#include "btDeformableVbdGpuApi.h"
#include "guard_ptx.h"
#include <cuda.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
namespace
{
struct Buffer
{
	CUdeviceptr ptr = 0;
	size_t capacity = 0;
};
struct ScalarSupport
{
	int node;
	double weight;
};
static_assert(sizeof(ScalarSupport) == 16, "CUDA scalar support layout");
struct State
{
	CUdevice device = 0;
	CUcontext context = nullptr;
	CUmodule module = nullptr;
	CUfunction kernel = nullptr, reduce = nullptr, evaluateKernel = nullptr;
	CUfunction scalarKernel = nullptr, scalarEvaluateKernel = nullptr;
	bool scalarMapping = false;
	Buffer maps, supports, x, p, r, reference, blocks, result, positions;
	std::vector<btVbdGpuMapping> hostMaps;
	std::vector<btVbdGpuSupport> hostSupports;
	bool mappingValid = false, referenceValid = false;
	int maximumNode = -1;
	std::string error;
	bool check(CUresult code)
	{
		if (code == CUDA_SUCCESS)
			return true;
		const char *name = nullptr;
		cuGetErrorName(code, &name);
		error = name ? name : "CUDA failure";
		return false;
	}
	bool reserve(Buffer &b, size_t bytes)
	{
		if (bytes <= b.capacity)
			return true;
		if (b.ptr)
		{
			cuMemFree(b.ptr);
			b = {};
		}
		if (!check(cuMemAlloc(&b.ptr, bytes)))
			return false;
		b.capacity = bytes;
		return true;
	}
	bool upload(Buffer &b, const void *data, size_t bytes)
	{
		return reserve(b, bytes) && (!bytes || check(cuMemcpyHtoD(b.ptr, data, bytes)));
	}
};
thread_local std::string creationError;
struct Current
{
	State &s;
	bool valid;
	Current(State &v) : s(v), valid(s.check(cuCtxPushCurrent(s.context))) {}
	~Current()
	{
		if (valid)
		{
			CUcontext c;
			cuCtxPopCurrent(&c);
		}
	}
};
void destroy(void *opaque)
{
	if (!opaque)
		return;
	auto &s = *static_cast<State *>(opaque);
	if (s.context)
	{
		Current current(s);
		if (current.valid)
		{
			for (auto *b : {&s.maps, &s.supports, &s.x, &s.p, &s.r, &s.reference, &s.blocks, &s.result, &s.positions})
				if (b->ptr)
					cuMemFree(b->ptr);
			if (s.module)
				cuModuleUnload(s.module);
		}
	}
	if (s.context)
		cuDevicePrimaryCtxRelease(s.device);
	delete &s;
}
void *create()
{
	auto *s = new State();
	if (!s->check(cuInit(0)) || !s->check(cuDeviceGet(&s->device, 0)) || !s->check(cuDevicePrimaryCtxRetain(&s->context, s->device)))
	{
		creationError = s->error;
		destroy(s);
		return nullptr;
	}
	bool ok = false;
	{
		Current current(*s);
		ok = current.valid && s->check(cuModuleLoadData(&s->module, vbdGuardPtx)) &&
			 s->check(cuModuleGetFunction(&s->kernel, s->module, "guard")) &&
			 s->check(cuModuleGetFunction(&s->reduce, s->module, "reduce")) &&
			 s->check(cuModuleGetFunction(&s->evaluateKernel, s->module, "evaluate")) &&
			 s->check(cuModuleGetFunction(&s->scalarKernel, s->module, "guardScalar")) &&
			 s->check(cuModuleGetFunction(&s->scalarEvaluateKernel, s->module, "evaluateScalar"));
	}
	if (!ok)
	{
		creationError = s->error;
		destroy(s);
		return nullptr;
	}
	return s;
}
int mapping(void *opaque, int count, const btVbdGpuMapping *maps, int supportCount, const btVbdGpuSupport *supports)
{
	auto &s = *static_cast<State *>(opaque);
	Current current(s);
	if (!current.valid || count < 0 || supportCount < 0)
		return 0;
	if (s.mappingValid && s.hostMaps.size() == size_t(count) && s.hostSupports.size() == size_t(supportCount) &&
		(!count || !std::memcmp(maps, s.hostMaps.data(), count * sizeof(*maps))) &&
		(!supportCount || !std::memcmp(supports, s.hostSupports.data(), supportCount * sizeof(*supports))))
		return 1;
	s.mappingValid = s.referenceValid = false;
	s.maximumNode = -1;
	for (int i = 0; i < count; ++i)
		if (maps[i].begin < 0 || maps[i].end < maps[i].begin || maps[i].end > supportCount)
		{
			s.error = "Invalid mapping support range";
			return 0;
		}
	for (int i = 0; i < supportCount; ++i)
	{
		if (supports[i].node < 0)
		{
			s.error = "Invalid mapping node";
			return 0;
		}
		s.maximumNode = std::max(s.maximumNode, supports[i].node);
	}
	s.scalarMapping = true;
	for (int i = 0; i < supportCount; ++i)
	{
		const double *j = supports[i].j;
		if (!std::isfinite(j[0]) || j[0] != j[4] || j[0] != j[8] || j[1] != 0 || j[2] != 0 || j[3] != 0 || j[5] != 0 || j[6] != 0 ||
			j[7] != 0)
		{
			s.scalarMapping = false;
			break;
		}
	}
	std::vector<ScalarSupport> compact;
	if (s.scalarMapping)
	{
		compact.reserve(supportCount);
		for (int i = 0; i < supportCount; ++i)
			compact.push_back({supports[i].node, supports[i].j[0]});
	}
	const void *supportData = s.scalarMapping ? static_cast<const void *>(compact.data()) : supports;
	const size_t supportBytes = supportCount * (s.scalarMapping ? sizeof(ScalarSupport) : sizeof(*supports));
	if (!s.upload(s.maps, maps, count * sizeof(*maps)) || !s.upload(s.supports, supportData, supportBytes) ||
		!s.reserve(s.reference, count * sizeof(btVbdGpuVec)))
		return 0;
	if (count)
		s.hostMaps.assign(maps, maps + count);
	else
		s.hostMaps.clear();
	if (supportCount)
		s.hostSupports.assign(supports, supports + supportCount);
	else
		s.hostSupports.clear();
	s.mappingValid = true;
	return 1;
}
int guard(void *opaque, int nodes, const btVbdGpuVec *x, const btVbdGpuVec *p, const btVbdGpuVec *r, int refresh, double bound, double *out)
{
	auto &s = *static_cast<State *>(opaque);
	Current current(s);
	if (!current.valid)
		return 0;
	if (!s.mappingValid || nodes <= s.maximumNode || nodes <= 0 || !(bound > 0))
	{
		s.error = "Invalid node count or displacement bound";
		return 0;
	}
	refresh = refresh || !s.referenceValid;
	int count = int(s.hostMaps.size()), blocks = (count + 255) / 256;
	if (!count)
	{
		*out = 1;
		return 1;
	}
	if (!s.upload(s.x, x, nodes * sizeof(*x)) || !s.upload(s.p, p, nodes * sizeof(*p)) ||
		(refresh && !s.upload(s.r, r, nodes * sizeof(*r))) || !s.reserve(s.blocks, blocks * sizeof(double)) ||
		!s.reserve(s.result, sizeof(double)))
		return 0;
	void *args[] = {&count, &s.maps.ptr, &s.supports.ptr, &s.x.ptr, &s.p.ptr, &s.r.ptr, &s.reference.ptr, &refresh, &bound, &s.blocks.ptr};
	void *reduceArgs[] = {&blocks, &s.blocks.ptr, &s.result.ptr};
	if (!s.check(cuLaunchKernel(s.scalarMapping ? s.scalarKernel : s.kernel, blocks, 1, 1, 256, 1, 1, 0, nullptr, args, nullptr)) ||
		!s.check(cuLaunchKernel(s.reduce, 1, 1, 1, 256, 1, 1, 0, nullptr, reduceArgs, nullptr)) ||
		!s.check(cuMemcpyDtoH(out, s.result.ptr, sizeof(double))))
		return 0;
	if (!std::isfinite(*out) || *out < 0 || *out > 1)
	{
		s.error = "Non-finite mapped guard input or output";
		return 0;
	}
	s.referenceValid = true;
	return 1;
}
int evaluate(void *opaque, int nodes, const btVbdGpuVec *x, int count, btVbdGpuVec *out)
{
	auto &s = *static_cast<State *>(opaque);
	Current current(s);
	if (!current.valid)
		return 0;
	if (!s.mappingValid || nodes <= s.maximumNode || nodes < 0 || count != int(s.hostMaps.size()))
	{
		s.error = "Invalid mapped evaluation dimensions";
		return 0;
	}
	if (!count)
		return 1;
	if (!s.upload(s.x, x, nodes * sizeof(*x)) || !s.reserve(s.positions, count * sizeof(*out)))
		return 0;
	void *args[] = {&count, &s.maps.ptr, &s.supports.ptr, &s.x.ptr, &s.positions.ptr};
	if (!s.check(cuLaunchKernel(s.scalarMapping ? s.scalarEvaluateKernel : s.evaluateKernel, (count + 255) / 256, 1, 1, 256, 1, 1, 0,
								nullptr, args, nullptr)) ||
		!s.check(cuMemcpyDtoH(out, s.positions.ptr, count * sizeof(*out))))
		return 0;
	for (int i = 0; i < count; ++i)
		if (!std::isfinite(out[i].x) || !std::isfinite(out[i].y) || !std::isfinite(out[i].z))
		{
			s.error = "Non-finite mapped surface position";
			return 0;
		}
	return 1;
}
const char *error(void *opaque)
{
	return opaque ? static_cast<State *>(opaque)->error.c_str() : creationError.c_str();
}
} // namespace
extern "C" __declspec(dllexport) const btVbdGpuApi *btVbdGpuGetApiV2(int version)
{
	static const btVbdGpuApi api{create, destroy, mapping, guard, error, evaluate};
	return version == 2 ? &api : nullptr;
}
