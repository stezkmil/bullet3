#ifndef BT_DEFORMABLE_VBD_GPU_H
#define BT_DEFORMABLE_VBD_GPU_H
#include "btDeformableVbdGpuApi.h"
#include <string>
#include <vector>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif
// Optional backend loaded only when explicitly requested; CPU operation needs no CUDA installation.
class btDeformableVbdGpu
{
	const btVbdGpuApi *api = nullptr;
	void *state = nullptr;
	bool attempted = false;
	bool mappingValid = false;
#ifdef _WIN32
	HMODULE library = nullptr;
#endif
	std::string failure;
	std::vector<btVbdGpuVec> surface;

  public:
	btDeformableVbdGpu() = default;
	btDeformableVbdGpu(const btDeformableVbdGpu &) = delete;
	btDeformableVbdGpu &operator=(const btDeformableVbdGpu &) = delete;
	~btDeformableVbdGpu()
	{
		if (state)
			api->destroy(state);
#ifdef _WIN32
		if (library)
			FreeLibrary(library);
#endif
	}
	void invalidateMapping()
	{
		mappingValid = false;
	}
	bool hasMapping() const
	{
		return mappingValid;
	}
	const std::string &error() const
	{
		return failure;
	}
	bool ready()
	{
		if (attempted)
			return state && failure.empty();
		attempted = true;
#ifdef _WIN32
		wchar_t path[32768];
		DWORD length = GetModuleFileNameW(nullptr, path, 32768);
		if (!length || length >= 32768)
		{
			failure = "Cannot resolve executable directory";
			return false;
		}
		std::wstring filename(path, length);
		filename.resize(filename.find_last_of(L"\\/") + 1);
		filename += L"BulletVbdGpu.dll";
		library = LoadLibraryExW(filename.c_str(), nullptr, LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR | LOAD_LIBRARY_SEARCH_SYSTEM32);
		if (!library)
		{
			failure = "BulletVbdGpu.dll or its CUDA driver dependency is unavailable";
			return false;
		}
		auto get = reinterpret_cast<btVbdGpuGetApi>(GetProcAddress(library, "btVbdGpuGetApiV2"));
		api = get ? get(2) : nullptr;
		if (!api)
		{
			failure = "Incompatible GPU guard backend";
			return false;
		}
		state = api->create();
		if (!state)
		{
			failure = api->error(nullptr);
			return false;
		}
		return true;
#else
		failure = "GPU guard prototype is currently available on Windows only";
		return false;
#endif
	}
	bool mapping(const std::vector<btVbdGpuMapping> &maps, const std::vector<btVbdGpuSupport> &supports)
	{
		if (!ready())
			return false;
		if (!api->mapping(state, int(maps.size()), maps.data(), int(supports.size()), supports.data()))
		{
			failure = api->error(state);
			return false;
		}
		mappingValid = true;
		return true;
	}
	const std::vector<btVbdGpuVec> &positions() const
	{
		return surface;
	}
	bool evaluate(const std::vector<btVbdGpuVec> &x, int count)
	{
		if (!ready())
			return false;
		surface.resize(count);
		if (!api->evaluate(state, int(x.size()), x.data(), count, surface.data()))
		{
			failure = api->error(state);
			return false;
		}
		return true;
	}
	bool guard(const std::vector<btVbdGpuVec> &x, const std::vector<btVbdGpuVec> &p, const std::vector<btVbdGpuVec> &r, bool refresh,
			   double bound, double &limit)
	{
		if (!ready())
			return false;
		if (!api->guard(state, int(x.size()), x.data(), p.data(), r.data(), refresh ? 1 : 0, bound, &limit))
		{
			failure = api->error(state);
			return false;
		}
		return true;
	}
};
#endif
