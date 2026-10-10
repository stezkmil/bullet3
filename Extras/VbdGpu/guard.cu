
struct Vec
{
	double x, y, z;
};
struct Mapping
{
	int begin, end;
	Vec offset;
};
struct Support
{
	int node;
	double j[9];
};
struct ScalarSupport
{
	int node;
	double weight;
};
__device__ Vec add(Vec a, Vec b)
{
	return {a.x + b.x, a.y + b.y, a.z + b.z};
}
__device__ Vec sub(Vec a, Vec b)
{
	return {a.x - b.x, a.y - b.y, a.z - b.z};
}
__device__ double dot(Vec a, Vec b)
{
	return (a.x * b.x + a.y * b.y) + a.z * b.z;
}
__device__ Vec position(Mapping m, const Support *s, const Vec *x)
{
	Vec p = m.offset;
	for (int k = m.begin; k < m.end; ++k)
	{
		Support a = s[k];
		Vec v = x[a.node];
		p = add(p, {(a.j[0] * v.x + a.j[1] * v.y) + a.j[2] * v.z, (a.j[3] * v.x + a.j[4] * v.y) + a.j[5] * v.z,
					(a.j[6] * v.x + a.j[7] * v.y) + a.j[8] * v.z});
	}
	return p;
}
// Scalar interpolation needs one coefficient rather than a full Jacobian.
__device__ Vec position(Mapping m, const ScalarSupport *s, const Vec *x)
{
	Vec p = m.offset;
	for (int k = m.begin; k < m.end; ++k)
	{
		ScalarSupport a = s[k];
		Vec v = x[a.node];
		p = add(p, {a.weight * v.x, a.weight * v.y, a.weight * v.z});
	}
	return p;
}
template <class T>
__device__ void guardImpl(int count, const Mapping *maps, const T *s, const Vec *x, const Vec *p, const Vec *r, Vec *reference, int refresh,
						  double bound, double *blocks)
{
	__shared__ double limits[256];
	int v = blockIdx.x * blockDim.x + threadIdx.x;
	double limit = 1;
	if (v < count)
	{
		Mapping m = maps[v];
		if (refresh)
			reference[v] = position(m, s, r);
		bool changed = false;
		for (int k = m.begin; k < m.end; ++k)
		{
			int i = s[k].node;
			changed = changed || x[i].x != p[i].x || x[i].y != p[i].y || x[i].z != p[i].z;
		}
		if (changed)
		{
			Vec start = position(m, s, x), d = sub(position(m, s, p), start), a0 = sub(start, reference[v]);
			if (dot(add(a0, d), add(a0, d)) > bound * bound && dot(d, d) > 0)
			{
				double a = dot(d, d), b = dot(a0, d), c = dot(a0, a0) - bound * bound;
				limit = fmin(1., fmax(0., (-b + sqrt(fmax(0., b * b - a * c))) / a));
			}
			if (!isfinite(dot(start, start)) || !isfinite(dot(d, d)) || !isfinite(dot(a0, a0)))
				limit = -1;
		}
	}
	limits[threadIdx.x] = limit;
	__syncthreads();
	for (int stride = 128; stride; stride /= 2)
	{
		if (threadIdx.x < stride)
			limits[threadIdx.x] = fmin(limits[threadIdx.x], limits[threadIdx.x + stride]);
		__syncthreads();
	}
	if (threadIdx.x == 0)
		blocks[blockIdx.x] = limits[0];
}
extern "C" __global__ void guard(int count, const Mapping *maps, const Support *s, const Vec *x, const Vec *p, const Vec *r, Vec *reference,
								 int refresh, double bound, double *blocks)
{
	guardImpl(count, maps, s, x, p, r, reference, refresh, bound, blocks);
}
extern "C" __global__ void guardScalar(int count, const Mapping *maps, const ScalarSupport *s, const Vec *x, const Vec *p, const Vec *r,
									   Vec *reference, int refresh, double bound, double *blocks)
{
	guardImpl(count, maps, s, x, p, r, reference, refresh, bound, blocks);
}
extern "C" __global__ void reduce(int count, const double *blocks, double *out)
{
	__shared__ double limits[256];
	double v = 1;
	for (int i = threadIdx.x; i < count; i += 256)
		v = fmin(v, blocks[i]);
	limits[threadIdx.x] = v;
	__syncthreads();
	for (int stride = 128; stride; stride /= 2)
	{
		if (threadIdx.x < stride)
			limits[threadIdx.x] = fmin(limits[threadIdx.x], limits[threadIdx.x + stride]);
		__syncthreads();
	}
	if (threadIdx.x == 0)
		*out = limits[0];
}

extern "C" __global__ void evaluate(int count, const Mapping *maps, const Support *s, const Vec *x, Vec *out)
{
	int v = blockIdx.x * blockDim.x + threadIdx.x;
	if (v < count)
		out[v] = position(maps[v], s, x);
}

extern "C" __global__ void evaluateScalar(int count, const Mapping *maps, const ScalarSupport *s, const Vec *x, Vec *out)
{
	int v = blockIdx.x * blockDim.x + threadIdx.x;
	if (v < count)
		out[v] = position(maps[v], s, x);
}
