#ifndef BT_DEFORMABLE_VBD_GPU_API_H
#define BT_DEFORMABLE_VBD_GPU_API_H
// Versioned POD boundary: the optional CUDA DLL does not depend on Bullet's scalar ABI.
struct btVbdGpuVec
{
	double x, y, z;
};
struct btVbdGpuMapping
{
	int begin, end;
	btVbdGpuVec offset;
};
struct btVbdGpuSupport
{
	int node;
	double j[9];
};
struct btVbdGpuApi
{
	void *(*create)();
	void (*destroy)(void *);
	int (*mapping)(void *, int, const btVbdGpuMapping *, int, const btVbdGpuSupport *);
	int (*guard)(void *, int, const btVbdGpuVec *, const btVbdGpuVec *, const btVbdGpuVec *, int, double, double *);
	const char *(*error)(void *);
	int (*evaluate)(void *, int, const btVbdGpuVec *, int, btVbdGpuVec *);
};
static_assert(sizeof(btVbdGpuVec) == 24 && sizeof(btVbdGpuMapping) == 32 && sizeof(btVbdGpuSupport) == 80, "GPU guard ABI layout mismatch");
using btVbdGpuGetApi = const btVbdGpuApi *(*)(int);
#endif
