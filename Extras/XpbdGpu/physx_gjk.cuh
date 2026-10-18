// Adapted from NVIDIA PhysX CUDA GJK/EPA (Apache-2.0; see LICENSE-PhysX.txt).
// Standalone math/storage adapter; the warp algorithms below retain their source structure.
#define FULL_MASK 0xffffffffu
#define PX_MAX_F32 3.402823466e38f
#define FLT_MAX PX_MAX_F32
#define PX_CUDA_CALLABLE __device__
#define PX_FORCE_INLINE __forceinline__
#define GPU_NP_CHECK(condition) do { if (!(condition)) asm("trap;"); } while(0)
namespace physx {
using PxReal=float; using PxU32=unsigned int; using PxU16=unsigned short; using PxU8=unsigned char;
__device__ float PxAbs(float x) { return fabsf(x); }
template<class T> __device__ T PxMin(T a,T b) { return a<b?a:b; }
__device__ unsigned lowestSetBit(unsigned x) { return x&(~x+1); }
__device__ unsigned lowestSetIndex(unsigned x) { return __ffs(x)-1; }
__device__ unsigned highestSetIndex(unsigned x) { return 31-__clz(x); }
struct PxVec3 {
 float x,y,z;
 __device__ PxVec3() {}
 __device__ explicit PxVec3(float a):x(a),y(a),z(a) {}
 __device__ PxVec3(float a,float b,float c):x(a),y(b),z(c) {}
 __device__ PxVec3 operator+(const PxVec3& b) const {return {x+b.x,y+b.y,z+b.z};}
 __device__ PxVec3 operator-(const PxVec3& b) const {return {x-b.x,y-b.y,z-b.z};}
 __device__ PxVec3 operator-() const {return {-x,-y,-z};}
 __device__ PxVec3 operator*(float a) const {return {x*a,y*a,z*a};}
 __device__ PxVec3 operator/(float a) const {return *this*(1.f/a);}
 __device__ float dot(const PxVec3& b) const {return x*b.x+y*b.y+z*b.z;}
 __device__ PxVec3 cross(const PxVec3& b) const {return {y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x};}
 __device__ float magnitudeSquared() const {return dot(*this);}
 __device__ float magnitude() const {return sqrtf(magnitudeSquared());}
 __device__ PxVec3 getNormalized() const {float m=magnitude();return m>0?*this/m:PxVec3(0);}
 __device__ PxVec3 abs() const {return {fabsf(x),fabsf(y),fabsf(z)};}
};
__device__ PxVec3 operator*(float a,const PxVec3& b) {return b*a;}
__device__ PxVec3 ldS(const volatile PxVec3& v) { return {v.x,v.y,v.z}; }
__device__ void stS(volatile PxVec3& v,const PxVec3& p) { v.x=p.x;v.y=p.y;v.z=p.z; }
}
// Copyright (c) 2001-2004 NovodeX AG. All rights reserved.
// Copyright (c) 2004-2008 AGEIA Technologies, Inc. All rights reserved.
// SPDX-FileCopyrightText: Copyright (c) 2008-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef GU_SCHLOCKSHARED_H
#define GU_SCHLOCKSHARED_H









#define GJK_LOGGING 0
#define EPA_REPLACE_TRIANGLE_LOGGING 0
#define EPA_LOGGING 0
#define EPA_LOG_DEGENERATE_CASES 0

#if GJK_LOGGING
#define GJK_LOG(...) if(gjkDebug) printf(__VA_ARGS__)
#else
#define GJK_LOG(...)
#endif


namespace schlock
{


using namespace physx;

template<typename T>
PX_CUDA_CALLABLE PX_FORCE_INLINE void swap(T& a, T& b)
{
	T tmp = a;
	a = b;
	b = tmp;
}

PX_FORCE_INLINE PX_CUDA_CALLABLE PxU32 encodeIndices(PxU8 a, PxU8 b)		{ return 0x3f800000u | b<<8 | a;}
PX_FORCE_INLINE PX_CUDA_CALLABLE PxU8 aIndex(PxU32 indices)					{ return indices&0xff; }
PX_FORCE_INLINE PX_CUDA_CALLABLE PxU8 bIndex(PxU32 indices)					{ return (indices>>8)&0xff; }

struct CachedVertex
{
	PxVec3 v;
	PxU32 indices;
	PX_FORCE_INLINE PX_CUDA_CALLABLE CachedVertex() {}
	PX_FORCE_INLINE PX_CUDA_CALLABLE CachedVertex(PxU8 iA, PxU8 iB, const PxVec3& iV): v(iV), indices(encodeIndices(iA,iB)) {}
	PX_FORCE_INLINE PX_CUDA_CALLABLE CachedVertex(const CachedVertex& other): v(other.v), indices(other.indices) {}
	PX_FORCE_INLINE PX_CUDA_CALLABLE void operator=(const CachedVertex& other)	{	v = other.v, indices = other.indices; }

	PX_FORCE_INLINE PX_CUDA_CALLABLE PxU8 a() const { return aIndex(indices); }
	PX_FORCE_INLINE PX_CUDA_CALLABLE PxU8 b() const { return bIndex(indices); }

	PX_FORCE_INLINE PX_CUDA_CALLABLE PxU8 a() const volatile { return aIndex(indices); }
	PX_FORCE_INLINE PX_CUDA_CALLABLE PxU8 b() const volatile { return bIndex(indices); }
	PX_FORCE_INLINE PX_CUDA_CALLABLE PxVec3 getV() const volatile { return PxVec3(v.x, v.y, v.z); }
};

struct GjkCachedData
{
	CachedVertex vertices[4];
	struct
	{
		PxVec3 n;
		PxReal d;
	} normals[4];

	PxU32 size;

	PX_FORCE_INLINE PX_CUDA_CALLABLE GjkCachedData(){}

	__device__ void Reset() { size = 0; __syncwarp(); }
};

struct GjkResult
{
	enum Enum
	{
		eDISTANT,
		eCLOSE,
		eOVERLAP
	};
};

struct GjkOutput
{
	PxVec3 direction;				// unit witness for lower bound
	PxReal lowerBound;				// provable lower bound
	PxReal upperBound;				// provable upper bound

	PxVec3 closestPointA;			// witnesses for upper bound
	PxVec3 closestPointB;

	bool degenerate;
	PxVec3 closestPointDir;
};

PX_CUDA_CALLABLE PX_FORCE_INLINE void barycentrics(const PxVec3& a, const PxVec3& b, const PxVec3& c, float &u, float &v, float &w, const PxVec3& p)
{
    PxVec3 m = (b - a).cross(c - a), ma = m.abs(), a1 = (b-p).cross(c-p), a2 = (c-p).cross(a-p);
    PxReal d;

    // area ratios in plane with largest normal component
    if (ma.x >= ma.y && ma.x >= ma.z)
		d = 1.0f/m.x, u = a1.x * d,  v = a2.x * d;
    else if (ma.y >= ma.z)
		d = 1.0f/m.y, u = a1.y * d,  v = a2.y * d;
	else
		d = 1.0f/m.z, u = a1.z * d,  v = a2.z * d;

    w = 1.0f - u - v;
}


// top-level functions

// unit-level functions for testing
static const PxU32 EPA_MAX_VERTS = 18, EPA_MAX_FACES = 2*EPA_MAX_VERTS-4;		// for N verts, we'll have 2(N-2) faces.

typedef int EPABitMap;

#define MAKE_Vec3I8(A, B, C) (physx::PxU32(A) | (physx::PxU32(B) << 8) | (physx::PxU32(C) << 16))

struct Vec3I8
{
	PxU32 a;
	PX_FORCE_INLINE PX_CUDA_CALLABLE PxU8 get(PxU32 i)	const			{ return (a >> (i << 3)) & 0xFFul; }
	PX_FORCE_INLINE PX_CUDA_CALLABLE void set(PxU32 i, PxU8 v)			{ a = (a & ~(0xFFul << (i << 3))) | (PxU32(v) << (i << 3)); }
	PX_FORCE_INLINE PX_CUDA_CALLABLE void setInt(int i)					{ a = i; }
	PX_FORCE_INLINE PX_CUDA_CALLABLE int getInt()						{ return a; }

	PX_FORCE_INLINE PX_CUDA_CALLABLE const Vec3I8& operator=(const Vec3I8& other)						{a = other.a; return *this;}

	__device__ bool operator==(const Vec3I8& other) { return (a & 0xFFffFFul) == (other.a & 0xFFffFFul); }
};

struct EPATriangle
{
	PxVec3 normal;			// unit normal
	PxReal distance;		// *plane* distance from origin
	Vec3I8 v;				// vertices
};

struct EPAHullData
{
	PxVec3				v[EPA_MAX_VERTS];
	PxU32				i[EPA_MAX_VERTS];
	EPATriangle			t[EPA_MAX_FACES];
	Vec3I8				a[EPA_MAX_FACES];

	PxU8				nbVertices;
	EPABitMap			liveMap;
};



}
#endif

// Copyright (c) 2001-2004 NovodeX AG. All rights reserved.
// Copyright (c) 2004-2008 AGEIA Technologies, Inc. All rights reserved.
// SPDX-FileCopyrightText: Copyright (c) 2008-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef __NPUTILS_CUH__
#define __NPUTILS_CUH__







namespace physx
{

static __device__ __forceinline__ int next3(int x)
{
	return x==2 ? 0 : x+1;
}

static __device__ inline int prev3(int x)
{
	return x==0 ? 2 : x-1;
}


// operators dealing with 8 4-vectors in a warp
static __device__ inline float dot(float a, float b)	// leaves answer in x
{
	float x = a*b;
	return x + __shfl_down_sync(FULL_MASK, x,1) + __shfl_down_sync(FULL_MASK, x,2);
}

static __device__ inline float cross(float a, float b, int cp0, int cp1)
{
	return __shfl_sync(FULL_MASK, a, cp0) * __shfl_sync(FULL_MASK, b, cp1) - __shfl_sync(FULL_MASK, a, cp1) * __shfl_sync(FULL_MASK, b, cp0);
}

static __device__ inline float splatV(int index, float v)	// splat a vector across all 8 vectors
{
	return __shfl_sync(FULL_MASK, v, (index << 2) | threadIdx.x & 3);
}

static __device__ inline float splatX(float a)				// splat the x-component of each vector into yzw
{
	return __shfl_sync(FULL_MASK, a, threadIdx.x & 28);
}

static __device__ inline float normalize(float x)			// normalize each vector
{
	return x * splatX(rsqrt(dot(x,x)));		// TODO: precision
}

static __device__ inline float magnitude(float x)			// normalize each vector
{
	return splatX(sqrt(dot(x,x)));		// TODO: precision
}

template<int N>
static __device__ inline int packXFlags(bool a)				// pack a boolean from each x into a byte
{
	return __ballot_sync(FULL_MASK, __shfl_sync(FULL_MASK, a, threadIdx.x<<2)) & ((1<<N)-1);
}

template<int N>
static __device__ inline bool allX(bool a)					// check some number initial set of the x-flags are all set
{
	return packXFlags<N>(a) == ((1<<N)-1);
}

template<int N>
static __device__ inline bool anyX(bool a)					// check if any of an initial segment of x-flags are all set
{
	return packXFlags<N>(a) != 0;
}

static __device__ inline float loadV3Unsafe(const PxVec3* v)	// NB requires padding (i.e. the 4th element is loaded)
{
	return reinterpret_cast<const float*>(v)[threadIdx.x&3];
}

static __device__ void storeV3(volatile PxVec3* v, int index, float value)
{
	int lane = index<<2;
	volatile float* dest = reinterpret_cast<volatile float*>(v);
	if(threadIdx.x>=lane && threadIdx.x<lane+3)
	{
		dest[threadIdx.x&3] = value;
	}
}

// experimentally, seems more register-efficient to coalesce this
static __device__ __forceinline__
PxReal shuffleDot(const PxU32 syncMask, const PxVec3& v0, int shuffle0, const PxVec3& v1)
{
	return __shfl_sync(syncMask, v0.x, shuffle0)*v1.x + __shfl_sync(syncMask, v0.y, shuffle0)*v1.y + __shfl_sync(syncMask, v0.z, shuffle0)*v1.z;
}

static __device__ __forceinline__
PxU32 maxIndex(PxReal v, PxU32 mask, PxReal& maxV)
{
	maxV = mask & (1 << (threadIdx.x&31)) ? v : -FLT_MAX;

	maxV = fmaxf(maxV, __shfl_xor_sync(FULL_MASK, maxV, 16));
	maxV = fmaxf(maxV, __shfl_xor_sync(FULL_MASK, maxV, 8));
	maxV = fmaxf(maxV, __shfl_xor_sync(FULL_MASK, maxV, 4));
	maxV = fmaxf(maxV, __shfl_xor_sync(FULL_MASK, maxV, 2));
	maxV = fmaxf(maxV, __shfl_xor_sync(FULL_MASK, maxV, 1));

	return lowestSetIndex(__ballot_sync(FULL_MASK, maxV == v)&mask);
}

static __device__ __forceinline__
PxU32 minIndex(PxReal v, PxU32 mask, PxReal& minV)
{
	minV = mask & (1 << (threadIdx.x & 31)) ? v : FLT_MAX;

	minV = fminf(minV, __shfl_xor_sync(FULL_MASK, minV, 16));
	minV = fminf(minV, __shfl_xor_sync(FULL_MASK, minV, 8));
	minV = fminf(minV, __shfl_xor_sync(FULL_MASK, minV, 4));
	minV = fminf(minV, __shfl_xor_sync(FULL_MASK, minV, 2));
	minV = fminf(minV, __shfl_xor_sync(FULL_MASK, minV, 1));

	return lowestSetIndex(__ballot_sync(FULL_MASK, minV == v)&mask);
}

// similar as above but only with blocks of 4 threads
// mask must have exactly 4 consecutive bits set and the corresponding
// threads must execute this function together
static __device__ __forceinline__
PxU32 minIndex4(PxReal v, PxU32 mask, PxU32 threadIndexInGroup)
{
	// sanity checks that function has been called as expected
	GPU_NP_CHECK(__popc(mask) == 4); // exactly 4 threads in mask
	GPU_NP_CHECK(mask & (1 << (threadIdx.x & 31))); // executing thread must be in the mask
	GPU_NP_CHECK(threadIndexInGroup == 0 || (mask & (1 << (threadIdx.x - 1 & 31)))); // threads consecutive

	PxReal minV =  v;

	PxReal minV_m1 = __shfl_sync(mask, minV, (threadIdx.x & 31) - 1);
	if(threadIndexInGroup > 0)
		minV = fminf(minV, minV_m1);

	PxReal minV_m2 = __shfl_sync(mask, minV, (threadIdx.x & 31) - 2);
	if(threadIndexInGroup > 2)
		minV = fminf(minV, minV_m2);

	// Now the last thread (idx 3) of the group knows the min.
	// Send it back to all the threads.
	minV = __shfl_sync(mask, minV, (threadIdx.x & 31) | 3);

	return lowestSetIndex(__ballot_sync(mask, minV == v)&mask);
}

} // namespace physx

#endif

// Copyright (c) 2001-2004 NovodeX AG. All rights reserved.
// Copyright (c) 2004-2008 AGEIA Technologies, Inc. All rights reserved.
// SPDX-FileCopyrightText: Copyright (c) 2008-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef __EPA_CUH__
#define __EPA_CUH__






namespace physx
{

namespace squawk
{

struct EpaScratch
{
	PxVec3 vertices[schlock::EPA_MAX_VERTS];
	PxVec3 bestDir;
	PxReal bestDist;
	PxU16 indices[schlock::EPA_MAX_VERTS];
	PxU8 loop[schlock::EPA_MAX_VERTS];
	PxU8 edges[schlock::EPA_MAX_VERTS];
	PxU8 next[schlock::EPA_MAX_VERTS];
};


__constant__ __device__ schlock::Vec3I8 sStartVertices4[4] =		{	MAKE_Vec3I8(0,1,2),	MAKE_Vec3I8(3,1,0),	MAKE_Vec3I8(3,2,1),	MAKE_Vec3I8(3,0,2) };
__constant__ __device__ schlock::Vec3I8 sStartAdjacencies4[4] =	{	MAKE_Vec3I8(2,3,1),	MAKE_Vec3I8(0,3,2),	MAKE_Vec3I8(0,1,3),	MAKE_Vec3I8(0,2,1) };

__constant__ __device__ schlock::Vec3I8 sStartVertices3[2] =		{	MAKE_Vec3I8(0,1,2),	MAKE_Vec3I8(1,0,2) };
__constant__ __device__ schlock::Vec3I8 sStartAdjacency3[2] =	{	MAKE_Vec3I8(1,1,1),	MAKE_Vec3I8(0,0,0)	};



typedef int VMapType;

__device__ inline bool replaceTriangle(EpaScratch& epaS, int nbVertices, int triIndex,
												PxReal lowerBound, PxVec3& normal, PxReal& distance,
												schlock::Vec3I8& v, schlock::Vec3I8& a, int& liveB, const PxReal releps)
{
	// to make it easier to keep track, we use a 'B' suffix for a bitmap where each bit represents a lane (and therefore a triangle)
	int tI = threadIdx.x, laneB = 1<<tI;
	//PxVec3 apex = ldS(S.vertices[nbVertices-1]);
	PxVec3 apex = epaS.vertices[nbVertices - 1];

	// triangles that are silhouetted by the new vertex are threatened, but not yet killed
	bool threatened = liveB&laneB && (normal.dot(apex) - distance >= /*-*/releps);
	schlock::EPABitMap threatenedB = __ballot_sync(FULL_MASK, threatened);

	if (threatenedB == 0 || threatenedB == liveB || !(threatenedB & (1<<triIndex)))
	{
#if EPA_LOG_DEGENERATE_CASES
		if(threadIdx.x == 0)
			printf("degenrate: threatenedB %d liveB %d < tri %d \n", threatenedB, liveB, tri);
#endif

		return false;
	}

	// the threatened set is not guaranteed connected, so flood fill to find killB, a bitmap of tris to remove
	schlock::EPABitMap b0B = 1<<a.get(0), b1B = 1<<a.get(1), b2B = 1<<a.get(2);
	// nbr = threatened neighbors of threatened triangles
	schlock::EPABitMap killB = 1<<triIndex, nbr = (b0B | b1B | b2B) & (threatened?threatenedB:0);

	// flood the threatened triangles to find the killed ones
	//for(schlock::EPABitMap pB = 0; pB != killB; pB = killB, killB |= __ballot(nbr&killB));

	for (schlock::EPABitMap pB = 0; pB != killB;)
	{
		pB = killB;
		killB |= __ballot_sync(FULL_MASK, nbr&killB);
	}

	//every triangles aren't alive and the triangles we want to kill in this iteration
	schlock::EPABitMap deadB = (~liveB) | killB;

	// d0, d1, d2 denote for each triangle whether the corresponding edge is a silhouette edge
	const int alive = laneB&(~deadB);
	const int d0 = alive && (b0B&deadB);
	const int d1 = alive && (b1B&deadB);
	const int d2 = alive && (b2B&deadB);

	// write the silhouette edges into the loop buffer
	int v0 = v.get(0), v1 = v.get(1), v2 = v.get(2), z = 0;

	GPU_NP_CHECK(!alive || v0 < schlock::EPA_MAX_VERTS);
	GPU_NP_CHECK(!alive || v1 < schlock::EPA_MAX_VERTS);
	GPU_NP_CHECK(!alive || v2 < schlock::EPA_MAX_VERTS);

	// S.edges is the adjacency for the new cone tri on this edge, and
	// also where the cone triangle will write back its index for later fixup
	if (d0)
		z = v2, epaS.next[v2] = v1, epaS.edges[v2] = tI;
	if (d1)
		z = v0, epaS.next[v0] = v2, epaS.edges[v0] = tI;
	if (d2)
		z = v1, epaS.next[v1] = v0, epaS.edges[v1] = tI;

	int d0d1d2B = __ballot_sync(FULL_MASK, d0|d1|d2);

	if (!d0d1d2B)
	{
#if EPA_LOG_DEGENERATE_CASES
		if (threadIdx.x == 0)
			printf("degenrate: d0d1d2B %d b0B %d b1B %d b2B %d deadB %d \n", d0d1d2B, b0B, b1B, b2B, deadB);
#endif

		return false;
	}

	int loopStart = __shfl_sync(FULL_MASK, z, lowestSetIndex(d0d1d2B));		// pick a loop vert from the lowest indexed triangle that has one
	int nbLoop = 0;

	__syncwarp();

	if(tI == 0)
	{
		GPU_NP_CHECK(loopStart < schlock::EPA_MAX_VERTS);

		epaS.loop[nbLoop++] = loopStart;									// flatten the linked list into a vertex array
		for(int vertex = epaS.next[loopStart]; vertex!=loopStart; vertex = epaS.next[vertex])
		{
			if(vertex >= schlock::EPA_MAX_VERTS)
			{
#if EPA_LOG_DEGENERATE_CASES
				printf("degenrate: vertex %d\n", vertex);
#endif

				nbLoop = -1;
				break;
			}

			if(nbLoop >= schlock::EPA_MAX_VERTS - 1)
			{
#if EPA_LOG_DEGENERATE_CASES
				printf("degenrate: nbLoop %d\n", nbLoop);
#endif

				nbLoop = -1;
				break;
			}

			epaS.loop[nbLoop++] = vertex;
		}

#if EPA_REPLACE_TRIANGLE_LOGGING
		for(int i=0;i<nbLoop;i++)
			printf("%d, ", S.loop[i]);
		printf("(starting at %d)\n",loopStart);
#endif
	}


	nbLoop = __shfl_sync(FULL_MASK, nbLoop, 0);

	if(nbLoop != __popc(__ballot_sync(FULL_MASK, d0)) + __popc(__ballot_sync(FULL_MASK, d1)) + __popc(__ballot_sync(FULL_MASK, d2)))
	{
#if EPA_LOG_DEGENERATE_CASES
		{
			printf("nbLoop %d __ballot(d0) %d __ballot(d1) %d __ballot(d2) %d tid %d", nbLoop, __ballot(d0), __ballot(d1), __ballot(d2), threadIdx.x);
		}
#endif
		return false;
	}

	__syncwarp();

	// make the cone with the nbLoop least non-survivors
	int coneIndex = __popc(deadB&(laneB-1)), coneB = __ballot_sync(FULL_MASK, deadB&laneB && coneIndex < nbLoop);
	bool fail = false;
	if(coneB & laneB)
	{
		GPU_NP_CHECK(coneIndex < nbLoop);
		GPU_NP_CHECK(nbLoop < schlock::EPA_MAX_VERTS);

		int v0 = nbVertices-1, v1 = epaS.loop[coneIndex], v2 = epaS.loop[coneIndex+1 == nbLoop ? 0 : coneIndex+1];
		v.setInt(MAKE_Vec3I8(v0, v1, v2));

		GPU_NP_CHECK(v0 < schlock::EPA_MAX_VERTS);
		GPU_NP_CHECK(v1 < schlock::EPA_MAX_VERTS);
		GPU_NP_CHECK(v2 < schlock::EPA_MAX_VERTS);

		a.set(0, epaS.edges[v1]);
		a.set(1, lowestSetIndex(coneB & ~(coneIndex+1 == nbLoop ? 0 : 2*laneB-1)));		// previous bit in the active mask, wrapping around
		a.set(2, highestSetIndex(coneB & ((coneIndex==0 ? 0 : laneB)-1)));				// next bit in the active mask, wrapping around


		epaS.edges[v1] = tI;																// write our index back into the buffer for later fixup

		//PxVec3 p0 = ldS(S.vertices[v0]), p1 = ldS(S.vertices[v1]), p2 = ldS(S.vertices[v2]);
		PxVec3 p0 = epaS.vertices[v0], p1 = epaS.vertices[v1], p2 = epaS.vertices[v2];
		PxVec3 e01 = p1 - p0, e02 = p2 - p0;
		PxVec3 n;

		if (e01.magnitudeSquared() < e02.magnitudeSquared())
		{
			n = e01.cross(p2 - p1);
		}
		else
		{
			n = e02.cross(p2 - p1);
		}

		PxReal nm = n.magnitude();
		normal = n/nm;
		distance = normal.dot(p0);
		/*fail = nm == 0 || (lowerBound > 0.0f && distance < lowerBound);
		if (fail)
		{
			printf("fail nm %f\ lowerBound %f distance %f \n", nm, lowerBound, distance);
		}*/
		fail = nm == 0 || (lowerBound > 0.0f && (distance - lowerBound) < -1e-4f);
	}

	__syncwarp();

	// ... to fix up the silhouette adjacency with the new tris
	if(d0) a.set(0, epaS.edges[v.get(2)]);
	if(d1) a.set(1, epaS.edges[v.get(0)]);
	if(d2) a.set(2, epaS.edges[v.get(1)]);

	liveB = (liveB & ~deadB) | coneB;
	return __ballot_sync(FULL_MASK, fail) == 0;
}


__device__ inline void epaStart4(EpaScratch& S,
										  const schlock::GjkCachedData* data,
										  schlock::Vec3I8& a, schlock::Vec3I8& v,
										  PxVec3& normal,
										  PxReal& distance)
{
	int tI = threadIdx.x;
	PxReal input = reinterpret_cast<const PxReal*>(data)[tI];
	if(tI<16)
	{
		if((tI&3)<3)
			(reinterpret_cast<volatile PxReal*>(S.vertices))[tI-tI/4] = input;
		else
			S.indices[tI>>2] = __float_as_int(input);
	}

	PxReal d = dot(input, input), m = __shfl_sync(FULL_MASK, 1/sqrt(d), tI*4+16);

	normal.x = __shfl_sync(FULL_MASK, input, tI*4+16) * m;
	normal.y = __shfl_sync(FULL_MASK, input, tI*4+17) * m;
	normal.z = __shfl_sync(FULL_MASK, input, tI*4+18) * m;
	distance = __shfl_sync(FULL_MASK, input, tI*4+19) * m;
	if(tI<4)
		a = sStartAdjacencies4[tI], v = sStartVertices4[tI];
}


__device__ inline PxU8 supportIndex(const PxVec3* v, PxU8 nbVertices, const PxVec3& dir)
{
	int bestN = 0xFFFFFFFF;
	PxReal bestDist = -FLT_MAX;
	for(PxU32 i=0;i<nbVertices;i+=32)
	{
		int n = i+threadIdx.x;
		PxReal d = n < nbVertices ? dir.dot(v[n]) : -FLT_MAX, e;

		e = fmaxf(d, __shfl_xor_sync(FULL_MASK, d, 16));
		e = fmaxf(e, __shfl_xor_sync(FULL_MASK, e, 8));
		e = fmaxf(e, __shfl_xor_sync(FULL_MASK, e, 4));
		e = fmaxf(e, __shfl_xor_sync(FULL_MASK, e, 2));
		e = fmaxf(e, __shfl_xor_sync(FULL_MASK, e, 1));

		if(e>bestDist)
			bestDist = e, bestN = i+__ffs(__ballot_sync(FULL_MASK, e == d))-1;
	}

	GPU_NP_CHECK(bestN < nbVertices);
	return PxU8(bestN);
}

// find a new support point in either the positive or negative direction furthest from the reference point
// we could be a bit more efficient here but this is a super-rare case
__device__ inline void newPoint(const PxVec3* vA, PxU32 nbA, const PxVec3* vB, PxU32 nbB, const PxVec3& dir, const PxVec3& ref, volatile PxVec3& vertex, volatile PxU16& indices)
{
	PxU8 indexAp = supportIndex(vA, nbA, dir), indexBp = supportIndex(vB, nbB, -dir);	// support point in +ve dir
	PxU8 indexAn = supportIndex(vA, nbA, -dir), indexBn = supportIndex(vB, nbB, dir);	// support point in +ve dir
	PxVec3 vp = vA[indexAp] - vB[indexBp], vn = vA[indexAn] - vB[indexBn];

	if(threadIdx.x == 0)
	{
		if(PxAbs((vp-ref).dot(dir)) > PxAbs((vn-ref).dot(dir)))
		{
			stS(vertex, vp);
			indices = indexBp<<8 | indexAp;
		}
		else
		{
			stS(vertex, vn);
			indices = indexBn<<8 | indexAn;
		}
	}
}

// find a new support point in either the positive or negative direction furthest from the reference point
// we could be a bit more efficient here but this is a super-rare case

__device__ inline void epaStart123(volatile EpaScratch& S,
											const schlock::GjkCachedData* data,
											const PxVec3* vA, PxU32 nbA,
											const PxVec3* vB, PxU32 nbB,
											schlock::Vec3I8& a, schlock::Vec3I8& v,
											PxVec3& normal,
											PxReal& distance,
											PxU32 warpMask)
{
	int tI = threadIdx.x;
	PxReal input = reinterpret_cast<const PxReal*>(data)[tI];
	if(tI<16)
	{
		if((tI&3)<3)
			(reinterpret_cast<volatile PxReal*>(S.vertices))[tI-tI/4] = input;
		else
			S.indices[tI>>2] = __float_as_int(input);
	}

	__syncwarp(warpMask); //S.vertices is written above and read below

	PxVec3 vertex0 = ldS(S.vertices[0]);

	if(data->size<2)
		newPoint(vA, nbA, vB, nbB, PxVec3(1,0,0), ldS(S.vertices[0]), S.vertices[1], S.indices[1]);

	__syncwarp(warpMask); //S.vertices is written (inside newPoint) above and read below
	PxVec3 vertex1 = ldS(S.vertices[1]);

	if(data->size<3)
	{
		PxVec3 dir = vertex1 - vertex0, s = (PxAbs(dir.x) < PxAbs(dir.y) ? PxVec3(1,0,0) : PxVec3(0,1,0)).cross(dir);
		newPoint(vA, nbA, vB, nbB, s, vertex0, S.vertices[2], S.indices[2]);
	}

	__syncwarp(warpMask); //S.vertices is written (inside newPoint) above and read below
	PxVec3 vertex2 = ldS(S.vertices[2]);

	PxVec3 e01 = vertex1 - vertex0, e12 = vertex2 - vertex1, e02 = vertex2 - vertex0;

	if (e01.magnitude()<=e02.magnitude() )
	{
		normal = e01.cross(e12).getNormalized();
	}
	else
	{
		normal = e02.cross(e12).getNormalized();
	}

	distance = normal.dot(vertex0);

	if(tI < 2)
		v = sStartVertices3[tI], a = sStartAdjacency3[tI];
	if(tI==1)
		normal = -normal, distance = -distance;
}


__device__ inline schlock::GjkResult::Enum epa(EpaScratch& epaS,
								const PxVec3* vA, PxU32 nbA,
								const PxVec3* vB, PxU32 nbB,
								const schlock::GjkCachedData* input,
								PxReal convergenceRatio,
								PxReal relEpsBase,
								schlock::GjkOutput& output)
{
	PxReal releps = relEpsBase * (1-convergenceRatio);

	schlock::Vec3I8 a, v;					// per-triangle adjancencies and vertex indices
	PxVec3 normal;					// per-triangle normal
	PxReal distance;				// per-triangle distance
	schlock::Vec3I8 bestVerts;

	PxReal lowerBound = -PX_MAX_F32, upperBound = PX_MAX_F32;
	int nbVertices, liveTriangles;

	bool b = input->size == 4;
	PxU32 warpMask = __ballot_sync(FULL_MASK, !b);

	if(b)
	{
		epaStart4(epaS, input, a, v, normal, distance);
		nbVertices = 4;
		//This is a tetrahedron so there are 4 triangles, hence liveTriangles = 0xf (1111)
		liveTriangles = 0xf;
	}
	else
	{
		epaStart123(epaS, input, vA, nbA, vB, nbB, a, v, normal, distance, warpMask);
		nbVertices = 3;
		//This is a double-sided triangle so there are 2 triangles, hence liveTriangles = 0x3 (11)
		liveTriangles = 0x3;
	}

	__syncwarp();

	// TODO: check here for vertex count is slighly conservative
	bool nondegenerate = true;
	while (nondegenerate && nbVertices < schlock::EPA_MAX_VERTS)
	{
		PxReal t = liveTriangles & 1 << threadIdx.x ? distance : FLT_MAX, e;
		//get the shortest distance triangle
		e = fminf(t, __shfl_xor_sync(FULL_MASK, t, 16));
		e = fminf(e, __shfl_xor_sync(FULL_MASK, e, 8));
		e = fminf(e, __shfl_xor_sync(FULL_MASK, e, 4));
		e = fminf(e, __shfl_xor_sync(FULL_MASK, e, 2));
		e = fminf(e, __shfl_xor_sync(FULL_MASK, e, 1));

		//try to find out which triangle index has the shortest distance
		int f = __ballot_sync(FULL_MASK, e == t);
		GPU_NP_CHECK(f); //Trigger assert if e == t failed on all lanes (indicates a nan)
		int triIndex = lowestSetIndex(f);

		//get the normal from the shortest distance triangle
		PxVec3 d(__shfl_sync(FULL_MASK, normal.x, triIndex), __shfl_sync(FULL_MASK, normal.y, triIndex), __shfl_sync(FULL_MASK, normal.z, triIndex));
		//get the support vertex
		PxU8 indexA = supportIndex(vA, nbA, d), indexB = supportIndex(vB, nbB, -d);
		PxVec3 supportVertex = vA[indexA] - vB[indexB];

		//calculate the distance from the origin to the support vertex
		PxReal thisFrameUb = supportVertex.dot(d);

		//found the MTD, epa terminate
		if (abs(thisFrameUb - e) <= releps)
		{
			upperBound = thisFrameUb, lowerBound = e;
			bestVerts.setInt(__shfl_sync(FULL_MASK, v.getInt(), triIndex));

			if (threadIdx.x == 0)
				epaS.bestDir = d;

			break;
		}

		// all remaining tris are worse than the best candidate, epa terminate
		if (e > upperBound)
			break;


		PxReal ub = PxMin(upperBound, thisFrameUb);

		if (ub < upperBound)
		{
			upperBound = ub, lowerBound = e;
			bestVerts.setInt(__shfl_sync(FULL_MASK, v.getInt(), triIndex));

			if (threadIdx.x == 0)
				epaS.bestDir = d;
		}

#if EPA_LOGGING
		if (threadIdx.x == 0)
			printf("[%f, %f], (%f, %f, %f)\n", lowerBound, upperBound, supportVertex.x, supportVertex.y, supportVertex.z);
#endif

		// save off the new vertex
		if (threadIdx.x == 0)
		{
			epaS.vertices[nbVertices] = supportVertex;
			epaS.indices[nbVertices] = PxU16(indexB) << 8 | PxU16(indexA);
		}

		__syncwarp();

		nondegenerate = replaceTriangle(epaS, ++nbVertices, triIndex, e, normal, distance, v, a, liveTriangles, releps);
	}

#if EPA_LOGGING
		if(!nondegenerate)
			printf("degenerate\n");
#endif

	output.degenerate = !nondegenerate;

	__syncwarp();

	schlock::GjkResult::Enum result;
	PxVec3 bestDir = epaS.bestDir;
	PxVec3 closestPoint = bestDir * lowerBound;

	if(threadIdx.x == 0)
	{
		output.direction = bestDir;
		output.lowerBound = lowerBound;
		output.upperBound = upperBound;

		int v0 = bestVerts.get(0), v1 = bestVerts.get(1), v2 = bestVerts.get(2);
		PxReal u, v, w;
		schlock::barycentrics(epaS.vertices[v0], epaS.vertices[v1], epaS.vertices[v2], u, v, w, closestPoint);
		PxVec3 cpA = vA[schlock::aIndex(epaS.indices[v0])] * u + vA[schlock::aIndex(epaS.indices[v1])] * v + vA[schlock::aIndex(epaS.indices[v2])] * w;
		output.closestPointA = cpA;
		output.closestPointB = cpA - closestPoint;
	}
	result = lowerBound < 0 ? schlock::GjkResult::eCLOSE : schlock::GjkResult::eOVERLAP;

	return result;
}
} // namespace squawk

} // namespace physx

#endif

// Copyright (c) 2001-2004 NovodeX AG. All rights reserved.
// Copyright (c) 2004-2008 AGEIA Technologies, Inc. All rights reserved.
// SPDX-FileCopyrightText: Copyright (c) 2008-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef __GJK_CUH__
#define __GJK_CUH__

namespace physx
{

namespace squawk
{
// The simplex is stored as 4 vectors and then 4 normals (so we can shoot it straight into EPA)
// The 3 component of each vector contains 3f80 in the high bits (so it's a valid float) and
// the vertex indices for a and b in the lower bits

typedef PxReal Simplex;

static __device__ inline Simplex setNormal0(Simplex z, PxReal nV, PxU32 index)
{
	PxReal tmpV = splatV(index, nV);
	return threadIdx.x < 16 || threadIdx.x > 19 ? z : tmpV;
}

static __device__ inline Simplex setNormal123(Simplex z, PxReal nV)
{
	PxReal tmpV = __shfl_up_sync(FULL_MASK, nV, 20);
	return threadIdx.x < 20 ? z : tmpV;
}

static __device__ inline Simplex loadSimplex(const schlock::GjkCachedData& data)
{
	return reinterpret_cast<const PxReal*>(&data)[threadIdx.x];
}

static __device__ inline void storeSimplex(schlock::GjkCachedData& data, Simplex z)
{
	reinterpret_cast<PxReal*>(&data)[threadIdx.x] = z;
}

static __device__ inline void storeSimplex(volatile schlock::GjkCachedData& data, Simplex z)
{
	reinterpret_cast<volatile PxReal*>(&data)[threadIdx.x] = z;
}

// suffixes for vector warps:
// * V = 3 vectors: (aX,aY,aZ,_,bX,bY,bZ,_,cX,cY,cZ,_)
// * D = dense scalar: (a,a,a,a,b,b,b,b,c,c,c,c)
// * S = sparse scalar: (a,_,_,_,b,_,_,_,c,_,_,_)
// * I = index vector for shuffle

#define DEGENERATE 0xffff

static __device__ inline PxU32 closestPointTri(Simplex &zV, PxReal cV, PxReal& closestPoint)
{
	PxU32 tI = threadIdx.x, cp0I = (tI&28)|next3(tI&3), cp1I = (tI&28)|prev3(tI&3);		// cross product indices

	closestPoint = 0;
	PxReal eV = cV - zV, lD = splatX(dot(eV, eV));										// edges to apex & squared lengths
	if(anyX<2>(lD==0))
		return DEGENERATE;																// FAIL: apex coincident with a base vert

	PxReal tmp = __shfl_xor_sync(FULL_MASK, eV, 4), shortEdgeV = __shfl_xor_sync(FULL_MASK, lD, 4) < lD ? tmp : eV;			// shorter edge to apex
	PxReal mV = cross(__shfl_xor_sync(FULL_MASK, zV,4)-zV, shortEdgeV, cp0I, cp1I);						// (normal, -normal)

	// lots of dot products here with cV, (plus proj below). TODO: optimize
	PxReal pS = dot(eV,cV);																// projection of origin onto edges
	bool uS = dot(cross(mV, eV, cp0I, cp1I), cV) <= 0;									// origin outside edge planes

	if(!anyX<2>(uS))																	// REGION: face - inside both edge planes
	{
		PxReal projD = __shfl_sync(FULL_MASK, dot(mV,cV),0);												// project onto face
		closestPoint = splatV(0, mV*projD/splatX(dot(mV,mV)));
		bool flip = projD<0;															// ensure the normal winding makes it outward-facing in the tet
		mV = (tI & 3) == 3 ? projD : mV;
		zV = setNormal0(zV, mV, flip);
		return flip ? 0x0123u : 0x1023u;
	}

	PxU32 eid = __ffs(packXFlags<2>(uS & bool(pS>0)));									// REGION: edges - outside the edge plane, and projects onto the edge
	if(eid--)																			// __ffs returns 0 for no bit, 1 for lsb
	{
		closestPoint = splatV(eid, cV - eV*splatX(pS)/lD);
		return eid<<8 | 0x22;
	}

	closestPoint = cV;																	// REGION: apex - not an edge, not a face
	return 0x0021;
}


static __device__ inline PxU32 closestPointTet(Simplex &zV, PxReal dV, PxReal& closestPoint)
{
	closestPoint = 0;

	PxReal n0V = splatV(4, zV);
	if(__shfl_sync(FULL_MASK, dot(n0V, dV) > dot(n0V, zV), 0))											// FAIL: new vertex is beyond old face (TODO: set up normal in UT)
		return DEGENERATE;

	PxU32 tI = threadIdx.x, cp0I = (tI&28)|next3(tI&3), cp1I = (tI&28)|prev3(tI&3);		// cross product indices
	PxU32 vNextI = tI>7 ? tI-8 : tI+4, vPrevI = tI<4 ? tI + 8 : tI-4;					// block rotate indices
	PxReal eV = dV - zV, lD = splatX(dot(eV, eV));										// edges to apex and squared lengths
	if(anyX<3>(lD==0))
		return DEGENERATE;																// FAIL: apex is coincident with a base vert

	PxReal tmp = __shfl_sync(FULL_MASK, eV, vNextI), shortEdgeV = __shfl_sync(FULL_MASK, lD, vNextI) < lD ? tmp : eV;		// shorter of edges to apex for each tri
	PxReal nV = cross(shortEdgeV, __shfl_sync(FULL_MASK, zV, vNextI) - zV, cp0I, cp1I);					// face normals, outward facing
	PxReal aS = dot(nV, dV);															// distances along normals of apex
	PxReal gS = dot(__shfl_sync(FULL_MASK, zV, vPrevI), nV);												// distances along normal of vert not in this triangle

	if(anyX<3>(aS<=gS))
		return DEGENERATE;																// FAIL: simplex inversion/collapse

	PxReal tmpV = __shfl_up_sync(FULL_MASK, aS, 3);
	nV = (tI & 3) == 3 ? tmpV : nV;

	if(allX<3>(aS>0))																	// REGION: inside simplex
	{
		zV = setNormal123(zV, nV);														// set normals for EPA start
		return 0x2104;
	}

	const PxReal pS = dot(eV,dV);														// projection of origin onto edges
	if(allX<3>(pS<=0))																	// REGION: apex
	{
		closestPoint = dV;
		return 0x2131;
	}

	bool uiS = dot(cross(eV, nV, cp0I, cp1I), dV) <= 0;									// origin outside one edge plane of each triangle
	bool ujS = dot(cross(__shfl_sync(FULL_MASK, eV,vNextI), nV, cp0I, cp1I), dV) > 0;						// origin outside other edge plane of each triangle
	int eid = __ffs(packXFlags<3>(uiS & __shfl_sync(FULL_MASK, ujS, vPrevI) & bool(pS>=0)));				// REGION: edges - origin outside edge planes and projects onto edge
	if(eid--)																			// __ffs returns 0 for no bit, 1 for lsb
	{
		closestPoint = splatV(eid, dV - eV * splatX(pS/lD));
		return eid<<8| 0x32;
	}

	PxU32 fid = __ffs(packXFlags<3>(!uiS & !ujS & bool(aS<=0)));						// REGION: faces - origin outside face plane, inside edge planes
	if(fid--)																			// __ffs returns 0 for no bit, 1 for lsb
	{
		closestPoint = splatV(fid, nV * splatX(aS / dot(nV, nV)));
		zV = setNormal0(zV, -nV, fid);													// flip normal to be outward facing for next iteration
		return next3(fid)<<12 | fid<<8| 0x33;
	}

	return DEGENERATE;																	// FAIL: couldn't find voronoi region. Chort vozmi!
}

static __device__ inline PxU8 supportIndex(const PxVec3* v, PxU8 nbVertices, PxReal dir)
{
	int bestN = 0xFFFFFFFF; //Ininitialize to a very illegal value to force a crash when we get nans
	PxReal bestDist = -FLT_MAX;
	PxReal x = __shfl_sync(FULL_MASK, dir, 0), y = __shfl_sync(FULL_MASK, dir, 1), z = __shfl_sync(FULL_MASK, dir, 2);
	for(PxU32 i=0;i<nbVertices;i+=32)
	{
		int n = i+threadIdx.x;
		PxReal d = n < nbVertices ? x*v[n].x + y*v[n].y + z*v[n].z : -FLT_MAX, e;

		e = fmaxf(d, __shfl_xor_sync(FULL_MASK, d, 16));
		e = fmaxf(e, __shfl_xor_sync(FULL_MASK, e, 8));
		e = fmaxf(e, __shfl_xor_sync(FULL_MASK, e, 4));
		e = fmaxf(e, __shfl_xor_sync(FULL_MASK, e, 2));
		e = fmaxf(e, __shfl_xor_sync(FULL_MASK, e, 1));

		unsigned mask_bestDist = __ballot_sync(FULL_MASK, e > bestDist);
		if(e>bestDist)
			bestDist = e, bestN = i+__ffs(__ballot_sync(mask_bestDist, e == d))-1;
	}

	GPU_NP_CHECK(bestN < nbVertices);
	return PxU8(bestN);
}

static __device__ inline PxU8 supportIndex(const PxVec3& v0, const PxVec3& v1, PxU8 nbVertices, PxReal dir)
{
	int bestN = 0xFFFFFFFF;
	//PxReal bestDist = -FLT_MAX;
	PxReal x = __shfl_sync(FULL_MASK, dir, 0), y = __shfl_sync(FULL_MASK, dir, 1), z = __shfl_sync(FULL_MASK, dir, 2);

	PxReal d = -FLT_MAX;
	PxU32 n = threadIdx.x;
	if (threadIdx.x < nbVertices)
	{
		d = x*v0.x + y*v0.y + z*v0.z;

		if ((threadIdx.x + 32) < nbVertices)
		{
			PxReal d1 = x*v0.x + y*v0.y + z*v0.z;
			if (d1 > d)
			{
				d = d1;
				n += 32;
			}
		}
	}

	PxReal e = fmaxf(d, __shfl_xor_sync(FULL_MASK, d, 16));
	e = fmaxf(e, __shfl_xor_sync(FULL_MASK, e, 8));
	e = fmaxf(e, __shfl_xor_sync(FULL_MASK, e, 4));
	e = fmaxf(e, __shfl_xor_sync(FULL_MASK, e, 2));
	e = fmaxf(e, __shfl_xor_sync(FULL_MASK, e, 1));

	//bestDist = e;
	bestN = __ffs(__ballot_sync(FULL_MASK, e == d)) - 1;
	bestN += __popc(__ballot_sync(FULL_MASK, n >= 32) & (1 << bestN)) * 32;

	GPU_NP_CHECK(bestN < nbVertices);

	return PxU8(bestN);
}


static __device__ inline schlock::GjkResult::Enum gjk(
										const PxVec3* vA, PxU32 nbA,
										const PxVec3* vB, PxU32 nbB,
										const PxVec3& initialUnitDir,
										PxReal minSep,
										PxReal convergenceRatio,
										schlock::GjkOutput& output, schlock::GjkCachedData& cacheData)
{
	int cacheMask = (1<< cacheData.size)-1;
	__syncwarp(); //cacheData.size is written further below
	Simplex simplexV = loadSimplex(cacheData);

restart:
	bool warmStarted = false;
	Simplex backupSimplexV;
	int size = 0, backupSize = 0;
	PxReal closestPointV, backupClosestV;

	if (cacheMask)
		closestPointV = splatV(0, simplexV);
	else
	{
		if ((threadIdx.x & 3) == 0)
			closestPointV = initialUnitDir.x;
		else if ((threadIdx.x & 3) == 1)
			closestPointV = initialUnitDir.y;
		else if ((threadIdx.x & 3) == 2)
			closestPointV = initialUnitDir.z;
		else
			closestPointV = PX_MAX_F32;
	}

	PxReal upperBound = PX_MAX_F32, lowerBound = -PX_MAX_F32, cm = magnitude(closestPointV); // if we hit a lower bound dir with distance exactly zero, -FLT_MIN means we keep it
	bool degenerate;

	do
	{
		PxReal dV = closestPointV/cm, supportVertexV;

#if GJK_LOGGING
		PxReal dx = dV, dy = shfl(dV, 1), dz = shfl(dV,2);
		PxReal cpx = closestPointV, cpy = shfl(closestPointV, 1), cpz = shfl(closestPointV,2);
		if(threadIdx.x == 0)
			printf("[%f, %f], dir (%f, %f, %f), closestPoint (%f, %f, %f), cm %f \n", lowerBound, upperBound, dx, dy, dz, cpx, cpy, cpz, cm);
#endif

		// find a warm start entry that can be added to the simplex because it's closer than the closest feature
		int candidates = 0;

		if (cacheMask)
			candidates = __ballot_sync(FULL_MASK, __shfl_sync(FULL_MASK, dot(dV, simplexV)<upperBound, threadIdx.x<<2)) & cacheMask;

		if(candidates)
		{
			warmStarted = true;
			cacheMask &= ~(lowestSetBit(candidates)*2-1);
			supportVertexV = splatV(lowestSetIndex(candidates), simplexV);
			if(threadIdx.x>>2 == size)
				simplexV = supportVertexV;
		}
		else				// if not, stop warm-starting and start using support vertices
		{
			cacheMask = 0;
			PxU8 indexA = supportIndex(vA, nbA, -dV), indexB = supportIndex(vB, nbB, dV);
			supportVertexV = loadV3Unsafe(vA+indexA) - loadV3Unsafe(vB+indexB);
			if(threadIdx.x>>2 == size)
				simplexV = (threadIdx.x&3) < 3 ? supportVertexV : __int_as_float(schlock::encodeIndices(indexA, indexB));
			PxReal lb = splatX(dot(dV,supportVertexV));
			if(lb>lowerBound)
			{
				storeV3(&output.direction, 0, -dV);
				lowerBound = lb;
				if(lowerBound > minSep || lowerBound>=upperBound*convergenceRatio)
					break;
			}
		}

#if GJK_LOGGING
		PxReal svx = supportVertexV, svy = shfl(supportVertexV,1), svz = shfl(supportVertexV,2);
		if(threadIdx.x == 0)
			printf("support (%f, %f, %f)\n", svx, svy, svz);
#endif

		PxU32 f = 0;											// f is feature index: (i2<<12 | i1<<8 | i0<<4 | size)
		PxReal pointV;
		switch(size)
		{
		case 0:		f = 0x2101, pointV = supportVertexV; break;	// have to be a bit for 0,1 cases not to trash the cached verts, which also live in the simplex
		case 1:
		{
			PxReal aV = splatV(0,simplexV), lV = supportVertexV - aV, m = splatX(dot(lV,lV)), x = -splatX(dot(aV,lV));
			if(m==0)		f = DEGENERATE, pointV = 0;							// FAIL: coincident vertex
			else if(x<m)	f = 0x2102,		pointV = aV + lV * x/m;				// REGION: line segment
			else 			f = 0x2111,		pointV = supportVertexV;			// REGION: new vertex
			break;
		}
		case 2:	f = closestPointTri(simplexV, supportVertexV, pointV);	break; // no worries for 2 & 3: only (0,1,2) are affected by permute, and they're already used
		case 3:	f = closestPointTet(simplexV, supportVertexV, pointV);	break;
		}

		cm = magnitude(pointV);
		int newSize = f&0xf;
		degenerate = f == DEGENERATE || (cm>=upperBound && newSize<=size);
		if(degenerate)
		{
			if(!warmStarted)	// hit a problem warm-starting: just restart cold
				break;
			cacheMask = 0;
			goto restart;
		}

		// local minima happen. Continue so long as either the simplex size increases or the upper bound decreases,
		backupSize = cm >= upperBound ? size : 0;
		if(backupSize)
			backupSimplexV = simplexV, backupClosestV = closestPointV;	// save the best simplex and size before updating state
		else
			upperBound = cm;

#if GJK_LOGGING
		PxReal px = pointV, py = shfl(pointV,1), pz = shfl(pointV,2);
		if(threadIdx.x == 0)
			printf("%04x (%f, %f, %f) %f\n", f, px, py, pz, cm);
#endif

		closestPointV = pointV, size = newSize;										// update closest point and simplex
		PxU32 sourceIndex = f>>((threadIdx.x&28)+4)&0xf;							// +4 because size is the first nibble
		PxReal tmp = __shfl_sync(FULL_MASK, simplexV, sourceIndex<<2 | (threadIdx.x&3));
		if(threadIdx.x<12)
			simplexV = tmp;
	}
	while(upperBound > 0 && lowerBound < upperBound * convergenceRatio);

	if(backupSize>0)
		simplexV = backupSimplexV, closestPointV = backupClosestV, size = backupSize;

	__syncwarp();
	storeSimplex(cacheData, simplexV);

	if(threadIdx.x == 0)
	{
		output.degenerate = degenerate;
		output.lowerBound = lowerBound;
		output.upperBound = upperBound;
		cacheData.size = size;
	}

	__syncwarp();

	if(upperBound == 0)
		return schlock::GjkResult::eOVERLAP;

	if(lowerBound>minSep)
		return schlock::GjkResult::eDISTANT;

	PxU8 a0 = cacheData.vertices[0].a(), a1 = cacheData.vertices[1].a(), a2 = cacheData.vertices[2].a();
	PxU8 b0 = cacheData.vertices[0].b(), b1 = cacheData.vertices[1].b(), b2 = cacheData.vertices[2].b();

	PxVec3 closestPoint(closestPointV, __shfl_sync(FULL_MASK, closestPointV,1), __shfl_sync(FULL_MASK, closestPointV,2));

	if(threadIdx.x == 0)
	{
		PxVec3 closestPointDir = -closestPoint.getNormalized();
		output.closestPointDir = closestPointDir;

		PxReal p, q, r;
		if(size == 3)
		{
			schlock::barycentrics(cacheData.vertices[0].getV(), cacheData.vertices[1].getV(), cacheData.vertices[2].getV(), p, q, r, closestPoint);
			output.closestPointA = p * vA[a0] + q * vA[a1] + r * vA[a2];
		}
		else if(size == 2)
		{
			PxVec3 d = cacheData.vertices[1].getV()- cacheData.vertices[0].getV();
			PxReal q = -cacheData.vertices[0].getV().dot(d)/d.magnitudeSquared();
			output.closestPointA = (1-q) * vA[a0] + q * vA[a1];
		}
		else
			output.closestPointA = vA[a0];

		output.closestPointB = output.closestPointA + upperBound * closestPointDir;
	}

	return schlock::GjkResult::eCLOSE;
}
} // namespace squawk

} // namespace physx

#undef DEGENERATE_SIZE_FLAG

#endif

#undef GPU_NP_CHECK
#undef PX_CUDA_CALLABLE
#undef PX_FORCE_INLINE
