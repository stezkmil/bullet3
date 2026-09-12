/*
 Written by Xuchen Han <xuchenhan2015@u.northwestern.edu>

 Bullet Continuous Collision Detection and Physics Library
 Copyright (c) 2019 Google Inc. http://bulletphysics.org
 This software is provided 'as-is', without any express or implied warranty.
 In no event will the authors be held liable for any damages arising from the use of this software.
 Permission is granted to anyone to use this software for any purpose,
 including commercial applications, and to alter it and redistribute it freely,
 subject to the following restrictions:
 1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
 2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
 3. This notice may not be removed or altered from any source distribution.
 */

#ifndef BT_KRYLOV_SOLVER_H
#define BT_KRYLOV_SOLVER_H

#include <iostream>
#include <cmath>
#include <limits>
#include <LinearMath/btAlignedObjectArray.h>
#include <LinearMath/btVector3.h>
#include <LinearMath/btScalar.h>
#include "LinearMath/btQuickprof.h"
#include "btDeformableOptimizationConfig.h"
#if defined(BT_USE_DOUBLE_PRECISION) && defined(__AVX2__)
#include <immintrin.h>
#endif

template <class MatrixX>
class btKrylovSolver
{
	typedef btAlignedObjectArray<btVector3> TVStack;

public:
	int m_maxIterations;
	btScalar m_tolerance;
	btKrylovSolver(int maxIterations, btScalar tolerance)
		: m_maxIterations(maxIterations), m_tolerance(tolerance)
	{
	}

	virtual ~btKrylovSolver() {}

	virtual int solve(MatrixX& A, TVStack& x, const TVStack& b, bool verbose = false) = 0;

	virtual void reinitialize(const TVStack& b) = 0;

	virtual SIMD_FORCE_INLINE TVStack sub(const TVStack& a, const TVStack& b)
	{
		// c = a-b
		btAssert(a.size() == b.size());
		TVStack c;
		c.resize(a.size());
		for (int i = 0; i < a.size(); ++i)
		{
			c[i] = a[i] - b[i];
		}
		return c;
	}

	// result = a - b. Reuses capacity; result may alias either input.
	SIMD_FORCE_INLINE void subtractInto(const TVStack& a, const TVStack& b, TVStack& result)
	{
		btAssert(a.size() == b.size());
		result.resize(a.size());
		for (int i = 0; i < a.size(); ++i)
			result[i] = a[i] - b[i];
	}

	virtual SIMD_FORCE_INLINE btScalar squaredNorm(const TVStack& a)
	{
		return dot(a, a);
	}

	virtual SIMD_FORCE_INLINE btScalar norm(const TVStack& a)
	{
		btScalar ret = 0;
		for (int i = 0; i < a.size(); ++i)
		{
			for (int d = 0; d < 3; ++d)
			{
				ret = btMax(ret, btFabs(a[i][d]));
			}
		}
		return ret;
	}

	virtual SIMD_FORCE_INLINE btScalar dot(const TVStack& a, const TVStack& b)
	{
		btScalar ans(0);
		for (int i = 0; i < a.size(); ++i)
			ans += a[i].dot(b[i]);
		return ans;
	}

	virtual SIMD_FORCE_INLINE void multAndAddTo(btScalar s, const TVStack& a, TVStack& result)
	{
#if (BT_DEFORMABLE_OPTIMIZATION_MASK & 32)

		//        result += s*a
		btAssert(a.size() == result.size());
		const int n = a.size();
		if (n == 0) return;
		const btVector3* const src = &a[0];
		btVector3* const dst = &result[0];
#if defined(BT_USE_DOUBLE_PRECISION) && defined(__AVX2__)
		const __m256d scale = _mm256_set1_pd(s);
		const __m256d zero = _mm256_setzero_pd();
		for (int i = 0; i < n; ++i)
		{
			const __m256d old = _mm256_loadu_pd(&dst[i][0]);
			const __m256d input = _mm256_blend_pd(_mm256_loadu_pd(&src[i][0]), zero, 8);
			const __m256d value = _mm256_add_pd(_mm256_blend_pd(old, zero, 8),
				_mm256_mul_pd(scale, input));
			// Unlike assignment arithmetic, operator+= preserves the fourth lane.
			_mm256_storeu_pd(&dst[i][0], _mm256_blend_pd(value, old, 8));
		}
#else
		for (int i = 0; i < n; ++i) dst[i] += s * src[i];
#endif

#else

		//        result += s*a
		btAssert(a.size() == result.size());
		for (int i = 0; i < a.size(); ++i)
			result[i] += s * a[i];

#endif
	}

	// result = s * result + b, reusing the destination storage.
	// Both arrays must already have equal sizes; existing result values are inputs.
	SIMD_FORCE_INLINE void scaleAndAddInPlace(btScalar s, const TVStack& b, TVStack& result)
	{
#if (BT_DEFORMABLE_OPTIMIZATION_MASK & 1)

		btAssert(b.size() == result.size());
		const int n = result.size();
		if (n == 0) return;
		btVector3* const dst = &result[0];
		const btVector3* const src = &b[0];
#if defined(BT_USE_DOUBLE_PRECISION) && defined(__AVX2__)
		const __m256d scale = _mm256_set1_pd(s);
		const __m256d zero = _mm256_setzero_pd();
		for (int i = 0; i < n; ++i)
		{
			// btVector3 guarantees only 16-byte alignment. Ignore its fourth
			// component on input and clear it on output, like vector arithmetic.
			const __m256d x = _mm256_blend_pd(_mm256_loadu_pd(&dst[i][0]), zero, 8);
			const __m256d y = _mm256_blend_pd(_mm256_loadu_pd(&src[i][0]), zero, 8);
			const __m256d value = _mm256_add_pd(_mm256_mul_pd(scale, x), y);
			_mm256_storeu_pd(&dst[i][0], _mm256_blend_pd(value, zero, 8));
		}
#else
		for (int i = 0; i < n; ++i)
			dst[i] = s * dst[i] + src[i];
#endif

#else

		btAssert(b.size() == result.size());
		for (int i = 0; i < result.size(); ++i)
			result[i] = s * result[i] + b[i];

#endif
	}

	virtual SIMD_FORCE_INLINE TVStack multAndAdd(btScalar s, const TVStack& a, const TVStack& b)
	{
		// result = a*s + b
		TVStack result;
		result.resize(a.size());
		for (int i = 0; i < a.size(); ++i)
			result[i] = s * a[i] + b[i];
		return result;
	}

	virtual SIMD_FORCE_INLINE void setTolerance(btScalar tolerance)
	{
		m_tolerance = tolerance;
	}
};
#endif /* BT_KRYLOV_SOLVER_H */
