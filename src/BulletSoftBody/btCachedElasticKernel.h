#ifndef BT_CACHED_ELASTIC_KERNEL_H
#define BT_CACHED_ELASTIC_KERNEL_H
#include "LinearMath/btMatrix3x3.h"
#if defined(BT_USE_DOUBLE_PRECISION) && defined(__AVX2__)
#include <immintrin.h>
// Apply one cached isotropic tetrahedron tangent. Input displacements and
// output forces are separate; output fourth components retain their values.
SIMD_FORCE_INLINE void btApplyCachedElasticTetra(const btMatrix3x3& gradients,
												 btScalar mu, btScalar lambda, const btVector3& x0, const btVector3& x1,
												 const btVector3& x2, const btVector3& x3,
												 btVector3& f0, btVector3& f1, btVector3& f2, btVector3& f3)
{
	const __m256d zero = _mm256_setzero_pd();
	const __m256d base = _mm256_blend_pd(_mm256_loadu_pd(&x0[0]), zero, 8);
	const __m256d d1 = _mm256_sub_pd(_mm256_blend_pd(_mm256_loadu_pd(&x1[0]), zero, 8), base);
	const __m256d d2 = _mm256_sub_pd(_mm256_blend_pd(_mm256_loadu_pd(&x2[0]), zero, 8), base);
	const __m256d d3 = _mm256_sub_pd(_mm256_blend_pd(_mm256_loadu_pd(&x3[0]), zero, 8), base);
	// Columns of G = sum_j (x_j-x_0) grad(N_j)^T; all XYZ lanes work together.
	const auto column = [&](int axis)
	{
		return _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(d1, _mm256_set1_pd(gradients[0][axis])),
										   _mm256_mul_pd(d2, _mm256_set1_pd(gradients[1][axis]))),
							 _mm256_mul_pd(d3, _mm256_set1_pd(gradients[2][axis])));
	};
	const __m256d g0 = column(0), g1 = column(1), g2 = column(2);
	// Transpose the 3x3 G, embedded in a zero-padded 4x4 register matrix.
	const __m256d lo01 = _mm256_unpacklo_pd(g0, g1), hi01 = _mm256_unpackhi_pd(g0, g1);
	const __m256d lo2 = _mm256_unpacklo_pd(g2, zero), hi2 = _mm256_unpackhi_pd(g2, zero);
	const __m256d row0 = _mm256_permute2f128_pd(lo01, lo2, 0x20);
	const __m256d row1 = _mm256_permute2f128_pd(hi01, hi2, 0x20);
	const __m256d row2 = _mm256_permute2f128_pd(lo01, lo2, 0x31);
	const btScalar trace = _mm256_cvtsd_f64(g0) +
						   _mm_cvtsd_f64(_mm_unpackhi_pd(_mm256_castpd256_pd128(g1), _mm256_castpd256_pd128(g1))) +
						   _mm_cvtsd_f64(_mm256_extractf128_pd(g2, 1));
	const btScalar diagonal = lambda * trace;
	const __m256d weight = _mm256_set1_pd(mu);
	const __m256d s0 = _mm256_add_pd(_mm256_mul_pd(_mm256_add_pd(g0, row0), weight), _mm256_setr_pd(diagonal, 0, 0, 0));
	const __m256d s1 = _mm256_add_pd(_mm256_mul_pd(_mm256_add_pd(g1, row1), weight), _mm256_setr_pd(0, diagonal, 0, 0));
	const __m256d s2 = _mm256_add_pd(_mm256_mul_pd(_mm256_add_pd(g2, row2), weight), _mm256_setr_pd(0, 0, diagonal, 0));
	const auto force = [&](int node)
	{
		return _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(s0, _mm256_set1_pd(gradients[node][0])),
										   _mm256_mul_pd(s1, _mm256_set1_pd(gradients[node][1]))),
							 _mm256_mul_pd(s2, _mm256_set1_pd(gradients[node][2])));
	};
	const __m256d c1 = force(0), c2 = force(1), c3 = force(2);
	const auto accumulate = [](btVector3& dst, __m256d contribution)
	{
		const __m256d old = _mm256_loadu_pd(&dst[0]);
		_mm256_storeu_pd(&dst[0], _mm256_blend_pd(_mm256_add_pd(old, contribution), old, 8));
	};
	accumulate(f0, _mm256_sub_pd(zero, _mm256_add_pd(_mm256_add_pd(c1, c2), c3)));
	accumulate(f1, c1);
	accumulate(f2, c2);
	accumulate(f3, c3);
}
#endif
#endif
