#ifndef BT_TRANSLATION_COUPLING_KERNEL_H
#define BT_TRANSLATION_COUPLING_KERNEL_H
#include "LinearMath/btVector3.h"
#if defined(BT_USE_DOUBLE_PRECISION) && defined(__AVX2__)
#include <immintrin.h>
#endif
// Reduce the three mode dot products in parallel, retaining node order.
SIMD_FORCE_INLINE btVector3 btComputeTranslationCoupling(int count, const btVector3* values,
														 const btVector3* az0, const btVector3* az1, const btVector3* az2)
{
#if defined(BT_USE_DOUBLE_PRECISION) && defined(__AVX2__)
	const __m256d zero = _mm256_setzero_pd();
	__m256d sum = zero;
	for (int n = 0; n < count; ++n)
	{
		const __m256d a0 = _mm256_loadu_pd(&az0[n][0]);
		const __m256d a1 = _mm256_loadu_pd(&az1[n][0]);
		const __m256d a2 = _mm256_loadu_pd(&az2[n][0]);
		// Transpose: SIMD lanes represent modes, not spatial components.
		const __m256d lo01 = _mm256_unpacklo_pd(a0, a1);
		const __m256d hi01 = _mm256_unpackhi_pd(a0, a1);
		const __m256d lo2 = _mm256_unpacklo_pd(a2, zero);
		const __m256d hi2 = _mm256_unpackhi_pd(a2, zero);
		const __m256d x = _mm256_permute2f128_pd(lo01, lo2, 0x20);
		const __m256d y = _mm256_permute2f128_pd(hi01, hi2, 0x20);
		const __m256d z = _mm256_permute2f128_pd(lo01, lo2, 0x31);
		const __m256d xy = _mm256_add_pd(_mm256_mul_pd(x, _mm256_set1_pd(values[n][0])),
										 _mm256_mul_pd(y, _mm256_set1_pd(values[n][1])));
		const __m256d dot = _mm256_add_pd(xy, _mm256_mul_pd(z, _mm256_set1_pd(values[n][2])));
		sum = _mm256_add_pd(sum, dot);
	}
	btVector3 result;
	_mm256_storeu_pd(&result[0], _mm256_blend_pd(sum, zero, 8));
	return result;
#else
	btVector3 result(0, 0, 0);
	for (int n = 0; n < count; ++n)
	{
		result[0] += az0[n].dot(values[n]);
		result[1] += az1[n].dot(values[n]);
		result[2] += az2[n].dot(values[n]);
	}
	return result;
#endif
}
#endif
