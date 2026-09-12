#ifndef BT_TRANSLATION_INPUT_KERNEL_H
#define BT_TRANSLATION_INPUT_KERNEL_H
#include "LinearMath/btVector3.h"
#if defined(BT_USE_DOUBLE_PRECISION) && defined(__AVX2__)
#include <immintrin.h>
#endif
// Subtract the three coarse-mode contributions in their original order.
// Work and mode arrays are distinct; preserve each work vector's fourth lane.
SIMD_FORCE_INLINE void btApplyTranslationInput(int count, btVector3* work,
											   const btVector3* az0, const btVector3* az1, const btVector3* az2, const btVector3& coarse)
{
#if defined(BT_USE_DOUBLE_PRECISION) && defined(__AVX2__)
	const __m256d c0 = _mm256_set1_pd(coarse[0]);
	const __m256d c1 = _mm256_set1_pd(coarse[1]);
	const __m256d c2 = _mm256_set1_pd(coarse[2]);
	const __m256d zero = _mm256_setzero_pd();
	for (int n = 0; n < count; ++n)
	{
		const __m256d old = _mm256_loadu_pd(&work[n][0]);
		__m256d value = _mm256_blend_pd(old, zero, 8);
		const __m256d a0 = _mm256_blend_pd(_mm256_loadu_pd(&az0[n][0]), zero, 8);
		const __m256d a1 = _mm256_blend_pd(_mm256_loadu_pd(&az1[n][0]), zero, 8);
		const __m256d a2 = _mm256_blend_pd(_mm256_loadu_pd(&az2[n][0]), zero, 8);
		value = _mm256_sub_pd(value, _mm256_mul_pd(a0, c0));
		value = _mm256_sub_pd(value, _mm256_mul_pd(a1, c1));
		value = _mm256_sub_pd(value, _mm256_mul_pd(a2, c2));
		_mm256_storeu_pd(&work[n][0], _mm256_blend_pd(value, old, 8));
	}
#else
	const btScalar c0 = coarse[0], c1 = coarse[1], c2 = coarse[2];
	for (int n = 0; n < count; ++n)
	{
		work[n] -= az0[n] * c0;
		work[n] -= az1[n] * c1;
		work[n] -= az2[n] * c2;
	}
#endif
}
#endif
