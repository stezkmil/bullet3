#ifndef BT_MASS_KERNEL_H
#define BT_MASS_KERNEL_H
#include "LinearMath/btVector3.h"
#if defined(BT_USE_DOUBLE_PRECISION) && defined(__AVX2__)
#include <immintrin.h>
#endif
// Node data is read on every application, preserving changes to inverse mass
// and frozen state. Input may equal output; vector padding is cleared.
template <class Node>
SIMD_FORCE_INLINE void btApplyNodeMass(int count, const Node* nodes,
									   const btVector3* input, btVector3* output)
{
	for (int n = 0; n < count; ++n)
	{
		if (nodes[n].m_frozen > 0)
		{
			output[n].setZero();
			continue;
		}
		const btScalar inverseMass = nodes[n].m_im;
		btFullAssert(inverseMass != btScalar(0));
		const btScalar mass = btScalar(1) / inverseMass;
#if defined(BT_USE_DOUBLE_PRECISION) && defined(__AVX2__)
		const __m256d zero = _mm256_setzero_pd();
		const __m256d value = _mm256_blend_pd(_mm256_loadu_pd(&input[n][0]), zero, 8);
		const __m256d product = _mm256_mul_pd(value, _mm256_set1_pd(mass));
		_mm256_storeu_pd(&output[n][0], _mm256_blend_pd(product, zero, 8));
#else
		output[n] = input[n] * mass;
#endif
	}
}
#endif
