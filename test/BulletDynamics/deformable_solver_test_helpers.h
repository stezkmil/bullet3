/* Independent numerical helpers for deformable solver regression tests. */
#ifndef BT_DEFORMABLE_SOLVER_TEST_HELPERS_H
#define BT_DEFORMABLE_SOLVER_TEST_HELPERS_H

#include "BulletSoftBody/btSoftBody.h"
#include <cmath>

namespace btDeformableTest
{
typedef btAlignedObjectArray<btVector3> Vectors;

inline btScalar dot(const Vectors& a, const Vectors& b)
{
	btScalar result = 0;
	for (int n = 0; n < a.size(); ++n) result += a[n].dot(b[n]);
	return result;
}

inline btVector3 sumOnBody(const btSoftBody& body, const Vectors& values, bool freeOnly)
{
	btVector3 sum(0, 0, 0);
	for (int n = 0; n < body.m_nodes.size(); ++n)
	{
		const btSoftBody::Node& node = body.m_nodes[n];
		if (!freeOnly || (node.m_frozen <= 0 && node.m_im > 0)) sum += values[node.index];
	}
	return sum;
}

struct OperatorProbe
{
	btScalar symmetryError, linearityError, repeatabilityError, zeroNorm;
	btScalar curvatureP, curvatureQ;
};

template<class Matrix>
OperatorProbe probeOperator(Matrix& A, int size)
{
	Vectors p, q, ap, aq, work, combined;
	p.resize(size); q.resize(size); ap.resize(size); aq.resize(size); work.resize(size); combined.resize(size);
	unsigned int seed = 137;
	for (int n = 0; n < size; ++n)
		for (int d = 0; d < 3; ++d)
		{
			seed = 1664525u * seed + 1013904223u;
			p[n][d] = btScalar((seed >> 8) & 65535u) / btScalar(32768) - btScalar(1);
			seed = 1664525u * seed + 1013904223u;
			q[n][d] = btScalar((seed >> 8) & 65535u) / btScalar(32768) - btScalar(1);
		}
	A.multiply(p, ap); A.multiply(q, aq);
	const btScalar pp = dot(p, p), qq = dot(q, q);
	const btScalar scale = btMax(btScalar(1e-30), btSqrt(pp * dot(aq, aq)) + btSqrt(qq * dot(ap, ap)));
	OperatorProbe result;
	result.symmetryError = btFabs(dot(p, aq) - dot(q, ap)) / scale;
	result.curvatureP = dot(p, ap) / btMax(btScalar(1e-30), pp);
	result.curvatureQ = dot(q, aq) / btMax(btScalar(1e-30), qq);
	A.multiply(p, work);
	for (int n = 0; n < size; ++n) work[n] -= ap[n];
	result.repeatabilityError = btSqrt(dot(work, work)) / btMax(btScalar(1e-30), btSqrt(dot(ap, ap)));
	for (int n = 0; n < size; ++n) combined[n] = p[n] + q[n];
	A.multiply(combined, work);
	for (int n = 0; n < size; ++n) work[n] -= ap[n] + aq[n];
	result.linearityError = btSqrt(dot(work, work)) / btMax(btScalar(1e-30), btSqrt(dot(ap, ap)) + btSqrt(dot(aq, aq)));
	for (int n = 0; n < size; ++n) combined[n].setZero();
	A.multiply(combined, work);
	result.zeroNorm = btSqrt(dot(work, work));
	return result;
}

struct PCGComparison
{
	int iterations;
	bool breakdown;
	btScalar initialResidual, finalResidual, target;
};

// Comparison only: checks positive curvature, keeps the best residual iterate,
// and reports a freshly recomputed physical residual. Never used by production solve().
template<class Matrix>
PCGComparison comparePCG(Matrix& A, Vectors& x, const Vectors& rhs, int maxIterations, btScalar target)
{
	PCGComparison result = {};
	result.target = target;
	Vectors r, z, p, ap, best;
	r.resize(rhs.size()); z.resize(rhs.size()); p.resize(rhs.size()); ap.resize(rhs.size());
	A.multiply(x, ap);
	for (int n = 0; n < rhs.size(); ++n) r[n] = rhs[n] - ap[n];
	btScalar bestNorm = btSqrt(dot(r, r));
	result.initialResidual = bestNorm;
	if (!std::isfinite((double)bestNorm)) { result.breakdown = true; result.finalResidual = bestNorm; return result; }
	if (bestNorm <= target) { result.finalResidual = bestNorm; return result; }
	best = x;
	A.precondition(r, z);
	p = z;
	btScalar rz = dot(r, z);
	for (int k = 0; k < maxIterations && bestNorm > target; ++k)
	{
		A.multiply(p, ap);
		const btScalar curvature = dot(p, ap);
		if (!(curvature > 0) || !(rz > 0) || !std::isfinite((double)curvature) || !std::isfinite((double)rz))
		{
			result.breakdown = true;
			break;
		}
		const btScalar alpha = rz / curvature;
		for (int n = 0; n < rhs.size(); ++n) { x[n] += alpha * p[n]; r[n] -= alpha * ap[n]; }
		result.iterations = k + 1;
		const btScalar residualNorm = btSqrt(dot(r, r));
		if (!std::isfinite((double)residualNorm)) { result.breakdown = true; break; }
		if (residualNorm < bestNorm) { bestNorm = residualNorm; best = x; }
		if (bestNorm <= target) break;
		A.precondition(r, z);
		const btScalar nextRz = dot(r, z);
		if (!(nextRz > 0) || !std::isfinite((double)nextRz)) { result.breakdown = true; break; }
		const btScalar beta = nextRz / rz;
		for (int n = 0; n < rhs.size(); ++n) p[n] = z[n] + beta * p[n];
		rz = nextRz;
	}
	x = best;
	A.multiply(x, ap);
	for (int n = 0; n < rhs.size(); ++n) r[n] = rhs[n] - ap[n];
	result.finalResidual = btSqrt(dot(r, r));
	return result;
}

}  // namespace btDeformableTest
#endif
