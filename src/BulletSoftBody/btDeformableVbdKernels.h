// SPDX-FileCopyrightText: Copyright (c) 2025 The Newton Developers
// SPDX-License-Identifier: Apache-2.0
// C++ adaptation of Newton 009158e62b86 particle/rigid_vbd_kernels.py.
// The determinant guard is adapted from our Newton reference's volume_guard.py.
#ifndef BT_DEFORMABLE_VBD_KERNELS_H
#define BT_DEFORMABLE_VBD_KERNELS_H
#include "LinearMath/btMatrix3x3.h"
#include <algorithm>
#include <cmath>

namespace btVbd
{
inline btMatrix3x3 columns(const btVector3 &a, const btVector3 &b, const btVector3 &c)
{
	return btMatrix3x3(a.x(), b.x(), c.x(), a.y(), b.y(), c.y(), a.z(), b.z(), c.z());
}
inline btMatrix3x3 outer(const btVector3 &a, const btVector3 &b)
{
	return columns(a * b.x(), a * b.y(), a * b.z());
}
inline btMatrix3x3 deformation(const btVector3 *x, const btMatrix3x3 &inverseRest)
{
	return columns(x[1] - x[0], x[2] - x[0], x[3] - x[0]) * inverseRest;
}
inline void tetraBlock(const btVector3 *x, const btVector3 *previous, const btMatrix3x3 &inverseRest, btScalar volume, int vertex,
					   btScalar mu, btScalar lambda, btScalar damping, btScalar dt, btVector3 &force, btMatrix3x3 &hessian)
{
	const btMatrix3x3 F = deformation(x, inverseRest);
	const btVector3 f[] = {F.getColumn(0), F.getColumn(1), F.getColumn(2)};
	const btMatrix3x3 cof = columns(f[1].cross(f[2]), f[2].cross(f[0]), f[0].cross(f[1]));
	const btVector3 m = vertex ? inverseRest[vertex - 1] : -(inverseRest[0] + inverseRest[1] + inverseRest[2]);
	const btScalar l = lambda + mu;
	const btScalar alpha = 1 + mu / btMax(l, btScalar(1e-6));
	const btVector3 g = cof * m;
	force -= volume * (mu * (F * m) + l * (F.determinant() - alpha) * g);
	hessian += (btMatrix3x3::getIdentity() * (mu * m.length2()) + outer(g, g) * l) * volume;
	if (damping == 0)
		return;
	const btMatrix3x3 F0 = deformation(previous, inverseRest);
	for (int i = 0; i < 3; ++i)
		for (int j = i; j < 3; ++j)
		{
			const btVector3 gradient = f[j] * m[i] + f[i] * m[j];
			const btScalar weight = i == j ? 1 : 2;
			const btScalar rate = (f[i].dot(f[j]) - F0.getColumn(i).dot(F0.getColumn(j))) / dt;
			force -= gradient * (volume * damping * weight * rate);
			hessian += outer(gradient, gradient) * (volume * damping * weight / dt);
		}
}
inline void contactBlock(btScalar distance, btScalar radius, const btVector3 &normal, const btVector3 &translation, btScalar ke,
						 btScalar kd, btScalar friction, btScalar epsilon, btScalar dt, btVector3 &force, btMatrix3x3 &hessian)
{
	if (distance >= radius)
		return;
	const btScalar load = ke * (radius - distance);
	const btMatrix3x3 nn = outer(normal, normal);
	force += normal * load;
	hessian += nn * ke;
	const btScalar normalMotion = normal.dot(translation);
	if (normalMotion < 0)
	{
		force -= normal * (kd / dt * normalMotion);
		hessian += nn * (kd / dt);
	}
	const btVector3 tangent = translation - normal * normalMotion;
	const btScalar length = tangent.length();
	if (length > 0)
	{
		const btScalar threshold = epsilon * dt;
		const btScalar scale = friction * load * (length > threshold ? 1 / length : (2 - length / threshold) / threshold);
		force -= tangent * scale;
		hessian += (btMatrix3x3::getIdentity() - nn) * scale;
	}
}
inline btScalar planarBound(const btVector3 &x, const btVector3 &displacement, const btVector3 &normal, const btVector3 &plane,
							btScalar epsilon)
{
	const btScalar s0 = normal.dot(x - plane) - epsilon, s1 = s0 + normal.dot(displacement);
	if (s0 < 0)
		return s1 >= s0 ? 1 : 0;
	if (s1 >= 0)
		return 1;
	return btClamped(btScalar(.85) * s0 / (s0 - s1), btScalar(0), btScalar(1));
}
inline btScalar polynomial(btScalar a, btScalar b, btScalar c, btScalar d, btScalar t)
{
	return ((a * t + b) * t + c) * t + d;
}
inline btScalar firstVolumeBound(btScalar a, btScalar b, btScalar c, btScalar d)
{
	if (!std::isfinite(double(a + b + c + d)))
		return 0;
	if (d <= 0)
	{
		d = 0;
		if (c < 0 || (c == 0 && b < 0) || (c == 0 && b == 0 && a < 0))
			return 0;
		if (a == 0 && b == 0 && c == 0)
			return 1;
	}
	btScalar cuts[] = {0, 1, 1, 1};
	if (btFabs(a) > 1e-30)
	{
		const btScalar disc = b * b - 3 * a * c;
		if (disc >= 0)
		{
			btScalar r0 = (-b - btSqrt(disc)) / (3 * a), r1 = (-b + btSqrt(disc)) / (3 * a);
			if (r0 > r1)
				std::swap(r0, r1);
			cuts[1] = btClamped(r0, btScalar(0), btScalar(1));
			cuts[2] = btClamped(r1, btScalar(0), btScalar(1));
		}
	}
	else if (btFabs(b) > 1e-30)
		cuts[1] = btClamped(-c / (2 * b), btScalar(0), btScalar(1));
	for (int k = 1; k < 4; ++k)
		if (cuts[k] > 0 && polynomial(a, b, c, d, cuts[k]) <= 0)
		{
			btScalar lo = cuts[k - 1], hi = cuts[k];
			for (int i = 0; i < 40; ++i)
			{
				const btScalar mid = (lo + hi) * btScalar(.5);
				if (polynomial(a, b, c, d, mid) > 0)
					lo = mid;
				else
					hi = mid;
			}
			return btScalar(.9) * lo;
		}
	return 1;
}
inline btScalar volumeBound(const btVector3 *old, const btVector3 *proposed, btScalar restDet, btScalar floor)
{
	const btVector3 e1 = old[1] - old[0], e2 = old[2] - old[0], e3 = old[3] - old[0];
	const btVector3 u1 = proposed[1] - proposed[0] - e1, u2 = proposed[2] - proposed[0] - e2, u3 = proposed[3] - proposed[0] - e3;
	const btScalar det = e1.dot(e2.cross(e3));
	return firstVolumeBound(u1.dot(u2.cross(u3)), u1.dot(u2.cross(e3)) + u1.dot(e2.cross(u3)) + e1.dot(u2.cross(u3)),
							u1.dot(e2.cross(e3)) + e1.dot(u2.cross(e3)) + e1.dot(e2.cross(u3)), det - btMin(floor * restDet, det));
}
struct ScalarMappingSupport
{
	int node;
	btScalar weight;
};
#ifdef _MSC_VER
#pragma float_control(precise, on, push)
#endif
// Keep scalar accumulation ordered like the matrix evaluation it replaces.
inline btVector3 scalarMappingPosition(const btVector3 &offset, const btVector3 *positions, const ScalarMappingSupport *begin,
									   const ScalarMappingSupport *end)
{
	btVector3 result = offset;
	for (auto support = begin; support != end; ++support)
	{
		const auto &p = positions[support->node];
		const btScalar weight = support->weight;
		result += btVector3(p.x() * weight, p.y() * weight, p.z() * weight);
	}
	return result;
}
#ifdef _MSC_VER
#pragma float_control(pop)
#endif

} // namespace btVbd
#endif
