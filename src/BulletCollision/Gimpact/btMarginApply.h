#ifndef BT_MARGIN_APPLY_H_INCLUDED
#define BT_MARGIN_APPLY_H_INCLUDED

#include "LinearMath/btVector3.h"

template <typename CLASS_POINT>
struct btMarginApply
{
	// This is the original way of how the margin was applied to a triangle. It has its advantages and disadvantages. Advantage is that the boxset creation is fast
	// and works for non-manifold meshes. For example it will work flawlessly for a simple plane. The plane can hit even with its backface and this uniform
	// AABB expansion will ensure that even the backface is covered by the AABB. Disadvantage is that it causes combinatorial explosions of tri tri
	// pairs when the margins are relatively big (even 1mm is typically a problem) during collision detection runtime.
	static void apply_margin_old(btScalar margin, CLASS_POINT &minPt, CLASS_POINT &maxPt)
	{
		minPt[0] -= margin;
		minPt[1] -= margin;
		minPt[2] -= margin;
		maxPt[0] += margin;
		maxPt[1] += margin;
		maxPt[2] += margin;
	}

	// The new method is not perfect for shell meshes. It can happen that the backfaces of two planes meet, and the triangles are not covered by the AABBs
	// because the AABBs are biased in the normal direcion. It is also slower in the boxset preparation time, but much faster in the collision detection runtime
	// for high tri count meshes. Further experimentation - if there are some unforeseen problems, we might try to apply this new method only for higher
	// margin values. For example fall back to the old method for margin <= 0.1 units.
	static void apply_margin(const CLASS_POINT &V1,
							 const CLASS_POINT &V2,
							 const CLASS_POINT &V3, btScalar margin, CLASS_POINT &minPt, CLASS_POINT &maxPt)
	{
		if (margin == btScalar(0))
			return;

		if (margin <= 0.2)  // So that the default value of 0.1 safely falls into the old method. It is only around 1.0 where things already start to go really slow,
							// so for those values we use the new method.
		{
			apply_margin_old(margin, minPt, maxPt);
			return;
		}

		// Compute geometric normal and shift the whole box by margin in the normal direction
		const btVector3 p1(V1[0], V1[1], V1[2]);
		const btVector3 p2(V2[0], V2[1], V2[2]);
		const btVector3 p3(V3[0], V3[1], V3[2]);

		btVector3 n = (p2 - p1).cross(p3 - p1);  // unnormalized normal
		const btScalar len2 = n.length2();

		if (len2 > btScalar(1e-24))
		{
			n *= btScalar(1) / btSqrt(len2);

			const btVector3 d = margin * n;  // translation vector

			minPt[0] += d.x();
			maxPt[0] += d.x();
			minPt[1] += d.y();
			maxPt[1] += d.y();
			minPt[2] += d.z();
			maxPt[2] += d.z();

			// This below makes sure that the AABB has a minimum size, to avoid an AABB with a zero dimension. It is
			// tempting to have the minimum to be the margin value, but that would completely remove the benefits of this new
			// method (not causing combinatorial explosions for large margins). So some minimum size independent of the margin is used
			// (based on maximum triangle edge length).

			const btVector3 a(V1[0], V1[1], V1[2]);
			const btVector3 b(V2[0], V2[1], V2[2]);
			const btVector3 c(V3[0], V3[1], V3[2]);

			const btScalar e1 = (b - a).length();
			const btScalar e2 = (c - b).length();
			const btScalar e3 = (a - c).length();
			const btScalar L = btMax(e1, btMax(e2, e3));  // Longest triangle edge

			// Smallest AABB dimension hard limit
			const btScalar kAbsEps = btScalar(1e-6);
			// Relative epsilon (fraction of the triangle size)
			const btScalar kRelEps = btScalar(1e-3);

			btScalar minExtentTarget = btMax(kAbsEps, kRelEps * L);
			minExtentTarget = btMin(minExtentTarget, margin);  // Not likely, but just to be sure we don't grow larger than margin

			for (int axis = 0; axis < 3; ++axis)
			{
				const btScalar extent = maxPt[axis] - minPt[axis];
				if (extent < minExtentTarget)
				{
					const btScalar grow = btScalar(0.5) * (minExtentTarget - extent);
					minPt[axis] -= grow;
					maxPt[axis] += grow;
				}
			}
		}
		else
		{
			// Degenerate triangle: be conservative
			apply_margin_old(margin, minPt, maxPt);
		}
	}
};
;

#endif
