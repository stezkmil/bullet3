#ifndef BT_DEFORMABLE_CONTACT_REFRESH_H
#define BT_DEFORMABLE_CONTACT_REFRESH_H

#include "btSoftBody.h"

// Keep distinct support points on an unchanged contact plane. Replace stale planes
// only; equivalent stencils are merged conservatively by the contact force.
inline int btRemoveRefreshedContactPatches(
	btAlignedObjectArray<btSoftBody::DeformableNodeNodeContact>& current,
	const btAlignedObjectArray<btSoftBody::DeformableNodeNodeContact>& refreshed)
{
	int kept = 0;
	const int original = current.size();
	for (int i = 0; i < original; ++i)
	{
		bool replaced = false;
		for (int j = 0; j < refreshed.size(); ++j)
			if (!current[i].m_surfaceInvalid && !refreshed[j].m_surfaceInvalid && current[i].sameSurfaceFeature(refreshed[j]))
			{
				const auto& old = current[i];
				const auto& fresh = refreshed[j];
				const bool sameOrder = old.m_surfaceObjects[0] == fresh.m_surfaceObjects[0] &&
									   old.m_surfaceParts[0] == fresh.m_surfaceParts[0] && old.m_surfaceTriangles[0] == fresh.m_surfaceTriangles[0];
				const btScalar orientation = sameOrder ? btScalar(1) : btScalar(-1);
				if (old.m_normal.length2() <= SIMD_EPSILON || fresh.m_normal.length2() <= SIMD_EPSILON) continue;
				const bool samePlane = old.m_normal.normalized().dot(fresh.m_normal.normalized()) * orientation > btScalar(0.999999);
				if (samePlane)
				{
					replaced = false;
					break;
				}
				replaced = true;
			}
		if (!replaced)
		{
			if (kept != i) current[kept] = current[i];
			++kept;
		}
	}
	current.resize(kept);
	return original - kept;
}

#endif
