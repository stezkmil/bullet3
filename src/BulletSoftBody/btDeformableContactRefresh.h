#ifndef BT_DEFORMABLE_CONTACT_REFRESH_H
#define BT_DEFORMABLE_CONTACT_REFRESH_H

#include "btSoftBody.h"

// Replace whole triangle-pair patches, preserving every point in the refreshed patch.
// Contacts without feature identity and patches absent from the refresh remain active.
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
			{ replaced = true; break; }
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
