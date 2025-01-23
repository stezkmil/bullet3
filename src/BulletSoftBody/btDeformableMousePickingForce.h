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

#ifndef BT_MOUSE_PICKING_FORCE_H
#define BT_MOUSE_PICKING_FORCE_H

#include "btDeformableLagrangianForce.h"
#include "btSoftBodyInternals.h"

class btDeformableMousePickingForce : public btDeformableLagrangianForce
{
	// If true, the damping force will be in the direction of the spring
	// If false, the damping force will be in the direction of the velocity
	btScalar m_elasticStiffness, m_dampingStiffness;
	const btSoftBody::Face* m_face;
	const btSoftBody::Tetra* m_tetra;
	btVector4 m_barycenter;
	btScalar m_maxForce;

	// We store the initial local positions of the nodes relative to the
	// element's centroid/orientation. We'll transform these to get
	// the current "desired" node positions in world space.
	btAlignedObjectArray<btVector3> m_initialLocalPositions_x, m_initialLocalPositions_q;

	// The current (desired) transform the user wants to impose (rotation + translation).
	// This is updated every frame (e.g., as the mouse is moved/rotated).
	btTransform m_desiredTransform;

	const btSoftBody::Node* getNode(int i) const
	{
		if (m_face)
			return m_face->m_n[i];
		else if (m_tetra)
			return m_tetra->m_n[i];
		return nullptr;
	}

	int getIndexCount() const
	{
		if (m_face)
			return 3;
		else if (m_tetra)
			return 4;
		return 0;
	}

	btVector3 getPickedPos(bool use_q) const
	{
		btVector3 pos;
		if (m_face)
		{
			btVector3 baryInFace(m_barycenter);
			if (use_q)
				return BaryEval(getNode(0)->m_q, getNode(1)->m_q, getNode(2)->m_q, baryInFace);
			else
				return BaryEval(getNode(0)->m_x, getNode(1)->m_x, getNode(2)->m_x, baryInFace);
		}
		else if (m_tetra)
		{
			if (use_q)
				return BaryEval(getNode(0)->m_q, getNode(1)->m_q, getNode(2)->m_q, getNode(3)->m_q, m_barycenter);
			else
				return BaryEval(getNode(0)->m_x, getNode(1)->m_x, getNode(2)->m_x, getNode(3)->m_x, m_barycenter);
		}
		return pos;
	}

public:
	typedef btAlignedObjectArray<btVector3> TVStack;
	btDeformableMousePickingForce(btScalar k, btScalar d, const btSoftBody::Face* face, const btSoftBody::Tetra* tetra, const btVector4& barycenter, const btTransform& initialTransform, btScalar maxForce = 0.3) : m_elasticStiffness(k), m_dampingStiffness(d), m_face(face), m_tetra(tetra), m_barycenter(barycenter), m_maxForce(maxForce), m_desiredTransform(initialTransform)
	{
		// Precompute the local positions of each node relative to the centroid
		// (and possibly orientation, if you want to do local frames).
		btVector3 picked_pos_x = getPickedPos(false);
		btVector3 picked_pos_q = getPickedPos(true);
		int count = getIndexCount();
		m_initialLocalPositions_x.resize(count);
		m_initialLocalPositions_q.resize(count);

		for (int i = 0; i < count; ++i)
		{
			btVector3 worldPos = getNode(i)->m_x;
			// local position relative to centroid
			m_initialLocalPositions_x[i] = worldPos - picked_pos_x;

			worldPos = getNode(i)->m_q;
			// local position relative to centroid
			m_initialLocalPositions_q[i] = worldPos - picked_pos_q;
		}
	}

	virtual void addScaledForces(btScalar scale, TVStack& force)
	{
		addScaledDampingForce(scale, force);
		addScaledElasticForce(scale, force);
	}

	virtual void addScaledExplicitForce(btScalar scale, TVStack& force)
	{
		addScaledElasticForce(scale, force);
	}

	virtual void addScaledDampingForce(btScalar scale, TVStack& force)
	{
		btScalar scaledDamp = m_dampingStiffness * scale;

		int count = getIndexCount();
		if (count == 0) return;

		// For each node, compute velocity damping in the direction
		// of the 'spring' to the target.
		for (int i = 0; i < count; ++i)
		{
			const btSoftBody::Node* node = getNode(i);
			btVector3 v = node->m_v;  // velocity
			// Desired (target) position for this node
			btVector3 desiredPos = m_desiredTransform * m_initialLocalPositions_x[i];
			// Current position
			btVector3 currPos = node->m_x;
			// Spring direction
			btVector3 diff = currPos - desiredPos;
			btScalar length = diff.length();
			btVector3 dir = length > SIMD_EPSILON ? diff.normalized() : btVector3(0, 0, 0);

			// Project velocity along the spring direction if you like,
			// or simply damp full velocity
			// This matches your original pattern: v_diff.dot(dir)*dir
			btVector3 dampingForce = scaledDamp * (v.dot(dir)) * dir;

			// Subtract from net force
			force[node->index] -= dampingForce;
		}
	}

	virtual void addScaledElasticForce(btScalar scale, TVStack& force)
	{
		btScalar scaled_stiffness = m_elasticStiffness * scale;

		int count = getIndexCount();
		if (count == 0) return;

		for (int i = 0; i < count; ++i)
		{
			const btSoftBody::Node* node = getNode(i);
			// Desired (target) position
			btVector3 desiredPos = m_desiredTransform * m_initialLocalPositions_x[i];
			// Current position
			btVector3 currPos = node->m_x;
			// Displacement
			btVector3 diff = currPos - desiredPos;

			// Force = -k * diff
			btVector3 springForce = -scaled_stiffness * diff;

			// Clamp to max force if desired
			if (m_maxForce > btScalar(0) && springForce.safeNorm() > m_maxForce)
			{
				springForce.safeNormalize();
				springForce *= m_maxForce;
			}

			// Add to the net force array
			force[node->index] += springForce;
		}
	}

	virtual void addScaledDampingForceDifferential(btScalar scale, const TVStack& dv, TVStack& df)
	{
		// Typically for a picking transform, you might do something similar
		// to the above. This is left as an exercise if you need full implicit
		// damping coupling. For simplicity, we can either omit or do a partial version:
		btScalar scaledDamp = m_dampingStiffness * scale;

		int count = getIndexCount();
		if (count == 0) return;

		for (int i = 0; i < count; ++i)
		{
			const btSoftBody::Node* node = getNode(i);
			// Get direction of spring at current state:
			btVector3 desiredPos = m_desiredTransform * m_initialLocalPositions_x[i];
			btVector3 currPos = node->m_x;
			btVector3 diff = currPos - desiredPos;
			btScalar length = diff.length();
			btVector3 dir = (length > SIMD_EPSILON) ? diff.normalized() : btVector3(0, 0, 0);

			// Project velocity differential along that direction
			btVector3 local_df = scaledDamp * (dv[node->index].dot(dir)) * dir;
			df[node->index] -= local_df;
		}
	}

	virtual void buildDampingForceDifferentialDiagonal(btScalar scale, TVStack& diagA) {}

	virtual double totalElasticEnergy(btScalar dt)
	{
		double energy = 0.0;
		int count = getIndexCount();
		if (count == 0) return energy;

		for (int i = 0; i < count; ++i)
		{
			const btSoftBody::Node* node = getNode(i);
			btVector3 desiredPos = m_desiredTransform * m_initialLocalPositions_q[i];
			btVector3 diff = node->m_q - desiredPos;
			// Force magnitude (if below max, otherwise clamp)
			btVector3 scaled_force = m_elasticStiffness * (-diff);
			if (m_maxForce > btScalar(0) && scaled_force.safeNorm() > m_maxForce)
			{
				scaled_force.safeNormalize();
				scaled_force *= m_maxForce;
			}
			// Energy = 0.5 * F • displacement
			energy += 0.5 * scaled_force.dot(diff);
		}
		return energy;
	}

	virtual double totalDampingEnergy(btScalar dt)
	{
		double energy = 0.0;
		int count = getIndexCount();
		if (count == 0) return energy;

		for (int i = 0; i < count; ++i)
		{
			const btSoftBody::Node* node = getNode(i);
			// Same pattern as addScaledDampingForce
			btVector3 v = node->m_v;
			btVector3 desiredPos = m_desiredTransform * m_initialLocalPositions_x[i];
			btVector3 diff = node->m_x - desiredPos;
			btVector3 dir = diff.length() > SIMD_EPSILON ? diff.normalized() : btVector3(0, 0, 0);
			btVector3 dampingForce = m_dampingStiffness * (v.dot(dir)) * dir;

			// Damping energy = - F • v * dt
			energy -= dampingForce.dot(v) / dt;
		}
		return energy;
	}

	virtual void addScaledElasticForceDifferential(btScalar scale, const TVStack& dx, TVStack& df)
	{
		btScalar scaled_stiffness = scale * m_elasticStiffness;
		for (int i = 0; i < getIndexCount(); ++i)
		{
			const btSoftBody::Node* node = getNode(i);
			btVector3 desiredPos = m_desiredTransform * m_initialLocalPositions_q[i];
			btVector3 diff = node->m_q - desiredPos;
			btScalar dir_norm = diff.norm();
			btVector3 dir_normalized = (dir_norm > SIMD_EPSILON) ? diff.normalized() : btVector3(0, 0, 0);

			int id = node->index;
			btVector3 dx_diff = dx[id];
			btScalar r = 0;  // rest length is 0 for picking spring
			btVector3 scaled_df = btVector3(0, 0, 0);
			if (dir_norm > SIMD_EPSILON)
			{
				scaled_df -= scaled_stiffness * dir_normalized.dot(dx_diff) * dir_normalized;
				scaled_df += scaled_stiffness * dir_normalized.dot(dx_diff) * ((dir_norm - r) / dir_norm) * dir_normalized;
				scaled_df -= scaled_stiffness * ((dir_norm - r) / dir_norm) * dx_diff;
			}
			df[id] += scaled_df;
		}
	}

	void setDesiredTransform(const btTransform& t)
	{
		m_desiredTransform = t;
	}

	void setMaxForce(btScalar maxForce)
	{
		m_maxForce = maxForce;
	}

	void setElasticStiffness(btScalar elasticStiffness)
	{
		m_elasticStiffness = elasticStiffness;
	}

	void setDampingStiffness(btScalar dampingStiffness)
	{
		m_dampingStiffness = dampingStiffness;
	}

	virtual btDeformableLagrangianForceType getForceType()
	{
		return BT_MOUSE_PICKING_FORCE;
	}
};

#endif
