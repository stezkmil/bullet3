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

class btDeformableMousePickingForce : public btDeformableLagrangianForce
{
	// If true, the damping force will be in the direction of the spring
	// If false, the damping force will be in the direction of the velocity
	btScalar m_elasticStiffness, m_dampingStiffness;
	const btSoftBody::Face* m_face;
	const btSoftBody::Tetra* m_tetra;
	btVector3 m_mouse_pos;
	btScalar m_maxForce;

	const btSoftBody::Node* getNode(int i) const
	{
		if (m_face)
			return m_face->m_n[i];
		if (m_tetra)
			return m_tetra->m_n[i];
		return nullptr;
	}

	int getIndexCount() const
	{
		if (m_face)
			return 3;
		if (m_tetra)
			return 4;
		return 0;
	}

public:
	typedef btAlignedObjectArray<btVector3> TVStack;
	btDeformableMousePickingForce(btScalar k, btScalar d, const btSoftBody::Face* face, const btSoftBody::Tetra* tetra, const btVector3& mouse_pos, btScalar maxForce = 0.3) : m_elasticStiffness(k), m_dampingStiffness(d), m_face(face), m_tetra(tetra), m_mouse_pos(mouse_pos), m_maxForce(maxForce)
	{
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
		for (int i = 0; i < getIndexCount(); ++i)
		{
			btVector3 v_diff = getNode(i)->m_v;
			btVector3 scaled_force = scale * m_dampingStiffness * v_diff;
			if ((getNode(i)->m_x - m_mouse_pos).norm() > SIMD_EPSILON)
			{
				btVector3 dir = (getNode(i)->m_x - m_mouse_pos).normalized();
				scaled_force = scale * m_dampingStiffness * v_diff.dot(dir) * dir;
			}
			force[getNode(i)->index] -= scaled_force;
		}
	}

	virtual void addScaledElasticForce(btScalar scale, TVStack& force)
	{
		btScalar scaled_stiffness = scale * m_elasticStiffness;
		btVector3 dir;
		for (int i = 0; i < getIndexCount(); ++i)
		{
			if (i == 0)
				dir = (getNode(i)->m_q - m_mouse_pos);
			btVector3 scaled_force = scaled_stiffness * dir;
			if (scaled_force.safeNorm() > m_maxForce)
			{
				scaled_force.safeNormalize();
				scaled_force *= m_maxForce;
			}
			force[getNode(i)->index] -= scaled_force;
		}
	}

	virtual void addScaledDampingForceDifferential(btScalar scale, const TVStack& dv, TVStack& df)
	{
		btScalar scaled_k_damp = m_dampingStiffness * scale;
		for (int i = 0; i < getIndexCount(); ++i)
		{
			btVector3 local_scaled_df = scaled_k_damp * dv[getNode(i)->index];
			if ((getNode(i)->m_x - m_mouse_pos).norm() > SIMD_EPSILON)
			{
				btVector3 dir = (getNode(i)->m_x - m_mouse_pos).normalized();
				local_scaled_df = scaled_k_damp * dv[getNode(i)->index].dot(dir) * dir;
			}
			df[getNode(i)->index] -= local_scaled_df;
		}
	}

	virtual void buildDampingForceDifferentialDiagonal(btScalar scale, TVStack& diagA) {}

	virtual double totalElasticEnergy(btScalar dt)
	{
		double energy = 0;
		for (int i = 0; i < getIndexCount(); ++i)
		{
			btVector3 dir = (getNode(i)->m_q - m_mouse_pos);
			btVector3 scaled_force = m_elasticStiffness * dir;
			if (scaled_force.safeNorm() > m_maxForce)
			{
				scaled_force.safeNormalize();
				scaled_force *= m_maxForce;
			}
			energy += 0.5 * scaled_force.dot(dir);
		}
		return energy;
	}

	virtual double totalDampingEnergy(btScalar dt)
	{
		double energy = 0;
		for (int i = 0; i < getIndexCount(); ++i)
		{
			btVector3 v_diff = getNode(i)->m_v;
			btVector3 scaled_force = m_dampingStiffness * v_diff;
			if ((getNode(i)->m_x - m_mouse_pos).norm() > SIMD_EPSILON)
			{
				btVector3 dir = (getNode(i)->m_x - m_mouse_pos).normalized();
				scaled_force = m_dampingStiffness * v_diff.dot(dir) * dir;
			}
			energy -= scaled_force.dot(getNode(i)->m_v) / dt;
		}
		return energy;
	}

	virtual void addScaledElasticForceDifferential(btScalar scale, const TVStack& dx, TVStack& df)
	{
		btScalar scaled_stiffness = scale * m_elasticStiffness;
		for (int i = 0; i < getIndexCount(); ++i)
		{
			btVector3 dir = (getNode(i)->m_q - m_mouse_pos);
			btScalar dir_norm = dir.norm();
			btVector3 dir_normalized = (dir_norm > SIMD_EPSILON) ? dir.normalized() : btVector3(0, 0, 0);
			int id = getNode(i)->index;
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

	void setMousePos(const btVector3& p)
	{
		m_mouse_pos = p;
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

#endif /* btMassSpring_h */
