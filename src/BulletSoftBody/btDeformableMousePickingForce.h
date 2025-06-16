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
	std::vector<const btSoftBody::Node*> m_nodes;
	// To grab not only by the positional DOFs, but also by the rotational DOFs, thus making the body retain its original orientation when it was first grabbed,
	// we use all the vertices of the grabbed primitive, be it triangle or a tetra. The m_vertex_to_barycenter is remembered when first grabbed and then when the
	// mouse is moving, these offsets are added to the new mouse pos. This creates forces which push the vertices of the primitive, so that the original orientation
	// is achieved. This has a side effect that the grabbed primitive resists deformation. Remains to be seen if that is a problem.
	// TODO it is theorised that if vertices from other more distant tetras were used to retain orientation, the stabilizing effect would be much stronger.
	btVector3 m_node_to_mouse_x[4], m_node_to_mouse_x_orig[4];
	btVector3 m_node_to_mouse_q[4], m_node_to_mouse_q_orig[4];
	btTransform m_mouse_transform, m_mouse_transform_orig, m_mouse_transform_orig_inv;
	btScalar m_maxForce;

	const btSoftBody::Node* getNode(int i) const
	{
		if (m_face)
			return m_face->m_n[i];
		else if (m_tetra)
			return m_tetra->m_n[i];
		else if (m_nodes[i])
			return m_nodes[i];
		return nullptr;
	}

	int getIndexCount() const
	{
		if (m_face)
			return 3;
		else if (m_tetra)
			return 4;
		else
			return m_nodes.size();
	}

	void calculateNodeToMouse()
	{
		for (int i = 0; i < getIndexCount(); ++i)
		{
			m_node_to_mouse_x[i] = getNode(i)->m_x - m_mouse_transform.getOrigin();
			m_node_to_mouse_x_orig[i] = m_node_to_mouse_x[i];
			m_node_to_mouse_q[i] = getNode(i)->m_q - m_mouse_transform.getOrigin();
			m_node_to_mouse_q_orig[i] = m_node_to_mouse_q[i];
		}
	}

	// This makes the soft body respect the controller orientation
	void rotateNodeToMouse()
	{
		// m_mouse_transform_orig * x = m_mouse_transform
		auto x = m_mouse_transform * m_mouse_transform_orig_inv;
		for (int i = 0; i < getIndexCount(); ++i)
		{
			m_node_to_mouse_x[i] = x.getBasis() * m_node_to_mouse_x_orig[i];
			m_node_to_mouse_q[i] = x.getBasis() * m_node_to_mouse_q_orig[i];
		}
	}

public:
	typedef btAlignedObjectArray<btVector3> TVStack;
	btDeformableMousePickingForce(btScalar k, btScalar d, const btSoftBody::Face* face, const btSoftBody::Tetra* tetra, std::vector<const btSoftBody::Node*>* nodes, const btTransform& mouse_transform, btScalar maxForce = 0.3) : m_elasticStiffness(k), m_dampingStiffness(d), m_face(face), m_tetra(tetra), m_mouse_transform(mouse_transform), m_mouse_transform_orig(mouse_transform), m_maxForce(maxForce)
	{
		if (nodes)
			m_nodes = *nodes;
		m_mouse_transform_orig_inv = m_mouse_transform_orig.inverse();
		calculateNodeToMouse();
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
			btVector3 diffForVert = (getNode(i)->m_x - m_mouse_transform.getOrigin()) - m_node_to_mouse_x[i];
			btVector3 v_diff = getNode(i)->m_v;
			btVector3 scaled_force = scale * m_dampingStiffness * v_diff;
			if (diffForVert.norm() > SIMD_EPSILON)
			{
				btVector3 dir = diffForVert.normalized();
				scaled_force = scale * m_dampingStiffness * v_diff.dot(dir) * dir;
			}
			force[getNode(i)->index] -= scaled_force;
		}
	}

	virtual void addScaledElasticForce(btScalar scale, TVStack& force)
	{
		btScalar scaled_stiffness = scale * m_elasticStiffness;
		for (int i = 0; i < getIndexCount(); ++i)
		{
			btVector3 diffForVert = (getNode(i)->m_x - m_mouse_transform.getOrigin()) - m_node_to_mouse_x[i];
			btVector3 scaled_force = scaled_stiffness * diffForVert;
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
			btVector3 diffForVert = (getNode(i)->m_x - m_mouse_transform.getOrigin()) - m_node_to_mouse_x[i];
			btVector3 local_scaled_df = scaled_k_damp * dv[getNode(i)->index];
			auto diffForVertLength = diffForVert.length();
			if (diffForVertLength > SIMD_EPSILON)
			{
				btVector3 dir = diffForVert / diffForVertLength;
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
			btVector3 diffForVert = (getNode(i)->m_q - m_mouse_transform.getOrigin()) - m_node_to_mouse_q[i];
			btVector3 scaled_force = m_elasticStiffness * diffForVert;
			if (scaled_force.safeNorm() > m_maxForce)
			{
				scaled_force.safeNormalize();
				scaled_force *= m_maxForce;
			}
			energy += 0.5 * scaled_force.dot(diffForVert);
		}
		return energy;
	}

	virtual double totalDampingEnergy(btScalar dt)
	{
		double energy = 0;
		for (int i = 0; i < getIndexCount(); ++i)
		{
			btVector3 diffForVert = (getNode(i)->m_x - m_mouse_transform.getOrigin()) - m_node_to_mouse_x[i];
			btVector3 v_diff = getNode(i)->m_v;
			btVector3 scaled_force = m_dampingStiffness * v_diff;
			if (diffForVert.norm() > SIMD_EPSILON)
			{
				btVector3 dir = diffForVert.normalized();
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
			btVector3 diffForVert = (getNode(i)->m_q - m_mouse_transform.getOrigin()) - m_node_to_mouse_q[i];
			btScalar dir_norm = diffForVert.norm();
			btVector3 dir_normalized = (dir_norm > SIMD_EPSILON) ? diffForVert.normalized() : btVector3(0, 0, 0);
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
		m_mouse_transform.setOrigin(p);
	}

	void setMouseTransform(const btTransform& t)
	{
		m_mouse_transform = t;
		rotateNodeToMouse();
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
