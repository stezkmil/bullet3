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
	// ------------------------------------------------------------------ params
	btScalar m_elasticStiffness;  //!< base Hooke constant (N/m)
	btScalar m_dampingStiffness;  //!< base damping constant (N·s/m)

	btScalar m_angleMax;   //!< angle at which the boost saturates (rad)
	btScalar m_scaleCeil;  //!< maximal multiplier (≥1)

	const btSoftBody::Face* m_face;
	const btSoftBody::Tetra* m_tetra;
	std::vector<const btSoftBody::Node*> m_nodes;

	btVector3 m_node_to_mouse_x[4], m_node_to_mouse_x_orig[4];
	btVector3 m_node_to_mouse_q[4], m_node_to_mouse_q_orig[4];

	btVector4 m_mouse_bary;

	btTransform m_mouse_transform, m_mouse_transform_orig, m_mouse_transform_orig_inv;
	btScalar m_maxForce;  //!< per‑spring clamping
	btScalar m_kRot = 1.0;

	// --------------------------------------------------------- helpers
	const btSoftBody::Node* getNode(int i) const
	{
		if (m_face) return m_face->m_n[i];
		if (m_tetra) return m_tetra->m_n[i];
		return m_nodes[i];
	}

	int getIndexCount() const
	{
		if (m_face) return 3;
		if (m_tetra) return 4;
		return static_cast<int>(m_nodes.size());
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
		m_mouse_bary = BaryCoord(getNode(0)->m_x, getNode(1)->m_x, getNode(2)->m_x, getNode(3)->m_x, m_mouse_transform.getOrigin());
	}

	void rotateNodeToMouse()
	{
		btTransform delta = m_mouse_transform * m_mouse_transform_orig_inv;
		for (int i = 0; i < getIndexCount(); ++i)
		{
			m_node_to_mouse_x[i] = delta.getBasis() * m_node_to_mouse_x_orig[i];
			m_node_to_mouse_q[i] = delta.getBasis() * m_node_to_mouse_q_orig[i];
		}
	}

	static inline int PolarDecomposition(const btMatrix3x3& m, btMatrix3x3& q, btMatrix3x3& s)
	{
		static const btPolarDecomposition polar;
		return polar.decompose(m, q, s);
	}

	void orientationError(btScalar& theta, btVector3& axis) const
	{
		int n = getIndexCount();
		btVector3 c_now(0, 0, 0), c_rest(0, 0, 0);
		for (int i = 0; i < n; ++i)
		{
			c_now += getNode(i)->m_x;
			c_rest += m_node_to_mouse_x[i];
		}
		c_now *= btScalar(1.0) / n;
		c_rest *= btScalar(1.0) / n;

		btVector3 p[4], q[4];
		for (int i = 0; i < n; ++i)
		{
			p[i] = getNode(i)->m_x - c_now;
			q[i] = m_node_to_mouse_x[i] - c_rest;
		}

		btMatrix3x3 H(btScalar(0), btScalar(0), btScalar(0),
					  btScalar(0), btScalar(0), btScalar(0),
					  btScalar(0), btScalar(0), btScalar(0));
		for (int i = 0; i < n; ++i)
		{
			H[0][0] += p[i].x() * q[i].x();
			H[0][1] += p[i].x() * q[i].y();
			H[0][2] += p[i].x() * q[i].z();
			H[1][0] += p[i].y() * q[i].x();
			H[1][1] += p[i].y() * q[i].y();
			H[1][2] += p[i].y() * q[i].z();
			H[2][0] += p[i].z() * q[i].x();
			H[2][1] += p[i].z() * q[i].y();
			H[2][2] += p[i].z() * q[i].z();
		}
		btMatrix3x3 R, S;
		static const btPolarDecomposition polar;
		polar.decompose(H, R, S);
		{
			btQuaternion q;
			R.getRotation(q);
			axis = q.getAxis();
			theta = q.getAngle();
		}
	}

	btScalar transBoost(btScalar theta) const
	{
		if (m_angleMax <= btScalar(0)) return btScalar(1);
		btScalar b = btScalar(1) + theta / m_angleMax;
		if (b > m_scaleCeil) b = m_scaleCeil;
		return b;
	}

	// ---------------------------------------------------------------------
	void applyTorqueSpring(const btVector3& T, TVStack& force)
	{
		int n = getIndexCount();
		if (n == 0) return;
		//btVector3 picked = /*m_mouse_transform*/ m_mouse_transform_orig.getOrigin();
		btScalar denom = btScalar(0);
		btVector3 r[4];
		btVector3 pickedFromBary = BaryEval(getNode(0)->m_x, getNode(1)->m_x, getNode(2)->m_x, getNode(3)->m_x, m_mouse_bary);
		for (int i = 0; i < n; ++i)
		{
			r[i] = getNode(i)->m_x - pickedFromBary;
			denom += r[i].length2();
		}
		if (denom < btScalar(1e-7)) return;
		for (int i = 0; i < n; ++i)
		{
			btVector3 dbg(0, 0, 1);
			btVector3 F = /*T.*/ dbg.cross(r[i].normalized()) * 10.0 /*/ denom*/;
			/*if (F.safeNorm() > m_maxForce)
			{
				F.safeNormalize();
				F *= m_maxForce;
			}*/
			force[getNode(i)->index] -= F;
		}
	}

public:
	typedef btAlignedObjectArray<btVector3> TVStack;

	btDeformableMousePickingForce(
		btScalar elasticK,
		btScalar dampingK,
		const btSoftBody::Face* face,
		const btSoftBody::Tetra* tetra,
		std::vector<const btSoftBody::Node*>* nodes,
		const btTransform& mouseTr,
		btScalar maxForce = btScalar(0.3),
		btScalar angleMax = 0.0,   // off by default
		btScalar scaleCeil = 1.0)  // no boost by default
		: m_elasticStiffness(elasticK),
		  m_dampingStiffness(dampingK),
		  m_angleMax(angleMax),
		  m_scaleCeil(scaleCeil),
		  m_face(face),
		  m_tetra(tetra),
		  m_mouse_transform(mouseTr),
		  m_mouse_transform_orig(mouseTr),
		  m_maxForce(maxForce)
	{
		if (nodes) m_nodes = *nodes;
		m_mouse_transform_orig_inv = m_mouse_transform_orig.inverse();
		calculateNodeToMouse();
	}

	// ------------------------------------------------- Force interface
	void addScaledForces(btScalar s, TVStack& f) override
	{
		addScaledDampingForce(s, f);
		addScaledElasticForce(s, f);
	}
	void addScaledExplicitForce(btScalar s, TVStack& f) override { addScaledElasticForce(s, f); }

	// ---------------- elastic ------------------------------------------------
	void addScaledElasticForce(btScalar scale, TVStack& force)
	{
		btScalar theta;
		btVector3 axis;
		orientationError(theta, axis);
		btScalar orientationBasedStiffnessBoost = transBoost(theta);
		btScalar k = scale * m_elasticStiffness * orientationBasedStiffnessBoost;
		for (int i = 0; i < getIndexCount(); ++i)
		{
			const btSoftBody::Node* n = getNode(i);
			btVector3 diff = (n->m_x - m_mouse_transform.getOrigin()) - m_node_to_mouse_x[i];
			btVector3 f = k * diff;
			if (f.safeNorm() > m_maxForce)
			{
				f.safeNormalize();
				f *= m_maxForce;
			}
			force[n->index] -= f;
		}
		// Should add a rotation to preserve the grab orientation, but I did not observe any positive effects even with high values of m_kRot
		if (m_kRot > btScalar(0) /*&& theta > btScalar(1e-4)*/)
		{
			btVector3 torque = m_kRot * theta * axis;
			applyTorqueSpring(torque, force);
		}
	}

	// --------------- damping -------------------------------------------------
	void addScaledDampingForce(btScalar scale, TVStack& force) override
	{
		btScalar theta;
		btVector3 axis;
		orientationError(theta, axis);
		btScalar orientationBasedStiffnessBoost = transBoost(theta);
		btScalar k = scale * m_dampingStiffness * orientationBasedStiffnessBoost;
		for (int i = 0; i < getIndexCount(); ++i)
		{
			const btSoftBody::Node* n = getNode(i);
			btVector3 diff = (n->m_x - m_mouse_transform.getOrigin()) - m_node_to_mouse_x[i];
			btVector3 v = n->m_v;
			btVector3 f = k * v;
			if (diff.norm() > SIMD_EPSILON)
			{
				btVector3 dir = diff.normalized();
				f = k * v.dot(dir) * dir;
			}
			force[n->index] -= f;
		}
		// Should add a rotation to preserve the grab orientation, but I did not observe any positive effects even with high values of m_kRot
		if (m_kRot > btScalar(0) /*&& theta > btScalar(1e-4)*/)
		{
			btVector3 torque = m_kRot * theta * axis;
			applyTorqueSpring(torque, force);
		}
	}

	// Jacobian‑vector products and energies are kept identical to the previous
	// version except that they use the *base* stiffness (orientation boosting
	// does not affect the linearised system for stability).

	void addScaledDampingForceDifferential(btScalar scale, const TVStack& dv, TVStack& df) override
	{
		btScalar theta;
		btVector3 axis;
		orientationError(theta, axis);
		btScalar orientationBasedStiffnessBoost = transBoost(theta);
		btScalar k = m_dampingStiffness * scale * orientationBasedStiffnessBoost;  // keep linear
		for (int i = 0; i < getIndexCount(); ++i)
		{
			const btSoftBody::Node* n = getNode(i);
			btVector3 diff = (n->m_x - m_mouse_transform.getOrigin()) - m_node_to_mouse_x[i];
			btVector3 local = k * dv[n->index];
			if (diff.length() > SIMD_EPSILON)
			{
				btVector3 dir = diff.normalized();
				local = k * dv[n->index].dot(dir) * dir;
			}
			df[n->index] -= local;
		}
	}
	void buildDampingForceDifferentialDiagonal(btScalar, TVStack&) override {}

	double totalElasticEnergy(btScalar) override
	{
		btScalar theta;
		btVector3 axis;
		orientationError(theta, axis);
		btScalar orientationBasedStiffnessBoost = transBoost(theta);
		double e = 0;
		for (int i = 0; i < getIndexCount(); ++i)
		{
			btVector3 diff = (getNode(i)->m_q - m_mouse_transform.getOrigin()) - m_node_to_mouse_q[i];
			btVector3 f = m_elasticStiffness * diff * orientationBasedStiffnessBoost;
			if (f.safeNorm() > m_maxForce)
			{
				f.safeNormalize();
				f *= m_maxForce;
			}
			e += 0.5 * f.dot(diff);
		}
		return e;
	}
	double totalDampingEnergy(btScalar dt) override
	{
		double e = 0;
		btScalar theta;
		btVector3 axis;
		orientationError(theta, axis);
		btScalar orientationBasedStiffnessBoost = transBoost(theta);
		for (int i = 0; i < getIndexCount(); ++i)
		{
			const btSoftBody::Node* n = getNode(i);
			btVector3 diff = (n->m_x - m_mouse_transform.getOrigin()) - m_node_to_mouse_x[i];
			btVector3 v = n->m_v;
			btVector3 f = m_dampingStiffness * v * orientationBasedStiffnessBoost;
			if (diff.norm() > SIMD_EPSILON)
			{
				btVector3 dir = diff.normalized();
				f = m_dampingStiffness * v.dot(dir) * dir * orientationBasedStiffnessBoost;
			}
			e -= f.dot(v) / dt;
		}
		return e;
	}

	void addScaledElasticForceDifferential(btScalar scale, const TVStack& dx, TVStack& df) override
	{
		btScalar theta;
		btVector3 axis;
		orientationError(theta, axis);
		btScalar orientationBasedStiffnessBoost = transBoost(theta);
		btScalar k = scale * m_elasticStiffness * orientationBasedStiffnessBoost;  // linear diff uses base k
		for (int i = 0; i < getIndexCount(); ++i)
		{
			const btSoftBody::Node* n = getNode(i);
			btVector3 diff = (n->m_q - m_mouse_transform.getOrigin()) - m_node_to_mouse_q[i];
			btScalar len = diff.length();
			btVector3 dir = (len > SIMD_EPSILON) ? diff / len : btVector3(0, 0, 0);
			int id = n->index;
			btVector3 dx_i = dx[id];
			btVector3 local(0, 0, 0);
			if (len > SIMD_EPSILON)
			{
				btScalar r = 0;
				local -= k * dir.dot(dx_i) * dir;
				local += k * dir.dot(dx_i) * ((len - r) / len) * dir;
				local -= k * ((len - r) / len) * dx_i;
			}
			df[id] += local;
		}
	}

	// ------------------------------------------------ setters ---------------
	void setMousePos(const btVector3& p) { m_mouse_transform.setOrigin(p); }
	void setMouseTransform(const btTransform& t)
	{
		m_mouse_transform = t;
		rotateNodeToMouse();
	}
	void setMaxForce(btScalar f) { m_maxForce = f; }
	void setElasticStiffness(btScalar k) { m_elasticStiffness = k; }
	void setDampingStiffness(btScalar k) { m_dampingStiffness = k; }

	btDeformableLagrangianForceType getForceType() override { return BT_MOUSE_PICKING_FORCE; }
};

#endif /* btMassSpring_h */
