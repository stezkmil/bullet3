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
	const btSoftBody::Node* getNode(int i) const;

	int getIndexCount() const;

	void calculateNodeToMouse();

	void rotateNodeToMouse();

	static inline int PolarDecomposition(const btMatrix3x3& m, btMatrix3x3& q, btMatrix3x3& s)
	{
		static const btPolarDecomposition polar;
		return polar.decompose(m, q, s);
	}

	void orientationError(btScalar& theta, btVector3& axis) const;

	btScalar transBoost(btScalar theta) const;

	// ---------------------------------------------------------------------
	void applyTorqueSpring(const btVector3& T, TVStack& force);

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
		btScalar angleMax = 0.0,  // off by default
		btScalar scaleCeil = 1.0);

	// ------------------------------------------------- Force interface
	void addScaledForces(btScalar s, TVStack& f) override;
	void addScaledExplicitForce(btScalar s, TVStack& f) override;

	// ---------------- elastic ------------------------------------------------
	void addScaledElasticForce(btScalar scale, TVStack& force);

	// --------------- damping -------------------------------------------------
	void addScaledDampingForce(btScalar scale, TVStack& force) override;

	// Jacobian‑vector products and energies are kept identical to the previous
	// version except that they use the *base* stiffness (orientation boosting
	// does not affect the linearised system for stability).

	void addScaledDampingForceDifferential(btScalar scale, const TVStack& dv, TVStack& df) override;
	void buildDampingForceDifferentialDiagonal(btScalar, TVStack&) override;

	double totalElasticEnergy(btScalar) override;
	double totalDampingEnergy(btScalar dt) override;

	void addScaledElasticForceDifferential(btScalar scale, const TVStack& dx, TVStack& df) override;

	// ------------------------------------------------ setters ---------------
	void setMousePos(const btVector3& p);
	void setMouseTransform(const btTransform& t);
	void setMaxForce(btScalar f);
	void setElasticStiffness(btScalar k);
	void setDampingStiffness(btScalar k);

	btDeformableLagrangianForceType getForceType() override;
};

#endif /* btMassSpring_h */
