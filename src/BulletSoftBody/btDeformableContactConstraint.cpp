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

#include "btDeformableContactConstraint.h"
/* ================   Deformable Node Anchor   =================== */
btDeformableNodeAnchorConstraint::btDeformableNodeAnchorConstraint(const btSoftBody::DeformableNodeRigidAnchor& a, const btContactSolverInfo& infoGlobal)
	: m_anchor(&a), btDeformableContactConstraint(a.m_cti.m_normal, infoGlobal)
{
}

btDeformableNodeAnchorConstraint::btDeformableNodeAnchorConstraint(const btDeformableNodeAnchorConstraint& other)
	: m_anchor(other.m_anchor), btDeformableContactConstraint(other)
{
}

btVector3 btDeformableNodeAnchorConstraint::getVa() const
{
	const btSoftBody::sCti& cti = m_anchor->m_cti;
	btVector3 va(0, 0, 0);
	if (cti.m_colObj->hasContactResponse())
	{
		btRigidBody* rigidCol = 0;
		btMultiBodyLinkCollider* multibodyLinkCol = 0;

		// grab the velocity of the rigid body
		if (cti.m_colObj->getInternalType() == btCollisionObject::CO_RIGID_BODY)
		{
			rigidCol = (btRigidBody*)btRigidBody::upcast(cti.m_colObj);
			va = rigidCol ? (rigidCol->getVelocityInLocalPoint(m_anchor->m_c1)) : btVector3(0, 0, 0);
		}
		else if (cti.m_colObj->getInternalType() == btCollisionObject::CO_FEATHERSTONE_LINK)
		{
			multibodyLinkCol = (btMultiBodyLinkCollider*)btMultiBodyLinkCollider::upcast(cti.m_colObj);
			if (multibodyLinkCol)
			{
				const int ndof = multibodyLinkCol->m_multiBody->getNumDofs() + 6;
				const btScalar* J_n = &m_anchor->jacobianData_normal.m_jacobians[0];
				const btScalar* J_t1 = &m_anchor->jacobianData_t1.m_jacobians[0];
				const btScalar* J_t2 = &m_anchor->jacobianData_t2.m_jacobians[0];
				const btScalar* local_v = multibodyLinkCol->m_multiBody->getVelocityVector();
				const btScalar* local_dv = multibodyLinkCol->m_multiBody->getDeltaVelocityVector();
				// add in the normal component of the va
				btScalar vel = 0.0;
				for (int k = 0; k < ndof; ++k)
				{
					vel += (local_v[k] + local_dv[k]) * J_n[k];
				}
				va = cti.m_normal * vel;
				// add in the tangential components of the va
				vel = 0.0;
				for (int k = 0; k < ndof; ++k)
				{
					vel += (local_v[k] + local_dv[k]) * J_t1[k];
				}
				va += m_anchor->t1 * vel;
				vel = 0.0;
				for (int k = 0; k < ndof; ++k)
				{
					vel += (local_v[k] + local_dv[k]) * J_t2[k];
				}
				va += m_anchor->t2 * vel;
			}
		}
	}
	return va;
}

// Velocity matching between the deformable node and the rigid/multibody it is anchored to. This is the main constraint solve function for the deformable node anchor constraint. Positional drift correction is handled in btDeformableNodeAnchorConstraint::solveSplitImpulse.
btScalar btDeformableNodeAnchorConstraint::solveConstraint(const btContactSolverInfo& infoGlobal)
{
	const btSoftBody::sCti& cti = m_anchor->m_cti;
	btScalar residualSquare = 0.0;

	if (!m_anchor->m_body->isStaticOrKinematicObject() && m_anchor->m_body->isActive())
	{
		btVector3 va = getVa();
		btVector3 vb = getVb();
		btVector3 vr = (vb - va);
		btScalar residualSquare = btDot(vr, vr);

		const btScalar misconvergenceRelaxationFactor = 0.5;
		const btScalar convergenceRelaxationFactor = 1.1;

		// If the current residual is larger than the last one, we are heading towards an explosion. We use this information as a hint that impulse magnitudes should be dampened.
		// This greatly improves convergence.
		if (m_previous_residual != -1.0)
			if (residualSquare > m_previous_residual)
				m_convergence_based_relaxation = std::max(0.1, m_convergence_based_relaxation * misconvergenceRelaxationFactor);
			else
				m_convergence_based_relaxation = std::min(1.0, m_convergence_based_relaxation * convergenceRelaxationFactor);
		m_previous_residual = residualSquare;

		// dn is the normal component of velocity diffrerence. Approximates the residual. // todo xuchenhan@: this prob needs to be scaled by dt
		btVector3 impulse = (m_anchor->m_c0 * vr) * m_convergence_based_relaxation;

		// TODO impulses applied in iteration 0 can be way off. Warmstart the first impulse. It should work similar to the way the m_appliedImpulse works. This should reduce the iteration count somewhat.
		// TODO still some instabilities were observed when the bodies have 3 orders of magnitude mass difference. To be investigated later. Maybe a mass ratio based clamping of the impulse is needed.

		// apply impulse to deformable nodes involved and change their velocities
		applyImpulse(impulse);

		// apply impulse to the rigid/multibodies involved and change their velocities
		if (cti.m_colObj->getInternalType() == btCollisionObject::CO_RIGID_BODY)
		{
			btRigidBody* rigidCol = 0;
			rigidCol = (btRigidBody*)btRigidBody::upcast(cti.m_colObj);
			if (rigidCol)
			{
				rigidCol->applyImpulse(impulse, m_anchor->m_c1);
			}
		}
		else if (cti.m_colObj->getInternalType() == btCollisionObject::CO_FEATHERSTONE_LINK)
		{
			btMultiBodyLinkCollider* multibodyLinkCol = 0;
			multibodyLinkCol = (btMultiBodyLinkCollider*)btMultiBodyLinkCollider::upcast(cti.m_colObj);
			if (multibodyLinkCol)
			{
				const btScalar* deltaV_normal = &m_anchor->jacobianData_normal.m_deltaVelocitiesUnitImpulse[0];
				// apply normal component of the impulse
				multibodyLinkCol->m_multiBody->applyDeltaVeeMultiDof2(deltaV_normal, impulse.dot(cti.m_normal));
				// apply tangential component of the impulse
				const btScalar* deltaV_t1 = &m_anchor->jacobianData_t1.m_deltaVelocitiesUnitImpulse[0];
				multibodyLinkCol->m_multiBody->applyDeltaVeeMultiDof2(deltaV_t1, impulse.dot(m_anchor->t1));
				const btScalar* deltaV_t2 = &m_anchor->jacobianData_t2.m_deltaVelocitiesUnitImpulse[0];
				multibodyLinkCol->m_multiBody->applyDeltaVeeMultiDof2(deltaV_t2, impulse.dot(m_anchor->t2));
			}
		}

		//fprintf(stderr, "btDeformableNodeAnchorConstraint::solveConstraint m_convergence_based_relaxation %f va %f %f %f vb %f %f %f residual %f\n", m_convergence_based_relaxation, va.x(), va.y(), va.z(), vb.x(), vb.y(), vb.z(),
		//		residualSquare);
	}
	// If the rigid is static, we freeze the soft node.

	return residualSquare;
}

btVector3 btDeformableNodeAnchorConstraint::getVb() const
{
	return m_anchor->m_node->m_v;
}

void btDeformableNodeAnchorConstraint::applyImpulse(const btVector3& impulse)
{
	btVector3 dv = impulse * m_anchor->m_c2;
	m_anchor->m_node->m_v -= dv;
}

btVector3 btDeformableNodeAnchorConstraint::getSplitVa() const

{
	const btSoftBody::sCti& cti = m_anchor->m_cti;
	btVector3 va(0, 0, 0);
	if (cti.m_colObj->hasContactResponse())
	{
		if (cti.m_colObj->getInternalType() == btCollisionObject::CO_RIGID_BODY)
		{
			btRigidBody* rigidCol = (btRigidBody*)btRigidBody::upcast(cti.m_colObj);
			va = rigidCol ? rigidCol->getPushVelocityInLocalPoint(m_anchor->m_c1) : btVector3(0, 0, 0);
		}
		// TODO feather
	}
	return va;
}

btVector3 btDeformableNodeAnchorConstraint::getSplitVb() const
{
	return m_anchor->m_node->m_splitv;
}

void btDeformableNodeAnchorConstraint::applySplitImpulse(const btVector3& impulse)
{
	btVector3 dv = impulse * m_anchor->m_c2;
	m_anchor->m_node->m_splitv -= dv;
}

// This function is to solve the positional drift of the deformable node anchor constraint by applying split impulse. Originally, the drift was
// removed by the n->m_x = a.m_cti.m_colObj->getWorldTransform() * a.m_local; hack in btDeformableBodySolver::applyTransforms. That cause massive destabilisation of the rigid body when the
// soft body stiffness was high, because it caused huge velocity spikes on that soft body node as it was yanked back into place by the elastic forces. These velocity spikes poisoned the rigid body velocity in
// btDeformableNodeAnchorConstraint::solveConstraint.
btScalar btDeformableNodeAnchorConstraint::solveSplitImpulse(const btContactSolverInfo& infoGlobal)
{
	btScalar residualSquare = 0.0;
	if (!m_anchor->m_body->isStaticOrKinematicObject() && m_anchor->m_body->isActive())
	{
		const btSoftBody::sCti& cti = m_anchor->m_cti;
		auto va = getSplitVa() + getVa();
		// Why is m_vn (previous velocity) here? m_v could not be used, because at this stage it contains only the explicit forces (like mouse drag) applied but not the implicit forces like the linear elasticity.
		// To predict the the position of the deformable node, we use vn as an approximation of the velocity at the end of the time step, which includes the effect of both explicit and implicit forces. It is not perfectly
		// accurate, but it works reasonably well due to time coherence. To get perfect velocity, this solve would have to be moved after the implicit forces are calculated. This could be a future TODO - not sure if easy.
		const auto vb = m_anchor->m_node->m_vn + m_anchor->m_node->m_splitv;

		//auto rigid_pt = (cti.m_colObj->getWorldTransform() * m_anchor->m_local);
		//fprintf(stderr, "m_anchor->m_node->m_x %f %f %f cti.m_colObj->getWorldTransform() * m_anchor->m_local %f %f %f diff %f %f %f\n", m_anchor->m_node->m_x.x(), m_anchor->m_node->m_x.y(), m_anchor->m_node->m_x.z(), rigid_pt.x(), rigid_pt.y(), rigid_pt.z(), trigger.x(), trigger.y(), trigger.z());

		// Predicted anchor poisitions for both the rigid and soft are calculated in this part.
		auto anchor_soft_predicted = (m_anchor->m_node->m_x + vb * infoGlobal.m_timeStep);

		btTransform rigid_tr = cti.m_colObj->getWorldTransform();
		btTransform tr_predicted;
		btTransformUtil::integrateTransform(
			rigid_tr,
			m_anchor->m_body->getLinearVelocity() + m_anchor->m_body->getPushVelocity(),
			m_anchor->m_body->getAngularVelocity() + m_anchor->m_body->getTurnVelocity() * infoGlobal.m_splitImpulseTurnErp,
			infoGlobal.m_timeStep,
			tr_predicted);
		btVector3 anchor_rigid_predicted = tr_predicted * m_anchor->m_local;

		btVector3 pos_diff = anchor_soft_predicted - anchor_rigid_predicted;

		residualSquare = btDot(pos_diff, pos_diff);

		// Maybe this early return would prevent some slight overshoot?
		//if (residualSquare < infoGlobal.m_leastSquaresResidualThreshold)
		//	return residualSquare;

		/*fprintf(stderr, "anchor_rigid_predicted %f %f %f anchor_soft_predicted %f %f %f pos_diff %f %f %f va %f %f %f vb %f %f %f residual %f\n",
			anchor_rigid_predicted.x(), anchor_rigid_predicted.y(), anchor_rigid_predicted.z(),
			anchor_soft_predicted.x(), anchor_soft_predicted.y(), anchor_soft_predicted.z(),
			pos_diff.x(), pos_diff.y(), pos_diff.z(),
			va.x(), va.y(), va.z(),
			vb.x(), vb.y(), vb.z(),
			residualSquare);*/

		auto rigid_impulse_fraction = 0.5;
		auto soft_impulse_fraction = 0.5;
		// About this heuristic. When there are penetrations, there is this fight between the solver body push velocity and the regular body push velocity (same applies for turn velocities naturally).
		// This results in bad penetration resolve when the anchored rigid is penetrating. One solution is to make this method be solver body based too. It was actually implemented on commit c619f118e060011efd3892b41f09851253a98f04
		// but there were problems. There was still a fight between solver body push velocities for unstuck and for anchor solve when they were acting against each other. For example the unstuck added positive X direction to the
		// push velocity, but the anchor solve subtracted from the X direction to maintain the anchorage. In worst case, this resulted in explosive growth of push velocities. One workaround would be to remeber the applied
		// push veloities from resolveSplitPenetrationImpulse and also add them to the soft anchor nodes here. This would ease the work for the solver tremendously. It still has to be tried if the solver body based solution
		// becomes unavoidable for some reason. Alternative workaround is the heuristic below. Instead of fighting against the solver body push velocity, we just give up and make the anchor solve be fully soft body based when
		// there are penetrations. This workaround can also be implemented on the solver body based solution, but obsoletes the solver body based solution because then there is no push velocities fight at all.
		if (m_anchor_rigid_penetration)
		{
			rigid_impulse_fraction = 0.0;
			soft_impulse_fraction = 1.0;
		}
		btVector3 rigid_velocity = (pos_diff * rigid_impulse_fraction) / infoGlobal.m_timeStep;
		btVector3 soft_velocity = (pos_diff * soft_impulse_fraction) / infoGlobal.m_timeStep;

		// Applies velocity to the rigid body, to minimize the gap between rigid and soft anchor positions.
		if (cti.m_colObj->getInternalType() == btCollisionObject::CO_RIGID_BODY)
		{
			if (m_anchor->m_body)
			{
				// Not applying real impulses but directly adding velocities significantly reduces the number of solver iterations needed to converge.
				// Remember - we do not need physically corerct impulse application, we just correct a drift.
				// One pitfall here is to reason that we could perhaps use applyCentralPushImpulse here and get away with it, because anchor is only a positional constraint. It does not converge when multiple
				// anchors are involved. Only when rotations are involved, the gap is reduced on the current anchor and at least partially preserved on the previous anchor.
				// UPDATE: If direct velocity manipulation turns out to be problematic for any reason, we can again try with regular impulses, but also with the impulse warmstarting
				// (see the comment about m_appliedImpulse in btDeformableNodeAnchorConstraint::solveConstraint) and similar m_previous_residual logic like in btDeformableNodeAnchorConstraint::solveConstraint

				// Since we are not applying impulses but directly changing velocities, we have to use a purely
				// kinematic way of applying velocity at a point, similar to the way the impulse application does it, but
				// without any inertia tensor or mass usage when updating rotation velocities.

				// Why this works:
				// Identity: if dvPar ⟂ r, then ω = (r×dvPar)/|r|^2 gives ω×r = dvPar.
				// This is a clean kinematic solve
				btVector3 dvPoint = rigid_velocity;
				btVector3 r = m_anchor->m_c1;
				btScalar r2 = r.length2();
				if (r2 > SIMD_EPSILON)
				{
					btVector3 dvPar = r * (dvPoint.dot(r) / r2);
					btVector3 dvPerp = dvPoint - dvPar;  // Remove component parallel to r

					m_anchor->m_body->setPushVelocity(m_anchor->m_body->getPushVelocity() + dvPar);

					// Choose omega so that omega x r = dvPerp
					btVector3 domega = r.cross(dvPerp) / r2;
					m_anchor->m_body->setTurnVelocity(m_anchor->m_body->getTurnVelocity() + domega);
				}
				else
				{
					m_anchor->m_body->setPushVelocity(m_anchor->m_body->getPushVelocity() + dvPoint);
				}

				// I thought that I could do it by only pushing the rigid body, but there are problems with that. If soft and rigid are connected by many anchors, then the soft anchored nodes might get slightly squashed
				// by elasticity and the rigid will then never be able to reach all anchor positions with sufficient accuracy, so it will never converge and will run full iteration count (slow).
				// To counter that, the soft must also contribute by being pushed.
				applySplitImpulse(soft_velocity / m_anchor->m_c2);

				//m_anchor->m_body->applyPushImpulse((rigid_velocity / m_anchor->m_body->getInvMass()) / m_anchor->m_body->getLinearFactor(), m_anchor->m_c1);
				fprintf(stderr, "btDeformableNodeAnchorConstraint::solveSplitImpulse real_push %f %f %f real_turn %f %f %f residual %f\n", m_anchor->m_body->getPushVelocity().x(), m_anchor->m_body->getPushVelocity().y(), m_anchor->m_body->getPushVelocity().z(),
						m_anchor->m_body->getTurnVelocity().x(), m_anchor->m_body->getTurnVelocity().y(), m_anchor->m_body->getTurnVelocity().z(), residualSquare);
				fprintf(stderr, "regular m_anchor->m_node->m_v %f %f %f\n", m_anchor->m_node->m_v.x(), m_anchor->m_node->m_v.y(), m_anchor->m_node->m_v.z());
			}
		}
	}
	// If the rigid is static, we freeze the soft node. In this case there is nothing to iterate to, which makes the normal drift correction unusable in this case.

	// Another pitfall here is to try to perform a final prediction here and return a residual based on that. Not only would that be a redundant code duplication, because that same functionality
	// would be done at the start of the next solveSplitImpulse, but it also causes partial mis-convergence.

	return residualSquare;
}

/* ================   Deformable vs. Rigid   =================== */
btDeformableRigidContactConstraint::btDeformableRigidContactConstraint(const btSoftBody::DeformableRigidContact& c, const btContactSolverInfo& infoGlobal)
	: m_contact(&c), btDeformableContactConstraint(c.m_cti.m_normal, infoGlobal)
{
	m_total_normal_dv.setZero();
	m_total_tangent_dv.setZero();
	// The magnitude of penetration is the depth of penetration.
	m_penetration = c.m_cti.m_offset;
	m_total_split_impulse = 0;
	m_binding = false;
}

btDeformableRigidContactConstraint::btDeformableRigidContactConstraint(const btDeformableRigidContactConstraint& other)
	: m_contact(other.m_contact), btDeformableContactConstraint(other), m_penetration(other.m_penetration), m_total_split_impulse(other.m_total_split_impulse), m_binding(other.m_binding)
{
	m_total_normal_dv = other.m_total_normal_dv;
	m_total_tangent_dv = other.m_total_tangent_dv;
}

btVector3 btDeformableRigidContactConstraint::getVa() const
{
	const btSoftBody::sCti& cti = m_contact->m_cti;
	btVector3 va(0, 0, 0);
	if (cti.m_colObj->hasContactResponse())
	{
		btRigidBody* rigidCol = 0;
		btMultiBodyLinkCollider* multibodyLinkCol = 0;

		// grab the velocity of the rigid body
		if (cti.m_colObj->getInternalType() == btCollisionObject::CO_RIGID_BODY)
		{
			rigidCol = (btRigidBody*)btRigidBody::upcast(cti.m_colObj);
			va = rigidCol ? (rigidCol->getVelocityInLocalPoint(m_contact->m_c1)) : btVector3(0, 0, 0);
		}
		else if (cti.m_colObj->getInternalType() == btCollisionObject::CO_FEATHERSTONE_LINK)
		{
			multibodyLinkCol = (btMultiBodyLinkCollider*)btMultiBodyLinkCollider::upcast(cti.m_colObj);
			if (multibodyLinkCol)
			{
				const int ndof = multibodyLinkCol->m_multiBody->getNumDofs() + 6;
				const btScalar* J_n = &m_contact->jacobianData_normal.m_jacobians[0];
				const btScalar* J_t1 = &m_contact->jacobianData_t1.m_jacobians[0];
				const btScalar* J_t2 = &m_contact->jacobianData_t2.m_jacobians[0];
				const btScalar* local_v = multibodyLinkCol->m_multiBody->getVelocityVector();
				const btScalar* local_dv = multibodyLinkCol->m_multiBody->getDeltaVelocityVector();
				// add in the normal component of the va
				btScalar vel = 0.0;
				for (int k = 0; k < ndof; ++k)
				{
					vel += (local_v[k] + local_dv[k]) * J_n[k];
				}
				va = cti.m_normal * vel;
				// add in the tangential components of the va
				vel = 0.0;
				for (int k = 0; k < ndof; ++k)
				{
					vel += (local_v[k] + local_dv[k]) * J_t1[k];
				}
				va += m_contact->t1 * vel;
				vel = 0.0;
				for (int k = 0; k < ndof; ++k)
				{
					vel += (local_v[k] + local_dv[k]) * J_t2[k];
				}
				va += m_contact->t2 * vel;
			}
		}
	}
	return va;
}

btVector3 btDeformableRigidContactConstraint::getSplitVa() const
{
	const btSoftBody::sCti& cti = m_contact->m_cti;
	btVector3 va(0, 0, 0);
	if (cti.m_colObj->hasContactResponse())
	{
		btRigidBody* rigidCol = 0;
		btMultiBodyLinkCollider* multibodyLinkCol = 0;

		// grab the velocity of the rigid body
		if (cti.m_colObj->getInternalType() == btCollisionObject::CO_RIGID_BODY)
		{
			rigidCol = (btRigidBody*)btRigidBody::upcast(cti.m_colObj);
			va = rigidCol ? (rigidCol->getPushVelocityInLocalPoint(m_contact->m_c1)) : btVector3(0, 0, 0);
		}
		else if (cti.m_colObj->getInternalType() == btCollisionObject::CO_FEATHERSTONE_LINK)
		{
			multibodyLinkCol = (btMultiBodyLinkCollider*)btMultiBodyLinkCollider::upcast(cti.m_colObj);
			if (multibodyLinkCol)
			{
				const int ndof = multibodyLinkCol->m_multiBody->getNumDofs() + 6;
				const btScalar* J_n = &m_contact->jacobianData_normal.m_jacobians[0];
				const btScalar* J_t1 = &m_contact->jacobianData_t1.m_jacobians[0];
				const btScalar* J_t2 = &m_contact->jacobianData_t2.m_jacobians[0];
				const btScalar* local_split_v = multibodyLinkCol->m_multiBody->getSplitVelocityVector();
				// add in the normal component of the va
				btScalar vel = 0.0;
				for (int k = 0; k < ndof; ++k)
				{
					vel += local_split_v[k] * J_n[k];
				}
				va = cti.m_normal * vel;
				// add in the tangential components of the va
				vel = 0.0;
				for (int k = 0; k < ndof; ++k)
				{
					vel += local_split_v[k] * J_t1[k];
				}
				va += m_contact->t1 * vel;
				vel = 0.0;
				for (int k = 0; k < ndof; ++k)
				{
					vel += local_split_v[k] * J_t2[k];
				}
				va += m_contact->t2 * vel;
			}
		}
	}
	return va;
}

btScalar btDeformableRigidContactConstraint::solveConstraint(const btContactSolverInfo& infoGlobal)
{
	const btSoftBody::sCti& cti = m_contact->m_cti;
	btVector3 va = getVa();
	btVector3 vb = getVb();
	btVector3 vr = vb - va;
	btScalar dn = btDot(vr, cti.m_normal) + m_total_normal_dv.dot(cti.m_normal) * infoGlobal.m_deformable_cfm;
	// TODO add some comments trying to describe the m_penetration logic. Maybe there should not be any m_penetration handling if infoGlobal.m_splitImpulse is true? Because then penetrations are completely handled in btDeformableRigidContactConstraint::solveSplitImpulse?
	if (m_penetration < 0)
	{
		dn += m_penetration / infoGlobal.m_timeStep;
	}
	if (!infoGlobal.m_splitImpulse && m_penetration < 0)
	{
		dn += m_penetration * infoGlobal.m_deformable_erp / infoGlobal.m_timeStep;
	}
	// dn is the normal component of velocity difference. Approximates the residual. // todo xuchenhan@: this prob needs to be scaled by dt
	btVector3 impulse = m_contact->m_c0 * (vr +
										   m_total_normal_dv * infoGlobal.m_deformable_cfm +
										   ((m_penetration < 0) ? m_penetration / infoGlobal.m_timeStep * cti.m_normal : btVector3(0, 0, 0)));

	if (!infoGlobal.m_splitImpulse && m_penetration < 0)
	{
		impulse += m_contact->m_c0 * (m_penetration * infoGlobal.m_deformable_erp / infoGlobal.m_timeStep * cti.m_normal);
	}
	btVector3 impulse_normal = m_contact->m_c0 * (cti.m_normal * dn);
	btVector3 impulse_tangent = impulse - impulse_normal;
	if (dn > 0)
	{
		return 0;
	}
	m_binding = true;
	btScalar residualSquare = dn * dn;
	btVector3 old_total_tangent_dv = m_total_tangent_dv;
	// m_c5 is the inverse mass of the deformable node/face
	m_total_normal_dv -= m_contact->m_c5 * impulse_normal;
	m_total_tangent_dv -= m_contact->m_c5 * impulse_tangent;

	if (m_total_normal_dv.dot(cti.m_normal) < 0)
	{
		// separating in the normal direction
		m_binding = false;
		m_static = false;
		impulse_tangent.setZero();
	}
	else
	{
		if (m_total_normal_dv.norm() * m_contact->m_c3 < m_total_tangent_dv.norm())
		{
			// dynamic friction
			// with dynamic friction, the impulse are still applied to the two objects colliding, however, it does not pose a constraint in the cg solve, hence the change to dv merely serves to update velocity in the contact iterations.
			m_static = false;
			if (m_total_tangent_dv.safeNorm() < SIMD_EPSILON)
			{
				m_total_tangent_dv = btVector3(0, 0, 0);
			}
			else
			{
				m_total_tangent_dv = m_total_tangent_dv.normalized() * m_total_normal_dv.safeNorm() * m_contact->m_c3;
			}
			//            impulse_tangent = -btScalar(1)/m_contact->m_c2 * (m_total_tangent_dv - old_total_tangent_dv);
			impulse_tangent = m_contact->m_c5.inverse() * (old_total_tangent_dv - m_total_tangent_dv);
		}
		else
		{
			// static friction
			m_static = true;
		}
	}
	impulse = impulse_normal + impulse_tangent;
	// apply impulse to deformable nodes involved and change their velocities
	//printf("impulse %f %f %f\n", impulse.x(), impulse.y(), impulse.z());
	applyImpulse(impulse);

	if (cti.m_contact_point_impulse_magnitude)
		*cti.m_contact_point_impulse_magnitude = impulse.length();

	// apply impulse to the rigid/multibodies involved and change their velocities
	if (cti.m_colObj->getInternalType() == btCollisionObject::CO_RIGID_BODY)
	{
		btRigidBody* rigidCol = (btRigidBody*)btRigidBody::upcast(cti.m_colObj);
		if (rigidCol)
		{
			rigidCol->applyImpulse(impulse, m_contact->m_c1);
		}
	}
	else if (cti.m_colObj->getInternalType() == btCollisionObject::CO_FEATHERSTONE_LINK)
	{
		btMultiBodyLinkCollider* multibodyLinkCol = (btMultiBodyLinkCollider*)btMultiBodyLinkCollider::upcast(cti.m_colObj);
		if (multibodyLinkCol)
		{
			const btScalar* deltaV_normal = &m_contact->jacobianData_normal.m_deltaVelocitiesUnitImpulse[0];
			// apply normal component of the impulse
			multibodyLinkCol->m_multiBody->applyDeltaVeeMultiDof2(deltaV_normal, impulse.dot(cti.m_normal));
			if (impulse_tangent.norm() > SIMD_EPSILON)
			{
				// apply tangential component of the impulse
				const btScalar* deltaV_t1 = &m_contact->jacobianData_t1.m_deltaVelocitiesUnitImpulse[0];
				multibodyLinkCol->m_multiBody->applyDeltaVeeMultiDof2(deltaV_t1, impulse.dot(m_contact->t1));
				const btScalar* deltaV_t2 = &m_contact->jacobianData_t2.m_deltaVelocitiesUnitImpulse[0];
				multibodyLinkCol->m_multiBody->applyDeltaVeeMultiDof2(deltaV_t2, impulse.dot(m_contact->t2));
			}
		}
	}
	return residualSquare;
}

btScalar btDeformableRigidContactConstraint::solveSplitImpulse(const btContactSolverInfo& infoGlobal)
{
	btScalar MAX_PENETRATION_CORRECTION = infoGlobal.m_deformable_maxErrorReduction;
	const btSoftBody::sCti& cti = m_contact->m_cti;
	btVector3 vb = getSplitVb();
	btVector3 va = getSplitVa();
	btScalar p = m_penetration;
	if (p > 0)
	{
		return 0;
	}
	btVector3 vr = vb - va;
	btScalar dn = btDot(vr, cti.m_normal) + p * infoGlobal.m_deformable_erp / infoGlobal.m_timeStep;
	if (dn > 0)
	{
		return 0;
	}
	if (m_total_split_impulse + dn > MAX_PENETRATION_CORRECTION)
	{
		dn = MAX_PENETRATION_CORRECTION - m_total_split_impulse;
	}
	if (m_total_split_impulse + dn < -MAX_PENETRATION_CORRECTION)
	{
		dn = -MAX_PENETRATION_CORRECTION - m_total_split_impulse;
	}
	m_total_split_impulse += dn;

	btScalar residualSquare = dn * dn;
	const btVector3 impulse = m_contact->m_c0 * (cti.m_normal * dn);
	applySplitImpulse(impulse);

	// apply split impulse to the rigid/multibodies involved and change their velocities
	if (cti.m_colObj->getInternalType() == btCollisionObject::CO_RIGID_BODY)
	{
		btRigidBody* rigidCol = 0;
		rigidCol = (btRigidBody*)btRigidBody::upcast(cti.m_colObj);
		if (rigidCol)
		{
			rigidCol->applyPushImpulse(impulse, m_contact->m_c1);
		}
	}
	else if (cti.m_colObj->getInternalType() == btCollisionObject::CO_FEATHERSTONE_LINK)
	{
		btMultiBodyLinkCollider* multibodyLinkCol = 0;
		multibodyLinkCol = (btMultiBodyLinkCollider*)btMultiBodyLinkCollider::upcast(cti.m_colObj);
		if (multibodyLinkCol)
		{
			const btScalar* deltaV_normal = &m_contact->jacobianData_normal.m_deltaVelocitiesUnitImpulse[0];
			// apply normal component of the impulse
			multibodyLinkCol->m_multiBody->applyDeltaSplitVeeMultiDof(deltaV_normal, impulse.dot(cti.m_normal));
		}
	}
	return residualSquare;
}
/* ================   Node vs. Rigid   =================== */
btDeformableNodeRigidContactConstraint::btDeformableNodeRigidContactConstraint(const btSoftBody::DeformableNodeRigidContact& contact, const btContactSolverInfo& infoGlobal)
	: m_node(contact.m_node), btDeformableRigidContactConstraint(contact, infoGlobal)
{
}

btDeformableNodeRigidContactConstraint::btDeformableNodeRigidContactConstraint(const btDeformableNodeRigidContactConstraint& other)
	: m_node(other.m_node), btDeformableRigidContactConstraint(other)
{
}

btVector3 btDeformableNodeRigidContactConstraint::getVb() const
{
	return m_node->m_v;
}

btVector3 btDeformableNodeRigidContactConstraint::getSplitVb() const
{
	return m_node->m_splitv;
}

btVector3 btDeformableNodeRigidContactConstraint::getDv(const btSoftBody::Node* node) const
{
	return m_total_normal_dv + m_total_tangent_dv;
}

void btDeformableNodeRigidContactConstraint::applyImpulse(const btVector3& impulse)
{
	const btSoftBody::DeformableNodeRigidContact* contact = getContact();
	btVector3 dv = contact->m_c5 * impulse;
	contact->m_node->m_v -= dv;
}

void btDeformableNodeRigidContactConstraint::applySplitImpulse(const btVector3& impulse)
{
	const btSoftBody::DeformableNodeRigidContact* contact = getContact();
	btVector3 dv = contact->m_c5 * impulse;
	contact->m_node->m_splitv -= dv;
}

/* ================   Face vs. Rigid   =================== */
btDeformableFaceRigidContactConstraint::btDeformableFaceRigidContactConstraint(const btSoftBody::DeformableFaceRigidContact& contact, const btContactSolverInfo& infoGlobal, bool useStrainLimiting)
	: m_face(contact.m_face), m_useStrainLimiting(useStrainLimiting), btDeformableRigidContactConstraint(contact, infoGlobal)
{
}

btDeformableFaceRigidContactConstraint::btDeformableFaceRigidContactConstraint(const btDeformableFaceRigidContactConstraint& other)
	: m_face(other.m_face), m_useStrainLimiting(other.m_useStrainLimiting), btDeformableRigidContactConstraint(other)
{
}

btVector3 btDeformableFaceRigidContactConstraint::getVb() const
{
	const btSoftBody::DeformableFaceRigidContact* contact = getContact();
	btVector3 vb = m_face->m_n[0]->m_v * contact->m_bary[0] + m_face->m_n[1]->m_v * contact->m_bary[1] + m_face->m_n[2]->m_v * contact->m_bary[2];
	return vb;
}

btVector3 btDeformableFaceRigidContactConstraint::getDv(const btSoftBody::Node* node) const
{
	btVector3 face_dv = m_total_normal_dv + m_total_tangent_dv;
	const btSoftBody::DeformableFaceRigidContact* contact = getContact();
	if (m_face->m_n[0] == node)
	{
		return face_dv * contact->m_weights[0];
	}
	if (m_face->m_n[1] == node)
	{
		return face_dv * contact->m_weights[1];
	}
	btAssert(node == m_face->m_n[2]);
	return face_dv * contact->m_weights[2];
}

void btDeformableFaceRigidContactConstraint::applyImpulse(const btVector3& impulse)
{
	const btSoftBody::DeformableFaceRigidContact* contact = getContact();
	btVector3 dv = impulse * contact->m_c2;
	btSoftBody::Face* face = contact->m_face;

	// save applied impulse
	contact->m_cti.m_impulse = impulse;

	btVector3& v0 = face->m_n[0]->m_v;
	btVector3& v1 = face->m_n[1]->m_v;
	btVector3& v2 = face->m_n[2]->m_v;
	const btScalar& im0 = face->m_n[0]->m_im;
	const btScalar& im1 = face->m_n[1]->m_im;
	const btScalar& im2 = face->m_n[2]->m_im;
	if (im0 > 0)
		v0 -= dv * contact->m_weights[0];
	if (im1 > 0)
		v1 -= dv * contact->m_weights[1];
	if (im2 > 0)
		v2 -= dv * contact->m_weights[2];
	if (m_useStrainLimiting)
	{
		btScalar relaxation = 1. / btScalar(m_infoGlobal->m_numIterations);
		btScalar m01 = (relaxation / (im0 + im1));
		btScalar m02 = (relaxation / (im0 + im2));
		btScalar m12 = (relaxation / (im1 + im2));
#ifdef USE_STRAIN_RATE_LIMITING
		// apply strain limiting to prevent the new velocity to change the current length of the edge by more than 1%.
		btScalar p = 0.01;
		btVector3& x0 = face->m_n[0]->m_x;
		btVector3& x1 = face->m_n[1]->m_x;
		btVector3& x2 = face->m_n[2]->m_x;
		const btVector3 x_diff[3] = {x1 - x0, x2 - x0, x2 - x1};
		const btVector3 v_diff[3] = {v1 - v0, v2 - v0, v2 - v1};
		btVector3 u[3];
		btScalar x_diff_dot_u, dn[3];
		btScalar dt = m_infoGlobal->m_timeStep;
		for (int i = 0; i < 3; ++i)
		{
			btScalar x_diff_norm = x_diff[i].safeNorm();
			btScalar x_diff_norm_new = (x_diff[i] + v_diff[i] * dt).safeNorm();
			btScalar strainRate = x_diff_norm_new / x_diff_norm;
			u[i] = v_diff[i];
			u[i].safeNormalize();
			if (x_diff_norm == 0 || (1 - p <= strainRate && strainRate <= 1 + p))
			{
				dn[i] = 0;
				continue;
			}
			x_diff_dot_u = btDot(x_diff[i], u[i]);
			btScalar s;
			if (1 - p > strainRate)
			{
				s = 1 / dt * (-x_diff_dot_u - btSqrt(x_diff_dot_u * x_diff_dot_u + (p * p - 2 * p) * x_diff_norm * x_diff_norm));
			}
			else
			{
				s = 1 / dt * (-x_diff_dot_u + btSqrt(x_diff_dot_u * x_diff_dot_u + (p * p + 2 * p) * x_diff_norm * x_diff_norm));
			}
			//		x_diff_norm_new = (x_diff[i] + s * u[i] * dt).safeNorm();
			//		strainRate = x_diff_norm_new/x_diff_norm;
			dn[i] = s - v_diff[i].safeNorm();
		}
		btVector3 dv0 = im0 * (m01 * u[0] * (-dn[0]) + m02 * u[1] * -(dn[1]));
		btVector3 dv1 = im1 * (m01 * u[0] * (dn[0]) + m12 * u[2] * (-dn[2]));
		btVector3 dv2 = im2 * (m12 * u[2] * (dn[2]) + m02 * u[1] * (dn[1]));
#else
		// apply strain limiting to prevent undamped modes
		btVector3 dv0 = im0 * (m01 * (v1 - v0) + m02 * (v2 - v0));
		btVector3 dv1 = im1 * (m01 * (v0 - v1) + m12 * (v2 - v1));
		btVector3 dv2 = im2 * (m12 * (v1 - v2) + m02 * (v0 - v2));
#endif
		v0 += dv0;
		v1 += dv1;
		v2 += dv2;
	}
}

btVector3 btDeformableFaceRigidContactConstraint::getSplitVb() const
{
	const btSoftBody::DeformableFaceRigidContact* contact = getContact();
	btVector3 vb = (m_face->m_n[0]->m_splitv) * contact->m_bary[0] + (m_face->m_n[1]->m_splitv) * contact->m_bary[1] + (m_face->m_n[2]->m_splitv) * contact->m_bary[2];
	return vb;
}

void btDeformableFaceRigidContactConstraint::applySplitImpulse(const btVector3& impulse)
{
	const btSoftBody::DeformableFaceRigidContact* contact = getContact();
	btVector3 dv = impulse * contact->m_c2;
	btSoftBody::Face* face = contact->m_face;
	btVector3& v0 = face->m_n[0]->m_splitv;
	btVector3& v1 = face->m_n[1]->m_splitv;
	btVector3& v2 = face->m_n[2]->m_splitv;
	const btScalar& im0 = face->m_n[0]->m_im;
	const btScalar& im1 = face->m_n[1]->m_im;
	const btScalar& im2 = face->m_n[2]->m_im;
	if (im0 > 0)
	{
		v0 -= dv * contact->m_weights[0];
	}
	if (im1 > 0)
	{
		v1 -= dv * contact->m_weights[1];
	}
	if (im2 > 0)
	{
		v2 -= dv * contact->m_weights[2];
	}
}

/* ================   Face vs. Node   =================== */
btDeformableFaceNodeContactConstraint::btDeformableFaceNodeContactConstraint(const btSoftBody::DeformableFaceNodeContact& contact, const btContactSolverInfo& infoGlobal)
	: m_node(contact.m_node), m_face(contact.m_face), m_contact(&contact), btDeformableContactConstraint(contact.m_normal, infoGlobal)
{
	m_total_normal_dv.setZero();
	m_total_tangent_dv.setZero();
}

btVector3 btDeformableFaceNodeContactConstraint::getVa() const
{
	return m_node->m_v;
}

btVector3 btDeformableFaceNodeContactConstraint::getVb() const
{
	const btSoftBody::DeformableFaceNodeContact* contact = getContact();
	btVector3 vb = m_face->m_n[0]->m_v * contact->m_bary[0] + m_face->m_n[1]->m_v * contact->m_bary[1] + m_face->m_n[2]->m_v * contact->m_bary[2];
	return vb;
}

btVector3 btDeformableFaceNodeContactConstraint::getDv(const btSoftBody::Node* n) const
{
	btVector3 dv = m_total_normal_dv + m_total_tangent_dv;
	if (n == m_node)
		return dv;
	const btSoftBody::DeformableFaceNodeContact* contact = getContact();
	if (m_face->m_n[0] == n)
	{
		return dv * contact->m_weights[0];
	}
	if (m_face->m_n[1] == n)
	{
		return dv * contact->m_weights[1];
	}
	btAssert(n == m_face->m_n[2]);
	return dv * contact->m_weights[2];
}

btScalar btDeformableFaceNodeContactConstraint::solveConstraint(const btContactSolverInfo& infoGlobal)
{
	btVector3 va = getVa();
	btVector3 vb = getVb();
	btVector3 vr = vb - va;
	const btScalar dn = btDot(vr, m_contact->m_normal);
	// dn is the normal component of velocity diffrerence. Approximates the residual. // todo xuchenhan@: this prob needs to be scaled by dt
	btScalar residualSquare = dn * dn;
	btVector3 impulse = m_contact->m_c0 * vr;
	const btVector3 impulse_normal = m_contact->m_c0 * (m_contact->m_normal * dn);
	btVector3 impulse_tangent = impulse - impulse_normal;

	btVector3 old_total_tangent_dv = m_total_tangent_dv;
	// m_c2 is the inverse mass of the deformable node/face
	if (m_node->m_im > 0)
	{
		m_total_normal_dv -= impulse_normal * m_node->m_im;
		m_total_tangent_dv -= impulse_tangent * m_node->m_im;
	}
	else
	{
		m_total_normal_dv -= impulse_normal * m_contact->m_imf;
		m_total_tangent_dv -= impulse_tangent * m_contact->m_imf;
	}

	if (m_total_normal_dv.dot(m_contact->m_normal) > 0)
	{
		// separating in the normal direction
		m_static = false;
		m_total_tangent_dv = btVector3(0, 0, 0);
		impulse_tangent.setZero();
	}
	else
	{
		if (m_total_normal_dv.norm() * m_contact->m_friction < m_total_tangent_dv.norm())
		{
			// dynamic friction
			// with dynamic friction, the impulse are still applied to the two objects colliding, however, it does not pose a constraint in the cg solve, hence the change to dv merely serves to update velocity in the contact iterations.
			m_static = false;
			if (m_total_tangent_dv.safeNorm() < SIMD_EPSILON)
			{
				m_total_tangent_dv = btVector3(0, 0, 0);
			}
			else
			{
				m_total_tangent_dv = m_total_tangent_dv.normalized() * m_total_normal_dv.safeNorm() * m_contact->m_friction;
			}
			impulse_tangent = -btScalar(1) / m_node->m_im * (m_total_tangent_dv - old_total_tangent_dv);
		}
		else
		{
			// static friction
			m_static = true;
		}
	}
	impulse = impulse_normal + impulse_tangent;
	// apply impulse to deformable nodes involved and change their velocities
	applyImpulse(impulse);
	return residualSquare;
}

void btDeformableFaceNodeContactConstraint::applyImpulse(const btVector3& impulse)
{
	const btSoftBody::DeformableFaceNodeContact* contact = getContact();
	btVector3 dva = impulse * contact->m_node->m_im;
	btVector3 dvb = impulse * contact->m_imf;
	if (contact->m_node->m_im > 0)
	{
		contact->m_node->m_v += dva;
	}

	btSoftBody::Face* face = contact->m_face;
	btVector3& v0 = face->m_n[0]->m_v;
	btVector3& v1 = face->m_n[1]->m_v;
	btVector3& v2 = face->m_n[2]->m_v;
	const btScalar& im0 = face->m_n[0]->m_im;
	const btScalar& im1 = face->m_n[1]->m_im;
	const btScalar& im2 = face->m_n[2]->m_im;
	if (im0 > 0)
	{
		v0 -= dvb * contact->m_weights[0];
	}
	if (im1 > 0)
	{
		v1 -= dvb * contact->m_weights[1];
	}
	if (im2 > 0)
	{
		v2 -= dvb * contact->m_weights[2];
	}
}
