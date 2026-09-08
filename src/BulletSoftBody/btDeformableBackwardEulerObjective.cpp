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

#include "btDeformableBackwardEulerObjective.h"
#include "btPreconditioner.h"
#include "LinearMath/btQuickprof.h"

btDeformableBackwardEulerObjective::btDeformableBackwardEulerObjective(btAlignedObjectArray<btSoftBody*>& softBodies, const TVStack& backup_v)
	: m_softBodies(softBodies), m_projection(softBodies), m_backupVelocity(backup_v), m_implicit(false)
{
	m_massPreconditioner = new MassPreconditioner(m_softBodies);
	m_KKTPreconditioner = new KKTPreconditioner(m_softBodies, m_projection, m_lf, m_dt, m_implicit);
	m_preconditioner = m_KKTPreconditioner;
}

btDeformableBackwardEulerObjective::~btDeformableBackwardEulerObjective()
{
	delete m_KKTPreconditioner;
	delete m_massPreconditioner;
}

void btDeformableBackwardEulerObjective::reinitialize(bool nodeUpdated, btScalar dt)
{
	BT_PROFILE("reinitialize");
	if (dt > 0)
	{
		setDt(dt);
	}
	if (nodeUpdated)
	{
		updateId();
	}
	for (int i = 0; i < m_lf.size(); ++i)
	{
		m_lf[i]->reinitialize(nodeUpdated);
	}
	btMatrix3x3 I;
	I.setIdentity();
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		for (int j = 0; j < psb->m_nodes.size(); ++j)
		{
			if (psb->m_nodes[j].m_frozen <= 0 && psb->m_nodes[j].m_im > 0)
				psb->m_nodes[j].m_effectiveMass = I * (1.0 / psb->m_nodes[j].m_im);
		}
	}
	m_projection.reinitialize(nodeUpdated);
	//    m_preconditioner->reinitialize(nodeUpdated);
}

void btDeformableBackwardEulerObjective::setDt(btScalar dt)
{
	m_dt = dt;
}

void btDeformableBackwardEulerObjective::multiply(const TVStack& x, TVStack& b) const
{
	BT_PROFILE("multiply");
	// add in the mass term
	size_t counter = 0;
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		for (int j = 0; j < psb->m_nodes.size(); ++j)
		{
			const btSoftBody::Node& node = psb->m_nodes[j];
			b[counter] = (node.m_frozen > 0) ? btVector3(0, 0, 0) : x[counter] / node.m_im;
			++counter;
		}
	}

	for (int i = 0; i < m_lf.size(); ++i)
	{
		// add damping matrix
		m_lf[i]->addScaledDampingForceDifferential(-m_dt, x, b);
		// Always integrate picking force implicitly for stability.
		if (m_implicit || m_lf[i]->getForceType() == BT_MOUSE_PICKING_FORCE)
		{
			m_lf[i]->addScaledElasticForceDifferential(-m_dt * m_dt, x, b);
		}
	}
	int offset = m_nodes.size();
	for (int i = offset; i < b.size(); ++i)
	{
		b[i].setZero();
	}
	// add in the lagrange multiplier terms

	for (int c = 0; c < m_projection.m_lagrangeMultipliers.size(); ++c)
	{
		// C^T * lambda
		const LagrangeMultiplier& lm = m_projection.m_lagrangeMultipliers[c];
		for (int i = 0; i < lm.m_num_nodes; ++i)
		{
			for (int j = 0; j < lm.m_num_constraints; ++j)
			{
				b[lm.m_indices[i]] += x[offset + c][j] * lm.m_weights[i] * lm.m_dirs[j];
			}
		}
		// C * x
		for (int d = 0; d < lm.m_num_constraints; ++d)
		{
			for (int i = 0; i < lm.m_num_nodes; ++i)
			{
				b[offset + c][d] += lm.m_weights[i] * x[lm.m_indices[i]].dot(lm.m_dirs[d]);
			}
		}
	}
}

void btDeformableBackwardEulerObjective::updateVelocity(const TVStack& dv)
{
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		for (int j = 0; j < psb->m_nodes.size(); ++j)
		{
			btSoftBody::Node& node = psb->m_nodes[j];
			node.m_v = m_backupVelocity[node.index] + dv[node.index];
		}
	}
}

void btDeformableBackwardEulerObjective::applyForce(TVStack& force, bool setZero)
{
	size_t counter = 0;
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		if (!psb->isActive() || psb->isStaticObject())
		{
			counter += psb->m_nodes.size();
			continue;
		}
		if (m_implicit)
		{
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				if (psb->m_nodes[j].m_frozen <= 0 && psb->m_nodes[j].m_im != 0)
				{
					psb->m_nodes[j].m_v += psb->m_nodes[j].m_effectiveMass_inv * force[counter++];
				}
			}
		}
		else
		{
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				btScalar one_over_mass = (psb->m_nodes[j].m_frozen > 0) ? 0 : psb->m_nodes[j].m_im;
				psb->m_nodes[j].m_v += one_over_mass * force[counter++];
			}
		}
	}
	if (setZero)
	{
		for (int i = 0; i < force.size(); ++i)
			force[i].setZero();
	}
}

void btDeformableBackwardEulerObjective::computeResidual(btScalar dt, TVStack& residual)
{
	BT_PROFILE("computeResidual");
	// add implicit force
	for (int i = 0; i < m_lf.size(); ++i)
	{
		// Always integrate picking force implicitly for stability.
		if (m_implicit || m_lf[i]->getForceType() == BT_MOUSE_PICKING_FORCE)
		{
			m_lf[i]->addScaledForces(dt, residual);
		}
		else
		{
			m_lf[i]->addScaledDampingForce(dt, residual);
		}
	}
	//    m_projection.project(residual);
}

btScalar btDeformableBackwardEulerObjective::computeNorm(const TVStack& residual) const
{
	btScalar mag = 0;
	for (int i = 0; i < residual.size(); ++i)
	{
		mag += residual[i].length2();
	}
	return std::sqrt(mag);
}

btScalar btDeformableBackwardEulerObjective::totalEnergy(btScalar dt)
{
	btScalar e = 0;
	for (int i = 0; i < m_lf.size(); ++i)
	{
		e += m_lf[i]->totalEnergy(dt);
	}
	return e;
}

void btDeformableBackwardEulerObjective::applyExplicitForce(TVStack& force)
{
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		m_softBodies[i]->advanceDeformation();
	}
	if (m_implicit)
	{
		// apply forces except gravity force
		btVector3 gravity;
		for (int i = 0; i < m_lf.size(); ++i)
		{
			if (m_lf[i]->getForceType() == BT_GRAVITY_FORCE)
			{
				gravity = static_cast<btDeformableGravityForce*>(m_lf[i])->m_gravity;
			}
			else
			{
				m_lf[i]->addScaledForces(m_dt, force);
			}
		}
		for (int i = 0; i < m_lf.size(); ++i)
		{
			m_lf[i]->addScaledHessian(m_dt);
		}
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (psb->isActive() && !psb->isStaticObject())
			{
				for (int j = 0; j < psb->m_nodes.size(); ++j)
				{
					btSoftBody::Node& node = psb->m_nodes[j];
					// Gravity is an acceleration, so skip frozen or zero-mass nodes
					// but do not scale by mass here.
					if (node.m_frozen <= 0 && node.m_im > 0)
					{
						node.m_v += m_dt * psb->m_gravityFactor * gravity;
					}
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < m_lf.size(); ++i)
		{
			m_lf[i]->addScaledExplicitForce(m_dt, force);
		}
	}
	// calculate inverse mass matrix for all nodes
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		if (psb->isActive() && !psb->isStaticObject())
		{
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				if (psb->m_nodes[j].m_frozen <= 0 && psb->m_nodes[j].m_im > 0)
				{
					psb->m_nodes[j].m_effectiveMass_inv = psb->m_nodes[j].m_effectiveMass.inverse();
				}
			}
		}
	}
	applyForce(force, true);
}

void btDeformableBackwardEulerObjective::initialGuess(TVStack& dv, const TVStack& residual)
{
	size_t counter = 0;
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		for (int j = 0; j < psb->m_nodes.size(); ++j)
		{
			dv[counter] = psb->m_nodes[j].m_im * residual[counter];
			++counter;
		}
	}
}

//set constraints as projections
void btDeformableBackwardEulerObjective::setConstraints(const btContactSolverInfo& infoGlobal)
{
	m_projection.setConstraints(infoGlobal);
}

void btDeformableBackwardEulerObjective::applyDynamicFriction(TVStack& r)
{
	m_projection.applyDynamicFriction(r);
}

bool btDeformableBackwardEulerObjective::setupTranslationCorrection()
{
	m_translationCorrection = false;
	m_translationBodies.clear();
	if (!m_implicit) return false;
	// These validated force operators act independently on each body. This
	// permits three shared A*Z products, rather than three per eligible body.
	for (int f = 0; f < m_lf.size(); ++f)
		if (m_lf[f]->getForceType() != BT_LINEAR_ELASTICITY_FORCE &&
			m_lf[f]->getForceType() != BT_GRAVITY_FORCE &&
			m_lf[f]->getForceType() != BT_NODAL_FORCE) return false;
	btAlignedObjectArray<int> constrained;
	constrained.resize(m_nodes.size(), 0);
	for (int c = 0; c < m_projection.m_lagrangeMultipliers.size(); ++c)
	{
		const LagrangeMultiplier& lm = m_projection.m_lagrangeMultipliers[c];
		for (int n = 0; n < lm.m_num_nodes; ++n) constrained[lm.m_indices[n]] = 1;
	}
	int offset = 0;
	for (int b = 0; b < m_softBodies.size(); ++b)
	{
		const btSoftBody& body = *m_softBodies[b];
		bool eligible = body.isActive() && !body.isStaticObject() && body.m_nodes.size() > 0;
		for (int n = 0; eligible && n < body.m_nodes.size(); ++n)
			eligible = body.m_nodes[n].m_frozen <= 0 && body.m_nodes[n].m_im > 0 && !constrained[offset + n];
		if (eligible)
		{
			TranslationBody entry;
			entry.offset = offset; entry.count = body.m_nodes.size();
			entry.inverse.setIdentity(); entry.coarse.setZero();
			m_translationBodies.push_back(entry);
		}
		offset += body.m_nodes.size();
	}
	if (m_translationBodies.size() == 0) return false;
	m_translationWork.resize(m_nodes.size() + m_projection.m_lagrangeMultipliers.size());
	for (int d = 0; d < 3; ++d)
	{
		for (int n = 0; n < m_translationWork.size(); ++n) m_translationWork[n].setZero();
		for (int b = 0; b < m_translationBodies.size(); ++b)
		{
			const TranslationBody& entry = m_translationBodies[b];
			for (int n = entry.offset; n < entry.offset + entry.count; ++n) m_translationWork[n][d] = 1;
		}
		m_translationAZ[d].resize(m_translationWork.size());
		multiply(m_translationWork, m_translationAZ[d]);
	}
	for (int b = m_translationBodies.size() - 1; b >= 0; --b)
	{
		TranslationBody& entry = m_translationBodies[b];
		btMatrix3x3 coarse;
		for (int d = 0; d < 3; ++d)
		{
			btVector3 sum(0, 0, 0);
			for (int n = entry.offset; n < entry.offset + entry.count; ++n) sum += m_translationAZ[d][n];
			for (int r = 0; r < 3; ++r) coarse[r][d] = sum[r];
		}
		bool regularized = false;
		// Do not substitute a regularized balance equation for this body.
		if (!KKTPreconditioner::invertPositiveBlock(coarse, entry.inverse, regularized) || regularized)
			m_translationBodies.removeAtIndex(b);
	}
	m_translationCorrection = m_translationBodies.size() > 0;
	return m_translationCorrection;
}

void btDeformableBackwardEulerObjective::precondition(const TVStack& x, TVStack& b)
{
	if (!m_translationCorrection) { m_preconditioner->operator()(x, b); return; }
	// B = Q + (I-Q*A)*D*(I-A*Q), with three coarse modes per free body.
	// Constrained bodies and multiplier entries keep the original D action.
	m_translationWork = x;
	for (int body = 0; body < m_translationBodies.size(); ++body)
	{
		TranslationBody& entry = m_translationBodies[body];
		btVector3 sum(0, 0, 0);
		for (int n = entry.offset; n < entry.offset + entry.count; ++n) sum += x[n];
		entry.coarse = entry.inverse * sum;
		for (int n = entry.offset; n < entry.offset + entry.count; ++n)
			for (int d = 0; d < 3; ++d) m_translationWork[n] -= m_translationAZ[d][n] * entry.coarse[d];
	}
	m_preconditioner->operator()(m_translationWork, b);
	for (int body = 0; body < m_translationBodies.size(); ++body)
	{
		const TranslationBody& entry = m_translationBodies[body];
		btVector3 coupling(0, 0, 0);
		for (int n = entry.offset; n < entry.offset + entry.count; ++n)
			for (int d = 0; d < 3; ++d) coupling[d] += m_translationAZ[d][n].dot(b[n]);
		const btVector3 translation = entry.coarse - entry.inverse * coupling;
		for (int n = entry.offset; n < entry.offset + entry.count; ++n) b[n] += translation;
	}
}

btScalar btDeformableBackwardEulerObjective::correctTranslation(TVStack& x, const TVStack& rhs)
{
	if (!m_translationCorrection) return 0;
	multiply(x, m_translationWork);
	btScalar largestCorrection = 0;
	for (int b = 0; b < m_translationBodies.size(); ++b)
	{
		const TranslationBody& entry = m_translationBodies[b];
		btVector3 sum(0, 0, 0);
		for (int n = entry.offset; n < entry.offset + entry.count; ++n) sum += rhs[n] - m_translationWork[n];
		const btVector3 correction = entry.inverse * sum;
		for (int n = entry.offset; n < entry.offset + entry.count; ++n) x[n] += correction;
		largestCorrection = btMax(largestCorrection, correction.length());
	}
	return largestCorrection;
}
