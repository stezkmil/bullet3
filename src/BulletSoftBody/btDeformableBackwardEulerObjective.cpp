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
#include "btTranslationInputKernel.h"
#include "btTranslationCouplingKernel.h"
#include "btMassKernel.h"
#include "btDeformableOptimizationConfig.h"
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
	{
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
#if (BT_DEFORMABLE_OPTIMIZATION_MASK & 16)
		const int count = psb->m_nodes.size();
		if (count > 0)
			btApplyNodeMass(count, &psb->m_nodes[0], &x[counter], &b[counter]);
		counter += count;
#else
		for (int j = 0; j < psb->m_nodes.size(); ++j)
		{
			const btSoftBody::Node& node = psb->m_nodes[j];
			b[counter] = (node.m_frozen > 0) ? btVector3(0, 0, 0) : x[counter] / node.m_im;
			++counter;
		}
#endif
	}

	}
	{
	for (int i = 0; i < m_lf.size(); ++i)
	{
		if (m_implicit)
		{
			m_lf[i]->addImplicitForceDifferential(m_dt, x, b);
		}
		else
		{
			m_lf[i]->addScaledDampingForceDifferential(-m_dt, x, b);
			// Always integrate picking force implicitly for stability.
			if (m_lf[i]->getForceType() == BT_MOUSE_PICKING_FORCE)
			{
				m_lf[i]->addScaledElasticForceDifferential(-m_dt * m_dt, x, b);
			}
		}
	}
	}
	{
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
		btVector3 gravity(0, 0, 0);
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

namespace
{
bool invertRigidCoarse(const btScalar matrix[6][6], btScalar inverse[6][6])
{
	btScalar scale[6], lower[6][6] = {};
	for (int i = 0; i < 6; ++i)
	{
		if (!(matrix[i][i] > 0) || !std::isfinite(double(matrix[i][i]))) return false;
		scale[i] = btSqrt(matrix[i][i]);
	}
	for (int i = 0; i < 6; ++i) for (int j = 0; j <= i; ++j)
	{
		btScalar value = ((matrix[i][j] + matrix[j][i]) * btScalar(.5)) / (scale[i] * scale[j]);
		for (int k = 0; k < j; ++k) value -= lower[i][k] * lower[j][k];
		if (!std::isfinite(double(value))) return false;
		if (i == j)
		{
			// Dependent rotation modes (e.g. collinear bodies) keep translation-only correction.
			if (value <= btScalar(256) * SIMD_EPSILON) return false;
			lower[i][j] = btSqrt(value);
		}
		else lower[i][j] = value / lower[j][j];
	}
	for (int col = 0; col < 6; ++col)
	{
		btScalar y[6], x[6];
		for (int i = 0; i < 6; ++i)
		{
			btScalar value = i == col ? btScalar(1) : btScalar(0);
			for (int k = 0; k < i; ++k) value -= lower[i][k] * y[k];
			y[i] = value / lower[i][i];
		}
		for (int i = 5; i >= 0; --i)
		{
			btScalar value = y[i];
			for (int k = i + 1; k < 6; ++k) value -= lower[k][i] * x[k];
			x[i] = value / lower[i][i];
			inverse[i][col] = x[i] / (scale[i] * scale[col]);
			if (!std::isfinite(double(inverse[i][col]))) return false;
		}
	}
	return true;
}
void applyRigidInverse(const btScalar inverse[6][6], const btScalar input[6], btScalar output[6])
{
	for (int i = 0; i < 6; ++i)
	{
		output[i] = 0;
		for (int j = 0; j < 6; ++j) output[i] += inverse[i][j] * input[j];
	}
}
}

bool btDeformableBackwardEulerObjective::setupContactCoarse()
{
	// Contact couples bodies: keep separate A*Z columns and one joint coarse solve.
	// Bound dense storage for worlds with many independent bodies.
	if (m_translationBodies.size() > 32) return false;
	for (int attempt = 0; attempt < 2; ++attempt)
	{
		m_contactZ.clear(); m_contactAZ.clear(); m_contactZRows.clear(); m_contactAZRows.clear();
		for (int b = 0; b < m_translationBodies.size(); ++b)
		{
			const TranslationBody& body = m_translationBodies[b];
			for (int d = 0; d < (attempt ? 3 : body.modes); ++d)
			{
				TVStack z, az; z.resize(m_translationWork.size(), btVector3(0, 0, 0)); az.resize(z.size());
				for (int n = body.offset; n < body.offset + body.count; ++n) z[n] = rigidMode(d, n);
				multiply(z, az); m_contactZ.push_back(z); m_contactAZ.push_back(az);
				std::vector<int> zr, ar;
				for (int n = 0; n < z.size(); ++n)
				{
					if (z[n].length2() > 0) zr.push_back(n);
					if (az[n].length2() > 0) ar.push_back(n);
				}
				m_contactZRows.push_back(zr); m_contactAZRows.push_back(ar);
			}
		}
		const int k = int(m_contactZ.size());
		std::vector<btScalar> matrix(k * k, 0);
		for (int r = 0; r < k; ++r) for (int c = 0; c < k; ++c)
			for (int n : m_contactZRows[r]) matrix[r*k+c] += m_contactZ[r][n].dot(m_contactAZ[c][n]);
		m_contactScale.resize(k); m_contactFactor.assign(k*k, 0);
		bool positive = true;
		for (int r = 0; r < k; ++r)
		{
			const btScalar diagonal = matrix[r*k+r];
			if (!(diagonal > 0) || !std::isfinite(double(diagonal))) { positive = false; break; }
			m_contactScale[r] = btSqrt(diagonal);
		}
		for (int r = 0; r < k && positive; ++r) for (int c = 0; c <= r; ++c)
		{
			btScalar value = (matrix[r*k+c] + matrix[c*k+r]) * btScalar(.5) / (m_contactScale[r] * m_contactScale[c]);
			for (int j = 0; j < c; ++j) value -= m_contactFactor[r*k+j] * m_contactFactor[c*k+j];
			if (!std::isfinite(double(value)) || (r == c && value <= btScalar(64)*SIMD_EPSILON)) { positive = false; break; }
			m_contactFactor[r*k+c] = r == c ? btSqrt(value) : value / m_contactFactor[c*k+c];
		}
		if (positive) { m_contactCoarse = m_translationCorrection = true; return true; }
	}
	// Degenerate rotations fall back to translations; nonpositive operators keep D.
	return false;
}

void btDeformableBackwardEulerObjective::solveContactCoarse(std::vector<btScalar>& values) const
{
	const int k = int(values.size());
	for (int r = 0; r < k; ++r)
	{
		values[r] /= m_contactScale[r];
		for (int c = 0; c < r; ++c) values[r] -= m_contactFactor[r*k+c] * values[c];
		values[r] /= m_contactFactor[r*k+r];
	}
	for (int r = k - 1; r >= 0; --r)
	{
		for (int c = r + 1; c < k; ++c) values[r] -= m_contactFactor[c*k+r] * values[c];
		values[r] /= m_contactFactor[r*k+r];
	}
	for (int r = 0; r < k; ++r) values[r] /= m_contactScale[r];
}

bool btDeformableBackwardEulerObjective::setupTranslationCorrection()
{
	m_translationCorrection = false;
	m_contactCoarse = false;
	m_translationBodies.clear();
	bool contact = false;
	for (int f = 0; f < m_lf.size(); ++f)
		if (m_lf[f]->getForceType() == BT_CONTACT_FORCE) contact = true;
	if (contact && !m_contactCoarseEnabled) return false;
	if (!m_implicit) return false;
	// Non-contact operators act independently on each body; only that path
	// can share A*Z columns across bodies.
	for (int f = 0; f < m_lf.size(); ++f)
		if (m_lf[f]->getForceType() != BT_LINEAR_ELASTICITY_FORCE &&
			m_lf[f]->getForceType() != BT_GRAVITY_FORCE &&
			m_lf[f]->getForceType() != BT_NODAL_FORCE &&
			m_lf[f]->getForceType() != BT_CONTACT_FORCE &&
			m_lf[f]->getForceType() != BT_VOLUME_BARRIER_FORCE) return false;
	btAlignedObjectArray<int> constrained;
	constrained.resize(m_nodes.size(), 0);
	for (int c = 0; c < m_projection.m_lagrangeMultipliers.size(); ++c)
	{
		const LagrangeMultiplier& lm = m_projection.m_lagrangeMultipliers[c];
		for (int n = 0; n < lm.m_num_nodes; ++n) constrained[lm.m_indices[n]] = 1;
	}
	if (m_rotationCorrectionEnabled)
		for (int d = 0; d < 3; ++d) m_rotationZ[d].resize(m_nodes.size());
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
			if (m_rotationCorrectionEnabled && entry.count >= 3)
			{
				const btVector3 origin = body.m_nodes[0].m_x;
				btVector3 centerOffset(0, 0, 0);
				btScalar totalMass = 0;
				for (int n = 0; n < entry.count; ++n)
				{
					const btScalar mass = 1 / body.m_nodes[n].m_im;
					totalMass += mass; centerOffset += (body.m_nodes[n].m_x - origin) * mass;
				}
				centerOffset /= totalMass;
				btScalar radius = 0;
				for (int n = 0; n < entry.count; ++n)
					radius = btMax(radius, (body.m_nodes[n].m_x - origin - centerOffset).length());
				if (radius > 0 && std::isfinite(double(radius)))
				{
					entry.modes = 6;
					for (int n = 0; n < entry.count; ++n)
					{
						const btVector3 relative = (body.m_nodes[n].m_x - origin - centerOffset) / radius;
						for (int d = 0; d < 3; ++d)
						{
							btVector3 axis(0, 0, 0); axis[d] = 1;
							m_rotationZ[d][offset + n] = axis.cross(relative);
						}
					}
				}
			}
			m_translationBodies.push_back(entry);
		}
		offset += body.m_nodes.size();
	}
	if (m_translationBodies.size() == 0) return false;
	m_translationWork.resize(m_nodes.size() + m_projection.m_lagrangeMultipliers.size());
	if (contact) return setupContactCoarse();
	int modes = 3;
	for (int b = 0; b < m_translationBodies.size(); ++b) modes = btMax(modes, m_translationBodies[b].modes);
	for (int d = 0; d < modes; ++d)
	{
		for (int n = 0; n < m_translationWork.size(); ++n) m_translationWork[n].setZero();
		for (int b = 0; b < m_translationBodies.size(); ++b)
		{
			const TranslationBody& entry = m_translationBodies[b];
			if (d < entry.modes)
				for (int n = entry.offset; n < entry.offset + entry.count; ++n) m_translationWork[n] = rigidMode(d, n);
		}
		m_translationAZ[d].resize(m_translationWork.size());
		multiply(m_translationWork, m_translationAZ[d]);
	}
	for (int b = m_translationBodies.size() - 1; b >= 0; --b)
	{
		TranslationBody& entry = m_translationBodies[b];
		if (entry.modes == 6)
		{
			btScalar coarse[6][6] = {};
			for (int n = entry.offset; n < entry.offset + entry.count; ++n)
				for (int r = 0; r < 6; ++r) for (int d = 0; d < 6; ++d)
					coarse[r][d] += rigidMode(r, n).dot(m_translationAZ[d][n]);
			if (invertRigidCoarse(coarse, entry.inverseRigid)) continue;
			entry.modes = 3;
		}
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
	if (m_translationCorrection && m_contactCoarse)
	{
		const int k = int(m_contactZ.size());
		std::vector<btScalar> coarse(k, 0), coupling(k, 0);
		m_translationWork = x;
		for (int d = 0; d < k; ++d)
			for (int n : m_contactZRows[d]) coarse[d] += m_contactZ[d][n].dot(x[n]);
		solveContactCoarse(coarse);
		for (int d = 0; d < k; ++d)
			for (int n : m_contactAZRows[d]) m_translationWork[n] -= m_contactAZ[d][n] * coarse[d];
		m_preconditioner->operator()(m_translationWork, b);
		for (int d = 0; d < k; ++d)
			for (int n : m_contactAZRows[d]) coupling[d] += m_contactAZ[d][n].dot(b[n]);
		solveContactCoarse(coupling);
		for (int d = 0; d < k; ++d)
			for (int n : m_contactZRows[d]) b[n] += m_contactZ[d][n] * (coarse[d] - coupling[d]);
		return;
	}
	if (!m_translationCorrection)
	{
		m_preconditioner->operator()(x, b);
		return;
	}
	// B = Q + (I-Q*A)*D*(I-A*Q), with rigid coarse modes per free body.
	// Constrained bodies and multiplier entries keep the original D action.
	{
		m_translationWork = x;
	}
	for (int body = 0; body < m_translationBodies.size(); ++body)
	{
		TranslationBody& entry = m_translationBodies[body];
		if (entry.modes == 6)
		{
			btScalar sum[6] = {};
			for (int n = entry.offset; n < entry.offset + entry.count; ++n)
				for (int d = 0; d < 6; ++d) sum[d] += rigidMode(d, n).dot(x[n]);
			applyRigidInverse(entry.inverseRigid, sum, entry.coarseRigid);
			for (int n = entry.offset; n < entry.offset + entry.count; ++n)
				for (int d = 0; d < 6; ++d) m_translationWork[n] -= m_translationAZ[d][n] * entry.coarseRigid[d];
			continue;
		}
		{
			btVector3 sum(0, 0, 0);
			for (int n = entry.offset; n < entry.offset + entry.count; ++n) sum += x[n];
			entry.coarse = entry.inverse * sum;
		}
		{
#if (BT_DEFORMABLE_OPTIMIZATION_MASK & 4)
			if (entry.count > 0)
				btApplyTranslationInput(entry.count, &m_translationWork[entry.offset],
					&m_translationAZ[0][entry.offset], &m_translationAZ[1][entry.offset],
					&m_translationAZ[2][entry.offset], entry.coarse);
#else
			for (int n = entry.offset; n < entry.offset + entry.count; ++n)
				for (int d = 0; d < 3; ++d) m_translationWork[n] -= m_translationAZ[d][n] * entry.coarse[d];
#endif
		}
	}
	{
		m_preconditioner->operator()(m_translationWork, b);
	}
	for (int body = 0; body < m_translationBodies.size(); ++body)
	{
		const TranslationBody& entry = m_translationBodies[body];
		if (entry.modes == 6)
		{
			btScalar coupling[6] = {}, correction[6];
			for (int n = entry.offset; n < entry.offset + entry.count; ++n)
				for (int d = 0; d < 6; ++d) coupling[d] += m_translationAZ[d][n].dot(b[n]);
			applyRigidInverse(entry.inverseRigid, coupling, correction);
			for (int n = entry.offset; n < entry.offset + entry.count; ++n)
				for (int d = 0; d < 6; ++d) b[n] += rigidMode(d, n) * (entry.coarseRigid[d] - correction[d]);
			continue;
		}
		btVector3 translation;
		{
			btVector3 coupling(0, 0, 0);
#if (BT_DEFORMABLE_OPTIMIZATION_MASK & 8)
			if (entry.count > 0)
				coupling = btComputeTranslationCoupling(entry.count, &b[entry.offset],
					&m_translationAZ[0][entry.offset], &m_translationAZ[1][entry.offset],
					&m_translationAZ[2][entry.offset]);
#else
			for (int n = entry.offset; n < entry.offset + entry.count; ++n)
				for (int d = 0; d < 3; ++d) coupling[d] += m_translationAZ[d][n].dot(b[n]);
#endif
			translation = entry.coarse - entry.inverse * coupling;
		}
		{
			for (int n = entry.offset; n < entry.offset + entry.count; ++n) b[n] += translation;
		}
	}
}

btScalar btDeformableBackwardEulerObjective::correctTranslation(TVStack& x, const TVStack& rhs)
{
	if (!m_translationCorrection) return 0;
	multiply(x, m_translationWork);
	if (m_contactCoarse)
	{
		std::vector<btScalar> coarse(m_contactZ.size(), 0);
		for (int d = 0; d < int(coarse.size()); ++d)
			for (int n = 0; n < x.size(); ++n) coarse[d] += m_contactZ[d][n].dot(rhs[n] - m_translationWork[n]);
		solveContactCoarse(coarse);
		btScalar largest = 0;
		for (int n = 0; n < x.size(); ++n)
		{
			btVector3 delta(0, 0, 0);
			for (int d = 0; d < int(coarse.size()); ++d) delta += m_contactZ[d][n] * coarse[d];
			x[n] += delta; largest = btMax(largest, delta.length());
		}
		return largest;
	}
	btScalar largestCorrection = 0;
	for (int b = 0; b < m_translationBodies.size(); ++b)
	{
		const TranslationBody& entry = m_translationBodies[b];
		if (entry.modes == 6)
		{
			btScalar sum[6] = {}, correction[6];
			for (int n = entry.offset; n < entry.offset + entry.count; ++n)
				for (int d = 0; d < 6; ++d) sum[d] += rigidMode(d, n).dot(rhs[n] - m_translationWork[n]);
			applyRigidInverse(entry.inverseRigid, sum, correction);
			for (int n = entry.offset; n < entry.offset + entry.count; ++n)
			{
				btVector3 delta(0, 0, 0);
				for (int d = 0; d < 6; ++d) delta += rigidMode(d, n) * correction[d];
				x[n] += delta; largestCorrection = btMax(largestCorrection, delta.length());
			}
			continue;
		}
		btVector3 sum(0, 0, 0);
		for (int n = entry.offset; n < entry.offset + entry.count; ++n) sum += rhs[n] - m_translationWork[n];
		const btVector3 correction = entry.inverse * sum;
		for (int n = entry.offset; n < entry.offset + entry.count; ++n) x[n] += correction;
		largestCorrection = btMax(largestCorrection, correction.length());
	}
	return largestCorrection;
}
