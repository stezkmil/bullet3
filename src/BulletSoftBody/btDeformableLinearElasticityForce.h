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

#ifndef BT_LINEAR_ELASTICITY_H
#define BT_LINEAR_ELASTICITY_H

#include "btCachedElasticKernel.h"
#include "btDeformableOptimizationConfig.h"
#include "btDeformableLagrangianForce.h"
#include "LinearMath/btQuickprof.h"
#include "LinearMath/btImplicitQRSVD.h"
#include "btSoftBodyInternals.h"
#define TETRA_FLAT_THRESHOLD 0.01
class btDeformableLinearElasticityForce : public btDeformableLagrangianForce
{
public:
	typedef btAlignedObjectArray<btVector3> TVStack;
	btScalar m_mu, m_lambda;
	btScalar m_E, m_nu;  // Young's modulus and Poisson ratio
	btScalar m_damping_alpha, m_damping_beta;
	btDeformableLinearElasticityForce() : m_mu(1), m_lambda(1), m_damping_alpha(0.01), m_damping_beta(0.01)
	{
		updateYoungsModulusAndPoissonRatio();
	}

	btDeformableLinearElasticityForce(btScalar mu, btScalar lambda, btScalar damping_alpha = 0.01, btScalar damping_beta = 0.01) : m_mu(mu), m_lambda(lambda), m_damping_alpha(damping_alpha), m_damping_beta(damping_beta)
	{
		updateYoungsModulusAndPoissonRatio();
	}

	void updateYoungsModulusAndPoissonRatio()
	{
		// conversion from Lame Parameters to Young's modulus and Poisson ratio
		// https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
		m_E = m_mu * (3 * m_lambda + 2 * m_mu) / (m_lambda + m_mu);
		m_nu = m_lambda * 0.5 / (m_mu + m_lambda);
	}

	void updateLameParameters()
	{
		// conversion from Young's modulus and Poisson ratio to Lame Parameters
		// https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
		m_mu = m_E * 0.5 / (1 + m_nu);
		m_lambda = m_E * m_nu / ((1 + m_nu) * (1 - 2 * m_nu));
	}

	void setYoungsModulus(btScalar E)
	{
		m_E = E;
		updateLameParameters();
	}

	void setPoissonRatio(btScalar nu)
	{
		m_nu = nu;
		updateLameParameters();
	}

	void setDamping(btScalar damping_alpha, btScalar damping_beta)
	{
		m_damping_alpha = damping_alpha;
		m_damping_beta = damping_beta;
	}

	void setLameParameters(btScalar mu, btScalar lambda)
	{
		m_mu = mu;
		m_lambda = lambda;
		updateYoungsModulusAndPoissonRatio();
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

	// The damping matrix is calculated using the time n state as described in https://www.math.ucla.edu/~jteran/papers/GSSJT15.pdf to allow line search
	virtual void addScaledDampingForce(btScalar scale, TVStack& force)
	{
		if (m_damping_alpha == 0 && m_damping_beta == 0)
			return;
		btScalar mu_damp = m_damping_beta * m_mu;
		btScalar lambda_damp = m_damping_beta * m_lambda;
		btVector3 grad_N_hat_1st_col = btVector3(-1, -1, -1);
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_tetras.size(); ++j)
			{
				bool close_to_flat = (psb->m_tetraScratches[j].m_J < TETRA_FLAT_THRESHOLD);
				btSoftBody::Tetra& tetra = psb->m_tetras[j];
				btSoftBody::Node* node0 = tetra.m_n[0];
				btSoftBody::Node* node1 = tetra.m_n[1];
				btSoftBody::Node* node2 = tetra.m_n[2];
				btSoftBody::Node* node3 = tetra.m_n[3];
				size_t id0 = node0->index;
				size_t id1 = node1->index;
				size_t id2 = node2->index;
				size_t id3 = node3->index;
				btMatrix3x3 dF = DsFromVelocity(node0, node1, node2, node3) * tetra.m_Dm_inverse;
				if (!close_to_flat)
				{
					dF = psb->m_tetraScratches[j].m_corotation.transpose() * dF;
				}
				btMatrix3x3 I;
				I.setIdentity();
				btMatrix3x3 dP = (dF + dF.transpose()) * mu_damp + I * ((dF[0][0] + dF[1][1] + dF[2][2]) * lambda_damp);
				btMatrix3x3 df_on_node123 = dP * tetra.m_Dm_inverse.transpose();
				if (!close_to_flat)
				{
					df_on_node123 = psb->m_tetraScratches[j].m_corotation * df_on_node123;
				}
				btVector3 df_on_node0 = df_on_node123 * grad_N_hat_1st_col;
				// damping force differential
				btScalar scale1 = scale * tetra.m_element_measure;
				force[id0] -= scale1 * df_on_node0;
				force[id1] -= scale1 * df_on_node123.getColumn(0);
				force[id2] -= scale1 * df_on_node123.getColumn(1);
				force[id3] -= scale1 * df_on_node123.getColumn(2);
			}
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				const btSoftBody::Node& node = psb->m_nodes[j];
				size_t id = node.index;
				if (node.m_frozen <= 0 && node.m_im > 0)
				{
					force[id] -= scale * node.m_v / node.m_im * m_damping_alpha;
				}
				//fprintf(stderr, "force[id] %d %f %f %f scale %f node.m_v %f %f %f node.m_im %f m_damping_alpha %f\n", j, force[id].x(), force[id].y(), force[id].z(), scale, node.m_v.x(), node.m_v.y(), node.m_v.z(), node.m_im, m_damping_alpha);
			}
		}
	}

	virtual double totalElasticEnergy(btScalar dt)
	{
		double energy = 0;
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_tetraScratches.size(); ++j)
			{
				btSoftBody::Tetra& tetra = psb->m_tetras[j];
				btSoftBody::TetraScratch& s = psb->m_tetraScratches[j];
				energy += tetra.m_element_measure * elasticEnergyDensity(s);
			}
		}
		return energy;
	}

	// The damping energy is formulated as in https://www.math.ucla.edu/~jteran/papers/GSSJT15.pdf to allow line search
	virtual double totalDampingEnergy(btScalar dt)
	{
		double energy = 0;
		int sz = 0;
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				sz = btMax(sz, psb->m_nodes[j].index);
			}
		}
		TVStack dampingForce;
		dampingForce.resize(sz + 1);
		for (int i = 0; i < dampingForce.size(); ++i)
			dampingForce[i].setZero();
		addScaledDampingForce(0.5, dampingForce);
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				const btSoftBody::Node& node = psb->m_nodes[j];
				energy -= dampingForce[node.index].dot(node.m_v) / dt;
			}
		}
		return energy;
	}

	double elasticEnergyDensity(const btSoftBody::TetraScratch& s)
	{
		double density = 0;
		// Keep the line-search merit function consistent with firstPiola():
		// rigid element rotation must not contribute elastic energy.
		const btMatrix3x3 corotated_F = s.m_corotation.transpose() * s.m_F;
		btMatrix3x3 epsilon = (corotated_F + corotated_F.transpose()) * 0.5 - btMatrix3x3::getIdentity();
		btScalar trace = epsilon[0][0] + epsilon[1][1] + epsilon[2][2];
		density += m_mu * (epsilon[0].length2() + epsilon[1].length2() + epsilon[2].length2());
		density += m_lambda * trace * trace * 0.5;
		return density;
	}

	virtual void addScaledElasticForce(btScalar scale, TVStack& force)
	{
		btVector3 grad_N_hat_1st_col = btVector3(-1, -1, -1);
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			btScalar max_p = psb->m_cfg.m_maxStress;
			btScalar averagePrincipalStress = 0.0;
			for (int j = 0; j < psb->m_tetras.size(); ++j)
			{
				btSoftBody::Tetra& tetra = psb->m_tetras[j];
				btMatrix3x3 P;
				firstPiola(psb->m_tetraScratches[j], P);
				btScalar trPTP = (P[0].length() + P[1].length() + P[2].length());
				averagePrincipalStress += trPTP;
#define USE_SVD 1
#if USE_SVD
				if (max_p > 0)
				{
					// since we want to clamp the principal stress to max_p, we only need to
					// calculate SVD when sigma_0^2 + sigma_1^2 + sigma_2^2 > max_p * max_p
					btScalar trPTP = (P[0].length2() + P[1].length2() + P[2].length2());
					if (trPTP > max_p * max_p)
					{
						btMatrix3x3 U, V;
						btVector3 sigma;
						singularValueDecomposition(P, U, sigma, V);
						sigma[0] = btMin(sigma[0], max_p);
						sigma[1] = btMin(sigma[1], max_p);
						sigma[2] = btMin(sigma[2], max_p);
						sigma[0] = btMax(sigma[0], -max_p);
						sigma[1] = btMax(sigma[1], -max_p);
						sigma[2] = btMax(sigma[2], -max_p);
						btMatrix3x3 Sigma;
						Sigma.setIdentity();
						Sigma[0][0] = sigma[0];
						Sigma[1][1] = sigma[1];
						Sigma[2][2] = sigma[2];
						P = U * Sigma * V.transpose();
					}
				}
#endif
				//                btVector3 force_on_node0 = P * (tetra.m_Dm_inverse.transpose()*grad_N_hat_1st_col);
				btMatrix3x3 force_on_node123 = psb->m_tetraScratches[j].m_corotation * P * tetra.m_Dm_inverse.transpose();
				btVector3 force_on_node0 = force_on_node123 * grad_N_hat_1st_col;

				btSoftBody::Node* node0 = tetra.m_n[0];
				btSoftBody::Node* node1 = tetra.m_n[1];
				btSoftBody::Node* node2 = tetra.m_n[2];
				btSoftBody::Node* node3 = tetra.m_n[3];
				size_t id0 = node0->index;
				size_t id1 = node1->index;
				size_t id2 = node2->index;
				size_t id3 = node3->index;

				// elastic force
				btScalar scale1 = scale * tetra.m_element_measure;
				force[id0] -= scale1 * force_on_node0;
				force[id1] -= scale1 * force_on_node123.getColumn(0);
				force[id2] -= scale1 * force_on_node123.getColumn(1);
				force[id3] -= scale1 * force_on_node123.getColumn(2);
			}
			averagePrincipalStress /= psb->m_tetras.size();
			psb->m_averagePrincipalStress = averagePrincipalStress;
		}
	}

	virtual void buildDampingForceDifferentialDiagonal(btScalar scale, TVStack& diagA) {}

	virtual bool addImplicitForceDifferentialBlocks(btScalar dt, btAlignedObjectArray<btMatrix3x3>& blocks)
	{
		BT_PROFILE("linearElasticityBlocks");
		const btMatrix3x3 identity = btMatrix3x3::getIdentity();
		for (int b = 0; b < m_softBodies.size(); ++b)
		{
			const btSoftBody* body = m_softBodies[b];
			if (!body->isActive() || body->isStaticObject())
				continue;
			for (int t = 0; t < body->m_tetras.size(); ++t)
			{
				const btSoftBody::Tetra& tetra = body->m_tetras[t];
				const btSoftBody::TetraScratch& scratch = body->m_tetraScratches[t];
				btVector3 gradient[4];
				gradient[1] = tetra.m_Dm_inverse[0];
				gradient[2] = tetra.m_Dm_inverse[1];
				gradient[3] = tetra.m_Dm_inverse[2];
				gradient[0] = -(gradient[1] + gradient[2] + gradient[3]);
				for (int n = 0; n < 4; ++n)
				{
					const btSoftBody::Node& node = *tetra.m_n[n];
					if (node.m_frozen > 0 || node.m_im <= 0)
						continue;
					const btVector3 g = gradient[n];
					const btVector3 rotated = scratch.m_corotation * g;
					// Frozen-rotation tangent used by addScaledElasticForceDifferential:
					// K_ii = V * (mu*|g|^2*I + (mu+lambda)*(R*g)*(R*g)^T).
					const btMatrix3x3 isotropic = identity * (m_mu * g.length2());
					const btMatrix3x3 elastic = isotropic + OuterProduct(rotated, rotated) * (m_mu + m_lambda);
					// Match the unrotated damping differential for nearly flat elements.
					const btVector3 dampingGradient = scratch.m_J < TETRA_FLAT_THRESHOLD ? g : rotated;
					const btMatrix3x3 damping = isotropic + OuterProduct(dampingGradient, dampingGradient) * (m_mu + m_lambda);
					blocks[node.index] += (elastic * (dt * dt) + damping * (dt * m_damping_beta)) * tetra.m_element_measure;
				}
			}
			for (int n = 0; n < body->m_nodes.size(); ++n)
			{
				const btSoftBody::Node& node = body->m_nodes[n];
				if (node.m_frozen <= 0 && node.m_im > 0)
					blocks[node.index] += identity * (dt * m_damping_alpha / node.m_im);
			}
		}
		return true;
	}


	// The damping matrix is calculated using the time n state as described in https://www.math.ucla.edu/~jteran/papers/GSSJT15.pdf to allow line search
	virtual void addScaledDampingForceDifferential(btScalar scale, const TVStack& dv, TVStack& df)
	{
		if (m_damping_alpha == 0 && m_damping_beta == 0)
			return;
		btScalar mu_damp = m_damping_beta * m_mu;
		btScalar lambda_damp = m_damping_beta * m_lambda;
		int numNodes = getNumNodes();
		btAssert(numNodes <= df.size());
		btVector3 grad_N_hat_1st_col = btVector3(-1, -1, -1);
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_tetras.size(); ++j)
			{
				bool close_to_flat = (psb->m_tetraScratches[j].m_J < TETRA_FLAT_THRESHOLD);
				btSoftBody::Tetra& tetra = psb->m_tetras[j];
				btSoftBody::Node* node0 = tetra.m_n[0];
				btSoftBody::Node* node1 = tetra.m_n[1];
				btSoftBody::Node* node2 = tetra.m_n[2];
				btSoftBody::Node* node3 = tetra.m_n[3];
				size_t id0 = node0->index;
				size_t id1 = node1->index;
				size_t id2 = node2->index;
				size_t id3 = node3->index;
				btMatrix3x3 dF = Ds(id0, id1, id2, id3, dv) * tetra.m_Dm_inverse;
				if (!close_to_flat)
				{
					dF = psb->m_tetraScratches[j].m_corotation.transpose() * dF;
				}
				btMatrix3x3 I;
				I.setIdentity();
				btMatrix3x3 dP = (dF + dF.transpose()) * mu_damp + I * ((dF[0][0] + dF[1][1] + dF[2][2]) * lambda_damp);
				btMatrix3x3 df_on_node123 = dP * tetra.m_Dm_inverse.transpose();
				if (!close_to_flat)
				{
					df_on_node123 = psb->m_tetraScratches[j].m_corotation * df_on_node123;
				}
				btVector3 df_on_node0 = df_on_node123 * grad_N_hat_1st_col;

				// damping force differential
				btScalar scale1 = scale * tetra.m_element_measure;
				df[id0] -= scale1 * df_on_node0;
				df[id1] -= scale1 * df_on_node123.getColumn(0);
				df[id2] -= scale1 * df_on_node123.getColumn(1);
				df[id3] -= scale1 * df_on_node123.getColumn(2);
			}
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				const btSoftBody::Node& node = psb->m_nodes[j];
				size_t id = node.index;
				if (node.m_frozen <= 0 && node.m_im > 0)
				{
					df[id] -= scale * dv[id] / node.m_im * m_damping_alpha;
				}
			}
		}
	}

	virtual void addScaledElasticForceDifferential(btScalar scale, const TVStack& dx, TVStack& df)
	{
		int numNodes = getNumNodes();
		btAssert(numNodes <= df.size());
		btVector3 grad_N_hat_1st_col = btVector3(-1, -1, -1);
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_tetras.size(); ++j)
			{
				btSoftBody::Tetra& tetra = psb->m_tetras[j];
				btSoftBody::Node* node0 = tetra.m_n[0];
				btSoftBody::Node* node1 = tetra.m_n[1];
				btSoftBody::Node* node2 = tetra.m_n[2];
				btSoftBody::Node* node3 = tetra.m_n[3];
				size_t id0 = node0->index;
				size_t id1 = node1->index;
				size_t id2 = node2->index;
				size_t id3 = node3->index;
				btMatrix3x3 dF = psb->m_tetraScratches[j].m_corotation.transpose() * Ds(id0, id1, id2, id3, dx) * tetra.m_Dm_inverse;
				btMatrix3x3 dP;
				firstPiolaDifferential(psb->m_tetraScratches[j], dF, dP);
				//                btVector3 df_on_node0 = dP * (tetra.m_Dm_inverse.transpose()*grad_N_hat_1st_col);
				btMatrix3x3 df_on_node123 = psb->m_tetraScratches[j].m_corotation * dP * tetra.m_Dm_inverse.transpose();
				btVector3 df_on_node0 = df_on_node123 * grad_N_hat_1st_col;

				// elastic force differential
				btScalar scale1 = scale * tetra.m_element_measure;
				df[id0] -= scale1 * df_on_node0;
				df[id1] -= scale1 * df_on_node123.getColumn(0);
				df[id2] -= scale1 * df_on_node123.getColumn(1);
				df[id3] -= scale1 * df_on_node123.getColumn(2);
			}
		}
	}

	struct CachedImplicitTetra
	{
		btMatrix3x3 gradients; // Rows are R * grad(N_1..3), expressed in world space.
		btScalar muWeight, lambdaWeight;
		bool usable;
	};
	btAlignedObjectArray<CachedImplicitTetra> m_implicitTetraCache;
	bool m_implicitCacheReady = false;
	btScalar m_implicitCacheDt = 0;

	virtual void prepareImplicitForceDifferential(btScalar dt)
	{
		m_implicitCacheReady = false;
		int count = 0;
		for (int b = 0; b < m_softBodies.size(); ++b)
			if (m_softBodies[b]->isActive() && !m_softBodies[b]->isStaticObject())
				count += m_softBodies[b]->m_tetras.size();
		m_implicitTetraCache.resize(count);
		int index = 0;
		for (int b = 0; b < m_softBodies.size(); ++b)
		{
			const btSoftBody& body = *m_softBodies[b];
			if (!body.isActive() || body.isStaticObject()) continue;
			for (int t = 0; t < body.m_tetras.size(); ++t, ++index)
			{
				CachedImplicitTetra& entry = m_implicitTetraCache[index];
				const btSoftBody::Tetra& tetra = body.m_tetras[t];
				const btSoftBody::TetraScratch& scratch = body.m_tetraScratches[t];
				// World-frame isotropic stress needs an orthogonal corotation.
				// Preserve the original path for flat or non-orthogonal scratch states.
				const btMatrix3x3 metric = scratch.m_corotation * scratch.m_corotation.transpose();
				entry.usable = scratch.m_J >= TETRA_FLAT_THRESHOLD;
				for (int i = 0; i < 3; ++i)
					for (int j = 0; j < 3; ++j)
						entry.usable = entry.usable && btFabs(metric[i][j] - (i == j ? btScalar(1) : btScalar(0))) <= btScalar(64) * SIMD_EPSILON;
				if (!entry.usable) continue;
				entry.gradients = tetra.m_Dm_inverse * scratch.m_corotation.transpose();
				const btScalar weight = (dt * dt + dt * m_damping_beta) * tetra.m_element_measure;
				entry.muWeight = weight * m_mu;
				entry.lambdaWeight = weight * m_lambda;
			}
		}
		m_implicitCacheDt = dt;
		m_implicitCacheReady = true;
	}

	virtual void finishImplicitForceDifferential()
	{
		m_implicitCacheReady = false;
	}

	// Non-flat stiffness damping is beta times the frozen-rotation elastic tangent.
	virtual void addImplicitForceDifferential(btScalar dt, const TVStack& dx, TVStack& df)
	{
		int cacheIndex = 0;
		const bool useCache = m_implicitCacheReady && dt == m_implicitCacheDt;
		const btScalar dampingScale = -dt * m_damping_beta;
		const btScalar elasticScale = -dt * dt;
		int numNodes = getNumNodes();
		btAssert(numNodes <= df.size());
		btVector3 grad_N_hat_1st_col = btVector3(-1, -1, -1);
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_tetras.size(); ++j, ++cacheIndex)
			{
				btSoftBody::Tetra& tetra = psb->m_tetras[j];
				btSoftBody::Node* node0 = tetra.m_n[0];
				btSoftBody::Node* node1 = tetra.m_n[1];
				btSoftBody::Node* node2 = tetra.m_n[2];
				btSoftBody::Node* node3 = tetra.m_n[3];
				size_t id0 = node0->index;
				size_t id1 = node1->index;
				size_t id2 = node2->index;
				size_t id3 = node3->index;
				if (useCache && cacheIndex < m_implicitTetraCache.size() && m_implicitTetraCache[cacheIndex].usable)
				{
					const CachedImplicitTetra& entry = m_implicitTetraCache[cacheIndex];
					// G = sum_j (dx_j-dx_0) (R*gradN_j)^T.
					// R*dP*R^T = mu*(G+G^T) + lambda*trace(G)*I.
#if (BT_DEFORMABLE_OPTIMIZATION_MASK & 2) && defined(BT_USE_DOUBLE_PRECISION) && defined(__AVX2__)
					btApplyCachedElasticTetra(entry.gradients, entry.muWeight, entry.lambdaWeight,
						dx[id0], dx[id1], dx[id2], dx[id3], df[id0], df[id1], df[id2], df[id3]);
#else
					const btMatrix3x3 G = Ds(id0, id1, id2, id3, dx) * entry.gradients;
					const btMatrix3x3 stress = (G + G.transpose()) * entry.muWeight +
						btMatrix3x3::getIdentity() * (entry.lambdaWeight * (G[0][0] + G[1][1] + G[2][2]));
					const btMatrix3x3 contributions = stress * entry.gradients.transpose();
					df[id0] += contributions * grad_N_hat_1st_col;
					df[id1] += contributions.getColumn(0);
					df[id2] += contributions.getColumn(1);
					df[id3] += contributions.getColumn(2);
#endif
					continue;
				}
				btMatrix3x3 dF = psb->m_tetraScratches[j].m_corotation.transpose() * Ds(id0, id1, id2, id3, dx) * tetra.m_Dm_inverse;
				btMatrix3x3 dP;
				firstPiolaDifferential(psb->m_tetraScratches[j], dF, dP);
				//                btVector3 df_on_node0 = dP * (tetra.m_Dm_inverse.transpose()*grad_N_hat_1st_col);
				btMatrix3x3 df_on_node123 = psb->m_tetraScratches[j].m_corotation * dP * tetra.m_Dm_inverse.transpose();
				btVector3 df_on_node0 = df_on_node123 * grad_N_hat_1st_col;

				// elastic force differential
				const bool flat = psb->m_tetraScratches[j].m_J < TETRA_FLAT_THRESHOLD;
				btScalar scale1 = (elasticScale + (flat ? btScalar(0) : dampingScale)) * tetra.m_element_measure;
				df[id0] -= scale1 * df_on_node0;
				df[id1] -= scale1 * df_on_node123.getColumn(0);
				df[id2] -= scale1 * df_on_node123.getColumn(1);
				df[id3] -= scale1 * df_on_node123.getColumn(2);
				if (flat && m_damping_beta != 0)
				{
					// Nearly flat elements retain the original unrotated damping tangent.
					const btMatrix3x3 dampingF = Ds(id0, id1, id2, id3, dx) * tetra.m_Dm_inverse;
					btMatrix3x3 dampingP;
					firstPiolaDifferential(psb->m_tetraScratches[j], dampingF, dampingP);
					const btMatrix3x3 dampingNodes = dampingP * tetra.m_Dm_inverse.transpose();
					const btScalar dampingWeight = dampingScale * tetra.m_element_measure;
					df[id0] -= dampingWeight * (dampingNodes * grad_N_hat_1st_col);
					df[id1] -= dampingWeight * dampingNodes.getColumn(0);
					df[id2] -= dampingWeight * dampingNodes.getColumn(1);
					df[id3] -= dampingWeight * dampingNodes.getColumn(2);
				}
			}
			if (m_damping_alpha != 0)
				for (int n = 0; n < psb->m_nodes.size(); ++n)
				{
					const btSoftBody::Node& node = psb->m_nodes[n];
					if (node.m_frozen <= 0 && node.m_im > 0)
						df[node.index] -= (-dt) * dx[node.index] / node.m_im * m_damping_alpha;
				}
		}
	}

	void firstPiola(const btSoftBody::TetraScratch& s, btMatrix3x3& P)
	{
		btMatrix3x3 corotated_F = s.m_corotation.transpose() * s.m_F;

		btMatrix3x3 epsilon = (corotated_F + corotated_F.transpose()) * 0.5 - btMatrix3x3::getIdentity();
		btScalar trace = epsilon[0][0] + epsilon[1][1] + epsilon[2][2];
		P = epsilon * btScalar(2) * m_mu + btMatrix3x3::getIdentity() * m_lambda * trace;
	}

	// Let P be the first piola stress.
	// This function calculates the dP = dP/dF * dF
	void firstPiolaDifferential(const btSoftBody::TetraScratch& s, const btMatrix3x3& dF, btMatrix3x3& dP)
	{
		btScalar trace = (dF[0][0] + dF[1][1] + dF[2][2]);
		dP = (dF + dF.transpose()) * m_mu + btMatrix3x3::getIdentity() * m_lambda * trace;
	}

	// Let Q be the damping stress.
	// This function calculates the dP = dQ/dF * dF
	void firstPiolaDampingDifferential(const btSoftBody::TetraScratch& s, const btMatrix3x3& dF, btMatrix3x3& dP)
	{
		btScalar mu_damp = m_damping_beta * m_mu;
		btScalar lambda_damp = m_damping_beta * m_lambda;
		btScalar trace = (dF[0][0] + dF[1][1] + dF[2][2]);
		dP = (dF + dF.transpose()) * mu_damp + btMatrix3x3::getIdentity() * lambda_damp * trace;
	}

	virtual void addScaledHessian(btScalar scale)
	{
		btVector3 grad_N_hat_1st_col = btVector3(-1, -1, -1);
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_tetras.size(); ++j)
			{
				btSoftBody::Tetra& tetra = psb->m_tetras[j];
				btMatrix3x3 P;
				firstPiola(psb->m_tetraScratches[j], P);  // make sure scratch is evaluated at x_n + dt * vn
				btMatrix3x3 force_on_node123 = psb->m_tetraScratches[j].m_corotation * P * tetra.m_Dm_inverse.transpose();
				btVector3 force_on_node0 = force_on_node123 * grad_N_hat_1st_col;
				btSoftBody::Node* node0 = tetra.m_n[0];
				btSoftBody::Node* node1 = tetra.m_n[1];
				btSoftBody::Node* node2 = tetra.m_n[2];
				btSoftBody::Node* node3 = tetra.m_n[3];
				btScalar scale1 = scale * (scale + m_damping_beta) * tetra.m_element_measure;  // stiff and stiffness-damping terms;
				node0->m_effectiveMass += OuterProduct(force_on_node0, force_on_node0) * scale1;
				node1->m_effectiveMass += OuterProduct(force_on_node123.getColumn(0), force_on_node123.getColumn(0)) * scale1;
				node2->m_effectiveMass += OuterProduct(force_on_node123.getColumn(1), force_on_node123.getColumn(1)) * scale1;
				node3->m_effectiveMass += OuterProduct(force_on_node123.getColumn(2), force_on_node123.getColumn(2)) * scale1;
			}
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				btSoftBody::Node& node = psb->m_nodes[j];
				if (node.m_frozen <= 0 && node.m_im > 0)
				{
					btMatrix3x3 I;
					I.setIdentity();
					node.m_effectiveMass += I * (scale * (1.0 / node.m_im) * m_damping_alpha);
				}
			}
		}
	}

	virtual btScalar getYoungsModulus() const
	{
		return m_E;
	}

	virtual btScalar getPoissonRatio() const
	{
		return m_nu;
	}

	virtual btDeformableLagrangianForceType getForceType()
	{
		return BT_LINEAR_ELASTICITY_FORCE;
	}
};
#endif /* BT_LINEAR_ELASTICITY_H */
