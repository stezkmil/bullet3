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

#ifndef BT_BACKWARD_EULER_OBJECTIVE_H
#define BT_BACKWARD_EULER_OBJECTIVE_H
//#include "btConjugateGradient.h"
#include "btDeformableLagrangianForce.h"
#include "btDeformableMassSpringForce.h"
#include "btDeformableGravityForce.h"
#include "btDeformableCorotatedForce.h"
#include "btDeformableMousePickingForce.h"
#include "btDeformableLinearElasticityForce.h"
#include "btDeformableNeoHookeanForce.h"
#include "btDeformableContactProjection.h"
#include "btPreconditioner.h"
// #include "btDeformableMultiBodyDynamicsWorld.h"
#include "LinearMath/btQuickprof.h"
#include <cstdlib>
#include <vector>

class btDeformableBackwardEulerObjective
{
public:
	enum _
	{
		Mass_preconditioner,
		KKT_preconditioner
	};

	typedef btAlignedObjectArray<btVector3> TVStack;
	btScalar m_dt;
	btAlignedObjectArray<btDeformableLagrangianForce*> m_lf;
	btAlignedObjectArray<btSoftBody*>& m_softBodies;
	Preconditioner* m_preconditioner;
	btDeformableContactProjection m_projection;
	const TVStack& m_backupVelocity;
	TVStack m_implicitConstraintDv; // Frozen post-contact velocity change for this timestep.
	btAlignedObjectArray<btSoftBody::Node*> m_nodes;
	bool m_implicit;
	MassPreconditioner* m_massPreconditioner;
	KKTPreconditioner* m_KKTPreconditioner;

	btDeformableBackwardEulerObjective(btAlignedObjectArray<btSoftBody*>& softBodies, const TVStack& backup_v);

	virtual ~btDeformableBackwardEulerObjective();

	void initialize() {}

	// compute the rhs for CG solve, i.e, add the dt scaled implicit force to residual
	void computeResidual(btScalar dt, TVStack& residual);

	// add explicit force to the velocity
	void applyExplicitForce(TVStack& force);

	// apply force to velocity and optionally reset the force to zero
	void applyForce(TVStack& force, bool setZero);

	// compute the norm of the residual
	btScalar computeNorm(const TVStack& residual) const;

	// compute one step of the solve (there is only one solve if the system is linear)
	void computeStep(TVStack& dv, const TVStack& residual, const btScalar& dt);

	// perform A*x = b
	void multiply(const TVStack& x, TVStack& b) const;

	// set initial guess for CG solve
	void initialGuess(TVStack& dv, const TVStack& residual);

	// reset data structure and reset dt
	void reinitialize(bool nodeUpdated, btScalar dt);

	void setDt(btScalar dt);

	// add friction force to residual
	void applyDynamicFriction(TVStack& r);

	// add dv to velocity
	void updateVelocity(const TVStack& dv);

	//set constraints as projections
	void setConstraints(const btContactSolverInfo& infoGlobal);

	// update the projections and project the residual
	void project(TVStack& r)
	{
		BT_PROFILE("project");
		m_projection.project(r);
	}

	// perform precondition M^(-1) x = b
	void precondition(const TVStack& x, TVStack& b);

	// Balanced two-level preconditioning with rigid translation/rotation modes Z.
	// Built from A*Z, so heterogeneous masses and damping need no special case.
	bool m_translationCorrection = false;
	bool m_contactCoarse = false;
	bool m_contactCoarseEnabled = []()
	{
		// TODO: Remove all environment-variable lookups across this feature before merging the feature branch.
		// The final implementation must not read environment variables.
		const char* value = std::getenv("BULLET_DEFORMABLE_CONTACT_COARSE");
		return value && value[0] == '1';
	}();
	std::vector<TVStack> m_contactZ, m_contactAZ;
	std::vector<std::vector<int> > m_contactZRows, m_contactAZRows;
	std::vector<btScalar> m_contactFactor, m_contactScale;
	bool setupContactCoarse();
	void solveContactCoarse(std::vector<btScalar>& values) const;
	TVStack m_translationAZ[6], m_rotationZ[3], m_translationWork;
	bool m_rotationCorrectionEnabled = []()
	{
		const char* value = std::getenv("BULLET_DEFORMABLE_ROTATION_CORRECTION");
		return !value || value[0] != '0';
	}();
	btVector3 rigidMode(int mode, int node) const
	{
		if (mode >= 3) return m_rotationZ[mode - 3][node];
		btVector3 axis(0, 0, 0); axis[mode] = 1; return axis;
	}
	struct TranslationBody
	{
		int offset, count;
		btMatrix3x3 inverse;
		btVector3 coarse;
		int modes = 3;
		btScalar inverseRigid[6][6], coarseRigid[6];
	};
	btAlignedObjectArray<TranslationBody> m_translationBodies;
	bool setupTranslationCorrection();
	btScalar correctTranslation(TVStack& x, const TVStack& rhs);

	// reindex all the vertices
	virtual void updateId()
	{
		size_t node_id = 0;
		size_t face_id = 0;
		m_nodes.clear();
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				psb->m_nodes[j].index = node_id;
				m_nodes.push_back(&psb->m_nodes[j]);
				++node_id;
			}
			for (int j = 0; j < psb->m_faces.size(); ++j)
			{
				psb->m_faces[j].m_index = face_id;
				++face_id;
			}
		}
	}

	const btAlignedObjectArray<btSoftBody::Node*>* getIndices() const
	{
		return &m_nodes;
	}

	void setImplicit(bool implicit)
	{
		m_implicit = implicit;
	}

	// Calculate the total potential energy in the system
	btScalar totalEnergy(btScalar dt);

	void addLagrangeMultiplier(const TVStack& vec, TVStack& extended_vec)
	{
		extended_vec.resize(vec.size() + m_projection.m_lagrangeMultipliers.size());
		for (int i = 0; i < vec.size(); ++i)
		{
			extended_vec[i] = vec[i];
		}
		int offset = vec.size();
		for (int i = 0; i < m_projection.m_lagrangeMultipliers.size(); ++i)
		{
			extended_vec[offset + i].setZero();
		}
	}

	void addLagrangeMultiplierRHS(const TVStack& residual, const TVStack& m_dv, TVStack& extended_residual)
	{
		extended_residual.resize(residual.size() + m_projection.m_lagrangeMultipliers.size());
		for (int i = 0; i < residual.size(); ++i)
		{
			extended_residual[i] = residual[i];
		}
		int offset = residual.size();
		for (int i = 0; i < m_projection.m_lagrangeMultipliers.size(); ++i)
		{
			const LagrangeMultiplier& lm = m_projection.m_lagrangeMultipliers[i];
			extended_residual[offset + i].setZero();
			for (int d = 0; d < lm.m_num_constraints; ++d)
			{
				for (int n = 0; n < lm.m_num_nodes; ++n)
				{
					// Newton adds ddv: C * ddv = C * (contactDv - m_dv).
					// Explicit integration replaces m_dv, so preserve the
					// velocity change already imposed by contacts/anchors.
					const btScalar rhsSign = m_implicit ? btScalar(-1) : btScalar(1);
					extended_residual[offset + i][d] += rhsSign * lm.m_weights[n] * m_dv[lm.m_indices[n]].dot(lm.m_dirs[d]);
					if (m_implicit && m_implicitConstraintDv.size() == m_dv.size())
						extended_residual[offset + i][d] += lm.m_weights[n] * m_implicitConstraintDv[lm.m_indices[n]].dot(lm.m_dirs[d]);
				}
			}
		}
	}

	void calculateContactForce(const TVStack& dv, const TVStack& rhs, TVStack& f)
	{
		size_t counter = 0;
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				const btSoftBody::Node& node = psb->m_nodes[j];
				f[counter] = (node.m_frozen > 0) ? btVector3(0, 0, 0) : dv[counter] / node.m_im;
				++counter;
			}
		}
		for (int i = 0; i < m_lf.size(); ++i)
		{
			// add damping matrix
			m_lf[i]->addScaledDampingForceDifferential(-m_dt, dv, f);
		}
		counter = 0;
		for (; counter < f.size(); ++counter)
		{
			f[counter] = rhs[counter] - f[counter];
		}
	}
};

#endif /* btBackwardEulerObjective_h */
