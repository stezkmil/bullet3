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

/* ====== Overview of the Deformable Algorithm ====== */

/*
A single step of the deformable body simulation contains the following main components:
Call internalStepSimulation multiple times, to achieve 240Hz (4 steps of 60Hz).
1. Deformable maintaintenance of rest lengths and volume preservation. Forces only depend on position: Update velocity to a temporary state v_{n+1}^* = v_n + explicit_force * dt / mass, where explicit forces include gravity and elastic forces.
2. Detect discrete collisions between rigid and deformable bodies at position x_{n+1}^* = x_n + dt * v_{n+1}^*.

3a. Solve all constraints, including LCP. Contact, position correction due to numerical drift, friction, and anchors for deformable.

3b. 5 Newton steps (multiple step). Conjugent Gradient solves linear system. Deformable Damping: Then velocities of deformable bodies v_{n+1} are solved in
        M(v_{n+1} - v_{n+1}^*) = damping_force * dt / mass,
   by a conjugate gradient solver, where the damping force is implicit and depends on v_{n+1}.
   Make sure contact constraints are not violated in step b by performing velocity projections as in the paper by Baraff and Witkin https://www.cs.cmu.edu/~baraff/papers/sig98.pdf. Dynamic frictions are treated as a force and added to the rhs of the CG solve, whereas static frictions are treated as constraints similar to contact.
4. Position is updated via x_{n+1} = x_n + dt * v_{n+1}.


The algorithm also closely resembles the one in http://physbam.stanford.edu/~fedkiw/papers/stanford2008-03.pdf
 */

#include "btDeformableVbdSolver.h"
#include "btDeformableMousePickingForce.h"
#include "btDeformableLinearElasticityForce.h"
#include "btDeformableDiagnostics.h"
#include "BulletCollision/Gimpact/btGImpactShape.h"
#include "BulletCollision/CollisionShapes/btCompoundShape.h"
#include <cmath>
#include <memory>
#include <stdio.h>
#include "btDeformableMultiBodyDynamicsWorld.h"
#include "DeformableBodyInplaceSolverIslandCallback.h"
#include "btDeformableBodySolver.h"
#include "LinearMath/btQuickprof.h"
#include "btSoftBodyInternals.h"

btDeformableMultiBodyDynamicsWorld::btDeformableMultiBodyDynamicsWorld(btDispatcher* dispatcher, btBroadphaseInterface* pairCache, btDeformableMultiBodyConstraintSolver* constraintSolver, btCollisionConfiguration* collisionConfiguration, btDeformableBodySolver* deformableBodySolver)
	: btMultiBodyDynamicsWorld(dispatcher, pairCache, (btMultiBodyConstraintSolver*)constraintSolver, collisionConfiguration),
	  m_deformableBodySolver(deformableBodySolver),
	  m_solverCallback(0)
{
	m_drawFlags = fDrawFlags::Std;
	m_drawNodeTree = true;
	m_drawFaceTree = false;
	m_drawClusterTree = false;
	m_sbi.m_broadphase = pairCache;
	m_sbi.m_dispatcher = dispatcher;
	m_sbi.m_sparsesdf.Initialize();
	m_sbi.m_sparsesdf.setDefaultVoxelsz(0.005);
	m_sbi.m_sparsesdf.Reset();

	m_sbi.air_density = (btScalar)1.2;
	m_sbi.water_density = 0;
	m_sbi.water_offset = 0;
	m_sbi.water_normal = btVector3(0, 0, 0);
	m_sbi.m_gravity.setValue(0, -9.8, 0);
	m_internalTime = 0.0;
	m_implicit = false;
	m_lineSearch = false;
	m_useProjection = false;
	m_ccdIterations = 5;
	m_solverDeformableBodyIslandCallback = new DeformableBodyInplaceSolverIslandCallback(constraintSolver, dispatcher);
}


void btDeformableMultiBodyDynamicsWorld::performDiscreteCollisionDetection()
{
	BT_PROFILE("performDiscreteCollisionDetection");

	btDispatcherInfo& dispatchInfo = getDispatchInfo();

	updateAabbs();

	computeOverlappingPairs();

	addSoftsWithSelfCollisionCheckToOverlappingPairs();

	btDispatcher* dispatcher = getDispatcher();
	{
		BT_PROFILE("dispatchAllCollisionPairs");
		if (dispatcher)
			dispatcher->dispatchAllCollisionPairs(m_broadphasePairCache->getOverlappingPairCache(), dispatchInfo, m_dispatcher1);
	}
}

void btDeformableMultiBodyDynamicsWorld::addSoftsWithSelfCollisionCheckToOverlappingPairs()
{
	btBroadphaseInterface* broadphase = getBroadphase();
	btOverlappingPairCache* pairCache = broadphase->getOverlappingPairCache();

	if (!pairCache)
		return;

	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		auto& soft = m_softBodies[i];
		if (soft->getCollisionShape()->getShapeType() != SOFTBODY_SHAPE_PROXYTYPE && (soft->m_cfg.collisions & btSoftBody::fCollision::CL_SELF))  // Do this only for "our" type of softs
		{
			if (!soft->getBroadphaseHandle())
				continue;

			pairCache->addOverlappingPair(soft->getBroadphaseHandle(), soft->getBroadphaseHandle());
		}
	}
}

void btDeformableMultiBodyDynamicsWorld::internalSingleStepSimulation(btScalar timeStep)
{
	if (m_gpuXpbdTgsEnabled && m_softBodies.size()) { gpuXpbdTgsSingleStepSimulation(timeStep); return; }
	if (m_vbdEnabled && m_softBodies.size()) { vbdSingleStepSimulation(timeStep); return; }
	BT_PROFILE("internalSingleStepSimulation");

	if (0 != m_internalPreTickCallback)
	{
		(*m_internalPreTickCallback)(this, timeStep);
	}
	reinitialize(timeStep);

	// add gravity to velocity of rigid and multi bodys
	applyRigidBodyGravity(timeStep);

	///apply gravity and explicit force to velocity, predict motion
	predictUnconstraintMotion(timeStep);

#ifdef BT_SAFE_UPDATE_DEBUG
	fprintf(stderr, "framestart()\n");
#endif
	///perform collision detection that involves rigid/multi bodies
	performDiscreteCollisionDetection();

	if (0 != m_internalPostDiscreteCollisionDetectionTickCallback)
	{
		(*m_internalPostDiscreteCollisionDetectionTickCallback)(this, timeStep);
	}

	btMultiBodyDynamicsWorld::calculateSimulationIslands();

	updateLastSafeTransforms();

#ifdef BT_SAFE_UPDATE_DEBUG
	fprintf(stderr, "frameend()\n");
	fprintf(stderr, "framestart()\n");
	btCollisionObject::gDebug = true;
	performDiscreteCollisionDetection();
	fprintf(stderr, "drawpoint \"VERIFY\" [0,0,0][1,1,1,1]\n");
	btCollisionObject::gDebug = false;

	fprintf(stderr, "frameend()\n");
#endif

	beforeSolverCallbacks(timeStep);

	// ///solve contact constraints and then deformable bodies momemtum equation
	solveConstraints(timeStep);

	afterSolverCallbacks(timeStep);

	performDeformableCollisionDetection();

	applyRepulsionForce(timeStep);

	performGeometricCollisions(timeStep);

	integrateTransforms(timeStep);

	///update vehicle simulation
	btMultiBodyDynamicsWorld::updateActions(timeStep);

	updateActivationState(timeStep);

	if (0 != m_internalTickCallback)
	{
		(*m_internalTickCallback)(this, timeStep);
	}

	// End solver-wise simulation step
	// ///////////////////////////////

	++m_executed_step_counter;
	m_dispatchInfo.m_stepCounter = m_executed_step_counter;
}

void btDeformableMultiBodyDynamicsWorld::performDeformableCollisionDetection()
{
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		m_softBodies[i]->m_softSoftCollision = true;
	}

	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		for (int j = i; j < m_softBodies.size(); ++j)
		{
			if (m_softBodies[i]->getCollisionShape()->getShapeType() != SOFTBODY_SHAPE_PROXYTYPE || m_softBodies[j]->getCollisionShape()->getShapeType() != SOFTBODY_SHAPE_PROXYTYPE)
			{
				// If any shape is not the default soft shape, then the collision is checked elsewhere
				continue;
			}
			m_softBodies[i]->defaultCollisionHandler(m_softBodies[j]);
		}
	}

	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		m_softBodies[i]->m_softSoftCollision = false;
	}
}

void btDeformableMultiBodyDynamicsWorld::updateActivationState(btScalar timeStep)
{
	for (int i = 0; i < m_softBodies.size(); i++)
	{
		btSoftBody* psb = m_softBodies[i];
		psb->updateDeactivation(timeStep);
		if (psb->wantsSleeping())
		{
			if (psb->getActivationState() == ACTIVE_TAG)
				psb->setActivationState(WANTS_DEACTIVATION);
			if (psb->getActivationState() == ISLAND_SLEEPING)
			{
				psb->setZeroVelocity();
			}
		}
		else
		{
			if (psb->getActivationState() != DISABLE_DEACTIVATION)
				psb->setActivationState(ACTIVE_TAG);
		}
	}
	btMultiBodyDynamicsWorld::updateActivationState(timeStep);
}

void btDeformableMultiBodyDynamicsWorld::applyRepulsionForce(btScalar timeStep)
{
	BT_PROFILE("btDeformableMultiBodyDynamicsWorld::applyRepulsionForce");
	for (int i = 0; i < m_softBodies.size(); i++)
	{
		btSoftBody* psb = m_softBodies[i];
		if (psb->isActive() && !psb->isStaticObject())
		{
			psb->applyRepulsionForce(timeStep, true);
		}
	}
}

void btDeformableMultiBodyDynamicsWorld::performGeometricCollisions(btScalar timeStep)
{
	BT_PROFILE("btDeformableMultiBodyDynamicsWorld::performGeometricCollisions");
	// refit the BVH tree for CCD
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		if (psb->isActive() && !psb->isStaticObject())
		{
			m_softBodies[i]->updateFaceTree(true, false);
			m_softBodies[i]->updateNodeTree(true, false);
			for (int j = 0; j < m_softBodies[i]->m_faces.size(); ++j)
			{
				btSoftBody::Face& f = m_softBodies[i]->m_faces[j];
				f.m_n0 = (f.m_n[1]->m_x - f.m_n[0]->m_x).cross(f.m_n[2]->m_x - f.m_n[0]->m_x);
			}
		}
	}

	// clear contact points & update DBVT
	for (int r = 0; r < m_ccdIterations; ++r)
	{
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (psb->isActive() && !psb->isStaticObject())
			{
				// clear contact points in the previous iteration
				psb->m_faceNodeContactsCCD.clear();

				// update m_q and normals for CCD calculation
				for (int j = 0; j < psb->m_nodes.size(); ++j)
				{
					psb->m_nodes[j].m_q = psb->m_nodes[j].m_x + timeStep * psb->m_nodes[j].m_v;
				}
				for (int j = 0; j < psb->m_faces.size(); ++j)
				{
					btSoftBody::Face& f = psb->m_faces[j];
					f.m_n1 = (f.m_n[1]->m_q - f.m_n[0]->m_q).cross(f.m_n[2]->m_q - f.m_n[0]->m_q);
					f.m_vn = (f.m_n[1]->m_v - f.m_n[0]->m_v).cross(f.m_n[2]->m_v - f.m_n[0]->m_v) * timeStep * timeStep;
				}
			}
		}

		// apply CCD to register new contact points
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			for (int j = i; j < m_softBodies.size(); ++j)
			{
				btSoftBody* psb1 = m_softBodies[i];
				btSoftBody* psb2 = m_softBodies[j];
				if (psb1->isActive() && !psb1->isStaticObject() && psb2->isActive() && !psb2->isStaticObject())
				{
					if (m_softBodies[i]->getCollisionShape()->getShapeType() != SOFTBODY_SHAPE_PROXYTYPE || m_softBodies[j]->getCollisionShape()->getShapeType() != SOFTBODY_SHAPE_PROXYTYPE)
						continue;
					m_softBodies[i]->geometricCollisionHandler(m_softBodies[j]);
				}
			}
		}

		int penetration_count = 0;
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (psb->isActive() && !psb->isStaticObject())
			{
				penetration_count += psb->m_faceNodeContactsCCD.size();
				;
			}
		}
		if (penetration_count == 0)
		{
			break;
		}

		// apply inelastic impulse
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (psb->isActive() && !psb->isStaticObject())
			{
				psb->applyRepulsionForce(timeStep, false);
			}
		}
	}
}

void btDeformableMultiBodyDynamicsWorld::softBodySelfCollision()
{
	BT_PROFILE("btDeformableMultiBodyDynamicsWorld::softBodySelfCollision");
	for (int i = 0; i < m_softBodies.size(); i++)
	{
		btSoftBody* psb = m_softBodies[i];
		if (psb->isActive() && !psb->isStaticObject())
		{
			psb->defaultCollisionHandler(psb);
		}
	}
}

void btDeformableMultiBodyDynamicsWorld::positionCorrection(btScalar timeStep)
{
	// correct the position of rigid bodies with temporary velocity generated from split impulse
	btContactSolverInfo infoGlobal;
	btVector3 zero(0, 0, 0);
	for (int i = 0; i < m_nonStaticRigidBodies.size(); ++i)
	{
		btRigidBody* rb = m_nonStaticRigidBodies[i];
		//correct the position/orientation based on push/turn recovery
		btTransform newTransform;
		btVector3 pushVelocity = rb->getPushVelocity();
		btVector3 turnVelocity = rb->getTurnVelocity();
		if (pushVelocity[0] != 0.f || pushVelocity[1] != 0 || pushVelocity[2] != 0 || turnVelocity[0] != 0.f || turnVelocity[1] != 0 || turnVelocity[2] != 0)
		{
			btTransformUtil::integrateTransform(rb->getWorldTransform(), pushVelocity, turnVelocity * infoGlobal.m_splitImpulseTurnErp, timeStep, newTransform);
			rb->setWorldTransform(newTransform);
			rb->setPushVelocity(zero);
			rb->setTurnVelocity(zero);
		}
	}
}

void btDeformableMultiBodyDynamicsWorld::integrateTransforms(btScalar timeStep)
{
	BT_PROFILE("integrateTransforms");
	positionCorrection(timeStep);
	btMultiBodyDynamicsWorld::integrateTransforms(timeStep);
	m_deformableBodySolver->applyTransforms(timeStep);
}

void btDeformableMultiBodyDynamicsWorld::solveConstraints(btScalar timeStep)
{
	BT_PROFILE("btDeformableMultiBodyDynamicsWorld::solveConstraints");
	// save v_{n+1}^* velocity after explicit forces
	m_deformableBodySolver->backupVelocity();

	// set up constraints among multibodies and between multibodies and deformable bodies
	setupConstraints();

	// solve contact constraints
	solveContactConstraints();

	// set up the directions in which the velocity does not change in the momentum solve
	if (m_useProjection)
		m_deformableBodySolver->setProjection();
	else
		m_deformableBodySolver->setLagrangeMultiplier();

	// for explicit scheme, m_backupVelocity = v_{n+1}^*
	// for implicit scheme, m_backupVelocity = v_n
	// Here, set dv = v_{n+1} - v_n for nodes in contact
	m_deformableBodySolver->setupDeformableSolve(m_implicit);

	// At this point, dv should be golden for nodes in contact
	// proceed to solve deformable momentum equation
	m_deformableBodySolver->solveDeformableConstraints(timeStep);
}

void btDeformableMultiBodyDynamicsWorld::setupConstraints()
{
	// set up constraints between multibody and deformable bodies
	m_deformableBodySolver->setConstraints(m_solverInfo);

	// set up constraints among multibodies
	{
		sortConstraints();
		// setup the solver callback
		btMultiBodyConstraint** sortedMultiBodyConstraints = m_sortedMultiBodyConstraints.size() ? &m_sortedMultiBodyConstraints[0] : 0;
		btTypedConstraint** constraintsPtr = getNumConstraints() ? &m_sortedConstraints[0] : 0;
		m_solverDeformableBodyIslandCallback->setup(&m_solverInfo, constraintsPtr, m_sortedConstraints.size(), sortedMultiBodyConstraints, m_sortedMultiBodyConstraints.size(), getDebugDrawer());

		// build islands
		m_islandManager->buildIslands(getCollisionWorld()->getDispatcher(), getCollisionWorld());
	}
}

void btDeformableMultiBodyDynamicsWorld::sortConstraints()
{
	m_sortedConstraints.resize(m_constraints.size());
	int i;
	for (i = 0; i < getNumConstraints(); i++)
	{
		m_sortedConstraints[i] = m_constraints[i];
	}
	m_sortedConstraints.quickSort(btSortConstraintOnIslandPredicate2());

	m_sortedMultiBodyConstraints.resize(m_multiBodyConstraints.size());
	for (i = 0; i < m_multiBodyConstraints.size(); i++)
	{
		m_sortedMultiBodyConstraints[i] = m_multiBodyConstraints[i];
	}
	m_sortedMultiBodyConstraints.quickSort(btSortMultiBodyConstraintOnIslandPredicate());
}

void btDeformableMultiBodyDynamicsWorld::solveContactConstraints()
{
	// process constraints on each island
	m_islandManager->processIslands(getCollisionWorld()->getDispatcher(), getCollisionWorld(), m_solverDeformableBodyIslandCallback);

	// process deferred
	m_solverDeformableBodyIslandCallback->processConstraints();
	m_constraintSolver->allSolved(m_solverInfo, m_debugDrawer);

	// write joint feedback
	{
		for (int i = 0; i < this->m_multiBodies.size(); i++)
		{
			btMultiBody* bod = m_multiBodies[i];

			bool isSleeping = false;

			if (bod->getBaseCollider() && bod->getBaseCollider()->getActivationState() == ISLAND_SLEEPING)
			{
				isSleeping = true;
			}
			for (int b = 0; b < bod->getNumLinks(); b++)
			{
				if (bod->getLink(b).m_collider && bod->getLink(b).m_collider->getActivationState() == ISLAND_SLEEPING)
					isSleeping = true;
			}

			if (!isSleeping)
			{
				//useless? they get resized in stepVelocities once again (AND DIFFERENTLY)
				m_scratch_r.resize(bod->getNumLinks() + 1);  //multidof? ("Y"s use it and it is used to store qdd)
				m_scratch_v.resize(bod->getNumLinks() + 1);
				m_scratch_m.resize(bod->getNumLinks() + 1);

				if (bod->internalNeedsJointFeedback())
				{
					if (!bod->isUsingRK4Integration())
					{
						if (bod->internalNeedsJointFeedback())
						{
							bool isConstraintPass = true;
							bod->computeAccelerationsArticulatedBodyAlgorithmMultiDof(m_solverInfo.m_timeStep, m_scratch_r, m_scratch_v, m_scratch_m, isConstraintPass,
																					  getSolverInfo().m_jointFeedbackInWorldSpace,
																					  getSolverInfo().m_jointFeedbackInJointFrame);
						}
					}
				}
			}
		}
	}

	for (int i = 0; i < this->m_multiBodies.size(); i++)
	{
		btMultiBody* bod = m_multiBodies[i];
		bod->processDeltaVeeMultiDof2();
	}
}

void btDeformableMultiBodyDynamicsWorld::addSoftBody(btSoftBody* body, int collisionFilterGroup, int collisionFilterMask)
{
	m_softBodies.push_back(body);

	// Set the soft body solver that will deal with this body
	// to be the world's solver
	body->setSoftBodySolver(m_deformableBodySolver);

	btCollisionWorld::addCollisionObject(body,
										 collisionFilterGroup,
										 collisionFilterMask);
}

void btDeformableMultiBodyDynamicsWorld::predictUnconstraintMotion(btScalar timeStep)
{
	BT_PROFILE("predictUnconstraintMotion");
	btMultiBodyDynamicsWorld::predictUnconstraintMotion(timeStep);
	m_deformableBodySolver->predictMotion(timeStep);
}

void btDeformableMultiBodyDynamicsWorld::setGravity(const btVector3& gravity)
{
	btDiscreteDynamicsWorld::setGravity(gravity);
	m_deformableBodySolver->setGravity(gravity);
}

void btDeformableMultiBodyDynamicsWorld::setMaxNewtonIterations(int maxNewtonIterations)
{
	m_deformableBodySolver->setMaxNewtonIterations(maxNewtonIterations);
}

void btDeformableMultiBodyDynamicsWorld::setNewtonTolerance(btScalar tolerance)
{
	m_deformableBodySolver->setNewtonTolerance(tolerance);
}

int btDeformableMultiBodyDynamicsWorld::getMaxNewtonIterations() const
{
	return m_deformableBodySolver->getMaxNewtonIterations();
}

void btDeformableMultiBodyDynamicsWorld::reinitialize(btScalar timeStep)
{
	m_internalTime += timeStep;
	m_solverInfo.m_deformable_implicit = m_implicit;
	m_deformableBodySolver->setImplicit(m_implicit);
	m_deformableBodySolver->setLineSearch(m_lineSearch);
	m_deformableBodySolver->reinitialize(m_softBodies, timeStep);
	btDispatcherInfo& dispatchInfo = btMultiBodyDynamicsWorld::getDispatchInfo();
	dispatchInfo.m_timeStep = timeStep;
	dispatchInfo.m_debugDraw = btMultiBodyDynamicsWorld::getDebugDrawer();
	btMultiBodyDynamicsWorld::getSolverInfo().m_timeStep = timeStep;
	if (m_useProjection)
	{
		m_deformableBodySolver->m_useProjection = true;
		m_deformableBodySolver->setStrainLimiting(true);
		m_deformableBodySolver->setPreconditioner(btDeformableBackwardEulerObjective::Mass_preconditioner);
	}
	else
	{
		m_deformableBodySolver->m_useProjection = false;
		m_deformableBodySolver->setStrainLimiting(false);
		m_deformableBodySolver->setPreconditioner(btDeformableBackwardEulerObjective::KKT_preconditioner);
	}
}

void btDeformableMultiBodyDynamicsWorld::debugDrawWorld()
{
	btMultiBodyDynamicsWorld::debugDrawWorld();

	for (int i = 0; i < getSoftBodyArray().size(); i++)
	{
		btSoftBody* psb = (btSoftBody*)getSoftBodyArray()[i];
		{
			btSoftBodyHelpers::DrawFrame(psb, getDebugDrawer());
			btSoftBodyHelpers::Draw(psb, getDebugDrawer(), getDrawFlags());
		}
	}
}

void btDeformableMultiBodyDynamicsWorld::applyRigidBodyGravity(btScalar timeStep)
{
	// Gravity is applied in stepSimulation and then cleared here and then applied here and then cleared here again
	// so that 1) gravity is applied to velocity before constraint solve and 2) gravity is applied in each substep
	// when there are multiple substeps
	btMultiBodyDynamicsWorld::applyGravity();
	// integrate rigid body gravity
	for (int i = 0; i < m_nonStaticRigidBodies.size(); ++i)
	{
		btRigidBody* rb = m_nonStaticRigidBodies[i];
		rb->integrateVelocities(timeStep);
	}

	// integrate multibody gravity
	{
		forwardKinematics();
		clearMultiBodyConstraintForces();
		{
			for (int i = 0; i < this->m_multiBodies.size(); i++)
			{
				btMultiBody* bod = m_multiBodies[i];

				bool isSleeping = false;

				if (bod->getBaseCollider() && bod->getBaseCollider()->getActivationState() == ISLAND_SLEEPING)
				{
					isSleeping = true;
				}
				for (int b = 0; b < bod->getNumLinks(); b++)
				{
					if (bod->getLink(b).m_collider && bod->getLink(b).m_collider->getActivationState() == ISLAND_SLEEPING)
						isSleeping = true;
				}

				if (!isSleeping)
				{
					m_scratch_r.resize(bod->getNumLinks() + 1);
					m_scratch_v.resize(bod->getNumLinks() + 1);
					m_scratch_m.resize(bod->getNumLinks() + 1);
					bool isConstraintPass = false;
					{
						if (!bod->isUsingRK4Integration())
						{
							bod->computeAccelerationsArticulatedBodyAlgorithmMultiDof(m_solverInfo.m_timeStep,
																					  m_scratch_r, m_scratch_v, m_scratch_m, isConstraintPass,
																					  getSolverInfo().m_jointFeedbackInWorldSpace,
																					  getSolverInfo().m_jointFeedbackInJointFrame);
						}
						else
						{
							btAssert(" RK4Integration is not supported");
						}
					}
				}
			}
		}
	}
	clearGravity();
}

void btDeformableMultiBodyDynamicsWorld::clearGravity()
{
	BT_PROFILE("btMultiBody clearGravity");
	// clear rigid body gravity
	for (int i = 0; i < m_nonStaticRigidBodies.size(); i++)
	{
		btRigidBody* body = m_nonStaticRigidBodies[i];
		if (body->isActive() && !body->isStaticObject())
		{
			body->clearGravity();
		}
	}
	// clear multibody gravity
	for (int i = 0; i < this->m_multiBodies.size(); i++)
	{
		btMultiBody* bod = m_multiBodies[i];

		bool isSleeping = false;

		if (bod->getBaseCollider() && bod->getBaseCollider()->getActivationState() == ISLAND_SLEEPING)
		{
			isSleeping = true;
		}
		for (int b = 0; b < bod->getNumLinks(); b++)
		{
			if (bod->getLink(b).m_collider && bod->getLink(b).m_collider->getActivationState() == ISLAND_SLEEPING)
				isSleeping = true;
		}

		if (!isSleeping)
		{
			bod->addBaseForce(-m_gravity * bod->getBaseMass());

			for (int j = 0; j < bod->getNumLinks(); ++j)
			{
				bod->addLinkForce(j, -m_gravity * bod->getLinkMass(j));
			}
		}
	}
}

void btDeformableMultiBodyDynamicsWorld::beforeSolverCallbacks(btScalar timeStep)
{
	if (0 != m_solverCallback)
	{
		(*m_solverCallback)(m_internalTime, this);
	}
}

void btDeformableMultiBodyDynamicsWorld::afterSolverCallbacks(btScalar timeStep)
{
	if (0 != m_solverCallback)
	{
		(*m_solverCallback)(m_internalTime, this);
	}
}

void btDeformableMultiBodyDynamicsWorld::addForce(btSoftBody* psb, btDeformableLagrangianForce* force)
{
	btAlignedObjectArray<btDeformableLagrangianForce*>& forces = *m_deformableBodySolver->getLagrangianForceArray();
	bool added = false;
	for (int i = 0; i < forces.size(); ++i)
	{
		if (forces[i]->getForceType() == force->getForceType())
		{
			forces[i]->addSoftBody(psb);
			added = true;
			break;
		}
	}
	if (!added)
	{
		force->addSoftBody(psb);
		force->setIndices(m_deformableBodySolver->getIndices());
		forces.push_back(force);
	}
}

void btDeformableMultiBodyDynamicsWorld::addForce(btDeformableLagrangianForce* force)
{
	btAlignedObjectArray<btDeformableLagrangianForce*>& forces = *m_deformableBodySolver->getLagrangianForceArray();
	force->setIndices(m_deformableBodySolver->getIndices());
	forces.push_back(force);
}

void btDeformableMultiBodyDynamicsWorld::removeForce(btSoftBody* psb, btDeformableLagrangianForce* force)
{
	btAlignedObjectArray<btDeformableLagrangianForce*>& forces = *m_deformableBodySolver->getLagrangianForceArray();
	int removed_index = -1;
	for (int i = 0; i < forces.size(); ++i)
	{
		if (forces[i]->getForceType() == force->getForceType())
		{
			forces[i]->removeSoftBody(psb);
			if (forces[i]->m_softBodies.size() == 0)
				removed_index = i;
			break;
		}
	}
	if (removed_index >= 0)
		forces.removeAtIndex(removed_index);
}

void btDeformableMultiBodyDynamicsWorld::removeForce(btDeformableLagrangianForce* force)
{
	btAlignedObjectArray<btDeformableLagrangianForce*>& forces = *m_deformableBodySolver->getLagrangianForceArray();
	int removed_index = -1;
	for (int i = 0; i < forces.size(); ++i)
	{
		if (forces[i] == force)
		{
			removed_index = i;
			break;
		}
	}
	if (removed_index >= 0)
		forces.removeAtIndex(removed_index);
}

void btDeformableMultiBodyDynamicsWorld::removeSoftBodyForce(btSoftBody* psb)
{
	btAlignedObjectArray<btDeformableLagrangianForce*>& forces = *m_deformableBodySolver->getLagrangianForceArray();
	for (int i = 0; i < forces.size(); ++i)
	{
		forces[i]->removeSoftBody(psb);
	}
}

void btDeformableMultiBodyDynamicsWorld::removeSoftBody(btSoftBody* body)
{
	clearGpuXpbdTgsKinematicTargets(body);
	if(m_gpuXpbdTgsSettlingBody==body) setGpuXpbdTgsSettling(nullptr,0,0);
	m_vbdAcceptedPositions.erase(body);
	removeSoftBodyForce(body);
	m_softBodies.remove(body);
	btCollisionWorld::removeCollisionObject(body);
	// force a reinitialize so that node indices get updated.
	m_deformableBodySolver->reinitialize(m_softBodies, btScalar(-1));
}

void btDeformableMultiBodyDynamicsWorld::removeCollisionObject(btCollisionObject* collisionObject)
{
	btSoftBody* body = btSoftBody::upcast(collisionObject);
	if (body)
		removeSoftBody(body);
	else
		btDiscreteDynamicsWorld::removeCollisionObject(collisionObject);
}

int btDeformableMultiBodyDynamicsWorld::stepSimulation(btScalar timeStep, int maxSubSteps, btScalar fixedTimeStep)
{
	const bool fixedStep = maxSubSteps > 0;
	// Failure describes this call; a later call may have changed forces or geometry.
	m_vbdStepFailed = false;
	startProfiling(timeStep);

	int numSimulationSubSteps = 0;

	if (maxSubSteps)
	{
		//fixed timestep with interpolation
		m_fixedTimeStep = fixedTimeStep;
		m_localTime += timeStep;
		if (m_localTime >= fixedTimeStep)
		{
			numSimulationSubSteps = int(m_localTime / fixedTimeStep);
			m_localTime -= numSimulationSubSteps * fixedTimeStep;
		}
	}
	else
	{
		//variable timestep
		fixedTimeStep = timeStep;
		m_localTime = m_latencyMotionStateInterpolation ? 0 : timeStep;
		m_fixedTimeStep = 0;
		if (btFuzzyZero(timeStep))
		{
			numSimulationSubSteps = 0;
			maxSubSteps = 0;
		}
		else
		{
			numSimulationSubSteps = 1;
			maxSubSteps = 1;
		}
	}

	//process some debugging flags
	if (getDebugDrawer())
	{
		btIDebugDraw* debugDrawer = getDebugDrawer();
		gDisableDeactivation = (debugDrawer->getDebugMode() & btIDebugDraw::DBG_NoDeactivation) != 0;
	}
	if (numSimulationSubSteps)
	{
		//clamp the number of substeps, to prevent simulation grinding spiralling down to a halt
		int clampedSimulationSteps = (numSimulationSubSteps > maxSubSteps) ? maxSubSteps : numSimulationSubSteps;

		saveKinematicState(fixedTimeStep * clampedSimulationSteps);

		for (int i = 0; i < clampedSimulationSteps; i++)
		{
			internalSingleStepSimulation(fixedTimeStep);
			if (m_vbdStepFailed)
			{
				// Discard this call's transient forces, as on a completed call.
				for (int b = 0; b < m_softBodies.size(); ++b)
					for (int n = 0; n < m_softBodies[b]->m_nodes.size(); ++n)
						m_softBodies[b]->m_nodes[n].m_f.setZero();
				clearForces();
				// Fixed-step wall time is consumed, retaining only its fractional remainder.
				// Physical time and callbacks advance only for completed substeps.
				if (!fixedStep) m_localTime = 0;
				return i;
			}
			synchronizeMotionStates();
		}
	}
	else
	{
		synchronizeMotionStates();
	}

	clearForces();

#ifndef BT_NO_PROFILE
	CProfileManager::Increment_Frame_Counter();
#endif  //BT_NO_PROFILE

	return numSimulationSubSteps;
}

void btDeformableMultiBodyDynamicsWorld::updateLastSafeTransforms()
{
	BT_PROFILE("updateLastSafeTransforms");

	processLastSafeTransforms(m_nonStaticRigidBodies.size() == 0 ? nullptr : reinterpret_cast<btCollisionObject**>(&m_nonStaticRigidBodies[0]), m_nonStaticRigidBodies.size(),
							  m_softBodies.size() == 0 ? nullptr : reinterpret_cast<btCollisionObject**>(&m_softBodies[0]), m_softBodies.size());
}

struct btDeformableVbdMappingCache
{
    btDeformableVbdSolver::SurfaceTopology topology;
	std::vector<unsigned char> key;
	std::vector<btDeformableVbdSolver::SurfaceVertex> vertices;
	std::vector<std::array<int, 3>> faces;
	std::vector<int> owners, faceEnds;
	// The mapping key invalidates these indices whenever mesh topology changes.
	std::vector<std::vector<std::pair<unsigned int, int>>> validationVertices;
	std::vector<btScalar> friction;
	bool valid = false;
};

namespace
{
// Exact input comparison also catches edits made in place without a revision counter.
template <class T> void vbdMappingKeyAppend(std::vector<unsigned char> &key, const T &value)
{
	const auto *bytes = reinterpret_cast<const unsigned char *>(&value);
	key.insert(key.end(), bytes, bytes + sizeof(T));
}
template <class Add> void vbdMappingFields(const std::vector<btSoftBody *> &bodies, btScalar scale, const Add &add)
{
	add(scale);
	add(bodies.size());
	for (auto *body : bodies)
	{
		add(body);
		add(body->m_nodes.size());
		add(body->m_tetras.size());
		add(body->getFriction());
		const auto &pose = body->getWorldTransform();
		for (int i = 0; i < 3; ++i)
		{
			add(pose.getOrigin()[i]);
			for (int j = 0; j < 3; ++j)
				add(pose.getBasis()[i][j]);
		}
		for (int t = 0; t < body->m_tetras.size(); ++t)
			for (int j = 0; j < 4; ++j)
			{
				const auto node = body->m_tetras[t].m_n[j] - &body->m_nodes[0];
				add(node);
			}
		const auto *mapping = body->getCollisionShapeVertexToSimTetra();
		add(mapping);
		if (!mapping)
			continue;
		add(mapping->size());
		for (const auto &entry : *mapping)
		{
			add(entry.vertexToTetra);
			static_assert(sizeof(btVector4) == 4 * sizeof(btScalar), "The mapping key requires four contiguous scalar components");
			add(entry.baryCoordInTetra);
		}
		auto *shape = body->getCollisionShape();
		add(shape);
		add(shape->getShapeType());
		if (shape->getShapeType() != GIMPACT_SHAPE_PROXYTYPE)
			continue;
		auto *gimpact = static_cast<btGImpactShapeInterface *>(shape);
		add(gimpact->getGImpactShapeType());
		if (gimpact->getGImpactShapeType() != CONST_GIMPACT_TRIMESH_SHAPE)
			continue;
		auto *mesh = static_cast<btGImpactMeshShape *>(shape);
		add(mesh->getMeshPartCount());
		for (int part = 0; part < mesh->getMeshPartCount(); ++part)
		{
			auto *piece = mesh->getMeshPart(part);
			for (int j = 0; j < 3; ++j)
				add(piece->getLocalScaling()[j]);
			piece->lockChildShapes();
			const auto *manager = piece->getPrimitiveManager();
			add(manager->get_primitive_count());
			for (int tri = 0; tri < manager->get_primitive_count(); ++tri)
			{
				unsigned int a, b, c;
				manager->get_primitive_indices(tri, a, b, c);
				add(a);
				add(b);
				add(c);
			}
			piece->unlockChildShapes();
		}
	}
}
std::vector<unsigned char> vbdMappingKey(const std::vector<btSoftBody *> &bodies, btScalar scale)
{
	std::vector<unsigned char> key;
	size_t capacity = 1024;
	for (auto *body : bodies)
		if (const auto *mapping = body->getCollisionShapeVertexToSimTetra())
			capacity += mapping->size() * 64;
	key.reserve(capacity);
	vbdMappingFields(bodies, scale, [&](const auto &value) { vbdMappingKeyAppend(key, value); });
	return key;
}
bool vbdMappingKeyMatches(const std::vector<unsigned char> &key, const std::vector<btSoftBody *> &bodies, btScalar scale)
{
	// Compare every serialized field exactly without constructing another copy of an unchanged key.
	size_t offset = 0;
	bool matches = true;
	vbdMappingFields(bodies, scale,
					 [&](const auto &value)
					 {
						 if (matches)
							 matches = offset <= key.size() && sizeof(value) <= key.size() - offset &&
									   std::memcmp(key.data() + offset, &value, sizeof(value)) == 0;
						 offset += sizeof(value);
					 });
	return matches && offset == key.size();
}
struct VbdMappingLease
{
	btDeformableVbdMappingCache &cache;
	btDeformableVbdSolver &solver;
	void swap()
	{
		cache.vertices.swap(solver.surfaceVertices);
		cache.faces.swap(solver.surface);
		cache.owners.swap(solver.surfaceOwners);
		cache.friction.swap(solver.surfaceFrictions);
	}
	~VbdMappingLease()
	{
		swap();
	}
};
} // namespace

namespace
{
struct VbdTriangleCollector : btTriangleCallback
{
	btDeformableVbdSolver &solver;
	btTransform transform;
	btScalar scale, friction;
	VbdTriangleCollector(btDeformableVbdSolver &s, const btTransform &t, btScalar units, btScalar mu)
		: solver(s), transform(t), scale(units), friction(mu)
	{
	}
	void processTriangle(btVector3 *triangle, int, int) override
	{
		btDeformableVbdSolver::Triangle result;
		for (int j = 0; j < 3; ++j)
			result.x[j] = (transform * triangle[j]) * scale;
		if ((result.x[1] - result.x[0]).cross(result.x[2] - result.x[0]).length2() < 1e-30)
			return;
		result.friction = friction;
		solver.barriers.push_back(result);
	}
};
} // namespace

void btDeformableMultiBodyDynamicsWorld::vbdSingleStepSimulation(btScalar timeStep)
{
	btDeformableDiagnostics::StepScope diagnostics(this, m_executed_step_counter, timeStep);
	auto fail = [&](const char *reason)
	{
		m_vbdStepFailed = true;
		btDeformableDiagnostics::write("STEP_FAILED", "solver=vbd reason=%s", reason);
	};
	if (!(m_vbdLengthScale > 0) || !std::isfinite(double(m_vbdLengthScale)) || !(timeStep > 0) || !std::isfinite(double(timeStep)) ||
		m_vbdIterations < 1)
	{
		fail("invalid_settings");
		return;
	}
	if (m_multiBodies.size() || m_solverCallback)
	{
		fail("unsupported_articulated_body_or_callback");
		return;
	}
	std::vector<btSoftBody *> bodies;
	std::vector<int> offsets;
	for (int b = 0; b < m_softBodies.size(); ++b)
		if (!m_softBodies[b]->isStaticObject())
		{
			auto *body = m_softBodies[b];
			if (!body->m_tetras.size())
			{
				fail("requires_tetra_bodies");
				return;
			}
			bodies.push_back(body);
		}
	if (bodies.empty())
	{
		fail("requires_tetra_body");
		return;
	}
	if (m_internalPreTickCallback)
		(*m_internalPreTickCallback)(this, timeStep);
	m_deformableBodySolver->reinitialize(m_softBodies, timeStep);
	const btScalar scale = m_vbdLengthScale;
	if (!m_vbdCollisionCache)
		m_vbdCollisionCache.reset(new btDeformableVbdCollisionCache());
	btDeformableVbdSolver vbd(m_vbdCollisionCache.get());
	if (!m_vbdMappingCache)
		m_vbdMappingCache.reset(new btDeformableVbdMappingCache());
	auto &mappingCache = *m_vbdMappingCache;
	const bool reuseMapping = mappingCache.valid && vbdMappingKeyMatches(mappingCache.key, bodies, scale);
	auto mappingKey = reuseMapping ? std::vector<unsigned char>() : vbdMappingKey(bodies, scale);
	if (!reuseMapping)
	{
		mappingCache.valid = false;
		mappingCache.vertices.clear();
		mappingCache.faces.clear();
		mappingCache.owners.clear();
		mappingCache.friction.clear();
		mappingCache.faceEnds.clear();
		mappingCache.validationVertices.clear();
	}
	VbdMappingLease mappingLease{mappingCache, vbd};
	mappingLease.swap();
	int cachedFace = 0, cachedPart = 0;

	vbd.settings = m_vbdSettings;
	vbd.settings.gap = btMax(vbd.settings.gap, vbd.settings.radius);
	vbd.settings.iterations = m_vbdIterations;
	vbd.settings.collisionUnitsPerMeter = 1 / scale;
	if (m_contactDiscoveryPadding >= 0)
		vbd.settings.gap = vbd.settings.radius + m_contactDiscoveryPadding * scale;
	auto &objective = *m_deformableBodySolver->m_objective;
	vbd.mappedSurface = true;
	for (int owner = 0; owner < int(bodies.size()); ++owner)
	{
		auto *body = bodies[owner];
		const int offset = int(vbd.x.size()), firstFace = reuseMapping ? cachedFace : int(vbd.surface.size()),
				  firstTet = int(vbd.tets.size());
		offsets.push_back(offset);
		btDeformableLinearElasticityForce *material = nullptr;
		TVStack external;
		external.resize(objective.m_nodes.size(), btVector3(0, 0, 0));
		for (int f = 0; f < objective.m_lf.size(); ++f)
		{
			auto *force = objective.m_lf[f];
			bool attached = false;
			for (int b = 0; b < force->m_softBodies.size(); ++b)
				attached = attached || force->m_softBodies[b] == body;
			if (!attached)
				continue;
			switch (force->getForceType())
			{
			case BT_LINEAR_ELASTICITY_FORCE:
				if (material)
				{
					fail("multiple_materials");
					return;
				}
				material = static_cast<btDeformableLinearElasticityForce *>(force);
				break;
			case BT_MOUSE_PICKING_FORCE:
			{
				std::vector<btDeformableMousePickingForce::NodeSpring> springs;
				static_cast<btDeformableMousePickingForce *>(force)->appendNodeSprings(springs);
				for (const auto &spring : springs)
				{
					int node = -1;
					for (int i = 0; i < body->m_nodes.size(); ++i)
						if (&body->m_nodes[i] == spring.node)
						{
							node = i;
							break;
						}
					if (node < 0)
					{
						fail("grab_node_not_in_body");
						return;
					}
					vbd.springs.push_back(
						{offset + node, spring.target * scale, spring.stiffness, spring.damping, spring.maxForce * scale});
				}
				break;
			}
			case BT_NODAL_FORCE:
				force->addScaledForces(1, external);
				break;
			case BT_GRAVITY_FORCE:
				break;
			default:
				fail("unsupported_force_type");
				return;
			}
		}
		if (!material)
		{
			fail("missing_tetra_material");
			return;
		}
		const btScalar mu = material->m_mu / scale, lambda = material->m_lambda / scale;
		for (int n = 0; n < body->m_nodes.size(); ++n)
		{
			const auto &node = body->m_nodes[n];
			int anchorFreeze = 0;
			for (int j = 0; j < body->m_deformableAnchors.size(); ++j)
				if (body->m_deformableAnchors[j].m_node == &node)
					anchorFreeze += body->m_deformableAnchors[j].m_freezeContribution;
			const btScalar mass = node.m_im > 0 && node.m_frozen <= anchorFreeze ? 1 / node.m_im : 0;
			vbd.x.push_back(node.m_x * scale);
			vbd.velocity.push_back(mass > 0 ? node.m_v * scale : btVector3(0, 0, 0));
			vbd.mass.push_back(mass);
			vbd.massDamping.push_back(material->m_damping_alpha);
			vbd.external.push_back((external[node.index] + getGravity() * (body->m_gravityFactor * mass)) * scale);
		}
		for (int k = 0; k < body->m_tetras.size(); ++k)
		{
			const auto &source = body->m_tetras[k];
			btDeformableVbdSolver::Tet tet;
			for (int j = 0; j < 4; ++j)
				tet.nodes[j] = offset + int(source.m_n[j] - &body->m_nodes[0]);
			tet.inverseRest = source.m_Dm_inverse * (1 / scale);
			tet.volume = source.m_element_measure * scale * scale * scale;
			tet.mu = mu;
			tet.lambda = lambda;
			tet.damping = material->m_damping_beta * mu / 4;
			vbd.tets.push_back(tet);
		}
		if (reuseMapping)
		{
			if (body->getCollisionShapeVertexToSimTetra())
			{
				auto *mesh = static_cast<btGImpactMeshShape *>(body->getCollisionShape());
				for (int part = 0; part < mesh->getMeshPartCount(); ++part)
				{
					auto *piece = mesh->getMeshPart(part);
					piece->lockChildShapes();
					const auto &vertices = mappingCache.validationVertices[cachedPart++];
					const int batch = 256, blocks = (int(vertices.size()) + batch - 1) / batch;
					std::vector<unsigned char> validBlocks(blocks, 1);
					// Read-only reconstruction is independent while the mesh stays locked.
					const int workers = vbd.settings.workers >= 1 && vbd.settings.workers <= 256 ? vbd.settings.workers : 1;
					btVbdParallelFor(blocks, workers,
									 [&](int block)
									 {
										 const int end = btMin(int(vertices.size()), (block + 1) * batch);
										 for (int i = block * batch; i < end; ++i)
										 {
											 const auto &vertex = vertices[i];
											 btVector3 original;
											 piece->getVertex(vertex.first, original);
											 const btVector3 expected = (body->getWorldTransform() * original) * scale;
											 const btVector3 actual = vbd.surfaceVertices[vertex.second].position(vbd.x);
											 if (!std::isfinite(double(expected.length2())) || (expected - actual).length2() > 1e-14)
											 {
												 validBlocks[block] = 0;
												 break;
											 }
										 }
									 });
					const bool valid = std::find(validBlocks.begin(), validBlocks.end(), 0) == validBlocks.end();
					piece->unlockChildShapes();
					if (!valid)
					{
						mappingCache.valid = false;
						fail("invalid_render_surface_mapping");
						return;
					}
				}
			}
			cachedFace = mappingCache.faceEnds[owner];
		}
		else
		{
			if (const auto *mapping = body->getCollisionShapeVertexToSimTetra())
			{
				if (body->getCollisionShape()->getShapeType() != GIMPACT_SHAPE_PROXYTYPE)
				{
					fail("mapped_surface_requires_gimpact");
					return;
				}
				auto *shape = static_cast<btGImpactShapeInterface *>(body->getCollisionShape());
				if (shape->getGImpactShapeType() != CONST_GIMPACT_TRIMESH_SHAPE)
				{
					fail("mapped_surface_requires_triangle_mesh");
					return;
				}
				auto *mesh = static_cast<btGImpactMeshShape *>(shape);
				vbd.mappedSurface = true;
				for (int part = 0; part < mesh->getMeshPartCount(); ++part)
				{
					auto *piece = mesh->getMeshPart(part);
					const btVector3 scaling = piece->getLocalScaling();
					const btMatrix3x3 transform =
						body->getWorldTransform().getBasis() * btMatrix3x3(scaling.x(), 0, 0, 0, scaling.y(), 0, 0, 0, scaling.z());
					piece->lockChildShapes();
					const auto *manager = piece->getPrimitiveManager();
					std::map<unsigned int, int> vertices;
					mappingCache.validationVertices.emplace_back();
					bool valid = true;
					for (int tri = 0; tri < manager->get_primitive_count() && valid; ++tri)
					{
						unsigned int indices[3];
						manager->get_primitive_indices(tri, indices[0], indices[1], indices[2]);
						std::array<int, 3> face;
						for (int corner = 0; corner < 3; ++corner)
						{
							const unsigned int index = indices[corner];
							if (index >= mapping->size())
							{
								valid = false;
								break;
							}
							auto found = vertices.find(index);
							if (found == vertices.end())
							{
								const auto &entry = (*mapping)[index];
								if (entry.vertexToTetra >= unsigned(body->m_tetras.size()))
								{
									valid = false;
									break;
								}
								btDeformableVbdSolver::SurfaceVertex vertex;
								vertex.offset = body->getWorldTransform().getOrigin() * scale;
								btScalar sum = 0;
								for (int j = 0; j < 4; ++j)
								{
									const btScalar w = entry.baryCoordInTetra[j];
									if (!std::isfinite(double(w)))
									{
										valid = false;
										break;
									}
									sum += w;
									if (w != 0)
										vertex.support.push_back(
											{offset + int(body->m_tetras[entry.vertexToTetra].m_n[j] - &body->m_nodes[0]), transform * w});
								}
								if (!valid || btFabs(sum - 1) > 1e-4)
								{
									valid = false;
									break;
								}
								const int id = int(vbd.surfaceVertices.size());
								vbd.surfaceVertices.push_back(vertex);
								found = vertices.emplace(index, id).first;
								mappingCache.validationVertices.back().emplace_back(index, id);
							}
							face[corner] = found->second;
						}
						if (valid)
						{
							btPrimitiveTriangle original;
							manager->get_primitive_triangle(tri, original, false);
							for (int corner = 0; corner < 3; ++corner)
							{
								const btVector3 expected = (body->getWorldTransform() * original.m_vertices[corner]) * scale;
								const btVector3 actual = vbd.surfaceVertices[face[corner]].position(vbd.x);
								// Catch frame/scaling mismatches instead of silently colliding a different surface.
								if (!std::isfinite(double(expected.length2())) || (expected - actual).length2() > 1e-14)
									valid = false;
							}
							if (valid)
								vbd.surface.push_back(face);
						}
					}
					piece->unlockChildShapes();
					if (!valid)
					{
						fail("invalid_render_surface_mapping");
						return;
					}
				}
			}
			else
			{
				const int firstVertex = int(vbd.surfaceVertices.size());
				for (int n = 0; n < body->m_nodes.size(); ++n)
				{
					btDeformableVbdSolver::SurfaceVertex vertex;
					vertex.support.push_back({offset + n, btMatrix3x3::getIdentity()});
					vbd.surfaceVertices.push_back(vertex);
				}
				std::map<std::array<int, 3>, int> faces;
				for (int k = firstTet; k < int(vbd.tets.size()); ++k)
					for (int omit = 0; omit < 4; ++omit)
					{
						std::array<int, 3> face;
						int j = 0;
						for (int n = 0; n < 4; ++n)
							if (n != omit)
								face[j++] = firstVertex + vbd.tets[k].nodes[n] - offset;
						std::sort(face.begin(), face.end());
						++faces[face];
					}
				for (const auto &face : faces)
					if (face.second == 1)
						vbd.surface.push_back(face.first);
			}
			mappingCache.faceEnds.push_back(int(vbd.surface.size()));
		}
		for (int face = firstFace; !reuseMapping && face < int(vbd.surface.size()); ++face)
		{
			vbd.surfaceOwners.push_back(owner);
			vbd.surfaceFrictions.push_back(body->getFriction());
		}
		vbd.selfContact.push_back(body->useSelfCollision() || (body->m_cfg.collisions & btSoftBody::fCollision::CL_SELF) != 0);
	}
	if (!reuseMapping)
		mappingCache.key = std::move(mappingKey);
	mappingCache.valid = true;
	std::vector<btCollisionObject *> objects(bodies.begin(), bodies.end());
	for (int b = 0; b < m_collisionObjects.size(); ++b)
		if (std::find(objects.begin(), objects.end(), m_collisionObjects[b]) == objects.end())
			objects.push_back(m_collisionObjects[b]);
	vbd.collisionAllowed.assign(objects.size(), std::vector<bool>(objects.size(), false));
	for (int i = 0; i < int(objects.size()); ++i)
		for (int j = 0; j < int(objects.size()); ++j)
			vbd.collisionAllowed[i][j] =
				i == j || (getPairCache()->needsBroadphaseCollision(objects[i]->getBroadphaseHandle(), objects[j]->getBroadphaseHandle()) &&
						   objects[i]->checkCollideWith(objects[j]) && objects[j]->checkCollideWith(objects[i]));
	std::vector<int> rigidIndices(objects.size(), -1);
	std::vector<btRigidBody *> rigidBodies;
	for (int object = 0; object < int(objects.size()); ++object)
	{
		auto *body = btRigidBody::upcast(objects[object]);
		if (!body)
			continue;
		btDeformableVbdSolver::Rigid rigid;
		rigid.pose = body->getWorldTransform();
		rigid.pose.setOrigin(rigid.pose.getOrigin() * scale);
		rigid.mass = body->isStaticOrKinematicObject() ? btScalar(0) : btScalar(1) / body->getInvMass();
		rigid.velocity = body->getLinearVelocity() * scale;
		rigid.angularVelocity = body->getAngularVelocity();
		rigid.linearFactor = body->getLinearFactor();
		rigid.angularFactor = body->getAngularFactor();
		rigid.linearDamping = body->getLinearDamping();
		rigid.angularDamping = body->getAngularDamping();
		rigid.force = body->getTotalForce() * scale + body->getGravity() * (rigid.mass * scale);
		rigid.torque = body->getTotalTorque() * (scale * scale);
		if (rigid.mass > 0)
		{
			if (body->getFlags() & BT_ENABLE_GYROSCOPIC_FORCE_EXPLICIT)
				rigid.torque -= body->computeGyroscopicForceExplicit(m_solverInfo.m_maxGyroscopicForce) * (scale * scale);
			if (body->getFlags() & BT_ENABLE_GYROSCOPIC_FORCE_IMPLICIT_WORLD)
				rigid.angularVelocity += body->computeGyroscopicImpulseImplicit_World(timeStep);
			if (body->getFlags() & BT_ENABLE_GYROSCOPIC_FORCE_IMPLICIT_BODY)
				rigid.angularVelocity += body->computeGyroscopicImpulseImplicit_Body(timeStep);
		}
		const auto inverse = body->getInvInertiaDiagLocal();
		for (int d = 0; d < 3; ++d)
			rigid.inertia[d] = inverse[d] > 0 ? scale * scale / inverse[d] : btScalar(1);
		btVector3 lo, hi;
		body->getCollisionShape()->getAabb(body->getWorldTransform(), lo, hi);
		rigid.radius =
			btMax((lo - body->getWorldTransform().getOrigin()).length(), (hi - body->getWorldTransform().getOrigin()).length()) * scale;
		rigidIndices[object] = int(vbd.rigids.size());
		vbd.rigids.push_back(rigid);
		rigidBodies.push_back(body);
	}
	for (int b = 0; b < int(bodies.size()); ++b)
	{
		auto *body = bodies[b];
		auto addAnchor = [&](btSoftBody::Node *node, btRigidBody *target, const btVector3 &local, btScalar influence)
		{
			if (influence <= 0)
				return true;
			if (!target)
				return false;
			auto found = std::find(objects.begin(), objects.end(), target);
			if (found == objects.end() || rigidIndices[found - objects.begin()] < 0)
				return false;
			const int index = int(node - &body->m_nodes[0]);
			if (index < 0 || index >= body->m_nodes.size())
				return false;
			vbd.attachments.push_back({offsets[b] + index, rigidIndices[found - objects.begin()], local * scale, btScalar(1e7) * influence,
									   btScalar(10) * influence});
			return true;
		};
		for (int k = 0; k < body->m_anchors.size(); ++k)
		{
			const auto &anchor = body->m_anchors[k];
			if (!addAnchor(anchor.m_node, anchor.m_body, anchor.m_local, anchor.m_influence))
			{
				fail("invalid_anchor");
				return;
			}
		}
		for (int k = 0; k < body->m_deformableAnchors.size(); ++k)
		{
			const auto &anchor = body->m_deformableAnchors[k];
			if (!addAnchor(anchor.m_node, anchor.m_body, anchor.m_local, 1))
			{
				fail("invalid_deformable_anchor");
				return;
			}
		}
	}
	struct VbdJointRange
	{
		btTypedConstraint *constraint;
		int first, count;
		std::vector<btScalar> scales;
	};
	std::vector<VbdJointRange> jointRanges;
	auto jointBody = [&](btRigidBody &body)
	{
		const auto found = std::find(rigidBodies.begin(), rigidBodies.end(), &body);
		return found == rigidBodies.end() ? -1 : int(found - rigidBodies.begin());
	};
	for (int c = 0; c < m_constraints.size(); ++c)
	{
		auto *constraint = m_constraints[c];
		if (!constraint->isEnabled())
			continue;
		const int a = jointBody(constraint->getRigidBodyA()), b = jointBody(constraint->getRigidBodyB());
		if ((a < 0 && !constraint->getRigidBodyA().isStaticObject()) || (b < 0 && !constraint->getRigidBodyB().isStaticObject()))
		{
			fail("joint_body_missing");
			return;
		}
		constraint->buildJacobian();
		btTypedConstraint::btConstraintInfo1 count{};
		constraint->getInfo1(&count);
		if (count.m_numConstraintRows <= 0)
			continue;
		// Some constraints write SIMD-padded btVector3 rows.
		const int size = count.m_numConstraintRows * 4;
		std::vector<btScalar> la(size, 0), aa(size, 0), lb(size, 0), ab(size, 0), rhs(size, 0), cfm(size, m_solverInfo.m_globalCfm),
			lower(size, -SIMD_INFINITY), upper(size, SIMD_INFINITY);
		btTypedConstraint::btConstraintInfo2 info{};
		info.fps = 1 / timeStep;
		info.erp = m_solverInfo.m_erp;
		info.rowskip = 4;
		info.m_J1linearAxis = la.data();
		info.m_J1angularAxis = aa.data();
		info.m_J2linearAxis = lb.data();
		info.m_J2angularAxis = ab.data();
		info.m_constraintError = rhs.data();
		info.cfm = cfm.data();
		info.m_lowerLimit = lower.data();
		info.m_upperLimit = upper.data();
		info.m_numIterations = vbd.settings.iterations;
		info.m_solverMode = m_solverInfo.m_solverMode;
		info.m_damping = m_solverInfo.m_damping;
		constraint->getInfo2(&info);
		VbdJointRange range{constraint, int(vbd.joints.size()), count.m_numConstraintRows, {}};
		for (int row = 0; row < count.m_numConstraintRows; ++row)
		{
			const int k = 4 * row;
			btDeformableVbdSolver::JointRow r;
			r.a = a;
			r.b = b;
			r.linearA = btVector3(la[k], la[k + 1], la[k + 2]);
			r.linearB = btVector3(lb[k], lb[k + 1], lb[k + 2]);
			r.angularA = btVector3(aa[k], aa[k + 1], aa[k + 2]);
			r.angularB = btVector3(ab[k], ab[k + 1], ab[k + 2]);
			const bool linear = r.linearA.length2() + r.linearB.length2() > 0;
			const btScalar impulseScale = linear ? scale : scale * scale;
			if (linear)
			{
				r.angularA *= scale;
				r.angularB *= scale;
			}
			r.rhs = rhs[k] * (linear ? scale : btScalar(1));
			r.cfm = cfm[k] / (linear ? btScalar(1) : scale * scale);
			const btScalar limit = constraint->getBreakingImpulseThreshold();
			r.lower = btMax(lower[k], -limit) * impulseScale;
			r.upper = btMin(upper[k], limit) * impulseScale;
			vbd.joints.push_back(r);
			range.scales.push_back(impulseScale);
		}
		jointRanges.push_back(range);
	}
	if (!vbd.initialize(reuseMapping ? &mappingCache.topology : nullptr))
	{
		mappingCache.valid = false;
		fail(vbd.error);
		return;
	}
	mappingCache.topology = vbd.surfaceTopology();
	btVector3 lo = vbd.x[0] / scale, hi = lo;
	for (const auto &point : vbd.x)
	{
		lo.setMin(point / scale);
		hi.setMax(point / scale);
	}
	vbd.includeSurfaceBounds(scale, lo, hi);
	for (auto *rigid : rigidBodies)
		if (!rigid->isStaticObject())
		{
			btVector3 rlo, rhi;
			rigid->getCollisionShape()->getAabb(rigid->getWorldTransform(), rlo, rhi);
			lo.setMin(rlo);
			hi.setMax(rhi);
		}
	const btVector3 padding(vbd.settings.gap * 2 / scale, vbd.settings.gap * 2 / scale, vbd.settings.gap * 2 / scale);
	lo -= padding;
	hi += padding;
	for (int b = int(bodies.size()); b < int(objects.size()); ++b)
	{
		btCollisionObject *other = objects[b];
		bool collides = false;
		for (int i = 0; i < int(objects.size()); ++i)
			if (!objects[i]->isStaticObject())
				collides = collides || vbd.collisionAllowed[i][b];
		if (!collides)
			continue;
		btVector3 blo, bhi;
		other->getCollisionShape()->getAabb(other->getWorldTransform(), blo, bhi);
		if (!TestAabbAgainstAabb2(lo, hi, blo, bhi))
			continue;
		if (!other->isStaticObject() && rigidIndices[b] < 0)
		{
			fail("non_rigid_moving_barrier");
			return;
		}
		const int collisionRigid = other->isStaticObject() ? -1 : rigidIndices[b];
		std::function<bool(btCollisionShape *, const btTransform &)> collect;
		collect = [&](btCollisionShape *shape, const btTransform &transform)
		{
			if (shape->isCompound())
			{
				auto *compound = static_cast<btCompoundShape *>(shape);
				for (int child = 0; child < compound->getNumChildShapes(); ++child)
					if (!collect(compound->getChildShape(child), transform * compound->getChildTransform(child)))
						return false;
				return true;
			}
			if (shape->isConvex())
			{
				vbd.convexBarriers.push_back({static_cast<btConvexShape *>(shape), transform, btMax(btScalar(0), other->getFriction()), b,
											  collisionRigid, other->getWorldTransform().inverse() * transform});
				return true;
			}
			if (collisionRigid < 0 && shape->getShapeType() == GIMPACT_SHAPE_PROXYTYPE &&
				static_cast<btGImpactShapeInterface *>(shape)->getGImpactShapeType() == CONST_GIMPACT_TRIMESH_SHAPE)
			{
				auto *mesh = static_cast<btGImpactMeshShape *>(shape);
				for (int part = 0; part < mesh->getMeshPartCount(); ++part)
					vbd.nativeBarriers.push_back({mesh->getMeshPart(part), transform, btMax(btScalar(0), other->getFriction()), b, {}});
				return true;
			}
			if (!shape->isConcave())
				return false;
			const btTransform inverse = transform.inverse();
			btVector3 localLo(SIMD_INFINITY, SIMD_INFINITY, SIMD_INFINITY), localHi = -localLo;
			for (int c = 0; c < 8; ++c)
			{
				const btVector3 p = inverse * btVector3(c & 1 ? hi.x() : lo.x(), c & 2 ? hi.y() : lo.y(), c & 4 ? hi.z() : lo.z());
				localLo.setMin(p);
				localHi.setMax(p);
			}
			const int first = int(vbd.barriers.size());
			VbdTriangleCollector collector(vbd, transform, scale, btMax(btScalar(0), other->getFriction()));
			static_cast<btConcaveShape *>(shape)->processAllTriangles(&collector, localLo, localHi);
			for (int tri = first; tri < int(vbd.barriers.size()); ++tri)
			{
				auto &triangle = vbd.barriers[tri];
				triangle.owner = b;
				triangle.rigid = collisionRigid;
				if (collisionRigid >= 0)
					for (int j = 0; j < 3; ++j)
						triangle.local[j] = vbd.rigids[collisionRigid].pose.invXform(triangle.x[j]);
			}
			return true;
		};
		if (!collect(other->getCollisionShape(), other->getWorldTransform()))
		{
			fail("unsupported_barrier_shape");
			return;
		}
	}
	for (auto *body : bodies)
	{
		const auto found = m_vbdAcceptedPositions.find(body);
		if (found == m_vbdAcceptedPositions.end() || found->second.size() != body->m_nodes.size())
		{
			vbd.recoveryPositions.clear();
			break;
		}
		vbd.recoveryPositions.insert(vbd.recoveryPositions.end(), found->second.begin(), found->second.end());
	}
	if (!vbd.step(timeStep))
	{
		fail(vbd.error);
		return;
	}
	for (int b = 0; b < int(bodies.size()); ++b)
	{
		auto *body = bodies[b];
		const int offset = offsets[b];
		// Commit only complete, finite steps; no legacy contact projection follows VBD.
		for (int n = 0; n < body->m_nodes.size(); ++n)
		{
			auto &node = body->m_nodes[n];
			node.m_x = vbd.x[offset + n] / scale;
			node.m_q = node.m_x;
			node.m_v = vbd.velocity[offset + n] / scale;
			node.m_vn = node.m_v;
			node.m_splitv.setZero();
		}
		body->updateDeformation();
		body->updateNormals();
		body->updateBounds();
		updateSingleAabb(body);
	}
	for (int r = 0; r < int(rigidBodies.size()); ++r)
		if (vbd.rigids[r].mass > 0)
		{
			auto *body = rigidBodies[r];
			const auto &rigid = vbd.rigids[r];
			btTransform transform = rigid.pose;
			transform.setOrigin(transform.getOrigin() / scale);
			body->setCenterOfMassTransform(transform);
			body->setInterpolationWorldTransform(transform);
			body->setLinearVelocity(rigid.velocity / scale);
			body->setAngularVelocity(rigid.angularVelocity);
			body->setInterpolationLinearVelocity(body->getLinearVelocity());
			body->setInterpolationAngularVelocity(body->getAngularVelocity());
			updateSingleAabb(body);
		}
	m_internalTime += timeStep;
	m_solverInfo.m_timeStep = timeStep;
	m_dispatchInfo.m_timeStep = timeStep;
	for (int b = 0; b < int(bodies.size()); ++b)
		m_vbdAcceptedPositions[bodies[b]] =
			std::vector<btVector3>(vbd.x.begin() + offsets[b], vbd.x.begin() + offsets[b] + bodies[b]->m_nodes.size());
	if (vbd.recoveredIntersections)
		btDeformableDiagnostics::write("VBD_RECOVERY", "bounded_initial_overlap=1 max_distance=%.9g",
									   double(vbd.settings.recoveryDistance));
	for (const auto &range : jointRanges)
	{
		btScalar maximum = 0;
		for (int row = 0; row < range.count; ++row)
			maximum = btMax(maximum, btFabs(vbd.joints[range.first + row].impulse) / range.scales[row]);
		if (auto *feedback = range.constraint->getJointFeedback())
		{
			feedback->m_appliedForceBodyA.setZero();
			feedback->m_appliedTorqueBodyA.setZero();
			feedback->m_appliedForceBodyB.setZero();
			feedback->m_appliedTorqueBodyB.setZero();
			for (int row = 0; row < range.count; ++row)
			{
				const auto &r = vbd.joints[range.first + row];
				feedback->m_appliedForceBodyA += r.linearA * (r.impulse / (timeStep * scale));
				feedback->m_appliedForceBodyB += r.linearB * (r.impulse / (timeStep * scale));
				feedback->m_appliedTorqueBodyA += r.angularA * (r.impulse / (timeStep * scale * scale));
				feedback->m_appliedTorqueBodyB += r.angularB * (r.impulse / (timeStep * scale * scale));
			}
		}
		range.constraint->internalSetAppliedImpulse(maximum);
		if (maximum >= range.constraint->getBreakingImpulseThreshold())
			range.constraint->setEnabled(false);
	}
	btDeformableDiagnostics::write(
		"VBD_STEP",
		"cd=gimpact soft_bodies=%d rigid_bodies=%d attachments=%d surface=%s surface_vertices=%d surface_triangles=%d grab_springs=%d "
		"iterations=%d workers=%d colors=%d contacts=%d planes=%d pairs=%d triangles=%d radius=%.9g discovery=%.9g minJ=%.9g gpu_guards=%d "
		"gpu_fallbacks=%d gpu_surfaces=%d gpu_surface_fallbacks=%d",
		int(bodies.size()), int(rigidBodies.size()), int(vbd.attachments.size()), vbd.mappedSurface ? "mapped" : "tet_boundary",
		int(vbd.surfaceVertices.size()), int(vbd.surface.size()), int(vbd.springs.size()), vbd.settings.iterations, vbd.settings.workers,
		vbd.colorCount, int(vbd.contacts.size()), int(vbd.planes.size()), vbd.broadphasePairs, int(vbd.barriers.size()),
		double(vbd.settings.radius), double(vbd.settings.gap), double(vbd.minimumJ), vbd.gpuGuardCalls, vbd.gpuGuardFallbacks,
		vbd.gpuSurfaceCalls, vbd.gpuSurfaceFallbacks);
	if (!vbd.gpuError.empty())
		btDeformableDiagnostics::write("VBD_GPU_FALLBACK", "reason=%s", vbd.gpuError.c_str());
	btVector3 center(0, 0, 0), minimum = vbd.x[0], maximum = vbd.x[0];
	btScalar totalMass = 0;
	for (int n = 0; n < int(vbd.x.size()); ++n)
	{
		center += vbd.x[n] * vbd.mass[n];
		totalMass += vbd.mass[n];
		minimum.setMin(vbd.x[n]);
		maximum.setMax(vbd.x[n]);
	}
	if (totalMass > 0)
		center /= totalMass;
	btDeformableDiagnostics::write("VBD_POSITION", "units=m cx=%.9g cy=%.9g cz=%.9g min_z=%.9g height=%.9g", double(center.x()),
								   double(center.y()), double(center.z()), double(minimum.z()), double(maximum.z() - minimum.z()));
	btDeformableDiagnostics::bodies("end", m_softBodies, true);
	if (m_internalTickCallback)
		(*m_internalTickCallback)(this, timeStep);
	++m_executed_step_counter;
	m_dispatchInfo.m_stepCounter = m_executed_step_counter;
}


btDeformableMultiBodyDynamicsWorld::~btDeformableMultiBodyDynamicsWorld()
{
	delete m_solverDeformableBodyIslandCallback;
}

void btDeformableMultiBodyDynamicsWorld::setVbdSolver(bool enabled, btScalar metersPerUnit, int iterations)
{
	if (m_vbdEnabled && !enabled)
	{
		// VBD never maintains legacy recovery history. Seed it from the accepted pose.
		for (int b = 0; b < m_softBodies.size(); ++b)
			m_softBodies[b]->updateLastSafeWorldTransform();
		for (int b = 0; b < m_nonStaticRigidBodies.size(); ++b)
			m_nonStaticRigidBodies[b]->updateLastSafeWorldTransform();
	}
	if (m_vbdEnabled != enabled || m_vbdLengthScale != metersPerUnit)
	{
		m_vbdStepFailed = false;
		m_vbdAcceptedPositions.clear();
		m_vbdMappingCache.reset();
		m_vbdCollisionCache.reset();
	}
	m_vbdEnabled = enabled;
	m_vbdLengthScale = metersPerUnit;
	m_vbdIterations = iterations;
}
