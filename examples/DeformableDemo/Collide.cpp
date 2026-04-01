/*
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

#include "Collide.h"
///btBulletDynamicsCommon.h is the main Bullet include file, contains most common include files.
#include "btBulletDynamicsCommon.h"
#include "BulletSoftBody/btDeformableMultiBodyDynamicsWorld.h"
#include "BulletSoftBody/btSoftBody.h"
#include "BulletSoftBody/btSoftBodyHelpers.h"
#include "BulletSoftBody/btDeformableBodySolver.h"
#include "BulletSoftBody/btSoftBodyRigidBodyCollisionConfiguration.h"
#include "BulletDynamics/Featherstone/btMultiBodyConstraintSolver.h"
#include "../CommonInterfaces/CommonParameterInterface.h"
#include "../OpenGLWindow/SimpleCamera.h"
#include <stdio.h>  //printf debugging

#include "../CommonInterfaces/CommonDeformableBodyBase.h"
#include "../Utils/b3ResourcePath.h"

///The Collide shows the contact between volumetric deformable objects and rigid objects.
static btScalar E = 68000000;
static btScalar nu = 0.3;
static btScalar damping_alpha = 0.1;
static btScalar damping_beta = 0.01;
static btScalar COLLIDING_VELOCITY = 15;
static const float CAMERA_FRUSTUM_NEAR = 1.f;
static const float CAMERA_FRUSTUM_FAR = 1000000.f;
static const bool kUseImplicitExperiment = true;
static const bool kUseProjectionExperiment = true;
static const bool kUseLineSearchExperiment = false;
static const float kInternalTimeStepExperiment = 0.001f;
static const int kMaxSubstepsExperiment = 1;
static const int kImplicitNewtonIterationsExperiment = 3;
static const float kImplicitNewtonToleranceExperiment = 1e-5f;
static const bool kEnableProgrammaticGrabExperiment = false;
static const float kProgrammaticGrabStartTime = 0.5f;
static const float kProgrammaticGrabDuration = 0.5f;
static const btVector3 kProgrammaticGrabOffset(0, 0, 150);

namespace FooSpace
{
struct TetraCube
{
#include "../SoftDemo/cantilever.inl"
};
}  // namespace FooSpace

class Collide : public CommonDeformableBodyBase
{
	btDeformableLinearElasticityForce* m_linearElasticity;
	btRigidBody* m_tipLoad;
	btScalar m_elapsedTime;
	bool m_scriptedRigidPickActive;
	btPoint2PointConstraint* m_scriptedRigidPickConstraint;
	btVector3 m_scriptedRigidPickTarget;

public:
	Collide(struct GUIHelperInterface* helper)
		: CommonDeformableBodyBase(helper)
	{
		m_linearElasticity = 0;
		m_tipLoad = 0;
		m_elapsedTime = 0;
		m_scriptedRigidPickActive = false;
		m_scriptedRigidPickConstraint = 0;
		m_scriptedRigidPickTarget = btVector3(0, 0, 0);
		m_pickingForceElasticStiffness = 50000;
		m_pickingForceDampingStiffness = 100;
		m_maxPickingForce = 1e10;
		m_rigidPickingImpulseClamp = 20000;
		m_rigidPickingTau = 0.2f;
		m_logPicking = true;
	}
	virtual ~Collide()
	{
	}
	void initPhysics();

	void exitPhysics();

	void resetCamera()
	{
		float dist = 1000;
		float pitch = 0;
		float yaw = 0;
		float targetPos[3] = {0, 0, 0};
		m_guiHelper->resetCamera(dist, yaw, pitch, targetPos[0], targetPos[1], targetPos[2]);

		if (m_guiHelper->getRenderInterface() && m_guiHelper->getRenderInterface()->getActiveCamera())
		{
			SimpleCamera* camera = static_cast<SimpleCamera*>(m_guiHelper->getRenderInterface()->getActiveCamera());
			camera->setCameraFrustumNear(CAMERA_FRUSTUM_NEAR);
			camera->setCameraFrustumFar(CAMERA_FRUSTUM_FAR);
		}
	}

	btRigidBody* Ctor_RbUpStack()
	{
		float mass = 100.0;
		btCollisionShape* shape = new btBoxShape(btVector3(10, 10, 10));
		btTransform startTransform;
		startTransform.setIdentity();
		startTransform.setOrigin(btVector3(-1020, 0, 0));
		btRigidBody* rb = createRigidBody(mass, startTransform, shape);
		//rb->setLinearVelocity(btVector3(0, +COLLIDING_VELOCITY, 0));
		return rb;
	}

	void stepSimulation(float deltaTime)
	{
		m_elapsedTime += deltaTime;
		m_linearElasticity->setPoissonRatio(nu);
		m_linearElasticity->setYoungsModulus(E);
		m_linearElasticity->setDamping(damping_alpha, damping_beta);
		if (kEnableProgrammaticGrabExperiment && m_tipLoad)
		{
			const bool shouldGrab = (m_elapsedTime >= kProgrammaticGrabStartTime &&
									 m_elapsedTime < kProgrammaticGrabStartTime + kProgrammaticGrabDuration);
			if (shouldGrab && !m_scriptedRigidPickActive)
			{
				const btVector3 bodyOrigin = m_tipLoad->getWorldTransform().getOrigin();
				const btVector3 localPivot = m_tipLoad->getCenterOfMassTransform().inverse() * bodyOrigin;
				m_scriptedRigidPickConstraint = new btPoint2PointConstraint(*m_tipLoad, localPivot);
				m_scriptedRigidPickConstraint->m_setting.m_impulseClamp = m_rigidPickingImpulseClamp;
				m_scriptedRigidPickConstraint->m_setting.m_tau = m_rigidPickingTau;
				m_dynamicsWorld->addConstraint(m_scriptedRigidPickConstraint, true);
				m_scriptedRigidPickTarget = bodyOrigin + kProgrammaticGrabOffset;
				m_scriptedRigidPickActive = true;
				fprintf(stderr, "scripted rigid pick begin t=%f bodyOrigin=(%f,%f,%f) targetOffset=(%f,%f,%f)\n",
						m_elapsedTime, bodyOrigin.x(), bodyOrigin.y(), bodyOrigin.z(),
						kProgrammaticGrabOffset.x(), kProgrammaticGrabOffset.y(), kProgrammaticGrabOffset.z());
			}
			else if (!shouldGrab && m_scriptedRigidPickActive)
			{
				m_dynamicsWorld->removeConstraint(m_scriptedRigidPickConstraint);
				delete m_scriptedRigidPickConstraint;
				m_scriptedRigidPickConstraint = 0;
				m_scriptedRigidPickActive = false;
				fprintf(stderr, "scripted rigid pick end t=%f\n", m_elapsedTime);
			}
			if (m_scriptedRigidPickActive && m_scriptedRigidPickConstraint)
			{
				const btVector3 currentOrigin = m_tipLoad->getWorldTransform().getOrigin();
				m_scriptedRigidPickConstraint->setPivotB(m_scriptedRigidPickTarget);
				fprintf(stderr, "scripted rigid pick move t=%f target=(%f,%f,%f) origin=(%f,%f,%f)\n",
						m_elapsedTime, m_scriptedRigidPickTarget.x(), m_scriptedRigidPickTarget.y(), m_scriptedRigidPickTarget.z(),
						currentOrigin.x(), currentOrigin.y(), currentOrigin.z());
			}
		}
		m_dynamicsWorld->stepSimulation(deltaTime, kMaxSubstepsExperiment, kInternalTimeStepExperiment);

		btDeformableMultiBodyDynamicsWorld* deformableWorld = getDeformableDynamicsWorld();
		btSoftBody* psb = (btSoftBody*)deformableWorld->getSoftBodyArray()[0];
		if (m_tipLoad)
		{
			const btVector3 bodyOrigin = m_tipLoad->getWorldTransform().getOrigin();
			fprintf(stderr, "collide frame t=%f tipLoadOrigin=(%f,%f,%f) tipLoadLinVel=(%f,%f,%f) freeNode3=(%f,%f,%f)\n",
					m_elapsedTime,
					bodyOrigin.x(), bodyOrigin.y(), bodyOrigin.z(),
					m_tipLoad->getLinearVelocity().x(), m_tipLoad->getLinearVelocity().y(), m_tipLoad->getLinearVelocity().z(),
					psb->m_nodes[3].m_x.x(), psb->m_nodes[3].m_x.y(), psb->m_nodes[3].m_x.z());
		}
	}

	virtual void renderScene()
	{
		CommonDeformableBodyBase::renderScene();
		btDeformableMultiBodyDynamicsWorld* deformableWorld = getDeformableDynamicsWorld();

		for (int i = 0; i < deformableWorld->getSoftBodyArray().size(); i++)
		{
			btSoftBody* psb = (btSoftBody*)deformableWorld->getSoftBodyArray()[i];
			{
				btSoftBodyHelpers::DrawFrame(psb, deformableWorld->getDebugDrawer());
				btSoftBodyHelpers::Draw(psb, deformableWorld->getDebugDrawer(), deformableWorld->getDrawFlags());
			}
		}
	}
};

void Collide::initPhysics()
{
	m_guiHelper->setUpAxis(2);

	///collision configuration contains default setup for memory, collision setup
	m_collisionConfiguration = new btSoftBodyRigidBodyCollisionConfiguration();

	///use the default collision dispatcher. For parallel processing you can use a diffent dispatcher (see Extras/BulletMultiThreaded)
	m_dispatcher = new btCollisionDispatcher(m_collisionConfiguration);

	m_broadphase = new btDbvtBroadphase();
	btDeformableBodySolver* deformableBodySolver = new btDeformableBodySolver();

	btDeformableMultiBodyConstraintSolver* sol = new btDeformableMultiBodyConstraintSolver();
	sol->setDeformableSolver(deformableBodySolver);
	m_solver = sol;

	m_dynamicsWorld = new btDeformableMultiBodyDynamicsWorld(m_dispatcher, m_broadphase, sol, m_collisionConfiguration, deformableBodySolver);
	btVector3 gravity = btVector3(0, 0, 0);
	m_dynamicsWorld->setGravity(gravity);
	m_guiHelper->createPhysicsDebugDrawer(m_dynamicsWorld);

	btSoftBody* psb = nullptr;

	// create volumetric soft body
	{
		psb = btSoftBodyHelpers::CreateFromTetGenData(getDeformableDynamicsWorld()->getWorldInfo(),
																  FooSpace::TetraCube::getElements(),
																  0,
																  FooSpace::TetraCube::getNodes(),
																  false, true, true);
		getDeformableDynamicsWorld()->addSoftBody(psb);
		//psb->scale(btVector3(2, 2, 2));
		//psb->translate(btVector3(0, 7, 0));
		psb->getCollisionShape()->setMargin(0.1);
		psb->setTotalMass(0.4);
		psb->m_cfg.kKHR = 1;  // collision hardness with kinematic objects
		psb->m_cfg.kCHR = 1;  // collision hardness with rigid body
		psb->m_cfg.kDF = 0;
		//psb->m_cfg.collisions = btSoftBody::fCollision::SDF_RD;
		//psb->m_cfg.collisions |= btSoftBody::fCollision::SDF_RDN;
		psb->m_sleepingThreshold = 0;
		
		psb->m_nodes[684].m_frozen = 1;
		psb->m_nodes[685].m_frozen = 1;
		psb->m_nodes[686].m_frozen = 1;
		psb->m_nodes[687].m_frozen = 1;
		psb->m_nodes[688].m_frozen = 1;
		psb->m_nodes[689].m_frozen = 1;
		psb->m_nodes[690].m_frozen = 1;
		psb->m_nodes[691].m_frozen = 1;
		psb->m_nodes[692].m_frozen = 1;

		btSoftBodyHelpers::generateBoundaryFaces(psb);

		//psb->setVelocity(btVector3(0, -COLLIDING_VELOCITY, 0));

		btDeformableLinearElasticityForce* linearElasticity = new btDeformableLinearElasticityForce(100, 100, 0.01);
		m_linearElasticity = linearElasticity;
		getDeformableDynamicsWorld()->addForce(psb, linearElasticity);
		m_forces.push_back(linearElasticity);

		btVector3 gravity = btVector3(0, 0, -10000);
		m_dynamicsWorld->setGravity(gravity);
		getDeformableDynamicsWorld()->getWorldInfo().m_gravity = gravity;

		/*btDeformableGravityForce* gravity_force = new btDeformableGravityForce(gravity);
		getDeformableDynamicsWorld()->addForce(psb, gravity_force);
		m_forces.push_back(gravity_force);*/
	}
	getDeformableDynamicsWorld()->setImplicit(kUseImplicitExperiment);
	getDeformableDynamicsWorld()->setLineSearch(kUseLineSearchExperiment);
	getDeformableDynamicsWorld()->setUseProjection(kUseProjectionExperiment);
	getDeformableDynamicsWorld()->setMaxNewtonIterations(kImplicitNewtonIterationsExperiment);
	getDeformableDynamicsWorld()->setNewtonTolerance(kImplicitNewtonToleranceExperiment);
	fprintf(stderr, "Collide config: implicit=%d projection=%d lineSearch=%d internalTimeStep=%f maxSubsteps=%d newtonIterations=%d newtonTolerance=%g\n",
		kUseImplicitExperiment ? 1 : 0,
		kUseProjectionExperiment ? 1 : 0,
		kUseLineSearchExperiment ? 1 : 0,
		kInternalTimeStepExperiment,
		kMaxSubstepsExperiment,
		kImplicitNewtonIterationsExperiment,
		(double)kImplicitNewtonToleranceExperiment);
	//getDeformableDynamicsWorld()->getSolverInfo().m_deformable_erp = 0.3;
	//getDeformableDynamicsWorld()->getSolverInfo().m_deformable_maxErrorReduction = btScalar(200);
	getDeformableDynamicsWorld()->getSolverInfo().m_leastSquaresResidualThreshold = 1e-3;
	//getDeformableDynamicsWorld()->getSolverInfo().m_splitImpulse = true;
	getDeformableDynamicsWorld()->getSolverInfo().m_numIterations = 500;
	// add a few rigid bodies
	auto rb = Ctor_RbUpStack();
	m_tipLoad = rb;
	fprintf(stderr, "Collide picking config: softElastic=%f softDamping=%f softMaxForce=%f rigidClamp=%f rigidTau=%f scriptedGrab=%d\n",
		m_pickingForceElasticStiffness,
		m_pickingForceDampingStiffness,
		m_maxPickingForce,
		m_rigidPickingImpulseClamp,
		m_rigidPickingTau,
		kEnableProgrammaticGrabExperiment ? 1 : 0);

	psb->setIgnoreCollisionCheck(rb, true);
	psb->appendDeformableAnchor(0, rb);
	psb->appendDeformableAnchor(2, rb);
	psb->appendDeformableAnchor(3, rb);
	psb->appendDeformableAnchor(6, rb);
	psb->appendDeformableAnchor(8, rb);
	psb->appendDeformableAnchor(10, rb);
	psb->appendDeformableAnchor(12, rb);
	psb->appendDeformableAnchor(14, rb);
	psb->appendDeformableAnchor(16, rb);
	/*psb->m_nodes[0].m_frozen = 1;
		psb->m_nodes[2].m_frozen = 1;
		psb->m_nodes[3].m_frozen = 1;
		psb->m_nodes[6].m_frozen = 1;
		psb->m_nodes[8].m_frozen = 1;
		psb->m_nodes[10].m_frozen = 1;
		psb->m_nodes[12].m_frozen = 1;
		psb->m_nodes[14].m_frozen = 1;
		psb->m_nodes[16].m_frozen = 1;*/

	m_guiHelper->autogenerateGraphicsObjects(m_dynamicsWorld);

	//	{
	//		SliderParams slider("Young's Modulus", &E);
	//		slider.m_minVal = 0;
	//		slider.m_maxVal = 2000;
	//		if (m_guiHelper->getParameterInterface())
	//			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(slider);
	//	}
	//	{
	//		SliderParams slider("Poisson Ratio", &nu);
	//		slider.m_minVal = 0.05;
	//		slider.m_maxVal = 0.49;
	//		if (m_guiHelper->getParameterInterface())
	//			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(slider);
	//	}
	//	{
	//		SliderParams slider("Mass Damping", &damping_alpha);
	//		slider.m_minVal = 0.001;
	//		slider.m_maxVal = 0.01;
	//		if (m_guiHelper->getParameterInterface())
	//			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(slider);
	//	}
	//    {
	//        SliderParams slider("Stiffness Damping", &damping_beta);
	//        slider.m_minVal = 0.001;
	//        slider.m_maxVal = 0.01;
	//        if (m_guiHelper->getParameterInterface())
	//            m_guiHelper->getParameterInterface()->registerSliderFloatParameter(slider);
	//    }
	{
		SliderParams slider("Young's Modulus", &E);
		slider.m_minVal = 0;
		slider.m_maxVal = 100'000'000;
		if (m_guiHelper->getParameterInterface())
			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(slider);
	}
	{
		SliderParams slider("Poisson Ratio", &nu);
		slider.m_minVal = 0.05;
		slider.m_maxVal = 0.49;
		if (m_guiHelper->getParameterInterface())
			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(slider);
	}
	{
		SliderParams slider("Mass Damping", &damping_alpha);
		slider.m_minVal = 0;
		slider.m_maxVal = 1;
		if (m_guiHelper->getParameterInterface())
			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(slider);
	}
	{
		SliderParams slider("Stiffness Damping", &damping_beta);
		slider.m_minVal = 0;
		slider.m_maxVal = 0.1;
		if (m_guiHelper->getParameterInterface())
			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(slider);
	}
}

void Collide::exitPhysics()
{
	//cleanup in the reverse order of creation/initialization
	if (m_scriptedRigidPickConstraint)
	{
		m_dynamicsWorld->removeConstraint(m_scriptedRigidPickConstraint);
		delete m_scriptedRigidPickConstraint;
		m_scriptedRigidPickConstraint = 0;
		m_scriptedRigidPickActive = false;
	}
	removePickingConstraint();
	//remove the rigidbodies from the dynamics world and delete them
	int i;
	for (i = m_dynamicsWorld->getNumCollisionObjects() - 1; i >= 0; i--)
	{
		btCollisionObject* obj = m_dynamicsWorld->getCollisionObjectArray()[i];
		btRigidBody* body = btRigidBody::upcast(obj);
		if (body && body->getMotionState())
		{
			delete body->getMotionState();
		}
		m_dynamicsWorld->removeCollisionObject(obj);
		delete obj;
	}
	// delete forces
	for (int j = 0; j < m_forces.size(); j++)
	{
		btDeformableLagrangianForce* force = m_forces[j];
		delete force;
	}
	m_forces.clear();

	//delete collision shapes
	for (int j = 0; j < m_collisionShapes.size(); j++)
	{
		btCollisionShape* shape = m_collisionShapes[j];
		delete shape;
	}
	m_collisionShapes.clear();

	delete m_dynamicsWorld;

	delete m_solver;

	delete m_broadphase;

	delete m_dispatcher;

	delete m_collisionConfiguration;
}

class CommonExampleInterface* CollideCreateFunc(struct CommonExampleOptions& options)
{
	return new Collide(options.m_guiHelper);
}
