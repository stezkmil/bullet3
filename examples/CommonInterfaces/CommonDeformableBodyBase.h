
#ifndef COMMON_DEFORMABLE_BODY_SETUP_H
#define COMMON_DEFORMABLE_BODY_SETUP_H
#include "btBulletDynamicsCommon.h"
#include "BulletDynamics/Featherstone/btMultiBodyDynamicsWorld.h"
#include "BulletDynamics/Featherstone/btMultiBodyConstraintSolver.h"
#include "BulletDynamics/Featherstone/btMultiBodyPoint2Point.h"
#include "BulletDynamics/Featherstone/btMultiBodyLinkCollider.h"
#include "btBulletDynamicsCommon.h"
#include "CommonExampleInterface.h"
#include "CommonGUIHelperInterface.h"
#include "CommonRenderInterface.h"
#include "CommonGraphicsAppInterface.h"
#include "CommonWindowInterface.h"
#include "CommonCameraInterface.h"
#include "CommonMultiBodyBase.h"
#include "BulletSoftBody/btSoftBody.h"

struct CommonDeformableBodyBase : public CommonMultiBodyBase
{
	btAlignedObjectArray<btDeformableLagrangianForce*> m_forces;
	btSoftBody* m_pickedSoftBody;
	btDeformableMousePickingForce* m_mouseForce;
	btScalar m_pickingForceElasticStiffness, m_pickingForceDampingStiffness, m_maxPickingForce;
	btScalar m_rigidPickingImpulseClamp, m_rigidPickingTau;
	bool m_logPicking;
	bool m_pickedRigidWasAnchoredCarrier;
	CommonDeformableBodyBase(GUIHelperInterface* helper)
		: CommonMultiBodyBase(helper),
		  m_pickedSoftBody(0),
		  m_mouseForce(0),
		  m_pickingForceElasticStiffness(100),
		  m_pickingForceDampingStiffness(0.0),
		  m_maxPickingForce(0.3),
		  m_rigidPickingImpulseClamp(30.f),
		  m_rigidPickingTau(0.001f),
		  m_logPicking(false),
		  m_pickedRigidWasAnchoredCarrier(false)
	{
	}

	virtual btDeformableMultiBodyDynamicsWorld* getDeformableDynamicsWorld()
	{
		return (btDeformableMultiBodyDynamicsWorld*)m_dynamicsWorld;
	}

	virtual const btDeformableMultiBodyDynamicsWorld* getDeformableDynamicsWorld() const
	{
		return (btDeformableMultiBodyDynamicsWorld*)m_dynamicsWorld;
	}

	struct ClosestRayResultCallbackWithInfo : public btCollisionWorld::ClosestRayResultCallback
	{
		ClosestRayResultCallbackWithInfo(const btVector3& rayFromWorld, const btVector3& rayToWorld)
			: ClosestRayResultCallback(rayFromWorld, rayToWorld)
		{
		}
		int m_faceId;

		virtual btScalar addSingleResult(btCollisionWorld::LocalRayResult& rayResult, bool normalInWorldSpace)
		{
			//caller already does the filter on the m_closestHitFraction
			btAssert(rayResult.m_hitFraction <= m_closestHitFraction);

			m_closestHitFraction = rayResult.m_hitFraction;
			m_collisionObject = rayResult.m_collisionObject;
			if (rayResult.m_localShapeInfo)
			{
				m_faceId = rayResult.m_localShapeInfo->m_triangleIndex;
			}
			else
			{
				m_faceId = -1;
			}
			if (normalInWorldSpace)
			{
				m_hitNormalWorld = rayResult.m_hitNormalLocal;
			}
			else
			{
				///need to transform normal into worldspace
				m_hitNormalWorld = m_collisionObject->getWorldTransform().getBasis() * rayResult.m_hitNormalLocal;
			}
			m_hitPointWorld.setInterpolate3(m_rayFromWorld, m_rayToWorld, rayResult.m_hitFraction);
			return rayResult.m_hitFraction;
		}
	};

	virtual bool pickBody(const btVector3& rayFromWorld, const btVector3& rayToWorld)
	{
		if (getDeformableDynamicsWorld() == 0)
			return false;
		ClosestRayResultCallbackWithInfo rayCallback(rayFromWorld, rayToWorld);
		getDeformableDynamicsWorld()->rayTest(rayFromWorld, rayToWorld, rayCallback);
		if (rayCallback.hasHit())
		{
			btVector3 pickPos = rayCallback.m_hitPointWorld;
			btRigidBody* body = (btRigidBody*)btRigidBody::upcast(rayCallback.m_collisionObject);
			btSoftBody* psb = (btSoftBody*)btSoftBody::upcast(rayCallback.m_collisionObject);
			m_oldPickingPos = rayToWorld;
			m_hitPos = pickPos;
			m_oldPickingDist = (pickPos - rayFromWorld).length();
			if (body)
			{
				if (!(body->isStaticObject() || body->isKinematicObject()))
				{
					m_pickedBody = body;
					m_pickedBody->setActivationState(DISABLE_DEACTIVATION);
					// Soft-anchor carrier bodies feel much better if the mouse pulls near the
					// center of mass instead of injecting mostly torque through an off-center pivot.
					const bool anchoredCarrier = body->hasAnchorRef();
					m_pickedRigidWasAnchoredCarrier = anchoredCarrier;
					btVector3 localPivot = anchoredCarrier ? btVector3(0, 0, 0) : (body->getCenterOfMassTransform().inverse() * pickPos);
					btPoint2PointConstraint* p2p = new btPoint2PointConstraint(*body, localPivot);
					m_dynamicsWorld->addConstraint(p2p, true);
					m_pickedConstraint = p2p;
					p2p->m_setting.m_impulseClamp = m_rigidPickingImpulseClamp;
					p2p->m_setting.m_tau = m_rigidPickingTau;
					if (m_logPicking)
					{
						fprintf(stderr, "pick rigid body=%p mass=%f anchoredCarrier=%d localPivot=(%f,%f,%f) worldHit=(%f,%f,%f) clamp=%f tau=%f\n",
								body,
								body->getInvMass() > 0 ? btScalar(1.0) / body->getInvMass() : btScalar(0.0),
								anchoredCarrier ? 1 : 0,
								localPivot.x(), localPivot.y(), localPivot.z(),
								pickPos.x(), pickPos.y(), pickPos.z(),
								m_rigidPickingImpulseClamp, m_rigidPickingTau);
					}
				}
			}
			else if (psb)
			{
				int face_id = rayCallback.m_faceId;
				if (face_id >= 0 && face_id < psb->m_faces.size())
				{
					m_pickedSoftBody = psb;
					psb->setActivationState(DISABLE_DEACTIVATION);
					const btSoftBody::Face& f = psb->m_faces[face_id];
					btDeformableMousePickingForce* mouse_force = new btDeformableMousePickingForce(m_pickingForceElasticStiffness, m_pickingForceDampingStiffness, &f, nullptr, nullptr, btTransform(btTransform::getIdentity().getBasis(), m_hitPos), m_maxPickingForce);
					m_mouseForce = mouse_force;
					getDeformableDynamicsWorld()->addForce(psb, mouse_force);
					if (m_logPicking)
					{
						fprintf(stderr, "pick soft body=%p face=%d hit=(%f,%f,%f) elastic=%f damping=%f maxForce=%f\n",
								psb, face_id, m_hitPos.x(), m_hitPos.y(), m_hitPos.z(),
								m_pickingForceElasticStiffness, m_pickingForceDampingStiffness, m_maxPickingForce);
					}
				}
			}
			else
			{
				btMultiBodyLinkCollider* multiCol = (btMultiBodyLinkCollider*)btMultiBodyLinkCollider::upcast(rayCallback.m_collisionObject);
				if (multiCol && multiCol->m_multiBody)
				{
					m_prevCanSleep = multiCol->m_multiBody->getCanSleep();
					multiCol->m_multiBody->setCanSleep(false);

					btVector3 pivotInA = multiCol->m_multiBody->worldPosToLocal(multiCol->m_link, pickPos);

					btMultiBodyPoint2Point* p2p = new btMultiBodyPoint2Point(multiCol->m_multiBody, multiCol->m_link, 0, pivotInA, pickPos);
					//if you add too much energy to the system, causing high angular velocities, simulation 'explodes'
					//see also http://www.bulletphysics.org/Bullet/phpBB3/viewtopic.php?f=4&t=949
					//so we try to avoid it by clamping the maximum impulse (force) that the mouse pick can apply
					//it is not satisfying, hopefully we find a better solution (higher order integrator, using joint friction using a zero-velocity target motor with limited force etc?)
					btScalar scaling = 1;
					p2p->setMaxAppliedImpulse(2 * scaling);
					btMultiBodyDynamicsWorld* world = (btMultiBodyDynamicsWorld*)m_dynamicsWorld;
					world->addMultiBodyConstraint(p2p);
					m_pickingMultiBodyPoint2Point = p2p;
				}
			}
		}
		return false;
	}

	virtual bool movePickedBody(const btVector3& rayFromWorld, const btVector3& rayToWorld)
	{
		if (m_pickedBody && m_pickedConstraint)
		{
			btPoint2PointConstraint* pickCon = static_cast<btPoint2PointConstraint*>(m_pickedConstraint);
			if (pickCon)
			{
				//keep it at the same picking distance
				btVector3 newPivotB;
				btVector3 dir = rayToWorld - rayFromWorld;
				dir.normalize();
				dir *= m_oldPickingDist;
				newPivotB = rayFromWorld + dir;
				pickCon->setPivotB(newPivotB);
				if (m_logPicking)
				{
					fprintf(stderr, "move rigid pick target=(%f,%f,%f) bodyOrigin=(%f,%f,%f)\n",
							newPivotB.x(), newPivotB.y(), newPivotB.z(),
							m_pickedBody->getWorldTransform().getOrigin().x(),
							m_pickedBody->getWorldTransform().getOrigin().y(),
							m_pickedBody->getWorldTransform().getOrigin().z());
				}
				return true;
			}
		}
		if (m_pickingMultiBodyPoint2Point)
		{
			//keep it at the same picking distance
			btVector3 dir = rayToWorld - rayFromWorld;
			dir.normalize();
			dir *= m_oldPickingDist;
			btVector3 newPivotB = rayFromWorld + dir;
			m_pickingMultiBodyPoint2Point->setPivotInB(newPivotB);
		}
		if (m_pickedSoftBody && m_mouseForce)
		{
			btVector3 newPivot;
			btVector3 dir = rayToWorld - rayFromWorld;
			dir.normalize();
			dir *= m_oldPickingDist;
			newPivot = rayFromWorld + dir;
			m_mouseForce->setMousePos(newPivot);
			if (m_logPicking)
			{
				fprintf(stderr, "move soft pick target=(%f,%f,%f)\n", newPivot.x(), newPivot.y(), newPivot.z());
			}
		}
		return false;
	}

	virtual void removePickingConstraint()
	{
		if (m_pickedConstraint)
		{
			m_dynamicsWorld->removeConstraint(m_pickedConstraint);

			if (m_pickedBody)
			{
				if (m_pickedRigidWasAnchoredCarrier)
				{
					const btVector3 linearVelocityBefore = m_pickedBody->getLinearVelocity();
					const btVector3 angularVelocityBefore = m_pickedBody->getAngularVelocity();
					m_pickedBody->setLinearVelocity(linearVelocityBefore * btScalar(0.1));
					m_pickedBody->setAngularVelocity(angularVelocityBefore * btScalar(0.1));
					if (m_logPicking)
					{
						fprintf(stderr,
								"release anchored rigid damping body=%p linBefore=(%f,%f,%f) angBefore=(%f,%f,%f) linAfter=(%f,%f,%f) angAfter=(%f,%f,%f)\n",
								m_pickedBody,
								linearVelocityBefore.x(), linearVelocityBefore.y(), linearVelocityBefore.z(),
								angularVelocityBefore.x(), angularVelocityBefore.y(), angularVelocityBefore.z(),
								m_pickedBody->getLinearVelocity().x(), m_pickedBody->getLinearVelocity().y(), m_pickedBody->getLinearVelocity().z(),
								m_pickedBody->getAngularVelocity().x(), m_pickedBody->getAngularVelocity().y(), m_pickedBody->getAngularVelocity().z());
					}
				}
				m_pickedBody->forceActivationState(ACTIVE_TAG);
				m_pickedBody->activate(true);
			}
			delete m_pickedConstraint;
			m_pickedConstraint = 0;
			m_pickedBody = 0;
			m_pickedRigidWasAnchoredCarrier = false;
		}
		if (m_pickingMultiBodyPoint2Point)
		{
			m_pickingMultiBodyPoint2Point->getMultiBodyA()->setCanSleep(m_prevCanSleep);
			btMultiBodyDynamicsWorld* world = (btMultiBodyDynamicsWorld*)m_dynamicsWorld;
			world->removeMultiBodyConstraint(m_pickingMultiBodyPoint2Point);
			delete m_pickingMultiBodyPoint2Point;
			m_pickingMultiBodyPoint2Point = 0;
		}
		if (m_pickedSoftBody)
		{
			getDeformableDynamicsWorld()->removeForce(m_pickedSoftBody, m_mouseForce);
			if (m_logPicking)
			{
				fprintf(stderr, "remove soft pick body=%p\n", m_pickedSoftBody);
			}
			delete m_mouseForce;
			m_mouseForce = 0;
			m_pickedSoftBody = 0;
		}
	}
};
#endif  //COMMON_MULTI_BODY_SETUP_H
