/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/

#include "btCollisionObject.h"
#include "LinearMath/btSerializer.h"
#include "BulletCollision/BroadphaseCollision/btBroadphaseProxy.h"

btScalar btCollisionObject::gFrictionOverride = -1.0;

btCollisionObject::btCollisionObject()
	: m_interpolationLinearVelocity(0.f, 0.f, 0.f),
	  m_interpolationAngularVelocity(0.f, 0.f, 0.f),
	  m_anisotropicFriction(1.f, 1.f, 1.f),
	  m_hasAnisotropicFriction(false),
	  m_contactProcessingThreshold(BT_LARGE_FLOAT),
	  m_broadphaseHandle(0),
	  m_collisionShape(0),
	  m_extensionPointer(0),
	  m_rootCollisionShape(0),
	  m_collisionFlags(btCollisionObject::CF_STATIC_OBJECT),
	  m_islandTag1(-1),
	  m_islandTag2(-1),
	  m_companionId(-1),
	  m_worldArrayIndex(-1),
	  m_activationState1(1),
	  m_deactivationTime(btScalar(0.)),
	  m_friction(btScalar(0.5)),
	  m_restitution(btScalar(0.)),
	  m_rollingFriction(0.0f),
	  m_spinningFriction(0.f),
	  m_contactDamping(.1),
	  m_contactStiffness(BT_LARGE_FLOAT),
	  m_internalType(CO_COLLISION_OBJECT),
	  m_userObjectPointer(0),
	  m_userIndex2(-1),
	  m_userIndex(-1),
	  m_userIndex3(-1),
	  m_hitFraction(btScalar(1.)),
	  m_ccdSweptSphereRadius(btScalar(0.)),
	  m_ccdMotionThreshold(btScalar(0.)),
	  m_checkCollideWith(false),
	  m_updateRevision(0),
	  m_lastSafeApplyCounter(0)
{
	m_worldTransform.setIdentity();
	m_lastSafeWorldTransform.setIdentity();
	m_interpolationWorldTransform.setIdentity();
}

btCollisionObject::~btCollisionObject()
{
	for (auto ref : m_anchorRefs)
		ref->resetColObjPtrsInAnchors(this);
}

void btCollisionObject::setActivationState(int newState) const
{
	if ((m_activationState1 != DISABLE_DEACTIVATION) && (m_activationState1 != DISABLE_SIMULATION))
		m_activationState1 = newState;
}

void btCollisionObject::forceActivationState(int newState) const
{
	m_activationState1 = newState;
}

void btCollisionObject::activate(bool forceActivation) const
{
	if (forceActivation || !(m_collisionFlags & (CF_STATIC_OBJECT | CF_KINEMATIC_OBJECT)))
	{
		setActivationState(ACTIVE_TAG);
		m_deactivationTime = btScalar(0.);
	}
}

const char* btCollisionObject::serialize(void* dataBuffer, btSerializer* serializer) const
{
	btCollisionObjectData* dataOut = (btCollisionObjectData*)dataBuffer;

	m_worldTransform.serialize(dataOut->m_worldTransform);
	m_interpolationWorldTransform.serialize(dataOut->m_interpolationWorldTransform);
	m_interpolationLinearVelocity.serialize(dataOut->m_interpolationLinearVelocity);
	m_interpolationAngularVelocity.serialize(dataOut->m_interpolationAngularVelocity);
	m_anisotropicFriction.serialize(dataOut->m_anisotropicFriction);
	dataOut->m_hasAnisotropicFriction = m_hasAnisotropicFriction;
	dataOut->m_contactProcessingThreshold = m_contactProcessingThreshold;
	dataOut->m_broadphaseHandle = 0;
	dataOut->m_collisionShape = serializer->getUniquePointer(m_collisionShape);
	dataOut->m_rootCollisionShape = 0;  //@todo
	dataOut->m_collisionFlags = m_collisionFlags;
	dataOut->m_islandTag1 = m_islandTag1;
	dataOut->m_companionId = m_companionId;
	dataOut->m_activationState1 = m_activationState1;
	dataOut->m_deactivationTime = m_deactivationTime;
	dataOut->m_friction = m_friction;
	dataOut->m_rollingFriction = m_rollingFriction;
	dataOut->m_contactDamping = m_contactDamping;
	dataOut->m_contactStiffness = m_contactStiffness;
	dataOut->m_restitution = m_restitution;
	dataOut->m_internalType = m_internalType;

	char* name = (char*)serializer->findNameForPointer(this);
	dataOut->m_name = (char*)serializer->getUniquePointer(name);
	if (dataOut->m_name)
	{
		serializer->serializeName(name);
	}
	dataOut->m_hitFraction = m_hitFraction;
	dataOut->m_ccdSweptSphereRadius = m_ccdSweptSphereRadius;
	dataOut->m_ccdMotionThreshold = m_ccdMotionThreshold;
	dataOut->m_checkCollideWith = m_checkCollideWith;
	if (m_broadphaseHandle)
	{
		dataOut->m_collisionFilterGroup = m_broadphaseHandle->m_collisionFilterGroup;
		dataOut->m_collisionFilterMask = m_broadphaseHandle->m_collisionFilterMask;
		dataOut->m_uniqueId = m_broadphaseHandle->m_uniqueId;
	}
	else
	{
		dataOut->m_collisionFilterGroup = 0;
		dataOut->m_collisionFilterMask = 0;
		dataOut->m_uniqueId = -1;
	}
	return btCollisionObjectDataName;
}

void btCollisionObject::serializeSingleObject(class btSerializer* serializer) const
{
	int len = calculateSerializeBufferSize();
	btChunk* chunk = serializer->allocate(len, 1);
	const char* structType = serialize(chunk->m_oldPtr, serializer);
	serializer->finalizeChunk(chunk, structType, BT_COLLISIONOBJECT_CODE, (void*)this);
}

void btCollisionObject::applyLastSafeWorldTransform(const std::map<int, btScalar>* partial)
{
	if (getCollisionFlags() & CF_APPLY_LAST_SAFE)
	{
		constexpr btScalar maxApplySteps = 10;
		// We sacrifice few iterations to move to the safe position only gradually. This significantly reduces the jitter of
		// jumping between the safe and stuck positions. The unstuck position will be much closer to the real point of contact.
		btScalar fraction = getLastSafeApplyCounter() / maxApplySteps;
		fraction = std::min(fraction, 1.0);

		btTransform dst = getLastSafeWorldTransform();
		btTransform src = getWorldTransform();

		btVector3 interpOrigin = src.getOrigin().lerp(dst.getOrigin(), fraction);
		btQuaternion interpRot = src.getRotation().slerp(dst.getRotation(), fraction);
		btTransform interp(interpRot, interpOrigin);
		setWorldTransform(interp);
		if (fraction < 1.0)
			++m_lastSafeApplyCounter;
		// Would it be possible to develop a velocity based unstuck for btRigidBody, where position would be overwritten, but also a linear and angular velocity would be added
		// for one simulation step, which would make it easier for the solvers (the deformable soft body constraint solver in particular) to cope without explosions?
	}
}
