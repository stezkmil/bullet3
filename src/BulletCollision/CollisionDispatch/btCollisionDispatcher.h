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

#ifndef BT_COLLISION__DISPATCHER_H
#define BT_COLLISION__DISPATCHER_H

#include "BulletCollision/BroadphaseCollision/btDispatcher.h"
#include "BulletCollision/NarrowPhaseCollision/btPersistentManifold.h"

#include "BulletCollision/CollisionDispatch/btManifoldResult.h"

#include "BulletCollision/BroadphaseCollision/btBroadphaseProxy.h"
#include "LinearMath/btAlignedObjectArray.h"

class btIDebugDraw;
class btOverlappingPairCache;
class btPoolAllocator;
class btCollisionConfiguration;

#include "btCollisionCreateFunc.h"

#define USE_DISPATCH_REGISTRY_ARRAY 1

class btCollisionDispatcher;
///user can override this nearcallback for collision filtering and more finegrained control over collision detection
typedef bool (*btNearCallback)(btBroadphasePair& collisionPair, btCollisionDispatcher& dispatcher, const btDispatcherInfo& dispatchInfo);

typedef void (*btProgressReportCallback)(void* ptr, int index, int total);

///btCollisionDispatcher supports algorithms that handle ConvexConvex and ConvexConcave collision pairs.
///Time of Impact, Closest Points and Penetration Depth.
class btCollisionDispatcher : public btDispatcher
{
protected:
	int m_dispatcherFlags;

	btAlignedObjectArray<btPersistentManifold*> m_manifoldsPtr;

	btNearCallback m_nearCallback;

	void* m_ptrForProgressReportCallback;
	btProgressReportCallback m_progressReportCallback;

	btPoolAllocator* m_collisionAlgorithmPoolAllocator;

	btPoolAllocator* m_persistentManifoldPoolAllocator;

	btCollisionAlgorithmCreateFunc* m_doubleDispatchContactPoints[MAX_BROADPHASE_COLLISION_TYPES][MAX_BROADPHASE_COLLISION_TYPES];

	btCollisionAlgorithmCreateFunc* m_doubleDispatchClosestPoints[MAX_BROADPHASE_COLLISION_TYPES][MAX_BROADPHASE_COLLISION_TYPES];

	btCollisionConfiguration* m_collisionConfiguration;

	btInitialCollisionParticipants initialCollisionParticipants;
	btInitialCollisionParticipantsSingle initialCollisionParticipants0, initialCollisionParticipants1;
	btPreviouslyFoundPairMap previouslyConsumedTime;

public:
	btScalar m_aabbSimilarityThreshold;

	enum DispatcherFlags
	{
		CD_STATIC_STATIC_REPORTED = 1,
		CD_USE_RELATIVE_CONTACT_BREAKING_THRESHOLD = 2,
		CD_DISABLE_CONTACTPOOL_DYNAMIC_ALLOCATION = 4,
		CD_USE_AABB_SIMILARITY_THRESHOLD = 8
	};

	int getDispatcherFlags() const
	{
		return m_dispatcherFlags;
	}

	void setDispatcherFlags(int flags)
	{
		m_dispatcherFlags = flags;
	}

	///registerCollisionCreateFunc allows registration of custom/alternative collision create functions
	void registerCollisionCreateFunc(int proxyType0, int proxyType1, btCollisionAlgorithmCreateFunc* createFunc);

	void registerClosestPointsCreateFunc(int proxyType0, int proxyType1, btCollisionAlgorithmCreateFunc* createFunc);

	void addInitialCollisionParticipant(btInitialCollisionParticipants::value_type pair) override
	{
		initialCollisionParticipants.insert(pair);
		initialCollisionParticipants0.insert(pair.first);
		initialCollisionParticipants1.insert(pair.second);
	}

	void addPreviouslyConsumedTime(btPreviouslyFoundPairMap::key_type key, btPreviouslyFoundPairMap::mapped_type value) override
	{
		previouslyConsumedTime.insert_or_assign(key, value);
	}

	const btPreviouslyFoundPairMap& getPreviouslyConsumedTime() const override
	{
		return previouslyConsumedTime;
	}

	const btInitialCollisionParticipants& getInitialCollisionParticipants() const override
	{
		return initialCollisionParticipants;
	}

	const btInitialCollisionParticipantsSingle& getInitialCollisionParticipants0() const override
	{
		return initialCollisionParticipants0;
	}

	const btInitialCollisionParticipantsSingle& getInitialCollisionParticipants1() const override
	{
		return initialCollisionParticipants1;
	}

	int getNumManifolds() const
	{
		return int(m_manifoldsPtr.size());
	}

	btPersistentManifold** getInternalManifoldPointer()
	{
		return m_manifoldsPtr.size() ? &m_manifoldsPtr[0] : 0;
	}

	btPersistentManifold* getManifoldByIndexInternal(int index)
	{
		btAssert(index >= 0);
		btAssert(index < m_manifoldsPtr.size());
		return m_manifoldsPtr[index];
	}

	const btPersistentManifold* getManifoldByIndexInternal(int index) const
	{
		btAssert(index >= 0);
		btAssert(index < m_manifoldsPtr.size());
		return m_manifoldsPtr[index];
	}

	btCollisionDispatcher(btCollisionConfiguration* collisionConfiguration);

	virtual ~btCollisionDispatcher();

	virtual btPersistentManifold* getNewManifold(const btCollisionObject* b0, const btCollisionObject* b1, size_t unlimitedSizeManifoldHint = 10);

	virtual void releaseManifold(btPersistentManifold* manifold);

	virtual void clearManifold(btPersistentManifold* manifold);

	btCollisionAlgorithm* findAlgorithm(const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, btPersistentManifold* sharedManifold, ebtDispatcherQueryType queryType);
	btCollisionAlgorithm* findAlgorithm(const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, const int body0ShapeType, const int body1ShapeType, btPersistentManifold* sharedManifold, ebtDispatcherQueryType algoType);

	virtual bool needsCollision(const btCollisionObject* body0, const btCollisionObject* body1);

	virtual bool needsCollisionUsingAABBSimilarity(btBroadphaseProxy* proxy0, btBroadphaseProxy* proxy1);

	virtual bool needsResponse(const btCollisionObject* body0, const btCollisionObject* body1);

	virtual void dispatchAllCollisionPairs(btOverlappingPairCache* pairCache, const btDispatcherInfo& dispatchInfo, btDispatcher* dispatcher);

	void setNearCallback(btNearCallback nearCallback)
	{
		m_nearCallback = nearCallback;
	}

	btNearCallback getNearCallback() const
	{
		return m_nearCallback;
	}

	void setProgressReportCallback(void* ptr, btProgressReportCallback progressReportCallback)
	{
		m_ptrForProgressReportCallback = ptr;
		m_progressReportCallback = progressReportCallback;
	}

	btProgressReportCallback getProgressReportCallback() const
	{
		return m_progressReportCallback;
	}

	void* getPtrForProgressReportCallback() const
	{
		return m_ptrForProgressReportCallback;
	}

	//by default, Bullet will use this near callback
	static bool defaultNearCallback(btBroadphasePair& collisionPair, btCollisionDispatcher& dispatcher, const btDispatcherInfo& dispatchInfo);

	virtual void* allocateCollisionAlgorithm(int size);

	virtual void freeCollisionAlgorithm(void* ptr);

	btCollisionConfiguration* getCollisionConfiguration()
	{
		return m_collisionConfiguration;
	}

	const btCollisionConfiguration* getCollisionConfiguration() const
	{
		return m_collisionConfiguration;
	}

	void setCollisionConfiguration(btCollisionConfiguration* config)
	{
		m_collisionConfiguration = config;
	}

	virtual btPoolAllocator* getInternalManifoldPool()
	{
		return m_persistentManifoldPoolAllocator;
	}

	virtual const btPoolAllocator* getInternalManifoldPool() const
	{
		return m_persistentManifoldPoolAllocator;
	}
};

#endif  //BT_COLLISION__DISPATCHER_H
