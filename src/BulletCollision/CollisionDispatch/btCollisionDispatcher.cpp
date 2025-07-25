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

#include "btCollisionDispatcher.h"
#include "LinearMath/btQuickprof.h"

#include "BulletCollision/BroadphaseCollision/btCollisionAlgorithm.h"

#include "BulletCollision/CollisionShapes/btCollisionShape.h"
#include "BulletCollision/CollisionDispatch/btCollisionObject.h"
#include "BulletCollision/BroadphaseCollision/btOverlappingPairCache.h"
#include "LinearMath/btPoolAllocator.h"
#include "BulletCollision/CollisionDispatch/btCollisionConfiguration.h"
#include "BulletCollision/CollisionDispatch/btCollisionObjectWrapper.h"

#ifdef BT_DEBUG
#include <stdio.h>
#endif

btCollisionDispatcher::btCollisionDispatcher(btCollisionConfiguration* collisionConfiguration) : m_dispatcherFlags(btCollisionDispatcher::CD_USE_RELATIVE_CONTACT_BREAKING_THRESHOLD),
																								 m_collisionConfiguration(collisionConfiguration),
																								 m_aabbSimilarityThreshold(0.9)
{
	int i;

	setNearCallback(defaultNearCallback);
	m_ptrForProgressReportCallback = nullptr;
	m_progressReportCallback = nullptr;

	m_collisionAlgorithmPoolAllocator = collisionConfiguration->getCollisionAlgorithmPool();

	m_persistentManifoldPoolAllocator = collisionConfiguration->getPersistentManifoldPool();

	for (i = 0; i < MAX_BROADPHASE_COLLISION_TYPES; i++)
	{
		for (int j = 0; j < MAX_BROADPHASE_COLLISION_TYPES; j++)
		{
			m_doubleDispatchContactPoints[i][j] = m_collisionConfiguration->getCollisionAlgorithmCreateFunc(i, j);
			btAssert(m_doubleDispatchContactPoints[i][j]);
			m_doubleDispatchClosestPoints[i][j] = m_collisionConfiguration->getClosestPointsAlgorithmCreateFunc(i, j);
		}
	}
}

void btCollisionDispatcher::registerCollisionCreateFunc(int proxyType0, int proxyType1, btCollisionAlgorithmCreateFunc* createFunc)
{
	m_doubleDispatchContactPoints[proxyType0][proxyType1] = createFunc;
}

void btCollisionDispatcher::registerClosestPointsCreateFunc(int proxyType0, int proxyType1, btCollisionAlgorithmCreateFunc* createFunc)
{
	m_doubleDispatchClosestPoints[proxyType0][proxyType1] = createFunc;
}

btCollisionDispatcher::~btCollisionDispatcher()
{
}

btPersistentManifold* btCollisionDispatcher::getNewManifold(const btCollisionObject* body0, const btCollisionObject* body1, size_t unlimitedSizeManifoldHint)
{
	//btAssert(gNumManifold < 65535);

	//optional relative contact breaking threshold, turned on by default (use setDispatcherFlags to switch off feature for improved performance)

	btScalar contactBreakingThreshold = (m_dispatcherFlags & btCollisionDispatcher::CD_USE_RELATIVE_CONTACT_BREAKING_THRESHOLD) ? btMin(body0->getCollisionShape()->getContactBreakingThreshold(gContactBreakingThreshold), body1->getCollisionShape()->getContactBreakingThreshold(gContactBreakingThreshold))
																																: gContactBreakingThreshold;

	btScalar contactProcessingThreshold = btMin(body0->getContactProcessingThreshold(), body1->getContactProcessingThreshold());

	void* mem = m_persistentManifoldPoolAllocator->allocate(sizeof(btPersistentManifold));
	if (NULL == mem)
	{
		//we got a pool memory overflow, by default we fallback to dynamically allocate memory. If we require a contiguous contact pool then assert.
		if ((m_dispatcherFlags & CD_DISABLE_CONTACTPOOL_DYNAMIC_ALLOCATION) == 0)
		{
			mem = btAlignedAlloc(sizeof(btPersistentManifold), 16);
		}
		else
		{
			btAssert(0);
			//make sure to increase the m_defaultMaxPersistentManifoldPoolSize in the btDefaultCollisionConstructionInfo/btDefaultCollisionConfiguration
			return 0;
		}
	}
	bool onlyGatherContactCounts0 = body0->getCollisionFlags() & btCollisionObject::CF_ONLY_GATHER_CONTACT_COUNTS;
	bool onlyGatherContactCounts1 = body1->getCollisionFlags() & btCollisionObject::CF_ONLY_GATHER_CONTACT_COUNTS;
	bool unlimitedSizeManifold = body0->getInternalType() == btCollisionObject::CO_SOFT_BODY || body1->getInternalType() == btCollisionObject::CO_SOFT_BODY;
	bool onlyGatherContactCounts = onlyGatherContactCounts0 || onlyGatherContactCounts1;
	btPersistentManifold* manifold = new (mem) btPersistentManifold(body0, body1, 0, contactBreakingThreshold, contactProcessingThreshold, unlimitedSizeManifold, unlimitedSizeManifoldHint, onlyGatherContactCounts);
	manifold->m_index1a = m_manifoldsPtr.size();
	m_manifoldsPtr.push_back(manifold);

	return manifold;
}

void btCollisionDispatcher::clearManifold(btPersistentManifold* manifold)
{
	manifold->clearManifold();
}

void btCollisionDispatcher::releaseManifold(btPersistentManifold* manifold)
{
	//printf("releaseManifold: gNumManifold %d\n",gNumManifold);
	clearManifold(manifold);

	int findIndex = manifold->m_index1a;
	btAssert(findIndex < m_manifoldsPtr.size());
	m_manifoldsPtr.swap(findIndex, m_manifoldsPtr.size() - 1);
	m_manifoldsPtr[findIndex]->m_index1a = findIndex;
	m_manifoldsPtr.pop_back();

	manifold->~btPersistentManifold();
	if (m_persistentManifoldPoolAllocator->validPtr(manifold))
	{
		m_persistentManifoldPoolAllocator->freeMemory(manifold);
	}
	else
	{
		btAlignedFree(manifold);
	}
}

btCollisionAlgorithm* btCollisionDispatcher::findAlgorithm(const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, btPersistentManifold* sharedManifold, ebtDispatcherQueryType algoType)
{
	btCollisionAlgorithmConstructionInfo ci;

	ci.m_dispatcher1 = this;
	ci.m_manifold = sharedManifold;
	btCollisionAlgorithm* algo = 0;
	if (algoType == BT_CONTACT_POINT_ALGORITHMS)
	{
		algo = m_doubleDispatchContactPoints[body0Wrap->getCollisionShape()->getShapeType()][body1Wrap->getCollisionShape()->getShapeType()]->CreateCollisionAlgorithm(ci, body0Wrap, body1Wrap);
	}
	else
	{
		algo = m_doubleDispatchClosestPoints[body0Wrap->getCollisionShape()->getShapeType()][body1Wrap->getCollisionShape()->getShapeType()]->CreateCollisionAlgorithm(ci, body0Wrap, body1Wrap);
	}

	return algo;
}

btCollisionAlgorithm* btCollisionDispatcher::findAlgorithm(const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, const int body0ShapeType, const int body1ShapeType, btPersistentManifold* sharedManifold, ebtDispatcherQueryType algoType)
{
	btCollisionAlgorithmConstructionInfo ci;

	ci.m_dispatcher1 = this;
	ci.m_manifold = sharedManifold;
	btCollisionAlgorithm* algo = 0;
	if (algoType == BT_CONTACT_POINT_ALGORITHMS)
	{
		algo = m_doubleDispatchContactPoints[body0ShapeType][body1ShapeType]->CreateCollisionAlgorithm(ci, body0Wrap, body1Wrap);
	}
	else
	{
		algo = m_doubleDispatchClosestPoints[body0ShapeType][body1ShapeType]->CreateCollisionAlgorithm(ci, body0Wrap, body1Wrap);
	}

	return algo;
}

bool btCollisionDispatcher::needsResponse(const btCollisionObject* body0, const btCollisionObject* body1)
{
	//here you can do filtering
	bool hasResponse =
		(body0->hasContactResponse() && body1->hasContactResponse());
	//no response between two static/kinematic bodies:
	hasResponse = hasResponse &&
				  ((!body0->isStaticOrKinematicObject()) || (!body1->isStaticOrKinematicObject()));
	return hasResponse;
}

bool btCollisionDispatcher::needsCollision(const btCollisionObject* body0, const btCollisionObject* body1)
{
	btAssert(body0);
	btAssert(body1);

	bool needsCollision = true;

#ifdef BT_DEBUG
	if (!(m_dispatcherFlags & btCollisionDispatcher::CD_STATIC_STATIC_REPORTED))
	{
		//broadphase filtering already deals with this
		if (body0->isStaticOrKinematicObject() && body1->isStaticOrKinematicObject())
		{
			m_dispatcherFlags |= btCollisionDispatcher::CD_STATIC_STATIC_REPORTED;
			printf("warning btCollisionDispatcher::needsCollision: static-static collision!\n");
		}
	}
#endif  //BT_DEBUG

	if ((!body0->isActive()) && (!body1->isActive()))
		needsCollision = false;
	else if ((!body0->checkCollideWith(body1)) || (!body1->checkCollideWith(body0)))
		needsCollision = false;

	return needsCollision;
}

bool btCollisionDispatcher::needsCollisionUsingAABBSimilarity(btBroadphaseProxy* proxy0, btBroadphaseProxy* proxy1)
{
	if (m_dispatcherFlags & btCollisionDispatcher::CD_USE_AABB_SIMILARITY_THRESHOLD)
	{
		// IoU. Intersection over Union metric.

		const btVector3& aabbMin0 = proxy0->m_aabbMin;
		const btVector3& aabbMax0 = proxy0->m_aabbMax;

		const btVector3& aabbMin1 = proxy1->m_aabbMin;
		const btVector3& aabbMax1 = proxy1->m_aabbMax;

		btVector3 interMin(btMax(aabbMin0.x(), aabbMin1.x()),
						   btMax(aabbMin0.y(), aabbMin1.y()),
						   btMax(aabbMin0.z(), aabbMin1.z()));

		btVector3 interMax(btMin(aabbMax0.x(), aabbMax1.x()),
						   btMin(aabbMax0.y(), aabbMax1.y()),
						   btMin(aabbMax0.z(), aabbMax1.z()));

		btVector3 interSize = interMax - interMin;
		btScalar interVol = interSize.x() * interSize.y() * interSize.z();

		btVector3 size1 = aabbMax0 - aabbMin0;
		btVector3 size2 = aabbMax1 - aabbMin1;
		btScalar vol1 = size1.x() * size1.y() * size1.z();
		btScalar vol2 = size2.x() * size2.y() * size2.z();

		btScalar iou = interVol / (vol1 + vol2 - interVol);

		return iou >= m_aabbSimilarityThreshold;
	}
	return true;
}

///interface for iterating all overlapping collision pairs, no matter how those pairs are stored (array, set, map etc)
///this is useful for the collision dispatcher.
class btCollisionPairCallback : public btOverlapCallback
{
	const btDispatcherInfo& m_dispatchInfo;
	btCollisionDispatcher* m_dispatcher;

public:
	btCollisionPairCallback(const btDispatcherInfo& dispatchInfo, btCollisionDispatcher* dispatcher)
		: m_dispatchInfo(dispatchInfo),
		  m_dispatcher(dispatcher)
	{
	}

	/*btCollisionPairCallback& operator=(btCollisionPairCallback& other)
	{
		m_dispatchInfo = other.m_dispatchInfo;
		m_dispatcher = other.m_dispatcher;
		return *this;
	}
	*/

	virtual ~btCollisionPairCallback() {}

	virtual bool processOverlap(btBroadphasePair& pair)
	{
		if (m_dispatcher->getProgressReportCallback())
			(*m_dispatcher->getProgressReportCallback())(m_dispatcher->getPtrForProgressReportCallback(), m_overlap_index, m_total_overlap_count);
		return (*m_dispatcher->getNearCallback())(pair, *m_dispatcher, m_dispatchInfo);
	}
};

void btCollisionDispatcher::dispatchAllCollisionPairs(btOverlappingPairCache* pairCache, const btDispatcherInfo& dispatchInfo, btDispatcher* dispatcher)
{
	//m_blockedForChanges = true;

	for (auto& time : previouslyConsumedTime)
	{
		std::get<1>(time.second) = false;  // out of date
	}

	btCollisionPairCallback collisionCallback(dispatchInfo, this);

	{
		BT_PROFILE("processAllOverlappingPairs");
		initialCollisionParticipants.clear();
		initialCollisionParticipants0.clear();
		initialCollisionParticipants1.clear();
		//printf("processAllOverlappingPairs START\n");
		pairCache->processAllOverlappingPairs(&collisionCallback, dispatcher, dispatchInfo);
		//printf("processAllOverlappingPairs END\n");
	}

	//m_blockedForChanges = false;
}

//by default, Bullet will use this near callback
bool btCollisionDispatcher::defaultNearCallback(btBroadphasePair& collisionPair, btCollisionDispatcher& dispatcher, const btDispatcherInfo& dispatchInfo)
{
	btCollisionObject* colObj0 = (btCollisionObject*)collisionPair.m_pProxy0->m_clientObject;
	btCollisionObject* colObj1 = (btCollisionObject*)collisionPair.m_pProxy1->m_clientObject;

	bool needsCollision = dispatcher.needsCollisionUsingAABBSimilarity(collisionPair.m_pProxy0, collisionPair.m_pProxy1) && dispatcher.needsCollision(colObj0, colObj1);
	if (needsCollision)
	{
		btCollisionObjectWrapper obj0Wrap(0, colObj0->getCollisionShape(), colObj0, colObj0->getWorldTransform(), -1, -1);
		btCollisionObjectWrapper obj1Wrap(0, colObj1->getCollisionShape(), colObj1, colObj1->getWorldTransform(), -1, -1);

		//dispatcher will keep algorithms persistent in the collision pair
		if (!collisionPair.m_algorithm)
		{
			collisionPair.m_algorithm = dispatcher.findAlgorithm(&obj0Wrap, &obj1Wrap, 0, BT_CONTACT_POINT_ALGORITHMS);
		}

		if (collisionPair.m_algorithm)
		{
			btManifoldResult contactPointResult(&obj0Wrap, &obj1Wrap);

			if (dispatchInfo.m_dispatchFunc == btDispatcherInfo::DISPATCH_DISCRETE)
			{
				//discrete collision detection query
				collisionPair.m_algorithm->processCollision(&obj0Wrap, &obj1Wrap, dispatchInfo, &contactPointResult);
			}
			else
			{
				//continuous collision detection query, time of impact (toi)
				btScalar toi = collisionPair.m_algorithm->calculateTimeOfImpact(colObj0, colObj1, dispatchInfo, &contactPointResult);
				if (dispatchInfo.m_timeOfImpact > toi)
					dispatchInfo.m_timeOfImpact = toi;
			}
		}
	}
	return !needsCollision;
}

void* btCollisionDispatcher::allocateCollisionAlgorithm(int size)
{
	void* mem = m_collisionAlgorithmPoolAllocator->allocate(size);
	if (NULL == mem)
	{
		//warn user for overflow?
		return btAlignedAlloc(static_cast<size_t>(size), 16);
	}
	return mem;
}

void btCollisionDispatcher::freeCollisionAlgorithm(void* ptr)
{
	if (m_collisionAlgorithmPoolAllocator->validPtr(ptr))
	{
		m_collisionAlgorithmPoolAllocator->freeMemory(ptr);
	}
	else
	{
		btAlignedFree(ptr);
	}
}
