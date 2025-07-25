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

#ifndef BT_DISPATCHER_H
#define BT_DISPATCHER_H
#include "LinearMath/btScalar.h"

#include <set>
#include <map>

class btCollisionAlgorithm;
struct btBroadphaseProxy;
class btRigidBody;
class btCollisionObject;
class btOverlappingPairCache;
struct btCollisionObjectWrapper;

class btPersistentManifold;
class btPoolAllocator;

struct btDispatcherInfo
{
	enum DispatchFunc
	{
		DISPATCH_DISCRETE = 1,
		DISPATCH_CONTINUOUS
	};
	btDispatcherInfo()
		: m_timeStep(btScalar(0.)),
		  m_stepCount(0),
		  m_dispatchFunc(DISPATCH_DISCRETE),
		  m_timeOfImpact(btScalar(1.)),
		  m_useContinuous(true),
		  m_debugDraw(0),
		  m_enableSatConvex(false),
		  m_enableSPU(true),
		  m_useEpa(true),
		  m_allowedCcdPenetration(btScalar(0.04)),
		  m_useConvexConservativeDistanceUtil(false),
		  m_convexConservativeDistanceThreshold(0.0f),
		  m_deterministicOverlappingPairs(false)
	{
	}
	btScalar m_timeStep;
	int m_stepCount;
	int m_dispatchFunc;
	mutable btScalar m_timeOfImpact;
	bool m_useContinuous;
	class btIDebugDraw* m_debugDraw;
	bool m_enableSatConvex;
	bool m_enableSPU;
	bool m_useEpa;
	btScalar m_allowedCcdPenetration;
	bool m_useConvexConservativeDistanceUtil;
	btScalar m_convexConservativeDistanceThreshold;
	bool m_deterministicOverlappingPairs;
};

enum ebtDispatcherQueryType
{
	BT_CONTACT_POINT_ALGORITHMS = 1,
	BT_CLOSEST_POINT_ALGORITHMS = 2
};

///The btDispatcher interface class can be used in combination with broadphase to dispatch calculations for overlapping pairs.
///For example for pairwise collision detection, calculating contact points stored in btPersistentManifold or user callbacks (game logic).
class btDispatcher
{
public:
	typedef std::set<std::pair<const btCollisionObject*, const btCollisionObject*>> btInitialCollisionParticipants;
	typedef std::set<const btCollisionObject*> btInitialCollisionParticipantsSingle;
	typedef std::map<std::pair<int, int>, std::tuple<int, bool>> btPreviouslyFoundPairMap;  // value int is the time itself, bool is whether the value is up to date

	virtual ~btDispatcher();

	virtual btCollisionAlgorithm* findAlgorithm(const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, btPersistentManifold* sharedManifold, ebtDispatcherQueryType queryType) = 0;
	virtual btCollisionAlgorithm* findAlgorithm(const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, const int body0ShapeType, const int body1ShapeType, btPersistentManifold* sharedManifold, ebtDispatcherQueryType queryType) = 0;

	virtual btPersistentManifold* getNewManifold(const btCollisionObject* b0, const btCollisionObject* b1, size_t unlimitedSizeManifoldHint = 10) = 0;

	virtual void releaseManifold(btPersistentManifold* manifold) = 0;

	virtual void clearManifold(btPersistentManifold* manifold) = 0;

	virtual bool needsCollisionUsingAABBSimilarity(btBroadphaseProxy* proxy0, btBroadphaseProxy* proxy1) = 0;

	virtual bool needsCollision(const btCollisionObject* body0, const btCollisionObject* body1) = 0;

	virtual bool needsResponse(const btCollisionObject* body0, const btCollisionObject* body1) = 0;

	virtual void dispatchAllCollisionPairs(btOverlappingPairCache* pairCache, const btDispatcherInfo& dispatchInfo, btDispatcher* dispatcher) = 0;

	virtual int getNumManifolds() const = 0;

	virtual btPersistentManifold* getManifoldByIndexInternal(int index) = 0;

	virtual btPersistentManifold** getInternalManifoldPointer() = 0;

	virtual btPoolAllocator* getInternalManifoldPool() = 0;

	virtual const btPoolAllocator* getInternalManifoldPool() const = 0;

	virtual void* allocateCollisionAlgorithm(int size) = 0;

	virtual void freeCollisionAlgorithm(void* ptr) = 0;

	virtual void addInitialCollisionParticipant(btInitialCollisionParticipants::value_type pair) = 0;

	virtual const btInitialCollisionParticipants& getInitialCollisionParticipants() const = 0;

	virtual const btInitialCollisionParticipantsSingle& getInitialCollisionParticipants0() const = 0;

	virtual const btInitialCollisionParticipantsSingle& getInitialCollisionParticipants1() const = 0;

	virtual void addPreviouslyConsumedTime(btPreviouslyFoundPairMap::key_type key, btPreviouslyFoundPairMap::mapped_type value) = 0;

	virtual const btPreviouslyFoundPairMap& getPreviouslyConsumedTime() const = 0;
};

#endif  //BT_DISPATCHER_H
