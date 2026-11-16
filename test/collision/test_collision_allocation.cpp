#include <cstdint>
#include <tuple>
#include <cstring>
#include <new>
#include "gtest/gtest.h"
#include "BulletCollision/CollisionDispatch/btCollisionDispatcher.h"
#include "BulletCollision/CollisionDispatch/btDefaultCollisionConfiguration.h"
#include "BulletCollision/Gimpact/btGImpactCollisionAlgorithm.h"
#include "LinearMath/btPoolAllocator.h"

class InspectableGImpactAlgorithm : public btGImpactCollisionAlgorithm
{
public:
	explicit InspectableGImpactAlgorithm(const btCollisionAlgorithmConstructionInfo& ci)
		: btGImpactCollisionAlgorithm(ci, NULL, NULL)
	{
	}

	bool softAlgorithmIsNull() const
	{
		// Inspect poisoned storage without evaluating an indeterminate pointer.
		btCollisionAlgorithm* nullAlgorithm = NULL;
		return std::memcmp(&m_algorithmForSoftVsSoft, &nullAlgorithm, sizeof(nullAlgorithm)) == 0;
	}

	void resetSoftAlgorithm()
	{
		// Allow a failing regression to clean up without dereferencing poisoned memory.
		m_algorithmForSoftVsSoft = NULL;
	}
};

TEST(CollisionAllocation, GImpactInitializesSoftAlgorithmInReusedPoolStorage)
{
	btDefaultCollisionConstructionInfo configInfo;
	configInfo.m_customCollisionAlgorithmMaxElementSize = sizeof(InspectableGImpactAlgorithm);
	btDefaultCollisionConfiguration config(configInfo);
	btCollisionDispatcher dispatcher(&config);
	btCollisionAlgorithmConstructionInfo ci;
	ci.m_dispatcher1 = &dispatcher;

	for (int iteration = 0; iteration < 2; ++iteration)
	{
		void* storage = dispatcher.allocateCollisionAlgorithm(sizeof(InspectableGImpactAlgorithm));
		ASSERT_TRUE(storage != NULL);
		EXPECT_TRUE(config.getCollisionAlgorithmPool()->validPtr(storage));
		std::memset(storage, 0xcd, sizeof(InspectableGImpactAlgorithm));
		InspectableGImpactAlgorithm* algorithm = new (storage) InspectableGImpactAlgorithm(ci);
		const bool initialized = algorithm->softAlgorithmIsNull();
		if (!initialized)
			algorithm->resetSoftAlgorithm();
		algorithm->~InspectableGImpactAlgorithm();
		dispatcher.freeCollisionAlgorithm(storage);
		EXPECT_TRUE(initialized);
		EXPECT_EQ(0, config.getCollisionAlgorithmPool()->getUsedCount());
	}
}

TEST(CollisionAllocation, OversizedRequestLeavesPoolUnchanged)
{
	btPoolAllocator pool(32, 2);
	void* oversized = pool.allocate(33);
	EXPECT_TRUE(oversized == NULL);
	EXPECT_EQ(2, pool.getFreeCount());
	if (oversized)
		pool.freeMemory(oversized);
	void* first = pool.allocate(32);
	void* second = pool.allocate(32);
	EXPECT_TRUE(first == pool.getPoolAddress());
	EXPECT_TRUE(second == pool.getPoolAddress() + 32);
	EXPECT_TRUE(pool.allocate(32) == NULL);
	pool.freeMemory(first);
	pool.freeMemory(second);
	EXPECT_EQ(2, pool.getFreeCount());
}

TEST(CollisionAllocation, DispatcherFallsBackForOversizedAlgorithm)
{
	btPoolAllocator pool(16, 2);
	btDefaultCollisionConstructionInfo configInfo;
	configInfo.m_collisionAlgorithmPool = &pool;
	btDefaultCollisionConfiguration config(configInfo);
	btCollisionDispatcher dispatcher(&config);
	const int size = sizeof(btGImpactCollisionAlgorithm);
	ASSERT_GT(size, pool.getElementSize());
	void* memory = dispatcher.allocateCollisionAlgorithm(size);
	ASSERT_TRUE(memory != NULL);
	const bool fromPool = pool.validPtr(memory);
	EXPECT_FALSE(fromPool);
	EXPECT_EQ(2, pool.getFreeCount());
	if (!fromPool)
		std::memset(memory, 0xcd, size);
	dispatcher.freeCollisionAlgorithm(memory);
	EXPECT_EQ(2, pool.getFreeCount());
	void* pooledMemory = dispatcher.allocateCollisionAlgorithm(16);
	EXPECT_TRUE(pool.validPtr(pooledMemory));
	dispatcher.freeCollisionAlgorithm(pooledMemory);
	EXPECT_EQ(2, pool.getFreeCount());
}

int main(int argc, char** argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
