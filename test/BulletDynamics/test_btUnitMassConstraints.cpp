#include <btBulletDynamicsCommon.h>
#include <BulletDynamics/ConstraintSolver/btGeneric6DofSpring2Constraint.h>
#include <gtest/gtest.h>
#include <memory>
#include <vector>

namespace
{
struct Chain
{
    btDefaultCollisionConfiguration config;
    btCollisionDispatcher dispatcher;
    btDbvtBroadphase broadphase;
    btSequentialImpulseConstraintSolver solver;
    btDiscreteDynamicsWorld world;
    btBoxShape shape;
    btBoxShape floorShape;
    btRigidBody floor;
    std::vector<std::unique_ptr<btRigidBody> > bodies;
    std::vector<std::unique_ptr<btGeneric6DofSpring2Constraint> > joints;

    Chain(bool unitMasses, bool workaround, bool contacts)
        : dispatcher(&config), world(&dispatcher, &broadphase, &solver, &config), shape(btVector3(.2, .5, .15)),
          floorShape(btVector3(20, .2, 20)), floor(btRigidBody::btRigidBodyConstructionInfo(0, 0, &floorShape))
    {
        world.setGravity(btVector3(0, -10, 0));
        world.getSolverInfo().m_numIterations = 100;
        if (workaround) world.getSolverInfo().m_solverMode |= SOLVER_AS_IF_UNIT_MASS;
        // Powers of two make inertia normalization exact, avoiding divergent contact
        // point choices from roundoff while testing a 1024:1 mass ratio.
        const btScalar masses[] = {1, .03125, 32, .125};
        for (int i = 0; i < 4; ++i)
        {
            btScalar mass = unitMasses ? 1 : masses[i];
            btVector3 inertia;
            shape.calculateLocalInertia(mass, inertia);
            btRigidBody::btRigidBodyConstructionInfo info(mass, 0, &shape, inertia);
            bodies.emplace_back(new btRigidBody(info));
            bodies.back()->setRollingFriction(.02);
            bodies.back()->setSpinningFriction(.02);
            btTransform transform;
            transform.setIdentity();
            transform.setOrigin(btVector3(.1 * i, 4 - .9 * i, 0));
            bodies.back()->setWorldTransform(transform);
            world.addRigidBody(bodies.back().get(), 1, contacts ? -1 : 0);
            if (i)
            {
                btTransform a, b;
                a.setIdentity(); b.setIdentity();
                a.setOrigin(btVector3(0, -.5, 0));
                b.setOrigin(btVector3(0, .5, 0));
                joints.emplace_back(new btGeneric6DofSpring2Constraint(*bodies[i-1], *bodies[i], a, b));
                joints.back()->setLinearLowerLimit(btVector3(0, 0, 0));
                joints.back()->setLinearUpperLimit(btVector3(0, 0, 0));
                // Mix a fixed joint with ball joints, as in the ragdoll's right arm.
                joints.back()->setAngularLowerLimit(btVector3(i == 1 ? 0 : 1, i == 1 ? 0 : 1, i == 1 ? 0 : 1));
                joints.back()->setAngularUpperLimit(btVector3(i == 1 ? 0 : -1, i == 1 ? 0 : -1, i == 1 ? 0 : -1));
                world.addConstraint(joints.back().get(), true);
            }
        }
        bodies.back()->setLinearVelocity(btVector3(2, 0, 1));
        floor.setRollingFriction(.02);
        floor.setSpinningFriction(.02);
        if (contacts) world.addRigidBody(&floor);
    }
    ~Chain()
    {
        if (floor.isInWorld()) world.removeRigidBody(&floor);
        for (auto& joint : joints) world.removeConstraint(joint.get());
        for (auto& body : bodies) world.removeRigidBody(body.get());
    }
};
}

TEST(UnitMassConstraints, MixedFixedAndBallJointsMatchActualUnitMasses)
{
    Chain reference(true, true, false);
    Chain candidate(false, true, false);
    for (int step = 0; step < 300; ++step)
    {
        reference.world.stepSimulation(.002, 0);
        candidate.world.stepSimulation(.002, 0);
        for (int i = 0; i < 4; ++i)
        {
            EXPECT_LT((reference.bodies[i]->getWorldTransform().getOrigin() - candidate.bodies[i]->getWorldTransform().getOrigin()).length(), 1e-4);
            EXPECT_LT((reference.bodies[i]->getAngularVelocity() - candidate.bodies[i]->getAngularVelocity()).length(), 1e-4);
        }
    }
}

TEST(UnitMassConstraints, FreeHandSpringRemainsExcluded)
{
    Chain reference(false, false, false);
    Chain candidate(false, true, false);
    for (int i = 0; i < 3; ++i)
    {
        for (Chain* chain : {&reference, &candidate})
        {
            auto& joint = chain->joints[i];
            joint->setLinearLowerLimit(btVector3(1, 1, 1));
            joint->setLinearUpperLimit(btVector3(-1, -1, -1));
            joint->setAngularLowerLimit(btVector3(1, 1, 1));
            joint->setAngularUpperLimit(btVector3(-1, -1, -1));
            for (int axis = 0; axis < 3; ++axis)
            {
                joint->enableSpring(axis, true);
                joint->setStiffness(axis, 10);
                joint->setDamping(axis, 1);
            }
        }
    }
    for (int step = 0; step < 100; ++step)
    {
        reference.world.stepSimulation(.002, 0);
        candidate.world.stepSimulation(.002, 0);
        for (int i = 0; i < 4; ++i)
            EXPECT_LT((reference.bodies[i]->getWorldTransform().getOrigin() - candidate.bodies[i]->getWorldTransform().getOrigin()).length(), 1e-6);
    }
}

TEST(UnitMassConstraints, GroundContactAndFrictionMatchActualUnitMasses)
{
    Chain reference(true, true, true);
    Chain candidate(false, true, true);
    bool hitFloor = false;
    for (int step = 0; step < 600; ++step)
    {
        reference.world.stepSimulation(.002, 0);
        candidate.world.stepSimulation(.002, 0);
        hitFloor |= reference.dispatcher.getNumManifolds() > 0;
        for (int i = 0; i < 4; ++i)
        {
            EXPECT_LT((reference.bodies[i]->getWorldTransform().getOrigin() - candidate.bodies[i]->getWorldTransform().getOrigin()).length(), 1e-4);
            ASSERT_LT((reference.bodies[i]->getAngularVelocity() - candidate.bodies[i]->getAngularVelocity()).length(), 1e-3) << "step=" << step << " body=" << i;
        }
    }
    EXPECT_TRUE(hitFloor);
}

TEST(UnitMassConstraints, UnconstrainedContactsRemainPhysical)
{
    Chain reference(false, false, true);
    Chain candidate(false, true, true);
    for (Chain* chain : {&reference, &candidate})
    {
        for (auto& joint : chain->joints) chain->world.removeConstraint(joint.get());
        chain->joints.clear();
    }
    for (int step = 0; step < 400; ++step)
    {
        reference.world.stepSimulation(.002, 0);
        candidate.world.stepSimulation(.002, 0);
        for (int i = 0; i < 4; ++i)
            ASSERT_LT((reference.bodies[i]->getWorldTransform().getOrigin() - candidate.bodies[i]->getWorldTransform().getOrigin()).length(), 1e-6);
    }
}
TEST(UnitMassConstraints, InactiveLimitsStillUseConsistentContactMass)
{
    Chain reference(true, true, true);
    Chain candidate(false, true, true);
    for (Chain* chain : {&reference, &candidate})
    {
        for (auto& joint : chain->joints)
        {
            joint->setLinearLowerLimit(btVector3(-100, -100, -100));
            joint->setLinearUpperLimit(btVector3(100, 100, 100));
            joint->setAngularLowerLimit(btVector3(1, 1, 1));
            joint->setAngularUpperLimit(btVector3(-1, -1, -1));
        }
    }
    for (int step = 0; step < 400; ++step)
    {
        reference.world.stepSimulation(.002, 0);
        candidate.world.stepSimulation(.002, 0);
        for (int i = 0; i < 4; ++i)
            ASSERT_LT((reference.bodies[i]->getWorldTransform().getOrigin() - candidate.bodies[i]->getWorldTransform().getOrigin()).length(), 1e-4);
    }
}
int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
