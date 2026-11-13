#include "BulletSoftBody/btDeformableContactConstraint.h"
#include "BulletCollision/CollisionShapes/btBoxShape.h"
#include <gtest/gtest.h>
#include "BulletSoftBody/btDeformableMultiBodyDynamicsWorld.h"
#include "BulletSoftBody/btDeformableMultiBodyConstraintSolver.h"
#include "BulletSoftBody/btDeformableBodySolver.h"
#include "BulletSoftBody/btDeformableLinearElasticityForce.h"
#include "BulletSoftBody/btSoftBodyRigidBodyCollisionConfiguration.h"
#include "BulletCollision/BroadphaseCollision/btDbvtBroadphase.h"
#include "BulletCollision/CollisionDispatch/btCollisionDispatcher.h"
#include <cmath>

class DeformableAnchor : public ::testing::Test
{
protected:
    btSoftBodyWorldInfo worldInfo;
    btBoxShape shape{btVector3(1, 2, 3)};
    btRigidBody rigid{btScalar(2), nullptr, &shape, btVector3(3, 5, 7)};
    btSoftBody::Node node{};
    btSoftBody::DeformableNodeRigidAnchor anchor{};
    btContactSolverInfo info;

    void SetUp() override
    {
        info.m_timeStep = btScalar(0.002);
        node.m_im = btScalar(0.25);
        node.m_x = btVector3(4, 2, 1);
        node.m_v = btVector3(3, -2, 1);
        node.m_vn = btVector3(-100, 200, 50);
        node.m_splitv.setZero();
        anchor.m_body = &rigid;
        anchor.m_cti.m_colObj = &rigid;
        anchor.m_cti.m_normal = btVector3(1, 0, 0);
        anchor.m_node = &node;
        anchor.m_c2 = node.m_im;
        anchor.m_c1 = node.m_x;
        anchor.m_local = node.m_x;
        // Deliberately stale: the solve must derive its response from the current inertia.
        anchor.m_c0 = btMatrix3x3::getIdentity();
        rigid.setLinearVelocity(btVector3(-2, 1, 3));
        rigid.setAngularVelocity(btVector3(1, 2, -1));
    }

    btVector3 momentum() const { return node.m_v / node.m_im + rigid.getLinearVelocity() / rigid.getInvMass(); }
    btVector3 angularMomentum() const
    {
        return node.m_x.cross(node.m_v / node.m_im) + rigid.getInvInertiaTensorWorld().inverse() * rigid.getAngularVelocity();
    }
    btScalar energy() const
    {
        return btScalar(0.5) * (node.m_v.length2() / node.m_im + rigid.getLinearVelocity().length2() / rigid.getInvMass()
            + rigid.getAngularVelocity().dot(rigid.getInvInertiaTensorWorld().inverse() * rigid.getAngularVelocity()));
    }
};

TEST_F(DeformableAnchor, RotatedOffCentreImpulseConservesMomentumAndDissipatesEnergy)
{
    rigid.setCenterOfMassTransform(btTransform(btQuaternion(btVector3(1, 2, 3).normalized(), btScalar(1.2)), btVector3(0, 0, 0)));
    anchor.m_local = rigid.getWorldTransform().invXform(node.m_x);
    const auto linear = momentum();
    const auto angular = angularMomentum();
    const auto kinetic = energy();
    btDeformableNodeAnchorConstraint constraint(anchor, info);
    constraint.solveConstraint(info);
    EXPECT_LT((constraint.getVa() - constraint.getVb()).length(), 1e-10);
    EXPECT_LT((momentum() - linear).length(), 1e-10);
    EXPECT_LT((angularMomentum() - angular).length(), 1e-10);
    EXPECT_LE(energy(), kinetic);
}

TEST_F(DeformableAnchor, SplitCorrectionUsesMassAndDoesNotChangePhysicalVelocity)
{
    for (const btScalar mass : {btScalar(0.001), btScalar(1), btScalar(1000)})
    {
        rigid.setMassProps(mass, btVector3(3, 5, 7));
        rigid.updateInertiaTensor();
        anchor.m_c1.setZero();
        anchor.m_local.setZero();
        node.m_x = btVector3(0.1, -0.2, 0.3);
        node.m_splitv.setZero();
        rigid.setPushVelocity(btVector3(0, 0, 0));
        rigid.setTurnVelocity(btVector3(0, 0, 0));
        const auto before = momentum();
        btDeformableNodeAnchorConstraint constraint(anchor, info);
        constraint.solveSplitImpulse(info);
        EXPECT_LT(constraint.solveSplitImpulse(info), 1e-18);
        EXPECT_LT((node.m_splitv / node.m_im + rigid.getPushVelocity() / rigid.getInvMass()).length(), 1e-10);
        EXPECT_LT((momentum() - before).length(), 1e-10);
        EXPECT_LT((node.m_x + info.m_timeStep * (node.m_splitv - rigid.getPushVelocity())
            - node.m_x * (1 - info.m_erp)).length(), 1e-10);
    }
}

TEST_F(DeformableAnchor, SplitCorrectionIgnoresPreviousPhysicalVelocity)
{
    node.m_x = rigid.getWorldTransform() * anchor.m_local;
    btDeformableNodeAnchorConstraint constraint(anchor, info);
    EXPECT_EQ(constraint.solveSplitImpulse(info), 0);
    EXPECT_EQ(node.m_splitv.length2(), 0);
    EXPECT_EQ(rigid.getPushVelocity().length2(), 0);
}

class DeformableAnchorWorld : public ::testing::Test
{
protected:
    void runDampedCable(btScalar rigidMass, bool multipleAnchors, bool groundContact = false, bool implicit = false, bool lineSearch = false, bool projection = false, bool rigidImpulse = false)
    {
        btSoftBodyRigidBodyCollisionConfiguration configuration;
        btCollisionDispatcher dispatcher(&configuration);
        btDbvtBroadphase broadphase;
        btDeformableMultiBodyConstraintSolver constraints;
        btDeformableBodySolver solver;
        btDeformableMultiBodyDynamicsWorld world(&dispatcher, &broadphase, &constraints, &configuration, &solver);
        constraints.setDeformableSolver(&solver);
        world.setImplicit(implicit);
        solver.setMaxNewtonIterations(10);
        solver.setNewtonTolerance(1e-6);
        solver.setLineSearch(lineSearch);
        world.setUseProjection(projection);
        world.setGravity(btVector3(0, 0, 0));
        world.getSolverInfo().m_numIterations = 100;
        world.getSolverInfo().m_leastSquaresResidualThreshold = 1e-16;
        const btVector3 positions[] = {btVector3(0, 0, 0), btVector3(1, 0, 0), btVector3(0, 1, 0), btVector3(0, 0, 1)};
        const btScalar masses[] = {1, 1, 1, 1};
        btSoftBody soft(&world.getWorldInfo(), 4, positions, masses);
        soft.appendTetra(0, 1, 2, 3);
        soft.initializeDmInverse();
        soft.m_tetraScratches.resize(1);
        soft.m_tetraScratchesTn.resize(1);
        soft.m_cfg.collisions = 0;
        soft.setActivationState(DISABLE_DEACTIVATION);
        world.addSoftBody(&soft, 1, 0);
        btDeformableLinearElasticityForce material(10000, 10000, 0, 1);
        world.addForce(&soft, &material);
        btBoxShape shape(btVector3(0.1, 0.1, 0.1));
        btVector3 inertia;
        shape.calculateLocalInertia(rigidMass, inertia);
        btRigidBody rigid(rigidMass, nullptr, &shape, inertia);
        rigid.setActivationState(DISABLE_DEACTIVATION);
        world.addRigidBody(&rigid, 1, groundContact ? 1 : 0);
        btBoxShape groundShape(btVector3(2, 2, 0.1));
        btRigidBody ground(0, nullptr, &groundShape);
        ground.setWorldTransform(btTransform(btQuaternion::getIdentity(), btVector3(0, 0, -0.15)));
        if (groundContact) world.addRigidBody(&ground, 1, 1);
        soft.appendDeformableAnchor(0, &rigid, 0);
        if (multipleAnchors) soft.appendDeformableAnchor(1, &rigid, 0);
        soft.m_nodes[3].m_v = btVector3(0, 0, (groundContact || rigidImpulse) ? 0 : 1);
        soft.m_nodes[3].m_vn = soft.m_nodes[3].m_v;
        if (rigidImpulse) rigid.applyCentralImpulse(btVector3(1, 0, 0));
        for (int step = 0; step < 200; ++step)
        {
            world.stepSimulation(btScalar(0.002), 0);
            btVector3 momentum = rigid.getLinearVelocity() * rigidMass;
            for (int n = 0; n < 4; ++n)
            {
                momentum += soft.m_nodes[n].m_v;
                EXPECT_TRUE(std::isfinite(double(soft.m_nodes[n].m_x.length2())));
                EXPECT_LT(soft.m_nodes[n].m_v.length(), groundContact ? 5 : 2);
                if (groundContact && step == 0) EXPECT_LT(soft.m_nodes[n].m_v.length(), 1e-8);
            }
            if (!groundContact) EXPECT_LT((momentum - (rigidImpulse ? btVector3(1, 0, 0) : btVector3(0, 0, 1))).length(), 2e-4);
            if (rigidImpulse && step == 0)
            {
                btAlignedObjectArray<btVector3> damping;
                damping.resize(4, btVector3(0, 0, 0));
                material.addScaledDampingForce(btScalar(0.002), damping);
                for (int n = 1; n < 4; ++n)
                    EXPECT_LT((soft.m_nodes[n].m_v - damping[n]).length(), 2e-6);
            }
            EXPECT_LT((soft.m_nodes[0].m_x - rigid.getWorldTransform().getOrigin()).length(), 1e-4);
        }
        if (multipleAnchors) soft.removeAnchor(1);
        soft.removeAnchor(0);
        if (groundContact) world.removeRigidBody(&ground);
        world.removeRigidBody(&rigid);
        world.removeForce(&material);
        world.removeSoftBody(&soft);
    }

};

TEST_F(DeformableAnchorWorld, DampedCableAndLightPlugConserveLinearMomentum)
{
    runDampedCable(btScalar(0.01), false);
}

TEST_F(DeformableAnchorWorld, MultipleAnchorsShareRigidInertia)
{
    runDampedCable(btScalar(0.01), true);
}

TEST_F(DeformableAnchorWorld, RigidImpulseIsIncludedInCableDampingSolve)
{
    runDampedCable(btScalar(0.01), false, false, false, false, false, true);
}

TEST_F(DeformableAnchorWorld, DampedCableAndHeavyPlugConserveLinearMomentum)
{
    runDampedCable(btScalar(1000), true);
}

TEST_F(DeformableAnchorWorld, ContactAndAnchorShareSplitCorrection)
{
    runDampedCable(btScalar(1), false, true);
}

TEST_F(DeformableAnchorWorld, ImplicitLightPlugConservesMomentum)
{
    runDampedCable(btScalar(0.01), true, false, true);
}

TEST_F(DeformableAnchorWorld, ImplicitHeavyPlugConservesMomentum)
{
    runDampedCable(btScalar(1000), true, false, true);
}

TEST_F(DeformableAnchorWorld, ImplicitLineSearchConservesMomentum)
{
    runDampedCable(btScalar(0.01), true, false, true, true);
}

TEST_F(DeformableAnchorWorld, ImplicitProjectionConservesMomentum)
{
    runDampedCable(btScalar(0.01), true, false, true, false, true);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
