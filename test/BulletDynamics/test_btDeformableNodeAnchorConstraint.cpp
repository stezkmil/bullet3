#include "BulletCollision/CollisionShapes/btSphereShape.h"
#include "BulletDynamics/Dynamics/btRigidBody.h"
#include "BulletSoftBody/btDeformableBodySolver.h"
#include "BulletSoftBody/btDeformableContactConstraint.h"
#include "BulletSoftBody/btDeformableContactProjection.h"
#include "gtest/gtest.h"

namespace
{
struct AnchorFixture
{
	btSphereShape shape;
	btRigidBody rigidBody;
	btSoftBody::Node node;
	btSoftBody::DeformableNodeRigidAnchor anchor;
	btContactSolverInfo solverInfo;

	AnchorFixture(btScalar softMass, btScalar rigidMass)
		: shape(1),
		  rigidBody(btRigidBody::btRigidBodyConstructionInfo(rigidMass, 0, &shape))
	{
		btVector3 localInertia;
		shape.calculateLocalInertia(rigidMass, localInertia);
		rigidBody.setMassProps(rigidMass, localInertia);
		rigidBody.updateInertiaTensor();

		node.m_im = btScalar(1) / softMass;
		node.m_v.setZero();
		node.m_vn.setZero();
		node.m_splitv.setZero();
		node.m_x.setZero();
		node.index = 0;

		anchor.m_node = &node;
		anchor.m_body = &rigidBody;
		anchor.m_userIndex = 0;
		anchor.m_cti.m_colObj = &rigidBody;
		anchor.m_cti.m_normal = btVector3(1, 0, 0);
		anchor.m_c1.setZero();
		anchor.m_c2 = node.m_im;
		anchor.m_local.setZero();

		solverInfo.m_sor = 1;
		solverInfo.m_timeStep = btScalar(1) / 100;
		solverInfo.m_deformable_erp = btScalar(0.2);
	}
};

void expectVectorNear(const btVector3& actual, const btVector3& expected, btScalar tolerance)
{
	EXPECT_NEAR(actual.x(), expected.x(), tolerance);
	EXPECT_NEAR(actual.y(), expected.y(), tolerance);
	EXPECT_NEAR(actual.z(), expected.z(), tolerance);
}
}  // namespace

GTEST_TEST(btDeformableNodeAnchorConstraint, VelocitySolveConservesMomentumAndMatchesVelocity)
{
	const btScalar softMass = 2;
	const btScalar rigidMass = 4;
	AnchorFixture fixture(softMass, rigidMass);
	fixture.node.m_v = btVector3(-3, 4, 1);
	fixture.rigidBody.setLinearVelocity(btVector3(1, -2, btScalar(0.5)));
	const btVector3 initialMomentum =
		fixture.node.m_v * softMass + fixture.rigidBody.getLinearVelocity() * rigidMass;

	btDeformableNodeAnchorConstraint constraint(fixture.anchor, fixture.solverInfo);
	EXPECT_GT(constraint.solveConstraint(fixture.solverInfo), 0);

	expectVectorNear(fixture.node.m_v, fixture.rigidBody.getLinearVelocity(), btScalar(1e-5));
	const btVector3 finalMomentum =
		fixture.node.m_v * softMass + fixture.rigidBody.getLinearVelocity() * rigidMass;
	expectVectorNear(finalMomentum, initialMomentum, btScalar(1e-5));
}

GTEST_TEST(btDeformableNodeAnchorConstraint, VelocitySolveIncludesLeverArmAndRigidFactors)
{
	AnchorFixture fixture(2, 4);
	fixture.anchor.m_c1 = btVector3(btScalar(0.75), btScalar(-0.4), btScalar(0.2));
	fixture.anchor.m_local = fixture.anchor.m_c1;
	fixture.node.m_x = fixture.anchor.m_c1;
	fixture.node.m_v = btVector3(-3, 4, 1);
	fixture.rigidBody.setLinearVelocity(btVector3(1, -2, btScalar(0.5)));
	fixture.rigidBody.setAngularVelocity(btVector3(btScalar(0.3), btScalar(0.7), btScalar(-0.6)));
	fixture.rigidBody.setLinearFactor(btVector3(btScalar(0.5), 1, btScalar(0.25)));
	fixture.rigidBody.setAngularFactor(btVector3(1, btScalar(0.3), btScalar(0.8)));

	btDeformableNodeAnchorConstraint constraint(fixture.anchor, fixture.solverInfo);
	EXPECT_GT(constraint.solveConstraint(fixture.solverInfo), 0);
	expectVectorNear(constraint.getVb(), constraint.getVa(), btScalar(1e-5));
}

GTEST_TEST(btDeformableNodeAnchorConstraint, SplitSolveUsesOnlyPositionAndPushVelocity)
{
	const btScalar softMass = 2;
	const btScalar rigidMass = 4;
	AnchorFixture fixture(softMass, rigidMass);
	const btVector3 positionError(btScalar(0.1), btScalar(-0.2), btScalar(0.05));
	fixture.node.m_x = positionError;
	fixture.node.m_v = btVector3(7, -8, 9);
	fixture.node.m_vn = btVector3(-100, 200, -300);
	fixture.rigidBody.setLinearVelocity(btVector3(-4, 5, -6));

	btDeformableNodeAnchorConstraint constraint(fixture.anchor, fixture.solverInfo);
	EXPECT_GT(constraint.solveSplitImpulse(fixture.solverInfo), 0);

	const btVector3 relativePushVelocity =
		fixture.node.m_splitv - fixture.rigidBody.getPushVelocity();
	const btVector3 expectedRelativePushVelocity =
		-positionError * (fixture.solverInfo.m_deformable_erp / fixture.solverInfo.m_timeStep);
	expectVectorNear(relativePushVelocity, expectedRelativePushVelocity, btScalar(1e-5));

	const btVector3 splitMomentum =
		fixture.node.m_splitv * softMass + fixture.rigidBody.getPushVelocity() * rigidMass;
	expectVectorNear(splitMomentum, btVector3(0, 0, 0), btScalar(1e-5));
	expectVectorNear(fixture.node.m_v, btVector3(7, -8, 9), btScalar(0));
}

GTEST_TEST(btDeformableNodeAnchorConstraint, ComplianceRegularizesVelocityConstraint)
{
	AnchorFixture fixture(1, 1);
	fixture.node.m_v = btVector3(2, 0, 0);
	fixture.node.m_x = btVector3(btScalar(0.1), 0, 0);
	fixture.anchor.m_compliance =
		btScalar(2) * fixture.solverInfo.m_timeStep * fixture.solverInfo.m_timeStep;

	btDeformableNodeAnchorConstraint constraint(fixture.anchor, fixture.solverInfo);
	EXPECT_GT(constraint.solveConstraint(fixture.solverInfo), 0);

	const btScalar regularization = fixture.anchor.m_compliance /
		(fixture.solverInfo.m_timeStep * fixture.solverInfo.m_timeStep);
	const btVector3 relativeVelocity = constraint.getVb() - constraint.getVa();
	const btVector3 residual = relativeVelocity +
		fixture.node.m_x / fixture.solverInfo.m_timeStep -
		constraint.m_totalImpulse * regularization;
	expectVectorNear(residual, btVector3(0, 0, 0), btScalar(1e-5));
	EXPECT_GT(relativeVelocity.length2(), btScalar(0));
	EXPECT_NEAR(constraint.solveConstraint(fixture.solverInfo), 0, btScalar(1e-5));

	const btVector3 momentum = fixture.node.m_v + fixture.rigidBody.getLinearVelocity();
	expectVectorNear(momentum, btVector3(2, 0, 0), btScalar(1e-5));
}

GTEST_TEST(btDeformableNodeAnchorConstraint, CompliantAnchorDoesNotDuplicatePositionCorrectionInSplitSolve)
{
	AnchorFixture fixture(1, 1);
	const btVector3 positionError(btScalar(0.1), 0, 0);
	fixture.node.m_x = positionError;
	fixture.anchor.m_compliance =
		btScalar(2) * fixture.solverInfo.m_timeStep * fixture.solverInfo.m_timeStep;

	btDeformableNodeAnchorConstraint constraint(fixture.anchor, fixture.solverInfo);
	EXPECT_EQ(constraint.solveSplitImpulse(fixture.solverInfo), 0);
	expectVectorNear(constraint.getSplitVb() - constraint.getSplitVa(), btVector3(0, 0, 0), btScalar(0));
	expectVectorNear(constraint.m_totalSplitImpulse, btVector3(0, 0, 0), btScalar(0));
}

GTEST_TEST(btDeformableNodeAnchorConstraint, CompliantAnchorSupportsStaticRigidBody)
{
	AnchorFixture fixture(1, 0);
	fixture.node.m_v = btVector3(2, 0, 0);
	fixture.anchor.m_compliance =
		fixture.solverInfo.m_timeStep * fixture.solverInfo.m_timeStep;

	btDeformableNodeAnchorConstraint constraint(fixture.anchor, fixture.solverInfo);
	EXPECT_GT(constraint.solveConstraint(fixture.solverInfo), 0);
	expectVectorNear(fixture.node.m_v, btVector3(1, 0, 0), btScalar(1e-5));
	expectVectorNear(fixture.rigidBody.getLinearVelocity(), btVector3(0, 0, 0), btScalar(0));
}

GTEST_TEST(btDeformableNodeAnchorConstraint, AnchorIgnoresCollisionResponseFlag)
{
	AnchorFixture fixture(1, 0);
	fixture.rigidBody.setCollisionFlags(
		fixture.rigidBody.getCollisionFlags() | btCollisionObject::CF_NO_CONTACT_RESPONSE);
	fixture.node.m_x = btVector3(btScalar(0.1), 0, 0);
	fixture.anchor.m_compliance =
		fixture.solverInfo.m_timeStep * fixture.solverInfo.m_timeStep;

	btDeformableNodeAnchorConstraint constraint(fixture.anchor, fixture.solverInfo);
	EXPECT_GT(constraint.solveConstraint(fixture.solverInfo), 0);
	EXPECT_LT(fixture.node.m_v.x(), 0);
}

GTEST_TEST(btDeformableNodeAnchorConstraint, CompliantAnchorVelocityIsPreservedByProjection)
{
	AnchorFixture fixture(1, 1);
	fixture.anchor.m_compliance = btScalar(1e-5);
	btSoftBodyWorldInfo worldInfo;
	btSoftBody softBody(&worldInfo);

	btAlignedObjectArray<btSoftBody*> softBodies;
	softBodies.push_back(&softBody);
	btDeformableContactProjection projection(softBodies);
	projection.reinitialize(true);
	projection.m_nodeAnchorConstraints[0].push_back(
		btDeformableNodeAnchorConstraint(fixture.anchor, fixture.solverInfo));
	projection.setProjection();

	btDeformableContactProjection::TVStack increment;
	increment.resize(1, btVector3(1, 2, 3));
	projection.project(increment);
	expectVectorNear(increment[0], btVector3(0, 0, 0), btScalar(0));
}

GTEST_TEST(btDeformableNodeAnchorConstraint, ZeroMassNodeIsNotIntegrated)
{
	btSoftBodyWorldInfo worldInfo;
	worldInfo.m_maxDisplacement = 100;
	btVector3 position(4, 5, 6);
	btScalar mass = 1;
	btSoftBody softBody(&worldInfo, 1, &position, &mass);
	softBody.m_nodes[0].m_v = btVector3(10, 20, 30);
	softBody.m_nodes[0].m_vn = btVector3(40, 50, 60);
	softBody.m_nodes[0].m_splitv = btVector3(70, 80, 90);
	softBody.setMass(0, 0);

	btAlignedObjectArray<btSoftBody*> softBodies;
	softBodies.push_back(&softBody);
	btDeformableBodySolver solver;
	solver.reinitialize(softBodies, btScalar(0.01));
	solver.applyTransforms(btScalar(0.01));

	expectVectorNear(softBody.m_nodes[0].m_x, position, btScalar(0));
	expectVectorNear(softBody.m_nodes[0].m_v, btVector3(0, 0, 0), btScalar(0));
	expectVectorNear(softBody.m_nodes[0].m_splitv, btVector3(0, 0, 0), btScalar(0));
}

int main(int argc, char** argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
