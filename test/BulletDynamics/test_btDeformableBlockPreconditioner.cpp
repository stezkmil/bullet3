#include "BulletSoftBody/btDeformableBodySolver.h"
#include "LinearMath/btQuaternion.h"
#include "deformable_solver_test_helpers.h"
#include <gtest/gtest.h>
#include <limits>

namespace
{
typedef btAlignedObjectArray<btVector3> Vectors;

class DeformableBlockPreconditioner : public ::testing::Test
{
protected:
	btSoftBodyWorldInfo info;
	btSoftBody* body;
	btDeformableLinearElasticityForce force;
	btAlignedObjectArray<btSoftBody*> bodies;
	btAlignedObjectArray<btDeformableLagrangianForce*> forces;

	DeformableBlockPreconditioner() : body(0), force(1200, 2400, btScalar(0.12), btScalar(0.04)) {}
	virtual void SetUp()
	{
		const btVector3 positions[] = {btVector3(0, 0, 0), btVector3(2, 0, 0), btVector3(0, 1, 0), btVector3(0, 0, 3)};
		const btScalar masses[] = {1, 2, 3, 4};
		body = new btSoftBody(&info, 4, positions, masses);
		body->appendTetra(0, 1, 2, 3);
		body->initializeDmInverse();
		body->m_tetraScratches.resize(1);
		body->m_tetraScratches[0].m_corotation = btMatrix3x3(btQuaternion(btVector3(1, 2, 3).normalized(), btScalar(0.73)));
		body->m_tetraScratches[0].m_J = 1;
		for (int i = 0; i < 4; ++i) body->m_nodes[i].index = i;
		force.addSoftBody(body);
		bodies.push_back(body);
		forces.push_back(&force);
	}
	virtual void TearDown() { delete body; }

	void compareBlocksToOperator(btScalar dt, bool flat)
	{
		body->m_tetraScratches[0].m_J = flat ? btScalar(0.001) : btScalar(1);
		btAlignedObjectArray<btMatrix3x3> blocks;
		blocks.resize(4);
		for (int n = 0; n < 4; ++n) blocks[n] = btMatrix3x3::getIdentity() * btScalar(0);
		ASSERT_TRUE(force.addImplicitForceDifferentialBlocks(dt, blocks));
		Vectors direction, product;
		direction.resize(4);
		product.resize(4);
		for (int n = 0; n < 4; ++n)
			for (int axis = 0; axis < 3; ++axis)
			{
				for (int j = 0; j < 4; ++j) { direction[j].setZero(); product[j].setZero(); }
				direction[n][axis] = 1;
				force.addScaledDampingForceDifferential(-dt, direction, product);
				force.addScaledElasticForceDifferential(-dt * dt, direction, product);
				for (int r = 0; r < 3; ++r)
					EXPECT_NEAR((double)product[n][r], (double)blocks[n][r][axis],
						2e-5 * btMax(btScalar(1), btFabs(product[n][r])));
			}
	}
};

TEST_F(DeformableBlockPreconditioner, MatchesForceDifferentialForRotatedAndFlatElements)
{
	compareBlocksToOperator(btScalar(0.0002), false);
	compareBlocksToOperator(btScalar(0.02), false);
	compareBlocksToOperator(btScalar(0.02), true);
	force.setDamping(0, 0);
	compareBlocksToOperator(btScalar(0.02), true);
	compareBlocksToOperator(0, false);
}


TEST_F(DeformableBlockPreconditioner, CombinedDifferentialMatchesSeparatePasses)
{
	Vectors x, reference, combined;
	x.resize(4); reference.resize(4); combined.resize(4);
	const btScalar timesteps[] = {0, btScalar(0.0002), btScalar(0.02)};
	for (int flat = 0; flat < 2; ++flat)
	for (int damping = 0; damping < 4; ++damping)
	for (int fixed = 0; fixed < 2; ++fixed)
	for (int inactive = 0; inactive < 2; ++inactive)
	for (int t = 0; t < 3; ++t)
	for (int basis = 0; basis < 12; ++basis)
	{
		body->m_tetraScratches[0].m_J = flat ? btScalar(0.001) : btScalar(1);
		force.setDamping(damping & 1 ? btScalar(0.12) : 0, damping & 2 ? btScalar(0.04) : 0);
		body->m_nodes[0].m_frozen = fixed;
		body->m_nodes[1].m_im = fixed ? btScalar(0) : btScalar(0.5);
		body->forceActivationState(inactive ? ISLAND_SLEEPING : ACTIVE_TAG);
		for (int n = 0; n < 4; ++n)
		{
			x[n].setZero();
			reference[n] = combined[n] = btVector3(btScalar(0.1), btScalar(-0.2), btScalar(0.3));
		}
		x[basis / 3][basis % 3] = 1;
		const btScalar dt = timesteps[t];
		force.addScaledDampingForceDifferential(-dt, x, reference);
		force.addScaledElasticForceDifferential(-dt * dt, x, reference);
		force.addImplicitForceDifferential(dt, x, combined);
		Vectors cached;
		cached.resize(4, btVector3(btScalar(0.1), btScalar(-0.2), btScalar(0.3)));
		force.prepareImplicitForceDifferential(dt);
		force.addImplicitForceDifferential(dt, x, cached);
		force.finishImplicitForceDifferential();
		for (int n = 0; n < 4; ++n)
		for (int d = 0; d < 3; ++d)
			EXPECT_NEAR(double(reference[n][d]), double(cached[n][d]),
				256 * SIMD_EPSILON * btMax(btScalar(1), btFabs(reference[n][d])));
		for (int n = 0; n < 4; ++n)
		for (int d = 0; d < 3; ++d)
			EXPECT_NEAR(double(reference[n][d]), double(combined[n][d]),
				128 * SIMD_EPSILON * btMax(btScalar(1), btFabs(reference[n][d])));
	}
}


TEST_F(DeformableBlockPreconditioner, CachedOperatorRefreshAndFallback)
{
	Vectors x, expected, actual;
	x.resize(4); expected.resize(4); actual.resize(4);
	for (int n = 0; n < 4; ++n) x[n] = btVector3(n - 2, n * n, 3 - n);
	for (int state = 0; state < 4; ++state)
	{
		body->m_tetraScratches[0].m_corotation = btMatrix3x3(btQuaternion(btVector3(2,1,3).normalized(), btScalar(0.3 * state)));
		if (state == 3) body->m_tetraScratches[0].m_corotation[0][0] += btScalar(0.1);
		force.setLameParameters(1000 + state * 700, 2000 + state * 900);
		const btScalar dt = btScalar(0.003) * (state + 1);
		for (int n = 0; n < 4; ++n) { expected[n].setZero(); actual[n].setZero(); }
		force.addImplicitForceDifferential(dt, x, expected);
		force.prepareImplicitForceDifferential(dt);
		ASSERT_EQ(1, force.m_implicitTetraCache.size());
		EXPECT_EQ(state != 3, force.m_implicitTetraCache[0].usable);
		force.addImplicitForceDifferential(dt, x, actual);
		force.finishImplicitForceDifferential();
		EXPECT_FALSE(force.m_implicitCacheReady);
		for (int n = 0; n < 4; ++n)
		for (int d = 0; d < 3; ++d)
			EXPECT_NEAR(double(expected[n][d]), double(actual[n][d]),
				256 * SIMD_EPSILON * btMax(btScalar(1), btFabs(expected[n][d])));
	}
}

TEST_F(DeformableBlockPreconditioner, CachedNewtonSolveClosesCache)
{
	btDeformableBodySolver solver;
	solver.setImplicit(true);
	solver.setMaxNewtonIterations(3);
	solver.m_objective->m_lf.push_back(&force);
	const btScalar dt = btScalar(0.0002);
	solver.reinitialize(bodies, dt);
	for (int n = 0; n < 4; ++n)
	{
		body->m_nodes[n].m_vn.setZero();
		body->m_nodes[n].m_v = btVector3(btScalar(0.1 * n), 0, 0);
	}
	solver.setupDeformableSolve(true);
	solver.solveDeformableConstraints(dt);
	EXPECT_FALSE(force.m_implicitCacheReady);

}

TEST_F(DeformableBlockPreconditioner, InvertsNodeBlocksAndScalesContactSchurDiagonal)
{
	btScalar dt = btScalar(0.02);
	bool implicit = true;
	btDeformableContactProjection projection(bodies);
	LagrangeMultiplier lm = {};
	lm.m_num_nodes = 2;
	lm.m_num_constraints = 1;
	lm.m_indices[0] = 0; lm.m_indices[1] = 2;
	lm.m_weights[0] = btScalar(0.3); lm.m_weights[1] = btScalar(0.7);
	lm.m_dirs[0] = btVector3(1, 2, -1).normalized();
	projection.m_lagrangeMultipliers.push_back(lm);
	KKTPreconditioner preconditioner(bodies, projection, forces, dt, implicit);
	preconditioner.reinitialize(true);
	btAlignedObjectArray<btMatrix3x3> blocks;
	blocks.resize(4);
	for (int n = 0; n < 4; ++n) blocks[n] = btMatrix3x3::getIdentity() * (1 / body->m_nodes[n].m_im);
	force.addImplicitForceDifferentialBlocks(dt, blocks);
	for (int n = 0; n < 4; ++n)
		for (int axis = 0; axis < 3; ++axis)
		{
			btVector3 unit(0, 0, 0); unit[axis] = 1;
			const btVector3 result = preconditioner.applyInverseNodeBlock(n, blocks[n] * unit);
			for (int d = 0; d < 3; ++d) EXPECT_NEAR((double)unit[d], (double)result[d], 2e-5);
		}
	btScalar schur = 0;
	for (int n = 0; n < 2; ++n)
		schur += lm.m_weights[n] * lm.m_weights[n] * lm.m_dirs[0].dot(blocks[lm.m_indices[n]].inverse() * lm.m_dirs[0]);
	Vectors input, output;
	input.resize(5); output.resize(5);
	for (int n = 0; n < 5; ++n) input[n].setZero();
	input[4][0] = 1;
	preconditioner(input, output);
	EXPECT_NEAR(1.0, (double)(output[4][0] * schur), 2e-5);
	EXPECT_EQ(btScalar(0), output[4][1]);
	EXPECT_EQ(btScalar(0), output[4][2]);

	// Switching back to explicit must recover the previous mass-only approximation.
	implicit = false;
	preconditioner.reinitialize(false);
	const btVector3 v(1, -2, 3);
	for (int n = 0; n < 4; ++n)
	{
		const btVector3 actual = preconditioner.applyInverseNodeBlock(n, v);
		const btVector3 expected = v * body->m_nodes[n].m_im;
		for (int d = 0; d < 3; ++d) EXPECT_EQ(expected[d], actual[d]);
	}
}

TEST_F(DeformableBlockPreconditioner, FrozenAndZeroMassNodesRemainZero)
{
	body->m_nodes[0].m_frozen = 1;
	body->m_nodes[1].m_im = 0;
	btScalar dt = btScalar(0.02);
	bool implicit = true;
	btDeformableContactProjection projection(bodies);
	KKTPreconditioner preconditioner(bodies, projection, forces, dt, implicit);
	preconditioner.reinitialize(true);
	for (int n = 0; n < 2; ++n)
		EXPECT_EQ(btScalar(0), preconditioner.applyInverseNodeBlock(n, btVector3(1, 2, 3)).length2());
}

TEST(DeformableBlockInverse, RegularizesSemidefiniteAndRejectsNonfiniteBlocks)
{
	btMatrix3x3 block(1, 0, 0, 0, 1, 0, 0, 0, 0), inverse;
	bool regularized = false;
	ASSERT_TRUE(KKTPreconditioner::invertPositiveBlock(block, inverse, regularized));
	EXPECT_TRUE(regularized);
	for (int d = 0; d < 3; ++d) EXPECT_GT(inverse[d][d], btScalar(0));
	block[0][0] = std::numeric_limits<btScalar>::quiet_NaN();
	EXPECT_FALSE(KKTPreconditioner::invertPositiveBlock(block, inverse, regularized));
}

TEST_F(DeformableBlockPreconditioner, ExplicitAndNewtonConstraintRhsHaveOppositeSigns)
{
	Vectors backup, residual, dv, rhs;
	backup.resize(4); residual.resize(4); dv.resize(4);
	for (int n = 0; n < 4; ++n) { backup[n].setZero(); residual[n].setZero(); dv[n] = btVector3(2, 3, 4); }
	btDeformableBackwardEulerObjective objective(bodies, backup);
	LagrangeMultiplier lm = {};
	lm.m_num_nodes = 1; lm.m_num_constraints = 1;
	lm.m_indices[0] = 0; lm.m_weights[0] = 1; lm.m_dirs[0] = btVector3(0, 1, 0);
	objective.m_projection.m_lagrangeMultipliers.push_back(lm);
	objective.setImplicit(false);
	objective.addLagrangeMultiplierRHS(residual, dv, rhs);
	EXPECT_EQ(btScalar(3), rhs[4][0]);
	objective.setImplicit(true);
	objective.addLagrangeMultiplierRHS(residual, dv, rhs);
	EXPECT_EQ(btScalar(-3), rhs[4][0]);
}
// Probe the residual checks without advancing or running a simulation.
class ImplicitResidualProbe : public btDeformableBodySolver
{
public:
	void configure(btSoftBody* body)
	{
		m_softBodies.push_back(body);
		m_numNodes = body->m_nodes.size();
		setImplicit(true);
		m_objective->setDt(0);
		m_objective->updateId();
	}
	void prepareContactStep()
	{
		m_dv.resize(m_numNodes, btVector3(0,0,0)); m_backupVelocity.resize(m_numNodes, btVector3(0,0,0));
		setupDeformableSolve(true);
	}
	btScalar constraintError(const Vectors& dv) const { return implicitConstraintError(dv); }
	void measure(const Vectors& x, const Vectors& rhs) { measureImplicitLinearResidual(x, rhs); }
	btScalar linearResidual() const { return m_lastLinearMomentumResidual; }
	btScalar stationarityResidual() const { return m_lastStationarityResidual; }
};

TEST_F(DeformableBlockPreconditioner, SmallAccurateCorrectionDoesNotImplyCurrentStateConverged)
{
	body->m_nodes[0].m_im = btScalar(1e-9);
	ImplicitResidualProbe solver;
	solver.configure(body);
	Vectors x, rhs;
	x.resize(4); rhs.resize(4);
	for (int n = 0; n < 4; ++n) { x[n].setZero(); rhs[n].setZero(); }
	x[0][0] = btScalar(1e-6);
	rhs[0] = x[0] / body->m_nodes[0].m_im;
	solver.measure(x, rhs);
	EXPECT_NEAR(0.0, (double)solver.linearResidual(), 1e-6);
	EXPECT_GT(solver.stationarityResidual(), btScalar(900));
}

TEST_F(DeformableBlockPreconditioner, StationarityAllowsBalancedContactReaction)
{
	ImplicitResidualProbe solver;
	solver.configure(body);
	LagrangeMultiplier lm = {};
	lm.m_num_nodes = 1; lm.m_num_constraints = 1;
	lm.m_indices[0] = 0; lm.m_weights[0] = 1; lm.m_dirs[0] = btVector3(0, 1, 0);
	solver.m_objective->m_projection.m_lagrangeMultipliers.push_back(lm);
	Vectors x, rhs;
	x.resize(5); rhs.resize(5);
	for (int n = 0; n < 5; ++n) { x[n].setZero(); rhs[n].setZero(); }
	rhs[0][1] = 100;
	x[4][0] = 100;
	solver.measure(x, rhs);
	EXPECT_EQ(btScalar(0), solver.linearResidual());
	EXPECT_EQ(btScalar(0), solver.stationarityResidual());
}

TEST_F(DeformableBlockPreconditioner, ElasticOperatorIsSymmetricAndMatchesKnownSolution)
{
	Vectors backup;
	backup.resize(4);
	for (int n = 0; n < 4; ++n) backup[n].setZero();
	btDeformableBackwardEulerObjective objective(bodies, backup);
	objective.updateId();
	objective.setDt(btScalar(0.02));
	objective.setImplicit(true);
	objective.m_lf.push_back(&force);
	objective.m_preconditioner->reinitialize(true);
	const btDeformableTest::OperatorProbe probe = btDeformableTest::probeOperator(objective, 4);
	EXPECT_LT(probe.symmetryError, btScalar(2e-6));
	EXPECT_LT(probe.linearityError, btScalar(2e-6));
	EXPECT_EQ(btScalar(0), probe.repeatabilityError);
	EXPECT_EQ(btScalar(0), probe.zeroNorm);
	EXPECT_GT(probe.curvatureP, btScalar(0));
	EXPECT_GT(probe.curvatureQ, btScalar(0));
	Vectors known, rhs, solution;
	known.resize(4); rhs.resize(4); solution.resize(4);
	for (int n = 0; n < 4; ++n) { known[n] = btVector3(n + 1, n - 2, 3 - n); solution[n].setZero(); }
	objective.multiply(known, rhs);
	const btScalar target = btMax(btScalar(1e-9), btScalar(200) * SIMD_EPSILON);
	const btDeformableTest::PCGComparison comparison = btDeformableTest::comparePCG(objective, solution, rhs, 100, target);
	EXPECT_FALSE(comparison.breakdown);
	EXPECT_LE(comparison.finalResidual, target * 10);
	for (int n = 0; n < 4; ++n)
		for (int d = 0; d < 3; ++d) EXPECT_NEAR((double)known[n][d], (double)solution[n][d], 1e-4);
}

TEST_F(DeformableBlockPreconditioner, InternalForcesBalanceApartFromMassDamping)
{
	const btScalar dt = btScalar(0.02);
	body->m_tetraScratches[0].m_F = body->m_tetraScratches[0].m_corotation * btMatrix3x3(btScalar(1.2), 0, 0, 0, btScalar(0.8), 0, 0, 0, btScalar(1.1));
	for (int n = 0; n < 4; ++n) body->m_nodes[n].m_v = btVector3(n + 1, n * 2, 1 - n);
	for (int frozen = 0; frozen < 2; ++frozen)
	{
		body->m_nodes[0].m_frozen = frozen;
		Vectors impulse;
		impulse.resize(4);
		for (int n = 0; n < 4; ++n) impulse[n].setZero();
		force.addScaledForces(dt, impulse);
		btVector3 expected(0, 0, 0);
		for (int n = 0; n < 4; ++n)
			if (body->m_nodes[n].m_frozen <= 0)
				expected -= dt * force.m_damping_alpha * body->m_nodes[n].m_v / body->m_nodes[n].m_im;
		const btVector3 actual = btDeformableTest::sumOnBody(*body, impulse, false);
		for (int d = 0; d < 3; ++d) EXPECT_NEAR((double)expected[d], (double)actual[d], 1e-4);
	}
}

struct NegativeOperator
{
	void multiply(const Vectors& x, Vectors& out) { for (int n = 0; n < x.size(); ++n) out[n] = -x[n]; }
	void precondition(const Vectors& x, Vectors& out) { out = x; }
};

TEST_F(DeformableBlockPreconditioner, TranslationPreservesMomentumWithUnequalMassesAndLocalizedLoad)
{
	Vectors backup, rhs, x, product;
	backup.resize(4, btVector3(0, 0, 0)); rhs = backup; x = backup; product = backup;
	// Apply an impulse to ONE node. The fixture masses are 1, 2, 3, 4.
	rhs[1] = btVector3(2, -3, 5);
	const Vectors originalRhs = rhs;
	btDeformableBackwardEulerObjective objective(bodies, backup);
	objective.updateId(); objective.setDt(btScalar(0.02)); objective.setImplicit(true);
	objective.m_lf.push_back(&force);
	objective.m_preconditioner->reinitialize(true);
	ASSERT_TRUE(objective.setupTranslationCorrection());
	objective.correctTranslation(x, rhs);
	// Even an intentionally incomplete fine solve must retain coarse balance.
	btConjugateResidual<btDeformableBackwardEulerObjective> cr(1);
	cr.solveWithConvergencePolicy(objective, x, rhs, false, true, true);
	objective.multiply(x, product);
	btVector3 defect(0, 0, 0), momentum(0, 0, 0);
	for (int n = 0; n < 4; ++n)
	{
		defect += rhs[n] - product[n];
		momentum += x[n] / body->m_nodes[n].m_im;
		EXPECT_EQ(btScalar(0), (rhs[n] - originalRhs[n]).length2());
	}
	EXPECT_LT(defect.length(), btScalar(2e-4));
	EXPECT_LT((momentum * (1 + btScalar(0.02) * btScalar(0.12)) - rhs[1]).length(), btScalar(2e-4));
	// A converged solution still responds locally, rather than moving rigidly.
	const auto result = btDeformableTest::comparePCG(objective, x, rhs, 100, btScalar(1e-6));
	EXPECT_FALSE(result.breakdown);
	EXPECT_LT(result.finalResidual, btScalar(2e-4));
	EXPECT_GT((x[1] - x[0]).length(), btScalar(0.01));
	// Compare against the original equations with the original preconditioner.
	objective.m_translationCorrection = false;
	Vectors reference = backup;
	const auto baseline = btDeformableTest::comparePCG(objective, reference, rhs, 100, btScalar(1e-6));
	EXPECT_FALSE(baseline.breakdown);
	EXPECT_LT(baseline.finalResidual, btScalar(2e-4));
	for (int n = 0; n < 4; ++n) EXPECT_LT((x[n] - reference[n]).length(), btScalar(2e-4));
}

TEST_F(DeformableBlockPreconditioner, TranslationPreconditionerIsSymmetricAndSolvesCoarseModes)
{
	Vectors backup, p, q, bp, bq, product;
	backup.resize(4, btVector3(0, 0, 0)); p = backup; q = backup; bp = backup; bq = backup; product = backup;
	btDeformableBackwardEulerObjective objective(bodies, backup);
	objective.updateId(); objective.setDt(btScalar(0.02)); objective.setImplicit(true);
	objective.m_lf.push_back(&force); objective.m_preconditioner->reinitialize(true);
	ASSERT_TRUE(objective.setupTranslationCorrection());
	for (int n = 0; n < 4; ++n) { p[n] = btVector3(n+1, 2-n, n*n); q[n] = btVector3(3-n, n+2, -n); }
	objective.precondition(p, bp); objective.precondition(q, bq);
	const btScalar pbq = btDeformableTest::dot(p, bq), qbp = btDeformableTest::dot(q, bp);
	EXPECT_NEAR((double)pbq, (double)qbp, 2e-5 * btMax(btScalar(1), btFabs(pbq)));
	EXPECT_GT(btDeformableTest::dot(p, bp), btScalar(0));
	for (int n = 0; n < 4; ++n) p[n] = btVector3(2, -1, 3);
	objective.multiply(p, product); objective.precondition(product, bp);
	for (int n = 0; n < 4; ++n) EXPECT_LT((bp[n] - p[n]).length(), btScalar(2e-4));
}

TEST_F(DeformableBlockPreconditioner, TranslationDisablesForExplicitFixedAndContact)
{
	Vectors backup; backup.resize(4, btVector3(0, 0, 0));
	btDeformableBackwardEulerObjective objective(bodies, backup);
	objective.updateId(); objective.setDt(btScalar(0.02));
	EXPECT_FALSE(objective.setupTranslationCorrection());
	objective.setImplicit(true);
	ASSERT_TRUE(objective.setupTranslationCorrection());
	body->m_nodes[0].m_frozen = 1;
	EXPECT_FALSE(objective.setupTranslationCorrection());
	EXPECT_FALSE(objective.m_translationCorrection);
	body->m_nodes[0].m_frozen = 0;
	const btScalar inverseMass = body->m_nodes[0].m_im;
	body->m_nodes[0].m_im = 0;
	EXPECT_FALSE(objective.setupTranslationCorrection());
	body->m_nodes[0].m_im = inverseMass;
	LagrangeMultiplier lm = {};
	lm.m_num_nodes = 1; lm.m_indices[0] = 0;
	objective.m_projection.m_lagrangeMultipliers.push_back(lm);
	EXPECT_FALSE(objective.setupTranslationCorrection());
	objective.m_projection.m_lagrangeMultipliers.clear();

}


TEST_F(DeformableBlockPreconditioner, TranslationBalanceSurvivesStiffTruncatedSolve)
{
	force.setLameParameters(btScalar(357142.856), btScalar(1428571.53));
	force.setDamping(btScalar(0.1), btScalar(0.01));
	for (int n = 0; n < 4; ++n) body->setMass(n, btScalar(n + 1) * btScalar(1e-5));
	Vectors backup, rhs, x, product;
	backup.resize(4, btVector3(0, 0, 0)); rhs = backup; x = backup; product = backup;
	rhs[1] = btVector3(0, btScalar(0.018), btScalar(-0.007));
	btDeformableBackwardEulerObjective objective(bodies, backup);
	objective.updateId(); objective.setDt(btScalar(0.0002)); objective.setImplicit(true);
	objective.m_lf.push_back(&force); objective.m_preconditioner->reinitialize(true);
	ASSERT_TRUE(objective.setupTranslationCorrection());
	objective.correctTranslation(x, rhs);
	btConjugateResidual<btDeformableBackwardEulerObjective> cr(3);
	cr.solveWithConvergencePolicy(objective, x, rhs, false, true, true);
	objective.correctTranslation(x, rhs);
	objective.multiply(x, product);
	btVector3 defect(0, 0, 0), momentum(0, 0, 0);
	for (int n = 0; n < 4; ++n) { defect += rhs[n] - product[n]; momentum += x[n] / body->m_nodes[n].m_im; }
	// Allow float cancellation in this ill-conditioned operator; double builds
	// must resolve the net impulse much more accurately than the local residual.
	const btScalar tolerance = sizeof(btScalar) == sizeof(double) ? btScalar(1e-10) : btScalar(1e-4);
	EXPECT_LT(defect.length(), tolerance);
	EXPECT_LT((momentum * (1 + btScalar(0.0002) * btScalar(0.1)) - rhs[1]).length(), tolerance);
}


TEST_F(DeformableBlockPreconditioner, MultipleBodiesBalanceSeparatelyAndSkipOnlyConstrainedBody)
{
	const btVector3 positions[] = {btVector3(5,0,0), btVector3(7,0,0), btVector3(5,1,0), btVector3(5,0,3)};
	const btScalar masses[] = {7, 3, 2, 5};
	btSoftBody other(&info, 4, positions, masses);
	other.appendTetra(0, 1, 2, 3); other.initializeDmInverse();
	other.m_tetraScratches.resize(1);
	other.m_tetraScratches[0].m_corotation.setIdentity(); other.m_tetraScratches[0].m_J = 1;
	bodies.push_back(&other); force.addSoftBody(&other);
	Vectors backup, rhs, x, product, baseline;
	backup.resize(8, btVector3(0,0,0)); rhs = backup; x = backup; product = backup; baseline = backup;
	btDeformableBackwardEulerObjective objective(bodies, backup);
	objective.updateId(); objective.setDt(btScalar(0.02)); objective.setImplicit(true);
	objective.m_lf.push_back(&force); objective.m_preconditioner->reinitialize(true);
	ASSERT_TRUE(objective.setupTranslationCorrection());
	ASSERT_EQ(2, objective.m_translationBodies.size());
	rhs[1] = btVector3(2,-3,5); rhs[6] = btVector3(-7,4,1);
	objective.correctTranslation(x, rhs);
	btConjugateResidual<btDeformableBackwardEulerObjective> cr(2);
	cr.solveWithConvergencePolicy(objective, x, rhs, false, true, true);
	objective.multiply(x, product);
	for (int b = 0; b < 2; ++b)
	{
		btVector3 defect(0,0,0);
		for (int n = b*4; n < b*4+4; ++n) defect += rhs[n] - product[n];
		EXPECT_LT(defect.length(), btScalar(2e-4));
	}
	// A multiplier on body 0 must not suppress body 1's correction. This also
	// exercises a nonzero node offset and the extended KKT vector length.
	LagrangeMultiplier lm = {};
	lm.m_num_nodes = 1; lm.m_indices[0] = 1; lm.m_num_constraints = 1;
	lm.m_weights[0] = 1; lm.m_dirs[0] = btVector3(1,0,0);
	objective.m_projection.m_lagrangeMultipliers.push_back(lm);
	objective.m_preconditioner->reinitialize(true);
	rhs.resize(9, btVector3(0,0,0)); x.resize(9); product.resize(9); baseline.resize(9);
	for (int n = 0; n < 9; ++n) x[n].setZero();
	ASSERT_TRUE(objective.setupTranslationCorrection());
	ASSERT_EQ(1, objective.m_translationBodies.size());
	EXPECT_EQ(4, objective.m_translationBodies[0].offset);
	objective.m_preconditioner->operator()(rhs, baseline);
	objective.precondition(rhs, product);
	for (int n = 0; n < 4; ++n) EXPECT_EQ(btScalar(0), (baseline[n] - product[n]).length2());
	EXPECT_EQ(btScalar(0), (baseline[8] - product[8]).length2());
	objective.correctTranslation(x, rhs);
	for (int n = 0; n < 4; ++n) EXPECT_EQ(btScalar(0), x[n].length2());
	cr.solveWithConvergencePolicy(objective, x, rhs, false, true, true);
	objective.multiply(x, product);
	btVector3 defect(0,0,0);
	for (int n = 4; n < 8; ++n) defect += rhs[n] - product[n];
	EXPECT_LT(defect.length(), btScalar(2e-4));
	// A contact joining the two bodies excludes both, then removal restores both.
	objective.m_projection.m_lagrangeMultipliers[0].m_num_nodes = 2;
	objective.m_projection.m_lagrangeMultipliers[0].m_indices[1] = 4;
	EXPECT_FALSE(objective.setupTranslationCorrection());
	objective.m_projection.m_lagrangeMultipliers.clear();
	ASSERT_TRUE(objective.setupTranslationCorrection());
	EXPECT_EQ(2, objective.m_translationBodies.size());
	body->m_nodes[0].m_frozen = 1;
	ASSERT_TRUE(objective.setupTranslationCorrection());
	ASSERT_EQ(1, objective.m_translationBodies.size());
	EXPECT_EQ(4, objective.m_translationBodies[0].offset);
	body->m_nodes[0].m_frozen = 0;
	force.removeSoftBody(&other); bodies.pop_back(); objective.updateId();
	ASSERT_TRUE(objective.setupTranslationCorrection());
	EXPECT_EQ(1, objective.m_translationBodies.size());
}


TEST(DeformableCRConvergence, WeightedResidualRetainsUsefulStepRejectedByInfinityNorm)
{
	struct Matrix
	{
		btMatrix3x3 a;
		void multiply(const Vectors& x, Vectors& y) { y[0] = a*x[0]; }
		void precondition(const Vectors& x, Vectors& y) { y=x; }
	} matrix = {btMatrix3x3(12,-2,0, -2,1,0, 0,0,1)};
	// SPD matrix, r=(1,1,0), A*r=(10,-1,0), alpha=9/101.
	// CR reduces ||r||_2^2 from 2 to 121/101, but raises ||r||_inf
	// from 1 to 110/101. The legacy best-iterate selection discards the step.
	Vectors rhs, oldX, newX, product;
	rhs.resize(1,btVector3(1,1,0)); oldX.resize(1,btVector3(0,0,0)); newX=oldX; product=oldX;
	btConjugateResidual<Matrix> legacy(1), weighted(1);
	legacy.solveWithConvergencePolicy(matrix,oldX,rhs,false,true,false);
	weighted.solveWithConvergencePolicy(matrix,newX,rhs,false,true,true);
	EXPECT_EQ(btScalar(0),oldX[0].length2());
	EXPECT_NEAR((double)newX[0].x(),9.0/101,1e-6);
	EXPECT_NEAR((double)newX[0].y(),9.0/101,1e-6);
	matrix.multiply(newX,product);
	EXPECT_NEAR((double)(rhs[0]-product[0]).length2(),121.0/101,1e-6);
	EXPECT_NEAR((double)weighted.getFinalResidual(),std::sqrt(121.0/101),1e-6);
	// Continuing from that nonzero guess solves the original equations.
	btConjugateResidual<Matrix> converged(8);
	converged.solveWithConvergencePolicy(matrix,newX,rhs,false,true,true);
	EXPECT_LT((newX[0]-btVector3(btScalar(0.375),btScalar(1.75),0)).length(),btScalar(1e-5));
}

TEST_F(DeformableBlockPreconditioner, WeightedCRSolvesElasticOperatorWithTranslationCorrection)
{
	Vectors backup, known, rhs, x;
	backup.resize(4,btVector3(0,0,0)); known=backup;rhs=backup;x=backup;
	btDeformableBackwardEulerObjective objective(bodies,backup);
	objective.updateId();objective.setDt(btScalar(0.02));objective.setImplicit(true);
	objective.m_lf.push_back(&force);objective.m_preconditioner->reinitialize(true);
	ASSERT_TRUE(objective.setupTranslationCorrection());
	for(int n=0;n<4;++n) known[n]=btVector3(n+1,2-n,3*n);
	objective.multiply(known,rhs);objective.correctTranslation(x,rhs);
	btConjugateResidual<btDeformableBackwardEulerObjective> cr(100);
	cr.solveWithConvergencePolicy(objective,x,rhs,false,false,true);
	for(int n=0;n<4;++n) EXPECT_LT((x[n]-known[n]).length(),btScalar(1e-4));
}

TEST(DeformableCRContinuation, PhysicalTargetIgnoresMisleadingPreconditionedTolerance)
{
	struct Matrix
	{
		void multiply(const Vectors& x, Vectors& y) { y[0]=x[0]*btVector3(1,10,100); }
		void precondition(const Vectors& x, Vectors& y) { y[0]=x[0]*btScalar(1e-20); }
	} matrix;
	Vectors x,rhs,product;
	x.resize(1,btVector3(0,0,0));rhs.resize(1,btVector3(1,1,1));product=x;
	btConjugateResidual<Matrix> cr(300);
	const int iterations=cr.solveWithConvergencePolicy(matrix,x,rhs,false,false,true,1200,btScalar(1e-5));
	matrix.multiply(x,product);
	EXPECT_GT(iterations,0);EXPECT_LT(iterations,300);
	EXPECT_LE((rhs[0]-product[0]).length(),btScalar(1e-5));
	EXPECT_EQ(300,cr.m_maxIterations); // Per-call limit must not leak to explicit solves.
	EXPECT_EQ(0,cr.solveWithConvergencePolicy(matrix,x,rhs,false,false,true,1200,btScalar(1e-5)));
}

TEST(DeformableCRContinuation, PreservesDirectionsBeyondDefaultBudget)
{
	struct Matrix
	{
		void multiply(const Vectors& x, Vectors& y)
		{
			const int size=3*x.size();
			for(int i=0;i<size;++i)
			{
				btScalar value=2*x[i/3][i%3];
				if(i>0)value-=x[(i-1)/3][(i-1)%3];
				if(i+1<size)value-=x[(i+1)/3][(i+1)%3];
				y[i/3][i%3]=value;
			}
		}
		void precondition(const Vectors& x,Vectors& y) { y=x; }
	} matrix;
	Vectors x,rhs,product;
	x.resize(128,btVector3(0,0,0));rhs=x;product=x;rhs[0][0]=1;
	btConjugateResidual<Matrix> cr(300);
	const int iterations=cr.solveWithConvergencePolicy(matrix,x,rhs,false,false,true,1200,btScalar(1e-5));
	matrix.multiply(x,product);
	btScalar error=0;
	for(int n=0;n<x.size();++n)error+=(rhs[n]-product[n]).length2();
	EXPECT_GT(iterations,300);EXPECT_LT(iterations,1200);
	EXPECT_LE(btSqrt(error),btScalar(1e-5));
}


TEST_F(DeformableBlockPreconditioner, NewtonPreservesContactVelocityInsteadOfPreviousFallingVelocity)
{
	ImplicitResidualProbe solver; solver.configure(body);
	for (int n=0;n<4;++n) { body->m_nodes[n].m_vn=btVector3(0,0,-10); body->m_nodes[n].m_v=btVector3(0,0,-10); }
	// A plane contact stopped node 0; previous-step velocity still points down.
	body->m_nodes[0].m_v=btVector3(0,0,0);
	LagrangeMultiplier lm = {};
	lm.m_num_nodes=1;lm.m_indices[0]=0;lm.m_weights[0]=1;lm.m_num_constraints=1;lm.m_dirs[0]=btVector3(0,0,1);
	solver.m_objective->m_projection.m_lagrangeMultipliers.push_back(lm);
	solver.prepareContactStep();
	Vectors dv=solver.m_objective->m_implicitConstraintDv, residual, rhs, correction;
	residual.resize(4,btVector3(0,0,0)); correction.resize(5,btVector3(0,0,0));
	EXPECT_EQ(btScalar(10),dv[0].z());
	EXPECT_EQ(btScalar(0),solver.constraintError(dv));
	// Include both gravity and a tangential load: only the normal is constrained.
	residual[0]=btVector3(1,0,-2);
	solver.m_objective->addLagrangeMultiplierRHS(residual,dv,rhs);
	EXPECT_EQ(btScalar(0),rhs[4][0]);
	solver.m_objective->m_preconditioner->reinitialize(true);
	btConjugateResidual<btDeformableBackwardEulerObjective> cr(100);
	cr.solveWithConvergencePolicy(*solver.m_objective,correction,rhs,false,false,true);
	const btVector3 finalVelocity=body->m_nodes[0].m_vn+dv[0]+correction[0];
	EXPECT_NEAR(0.0,(double)finalVelocity.z(),1e-5);
	EXPECT_GT(finalVelocity.x(),btScalar(0));
	// Later Newton errors are corrected toward the saved contact velocity.
	dv[0][2]-=3;
	solver.m_objective->addLagrangeMultiplierRHS(residual,dv,rhs);
	EXPECT_EQ(btScalar(3),rhs[4][0]);
	EXPECT_EQ(btScalar(3),solver.constraintError(dv));
	// A moving contact target is preserved too, and rebuilt next timestep.
	body->m_nodes[0].m_v=btVector3(0,0,2);
	solver.prepareContactStep();
	EXPECT_EQ(btScalar(12),solver.m_objective->m_implicitConstraintDv[0].z());
}

}  // namespace

int main(int argc, char** argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
