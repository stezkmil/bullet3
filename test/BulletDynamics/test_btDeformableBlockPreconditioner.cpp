#include "BulletSoftBody/btDeformableBodySolver.h"
#include "BulletSoftBody/btDeformableDiagnostics.h"
#include "LinearMath/btQuaternion.h"
#include "deformable_solver_test_helpers.h"
#include <gtest/gtest.h>
#include <limits>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <thread>
#include "BulletSoftBody/btDeformableContactForce.h"
#include "BulletSoftBody/btDeformableContactRefresh.h"
#include "BulletSoftBody/btDeformableEnergyChange.h"
#include "BulletSoftBody/btDeformableVolumeBarrierForce.h"
#include "BulletSoftBody/btDeformableNodalForce.h"
#include "BulletSoftBody/btDeformableNewtonSnapshot.h"
#include "BulletCollision/Gimpact/btGImpactShape.h"
#include "BulletCollision/Gimpact/btGImpactVertexCache.h"
#include "BulletCollision/CollisionShapes/btTriangleIndexVertexArray.h"
#include "BulletSoftBody/btSoftBodyRigidBodyCollisionConfiguration.h"
#include "BulletCollision/BroadphaseCollision/btDbvtBroadphase.h"
#include "BulletCollision/CollisionDispatch/btCollisionDispatcher.h"

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

TEST_F(DeformableBlockPreconditioner, FinalContactNewtonUpdateIsChecked)
{
	btDeformableBodySolver solver;solver.setImplicit(true);solver.setLineSearch(true);solver.setMaxNewtonIterations(1);
	solver.reinitialize(bodies,.01);solver.setPreconditioner(btDeformableBackwardEulerObjective::KKT_preconditioner);
	btDeformableContactForce contact(.01);btSoftBody::DeformableNodeNodeContact c={};
	c.m_normal=btVector3(1,0,0);btSoftBody::ContactNode node={&body->m_nodes[0],btMatrix3x3::getIdentity()};c.m_surfaceNodes.push_back(node);
	ASSERT_TRUE(contact.add(c));contact.contacts[0].normalImpulse=.1;
	solver.m_objective->m_lf.push_back(&contact);solver.setupDeformableSolve(true);
	solver.solveDeformableConstraints(.01);EXPECT_TRUE(solver.m_lastSolveConverged);
	EXPECT_NEAR(double(body->m_nodes[0].m_v.x()),.1/11,1e-9);
	solver.m_objective->m_lf.clear();
}

TEST_F(DeformableBlockPreconditioner, PreparedBlocksMatchFreshSolveAndDoNotPersist)
{
	btDeformableBodySolver solver;
	solver.setImplicit(true); solver.m_useProjection = false;
	solver.reinitialize(bodies, btScalar(.01));
	solver.setPreconditioner(btDeformableBackwardEulerObjective::KKT_preconditioner);
	solver.m_objective->m_lf.push_back(&force);
	solver.setupDeformableSolve(true); solver.updateState();
	Vectors rhs, reused, fresh, changed;
	rhs.resize(4, btVector3(0,0,0)); rhs[0] = btVector3(1,-2,3); rhs[1] = -rhs[0];
	reused.resize(4, btVector3(0,0,0)); fresh = reused; changed = reused;
	solver.m_objective->m_preconditioner->reinitialize(true);
	solver.computeStep(reused, rhs, true);
	solver.computeStep(fresh, rhs);
	for (int n = 0; n < 4; ++n) EXPECT_LT((reused[n]-fresh[n]).length(), btScalar(1e-8));
	force.setYoungsModulus(btScalar(1e7));
	solver.computeStep(changed, rhs);
	btScalar difference = 0;
	for (int n = 0; n < 4; ++n) difference += (changed[n]-fresh[n]).length2();
	EXPECT_GT(difference, btScalar(1e-6));
	solver.m_objective->m_lf.clear();
}

TEST_F(DeformableBlockPreconditioner, EnergyChangeRemovesLargeLinearPotentialOffsets)
{
	Vectors zero,trial;zero.resize(4,btVector3(0,0,0));trial=zero;trial[0]=btVector3(1e-6,0,0);
	btDeformableBackwardEulerObjective o(bodies,zero);o.updateId();
	body->m_nodes[0].m_q.setX(btScalar(1e12));
	btAlignedObjectArray<int> indices;indices.push_back(0);
	btDeformableNodalForce nodal(body,indices,btVector3(100,0,0));o.m_lf.push_back(&nodal);
	btDeformableEnergyChange change(o,zero,.01);
	EXPECT_NEAR(change.difference(trial),-1e-6+0.5e-12,1e-18);
	o.m_lf.clear();btDeformableGravityForce gravity(btVector3(100,0,0));gravity.addSoftBody(body);o.m_lf.push_back(&gravity);
	btDeformableEnergyChange gravityChange(o,zero,.01);
	EXPECT_NEAR(gravityChange.difference(trial),-1e-6+0.5e-12,1e-18);
	EXPECT_EQ(0.,gravityChange.difference(zero));
}

TEST_F(DeformableBlockPreconditioner, VelocityGuessPreservesMomentumReferenceAndFixedNodes)
{
	btDeformableBodySolver solver;solver.setImplicit(true);solver.setMaxNewtonIterations(5);
	solver.reinitialize(bodies,.01);solver.setPreconditioner(btDeformableBackwardEulerObjective::KKT_preconditioner);
	for(int n=0;n<4;++n)body->m_nodes[n].m_v=body->m_nodes[n].m_vn=btVector3(2,0,0);
	solver.setupDeformableSolve(true);solver.updateState();
	const btVector3 position=body->m_nodes[0].m_x;
	body->m_nodes[3].m_frozen=1;
	Vectors guess;guess.resize(4,btVector3(7,1,0));
	ASSERT_TRUE(solver.setImplicitVelocityGuess(guess));
	EXPECT_EQ(guess[0],body->m_nodes[0].m_v);EXPECT_EQ(btVector3(2,0,0),body->m_nodes[3].m_v);
	EXPECT_EQ(position,body->m_nodes[0].m_x);EXPECT_EQ(btVector3(2,0,0),body->m_nodes[0].m_vn);
	guess[0].setX(std::numeric_limits<btScalar>::quiet_NaN());
	EXPECT_FALSE(solver.setImplicitVelocityGuess(guess));EXPECT_EQ(btVector3(7,1,0),body->m_nodes[0].m_v);
	solver.solveDeformableConstraints(.01);EXPECT_TRUE(solver.m_lastSolveConverged);
	EXPECT_LT((body->m_nodes[0].m_v-btVector3(2,0,0)).length(),btScalar(1e-10));
}

TEST_F(DeformableBlockPreconditioner, VelocityGuessRejectsInversionAndRestoresColdGuess)
{
	btDeformableBodySolver solver;solver.setImplicit(true);solver.reinitialize(bodies,.01);
	solver.setupDeformableSolve(true);solver.updateState();
	btDeformableVolumeBarrierForce barrier;btDeformableVolumeBarrierForce::Material material={body,200};barrier.materials.push_back(material);
	solver.m_objective->m_lf.push_back(&barrier);
	const btVector3 q=body->m_nodes[3].m_q,v=body->m_nodes[3].m_v;
	Vectors guess;guess.resize(4,btVector3(0,0,0));guess[3]=btVector3(0,0,-1000);
	EXPECT_FALSE(solver.setImplicitVelocityGuess(guess));EXPECT_EQ(q,body->m_nodes[3].m_q);EXPECT_EQ(v,body->m_nodes[3].m_v);
	EXPECT_TRUE(barrier.admissible());solver.m_useProjection=true;EXPECT_FALSE(solver.setImplicitVelocityGuess(guess));
	solver.m_objective->m_lf.clear();
}

TEST_F(DeformableBlockPreconditioner, NewtonSnapshotReplaysEnergiesAndLinearCorrection)
{
	for (int frozenSnapshot = 0; frozenSnapshot < 2; ++frozenSnapshot)
	{
		class Probe : public btDeformableBodySolver
		{
		  public:
			void prepare()
			{
				setupDeformableSolve(true);
				updateState();
				m_objective->computeResidual(m_dt, m_residual);
				computeStep(m_ddv, m_residual);
			}
			btScalar slope()
			{
				return m_cg.dot(m_ddv, m_residual);
			}
			const Vectors &direction() const
			{
				return m_ddv;
			}
		} source;
		const btScalar dt = btScalar(.0002);
		body->m_nodes[3].m_x.setZ(btScalar(.9));
		for (int n = 0; n < 4; ++n)
			body->m_nodes[n].m_v = body->m_nodes[n].m_vn = btVector3(btScalar(.1) * n, btScalar(.2), 0);
		source.setImplicit(true);
		source.setLineSearch(true);
		source.reinitialize(bodies, dt);
		source.setPreconditioner(btDeformableBackwardEulerObjective::KKT_preconditioner);
		source.m_objective->m_lf.push_back(&force);
		btDeformableVolumeBarrierForce barrier;
		btDeformableVolumeBarrierForce::Material material = {body, btScalar(200)};
		barrier.materials.push_back(material);
		source.m_objective->m_lf.push_back(&barrier);
		btAlignedObjectArray<int> indices;
		indices.push_back(1);
		indices.push_back(3);
		btDeformableNodalForce nodal(body, indices, btVector3(1, 2, 3));
		source.m_objective->m_lf.push_back(&nodal);
		btDeformableGravityForce gravity(btVector3(0, -1, 0));
		gravity.addSoftBody(body);
		source.m_objective->m_lf.push_back(&gravity);
		btDeformableContactForce contact(dt);
		btSoftBody::DeformableNodeNodeContact c = {};
		c.m_normal = btVector3(1, 0, 0);
		c.m_friction = btScalar(.5);
		btSoftBody::ContactNode cn = {&body->m_nodes[0], btMatrix3x3::getIdentity()};
		c.m_surfaceNodes.push_back(cn);
		ASSERT_TRUE(contact.add(c));
		contact.contacts[0].normalImpulse = btScalar(.01);
		contact.contacts[0].tangentImpulse = btVector3(0, btScalar(.001), 0);
		source.m_objective->m_lf.push_back(&contact);
		body->m_tetraScratchesTn.clear();
		if (frozenSnapshot)
		{
			body->m_tetraScratchesTn.resize(1);
			body->advanceDeformation();
		}
		source.prepare();
		btDeformableNewtonSnapshot snapshot;
		snapshot.dt = dt;
		snapshot.slope = source.slope();
		snapshot.baselineEnergy = source.m_objective->totalEnergy(dt) + source.kineticEnergy();
		source.backupDv();
		for (btScalar scale : {btScalar(1), btScalar(.25), btScalar(1e-6)})
		{
			source.updateEnergy(scale);
			snapshot.scales.push_back(scale);
			snapshot.energies.push_back(source.m_objective->totalEnergy(dt) + source.kineticEnergy());
		}
		source.revertDv();
		source.updateState();
		const char *temp = std::getenv("TEMP");
		const std::string path = std::string(temp ? temp : ".") + "/bullet-newton-" +
								 std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()) + ".bin";
		ASSERT_TRUE(snapshot.save(path.c_str(), source));
		btDeformableNewtonSnapshot replay;
		ASSERT_TRUE(replay.load(path.c_str()));
		Vectors product, conditioned;
		btDeformableNewtonSnapshot::probe(replay.solver, product, conditioned);
		ASSERT_EQ(product.size(), replay.operatorProbe.size());
		ASSERT_EQ(conditioned.size(), replay.preconditionerProbe.size());
		for (int n = 0; n < product.size(); ++n)
			for (int d = 0; d < 3; ++d)
			{
				EXPECT_EQ(product[n][d], replay.operatorProbe[n][d]);
				EXPECT_EQ(conditioned[n][d], replay.preconditionerProbe[n][d]);
			}
		EXPECT_NEAR(double(snapshot.baselineEnergy), double(replay.energy()), 1e-10);
		replay.solver.backupDv();
		for (int i = 0; i < snapshot.scales.size(); ++i)
		{
			replay.solver.updateEnergy(snapshot.scales[i]);
			EXPECT_NEAR(double(snapshot.energies[i]), double(replay.energy()), 1e-9);
		}
		replay.solver.revertDv();
		replay.solver.updateState();
		replay.recomputeDirection();
		for (int n = 0; n < 4; ++n)
			EXPECT_LT((source.direction()[n] - replay.direction()[n]).length(), btScalar(1e-9));
		// A truncated capture must be rejected, not interpreted as a partial scene.
		FILE *invalid = std::fopen(path.c_str(), "wb");
		ASSERT_TRUE(invalid != nullptr);
		std::fputc(0, invalid);
		std::fclose(invalid);
		btDeformableNewtonSnapshot truncated;
		EXPECT_FALSE(truncated.load(path.c_str()));
		std::remove(path.c_str());
		source.m_objective->m_lf.clear();
	}
}

TEST(NewtonReplay, DISABLED_CapturedLineSearch)
{
	const char* path=std::getenv("BULLET_NEWTON_REPLAY");ASSERT_TRUE(path && *path);
	btDeformableNewtonSnapshot replay;ASSERT_TRUE(replay.load(path));
	if (replay.operatorProbe.size())
	{
		Vectors product,conditioned;btDeformableNewtonSnapshot::probe(replay.solver,product,conditioned);
		for(int kind=0;kind<2;++kind)
		{
			const Vectors& expected=kind?replay.preconditionerProbe:replay.operatorProbe;
			const Vectors& actual=kind?conditioned:product;btScalar difference=0,norm=0;
			for(int n=0;n<actual.size();++n){difference+=(actual[n]-expected[n]).length2();norm+=expected[n].length2();}
			printf("OPERATOR_AUDIT kind=%s relative_error=%.17g\n",kind?"preconditioner":"multiply",double(btSqrt(difference/btMax(norm,btScalar(1e-30)))));
			EXPECT_LE(btSqrt(difference),btScalar(1e-12)*btMax(btSqrt(norm),btScalar(1e-30)));
		}
	}
	const char* assembly=std::getenv("BULLET_DEFORMABLE_ASSEMBLED_ELASTIC");
	ASSERT_EQ(replay.assembledElastic, int(!(assembly && assembly[0]=='0')));
	printf("REPLAY step=%lld h=%.17g baseline=%.17g captured=%.17g slope=%.17g trials=%d\n",replay.step,double(replay.dt),double(replay.energy()),double(replay.baselineEnergy),double(replay.slope),replay.scales.size());
	replay.solver.m_objective->m_KKTPreconditioner->reinitialize(true);
	printf("BASE_RESIDUAL %.17g\n",double(replay.weightedResidual()));
	EXPECT_NEAR(double(replay.baselineEnergy),double(replay.energy()),1e-11*btMax(1.,std::abs(double(replay.baselineEnergy))));
	replay.solver.backupDv();
	for(int i=0;i<replay.scales.size();++i)
	{
		replay.solver.updateEnergy(replay.scales[i]);const btScalar energy=replay.energy();
		const btScalar limit=replay.baselineEnergy-btScalar(.01)*replay.scales[i]*replay.slope+SIMD_EPSILON;
		printf("TRIAL scale=%.17g energy=%.17g captured=%.17g delta=%.17g accepted=%d\n",double(replay.scales[i]),double(energy),double(replay.energies[i]),double(energy-replay.baselineEnergy),int(energy<limit));
		if(i<3)printf("TRIAL_RESIDUAL %.17g\n",double(replay.weightedResidual()));
		EXPECT_NEAR(double(replay.energies[i]),double(energy),1e-11*btMax(1.,std::abs(double(energy))));
		EXPECT_EQ(replay.energies[i]<limit,energy<limit);
	}
	replay.solver.revertDv();replay.solver.updateState();
	const Vectors saved=replay.direction();replay.recomputeDirection();btScalar error=0,length=0;
	for(int n=0;n<saved.size();++n){error+=(saved[n]-replay.direction()[n]).length2();length+=saved[n].length2();}
	printf("LINEAR_REPLAY relative_direction_error=%.17g\n",double(btSqrt(error/btMax(length,btScalar(1e-30)))));
	EXPECT_LE(btSqrt(error),btScalar(1e-6)*btMax(btScalar(1),btSqrt(length)));
	replay.solver.solveDeformableConstraints(replay.dt);
	printf("RESUMED_SOLVE converged=%d\n",int(replay.solver.m_lastSolveConverged));
	if (replay.scales.size()) EXPECT_TRUE(replay.solver.m_lastSolveConverged);
}

TEST(NewtonReplay, DISABLED_CapturedRoundTrip)
{
	const char* path=std::getenv("BULLET_NEWTON_REPLAY");ASSERT_TRUE(path && *path);
	btDeformableNewtonSnapshot first;ASSERT_TRUE(first.load(path));first.recomputeDirection();
	const Vectors expected=first.direction();
	first.recomputeDirection();
	for(int n=0;n<expected.size();++n)for(int d=0;d<3;++d)EXPECT_EQ(expected[n][d],first.direction()[n][d]);
	const std::string copy=std::string(path)+".roundtrip.bin";
	ASSERT_TRUE(first.save(copy.c_str(),first.solver));
	btDeformableNewtonSnapshot second;ASSERT_TRUE(second.load(copy.c_str()));second.recomputeDirection();
	for(int n=0;n<expected.size();++n)for(int d=0;d<3;++d)EXPECT_EQ(expected[n][d],second.direction()[n][d]);
	std::remove(copy.c_str());
}

TEST(NewtonReplay, DISABLED_ContactCoarseBenchmark)
{
	const char* path=std::getenv("BULLET_NEWTON_REPLAY");ASSERT_TRUE(path && *path);
	btDeformableNewtonSnapshot replay;ASSERT_TRUE(replay.load(path));
	const char* coarse=std::getenv("BULLET_DEFORMABLE_CONTACT_COARSE");
	if(coarse)replay.solver.m_objective->m_contactCoarseEnabled=coarse[0]=='1';
	const auto start=std::chrono::steady_clock::now();replay.recomputeDirection();
	const double ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
	auto& o=*replay.solver.m_objective;
	for(int f=0;f<o.m_lf.size();++f)o.m_lf[f]->prepareImplicitForceDifferential(replay.dt);
	o.m_KKTPreconditioner->reinitialize(true);
	Vectors product;product.resize(replay.direction().size());o.multiply(replay.direction(),product);
	btScalar raw=0,weighted=0;
	for(int n=0;n<product.size();++n)
	{
		const btVector3 r=replay.residual()[n]-product[n];raw+=r.length2();
		weighted+=r.dot(o.m_KKTPreconditioner->applyInverseNodeBlock(n,r));
	}
	printf("CONTACT_COARSE_BENCH enabled=%d modes=%d ms=%.6f raw=%.17g weighted=%.17g\n",int(o.m_contactCoarse),int(o.m_contactZ.size()),ms,double(btSqrt(raw)),double(btSqrt(weighted)));
	for(int f=0;f<o.m_lf.size();++f)o.m_lf[f]->finishImplicitForceDifferential();
	const auto resumed=std::chrono::steady_clock::now();replay.solver.solveDeformableConstraints(replay.dt);
	printf("CONTACT_COARSE_RESUMED ms=%.6f converged=%d\n",std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-resumed).count(),int(replay.solver.m_lastSolveConverged));
}

TEST(NewtonReplay, DISABLED_NewtonBudgetSweep)
{
	const char* path=std::getenv("BULLET_NEWTON_REPLAY");ASSERT_TRUE(path && *path);
	for(int budget : {5,10,20,40,80})
	{
		btDeformableNewtonSnapshot replay;ASSERT_TRUE(replay.load(path));
		replay.solver.setMaxNewtonIterations(budget);
		const auto start=std::chrono::steady_clock::now();
		replay.solver.solveDeformableConstraints(replay.dt);
		const double ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
		replay.solver.m_objective->m_KKTPreconditioner->reinitialize(true);
		printf("NEWTON_BUDGET budget=%d ms=%.6f converged=%d residual=%.17g energy=%.17g\n",budget,ms,int(replay.solver.m_lastSolveConverged),double(replay.weightedResidual()),double(replay.energy()));
	}
}

TEST(NewtonReplay, DISABLED_ContactNewtonBudgetSweep)
{
	const char* path=std::getenv("BULLET_NEWTON_REPLAY");ASSERT_TRUE(path && *path);
	for(int budget : {5,10,15,20,30})
	{
		btDeformableNewtonSnapshot replay;ASSERT_TRUE(replay.load(path));
		replay.solver.setMaxNewtonIterations(budget);
		btDeformableContactForce* contact=nullptr;
		for(int f=0;f<replay.solver.m_objective->m_lf.size();++f)
			if(replay.solver.m_objective->m_lf[f]->getForceType()==BT_CONTACT_FORCE)
				contact=static_cast<btDeformableContactForce*>(replay.solver.m_objective->m_lf[f]);
		ASSERT_TRUE(contact!=nullptr);
		const auto start=std::chrono::steady_clock::now();bool converged=false;btScalar error=0;int calls=0;
		for(int outer=0;outer<20;++outer)
		{
			replay.solver.solveDeformableConstraints(replay.dt);++calls;
			if(replay.solver.m_lastSolveInvalidPredictor)break;
			error=contact->updateMultipliers();
			printf("CONTACT_BUDGET_ITER budget=%d outer=%d inner_converged=%d error=%.17g\n",budget,outer,int(replay.solver.m_lastSolveConverged),double(error));
			if(replay.solver.m_lastSolveConverged && error<=btScalar(.001)){converged=true;break;}
		}
		printf("CONTACT_BUDGET budget=%d calls=%d ms=%.6f converged=%d error=%.17g\n",budget,calls,std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count(),int(converged),double(error));
	}
}

TEST(NewtonReplay, DISABLED_RequestedCapturePreservesSolve)
{
	ASSERT_TRUE(std::getenv("BULLET_DEFORMABLE_DIAGNOSTICS")!=nullptr);
	ASSERT_TRUE(std::getenv("BULLET_DEFORMABLE_CAPTURE_SOLVER_STEP")!=nullptr);
	const long long requested=std::strtoll(std::getenv("BULLET_DEFORMABLE_CAPTURE_SOLVER_STEP"),nullptr,10);
	btVector3 reference;
	for(int capture=0;capture<2;++capture)
	{
		btSoftBodyWorldInfo info;const btVector3 p(0,0,0);const btScalar mass=1,dt=btScalar(.01);
		btSoftBody body(&info,1,&p,&mass);body.m_nodes[0].m_v=body.m_nodes[0].m_vn=btVector3(0,0,0);
		btAlignedObjectArray<btSoftBody*> bodies;bodies.push_back(&body);
		btDeformableBodySolver solver;solver.setImplicit(true);solver.setLineSearch(true);solver.setMaxNewtonIterations(1);
		solver.reinitialize(bodies,dt);solver.setPreconditioner(btDeformableBackwardEulerObjective::KKT_preconditioner);
		btAlignedObjectArray<int> indices;indices.push_back(0);
		btDeformableNodalForce force(&body,indices,btVector3(1,2,3));solver.m_objective->m_lf.push_back(&force);
		solver.setupDeformableSolve(true);
		{btDeformableDiagnostics::StepScope scope(&solver,requested+(capture?0:1),dt);solver.solveDeformableConstraints(dt);}
		if(!capture)reference=body.m_nodes[0].m_v;else EXPECT_EQ(reference,body.m_nodes[0].m_v);
		solver.m_objective->m_lf.clear();
	}
	const std::string path=std::string(std::getenv("BULLET_DEFORMABLE_DIAGNOSTICS"))+".step-"+std::to_string(requested)+".newton.bin";
	btDeformableNewtonSnapshot replay;ASSERT_TRUE(replay.load(path.c_str()));EXPECT_EQ(0,replay.scales.size());
	EXPECT_NEAR(double(replay.baselineEnergy),double(replay.energy()),1e-12);
	const Vectors saved=replay.direction();replay.recomputeDirection();
	EXPECT_LT((saved[0]-replay.direction()[0]).length(),btScalar(1e-12));
}

TEST(NewtonReplay, DISABLED_AutomaticFailureCapture)
{
	ASSERT_TRUE(std::getenv("BULLET_DEFORMABLE_DIAGNOSTICS")!=nullptr);
	btSoftBodyWorldInfo info;const btVector3 p(0,0,0);const btScalar mass=1,dt=btScalar(.01);
	btSoftBody body(&info,1,&p,&mass);body.m_gravityFactor=2;
	body.m_nodes[0].m_v=body.m_nodes[0].m_vn=btVector3(0,0,0);
	btAlignedObjectArray<btSoftBody*> bodies;bodies.push_back(&body);
	btDeformableBodySolver solver;solver.setImplicit(true);solver.setLineSearch(true);solver.setMaxNewtonIterations(5);
	solver.reinitialize(bodies,dt);solver.setPreconditioner(btDeformableBackwardEulerObjective::KKT_preconditioner);
	btDeformableGravityForce gravity(btVector3(1e8,0,0));gravity.addSoftBody(&body);
	btAlignedObjectArray<int> indices;indices.push_back(0);
	// Gravity's legacy energy omits gravityFactor, making this a deliberate
	// force/energy inconsistency that rejects every line-search trial.
	btDeformableNodalForce nodal(&body,indices,btVector3(-1.5e8,0,0));
	solver.m_objective->m_lf.push_back(&gravity);solver.m_objective->m_lf.push_back(&nodal);
	solver.setupDeformableSolve(true);
	{ btDeformableDiagnostics::StepScope scope(&solver,987654,dt*256);solver.solveDeformableConstraints(dt); }
	EXPECT_FALSE(solver.m_lastSolveConverged);
	const std::string path=std::string(std::getenv("BULLET_DEFORMABLE_DIAGNOSTICS"))+".step-987654.newton.bin";
	btDeformableNewtonSnapshot replay;ASSERT_TRUE(replay.load(path.c_str()));EXPECT_GT(replay.scales.size(),20);
	replay.solver.backupDv();
	for(int i=0;i<replay.scales.size();++i)
	{
		replay.solver.updateEnergy(replay.scales[i]);
		EXPECT_NEAR(double(replay.energies[i]),double(replay.energy()),1e-10*btMax(1.,std::abs(double(replay.energies[i]))));
		EXPECT_FALSE(replay.energy()<replay.baselineEnergy-btScalar(.01)*replay.scales[i]*replay.slope+SIMD_EPSILON);
	}
	solver.m_objective->m_lf.clear();
}

TEST_F(DeformableBlockPreconditioner, MatchesForceDifferentialForRotatedAndFlatElements)
{
	compareBlocksToOperator(btScalar(0.0002), false);
	compareBlocksToOperator(btScalar(0.02), false);
	compareBlocksToOperator(btScalar(0.02), true);
	force.setDamping(0, 0);
	compareBlocksToOperator(btScalar(0.02), true);
	compareBlocksToOperator(0, false);
}


TEST_F(DeformableBlockPreconditioner, DampingEnergyGradientMatchesVelocityResidual)
{
	Vectors velocities, direction, impulse;
	velocities.resize(4); direction.resize(4); impulse.resize(4);
	for (int n = 0; n < 4; ++n)
	{
		velocities[n] = btVector3(n + 1, btScalar(.3)*n, btScalar(-.2));
		direction[n] = btVector3(btScalar(.2), n - 1, btScalar(.4)*n);
	}
	for (btScalar dt : {btScalar(.02), btScalar(.0002), btScalar(.0002/256)})
	{
		for (int n = 0; n < 4; ++n) { body->m_nodes[n].m_v = velocities[n]; impulse[n].setZero(); }
		force.addScaledDampingForce(dt, impulse);
		double expected = 0;
		for (int n = 0; n < 4; ++n) expected -= double(impulse[n].dot(direction[n]));
		const btScalar epsilon = btScalar(1e-4);
		for (int n = 0; n < 4; ++n) body->m_nodes[n].m_v = velocities[n] + epsilon * direction[n];
		const double plus = force.totalDampingEnergy(dt);
		for (int n = 0; n < 4; ++n) body->m_nodes[n].m_v = velocities[n] - epsilon * direction[n];
		const double minus = force.totalDampingEnergy(dt);
		EXPECT_NEAR(expected, (plus - minus) / (2 * epsilon), 1e-5 * btMax(1., std::abs(expected)));
	}
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
	btDeformableDiagnostics::StepScope diagnostics(&solver, 0, dt);
	solver.solveDeformableConstraints(dt);
	btDeformableDiagnostics::bodies("newton_test", bodies, true);
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
	void solveContact(Vectors& x, const Vectors& rhs, int newton)
	{
		m_newtonIteration = newton;
		solveImplicitKKT(x, rhs);
	}
	btScalar linearConstraintResidual() const { return m_lastLinearConstraintResidual; }
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

TEST(DeformableCRContinuation, WeightedTargetCanTightenButNotRelaxDefault)
{
	struct Matrix
	{
		void multiply(const Vectors& x, Vectors& y) { y[0] = x[0] * btScalar(1e-6); }
		void precondition(const Vectors& x, Vectors& y) { y[0] = x[0] * btScalar(1e6); }
	} matrix;
	Vectors x, rhs, product;
	x.resize(1, btVector3(0,0,0)); rhs.resize(1, btVector3(btScalar(1e-12),0,0)); product = x;
	btConjugateResidual<Matrix> cr(20);
	EXPECT_EQ(0, cr.solveWithConvergencePolicy(matrix, x, rhs, false, true, true));
	EXPECT_GT(cr.solveWithConvergencePolicy(matrix, x, rhs, false, true, true, 0, 0, btScalar(1e-11)), 0);
	matrix.multiply(x, product);
	EXPECT_LE((rhs[0] - product[0]).length() * btScalar(1e3), btScalar(1e-11));
	x[0].setZero(); rhs[0] = btVector3(btScalar(1e-9),0,0);
	EXPECT_GT(cr.solveWithConvergencePolicy(matrix, x, rhs, false, true, true, 0, 0, btScalar(1e-3)), 0);
	EXPECT_LE(cr.getTargetResidual(), btScalar(1e-8));
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
	// Exercise the actual first and subsequent Newton solves, preserving the
	// contact velocity while later corrections meet the strict physical target.
	for (int newton = 0; newton < 3; ++newton)
	{
		for (int n = 0; n < correction.size(); ++n) correction[n].setZero();
		solver.solveContact(correction, rhs, newton);
		if (newton > 0)
		{
			EXPECT_LE(solver.linearResidual(), btScalar(5e-5));
			EXPECT_LE(solver.linearConstraintResidual(), btScalar(5e-5));
		}
	}
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


class IntegrationRecoveryProbe : public btDeformableRigidContactConstraint
{
public:
 btVector3 velocity=btVector3(0,0,0), split=btVector3(0,0,0);
 IntegrationRecoveryProbe(const btSoftBody::DeformableRigidContact& c,const btContactSolverInfo& info):btDeformableRigidContactConstraint(c,info){}
 btVector3 getVb() const override{return velocity;}
 btVector3 getSplitVb() const override{return split;}
 btVector3 getDv(const btSoftBody::Node*) const override{return velocity;}
 void applyImpulse(const btVector3& impulse) override{velocity-=impulse;}
 void applySplitImpulse(const btVector3& impulse) override{split-=impulse;}
};
TEST(ContactRecovery, IntegrationModePreservesExplicitLegacyResponse)
{
 btCollisionObject rigid;rigid.setCollisionFlags(btCollisionObject::CF_STATIC_OBJECT);
 btSoftBody::DeformableRigidContact c;
 c.m_cti.m_colObj=&rigid;c.m_cti.m_normal=btVector3(1,0,0);c.m_cti.m_offset=-2;c.m_cti.m_contact_point_impulse_magnitude=nullptr;
 c.m_c0.setIdentity();c.m_c5.setIdentity();c.m_c3=0;
 btContactSolverInfo info;EXPECT_FALSE(info.m_deformable_implicit);
 info.m_timeStep=btScalar(.002);info.m_deformable_cfm=0;
 for(int implicit=0;implicit<2;++implicit)
  for(int split=0;split<2;++split)
  {
   info.m_deformable_implicit=implicit!=0;info.m_splitImpulse=split!=0;
   IntegrationRecoveryProbe probe(c,info);probe.solveConstraint(info);
   const btScalar expected=implicit&&split?btScalar(0):btScalar(2)/info.m_timeStep*(split?btScalar(1):btScalar(1)+info.m_deformable_erp);
   EXPECT_NEAR(double(probe.velocity.x()),double(expected),1e-4);
   if(split){probe.solveSplitImpulse(info);EXPECT_GT(probe.split.x(),btScalar(0));EXPECT_NEAR(double(probe.velocity.x()),double(expected),1e-4);}
  }
}
}  // namespace

// Characterization of the current custom contact response, not desired physics.
TEST(SoftContactDiagnostics, RecordsConstraintBoostWithoutChangingResponse)
{
	btSoftBodyWorldInfo info;
	const btVector3 position(0, 0, 0);
	const btScalar mass = 1;
	btSoftBody a(&info, 1, &position, &mass), b(&info, 1, &position, &mass);
	a.setUserIndex(101);
	b.setUserIndex(102);
	a.m_nodes[0].local_index = b.m_nodes[0].local_index = 0;
	btSoftBody::DeformableNodeNodeContact contact;
	contact.m_node0 = &a.m_nodes[0];
	contact.m_node1 = &b.m_nodes[0];
	contact.m_colObj = &b;
	contact.m_normal = btVector3(1, 0, 0);
	contact.m_offset = 0;
	contact.m_friction = 0;
	contact.m_contact_point_impulse_magnitude = nullptr;
	a.m_nodeNodeContacts.push_back(contact);
	btAlignedObjectArray<btSoftBody*> bodies;
	bodies.push_back(&a);
	bodies.push_back(&b);
	const double expectedEnergy[] = {0.25, 8.125, 12.625, 20.5};
	const double expectedSeparation[] = {0, 4.5, 4.5, 9};
	for (int constrained = 0; constrained < 4; ++constrained)
	{
		btVector3 referenceA, referenceB;
		for (int logging = 0; logging < 2; ++logging)
		{
			a.m_nodes[0].m_v = btVector3(-1, 0, 0);
			b.m_nodes[0].m_v.setZero();
			a.m_nodes[0].m_constrained = (constrained & 1) != 0;
			b.m_nodes[0].m_constrained = (constrained & 2) != 0;
			if (logging)
			{
				btDeformableDiagnostics::StepScope scope(&a, constrained, btScalar(0.0002));
				btDeformableDiagnostics::bodies("before", bodies);
				a.applyRepulsionForce(btScalar(0.0002), true);
				btDeformableDiagnostics::bodies("after", bodies, true);
				EXPECT_EQ(referenceA, a.m_nodes[0].m_v);
				EXPECT_EQ(referenceB, b.m_nodes[0].m_v);
			}
			else
			{
				a.applyRepulsionForce(btScalar(0.0002), true);
				referenceA = a.m_nodes[0].m_v;
				referenceB = b.m_nodes[0].m_v;
			}
			const double energy = 0.5 * double(a.m_nodes[0].m_v.length2() + b.m_nodes[0].m_v.length2());
			EXPECT_DOUBLE_EQ(expectedEnergy[constrained], energy);
			EXPECT_DOUBLE_EQ(expectedSeparation[constrained], double((a.m_nodes[0].m_v - b.m_nodes[0].m_v).x()));
		}
	}
}

int main(int argc, char** argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}


TEST(CoupledContact, ImplicitMomentumConservationFrictionAndSeparation)
{
	for (int separating = 0; separating < 2; ++separating)
	for (int fixed = 0; fixed < 2; ++fixed)
	{
		btSoftBodyWorldInfo info; const btVector3 p(0,0,0); const btScalar mass = 1;
		btSoftBody a(&info,1,&p,&mass), b(&info,1,&p,&mass);
		btAlignedObjectArray<btSoftBody*> bodies; bodies.push_back(&a); bodies.push_back(&b);
		a.m_nodes[0].m_v = a.m_nodes[0].m_vn = btVector3(separating ? 1 : -1, 1, 0);
		b.m_nodes[0].m_v = b.m_nodes[0].m_vn = btVector3(0,0,0);
		b.m_nodes[0].m_frozen = fixed;
		btDeformableBodySolver solver; solver.setImplicit(true); solver.m_useProjection = false;
		solver.setMaxNewtonIterations(8); solver.setNewtonTolerance(btScalar(1e-6)); solver.reinitialize(bodies, btScalar(.01));
		solver.setPreconditioner(btDeformableBackwardEulerObjective::KKT_preconditioner);
		btContactSolverInfo si; si.m_timeStep=btScalar(.01); solver.setConstraints(si); solver.setLagrangeMultiplier();
		solver.setupDeformableSolve(true);
		btSoftBody::DeformableNodeNodeContact c = {}; c.m_node0=&a.m_nodes[0]; c.m_node1=&b.m_nodes[0]; c.m_normal=btVector3(1,0,0);c.m_friction=btScalar(.5);
		btDeformableContactForce contact(btScalar(.01)); contact.add(c); solver.m_objective->m_lf.push_back(&contact);
		btScalar error = 1;
		for(int i=0;i<30;++i) { solver.solveDeformableConstraints(btScalar(.01)); error=contact.updateMultipliers(); if(error<btScalar(1e-5)&&solver.m_lastSolveConverged)break; }
		EXPECT_LT(error,btScalar(1e-4));
		const btVector3 va=a.m_nodes[0].m_v,vb=b.m_nodes[0].m_v;
		if(separating) { EXPECT_NEAR(1.,double(va.x()),1e-5);EXPECT_NEAR(1.,double(va.y()),1e-5); }
		else
		{
			EXPECT_NEAR(0.,double((va-vb).x()),1e-4);
			EXPECT_NEAR(fixed ? .5 : .75,double(va.y()),1e-4);
			EXPECT_LE(double(va.length2()+vb.length2()),2.00001);
		}
		if(fixed) EXPECT_EQ(btVector3(0,0,0),vb);
		else EXPECT_NEAR(0.,double((va+vb-btVector3(separating?1:-1,1,0)).length()),1e-4);
		solver.m_objective->m_lf.clear();
	}
}

TEST(CoupledContact, LineSearchIgnoresConstantEnergyOffsetWhenResidualDecreases)
{
	class OffsetForce : public btDeformableNodalForce
	{
		const double offset;
	public:
		explicit OffsetForce(double value) : btDeformableNodalForce(nullptr,btAlignedObjectArray<int>(),btVector3(0,0,0)),offset(value) {}
		double totalElasticEnergy(btScalar) override { return offset; }
	};
	btVector3 reference;
	for(int shifted=0;shifted<2;++shifted)
	{
		btSoftBodyWorldInfo info;const btVector3 p(0,0,0);const btScalar mass=btScalar(.001/1943),dt=btScalar(.0002);
		btSoftBody body(&info,1,&p,&mass);body.m_nodes[0].m_v=body.m_nodes[0].m_vn=btVector3(-1,btScalar(.2),btScalar(.1));
		btAlignedObjectArray<btSoftBody*> bodies;bodies.push_back(&body);
		btDeformableBodySolver solver;solver.setImplicit(true);solver.setLineSearch(true);solver.setMaxNewtonIterations(5);solver.setNewtonTolerance(btScalar(.01));
		solver.reinitialize(bodies,dt);solver.setPreconditioner(btDeformableBackwardEulerObjective::KKT_preconditioner);solver.setupDeformableSolve(true);
		OffsetForce offset(shifted?1e12:0);solver.m_objective->m_lf.push_back(&offset);
		btSoftBody::DeformableNodeNodeContact c={};c.m_normal=btVector3(1,0,0);c.m_friction=btScalar(.5);
		btSoftBody::ContactNode node={&body.m_nodes[0],btMatrix3x3::getIdentity()};c.m_surfaceNodes.push_back(node);
		btDeformableContactForce contact(dt);ASSERT_TRUE(contact.add(c));solver.m_objective->m_lf.push_back(&contact);
		btScalar error=1;int iteration=0;
		for(;iteration<20;++iteration)
		{
			solver.solveDeformableConstraints(dt);error=contact.updateMultipliers();
			if(solver.m_lastSolveConverged && error<=btScalar(.001))break;
		}
		EXPECT_LT(iteration,20);EXPECT_TRUE(solver.m_lastSolveConverged);EXPECT_LE(error,btScalar(.001));
		EXPECT_LT(body.m_nodes[0].m_v.length(),btScalar(.001));
		if(!shifted)reference=body.m_nodes[0].m_v;
		else EXPECT_LT((body.m_nodes[0].m_v-reference).length(),btScalar(.001));
		solver.m_objective->m_lf.clear();
	}
}

TEST(CoupledContact, LineSearchResolvesSlidingToStickingAtClipMass)
{
	for (int tight = 0; tight < 2; ++tight)
	for (int reducedStep = 0; reducedStep < 2; ++reducedStep)
	for (int speed = 1; speed <= 100; speed *= 10)
	{
		btSoftBodyWorldInfo info;
		const btVector3 p(0,0,0);
		const btScalar mass = btScalar(.001 / 1943);
		btSoftBody a(&info,1,&p,&mass), b(&info,1,&p,&mass);
		b.m_nodes[0].m_frozen = 1;
		a.m_nodes[0].m_v = a.m_nodes[0].m_vn = btVector3(-speed, btScalar(.2)*speed, btScalar(.1)*speed);
		btAlignedObjectArray<btSoftBody*> bodies;
		bodies.push_back(&a); bodies.push_back(&b);
		btDeformableBodySolver solver;
		solver.setImplicit(true); solver.setLineSearch(true); solver.m_useProjection = false;
		solver.setMaxNewtonIterations(5);
		solver.setNewtonTolerance(tight ? btScalar(1e-12) : btScalar(.01));
		const btScalar dt = btScalar(.0002) / (reducedStep ? 256 : 1);
		solver.reinitialize(bodies, dt);
		solver.setPreconditioner(btDeformableBackwardEulerObjective::KKT_preconditioner);
		solver.setupDeformableSolve(true);
		btSoftBody::DeformableNodeNodeContact c = {};
		c.m_node0 = &a.m_nodes[0]; c.m_node1 = &b.m_nodes[0];
		c.m_normal = btVector3(1,0,0); c.m_friction = btScalar(.5);
		btDeformableContactForce contact(dt);
		ASSERT_TRUE(contact.add(c)); solver.m_objective->m_lf.push_back(&contact);
		btScalar error = 0;
		int iteration = 0;
		for (; iteration < 20; ++iteration)
		{
			solver.solveDeformableConstraints(dt);
			error = contact.updateMultipliers();
			if (error < btScalar(.001) && solver.m_lastSolveConverged) break;
		}
		// The incoming tangential impulse lies inside the final friction cone.
		// Full Newton steps overshoot its sticking region and alternate directions.
		EXPECT_LT(iteration, 20);
		EXPECT_LT(error, btScalar(.001));
		EXPECT_LT(a.m_nodes[0].m_v.length(), btScalar(.001));
		EXPECT_TRUE(solver.m_lastSolveConverged);
		EXPECT_EQ(error, btMax(contact.lastNormalError, contact.lastTangentError));
		solver.m_objective->m_lf.clear();
	}
}

TEST(CoupledContact, ContactTangentMatchesFiniteDifferenceAndBlocks)
{
	btSoftBodyWorldInfo info; const btVector3 p(0,0,0); const btScalar mass=1;
	btSoftBody a(&info,1,&p,&mass),b(&info,1,&p,&mass);a.m_nodes[0].index=0;b.m_nodes[0].index=1;
	btSoftBody::DeformableNodeNodeContact c={};c.m_node0=&a.m_nodes[0];c.m_node1=&b.m_nodes[0];c.m_normal=btVector3(1,0,0);c.m_friction=btScalar(.5);
	btDeformableContactForce contact(btScalar(.01));contact.add(c);contact.contacts[0].normalImpulse=1;
	for(int slide=0;slide<2;++slide)
	{
		a.m_nodes[0].m_v=btVector3(-1,slide?2:btScalar(.01),btScalar(.02));b.m_nodes[0].m_v.setZero();
		const btVector3 original=a.m_nodes[0].m_v;
		for(int axis=0;axis<3;++axis)
		{
			Vectors x,ax,plus,minus;x.resize(2,btVector3(0,0,0));ax=x;plus=x;minus=x;x[0][axis]=1;
			contact.addImplicitForceDifferential(btScalar(.01),x,ax);
			btAlignedObjectArray<btMatrix3x3> blocks;blocks.resize(2,btMatrix3x3::getIdentity()*btScalar(0));
			contact.addImplicitForceDifferentialBlocks(btScalar(.01),blocks);
			const btScalar eps=btScalar(.0001);
			a.m_nodes[0].m_v=original+eps*x[0];contact.addScaledForces(btScalar(.01),plus);
			a.m_nodes[0].m_v=original-eps*x[0];contact.addScaledForces(btScalar(.01),minus);
			a.m_nodes[0].m_v=original;
			for(int d=0;d<3;++d)
			{
				EXPECT_NEAR(double(ax[0][d]),double(-(plus[0][d]-minus[0][d])/(2*eps)),.004);
				EXPECT_NEAR(double(ax[0][d]),double(blocks[0][d][axis]),1e-6);
				EXPECT_NEAR(double(ax[0][d]+ax[1][d]),0.,1e-6);
			}
		}
	}
}

namespace
{
class ContactTestWorld : public btDeformableMultiBodyDynamicsWorld
{
public:
	int detections=0;
	bool forceFailure=false;
	bool weighted=false;
	bool staticSurface=false;
	ContactTestWorld(btDispatcher* d,btBroadphaseInterface* b,btDeformableMultiBodyConstraintSolver* c,btCollisionConfiguration* config,btDeformableBodySolver* s)
		:btDeformableMultiBodyDynamicsWorld(d,b,c,config,s){}
	void performDiscreteCollisionDetection() override
	{
		++detections;
		auto& bodies=getSoftBodyArray();if(bodies.size()!=2)return;
		auto& a=bodies[0]->m_nodes[0];auto& b=bodies[1]->m_nodes[0];
		const btScalar gap=(weighted ? btScalar(.25)*a.m_x.x()+btScalar(.75)*bodies[0]->m_nodes[1].m_x.x() : a.m_x.x())-b.m_x.x();
		if(gap<=0 || forceFailure)
		{
			btSoftBody::DeformableNodeNodeContact c={};c.m_node0=&a;c.m_node1=&b;c.m_normal=btVector3(1,0,0);c.m_colObj=bodies[1];c.m_offset=forceFailure?btScalar(-1):gap;
			if (weighted)
			{
				btSoftBody::ContactNode n0={&a,btMatrix3x3::getIdentity()*btScalar(.25)};
				btSoftBody::ContactNode n1={&bodies[0]->m_nodes[1],btMatrix3x3::getIdentity()*btScalar(.75)};
				btSoftBody::ContactNode n2={&b,btMatrix3x3::getIdentity()*btScalar(-1)};
				c.m_surfaceNodes.push_back(n0);c.m_surfaceNodes.push_back(n1);
                if(staticSurface)c.m_node1=nullptr;else c.m_surfaceNodes.push_back(n2);
			}
			bodies[0]->m_nodeNodeContacts.push_back(c);
		}
	}
};
int acceptedCallbacks=0;
void countAccepted(btDynamicsWorld*,btScalar){++acceptedCallbacks;}
}

TEST(CoupledContact, NewContactRefreshAndFailedStepAreTransactional)
{
	for(int fail=0;fail<2;++fail)
	{
		btSoftBodyRigidBodyCollisionConfiguration config;
		btCollisionDispatcher dispatcher(&config);btDbvtBroadphase broadphase;
		btDeformableBodySolver solver;btDeformableMultiBodyConstraintSolver constraints;constraints.setDeformableSolver(&solver);
		ContactTestWorld world(&dispatcher,&broadphase,&constraints,&config,&solver);
		world.setImplicit(true);world.setCoupledContact(true);world.setMaxNewtonIterations(8);world.setNewtonTolerance(btScalar(1e-6));
		world.setGravity(btVector3(0,0,0));world.setInternalTickCallback(countAccepted);
		const btVector3 pa(btScalar(.005),0,0),pb(0,0,0);const btScalar mass=1;
		btSoftBody a(&world.getWorldInfo(),1,&pa,&mass),b(&world.getWorldInfo(),1,&pb,&mass);
		a.m_cfg.drag=b.m_cfg.drag=0;a.m_cfg.collisions=b.m_cfg.collisions=0;
		a.m_nodes[0].m_v=a.m_nodes[0].m_vn=btVector3(-1,0,0);b.m_nodes[0].m_v=b.m_nodes[0].m_vn=btVector3(0,0,0);
		world.addSoftBody(&a);world.addSoftBody(&b);world.forceFailure=fail!=0;acceptedCallbacks=0;
		const int steps=world.stepSimulation(btScalar(.01),0);
		if(fail)
		{
			EXPECT_EQ(0,steps);EXPECT_TRUE(world.hasCoupledStepFailed());EXPECT_EQ(0,acceptedCallbacks);
			EXPECT_EQ(pa,a.m_nodes[0].m_x);EXPECT_EQ(pb,b.m_nodes[0].m_x);
			EXPECT_EQ(btVector3(-1,0,0),a.m_nodes[0].m_v);EXPECT_EQ(btVector3(0,0,0),b.m_nodes[0].m_v);
			EXPECT_EQ(0,world.stepSimulation(btScalar(.01),0));
		}
		else
		{
			EXPECT_EQ(1,steps);EXPECT_FALSE(world.hasCoupledStepFailed());EXPECT_EQ(1,acceptedCallbacks);
			EXPECT_GE(double(a.m_nodes[0].m_x.x()-b.m_nodes[0].m_x.x()),-1e-5);
			EXPECT_NEAR(-1.,double(a.m_nodes[0].m_v.x()+b.m_nodes[0].m_v.x()),1e-4);
			EXPECT_GT(world.detections,3);
		}
		world.removeSoftBody(&a);world.removeSoftBody(&b);
	}
}


TEST(CoupledContact, InvalidPredictorSkipsMultiplierRetriesAndPreservesRollback)
{
	class CountingSolver : public btDeformableBodySolver
	{
	public:
		int calls=0;
		void solveDeformableConstraints(btScalar dt) override
		{ ++calls;btDeformableBodySolver::solveDeformableConstraints(dt); }
	} solver;
	btSoftBodyRigidBodyCollisionConfiguration config;
	btCollisionDispatcher dispatcher(&config);btDbvtBroadphase broadphase;
	btDeformableMultiBodyConstraintSolver constraints;constraints.setDeformableSolver(&solver);
	ContactTestWorld world(&dispatcher,&broadphase,&constraints,&config,&solver);
	world.setImplicit(true);world.setCoupledContact(true);world.setMaxNewtonIterations(5);world.setGravity(btVector3(0,0,0));
	const btVector3 p[]={btVector3(0,0,0),btVector3(1,0,0),btVector3(0,1,0),btVector3(0,0,1)};
	const btScalar masses[]={1,1,1,1};
	btSoftBody body(&world.getWorldInfo(),4,p,masses);
	body.appendTetra(0,1,2,3);body.initializeDmInverse();body.m_tetraScratches.resize(1);body.m_tetraScratchesTn.resize(1);
	body.m_cfg.collisions=0;body.m_cfg.drag=0;
	for(int n=0;n<4;++n)body.m_nodes[n].m_v=body.m_nodes[n].m_vn=btVector3(0,0,0);
	body.m_nodes[3].m_x.setZ(btScalar(-.1));body.updateDeformation();
	btDeformableLinearElasticityForce elastic(1,1,0,0);elastic.addSoftBody(&body);
	world.addSoftBody(&body);world.addForce(&elastic);
	EXPECT_EQ(0,world.stepSimulation(btScalar(.0002),0));
	EXPECT_TRUE(world.hasCoupledStepFailed());EXPECT_TRUE(solver.m_lastSolveInvalidPredictor);
	EXPECT_EQ(9,solver.calls); // dt through dt/256, one solve per rejected trial.
	EXPECT_EQ(btScalar(-.1),body.m_nodes[3].m_x.z());
	world.removeForce(&elastic);world.removeSoftBody(&body);
	// The status is per call and must not poison a later valid solve.
	body.m_nodes[3].m_x=p[3];
	btAlignedObjectArray<btSoftBody*> bodies;bodies.push_back(&body);
	solver.reinitialize(bodies,btScalar(.0002));solver.setupDeformableSolve(true);
	solver.solveDeformableConstraints(btScalar(.0002));
	EXPECT_FALSE(solver.m_lastSolveInvalidPredictor);EXPECT_TRUE(solver.m_lastSolveConverged);
}

TEST(CoupledContact, SuccessfulSubstepsDoNotExhaustRetryBudget)
{
	class SmallStepSolver : public btDeformableBodySolver
	{
	public:
		btScalar largestStep;
		int successfulSolves = 0;
		bool failLate = false;
		void solveDeformableConstraints(btScalar dt) override
		{
			if (dt > largestStep || (failLate && successfulSolves >= 50))
			{ m_lastSolveConverged = false; return; }
			btDeformableBodySolver::solveDeformableConstraints(dt);
			if (m_lastSolveConverged) ++successfulSolves;
		}
	};
	for (int subdivisions : {64, 256})
	for (int failLate = 0; failLate < 2; ++failLate)
	{
		btSoftBodyRigidBodyCollisionConfiguration config;
		btCollisionDispatcher dispatcher(&config); btDbvtBroadphase broadphase;
		SmallStepSolver solver;
		const btScalar dt = btScalar(.01);
		solver.largestStep = dt / subdivisions; solver.failLate = failLate != 0;
		btDeformableMultiBodyConstraintSolver constraints; constraints.setDeformableSolver(&solver);
		ContactTestWorld world(&dispatcher, &broadphase, &constraints, &config, &solver);
		world.setImplicit(true); world.setCoupledContact(true); world.setMaxNewtonIterations(5);
		world.setGravity(btVector3(0,0,0)); world.setInternalTickCallback(countAccepted);
		const btVector3 position(0,0,0), velocity(1,0,0); const btScalar mass = 1;
		btSoftBody body(&world.getWorldInfo(), 1, &position, &mass);
		body.m_cfg.drag = 0; body.m_cfg.collisions = 0;
		body.m_nodes[0].m_v = body.m_nodes[0].m_vn = velocity;
		world.addSoftBody(&body); acceptedCallbacks = 0;
		EXPECT_EQ(failLate ? 0 : 1, world.stepSimulation(dt, 0));
		EXPECT_EQ(failLate != 0, world.hasCoupledStepFailed());
		EXPECT_EQ(failLate ? 0 : 1, acceptedCallbacks);
		EXPECT_NEAR(failLate ? 0. : double(dt), double(body.m_nodes[0].m_x.x()), 1e-8);
		EXPECT_EQ(velocity, body.m_nodes[0].m_v);
		EXPECT_GE(solver.successfulSolves, failLate ? 50 : subdivisions);
		world.removeSoftBody(&body);
	}
}

TEST(ContactRefresh, ReplacesChangedPlaneWithoutMergingDistinctPoints)
{
	btCollisionObject a,b,c;
	btSoftBody::DeformableNodeNodeContact old={};
	old.m_surfaceObjects[0]=&a;old.m_surfaceObjects[1]=&b;
	old.m_surfaceParts[0]=1;old.m_surfaceParts[1]=2;
	old.m_surfaceTriangles[0]=3;old.m_surfaceTriangles[1]=4;
	old.m_normal=btVector3(1,0,0);old.m_offset=.1;
	auto second=old;second.m_offset=.2;
	auto distinct=old;distinct.m_surfaceTriangles[1]=5;
	auto differentBody=old;differentBody.m_surfaceObjects[1]=&c;
	auto differentPart=old;differentPart.m_surfaceParts[1]=6;
	btSoftBody::DeformableNodeNodeContact unknown={};
	auto invalid=old;invalid.m_surfaceInvalid=true;
	btAlignedObjectArray<btSoftBody::DeformableNodeNodeContact> current,updated;
	for(const auto& contact : {old,second,distinct,differentBody,differentPart,unknown,invalid})current.push_back(contact);
	auto fresh=old;fresh.m_normal=btVector3(0,1,0);fresh.m_offset=-.02;
	updated.push_back(fresh);fresh.m_offset=-.03;updated.push_back(fresh);
	EXPECT_EQ(2,btRemoveRefreshedContactPatches(current,updated));
	ASSERT_EQ(5,current.size());EXPECT_EQ(5,current[0].m_surfaceTriangles[1]);
	for(int i=0;i<updated.size();++i)current.push_back(updated[i]);
	ASSERT_EQ(7,current.size());EXPECT_EQ(btScalar(-.02),current[5].m_offset);EXPECT_EQ(btScalar(-.03),current[6].m_offset);
	auto reversed=old;
	for(int i=0;i<2;++i){reversed.m_surfaceObjects[i]=old.m_surfaceObjects[1-i];reversed.m_surfaceParts[i]=old.m_surfaceParts[1-i];reversed.m_surfaceTriangles[i]=old.m_surfaceTriangles[1-i];}
	EXPECT_TRUE(old.sameSurfaceFeature(reversed));EXPECT_FALSE(old.sameSurfaceFeature(unknown));
	updated.clear();reversed.m_surfaceInvalid=true;updated.push_back(reversed);
	EXPECT_EQ(0,btRemoveRefreshedContactPatches(current,updated));
}

TEST(CoupledContact, AdaptiveStartingStepGrowsResetsAndRollsBack)
{
	class LimitedSolver : public btDeformableBodySolver
	{
	public:
		btScalar limit;std::vector<btScalar> calls;
		void solveDeformableConstraints(btScalar dt) override
		{
			calls.push_back(dt);
			if(dt>limit){m_lastSolveConverged=false;return;}
			btDeformableBodySolver::solveDeformableConstraints(dt);
		}
	};
	btSoftBodyRigidBodyCollisionConfiguration config;btCollisionDispatcher dispatcher(&config);btDbvtBroadphase broadphase;
	LimitedSolver solver;btDeformableMultiBodyConstraintSolver constraints;constraints.setDeformableSolver(&solver);
	ContactTestWorld world(&dispatcher,&broadphase,&constraints,&config,&solver);
	world.setImplicit(true);world.setCoupledContact(true);world.setAdaptiveCoupledTimesteps(true);world.setMaxNewtonIterations(5);world.setGravity(btVector3(0,0,0));
	const btScalar dt=.01,mass=1;const btVector3 position(0,0,0),velocity(1,0,0);
	btSoftBody body(&world.getWorldInfo(),1,&position,&mass);body.m_cfg.drag=0;body.m_cfg.collisions=0;
	body.m_nodes[0].m_v=body.m_nodes[0].m_vn=velocity;world.addSoftBody(&body);
	solver.limit=dt/8;
	ASSERT_EQ(1,world.stepSimulation(dt,0));EXPECT_EQ(dt,solver.calls.front());
	solver.calls.clear();ASSERT_EQ(1,world.stepSimulation(dt,0));EXPECT_EQ(dt/4,solver.calls.front());
	EXPECT_NEAR(double(body.m_nodes[0].m_x.x()),double(2*dt),1e-8);
	solver.limit=dt;
	for(int divisor : {4,1,1})
	{
		solver.calls.clear();ASSERT_EQ(1,world.stepSimulation(dt,0));EXPECT_EQ(dt/divisor,solver.calls.front());
		if(divisor==4)
		{
			ASSERT_EQ(3u,solver.calls.size());EXPECT_NEAR(double(dt/4),double(solver.calls[1]),1e-14);EXPECT_NEAR(double(dt/2),double(solver.calls[2]),1e-14);
		}
	}
	// A failed growth probe is rolled back and is not attempted repeatedly.
	solver.limit=dt/8;ASSERT_EQ(1,world.stepSimulation(dt,0));solver.limit=dt/4;solver.calls.clear();
	const btScalar beforeGrowth=body.m_nodes[0].m_x.x();
	ASSERT_EQ(1,world.stepSimulation(dt,0));
	EXPECT_EQ(20,std::count_if(solver.calls.begin(),solver.calls.end(),[&](btScalar h){return btFabs(h-dt/2)<btScalar(1e-14);}));
	EXPECT_EQ(4,std::count_if(solver.calls.begin(),solver.calls.end(),[&](btScalar h){return btFabs(h-dt/4)<btScalar(1e-14);}));
	EXPECT_NEAR(double(beforeGrowth+dt),double(body.m_nodes[0].m_x.x()),1e-8);
	// Establish subdivision again, then verify a changed outer timestep starts fresh.
	solver.limit=dt/8;ASSERT_EQ(1,world.stepSimulation(dt,0));solver.limit=dt;
	solver.calls.clear();ASSERT_EQ(1,world.stepSimulation(dt/2,0));EXPECT_EQ(dt/2,solver.calls.front());
	solver.limit=dt/8;ASSERT_EQ(1,world.stepSimulation(dt,0));
	const btVector3 before=body.m_nodes[0].m_x;solver.limit=0;
	EXPECT_EQ(0,world.stepSimulation(dt,0));EXPECT_TRUE(world.hasCoupledStepFailed());EXPECT_EQ(before,body.m_nodes[0].m_x);
	world.setCoupledContact(true);solver.limit=dt;solver.calls.clear();
	ASSERT_EQ(1,world.stepSimulation(dt,0));EXPECT_EQ(dt,solver.calls.front());
	solver.limit=dt/8;ASSERT_EQ(1,world.stepSimulation(dt,0));world.setAdaptiveCoupledTimesteps(false);
	solver.calls.clear();ASSERT_EQ(1,world.stepSimulation(dt,0));EXPECT_EQ(dt,solver.calls.front());
	world.removeSoftBody(&body);
}

TEST(CoupledContact, StiffTetraAgainstSoftTetraDoesNotInjectEnergy)
{
	btSoftBodyRigidBodyCollisionConfiguration config;btCollisionDispatcher dispatcher(&config);btDbvtBroadphase broadphase;
	btDeformableBodySolver solver;btDeformableMultiBodyConstraintSolver constraints;constraints.setDeformableSolver(&solver);
	ContactTestWorld world(&dispatcher,&broadphase,&constraints,&config,&solver);
	world.setImplicit(true);world.setCoupledContact(true);world.setMaxNewtonIterations(8);world.setNewtonTolerance(btScalar(1e-5));world.setGravity(btVector3(0,0,0));
	const btVector3 pa[]={btVector3(.005,0,0),btVector3(1.005,0,0),btVector3(.005,1,0),btVector3(.005,0,1)};
	const btVector3 pb[]={btVector3(0,0,0),btVector3(1,0,0),btVector3(0,1,0),btVector3(0,0,1)};
	const btScalar masses[]={1,1,1,1};
	btSoftBody a(&world.getWorldInfo(),4,pa,masses),b(&world.getWorldInfo(),4,pb,masses);
	a.appendTetra(0,1,2,3);b.appendTetra(0,1,2,3);a.initializeDmInverse();b.initializeDmInverse();
	a.m_tetraScratches.resize(1);b.m_tetraScratches.resize(1);a.m_tetraScratchesTn.resize(1);b.m_tetraScratchesTn.resize(1);a.updateDeformation();b.updateDeformation();
	a.m_cfg.drag=b.m_cfg.drag=0;a.m_cfg.collisions=b.m_cfg.collisions=0;
	for(int n=0;n<4;++n){a.m_nodes[n].m_v=a.m_nodes[n].m_vn=btVector3(-1,0,0);b.m_nodes[n].m_v=b.m_nodes[n].m_vn=btVector3(0,0,0);}
	btDeformableLinearElasticityForce stiff(1e6,1e6,0,0),soft(1,1,0,0);
	stiff.addSoftBody(&a);soft.addSoftBody(&b);world.addSoftBody(&a);world.addSoftBody(&b);world.addForce(&stiff);world.addForce(&soft);
	for(int step=0;step<20;++step)
	{
		ASSERT_EQ(1,world.stepSimulation(btScalar(.001),0));
		a.updateDeformation();b.updateDeformation();
		double energy=stiff.totalElasticEnergy(0)+soft.totalElasticEnergy(0);
		btVector3 momentum(0,0,0);
		for(int n=0;n<4;++n){energy+=.5*(a.m_nodes[n].m_v.length2()+b.m_nodes[n].m_v.length2());momentum+=a.m_nodes[n].m_v+b.m_nodes[n].m_v;}
		EXPECT_LE(energy,2.001);EXPECT_NEAR(-4.,double(momentum.x()),.002);
		EXPECT_GE(double(a.m_nodes[0].m_x.x()-b.m_nodes[0].m_x.x()),-1e-5);
	}
	world.removeForce(&stiff);world.removeForce(&soft);world.removeSoftBody(&a);world.removeSoftBody(&b);
}


TEST(CoupledContact, ComplianceWeightedAccuracyResolvesSoftFriction)
{
	for (int tight = 0; tight < 2; ++tight)
	for (int trial = 0; trial < 10; ++trial)
	{
		btSoftBodyWorldInfo info;
		const btVector3 p[] = {btVector3(0,0,0),btVector3(1,0,0),btVector3(0,1,0),btVector3(0,0,1)};
		const btScalar m = btScalar(.001/1943), masses[] = {m,m,m,m}, dt = btScalar(.0002);
		btSoftBody body(&info,4,p,masses);
		body.appendTetra(0,1,2,3); body.initializeDmInverse();
		body.m_tetraScratches.resize(1); body.m_tetraScratchesTn.resize(1); body.updateDeformation();
		for (int n=0;n<4;++n) body.m_nodes[n].m_v=body.m_nodes[n].m_vn=btVector3(-1,btScalar(.1)*(trial+1),btScalar(.15)*n);
		btAlignedObjectArray<btSoftBody*> bodies; bodies.push_back(&body);
		btDeformableLinearElasticityForce elastic(1,1,0,btScalar(.01)); elastic.addSoftBody(&body);
		btDeformableBodySolver solver; solver.setImplicit(true); solver.setLineSearch(true); solver.m_useProjection=false;
		solver.setMaxNewtonIterations(5); solver.setNewtonTolerance(tight ? btScalar(1e-8) : btScalar(.01));
		solver.reinitialize(bodies,dt); solver.setPreconditioner(btDeformableBackwardEulerObjective::KKT_preconditioner);
		solver.m_objective->m_lf.push_back(&elastic); solver.setupDeformableSolve(true);
		btDeformableContactForce contact(dt);
		for (int k=0;k<3;++k)
		{
			btSoftBody::DeformableNodeNodeContact c={}; c.m_normal=btVector3(1,0,0); c.m_friction=btScalar(.5);
			btSoftBody::ContactNode a={&body.m_nodes[k],btMatrix3x3::getIdentity()*btScalar(.8)};
			btSoftBody::ContactNode b={&body.m_nodes[3],btMatrix3x3::getIdentity()*btScalar(.2)};
			c.m_surfaceNodes.push_back(a);c.m_surfaceNodes.push_back(b);contact.add(c);
		}
		solver.m_objective->m_KKTPreconditioner->reinitialize(true);
		contact.configureScaling(*solver.m_objective->m_KKTPreconditioner,solver.m_objective->m_projection.m_lagrangeMultipliers,4);
		solver.m_objective->m_lf.push_back(&contact);
		btScalar error=0; int iteration=0;
		for (;iteration<20;++iteration) {solver.solveDeformableConstraints(dt);error=contact.updateMultipliers();if(error<=btScalar(.001)&&solver.m_lastSolveConverged)break;}
		EXPECT_LT(iteration, 20); EXPECT_LE(error, btScalar(.001)); EXPECT_TRUE(solver.m_lastSolveConverged);
		solver.m_objective->m_lf.clear();
	}
}

TEST(CoupledContact, SheetConstraintScalingConvergesAtSceneMasses)
{
	btSoftBodyWorldInfo info;
	const btVector3 p(0,0,0);
	const btScalar clipMass = btScalar(.001/1943), triangleMass = btScalar(.0836468/726), dt = btScalar(.0002);
	btSoftBody clip(&info,1,&p,&clipMass), triangle(&info,1,&p,&triangleMass);
	btAlignedObjectArray<btSoftBody*> bodies; bodies.push_back(&clip); bodies.push_back(&triangle);
	clip.m_nodes[0].m_v = clip.m_nodes[0].m_vn = btVector3(0,0,0);
	triangle.m_nodes[0].m_v = triangle.m_nodes[0].m_vn = btVector3(2,0,0);
	btDeformableBodySolver solver;
	solver.setImplicit(true); solver.m_useProjection = false;
	solver.setMaxNewtonIterations(5); solver.setNewtonTolerance(btScalar(.01)); solver.reinitialize(bodies,dt);
	solver.setPreconditioner(btDeformableBackwardEulerObjective::KKT_preconditioner);
	LagrangeMultiplier lm = {};
	lm.m_num_nodes=1; lm.m_indices[0]=0; lm.m_weights[0]=1; lm.m_num_constraints=1; lm.m_dirs[0]=btVector3(1,0,0);
	solver.m_objective->m_projection.m_lagrangeMultipliers.push_back(lm);
	solver.setupDeformableSolve(true);
	btSoftBody::DeformableNodeNodeContact c = {};
	c.m_node0=&clip.m_nodes[0]; c.m_node1=&triangle.m_nodes[0]; c.m_normal=btVector3(1,0,0);
	btDeformableContactForce contact(dt); ASSERT_TRUE(contact.add(c));
	auto& objective = *solver.m_objective;
	objective.m_KKTPreconditioner->reinitialize(true);
	contact.configureScaling(*objective.m_KKTPreconditioner,objective.m_projection.m_lagrangeMultipliers,objective.m_nodes.size());
	EXPECT_NEAR(double(10*triangleMass),double(contact.contacts[0].rho),1e-8);
	EXPECT_LT(contact.contacts[0].tangentRho,contact.contacts[0].rho/btScalar(100));
	objective.m_lf.push_back(&contact);
	btScalar error=1;
	for(int iteration=0;iteration<20;++iteration)
	{
		solver.solveDeformableConstraints(dt); error=contact.updateMultipliers();
		if(error<=btScalar(.001) && solver.m_lastSolveConverged) break;
	}
	EXPECT_LE(error,btScalar(.001));
	EXPECT_NEAR(0.,double(clip.m_nodes[0].m_v.x()),.001);
	EXPECT_NEAR(0.,double(triangle.m_nodes[0].m_v.x()),.001);
	objective.m_lf.clear();
}

TEST(CoupledContact, WeightedSurfaceStopsWhenNearestNodeWouldNot)
{
	btSoftBodyWorldInfo info;
	const btVector3 p[2]={btVector3(0,0,0),btVector3(0,1,0)}; const btScalar masses[2]={1,1};
	btSoftBody a(&info,2,p,masses),b(&info,1,p,masses);
	btAlignedObjectArray<btSoftBody*> bodies; bodies.push_back(&a); bodies.push_back(&b);
	for(int n=0;n<2;++n) a.m_nodes[n].m_v=a.m_nodes[n].m_vn=btVector3(-1,0,0);
	b.m_nodes[0].m_v=b.m_nodes[0].m_vn=btVector3(0,0,0);
	btDeformableBodySolver solver;
	solver.setImplicit(true); solver.m_useProjection=false; solver.setMaxNewtonIterations(5); solver.setNewtonTolerance(btScalar(1e-8));
	solver.reinitialize(bodies,btScalar(.0002)); solver.setPreconditioner(btDeformableBackwardEulerObjective::KKT_preconditioner); solver.setupDeformableSolve(true);
	btSoftBody::DeformableNodeNodeContact c={}; c.m_node0=&a.m_nodes[0]; c.m_node1=&b.m_nodes[0]; c.m_normal=btVector3(1,0,0);
	const btScalar weights[]={btScalar(.25),btScalar(.75),-1};
	for(int n=0;n<3;++n)
	{
		btSoftBody::ContactNode node={n<2?&a.m_nodes[n]:&b.m_nodes[0],btMatrix3x3::getIdentity()*weights[n]};
		c.m_surfaceNodes.push_back(node);
	}
	btDeformableContactForce contact(btScalar(.0002)); ASSERT_TRUE(contact.add(c)); solver.m_objective->m_lf.push_back(&contact);
	btScalar error=1;
	for(int iteration=0;iteration<20;++iteration){solver.solveDeformableConstraints(btScalar(.0002));error=contact.updateMultipliers();}
	const btScalar surfaceVelocity=btScalar(.25)*a.m_nodes[0].m_v.x()+btScalar(.75)*a.m_nodes[1].m_v.x()-b.m_nodes[0].m_v.x();
	EXPECT_LT(error,btScalar(.001)); EXPECT_NEAR(0.,double(surfaceVelocity),.001);
	EXPECT_NEAR(-2.,double(a.m_nodes[0].m_v.x()+a.m_nodes[1].m_v.x()+b.m_nodes[0].m_v.x()),1e-6);
	// The contact impulse must follow the surface weights, including its moment.
	EXPECT_NEAR(3.,double((a.m_nodes[1].m_v.x()+1)/(a.m_nodes[0].m_v.x()+1)),1e-5);
	solver.m_objective->m_lf.clear();
}

TEST(CoupledContact, ScalingProjectsWeightedFaceAndRedundantConstraints)
{
	btSoftBodyWorldInfo info;const btVector3 p[3]={btVector3(0,0,0),btVector3(1,0,0),btVector3(0,1,0)};const btScalar masses[3]={1,1,1};
	btSoftBody body(&info,3,p,masses);btAlignedObjectArray<btSoftBody*> bodies;bodies.push_back(&body);
	btDeformableBodySolver solver;solver.setImplicit(true);solver.reinitialize(bodies,btScalar(.01));
	LagrangeMultiplier row={};row.m_num_nodes=2;row.m_indices[0]=0;row.m_indices[1]=1;row.m_weights[0]=row.m_weights[1]=btScalar(.5);row.m_num_constraints=1;row.m_dirs[0]=btVector3(1,0,0);
	auto& objective=*solver.m_objective;objective.m_projection.m_lagrangeMultipliers.push_back(row);objective.m_projection.m_lagrangeMultipliers.push_back(row);
	btSoftBody::DeformableNodeNodeContact c={};c.m_node0=&body.m_nodes[0];c.m_node1=&body.m_nodes[2];c.m_normal=btVector3(1,0,0);
	btDeformableContactForce contact(btScalar(.01));ASSERT_TRUE(contact.add(c));objective.m_KKTPreconditioner->reinitialize(true);
	contact.configureScaling(*objective.m_KKTPreconditioner,objective.m_projection.m_lagrangeMultipliers,objective.m_nodes.size());
	EXPECT_NEAR(10./1.5,double(contact.contacts[0].rho),1e-6);
	EXPECT_NEAR(5.,double(contact.contacts[0].tangentRho),1e-6);
}

namespace
{
class MappedContactBody : public btSoftBody
{
public:
	std::vector<btVertexToTetraMapping> mapping;
	MappedContactBody(btSoftBodyWorldInfo* info,const btVector3* p,const btScalar* m):btSoftBody(info,4,p,m){}
	const std::vector<btVertexToTetraMapping>* getCollisionShapeVertexToSimTetra() const override {return &mapping;}
};
class MappedContactManager : public btGImpactMeshShapePart::TrimeshPrimitiveManager
{
public:
	MappedContactBody* body;
	explicit MappedContactManager(MappedContactBody* b):body(b){}
	void get_vertex(unsigned int index,btVector3& v,bool) const override
	{
		v.setZero();const auto& entry=body->mapping[index];
		for(int n=0;n<4;++n)v+=entry.baryCoordInTetra[n]*body->m_tetras[entry.vertexToTetra].m_n[n]->m_x;
		v*=m_scale;
	}
};
}

TEST(CoupledContact, GImpactSurfaceJacobianMatchesActualMappedGeometry)
{
	btSoftBodyWorldInfo info;const btVector3 positions[]={btVector3(0,0,0),btVector3(2,0,0),btVector3(0,2,0),btVector3(0,0,2)};const btScalar masses[]={1,1,1,1};
	int indices[]={0,1,2};btScalar vertices[]={0,0,0,1,0,0,0,1,0};
	btTriangleIndexVertexArray mesh(1,indices,3*sizeof(int),3,vertices,3*sizeof(btScalar));
	MappedContactBody body(&info,positions,masses);body.appendTetra(0,1,2,3);body.mapping.resize(3);
	const btVector4 weights[]={btVector4(.5,.5,0,0),btVector4(0,.25,.75,0),btVector4(.2,0,.2,.6)};
	for(int v=0;v<3;++v){body.mapping[v].vertexToTetra=0;body.mapping[v].baryCoordInTetra=weights[v];}
	auto* manager=new MappedContactManager(&body);auto* shape=new btGImpactMeshShape(&mesh,manager);
	delete body.getCollisionShape();body.setCollisionShape(shape);
	shape->setLocalScaling(btVector3(2,3,4));shape->updateBound();
	btTransform transform;transform.setIdentity();transform.setRotation(btQuaternion(btVector3(1,2,3).normalized(),btScalar(.3)));transform.setOrigin(btVector3(4,5,6));body.setWorldTransform(transform);
	const btVector3 bary(btScalar(.2),btScalar(.3),btScalar(.5));
	auto point=[&](){btVector3 result(0,0,0);for(int v=0;v<3;++v){btVector3 vertex;manager->get_vertex(v,vertex,false);result+=bary[v]*vertex;}return transform*result;};
	const btVector3 originalPoint=point();btAlignedObjectArray<btSoftBody::ContactNode> stencil;
	ASSERT_TRUE(body.appendSurfaceContactNodes(0,0,originalPoint,1,stencil));ASSERT_EQ(4,stencil.size());
	for(int n=0;n<stencil.size();++n)for(int d=0;d<3;++d)
	{
		const btScalar eps=btScalar(1e-5);const btVector3 original=stencil[n].node->m_x;
		stencil[n].node->m_x[d]+=eps;const btVector3 plus=point();stencil[n].node->m_x=original;stencil[n].node->m_x[d]-=eps;const btVector3 minus=point();stencil[n].node->m_x=original;
		btVector3 direction(0,0,0);direction[d]=1;
		EXPECT_NEAR(0.,double(((plus-minus)/(2*eps)-stencil[n].jacobian*direction).length()),1e-6);
	}
	btAlignedObjectArray<btSoftBody::ContactNode> invalid;
	EXPECT_FALSE(body.appendSurfaceContactNodes(0,4,originalPoint,1,invalid));
	btSoftBody::DeformableNodeNodeContact c={};c.m_surfaceInvalid=true;
	btDeformableContactForce force(btScalar(.01));EXPECT_FALSE(force.add(c));
}


TEST(CoupledContact, WeightedContactDifferentialAndEnergyAreConsistent)
{
	btSoftBodyWorldInfo info;const btVector3 p[3]={btVector3(0,0,0),btVector3(0,1,0),btVector3(0,2,0)};const btScalar masses[3]={1,1,1};
	btSoftBody body(&info,3,p,masses);btSoftBody::DeformableNodeNodeContact c={};c.m_normal=btVector3(1,0,0);c.m_friction=btScalar(.5);
	const btScalar weights[]={btScalar(.25),btScalar(.75),-1};
	for(int n=0;n<3;++n){body.m_nodes[n].index=n;btSoftBody::ContactNode entry={&body.m_nodes[n],btMatrix3x3::getIdentity()*weights[n]};c.m_surfaceNodes.push_back(entry);}
	btDeformableContactForce contact(btScalar(.01));ASSERT_TRUE(contact.add(c));ASSERT_TRUE(contact.add(c));EXPECT_EQ(1,contact.contacts.size());
	contact.contacts[0].normalImpulse=1;contact.contacts[0].rho=9;contact.contacts[0].tangentRho=3;
	for(int slide=0;slide<2;++slide)
	{
		for(int n=0;n<3;++n)body.m_nodes[n].m_v=btVector3(n==2?0:-1, n==2?0:(slide?2:btScalar(.01)),0);
		for(int node=0;node<3;++node)for(int axis=0;axis<3;++axis)
		{
			Vectors direction,product,plus,minus,impulse;direction.resize(3,btVector3(0,0,0));product=plus=minus=impulse=direction;direction[node][axis]=1;
			contact.addImplicitForceDifferential(btScalar(.01),direction,product);contact.addScaledForces(btScalar(.01),impulse);
			btAlignedObjectArray<btMatrix3x3> blocks;blocks.resize(3,btMatrix3x3::getIdentity()*btScalar(0));contact.addImplicitForceDifferentialBlocks(btScalar(.01),blocks);
			const btScalar eps=btScalar(1e-4);const btVector3 original=body.m_nodes[node].m_v;
			body.m_nodes[node].m_v=original+eps*direction[node];contact.addScaledForces(btScalar(.01),plus);const double energyPlus=contact.totalEnergy(btScalar(.01));
			body.m_nodes[node].m_v=original-eps*direction[node];contact.addScaledForces(btScalar(.01),minus);const double energyMinus=contact.totalEnergy(btScalar(.01));body.m_nodes[node].m_v=original;
			for(int n=0;n<3;++n)for(int d=0;d<3;++d)EXPECT_NEAR(double(product[n][d]),double(-(plus[n][d]-minus[n][d])/(2*eps)),.004);
			for(int d=0;d<3;++d)EXPECT_NEAR(double(blocks[node][d][axis]),double(product[node][d]),1e-6);
			EXPECT_NEAR(double(-impulse[node][axis]),(energyPlus-energyMinus)/(2*eps),.004);
			EXPECT_NEAR(0.,double((impulse[0]+impulse[1]+impulse[2]).length()),1e-6);
		}
	}
	c.m_surfaceNodes[0].jacobian=btMatrix3x3::getIdentity()*btScalar(.3);c.m_surfaceNodes[1].jacobian=btMatrix3x3::getIdentity()*btScalar(.7);
	ASSERT_TRUE(contact.add(c));EXPECT_EQ(2,contact.contacts.size());
}


TEST(CoupledContact, TrialRefreshUsesWeightedSurfaceDisplacement)
{
	btSoftBodyRigidBodyCollisionConfiguration config;btCollisionDispatcher dispatcher(&config);btDbvtBroadphase broadphase;
	btDeformableBodySolver solver;btDeformableMultiBodyConstraintSolver constraints;constraints.setDeformableSolver(&solver);
	ContactTestWorld world(&dispatcher,&broadphase,&constraints,&config,&solver);world.weighted=true;
	world.setImplicit(true);world.setCoupledContact(true);world.setMaxNewtonIterations(8);world.setNewtonTolerance(btScalar(1e-6));world.setGravity(btVector3(0,0,0));
	const btVector3 pa[]={btVector3(.005,0,0),btVector3(.015,1,0)},pb(0,0,0);const btScalar masses[]={1,1};
	btSoftBody a(&world.getWorldInfo(),2,pa,masses),b(&world.getWorldInfo(),1,&pb,masses);
	a.m_cfg.drag=b.m_cfg.drag=0;a.m_cfg.collisions=b.m_cfg.collisions=0;
	for(int n=0;n<2;++n)a.m_nodes[n].m_v=a.m_nodes[n].m_vn=btVector3(-1,0,0);
	b.m_nodes[0].m_v=b.m_nodes[0].m_vn=btVector3(0,0,0);
	world.addSoftBody(&a);world.addSoftBody(&b);
	EXPECT_EQ(1,world.stepSimulation(btScalar(.02),0));EXPECT_FALSE(world.hasCoupledStepFailed());
	const btScalar gap=btScalar(.25)*a.m_nodes[0].m_x.x()+btScalar(.75)*a.m_nodes[1].m_x.x()-b.m_nodes[0].m_x.x();
	EXPECT_GE(double(gap),-1e-5);
	EXPECT_NEAR(-2.,double(a.m_nodes[0].m_v.x()+a.m_nodes[1].m_v.x()+b.m_nodes[0].m_v.x()),1e-5);
	EXPECT_GT(world.detections,3);
	world.removeSoftBody(&a);world.removeSoftBody(&b);
}

TEST(CoupledContact, GeometricGapIgnoresLegacyPushbackAndRejectsRecoveryDistance)
{
	btSoftBodyWorldInfo info;
	const btVector3 positions[]={btVector3(0,0,0),btVector3(2,0,0),btVector3(0,2,0),btVector3(0,0,2)};
	const btScalar masses[]={1,1,1,1};
	int indices[]={0,1,2};btScalar vertices[]={0,0,0,2,0,0,0,2,0};
	btTriangleIndexVertexArray mesh(1,indices,3*sizeof(int),3,vertices,3*sizeof(btScalar));
	MappedContactBody a(&info,positions,masses),b(&info,positions,masses);
	for(auto* body : {&a,&b})
	{
		body->appendTetra(0,1,2,3);body->mapping.resize(3);
		for(int v=0;v<3;++v){body->mapping[v].vertexToTetra=0;body->mapping[v].baryCoordInTetra=btVector4(0,0,0,0);body->mapping[v].baryCoordInTetra[v]=1;}
		auto* shape=new btGImpactMeshShape(&mesh,new MappedContactManager(body));
		delete body->getCollisionShape();body->setCollisionShape(shape);shape->setMargin(btScalar(.1));shape->updateBound();body->m_useSurfaceContact=true;
	}
	btTransform transform;transform.setIdentity();transform.setOrigin(btVector3(0,0,btScalar(.15)));a.setWorldTransform(transform);
	const btVector3 point(btScalar(.4),btScalar(.6),btScalar(.15));
	auto contact=[&](btScalar legacy,btScalar separation,bool penetrating)
	{
		a.m_nodeNodeContacts.clear();
		a.skinSoftSoftCollisionHandler(&b,0,0,0,0,point,btVector3(0,0,-1),legacy,penetrating,nullptr,separation);
		return a.m_nodeNodeContacts[0];
	};
	btPrimitiveTriangle faceA,faceB;
	for(int v=0;v<3;++v){faceB.m_vertices[v]=positions[v];faceA.m_vertices[v]=transform*positions[v];}
	faceA.m_margin=faceB.m_margin=btScalar(.1);faceA.buildTriPlane();faceB.buildTriPlane();
	GIM_TRIANGLE_CONTACT generated;
	btTransform identity;identity.setIdentity();
	ASSERT_TRUE(faceA.find_triangle_collision_alt_method_outer(faceB,generated,btScalar(.2),identity,identity,faceA,faceB,false,false));
	EXPECT_NEAR(.15,double(generated.m_unmodified_depth),1e-7);
	EXPECT_GT(btFabs(generated.m_unmodified_depth-generated.m_penetration_depth),btScalar(.1));
	a.skinSoftSoftCollisionHandler(&b,0,0,0,0,generated.m_points[0],-generated.m_separating_normal,
		-generated.m_penetration_depth,false,nullptr,generated.m_unmodified_depth);
	ASSERT_FALSE(a.m_nodeNodeContacts[0].m_surfaceInvalid);
	EXPECT_NEAR(-.045,double(a.m_nodeNodeContacts[0].m_offset),1e-7);
	const btVector3 pointB=generated.m_points[0]-generated.m_separating_normal*generated.m_unmodified_depth;
	b.skinSoftSoftCollisionHandler(&a,0,0,0,0,pointB,generated.m_separating_normal,
		-generated.m_penetration_depth,false,nullptr,generated.m_unmodified_depth);
	ASSERT_FALSE(b.m_nodeNodeContacts[0].m_surfaceInvalid);
	EXPECT_EQ(a.m_nodeNodeContacts[0].m_offset,b.m_nodeNodeContacts[0].m_offset);
	const auto c=contact(btScalar(-.009),btScalar(.15),false);
	ASSERT_FALSE(c.m_surfaceInvalid);EXPECT_NEAR(-.045,double(c.m_offset),1e-7);
	const auto rescaled=contact(btScalar(-9),btScalar(.15),false);
	EXPECT_EQ(c.m_offset,rescaled.m_offset);
	btVector3 relative(0,0,0);
	for(int i=0;i<c.m_surfaceNodes.size();++i)
	{
		const auto& node=c.m_surfaceNodes[i];bool onA=false;
		for(int n=0;n<a.m_nodes.size();++n)onA|=node.node==&a.m_nodes[n];
		if(onA)relative+=node.jacobian*btVector3(0,0,1);
	}
	const btScalar eps=btScalar(1e-4);
	const auto shifted=contact(btScalar(-.009),btScalar(.15)+eps,false);
	EXPECT_NEAR(double(c.m_normal.dot(relative)),double((shifted.m_offset-c.m_offset)/eps),1e-6);
	EXPECT_TRUE(contact(btScalar(-.1),btScalar(.15),true).m_surfaceInvalid);
	EXPECT_TRUE(contact(btScalar(-.1),btScalar(-1),false).m_surfaceInvalid);
    btGImpactMeshShape rigidShape(&mesh);rigidShape.setMargin(btScalar(.1));rigidShape.updateBound();
    btCollisionObject rigid;rigid.setCollisionShape(&rigidShape);rigid.setCollisionFlags(btCollisionObject::CF_STATIC_OBJECT);
    btCollisionObjectWrapper wrapper(nullptr,&rigidShape,&rigid,identity,-1,-1);
    a.m_nodeNodeContacts.clear();
    a.skinSoftStaticCollisionHandler(&wrapper,0,0,0,0,point,btVector3(0,0,1),btScalar(.15),false,nullptr);
    ASSERT_EQ(4,a.m_nodeNodeContacts.size());const auto& fixedContact=a.m_nodeNodeContacts[0];
    EXPECT_EQ(&a,fixedContact.m_surfaceObjects[0]);EXPECT_EQ(&rigid,fixedContact.m_surfaceObjects[1]);
    EXPECT_EQ(0,fixedContact.m_surfaceTriangles[0]);EXPECT_EQ(0,fixedContact.m_surfaceTriangles[1]);
    ASSERT_FALSE(fixedContact.m_surfaceInvalid);EXPECT_EQ(nullptr,fixedContact.m_node1);EXPECT_EQ(0,a.m_nodeRigidContacts.size());
    EXPECT_NEAR(-.045,double(fixedContact.m_offset),1e-7);
    btVector3 translation(0,0,0);
    for(int n=0;n<fixedContact.m_surfaceNodes.size();++n)translation+=fixedContact.m_surfaceNodes[n].jacobian*btVector3(1,2,3);
    EXPECT_LT((translation-btVector3(1,2,3)).length(),btScalar(1e-10));
    btDeformableContactForce fixedForce(btScalar(.0002));EXPECT_TRUE(fixedForce.add(fixedContact));
    a.skinSoftStaticCollisionHandler(&wrapper,0,0,0,0,point,btVector3(0,0,1),btScalar(.15),true,nullptr);
    EXPECT_TRUE(a.m_nodeNodeContacts[4].m_surfaceInvalid);
    rigid.setCollisionFlags(btCollisionObject::CF_KINEMATIC_OBJECT);
    a.skinSoftStaticCollisionHandler(&wrapper,0,0,0,0,point,btVector3(0,0,1),btScalar(.15),false,nullptr);
    EXPECT_TRUE(a.m_nodeNodeContacts[5].m_surfaceInvalid);
    rigid.setCollisionFlags(btCollisionObject::CF_STATIC_OBJECT);
    a.skinSoftStaticCollisionHandler(&wrapper,0,0,0,0,point,btVector3(0,0,0),btScalar(.15),false,nullptr);
    EXPECT_TRUE(a.m_nodeNodeContacts[6].m_surfaceInvalid);

}


TEST_F(DeformableBlockPreconditioner, VolumeBarrierEnergyForceAndTangentAgree)
{
    btDeformableVolumeBarrierForce barrier;
    btDeformableVolumeBarrierForce::Material material={body,btScalar(2)};barrier.materials.push_back(material);
    Vectors original,direction,plus,minus,analytic,load;original.resize(4);direction.resize(4);plus.resize(4);minus.resize(4);analytic.resize(4);load.resize(4);
    for(int n=0;n<4;++n)
    {
        original[n]=body->m_nodes[n].m_x;original[n].setZ(original[n].z()*btScalar(.2));body->m_nodes[n].m_q=original[n];
        direction[n]=btVector3(btScalar(.13*(n+1)),btScalar(.1-.07*n),btScalar(.03*n));analytic[n].setZero();load[n].setZero();
    }
    barrier.addScaledForces(1,load);barrier.addScaledElasticForceDifferential(1,direction,analytic);
    btVector3 sum(0,0,0),torque(0,0,0);btScalar work=0;
    for(int n=0;n<4;++n){sum+=load[n];torque+=original[n].cross(load[n]);work+=load[n].dot(direction[n]);}
    EXPECT_LT(sum.length(),1e-10);EXPECT_LT(torque.length(),1e-10);
    const btScalar eps=btScalar(1e-6);
    for(int n=0;n<4;++n){body->m_nodes[n].m_q=original[n]+eps*direction[n];plus[n].setZero();}
    const double ep=barrier.totalEnergy(0);barrier.addScaledForces(1,plus);
    for(int n=0;n<4;++n){body->m_nodes[n].m_q=original[n]-eps*direction[n];minus[n].setZero();}
    const double em=barrier.totalEnergy(0);barrier.addScaledForces(1,minus);
    EXPECT_NEAR(double(-work),(ep-em)/(2*eps),1e-7);
    for(int n=0;n<4;++n){EXPECT_LT(((plus[n]-minus[n])/(2*eps)-analytic[n]).length(),1e-6);body->m_nodes[n].m_q=original[n];}
    btAlignedObjectArray<btMatrix3x3> blocks;blocks.resize(4);
    for(int n=0;n<4;++n)blocks[n]=btMatrix3x3::getIdentity()*btScalar(0);
    ASSERT_TRUE(barrier.addImplicitForceDifferentialBlocks(btScalar(.01),blocks));
    for(int n=0;n<4;++n)for(int axis=0;axis<3;++axis)
    {
        for(int k=0;k<4;++k){direction[k].setZero();analytic[k].setZero();}direction[n][axis]=1;
        barrier.addImplicitForceDifferential(btScalar(.01),direction,analytic);
        EXPECT_LT((analytic[n]-blocks[n].getColumn(axis)).length(),1e-10);
    }
    btScalar energy,first,second;
    for(btScalar j : {btScalar(.5),btScalar(1),btScalar(2)})
    {barrier.density(j,energy,first,second);EXPECT_EQ(0,energy);EXPECT_EQ(0,first);EXPECT_EQ(0,second);}
    barrier.density(btScalar(.050001),energy,first,second);EXPECT_GT(energy,1);EXPECT_LT(first,-100000);EXPECT_GT(second,1e10);
    barrier.density(btScalar(.05),energy,first,second);EXPECT_FALSE(std::isfinite(double(energy)));
}

TEST_F(DeformableBlockPreconditioner, VolumeStepLimitProtectsEntirePath)
{
    btDeformableVolumeBarrierForce barrier;btDeformableVolumeBarrierForce::Material material={body,1};barrier.materials.push_back(material);
    Vectors delta;delta.resize(4);
    for(int n=0;n<4;++n){body->m_nodes[n].m_q=body->m_nodes[n].m_x;delta[n]=btVector3(-2*body->m_nodes[n].m_x.x(),-2*body->m_nodes[n].m_x.y(),0);}
    // Both endpoints have positive determinant, but the unrestricted path crosses zero.
    const btScalar scale=barrier.safeStep(1,delta);EXPECT_GT(scale,0);EXPECT_LT(scale,btScalar(.5));
    for(int sample=0;sample<=100;++sample)
    {
        for(int n=0;n<4;++n)body->m_nodes[n].m_q=body->m_nodes[n].m_x+delta[n]*scale*btScalar(sample/100.0);
        EXPECT_TRUE(barrier.admissible());
    }
    body->m_nodes[3].m_q.setZ(0);EXPECT_FALSE(barrier.admissible());EXPECT_EQ(0,barrier.safeStep(1,delta));
}

TEST(VolumeBarrier, CompressionSolveRemainsValidWithSoftMaterial)
{
    for(bool lineSearch : {false,true})
    {
    btSoftBodyWorldInfo info;const btVector3 p[]={btVector3(0,0,0),btVector3(1,0,0),btVector3(0,1,0),btVector3(0,0,1)};
    const btScalar masses[]={1,1,1,btScalar(.001/1943)};btSoftBody body(&info,4,p,masses);body.appendTetra(0,1,2,3);body.initializeDmInverse();
    body.m_tetraScratches.resize(1);body.m_tetraScratchesTn.resize(1);body.updateDeformation();
    for(int k=0;k<3;++k){body.m_nodes[k].m_v=body.m_nodes[k].m_vn=btVector3(0,0,0);}
    body.m_nodes[3].m_x.setZ(btScalar(.1));body.m_nodes[3].m_v=body.m_nodes[3].m_vn=btVector3(0,0,-100);
    btAlignedObjectArray<btSoftBody*> bodies;bodies.push_back(&body);
    btDeformableLinearElasticityForce elastic;elastic.setYoungsModulus(1);elastic.setPoissonRatio(btScalar(.4));elastic.setDamping(0,0);elastic.addSoftBody(&body);
    btDeformableVolumeBarrierForce barrier;btDeformableVolumeBarrierForce::Material material={&body,btScalar(1/.6)};barrier.materials.push_back(material);
    btDeformableBodySolver solver;solver.setImplicit(true);solver.setLineSearch(lineSearch);solver.m_useProjection=false;solver.setMaxNewtonIterations(30);solver.setNewtonTolerance(btScalar(1e-4));
    btAlignedObjectArray<int> loaded;loaded.push_back(3);btDeformableNodalForce pressure(&body,loaded,btVector3(0,0,-5));
    solver.m_objective->m_lf.push_back(&pressure);
    solver.m_objective->m_lf.push_back(&elastic);solver.m_objective->m_lf.push_back(&barrier);
    const btScalar dt=btScalar(.0002);solver.reinitialize(bodies,dt);solver.setPreconditioner(btDeformableBackwardEulerObjective::KKT_preconditioner);
    for(int n=0;n<3;++n){LagrangeMultiplier lm={};lm.m_num_nodes=1;lm.m_indices[0]=n;lm.m_weights[0]=1;lm.m_num_constraints=3;for(int d=0;d<3;++d){lm.m_dirs[d].setZero();lm.m_dirs[d][d]=1;}solver.m_objective->m_projection.m_lagrangeMultipliers.push_back(lm);}
    solver.setupDeformableSolve(true);
    btDeformableDiagnostics::StepScope diagnostic(&solver,0,dt);
    solver.solveDeformableConstraints(dt);solver.updateState();
    EXPECT_TRUE(solver.m_lastSolveConverged);EXPECT_TRUE(barrier.admissible());EXPECT_GT(body.m_nodes[3].m_q.z(),btScalar(.05));EXPECT_LT(body.m_nodes[3].m_q.z(),btScalar(.08));
    Vectors balance;balance.resize(4,btVector3(0,0,0));elastic.addScaledForces(dt,balance);barrier.addScaledForces(dt,balance);pressure.addScaledForces(dt,balance);
    balance[3]-=(body.m_nodes[3].m_v-btVector3(0,0,-100))/body.m_nodes[3].m_im;
    EXPECT_LT(balance[3].length(),btScalar(1e-6));
    for(int n=0;n<3;++n)EXPECT_LT(body.m_nodes[n].m_v.length(),btScalar(1e-6));
    solver.m_objective->m_lf.clear();
    }
}


TEST_F(DeformableBlockPreconditioner, VolumeBarrierPreservesTranslationCorrection)
{
    btDeformableVolumeBarrierForce barrier;
    btDeformableVolumeBarrierForce::Material material={body,btScalar(2)};barrier.materials.push_back(material);
    Vectors backup,translation,product,result;
    backup.resize(4,btVector3(0,0,0));translation.resize(4,btVector3(2,-1,3));product=backup;result=backup;
    btDeformableBackwardEulerObjective objective(bodies,backup);
    objective.updateId();objective.setDt(btScalar(.02));objective.setImplicit(true);
    objective.m_lf.push_back(&force);objective.m_lf.push_back(&barrier);
    for(btScalar compression : {btScalar(1),btScalar(.2)})
    {
        for(int n=0;n<4;++n){body->m_nodes[n].m_q=body->m_nodes[n].m_x;body->m_nodes[n].m_q[2]*=compression;}
        body->updateDeformation();
        barrier.prepareImplicitForceDifferential(btScalar(.02));
        objective.m_preconditioner->reinitialize(true);
        ASSERT_TRUE(objective.setupTranslationCorrection());
        objective.multiply(translation,product);objective.precondition(product,result);
        for(int n=0;n<4;++n)EXPECT_LT((result[n]-translation[n]).length(),btScalar(2e-4));
        barrier.finishImplicitForceDifferential();
    }
}


TEST(CoupledContact, StaticSurfaceStopsWhileNearestNodeIsAlreadyConstrained)
{
    btSoftBodyWorldInfo info;const btVector3 p[]={btVector3(0,0,0),btVector3(0,1,0)};const btScalar masses[]={1,1};
    btSoftBody body(&info,2,p,masses);btAlignedObjectArray<btSoftBody*> bodies;bodies.push_back(&body);
    body.m_nodes[0].m_v=body.m_nodes[0].m_vn=btVector3(0,0,0);
    body.m_nodes[1].m_v=body.m_nodes[1].m_vn=btVector3(-1,0,0);
    const btScalar dt=btScalar(.0002);
    btDeformableBodySolver solver;solver.setImplicit(true);solver.setMaxNewtonIterations(5);solver.setNewtonTolerance(btScalar(1e-8));
    solver.reinitialize(bodies,dt);solver.setPreconditioner(btDeformableBackwardEulerObjective::KKT_preconditioner);
    LagrangeMultiplier lm={};lm.m_num_nodes=1;lm.m_indices[0]=0;lm.m_weights[0]=1;lm.m_num_constraints=1;lm.m_dirs[0]=btVector3(1,0,0);
    solver.m_objective->m_projection.m_lagrangeMultipliers.push_back(lm);solver.setupDeformableSolve(true);
    btSoftBody::DeformableNodeNodeContact c={};c.m_normal=btVector3(1,0,0);
    for(int n=0;n<2;++n){btSoftBody::ContactNode entry={&body.m_nodes[n],btMatrix3x3::getIdentity()*btScalar(n? .75:.25)};c.m_surfaceNodes.push_back(entry);}
    btDeformableContactForce contact(dt);ASSERT_TRUE(contact.add(c));
    solver.m_objective->m_KKTPreconditioner->reinitialize(true);
    contact.configureScaling(*solver.m_objective->m_KKTPreconditioner,solver.m_objective->m_projection.m_lagrangeMultipliers,2);
    solver.m_objective->m_lf.push_back(&contact);btScalar error=1;
    for(int i=0;i<20;++i){solver.solveDeformableConstraints(dt);error=contact.updateMultipliers();if(solver.m_lastSolveConverged && error<btScalar(.001))break;}
    EXPECT_TRUE(solver.m_lastSolveConverged);EXPECT_LT(error,btScalar(.001));
    EXPECT_NEAR(0.,double(body.m_nodes[0].m_v.x()),1e-6);
    EXPECT_NEAR(0.,double(btScalar(.25)*body.m_nodes[0].m_v.x()+btScalar(.75)*body.m_nodes[1].m_v.x()),.001);
    EXPECT_LE(body.m_nodes[1].m_v.length2(),btScalar(1));
    solver.m_objective->m_lf.clear();
}


TEST(CoupledContact, StaticSurfaceTrialRefreshUsesOneSidedDisplacement)
{
    btSoftBodyRigidBodyCollisionConfiguration config;btCollisionDispatcher dispatcher(&config);btDbvtBroadphase broadphase;
    btDeformableBodySolver solver;btDeformableMultiBodyConstraintSolver constraints;constraints.setDeformableSolver(&solver);
    ContactTestWorld world(&dispatcher,&broadphase,&constraints,&config,&solver);world.weighted=true;world.staticSurface=true;
    world.setImplicit(true);world.setCoupledContact(true);world.setMaxNewtonIterations(8);world.setNewtonTolerance(btScalar(1e-6));world.setGravity(btVector3(0,0,0));
    const btVector3 pa[]={btVector3(.005,0,0),btVector3(.015,1,0)},pb(0,0,0);const btScalar masses[]={1,1};
    btSoftBody a(&world.getWorldInfo(),2,pa,masses),reference(&world.getWorldInfo(),1,&pb,masses);
    a.m_cfg.drag=reference.m_cfg.drag=0;a.m_cfg.collisions=reference.m_cfg.collisions=0;
    for(int n=0;n<2;++n)a.m_nodes[n].m_v=a.m_nodes[n].m_vn=btVector3(-1,0,0);
    reference.m_nodes[0].m_v=reference.m_nodes[0].m_vn=btVector3(0,0,0);
    world.addSoftBody(&a);world.addSoftBody(&reference);
    ASSERT_EQ(1,world.stepSimulation(btScalar(.02),0));EXPECT_FALSE(world.hasCoupledStepFailed());
    EXPECT_GE(double(btScalar(.25)*a.m_nodes[0].m_x.x()+btScalar(.75)*a.m_nodes[1].m_x.x()),-1e-5);
    EXPECT_EQ(pb,reference.m_nodes[0].m_x);EXPECT_LT(a.m_nodes[0].m_v.length2()+a.m_nodes[1].m_v.length2(),btScalar(2));
    EXPECT_GT(world.detections,3);EXPECT_EQ(0,a.m_nodeRigidContacts.size());
    world.removeSoftBody(&a);world.removeSoftBody(&reference);
}


TEST(GImpactVertexCache, NestedQueriesRefreshAfterRollbackScalingAndSafeUpdate)
{
	btGImpactVertexCache cache;
	btVector3 current(1, 2, 3), safe(4, 5, 6), scale(2, 3, 4);
	int evaluations = 0;
	auto reconstruct = [&](int, btVector3& c, btVector3& s)
	{
		++evaluations;
		c = current * scale;
		s = safe * scale;
	};
	EXPECT_EQ(static_cast<const btVector3*>(nullptr), cache.current(0));
	for (int trial = 0; trial < 3; ++trial)
	{
		cache.begin(1, reconstruct);
		cache.begin(1, reconstruct);
		EXPECT_EQ(trial + 1, evaluations);
		EXPECT_EQ(current * scale, *cache.current(0));
		EXPECT_EQ(safe * scale, *cache.safe(0));
		btGImpactVertexCache clone(cache);
		EXPECT_EQ(static_cast<const btVector3*>(nullptr), clone.current(0));
		cache.end();
		EXPECT_NE(static_cast<const btVector3*>(nullptr), cache.current(0));
		cache.end();
		EXPECT_EQ(static_cast<const btVector3*>(nullptr), cache.current(0));
		EXPECT_EQ(static_cast<const btVector3*>(nullptr), cache.safe(0));
		current = trial == 0 ? btVector3(10, 20, 30) : btVector3(1, 2, 3);
		safe += btVector3(1, 0, 0);
		scale = btVector3(3, 2, 1);
	}
}

namespace
{
class CachedMappedContactManager : public MappedContactManager
{
public:
	mutable btGImpactVertexCache cache;
	mutable int evaluations = 0;
	explicit CachedMappedContactManager(MappedContactBody* body) : MappedContactManager(body) {}
	void begin_geometry_query() const override
	{
		cache.begin(static_cast<int>(body->mapping.size()), [this](int i, btVector3& c, btVector3& s)
		{
			++evaluations;
			MappedContactManager::get_vertex(i, c, false);
			s = c;
		});
	}
	void end_geometry_query() const override { cache.end(); }
	void get_vertex(unsigned int i, btVector3& v, bool original) const override
	{
		if (const auto* cached = cache.current(i)) v = *cached;
		else MappedContactManager::get_vertex(i, v, original);
	}
};
}

TEST(GImpactVertexCache, BvhRefitMatchesDirectGeometryAfterPositionAndMappingChanges)
{
	btSoftBodyWorldInfo info;
	const btVector3 positions[] = {btVector3(0,0,0), btVector3(2,0,0), btVector3(0,2,0), btVector3(0,0,2)};
	const btScalar masses[] = {1,1,1,1};
	MappedContactBody body(&info, positions, masses);
	body.appendTetra(0,1,2,3);
	body.mapping.resize(3);
	for (int i = 0; i < 3; ++i)
	{
		body.mapping[i].vertexToTetra = 0;
		body.mapping[i].baryCoordInTetra = btVector4(0,0,0,0);
		body.mapping[i].baryCoordInTetra[i] = 1;
	}
	int indices[] = {0,1,2, 2,1,0};
	btScalar vertices[] = {0,0,0, 2,0,0, 0,2,0};
	btTriangleIndexVertexArray mesh(2, indices, 3*sizeof(int), 3, vertices, 3*sizeof(btScalar));
	CachedMappedContactManager cached(&body);
	MappedContactManager direct(&body);
	cached.m_meshInterface = direct.m_meshInterface = &mesh;
	cached.m_part = direct.m_part = 0;
	cached.lock(); direct.lock();
	btGImpactQuantizedBvh tree(&cached), reference(&direct);
	tree.setStoreIndicesPerLevel();
	tree.buildSet();
	reference.buildSet();
	EXPECT_EQ(3, cached.evaluations);
	EXPECT_EQ(static_cast<const btVector3*>(nullptr), cached.cache.current(0));
	for (int trial = 0; trial < 3; ++trial)
	{
		body.m_nodes[0].m_x = trial == 1 ? btVector3(0,0,0) : btVector3(btScalar(.1),btScalar(.2),btScalar(.3));
		body.mapping[1].baryCoordInTetra = btVector4(btScalar(.25),btScalar(.75),0,0);
		cached.m_scale = direct.m_scale = btVector3(2,1,3);
		btAABB expected;
		direct.get_primitive_box(0, expected);
		tree.setGlobalBoundHint(expected);
		reference.setGlobalBoundHint(expected);
		tree.update();
		reference.update();
		EXPECT_EQ(3 * (trial + 2), cached.evaluations);
		EXPECT_EQ(static_cast<const btVector3*>(nullptr), cached.cache.current(0));
		const btAABB actual = tree.getGlobalBox();
		EXPECT_EQ(reference.getGlobalBox().m_min, actual.m_min);
		EXPECT_EQ(reference.getGlobalBox().m_max, actual.m_max);
		// Quantized bounds enclose the exact triangle, to within quantization resolution.
		for (int axis = 0; axis < 3; ++axis)
		{
			EXPECT_LE(actual.m_min[axis], expected.m_min[axis] + btScalar(1e-4));
			EXPECT_GE(actual.m_max[axis], expected.m_max[axis] - btScalar(1e-4));
		}
	}
	direct.unlock(); cached.unlock();
}


namespace
{
class GeometryRefreshTestWorld : public btDeformableMultiBodyDynamicsWorld
{
public:
	CachedMappedContactManager* manager = nullptr;
	int detections = 0;
	using btDeformableMultiBodyDynamicsWorld::btDeformableMultiBodyDynamicsWorld;
	void performDiscreteCollisionDetection() override
	{
		++detections;
		ASSERT_TRUE(manager != nullptr);
		EXPECT_TRUE(manager->cache.current(0) != nullptr);
		// Bounds and BVH work must not have reconstructed the same vertices again.
		EXPECT_EQ(3 * detections, manager->evaluations);
		btPrimitiveGeometryQuery nestedQuery(manager);
		for (int i = 0; i < 3; ++i)
		{
			btVector3 expected, actual;
			manager->MappedContactManager::get_vertex(i, expected, false);
			manager->get_vertex(i, actual, false);
			EXPECT_EQ(expected, actual);
		}
		EXPECT_EQ(3 * detections, manager->evaluations);
	}
};
}

TEST(GImpactVertexCache, EntireRefreshSharesCacheAndRebuildsAfterIntegration)
{
	btSoftBodyRigidBodyCollisionConfiguration config;
	btCollisionDispatcher dispatcher(&config);
	btDbvtBroadphase broadphase;
	btDeformableBodySolver solver;
	btDeformableMultiBodyConstraintSolver constraints;
	constraints.setDeformableSolver(&solver);
	GeometryRefreshTestWorld world(&dispatcher, &broadphase, &constraints, &config, &solver);
	world.setImplicit(true);
	world.setCoupledContact(true);
	world.setMaxNewtonIterations(8);
	world.setGravity(btVector3(0,0,0));
	const btVector3 positions[] = {btVector3(0,0,0),btVector3(2,0,0),btVector3(0,2,0),btVector3(0,0,2)};
	const btScalar masses[] = {1,1,1,1};
	MappedContactBody body(&world.getWorldInfo(), positions, masses);
	body.appendTetra(0,1,2,3);
	body.initializeDmInverse();
	body.m_tetraScratches.resize(1);
	body.m_tetraScratchesTn.resize(1);
	body.updateDeformation();
	body.mapping.resize(3);
	for (int i = 0; i < 3; ++i)
	{
		body.mapping[i].vertexToTetra = 0;
		body.mapping[i].baryCoordInTetra = btVector4(0,0,0,0);
		body.mapping[i].baryCoordInTetra[i] = 1;
	}
	int indices[] = {0,1,2};
	btScalar vertices[] = {0,0,0,2,0,0,0,2,0};
	btTriangleIndexVertexArray mesh(1, indices, 3*sizeof(int), 3, vertices, 3*sizeof(btScalar));
	auto* manager = new CachedMappedContactManager(&body);
	auto* shape = new btGImpactMeshShape(&mesh, manager);
	delete body.getCollisionShape();
	body.setCollisionShape(shape);
	shape->updateBound();
	body.m_cfg.drag = 0;
	body.m_cfg.collisions = 0;
	for (int n = 0; n < body.m_nodes.size(); ++n)
		body.m_nodes[n].m_v = body.m_nodes[n].m_vn = btVector3(1,0,0);
	world.addSoftBody(&body);
	world.manager = manager;
	manager->evaluations = 0;
	for (int step = 0; step < 2; ++step)
	{
		EXPECT_EQ(1, world.stepSimulation(btScalar(.001), 0));
		EXPECT_FALSE(world.hasCoupledStepFailed());
		EXPECT_TRUE(manager->cache.current(0) == nullptr);
		EXPECT_EQ(3 * world.detections, manager->evaluations);
	}
	EXPECT_GE(world.detections, 6);
	EXPECT_GT(body.m_nodes[0].m_x.x(), btScalar(0));
	world.removeSoftBody(&body);
}


TEST_F(DeformableBlockPreconditioner, AssembledSharedNodesMixedFallbackAndTimeStepMismatch)
{
	body->appendTetra(0,1,2,3);
	body->initializeDmInverse();
	body->m_tetraScratches.resize(2);
	for (int i = 0; i < 4; ++i) body->m_nodes[i].index = 3 - i;
	Vectors x, reference, actual;
	x.resize(4); reference.resize(4); actual.resize(4);
	for (int i = 0; i < 4; ++i) x[i] = btVector3(i+1, 2-i*i, i*3);
	force.m_useAssembledImplicit = true;
	for (int state = 0; state < 4; ++state)
	{
		for (int t = 0; t < 2; ++t)
		{
			body->m_tetraScratches[t].m_J = state == 1 && t == 1 ? btScalar(.001) : btScalar(1);
			body->m_tetraScratches[t].m_corotation = btMatrix3x3(btQuaternion(btVector3(1,2,3).normalized(), btScalar(.3 + t*.1)));
		}
		if (state == 3)
			for (int i = 0; i < 4; ++i) body->m_nodes[i].index = i;
		const btScalar dt = btScalar(.003);
		for (int i = 0; i < 4; ++i) reference[i] = actual[i] = btVector3(1,2,3);
		force.addImplicitForceDifferential(dt, x, reference);
		force.prepareImplicitForceDifferential(state == 2 ? dt*2 : dt);
		ASSERT_TRUE(force.m_assembledImplicitReady);
		EXPECT_EQ(size_t(12), force.m_implicitBlocks.size());
		force.addImplicitForceDifferential(dt, x, actual);
		force.finishImplicitForceDifferential();
		EXPECT_FALSE(force.m_assembledImplicitReady);
		for (int i = 0; i < 4; ++i) for (int axis = 0; axis < 3; ++axis)
			EXPECT_NEAR(double(reference[i][axis]), double(actual[i][axis]),
				1024 * SIMD_EPSILON * btMax(btScalar(1), btFabs(reference[i][axis])));
	}
	force.setDamping(0,0);
	force.prepareImplicitForceDifferential(btScalar(.003));
	for (int i = 0; i < 4; ++i) { x[i] = btVector3(1e8,-2e8,3e8); actual[i].setZero(); }
	force.addImplicitForceDifferential(btScalar(.003), x, actual);
	for (int i = 0; i < 4; ++i) EXPECT_EQ(btVector3(0,0,0), actual[i]);
	force.finishImplicitForceDifferential();
}

TEST(ElasticOperatorBenchmark, DISABLED_SharedGrid)
{
	const int side = 12, stride = side + 1;
	const int count = stride * stride * stride;
	std::vector<btVector3> positions;
	std::vector<btScalar> masses(count, 1);
	for (int z = 0; z <= side; ++z) for (int y = 0; y <= side; ++y) for (int x = 0; x <= side; ++x)
		positions.push_back(btVector3(x,y,z));
	btSoftBodyWorldInfo info;
	btSoftBody body(&info, count, positions.data(), masses.data());
	const int permutations[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
	const int offsets[3] = {1,stride,stride*stride};
	for (int z = 0; z < side; ++z) for (int y = 0; y < side; ++y) for (int x = 0; x < side; ++x)
	{
		const int a = x + stride*y + stride*stride*z;
		for (const auto& permutation : permutations)
		{
			int b = a + offsets[permutation[0]], c = b + offsets[permutation[1]];
			body.appendTetra(a,b,c,a+1+stride+stride*stride);
		}
	}
	body.initializeDmInverse();
	body.m_tetraScratches.resize(body.m_tetras.size());
	for (int t = 0; t < body.m_tetras.size(); ++t)
	{
		body.m_tetras[t].m_element_measure = btFabs(body.m_tetras[t].m_element_measure);
		body.m_tetraScratches[t].m_J = 1;
		body.m_tetraScratches[t].m_corotation = btMatrix3x3(btQuaternion(btVector3(1,2,3).normalized(),btScalar(.1*(t%7))));
	}
	btDeformableLinearElasticityForce force(1e6,2e6,btScalar(.01),btScalar(.01));
	force.addSoftBody(&body);
	Vectors input, output, reference;
	input.resize(count); output.resize(count); reference.resize(count);
	for (int n = 0; n < count; ++n)
	{
		body.m_nodes[n].index = n;
		input[n] = btVector3(btScalar(n%13)/13,btScalar(n%17)/17,btScalar(n%23)/23);
	}
	const btScalar dt = btScalar(.0002);
	for (int mode = 0; mode < 3; ++mode)
	{
		force.m_useAssembledImplicit = mode != 0;
		const auto start = std::chrono::steady_clock::now();
		force.prepareImplicitForceDifferential(dt);
		const auto prepared = std::chrono::steady_clock::now();
		for (int iteration = 0; iteration < 200; ++iteration)
		{
			for (int n = 0; n < count; ++n) output[n].setZero();
			force.addImplicitForceDifferential(dt,input,output);
		}
		const auto end = std::chrono::steady_clock::now();
		if (!mode) reference = output;
		else for (int n = 0; n < count; ++n) for (int axis = 0; axis < 3; ++axis)
			EXPECT_NEAR(double(reference[n][axis]),double(output[n][axis]),
				2048*SIMD_EPSILON*btMax(btScalar(1),btFabs(reference[n][axis])));
		printf("elastic mode=%s nodes=%d tets=%d blocks=%zu setup_ms=%.3f apply_200_ms=%.3f total_ms=%.3f\n",
			mode == 2 ? "assembled_reuse" : (mode ? "assembled_first" : "tetra"), count,body.m_tetras.size(),force.m_implicitBlocks.size(),
			std::chrono::duration<double,std::milli>(prepared-start).count(),
			std::chrono::duration<double,std::milli>(end-prepared).count(),
			std::chrono::duration<double,std::milli>(end-start).count());
		force.finishImplicitForceDifferential();
	}
}


TEST(ElasticOperatorInvestigation, DISABLED_SavedTriangle)
{
	const char* path = std::getenv("BULLET_ELASTIC_SCENE_PROBE");
	ASSERT_TRUE(path != nullptr);
	std::ifstream input(path);
	int nodes = 0, tets = 0;
	input >> nodes >> tets;
	ASSERT_GT(nodes, 0); ASSERT_GT(tets, 0);
	std::vector<btVector3> rest(nodes), current(nodes);
	std::vector<btScalar> mass(nodes, btScalar(.0836468325) / nodes);
	for (int n = 0; n < nodes; ++n)
		for (int c = 0; c < 6; ++c) input >> (c < 3 ? rest[n][c] : current[n][c-3]);
	btSoftBodyWorldInfo info;
	btSoftBody body(&info, nodes, rest.data(), mass.data());
	for (int t = 0; t < tets; ++t)
	{
		int a,b,c,d; input >> a >> b >> c >> d; body.appendTetra(a,b,c,d);
	}
	ASSERT_TRUE(bool(input));
	body.initializeDmInverse(); body.m_tetraScratches.resize(tets); body.m_tetraScratchesTn.resize(tets);
	btVector3 center(0,0,0);
	for (int n = 0; n < nodes; ++n) { body.m_nodes[n].m_x = current[n]; body.m_nodes[n].index = n; center += current[n]; }
	center /= nodes; body.updateDeformation();
	btDeformableLinearElasticityForce force;
	force.setYoungsModulus(btScalar(1e6)); force.setPoissonRatio(btScalar(.4)); force.setDamping(0,btScalar(.01)); force.addSoftBody(&body);
	const btScalar dt = btScalar(.0002);
	Vectors x, off, on, zero;
	zero.resize(nodes,btVector3(0,0,0)); x=off=on=zero;
	for (int direction = 0; direction < 4; ++direction)
	{
		for (int n = 0; n < nodes; ++n)
			x[n] = direction == 0 ? btVector3(btScalar(n%13)/13,btScalar(n%17)/17,btScalar(n%23)/23) :
				(direction == 1 ? btVector3(1,2,3) : btVector3(.3,.7,.2).cross((direction == 2 ? rest[n] : current[n])-center));
		off=on=zero;
		force.m_useAssembledImplicit=false; force.prepareImplicitForceDifferential(dt); force.addImplicitForceDifferential(dt,x,off); force.finishImplicitForceDifferential();
		force.m_useAssembledImplicit=true; force.prepareImplicitForceDifferential(dt); force.addImplicitForceDifferential(dt,x,on); force.finishImplicitForceDifferential();
		double difference=0,norm=0; btVector3 netOff(0,0,0),netOn(0,0,0);
		for(int n=0;n<nodes;++n){difference+=(off[n]-on[n]).length2();norm+=off[n].length2();netOff+=off[n];netOn+=on[n];}
		printf("SCENE_OPERATOR direction=%d off_norm=%.12g diff_norm=%.12g relative=%.12g net_off=%.12g net_on=%.12g\n",direction,sqrt(norm),sqrt(difference),sqrt(difference)/btMax(1.,sqrt(norm)),double(netOff.length()),double(netOn.length()));
		EXPECT_LT(sqrt(difference)/btMax(1.,sqrt(norm)),1e-8);
	}
	btAlignedObjectArray<btSoftBody*> bodies; bodies.push_back(&body);
	btDeformableBackwardEulerObjective objective(bodies, zero);
	objective.updateId(); objective.setDt(dt); objective.setImplicit(true); objective.m_lf.push_back(&force);
	Vectors rhs=zero, exact=zero, product=zero;
	for(int n=0;n<nodes;++n) exact[n]=btVector3(.3,.7,.2).cross(rest[n]-center);
	force.m_useAssembledImplicit=false; force.prepareImplicitForceDifferential(dt); objective.multiply(exact,rhs); force.finishImplicitForceDifferential();
	Vectors solved[2][2];
	for(int mode=0;mode<2;++mode) for(int tight=0;tight<2;++tight)
	{
		force.m_useAssembledImplicit=mode!=0; force.prepareImplicitForceDifferential(dt);
		objective.m_preconditioner->reinitialize(true); objective.setupTranslationCorrection();
		x=zero; objective.correctTranslation(x,rhs);
		btConjugateResidual<btDeformableBackwardEulerObjective> cr(2400);
		int iterations=cr.solveWithConvergencePolicy(objective,x,rhs,false,false,true,2400,tight?btScalar(1e-6):btScalar(.005));
		objective.correctTranslation(x,rhs); objective.multiply(x,product);
		double residual=0,error=0,solution=0;
		for(int n=0;n<nodes;++n){residual+=(rhs[n]-product[n]).length2();error+=(x[n]-exact[n]).length2();solution+=exact[n].length2();}
		printf("SCENE_SOLVE assembled=%d tight=%d iterations=%d residual=%.12g velocity_relative_error=%.12g\n",mode,tight,iterations,sqrt(residual),sqrt(error/solution));
		solved[mode][tight] = x;
		if (tight) { EXPECT_LT(sqrt(residual),1e-6); EXPECT_LT(sqrt(error/solution),1e-5); }
		force.finishImplicitForceDifferential(); objective.m_translationCorrection=false;
	}
	for (int tight = 0; tight < 2; ++tight)
	{
		double difference = 0, norm = 0;
		for (int n = 0; n < nodes; ++n) { difference += (solved[0][tight][n]-solved[1][tight][n]).length2(); norm += solved[0][tight][n].length2(); }
		printf("SCENE_SOLVE_COMPARE tight=%d velocity_relative_difference=%.12g\n",tight,sqrt(difference/norm));
		if (tight) EXPECT_LT(sqrt(difference/norm),1e-7);
	}
	objective.m_lf.clear();
}


TEST_F(DeformableBlockPreconditioner, RigidCorrectionBalancesForceTorqueAndIsSymmetric)
{
	Vectors zero, rhs, x, product, p, q, bp, bq;
	zero.resize(4,btVector3(0,0,0)); rhs=x=product=p=q=bp=bq=zero;
	btDeformableBackwardEulerObjective objective(bodies,zero);
	objective.m_rotationCorrectionEnabled=true;
	objective.updateId(); objective.setDt(btScalar(.0002)); objective.setImplicit(true);
	objective.m_lf.push_back(&force);
	force.prepareImplicitForceDifferential(btScalar(.0002));
	objective.m_preconditioner->reinitialize(true);
	ASSERT_TRUE(objective.setupTranslationCorrection());
	ASSERT_EQ(6,objective.m_translationBodies[0].modes);
	rhs[1]=btVector3(2,-3,5);
	objective.correctTranslation(x,rhs); objective.multiply(x,product);
	btVector3 linear(0,0,0), angular(0,0,0);
	for(int n=0;n<4;++n)
	{
		const btVector3 residual=rhs[n]-product[n];
		linear+=residual; angular+=body->m_nodes[n].m_x.cross(residual);
		p[n]=btVector3(n+1,2-n,n*n); q[n]=btVector3(3-n,n+2,-n);
	}
	EXPECT_LT(linear.length(),btScalar(1e-10));
	EXPECT_LT(angular.length(),btScalar(1e-10));
	objective.precondition(p,bp); objective.precondition(q,bq);
	const btScalar pbq=btDeformableTest::dot(p,bq),qbp=btDeformableTest::dot(q,bp);
	EXPECT_NEAR(double(pbq),double(qbp),1e-10*btMax(btScalar(1),btFabs(pbq)));
	EXPECT_GT(btDeformableTest::dot(p,bp),btScalar(0));
	for(int d=0;d<6;++d)
	{
		for(int n=0;n<4;++n) p[n]=objective.rigidMode(d,n);
		objective.multiply(p,product); objective.precondition(product,bp);
		for(int n=0;n<4;++n) EXPECT_LT((bp[n]-p[n]).length(),btScalar(1e-10));
	}
	force.finishImplicitForceDifferential();
	// A fresh geometry state rebuilds the rotation basis, without retained state flags.
	body->m_nodes[3].m_x+=btVector3(1,2,3);
	ASSERT_TRUE(objective.setupTranslationCorrection());
	EXPECT_EQ(6,objective.m_translationBodies[0].modes);
}

TEST_F(DeformableBlockPreconditioner, ContactCoarseCouplesBodiesAndPreservesSymmetry)
{
	const btVector3 position(3,2,1);const btScalar mass=2;
	btSoftBody other(&info,1,&position,&mass);bodies.push_back(&other);
	Vectors zero,p,q,bp,bq,az;zero.resize(5,btVector3(0,0,0));p=q=bp=bq=az=zero;
	btDeformableBackwardEulerObjective o(bodies,zero);o.updateId();o.setDt(.0002);o.setImplicit(true);
	o.m_rotationCorrectionEnabled=true;o.m_contactCoarseEnabled=true;
	btDeformableContactForce contact(.0002);btSoftBody::DeformableNodeNodeContact c={};
	c.m_node0=&body->m_nodes[1];c.m_node1=&other.m_nodes[0];c.m_normal=btVector3(1,0,0);c.m_offset=-.01;c.m_friction=.5;
	ASSERT_TRUE(contact.add(c));contact.contacts[0].normalImpulse=1;
	o.m_lf.push_back(&force);o.m_lf.push_back(&contact);force.prepareImplicitForceDifferential(.0002);
	o.m_preconditioner->reinitialize(true);ASSERT_TRUE(o.setupTranslationCorrection());ASSERT_TRUE(o.m_contactCoarse);
	ASSERT_EQ(9u,o.m_contactZ.size());
	EXPECT_GT(o.m_contactAZ[0][4].length2(),btScalar(0));
	for(int n=0;n<5;++n){p[n]=btVector3(n+1,2-n,n*n);q[n]=btVector3(3-n,n+2,-n);}
	o.precondition(p,bp);o.precondition(q,bq);
	EXPECT_NEAR(double(btDeformableTest::dot(p,bq)),double(btDeformableTest::dot(q,bp)),1e-9);
	EXPECT_GT(btDeformableTest::dot(p,bp),btScalar(0));
	for(int d=0;d<int(o.m_contactZ.size());++d)
	{
		o.multiply(o.m_contactZ[d],az);o.precondition(az,bp);
		for(int n=0;n<5;++n)EXPECT_LT((bp[n]-o.m_contactZ[d][n]).length(),btScalar(1e-9));
	}
	bp=zero;o.correctTranslation(bp,p);o.multiply(bp,az);
	for(int d=0;d<int(o.m_contactZ.size());++d)
	{
		btScalar balance=0;for(int n=0;n<5;++n)balance+=o.m_contactZ[d][n].dot(p[n]-az[n]);
		EXPECT_NEAR(double(balance),0.,1e-9);
	}
	o.m_contactCoarseEnabled=false;EXPECT_FALSE(o.setupTranslationCorrection());
	force.finishImplicitForceDifferential();bodies.pop_back();
}

TEST(RigidCorrection, CollinearGeometryRetainsTranslationAndSwitchDisablesRotation)
{
	btSoftBodyWorldInfo info;
	const btVector3 positions[]={btVector3(0,0,0),btVector3(1,0,0),btVector3(2,0,0)};
	const btScalar masses[]={1,2,3};
	btSoftBody body(&info,3,positions,masses);
	btAlignedObjectArray<btSoftBody*> bodies;bodies.push_back(&body);
	Vectors zero;zero.resize(3,btVector3(0,0,0));
	btDeformableBackwardEulerObjective objective(bodies,zero);
	objective.m_rotationCorrectionEnabled=true;
	objective.updateId();objective.setDt(btScalar(.01));objective.setImplicit(true);
	objective.m_preconditioner->reinitialize(true);
	ASSERT_TRUE(objective.setupTranslationCorrection());
	EXPECT_EQ(3,objective.m_translationBodies[0].modes);
	body.m_nodes[2].m_x=btVector3(0,1,0);
	ASSERT_TRUE(objective.setupTranslationCorrection());
	EXPECT_EQ(6,objective.m_translationBodies[0].modes);
	objective.m_rotationCorrectionEnabled=false;
	ASSERT_TRUE(objective.setupTranslationCorrection());
	EXPECT_EQ(3,objective.m_translationBodies[0].modes);
}
TEST_F(DeformableBlockPreconditioner, DampingGradientIncludesTrialGeometry)
{
	body->m_tetraScratchesTn.resize(1);
	body->advanceDeformation();
	Vectors v, d, impulse;
	v.resize(4);
	d.resize(4);
	impulse.resize(4, btVector3(0, 0, 0));
	const btScalar dt = .2, eps = 1e-5;
	for (int n = 0; n < 4; ++n)
	{
		v[n] = btVector3(n + 1, .3 * n, -.2);
		d[n] = btVector3(.2, n - 1, .4 * n);
	}
	auto trial = [&](btScalar a) {
		for (int n = 0; n < 4; ++n)
		{
			body->m_nodes[n].m_v = v[n] + a * d[n];
			body->m_nodes[n].m_q = body->m_nodes[n].m_x + dt * body->m_nodes[n].m_v;
		}
		body->updateDeformation();
	};
	trial(0);
	force.addScaledDampingForce(dt, impulse);
	double expected = 0;
	for (int n = 0; n < 4; ++n)
		expected -= impulse[n].dot(d[n]);
	trial(eps);
	double plus = force.totalDampingEnergy(dt);
	trial(-eps);
	double minus = force.totalDampingEnergy(dt);
	EXPECT_NEAR(expected, (plus - minus) / (2 * eps), 1e-6 * btMax(1., std::abs(expected)));
}

TEST_F(DeformableBlockPreconditioner, ElasticGradientWithRotatedShear)
{
	const btMatrix3x3 rotation(btQuaternion(btVector3(1, 2, 3).normalized(), .8));
	const btMatrix3x3 deformation = rotation * btMatrix3x3(1.2, .4, .1, 0, .8, .2, 0, 0, 1.1);
	Vectors q, d, f;
	q.resize(4);
	d.resize(4);
	f.resize(4, btVector3(0, 0, 0));
	for (int n = 0; n < 4; ++n)
	{
		q[n] = deformation * body->m_nodes[n].m_x;
		d[n] = btVector3(.2, n - 1, .4 * n);
	}
	auto trial = [&](btScalar a) {
		for (int n = 0; n < 4; ++n)
			body->m_nodes[n].m_q = q[n] + a * d[n];
		body->updateDeformation();
	};
	trial(0);
	force.addScaledElasticForce(1, f);
	double expected = 0;
	for (int n = 0; n < 4; ++n)
		expected -= f[n].dot(d[n]);
	const btScalar eps = 1e-5;
	trial(eps);
	double plus = force.totalElasticEnergy(1);
	trial(-eps);
	double minus = force.totalElasticEnergy(1);
	EXPECT_NEAR(expected, (plus - minus) / (2 * eps), 1e-6 * btMax(1., std::abs(expected)));
}

TEST_F(DeformableBlockPreconditioner, FrozenDampingOperatorsMatchWithDifferentRotations)
{
	body->m_tetraScratchesTn.resize(1);
	body->m_tetraScratchesTn[0] = body->m_tetraScratches[0];
	body->m_tetraScratchesTn[0].m_corotation = btMatrix3x3(btQuaternion(btVector3(3, 1, 2).normalized(), -.6));
	Vectors x, expected, actual;
	x.resize(4);
	expected.resize(4);
	actual.resize(4);
	for (int n = 0; n < 4; ++n)
		x[n] = btVector3(n - 2, n * n, 3 - n);
	for (int flat = 0; flat < 2; ++flat)
		for (int assemble = 0; assemble < 2; ++assemble)
		for (int currentFlat = 0; currentFlat < 2; ++currentFlat)
		{
			body->m_tetraScratchesTn[0].m_J = flat ? .001 : 1.;
			force.m_useAssembledImplicit = assemble != 0;
			compareBlocksToOperator(.02, currentFlat != 0);
			for (int n = 0; n < 4; ++n)
			{
				expected[n].setZero();
				actual[n].setZero();
			}
			force.addScaledDampingForceDifferential(-.02, x, expected);
			force.addScaledElasticForceDifferential(-.0004, x, expected);
			force.prepareImplicitForceDifferential(.02);
			force.addImplicitForceDifferential(.02, x, actual);
			force.finishImplicitForceDifferential();
			for (int n = 0; n < 4; ++n)
				EXPECT_LT((actual[n] - expected[n]).length(), 1e-9);
		}
}

TEST(NewtonReplay, DISABLED_CorrectedCapturedSolve)
{
	const char* path = std::getenv("BULLET_NEWTON_REPLAY");
	ASSERT_TRUE(path);
	btDeformableNewtonSnapshot r;
	ASSERT_TRUE(r.load(path));
	r.solver.updateState();
	r.solver.setMaxNewtonIterations(20);
	r.solver.solveDeformableConstraints(r.dt);
	printf("CORRECTED_CAPTURE converged=%d weighted=%.17g\n", int(r.solver.m_lastSolveConverged), double(r.weightedResidual()));
	EXPECT_TRUE(r.solver.m_lastSolveConverged);
}

TEST_F(DeformableBlockPreconditioner, PolarRotationIsProperForDegenerateTetrahedra)
{
	const btMatrix3x3 rotation(btQuaternion(btVector3(1, 2, 3).normalized(), 3.13));
	for (btScalar stretch : {btScalar(1), btScalar(1e-8), btScalar(0), btScalar(-.1)})
	{
		const btMatrix3x3 f = rotation * btMatrix3x3(1, .3, 0, 0, .8, .1, 0, 0, stretch);
		for (int n = 0; n < 4; ++n)
			body->m_nodes[n].m_q = f * body->m_nodes[n].m_x;
		body->updateDeformation();
		const auto &r = body->m_tetraScratches[0].m_corotation;
		EXPECT_NEAR(double(r.determinant()), 1., 1e-10);
		const btMatrix3x3 metric = r.transpose() * r;
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				EXPECT_NEAR(double(metric[i][j]), i == j ? 1. : 0., 1e-10);
	}
}

TEST_F(DeformableBlockPreconditioner, ContactBudgetFinishesSlowMultiplierConvergence)
{
	btDeformableContactForce contact(.01);
	btSoftBody::DeformableNodeNodeContact input = {};
	input.m_normal = btVector3(1, 0, 0);
	btSoftBody::ContactNode node = {&body->m_nodes[0], btMatrix3x3::getIdentity()};
	input.m_surfaceNodes.push_back(node);
	ASSERT_TRUE(contact.add(input));
	auto &c = contact.contacts[0];
	c.rho = .1;
	btDeformableContactConvergence budget;
	btScalar error = 1, errorAt20 = 0;
	int iterations = 0;
	while (iterations < budget.limit() && error > .001)
	{
		// Exact minimizer of the one-node inertial plus augmented-contact objective.
		body->m_nodes[0].m_v = btVector3((-1 + c.normalImpulse) / (1 + c.rho), 0, 0);
		error = contact.updateMultipliers();
		++iterations;
		if (iterations == 20)
			errorAt20 = error;
		budget.observe(iterations, error, true);
	}
	EXPECT_GT(errorAt20, .001);
	EXPECT_GT(iterations, 20);
	EXPECT_LT(iterations, 100);
	EXPECT_LE(error, .001);
	EXPECT_GE(c.normalImpulse, 0);
	EXPECT_NEAR(double(body->m_nodes[0].m_v.x()), 0., .001);
}
TEST(ContactConvergence, StopsStagnatingOrUnconvergedSolvesAndBoundsWork)
{
	for (int mode = 0; mode < 4; ++mode)
	{
		btDeformableContactConvergence budget;
		budget.observe(1, 1, true);
		EXPECT_FALSE(budget.observe(20,
			mode == 0	? 1
			: mode == 1 ? .1
			: mode == 2 ? SIMD_INFINITY
						: std::numeric_limits<btScalar>::quiet_NaN(),
			mode != 1));
		EXPECT_EQ(20, budget.limit());
	}
	btDeformableContactConvergence budget;
	budget.observe(1, 1, true);
	btScalar error = 1;
	for (int completed = 20; completed <= 100; completed += 20)
	{
		error *= .1;
		budget.observe(completed, error, true);
	}
	EXPECT_EQ(100, budget.limit());
}

TEST(ContactConvergence, PenaltyGrowthRequiresStagnationAndConvergedInnerSolve)
{
	btDeformableContactConvergence progress;
	EXPECT_FALSE(progress.increaseNormalPenalty(1, 1, true));
	EXPECT_FALSE(progress.increaseNormalPenalty(5, .1, true));
	EXPECT_FALSE(progress.increaseNormalPenalty(10, .1, false));
	EXPECT_TRUE(progress.increaseNormalPenalty(15, .1, true));
	EXPECT_EQ(4, progress.penaltyScale());
	EXPECT_TRUE(progress.increaseNormalPenalty(20, .1, true));
	EXPECT_TRUE(progress.increaseNormalPenalty(25, .1, true));
	EXPECT_FALSE(progress.increaseNormalPenalty(30, .1, true));
	EXPECT_EQ(64, progress.penaltyScale());
	EXPECT_NEAR(double(progress.normalError(.001)), .064, 1e-12);
	btDeformableContactConvergence accepted;
	accepted.increaseNormalPenalty(1, .0001, true);
	EXPECT_FALSE(accepted.increaseNormalPenalty(5, .0001, true));
}

TEST(NewtonReplay, DISABLED_AdaptiveCapturedContactSolve)
{
	const char *path = std::getenv("BULLET_NEWTON_REPLAY");
	ASSERT_TRUE(path);
	btDeformableNewtonSnapshot r;
	ASSERT_TRUE(r.load(path));
	btDeformableContactForce *c = nullptr;
	for (int f = 0; f < r.solver.m_objective->m_lf.size(); ++f)
		if (r.solver.m_objective->m_lf[f]->getForceType() == BT_CONTACT_FORCE)
			c = static_cast<btDeformableContactForce *>(r.solver.m_objective->m_lf[f]);
	ASSERT_TRUE(c);
	btDeformableContactConvergence progress;
	bool converged = false;
	for (int k = 0; k < progress.limit(); ++k)
	{
		r.solver.solveDeformableConstraints(r.dt);
		c->updateMultipliers();
		const auto normalError = progress.normalError(c->lastNormalError);
		const auto error = btMax(normalError, c->lastTangentError);
		if (error <= .001 && r.solver.m_lastSolveConverged)
		{
			converged = true;
			printf("ADAPTIVE_CAPTURE iterations=%d error=%.12g scale=%d\n", k + 1, double(error), progress.penaltyScale());
			break;
		}
		const bool increase = progress.increaseNormalPenalty(k + 1, normalError, r.solver.m_lastSolveConverged);
		if (increase)
			for (int j = 0; j < c->contacts.size(); ++j)
				c->contacts[j].rho *= 4;
		progress.observe(k + 1, error, r.solver.m_lastSolveConverged, increase);
	}
	EXPECT_TRUE(converged);
}

TEST(CoupledContact, FailedCallCanRetryWithoutResetOrAccumulatedForces)
{
	for (int fixed = 0; fixed < 2; ++fixed)
	{
		btSoftBodyRigidBodyCollisionConfiguration config;
		btCollisionDispatcher dispatcher(&config);
		btDbvtBroadphase broadphase;
		btDeformableBodySolver solver;
		btDeformableMultiBodyConstraintSolver constraints;
		constraints.setDeformableSolver(&solver);
		ContactTestWorld world(&dispatcher, &broadphase, &constraints, &config, &solver);
		world.setImplicit(true);
		world.setCoupledContact(true);
		world.setMaxNewtonIterations(8);
		world.setGravity(btVector3(0, 0, 0));
		world.setInternalTickCallback(countAccepted);
		world.setLatencyMotionStateInterpolation(false);
		const btVector3 pa(.005, 0, 0), pb(0, 0, 0);
		const btScalar mass = 1, dt = .125;
		btSoftBody a(&world.getWorldInfo(), 1, &pa, &mass), b(&world.getWorldInfo(), 1, &pb, &mass);
		a.m_cfg.drag = b.m_cfg.drag = 0;
		a.m_cfg.collisions = b.m_cfg.collisions = 0;
		a.m_nodes[0].m_v = a.m_nodes[0].m_vn = btVector3(-.01, 0, 0);
		b.m_nodes[0].m_v = b.m_nodes[0].m_vn = btVector3(0, 0, 0);
		world.addSoftBody(&a);
		world.addSoftBody(&b);
		acceptedCallbacks = 0;
		world.forceFailure = true;
		for (int attempt = 0; attempt < 2; ++attempt)
		{
			a.addForce(btVector3(0, 1, 0), 0);
			const int detections = world.detections;
			EXPECT_EQ(0, world.stepSimulation(dt, fixed, dt));
			EXPECT_TRUE(world.hasCoupledStepFailed());
			EXPECT_GT(world.detections, detections);
			EXPECT_EQ(pa, a.m_nodes[0].m_x);
			EXPECT_EQ(pb, b.m_nodes[0].m_x);
			EXPECT_EQ(btVector3(0, 0, 0), a.m_nodes[0].m_f);
			EXPECT_EQ(0, acceptedCallbacks);
			EXPECT_EQ(0, world.getLocalTime());
		}
		world.forceFailure = false;
		btAlignedObjectArray<int> indices;
		indices.push_back(0);
		btDeformableNodalForce persistent(&a, indices, btVector3(0, 1, 0));
		world.addForce(&persistent);
		EXPECT_EQ(1, world.stepSimulation(dt, fixed, dt));
		EXPECT_FALSE(world.hasCoupledStepFailed());
		EXPECT_EQ(1, acceptedCallbacks);
		EXPECT_NEAR(double(a.m_nodes[0].m_v.y()), double(dt), 1e-9);
		EXPECT_NEAR(double(a.m_nodes[0].m_x.y()), double(dt * dt), 1e-9);
		world.removeForce(&persistent);
		world.removeSoftBody(&a);
		world.removeSoftBody(&b);
	}
}

TEST(CoupledContact, LaterFailedSubstepRetainsCompletedTimeAndFractionalRemainder)
{
	class FailAfterFirst : public btDeformableBodySolver
	{
	public:
		bool fail = true;
		void solveDeformableConstraints(btScalar dt) override
		{
			if (fail && acceptedCallbacks == 1)
			{
				m_lastSolveConverged = false;
				return;
			}
			btDeformableBodySolver::solveDeformableConstraints(dt);
		}
	} solver;
	btSoftBodyRigidBodyCollisionConfiguration config;
	btCollisionDispatcher dispatcher(&config);
	btDbvtBroadphase broadphase;
	btDeformableMultiBodyConstraintSolver constraints;
	constraints.setDeformableSolver(&solver);
	ContactTestWorld world(&dispatcher, &broadphase, &constraints, &config, &solver);
	world.setImplicit(true);
	world.setCoupledContact(true);
	world.setMaxNewtonIterations(5);
	world.setGravity(btVector3(0, 0, 0));
	world.setInternalTickCallback(countAccepted);
	const btVector3 p(0, 0, 0);
	const btScalar mass = 1, dt = .125;
	btSoftBody body(&world.getWorldInfo(), 1, &p, &mass);
	body.m_cfg.drag = 0;
	body.m_cfg.collisions = 0;
	body.m_nodes[0].m_v = body.m_nodes[0].m_vn = btVector3(1, 0, 0);
	world.addSoftBody(&body);
	acceptedCallbacks = 0;
	EXPECT_EQ(1, world.stepSimulation(btScalar(2.5) * dt, 2, dt));
	EXPECT_TRUE(world.hasCoupledStepFailed());
	EXPECT_EQ(1, acceptedCallbacks);
	EXPECT_NEAR(double(body.m_nodes[0].m_x.x()), double(dt), 1e-9);
	EXPECT_EQ(btScalar(.5) * dt, world.getLocalTime());
	solver.fail = false;
	EXPECT_EQ(1, world.stepSimulation(btScalar(.5) * dt, 2, dt));
	EXPECT_FALSE(world.hasCoupledStepFailed());
	EXPECT_EQ(2, acceptedCallbacks);
	EXPECT_NEAR(double(body.m_nodes[0].m_x.x()), double(2 * dt), 1e-9);
	EXPECT_EQ(0, world.getLocalTime());
	world.removeSoftBody(&body);
}

TEST_F(DeformableBlockPreconditioner, WarmStartScalesImpulsesAndProjectsFrictionWithoutChangingPenalty)
{
	btCollisionObject a, b;
	btSoftBody::DeformableNodeNodeContact c = {};
	c.m_normal = btVector3(1, 0, 0);
	c.m_friction = .25;
	c.m_surfaceObjects[0] = &a;
	c.m_surfaceObjects[1] = &b;
	c.m_surfaceTriangles[0] = 1;
	c.m_surfaceTriangles[1] = 2;
	btSoftBody::ContactNode n = {&body->m_nodes[0], btMatrix3x3::getIdentity()};
	c.m_surfaceNodes.push_back(n);
	btDeformableContactForce previous(.01), current(.02);
	ASSERT_TRUE(previous.add(c));
	previous.contacts[0].normalImpulse = 2;
	previous.contacts[0].tangentImpulse = btVector3(0, 1, 0);
	c.m_normal = btVector3(1, .001, 0).normalized();
	ASSERT_TRUE(current.add(c));
	const auto rho = current.contacts[0].rho;
	EXPECT_EQ(1, current.warmStart(previous.contacts, previous.dt));
	EXPECT_EQ(4, current.contacts[0].normalImpulse);
	EXPECT_NEAR(0., double(current.contacts[0].normal.dot(current.contacts[0].tangentImpulse)), 1e-12);
	EXPECT_NEAR(1., double(current.contacts[0].tangentImpulse.length()), 1e-12);
	EXPECT_EQ(rho, current.contacts[0].rho);
}
TEST_F(DeformableBlockPreconditioner, WarmStartRejectsChangedIdentityStencilAndLargeTimestepRatio)
{
	btCollisionObject a, b;
	btSoftBody::DeformableNodeNodeContact c = {};
	c.m_normal = btVector3(1, 0, 0);
	c.m_friction = .5;
	c.m_surfaceObjects[0] = &a;
	c.m_surfaceObjects[1] = &b;
	c.m_surfaceTriangles[0] = 1;
	c.m_surfaceTriangles[1] = 2;
	btSoftBody::ContactNode n = {&body->m_nodes[0], btMatrix3x3::getIdentity()};
	c.m_surfaceNodes.push_back(n);
	btDeformableContactForce previous(.01);
	ASSERT_TRUE(previous.add(c));
	previous.contacts[0].normalImpulse = 2;
	for (int mode = 0; mode < 4; ++mode)
	{
		auto changed = c;
		if (mode == 0) changed.m_surfaceTriangles[0]++;
		if (mode == 1) changed.m_surfaceNodes[0].node = &body->m_nodes[1];
		if (mode == 2) changed.m_surfaceNodes[0].jacobian = changed.m_surfaceNodes[0].jacobian * .5;
		btDeformableContactForce current(mode == 3 ? .1 : .01);
		ASSERT_TRUE(current.add(changed));
		EXPECT_EQ(0, current.warmStart(previous.contacts, previous.dt));
		EXPECT_EQ(0, current.contacts[0].normalImpulse);
	}
	btDeformableContactForce current(.01);
	ASSERT_TRUE(current.add(c));
	c.m_surfaceNodes[0].jacobian = c.m_surfaceNodes[0].jacobian * 1.001;
	ASSERT_TRUE(current.add(c));
	ASSERT_EQ(2, current.contacts.size());
	EXPECT_EQ(1, current.warmStart(previous.contacts, previous.dt));
	EXPECT_EQ(2, current.contacts[0].normalImpulse);
	EXPECT_EQ(0, current.contacts[1].normalImpulse);
}

TEST(CoupledContact, WarmStartHistoryIsDiscardedAfterFailedOuterStep)
{
	class ProbeSolver : public btDeformableBodySolver
	{
	public:
		bool fail = false;
		btScalar firstSeed = -1;
		void solveDeformableConstraints(btScalar dt) override
		{
			if (firstSeed < 0)
			{
				firstSeed = 0;
				for (int f = 0; f < m_objective->m_lf.size(); ++f)
					if (m_objective->m_lf[f]->getForceType() == BT_CONTACT_FORCE)
					{
						auto* contact = static_cast<btDeformableContactForce*>(m_objective->m_lf[f]);
						for (int c = 0; c < contact->contacts.size(); ++c) firstSeed += contact->contacts[c].normalImpulse;
					}
			}
			if (fail)
			{
				m_lastSolveConverged = false;
				return;
			}
			btDeformableBodySolver::solveDeformableConstraints(dt);
		}
	} solver;
	class WarmWorld : public ContactTestWorld
	{
	public:
		using ContactTestWorld::ContactTestWorld;
		void performDiscreteCollisionDetection() override
		{
			auto& bodies = getSoftBodyArray();
			if (bodies.size() != 2) return;
			btSoftBody::DeformableNodeNodeContact c = {};
			c.m_node0 = &bodies[0]->m_nodes[0];
			c.m_node1 = &bodies[1]->m_nodes[0];
			c.m_normal = btVector3(1, 0, 0);
			c.m_offset = c.m_node0->m_x.x() - c.m_node1->m_x.x();
			c.m_surfaceObjects[0] = bodies[0];
			c.m_surfaceObjects[1] = bodies[1];
			c.m_surfaceTriangles[0] = c.m_surfaceTriangles[1] = 0;
			bodies[0]->m_nodeNodeContacts.push_back(c);
		}
	};
	btSoftBodyRigidBodyCollisionConfiguration config;
	btCollisionDispatcher dispatcher(&config);
	btDbvtBroadphase broadphase;
	btDeformableMultiBodyConstraintSolver constraints;
	constraints.setDeformableSolver(&solver);
	WarmWorld world(&dispatcher, &broadphase, &constraints, &config, &solver);
	world.setImplicit(true);
	world.setCoupledContact(true);
	world.setMaxNewtonIterations(8);
	world.setGravity(btVector3(0, 0, 0));
	const btVector3 p(0, 0, 0);
	const btScalar mass = 1, dt = .01;
	btSoftBody a(&world.getWorldInfo(), 1, &p, &mass), b(&world.getWorldInfo(), 1, &p, &mass);
	a.m_cfg.drag = b.m_cfg.drag = 0;
	a.m_cfg.collisions = b.m_cfg.collisions = 0;
	a.m_nodes[0].m_v = a.m_nodes[0].m_vn = btVector3(-1, 0, 0);
	b.m_nodes[0].m_v = b.m_nodes[0].m_vn = btVector3(0, 0, 0);
	world.addSoftBody(&a);
	world.addSoftBody(&b);
	ASSERT_EQ(1, world.stepSimulation(dt, 0));
	const auto position = a.m_nodes[0].m_x;
	solver.firstSeed = -1;
	solver.fail = true;
	EXPECT_EQ(0, world.stepSimulation(dt, 0));
	EXPECT_GT(solver.firstSeed, 0);
	EXPECT_EQ(position, a.m_nodes[0].m_x);
	solver.firstSeed = -1;
	solver.fail = false;
	EXPECT_EQ(1, world.stepSimulation(dt, 0));
	EXPECT_EQ(0, solver.firstSeed);
	world.removeSoftBody(&a);
	world.removeSoftBody(&b);
}

TEST(ElasticOperator, ParallelRowsMatchSerialAcrossRebuildAndFallback)
{
	const int side = 8, stride = side + 1;
	const int count = stride * stride * stride;
	std::vector<btVector3> positions;
	std::vector<btScalar> masses(count, 1);
	for (int z = 0; z <= side; ++z)
		for (int y = 0; y <= side; ++y)
			for (int x = 0; x <= side; ++x)
				positions.push_back(btVector3(x, y, z));
	btSoftBodyWorldInfo info;
	btSoftBody body(&info, count, positions.data(), masses.data());
	const int permutations[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
	const int offsets[3] = {1, stride, stride * stride};
	for (int z = 0; z < side; ++z)
		for (int y = 0; y < side; ++y)
			for (int x = 0; x < side; ++x)
			{
				const int a = x + stride * y + stride * stride * z;
				for (const auto& permutation : permutations)
				{
					int b = a + offsets[permutation[0]], c = b + offsets[permutation[1]];
					body.appendTetra(a, b, c, a + 1 + stride + stride * stride);
				}
			}
	body.initializeDmInverse();
	body.m_tetraScratches.resize(body.m_tetras.size());
	for (int t = 0; t < body.m_tetras.size(); ++t)
	{
		body.m_tetras[t].m_element_measure = btFabs(body.m_tetras[t].m_element_measure);
		body.m_tetraScratches[t].m_J = 1;
		body.m_tetraScratches[t].m_corotation = btMatrix3x3(btQuaternion(btVector3(1, 2, 3).normalized(), btScalar(.1 * (t % 7))));
	}
	btDeformableLinearElasticityForce force(1e6, 2e6, btScalar(.01), btScalar(.01));
	force.addSoftBody(&body);
	Vectors input, output, reference;
	input.resize(count);
	output.resize(count);
	reference.resize(count);
	for (int n = 0; n < count; ++n)
	{
		body.m_nodes[n].index = n;
		input[n] = btVector3(btScalar(n % 13) / 13, btScalar(n % 17) / 17, btScalar(n % 23) / 23);
	}
	const btScalar dt = btScalar(.0002);
	force.m_useAssembledImplicit = true;
	for (int state = 0; state < 4; ++state)
	{
		for (int n = 0; n < count; ++n) body.m_nodes[n].index = state == 1 ? count - 1 - n : n;
		for (int t = 0; t < body.m_tetras.size(); ++t) body.m_tetraScratches[t].m_J = state == 2 && t % 19 == 0 ? btScalar(.001) : btScalar(1);
		if (state == 3) body.setActivationState(ISLAND_SLEEPING);
		force.prepareImplicitForceDifferential(dt);
		force.setImplicitRowDispatcher({});
		for (int n = 0; n < count; ++n) reference[n] = btVector3(1, 2, 3);
		force.addImplicitForceDifferential(dt, input, reference);
		for (int workers : {1, 2, 4, 8})
		{
			int calls = 0;
			force.setImplicitRowDispatcher([&](int rows, const btDeformableLinearElasticityForce::ImplicitRowBody& work)
										   {
				++calls;
				std::vector<std::thread> threads;
				for (int w = 0; w < workers; ++w)
					threads.emplace_back([&,w]() { work(rows*w/workers, rows*(w+1)/workers); });
				for (auto& thread : threads) thread.join(); });
			for (int n = 0; n < count; ++n) output[n] = btVector3(1, 2, 3);
			force.addImplicitForceDifferential(dt, input, output);
			EXPECT_EQ(state == 3 ? 0 : 1, calls);
			for (int n = 0; n < count; ++n) EXPECT_EQ(reference[n], output[n]);
		}
		const btScalar aliasDt = btScalar(1e-8);
		force.prepareImplicitForceDifferential(aliasDt);
		force.setImplicitRowDispatcher({});
		reference = input;
		force.addImplicitForceDifferential(aliasDt, reference, reference);
		force.setImplicitRowDispatcher([](int, const btDeformableLinearElasticityForce::ImplicitRowBody&)
									   { ADD_FAILURE() << "Aliased product dispatched"; });
		output = input;
		force.addImplicitForceDifferential(aliasDt, output, output);
		for (int n = 0; n < count; ++n) EXPECT_EQ(reference[n], output[n]);

		force.finishImplicitForceDifferential();
	}
}

TEST_F(DeformableBlockPreconditioner, SmallAssembledProductDoesNotDispatch)
{
	force.m_useAssembledImplicit = true;
	force.setDamping(0, 0);
	force.setImplicitRowDispatcher([](int, const btDeformableLinearElasticityForce::ImplicitRowBody&)
								   { ADD_FAILURE() << "Small product dispatched"; });
	force.prepareImplicitForceDifferential(.01);
	Vectors input, output;
	input.resize(4, btVector3(1, 2, 3));
	output.resize(4, btVector3(0, 0, 0));
	force.addImplicitForceDifferential(.01, input, output);
	for (int i = 0; i < 4; ++i) EXPECT_EQ(btVector3(0, 0, 0), output[i]);
	force.finishImplicitForceDifferential();
}

TEST_F(DeformableBlockPreconditioner, ParallelDeformationMatchesSerialAndJoinsBeforeAdvance)
{
	for (int t = 1; t < 769; ++t) body->appendTetra(0, 1, 2, 3);
	body->initializeDmInverse();
	body->m_tetraScratches.resize(769);
	body->m_tetraScratchesTn.resize(769);
	for (int state = 0; state < 5; ++state)
	{
		btMatrix3x3 transform(btQuaternion(btVector3(1, 2, 3).normalized(), btScalar(.73)));
		if (state == 1) transform = btMatrix3x3(1, 0, 0, 0, 1, 0, 0, 0, 0);
		if (state == 2) transform = btMatrix3x3(-1, 0, 0, 0, 1, 0, 0, 0, 1);
		if (state == 3) transform = btMatrix3x3::getIdentity() * btScalar(0);
		if (state == 4) transform = btMatrix3x3(1, .2, 0, 0, .7, 0, 0, 0, 1e-12);
		for (int n = 0; n < 4; ++n) body->m_nodes[n].m_q = transform * body->m_nodes[n].m_x;
		body->setDeformationDispatcher({});
		body->updateDeformation();
		const auto expected = body->m_tetraScratches;
		for (int workers : {1, 2, 4, 8})
		{
			int calls = 0;
			body->setDeformationDispatcher([&](int count, const btSoftBody::DeformationRange& work)
										   {
				++calls;
				std::vector<std::thread> threads;
				for (int w = 0; w < workers; ++w)
					threads.emplace_back([&,w]() { work(count*w/workers, count*(w+1)/workers); });
				for (auto& thread : threads) thread.join(); });
			for (int t = 0; t < 769; ++t) body->m_tetraScratches[t].m_J = -123;
			body->advanceDeformation();
			EXPECT_EQ(1, calls);
			for (int t = 0; t < 769; ++t)
			{
				const auto& actual = body->m_tetraScratches[t];
				EXPECT_EQ(expected[t].m_J, actual.m_J);
				EXPECT_EQ(expected[t].m_trace, actual.m_trace);
				EXPECT_EQ(expected[t].m_J, body->m_tetraScratchesTn[t].m_J);
				for (int r = 0; r < 3; ++r)
				{
					EXPECT_EQ(expected[t].m_F[r], actual.m_F[r]);
					EXPECT_EQ(expected[t].m_F[r], body->m_tetras[t].m_F[r]);
					EXPECT_EQ(expected[t].m_cofF[r], actual.m_cofF[r]);
					EXPECT_EQ(expected[t].m_corotation[r], actual.m_corotation[r]);
					EXPECT_EQ(expected[t].m_corotation[r], body->m_tetraScratchesTn[t].m_corotation[r]);
				}
			}
		}
	}
	body->setDeformationDispatcher({});
}

TEST_F(DeformableBlockPreconditioner, SmallDeformationUpdateStaysSerial)
{
	body->setDeformationDispatcher([](int, const btSoftBody::DeformationRange&)
								   { ADD_FAILURE() << "Small body dispatched"; });
	body->updateDeformation();
	EXPECT_TRUE(std::isfinite(double(body->m_tetraScratches[0].m_J)));
	body->m_tetras.clear();
	body->updateDeformation();
}

TEST(ElasticOperator, ParallelAssemblyMatchesSerialAcrossRebuildAndFallback)
{
	const int side = 8, stride = side + 1;
	const int count = stride * stride * stride;
	std::vector<btVector3> positions;
	std::vector<btScalar> masses(count, 1);
	for (int z = 0; z <= side; ++z)
		for (int y = 0; y <= side; ++y)
			for (int x = 0; x <= side; ++x)
				positions.push_back(btVector3(x, y, z));
	btSoftBodyWorldInfo info;
	btSoftBody body(&info, count, positions.data(), masses.data());
	const int permutations[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
	const int offsets[3] = {1, stride, stride * stride};
	for (int z = 0; z < side; ++z)
		for (int y = 0; y < side; ++y)
			for (int x = 0; x < side; ++x)
			{
				const int a = x + stride * y + stride * stride * z;
				for (const auto& permutation : permutations)
				{
					int b = a + offsets[permutation[0]], c = b + offsets[permutation[1]];
					body.appendTetra(a, b, c, a + 1 + stride + stride * stride);
				}
			}
	body.initializeDmInverse();
	body.m_tetraScratches.resize(body.m_tetras.size());
	body.m_tetraScratchesTn.resize(body.m_tetras.size());
	for (int t = 0; t < body.m_tetras.size(); ++t)
	{
		body.m_tetras[t].m_element_measure = btFabs(body.m_tetras[t].m_element_measure);
		body.m_tetraScratches[t].m_J = 1;
		body.m_tetraScratches[t].m_corotation = btMatrix3x3(btQuaternion(btVector3(1, 2, 3).normalized(), btScalar(.1 * (t % 7))));
		body.m_tetraScratchesTn[t].m_J = t % 17 == 0 ? btScalar(.001) : btScalar(1);
		body.m_tetraScratchesTn[t].m_corotation = btMatrix3x3(btQuaternion(btVector3(3, 1, 2).normalized(), btScalar(.07 * (t % 11))));
	}
	btDeformableLinearElasticityForce force(1e6, 2e6, btScalar(.01), btScalar(.01));
	force.addSoftBody(&body);
	Vectors input, output, reference;
	input.resize(count);
	output.resize(count);
	reference.resize(count);
	for (int n = 0; n < count; ++n)
	{
		body.m_nodes[n].index = n;
		input[n] = btVector3(btScalar(n % 13) / 13, btScalar(n % 17) / 17, btScalar(n % 23) / 23);
	}
	const btScalar dt = btScalar(.0002);
	force.m_useAssembledImplicit = true;
	for (int state = 0; state < 4; ++state)
	{
		for (int n = 0; n < count; ++n) body.m_nodes[n].index = state == 1 ? count - 1 - n : n;
		for (int t = 0; t < body.m_tetras.size(); ++t) body.m_tetraScratches[t].m_J = state == 2 && t % 19 == 0 ? btScalar(.001) : btScalar(1);
		if (state == 3) body.setActivationState(ISLAND_SLEEPING);
		force.setImplicitAssemblyDispatcher({});
		force.prepareImplicitForceDifferential(dt);
		const auto serialBlocks = force.m_implicitBlocks;
		btAlignedObjectArray<btMatrix3x3> diagonalReference, diagonalParallel;
		diagonalReference.resize(count, btMatrix3x3::getIdentity());
		for (int n = 0; n < count; ++n)
		{
			body.m_nodes[n].m_im = n % 11 == 0 ? 0 : 1;
			body.m_nodes[n].m_frozen = n % 13 == 0 ? 1 : 0;
		}
		force.addImplicitForceDifferentialBlocks(dt, diagonalReference);
		for (int n = 0; n < count; ++n) reference[n] = btVector3(1, 2, 3);
		force.addImplicitForceDifferential(dt, input, reference);
		for (int workers : {1, 2, 4, 8, 31, 32})
		{
			int calls = 0;
			force.setImplicitAssemblyDispatcher([&](int rows, const btDeformableLinearElasticityForce::ImplicitRowBody& work)
												{
				++calls;
				std::vector<std::thread> threads;
				for (int w = 0; w < workers; ++w)
					threads.emplace_back([&,w]() { work(rows*w/workers, rows*(w+1)/workers); });
				for (auto& thread : threads) thread.join(); });
			force.prepareImplicitForceDifferential(dt);
			ASSERT_EQ(serialBlocks.size(), force.m_implicitBlocks.size());
			for (int i = 0; i < static_cast<int>(serialBlocks.size()); ++i)
				for (int r = 0; r < 3; ++r) EXPECT_EQ(serialBlocks[i].value[r], force.m_implicitBlocks[i].value[r]);
			for (int n = 0; n < count; ++n) output[n] = btVector3(1, 2, 3);
			force.addImplicitForceDifferential(dt, input, output);
			EXPECT_EQ(state == 3 ? 0 : 1, calls);
			diagonalParallel.resize(count);
			for (int n = 0; n < count; ++n) diagonalParallel[n] = btMatrix3x3::getIdentity();
			force.addImplicitForceDifferentialBlocks(dt, diagonalParallel);
			EXPECT_EQ(state == 3 ? 0 : 2, calls);
			for (int n = 0; n < count; ++n)
				for (int r = 0; r < 3; ++r) EXPECT_EQ(diagonalReference[n][r], diagonalParallel[n][r]);
			for (int n = 0; n < count; ++n) EXPECT_EQ(reference[n], output[n]);
		}
		force.finishImplicitForceDifferential();
	}
}

TEST_F(DeformableBlockPreconditioner, SmallAssemblyDoesNotDispatch)
{
	force.setImplicitAssemblyDispatcher([](int, const btDeformableLinearElasticityForce::ImplicitRowBody&)
										{ ADD_FAILURE() << "Small assembly dispatched"; });
	force.m_useAssembledImplicit = true;
	force.prepareImplicitForceDifferential(btScalar(.0002));
	EXPECT_TRUE(force.m_assembledImplicitReady);
	force.finishImplicitForceDifferential();
	force.m_useAssembledImplicit = false;
	force.prepareImplicitForceDifferential(btScalar(.0002));
	EXPECT_FALSE(force.m_assembledImplicitReady);
	force.finishImplicitForceDifferential();
}

TEST(ContactRefresh, PreservesCoplanarSupportIncludingStricterEquivalentStencil)
{
	btCollisionObject a, b;
	btSoftBody::Node nodes[2] = {};
	btSoftBody::DeformableNodeNodeContact old = {};
	old.m_surfaceObjects[0] = &a;
	old.m_surfaceObjects[1] = &b;
	old.m_surfaceParts[0] = old.m_surfaceParts[1] = 0;
	old.m_surfaceTriangles[0] = old.m_surfaceTriangles[1] = 0;
	old.m_normal = btVector3(0, 0, 1);
	old.m_offset = btScalar(-.1);
	old.m_surfaceNodes.push_back({&nodes[0], btMatrix3x3::getIdentity()});
	auto fresh = old;
	fresh.m_surfaceNodes[0].node = &nodes[1];
	btAlignedObjectArray<btSoftBody::DeformableNodeNodeContact> current, updated;
	current.push_back(old);
	updated.push_back(fresh);
	EXPECT_EQ(0, btRemoveRefreshedContactPatches(current, updated));
	ASSERT_EQ(1, current.size());
	EXPECT_EQ(&nodes[0], current[0].m_surfaceNodes[0].node);
	// Same nodes with different interpolation weights still provide distinct support.
	updated[0] = old;
	updated[0].m_surfaceNodes[0].jacobian = updated[0].m_surfaceNodes[0].jacobian * btScalar(.5);
	EXPECT_EQ(0, btRemoveRefreshedContactPatches(current, updated));
	updated[0] = old;
	updated[0].m_offset = btScalar(.2);
	EXPECT_EQ(0, btRemoveRefreshedContactPatches(current, updated));
	btDeformableContactForce force(btScalar(.002));
	EXPECT_TRUE(force.add(current[0]));
	EXPECT_TRUE(force.add(updated[0]));
	ASSERT_EQ(1, force.contacts.size());
	EXPECT_EQ(btScalar(-.1), force.contacts[0].gap);
	updated[0].m_surfaceObjects[0] = &b;
	updated[0].m_surfaceObjects[1] = &a;
	updated[0].m_normal *= -1;
	updated[0].m_surfaceNodes[0].jacobian = updated[0].m_surfaceNodes[0].jacobian * btScalar(-1);
	EXPECT_EQ(0, btRemoveRefreshedContactPatches(current, updated));
	updated[0].m_surfaceInvalid = true;
	EXPECT_EQ(0, btRemoveRefreshedContactPatches(current, updated));
}

TEST(ContactRefresh, CompatibleSupportSurvivesMixedNormalsInEitherOrder)
{
	btCollisionObject a, b;
	btSoftBody::DeformableNodeNodeContact old = {};
	old.m_surfaceObjects[0] = &a;
	old.m_surfaceObjects[1] = &b;
	old.m_surfaceParts[0] = old.m_surfaceParts[1] = 0;
	old.m_surfaceTriangles[0] = old.m_surfaceTriangles[1] = 0;
	old.m_normal = btVector3(0, 0, 1);
	auto rotated = old;
	rotated.m_normal = btVector3(1, 0, 0);
	for (int order = 0; order < 2; ++order)
	{
		btAlignedObjectArray<btSoftBody::DeformableNodeNodeContact> current, updated;
		current.push_back(old);
		updated.push_back(order ? old : rotated);
		updated.push_back(order ? rotated : old);
		EXPECT_EQ(0, btRemoveRefreshedContactPatches(current, updated));
		ASSERT_EQ(1, current.size());
	}
}

TEST(ContactRefresh, ReversedSelfContactKeepsCompatibleSupport)
{
	btCollisionObject body;
	btSoftBody::DeformableNodeNodeContact old = {};
	old.m_surfaceObjects[0] = old.m_surfaceObjects[1] = &body;
	old.m_surfaceParts[0] = old.m_surfaceParts[1] = 0;
	old.m_surfaceTriangles[0] = 3;
	old.m_surfaceTriangles[1] = 7;
	old.m_normal = btVector3(0, 0, 1);
	auto reversed = old;
	reversed.m_surfaceTriangles[0] = 7;
	reversed.m_surfaceTriangles[1] = 3;
	reversed.m_normal *= -1;
	btAlignedObjectArray<btSoftBody::DeformableNodeNodeContact> current, updated;
	current.push_back(old);
	updated.push_back(reversed);
	EXPECT_EQ(0, btRemoveRefreshedContactPatches(current, updated));
	ASSERT_EQ(1, current.size());
}

TEST(CoupledContact, FeasiblePredictorPreservesIncomingMomentum)
{
    btSoftBodyWorldInfo info;
    const btVector3 p[] = {btVector3(0,0,0),btVector3(1,0,0),btVector3(0,1,0),btVector3(0,0,1)};
    const btScalar masses[] = {btScalar(.001),btScalar(.001),btScalar(.001),btScalar(.001)};
    btSoftBody body(&info,4,p,masses);
    body.appendTetra(0,1,2,3);body.initializeDmInverse();
    body.m_tetraScratches.resize(1);body.m_tetraScratchesTn.resize(1);
    for(int n=0;n<4;++n)body.m_nodes[n].m_v=body.m_nodes[n].m_vn=btVector3(0,0,0);
    const btVector3 incoming(0,0,-12);
    body.m_nodes[3].m_v=body.m_nodes[3].m_vn=incoming;
    body.updateDeformation();
    btDeformableLinearElasticityForce elastic;
    elastic.setYoungsModulus(100);elastic.setPoissonRatio(btScalar(.4));elastic.setDamping(0,0);elastic.addSoftBody(&body);
    btDeformableVolumeBarrierForce barrier;
    btDeformableVolumeBarrierForce::Material material={&body,btScalar(100/.6)};barrier.materials.push_back(material);
    btAlignedObjectArray<btSoftBody*> bodies;bodies.push_back(&body);
    btDeformableBodySolver solver;solver.setImplicit(true);solver.setLineSearch(true);solver.m_useProjection=false;
    solver.setMaxNewtonIterations(50);solver.setNewtonTolerance(btScalar(1e-6));
    const btScalar dt=btScalar(.1);
    solver.reinitialize(bodies,dt);solver.setPreconditioner(btDeformableBackwardEulerObjective::KKT_preconditioner);
    solver.m_objective->m_lf.push_back(&elastic);solver.m_objective->m_lf.push_back(&barrier);
    solver.setupDeformableSolve(true);solver.updateState();
    ASSERT_FALSE(barrier.admissible());
    btDeformableDiagnostics::StepScope diagnostic(&solver,0,dt);
    solver.solveDeformableConstraints(dt);solver.updateState();
    EXPECT_TRUE(solver.m_lastSolveConverged);EXPECT_FALSE(solver.m_lastSolveInvalidPredictor);EXPECT_TRUE(barrier.admissible());
    Vectors balance;balance.resize(4,btVector3(0,0,0));elastic.addScaledForces(dt,balance);barrier.addScaledForces(dt,balance);
    balance[3]-=(body.m_nodes[3].m_v-incoming)/body.m_nodes[3].m_im;
    for(int n=0;n<3;++n)balance[n]-=body.m_nodes[n].m_v/body.m_nodes[n].m_im;
    for(int n=0;n<4;++n)EXPECT_LT(balance[n].length(),btScalar(1e-6));
    btVector3 momentum(0,0,0);for(int n=0;n<4;++n)momentum+=body.m_nodes[n].m_v/body.m_nodes[n].m_im;
    EXPECT_LT((momentum-incoming*btScalar(.001)).length(),btScalar(1e-6));
    EXPECT_GT((body.m_nodes[3].m_v-incoming).length(),btScalar(1));
    EXPECT_EQ(p[3],body.m_nodes[3].m_x);
    solver.m_objective->m_lf.clear();
}

TEST(CoupledContact, StaticPlaneSupportsUseLocalClearanceAndFiniteTriangle)
{
    btSoftBodyWorldInfo info;
    const btVector3 p[]={btVector3(0,0,.1),btVector3(2,0,.15),btVector3(0,2,.18),btVector3(0,0,2)};
    const btScalar masses[]={1,1,1,1};
    int indices[]={0,1,2};btScalar vertices[]={0,0,0,2,0,0,0,2,0};
    btTriangleIndexVertexArray mesh(1,indices,3*sizeof(int),3,vertices,3*sizeof(btScalar));
    MappedContactBody body(&info,p,masses);body.appendTetra(0,1,2,3);body.mapping.resize(3);
    for(int v=0;v<3;++v){body.mapping[v].vertexToTetra=0;body.mapping[v].baryCoordInTetra=btVector4(0,0,0,0);body.mapping[v].baryCoordInTetra[v]=1;}
    auto* shape=new btGImpactMeshShape(&mesh,new MappedContactManager(&body));delete body.getCollisionShape();body.setCollisionShape(shape);
    shape->setMargin(btScalar(.1));shape->updateBound();
    btGImpactMeshShape rigidShape(&mesh);rigidShape.setMargin(btScalar(.1));rigidShape.updateBound();
    btCollisionObject rigid;rigid.setCollisionShape(&rigidShape);rigid.setCollisionFlags(btCollisionObject::CF_STATIC_OBJECT);
    for(int rotated=0;rotated<2;++rotated)
    {
        btTransform transform;transform.setIdentity();
        if(rotated){transform.setRotation(btQuaternion(btVector3(1,2,3).normalized(),btScalar(.7)));transform.setOrigin(btVector3(4,-3,2));}
        body.setWorldTransform(transform);rigid.setWorldTransform(transform);
        btCollisionObjectWrapper wrapper(nullptr,&rigidShape,&rigid,transform,-1,-1);
        const btVector3 normal=transform.getBasis()*btVector3(0,0,1);
        body.m_nodeNodeContacts.clear();
        body.skinSoftStaticCollisionHandler(&wrapper,0,0,0,0,transform*p[0],normal,btScalar(.1),false,nullptr);
        ASSERT_EQ(4,body.m_nodeNodeContacts.size());
        for(int v=0;v<3;++v)
        {
            const auto& c=body.m_nodeNodeContacts[v+1];
            EXPECT_NEAR(double(p[v].z()-.195),double(c.m_offset),1e-7);
            btVector3 displacement(0,0,0);
            for(int n=0;n<c.m_surfaceNodes.size();++n)
                if(c.m_surfaceNodes[n].node==&body.m_nodes[v])displacement+=c.m_surfaceNodes[n].jacobian*btVector3(0,0,1);
            EXPECT_NEAR(1.,double(normal.dot(displacement)),1e-8);
        }
    }
    btTransform identity;identity.setIdentity();body.setWorldTransform(identity);
    btTransform translated=identity;translated.setOrigin(btVector3(1,0,0));rigid.setWorldTransform(translated);
    btCollisionObjectWrapper finite(nullptr,&rigidShape,&rigid,translated,-1,-1);
    body.m_nodeNodeContacts.clear();
    body.skinSoftStaticCollisionHandler(&finite,0,0,0,0,p[1],btVector3(0,0,1),btScalar(.15),false,nullptr);
    ASSERT_EQ(2,body.m_nodeNodeContacts.size());
    EXPECT_NEAR(-.045,double(body.m_nodeNodeContacts[1].m_offset),1e-7);
}

TEST(DeformableCRContinuation, InexactCorrectionVerifiesResidualAndDoesNotLeakToStrictSolve)
{
    struct Matrix
    {
        void multiply(const Vectors& x,Vectors& y){y[0]=x[0]*btVector3(1,10,100);}
        void precondition(const Vectors& x,Vectors& y){y=x;}
    } matrix;
    Vectors rhs,x,product;
    rhs.resize(1,btVector3(10000,1,1));x.resize(1,btVector3(0,0,0));product=x;
    btConjugateResidual<Matrix> cr(100);
    const int coarse=cr.solveWithConvergencePolicy(matrix,x,rhs,false,false,true,100,0,btScalar(1e-6),btScalar(.01));
    matrix.multiply(x,product);
    EXPECT_LE((rhs[0]-product[0]).length(),cr.getTargetResidual());
    EXPECT_NEAR(double((rhs[0]-product[0]).length()),double(cr.getFinalResidual()),1e-8);
    EXPECT_GT(cr.getFinalResidual(),btScalar(1));
    EXPECT_LT(coarse,3);
    cr.solveWithConvergencePolicy(matrix,x,rhs,false,false,true,100,0,btScalar(1e-6));
    matrix.multiply(x,product);
    EXPECT_LE(cr.getTargetResidual(),btScalar(1e-8));
    EXPECT_LT((rhs[0]-product[0]).length(),btScalar(1e-8));
}

TEST(CoupledContact, RoundoffStencilEntriesDoNotCreateDuplicateMultipliers)
{
    btSoftBodyWorldInfo info;const btVector3 p[]={btVector3(0,0,0),btVector3(1,0,0)};const btScalar masses[]={1,1};
    btSoftBody body(&info,2,p,masses);
    btSoftBody::DeformableNodeNodeContact c={};c.m_normal=btVector3(0,0,1);c.m_offset=btScalar(.1);
    btSoftBody::ContactNode primary={&body.m_nodes[0],btMatrix3x3::getIdentity()};c.m_surfaceNodes.push_back(primary);
    btDeformableContactForce force(btScalar(.01));ASSERT_TRUE(force.add(c));
    btSoftBody::ContactNode roundoff={&body.m_nodes[1],btMatrix3x3::getIdentity()*btScalar(1e-14)};c.m_surfaceNodes.push_back(roundoff);
    c.m_offset=btScalar(.05);ASSERT_TRUE(force.add(c));
    ASSERT_EQ(1,force.contacts.size());EXPECT_EQ(btScalar(.05),force.contacts[0].gap);
    c.m_surfaceNodes[1].jacobian=btMatrix3x3::getIdentity()*btScalar(1e-6);
    ASSERT_TRUE(force.add(c));EXPECT_EQ(2,force.contacts.size());
    c.m_surfaceNodes[0].jacobian=btMatrix3x3::getIdentity()*btScalar(0);
    c.m_surfaceNodes[1].jacobian=btMatrix3x3::getIdentity()*btScalar(0);
    EXPECT_FALSE(force.add(c));
    c.m_surfaceNodes.resize(1);c.m_surfaceNodes[0].jacobian=btMatrix3x3::getIdentity()*btScalar(1e-15);
    ASSERT_TRUE(force.add(c));EXPECT_EQ(1,force.contacts[force.contacts.size()-1].nodes.size());
}

TEST(CoupledContact, StalledWarmStartRetriesColdBeforeSubdividing)
{
	class ProbeSolver : public btDeformableBodySolver
	{
	public:
		bool fail = false;
		int rejectedWarmCalls = 0;
		bool recoveredCold = false;
		btScalar firstSeed = -1;
		void solveDeformableConstraints(btScalar dt) override
		{
			if (firstSeed < 0)
			{
				firstSeed = 0;
				for (int f = 0; f < m_objective->m_lf.size(); ++f)
					if (m_objective->m_lf[f]->getForceType() == BT_CONTACT_FORCE)
					{
						auto* contact = static_cast<btDeformableContactForce*>(m_objective->m_lf[f]);
						for (int c = 0; c < contact->contacts.size(); ++c) firstSeed += contact->contacts[c].normalImpulse;
					}
			}
			if (fail)
			{
				btScalar seed = 0;
				for (int f = 0; f < m_objective->m_lf.size(); ++f)
					if (m_objective->m_lf[f]->getForceType() == BT_CONTACT_FORCE)
						{
							const auto& contacts = static_cast<btDeformableContactForce*>(m_objective->m_lf[f])->contacts;
							for (int c = 0; c < contacts.size(); ++c) seed += contacts[c].normalImpulse;
						}
				if (seed > 0)
				{
					++rejectedWarmCalls;
					m_lastSolveConverged = false;
					return;
				}
				recoveredCold = true;
				fail = false;
			}
			btDeformableBodySolver::solveDeformableConstraints(dt);
		}
	} solver;
	class WarmWorld : public ContactTestWorld
	{
	public:
		using ContactTestWorld::ContactTestWorld;
		void performDiscreteCollisionDetection() override
		{
			auto& bodies = getSoftBodyArray();
			if (bodies.size() != 2) return;
			btSoftBody::DeformableNodeNodeContact c = {};
			c.m_node0 = &bodies[0]->m_nodes[0];
			c.m_node1 = &bodies[1]->m_nodes[0];
			c.m_normal = btVector3(1, 0, 0);
			c.m_offset = c.m_node0->m_x.x() - c.m_node1->m_x.x();
			c.m_surfaceObjects[0] = bodies[0];
			c.m_surfaceObjects[1] = bodies[1];
			c.m_surfaceTriangles[0] = c.m_surfaceTriangles[1] = 0;
			bodies[0]->m_nodeNodeContacts.push_back(c);
		}
	};
	btSoftBodyRigidBodyCollisionConfiguration config;
	btCollisionDispatcher dispatcher(&config);
	btDbvtBroadphase broadphase;
	btDeformableMultiBodyConstraintSolver constraints;
	constraints.setDeformableSolver(&solver);
	WarmWorld world(&dispatcher, &broadphase, &constraints, &config, &solver);
	world.setImplicit(true);
	world.setCoupledContact(true);
	world.setMaxNewtonIterations(8);
	world.setGravity(btVector3(0, 0, 0));
	const btVector3 p(0, 0, 0);
	const btScalar mass = 1, dt = .01;
	btSoftBody a(&world.getWorldInfo(), 1, &p, &mass), b(&world.getWorldInfo(), 1, &p, &mass);
	a.m_cfg.drag = b.m_cfg.drag = 0;
	a.m_cfg.collisions = b.m_cfg.collisions = 0;
	a.m_nodes[0].m_v = a.m_nodes[0].m_vn = btVector3(-1, 0, 0);
	b.m_nodes[0].m_v = b.m_nodes[0].m_vn = btVector3(0, 0, 0);
	world.addSoftBody(&a);
	world.addSoftBody(&b);
	ASSERT_EQ(1, world.stepSimulation(dt, 0));
	solver.firstSeed = -1;
	solver.fail = true;
	EXPECT_EQ(1, world.stepSimulation(dt, 0));
	EXPECT_GT(solver.firstSeed, 0);
	EXPECT_EQ(20, solver.rejectedWarmCalls);
	EXPECT_TRUE(solver.recoveredCold);
	EXPECT_FALSE(world.hasCoupledStepFailed());
	EXPECT_GE(double(a.m_nodes[0].m_x.x()-b.m_nodes[0].m_x.x()), -1e-5);
	world.removeSoftBody(&a);
	world.removeSoftBody(&b);
}

#include "deformable_vbd_tests.h"
