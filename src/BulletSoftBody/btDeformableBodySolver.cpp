/*
 Written by Xuchen Han <xuchenhan2015@u.northwestern.edu>
 
 Bullet Continuous Collision Detection and Physics Library
 Copyright (c) 2019 Google Inc. http://bulletphysics.org
 This software is provided 'as-is', without any express or implied warranty.
 In no event will the authors be held liable for any damages arising from the use of this software.
 Permission is granted to anyone to use this software for any purpose,
 including commercial applications, and to alter it and redistribute it freely,
 subject to the following restrictions:
 1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
 2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
 3. This notice may not be removed or altered from any source distribution.
 */

#include <limits>
#include <stdio.h>
#include "btDeformableBodySolver.h"
#include "btSoftBodyInternals.h"
#include "LinearMath/btQuickprof.h"
static const int kMaxConjugateGradientIterations = 300;

namespace
{
class ImplicitOperatorCacheScope
{
	btDeformableBackwardEulerObjective& m_objective;
	bool m_enabled;
public:
	explicit ImplicitOperatorCacheScope(btDeformableBackwardEulerObjective& objective)
		: m_objective(objective), m_enabled(BT_DEFORMABLE_USE_CACHED_OPERATOR && objective.m_implicit)
	{
		if (!m_enabled) return;
		btDeformablePerformanceScope perf(objective.m_performance.cacheSetup, objective.m_performance.active);
		for (int f = 0; f < objective.m_lf.size(); ++f)
			objective.m_lf[f]->prepareImplicitForceDifferential(objective.m_dt);
	}
	~ImplicitOperatorCacheScope()
	{
		if (!m_enabled) return;
		for (int f = 0; f < m_objective.m_lf.size(); ++f)
			m_objective.m_lf[f]->finishImplicitForceDifferential();
	}
};
}


btDeformableBodySolver::btDeformableBodySolver()
	: m_numNodes(0), m_cg(kMaxConjugateGradientIterations), m_cr(kMaxConjugateGradientIterations), m_maxNewtonIterations(1), m_newtonTolerance(1e-4), m_lineSearch(false), m_useProjection(false)
{
	m_objective = new btDeformableBackwardEulerObjective(m_softBodies, m_backupVelocity);
	m_reducedSolver = false;
}

btDeformableBodySolver::~btDeformableBodySolver()
{
	delete m_objective;
}

void btDeformableBodySolver::solveDeformableConstraints(btScalar solverdt)
{
	BT_PROFILE("solveDeformableConstraints");
	btDeformableSolverPerformance& perf = m_objective->m_performance;
	perf = btDeformableSolverPerformance();
	perf.active = m_implicit;
	++m_performanceStep;
	const std::chrono::steady_clock::time_point perfStart = std::chrono::steady_clock::now();
	if (!m_implicit)
	{
		m_objective->computeResidual(solverdt, m_residual);
		m_objective->applyDynamicFriction(m_residual);
		if (m_useProjection)
		{
			computeStep(m_dv, m_residual);
		}
		else
		{
			TVStack rhs, x;
			m_objective->addLagrangeMultiplierRHS(m_residual, m_dv, rhs);
			m_objective->addLagrangeMultiplier(m_dv, x);
			{
				btDeformablePerformanceScope blockPerf(m_objective->m_performance.blockSetup, m_objective->m_performance.active);
				m_objective->m_preconditioner->reinitialize(true);
			}
			// Explicit damping needs the previous absolute-tolerance/full-budget
			// policy: a small relative residual can still leave large velocity errors.
			{
				btDeformablePerformanceScope linearPerf(m_objective->m_performance.linear, m_objective->m_performance.active);
				m_objective->m_performance.krylovIterations += m_cr.solveWithConvergencePolicy(*m_objective, x, rhs, false, false);
			}
			for (int i = 0; i < m_dv.size(); ++i)
			{
				m_dv[i] = x[i];
			}
		}
		updateVelocity();
	}
	else
	{
		m_implicitRecoveryUsed = false;
		for (int i = 0; i < m_maxNewtonIterations; ++i)
		{
			m_newtonIteration = i;
			++perf.newton;
			btDeformableSolverPerformance::NewtonRecord newton;
			newton.iteration = i + 1;
			perf.newtonRecords.push_back(newton);
			btDeformableSolverPerformance::NewtonRecord& nr = perf.newtonRecords[perf.newtonRecords.size() - 1];
			updateState();
			// add the inertia term in the residual
			int counter = 0;
			for (int k = 0; k < m_softBodies.size(); ++k)
			{
				btSoftBody* psb = m_softBodies[k];
				for (int j = 0; j < psb->m_nodes.size(); ++j)
				{
					if (psb->m_nodes[j].m_frozen <= 0 && psb->m_nodes[j].m_im > 0)
					{
						m_residual[counter] = (-1. / psb->m_nodes[j].m_im) * m_dv[counter];
					}
					++counter;
				}
			}

			m_objective->computeResidual(solverdt, m_residual);
			const btScalar residualNorm = m_objective->computeNorm(m_residual);
			const btScalar constraintError = implicitConstraintError(m_dv);
			nr.inputForceL2 = residualNorm;
			nr.inputConstraintInf = constraintError;
			nr.outputConstraintInf = constraintError;
			if (residualNorm < m_newtonTolerance && constraintError <= m_newtonTolerance && i > 0)
			{
				nr.outcome = "input_residual_target";
				break;
			}
			// todo xuchenhan@: this really only needs to be calculated once
			m_objective->applyDynamicFriction(m_residual);
			if (m_lineSearch)
			{
				btScalar inner_product = computeDescentStep(m_ddv, m_residual);
				const btScalar stepNorm = m_objective->computeNorm(m_ddv);
				const btScalar relativeStepTolerance = m_newtonTolerance * btMax(btScalar(1), m_objective->computeNorm(m_dv));
				nr.solved = true; nr.stepL2 = stepNorm; nr.stepTolerance = relativeStepTolerance;
				nr.stationarityL2 = m_lastStationarityResidual;
				nr.linearMomentumL2 = m_lastLinearMomentumResidual;
				nr.linearConstraintInf = m_lastLinearConstraintResidual;
				const bool accurateLinearStep = !m_useProjection &&
					m_lastLinearMomentumResidual <= m_newtonTolerance &&
					m_lastLinearConstraintResidual <= m_newtonTolerance;
				const bool stepConverged = accurateLinearStep &&
					m_lastStationarityResidual <= m_newtonTolerance && constraintError <= m_newtonTolerance;
				if (i > 0 && stepNorm <= relativeStepTolerance && (stepConverged || !accurateLinearStep))
				{
					nr.outcome = stepConverged ? "small_stationary_step" : "small_inaccurate_step";
					// Check stationarity at the accepted state, not just the solved
					// linear system. Otherwise apply a small but accurate correction.
					break;
				}
				btScalar alpha = 0.01, beta = 0.5;  // Boyd & Vandenberghe suggested alpha between 0.01 and 0.3, beta between 0.1 to 0.8
				btScalar scale = 2;
				btScalar f0 = m_objective->totalEnergy(solverdt) + kineticEnergy(), f1, f2;
				bool lineSearchFailed = false;
				backupDv();
				do
				{
					scale *= beta;
					if (scale < 1e-8)
					{
						lineSearchFailed = true;
						++perf.lineFailures;
						break;
					}
					++perf.lineTrials;
					updateEnergy(scale);
					f1 = m_objective->totalEnergy(solverdt) + kineticEnergy();
					f2 = f0 - alpha * scale * inner_product;
				} while (!(f1 < f2 + SIMD_EPSILON));  // if anything here is nan then the search continues
				if (lineSearchFailed)
				{
					nr.outcome = "line_search_failed";
					// The trial evaluations modified m_dv, node velocities, temporary
					// positions, and deformation scratch data. Restore the accepted
					// iterate before terminating this Newton solve.
					revertDv();
					updateState();
					break;
				}
				revertDv();
				updateDv(scale);
				nr.appliedScale = scale;
			}
			else
			{
				computeStep(m_ddv, m_residual);
				const btScalar stepNorm = m_objective->computeNorm(m_ddv);
				const btScalar relativeStepTolerance = m_newtonTolerance * btMax(btScalar(1), m_objective->computeNorm(m_dv));
				nr.solved = true; nr.stepL2 = stepNorm; nr.stepTolerance = relativeStepTolerance;
				nr.stationarityL2 = m_lastStationarityResidual;
				nr.linearMomentumL2 = m_lastLinearMomentumResidual;
				nr.linearConstraintInf = m_lastLinearConstraintResidual;
				const bool accurateLinearStep = !m_useProjection &&
					m_lastLinearMomentumResidual <= m_newtonTolerance &&
					m_lastLinearConstraintResidual <= m_newtonTolerance;
				const bool stepConverged = accurateLinearStep &&
					m_lastStationarityResidual <= m_newtonTolerance && constraintError <= m_newtonTolerance;
				if (i > 0 && stepNorm <= relativeStepTolerance && (stepConverged || !accurateLinearStep))
				{
					nr.outcome = stepConverged ? "small_stationary_step" : "small_inaccurate_step";
					// Check stationarity at the accepted state, not just the solved
					// linear system. Otherwise apply a small but accurate correction.
					break;
				}
				updateDv();
				nr.appliedScale = 1;
			}
			nr.outputConstraintInf = implicitConstraintError(m_dv);
			nr.outcome = i + 1 == m_maxNewtonIterations ? "iteration_limit" : "applied";
			for (int j = 0; j < m_numNodes; ++j)
			{
				m_ddv[j].setZero();
				m_residual[j].setZero();
			}
		}
		updateVelocity();
	}
	if (perf.active)
	{
		const double totalMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - perfStart).count();
		perf.active = false; // Exclude stderr I/O from all measured intervals.
		int tetras = 0;
		for (int b = 0; b < m_softBodies.size(); ++b) tetras += m_softBodies[b]->m_tetras.size();
		// Each timing is milliseconds/calls; nested timings must not be summed.
		fprintf(stderr, "[deform-perf] solver=%p step=%llu inplace=%d cached=%d dt=%.6g nodes=%d tetras=%d constraints=%d projection=%d linesearch=%d total_ms=%.3f newton=%d krylov_iters=%d trials=%d line_failures=%d recovery=%d linear=%.3f/%d multiply=%.3f/%d combined=%.3f/%d damping=%.3f/%d elastic=%.3f/%d precondition=%.3f/%d blocks=%.3f/%d translation_setup=%.3f/%d translation_correct=%.3f/%d verify=%.3f/%d state=%.3f/%d energy=%.3f/%d residual=%.3f/%d cache_setup=%.3f/%d\n",
			static_cast<void*>(this), m_performanceStep, int(BT_KRYLOV_USE_INPLACE_UPDATES), int(BT_DEFORMABLE_USE_CACHED_OPERATOR), double(solverdt), m_numNodes, tetras, m_objective->m_projection.m_lagrangeMultipliers.size(), int(m_useProjection), int(m_lineSearch), totalMs, perf.newton, perf.krylovIterations, perf.lineTrials, perf.lineFailures, int(m_implicitRecoveryUsed),
			perf.linear.ms, perf.linear.calls,
			perf.multiply.ms, perf.multiply.calls,
			perf.combined.ms, perf.combined.calls,
			perf.damping.ms, perf.damping.calls,
			perf.elastic.ms, perf.elastic.calls,
			perf.precondition.ms, perf.precondition.calls,
			perf.blockSetup.ms, perf.blockSetup.calls,
			perf.translationSetup.ms, perf.translationSetup.calls,
			perf.translationCorrect.ms, perf.translationCorrect.calls,
			perf.verification.ms, perf.verification.calls,
			perf.state.ms, perf.state.calls,
			perf.energy.ms, perf.energy.calls,
			perf.residual.ms, perf.residual.calls, perf.cacheSetup.ms, perf.cacheSetup.calls);
		for (int i = 0; i < perf.newtonRecords.size(); ++i)
		{
			const auto& nr = perf.newtonRecords[i];
			fprintf(stderr, "[deform-newton] solver=%p step=%llu newton=%d outcome=%s solved=%d input_force_l2=%.9g input_constraint_inf=%.9g correction_l2=%.9g correction_tolerance=%.9g stationarity_l2=%.9g linear_momentum_l2=%.9g linear_constraint_inf=%.9g applied_scale=%.9g output_constraint_inf=%.9g\n",
				static_cast<void*>(this), m_performanceStep, nr.iteration, nr.outcome, int(nr.solved), nr.inputForceL2, nr.inputConstraintInf, nr.stepL2, nr.stepTolerance, nr.stationarityL2, nr.linearMomentumL2, nr.linearConstraintInf, nr.appliedScale, nr.outputConstraintInf);
		}
		for (int i = 0; i < perf.linearRecords.size(); ++i)
		{
			const btDeformableSolverPerformance::LinearRecord& record = perf.linearRecords[i];
			// Weighted CR values precede translation correction. Physical values
			// reuse the existing verification AFTER that correction, with no extra A*x.
			fprintf(stderr, "[deform-linear] solver=%p step=%llu newton=%d recovery=%d translation=%d relative=%d iterations=%d budget=%d stop=%s linear_ms=%.3f weighted_initial=%.9g weighted_final=%.9g weighted_target=%.9g physical_target=%.9g post_translation_momentum_l2=%.9g post_translation_constraint_inf=%.9g stationarity_l2=%.9g rhs_momentum_l2=%.9g rhs_constraint_inf=%.9g step_l2=%.9g checkpoint_ms=%.3f\n",
				static_cast<void*>(this), m_performanceStep, record.newton, int(record.recovery), int(record.translation), int(record.relative), record.iterations, record.budget, record.stop, record.ms, record.initial, record.final, record.target, record.physicalTarget, record.momentum, record.constraint, record.stationarity, record.rhsMomentum, record.rhsConstraint, record.stepL2, record.checkpointMs);
			for (int k = record.progressBegin; k < record.progressEnd; ++k)
			{
				const auto& progress = perf.progressRecords[k];
				fprintf(stderr, "[deform-cr-progress] solver=%p step=%llu newton=%d recovery=%d iteration=%d recurrence_weighted=%.9g best_weighted=%.9g recurrence_l2=%.9g true_momentum_l2=%.9g true_constraint_inf=%.9g recurrence_gap_l2=%.9g\n",
					static_cast<void*>(this), m_performanceStep, record.newton, int(record.recovery), progress.iteration,
					progress.recurrenceWeighted, progress.weightedBest, progress.recurrenceL2, progress.trueMomentumL2, progress.trueConstraintInf, progress.residualGapL2);
			}
		}
	}
}

btScalar btDeformableBodySolver::kineticEnergy()
{
	btScalar ke = 0;
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		for (int j = 0; j < psb->m_nodes.size(); ++j)
		{
			btSoftBody::Node& node = psb->m_nodes[j];
			if (node.m_frozen <= 0 && node.m_im > 0)
			{
				ke += m_dv[node.index].length2() * 0.5 / node.m_im;
			}
		}
	}
	return ke;
}

void btDeformableBodySolver::backupDv()
{
	m_backup_dv.resize(m_dv.size());
	for (int i = 0; i < m_backup_dv.size(); ++i)
	{
		m_backup_dv[i] = m_dv[i];
	}
}

void btDeformableBodySolver::revertDv()
{
	for (int i = 0; i < m_backup_dv.size(); ++i)
	{
		m_dv[i] = m_backup_dv[i];
	}
}

void btDeformableBodySolver::updateEnergy(btScalar scale)
{
	for (int i = 0; i < m_dv.size(); ++i)
	{
		m_dv[i] = m_backup_dv[i] + scale * m_ddv[i];
	}
	updateState();
}

btScalar btDeformableBodySolver::computeDescentStep(TVStack& ddv, const TVStack& residual, bool verbose)
{
	btScalar inner_product = 0;
	if (m_useProjection)
	{
		{
			btDeformablePerformanceScope linearPerf(m_objective->m_performance.linear, m_objective->m_performance.active);
			ImplicitOperatorCacheScope cacheScope(*m_objective);
			m_objective->m_performance.krylovIterations += m_cg.solve(*m_objective, ddv, residual, false);
		}
		inner_product = m_cg.dot(residual, m_ddv);
	}
	else
	{
		TVStack rhs, x;
		m_objective->addLagrangeMultiplierRHS(residual, m_dv, rhs);
		m_objective->addLagrangeMultiplier(ddv, x);
		solveImplicitKKT(x, rhs);
		for (int i = 0; i < ddv.size(); ++i)
		{
			ddv[i] = x[i];
		}
		inner_product = m_cg.dot(residual, ddv);
	}
	btScalar res_norm = m_objective->computeNorm(residual);
	btScalar tol = 1e-5 * res_norm * m_objective->computeNorm(m_ddv);
	if (inner_product < -tol)
	{
		if (verbose)
		{
			std::cout << "Looking backwards!" << std::endl;
		}
		for (int i = 0; i < m_ddv.size(); ++i)
		{
			m_ddv[i] = -m_ddv[i];
		}
		inner_product = -inner_product;
	}
	else if (std::abs(inner_product) < tol)
	{
		if (verbose)
		{
			std::cout << "Gradient Descent!" << std::endl;
		}
		btScalar scale = m_objective->computeNorm(m_ddv) / res_norm;
		for (int i = 0; i < m_ddv.size(); ++i)
		{
			m_ddv[i] = scale * residual[i];
		}
		inner_product = scale * res_norm * res_norm;
	}
	return inner_product;
}

void btDeformableBodySolver::updateState()
{
	btDeformablePerformanceScope perf(m_objective->m_performance.state, m_objective->m_performance.active);
	updateVelocity();
	updateTempPosition();
}

void btDeformableBodySolver::updateDv(btScalar scale)
{
	for (int i = 0; i < m_numNodes; ++i)
	{
		m_dv[i] += scale * m_ddv[i];
	}
}


btScalar btDeformableBodySolver::implicitConstraintError(const TVStack& dv) const
{
	btScalar error = 0;
	const btAlignedObjectArray<LagrangeMultiplier>& multipliers = m_objective->m_projection.m_lagrangeMultipliers;
	for (int c = 0; c < multipliers.size(); ++c)
	{
		const LagrangeMultiplier& lm = multipliers[c];
		for (int d = 0; d < lm.m_num_constraints; ++d)
		{
			btScalar value = 0;
			for (int n = 0; n < lm.m_num_nodes; ++n)
				{
				btVector3 difference = dv[lm.m_indices[n]];
				if (m_objective->m_implicitConstraintDv.size() == dv.size())
					difference -= m_objective->m_implicitConstraintDv[lm.m_indices[n]];
				value += lm.m_weights[n] * difference.dot(lm.m_dirs[d]);
			}
			if (!std::isfinite((double)value)) return SIMD_INFINITY;
			error = btMax(error, btFabs(value));
		}
	}
	return error;
}

void btDeformableBodySolver::measureImplicitLinearResidual(const TVStack& x, const TVStack& rhs)
{
	btDeformablePerformanceScope perf(m_objective->m_performance.verification, m_objective->m_performance.active);
	TVStack product;
	product.resize(rhs.size());
	m_objective->multiply(x, product);
	btScalar momentumSquared = 0;
	m_lastLinearConstraintResidual = 0;
	for (int i = 0; i < rhs.size(); ++i)
	{
		const btVector3 error = rhs[i] - product[i];
		if (!std::isfinite((double)error.x()) || !std::isfinite((double)error.y()) || !std::isfinite((double)error.z()))
		{
			m_lastStationarityResidual = m_lastLinearMomentumResidual = m_lastLinearConstraintResidual = SIMD_INFINITY;
			return;
		}
		if (i < m_numNodes)
			momentumSquared += error.length2();
		else
			for (int d = 0; d < 3; ++d)
				m_lastLinearConstraintResidual = btMax(m_lastLinearConstraintResidual, btFabs(error[d]));
	}
	m_lastLinearMomentumResidual = btSqrt(momentumSquared);
	// Momentum residual at the CURRENT nonlinear state, allowing reaction
	// forces C^T*lambda. Do not include A*ddv: a small ddv can still have a
	// large effect in a stiff system, so it must not justify discarding ddv.
	TVStack stationarity = rhs;
	const btAlignedObjectArray<LagrangeMultiplier>& multipliers = m_objective->m_projection.m_lagrangeMultipliers;
	for (int c = 0; c < multipliers.size(); ++c)
	{
		const LagrangeMultiplier& lm = multipliers[c];
		for (int n = 0; n < lm.m_num_nodes; ++n)
			for (int d = 0; d < lm.m_num_constraints; ++d)
				stationarity[lm.m_indices[n]] -= lm.m_weights[n] * x[m_numNodes + c][d] * lm.m_dirs[d];
	}
	btScalar stationaritySquared = 0;
	for (int n = 0; n < m_numNodes; ++n) stationaritySquared += stationarity[n].length2();
	m_lastStationarityResidual = std::isfinite((double)stationaritySquared) ? btSqrt(stationaritySquared) : SIMD_INFINITY;
}

void btDeformableBodySolver::solveImplicitKKT(TVStack& x, const TVStack& rhs)
{
	ImplicitOperatorCacheScope cacheScope(*m_objective);
	{
		btDeformablePerformanceScope blockPerf(m_objective->m_performance.blockSetup, m_objective->m_performance.active);
		m_objective->m_preconditioner->reinitialize(true);
	}
	const bool translation = m_objective->setupTranslationCorrection();
	m_cr.configureProgress(m_objective->m_performance.active &&
		m_objective->m_projection.m_lagrangeMultipliers.size() > 0, m_numNodes);
	m_objective->correctTranslation(x, rhs);
	// Preserve the cheap first Newton solve. Continue difficult free-body
	// corrections without discarding conjugate directions at iteration 300.
	const bool continued = translation && m_objective->m_projection.m_lagrangeMultipliers.size() == 0 && m_newtonIteration > 0;
	const int linearBudget = continued ? 1200 : kMaxConjugateGradientIterations;
	const btScalar physicalTarget = continued ? btScalar(0.5) * m_newtonTolerance : btScalar(0);
	const auto recordSolve = [&](int iterations, int budget, bool recovery, bool relative, btScalar target, double ms)
	{
		if (!m_objective->m_performance.active) return;
		btDeformableSolverPerformance::LinearRecord record;
		record.newton = m_newtonIteration + 1;
		record.iterations = iterations; record.budget = budget;
		record.recovery = recovery; record.translation = translation; record.relative = relative;
		record.stop = m_cr.getStopReason(); record.ms = ms;
		record.initial = m_cr.getInitialResidual(); record.final = m_cr.getFinalResidual();
		record.target = m_cr.getTargetResidual(); record.physicalTarget = target;
		record.momentum = m_lastLinearMomentumResidual; record.constraint = m_lastLinearConstraintResidual;
		record.stationarity = m_lastStationarityResidual;
		double rhs2 = 0, step2 = 0, constraint = 0;
		for (int n = 0; n < rhs.size(); ++n)
		{
			if (n < m_numNodes) { rhs2 += double(rhs[n].length2()); step2 += double(x[n].length2()); }
			else for (int d = 0; d < 3; ++d) constraint = btMax(constraint, double(btFabs(rhs[n][d])));
		}
		record.rhsMomentum = std::sqrt(rhs2); record.rhsConstraint = constraint;
		record.stepL2 = std::sqrt(step2); record.checkpointMs = m_cr.m_progressMs;
		record.progressBegin = m_objective->m_performance.progressRecords.size();
		for (int k = 0; k < m_cr.m_progress.size(); ++k)
		{
			const auto& source = m_cr.m_progress[k];
			btDeformableSolverPerformance::ProgressRecord entry;
			entry.iteration = source.iteration; entry.recurrenceWeighted = source.recurrenceWeighted;
			entry.recurrenceL2 = source.recurrenceL2; entry.trueMomentumL2 = source.trueMomentumL2;
			entry.trueConstraintInf = source.trueConstraintInf; entry.residualGapL2 = source.residualGapL2;
			entry.weightedBest = source.weightedBest;
			m_objective->m_performance.progressRecords.push_back(entry);
		}
		record.progressEnd = m_objective->m_performance.progressRecords.size();
		m_objective->m_performance.linearRecords.push_back(record);
	};
	int iterations;
	double linearStartMs = m_objective->m_performance.linear.ms;
	{
		btDeformablePerformanceScope linearPerf(m_objective->m_performance.linear, m_objective->m_performance.active);
		iterations = m_cr.solveWithConvergencePolicy(*m_objective, x, rhs, false, !continued, true, linearBudget, physicalTarget);
		m_objective->m_performance.krylovIterations += iterations;
	}
	m_objective->correctTranslation(x, rhs);
	measureImplicitLinearResidual(x, rhs);
	recordSolve(iterations, linearBudget, false, !continued, physicalTarget, m_objective->m_performance.linear.ms - linearStartMs);
	btScalar stepSquared = 0;
	for (int n = 0; n < m_numNodes; ++n) stepSquared += x[n].length2();
	const btScalar stepTolerance = m_newtonTolerance * btMax(btScalar(1), m_objective->computeNorm(m_dv));
	if (!continued && !m_implicitRecoveryUsed && btSqrt(stepSquared) <= stepTolerance &&
		(m_lastLinearMomentumResidual > m_newtonTolerance || m_lastLinearConstraintResidual > m_newtonTolerance))
	{
		// One full-budget restart per timestep for an inaccurate, tiny step.
		m_implicitRecoveryUsed = true;
		linearStartMs = m_objective->m_performance.linear.ms;
		{
			btDeformablePerformanceScope linearPerf(m_objective->m_performance.linear, m_objective->m_performance.active);
			iterations = m_cr.solveWithConvergencePolicy(*m_objective, x, rhs, false, false, true);
			m_objective->m_performance.krylovIterations += iterations;
		}
		m_objective->correctTranslation(x, rhs);
		measureImplicitLinearResidual(x, rhs);
		recordSolve(iterations, kMaxConjugateGradientIterations, true, false, 0, m_objective->m_performance.linear.ms - linearStartMs);
	}
	m_cr.configureProgress(false, 0);
	m_objective->m_translationCorrection = false;
}

void btDeformableBodySolver::computeStep(TVStack& ddv, const TVStack& residual)
{
	if (m_useProjection)
	{
		{
			btDeformablePerformanceScope linearPerf(m_objective->m_performance.linear, m_objective->m_performance.active);
			ImplicitOperatorCacheScope cacheScope(*m_objective);
			m_objective->m_performance.krylovIterations += m_cg.solve(*m_objective, ddv, residual, false);
		}
	}
	else
	{
		TVStack rhs, x;
		m_objective->addLagrangeMultiplierRHS(residual, m_dv, rhs);
		m_objective->addLagrangeMultiplier(ddv, x);
		solveImplicitKKT(x, rhs);
		for (int i = 0; i < ddv.size(); ++i)
		{
			ddv[i] = x[i];
		}
	}
}

void btDeformableBodySolver::reinitialize(const btAlignedObjectArray<btSoftBody*>& softBodies, btScalar dt)
{
	m_softBodies.copyFromArray(softBodies);
	bool nodeUpdated = updateNodes();

	if (nodeUpdated)
	{
		m_dv.resize(m_numNodes, btVector3(0, 0, 0));
		m_ddv.resize(m_numNodes, btVector3(0, 0, 0));
		m_residual.resize(m_numNodes, btVector3(0, 0, 0));
		m_backupVelocity.resize(m_numNodes, btVector3(0, 0, 0));
	}

	// need to setZero here as resize only set value for newly allocated items
	for (int i = 0; i < m_numNodes; ++i)
	{
		m_dv[i].setZero();
		m_ddv[i].setZero();
		m_residual[i].setZero();
	}

	if (dt > 0)
	{
		m_dt = dt;
	}
	m_objective->reinitialize(nodeUpdated, dt);
	updateSoftBodies();
}

void btDeformableBodySolver::setConstraints(const btContactSolverInfo& infoGlobal)
{
	BT_PROFILE("setConstraint");
	m_objective->setConstraints(infoGlobal);
}

btScalar btDeformableBodySolver::solveContactConstraints(btCollisionObject** deformableBodies, int numDeformableBodies, const btContactSolverInfo& infoGlobal)
{
	BT_PROFILE("solveContactConstraints");
	btScalar maxSquaredResidual = m_objective->m_projection.update(deformableBodies, numDeformableBodies, infoGlobal);
	return maxSquaredResidual;
}

void btDeformableBodySolver::updateVelocity()
{
	int counter = 0;
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		psb->m_maxSpeedSquared = 0;
		if (!psb->isActive() || psb->isStaticObject())
		{
			counter += psb->m_nodes.size();
			continue;
		}
		for (int j = 0; j < psb->m_nodes.size(); ++j)
		{
			// set NaN to zero;
			if (m_dv[counter] != m_dv[counter])
			{
				m_dv[counter].setZero();
			}
			if (m_implicit)
			{
				psb->m_nodes[j].m_v = m_backupVelocity[counter] + m_dv[counter];
			}
			else
			{
				psb->m_nodes[j].m_v = m_backupVelocity[counter] + m_dv[counter] - psb->m_nodes[j].m_splitv;
			}
			psb->m_maxSpeedSquared = btMax(psb->m_maxSpeedSquared, psb->m_nodes[j].m_v.length2());
			++counter;
		}
	}
}

void btDeformableBodySolver::updateTempPosition()
{
	int counter = 0;
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		if (!psb->isActive() || psb->isStaticObject())
		{
			counter += psb->m_nodes.size();
			continue;
		}
		for (int j = 0; j < psb->m_nodes.size(); ++j)
		{
			psb->m_nodes[j].m_q = psb->m_nodes[j].m_x + m_dt * (psb->m_nodes[j].m_v + psb->m_nodes[j].m_splitv);
			++counter;
		}
		psb->updateDeformation();
	}
}

void btDeformableBodySolver::backupVelocity()
{
	int counter = 0;
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		for (int j = 0; j < psb->m_nodes.size(); ++j)
		{
			m_backupVelocity[counter++] = psb->m_nodes[j].m_v;
		}
	}
}

void btDeformableBodySolver::setupDeformableSolve(bool implicit)
{
	int counter = 0;
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		if (!psb->isActive() || psb->isStaticObject())
		{
			counter += psb->m_nodes.size();
			continue;
		}
		for (int j = 0; j < psb->m_nodes.size(); ++j)
		{
			if (implicit)
			{
				// Use the current post-explicit/post-contact velocity as the initial guess for the implicit solve.
				// Zeroing unconstrained nodes here throws away the explicit gravity update in free flight.
				m_dv[counter] = psb->m_nodes[j].m_v - psb->m_nodes[j].m_vn;
				m_backupVelocity[counter] = psb->m_nodes[j].m_vn;
			}
			else
			{
				m_dv[counter] = psb->m_nodes[j].m_v + psb->m_nodes[j].m_splitv - m_backupVelocity[counter];
			}
			psb->m_nodes[j].m_v = m_backupVelocity[counter];
			++counter;
		}
	}
	// Preserve the constrained components established by the contact solver.
	// Newton changes deformation velocities, not the contact's target velocity.
	if (implicit) m_objective->m_implicitConstraintDv = m_dv;
	else m_objective->m_implicitConstraintDv.clear();

}

void btDeformableBodySolver::revertVelocity()
{
	int counter = 0;
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		for (int j = 0; j < psb->m_nodes.size(); ++j)
		{
			psb->m_nodes[j].m_v = m_backupVelocity[counter++];
		}
	}
}

bool btDeformableBodySolver::updateNodes()
{
	int numNodes = 0;
	for (int i = 0; i < m_softBodies.size(); ++i)
		numNodes += m_softBodies[i]->m_nodes.size();
	if (numNodes != m_numNodes)
	{
		m_numNodes = numNodes;
		return true;
	}
	return false;
}

void btDeformableBodySolver::predictMotion(btScalar solverdt)
{
	// apply explicit forces to velocity
	if (m_implicit)
	{
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (psb->isActive() && !psb->isStaticObject())
			{
				for (int j = 0; j < psb->m_nodes.size(); ++j)
				{
					psb->m_nodes[j].m_q = psb->m_nodes[j].m_x + psb->m_nodes[j].m_v * solverdt;
				}
			}
		}
	}
	applyExplicitForce();
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		if (psb->isActive() && !psb->isStaticObject())
		{
			/* Clear contacts when softbody is active*/
			psb->m_nodeRigidContacts.resize(0);
			psb->m_faceRigidContacts.resize(0);
			psb->m_faceNodeContacts.resize(0);
			psb->m_nodeNodeContacts.resize(0);
			psb->m_faceNodeContactsCCD.resize(0);
			// predict motion for collision detection
			predictDeformableMotion(psb, solverdt);
		}
	}
}

void btDeformableBodySolver::predictDeformableMotion(btSoftBody* psb, btScalar dt)
{
	BT_PROFILE("btDeformableBodySolver::predictDeformableMotion");
	int i, ni;

	/* Update                */
	if (psb->m_bUpdateRtCst)
	{
		psb->m_bUpdateRtCst = false;
		psb->updateConstants();
		psb->m_fdbvt.clear();
		if (psb->m_cfg.collisions & btSoftBody::fCollision::SDF_RD)
		{
			psb->initializeFaceTree();
		}
	}

	/* Prepare                */
	psb->m_sst.sdt = dt * psb->m_cfg.timescale;
	psb->m_sst.isdt = 1 / psb->m_sst.sdt;
	psb->m_sst.velmrg = psb->m_sst.sdt * 3;
	psb->m_sst.radmrg = psb->getCollisionShape()->getMargin();
	psb->m_sst.updmrg = psb->m_sst.radmrg * (btScalar)0.25;
	/* Bounds                */
	psb->updateBounds();

	/* Integrate            */
	// do not allow particles to move more than the bounding box size
	btScalar max_v = (psb->m_bounds[1] - psb->m_bounds[0]).norm() / dt;
	for (i = 0, ni = psb->m_nodes.size(); i < ni; ++i)
	{
		btSoftBody::Node& n = psb->m_nodes[i];
		// apply drag
		n.m_v *= (1 - psb->m_cfg.drag);
		// scale velocity back
		if (m_implicit)
		{
			n.m_q = n.m_x;
		}
		else
		{
			if (n.m_v.norm() > max_v)
			{
				n.m_v.safeNormalize();
				n.m_v *= max_v;
			}
			n.m_q = n.m_x + n.m_v * dt;
		}
		n.m_splitv.setZero();
		n.m_constrained = false;
	}

	/* Nodes                */
	psb->updateNodeTree(true, true);
	if (!psb->m_fdbvt.empty())
	{
		psb->updateFaceTree(true, true);
	}
	/* Optimize dbvt's        */
	//    psb->m_ndbvt.optimizeIncremental(1);
	//    psb->m_fdbvt.optimizeIncremental(1);
}

void btDeformableBodySolver::updateSoftBodies()
{
	BT_PROFILE("updateSoftBodies");
	for (int i = 0; i < m_softBodies.size(); i++)
	{
		btSoftBody* psb = (btSoftBody*)m_softBodies[i];
		if (psb->isActive() && !psb->isStaticObject())
		{
			psb->updateNormals();
		}
	}
}

void btDeformableBodySolver::setImplicit(bool implicit)
{
	m_implicit = implicit;
	m_objective->setImplicit(implicit);
}

void btDeformableBodySolver::setLineSearch(bool lineSearch)
{
	m_lineSearch = lineSearch;
}

void btDeformableBodySolver::setMaxNewtonIterations(int maxNewtonIterations)
{
	m_maxNewtonIterations = btMax(1, maxNewtonIterations);
}

void btDeformableBodySolver::applyExplicitForce()
{
	m_objective->applyExplicitForce(m_residual);
}

void btDeformableBodySolver::applyTransforms(btScalar timeStep)
{
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		//fprintf(stderr, "framestart()\n");
		for (int j = 0; j < psb->m_nodes.size(); ++j)
		{
			btSoftBody::Node& node = psb->m_nodes[j];
			btScalar maxDisplacement = psb->getWorldInfo()->m_maxDisplacement;
			btScalar clampDeltaV = maxDisplacement / timeStep;
			//fprintf(stderr, "v %d %f %f %f\n", j, node.m_v.x(), node.m_v.y(), node.m_v.z());
			for (int c = 0; c < 3; c++)
			{
				if (node.m_v[c] > clampDeltaV)
				{
					node.m_v[c] = clampDeltaV;
				}
				if (node.m_v[c] < -clampDeltaV)
				{
					node.m_v[c] = -clampDeltaV;
				}
			}
			node.m_x = node.m_x + timeStep * (node.m_v + node.m_splitv);
			node.m_q = node.m_x;
			node.m_vn = node.m_v;
		}
		// enforce anchor constraints
		for (int j = 0; j < psb->m_deformableAnchors.size(); ++j)
		{
			btSoftBody::DeformableNodeRigidAnchor& a = psb->m_deformableAnchors[j];
			btSoftBody::Node* n = a.m_node;
			//fprintf(stderr, "a.m_local %d %f %f %f\n", j, a.m_local.x(), a.m_local.y(), a.m_local.z());
			//fprintf(stderr, "n->m_x orig %d %f %f %f\n", j, n->m_x.x(), n->m_x.y(), n->m_x.z());
			//fprintf(stderr, "getWorldTransform %d %f %f %f\n", j, a.m_cti.m_colObj->getWorldTransform().getOrigin().x(), a.m_cti.m_colObj->getWorldTransform().getOrigin().y(), a.m_cti.m_colObj->getWorldTransform().getOrigin().z());
			//n->m_x = a.m_cti.m_colObj->getWorldTransform() * a.m_local;
			//fprintf(stderr, "drawpoint \"pt\" [%f,%f,%f]\n", j, n->m_x.x(), n->m_x.y(), n->m_x.z());

			// update multibody anchor info
			if (a.m_cti.m_colObj->getInternalType() == btCollisionObject::CO_FEATHERSTONE_LINK)
			{
				btMultiBodyLinkCollider* multibodyLinkCol = (btMultiBodyLinkCollider*)btMultiBodyLinkCollider::upcast(a.m_cti.m_colObj);
				if (multibodyLinkCol)
				{
					btVector3 nrm;
					const btCollisionShape* shp = multibodyLinkCol->getCollisionShape();
					const btTransform& wtr = multibodyLinkCol->getWorldTransform();
					psb->m_worldInfo->m_sparsesdf.Evaluate(
						wtr.invXform(n->m_x),
						shp,
						nrm,
						0);
					a.m_cti.m_normal = wtr.getBasis() * nrm;
					btVector3 normal = a.m_cti.m_normal;
					btVector3 t1 = generateUnitOrthogonalVector(normal);
					btVector3 t2 = btCross(normal, t1);
					btMultiBodyJacobianData jacobianData_normal, jacobianData_t1, jacobianData_t2;
					findJacobian(multibodyLinkCol, jacobianData_normal, a.m_node->m_x, normal);
					findJacobian(multibodyLinkCol, jacobianData_t1, a.m_node->m_x, t1);
					findJacobian(multibodyLinkCol, jacobianData_t2, a.m_node->m_x, t2);

					btScalar* J_n = &jacobianData_normal.m_jacobians[0];
					btScalar* J_t1 = &jacobianData_t1.m_jacobians[0];
					btScalar* J_t2 = &jacobianData_t2.m_jacobians[0];

					btScalar* u_n = &jacobianData_normal.m_deltaVelocitiesUnitImpulse[0];
					btScalar* u_t1 = &jacobianData_t1.m_deltaVelocitiesUnitImpulse[0];
					btScalar* u_t2 = &jacobianData_t2.m_deltaVelocitiesUnitImpulse[0];

					btMatrix3x3 rot(normal.getX(), normal.getY(), normal.getZ(),
									t1.getX(), t1.getY(), t1.getZ(),
									t2.getX(), t2.getY(), t2.getZ());  // world frame to local frame
					const int ndof = multibodyLinkCol->m_multiBody->getNumDofs() + 6;
					btMatrix3x3 local_impulse_matrix = (Diagonal(n->m_im) + OuterProduct(J_n, J_t1, J_t2, u_n, u_t1, u_t2, ndof)).inverse();
					a.m_c0 = rot.transpose() * local_impulse_matrix * rot;
					a.jacobianData_normal = jacobianData_normal;
					a.jacobianData_t1 = jacobianData_t1;
					a.jacobianData_t2 = jacobianData_t2;
					a.t1 = t1;
					a.t2 = t2;
				}
			}
		}
		//fprintf(stderr, "frameend()\n");
		psb->interpolateRenderMesh();
	}
}

void btDeformableBodySolver::processCollision(btSoftBody* softBody, const btCollisionObjectWrapper* collisionObjectWrap, btManifoldResultForSkin* resultOut)
{
	if (softBody->getCollisionShape()->getShapeType() == SOFTBODY_SHAPE_PROXYTYPE)
		softBody->defaultCollisionHandler(collisionObjectWrap);
	else
	{
		resultOut->getPersistentManifold()->m_responseProcessedEarly = true;
		auto& cp = resultOut->getPersistentManifold()->getContactPoint(resultOut->contactIndex);
		softBody->skinSoftRigidCollisionHandler(collisionObjectWrap, resultOut->getPartId0(), resultOut->getIndex0(),
												resultOut->swapped ? (/*not using cp.getPositionWorldOnA on purpose because it is calculated using wrong depth at the moment.
                                                    See the comment in btManifoldResult::addContactPoint (the one which starts "Ideally there should be this commented out...") */
																	  cp.getPositionWorldOnB() - cp.m_normalWorldOnB * cp.getUnmodifiedDistance())
																   : cp.getPositionWorldOnB(),  // Not sure that this is correct. I am sure that I have seen it swapped once, but was not able to reproduce it since.
												resultOut->swapped ? cp.m_normalWorldOnB : -cp.m_normalWorldOnB,
												cp.getDistance(), cp.m_contactPointFlags & BT_CONTACT_FLAG_PENETRATING, &cp.m_appliedImpulse);
	}
}

void btDeformableBodySolver::processCollision(btSoftBody* softBody, btSoftBody* otherSoftBody)
{
	softBody->defaultCollisionHandler(otherSoftBody);
}

void btDeformableBodySolver::processCollision(btSoftBody* softBody, btSoftBody* otherSoftBody, btManifoldResultForSkin* resultOut)
{
	resultOut->getPersistentManifold()->m_responseProcessedEarly = true;
	auto& cp = resultOut->getPersistentManifold()->getContactPoint(resultOut->contactIndex);
	auto contactPoint = resultOut->swapped ? (/*not using cp.getPositionWorldOnA on purpose because it is calculated using wrong depth at the moment.
                                                    See the comment in btManifoldResult::addContactPoint (the one which starts "Ideally there should be this commented out...") */
											  cp.getPositionWorldOnB() - cp.m_normalWorldOnB * cp.getUnmodifiedDistance())
										   : cp.getPositionWorldOnB();  // Not sure that this is correct. I am sure that I have seen it swapped once, but was not able to reproduce it since.
	auto normal = resultOut->swapped ? cp.m_normalWorldOnB : -cp.m_normalWorldOnB;
	softBody->skinSoftSoftCollisionHandler(otherSoftBody, resultOut->getPartId0(), resultOut->getIndex0(), resultOut->getPartId1(), resultOut->getIndex1(), contactPoint, normal, cp.getDistance(), cp.m_contactPointFlags & BT_CONTACT_FLAG_PENETRATING,
										   &cp.m_appliedImpulse);
}
