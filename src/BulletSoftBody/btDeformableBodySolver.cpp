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
#include "btDeformableBodySolver.h"
#include "btSoftBodyInternals.h"
#include "btDeformableDiagnostics.h"
#include "btDeformableVolumeBarrierForce.h"
#include "btDeformableContactForce.h"
#include "btDeformableNewtonSnapshot.h"
#include "btDeformableEnergyChange.h"
#include <string>
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
		: m_objective(objective), m_enabled(objective.m_implicit)
	{
		if (!m_enabled) return;
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
	m_lastSolveConverged = false;
	m_lastSolveInvalidPredictor = false;
	for (int j = 0; j < m_numNodes; ++j) { m_ddv[j].setZero(); m_residual[j].setZero(); }
	const bool diagnose = btDeformableDiagnostics::enabled();
	btDeformableVolumeBarrierForce* barrier = nullptr;
	m_contactWeightedTarget = 0;
	if (m_implicit && !m_useProjection && m_objective->m_preconditioner == m_objective->m_KKTPreconditioner &&
		!m_objective->m_projection.m_lagrangeMultipliers.size())
		for (int f = 0; f < m_objective->m_lf.size(); ++f)
			if (m_objective->m_lf[f]->getForceType() == BT_CONTACT_FORCE)
			{
				const auto& contacts = static_cast<btDeformableContactForce*>(m_objective->m_lf[f])->contacts;
				for (int c = 0; c < contacts.size(); ++c)
				{
					// rho = 10/compliance. Reserve a tenth of the outer velocity tolerance.
					const btScalar target = btScalar(0.0001) * btSqrt(btMin(contacts[c].rho, contacts[c].tangentRho) / 10);
					if (target > 0 && (m_contactWeightedTarget == 0 || target < m_contactWeightedTarget)) m_contactWeightedTarget = target;
				}
			}
	for (int f = 0; f < m_objective->m_lf.size(); ++f)
		if (m_objective->m_lf[f]->getForceType() == BT_VOLUME_BARRIER_FORCE)
			barrier = static_cast<btDeformableVolumeBarrierForce*>(m_objective->m_lf[f]);
	if (diagnose)
	{
		btDeformableDiagnostics::write("SOLVER", "implicit=%d projection=%d line_search=%d max_newton=%d tolerance=%.9g multipliers=%d",
			int(m_implicit), int(m_useProjection), int(m_lineSearch), m_maxNewtonIterations, double(m_newtonTolerance), m_objective->m_projection.m_lagrangeMultipliers.size());
		if (btDeformableDiagnostics::current()->step % 1000 == 0)
			for (int f = 0; f < m_objective->m_lf.size(); ++f)
			{
				auto& force = *m_objective->m_lf[f];
				for (int b = 0; b < force.m_softBodies.size(); ++b)
					btDeformableDiagnostics::write("MATERIAL", "id=%d force_type=%d young=%.9g poisson=%.9g",
						force.m_softBodies[b]->getUserIndex(), int(force.getForceType()), double(force.getYoungsModulus()), double(force.getPoissonRatio()));
			}
	}
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
			m_objective->m_preconditioner->reinitialize(true);
			// Explicit damping needs the previous absolute-tolerance/full-budget
			// policy: a small relative residual can still leave large velocity errors.
			m_cr.solveWithConvergencePolicy(*m_objective, x, rhs, false, false);
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
		const char* diagnosticExit = "iteration_limit";
		int diagnosticIterations = 0;
		for (int i = 0; i < m_maxNewtonIterations; ++i)
		{
			diagnosticIterations = i + 1;
			m_newtonIteration = i;
			updateState();
			if (barrier && !barrier->admissible())
			{
				m_lastSolveInvalidPredictor = true;
				diagnosticExit = "volume_predictor_invalid";
				break;
			}
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
			bool contactAccurate = true;
			if (m_contactWeightedTarget > 0)
			{
				m_objective->m_KKTPreconditioner->reinitialize(true);
				btScalar squared = 0;
				for (int n = 0; n < m_numNodes; ++n)
					squared += m_residual[n].dot(m_objective->m_KKTPreconditioner->applyInverseNodeBlock(n, m_residual[n]));
				const btScalar weighted = squared >= 0 && std::isfinite(double(squared)) ? btSqrt(squared) : SIMD_INFINITY;
				contactAccurate = weighted <= m_contactWeightedTarget;
				if (diagnose) btDeformableDiagnostics::write("CONTACT_ACCURACY", "iteration=%d weighted_residual=%.9g target=%.9g accurate=%d",
					i, double(weighted), double(m_contactWeightedTarget), int(contactAccurate));
			}
			if (diagnose) btDeformableDiagnostics::write("NEWTON_STATE", "iteration=%d residual=%.9g constraint_error=%.9g",
				i, double(residualNorm), double(constraintError));
			if (contactAccurate && residualNorm < m_newtonTolerance && constraintError <= m_newtonTolerance && i > 0)
			{
				m_lastSolveConverged = true;
				diagnosticExit = "residual_tolerance";
				break;
			}
			// todo xuchenhan@: this really only needs to be calculated once
			m_objective->applyDynamicFriction(m_residual);
			if (m_lineSearch)
			{
				const bool linearRecoveryBefore = m_implicitRecoveryUsed;
				btScalar inner_product;
				if (barrier) { computeStep(m_ddv, m_residual); inner_product = m_cg.dot(m_residual, m_ddv); }
				else inner_product = computeDescentStep(m_ddv, m_residual);
				const btScalar stepNorm = m_objective->computeNorm(m_ddv);
				const btScalar relativeStepTolerance = m_newtonTolerance * btMax(btScalar(1), m_objective->computeNorm(m_dv));
				const bool accurateLinearStep = !m_useProjection &&
					m_lastLinearMomentumResidual <= m_newtonTolerance &&
					m_lastLinearConstraintResidual <= m_newtonTolerance;
				const bool stepConverged = contactAccurate && accurateLinearStep &&
					m_lastStationarityResidual <= m_newtonTolerance && constraintError <= m_newtonTolerance;
				if (diagnose) btDeformableDiagnostics::write("NEWTON_LINEAR", "iteration=%d step_norm=%.9g step_tolerance=%.9g accurate=%d converged=%d momentum_residual=%.9g constraint_residual=%.9g stationarity=%.9g recovery=%d",
					i, double(stepNorm), double(relativeStepTolerance), int(accurateLinearStep), int(stepConverged), double(m_lastLinearMomentumResidual), double(m_lastLinearConstraintResidual), double(m_lastStationarityResidual), int(m_implicitRecoveryUsed));
				if (i > 0 && stepNorm <= relativeStepTolerance && (stepConverged || !accurateLinearStep))
				{
					m_lastSolveConverged = stepConverged;
					diagnosticExit = stepConverged ? "small_step_verified" : "small_step_unverified";
					// Check stationarity at the accepted state, not just the solved
					// linear system. Otherwise apply a small but accurate correction.
					break;
				}
				btScalar alpha = 0.01, beta = 0.5;  // Boyd & Vandenberghe suggested alpha between 0.01 and 0.3, beta between 0.1 to 0.8
				btScalar scale = 2 * (barrier ? barrier->safeStep(solverdt, m_ddv) : btScalar(1));
				const btScalar initialScale = scale * btScalar(.5);
				const bool captureCandidate = diagnose && btDeformableDiagnostics::current()->step != m_lastCapturedStep &&
					double(solverdt) <= btDeformableDiagnostics::current()->dt / 256 * 1.001;
				btAlignedObjectArray<btScalar> trialScales, trialEnergies;
				btScalar f0 = m_objective->totalEnergy(solverdt) + kineticEnergy(), f1, f2;
				// TODO: Remove all environment-variable lookups across this feature before merging the feature branch.
				// The final implementation must not read environment variables.
				const char* requestedCapture = std::getenv("BULLET_DEFORMABLE_CAPTURE_SOLVER_STEP");
				if (diagnose && requestedCapture && i == m_maxNewtonIterations - 1 &&
					btDeformableDiagnostics::current()->step != m_lastCapturedStep)
				{
					char* end = nullptr;
					const long long requestedStep = std::strtoll(requestedCapture, &end, 10);
					const char* logPath = std::getenv("BULLET_DEFORMABLE_DIAGNOSTICS");
					if (end != requestedCapture && *end == '\0' && requestedStep >= 0 &&
						requestedStep == btDeformableDiagnostics::current()->step && logPath && *logPath)
					{
						m_lastCapturedStep = requestedStep;
						const std::string path = std::string(logPath) + ".step-" + std::to_string(requestedStep) + ".newton.bin";
						btDeformableNewtonSnapshot snapshot;
						snapshot.dt=solverdt;snapshot.step=requestedStep;snapshot.baselineEnergy=f0;
						snapshot.slope=inner_product;snapshot.initialScale=initialScale;snapshot.linearRecoveryBefore=linearRecoveryBefore;
						const bool saved=snapshot.save(path.c_str(),*this);
						btDeformableDiagnostics::write("NEWTON_CAPTURE", "success=%d kind=iteration_limit_probe path=%s h=%.17g",
							int(saved), path.c_str(), double(solverdt));
					}
				}
				bool lineSearchFailed = false;
				bool sufficientDecrease = false;
				backupDv();
				const char* stableValue = std::getenv("BULLET_DEFORMABLE_STABLE_ENERGY");
				std::unique_ptr<btDeformableEnergyChange> energyChange;
				if (m_contactWeightedTarget > 0 && stableValue && stableValue[0] == '1')
					energyChange.reset(new btDeformableEnergyChange(*m_objective, m_dv, solverdt));
				do
				{
					scale *= beta;
					if (scale < 1e-8)
					{
						lineSearchFailed = true;
						break;
					}
					updateEnergy(scale);
					if (!energyChange) f1 = m_objective->totalEnergy(solverdt) + kineticEnergy();
					const double energyDelta = energyChange ? energyChange->difference(m_dv) : double(f1 - f0);
					if (energyChange) f1 = f0 + energyDelta;
					f2 = f0 - alpha * scale * inner_product;
					sufficientDecrease = energyChange ? energyDelta < -double(alpha * scale * inner_product) : f1 < f2 + SIMD_EPSILON;
					const btScalar roundoff = btScalar(8) * std::numeric_limits<btScalar>::epsilon() * btMax(btScalar(1), btFabs(f0));
					if (!sufficientDecrease && m_contactWeightedTarget > 0 && inner_product > 0 &&
						alpha * scale * inner_product <= roundoff && std::abs(energyDelta) <= roundoff)
					{
						// Absolute potential offsets can hide a valid decrease. In that
						// roundoff band, require progress in the fixed weighted residual.
						TVStack trialResidual; trialResidual.resize(m_numNodes, btVector3(0,0,0));
						for (int n = 0; n < m_numNodes; ++n)
						{
							const auto* node = m_objective->m_nodes[n];
							if (node->m_im > 0 && node->m_frozen <= 0) trialResidual[n] = -m_dv[n] / node->m_im;
						}
						m_objective->computeResidual(solverdt, trialResidual);
						btScalar before = 0, after = 0;
						for (int n = 0; n < m_numNodes; ++n)
						{
							before += m_residual[n].dot(m_objective->m_KKTPreconditioner->applyInverseNodeBlock(n, m_residual[n]));
							after += trialResidual[n].dot(m_objective->m_KKTPreconditioner->applyInverseNodeBlock(n, trialResidual[n]));
						}
						sufficientDecrease = std::isfinite(double(before)) && std::isfinite(double(after)) &&
							after >= 0 && after < before * (1 - btScalar(.0001) * scale);
						if (diagnose && sufficientDecrease) btDeformableDiagnostics::write("LINE_SEARCH_ROUNDOFF",
							"scale=%.9g energy_delta=%.17g allowance=%.9g residual_before=%.9g residual_after=%.9g",
							double(scale), energyDelta, double(roundoff), double(btSqrt(before)), double(btSqrt(after)));
					}
					if (captureCandidate) { trialScales.push_back(scale); trialEnergies.push_back(f1); }
				} while (!sufficientDecrease);
				if (lineSearchFailed)
				{
					diagnosticExit = "line_search_failed";
					// The trial evaluations modified m_dv, node velocities, temporary
					// positions, and deformation scratch data. Restore the accepted
					// iterate before terminating this Newton solve.
					revertDv();
					updateState();
					if (captureCandidate)
					{
						m_lastCapturedStep = btDeformableDiagnostics::current()->step;
						const char* logPath = std::getenv("BULLET_DEFORMABLE_DIAGNOSTICS");
						if (logPath && *logPath)
						{
							const std::string path = std::string(logPath) + ".step-" + std::to_string(m_lastCapturedStep) + ".newton.bin";
							btDeformableNewtonSnapshot snapshot;
							snapshot.dt=solverdt;snapshot.step=m_lastCapturedStep;snapshot.baselineEnergy=f0;
							snapshot.slope=inner_product;snapshot.initialScale=initialScale;
							snapshot.linearRecoveryBefore=linearRecoveryBefore;
							snapshot.scales=trialScales;snapshot.energies=trialEnergies;
							const bool saved=snapshot.save(path.c_str(),*this);
							btDeformableDiagnostics::write("NEWTON_CAPTURE", "success=%d path=%s trials=%d h=%.17g", int(saved), path.c_str(), trialScales.size(), double(solverdt));
						}
					}
					break;
				}
				revertDv();
				updateDv(scale);
			}
			else
			{
				computeStep(m_ddv, m_residual);
				const btScalar stepNorm = m_objective->computeNorm(m_ddv);
				const btScalar relativeStepTolerance = m_newtonTolerance * btMax(btScalar(1), m_objective->computeNorm(m_dv));
				const bool accurateLinearStep = !m_useProjection &&
					m_lastLinearMomentumResidual <= m_newtonTolerance &&
					m_lastLinearConstraintResidual <= m_newtonTolerance;
				const bool stepConverged = contactAccurate && accurateLinearStep &&
					m_lastStationarityResidual <= m_newtonTolerance && constraintError <= m_newtonTolerance;
				if (diagnose) btDeformableDiagnostics::write("NEWTON_LINEAR", "iteration=%d step_norm=%.9g step_tolerance=%.9g accurate=%d converged=%d momentum_residual=%.9g constraint_residual=%.9g stationarity=%.9g recovery=%d",
					i, double(stepNorm), double(relativeStepTolerance), int(accurateLinearStep), int(stepConverged), double(m_lastLinearMomentumResidual), double(m_lastLinearConstraintResidual), double(m_lastStationarityResidual), int(m_implicitRecoveryUsed));
				if (i > 0 && stepNorm <= relativeStepTolerance && (stepConverged || !accurateLinearStep))
				{
					m_lastSolveConverged = stepConverged;
					diagnosticExit = stepConverged ? "small_step_verified" : "small_step_unverified";
					// Check stationarity at the accepted state, not just the solved
					// linear system. Otherwise apply a small but accurate correction.
					break;
				}
				const btScalar scale = barrier ? barrier->safeStep(solverdt, m_ddv) : btScalar(1);
				if (barrier && scale < 1)
					btDeformableDiagnostics::write("VOLUME_STEP_LIMIT", "iteration=%d scale=%.12g", i, double(scale));
				if (scale < btScalar(1e-8)) { diagnosticExit = "volume_step_stalled"; break; }
				updateDv(scale);
			}
			for (int j = 0; j < m_numNodes; ++j)
			{
				m_ddv[j].setZero();
				m_residual[j].setZero();
			}
		}
		if (m_contactWeightedTarget > 0 && diagnosticExit == std::string("iteration_limit"))
		{
			// The last allowed update deserves the same residual check as every other iterate.
			updateState();
			bool admissible = !barrier || barrier->admissible();
			for (int n = 0; n < m_numNodes; ++n)
			{
				const auto* node = m_objective->m_nodes[n];
				m_residual[n] = node->m_im > 0 && node->m_frozen <= 0 ? -m_dv[n] / node->m_im : btVector3(0,0,0);
			}
			m_objective->computeResidual(solverdt, m_residual);
			m_objective->m_KKTPreconditioner->reinitialize(true);
			btScalar squared = 0;
			for (int n = 0; n < m_numNodes; ++n) squared += m_residual[n].dot(m_objective->m_KKTPreconditioner->applyInverseNodeBlock(n, m_residual[n]));
			if (admissible && std::isfinite(double(squared)) && squared >= 0 && btSqrt(squared) <= m_contactWeightedTarget &&
				m_objective->computeNorm(m_residual) < m_newtonTolerance && implicitConstraintError(m_dv) <= m_newtonTolerance)
			{ m_lastSolveConverged = true; diagnosticExit = "final_residual_tolerance"; }
		}
		updateVelocity();
		if (diagnose) btDeformableDiagnostics::write("NEWTON_EXIT", "reason=%s iterations=%d", diagnosticExit, diagnosticIterations);
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
		ImplicitOperatorCacheScope cacheScope(*m_objective);
		m_cg.solve(*m_objective, ddv, residual, false);
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

bool btDeformableBodySolver::setImplicitVelocityGuess(const btAlignedObjectArray<btVector3>& velocity)
{
	if (!m_implicit || m_useProjection || m_objective->m_projection.m_lagrangeMultipliers.size() || velocity.size() != m_numNodes) return false;
	TVStack guess = m_dv;
	int index = 0;
	for (int b = 0; b < m_softBodies.size(); ++b)
	{
		const btSoftBody& body = *m_softBodies[b];
		for (int n = 0; n < body.m_nodes.size(); ++n, ++index)
		{
			if (!body.isActive() || body.isStaticObject() || body.m_nodes[n].m_im <= 0 || body.m_nodes[n].m_frozen > 0) continue;
			guess[index] = velocity[index] - m_backupVelocity[index];
			for (int d = 0; d < 3; ++d) if (!std::isfinite(double(guess[index][d]))) return false;
		}
	}
	const TVStack original = m_dv;
	m_dv = guess; updateState();
	for (int f = 0; f < m_objective->m_lf.size(); ++f)
		if (m_objective->m_lf[f]->getForceType() == BT_VOLUME_BARRIER_FORCE &&
			!static_cast<btDeformableVolumeBarrierForce*>(m_objective->m_lf[f])->admissible())
		{
			m_dv = original; updateState(); return false;
		}
	// The beginning-of-step velocity and constraint targets remain unchanged.
	return true;
}

void btDeformableBodySolver::updateState()
{
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
	m_objective->m_preconditioner->reinitialize(true);
	const bool translation = m_objective->setupTranslationCorrection();
	m_objective->correctTranslation(x, rhs);
	// Keep the first Newton correction cheap. Later corrections retain Krylov
	// directions until a verified physical target or the per-call budget.
	const bool constrained = m_objective->m_projection.m_lagrangeMultipliers.size() > 0;
	const bool continued = m_newtonIteration > 0 && (constrained || translation || m_contactWeightedTarget > 0);
	const int linearBudget = continued ? (constrained ? 2400 : 1200) : kMaxConjugateGradientIterations;
	const btScalar physicalTarget = continued && m_contactWeightedTarget == 0 ? btScalar(0.5) * m_newtonTolerance : btScalar(0);
	const int iterations = m_cr.solveWithConvergencePolicy(*m_objective, x, rhs, false, !continued, true, linearBudget, physicalTarget, btScalar(.5) * m_contactWeightedTarget);
	btDeformableDiagnostics::write("KRYLOV_SOLVE", "newton=%d iterations=%d budget=%d initial=%.9g final=%.9g target=%.9g stagnated=%d",
		m_newtonIteration, iterations, linearBudget, double(m_cr.getInitialResidual()), double(m_cr.getFinalResidual()),
		double(m_cr.getTargetResidual()), int(m_cr.getStagnated()));
	m_objective->correctTranslation(x, rhs);
	measureImplicitLinearResidual(x, rhs);
	btScalar stepSquared = 0;
	for (int n = 0; n < m_numNodes; ++n) stepSquared += x[n].length2();
	const btScalar stepTolerance = m_newtonTolerance * btMax(btScalar(1), m_objective->computeNorm(m_dv));
	if (!continued && !m_implicitRecoveryUsed && btSqrt(stepSquared) <= stepTolerance &&
		(m_lastLinearMomentumResidual > m_newtonTolerance || m_lastLinearConstraintResidual > m_newtonTolerance))
	{
		// One full-budget restart per timestep for an inaccurate, tiny step.
		m_implicitRecoveryUsed = true;
		m_cr.solveWithConvergencePolicy(*m_objective, x, rhs, false, false, true, 0, 0, btScalar(.5) * m_contactWeightedTarget);
		m_objective->correctTranslation(x, rhs);
		measureImplicitLinearResidual(x, rhs);
	}
	m_objective->m_translationCorrection = false;
}

void btDeformableBodySolver::computeStep(TVStack& ddv, const TVStack& residual)
{
	if (m_useProjection)
	{
		ImplicitOperatorCacheScope cacheScope(*m_objective);
		m_cg.solve(*m_objective, ddv, residual, false);
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
		if (softBody->m_useSurfaceContact && btDeformableDiagnostics::enabled() &&
			(cp.m_contactPointFlags & BT_CONTACT_FLAG_PENETRATING) && btDeformableDiagnostics::current()->rigidSurfaceSamples++ < 8)
		{
			const auto point = resultOut->swapped ? cp.getPositionWorldOnB() - cp.m_normalWorldOnB * cp.getUnmodifiedDistance() : cp.getPositionWorldOnB();
			const auto normal = resultOut->swapped ? -cp.m_normalWorldOnB : cp.m_normalWorldOnB;
			btAlignedObjectArray<btSoftBody::ContactNode> stencil;
			const bool mapped = softBody->appendSurfaceContactNodes(resultOut->getPartId0(), resultOut->getIndex0(), point, 1, stencil);
			btVector3 surfaceVelocity(0,0,0), surfaceSplit(0,0,0);
			for (int n = 0; n < stencil.size(); ++n)
			{
				surfaceVelocity += stencil[n].jacobian * stencil[n].node->m_v;
				surfaceSplit += stencil[n].jacobian * stencil[n].node->m_splitv;
			}
			btDeformableDiagnostics::write("RIGID_SURFACE_SAMPLE", "body=%d other=%d part=%d triangle=%d swapped=%d mapped=%d nodes=%d point_source=legacy_recovery modified_distance=%.12g recovery_distance=%.12g normal_velocity=%.12g split_normal_velocity=%.12g point=%.12g,%.12g,%.12g normal=%.12g,%.12g,%.12g",
				softBody->getUserIndex(), collisionObjectWrap->getCollisionObject()->getUserIndex(), resultOut->getPartId0(), resultOut->getIndex0(), int(resultOut->swapped), int(mapped), stencil.size(),
				double(cp.getDistance()), double(cp.getUnmodifiedDistance()), double(normal.dot(surfaceVelocity)), double(normal.dot(surfaceSplit)),
				double(point.x()), double(point.y()), double(point.z()), double(normal.x()), double(normal.y()), double(normal.z()));
		}
		if (softBody->m_useSurfaceContact)
		{
			const auto point = resultOut->swapped ? cp.getPositionWorldOnB() - cp.m_normalWorldOnB * cp.getUnmodifiedDistance() : cp.getPositionWorldOnB();
			const auto normal = resultOut->swapped ? -cp.m_normalWorldOnB : cp.m_normalWorldOnB;
			softBody->skinSoftStaticCollisionHandler(collisionObjectWrap, resultOut->getPartId0(), resultOut->getIndex0(), resultOut->getPartId1(), resultOut->getIndex1(),
				point, normal, cp.getUnmodifiedDistance(), (cp.m_contactPointFlags & BT_CONTACT_FLAG_PENETRATING) != 0, &cp.m_appliedImpulse);
			return;
		}
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
										   &cp.m_appliedImpulse, cp.getUnmodifiedDistance());
}
