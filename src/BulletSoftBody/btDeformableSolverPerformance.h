// Temporary baseline instrumentation. Timings are inclusive and overlap.
#ifndef BT_DEFORMABLE_SOLVER_PERFORMANCE_H
#define BT_DEFORMABLE_SOLVER_PERFORMANCE_H
#include <chrono>
#include "LinearMath/btAlignedObjectArray.h"
struct btDeformableSolverPerformance
{
    struct Sample { double ms = 0; int calls = 0; };
    bool active = false;
    int newton = 0, krylovIterations = 0, lineTrials = 0, lineFailures = 0;
    struct LinearRecord
    {
        int newton, iterations, budget;
        bool recovery, translation, relative;
        const char* stop;
        double ms, initial, final, target, physicalTarget;
        double momentum, constraint, stationarity;
        double rhsMomentum, rhsConstraint, stepL2, checkpointMs;
        int progressBegin, progressEnd;
    };
    struct ProgressRecord
    {
        int iteration;
        double recurrenceWeighted, recurrenceL2, trueMomentumL2, trueConstraintInf, residualGapL2, weightedBest;
    };
    struct NewtonRecord
    {
        int iteration;
        double inputForceL2 = 0, inputConstraintInf = 0, stepL2 = 0, stepTolerance = 0;
        double stationarityL2 = 0, linearMomentumL2 = 0, linearConstraintInf = 0;
        double appliedScale = 0, outputConstraintInf = 0;
        bool solved = false;
        const char* outcome = "pending";
    };
    btAlignedObjectArray<NewtonRecord> newtonRecords;
    btAlignedObjectArray<ProgressRecord> progressRecords;
    btAlignedObjectArray<LinearRecord> linearRecords;
    Sample cacheSetup;
    Sample combined, multiply, damping, elastic, precondition, blockSetup, translationSetup;
    Sample translationCorrect, linear, verification, state, energy, residual;
};
class btDeformablePerformanceScope
{
    typedef std::chrono::steady_clock Clock;
    btDeformableSolverPerformance::Sample& m_sample;
    bool m_active;
    Clock::time_point m_start;
public:
    btDeformablePerformanceScope(btDeformableSolverPerformance::Sample& sample, bool active)
        : m_sample(sample), m_active(active)
    {
        if (m_active) m_start = Clock::now();
    }
    ~btDeformablePerformanceScope()
    {
        if (m_active)
        {
            m_sample.ms += std::chrono::duration<double, std::milli>(Clock::now() - m_start).count();
            ++m_sample.calls;
        }
    }
};
#endif
