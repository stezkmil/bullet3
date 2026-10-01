#ifndef BT_DEFORMABLE_ENERGY_CHANGE_H
#define BT_DEFORMABLE_ENERGY_CHANGE_H
#include "btDeformableBackwardEulerObjective.h"

// Subtract nonlinear force energies separately; linear potentials and inertia
// use algebraic differences so large absolute offsets do not enter the sum.
class btDeformableEnergyChange
{
	btDeformableBackwardEulerObjective& objective;
	btAlignedObjectArray<btVector3> baseline;
	btAlignedObjectArray<double> energies;
	btScalar dt;
public:
	btDeformableEnergyChange(btDeformableBackwardEulerObjective& o, const btAlignedObjectArray<btVector3>& dv, btScalar step)
		: objective(o), baseline(dv), dt(step)
	{
		for (int f = 0; f < o.m_lf.size(); ++f) energies.push_back(o.m_lf[f]->totalEnergy(dt));
	}
	double difference(const btAlignedObjectArray<btVector3>& dv)
	{
		btAlignedObjectArray<btVector3> delta;delta.resize(dv.size());
		double change = 0;
		for (int n = 0; n < dv.size(); ++n)
		{
			delta[n] = dv[n] - baseline[n];
			const auto* node = objective.m_nodes[n];
			if (node->m_im > 0 && node->m_frozen <= 0)
				change += double(delta[n].dot(dv[n] + baseline[n]) * btScalar(.5) / node->m_im);
		}
		for (int f = 0; f < objective.m_lf.size(); ++f)
		{
			double term;
			if (!objective.m_lf[f]->linearEnergyChange(dt, delta, term)) term = objective.m_lf[f]->totalEnergy(dt) - energies[f];
			change += term;
		}
		return change;
	}
};
#endif
