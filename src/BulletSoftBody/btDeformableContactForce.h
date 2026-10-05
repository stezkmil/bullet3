#ifndef BT_DEFORMABLE_CONTACT_FORCE_H
#define BT_DEFORMABLE_CONTACT_FORCE_H

#include "btDeformableLagrangianForce.h"
#include "btDeformableContactProjection.h"
#include "btPreconditioner.h"
#include <cmath>

// Extend useful multiplier iterations before resorting to a smaller physical timestep.
class btDeformableContactConvergence
{
	btScalar m_checkpointError = SIMD_INFINITY;
	int m_limit = 20;
	btScalar m_normalCheckpoint = SIMD_INFINITY;
	int m_penaltyScale = 1;

public:
	int limit() const
	{
		return m_limit;
	}
	bool increaseNormalPenalty(int completed, btScalar error, bool newtonConverged)
	{
		if (completed == 1)
			m_normalCheckpoint = error;
		if (completed % 5 != 0)
			return false;
		const bool increase = newtonConverged && std::isfinite(double(error)) && std::isfinite(double(m_normalCheckpoint)) &&
							  error > btScalar(.001) && error > m_normalCheckpoint * btScalar(.5) && m_penaltyScale < 64;
		m_normalCheckpoint = error;
		if (increase)
			m_penaltyScale *= 4;
		return increase;
	}
	int penaltyScale() const
	{
		return m_penaltyScale;
	}
	// Preserve the original impulse-to-velocity error scale after penalty growth.
	btScalar normalError(btScalar error) const
	{
		return error * m_penaltyScale;
	}
	bool observe(int completed, btScalar error, bool newtonConverged, bool penaltyIncreased = false)
	{
		if (completed == 1)
			m_checkpointError = error;
		if (!std::isfinite(double(error)) || !std::isfinite(double(m_checkpointError)) || completed != m_limit || m_limit >= 100 ||
			!newtonConverged || !(error >= 0 && (error < m_checkpointError * btScalar(.5) || penaltyIncreased)))
			return false;
		m_checkpointError = error;
		m_limit += 20;
		return true;
	}
};

// Augmented contact impulses, solved together with elasticity. Geometry and
// the friction radius stay fixed during each Newton solve.
class btDeformableContactForce : public btDeformableLagrangianForce
{
public:
	struct Contact
	{
		btAlignedObjectArray<btSoftBody::ContactNode> nodes;
		btVector3 normal, tangentImpulse;
		btScalar gap, friction, rho, tangentRho, normalImpulse;
	};
	btAlignedObjectArray<Contact> contacts;
	btScalar dt;
	btScalar lastNormalError = 0, lastTangentError = 0;
	int worstNormalContact = -1, worstTangentContact = -1;

	explicit btDeformableContactForce(btScalar step) : dt(step) {}
	static bool movable(const btSoftBody::Node* n) { return n->m_im > 0 && n->m_frozen <= 0; }
	bool add(const btSoftBody::DeformableNodeNodeContact& c)
	{
		if (c.m_surfaceInvalid) return false;
		if (c.m_normal.length2() < SIMD_EPSILON) return false;
		Contact v;
		v.nodes = c.m_surfaceNodes;
		if (!v.nodes.size())
		{
			if (!c.m_node0 || !c.m_node1) return false;
			if (c.m_node0 == c.m_node1) return true;
			btSoftBody::ContactNode a = {c.m_node0, btMatrix3x3::getIdentity()};
			btSoftBody::ContactNode b = {c.m_node1, btMatrix3x3::getIdentity() * btScalar(-1)};
			v.nodes.push_back(a); v.nodes.push_back(b);
		}
		v.normal = c.m_normal.normalized();
		v.gap = c.m_offset; v.friction = btMax(btScalar(0), c.m_friction);
		btScalar inverseMass = 0;
		for (int n = 0; n < v.nodes.size(); ++n)
		{
			const auto& entry = v.nodes[n];
			if (movable(entry.node)) inverseMass += entry.node->m_im * (entry.jacobian.transpose() * v.normal).length2();
		}
		// Keep immovable contacts in validation; do not turn them into free motion.
		v.rho = v.tangentRho = inverseMass > 0 ? btScalar(10) / inverseMass : btScalar(1);
		v.normalImpulse = 0; v.tangentImpulse.setZero();
		for (int i = 0; i < contacts.size(); ++i)
		{
			Contact& old = contacts[i];
			const btScalar orientation = old.normal.dot(v.normal) >= 0 ? 1 : -1;
			bool same = old.nodes.size() == v.nodes.size() && old.normal.dot(v.normal) * orientation > btScalar(0.999999);
			for (int n = 0; same && n < v.nodes.size(); ++n)
			{
				bool found = false;
				for (int k = 0; k < old.nodes.size(); ++k)
					if (old.nodes[k].node == v.nodes[n].node)
					{
						const btMatrix3x3 diff = old.nodes[k].jacobian - v.nodes[n].jacobian * orientation;
						found = diff[0].length2() + diff[1].length2() + diff[2].length2() < btScalar(1e-18);
						break;
					}
				same = found;
			}
			if (same) { old.gap = btMin(old.gap, v.gap); return true; }
		}
		contacts.push_back(v);
		return true;
	}
	btVector3 velocity(const Contact& c, bool split = false) const
	{
		btVector3 result(0,0,0);
		for (int n = 0; n < c.nodes.size(); ++n)
			result += c.nodes[n].jacobian * (split ? c.nodes[n].node->m_splitv : c.nodes[n].node->m_v);
		return result;
	}

	// Approximate constrained compliance D-D*C^T*(C*D*C^T)^-1*C*D.
	// D includes elastic/damping blocks; row projections also handle face constraints.
	static btScalar compliance(const Contact& c, const btVector3& direction, const KKTPreconditioner& preconditioner,
		const btAlignedObjectArray<LagrangeMultiplier>& constraints, int nodeCount)
	{
		TVStack response;
		response.resize(nodeCount, btVector3(0,0,0));
		for (int n = 0; n < c.nodes.size(); ++n)
		{
			const auto& entry = c.nodes[n];
			if (movable(entry.node)) response[entry.node->index] += preconditioner.applyInverseNodeBlock(entry.node->index, entry.jacobian.transpose() * direction);
		}
		for (int sweep = 0; sweep < 64; ++sweep)
		{
			btScalar largest = 0;
			for (int k = 0; k < constraints.size(); ++k)
			{
				const LagrangeMultiplier& row = constraints[k];
				for (int d = 0; d < row.m_num_constraints; ++d)
				{
					btScalar residual = 0, diagonal = 0;
					for (int n = 0; n < row.m_num_nodes; ++n)
					{
						const btVector3 gradient = row.m_weights[n] * row.m_dirs[d];
						residual += gradient.dot(response[row.m_indices[n]]);
						diagonal += gradient.dot(preconditioner.applyInverseNodeBlock(row.m_indices[n], gradient));
					}
					if (diagonal <= 0) continue;
					largest = btMax(largest, btFabs(residual));
					for (int n = 0; n < row.m_num_nodes; ++n)
						response[row.m_indices[n]] -= preconditioner.applyInverseNodeBlock(row.m_indices[n], row.m_weights[n] * row.m_dirs[d]) * (residual / diagonal);
				}
			}
			if (largest <= btScalar(1e-10)) break;
		}
		btScalar result = 0;
		for (int n = 0; n < c.nodes.size(); ++n) result += direction.dot(c.nodes[n].jacobian * response[c.nodes[n].node->index]);
		return btMax(btScalar(0), result);
	}
	void configureScaling(const KKTPreconditioner& preconditioner, const btAlignedObjectArray<LagrangeMultiplier>& constraints, int nodeCount)
	{
		for (int i = 0; i < contacts.size(); ++i)
		{
			Contact& c = contacts[i];
			btVector3 t1, t2; btPlaneSpace1(c.normal, t1, t2);
			const btScalar normalResponse = compliance(c, c.normal, preconditioner, constraints, nodeCount);
			const btScalar tangentResponse = (compliance(c, t1, preconditioner, constraints, nodeCount) + compliance(c, t2, preconditioner, constraints, nodeCount)) * btScalar(0.5);
			const btScalar minimumResponse = btScalar(10) / (c.rho * btScalar(1e6));
			c.rho = btScalar(10) / btMax(normalResponse, minimumResponse);
			c.tangentRho = btScalar(10) / btMax(tangentResponse, minimumResponse);
		}
	}
	btScalar normalImpulse(const Contact& c) const
	{
		return btMax(btScalar(0), c.normalImpulse - c.rho * (c.gap / dt + c.normal.dot(velocity(c) + velocity(c, true))));
	}
	btVector3 tangentTrial(const Contact& c) const
	{
		const btVector3 v = velocity(c);
		return c.tangentImpulse - c.tangentRho * (v - c.normal * c.normal.dot(v));
	}
	btVector3 tangentImpulse(const Contact& c) const
	{
		const btVector3 trial = tangentTrial(c);
		const btScalar length = trial.length(), radius = c.friction * c.normalImpulse;
		return length > radius && length > 0 ? trial * (radius / length) : trial;
	}
	// Returns a velocity-equivalent complementarity / friction residual.
	btScalar updateMultipliers()
	{
		btScalar error = 0;
		lastNormalError = lastTangentError = 0;
		worstNormalContact = worstTangentContact = -1;
		for (int i = 0; i < contacts.size(); ++i)
		{
			Contact& c = contacts[i];
			const btScalar j = normalImpulse(c);
			const btVector3 t = tangentImpulse(c);
			const btScalar normalError = btFabs(j - c.normalImpulse) / c.rho;
			const btScalar tangentError = (t - c.tangentImpulse).length() / c.tangentRho;
			if (normalError > lastNormalError) { lastNormalError = normalError; worstNormalContact = i; }
			if (tangentError > lastTangentError) { lastTangentError = tangentError; worstTangentContact = i; }
			error = btMax(error, btMax(normalError, tangentError));
			c.normalImpulse = j; c.tangentImpulse = t;
		}
		return error;
	}
	virtual void addScaledForces(btScalar scale, TVStack& f)
	{
		for (int i = 0; i < contacts.size(); ++i)
		{
			const Contact& c = contacts[i];
			const btVector3 impulse = (scale / dt) * (normalImpulse(c) * c.normal + tangentImpulse(c));
			for (int n = 0; n < c.nodes.size(); ++n)
				if (movable(c.nodes[n].node)) f[c.nodes[n].node->index] += c.nodes[n].jacobian.transpose() * impulse;
		}
	}
	btMatrix3x3 hessian(const Contact& c) const
	{
		btMatrix3x3 nn, tt;
		for (int r = 0; r < 3; ++r)
			for (int d = 0; d < 3; ++d) nn[r][d] = c.normal[r] * c.normal[d];
		tt = btMatrix3x3::getIdentity() - nn;
		const btVector3 trial = tangentTrial(c);
		const btScalar length = trial.length(), radius = c.friction * c.normalImpulse;
		if (radius <= 0) tt = tt * btScalar(0);
		else if (length > radius)
		{
			const btVector3 direction = trial / length;
			for (int r = 0; r < 3; ++r)
				for (int d = 0; d < 3; ++d) tt[r][d] -= direction[r] * direction[d];
			tt = tt * (radius / length);
		}
		if (normalImpulse(c) <= 0) nn = nn * btScalar(0);
		return nn * c.rho + tt * c.tangentRho;
	}
	virtual void addImplicitForceDifferential(btScalar, const TVStack& x, TVStack& out)
	{
		for (int i = 0; i < contacts.size(); ++i)
		{
			const Contact& c = contacts[i];
			btVector3 relative(0,0,0);
			for (int n = 0; n < c.nodes.size(); ++n)
				if (movable(c.nodes[n].node)) relative += c.nodes[n].jacobian * x[c.nodes[n].node->index];
			const btVector3 d = hessian(c) * relative;
			for (int n = 0; n < c.nodes.size(); ++n)
				if (movable(c.nodes[n].node)) out[c.nodes[n].node->index] += c.nodes[n].jacobian.transpose() * d;
		}
	}
	virtual bool addImplicitForceDifferentialBlocks(btScalar, btAlignedObjectArray<btMatrix3x3>& blocks)
	{
		for (int i = 0; i < contacts.size(); ++i)
		{
			const Contact& c = contacts[i]; const btMatrix3x3 h = hessian(c);
			for (int n = 0; n < c.nodes.size(); ++n)
				if (movable(c.nodes[n].node)) blocks[c.nodes[n].node->index] += c.nodes[n].jacobian.transpose() * h * c.nodes[n].jacobian;
		}
		return true;
	}
	virtual double totalEnergy(btScalar)
	{
		double energy = 0;
		for (int i = 0; i < contacts.size(); ++i)
		{
			const Contact& c = contacts[i];
			const btScalar j = normalImpulse(c), t = tangentTrial(c).length(), radius = c.friction * c.normalImpulse;
			energy += double(j) * j / (2 * c.rho);
			energy += (t <= radius ? double(t)*t/2 : double(radius)*(t-radius/2)) / c.tangentRho;
		}
		return energy;
	}
	virtual void addScaledDampingForceDifferential(btScalar, const TVStack&, TVStack&) {}
	virtual void addScaledElasticForceDifferential(btScalar, const TVStack&, TVStack&) {}
	virtual void buildDampingForceDifferentialDiagonal(btScalar, TVStack&) {}
	virtual void addScaledExplicitForce(btScalar, TVStack&) {}
	virtual void addScaledDampingForce(btScalar, TVStack&) {}
	virtual btDeformableLagrangianForceType getForceType() { return BT_CONTACT_FORCE; }
};
#endif
