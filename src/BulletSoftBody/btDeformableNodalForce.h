/*
 Bullet Continuous Collision Detection and Physics Library
 Copyright (c) 2019 Google Inc. http://bulletphysics.org
 This source is provided under the same terms as the rest of BulletSoftBody.
 Altered source: constant force support for selected deformable nodes.
*/

#ifndef BT_DEFORMABLE_NODAL_FORCE_H
#define BT_DEFORMABLE_NODAL_FORCE_H

#include "btDeformableLagrangianForce.h"

class btDeformableNodalForce : public btDeformableLagrangianForce
{
	btAlignedObjectArray<int> m_nodeIndices;
	btVector3 m_totalForce;

	void activateSoftBodies()
	{
		for (int bodyIndex = 0; bodyIndex < m_softBodies.size(); ++bodyIndex)
		{
			if (m_softBodies[bodyIndex])
				m_softBodies[bodyIndex]->activate();
		}
	}

	bool canApplyForce(const btSoftBody::Node& node) const
	{
		return node.m_frozen <= 0 && node.m_im > 0;
	}

	int getAffectedNodeCount(const btSoftBody& softBody) const
	{
		int count = 0;
		for (int i = 0; i < m_nodeIndices.size(); ++i)
		{
			const int nodeIndex = m_nodeIndices[i];
			if (nodeIndex >= 0 && nodeIndex < softBody.m_nodes.size() && canApplyForce(softBody.m_nodes[nodeIndex]))
				++count;
		}
		return count;
	}

	void addScaledNodalForce(btScalar scale, TVStack& force)
	{
		for (int bodyIndex = 0; bodyIndex < m_softBodies.size(); ++bodyIndex)
		{
			btSoftBody* softBody = m_softBodies[bodyIndex];
			if (!softBody || !softBody->isActive() || softBody->isStaticObject())
				continue;

			const int affectedNodeCount = getAffectedNodeCount(*softBody);
			if (affectedNodeCount == 0)
				continue;
			const btVector3 forcePerNode = m_totalForce / btScalar(affectedNodeCount);

			for (int i = 0; i < m_nodeIndices.size(); ++i)
			{
				const int nodeIndex = m_nodeIndices[i];
				if (nodeIndex >= 0 && nodeIndex < softBody->m_nodes.size())
				{
					const btSoftBody::Node& node = softBody->m_nodes[nodeIndex];
					if (canApplyForce(node))
						force[node.index] += scale * forcePerNode;
				}
			}
		}
	}

public:
	const btAlignedObjectArray<int>& nodeIndices() const { return m_nodeIndices; }
	const btVector3& totalForce() const { return m_totalForce; }
	btDeformableNodalForce(btSoftBody* softBody, const btAlignedObjectArray<int>& nodeIndices, const btVector3& totalForce)
		: m_nodeIndices(nodeIndices), m_totalForce(totalForce)
	{
		if (softBody)
			addSoftBody(softBody);
		activateSoftBodies();
	}

	void setForce(const btVector3& totalForce)
	{
		m_totalForce = totalForce;
		activateSoftBodies();
	}
	void setNodeIndices(const btAlignedObjectArray<int>& nodeIndices) { m_nodeIndices = nodeIndices; }

	void addScaledForces(btScalar scale, TVStack& force) override { addScaledNodalForce(scale, force); }
	void addScaledExplicitForce(btScalar scale, TVStack& force) override { addScaledNodalForce(scale, force); }
	void addScaledDampingForce(btScalar, TVStack&) override {}
	void addScaledDampingForceDifferential(btScalar, const TVStack&, TVStack&) override {}
	void buildDampingForceDifferentialDiagonal(btScalar, TVStack&) override {}
	void addScaledElasticForceDifferential(btScalar, const TVStack&, TVStack&) override {}

	bool linearEnergyChange(btScalar dt, const TVStack& deltaVelocity, double& change) override
	{
		change = 0;
		for (int b = 0; b < m_softBodies.size(); ++b)
		{
			const auto* body = m_softBodies[b];
			if (!body || !body->isActive() || body->isStaticObject()) continue;
			const int count = getAffectedNodeCount(*body); if (!count) continue;
			const btVector3 force = m_totalForce / btScalar(count);
			for (int i = 0; i < m_nodeIndices.size(); ++i)
			{
				const int n = m_nodeIndices[i];
				if (n >= 0 && n < body->m_nodes.size() && canApplyForce(body->m_nodes[n]))
					change -= double(dt * force.dot(deltaVelocity[body->m_nodes[n].index]));
			}
		}
		return true;
	}

	double totalElasticEnergy(btScalar) override
	{
		double energy = 0;
		for (int bodyIndex = 0; bodyIndex < m_softBodies.size(); ++bodyIndex)
		{
			const btSoftBody* softBody = m_softBodies[bodyIndex];
			if (!softBody || !softBody->isActive() || softBody->isStaticObject())
				continue;
			const int affectedNodeCount = getAffectedNodeCount(*softBody);
			if (affectedNodeCount == 0)
				continue;
			const btVector3 forcePerNode = m_totalForce / btScalar(affectedNodeCount);
			for (int i = 0; i < m_nodeIndices.size(); ++i)
			{
				const int nodeIndex = m_nodeIndices[i];
				if (nodeIndex >= 0 && nodeIndex < softBody->m_nodes.size() && canApplyForce(softBody->m_nodes[nodeIndex]))
					energy -= forcePerNode.dot(softBody->m_nodes[nodeIndex].m_q);
			}
		}
		return energy;
	}

	btDeformableLagrangianForceType getForceType() override { return BT_NODAL_FORCE; }
};

#endif /* BT_DEFORMABLE_NODAL_FORCE_H */
