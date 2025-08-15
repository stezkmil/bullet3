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

#include "btDeformableMousePickingForce.h"

const btSoftBody::Node* btDeformableMousePickingForce::getNode(int i) const
{
	if (m_face) return m_face->m_n[i];
	if (m_tetra) return m_tetra->m_n[i];
	return m_nodes[i];
}

int btDeformableMousePickingForce::getIndexCount() const
{
	if (m_face) return 3;
	if (m_tetra) return 4;
	return static_cast<int>(m_nodes.size());
}

void btDeformableMousePickingForce::calculateNodeToMouse()
{
	for (int i = 0; i < getIndexCount(); ++i)
	{
		m_node_to_mouse_x[i] = getNode(i)->m_x - m_mouse_transform.getOrigin();
		m_node_to_mouse_x_orig[i] = m_node_to_mouse_x[i];
		m_node_to_mouse_q[i] = getNode(i)->m_q - m_mouse_transform.getOrigin();
		m_node_to_mouse_q_orig[i] = m_node_to_mouse_q[i];
	}
	m_mouse_bary = BaryCoord(getNode(0)->m_x, getNode(1)->m_x, getNode(2)->m_x, getNode(3)->m_x, m_mouse_transform.getOrigin());
}

void btDeformableMousePickingForce::rotateNodeToMouse()
{
	btTransform delta = m_mouse_transform * m_mouse_transform_orig_inv;
	for (int i = 0; i < getIndexCount(); ++i)
	{
		m_node_to_mouse_x[i] = delta.getBasis() * m_node_to_mouse_x_orig[i];
		m_node_to_mouse_q[i] = delta.getBasis() * m_node_to_mouse_q_orig[i];
	}
}

void btDeformableMousePickingForce::orientationError(btScalar& theta, btVector3& axis) const
{
	int n = getIndexCount();
	btVector3 c_now(0, 0, 0), c_rest(0, 0, 0);
	for (int i = 0; i < n; ++i)
	{
		c_now += getNode(i)->m_x;
		c_rest += m_node_to_mouse_x[i];
	}
	c_now *= btScalar(1.0) / n;
	c_rest *= btScalar(1.0) / n;

	btVector3 p[4], q[4];
	for (int i = 0; i < n; ++i)
	{
		p[i] = getNode(i)->m_x - c_now;
		q[i] = m_node_to_mouse_x[i] - c_rest;
	}

	btMatrix3x3 H(btScalar(0), btScalar(0), btScalar(0),
				  btScalar(0), btScalar(0), btScalar(0),
				  btScalar(0), btScalar(0), btScalar(0));
	for (int i = 0; i < n; ++i)
	{
		H[0][0] += p[i].x() * q[i].x();
		H[0][1] += p[i].x() * q[i].y();
		H[0][2] += p[i].x() * q[i].z();
		H[1][0] += p[i].y() * q[i].x();
		H[1][1] += p[i].y() * q[i].y();
		H[1][2] += p[i].y() * q[i].z();
		H[2][0] += p[i].z() * q[i].x();
		H[2][1] += p[i].z() * q[i].y();
		H[2][2] += p[i].z() * q[i].z();
	}
	btMatrix3x3 R, S;
	static const btPolarDecomposition polar;
	polar.decompose(H, R, S);
	{
		btQuaternion q;
		R.getRotation(q);
		axis = q.getAxis();
		theta = q.getAngle();
	}
}

btScalar btDeformableMousePickingForce::transBoost(btScalar theta) const
{
	if (m_angleMax <= btScalar(0)) return btScalar(1);
	btScalar b = btScalar(1) + theta / m_angleMax;
	if (b > m_scaleCeil) b = m_scaleCeil;
	return b;
}

// ---------------------------------------------------------------------
void btDeformableMousePickingForce::applyTorqueSpring(const btVector3& T, TVStack& force)
{
	int n = getIndexCount();
	if (n == 0) return;
	//btVector3 picked = /*m_mouse_transform*/ m_mouse_transform_orig.getOrigin();
	btScalar denom = btScalar(0);
	btVector3 r[4];
	btVector3 pickedFromBary = BaryEval(getNode(0)->m_x, getNode(1)->m_x, getNode(2)->m_x, getNode(3)->m_x, m_mouse_bary);

	const auto a = getNode(0);
	const auto b = getNode(1);
	const auto c = getNode(2);
	const auto d = getNode(3);
	const btSoftBody::Node* nodeArr[4] = {a, b, c, d};
	fprintf(stderr, "drawline \"line\" [%f,%f,%f][%f,%f,%f][0,1,0,1] \n", a->m_x.x(), a->m_x.y(), a->m_x.z(),
			b->m_x.x(), b->m_x.y(), b->m_x.z());
	fprintf(stderr, "drawline \"line\" [%f,%f,%f][%f,%f,%f][0,1,0,1] \n", a->m_x.x(), a->m_x.y(), a->m_x.z(),
			c->m_x.x(), c->m_x.y(), c->m_x.z());
	fprintf(stderr, "drawline \"line\" [%f,%f,%f][%f,%f,%f][0,1,0,1] \n", a->m_x.x(), a->m_x.y(), a->m_x.z(),
			d->m_x.x(), d->m_x.y(), d->m_x.z());
	fprintf(stderr, "drawline \"line\" [%f,%f,%f][%f,%f,%f][0,1,0,1] \n", b->m_x.x(), b->m_x.y(), b->m_x.z(),
			d->m_x.x(), d->m_x.y(), d->m_x.z());
	fprintf(stderr, "drawline \"line\" [%f,%f,%f][%f,%f,%f][0,1,0,1] \n", b->m_x.x(), b->m_x.y(), b->m_x.z(),
			c->m_x.x(), c->m_x.y(), c->m_x.z());
	fprintf(stderr, "drawline \"line\" [%f,%f,%f][%f,%f,%f][0,1,0,1] \n", d->m_x.x(), d->m_x.y(), d->m_x.z(),
			c->m_x.x(), c->m_x.y(), c->m_x.z());

	fprintf(stderr, "drawline \"axisX\" [%f,%f,%f][%f,%f,%f][1,0,0,1] \n", a->m_x.x(), a->m_x.y(), a->m_x.z(),
			a->m_x.x() + 10.0, a->m_x.y(), a->m_x.z());
	fprintf(stderr, "drawline \"axisY\" [%f,%f,%f][%f,%f,%f][0,1,0,1] \n", a->m_x.x(), a->m_x.y(), a->m_x.z(),
			a->m_x.x(), a->m_x.y() + 10.0, a->m_x.z());
	fprintf(stderr, "drawline \"axisZ\" [%f,%f,%f][%f,%f,%f][0,0,1,1] \n", a->m_x.x(), a->m_x.y(), a->m_x.z(),
			a->m_x.x(), a->m_x.y(), a->m_x.z() + 10.0);

	fprintf(stderr, "drawpoint \"picked\" [%f,%f,%f][1,0,0,1] \n", pickedFromBary.x(), pickedFromBary.y(), pickedFromBary.z());
	fprintf(stderr, "drawline \"T\" [%f,%f,%f][%f,%f,%f][1,0,0,1] \n", pickedFromBary.x(), pickedFromBary.y(), pickedFromBary.z(),
			pickedFromBary.x() + (T.x() * 10000.0), pickedFromBary.y() + (T.y() * 10000.0), pickedFromBary.z() + (T.z() * 10000.0));

	for (int i = 0; i < n; ++i)
	{
		r[i] = getNode(i)->m_x - pickedFromBary;
		denom += r[i].length2();
	}
	if (denom < btScalar(1e-7)) return;

	std::string str;
	for (int i = 0; i < n; ++i)
	{
		btVector3 dbg(0, 0, 1);
		if (T./*dbg.*/ cross(r[i].normalized()).norm() > SIMD_EPSILON)
		{
			btVector3 F = (T./*dbg.*/ cross(r[i].normalized()) /** 30.0*/ /*/ denom*/) * 100.0;
			str += " F " + std::to_string(F.x()) + " " + std::to_string(F.y()) + " " + std::to_string(F.z());
			/*if (F.safeNorm() > m_maxForce)
			{
				F.safeNormalize();
				F *= m_maxForce;
			}*/
			force[getNode(i)->index] -= F;
			fprintf(stderr, "drawline \"line\" [%f,%f,%f][%f,%f,%f][1,0,0,1] \n", nodeArr[i]->m_x.x(), nodeArr[i]->m_x.y(), nodeArr[i]->m_x.z(),
					nodeArr[i]->m_x.x() + F.x(), nodeArr[i]->m_x.y() + F.y(), nodeArr[i]->m_x.z() + F.z());
		}
	}
	//fprintf(stderr, "%s\n", str.c_str());
}

btDeformableMousePickingForce::btDeformableMousePickingForce(
	btScalar elasticK,
	btScalar dampingK,
	const btSoftBody::Face* face,
	const btSoftBody::Tetra* tetra,
	std::vector<const btSoftBody::Node*>* nodes,
	const btTransform& mouseTr,
	btScalar maxForce,
	btScalar angleMax,   // off by default
	btScalar scaleCeil)  // no boost by default
	: m_elasticStiffness(elasticK),
	  m_dampingStiffness(dampingK),
	  m_angleMax(angleMax),
	  m_scaleCeil(scaleCeil),
	  m_face(face),
	  m_tetra(tetra),
	  m_mouse_transform(mouseTr),
	  m_mouse_transform_orig(mouseTr),
	  m_maxForce(maxForce)
{
	if (nodes) m_nodes = *nodes;
	m_mouse_transform_orig_inv = m_mouse_transform_orig.inverse();
	calculateNodeToMouse();
}

// ------------------------------------------------- Force interface
void btDeformableMousePickingForce::addScaledForces(btScalar s, TVStack& f)
{
	addScaledDampingForce(s, f);
	addScaledElasticForce(s, f);
}
void btDeformableMousePickingForce::addScaledExplicitForce(btScalar s, TVStack& f) { addScaledElasticForce(s, f); }

// ---------------- elastic ------------------------------------------------
void btDeformableMousePickingForce::addScaledElasticForce(btScalar scale, TVStack& force)
{
	//fprintf(stderr, "framestart()\n");
	//btScalar theta;
	//btVector3 axis;
	//orientationError(theta, axis);
	//btScalar orientationBasedStiffnessBoost = transBoost(theta);
	btScalar k = scale * m_elasticStiffness /* * orientationBasedStiffnessBoost*/;
	for (int i = 0; i < getIndexCount(); ++i)
	{
		const btSoftBody::Node* n = getNode(i);
		btVector3 diff = (n->m_x - m_mouse_transform.getOrigin()) - m_node_to_mouse_x[i];
		btVector3 f = k * diff;
		if (f.safeNorm() > m_maxForce)
		{
			f.safeNormalize();
			f *= m_maxForce;
		}
		force[n->index] -= f;
	}
	// Should add a rotation to preserve the grab orientation, but I did not observe any positive effects even with high values of m_kRot
	//if (m_kRot > btScalar(0) /*&& theta > btScalar(1e-4)*/)
	//{
	//	btVector3 torque = m_kRot * theta * axis;
	//	applyTorqueSpring(torque, force);
	//}
	//fprintf(stderr, "frameend()\n");
}

// --------------- damping -------------------------------------------------
void btDeformableMousePickingForce::addScaledDampingForce(btScalar scale, TVStack& force)
{
	//btScalar theta;
	//btVector3 axis;
	//orientationError(theta, axis);
	//btScalar orientationBasedStiffnessBoost = transBoost(theta);
	btScalar k = scale * m_dampingStiffness /** orientationBasedStiffnessBoost*/;
	for (int i = 0; i < getIndexCount(); ++i)
	{
		const btSoftBody::Node* n = getNode(i);
		btVector3 diff = (n->m_x - m_mouse_transform.getOrigin()) - m_node_to_mouse_x[i];
		btVector3 v = n->m_v;
		btVector3 f = k * v;
		if (diff.norm() > SIMD_EPSILON)
		{
			btVector3 dir = diff.normalized();
			f = k * v.dot(dir) * dir;
		}
		force[n->index] -= f;
	}
	// Should add a rotation to preserve the grab orientation, but I did not observe any positive effects even with high values of m_kRot
	//if (m_kRot > btScalar(0) /*&& theta > btScalar(1e-4)*/)
	//{
	//	btVector3 torque = m_kRot * theta * axis;
	//	applyTorqueSpring(torque, force);
	//}
}

// Jacobian‑vector products and energies are kept identical to the previous
// version except that they use the *base* stiffness (orientation boosting
// does not affect the linearised system for stability).

void btDeformableMousePickingForce::addScaledDampingForceDifferential(btScalar scale, const TVStack& dv, TVStack& df)
{
	//btScalar theta;
	//btVector3 axis;
	//orientationError(theta, axis);
	//btScalar orientationBasedStiffnessBoost = transBoost(theta);
	btScalar k = m_dampingStiffness * scale /** orientationBasedStiffnessBoost*/;  // keep linear
	for (int i = 0; i < getIndexCount(); ++i)
	{
		const btSoftBody::Node* n = getNode(i);
		btVector3 diff = (n->m_x - m_mouse_transform.getOrigin()) - m_node_to_mouse_x[i];
		btVector3 local = k * dv[n->index];
		if (diff.length() > SIMD_EPSILON)
		{
			btVector3 dir = diff.normalized();
			local = k * dv[n->index].dot(dir) * dir;
		}
		df[n->index] -= local;
	}
}
void btDeformableMousePickingForce::buildDampingForceDifferentialDiagonal(btScalar, TVStack&) {}

double btDeformableMousePickingForce::totalElasticEnergy(btScalar)
{
	//btScalar theta;
	//btVector3 axis;
	//orientationError(theta, axis);
	//btScalar orientationBasedStiffnessBoost = transBoost(theta);
	double e = 0;
	for (int i = 0; i < getIndexCount(); ++i)
	{
		btVector3 diff = (getNode(i)->m_q - m_mouse_transform.getOrigin()) - m_node_to_mouse_q[i];
		btVector3 f = m_elasticStiffness * diff /** orientationBasedStiffnessBoost*/;
		if (f.safeNorm() > m_maxForce)
		{
			f.safeNormalize();
			f *= m_maxForce;
		}
		e += 0.5 * f.dot(diff);
	}
	return e;
}
double btDeformableMousePickingForce::totalDampingEnergy(btScalar dt)
{
	double e = 0;
	//btScalar theta;
	//btVector3 axis;
	//orientationError(theta, axis);
	//btScalar orientationBasedStiffnessBoost = transBoost(theta);
	for (int i = 0; i < getIndexCount(); ++i)
	{
		const btSoftBody::Node* n = getNode(i);
		btVector3 diff = (n->m_x - m_mouse_transform.getOrigin()) - m_node_to_mouse_x[i];
		btVector3 v = n->m_v;
		btVector3 f = m_dampingStiffness * v /** orientationBasedStiffnessBoost*/;
		if (diff.norm() > SIMD_EPSILON)
		{
			btVector3 dir = diff.normalized();
			f = m_dampingStiffness * v.dot(dir) * dir /** orientationBasedStiffnessBoost*/;
		}
		e -= f.dot(v) / dt;
	}
	return e;
}

void btDeformableMousePickingForce::addScaledElasticForceDifferential(btScalar scale, const TVStack& dx, TVStack& df)
{
	//btScalar theta;
	//btVector3 axis;
	//orientationError(theta, axis);
	//btScalar orientationBasedStiffnessBoost = transBoost(theta);
	btScalar k = scale * m_elasticStiffness /** orientationBasedStiffnessBoost*/;  // linear diff uses base k
	for (int i = 0; i < getIndexCount(); ++i)
	{
		const btSoftBody::Node* n = getNode(i);
		btVector3 diff = (n->m_q - m_mouse_transform.getOrigin()) - m_node_to_mouse_q[i];
		btScalar len = diff.length();
		btVector3 dir = (len > SIMD_EPSILON) ? diff / len : btVector3(0, 0, 0);
		int id = n->index;
		btVector3 dx_i = dx[id];
		btVector3 local(0, 0, 0);
		if (len > SIMD_EPSILON)
		{
			btScalar r = 0;
			local -= k * dir.dot(dx_i) * dir;
			local += k * dir.dot(dx_i) * ((len - r) / len) * dir;
			local -= k * ((len - r) / len) * dx_i;
		}
		df[id] += local;
	}
}

// ------------------------------------------------ setters ---------------
void btDeformableMousePickingForce::setMousePos(const btVector3& p) { m_mouse_transform.setOrigin(p); }
void btDeformableMousePickingForce::setMouseTransform(const btTransform& t)
{
	m_mouse_transform = t;
	rotateNodeToMouse();
}
void btDeformableMousePickingForce::setMaxForce(btScalar f) { m_maxForce = f; }
void btDeformableMousePickingForce::setElasticStiffness(btScalar k) { m_elasticStiffness = k; }
void btDeformableMousePickingForce::setDampingStiffness(btScalar k) { m_dampingStiffness = k; }

btDeformableLagrangianForceType btDeformableMousePickingForce::getForceType() { return BT_MOUSE_PICKING_FORCE; }
