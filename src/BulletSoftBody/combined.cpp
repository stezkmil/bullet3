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

#ifndef BT_CG_PROJECTION_H
#define BT_CG_PROJECTION_H

#include "btSoftBody.h"
#include "BulletDynamics/Featherstone/btMultiBodyLinkCollider.h"
#include "BulletDynamics/Featherstone/btMultiBodyConstraint.h"

struct DeformableContactConstraint
{
	const btSoftBody::Node* m_node;
	btAlignedObjectArray<const btSoftBody::RContact*> m_contact;
	btAlignedObjectArray<btVector3> m_total_normal_dv;
	btAlignedObjectArray<btVector3> m_total_tangent_dv;
	btAlignedObjectArray<bool> m_static;
	btAlignedObjectArray<bool> m_can_be_dynamic;

	DeformableContactConstraint(const btSoftBody::RContact& rcontact) : m_node(rcontact.m_node)
	{
		append(rcontact);
	}

	DeformableContactConstraint() : m_node(NULL)
	{
		m_contact.push_back(NULL);
	}

	void append(const btSoftBody::RContact& rcontact)
	{
		m_contact.push_back(&rcontact);
		m_total_normal_dv.push_back(btVector3(0, 0, 0));
		m_total_tangent_dv.push_back(btVector3(0, 0, 0));
		m_static.push_back(false);
		m_can_be_dynamic.push_back(true);
	}

	void replace(const btSoftBody::RContact& rcontact)
	{
		m_contact.clear();
		m_total_normal_dv.clear();
		m_total_tangent_dv.clear();
		m_static.clear();
		m_can_be_dynamic.clear();
		append(rcontact);
	}

	~DeformableContactConstraint()
	{
	}
};

class btCGProjection
{
public:
	typedef btAlignedObjectArray<btVector3> TVStack;
	typedef btAlignedObjectArray<btAlignedObjectArray<btVector3> > TVArrayStack;
	typedef btAlignedObjectArray<btAlignedObjectArray<btScalar> > TArrayStack;
	btAlignedObjectArray<btSoftBody*>& m_softBodies;
	const btScalar& m_dt;
	// map from node indices to node pointers
	const btAlignedObjectArray<btSoftBody::Node*>* m_nodes;

	btCGProjection(btAlignedObjectArray<btSoftBody*>& softBodies, const btScalar& dt)
		: m_softBodies(softBodies), m_dt(dt)
	{
	}

	virtual ~btCGProjection()
	{
	}

	// apply the constraints
	virtual void project(TVStack& x) = 0;

	virtual void setConstraints() = 0;

	// update the constraints
	virtual btScalar update() = 0;

	virtual void reinitialize(bool nodeUpdated)
	{
	}

	virtual void setIndices(const btAlignedObjectArray<btSoftBody::Node*>* nodes)
	{
		m_nodes = nodes;
	}
};

#endif /* btCGProjection_h */
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

#ifndef BT_CONJUGATE_GRADIENT_H
#define BT_CONJUGATE_GRADIENT_H
#include "btKrylovSolver.h"
template <class MatrixX>
class btConjugateGradient : public btKrylovSolver<MatrixX>
{
	typedef btAlignedObjectArray<btVector3> TVStack;
	typedef btKrylovSolver<MatrixX> Base;
	TVStack r, p, z, temp;

public:
	btConjugateGradient(const int max_it_in)
		: btKrylovSolver<MatrixX>(max_it_in, SIMD_EPSILON)
	{
	}

	virtual ~btConjugateGradient() {}

	// return the number of iterations taken
	int solve(MatrixX& A, TVStack& x, const TVStack& b, bool verbose = false)
	{
		BT_PROFILE("CGSolve");
		btAssert(x.size() == b.size());
		reinitialize(b);
		temp = b;
		A.project(temp);
		p = temp;
		A.precondition(p, z);
		btScalar d0 = this->dot(z, temp);
		d0 = btMin(btScalar(1), d0);
		// r = b - A * x --with assigned dof zeroed out
		A.multiply(x, temp);
		r = this->sub(b, temp);
		A.project(r);
		// z = M^(-1) * r
		A.precondition(r, z);
		A.project(z);
		btScalar r_dot_z = this->dot(z, r);
		if (r_dot_z <= Base::m_tolerance * d0)
		{
			if (verbose)
			{
				std::cout << "Iteration = 0" << std::endl;
				std::cout << "Two norm of the residual = " << r_dot_z << std::endl;
			}
			return 0;
		}
		p = z;
		btScalar r_dot_z_new = r_dot_z;
		for (int k = 1; k <= Base::m_maxIterations; k++)
		{
			// temp = A*p
			A.multiply(p, temp);
			A.project(temp);
			if (this->dot(p, temp) < 0)
			{
				if (verbose)
					std::cout << "Encountered negative direction in CG!" << std::endl;
				if (k == 1)
				{
					x = b;
				}
				return k;
			}
			// alpha = r^T * z / (p^T * A * p)
			btScalar alpha = r_dot_z_new / this->dot(p, temp);
			//  x += alpha * p;
			this->multAndAddTo(alpha, p, x);
			//  r -= alpha * temp;
			this->multAndAddTo(-alpha, temp, r);
			// z = M^(-1) * r
			A.precondition(r, z);
			r_dot_z = r_dot_z_new;
			r_dot_z_new = this->dot(r, z);
			if (r_dot_z_new < Base::m_tolerance * d0)
			{
				if (verbose)
				{
					std::cout << "ConjugateGradient iterations " << k << " residual = " << r_dot_z_new << std::endl;
				}
				return k;
			}

			btScalar beta = r_dot_z_new / r_dot_z;
			p = this->multAndAdd(beta, p, z);
		}
		if (verbose)
		{
			std::cout << "ConjugateGradient max iterations reached " << Base::m_maxIterations << " error = " << r_dot_z_new << std::endl;
		}
		return Base::m_maxIterations;
	}

	void reinitialize(const TVStack& b)
	{
		r.resize(b.size());
		p.resize(b.size());
		z.resize(b.size());
		temp.resize(b.size());
	}
};
#endif /* btConjugateGradient_h */
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

#ifndef BT_CONJUGATE_RESIDUAL_H
#define BT_CONJUGATE_RESIDUAL_H
#include "btKrylovSolver.h"

template <class MatrixX>
class btConjugateResidual : public btKrylovSolver<MatrixX>
{
	typedef btAlignedObjectArray<btVector3> TVStack;
	typedef btKrylovSolver<MatrixX> Base;
	TVStack r, p, z, temp_p, temp_r, best_x;
	// temp_r = A*r
	// temp_p = A*p
	// z = M^(-1) * temp_p = M^(-1) * A * p
	btScalar best_r;

public:
	btConjugateResidual(const int max_it_in)
		: Base(max_it_in, 1e-8)
	{
	}

	virtual ~btConjugateResidual() {}

	// return the number of iterations taken
	int solve(MatrixX& A, TVStack& x, const TVStack& b, bool verbose = false)
	{
		BT_PROFILE("CRSolve");
		btAssert(x.size() == b.size());
		reinitialize(b);
		// r = b - A * x --with assigned dof zeroed out
		A.multiply(x, temp_r);  // borrow temp_r here to store A*x
		r = this->sub(b, temp_r);
		// z = M^(-1) * r
		A.precondition(r, z);  // borrow z to store preconditioned r
		r = z;
		btScalar residual_norm = this->norm(r);
		if (residual_norm <= Base::m_tolerance)
		{
			return 0;
		}
		p = r;
		btScalar r_dot_Ar, r_dot_Ar_new;
		// temp_p = A*p
		A.multiply(p, temp_p);
		// temp_r = A*r
		temp_r = temp_p;
		r_dot_Ar = this->dot(r, temp_r);
		for (int k = 1; k <= Base::m_maxIterations; k++)
		{
			// z = M^(-1) * Ap
			A.precondition(temp_p, z);
			// alpha = r^T * A * r / (Ap)^T * M^-1 * Ap)
			btScalar alpha = r_dot_Ar / this->dot(temp_p, z);
			//  x += alpha * p;
			this->multAndAddTo(alpha, p, x);
			//  r -= alpha * z;
			this->multAndAddTo(-alpha, z, r);
			btScalar norm_r = this->norm(r);
			if (norm_r < best_r)
			{
				best_x = x;
				best_r = norm_r;
				if (norm_r < Base::m_tolerance)
				{
					return k;
				}
			}
			// temp_r = A * r;
			A.multiply(r, temp_r);
			r_dot_Ar_new = this->dot(r, temp_r);
			btScalar beta = r_dot_Ar_new / r_dot_Ar;
			r_dot_Ar = r_dot_Ar_new;
			// p = beta*p + r;
			p = this->multAndAdd(beta, p, r);
			// temp_p = beta*temp_p + temp_r;
			temp_p = this->multAndAdd(beta, temp_p, temp_r);
		}
		if (verbose)
		{
			std::cout << "ConjugateResidual max iterations reached, residual = " << best_r << std::endl;
		}
		x = best_x;
		return Base::m_maxIterations;
	}

	void reinitialize(const TVStack& b)
	{
		r.resize(b.size());
		p.resize(b.size());
		z.resize(b.size());
		temp_p.resize(b.size());
		temp_r.resize(b.size());
		best_x.resize(b.size());
		best_r = SIMD_INFINITY;
	}
};
#endif /* btConjugateResidual_h */
/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/

#include "BulletCollision/CollisionShapes/btTriangleIndexVertexArray.h"
#include "BulletCollision/CollisionDispatch/btCollisionObject.h"
#include "BulletCollision/CollisionShapes/btCollisionShape.h"

#include "btDefaultSoftBodySolver.h"
#include "BulletCollision/CollisionShapes/btCapsuleShape.h"
#include "BulletSoftBody/btSoftBody.h"

btDefaultSoftBodySolver::btDefaultSoftBodySolver()
{
	// Initial we will clearly need to update solver constants
	// For now this is global for the cloths linked with this solver - we should probably make this body specific
	// for performance in future once we understand more clearly when constants need to be updated
	m_updateSolverConstants = true;
}

btDefaultSoftBodySolver::~btDefaultSoftBodySolver()
{
}

// In this case the data is already in the soft bodies so there is no need for us to do anything
void btDefaultSoftBodySolver::copyBackToSoftBodies(bool bMove)
{
}

void btDefaultSoftBodySolver::optimize(btAlignedObjectArray<btSoftBody *> &softBodies, bool forceUpdate)
{
	m_softBodySet.copyFromArray(softBodies);
}

void btDefaultSoftBodySolver::updateSoftBodies()
{
	for (int i = 0; i < m_softBodySet.size(); i++)
	{
		btSoftBody *psb = (btSoftBody *)m_softBodySet[i];
		if (psb->isActive() && !psb->isStaticObject())
		{
			psb->integrateMotion();
		}
	}
}  // updateSoftBodies

bool btDefaultSoftBodySolver::checkInitialized()
{
	return true;
}

void btDefaultSoftBodySolver::solveConstraints(btScalar solverdt)
{
	// Solve constraints for non-solver softbodies
	for (int i = 0; i < m_softBodySet.size(); ++i)
	{
		btSoftBody *psb = static_cast<btSoftBody *>(m_softBodySet[i]);
		if (psb->isActive() && !psb->isStaticObject())
		{
			psb->solveConstraints();
		}
	}
}  // btDefaultSoftBodySolver::solveConstraints

void btDefaultSoftBodySolver::copySoftBodyToVertexBuffer(const btSoftBody *const softBody, btVertexBufferDescriptor *vertexBuffer)
{
	// Currently only support CPU output buffers
	// TODO: check for DX11 buffers. Take all offsets into the same DX11 buffer
	// and use them together on a single kernel call if possible by setting up a
	// per-cloth target buffer array for the copy kernel.

	if (vertexBuffer->getBufferType() == btVertexBufferDescriptor::CPU_BUFFER)
	{
		const btAlignedObjectArray<btSoftBody::Node> &clothVertices(softBody->m_nodes);
		int numVertices = clothVertices.size();

		const btCPUVertexBufferDescriptor *cpuVertexBuffer = static_cast<btCPUVertexBufferDescriptor *>(vertexBuffer);
		float *basePointer = cpuVertexBuffer->getBasePointer();

		if (vertexBuffer->hasVertexPositions())
		{
			const int vertexOffset = cpuVertexBuffer->getVertexOffset();
			const int vertexStride = cpuVertexBuffer->getVertexStride();
			float *vertexPointer = basePointer + vertexOffset;

			for (int vertexIndex = 0; vertexIndex < numVertices; ++vertexIndex)
			{
				btVector3 position = clothVertices[vertexIndex].m_x;
				*(vertexPointer + 0) = (float)position.getX();
				*(vertexPointer + 1) = (float)position.getY();
				*(vertexPointer + 2) = (float)position.getZ();
				vertexPointer += vertexStride;
			}
		}
		if (vertexBuffer->hasNormals())
		{
			const int normalOffset = cpuVertexBuffer->getNormalOffset();
			const int normalStride = cpuVertexBuffer->getNormalStride();
			float *normalPointer = basePointer + normalOffset;

			for (int vertexIndex = 0; vertexIndex < numVertices; ++vertexIndex)
			{
				btVector3 normal = clothVertices[vertexIndex].m_n;
				*(normalPointer + 0) = (float)normal.getX();
				*(normalPointer + 1) = (float)normal.getY();
				*(normalPointer + 2) = (float)normal.getZ();
				normalPointer += normalStride;
			}
		}
	}
}  // btDefaultSoftBodySolver::copySoftBodyToVertexBuffer

void btDefaultSoftBodySolver::processCollision(btSoftBody *softBody, btSoftBody *otherSoftBody)
{
	softBody->defaultCollisionHandler(otherSoftBody);
}

void btDefaultSoftBodySolver::processCollision(btSoftBody *softBody, btSoftBody *otherSoftBody, btManifoldResultForSkin *)
{
	softBody->defaultCollisionHandler(otherSoftBody);
}

// For the default solver just leave the soft body to do its collision processing
void btDefaultSoftBodySolver::processCollision(btSoftBody *softBody, const btCollisionObjectWrapper *collisionObjectWrap, btManifoldResultForSkin *resultOut)
{
	softBody->defaultCollisionHandler(collisionObjectWrap);
}  // btDefaultSoftBodySolver::processCollision

void btDefaultSoftBodySolver::predictMotion(btScalar timeStep)
{
	for (int i = 0; i < m_softBodySet.size(); ++i)
	{
		btSoftBody *psb = m_softBodySet[i];

		if (psb->isActive() && !psb->isStaticObject())
		{
			psb->predictMotion(timeStep);
		}
	}
}
/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/

#ifndef BT_SOFT_BODY_DEFAULT_SOLVER_H
#define BT_SOFT_BODY_DEFAULT_SOLVER_H

#include "BulletSoftBody/btSoftBodySolvers.h"
#include "btSoftBodySolverVertexBuffer.h"
struct btCollisionObjectWrapper;

class btDefaultSoftBodySolver : public btSoftBodySolver
{
protected:
	/** Variable to define whether we need to update solver constants on the next iteration */
	bool m_updateSolverConstants;

	btAlignedObjectArray<btSoftBody *> m_softBodySet;

public:
	btDefaultSoftBodySolver();

	virtual ~btDefaultSoftBodySolver();

	virtual SolverTypes getSolverType() const
	{
		return DEFAULT_SOLVER;
	}

	virtual bool checkInitialized();

	virtual void updateSoftBodies();

	virtual void optimize(btAlignedObjectArray<btSoftBody *> &softBodies, bool forceUpdate = false);

	virtual void copyBackToSoftBodies(bool bMove = true);

	virtual void solveConstraints(btScalar solverdt);

	virtual void predictMotion(btScalar solverdt);

	virtual void copySoftBodyToVertexBuffer(const btSoftBody *const softBody, btVertexBufferDescriptor *vertexBuffer);

	virtual void processCollision(btSoftBody *, const btCollisionObjectWrapper *, btManifoldResultForSkin *);

	virtual void processCollision(btSoftBody *, btSoftBody *);

	virtual void processCollision(btSoftBody *, btSoftBody *, btManifoldResultForSkin *);
};

#endif  // #ifndef BT_ACCELERATED_SOFT_BODY_CPU_SOLVER_H
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

#include "btDeformableBackwardEulerObjective.h"
#include "btPreconditioner.h"
#include "LinearMath/btQuickprof.h"

btDeformableBackwardEulerObjective::btDeformableBackwardEulerObjective(btAlignedObjectArray<btSoftBody*>& softBodies, const TVStack& backup_v)
	: m_softBodies(softBodies), m_projection(softBodies), m_backupVelocity(backup_v), m_implicit(false)
{
	m_massPreconditioner = new MassPreconditioner(m_softBodies);
	m_KKTPreconditioner = new KKTPreconditioner(m_softBodies, m_projection, m_lf, m_dt, m_implicit);
	m_preconditioner = m_KKTPreconditioner;
}

btDeformableBackwardEulerObjective::~btDeformableBackwardEulerObjective()
{
	delete m_KKTPreconditioner;
	delete m_massPreconditioner;
}

void btDeformableBackwardEulerObjective::reinitialize(bool nodeUpdated, btScalar dt)
{
	BT_PROFILE("reinitialize");
	if (dt > 0)
	{
		setDt(dt);
	}
	if (nodeUpdated)
	{
		updateId();
	}
	for (int i = 0; i < m_lf.size(); ++i)
	{
		m_lf[i]->reinitialize(nodeUpdated);
	}
	btMatrix3x3 I;
	I.setIdentity();
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		for (int j = 0; j < psb->m_nodes.size(); ++j)
		{
			if (psb->m_nodes[j].m_im > 0)
				psb->m_nodes[j].m_effectiveMass = I * (1.0 / psb->m_nodes[j].m_im);
		}
	}
	m_projection.reinitialize(nodeUpdated);
	//    m_preconditioner->reinitialize(nodeUpdated);
}

void btDeformableBackwardEulerObjective::setDt(btScalar dt)
{
	m_dt = dt;
}

void btDeformableBackwardEulerObjective::multiply(const TVStack& x, TVStack& b) const
{
	BT_PROFILE("multiply");
	// add in the mass term
	size_t counter = 0;
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		for (int j = 0; j < psb->m_nodes.size(); ++j)
		{
			const btSoftBody::Node& node = psb->m_nodes[j];
			b[counter] = (node.m_im == 0) ? btVector3(0, 0, 0) : x[counter] / node.m_im;
			++counter;
		}
	}

	for (int i = 0; i < m_lf.size(); ++i)
	{
		// add damping matrix
		m_lf[i]->addScaledDampingForceDifferential(-m_dt, x, b);
		// Always integrate picking force implicitly for stability.
		if (m_implicit || m_lf[i]->getForceType() == BT_MOUSE_PICKING_FORCE)
		{
			m_lf[i]->addScaledElasticForceDifferential(-m_dt * m_dt, x, b);
		}
	}
	int offset = m_nodes.size();
	for (int i = offset; i < b.size(); ++i)
	{
		b[i].setZero();
	}
	// add in the lagrange multiplier terms

	for (int c = 0; c < m_projection.m_lagrangeMultipliers.size(); ++c)
	{
		// C^T * lambda
		const LagrangeMultiplier& lm = m_projection.m_lagrangeMultipliers[c];
		for (int i = 0; i < lm.m_num_nodes; ++i)
		{
			for (int j = 0; j < lm.m_num_constraints; ++j)
			{
				b[lm.m_indices[i]] += x[offset + c][j] * lm.m_weights[i] * lm.m_dirs[j];
			}
		}
		// C * x
		for (int d = 0; d < lm.m_num_constraints; ++d)
		{
			for (int i = 0; i < lm.m_num_nodes; ++i)
			{
				b[offset + c][d] += lm.m_weights[i] * x[lm.m_indices[i]].dot(lm.m_dirs[d]);
			}
		}
	}
}

void btDeformableBackwardEulerObjective::updateVelocity(const TVStack& dv)
{
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		for (int j = 0; j < psb->m_nodes.size(); ++j)
		{
			btSoftBody::Node& node = psb->m_nodes[j];
			node.m_v = m_backupVelocity[node.index] + dv[node.index];
		}
	}
}

void btDeformableBackwardEulerObjective::applyForce(TVStack& force, bool setZero)
{
	size_t counter = 0;
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		if (!psb->isActive() || psb->isStaticObject())
		{
			counter += psb->m_nodes.size();
			continue;
		}
		if (m_implicit)
		{
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				if (psb->m_nodes[j].m_im != 0)
				{
					psb->m_nodes[j].m_v += psb->m_nodes[j].m_effectiveMass_inv * force[counter++];
				}
			}
		}
		else
		{
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				btScalar one_over_mass = (psb->m_nodes[j].m_im == 0) ? 0 : psb->m_nodes[j].m_im;
				psb->m_nodes[j].m_v += one_over_mass * force[counter++];
			}
		}
	}
	if (setZero)
	{
		for (int i = 0; i < force.size(); ++i)
			force[i].setZero();
	}
}

void btDeformableBackwardEulerObjective::computeResidual(btScalar dt, TVStack& residual)
{
	BT_PROFILE("computeResidual");
	// add implicit force
	for (int i = 0; i < m_lf.size(); ++i)
	{
		// Always integrate picking force implicitly for stability.
		if (m_implicit || m_lf[i]->getForceType() == BT_MOUSE_PICKING_FORCE)
		{
			m_lf[i]->addScaledForces(dt, residual);
		}
		else
		{
			m_lf[i]->addScaledDampingForce(dt, residual);
		}
	}
	//    m_projection.project(residual);
}

btScalar btDeformableBackwardEulerObjective::computeNorm(const TVStack& residual) const
{
	btScalar mag = 0;
	for (int i = 0; i < residual.size(); ++i)
	{
		mag += residual[i].length2();
	}
	return std::sqrt(mag);
}

btScalar btDeformableBackwardEulerObjective::totalEnergy(btScalar dt)
{
	btScalar e = 0;
	for (int i = 0; i < m_lf.size(); ++i)
	{
		e += m_lf[i]->totalEnergy(dt);
	}
	return e;
}

void btDeformableBackwardEulerObjective::applyExplicitForce(TVStack& force)
{
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		m_softBodies[i]->advanceDeformation();
	}
	if (m_implicit)
	{
		// apply forces except gravity force
		btVector3 gravity;
		for (int i = 0; i < m_lf.size(); ++i)
		{
			if (m_lf[i]->getForceType() == BT_GRAVITY_FORCE)
			{
				gravity = static_cast<btDeformableGravityForce*>(m_lf[i])->m_gravity;
			}
			else
			{
				m_lf[i]->addScaledForces(m_dt, force);
			}
		}
		for (int i = 0; i < m_lf.size(); ++i)
		{
			m_lf[i]->addScaledHessian(m_dt);
		}
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (psb->isActive() && !psb->isStaticObject())
			{
				for (int j = 0; j < psb->m_nodes.size(); ++j)
				{
					// add gravity explicitly
					psb->m_nodes[j].m_v += m_dt * psb->m_gravityFactor * gravity;
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < m_lf.size(); ++i)
		{
			m_lf[i]->addScaledExplicitForce(m_dt, force);
		}
	}
	// calculate inverse mass matrix for all nodes
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		if (psb->isActive() && !psb->isStaticObject())
		{
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				if (psb->m_nodes[j].m_im > 0)
				{
					psb->m_nodes[j].m_effectiveMass_inv = psb->m_nodes[j].m_effectiveMass.inverse();
				}
			}
		}
	}
	applyForce(force, true);
}

void btDeformableBackwardEulerObjective::initialGuess(TVStack& dv, const TVStack& residual)
{
	size_t counter = 0;
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		for (int j = 0; j < psb->m_nodes.size(); ++j)
		{
			dv[counter] = psb->m_nodes[j].m_im * residual[counter];
			++counter;
		}
	}
}

//set constraints as projections
void btDeformableBackwardEulerObjective::setConstraints(const btContactSolverInfo& infoGlobal)
{
	m_projection.setConstraints(infoGlobal);
}

void btDeformableBackwardEulerObjective::applyDynamicFriction(TVStack& r)
{
	m_projection.applyDynamicFriction(r);
}
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

#ifndef BT_BACKWARD_EULER_OBJECTIVE_H
#define BT_BACKWARD_EULER_OBJECTIVE_H
//#include "btConjugateGradient.h"
#include "btDeformableLagrangianForce.h"
#include "btDeformableMassSpringForce.h"
#include "btDeformableGravityForce.h"
#include "btDeformableCorotatedForce.h"
#include "btDeformableMousePickingForce.h"
#include "btDeformableLinearElasticityForce.h"
#include "btDeformableNeoHookeanForce.h"
#include "btDeformableContactProjection.h"
#include "btPreconditioner.h"
// #include "btDeformableMultiBodyDynamicsWorld.h"
#include "LinearMath/btQuickprof.h"

class btDeformableBackwardEulerObjective
{
public:
	enum _
	{
		Mass_preconditioner,
		KKT_preconditioner
	};

	typedef btAlignedObjectArray<btVector3> TVStack;
	btScalar m_dt;
	btAlignedObjectArray<btDeformableLagrangianForce*> m_lf;
	btAlignedObjectArray<btSoftBody*>& m_softBodies;
	Preconditioner* m_preconditioner;
	btDeformableContactProjection m_projection;
	const TVStack& m_backupVelocity;
	btAlignedObjectArray<btSoftBody::Node*> m_nodes;
	bool m_implicit;
	MassPreconditioner* m_massPreconditioner;
	KKTPreconditioner* m_KKTPreconditioner;

	btDeformableBackwardEulerObjective(btAlignedObjectArray<btSoftBody*>& softBodies, const TVStack& backup_v);

	virtual ~btDeformableBackwardEulerObjective();

	void initialize() {}

	// compute the rhs for CG solve, i.e, add the dt scaled implicit force to residual
	void computeResidual(btScalar dt, TVStack& residual);

	// add explicit force to the velocity
	void applyExplicitForce(TVStack& force);

	// apply force to velocity and optionally reset the force to zero
	void applyForce(TVStack& force, bool setZero);

	// compute the norm of the residual
	btScalar computeNorm(const TVStack& residual) const;

	// compute one step of the solve (there is only one solve if the system is linear)
	void computeStep(TVStack& dv, const TVStack& residual, const btScalar& dt);

	// perform A*x = b
	void multiply(const TVStack& x, TVStack& b) const;

	// set initial guess for CG solve
	void initialGuess(TVStack& dv, const TVStack& residual);

	// reset data structure and reset dt
	void reinitialize(bool nodeUpdated, btScalar dt);

	void setDt(btScalar dt);

	// add friction force to residual
	void applyDynamicFriction(TVStack& r);

	// add dv to velocity
	void updateVelocity(const TVStack& dv);

	//set constraints as projections
	void setConstraints(const btContactSolverInfo& infoGlobal);

	// update the projections and project the residual
	void project(TVStack& r)
	{
		BT_PROFILE("project");
		m_projection.project(r);
	}

	// perform precondition M^(-1) x = b
	void precondition(const TVStack& x, TVStack& b)
	{
		m_preconditioner->operator()(x, b);
	}

	// reindex all the vertices
	virtual void updateId()
	{
		size_t node_id = 0;
		size_t face_id = 0;
		m_nodes.clear();
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				psb->m_nodes[j].index = node_id;
				m_nodes.push_back(&psb->m_nodes[j]);
				++node_id;
			}
			for (int j = 0; j < psb->m_faces.size(); ++j)
			{
				psb->m_faces[j].m_index = face_id;
				++face_id;
			}
		}
	}

	const btAlignedObjectArray<btSoftBody::Node*>* getIndices() const
	{
		return &m_nodes;
	}

	void setImplicit(bool implicit)
	{
		m_implicit = implicit;
	}

	// Calculate the total potential energy in the system
	btScalar totalEnergy(btScalar dt);

	void addLagrangeMultiplier(const TVStack& vec, TVStack& extended_vec)
	{
		extended_vec.resize(vec.size() + m_projection.m_lagrangeMultipliers.size());
		for (int i = 0; i < vec.size(); ++i)
		{
			extended_vec[i] = vec[i];
		}
		int offset = vec.size();
		for (int i = 0; i < m_projection.m_lagrangeMultipliers.size(); ++i)
		{
			extended_vec[offset + i].setZero();
		}
	}

	void addLagrangeMultiplierRHS(const TVStack& residual, const TVStack& m_dv, TVStack& extended_residual)
	{
		extended_residual.resize(residual.size() + m_projection.m_lagrangeMultipliers.size());
		for (int i = 0; i < residual.size(); ++i)
		{
			extended_residual[i] = residual[i];
		}
		int offset = residual.size();
		for (int i = 0; i < m_projection.m_lagrangeMultipliers.size(); ++i)
		{
			const LagrangeMultiplier& lm = m_projection.m_lagrangeMultipliers[i];
			extended_residual[offset + i].setZero();
			for (int d = 0; d < lm.m_num_constraints; ++d)
			{
				for (int n = 0; n < lm.m_num_nodes; ++n)
				{
					extended_residual[offset + i][d] += lm.m_weights[n] * m_dv[lm.m_indices[n]].dot(lm.m_dirs[d]);
				}
			}
		}
	}

	void calculateContactForce(const TVStack& dv, const TVStack& rhs, TVStack& f)
	{
		size_t counter = 0;
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				const btSoftBody::Node& node = psb->m_nodes[j];
				f[counter] = (node.m_im == 0) ? btVector3(0, 0, 0) : dv[counter] / node.m_im;
				++counter;
			}
		}
		for (int i = 0; i < m_lf.size(); ++i)
		{
			// add damping matrix
			m_lf[i]->addScaledDampingForceDifferential(-m_dt, dv, f);
		}
		counter = 0;
		for (; counter < f.size(); ++counter)
		{
			f[counter] = rhs[counter] - f[counter];
		}
	}
};

#endif /* btBackwardEulerObjective_h */
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

#include <stdio.h>
#include <limits>
#include "btDeformableBodySolver.h"
#include "btSoftBodyInternals.h"
#include "LinearMath/btQuickprof.h"
static const int kMaxConjugateGradientIterations = 300;
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
			computeStep(x, rhs);
			for (int i = 0; i < m_dv.size(); ++i)
			{
				m_dv[i] = x[i];
			}
		}
		updateVelocity();
	}
	else
	{
		for (int i = 0; i < m_maxNewtonIterations; ++i)
		{
			updateState();
			// add the inertia term in the residual
			int counter = 0;
			for (int k = 0; k < m_softBodies.size(); ++k)
			{
				btSoftBody* psb = m_softBodies[k];
				for (int j = 0; j < psb->m_nodes.size(); ++j)
				{
					if (psb->m_nodes[j].m_im > 0)
					{
						m_residual[counter] = (-1. / psb->m_nodes[j].m_im) * m_dv[counter];
					}
					++counter;
				}
			}

			m_objective->computeResidual(solverdt, m_residual);
			if (m_objective->computeNorm(m_residual) < m_newtonTolerance && i > 0)
			{
				break;
			}
			// todo xuchenhan@: this really only needs to be calculated once
			m_objective->applyDynamicFriction(m_residual);
			if (m_lineSearch)
			{
				btScalar inner_product = computeDescentStep(m_ddv, m_residual);
				btScalar alpha = 0.01, beta = 0.5;  // Boyd & Vandenberghe suggested alpha between 0.01 and 0.3, beta between 0.1 to 0.8
				btScalar scale = 2;
				btScalar f0 = m_objective->totalEnergy(solverdt) + kineticEnergy(), f1, f2;
				backupDv();
				do
				{
					scale *= beta;
					if (scale < 1e-8)
					{
						return;
					}
					updateEnergy(scale);
					f1 = m_objective->totalEnergy(solverdt) + kineticEnergy();
					f2 = f0 - alpha * scale * inner_product;
				} while (!(f1 < f2 + SIMD_EPSILON));  // if anything here is nan then the search continues
				revertDv();
				updateDv(scale);
			}
			else
			{
				computeStep(m_ddv, m_residual);
				updateDv();
			}
			for (int j = 0; j < m_numNodes; ++j)
			{
				m_ddv[j].setZero();
				m_residual[j].setZero();
			}
		}
		updateVelocity();
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
			if (node.m_im > 0)
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
	m_cg.solve(*m_objective, ddv, residual, false);
	btScalar inner_product = m_cg.dot(residual, m_ddv);
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

void btDeformableBodySolver::computeStep(TVStack& ddv, const TVStack& residual)
{
	if (m_useProjection)
		m_cg.solve(*m_objective, ddv, residual, false);
	else
		m_cr.solve(*m_objective, ddv, residual, false);
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
				// setting the initial guess for newton, need m_dv = v_{n+1} - v_n for dofs that are in constraint.
				if (psb->m_nodes[j].m_v == m_backupVelocity[counter])
					m_dv[counter].setZero();
				else
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
			if (j == 866)
				fprintf(stderr, "b real pos %f %f %f node.m_v %f %f %f node.m_splitv %f %f %f\n",
                    node.m_x.x(), node.m_x.y(), node.m_x.z(),
                    node.m_v.x(), node.m_v.y(), node.m_v.z(),
                    node.m_splitv.x(), node.m_splitv.y(), node.m_splitv.z());
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

#ifndef BT_DEFORMABLE_BODY_SOLVERS_H
#define BT_DEFORMABLE_BODY_SOLVERS_H

#include "btSoftBodySolvers.h"
#include "btDeformableBackwardEulerObjective.h"
#include "btDeformableMultiBodyDynamicsWorld.h"
#include "BulletDynamics/Featherstone/btMultiBodyLinkCollider.h"
#include "BulletDynamics/Featherstone/btMultiBodyConstraint.h"
#include "btConjugateResidual.h"
#include "btConjugateGradient.h"
struct btCollisionObjectWrapper;
// class btDeformableBackwardEulerObjective;
// class btDeformableMultiBodyDynamicsWorld;

class btDeformableBodySolver : public btSoftBodySolver
{
	typedef btAlignedObjectArray<btVector3> TVStack;

protected:
	int m_numNodes;                                                // total number of deformable body nodes
	TVStack m_dv;                                                  // v_{n+1} - v_n
	TVStack m_backup_dv;                                           // backed up dv
	TVStack m_ddv;                                                 // incremental dv
	TVStack m_residual;                                            // rhs of the linear solve
	btAlignedObjectArray<btSoftBody*> m_softBodies;                // all deformable bodies
	TVStack m_backupVelocity;                                      // backed up v, equals v_n for implicit, equals v_{n+1}^* for explicit
	btScalar m_dt;                                                 // dt
	btConjugateGradient<btDeformableBackwardEulerObjective> m_cg;  // CG solver
	btConjugateResidual<btDeformableBackwardEulerObjective> m_cr;  // CR solver
	bool m_implicit;                                               // use implicit scheme if true, explicit scheme if false
	int m_maxNewtonIterations;                                     // max number of newton iterations
	btScalar m_newtonTolerance;                                    // stop newton iterations if f(x) < m_newtonTolerance
	bool m_lineSearch;                                             // If true, use newton's method with line search under implicit scheme
	bool m_reducedSolver;                                          // flag for reduced soft body solver
public:
	// handles data related to objective function
	btDeformableBackwardEulerObjective* m_objective;
	bool m_useProjection;

	btDeformableBodySolver();

	virtual ~btDeformableBodySolver();

	virtual SolverTypes getSolverType() const
	{
		return DEFORMABLE_SOLVER;
	}

	// update soft body normals
	virtual void updateSoftBodies();

	virtual btScalar solveContactConstraints(btCollisionObject** deformableBodies, int numDeformableBodies, const btContactSolverInfo& infoGlobal);

	// solve the momentum equation
	virtual void solveDeformableConstraints(btScalar solverdt);

	// set gravity (get from deformable world)
	virtual void setGravity(const btVector3& gravity)
	{
		// for full deformable object, we don't store gravity in the solver
		// this function is overriden in the reduced deformable object
	}

	// resize/clear data structures
	virtual void reinitialize(const btAlignedObjectArray<btSoftBody*>& softBodies, btScalar dt);

	// set up contact constraints
	virtual void setConstraints(const btContactSolverInfo& infoGlobal);

	// add in elastic forces and gravity to obtain v_{n+1}^* and calls predictDeformableMotion
	virtual void predictMotion(btScalar solverdt);

	// move to temporary position x_{n+1}^* = x_n + dt * v_{n+1}^*
	// x_{n+1}^* is stored in m_q
	void predictDeformableMotion(btSoftBody* psb, btScalar dt);

	// save the current velocity to m_backupVelocity
	void backupVelocity();

	// set m_dv and m_backupVelocity to desired value to prepare for momentum solve
	virtual void setupDeformableSolve(bool implicit);

	// set the current velocity to that backed up in m_backupVelocity
	void revertVelocity();

	// set velocity to m_dv + m_backupVelocity
	void updateVelocity();

	// update the node count
	bool updateNodes();

	// calculate the change in dv resulting from the momentum solve
	void computeStep(TVStack& ddv, const TVStack& residual);

	// calculate the change in dv resulting from the momentum solve when line search is turned on
	btScalar computeDescentStep(TVStack& ddv, const TVStack& residual, bool verbose = false);

	virtual void copySoftBodyToVertexBuffer(const btSoftBody* const softBody, btVertexBufferDescriptor* vertexBuffer) {}

	// process collision between deformable and rigid
	virtual void processCollision(btSoftBody* softBody, const btCollisionObjectWrapper* collisionObjectWrap, btManifoldResultForSkin* resultOut);

	// process collision between deformable and deformable
	virtual void processCollision(btSoftBody* softBody, btSoftBody* otherSoftBody);

	virtual void processCollision(btSoftBody* softBody, btSoftBody* otherSoftBody, btManifoldResultForSkin* resultOut);

	// If true, implicit time stepping scheme is used.
	// Otherwise, explicit time stepping scheme is used
	void setImplicit(bool implicit);

	// If true, newton's method with line search is used when implicit time stepping scheme is turned on
	void setLineSearch(bool lineSearch);

	// set temporary position x^* = x_n + dt * v
	// update the deformation gradient at position x^*
	void updateState();

	// set dv = dv + scale * ddv
	void updateDv(btScalar scale = 1);

	// set temporary position x^* = x_n + dt * v^*
	void updateTempPosition();

	// save the current dv to m_backup_dv;
	void backupDv();

	// set dv to the backed-up value
	void revertDv();

	// set dv = dv + scale * ddv
	// set v^* = v_n + dv
	// set temporary position x^* = x_n + dt * v^*
	// update the deformation gradient at position x^*
	void updateEnergy(btScalar scale);

	// calculates the appropriately scaled kinetic energy in the system, which is
	// 1/2 * dv^T * M * dv
	// used in line search
	btScalar kineticEnergy();

	// add explicit force to the velocity in the objective class
	virtual void applyExplicitForce();

	// execute position/velocity update and apply anchor constraints in the integrateTransforms from the Dynamics world
	virtual void applyTransforms(btScalar timeStep);

	virtual void setStrainLimiting(bool opt)
	{
		m_objective->m_projection.m_useStrainLimiting = opt;
	}

	virtual void setPreconditioner(int opt)
	{
		switch (opt)
		{
			case btDeformableBackwardEulerObjective::Mass_preconditioner:
				m_objective->m_preconditioner = m_objective->m_massPreconditioner;
				break;

			case btDeformableBackwardEulerObjective::KKT_preconditioner:
				m_objective->m_preconditioner = m_objective->m_KKTPreconditioner;
				break;

			default:
				btAssert(false);
				break;
		}
	}

	virtual btAlignedObjectArray<btDeformableLagrangianForce*>* getLagrangianForceArray()
	{
		return &(m_objective->m_lf);
	}

	virtual const btAlignedObjectArray<btSoftBody::Node*>* getIndices()
	{
		return m_objective->getIndices();
	}

	virtual void setProjection()
	{
		m_objective->m_projection.setProjection();
	}

	virtual void setLagrangeMultiplier()
	{
		m_objective->m_projection.setLagrangeMultiplier();
	}

	virtual bool isReducedSolver()
	{
		return m_reducedSolver;
	}

	virtual void deformableBodyInternalWriteBack() {}

	// unused functions
	virtual void optimize(btAlignedObjectArray<btSoftBody*>& softBodies, bool forceUpdate = false) {}
	virtual void solveConstraints(btScalar dt) {}
	virtual bool checkInitialized() { return true; }
	virtual void copyBackToSoftBodies(bool bMove = true) {}
};

#endif /* btDeformableBodySolver_h */
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

#include "btDeformableContactConstraint.h"
/* ================   Deformable Node Anchor   =================== */
btDeformableNodeAnchorConstraint::btDeformableNodeAnchorConstraint(const btSoftBody::DeformableNodeRigidAnchor& a, const btContactSolverInfo& infoGlobal)
	: m_anchor(&a), btDeformableContactConstraint(a.m_cti.m_normal, infoGlobal)
{
}

btDeformableNodeAnchorConstraint::btDeformableNodeAnchorConstraint(const btDeformableNodeAnchorConstraint& other)
	: m_anchor(other.m_anchor), btDeformableContactConstraint(other)
{
}

btVector3 btDeformableNodeAnchorConstraint::getVa() const
{
	const btSoftBody::sCti& cti = m_anchor->m_cti;
	btVector3 va(0, 0, 0);
	if (cti.m_colObj->hasContactResponse())
	{
		btRigidBody* rigidCol = 0;
		btMultiBodyLinkCollider* multibodyLinkCol = 0;

		// grab the velocity of the rigid body
		if (cti.m_colObj->getInternalType() == btCollisionObject::CO_RIGID_BODY)
		{
			rigidCol = (btRigidBody*)btRigidBody::upcast(cti.m_colObj);
			va = rigidCol ? (rigidCol->getVelocityInLocalPoint(m_anchor->m_c1)) : btVector3(0, 0, 0);
		}
		else if (cti.m_colObj->getInternalType() == btCollisionObject::CO_FEATHERSTONE_LINK)
		{
			multibodyLinkCol = (btMultiBodyLinkCollider*)btMultiBodyLinkCollider::upcast(cti.m_colObj);
			if (multibodyLinkCol)
			{
				const int ndof = multibodyLinkCol->m_multiBody->getNumDofs() + 6;
				const btScalar* J_n = &m_anchor->jacobianData_normal.m_jacobians[0];
				const btScalar* J_t1 = &m_anchor->jacobianData_t1.m_jacobians[0];
				const btScalar* J_t2 = &m_anchor->jacobianData_t2.m_jacobians[0];
				const btScalar* local_v = multibodyLinkCol->m_multiBody->getVelocityVector();
				const btScalar* local_dv = multibodyLinkCol->m_multiBody->getDeltaVelocityVector();
				// add in the normal component of the va
				btScalar vel = 0.0;
				for (int k = 0; k < ndof; ++k)
				{
					vel += (local_v[k] + local_dv[k]) * J_n[k];
				}
				va = cti.m_normal * vel;
				// add in the tangential components of the va
				vel = 0.0;
				for (int k = 0; k < ndof; ++k)
				{
					vel += (local_v[k] + local_dv[k]) * J_t1[k];
				}
				va += m_anchor->t1 * vel;
				vel = 0.0;
				for (int k = 0; k < ndof; ++k)
				{
					vel += (local_v[k] + local_dv[k]) * J_t2[k];
				}
				va += m_anchor->t2 * vel;
			}
		}
	}
	return va;
}

// Velocity matching between the deformable node and the rigid/multibody it is anchored to. This is the main constraint solve function for the deformable node anchor constraint. Positional drift correction is handled in btDeformableNodeAnchorConstraint::solveSplitImpulse.
btScalar btDeformableNodeAnchorConstraint::solveConstraint(const btContactSolverInfo& infoGlobal)
{
	const btSoftBody::sCti& cti = m_anchor->m_cti;
	btVector3 va = getVa();
	btVector3 vb = getVb();
	btVector3 vr = (vb - va);
	const btScalar dn = btDot(vr, vr);
	// dn is the normal component of velocity diffrerence. Approximates the residual. // todo xuchenhan@: this prob needs to be scaled by dt
	btScalar residualSquare = dn * dn;
	btVector3 impulse = m_anchor->m_c0 * vr;
	// apply impulse to deformable nodes involved and change their velocities
	applyImpulse(impulse);

	// apply impulse to the rigid/multibodies involved and change their velocities
	if (cti.m_colObj->getInternalType() == btCollisionObject::CO_RIGID_BODY)
	{
		btRigidBody* rigidCol = 0;
		rigidCol = (btRigidBody*)btRigidBody::upcast(cti.m_colObj);
		if (rigidCol && rigidCol->isActive())
		{
			rigidCol->applyImpulse(impulse, m_anchor->m_c1);
		}
	}
	else if (cti.m_colObj->getInternalType() == btCollisionObject::CO_FEATHERSTONE_LINK)
	{
		btMultiBodyLinkCollider* multibodyLinkCol = 0;
		multibodyLinkCol = (btMultiBodyLinkCollider*)btMultiBodyLinkCollider::upcast(cti.m_colObj);
		if (multibodyLinkCol)
		{
			const btScalar* deltaV_normal = &m_anchor->jacobianData_normal.m_deltaVelocitiesUnitImpulse[0];
			// apply normal component of the impulse
			multibodyLinkCol->m_multiBody->applyDeltaVeeMultiDof2(deltaV_normal, impulse.dot(cti.m_normal));
			// apply tangential component of the impulse
			const btScalar* deltaV_t1 = &m_anchor->jacobianData_t1.m_deltaVelocitiesUnitImpulse[0];
			multibodyLinkCol->m_multiBody->applyDeltaVeeMultiDof2(deltaV_t1, impulse.dot(m_anchor->t1));
			const btScalar* deltaV_t2 = &m_anchor->jacobianData_t2.m_deltaVelocitiesUnitImpulse[0];
			multibodyLinkCol->m_multiBody->applyDeltaVeeMultiDof2(deltaV_t2, impulse.dot(m_anchor->t2));
		}
	}

	return residualSquare;
}

btVector3 btDeformableNodeAnchorConstraint::getVb() const
{
	return m_anchor->m_node->m_v;
}

void btDeformableNodeAnchorConstraint::applyImpulse(const btVector3& impulse)
{
	btVector3 dv = impulse * m_anchor->m_c2;
	m_anchor->m_node->m_v -= dv;
}

btVector3 btDeformableNodeAnchorConstraint::getSplitVa() const

{
	const btSoftBody::sCti& cti = m_anchor->m_cti;
	btVector3 va(0, 0, 0);
	if (cti.m_colObj->hasContactResponse())
	{
		if (cti.m_colObj->getInternalType() == btCollisionObject::CO_RIGID_BODY)
		{
			btRigidBody* rigidCol = (btRigidBody*)btRigidBody::upcast(cti.m_colObj);
			va = rigidCol ? rigidCol->getPushVelocityInLocalPoint(m_anchor->m_c1) : btVector3(0, 0, 0);
		}
		// TODO feather
	}
	return va;
}

btVector3 btDeformableNodeAnchorConstraint::getSplitVb() const
{
	return m_anchor->m_node->m_splitv;
}

void btDeformableNodeAnchorConstraint::applySplitImpulse(const btVector3& impulse)
{
	btVector3 dv = impulse * m_anchor->m_c2;
	m_anchor->m_node->m_splitv -= dv;
}

// This function is to solve the positional drift of the deformable node anchor constraint by applying split impulse. Originally, the drift was
// removed by the n->m_x = a.m_cti.m_colObj->getWorldTransform() * a.m_local; hack in btDeformableBodySolver::applyTransforms. That cause massive destabilisation of the rigid body when the
// soft body stiffness was high, because it caused huge velocity spikes on that soft body node as it was yanked back into place by the elastic forces. These velocity spikes poisoned the rigid body velocity in
// btDeformableNodeAnchorConstraint::solveConstraint.
btScalar btDeformableNodeAnchorConstraint::solveSplitImpulse(const btContactSolverInfo& infoGlobal)
{
	const btSoftBody::sCti& cti = m_anchor->m_cti;
	auto va = getSplitVa() + getVa();
	// Why is m_vn (previous velocity) here? m_v could not be used, because at this stage it contains only the explicit forces (like mouse drag) applied but not the implicit forces like the linear elasticity.
	// To predict the the position of the deformable node, we use vn as an approximation of the velocity at the end of the time step, which includes the effect of both explicit and implicit forces. It is not perfectly
	// accurate, but it works reasonably well due to time coherence. To get perfect velocity, this solve would have to be moved after the implicit forces are calculated. This could be a future TODO - not sure if easy.
	const auto vb = m_anchor->m_node->m_vn + m_anchor->m_node->m_splitv;

	//auto rigid_pt = (cti.m_colObj->getWorldTransform() * m_anchor->m_local);
	//fprintf(stderr, "m_anchor->m_node->m_x %f %f %f cti.m_colObj->getWorldTransform() * m_anchor->m_local %f %f %f diff %f %f %f\n", m_anchor->m_node->m_x.x(), m_anchor->m_node->m_x.y(), m_anchor->m_node->m_x.z(), rigid_pt.x(), rigid_pt.y(), rigid_pt.z(), trigger.x(), trigger.y(), trigger.z());

	// Predicted anchor poisitions for both the rigid and soft are calculated in this part.
	auto anchor_soft_predicted = (m_anchor->m_node->m_x + vb * infoGlobal.m_timeStep);

	btTransform rigid_tr = cti.m_colObj->getWorldTransform();
	btTransform tr_predicted;
	btTransformUtil::integrateTransform(
		rigid_tr,
		m_anchor->m_body->getLinearVelocity() + m_anchor->m_body->getPushVelocity(),
		m_anchor->m_body->getAngularVelocity() + m_anchor->m_body->getTurnVelocity() * infoGlobal.m_splitImpulseTurnErp,
		infoGlobal.m_timeStep,
		tr_predicted);
	btVector3 anchor_rigid_predicted = tr_predicted * m_anchor->m_local;

	btVector3 pos_diff = anchor_soft_predicted - anchor_rigid_predicted;

	const btScalar diff_len_squared = btDot(pos_diff, pos_diff);
	btScalar residualSquare = diff_len_squared * diff_len_squared;

	// Maybe this early return would prevent some slight overshoot?
	//if (residualSquare < infoGlobal.m_leastSquaresResidualThreshold)
	//	return residualSquare;

	/*fprintf(stderr, "anchor_rigid_predicted %f %f %f anchor_soft_predicted %f %f %f pos_diff %f %f %f va %f %f %f vb %f %f %f residual %f\n",
			anchor_rigid_predicted.x(), anchor_rigid_predicted.y(), anchor_rigid_predicted.z(),
			anchor_soft_predicted.x(), anchor_soft_predicted.y(), anchor_soft_predicted.z(),
			pos_diff.x(), pos_diff.y(), pos_diff.z(),
			va.x(), va.y(), va.z(),
			vb.x(), vb.y(), vb.z(),
			residualSquare);*/

	// 0.5 so that the rigid goes halfway and the soft also goes halfway.
	btVector3 impulse = (pos_diff * 0.5) / infoGlobal.m_timeStep;

	// Applies impulse to the rigid body, to minimize the gap between rigid and soft anchor positions.
	if (cti.m_colObj->getInternalType() == btCollisionObject::CO_RIGID_BODY)
	{
		if (m_anchor->m_body)
		{
			// TODO the question of activity. If the activate line is here, the soft never falls asleep. If it isn't there, the rigid sometimes falls asleep even if the soft is moving.
			m_anchor->m_body->activate();
			// Inv mass and linear factor applied to counter their application in applyCentralPushImpulse. Remember - we do not need physically corerct impulse application, we just correct a drift.
			// One pitfall here is to reason that we could perhaps use applyCentralPushImpulse here and get away with it, because anchor is only a positional constraint. It does not converge when multiple
			// anchors are involved. Only when rotations are involved, the gap is reduced on the current anchor and at least partially preserved on the previous anchor.
			m_anchor->m_body->applyPushImpulse((impulse / m_anchor->m_body->getInvMass()) / m_anchor->m_body->getLinearFactor(), m_anchor->m_c1);
		}
	}

	// I thought that I could do it by only pushing the rigid body, but there are problems with that. If soft and rigid are connected by many anchors, then the soft anchored nodes might get slightly squashed
	// by elasticity and the rigid will then never be able to reach all anchor positions with sufficient accuracy, so it will never converge and run full iteration count. To counter that, the soft must also contribute
	// by being impulsed.
	applySplitImpulse(impulse / m_anchor->m_c2);

	// Another pitfall here is to try to perform a final prediction here and return a residual based on that. Not only would that be a redundant code duplication, because that same functionality
	// would be done at the start of the next solveSplitImpulse, but it also causes partial mis-convergence.

	return residualSquare;
}

/* ================   Deformable vs. Rigid   =================== */
btDeformableRigidContactConstraint::btDeformableRigidContactConstraint(const btSoftBody::DeformableRigidContact& c, const btContactSolverInfo& infoGlobal)
	: m_contact(&c), btDeformableContactConstraint(c.m_cti.m_normal, infoGlobal)
{
	m_total_normal_dv.setZero();
	m_total_tangent_dv.setZero();
	// The magnitude of penetration is the depth of penetration.
	m_penetration = c.m_cti.m_offset;
	m_total_split_impulse = 0;
	m_binding = false;
}

btDeformableRigidContactConstraint::btDeformableRigidContactConstraint(const btDeformableRigidContactConstraint& other)
	: m_contact(other.m_contact), btDeformableContactConstraint(other), m_penetration(other.m_penetration), m_total_split_impulse(other.m_total_split_impulse), m_binding(other.m_binding)
{
	m_total_normal_dv = other.m_total_normal_dv;
	m_total_tangent_dv = other.m_total_tangent_dv;
}

btVector3 btDeformableRigidContactConstraint::getVa() const
{
	const btSoftBody::sCti& cti = m_contact->m_cti;
	btVector3 va(0, 0, 0);
	if (cti.m_colObj->hasContactResponse())
	{
		btRigidBody* rigidCol = 0;
		btMultiBodyLinkCollider* multibodyLinkCol = 0;

		// grab the velocity of the rigid body
		if (cti.m_colObj->getInternalType() == btCollisionObject::CO_RIGID_BODY)
		{
			rigidCol = (btRigidBody*)btRigidBody::upcast(cti.m_colObj);
			va = rigidCol ? (rigidCol->getVelocityInLocalPoint(m_contact->m_c1)) : btVector3(0, 0, 0);
		}
		else if (cti.m_colObj->getInternalType() == btCollisionObject::CO_FEATHERSTONE_LINK)
		{
			multibodyLinkCol = (btMultiBodyLinkCollider*)btMultiBodyLinkCollider::upcast(cti.m_colObj);
			if (multibodyLinkCol)
			{
				const int ndof = multibodyLinkCol->m_multiBody->getNumDofs() + 6;
				const btScalar* J_n = &m_contact->jacobianData_normal.m_jacobians[0];
				const btScalar* J_t1 = &m_contact->jacobianData_t1.m_jacobians[0];
				const btScalar* J_t2 = &m_contact->jacobianData_t2.m_jacobians[0];
				const btScalar* local_v = multibodyLinkCol->m_multiBody->getVelocityVector();
				const btScalar* local_dv = multibodyLinkCol->m_multiBody->getDeltaVelocityVector();
				// add in the normal component of the va
				btScalar vel = 0.0;
				for (int k = 0; k < ndof; ++k)
				{
					vel += (local_v[k] + local_dv[k]) * J_n[k];
				}
				va = cti.m_normal * vel;
				// add in the tangential components of the va
				vel = 0.0;
				for (int k = 0; k < ndof; ++k)
				{
					vel += (local_v[k] + local_dv[k]) * J_t1[k];
				}
				va += m_contact->t1 * vel;
				vel = 0.0;
				for (int k = 0; k < ndof; ++k)
				{
					vel += (local_v[k] + local_dv[k]) * J_t2[k];
				}
				va += m_contact->t2 * vel;
			}
		}
	}
	return va;
}

btVector3 btDeformableRigidContactConstraint::getSplitVa() const
{
	const btSoftBody::sCti& cti = m_contact->m_cti;
	btVector3 va(0, 0, 0);
	if (cti.m_colObj->hasContactResponse())
	{
		btRigidBody* rigidCol = 0;
		btMultiBodyLinkCollider* multibodyLinkCol = 0;

		// grab the velocity of the rigid body
		if (cti.m_colObj->getInternalType() == btCollisionObject::CO_RIGID_BODY)
		{
			rigidCol = (btRigidBody*)btRigidBody::upcast(cti.m_colObj);
			va = rigidCol ? (rigidCol->getPushVelocityInLocalPoint(m_contact->m_c1)) : btVector3(0, 0, 0);
		}
		else if (cti.m_colObj->getInternalType() == btCollisionObject::CO_FEATHERSTONE_LINK)
		{
			multibodyLinkCol = (btMultiBodyLinkCollider*)btMultiBodyLinkCollider::upcast(cti.m_colObj);
			if (multibodyLinkCol)
			{
				const int ndof = multibodyLinkCol->m_multiBody->getNumDofs() + 6;
				const btScalar* J_n = &m_contact->jacobianData_normal.m_jacobians[0];
				const btScalar* J_t1 = &m_contact->jacobianData_t1.m_jacobians[0];
				const btScalar* J_t2 = &m_contact->jacobianData_t2.m_jacobians[0];
				const btScalar* local_split_v = multibodyLinkCol->m_multiBody->getSplitVelocityVector();
				// add in the normal component of the va
				btScalar vel = 0.0;
				for (int k = 0; k < ndof; ++k)
				{
					vel += local_split_v[k] * J_n[k];
				}
				va = cti.m_normal * vel;
				// add in the tangential components of the va
				vel = 0.0;
				for (int k = 0; k < ndof; ++k)
				{
					vel += local_split_v[k] * J_t1[k];
				}
				va += m_contact->t1 * vel;
				vel = 0.0;
				for (int k = 0; k < ndof; ++k)
				{
					vel += local_split_v[k] * J_t2[k];
				}
				va += m_contact->t2 * vel;
			}
		}
	}
	return va;
}

btScalar btDeformableRigidContactConstraint::solveConstraint(const btContactSolverInfo& infoGlobal)
{
	const btSoftBody::sCti& cti = m_contact->m_cti;
	btVector3 va = getVa();
	btVector3 vb = getVb();
	btVector3 vr = vb - va;
	btScalar dn = btDot(vr, cti.m_normal) + m_total_normal_dv.dot(cti.m_normal) * infoGlobal.m_deformable_cfm;
	// TODO add some comments trying to describe the m_penetration logic. Maybe there should not be any m_penetration handling if infoGlobal.m_splitImpulse is true? Because then penetrations are completely handled in btDeformableRigidContactConstraint::solveSplitImpulse?
	if (m_penetration < 0)
	{
		dn += m_penetration / infoGlobal.m_timeStep;
	}
	if (!infoGlobal.m_splitImpulse && m_penetration < 0)
	{
		dn += m_penetration * infoGlobal.m_deformable_erp / infoGlobal.m_timeStep;
	}
	// dn is the normal component of velocity difference. Approximates the residual. // todo xuchenhan@: this prob needs to be scaled by dt
	btVector3 impulse = m_contact->m_c0 * (vr +
										   m_total_normal_dv * infoGlobal.m_deformable_cfm +
										   ((m_penetration < 0) ? m_penetration / infoGlobal.m_timeStep * cti.m_normal : btVector3(0, 0, 0)));

	if (!infoGlobal.m_splitImpulse && m_penetration < 0)
	{
		impulse += m_contact->m_c0 * (m_penetration * infoGlobal.m_deformable_erp / infoGlobal.m_timeStep * cti.m_normal);
	}
	btVector3 impulse_normal = m_contact->m_c0 * (cti.m_normal * dn);
	btVector3 impulse_tangent = impulse - impulse_normal;
	if (dn > 0)
	{
		return 0;
	}
	m_binding = true;
	btScalar residualSquare = dn * dn;
	btVector3 old_total_tangent_dv = m_total_tangent_dv;
	// m_c5 is the inverse mass of the deformable node/face
	m_total_normal_dv -= m_contact->m_c5 * impulse_normal;
	m_total_tangent_dv -= m_contact->m_c5 * impulse_tangent;

	if (m_total_normal_dv.dot(cti.m_normal) < 0)
	{
		// separating in the normal direction
		m_binding = false;
		m_static = false;
		impulse_tangent.setZero();
	}
	else
	{
		if (m_total_normal_dv.norm() * m_contact->m_c3 < m_total_tangent_dv.norm())
		{
			// dynamic friction
			// with dynamic friction, the impulse are still applied to the two objects colliding, however, it does not pose a constraint in the cg solve, hence the change to dv merely serves to update velocity in the contact iterations.
			m_static = false;
			if (m_total_tangent_dv.safeNorm() < SIMD_EPSILON)
			{
				m_total_tangent_dv = btVector3(0, 0, 0);
			}
			else
			{
				m_total_tangent_dv = m_total_tangent_dv.normalized() * m_total_normal_dv.safeNorm() * m_contact->m_c3;
			}
			//            impulse_tangent = -btScalar(1)/m_contact->m_c2 * (m_total_tangent_dv - old_total_tangent_dv);
			impulse_tangent = m_contact->m_c5.inverse() * (old_total_tangent_dv - m_total_tangent_dv);
		}
		else
		{
			// static friction
			m_static = true;
		}
	}
	impulse = impulse_normal + impulse_tangent;
	// apply impulse to deformable nodes involved and change their velocities
	//printf("impulse %f %f %f\n", impulse.x(), impulse.y(), impulse.z());
	applyImpulse(impulse);

	if (cti.m_contact_point_impulse_magnitude)
		*cti.m_contact_point_impulse_magnitude = impulse.length();

	// apply impulse to the rigid/multibodies involved and change their velocities
	if (cti.m_colObj->getInternalType() == btCollisionObject::CO_RIGID_BODY)
	{
		btRigidBody* rigidCol = (btRigidBody*)btRigidBody::upcast(cti.m_colObj);
		if (rigidCol)
		{
			rigidCol->applyImpulse(impulse, m_contact->m_c1);
		}
	}
	else if (cti.m_colObj->getInternalType() == btCollisionObject::CO_FEATHERSTONE_LINK)
	{
		btMultiBodyLinkCollider* multibodyLinkCol = (btMultiBodyLinkCollider*)btMultiBodyLinkCollider::upcast(cti.m_colObj);
		if (multibodyLinkCol)
		{
			const btScalar* deltaV_normal = &m_contact->jacobianData_normal.m_deltaVelocitiesUnitImpulse[0];
			// apply normal component of the impulse
			multibodyLinkCol->m_multiBody->applyDeltaVeeMultiDof2(deltaV_normal, impulse.dot(cti.m_normal));
			if (impulse_tangent.norm() > SIMD_EPSILON)
			{
				// apply tangential component of the impulse
				const btScalar* deltaV_t1 = &m_contact->jacobianData_t1.m_deltaVelocitiesUnitImpulse[0];
				multibodyLinkCol->m_multiBody->applyDeltaVeeMultiDof2(deltaV_t1, impulse.dot(m_contact->t1));
				const btScalar* deltaV_t2 = &m_contact->jacobianData_t2.m_deltaVelocitiesUnitImpulse[0];
				multibodyLinkCol->m_multiBody->applyDeltaVeeMultiDof2(deltaV_t2, impulse.dot(m_contact->t2));
			}
		}
	}
	return residualSquare;
}

btScalar btDeformableRigidContactConstraint::solveSplitImpulse(const btContactSolverInfo& infoGlobal)
{
	btScalar MAX_PENETRATION_CORRECTION = infoGlobal.m_deformable_maxErrorReduction;
	const btSoftBody::sCti& cti = m_contact->m_cti;
	btVector3 vb = getSplitVb();
	btVector3 va = getSplitVa();
	btScalar p = m_penetration;
	if (p > 0)
	{
		return 0;
	}
	btVector3 vr = vb - va;
	btScalar dn = btDot(vr, cti.m_normal) + p * infoGlobal.m_deformable_erp / infoGlobal.m_timeStep;
	if (dn > 0)
	{
		return 0;
	}
	if (m_total_split_impulse + dn > MAX_PENETRATION_CORRECTION)
	{
		dn = MAX_PENETRATION_CORRECTION - m_total_split_impulse;
	}
	if (m_total_split_impulse + dn < -MAX_PENETRATION_CORRECTION)
	{
		dn = -MAX_PENETRATION_CORRECTION - m_total_split_impulse;
	}
	m_total_split_impulse += dn;

	btScalar residualSquare = dn * dn;
	const btVector3 impulse = m_contact->m_c0 * (cti.m_normal * dn);
	applySplitImpulse(impulse);

	// apply split impulse to the rigid/multibodies involved and change their velocities
	if (cti.m_colObj->getInternalType() == btCollisionObject::CO_RIGID_BODY)
	{
		btRigidBody* rigidCol = 0;
		rigidCol = (btRigidBody*)btRigidBody::upcast(cti.m_colObj);
		if (rigidCol)
		{
			rigidCol->applyPushImpulse(impulse, m_contact->m_c1);
		}
	}
	else if (cti.m_colObj->getInternalType() == btCollisionObject::CO_FEATHERSTONE_LINK)
	{
		btMultiBodyLinkCollider* multibodyLinkCol = 0;
		multibodyLinkCol = (btMultiBodyLinkCollider*)btMultiBodyLinkCollider::upcast(cti.m_colObj);
		if (multibodyLinkCol)
		{
			const btScalar* deltaV_normal = &m_contact->jacobianData_normal.m_deltaVelocitiesUnitImpulse[0];
			// apply normal component of the impulse
			multibodyLinkCol->m_multiBody->applyDeltaSplitVeeMultiDof(deltaV_normal, impulse.dot(cti.m_normal));
		}
	}
	return residualSquare;
}
/* ================   Node vs. Rigid   =================== */
btDeformableNodeRigidContactConstraint::btDeformableNodeRigidContactConstraint(const btSoftBody::DeformableNodeRigidContact& contact, const btContactSolverInfo& infoGlobal)
	: m_node(contact.m_node), btDeformableRigidContactConstraint(contact, infoGlobal)
{
}

btDeformableNodeRigidContactConstraint::btDeformableNodeRigidContactConstraint(const btDeformableNodeRigidContactConstraint& other)
	: m_node(other.m_node), btDeformableRigidContactConstraint(other)
{
}

btVector3 btDeformableNodeRigidContactConstraint::getVb() const
{
	return m_node->m_v;
}

btVector3 btDeformableNodeRigidContactConstraint::getSplitVb() const
{
	return m_node->m_splitv;
}

btVector3 btDeformableNodeRigidContactConstraint::getDv(const btSoftBody::Node* node) const
{
	return m_total_normal_dv + m_total_tangent_dv;
}

void btDeformableNodeRigidContactConstraint::applyImpulse(const btVector3& impulse)
{
	const btSoftBody::DeformableNodeRigidContact* contact = getContact();
	btVector3 dv = contact->m_c5 * impulse;
	contact->m_node->m_v -= dv;
}

void btDeformableNodeRigidContactConstraint::applySplitImpulse(const btVector3& impulse)
{
	const btSoftBody::DeformableNodeRigidContact* contact = getContact();
	btVector3 dv = contact->m_c5 * impulse;
	contact->m_node->m_splitv -= dv;
}

/* ================   Face vs. Rigid   =================== */
btDeformableFaceRigidContactConstraint::btDeformableFaceRigidContactConstraint(const btSoftBody::DeformableFaceRigidContact& contact, const btContactSolverInfo& infoGlobal, bool useStrainLimiting)
	: m_face(contact.m_face), m_useStrainLimiting(useStrainLimiting), btDeformableRigidContactConstraint(contact, infoGlobal)
{
}

btDeformableFaceRigidContactConstraint::btDeformableFaceRigidContactConstraint(const btDeformableFaceRigidContactConstraint& other)
	: m_face(other.m_face), m_useStrainLimiting(other.m_useStrainLimiting), btDeformableRigidContactConstraint(other)
{
}

btVector3 btDeformableFaceRigidContactConstraint::getVb() const
{
	const btSoftBody::DeformableFaceRigidContact* contact = getContact();
	btVector3 vb = m_face->m_n[0]->m_v * contact->m_bary[0] + m_face->m_n[1]->m_v * contact->m_bary[1] + m_face->m_n[2]->m_v * contact->m_bary[2];
	return vb;
}

btVector3 btDeformableFaceRigidContactConstraint::getDv(const btSoftBody::Node* node) const
{
	btVector3 face_dv = m_total_normal_dv + m_total_tangent_dv;
	const btSoftBody::DeformableFaceRigidContact* contact = getContact();
	if (m_face->m_n[0] == node)
	{
		return face_dv * contact->m_weights[0];
	}
	if (m_face->m_n[1] == node)
	{
		return face_dv * contact->m_weights[1];
	}
	btAssert(node == m_face->m_n[2]);
	return face_dv * contact->m_weights[2];
}

void btDeformableFaceRigidContactConstraint::applyImpulse(const btVector3& impulse)
{
	const btSoftBody::DeformableFaceRigidContact* contact = getContact();
	btVector3 dv = impulse * contact->m_c2;
	btSoftBody::Face* face = contact->m_face;

	// save applied impulse
	contact->m_cti.m_impulse = impulse;

	btVector3& v0 = face->m_n[0]->m_v;
	btVector3& v1 = face->m_n[1]->m_v;
	btVector3& v2 = face->m_n[2]->m_v;
	const btScalar& im0 = face->m_n[0]->m_im;
	const btScalar& im1 = face->m_n[1]->m_im;
	const btScalar& im2 = face->m_n[2]->m_im;
	if (im0 > 0)
		v0 -= dv * contact->m_weights[0];
	if (im1 > 0)
		v1 -= dv * contact->m_weights[1];
	if (im2 > 0)
		v2 -= dv * contact->m_weights[2];
	if (m_useStrainLimiting)
	{
		btScalar relaxation = 1. / btScalar(m_infoGlobal->m_numIterations);
		btScalar m01 = (relaxation / (im0 + im1));
		btScalar m02 = (relaxation / (im0 + im2));
		btScalar m12 = (relaxation / (im1 + im2));
#ifdef USE_STRAIN_RATE_LIMITING
		// apply strain limiting to prevent the new velocity to change the current length of the edge by more than 1%.
		btScalar p = 0.01;
		btVector3& x0 = face->m_n[0]->m_x;
		btVector3& x1 = face->m_n[1]->m_x;
		btVector3& x2 = face->m_n[2]->m_x;
		const btVector3 x_diff[3] = {x1 - x0, x2 - x0, x2 - x1};
		const btVector3 v_diff[3] = {v1 - v0, v2 - v0, v2 - v1};
		btVector3 u[3];
		btScalar x_diff_dot_u, dn[3];
		btScalar dt = m_infoGlobal->m_timeStep;
		for (int i = 0; i < 3; ++i)
		{
			btScalar x_diff_norm = x_diff[i].safeNorm();
			btScalar x_diff_norm_new = (x_diff[i] + v_diff[i] * dt).safeNorm();
			btScalar strainRate = x_diff_norm_new / x_diff_norm;
			u[i] = v_diff[i];
			u[i].safeNormalize();
			if (x_diff_norm == 0 || (1 - p <= strainRate && strainRate <= 1 + p))
			{
				dn[i] = 0;
				continue;
			}
			x_diff_dot_u = btDot(x_diff[i], u[i]);
			btScalar s;
			if (1 - p > strainRate)
			{
				s = 1 / dt * (-x_diff_dot_u - btSqrt(x_diff_dot_u * x_diff_dot_u + (p * p - 2 * p) * x_diff_norm * x_diff_norm));
			}
			else
			{
				s = 1 / dt * (-x_diff_dot_u + btSqrt(x_diff_dot_u * x_diff_dot_u + (p * p + 2 * p) * x_diff_norm * x_diff_norm));
			}
			//		x_diff_norm_new = (x_diff[i] + s * u[i] * dt).safeNorm();
			//		strainRate = x_diff_norm_new/x_diff_norm;
			dn[i] = s - v_diff[i].safeNorm();
		}
		btVector3 dv0 = im0 * (m01 * u[0] * (-dn[0]) + m02 * u[1] * -(dn[1]));
		btVector3 dv1 = im1 * (m01 * u[0] * (dn[0]) + m12 * u[2] * (-dn[2]));
		btVector3 dv2 = im2 * (m12 * u[2] * (dn[2]) + m02 * u[1] * (dn[1]));
#else
		// apply strain limiting to prevent undamped modes
		btVector3 dv0 = im0 * (m01 * (v1 - v0) + m02 * (v2 - v0));
		btVector3 dv1 = im1 * (m01 * (v0 - v1) + m12 * (v2 - v1));
		btVector3 dv2 = im2 * (m12 * (v1 - v2) + m02 * (v0 - v2));
#endif
		v0 += dv0;
		v1 += dv1;
		v2 += dv2;
	}
}

btVector3 btDeformableFaceRigidContactConstraint::getSplitVb() const
{
	const btSoftBody::DeformableFaceRigidContact* contact = getContact();
	btVector3 vb = (m_face->m_n[0]->m_splitv) * contact->m_bary[0] + (m_face->m_n[1]->m_splitv) * contact->m_bary[1] + (m_face->m_n[2]->m_splitv) * contact->m_bary[2];
	return vb;
}

void btDeformableFaceRigidContactConstraint::applySplitImpulse(const btVector3& impulse)
{
	const btSoftBody::DeformableFaceRigidContact* contact = getContact();
	btVector3 dv = impulse * contact->m_c2;
	btSoftBody::Face* face = contact->m_face;
	btVector3& v0 = face->m_n[0]->m_splitv;
	btVector3& v1 = face->m_n[1]->m_splitv;
	btVector3& v2 = face->m_n[2]->m_splitv;
	const btScalar& im0 = face->m_n[0]->m_im;
	const btScalar& im1 = face->m_n[1]->m_im;
	const btScalar& im2 = face->m_n[2]->m_im;
	if (im0 > 0)
	{
		v0 -= dv * contact->m_weights[0];
	}
	if (im1 > 0)
	{
		v1 -= dv * contact->m_weights[1];
	}
	if (im2 > 0)
	{
		v2 -= dv * contact->m_weights[2];
	}
}

/* ================   Face vs. Node   =================== */
btDeformableFaceNodeContactConstraint::btDeformableFaceNodeContactConstraint(const btSoftBody::DeformableFaceNodeContact& contact, const btContactSolverInfo& infoGlobal)
	: m_node(contact.m_node), m_face(contact.m_face), m_contact(&contact), btDeformableContactConstraint(contact.m_normal, infoGlobal)
{
	m_total_normal_dv.setZero();
	m_total_tangent_dv.setZero();
}

btVector3 btDeformableFaceNodeContactConstraint::getVa() const
{
	return m_node->m_v;
}

btVector3 btDeformableFaceNodeContactConstraint::getVb() const
{
	const btSoftBody::DeformableFaceNodeContact* contact = getContact();
	btVector3 vb = m_face->m_n[0]->m_v * contact->m_bary[0] + m_face->m_n[1]->m_v * contact->m_bary[1] + m_face->m_n[2]->m_v * contact->m_bary[2];
	return vb;
}

btVector3 btDeformableFaceNodeContactConstraint::getDv(const btSoftBody::Node* n) const
{
	btVector3 dv = m_total_normal_dv + m_total_tangent_dv;
	if (n == m_node)
		return dv;
	const btSoftBody::DeformableFaceNodeContact* contact = getContact();
	if (m_face->m_n[0] == n)
	{
		return dv * contact->m_weights[0];
	}
	if (m_face->m_n[1] == n)
	{
		return dv * contact->m_weights[1];
	}
	btAssert(n == m_face->m_n[2]);
	return dv * contact->m_weights[2];
}

btScalar btDeformableFaceNodeContactConstraint::solveConstraint(const btContactSolverInfo& infoGlobal)
{
	btVector3 va = getVa();
	btVector3 vb = getVb();
	btVector3 vr = vb - va;
	const btScalar dn = btDot(vr, m_contact->m_normal);
	// dn is the normal component of velocity diffrerence. Approximates the residual. // todo xuchenhan@: this prob needs to be scaled by dt
	btScalar residualSquare = dn * dn;
	btVector3 impulse = m_contact->m_c0 * vr;
	const btVector3 impulse_normal = m_contact->m_c0 * (m_contact->m_normal * dn);
	btVector3 impulse_tangent = impulse - impulse_normal;

	btVector3 old_total_tangent_dv = m_total_tangent_dv;
	// m_c2 is the inverse mass of the deformable node/face
	if (m_node->m_im > 0)
	{
		m_total_normal_dv -= impulse_normal * m_node->m_im;
		m_total_tangent_dv -= impulse_tangent * m_node->m_im;
	}
	else
	{
		m_total_normal_dv -= impulse_normal * m_contact->m_imf;
		m_total_tangent_dv -= impulse_tangent * m_contact->m_imf;
	}

	if (m_total_normal_dv.dot(m_contact->m_normal) > 0)
	{
		// separating in the normal direction
		m_static = false;
		m_total_tangent_dv = btVector3(0, 0, 0);
		impulse_tangent.setZero();
	}
	else
	{
		if (m_total_normal_dv.norm() * m_contact->m_friction < m_total_tangent_dv.norm())
		{
			// dynamic friction
			// with dynamic friction, the impulse are still applied to the two objects colliding, however, it does not pose a constraint in the cg solve, hence the change to dv merely serves to update velocity in the contact iterations.
			m_static = false;
			if (m_total_tangent_dv.safeNorm() < SIMD_EPSILON)
			{
				m_total_tangent_dv = btVector3(0, 0, 0);
			}
			else
			{
				m_total_tangent_dv = m_total_tangent_dv.normalized() * m_total_normal_dv.safeNorm() * m_contact->m_friction;
			}
			impulse_tangent = -btScalar(1) / m_node->m_im * (m_total_tangent_dv - old_total_tangent_dv);
		}
		else
		{
			// static friction
			m_static = true;
		}
	}
	impulse = impulse_normal + impulse_tangent;
	// apply impulse to deformable nodes involved and change their velocities
	applyImpulse(impulse);
	return residualSquare;
}

void btDeformableFaceNodeContactConstraint::applyImpulse(const btVector3& impulse)
{
	const btSoftBody::DeformableFaceNodeContact* contact = getContact();
	btVector3 dva = impulse * contact->m_node->m_im;
	btVector3 dvb = impulse * contact->m_imf;
	if (contact->m_node->m_im > 0)
	{
		contact->m_node->m_v += dva;
	}

	btSoftBody::Face* face = contact->m_face;
	btVector3& v0 = face->m_n[0]->m_v;
	btVector3& v1 = face->m_n[1]->m_v;
	btVector3& v2 = face->m_n[2]->m_v;
	const btScalar& im0 = face->m_n[0]->m_im;
	const btScalar& im1 = face->m_n[1]->m_im;
	const btScalar& im2 = face->m_n[2]->m_im;
	if (im0 > 0)
	{
		v0 -= dvb * contact->m_weights[0];
	}
	if (im1 > 0)
	{
		v1 -= dvb * contact->m_weights[1];
	}
	if (im2 > 0)
	{
		v2 -= dvb * contact->m_weights[2];
	}
}
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

#ifndef BT_DEFORMABLE_CONTACT_CONSTRAINT_H
#define BT_DEFORMABLE_CONTACT_CONSTRAINT_H
#include "btSoftBody.h"

// btDeformableContactConstraint is an abstract class specifying the method that each type of contact constraint needs to implement
class btDeformableContactConstraint
{
public:
	// True if the friction is static
	// False if the friction is dynamic
	bool m_static;
	const btContactSolverInfo* m_infoGlobal;

	// normal of the contact
	btVector3 m_normal;

	btDeformableContactConstraint(const btVector3& normal, const btContactSolverInfo& infoGlobal) : m_static(false), m_normal(normal), m_infoGlobal(&infoGlobal)
	{
	}

	btDeformableContactConstraint(bool isStatic, const btVector3& normal, const btContactSolverInfo& infoGlobal) : m_static(isStatic), m_normal(normal), m_infoGlobal(&infoGlobal)
	{
	}

	btDeformableContactConstraint() : m_static(false) {}

	btDeformableContactConstraint(const btDeformableContactConstraint& other)
		: m_static(other.m_static), m_normal(other.m_normal), m_infoGlobal(other.m_infoGlobal)
	{
	}

	virtual ~btDeformableContactConstraint() {}

	// solve the constraint with inelastic impulse and return the error, which is the square of normal component of velocity diffrerence
	// the constraint is solved by calculating the impulse between object A and B in the contact and apply the impulse to both objects involved in the contact
	virtual btScalar solveConstraint(const btContactSolverInfo& infoGlobal) = 0;

	// get the velocity of the object A in the contact
	virtual btVector3 getVa() const = 0;

	// get the velocity of the object B in the contact
	virtual btVector3 getVb() const = 0;

	// get the velocity change of the soft body node in the constraint
	virtual btVector3 getDv(const btSoftBody::Node*) const = 0;

	// apply impulse to the soft body node and/or face involved
	virtual void applyImpulse(const btVector3& impulse) = 0;

	// scale the penetration depth by erp
	virtual void setPenetrationScale(btScalar scale) = 0;
};

//
// Constraint that a certain node in the deformable objects cannot move
class btDeformableStaticConstraint : public btDeformableContactConstraint
{
public:
	btSoftBody::Node* m_node;

	btDeformableStaticConstraint(btSoftBody::Node* node, const btContactSolverInfo& infoGlobal) : m_node(node), btDeformableContactConstraint(false, btVector3(0, 0, 0), infoGlobal)
	{
	}
	btDeformableStaticConstraint() {}
	btDeformableStaticConstraint(const btDeformableStaticConstraint& other)
		: m_node(other.m_node), btDeformableContactConstraint(other)
	{
	}

	virtual ~btDeformableStaticConstraint() {}

	virtual btScalar solveConstraint(const btContactSolverInfo& infoGlobal)
	{
		return 0;
	}

	virtual btVector3 getVa() const
	{
		return btVector3(0, 0, 0);
	}

	virtual btVector3 getVb() const
	{
		return btVector3(0, 0, 0);
	}

	virtual btVector3 getDv(const btSoftBody::Node* n) const
	{
		return btVector3(0, 0, 0);
	}

	virtual void applyImpulse(const btVector3& impulse) {}
	virtual void setPenetrationScale(btScalar scale) {}
};

//
// Anchor Constraint between rigid and deformable node
class btDeformableNodeAnchorConstraint : public btDeformableContactConstraint
{
public:
	const btSoftBody::DeformableNodeRigidAnchor* m_anchor;

	btDeformableNodeAnchorConstraint(const btSoftBody::DeformableNodeRigidAnchor& c, const btContactSolverInfo& infoGlobal);
	btDeformableNodeAnchorConstraint(const btDeformableNodeAnchorConstraint& other);
	btDeformableNodeAnchorConstraint() {}
	virtual ~btDeformableNodeAnchorConstraint()
	{
	}
	virtual btScalar solveConstraint(const btContactSolverInfo& infoGlobal);

	// object A is the rigid/multi body, and object B is the deformable node/face
	virtual btVector3 getVa() const;
	// get the velocity of the deformable node in contact
	virtual btVector3 getVb() const;

	// --- hard positional (split-impulse) anchor solve ---
	btVector3 getSplitVa() const;
	btVector3 getSplitVb() const;
	void applySplitImpulse(const btVector3& impulse);
	btScalar solveSplitImpulse(const btContactSolverInfo& infoGlobal);

	virtual btVector3 getDv(const btSoftBody::Node* n) const
	{
		return btVector3(0, 0, 0);
	}
	virtual void applyImpulse(const btVector3& impulse);

	virtual void setPenetrationScale(btScalar scale) {}
};

//
// Constraint between rigid/multi body and deformable objects
class btDeformableRigidContactConstraint : public btDeformableContactConstraint
{
public:
	btVector3 m_total_normal_dv;
	btVector3 m_total_tangent_dv;
	btScalar m_penetration;
	btScalar m_total_split_impulse;
	bool m_binding;
	const btSoftBody::DeformableRigidContact* m_contact;

	btDeformableRigidContactConstraint(const btSoftBody::DeformableRigidContact& c, const btContactSolverInfo& infoGlobal);
	btDeformableRigidContactConstraint(const btDeformableRigidContactConstraint& other);
	btDeformableRigidContactConstraint() : m_binding(false) {}
	virtual ~btDeformableRigidContactConstraint()
	{
	}

	// object A is the rigid/multi body, and object B is the deformable node/face
	virtual btVector3 getVa() const;

	// get the split impulse velocity of the deformable face at the contact point
	virtual btVector3 getSplitVb() const = 0;

	// get the split impulse velocity of the rigid/multibdoy at the contaft
	virtual btVector3 getSplitVa() const;

	virtual btScalar solveConstraint(const btContactSolverInfo& infoGlobal);

	virtual void setPenetrationScale(btScalar scale)
	{
		m_penetration *= scale;
	}

	btScalar solveSplitImpulse(const btContactSolverInfo& infoGlobal);

	virtual void applySplitImpulse(const btVector3& impulse) = 0;
};

//
// Constraint between rigid/multi body and deformable objects nodes
class btDeformableNodeRigidContactConstraint : public btDeformableRigidContactConstraint
{
public:
	// the deformable node in contact
	btSoftBody::Node* m_node;

	btDeformableNodeRigidContactConstraint(const btSoftBody::DeformableNodeRigidContact& contact, const btContactSolverInfo& infoGlobal);
	btDeformableNodeRigidContactConstraint(const btDeformableNodeRigidContactConstraint& other);
	btDeformableNodeRigidContactConstraint() {}
	virtual ~btDeformableNodeRigidContactConstraint()
	{
	}

	// get the velocity of the deformable node in contact
	virtual btVector3 getVb() const;

	// get the split impulse velocity of the deformable face at the contact point
	virtual btVector3 getSplitVb() const;

	// get the velocity change of the input soft body node in the constraint
	virtual btVector3 getDv(const btSoftBody::Node*) const;

	// cast the contact to the desired type
	const btSoftBody::DeformableNodeRigidContact* getContact() const
	{
		return static_cast<const btSoftBody::DeformableNodeRigidContact*>(m_contact);
	}

	virtual void applyImpulse(const btVector3& impulse);

	virtual void applySplitImpulse(const btVector3& impulse);
};

//
// Constraint between rigid/multi body and deformable objects faces
class btDeformableFaceRigidContactConstraint : public btDeformableRigidContactConstraint
{
public:
	btSoftBody::Face* m_face;
	bool m_useStrainLimiting;
	btDeformableFaceRigidContactConstraint(const btSoftBody::DeformableFaceRigidContact& contact, const btContactSolverInfo& infoGlobal, bool useStrainLimiting);
	btDeformableFaceRigidContactConstraint(const btDeformableFaceRigidContactConstraint& other);
	btDeformableFaceRigidContactConstraint() : m_useStrainLimiting(false) {}
	virtual ~btDeformableFaceRigidContactConstraint()
	{
	}

	// get the velocity of the deformable face at the contact point
	virtual btVector3 getVb() const;

	// get the split impulse velocity of the deformable face at the contact point
	virtual btVector3 getSplitVb() const;

	// get the velocity change of the input soft body node in the constraint
	virtual btVector3 getDv(const btSoftBody::Node*) const;

	// cast the contact to the desired type
	const btSoftBody::DeformableFaceRigidContact* getContact() const
	{
		return static_cast<const btSoftBody::DeformableFaceRigidContact*>(m_contact);
	}

	virtual void applyImpulse(const btVector3& impulse);

	virtual void applySplitImpulse(const btVector3& impulse);
};

//
// Constraint between  deformable objects faces and deformable objects nodes
class btDeformableFaceNodeContactConstraint : public btDeformableContactConstraint
{
public:
	btSoftBody::Node* m_node;
	btSoftBody::Face* m_face;
	const btSoftBody::DeformableFaceNodeContact* m_contact;
	btVector3 m_total_normal_dv;
	btVector3 m_total_tangent_dv;

	btDeformableFaceNodeContactConstraint(const btSoftBody::DeformableFaceNodeContact& contact, const btContactSolverInfo& infoGlobal);
	btDeformableFaceNodeContactConstraint() {}
	virtual ~btDeformableFaceNodeContactConstraint() {}

	virtual btScalar solveConstraint(const btContactSolverInfo& infoGlobal);

	// get the velocity of the object A in the contact
	virtual btVector3 getVa() const;

	// get the velocity of the object B in the contact
	virtual btVector3 getVb() const;

	// get the velocity change of the input soft body node in the constraint
	virtual btVector3 getDv(const btSoftBody::Node*) const;

	// cast the contact to the desired type
	const btSoftBody::DeformableFaceNodeContact* getContact() const
	{
		return static_cast<const btSoftBody::DeformableFaceNodeContact*>(m_contact);
	}

	virtual void applyImpulse(const btVector3& impulse);

	virtual void setPenetrationScale(btScalar scale) {}
};
#endif /* BT_DEFORMABLE_CONTACT_CONSTRAINT_H */
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

#include "btDeformableContactProjection.h"
#include "btDeformableMultiBodyDynamicsWorld.h"
#include <algorithm>
#include <cmath>
btScalar btDeformableContactProjection::update(btCollisionObject** deformableBodies, int numDeformableBodies, const btContactSolverInfo& infoGlobal)
{
	btScalar residualSquare = 0;
	for (int i = 0; i < numDeformableBodies; ++i)
	{
		for (int j = 0; j < m_softBodies.size(); ++j)
		{
			btCollisionObject* psb = m_softBodies[j];
			if (psb != deformableBodies[i])
			{
				continue;
			}
			for (int k = 0; k < m_nodeRigidConstraints[j].size(); ++k)
			{
				btDeformableNodeRigidContactConstraint& constraint = m_nodeRigidConstraints[j][k];
				btScalar localResidualSquare = constraint.solveConstraint(infoGlobal);
				residualSquare = btMax(residualSquare, localResidualSquare);
			}
			for (int k = 0; k < m_nodeAnchorConstraints[j].size(); ++k)
			{
				btDeformableNodeAnchorConstraint& constraint = m_nodeAnchorConstraints[j][k];
				btScalar localResidualSquare = constraint.solveConstraint(infoGlobal);
				residualSquare = btMax(residualSquare, localResidualSquare);
			}
			for (int k = 0; k < m_faceRigidConstraints[j].size(); ++k)
			{
				btDeformableFaceRigidContactConstraint& constraint = m_faceRigidConstraints[j][k];
				btScalar localResidualSquare = constraint.solveConstraint(infoGlobal);
				residualSquare = btMax(residualSquare, localResidualSquare);
			}
			for (int k = 0; k < m_deformableConstraints[j].size(); ++k)
			{
				btDeformableFaceNodeContactConstraint& constraint = m_deformableConstraints[j][k];
				btScalar localResidualSquare = constraint.solveConstraint(infoGlobal);
				residualSquare = btMax(residualSquare, localResidualSquare);
			}
		}
	}
	return residualSquare;
}

btScalar btDeformableContactProjection::solveSplitImpulse(btCollisionObject** deformableBodies, int numDeformableBodies, const btContactSolverInfo& infoGlobal)
{
	btScalar residualSquare = 0;
	for (int i = 0; i < numDeformableBodies; ++i)
	{
		for (int j = 0; j < m_softBodies.size(); ++j)
		{
			btCollisionObject* psb = m_softBodies[j];
			if (psb != deformableBodies[i])
			{
				continue;
			}
			for (int k = 0; k < m_nodeRigidConstraints[j].size(); ++k)
			{
				btDeformableNodeRigidContactConstraint& constraint = m_nodeRigidConstraints[j][k];
				btScalar localResidualSquare = constraint.solveSplitImpulse(infoGlobal);
				residualSquare = btMax(residualSquare, localResidualSquare);
			}
			for (int k = 0; k < m_faceRigidConstraints[j].size(); ++k)
			{
				btDeformableFaceRigidContactConstraint& constraint = m_faceRigidConstraints[j][k];
				btScalar localResidualSquare = constraint.solveSplitImpulse(infoGlobal);
				residualSquare = btMax(residualSquare, localResidualSquare);
			}
			for (int k = 0; k < m_nodeAnchorConstraints[j].size(); ++k)
			{
				btDeformableNodeAnchorConstraint& constraint = m_nodeAnchorConstraints[j][k];
				btScalar localResidualSquare = constraint.solveSplitImpulse(infoGlobal);
				residualSquare = btMax(residualSquare, localResidualSquare);
			}
		}
	}
	return residualSquare;
}

void btDeformableContactProjection::setConstraints(const btContactSolverInfo& infoGlobal)
{
	BT_PROFILE("setConstraints");
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		if (!psb->isActive() || psb->isStaticObject())
		{
			continue;
		}

		// set Dirichlet constraint
		for (int j = 0; j < psb->m_nodes.size(); ++j)
		{
			if (psb->m_nodes[j].m_im == 0)
			{
				btDeformableStaticConstraint static_constraint(&psb->m_nodes[j], infoGlobal);
				m_staticConstraints[i].push_back(static_constraint);
			}
		}

		// set up deformable anchors
		for (int j = 0; j < psb->m_deformableAnchors.size(); ++j)
		{
			btSoftBody::DeformableNodeRigidAnchor& anchor = psb->m_deformableAnchors[j];
			// skip fixed points
			if (anchor.m_node->m_im == 0)
			{
				continue;
			}
			anchor.m_c1 = anchor.m_cti.m_colObj->getWorldTransform().getBasis() * anchor.m_local;
			btDeformableNodeAnchorConstraint constraint(anchor, infoGlobal);
			m_nodeAnchorConstraints[i].push_back(constraint);
		}

		// set Deformable Node vs. Rigid constraint
		//fprintf(stderr, "framestart()\n");
		/*if (psb->m_nodeRigidContacts.size() > 0)
			for (auto fi = 0; fi < psb->m_faces.size(); ++fi)
			{
				auto& f = psb->m_faces[fi];
				fprintf(stderr, "drawtriangle \"tri\" [%f,%f,%f][%f,%f,%f][%f,%f,%f]\n", f.m_n[0]->m_x.x(), f.m_n[0]->m_x.y(), f.m_n[0]->m_x.z(),
						f.m_n[1]->m_x.x(), f.m_n[1]->m_x.y(), f.m_n[1]->m_x.z(),
						f.m_n[2]->m_x.x(), f.m_n[2]->m_x.y(), f.m_n[2]->m_x.z());
			}*/
		//fprintf(stderr, "psb->m_nodeRigidContacts %d\n", psb->m_nodeRigidContacts.size());
		for (int j = 0; j < psb->m_nodeRigidContacts.size(); ++j)
		{
			const btSoftBody::DeformableNodeRigidContact& contact = psb->m_nodeRigidContacts[j];
			// skip fixed points
			if (contact.m_node->m_im == 0)
			{
				continue;
			}

			auto lineStart = contact.m_node->m_x;
			auto lineEnd = contact.m_node->m_x + contact.m_cti.m_normal * 100.0;
			//fprintf(stderr, "drawline \"line%d\" [%f,%f,%f][%f,%f,%f] \n", contact.m_cti.m_count, lineStart.x(), lineStart.y(), lineStart.z(),
			//		lineEnd.x(), lineEnd.y(), lineEnd.z());
			//fprintf(stderr, "normal %d %f %f %f\n", j, contact.m_cti.m_normal.x(), contact.m_cti.m_normal.y(), contact.m_cti.m_normal.z());
			btDeformableNodeRigidContactConstraint constraint(contact, infoGlobal);
			m_nodeRigidConstraints[i].push_back(constraint);
		}

		// set Deformable Face vs. Rigid constraint
		for (int j = 0; j < psb->m_faceRigidContacts.size(); ++j)
		{
			const btSoftBody::DeformableFaceRigidContact& contact = psb->m_faceRigidContacts[j];
			// skip fixed faces
			if (contact.m_c2 == 0)
			{
				continue;
			}
			btDeformableFaceRigidContactConstraint constraint(contact, infoGlobal, m_useStrainLimiting);
			m_faceRigidConstraints[i].push_back(constraint);
		}
	}
	//fprintf(stderr, "frameend()\n");
}

void btDeformableContactProjection::project(TVStack& x)
{
#ifndef USE_MGS
	const int dim = 3;
	for (int index = 0; index < m_projectionsDict.size(); ++index)
	{
		btAlignedObjectArray<btVector3>& projectionDirs = *m_projectionsDict.getAtIndex(index);
		size_t i = m_projectionsDict.getKeyAtIndex(index).getUid1();
		if (projectionDirs.size() >= dim)
		{
			// static node
			x[i].setZero();
			continue;
		}
		else if (projectionDirs.size() == 2)
		{
			btVector3 dir0 = projectionDirs[0];
			btVector3 dir1 = projectionDirs[1];
			btVector3 free_dir = btCross(dir0, dir1);
			if (free_dir.safeNorm() < SIMD_EPSILON)
			{
				x[i] -= x[i].dot(dir0) * dir0;
			}
			else
			{
				free_dir.normalize();
				x[i] = x[i].dot(free_dir) * free_dir;
			}
		}
		else
		{
			btAssert(projectionDirs.size() == 1);
			btVector3 dir0 = projectionDirs[0];
			x[i] -= x[i].dot(dir0) * dir0;
		}
	}
#else
	btReducedVector p(x.size());
	for (int i = 0; i < m_projections.size(); ++i)
	{
		p += (m_projections[i].dot(x) * m_projections[i]);
	}
	for (int i = 0; i < p.m_indices.size(); ++i)
	{
		x[p.m_indices[i]] -= p.m_vecs[i];
	}
#endif
}

void btDeformableContactProjection::setProjection()
{
#ifndef USE_MGS
	BT_PROFILE("btDeformableContactProjection::setProjection");
	btAlignedObjectArray<btVector3> units;
	units.push_back(btVector3(1, 0, 0));
	units.push_back(btVector3(0, 1, 0));
	units.push_back(btVector3(0, 0, 1));
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		if (!psb->isActive() || psb->isStaticObject())
		{
			continue;
		}
		for (int j = 0; j < m_staticConstraints[i].size(); ++j)
		{
			int index = m_staticConstraints[i][j].m_node->index;
			m_staticConstraints[i][j].m_node->m_constrained = true;
			if (m_projectionsDict.find(index) == NULL)
			{
				m_projectionsDict.insert(index, units);
			}
			else
			{
				btAlignedObjectArray<btVector3>& projections = *m_projectionsDict[index];
				for (int k = 0; k < 3; ++k)
				{
					projections.push_back(units[k]);
				}
			}
		}
		for (int j = 0; j < m_nodeAnchorConstraints[i].size(); ++j)
		{
			int index = m_nodeAnchorConstraints[i][j].m_anchor->m_node->index;
			m_nodeAnchorConstraints[i][j].m_anchor->m_node->m_constrained = true;
			if (m_projectionsDict.find(index) == NULL)
			{
				m_projectionsDict.insert(index, units);
			}
			else
			{
				btAlignedObjectArray<btVector3>& projections = *m_projectionsDict[index];
				for (int k = 0; k < 3; ++k)
				{
					projections.push_back(units[k]);
				}
			}
		}
		for (int j = 0; j < m_nodeRigidConstraints[i].size(); ++j)
		{
			int index = m_nodeRigidConstraints[i][j].m_node->index;
			m_nodeRigidConstraints[i][j].m_node->m_constrained = true;
			if (m_nodeRigidConstraints[i][j].m_binding)
			{
				if (m_nodeRigidConstraints[i][j].m_static)
				{
					if (m_projectionsDict.find(index) == NULL)
					{
						m_projectionsDict.insert(index, units);
					}
					else
					{
						btAlignedObjectArray<btVector3>& projections = *m_projectionsDict[index];
						for (int k = 0; k < 3; ++k)
						{
							projections.push_back(units[k]);
						}
					}
				}
				else
				{
					if (m_projectionsDict.find(index) == NULL)
					{
						btAlignedObjectArray<btVector3> projections;
						projections.push_back(m_nodeRigidConstraints[i][j].m_normal);
						m_projectionsDict.insert(index, projections);
					}
					else
					{
						btAlignedObjectArray<btVector3>& projections = *m_projectionsDict[index];
						projections.push_back(m_nodeRigidConstraints[i][j].m_normal);
					}
				}
			}
		}
		for (int j = 0; j < m_faceRigidConstraints[i].size(); ++j)
		{
			const btSoftBody::Face* face = m_faceRigidConstraints[i][j].m_face;
			if (m_faceRigidConstraints[i][j].m_binding)
			{
				for (int k = 0; k < 3; ++k)
				{
					face->m_n[k]->m_constrained = true;
				}
			}
			for (int k = 0; k < 3; ++k)
			{
				btSoftBody::Node* node = face->m_n[k];
				int index = node->index;
				if (m_faceRigidConstraints[i][j].m_static)
				{
					if (m_projectionsDict.find(index) == NULL)
					{
						m_projectionsDict.insert(index, units);
					}
					else
					{
						btAlignedObjectArray<btVector3>& projections = *m_projectionsDict[index];
						for (int l = 0; l < 3; ++l)
						{
							projections.push_back(units[l]);
						}
					}
				}
				else
				{
					if (m_projectionsDict.find(index) == NULL)
					{
						btAlignedObjectArray<btVector3> projections;
						projections.push_back(m_faceRigidConstraints[i][j].m_normal);
						m_projectionsDict.insert(index, projections);
					}
					else
					{
						btAlignedObjectArray<btVector3>& projections = *m_projectionsDict[index];
						projections.push_back(m_faceRigidConstraints[i][j].m_normal);
					}
				}
			}
		}
	}
#else
	int dof = 0;
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		dof += m_softBodies[i]->m_nodes.size();
	}
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		if (!psb->isActive() || psb->isStaticObject())
		{
			continue;
		}
		for (int j = 0; j < m_staticConstraints[i].size(); ++j)
		{
			int index = m_staticConstraints[i][j].m_node->index;
			m_staticConstraints[i][j].m_node->m_penetration = SIMD_INFINITY;
			btAlignedObjectArray<int> indices;
			btAlignedObjectArray<btVector3> vecs1, vecs2, vecs3;
			indices.push_back(index);
			vecs1.push_back(btVector3(1, 0, 0));
			vecs2.push_back(btVector3(0, 1, 0));
			vecs3.push_back(btVector3(0, 0, 1));
			m_projections.push_back(btReducedVector(dof, indices, vecs1));
			m_projections.push_back(btReducedVector(dof, indices, vecs2));
			m_projections.push_back(btReducedVector(dof, indices, vecs3));
		}

		for (int j = 0; j < m_nodeAnchorConstraints[i].size(); ++j)
		{
			int index = m_nodeAnchorConstraints[i][j].m_anchor->m_node->index;
			m_nodeAnchorConstraints[i][j].m_anchor->m_node->m_penetration = SIMD_INFINITY;
			btAlignedObjectArray<int> indices;
			btAlignedObjectArray<btVector3> vecs1, vecs2, vecs3;
			indices.push_back(index);
			vecs1.push_back(btVector3(1, 0, 0));
			vecs2.push_back(btVector3(0, 1, 0));
			vecs3.push_back(btVector3(0, 0, 1));
			m_projections.push_back(btReducedVector(dof, indices, vecs1));
			m_projections.push_back(btReducedVector(dof, indices, vecs2));
			m_projections.push_back(btReducedVector(dof, indices, vecs3));
		}
		for (int j = 0; j < m_nodeRigidConstraints[i].size(); ++j)
		{
			int index = m_nodeRigidConstraints[i][j].m_node->index;
			m_nodeRigidConstraints[i][j].m_node->m_penetration = -m_nodeRigidConstraints[i][j].getContact()->m_cti.m_offset;
			btAlignedObjectArray<int> indices;
			indices.push_back(index);
			btAlignedObjectArray<btVector3> vecs1, vecs2, vecs3;
			if (m_nodeRigidConstraints[i][j].m_static)
			{
				vecs1.push_back(btVector3(1, 0, 0));
				vecs2.push_back(btVector3(0, 1, 0));
				vecs3.push_back(btVector3(0, 0, 1));
				m_projections.push_back(btReducedVector(dof, indices, vecs1));
				m_projections.push_back(btReducedVector(dof, indices, vecs2));
				m_projections.push_back(btReducedVector(dof, indices, vecs3));
			}
			else
			{
				vecs1.push_back(m_nodeRigidConstraints[i][j].m_normal);
				m_projections.push_back(btReducedVector(dof, indices, vecs1));
			}
		}
		for (int j = 0; j < m_faceRigidConstraints[i].size(); ++j)
		{
			const btSoftBody::Face* face = m_faceRigidConstraints[i][j].m_face;
			btVector3 bary = m_faceRigidConstraints[i][j].getContact()->m_bary;
			btScalar penetration = -m_faceRigidConstraints[i][j].getContact()->m_cti.m_offset;
			for (int k = 0; k < 3; ++k)
			{
				face->m_n[k]->m_penetration = btMax(face->m_n[k]->m_penetration, penetration);
			}
			if (m_faceRigidConstraints[i][j].m_static)
			{
				for (int l = 0; l < 3; ++l)
				{
					btReducedVector rv(dof);
					for (int k = 0; k < 3; ++k)
					{
						rv.m_indices.push_back(face->m_n[k]->index);
						btVector3 v(0, 0, 0);
						v[l] = bary[k];
						rv.m_vecs.push_back(v);
						rv.sort();
					}
					m_projections.push_back(rv);
				}
			}
			else
			{
				btReducedVector rv(dof);
				for (int k = 0; k < 3; ++k)
				{
					rv.m_indices.push_back(face->m_n[k]->index);
					rv.m_vecs.push_back(bary[k] * m_faceRigidConstraints[i][j].m_normal);
					rv.sort();
				}
				m_projections.push_back(rv);
			}
		}
	}
	btModifiedGramSchmidt<btReducedVector> mgs(m_projections);
	mgs.solve();
	m_projections = mgs.m_out;
#endif
}

void btDeformableContactProjection::checkConstraints(const TVStack& x)
{
	for (int i = 0; i < m_lagrangeMultipliers.size(); ++i)
	{
		btVector3 d(0, 0, 0);
		const LagrangeMultiplier& lm = m_lagrangeMultipliers[i];
		for (int j = 0; j < lm.m_num_constraints; ++j)
		{
			for (int k = 0; k < lm.m_num_nodes; ++k)
			{
				d[j] += lm.m_weights[k] * x[lm.m_indices[k]].dot(lm.m_dirs[j]);
			}
		}
		//		printf("d = %f, %f, %f\n", d[0], d[1], d[2]);
		//        printf("val = %f, %f, %f\n", lm.m_vals[0], lm.m_vals[1], lm.m_vals[2]);
	}
}

void btDeformableContactProjection::setLagrangeMultiplier()
{
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		if (!psb->isActive() || psb->isStaticObject())
		{
			continue;
		}
		for (int j = 0; j < m_staticConstraints[i].size(); ++j)
		{
			int index = m_staticConstraints[i][j].m_node->index;
			m_staticConstraints[i][j].m_node->m_constrained = true;
			LagrangeMultiplier lm;
			lm.m_num_nodes = 1;
			lm.m_indices[0] = index;
			lm.m_weights[0] = 1.0;
			lm.m_num_constraints = 3;
			lm.m_dirs[0] = btVector3(1, 0, 0);
			lm.m_dirs[1] = btVector3(0, 1, 0);
			lm.m_dirs[2] = btVector3(0, 0, 1);
			m_lagrangeMultipliers.push_back(lm);
		}
		for (int j = 0; j < m_nodeAnchorConstraints[i].size(); ++j)
		{
			int index = m_nodeAnchorConstraints[i][j].m_anchor->m_node->index;
			m_nodeAnchorConstraints[i][j].m_anchor->m_node->m_constrained = true;
			LagrangeMultiplier lm;
			lm.m_num_nodes = 1;
			lm.m_indices[0] = index;
			lm.m_weights[0] = 1.0;
			lm.m_num_constraints = 3;
			lm.m_dirs[0] = btVector3(1, 0, 0);
			lm.m_dirs[1] = btVector3(0, 1, 0);
			lm.m_dirs[2] = btVector3(0, 0, 1);
			m_lagrangeMultipliers.push_back(lm);
		}

		for (int j = 0; j < m_nodeRigidConstraints[i].size(); ++j)
		{
			if (!m_nodeRigidConstraints[i][j].m_binding)
			{
				continue;
			}
			int index = m_nodeRigidConstraints[i][j].m_node->index;
			m_nodeRigidConstraints[i][j].m_node->m_constrained = true;
			LagrangeMultiplier lm;
			lm.m_num_nodes = 1;
			lm.m_indices[0] = index;
			lm.m_weights[0] = 1.0;
			if (m_nodeRigidConstraints[i][j].m_static)
			{
				lm.m_num_constraints = 3;
				lm.m_dirs[0] = btVector3(1, 0, 0);
				lm.m_dirs[1] = btVector3(0, 1, 0);
				lm.m_dirs[2] = btVector3(0, 0, 1);
			}
			else
			{
				lm.m_num_constraints = 1;
				lm.m_dirs[0] = m_nodeRigidConstraints[i][j].m_normal;
			}
			m_lagrangeMultipliers.push_back(lm);
		}

		for (int j = 0; j < m_faceRigidConstraints[i].size(); ++j)
		{
			if (!m_faceRigidConstraints[i][j].m_binding)
			{
				continue;
			}
			btSoftBody::Face* face = m_faceRigidConstraints[i][j].m_face;

			btVector3 bary = m_faceRigidConstraints[i][j].getContact()->m_bary;
			LagrangeMultiplier lm;
			lm.m_num_nodes = 3;

			for (int k = 0; k < 3; ++k)
			{
				face->m_n[k]->m_constrained = true;
				lm.m_indices[k] = face->m_n[k]->index;
				lm.m_weights[k] = bary[k];
			}
			if (m_faceRigidConstraints[i][j].m_static)
			{
				face->m_pcontact[3] = 1;
				lm.m_num_constraints = 3;
				lm.m_dirs[0] = btVector3(1, 0, 0);
				lm.m_dirs[1] = btVector3(0, 1, 0);
				lm.m_dirs[2] = btVector3(0, 0, 1);
			}
			else
			{
				face->m_pcontact[3] = 0;
				lm.m_num_constraints = 1;
				lm.m_dirs[0] = m_faceRigidConstraints[i][j].m_normal;
			}
			m_lagrangeMultipliers.push_back(lm);
		}
	}
}

//
void btDeformableContactProjection::applyDynamicFriction(TVStack& f)
{
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		for (int j = 0; j < m_nodeRigidConstraints[i].size(); ++j)
		{
			const btDeformableNodeRigidContactConstraint& constraint = m_nodeRigidConstraints[i][j];
			const btSoftBody::Node* node = constraint.m_node;
			if (node->m_im != 0)
			{
				int index = node->index;
				f[index] += constraint.getDv(node) * (1. / node->m_im);
			}
		}
		for (int j = 0; j < m_faceRigidConstraints[i].size(); ++j)
		{
			const btDeformableFaceRigidContactConstraint& constraint = m_faceRigidConstraints[i][j];
			const btSoftBody::Face* face = constraint.getContact()->m_face;
			for (int k = 0; k < 3; ++k)
			{
				const btSoftBody::Node* node = face->m_n[k];
				if (node->m_im != 0)
				{
					int index = node->index;
					f[index] += constraint.getDv(node) * (1. / node->m_im);
				}
			}
		}
		for (int j = 0; j < m_deformableConstraints[i].size(); ++j)
		{
			const btDeformableFaceNodeContactConstraint& constraint = m_deformableConstraints[i][j];
			const btSoftBody::Face* face = constraint.getContact()->m_face;
			const btSoftBody::Node* node = constraint.getContact()->m_node;
			if (node->m_im != 0)
			{
				int index = node->index;
				f[index] += constraint.getDv(node) * (1. / node->m_im);
			}
			for (int k = 0; k < 3; ++k)
			{
				const btSoftBody::Node* node = face->m_n[k];
				if (node->m_im != 0)
				{
					int index = node->index;
					f[index] += constraint.getDv(node) * (1. / node->m_im);
				}
			}
		}
	}
}

void btDeformableContactProjection::reinitialize(bool nodeUpdated)
{
	int N = m_softBodies.size();
	if (nodeUpdated)
	{
		m_staticConstraints.resize(N);
		m_nodeAnchorConstraints.resize(N);
		m_nodeRigidConstraints.resize(N);
		m_faceRigidConstraints.resize(N);
		m_deformableConstraints.resize(N);
	}
	for (int i = 0; i < N; ++i)
	{
		m_staticConstraints[i].clear();
		m_nodeAnchorConstraints[i].clear();
		m_nodeRigidConstraints[i].clear();
		m_faceRigidConstraints[i].clear();
		m_deformableConstraints[i].clear();
	}
#ifndef USE_MGS
	m_projectionsDict.clear();
#else
	m_projections.clear();
#endif
	m_lagrangeMultipliers.clear();
}
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

#ifndef BT_CONTACT_PROJECTION_H
#define BT_CONTACT_PROJECTION_H
#include "btCGProjection.h"
#include "btSoftBody.h"
#include "BulletDynamics/Featherstone/btMultiBodyLinkCollider.h"
#include "BulletDynamics/Featherstone/btMultiBodyConstraint.h"
#include "btDeformableContactConstraint.h"
#include "LinearMath/btHashMap.h"
#include "LinearMath/btReducedVector.h"
#include "LinearMath/btModifiedGramSchmidt.h"
#include <vector>

struct LagrangeMultiplier
{
	int m_num_constraints;  // Number of constraints
	int m_num_nodes;        // Number of nodes in these constraints
	btScalar m_weights[3];  // weights of the nodes involved, same size as m_num_nodes
	btVector3 m_dirs[3];    // Constraint directions, same size of m_num_constraints;
	int m_indices[3];       // indices of the nodes involved, same size as m_num_nodes;
};

class btDeformableContactProjection
{
public:
	typedef btAlignedObjectArray<btVector3> TVStack;
	btAlignedObjectArray<btSoftBody*>& m_softBodies;

	// all constraints involving face
	btAlignedObjectArray<btDeformableContactConstraint*> m_allFaceConstraints;
#ifndef USE_MGS
	// map from node index to projection directions
	btHashMap<btHashInt, btAlignedObjectArray<btVector3> > m_projectionsDict;
#else
	btAlignedObjectArray<btReducedVector> m_projections;
#endif

	btAlignedObjectArray<LagrangeMultiplier> m_lagrangeMultipliers;

	// map from node index to static constraint
	btAlignedObjectArray<btAlignedObjectArray<btDeformableStaticConstraint> > m_staticConstraints;
	// map from node index to node rigid constraint
	btAlignedObjectArray<btAlignedObjectArray<btDeformableNodeRigidContactConstraint> > m_nodeRigidConstraints;
	// map from node index to face rigid constraint
	btAlignedObjectArray<btAlignedObjectArray<btDeformableFaceRigidContactConstraint> > m_faceRigidConstraints;
	// map from node index to deformable constraint
	btAlignedObjectArray<btAlignedObjectArray<btDeformableFaceNodeContactConstraint> > m_deformableConstraints;
	// map from node index to node anchor constraint
	btAlignedObjectArray<btAlignedObjectArray<btDeformableNodeAnchorConstraint> > m_nodeAnchorConstraints;

	bool m_useStrainLimiting;

	btDeformableContactProjection(btAlignedObjectArray<btSoftBody*>& softBodies)
		: m_softBodies(softBodies)
	{
	}

	virtual ~btDeformableContactProjection()
	{
	}

	// apply the constraints to the rhs of the linear solve
	virtual void project(TVStack& x);

	// add friction force to the rhs of the linear solve
	virtual void applyDynamicFriction(TVStack& f);

	// update and solve the constraints
	virtual btScalar update(btCollisionObject** deformableBodies, int numDeformableBodies, const btContactSolverInfo& infoGlobal);

	// Add constraints to m_constraints. In addition, the constraints that each vertex own are recorded in m_constraintsDict.
	virtual void setConstraints(const btContactSolverInfo& infoGlobal);

	// Set up projections for each vertex by adding the projection direction to
	virtual void setProjection();

	virtual void reinitialize(bool nodeUpdated);

	btScalar solveSplitImpulse(btCollisionObject** deformableBodies, int numDeformableBodies, const btContactSolverInfo& infoGlobal);

	virtual void setLagrangeMultiplier();

	void checkConstraints(const TVStack& x);
};
#endif /* btDeformableContactProjection_h */
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

#ifndef BT_COROTATED_H
#define BT_COROTATED_H

#include "btDeformableLagrangianForce.h"
#include "LinearMath/btPolarDecomposition.h"

static inline int PolarDecomposition(const btMatrix3x3& m, btMatrix3x3& q, btMatrix3x3& s)
{
	static const btPolarDecomposition polar;
	return polar.decompose(m, q, s);
}

class btDeformableCorotatedForce : public btDeformableLagrangianForce
{
public:
	typedef btAlignedObjectArray<btVector3> TVStack;
	btScalar m_mu, m_lambda;
	btDeformableCorotatedForce() : m_mu(1), m_lambda(1)
	{
	}

	btDeformableCorotatedForce(btScalar mu, btScalar lambda) : m_mu(mu), m_lambda(lambda)
	{
	}

	virtual void addScaledForces(btScalar scale, TVStack& force)
	{
		addScaledElasticForce(scale, force);
	}

	virtual void addScaledExplicitForce(btScalar scale, TVStack& force)
	{
		addScaledElasticForce(scale, force);
	}

	virtual void addScaledDampingForce(btScalar scale, TVStack& force)
	{
	}

	virtual void addScaledElasticForce(btScalar scale, TVStack& force)
	{
		int numNodes = getNumNodes();
		btAssert(numNodes <= force.size());
		btVector3 grad_N_hat_1st_col = btVector3(-1, -1, -1);
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			for (int j = 0; j < psb->m_tetras.size(); ++j)
			{
				btSoftBody::Tetra& tetra = psb->m_tetras[j];
				btMatrix3x3 P;
				firstPiola(tetra.m_F, P);
				btVector3 force_on_node0 = P * (tetra.m_Dm_inverse.transpose() * grad_N_hat_1st_col);
				btMatrix3x3 force_on_node123 = P * tetra.m_Dm_inverse.transpose();

				btSoftBody::Node* node0 = tetra.m_n[0];
				btSoftBody::Node* node1 = tetra.m_n[1];
				btSoftBody::Node* node2 = tetra.m_n[2];
				btSoftBody::Node* node3 = tetra.m_n[3];
				size_t id0 = node0->index;
				size_t id1 = node1->index;
				size_t id2 = node2->index;
				size_t id3 = node3->index;

				// elastic force
				// explicit elastic force
				btScalar scale1 = scale * tetra.m_element_measure;
				force[id0] -= scale1 * force_on_node0;
				force[id1] -= scale1 * force_on_node123.getColumn(0);
				force[id2] -= scale1 * force_on_node123.getColumn(1);
				force[id3] -= scale1 * force_on_node123.getColumn(2);
			}
		}
	}

	void firstPiola(const btMatrix3x3& F, btMatrix3x3& P)
	{
		// btMatrix3x3 JFinvT = F.adjoint();
		btScalar J = F.determinant();
		P = F.adjoint().transpose() * (m_lambda * (J - 1));
		if (m_mu > SIMD_EPSILON)
		{
			btMatrix3x3 R, S;
			if (J < 1024 * SIMD_EPSILON)
				R.setIdentity();
			else
				PolarDecomposition(F, R, S);  // this QR is not robust, consider using implicit shift svd
			/*https://fuchuyuan.github.io/research/svd/paper.pdf*/
			P += (F - R) * 2 * m_mu;
		}
	}

	virtual void addScaledElasticForceDifferential(btScalar scale, const TVStack& dx, TVStack& df)
	{
	}

	virtual void addScaledDampingForceDifferential(btScalar scale, const TVStack& dv, TVStack& df)
	{
	}

	virtual void buildDampingForceDifferentialDiagonal(btScalar scale, TVStack& diagA) {}

	virtual btDeformableLagrangianForceType getForceType()
	{
		return BT_COROTATED_FORCE;
	}
};

#endif /* btCorotated_h */
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

#ifndef BT_DEFORMABLE_GRAVITY_FORCE_H
#define BT_DEFORMABLE_GRAVITY_FORCE_H

#include "btDeformableLagrangianForce.h"

class btDeformableGravityForce : public btDeformableLagrangianForce
{
public:
	typedef btAlignedObjectArray<btVector3> TVStack;
	btVector3 m_gravity;

	btDeformableGravityForce(const btVector3& g) : m_gravity(g)
	{
	}

	virtual void addScaledForces(btScalar scale, TVStack& force)
	{
		addScaledGravityForce(scale, force);
	}

	virtual void addScaledExplicitForce(btScalar scale, TVStack& force)
	{
		addScaledGravityForce(scale, force);
	}

	virtual void addScaledDampingForce(btScalar scale, TVStack& force)
	{
	}

	virtual void addScaledElasticForceDifferential(btScalar scale, const TVStack& dx, TVStack& df)
	{
	}

	virtual void addScaledDampingForceDifferential(btScalar scale, const TVStack& dv, TVStack& df)
	{
	}

	virtual void buildDampingForceDifferentialDiagonal(btScalar scale, TVStack& diagA) {}

	virtual void addScaledGravityForce(btScalar scale, TVStack& force)
	{
		int numNodes = getNumNodes();
		btAssert(numNodes <= force.size());
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				btSoftBody::Node& n = psb->m_nodes[j];
				size_t id = n.index;
				btScalar mass = (n.m_im == 0) ? 0 : 1. / n.m_im;
				btVector3 scaled_force = scale * m_gravity * mass * m_softBodies[i]->m_gravityFactor;
				force[id] += scaled_force;
			}
		}
	}

	virtual btDeformableLagrangianForceType getForceType()
	{
		return BT_GRAVITY_FORCE;
	}

	// the gravitational potential energy
	virtual double totalEnergy(btScalar dt)
	{
		double e = 0;
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				const btSoftBody::Node& node = psb->m_nodes[j];
				if (node.m_im > 0)
				{
					e -= m_gravity.dot(node.m_q) / node.m_im;
				}
			}
		}
		return e;
	}
};
#endif /* BT_DEFORMABLE_GRAVITY_FORCE_H */
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

#ifndef BT_DEFORMABLE_LAGRANGIAN_FORCE_H
#define BT_DEFORMABLE_LAGRANGIAN_FORCE_H

#include "btSoftBody.h"
#include <LinearMath/btHashMap.h>
#include <iostream>

enum btDeformableLagrangianForceType
{
	BT_GRAVITY_FORCE = 1,
	BT_MASSSPRING_FORCE = 2,
	BT_COROTATED_FORCE = 3,
	BT_NEOHOOKEAN_FORCE = 4,
	BT_LINEAR_ELASTICITY_FORCE = 5,
	BT_MOUSE_PICKING_FORCE = 6
};

static inline double randomDouble(double low, double high)
{
	return low + static_cast<double>(rand()) / RAND_MAX * (high - low);
}

class btDeformableLagrangianForce
{
public:
	typedef btAlignedObjectArray<btVector3> TVStack;
	btAlignedObjectArray<btSoftBody*> m_softBodies;
	const btAlignedObjectArray<btSoftBody::Node*>* m_nodes;

	btDeformableLagrangianForce()
	{
	}

	virtual ~btDeformableLagrangianForce() {}

	// add all forces
	virtual void addScaledForces(btScalar scale, TVStack& force) = 0;

	// add damping df
	virtual void addScaledDampingForceDifferential(btScalar scale, const TVStack& dv, TVStack& df) = 0;

	// build diagonal of A matrix
	virtual void buildDampingForceDifferentialDiagonal(btScalar scale, TVStack& diagA) = 0;

	// add elastic df
	virtual void addScaledElasticForceDifferential(btScalar scale, const TVStack& dx, TVStack& df) = 0;

	// add all forces that are explicit in explicit solve
	virtual void addScaledExplicitForce(btScalar scale, TVStack& force) = 0;

	// add all damping forces
	virtual void addScaledDampingForce(btScalar scale, TVStack& force) = 0;

	virtual void addScaledHessian(btScalar scale) {}

	virtual btDeformableLagrangianForceType getForceType() = 0;

	virtual void reinitialize(bool nodeUpdated)
	{
	}

	// get number of nodes that have the force
	virtual int getNumNodes()
	{
		int numNodes = 0;
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			numNodes += m_softBodies[i]->m_nodes.size();
		}
		return numNodes;
	}

	// add a soft body to be affected by the particular lagrangian force
	virtual void addSoftBody(btSoftBody* psb)
	{
		m_softBodies.push_back(psb);
	}

	virtual void removeSoftBody(btSoftBody* psb)
	{
		m_softBodies.remove(psb);
	}

	virtual const btAlignedObjectArray<btSoftBody*>& getSoftBodies() const
	{
		return m_softBodies;
	}

	virtual void setIndices(const btAlignedObjectArray<btSoftBody::Node*>* nodes)
	{
		m_nodes = nodes;
	}

	// Calculate the incremental deformable generated from the input dx
	virtual btMatrix3x3 Ds(int id0, int id1, int id2, int id3, const TVStack& dx)
	{
		btVector3 c1 = dx[id1] - dx[id0];
		btVector3 c2 = dx[id2] - dx[id0];
		btVector3 c3 = dx[id3] - dx[id0];
		return btMatrix3x3(c1, c2, c3).transpose();
	}

	// Calculate the incremental deformable generated from the current velocity
	virtual btMatrix3x3 DsFromVelocity(const btSoftBody::Node* n0, const btSoftBody::Node* n1, const btSoftBody::Node* n2, const btSoftBody::Node* n3)
	{
		btVector3 c1 = n1->m_v - n0->m_v;
		btVector3 c2 = n2->m_v - n0->m_v;
		btVector3 c3 = n3->m_v - n0->m_v;
		return btMatrix3x3(c1, c2, c3).transpose();
	}

	// test for addScaledElasticForce function
	virtual void testDerivative()
	{
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				psb->m_nodes[j].m_q += btVector3(randomDouble(-.1, .1), randomDouble(-.1, .1), randomDouble(-.1, .1));
			}
			psb->updateDeformation();
		}

		TVStack dx;
		dx.resize(getNumNodes());
		TVStack dphi_dx;
		dphi_dx.resize(dx.size());
		for (int i = 0; i < dphi_dx.size(); ++i)
		{
			dphi_dx[i].setZero();
		}
		addScaledForces(-1, dphi_dx);

		// write down the current position
		TVStack x;
		x.resize(dx.size());
		int counter = 0;
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				x[counter] = psb->m_nodes[j].m_q;
				counter++;
			}
		}
		counter = 0;

		// populate dx with random vectors
		for (int i = 0; i < dx.size(); ++i)
		{
			dx[i].setX(randomDouble(-1, 1));
			dx[i].setY(randomDouble(-1, 1));
			dx[i].setZ(randomDouble(-1, 1));
		}

		btAlignedObjectArray<double> errors;
		for (int it = 0; it < 10; ++it)
		{
			for (int i = 0; i < dx.size(); ++i)
			{
				dx[i] *= 0.5;
			}

			// get dphi/dx * dx
			double dphi = 0;
			for (int i = 0; i < dx.size(); ++i)
			{
				dphi += dphi_dx[i].dot(dx[i]);
			}

			for (int i = 0; i < m_softBodies.size(); ++i)
			{
				btSoftBody* psb = m_softBodies[i];
				for (int j = 0; j < psb->m_nodes.size(); ++j)
				{
					psb->m_nodes[j].m_q = x[counter] + dx[counter];
					counter++;
				}
				psb->updateDeformation();
			}
			counter = 0;
			double f1 = totalElasticEnergy(0);

			for (int i = 0; i < m_softBodies.size(); ++i)
			{
				btSoftBody* psb = m_softBodies[i];
				for (int j = 0; j < psb->m_nodes.size(); ++j)
				{
					psb->m_nodes[j].m_q = x[counter] - dx[counter];
					counter++;
				}
				psb->updateDeformation();
			}
			counter = 0;

			double f2 = totalElasticEnergy(0);

			//restore m_q
			for (int i = 0; i < m_softBodies.size(); ++i)
			{
				btSoftBody* psb = m_softBodies[i];
				for (int j = 0; j < psb->m_nodes.size(); ++j)
				{
					psb->m_nodes[j].m_q = x[counter];
					counter++;
				}
				psb->updateDeformation();
			}
			counter = 0;
			double error = f1 - f2 - 2 * dphi;
			errors.push_back(error);
			std::cout << "Iteration = " << it << ", f1 = " << f1 << ", f2 = " << f2 << ", error = " << error << std::endl;
		}
		for (int i = 1; i < errors.size(); ++i)
		{
			std::cout << "Iteration = " << i << ", ratio = " << errors[i - 1] / errors[i] << std::endl;
		}
	}

	// test for addScaledElasticForce function
	virtual void testHessian()
	{
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				psb->m_nodes[j].m_q += btVector3(randomDouble(-.1, .1), randomDouble(-.1, .1), randomDouble(-.1, .1));
			}
			psb->updateDeformation();
		}

		TVStack dx;
		dx.resize(getNumNodes());
		TVStack df;
		df.resize(dx.size());
		TVStack f1;
		f1.resize(dx.size());
		TVStack f2;
		f2.resize(dx.size());

		// write down the current position
		TVStack x;
		x.resize(dx.size());
		int counter = 0;
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				x[counter] = psb->m_nodes[j].m_q;
				counter++;
			}
		}
		counter = 0;

		// populate dx with random vectors
		for (int i = 0; i < dx.size(); ++i)
		{
			dx[i].setX(randomDouble(-1, 1));
			dx[i].setY(randomDouble(-1, 1));
			dx[i].setZ(randomDouble(-1, 1));
		}

		btAlignedObjectArray<double> errors;
		for (int it = 0; it < 10; ++it)
		{
			for (int i = 0; i < dx.size(); ++i)
			{
				dx[i] *= 0.5;
			}

			// get df
			for (int i = 0; i < df.size(); ++i)
			{
				df[i].setZero();
				f1[i].setZero();
				f2[i].setZero();
			}

			//set df
			addScaledElasticForceDifferential(-1, dx, df);

			for (int i = 0; i < m_softBodies.size(); ++i)
			{
				btSoftBody* psb = m_softBodies[i];
				for (int j = 0; j < psb->m_nodes.size(); ++j)
				{
					psb->m_nodes[j].m_q = x[counter] + dx[counter];
					counter++;
				}
				psb->updateDeformation();
			}
			counter = 0;

			//set f1
			addScaledForces(-1, f1);

			for (int i = 0; i < m_softBodies.size(); ++i)
			{
				btSoftBody* psb = m_softBodies[i];
				for (int j = 0; j < psb->m_nodes.size(); ++j)
				{
					psb->m_nodes[j].m_q = x[counter] - dx[counter];
					counter++;
				}
				psb->updateDeformation();
			}
			counter = 0;

			//set f2
			addScaledForces(-1, f2);

			//restore m_q
			for (int i = 0; i < m_softBodies.size(); ++i)
			{
				btSoftBody* psb = m_softBodies[i];
				for (int j = 0; j < psb->m_nodes.size(); ++j)
				{
					psb->m_nodes[j].m_q = x[counter];
					counter++;
				}
				psb->updateDeformation();
			}
			counter = 0;
			double error = 0;
			for (int i = 0; i < df.size(); ++i)
			{
				btVector3 error_vector = f1[i] - f2[i] - 2 * df[i];
				error += error_vector.length2();
			}
			error = btSqrt(error);
			errors.push_back(error);
			std::cout << "Iteration = " << it << ", error = " << error << std::endl;
		}
		for (int i = 1; i < errors.size(); ++i)
		{
			std::cout << "Iteration = " << i << ", ratio = " << errors[i - 1] / errors[i] << std::endl;
		}
	}

	//
	virtual double totalElasticEnergy(btScalar dt)
	{
		return 0;
	}

	//
	virtual double totalDampingEnergy(btScalar dt)
	{
		return 0;
	}

	// total Energy takes dt as input because certain energies depend on dt
	virtual double totalEnergy(btScalar dt)
	{
		return totalElasticEnergy(dt) + totalDampingEnergy(dt);
	}

	virtual btScalar getYoungsModulus() const
	{
		return -1.0;
	}

	virtual btScalar getPoissonRatio() const
	{
		return -1.0;
	}
};
#endif /* BT_DEFORMABLE_LAGRANGIAN_FORCE */
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

#ifndef BT_LINEAR_ELASTICITY_H
#define BT_LINEAR_ELASTICITY_H

#include "btDeformableLagrangianForce.h"
#include "LinearMath/btQuickprof.h"
#include "LinearMath/btImplicitQRSVD.h"
#include "btSoftBodyInternals.h"
#define TETRA_FLAT_THRESHOLD 0.01
class btDeformableLinearElasticityForce : public btDeformableLagrangianForce
{
public:
	typedef btAlignedObjectArray<btVector3> TVStack;
	btScalar m_mu, m_lambda;
	btScalar m_E, m_nu;  // Young's modulus and Poisson ratio
	btScalar m_damping_alpha, m_damping_beta;
	btDeformableLinearElasticityForce() : m_mu(1), m_lambda(1), m_damping_alpha(0.01), m_damping_beta(0.01)
	{
		updateYoungsModulusAndPoissonRatio();
	}

	btDeformableLinearElasticityForce(btScalar mu, btScalar lambda, btScalar damping_alpha = 0.01, btScalar damping_beta = 0.01) : m_mu(mu), m_lambda(lambda), m_damping_alpha(damping_alpha), m_damping_beta(damping_beta)
	{
		updateYoungsModulusAndPoissonRatio();
	}

	void updateYoungsModulusAndPoissonRatio()
	{
		// conversion from Lame Parameters to Young's modulus and Poisson ratio
		// https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
		m_E = m_mu * (3 * m_lambda + 2 * m_mu) / (m_lambda + m_mu);
		m_nu = m_lambda * 0.5 / (m_mu + m_lambda);
	}

	void updateLameParameters()
	{
		// conversion from Young's modulus and Poisson ratio to Lame Parameters
		// https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
		m_mu = m_E * 0.5 / (1 + m_nu);
		m_lambda = m_E * m_nu / ((1 + m_nu) * (1 - 2 * m_nu));
	}

	void setYoungsModulus(btScalar E)
	{
		m_E = E;
		updateLameParameters();
	}

	void setPoissonRatio(btScalar nu)
	{
		m_nu = nu;
		updateLameParameters();
	}

	void setDamping(btScalar damping_alpha, btScalar damping_beta)
	{
		m_damping_alpha = damping_alpha;
		m_damping_beta = damping_beta;
	}

	void setLameParameters(btScalar mu, btScalar lambda)
	{
		m_mu = mu;
		m_lambda = lambda;
		updateYoungsModulusAndPoissonRatio();
	}

	virtual void addScaledForces(btScalar scale, TVStack& force)
	{
		addScaledDampingForce(scale, force);
		addScaledElasticForce(scale, force);
	}

	virtual void addScaledExplicitForce(btScalar scale, TVStack& force)
	{
		addScaledElasticForce(scale, force);
	}

	// The damping matrix is calculated using the time n state as described in https://www.math.ucla.edu/~jteran/papers/GSSJT15.pdf to allow line search
	virtual void addScaledDampingForce(btScalar scale, TVStack& force)
	{
		if (m_damping_alpha == 0 && m_damping_beta == 0)
			return;
		btScalar mu_damp = m_damping_beta * m_mu;
		btScalar lambda_damp = m_damping_beta * m_lambda;
		int numNodes = getNumNodes();
		btAssert(numNodes <= force.size());
		btVector3 grad_N_hat_1st_col = btVector3(-1, -1, -1);
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_tetras.size(); ++j)
			{
				bool close_to_flat = (psb->m_tetraScratches[j].m_J < TETRA_FLAT_THRESHOLD);
				btSoftBody::Tetra& tetra = psb->m_tetras[j];
				btSoftBody::Node* node0 = tetra.m_n[0];
				btSoftBody::Node* node1 = tetra.m_n[1];
				btSoftBody::Node* node2 = tetra.m_n[2];
				btSoftBody::Node* node3 = tetra.m_n[3];
				size_t id0 = node0->index;
				size_t id1 = node1->index;
				size_t id2 = node2->index;
				size_t id3 = node3->index;
				btMatrix3x3 dF = DsFromVelocity(node0, node1, node2, node3) * tetra.m_Dm_inverse;
				if (!close_to_flat)
				{
					dF = psb->m_tetraScratches[j].m_corotation.transpose() * dF;
				}
				btMatrix3x3 I;
				I.setIdentity();
				btMatrix3x3 dP = (dF + dF.transpose()) * mu_damp + I * ((dF[0][0] + dF[1][1] + dF[2][2]) * lambda_damp);
				btMatrix3x3 df_on_node123 = dP * tetra.m_Dm_inverse.transpose();
				if (!close_to_flat)
				{
					df_on_node123 = psb->m_tetraScratches[j].m_corotation * df_on_node123;
				}
				btVector3 df_on_node0 = df_on_node123 * grad_N_hat_1st_col;
				// damping force differential
				btScalar scale1 = scale * tetra.m_element_measure;
				force[id0] -= scale1 * df_on_node0;
				force[id1] -= scale1 * df_on_node123.getColumn(0);
				force[id2] -= scale1 * df_on_node123.getColumn(1);
				force[id3] -= scale1 * df_on_node123.getColumn(2);
			}
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				const btSoftBody::Node& node = psb->m_nodes[j];
				size_t id = node.index;
				if (node.m_im > 0)
				{
					force[id] -= scale * node.m_v / node.m_im * m_damping_alpha;
				}
				//fprintf(stderr, "force[id] %d %f %f %f scale %f node.m_v %f %f %f node.m_im %f m_damping_alpha %f\n", j, force[id].x(), force[id].y(), force[id].z(), scale, node.m_v.x(), node.m_v.y(), node.m_v.z(), node.m_im, m_damping_alpha);
			}
		}
	}

	virtual double totalElasticEnergy(btScalar dt)
	{
		double energy = 0;
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_tetraScratches.size(); ++j)
			{
				btSoftBody::Tetra& tetra = psb->m_tetras[j];
				btSoftBody::TetraScratch& s = psb->m_tetraScratches[j];
				energy += tetra.m_element_measure * elasticEnergyDensity(s);
			}
		}
		return energy;
	}

	// The damping energy is formulated as in https://www.math.ucla.edu/~jteran/papers/GSSJT15.pdf to allow line search
	virtual double totalDampingEnergy(btScalar dt)
	{
		double energy = 0;
		int sz = 0;
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				sz = btMax(sz, psb->m_nodes[j].index);
			}
		}
		TVStack dampingForce;
		dampingForce.resize(sz + 1);
		for (int i = 0; i < dampingForce.size(); ++i)
			dampingForce[i].setZero();
		addScaledDampingForce(0.5, dampingForce);
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				const btSoftBody::Node& node = psb->m_nodes[j];
				energy -= dampingForce[node.index].dot(node.m_v) / dt;
			}
		}
		return energy;
	}

	double elasticEnergyDensity(const btSoftBody::TetraScratch& s)
	{
		double density = 0;
		btMatrix3x3 epsilon = (s.m_F + s.m_F.transpose()) * 0.5 - btMatrix3x3::getIdentity();
		btScalar trace = epsilon[0][0] + epsilon[1][1] + epsilon[2][2];
		density += m_mu * (epsilon[0].length2() + epsilon[1].length2() + epsilon[2].length2());
		density += m_lambda * trace * trace * 0.5;
		return density;
	}

	virtual void addScaledElasticForce(btScalar scale, TVStack& force)
	{
		int numNodes = getNumNodes();
		btAssert(numNodes <= force.size());
		btVector3 grad_N_hat_1st_col = btVector3(-1, -1, -1);
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			btScalar max_p = psb->m_cfg.m_maxStress;
			btScalar averagePrincipalStress = 0.0;
			for (int j = 0; j < psb->m_tetras.size(); ++j)
			{
				btSoftBody::Tetra& tetra = psb->m_tetras[j];
				btMatrix3x3 P;
				firstPiola(psb->m_tetraScratches[j], P);
				btScalar trPTP = (P[0].length() + P[1].length() + P[2].length());
				averagePrincipalStress += trPTP;
#define USE_SVD 1
#if USE_SVD
				if (max_p > 0)
				{
					// since we want to clamp the principal stress to max_p, we only need to
					// calculate SVD when sigma_0^2 + sigma_1^2 + sigma_2^2 > max_p * max_p
					btScalar trPTP = (P[0].length2() + P[1].length2() + P[2].length2());
					if (trPTP > max_p * max_p)
					{
						btMatrix3x3 U, V;
						btVector3 sigma;
						singularValueDecomposition(P, U, sigma, V);
						sigma[0] = btMin(sigma[0], max_p);
						sigma[1] = btMin(sigma[1], max_p);
						sigma[2] = btMin(sigma[2], max_p);
						sigma[0] = btMax(sigma[0], -max_p);
						sigma[1] = btMax(sigma[1], -max_p);
						sigma[2] = btMax(sigma[2], -max_p);
						btMatrix3x3 Sigma;
						Sigma.setIdentity();
						Sigma[0][0] = sigma[0];
						Sigma[1][1] = sigma[1];
						Sigma[2][2] = sigma[2];
						P = U * Sigma * V.transpose();
					}
				}
#endif
				//                btVector3 force_on_node0 = P * (tetra.m_Dm_inverse.transpose()*grad_N_hat_1st_col);
				btMatrix3x3 force_on_node123 = psb->m_tetraScratches[j].m_corotation * P * tetra.m_Dm_inverse.transpose();
				btVector3 force_on_node0 = force_on_node123 * grad_N_hat_1st_col;

				btSoftBody::Node* node0 = tetra.m_n[0];
				btSoftBody::Node* node1 = tetra.m_n[1];
				btSoftBody::Node* node2 = tetra.m_n[2];
				btSoftBody::Node* node3 = tetra.m_n[3];
				size_t id0 = node0->index;
				size_t id1 = node1->index;
				size_t id2 = node2->index;
				size_t id3 = node3->index;

				// elastic force
				btScalar scale1 = scale * tetra.m_element_measure;
				force[id0] -= scale1 * force_on_node0;
				force[id1] -= scale1 * force_on_node123.getColumn(0);
				force[id2] -= scale1 * force_on_node123.getColumn(1);
				force[id3] -= scale1 * force_on_node123.getColumn(2);
			}
			averagePrincipalStress /= psb->m_tetras.size();
			psb->m_averagePrincipalStress = averagePrincipalStress;
		}
	}

	virtual void buildDampingForceDifferentialDiagonal(btScalar scale, TVStack& diagA) {}

	// The damping matrix is calculated using the time n state as described in https://www.math.ucla.edu/~jteran/papers/GSSJT15.pdf to allow line search
	virtual void addScaledDampingForceDifferential(btScalar scale, const TVStack& dv, TVStack& df)
	{
		if (m_damping_alpha == 0 && m_damping_beta == 0)
			return;
		btScalar mu_damp = m_damping_beta * m_mu;
		btScalar lambda_damp = m_damping_beta * m_lambda;
		int numNodes = getNumNodes();
		btAssert(numNodes <= df.size());
		btVector3 grad_N_hat_1st_col = btVector3(-1, -1, -1);
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_tetras.size(); ++j)
			{
				bool close_to_flat = (psb->m_tetraScratches[j].m_J < TETRA_FLAT_THRESHOLD);
				btSoftBody::Tetra& tetra = psb->m_tetras[j];
				btSoftBody::Node* node0 = tetra.m_n[0];
				btSoftBody::Node* node1 = tetra.m_n[1];
				btSoftBody::Node* node2 = tetra.m_n[2];
				btSoftBody::Node* node3 = tetra.m_n[3];
				size_t id0 = node0->index;
				size_t id1 = node1->index;
				size_t id2 = node2->index;
				size_t id3 = node3->index;
				btMatrix3x3 dF = Ds(id0, id1, id2, id3, dv) * tetra.m_Dm_inverse;
				if (!close_to_flat)
				{
					dF = psb->m_tetraScratches[j].m_corotation.transpose() * dF;
				}
				btMatrix3x3 I;
				I.setIdentity();
				btMatrix3x3 dP = (dF + dF.transpose()) * mu_damp + I * ((dF[0][0] + dF[1][1] + dF[2][2]) * lambda_damp);
				btMatrix3x3 df_on_node123 = dP * tetra.m_Dm_inverse.transpose();
				if (!close_to_flat)
				{
					df_on_node123 = psb->m_tetraScratches[j].m_corotation * df_on_node123;
				}
				btVector3 df_on_node0 = df_on_node123 * grad_N_hat_1st_col;

				// damping force differential
				btScalar scale1 = scale * tetra.m_element_measure;
				df[id0] -= scale1 * df_on_node0;
				df[id1] -= scale1 * df_on_node123.getColumn(0);
				df[id2] -= scale1 * df_on_node123.getColumn(1);
				df[id3] -= scale1 * df_on_node123.getColumn(2);
			}
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				const btSoftBody::Node& node = psb->m_nodes[j];
				size_t id = node.index;
				if (node.m_im > 0)
				{
					df[id] -= scale * dv[id] / node.m_im * m_damping_alpha;
				}
			}
		}
	}

	virtual void addScaledElasticForceDifferential(btScalar scale, const TVStack& dx, TVStack& df)
	{
		int numNodes = getNumNodes();
		btAssert(numNodes <= df.size());
		btVector3 grad_N_hat_1st_col = btVector3(-1, -1, -1);
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_tetras.size(); ++j)
			{
				btSoftBody::Tetra& tetra = psb->m_tetras[j];
				btSoftBody::Node* node0 = tetra.m_n[0];
				btSoftBody::Node* node1 = tetra.m_n[1];
				btSoftBody::Node* node2 = tetra.m_n[2];
				btSoftBody::Node* node3 = tetra.m_n[3];
				size_t id0 = node0->index;
				size_t id1 = node1->index;
				size_t id2 = node2->index;
				size_t id3 = node3->index;
				btMatrix3x3 dF = psb->m_tetraScratches[j].m_corotation.transpose() * Ds(id0, id1, id2, id3, dx) * tetra.m_Dm_inverse;
				btMatrix3x3 dP;
				firstPiolaDifferential(psb->m_tetraScratches[j], dF, dP);
				//                btVector3 df_on_node0 = dP * (tetra.m_Dm_inverse.transpose()*grad_N_hat_1st_col);
				btMatrix3x3 df_on_node123 = psb->m_tetraScratches[j].m_corotation * dP * tetra.m_Dm_inverse.transpose();
				btVector3 df_on_node0 = df_on_node123 * grad_N_hat_1st_col;

				// elastic force differential
				btScalar scale1 = scale * tetra.m_element_measure;
				df[id0] -= scale1 * df_on_node0;
				df[id1] -= scale1 * df_on_node123.getColumn(0);
				df[id2] -= scale1 * df_on_node123.getColumn(1);
				df[id3] -= scale1 * df_on_node123.getColumn(2);
			}
		}
	}

	void firstPiola(const btSoftBody::TetraScratch& s, btMatrix3x3& P)
	{
		btMatrix3x3 corotated_F = s.m_corotation.transpose() * s.m_F;

		btMatrix3x3 epsilon = (corotated_F + corotated_F.transpose()) * 0.5 - btMatrix3x3::getIdentity();
		btScalar trace = epsilon[0][0] + epsilon[1][1] + epsilon[2][2];
		P = epsilon * btScalar(2) * m_mu + btMatrix3x3::getIdentity() * m_lambda * trace;
	}

	// Let P be the first piola stress.
	// This function calculates the dP = dP/dF * dF
	void firstPiolaDifferential(const btSoftBody::TetraScratch& s, const btMatrix3x3& dF, btMatrix3x3& dP)
	{
		btScalar trace = (dF[0][0] + dF[1][1] + dF[2][2]);
		dP = (dF + dF.transpose()) * m_mu + btMatrix3x3::getIdentity() * m_lambda * trace;
	}

	// Let Q be the damping stress.
	// This function calculates the dP = dQ/dF * dF
	void firstPiolaDampingDifferential(const btSoftBody::TetraScratch& s, const btMatrix3x3& dF, btMatrix3x3& dP)
	{
		btScalar mu_damp = m_damping_beta * m_mu;
		btScalar lambda_damp = m_damping_beta * m_lambda;
		btScalar trace = (dF[0][0] + dF[1][1] + dF[2][2]);
		dP = (dF + dF.transpose()) * mu_damp + btMatrix3x3::getIdentity() * lambda_damp * trace;
	}

	virtual void addScaledHessian(btScalar scale)
	{
		btVector3 grad_N_hat_1st_col = btVector3(-1, -1, -1);
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_tetras.size(); ++j)
			{
				btSoftBody::Tetra& tetra = psb->m_tetras[j];
				btMatrix3x3 P;
				firstPiola(psb->m_tetraScratches[j], P);  // make sure scratch is evaluated at x_n + dt * vn
				btMatrix3x3 force_on_node123 = psb->m_tetraScratches[j].m_corotation * P * tetra.m_Dm_inverse.transpose();
				btVector3 force_on_node0 = force_on_node123 * grad_N_hat_1st_col;
				btSoftBody::Node* node0 = tetra.m_n[0];
				btSoftBody::Node* node1 = tetra.m_n[1];
				btSoftBody::Node* node2 = tetra.m_n[2];
				btSoftBody::Node* node3 = tetra.m_n[3];
				btScalar scale1 = scale * (scale + m_damping_beta) * tetra.m_element_measure;  // stiff and stiffness-damping terms;
				node0->m_effectiveMass += OuterProduct(force_on_node0, force_on_node0) * scale1;
				node1->m_effectiveMass += OuterProduct(force_on_node123.getColumn(0), force_on_node123.getColumn(0)) * scale1;
				node2->m_effectiveMass += OuterProduct(force_on_node123.getColumn(1), force_on_node123.getColumn(1)) * scale1;
				node3->m_effectiveMass += OuterProduct(force_on_node123.getColumn(2), force_on_node123.getColumn(2)) * scale1;
			}
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				btSoftBody::Node& node = psb->m_nodes[j];
				if (node.m_im > 0)
				{
					btMatrix3x3 I;
					I.setIdentity();
					node.m_effectiveMass += I * (scale * (1.0 / node.m_im) * m_damping_alpha);
				}
			}
		}
	}

	virtual btScalar getYoungsModulus() const
	{
		return m_E;
	}

	virtual btScalar getPoissonRatio() const
	{
		return m_nu;
	}

	virtual btDeformableLagrangianForceType getForceType()
	{
		return BT_LINEAR_ELASTICITY_FORCE;
	}
};
#endif /* BT_LINEAR_ELASTICITY_H */
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

#ifndef BT_MASS_SPRING_H
#define BT_MASS_SPRING_H

#include "btDeformableLagrangianForce.h"

class btDeformableMassSpringForce : public btDeformableLagrangianForce
{
	// If true, the damping force will be in the direction of the spring
	// If false, the damping force will be in the direction of the velocity
	bool m_momentum_conserving;
	btScalar m_elasticStiffness, m_dampingStiffness, m_bendingStiffness;

public:
	typedef btAlignedObjectArray<btVector3> TVStack;
	btDeformableMassSpringForce() : m_momentum_conserving(false), m_elasticStiffness(1), m_dampingStiffness(0.05)
	{
	}
	btDeformableMassSpringForce(btScalar k, btScalar d, bool conserve_angular = true, double bending_k = -1) : m_momentum_conserving(conserve_angular), m_elasticStiffness(k), m_dampingStiffness(d), m_bendingStiffness(bending_k)
	{
		if (m_bendingStiffness < btScalar(0))
		{
			m_bendingStiffness = m_elasticStiffness;
		}
	}

	virtual void addScaledForces(btScalar scale, TVStack& force)
	{
		addScaledDampingForce(scale, force);
		addScaledElasticForce(scale, force);
	}

	virtual void addScaledExplicitForce(btScalar scale, TVStack& force)
	{
		addScaledElasticForce(scale, force);
	}

	virtual void addScaledDampingForce(btScalar scale, TVStack& force)
	{
		int numNodes = getNumNodes();
		btAssert(numNodes <= force.size());
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			const btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_links.size(); ++j)
			{
				const btSoftBody::Link& link = psb->m_links[j];
				btSoftBody::Node* node1 = link.m_n[0];
				btSoftBody::Node* node2 = link.m_n[1];
				size_t id1 = node1->index;
				size_t id2 = node2->index;

				// damping force
				btVector3 v_diff = (node2->m_v - node1->m_v);
				btVector3 scaled_force = scale * m_dampingStiffness * v_diff;
				if (m_momentum_conserving)
				{
					if ((node2->m_x - node1->m_x).norm() > SIMD_EPSILON)
					{
						btVector3 dir = (node2->m_x - node1->m_x).normalized();
						scaled_force = scale * m_dampingStiffness * v_diff.dot(dir) * dir;
					}
				}
				force[id1] += scaled_force;
				force[id2] -= scaled_force;
			}
		}
	}

	virtual void addScaledElasticForce(btScalar scale, TVStack& force)
	{
		int numNodes = getNumNodes();
		btAssert(numNodes <= force.size());
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			const btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_links.size(); ++j)
			{
				const btSoftBody::Link& link = psb->m_links[j];
				btSoftBody::Node* node1 = link.m_n[0];
				btSoftBody::Node* node2 = link.m_n[1];
				btScalar r = link.m_rl;
				size_t id1 = node1->index;
				size_t id2 = node2->index;

				// elastic force
				btVector3 dir = (node2->m_q - node1->m_q);
				btVector3 dir_normalized = (dir.norm() > SIMD_EPSILON) ? dir.normalized() : btVector3(0, 0, 0);
				btScalar scaled_stiffness = scale * (link.m_bbending ? m_bendingStiffness : m_elasticStiffness);
				btVector3 scaled_force = scaled_stiffness * (dir - dir_normalized * r);
				force[id1] += scaled_force;
				force[id2] -= scaled_force;
			}
		}
	}

	virtual void addScaledDampingForceDifferential(btScalar scale, const TVStack& dv, TVStack& df)
	{
		// implicit damping force differential
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			btScalar scaled_k_damp = m_dampingStiffness * scale;
			for (int j = 0; j < psb->m_links.size(); ++j)
			{
				const btSoftBody::Link& link = psb->m_links[j];
				btSoftBody::Node* node1 = link.m_n[0];
				btSoftBody::Node* node2 = link.m_n[1];
				size_t id1 = node1->index;
				size_t id2 = node2->index;

				btVector3 local_scaled_df = scaled_k_damp * (dv[id2] - dv[id1]);
				if (m_momentum_conserving)
				{
					if ((node2->m_x - node1->m_x).norm() > SIMD_EPSILON)
					{
						btVector3 dir = (node2->m_x - node1->m_x).normalized();
						local_scaled_df = scaled_k_damp * (dv[id2] - dv[id1]).dot(dir) * dir;
					}
				}
				df[id1] += local_scaled_df;
				df[id2] -= local_scaled_df;
			}
		}
	}

	virtual void buildDampingForceDifferentialDiagonal(btScalar scale, TVStack& diagA)
	{
		// implicit damping force differential
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			btScalar scaled_k_damp = m_dampingStiffness * scale;
			for (int j = 0; j < psb->m_links.size(); ++j)
			{
				const btSoftBody::Link& link = psb->m_links[j];
				btSoftBody::Node* node1 = link.m_n[0];
				btSoftBody::Node* node2 = link.m_n[1];
				size_t id1 = node1->index;
				size_t id2 = node2->index;
				if (m_momentum_conserving)
				{
					if ((node2->m_x - node1->m_x).norm() > SIMD_EPSILON)
					{
						btVector3 dir = (node2->m_x - node1->m_x).normalized();
						for (int d = 0; d < 3; ++d)
						{
							if (node1->m_im > 0)
								diagA[id1][d] -= scaled_k_damp * dir[d] * dir[d];
							if (node2->m_im > 0)
								diagA[id2][d] -= scaled_k_damp * dir[d] * dir[d];
						}
					}
				}
				else
				{
					for (int d = 0; d < 3; ++d)
					{
						if (node1->m_im > 0)
							diagA[id1][d] -= scaled_k_damp;
						if (node2->m_im > 0)
							diagA[id2][d] -= scaled_k_damp;
					}
				}
			}
		}
	}

	virtual double totalElasticEnergy(btScalar dt)
	{
		double energy = 0;
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			const btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_links.size(); ++j)
			{
				const btSoftBody::Link& link = psb->m_links[j];
				btSoftBody::Node* node1 = link.m_n[0];
				btSoftBody::Node* node2 = link.m_n[1];
				btScalar r = link.m_rl;

				// elastic force
				btVector3 dir = (node2->m_q - node1->m_q);
				energy += 0.5 * m_elasticStiffness * (dir.norm() - r) * (dir.norm() - r);
			}
		}
		return energy;
	}

	virtual double totalDampingEnergy(btScalar dt)
	{
		double energy = 0;
		int sz = 0;
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				sz = btMax(sz, psb->m_nodes[j].index);
			}
		}
		TVStack dampingForce;
		dampingForce.resize(sz + 1);
		for (int i = 0; i < dampingForce.size(); ++i)
			dampingForce[i].setZero();
		addScaledDampingForce(0.5, dampingForce);
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				const btSoftBody::Node& node = psb->m_nodes[j];
				energy -= dampingForce[node.index].dot(node.m_v) / dt;
			}
		}
		return energy;
	}

	virtual void addScaledElasticForceDifferential(btScalar scale, const TVStack& dx, TVStack& df)
	{
		// implicit damping force differential
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			const btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_links.size(); ++j)
			{
				const btSoftBody::Link& link = psb->m_links[j];
				btSoftBody::Node* node1 = link.m_n[0];
				btSoftBody::Node* node2 = link.m_n[1];
				size_t id1 = node1->index;
				size_t id2 = node2->index;
				btScalar r = link.m_rl;

				btVector3 dir = (node1->m_q - node2->m_q);
				btScalar dir_norm = dir.norm();
				btVector3 dir_normalized = (dir_norm > SIMD_EPSILON) ? dir.normalized() : btVector3(0, 0, 0);
				btVector3 dx_diff = dx[id1] - dx[id2];
				btVector3 scaled_df = btVector3(0, 0, 0);
				btScalar scaled_k = scale * (link.m_bbending ? m_bendingStiffness : m_elasticStiffness);
				if (dir_norm > SIMD_EPSILON)
				{
					scaled_df -= scaled_k * dir_normalized.dot(dx_diff) * dir_normalized;
					scaled_df += scaled_k * dir_normalized.dot(dx_diff) * ((dir_norm - r) / dir_norm) * dir_normalized;
					scaled_df -= scaled_k * ((dir_norm - r) / dir_norm) * dx_diff;
				}

				df[id1] += scaled_df;
				df[id2] -= scaled_df;
			}
		}
	}

	virtual btDeformableLagrangianForceType getForceType()
	{
		return BT_MASSSPRING_FORCE;
	}
};

#endif /* btMassSpring_h */
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

// Jacobian-vector products and energies are kept identical to the previous
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

#ifndef BT_MOUSE_PICKING_FORCE_H
#define BT_MOUSE_PICKING_FORCE_H

#include "btDeformableLagrangianForce.h"
#include "btSoftBodyInternals.h"

class btDeformableMousePickingForce : public btDeformableLagrangianForce
{
	// ------------------------------------------------------------------ params
	btScalar m_elasticStiffness;  //!< base Hooke constant (N/m)
	btScalar m_dampingStiffness;  //!< base damping constant (N�s/m)

	btScalar m_angleMax;   //!< angle at which the boost saturates (rad)
	btScalar m_scaleCeil;  //!< maximal multiplier (?1)

	const btSoftBody::Face* m_face;
	const btSoftBody::Tetra* m_tetra;
	std::vector<const btSoftBody::Node*> m_nodes;

	btVector3 m_node_to_mouse_x[4], m_node_to_mouse_x_orig[4];
	btVector3 m_node_to_mouse_q[4], m_node_to_mouse_q_orig[4];

	btVector4 m_mouse_bary;

	btTransform m_mouse_transform, m_mouse_transform_orig, m_mouse_transform_orig_inv;
	btScalar m_maxForce;  //!< per-spring clamping
	btScalar m_kRot = 1.0;

	// --------------------------------------------------------- helpers
	const btSoftBody::Node* getNode(int i) const;

	int getIndexCount() const;

	void calculateNodeToMouse();

	void rotateNodeToMouse();

	static inline int PolarDecomposition(const btMatrix3x3& m, btMatrix3x3& q, btMatrix3x3& s)
	{
		static const btPolarDecomposition polar;
		return polar.decompose(m, q, s);
	}

	void orientationError(btScalar& theta, btVector3& axis) const;

	btScalar transBoost(btScalar theta) const;

	// ---------------------------------------------------------------------
	void applyTorqueSpring(const btVector3& T, TVStack& force);

public:
	typedef btAlignedObjectArray<btVector3> TVStack;

	btDeformableMousePickingForce(
		btScalar elasticK,
		btScalar dampingK,
		const btSoftBody::Face* face,
		const btSoftBody::Tetra* tetra,
		std::vector<const btSoftBody::Node*>* nodes,
		const btTransform& mouseTr,
		btScalar maxForce = btScalar(0.3),
		btScalar angleMax = 0.0,  // off by default
		btScalar scaleCeil = 1.0);

	// ------------------------------------------------- Force interface
	void addScaledForces(btScalar s, TVStack& f) override;
	void addScaledExplicitForce(btScalar s, TVStack& f) override;

	// ---------------- elastic ------------------------------------------------
	void addScaledElasticForce(btScalar scale, TVStack& force);

	// --------------- damping -------------------------------------------------
	void addScaledDampingForce(btScalar scale, TVStack& force) override;

	// Jacobian-vector products and energies are kept identical to the previous
	// version except that they use the *base* stiffness (orientation boosting
	// does not affect the linearised system for stability).

	void addScaledDampingForceDifferential(btScalar scale, const TVStack& dv, TVStack& df) override;
	void buildDampingForceDifferentialDiagonal(btScalar, TVStack&) override;

	double totalElasticEnergy(btScalar) override;
	double totalDampingEnergy(btScalar dt) override;

	void addScaledElasticForceDifferential(btScalar scale, const TVStack& dx, TVStack& df) override;

	// ------------------------------------------------ setters ---------------
	void setMousePos(const btVector3& p);
	void setMouseTransform(const btTransform& t);
	void setMaxForce(btScalar f);
	void setElasticStiffness(btScalar k);
	void setDampingStiffness(btScalar k);

	btDeformableLagrangianForceType getForceType() override;
};

#endif /* btMassSpring_h */
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

#include "btDeformableMultiBodyConstraintSolver.h"
#include "BulletReducedDeformableBody/btReducedDeformableBodySolver.h"
#include <iostream>

// override the iterations method to include deformable/multibody contact
btScalar btDeformableMultiBodyConstraintSolver::solveDeformableGroupIterations(btCollisionObject** bodies, int numBodies, btCollisionObject** deformableBodies, int numDeformableBodies, btPersistentManifold** manifoldPtr, int numManifolds, btTypedConstraint** constraints, int numConstraints, const btContactSolverInfo& infoGlobal, btIDebugDraw* debugDrawer)
{
	{
		// pair deformable body with solver body
		pairDeformableAndSolverBody(bodies, numBodies, numDeformableBodies, infoGlobal);

		///this is a special step to resolve penetrations (just for contacts)
		solveGroupCacheFriendlySplitImpulseIterations(bodies, numBodies, deformableBodies, numDeformableBodies, manifoldPtr, numManifolds, constraints, numConstraints, infoGlobal, debugDrawer);

		int maxIterations = m_maxOverrideNumSolverIterations > infoGlobal.m_numIterations ? m_maxOverrideNumSolverIterations : infoGlobal.m_numIterations;
		for (int iteration = 0; iteration < maxIterations; iteration++)
		{
			//fprintf(stderr, "iteration %d\n", iteration);
			// rigid bodies are solved using solver body velocity, but rigid/deformable contact directly uses the velocity of the actual rigid body. So we have to do the following: Solve one iteration of the rigid/rigid contact, get the updated velocity in the solver body and update the velocity of the underlying rigid body. Then solve the rigid/deformable contact. Finally, grab the (once again) updated rigid velocity and update the velocity of the wrapping solver body

			// solve rigid/rigid in solver body
			m_leastSquaresResidual = solveSingleIteration(iteration, bodies, numBodies, manifoldPtr, numManifolds, constraints, numConstraints, infoGlobal, debugDrawer);
			// solver body velocity -> rigid body velocity
			solverBodyWriteBack(infoGlobal);
			btScalar deformableResidual = m_deformableSolver->solveContactConstraints(deformableBodies, numDeformableBodies, infoGlobal);
			// update rigid body velocity in rigid/deformable contact
			m_leastSquaresResidual = btMax(m_leastSquaresResidual, deformableResidual);
			// solver body velocity <- rigid body velocity
			writeToSolverBody(bodies, numBodies, infoGlobal);

			// std::cout << "------------Iteration " << iteration << "------------\n";
			// std::cout << "m_leastSquaresResidual: " << m_leastSquaresResidual << "\n";

			if (m_leastSquaresResidual <= infoGlobal.m_leastSquaresResidualThreshold || (iteration >= (maxIterations - 1)))
			{
#ifdef VERBOSE_RESIDUAL_PRINTF
				if (iteration >= (maxIterations - 1))
					printf("residual = %f at iteration #%d\n", m_leastSquaresResidual, iteration);
#endif
				m_analyticsData.m_numSolverCalls++;
				m_analyticsData.m_numIterationsUsed = iteration + 1;
				m_analyticsData.m_islandId = -2;
				if (numBodies > 0)
					m_analyticsData.m_islandId = bodies[0]->getCompanionId();
				m_analyticsData.m_numBodies = numBodies;
				m_analyticsData.m_numContactManifolds = numManifolds;
				m_analyticsData.m_remainingLeastSquaresResidual = m_leastSquaresResidual;

				m_deformableSolver->deformableBodyInternalWriteBack();
				// std::cout << "[===================Next Step===================]\n";
				break;
			}
		}
	}
	return 0.f;
}

void btDeformableMultiBodyConstraintSolver::solveDeformableBodyGroup(btCollisionObject** bodies, int numBodies, btCollisionObject** deformableBodies, int numDeformableBodies, btPersistentManifold** manifold, int numManifolds, btTypedConstraint** constraints, int numConstraints, btMultiBodyConstraint** multiBodyConstraints, int numMultiBodyConstraints, const btContactSolverInfo& info, btIDebugDraw* debugDrawer, btDispatcher* dispatcher)
{
	m_tmpMultiBodyConstraints = multiBodyConstraints;
	m_tmpNumMultiBodyConstraints = numMultiBodyConstraints;

	// inherited from MultiBodyConstraintSolver
	solveGroupCacheFriendlySetup(bodies, numBodies, manifold, numManifolds, constraints, numConstraints, info, debugDrawer);

	// overriden
	solveDeformableGroupIterations(bodies, numBodies, deformableBodies, numDeformableBodies, manifold, numManifolds, constraints, numConstraints, info, debugDrawer);

	// inherited from MultiBodyConstraintSolver
	solveGroupCacheFriendlyFinish(bodies, numBodies, info);

	m_tmpMultiBodyConstraints = 0;
	m_tmpNumMultiBodyConstraints = 0;
}

void btDeformableMultiBodyConstraintSolver::synchronizeSolverBodyWithRigidBody(btSolverBody* solverBody, btRigidBody* rigidBody)
{
	// Compute the total velocity change
	btVector3 totalDeltaLinearVelocity = rigidBody->getLinearVelocity() - (solverBody->m_linearVelocity + solverBody->m_deltaLinearVelocity);
	btVector3 totalDeltaAngularVelocity = rigidBody->getAngularVelocity() - (solverBody->m_angularVelocity + solverBody->m_deltaAngularVelocity);

	// Update the delta velocities
	solverBody->m_deltaLinearVelocity += totalDeltaLinearVelocity;
	solverBody->m_deltaAngularVelocity += totalDeltaAngularVelocity;

	// Adjust the solver body's base velocities
	solverBody->m_linearVelocity = rigidBody->getLinearVelocity() - solverBody->m_deltaLinearVelocity;
	solverBody->m_angularVelocity = rigidBody->getAngularVelocity() - solverBody->m_deltaAngularVelocity;
}

void btDeformableMultiBodyConstraintSolver::writeToSolverBody(btCollisionObject** bodies, int numBodies, const btContactSolverInfo& infoGlobal)
{
	// reduced soft body solver directly modifies the solver body
	if (m_deformableSolver->isReducedSolver())
	{
		return;
	}

	for (int i = 0; i < numBodies; i++)
	{
		int bodyId = getOrInitSolverBody(*bodies[i], infoGlobal.m_timeStep);

		btRigidBody* body = btRigidBody::upcast(bodies[i]);
		if (body && body->getInvMass())
		{
			btSolverBody& solverBody = m_tmpSolverBodyPool[bodyId];
			synchronizeSolverBodyWithRigidBody(&solverBody, body);
		}
	}
}

void btDeformableMultiBodyConstraintSolver::solverBodyWriteBack(const btContactSolverInfo& infoGlobal)
{
	// reduced soft body solver directly modifies the solver body
	if (m_deformableSolver->isReducedSolver())
	{
		return;
	}

	for (int i = 0; i < m_tmpSolverBodyPool.size(); i++)
	{
		btRigidBody* body = m_tmpSolverBodyPool[i].m_originalBody;
		if (body)
		{
			m_tmpSolverBodyPool[i].m_originalBody->setLinearVelocity(m_tmpSolverBodyPool[i].m_linearVelocity + m_tmpSolverBodyPool[i].m_deltaLinearVelocity);
			m_tmpSolverBodyPool[i].m_originalBody->setAngularVelocity(m_tmpSolverBodyPool[i].m_angularVelocity + m_tmpSolverBodyPool[i].m_deltaAngularVelocity);
		}
	}
}

void btDeformableMultiBodyConstraintSolver::pairDeformableAndSolverBody(btCollisionObject** bodies, int numBodies, int numDeformableBodies, const btContactSolverInfo& infoGlobal)
{
	if (!m_deformableSolver->isReducedSolver())
	{
		return;
	}

	btReducedDeformableBodySolver* solver = static_cast<btReducedDeformableBodySolver*>(m_deformableSolver);

	for (int i = 0; i < numDeformableBodies; ++i)
	{
		for (int k = 0; k < solver->m_nodeRigidConstraints[i].size(); ++k)
		{
			btReducedDeformableNodeRigidContactConstraint& constraint = solver->m_nodeRigidConstraints[i][k];

			if (!constraint.m_contact->m_cti.m_colObj->isStaticObject())
			{
				btCollisionObject& col_obj = const_cast<btCollisionObject&>(*constraint.m_contact->m_cti.m_colObj);

				// object index in the solver body pool
				int bodyId = getOrInitSolverBody(col_obj, infoGlobal.m_timeStep);

				const btRigidBody* body = btRigidBody::upcast(bodies[bodyId]);
				if (body && body->getInvMass())
				{
					// std::cout << "Node: " << constraint.m_node->index << ", body: " << bodyId << "\n";
					btSolverBody& solverBody = m_tmpSolverBodyPool[bodyId];
					constraint.setSolverBody(bodyId, solverBody);
				}
			}
		}

		// for (int j = 0; j < numBodies; j++)
		// {
		// 	int bodyId = getOrInitSolverBody(*bodies[j], infoGlobal.m_timeStep);

		// 	btRigidBody* body = btRigidBody::upcast(bodies[j]);
		// 	if (body && body->getInvMass())
		// 	{
		// 		btSolverBody& solverBody = m_tmpSolverBodyPool[bodyId];
		// 		m_deformableSolver->pairConstraintWithSolverBody(i, bodyId, solverBody);
		// 	}
		// }
	}
}

void btDeformableMultiBodyConstraintSolver::solveGroupCacheFriendlySplitImpulseIterations(btCollisionObject** bodies, int numBodies, btCollisionObject** deformableBodies, int numDeformableBodies, btPersistentManifold** manifoldPtr, int numManifolds, btTypedConstraint** constraints, int numConstraints, const btContactSolverInfo& infoGlobal, btIDebugDraw* debugDrawer)
{
	BT_PROFILE("solveGroupCacheFriendlySplitImpulseIterations");
	int iteration;
	if (infoGlobal.m_splitImpulse)
	{
		{
			for (iteration = 0; iteration < infoGlobal.m_numIterations; iteration++)
			{
				fprintf(stderr, "iteration %d\n", iteration);
				btScalar leastSquaresResidual = 0.f;
				{
					int numPoolConstraints = m_tmpSolverContactConstraintPool.size();
					int j;
					for (j = 0; j < numPoolConstraints; j++)
					{
						const btSolverConstraint& solveManifold = m_tmpSolverContactConstraintPool[m_orderTmpConstraintPool[j]];

						btScalar residual = resolveSplitPenetrationImpulse(m_tmpSolverBodyPool[solveManifold.m_solverBodyIdA], m_tmpSolverBodyPool[solveManifold.m_solverBodyIdB], solveManifold);
						leastSquaresResidual = btMax(leastSquaresResidual, residual * residual);
					}
					// solve the position correction between deformable and rigid/multibody
					//                    btScalar residual = m_deformableSolver->solveSplitImpulse(infoGlobal);
					btScalar residual = m_deformableSolver->m_objective->m_projection.solveSplitImpulse(deformableBodies, numDeformableBodies, infoGlobal);
					leastSquaresResidual = btMax(leastSquaresResidual, residual * residual);
				}
				if (leastSquaresResidual <= infoGlobal.m_leastSquaresResidualThreshold || iteration >= (infoGlobal.m_numIterations - 1))
				{
#ifdef VERBOSE_RESIDUAL_PRINTF
					if (iteration >= (infoGlobal.m_numIterations - 1))
						printf("split impulse residual = %f at iteration #%d\n", leastSquaresResidual, iteration);
#endif
					break;
				}
			}
		}
	}
}
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

#ifndef BT_DEFORMABLE_MULTIBODY_CONSTRAINT_SOLVER_H
#define BT_DEFORMABLE_MULTIBODY_CONSTRAINT_SOLVER_H

#include "btDeformableBodySolver.h"
#include "BulletDynamics/Featherstone/btMultiBodyConstraintSolver.h"

class btDeformableBodySolver;

// btDeformableMultiBodyConstraintSolver extendsn btMultiBodyConstraintSolver to solve for the contact among rigid/multibody and deformable bodies. Notice that the following constraints
// 1. rigid/multibody against rigid/multibody
// 2. rigid/multibody against deforamble
// 3. deformable against deformable
// 4. deformable self collision
// 5. joint constraints
// are all coupled in this solve.
ATTRIBUTE_ALIGNED16(class)
btDeformableMultiBodyConstraintSolver : public btMultiBodyConstraintSolver
{
	btDeformableBodySolver* m_deformableSolver;

protected:
	// override the iterations method to include deformable/multibody contact
	//    virtual btScalar solveGroupCacheFriendlyIterations(btCollisionObject** bodies,int numBodies,btPersistentManifold** manifoldPtr, int numManifolds,btTypedConstraint** constraints,int numConstraints,const btContactSolverInfo& infoGlobal,btIDebugDraw* debugDrawer);

	// write the velocity of the the solver body to the underlying rigid body
	void solverBodyWriteBack(const btContactSolverInfo& infoGlobal);

	void synchronizeSolverBodyWithRigidBody(btSolverBody * solverBody, btRigidBody * rigidBody);

	// write the velocity of the underlying rigid body to the the the solver body
	void writeToSolverBody(btCollisionObject * *bodies, int numBodies, const btContactSolverInfo& infoGlobal);

	// let each deformable body knows which solver body is in constact
	void pairDeformableAndSolverBody(btCollisionObject * *bodies, int numBodies, int numDeformableBodies, const btContactSolverInfo& infoGlobal);

	virtual void solveGroupCacheFriendlySplitImpulseIterations(btCollisionObject * *bodies, int numBodies, btCollisionObject** deformableBodies, int numDeformableBodies, btPersistentManifold** manifoldPtr, int numManifolds, btTypedConstraint** constraints, int numConstraints, const btContactSolverInfo& infoGlobal, btIDebugDraw* debugDrawer);

	virtual btScalar solveDeformableGroupIterations(btCollisionObject * *bodies, int numBodies, btCollisionObject** deformableBodies, int numDeformableBodies, btPersistentManifold** manifoldPtr, int numManifolds, btTypedConstraint** constraints, int numConstraints, const btContactSolverInfo& infoGlobal, btIDebugDraw* debugDrawer);

public:
	BT_DECLARE_ALIGNED_ALLOCATOR();

	void setDeformableSolver(btDeformableBodySolver * deformableSolver)
	{
		m_deformableSolver = deformableSolver;
	}

	virtual void solveDeformableBodyGroup(btCollisionObject * *bodies, int numBodies, btCollisionObject** deformableBodies, int numDeformableBodies, btPersistentManifold** manifold, int numManifolds, btTypedConstraint** constraints, int numConstraints, btMultiBodyConstraint** multiBodyConstraints, int numMultiBodyConstraints, const btContactSolverInfo& info, btIDebugDraw* debugDrawer, btDispatcher* dispatcher);
};

#endif /* BT_DEFORMABLE_MULTIBODY_CONSTRAINT_SOLVER_H */
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

/* ====== Overview of the Deformable Algorithm ====== */

/*
A single step of the deformable body simulation contains the following main components:
Call internalStepSimulation multiple times, to achieve 240Hz (4 steps of 60Hz).
1. Deformable maintaintenance of rest lengths and volume preservation. Forces only depend on position: Update velocity to a temporary state v_{n+1}^* = v_n + explicit_force * dt / mass, where explicit forces include gravity and elastic forces.
2. Detect discrete collisions between rigid and deformable bodies at position x_{n+1}^* = x_n + dt * v_{n+1}^*.

3a. Solve all constraints, including LCP. Contact, position correction due to numerical drift, friction, and anchors for deformable.

3b. 5 Newton steps (multiple step). Conjugent Gradient solves linear system. Deformable Damping: Then velocities of deformable bodies v_{n+1} are solved in
        M(v_{n+1} - v_{n+1}^*) = damping_force * dt / mass,
   by a conjugate gradient solver, where the damping force is implicit and depends on v_{n+1}.
   Make sure contact constraints are not violated in step b by performing velocity projections as in the paper by Baraff and Witkin https://www.cs.cmu.edu/~baraff/papers/sig98.pdf. Dynamic frictions are treated as a force and added to the rhs of the CG solve, whereas static frictions are treated as constraints similar to contact.
4. Position is updated via x_{n+1} = x_n + dt * v_{n+1}.


The algorithm also closely resembles the one in http://physbam.stanford.edu/~fedkiw/papers/stanford2008-03.pdf
 */

#include <stdio.h>
#include "btDeformableMultiBodyDynamicsWorld.h"
#include "DeformableBodyInplaceSolverIslandCallback.h"
#include "btDeformableBodySolver.h"
#include "LinearMath/btQuickprof.h"
#include "btSoftBodyInternals.h"

btDeformableMultiBodyDynamicsWorld::btDeformableMultiBodyDynamicsWorld(btDispatcher* dispatcher, btBroadphaseInterface* pairCache, btDeformableMultiBodyConstraintSolver* constraintSolver, btCollisionConfiguration* collisionConfiguration, btDeformableBodySolver* deformableBodySolver)
	: btMultiBodyDynamicsWorld(dispatcher, pairCache, (btMultiBodyConstraintSolver*)constraintSolver, collisionConfiguration),
	  m_deformableBodySolver(deformableBodySolver),
	  m_solverCallback(0)
{
	m_drawFlags = fDrawFlags::Std;
	m_drawNodeTree = true;
	m_drawFaceTree = false;
	m_drawClusterTree = false;
	m_sbi.m_broadphase = pairCache;
	m_sbi.m_dispatcher = dispatcher;
	m_sbi.m_sparsesdf.Initialize();
	m_sbi.m_sparsesdf.setDefaultVoxelsz(0.005);
	m_sbi.m_sparsesdf.Reset();

	m_sbi.air_density = (btScalar)1.2;
	m_sbi.water_density = 0;
	m_sbi.water_offset = 0;
	m_sbi.water_normal = btVector3(0, 0, 0);
	m_sbi.m_gravity.setValue(0, -9.8, 0);
	m_internalTime = 0.0;
	m_implicit = false;
	m_lineSearch = false;
	m_useProjection = false;
	m_ccdIterations = 5;
	m_solverDeformableBodyIslandCallback = new DeformableBodyInplaceSolverIslandCallback(constraintSolver, dispatcher);
}

btDeformableMultiBodyDynamicsWorld::~btDeformableMultiBodyDynamicsWorld()
{
	delete m_solverDeformableBodyIslandCallback;
}

void btDeformableMultiBodyDynamicsWorld::performDiscreteCollisionDetection()
{
	BT_PROFILE("performDiscreteCollisionDetection");

	btDispatcherInfo& dispatchInfo = getDispatchInfo();

	updateAabbs();

	computeOverlappingPairs();

	addSoftsWithSelfCollisionCheckToOverlappingPairs();

	btDispatcher* dispatcher = getDispatcher();
	{
		BT_PROFILE("dispatchAllCollisionPairs");
		if (dispatcher)
			dispatcher->dispatchAllCollisionPairs(m_broadphasePairCache->getOverlappingPairCache(), dispatchInfo, m_dispatcher1);
	}
}

void btDeformableMultiBodyDynamicsWorld::addSoftsWithSelfCollisionCheckToOverlappingPairs()
{
	btBroadphaseInterface* broadphase = getBroadphase();
	btOverlappingPairCache* pairCache = broadphase->getOverlappingPairCache();

	if (!pairCache)
		return;

	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		auto& soft = m_softBodies[i];
		if (soft->getCollisionShape()->getShapeType() != SOFTBODY_SHAPE_PROXYTYPE && (soft->m_cfg.collisions & btSoftBody::fCollision::CL_SELF))  // Do this only for "our" type of softs
		{
			if (!soft->getBroadphaseHandle())
				continue;

			pairCache->addOverlappingPair(soft->getBroadphaseHandle(), soft->getBroadphaseHandle());
		}
	}
}

void btDeformableMultiBodyDynamicsWorld::internalSingleStepSimulation(btScalar timeStep)
{
	BT_PROFILE("internalSingleStepSimulation");

	if (0 != m_internalPreTickCallback)
	{
		(*m_internalPreTickCallback)(this, timeStep);
	}
	reinitialize(timeStep);

	// add gravity to velocity of rigid and multi bodys
	applyRigidBodyGravity(timeStep);

	///apply gravity and explicit force to velocity, predict motion
	predictUnconstraintMotion(timeStep);

#ifdef BT_SAFE_UPDATE_DEBUG
	fprintf(stderr, "framestart()\n");
#endif
	///perform collision detection that involves rigid/multi bodies
	performDiscreteCollisionDetection();

	if (0 != m_internalPostDiscreteCollisionDetectionTickCallback)
	{
		(*m_internalPostDiscreteCollisionDetectionTickCallback)(this, timeStep);
	}

	btMultiBodyDynamicsWorld::calculateSimulationIslands();

	updateLastSafeTransforms();

#ifdef BT_SAFE_UPDATE_DEBUG
	fprintf(stderr, "frameend()\n");
	fprintf(stderr, "framestart()\n");
	btCollisionObject::gDebug = true;
	performDiscreteCollisionDetection();
	fprintf(stderr, "drawpoint \"VERIFY\" [0,0,0][1,1,1,1]\n");
	btCollisionObject::gDebug = false;

	fprintf(stderr, "frameend()\n");
#endif

	beforeSolverCallbacks(timeStep);

	// ///solve contact constraints and then deformable bodies momemtum equation
	solveConstraints(timeStep);

	afterSolverCallbacks(timeStep);

	performDeformableCollisionDetection();

	applyRepulsionForce(timeStep);

	performGeometricCollisions(timeStep);

	integrateTransforms(timeStep);

	///update vehicle simulation
	btMultiBodyDynamicsWorld::updateActions(timeStep);

	updateActivationState(timeStep);

	if (0 != m_internalTickCallback)
	{
		(*m_internalTickCallback)(this, timeStep);
	}

	// End solver-wise simulation step
	// ///////////////////////////////

	++m_executed_step_counter;
	m_dispatchInfo.m_stepCounter = m_executed_step_counter;
}

void btDeformableMultiBodyDynamicsWorld::performDeformableCollisionDetection()
{
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		m_softBodies[i]->m_softSoftCollision = true;
	}

	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		for (int j = i; j < m_softBodies.size(); ++j)
		{
			if (m_softBodies[i]->getCollisionShape()->getShapeType() != SOFTBODY_SHAPE_PROXYTYPE || m_softBodies[j]->getCollisionShape()->getShapeType() != SOFTBODY_SHAPE_PROXYTYPE)
			{
				// If any shape is not the default soft shape, then the collision is checked elsewhere
				continue;
			}
			m_softBodies[i]->defaultCollisionHandler(m_softBodies[j]);
		}
	}

	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		m_softBodies[i]->m_softSoftCollision = false;
	}
}

void btDeformableMultiBodyDynamicsWorld::updateActivationState(btScalar timeStep)
{
	for (int i = 0; i < m_softBodies.size(); i++)
	{
		btSoftBody* psb = m_softBodies[i];
		psb->updateDeactivation(timeStep);
		if (psb->wantsSleeping())
		{
			if (psb->getActivationState() == ACTIVE_TAG)
				psb->setActivationState(WANTS_DEACTIVATION);
			if (psb->getActivationState() == ISLAND_SLEEPING)
			{
				psb->setZeroVelocity();
			}
		}
		else
		{
			if (psb->getActivationState() != DISABLE_DEACTIVATION)
				psb->setActivationState(ACTIVE_TAG);
		}
	}
	btMultiBodyDynamicsWorld::updateActivationState(timeStep);
}

void btDeformableMultiBodyDynamicsWorld::applyRepulsionForce(btScalar timeStep)
{
	BT_PROFILE("btDeformableMultiBodyDynamicsWorld::applyRepulsionForce");
	for (int i = 0; i < m_softBodies.size(); i++)
	{
		btSoftBody* psb = m_softBodies[i];
		if (psb->isActive() && !psb->isStaticObject())
		{
			psb->applyRepulsionForce(timeStep, true);
		}
	}
}

void btDeformableMultiBodyDynamicsWorld::performGeometricCollisions(btScalar timeStep)
{
	BT_PROFILE("btDeformableMultiBodyDynamicsWorld::performGeometricCollisions");
	// refit the BVH tree for CCD
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btSoftBody* psb = m_softBodies[i];
		if (psb->isActive() && !psb->isStaticObject())
		{
			m_softBodies[i]->updateFaceTree(true, false);
			m_softBodies[i]->updateNodeTree(true, false);
			for (int j = 0; j < m_softBodies[i]->m_faces.size(); ++j)
			{
				btSoftBody::Face& f = m_softBodies[i]->m_faces[j];
				f.m_n0 = (f.m_n[1]->m_x - f.m_n[0]->m_x).cross(f.m_n[2]->m_x - f.m_n[0]->m_x);
			}
		}
	}

	// clear contact points & update DBVT
	for (int r = 0; r < m_ccdIterations; ++r)
	{
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (psb->isActive() && !psb->isStaticObject())
			{
				// clear contact points in the previous iteration
				psb->m_faceNodeContactsCCD.clear();

				// update m_q and normals for CCD calculation
				for (int j = 0; j < psb->m_nodes.size(); ++j)
				{
					psb->m_nodes[j].m_q = psb->m_nodes[j].m_x + timeStep * psb->m_nodes[j].m_v;
				}
				for (int j = 0; j < psb->m_faces.size(); ++j)
				{
					btSoftBody::Face& f = psb->m_faces[j];
					f.m_n1 = (f.m_n[1]->m_q - f.m_n[0]->m_q).cross(f.m_n[2]->m_q - f.m_n[0]->m_q);
					f.m_vn = (f.m_n[1]->m_v - f.m_n[0]->m_v).cross(f.m_n[2]->m_v - f.m_n[0]->m_v) * timeStep * timeStep;
				}
			}
		}

		// apply CCD to register new contact points
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			for (int j = i; j < m_softBodies.size(); ++j)
			{
				btSoftBody* psb1 = m_softBodies[i];
				btSoftBody* psb2 = m_softBodies[j];
				if (psb1->isActive() && !psb1->isStaticObject() && psb2->isActive() && !psb2->isStaticObject())
				{
					if (m_softBodies[i]->getCollisionShape()->getShapeType() != SOFTBODY_SHAPE_PROXYTYPE || m_softBodies[j]->getCollisionShape()->getShapeType() != SOFTBODY_SHAPE_PROXYTYPE)
						continue;
					m_softBodies[i]->geometricCollisionHandler(m_softBodies[j]);
				}
			}
		}

		int penetration_count = 0;
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (psb->isActive() && !psb->isStaticObject())
			{
				penetration_count += psb->m_faceNodeContactsCCD.size();
				;
			}
		}
		if (penetration_count == 0)
		{
			break;
		}

		// apply inelastic impulse
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (psb->isActive() && !psb->isStaticObject())
			{
				psb->applyRepulsionForce(timeStep, false);
			}
		}
	}
}

void btDeformableMultiBodyDynamicsWorld::softBodySelfCollision()
{
	BT_PROFILE("btDeformableMultiBodyDynamicsWorld::softBodySelfCollision");
	for (int i = 0; i < m_softBodies.size(); i++)
	{
		btSoftBody* psb = m_softBodies[i];
		if (psb->isActive() && !psb->isStaticObject())
		{
			psb->defaultCollisionHandler(psb);
		}
	}
}

void btDeformableMultiBodyDynamicsWorld::positionCorrection(btScalar timeStep)
{
	// correct the position of rigid bodies with temporary velocity generated from split impulse
	btContactSolverInfo infoGlobal;
	btVector3 zero(0, 0, 0);
	for (int i = 0; i < m_nonStaticRigidBodies.size(); ++i)
	{
		btRigidBody* rb = m_nonStaticRigidBodies[i];
		//correct the position/orientation based on push/turn recovery
		btTransform newTransform;
		btVector3 pushVelocity = rb->getPushVelocity();
		btVector3 turnVelocity = rb->getTurnVelocity();
		if (pushVelocity[0] != 0.f || pushVelocity[1] != 0 || pushVelocity[2] != 0 || turnVelocity[0] != 0.f || turnVelocity[1] != 0 || turnVelocity[2] != 0)
		{
			btTransformUtil::integrateTransform(rb->getWorldTransform(), pushVelocity, turnVelocity * infoGlobal.m_splitImpulseTurnErp, timeStep, newTransform);
			fprintf(stderr, "btDeformableMultiBodyDynamicsWorld::positionCorrection from %f %f %f to %f %f %f pushVelocity %f %f %f\n",
					rb->getWorldTransform().getOrigin().x(), rb->getWorldTransform().getOrigin().y(), rb->getWorldTransform().getOrigin().z(),
					newTransform.getOrigin().x(), newTransform.getOrigin().y(), newTransform.getOrigin().z(),
					pushVelocity.x(), pushVelocity.y(), pushVelocity.z());
			rb->setWorldTransform(newTransform);
			rb->setPushVelocity(zero);
			rb->setTurnVelocity(zero);
		}
	}
}

void btDeformableMultiBodyDynamicsWorld::integrateTransforms(btScalar timeStep)
{
	BT_PROFILE("integrateTransforms");
	positionCorrection(timeStep);
	btMultiBodyDynamicsWorld::integrateTransforms(timeStep);
	m_deformableBodySolver->applyTransforms(timeStep);
}

void btDeformableMultiBodyDynamicsWorld::solveConstraints(btScalar timeStep)
{
	BT_PROFILE("btDeformableMultiBodyDynamicsWorld::solveConstraints");
	// save v_{n+1}^* velocity after explicit forces
	m_deformableBodySolver->backupVelocity();

	// set up constraints among multibodies and between multibodies and deformable bodies
	setupConstraints();

	// solve contact constraints
	solveContactConstraints();

	// set up the directions in which the velocity does not change in the momentum solve
	if (m_useProjection)
		m_deformableBodySolver->setProjection();
	else
		m_deformableBodySolver->setLagrangeMultiplier();

	// for explicit scheme, m_backupVelocity = v_{n+1}^*
	// for implicit scheme, m_backupVelocity = v_n
	// Here, set dv = v_{n+1} - v_n for nodes in contact
	m_deformableBodySolver->setupDeformableSolve(m_implicit);

	// At this point, dv should be golden for nodes in contact
	// proceed to solve deformable momentum equation
	m_deformableBodySolver->solveDeformableConstraints(timeStep);
}

void btDeformableMultiBodyDynamicsWorld::setupConstraints()
{
	// set up constraints between multibody and deformable bodies
	m_deformableBodySolver->setConstraints(m_solverInfo);

	// set up constraints among multibodies
	{
		sortConstraints();
		// setup the solver callback
		btMultiBodyConstraint** sortedMultiBodyConstraints = m_sortedMultiBodyConstraints.size() ? &m_sortedMultiBodyConstraints[0] : 0;
		btTypedConstraint** constraintsPtr = getNumConstraints() ? &m_sortedConstraints[0] : 0;
		m_solverDeformableBodyIslandCallback->setup(&m_solverInfo, constraintsPtr, m_sortedConstraints.size(), sortedMultiBodyConstraints, m_sortedMultiBodyConstraints.size(), getDebugDrawer());

		// build islands
		m_islandManager->buildIslands(getCollisionWorld()->getDispatcher(), getCollisionWorld());
	}
}

void btDeformableMultiBodyDynamicsWorld::sortConstraints()
{
	m_sortedConstraints.resize(m_constraints.size());
	int i;
	for (i = 0; i < getNumConstraints(); i++)
	{
		m_sortedConstraints[i] = m_constraints[i];
	}
	m_sortedConstraints.quickSort(btSortConstraintOnIslandPredicate2());

	m_sortedMultiBodyConstraints.resize(m_multiBodyConstraints.size());
	for (i = 0; i < m_multiBodyConstraints.size(); i++)
	{
		m_sortedMultiBodyConstraints[i] = m_multiBodyConstraints[i];
	}
	m_sortedMultiBodyConstraints.quickSort(btSortMultiBodyConstraintOnIslandPredicate());
}

void btDeformableMultiBodyDynamicsWorld::solveContactConstraints()
{
	// process constraints on each island
	m_islandManager->processIslands(getCollisionWorld()->getDispatcher(), getCollisionWorld(), m_solverDeformableBodyIslandCallback);

	// process deferred
	m_solverDeformableBodyIslandCallback->processConstraints();
	m_constraintSolver->allSolved(m_solverInfo, m_debugDrawer);

	// write joint feedback
	{
		for (int i = 0; i < this->m_multiBodies.size(); i++)
		{
			btMultiBody* bod = m_multiBodies[i];

			bool isSleeping = false;

			if (bod->getBaseCollider() && bod->getBaseCollider()->getActivationState() == ISLAND_SLEEPING)
			{
				isSleeping = true;
			}
			for (int b = 0; b < bod->getNumLinks(); b++)
			{
				if (bod->getLink(b).m_collider && bod->getLink(b).m_collider->getActivationState() == ISLAND_SLEEPING)
					isSleeping = true;
			}

			if (!isSleeping)
			{
				//useless? they get resized in stepVelocities once again (AND DIFFERENTLY)
				m_scratch_r.resize(bod->getNumLinks() + 1);  //multidof? ("Y"s use it and it is used to store qdd)
				m_scratch_v.resize(bod->getNumLinks() + 1);
				m_scratch_m.resize(bod->getNumLinks() + 1);

				if (bod->internalNeedsJointFeedback())
				{
					if (!bod->isUsingRK4Integration())
					{
						if (bod->internalNeedsJointFeedback())
						{
							bool isConstraintPass = true;
							bod->computeAccelerationsArticulatedBodyAlgorithmMultiDof(m_solverInfo.m_timeStep, m_scratch_r, m_scratch_v, m_scratch_m, isConstraintPass,
																					  getSolverInfo().m_jointFeedbackInWorldSpace,
																					  getSolverInfo().m_jointFeedbackInJointFrame);
						}
					}
				}
			}
		}
	}

	for (int i = 0; i < this->m_multiBodies.size(); i++)
	{
		btMultiBody* bod = m_multiBodies[i];
		bod->processDeltaVeeMultiDof2();
	}
}

void btDeformableMultiBodyDynamicsWorld::addSoftBody(btSoftBody* body, int collisionFilterGroup, int collisionFilterMask)
{
	m_softBodies.push_back(body);

	// Set the soft body solver that will deal with this body
	// to be the world's solver
	body->setSoftBodySolver(m_deformableBodySolver);

	btCollisionWorld::addCollisionObject(body,
										 collisionFilterGroup,
										 collisionFilterMask);
}

void btDeformableMultiBodyDynamicsWorld::predictUnconstraintMotion(btScalar timeStep)
{
	BT_PROFILE("predictUnconstraintMotion");
	btMultiBodyDynamicsWorld::predictUnconstraintMotion(timeStep);
	m_deformableBodySolver->predictMotion(timeStep);
}

void btDeformableMultiBodyDynamicsWorld::setGravity(const btVector3& gravity)
{
	btDiscreteDynamicsWorld::setGravity(gravity);
	m_deformableBodySolver->setGravity(gravity);
}

void btDeformableMultiBodyDynamicsWorld::reinitialize(btScalar timeStep)
{
	m_internalTime += timeStep;
	m_deformableBodySolver->setImplicit(m_implicit);
	m_deformableBodySolver->setLineSearch(m_lineSearch);
	m_deformableBodySolver->reinitialize(m_softBodies, timeStep);
	btDispatcherInfo& dispatchInfo = btMultiBodyDynamicsWorld::getDispatchInfo();
	dispatchInfo.m_timeStep = timeStep;
	dispatchInfo.m_debugDraw = btMultiBodyDynamicsWorld::getDebugDrawer();
	btMultiBodyDynamicsWorld::getSolverInfo().m_timeStep = timeStep;
	if (m_useProjection)
	{
		m_deformableBodySolver->m_useProjection = true;
		m_deformableBodySolver->setStrainLimiting(true);
		m_deformableBodySolver->setPreconditioner(btDeformableBackwardEulerObjective::Mass_preconditioner);
	}
	else
	{
		m_deformableBodySolver->m_useProjection = false;
		m_deformableBodySolver->setStrainLimiting(false);
		m_deformableBodySolver->setPreconditioner(btDeformableBackwardEulerObjective::KKT_preconditioner);
	}
}

void btDeformableMultiBodyDynamicsWorld::debugDrawWorld()
{
	btMultiBodyDynamicsWorld::debugDrawWorld();

	for (int i = 0; i < getSoftBodyArray().size(); i++)
	{
		btSoftBody* psb = (btSoftBody*)getSoftBodyArray()[i];
		{
			btSoftBodyHelpers::DrawFrame(psb, getDebugDrawer());
			btSoftBodyHelpers::Draw(psb, getDebugDrawer(), getDrawFlags());
		}
	}
}

void btDeformableMultiBodyDynamicsWorld::applyRigidBodyGravity(btScalar timeStep)
{
	// Gravity is applied in stepSimulation and then cleared here and then applied here and then cleared here again
	// so that 1) gravity is applied to velocity before constraint solve and 2) gravity is applied in each substep
	// when there are multiple substeps
	btMultiBodyDynamicsWorld::applyGravity();
	// integrate rigid body gravity
	for (int i = 0; i < m_nonStaticRigidBodies.size(); ++i)
	{
		btRigidBody* rb = m_nonStaticRigidBodies[i];
		rb->integrateVelocities(timeStep);
	}

	// integrate multibody gravity
	{
		forwardKinematics();
		clearMultiBodyConstraintForces();
		{
			for (int i = 0; i < this->m_multiBodies.size(); i++)
			{
				btMultiBody* bod = m_multiBodies[i];

				bool isSleeping = false;

				if (bod->getBaseCollider() && bod->getBaseCollider()->getActivationState() == ISLAND_SLEEPING)
				{
					isSleeping = true;
				}
				for (int b = 0; b < bod->getNumLinks(); b++)
				{
					if (bod->getLink(b).m_collider && bod->getLink(b).m_collider->getActivationState() == ISLAND_SLEEPING)
						isSleeping = true;
				}

				if (!isSleeping)
				{
					m_scratch_r.resize(bod->getNumLinks() + 1);
					m_scratch_v.resize(bod->getNumLinks() + 1);
					m_scratch_m.resize(bod->getNumLinks() + 1);
					bool isConstraintPass = false;
					{
						if (!bod->isUsingRK4Integration())
						{
							bod->computeAccelerationsArticulatedBodyAlgorithmMultiDof(m_solverInfo.m_timeStep,
																					  m_scratch_r, m_scratch_v, m_scratch_m, isConstraintPass,
																					  getSolverInfo().m_jointFeedbackInWorldSpace,
																					  getSolverInfo().m_jointFeedbackInJointFrame);
						}
						else
						{
							btAssert(" RK4Integration is not supported");
						}
					}
				}
			}
		}
	}
	clearGravity();
}

void btDeformableMultiBodyDynamicsWorld::clearGravity()
{
	BT_PROFILE("btMultiBody clearGravity");
	// clear rigid body gravity
	for (int i = 0; i < m_nonStaticRigidBodies.size(); i++)
	{
		btRigidBody* body = m_nonStaticRigidBodies[i];
		if (body->isActive() && !body->isStaticObject())
		{
			body->clearGravity();
		}
	}
	// clear multibody gravity
	for (int i = 0; i < this->m_multiBodies.size(); i++)
	{
		btMultiBody* bod = m_multiBodies[i];

		bool isSleeping = false;

		if (bod->getBaseCollider() && bod->getBaseCollider()->getActivationState() == ISLAND_SLEEPING)
		{
			isSleeping = true;
		}
		for (int b = 0; b < bod->getNumLinks(); b++)
		{
			if (bod->getLink(b).m_collider && bod->getLink(b).m_collider->getActivationState() == ISLAND_SLEEPING)
				isSleeping = true;
		}

		if (!isSleeping)
		{
			bod->addBaseForce(-m_gravity * bod->getBaseMass());

			for (int j = 0; j < bod->getNumLinks(); ++j)
			{
				bod->addLinkForce(j, -m_gravity * bod->getLinkMass(j));
			}
		}
	}
}

void btDeformableMultiBodyDynamicsWorld::beforeSolverCallbacks(btScalar timeStep)
{
	if (0 != m_solverCallback)
	{
		(*m_solverCallback)(m_internalTime, this);
	}
}

void btDeformableMultiBodyDynamicsWorld::afterSolverCallbacks(btScalar timeStep)
{
	if (0 != m_solverCallback)
	{
		(*m_solverCallback)(m_internalTime, this);
	}
}

void btDeformableMultiBodyDynamicsWorld::addForce(btSoftBody* psb, btDeformableLagrangianForce* force)
{
	btAlignedObjectArray<btDeformableLagrangianForce*>& forces = *m_deformableBodySolver->getLagrangianForceArray();
	bool added = false;
	for (int i = 0; i < forces.size(); ++i)
	{
		if (forces[i]->getForceType() == force->getForceType())
		{
			forces[i]->addSoftBody(psb);
			added = true;
			break;
		}
	}
	if (!added)
	{
		force->addSoftBody(psb);
		force->setIndices(m_deformableBodySolver->getIndices());
		forces.push_back(force);
	}
}

void btDeformableMultiBodyDynamicsWorld::addForce(btDeformableLagrangianForce* force)
{
	btAlignedObjectArray<btDeformableLagrangianForce*>& forces = *m_deformableBodySolver->getLagrangianForceArray();
	force->setIndices(m_deformableBodySolver->getIndices());
	forces.push_back(force);
}

void btDeformableMultiBodyDynamicsWorld::removeForce(btSoftBody* psb, btDeformableLagrangianForce* force)
{
	btAlignedObjectArray<btDeformableLagrangianForce*>& forces = *m_deformableBodySolver->getLagrangianForceArray();
	int removed_index = -1;
	for (int i = 0; i < forces.size(); ++i)
	{
		if (forces[i]->getForceType() == force->getForceType())
		{
			forces[i]->removeSoftBody(psb);
			if (forces[i]->m_softBodies.size() == 0)
				removed_index = i;
			break;
		}
	}
	if (removed_index >= 0)
		forces.removeAtIndex(removed_index);
}

void btDeformableMultiBodyDynamicsWorld::removeForce(btDeformableLagrangianForce* force)
{
	btAlignedObjectArray<btDeformableLagrangianForce*>& forces = *m_deformableBodySolver->getLagrangianForceArray();
	int removed_index = -1;
	for (int i = 0; i < forces.size(); ++i)
	{
		if (forces[i] == force)
		{
			removed_index = i;
			break;
		}
	}
	if (removed_index >= 0)
		forces.removeAtIndex(removed_index);
}

void btDeformableMultiBodyDynamicsWorld::removeSoftBodyForce(btSoftBody* psb)
{
	btAlignedObjectArray<btDeformableLagrangianForce*>& forces = *m_deformableBodySolver->getLagrangianForceArray();
	for (int i = 0; i < forces.size(); ++i)
	{
		forces[i]->removeSoftBody(psb);
	}
}

void btDeformableMultiBodyDynamicsWorld::removeSoftBody(btSoftBody* body)
{
	removeSoftBodyForce(body);
	m_softBodies.remove(body);
	btCollisionWorld::removeCollisionObject(body);
	// force a reinitialize so that node indices get updated.
	m_deformableBodySolver->reinitialize(m_softBodies, btScalar(-1));
}

void btDeformableMultiBodyDynamicsWorld::removeCollisionObject(btCollisionObject* collisionObject)
{
	btSoftBody* body = btSoftBody::upcast(collisionObject);
	if (body)
		removeSoftBody(body);
	else
		btDiscreteDynamicsWorld::removeCollisionObject(collisionObject);
}

int btDeformableMultiBodyDynamicsWorld::stepSimulation(btScalar timeStep, int maxSubSteps, btScalar fixedTimeStep)
{
	startProfiling(timeStep);

	int numSimulationSubSteps = 0;

	if (maxSubSteps)
	{
		//fixed timestep with interpolation
		m_fixedTimeStep = fixedTimeStep;
		m_localTime += timeStep;
		if (m_localTime >= fixedTimeStep)
		{
			numSimulationSubSteps = int(m_localTime / fixedTimeStep);
			m_localTime -= numSimulationSubSteps * fixedTimeStep;
		}
	}
	else
	{
		//variable timestep
		fixedTimeStep = timeStep;
		m_localTime = m_latencyMotionStateInterpolation ? 0 : timeStep;
		m_fixedTimeStep = 0;
		if (btFuzzyZero(timeStep))
		{
			numSimulationSubSteps = 0;
			maxSubSteps = 0;
		}
		else
		{
			numSimulationSubSteps = 1;
			maxSubSteps = 1;
		}
	}

	//process some debugging flags
	if (getDebugDrawer())
	{
		btIDebugDraw* debugDrawer = getDebugDrawer();
		gDisableDeactivation = (debugDrawer->getDebugMode() & btIDebugDraw::DBG_NoDeactivation) != 0;
	}
	if (numSimulationSubSteps)
	{
		//clamp the number of substeps, to prevent simulation grinding spiralling down to a halt
		int clampedSimulationSteps = (numSimulationSubSteps > maxSubSteps) ? maxSubSteps : numSimulationSubSteps;

		saveKinematicState(fixedTimeStep * clampedSimulationSteps);

		for (int i = 0; i < clampedSimulationSteps; i++)
		{
			internalSingleStepSimulation(fixedTimeStep);
			synchronizeMotionStates();
		}
	}
	else
	{
		synchronizeMotionStates();
	}

	clearForces();

#ifndef BT_NO_PROFILE
	CProfileManager::Increment_Frame_Counter();
#endif  //BT_NO_PROFILE

	return numSimulationSubSteps;
}

void btDeformableMultiBodyDynamicsWorld::updateLastSafeTransforms()
{
	BT_PROFILE("updateLastSafeTransforms");

	processLastSafeTransforms(m_nonStaticRigidBodies.size() == 0 ? nullptr : reinterpret_cast<btCollisionObject**>(&m_nonStaticRigidBodies[0]), m_nonStaticRigidBodies.size(),
							  m_softBodies.size() == 0 ? nullptr : reinterpret_cast<btCollisionObject**>(&m_softBodies[0]), m_softBodies.size());
}
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

#ifndef BT_DEFORMABLE_MULTIBODY_DYNAMICS_WORLD_H
#define BT_DEFORMABLE_MULTIBODY_DYNAMICS_WORLD_H

#include "btSoftMultiBodyDynamicsWorld.h"
#include "btDeformableLagrangianForce.h"
#include "btDeformableMassSpringForce.h"
// #include "btDeformableBodySolver.h"
#include "btDeformableMultiBodyConstraintSolver.h"
#include "btSoftBodyHelpers.h"
#include "BulletCollision/CollisionDispatch/btSimulationIslandManager.h"
#include <functional>
typedef btAlignedObjectArray<btSoftBody*> btSoftBodyArray;

class btDeformableBodySolver;
class btDeformableLagrangianForce;
struct MultiBodyInplaceSolverIslandCallback;
struct DeformableBodyInplaceSolverIslandCallback;
class btDeformableMultiBodyConstraintSolver;

typedef btAlignedObjectArray<btSoftBody*> btSoftBodyArray;

class btDeformableMultiBodyDynamicsWorld : public btMultiBodyDynamicsWorld
{
	typedef btAlignedObjectArray<btVector3> TVStack;
	///Solver classes that encapsulate multiple deformable bodies for solving
	btDeformableBodySolver* m_deformableBodySolver;
	btSoftBodyArray m_softBodies;
	int m_drawFlags;
	bool m_drawNodeTree;
	bool m_drawFaceTree;
	bool m_drawClusterTree;
	btSoftBodyWorldInfo m_sbi;
	btScalar m_internalTime;
	int m_ccdIterations;
	bool m_implicit;
	bool m_lineSearch;
	bool m_useProjection;
	DeformableBodyInplaceSolverIslandCallback* m_solverDeformableBodyIslandCallback;

	typedef void (*btSolverCallback)(btScalar time, btDeformableMultiBodyDynamicsWorld* world);
	btSolverCallback m_solverCallback;

protected:
	virtual void internalSingleStepSimulation(btScalar timeStep);

	virtual void integrateTransforms(btScalar timeStep);

	void positionCorrection(btScalar timeStep);

	void solveConstraints(btScalar timeStep);

	void updateActivationState(btScalar timeStep);

	void clearGravity();

	void updateLastSafeTransforms();

	void addSoftsWithSelfCollisionCheckToOverlappingPairs();

public:
	btDeformableMultiBodyDynamicsWorld(btDispatcher* dispatcher, btBroadphaseInterface* pairCache, btDeformableMultiBodyConstraintSolver* constraintSolver, btCollisionConfiguration* collisionConfiguration, btDeformableBodySolver* deformableBodySolver = 0);

	virtual int stepSimulation(btScalar timeStep, int maxSubSteps = 1, btScalar fixedTimeStep = btScalar(1.) / btScalar(60.));

	virtual void debugDrawWorld();

	void setSolverCallback(btSolverCallback cb)
	{
		m_solverCallback = cb;
	}

	virtual ~btDeformableMultiBodyDynamicsWorld();

	virtual btMultiBodyDynamicsWorld* getMultiBodyDynamicsWorld()
	{
		return (btMultiBodyDynamicsWorld*)(this);
	}

	virtual const btMultiBodyDynamicsWorld* getMultiBodyDynamicsWorld() const
	{
		return (const btMultiBodyDynamicsWorld*)(this);
	}

	virtual btDynamicsWorldType getWorldType() const
	{
		return BT_DEFORMABLE_MULTIBODY_DYNAMICS_WORLD;
	}

	virtual void predictUnconstraintMotion(btScalar timeStep);

	virtual void addSoftBody(btSoftBody* body, int collisionFilterGroup = btBroadphaseProxy::DefaultFilter, int collisionFilterMask = btBroadphaseProxy::AllFilter);

	btSoftBodyArray& getSoftBodyArray()
	{
		return m_softBodies;
	}

	const btSoftBodyArray& getSoftBodyArray() const
	{
		return m_softBodies;
	}

	btSoftBodyWorldInfo& getWorldInfo()
	{
		return m_sbi;
	}

	const btSoftBodyWorldInfo& getWorldInfo() const
	{
		return m_sbi;
	}

	virtual void setGravity(const btVector3& gravity);

	void performDiscreteCollisionDetection() override;

	void reinitialize(btScalar timeStep);

	void applyRigidBodyGravity(btScalar timeStep);

	void beforeSolverCallbacks(btScalar timeStep);

	void afterSolverCallbacks(btScalar timeStep);

	void addForce(btSoftBody* psb, btDeformableLagrangianForce* force);

	void addForce(btDeformableLagrangianForce* force);

	void removeForce(btSoftBody* psb, btDeformableLagrangianForce* force);

	void removeForce(btDeformableLagrangianForce* force);

	void removeSoftBodyForce(btSoftBody* psb);

	void removeSoftBody(btSoftBody* body);

	void removeCollisionObject(btCollisionObject* collisionObject);

	int getDrawFlags() const { return (m_drawFlags); }
	void setDrawFlags(int f) { m_drawFlags = f; }

	void setupConstraints();

	void performDeformableCollisionDetection();

	void solveMultiBodyConstraints();

	void solveContactConstraints();

	void sortConstraints();

	void softBodySelfCollision();

	void setImplicit(bool implicit)
	{
		m_implicit = implicit;
	}

	void setLineSearch(bool lineSearch)
	{
		m_lineSearch = lineSearch;
	}

	void setUseProjection(bool useProjection)
	{
		m_useProjection = useProjection;
	}

	void applyRepulsionForce(btScalar timeStep);

	void performGeometricCollisions(btScalar timeStep);

	struct btDeformableSingleRayCallback : public btBroadphaseRayCallback
	{
		btVector3 m_rayFromWorld;
		btVector3 m_rayToWorld;
		btTransform m_rayFromTrans;
		btTransform m_rayToTrans;
		btVector3 m_hitNormal;

		const btDeformableMultiBodyDynamicsWorld* m_world;
		btCollisionWorld::RayResultCallback& m_resultCallback;

		btDeformableSingleRayCallback(const btVector3& rayFromWorld, const btVector3& rayToWorld, const btDeformableMultiBodyDynamicsWorld* world, btCollisionWorld::RayResultCallback& resultCallback)
			: m_rayFromWorld(rayFromWorld),
			  m_rayToWorld(rayToWorld),
			  m_world(world),
			  m_resultCallback(resultCallback)
		{
			m_rayFromTrans.setIdentity();
			m_rayFromTrans.setOrigin(m_rayFromWorld);
			m_rayToTrans.setIdentity();
			m_rayToTrans.setOrigin(m_rayToWorld);

			btVector3 rayDir = (rayToWorld - rayFromWorld);

			rayDir.normalize();
			///what about division by zero? --> just set rayDirection[i] to INF/1e30
			m_rayDirectionInverse[0] = rayDir[0] == btScalar(0.0) ? btScalar(1e30) : btScalar(1.0) / rayDir[0];
			m_rayDirectionInverse[1] = rayDir[1] == btScalar(0.0) ? btScalar(1e30) : btScalar(1.0) / rayDir[1];
			m_rayDirectionInverse[2] = rayDir[2] == btScalar(0.0) ? btScalar(1e30) : btScalar(1.0) / rayDir[2];
			m_signs[0] = m_rayDirectionInverse[0] < 0.0;
			m_signs[1] = m_rayDirectionInverse[1] < 0.0;
			m_signs[2] = m_rayDirectionInverse[2] < 0.0;

			m_lambda_max = rayDir.dot(m_rayToWorld - m_rayFromWorld);
		}

		virtual bool process(const btBroadphaseProxy* proxy)
		{
			///terminate further ray tests, once the closestHitFraction reached zero
			if (m_resultCallback.m_closestHitFraction == btScalar(0.f))
				return false;

			btCollisionObject* collisionObject = (btCollisionObject*)proxy->m_clientObject;

			//only perform raycast if filterMask matches
			if (m_resultCallback.needsCollision(collisionObject->getBroadphaseHandle()))
			{
				//RigidcollisionObject* collisionObject = ctrl->GetRigidcollisionObject();
				//btVector3 collisionObjectAabbMin,collisionObjectAabbMax;
#if 0
#ifdef RECALCULATE_AABB
                btVector3 collisionObjectAabbMin,collisionObjectAabbMax;
                collisionObject->getCollisionShape()->getAabb(collisionObject->getWorldTransform(),collisionObjectAabbMin,collisionObjectAabbMax);
#else
                //getBroadphase()->getAabb(collisionObject->getBroadphaseHandle(),collisionObjectAabbMin,collisionObjectAabbMax);
                const btVector3& collisionObjectAabbMin = collisionObject->getBroadphaseHandle()->m_aabbMin;
                const btVector3& collisionObjectAabbMax = collisionObject->getBroadphaseHandle()->m_aabbMax;
#endif
#endif
				//btScalar hitLambda = m_resultCallback.m_closestHitFraction;
				//culling already done by broadphase
				//if (btRayAabb(m_rayFromWorld,m_rayToWorld,collisionObjectAabbMin,collisionObjectAabbMax,hitLambda,m_hitNormal))
				{
					m_world->rayTestSingle(m_rayFromTrans, m_rayToTrans,
										   collisionObject,
										   collisionObject->getCollisionShape(),
										   collisionObject->getWorldTransform(),
										   m_resultCallback);
				}
			}
			return true;
		}
	};

	void rayTest(const btVector3& rayFromWorld, const btVector3& rayToWorld, RayResultCallback& resultCallback) const
	{
		BT_PROFILE("rayTest");
		/// use the broadphase to accelerate the search for objects, based on their aabb
		/// and for each object with ray-aabb overlap, perform an exact ray test
		btDeformableSingleRayCallback rayCB(rayFromWorld, rayToWorld, this, resultCallback);

#ifndef USE_BRUTEFORCE_RAYBROADPHASE
		m_broadphasePairCache->rayTest(rayFromWorld, rayToWorld, rayCB);
#else
		for (int i = 0; i < this->getNumCollisionObjects(); i++)
		{
			rayCB.process(m_collisionObjects[i]->getBroadphaseHandle());
		}
#endif  //USE_BRUTEFORCE_RAYBROADPHASE
	}

	void rayTestSingle(const btTransform& rayFromTrans, const btTransform& rayToTrans,
					   btCollisionObject* collisionObject,
					   const btCollisionShape* collisionShape,
					   const btTransform& colObjWorldTransform,
					   RayResultCallback& resultCallback) const
	{
		if (collisionShape->isSoftBody())
		{
			btSoftBody* softBody = btSoftBody::upcast(collisionObject);
			if (softBody)
			{
				btSoftBody::sRayCast softResult;
				if (softBody->rayFaceTest(rayFromTrans.getOrigin(), rayToTrans.getOrigin(), softResult))
				{
					if (softResult.fraction <= resultCallback.m_closestHitFraction)
					{
						btCollisionWorld::LocalShapeInfo shapeInfo;
						shapeInfo.m_shapePart = 0;
						shapeInfo.m_triangleIndex = softResult.index;
						// get the normal
						btVector3 rayDir = rayToTrans.getOrigin() - rayFromTrans.getOrigin();
						btVector3 normal = -rayDir;
						normal.normalize();
						{
							normal = softBody->m_faces[softResult.index].m_normal;
							if (normal.dot(rayDir) > 0)
							{
								// normal always point toward origin of the ray
								normal = -normal;
							}
						}

						btCollisionWorld::LocalRayResult rayResult(collisionObject,
																   &shapeInfo,
																   normal,
																   softResult.fraction);
						bool normalInWorldSpace = true;
						resultCallback.addSingleResult(rayResult, normalInWorldSpace);
					}
				}
			}
		}
		else
		{
			btCollisionWorld::rayTestSingle(rayFromTrans, rayToTrans, collisionObject, collisionShape, colObjWorldTransform, resultCallback);
		}
	}
};

#endif  //BT_DEFORMABLE_MULTIBODY_DYNAMICS_WORLD_H
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

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/

#ifndef BT_NEOHOOKEAN_H
#define BT_NEOHOOKEAN_H

#include "btDeformableLagrangianForce.h"
#include "LinearMath/btQuickprof.h"
#include "LinearMath/btImplicitQRSVD.h"
// This energy is as described in https://graphics.pixar.com/library/StableElasticity/paper.pdf
class btDeformableNeoHookeanForce : public btDeformableLagrangianForce
{
public:
	typedef btAlignedObjectArray<btVector3> TVStack;
	btScalar m_mu, m_lambda;  // Lame Parameters
	btScalar m_E, m_nu;       // Young's modulus and Poisson ratio
	btScalar m_mu_damp, m_lambda_damp;
	btDeformableNeoHookeanForce() : m_mu(1), m_lambda(1)
	{
		btScalar damping = 0.05;
		m_mu_damp = damping * m_mu;
		m_lambda_damp = damping * m_lambda;
		updateYoungsModulusAndPoissonRatio();
	}

	btDeformableNeoHookeanForce(btScalar mu, btScalar lambda, btScalar damping = 0.05) : m_mu(mu), m_lambda(lambda)
	{
		m_mu_damp = damping * m_mu;
		m_lambda_damp = damping * m_lambda;
		updateYoungsModulusAndPoissonRatio();
	}

	void updateYoungsModulusAndPoissonRatio()
	{
		// conversion from Lame Parameters to Young's modulus and Poisson ratio
		// https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
		m_E = m_mu * (3 * m_lambda + 2 * m_mu) / (m_lambda + m_mu);
		m_nu = m_lambda * 0.5 / (m_mu + m_lambda);
	}

	void updateLameParameters()
	{
		// conversion from Young's modulus and Poisson ratio to Lame Parameters
		// https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
		m_mu = m_E * 0.5 / (1 + m_nu);
		m_lambda = m_E * m_nu / ((1 + m_nu) * (1 - 2 * m_nu));
	}

	void setYoungsModulus(btScalar E)
	{
		m_E = E;
		updateLameParameters();
	}

	void setPoissonRatio(btScalar nu)
	{
		m_nu = nu;
		updateLameParameters();
	}

	void setDamping(btScalar damping)
	{
		m_mu_damp = damping * m_mu;
		m_lambda_damp = damping * m_lambda;
	}

	void setLameParameters(btScalar mu, btScalar lambda)
	{
		m_mu = mu;
		m_lambda = lambda;
		updateYoungsModulusAndPoissonRatio();
	}

	virtual void addScaledForces(btScalar scale, TVStack& force)
	{
		addScaledDampingForce(scale, force);
		addScaledElasticForce(scale, force);
	}

	virtual void addScaledExplicitForce(btScalar scale, TVStack& force)
	{
		addScaledElasticForce(scale, force);
	}

	// The damping matrix is calculated using the time n state as described in https://www.math.ucla.edu/~jteran/papers/GSSJT15.pdf to allow line search
	virtual void addScaledDampingForce(btScalar scale, TVStack& force)
	{
		if (m_mu_damp == 0 && m_lambda_damp == 0)
			return;
		int numNodes = getNumNodes();
		btAssert(numNodes <= force.size());
		btVector3 grad_N_hat_1st_col = btVector3(-1, -1, -1);
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_tetras.size(); ++j)
			{
				btSoftBody::Tetra& tetra = psb->m_tetras[j];
				btSoftBody::Node* node0 = tetra.m_n[0];
				btSoftBody::Node* node1 = tetra.m_n[1];
				btSoftBody::Node* node2 = tetra.m_n[2];
				btSoftBody::Node* node3 = tetra.m_n[3];
				size_t id0 = node0->index;
				size_t id1 = node1->index;
				size_t id2 = node2->index;
				size_t id3 = node3->index;
				btMatrix3x3 dF = DsFromVelocity(node0, node1, node2, node3) * tetra.m_Dm_inverse;
				btMatrix3x3 I;
				I.setIdentity();
				btMatrix3x3 dP = (dF + dF.transpose()) * m_mu_damp + I * (dF[0][0] + dF[1][1] + dF[2][2]) * m_lambda_damp;
				//                firstPiolaDampingDifferential(psb->m_tetraScratchesTn[j], dF, dP);
				btVector3 df_on_node0 = dP * (tetra.m_Dm_inverse.transpose() * grad_N_hat_1st_col);
				btMatrix3x3 df_on_node123 = dP * tetra.m_Dm_inverse.transpose();

				// damping force differential
				btScalar scale1 = scale * tetra.m_element_measure;
				force[id0] -= scale1 * df_on_node0;
				force[id1] -= scale1 * df_on_node123.getColumn(0);
				force[id2] -= scale1 * df_on_node123.getColumn(1);
				force[id3] -= scale1 * df_on_node123.getColumn(2);
			}
		}
	}

	virtual double totalElasticEnergy(btScalar dt)
	{
		double energy = 0;
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_tetraScratches.size(); ++j)
			{
				btSoftBody::Tetra& tetra = psb->m_tetras[j];
				btSoftBody::TetraScratch& s = psb->m_tetraScratches[j];
				energy += tetra.m_element_measure * elasticEnergyDensity(s);
			}
		}
		return energy;
	}

	// The damping energy is formulated as in https://www.math.ucla.edu/~jteran/papers/GSSJT15.pdf to allow line search
	virtual double totalDampingEnergy(btScalar dt)
	{
		double energy = 0;
		int sz = 0;
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				sz = btMax(sz, psb->m_nodes[j].index);
			}
		}
		TVStack dampingForce;
		dampingForce.resize(sz + 1);
		for (int i = 0; i < dampingForce.size(); ++i)
			dampingForce[i].setZero();
		addScaledDampingForce(0.5, dampingForce);
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				const btSoftBody::Node& node = psb->m_nodes[j];
				energy -= dampingForce[node.index].dot(node.m_v) / dt;
			}
		}
		return energy;
	}

	double elasticEnergyDensity(const btSoftBody::TetraScratch& s)
	{
		double density = 0;
		density += m_mu * 0.5 * (s.m_trace - 3.);
		density += m_lambda * 0.5 * (s.m_J - 1. - 0.75 * m_mu / m_lambda) * (s.m_J - 1. - 0.75 * m_mu / m_lambda);
		density -= m_mu * 0.5 * log(s.m_trace + 1);
		return density;
	}

	virtual void addScaledElasticForce(btScalar scale, TVStack& force)
	{
		int numNodes = getNumNodes();
		btAssert(numNodes <= force.size());
		btVector3 grad_N_hat_1st_col = btVector3(-1, -1, -1);
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			btScalar max_p = psb->m_cfg.m_maxStress;
			btScalar averagePrincipalStress = 0.0;
			for (int j = 0; j < psb->m_tetras.size(); ++j)
			{
				btSoftBody::Tetra& tetra = psb->m_tetras[j];
				btMatrix3x3 P;
				firstPiola(psb->m_tetraScratches[j], P);
				btScalar trPTP = (P[0].length() + P[1].length() + P[2].length());
				averagePrincipalStress += trPTP;
#define USE_SVD 1
#ifdef USE_SVD
				if (max_p > 0)
				{
					// since we want to clamp the principal stress to max_p, we only need to
					// calculate SVD when sigma_0^2 + sigma_1^2 + sigma_2^2 > max_p * max_p
					btScalar trPTP = (P[0].length2() + P[1].length2() + P[2].length2());
					if (trPTP > max_p * max_p)
					{
						btMatrix3x3 U, V;
						btVector3 sigma;
						singularValueDecomposition(P, U, sigma, V);
						sigma[0] = btMin(sigma[0], max_p);
						sigma[1] = btMin(sigma[1], max_p);
						sigma[2] = btMin(sigma[2], max_p);
						sigma[0] = btMax(sigma[0], -max_p);
						sigma[1] = btMax(sigma[1], -max_p);
						sigma[2] = btMax(sigma[2], -max_p);
						btMatrix3x3 Sigma;
						Sigma.setIdentity();
						Sigma[0][0] = sigma[0];
						Sigma[1][1] = sigma[1];
						Sigma[2][2] = sigma[2];
						P = U * Sigma * V.transpose();
					}
				}
#endif
				//                btVector3 force_on_node0 = P * (tetra.m_Dm_inverse.transpose()*grad_N_hat_1st_col);
				btMatrix3x3 force_on_node123 = P * tetra.m_Dm_inverse.transpose();
				btVector3 force_on_node0 = force_on_node123 * grad_N_hat_1st_col;

				btSoftBody::Node* node0 = tetra.m_n[0];
				btSoftBody::Node* node1 = tetra.m_n[1];
				btSoftBody::Node* node2 = tetra.m_n[2];
				btSoftBody::Node* node3 = tetra.m_n[3];
				size_t id0 = node0->index;
				size_t id1 = node1->index;
				size_t id2 = node2->index;
				size_t id3 = node3->index;

				// elastic force
				btScalar scale1 = scale * tetra.m_element_measure;
				force[id0] -= scale1 * force_on_node0;
				force[id1] -= scale1 * force_on_node123.getColumn(0);
				force[id2] -= scale1 * force_on_node123.getColumn(1);
				force[id3] -= scale1 * force_on_node123.getColumn(2);
			}
			averagePrincipalStress /= psb->m_tetras.size();
			psb->m_averagePrincipalStress = averagePrincipalStress;
			//fprintf(stderr, "m_averagePrincipalStress %f\n", psb->m_averagePrincipalStress);
		}
	}

	// The damping matrix is calculated using the time n state as described in https://www.math.ucla.edu/~jteran/papers/GSSJT15.pdf to allow line search
	virtual void addScaledDampingForceDifferential(btScalar scale, const TVStack& dv, TVStack& df)
	{
		if (m_mu_damp == 0 && m_lambda_damp == 0)
			return;
		int numNodes = getNumNodes();
		btAssert(numNodes <= df.size());
		btVector3 grad_N_hat_1st_col = btVector3(-1, -1, -1);
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_tetras.size(); ++j)
			{
				btSoftBody::Tetra& tetra = psb->m_tetras[j];
				btSoftBody::Node* node0 = tetra.m_n[0];
				btSoftBody::Node* node1 = tetra.m_n[1];
				btSoftBody::Node* node2 = tetra.m_n[2];
				btSoftBody::Node* node3 = tetra.m_n[3];
				size_t id0 = node0->index;
				size_t id1 = node1->index;
				size_t id2 = node2->index;
				size_t id3 = node3->index;
				btMatrix3x3 dF = Ds(id0, id1, id2, id3, dv) * tetra.m_Dm_inverse;
				btMatrix3x3 I;
				I.setIdentity();
				btMatrix3x3 dP = (dF + dF.transpose()) * m_mu_damp + I * (dF[0][0] + dF[1][1] + dF[2][2]) * m_lambda_damp;
				//                firstPiolaDampingDifferential(psb->m_tetraScratchesTn[j], dF, dP);
				//                btVector3 df_on_node0 = dP * (tetra.m_Dm_inverse.transpose()*grad_N_hat_1st_col);
				btMatrix3x3 df_on_node123 = dP * tetra.m_Dm_inverse.transpose();
				btVector3 df_on_node0 = df_on_node123 * grad_N_hat_1st_col;

				// damping force differential
				btScalar scale1 = scale * tetra.m_element_measure;
				df[id0] -= scale1 * df_on_node0;
				df[id1] -= scale1 * df_on_node123.getColumn(0);
				df[id2] -= scale1 * df_on_node123.getColumn(1);
				df[id3] -= scale1 * df_on_node123.getColumn(2);
			}
		}
	}

	virtual void buildDampingForceDifferentialDiagonal(btScalar scale, TVStack& diagA) {}

	virtual void addScaledElasticForceDifferential(btScalar scale, const TVStack& dx, TVStack& df)
	{
		int numNodes = getNumNodes();
		btAssert(numNodes <= df.size());
		btVector3 grad_N_hat_1st_col = btVector3(-1, -1, -1);
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			if (!psb->isActive() || psb->isStaticObject())
			{
				continue;
			}
			for (int j = 0; j < psb->m_tetras.size(); ++j)
			{
				btSoftBody::Tetra& tetra = psb->m_tetras[j];
				btSoftBody::Node* node0 = tetra.m_n[0];
				btSoftBody::Node* node1 = tetra.m_n[1];
				btSoftBody::Node* node2 = tetra.m_n[2];
				btSoftBody::Node* node3 = tetra.m_n[3];
				size_t id0 = node0->index;
				size_t id1 = node1->index;
				size_t id2 = node2->index;
				size_t id3 = node3->index;
				btMatrix3x3 dF = Ds(id0, id1, id2, id3, dx) * tetra.m_Dm_inverse;
				btMatrix3x3 dP;
				firstPiolaDifferential(psb->m_tetraScratches[j], dF, dP);
				//                btVector3 df_on_node0 = dP * (tetra.m_Dm_inverse.transpose()*grad_N_hat_1st_col);
				btMatrix3x3 df_on_node123 = dP * tetra.m_Dm_inverse.transpose();
				btVector3 df_on_node0 = df_on_node123 * grad_N_hat_1st_col;

				// elastic force differential
				btScalar scale1 = scale * tetra.m_element_measure;
				df[id0] -= scale1 * df_on_node0;
				df[id1] -= scale1 * df_on_node123.getColumn(0);
				df[id2] -= scale1 * df_on_node123.getColumn(1);
				df[id3] -= scale1 * df_on_node123.getColumn(2);
			}
		}
	}

	void firstPiola(const btSoftBody::TetraScratch& s, btMatrix3x3& P)
	{
		btScalar c1 = (m_mu * (1. - 1. / (s.m_trace + 1.)));
		btScalar c2 = (m_lambda * (s.m_J - 1.) - 0.75 * m_mu);
		P = s.m_F * c1 + s.m_cofF * c2;
	}

	// Let P be the first piola stress.
	// This function calculates the dP = dP/dF * dF
	void firstPiolaDifferential(const btSoftBody::TetraScratch& s, const btMatrix3x3& dF, btMatrix3x3& dP)
	{
		btScalar c1 = m_mu * (1. - 1. / (s.m_trace + 1.));
		btScalar c2 = (2. * m_mu) * DotProduct(s.m_F, dF) * (1. / ((1. + s.m_trace) * (1. + s.m_trace)));
		btScalar c3 = (m_lambda * DotProduct(s.m_cofF, dF));
		dP = dF * c1 + s.m_F * c2;
		addScaledCofactorMatrixDifferential(s.m_F, dF, m_lambda * (s.m_J - 1.) - 0.75 * m_mu, dP);
		dP += s.m_cofF * c3;
	}

	// Let Q be the damping stress.
	// This function calculates the dP = dQ/dF * dF
	void firstPiolaDampingDifferential(const btSoftBody::TetraScratch& s, const btMatrix3x3& dF, btMatrix3x3& dP)
	{
		btScalar c1 = (m_mu_damp * (1. - 1. / (s.m_trace + 1.)));
		btScalar c2 = ((2. * m_mu_damp) * DotProduct(s.m_F, dF) * (1. / ((1. + s.m_trace) * (1. + s.m_trace))));
		btScalar c3 = (m_lambda_damp * DotProduct(s.m_cofF, dF));
		dP = dF * c1 + s.m_F * c2;
		addScaledCofactorMatrixDifferential(s.m_F, dF, m_lambda_damp * (s.m_J - 1.) - 0.75 * m_mu_damp, dP);
		dP += s.m_cofF * c3;
	}

	btScalar DotProduct(const btMatrix3x3& A, const btMatrix3x3& B)
	{
		btScalar ans = 0;
		for (int i = 0; i < 3; ++i)
		{
			ans += A[i].dot(B[i]);
		}
		return ans;
	}

	// Let C(A) be the cofactor of the matrix A
	// Let H = the derivative of C(A) with respect to A evaluated at F = A
	// This function calculates H*dF
	void addScaledCofactorMatrixDifferential(const btMatrix3x3& F, const btMatrix3x3& dF, btScalar scale, btMatrix3x3& M)
	{
		M[0][0] += scale * (dF[1][1] * F[2][2] + F[1][1] * dF[2][2] - dF[2][1] * F[1][2] - F[2][1] * dF[1][2]);
		M[1][0] += scale * (dF[2][1] * F[0][2] + F[2][1] * dF[0][2] - dF[0][1] * F[2][2] - F[0][1] * dF[2][2]);
		M[2][0] += scale * (dF[0][1] * F[1][2] + F[0][1] * dF[1][2] - dF[1][1] * F[0][2] - F[1][1] * dF[0][2]);
		M[0][1] += scale * (dF[2][0] * F[1][2] + F[2][0] * dF[1][2] - dF[1][0] * F[2][2] - F[1][0] * dF[2][2]);
		M[1][1] += scale * (dF[0][0] * F[2][2] + F[0][0] * dF[2][2] - dF[2][0] * F[0][2] - F[2][0] * dF[0][2]);
		M[2][1] += scale * (dF[1][0] * F[0][2] + F[1][0] * dF[0][2] - dF[0][0] * F[1][2] - F[0][0] * dF[1][2]);
		M[0][2] += scale * (dF[1][0] * F[2][1] + F[1][0] * dF[2][1] - dF[2][0] * F[1][1] - F[2][0] * dF[1][1]);
		M[1][2] += scale * (dF[2][0] * F[0][1] + F[2][0] * dF[0][1] - dF[0][0] * F[2][1] - F[0][0] * dF[2][1]);
		M[2][2] += scale * (dF[0][0] * F[1][1] + F[0][0] * dF[1][1] - dF[1][0] * F[0][1] - F[1][0] * dF[0][1]);
	}

	virtual btScalar getYoungsModulus() const
	{
		return m_E;
	}

	virtual btScalar getPoissonRatio() const
	{
		return m_nu;
	}

	virtual btDeformableLagrangianForceType getForceType()
	{
		return BT_NEOHOOKEAN_FORCE;
	}
};
#endif /* BT_NEOHOOKEAN_H */
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

#ifndef BT_KRYLOV_SOLVER_H
#define BT_KRYLOV_SOLVER_H
#include <iostream>
#include <cmath>
#include <limits>
#include <LinearMath/btAlignedObjectArray.h>
#include <LinearMath/btVector3.h>
#include <LinearMath/btScalar.h>
#include "LinearMath/btQuickprof.h"

template <class MatrixX>
class btKrylovSolver
{
	typedef btAlignedObjectArray<btVector3> TVStack;

public:
	int m_maxIterations;
	btScalar m_tolerance;
	btKrylovSolver(int maxIterations, btScalar tolerance)
		: m_maxIterations(maxIterations), m_tolerance(tolerance)
	{
	}

	virtual ~btKrylovSolver() {}

	virtual int solve(MatrixX& A, TVStack& x, const TVStack& b, bool verbose = false) = 0;

	virtual void reinitialize(const TVStack& b) = 0;

	virtual SIMD_FORCE_INLINE TVStack sub(const TVStack& a, const TVStack& b)
	{
		// c = a-b
		btAssert(a.size() == b.size());
		TVStack c;
		c.resize(a.size());
		for (int i = 0; i < a.size(); ++i)
		{
			c[i] = a[i] - b[i];
		}
		return c;
	}

	virtual SIMD_FORCE_INLINE btScalar squaredNorm(const TVStack& a)
	{
		return dot(a, a);
	}

	virtual SIMD_FORCE_INLINE btScalar norm(const TVStack& a)
	{
		btScalar ret = 0;
		for (int i = 0; i < a.size(); ++i)
		{
			for (int d = 0; d < 3; ++d)
			{
				ret = btMax(ret, btFabs(a[i][d]));
			}
		}
		return ret;
	}

	virtual SIMD_FORCE_INLINE btScalar dot(const TVStack& a, const TVStack& b)
	{
		btScalar ans(0);
		for (int i = 0; i < a.size(); ++i)
			ans += a[i].dot(b[i]);
		return ans;
	}

	virtual SIMD_FORCE_INLINE void multAndAddTo(btScalar s, const TVStack& a, TVStack& result)
	{
		//        result += s*a
		btAssert(a.size() == result.size());
		for (int i = 0; i < a.size(); ++i)
			result[i] += s * a[i];
	}

	virtual SIMD_FORCE_INLINE TVStack multAndAdd(btScalar s, const TVStack& a, const TVStack& b)
	{
		// result = a*s + b
		TVStack result;
		result.resize(a.size());
		for (int i = 0; i < a.size(); ++i)
			result[i] = s * a[i] + b[i];
		return result;
	}

	virtual SIMD_FORCE_INLINE void setTolerance(btScalar tolerance)
	{
		m_tolerance = tolerance;
	}
};
#endif /* BT_KRYLOV_SOLVER_H */
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

#ifndef BT_PRECONDITIONER_H
#define BT_PRECONDITIONER_H

class Preconditioner
{
public:
	typedef btAlignedObjectArray<btVector3> TVStack;
	virtual void operator()(const TVStack& x, TVStack& b) = 0;
	virtual void reinitialize(bool nodeUpdated) = 0;
	virtual ~Preconditioner() {}
};

class DefaultPreconditioner : public Preconditioner
{
public:
	virtual void operator()(const TVStack& x, TVStack& b)
	{
		btAssert(b.size() == x.size());
		for (int i = 0; i < b.size(); ++i)
			b[i] = x[i];
	}
	virtual void reinitialize(bool nodeUpdated)
	{
	}

	virtual ~DefaultPreconditioner() {}
};

class MassPreconditioner : public Preconditioner
{
	btAlignedObjectArray<btScalar> m_inv_mass;
	const btAlignedObjectArray<btSoftBody*>& m_softBodies;

public:
	MassPreconditioner(const btAlignedObjectArray<btSoftBody*>& softBodies)
		: m_softBodies(softBodies)
	{
	}

	virtual void reinitialize(bool nodeUpdated)
	{
		if (nodeUpdated)
		{
			m_inv_mass.clear();
			for (int i = 0; i < m_softBodies.size(); ++i)
			{
				btSoftBody* psb = m_softBodies[i];
				for (int j = 0; j < psb->m_nodes.size(); ++j)
					m_inv_mass.push_back(psb->m_nodes[j].m_im);
			}
		}
	}

	virtual void operator()(const TVStack& x, TVStack& b)
	{
		btAssert(b.size() == x.size());
		btAssert(m_inv_mass.size() <= x.size());
		for (int i = 0; i < m_inv_mass.size(); ++i)
		{
			b[i] = x[i] * m_inv_mass[i];
		}
		for (int i = m_inv_mass.size(); i < b.size(); ++i)
		{
			b[i] = x[i];
		}
	}
};

class KKTPreconditioner : public Preconditioner
{
	const btAlignedObjectArray<btSoftBody*>& m_softBodies;
	const btDeformableContactProjection& m_projections;
	const btAlignedObjectArray<btDeformableLagrangianForce*>& m_lf;
	TVStack m_inv_A, m_inv_S;
	const btScalar& m_dt;
	const bool& m_implicit;

public:
	KKTPreconditioner(const btAlignedObjectArray<btSoftBody*>& softBodies, const btDeformableContactProjection& projections, const btAlignedObjectArray<btDeformableLagrangianForce*>& lf, const btScalar& dt, const bool& implicit)
		: m_softBodies(softBodies), m_projections(projections), m_lf(lf), m_dt(dt), m_implicit(implicit)
	{
	}

	virtual void reinitialize(bool nodeUpdated)
	{
		if (nodeUpdated)
		{
			int num_nodes = 0;
			for (int i = 0; i < m_softBodies.size(); ++i)
			{
				btSoftBody* psb = m_softBodies[i];
				num_nodes += psb->m_nodes.size();
			}
			m_inv_A.resize(num_nodes);
		}
		buildDiagonalA(m_inv_A);
		for (int i = 0; i < m_inv_A.size(); ++i)
		{
			//            printf("A[%d] = %f, %f, %f \n", i, m_inv_A[i][0], m_inv_A[i][1], m_inv_A[i][2]);
			for (int d = 0; d < 3; ++d)
			{
				m_inv_A[i][d] = (m_inv_A[i][d] == 0) ? 0.0 : 1.0 / m_inv_A[i][d];
			}
		}
		m_inv_S.resize(m_projections.m_lagrangeMultipliers.size());
		//        printf("S.size() = %d \n", m_inv_S.size());
		buildDiagonalS(m_inv_A, m_inv_S);
		for (int i = 0; i < m_inv_S.size(); ++i)
		{
			//            printf("S[%d] = %f, %f, %f \n", i, m_inv_S[i][0], m_inv_S[i][1], m_inv_S[i][2]);
			for (int d = 0; d < 3; ++d)
			{
				m_inv_S[i][d] = (m_inv_S[i][d] == 0) ? 0.0 : 1.0 / m_inv_S[i][d];
			}
		}
	}

	void buildDiagonalA(TVStack& diagA) const
	{
		size_t counter = 0;
		for (int i = 0; i < m_softBodies.size(); ++i)
		{
			btSoftBody* psb = m_softBodies[i];
			for (int j = 0; j < psb->m_nodes.size(); ++j)
			{
				const btSoftBody::Node& node = psb->m_nodes[j];
				diagA[counter] = (node.m_im == 0) ? btVector3(0, 0, 0) : btVector3(1.0 / node.m_im, 1.0 / node.m_im, 1.0 / node.m_im);
				++counter;
			}
		}
		if (m_implicit)
		{
			printf("implicit not implemented\n");
			btAssert(false);
		}
		for (int i = 0; i < m_lf.size(); ++i)
		{
			// add damping matrix
			m_lf[i]->buildDampingForceDifferentialDiagonal(-m_dt, diagA);
		}
	}

	void buildDiagonalS(const TVStack& inv_A, TVStack& diagS)
	{
		for (int c = 0; c < m_projections.m_lagrangeMultipliers.size(); ++c)
		{
			// S[k,k] = e_k^T * C A_d^-1 C^T * e_k
			const LagrangeMultiplier& lm = m_projections.m_lagrangeMultipliers[c];
			btVector3& t = diagS[c];
			t.setZero();
			for (int j = 0; j < lm.m_num_constraints; ++j)
			{
				for (int i = 0; i < lm.m_num_nodes; ++i)
				{
					for (int d = 0; d < 3; ++d)
					{
						t[j] += inv_A[lm.m_indices[i]][d] * lm.m_dirs[j][d] * lm.m_dirs[j][d] * lm.m_weights[i] * lm.m_weights[i];
					}
				}
			}
		}
	}
//#define USE_FULL_PRECONDITIONER
#ifndef USE_FULL_PRECONDITIONER
	virtual void operator()(const TVStack& x, TVStack& b)
	{
		btAssert(b.size() == x.size());
		for (int i = 0; i < m_inv_A.size(); ++i)
		{
			b[i] = x[i] * m_inv_A[i];
		}
		int offset = m_inv_A.size();
		for (int i = 0; i < m_inv_S.size(); ++i)
		{
			b[i + offset] = x[i + offset] * m_inv_S[i];
		}
	}
#else
	virtual void operator()(const TVStack& x, TVStack& b)
	{
		btAssert(b.size() == x.size());
		int offset = m_inv_A.size();

		for (int i = 0; i < m_inv_A.size(); ++i)
		{
			b[i] = x[i] * m_inv_A[i];
		}

		for (int i = 0; i < m_inv_S.size(); ++i)
		{
			b[i + offset].setZero();
		}

		for (int c = 0; c < m_projections.m_lagrangeMultipliers.size(); ++c)
		{
			const LagrangeMultiplier& lm = m_projections.m_lagrangeMultipliers[c];
			// C * x
			for (int d = 0; d < lm.m_num_constraints; ++d)
			{
				for (int i = 0; i < lm.m_num_nodes; ++i)
				{
					b[offset + c][d] += lm.m_weights[i] * b[lm.m_indices[i]].dot(lm.m_dirs[d]);
				}
			}
		}

		for (int i = 0; i < m_inv_S.size(); ++i)
		{
			b[i + offset] = b[i + offset] * m_inv_S[i];
		}

		for (int i = 0; i < m_inv_A.size(); ++i)
		{
			b[i].setZero();
		}

		for (int c = 0; c < m_projections.m_lagrangeMultipliers.size(); ++c)
		{
			// C^T * lambda
			const LagrangeMultiplier& lm = m_projections.m_lagrangeMultipliers[c];
			for (int i = 0; i < lm.m_num_nodes; ++i)
			{
				for (int j = 0; j < lm.m_num_constraints; ++j)
				{
					b[lm.m_indices[i]] += b[offset + c][j] * lm.m_weights[i] * lm.m_dirs[j];
				}
			}
		}

		for (int i = 0; i < m_inv_A.size(); ++i)
		{
			b[i] = (x[i] - b[i]) * m_inv_A[i];
		}

		TVStack t;
		t.resize(b.size());
		for (int i = 0; i < m_inv_S.size(); ++i)
		{
			t[i + offset] = x[i + offset] * m_inv_S[i];
		}
		for (int i = 0; i < m_inv_A.size(); ++i)
		{
			t[i].setZero();
		}
		for (int c = 0; c < m_projections.m_lagrangeMultipliers.size(); ++c)
		{
			// C^T * lambda
			const LagrangeMultiplier& lm = m_projections.m_lagrangeMultipliers[c];
			for (int i = 0; i < lm.m_num_nodes; ++i)
			{
				for (int j = 0; j < lm.m_num_constraints; ++j)
				{
					t[lm.m_indices[i]] += t[offset + c][j] * lm.m_weights[i] * lm.m_dirs[j];
				}
			}
		}
		for (int i = 0; i < m_inv_A.size(); ++i)
		{
			b[i] += t[i] * m_inv_A[i];
		}

		for (int i = 0; i < m_inv_S.size(); ++i)
		{
			b[i + offset] -= x[i + offset] * m_inv_S[i];
		}
	}
#endif
};

#endif /* BT_PRECONDITIONER_H */
/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/

///btSoftBody implementation by Nathanael Presson

#include "btSoftBodyInternals.h"
#include "BulletSoftBody/btSoftBodySolvers.h"
#include "btSoftBodyData.h"
#include "LinearMath/btSerializer.h"
#include "LinearMath/btImplicitQRSVD.h"
#include "LinearMath/btAlignedAllocator.h"
#include "BulletDynamics/Featherstone/btMultiBodyLinkCollider.h"
#include "BulletDynamics/Featherstone/btMultiBodyConstraint.h"
#include "BulletCollision/NarrowPhaseCollision/btGjkEpa2.h"
#include "BulletCollision/CollisionShapes/btTriangleShape.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <list>
//
static inline btDbvtNode* buildTreeBottomUp(btAlignedObjectArray<btDbvtNode*>& leafNodes, btAlignedObjectArray<btAlignedObjectArray<int>>& adj)
{
	int N = leafNodes.size();
	if (N == 0)
	{
		return NULL;
	}
	while (N > 1)
	{
		btAlignedObjectArray<bool> marked;
		btAlignedObjectArray<btDbvtNode*> newLeafNodes;
		btAlignedObjectArray<std::pair<int, int>> childIds;
		btAlignedObjectArray<btAlignedObjectArray<int>> newAdj;
		marked.resize(N);
		for (int i = 0; i < N; ++i)
			marked[i] = false;

		// pair adjacent nodes into new(parent) node
		for (int i = 0; i < N; ++i)
		{
			if (marked[i])
				continue;
			bool merged = false;
			for (int j = 0; j < adj[i].size(); ++j)
			{
				int n = adj[i][j];
				if (!marked[adj[i][j]])
				{
					btDbvtNode* node = new (btAlignedAlloc(sizeof(btDbvtNode), 16)) btDbvtNode();
					node->parent = NULL;
					node->childs[0] = leafNodes[i];
					node->childs[1] = leafNodes[n];
					leafNodes[i]->parent = node;
					leafNodes[n]->parent = node;
					newLeafNodes.push_back(node);
					childIds.push_back(std::make_pair(i, n));
					merged = true;
					marked[n] = true;
					break;
				}
			}
			if (!merged)
			{
				newLeafNodes.push_back(leafNodes[i]);
				childIds.push_back(std::make_pair(i, -1));
			}
			marked[i] = true;
		}
		// update adjacency matrix
		newAdj.resize(newLeafNodes.size());
		for (int i = 0; i < newLeafNodes.size(); ++i)
		{
			for (int j = i + 1; j < newLeafNodes.size(); ++j)
			{
				bool neighbor = false;
				const btAlignedObjectArray<int>& leftChildNeighbors = adj[childIds[i].first];
				for (int k = 0; k < leftChildNeighbors.size(); ++k)
				{
					if (leftChildNeighbors[k] == childIds[j].first || leftChildNeighbors[k] == childIds[j].second)
					{
						neighbor = true;
						break;
					}
				}
				if (!neighbor && childIds[i].second != -1)
				{
					const btAlignedObjectArray<int>& rightChildNeighbors = adj[childIds[i].second];
					for (int k = 0; k < rightChildNeighbors.size(); ++k)
					{
						if (rightChildNeighbors[k] == childIds[j].first || rightChildNeighbors[k] == childIds[j].second)
						{
							neighbor = true;
							break;
						}
					}
				}
				if (neighbor)
				{
					newAdj[i].push_back(j);
					newAdj[j].push_back(i);
				}
			}
		}
		leafNodes = newLeafNodes;
		//this assignment leaks memory, the assignment doesn't do a deep copy, for now a manual copy
		//adj = newAdj;
		adj.clear();
		adj.resize(newAdj.size());
		for (int i = 0; i < newAdj.size(); i++)
		{
			for (int j = 0; j < newAdj[i].size(); j++)
			{
				adj[i].push_back(newAdj[i][j]);
			}
		}
		N = leafNodes.size();
	}
	return leafNodes[0];
}

//
btSoftBody::btSoftBody(btSoftBodyWorldInfo* worldInfo, int node_count, const btVector3* x, const btScalar* m, btCollisionShape* collisionShape)
	: m_softBodySolver(0), m_worldInfo(worldInfo)
{
	/* Init		*/
	initDefaults(collisionShape);

	/* Default material	*/
	Material* pm = appendMaterial();
	pm->m_kLST = 1;
	pm->m_kAST = 1;
	pm->m_kVST = 1;
	pm->m_flags = fMaterial::Default;

	/* Nodes			*/
	const btScalar margin = getCollisionShape()->getMargin();
	m_nodes.resize(node_count);
	m_X.resize(node_count);
	for (int i = 0, ni = node_count; i < ni; ++i)
	{
		Node& n = m_nodes[i];
		ZeroInitialize(n);
		n.m_x = x ? *x++ : btVector3(0, 0, 0);
		n.m_q = n.m_x;
		n.m_im = m ? *m++ : 1;
		n.m_im = n.m_im > 0 ? 1 / n.m_im : 0;
		n.m_leaf = m_ndbvt.insert(btDbvtVolume::FromCR(n.m_x, margin), &n);
		n.m_material = pm;
		n.m_safe.copy_from_node(n);
		m_X[i] = n.m_x;
	}
	updateBounds();
	setCollisionQuadrature(3);
	m_fdbvnt = 0;
}

btSoftBody::btSoftBody(btSoftBodyWorldInfo* worldInfo)
	: m_worldInfo(worldInfo)
{
	initDefaults(nullptr);
}

void btSoftBody::initDefaults(btCollisionShape* collisionShape)
{
	m_internalType = CO_SOFT_BODY;
	m_cfg.aeromodel = eAeroModel::V_Point;
	m_cfg.kVCF = 1;
	m_cfg.kDG = 0;
	m_cfg.kLF = 0;
	m_cfg.kDP = 0;
	m_cfg.kPR = 0;
	m_cfg.kVC = 0;
	m_cfg.kDF = (btScalar)0.2;
	m_cfg.kMT = 0;
	m_cfg.kCHR = (btScalar)1.0;
	m_cfg.kKHR = (btScalar)0.1;
	m_cfg.kSHR = (btScalar)1.0;
	m_cfg.kAHR = (btScalar)0.7;
	m_cfg.kSRHR_CL = (btScalar)0.1;
	m_cfg.kSKHR_CL = (btScalar)1;
	m_cfg.kSSHR_CL = (btScalar)0.5;
	m_cfg.kSR_SPLT_CL = (btScalar)0.5;
	m_cfg.kSK_SPLT_CL = (btScalar)0.5;
	m_cfg.kSS_SPLT_CL = (btScalar)0.5;
	m_cfg.maxvolume = (btScalar)1;
	m_cfg.timescale = 1;
	m_cfg.viterations = 0;
	m_cfg.piterations = 1;
	m_cfg.diterations = 0;
	m_cfg.citerations = 4;
	m_cfg.drag = 0;
	m_cfg.m_maxStress = 0;
	m_cfg.collisions = fCollision::Default;
	m_pose.m_bvolume = false;
	m_pose.m_bframe = false;
	m_pose.m_volume = 0;
	m_pose.m_com = btVector3(0, 0, 0);
	m_pose.m_rot.setIdentity();
	m_pose.m_scl.setIdentity();
	m_tag = 0;
	m_timeacc = 0;
	m_bUpdateRtCst = true;
	m_bounds[0] = btVector3(0, 0, 0);
	m_bounds[1] = btVector3(0, 0, 0);
	m_worldTransform.setIdentity();
	setSolver(eSolverPresets::Positions);

	/* Collision shape	*/
	///for now, create a collision shape internally
	if (collisionShape)
		m_collisionShape = collisionShape;
	else
	{
		m_collisionShape = new btSoftBodyCollisionShape(this);
		m_collisionShape->setMargin(0.25f);
	}

	m_worldTransform.setIdentity();

	m_windVelocity = btVector3(0, 0, 0);
	m_restLengthScale = btScalar(1.0);
	m_dampingCoefficient = 1.0;
	m_sleepingThreshold = .04;
	m_useSelfCollision = false;
	m_collisionFlags = 0;
	m_softSoftCollision = false;
	m_maxSpeedSquared = 0;
	m_repulsionStiffness = 0.5;
	m_gravityFactor = 1;
	m_cacheBarycenter = false;
	m_fdbvnt = 0;

	// reduced flag
	m_reducedModel = false;

	m_averagePrincipalStress = 0.0;
	m_lastSafeApplyDepthThreshold = 5.0;
	m_lastSafeApplyVelocityDamping = 0.5;
	m_softVsSoftContactStiffness = 100.0;
}

//
btSoftBody::~btSoftBody()
{
	//for now, delete the internal shape
	delete m_collisionShape;
	int i;

	releaseClusters();
	for (i = 0; i < m_materials.size(); ++i)
		btAlignedFree(m_materials[i]);
	for (i = 0; i < m_joints.size(); ++i)
		btAlignedFree(m_joints[i]);
	if (m_fdbvnt)
		delete m_fdbvnt;

	for (int i = 0; i < m_anchors.size(); ++i)
	{
		const Anchor& c = m_anchors[i];
		if (c.m_body)
			c.m_body->removeAnchorRef(this);
	}

	for (int i = 0; i < m_deformableAnchors.size(); ++i)
	{
		const DeformableNodeRigidAnchor& c = m_deformableAnchors[i];
		if (c.m_body)
			c.m_body->removeAnchorRef(this);
	}
}

//
bool btSoftBody::checkLink(int node0, int node1) const
{
	return (checkLink(&m_nodes[node0], &m_nodes[node1]));
}

//
bool btSoftBody::checkLink(const Node* node0, const Node* node1) const
{
	const Node* n[] = {node0, node1};
	for (int i = 0, ni = m_links.size(); i < ni; ++i)
	{
		const Link& l = m_links[i];
		if ((l.m_n[0] == n[0] && l.m_n[1] == n[1]) ||
			(l.m_n[0] == n[1] && l.m_n[1] == n[0]))
		{
			return (true);
		}
	}
	return (false);
}

//
bool btSoftBody::checkFace(int node0, int node1, int node2) const
{
	const Node* n[] = {&m_nodes[node0],
					   &m_nodes[node1],
					   &m_nodes[node2]};
	for (int i = 0, ni = m_faces.size(); i < ni; ++i)
	{
		const Face& f = m_faces[i];
		int c = 0;
		for (int j = 0; j < 3; ++j)
		{
			if ((f.m_n[j] == n[0]) ||
				(f.m_n[j] == n[1]) ||
				(f.m_n[j] == n[2]))
				c |= 1 << j;
			else
				break;
		}
		if (c == 7) return (true);
	}
	return (false);
}

//
btSoftBody::Material* btSoftBody::appendMaterial()
{
	Material* pm = new (btAlignedAlloc(sizeof(Material), 16)) Material();
	if (m_materials.size() > 0)
		*pm = *m_materials[0];
	else
		ZeroInitialize(*pm);
	m_materials.push_back(pm);
	return (pm);
}

//
void btSoftBody::appendNote(const char* text,
							const btVector3& o,
							const btVector4& c,
							Node* n0,
							Node* n1,
							Node* n2,
							Node* n3)
{
	Note n;
	ZeroInitialize(n);
	n.m_rank = 0;
	n.m_text = text;
	n.m_offset = o;
	n.m_coords[0] = c.x();
	n.m_coords[1] = c.y();
	n.m_coords[2] = c.z();
	n.m_coords[3] = c.w();
	n.m_nodes[0] = n0;
	n.m_rank += n0 ? 1 : 0;
	n.m_nodes[1] = n1;
	n.m_rank += n1 ? 1 : 0;
	n.m_nodes[2] = n2;
	n.m_rank += n2 ? 1 : 0;
	n.m_nodes[3] = n3;
	n.m_rank += n3 ? 1 : 0;
	m_notes.push_back(n);
}

//
void btSoftBody::appendNote(const char* text,
							const btVector3& o,
							Node* feature)
{
	appendNote(text, o, btVector4(1, 0, 0, 0), feature);
}

//
void btSoftBody::appendNote(const char* text,
							const btVector3& o,
							Link* feature)
{
	static const btScalar w = 1 / (btScalar)2;
	appendNote(text, o, btVector4(w, w, 0, 0), feature->m_n[0],
			   feature->m_n[1]);
}

//
void btSoftBody::appendNote(const char* text,
							const btVector3& o,
							Face* feature)
{
	static const btScalar w = 1 / (btScalar)3;
	appendNote(text, o, btVector4(w, w, w, 0), feature->m_n[0],
			   feature->m_n[1],
			   feature->m_n[2]);
}

//
void btSoftBody::appendNode(const btVector3& x, btScalar m)
{
	if (m_nodes.capacity() == m_nodes.size())
	{
		pointersToIndices();
		m_nodes.reserve(m_nodes.size() * 2 + 1);
		indicesToPointers();
	}
	const btScalar margin = getCollisionShape()->getMargin();
	m_nodes.push_back(Node());
	Node& n = m_nodes[m_nodes.size() - 1];
	ZeroInitialize(n);
	n.m_x = x;
	n.m_q = n.m_x;
	n.m_im = m > 0 ? 1 / m : 0;
	n.m_material = m_materials[0];
	n.m_safe.copy_from_node(n);
	n.m_leaf = m_ndbvt.insert(btDbvtVolume::FromCR(n.m_x, margin), &n);
}

//
void btSoftBody::appendLink(int model, Material* mat)
{
	Link l;
	if (model >= 0)
		l = m_links[model];
	else
	{
		ZeroInitialize(l);
		l.m_material = mat ? mat : m_materials[0];
	}
	m_links.push_back(l);
}

//
void btSoftBody::appendLink(int node0,
							int node1,
							Material* mat,
							bool bcheckexist)
{
	appendLink(&m_nodes[node0], &m_nodes[node1], mat, bcheckexist);
}

//
void btSoftBody::appendLink(Node* node0,
							Node* node1,
							Material* mat,
							bool bcheckexist)
{
	if ((!bcheckexist) || (!checkLink(node0, node1)))
	{
		appendLink(-1, mat);
		Link& l = m_links[m_links.size() - 1];
		l.m_n[0] = node0;
		l.m_n[1] = node1;
		l.m_rl = (l.m_n[0]->m_x - l.m_n[1]->m_x).length();
		m_bUpdateRtCst = true;
	}
}

//
void btSoftBody::appendFace(int model, Material* mat)
{
	Face f;
	if (model >= 0)
	{
		f = m_faces[model];
	}
	else
	{
		ZeroInitialize(f);
		f.m_material = mat ? mat : m_materials[0];
	}
	m_faces.push_back(f);
}

//
void btSoftBody::appendFace(int node0, int node1, int node2, Material* mat)
{
	if (node0 == node1)
		return;
	if (node1 == node2)
		return;
	if (node2 == node0)
		return;

	appendFace(-1, mat);
	Face& f = m_faces[m_faces.size() - 1];
	btAssert(node0 != node1);
	btAssert(node1 != node2);
	btAssert(node2 != node0);
	f.m_n[0] = &m_nodes[node0];
	f.m_n[1] = &m_nodes[node1];
	f.m_n[2] = &m_nodes[node2];
	f.m_ra = AreaOf(f.m_n[0]->m_x,
					f.m_n[1]->m_x,
					f.m_n[2]->m_x);
	m_bUpdateRtCst = true;
}

//
void btSoftBody::appendTetra(int model, Material* mat)
{
	Tetra t;
	if (model >= 0)
		t = m_tetras[model];
	else
	{
		ZeroInitialize(t);
		t.m_material = mat ? mat : m_materials[0];
	}
	m_tetras.push_back(t);
}

//
void btSoftBody::appendTetra(int node0,
							 int node1,
							 int node2,
							 int node3,
							 Material* mat)
{
	appendTetra(-1, mat);
	auto tetraIndex = m_tetras.size() - 1;
	Tetra& t = m_tetras[tetraIndex];
	t.m_n[0] = &m_nodes[node0];
	t.m_n[1] = &m_nodes[node1];
	t.m_n[2] = &m_nodes[node2];
	t.m_n[3] = &m_nodes[node3];
	btAssert(m_nodes[node0].m_tetraMembershipCount < Node::kMaxTetraMembershipPerNode);
	btAssert(m_nodes[node1].m_tetraMembershipCount < Node::kMaxTetraMembershipPerNode);
	btAssert(m_nodes[node2].m_tetraMembershipCount < Node::kMaxTetraMembershipPerNode);
	btAssert(m_nodes[node3].m_tetraMembershipCount < Node::kMaxTetraMembershipPerNode);
	m_nodes[node0].m_tetraMembership[m_nodes[node0].m_tetraMembershipCount++] = tetraIndex;
	m_nodes[node1].m_tetraMembership[m_nodes[node1].m_tetraMembershipCount++] = tetraIndex;
	m_nodes[node2].m_tetraMembership[m_nodes[node2].m_tetraMembershipCount++] = tetraIndex;
	m_nodes[node3].m_tetraMembership[m_nodes[node3].m_tetraMembershipCount++] = tetraIndex;
	t.m_rv = VolumeOf(t.m_n[0]->m_x, t.m_n[1]->m_x, t.m_n[2]->m_x, t.m_n[3]->m_x);
	m_bUpdateRtCst = true;
}

//

void btSoftBody::appendAnchor(int node, btRigidBody* body, bool disableCollisionBetweenLinkedBodies, btScalar influence)
{
	btVector3 local = body->getWorldTransform().inverse() * m_nodes[node].m_x;
	appendAnchor(node, body, local, disableCollisionBetweenLinkedBodies, influence);
}

//
void btSoftBody::appendAnchor(int node, btRigidBody* body, const btVector3& localPivot, bool disableCollisionBetweenLinkedBodies, btScalar influence)
{
	if (disableCollisionBetweenLinkedBodies)
	{
		if (m_collisionDisabledObjects.findLinearSearch(body) == m_collisionDisabledObjects.size())
		{
			m_collisionDisabledObjects.push_back(body);
		}
	}

	Anchor a;
	a.m_node = &m_nodes[node];
	a.m_body = body;
	a.m_local = localPivot;
	a.m_node->m_battach = 1;
	a.m_influence = influence;
	m_anchors.push_back(a);
	body->addAnchorRef(this);
}

//
void btSoftBody::appendDeformableAnchor(int node, btRigidBody* body, uint32_t userIndex)
{
	DeformableNodeRigidAnchor c;
	btSoftBody::Node& n = m_nodes[node];
	const btScalar ima = n.m_im;
	const btScalar imb = body->getInvMass();
	btVector3 nrm;
	const btCollisionShape* shp = body->getCollisionShape();
	const btTransform& wtr = body->getWorldTransform();
	btScalar dst =
		m_worldInfo->m_sparsesdf.Evaluate(
			wtr.invXform(m_nodes[node].m_x),
			shp,
			nrm,
			0);

	c.m_cti.m_colObj = body;
	c.m_body = body;
	c.m_cti.m_normal = wtr.getBasis() * nrm;
	c.m_cti.m_offset = dst;
	c.m_node = &m_nodes[node];
	const btScalar fc = m_cfg.kDF * body->getFriction();
	c.m_c2 = ima;
	c.m_c3 = fc;
	c.m_c4 = body->isStaticOrKinematicObject() ? m_cfg.kKHR : m_cfg.kCHR;
	static const btMatrix3x3 iwiStatic(0, 0, 0, 0, 0, 0, 0, 0, 0);
	const btMatrix3x3& iwi = body->getInvInertiaTensorWorld();
	const btVector3 ra = n.m_x - wtr.getOrigin();

	c.m_c0 = ImpulseMatrix(1, ima, imb, iwi, ra);
	c.m_c1 = ra;
	c.m_local = body->getWorldTransform().inverse() * m_nodes[node].m_x;
	c.m_node->m_battach = 1;
	c.m_userIndex = userIndex;
	m_deformableAnchors.push_back(c);
	body->addAnchorRef(this);
}

void btSoftBody::removeAnchor(int node)
{
	const btSoftBody::Node& n = m_nodes[node];
	for (int i = 0; i < m_anchors.size();)
	{
		const Anchor& c = m_anchors[i];
		if (c.m_node == &n)
		{
			if (c.m_body)
				c.m_body->removeAnchorRef(this);
			m_anchors.removeAtIndex(i);
		}
		else
		{
			i++;
		}
	}
}

void btSoftBody::removeDeformableAnchor(int node)
{
	const btSoftBody::Node& n = m_nodes[node];
	for (int i = 0; i < m_deformableAnchors.size();)
	{
		const DeformableNodeRigidAnchor& c = m_deformableAnchors[i];
		if (c.m_node == &n)
		{
			if (c.m_body)
				c.m_body->removeAnchorRef(this);
			m_deformableAnchors.removeAtIndex(i);
		}
		else
		{
			i++;
		}
	}
}

int btSoftBody::removeDeformableAnchorByUserIndex(int userIndex)
{
	int removedCount = 0;
	for (int i = 0; i < m_deformableAnchors.size();)
	{
		const DeformableNodeRigidAnchor& c = m_deformableAnchors[i];
		if (c.m_userIndex == userIndex)
		{
			if (c.m_body)
				c.m_body->removeAnchorRef(this);
			m_deformableAnchors.removeAtIndex(i);
			++removedCount;
		}
		else
		{
			i++;
		}
	}
	return removedCount;
}

std::vector<int> btSoftBody::getDeformableAnchorByUserIndex(int userIndex) const
{
	std::vector<int> result;
	for (int i = 0; i < m_deformableAnchors.size(); ++i)
	{
		const DeformableNodeRigidAnchor& c = m_deformableAnchors[i];
		if (c.m_userIndex == userIndex)
			result.emplace_back(i);
	}
	return result;
}

//
void btSoftBody::appendDeformableAnchor(int node, btMultiBodyLinkCollider* link, uint32_t userIndex)
{
	DeformableNodeRigidAnchor c;
	btSoftBody::Node& n = m_nodes[node];
	const btScalar ima = n.m_im;
	btVector3 nrm;
	const btCollisionShape* shp = link->getCollisionShape();
	const btTransform& wtr = link->getWorldTransform();
	btScalar dst =
		m_worldInfo->m_sparsesdf.Evaluate(
			wtr.invXform(m_nodes[node].m_x),
			shp,
			nrm,
			0);
	c.m_cti.m_colObj = link;
	c.m_cti.m_normal = wtr.getBasis() * nrm;
	c.m_cti.m_offset = dst;
	c.m_node = &m_nodes[node];
	const btScalar fc = m_cfg.kDF * link->getFriction();
	c.m_c2 = ima;
	c.m_c3 = fc;
	c.m_c4 = link->isStaticOrKinematicObject() ? m_cfg.kKHR : m_cfg.kCHR;
	btVector3 normal = c.m_cti.m_normal;
	btVector3 t1 = generateUnitOrthogonalVector(normal);
	btVector3 t2 = btCross(normal, t1);
	btMultiBodyJacobianData jacobianData_normal, jacobianData_t1, jacobianData_t2;
	findJacobian(link, jacobianData_normal, c.m_node->m_x, normal);
	findJacobian(link, jacobianData_t1, c.m_node->m_x, t1);
	findJacobian(link, jacobianData_t2, c.m_node->m_x, t2);

	btScalar* J_n = &jacobianData_normal.m_jacobians[0];
	btScalar* J_t1 = &jacobianData_t1.m_jacobians[0];
	btScalar* J_t2 = &jacobianData_t2.m_jacobians[0];

	btScalar* u_n = &jacobianData_normal.m_deltaVelocitiesUnitImpulse[0];
	btScalar* u_t1 = &jacobianData_t1.m_deltaVelocitiesUnitImpulse[0];
	btScalar* u_t2 = &jacobianData_t2.m_deltaVelocitiesUnitImpulse[0];

	btMatrix3x3 rot(normal.getX(), normal.getY(), normal.getZ(),
					t1.getX(), t1.getY(), t1.getZ(),
					t2.getX(), t2.getY(), t2.getZ());  // world frame to local frame
	const int ndof = link->m_multiBody->getNumDofs() + 6;
	btMatrix3x3 local_impulse_matrix = (Diagonal(n.m_im) + OuterProduct(J_n, J_t1, J_t2, u_n, u_t1, u_t2, ndof)).inverse();
	c.m_c0 = rot.transpose() * local_impulse_matrix * rot;
	c.jacobianData_normal = jacobianData_normal;
	c.jacobianData_t1 = jacobianData_t1;
	c.jacobianData_t2 = jacobianData_t2;
	c.t1 = t1;
	c.t2 = t2;
	const btVector3 ra = n.m_x - wtr.getOrigin();
	c.m_c1 = ra;
	c.m_local = link->getWorldTransform().inverse() * m_nodes[node].m_x;
	c.m_node->m_battach = 1;
	c.m_userIndex = userIndex;
	m_deformableAnchors.push_back(c);
}
//
void btSoftBody::appendLinearJoint(const LJoint::Specs& specs, Cluster* body0, Body body1)
{
	LJoint* pj = new (btAlignedAlloc(sizeof(LJoint), 16)) LJoint();
	pj->m_bodies[0] = body0;
	pj->m_bodies[1] = body1;
	pj->m_refs[0] = pj->m_bodies[0].xform().inverse() * specs.position;
	pj->m_refs[1] = pj->m_bodies[1].xform().inverse() * specs.position;
	pj->m_cfm = specs.cfm;
	pj->m_erp = specs.erp;
	pj->m_split = specs.split;
	m_joints.push_back(pj);
}

//
void btSoftBody::appendLinearJoint(const LJoint::Specs& specs, Body body)
{
	appendLinearJoint(specs, m_clusters[0], body);
}

//
void btSoftBody::appendLinearJoint(const LJoint::Specs& specs, btSoftBody* body)
{
	appendLinearJoint(specs, m_clusters[0], body->m_clusters[0]);
}

//
void btSoftBody::appendAngularJoint(const AJoint::Specs& specs, Cluster* body0, Body body1)
{
	AJoint* pj = new (btAlignedAlloc(sizeof(AJoint), 16)) AJoint();
	pj->m_bodies[0] = body0;
	pj->m_bodies[1] = body1;
	pj->m_refs[0] = pj->m_bodies[0].xform().inverse().getBasis() * specs.axis;
	pj->m_refs[1] = pj->m_bodies[1].xform().inverse().getBasis() * specs.axis;
	pj->m_cfm = specs.cfm;
	pj->m_erp = specs.erp;
	pj->m_split = specs.split;
	pj->m_icontrol = specs.icontrol;
	m_joints.push_back(pj);
}

//
void btSoftBody::appendAngularJoint(const AJoint::Specs& specs, Body body)
{
	appendAngularJoint(specs, m_clusters[0], body);
}

//
void btSoftBody::appendAngularJoint(const AJoint::Specs& specs, btSoftBody* body)
{
	appendAngularJoint(specs, m_clusters[0], body->m_clusters[0]);
}

//
void btSoftBody::addForce(const btVector3& force)
{
	for (int i = 0, ni = m_nodes.size(); i < ni; ++i) addForce(force, i);
}

//
void btSoftBody::addForce(const btVector3& force, int node)
{
	Node& n = m_nodes[node];
	if (n.m_im > 0)
	{
		n.m_f += force;
	}
}

void btSoftBody::addAeroForceToNode(const btVector3& windVelocity, int nodeIndex)
{
	btAssert(nodeIndex >= 0 && nodeIndex < m_nodes.size());

	const btScalar dt = m_sst.sdt;
	const btScalar kLF = m_cfg.kLF;
	const btScalar kDG = m_cfg.kDG;
	//const btScalar kPR = m_cfg.kPR;
	//const btScalar kVC = m_cfg.kVC;
	const bool as_lift = kLF > 0;
	const bool as_drag = kDG > 0;
	const bool as_aero = as_lift || as_drag;
	const bool as_vaero = as_aero && (m_cfg.aeromodel < btSoftBody::eAeroModel::F_TwoSided);

	Node& n = m_nodes[nodeIndex];

	if (n.m_im > 0)
	{
		btSoftBody::sMedium medium;

		EvaluateMedium(m_worldInfo, n.m_x, medium);
		medium.m_velocity = windVelocity;
		medium.m_density = m_worldInfo->air_density;

		/* Aerodynamics			*/
		if (as_vaero)
		{
			const btVector3 rel_v = n.m_v - medium.m_velocity;
			const btScalar rel_v_len = rel_v.length();
			const btScalar rel_v2 = rel_v.length2();

			if (rel_v2 > SIMD_EPSILON)
			{
				const btVector3 rel_v_nrm = rel_v.normalized();
				btVector3 nrm = n.m_n;

				if (m_cfg.aeromodel == btSoftBody::eAeroModel::V_TwoSidedLiftDrag)
				{
					nrm *= (btScalar)((btDot(nrm, rel_v) < 0) ? -1 : +1);
					btVector3 fDrag(0, 0, 0);
					btVector3 fLift(0, 0, 0);

					btScalar n_dot_v = nrm.dot(rel_v_nrm);
					btScalar tri_area = 0.5f * n.m_area;

					fDrag = 0.5f * kDG * medium.m_density * rel_v2 * tri_area * n_dot_v * (-rel_v_nrm);

					// Check angle of attack
					// cos(10º) = 0.98480
					if (0 < n_dot_v && n_dot_v < 0.98480f)
						fLift = 0.5f * kLF * medium.m_density * rel_v_len * tri_area * btSqrt(1.0f - n_dot_v * n_dot_v) * (nrm.cross(rel_v_nrm).cross(rel_v_nrm));

					// Check if the velocity change resulted by aero drag force exceeds the current velocity of the node.
					btVector3 del_v_by_fDrag = fDrag * n.m_im * m_sst.sdt;
					btScalar del_v_by_fDrag_len2 = del_v_by_fDrag.length2();
					btScalar v_len2 = n.m_v.length2();

					if (del_v_by_fDrag_len2 >= v_len2 && del_v_by_fDrag_len2 > 0)
					{
						btScalar del_v_by_fDrag_len = del_v_by_fDrag.length();
						btScalar v_len = n.m_v.length();
						fDrag *= btScalar(0.8) * (v_len / del_v_by_fDrag_len);
					}

					n.m_f += fDrag;
					n.m_f += fLift;
				}
				else if (m_cfg.aeromodel == btSoftBody::eAeroModel::V_Point || m_cfg.aeromodel == btSoftBody::eAeroModel::V_OneSided || m_cfg.aeromodel == btSoftBody::eAeroModel::V_TwoSided)
				{
					if (m_cfg.aeromodel == btSoftBody::eAeroModel::V_TwoSided)
						nrm *= (btScalar)((btDot(nrm, rel_v) < 0) ? -1 : +1);

					const btScalar dvn = btDot(rel_v, nrm);
					/* Compute forces	*/
					if (dvn > 0)
					{
						btVector3 force(0, 0, 0);
						const btScalar c0 = n.m_area * dvn * rel_v2 / 2;
						const btScalar c1 = c0 * medium.m_density;
						force += nrm * (-c1 * kLF);
						force += rel_v.normalized() * (-c1 * kDG);
						ApplyClampedForce(n, force, dt);
					}
				}
			}
		}
	}
}

void btSoftBody::addAeroForceToFace(const btVector3& windVelocity, int faceIndex)
{
	const btScalar dt = m_sst.sdt;
	const btScalar kLF = m_cfg.kLF;
	const btScalar kDG = m_cfg.kDG;
	//	const btScalar kPR = m_cfg.kPR;
	//	const btScalar kVC = m_cfg.kVC;
	const bool as_lift = kLF > 0;
	const bool as_drag = kDG > 0;
	const bool as_aero = as_lift || as_drag;
	const bool as_faero = as_aero && (m_cfg.aeromodel >= btSoftBody::eAeroModel::F_TwoSided);

	if (as_faero)
	{
		btSoftBody::Face& f = m_faces[faceIndex];

		btSoftBody::sMedium medium;

		const btVector3 v = (f.m_n[0]->m_v + f.m_n[1]->m_v + f.m_n[2]->m_v) / 3;
		const btVector3 x = (f.m_n[0]->m_x + f.m_n[1]->m_x + f.m_n[2]->m_x) / 3;
		EvaluateMedium(m_worldInfo, x, medium);
		medium.m_velocity = windVelocity;
		medium.m_density = m_worldInfo->air_density;
		const btVector3 rel_v = v - medium.m_velocity;
		const btScalar rel_v_len = rel_v.length();
		const btScalar rel_v2 = rel_v.length2();

		if (rel_v2 > SIMD_EPSILON)
		{
			const btVector3 rel_v_nrm = rel_v.normalized();
			btVector3 nrm = f.m_normal;

			if (m_cfg.aeromodel == btSoftBody::eAeroModel::F_TwoSidedLiftDrag)
			{
				nrm *= (btScalar)((btDot(nrm, rel_v) < 0) ? -1 : +1);

				btVector3 fDrag(0, 0, 0);
				btVector3 fLift(0, 0, 0);

				btScalar n_dot_v = nrm.dot(rel_v_nrm);
				btScalar tri_area = 0.5f * f.m_ra;

				fDrag = 0.5f * kDG * medium.m_density * rel_v2 * tri_area * n_dot_v * (-rel_v_nrm);

				// Check angle of attack
				// cos(10º) = 0.98480
				if (0 < n_dot_v && n_dot_v < 0.98480f)
					fLift = 0.5f * kLF * medium.m_density * rel_v_len * tri_area * btSqrt(1.0f - n_dot_v * n_dot_v) * (nrm.cross(rel_v_nrm).cross(rel_v_nrm));

				fDrag /= 3;
				fLift /= 3;

				for (int j = 0; j < 3; ++j)
				{
					if (f.m_n[j]->m_im > 0)
					{
						// Check if the velocity change resulted by aero drag force exceeds the current velocity of the node.
						btVector3 del_v_by_fDrag = fDrag * f.m_n[j]->m_im * m_sst.sdt;
						btScalar del_v_by_fDrag_len2 = del_v_by_fDrag.length2();
						btScalar v_len2 = f.m_n[j]->m_v.length2();

						if (del_v_by_fDrag_len2 >= v_len2 && del_v_by_fDrag_len2 > 0)
						{
							btScalar del_v_by_fDrag_len = del_v_by_fDrag.length();
							btScalar v_len = f.m_n[j]->m_v.length();
							fDrag *= btScalar(0.8) * (v_len / del_v_by_fDrag_len);
						}

						f.m_n[j]->m_f += fDrag;
						f.m_n[j]->m_f += fLift;
					}
				}
			}
			else if (m_cfg.aeromodel == btSoftBody::eAeroModel::F_OneSided || m_cfg.aeromodel == btSoftBody::eAeroModel::F_TwoSided)
			{
				if (m_cfg.aeromodel == btSoftBody::eAeroModel::F_TwoSided)
					nrm *= (btScalar)((btDot(nrm, rel_v) < 0) ? -1 : +1);

				const btScalar dvn = btDot(rel_v, nrm);
				/* Compute forces	*/
				if (dvn > 0)
				{
					btVector3 force(0, 0, 0);
					const btScalar c0 = f.m_ra * dvn * rel_v2;
					const btScalar c1 = c0 * medium.m_density;
					force += nrm * (-c1 * kLF);
					force += rel_v.normalized() * (-c1 * kDG);
					force /= 3;
					for (int j = 0; j < 3; ++j) ApplyClampedForce(*f.m_n[j], force, dt);
				}
			}
		}
	}
}

//
void btSoftBody::addVelocity(const btVector3& velocity)
{
	for (int i = 0, ni = m_nodes.size(); i < ni; ++i) addVelocity(velocity, i);
}

/* Set velocity for the entire body										*/
void btSoftBody::setVelocity(const btVector3& velocity)
{
	for (int i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		Node& n = m_nodes[i];
		if (n.m_im > 0)
		{
			n.m_v = velocity;
			n.m_vn = velocity;
		}
	}
}

//
void btSoftBody::addVelocity(const btVector3& velocity, int node)
{
	Node& n = m_nodes[node];
	if (n.m_im > 0)
	{
		n.m_v += velocity;
	}
}

//
void btSoftBody::setMass(int node, btScalar mass)
{
	m_nodes[node].m_im = mass > 0 ? 1 / mass : 0;
	m_bUpdateRtCst = true;
}

//
btScalar btSoftBody::getMass(int node) const
{
	return (m_nodes[node].m_im > 0 ? 1 / m_nodes[node].m_im : 0);
}

//
btScalar btSoftBody::getTotalMass() const
{
	btScalar mass = 0;
	for (int i = 0; i < m_nodes.size(); ++i)
	{
		mass += getMass(i);
	}
	return (mass);
}

//
void btSoftBody::setTotalMass(btScalar mass, bool fromfaces)
{
	int i;
	std::list<int> anchors;

	if (fromfaces)
	{
		for (i = 0; i < m_nodes.size(); ++i)
		{
			if (m_nodes[i].m_im == 0)
				anchors.emplace_back(i);
			else
				m_nodes[i].m_im = 0;
		}
		for (i = 0; i < m_faces.size(); ++i)
		{
			const Face& f = m_faces[i];
			const btScalar twicearea = AreaOf(f.m_n[0]->m_x,
											  f.m_n[1]->m_x,
											  f.m_n[2]->m_x);
			for (int j = 0; j < 3; ++j)
			{
				f.m_n[j]->m_im += twicearea;
			}
		}
		for (i = 0; i < m_nodes.size(); ++i)
		{
			m_nodes[i].m_im = 1 / m_nodes[i].m_im;
		}
	}
	const btScalar tm = getTotalMass();
	const btScalar itm = 1 / tm;
	for (i = 0; i < m_nodes.size(); ++i)
	{
		m_nodes[i].m_im /= itm * mass;
	}
	for (auto anchor : anchors)
		m_nodes[anchor].m_im = 0;

	m_bUpdateRtCst = true;
}

//
void btSoftBody::setTotalDensity(btScalar density)
{
	setTotalMass(getVolume() * density, true);
}

//
void btSoftBody::setVolumeMass(btScalar mass)
{
	btAlignedObjectArray<btScalar> ranks;
	std::list<int> anchors;
	ranks.resize(m_nodes.size(), 0);
	int i;

	for (i = 0; i < m_nodes.size(); ++i)
	{
		if (m_nodes[i].m_im == 0)
			anchors.emplace_back(i);
		else
			m_nodes[i].m_im = 0;
	}
	for (i = 0; i < m_tetras.size(); ++i)
	{
		const Tetra& t = m_tetras[i];
		for (int j = 0; j < 4; ++j)
		{
			t.m_n[j]->m_im += btFabs(t.m_rv);
			ranks[int(t.m_n[j] - &m_nodes[0])] += 1;
		}
	}
	for (i = 0; i < m_nodes.size(); ++i)
	{
		if (m_nodes[i].m_im > 0)
		{
			m_nodes[i].m_im = ranks[i] / m_nodes[i].m_im;
		}
	}
	for (auto anchor : anchors)
		m_nodes[anchor].m_im = 0;

	setTotalMass(mass, false);
}

//
void btSoftBody::setVolumeDensity(btScalar density)
{
	btScalar volume = 0;
	for (int i = 0; i < m_tetras.size(); ++i)
	{
		const Tetra& t = m_tetras[i];
		volume += btFabs(t.m_rv);
	}
	setVolumeMass(volume * density / 6);
}

//
btVector3 btSoftBody::getLinearVelocity()
{
	btVector3 total_momentum = btVector3(0, 0, 0);
	for (int i = 0; i < m_nodes.size(); ++i)
	{
		btScalar mass = m_nodes[i].m_im == 0 ? 0 : 1.0 / m_nodes[i].m_im;
		total_momentum += mass * m_nodes[i].m_v;
	}
	btScalar total_mass = getTotalMass();
	return total_mass == 0 ? total_momentum : total_momentum / total_mass;
}

//
void btSoftBody::setLinearVelocity(const btVector3& linVel)
{
	btVector3 old_vel = getLinearVelocity();
	btVector3 diff = linVel - old_vel;
	for (int i = 0; i < m_nodes.size(); ++i)
		m_nodes[i].m_v += diff;
}

//
void btSoftBody::setAngularVelocity(const btVector3& angVel)
{
	btVector3 old_vel = getLinearVelocity();
	btVector3 com = getCenterOfMass();
	for (int i = 0; i < m_nodes.size(); ++i)
	{
		m_nodes[i].m_v = angVel.cross(m_nodes[i].m_x - com) + old_vel;
	}
}

//
btTransform btSoftBody::getRigidTransform()
{
	btVector3 t = getCenterOfMass();
	btMatrix3x3 S;
	S.setZero();
	// Get rotation that minimizes L2 difference: \sum_i || RX_i + t - x_i ||
	// It's important to make sure that S has the correct signs.
	// SVD is only unique up to the ordering of singular values.
	// SVD will manipulate U and V to ensure the ordering of singular values. If all three singular
	// vaues are negative, SVD will permute colums of U to make two of them positive.
	for (int i = 0; i < m_nodes.size(); ++i)
	{
		S -= OuterProduct(m_X[i], t - m_nodes[i].m_x);
	}
	btVector3 sigma;
	btMatrix3x3 U, V;
	singularValueDecomposition(S, U, sigma, V);
	btMatrix3x3 R = V * U.transpose();
	btTransform trs;
	trs.setIdentity();
	trs.setOrigin(t);
	trs.setBasis(R);
	return trs;
}

//
void btSoftBody::transformTo(const btTransform& trs)
{
	// get the current best rigid fit
	btTransform current_transform = getRigidTransform();
	// apply transform in material space
	btTransform new_transform = trs * current_transform.inverse();
	transform(new_transform);
}

//
void btSoftBody::transform(const btTransform& trs)
{
	const btScalar margin = getCollisionShape()->getMargin();
	ATTRIBUTE_ALIGNED16(btDbvtVolume)
	vol;

	for (int i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		Node& n = m_nodes[i];
		n.m_x = trs * n.m_x;
		n.m_q = trs * n.m_q;
		n.m_n = trs.getBasis() * n.m_n;
		vol = btDbvtVolume::FromCR(n.m_x, margin);

		m_ndbvt.update(n.m_leaf, vol);
	}
	updateNormals();
	updateBounds();
	updateConstants();
}

//
void btSoftBody::translate(const btVector3& trs)
{
	btTransform t;
	t.setIdentity();
	t.setOrigin(trs);
	transform(t);
}

//
void btSoftBody::rotate(const btQuaternion& rot)
{
	btTransform t;
	t.setIdentity();
	t.setRotation(rot);
	transform(t);
}

//
void btSoftBody::scale(const btVector3& scl)
{
	const btScalar margin = getCollisionShape()->getMargin();
	ATTRIBUTE_ALIGNED16(btDbvtVolume)
	vol;

	for (int i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		Node& n = m_nodes[i];
		n.m_x *= scl;
		n.m_q *= scl;
		vol = btDbvtVolume::FromCR(n.m_x, margin);
		m_ndbvt.update(n.m_leaf, vol);
	}
	updateNormals();
	updateBounds();
	updateConstants();
	initializeDmInverse();
}

//
btScalar btSoftBody::getRestLengthScale()
{
	return m_restLengthScale;
}

//
void btSoftBody::setRestLengthScale(btScalar restLengthScale)
{
	for (int i = 0, ni = m_links.size(); i < ni; ++i)
	{
		Link& l = m_links[i];
		l.m_rl = l.m_rl / m_restLengthScale * restLengthScale;
		l.m_c1 = l.m_rl * l.m_rl;
	}
	m_restLengthScale = restLengthScale;

	if (getActivationState() == ISLAND_SLEEPING)
		activate();
}

//
void btSoftBody::setPose(bool bvolume, bool bframe)
{
	m_pose.m_bvolume = bvolume;
	m_pose.m_bframe = bframe;
	int i, ni;

	/* Weights		*/
	const btScalar omass = getTotalMass();
	const btScalar kmass = omass * m_nodes.size() * 1000;
	btScalar tmass = omass;
	m_pose.m_wgh.resize(m_nodes.size());
	for (i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		if (m_nodes[i].m_im <= 0) tmass += kmass;
	}
	for (i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		Node& n = m_nodes[i];
		m_pose.m_wgh[i] = n.m_im > 0 ? 1 / (m_nodes[i].m_im * tmass) : kmass / tmass;
	}
	/* Pos		*/
	const btVector3 com = evaluateCom();
	m_pose.m_pos.resize(m_nodes.size());
	for (i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		m_pose.m_pos[i] = m_nodes[i].m_x - com;
	}
	m_pose.m_volume = bvolume ? getVolume() : 0;
	m_pose.m_com = com;
	m_pose.m_rot.setIdentity();
	m_pose.m_scl.setIdentity();
	/* Aqq		*/
	m_pose.m_aqq[0] =
		m_pose.m_aqq[1] =
			m_pose.m_aqq[2] = btVector3(0, 0, 0);
	for (i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		const btVector3& q = m_pose.m_pos[i];
		const btVector3 mq = m_pose.m_wgh[i] * q;
		m_pose.m_aqq[0] += mq.x() * q;
		m_pose.m_aqq[1] += mq.y() * q;
		m_pose.m_aqq[2] += mq.z() * q;
	}
	m_pose.m_aqq = m_pose.m_aqq.inverse();

	updateConstants();
}

void btSoftBody::resetLinkRestLengths()
{
	for (int i = 0, ni = m_links.size(); i < ni; ++i)
	{
		Link& l = m_links[i];
		l.m_rl = (l.m_n[0]->m_x - l.m_n[1]->m_x).length();
		l.m_c1 = l.m_rl * l.m_rl;
	}
}

//
btScalar btSoftBody::getVolume() const
{
	btScalar vol = 0;
	if (m_nodes.size() > 0)
	{
		int i, ni;

		const btVector3 org = m_nodes[0].m_x;
		for (i = 0, ni = m_faces.size(); i < ni; ++i)
		{
			const Face& f = m_faces[i];
			vol += btDot(f.m_n[0]->m_x - org, btCross(f.m_n[1]->m_x - org, f.m_n[2]->m_x - org));
		}
		vol /= (btScalar)6;
	}
	return (vol);
}

//
int btSoftBody::clusterCount() const
{
	return (m_clusters.size());
}

//
btVector3 btSoftBody::clusterCom(const Cluster* cluster)
{
	btVector3 com(0, 0, 0);
	for (int i = 0, ni = cluster->m_nodes.size(); i < ni; ++i)
	{
		com += cluster->m_nodes[i]->m_x * cluster->m_masses[i];
	}
	return (com * cluster->m_imass);
}

//
btVector3 btSoftBody::clusterCom(int cluster) const
{
	return (clusterCom(m_clusters[cluster]));
}

//
btVector3 btSoftBody::clusterVelocity(const Cluster* cluster, const btVector3& rpos)
{
	return (cluster->m_lv + btCross(cluster->m_av, rpos));
}

//
void btSoftBody::clusterVImpulse(Cluster* cluster, const btVector3& rpos, const btVector3& impulse)
{
	const btVector3 li = cluster->m_imass * impulse;
	const btVector3 ai = cluster->m_invwi * btCross(rpos, impulse);
	cluster->m_vimpulses[0] += li;
	cluster->m_lv += li;
	cluster->m_vimpulses[1] += ai;
	cluster->m_av += ai;
	cluster->m_nvimpulses++;
}

//
void btSoftBody::clusterDImpulse(Cluster* cluster, const btVector3& rpos, const btVector3& impulse)
{
	const btVector3 li = cluster->m_imass * impulse;
	const btVector3 ai = cluster->m_invwi * btCross(rpos, impulse);
	cluster->m_dimpulses[0] += li;
	cluster->m_dimpulses[1] += ai;
	cluster->m_ndimpulses++;
}

//
void btSoftBody::clusterImpulse(Cluster* cluster, const btVector3& rpos, const Impulse& impulse)
{
	if (impulse.m_asVelocity) clusterVImpulse(cluster, rpos, impulse.m_velocity);
	if (impulse.m_asDrift) clusterDImpulse(cluster, rpos, impulse.m_drift);
}

//
void btSoftBody::clusterVAImpulse(Cluster* cluster, const btVector3& impulse)
{
	const btVector3 ai = cluster->m_invwi * impulse;
	cluster->m_vimpulses[1] += ai;
	cluster->m_av += ai;
	cluster->m_nvimpulses++;
}

//
void btSoftBody::clusterDAImpulse(Cluster* cluster, const btVector3& impulse)
{
	const btVector3 ai = cluster->m_invwi * impulse;
	cluster->m_dimpulses[1] += ai;
	cluster->m_ndimpulses++;
}

//
void btSoftBody::clusterAImpulse(Cluster* cluster, const Impulse& impulse)
{
	if (impulse.m_asVelocity) clusterVAImpulse(cluster, impulse.m_velocity);
	if (impulse.m_asDrift) clusterDAImpulse(cluster, impulse.m_drift);
}

//
void btSoftBody::clusterDCImpulse(Cluster* cluster, const btVector3& impulse)
{
	cluster->m_dimpulses[0] += impulse * cluster->m_imass;
	cluster->m_ndimpulses++;
}

struct NodeLinks
{
	btAlignedObjectArray<int> m_links;
};

//
int btSoftBody::generateBendingConstraints(int distance, Material* mat)
{
	int i, j;

	if (distance > 1)
	{
		/* Build graph	*/
		const int n = m_nodes.size();
		const unsigned inf = (~(unsigned)0) >> 1;
		unsigned* adj = new unsigned[n * n];

#define IDX(_x_, _y_) ((_y_) * n + (_x_))
		for (j = 0; j < n; ++j)
		{
			for (i = 0; i < n; ++i)
			{
				if (i != j)
				{
					adj[IDX(i, j)] = adj[IDX(j, i)] = inf;
				}
				else
				{
					adj[IDX(i, j)] = adj[IDX(j, i)] = 0;
				}
			}
		}
		for (i = 0; i < m_links.size(); ++i)
		{
			const int ia = (int)(m_links[i].m_n[0] - &m_nodes[0]);
			const int ib = (int)(m_links[i].m_n[1] - &m_nodes[0]);
			adj[IDX(ia, ib)] = 1;
			adj[IDX(ib, ia)] = 1;
		}

		//special optimized case for distance == 2
		if (distance == 2)
		{
			btAlignedObjectArray<NodeLinks> nodeLinks;

			/* Build node links */
			nodeLinks.resize(m_nodes.size());

			for (i = 0; i < m_links.size(); ++i)
			{
				const int ia = (int)(m_links[i].m_n[0] - &m_nodes[0]);
				const int ib = (int)(m_links[i].m_n[1] - &m_nodes[0]);
				if (nodeLinks[ia].m_links.findLinearSearch(ib) == nodeLinks[ia].m_links.size())
					nodeLinks[ia].m_links.push_back(ib);

				if (nodeLinks[ib].m_links.findLinearSearch(ia) == nodeLinks[ib].m_links.size())
					nodeLinks[ib].m_links.push_back(ia);
			}
			for (int ii = 0; ii < nodeLinks.size(); ii++)
			{
				int i = ii;

				for (int jj = 0; jj < nodeLinks[ii].m_links.size(); jj++)
				{
					int k = nodeLinks[ii].m_links[jj];
					for (int kk = 0; kk < nodeLinks[k].m_links.size(); kk++)
					{
						int j = nodeLinks[k].m_links[kk];
						if (i != j)
						{
							const unsigned sum = adj[IDX(i, k)] + adj[IDX(k, j)];
							btAssert(sum == 2);
							if (adj[IDX(i, j)] > sum)
							{
								adj[IDX(i, j)] = adj[IDX(j, i)] = sum;
							}
						}
					}
				}
			}
		}
		else
		{
			///generic Floyd's algorithm
			for (int k = 0; k < n; ++k)
			{
				for (j = 0; j < n; ++j)
				{
					for (i = j + 1; i < n; ++i)
					{
						const unsigned sum = adj[IDX(i, k)] + adj[IDX(k, j)];
						if (adj[IDX(i, j)] > sum)
						{
							adj[IDX(i, j)] = adj[IDX(j, i)] = sum;
						}
					}
				}
			}
		}

		/* Build links	*/
		int nlinks = 0;
		for (j = 0; j < n; ++j)
		{
			for (i = j + 1; i < n; ++i)
			{
				if (adj[IDX(i, j)] == (unsigned)distance)
				{
					appendLink(i, j, mat);
					m_links[m_links.size() - 1].m_bbending = 1;
					++nlinks;
				}
			}
		}
		delete[] adj;
		return (nlinks);
	}
	return (0);
}

//
void btSoftBody::randomizeConstraints()
{
	unsigned long seed = 243703;
#define NEXTRAND (seed = (1664525L * seed + 1013904223L) & 0xffffffff)
	int i, ni;

	for (i = 0, ni = m_links.size(); i < ni; ++i)
	{
		btSwap(m_links[i], m_links[NEXTRAND % ni]);
	}
	for (i = 0, ni = m_faces.size(); i < ni; ++i)
	{
		btSwap(m_faces[i], m_faces[NEXTRAND % ni]);
	}
#undef NEXTRAND
}

void btSoftBody::updateState(const btAlignedObjectArray<btVector3>& q, const btAlignedObjectArray<btVector3>& v)
{
	int node_count = m_nodes.size();
	btAssert(node_count == q.size());
	btAssert(node_count == v.size());
	for (int i = 0; i < node_count; i++)
	{
		Node& n = m_nodes[i];
		n.m_x = q[i];
		n.m_q = q[i];
		n.m_v = v[i];
		n.m_vn = v[i];
	}
}

//
void btSoftBody::releaseCluster(int index)
{
	Cluster* c = m_clusters[index];
	if (c->m_leaf) m_cdbvt.remove(c->m_leaf);
	c->~Cluster();
	btAlignedFree(c);
	m_clusters.remove(c);
}

//
void btSoftBody::releaseClusters()
{
	while (m_clusters.size() > 0) releaseCluster(0);
}

//
int btSoftBody::generateClusters(int k, int maxiterations)
{
	int i;
	releaseClusters();
	m_clusters.resize(btMin(k, m_nodes.size()));
	for (i = 0; i < m_clusters.size(); ++i)
	{
		m_clusters[i] = new (btAlignedAlloc(sizeof(Cluster), 16)) Cluster();
		m_clusters[i]->m_collide = true;
	}
	k = m_clusters.size();
	if (k > 0)
	{
		/* Initialize		*/
		btAlignedObjectArray<btVector3> centers;
		btVector3 cog(0, 0, 0);
		int i;
		for (i = 0; i < m_nodes.size(); ++i)
		{
			cog += m_nodes[i].m_x;
			m_clusters[(i * 29873) % m_clusters.size()]->m_nodes.push_back(&m_nodes[i]);
		}
		cog /= (btScalar)m_nodes.size();
		centers.resize(k, cog);
		/* Iterate			*/
		const btScalar slope = 16;
		bool changed;
		int iterations = 0;
		do
		{
			const btScalar w = 2 - btMin<btScalar>(1, iterations / slope);
			changed = false;
			iterations++;
			int i;

			for (i = 0; i < k; ++i)
			{
				btVector3 c(0, 0, 0);
				for (int j = 0; j < m_clusters[i]->m_nodes.size(); ++j)
				{
					c += m_clusters[i]->m_nodes[j]->m_x;
				}
				if (m_clusters[i]->m_nodes.size())
				{
					c /= (btScalar)m_clusters[i]->m_nodes.size();
					c = centers[i] + (c - centers[i]) * w;
					changed |= ((c - centers[i]).length2() > SIMD_EPSILON);
					centers[i] = c;
					m_clusters[i]->m_nodes.resize(0);
				}
			}
			for (i = 0; i < m_nodes.size(); ++i)
			{
				const btVector3 nx = m_nodes[i].m_x;
				int kbest = 0;
				btScalar kdist = ClusterMetric(centers[0], nx);
				for (int j = 1; j < k; ++j)
				{
					const btScalar d = ClusterMetric(centers[j], nx);
					if (d < kdist)
					{
						kbest = j;
						kdist = d;
					}
				}
				m_clusters[kbest]->m_nodes.push_back(&m_nodes[i]);
			}
		} while (changed && (iterations < maxiterations));
		/* Merge		*/
		btAlignedObjectArray<int> cids;
		cids.resize(m_nodes.size(), -1);
		for (i = 0; i < m_clusters.size(); ++i)
		{
			for (int j = 0; j < m_clusters[i]->m_nodes.size(); ++j)
			{
				cids[int(m_clusters[i]->m_nodes[j] - &m_nodes[0])] = i;
			}
		}
		for (i = 0; i < m_faces.size(); ++i)
		{
			const int idx[] = {int(m_faces[i].m_n[0] - &m_nodes[0]),
							   int(m_faces[i].m_n[1] - &m_nodes[0]),
							   int(m_faces[i].m_n[2] - &m_nodes[0])};
			for (int j = 0; j < 3; ++j)
			{
				const int cid = cids[idx[j]];
				for (int q = 1; q < 3; ++q)
				{
					const int kid = idx[(j + q) % 3];
					if (cids[kid] != cid)
					{
						if (m_clusters[cid]->m_nodes.findLinearSearch(&m_nodes[kid]) == m_clusters[cid]->m_nodes.size())
						{
							m_clusters[cid]->m_nodes.push_back(&m_nodes[kid]);
						}
					}
				}
			}
		}
		/* Master		*/
		if (m_clusters.size() > 1)
		{
			Cluster* pmaster = new (btAlignedAlloc(sizeof(Cluster), 16)) Cluster();
			pmaster->m_collide = false;
			pmaster->m_nodes.reserve(m_nodes.size());
			for (int i = 0; i < m_nodes.size(); ++i) pmaster->m_nodes.push_back(&m_nodes[i]);
			m_clusters.push_back(pmaster);
			btSwap(m_clusters[0], m_clusters[m_clusters.size() - 1]);
		}
		/* Terminate	*/
		for (i = 0; i < m_clusters.size(); ++i)
		{
			if (m_clusters[i]->m_nodes.size() == 0)
			{
				releaseCluster(i--);
			}
		}
	}
	else
	{
		//create a cluster for each tetrahedron (if tetrahedra exist) or each face
		if (m_tetras.size())
		{
			m_clusters.resize(m_tetras.size());
			for (i = 0; i < m_clusters.size(); ++i)
			{
				m_clusters[i] = new (btAlignedAlloc(sizeof(Cluster), 16)) Cluster();
				m_clusters[i]->m_collide = true;
			}
			for (i = 0; i < m_tetras.size(); i++)
			{
				for (int j = 0; j < 4; j++)
				{
					m_clusters[i]->m_nodes.push_back(m_tetras[i].m_n[j]);
				}
			}
		}
		else
		{
			m_clusters.resize(m_faces.size());
			for (i = 0; i < m_clusters.size(); ++i)
			{
				m_clusters[i] = new (btAlignedAlloc(sizeof(Cluster), 16)) Cluster();
				m_clusters[i]->m_collide = true;
			}

			for (i = 0; i < m_faces.size(); ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					m_clusters[i]->m_nodes.push_back(m_faces[i].m_n[j]);
				}
			}
		}
	}

	if (m_clusters.size())
	{
		initializeClusters();
		updateClusters();

		//for self-collision
		m_clusterConnectivity.resize(m_clusters.size() * m_clusters.size());
		{
			for (int c0 = 0; c0 < m_clusters.size(); c0++)
			{
				m_clusters[c0]->m_clusterIndex = c0;
				for (int c1 = 0; c1 < m_clusters.size(); c1++)
				{
					bool connected = false;
					Cluster* cla = m_clusters[c0];
					Cluster* clb = m_clusters[c1];
					for (int i = 0; !connected && i < cla->m_nodes.size(); i++)
					{
						for (int j = 0; j < clb->m_nodes.size(); j++)
						{
							if (cla->m_nodes[i] == clb->m_nodes[j])
							{
								connected = true;
								break;
							}
						}
					}
					m_clusterConnectivity[c0 + c1 * m_clusters.size()] = connected;
				}
			}
		}
	}

	return (m_clusters.size());
}

//
void btSoftBody::refine(ImplicitFn* ifn, btScalar accurary, bool cut)
{
	const Node* nbase = &m_nodes[0];
	int ncount = m_nodes.size();
	btSymMatrix<int> edges(ncount, -2);
	int newnodes = 0;
	int i, j, k, ni;

	/* Filter out		*/
	for (i = 0; i < m_links.size(); ++i)
	{
		Link& l = m_links[i];
		if (l.m_bbending)
		{
			if (!SameSign(ifn->Eval(l.m_n[0]->m_x), ifn->Eval(l.m_n[1]->m_x)))
			{
				btSwap(m_links[i], m_links[m_links.size() - 1]);
				m_links.pop_back();
				--i;
			}
		}
	}
	/* Fill edges		*/
	for (i = 0; i < m_links.size(); ++i)
	{
		Link& l = m_links[i];
		edges(int(l.m_n[0] - nbase), int(l.m_n[1] - nbase)) = -1;
	}
	for (i = 0; i < m_faces.size(); ++i)
	{
		Face& f = m_faces[i];
		edges(int(f.m_n[0] - nbase), int(f.m_n[1] - nbase)) = -1;
		edges(int(f.m_n[1] - nbase), int(f.m_n[2] - nbase)) = -1;
		edges(int(f.m_n[2] - nbase), int(f.m_n[0] - nbase)) = -1;
	}
	/* Intersect		*/
	for (i = 0; i < ncount; ++i)
	{
		for (j = i + 1; j < ncount; ++j)
		{
			if (edges(i, j) == -1)
			{
				Node& a = m_nodes[i];
				Node& b = m_nodes[j];
				const btScalar t = ImplicitSolve(ifn, a.m_x, b.m_x, accurary);
				if (t > 0)
				{
					const btVector3 x = Lerp(a.m_x, b.m_x, t);
					const btVector3 v = Lerp(a.m_v, b.m_v, t);
					btScalar m = 0;
					if (a.m_im > 0)
					{
						if (b.m_im > 0)
						{
							const btScalar ma = 1 / a.m_im;
							const btScalar mb = 1 / b.m_im;
							const btScalar mc = Lerp(ma, mb, t);
							const btScalar f = (ma + mb) / (ma + mb + mc);
							a.m_im = 1 / (ma * f);
							b.m_im = 1 / (mb * f);
							m = mc * f;
						}
						else
						{
							a.m_im /= 0.5f;
							m = 1 / a.m_im;
						}
					}
					else
					{
						if (b.m_im > 0)
						{
							b.m_im /= 0.5f;
							m = 1 / b.m_im;
						}
						else
							m = 0;
					}
					appendNode(x, m);
					edges(i, j) = m_nodes.size() - 1;
					m_nodes[edges(i, j)].m_v = v;
					++newnodes;
				}
			}
		}
	}
	nbase = &m_nodes[0];
	/* Refine links		*/
	for (i = 0, ni = m_links.size(); i < ni; ++i)
	{
		Link& feat = m_links[i];
		const int idx[] = {int(feat.m_n[0] - nbase),
						   int(feat.m_n[1] - nbase)};
		if ((idx[0] < ncount) && (idx[1] < ncount))
		{
			const int ni = edges(idx[0], idx[1]);
			if (ni > 0)
			{
				appendLink(i);
				Link* pft[] = {&m_links[i],
							   &m_links[m_links.size() - 1]};
				pft[0]->m_n[0] = &m_nodes[idx[0]];
				pft[0]->m_n[1] = &m_nodes[ni];
				pft[1]->m_n[0] = &m_nodes[ni];
				pft[1]->m_n[1] = &m_nodes[idx[1]];
			}
		}
	}
	/* Refine faces		*/
	for (i = 0; i < m_faces.size(); ++i)
	{
		const Face& feat = m_faces[i];
		const int idx[] = {int(feat.m_n[0] - nbase),
						   int(feat.m_n[1] - nbase),
						   int(feat.m_n[2] - nbase)};
		for (j = 2, k = 0; k < 3; j = k++)
		{
			if ((idx[j] < ncount) && (idx[k] < ncount))
			{
				const int ni = edges(idx[j], idx[k]);
				if (ni > 0)
				{
					appendFace(i);
					const int l = (k + 1) % 3;
					Face* pft[] = {&m_faces[i],
								   &m_faces[m_faces.size() - 1]};
					pft[0]->m_n[0] = &m_nodes[idx[l]];
					pft[0]->m_n[1] = &m_nodes[idx[j]];
					pft[0]->m_n[2] = &m_nodes[ni];
					pft[1]->m_n[0] = &m_nodes[ni];
					pft[1]->m_n[1] = &m_nodes[idx[k]];
					pft[1]->m_n[2] = &m_nodes[idx[l]];
					appendLink(ni, idx[l], pft[0]->m_material);
					--i;
					break;
				}
			}
		}
	}
	/* Cut				*/
	if (cut)
	{
		btAlignedObjectArray<int> cnodes;
		const int pcount = ncount;
		int i;
		ncount = m_nodes.size();
		cnodes.resize(ncount, 0);
		/* Nodes		*/
		for (i = 0; i < ncount; ++i)
		{
			const btVector3 x = m_nodes[i].m_x;
			if ((i >= pcount) || (btFabs(ifn->Eval(x)) < accurary))
			{
				const btVector3 v = m_nodes[i].m_v;
				btScalar m = getMass(i);
				if (m > 0)
				{
					m *= 0.5f;
					m_nodes[i].m_im /= 0.5f;
				}
				appendNode(x, m);
				cnodes[i] = m_nodes.size() - 1;
				m_nodes[cnodes[i]].m_v = v;
			}
		}
		nbase = &m_nodes[0];
		/* Links		*/
		for (i = 0, ni = m_links.size(); i < ni; ++i)
		{
			const int id[] = {int(m_links[i].m_n[0] - nbase),
							  int(m_links[i].m_n[1] - nbase)};
			int todetach = 0;
			if (cnodes[id[0]] && cnodes[id[1]])
			{
				appendLink(i);
				todetach = m_links.size() - 1;
			}
			else
			{
				if (((ifn->Eval(m_nodes[id[0]].m_x) < accurary) &&
					 (ifn->Eval(m_nodes[id[1]].m_x) < accurary)))
					todetach = i;
			}
			if (todetach)
			{
				Link& l = m_links[todetach];
				for (int j = 0; j < 2; ++j)
				{
					int cn = cnodes[int(l.m_n[j] - nbase)];
					if (cn) l.m_n[j] = &m_nodes[cn];
				}
			}
		}
		/* Faces		*/
		for (i = 0, ni = m_faces.size(); i < ni; ++i)
		{
			Node** n = m_faces[i].m_n;
			if ((ifn->Eval(n[0]->m_x) < accurary) &&
				(ifn->Eval(n[1]->m_x) < accurary) &&
				(ifn->Eval(n[2]->m_x) < accurary))
			{
				for (int j = 0; j < 3; ++j)
				{
					int cn = cnodes[int(n[j] - nbase)];
					if (cn) n[j] = &m_nodes[cn];
				}
			}
		}
		/* Clean orphans	*/
		int nnodes = m_nodes.size();
		btAlignedObjectArray<int> ranks;
		btAlignedObjectArray<int> todelete;
		ranks.resize(nnodes, 0);
		for (i = 0, ni = m_links.size(); i < ni; ++i)
		{
			for (int j = 0; j < 2; ++j) ranks[int(m_links[i].m_n[j] - nbase)]++;
		}
		for (i = 0, ni = m_faces.size(); i < ni; ++i)
		{
			for (int j = 0; j < 3; ++j) ranks[int(m_faces[i].m_n[j] - nbase)]++;
		}
		for (i = 0; i < m_links.size(); ++i)
		{
			const int id[] = {int(m_links[i].m_n[0] - nbase),
							  int(m_links[i].m_n[1] - nbase)};
			const bool sg[] = {ranks[id[0]] == 1,
							   ranks[id[1]] == 1};
			if (sg[0] || sg[1])
			{
				--ranks[id[0]];
				--ranks[id[1]];
				btSwap(m_links[i], m_links[m_links.size() - 1]);
				m_links.pop_back();
				--i;
			}
		}
#if 0	
		for(i=nnodes-1;i>=0;--i)
		{
			if(!ranks[i]) todelete.push_back(i);
		}	
		if(todelete.size())
		{		
			btAlignedObjectArray<int>&	map=ranks;
			for(int i=0;i<nnodes;++i) map[i]=i;
			PointersToIndices(this);
			for(int i=0,ni=todelete.size();i<ni;++i)
			{
				int		j=todelete[i];
				int&	a=map[j];
				int&	b=map[--nnodes];
				m_ndbvt.remove(m_nodes[a].m_leaf);m_nodes[a].m_leaf=0;
				btSwap(m_nodes[a],m_nodes[b]);
				j=a;a=b;b=j;			
			}
			IndicesToPointers(this,&map[0]);
			m_nodes.resize(nnodes);
		}
#endif
	}
	m_bUpdateRtCst = true;
}

//
bool btSoftBody::cutLink(const Node* node0, const Node* node1, btScalar position)
{
	return (cutLink(int(node0 - &m_nodes[0]), int(node1 - &m_nodes[0]), position));
}

//
bool btSoftBody::cutLink(int node0, int node1, btScalar position)
{
	bool done = false;
	int i, ni;
	//	const btVector3	d=m_nodes[node0].m_x-m_nodes[node1].m_x;
	const btVector3 x = Lerp(m_nodes[node0].m_x, m_nodes[node1].m_x, position);
	const btVector3 v = Lerp(m_nodes[node0].m_v, m_nodes[node1].m_v, position);
	const btScalar m = 1;
	appendNode(x, m);
	appendNode(x, m);
	Node* pa = &m_nodes[node0];
	Node* pb = &m_nodes[node1];
	Node* pn[2] = {&m_nodes[m_nodes.size() - 2],
				   &m_nodes[m_nodes.size() - 1]};
	pn[0]->m_v = v;
	pn[1]->m_v = v;
	for (i = 0, ni = m_links.size(); i < ni; ++i)
	{
		const int mtch = MatchEdge(m_links[i].m_n[0], m_links[i].m_n[1], pa, pb);
		if (mtch != -1)
		{
			appendLink(i);
			Link* pft[] = {&m_links[i], &m_links[m_links.size() - 1]};
			pft[0]->m_n[1] = pn[mtch];
			pft[1]->m_n[0] = pn[1 - mtch];
			done = true;
		}
	}
	for (i = 0, ni = m_faces.size(); i < ni; ++i)
	{
		for (int k = 2, l = 0; l < 3; k = l++)
		{
			const int mtch = MatchEdge(m_faces[i].m_n[k], m_faces[i].m_n[l], pa, pb);
			if (mtch != -1)
			{
				appendFace(i);
				Face* pft[] = {&m_faces[i], &m_faces[m_faces.size() - 1]};
				pft[0]->m_n[l] = pn[mtch];
				pft[1]->m_n[k] = pn[1 - mtch];
				appendLink(pn[0], pft[0]->m_n[(l + 1) % 3], pft[0]->m_material, true);
				appendLink(pn[1], pft[0]->m_n[(l + 1) % 3], pft[0]->m_material, true);
			}
		}
	}
	if (!done)
	{
		m_ndbvt.remove(pn[0]->m_leaf);
		m_ndbvt.remove(pn[1]->m_leaf);
		m_nodes.pop_back();
		m_nodes.pop_back();
	}
	return (done);
}

//
bool btSoftBody::rayTest(const btVector3& rayFrom,
						 const btVector3& rayTo,
						 sRayCast& results)
{
	if (m_faces.size() && m_fdbvt.empty())
		initializeFaceTree();

	results.body = this;
	results.fraction = 1.f;
	results.feature = eFeature::None;
	results.index = -1;

	return (rayTest(rayFrom, rayTo, results.fraction, results.feature, results.index, false) != 0);
}

bool btSoftBody::rayFaceTest(const btVector3& rayFrom,
							 const btVector3& rayTo,
							 sRayCast& results)
{
	if (m_faces.size() == 0)
		return false;
	else
	{
		if (m_fdbvt.empty())
			initializeFaceTree();
	}

	results.body = this;
	results.fraction = 1.f;
	results.index = -1;

	return (rayFaceTest(rayFrom, rayTo, results.fraction, results.index) != 0);
}

//
void btSoftBody::setSolver(eSolverPresets::_ preset)
{
	m_cfg.m_vsequence.clear();
	m_cfg.m_psequence.clear();
	m_cfg.m_dsequence.clear();
	switch (preset)
	{
		case eSolverPresets::Positions:
			m_cfg.m_psequence.push_back(ePSolver::Anchors);
			m_cfg.m_psequence.push_back(ePSolver::RContacts);
			m_cfg.m_psequence.push_back(ePSolver::SContacts);
			m_cfg.m_psequence.push_back(ePSolver::Linear);
			break;
		case eSolverPresets::Velocities:
			m_cfg.m_vsequence.push_back(eVSolver::Linear);

			m_cfg.m_psequence.push_back(ePSolver::Anchors);
			m_cfg.m_psequence.push_back(ePSolver::RContacts);
			m_cfg.m_psequence.push_back(ePSolver::SContacts);

			m_cfg.m_dsequence.push_back(ePSolver::Linear);
			break;
	}
}

void btSoftBody::predictMotion(btScalar dt)
{
	int i, ni;

	/* Update                */
	if (m_bUpdateRtCst)
	{
		m_bUpdateRtCst = false;
		updateConstants();
		m_fdbvt.clear();
		if (m_cfg.collisions & fCollision::VF_SS)
		{
			initializeFaceTree();
		}
	}

	/* Prepare                */
	m_sst.sdt = dt * m_cfg.timescale;
	m_sst.isdt = 1 / m_sst.sdt;
	m_sst.velmrg = m_sst.sdt * 3;
	m_sst.radmrg = getCollisionShape()->getMargin();
	m_sst.updmrg = m_sst.radmrg * (btScalar)0.25;
	/* Forces                */
	addVelocity(m_worldInfo->m_gravity * m_sst.sdt);
	applyForces();
	/* Integrate            */
	for (i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		Node& n = m_nodes[i];
		n.m_q = n.m_x;
		btVector3 deltaV = n.m_f * n.m_im * m_sst.sdt;
		{
			btScalar maxDisplacement = m_worldInfo->m_maxDisplacement;
			btScalar clampDeltaV = maxDisplacement / m_sst.sdt;
			for (int c = 0; c < 3; c++)
			{
				if (deltaV[c] > clampDeltaV)
				{
					deltaV[c] = clampDeltaV;
				}
				if (deltaV[c] < -clampDeltaV)
				{
					deltaV[c] = -clampDeltaV;
				}
			}
		}
		n.m_v += deltaV;
		n.m_x += n.m_v * m_sst.sdt;
		n.m_f = btVector3(0, 0, 0);
	}
	/* Clusters                */
	updateClusters();
	/* Bounds                */
	updateBounds();
	/* Nodes                */
	ATTRIBUTE_ALIGNED16(btDbvtVolume)
	vol;
	for (i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		Node& n = m_nodes[i];
		vol = btDbvtVolume::FromCR(n.m_x, m_sst.radmrg);
		m_ndbvt.update(n.m_leaf,
					   vol,
					   n.m_v * m_sst.velmrg,
					   m_sst.updmrg);
	}
	/* Faces                */
	if (!m_fdbvt.empty())
	{
		for (int i = 0; i < m_faces.size(); ++i)
		{
			Face& f = m_faces[i];
			const btVector3 v = (f.m_n[0]->m_v +
								 f.m_n[1]->m_v +
								 f.m_n[2]->m_v) /
								3;
			vol = VolumeOf(f, m_sst.radmrg);
			m_fdbvt.update(f.m_leaf,
						   vol,
						   v * m_sst.velmrg,
						   m_sst.updmrg);
		}
	}
	/* Pose                    */
	updatePose();
	/* Match                */
	if (m_pose.m_bframe && (m_cfg.kMT > 0))
	{
		const btMatrix3x3 posetrs = m_pose.m_rot;
		for (int i = 0, ni = m_nodes.size(); i < ni; ++i)
		{
			Node& n = m_nodes[i];
			if (n.m_im > 0)
			{
				const btVector3 x = posetrs * m_pose.m_pos[i] + m_pose.m_com;
				n.m_x = Lerp(n.m_x, x, m_cfg.kMT);
			}
		}
	}
	/* Clear contacts        */
	m_rcontacts.resize(0);
	m_scontacts.resize(0);
	/* Optimize dbvt's        */
	m_ndbvt.optimizeIncremental(1);
	m_fdbvt.optimizeIncremental(1);
	m_cdbvt.optimizeIncremental(1);
}

//
void btSoftBody::solveConstraints()
{
	/* Apply clusters		*/
	applyClusters(false);
	/* Prepare links		*/

	int i, ni;

	for (i = 0, ni = m_links.size(); i < ni; ++i)
	{
		Link& l = m_links[i];
		l.m_c3 = l.m_n[1]->m_q - l.m_n[0]->m_q;
		l.m_c2 = 1 / (l.m_c3.length2() * l.m_c0);
	}
	/* Prepare anchors		*/
	for (i = 0, ni = m_anchors.size(); i < ni; ++i)
	{
		Anchor& a = m_anchors[i];
		const btVector3 ra = a.m_body->getWorldTransform().getBasis() * a.m_local;
		a.m_c0 = ImpulseMatrix(m_sst.sdt,
							   a.m_node->m_im,
							   a.m_body->getInvMass(),
							   a.m_body->getInvInertiaTensorWorld(),
							   ra);
		a.m_c1 = ra;
		a.m_c2 = m_sst.sdt * a.m_node->m_im;
		a.m_body->activate();
	}
	/* Solve velocities		*/
	if (m_cfg.viterations > 0)
	{
		/* Solve			*/
		for (int isolve = 0; isolve < m_cfg.viterations; ++isolve)
		{
			for (int iseq = 0; iseq < m_cfg.m_vsequence.size(); ++iseq)
			{
				getSolver(m_cfg.m_vsequence[iseq])(this, 1);
			}
		}
		/* Update			*/
		for (i = 0, ni = m_nodes.size(); i < ni; ++i)
		{
			Node& n = m_nodes[i];
			n.m_x = n.m_q + n.m_v * m_sst.sdt;
		}
	}
	/* Solve positions		*/
	if (m_cfg.piterations > 0)
	{
		for (int isolve = 0; isolve < m_cfg.piterations; ++isolve)
		{
			const btScalar ti = isolve / (btScalar)m_cfg.piterations;
			for (int iseq = 0; iseq < m_cfg.m_psequence.size(); ++iseq)
			{
				getSolver(m_cfg.m_psequence[iseq])(this, 1, ti);
			}
		}
		const btScalar vc = m_sst.isdt * (1 - m_cfg.kDP);
		for (i = 0, ni = m_nodes.size(); i < ni; ++i)
		{
			Node& n = m_nodes[i];
			n.m_v = (n.m_x - n.m_q) * vc;
			n.m_f = btVector3(0, 0, 0);
		}
	}
	/* Solve drift			*/
	if (m_cfg.diterations > 0)
	{
		const btScalar vcf = m_cfg.kVCF * m_sst.isdt;
		for (i = 0, ni = m_nodes.size(); i < ni; ++i)
		{
			Node& n = m_nodes[i];
			n.m_q = n.m_x;
		}
		for (int idrift = 0; idrift < m_cfg.diterations; ++idrift)
		{
			for (int iseq = 0; iseq < m_cfg.m_dsequence.size(); ++iseq)
			{
				getSolver(m_cfg.m_dsequence[iseq])(this, 1, 0);
			}
		}
		for (int i = 0, ni = m_nodes.size(); i < ni; ++i)
		{
			Node& n = m_nodes[i];
			n.m_v += (n.m_x - n.m_q) * vcf;
		}
	}
	/* Apply clusters		*/
	dampClusters();
	applyClusters(true);
}

//
void btSoftBody::staticSolve(int iterations)
{
	for (int isolve = 0; isolve < iterations; ++isolve)
	{
		for (int iseq = 0; iseq < m_cfg.m_psequence.size(); ++iseq)
		{
			getSolver(m_cfg.m_psequence[iseq])(this, 1, 0);
		}
	}
}

//
void btSoftBody::solveCommonConstraints(btSoftBody** /*bodies*/, int /*count*/, int /*iterations*/)
{
	/// placeholder
}

//
void btSoftBody::solveClusters(const btAlignedObjectArray<btSoftBody*>& bodies)
{
	const int nb = bodies.size();
	int iterations = 0;
	int i;

	for (i = 0; i < nb; ++i)
	{
		iterations = btMax(iterations, bodies[i]->m_cfg.citerations);
	}
	for (i = 0; i < nb; ++i)
	{
		bodies[i]->prepareClusters(iterations);
	}
	for (i = 0; i < iterations; ++i)
	{
		const btScalar sor = 1;
		for (int j = 0; j < nb; ++j)
		{
			bodies[j]->solveClusters(sor);
		}
	}
	for (i = 0; i < nb; ++i)
	{
		bodies[i]->cleanupClusters();
	}
}

//
void btSoftBody::integrateMotion()
{
	/* Update			*/
	updateNormals();
}

//
btSoftBody::RayFromToCaster::RayFromToCaster(const btVector3& rayFrom, const btVector3& rayTo, btScalar mxt)
{
	m_rayFrom = rayFrom;
	m_rayNormalizedDirection = (rayTo - rayFrom);
	m_rayTo = rayTo;
	m_mint = mxt;
	m_face = 0;
	m_tests = 0;
}

//
void btSoftBody::RayFromToCaster::Process(const btDbvtNode* leaf)
{
	btSoftBody::Face& f = *(btSoftBody::Face*)leaf->data;
	const btScalar t = rayFromToTriangle(m_rayFrom, m_rayTo, m_rayNormalizedDirection,
										 f.m_n[0]->m_x,
										 f.m_n[1]->m_x,
										 f.m_n[2]->m_x,
										 m_mint);
	if ((t > 0) && (t < m_mint))
	{
		m_mint = t;
		m_face = &f;
	}
	++m_tests;
}

//
btScalar btSoftBody::RayFromToCaster::rayFromToTriangle(const btVector3& rayFrom,
														const btVector3& rayTo,
														const btVector3& rayNormalizedDirection,
														const btVector3& a,
														const btVector3& b,
														const btVector3& c,
														btScalar maxt)
{
	static const btScalar ceps = -SIMD_EPSILON * 10;
	static const btScalar teps = SIMD_EPSILON * 10;

	const btVector3 n = btCross(b - a, c - a);
	const btScalar d = btDot(a, n);
	const btScalar den = btDot(rayNormalizedDirection, n);
	if (!btFuzzyZero(den))
	{
		const btScalar num = btDot(rayFrom, n) - d;
		const btScalar t = -num / den;
		if ((t > teps) && (t < maxt))
		{
			const btVector3 hit = rayFrom + rayNormalizedDirection * t;
			if ((btDot(n, btCross(a - hit, b - hit)) > ceps) &&
				(btDot(n, btCross(b - hit, c - hit)) > ceps) &&
				(btDot(n, btCross(c - hit, a - hit)) > ceps))
			{
				return (t);
			}
		}
	}
	return (-1);
}

//
void btSoftBody::pointersToIndices()
{
#define PTR2IDX(_p_, _b_) reinterpret_cast<btSoftBody::Node*>((_p_) - (_b_))
	btSoftBody::Node* base = m_nodes.size() ? &m_nodes[0] : 0;
	int i, ni;

	for (i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		if (m_nodes[i].m_leaf)
		{
			m_nodes[i].m_leaf->data = *(void**)&i;
		}
	}
	for (i = 0, ni = m_links.size(); i < ni; ++i)
	{
		m_links[i].m_n[0] = PTR2IDX(m_links[i].m_n[0], base);
		m_links[i].m_n[1] = PTR2IDX(m_links[i].m_n[1], base);
	}
	for (i = 0, ni = m_faces.size(); i < ni; ++i)
	{
		m_faces[i].m_n[0] = PTR2IDX(m_faces[i].m_n[0], base);
		m_faces[i].m_n[1] = PTR2IDX(m_faces[i].m_n[1], base);
		m_faces[i].m_n[2] = PTR2IDX(m_faces[i].m_n[2], base);
		if (m_faces[i].m_leaf)
		{
			m_faces[i].m_leaf->data = *(void**)&i;
		}
	}
	for (i = 0, ni = m_anchors.size(); i < ni; ++i)
	{
		m_anchors[i].m_node = PTR2IDX(m_anchors[i].m_node, base);
	}
	for (i = 0, ni = m_notes.size(); i < ni; ++i)
	{
		for (int j = 0; j < m_notes[i].m_rank; ++j)
		{
			m_notes[i].m_nodes[j] = PTR2IDX(m_notes[i].m_nodes[j], base);
		}
	}
#undef PTR2IDX
}

//
void btSoftBody::indicesToPointers(const int* map)
{
#define IDX2PTR(_p_, _b_) map ? (&(_b_)[map[(((char*)_p_) - (char*)0)]]) : (&(_b_)[(((char*)_p_) - (char*)0)])
	btSoftBody::Node* base = m_nodes.size() ? &m_nodes[0] : 0;
	int i, ni;

	for (i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		if (m_nodes[i].m_leaf)
		{
			m_nodes[i].m_leaf->data = &m_nodes[i];
		}
	}
	for (i = 0, ni = m_links.size(); i < ni; ++i)
	{
		m_links[i].m_n[0] = IDX2PTR(m_links[i].m_n[0], base);
		m_links[i].m_n[1] = IDX2PTR(m_links[i].m_n[1], base);
	}
	for (i = 0, ni = m_faces.size(); i < ni; ++i)
	{
		m_faces[i].m_n[0] = IDX2PTR(m_faces[i].m_n[0], base);
		m_faces[i].m_n[1] = IDX2PTR(m_faces[i].m_n[1], base);
		m_faces[i].m_n[2] = IDX2PTR(m_faces[i].m_n[2], base);
		if (m_faces[i].m_leaf)
		{
			m_faces[i].m_leaf->data = &m_faces[i];
		}
	}
	for (i = 0, ni = m_anchors.size(); i < ni; ++i)
	{
		m_anchors[i].m_node = IDX2PTR(m_anchors[i].m_node, base);
	}
	for (i = 0, ni = m_notes.size(); i < ni; ++i)
	{
		for (int j = 0; j < m_notes[i].m_rank; ++j)
		{
			m_notes[i].m_nodes[j] = IDX2PTR(m_notes[i].m_nodes[j], base);
		}
	}
#undef IDX2PTR
}

//
int btSoftBody::rayTest(const btVector3& rayFrom, const btVector3& rayTo,
						btScalar& mint, eFeature::_& feature, int& index, bool bcountonly) const
{
	int cnt = 0;
	btVector3 dir = rayTo - rayFrom;

	if (bcountonly || m_fdbvt.empty())
	{ /* Full search	*/

		for (int i = 0, ni = m_faces.size(); i < ni; ++i)
		{
			const btSoftBody::Face& f = m_faces[i];

			const btScalar t = RayFromToCaster::rayFromToTriangle(rayFrom, rayTo, dir,
																  f.m_n[0]->m_x,
																  f.m_n[1]->m_x,
																  f.m_n[2]->m_x,
																  mint);
			if (t > 0)
			{
				++cnt;
				if (!bcountonly)
				{
					feature = btSoftBody::eFeature::Face;
					index = i;
					mint = t;
				}
			}
		}
	}
	else
	{ /* Use dbvt	*/
		RayFromToCaster collider(rayFrom, rayTo, mint);

		btDbvt::rayTest(m_fdbvt.m_root, rayFrom, rayTo, collider);
		if (collider.m_face)
		{
			mint = collider.m_mint;
			feature = btSoftBody::eFeature::Face;
			index = (int)(collider.m_face - &m_faces[0]);
			cnt = 1;
		}
	}

	for (int i = 0; i < m_tetras.size(); i++)
	{
		const btSoftBody::Tetra& tet = m_tetras[i];
		int tetfaces[4][3] = {{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}};
		for (int f = 0; f < 4; f++)
		{
			int index0 = tetfaces[f][0];
			int index1 = tetfaces[f][1];
			int index2 = tetfaces[f][2];
			btVector3 v0 = tet.m_n[index0]->m_x;
			btVector3 v1 = tet.m_n[index1]->m_x;
			btVector3 v2 = tet.m_n[index2]->m_x;

			const btScalar t = RayFromToCaster::rayFromToTriangle(rayFrom, rayTo, dir,
																  v0, v1, v2,
																  mint);
			if (t > 0)
			{
				++cnt;
				if (!bcountonly)
				{
					feature = btSoftBody::eFeature::Tetra;
					index = i;
					mint = t;
				}
			}
		}
	}
	return (cnt);
}

int btSoftBody::rayFaceTest(const btVector3& rayFrom, const btVector3& rayTo,
							btScalar& mint, int& index) const
{
	int cnt = 0;
	{ /* Use dbvt	*/
		RayFromToCaster collider(rayFrom, rayTo, mint);

		btDbvt::rayTest(m_fdbvt.m_root, rayFrom, rayTo, collider);
		if (collider.m_face)
		{
			mint = collider.m_mint;
			index = (int)(collider.m_face - &m_faces[0]);
			cnt = 1;
		}
	}
	return (cnt);
}

//
static inline btDbvntNode* copyToDbvnt(const btDbvtNode* n)
{
	if (n == 0)
		return 0;
	btDbvntNode* root = new btDbvntNode(n);
	if (n->isinternal())
	{
		btDbvntNode* c0 = copyToDbvnt(n->childs[0]);
		root->childs[0] = c0;
		btDbvntNode* c1 = copyToDbvnt(n->childs[1]);
		root->childs[1] = c1;
	}
	return root;
}

static inline void calculateNormalCone(btDbvntNode* root)
{
	if (!root)
		return;
	if (root->isleaf())
	{
		const btSoftBody::Face* face = (btSoftBody::Face*)root->data;
		root->normal = face->m_normal;
		root->angle = 0;
	}
	else
	{
		btVector3 n0(0, 0, 0), n1(0, 0, 0);
		btScalar a0 = 0, a1 = 0;
		if (root->childs[0])
		{
			calculateNormalCone(root->childs[0]);
			n0 = root->childs[0]->normal;
			a0 = root->childs[0]->angle;
		}
		if (root->childs[1])
		{
			calculateNormalCone(root->childs[1]);
			n1 = root->childs[1]->normal;
			a1 = root->childs[1]->angle;
		}
		root->normal = (n0 + n1).safeNormalize();
		root->angle = btMax(a0, a1) + btAngle(n0, n1) * 0.5;
	}
}

void btSoftBody::initializeFaceTree()
{
	BT_PROFILE("btSoftBody::initializeFaceTree");
	m_fdbvt.clear();
	// create leaf nodes;
	btAlignedObjectArray<btDbvtNode*> leafNodes;
	leafNodes.resize(m_faces.size());
	for (int i = 0; i < m_faces.size(); ++i)
	{
		Face& f = m_faces[i];
		ATTRIBUTE_ALIGNED16(btDbvtVolume)
		vol = VolumeOf(f, 0);
		btDbvtNode* node = new (btAlignedAlloc(sizeof(btDbvtNode), 16)) btDbvtNode();
		node->parent = NULL;
		node->data = &f;
		node->childs[1] = 0;
		node->volume = vol;
		leafNodes[i] = node;
		f.m_leaf = node;
	}
	btAlignedObjectArray<btAlignedObjectArray<int>> adj;
	adj.resize(m_faces.size());
	// construct the adjacency list for triangles
	for (int i = 0; i < adj.size(); ++i)
	{
		for (int j = i + 1; j < adj.size(); ++j)
		{
			int dup = 0;
			for (int k = 0; k < 3; ++k)
			{
				for (int l = 0; l < 3; ++l)
				{
					if (m_faces[i].m_n[k] == m_faces[j].m_n[l])
					{
						++dup;
						break;
					}
				}
				if (dup == 2)
				{
					adj[i].push_back(j);
					adj[j].push_back(i);
				}
			}
		}
	}
	m_fdbvt.m_root = buildTreeBottomUp(leafNodes, adj);
	if (m_fdbvnt)
		delete m_fdbvnt;
	m_fdbvnt = copyToDbvnt(m_fdbvt.m_root);
	updateFaceTree(false, false);
	rebuildNodeTree();
}

//
void btSoftBody::rebuildNodeTree()
{
	m_ndbvt.clear();
	btAlignedObjectArray<btDbvtNode*> leafNodes;
	leafNodes.resize(m_nodes.size());
	for (int i = 0; i < m_nodes.size(); ++i)
	{
		Node& n = m_nodes[i];
		ATTRIBUTE_ALIGNED16(btDbvtVolume)
		vol = btDbvtVolume::FromCR(n.m_x, 0);
		btDbvtNode* node = new (btAlignedAlloc(sizeof(btDbvtNode), 16)) btDbvtNode();
		node->parent = NULL;
		node->data = &n;
		node->childs[1] = 0;
		node->volume = vol;
		leafNodes[i] = node;
		n.m_leaf = node;
	}
	btAlignedObjectArray<btAlignedObjectArray<int>> adj;
	adj.resize(m_nodes.size());
	btAlignedObjectArray<int> old_id;
	old_id.resize(m_nodes.size());
	for (int i = 0; i < m_nodes.size(); ++i)
		old_id[i] = m_nodes[i].index;
	for (int i = 0; i < m_nodes.size(); ++i)
		m_nodes[i].index = i;
	for (int i = 0; i < m_links.size(); ++i)
	{
		Link& l = m_links[i];
		adj[l.m_n[0]->index].push_back(l.m_n[1]->index);
		adj[l.m_n[1]->index].push_back(l.m_n[0]->index);
	}
	m_ndbvt.m_root = buildTreeBottomUp(leafNodes, adj);
	for (int i = 0; i < m_nodes.size(); ++i)
		m_nodes[i].index = old_id[i];
}

//
btVector3 btSoftBody::evaluateCom() const
{
	btVector3 com(0, 0, 0);
	if (m_pose.m_bframe)
	{
		for (int i = 0, ni = m_nodes.size(); i < ni; ++i)
		{
			com += m_nodes[i].m_x * m_pose.m_wgh[i];
		}
	}
	return (com);
}

bool btSoftBody::checkContact(const btCollisionObjectWrapper* colObjWrap,
							  const btVector3& x,
							  btScalar margin,
							  btSoftBody::sCti& cti) const
{
	btVector3 nrm;
	const btCollisionShape* shp = colObjWrap->getCollisionShape();
	//    const btRigidBody *tmpRigid = btRigidBody::upcast(colObjWrap->getCollisionObject());
	//const btTransform &wtr = tmpRigid ? tmpRigid->getWorldTransform() : colObjWrap->getWorldTransform();
	const btTransform& wtr = colObjWrap->getWorldTransform();
	//todo: check which transform is needed here

	btScalar dst =
		m_worldInfo->m_sparsesdf.Evaluate(
			wtr.invXform(x),
			shp,
			nrm,
			margin);
	if (dst < 0)
	{
		cti.m_colObj = colObjWrap->getCollisionObject();
		cti.m_normal = wtr.getBasis() * nrm;
		cti.m_offset = -btDot(cti.m_normal, x - cti.m_normal * dst);
		return (true);
	}
	return (false);
}

//
bool btSoftBody::checkDeformableContact(const btCollisionObjectWrapper* colObjWrap,
										const btVector3& x,
										btScalar margin,
										btSoftBody::sCti& cti, bool predict) const
{
	btVector3 nrm;
	const btCollisionShape* shp = colObjWrap->getCollisionShape();
	const btCollisionObject* tmpCollisionObj = colObjWrap->getCollisionObject();
	// use the position x_{n+1}^* = x_n + dt * v_{n+1}^* where v_{n+1}^* = v_n + dtg for collision detect
	// but resolve contact at x_n
	btTransform wtr = (predict) ? (colObjWrap->m_preTransform != NULL ? tmpCollisionObj->getInterpolationWorldTransform() * (*colObjWrap->m_preTransform) : tmpCollisionObj->getInterpolationWorldTransform())
								: colObjWrap->getWorldTransform();
	btScalar dst =
		m_worldInfo->m_sparsesdf.Evaluate(
			wtr.invXform(x),
			shp,
			nrm,
			margin);

	if (!predict)
	{
		cti.m_colObj = colObjWrap->getCollisionObject();
		cti.m_normal = wtr.getBasis() * nrm;
		cti.m_offset = dst;
	}
	if (dst < 0)
		return true;
	return (false);
}

//
// Compute barycentric coordinates (u, v, w) for
// point p with respect to triangle (a, b, c)
static void getBarycentric(const btVector3& p, const btVector3& a, const btVector3& b, const btVector3& c, btVector3& bary)
{
	btVector3 v0 = b - a, v1 = c - a, v2 = p - a;
	btScalar d00 = v0.dot(v0);
	btScalar d01 = v0.dot(v1);
	btScalar d11 = v1.dot(v1);
	btScalar d20 = v2.dot(v0);
	btScalar d21 = v2.dot(v1);
	btScalar denom = d00 * d11 - d01 * d01;
	// In the case of a degenerate triangle, pick a vertex.
	if (btFabs(denom) < SIMD_EPSILON)
	{
		bary.setY(btScalar(0.0));
		bary.setZ(btScalar(0.0));
	}
	else
	{
		bary.setY((d11 * d20 - d01 * d21) / denom);
		bary.setZ((d00 * d21 - d01 * d20) / denom);
	}
	bary.setX(btScalar(1) - bary.getY() - bary.getZ());
}

//
bool btSoftBody::checkDeformableFaceContact(const btCollisionObjectWrapper* colObjWrap,
											Face& f,
											btVector3& contact_point,
											btVector3& bary,
											btScalar margin,
											btSoftBody::sCti& cti, bool predict) const
{
	btVector3 nrm;
	const btCollisionShape* shp = colObjWrap->getCollisionShape();
	const btCollisionObject* tmpCollisionObj = colObjWrap->getCollisionObject();
	// use the position x_{n+1}^* = x_n + dt * v_{n+1}^* where v_{n+1}^* = v_n + dtg for collision detect
	// but resolve contact at x_n
	btTransform wtr = (predict) ? (colObjWrap->m_preTransform != NULL ? tmpCollisionObj->getInterpolationWorldTransform() * (*colObjWrap->m_preTransform) : tmpCollisionObj->getInterpolationWorldTransform())
								: colObjWrap->getWorldTransform();
	btScalar dst;
	btGjkEpaSolver2::sResults results;

	//	#define USE_QUADRATURE 1

	// use collision quadrature point
#ifdef USE_QUADRATURE
	{
		dst = SIMD_INFINITY;
		btVector3 local_nrm;
		for (int q = 0; q < m_quads.size(); ++q)
		{
			btVector3 p;
			if (predict)
				p = BaryEval(f.m_n[0]->m_q, f.m_n[1]->m_q, f.m_n[2]->m_q, m_quads[q]);
			else
				p = BaryEval(f.m_n[0]->m_x, f.m_n[1]->m_x, f.m_n[2]->m_x, m_quads[q]);
			btScalar local_dst = m_worldInfo->m_sparsesdf.Evaluate(
				wtr.invXform(p),
				shp,
				local_nrm,
				margin);
			if (local_dst < dst)
			{
				if (local_dst < 0 && predict)
					return true;
				dst = local_dst;
				contact_point = p;
				bary = m_quads[q];
				nrm = local_nrm;
			}
			if (!predict)
			{
				cti.m_colObj = colObjWrap->getCollisionObject();
				cti.m_normal = wtr.getBasis() * nrm;
				cti.m_offset = dst;
			}
		}
		return (dst < 0);
	}
#endif

	// collision detection using x*
	btTransform triangle_transform;
	triangle_transform.setIdentity();
	triangle_transform.setOrigin(f.m_n[0]->m_q);
	btTriangleShape triangle(btVector3(0, 0, 0), f.m_n[1]->m_q - f.m_n[0]->m_q, f.m_n[2]->m_q - f.m_n[0]->m_q);
	btVector3 guess(0, 0, 0);
	const btConvexShape* csh = static_cast<const btConvexShape*>(shp);
	btGjkEpaSolver2::SignedDistance(&triangle, triangle_transform, csh, wtr, guess, results);
	dst = results.distance - 2.0 * csh->getMargin() - margin;  // margin padding so that the distance = the actual distance between face and rigid - margin of rigid - margin of deformable
	if (dst >= 0)
		return false;

	// Use consistent barycenter to recalculate distance.
	if (this->m_cacheBarycenter)
	{
		if (f.m_pcontact[3] != 0)
		{
			for (int i = 0; i < 3; ++i)
				bary[i] = f.m_pcontact[i];
			contact_point = BaryEval(f.m_n[0]->m_x, f.m_n[1]->m_x, f.m_n[2]->m_x, bary);
			const btConvexShape* csh = static_cast<const btConvexShape*>(shp);
			btGjkEpaSolver2::SignedDistance(contact_point, margin, csh, wtr, results);
			cti.m_colObj = colObjWrap->getCollisionObject();
			dst = results.distance;
			cti.m_normal = results.normal;
			cti.m_offset = dst;

			//point-convex CD
			wtr = colObjWrap->getWorldTransform();
			btTriangleShape triangle2(btVector3(0, 0, 0), f.m_n[1]->m_x - f.m_n[0]->m_x, f.m_n[2]->m_x - f.m_n[0]->m_x);
			triangle_transform.setOrigin(f.m_n[0]->m_x);
			btGjkEpaSolver2::SignedDistance(&triangle2, triangle_transform, csh, wtr, guess, results);

			dst = results.distance - csh->getMargin() - margin;
			return true;
		}
	}

	// Use triangle-convex CD.
	wtr = colObjWrap->getWorldTransform();
	btTriangleShape triangle2(btVector3(0, 0, 0), f.m_n[1]->m_x - f.m_n[0]->m_x, f.m_n[2]->m_x - f.m_n[0]->m_x);
	triangle_transform.setOrigin(f.m_n[0]->m_x);
	btGjkEpaSolver2::SignedDistance(&triangle2, triangle_transform, csh, wtr, guess, results);
	contact_point = results.witnesses[0];
	getBarycentric(contact_point, f.m_n[0]->m_x, f.m_n[1]->m_x, f.m_n[2]->m_x, bary);

	for (int i = 0; i < 3; ++i)
		f.m_pcontact[i] = bary[i];

	dst = results.distance - csh->getMargin() - margin;
	cti.m_colObj = colObjWrap->getCollisionObject();
	cti.m_normal = results.normal;
	cti.m_offset = dst;
	return true;
}

void btSoftBody::updateNormals()
{
	const btVector3 zv(0, 0, 0);
	int i, ni;

	for (i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		m_nodes[i].m_n = zv;
	}
	for (i = 0, ni = m_faces.size(); i < ni; ++i)
	{
		btSoftBody::Face& f = m_faces[i];
		const btVector3 n = btCross(f.m_n[1]->m_x - f.m_n[0]->m_x,
									f.m_n[2]->m_x - f.m_n[0]->m_x);
		f.m_normal = n;
		f.m_normal.safeNormalize();
		f.m_n[0]->m_n += n;
		f.m_n[1]->m_n += n;
		f.m_n[2]->m_n += n;
	}
	for (i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		btScalar len = m_nodes[i].m_n.length();
		if (len > SIMD_EPSILON)
			m_nodes[i].m_n /= len;
	}
}

//
void btSoftBody::updateBounds()
{
	/*if( m_acceleratedSoftBody )
	{
		// If we have an accelerated softbody we need to obtain the bounds correctly
		// For now (slightly hackily) just have a very large AABB
		// TODO: Write get bounds kernel
		// If that is updating in place, atomic collisions might be low (when the cloth isn't perfectly aligned to an axis) and we could
		// probably do a test and exchange reasonably efficiently.

		m_bounds[0] = btVector3(-1000, -1000, -1000);
		m_bounds[1] = btVector3(1000, 1000, 1000);

	} else {*/
	//    if (m_ndbvt.m_root)
	//    {
	//        const btVector3& mins = m_ndbvt.m_root->volume.Mins();
	//        const btVector3& maxs = m_ndbvt.m_root->volume.Maxs();
	//        const btScalar csm = getCollisionShape()->getMargin();
	//        const btVector3 mrg = btVector3(csm,
	//                                        csm,
	//                                        csm) *
	//                              1;  // ??? to investigate...
	//        m_bounds[0] = mins - mrg;
	//        m_bounds[1] = maxs + mrg;
	//        if (0 != getBroadphaseHandle())
	//        {
	//            m_worldInfo->m_broadphase->setAabb(getBroadphaseHandle(),
	//                                               m_bounds[0],
	//                                               m_bounds[1],
	//                                               m_worldInfo->m_dispatcher);
	//        }
	//    }
	//    else
	//    {
	//        m_bounds[0] =
	//            m_bounds[1] = btVector3(0, 0, 0);
	//    }
	if (m_nodes.size())
	{
		btVector3 mins = m_nodes[0].m_x;
		btVector3 maxs = m_nodes[0].m_x;
		for (int i = 1; i < m_nodes.size(); ++i)
		{
			for (int d = 0; d < 3; ++d)
			{
				if (m_nodes[i].m_x[d] > maxs[d])
					maxs[d] = m_nodes[i].m_x[d];
				if (m_nodes[i].m_x[d] < mins[d])
					mins[d] = m_nodes[i].m_x[d];
			}
		}
		const btScalar csm = getCollisionShape()->getMargin();
		const btVector3 mrg = btVector3(csm,
										csm,
										csm);
		m_bounds[0] = mins - mrg;
		m_bounds[1] = maxs + mrg;
		if (0 != getBroadphaseHandle())
		{
			m_worldInfo->m_broadphase->setAabb(getBroadphaseHandle(),
											   m_bounds[0],
											   m_bounds[1],
											   m_worldInfo->m_dispatcher);
		}
	}
	else
	{
		m_bounds[0] =
			m_bounds[1] = btVector3(0, 0, 0);
	}
}

//
void btSoftBody::updatePose()
{
	if (m_pose.m_bframe)
	{
		btSoftBody::Pose& pose = m_pose;
		const btVector3 com = evaluateCom();
		/* Com			*/
		pose.m_com = com;
		/* Rotation		*/
		btMatrix3x3 Apq;
		const btScalar eps = SIMD_EPSILON;
		Apq[0] = Apq[1] = Apq[2] = btVector3(0, 0, 0);
		Apq[0].setX(eps);
		Apq[1].setY(eps * 2);
		Apq[2].setZ(eps * 3);
		for (int i = 0, ni = m_nodes.size(); i < ni; ++i)
		{
			const btVector3 a = pose.m_wgh[i] * (m_nodes[i].m_x - com);
			const btVector3& b = pose.m_pos[i];
			Apq[0] += a.x() * b;
			Apq[1] += a.y() * b;
			Apq[2] += a.z() * b;
		}
		btMatrix3x3 r, s;
		PolarDecompose(Apq, r, s);
		pose.m_rot = r;
		pose.m_scl = pose.m_aqq * r.transpose() * Apq;
		if (m_cfg.maxvolume > 1)
		{
			const btScalar idet = Clamp<btScalar>(1 / pose.m_scl.determinant(),
												  1, m_cfg.maxvolume);
			pose.m_scl = Mul(pose.m_scl, idet);
		}
	}
}

//
void btSoftBody::updateArea(bool averageArea)
{
	int i, ni;

	/* Face area		*/
	for (i = 0, ni = m_faces.size(); i < ni; ++i)
	{
		Face& f = m_faces[i];
		f.m_ra = AreaOf(f.m_n[0]->m_x, f.m_n[1]->m_x, f.m_n[2]->m_x);
	}

	/* Node area		*/

	if (averageArea)
	{
		btAlignedObjectArray<int> counts;
		counts.resize(m_nodes.size(), 0);
		for (i = 0, ni = m_nodes.size(); i < ni; ++i)
		{
			m_nodes[i].m_area = 0;
		}
		for (i = 0, ni = m_faces.size(); i < ni; ++i)
		{
			btSoftBody::Face& f = m_faces[i];
			for (int j = 0; j < 3; ++j)
			{
				const int index = (int)(f.m_n[j] - &m_nodes[0]);
				counts[index]++;
				f.m_n[j]->m_area += btFabs(f.m_ra);
			}
		}
		for (i = 0, ni = m_nodes.size(); i < ni; ++i)
		{
			if (counts[i] > 0)
				m_nodes[i].m_area /= (btScalar)counts[i];
			else
				m_nodes[i].m_area = 0;
		}
	}
	else
	{
		// initialize node area as zero
		for (i = 0, ni = m_nodes.size(); i < ni; ++i)
		{
			m_nodes[i].m_area = 0;
		}

		for (i = 0, ni = m_faces.size(); i < ni; ++i)
		{
			btSoftBody::Face& f = m_faces[i];

			for (int j = 0; j < 3; ++j)
			{
				f.m_n[j]->m_area += f.m_ra;
			}
		}

		for (i = 0, ni = m_nodes.size(); i < ni; ++i)
		{
			m_nodes[i].m_area *= 0.3333333f;
		}
	}
}

void btSoftBody::updateLinkConstants()
{
	int i, ni;

	/* Links		*/
	for (i = 0, ni = m_links.size(); i < ni; ++i)
	{
		Link& l = m_links[i];
		Material& m = *l.m_material;
		l.m_c0 = (l.m_n[0]->m_im + l.m_n[1]->m_im) / m.m_kLST;
	}
}

void btSoftBody::updateConstants()
{
	resetLinkRestLengths();
	updateLinkConstants();
	updateArea();
}

//
void btSoftBody::initializeClusters()
{
	int i;

	for (i = 0; i < m_clusters.size(); ++i)
	{
		Cluster& c = *m_clusters[i];
		c.m_imass = 0;
		c.m_masses.resize(c.m_nodes.size());
		for (int j = 0; j < c.m_nodes.size(); ++j)
		{
			if (c.m_nodes[j]->m_im == 0)
			{
				c.m_containsAnchor = true;
				c.m_masses[j] = BT_LARGE_FLOAT;
			}
			else
			{
				c.m_masses[j] = btScalar(1.) / c.m_nodes[j]->m_im;
			}
			c.m_imass += c.m_masses[j];
		}
		c.m_imass = btScalar(1.) / c.m_imass;
		c.m_com = btSoftBody::clusterCom(&c);
		c.m_lv = btVector3(0, 0, 0);
		c.m_av = btVector3(0, 0, 0);
		c.m_leaf = 0;
		/* Inertia	*/
		btMatrix3x3& ii = c.m_locii;
		ii[0] = ii[1] = ii[2] = btVector3(0, 0, 0);
		{
			int i, ni;

			for (i = 0, ni = c.m_nodes.size(); i < ni; ++i)
			{
				const btVector3 k = c.m_nodes[i]->m_x - c.m_com;
				const btVector3 q = k * k;
				const btScalar m = c.m_masses[i];
				ii[0][0] += m * (q[1] + q[2]);
				ii[1][1] += m * (q[0] + q[2]);
				ii[2][2] += m * (q[0] + q[1]);
				ii[0][1] -= m * k[0] * k[1];
				ii[0][2] -= m * k[0] * k[2];
				ii[1][2] -= m * k[1] * k[2];
			}
		}
		ii[1][0] = ii[0][1];
		ii[2][0] = ii[0][2];
		ii[2][1] = ii[1][2];

		ii = ii.inverse();

		/* Frame	*/
		c.m_framexform.setIdentity();
		c.m_framexform.setOrigin(c.m_com);
		c.m_framerefs.resize(c.m_nodes.size());
		{
			int i;
			for (i = 0; i < c.m_framerefs.size(); ++i)
			{
				c.m_framerefs[i] = c.m_nodes[i]->m_x - c.m_com;
			}
		}
	}
}

//
void btSoftBody::updateClusters()
{
	BT_PROFILE("UpdateClusters");
	int i;

	for (i = 0; i < m_clusters.size(); ++i)
	{
		btSoftBody::Cluster& c = *m_clusters[i];
		const int n = c.m_nodes.size();
		//const btScalar			invn=1/(btScalar)n;
		if (n)
		{
			/* Frame				*/
			const btScalar eps = btScalar(0.0001);
			btMatrix3x3 m, r, s;
			m[0] = m[1] = m[2] = btVector3(0, 0, 0);
			m[0][0] = eps * 1;
			m[1][1] = eps * 2;
			m[2][2] = eps * 3;
			c.m_com = clusterCom(&c);
			for (int i = 0; i < c.m_nodes.size(); ++i)
			{
				const btVector3 a = c.m_nodes[i]->m_x - c.m_com;
				const btVector3& b = c.m_framerefs[i];
				m[0] += a[0] * b;
				m[1] += a[1] * b;
				m[2] += a[2] * b;
			}
			PolarDecompose(m, r, s);
			c.m_framexform.setOrigin(c.m_com);
			c.m_framexform.setBasis(r);
			/* Inertia			*/
#if 1 /* Constant	*/
			c.m_invwi = c.m_framexform.getBasis() * c.m_locii * c.m_framexform.getBasis().transpose();
#else
#if 0 /* Sphere	*/ 
			const btScalar	rk=(2*c.m_extents.length2())/(5*c.m_imass);
			const btVector3	inertia(rk,rk,rk);
			const btVector3	iin(btFabs(inertia[0])>SIMD_EPSILON?1/inertia[0]:0,
				btFabs(inertia[1])>SIMD_EPSILON?1/inertia[1]:0,
				btFabs(inertia[2])>SIMD_EPSILON?1/inertia[2]:0);

			c.m_invwi=c.m_xform.getBasis().scaled(iin)*c.m_xform.getBasis().transpose();
#else /* Actual	*/
			c.m_invwi[0] = c.m_invwi[1] = c.m_invwi[2] = btVector3(0, 0, 0);
			for (int i = 0; i < n; ++i)
			{
				const btVector3 k = c.m_nodes[i]->m_x - c.m_com;
				const btVector3 q = k * k;
				const btScalar m = 1 / c.m_nodes[i]->m_im;
				c.m_invwi[0][0] += m * (q[1] + q[2]);
				c.m_invwi[1][1] += m * (q[0] + q[2]);
				c.m_invwi[2][2] += m * (q[0] + q[1]);
				c.m_invwi[0][1] -= m * k[0] * k[1];
				c.m_invwi[0][2] -= m * k[0] * k[2];
				c.m_invwi[1][2] -= m * k[1] * k[2];
			}
			c.m_invwi[1][0] = c.m_invwi[0][1];
			c.m_invwi[2][0] = c.m_invwi[0][2];
			c.m_invwi[2][1] = c.m_invwi[1][2];
			c.m_invwi = c.m_invwi.inverse();
#endif
#endif
			/* Velocities			*/
			c.m_lv = btVector3(0, 0, 0);
			c.m_av = btVector3(0, 0, 0);
			{
				int i;

				for (i = 0; i < n; ++i)
				{
					const btVector3 v = c.m_nodes[i]->m_v * c.m_masses[i];
					c.m_lv += v;
					c.m_av += btCross(c.m_nodes[i]->m_x - c.m_com, v);
				}
			}
			c.m_lv = c.m_imass * c.m_lv * (1 - c.m_ldamping);
			c.m_av = c.m_invwi * c.m_av * (1 - c.m_adamping);
			c.m_vimpulses[0] =
				c.m_vimpulses[1] = btVector3(0, 0, 0);
			c.m_dimpulses[0] =
				c.m_dimpulses[1] = btVector3(0, 0, 0);
			c.m_nvimpulses = 0;
			c.m_ndimpulses = 0;
			/* Matching				*/
			if (c.m_matching > 0)
			{
				for (int j = 0; j < c.m_nodes.size(); ++j)
				{
					Node& n = *c.m_nodes[j];
					const btVector3 x = c.m_framexform * c.m_framerefs[j];
					n.m_x = Lerp(n.m_x, x, c.m_matching);
				}
			}
			/* Dbvt					*/
			if (c.m_collide)
			{
				btVector3 mi = c.m_nodes[0]->m_x;
				btVector3 mx = mi;
				for (int j = 1; j < n; ++j)
				{
					mi.setMin(c.m_nodes[j]->m_x);
					mx.setMax(c.m_nodes[j]->m_x);
				}
				ATTRIBUTE_ALIGNED16(btDbvtVolume)
				bounds = btDbvtVolume::FromMM(mi, mx);
				if (c.m_leaf)
					m_cdbvt.update(c.m_leaf, bounds, c.m_lv * m_sst.sdt * 3, m_sst.radmrg);
				else
					c.m_leaf = m_cdbvt.insert(bounds, &c);
			}
		}
	}
}

//
void btSoftBody::cleanupClusters()
{
	for (int i = 0; i < m_joints.size(); ++i)
	{
		m_joints[i]->Terminate(m_sst.sdt);
		if (m_joints[i]->m_delete)
		{
			btAlignedFree(m_joints[i]);
			m_joints.remove(m_joints[i--]);
		}
	}
}

//
void btSoftBody::prepareClusters(int iterations)
{
	for (int i = 0; i < m_joints.size(); ++i)
	{
		m_joints[i]->Prepare(m_sst.sdt, iterations);
	}
}

//
void btSoftBody::solveClusters(btScalar sor)
{
	for (int i = 0, ni = m_joints.size(); i < ni; ++i)
	{
		m_joints[i]->Solve(m_sst.sdt, sor);
	}
}

//
void btSoftBody::applyClusters(bool drift)
{
	BT_PROFILE("ApplyClusters");
	//	const btScalar					f0=m_sst.sdt;
	//const btScalar					f1=f0/2;
	btAlignedObjectArray<btVector3> deltas;
	btAlignedObjectArray<btScalar> weights;
	deltas.resize(m_nodes.size(), btVector3(0, 0, 0));
	weights.resize(m_nodes.size(), 0);
	int i;

	if (drift)
	{
		for (i = 0; i < m_clusters.size(); ++i)
		{
			Cluster& c = *m_clusters[i];
			if (c.m_ndimpulses)
			{
				c.m_dimpulses[0] /= (btScalar)c.m_ndimpulses;
				c.m_dimpulses[1] /= (btScalar)c.m_ndimpulses;
			}
		}
	}

	for (i = 0; i < m_clusters.size(); ++i)
	{
		Cluster& c = *m_clusters[i];
		if (0 < (drift ? c.m_ndimpulses : c.m_nvimpulses))
		{
			const btVector3 v = (drift ? c.m_dimpulses[0] : c.m_vimpulses[0]) * m_sst.sdt;
			const btVector3 w = (drift ? c.m_dimpulses[1] : c.m_vimpulses[1]) * m_sst.sdt;
			for (int j = 0; j < c.m_nodes.size(); ++j)
			{
				const int idx = int(c.m_nodes[j] - &m_nodes[0]);
				const btVector3& x = c.m_nodes[j]->m_x;
				const btScalar q = c.m_masses[j];
				deltas[idx] += (v + btCross(w, x - c.m_com)) * q;
				weights[idx] += q;
			}
		}
	}
	for (i = 0; i < deltas.size(); ++i)
	{
		if (weights[i] > 0)
		{
			m_nodes[i].m_x += deltas[i] / weights[i];
		}
	}
}

//
void btSoftBody::dampClusters()
{
	int i;

	for (i = 0; i < m_clusters.size(); ++i)
	{
		Cluster& c = *m_clusters[i];
		if (c.m_ndamping > 0)
		{
			for (int j = 0; j < c.m_nodes.size(); ++j)
			{
				Node& n = *c.m_nodes[j];
				if (n.m_im > 0)
				{
					const btVector3 vx = c.m_lv + btCross(c.m_av, c.m_nodes[j]->m_q - c.m_com);
					if (vx.length2() <= n.m_v.length2())
					{
						n.m_v += c.m_ndamping * (vx - n.m_v);
					}
				}
			}
		}
	}
}

void btSoftBody::setSpringStiffness(btScalar k)
{
	for (int i = 0; i < m_links.size(); ++i)
	{
		m_links[i].Feature::m_material->m_kLST = k;
	}
	m_repulsionStiffness = k;
}

void btSoftBody::setGravityFactor(btScalar gravFactor)
{
	m_gravityFactor = gravFactor;
}

void btSoftBody::setCacheBarycenter(bool cacheBarycenter)
{
	m_cacheBarycenter = cacheBarycenter;
}

void btSoftBody::initializeDmInverse()
{
	btScalar unit_simplex_measure = 1. / 6.;

	for (int i = 0; i < m_tetras.size(); ++i)
	{
		Tetra& t = m_tetras[i];
		btVector3 c1 = t.m_n[1]->m_x - t.m_n[0]->m_x;
		btVector3 c2 = t.m_n[2]->m_x - t.m_n[0]->m_x;
		btVector3 c3 = t.m_n[3]->m_x - t.m_n[0]->m_x;
		btMatrix3x3 Dm(c1.getX(), c2.getX(), c3.getX(),
					   c1.getY(), c2.getY(), c3.getY(),
					   c1.getZ(), c2.getZ(), c3.getZ());
		t.m_element_measure = Dm.determinant() * unit_simplex_measure;
		t.m_Dm_inverse = Dm.inverse();

		// calculate the first three columns of P^{-1}
		btVector3 a = t.m_n[0]->m_x;
		btVector3 b = t.m_n[1]->m_x;
		btVector3 c = t.m_n[2]->m_x;
		btVector3 d = t.m_n[3]->m_x;

		btScalar det = 1 / (a[0] * b[1] * c[2] - a[0] * b[1] * d[2] - a[0] * b[2] * c[1] + a[0] * b[2] * d[1] + a[0] * c[1] * d[2] - a[0] * c[2] * d[1] + a[1] * (-b[0] * c[2] + b[0] * d[2] + b[2] * c[0] - b[2] * d[0] - c[0] * d[2] + c[2] * d[0]) + a[2] * (b[0] * c[1] - b[0] * d[1] + b[1] * (d[0] - c[0]) + c[0] * d[1] - c[1] * d[0]) - b[0] * c[1] * d[2] + b[0] * c[2] * d[1] + b[1] * c[0] * d[2] - b[1] * c[2] * d[0] - b[2] * c[0] * d[1] + b[2] * c[1] * d[0]);

		btScalar P11 = -b[2] * c[1] + d[2] * c[1] + b[1] * c[2] + b[2] * d[1] - c[2] * d[1] - b[1] * d[2];
		btScalar P12 = b[2] * c[0] - d[2] * c[0] - b[0] * c[2] - b[2] * d[0] + c[2] * d[0] + b[0] * d[2];
		btScalar P13 = -b[1] * c[0] + d[1] * c[0] + b[0] * c[1] + b[1] * d[0] - c[1] * d[0] - b[0] * d[1];
		btScalar P21 = a[2] * c[1] - d[2] * c[1] - a[1] * c[2] - a[2] * d[1] + c[2] * d[1] + a[1] * d[2];
		btScalar P22 = -a[2] * c[0] + d[2] * c[0] + a[0] * c[2] + a[2] * d[0] - c[2] * d[0] - a[0] * d[2];
		btScalar P23 = a[1] * c[0] - d[1] * c[0] - a[0] * c[1] - a[1] * d[0] + c[1] * d[0] + a[0] * d[1];
		btScalar P31 = -a[2] * b[1] + d[2] * b[1] + a[1] * b[2] + a[2] * d[1] - b[2] * d[1] - a[1] * d[2];
		btScalar P32 = a[2] * b[0] - d[2] * b[0] - a[0] * b[2] - a[2] * d[0] + b[2] * d[0] + a[0] * d[2];
		btScalar P33 = -a[1] * b[0] + d[1] * b[0] + a[0] * b[1] + a[1] * d[0] - b[1] * d[0] - a[0] * d[1];
		btScalar P41 = a[2] * b[1] - c[2] * b[1] - a[1] * b[2] - a[2] * c[1] + b[2] * c[1] + a[1] * c[2];
		btScalar P42 = -a[2] * b[0] + c[2] * b[0] + a[0] * b[2] + a[2] * c[0] - b[2] * c[0] - a[0] * c[2];
		btScalar P43 = a[1] * b[0] - c[1] * b[0] - a[0] * b[1] - a[1] * c[0] + b[1] * c[0] + a[0] * c[1];

		btVector4 p1(P11 * det, P21 * det, P31 * det, P41 * det);
		btVector4 p2(P12 * det, P22 * det, P32 * det, P42 * det);
		btVector4 p3(P13 * det, P23 * det, P33 * det, P43 * det);

		t.m_P_inv[0] = p1;
		t.m_P_inv[1] = p2;
		t.m_P_inv[2] = p3;

		t.m_rv = VolumeOf(a, b, c, d);
	}
}

static btScalar Dot4(const btVector4& a, const btVector4& b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
}

void btSoftBody::updateDeformation()
{
	btQuaternion q;
	for (int i = 0; i < m_tetras.size(); ++i)
	{
		btSoftBody::Tetra& t = m_tetras[i];
		btVector3 c1 = t.m_n[1]->m_q - t.m_n[0]->m_q;
		btVector3 c2 = t.m_n[2]->m_q - t.m_n[0]->m_q;
		btVector3 c3 = t.m_n[3]->m_q - t.m_n[0]->m_q;
		btMatrix3x3 Ds(c1.getX(), c2.getX(), c3.getX(),
					   c1.getY(), c2.getY(), c3.getY(),
					   c1.getZ(), c2.getZ(), c3.getZ());
		t.m_F = Ds * t.m_Dm_inverse;

		btSoftBody::TetraScratch& s = m_tetraScratches[i];
		s.m_F = t.m_F;
		s.m_J = t.m_F.determinant();
		btMatrix3x3 C = t.m_F.transpose() * t.m_F;
		s.m_trace = C[0].getX() + C[1].getY() + C[2].getZ();
		s.m_cofF = t.m_F.adjoint().transpose();

		btVector3 a = t.m_n[0]->m_q;
		btVector3 b = t.m_n[1]->m_q;
		btVector3 c = t.m_n[2]->m_q;
		btVector3 d = t.m_n[3]->m_q;
		btVector4 q1(a[0], b[0], c[0], d[0]);
		btVector4 q2(a[1], b[1], c[1], d[1]);
		btVector4 q3(a[2], b[2], c[2], d[2]);
		btMatrix3x3 B(Dot4(q1, t.m_P_inv[0]), Dot4(q1, t.m_P_inv[1]), Dot4(q1, t.m_P_inv[2]),
					  Dot4(q2, t.m_P_inv[0]), Dot4(q2, t.m_P_inv[1]), Dot4(q2, t.m_P_inv[2]),
					  Dot4(q3, t.m_P_inv[0]), Dot4(q3, t.m_P_inv[1]), Dot4(q3, t.m_P_inv[2]));
		q.setRotation(btVector3(0, 0, 1), 0);
		B.getRotation(q);
		btMatrix3x3 Q(q);
		s.m_corotation = Q;
	}
}

void btSoftBody::advanceDeformation()
{
	updateDeformation();
	for (int i = 0; i < m_tetras.size(); ++i)
	{
		m_tetraScratchesTn[i] = m_tetraScratches[i];
	}
}
//
void btSoftBody::Joint::Prepare(btScalar dt, int)
{
	m_bodies[0].activate();
	m_bodies[1].activate();
}

//
void btSoftBody::LJoint::Prepare(btScalar dt, int iterations)
{
	static const btScalar maxdrift = 4;
	Joint::Prepare(dt, iterations);
	m_rpos[0] = m_bodies[0].xform() * m_refs[0];
	m_rpos[1] = m_bodies[1].xform() * m_refs[1];
	m_drift = Clamp(m_rpos[0] - m_rpos[1], maxdrift) * m_erp / dt;
	m_rpos[0] -= m_bodies[0].xform().getOrigin();
	m_rpos[1] -= m_bodies[1].xform().getOrigin();
	m_massmatrix = ImpulseMatrix(m_bodies[0].invMass(), m_bodies[0].invWorldInertia(), m_rpos[0],
								 m_bodies[1].invMass(), m_bodies[1].invWorldInertia(), m_rpos[1]);
	if (m_split > 0)
	{
		m_sdrift = m_massmatrix * (m_drift * m_split);
		m_drift *= 1 - m_split;
	}
	m_drift /= (btScalar)iterations;
}

//
void btSoftBody::LJoint::Solve(btScalar dt, btScalar sor)
{
	const btVector3 va = m_bodies[0].velocity(m_rpos[0]);
	const btVector3 vb = m_bodies[1].velocity(m_rpos[1]);
	const btVector3 vr = va - vb;
	btSoftBody::Impulse impulse;
	impulse.m_asVelocity = 1;
	impulse.m_velocity = m_massmatrix * (m_drift + vr * m_cfm) * sor;
	m_bodies[0].applyImpulse(-impulse, m_rpos[0]);
	m_bodies[1].applyImpulse(impulse, m_rpos[1]);
}

//
void btSoftBody::LJoint::Terminate(btScalar dt)
{
	if (m_split > 0)
	{
		m_bodies[0].applyDImpulse(-m_sdrift, m_rpos[0]);
		m_bodies[1].applyDImpulse(m_sdrift, m_rpos[1]);
	}
}

//
void btSoftBody::AJoint::Prepare(btScalar dt, int iterations)
{
	static const btScalar maxdrift = SIMD_PI / 16;
	m_icontrol->Prepare(this);
	Joint::Prepare(dt, iterations);
	m_axis[0] = m_bodies[0].xform().getBasis() * m_refs[0];
	m_axis[1] = m_bodies[1].xform().getBasis() * m_refs[1];
	m_drift = NormalizeAny(btCross(m_axis[1], m_axis[0]));
	m_drift *= btMin(maxdrift, btAcos(Clamp<btScalar>(btDot(m_axis[0], m_axis[1]), -1, +1)));
	m_drift *= m_erp / dt;
	m_massmatrix = AngularImpulseMatrix(m_bodies[0].invWorldInertia(), m_bodies[1].invWorldInertia());
	if (m_split > 0)
	{
		m_sdrift = m_massmatrix * (m_drift * m_split);
		m_drift *= 1 - m_split;
	}
	m_drift /= (btScalar)iterations;
}

//
void btSoftBody::AJoint::Solve(btScalar dt, btScalar sor)
{
	const btVector3 va = m_bodies[0].angularVelocity();
	const btVector3 vb = m_bodies[1].angularVelocity();
	const btVector3 vr = va - vb;
	const btScalar sp = btDot(vr, m_axis[0]);
	const btVector3 vc = vr - m_axis[0] * m_icontrol->Speed(this, sp);
	btSoftBody::Impulse impulse;
	impulse.m_asVelocity = 1;
	impulse.m_velocity = m_massmatrix * (m_drift + vc * m_cfm) * sor;
	m_bodies[0].applyAImpulse(-impulse);
	m_bodies[1].applyAImpulse(impulse);
}

//
void btSoftBody::AJoint::Terminate(btScalar dt)
{
	if (m_split > 0)
	{
		m_bodies[0].applyDAImpulse(-m_sdrift);
		m_bodies[1].applyDAImpulse(m_sdrift);
	}
}

//
void btSoftBody::CJoint::Prepare(btScalar dt, int iterations)
{
	Joint::Prepare(dt, iterations);
	const bool dodrift = (m_life == 0);
	m_delete = (++m_life) > m_maxlife;
	if (dodrift)
	{
		m_drift = m_drift * m_erp / dt;
		if (m_split > 0)
		{
			m_sdrift = m_massmatrix * (m_drift * m_split);
			m_drift *= 1 - m_split;
		}
		m_drift /= (btScalar)iterations;
	}
	else
	{
		m_drift = m_sdrift = btVector3(0, 0, 0);
	}
}

//
void btSoftBody::CJoint::Solve(btScalar dt, btScalar sor)
{
	const btVector3 va = m_bodies[0].velocity(m_rpos[0]);
	const btVector3 vb = m_bodies[1].velocity(m_rpos[1]);
	const btVector3 vrel = va - vb;
	const btScalar rvac = btDot(vrel, m_normal);
	btSoftBody::Impulse impulse;
	impulse.m_asVelocity = 1;
	impulse.m_velocity = m_drift;
	if (rvac < 0)
	{
		const btVector3 iv = m_normal * rvac;
		const btVector3 fv = vrel - iv;
		impulse.m_velocity += iv + fv * m_friction;
	}
	impulse.m_velocity = m_massmatrix * impulse.m_velocity * sor;

	if (m_bodies[0].m_soft == m_bodies[1].m_soft)
	{
		if ((impulse.m_velocity.getX() == impulse.m_velocity.getX()) && (impulse.m_velocity.getY() == impulse.m_velocity.getY()) &&
			(impulse.m_velocity.getZ() == impulse.m_velocity.getZ()))
		{
			if (impulse.m_asVelocity)
			{
				if (impulse.m_velocity.length() < m_bodies[0].m_soft->m_maxSelfCollisionImpulse)
				{
				}
				else
				{
					m_bodies[0].applyImpulse(-impulse * m_bodies[0].m_soft->m_selfCollisionImpulseFactor, m_rpos[0]);
					m_bodies[1].applyImpulse(impulse * m_bodies[0].m_soft->m_selfCollisionImpulseFactor, m_rpos[1]);
				}
			}
		}
	}
	else
	{
		m_bodies[0].applyImpulse(-impulse, m_rpos[0]);
		m_bodies[1].applyImpulse(impulse, m_rpos[1]);
	}
}

//
void btSoftBody::CJoint::Terminate(btScalar dt)
{
	if (m_split > 0)
	{
		m_bodies[0].applyDImpulse(-m_sdrift, m_rpos[0]);
		m_bodies[1].applyDImpulse(m_sdrift, m_rpos[1]);
	}
}

//
void btSoftBody::applyForces()
{
	BT_PROFILE("SoftBody applyForces");
	//	const btScalar					dt =			m_sst.sdt;
	const btScalar kLF = m_cfg.kLF;
	const btScalar kDG = m_cfg.kDG;
	const btScalar kPR = m_cfg.kPR;
	const btScalar kVC = m_cfg.kVC;
	const bool as_lift = kLF > 0;
	const bool as_drag = kDG > 0;
	const bool as_pressure = kPR != 0;
	const bool as_volume = kVC > 0;
	const bool as_aero = as_lift ||
						 as_drag;
	//const bool						as_vaero =		as_aero	&&
	//												(m_cfg.aeromodel < btSoftBody::eAeroModel::F_TwoSided);
	//const bool						as_faero =		as_aero	&&
	//												(m_cfg.aeromodel >= btSoftBody::eAeroModel::F_TwoSided);
	const bool use_medium = as_aero;
	const bool use_volume = as_pressure ||
							as_volume;
	btScalar volume = 0;
	btScalar ivolumetp = 0;
	btScalar dvolumetv = 0;
	btSoftBody::sMedium medium;
	if (use_volume)
	{
		volume = getVolume();
		ivolumetp = 1 / btFabs(volume) * kPR;
		dvolumetv = (m_pose.m_volume - volume) * kVC;
	}
	/* Per vertex forces			*/
	int i, ni;

	for (i = 0, ni = m_nodes.size(); i < ni; ++i)
	{
		btSoftBody::Node& n = m_nodes[i];
		if (n.m_im > 0)
		{
			if (use_medium)
			{
				/* Aerodynamics			*/
				addAeroForceToNode(m_windVelocity, i);
			}
			/* Pressure				*/
			if (as_pressure)
			{
				n.m_f += n.m_n * (n.m_area * ivolumetp);
			}
			/* Volume				*/
			if (as_volume)
			{
				n.m_f += n.m_n * (n.m_area * dvolumetv);
			}
		}
	}

	/* Per face forces				*/
	for (i = 0, ni = m_faces.size(); i < ni; ++i)
	{
		//	btSoftBody::Face&	f=m_faces[i];

		/* Aerodynamics			*/
		addAeroForceToFace(m_windVelocity, i);
	}
}

//
void btSoftBody::setMaxStress(btScalar maxStress)
{
	m_cfg.m_maxStress = maxStress;
}

//
void btSoftBody::interpolateRenderMesh()
{
	if (m_z.size() > 0)
	{
		for (int i = 0; i < m_renderNodes.size(); ++i)
		{
			const Node* p0 = m_renderNodesParents[i][0];
			const Node* p1 = m_renderNodesParents[i][1];
			const Node* p2 = m_renderNodesParents[i][2];
			btVector3 normal = btCross(p1->m_x - p0->m_x, p2->m_x - p0->m_x);
			btVector3 unit_normal = normal.normalized();
			RenderNode& n = m_renderNodes[i];
			n.m_x.setZero();
			for (int j = 0; j < 3; ++j)
			{
				n.m_x += m_renderNodesParents[i][j]->m_x * m_renderNodesInterpolationWeights[i][j];
			}
			n.m_x += m_z[i] * unit_normal;
		}
	}
	else
	{
		for (int i = 0; i < m_renderNodes.size(); ++i)
		{
			RenderNode& n = m_renderNodes[i];
			n.m_x.setZero();
			for (int j = 0; j < 4; ++j)
			{
				if (m_renderNodesParents[i].size())
				{
					n.m_x += m_renderNodesParents[i][j]->m_x * m_renderNodesInterpolationWeights[i][j];
				}
			}
		}
	}
}

void btSoftBody::setCollisionQuadrature(int N)
{
	for (int i = 0; i <= N; ++i)
	{
		for (int j = 0; i + j <= N; ++j)
		{
			m_quads.push_back(btVector3(btScalar(i) / btScalar(N), btScalar(j) / btScalar(N), btScalar(N - i - j) / btScalar(N)));
		}
	}
}

//
void btSoftBody::PSolve_Anchors(btSoftBody* psb, btScalar kst, btScalar ti)
{
	BT_PROFILE("PSolve_Anchors");
	const btScalar kAHR = psb->m_cfg.kAHR * kst;
	const btScalar dt = psb->m_sst.sdt;
	for (int i = 0, ni = psb->m_anchors.size(); i < ni; ++i)
	{
		const Anchor& a = psb->m_anchors[i];
		const btTransform& t = a.m_body->getWorldTransform();
		Node& n = *a.m_node;
		const btVector3 wa = t * a.m_local;
		const btVector3 va = a.m_body->getVelocityInLocalPoint(a.m_c1) * dt;
		const btVector3 vb = n.m_x - n.m_q;
		const btVector3 vr = (va - vb) + (wa - n.m_x) * kAHR;
		const btVector3 impulse = a.m_c0 * vr * a.m_influence;
		n.m_x += impulse * a.m_c2;
		a.m_body->applyImpulse(-impulse, a.m_c1);
	}
}

//
void btSoftBody::PSolve_RContacts(btSoftBody* psb, btScalar kst, btScalar ti)
{
	BT_PROFILE("PSolve_RContacts");
	const btScalar dt = psb->m_sst.sdt;
	const btScalar mrg = psb->getCollisionShape()->getMargin();
	btMultiBodyJacobianData jacobianData;
	for (int i = 0, ni = psb->m_rcontacts.size(); i < ni; ++i)
	{
		const RContact& c = psb->m_rcontacts[i];
		const sCti& cti = c.m_cti;
		if (cti.m_colObj->hasContactResponse())
		{
			btVector3 va(0, 0, 0);
			btRigidBody* rigidCol = 0;
			btMultiBodyLinkCollider* multibodyLinkCol = 0;
			btScalar* deltaV = NULL;

			if (cti.m_colObj->getInternalType() == btCollisionObject::CO_RIGID_BODY)
			{
				rigidCol = (btRigidBody*)btRigidBody::upcast(cti.m_colObj);
				va = rigidCol ? rigidCol->getVelocityInLocalPoint(c.m_c1) * dt : btVector3(0, 0, 0);
			}
			else if (cti.m_colObj->getInternalType() == btCollisionObject::CO_FEATHERSTONE_LINK)
			{
				multibodyLinkCol = (btMultiBodyLinkCollider*)btMultiBodyLinkCollider::upcast(cti.m_colObj);
				if (multibodyLinkCol)
				{
					const int ndof = multibodyLinkCol->m_multiBody->getNumDofs() + 6;
					jacobianData.m_jacobians.resize(ndof);
					jacobianData.m_deltaVelocitiesUnitImpulse.resize(ndof);
					btScalar* jac = &jacobianData.m_jacobians[0];

					multibodyLinkCol->m_multiBody->fillContactJacobianMultiDof(multibodyLinkCol->m_link, c.m_node->m_x, cti.m_normal, jac, jacobianData.scratch_r, jacobianData.scratch_v, jacobianData.scratch_m);
					deltaV = &jacobianData.m_deltaVelocitiesUnitImpulse[0];
					multibodyLinkCol->m_multiBody->calcAccelerationDeltasMultiDof(&jacobianData.m_jacobians[0], deltaV, jacobianData.scratch_r, jacobianData.scratch_v);

					btScalar vel = 0.0;
					for (int j = 0; j < ndof; ++j)
					{
						vel += multibodyLinkCol->m_multiBody->getVelocityVector()[j] * jac[j];
					}
					va = cti.m_normal * vel * dt;
				}
			}

			const btVector3 vb = c.m_node->m_x - c.m_node->m_q;
			const btVector3 vr = vb - va;
			const btScalar dn = btDot(vr, cti.m_normal);
			if (dn <= SIMD_EPSILON)
			{
				const btScalar dp = btMin((btDot(c.m_node->m_x, cti.m_normal) + cti.m_offset), mrg);
				const btVector3 fv = vr - (cti.m_normal * dn);
				// c0 is the impulse matrix, c3 is 1 - the friction coefficient or 0, c4 is the contact hardness coefficient
				const btVector3 impulse = c.m_c0 * ((vr - (fv * c.m_c3) + (cti.m_normal * (dp * c.m_c4))) * kst);
				c.m_node->m_x -= impulse * c.m_c2;

				if (cti.m_colObj->getInternalType() == btCollisionObject::CO_RIGID_BODY)
				{
					if (rigidCol)
						rigidCol->applyImpulse(impulse, c.m_c1);
				}
				else if (cti.m_colObj->getInternalType() == btCollisionObject::CO_FEATHERSTONE_LINK)
				{
					if (multibodyLinkCol)
					{
						double multiplier = 0.5;
						multibodyLinkCol->m_multiBody->applyDeltaVeeMultiDof(deltaV, -impulse.length() * multiplier);
					}
				}
			}
		}
	}
}

//
void btSoftBody::PSolve_SContacts(btSoftBody* psb, btScalar, btScalar ti)
{
	BT_PROFILE("PSolve_SContacts");

	for (int i = 0, ni = psb->m_scontacts.size(); i < ni; ++i)
	{
		const SContact& c = psb->m_scontacts[i];
		const btVector3& nr = c.m_normal;
		Node& n = *c.m_node;
		Face& f = *c.m_face;
		const btVector3 p = BaryEval(f.m_n[0]->m_x,
									 f.m_n[1]->m_x,
									 f.m_n[2]->m_x,
									 c.m_weights);
		const btVector3 q = BaryEval(f.m_n[0]->m_q,
									 f.m_n[1]->m_q,
									 f.m_n[2]->m_q,
									 c.m_weights);
		const btVector3 vr = (n.m_x - n.m_q) - (p - q);
		btVector3 corr(0, 0, 0);
		btScalar dot = btDot(vr, nr);
		if (dot < 0)
		{
			const btScalar j = c.m_margin - (btDot(nr, n.m_x) - btDot(nr, p));
			corr += c.m_normal * j;
		}
		corr -= ProjectOnPlane(vr, nr) * c.m_friction;
		n.m_x += corr * c.m_cfm[0];
		f.m_n[0]->m_x -= corr * (c.m_cfm[1] * c.m_weights.x());
		f.m_n[1]->m_x -= corr * (c.m_cfm[1] * c.m_weights.y());
		f.m_n[2]->m_x -= corr * (c.m_cfm[1] * c.m_weights.z());
	}
}

//
void btSoftBody::PSolve_Links(btSoftBody* psb, btScalar kst, btScalar ti)
{
	BT_PROFILE("PSolve_Links");
	for (int i = 0, ni = psb->m_links.size(); i < ni; ++i)
	{
		Link& l = psb->m_links[i];
		if (l.m_c0 > 0)
		{
			Node& a = *l.m_n[0];
			Node& b = *l.m_n[1];
			const btVector3 del = b.m_x - a.m_x;
			const btScalar len = del.length2();
			if (l.m_c1 + len > SIMD_EPSILON)
			{
				const btScalar k = ((l.m_c1 - len) / (l.m_c0 * (l.m_c1 + len))) * kst;
				a.m_x -= del * (k * a.m_im);
				b.m_x += del * (k * b.m_im);
			}
		}
	}
}

//
void btSoftBody::VSolve_Links(btSoftBody* psb, btScalar kst)
{
	BT_PROFILE("VSolve_Links");
	for (int i = 0, ni = psb->m_links.size(); i < ni; ++i)
	{
		Link& l = psb->m_links[i];
		Node** n = l.m_n;
		const btScalar j = -btDot(l.m_c3, n[0]->m_v - n[1]->m_v) * l.m_c2 * kst;
		n[0]->m_v += l.m_c3 * (j * n[0]->m_im);
		n[1]->m_v -= l.m_c3 * (j * n[1]->m_im);
	}
}

//
btSoftBody::psolver_t btSoftBody::getSolver(ePSolver::_ solver)
{
	switch (solver)
	{
		case ePSolver::Anchors:
			return (&btSoftBody::PSolve_Anchors);
		case ePSolver::Linear:
			return (&btSoftBody::PSolve_Links);
		case ePSolver::RContacts:
			return (&btSoftBody::PSolve_RContacts);
		case ePSolver::SContacts:
			return (&btSoftBody::PSolve_SContacts);
		default:
		{
		}
	}
	return (0);
}

//
btSoftBody::vsolver_t btSoftBody::getSolver(eVSolver::_ solver)
{
	switch (solver)
	{
		case eVSolver::Linear:
			return (&btSoftBody::VSolve_Links);
		default:
		{
		}
	}
	return (0);
}

void btSoftBody::setSelfCollision(bool useSelfCollision)
{
	m_useSelfCollision = useSelfCollision;
}

bool btSoftBody::useSelfCollision()
{
	return m_useSelfCollision;
}

//
void btSoftBody::defaultCollisionHandler(const btCollisionObjectWrapper* pcoWrap)
{
	switch (m_cfg.collisions & fCollision::RVSmask)
	{
		case fCollision::SDF_RS:
		{
			btSoftColliders::CollideSDF_RS docollide;
			btRigidBody* prb1 = (btRigidBody*)btRigidBody::upcast(pcoWrap->getCollisionObject());
			btTransform wtr = pcoWrap->getWorldTransform();

			const btTransform ctr = pcoWrap->getWorldTransform();
			const btScalar timemargin = (wtr.getOrigin() - ctr.getOrigin()).length();
			const btScalar basemargin = getCollisionShape()->getMargin();
			btVector3 mins;
			btVector3 maxs;
			ATTRIBUTE_ALIGNED16(btDbvtVolume)
			volume;
			pcoWrap->getCollisionShape()->getAabb(pcoWrap->getWorldTransform(),
												  mins,
												  maxs);
			volume = btDbvtVolume::FromMM(mins, maxs);
			volume.Expand(btVector3(basemargin, basemargin, basemargin));
			docollide.psb = this;
			docollide.m_colObj1Wrap = pcoWrap;
			docollide.m_rigidBody = prb1;

			docollide.dynmargin = basemargin + timemargin;
			docollide.stamargin = basemargin;
			m_ndbvt.collideTV(m_ndbvt.m_root, volume, docollide);
		}
		break;
		case fCollision::CL_RS:
		{
			btSoftColliders::CollideCL_RS collider;
			collider.ProcessColObj(this, pcoWrap);
		}
		break;
		case fCollision::SDF_RD:
		{
			btRigidBody* prb1 = (btRigidBody*)btRigidBody::upcast(pcoWrap->getCollisionObject());
			if (this->isActive() && !this->isStaticObject())
			{
				const btTransform wtr = pcoWrap->getWorldTransform();
				const btScalar timemargin = 0;
				const btScalar basemargin = getCollisionShape()->getMargin();
				btVector3 mins;
				btVector3 maxs;
				ATTRIBUTE_ALIGNED16(btDbvtVolume)
				volume;
				pcoWrap->getCollisionShape()->getAabb(wtr,
													  mins,
													  maxs);
				volume = btDbvtVolume::FromMM(mins, maxs);
				volume.Expand(btVector3(basemargin, basemargin, basemargin));
				if (m_cfg.collisions & fCollision::SDF_RDN)
				{
					btSoftColliders::CollideSDF_RD docollideNode;
					docollideNode.psb = this;
					docollideNode.m_colObj1Wrap = pcoWrap;
					docollideNode.m_rigidBody = prb1;
					docollideNode.dynmargin = basemargin + timemargin;
					docollideNode.stamargin = basemargin;
					m_ndbvt.collideTV(m_ndbvt.m_root, volume, docollideNode);
				}

				if (((pcoWrap->getCollisionObject()->getInternalType() == CO_RIGID_BODY) && (m_cfg.collisions & fCollision::SDF_RDF)) || ((pcoWrap->getCollisionObject()->getInternalType() == CO_FEATHERSTONE_LINK) && (m_cfg.collisions & fCollision::SDF_MDF)))
				{
					btSoftColliders::CollideSDF_RDF docollideFace;
					docollideFace.psb = this;
					docollideFace.m_colObj1Wrap = pcoWrap;
					docollideFace.m_rigidBody = prb1;
					docollideFace.dynmargin = basemargin + timemargin;
					docollideFace.stamargin = basemargin;
					m_fdbvt.collideTV(m_fdbvt.m_root, volume, docollideFace);
				}
			}
		}
		break;
	}
}

//
void btSoftBody::defaultCollisionHandler(btSoftBody* psb)
{
	BT_PROFILE("Deformable Collision");
	const int cf = m_cfg.collisions & psb->m_cfg.collisions;
	switch (cf & fCollision::SVSmask)
	{
		case fCollision::CL_SS:
		{
			//support self-collision if CL_SELF flag set
			if (this != psb || psb->m_cfg.collisions & fCollision::CL_SELF)
			{
				btSoftColliders::CollideCL_SS docollide;
				docollide.ProcessSoftSoft(this, psb);
			}
		}
		break;
		case fCollision::VF_SS:
		{
			//only self-collision for Cluster, not Vertex-Face yet
			if (this != psb)
			{
				btSoftColliders::CollideVF_SS docollide;
				/* common					*/
				docollide.mrg = getCollisionShape()->getMargin() +
								psb->getCollisionShape()->getMargin();
				/* psb0 nodes vs psb1 faces	*/
				docollide.psb[0] = this;
				docollide.psb[1] = psb;
				docollide.psb[0]->m_ndbvt.collideTT(docollide.psb[0]->m_ndbvt.m_root,
													docollide.psb[1]->m_fdbvt.m_root,
													docollide);
				/* psb1 nodes vs psb0 faces	*/
				docollide.psb[0] = psb;
				docollide.psb[1] = this;
				docollide.psb[0]->m_ndbvt.collideTT(docollide.psb[0]->m_ndbvt.m_root,
													docollide.psb[1]->m_fdbvt.m_root,
													docollide);
			}
		}
		break;
		case fCollision::VF_DD:
		{
			if (!psb->m_softSoftCollision)
				return;
			if ((psb->isActive() && !psb->isStaticObject()) || (this->isActive() && !this->isStaticObject()))
			{
				if (this != psb)
				{
					btSoftColliders::CollideVF_DD docollide;
					/* common                    */
					docollide.mrg = getCollisionShape()->getMargin() +
									psb->getCollisionShape()->getMargin();
					/* psb0 nodes vs psb1 faces    */
					if (psb->m_tetras.size() > 0)
						docollide.useFaceNormal = true;
					else
						docollide.useFaceNormal = false;
					docollide.psb[0] = this;
					docollide.psb[1] = psb;
					docollide.psb[0]->m_ndbvt.collideTT(docollide.psb[0]->m_ndbvt.m_root,
														docollide.psb[1]->m_fdbvt.m_root,
														docollide);

					/* psb1 nodes vs psb0 faces    */
					if (this->m_tetras.size() > 0)
						docollide.useFaceNormal = true;
					else
						docollide.useFaceNormal = false;
					docollide.psb[0] = psb;
					docollide.psb[1] = this;
					docollide.psb[0]->m_ndbvt.collideTT(docollide.psb[0]->m_ndbvt.m_root,
														docollide.psb[1]->m_fdbvt.m_root,
														docollide);
				}
				else
				{
					if (psb->useSelfCollision())
					{
						btSoftColliders::CollideFF_DD docollide;
						docollide.mrg = 2 * getCollisionShape()->getMargin();
						docollide.psb[0] = this;
						docollide.psb[1] = psb;
						if (this->m_tetras.size() > 0)
							docollide.useFaceNormal = true;
						else
							docollide.useFaceNormal = false;
						/* psb0 faces vs psb0 faces    */
						calculateNormalCone(this->m_fdbvnt);
						this->m_fdbvt.selfCollideT(m_fdbvnt, docollide);
					}
				}
			}
		}
		break;
		default:
		{
		}
	}
}

int btSoftBody::findClosestNodeByMapping(int part, int triIndex, const btVector3& p) const
{
	auto tetraIndices = getCollisionShape()->getMappingForTri(part, triIndex);
	std::map<btScalar, int> nodeDistances;
	for (auto tetraIndex : tetraIndices)
	{
		const auto& tetra = m_tetras[tetraIndex];
		for (auto n = 0; n < 4; ++n)
		{
			auto& node = tetra.m_n[n];
			nodeDistances.insert({(p - node->m_x).length2(), node->local_index});
		}
	}
	if (!nodeDistances.empty())
		return nodeDistances.begin()->second;
	else
		return -1;
}

std::vector<int> btSoftBody::findNClosestNodesLinearComplexity(const btVector3& p, int N) const
{
	std::vector<std::tuple<int, btScalar>> nodeDistances;

	for (int i = 0; i < m_nodes.size(); ++i)
	{
		auto& n = m_nodes[i];
		if (n.m_im <= 0.0)
			continue;
		btScalar distSq = (p - n.m_x).length2();
		nodeDistances.emplace_back(i, distSq);
	}

	std::nth_element(nodeDistances.begin(), nodeDistances.begin() + N, nodeDistances.end(),
					 [](const std::tuple<int, btScalar>& a, const std::tuple<int, btScalar>& b)
					 {
						 return std::get<1>(a) < std::get<1>(b);
					 });

	// Prepare the result vector for N closest faces
	std::vector<int> closestNodes;
	for (int i = 0; i < std::min(N, static_cast<int>(nodeDistances.size())); ++i)
	{
		closestNodes.emplace_back(std::get<0>(nodeDistances[i]));
	}

	return closestNodes;
}

// 60 degrees to make sure we cover the whole sphere with kMaxBucketsPerNode = 6 buckets
static const btScalar kMergeCos = btCos(btRadians(60.f));
// One for each cardinal direction (-x, x, -y, y, -z, z)
static const int kMaxBucketsPerNode = 6;

template <typename T, typename U>
bool mergeContactIntoBucket(const btCollisionObject* body, T& contacts, const btVector3& contactNormalOnSoftCollisionMesh, const btScalar& penetrationDepth, btSoftBody::Node& n, btSoftBody::Node& nOther)
{
	auto sameDirection = [](const btVector3& a, const btVector3& b)
	{
		// a·b ≥ |a||b|cosθ  (both already unit-length here)
		return a.dot(b) > kMergeCos;
	};

	bool merged = false;
	int bucketCount = 0;
	for (int idx = 0, cnt = contacts.size(); idx < cnt; ++idx)
	{
		auto& c = contacts[idx];
		U* cCti;
		if constexpr (std::is_same_v<U, btSoftBody::sCti>)
			cCti = &contacts[idx].m_cti;
		else
			cCti = &contacts[idx];

		bool irrelevantContact;
		if constexpr (std::is_same_v<U, btSoftBody::sCti>)
			irrelevantContact = cCti->m_colObj != body || c.m_node != &n;
		else
			irrelevantContact = c.m_colObj != body || c.m_node0 != &n || c.m_node1 != &nOther;

		if (irrelevantContact)
			continue;

		++bucketCount;

		if (sameDirection(cCti->m_normal, contactNormalOnSoftCollisionMesh))
		{
			// Merge into this bucket
			const int newCount = cCti->m_count + 1;
			const btVector3 avgN = (cCti->m_normal * cCti->m_count +
									contactNormalOnSoftCollisionMesh) /
								   newCount;

			cCti->m_normal = avgN.normalized();
			cCti->m_offset = btMin(cCti->m_offset, penetrationDepth);
			cCti->m_count = newCount;

			merged = true;
			break;
		}
	}

	// Whole sphere should be covered by these conical buckets, so if we now have number of bins which is equal or larger than kMaxBucketsPerNode
	// it indicates a bug somewhere
	btAssert(bucketCount <= kMaxBucketsPerNode &&
			 "Directional coverage theory broken – investigate!");

	if (bucketCount > kMaxBucketsPerNode)
		return true;  // safety valve – should never fire in practice

	return merged;
}

void btSoftBody::skinSoftRigidCollisionHandler(const btCollisionObjectWrapper* rigidWrap, int part0, int index0, const btVector3& contactPointOnSoftCollisionMesh, btVector3 contactNormalOnSoftCollisionMesh,
											   btScalar penetrationDepth, const bool penetrating, btScalar* contactPointImpulseMagnitude)
{
	contactNormalOnSoftCollisionMesh = -contactNormalOnSoftCollisionMesh;
	const auto rigidBody = static_cast<const btRigidBody*>(rigidWrap->getCollisionObject());

	// Uniform distribution for [0.0, 1.0]

	/*float r = dist(gen);
	float g = dist(gen);
	float b = dist(gen);
	fprintf(stderr, "drawpoint \"pt\" [%f,%f,%f][%f,%f,%f,1]\n", contactPointOnSoftCollisionMesh.x(), contactPointOnSoftCollisionMesh.y(), contactPointOnSoftCollisionMesh.z(), r, g, b);
	auto lineStart = contactPointOnSoftCollisionMesh;
	auto lineEnd = contactPointOnSoftCollisionMesh + contactNormalOnSoftCollisionMesh * 10.0;
	fprintf(stderr, "drawline \"line\" [%f,%f,%f][%f,%f,%f][%f,%f,%f,1] \n", lineStart.x(), lineStart.y(), lineStart.z(),
			lineEnd.x(), lineEnd.y(), lineEnd.z(), r, g, b);*/

	auto nodeIndex = findClosestNodeByMapping(part0, index0, contactPointOnSoftCollisionMesh);

	if (nodeIndex == -1)
		return;

	btSoftBody::Node& n = m_nodes[nodeIndex];

	bool merged = mergeContactIntoBucket<btAlignedObjectArray<btSoftBody::DeformableNodeRigidContact>, btSoftBody::sCti>(rigidBody, m_nodeRigidContacts, contactNormalOnSoftCollisionMesh, penetrationDepth, n, n);
	if (merged)
		return;

	btSoftBody::DeformableNodeRigidContact c;

	c.m_cti.m_colObj = rigidBody;
	c.m_cti.m_normal = contactNormalOnSoftCollisionMesh;
	c.m_cti.m_offset = penetrationDepth;
	c.m_cti.m_contact_point_impulse_magnitude = contactPointImpulseMagnitude;

	btScalar ima = n.m_im;
	const btScalar imb = rigidBody ? rigidBody->getInvMass() : 0.f;

	const btScalar ms = ima + imb;
	if (ms > 0)
	{
		btSoftBody::sCti& cti = c.m_cti;

		c.m_node = &n;
		// friction is handled by the nodes to prevent sticking
		//                    const btScalar fc = 0;
		const btScalar fc = m_cfg.kDF * rigidWrap->getCollisionObject()->getFriction();

		c.m_c2 = ima;
		c.m_c3 = fc;
		c.m_c4 = rigidWrap->getCollisionObject()->isStaticOrKinematicObject() ? m_cfg.kKHR : m_cfg.kCHR;
		c.m_c5 = Diagonal(ima);
		if (cti.m_colObj->getInternalType() == btCollisionObject::CO_RIGID_BODY)
		{
			const btTransform& wtr = rigidBody ? rigidBody->getWorldTransform() : rigidWrap->getCollisionObject()->getWorldTransform();
			static const btMatrix3x3 iwiStatic(0, 0, 0, 0, 0, 0, 0, 0, 0);
			const btMatrix3x3& iwi = rigidBody ? rigidBody->getInvInertiaTensorWorld() : iwiStatic;
			const btVector3 ra = contactPointOnSoftCollisionMesh - wtr.getOrigin();

			// we do not scale the impulse matrix by dt
			c.m_c0 = ImpulseMatrix(1, ima, imb, iwi, ra);
			c.m_c1 = ra;
		}
		else if (cti.m_colObj->getInternalType() == btCollisionObject::CO_FEATHERSTONE_LINK)
		{
			btMultiBodyLinkCollider* multibodyLinkCol = (btMultiBodyLinkCollider*)btMultiBodyLinkCollider::upcast(cti.m_colObj);
			if (multibodyLinkCol)
			{
				btVector3 normal = cti.m_normal;
				btVector3 t1 = generateUnitOrthogonalVector(normal);
				btVector3 t2 = btCross(normal, t1);
				btMultiBodyJacobianData jacobianData_normal, jacobianData_t1, jacobianData_t2;
				findJacobian(multibodyLinkCol, jacobianData_normal, contactPointOnSoftCollisionMesh, normal);
				findJacobian(multibodyLinkCol, jacobianData_t1, contactPointOnSoftCollisionMesh, t1);
				findJacobian(multibodyLinkCol, jacobianData_t2, contactPointOnSoftCollisionMesh, t2);

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
				btMatrix3x3 local_impulse_matrix = (Diagonal(ima) + OuterProduct(J_n, J_t1, J_t2, u_n, u_t1, u_t2, ndof)).inverse();
				c.m_c0 = rot.transpose() * local_impulse_matrix * rot;
				c.jacobianData_normal = jacobianData_normal;
				c.jacobianData_t1 = jacobianData_t1;
				c.jacobianData_t2 = jacobianData_t2;
				c.t1 = t1;
				c.t2 = t2;
			}
		}
		m_nodeRigidContacts.push_back(c);
	}
}

//int cnt = 0;
//int origsDrawn = false;

void btSoftBody::skinSoftSoftCollisionHandler(btSoftBody* otherSoft, int part0, int index0, int part1, int index1, const btVector3& contactPointOnSoftCollisionMesh, btVector3 contactNormalOnSoftCollisionMesh, btScalar distance, const bool penetrating, btScalar* contactPointImpulseMagnitude)
{
	contactNormalOnSoftCollisionMesh = -contactNormalOnSoftCollisionMesh;

	/*fprintf(stderr, "drawpoint \"pt\" [%f,%f,%f]\n", contactPointOnSoftCollisionMesh.x(), contactPointOnSoftCollisionMesh.y(), contactPointOnSoftCollisionMesh.z());
	auto lineStart = contactPointOnSoftCollisionMesh;
	auto lineEnd = contactPointOnSoftCollisionMesh + contactNormalOnSoftCollisionMesh * 10.0;
	fprintf(stderr, "drawline \"line\" [%f,%f,%f][%f,%f,%f] \n", lineStart.x(), lineStart.y(), lineStart.z(),
			lineEnd.x(), lineEnd.y(), lineEnd.z());*/

	auto nodeIndex = findClosestNodeByMapping(part0, index0, contactPointOnSoftCollisionMesh);
	auto nodeIndexOther = otherSoft->findClosestNodeByMapping(part1, index1, contactPointOnSoftCollisionMesh);

	if (nodeIndex == -1 || nodeIndexOther == -1)
		return;

	btSoftBody::Node& n = m_nodes[nodeIndex];
	btSoftBody::Node& nOther = otherSoft->m_nodes[nodeIndexOther];

	bool merged = mergeContactIntoBucket<btAlignedObjectArray<DeformableNodeNodeContact>, DeformableNodeNodeContact>(otherSoft, m_nodeNodeContacts, contactNormalOnSoftCollisionMesh, distance, n, nOther);
	if (merged)
		return;

	//if (!origsDrawn)
	{
		/*for (auto i = 0; i < m_faces.size(); ++i)
		{
			auto& f = m_faces[i];
			fprintf(stderr, "drawline \"lin%d\" [%f, %f, %f] [%f, %f, %f]\n", i, f.m_n[0]->m_x.x(), f.m_n[0]->m_x.y(), f.m_n[0]->m_x.z(),
					f.m_n[1]->m_x.x(), f.m_n[1]->m_x.y(), f.m_n[1]->m_x.z());
			fprintf(stderr, "drawline \"lin%d\" [%f, %f, %f] [%f, %f, %f]\n", i, f.m_n[0]->m_x.x(), f.m_n[0]->m_x.y(), f.m_n[0]->m_x.z(),
					f.m_n[2]->m_x.x(), f.m_n[2]->m_x.y(), f.m_n[2]->m_x.z());
			fprintf(stderr, "drawline \"lin%d\" [%f, %f, %f] [%f, %f, %f]\n", i, f.m_n[1]->m_x.x(), f.m_n[1]->m_x.y(), f.m_n[1]->m_x.z(),
					f.m_n[2]->m_x.x(), f.m_n[2]->m_x.y(), f.m_n[2]->m_x.z());
		}
		for (auto i = 0; i < otherSoft->m_faces.size(); ++i)
		{
			auto& f = otherSoft->m_faces[i];
			fprintf(stderr, "drawline \"otherlin%d\" [%f, %f, %f] [%f, %f, %f]\n", i, f.m_n[0]->m_x.x(), f.m_n[0]->m_x.y(), f.m_n[0]->m_x.z(),
					f.m_n[1]->m_x.x(), f.m_n[1]->m_x.y(), f.m_n[1]->m_x.z());
			fprintf(stderr, "drawline \"otherlin%d\" [%f, %f, %f] [%f, %f, %f]\n", i, f.m_n[0]->m_x.x(), f.m_n[0]->m_x.y(), f.m_n[0]->m_x.z(),
					f.m_n[2]->m_x.x(), f.m_n[2]->m_x.y(), f.m_n[2]->m_x.z());
			fprintf(stderr, "drawline \"otherlin%d\" [%f, %f, %f] [%f, %f, %f]\n", i, f.m_n[1]->m_x.x(), f.m_n[1]->m_x.y(), f.m_n[1]->m_x.z(),
					f.m_n[2]->m_x.x(), f.m_n[2]->m_x.y(), f.m_n[2]->m_x.z());
		}*/
		//origsDrawn = true;
	}

	/*fprintf(stderr, "createPointHelperObject \"contactpoint\" [%f, %f, %f] \n", contactPointOnSoftCollisionMesh.x(), contactPointOnSoftCollisionMesh.y(), contactPointOnSoftCollisionMesh.z());
		fprintf(stderr, "createPointHelperObject \"n%d\" [%f, %f, %f] \n", getUserIndex(), n.m_x.x(), n.m_x.y(), n.m_x.z());
		fprintf(stderr, "createPointHelperObject \"nOther%d\" [%f, %f, %f] \n", otherSoft->getUserIndex(), nOther.m_x.x(), nOther.m_x.y(), nOther.m_x.z());*/

	{
		btSoftBody::DeformableNodeNodeContact c;
		c.m_normal = contactNormalOnSoftCollisionMesh;
		c.m_offset = distance;
		c.m_node0 = &n;
		c.m_node1 = &nOther;
		c.m_colObj = otherSoft;
		c.m_friction = m_cfg.kDF * otherSoft->m_cfg.kDF;
		c.m_contact_point_impulse_magnitude = contactPointImpulseMagnitude;
		m_nodeNodeContacts.push_back(c);
	}
}

void btSoftBody::applyRepulsionForce(btScalar timeStep, bool applySpringForce)
{
	btAlignedObjectArray<int> indices;
	{
		// randomize the order of repulsive force
		indices.resize(m_faceNodeContacts.size());
		for (int i = 0; i < m_faceNodeContacts.size(); ++i)
			indices[i] = i;
#define NEXTRAND (seed = (1664525L * seed + 1013904223L) & 0xffffffff)
		int i, ni;

		for (i = 0, ni = indices.size(); i < ni; ++i)
		{
			btSwap(indices[i], indices[NEXTRAND % ni]);
		}
	}
	for (int k = 0; k < m_faceNodeContacts.size(); ++k)
	{
		int idx = indices[k];
		btSoftBody::DeformableFaceNodeContact& c = m_faceNodeContacts[idx];
		btSoftBody::Node* node = c.m_node;
		btSoftBody::Face* face = c.m_face;
		const btVector3& w = c.m_bary;
		const btVector3& n = c.m_normal;
		btVector3 l = node->m_x - BaryEval(face->m_n[0]->m_x, face->m_n[1]->m_x, face->m_n[2]->m_x, w);
		btScalar d = c.m_margin - n.dot(l);
		d = btMax(btScalar(0), d);

		const btVector3& va = node->m_v;
		btVector3 vb = BaryEval(face->m_n[0]->m_v, face->m_n[1]->m_v, face->m_n[2]->m_v, w);
		btVector3 vr = va - vb;
		const btScalar vn = btDot(vr, n);  // dn < 0 <==> opposing
		if (vn > OVERLAP_REDUCTION_FACTOR * d / timeStep)
			continue;
		btVector3 vt = vr - vn * n;
		btScalar I = 0;
		btScalar mass = node->m_im == 0 ? 0 : btScalar(1) / node->m_im;
		if (applySpringForce)
			I = -btMin(m_repulsionStiffness * timeStep * d, mass * (OVERLAP_REDUCTION_FACTOR * d / timeStep - vn));
		if (vn < 0)
			I += 0.5 * mass * vn;
		int face_penetration = 0, node_penetration = node->m_constrained;
		for (int i = 0; i < 3; ++i)
			face_penetration |= face->m_n[i]->m_constrained;
		btScalar I_tilde = 2.0 * I / (1.0 + w.length2());

		//             double the impulse if node or face is constrained.
		if (face_penetration > 0 || node_penetration > 0)
		{
			I_tilde *= 2.0;
		}
		if (face_penetration <= 0)
		{
			for (int j = 0; j < 3; ++j)
				face->m_n[j]->m_v += w[j] * n * I_tilde * node->m_im;
		}
		if (node_penetration <= 0)
		{
			node->m_v -= I_tilde * node->m_im * n;
		}

		// apply frictional impulse
		btScalar vt_norm = vt.safeNorm();
		if (vt_norm > SIMD_EPSILON)
		{
			btScalar delta_vn = -2 * I * node->m_im;
			btScalar mu = c.m_friction;
			btScalar vt_new = btMax(btScalar(1) - mu * delta_vn / (vt_norm + SIMD_EPSILON), btScalar(0)) * vt_norm;
			I = 0.5 * mass * (vt_norm - vt_new);
			vt.safeNormalize();
			I_tilde = 2.0 * I / (1.0 + w.length2());
			//                 double the impulse if node or face is constrained.
			if (face_penetration > 0 || node_penetration > 0)
				I_tilde *= 2.0;
			if (face_penetration <= 0)
			{
				for (int j = 0; j < 3; ++j)
					face->m_n[j]->m_v += w[j] * vt * I_tilde * (face->m_n[j])->m_im;
			}
			if (node_penetration <= 0)
			{
				node->m_v -= I_tilde * node->m_im * vt;
			}
		}
	}

	//fprintf(stderr, "m_nodeNodeContacts %d idx %d\n", m_nodeNodeContacts.size(), getUserIndex());
	for (int k = 0; k < m_nodeNodeContacts.size(); ++k)
	{
		//int idx = indices[k];
		btSoftBody::DeformableNodeNodeContact& c = m_nodeNodeContacts[k];
		btSoftBody::Node* node0 = c.m_node0;
		btSoftBody::Node* node1 = c.m_node1;
		const btVector3& n = c.m_normal;
		//fprintf(stderr, "c.m_normal %f %f %f\n", c.m_normal.x(), c.m_normal.y(), c.m_normal.z());

		const btVector3& va = node0->m_v;
		btVector3 vb = node1->m_v;
		btVector3 vr = va - vb;
		const btScalar vn = btDot(vr, n);  // dn < 0 <==> opposing
		btVector3 vt = vr - vn * n;
		btScalar vtNorm = vt.safeNorm();

		//fprintf(stderr, "vn %f\n", vn);
		if (vn >= 0.0)
		{
			//fprintf(stderr, "CONTINUE\n");
			continue;
		}

		btScalar invMass0 = node0->m_im;
		btScalar invMass1 = node1->m_im;
		btScalar invMassSum = invMass0 + invMass1;
		if (invMassSum < SIMD_EPSILON)
			continue;

		// Typical "collision" style impulse
		btScalar penetration = 0.0;
		if (c.m_offset < 0.0)
			penetration = c.m_offset * m_softVsSoftContactStiffness;
		btScalar restitution = 0.0;  // TODO coefficient of restitution for soft bodies by mapping the value of btCollisionObject::m_restitution and using it here. Will have to be also done in btSoftBody::skinSoftRigidCollisionHandler.
		btScalar jCollision = -(1.0 + restitution) * (vn + penetration) / invMassSum;
		btScalar jTotal = jCollision;

		if (c.m_contact_point_impulse_magnitude)
			*c.m_contact_point_impulse_magnitude = std::abs(jTotal);

		{
			auto delta = (jTotal * invMass0) * n;
			if (node0->m_constrained != 0 && !delta.fuzzyZero())
			{
				auto len = delta.length();
				delta = (delta / len) * (len * 10.0);
			}
			//fprintf(stderr, "node0->m_v delta %d %f %f %f\n", node0->local_index, delta.x(), delta.y(), delta.z());
			if (invMass0 != 0.0)
				node0->m_v += delta;
			//fprintf(stderr, "node0->m_v %f %f %f\n", node0->m_v.x(), node0->m_v.y(), node0->m_v.z());
		}

		{
			auto delta = (jTotal * invMass1) * n;
			if (node1->m_constrained != 0 && !delta.fuzzyZero())
			{
				auto len = delta.length();
				delta = (delta / len) * (len * 10.0);
			}
			//fprintf(stderr, "node1->m_v delta %d %f %f %f\n", node1->local_index, delta.x(), delta.y(), delta.z());
			if (invMass1 != 0.0)
				node1->m_v -= delta;
			//fprintf(stderr, "node1->m_v %f %f %f\n", node1->m_v.x(), node1->m_v.y(), node1->m_v.z());
		}

		// apply frictional impulse
		btScalar vt_norm = vt.safeNorm();
		if (vt_norm > SIMD_EPSILON)
		{
			// friction coefficient
			btScalar mu = c.m_friction;

			// ideal friction impulse (no clamp)
			// note: we re-check relative velocity in tangent after normal impulse if you prefer
			btVector3 vtDir = vt / vtNorm;
			// friction tries to cancel tangential velocity
			// jFric = - vtNorm / invMassSum
			btScalar jFric = -vtNorm / invMassSum;

			// Coulomb clamp: limit friction by mu * normal impulse
			btScalar frictionLimit = mu * jTotal;  // if j is positive normal impulse
			btClamp(jFric, -frictionLimit, frictionLimit);

			// Apply friction impulse
			if (!node0->m_constrained)
			{
				auto delta = (jFric * invMass0) * vtDir;
				node0->m_v += delta;
			}
			if (!node1->m_constrained)
			{
				auto delta = (jFric * invMass1) * vtDir;
				node1->m_v -= delta;
			}
		}
	}
}

void btSoftBody::geometricCollisionHandler(btSoftBody* psb)
{
	if ((psb->isActive() && !psb->isStaticObject()) || (this->isActive() && !this->isStaticObject()))
	{
		if (this != psb)
		{
			btSoftColliders::CollideCCD docollide;
			/* common                    */
			docollide.mrg = SAFE_EPSILON;  // for rounding error instead of actual margin
			docollide.dt = psb->m_sst.sdt;
			/* psb0 nodes vs psb1 faces    */
			if (psb->m_tetras.size() > 0)
				docollide.useFaceNormal = true;
			else
				docollide.useFaceNormal = false;
			docollide.psb[0] = this;
			docollide.psb[1] = psb;
			docollide.psb[0]->m_ndbvt.collideTT(docollide.psb[0]->m_ndbvt.m_root,
												docollide.psb[1]->m_fdbvt.m_root,
												docollide);
			/* psb1 nodes vs psb0 faces    */
			if (this->m_tetras.size() > 0)
				docollide.useFaceNormal = true;
			else
				docollide.useFaceNormal = false;
			docollide.psb[0] = psb;
			docollide.psb[1] = this;
			docollide.psb[0]->m_ndbvt.collideTT(docollide.psb[0]->m_ndbvt.m_root,
												docollide.psb[1]->m_fdbvt.m_root,
												docollide);
		}
		else
		{
			if (psb->useSelfCollision())
			{
				btSoftColliders::CollideCCD docollide;
				docollide.mrg = SAFE_EPSILON;
				docollide.psb[0] = this;
				docollide.psb[1] = psb;
				docollide.dt = psb->m_sst.sdt;
				if (this->m_tetras.size() > 0)
					docollide.useFaceNormal = true;
				else
					docollide.useFaceNormal = false;
				/* psb0 faces vs psb0 faces    */
				calculateNormalCone(this->m_fdbvnt);  // should compute this outside of this scope
				this->m_fdbvt.selfCollideT(m_fdbvnt, docollide);
			}
		}
	}
}

void btSoftBody::setWindVelocity(const btVector3& velocity)
{
	m_windVelocity = velocity;
}

const btVector3& btSoftBody::getWindVelocity()
{
	return m_windVelocity;
}

int btSoftBody::calculateSerializeBufferSize() const
{
	int sz = sizeof(btSoftBodyData);
	return sz;
}

///fills the dataBuffer and returns the struct name (and 0 on failure)
const char* btSoftBody::serialize(void* dataBuffer, class btSerializer* serializer) const
{
	btSoftBodyData* sbd = (btSoftBodyData*)dataBuffer;

	btCollisionObject::serialize(&sbd->m_collisionObjectData, serializer);

	btHashMap<btHashPtr, int> m_nodeIndexMap;

	sbd->m_numMaterials = m_materials.size();
	sbd->m_materials = sbd->m_numMaterials ? (SoftBodyMaterialData**)serializer->getUniquePointer((void*)&m_materials) : 0;

	if (sbd->m_materials)
	{
		int sz = sizeof(SoftBodyMaterialData*);
		int numElem = sbd->m_numMaterials;
		btChunk* chunk = serializer->allocate(sz, numElem);
		//SoftBodyMaterialData** memPtr = chunk->m_oldPtr;
		SoftBodyMaterialData** memPtr = (SoftBodyMaterialData**)chunk->m_oldPtr;
		for (int i = 0; i < numElem; i++, memPtr++)
		{
			btSoftBody::Material* mat = m_materials[i];
			*memPtr = mat ? (SoftBodyMaterialData*)serializer->getUniquePointer((void*)mat) : 0;
			if (!serializer->findPointer(mat))
			{
				//serialize it here
				btChunk* chunk = serializer->allocate(sizeof(SoftBodyMaterialData), 1);
				SoftBodyMaterialData* memPtr = (SoftBodyMaterialData*)chunk->m_oldPtr;
				memPtr->m_flags = mat->m_flags;
				memPtr->m_angularStiffness = mat->m_kAST;
				memPtr->m_linearStiffness = mat->m_kLST;
				memPtr->m_volumeStiffness = mat->m_kVST;
				serializer->finalizeChunk(chunk, "SoftBodyMaterialData", BT_SBMATERIAL_CODE, mat);
			}
		}
		serializer->finalizeChunk(chunk, "SoftBodyMaterialData", BT_ARRAY_CODE, (void*)&m_materials);
	}

	sbd->m_numNodes = m_nodes.size();
	sbd->m_nodes = sbd->m_numNodes ? (SoftBodyNodeData*)serializer->getUniquePointer((void*)&m_nodes) : 0;
	if (sbd->m_nodes)
	{
		int sz = sizeof(SoftBodyNodeData);
		int numElem = sbd->m_numNodes;
		btChunk* chunk = serializer->allocate(sz, numElem);
		SoftBodyNodeData* memPtr = (SoftBodyNodeData*)chunk->m_oldPtr;
		for (int i = 0; i < numElem; i++, memPtr++)
		{
			m_nodes[i].m_f.serializeFloat(memPtr->m_accumulatedForce);
			memPtr->m_area = m_nodes[i].m_area;
			memPtr->m_attach = m_nodes[i].m_battach;
			memPtr->m_inverseMass = m_nodes[i].m_im;
			memPtr->m_material = m_nodes[i].m_material ? (SoftBodyMaterialData*)serializer->getUniquePointer((void*)m_nodes[i].m_material) : 0;
			m_nodes[i].m_n.serializeFloat(memPtr->m_normal);
			m_nodes[i].m_x.serializeFloat(memPtr->m_position);
			m_nodes[i].m_q.serializeFloat(memPtr->m_previousPosition);
			m_nodes[i].m_v.serializeFloat(memPtr->m_velocity);
			m_nodeIndexMap.insert(&m_nodes[i], i);
		}
		serializer->finalizeChunk(chunk, "SoftBodyNodeData", BT_SBNODE_CODE, (void*)&m_nodes);
	}

	sbd->m_numLinks = m_links.size();
	sbd->m_links = sbd->m_numLinks ? (SoftBodyLinkData*)serializer->getUniquePointer((void*)&m_links[0]) : 0;
	if (sbd->m_links)
	{
		int sz = sizeof(SoftBodyLinkData);
		int numElem = sbd->m_numLinks;
		btChunk* chunk = serializer->allocate(sz, numElem);
		SoftBodyLinkData* memPtr = (SoftBodyLinkData*)chunk->m_oldPtr;
		for (int i = 0; i < numElem; i++, memPtr++)
		{
			memPtr->m_bbending = m_links[i].m_bbending;
			memPtr->m_material = m_links[i].m_material ? (SoftBodyMaterialData*)serializer->getUniquePointer((void*)m_links[i].m_material) : 0;
			memPtr->m_nodeIndices[0] = m_links[i].m_n[0] ? m_links[i].m_n[0] - &m_nodes[0] : -1;
			memPtr->m_nodeIndices[1] = m_links[i].m_n[1] ? m_links[i].m_n[1] - &m_nodes[0] : -1;
			btAssert(memPtr->m_nodeIndices[0] < m_nodes.size());
			btAssert(memPtr->m_nodeIndices[1] < m_nodes.size());
			memPtr->m_restLength = m_links[i].m_rl;
		}
		serializer->finalizeChunk(chunk, "SoftBodyLinkData", BT_ARRAY_CODE, (void*)&m_links[0]);
	}

	sbd->m_numFaces = m_faces.size();
	sbd->m_faces = sbd->m_numFaces ? (SoftBodyFaceData*)serializer->getUniquePointer((void*)&m_faces[0]) : 0;
	if (sbd->m_faces)
	{
		int sz = sizeof(SoftBodyFaceData);
		int numElem = sbd->m_numFaces;
		btChunk* chunk = serializer->allocate(sz, numElem);
		SoftBodyFaceData* memPtr = (SoftBodyFaceData*)chunk->m_oldPtr;
		for (int i = 0; i < numElem; i++, memPtr++)
		{
			memPtr->m_material = m_faces[i].m_material ? (SoftBodyMaterialData*)serializer->getUniquePointer((void*)m_faces[i].m_material) : 0;
			m_faces[i].m_normal.serializeFloat(memPtr->m_normal);
			for (int j = 0; j < 3; j++)
			{
				memPtr->m_nodeIndices[j] = m_faces[i].m_n[j] ? m_faces[i].m_n[j] - &m_nodes[0] : -1;
			}
			memPtr->m_restArea = m_faces[i].m_ra;
		}
		serializer->finalizeChunk(chunk, "SoftBodyFaceData", BT_ARRAY_CODE, (void*)&m_faces[0]);
	}

	sbd->m_numTetrahedra = m_tetras.size();
	sbd->m_tetrahedra = sbd->m_numTetrahedra ? (SoftBodyTetraData*)serializer->getUniquePointer((void*)&m_tetras[0]) : 0;
	if (sbd->m_tetrahedra)
	{
		int sz = sizeof(SoftBodyTetraData);
		int numElem = sbd->m_numTetrahedra;
		btChunk* chunk = serializer->allocate(sz, numElem);
		SoftBodyTetraData* memPtr = (SoftBodyTetraData*)chunk->m_oldPtr;
		for (int i = 0; i < numElem; i++, memPtr++)
		{
			for (int j = 0; j < 4; j++)
			{
				m_tetras[i].m_c0[j].serializeFloat(memPtr->m_c0[j]);
				memPtr->m_nodeIndices[j] = m_tetras[i].m_n[j] ? m_tetras[i].m_n[j] - &m_nodes[0] : -1;
			}
			memPtr->m_c1 = m_tetras[i].m_c1;
			memPtr->m_c2 = m_tetras[i].m_c2;
			memPtr->m_material = m_tetras[i].m_material ? (SoftBodyMaterialData*)serializer->getUniquePointer((void*)m_tetras[i].m_material) : 0;
			memPtr->m_restVolume = m_tetras[i].m_rv;
		}
		serializer->finalizeChunk(chunk, "SoftBodyTetraData", BT_ARRAY_CODE, (void*)&m_tetras[0]);
	}

	sbd->m_numAnchors = m_anchors.size();
	sbd->m_anchors = sbd->m_numAnchors ? (SoftRigidAnchorData*)serializer->getUniquePointer((void*)&m_anchors[0]) : 0;
	if (sbd->m_anchors)
	{
		int sz = sizeof(SoftRigidAnchorData);
		int numElem = sbd->m_numAnchors;
		btChunk* chunk = serializer->allocate(sz, numElem);
		SoftRigidAnchorData* memPtr = (SoftRigidAnchorData*)chunk->m_oldPtr;
		for (int i = 0; i < numElem; i++, memPtr++)
		{
			m_anchors[i].m_c0.serializeFloat(memPtr->m_c0);
			m_anchors[i].m_c1.serializeFloat(memPtr->m_c1);
			memPtr->m_c2 = m_anchors[i].m_c2;
			m_anchors[i].m_local.serializeFloat(memPtr->m_localFrame);
			memPtr->m_nodeIndex = m_anchors[i].m_node ? m_anchors[i].m_node - &m_nodes[0] : -1;

			memPtr->m_rigidBody = m_anchors[i].m_body ? (btRigidBodyData*)serializer->getUniquePointer((void*)m_anchors[i].m_body) : 0;
			btAssert(memPtr->m_nodeIndex < m_nodes.size());
		}
		serializer->finalizeChunk(chunk, "SoftRigidAnchorData", BT_ARRAY_CODE, (void*)&m_anchors[0]);
	}

	sbd->m_config.m_dynamicFriction = m_cfg.kDF;
	sbd->m_config.m_baumgarte = m_cfg.kVCF;
	sbd->m_config.m_pressure = m_cfg.kPR;
	sbd->m_config.m_aeroModel = this->m_cfg.aeromodel;
	sbd->m_config.m_lift = m_cfg.kLF;
	sbd->m_config.m_drag = m_cfg.kDG;
	sbd->m_config.m_positionIterations = m_cfg.piterations;
	sbd->m_config.m_driftIterations = m_cfg.diterations;
	sbd->m_config.m_clusterIterations = m_cfg.citerations;
	sbd->m_config.m_velocityIterations = m_cfg.viterations;
	sbd->m_config.m_maxVolume = m_cfg.maxvolume;
	sbd->m_config.m_damping = m_cfg.kDP;
	sbd->m_config.m_poseMatch = m_cfg.kMT;
	sbd->m_config.m_collisionFlags = m_cfg.collisions;
	sbd->m_config.m_volume = m_cfg.kVC;
	sbd->m_config.m_rigidContactHardness = m_cfg.kCHR;
	sbd->m_config.m_kineticContactHardness = m_cfg.kKHR;
	sbd->m_config.m_softContactHardness = m_cfg.kSHR;
	sbd->m_config.m_anchorHardness = m_cfg.kAHR;
	sbd->m_config.m_timeScale = m_cfg.timescale;
	sbd->m_config.m_maxVolume = m_cfg.maxvolume;
	sbd->m_config.m_softRigidClusterHardness = m_cfg.kSRHR_CL;
	sbd->m_config.m_softKineticClusterHardness = m_cfg.kSKHR_CL;
	sbd->m_config.m_softSoftClusterHardness = m_cfg.kSSHR_CL;
	sbd->m_config.m_softRigidClusterImpulseSplit = m_cfg.kSR_SPLT_CL;
	sbd->m_config.m_softKineticClusterImpulseSplit = m_cfg.kSK_SPLT_CL;
	sbd->m_config.m_softSoftClusterImpulseSplit = m_cfg.kSS_SPLT_CL;

	//pose for shape matching
	{
		sbd->m_pose = (SoftBodyPoseData*)serializer->getUniquePointer((void*)&m_pose);

		int sz = sizeof(SoftBodyPoseData);
		btChunk* chunk = serializer->allocate(sz, 1);
		SoftBodyPoseData* memPtr = (SoftBodyPoseData*)chunk->m_oldPtr;

		m_pose.m_aqq.serializeFloat(memPtr->m_aqq);
		memPtr->m_bframe = m_pose.m_bframe;
		memPtr->m_bvolume = m_pose.m_bvolume;
		m_pose.m_com.serializeFloat(memPtr->m_com);

		memPtr->m_numPositions = m_pose.m_pos.size();
		memPtr->m_positions = memPtr->m_numPositions ? (btVector3FloatData*)serializer->getUniquePointer((void*)&m_pose.m_pos[0]) : 0;
		if (memPtr->m_numPositions)
		{
			int numElem = memPtr->m_numPositions;
			int sz = sizeof(btVector3Data);
			btChunk* chunk = serializer->allocate(sz, numElem);
			btVector3FloatData* memPtr = (btVector3FloatData*)chunk->m_oldPtr;
			for (int i = 0; i < numElem; i++, memPtr++)
			{
				m_pose.m_pos[i].serializeFloat(*memPtr);
			}
			serializer->finalizeChunk(chunk, "btVector3FloatData", BT_ARRAY_CODE, (void*)&m_pose.m_pos[0]);
		}
		memPtr->m_restVolume = m_pose.m_volume;
		m_pose.m_rot.serializeFloat(memPtr->m_rot);
		m_pose.m_scl.serializeFloat(memPtr->m_scale);

		memPtr->m_numWeigts = m_pose.m_wgh.size();
		memPtr->m_weights = memPtr->m_numWeigts ? (float*)serializer->getUniquePointer((void*)&m_pose.m_wgh[0]) : 0;
		if (memPtr->m_numWeigts)
		{
			int numElem = memPtr->m_numWeigts;
			int sz = sizeof(float);
			btChunk* chunk = serializer->allocate(sz, numElem);
			float* memPtr = (float*)chunk->m_oldPtr;
			for (int i = 0; i < numElem; i++, memPtr++)
			{
				*memPtr = m_pose.m_wgh[i];
			}
			serializer->finalizeChunk(chunk, "float", BT_ARRAY_CODE, (void*)&m_pose.m_wgh[0]);
		}

		serializer->finalizeChunk(chunk, "SoftBodyPoseData", BT_ARRAY_CODE, (void*)&m_pose);
	}

	//clusters for convex-cluster collision detection

	sbd->m_numClusters = m_clusters.size();
	sbd->m_clusters = sbd->m_numClusters ? (SoftBodyClusterData*)serializer->getUniquePointer((void*)m_clusters[0]) : 0;
	if (sbd->m_numClusters)
	{
		int numElem = sbd->m_numClusters;
		int sz = sizeof(SoftBodyClusterData);
		btChunk* chunk = serializer->allocate(sz, numElem);
		SoftBodyClusterData* memPtr = (SoftBodyClusterData*)chunk->m_oldPtr;
		for (int i = 0; i < numElem; i++, memPtr++)
		{
			memPtr->m_adamping = m_clusters[i]->m_adamping;
			m_clusters[i]->m_av.serializeFloat(memPtr->m_av);
			memPtr->m_clusterIndex = m_clusters[i]->m_clusterIndex;
			memPtr->m_collide = m_clusters[i]->m_collide;
			m_clusters[i]->m_com.serializeFloat(memPtr->m_com);
			memPtr->m_containsAnchor = m_clusters[i]->m_containsAnchor;
			m_clusters[i]->m_dimpulses[0].serializeFloat(memPtr->m_dimpulses[0]);
			m_clusters[i]->m_dimpulses[1].serializeFloat(memPtr->m_dimpulses[1]);
			m_clusters[i]->m_framexform.serializeFloat(memPtr->m_framexform);
			memPtr->m_idmass = m_clusters[i]->m_idmass;
			memPtr->m_imass = m_clusters[i]->m_imass;
			m_clusters[i]->m_invwi.serializeFloat(memPtr->m_invwi);
			memPtr->m_ldamping = m_clusters[i]->m_ldamping;
			m_clusters[i]->m_locii.serializeFloat(memPtr->m_locii);
			m_clusters[i]->m_lv.serializeFloat(memPtr->m_lv);
			memPtr->m_matching = m_clusters[i]->m_matching;
			memPtr->m_maxSelfCollisionImpulse = m_clusters[i]->m_maxSelfCollisionImpulse;
			memPtr->m_ndamping = m_clusters[i]->m_ndamping;
			memPtr->m_ldamping = m_clusters[i]->m_ldamping;
			memPtr->m_adamping = m_clusters[i]->m_adamping;
			memPtr->m_selfCollisionImpulseFactor = m_clusters[i]->m_selfCollisionImpulseFactor;

			memPtr->m_numFrameRefs = m_clusters[i]->m_framerefs.size();
			memPtr->m_numMasses = m_clusters[i]->m_masses.size();
			memPtr->m_numNodes = m_clusters[i]->m_nodes.size();

			memPtr->m_nvimpulses = m_clusters[i]->m_nvimpulses;
			m_clusters[i]->m_vimpulses[0].serializeFloat(memPtr->m_vimpulses[0]);
			m_clusters[i]->m_vimpulses[1].serializeFloat(memPtr->m_vimpulses[1]);
			memPtr->m_ndimpulses = m_clusters[i]->m_ndimpulses;

			memPtr->m_framerefs = memPtr->m_numFrameRefs ? (btVector3FloatData*)serializer->getUniquePointer((void*)&m_clusters[i]->m_framerefs[0]) : 0;
			if (memPtr->m_framerefs)
			{
				int numElem = memPtr->m_numFrameRefs;
				int sz = sizeof(btVector3FloatData);
				btChunk* chunk = serializer->allocate(sz, numElem);
				btVector3FloatData* memPtr = (btVector3FloatData*)chunk->m_oldPtr;
				for (int j = 0; j < numElem; j++, memPtr++)
				{
					m_clusters[i]->m_framerefs[j].serializeFloat(*memPtr);
				}
				serializer->finalizeChunk(chunk, "btVector3FloatData", BT_ARRAY_CODE, (void*)&m_clusters[i]->m_framerefs[0]);
			}

			memPtr->m_masses = memPtr->m_numMasses ? (float*)serializer->getUniquePointer((void*)&m_clusters[i]->m_masses[0]) : 0;
			if (memPtr->m_masses)
			{
				int numElem = memPtr->m_numMasses;
				int sz = sizeof(float);
				btChunk* chunk = serializer->allocate(sz, numElem);
				float* memPtr = (float*)chunk->m_oldPtr;
				for (int j = 0; j < numElem; j++, memPtr++)
				{
					*memPtr = m_clusters[i]->m_masses[j];
				}
				serializer->finalizeChunk(chunk, "float", BT_ARRAY_CODE, (void*)&m_clusters[i]->m_masses[0]);
			}

			memPtr->m_nodeIndices = memPtr->m_numNodes ? (int*)serializer->getUniquePointer((void*)&m_clusters[i]->m_nodes) : 0;
			if (memPtr->m_nodeIndices)
			{
				int numElem = memPtr->m_numMasses;
				int sz = sizeof(int);
				btChunk* chunk = serializer->allocate(sz, numElem);
				int* memPtr = (int*)chunk->m_oldPtr;
				for (int j = 0; j < numElem; j++, memPtr++)
				{
					int* indexPtr = m_nodeIndexMap.find(m_clusters[i]->m_nodes[j]);
					btAssert(indexPtr);
					*memPtr = *indexPtr;
				}
				serializer->finalizeChunk(chunk, "int", BT_ARRAY_CODE, (void*)&m_clusters[i]->m_nodes);
			}
		}
		serializer->finalizeChunk(chunk, "SoftBodyClusterData", BT_ARRAY_CODE, (void*)m_clusters[0]);
	}

	sbd->m_numJoints = m_joints.size();
	sbd->m_joints = m_joints.size() ? (btSoftBodyJointData*)serializer->getUniquePointer((void*)&m_joints[0]) : 0;

	if (sbd->m_joints)
	{
		int sz = sizeof(btSoftBodyJointData);
		int numElem = m_joints.size();
		btChunk* chunk = serializer->allocate(sz, numElem);
		btSoftBodyJointData* memPtr = (btSoftBodyJointData*)chunk->m_oldPtr;

		for (int i = 0; i < numElem; i++, memPtr++)
		{
			memPtr->m_jointType = (int)m_joints[i]->Type();
			m_joints[i]->m_refs[0].serializeFloat(memPtr->m_refs[0]);
			m_joints[i]->m_refs[1].serializeFloat(memPtr->m_refs[1]);
			memPtr->m_cfm = m_joints[i]->m_cfm;
			memPtr->m_erp = float(m_joints[i]->m_erp);
			memPtr->m_split = float(m_joints[i]->m_split);
			memPtr->m_delete = m_joints[i]->m_delete;

			for (int j = 0; j < 4; j++)
			{
				memPtr->m_relPosition[0].m_floats[j] = 0.f;
				memPtr->m_relPosition[1].m_floats[j] = 0.f;
			}
			memPtr->m_bodyA = 0;
			memPtr->m_bodyB = 0;
			if (m_joints[i]->m_bodies[0].m_soft)
			{
				memPtr->m_bodyAtype = BT_JOINT_SOFT_BODY_CLUSTER;
				memPtr->m_bodyA = serializer->getUniquePointer((void*)m_joints[i]->m_bodies[0].m_soft);
			}
			if (m_joints[i]->m_bodies[0].m_collisionObject)
			{
				memPtr->m_bodyAtype = BT_JOINT_COLLISION_OBJECT;
				memPtr->m_bodyA = serializer->getUniquePointer((void*)m_joints[i]->m_bodies[0].m_collisionObject);
			}
			if (m_joints[i]->m_bodies[0].m_rigid)
			{
				memPtr->m_bodyAtype = BT_JOINT_RIGID_BODY;
				memPtr->m_bodyA = serializer->getUniquePointer((void*)m_joints[i]->m_bodies[0].m_rigid);
			}

			if (m_joints[i]->m_bodies[1].m_soft)
			{
				memPtr->m_bodyBtype = BT_JOINT_SOFT_BODY_CLUSTER;
				memPtr->m_bodyB = serializer->getUniquePointer((void*)m_joints[i]->m_bodies[1].m_soft);
			}
			if (m_joints[i]->m_bodies[1].m_collisionObject)
			{
				memPtr->m_bodyBtype = BT_JOINT_COLLISION_OBJECT;
				memPtr->m_bodyB = serializer->getUniquePointer((void*)m_joints[i]->m_bodies[1].m_collisionObject);
			}
			if (m_joints[i]->m_bodies[1].m_rigid)
			{
				memPtr->m_bodyBtype = BT_JOINT_RIGID_BODY;
				memPtr->m_bodyB = serializer->getUniquePointer((void*)m_joints[i]->m_bodies[1].m_rigid);
			}
		}
		serializer->finalizeChunk(chunk, "btSoftBodyJointData", BT_ARRAY_CODE, (void*)&m_joints[0]);
	}

	return btSoftBodyDataName;
}

void btSoftBody::updateDeactivation(btScalar timeStep)
{
	if ((getActivationState() == ISLAND_SLEEPING) || (getActivationState() == DISABLE_DEACTIVATION))
		return;

	// If the soft body is attached to a rigid body which is active, we keep this soft also active, otherwise the anchors would
	// not be considered
	for (auto i = 0; i < m_deformableAnchors.size(); ++i)
	{
		const auto& rigid_body = m_deformableAnchors[i].m_body;
		if (rigid_body && rigid_body->isActive())
		{
			return;
		}
	}

	if (m_maxSpeedSquared < m_sleepingThreshold * m_sleepingThreshold)
	{
		m_deactivationTime += timeStep;
	}
	else
	{
		m_deactivationTime = btScalar(0.);
		setActivationState(0);
	}
}

void btSoftBody::setZeroVelocity()
{
	for (int i = 0; i < m_nodes.size(); ++i)
	{
		m_nodes[i].m_v.setZero();
	}
}

bool btSoftBody::wantsSleeping()
{
	if (getActivationState() == DISABLE_DEACTIVATION)
		return false;

	//disable deactivation
	if (gDisableDeactivation || (gDeactivationTime == btScalar(0.)))
		return false;

	if ((getActivationState() == ISLAND_SLEEPING) || (getActivationState() == WANTS_DEACTIVATION))
		return true;

	if (m_deactivationTime > gDeactivationTime)
	{
		return true;
	}
	return false;
}

// Surrounding nodes also, because my CD phase is not bulletproof (a fully tunnelled through triangle could cause an update here which should not happen).
// So this is a penetrating node cloud border growth. It is better to be conservative here and update only the nodes I am sure they do not penetrate,
// otherwise a position could be marked as safe even if it isn't - that is very bad. Border growth of 1 should be enough. That will cover the riskiest
// zone. Further nodes will be sanitized using the normal plane projection in updateLastSafeWorldTransform.
void btSoftBody::lastSafeBorderGrow(int growth, std::map<btSoftBody::Node*, StuckTetraIndicesMapped>& nodesInCollision)
{
	for (auto pass = 0; pass < growth; ++pass)
	{
		std::map<btSoftBody::Node*, StuckTetraIndicesMapped> grownNodesInCollision = nodesInCollision;
		for (auto& [node, mapped] : nodesInCollision)
			if (mapped.penetrating)
			{
				for (auto m = 0; m < node->m_tetraMembershipCount; ++m)
				{
					auto tetraMembershipIndex = node->m_tetraMembership[m];
					for (auto j = 0; j < 4; ++j)
					{
						auto iter = grownNodesInCollision.find(m_tetras[tetraMembershipIndex].m_n[j]);
						if (iter == grownNodesInCollision.end())
							grownNodesInCollision.insert({m_tetras[tetraMembershipIndex].m_n[j], mapped});
						else
							iter->second.penetrating = true;
					}
				}
			}
		nodesInCollision = grownNodesInCollision;
	}
}

int kLastSafeGrowth = 1;

// Notes about the OGC paper (Offset Geometric Contact - https://dl.acm.org/doi/10.1145/3731205). I was trying to find if the ideas in the paper were applicable to the
// current soft body simulation implemented here. There are some key differences which make it hard. The most useful idea for me is the penetration free simulation.
// They calculate a safe circular bound for every simulation vertex (they do not have a separate collision mesh - they are mostly dealing with cloths and strings, so it is ok for them)
// where the vertex can move in a simulation iteration. To find this bound for one vertex, they do not use an overlap query, but they use a closest distance query, so
// this enables them to "see" obstacles no matter how far away from them they are, so there can not be any tunnelling problems. I would have to do (because I have a separate collision mesh)
// a closest distance query for every triangle of the collision mesh. Then, for one simulation vertex, I would go through all tetras which share that vertex and through all triangles which
// pass through those tetras. I would gather all the min distances of those triangles and assign the bound for the sim vertex to be the minimum of those min distances.
// This would have to be done in multiple simulation iterations (!) as in the paper, they also do the collision detection per iteration (not always - see the collisionDetectionRequired
// variable in their pseudocode). I have a strong feeling that it would be a very hard performance hit given the kind of hi poly meshes I use. On GPU, I hazard a guess
// that it could be doable because in the paper, they have this 50 cloth layers case with 1M vertices (1.96M triangles) and the timestep takes 6-11 ms on GPU. My
// recently used hi poly mesh has 300K triangles, so it should not be so bad, but we have to account for the fact that triangle-triangle nearest query would be more expensive
// than theirs vertex-triangle nearest query. Problem is that running on GPU would be very time consuming to implement now.
// Another interesting idea are the offset meshes constructed by replacing every triangle by a prism, every edge by a capsule and every vertex by a sphere. I currently only do
// "replacement" of every triangle by a capsule-like triangle (simply a tri-tri distance query with some margin tolerance). This does not give a perfect contact normal in all cases.

void btSoftBody::updateLastSafeWorldTransform()
{
	// Note that unlike btSoftBody::applyLastSafeWorldTransform, this can not be disabled,
	// because the penetration contact generation in btPrimitiveTriangle::find_triangle_collision_alt_method_outer
	// also depends on this last safe data.

	// updateLastSafeWorldTransform also used to be partially updated, but it turned out to be a nightmare to keep the safe positions truly safe like that.
	// Miniscule (and hard to diagnose the reason - this is why the BT_SAFE_UPDATE_DEBUG exists) errors creeped in causing the safe positions to be unsafe rendering them unusable and causing
	// explosions. So now we always do a full update of all nodes.

	for (auto i = 0; i < m_nodes.size(); ++i)
	{
		const auto& node = m_nodes[i];
		auto& safe = m_nodes[i].m_safe;
#ifdef BT_SAFE_UPDATE_DEBUG
		fprintf(stderr, "drawpoint \"whole pt old safe\" [%f,%f,%f][1,1,1,1] \n", safe.m_x.x(), safe.m_x.y(), safe.m_x.z());
#endif
		safe.m_x = node.m_x;
		safe.m_q = node.m_q;

#ifdef BT_SAFE_UPDATE_DEBUG
		fprintf(stderr, "drawpoint \"whole pt\" [%f,%f,%f][1,1,1,1] \n", node.m_x.x(), node.m_x.y(), node.m_x.z());
#endif
		/*dst.m_q = src.m_q;
		dst.m_v = src.m_v;
		dst.m_vn = src.m_vn;
		dst.m_f = src.m_f;
		dst.m_n = src.m_n;
		dst.m_im = src.m_im;
		dst.m_area = src.m_area;
		dst.m_splitv = src.m_splitv;
		dst.m_effectiveMass = src.m_effectiveMass;
		dst.m_effectiveMass_inv = src.m_effectiveMass_inv;*/
	}
}

void btSoftBody::applyLastSafeWorldTransform(const std::map<int, StuckTetraIndicesMapped>* partial)
{
	// If it turns out that there are some situations where penetrations have to be resolved in a single step:
	// it is implemented on this branch in bullet3 repo recursive_safe_apply and on VRUT-6981_recursive_safe_apply

	if (getCollisionFlags() & CF_APPLY_LAST_SAFE)
	{
		std::map<btSoftBody::Node*, std::vector<btVector3>> nodesInCollision;
		if (partial)
		{
			for (auto& [partTetraIndex, partMapped] : *partial)
			{
				if (partMapped.penetrating)
					for (auto i = 0; i < 4; ++i)
					{
						auto iter = nodesInCollision.find(m_tetras[partTetraIndex].m_n[i]);
						if (iter == nodesInCollision.end())
							nodesInCollision.insert({m_tetras[partTetraIndex].m_n[i], partMapped.opposingNormals});
						else
							iter->second.insert(iter->second.end(), partMapped.opposingNormals.begin(), partMapped.opposingNormals.end());
					}
			}
		}
		else
		{
			btAssert(false);  // Softs should always be partial now
			for (int i = 0; i < m_nodes.size(); ++i)
			{
				nodesInCollision.insert({&m_nodes[i], {}});
			}
		}

		for (auto& [nodeInCollision, normals] : nodesInCollision)
		{
			const auto& src = nodeInCollision->m_safe;
			auto& dst = nodeInCollision;

#ifdef BT_SAFE_UPDATE_DEBUG
			fprintf(stderr, "drawpoint \"apply pt\" [%f,%f,%f][0,0,1,1] \n", src.m_x.x(), src.m_x.y(), src.m_x.z());
			fprintf(stderr, "drawline \"apply ln\" [%f,%f,%f][%f,%f,%f][0,0,1,1] \n", dst->m_x.x(), dst->m_x.y(), dst->m_x.z(), src.m_x.x(), src.m_x.y(), src.m_x.z());
#endif

			dst->m_x = src.m_x;
			dst->m_q = src.m_q;

			// Velocities are modified (projected), so that they do not cause the penetration again
			for (const auto& normal : normals)
			{
				if (dst->m_v.dot(normal) < 0.0)
					dst->m_v = dst->m_v.rejectFrom(normal);

				if (dst->m_vn.dot(normal) < 0.0)
					dst->m_vn = dst->m_vn.rejectFrom(normal);
			}
		}
	}
}
/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/

///btSoftBody implementation by Nathanael Presson

#ifndef _BT_SOFT_BODY_H
#define _BT_SOFT_BODY_H

#include "LinearMath/btAlignedObjectArray.h"
#include "LinearMath/btTransform.h"
#include "LinearMath/btIDebugDraw.h"
#include "LinearMath/btVector3.h"
#include "BulletDynamics/Dynamics/btRigidBody.h"

#include "BulletCollision/CollisionShapes/btConcaveShape.h"
#include "BulletCollision/CollisionDispatch/btCollisionCreateFunc.h"
#include "btSparseSDF.h"
#include "BulletCollision/BroadphaseCollision/btDbvt.h"
#include "BulletDynamics/Featherstone/btMultiBodyLinkCollider.h"
#include "BulletDynamics/Featherstone/btMultiBodyConstraint.h"

#include <memory>
#include <vector>

//#ifdef BT_USE_DOUBLE_PRECISION
//#define btRigidBodyData	btRigidBodyDoubleData
//#define btRigidBodyDataName	"btRigidBodyDoubleData"
//#else
#define btSoftBodyData btSoftBodyFloatData
#define btSoftBodyDataName "btSoftBodyFloatData"
static const btScalar OVERLAP_REDUCTION_FACTOR = 0.1;
static unsigned long seed = 243703;
//#endif //BT_USE_DOUBLE_PRECISION

class btBroadphaseInterface;
class btDispatcher;
class btSoftBodySolver;

/* btSoftBodyWorldInfo	*/
struct btSoftBodyWorldInfo
{
	btScalar air_density;
	btScalar water_density;
	btScalar water_offset;
	btScalar m_maxDisplacement;
	btVector3 water_normal;
	btBroadphaseInterface* m_broadphase;
	btDispatcher* m_dispatcher;
	btVector3 m_gravity;
	btSparseSdf<3> m_sparsesdf;

	btSoftBodyWorldInfo()
		: air_density((btScalar)1.2),
		  water_density(0),
		  water_offset(0),
		  m_maxDisplacement(1000.f),  //avoid soft body from 'exploding' so use some upper threshold of maximum motion that a node can travel per frame
		  water_normal(0, 0, 0),
		  m_broadphase(0),
		  m_dispatcher(0),
		  m_gravity(0, -10, 0)
	{
	}
};

///The btSoftBody is an class to simulate cloth and volumetric soft bodies.
///There is two-way interaction between btSoftBody and btRigidBody/btCollisionObject.
class btSoftBody : public btCollisionObject
{
	static constexpr btScalar influencedNodesFactor = 0.0;  // At least 0.05 for neohookean. It seems that linear elasticity is ok with just one node being inluenced.

public:
	btAlignedObjectArray<const class btCollisionObject*> m_collisionDisabledObjects;

	// The solver object that handles this soft body
	btSoftBodySolver* m_softBodySolver;

	//
	// Enumerations
	//

	///eAeroModel
	struct eAeroModel
	{
		enum _
		{
			V_Point,             ///Vertex normals are oriented toward velocity
			V_TwoSided,          ///Vertex normals are flipped to match velocity
			V_TwoSidedLiftDrag,  ///Vertex normals are flipped to match velocity and lift and drag forces are applied
			V_OneSided,          ///Vertex normals are taken as it is
			F_TwoSided,          ///Face normals are flipped to match velocity
			F_TwoSidedLiftDrag,  ///Face normals are flipped to match velocity and lift and drag forces are applied
			F_OneSided,          ///Face normals are taken as it is
			END
		};
	};

	///eVSolver : velocities solvers
	struct eVSolver
	{
		enum _
		{
			Linear,  ///Linear solver
			END
		};
	};

	///ePSolver : positions solvers
	struct ePSolver
	{
		enum _
		{
			Linear,     ///Linear solver
			Anchors,    ///Anchor solver
			RContacts,  ///Rigid contacts solver
			SContacts,  ///Soft contacts solver
			END
		};
	};

	///eSolverPresets
	struct eSolverPresets
	{
		enum _
		{
			Positions,
			Velocities,
			Default = Positions,
			END
		};
	};

	///eFeature
	struct eFeature
	{
		enum _
		{
			None,
			Node,
			Link,
			Face,
			Tetra,
			END
		};
	};

	typedef btAlignedObjectArray<eVSolver::_> tVSolverArray;
	typedef btAlignedObjectArray<ePSolver::_> tPSolverArray;

	//
	// Flags
	//

	///fCollision
	struct fCollision
	{
		enum _
		{
			RVSmask = 0x000f,  ///Rigid versus soft mask
			SDF_RS = 0x0001,   ///SDF based rigid vs soft
			CL_RS = 0x0002,    ///Cluster vs convex rigid vs soft
			SDF_RD = 0x0004,   ///rigid vs deformable

			SVSmask = 0x00f0,  ///Rigid versus soft mask
			VF_SS = 0x0010,    ///Vertex vs face soft vs soft handling
			CL_SS = 0x0020,    ///Cluster vs cluster soft vs soft handling
			CL_SELF = 0x0040,  ///Cluster soft body self collision
			VF_DD = 0x0080,    ///Vertex vs face soft vs soft handling

			RVDFmask = 0x0f00,  /// Rigid versus deformable face mask
			SDF_RDF = 0x0100,   /// GJK based Rigid vs. deformable face
			SDF_MDF = 0x0200,   /// GJK based Multibody vs. deformable face
			SDF_RDN = 0x0400,   /// SDF based Rigid vs. deformable node
			/* presets	*/
			Default = SDF_RS,
			END
		};
	};

	///fMaterial
	struct fMaterial
	{
		enum _
		{
			DebugDraw = 0x0001,  /// Enable debug draw
			/* presets	*/
			Default = DebugDraw,
			END
		};
	};

	//
	// API Types
	//

	/* sRayCast		*/
	struct sRayCast
	{
		btSoftBody* body;     /// soft body
		eFeature::_ feature;  /// feature type
		int index;            /// feature index
		btScalar fraction;    /// time of impact fraction (rayorg+(rayto-rayfrom)*fraction)
	};

	/* ImplicitFn	*/
	struct ImplicitFn
	{
		virtual ~ImplicitFn() {}
		virtual btScalar Eval(const btVector3& x) = 0;
	};

	//
	// Internal types
	//

	typedef btAlignedObjectArray<btScalar> tScalarArray;
	typedef btAlignedObjectArray<btVector3> tVector3Array;

	/* sCti is Softbody contact info	*/
	struct sCti
	{
		const btCollisionObject* m_colObj; /* Rigid body			        */
		btVector3 m_normal;                /* Outward normal		        */
		mutable btVector3 m_impulse;       /* Applied impulse        	    */
		btScalar m_offset;                 /* Offset from origin	        */
		btVector3 m_bary;                  /* Barycentric weights for faces */
		int m_count = 0;
		mutable btScalar* m_contact_point_impulse_magnitude = nullptr;
		sCti() : m_impulse(0, 0, 0) {}
	};

	/* sMedium		*/
	struct sMedium
	{
		btVector3 m_velocity; /* Velocity				*/
		btScalar m_pressure;  /* Pressure				*/
		btScalar m_density;   /* Density				*/
	};

	/* Base type	*/
	struct Element
	{
		void* m_tag;  // User data
		Element() : m_tag(0) {}
	};
	/* Material		*/
	struct Material : Element
	{
		btScalar m_kLST;  // Linear stiffness coefficient [0,1]
		btScalar m_kAST;  // Area/Angular stiffness coefficient [0,1]
		btScalar m_kVST;  // Volume stiffness coefficient [0,1]
		int m_flags;      // Flags
	};

	/* Feature		*/
	struct Feature : Element
	{
		Material* m_material;  // Material
	};
	/* Node			*/
	struct RenderNode
	{
		btVector3 m_x;
		btVector3 m_uv1;
		btVector3 m_normal;
	};

	struct Node;

	struct NodeSafe
	{
		btVector3 m_x;                    // Safe position when there was no penetration. Used in TrimeshDeformedPrimitiveManager.
		btVector3 m_q;                    // Previous step position/Test position
		btVector3 m_v;                    // Velocity
		btVector3 m_vn;                   // Previous step velocity
		btVector3 m_f;                    // Force accumulator
		btVector3 m_n;                    // Normal
		btScalar m_im;                    // 1/mass
		btScalar m_area;                  // Area
		btVector3 m_splitv;               // velocity associated with split impulse
		btMatrix3x3 m_effectiveMass;      // effective mass in contact
		btMatrix3x3 m_effectiveMass_inv;  // inverse of effective mass

		void copy_from_node(Node& node)
		{
			m_x = node.m_x;
			m_q = node.m_q;
			m_v = node.m_v;
			m_vn = node.m_vn;
			m_f = node.m_f;
			m_n = node.m_n;
			m_im = node.m_im;
			m_area = node.m_area;
			m_splitv = node.m_splitv;
			m_effectiveMass = node.m_effectiveMass;
			m_effectiveMass_inv = node.m_effectiveMass_inv;
		}
	};

	struct Node : Feature
	{
		btVector3 m_x;       // Position
		btVector3 m_q;       // Previous step position/Test position
		btVector3 m_v;       // Velocity
		btVector3 m_vn;      // Previous step velocity
		btVector3 m_f;       // Force accumulator
		btVector3 m_n;       // Normal
		btScalar m_im;       // 1/mass
		btScalar m_area;     // Area
		btDbvtNode* m_leaf;  // Leaf data
		int m_constrained;   // depth of penetration
		int m_battach : 1;   // Attached
		int index;           // Watch out - there is a certain design smell here. Index is sometimes an index into the global m_nodes array (when using the deformable soft bodies) and sometimes into the soft body's own m_nodes array. This is why the local_index was introduced. To reliably have an index into the local m_nodes.
		int local_index;
		btVector3 m_splitv;               // velocity associated with split impulse
		btMatrix3x3 m_effectiveMass;      // effective mass in contact
		btMatrix3x3 m_effectiveMass_inv;  // inverse of effective mass
		NodeSafe m_safe;
		static constexpr int kMaxTetraMembershipPerNode = 32;
		int m_tetraMembership[kMaxTetraMembershipPerNode];
		int m_tetraMembershipCount;
	};
	/* Link			*/
	ATTRIBUTE_ALIGNED16(struct)
	Link : Feature
	{
		btVector3 m_c3;      // gradient
		Node* m_n[2];        // Node pointers
		btScalar m_rl;       // Rest length
		int m_bbending : 1;  // Bending link
		btScalar m_c0;       // (ima+imb)*kLST
		btScalar m_c1;       // rl^2
		btScalar m_c2;       // |gradient|^2/c0

		BT_DECLARE_ALIGNED_ALLOCATOR();
	};
	struct RenderFace
	{
		RenderNode* m_n[3];  // Node pointers
	};

	/* Face			*/
	struct Face : Feature
	{
		Node* m_n[3];          // Node pointers
		btVector3 m_normal;    // Normal
		btScalar m_ra;         // Rest area
		btDbvtNode* m_leaf;    // Leaf data
		btVector4 m_pcontact;  // barycentric weights of the persistent contact
		btVector3 m_n0, m_n1, m_vn;
		int m_index;
	};
	/* Tetra		*/
	struct Tetra : Feature
	{
		Node* m_n[4];              // Node pointers
		btScalar m_rv;             // Rest volume
		btDbvtNode* m_leaf;        // Leaf data
		btVector3 m_c0[4];         // gradients
		btScalar m_c1;             // (4*kVST)/(im0+im1+im2+im3)
		btScalar m_c2;             // m_c1/sum(|g0..3|^2)
		btMatrix3x3 m_Dm_inverse;  // rest Dm^-1
		btMatrix3x3 m_F;
		btScalar m_element_measure;
		btVector4 m_P_inv[3];  // first three columns of P_inv matrix
	};

	/*  TetraScratch  */
	struct TetraScratch
	{
		btMatrix3x3 m_F;           // deformation gradient F
		btScalar m_trace;          // trace of F^T * F
		btScalar m_J;              // det(F)
		btMatrix3x3 m_cofF;        // cofactor of F
		btMatrix3x3 m_corotation;  // corotatio of the tetra
	};

	/* RContact		*/
	struct RContact
	{
		sCti m_cti;        // Contact infos
		Node* m_node;      // Owner node
		btMatrix3x3 m_c0;  // Impulse matrix
		btVector3 m_c1;    // Relative anchor
		btScalar m_c2;     // ima*dt
		btScalar m_c3;     // Friction
		btScalar m_c4;     // Hardness

		// jacobians and unit impulse responses for multibody
		btMultiBodyJacobianData jacobianData_normal;
		btMultiBodyJacobianData jacobianData_t1;
		btMultiBodyJacobianData jacobianData_t2;
		btVector3 t1;
		btVector3 t2;
	};

	class DeformableRigidContact
	{
	public:
		sCti m_cti;        // Contact infos
		btMatrix3x3 m_c0;  // Impulse matrix
		btVector3 m_c1;    // Relative anchor
		btScalar m_c2;     // inverse mass of node/face
		btScalar m_c3;     // Friction
		btScalar m_c4;     // Hardness
		btMatrix3x3 m_c5;  // inverse effective mass

		// jacobians and unit impulse responses for multibody
		btMultiBodyJacobianData jacobianData_normal;
		btMultiBodyJacobianData jacobianData_t1;
		btMultiBodyJacobianData jacobianData_t2;
		btVector3 t1;
		btVector3 t2;
	};

	class DeformableNodeRigidContact : public DeformableRigidContact
	{
	public:
		Node* m_node;  // Owner node
	};

	class DeformableNodeRigidAnchor : public DeformableNodeRigidContact
	{
	public:
		btVector3 m_local;  // Anchor position in body space
		uint32_t m_userIndex;
		btRigidBody* m_body = nullptr;  // Body
	};

	class DeformableFaceRigidContact : public DeformableRigidContact
	{
	public:
		Face* m_face;              // Owner face
		btVector3 m_contactPoint;  // Contact point
		btVector3 m_bary;          // Barycentric weights
		btVector3 m_weights;       // v_contactPoint * m_weights[i] = m_face->m_node[i]->m_v;
	};

	struct DeformableFaceNodeContact
	{
		Node* m_node;                       // Node
		Face* m_face;                       // Face
		btVector3 m_bary;                   // Barycentric weights
		btVector3 m_weights;                // v_contactPoint * m_weights[i] = m_face->m_node[i]->m_v;
		btVector3 m_normal;                 // Normal
		btScalar m_margin;                  // Margin
		btScalar m_friction;                // Friction
		btScalar m_imf;                     // inverse mass of the face at contact point
		btScalar m_c0;                      // scale of the impulse matrix;
		const btCollisionObject* m_colObj;  // Collision object to collide with.
	};

	struct DeformableNodeNodeContact
	{
		Node* m_node0;       // Node0
		Node* m_node1;       // Node1
		btVector3 m_normal;  // Normal
		btScalar m_offset;
		btScalar m_friction;                // Friction
		const btCollisionObject* m_colObj;  // Collision object to collide with.
		btScalar* m_contact_point_impulse_magnitude = nullptr;
		int m_count = 0;
	};

	/* SContact		*/
	struct SContact
	{
		Node* m_node;         // Node
		Face* m_face;         // Face
		btVector3 m_weights;  // Weigths
		btVector3 m_normal;   // Normal
		btScalar m_margin;    // Margin
		btScalar m_friction;  // Friction
		btScalar m_cfm[2];    // Constraint force mixing
	};
	/* Anchor		*/
	struct Anchor
	{
		Node* m_node;         // Node pointer
		btVector3 m_local;    // Anchor position in body space
		btRigidBody* m_body;  // Body
		btScalar m_influence;
		btMatrix3x3 m_c0;  // Impulse matrix
		btVector3 m_c1;    // Relative anchor
		btScalar m_c2;     // ima*dt
	};
	/* Note			*/
	struct Note : Element
	{
		const char* m_text;    // Text
		btVector3 m_offset;    // Offset
		int m_rank;            // Rank
		Node* m_nodes[4];      // Nodes
		btScalar m_coords[4];  // Coordinates
	};
	/* Pose			*/
	struct Pose
	{
		bool m_bvolume;       // Is valid
		bool m_bframe;        // Is frame
		btScalar m_volume;    // Rest volume
		tVector3Array m_pos;  // Reference positions
		tScalarArray m_wgh;   // Weights
		btVector3 m_com;      // COM
		btMatrix3x3 m_rot;    // Rotation
		btMatrix3x3 m_scl;    // Scale
		btMatrix3x3 m_aqq;    // Base scaling
	};
	/* Cluster		*/
	struct Cluster
	{
		tScalarArray m_masses;
		btAlignedObjectArray<Node*> m_nodes;
		tVector3Array m_framerefs;
		btTransform m_framexform;
		btScalar m_idmass;
		btScalar m_imass;
		btMatrix3x3 m_locii;
		btMatrix3x3 m_invwi;
		btVector3 m_com;
		btVector3 m_vimpulses[2];
		btVector3 m_dimpulses[2];
		int m_nvimpulses;
		int m_ndimpulses;
		btVector3 m_lv;
		btVector3 m_av;
		btDbvtNode* m_leaf;
		btScalar m_ndamping; /* Node damping		*/
		btScalar m_ldamping; /* Linear damping	*/
		btScalar m_adamping; /* Angular damping	*/
		btScalar m_matching;
		btScalar m_maxSelfCollisionImpulse;
		btScalar m_selfCollisionImpulseFactor;
		bool m_containsAnchor;
		bool m_collide;
		int m_clusterIndex;
		Cluster() : m_leaf(0), m_ndamping(0), m_ldamping(0), m_adamping(0), m_matching(0), m_maxSelfCollisionImpulse(100.f), m_selfCollisionImpulseFactor(0.01f), m_containsAnchor(false)
		{
		}
	};
	/* Impulse		*/
	struct Impulse
	{
		btVector3 m_velocity;
		btVector3 m_drift;
		int m_asVelocity : 1;
		int m_asDrift : 1;
		Impulse() : m_velocity(0, 0, 0), m_drift(0, 0, 0), m_asVelocity(0), m_asDrift(0) {}
		Impulse operator-() const
		{
			Impulse i = *this;
			i.m_velocity = -i.m_velocity;
			i.m_drift = -i.m_drift;
			return (i);
		}
		Impulse operator*(btScalar x) const
		{
			Impulse i = *this;
			i.m_velocity *= x;
			i.m_drift *= x;
			return (i);
		}
	};
	/* Body			*/
	struct Body
	{
		Cluster* m_soft;
		btRigidBody* m_rigid;
		const btCollisionObject* m_collisionObject;

		Body() : m_soft(0), m_rigid(0), m_collisionObject(0) {}
		Body(Cluster* p) : m_soft(p), m_rigid(0), m_collisionObject(0) {}
		Body(const btCollisionObject* colObj) : m_soft(0), m_collisionObject(colObj)
		{
			m_rigid = (btRigidBody*)btRigidBody::upcast(m_collisionObject);
		}

		void activate() const
		{
			if (m_rigid)
				m_rigid->activate();
			if (m_collisionObject)
				m_collisionObject->activate();
		}
		const btMatrix3x3& invWorldInertia() const
		{
			static const btMatrix3x3 iwi(0, 0, 0, 0, 0, 0, 0, 0, 0);
			if (m_rigid) return (m_rigid->getInvInertiaTensorWorld());
			if (m_soft) return (m_soft->m_invwi);
			return (iwi);
		}
		btScalar invMass() const
		{
			if (m_rigid) return (m_rigid->getInvMass());
			if (m_soft) return (m_soft->m_imass);
			return (0);
		}
		const btTransform& xform() const
		{
			static const btTransform identity = btTransform::getIdentity();
			if (m_collisionObject) return (m_collisionObject->getWorldTransform());
			if (m_soft) return (m_soft->m_framexform);
			return (identity);
		}
		btVector3 linearVelocity() const
		{
			if (m_rigid) return (m_rigid->getLinearVelocity());
			if (m_soft) return (m_soft->m_lv);
			return (btVector3(0, 0, 0));
		}
		btVector3 angularVelocity(const btVector3& rpos) const
		{
			if (m_rigid) return (btCross(m_rigid->getAngularVelocity(), rpos));
			if (m_soft) return (btCross(m_soft->m_av, rpos));
			return (btVector3(0, 0, 0));
		}
		btVector3 angularVelocity() const
		{
			if (m_rigid) return (m_rigid->getAngularVelocity());
			if (m_soft) return (m_soft->m_av);
			return (btVector3(0, 0, 0));
		}
		btVector3 velocity(const btVector3& rpos) const
		{
			return (linearVelocity() + angularVelocity(rpos));
		}
		void applyVImpulse(const btVector3& impulse, const btVector3& rpos) const
		{
			if (m_rigid) m_rigid->applyImpulse(impulse, rpos);
			if (m_soft) btSoftBody::clusterVImpulse(m_soft, rpos, impulse);
		}
		void applyDImpulse(const btVector3& impulse, const btVector3& rpos) const
		{
			if (m_rigid) m_rigid->applyImpulse(impulse, rpos);
			if (m_soft) btSoftBody::clusterDImpulse(m_soft, rpos, impulse);
		}
		void applyImpulse(const Impulse& impulse, const btVector3& rpos) const
		{
			if (impulse.m_asVelocity)
			{
				//				printf("impulse.m_velocity = %f,%f,%f\n",impulse.m_velocity.getX(),impulse.m_velocity.getY(),impulse.m_velocity.getZ());
				applyVImpulse(impulse.m_velocity, rpos);
			}
			if (impulse.m_asDrift)
			{
				//				printf("impulse.m_drift = %f,%f,%f\n",impulse.m_drift.getX(),impulse.m_drift.getY(),impulse.m_drift.getZ());
				applyDImpulse(impulse.m_drift, rpos);
			}
		}
		void applyVAImpulse(const btVector3& impulse) const
		{
			if (m_rigid) m_rigid->applyTorqueImpulse(impulse);
			if (m_soft) btSoftBody::clusterVAImpulse(m_soft, impulse);
		}
		void applyDAImpulse(const btVector3& impulse) const
		{
			if (m_rigid) m_rigid->applyTorqueImpulse(impulse);
			if (m_soft) btSoftBody::clusterDAImpulse(m_soft, impulse);
		}
		void applyAImpulse(const Impulse& impulse) const
		{
			if (impulse.m_asVelocity) applyVAImpulse(impulse.m_velocity);
			if (impulse.m_asDrift) applyDAImpulse(impulse.m_drift);
		}
		void applyDCImpulse(const btVector3& impulse) const
		{
			if (m_rigid) m_rigid->applyCentralImpulse(impulse);
			if (m_soft) btSoftBody::clusterDCImpulse(m_soft, impulse);
		}
	};
	/* Joint		*/
	struct Joint
	{
		struct eType
		{
			enum _
			{
				Linear = 0,
				Angular,
				Contact
			};
		};
		struct Specs
		{
			Specs() : erp(1), cfm(1), split(1) {}
			btScalar erp;
			btScalar cfm;
			btScalar split;
		};
		Body m_bodies[2];
		btVector3 m_refs[2];
		btScalar m_cfm;
		btScalar m_erp;
		btScalar m_split;
		btVector3 m_drift;
		btVector3 m_sdrift;
		btMatrix3x3 m_massmatrix;
		bool m_delete;
		virtual ~Joint() {}
		Joint() : m_delete(false) {}
		virtual void Prepare(btScalar dt, int iterations);
		virtual void Solve(btScalar dt, btScalar sor) = 0;
		virtual void Terminate(btScalar dt) = 0;
		virtual eType::_ Type() const = 0;
	};
	/* LJoint		*/
	struct LJoint : Joint
	{
		struct Specs : Joint::Specs
		{
			btVector3 position;
		};
		btVector3 m_rpos[2];
		void Prepare(btScalar dt, int iterations);
		void Solve(btScalar dt, btScalar sor);
		void Terminate(btScalar dt);
		eType::_ Type() const { return (eType::Linear); }
	};
	/* AJoint		*/
	struct AJoint : Joint
	{
		struct IControl
		{
			virtual ~IControl() {}
			virtual void Prepare(AJoint*) {}
			virtual btScalar Speed(AJoint*, btScalar current) { return (current); }
			static IControl* Default()
			{
				static IControl def;
				return (&def);
			}
		};
		struct Specs : Joint::Specs
		{
			Specs() : icontrol(IControl::Default()) {}
			btVector3 axis;
			IControl* icontrol;
		};
		btVector3 m_axis[2];
		IControl* m_icontrol;
		void Prepare(btScalar dt, int iterations);
		void Solve(btScalar dt, btScalar sor);
		void Terminate(btScalar dt);
		eType::_ Type() const { return (eType::Angular); }
	};
	/* CJoint		*/
	struct CJoint : Joint
	{
		int m_life;
		int m_maxlife;
		btVector3 m_rpos[2];
		btVector3 m_normal;
		btScalar m_friction;
		void Prepare(btScalar dt, int iterations);
		void Solve(btScalar dt, btScalar sor);
		void Terminate(btScalar dt);
		eType::_ Type() const { return (eType::Contact); }
	};
	/* Config		*/
	struct Config
	{
		eAeroModel::_ aeromodel;    // Aerodynamic model (default: V_Point)
		btScalar kVCF;              // Velocities correction factor (Baumgarte)
		btScalar kDP;               // Damping coefficient [0,1]
		btScalar kDG;               // Drag coefficient [0,+inf]
		btScalar kLF;               // Lift coefficient [0,+inf]
		btScalar kPR;               // Pressure coefficient [-inf,+inf]
		btScalar kVC;               // Volume conversation coefficient [0,+inf]
		btScalar kDF;               // Dynamic friction coefficient [0,1]
		btScalar kMT;               // Pose matching coefficient [0,1]
		btScalar kCHR;              // Rigid contacts hardness [0,1]
		btScalar kKHR;              // Kinetic contacts hardness [0,1]
		btScalar kSHR;              // Soft contacts hardness [0,1]
		btScalar kAHR;              // Anchors hardness [0,1]
		btScalar kSRHR_CL;          // Soft vs rigid hardness [0,1] (cluster only)
		btScalar kSKHR_CL;          // Soft vs kinetic hardness [0,1] (cluster only)
		btScalar kSSHR_CL;          // Soft vs soft hardness [0,1] (cluster only)
		btScalar kSR_SPLT_CL;       // Soft vs rigid impulse split [0,1] (cluster only)
		btScalar kSK_SPLT_CL;       // Soft vs rigid impulse split [0,1] (cluster only)
		btScalar kSS_SPLT_CL;       // Soft vs rigid impulse split [0,1] (cluster only)
		btScalar maxvolume;         // Maximum volume ratio for pose
		btScalar timescale;         // Time scale
		int viterations;            // Velocities solver iterations
		int piterations;            // Positions solver iterations
		int diterations;            // Drift solver iterations
		int citerations;            // Cluster solver iterations
		int collisions;             // Collisions flags
		tVSolverArray m_vsequence;  // Velocity solvers sequence
		tPSolverArray m_psequence;  // Position solvers sequence
		tPSolverArray m_dsequence;  // Drift solvers sequence
		btScalar drag;              // deformable air drag
		btScalar m_maxStress;       // Maximum principle first Piola stress
	};
	/* SolverState	*/
	struct SolverState
	{
		//if you add new variables, always initialize them!
		SolverState()
			: sdt(0),
			  isdt(0),
			  velmrg(0),
			  radmrg(0),
			  updmrg(0)
		{
		}
		btScalar sdt;     // dt*timescale
		btScalar isdt;    // 1/sdt
		btScalar velmrg;  // velocity margin
		btScalar radmrg;  // radial margin
		btScalar updmrg;  // Update margin
	};
	/// RayFromToCaster takes a ray from, ray to (instead of direction!)
	struct RayFromToCaster : btDbvt::ICollide
	{
		btVector3 m_rayFrom;
		btVector3 m_rayTo;
		btVector3 m_rayNormalizedDirection;
		btScalar m_mint;
		Face* m_face;
		int m_tests;
		RayFromToCaster(const btVector3& rayFrom, const btVector3& rayTo, btScalar mxt);
		void Process(const btDbvtNode* leaf);

		static /*inline*/ btScalar rayFromToTriangle(const btVector3& rayFrom,
													 const btVector3& rayTo,
													 const btVector3& rayNormalizedDirection,
													 const btVector3& a,
													 const btVector3& b,
													 const btVector3& c,
													 btScalar maxt = SIMD_INFINITY);
	};

	struct btVertexToTetraMapping
	{
		unsigned vertexToTetra;               // Index of first tetra index to which the collision mesh vertex is attached
		btVector4 baryCoordInTetra;           // Barycentric coordiante of collision shape vertex in that tetra
		btVector4 baryCoordNormalEndInTetra;  // Barycentric coordiante of collision shape vertex normal end in that tetra
	};

	//
	// Typedefs
	//

	typedef void (*psolver_t)(btSoftBody*, btScalar, btScalar);
	typedef void (*vsolver_t)(btSoftBody*, btScalar);
	typedef btAlignedObjectArray<Cluster*> tClusterArray;
	typedef btAlignedObjectArray<Note> tNoteArray;
	typedef btAlignedObjectArray<Node> tNodeArray;
	typedef btAlignedObjectArray<RenderNode> tRenderNodeArray;
	typedef btAlignedObjectArray<btDbvtNode*> tLeafArray;
	typedef btAlignedObjectArray<Link> tLinkArray;
	typedef btAlignedObjectArray<Face> tFaceArray;
	typedef btAlignedObjectArray<RenderFace> tRenderFaceArray;
	typedef btAlignedObjectArray<Tetra> tTetraArray;
	typedef btAlignedObjectArray<Anchor> tAnchorArray;
	typedef btAlignedObjectArray<RContact> tRContactArray;
	typedef btAlignedObjectArray<SContact> tSContactArray;
	typedef btAlignedObjectArray<Material*> tMaterialArray;
	typedef btAlignedObjectArray<Joint*> tJointArray;
	typedef btAlignedObjectArray<btSoftBody*> tSoftBodyArray;
	typedef btAlignedObjectArray<btAlignedObjectArray<btScalar>> tDenseMatrix;

	//
	// Fields
	//

	Config m_cfg;                      // Configuration
	SolverState m_sst;                 // Solver state
	Pose m_pose;                       // Pose
	void* m_tag;                       // User data
	btSoftBodyWorldInfo* m_worldInfo;  // World info
	tNoteArray m_notes;                // Notes
	tNodeArray m_nodes;                // Nodes
	tRenderNodeArray m_renderNodes;    // Render Nodes
	tLinkArray m_links;                // Links
	tFaceArray m_faces;                // Faces
	tRenderFaceArray m_renderFaces;    // Faces
	tTetraArray m_tetras;              // Tetras
	btAlignedObjectArray<TetraScratch> m_tetraScratches;
	btAlignedObjectArray<TetraScratch> m_tetraScratchesTn;
	tAnchorArray m_anchors;  // Anchors
	btAlignedObjectArray<DeformableNodeRigidAnchor> m_deformableAnchors;
	tRContactArray m_rcontacts;  // Rigid contacts
	btAlignedObjectArray<DeformableNodeRigidContact> m_nodeRigidContacts;
	btAlignedObjectArray<DeformableFaceNodeContact> m_faceNodeContacts;
	btAlignedObjectArray<DeformableNodeNodeContact> m_nodeNodeContacts;
	btAlignedObjectArray<DeformableFaceRigidContact> m_faceRigidContacts;
	btAlignedObjectArray<DeformableFaceNodeContact> m_faceNodeContactsCCD;
	tSContactArray m_scontacts;     // Soft contacts
	tJointArray m_joints;           // Joints
	tMaterialArray m_materials;     // Materials
	btScalar m_timeacc;             // Time accumulator
	btVector3 m_bounds[2];          // Spatial bounds
	bool m_bUpdateRtCst;            // Update runtime constants
	btDbvt m_ndbvt;                 // Nodes tree
	btDbvt m_fdbvt;                 // Faces tree
	btDbvntNode* m_fdbvnt;          // Faces tree with normals
	btDbvt m_cdbvt;                 // Clusters tree
	tClusterArray m_clusters;       // Clusters
	btScalar m_dampingCoefficient;  // Damping Coefficient
	btScalar m_sleepingThreshold;
	btScalar m_maxSpeedSquared;
	btAlignedObjectArray<btVector3> m_quads;  // quadrature points for collision detection
	btScalar m_repulsionStiffness;
	btScalar m_gravityFactor;
	bool m_cacheBarycenter;
	btAlignedObjectArray<btVector3> m_X;  // initial positions

	btAlignedObjectArray<btVector4> m_renderNodesInterpolationWeights;
	btAlignedObjectArray<btAlignedObjectArray<const btSoftBody::Node*>> m_renderNodesParents;
	btAlignedObjectArray<btScalar> m_z;  // vertical distance used in extrapolation
	bool m_useSelfCollision;
	bool m_softSoftCollision;

	btAlignedObjectArray<bool> m_clusterConnectivity;  //cluster connectivity, for self-collision

	btVector3 m_windVelocity;

	btScalar m_restLengthScale;
	btScalar m_averagePrincipalStress;
	btScalar m_lastSafeApplyDepthThreshold;
	btScalar m_lastSafeApplyVelocityDamping;
	btScalar m_softVsSoftContactStiffness;

	bool m_reducedModel;  // Reduced deformable model flag

	//static SimpleOBJDrawer drawer;

	//
	// Api
	//

	/* ctor																	*/
	btSoftBody(btSoftBodyWorldInfo* worldInfo, int node_count, const btVector3* x, const btScalar* m, btCollisionShape* collisionShape = nullptr);

	/* ctor																	*/
	btSoftBody(btSoftBodyWorldInfo* worldInfo);

	void initDefaults(btCollisionShape* collisionShape);

	/* dtor																	*/
	virtual ~btSoftBody();
	/* Check for existing link												*/

	btAlignedObjectArray<int> m_userIndexMapping;

	btSoftBodyWorldInfo* getWorldInfo()
	{
		return m_worldInfo;
	}

	void setDampingCoefficient(btScalar damping_coeff)
	{
		m_dampingCoefficient = damping_coeff;
	}

	virtual void setCollisionShape(btCollisionShape* collisionShape)
	{
		m_collisionShape = collisionShape;
	}

	bool checkLink(int node0,
				   int node1) const;
	bool checkLink(const Node* node0,
				   const Node* node1) const;
	/* Check for existing face												*/
	bool checkFace(int node0,
				   int node1,
				   int node2) const;
	/* Append material														*/
	Material* appendMaterial();
	/* Append note															*/
	void appendNote(const char* text,
					const btVector3& o,
					const btVector4& c = btVector4(1, 0, 0, 0),
					Node* n0 = 0,
					Node* n1 = 0,
					Node* n2 = 0,
					Node* n3 = 0);
	void appendNote(const char* text,
					const btVector3& o,
					Node* feature);
	void appendNote(const char* text,
					const btVector3& o,
					Link* feature);
	void appendNote(const char* text,
					const btVector3& o,
					Face* feature);
	/* Append node															*/
	void appendNode(const btVector3& x, btScalar m);
	/* Append link															*/
	void appendLink(int model = -1, Material* mat = 0);
	void appendLink(int node0,
					int node1,
					Material* mat = 0,
					bool bcheckexist = false);
	void appendLink(Node* node0,
					Node* node1,
					Material* mat = 0,
					bool bcheckexist = false);
	/* Append face															*/
	void appendFace(int model = -1, Material* mat = 0);
	void appendFace(int node0,
					int node1,
					int node2,
					Material* mat = 0);
	void appendTetra(int model, Material* mat);
	//
	void appendTetra(int node0,
					 int node1,
					 int node2,
					 int node3,
					 Material* mat = 0);

	/* Append anchor														*/
	void appendDeformableAnchor(int node, btRigidBody* body, uint32_t userIndex = -1);
	void appendDeformableAnchor(int node, btMultiBodyLinkCollider* link, uint32_t userIndex = -1);
	void appendAnchor(int node,
					  btRigidBody* body, bool disableCollisionBetweenLinkedBodies = false, btScalar influence = 1);
	void appendAnchor(int node, btRigidBody* body, const btVector3& localPivot, bool disableCollisionBetweenLinkedBodies = false, btScalar influence = 1);
	void removeAnchor(int node);
	void removeDeformableAnchor(int node);
	int removeDeformableAnchorByUserIndex(int userIndex);
	std::vector<int> getDeformableAnchorByUserIndex(int userIndex) const;
	/* Append linear joint													*/
	void appendLinearJoint(const LJoint::Specs& specs, Cluster* body0, Body body1);
	void appendLinearJoint(const LJoint::Specs& specs, Body body = Body());
	void appendLinearJoint(const LJoint::Specs& specs, btSoftBody* body);
	/* Append linear joint													*/
	void appendAngularJoint(const AJoint::Specs& specs, Cluster* body0, Body body1);
	void appendAngularJoint(const AJoint::Specs& specs, Body body = Body());
	void appendAngularJoint(const AJoint::Specs& specs, btSoftBody* body);
	/* Add force (or gravity) to the entire body							*/
	void addForce(const btVector3& force);
	/* Add force (or gravity) to a node of the body							*/
	void addForce(const btVector3& force,
				  int node);
	/* Add aero force to a node of the body */
	void addAeroForceToNode(const btVector3& windVelocity, int nodeIndex);

	/* Add aero force to a face of the body */
	void addAeroForceToFace(const btVector3& windVelocity, int faceIndex);

	/* Add velocity to the entire body										*/
	void addVelocity(const btVector3& velocity);

	/* Set velocity for the entire body										*/
	void setVelocity(const btVector3& velocity);

	/* Add velocity to a node of the body									*/
	void addVelocity(const btVector3& velocity,
					 int node);
	/* Set mass																*/
	void setMass(int node,
				 btScalar mass);
	/* Get mass																*/
	btScalar getMass(int node) const;
	/* Get total mass														*/
	btScalar getTotalMass() const;
	/* Set total mass (weighted by previous masses)							*/
	void setTotalMass(btScalar mass,
					  bool fromfaces = false);
	/* Set total density													*/
	void setTotalDensity(btScalar density);
	/* Set volume mass (using tetrahedrons)									*/
	void setVolumeMass(btScalar mass);
	/* Set volume density (using tetrahedrons)								*/
	void setVolumeDensity(btScalar density);
	/* Get the linear velocity of the center of mass                        */
	btVector3 getLinearVelocity();
	/* Set the linear velocity of the center of mass                        */
	void setLinearVelocity(const btVector3& linVel);
	/* Set the angular velocity of the center of mass                       */
	void setAngularVelocity(const btVector3& angVel);
	/* Get best fit rigid transform                                         */
	btTransform getRigidTransform();
	/* Transform to given pose                                              */
	virtual void transformTo(const btTransform& trs);
	/* Transform															*/
	virtual void transform(const btTransform& trs);
	/* Translate															*/
	virtual void translate(const btVector3& trs);
	/* Rotate															*/
	virtual void rotate(const btQuaternion& rot);
	/* Scale																*/
	virtual void scale(const btVector3& scl);
	/* Get link resting lengths scale										*/
	btScalar getRestLengthScale();
	/* Scale resting length of all springs									*/
	void setRestLengthScale(btScalar restLength);
	/* Set current state as pose											*/
	void setPose(bool bvolume,
				 bool bframe);
	/* Set current link lengths as resting lengths							*/
	void resetLinkRestLengths();
	/* Return the volume													*/
	btScalar getVolume() const;
	/* Cluster count														*/
	btVector3 getCenterOfMass() const
	{
		btVector3 com(0, 0, 0);
		for (int i = 0; i < m_nodes.size(); i++)
		{
			com += (m_nodes[i].m_x * this->getMass(i));
		}
		com /= this->getTotalMass();
		return com;
	}
	int clusterCount() const;
	/* Cluster center of mass												*/
	static btVector3 clusterCom(const Cluster* cluster);
	btVector3 clusterCom(int cluster) const;
	/* Cluster velocity at rpos												*/
	static btVector3 clusterVelocity(const Cluster* cluster, const btVector3& rpos);
	/* Cluster impulse														*/
	static void clusterVImpulse(Cluster* cluster, const btVector3& rpos, const btVector3& impulse);
	static void clusterDImpulse(Cluster* cluster, const btVector3& rpos, const btVector3& impulse);
	static void clusterImpulse(Cluster* cluster, const btVector3& rpos, const Impulse& impulse);
	static void clusterVAImpulse(Cluster* cluster, const btVector3& impulse);
	static void clusterDAImpulse(Cluster* cluster, const btVector3& impulse);
	static void clusterAImpulse(Cluster* cluster, const Impulse& impulse);
	static void clusterDCImpulse(Cluster* cluster, const btVector3& impulse);
	/* Generate bending constraints based on distance in the adjency graph	*/
	int generateBendingConstraints(int distance,
								   Material* mat = 0);
	/* Randomize constraints to reduce solver bias							*/
	void randomizeConstraints();

	void updateState(const btAlignedObjectArray<btVector3>& qs, const btAlignedObjectArray<btVector3>& vs);

	/* Release clusters														*/
	void releaseCluster(int index);
	void releaseClusters();
	/* Generate clusters (K-mean)											*/
	///generateClusters with k=0 will create a convex cluster for each tetrahedron or triangle
	///otherwise an approximation will be used (better performance)
	int generateClusters(int k, int maxiterations = 8192);
	/* Refine																*/
	void refine(ImplicitFn* ifn, btScalar accurary, bool cut);
	/* CutLink																*/
	bool cutLink(int node0, int node1, btScalar position);
	bool cutLink(const Node* node0, const Node* node1, btScalar position);

	///Ray casting using rayFrom and rayTo in worldspace, (not direction!)
	bool rayTest(const btVector3& rayFrom,
				 const btVector3& rayTo,
				 sRayCast& results);
	bool rayFaceTest(const btVector3& rayFrom,
					 const btVector3& rayTo,
					 sRayCast& results);
	int rayFaceTest(const btVector3& rayFrom, const btVector3& rayTo,
					btScalar& mint, int& index) const;
	/* Solver presets														*/
	void setSolver(eSolverPresets::_ preset);
	/* predictMotion														*/
	void predictMotion(btScalar dt);
	/* solveConstraints														*/
	void solveConstraints();
	/* staticSolve															*/
	void staticSolve(int iterations);
	/* solveCommonConstraints												*/
	static void solveCommonConstraints(btSoftBody** bodies, int count, int iterations);
	/* solveClusters														*/
	static void solveClusters(const btAlignedObjectArray<btSoftBody*>& bodies);
	/* integrateMotion														*/
	void integrateMotion();
	/* defaultCollisionHandlers												*/
	void defaultCollisionHandler(const btCollisionObjectWrapper* pcoWrap);
	void defaultCollisionHandler(btSoftBody* psb);
	void skinSoftRigidCollisionHandler(const btCollisionObjectWrapper* pcoWrap, int part0, int index0, const btVector3& contactPointOnSoftCollisionMesh, btVector3 contactNormalOnSoftCollisionMesh, btScalar penetrationDepth, const bool penetrating, btScalar* contactPointImpulseMagnitude);
	void skinSoftSoftCollisionHandler(btSoftBody* otherSoft, int part0, int index0, int part1, int index1, const btVector3& contactPointOnSoftCollisionMesh, btVector3 contactNormalOnSoftCollisionMesh, btScalar distance, const bool penetrating, btScalar* contactPointImpulseMagnitude);
	std::vector<int> findNClosestNodesLinearComplexity(const btVector3& p, int N) const;
	int findClosestNodeByMapping(int part, int triIndex, const btVector3& p) const;
	void setSelfCollision(bool useSelfCollision);
	bool useSelfCollision();
	void updateDeactivation(btScalar timeStep);
	void setZeroVelocity();
	bool wantsSleeping();

	virtual btMatrix3x3 getImpulseFactor(int n_node)
	{
		btMatrix3x3 tmp;
		tmp.setIdentity();
		return tmp;
	}

	//
	// Functionality to deal with new accelerated solvers.
	//

	/**
	 * Set a wind velocity for interaction with the air.
	 */
	void setWindVelocity(const btVector3& velocity);

	/**
	 * Return the wind velocity for interaction with the air.
	 */
	const btVector3& getWindVelocity();

	//
	// Set the solver that handles this soft body
	// Should not be allowed to get out of sync with reality
	// Currently called internally on addition to the world
	void setSoftBodySolver(btSoftBodySolver* softBodySolver)
	{
		m_softBodySolver = softBodySolver;
	}

	//
	// Return the solver that handles this soft body
	//
	btSoftBodySolver* getSoftBodySolver()
	{
		return m_softBodySolver;
	}

	//
	// Return the solver that handles this soft body
	//
	btSoftBodySolver* getSoftBodySolver() const
	{
		return m_softBodySolver;
	}

	//
	// Cast
	//

	static const btSoftBody* upcast(const btCollisionObject* colObj)
	{
		if (colObj->getInternalType() == CO_SOFT_BODY)
			return (const btSoftBody*)colObj;
		return 0;
	}
	static btSoftBody* upcast(btCollisionObject* colObj)
	{
		if (colObj->getInternalType() == CO_SOFT_BODY)
			return (btSoftBody*)colObj;
		return 0;
	}

	//
	// ::btCollisionObject
	//

	virtual void getAabb(btVector3& aabbMin, btVector3& aabbMax) const
	{
		aabbMin = m_bounds[0];
		aabbMax = m_bounds[1];
	}
	//
	// Private
	//
	void pointersToIndices();
	void indicesToPointers(const int* map = 0);

	int rayTest(const btVector3& rayFrom, const btVector3& rayTo,
				btScalar& mint, eFeature::_& feature, int& index, bool bcountonly) const;
	void initializeFaceTree();
	void rebuildNodeTree();
	btVector3 evaluateCom() const;
	bool checkDeformableContact(const btCollisionObjectWrapper* colObjWrap, const btVector3& x, btScalar margin, btSoftBody::sCti& cti, bool predict = false) const;
	bool checkDeformableFaceContact(const btCollisionObjectWrapper* colObjWrap, Face& f, btVector3& contact_point, btVector3& bary, btScalar margin, btSoftBody::sCti& cti, bool predict = false) const;
	bool checkContact(const btCollisionObjectWrapper* colObjWrap, const btVector3& x, btScalar margin, btSoftBody::sCti& cti) const;
	void updateNormals();
	void updateBounds();
	void updatePose();
	void updateConstants();
	void updateLinkConstants();
	void updateArea(bool averageArea = true);
	void initializeClusters();
	void updateClusters();
	void cleanupClusters();
	void prepareClusters(int iterations);
	void solveClusters(btScalar sor);
	void applyClusters(bool drift);
	void dampClusters();
	void setSpringStiffness(btScalar k);
	void setGravityFactor(btScalar gravFactor);
	void setCacheBarycenter(bool cacheBarycenter);
	void initializeDmInverse();
	void updateDeformation();
	void advanceDeformation();
	void applyForces();
	void setMaxStress(btScalar maxStress);
	void interpolateRenderMesh();
	void setCollisionQuadrature(int N);
	static void PSolve_Anchors(btSoftBody* psb, btScalar kst, btScalar ti);
	static void PSolve_RContacts(btSoftBody* psb, btScalar kst, btScalar ti);
	static void PSolve_SContacts(btSoftBody* psb, btScalar, btScalar ti);
	static void PSolve_Links(btSoftBody* psb, btScalar kst, btScalar ti);
	static void VSolve_Links(btSoftBody* psb, btScalar kst);
	static psolver_t getSolver(ePSolver::_ solver);
	static vsolver_t getSolver(eVSolver::_ solver);
	void geometricCollisionHandler(btSoftBody* psb);
#define SAFE_EPSILON SIMD_EPSILON * 100.0
	void updateNode(btDbvtNode* node, bool use_velocity, bool margin)
	{
		if (node->isleaf())
		{
			btSoftBody::Node* n = (btSoftBody::Node*)(node->data);
			ATTRIBUTE_ALIGNED16(btDbvtVolume)
			vol;
			btScalar pad = margin ? m_sst.radmrg : SAFE_EPSILON;  // use user defined margin or margin for floating point precision
			if (use_velocity)
			{
				btVector3 points[2] = {n->m_x, n->m_x + m_sst.sdt * n->m_v};
				vol = btDbvtVolume::FromPoints(points, 2);
				vol.Expand(btVector3(pad, pad, pad));
			}
			else
			{
				vol = btDbvtVolume::FromCR(n->m_x, pad);
			}
			node->volume = vol;
			return;
		}
		else
		{
			updateNode(node->childs[0], use_velocity, margin);
			updateNode(node->childs[1], use_velocity, margin);
			ATTRIBUTE_ALIGNED16(btDbvtVolume)
			vol;
			Merge(node->childs[0]->volume, node->childs[1]->volume, vol);
			node->volume = vol;
		}
	}

	void updateNodeTree(bool use_velocity, bool margin)
	{
		if (m_ndbvt.m_root)
			updateNode(m_ndbvt.m_root, use_velocity, margin);
	}

	template <class DBVTNODE>  // btDbvtNode or btDbvntNode
	void updateFace(DBVTNODE* node, bool use_velocity, bool margin)
	{
		if (node->isleaf())
		{
			btSoftBody::Face* f = (btSoftBody::Face*)(node->data);
			btScalar pad = margin ? m_sst.radmrg : SAFE_EPSILON;  // use user defined margin or margin for floating point precision
			ATTRIBUTE_ALIGNED16(btDbvtVolume)
			vol;
			if (use_velocity)
			{
				btVector3 points[6] = {f->m_n[0]->m_x, f->m_n[0]->m_x + m_sst.sdt * f->m_n[0]->m_v,
									   f->m_n[1]->m_x, f->m_n[1]->m_x + m_sst.sdt * f->m_n[1]->m_v,
									   f->m_n[2]->m_x, f->m_n[2]->m_x + m_sst.sdt * f->m_n[2]->m_v};
				vol = btDbvtVolume::FromPoints(points, 6);
			}
			else
			{
				btVector3 points[3] = {f->m_n[0]->m_x,
									   f->m_n[1]->m_x,
									   f->m_n[2]->m_x};
				vol = btDbvtVolume::FromPoints(points, 3);
			}
			vol.Expand(btVector3(pad, pad, pad));
			node->volume = vol;
			return;
		}
		else
		{
			updateFace(node->childs[0], use_velocity, margin);
			updateFace(node->childs[1], use_velocity, margin);
			ATTRIBUTE_ALIGNED16(btDbvtVolume)
			vol;
			Merge(node->childs[0]->volume, node->childs[1]->volume, vol);
			node->volume = vol;
		}
	}
	void updateFaceTree(bool use_velocity, bool margin)
	{
		if (m_fdbvt.m_root)
			updateFace(m_fdbvt.m_root, use_velocity, margin);
		if (m_fdbvnt)
			updateFace(m_fdbvnt, use_velocity, margin);
	}

	template <typename T>
	static inline T BaryEval(const T& a,
							 const T& b,
							 const T& c,
							 const btVector3& coord)
	{
		return (a * coord.x() + b * coord.y() + c * coord.z());
	}

	void applyRepulsionForce(btScalar timeStep, bool applySpringForce);
	virtual int calculateSerializeBufferSize() const;

	///fills the dataBuffer and returns the struct name (and 0 on failure)
	virtual const char* serialize(void* dataBuffer, class btSerializer* serializer) const;

	virtual const std::vector<btVertexToTetraMapping>* getCollisionShapeVertexToSimTetra() const
	{
		return nullptr;
	}

	const tAnchorArray& getAnchors() const
	{
		return m_anchors;
	}

	const btAlignedObjectArray<DeformableNodeRigidAnchor>& getDeformableAnchors() const
	{
		return m_deformableAnchors;
	}

	virtual void resetColObjPtrsInAnchors(btCollisionObject* co)
	{
		for (auto i = 0; i < m_anchors.size(); ++i)
			if (m_anchors[i].m_body == co)
				m_anchors[i].m_body = nullptr;
		for (auto i = 0; i < m_deformableAnchors.size(); ++i)
		{
			if (m_deformableAnchors[i].m_body == co)
				m_deformableAnchors[i].m_body = nullptr;
			if (m_deformableAnchors[i].m_cti.m_colObj == co)
				m_deformableAnchors[i].m_cti.m_colObj = nullptr;
		}
	}

	void lastSafeBorderGrow(int growth, std::map<btSoftBody::Node*, StuckTetraIndicesMapped>& nodesInCollision);
	virtual void updateLastSafeWorldTransform() override;
	virtual void applyLastSafeWorldTransform(const std::map<int, StuckTetraIndicesMapped>* partial) override;
};

#endif  //_BT_SOFT_BODY_H
/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/


#include "btSoftBodyConcaveCollisionAlgorithm.h"
#include "BulletCollision/CollisionDispatch/btCollisionObject.h"
#include "BulletCollision/CollisionShapes/btMultiSphereShape.h"
#include "BulletCollision/BroadphaseCollision/btBroadphaseProxy.h"
#include "BulletCollision/CollisionShapes/btConcaveShape.h"
#include "BulletCollision/CollisionDispatch/btManifoldResult.h"
#include "BulletCollision/NarrowPhaseCollision/btRaycastCallback.h"
#include "BulletCollision/CollisionShapes/btTriangleShape.h"
#include "BulletCollision/CollisionShapes/btSphereShape.h"
#include "BulletCollision/CollisionShapes/btTetrahedronShape.h"
#include "BulletCollision/CollisionShapes/btConvexHullShape.h"
#include "BulletCollision/CollisionDispatch/btCollisionObjectWrapper.h"

#include "LinearMath/btIDebugDraw.h"
#include "BulletCollision/NarrowPhaseCollision/btSubSimplexConvexCast.h"
#include "BulletSoftBody/btSoftBody.h"

#define BT_SOFTBODY_TRIANGLE_EXTRUSION btScalar(0.06)  //make this configurable

btSoftBodyConcaveCollisionAlgorithm::btSoftBodyConcaveCollisionAlgorithm(const btCollisionAlgorithmConstructionInfo& ci, const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, bool isSwapped)
	: btCollisionAlgorithm(ci),
	  m_isSwapped(isSwapped),
	  m_btSoftBodyTriangleCallback(ci.m_dispatcher1, body0Wrap, body1Wrap, isSwapped)
{
}

btSoftBodyConcaveCollisionAlgorithm::~btSoftBodyConcaveCollisionAlgorithm()
{
}

btSoftBodyTriangleCallback::btSoftBodyTriangleCallback(btDispatcher* dispatcher, const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, bool isSwapped) : m_dispatcher(dispatcher),
																																														 m_dispatchInfoPtr(0)
{
	m_softBody = (isSwapped ? (btSoftBody*)body1Wrap->getCollisionObject() : (btSoftBody*)body0Wrap->getCollisionObject());
	m_triBody = isSwapped ? body0Wrap->getCollisionObject() : body1Wrap->getCollisionObject();

	//
	// create the manifold from the dispatcher 'manifold pool'
	//
	//	  m_manifoldPtr = m_dispatcher->getNewManifold(m_convexBody,m_triBody);

	clearCache();
}

btSoftBodyTriangleCallback::~btSoftBodyTriangleCallback()
{
	clearCache();
	//	m_dispatcher->releaseManifold( m_manifoldPtr );
}

void btSoftBodyTriangleCallback::clearCache()
{
	for (int i = 0; i < m_shapeCache.size(); i++)
	{
		btTriIndex* tmp = m_shapeCache.getAtIndex(i);
		btAssert(tmp);
		btAssert(tmp->m_childShape);
		m_softBody->getWorldInfo()->m_sparsesdf.RemoveReferences(tmp->m_childShape);  //necessary?
		delete tmp->m_childShape;
	}
	m_shapeCache.clear();
}

void btSoftBodyTriangleCallback::processTriangle(btVector3* triangle, int partId, int triangleIndex)
{
	//just for debugging purposes
	//printf("triangle %d",m_triangleCount++);

	btCollisionAlgorithmConstructionInfo ci;
	ci.m_dispatcher1 = m_dispatcher;

	///debug drawing of the overlapping triangles
	if (m_dispatchInfoPtr && m_dispatchInfoPtr->m_debugDraw && (m_dispatchInfoPtr->m_debugDraw->getDebugMode() & btIDebugDraw::DBG_DrawWireframe))
	{
		btVector3 color(1, 1, 0);
		const btTransform& tr = m_triBody->getWorldTransform();
		m_dispatchInfoPtr->m_debugDraw->drawLine(tr(triangle[0]), tr(triangle[1]), color);
		m_dispatchInfoPtr->m_debugDraw->drawLine(tr(triangle[1]), tr(triangle[2]), color);
		m_dispatchInfoPtr->m_debugDraw->drawLine(tr(triangle[2]), tr(triangle[0]), color);
	}

	btTriIndex triIndex(partId, triangleIndex, 0);
	btHashKey<btTriIndex> triKey(triIndex.getUid());

	btTriIndex* shapeIndex = m_shapeCache[triKey];
	if (shapeIndex)
	{
		btCollisionShape* tm = shapeIndex->m_childShape;
		btAssert(tm);

		//copy over user pointers to temporary shape
		tm->setUserPointer(m_triBody->getCollisionShape()->getUserPointer());

		btCollisionObjectWrapper softBody(0, m_softBody->getCollisionShape(), m_softBody, m_softBody->getWorldTransform(), -1, -1);
		//btCollisionObjectWrapper triBody(0,tm, ob, btTransform::getIdentity());//ob->getWorldTransform());//??
		btCollisionObjectWrapper triBody(0, tm, m_triBody, m_triBody->getWorldTransform(), partId, triangleIndex);
		ebtDispatcherQueryType algoType = m_resultOut->m_closestPointDistanceThreshold > 0 ? BT_CLOSEST_POINT_ALGORITHMS : BT_CONTACT_POINT_ALGORITHMS;
		btCollisionAlgorithm* colAlgo = ci.m_dispatcher1->findAlgorithm(&softBody, &triBody, 0, algoType);  //m_manifoldPtr);

		colAlgo->processCollision(&softBody, &triBody, *m_dispatchInfoPtr, m_resultOut);
		colAlgo->~btCollisionAlgorithm();
		ci.m_dispatcher1->freeCollisionAlgorithm(colAlgo);

		return;
	}

	//aabb filter is already applied!

	//btCollisionObject* colObj = static_cast<btCollisionObject*>(m_convexProxy->m_clientObject);

	//	if (m_softBody->getCollisionShape()->getShapeType()==
	{
		//		btVector3 other;
		btVector3 normal = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
		normal.normalize();
		normal *= BT_SOFTBODY_TRIANGLE_EXTRUSION;
		//		other=(triangle[0]+triangle[1]+triangle[2])*0.333333f;
		//		other+=normal*22.f;
		btVector3 pts[6] = {triangle[0] + normal,
							triangle[1] + normal,
							triangle[2] + normal,
							triangle[0] - normal,
							triangle[1] - normal,
							triangle[2] - normal};

		btConvexHullShape* tm = new btConvexHullShape(&pts[0].getX(), 6);

		//		btBU_Simplex1to4 tm(triangle[0],triangle[1],triangle[2],other);

		//btTriangleShape tm(triangle[0],triangle[1],triangle[2]);
		//	tm.setMargin(m_collisionMarginTriangle);

		//copy over user pointers to temporary shape
		tm->setUserPointer(m_triBody->getCollisionShape()->getUserPointer());

		btCollisionObjectWrapper softBody(0, m_softBody->getCollisionShape(), m_softBody, m_softBody->getWorldTransform(), -1, -1);
		btCollisionObjectWrapper triBody(0, tm, m_triBody, m_triBody->getWorldTransform(), partId, triangleIndex);  //btTransform::getIdentity());//??

		ebtDispatcherQueryType algoType = m_resultOut->m_closestPointDistanceThreshold > 0 ? BT_CLOSEST_POINT_ALGORITHMS : BT_CONTACT_POINT_ALGORITHMS;
		btCollisionAlgorithm* colAlgo = ci.m_dispatcher1->findAlgorithm(&softBody, &triBody, 0, algoType);  //m_manifoldPtr);

		colAlgo->processCollision(&softBody, &triBody, *m_dispatchInfoPtr, m_resultOut);
		colAlgo->~btCollisionAlgorithm();
		ci.m_dispatcher1->freeCollisionAlgorithm(colAlgo);

		triIndex.m_childShape = tm;
		m_shapeCache.insert(triKey, triIndex);
	}
}

void btSoftBodyTriangleCallback::setTimeStepAndCounters(btScalar collisionMarginTriangle, const btCollisionObjectWrapper* triBodyWrap, const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
{
	m_dispatchInfoPtr = &dispatchInfo;
	m_collisionMarginTriangle = collisionMarginTriangle + btScalar(BT_SOFTBODY_TRIANGLE_EXTRUSION);
	m_resultOut = resultOut;

	btVector3 aabbWorldSpaceMin, aabbWorldSpaceMax;
	m_softBody->getAabb(aabbWorldSpaceMin, aabbWorldSpaceMax);
	btVector3 halfExtents = (aabbWorldSpaceMax - aabbWorldSpaceMin) * btScalar(0.5);
	btVector3 softBodyCenter = (aabbWorldSpaceMax + aabbWorldSpaceMin) * btScalar(0.5);

	btTransform softTransform;
	softTransform.setIdentity();
	softTransform.setOrigin(softBodyCenter);

	btTransform convexInTriangleSpace;
	convexInTriangleSpace = triBodyWrap->getWorldTransform().inverse() * softTransform;
	btTransformAabb(halfExtents, m_collisionMarginTriangle, convexInTriangleSpace, m_aabbMin, m_aabbMax);
}

void btSoftBodyConcaveCollisionAlgorithm::clearCache()
{
	m_btSoftBodyTriangleCallback.clearCache();
}

void btSoftBodyConcaveCollisionAlgorithm::processCollision(const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
{
	//btCollisionObject* convexBody = m_isSwapped ? body1 : body0;
	const btCollisionObjectWrapper* triBody = m_isSwapped ? body0Wrap : body1Wrap;

	if (triBody->getCollisionShape()->isConcave())
	{
		const btConcaveShape* concaveShape = static_cast<const btConcaveShape*>(triBody->getCollisionShape());

		//	if (convexBody->getCollisionShape()->isConvex())
		{
			btScalar collisionMarginTriangle = concaveShape->getMargin();

			//			resultOut->setPersistentManifold(m_btSoftBodyTriangleCallback.m_manifoldPtr);
			m_btSoftBodyTriangleCallback.setTimeStepAndCounters(collisionMarginTriangle, triBody, dispatchInfo, resultOut);

			concaveShape->processAllTriangles(&m_btSoftBodyTriangleCallback, m_btSoftBodyTriangleCallback.getAabbMin(), m_btSoftBodyTriangleCallback.getAabbMax());

			//	resultOut->refreshContactPoints();
		}
	}
}

btScalar btSoftBodyConcaveCollisionAlgorithm::calculateTimeOfImpact(btCollisionObject* body0, btCollisionObject* body1, const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
{
	(void)resultOut;
	(void)dispatchInfo;
	btCollisionObject* convexbody = m_isSwapped ? body1 : body0;
	btCollisionObject* triBody = m_isSwapped ? body0 : body1;

	//quick approximation using raycast, todo: hook up to the continuous collision detection (one of the btConvexCast)

	//only perform CCD above a certain threshold, this prevents blocking on the long run
	//because object in a blocked ccd state (hitfraction<1) get their linear velocity halved each frame...
	btScalar squareMot0 = (convexbody->getInterpolationWorldTransform().getOrigin() - convexbody->getWorldTransform().getOrigin()).length2();
	if (squareMot0 < convexbody->getCcdSquareMotionThreshold())
	{
		return btScalar(1.);
	}

	//const btVector3& from = convexbody->m_worldTransform.getOrigin();
	//btVector3 to = convexbody->m_interpolationWorldTransform.getOrigin();
	//todo: only do if the motion exceeds the 'radius'

	btTransform triInv = triBody->getWorldTransform().inverse();
	btTransform convexFromLocal = triInv * convexbody->getWorldTransform();
	btTransform convexToLocal = triInv * convexbody->getInterpolationWorldTransform();

	struct LocalTriangleSphereCastCallback : public btTriangleCallback
	{
		btTransform m_ccdSphereFromTrans;
		btTransform m_ccdSphereToTrans;
		btTransform m_meshTransform;

		btScalar m_ccdSphereRadius;
		btScalar m_hitFraction;

		LocalTriangleSphereCastCallback(const btTransform& from, const btTransform& to, btScalar ccdSphereRadius, btScalar hitFraction)
			: m_ccdSphereFromTrans(from),
			  m_ccdSphereToTrans(to),
			  m_ccdSphereRadius(ccdSphereRadius),
			  m_hitFraction(hitFraction)
		{
		}

		virtual void processTriangle(btVector3* triangle, int partId, int triangleIndex)
		{
			(void)partId;
			(void)triangleIndex;
			//do a swept sphere for now
			btTransform ident;
			ident.setIdentity();
			btConvexCast::CastResult castResult;
			castResult.m_fraction = m_hitFraction;
			btSphereShape pointShape(m_ccdSphereRadius);
			btTriangleShape triShape(triangle[0], triangle[1], triangle[2]);
			btVoronoiSimplexSolver simplexSolver;
			btSubsimplexConvexCast convexCaster(&pointShape, &triShape, &simplexSolver);
			//GjkConvexCast	convexCaster(&pointShape,convexShape,&simplexSolver);
			//ContinuousConvexCollision convexCaster(&pointShape,convexShape,&simplexSolver,0);
			//local space?

			if (convexCaster.calcTimeOfImpact(m_ccdSphereFromTrans, m_ccdSphereToTrans,
											  ident, ident, castResult))
			{
				if (m_hitFraction > castResult.m_fraction)
					m_hitFraction = castResult.m_fraction;
			}
		}
	};

	if (triBody->getCollisionShape()->isConcave())
	{
		btVector3 rayAabbMin = convexFromLocal.getOrigin();
		rayAabbMin.setMin(convexToLocal.getOrigin());
		btVector3 rayAabbMax = convexFromLocal.getOrigin();
		rayAabbMax.setMax(convexToLocal.getOrigin());
		btScalar ccdRadius0 = convexbody->getCcdSweptSphereRadius();
		rayAabbMin -= btVector3(ccdRadius0, ccdRadius0, ccdRadius0);
		rayAabbMax += btVector3(ccdRadius0, ccdRadius0, ccdRadius0);

		btScalar curHitFraction = btScalar(1.);  //is this available?
		LocalTriangleSphereCastCallback raycastCallback(convexFromLocal, convexToLocal,
														convexbody->getCcdSweptSphereRadius(), curHitFraction);

		raycastCallback.m_hitFraction = convexbody->getHitFraction();

		btCollisionObject* concavebody = triBody;

		btConcaveShape* triangleMesh = (btConcaveShape*)concavebody->getCollisionShape();

		if (triangleMesh)
		{
			triangleMesh->processAllTriangles(&raycastCallback, rayAabbMin, rayAabbMax);
		}

		if (raycastCallback.m_hitFraction < convexbody->getHitFraction())
		{
			convexbody->setHitFraction(raycastCallback.m_hitFraction);
			return raycastCallback.m_hitFraction;
		}
	}

	return btScalar(1.);
}
/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/


#ifndef BT_SOFT_BODY_CONCAVE_COLLISION_ALGORITHM_H
#define BT_SOFT_BODY_CONCAVE_COLLISION_ALGORITHM_H

#include "BulletCollision/BroadphaseCollision/btCollisionAlgorithm.h"
#include "BulletCollision/BroadphaseCollision/btDispatcher.h"
#include "BulletCollision/BroadphaseCollision/btBroadphaseInterface.h"
#include "BulletCollision/CollisionShapes/btTriangleCallback.h"
#include "BulletCollision/NarrowPhaseCollision/btPersistentManifold.h"
class btDispatcher;
#include "BulletCollision/BroadphaseCollision/btBroadphaseProxy.h"
#include "BulletCollision/CollisionDispatch/btCollisionCreateFunc.h"
class btSoftBody;
class btCollisionShape;

#include "LinearMath/btHashMap.h"

#include "BulletCollision/BroadphaseCollision/btQuantizedBvh.h"  //for definition of MAX_NUM_PARTS_IN_BITS

struct btTriIndex
{
	int m_PartIdTriangleIndex;
	class btCollisionShape* m_childShape;

	btTriIndex(int partId, int triangleIndex, btCollisionShape* shape)
	{
		m_PartIdTriangleIndex = (partId << (31 - MAX_NUM_PARTS_IN_BITS)) | triangleIndex;
		m_childShape = shape;
	}

	int getTriangleIndex() const
	{
		// Get only the lower bits where the triangle index is stored
		unsigned int x = 0;
		unsigned int y = (~(x & 0)) << (31 - MAX_NUM_PARTS_IN_BITS);
		return (m_PartIdTriangleIndex & ~(y));
	}
	int getPartId() const
	{
		// Get only the highest bits where the part index is stored
		return (m_PartIdTriangleIndex >> (31 - MAX_NUM_PARTS_IN_BITS));
	}
	int getUid() const
	{
		return m_PartIdTriangleIndex;
	}
};

///For each triangle in the concave mesh that overlaps with the AABB of a soft body (m_softBody), processTriangle is called.
class btSoftBodyTriangleCallback : public btTriangleCallback
{
	btSoftBody* m_softBody;
	const btCollisionObject* m_triBody;

	btVector3 m_aabbMin;
	btVector3 m_aabbMax;

	btManifoldResult* m_resultOut;

	btDispatcher* m_dispatcher;
	const btDispatcherInfo* m_dispatchInfoPtr;
	btScalar m_collisionMarginTriangle;

	btHashMap<btHashKey<btTriIndex>, btTriIndex> m_shapeCache;

public:
	int m_triangleCount;

	//	btPersistentManifold*	m_manifoldPtr;

	btSoftBodyTriangleCallback(btDispatcher* dispatcher, const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, bool isSwapped);

	void setTimeStepAndCounters(btScalar collisionMarginTriangle, const btCollisionObjectWrapper* triObjWrap, const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut);

	virtual ~btSoftBodyTriangleCallback();

	virtual void processTriangle(btVector3* triangle, int partId, int triangleIndex);

	void clearCache();

	SIMD_FORCE_INLINE const btVector3& getAabbMin() const
	{
		return m_aabbMin;
	}
	SIMD_FORCE_INLINE const btVector3& getAabbMax() const
	{
		return m_aabbMax;
	}
};

/// btSoftBodyConcaveCollisionAlgorithm  supports collision between soft body shapes and (concave) trianges meshes.
class btSoftBodyConcaveCollisionAlgorithm : public btCollisionAlgorithm
{
	bool m_isSwapped;

	btSoftBodyTriangleCallback m_btSoftBodyTriangleCallback;

public:
	btSoftBodyConcaveCollisionAlgorithm(const btCollisionAlgorithmConstructionInfo& ci, const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, bool isSwapped);

	virtual ~btSoftBodyConcaveCollisionAlgorithm();

	virtual void processCollision(const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut);

	btScalar calculateTimeOfImpact(btCollisionObject* body0, btCollisionObject* body1, const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut);

	virtual void getAllContactManifolds(btManifoldArray& manifoldArray)
	{
		//we don't add any manifolds
	}

	void clearCache();

	struct CreateFunc : public btCollisionAlgorithmCreateFunc
	{
		virtual btCollisionAlgorithm* CreateCollisionAlgorithm(btCollisionAlgorithmConstructionInfo& ci, const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap)
		{
			void* mem = ci.m_dispatcher1->allocateCollisionAlgorithm(sizeof(btSoftBodyConcaveCollisionAlgorithm));
			return new (mem) btSoftBodyConcaveCollisionAlgorithm(ci, body0Wrap, body1Wrap, false);
		}
	};

	struct SwappedCreateFunc : public btCollisionAlgorithmCreateFunc
	{
		virtual btCollisionAlgorithm* CreateCollisionAlgorithm(btCollisionAlgorithmConstructionInfo& ci, const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap)
		{
			void* mem = ci.m_dispatcher1->allocateCollisionAlgorithm(sizeof(btSoftBodyConcaveCollisionAlgorithm));
			return new (mem) btSoftBodyConcaveCollisionAlgorithm(ci, body0Wrap, body1Wrap, true);
		}
	};
};

#endif  //BT_SOFT_BODY_CONCAVE_COLLISION_ALGORITHM_H
/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/


#ifndef BT_SOFTBODY_FLOAT_DATA
#define BT_SOFTBODY_FLOAT_DATA

#include "BulletCollision/CollisionDispatch/btCollisionObject.h"
#include "BulletDynamics/Dynamics/btRigidBody.h"

struct SoftBodyMaterialData
{
	float m_linearStiffness;
	float m_angularStiffness;
	float m_volumeStiffness;
	int m_flags;
};

struct SoftBodyNodeData
{
	SoftBodyMaterialData *m_material;
	btVector3FloatData m_position;
	btVector3FloatData m_previousPosition;
	btVector3FloatData m_velocity;
	btVector3FloatData m_accumulatedForce;
	btVector3FloatData m_normal;
	float m_inverseMass;
	float m_area;
	int m_attach;
	int m_pad;
};

struct SoftBodyLinkData
{
	SoftBodyMaterialData *m_material;
	int m_nodeIndices[2];  // Node pointers
	float m_restLength;    // Rest length
	int m_bbending;        // Bending link
};

struct SoftBodyFaceData
{
	btVector3FloatData m_normal;  // Normal
	SoftBodyMaterialData *m_material;
	int m_nodeIndices[3];  // Node pointers
	float m_restArea;      // Rest area
};

struct SoftBodyTetraData
{
	btVector3FloatData m_c0[4];  // gradients
	SoftBodyMaterialData *m_material;
	int m_nodeIndices[4];  // Node pointers
	float m_restVolume;    // Rest volume
	float m_c1;            // (4*kVST)/(im0+im1+im2+im3)
	float m_c2;            // m_c1/sum(|g0..3|^2)
	int m_pad;
};

struct SoftRigidAnchorData
{
	btMatrix3x3FloatData m_c0;        // Impulse matrix
	btVector3FloatData m_c1;          // Relative anchor
	btVector3FloatData m_localFrame;  // Anchor position in body space
	btRigidBodyData *m_rigidBody;
	int m_nodeIndex;  // Node pointer
	float m_c2;       // ima*dt
};

struct SoftBodyConfigData
{
	int m_aeroModel;                         // Aerodynamic model (default: V_Point)
	float m_baumgarte;                       // Velocities correction factor (Baumgarte)
	float m_damping;                         // Damping coefficient [0,1]
	float m_drag;                            // Drag coefficient [0,+inf]
	float m_lift;                            // Lift coefficient [0,+inf]
	float m_pressure;                        // Pressure coefficient [-inf,+inf]
	float m_volume;                          // Volume conversation coefficient [0,+inf]
	float m_dynamicFriction;                 // Dynamic friction coefficient [0,1]
	float m_poseMatch;                       // Pose matching coefficient [0,1]
	float m_rigidContactHardness;            // Rigid contacts hardness [0,1]
	float m_kineticContactHardness;          // Kinetic contacts hardness [0,1]
	float m_softContactHardness;             // Soft contacts hardness [0,1]
	float m_anchorHardness;                  // Anchors hardness [0,1]
	float m_softRigidClusterHardness;        // Soft vs rigid hardness [0,1] (cluster only)
	float m_softKineticClusterHardness;      // Soft vs kinetic hardness [0,1] (cluster only)
	float m_softSoftClusterHardness;         // Soft vs soft hardness [0,1] (cluster only)
	float m_softRigidClusterImpulseSplit;    // Soft vs rigid impulse split [0,1] (cluster only)
	float m_softKineticClusterImpulseSplit;  // Soft vs rigid impulse split [0,1] (cluster only)
	float m_softSoftClusterImpulseSplit;     // Soft vs rigid impulse split [0,1] (cluster only)
	float m_maxVolume;                       // Maximum volume ratio for pose
	float m_timeScale;                       // Time scale
	int m_velocityIterations;                // Velocities solver iterations
	int m_positionIterations;                // Positions solver iterations
	int m_driftIterations;                   // Drift solver iterations
	int m_clusterIterations;                 // Cluster solver iterations
	int m_collisionFlags;                    // Collisions flags
};

struct SoftBodyPoseData
{
	btMatrix3x3FloatData m_rot;    // Rotation
	btMatrix3x3FloatData m_scale;  // Scale
	btMatrix3x3FloatData m_aqq;    // Base scaling
	btVector3FloatData m_com;      // COM

	btVector3FloatData *m_positions;  // Reference positions
	float *m_weights;                 // Weights
	int m_numPositions;
	int m_numWeigts;

	int m_bvolume;       // Is valid
	int m_bframe;        // Is frame
	float m_restVolume;  // Rest volume
	int m_pad;
};

struct SoftBodyClusterData
{
	btTransformFloatData m_framexform;
	btMatrix3x3FloatData m_locii;
	btMatrix3x3FloatData m_invwi;
	btVector3FloatData m_com;
	btVector3FloatData m_vimpulses[2];
	btVector3FloatData m_dimpulses[2];
	btVector3FloatData m_lv;
	btVector3FloatData m_av;

	btVector3FloatData *m_framerefs;
	int *m_nodeIndices;
	float *m_masses;

	int m_numFrameRefs;
	int m_numNodes;
	int m_numMasses;

	float m_idmass;
	float m_imass;
	int m_nvimpulses;
	int m_ndimpulses;
	float m_ndamping;
	float m_ldamping;
	float m_adamping;
	float m_matching;
	float m_maxSelfCollisionImpulse;
	float m_selfCollisionImpulseFactor;
	int m_containsAnchor;
	int m_collide;
	int m_clusterIndex;
};

enum btSoftJointBodyType
{
	BT_JOINT_SOFT_BODY_CLUSTER = 1,
	BT_JOINT_RIGID_BODY,
	BT_JOINT_COLLISION_OBJECT
};

struct btSoftBodyJointData
{
	void *m_bodyA;
	void *m_bodyB;
	btVector3FloatData m_refs[2];
	float m_cfm;
	float m_erp;
	float m_split;
	int m_delete;
	btVector3FloatData m_relPosition[2];  //linear
	int m_bodyAtype;
	int m_bodyBtype;
	int m_jointType;
	int m_pad;
};

///do not change those serialization structures, it requires an updated sBulletDNAstr/sBulletDNAstr64
struct btSoftBodyFloatData
{
	btCollisionObjectFloatData m_collisionObjectData;

	SoftBodyPoseData *m_pose;
	SoftBodyMaterialData **m_materials;
	SoftBodyNodeData *m_nodes;
	SoftBodyLinkData *m_links;
	SoftBodyFaceData *m_faces;
	SoftBodyTetraData *m_tetrahedra;
	SoftRigidAnchorData *m_anchors;
	SoftBodyClusterData *m_clusters;
	btSoftBodyJointData *m_joints;

	int m_numMaterials;
	int m_numNodes;
	int m_numLinks;
	int m_numFaces;
	int m_numTetrahedra;
	int m_numAnchors;
	int m_numClusters;
	int m_numJoints;
	SoftBodyConfigData m_config;
};

#endif  //BT_SOFTBODY_FLOAT_DATA
/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/

///btSoftBodyHelpers.cpp by Nathanael Presson

#include "btSoftBodyInternals.h"
#include <stdio.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string.h>
#include <algorithm>
#include "btSoftBodyHelpers.h"
#include "LinearMath/btConvexHull.h"
#include "LinearMath/btConvexHullComputer.h"

#include <map>
#include <vector>

static void drawVertex(btIDebugDraw* idraw,
					   const btVector3& x, btScalar s, const btVector3& c)
{
	idraw->drawLine(x - btVector3(s, 0, 0), x + btVector3(s, 0, 0), c);
	idraw->drawLine(x - btVector3(0, s, 0), x + btVector3(0, s, 0), c);
	idraw->drawLine(x - btVector3(0, 0, s), x + btVector3(0, 0, s), c);
}

//
static void drawBox(btIDebugDraw* idraw,
					const btVector3& mins,
					const btVector3& maxs,
					const btVector3& color)
{
	const btVector3 c[] = {btVector3(mins.x(), mins.y(), mins.z()),
						   btVector3(maxs.x(), mins.y(), mins.z()),
						   btVector3(maxs.x(), maxs.y(), mins.z()),
						   btVector3(mins.x(), maxs.y(), mins.z()),
						   btVector3(mins.x(), mins.y(), maxs.z()),
						   btVector3(maxs.x(), mins.y(), maxs.z()),
						   btVector3(maxs.x(), maxs.y(), maxs.z()),
						   btVector3(mins.x(), maxs.y(), maxs.z())};
	idraw->drawLine(c[0], c[1], color);
	idraw->drawLine(c[1], c[2], color);
	idraw->drawLine(c[2], c[3], color);
	idraw->drawLine(c[3], c[0], color);
	idraw->drawLine(c[4], c[5], color);
	idraw->drawLine(c[5], c[6], color);
	idraw->drawLine(c[6], c[7], color);
	idraw->drawLine(c[7], c[4], color);
	idraw->drawLine(c[0], c[4], color);
	idraw->drawLine(c[1], c[5], color);
	idraw->drawLine(c[2], c[6], color);
	idraw->drawLine(c[3], c[7], color);
}

//
static void drawTree(btIDebugDraw* idraw,
					 const btDbvtNode* node,
					 int depth,
					 const btVector3& ncolor,
					 const btVector3& lcolor,
					 int mindepth,
					 int maxdepth)
{
	if (node)
	{
		if (node->isinternal() && ((depth < maxdepth) || (maxdepth < 0)))
		{
			drawTree(idraw, node->childs[0], depth + 1, ncolor, lcolor, mindepth, maxdepth);
			drawTree(idraw, node->childs[1], depth + 1, ncolor, lcolor, mindepth, maxdepth);
		}
		if (depth >= mindepth)
		{
			const btScalar scl = (btScalar)(node->isinternal() ? 1 : 1);
			const btVector3 mi = node->volume.Center() - node->volume.Extents() * scl;
			const btVector3 mx = node->volume.Center() + node->volume.Extents() * scl;
			drawBox(idraw, mi, mx, node->isleaf() ? lcolor : ncolor);
		}
	}
}

//
template <typename T>
static inline T sum(const btAlignedObjectArray<T>& items)
{
	T v;
	if (items.size())
	{
		v = items[0];
		for (int i = 1, ni = items.size(); i < ni; ++i)
		{
			v += items[i];
		}
	}
	return (v);
}

//
template <typename T, typename Q>
static inline void add(btAlignedObjectArray<T>& items, const Q& value)
{
	for (int i = 0, ni = items.size(); i < ni; ++i)
	{
		items[i] += value;
	}
}

//
template <typename T, typename Q>
static inline void mul(btAlignedObjectArray<T>& items, const Q& value)
{
	for (int i = 0, ni = items.size(); i < ni; ++i)
	{
		items[i] *= value;
	}
}

//
template <typename T>
static inline T average(const btAlignedObjectArray<T>& items)
{
	const btScalar n = (btScalar)(items.size() > 0 ? items.size() : 1);
	return (sum(items) / n);
}

#if 0
//
 inline static btScalar		tetravolume(const btVector3& x0,
										const btVector3& x1,
										const btVector3& x2,
										const btVector3& x3)
{
	const btVector3	a=x1-x0;
	const btVector3	b=x2-x0;
	const btVector3	c=x3-x0;
	return(btDot(a,btCross(b,c)));
}
#endif

//
#if 0
static btVector3		stresscolor(btScalar stress)
{
	static const btVector3	spectrum[]=	{	btVector3(1,0,1),
		btVector3(0,0,1),
		btVector3(0,1,1),
		btVector3(0,1,0),
		btVector3(1,1,0),
		btVector3(1,0,0),
		btVector3(1,0,0)};
	static const int		ncolors=sizeof(spectrum)/sizeof(spectrum[0])-1;
	static const btScalar	one=1;
	stress=btMax<btScalar>(0,btMin<btScalar>(1,stress))*ncolors;
	const int				sel=(int)stress;
	const btScalar			frc=stress-sel;
	return(spectrum[sel]+(spectrum[sel+1]-spectrum[sel])*frc);
}
#endif

//
void btSoftBodyHelpers::Draw(btSoftBody* psb,
							 btIDebugDraw* idraw,
							 int drawflags)
{
	const btScalar scl = (btScalar)0.1;
	const btScalar nscl = scl * 5;
	const btVector3 lcolor = btVector3(0, 0, 0);
	const btVector3 ncolor = btVector3(1, 1, 1);
	const btVector3 ccolor = btVector3(1, 0, 0);
	int i, j, nj;

	/* Clusters	*/
	if (0 != (drawflags & fDrawFlags::Clusters))
	{
		srand(1806);
		for (i = 0; i < psb->m_clusters.size(); ++i)
		{
			if (psb->m_clusters[i]->m_collide)
			{
				btVector3 color(rand() / (btScalar)RAND_MAX,
								rand() / (btScalar)RAND_MAX,
								rand() / (btScalar)RAND_MAX);
				color = color.normalized() * 0.75;
				btAlignedObjectArray<btVector3> vertices;
				vertices.resize(psb->m_clusters[i]->m_nodes.size());
				for (j = 0, nj = vertices.size(); j < nj; ++j)
				{
					vertices[j] = psb->m_clusters[i]->m_nodes[j]->m_x;
				}
#define USE_NEW_CONVEX_HULL_COMPUTER
#ifdef USE_NEW_CONVEX_HULL_COMPUTER
				btConvexHullComputer computer;
				int stride = sizeof(btVector3);
				int count = vertices.size();
				btScalar shrink = 0.f;
				btScalar shrinkClamp = 0.f;
				computer.compute(&vertices[0].getX(), stride, count, shrink, shrinkClamp);
				for (int i = 0; i < computer.faces.size(); i++)
				{
					int face = computer.faces[i];
					//printf("face=%d\n",face);
					const btConvexHullComputer::Edge* firstEdge = &computer.edges[face];
					const btConvexHullComputer::Edge* edge = firstEdge->getNextEdgeOfFace();

					int v0 = firstEdge->getSourceVertex();
					int v1 = firstEdge->getTargetVertex();
					while (edge != firstEdge)
					{
						int v2 = edge->getTargetVertex();
						idraw->drawTriangle(computer.vertices[v0], computer.vertices[v1], computer.vertices[v2], color, 1);
						edge = edge->getNextEdgeOfFace();
						v0 = v1;
						v1 = v2;
					};
				}
#else

				HullDesc hdsc(QF_TRIANGLES, vertices.size(), &vertices[0]);
				HullResult hres;
				HullLibrary hlib;
				hdsc.mMaxVertices = vertices.size();
				hlib.CreateConvexHull(hdsc, hres);
				const btVector3 center = average(hres.m_OutputVertices);
				add(hres.m_OutputVertices, -center);
				mul(hres.m_OutputVertices, (btScalar)1);
				add(hres.m_OutputVertices, center);
				for (j = 0; j < (int)hres.mNumFaces; ++j)
				{
					const int idx[] = {hres.m_Indices[j * 3 + 0], hres.m_Indices[j * 3 + 1], hres.m_Indices[j * 3 + 2]};
					idraw->drawTriangle(hres.m_OutputVertices[idx[0]],
										hres.m_OutputVertices[idx[1]],
										hres.m_OutputVertices[idx[2]],
										color, 1);
				}
				hlib.ReleaseResult(hres);
#endif
			}
			/* Velocities	*/
#if 0
			for(int j=0;j<psb->m_clusters[i].m_nodes.size();++j)
			{
				const btSoftBody::Cluster&	c=psb->m_clusters[i];
				const btVector3				r=c.m_nodes[j]->m_x-c.m_com;
				const btVector3				v=c.m_lv+btCross(c.m_av,r);
				idraw->drawLine(c.m_nodes[j]->m_x,c.m_nodes[j]->m_x+v,btVector3(1,0,0));
			}
#endif
			/* Frame		*/
			//		btSoftBody::Cluster& c=*psb->m_clusters[i];
			//		idraw->drawLine(c.m_com,c.m_framexform*btVector3(10,0,0),btVector3(1,0,0));
			//		idraw->drawLine(c.m_com,c.m_framexform*btVector3(0,10,0),btVector3(0,1,0));
			//		idraw->drawLine(c.m_com,c.m_framexform*btVector3(0,0,10),btVector3(0,0,1));
		}
	}
	else
	{
		/* Nodes	*/
		if (0 != (drawflags & fDrawFlags::Nodes))
		{
			for (i = 0; i < psb->m_nodes.size(); ++i)
			{
				const btSoftBody::Node& n = psb->m_nodes[i];
				if (0 == (n.m_material->m_flags & btSoftBody::fMaterial::DebugDraw)) continue;
				idraw->drawLine(n.m_x - btVector3(scl, 0, 0), n.m_x + btVector3(scl, 0, 0), btVector3(1, 0, 0));
				idraw->drawLine(n.m_x - btVector3(0, scl, 0), n.m_x + btVector3(0, scl, 0), btVector3(0, 1, 0));
				idraw->drawLine(n.m_x - btVector3(0, 0, scl), n.m_x + btVector3(0, 0, scl), btVector3(0, 0, 1));
			}
		}
		/* Links	*/
		if (0 != (drawflags & fDrawFlags::Links))
		{
			for (i = 0; i < psb->m_links.size(); ++i)
			{
				const btSoftBody::Link& l = psb->m_links[i];
				if (0 == (l.m_material->m_flags & btSoftBody::fMaterial::DebugDraw)) continue;
				idraw->drawLine(l.m_n[0]->m_x, l.m_n[1]->m_x, lcolor);
			}
		}
		/* Normals	*/
		if (0 != (drawflags & fDrawFlags::Normals))
		{
			for (i = 0; i < psb->m_nodes.size(); ++i)
			{
				const btSoftBody::Node& n = psb->m_nodes[i];
				if (0 == (n.m_material->m_flags & btSoftBody::fMaterial::DebugDraw)) continue;
				const btVector3 d = n.m_n * nscl;
				idraw->drawLine(n.m_x, n.m_x + d, ncolor);
				idraw->drawLine(n.m_x, n.m_x - d, ncolor * 0.5);
			}
		}
		/* Contacts	*/
		if (0 != (drawflags & fDrawFlags::Contacts))
		{
			static const btVector3 axis[] = {btVector3(1, 0, 0),
											 btVector3(0, 1, 0),
											 btVector3(0, 0, 1)};
			for (i = 0; i < psb->m_rcontacts.size(); ++i)
			{
				const btSoftBody::RContact& c = psb->m_rcontacts[i];
				const btVector3 o = c.m_node->m_x - c.m_cti.m_normal *
														(btDot(c.m_node->m_x, c.m_cti.m_normal) + c.m_cti.m_offset);
				const btVector3 x = btCross(c.m_cti.m_normal, axis[c.m_cti.m_normal.minAxis()]).normalized();
				const btVector3 y = btCross(x, c.m_cti.m_normal).normalized();
				idraw->drawLine(o - x * nscl, o + x * nscl, ccolor);
				idraw->drawLine(o - y * nscl, o + y * nscl, ccolor);
				idraw->drawLine(o, o + c.m_cti.m_normal * nscl * 3, btVector3(1, 1, 0));
			}
		}
		/* Faces	*/
		if (0 != (drawflags & fDrawFlags::Faces))
		{
			const btScalar scl = (btScalar)0.8;
			const btScalar alp = (btScalar)1;
			const btVector3 col(0, (btScalar)0.7, 0);
			for (i = 0; i < psb->m_faces.size(); ++i)
			{
				const btSoftBody::Face& f = psb->m_faces[i];
				if (0 == (f.m_material->m_flags & btSoftBody::fMaterial::DebugDraw)) continue;
				const btVector3 x[] = {f.m_n[0]->m_x, f.m_n[1]->m_x, f.m_n[2]->m_x};
				const btVector3 c = (x[0] + x[1] + x[2]) / 3;
				idraw->drawTriangle((x[0] - c) * scl + c,
									(x[1] - c) * scl + c,
									(x[2] - c) * scl + c,
									col, alp);
			}
		}
		/* Tetras	*/
		if (0 != (drawflags & fDrawFlags::Tetras))
		{
			const btScalar scl = (btScalar)0.8;
			const btScalar alp = (btScalar)1;
			const btVector3 col((btScalar)0.3, (btScalar)0.3, (btScalar)0.7);
			for (int i = 0; i < psb->m_tetras.size(); ++i)
			{
				const btSoftBody::Tetra& t = psb->m_tetras[i];
				if (0 == (t.m_material->m_flags & btSoftBody::fMaterial::DebugDraw)) continue;
				const btVector3 x[] = {t.m_n[0]->m_x, t.m_n[1]->m_x, t.m_n[2]->m_x, t.m_n[3]->m_x};
				const btVector3 c = (x[0] + x[1] + x[2] + x[3]) / 4;
				idraw->drawTriangle((x[0] - c) * scl + c, (x[1] - c) * scl + c, (x[2] - c) * scl + c, col, alp);
				idraw->drawTriangle((x[0] - c) * scl + c, (x[1] - c) * scl + c, (x[3] - c) * scl + c, col, alp);
				idraw->drawTriangle((x[1] - c) * scl + c, (x[2] - c) * scl + c, (x[3] - c) * scl + c, col, alp);
				idraw->drawTriangle((x[2] - c) * scl + c, (x[0] - c) * scl + c, (x[3] - c) * scl + c, col, alp);
			}
		}
	}
	/* Anchors	*/
	if (0 != (drawflags & fDrawFlags::Anchors))
	{
		for (i = 0; i < psb->m_anchors.size(); ++i)
		{
			const btSoftBody::Anchor& a = psb->m_anchors[i];
			const btVector3 q = a.m_body->getWorldTransform() * a.m_local;
			drawVertex(idraw, a.m_node->m_x, 0.25, btVector3(1, 0, 0));
			drawVertex(idraw, q, 0.25, btVector3(0, 1, 0));
			idraw->drawLine(a.m_node->m_x, q, btVector3(1, 1, 1));
		}
		for (i = 0; i < psb->m_nodes.size(); ++i)
		{
			const btSoftBody::Node& n = psb->m_nodes[i];
			if (0 == (n.m_material->m_flags & btSoftBody::fMaterial::DebugDraw)) continue;
			if (n.m_im <= 0)
			{
				drawVertex(idraw, n.m_x, 0.25, btVector3(1, 0, 0));
			}
		}
	}

	/* Notes	*/
	if (0 != (drawflags & fDrawFlags::Notes))
	{
		for (i = 0; i < psb->m_notes.size(); ++i)
		{
			const btSoftBody::Note& n = psb->m_notes[i];
			btVector3 p = n.m_offset;
			for (int j = 0; j < n.m_rank; ++j)
			{
				p += n.m_nodes[j]->m_x * n.m_coords[j];
			}
			idraw->draw3dText(p, n.m_text);
		}
	}
	/* Node tree	*/
	if (0 != (drawflags & fDrawFlags::NodeTree)) DrawNodeTree(psb, idraw);
	/* Face tree	*/
	if (0 != (drawflags & fDrawFlags::FaceTree)) DrawFaceTree(psb, idraw);
	/* Cluster tree	*/
	if (0 != (drawflags & fDrawFlags::ClusterTree)) DrawClusterTree(psb, idraw);
	/* Joints		*/
	if (0 != (drawflags & fDrawFlags::Joints))
	{
		for (i = 0; i < psb->m_joints.size(); ++i)
		{
			const btSoftBody::Joint* pj = psb->m_joints[i];
			switch (pj->Type())
			{
				case btSoftBody::Joint::eType::Linear:
				{
					const btSoftBody::LJoint* pjl = (const btSoftBody::LJoint*)pj;
					const btVector3 a0 = pj->m_bodies[0].xform() * pjl->m_refs[0];
					const btVector3 a1 = pj->m_bodies[1].xform() * pjl->m_refs[1];
					idraw->drawLine(pj->m_bodies[0].xform().getOrigin(), a0, btVector3(1, 1, 0));
					idraw->drawLine(pj->m_bodies[1].xform().getOrigin(), a1, btVector3(0, 1, 1));
					drawVertex(idraw, a0, 0.25, btVector3(1, 1, 0));
					drawVertex(idraw, a1, 0.25, btVector3(0, 1, 1));
				}
				break;
				case btSoftBody::Joint::eType::Angular:
				{
					//const btSoftBody::AJoint*	pja=(const btSoftBody::AJoint*)pj;
					const btVector3 o0 = pj->m_bodies[0].xform().getOrigin();
					const btVector3 o1 = pj->m_bodies[1].xform().getOrigin();
					const btVector3 a0 = pj->m_bodies[0].xform().getBasis() * pj->m_refs[0];
					const btVector3 a1 = pj->m_bodies[1].xform().getBasis() * pj->m_refs[1];
					idraw->drawLine(o0, o0 + a0 * 10, btVector3(1, 1, 0));
					idraw->drawLine(o0, o0 + a1 * 10, btVector3(1, 1, 0));
					idraw->drawLine(o1, o1 + a0 * 10, btVector3(0, 1, 1));
					idraw->drawLine(o1, o1 + a1 * 10, btVector3(0, 1, 1));
					break;
				}
				default:
				{
				}
			}
		}
	}
}

//
void btSoftBodyHelpers::DrawInfos(btSoftBody* psb,
								  btIDebugDraw* idraw,
								  bool masses,
								  bool areas,
								  bool /*stress*/)
{
	for (int i = 0; i < psb->m_nodes.size(); ++i)
	{
		const btSoftBody::Node& n = psb->m_nodes[i];
		char text[2048] = {0};
		char buff[1024];
		if (masses)
		{
			sprintf(buff, " M(%.2f)", 1 / n.m_im);
			strcat(text, buff);
		}
		if (areas)
		{
			sprintf(buff, " A(%.2f)", n.m_area);
			strcat(text, buff);
		}
		if (text[0]) idraw->draw3dText(n.m_x, text);
	}
}

//
void btSoftBodyHelpers::DrawNodeTree(btSoftBody* psb,
									 btIDebugDraw* idraw,
									 int mindepth,
									 int maxdepth)
{
	drawTree(idraw, psb->m_ndbvt.m_root, 0, btVector3(1, 0, 1), btVector3(1, 1, 1), mindepth, maxdepth);
}

//
void btSoftBodyHelpers::DrawFaceTree(btSoftBody* psb,
									 btIDebugDraw* idraw,
									 int mindepth,
									 int maxdepth)
{
	drawTree(idraw, psb->m_fdbvt.m_root, 0, btVector3(0, 1, 0), btVector3(1, 0, 0), mindepth, maxdepth);
}

//
void btSoftBodyHelpers::DrawClusterTree(btSoftBody* psb,
										btIDebugDraw* idraw,
										int mindepth,
										int maxdepth)
{
	drawTree(idraw, psb->m_cdbvt.m_root, 0, btVector3(0, 1, 1), btVector3(1, 0, 0), mindepth, maxdepth);
}

//The btSoftBody object from the BulletSDK includes an array of Nodes and Links. These links appear
// to be first set up to connect a node to between 5 and 6 of its neighbors [480 links],
//and then to the rest of the nodes after the execution of the Floyd-Warshall graph algorithm
//[another 930 links].
//The way the links are stored by default, we have a number of cases where adjacent links share a node in common
// - this leads to the creation of a data dependency through memory.
//The PSolve_Links() function reads and writes nodes as it iterates over each link.
//So, we now have the possibility of a data dependency between iteration X
//that processes link L with iteration X+1 that processes link L+1
//because L and L+1 have one node in common, and iteration X updates the positions of that node,
//and iteration X+1 reads in the position of that shared node.
//
//Such a memory dependency limits the ability of a modern CPU to speculate beyond
//a certain point because it has to respect a possible dependency
//- this prevents the CPU from making full use of its out-of-order resources.
//If we re-order the links such that we minimize the cases where a link L and L+1 share a common node,
//we create a temporal gap between when the node position is written,
//and when it is subsequently read. This in turn allows the CPU to continue execution without
//risking a dependency violation. Such a reordering would result in significant speedups on
//modern CPUs with lots of execution resources.
//In our testing, we see it have a tremendous impact not only on the A7,
//but also on all x86 cores that ship with modern Macs.
//The attached source file includes a single function (ReoptimizeLinkOrder) which can be called on a
//btSoftBody object in the solveConstraints() function before the actual solver is invoked,
//or right after generateBendingConstraints() once we have all 1410 links.

//===================================================================
//
//
// This function takes in a list of interdependent Links and tries
// to maximize the distance between calculation
// of dependent links.  This increases the amount of parallelism that can
// be exploited by out-of-order instruction processors with large but
// (inevitably) finite instruction windows.
//
//===================================================================

// A small structure to track lists of dependent link calculations
class LinkDeps_t
{
public:
	int value;         // A link calculation that is dependent on this one
					   // Positive values = "input A" while negative values = "input B"
	LinkDeps_t* next;  // Next dependence in the list
};
typedef LinkDeps_t* LinkDepsPtr_t;

// Dependency list constants
#define REOP_NOT_DEPENDENT -1
#define REOP_NODE_COMPLETE -2  // Must be less than REOP_NOT_DEPENDENT

void btSoftBodyHelpers::ReoptimizeLinkOrder(btSoftBody* psb /* This can be replaced by a btSoftBody pointer */)
{
	int i, nLinks = psb->m_links.size(), nNodes = psb->m_nodes.size();
	btSoftBody::Link* lr;
	int ar, br;
	btSoftBody::Node* node0 = &(psb->m_nodes[0]);
	btSoftBody::Node* node1 = &(psb->m_nodes[1]);
	LinkDepsPtr_t linkDep;
	int readyListHead, readyListTail, linkNum, linkDepFrees, depLink;

	// Allocate temporary buffers
	int* nodeWrittenAt = new int[nNodes + 1];  // What link calculation produced this node's current values?
	int* linkDepA = new int[nLinks];           // Link calculation input is dependent upon prior calculation #N
	int* linkDepB = new int[nLinks];
	int* readyList = new int[nLinks];                              // List of ready-to-process link calculations (# of links, maximum)
	LinkDeps_t* linkDepFreeList = new LinkDeps_t[2 * nLinks];      // Dependent-on-me list elements (2x# of links, maximum)
	LinkDepsPtr_t* linkDepListStarts = new LinkDepsPtr_t[nLinks];  // Start nodes of dependent-on-me lists, one for each link

	// Copy the original, unsorted links to a side buffer
	btSoftBody::Link* linkBuffer = new btSoftBody::Link[nLinks];
	memcpy(linkBuffer, &(psb->m_links[0]), sizeof(btSoftBody::Link) * nLinks);

	// Clear out the node setup and ready list
	for (i = 0; i < nNodes + 1; i++)
	{
		nodeWrittenAt[i] = REOP_NOT_DEPENDENT;
	}
	for (i = 0; i < nLinks; i++)
	{
		linkDepListStarts[i] = NULL;
	}
	readyListHead = readyListTail = linkDepFrees = 0;

	// Initial link analysis to set up data structures
	for (i = 0; i < nLinks; i++)
	{
		// Note which prior link calculations we are dependent upon & build up dependence lists
		lr = &(psb->m_links[i]);
		ar = (lr->m_n[0] - node0) / (node1 - node0);
		br = (lr->m_n[1] - node0) / (node1 - node0);
		if (nodeWrittenAt[ar] > REOP_NOT_DEPENDENT)
		{
			linkDepA[i] = nodeWrittenAt[ar];
			linkDep = &linkDepFreeList[linkDepFrees++];
			linkDep->value = i;
			linkDep->next = linkDepListStarts[nodeWrittenAt[ar]];
			linkDepListStarts[nodeWrittenAt[ar]] = linkDep;
		}
		else
		{
			linkDepA[i] = REOP_NOT_DEPENDENT;
		}
		if (nodeWrittenAt[br] > REOP_NOT_DEPENDENT)
		{
			linkDepB[i] = nodeWrittenAt[br];
			linkDep = &linkDepFreeList[linkDepFrees++];
			linkDep->value = -(i + 1);
			linkDep->next = linkDepListStarts[nodeWrittenAt[br]];
			linkDepListStarts[nodeWrittenAt[br]] = linkDep;
		}
		else
		{
			linkDepB[i] = REOP_NOT_DEPENDENT;
		}

		// Add this link to the initial ready list, if it is not dependent on any other links
		if ((linkDepA[i] == REOP_NOT_DEPENDENT) && (linkDepB[i] == REOP_NOT_DEPENDENT))
		{
			readyList[readyListTail++] = i;
			linkDepA[i] = linkDepB[i] = REOP_NODE_COMPLETE;  // Probably not needed now
		}

		// Update the nodes to mark which ones are calculated by this link
		nodeWrittenAt[ar] = nodeWrittenAt[br] = i;
	}

	// Process the ready list and create the sorted list of links
	// -- By treating the ready list as a queue, we maximize the distance between any
	//    inter-dependent node calculations
	// -- All other (non-related) nodes in the ready list will automatically be inserted
	//    in between each set of inter-dependent link calculations by this loop
	i = 0;
	while (readyListHead != readyListTail)
	{
		// Use ready list to select the next link to process
		linkNum = readyList[readyListHead++];
		// Copy the next-to-calculate link back into the original link array
		psb->m_links[i++] = linkBuffer[linkNum];

		// Free up any link inputs that are dependent on this one
		linkDep = linkDepListStarts[linkNum];
		while (linkDep)
		{
			depLink = linkDep->value;
			if (depLink >= 0)
			{
				linkDepA[depLink] = REOP_NOT_DEPENDENT;
			}
			else
			{
				depLink = -depLink - 1;
				linkDepB[depLink] = REOP_NOT_DEPENDENT;
			}
			// Add this dependent link calculation to the ready list if *both* inputs are clear
			if ((linkDepA[depLink] == REOP_NOT_DEPENDENT) && (linkDepB[depLink] == REOP_NOT_DEPENDENT))
			{
				readyList[readyListTail++] = depLink;
				linkDepA[depLink] = linkDepB[depLink] = REOP_NODE_COMPLETE;  // Probably not needed now
			}
			linkDep = linkDep->next;
		}
	}

	// Delete the temporary buffers
	delete[] nodeWrittenAt;
	delete[] linkDepA;
	delete[] linkDepB;
	delete[] readyList;
	delete[] linkDepFreeList;
	delete[] linkDepListStarts;
	delete[] linkBuffer;
}

//
void btSoftBodyHelpers::DrawFrame(btSoftBody* psb,
								  btIDebugDraw* idraw)
{
	if (psb->m_pose.m_bframe)
	{
		static const btScalar ascl = 10;
		static const btScalar nscl = (btScalar)0.1;
		const btVector3 com = psb->m_pose.m_com;
		const btMatrix3x3 trs = psb->m_pose.m_rot * psb->m_pose.m_scl;
		const btVector3 Xaxis = (trs * btVector3(1, 0, 0)).normalized();
		const btVector3 Yaxis = (trs * btVector3(0, 1, 0)).normalized();
		const btVector3 Zaxis = (trs * btVector3(0, 0, 1)).normalized();
		idraw->drawLine(com, com + Xaxis * ascl, btVector3(1, 0, 0));
		idraw->drawLine(com, com + Yaxis * ascl, btVector3(0, 1, 0));
		idraw->drawLine(com, com + Zaxis * ascl, btVector3(0, 0, 1));
		for (int i = 0; i < psb->m_pose.m_pos.size(); ++i)
		{
			const btVector3 x = com + trs * psb->m_pose.m_pos[i];
			drawVertex(idraw, x, nscl, btVector3(1, 0, 1));
		}
	}
}

//
btSoftBody* btSoftBodyHelpers::CreateRope(btSoftBodyWorldInfo& worldInfo, const btVector3& from,
										  const btVector3& to,
										  int res,
										  int fixeds)
{
	/* Create nodes	*/
	const int r = res + 2;
	btVector3* x = new btVector3[r];
	btScalar* m = new btScalar[r];
	int i;

	for (i = 0; i < r; ++i)
	{
		const btScalar t = i / (btScalar)(r - 1);
		x[i] = lerp(from, to, t);
		m[i] = 1;
	}
	btSoftBody* psb = new btSoftBody(&worldInfo, r, x, m);
	if (fixeds & 1) psb->setMass(0, 0);
	if (fixeds & 2) psb->setMass(r - 1, 0);
	delete[] x;
	delete[] m;
	/* Create links	*/
	for (i = 1; i < r; ++i)
	{
		psb->appendLink(i - 1, i);
	}
	/* Finished		*/
	return (psb);
}

//
btSoftBody* btSoftBodyHelpers::CreatePatch(btSoftBodyWorldInfo& worldInfo, const btVector3& corner00,
										   const btVector3& corner10,
										   const btVector3& corner01,
										   const btVector3& corner11,
										   int resx,
										   int resy,
										   int fixeds,
										   bool gendiags,
										   btScalar perturbation)
{
#define IDX(_x_, _y_) ((_y_) * rx + (_x_))
	/* Create nodes	*/
	if ((resx < 2) || (resy < 2)) return (0);
	const int rx = resx;
	const int ry = resy;
	const int tot = rx * ry;
	btVector3* x = new btVector3[tot];
	btScalar* m = new btScalar[tot];
	int iy;

	for (iy = 0; iy < ry; ++iy)
	{
		const btScalar ty = iy / (btScalar)(ry - 1);
		const btVector3 py0 = lerp(corner00, corner01, ty);
		const btVector3 py1 = lerp(corner10, corner11, ty);
		for (int ix = 0; ix < rx; ++ix)
		{
			const btScalar tx = ix / (btScalar)(rx - 1);
			btScalar pert = perturbation * btScalar(rand()) / RAND_MAX;
			btVector3 temp1 = py1;
			temp1.setY(py1.getY() + pert);
			btVector3 temp = py0;
			pert = perturbation * btScalar(rand()) / RAND_MAX;
			temp.setY(py0.getY() + pert);
			x[IDX(ix, iy)] = lerp(temp, temp1, tx);
			m[IDX(ix, iy)] = 1;
		}
	}
	btSoftBody* psb = new btSoftBody(&worldInfo, tot, x, m);
	if (fixeds & 1) psb->setMass(IDX(0, 0), 0);
	if (fixeds & 2) psb->setMass(IDX(rx - 1, 0), 0);
	if (fixeds & 4) psb->setMass(IDX(0, ry - 1), 0);
	if (fixeds & 8) psb->setMass(IDX(rx - 1, ry - 1), 0);
	delete[] x;
	delete[] m;
	/* Create links	and faces */
	for (iy = 0; iy < ry; ++iy)
	{
		for (int ix = 0; ix < rx; ++ix)
		{
			const int idx = IDX(ix, iy);
			const bool mdx = (ix + 1) < rx;
			const bool mdy = (iy + 1) < ry;
			if (mdx) psb->appendLink(idx, IDX(ix + 1, iy));
			if (mdy) psb->appendLink(idx, IDX(ix, iy + 1));
			if (mdx && mdy)
			{
				if ((ix + iy) & 1)
				{
					psb->appendFace(IDX(ix, iy), IDX(ix + 1, iy), IDX(ix + 1, iy + 1));
					psb->appendFace(IDX(ix, iy), IDX(ix + 1, iy + 1), IDX(ix, iy + 1));
					if (gendiags)
					{
						psb->appendLink(IDX(ix, iy), IDX(ix + 1, iy + 1));
					}
				}
				else
				{
					psb->appendFace(IDX(ix, iy + 1), IDX(ix, iy), IDX(ix + 1, iy));
					psb->appendFace(IDX(ix, iy + 1), IDX(ix + 1, iy), IDX(ix + 1, iy + 1));
					if (gendiags)
					{
						psb->appendLink(IDX(ix + 1, iy), IDX(ix, iy + 1));
					}
				}
			}
		}
	}
	/* Finished		*/
#undef IDX
	return (psb);
}

//
btSoftBody* btSoftBodyHelpers::CreatePatchUV(btSoftBodyWorldInfo& worldInfo,
											 const btVector3& corner00,
											 const btVector3& corner10,
											 const btVector3& corner01,
											 const btVector3& corner11,
											 int resx,
											 int resy,
											 int fixeds,
											 bool gendiags,
											 float* tex_coords)
{
	/*
	*
	*  corners:
	*
	*  [0][0]     corner00 ------- corner01   [resx][0]
	*                |                |
	*                |                |
	*  [0][resy]  corner10 -------- corner11  [resx][resy]
	*
	*
	*
	*
	*
	*
	*   "fixedgs" map:
	*
	*  corner00     -->   +1
	*  corner01     -->   +2
	*  corner10     -->   +4
	*  corner11     -->   +8
	*  upper middle -->  +16
	*  left middle  -->  +32
	*  right middle -->  +64
	*  lower middle --> +128
	*  center       --> +256
	*
	*
	*   tex_coords size   (resx-1)*(resy-1)*12
	*
	*
	*
	*     SINGLE QUAD INTERNALS
	*
	*  1) btSoftBody's nodes and links,
	*     diagonal link is optional ("gendiags")
	*
	*
	*    node00 ------ node01
	*      | .              
	*      |   .            
	*      |     .          
	*      |       .        
	*      |         .      
	*    node10        node11
	*
	*
	*
	*   2) Faces:
	*      two triangles,
	*      UV Coordinates (hier example for single quad)
	*      
	*     (0,1)          (0,1)  (1,1)
	*     1 |\            3 \-----| 2
	*       | \              \    |
	*       |  \              \   |
	*       |   \              \  |
	*       |    \              \ |
	*     2 |-----\ 3            \| 1
	*     (0,0)    (1,0)       (1,0)
	*
	*
	*
	*
	*
	*
	*/

#define IDX(_x_, _y_) ((_y_) * rx + (_x_))
	/* Create nodes		*/
	if ((resx < 2) || (resy < 2)) return (0);
	const int rx = resx;
	const int ry = resy;
	const int tot = rx * ry;
	btVector3* x = new btVector3[tot];
	btScalar* m = new btScalar[tot];

	int iy;

	for (iy = 0; iy < ry; ++iy)
	{
		const btScalar ty = iy / (btScalar)(ry - 1);
		const btVector3 py0 = lerp(corner00, corner01, ty);
		const btVector3 py1 = lerp(corner10, corner11, ty);
		for (int ix = 0; ix < rx; ++ix)
		{
			const btScalar tx = ix / (btScalar)(rx - 1);
			x[IDX(ix, iy)] = lerp(py0, py1, tx);
			m[IDX(ix, iy)] = 1;
		}
	}
	btSoftBody* psb = new btSoftBody(&worldInfo, tot, x, m);
	if (fixeds & 1) psb->setMass(IDX(0, 0), 0);
	if (fixeds & 2) psb->setMass(IDX(rx - 1, 0), 0);
	if (fixeds & 4) psb->setMass(IDX(0, ry - 1), 0);
	if (fixeds & 8) psb->setMass(IDX(rx - 1, ry - 1), 0);
	if (fixeds & 16) psb->setMass(IDX((rx - 1) / 2, 0), 0);
	if (fixeds & 32) psb->setMass(IDX(0, (ry - 1) / 2), 0);
	if (fixeds & 64) psb->setMass(IDX(rx - 1, (ry - 1) / 2), 0);
	if (fixeds & 128) psb->setMass(IDX((rx - 1) / 2, ry - 1), 0);
	if (fixeds & 256) psb->setMass(IDX((rx - 1) / 2, (ry - 1) / 2), 0);
	delete[] x;
	delete[] m;

	int z = 0;
	/* Create links	and faces	*/
	for (iy = 0; iy < ry; ++iy)
	{
		for (int ix = 0; ix < rx; ++ix)
		{
			const bool mdx = (ix + 1) < rx;
			const bool mdy = (iy + 1) < ry;

			int node00 = IDX(ix, iy);
			int node01 = IDX(ix + 1, iy);
			int node10 = IDX(ix, iy + 1);
			int node11 = IDX(ix + 1, iy + 1);

			if (mdx) psb->appendLink(node00, node01);
			if (mdy) psb->appendLink(node00, node10);
			if (mdx && mdy)
			{
				psb->appendFace(node00, node10, node11);
				if (tex_coords)
				{
					tex_coords[z + 0] = CalculateUV(resx, resy, ix, iy, 0);
					tex_coords[z + 1] = CalculateUV(resx, resy, ix, iy, 1);
					tex_coords[z + 2] = CalculateUV(resx, resy, ix, iy, 0);
					tex_coords[z + 3] = CalculateUV(resx, resy, ix, iy, 2);
					tex_coords[z + 4] = CalculateUV(resx, resy, ix, iy, 3);
					tex_coords[z + 5] = CalculateUV(resx, resy, ix, iy, 2);
				}
				psb->appendFace(node11, node01, node00);
				if (tex_coords)
				{
					tex_coords[z + 6] = CalculateUV(resx, resy, ix, iy, 3);
					tex_coords[z + 7] = CalculateUV(resx, resy, ix, iy, 2);
					tex_coords[z + 8] = CalculateUV(resx, resy, ix, iy, 3);
					tex_coords[z + 9] = CalculateUV(resx, resy, ix, iy, 1);
					tex_coords[z + 10] = CalculateUV(resx, resy, ix, iy, 0);
					tex_coords[z + 11] = CalculateUV(resx, resy, ix, iy, 1);
				}
				if (gendiags) psb->appendLink(node00, node11);
				z += 12;
			}
		}
	}
	/* Finished	*/
#undef IDX
	return (psb);
}

float btSoftBodyHelpers::CalculateUV(int resx, int resy, int ix, int iy, int id)
{
	/*
	*
	*
	*    node00 --- node01
	*      |          |
	*    node10 --- node11
	*
	*
	*   ID map:
	*
	*   node00 s --> 0
	*   node00 t --> 1
	*
	*   node01 s --> 3
	*   node01 t --> 1
	*
	*   node10 s --> 0
	*   node10 t --> 2
	*
	*   node11 s --> 3
	*   node11 t --> 2
	*
	*
	*/

	float tc = 0.0f;
	if (id == 0)
	{
		tc = (1.0f / ((resx - 1)) * ix);
	}
	else if (id == 1)
	{
		tc = (1.0f / ((resy - 1)) * (resy - 1 - iy));
	}
	else if (id == 2)
	{
		tc = (1.0f / ((resy - 1)) * (resy - 1 - iy - 1));
	}
	else if (id == 3)
	{
		tc = (1.0f / ((resx - 1)) * (ix + 1));
	}
	return tc;
}
//
btSoftBody* btSoftBodyHelpers::CreateEllipsoid(btSoftBodyWorldInfo& worldInfo, const btVector3& center,
											   const btVector3& radius,
											   int res)
{
	struct Hammersley
	{
		static void Generate(btVector3* x, int n)
		{
			for (int i = 0; i < n; i++)
			{
				btScalar p = 0.5, t = 0;
				for (int j = i; j; p *= 0.5, j >>= 1)
					if (j & 1) t += p;
				btScalar w = 2 * t - 1;
				btScalar a = (SIMD_PI + 2 * i * SIMD_PI) / n;
				btScalar s = btSqrt(1 - w * w);
				*x++ = btVector3(s * btCos(a), s * btSin(a), w);
			}
		}
	};
	btAlignedObjectArray<btVector3> vtx;
	vtx.resize(3 + res);
	Hammersley::Generate(&vtx[0], vtx.size());
	for (int i = 0; i < vtx.size(); ++i)
	{
		vtx[i] = vtx[i] * radius + center;
	}
	return (CreateFromConvexHull(worldInfo, &vtx[0], vtx.size()));
}

//
btSoftBody* btSoftBodyHelpers::CreateFromTriMesh(btSoftBodyWorldInfo& worldInfo, const btScalar* vertices,
												 const int* triangles,
												 int ntriangles, bool randomizeConstraints)
{
	int maxidx = 0;
	int i, j, ni;

	for (i = 0, ni = ntriangles * 3; i < ni; ++i)
	{
		maxidx = btMax(triangles[i], maxidx);
	}
	++maxidx;
	btAlignedObjectArray<bool> chks;
	btAlignedObjectArray<btVector3> vtx;
	chks.resize(maxidx * maxidx, false);
	vtx.resize(maxidx);
	for (i = 0, j = 0, ni = maxidx * 3; i < ni; ++j, i += 3)
	{
		vtx[j] = btVector3(vertices[i], vertices[i + 1], vertices[i + 2]);
	}
	btSoftBody* psb = new btSoftBody(&worldInfo, vtx.size(), &vtx[0], 0);
	for (i = 0, ni = ntriangles * 3; i < ni; i += 3)
	{
		const int idx[] = {triangles[i], triangles[i + 1], triangles[i + 2]};
#define IDX(_x_, _y_) ((_y_) * maxidx + (_x_))
		for (int j = 2, k = 0; k < 3; j = k++)
		{
			if (!chks[IDX(idx[j], idx[k])])
			{
				chks[IDX(idx[j], idx[k])] = true;
				chks[IDX(idx[k], idx[j])] = true;
				psb->appendLink(idx[j], idx[k]);
			}
		}
#undef IDX
		psb->appendFace(idx[0], idx[1], idx[2]);
	}

	if (randomizeConstraints)
	{
		psb->randomizeConstraints();
	}

	return (psb);
}

//
btSoftBody* btSoftBodyHelpers::CreateFromConvexHull(btSoftBodyWorldInfo& worldInfo, const btVector3* vertices,
													int nvertices, bool randomizeConstraints)
{
	HullDesc hdsc(QF_TRIANGLES, nvertices, vertices);
	HullResult hres;
	HullLibrary hlib; /*??*/
	hdsc.mMaxVertices = nvertices;
	hlib.CreateConvexHull(hdsc, hres);
	btSoftBody* psb = new btSoftBody(&worldInfo, (int)hres.mNumOutputVertices,
									 &hres.m_OutputVertices[0], 0);
	for (int i = 0; i < (int)hres.mNumFaces; ++i)
	{
		const int idx[] = {static_cast<int>(hres.m_Indices[i * 3 + 0]),
						   static_cast<int>(hres.m_Indices[i * 3 + 1]),
						   static_cast<int>(hres.m_Indices[i * 3 + 2])};
		if (idx[0] < idx[1]) psb->appendLink(idx[0], idx[1]);
		if (idx[1] < idx[2]) psb->appendLink(idx[1], idx[2]);
		if (idx[2] < idx[0]) psb->appendLink(idx[2], idx[0]);
		psb->appendFace(idx[0], idx[1], idx[2]);
	}
	hlib.ReleaseResult(hres);
	if (randomizeConstraints)
	{
		psb->randomizeConstraints();
	}
	return (psb);
}

void btSoftBodyHelpers::PopulateTetras(btSoftBody* psb, const std::vector<std::array<int, 4>>& tetras, bool createLinks)
{
	for (auto& tetra : tetras)
	{
		psb->appendTetra(tetra[0], tetra[1], tetra[2], tetra[3]);
		if (createLinks)
		{
			psb->appendLink(tetra[0], tetra[1], 0, true);
			psb->appendLink(tetra[1], tetra[2], 0, true);
			psb->appendLink(tetra[2], tetra[0], 0, true);
			psb->appendLink(tetra[0], tetra[3], 0, true);
			psb->appendLink(tetra[1], tetra[3], 0, true);
			psb->appendLink(tetra[2], tetra[3], 0, true);
		}
	}

	btSoftBodyHelpers::generateBoundaryFaces(psb);
	psb->initializeDmInverse();
	psb->m_tetraScratches.resize(psb->m_tetras.size());
	psb->m_tetraScratchesTn.resize(psb->m_tetras.size());
}

static int nextLine(const char* buffer)
{
	int numBytesRead = 0;

	while (*buffer != '\n')
	{
		buffer++;
		numBytesRead++;
	}

	if (buffer[0] == 0x0a)
	{
		buffer++;
		numBytesRead++;
	}
	return numBytesRead;
}

/* Create from TetGen .ele, .face, .node data							*/
btSoftBody* btSoftBodyHelpers::CreateFromTetGenData(btSoftBodyWorldInfo& worldInfo,
													const char* ele,
													const char* face,
													const char* node,
													bool bfacelinks,
													bool btetralinks,
													bool bfacesfromtetras)
{
	btAlignedObjectArray<btVector3> pos;
	int nnode = 0;
	int ndims = 0;
	int nattrb = 0;
	int hasbounds = 0;
	int result = sscanf(node, "%d %d %d %d", &nnode, &ndims, &nattrb, &hasbounds);
	result = sscanf(node, "%d %d %d %d", &nnode, &ndims, &nattrb, &hasbounds);
	node += nextLine(node);

	pos.resize(nnode);
	for (int i = 0; i < pos.size(); ++i)
	{
		int index = 0;
		//int			bound=0;
		float x, y, z;
		sscanf(node, "%d %f %f %f", &index, &x, &y, &z);

		//	sn>>index;
		//	sn>>x;sn>>y;sn>>z;
		node += nextLine(node);

		//for(int j=0;j<nattrb;++j)
		//	sn>>a;

		//if(hasbounds)
		//	sn>>bound;

		pos[index].setX(btScalar(x));
		pos[index].setY(btScalar(y));
		pos[index].setZ(btScalar(z));
	}
	btSoftBody* psb = new btSoftBody(&worldInfo, nnode, &pos[0], 0);
#if 0
if(face&&face[0])
	{
	int								nface=0;
	sf>>nface;sf>>hasbounds;
	for(int i=0;i<nface;++i)
		{
		int			index=0;
		int			bound=0;
		int			ni[3];
		sf>>index;
		sf>>ni[0];sf>>ni[1];sf>>ni[2];
		sf>>bound;
		psb->appendFace(ni[0],ni[1],ni[2]);	
		if(btetralinks)
			{
			psb->appendLink(ni[0],ni[1],0,true);
			psb->appendLink(ni[1],ni[2],0,true);
			psb->appendLink(ni[2],ni[0],0,true);
			}
		}
	}
#endif

	if (ele && ele[0])
	{
		int ntetra = 0;
		int ncorner = 0;
		int neattrb = 0;
		sscanf(ele, "%d %d %d", &ntetra, &ncorner, &neattrb);
		ele += nextLine(ele);

		//se>>ntetra;se>>ncorner;se>>neattrb;
		for (int i = 0; i < ntetra; ++i)
		{
			int index = 0;
			int ni[4];

			//se>>index;
			//se>>ni[0];se>>ni[1];se>>ni[2];se>>ni[3];
			sscanf(ele, "%d %d %d %d %d", &index, &ni[0], &ni[1], &ni[2], &ni[3]);
			ele += nextLine(ele);
			//for(int j=0;j<neattrb;++j)
			//	se>>a;
			psb->appendTetra(ni[0], ni[1], ni[2], ni[3]);
			if (btetralinks)
			{
				psb->appendLink(ni[0], ni[1], 0, true);
				psb->appendLink(ni[1], ni[2], 0, true);
				psb->appendLink(ni[2], ni[0], 0, true);
				psb->appendLink(ni[0], ni[3], 0, true);
				psb->appendLink(ni[1], ni[3], 0, true);
				psb->appendLink(ni[2], ni[3], 0, true);
			}
		}
	}
	psb->initializeDmInverse();
	psb->m_tetraScratches.resize(psb->m_tetras.size());
	psb->m_tetraScratchesTn.resize(psb->m_tetras.size());
	printf("Nodes:  %u\r\n", psb->m_nodes.size());
	printf("Links:  %u\r\n", psb->m_links.size());
	printf("Faces:  %u\r\n", psb->m_faces.size());
	printf("Tetras: %u\r\n", psb->m_tetras.size());
	return (psb);
}

btSoftBody* btSoftBodyHelpers::CreateFromVtkFile(btSoftBodyWorldInfo& worldInfo, const char* vtk_file)
{
	std::ifstream fs;
	fs.open(vtk_file);
	btAssert(fs);

	typedef btAlignedObjectArray<int> Index;
	std::string line;
	btAlignedObjectArray<btVector3> X;
	btVector3 position;
	btAlignedObjectArray<Index> indices;
	bool reading_points = false;
	bool reading_tets = false;
	size_t n_points = 0;
	size_t n_tets = 0;
	size_t x_count = 0;
	size_t indices_count = 0;
	while (std::getline(fs, line))
	{
		std::stringstream ss(line);
		if (line.size() == (size_t)(0))
		{
		}
		else if (line.substr(0, 6) == "POINTS")
		{
			reading_points = true;
			reading_tets = false;
			ss.ignore(128, ' ');  // ignore "POINTS"
			ss >> n_points;
			X.resize(n_points);
		}
		else if (line.substr(0, 5) == "CELLS")
		{
			reading_points = false;
			reading_tets = true;
			ss.ignore(128, ' ');  // ignore "CELLS"
			ss >> n_tets;
			indices.resize(n_tets);
		}
		else if (line.substr(0, 10) == "CELL_TYPES")
		{
			reading_points = false;
			reading_tets = false;
		}
		else if (reading_points)
		{
			btScalar p;
			ss >> p;
			position.setX(p);
			ss >> p;
			position.setY(p);
			ss >> p;
			position.setZ(p);
			//printf("v %f %f %f\n", position.getX(), position.getY(), position.getZ());
			X[x_count++] = position;
		}
		else if (reading_tets)
		{
			int d;
			ss >> d;
			if (d != 4)
			{
				printf("Load deformable failed: Only Tetrahedra are supported in VTK file.\n");
				fs.close();
				return 0;
			}
			ss.ignore(128, ' ');  // ignore "4"
			Index tet;
			tet.resize(4);
			for (size_t i = 0; i < 4; i++)
			{
				ss >> tet[i];
				//printf("%d ", tet[i]);
			}
			//printf("\n");
			indices[indices_count++] = tet;
		}
	}
	btSoftBody* psb = new btSoftBody(&worldInfo, n_points, &X[0], 0);

	for (int i = 0; i < n_tets; ++i)
	{
		const Index& ni = indices[i];
		psb->appendTetra(ni[0], ni[1], ni[2], ni[3]);
		{
			psb->appendLink(ni[0], ni[1], 0, true);
			psb->appendLink(ni[1], ni[2], 0, true);
			psb->appendLink(ni[2], ni[0], 0, true);
			psb->appendLink(ni[0], ni[3], 0, true);
			psb->appendLink(ni[1], ni[3], 0, true);
			psb->appendLink(ni[2], ni[3], 0, true);
		}
	}

	generateBoundaryFaces(psb);
	psb->initializeDmInverse();
	psb->m_tetraScratches.resize(psb->m_tetras.size());
	psb->m_tetraScratchesTn.resize(psb->m_tetras.size());
	printf("Nodes:  %u\r\n", psb->m_nodes.size());
	printf("Links:  %u\r\n", psb->m_links.size());
	printf("Faces:  %u\r\n", psb->m_faces.size());
	printf("Tetras: %u\r\n", psb->m_tetras.size());

	fs.close();
	return psb;
}

void btSoftBodyHelpers::generateBoundaryFaces(btSoftBody* psb)
{
	int counter = 0;
	for (int i = 0; i < psb->m_nodes.size(); ++i)
	{
		psb->m_nodes[i].index = counter++;
		psb->m_nodes[i].local_index = psb->m_nodes[i].index;
	}
	typedef btAlignedObjectArray<int> Index;
	btAlignedObjectArray<Index> indices;
	indices.resize(psb->m_tetras.size());
	for (int i = 0; i < indices.size(); ++i)
	{
		Index index;
		index.push_back(psb->m_tetras[i].m_n[0]->index);
		index.push_back(psb->m_tetras[i].m_n[1]->index);
		index.push_back(psb->m_tetras[i].m_n[2]->index);
		index.push_back(psb->m_tetras[i].m_n[3]->index);
		indices[i] = index;
	}

	// Originally there was simply a std::vector as the map key. I worte this for a signle reason - to dodge one MSVC linker error (caused by MSVC's operator< for comparing two arrays or vectors) when using multiple MSVC versions at once. So, this can be safely deleted over time...
	class FaceIndices
	{
		std::array<int, 3> arr;

	public:
		FaceIndices() {}
		FaceIndices(const std::array<int, 3>& arr) : arr(arr) {}

		// Overload the less-than operator
		bool operator<(const FaceIndices& other) const
		{
			for (int i = 0; i < 3; ++i)
			{
				if (arr[i] < other.arr[i])
				{
					return true;
				}
				if (arr[i] > other.arr[i])
				{
					return false;
				}
			}
			return false;  // If all elements are equal, return false
		}

		void Set(int idx, int val)
		{
			arr[idx] = val;
		}

		void Sort()
		{
			std::sort(arr.begin(), arr.end());
		}

		int Get(int idx) const
		{
			return arr[idx];
		}
	};

	std::map<FaceIndices, FaceIndices> dict;
	for (int i = 0; i < indices.size(); ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			FaceIndices f;
			if (j == 0)
			{
				f.Set(0, indices[i][1]);
				f.Set(1, indices[i][0]);
				f.Set(2, indices[i][2]);
			}
			if (j == 1)
			{
				f.Set(0, indices[i][3]);
				f.Set(1, indices[i][0]);
				f.Set(2, indices[i][1]);
			}
			if (j == 2)
			{
				f.Set(0, indices[i][3]);
				f.Set(1, indices[i][1]);
				f.Set(2, indices[i][2]);
			}
			if (j == 3)
			{
				f.Set(0, indices[i][2]);
				f.Set(1, indices[i][0]);
				f.Set(2, indices[i][3]);
			}
			FaceIndices f_sorted = f;
			f_sorted.Sort();
			if (dict.find(f_sorted) != dict.end())
			{
				dict.erase(f_sorted);
			}
			else
			{
				dict.insert(std::make_pair(f_sorted, f));
			}
		}
	}

	for (auto it = dict.begin(); it != dict.end(); ++it)
	{
		FaceIndices f = it->second;
		psb->appendFace(f.Get(0), f.Get(1), f.Get(2));
		//printf("f %d %d %d\n", f[0] + 1, f[1] + 1, f[2] + 1);
	}
}

//Write the surface mesh to an obj file.
void btSoftBodyHelpers::writeObj(const char* filename, const btSoftBody* psb)
{
	std::ofstream fs;
	fs.open(filename);
	btAssert(fs);

	if (psb->m_tetras.size() > 0)
	{
		// For tetrahedron mesh, we need to re-index the surface mesh for it to be in obj file/
		std::map<int, int> dict;
		for (int i = 0; i < psb->m_faces.size(); i++)
		{
			for (int d = 0; d < 3; d++)
			{
				int index = psb->m_faces[i].m_n[d]->index;
				if (dict.find(index) == dict.end())
				{
					int dict_size = dict.size();
					dict[index] = dict_size;
					fs << "v";
					for (int k = 0; k < 3; k++)
					{
						fs << " " << psb->m_nodes[index].m_x[k];
					}
					fs << "\n";
				}
			}
		}
		// Write surface mesh.
		for (int i = 0; i < psb->m_faces.size(); ++i)
		{
			fs << "f";
			for (int n = 0; n < 3; n++)
			{
				fs << " " << dict[psb->m_faces[i].m_n[n]->index] + 1;
			}
			fs << "\n";
		}
	}
	else
	{
		// For trimesh, directly write out all the nodes and faces.xs
		for (int i = 0; i < psb->m_nodes.size(); ++i)
		{
			fs << "v";
			for (int d = 0; d < 3; d++)
			{
				fs << " " << psb->m_nodes[i].m_x[d];
			}
			fs << "\n";
		}

		for (int i = 0; i < psb->m_faces.size(); ++i)
		{
			fs << "f";
			for (int n = 0; n < 3; n++)
			{
				fs << " " << psb->m_faces[i].m_n[n]->index + 1;
			}
			fs << "\n";
		}
	}
	fs.close();
}

void btSoftBodyHelpers::writeState(const char* file, const btSoftBody* psb)
{
	std::ofstream fs;
	fs.open(file);
	btAssert(fs);
	fs << std::scientific << std::setprecision(16);

	// Only write out for trimesh, directly write out all the nodes and faces.xs
	for (int i = 0; i < psb->m_nodes.size(); ++i)
	{
		fs << "q";
		for (int d = 0; d < 3; d++)
		{
			fs << " " << psb->m_nodes[i].m_q[d];
		}
		fs << "\n";
	}

	for (int i = 0; i < psb->m_nodes.size(); ++i)
	{
		fs << "v";
		for (int d = 0; d < 3; d++)
		{
			fs << " " << psb->m_nodes[i].m_v[d];
		}
		fs << "\n";
	}
	fs.close();
}

void btSoftBodyHelpers::duplicateFaces(const char* filename, const btSoftBody* psb)
{
	std::ifstream fs_read;
	fs_read.open(filename);
	std::string line;
	btVector3 pos;
	btAlignedObjectArray<btAlignedObjectArray<int>> additional_faces;
	while (std::getline(fs_read, line))
	{
		std::stringstream ss(line);
		if (line[0] == 'v')
		{
		}
		else if (line[0] == 'f')
		{
			ss.ignore();
			int id0, id1, id2;
			ss >> id0;
			ss >> id1;
			ss >> id2;
			btAlignedObjectArray<int> new_face;
			new_face.push_back(id1);
			new_face.push_back(id0);
			new_face.push_back(id2);
			additional_faces.push_back(new_face);
		}
	}
	fs_read.close();

	std::ofstream fs_write;
	fs_write.open(filename, std::ios_base::app);
	for (int i = 0; i < additional_faces.size(); ++i)
	{
		fs_write << "f";
		for (int n = 0; n < 3; n++)
		{
			fs_write << " " << additional_faces[i][n];
		}
		fs_write << "\n";
	}
	fs_write.close();
}

// Given a simplex with vertices a,b,c,d, find the barycentric weights of p in this simplex
void btSoftBodyHelpers::getBarycentricWeights(const btVector3& a, const btVector3& b, const btVector3& c, const btVector3& d, const btVector3& p, btVector4& bary)
{
	btVector3 vap = p - a;
	btVector3 vbp = p - b;

	btVector3 vab = b - a;
	btVector3 vac = c - a;
	btVector3 vad = d - a;

	btVector3 vbc = c - b;
	btVector3 vbd = d - b;
	btScalar va6 = (vbp.cross(vbd)).dot(vbc);
	btScalar vb6 = (vap.cross(vac)).dot(vad);
	btScalar vc6 = (vap.cross(vad)).dot(vab);
	btScalar vd6 = (vap.cross(vab)).dot(vac);
	btScalar v6 = btScalar(1) / (vab.cross(vac).dot(vad));
	bary = btVector4(va6 * v6, vb6 * v6, vc6 * v6, vd6 * v6);
}

// Given a simplex with vertices a,b,c, find the barycentric weights of p in this simplex. bary[3] = 0.
void btSoftBodyHelpers::getBarycentricWeights(const btVector3& a, const btVector3& b, const btVector3& c, const btVector3& p, btVector4& bary)
{
	btVector3 v0 = b - a, v1 = c - a, v2 = p - a;
	btScalar d00 = btDot(v0, v0);
	btScalar d01 = btDot(v0, v1);
	btScalar d11 = btDot(v1, v1);
	btScalar d20 = btDot(v2, v0);
	btScalar d21 = btDot(v2, v1);
	btScalar invDenom = 1.0 / (d00 * d11 - d01 * d01);
	bary[1] = (d11 * d20 - d01 * d21) * invDenom;
	bary[2] = (d00 * d21 - d01 * d20) * invDenom;
	bary[0] = 1.0 - bary[1] - bary[2];
	bary[3] = 0;
}

// Iterate through all render nodes to find the simulation tetrahedron that contains the render node and record the barycentric weights
// If the node is not inside any tetrahedron, assign it to the tetrahedron in which the node has the least negative barycentric weight
void btSoftBodyHelpers::interpolateBarycentricWeights(btSoftBody* psb)
{
	psb->m_z.resize(0);
	psb->m_renderNodesInterpolationWeights.resize(psb->m_renderNodes.size());
	psb->m_renderNodesParents.resize(psb->m_renderNodes.size());
	for (int i = 0; i < psb->m_renderNodes.size(); ++i)
	{
		const btVector3& p = psb->m_renderNodes[i].m_x;
		btVector4 bary;
		btVector4 optimal_bary;
		btScalar min_bary_weight = -1e3;
		btAlignedObjectArray<const btSoftBody::Node*> optimal_parents;
		for (int j = 0; j < psb->m_tetras.size(); ++j)
		{
			const btSoftBody::Tetra& t = psb->m_tetras[j];
			getBarycentricWeights(t.m_n[0]->m_x, t.m_n[1]->m_x, t.m_n[2]->m_x, t.m_n[3]->m_x, p, bary);
			btScalar new_min_bary_weight = bary[0];
			for (int k = 1; k < 4; ++k)
			{
				new_min_bary_weight = btMin(new_min_bary_weight, bary[k]);
			}
			if (new_min_bary_weight > min_bary_weight)
			{
				btAlignedObjectArray<const btSoftBody::Node*> parents;
				parents.push_back(t.m_n[0]);
				parents.push_back(t.m_n[1]);
				parents.push_back(t.m_n[2]);
				parents.push_back(t.m_n[3]);
				optimal_parents = parents;
				optimal_bary = bary;
				min_bary_weight = new_min_bary_weight;
				// stop searching if p is inside the tetrahedron at hand
				if (bary[0] >= 0. && bary[1] >= 0. && bary[2] >= 0. && bary[3] >= 0.)
				{
					break;
				}
			}
		}
		psb->m_renderNodesInterpolationWeights[i] = optimal_bary;
		psb->m_renderNodesParents[i] = optimal_parents;
	}
}

// Iterate through all render nodes to find the simulation triangle that's closest to the node in the barycentric sense.
void btSoftBodyHelpers::extrapolateBarycentricWeights(btSoftBody* psb)
{
	psb->m_renderNodesInterpolationWeights.resize(psb->m_renderNodes.size());
	psb->m_renderNodesParents.resize(psb->m_renderNodes.size());
	psb->m_z.resize(psb->m_renderNodes.size());
	for (int i = 0; i < psb->m_renderNodes.size(); ++i)
	{
		const btVector3& p = psb->m_renderNodes[i].m_x;
		btVector4 bary;
		btVector4 optimal_bary;
		btScalar min_bary_weight = -SIMD_INFINITY;
		btAlignedObjectArray<const btSoftBody::Node*> optimal_parents;
		btScalar dist = 0, optimal_dist = 0;
		for (int j = 0; j < psb->m_faces.size(); ++j)
		{
			const btSoftBody::Face& f = psb->m_faces[j];
			btVector3 n = btCross(f.m_n[1]->m_x - f.m_n[0]->m_x, f.m_n[2]->m_x - f.m_n[0]->m_x);
			btVector3 unit_n = n.normalized();
			dist = (p - f.m_n[0]->m_x).dot(unit_n);
			btVector3 proj_p = p - dist * unit_n;
			getBarycentricWeights(f.m_n[0]->m_x, f.m_n[1]->m_x, f.m_n[2]->m_x, proj_p, bary);
			btScalar new_min_bary_weight = bary[0];
			for (int k = 1; k < 3; ++k)
			{
				new_min_bary_weight = btMin(new_min_bary_weight, bary[k]);
			}

			// p is out of the current best triangle, we found a traingle that's better
			bool better_than_closest_outisde = (new_min_bary_weight > min_bary_weight && min_bary_weight < 0.);
			// p is inside of the current best triangle, we found a triangle that's better
			bool better_than_best_inside = (new_min_bary_weight >= 0 && min_bary_weight >= 0 && btFabs(dist) < btFabs(optimal_dist));

			if (better_than_closest_outisde || better_than_best_inside)
			{
				btAlignedObjectArray<const btSoftBody::Node*> parents;
				parents.push_back(f.m_n[0]);
				parents.push_back(f.m_n[1]);
				parents.push_back(f.m_n[2]);
				optimal_parents = parents;
				optimal_bary = bary;
				optimal_dist = dist;
				min_bary_weight = new_min_bary_weight;
			}
		}
		psb->m_renderNodesInterpolationWeights[i] = optimal_bary;
		psb->m_renderNodesParents[i] = optimal_parents;
		psb->m_z[i] = optimal_dist;
	}
}
/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2008 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/


#ifndef BT_SOFT_BODY_HELPERS_H
#define BT_SOFT_BODY_HELPERS_H

#include "btSoftBody.h"
#include <fstream>
#include <string>

//
// Helpers
//

/* fDrawFlags															*/
struct fDrawFlags
{
	enum _
	{
		Nodes = 0x0001,
		Links = 0x0002,
		Faces = 0x0004,
		Tetras = 0x0008,
		Normals = 0x0010,
		Contacts = 0x0020,
		Anchors = 0x0040,
		Notes = 0x0080,
		Clusters = 0x0100,
		NodeTree = 0x0200,
		FaceTree = 0x0400,
		ClusterTree = 0x0800,
		Joints = 0x1000,
		/* presets	*/
		Std = Links + Faces + Tetras + Anchors + Notes + Joints,
		StdTetra = Std - Faces + Tetras
	};
};

struct btSoftBodyHelpers
{
	/* Draw body															*/
	static void Draw(btSoftBody* psb,
					 btIDebugDraw* idraw,
					 int drawflags = fDrawFlags::Std);
	/* Draw body infos														*/
	static void DrawInfos(btSoftBody* psb,
						  btIDebugDraw* idraw,
						  bool masses,
						  bool areas,
						  bool stress);
	/* Draw node tree														*/
	static void DrawNodeTree(btSoftBody* psb,
							 btIDebugDraw* idraw,
							 int mindepth = 0,
							 int maxdepth = -1);
	/* Draw face tree														*/
	static void DrawFaceTree(btSoftBody* psb,
							 btIDebugDraw* idraw,
							 int mindepth = 0,
							 int maxdepth = -1);
	/* Draw cluster tree													*/
	static void DrawClusterTree(btSoftBody* psb,
								btIDebugDraw* idraw,
								int mindepth = 0,
								int maxdepth = -1);
	/* Draw rigid frame														*/
	static void DrawFrame(btSoftBody* psb,
						  btIDebugDraw* idraw);
	/* Create a rope														*/
	static btSoftBody* CreateRope(btSoftBodyWorldInfo& worldInfo,
								  const btVector3& from,
								  const btVector3& to,
								  int res,
								  int fixeds);
	/* Create a patch														*/
	static btSoftBody* CreatePatch(btSoftBodyWorldInfo& worldInfo,
								   const btVector3& corner00,
								   const btVector3& corner10,
								   const btVector3& corner01,
								   const btVector3& corner11,
								   int resx,
								   int resy,
								   int fixeds,
								   bool gendiags,
								   btScalar perturbation = 0.);
	/* Create a patch with UV Texture Coordinates	*/
	static btSoftBody* CreatePatchUV(btSoftBodyWorldInfo& worldInfo,
									 const btVector3& corner00,
									 const btVector3& corner10,
									 const btVector3& corner01,
									 const btVector3& corner11,
									 int resx,
									 int resy,
									 int fixeds,
									 bool gendiags,
									 float* tex_coords = 0);
	static float CalculateUV(int resx, int resy, int ix, int iy, int id);
	/* Create an ellipsoid													*/
	static btSoftBody* CreateEllipsoid(btSoftBodyWorldInfo& worldInfo,
									   const btVector3& center,
									   const btVector3& radius,
									   int res);
	/* Create from trimesh													*/
	static btSoftBody* CreateFromTriMesh(btSoftBodyWorldInfo& worldInfo,
										 const btScalar* vertices,
										 const int* triangles,
										 int ntriangles,
										 bool randomizeConstraints = true);
	/* Create from convex-hull												*/
	static btSoftBody* CreateFromConvexHull(btSoftBodyWorldInfo& worldInfo,
											const btVector3* vertices,
											int nvertices,
											bool randomizeConstraints = true);

	static void PopulateTetras(btSoftBody* psb, const std::vector<std::array<int, 4>>& tetras, bool createLinks);

	/* Export TetGen compatible .smesh file									*/
	//	static void				ExportAsSMeshFile(	btSoftBody* psb,
	//												const char* filename);
	/* Create from TetGen .ele, .face, .node files							*/
	//	static btSoftBody*		CreateFromTetGenFile(	btSoftBodyWorldInfo& worldInfo,
	//													const char* ele,
	//													const char* face,
	//													const char* node,
	//													bool bfacelinks,
	//													bool btetralinks,
	//													bool bfacesfromtetras);
	/* Create from TetGen .ele, .face, .node data							*/
	static btSoftBody* CreateFromTetGenData(btSoftBodyWorldInfo& worldInfo,
											const char* ele,
											const char* face,
											const char* node,
											bool bfacelinks,
											bool btetralinks,
											bool bfacesfromtetras);
	static btSoftBody* CreateFromVtkFile(btSoftBodyWorldInfo& worldInfo, const char* vtk_file);

	static void writeObj(const char* file, const btSoftBody* psb);

	static void writeState(const char* file, const btSoftBody* psb);

  //this code cannot be here, dependency on example code are not allowed
	//static std::string loadDeformableState(btAlignedObjectArray<btVector3>& qs, btAlignedObjectArray<btVector3>& vs, const char* filename, CommonFileIOInterface* fileIO);

	static void getBarycentricWeights(const btVector3& a, const btVector3& b, const btVector3& c, const btVector3& d, const btVector3& p, btVector4& bary);

	static void getBarycentricWeights(const btVector3& a, const btVector3& b, const btVector3& c, const btVector3& p, btVector4& bary);

	static void interpolateBarycentricWeights(btSoftBody* psb);

	static void extrapolateBarycentricWeights(btSoftBody* psb);

	static void generateBoundaryFaces(btSoftBody* psb);

	static void duplicateFaces(const char* filename, const btSoftBody* psb);
	/// Sort the list of links to move link calculations that are dependent upon earlier
	/// ones as far as possible away from the calculation of those values
	/// This tends to make adjacent loop iterations not dependent upon one another,
	/// so out-of-order processors can execute instructions from multiple iterations at once
	static void ReoptimizeLinkOrder(btSoftBody* psb);
};

#endif  //BT_SOFT_BODY_HELPERS_H
/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
///btSoftBody implementation by Nathanael Presson

#ifndef _BT_SOFT_BODY_INTERNALS_H
#define _BT_SOFT_BODY_INTERNALS_H

#include "btSoftBody.h"
#include "LinearMath/btQuickprof.h"
#include "LinearMath/btPolarDecomposition.h"
#include "BulletCollision/BroadphaseCollision/btBroadphaseInterface.h"
#include "BulletCollision/CollisionDispatch/btCollisionDispatcher.h"
#include "BulletCollision/CollisionShapes/btConvexInternalShape.h"
#include "BulletCollision/NarrowPhaseCollision/btGjkEpa2.h"
#include "BulletDynamics/Featherstone/btMultiBodyLinkCollider.h"
#include "BulletDynamics/Featherstone/btMultiBodyConstraint.h"
#include <string.h>  //for memset
#include <cmath>
#include "poly34.h"

// Given a multibody link, a contact point and a contact direction, fill in the jacobian data needed to calculate the velocity change given an impulse in the contact direction
static SIMD_FORCE_INLINE void findJacobian(const btMultiBodyLinkCollider* multibodyLinkCol,
										   btMultiBodyJacobianData& jacobianData,
										   const btVector3& contact_point,
										   const btVector3& dir)
{
	const int ndof = multibodyLinkCol->m_multiBody->getNumDofs() + 6;
	jacobianData.m_jacobians.resize(ndof);
	jacobianData.m_deltaVelocitiesUnitImpulse.resize(ndof);
	btScalar* jac = &jacobianData.m_jacobians[0];

	multibodyLinkCol->m_multiBody->fillContactJacobianMultiDof(multibodyLinkCol->m_link, contact_point, dir, jac, jacobianData.scratch_r, jacobianData.scratch_v, jacobianData.scratch_m);
	multibodyLinkCol->m_multiBody->calcAccelerationDeltasMultiDof(&jacobianData.m_jacobians[0], &jacobianData.m_deltaVelocitiesUnitImpulse[0], jacobianData.scratch_r, jacobianData.scratch_v);
}
static SIMD_FORCE_INLINE btVector3 generateUnitOrthogonalVector(const btVector3& u)
{
	btScalar ux = u.getX();
	btScalar uy = u.getY();
	btScalar uz = u.getZ();
	btScalar ax = std::abs(ux);
	btScalar ay = std::abs(uy);
	btScalar az = std::abs(uz);
	btVector3 v;
	if (ax <= ay && ax <= az)
		v = btVector3(0, -uz, uy);
	else if (ay <= ax && ay <= az)
		v = btVector3(-uz, 0, ux);
	else
		v = btVector3(-uy, ux, 0);
	v.normalize();
	return v;
}

static SIMD_FORCE_INLINE bool proximityTest(const btVector3& x1, const btVector3& x2, const btVector3& x3, const btVector3& x4, const btVector3& normal, const btScalar& mrg, btVector3& bary)
{
	btVector3 x43 = x4 - x3;
	if (std::abs(x43.dot(normal)) > mrg)
		return false;
	btVector3 x13 = x1 - x3;
	btVector3 x23 = x2 - x3;
	btScalar a11 = x13.length2();
	btScalar a22 = x23.length2();
	btScalar a12 = x13.dot(x23);
	btScalar b1 = x13.dot(x43);
	btScalar b2 = x23.dot(x43);
	btScalar det = a11 * a22 - a12 * a12;
	if (det < SIMD_EPSILON)
		return false;
	btScalar w1 = (b1 * a22 - b2 * a12) / det;
	btScalar w2 = (b2 * a11 - b1 * a12) / det;
	btScalar w3 = 1 - w1 - w2;
	btScalar delta = 0.01 /*mrg / std::sqrt(0.5 * std::abs(x13.cross(x23).safeNorm()))*/;
	bary = btVector3(w1, w2, w3);
	for (int i = 0; i < 3; ++i)
	{
		if (bary[i] < -delta || bary[i] > 1 + delta)
			return false;
	}
	return true;
}
static const int KDOP_COUNT = 13;
static btVector3 dop[KDOP_COUNT] = {btVector3(1, 0, 0),
									btVector3(0, 1, 0),
									btVector3(0, 0, 1),
									btVector3(1, 1, 0),
									btVector3(1, 0, 1),
									btVector3(0, 1, 1),
									btVector3(1, -1, 0),
									btVector3(1, 0, -1),
									btVector3(0, 1, -1),
									btVector3(1, 1, 1),
									btVector3(1, -1, 1),
									btVector3(1, 1, -1),
									btVector3(1, -1, -1)};

static inline int getSign(const btVector3& n, const btVector3& x)
{
	btScalar d = n.dot(x);
	if (d > SIMD_EPSILON)
		return 1;
	if (d < -SIMD_EPSILON)
		return -1;
	return 0;
}

static SIMD_FORCE_INLINE bool hasSeparatingPlane(const btSoftBody::Face* face, const btSoftBody::Node* node, const btScalar& dt)
{
	btVector3 hex[6] = {face->m_n[0]->m_x - node->m_x,
						face->m_n[1]->m_x - node->m_x,
						face->m_n[2]->m_x - node->m_x,
						face->m_n[0]->m_x + dt * face->m_n[0]->m_v - node->m_x,
						face->m_n[1]->m_x + dt * face->m_n[1]->m_v - node->m_x,
						face->m_n[2]->m_x + dt * face->m_n[2]->m_v - node->m_x};
	btVector3 segment = dt * node->m_v;
	for (int i = 0; i < KDOP_COUNT; ++i)
	{
		int s = getSign(dop[i], segment);
		int j = 0;
		for (; j < 6; ++j)
		{
			if (getSign(dop[i], hex[j]) == s)
				break;
		}
		if (j == 6)
			return true;
	}
	return false;
}

static SIMD_FORCE_INLINE bool nearZero(const btScalar& a)
{
	return (a > -SAFE_EPSILON && a < SAFE_EPSILON);
}
static SIMD_FORCE_INLINE bool sameSign(const btScalar& a, const btScalar& b)
{
	return (nearZero(a) || nearZero(b) || (a > SAFE_EPSILON && b > SAFE_EPSILON) || (a < -SAFE_EPSILON && b < -SAFE_EPSILON));
}
static SIMD_FORCE_INLINE bool diffSign(const btScalar& a, const btScalar& b)
{
	return !sameSign(a, b);
}
inline btScalar evaluateBezier2(const btScalar& p0, const btScalar& p1, const btScalar& p2, const btScalar& t, const btScalar& s)
{
	btScalar s2 = s * s;
	btScalar t2 = t * t;

	return p0 * s2 + p1 * btScalar(2.0) * s * t + p2 * t2;
}
inline btScalar evaluateBezier(const btScalar& p0, const btScalar& p1, const btScalar& p2, const btScalar& p3, const btScalar& t, const btScalar& s)
{
	btScalar s2 = s * s;
	btScalar s3 = s2 * s;
	btScalar t2 = t * t;
	btScalar t3 = t2 * t;

	return p0 * s3 + p1 * btScalar(3.0) * s2 * t + p2 * btScalar(3.0) * s * t2 + p3 * t3;
}
static SIMD_FORCE_INLINE bool getSigns(bool type_c, const btScalar& k0, const btScalar& k1, const btScalar& k2, const btScalar& k3, const btScalar& t0, const btScalar& t1, btScalar& lt0, btScalar& lt1)
{
	if (sameSign(t0, t1))
	{
		lt0 = t0;
		lt1 = t0;
		return true;
	}

	if (type_c || diffSign(k0, k3))
	{
		btScalar ft = evaluateBezier(k0, k1, k2, k3, t0, -t1);
		if (t0 < -0)
			ft = -ft;

		if (sameSign(ft, k0))
		{
			lt0 = t1;
			lt1 = t1;
		}
		else
		{
			lt0 = t0;
			lt1 = t0;
		}
		return true;
	}

	if (!type_c)
	{
		btScalar ft = evaluateBezier(k0, k1, k2, k3, t0, -t1);
		if (t0 < -0)
			ft = -ft;

		if (diffSign(ft, k0))
		{
			lt0 = t0;
			lt1 = t1;
			return true;
		}

		btScalar fk = evaluateBezier2(k1 - k0, k2 - k1, k3 - k2, t0, -t1);

		if (sameSign(fk, k1 - k0))
			lt0 = lt1 = t1;
		else
			lt0 = lt1 = t0;

		return true;
	}
	return false;
}

static SIMD_FORCE_INLINE void getBernsteinCoeff(const btSoftBody::Face* face, const btSoftBody::Node* node, const btScalar& dt, btScalar& k0, btScalar& k1, btScalar& k2, btScalar& k3)
{
	const btVector3& n0 = face->m_n0;
	const btVector3& n1 = face->m_n1;
	btVector3 n_hat = n0 + n1 - face->m_vn;
	btVector3 p0ma0 = node->m_x - face->m_n[0]->m_x;
	btVector3 p1ma1 = node->m_q - face->m_n[0]->m_q;
	k0 = (p0ma0).dot(n0) * 3.0;
	k1 = (p0ma0).dot(n_hat) + (p1ma1).dot(n0);
	k2 = (p1ma1).dot(n_hat) + (p0ma0).dot(n1);
	k3 = (p1ma1).dot(n1) * 3.0;
}

static SIMD_FORCE_INLINE void polyDecomposition(const btScalar& k0, const btScalar& k1, const btScalar& k2, const btScalar& k3, const btScalar& j0, const btScalar& j1, const btScalar& j2, btScalar& u0, btScalar& u1, btScalar& v0, btScalar& v1)
{
	btScalar denom = 4.0 * (j1 - j2) * (j1 - j0) + (j2 - j0) * (j2 - j0);
	u0 = (2.0 * (j1 - j2) * (3.0 * k1 - 2.0 * k0 - k3) - (j0 - j2) * (3.0 * k2 - 2.0 * k3 - k0)) / denom;
	u1 = (2.0 * (j1 - j0) * (3.0 * k2 - 2.0 * k3 - k0) - (j2 - j0) * (3.0 * k1 - 2.0 * k0 - k3)) / denom;
	v0 = k0 - u0 * j0;
	v1 = k3 - u1 * j2;
}

static SIMD_FORCE_INLINE bool rootFindingLemma(const btScalar& k0, const btScalar& k1, const btScalar& k2, const btScalar& k3)
{
	btScalar u0, u1, v0, v1;
	btScalar j0 = 3.0 * (k1 - k0);
	btScalar j1 = 3.0 * (k2 - k1);
	btScalar j2 = 3.0 * (k3 - k2);
	polyDecomposition(k0, k1, k2, k3, j0, j1, j2, u0, u1, v0, v1);
	if (sameSign(v0, v1))
	{
		btScalar Ypa = j0 * (1.0 - v0) * (1.0 - v0) + 2.0 * j1 * v0 * (1.0 - v0) + j2 * v0 * v0;  // Y'(v0)
		if (sameSign(Ypa, j0))
		{
			return (diffSign(k0, v1));
		}
	}
	return diffSign(k0, v0);
}

static SIMD_FORCE_INLINE void getJs(const btScalar& k0, const btScalar& k1, const btScalar& k2, const btScalar& k3, const btSoftBody::Node* a, const btSoftBody::Node* b, const btSoftBody::Node* c, const btSoftBody::Node* p, const btScalar& dt, btScalar& j0, btScalar& j1, btScalar& j2)
{
	const btVector3& a0 = a->m_x;
	const btVector3& b0 = b->m_x;
	const btVector3& c0 = c->m_x;
	const btVector3& va = a->m_v;
	const btVector3& vb = b->m_v;
	const btVector3& vc = c->m_v;
	const btVector3 a1 = a0 + dt * va;
	const btVector3 b1 = b0 + dt * vb;
	const btVector3 c1 = c0 + dt * vc;
	btVector3 n0 = (b0 - a0).cross(c0 - a0);
	btVector3 n1 = (b1 - a1).cross(c1 - a1);
	btVector3 n_hat = n0 + n1 - dt * dt * (vb - va).cross(vc - va);
	const btVector3& p0 = p->m_x;
	const btVector3& vp = p->m_v;
	btVector3 p1 = p0 + dt * vp;
	btVector3 m0 = (b0 - p0).cross(c0 - p0);
	btVector3 m1 = (b1 - p1).cross(c1 - p1);
	btVector3 m_hat = m0 + m1 - dt * dt * (vb - vp).cross(vc - vp);
	btScalar l0 = m0.dot(n0);
	btScalar l1 = 0.25 * (m0.dot(n_hat) + m_hat.dot(n0));
	btScalar l2 = btScalar(1) / btScalar(6) * (m0.dot(n1) + m_hat.dot(n_hat) + m1.dot(n0));
	btScalar l3 = 0.25 * (m_hat.dot(n1) + m1.dot(n_hat));
	btScalar l4 = m1.dot(n1);

	btScalar k1p = 0.25 * k0 + 0.75 * k1;
	btScalar k2p = 0.5 * k1 + 0.5 * k2;
	btScalar k3p = 0.75 * k2 + 0.25 * k3;

	btScalar s0 = (l1 * k0 - l0 * k1p) * 4.0;
	btScalar s1 = (l2 * k0 - l0 * k2p) * 2.0;
	btScalar s2 = (l3 * k0 - l0 * k3p) * btScalar(4) / btScalar(3);
	btScalar s3 = l4 * k0 - l0 * k3;

	j0 = (s1 * k0 - s0 * k1) * 3.0;
	j1 = (s2 * k0 - s0 * k2) * 1.5;
	j2 = (s3 * k0 - s0 * k3);
}

static SIMD_FORCE_INLINE bool signDetermination1Internal(const btScalar& k0, const btScalar& k1, const btScalar& k2, const btScalar& k3, const btScalar& u0, const btScalar& u1, const btScalar& v0, const btScalar& v1)
{
	btScalar Yu0 = k0 * (1.0 - u0) * (1.0 - u0) * (1.0 - u0) + 3.0 * k1 * u0 * (1.0 - u0) * (1.0 - u0) + 3.0 * k2 * u0 * u0 * (1.0 - u0) + k3 * u0 * u0 * u0;  // Y(u0)
	btScalar Yv0 = k0 * (1.0 - v0) * (1.0 - v0) * (1.0 - v0) + 3.0 * k1 * v0 * (1.0 - v0) * (1.0 - v0) + 3.0 * k2 * v0 * v0 * (1.0 - v0) + k3 * v0 * v0 * v0;  // Y(v0)

	btScalar sign_Ytp = (u0 > u1) ? Yu0 : -Yu0;
	btScalar L = sameSign(sign_Ytp, k0) ? u1 : u0;
	sign_Ytp = (v0 > v1) ? Yv0 : -Yv0;
	btScalar K = (sameSign(sign_Ytp, k0)) ? v1 : v0;
	return diffSign(L, K);
}

static SIMD_FORCE_INLINE bool signDetermination2Internal(const btScalar& k0, const btScalar& k1, const btScalar& k2, const btScalar& k3, const btScalar& j0, const btScalar& j1, const btScalar& j2, const btScalar& u0, const btScalar& u1, const btScalar& v0, const btScalar& v1)
{
	btScalar Yu0 = k0 * (1.0 - u0) * (1.0 - u0) * (1.0 - u0) + 3.0 * k1 * u0 * (1.0 - u0) * (1.0 - u0) + 3.0 * k2 * u0 * u0 * (1.0 - u0) + k3 * u0 * u0 * u0;  // Y(u0)
	btScalar sign_Ytp = (u0 > u1) ? Yu0 : -Yu0, L1, L2;
	if (diffSign(sign_Ytp, k0))
	{
		L1 = u0;
		L2 = u1;
	}
	else
	{
		btScalar Yp_u0 = j0 * (1.0 - u0) * (1.0 - u0) + 2.0 * j1 * (1.0 - u0) * u0 + j2 * u0 * u0;
		if (sameSign(Yp_u0, j0))
		{
			L1 = u1;
			L2 = u1;
		}
		else
		{
			L1 = u0;
			L2 = u0;
		}
	}
	btScalar Yv0 = k0 * (1.0 - v0) * (1.0 - v0) * (1.0 - v0) + 3.0 * k1 * v0 * (1.0 - v0) * (1.0 - v0) + 3.0 * k2 * v0 * v0 * (1.0 - v0) + k3 * v0 * v0 * v0;  // Y(uv0)
	sign_Ytp = (v0 > v1) ? Yv0 : -Yv0;
	btScalar K1, K2;
	if (diffSign(sign_Ytp, k0))
	{
		K1 = v0;
		K2 = v1;
	}
	else
	{
		btScalar Yp_v0 = j0 * (1.0 - v0) * (1.0 - v0) + 2.0 * j1 * (1.0 - v0) * v0 + j2 * v0 * v0;
		if (sameSign(Yp_v0, j0))
		{
			K1 = v1;
			K2 = v1;
		}
		else
		{
			K1 = v0;
			K2 = v0;
		}
	}
	return (diffSign(K1, L1) || diffSign(L2, K2));
}

static SIMD_FORCE_INLINE bool signDetermination1(const btScalar& k0, const btScalar& k1, const btScalar& k2, const btScalar& k3, const btSoftBody::Face* face, const btSoftBody::Node* node, const btScalar& dt)
{
	btScalar j0, j1, j2, u0, u1, v0, v1;
	// p1
	getJs(k0, k1, k2, k3, face->m_n[0], face->m_n[1], face->m_n[2], node, dt, j0, j1, j2);
	if (nearZero(j0 + j2 - j1 * 2.0))
	{
		btScalar lt0, lt1;
		getSigns(true, k0, k1, k2, k3, j0, j2, lt0, lt1);
		if (lt0 < -SAFE_EPSILON)
			return false;
	}
	else
	{
		polyDecomposition(k0, k1, k2, k3, j0, j1, j2, u0, u1, v0, v1);
		if (!signDetermination1Internal(k0, k1, k2, k3, u0, u1, v0, v1))
			return false;
	}
	// p2
	getJs(k0, k1, k2, k3, face->m_n[1], face->m_n[2], face->m_n[0], node, dt, j0, j1, j2);
	if (nearZero(j0 + j2 - j1 * 2.0))
	{
		btScalar lt0, lt1;
		getSigns(true, k0, k1, k2, k3, j0, j2, lt0, lt1);
		if (lt0 < -SAFE_EPSILON)
			return false;
	}
	else
	{
		polyDecomposition(k0, k1, k2, k3, j0, j1, j2, u0, u1, v0, v1);
		if (!signDetermination1Internal(k0, k1, k2, k3, u0, u1, v0, v1))
			return false;
	}
	// p3
	getJs(k0, k1, k2, k3, face->m_n[2], face->m_n[0], face->m_n[1], node, dt, j0, j1, j2);
	if (nearZero(j0 + j2 - j1 * 2.0))
	{
		btScalar lt0, lt1;
		getSigns(true, k0, k1, k2, k3, j0, j2, lt0, lt1);
		if (lt0 < -SAFE_EPSILON)
			return false;
	}
	else
	{
		polyDecomposition(k0, k1, k2, k3, j0, j1, j2, u0, u1, v0, v1);
		if (!signDetermination1Internal(k0, k1, k2, k3, u0, u1, v0, v1))
			return false;
	}
	return true;
}

static SIMD_FORCE_INLINE bool signDetermination2(const btScalar& k0, const btScalar& k1, const btScalar& k2, const btScalar& k3, const btSoftBody::Face* face, const btSoftBody::Node* node, const btScalar& dt)
{
	btScalar j0, j1, j2, u0, u1, v0, v1;
	// p1
	getJs(k0, k1, k2, k3, face->m_n[0], face->m_n[1], face->m_n[2], node, dt, j0, j1, j2);
	if (nearZero(j0 + j2 - j1 * 2.0))
	{
		btScalar lt0, lt1;
		bool bt0 = true, bt1 = true;
		getSigns(false, k0, k1, k2, k3, j0, j2, lt0, lt1);
		if (lt0 < -SAFE_EPSILON)
			bt0 = false;
		if (lt1 < -SAFE_EPSILON)
			bt1 = false;
		if (!bt0 && !bt1)
			return false;
	}
	else
	{
		polyDecomposition(k0, k1, k2, k3, j0, j1, j2, u0, u1, v0, v1);
		if (!signDetermination2Internal(k0, k1, k2, k3, j0, j1, j2, u0, u1, v0, v1))
			return false;
	}
	// p2
	getJs(k0, k1, k2, k3, face->m_n[1], face->m_n[2], face->m_n[0], node, dt, j0, j1, j2);
	if (nearZero(j0 + j2 - j1 * 2.0))
	{
		btScalar lt0, lt1;
		bool bt0 = true, bt1 = true;
		getSigns(false, k0, k1, k2, k3, j0, j2, lt0, lt1);
		if (lt0 < -SAFE_EPSILON)
			bt0 = false;
		if (lt1 < -SAFE_EPSILON)
			bt1 = false;
		if (!bt0 && !bt1)
			return false;
	}
	else
	{
		polyDecomposition(k0, k1, k2, k3, j0, j1, j2, u0, u1, v0, v1);
		if (!signDetermination2Internal(k0, k1, k2, k3, j0, j1, j2, u0, u1, v0, v1))
			return false;
	}
	// p3
	getJs(k0, k1, k2, k3, face->m_n[2], face->m_n[0], face->m_n[1], node, dt, j0, j1, j2);
	if (nearZero(j0 + j2 - j1 * 2.0))
	{
		btScalar lt0, lt1;
		bool bt0 = true, bt1 = true;
		getSigns(false, k0, k1, k2, k3, j0, j2, lt0, lt1);
		if (lt0 < -SAFE_EPSILON)
			bt0 = false;
		if (lt1 < -SAFE_EPSILON)
			bt1 = false;
		if (!bt0 && !bt1)
			return false;
	}
	else
	{
		polyDecomposition(k0, k1, k2, k3, j0, j1, j2, u0, u1, v0, v1);
		if (!signDetermination2Internal(k0, k1, k2, k3, j0, j1, j2, u0, u1, v0, v1))
			return false;
	}
	return true;
}

static SIMD_FORCE_INLINE bool coplanarAndInsideTest(const btScalar& k0, const btScalar& k1, const btScalar& k2, const btScalar& k3, const btSoftBody::Face* face, const btSoftBody::Node* node, const btScalar& dt)
{
	// Coplanar test
	if (diffSign(k1 - k0, k3 - k2))
	{
		// Case b:
		if (sameSign(k0, k3) && !rootFindingLemma(k0, k1, k2, k3))
			return false;
		// inside test
		return signDetermination2(k0, k1, k2, k3, face, node, dt);
	}
	else
	{
		// Case c:
		if (sameSign(k0, k3))
			return false;
		// inside test
		return signDetermination1(k0, k1, k2, k3, face, node, dt);
	}
	return false;
}
static SIMD_FORCE_INLINE bool conservativeCulling(const btScalar& k0, const btScalar& k1, const btScalar& k2, const btScalar& k3, const btScalar& mrg)
{
	if (k0 > mrg && k1 > mrg && k2 > mrg && k3 > mrg)
		return true;
	if (k0 < -mrg && k1 < -mrg && k2 < -mrg && k3 < -mrg)
		return true;
	return false;
}

static SIMD_FORCE_INLINE bool bernsteinVFTest(const btScalar& k0, const btScalar& k1, const btScalar& k2, const btScalar& k3, const btScalar& mrg, const btSoftBody::Face* face, const btSoftBody::Node* node, const btScalar& dt)
{
	if (conservativeCulling(k0, k1, k2, k3, mrg))
		return false;
	return coplanarAndInsideTest(k0, k1, k2, k3, face, node, dt);
}

static SIMD_FORCE_INLINE void deCasteljau(const btScalar& k0, const btScalar& k1, const btScalar& k2, const btScalar& k3, const btScalar& t0, btScalar& k10, btScalar& k20, btScalar& k30, btScalar& k21, btScalar& k12)
{
	k10 = k0 * (1.0 - t0) + k1 * t0;
	btScalar k11 = k1 * (1.0 - t0) + k2 * t0;
	k12 = k2 * (1.0 - t0) + k3 * t0;
	k20 = k10 * (1.0 - t0) + k11 * t0;
	k21 = k11 * (1.0 - t0) + k12 * t0;
	k30 = k20 * (1.0 - t0) + k21 * t0;
}
static SIMD_FORCE_INLINE bool bernsteinVFTest(const btSoftBody::Face* face, const btSoftBody::Node* node, const btScalar& dt, const btScalar& mrg)
{
	btScalar k0, k1, k2, k3;
	getBernsteinCoeff(face, node, dt, k0, k1, k2, k3);
	if (conservativeCulling(k0, k1, k2, k3, mrg))
		return false;
	return true;
	if (diffSign(k2 - 2.0 * k1 + k0, k3 - 2.0 * k2 + k1))
	{
		btScalar k10, k20, k30, k21, k12;
		btScalar t0 = (k2 - 2.0 * k1 + k0) / (k0 - 3.0 * k1 + 3.0 * k2 - k3);
		deCasteljau(k0, k1, k2, k3, t0, k10, k20, k30, k21, k12);
		return bernsteinVFTest(k0, k10, k20, k30, mrg, face, node, dt) || bernsteinVFTest(k30, k21, k12, k3, mrg, face, node, dt);
	}
	return coplanarAndInsideTest(k0, k1, k2, k3, face, node, dt);
}

static SIMD_FORCE_INLINE bool continuousCollisionDetection(const btSoftBody::Face* face, const btSoftBody::Node* node, const btScalar& dt, const btScalar& mrg, btVector3& bary)
{
	if (hasSeparatingPlane(face, node, dt))
		return false;
	btVector3 x21 = face->m_n[1]->m_x - face->m_n[0]->m_x;
	btVector3 x31 = face->m_n[2]->m_x - face->m_n[0]->m_x;
	btVector3 x41 = node->m_x - face->m_n[0]->m_x;
	btVector3 v21 = face->m_n[1]->m_v - face->m_n[0]->m_v;
	btVector3 v31 = face->m_n[2]->m_v - face->m_n[0]->m_v;
	btVector3 v41 = node->m_v - face->m_n[0]->m_v;
	btVector3 a = x21.cross(x31);
	btVector3 b = x21.cross(v31) + v21.cross(x31);
	btVector3 c = v21.cross(v31);
	btVector3 d = x41;
	btVector3 e = v41;
	btScalar a0 = a.dot(d);
	btScalar a1 = a.dot(e) + b.dot(d);
	btScalar a2 = c.dot(d) + b.dot(e);
	btScalar a3 = c.dot(e);
	btScalar eps = SAFE_EPSILON;
	int num_roots = 0;
	btScalar roots[3];
	if (std::abs(a3) < eps)
	{
		// cubic term is zero
		if (std::abs(a2) < eps)
		{
			if (std::abs(a1) < eps)
			{
				if (std::abs(a0) < eps)
				{
					num_roots = 2;
					roots[0] = 0;
					roots[1] = dt;
				}
			}
			else
			{
				num_roots = 1;
				roots[0] = -a0 / a1;
			}
		}
		else
		{
			num_roots = SolveP2(roots, a1 / a2, a0 / a2);
		}
	}
	else
	{
		num_roots = SolveP3(roots, a2 / a3, a1 / a3, a0 / a3);
	}
	//    std::sort(roots, roots+num_roots);
	if (num_roots > 1)
	{
		if (roots[0] > roots[1])
			btSwap(roots[0], roots[1]);
	}
	if (num_roots > 2)
	{
		if (roots[0] > roots[2])
			btSwap(roots[0], roots[2]);
		if (roots[1] > roots[2])
			btSwap(roots[1], roots[2]);
	}
	for (int r = 0; r < num_roots; ++r)
	{
		double root = roots[r];
		if (root <= 0)
			continue;
		if (root > dt + SIMD_EPSILON)
			return false;
		btVector3 x1 = face->m_n[0]->m_x + root * face->m_n[0]->m_v;
		btVector3 x2 = face->m_n[1]->m_x + root * face->m_n[1]->m_v;
		btVector3 x3 = face->m_n[2]->m_x + root * face->m_n[2]->m_v;
		btVector3 x4 = node->m_x + root * node->m_v;
		btVector3 normal = (x2 - x1).cross(x3 - x1);
		normal.safeNormalize();
		if (proximityTest(x1, x2, x3, x4, normal, mrg, bary))
			return true;
	}
	return false;
}
static SIMD_FORCE_INLINE bool bernsteinCCD(const btSoftBody::Face* face, const btSoftBody::Node* node, const btScalar& dt, const btScalar& mrg, btVector3& bary)
{
	if (!bernsteinVFTest(face, node, dt, mrg))
		return false;
	if (!continuousCollisionDetection(face, node, dt, 1e-6, bary))
		return false;
	return true;
}

//
// btSymMatrix
//
template <typename T>
struct btSymMatrix
{
	btSymMatrix() : dim(0) {}
	btSymMatrix(int n, const T& init = T()) { resize(n, init); }
	void resize(int n, const T& init = T())
	{
		dim = n;
		store.resize((n * (n + 1)) / 2, init);
	}
	int index(int c, int r) const
	{
		if (c > r) btSwap(c, r);
		btAssert(r < dim);
		return ((r * (r + 1)) / 2 + c);
	}
	T& operator()(int c, int r) { return (store[index(c, r)]); }
	const T& operator()(int c, int r) const { return (store[index(c, r)]); }
	btAlignedObjectArray<T> store;
	int dim;
};

//
// btSoftBodyCollisionShape
//
class btSoftBodyCollisionShape : public btConcaveShape
{
public:
	btSoftBody* m_body;

	btSoftBodyCollisionShape(btSoftBody* backptr)
	{
		m_shapeType = SOFTBODY_SHAPE_PROXYTYPE;
		m_body = backptr;
	}

	virtual ~btSoftBodyCollisionShape()
	{
	}

	void processAllTriangles(btTriangleCallback* /*callback*/, const btVector3& /*aabbMin*/, const btVector3& /*aabbMax*/) const
	{
		//not yet
		btAssert(0);
	}

	///getAabb returns the axis aligned bounding box in the coordinate frame of the given transform t.
	virtual void getAabb(const btTransform& t, btVector3& aabbMin, btVector3& aabbMax) const
	{
		/* t is usually identity, except when colliding against btCompoundShape. See Issue 512 */
		const btVector3 mins = m_body->m_bounds[0];
		const btVector3 maxs = m_body->m_bounds[1];
		const btVector3 crns[] = {t * btVector3(mins.x(), mins.y(), mins.z()),
								  t * btVector3(maxs.x(), mins.y(), mins.z()),
								  t * btVector3(maxs.x(), maxs.y(), mins.z()),
								  t * btVector3(mins.x(), maxs.y(), mins.z()),
								  t * btVector3(mins.x(), mins.y(), maxs.z()),
								  t * btVector3(maxs.x(), mins.y(), maxs.z()),
								  t * btVector3(maxs.x(), maxs.y(), maxs.z()),
								  t * btVector3(mins.x(), maxs.y(), maxs.z())};
		aabbMin = aabbMax = crns[0];
		for (int i = 1; i < 8; ++i)
		{
			aabbMin.setMin(crns[i]);
			aabbMax.setMax(crns[i]);
		}
	}

	virtual void setLocalScaling(const btVector3& /*scaling*/)
	{
		///na
	}
	virtual const btVector3& getLocalScaling() const
	{
		static const btVector3 dummy(1, 1, 1);
		return dummy;
	}
	virtual void calculateLocalInertia(btScalar /*mass*/, btVector3& /*inertia*/) const
	{
		///not yet
		btAssert(0);
	}
	virtual const char* getName() const
	{
		return "SoftBody";
	}
};

//
// btSoftClusterCollisionShape
//
class btSoftClusterCollisionShape : public btConvexInternalShape
{
public:
	const btSoftBody::Cluster* m_cluster;

	btSoftClusterCollisionShape(const btSoftBody::Cluster* cluster) : m_cluster(cluster) { setMargin(0); }

	virtual btVector3 localGetSupportingVertex(const btVector3& vec) const
	{
		btSoftBody::Node* const* n = &m_cluster->m_nodes[0];
		btScalar d = btDot(vec, n[0]->m_x);
		int j = 0;
		for (int i = 1, ni = m_cluster->m_nodes.size(); i < ni; ++i)
		{
			const btScalar k = btDot(vec, n[i]->m_x);
			if (k > d)
			{
				d = k;
				j = i;
			}
		}
		return (n[j]->m_x);
	}
	virtual btVector3 localGetSupportingVertexWithoutMargin(const btVector3& vec) const
	{
		return (localGetSupportingVertex(vec));
	}
	//notice that the vectors should be unit length
	virtual void batchedUnitVectorGetSupportingVertexWithoutMargin(const btVector3* vectors, btVector3* supportVerticesOut, int numVectors) const
	{
	}

	virtual void calculateLocalInertia(btScalar mass, btVector3& inertia) const
	{
	}

	virtual void getAabb(const btTransform& t, btVector3& aabbMin, btVector3& aabbMax) const
	{
	}

	virtual int getShapeType() const { return SOFTBODY_SHAPE_PROXYTYPE; }

	//debugging
	virtual const char* getName() const { return "SOFTCLUSTER"; }

	virtual void setMargin(btScalar margin)
	{
		btConvexInternalShape::setMargin(margin);
	}
	virtual btScalar getMargin() const
	{
		return btConvexInternalShape::getMargin();
	}
};

//
// Inline's
//

//
template <typename T>
static inline void ZeroInitialize(T& value)
{
	memset(&value, 0, sizeof(T));
}
//
template <typename T>
static inline bool CompLess(const T& a, const T& b)
{
	return (a < b);
}
//
template <typename T>
static inline bool CompGreater(const T& a, const T& b)
{
	return (a > b);
}
//
template <typename T>
static inline T Lerp(const T& a, const T& b, btScalar t)
{
	return (a + (b - a) * t);
}
//
template <typename T>
static inline T InvLerp(const T& a, const T& b, btScalar t)
{
	return ((b + a * t - b * t) / (a * b));
}
//
static inline btMatrix3x3 Lerp(const btMatrix3x3& a,
							   const btMatrix3x3& b,
							   btScalar t)
{
	btMatrix3x3 r;
	r[0] = Lerp(a[0], b[0], t);
	r[1] = Lerp(a[1], b[1], t);
	r[2] = Lerp(a[2], b[2], t);
	return (r);
}
//
static inline btVector3 Clamp(const btVector3& v, btScalar maxlength)
{
	const btScalar sql = v.length2();
	if (sql > (maxlength * maxlength))
		return ((v * maxlength) / btSqrt(sql));
	else
		return (v);
}
//
template <typename T>
static inline T Clamp(const T& x, const T& l, const T& h)
{
	return (x < l ? l : x > h ? h
							  : x);
}
//
template <typename T>
static inline T Sq(const T& x)
{
	return (x * x);
}
//
template <typename T>
static inline T Cube(const T& x)
{
	return (x * x * x);
}
//
template <typename T>
static inline T Sign(const T& x)
{
	return ((T)(x < 0 ? -1 : +1));
}
//
template <typename T>
static inline bool SameSign(const T& x, const T& y)
{
	return ((x * y) > 0);
}
//
static inline btScalar ClusterMetric(const btVector3& x, const btVector3& y)
{
	const btVector3 d = x - y;
	return (btFabs(d[0]) + btFabs(d[1]) + btFabs(d[2]));
}
//
static inline btMatrix3x3 ScaleAlongAxis(const btVector3& a, btScalar s)
{
	const btScalar xx = a.x() * a.x();
	const btScalar yy = a.y() * a.y();
	const btScalar zz = a.z() * a.z();
	const btScalar xy = a.x() * a.y();
	const btScalar yz = a.y() * a.z();
	const btScalar zx = a.z() * a.x();
	btMatrix3x3 m;
	m[0] = btVector3(1 - xx + xx * s, xy * s - xy, zx * s - zx);
	m[1] = btVector3(xy * s - xy, 1 - yy + yy * s, yz * s - yz);
	m[2] = btVector3(zx * s - zx, yz * s - yz, 1 - zz + zz * s);
	return (m);
}
//
static inline btMatrix3x3 Cross(const btVector3& v)
{
	btMatrix3x3 m;
	m[0] = btVector3(0, -v.z(), +v.y());
	m[1] = btVector3(+v.z(), 0, -v.x());
	m[2] = btVector3(-v.y(), +v.x(), 0);
	return (m);
}
//
static inline btMatrix3x3 Diagonal(btScalar x)
{
	btMatrix3x3 m;
	m[0] = btVector3(x, 0, 0);
	m[1] = btVector3(0, x, 0);
	m[2] = btVector3(0, 0, x);
	return (m);
}

static inline btMatrix3x3 Diagonal(const btVector3& v)
{
	btMatrix3x3 m;
	m[0] = btVector3(v.getX(), 0, 0);
	m[1] = btVector3(0, v.getY(), 0);
	m[2] = btVector3(0, 0, v.getZ());
	return (m);
}

static inline btScalar Dot(const btScalar* a, const btScalar* b, int ndof)
{
	btScalar result = 0;
	for (int i = 0; i < ndof; ++i)
		result += a[i] * b[i];
	return result;
}

static inline btMatrix3x3 OuterProduct(const btScalar* v1, const btScalar* v2, const btScalar* v3,
									   const btScalar* u1, const btScalar* u2, const btScalar* u3, int ndof)
{
	btMatrix3x3 m;
	btScalar a11 = Dot(v1, u1, ndof);
	btScalar a12 = Dot(v1, u2, ndof);
	btScalar a13 = Dot(v1, u3, ndof);

	btScalar a21 = Dot(v2, u1, ndof);
	btScalar a22 = Dot(v2, u2, ndof);
	btScalar a23 = Dot(v2, u3, ndof);

	btScalar a31 = Dot(v3, u1, ndof);
	btScalar a32 = Dot(v3, u2, ndof);
	btScalar a33 = Dot(v3, u3, ndof);
	m[0] = btVector3(a11, a12, a13);
	m[1] = btVector3(a21, a22, a23);
	m[2] = btVector3(a31, a32, a33);
	return (m);
}

static inline btMatrix3x3 OuterProduct(const btVector3& v1, const btVector3& v2)
{
	btMatrix3x3 m;
	btScalar a11 = v1[0] * v2[0];
	btScalar a12 = v1[0] * v2[1];
	btScalar a13 = v1[0] * v2[2];

	btScalar a21 = v1[1] * v2[0];
	btScalar a22 = v1[1] * v2[1];
	btScalar a23 = v1[1] * v2[2];

	btScalar a31 = v1[2] * v2[0];
	btScalar a32 = v1[2] * v2[1];
	btScalar a33 = v1[2] * v2[2];
	m[0] = btVector3(a11, a12, a13);
	m[1] = btVector3(a21, a22, a23);
	m[2] = btVector3(a31, a32, a33);
	return (m);
}

//
static inline btMatrix3x3 Add(const btMatrix3x3& a,
							  const btMatrix3x3& b)
{
	btMatrix3x3 r;
	for (int i = 0; i < 3; ++i) r[i] = a[i] + b[i];
	return (r);
}
//
static inline btMatrix3x3 Sub(const btMatrix3x3& a,
							  const btMatrix3x3& b)
{
	btMatrix3x3 r;
	for (int i = 0; i < 3; ++i) r[i] = a[i] - b[i];
	return (r);
}
//
static inline btMatrix3x3 Mul(const btMatrix3x3& a,
							  btScalar b)
{
	btMatrix3x3 r;
	for (int i = 0; i < 3; ++i) r[i] = a[i] * b;
	return (r);
}
//
static inline void Orthogonalize(btMatrix3x3& m)
{
	m[2] = btCross(m[0], m[1]).normalized();
	m[1] = btCross(m[2], m[0]).normalized();
	m[0] = btCross(m[1], m[2]).normalized();
}
//
static inline btMatrix3x3 MassMatrix(btScalar im, const btMatrix3x3& iwi, const btVector3& r)
{
	const btMatrix3x3 cr = Cross(r);
	return (Sub(Diagonal(im), cr * iwi * cr));
}

//
static inline btMatrix3x3 ImpulseMatrix(btScalar dt,
										btScalar ima,
										btScalar imb,
										const btMatrix3x3& iwi,
										const btVector3& r)
{
	return (Diagonal(1 / dt) * Add(Diagonal(ima), MassMatrix(imb, iwi, r)).inverse());
}

//
static inline btMatrix3x3 ImpulseMatrix(btScalar dt,
										const btMatrix3x3& effective_mass_inv,
										btScalar imb,
										const btMatrix3x3& iwi,
										const btVector3& r)
{
	return (Diagonal(1 / dt) * Add(effective_mass_inv, MassMatrix(imb, iwi, r)).inverse());
	//    btMatrix3x3 iimb = MassMatrix(imb, iwi, r);
	//    if (iimb.determinant() == 0)
	//        return effective_mass_inv.inverse();
	//    return effective_mass_inv.inverse() *  Add(effective_mass_inv.inverse(), iimb.inverse()).inverse() * iimb.inverse();
}

//
static inline btMatrix3x3 ImpulseMatrix(btScalar ima, const btMatrix3x3& iia, const btVector3& ra,
										btScalar imb, const btMatrix3x3& iib, const btVector3& rb)
{
	return (Add(MassMatrix(ima, iia, ra), MassMatrix(imb, iib, rb)).inverse());
}

//
static inline btMatrix3x3 AngularImpulseMatrix(const btMatrix3x3& iia,
											   const btMatrix3x3& iib)
{
	return (Add(iia, iib).inverse());
}

//
static inline btVector3 ProjectOnAxis(const btVector3& v,
									  const btVector3& a)
{
	return (a * btDot(v, a));
}
//
static inline btVector3 ProjectOnPlane(const btVector3& v,
									   const btVector3& a)
{
	return (v - ProjectOnAxis(v, a));
}

//
static inline void ProjectOrigin(const btVector3& a,
								 const btVector3& b,
								 btVector3& prj,
								 btScalar& sqd)
{
	const btVector3 d = b - a;
	const btScalar m2 = d.length2();
	if (m2 > SIMD_EPSILON)
	{
		const btScalar t = Clamp<btScalar>(-btDot(a, d) / m2, 0, 1);
		const btVector3 p = a + d * t;
		const btScalar l2 = p.length2();
		if (l2 < sqd)
		{
			prj = p;
			sqd = l2;
		}
	}
}
//
static inline void ProjectOrigin(const btVector3& a,
								 const btVector3& b,
								 const btVector3& c,
								 btVector3& prj,
								 btScalar& sqd)
{
	const btVector3& q = btCross(b - a, c - a);
	const btScalar m2 = q.length2();
	if (m2 > SIMD_EPSILON)
	{
		const btVector3 n = q / btSqrt(m2);
		const btScalar k = btDot(a, n);
		const btScalar k2 = k * k;
		if (k2 < sqd)
		{
			const btVector3 p = n * k;
			if ((btDot(btCross(a - p, b - p), q) > 0) &&
				(btDot(btCross(b - p, c - p), q) > 0) &&
				(btDot(btCross(c - p, a - p), q) > 0))
			{
				prj = p;
				sqd = k2;
			}
			else
			{
				ProjectOrigin(a, b, prj, sqd);
				ProjectOrigin(b, c, prj, sqd);
				ProjectOrigin(c, a, prj, sqd);
			}
		}
	}
}

//
static inline bool rayIntersectsTriangle(const btVector3& origin, const btVector3& dir, const btVector3& v0, const btVector3& v1, const btVector3& v2, btScalar& t)
{
	btScalar a, f, u, v;

	btVector3 e1 = v1 - v0;
	btVector3 e2 = v2 - v0;
	btVector3 h = dir.cross(e2);
	a = e1.dot(h);

	if (a > -0.00001 && a < 0.00001)
		return (false);

	f = btScalar(1) / a;
	btVector3 s = origin - v0;
	u = f * s.dot(h);

	if (u < 0.0 || u > 1.0)
		return (false);

	btVector3 q = s.cross(e1);
	v = f * dir.dot(q);
	if (v < 0.0 || u + v > 1.0)
		return (false);
	// at this stage we can compute t to find out where
	// the intersection point is on the line
	t = f * e2.dot(q);
	if (t > 0)  // ray intersection
		return (true);
	else  // this means that there is a line intersection
		// but not a ray intersection
		return (false);
}

static inline bool lineIntersectsTriangle(const btVector3& rayStart, const btVector3& rayEnd, const btVector3& p1, const btVector3& p2, const btVector3& p3, btVector3& sect, btVector3& normal)
{
	btVector3 dir = rayEnd - rayStart;
	btScalar dir_norm = dir.norm();
	if (dir_norm < SIMD_EPSILON)
		return false;
	dir.normalize();
	btScalar t;
	bool ret = rayIntersectsTriangle(rayStart, dir, p1, p2, p3, t);

	if (ret)
	{
		if (t <= dir_norm)
		{
			sect = rayStart + dir * t;
		}
		else
		{
			ret = false;
		}
	}

	if (ret)
	{
		btVector3 n = (p3 - p1).cross(p2 - p1);
		n.safeNormalize();
		if (n.dot(dir) < 0)
			normal = n;
		else
			normal = -n;
	}
	return ret;
}

//
template <typename T>
static inline T BaryEval(const T& a,
						 const T& b,
						 const T& c,
						 const btVector3& coord)
{
	return (a * coord.x() + b * coord.y() + c * coord.z());
}

template <typename T>
static inline T BaryEval(const T& a,
						 const T& b,
						 const T& c,
						 const T& d,
						 const btVector4& coord)
{
	return (a * coord.x() + b * coord.y() + c * coord.z() + d * coord.w());
}

//
static inline btVector3 BaryCoord(const btVector3& a,
								  const btVector3& b,
								  const btVector3& c,
								  const btVector3& p)
{
	const btScalar w[] = {btCross(a - p, b - p).length(),
						  btCross(b - p, c - p).length(),
						  btCross(c - p, a - p).length()};
	const btScalar isum = 1 / (w[0] + w[1] + w[2]);
	return (btVector3(w[1] * isum, w[2] * isum, w[0] * isum));
}

static inline btVector4 BaryCoord(const btVector3& a,
								  const btVector3& b,
								  const btVector3& c,
								  const btVector3& d,
								  const btVector3& p)
{
	// signed volume of a tetrahedron (factor 6 cancels out later)
	auto vol = [](const btVector3& p0,
				  const btVector3& p1,
				  const btVector3& p2,
				  const btVector3& p3) -> btScalar
	{
		return (p1 - p0).dot((p2 - p0).cross(p3 - p0));
	};

	const btScalar V = vol(a, b, c, d);  // total volume
	btAssert(btFabs(V) > SIMD_EPSILON);  // degenerate check

	const btScalar invV = 1.0f / V;
	const btScalar l0 = vol(p, b, c, d) * invV;  // opposite vertex a
	const btScalar l1 = vol(a, p, c, d) * invV;  // opposite vertex b
	const btScalar l2 = vol(a, b, p, d) * invV;  // opposite vertex c
	const btScalar l3 = vol(a, b, c, p) * invV;  // opposite vertex d

	return btVector4(l0, l1, l2, l3);  // ?0+?1+?2+?3 == 1
}

//
inline static btScalar ImplicitSolve(btSoftBody::ImplicitFn* fn,
									 const btVector3& a,
									 const btVector3& b,
									 const btScalar accuracy,
									 const int maxiterations = 256)
{
	btScalar span[2] = {0, 1};
	btScalar values[2] = {fn->Eval(a), fn->Eval(b)};
	if (values[0] > values[1])
	{
		btSwap(span[0], span[1]);
		btSwap(values[0], values[1]);
	}
	if (values[0] > -accuracy) return (-1);
	if (values[1] < +accuracy) return (-1);
	for (int i = 0; i < maxiterations; ++i)
	{
		const btScalar t = Lerp(span[0], span[1], values[0] / (values[0] - values[1]));
		const btScalar v = fn->Eval(Lerp(a, b, t));
		if ((t <= 0) || (t >= 1)) break;
		if (btFabs(v) < accuracy) return (t);
		if (v < 0)
		{
			span[0] = t;
			values[0] = v;
		}
		else
		{
			span[1] = t;
			values[1] = v;
		}
	}
	return (-1);
}

inline static void EvaluateMedium(const btSoftBodyWorldInfo* wfi,
								  const btVector3& x,
								  btSoftBody::sMedium& medium)
{
	medium.m_velocity = btVector3(0, 0, 0);
	medium.m_pressure = 0;
	medium.m_density = wfi->air_density;
	if (wfi->water_density > 0)
	{
		const btScalar depth = -(btDot(x, wfi->water_normal) + wfi->water_offset);
		if (depth > 0)
		{
			medium.m_density = wfi->water_density;
			medium.m_pressure = depth * wfi->water_density * wfi->m_gravity.length();
		}
	}
}

//
static inline btVector3 NormalizeAny(const btVector3& v)
{
	const btScalar l = v.length();
	if (l > SIMD_EPSILON)
		return (v / l);
	else
		return (btVector3(0, 0, 0));
}

//
static inline btDbvtVolume VolumeOf(const btSoftBody::Face& f,
									btScalar margin)
{
	const btVector3* pts[] = {&f.m_n[0]->m_x,
							  &f.m_n[1]->m_x,
							  &f.m_n[2]->m_x};
	btDbvtVolume vol = btDbvtVolume::FromPoints(pts, 3);
	vol.Expand(btVector3(margin, margin, margin));
	return (vol);
}

//
static inline btVector3 CenterOf(const btSoftBody::Face& f)
{
	return ((f.m_n[0]->m_x + f.m_n[1]->m_x + f.m_n[2]->m_x) / 3);
}

// WARNING: The return value has to be divided by 2.0 to obtain true area value. The division is probably omitted as an optimization.
static inline btScalar AreaOf(const btVector3& x0,
							  const btVector3& x1,
							  const btVector3& x2)
{
	const btVector3 a = x1 - x0;
	const btVector3 b = x2 - x0;
	const btVector3 cr = btCross(a, b);
	const btScalar area = cr.length();
	return (area);
}

//
static inline btScalar VolumeOf(const btVector3& x0,
								const btVector3& x1,
								const btVector3& x2,
								const btVector3& x3)
{
	const btVector3 a = x1 - x0;
	const btVector3 b = x2 - x0;
	const btVector3 c = x3 - x0;
	return (btDot(a, btCross(b, c)));
}

//

//
static inline void ApplyClampedForce(btSoftBody::Node& n,
									 const btVector3& f,
									 btScalar dt)
{
	const btScalar dtim = dt * n.m_im;
	if ((f * dtim).length2() > n.m_v.length2())
	{ /* Clamp	*/
		n.m_f -= ProjectOnAxis(n.m_v, f.normalized()) / dtim;
	}
	else
	{ /* Apply	*/
		n.m_f += f;
	}
}

//
static inline int MatchEdge(const btSoftBody::Node* a,
							const btSoftBody::Node* b,
							const btSoftBody::Node* ma,
							const btSoftBody::Node* mb)
{
	if ((a == ma) && (b == mb)) return (0);
	if ((a == mb) && (b == ma)) return (1);
	return (-1);
}

//
// btEigen : Extract eigen system,
// straitforward implementation of http://math.fullerton.edu/mathews/n2003/JacobiMethodMod.html
// outputs are NOT sorted.
//
struct btEigen
{
	static int system(btMatrix3x3& a, btMatrix3x3* vectors, btVector3* values = 0)
	{
		static const int maxiterations = 16;
		static const btScalar accuracy = (btScalar)0.0001;
		btMatrix3x3& v = *vectors;
		int iterations = 0;
		vectors->setIdentity();
		do
		{
			int p = 0, q = 1;
			if (btFabs(a[p][q]) < btFabs(a[0][2]))
			{
				p = 0;
				q = 2;
			}
			if (btFabs(a[p][q]) < btFabs(a[1][2]))
			{
				p = 1;
				q = 2;
			}
			if (btFabs(a[p][q]) > accuracy)
			{
				const btScalar w = (a[q][q] - a[p][p]) / (2 * a[p][q]);
				const btScalar z = btFabs(w);
				const btScalar t = w / (z * (btSqrt(1 + w * w) + z));
				if (t == t) /* [WARNING] let hope that one does not get thrown aways by some compilers... */
				{
					const btScalar c = 1 / btSqrt(t * t + 1);
					const btScalar s = c * t;
					mulPQ(a, c, s, p, q);
					mulTPQ(a, c, s, p, q);
					mulPQ(v, c, s, p, q);
				}
				else
					break;
			}
			else
				break;
		} while ((++iterations) < maxiterations);
		if (values)
		{
			*values = btVector3(a[0][0], a[1][1], a[2][2]);
		}
		return (iterations);
	}

private:
	static inline void mulTPQ(btMatrix3x3& a, btScalar c, btScalar s, int p, int q)
	{
		const btScalar m[2][3] = {{a[p][0], a[p][1], a[p][2]},
								  {a[q][0], a[q][1], a[q][2]}};
		int i;

		for (i = 0; i < 3; ++i) a[p][i] = c * m[0][i] - s * m[1][i];
		for (i = 0; i < 3; ++i) a[q][i] = c * m[1][i] + s * m[0][i];
	}
	static inline void mulPQ(btMatrix3x3& a, btScalar c, btScalar s, int p, int q)
	{
		const btScalar m[2][3] = {{a[0][p], a[1][p], a[2][p]},
								  {a[0][q], a[1][q], a[2][q]}};
		int i;

		for (i = 0; i < 3; ++i) a[i][p] = c * m[0][i] - s * m[1][i];
		for (i = 0; i < 3; ++i) a[i][q] = c * m[1][i] + s * m[0][i];
	}
};

//
// Polar decomposition,
// "Computing the Polar Decomposition with Applications", Nicholas J. Higham, 1986.
//
static inline int PolarDecompose(const btMatrix3x3& m, btMatrix3x3& q, btMatrix3x3& s)
{
	static const btPolarDecomposition polar;
	return polar.decompose(m, q, s);
}

//
// btSoftColliders
//
struct btSoftColliders
{
	//
	// ClusterBase
	//
	struct ClusterBase : btDbvt::ICollide
	{
		btScalar erp;
		btScalar idt;
		btScalar m_margin;
		btScalar friction;
		btScalar threshold;
		ClusterBase()
		{
			erp = (btScalar)1;
			idt = 0;
			m_margin = 0;
			friction = 0;
			threshold = (btScalar)0;
		}
		bool SolveContact(const btGjkEpaSolver2::sResults& res,
						  btSoftBody::Body ba, const btSoftBody::Body bb,
						  btSoftBody::CJoint& joint)
		{
			if (res.distance < m_margin)
			{
				btVector3 norm = res.normal;
				norm.normalize();  //is it necessary?

				const btVector3 ra = res.witnesses[0] - ba.xform().getOrigin();
				const btVector3 rb = res.witnesses[1] - bb.xform().getOrigin();
				const btVector3 va = ba.velocity(ra);
				const btVector3 vb = bb.velocity(rb);
				const btVector3 vrel = va - vb;
				const btScalar rvac = btDot(vrel, norm);
				btScalar depth = res.distance - m_margin;

				//				printf("depth=%f\n",depth);
				const btVector3 iv = norm * rvac;
				const btVector3 fv = vrel - iv;
				joint.m_bodies[0] = ba;
				joint.m_bodies[1] = bb;
				joint.m_refs[0] = ra * ba.xform().getBasis();
				joint.m_refs[1] = rb * bb.xform().getBasis();
				joint.m_rpos[0] = ra;
				joint.m_rpos[1] = rb;
				joint.m_cfm = 1;
				joint.m_erp = 1;
				joint.m_life = 0;
				joint.m_maxlife = 0;
				joint.m_split = 1;

				joint.m_drift = depth * norm;

				joint.m_normal = norm;
				//				printf("normal=%f,%f,%f\n",res.normal.getX(),res.normal.getY(),res.normal.getZ());
				joint.m_delete = false;
				joint.m_friction = fv.length2() < (rvac * friction * rvac * friction) ? 1 : friction;
				joint.m_massmatrix = ImpulseMatrix(ba.invMass(), ba.invWorldInertia(), joint.m_rpos[0],
												   bb.invMass(), bb.invWorldInertia(), joint.m_rpos[1]);

				return (true);
			}
			return (false);
		}
	};
	//
	// CollideCL_RS
	//
	struct CollideCL_RS : ClusterBase
	{
		btSoftBody* psb;
		const btCollisionObjectWrapper* m_colObjWrap;

		void Process(const btDbvtNode* leaf)
		{
			btSoftBody::Cluster* cluster = (btSoftBody::Cluster*)leaf->data;
			btSoftClusterCollisionShape cshape(cluster);

			const btConvexShape* rshape = (const btConvexShape*)m_colObjWrap->getCollisionShape();

			///don't collide an anchored cluster with a static/kinematic object
			if (m_colObjWrap->getCollisionObject()->isStaticOrKinematicObject() && cluster->m_containsAnchor)
				return;

			btGjkEpaSolver2::sResults res;
			if (btGjkEpaSolver2::SignedDistance(&cshape, btTransform::getIdentity(),
												rshape, m_colObjWrap->getWorldTransform(),
												btVector3(1, 0, 0), res))
			{
				btSoftBody::CJoint joint;
				if (SolveContact(res, cluster, m_colObjWrap->getCollisionObject(), joint))  //prb,joint))
				{
					btSoftBody::CJoint* pj = new (btAlignedAlloc(sizeof(btSoftBody::CJoint), 16)) btSoftBody::CJoint();
					*pj = joint;
					psb->m_joints.push_back(pj);
					if (m_colObjWrap->getCollisionObject()->isStaticOrKinematicObject())
					{
						pj->m_erp *= psb->m_cfg.kSKHR_CL;
						pj->m_split *= psb->m_cfg.kSK_SPLT_CL;
					}
					else
					{
						pj->m_erp *= psb->m_cfg.kSRHR_CL;
						pj->m_split *= psb->m_cfg.kSR_SPLT_CL;
					}
				}
			}
		}
		void ProcessColObj(btSoftBody* ps, const btCollisionObjectWrapper* colObWrap)
		{
			psb = ps;
			m_colObjWrap = colObWrap;
			idt = ps->m_sst.isdt;
			m_margin = m_colObjWrap->getCollisionShape()->getMargin() + psb->getCollisionShape()->getMargin();
			///Bullet rigid body uses multiply instead of minimum to determine combined friction. Some customization would be useful.
			friction = btMin(psb->m_cfg.kDF, m_colObjWrap->getCollisionObject()->getFriction());
			btVector3 mins;
			btVector3 maxs;

			ATTRIBUTE_ALIGNED16(btDbvtVolume)
			volume;
			colObWrap->getCollisionShape()->getAabb(colObWrap->getWorldTransform(), mins, maxs);
			volume = btDbvtVolume::FromMM(mins, maxs);
			volume.Expand(btVector3(1, 1, 1) * m_margin);
			ps->m_cdbvt.collideTV(ps->m_cdbvt.m_root, volume, *this);
		}
	};
	//
	// CollideCL_SS
	//
	struct CollideCL_SS : ClusterBase
	{
		btSoftBody* bodies[2];
		void Process(const btDbvtNode* la, const btDbvtNode* lb)
		{
			btSoftBody::Cluster* cla = (btSoftBody::Cluster*)la->data;
			btSoftBody::Cluster* clb = (btSoftBody::Cluster*)lb->data;

			bool connected = false;
			if ((bodies[0] == bodies[1]) && (bodies[0]->m_clusterConnectivity.size()))
			{
				connected = bodies[0]->m_clusterConnectivity[cla->m_clusterIndex + bodies[0]->m_clusters.size() * clb->m_clusterIndex];
			}

			if (!connected)
			{
				btSoftClusterCollisionShape csa(cla);
				btSoftClusterCollisionShape csb(clb);
				btGjkEpaSolver2::sResults res;
				if (btGjkEpaSolver2::SignedDistance(&csa, btTransform::getIdentity(),
													&csb, btTransform::getIdentity(),
													cla->m_com - clb->m_com, res))
				{
					btSoftBody::CJoint joint;
					if (SolveContact(res, cla, clb, joint))
					{
						btSoftBody::CJoint* pj = new (btAlignedAlloc(sizeof(btSoftBody::CJoint), 16)) btSoftBody::CJoint();
						*pj = joint;
						bodies[0]->m_joints.push_back(pj);
						pj->m_erp *= btMax(bodies[0]->m_cfg.kSSHR_CL, bodies[1]->m_cfg.kSSHR_CL);
						pj->m_split *= (bodies[0]->m_cfg.kSS_SPLT_CL + bodies[1]->m_cfg.kSS_SPLT_CL) / 2;
					}
				}
			}
			else
			{
				static int count = 0;
				count++;
				//printf("count=%d\n",count);
			}
		}
		void ProcessSoftSoft(btSoftBody* psa, btSoftBody* psb)
		{
			idt = psa->m_sst.isdt;
			//m_margin		=	(psa->getCollisionShape()->getMargin()+psb->getCollisionShape()->getMargin())/2;
			m_margin = (psa->getCollisionShape()->getMargin() + psb->getCollisionShape()->getMargin());
			friction = btMin(psa->m_cfg.kDF, psb->m_cfg.kDF);
			bodies[0] = psa;
			bodies[1] = psb;
			psa->m_cdbvt.collideTT(psa->m_cdbvt.m_root, psb->m_cdbvt.m_root, *this);
		}
	};
	//
	// CollideSDF_RS
	//
	struct CollideSDF_RS : btDbvt::ICollide
	{
		void Process(const btDbvtNode* leaf)
		{
			btSoftBody::Node* node = (btSoftBody::Node*)leaf->data;
			DoNode(*node);
		}
		void DoNode(btSoftBody::Node& n) const
		{
			const btScalar m = n.m_im > 0 ? dynmargin : stamargin;
			btSoftBody::RContact c;

			if ((!n.m_battach) &&
				psb->checkContact(m_colObj1Wrap, n.m_x, m, c.m_cti))
			{
				const btScalar ima = n.m_im;
				const btScalar imb = m_rigidBody ? m_rigidBody->getInvMass() : 0.f;
				const btScalar ms = ima + imb;
				if (ms > 0)
				{
					const btTransform& wtr = m_rigidBody ? m_rigidBody->getWorldTransform() : m_colObj1Wrap->getCollisionObject()->getWorldTransform();
					static const btMatrix3x3 iwiStatic(0, 0, 0, 0, 0, 0, 0, 0, 0);
					const btMatrix3x3& iwi = m_rigidBody ? m_rigidBody->getInvInertiaTensorWorld() : iwiStatic;
					const btVector3 ra = n.m_x - wtr.getOrigin();
					const btVector3 va = m_rigidBody ? m_rigidBody->getVelocityInLocalPoint(ra) * psb->m_sst.sdt : btVector3(0, 0, 0);
					const btVector3 vb = n.m_x - n.m_q;
					const btVector3 vr = vb - va;
					const btScalar dn = btDot(vr, c.m_cti.m_normal);
					const btVector3 fv = vr - c.m_cti.m_normal * dn;
					const btScalar fc = psb->m_cfg.kDF * m_colObj1Wrap->getCollisionObject()->getFriction();
					c.m_node = &n;
					c.m_c0 = ImpulseMatrix(psb->m_sst.sdt, ima, imb, iwi, ra);
					c.m_c1 = ra;
					c.m_c2 = ima * psb->m_sst.sdt;
					c.m_c3 = fv.length2() < (dn * fc * dn * fc) ? 0 : 1 - fc;
					c.m_c4 = m_colObj1Wrap->getCollisionObject()->isStaticOrKinematicObject() ? psb->m_cfg.kKHR : psb->m_cfg.kCHR;
					psb->m_rcontacts.push_back(c);
					if (m_rigidBody)
						m_rigidBody->activate();
				}
			}
		}
		btSoftBody* psb;
		const btCollisionObjectWrapper* m_colObj1Wrap;
		btRigidBody* m_rigidBody;
		btScalar dynmargin;
		btScalar stamargin;
	};

	//
	// CollideSDF_RD
	//
	struct CollideSDF_RD : btDbvt::ICollide
	{
		void Process(const btDbvtNode* leaf)
		{
			btSoftBody::Node* node = (btSoftBody::Node*)leaf->data;
			DoNode(*node);
		}
		void DoNode(btSoftBody::Node& n) const
		{
			const btScalar m = n.m_im > 0 ? dynmargin : stamargin;
			btSoftBody::DeformableNodeRigidContact c;

			if (!n.m_battach)
			{
				// check for collision at x_{n+1}^*
				if (psb->checkDeformableContact(m_colObj1Wrap, n.m_q, m, c.m_cti, /*predict = */ true))
				{
					const btScalar ima = n.m_im;
					// todo: collision between multibody and fixed deformable node will be missed.
					const btScalar imb = m_rigidBody ? m_rigidBody->getInvMass() : 0.f;
					const btScalar ms = ima + imb;
					if (ms > 0)
					{
						// resolve contact at x_n
						psb->checkDeformableContact(m_colObj1Wrap, n.m_x, m, c.m_cti, /*predict = */ false);
						btSoftBody::sCti& cti = c.m_cti;
						c.m_node = &n;
						const btScalar fc = psb->m_cfg.kDF * m_colObj1Wrap->getCollisionObject()->getFriction();
						c.m_c2 = ima;
						c.m_c3 = fc;
						c.m_c4 = m_colObj1Wrap->getCollisionObject()->isStaticOrKinematicObject() ? psb->m_cfg.kKHR : psb->m_cfg.kCHR;
						c.m_c5 = n.m_effectiveMass_inv;

						if (cti.m_colObj->getInternalType() == btCollisionObject::CO_RIGID_BODY)
						{
							const btTransform& wtr = m_rigidBody ? m_rigidBody->getWorldTransform() : m_colObj1Wrap->getCollisionObject()->getWorldTransform();
							const btVector3 ra = n.m_x - wtr.getOrigin();

							static const btMatrix3x3 iwiStatic(0, 0, 0, 0, 0, 0, 0, 0, 0);
							const btMatrix3x3& iwi = m_rigidBody ? m_rigidBody->getInvInertiaTensorWorld() : iwiStatic;
							if (psb->m_reducedModel)
							{
								c.m_c0 = MassMatrix(imb, iwi, ra);  //impulse factor K of the rigid body only (not the inverse)
							}
							else
							{
								c.m_c0 = ImpulseMatrix(1, n.m_effectiveMass_inv, imb, iwi, ra);
								//                            c.m_c0 = ImpulseMatrix(1, ima, imb, iwi, ra);
							}
							c.m_c1 = ra;
						}
						else if (cti.m_colObj->getInternalType() == btCollisionObject::CO_FEATHERSTONE_LINK)
						{
							btMultiBodyLinkCollider* multibodyLinkCol = (btMultiBodyLinkCollider*)btMultiBodyLinkCollider::upcast(cti.m_colObj);
							if (multibodyLinkCol)
							{
								btVector3 normal = cti.m_normal;
								btVector3 t1 = generateUnitOrthogonalVector(normal);
								btVector3 t2 = btCross(normal, t1);
								btMultiBodyJacobianData jacobianData_normal, jacobianData_t1, jacobianData_t2;
								findJacobian(multibodyLinkCol, jacobianData_normal, c.m_node->m_x, normal);
								findJacobian(multibodyLinkCol, jacobianData_t1, c.m_node->m_x, t1);
								findJacobian(multibodyLinkCol, jacobianData_t2, c.m_node->m_x, t2);

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

								btMatrix3x3 local_impulse_matrix;
								if (psb->m_reducedModel)
								{
									local_impulse_matrix = OuterProduct(J_n, J_t1, J_t2, u_n, u_t1, u_t2, ndof);
								}
								else
								{
									local_impulse_matrix = (n.m_effectiveMass_inv + OuterProduct(J_n, J_t1, J_t2, u_n, u_t1, u_t2, ndof)).inverse();
								}
								c.m_c0 = rot.transpose() * local_impulse_matrix * rot;
								c.jacobianData_normal = jacobianData_normal;
								c.jacobianData_t1 = jacobianData_t1;
								c.jacobianData_t2 = jacobianData_t2;
								c.t1 = t1;
								c.t2 = t2;
							}
						}
						psb->m_nodeRigidContacts.push_back(c);
					}
				}
			}
		}
		btSoftBody* psb;
		const btCollisionObjectWrapper* m_colObj1Wrap;
		btRigidBody* m_rigidBody;
		btScalar dynmargin;
		btScalar stamargin;
	};

	//
	// CollideSDF_RDF
	//
	struct CollideSDF_RDF : btDbvt::ICollide
	{
		void Process(const btDbvtNode* leaf)
		{
			btSoftBody::Face* face = (btSoftBody::Face*)leaf->data;
			DoNode(*face);
		}
		void DoNode(btSoftBody::Face& f) const
		{
			btSoftBody::Node* n0 = f.m_n[0];
			btSoftBody::Node* n1 = f.m_n[1];
			btSoftBody::Node* n2 = f.m_n[2];
			const btScalar m = (n0->m_im > 0 && n1->m_im > 0 && n2->m_im > 0) ? dynmargin : stamargin;
			btSoftBody::DeformableFaceRigidContact c;
			btVector3 contact_point;
			btVector3 bary;
			if (psb->checkDeformableFaceContact(m_colObj1Wrap, f, contact_point, bary, m, c.m_cti, true))
			{
				btScalar ima = n0->m_im + n1->m_im + n2->m_im;
				const btScalar imb = m_rigidBody ? m_rigidBody->getInvMass() : 0.f;
				// todo: collision between multibody and fixed deformable face will be missed.
				const btScalar ms = ima + imb;
				if (ms > 0)
				{
					// resolve contact at x_n
					//                    psb->checkDeformableFaceContact(m_colObj1Wrap, f, contact_point, bary, m, c.m_cti, /*predict = */ false);
					btSoftBody::sCti& cti = c.m_cti;
					c.m_contactPoint = contact_point;

					c.m_bary = bary;
					// todo xuchenhan@: this is assuming mass of all vertices are the same. Need to modify if mass are different for distinct vertices
					c.m_weights = btScalar(2) / (btScalar(1) + bary.length2()) * bary;
					c.m_face = &f;
					// friction is handled by the nodes to prevent sticking
					//                    const btScalar fc = 0;
					const btScalar fc = psb->m_cfg.kDF * m_colObj1Wrap->getCollisionObject()->getFriction();

					// the effective inverse mass of the face as in https://graphics.stanford.edu/papers/cloth-sig02/cloth.pdf
					ima = bary.getX() * c.m_weights.getX() * n0->m_im + bary.getY() * c.m_weights.getY() * n1->m_im + bary.getZ() * c.m_weights.getZ() * n2->m_im;
					c.m_c2 = ima;
					c.m_c3 = fc;
					c.m_c4 = m_colObj1Wrap->getCollisionObject()->isStaticOrKinematicObject() ? psb->m_cfg.kKHR : psb->m_cfg.kCHR;
					c.m_c5 = Diagonal(ima);
					if (cti.m_colObj->getInternalType() == btCollisionObject::CO_RIGID_BODY)
					{
						const btTransform& wtr = m_rigidBody ? m_rigidBody->getWorldTransform() : m_colObj1Wrap->getCollisionObject()->getWorldTransform();
						static const btMatrix3x3 iwiStatic(0, 0, 0, 0, 0, 0, 0, 0, 0);
						const btMatrix3x3& iwi = m_rigidBody ? m_rigidBody->getInvInertiaTensorWorld() : iwiStatic;
						const btVector3 ra = contact_point - wtr.getOrigin();

						// we do not scale the impulse matrix by dt
						c.m_c0 = ImpulseMatrix(1, ima, imb, iwi, ra);
						c.m_c1 = ra;
					}
					else if (cti.m_colObj->getInternalType() == btCollisionObject::CO_FEATHERSTONE_LINK)
					{
						btMultiBodyLinkCollider* multibodyLinkCol = (btMultiBodyLinkCollider*)btMultiBodyLinkCollider::upcast(cti.m_colObj);
						if (multibodyLinkCol)
						{
							btVector3 normal = cti.m_normal;
							btVector3 t1 = generateUnitOrthogonalVector(normal);
							btVector3 t2 = btCross(normal, t1);
							btMultiBodyJacobianData jacobianData_normal, jacobianData_t1, jacobianData_t2;
							findJacobian(multibodyLinkCol, jacobianData_normal, contact_point, normal);
							findJacobian(multibodyLinkCol, jacobianData_t1, contact_point, t1);
							findJacobian(multibodyLinkCol, jacobianData_t2, contact_point, t2);

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
							btMatrix3x3 local_impulse_matrix = (Diagonal(ima) + OuterProduct(J_n, J_t1, J_t2, u_n, u_t1, u_t2, ndof)).inverse();
							c.m_c0 = rot.transpose() * local_impulse_matrix * rot;
							c.jacobianData_normal = jacobianData_normal;
							c.jacobianData_t1 = jacobianData_t1;
							c.jacobianData_t2 = jacobianData_t2;
							c.t1 = t1;
							c.t2 = t2;
						}
					}
					psb->m_faceRigidContacts.push_back(c);
				}
			}
			// Set caching barycenters to be false after collision detection.
			// Only turn on when contact is static.
			f.m_pcontact[3] = 0;
		}
		btSoftBody* psb;
		const btCollisionObjectWrapper* m_colObj1Wrap;
		btRigidBody* m_rigidBody;
		btScalar dynmargin;
		btScalar stamargin;
	};

	//
	// CollideVF_SS
	//
	struct CollideVF_SS : btDbvt::ICollide
	{
		void Process(const btDbvtNode* lnode,
					 const btDbvtNode* lface)
		{
			btSoftBody::Node* node = (btSoftBody::Node*)lnode->data;
			btSoftBody::Face* face = (btSoftBody::Face*)lface->data;
			for (int i = 0; i < 3; ++i)
			{
				if (face->m_n[i] == node)
					continue;
			}

			btVector3 o = node->m_x;
			btVector3 p;
			btScalar d = SIMD_INFINITY;
			ProjectOrigin(face->m_n[0]->m_x - o,
						  face->m_n[1]->m_x - o,
						  face->m_n[2]->m_x - o,
						  p, d);
			const btScalar m = mrg + (o - node->m_q).length() * 2;
			if (d < (m * m))
			{
				const btSoftBody::Node* n[] = {face->m_n[0], face->m_n[1], face->m_n[2]};
				const btVector3 w = BaryCoord(n[0]->m_x, n[1]->m_x, n[2]->m_x, p + o);
				const btScalar ma = node->m_im;
				btScalar mb = BaryEval(n[0]->m_im, n[1]->m_im, n[2]->m_im, w);
				if ((n[0]->m_im <= 0) ||
					(n[1]->m_im <= 0) ||
					(n[2]->m_im <= 0))
				{
					mb = 0;
				}
				const btScalar ms = ma + mb;
				if (ms > 0)
				{
					btSoftBody::SContact c;
					c.m_normal = p / -btSqrt(d);
					c.m_margin = m;
					c.m_node = node;
					c.m_face = face;
					c.m_weights = w;
					c.m_friction = btMax(psb[0]->m_cfg.kDF, psb[1]->m_cfg.kDF);
					c.m_cfm[0] = ma / ms * psb[0]->m_cfg.kSHR;
					c.m_cfm[1] = mb / ms * psb[1]->m_cfg.kSHR;
					psb[0]->m_scontacts.push_back(c);
				}
			}
		}
		btSoftBody* psb[2];
		btScalar mrg;
	};

	//
	// CollideVF_DD
	//
	struct CollideVF_DD : btDbvt::ICollide
	{
		void Process(const btDbvtNode* lnode,
					 const btDbvtNode* lface)
		{
			btSoftBody::Node* node = (btSoftBody::Node*)lnode->data;
			btSoftBody::Face* face = (btSoftBody::Face*)lface->data;
			btVector3 bary;
			if (proximityTest(face->m_n[0]->m_x, face->m_n[1]->m_x, face->m_n[2]->m_x, node->m_x, face->m_normal, mrg, bary))
			{
				const btSoftBody::Node* n[] = {face->m_n[0], face->m_n[1], face->m_n[2]};
				const btVector3 w = bary;
				const btScalar ma = node->m_im;
				btScalar mb = BaryEval(n[0]->m_im, n[1]->m_im, n[2]->m_im, w);
				if ((n[0]->m_im <= 0) ||
					(n[1]->m_im <= 0) ||
					(n[2]->m_im <= 0))
				{
					mb = 0;
				}
				const btScalar ms = ma + mb;
				if (ms > 0)
				{
					btSoftBody::DeformableFaceNodeContact c;
					c.m_normal = face->m_normal;
					if (!useFaceNormal && c.m_normal.dot(node->m_x - face->m_n[2]->m_x) < 0)
						c.m_normal = -face->m_normal;
					c.m_margin = mrg;
					c.m_node = node;
					c.m_face = face;
					c.m_bary = w;
					c.m_friction = psb[0]->m_cfg.kDF * psb[1]->m_cfg.kDF;
					// Initialize unused fields.
					c.m_weights = btVector3(0, 0, 0);
					c.m_imf = 0;
					c.m_c0 = 0;
					c.m_colObj = psb[1];
					psb[0]->m_faceNodeContacts.push_back(c);
				}
			}
		}
		btSoftBody* psb[2];
		btScalar mrg;
		bool useFaceNormal;
	};

	//
	// CollideFF_DD
	//
	struct CollideFF_DD : btDbvt::ICollide
	{
		void Process(const btDbvntNode* lface1,
					 const btDbvntNode* lface2)
		{
			btSoftBody::Face* f1 = (btSoftBody::Face*)lface1->data;
			btSoftBody::Face* f2 = (btSoftBody::Face*)lface2->data;
			if (f1 != f2)
			{
				Repel(f1, f2);
				Repel(f2, f1);
			}
		}
		void Repel(btSoftBody::Face* f1, btSoftBody::Face* f2)
		{
			//#define REPEL_NEIGHBOR 1
#ifndef REPEL_NEIGHBOR
			for (int node_id = 0; node_id < 3; ++node_id)
			{
				btSoftBody::Node* node = f1->m_n[node_id];
				for (int i = 0; i < 3; ++i)
				{
					if (f2->m_n[i] == node)
						return;
				}
			}
#endif
			bool skip = false;
			for (int node_id = 0; node_id < 3; ++node_id)
			{
				btSoftBody::Node* node = f1->m_n[node_id];
#ifdef REPEL_NEIGHBOR
				for (int i = 0; i < 3; ++i)
				{
					if (f2->m_n[i] == node)
					{
						skip = true;
						break;
					}
				}
				if (skip)
				{
					skip = false;
					continue;
				}
#endif
				btSoftBody::Face* face = f2;
				btVector3 bary;
				if (!proximityTest(face->m_n[0]->m_x, face->m_n[1]->m_x, face->m_n[2]->m_x, node->m_x, face->m_normal, mrg, bary))
					continue;
				btSoftBody::DeformableFaceNodeContact c;
				c.m_normal = face->m_normal;
				if (!useFaceNormal && c.m_normal.dot(node->m_x - face->m_n[2]->m_x) < 0)
					c.m_normal = -face->m_normal;
				c.m_margin = mrg;
				c.m_node = node;
				c.m_face = face;
				c.m_bary = bary;
				c.m_friction = psb[0]->m_cfg.kDF * psb[1]->m_cfg.kDF;
				// Initialize unused fields.
				c.m_weights = btVector3(0, 0, 0);
				c.m_imf = 0;
				c.m_c0 = 0;
				c.m_colObj = psb[1];
				psb[0]->m_faceNodeContacts.push_back(c);
			}
		}

		btSoftBody* psb[2];
		btScalar mrg;
		bool useFaceNormal;
	};

	struct CollideCCD : btDbvt::ICollide
	{
		void Process(const btDbvtNode* lnode,
					 const btDbvtNode* lface)
		{
			btSoftBody::Node* node = (btSoftBody::Node*)lnode->data;
			btSoftBody::Face* face = (btSoftBody::Face*)lface->data;
			btVector3 bary;
			if (bernsteinCCD(face, node, dt, SAFE_EPSILON, bary))
			{
				btSoftBody::DeformableFaceNodeContact c;
				c.m_normal = face->m_normal;
				if (!useFaceNormal && c.m_normal.dot(node->m_x - face->m_n[2]->m_x) < 0)
					c.m_normal = -face->m_normal;
				c.m_node = node;
				c.m_face = face;
				c.m_bary = bary;
				c.m_friction = psb[0]->m_cfg.kDF * psb[1]->m_cfg.kDF;
				// Initialize unused fields.
				c.m_weights = btVector3(0, 0, 0);
				c.m_margin = mrg;
				c.m_imf = 0;
				c.m_c0 = 0;
				c.m_colObj = psb[1];
				psb[0]->m_faceNodeContactsCCD.push_back(c);
			}
		}
		void Process(const btDbvntNode* lface1,
					 const btDbvntNode* lface2)
		{
			btSoftBody::Face* f1 = (btSoftBody::Face*)lface1->data;
			btSoftBody::Face* f2 = (btSoftBody::Face*)lface2->data;
			if (f1 != f2)
			{
				Repel(f1, f2);
				Repel(f2, f1);
			}
		}
		void Repel(btSoftBody::Face* f1, btSoftBody::Face* f2)
		{
			//#define REPEL_NEIGHBOR 1
#ifndef REPEL_NEIGHBOR
			for (int node_id = 0; node_id < 3; ++node_id)
			{
				btSoftBody::Node* node = f1->m_n[node_id];
				for (int i = 0; i < 3; ++i)
				{
					if (f2->m_n[i] == node)
						return;
				}
			}
#endif
			bool skip = false;
			for (int node_id = 0; node_id < 3; ++node_id)
			{
				btSoftBody::Node* node = f1->m_n[node_id];
#ifdef REPEL_NEIGHBOR
				for (int i = 0; i < 3; ++i)
				{
					if (f2->m_n[i] == node)
					{
						skip = true;
						break;
					}
				}
				if (skip)
				{
					skip = false;
					continue;
				}
#endif
				btSoftBody::Face* face = f2;
				btVector3 bary;
				if (bernsteinCCD(face, node, dt, SAFE_EPSILON, bary))
				{
					btSoftBody::DeformableFaceNodeContact c;
					c.m_normal = face->m_normal;
					if (!useFaceNormal && c.m_normal.dot(node->m_x - face->m_n[2]->m_x) < 0)
						c.m_normal = -face->m_normal;
					c.m_node = node;
					c.m_face = face;
					c.m_bary = bary;
					c.m_friction = psb[0]->m_cfg.kDF * psb[1]->m_cfg.kDF;
					// Initialize unused fields.
					c.m_weights = btVector3(0, 0, 0);
					c.m_margin = mrg;
					c.m_imf = 0;
					c.m_c0 = 0;
					c.m_colObj = psb[1];
					psb[0]->m_faceNodeContactsCCD.push_back(c);
				}
			}
		}
		btSoftBody* psb[2];
		btScalar dt, mrg;
		bool useFaceNormal;
	};
};
#endif  //_BT_SOFT_BODY_INTERNALS_H
/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/


#include "btSoftBodyRigidBodyCollisionConfiguration.h"
#include "btSoftRigidCollisionAlgorithm.h"
#include "btSoftBodyConcaveCollisionAlgorithm.h"
#include "btSoftSoftCollisionAlgorithm.h"

#include "LinearMath/btPoolAllocator.h"

#define ENABLE_SOFTBODY_CONCAVE_COLLISIONS 1

btSoftBodyRigidBodyCollisionConfiguration::btSoftBodyRigidBodyCollisionConfiguration(const btDefaultCollisionConstructionInfo& constructionInfo)
	: btDefaultCollisionConfiguration(constructionInfo)
{
	void* mem;

	mem = btAlignedAlloc(sizeof(btSoftSoftCollisionAlgorithm::CreateFunc), 16);
	m_softSoftCreateFunc = new (mem) btSoftSoftCollisionAlgorithm::CreateFunc;

	mem = btAlignedAlloc(sizeof(btSoftRigidCollisionAlgorithm::CreateFunc), 16);
	m_softRigidConvexCreateFunc = new (mem) btSoftRigidCollisionAlgorithm::CreateFunc;

	mem = btAlignedAlloc(sizeof(btSoftRigidCollisionAlgorithm::CreateFunc), 16);
	m_swappedSoftRigidConvexCreateFunc = new (mem) btSoftRigidCollisionAlgorithm::CreateFunc;
	m_swappedSoftRigidConvexCreateFunc->m_swapped = true;

#ifdef ENABLE_SOFTBODY_CONCAVE_COLLISIONS
	mem = btAlignedAlloc(sizeof(btSoftBodyConcaveCollisionAlgorithm::CreateFunc), 16);
	m_softRigidConcaveCreateFunc = new (mem) btSoftBodyConcaveCollisionAlgorithm::CreateFunc;

	mem = btAlignedAlloc(sizeof(btSoftBodyConcaveCollisionAlgorithm::CreateFunc), 16);
	m_swappedSoftRigidConcaveCreateFunc = new (mem) btSoftBodyConcaveCollisionAlgorithm::SwappedCreateFunc;
	m_swappedSoftRigidConcaveCreateFunc->m_swapped = true;
#endif

	//replace pool by a new one, with potential larger size

	if (m_ownsCollisionAlgorithmPool && m_collisionAlgorithmPool)
	{
		int curElemSize = m_collisionAlgorithmPool->getElementSize();
		///calculate maximum element size, big enough to fit any collision algorithm in the memory pool

		int maxSize0 = sizeof(btSoftSoftCollisionAlgorithm);
		int maxSize1 = sizeof(btSoftRigidCollisionAlgorithm);
		int maxSize2 = sizeof(btSoftBodyConcaveCollisionAlgorithm);

		int collisionAlgorithmMaxElementSize = btMax(maxSize0, maxSize1);
		collisionAlgorithmMaxElementSize = btMax(collisionAlgorithmMaxElementSize, maxSize2);

		if (collisionAlgorithmMaxElementSize > curElemSize)
		{
			m_collisionAlgorithmPool->~btPoolAllocator();
			btAlignedFree(m_collisionAlgorithmPool);
			void* mem = btAlignedAlloc(sizeof(btPoolAllocator), 16);
			m_collisionAlgorithmPool = new (mem) btPoolAllocator(collisionAlgorithmMaxElementSize, constructionInfo.m_defaultMaxCollisionAlgorithmPoolSize);
		}
	}
}

btSoftBodyRigidBodyCollisionConfiguration::~btSoftBodyRigidBodyCollisionConfiguration()
{
	m_softSoftCreateFunc->~btCollisionAlgorithmCreateFunc();
	btAlignedFree(m_softSoftCreateFunc);

	m_softRigidConvexCreateFunc->~btCollisionAlgorithmCreateFunc();
	btAlignedFree(m_softRigidConvexCreateFunc);

	m_swappedSoftRigidConvexCreateFunc->~btCollisionAlgorithmCreateFunc();
	btAlignedFree(m_swappedSoftRigidConvexCreateFunc);

#ifdef ENABLE_SOFTBODY_CONCAVE_COLLISIONS
	m_softRigidConcaveCreateFunc->~btCollisionAlgorithmCreateFunc();
	btAlignedFree(m_softRigidConcaveCreateFunc);

	m_swappedSoftRigidConcaveCreateFunc->~btCollisionAlgorithmCreateFunc();
	btAlignedFree(m_swappedSoftRigidConcaveCreateFunc);
#endif
}

///creation of soft-soft and soft-rigid, and otherwise fallback to base class implementation
btCollisionAlgorithmCreateFunc* btSoftBodyRigidBodyCollisionConfiguration::getCollisionAlgorithmCreateFunc(int proxyType0, int proxyType1)
{
	///try to handle the softbody interactions first

	if ((proxyType0 == SOFTBODY_SHAPE_PROXYTYPE) && (proxyType1 == SOFTBODY_SHAPE_PROXYTYPE))
	{
		return m_softSoftCreateFunc;
	}

	///softbody versus convex
	if (proxyType0 == SOFTBODY_SHAPE_PROXYTYPE && btBroadphaseProxy::isConvex(proxyType1))
	{
		return m_softRigidConvexCreateFunc;
	}

	///convex versus soft body
	if (btBroadphaseProxy::isConvex(proxyType0) && proxyType1 == SOFTBODY_SHAPE_PROXYTYPE)
	{
		return m_swappedSoftRigidConvexCreateFunc;
	}

#ifdef ENABLE_SOFTBODY_CONCAVE_COLLISIONS
	///softbody versus convex
	if (proxyType0 == SOFTBODY_SHAPE_PROXYTYPE && btBroadphaseProxy::isConcave(proxyType1))
	{
		return m_softRigidConcaveCreateFunc;
	}

	///convex versus soft body
	if (btBroadphaseProxy::isConcave(proxyType0) && proxyType1 == SOFTBODY_SHAPE_PROXYTYPE)
	{
		return m_swappedSoftRigidConcaveCreateFunc;
	}
#endif

	///fallback to the regular rigid collision shape
	return btDefaultCollisionConfiguration::getCollisionAlgorithmCreateFunc(proxyType0, proxyType1);
}
/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/


#ifndef BT_SOFTBODY_RIGIDBODY_COLLISION_CONFIGURATION
#define BT_SOFTBODY_RIGIDBODY_COLLISION_CONFIGURATION

#include "BulletCollision/CollisionDispatch/btDefaultCollisionConfiguration.h"

class btVoronoiSimplexSolver;
class btGjkEpaPenetrationDepthSolver;

///btSoftBodyRigidBodyCollisionConfiguration add softbody interaction on top of btDefaultCollisionConfiguration
class btSoftBodyRigidBodyCollisionConfiguration : public btDefaultCollisionConfiguration
{
	//default CreationFunctions, filling the m_doubleDispatch table
	btCollisionAlgorithmCreateFunc* m_softSoftCreateFunc;
	btCollisionAlgorithmCreateFunc* m_softRigidConvexCreateFunc;
	btCollisionAlgorithmCreateFunc* m_swappedSoftRigidConvexCreateFunc;
	btCollisionAlgorithmCreateFunc* m_softRigidConcaveCreateFunc;
	btCollisionAlgorithmCreateFunc* m_swappedSoftRigidConcaveCreateFunc;

public:
	btSoftBodyRigidBodyCollisionConfiguration(const btDefaultCollisionConstructionInfo& constructionInfo = btDefaultCollisionConstructionInfo());

	virtual ~btSoftBodyRigidBodyCollisionConfiguration();

	///creation of soft-soft and soft-rigid, and otherwise fallback to base class implementation
	virtual btCollisionAlgorithmCreateFunc* getCollisionAlgorithmCreateFunc(int proxyType0, int proxyType1);
};

#endif  //BT_SOFTBODY_RIGIDBODY_COLLISION_CONFIGURATION
/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/

#ifndef BT_SOFT_BODY_SOLVERS_H
#define BT_SOFT_BODY_SOLVERS_H

#include "BulletCollision/CollisionShapes/btTriangleIndexVertexArray.h"

class btSoftBodyTriangleData;
class btSoftBodyLinkData;
class btSoftBodyVertexData;
class btVertexBufferDescriptor;
class btCollisionObject;
class btSoftBody;
struct btDispatcherInfo;
class btManifoldResult;
class btManifoldResultForSkin;

class btSoftBodySolver
{
public:
	enum SolverTypes
	{
		DEFAULT_SOLVER,
		CPU_SOLVER,
		CL_SOLVER,
		CL_SIMD_SOLVER,
		DX_SOLVER,
		DX_SIMD_SOLVER,
		DEFORMABLE_SOLVER,
		REDUCED_DEFORMABLE_SOLVER
	};

protected:
	int m_numberOfPositionIterations;
	int m_numberOfVelocityIterations;
	// Simulation timescale
	float m_timeScale;

public:
	btSoftBodySolver() : m_numberOfPositionIterations(10),
						 m_timeScale(1)
	{
		m_numberOfVelocityIterations = 0;
		m_numberOfPositionIterations = 5;
	}

	virtual ~btSoftBodySolver()
	{
	}

	/**
	 * Return the type of the solver.
	 */
	virtual SolverTypes getSolverType() const = 0;

	/** Ensure that this solver is initialized. */
	virtual bool checkInitialized() = 0;

	/** Optimize soft bodies in this solver. */
	virtual void optimize(btAlignedObjectArray<btSoftBody *> &softBodies, bool forceUpdate = false) = 0;

	/** Copy necessary data back to the original soft body source objects. */
	virtual void copyBackToSoftBodies(bool bMove = true) = 0;

	/** Predict motion of soft bodies into next timestep */
	virtual void predictMotion(btScalar solverdt) = 0;

	/** Solve constraints for a set of soft bodies */
	virtual void solveConstraints(btScalar solverdt) = 0;

	/** Perform necessary per-step updates of soft bodies such as recomputing normals and bounding boxes */
	virtual void updateSoftBodies() = 0;

	/** Process a collision between one of the world's soft bodies and another collision object */
	virtual void processCollision(btSoftBody *, const struct btCollisionObjectWrapper *, btManifoldResultForSkin *) = 0;

	/** Process a collision between two soft bodies */
	virtual void processCollision(btSoftBody *, btSoftBody *) = 0;

	virtual void processCollision(btSoftBody *, btSoftBody *, btManifoldResultForSkin *) = 0;

	/** Set the number of velocity constraint solver iterations this solver uses. */
	virtual void setNumberOfPositionIterations(int iterations)
	{
		m_numberOfPositionIterations = iterations;
	}

	/** Get the number of velocity constraint solver iterations this solver uses. */
	virtual int getNumberOfPositionIterations()
	{
		return m_numberOfPositionIterations;
	}

	/** Set the number of velocity constraint solver iterations this solver uses. */
	virtual void setNumberOfVelocityIterations(int iterations)
	{
		m_numberOfVelocityIterations = iterations;
	}

	/** Get the number of velocity constraint solver iterations this solver uses. */
	virtual int getNumberOfVelocityIterations()
	{
		return m_numberOfVelocityIterations;
	}

	/** Return the timescale that the simulation is using */
	float getTimeScale()
	{
		return m_timeScale;
	}

#if 0
	/**
	 * Add a collision object to be used by the indicated softbody.
	 */
	virtual void addCollisionObjectForSoftBody( int clothIdentifier, btCollisionObject *collisionObject ) = 0;
#endif
};

/** 
 * Class to manage movement of data from a solver to a given target.
 * This version is abstract. Subclasses will have custom pairings for different combinations.
 */
class btSoftBodySolverOutput
{
protected:
public:
	btSoftBodySolverOutput()
	{
	}

	virtual ~btSoftBodySolverOutput()
	{
	}

	/** Output current computed vertex data to the vertex buffers for all cloths in the solver. */
	virtual void copySoftBodyToVertexBuffer(const btSoftBody *const softBody, btVertexBufferDescriptor *vertexBuffer) = 0;
};

#endif  // #ifndef BT_SOFT_BODY_SOLVERS_H
/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/


#ifndef BT_SOFT_BODY_SOLVER_VERTEX_BUFFER_H
#define BT_SOFT_BODY_SOLVER_VERTEX_BUFFER_H

class btVertexBufferDescriptor
{
public:
	enum BufferTypes
	{
		CPU_BUFFER,
		DX11_BUFFER,
		OPENGL_BUFFER
	};

protected:
	bool m_hasVertexPositions;
	bool m_hasNormals;

	int m_vertexOffset;
	int m_vertexStride;

	int m_normalOffset;
	int m_normalStride;

public:
	btVertexBufferDescriptor()
	{
		m_hasVertexPositions = false;
		m_hasNormals = false;
		m_vertexOffset = 0;
		m_vertexStride = 0;
		m_normalOffset = 0;
		m_normalStride = 0;
	}

	virtual ~btVertexBufferDescriptor()
	{
	}

	virtual bool hasVertexPositions() const
	{
		return m_hasVertexPositions;
	}

	virtual bool hasNormals() const
	{
		return m_hasNormals;
	}

	/**
	 * Return the type of the vertex buffer descriptor.
	 */
	virtual BufferTypes getBufferType() const = 0;

	/**
	 * Return the vertex offset in floats from the base pointer.
	 */
	virtual int getVertexOffset() const
	{
		return m_vertexOffset;
	}

	/**
	 * Return the vertex stride in number of floats between vertices.
	 */
	virtual int getVertexStride() const
	{
		return m_vertexStride;
	}

	/**
	 * Return the vertex offset in floats from the base pointer.
	 */
	virtual int getNormalOffset() const
	{
		return m_normalOffset;
	}

	/**
	 * Return the vertex stride in number of floats between vertices.
	 */
	virtual int getNormalStride() const
	{
		return m_normalStride;
	}
};

class btCPUVertexBufferDescriptor : public btVertexBufferDescriptor
{
protected:
	float *m_basePointer;

public:
	/**
	 * vertexBasePointer is pointer to beginning of the buffer.
	 * vertexOffset is the offset in floats to the first vertex.
	 * vertexStride is the stride in floats between vertices.
	 */
	btCPUVertexBufferDescriptor(float *basePointer, int vertexOffset, int vertexStride)
	{
		m_basePointer = basePointer;
		m_vertexOffset = vertexOffset;
		m_vertexStride = vertexStride;
		m_hasVertexPositions = true;
	}

	/**
	 * vertexBasePointer is pointer to beginning of the buffer.
	 * vertexOffset is the offset in floats to the first vertex.
	 * vertexStride is the stride in floats between vertices.
	 */
	btCPUVertexBufferDescriptor(float *basePointer, int vertexOffset, int vertexStride, int normalOffset, int normalStride)
	{
		m_basePointer = basePointer;

		m_vertexOffset = vertexOffset;
		m_vertexStride = vertexStride;
		m_hasVertexPositions = true;

		m_normalOffset = normalOffset;
		m_normalStride = normalStride;
		m_hasNormals = true;
	}

	virtual ~btCPUVertexBufferDescriptor()
	{
	}

	/**
	 * Return the type of the vertex buffer descriptor.
	 */
	virtual BufferTypes getBufferType() const
	{
		return CPU_BUFFER;
	}

	/**
	 * Return the base pointer in memory to the first vertex.
	 */
	virtual float *getBasePointer() const
	{
		return m_basePointer;
	}
};

#endif  // #ifndef BT_SOFT_BODY_SOLVER_VERTEX_BUFFER_H
/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/


#include "btSoftMultiBodyDynamicsWorld.h"
#include "LinearMath/btQuickprof.h"

//softbody & helpers
#include "BulletSoftBody/btSoftBody.h"
#include "BulletSoftBody/btSoftBodyHelpers.h"
#include "BulletSoftBody/btSoftBodySolvers.h"
#include "BulletSoftBody/btDefaultSoftBodySolver.h"
#include "LinearMath/btSerializer.h"

btSoftMultiBodyDynamicsWorld::btSoftMultiBodyDynamicsWorld(
	btDispatcher* dispatcher,
	btBroadphaseInterface* pairCache,
	btMultiBodyConstraintSolver* constraintSolver,
	btCollisionConfiguration* collisionConfiguration,
	btSoftBodySolver* softBodySolver) : btMultiBodyDynamicsWorld(dispatcher, pairCache, constraintSolver, collisionConfiguration),
										m_softBodySolver(softBodySolver),
										m_ownsSolver(false)
{
	if (!m_softBodySolver)
	{
		void* ptr = btAlignedAlloc(sizeof(btDefaultSoftBodySolver), 16);
		m_softBodySolver = new (ptr) btDefaultSoftBodySolver();
		m_ownsSolver = true;
	}

	m_drawFlags = fDrawFlags::Std;
	m_drawNodeTree = true;
	m_drawFaceTree = false;
	m_drawClusterTree = false;
	m_sbi.m_broadphase = pairCache;
	m_sbi.m_dispatcher = dispatcher;
	m_sbi.m_sparsesdf.Initialize();
	m_sbi.m_sparsesdf.Reset();

	m_sbi.air_density = (btScalar)1.2;
	m_sbi.water_density = 0;
	m_sbi.water_offset = 0;
	m_sbi.water_normal = btVector3(0, 0, 0);
	m_sbi.m_gravity.setValue(0, -10, 0);

	m_sbi.m_sparsesdf.Initialize();
}

btSoftMultiBodyDynamicsWorld::~btSoftMultiBodyDynamicsWorld()
{
	if (m_ownsSolver)
	{
		m_softBodySolver->~btSoftBodySolver();
		btAlignedFree(m_softBodySolver);
	}
}

void btSoftMultiBodyDynamicsWorld::predictUnconstraintMotion(btScalar timeStep)
{
	btDiscreteDynamicsWorld::predictUnconstraintMotion(timeStep);
	{
		BT_PROFILE("predictUnconstraintMotionSoftBody");
		m_softBodySolver->predictMotion(float(timeStep));
	}
}

void btSoftMultiBodyDynamicsWorld::internalSingleStepSimulation(btScalar timeStep)
{
	// Let the solver grab the soft bodies and if necessary optimize for it
	m_softBodySolver->optimize(getSoftBodyArray());

	if (!m_softBodySolver->checkInitialized())
	{
		btAssert("Solver initialization failed\n");
	}

	btDiscreteDynamicsWorld::internalSingleStepSimulation(timeStep);

	///solve soft bodies constraints
	solveSoftBodiesConstraints(timeStep);

	//self collisions
	for (int i = 0; i < m_softBodies.size(); i++)
	{
		btSoftBody* psb = (btSoftBody*)m_softBodies[i];
		psb->defaultCollisionHandler(psb);
	}

	///update soft bodies
	m_softBodySolver->updateSoftBodies();

	for (int i = 0; i < m_softBodies.size(); i++)
	{
		btSoftBody* psb = (btSoftBody*)m_softBodies[i];
		psb->interpolateRenderMesh();
	}
	// End solver-wise simulation step
	// ///////////////////////////////
}

void btSoftMultiBodyDynamicsWorld::solveSoftBodiesConstraints(btScalar timeStep)
{
	BT_PROFILE("solveSoftConstraints");

	if (m_softBodies.size())
	{
		btSoftBody::solveClusters(m_softBodies);
	}

	// Solve constraints solver-wise
	m_softBodySolver->solveConstraints(timeStep * m_softBodySolver->getTimeScale());
}

void btSoftMultiBodyDynamicsWorld::addSoftBody(btSoftBody* body, int collisionFilterGroup, int collisionFilterMask)
{
	m_softBodies.push_back(body);

	// Set the soft body solver that will deal with this body
	// to be the world's solver
	body->setSoftBodySolver(m_softBodySolver);

	btCollisionWorld::addCollisionObject(body,
										 collisionFilterGroup,
										 collisionFilterMask);
}

void btSoftMultiBodyDynamicsWorld::removeSoftBody(btSoftBody* body)
{
	m_softBodies.remove(body);

	btCollisionWorld::removeCollisionObject(body);
}

void btSoftMultiBodyDynamicsWorld::removeCollisionObject(btCollisionObject* collisionObject)
{
	btSoftBody* body = btSoftBody::upcast(collisionObject);
	if (body)
		removeSoftBody(body);
	else
		btDiscreteDynamicsWorld::removeCollisionObject(collisionObject);
}

void btSoftMultiBodyDynamicsWorld::debugDrawWorld()
{
	btMultiBodyDynamicsWorld::debugDrawWorld();

	if (getDebugDrawer())
	{
		int i;
		for (i = 0; i < this->m_softBodies.size(); i++)
		{
			btSoftBody* psb = (btSoftBody*)this->m_softBodies[i];
			if (getDebugDrawer() && (getDebugDrawer()->getDebugMode() & (btIDebugDraw::DBG_DrawWireframe)))
			{
				btSoftBodyHelpers::DrawFrame(psb, m_debugDrawer);
				btSoftBodyHelpers::Draw(psb, m_debugDrawer, m_drawFlags);
			}

			if (m_debugDrawer && (m_debugDrawer->getDebugMode() & btIDebugDraw::DBG_DrawAabb))
			{
				if (m_drawNodeTree) btSoftBodyHelpers::DrawNodeTree(psb, m_debugDrawer);
				if (m_drawFaceTree) btSoftBodyHelpers::DrawFaceTree(psb, m_debugDrawer);
				if (m_drawClusterTree) btSoftBodyHelpers::DrawClusterTree(psb, m_debugDrawer);
			}
		}
	}
}

struct btSoftSingleRayCallback : public btBroadphaseRayCallback
{
	btVector3 m_rayFromWorld;
	btVector3 m_rayToWorld;
	btTransform m_rayFromTrans;
	btTransform m_rayToTrans;
	btVector3 m_hitNormal;

	const btSoftMultiBodyDynamicsWorld* m_world;
	btCollisionWorld::RayResultCallback& m_resultCallback;

	btSoftSingleRayCallback(const btVector3& rayFromWorld, const btVector3& rayToWorld, const btSoftMultiBodyDynamicsWorld* world, btCollisionWorld::RayResultCallback& resultCallback)
		: m_rayFromWorld(rayFromWorld),
		  m_rayToWorld(rayToWorld),
		  m_world(world),
		  m_resultCallback(resultCallback)
	{
		m_rayFromTrans.setIdentity();
		m_rayFromTrans.setOrigin(m_rayFromWorld);
		m_rayToTrans.setIdentity();
		m_rayToTrans.setOrigin(m_rayToWorld);

		btVector3 rayDir = (rayToWorld - rayFromWorld);

		rayDir.normalize();
		///what about division by zero? --> just set rayDirection[i] to INF/1e30
		m_rayDirectionInverse[0] = rayDir[0] == btScalar(0.0) ? btScalar(1e30) : btScalar(1.0) / rayDir[0];
		m_rayDirectionInverse[1] = rayDir[1] == btScalar(0.0) ? btScalar(1e30) : btScalar(1.0) / rayDir[1];
		m_rayDirectionInverse[2] = rayDir[2] == btScalar(0.0) ? btScalar(1e30) : btScalar(1.0) / rayDir[2];
		m_signs[0] = m_rayDirectionInverse[0] < 0.0;
		m_signs[1] = m_rayDirectionInverse[1] < 0.0;
		m_signs[2] = m_rayDirectionInverse[2] < 0.0;

		m_lambda_max = rayDir.dot(m_rayToWorld - m_rayFromWorld);
	}

	virtual bool process(const btBroadphaseProxy* proxy)
	{
		///terminate further ray tests, once the closestHitFraction reached zero
		if (m_resultCallback.m_closestHitFraction == btScalar(0.f))
			return false;

		btCollisionObject* collisionObject = (btCollisionObject*)proxy->m_clientObject;

		//only perform raycast if filterMask matches
		if (m_resultCallback.needsCollision(collisionObject->getBroadphaseHandle()))
		{
			//RigidcollisionObject* collisionObject = ctrl->GetRigidcollisionObject();
			//btVector3 collisionObjectAabbMin,collisionObjectAabbMax;
#if 0
#ifdef RECALCULATE_AABB
			btVector3 collisionObjectAabbMin,collisionObjectAabbMax;
			collisionObject->getCollisionShape()->getAabb(collisionObject->getWorldTransform(),collisionObjectAabbMin,collisionObjectAabbMax);
#else
			//getBroadphase()->getAabb(collisionObject->getBroadphaseHandle(),collisionObjectAabbMin,collisionObjectAabbMax);
			const btVector3& collisionObjectAabbMin = collisionObject->getBroadphaseHandle()->m_aabbMin;
			const btVector3& collisionObjectAabbMax = collisionObject->getBroadphaseHandle()->m_aabbMax;
#endif
#endif
			//btScalar hitLambda = m_resultCallback.m_closestHitFraction;
			//culling already done by broadphase
			//if (btRayAabb(m_rayFromWorld,m_rayToWorld,collisionObjectAabbMin,collisionObjectAabbMax,hitLambda,m_hitNormal))
			{
				m_world->rayTestSingle(m_rayFromTrans, m_rayToTrans,
									   collisionObject,
									   collisionObject->getCollisionShape(),
									   collisionObject->getWorldTransform(),
									   m_resultCallback);
			}
		}
		return true;
	}
};

void btSoftMultiBodyDynamicsWorld::rayTest(const btVector3& rayFromWorld, const btVector3& rayToWorld, RayResultCallback& resultCallback) const
{
	BT_PROFILE("rayTest");
	/// use the broadphase to accelerate the search for objects, based on their aabb
	/// and for each object with ray-aabb overlap, perform an exact ray test
	btSoftSingleRayCallback rayCB(rayFromWorld, rayToWorld, this, resultCallback);

#ifndef USE_BRUTEFORCE_RAYBROADPHASE
	m_broadphasePairCache->rayTest(rayFromWorld, rayToWorld, rayCB);
#else
	for (int i = 0; i < this->getNumCollisionObjects(); i++)
	{
		rayCB.process(m_collisionObjects[i]->getBroadphaseHandle());
	}
#endif  //USE_BRUTEFORCE_RAYBROADPHASE
}

void btSoftMultiBodyDynamicsWorld::rayTestSingle(const btTransform& rayFromTrans, const btTransform& rayToTrans,
												 btCollisionObject* collisionObject,
												 const btCollisionShape* collisionShape,
												 const btTransform& colObjWorldTransform,
												 RayResultCallback& resultCallback)
{
	if (collisionShape->isSoftBody())
	{
		btSoftBody* softBody = btSoftBody::upcast(collisionObject);
		if (softBody)
		{
			btSoftBody::sRayCast softResult;
			if (softBody->rayTest(rayFromTrans.getOrigin(), rayToTrans.getOrigin(), softResult))
			{
				if (softResult.fraction <= resultCallback.m_closestHitFraction)
				{
					btCollisionWorld::LocalShapeInfo shapeInfo;
					shapeInfo.m_shapePart = 0;
					shapeInfo.m_triangleIndex = softResult.index;
					// get the normal
					btVector3 rayDir = rayToTrans.getOrigin() - rayFromTrans.getOrigin();
					btVector3 normal = -rayDir;
					normal.normalize();

					if (softResult.feature == btSoftBody::eFeature::Face)
					{
						normal = softBody->m_faces[softResult.index].m_normal;
						if (normal.dot(rayDir) > 0)
						{
							// normal always point toward origin of the ray
							normal = -normal;
						}
					}

					btCollisionWorld::LocalRayResult rayResult(collisionObject,
															   &shapeInfo,
															   normal,
															   softResult.fraction);
					bool normalInWorldSpace = true;
					resultCallback.addSingleResult(rayResult, normalInWorldSpace);
				}
			}
		}
	}
	else
	{
		btCollisionWorld::rayTestSingle(rayFromTrans, rayToTrans, collisionObject, collisionShape, colObjWorldTransform, resultCallback);
	}
}

void btSoftMultiBodyDynamicsWorld::serializeSoftBodies(btSerializer* serializer)
{
	int i;
	//serialize all collision objects
	for (i = 0; i < m_collisionObjects.size(); i++)
	{
		btCollisionObject* colObj = m_collisionObjects[i];
		if (colObj->getInternalType() & btCollisionObject::CO_SOFT_BODY)
		{
			int len = colObj->calculateSerializeBufferSize();
			btChunk* chunk = serializer->allocate(len, 1);
			const char* structType = colObj->serialize(chunk->m_oldPtr, serializer);
			serializer->finalizeChunk(chunk, structType, BT_SOFTBODY_CODE, colObj);
		}
	}
}

void btSoftMultiBodyDynamicsWorld::serialize(btSerializer* serializer)
{
	serializer->startSerialization();

	serializeDynamicsWorldInfo(serializer);

	serializeSoftBodies(serializer);

	serializeMultiBodies(serializer);

	serializeRigidBodies(serializer);

	serializeCollisionObjects(serializer);

	serializeContactManifolds(serializer);

	serializer->finishSerialization();
}
/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/


#ifndef BT_SOFT_MULTIBODY_DYNAMICS_WORLD_H
#define BT_SOFT_MULTIBODY_DYNAMICS_WORLD_H

#include "BulletDynamics/Dynamics/btDiscreteDynamicsWorld.h"
#include "BulletDynamics/Featherstone/btMultiBodyDynamicsWorld.h"
#include "BulletSoftBody/btSoftBody.h"

#ifndef BT_SOFT_RIGID_DYNAMICS_WORLD_H
typedef btAlignedObjectArray<btSoftBody*> btSoftBodyArray;
#endif

class btSoftBodySolver;

class btSoftMultiBodyDynamicsWorld : public btMultiBodyDynamicsWorld
{
	btSoftBodyArray m_softBodies;
	int m_drawFlags;
	bool m_drawNodeTree;
	bool m_drawFaceTree;
	bool m_drawClusterTree;
	btSoftBodyWorldInfo m_sbi;
	///Solver classes that encapsulate multiple soft bodies for solving
	btSoftBodySolver* m_softBodySolver;
	bool m_ownsSolver;

protected:
	virtual void predictUnconstraintMotion(btScalar timeStep);

	virtual void internalSingleStepSimulation(btScalar timeStep);

	void solveSoftBodiesConstraints(btScalar timeStep);

	void serializeSoftBodies(btSerializer* serializer);

public:
	btSoftMultiBodyDynamicsWorld(btDispatcher* dispatcher, btBroadphaseInterface* pairCache, btMultiBodyConstraintSolver* constraintSolver, btCollisionConfiguration* collisionConfiguration, btSoftBodySolver* softBodySolver = 0);

	virtual ~btSoftMultiBodyDynamicsWorld();

	virtual void debugDrawWorld();

	void addSoftBody(btSoftBody* body, int collisionFilterGroup = btBroadphaseProxy::DefaultFilter, int collisionFilterMask = btBroadphaseProxy::AllFilter);

	void removeSoftBody(btSoftBody* body);

	///removeCollisionObject will first check if it is a rigid body, if so call removeRigidBody otherwise call btDiscreteDynamicsWorld::removeCollisionObject
	virtual void removeCollisionObject(btCollisionObject* collisionObject);

	int getDrawFlags() const { return (m_drawFlags); }
	void setDrawFlags(int f) { m_drawFlags = f; }

	btSoftBodyWorldInfo& getWorldInfo()
	{
		return m_sbi;
	}
	const btSoftBodyWorldInfo& getWorldInfo() const
	{
		return m_sbi;
	}

	virtual btDynamicsWorldType getWorldType() const
	{
		return BT_SOFT_MULTIBODY_DYNAMICS_WORLD;
	}

	btSoftBodyArray& getSoftBodyArray()
	{
		return m_softBodies;
	}

	const btSoftBodyArray& getSoftBodyArray() const
	{
		return m_softBodies;
	}

	virtual void rayTest(const btVector3& rayFromWorld, const btVector3& rayToWorld, RayResultCallback& resultCallback) const;

	/// rayTestSingle performs a raycast call and calls the resultCallback. It is used internally by rayTest.
	/// In a future implementation, we consider moving the ray test as a virtual method in btCollisionShape.
	/// This allows more customization.
	static void rayTestSingle(const btTransform& rayFromTrans, const btTransform& rayToTrans,
							  btCollisionObject* collisionObject,
							  const btCollisionShape* collisionShape,
							  const btTransform& colObjWorldTransform,
							  RayResultCallback& resultCallback);

	virtual void serialize(btSerializer* serializer);
};

#endif  //BT_SOFT_MULTIBODY_DYNAMICS_WORLD_H
/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/

#include "btSoftRigidCollisionAlgorithm.h"
#include "BulletCollision/CollisionDispatch/btCollisionDispatcher.h"
#include "BulletCollision/CollisionShapes/btSphereShape.h"
#include "BulletCollision/CollisionShapes/btBoxShape.h"
#include "BulletCollision/CollisionDispatch/btCollisionObject.h"
#include "btSoftBody.h"
#include "BulletSoftBody/btSoftBodySolvers.h"
#include "BulletCollision/CollisionDispatch/btCollisionObjectWrapper.h"

///TODO: include all the shapes that the softbody can collide with
///alternatively, implement special case collision algorithms (just like for rigid collision shapes)

//#include <stdio.h>

btSoftRigidCollisionAlgorithm::btSoftRigidCollisionAlgorithm(btPersistentManifold* /*mf*/, const btCollisionAlgorithmConstructionInfo& ci, const btCollisionObjectWrapper*, const btCollisionObjectWrapper*, bool isSwapped)
	: btCollisionAlgorithm(ci),
	  //m_ownManifold(false),
	  //m_manifoldPtr(mf),
	  m_isSwapped(isSwapped)
{
}

btSoftRigidCollisionAlgorithm::~btSoftRigidCollisionAlgorithm()
{
	//m_softBody->m_overlappingRigidBodies.remove(m_rigidCollisionObject);

	/*if (m_ownManifold)
	{
	if (m_manifoldPtr)
	m_dispatcher->releaseManifold(m_manifoldPtr);
	}
	*/
}

#include <stdio.h>
#include "LinearMath/btQuickprof.h"
void btSoftRigidCollisionAlgorithm::processCollision(const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
{
	BT_PROFILE("btSoftRigidCollisionAlgorithm::processCollision");
	(void)dispatchInfo;
	//printf("btSoftRigidCollisionAlgorithm\n");
	//	const btCollisionObjectWrapper* softWrap = m_isSwapped?body1Wrap:body0Wrap;
	//	const btCollisionObjectWrapper* rigidWrap = m_isSwapped?body0Wrap:body1Wrap;
	btSoftBody* softBody = m_isSwapped ? (btSoftBody*)body1Wrap->getCollisionObject() : (btSoftBody*)body0Wrap->getCollisionObject();
	const btCollisionObjectWrapper* rigidCollisionObjectWrap = m_isSwapped ? body0Wrap : body1Wrap;

	if (softBody->m_collisionDisabledObjects.findLinearSearch(rigidCollisionObjectWrap->getCollisionObject()) == softBody->m_collisionDisabledObjects.size())
	{
		softBody->getSoftBodySolver()->processCollision(softBody, rigidCollisionObjectWrap, static_cast<btManifoldResultForSkin*>(resultOut));
	}
}

btScalar btSoftRigidCollisionAlgorithm::calculateTimeOfImpact(btCollisionObject* col0, btCollisionObject* col1, const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
{
	(void)resultOut;
	(void)dispatchInfo;
	(void)col0;
	(void)col1;

	//not yet
	return btScalar(1.);
}
/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/


#ifndef BT_SOFT_RIGID_COLLISION_ALGORITHM_H
#define BT_SOFT_RIGID_COLLISION_ALGORITHM_H

#include "BulletCollision/BroadphaseCollision/btCollisionAlgorithm.h"
#include "BulletCollision/BroadphaseCollision/btBroadphaseProxy.h"
#include "BulletCollision/CollisionDispatch/btCollisionCreateFunc.h"
class btPersistentManifold;
#include "BulletCollision/CollisionDispatch/btCollisionDispatcher.h"

#include "LinearMath/btVector3.h"
class btSoftBody;

/// btSoftRigidCollisionAlgorithm  provides collision detection between btSoftBody and btRigidBody
class btSoftRigidCollisionAlgorithm : public btCollisionAlgorithm
{
	//	bool	m_ownManifold;
	//	btPersistentManifold*	m_manifoldPtr;

	//btSoftBody*				m_softBody;
	//btCollisionObject*		m_rigidCollisionObject;

	///for rigid versus soft (instead of soft versus rigid), we use this swapped boolean
	bool m_isSwapped;

public:
	btSoftRigidCollisionAlgorithm(btPersistentManifold* mf, const btCollisionAlgorithmConstructionInfo& ci, const btCollisionObjectWrapper* col0, const btCollisionObjectWrapper* col1Wrap, bool isSwapped);

	virtual ~btSoftRigidCollisionAlgorithm();

	virtual void processCollision(const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut);

	virtual btScalar calculateTimeOfImpact(btCollisionObject* body0, btCollisionObject* body1, const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut);

	virtual void getAllContactManifolds(btManifoldArray& manifoldArray)
	{
		//we don't add any manifolds
	}

	struct CreateFunc : public btCollisionAlgorithmCreateFunc
	{
		virtual btCollisionAlgorithm* CreateCollisionAlgorithm(btCollisionAlgorithmConstructionInfo& ci, const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap)
		{
			void* mem = ci.m_dispatcher1->allocateCollisionAlgorithm(sizeof(btSoftRigidCollisionAlgorithm));
			if (!m_swapped)
			{
				return new (mem) btSoftRigidCollisionAlgorithm(0, ci, body0Wrap, body1Wrap, false);
			}
			else
			{
				return new (mem) btSoftRigidCollisionAlgorithm(0, ci, body0Wrap, body1Wrap, true);
			}
		}
	};
};

#endif  //BT_SOFT_RIGID_COLLISION_ALGORITHM_H
/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/


#include "btSoftRigidDynamicsWorld.h"
#include "LinearMath/btQuickprof.h"

//softbody & helpers
#include "btSoftBody.h"
#include "btSoftBodyHelpers.h"
#include "btSoftBodySolvers.h"
#include "btDefaultSoftBodySolver.h"
#include "LinearMath/btSerializer.h"

btSoftRigidDynamicsWorld::btSoftRigidDynamicsWorld(
	btDispatcher* dispatcher,
	btBroadphaseInterface* pairCache,
	btConstraintSolver* constraintSolver,
	btCollisionConfiguration* collisionConfiguration,
	btSoftBodySolver* softBodySolver) : btDiscreteDynamicsWorld(dispatcher, pairCache, constraintSolver, collisionConfiguration),
										m_softBodySolver(softBodySolver),
										m_ownsSolver(false)
{
	if (!m_softBodySolver)
	{
		void* ptr = btAlignedAlloc(sizeof(btDefaultSoftBodySolver), 16);
		m_softBodySolver = new (ptr) btDefaultSoftBodySolver();
		m_ownsSolver = true;
	}

	m_drawFlags = fDrawFlags::Std;
	m_drawNodeTree = true;
	m_drawFaceTree = false;
	m_drawClusterTree = false;
	m_sbi.m_broadphase = pairCache;
	m_sbi.m_dispatcher = dispatcher;
	m_sbi.m_sparsesdf.Initialize();
	m_sbi.m_sparsesdf.Reset();

	m_sbi.air_density = (btScalar)1.2;
	m_sbi.water_density = 0;
	m_sbi.water_offset = 0;
	m_sbi.water_normal = btVector3(0, 0, 0);
	m_sbi.m_gravity.setValue(0, -10, 0);

	m_sbi.m_sparsesdf.Initialize();
}

btSoftRigidDynamicsWorld::~btSoftRigidDynamicsWorld()
{
	if (m_ownsSolver)
	{
		m_softBodySolver->~btSoftBodySolver();
		btAlignedFree(m_softBodySolver);
	}
}

void btSoftRigidDynamicsWorld::predictUnconstraintMotion(btScalar timeStep)
{
	btDiscreteDynamicsWorld::predictUnconstraintMotion(timeStep);
	{
		BT_PROFILE("predictUnconstraintMotionSoftBody");
		m_softBodySolver->predictMotion(float(timeStep));
	}
}

void btSoftRigidDynamicsWorld::internalSingleStepSimulation(btScalar timeStep)
{
	// Let the solver grab the soft bodies and if necessary optimize for it
	m_softBodySolver->optimize(getSoftBodyArray());

	if (!m_softBodySolver->checkInitialized())
	{
		btAssert("Solver initialization failed\n");
	}

	btDiscreteDynamicsWorld::internalSingleStepSimulation(timeStep);

	///solve soft bodies constraints
	solveSoftBodiesConstraints(timeStep);

	//self collisions
	for (int i = 0; i < m_softBodies.size(); i++)
	{
		btSoftBody* psb = (btSoftBody*)m_softBodies[i];
		psb->defaultCollisionHandler(psb);
	}

	///update soft bodies
	m_softBodySolver->updateSoftBodies();

	// End solver-wise simulation step
	// ///////////////////////////////
}

void btSoftRigidDynamicsWorld::solveSoftBodiesConstraints(btScalar timeStep)
{
	BT_PROFILE("solveSoftConstraints");

	if (m_softBodies.size())
	{
		btSoftBody::solveClusters(m_softBodies);
	}

	// Solve constraints solver-wise
	m_softBodySolver->solveConstraints(timeStep * m_softBodySolver->getTimeScale());
}

void btSoftRigidDynamicsWorld::addSoftBody(btSoftBody* body, int collisionFilterGroup, int collisionFilterMask)
{
	m_softBodies.push_back(body);

	// Set the soft body solver that will deal with this body
	// to be the world's solver
	body->setSoftBodySolver(m_softBodySolver);

	btCollisionWorld::addCollisionObject(body,
										 collisionFilterGroup,
										 collisionFilterMask);
}

void btSoftRigidDynamicsWorld::removeSoftBody(btSoftBody* body)
{
	m_softBodies.remove(body);

	btCollisionWorld::removeCollisionObject(body);
}

void btSoftRigidDynamicsWorld::removeCollisionObject(btCollisionObject* collisionObject)
{
	btSoftBody* body = btSoftBody::upcast(collisionObject);
	if (body)
		removeSoftBody(body);
	else
		btDiscreteDynamicsWorld::removeCollisionObject(collisionObject);
}

void btSoftRigidDynamicsWorld::debugDrawWorld()
{
	btDiscreteDynamicsWorld::debugDrawWorld();

	if (getDebugDrawer())
	{
		int i;
		for (i = 0; i < this->m_softBodies.size(); i++)
		{
			btSoftBody* psb = (btSoftBody*)this->m_softBodies[i];
			if (getDebugDrawer() && (getDebugDrawer()->getDebugMode() & (btIDebugDraw::DBG_DrawWireframe)))
			{
				btSoftBodyHelpers::DrawFrame(psb, m_debugDrawer);
				btSoftBodyHelpers::Draw(psb, m_debugDrawer, m_drawFlags);
			}

			if (m_debugDrawer && (m_debugDrawer->getDebugMode() & btIDebugDraw::DBG_DrawAabb))
			{
				if (m_drawNodeTree) btSoftBodyHelpers::DrawNodeTree(psb, m_debugDrawer);
				if (m_drawFaceTree) btSoftBodyHelpers::DrawFaceTree(psb, m_debugDrawer);
				if (m_drawClusterTree) btSoftBodyHelpers::DrawClusterTree(psb, m_debugDrawer);
			}
		}
	}
}

struct btSoftSingleRayCallback : public btBroadphaseRayCallback
{
	btVector3 m_rayFromWorld;
	btVector3 m_rayToWorld;
	btTransform m_rayFromTrans;
	btTransform m_rayToTrans;
	btVector3 m_hitNormal;

	const btSoftRigidDynamicsWorld* m_world;
	btCollisionWorld::RayResultCallback& m_resultCallback;

	btSoftSingleRayCallback(const btVector3& rayFromWorld, const btVector3& rayToWorld, const btSoftRigidDynamicsWorld* world, btCollisionWorld::RayResultCallback& resultCallback)
		: m_rayFromWorld(rayFromWorld),
		  m_rayToWorld(rayToWorld),
		  m_world(world),
		  m_resultCallback(resultCallback)
	{
		m_rayFromTrans.setIdentity();
		m_rayFromTrans.setOrigin(m_rayFromWorld);
		m_rayToTrans.setIdentity();
		m_rayToTrans.setOrigin(m_rayToWorld);

		btVector3 rayDir = (rayToWorld - rayFromWorld);

		rayDir.normalize();
		///what about division by zero? --> just set rayDirection[i] to INF/1e30
		m_rayDirectionInverse[0] = rayDir[0] == btScalar(0.0) ? btScalar(1e30) : btScalar(1.0) / rayDir[0];
		m_rayDirectionInverse[1] = rayDir[1] == btScalar(0.0) ? btScalar(1e30) : btScalar(1.0) / rayDir[1];
		m_rayDirectionInverse[2] = rayDir[2] == btScalar(0.0) ? btScalar(1e30) : btScalar(1.0) / rayDir[2];
		m_signs[0] = m_rayDirectionInverse[0] < 0.0;
		m_signs[1] = m_rayDirectionInverse[1] < 0.0;
		m_signs[2] = m_rayDirectionInverse[2] < 0.0;

		m_lambda_max = rayDir.dot(m_rayToWorld - m_rayFromWorld);
	}

	virtual bool process(const btBroadphaseProxy* proxy)
	{
		///terminate further ray tests, once the closestHitFraction reached zero
		if (m_resultCallback.m_closestHitFraction == btScalar(0.f))
			return false;

		btCollisionObject* collisionObject = (btCollisionObject*)proxy->m_clientObject;

		//only perform raycast if filterMask matches
		if (m_resultCallback.needsCollision(collisionObject->getBroadphaseHandle()))
		{
			//RigidcollisionObject* collisionObject = ctrl->GetRigidcollisionObject();
			//btVector3 collisionObjectAabbMin,collisionObjectAabbMax;
#if 0
#ifdef RECALCULATE_AABB
			btVector3 collisionObjectAabbMin,collisionObjectAabbMax;
			collisionObject->getCollisionShape()->getAabb(collisionObject->getWorldTransform(),collisionObjectAabbMin,collisionObjectAabbMax);
#else
			//getBroadphase()->getAabb(collisionObject->getBroadphaseHandle(),collisionObjectAabbMin,collisionObjectAabbMax);
			const btVector3& collisionObjectAabbMin = collisionObject->getBroadphaseHandle()->m_aabbMin;
			const btVector3& collisionObjectAabbMax = collisionObject->getBroadphaseHandle()->m_aabbMax;
#endif
#endif
			//btScalar hitLambda = m_resultCallback.m_closestHitFraction;
			//culling already done by broadphase
			//if (btRayAabb(m_rayFromWorld,m_rayToWorld,collisionObjectAabbMin,collisionObjectAabbMax,hitLambda,m_hitNormal))
			{
				m_world->rayTestSingle(m_rayFromTrans, m_rayToTrans,
									   collisionObject,
									   collisionObject->getCollisionShape(),
									   collisionObject->getWorldTransform(),
									   m_resultCallback);
			}
		}
		return true;
	}
};

void btSoftRigidDynamicsWorld::rayTest(const btVector3& rayFromWorld, const btVector3& rayToWorld, RayResultCallback& resultCallback) const
{
	BT_PROFILE("rayTest");
	/// use the broadphase to accelerate the search for objects, based on their aabb
	/// and for each object with ray-aabb overlap, perform an exact ray test
	btSoftSingleRayCallback rayCB(rayFromWorld, rayToWorld, this, resultCallback);

#ifndef USE_BRUTEFORCE_RAYBROADPHASE
	m_broadphasePairCache->rayTest(rayFromWorld, rayToWorld, rayCB);
#else
	for (int i = 0; i < this->getNumCollisionObjects(); i++)
	{
		rayCB.process(m_collisionObjects[i]->getBroadphaseHandle());
	}
#endif  //USE_BRUTEFORCE_RAYBROADPHASE
}

void btSoftRigidDynamicsWorld::rayTestSingle(const btTransform& rayFromTrans, const btTransform& rayToTrans,
											 btCollisionObject* collisionObject,
											 const btCollisionShape* collisionShape,
											 const btTransform& colObjWorldTransform,
											 RayResultCallback& resultCallback)
{
	if (collisionShape->isSoftBody())
	{
		btSoftBody* softBody = btSoftBody::upcast(collisionObject);
		if (softBody)
		{
			btSoftBody::sRayCast softResult;
			if (softBody->rayTest(rayFromTrans.getOrigin(), rayToTrans.getOrigin(), softResult))
			{
				if (softResult.fraction <= resultCallback.m_closestHitFraction)
				{
					btCollisionWorld::LocalShapeInfo shapeInfo;
					shapeInfo.m_shapePart = 0;
					shapeInfo.m_triangleIndex = softResult.index;
					// get the normal
					btVector3 rayDir = rayToTrans.getOrigin() - rayFromTrans.getOrigin();
					btVector3 normal = -rayDir;
					normal.normalize();

					if (softResult.feature == btSoftBody::eFeature::Face)
					{
						normal = softBody->m_faces[softResult.index].m_normal;
						if (normal.dot(rayDir) > 0)
						{
							// normal always point toward origin of the ray
							normal = -normal;
						}
					}

					btCollisionWorld::LocalRayResult rayResult(collisionObject,
															   &shapeInfo,
															   normal,
															   softResult.fraction);
					bool normalInWorldSpace = true;
					resultCallback.addSingleResult(rayResult, normalInWorldSpace);
				}
			}
		}
	}
	else
	{
		btCollisionWorld::rayTestSingle(rayFromTrans, rayToTrans, collisionObject, collisionShape, colObjWorldTransform, resultCallback);
	}
}

void btSoftRigidDynamicsWorld::serializeSoftBodies(btSerializer* serializer)
{
	int i;
	//serialize all collision objects
	for (i = 0; i < m_collisionObjects.size(); i++)
	{
		btCollisionObject* colObj = m_collisionObjects[i];
		if (colObj->getInternalType() & btCollisionObject::CO_SOFT_BODY)
		{
			int len = colObj->calculateSerializeBufferSize();
			btChunk* chunk = serializer->allocate(len, 1);
			const char* structType = colObj->serialize(chunk->m_oldPtr, serializer);
			serializer->finalizeChunk(chunk, structType, BT_SOFTBODY_CODE, colObj);
		}
	}
}

void btSoftRigidDynamicsWorld::serialize(btSerializer* serializer)
{
	serializer->startSerialization();

	serializeDynamicsWorldInfo(serializer);

	serializeSoftBodies(serializer);

	serializeRigidBodies(serializer);

	serializeCollisionObjects(serializer);

	serializer->finishSerialization();
}
/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/


#ifndef BT_SOFT_RIGID_DYNAMICS_WORLD_H
#define BT_SOFT_RIGID_DYNAMICS_WORLD_H

#include "BulletDynamics/Dynamics/btDiscreteDynamicsWorld.h"
#include "btSoftBody.h"

typedef btAlignedObjectArray<btSoftBody*> btSoftBodyArray;

class btSoftBodySolver;

class btSoftRigidDynamicsWorld : public btDiscreteDynamicsWorld
{
	btSoftBodyArray m_softBodies;
	int m_drawFlags;
	bool m_drawNodeTree;
	bool m_drawFaceTree;
	bool m_drawClusterTree;
	btSoftBodyWorldInfo m_sbi;
	///Solver classes that encapsulate multiple soft bodies for solving
	btSoftBodySolver* m_softBodySolver;
	bool m_ownsSolver;

protected:
	virtual void predictUnconstraintMotion(btScalar timeStep);

	virtual void internalSingleStepSimulation(btScalar timeStep);

	void solveSoftBodiesConstraints(btScalar timeStep);

	void serializeSoftBodies(btSerializer* serializer);

public:
	btSoftRigidDynamicsWorld(btDispatcher* dispatcher, btBroadphaseInterface* pairCache, btConstraintSolver* constraintSolver, btCollisionConfiguration* collisionConfiguration, btSoftBodySolver* softBodySolver = 0);

	virtual ~btSoftRigidDynamicsWorld();

	virtual void debugDrawWorld();

	void addSoftBody(btSoftBody* body, int collisionFilterGroup = btBroadphaseProxy::DefaultFilter, int collisionFilterMask = btBroadphaseProxy::AllFilter);

	void removeSoftBody(btSoftBody* body);

	///removeCollisionObject will first check if it is a rigid body, if so call removeRigidBody otherwise call btDiscreteDynamicsWorld::removeCollisionObject
	virtual void removeCollisionObject(btCollisionObject* collisionObject);

	int getDrawFlags() const { return (m_drawFlags); }
	void setDrawFlags(int f) { m_drawFlags = f; }

	btSoftBodyWorldInfo& getWorldInfo()
	{
		return m_sbi;
	}
	const btSoftBodyWorldInfo& getWorldInfo() const
	{
		return m_sbi;
	}

	virtual btDynamicsWorldType getWorldType() const
	{
		return BT_SOFT_RIGID_DYNAMICS_WORLD;
	}

	btSoftBodyArray& getSoftBodyArray()
	{
		return m_softBodies;
	}

	const btSoftBodyArray& getSoftBodyArray() const
	{
		return m_softBodies;
	}

	virtual void rayTest(const btVector3& rayFromWorld, const btVector3& rayToWorld, RayResultCallback& resultCallback) const;

	/// rayTestSingle performs a raycast call and calls the resultCallback. It is used internally by rayTest.
	/// In a future implementation, we consider moving the ray test as a virtual method in btCollisionShape.
	/// This allows more customization.
	static void rayTestSingle(const btTransform& rayFromTrans, const btTransform& rayToTrans,
							  btCollisionObject* collisionObject,
							  const btCollisionShape* collisionShape,
							  const btTransform& colObjWorldTransform,
							  RayResultCallback& resultCallback);

	virtual void serialize(btSerializer* serializer);
};

#endif  //BT_SOFT_RIGID_DYNAMICS_WORLD_H
/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/

#include "btSoftSoftCollisionAlgorithm.h"
#include "BulletCollision/CollisionDispatch/btCollisionDispatcher.h"
#include "BulletCollision/CollisionShapes/btBoxShape.h"
#include "BulletCollision/CollisionDispatch/btCollisionObject.h"
#include "BulletSoftBody/btSoftBodySolvers.h"
#include "btSoftBody.h"
#include "BulletCollision/CollisionDispatch/btCollisionObjectWrapper.h"

#define USE_PERSISTENT_CONTACTS 1

btSoftSoftCollisionAlgorithm::btSoftSoftCollisionAlgorithm(btPersistentManifold* /*mf*/, const btCollisionAlgorithmConstructionInfo& ci, const btCollisionObjectWrapper* /*obj0*/, const btCollisionObjectWrapper* /*obj1*/)
	: btCollisionAlgorithm(ci)
//m_ownManifold(false),
//m_manifoldPtr(mf)
{
}

btSoftSoftCollisionAlgorithm::~btSoftSoftCollisionAlgorithm()
{
}

void btSoftSoftCollisionAlgorithm::processCollision(const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, const btDispatcherInfo& /*dispatchInfo*/, btManifoldResult* resultOut)
{
	btSoftBody* soft0 = (btSoftBody*)body0Wrap->getCollisionObject();
	btSoftBody* soft1 = (btSoftBody*)body1Wrap->getCollisionObject();
	if ((soft0->getCollisionShape()->getShapeType() == SOFTBODY_SHAPE_PROXYTYPE && soft1->getCollisionShape()->getShapeType() == SOFTBODY_SHAPE_PROXYTYPE))
		soft0->getSoftBodySolver()->processCollision(soft0, soft1);
	else
		soft0->getSoftBodySolver()->processCollision(soft0, soft1, static_cast<btManifoldResultForSkin*>(resultOut));
}

btScalar btSoftSoftCollisionAlgorithm::calculateTimeOfImpact(btCollisionObject* /*body0*/, btCollisionObject* /*body1*/, const btDispatcherInfo& /*dispatchInfo*/, btManifoldResult* /*resultOut*/)
{
	//not yet
	return 1.f;
}
/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/


#ifndef BT_SOFT_SOFT_COLLISION_ALGORITHM_H
#define BT_SOFT_SOFT_COLLISION_ALGORITHM_H

#include "BulletCollision/BroadphaseCollision/btCollisionAlgorithm.h"
#include "BulletCollision/BroadphaseCollision/btBroadphaseProxy.h"
#include "BulletCollision/BroadphaseCollision/btDispatcher.h"
#include "BulletCollision/CollisionDispatch/btCollisionCreateFunc.h"

class btPersistentManifold;
class btSoftBody;

///collision detection between two btSoftBody shapes
class btSoftSoftCollisionAlgorithm : public btCollisionAlgorithm
{
	bool m_ownManifold;
	btPersistentManifold* m_manifoldPtr;

	//	btSoftBody*	m_softBody0;
	//	btSoftBody*	m_softBody1;

public:
	btSoftSoftCollisionAlgorithm(const btCollisionAlgorithmConstructionInfo& ci)
		: btCollisionAlgorithm(ci) {}

	virtual void processCollision(const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap, const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut);

	virtual btScalar calculateTimeOfImpact(btCollisionObject* body0, btCollisionObject* body1, const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut);

	virtual void getAllContactManifolds(btManifoldArray& manifoldArray)
	{
		if (m_manifoldPtr && m_ownManifold)
			manifoldArray.push_back(m_manifoldPtr);
	}

	btSoftSoftCollisionAlgorithm(btPersistentManifold* mf, const btCollisionAlgorithmConstructionInfo& ci, const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap);

	virtual ~btSoftSoftCollisionAlgorithm();

	struct CreateFunc : public btCollisionAlgorithmCreateFunc
	{
		virtual btCollisionAlgorithm* CreateCollisionAlgorithm(btCollisionAlgorithmConstructionInfo& ci, const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap)
		{
			int bbsize = sizeof(btSoftSoftCollisionAlgorithm);
			void* ptr = ci.m_dispatcher1->allocateCollisionAlgorithm(bbsize);
			return new (ptr) btSoftSoftCollisionAlgorithm(0, ci, body0Wrap, body1Wrap);
		}
	};
};

#endif  //BT_SOFT_SOFT_COLLISION_ALGORITHM_H
/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/*
This is a modified version of the Bullet Continuous Collision Detection and Physics Library
*/

///btSparseSdf implementation by Nathanael Presson

#ifndef BT_SPARSE_SDF_H
#define BT_SPARSE_SDF_H

#include "BulletCollision/CollisionDispatch/btCollisionObject.h"
#include "BulletCollision/NarrowPhaseCollision/btGjkEpa2.h"

// Fast Hash

#if !defined(get16bits)
#define get16bits(d) ((((unsigned int)(((const unsigned char*)(d))[1])) << 8) + (unsigned int)(((const unsigned char*)(d))[0]))
#endif
//
// super hash function by Paul Hsieh
//
inline unsigned int HsiehHash(const char* data, int len)
{
	unsigned int hash = len, tmp;
	len >>= 2;

	/* Main loop */
	for (; len > 0; len--)
	{
		hash += get16bits(data);
		tmp = (get16bits(data + 2) << 11) ^ hash;
		hash = (hash << 16) ^ tmp;
		data += 2 * sizeof(unsigned short);
		hash += hash >> 11;
	}

	/* Force "avalanching" of final 127 bits */
	hash ^= hash << 3;
	hash += hash >> 5;
	hash ^= hash << 4;
	hash += hash >> 17;
	hash ^= hash << 25;
	hash += hash >> 6;

	return hash;
}

template <const int CELLSIZE>
struct btSparseSdf
{
	//
	// Inner types
	//
	struct IntFrac
	{
		int b;
		int i;
		btScalar f;
	};
	struct Cell
	{
		btScalar d[CELLSIZE + 1][CELLSIZE + 1][CELLSIZE + 1];
		int c[3];
		int puid;
		unsigned hash;
		const btCollisionShape* pclient;
		Cell* next;
	};
	//
	// Fields
	//

	btAlignedObjectArray<Cell*> cells;
	btScalar voxelsz;
	btScalar m_defaultVoxelsz;
	int puid;
	int ncells;
	int m_clampCells;
	int nprobes;
	int nqueries;

	~btSparseSdf()
	{
		Reset();
	}
	//
	// Methods
	//

	//
	void Initialize(int hashsize = 2383, int clampCells = 256 * 1024)
	{
		//avoid a crash due to running out of memory, so clamp the maximum number of cells allocated
		//if this limit is reached, the SDF is reset (at the cost of some performance during the reset)
		m_clampCells = clampCells;
		cells.resize(hashsize, 0);
		m_defaultVoxelsz = 0.25;
		Reset();
	}
	//

	void setDefaultVoxelsz(btScalar sz)
	{
		m_defaultVoxelsz = sz;
	}

	void Reset()
	{
		for (int i = 0, ni = cells.size(); i < ni; ++i)
		{
			Cell* pc = cells[i];
			cells[i] = 0;
			while (pc)
			{
				Cell* pn = pc->next;
				delete pc;
				pc = pn;
			}
		}
		voxelsz = m_defaultVoxelsz;
		puid = 0;
		ncells = 0;
		nprobes = 1;
		nqueries = 1;
	}
	//
	void GarbageCollect(int lifetime = 256)
	{
		const int life = puid - lifetime;
		for (int i = 0; i < cells.size(); ++i)
		{
			Cell*& root = cells[i];
			Cell* pp = 0;
			Cell* pc = root;
			while (pc)
			{
				Cell* pn = pc->next;
				if (pc->puid < life)
				{
					if (pp)
						pp->next = pn;
					else
						root = pn;
					delete pc;
					pc = pp;
					--ncells;
				}
				pp = pc;
				pc = pn;
			}
		}
		//printf("GC[%d]: %d cells, PpQ: %f\r\n",puid,ncells,nprobes/(btScalar)nqueries);
		nqueries = 1;
		nprobes = 1;
		++puid;  ///@todo: Reset puid's when int range limit is reached	*/
		/* else setup a priority list...						*/
	}
	//
	int RemoveReferences(btCollisionShape* pcs)
	{
		int refcount = 0;
		for (int i = 0; i < cells.size(); ++i)
		{
			Cell*& root = cells[i];
			Cell* pp = 0;
			Cell* pc = root;
			while (pc)
			{
				Cell* pn = pc->next;
				if (pc->pclient == pcs)
				{
					if (pp)
						pp->next = pn;
					else
						root = pn;
					delete pc;
					pc = pp;
					++refcount;
				}
				pp = pc;
				pc = pn;
			}
		}
		return (refcount);
	}
	//
	btScalar Evaluate(const btVector3& x,
					  const btCollisionShape* shape,
					  btVector3& normal,
					  btScalar margin)
	{
		/* Lookup cell			*/
		const btVector3 scx = x / voxelsz;
		const IntFrac ix = Decompose(scx.x());
		const IntFrac iy = Decompose(scx.y());
		const IntFrac iz = Decompose(scx.z());
		const unsigned h = Hash(ix.b, iy.b, iz.b, shape);
		Cell*& root = cells[static_cast<int>(h % cells.size())];
		Cell* c = root;
		++nqueries;
		while (c)
		{
			++nprobes;
			if ((c->hash == h) &&
				(c->c[0] == ix.b) &&
				(c->c[1] == iy.b) &&
				(c->c[2] == iz.b) &&
				(c->pclient == shape))
			{
				break;
			}
			else
			{
				// printf("c->hash/c[0][1][2]=%d,%d,%d,%d\n", c->hash, c->c[0], c->c[1],c->c[2]);
				//printf("h,ixb,iyb,izb=%d,%d,%d,%d\n", h,ix.b, iy.b, iz.b);

				c = c->next;
			}
		}
		if (!c)
		{
			++nprobes;
			++ncells;
			//int sz = sizeof(Cell);
			if (ncells > m_clampCells)
			{
				//static int numResets = 0;
				//numResets++;
				//printf("numResets=%d\n",numResets);
				Reset();
			}

			c = new Cell();
			c->next = root;
			root = c;
			c->pclient = shape;
			c->hash = h;
			c->c[0] = ix.b;
			c->c[1] = iy.b;
			c->c[2] = iz.b;
			BuildCell(*c);
		}
		c->puid = puid;
		/* Extract infos		*/
		const int o[] = {ix.i, iy.i, iz.i};
		const btScalar d[] = {c->d[o[0] + 0][o[1] + 0][o[2] + 0],
							  c->d[o[0] + 1][o[1] + 0][o[2] + 0],
							  c->d[o[0] + 1][o[1] + 1][o[2] + 0],
							  c->d[o[0] + 0][o[1] + 1][o[2] + 0],
							  c->d[o[0] + 0][o[1] + 0][o[2] + 1],
							  c->d[o[0] + 1][o[1] + 0][o[2] + 1],
							  c->d[o[0] + 1][o[1] + 1][o[2] + 1],
							  c->d[o[0] + 0][o[1] + 1][o[2] + 1]};
		/* Normal	*/
#if 1
		const btScalar gx[] = {d[1] - d[0], d[2] - d[3],
							   d[5] - d[4], d[6] - d[7]};
		const btScalar gy[] = {d[3] - d[0], d[2] - d[1],
							   d[7] - d[4], d[6] - d[5]};
		const btScalar gz[] = {d[4] - d[0], d[5] - d[1],
							   d[7] - d[3], d[6] - d[2]};
		normal.setX(Lerp(Lerp(gx[0], gx[1], iy.f),
						 Lerp(gx[2], gx[3], iy.f), iz.f));
		normal.setY(Lerp(Lerp(gy[0], gy[1], ix.f),
						 Lerp(gy[2], gy[3], ix.f), iz.f));
		normal.setZ(Lerp(Lerp(gz[0], gz[1], ix.f),
						 Lerp(gz[2], gz[3], ix.f), iy.f));
		normal.safeNormalize();
#else
		normal = btVector3(d[1] - d[0], d[3] - d[0], d[4] - d[0]).normalized();
#endif
		/* Distance	*/
		const btScalar d0 = Lerp(Lerp(d[0], d[1], ix.f),
								 Lerp(d[3], d[2], ix.f), iy.f);
		const btScalar d1 = Lerp(Lerp(d[4], d[5], ix.f),
								 Lerp(d[7], d[6], ix.f), iy.f);
		return (Lerp(d0, d1, iz.f) - margin);
	}
	//
	void BuildCell(Cell& c)
	{
		const btVector3 org = btVector3((btScalar)c.c[0],
										(btScalar)c.c[1],
										(btScalar)c.c[2]) *
							  CELLSIZE * voxelsz;
		for (int k = 0; k <= CELLSIZE; ++k)
		{
			const btScalar z = voxelsz * k + org.z();
			for (int j = 0; j <= CELLSIZE; ++j)
			{
				const btScalar y = voxelsz * j + org.y();
				for (int i = 0; i <= CELLSIZE; ++i)
				{
					const btScalar x = voxelsz * i + org.x();
					c.d[i][j][k] = DistanceToShape(btVector3(x, y, z),
												   c.pclient);
				}
			}
		}
	}
	//
	static inline btScalar DistanceToShape(const btVector3& x,
										   const btCollisionShape* shape)
	{
		btTransform unit;
		unit.setIdentity();
		if (shape->isConvex())
		{
			btGjkEpaSolver2::sResults res;
			const btConvexShape* csh = static_cast<const btConvexShape*>(shape);
			return (btGjkEpaSolver2::SignedDistance(x, 0, csh, unit, res));
		}
		return (0);
	}
	//
	static inline IntFrac Decompose(btScalar x)
	{
		/* That one need a lot of improvements...	*/
		/* Remove test, faster floor...				*/
		IntFrac r;
		x /= CELLSIZE;
		const int o = x < 0 ? (int)(-x + 1) : 0;
		x += o;
		r.b = (int)x;
		const btScalar k = (x - r.b) * CELLSIZE;
		r.i = (int)k;
		r.f = k - r.i;
		r.b -= o;
		return (r);
	}
	//
	static inline btScalar Lerp(btScalar a, btScalar b, btScalar t)
	{
		return (a + (b - a) * t);
	}

	//
	static inline unsigned int Hash(int x, int y, int z, const btCollisionShape* shape)
	{
		struct btS
		{
			int x, y, z, w;
			void* p;
		};

		btS myset;
		//memset may be needed in case of additional (uninitialized) padding!
		//memset(&myset, 0, sizeof(btS));

		myset.x = x;
		myset.y = y;
		myset.z = z;
		myset.w = 0;
		myset.p = (void*)shape;
		const char* ptr = (const char*)&myset;

		unsigned int result = HsiehHash(ptr, sizeof(btS));

		return result;
	}
};

#endif  //BT_SPARSE_SDF_H
#include "btReducedDeformableBody.h"
#include "../btSoftBodyInternals.h"
#include "btReducedDeformableBodyHelpers.h"
#include "LinearMath/btTransformUtil.h"
#include <iostream>
#include <fstream>

btReducedDeformableBody::btReducedDeformableBody(btSoftBodyWorldInfo* worldInfo, int node_count, const btVector3* x, const btScalar* m)
 : btSoftBody(worldInfo, node_count, x, m), m_rigidOnly(false)
{
  // reduced deformable
  m_reducedModel = true;
  m_nReduced = 0;
  m_nFull = 0;
  m_nodeIndexOffset = 0;

  m_transform_lock = false;
  m_ksScale = 1.0;
  m_rhoScale = 1.0;

  // rigid motion
  m_linearVelocity.setZero();
	m_angularVelocity.setZero();
  m_internalDeltaLinearVelocity.setZero();
  m_internalDeltaAngularVelocity.setZero();
  m_angularVelocityFromReduced.setZero();
  m_internalDeltaAngularVelocityFromReduced.setZero();
	m_angularFactor.setValue(1, 1, 1);
	m_linearFactor.setValue(1, 1, 1);
  // m_invInertiaLocal.setValue(1, 1, 1);
  m_invInertiaLocal.setIdentity();
  m_mass = 0.0;
  m_inverseMass = 0.0;

  m_linearDamping = 0;
  m_angularDamping = 0;

  // Rayleigh damping
  m_dampingAlpha = 0;
  m_dampingBeta = 0;

  m_rigidTransformWorld.setIdentity();
}

void btReducedDeformableBody::setReducedModes(int num_modes, int full_size)
{
  m_nReduced = num_modes;
  m_nFull = full_size;
  m_reducedDofs.resize(m_nReduced, 0);
  m_reducedDofsBuffer.resize(m_nReduced, 0);
  m_reducedVelocity.resize(m_nReduced, 0);
  m_reducedVelocityBuffer.resize(m_nReduced, 0);
  m_reducedForceElastic.resize(m_nReduced, 0);
  m_reducedForceDamping.resize(m_nReduced, 0);
  m_reducedForceExternal.resize(m_nReduced, 0);
  m_internalDeltaReducedVelocity.resize(m_nReduced, 0);
  m_nodalMass.resize(full_size, 0);
  m_localMomentArm.resize(m_nFull);
}

void btReducedDeformableBody::setMassProps(const tDenseArray& mass_array)
{
  btScalar total_mass = 0;
  btVector3 CoM(0, 0, 0);
	for (int i = 0; i < m_nFull; ++i)
	{
		m_nodalMass[i] = m_rhoScale * mass_array[i];
		m_nodes[i].m_im = mass_array[i] > 0 ? 1.0 / (m_rhoScale * mass_array[i]) : 0;
		total_mass += m_rhoScale * mass_array[i];

    CoM += m_nodalMass[i] * m_nodes[i].m_x;
	}
  // total rigid body mass
  m_mass = total_mass;
  m_inverseMass = total_mass > 0 ? 1.0 / total_mass : 0;
  // original CoM
  m_initialCoM = CoM / total_mass;
}

void btReducedDeformableBody::setInertiaProps()
{
  // make sure the initial CoM is at the origin (0,0,0)
  // for (int i = 0; i < m_nFull; ++i)
  // {
  //   m_nodes[i].m_x -= m_initialCoM;
  // }
  // m_initialCoM.setZero();
  m_rigidTransformWorld.setOrigin(m_initialCoM);
  m_interpolationWorldTransform = m_rigidTransformWorld;
  
  updateLocalInertiaTensorFromNodes();

  // update world inertia tensor
  btMatrix3x3 rotation;
  rotation.setIdentity();
  updateInitialInertiaTensor(rotation);
  updateInertiaTensor();
  m_interpolateInvInertiaTensorWorld = m_invInertiaTensorWorld;
}

void btReducedDeformableBody::setRigidVelocity(const btVector3& v)
{
  m_linearVelocity = v;
}

void btReducedDeformableBody::setRigidAngularVelocity(const btVector3& omega)
{
  m_angularVelocity = omega;
}

void btReducedDeformableBody::setStiffnessScale(const btScalar ks)
{
  m_ksScale = ks;
}

void btReducedDeformableBody::setMassScale(const btScalar rho)
{
  m_rhoScale = rho;
}

void btReducedDeformableBody::setFixedNodes(const int n_node)
{
  m_fixedNodes.push_back(n_node);
  m_nodes[n_node].m_im = 0;   // set inverse mass to be zero for the constraint solver.
}

void btReducedDeformableBody::setDamping(const btScalar alpha, const btScalar beta)
{
  m_dampingAlpha = alpha;
  m_dampingBeta = beta;
}

void btReducedDeformableBody::internalInitialization()
{
  // zeroing
  endOfTimeStepZeroing();
  // initialize rest position
  updateRestNodalPositions();
  // initialize local nodal moment arm form the CoM
  updateLocalMomentArm();
  // initialize projection matrix
  updateExternalForceProjectMatrix(false);
}

void btReducedDeformableBody::updateLocalMomentArm()
{
  TVStack delta_x;
  delta_x.resize(m_nFull);

  for (int i = 0; i < m_nFull; ++i)
  {
    for (int k = 0; k < 3; ++k)
    {
      // compute displacement
      delta_x[i][k] = 0;
      for (int j = 0; j < m_nReduced; ++j) 
      {
        delta_x[i][k] += m_modes[j][3 * i + k] * m_reducedDofs[j];
      }
    }
    // get new moment arm Sq + x0
    m_localMomentArm[i] = m_x0[i] - m_initialCoM + delta_x[i];
  }
}

void btReducedDeformableBody::updateExternalForceProjectMatrix(bool initialized)
{
  // if not initialized, need to compute both P_A and Cq
  // otherwise, only need to udpate Cq
  if (!initialized)
  {
    // resize
    m_projPA.resize(m_nReduced);
    m_projCq.resize(m_nReduced);

    m_STP.resize(m_nReduced);
    m_MrInvSTP.resize(m_nReduced);

    // P_A
    for (int r = 0; r < m_nReduced; ++r)
    {
      m_projPA[r].resize(3 * m_nFull, 0);
      for (int i = 0; i < m_nFull; ++i)
      {
        btMatrix3x3 mass_scaled_i = Diagonal(1) - Diagonal(m_nodalMass[i] / m_mass);
        btVector3 s_ri(m_modes[r][3 * i], m_modes[r][3 * i + 1], m_modes[r][3 * i + 2]);
        btVector3 prod_i = mass_scaled_i * s_ri;

        for (int k = 0; k < 3; ++k)
          m_projPA[r][3 * i + k] = prod_i[k];

        // btScalar ratio = m_nodalMass[i] / m_mass;
        // m_projPA[r] += btVector3(- m_modes[r][3 * i] * ratio,
        //                          - m_modes[r][3 * i + 1] * ratio,
        //                          - m_modes[r][3 * i + 2] * ratio);
      }
    }
  }

  // C(q) is updated once per position update
  for (int r = 0; r < m_nReduced; ++r)
  {
  	m_projCq[r].resize(3 * m_nFull, 0);
    for (int i = 0; i < m_nFull; ++i)
    {
      btMatrix3x3 r_star = Cross(m_localMomentArm[i]);
      btVector3 s_ri(m_modes[r][3 * i], m_modes[r][3 * i + 1], m_modes[r][3 * i + 2]);
      btVector3 prod_i = r_star * m_invInertiaTensorWorld * r_star * s_ri;

      for (int k = 0; k < 3; ++k)
        m_projCq[r][3 * i + k] = m_nodalMass[i] * prod_i[k];

      // btVector3 si(m_modes[r][3 * i], m_modes[r][3 * i + 1], m_modes[r][3 * i + 2]);
      // m_projCq[r] += m_nodalMass[i] * si.cross(m_localMomentArm[i]);
    }
  }
}

void btReducedDeformableBody::endOfTimeStepZeroing()
{
  for (int i = 0; i < m_nReduced; ++i)
  {
    m_reducedForceElastic[i] = 0;
    m_reducedForceDamping[i] = 0;
    m_reducedForceExternal[i] = 0;
    m_internalDeltaReducedVelocity[i] = 0;
    m_reducedDofsBuffer[i] = m_reducedDofs[i];
    m_reducedVelocityBuffer[i] = m_reducedVelocity[i];
  }
  // std::cout << "zeroed!\n";
}

void btReducedDeformableBody::applyInternalVelocityChanges()
{
  m_linearVelocity += m_internalDeltaLinearVelocity;
  m_angularVelocity += m_internalDeltaAngularVelocity;
  m_internalDeltaLinearVelocity.setZero();
  m_internalDeltaAngularVelocity.setZero();
  for (int r = 0; r < m_nReduced; ++r)
  {
    m_reducedVelocity[r] += m_internalDeltaReducedVelocity[r];
    m_internalDeltaReducedVelocity[r] = 0;
  }
}

void btReducedDeformableBody::predictIntegratedTransform(btScalar dt, btTransform& predictedTransform)
{
	btTransformUtil::integrateTransform(m_rigidTransformWorld, m_linearVelocity, m_angularVelocity, dt, predictedTransform);
}

void btReducedDeformableBody::updateReducedDofs(btScalar solverdt)
{
  for (int r = 0; r < m_nReduced; ++r)
  { 
    m_reducedDofs[r] = m_reducedDofsBuffer[r] + solverdt * m_reducedVelocity[r];
  }
}

void btReducedDeformableBody::mapToFullPosition(const btTransform& ref_trans)
{
  btVector3 origin = ref_trans.getOrigin();
  btMatrix3x3 rotation = ref_trans.getBasis();
  

  for (int i = 0; i < m_nFull; ++i)
  {
    m_nodes[i].m_x = rotation * m_localMomentArm[i] + origin;
    m_nodes[i].m_q = m_nodes[i].m_x;
  }
}

void btReducedDeformableBody::updateReducedVelocity(btScalar solverdt)
{
  // update reduced velocity
  for (int r = 0; r < m_nReduced; ++r)
  {
    // the reduced mass is always identity!
    btScalar delta_v = 0;
    delta_v = solverdt * (m_reducedForceElastic[r] + m_reducedForceDamping[r]);
    // delta_v = solverdt * (m_reducedForceElastic[r] + m_reducedForceDamping[r] + m_reducedForceExternal[r]);
    m_reducedVelocity[r] = m_reducedVelocityBuffer[r] + delta_v;
  }
}

void btReducedDeformableBody::mapToFullVelocity(const btTransform& ref_trans)
{
  // compute the reduced contribution to the angular velocity
  // btVector3 sum_linear(0, 0, 0);
  // btVector3 sum_angular(0, 0, 0);
  // m_linearVelocityFromReduced.setZero();
  // m_angularVelocityFromReduced.setZero();
  // for (int i = 0; i < m_nFull; ++i)
  // {
  //   btVector3 r_com = ref_trans.getBasis() * m_localMomentArm[i];
  //   btMatrix3x3 r_star = Cross(r_com);

  //   btVector3 v_from_reduced(0, 0, 0);
  //   for (int k = 0; k < 3; ++k)
  //   {
  //     for (int r = 0; r < m_nReduced; ++r)
  //     {
  //       v_from_reduced[k] += m_modes[r][3 * i + k] * m_reducedVelocity[r];
  //     }
  //   }

  //   btVector3 delta_linear = m_nodalMass[i] * v_from_reduced;
  //   btVector3 delta_angular = m_nodalMass[i] * (r_star * ref_trans.getBasis() * v_from_reduced);
  //   sum_linear += delta_linear;
  //   sum_angular += delta_angular;
  //   // std::cout << "delta_linear: " << delta_linear[0] << "\t" << delta_linear[1] << "\t" << delta_linear[2] << "\n";
  //   // std::cout << "delta_angular: " << delta_angular[0] << "\t" << delta_angular[1] << "\t" << delta_angular[2] << "\n";
  //   // std::cout << "sum_linear: " << sum_linear[0] << "\t" << sum_linear[1] << "\t" << sum_linear[2] << "\n";
  //   // std::cout << "sum_angular: " << sum_angular[0] << "\t" << sum_angular[1] << "\t" << sum_angular[2] << "\n";
  // }
  // m_linearVelocityFromReduced = 1.0 / m_mass * (ref_trans.getBasis() * sum_linear);
  // m_angularVelocityFromReduced = m_interpolateInvInertiaTensorWorld * sum_angular;

  // m_linearVelocity -= m_linearVelocityFromReduced;
  // m_angularVelocity -= m_angularVelocityFromReduced;

  for (int i = 0; i < m_nFull; ++i)
  {
    m_nodes[i].m_v = computeNodeFullVelocity(ref_trans, i);
  }
}

const btVector3 btReducedDeformableBody::computeTotalAngularMomentum() const
{
  btVector3 L_rigid = m_invInertiaTensorWorld.inverse() * m_angularVelocity;
  btVector3 L_reduced(0, 0, 0);
  btMatrix3x3 omega_prime_star = Cross(m_angularVelocityFromReduced);

  for (int i = 0; i < m_nFull; ++i)
  {
    btVector3 r_com = m_rigidTransformWorld.getBasis() * m_localMomentArm[i];
    btMatrix3x3 r_star = Cross(r_com);

    btVector3 v_from_reduced(0, 0, 0);
    for (int k = 0; k < 3; ++k)
    {
      for (int r = 0; r < m_nReduced; ++r)
      {
        v_from_reduced[k] += m_modes[r][3 * i + k] * m_reducedVelocity[r];
      }
    }

    L_reduced += m_nodalMass[i] * (r_star * (m_rigidTransformWorld.getBasis() * v_from_reduced - omega_prime_star * r_com));
    // L_reduced += m_nodalMass[i] * (r_star * (m_rigidTransformWorld.getBasis() * v_from_reduced));
  }
  return L_rigid + L_reduced;
}

const btVector3 btReducedDeformableBody::computeNodeFullVelocity(const btTransform& ref_trans, int n_node) const
{
  btVector3 v_from_reduced(0, 0, 0);
  btVector3 r_com = ref_trans.getBasis() * m_localMomentArm[n_node];
  // compute velocity contributed by the reduced velocity
  for (int k = 0; k < 3; ++k)
  {
    for (int r = 0; r < m_nReduced; ++r)
    {
      v_from_reduced[k] += m_modes[r][3 * n_node + k] * m_reducedVelocity[r];
    }
  }
  // get new velocity
  btVector3 vel = m_angularVelocity.cross(r_com) + 
                  ref_trans.getBasis() * v_from_reduced +
                  m_linearVelocity;
  return vel;
}

const btVector3 btReducedDeformableBody::internalComputeNodeDeltaVelocity(const btTransform& ref_trans, int n_node) const
{
  btVector3 deltaV_from_reduced(0, 0, 0);
  btVector3 r_com = ref_trans.getBasis() * m_localMomentArm[n_node];

  // compute velocity contributed by the reduced velocity
  for (int k = 0; k < 3; ++k)
  {
    for (int r = 0; r < m_nReduced; ++r)
    {
      deltaV_from_reduced[k] += m_modes[r][3 * n_node + k] * m_internalDeltaReducedVelocity[r];
    }
  }

  // get delta velocity
  btVector3 deltaV = m_internalDeltaAngularVelocity.cross(r_com) + 
                     ref_trans.getBasis() * deltaV_from_reduced +
                     m_internalDeltaLinearVelocity;
  return deltaV;
}

void btReducedDeformableBody::proceedToTransform(btScalar dt, bool end_of_time_step)
{
  btTransformUtil::integrateTransform(m_rigidTransformWorld, m_linearVelocity, m_angularVelocity, dt, m_interpolationWorldTransform);
  updateInertiaTensor();
  // m_interpolateInvInertiaTensorWorld = m_interpolationWorldTransform.getBasis().scaled(m_invInertiaLocal) * m_interpolationWorldTransform.getBasis().transpose();
  m_rigidTransformWorld = m_interpolationWorldTransform;
  m_invInertiaTensorWorld = m_interpolateInvInertiaTensorWorld;
}

void btReducedDeformableBody::transformTo(const btTransform& trs)
{
	btTransform current_transform = getRigidTransform();
	btTransform new_transform(trs.getBasis() * current_transform.getBasis().transpose(),
                            trs.getOrigin() - current_transform.getOrigin());
  transform(new_transform);
}

void btReducedDeformableBody::transform(const btTransform& trs)
{
  m_transform_lock = true;

  // transform mesh
  {
    const btScalar margin = getCollisionShape()->getMargin();
    ATTRIBUTE_ALIGNED16(btDbvtVolume)
    vol;

    btVector3 CoM = m_rigidTransformWorld.getOrigin();
    btVector3 translation = trs.getOrigin();
    btMatrix3x3 rotation = trs.getBasis();

    for (int i = 0; i < m_nodes.size(); ++i)
    {
      Node& n = m_nodes[i];
      n.m_x = rotation * (n.m_x - CoM) + CoM + translation;
      n.m_q = rotation * (n.m_q - CoM) + CoM + translation;
      n.m_n = rotation * n.m_n;
      vol = btDbvtVolume::FromCR(n.m_x, margin);

      m_ndbvt.update(n.m_leaf, vol);
    }
    updateNormals();
    updateBounds();
    updateConstants();
  }

  // update modes
  updateModesByRotation(trs.getBasis());

  // update inertia tensor
  updateInitialInertiaTensor(trs.getBasis());
  updateInertiaTensor();
  m_interpolateInvInertiaTensorWorld = m_invInertiaTensorWorld;
  
  // update rigid frame (No need to update the rotation. Nodes have already been updated.)
  m_rigidTransformWorld.setOrigin(m_initialCoM + trs.getOrigin());
  m_interpolationWorldTransform = m_rigidTransformWorld;
  m_initialCoM = m_rigidTransformWorld.getOrigin();

  internalInitialization();
}

void btReducedDeformableBody::scale(const btVector3& scl)
{
  // Scaling the mesh after transform is applied is not allowed
  btAssert(!m_transform_lock);

  // scale the mesh
  {
    const btScalar margin = getCollisionShape()->getMargin();
    ATTRIBUTE_ALIGNED16(btDbvtVolume)
    vol;

    btVector3 CoM = m_rigidTransformWorld.getOrigin();

    for (int i = 0; i < m_nodes.size(); ++i)
    {
      Node& n = m_nodes[i];
      n.m_x = (n.m_x - CoM) * scl + CoM;
      n.m_q = (n.m_q - CoM) * scl + CoM;
      vol = btDbvtVolume::FromCR(n.m_x, margin);
      m_ndbvt.update(n.m_leaf, vol);
    }
    updateNormals();
    updateBounds();
    updateConstants();
    initializeDmInverse();
  }

  // update inertia tensor
  updateLocalInertiaTensorFromNodes();

  btMatrix3x3 id;
  id.setIdentity();
  updateInitialInertiaTensor(id);   // there is no rotation, but the local inertia tensor has changed
  updateInertiaTensor();
  m_interpolateInvInertiaTensorWorld = m_invInertiaTensorWorld;

  internalInitialization();
}

void btReducedDeformableBody::setTotalMass(btScalar mass, bool fromfaces)
{
  // Changing the total mass after transform is applied is not allowed
  btAssert(!m_transform_lock);

  btScalar scale_ratio = mass / m_mass;

  // update nodal mass
  for (int i = 0; i < m_nFull; ++i)
  {
    m_nodalMass[i] *= scale_ratio;
  }
  m_mass = mass;
  m_inverseMass = mass > 0 ? 1.0 / mass : 0;

  // update inertia tensors
  updateLocalInertiaTensorFromNodes();

  btMatrix3x3 id;
  id.setIdentity();
  updateInitialInertiaTensor(id);   // there is no rotation, but the local inertia tensor has changed
  updateInertiaTensor();
  m_interpolateInvInertiaTensorWorld = m_invInertiaTensorWorld;

  internalInitialization();
}

void btReducedDeformableBody::updateRestNodalPositions()
{
  // update reset nodal position
  m_x0.resize(m_nFull);
  for (int i = 0; i < m_nFull; ++i)
  {
    m_x0[i] = m_nodes[i].m_x;
  }
}

// reference notes:
// https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-07-dynamics-fall-2009/lecture-notes/MIT16_07F09_Lec26.pdf
void btReducedDeformableBody::updateLocalInertiaTensorFromNodes()
{
  btMatrix3x3 inertia_tensor;
  inertia_tensor.setZero();

  for (int p = 0; p < m_nFull; ++p)
  {
    btMatrix3x3 particle_inertia;
    particle_inertia.setZero();

    btVector3 r = m_nodes[p].m_x - m_initialCoM;

    particle_inertia[0][0] = m_nodalMass[p] * (r[1] * r[1] + r[2] * r[2]);
    particle_inertia[1][1] = m_nodalMass[p] * (r[0] * r[0] + r[2] * r[2]);
    particle_inertia[2][2] = m_nodalMass[p] * (r[0] * r[0] + r[1] * r[1]);

    particle_inertia[0][1] = - m_nodalMass[p] * (r[0] * r[1]);
    particle_inertia[0][2] = - m_nodalMass[p] * (r[0] * r[2]);
    particle_inertia[1][2] = - m_nodalMass[p] * (r[1] * r[2]);

    particle_inertia[1][0] = particle_inertia[0][1];
    particle_inertia[2][0] = particle_inertia[0][2];
    particle_inertia[2][1] = particle_inertia[1][2];

    inertia_tensor += particle_inertia;
  }
  m_invInertiaLocal = inertia_tensor.inverse();
}

void btReducedDeformableBody::updateInitialInertiaTensor(const btMatrix3x3& rotation)
{
  // m_invInertiaTensorWorldInitial = rotation.scaled(m_invInertiaLocal) * rotation.transpose();
  m_invInertiaTensorWorldInitial = rotation * m_invInertiaLocal * rotation.transpose();
}

void btReducedDeformableBody::updateModesByRotation(const btMatrix3x3& rotation)
{
  for (int r = 0; r < m_nReduced; ++r)
  {
    for (int i = 0; i < m_nFull; ++i)
    {
      btVector3 nodal_disp(m_modes[r][3 * i], m_modes[r][3 * i + 1], m_modes[r][3 * i + 2]);
      nodal_disp = rotation * nodal_disp;

      for (int k = 0; k < 3; ++k)
      {
        m_modes[r][3 * i + k] = nodal_disp[k];
      }
    }
  }
}

void btReducedDeformableBody::updateInertiaTensor()
{
	m_invInertiaTensorWorld = m_rigidTransformWorld.getBasis() * m_invInertiaTensorWorldInitial * m_rigidTransformWorld.getBasis().transpose();
}

void btReducedDeformableBody::applyDamping(btScalar timeStep)
{
  m_linearVelocity *= btScalar(1) - m_linearDamping;
  m_angularDamping *= btScalar(1) - m_angularDamping;
}

void btReducedDeformableBody::applyCentralImpulse(const btVector3& impulse)
{
  m_linearVelocity += impulse * m_linearFactor * m_inverseMass;
  #if defined(BT_CLAMP_VELOCITY_TO) && BT_CLAMP_VELOCITY_TO > 0
  clampVelocity(m_linearVelocity);
  #endif
}

void btReducedDeformableBody::applyTorqueImpulse(const btVector3& torque)
{
  m_angularVelocity += m_interpolateInvInertiaTensorWorld * torque * m_angularFactor;
  #if defined(BT_CLAMP_VELOCITY_TO) && BT_CLAMP_VELOCITY_TO > 0
  clampVelocity(m_angularVelocity);
  #endif
}

void btReducedDeformableBody::internalApplyRigidImpulse(const btVector3& impulse, const btVector3& rel_pos)
{
  if (m_inverseMass == btScalar(0.))
  {
    std::cout << "something went wrong...probably didn't initialize?\n";
    btAssert(false);
  }
  // delta linear velocity
  m_internalDeltaLinearVelocity += impulse * m_linearFactor * m_inverseMass;
  // delta angular velocity
  btVector3 torque = rel_pos.cross(impulse * m_linearFactor);
  m_internalDeltaAngularVelocity += m_interpolateInvInertiaTensorWorld * torque * m_angularFactor;
}

btVector3 btReducedDeformableBody::getRelativePos(int n_node)
{
  btMatrix3x3 rotation = m_interpolationWorldTransform.getBasis();
  btVector3 ri = rotation * m_localMomentArm[n_node];
  return ri;
}

btMatrix3x3 btReducedDeformableBody::getImpulseFactor(int n_node)
{
  // relative position
  btMatrix3x3 rotation = m_interpolationWorldTransform.getBasis();
  btVector3 ri = rotation * m_localMomentArm[n_node];
  btMatrix3x3 ri_skew = Cross(ri);

  // calculate impulse factor
  // rigid part
  btScalar inv_mass = m_nodalMass[n_node] > btScalar(0) ? btScalar(1) / m_mass : btScalar(0);
  btMatrix3x3 K1 = Diagonal(inv_mass);
  K1 -= ri_skew * m_interpolateInvInertiaTensorWorld * ri_skew;

  // reduced deformable part
  btMatrix3x3 SA;
  SA.setZero();
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      for (int r = 0; r < m_nReduced; ++r)
      {
        SA[i][j] += m_modes[r][3 * n_node + i] * (m_projPA[r][3 * n_node + j] + m_projCq[r][3 * n_node + j]);
      }
    }
  }
  btMatrix3x3 RSARinv = rotation * SA * rotation.transpose();


  TVStack omega_helper; // Sum_i m_i r*_i R S_i
  omega_helper.resize(m_nReduced);
  for (int r = 0; r < m_nReduced; ++r)
  {
    omega_helper[r].setZero();
    for (int i = 0; i < m_nFull; ++i)
    {
      btMatrix3x3 mi_rstar_i = rotation * Cross(m_localMomentArm[i]) * m_nodalMass[i];
      btVector3 s_ri(m_modes[r][3 * i], m_modes[r][3 * i + 1], m_modes[r][3 * i + 2]);
      omega_helper[r] += mi_rstar_i * rotation * s_ri;
    }
  }

  btMatrix3x3 sum_multiply_A;
  sum_multiply_A.setZero();
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      for (int r = 0; r < m_nReduced; ++r)
      {
        sum_multiply_A[i][j] += omega_helper[r][i] * (m_projPA[r][3 * n_node + j] + m_projCq[r][3 * n_node + j]);
      }
    }
  }

  btMatrix3x3 K2 = RSARinv + ri_skew * m_interpolateInvInertiaTensorWorld * sum_multiply_A * rotation.transpose();

  return m_rigidOnly ? K1 : K1 + K2;
}

void btReducedDeformableBody::internalApplyFullSpaceImpulse(const btVector3& impulse, const btVector3& rel_pos, int n_node, btScalar dt)
{
  if (!m_rigidOnly)
  {
    // apply impulse force
    applyFullSpaceNodalForce(impulse / dt, n_node);

    // update delta damping force
    tDenseArray reduced_vel_tmp;
    reduced_vel_tmp.resize(m_nReduced);
    for (int r = 0; r < m_nReduced; ++r)
    {
      reduced_vel_tmp[r] = m_reducedVelocity[r] + m_internalDeltaReducedVelocity[r];
    }
    applyReducedDampingForce(reduced_vel_tmp);
    // applyReducedDampingForce(m_internalDeltaReducedVelocity);

    // delta reduced velocity
    for (int r = 0; r < m_nReduced; ++r)
    {
      // The reduced mass is always identity!
      m_internalDeltaReducedVelocity[r] += dt * (m_reducedForceDamping[r] + m_reducedForceExternal[r]);
    }
  }

  internalApplyRigidImpulse(impulse, rel_pos);
}

void btReducedDeformableBody::applyFullSpaceNodalForce(const btVector3& f_ext, int n_node)
{
  // f_local = R^-1 * f_ext //TODO: interpoalted transfrom
  // btVector3 f_local = m_rigidTransformWorld.getBasis().transpose() * f_ext;
  btVector3 f_local = m_interpolationWorldTransform.getBasis().transpose() * f_ext;

  // f_ext_r = [S^T * P]_{n_node} * f_local
  tDenseArray f_ext_r;
  f_ext_r.resize(m_nReduced, 0);
  for (int r = 0; r < m_nReduced; ++r)
  {
    m_reducedForceExternal[r] = 0;
    for (int k = 0; k < 3; ++k)
    {
      f_ext_r[r] += (m_projPA[r][3 * n_node + k] + m_projCq[r][3 * n_node + k]) * f_local[k];
    }

    m_reducedForceExternal[r] += f_ext_r[r];
  }
}

void btReducedDeformableBody::applyRigidGravity(const btVector3& gravity, btScalar dt)
{
  // update rigid frame velocity
  m_linearVelocity += dt * gravity;
}

void btReducedDeformableBody::applyReducedElasticForce(const tDenseArray& reduce_dofs)
{
  for (int r = 0; r < m_nReduced; ++r) 
  {
    m_reducedForceElastic[r] = - m_ksScale * m_Kr[r] * reduce_dofs[r];
  }
}

void btReducedDeformableBody::applyReducedDampingForce(const tDenseArray& reduce_vel)
{
  for (int r = 0; r < m_nReduced; ++r) 
  {
    m_reducedForceDamping[r] = - m_dampingBeta * m_ksScale * m_Kr[r] * reduce_vel[r];
  }
}

btScalar btReducedDeformableBody::getTotalMass() const
{
  return m_mass;
}

btTransform& btReducedDeformableBody::getRigidTransform()
{
  return m_rigidTransformWorld;
}

const btVector3& btReducedDeformableBody::getLinearVelocity() const
{
  return m_linearVelocity;
}

const btVector3& btReducedDeformableBody::getAngularVelocity() const
{
  return m_angularVelocity;
}

void btReducedDeformableBody::disableReducedModes(const bool rigid_only)
{
  m_rigidOnly = rigid_only;
}

bool btReducedDeformableBody::isReducedModesOFF() const
{
  return m_rigidOnly;
}
#ifndef BT_REDUCED_SOFT_BODY_H
#define BT_REDUCED_SOFT_BODY_H

#include "../btSoftBody.h"
#include "LinearMath/btAlignedObjectArray.h"
#include "LinearMath/btVector3.h"
#include "LinearMath/btMatrix3x3.h"
#include "LinearMath/btTransform.h"

// Reduced deformable body is a simplified deformable object embedded in a rigid frame.
class btReducedDeformableBody : public btSoftBody
{
 public:
  //
  //  Typedefs
  //
  typedef btAlignedObjectArray<btVector3> TVStack;
  // typedef btAlignedObjectArray<btMatrix3x3> tBlockDiagMatrix;
  typedef btAlignedObjectArray<btScalar> tDenseArray;
  typedef btAlignedObjectArray<btAlignedObjectArray<btScalar> > tDenseMatrix;

 private:
  // flag to turn off the reduced modes
  bool m_rigidOnly;

  // Flags for transform. Once transform is applied, users cannot scale the mesh or change its total mass.
  bool m_transform_lock;

  // scaling factors
  btScalar m_rhoScale;         // mass density scale
  btScalar m_ksScale;          // stiffness scale

  // projection matrix
  tDenseMatrix m_projPA;        // Eqn. 4.11 from Rahul Sheth's thesis
  tDenseMatrix m_projCq;
  tDenseArray m_STP;
  tDenseArray m_MrInvSTP;

  TVStack m_localMomentArm; // Sq + x0

  btVector3 m_internalDeltaLinearVelocity;
  btVector3 m_internalDeltaAngularVelocity;
  tDenseArray m_internalDeltaReducedVelocity;
  
  btVector3 m_linearVelocityFromReduced;  // contribution to the linear velocity from reduced velocity
  btVector3 m_angularVelocityFromReduced; // contribution to the angular velocity from reduced velocity
  btVector3 m_internalDeltaAngularVelocityFromReduced;

 protected:
  // rigid frame
  btScalar m_mass;          // total mass of the rigid frame
  btScalar m_inverseMass;   // inverse of the total mass of the rigid frame
  btVector3 m_linearVelocity;
  btVector3 m_angularVelocity;
  btScalar m_linearDamping;    // linear damping coefficient
  btScalar m_angularDamping;    // angular damping coefficient
  btVector3 m_linearFactor;
  btVector3 m_angularFactor;
  // btVector3 m_invInertiaLocal;
  btMatrix3x3 m_invInertiaLocal;
  btTransform m_rigidTransformWorld;
  btMatrix3x3 m_invInertiaTensorWorldInitial;
  btMatrix3x3 m_invInertiaTensorWorld;
  btMatrix3x3 m_interpolateInvInertiaTensorWorld;
  btVector3 m_initialCoM;  // initial center of mass (original of the m_rigidTransformWorld)

  // damping
  btScalar m_dampingAlpha;
  btScalar m_dampingBeta;

 public:
  //
  //  Fields
  //

  // reduced space
  int m_nReduced;
  int m_nFull;
  tDenseMatrix m_modes;														// modes of the reduced deformable model. Each inner array is a mode, outer array size = n_modes
  tDenseArray m_reducedDofs;				   // Reduced degree of freedom
  tDenseArray m_reducedDofsBuffer;     // Reduced degree of freedom at t^n
  tDenseArray m_reducedVelocity;		   // Reduced velocity array
  tDenseArray m_reducedVelocityBuffer; // Reduced velocity array at t^n
  tDenseArray m_reducedForceExternal;          // reduced external force
  tDenseArray m_reducedForceElastic;           // reduced internal elastic force
  tDenseArray m_reducedForceDamping;           // reduced internal damping force
  tDenseArray m_eigenvalues;		// eigenvalues of the reduce deformable model
  tDenseArray m_Kr;	// reduced stiffness matrix
  
  // full space
  TVStack m_x0;					     				 // Rest position
  tDenseArray m_nodalMass;           // Mass on each node
  btAlignedObjectArray<int> m_fixedNodes; // index of the fixed nodes
  int m_nodeIndexOffset;             // offset of the node index needed for contact solver when there are multiple reduced deformable body in the world.

  // contacts
  btAlignedObjectArray<int> m_contactNodesList;

  //
  // Api
  //
  btReducedDeformableBody(btSoftBodyWorldInfo* worldInfo, int node_count, const btVector3* x, const btScalar* m);

  ~btReducedDeformableBody() {}

  //
  // initializing helpers
  //
  void internalInitialization();

  void setReducedModes(int num_modes, int full_size);

  void setMassProps(const tDenseArray& mass_array);

  void setInertiaProps();

  void setRigidVelocity(const btVector3& v);

  void setRigidAngularVelocity(const btVector3& omega);

  void setStiffnessScale(const btScalar ks);

  void setMassScale(const btScalar rho);

  void setFixedNodes(const int n_node);

  void setDamping(const btScalar alpha, const btScalar beta);

  void disableReducedModes(const bool rigid_only);

  virtual void setTotalMass(btScalar mass, bool fromfaces = false);

  //
  // various internal updates
  //
  virtual void transformTo(const btTransform& trs);
  virtual void transform(const btTransform& trs);
  // caution: 
  // need to use scale before using transform, because the scale is performed in the local frame 
  // (i.e., may have some rotation already, but the m_rigidTransformWorld doesn't have this info)
  virtual void scale(const btVector3& scl);

 private:
  void updateRestNodalPositions();

  void updateInitialInertiaTensor(const btMatrix3x3& rotation);

  void updateLocalInertiaTensorFromNodes();

  void updateInertiaTensor();

  void updateModesByRotation(const btMatrix3x3& rotation);
 
 public:
  void updateLocalMomentArm();

  void predictIntegratedTransform(btScalar dt, btTransform& predictedTransform);

  // update the external force projection matrix 
  void updateExternalForceProjectMatrix(bool initialized);

  void endOfTimeStepZeroing();

  void applyInternalVelocityChanges();

  //
  // position and velocity update related
  //

  // compute reduced degree of freedoms
  void updateReducedDofs(btScalar solverdt);

  // compute reduced velocity update (for explicit time stepping)
  void updateReducedVelocity(btScalar solverdt);

  // map to full degree of freedoms
  void mapToFullPosition(const btTransform& ref_trans);

  // compute full space velocity from the reduced velocity
  void mapToFullVelocity(const btTransform& ref_trans);

  // compute total angular momentum
  const btVector3 computeTotalAngularMomentum() const;

  // get a single node's full space velocity from the reduced velocity
  const btVector3 computeNodeFullVelocity(const btTransform& ref_trans, int n_node) const;

  // get a single node's all delta velocity
  const btVector3 internalComputeNodeDeltaVelocity(const btTransform& ref_trans, int n_node) const;

  //
  // rigid motion related
  //
  void applyDamping(btScalar timeStep);

  void applyCentralImpulse(const btVector3& impulse);

  void applyTorqueImpulse(const btVector3& torque);

  void proceedToTransform(btScalar dt, bool end_of_time_step);

  //
  // force related
  //

  // apply impulse to the rigid frame
  void internalApplyRigidImpulse(const btVector3& impulse, const btVector3& rel_pos);

  // apply impulse to nodes in the full space
  void internalApplyFullSpaceImpulse(const btVector3& impulse, const btVector3& rel_pos, int n_node, btScalar dt);

  // apply nodal external force in the full space
  void applyFullSpaceNodalForce(const btVector3& f_ext, int n_node);

  // apply gravity to the rigid frame
  void applyRigidGravity(const btVector3& gravity, btScalar dt);

  // apply reduced elastic force
  void applyReducedElasticForce(const tDenseArray& reduce_dofs);

  // apply reduced damping force
  void applyReducedDampingForce(const tDenseArray& reduce_vel);

  // calculate the impulse factor
  virtual btMatrix3x3 getImpulseFactor(int n_node);

  // get relative position from a node to the CoM of the rigid frame
  btVector3 getRelativePos(int n_node);

  //
  // accessors
  //
  bool isReducedModesOFF() const;
  btScalar getTotalMass() const;
  btTransform& getRigidTransform();
  const btVector3& getLinearVelocity() const;
	const btVector3& getAngularVelocity() const;

  #if defined(BT_CLAMP_VELOCITY_TO) && BT_CLAMP_VELOCITY_TO > 0
  void clampVelocity(btVector3& v) const {
      v.setX(
          fmax(-BT_CLAMP_VELOCITY_TO,
                fmin(BT_CLAMP_VELOCITY_TO, v.getX()))
      );
      v.setY(
          fmax(-BT_CLAMP_VELOCITY_TO,
                fmin(BT_CLAMP_VELOCITY_TO, v.getY()))
      );
      v.setZ(
          fmax(-BT_CLAMP_VELOCITY_TO,
                fmin(BT_CLAMP_VELOCITY_TO, v.getZ()))
      );
  }
  #endif
};

#endif // BT_REDUCED_SOFT_BODY_H
#include "btReducedDeformableBodyHelpers.h"
#include "../btSoftBodyHelpers.h"
#include <iostream>
#include <string>
#include <sstream>

btReducedDeformableBody* btReducedDeformableBodyHelpers::createReducedDeformableObject(btSoftBodyWorldInfo& worldInfo, const std::string& file_path, const std::string& vtk_file, const int num_modes, bool rigid_only) {
	std::string filename = file_path + vtk_file;
	btReducedDeformableBody* rsb = btReducedDeformableBodyHelpers::createFromVtkFile(worldInfo, filename.c_str());
	
	rsb->setReducedModes(num_modes, rsb->m_nodes.size());
	btReducedDeformableBodyHelpers::readReducedDeformableInfoFromFiles(rsb, file_path.c_str());
	
	rsb->disableReducedModes(rigid_only);

	return rsb;
}

btReducedDeformableBody* btReducedDeformableBodyHelpers::createFromVtkFile(btSoftBodyWorldInfo& worldInfo, const char* vtk_file)
{
	std::ifstream fs;
	fs.open(vtk_file);
	btAssert(fs);

	typedef btAlignedObjectArray<int> Index;
	std::string line;
	btAlignedObjectArray<btVector3> X;
	btVector3 position;
	btAlignedObjectArray<Index> indices;
	bool reading_points = false;
	bool reading_tets = false;
	size_t n_points = 0;
	size_t n_tets = 0;
	size_t x_count = 0;
	size_t indices_count = 0;
	while (std::getline(fs, line))
	{
		std::stringstream ss(line);
		if (line.size() == (size_t)(0))
		{
		}
		else if (line.substr(0, 6) == "POINTS")
		{
			reading_points = true;
			reading_tets = false;
			ss.ignore(128, ' ');  // ignore "POINTS"
			ss >> n_points;
			X.resize(n_points);
		}
		else if (line.substr(0, 5) == "CELLS")
		{
			reading_points = false;
			reading_tets = true;
			ss.ignore(128, ' ');  // ignore "CELLS"
			ss >> n_tets;
			indices.resize(n_tets);
		}
		else if (line.substr(0, 10) == "CELL_TYPES")
		{
			reading_points = false;
			reading_tets = false;
		}
		else if (reading_points)
		{
			btScalar p;
			ss >> p;
			position.setX(p);
			ss >> p;
			position.setY(p);
			ss >> p;
			position.setZ(p);
			//printf("v %f %f %f\n", position.getX(), position.getY(), position.getZ());
			X[x_count++] = position;
		}
		else if (reading_tets)
		{
			int d;
			ss >> d;
			if (d != 4)
			{
				printf("Load deformable failed: Only Tetrahedra are supported in VTK file.\n");
				fs.close();
				return 0;
			}
			ss.ignore(128, ' ');  // ignore "4"
			Index tet;
			tet.resize(4);
			for (size_t i = 0; i < 4; i++)
			{
				ss >> tet[i];
				//printf("%d ", tet[i]);
			}
			//printf("\n");
			indices[indices_count++] = tet;
		}
	}
	btReducedDeformableBody* rsb = new btReducedDeformableBody(&worldInfo, n_points, &X[0], 0);

	for (int i = 0; i < n_tets; ++i)
	{
		const Index& ni = indices[i];
		rsb->appendTetra(ni[0], ni[1], ni[2], ni[3]);
		{
			rsb->appendLink(ni[0], ni[1], 0, true);
			rsb->appendLink(ni[1], ni[2], 0, true);
			rsb->appendLink(ni[2], ni[0], 0, true);
			rsb->appendLink(ni[0], ni[3], 0, true);
			rsb->appendLink(ni[1], ni[3], 0, true);
			rsb->appendLink(ni[2], ni[3], 0, true);
		}
	}

	btSoftBodyHelpers::generateBoundaryFaces(rsb);
	rsb->initializeDmInverse();
	rsb->m_tetraScratches.resize(rsb->m_tetras.size());
	rsb->m_tetraScratchesTn.resize(rsb->m_tetras.size());
	printf("Nodes:  %u\r\n", rsb->m_nodes.size());
	printf("Links:  %u\r\n", rsb->m_links.size());
	printf("Faces:  %u\r\n", rsb->m_faces.size());
	printf("Tetras: %u\r\n", rsb->m_tetras.size());

	fs.close();

	return rsb;
}

void btReducedDeformableBodyHelpers::readReducedDeformableInfoFromFiles(btReducedDeformableBody* rsb, const char* file_path)
{
	// read in eigenmodes, stiffness and mass matrices
	std::string eigenvalues_file = std::string(file_path) + "eigenvalues.bin";
	btReducedDeformableBodyHelpers::readBinaryVec(rsb->m_eigenvalues, rsb->m_nReduced, eigenvalues_file.c_str());

	std::string Kr_file = std::string(file_path) + "K_r_diag_mat.bin";
	btReducedDeformableBodyHelpers::readBinaryVec(rsb->m_Kr,  rsb->m_nReduced, Kr_file.c_str());

	// std::string Mr_file = std::string(file_path) + "M_r_diag_mat.bin";
	// btReducedDeformableBodyHelpers::readBinaryVec(rsb->m_Mr, rsb->m_nReduced, Mr_file.c_str());

	std::string modes_file = std::string(file_path) + "modes.bin";
	btReducedDeformableBodyHelpers::readBinaryMat(rsb->m_modes, rsb->m_nReduced, 3 * rsb->m_nFull, modes_file.c_str());
	
	// read in full nodal mass
	std::string M_file = std::string(file_path) + "M_diag_mat.bin";
	btAlignedObjectArray<btScalar> mass_array;
	btReducedDeformableBodyHelpers::readBinaryVec(mass_array, rsb->m_nFull, M_file.c_str());
	rsb->setMassProps(mass_array);
	
	// calculate the inertia tensor in the local frame 
 	rsb->setInertiaProps();

	// other internal initialization
	rsb->internalInitialization();
}

// read in a vector from the binary file
void btReducedDeformableBodyHelpers::readBinaryVec(btReducedDeformableBody::tDenseArray& vec, 
																				  	 const unsigned int n_size, 				// #entries read
																						 const char* file)
{
	std::ifstream f_in(file, std::ios::in | std::ios::binary);
	// first get size
	unsigned int size=0;
	f_in.read((char*)&size, 4);//sizeof(unsigned int));
	btAssert(size >= n_size); 	// make sure the #requested mode is smaller than the #available modes

	// read data
	vec.resize(n_size);
	double temp;
	for (unsigned int i = 0; i < n_size; ++i)
	{
		f_in.read((char*)&temp, sizeof(double));
		vec[i] = btScalar(temp);
	}
  f_in.close();
}

// read in a matrix from the binary file
void btReducedDeformableBodyHelpers::readBinaryMat(btReducedDeformableBody::tDenseMatrix& mat, 
																						 const unsigned int n_modes, 		// #modes, outer array size
																						 const unsigned int n_full, 		// inner array size
																						 const char* file)
{
	std::ifstream f_in(file, std::ios::in | std::ios::binary);
	// first get size
	unsigned int v_size=0;
	f_in.read((char*)&v_size, 4);//sizeof(unsigned int));
	btAssert(v_size >= n_modes * n_full); 	// make sure the #requested mode is smaller than the #available modes

	// read data
	mat.resize(n_modes);
	for (int i = 0; i < n_modes; ++i) 
	{
		for (int j = 0; j < n_full; ++j)
		{
			double temp;
			f_in.read((char*)&temp, sizeof(double));

			if (mat[i].size() != n_modes)
				mat[i].resize(n_full);
			mat[i][j] = btScalar(temp);
		}
	}
  f_in.close();
}

void btReducedDeformableBodyHelpers::calculateLocalInertia(btVector3& inertia, const btScalar mass, const btVector3& half_extents, const btVector3& margin)
{
	btScalar lx = btScalar(2.) * (half_extents[0] + margin[0]);
	btScalar ly = btScalar(2.) * (half_extents[1] + margin[1]);
	btScalar lz = btScalar(2.) * (half_extents[2] + margin[2]);

	inertia.setValue(mass / (btScalar(12.0)) * (ly * ly + lz * lz),
								   mass / (btScalar(12.0)) * (lx * lx + lz * lz),
								   mass / (btScalar(12.0)) * (lx * lx + ly * ly));
}
#ifndef BT_REDUCED_SOFT_BODY_HELPERS_H
#define BT_REDUCED_SOFT_BODY_HELPERS_H

#include "btReducedDeformableBody.h"
#include <string>

struct btReducedDeformableBodyHelpers
{
	// create a reduced deformable object
	static btReducedDeformableBody* createReducedDeformableObject(btSoftBodyWorldInfo& worldInfo, const std::string& file_path, const std::string& vtk_file, const int num_modes, bool rigid_only);
	// read in geometry info from Vtk file
  static btReducedDeformableBody* createFromVtkFile(btSoftBodyWorldInfo& worldInfo, const char* vtk_file);
	// read in all reduced files
	static void readReducedDeformableInfoFromFiles(btReducedDeformableBody* rsb, const char* file_path);
	// read in a binary vector
	static void readBinaryVec(btReducedDeformableBody::tDenseArray& vec, const unsigned int n_size, const char* file);
	// read in a binary matrix
	static void readBinaryMat(btReducedDeformableBody::tDenseMatrix& mat, const unsigned int n_modes, const unsigned int n_full, const char* file);
	
	// calculate the local inertia tensor for a box shape reduced deformable object
	static void calculateLocalInertia(btVector3& inertia, const btScalar mass, const btVector3& half_extents, const btVector3& margin);
};


#endif // BT_REDUCED_SOFT_BODY_HELPERS_H
#include "btReducedDeformableBodySolver.h"
#include "btReducedDeformableBody.h"

btReducedDeformableBodySolver::btReducedDeformableBodySolver()
{
	m_ascendOrder = true;
	m_reducedSolver = true;
	m_dampingAlpha = 0;
	m_dampingBeta = 0;
	m_gravity = btVector3(0, 0, 0);
}

void btReducedDeformableBodySolver::setGravity(const btVector3& gravity)
{
	m_gravity = gravity;
}

void btReducedDeformableBodySolver::reinitialize(const btAlignedObjectArray<btSoftBody*>& bodies, btScalar dt)
{
	m_softBodies.copyFromArray(bodies);
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

	int N = bodies.size();
	if (nodeUpdated)
	{
		m_staticConstraints.resize(N);
		m_nodeRigidConstraints.resize(N);
		// m_faceRigidConstraints.resize(N);
	}
	for (int i = 0; i < N; ++i)
	{
		m_staticConstraints[i].clear();
		m_nodeRigidConstraints[i].clear();
		// m_faceRigidConstraints[i].clear();
	}

	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btReducedDeformableBody* rsb = static_cast<btReducedDeformableBody*>(m_softBodies[i]);
		rsb->m_contactNodesList.clear();
	}

	// set node index offsets
	int sum = 0;
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btReducedDeformableBody* rsb = static_cast<btReducedDeformableBody*>(m_softBodies[i]);
		rsb->m_nodeIndexOffset = sum;
		sum += rsb->m_nodes.size();
	}

	btDeformableBodySolver::updateSoftBodies();
}

void btReducedDeformableBodySolver::predictMotion(btScalar solverdt)
{
	applyExplicitForce(solverdt);

	// predict new mesh location
	predictReduceDeformableMotion(solverdt);

	//TODO: check if there is anything missed from btDeformableBodySolver::predictDeformableMotion
}

void btReducedDeformableBodySolver::predictReduceDeformableMotion(btScalar solverdt)
{
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btReducedDeformableBody* rsb = static_cast<btReducedDeformableBody*>(m_softBodies[i]);
		if (!rsb->isActive() || rsb->isStaticObject())
		{
			continue;
		}

		// clear contacts variables
		rsb->m_nodeRigidContacts.resize(0);
		rsb->m_faceRigidContacts.resize(0);
		rsb->m_faceNodeContacts.resize(0);
		rsb->m_nodeNodeContacts.resize(0);

		// calculate inverse mass matrix for all nodes
		for (int j = 0; j < rsb->m_nodes.size(); ++j)
		{
			if (rsb->m_nodes[j].m_im > 0)
			{
				rsb->m_nodes[j].m_effectiveMass_inv = rsb->m_nodes[j].m_effectiveMass.inverse();
			}
		}

		// rigid motion: t, R at time^*
		rsb->predictIntegratedTransform(solverdt, rsb->getInterpolationWorldTransform());

		// update reduced dofs at time^*
		// rsb->updateReducedDofs(solverdt);

		// update local moment arm at time^*
		// rsb->updateLocalMomentArm();
		// rsb->updateExternalForceProjectMatrix(true);

		// predict full space velocity at time^* (needed for constraints)
		rsb->mapToFullVelocity(rsb->getInterpolationWorldTransform());

		// update full space nodal position at time^*
		rsb->mapToFullPosition(rsb->getInterpolationWorldTransform());

		// update bounding box
		rsb->updateBounds();

		// update tree
		rsb->updateNodeTree(true, true);
		if (!rsb->m_fdbvt.empty())
		{
			rsb->updateFaceTree(true, true);
		}
	}
}

void btReducedDeformableBodySolver::applyExplicitForce(btScalar solverdt)
{
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btReducedDeformableBody* rsb = static_cast<btReducedDeformableBody*>(m_softBodies[i]);

		// apply gravity to the rigid frame, get m_linearVelocity at time^*
		rsb->applyRigidGravity(m_gravity, solverdt);

		if (!rsb->isReducedModesOFF())
		{
			// add internal force (elastic force & damping force)
			rsb->applyReducedElasticForce(rsb->m_reducedDofsBuffer);
			rsb->applyReducedDampingForce(rsb->m_reducedVelocityBuffer);

			// get reduced velocity at time^*
			rsb->updateReducedVelocity(solverdt);
		}

		// apply damping (no need at this point)
		// rsb->applyDamping(solverdt);
	}
}

void btReducedDeformableBodySolver::applyTransforms(btScalar timeStep)
{
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btReducedDeformableBody* rsb = static_cast<btReducedDeformableBody*>(m_softBodies[i]);

		// rigid motion
		rsb->proceedToTransform(timeStep, true);

		if (!rsb->isReducedModesOFF())
		{
			// update reduced dofs for time^n+1
			rsb->updateReducedDofs(timeStep);

			// update local moment arm for time^n+1
			rsb->updateLocalMomentArm();
			rsb->updateExternalForceProjectMatrix(true);
		}

		// update mesh nodal positions for time^n+1
		rsb->mapToFullPosition(rsb->getRigidTransform());

		// update mesh nodal velocity
		rsb->mapToFullVelocity(rsb->getRigidTransform());

		// end of time step clean up and update
		rsb->endOfTimeStepZeroing();

		// update the rendering mesh
		rsb->interpolateRenderMesh();
	}
}

void btReducedDeformableBodySolver::setConstraints(const btContactSolverInfo& infoGlobal)
{
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btReducedDeformableBody* rsb = static_cast<btReducedDeformableBody*>(m_softBodies[i]);
		if (!rsb->isActive() || rsb->isStaticObject())
		{
			continue;
		}

		// set fixed constraints
		for (int j = 0; j < rsb->m_fixedNodes.size(); ++j)
		{
			int i_node = rsb->m_fixedNodes[j];
			if (rsb->m_nodes[i_node].m_im == 0)
			{
				for (int k = 0; k < 3; ++k)
				{
					btVector3 dir(0, 0, 0);
					dir[k] = 1;
					btReducedDeformableStaticConstraint static_constraint(rsb, &rsb->m_nodes[i_node], rsb->getRelativePos(i_node), rsb->m_x0[i_node], dir, infoGlobal, m_dt);
					m_staticConstraints[i].push_back(static_constraint);
				}
			}
		}
		btAssert(rsb->m_fixedNodes.size() * 3 == m_staticConstraints[i].size());

		// set Deformable Node vs. Rigid constraint
		for (int j = 0; j < rsb->m_nodeRigidContacts.size(); ++j)
		{
			const btSoftBody::DeformableNodeRigidContact& contact = rsb->m_nodeRigidContacts[j];
			// skip fixed points
			if (contact.m_node->m_im == 0)
			{
				continue;
			}
			btReducedDeformableNodeRigidContactConstraint constraint(rsb, contact, infoGlobal, m_dt);
			m_nodeRigidConstraints[i].push_back(constraint);
			rsb->m_contactNodesList.push_back(contact.m_node->index - rsb->m_nodeIndexOffset);
		}
		// std::cout << "contact node list size: " << rsb->m_contactNodesList.size() << "\n";
		// std::cout << "#contact nodes: " << m_nodeRigidConstraints[i].size() << "\n";
	}
}

btScalar btReducedDeformableBodySolver::solveContactConstraints(btCollisionObject** deformableBodies, int numDeformableBodies, const btContactSolverInfo& infoGlobal)
{
	btScalar residualSquare = 0;

	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btAlignedObjectArray<int> m_orderNonContactConstraintPool;
		btAlignedObjectArray<int> m_orderContactConstraintPool;

		btReducedDeformableBody* rsb = static_cast<btReducedDeformableBody*>(m_softBodies[i]);

		// shuffle the order of applying constraint
		m_orderNonContactConstraintPool.resize(m_staticConstraints[i].size());
		m_orderContactConstraintPool.resize(m_nodeRigidConstraints[i].size());
		if (infoGlobal.m_solverMode & SOLVER_RANDMIZE_ORDER)
		{
			// fixed constraint order
			for (int j = 0; j < m_staticConstraints[i].size(); ++j)
			{
				m_orderNonContactConstraintPool[j] = m_ascendOrder ? j : m_staticConstraints[i].size() - 1 - j;
			}
			// contact constraint order
			for (int j = 0; j < m_nodeRigidConstraints[i].size(); ++j)
			{
				m_orderContactConstraintPool[j] = m_ascendOrder ? j : m_nodeRigidConstraints[i].size() - 1 - j;
			}

			m_ascendOrder = m_ascendOrder ? false : true;
		}
		else
		{
			for (int j = 0; j < m_staticConstraints[i].size(); ++j)
			{
				m_orderNonContactConstraintPool[j] = j;
			}
			// contact constraint order
			for (int j = 0; j < m_nodeRigidConstraints[i].size(); ++j)
			{
				m_orderContactConstraintPool[j] = j;
			}
		}

		// handle fixed constraint
		for (int k = 0; k < m_staticConstraints[i].size(); ++k)
		{
			btReducedDeformableStaticConstraint& constraint = m_staticConstraints[i][m_orderNonContactConstraintPool[k]];
			btScalar localResidualSquare = constraint.solveConstraint(infoGlobal);
			residualSquare = btMax(residualSquare, localResidualSquare);
		}

		// handle contact constraint

		// node vs rigid contact
		// std::cout << "!!#contact_nodes: " << m_nodeRigidConstraints[i].size() << '\n';
		for (int k = 0; k < m_nodeRigidConstraints[i].size(); ++k)
		{
			btReducedDeformableNodeRigidContactConstraint& constraint = m_nodeRigidConstraints[i][m_orderContactConstraintPool[k]];
			btScalar localResidualSquare = constraint.solveConstraint(infoGlobal);
			residualSquare = btMax(residualSquare, localResidualSquare);
		}

		// face vs rigid contact
		// for (int k = 0; k < m_faceRigidConstraints[i].size(); ++k)
		// {
		// 	btReducedDeformableFaceRigidContactConstraint& constraint = m_faceRigidConstraints[i][k];
		// 	btScalar localResidualSquare = constraint.solveConstraint(infoGlobal);
		// 	residualSquare = btMax(residualSquare, localResidualSquare);
		// }
	}

	return residualSquare;
}

void btReducedDeformableBodySolver::deformableBodyInternalWriteBack()
{
	// reduced deformable update
	for (int i = 0; i < m_softBodies.size(); ++i)
	{
		btReducedDeformableBody* rsb = static_cast<btReducedDeformableBody*>(m_softBodies[i]);
		rsb->applyInternalVelocityChanges();
	}
	m_ascendOrder = true;
}
#ifndef BT_REDUCED_DEFORMABLE_BODY_DYNAMICS_WORLD_H
#define BT_REDUCED_DEFORMABLE_BODY_DYNAMICS_WORLD_H

#include "BulletSoftBody/btDeformableBodySolver.h"
#include "btReducedDeformableContactConstraint.h"

class btReducedDeformableBody;

class btReducedDeformableBodySolver : public btDeformableBodySolver
{
 protected:
  bool m_ascendOrder;
  btScalar m_dampingAlpha;
  btScalar m_dampingBeta;

  btVector3 m_gravity;

  void predictReduceDeformableMotion(btScalar solverdt);

  void applyExplicitForce(btScalar solverdt);

 public:
  btAlignedObjectArray<btAlignedObjectArray<btReducedDeformableStaticConstraint> > m_staticConstraints;
  btAlignedObjectArray<btAlignedObjectArray<btReducedDeformableNodeRigidContactConstraint> > m_nodeRigidConstraints;
  btAlignedObjectArray<btAlignedObjectArray<btReducedDeformableFaceRigidContactConstraint> > m_faceRigidConstraints;
  
  btReducedDeformableBodySolver();
  ~btReducedDeformableBodySolver() {}

  virtual void setGravity(const btVector3& gravity);

  virtual SolverTypes getSolverType() const
  {
    return REDUCED_DEFORMABLE_SOLVER;
  }

  // resize/clear data structures
	virtual void reinitialize(const btAlignedObjectArray<btSoftBody*>& bodies, btScalar dt);

  virtual void predictMotion(btScalar solverdt);

  virtual void applyTransforms(btScalar timeStep);

  // set up contact constraints
	virtual void setConstraints(const btContactSolverInfo& infoGlobal);

  // solve all constraints (fixed and contact)
  virtual btScalar solveContactConstraints(btCollisionObject** deformableBodies, int numDeformableBodies, const btContactSolverInfo& infoGlobal);

  // apply all the delta velocities
  virtual void deformableBodyInternalWriteBack();

  // virtual void setProjection() {}

  // virtual void setLagrangeMultiplier() {}

  // virtual void setupDeformableSolve(bool implicit);

};

#endif // BT_REDUCED_DEFORMABLE_BODY_DYNAMICS_WORLD_H
#include "btReducedDeformableContactConstraint.h"
#include <iostream>

// ================= static constraints ===================
btReducedDeformableStaticConstraint::btReducedDeformableStaticConstraint(
  btReducedDeformableBody* rsb, 
  btSoftBody::Node* node,
	const btVector3& ri,
	const btVector3& x0,
	const btVector3& dir,
  const btContactSolverInfo& infoGlobal,
	btScalar dt)
  : m_rsb(rsb), m_ri(ri), m_targetPos(x0), m_impulseDirection(dir), m_dt(dt), btDeformableStaticConstraint(node, infoGlobal)
{
	m_erp = 0.2;
	m_appliedImpulse = 0;

	// get impulse factor
  m_impulseFactorMatrix = rsb->getImpulseFactor(m_node->index);
	m_impulseFactor = (m_impulseFactorMatrix * m_impulseDirection).dot(m_impulseDirection);

	btScalar vel_error = btDot(-m_node->m_v, m_impulseDirection);
	btScalar pos_error = btDot(m_targetPos - m_node->m_x, m_impulseDirection);

	m_rhs = (vel_error + m_erp * pos_error / m_dt) / m_impulseFactor;
}

btScalar btReducedDeformableStaticConstraint::solveConstraint(const btContactSolverInfo& infoGlobal)
{
	// target velocity of fixed constraint is 0
	btVector3 deltaVa = getDeltaVa();
	btScalar deltaV_rel = btDot(deltaVa, m_impulseDirection);
  btScalar deltaImpulse = m_rhs - deltaV_rel / m_impulseFactor;
	m_appliedImpulse = m_appliedImpulse + deltaImpulse;

	btVector3 impulse = deltaImpulse * m_impulseDirection;
	applyImpulse(impulse);

	// calculate residual
	btScalar residualSquare = m_impulseFactor * deltaImpulse;
	residualSquare *= residualSquare;

	return residualSquare;
}
  
// this calls reduced deformable body's internalApplyFullSpaceImpulse
void btReducedDeformableStaticConstraint::applyImpulse(const btVector3& impulse)
{
	// apply full space impulse
	m_rsb->internalApplyFullSpaceImpulse(impulse, m_ri, m_node->index, m_dt);
}

btVector3 btReducedDeformableStaticConstraint::getDeltaVa() const
{
	return m_rsb->internalComputeNodeDeltaVelocity(m_rsb->getInterpolationWorldTransform(), m_node->index);
}

// ================= base contact constraints ===================
btReducedDeformableRigidContactConstraint::btReducedDeformableRigidContactConstraint(
  btReducedDeformableBody* rsb, 
  const btSoftBody::DeformableRigidContact& c, 
  const btContactSolverInfo& infoGlobal,
	btScalar dt)
  : m_rsb(rsb), m_dt(dt), btDeformableRigidContactConstraint(c, infoGlobal)
{
	m_nodeQueryIndex = 0;
	m_appliedNormalImpulse = 0;
  m_appliedTangentImpulse = 0;
	m_rhs = 0;
	m_rhs_tangent = 0;
	m_cfm = infoGlobal.m_deformable_cfm;
	m_cfm_friction = 0;
	m_erp = infoGlobal.m_deformable_erp;
	m_erp_friction = infoGlobal.m_deformable_erp;
	m_friction = infoGlobal.m_friction;

	m_collideStatic = m_contact->m_cti.m_colObj->isStaticObject();
	m_collideMultibody = (m_contact->m_cti.m_colObj->getInternalType() == btCollisionObject::CO_FEATHERSTONE_LINK);
}

void btReducedDeformableRigidContactConstraint::setSolverBody(const int bodyId, btSolverBody& solver_body)
{
	if (!m_collideMultibody)
	{
		m_solverBodyId = bodyId;
		m_solverBody = &solver_body;
		m_linearComponentNormal = -m_contactNormalA * m_solverBody->internalGetInvMass();
		btVector3	torqueAxis = -m_relPosA.cross(m_contactNormalA);
		m_angularComponentNormal = m_solverBody->m_originalBody->getInvInertiaTensorWorld() * torqueAxis;
		
		m_linearComponentTangent = m_contactTangent * m_solverBody->internalGetInvMass();
		btVector3 torqueAxisTangent = m_relPosA.cross(m_contactTangent);
		m_angularComponentTangent = m_solverBody->m_originalBody->getInvInertiaTensorWorld() * torqueAxisTangent;
	}
}

btVector3 btReducedDeformableRigidContactConstraint::getVa() const
{
	btVector3 Va(0, 0, 0);
	if (!m_collideStatic)
	{
		Va = btDeformableRigidContactConstraint::getVa();
	}
	return Va;
}

btScalar btReducedDeformableRigidContactConstraint::solveConstraint(const btContactSolverInfo& infoGlobal)
{
	// btVector3 Va = getVa();
	// btVector3 deltaVa = Va - m_bufferVelocityA;
	// if (!m_collideStatic)
	// {
		// std::cout << "moving collision!!!\n";
		// std::cout << "relPosA: " << m_relPosA[0] << "\t" << m_relPosA[1] << "\t" << m_relPosA[2] << "\n";
		// std::cout << "moving rigid linear_vel: " << m_solverBody->m_originalBody->getLinearVelocity()[0] << '\t'
		//  << m_solverBody->m_originalBody->getLinearVelocity()[1] << '\t'
		//   << m_solverBody->m_originalBody->getLinearVelocity()[2] << '\n';
	// }
	btVector3 deltaVa = getDeltaVa();
	btVector3 deltaVb = getDeltaVb();

	// if (!m_collideStatic)
	// {
	// 	std::cout << "deltaVa: " << deltaVa[0] << '\t' << deltaVa[1] << '\t' << deltaVa[2] << '\n';
	// 	std::cout << "deltaVb: " << deltaVb[0] << '\t' << deltaVb[1] << '\t' << deltaVb[2] << '\n';
	// }

	// get delta relative velocity and magnitude (i.e., how much impulse has been applied?)
	btVector3 deltaV_rel = deltaVa - deltaVb;
	btScalar deltaV_rel_normal = -btDot(deltaV_rel, m_contactNormalA);

	// if (!m_collideStatic)
	// {
	// 	std::cout << "deltaV_rel: " << deltaV_rel[0] << '\t' << deltaV_rel[1] << '\t' << deltaV_rel[2] << "\n";
	// 	std::cout << "deltaV_rel_normal: " << deltaV_rel_normal << "\n";
	// 	std::cout << "normal_A: " << m_contactNormalA[0] << '\t' << m_contactNormalA[1] << '\t' << m_contactNormalA[2] << '\n';
	// }
	
	// get the normal impulse to be applied
	btScalar deltaImpulse = m_rhs - m_appliedNormalImpulse * m_cfm - deltaV_rel_normal / m_normalImpulseFactor;
	// if (!m_collideStatic)
	// {
	// 	std::cout << "m_rhs: " << m_rhs << '\t' << "m_appliedNormalImpulse: "  << m_appliedNormalImpulse << "\n";
	// 	std::cout << "m_normalImpulseFactor: " << m_normalImpulseFactor << '\n';
	// }

	{
		// cumulative impulse that has been applied
		btScalar sum = m_appliedNormalImpulse + deltaImpulse;
		// if the cumulative impulse is pushing the object into the rigid body, set it zero
		if (sum < 0)
		{
			deltaImpulse = -m_appliedNormalImpulse;
			m_appliedNormalImpulse = 0;
		}
		else
		{
			m_appliedNormalImpulse = sum;
		}	
	}

	// if (!m_collideStatic)
	// {
	// 	std::cout << "m_appliedNormalImpulse: " << m_appliedNormalImpulse << '\n';
	// 	std::cout << "deltaImpulse: " << deltaImpulse << '\n';
	// }

	// residual is the nodal normal velocity change in current iteration
	btScalar residualSquare = deltaImpulse * m_normalImpulseFactor;	// get residual
	residualSquare *= residualSquare;

	
	// apply Coulomb friction (based on delta velocity, |dv_t| = |dv_n * friction|)
	btScalar deltaImpulse_tangent = 0;
	btScalar deltaImpulse_tangent2 = 0;
	{
		// calculate how much impulse is needed
		// btScalar deltaV_rel_tangent = btDot(deltaV_rel, m_contactTangent);
		// btScalar impulse_changed = deltaV_rel_tangent * m_tangentImpulseFactorInv;
		// deltaImpulse_tangent = m_rhs_tangent - impulse_changed;

		// btScalar sum = m_appliedTangentImpulse + deltaImpulse_tangent;
		btScalar lower_limit = - m_appliedNormalImpulse * m_friction;
		btScalar upper_limit = m_appliedNormalImpulse * m_friction;
		// if (sum > upper_limit)
		// {
		// 	deltaImpulse_tangent = upper_limit - m_appliedTangentImpulse;
		// 	m_appliedTangentImpulse = upper_limit;
		// }
		// else if (sum < lower_limit)
		// {
		// 	deltaImpulse_tangent = lower_limit - m_appliedTangentImpulse;
		// 	m_appliedTangentImpulse = lower_limit;
		// }
		// else
		// {
		// 	m_appliedTangentImpulse = sum;
		// }
		// 
		calculateTangentialImpulse(deltaImpulse_tangent, m_appliedTangentImpulse, m_rhs_tangent,
															 m_tangentImpulseFactorInv, m_contactTangent, lower_limit, upper_limit, deltaV_rel);
		
		if (m_collideMultibody)
		{
			calculateTangentialImpulse(deltaImpulse_tangent2, m_appliedTangentImpulse2, m_rhs_tangent2,
															   m_tangentImpulseFactorInv2, m_contactTangent2, lower_limit, upper_limit, deltaV_rel);
		}
															 

		if (!m_collideStatic)
		{
			// std::cout << "m_contactTangent: " << m_contactTangent[0] << "\t"  << m_contactTangent[1] << "\t"  << m_contactTangent[2] << "\n";
			// std::cout << "deltaV_rel_tangent: " << deltaV_rel_tangent << '\n';
			// std::cout << "deltaImpulseTangent: " << deltaImpulse_tangent << '\n';
			// std::cout << "m_appliedTangentImpulse: " << m_appliedTangentImpulse << '\n';
		}
	}

	// get the total impulse vector
	btVector3 impulse_normal = deltaImpulse * m_contactNormalA;
	btVector3 impulse_tangent = deltaImpulse_tangent * (-m_contactTangent);
	btVector3 impulse_tangent2 = deltaImpulse_tangent2 * (-m_contactTangent2);
	btVector3 impulse = impulse_normal + impulse_tangent + impulse_tangent2;

	applyImpulse(impulse);
	
	// apply impulse to the rigid/multibodies involved and change their velocities
	if (!m_collideStatic)
	{
		// std::cout << "linear_component: " << m_linearComponentNormal[0] << '\t'
		// 																	<< m_linearComponentNormal[1] << '\t'
		// 																	<< m_linearComponentNormal[2] << '\n';
		// std::cout << "angular_component: " << m_angularComponentNormal[0] << '\t'
		// 																	<< m_angularComponentNormal[1] << '\t'
		// 																	<< m_angularComponentNormal[2] << '\n';

		if (!m_collideMultibody)		// collision with rigid body
		{
			// std::cout << "rigid impulse applied!!\n";
			// std::cout << "delta Linear: " << m_solverBody->getDeltaLinearVelocity()[0] << '\t'
			// << m_solverBody->getDeltaLinearVelocity()[1] << '\t'
			// 	<< m_solverBody->getDeltaLinearVelocity()[2] << '\n';
			// std::cout << "delta Angular: " << m_solverBody->getDeltaAngularVelocity()[0] << '\t'
			// << m_solverBody->getDeltaAngularVelocity()[1] << '\t'
			// 	<< m_solverBody->getDeltaAngularVelocity()[2] << '\n';

			m_solverBody->internalApplyImpulse(m_linearComponentNormal, m_angularComponentNormal, deltaImpulse);
			m_solverBody->internalApplyImpulse(m_linearComponentTangent, m_angularComponentTangent, deltaImpulse_tangent);

			// std::cout << "after\n";
			// std::cout << "rigid impulse applied!!\n";
			// std::cout << "delta Linear: " << m_solverBody->getDeltaLinearVelocity()[0] << '\t'
			// << m_solverBody->getDeltaLinearVelocity()[1] << '\t'
			// 	<< m_solverBody->getDeltaLinearVelocity()[2] << '\n';
			// std::cout << "delta Angular: " << m_solverBody->getDeltaAngularVelocity()[0] << '\t'
			// << m_solverBody->getDeltaAngularVelocity()[1] << '\t'
			// 	<< m_solverBody->getDeltaAngularVelocity()[2] << '\n';
		}
		else		// collision with multibody
		{
			btMultiBodyLinkCollider* multibodyLinkCol = 0;
			multibodyLinkCol = (btMultiBodyLinkCollider*)btMultiBodyLinkCollider::upcast(m_contact->m_cti.m_colObj);
			if (multibodyLinkCol)
			{
				const btScalar* deltaV_normal = &m_contact->jacobianData_normal.m_deltaVelocitiesUnitImpulse[0];
				// apply normal component of the impulse
				multibodyLinkCol->m_multiBody->applyDeltaVeeMultiDof2(deltaV_normal, -deltaImpulse);
				
				// const int ndof = multibodyLinkCol->m_multiBody->getNumDofs() + 6;
				// std::cout << "deltaV_normal: ";
				// for (int i = 0; i < ndof; ++i)
				// {
				// 	std::cout << i << "\t" << deltaV_normal[i] << '\n';
				// }

				if (impulse_tangent.norm() > SIMD_EPSILON)
				{
					// apply tangential component of the impulse
					const btScalar* deltaV_t1 = &m_contact->jacobianData_t1.m_deltaVelocitiesUnitImpulse[0];
					multibodyLinkCol->m_multiBody->applyDeltaVeeMultiDof2(deltaV_t1, deltaImpulse_tangent);
					const btScalar* deltaV_t2 = &m_contact->jacobianData_t2.m_deltaVelocitiesUnitImpulse[0];
					multibodyLinkCol->m_multiBody->applyDeltaVeeMultiDof2(deltaV_t2, deltaImpulse_tangent2);
				}
			}
		}
	}
	return residualSquare;
}

void btReducedDeformableRigidContactConstraint::calculateTangentialImpulse(btScalar& deltaImpulse_tangent, 
																		 																			 btScalar& appliedImpulse, 
																																					 const btScalar rhs_tangent,
																																					 const btScalar tangentImpulseFactorInv,
																																					 const btVector3& tangent,
																		 																			 const btScalar lower_limit,
																																					 const btScalar upper_limit,
																																					 const btVector3& deltaV_rel)
{
	btScalar deltaV_rel_tangent = btDot(deltaV_rel, tangent);
	btScalar impulse_changed = deltaV_rel_tangent * tangentImpulseFactorInv;
	deltaImpulse_tangent = rhs_tangent - m_cfm_friction * appliedImpulse - impulse_changed;

	btScalar sum = appliedImpulse + deltaImpulse_tangent;
	if (sum > upper_limit)
	{
		deltaImpulse_tangent = upper_limit - appliedImpulse;
		appliedImpulse = upper_limit;
	}
	else if (sum < lower_limit)
	{
		deltaImpulse_tangent = lower_limit - appliedImpulse;
		appliedImpulse = lower_limit;
	}
	else
	{
		appliedImpulse = sum;
	}
}

// ================= node vs rigid constraints ===================
btReducedDeformableNodeRigidContactConstraint::btReducedDeformableNodeRigidContactConstraint(
  btReducedDeformableBody* rsb, 
  const btSoftBody::DeformableNodeRigidContact& contact, 
  const btContactSolverInfo& infoGlobal,
	btScalar dt)
  : m_node(contact.m_node), btReducedDeformableRigidContactConstraint(rsb, contact, infoGlobal, dt)
{
	m_contactNormalA = contact.m_cti.m_normal;
  m_contactNormalB = -contact.m_cti.m_normal;

	if (contact.m_node->index < rsb->m_nodes.size())
	{
		m_nodeQueryIndex = contact.m_node->index;
	}
	else
	{
		m_nodeQueryIndex = m_node->index - rsb->m_nodeIndexOffset;
	}

	if (m_contact->m_cti.m_colObj->getInternalType() == btCollisionObject::CO_RIGID_BODY)
	{
		m_relPosA = contact.m_c1;
	}
	else
	{
		m_relPosA = btVector3(0,0,0);
	}
	m_relPosB = m_node->m_x - m_rsb->getRigidTransform().getOrigin();

	if (m_collideStatic)		// colliding with static object, only consider reduced deformable body's impulse factor
	{
		m_impulseFactor = m_rsb->getImpulseFactor(m_nodeQueryIndex);
	}
	else		// colliding with dynamic object, consider both reduced deformable and rigid body's impulse factors
	{
		m_impulseFactor = m_rsb->getImpulseFactor(m_nodeQueryIndex) + contact.m_c0;
	}

	m_normalImpulseFactor = (m_impulseFactor * m_contactNormalA).dot(m_contactNormalA);
	m_tangentImpulseFactor = 0;

	warmStarting();
}

void btReducedDeformableNodeRigidContactConstraint::warmStarting()
{
	btVector3 va = getVa();
	btVector3 vb = getVb();
	m_bufferVelocityA = va;
	m_bufferVelocityB = vb;

	// we define the (+) direction of errors to be the outward surface normal of the rigid object
	btVector3 v_rel = vb - va;
	// get tangent direction of the relative velocity
	btVector3 v_tangent = v_rel - v_rel.dot(m_contactNormalA) * m_contactNormalA;
	if (v_tangent.norm() < SIMD_EPSILON)
	{
		m_contactTangent = btVector3(0, 0, 0);
		// tangent impulse factor
		m_tangentImpulseFactor = 0;
		m_tangentImpulseFactorInv = 0;
	}
	else
	{
		if (!m_collideMultibody)
		{
			m_contactTangent = v_tangent.normalized();
			m_contactTangent2.setZero();
			// tangent impulse factor 1
			m_tangentImpulseFactor = (m_impulseFactor * m_contactTangent).dot(m_contactTangent);
			m_tangentImpulseFactorInv = btScalar(1) / m_tangentImpulseFactor;
			// tangent impulse factor 2
			m_tangentImpulseFactor2 = 0;
			m_tangentImpulseFactorInv2 = 0;
		}
		else	// multibody requires 2 contact directions
		{
			m_contactTangent = m_contact->t1;
			m_contactTangent2 = m_contact->t2;

			// tangent impulse factor 1
			m_tangentImpulseFactor = (m_impulseFactor * m_contactTangent).dot(m_contactTangent);
			m_tangentImpulseFactorInv = btScalar(1) / m_tangentImpulseFactor;
			// tangent impulse factor 2
			m_tangentImpulseFactor2 = (m_impulseFactor * m_contactTangent2).dot(m_contactTangent2);
			m_tangentImpulseFactorInv2 = btScalar(1) / m_tangentImpulseFactor2;
		}
	}


	// initial guess for normal impulse
	{
		btScalar velocity_error = btDot(v_rel, m_contactNormalA);	// magnitude of relative velocity
		btScalar position_error = 0;
		if (m_penetration > 0)
		{
			velocity_error += m_penetration / m_dt;
		}
		else
		{
			// add penetration correction vel
			position_error = m_penetration * m_erp / m_dt;
		}
		// get the initial estimate of impulse magnitude to be applied
		m_rhs = -(velocity_error + position_error) / m_normalImpulseFactor;
	}

	// initial guess for tangential impulse
	{
		btScalar velocity_error = btDot(v_rel, m_contactTangent);
		m_rhs_tangent = velocity_error * m_tangentImpulseFactorInv;

		if (m_collideMultibody)
		{
			btScalar velocity_error2 = btDot(v_rel, m_contactTangent2);
			m_rhs_tangent2 = velocity_error2 * m_tangentImpulseFactorInv2;
		}
	}

	// warm starting
	// applyImpulse(m_rhs * m_contactNormalA);
	// if (!m_collideStatic)
	// {
	// 	const btSoftBody::sCti& cti = m_contact->m_cti;
	// 	if (cti.m_colObj->getInternalType() == btCollisionObject::CO_RIGID_BODY)
	// 	{
	// 		m_solverBody->internalApplyImpulse(m_linearComponentNormal, m_angularComponentNormal, -m_rhs);
	// 	}
	// }
}

btVector3 btReducedDeformableNodeRigidContactConstraint::getVb() const
{
	return m_node->m_v;
}

btVector3 btReducedDeformableNodeRigidContactConstraint::getDeltaVa() const
{
	btVector3 deltaVa(0, 0, 0);
	if (!m_collideStatic)
	{
		if (!m_collideMultibody)		// for rigid body
		{
			deltaVa = m_solverBody->internalGetDeltaLinearVelocity() + m_solverBody->internalGetDeltaAngularVelocity().cross(m_relPosA);
		}
		else		// for multibody
		{
			btMultiBodyLinkCollider* multibodyLinkCol = 0;
			multibodyLinkCol = (btMultiBodyLinkCollider*)btMultiBodyLinkCollider::upcast(m_contact->m_cti.m_colObj);
			if (multibodyLinkCol)
			{
				const int ndof = multibodyLinkCol->m_multiBody->getNumDofs() + 6;
				const btScalar* J_n = &m_contact->jacobianData_normal.m_jacobians[0];
				const btScalar* J_t1 = &m_contact->jacobianData_t1.m_jacobians[0];
				const btScalar* J_t2 = &m_contact->jacobianData_t2.m_jacobians[0];
				const btScalar* local_dv = multibodyLinkCol->m_multiBody->getDeltaVelocityVector();
				// add in the normal component of the va
				btScalar vel = 0;
				for (int k = 0; k < ndof; ++k)
				{
					vel += local_dv[k] * J_n[k];
				}
				deltaVa = m_contact->m_cti.m_normal * vel;
				
				// add in the tangential components of the va
				vel = 0;
				for (int k = 0; k < ndof; ++k)
				{
					vel += local_dv[k] * J_t1[k];
				}
				deltaVa += m_contact->t1 * vel;

				vel = 0;
				for (int k = 0; k < ndof; ++k)
				{
					vel += local_dv[k] * J_t2[k];
				}
				deltaVa += m_contact->t2 * vel;
			}
		}
	}
	return deltaVa;
}

btVector3 btReducedDeformableNodeRigidContactConstraint::getDeltaVb() const
{	
	// std::cout << "node: " << m_node->index << '\n';
	return m_rsb->internalComputeNodeDeltaVelocity(m_rsb->getInterpolationWorldTransform(), m_nodeQueryIndex);
}

btVector3 btReducedDeformableNodeRigidContactConstraint::getSplitVb() const
{
	return m_node->m_splitv;
}

btVector3 btReducedDeformableNodeRigidContactConstraint::getDv(const btSoftBody::Node* node) const
{
	return m_total_normal_dv + m_total_tangent_dv;
}

void btReducedDeformableNodeRigidContactConstraint::applyImpulse(const btVector3& impulse)
{
  m_rsb->internalApplyFullSpaceImpulse(impulse, m_relPosB, m_nodeQueryIndex, m_dt);
	// m_rsb->applyFullSpaceImpulse(impulse, m_relPosB, m_node->index, m_dt);
	// m_rsb->mapToFullVelocity(m_rsb->getInterpolationWorldTransform());
	// if (!m_collideStatic)
	// {
	// 	// std::cout << "impulse applied: " << impulse[0] << '\t' << impulse[1] << '\t' << impulse[2] << '\n';
	// 	// std::cout << "node: " << m_node->index << " vel: " << m_node->m_v[0] << '\t' << m_node->m_v[1] << '\t' << m_node->m_v[2] << '\n';
	// 	btVector3 v_after = getDeltaVb() + m_node->m_v;
	// 	// std::cout << "vel after: " << v_after[0] << '\t' << v_after[1] << '\t' << v_after[2] << '\n';
	// }
	// std::cout << "node: " << m_node->index << " pos: " << m_node->m_x[0] << '\t' << m_node->m_x[1] << '\t' << m_node->m_x[2] << '\n';
}

// ================= face vs rigid constraints ===================
btReducedDeformableFaceRigidContactConstraint::btReducedDeformableFaceRigidContactConstraint(
  btReducedDeformableBody* rsb, 
  const btSoftBody::DeformableFaceRigidContact& contact, 
  const btContactSolverInfo& infoGlobal,
	btScalar dt, 
  bool useStrainLimiting)
  : m_face(contact.m_face), m_useStrainLimiting(useStrainLimiting), btReducedDeformableRigidContactConstraint(rsb, contact, infoGlobal, dt)
{}

btVector3 btReducedDeformableFaceRigidContactConstraint::getVb() const
{
	const btSoftBody::DeformableFaceRigidContact* contact = getContact();
	btVector3 vb = m_face->m_n[0]->m_v * contact->m_bary[0] + m_face->m_n[1]->m_v * contact->m_bary[1] + m_face->m_n[2]->m_v * contact->m_bary[2];
	return vb;
}

btVector3 btReducedDeformableFaceRigidContactConstraint::getSplitVb() const
{
	const btSoftBody::DeformableFaceRigidContact* contact = getContact();
	btVector3 vb = (m_face->m_n[0]->m_splitv) * contact->m_bary[0] + (m_face->m_n[1]->m_splitv) * contact->m_bary[1] + (m_face->m_n[2]->m_splitv) * contact->m_bary[2];
	return vb;
}

btVector3 btReducedDeformableFaceRigidContactConstraint::getDv(const btSoftBody::Node* node) const
{
	btVector3 face_dv = m_total_normal_dv + m_total_tangent_dv;
	const btSoftBody::DeformableFaceRigidContact* contact = getContact();
	if (m_face->m_n[0] == node)
	{
		return face_dv * contact->m_weights[0];
	}
	if (m_face->m_n[1] == node)
	{
		return face_dv * contact->m_weights[1];
	}
	btAssert(node == m_face->m_n[2]);
	return face_dv * contact->m_weights[2];
}

void btReducedDeformableFaceRigidContactConstraint::applyImpulse(const btVector3& impulse)
{
  //
}
#include "../btDeformableContactConstraint.h"
#include "btReducedDeformableBody.h"

// ================= static constraints ===================
class btReducedDeformableStaticConstraint : public btDeformableStaticConstraint
{
 public:
  btReducedDeformableBody* m_rsb;
  btScalar m_dt;
  btVector3 m_ri;
  btVector3 m_targetPos;
  btVector3 m_impulseDirection;
  btMatrix3x3 m_impulseFactorMatrix;
  btScalar m_impulseFactor;
  btScalar m_rhs;
  btScalar m_appliedImpulse;
  btScalar m_erp;

  btReducedDeformableStaticConstraint(btReducedDeformableBody* rsb, 
                                      btSoftBody::Node* node,
                                      const btVector3& ri,
                                      const btVector3& x0,
                                      const btVector3& dir,
                                      const btContactSolverInfo& infoGlobal,
                                      btScalar dt);
	// btReducedDeformableStaticConstraint(const btReducedDeformableStaticConstraint& other);
  btReducedDeformableStaticConstraint() {}
  virtual ~btReducedDeformableStaticConstraint() {}

  virtual btScalar solveConstraint(const btContactSolverInfo& infoGlobal);
  
  // this calls reduced deformable body's applyFullSpaceImpulse
  virtual void applyImpulse(const btVector3& impulse);

  btVector3 getDeltaVa() const;

  // virtual void applySplitImpulse(const btVector3& impulse) {}
};

// ================= base contact constraints ===================
class btReducedDeformableRigidContactConstraint : public btDeformableRigidContactConstraint
{
 public:
  bool m_collideStatic;     // flag for collision with static object
  bool m_collideMultibody;  // flag for collision with multibody

  int m_nodeQueryIndex;
  int m_solverBodyId;       // for debugging

  btReducedDeformableBody* m_rsb;
  btSolverBody* m_solverBody;
  btScalar m_dt;

  btScalar m_appliedNormalImpulse;
  btScalar m_appliedTangentImpulse;
  btScalar m_appliedTangentImpulse2;
  btScalar m_normalImpulseFactor;
  btScalar m_tangentImpulseFactor;
  btScalar m_tangentImpulseFactor2;
  btScalar m_tangentImpulseFactorInv;
  btScalar m_tangentImpulseFactorInv2;
  btScalar m_rhs;
  btScalar m_rhs_tangent;
  btScalar m_rhs_tangent2;
  
  btScalar m_cfm;
  btScalar m_cfm_friction;
  btScalar m_erp;
  btScalar m_erp_friction;
  btScalar m_friction;

  btVector3 m_contactNormalA;     // surface normal for rigid body (opposite direction as impulse)
  btVector3 m_contactNormalB;     // surface normal for reduced deformable body (opposite direction as impulse)
  btVector3 m_contactTangent;     // tangential direction of the relative velocity
  btVector3 m_contactTangent2;    // 2nd tangential direction of the relative velocity
  btVector3 m_relPosA;            // relative position of the contact point for A (rigid)
  btVector3 m_relPosB;            // relative position of the contact point for B
  btMatrix3x3 m_impulseFactor;    // total impulse matrix

  btVector3 m_bufferVelocityA;    // velocity at the beginning of the iteration
  btVector3 m_bufferVelocityB;
  btVector3 m_linearComponentNormal;    // linear components for the solver body
  btVector3 m_angularComponentNormal;   // angular components for the solver body
  // since 2nd contact direction only applies to multibody, these components will never be used
  btVector3 m_linearComponentTangent;
  btVector3 m_angularComponentTangent;

  btReducedDeformableRigidContactConstraint(btReducedDeformableBody* rsb, 
                                            const btSoftBody::DeformableRigidContact& c, 
                                            const btContactSolverInfo& infoGlobal,
                                            btScalar dt);
	// btReducedDeformableRigidContactConstraint(const btReducedDeformableRigidContactConstraint& other);
  btReducedDeformableRigidContactConstraint() {}
  virtual ~btReducedDeformableRigidContactConstraint() {}

  void setSolverBody(const int bodyId, btSolverBody& solver_body);

  virtual void warmStarting() {}

  virtual btScalar solveConstraint(const btContactSolverInfo& infoGlobal);

  void calculateTangentialImpulse(btScalar& deltaImpulse_tangent, 
                                  btScalar& appliedImpulse, 
                                  const btScalar rhs_tangent,
                                  const btScalar tangentImpulseFactorInv,
                                  const btVector3& tangent,
                                  const btScalar lower_limit,
                                  const btScalar upper_limit,
                                  const btVector3& deltaV_rel);

  virtual void applyImpulse(const btVector3& impulse) {}

  virtual void applySplitImpulse(const btVector3& impulse) {} // TODO: may need later

  virtual btVector3 getVa() const;
  virtual btVector3 getDeltaVa() const = 0;
  virtual btVector3 getDeltaVb() const = 0;
};

// ================= node vs rigid constraints ===================
class btReducedDeformableNodeRigidContactConstraint : public btReducedDeformableRigidContactConstraint
{
 public:
  btSoftBody::Node* m_node;

  btReducedDeformableNodeRigidContactConstraint(btReducedDeformableBody* rsb, 
                                                const btSoftBody::DeformableNodeRigidContact& contact, 
                                                const btContactSolverInfo& infoGlobal,
                                                btScalar dt);
	// btReducedDeformableNodeRigidContactConstraint(const btReducedDeformableNodeRigidContactConstraint& other);
  btReducedDeformableNodeRigidContactConstraint() {}
  virtual ~btReducedDeformableNodeRigidContactConstraint() {}

  virtual void warmStarting();

  // get the velocity of the deformable node in contact
	virtual btVector3 getVb() const;

  // get the velocity change of the rigid body
  virtual btVector3 getDeltaVa() const;

  // get velocity change of the node in contat
  virtual btVector3 getDeltaVb() const;

	// get the split impulse velocity of the deformable face at the contact point
	virtual btVector3 getSplitVb() const;

	// get the velocity change of the input soft body node in the constraint
	virtual btVector3 getDv(const btSoftBody::Node*) const;

	// cast the contact to the desired type
	const btSoftBody::DeformableNodeRigidContact* getContact() const
	{
		return static_cast<const btSoftBody::DeformableNodeRigidContact*>(m_contact);
	}
  
  // this calls reduced deformable body's applyFullSpaceImpulse
  virtual void applyImpulse(const btVector3& impulse);
};

// ================= face vs rigid constraints ===================
class btReducedDeformableFaceRigidContactConstraint : public btReducedDeformableRigidContactConstraint
{
 public:
  btSoftBody::Face* m_face;
	bool m_useStrainLimiting;

  btReducedDeformableFaceRigidContactConstraint(btReducedDeformableBody* rsb, 
                                                const btSoftBody::DeformableFaceRigidContact& contact, 
                                                const btContactSolverInfo& infoGlobal,
                                                btScalar dt, 
                                                bool useStrainLimiting);
	// btReducedDeformableFaceRigidContactConstraint(const btReducedDeformableFaceRigidContactConstraint& other);
  btReducedDeformableFaceRigidContactConstraint() {}
  virtual ~btReducedDeformableFaceRigidContactConstraint() {}

  // get the velocity of the deformable face at the contact point
	virtual btVector3 getVb() const;

	// get the split impulse velocity of the deformable face at the contact point
	virtual btVector3 getSplitVb() const;

	// get the velocity change of the input soft body node in the constraint
	virtual btVector3 getDv(const btSoftBody::Node*) const;

	// cast the contact to the desired type
	const btSoftBody::DeformableFaceRigidContact* getContact() const
	{
		return static_cast<const btSoftBody::DeformableFaceRigidContact*>(m_contact);
	}

  // this calls reduced deformable body's applyFullSpaceImpulse
  virtual void applyImpulse(const btVector3& impulse);
};
//
//  DeformableBodyInplaceSolverIslandCallback.h
//  BulletSoftBody
//
//  Created by Xuchen Han on 12/16/19.
//

#ifndef DeformableBodyInplaceSolverIslandCallback_h
#define DeformableBodyInplaceSolverIslandCallback_h

struct DeformableBodyInplaceSolverIslandCallback : public MultiBodyInplaceSolverIslandCallback
{
	btDeformableMultiBodyConstraintSolver* m_deformableSolver;

	DeformableBodyInplaceSolverIslandCallback(btDeformableMultiBodyConstraintSolver* solver,
											  btDispatcher* dispatcher)
		: MultiBodyInplaceSolverIslandCallback(solver, dispatcher), m_deformableSolver(solver)
	{
	}

	virtual void processConstraints(int islandId = -1)
	{
		btCollisionObject** bodies = m_bodies.size() ? &m_bodies[0] : 0;
		btCollisionObject** softBodies = m_softBodies.size() ? &m_softBodies[0] : 0;
		btPersistentManifold** manifold = m_manifolds.size() ? &m_manifolds[0] : 0;
		btTypedConstraint** constraints = m_constraints.size() ? &m_constraints[0] : 0;
		btMultiBodyConstraint** multiBodyConstraints = m_multiBodyConstraints.size() ? &m_multiBodyConstraints[0] : 0;

		//printf("mb contacts = %d, mb constraints = %d\n", mbContacts, m_multiBodyConstraints.size());

		m_deformableSolver->solveDeformableBodyGroup(bodies, m_bodies.size(), softBodies, m_softBodies.size(), manifold, m_manifolds.size(), constraints, m_constraints.size(), multiBodyConstraints, m_multiBodyConstraints.size(), *m_solverInfo, m_debugDrawer, m_dispatcher);
		if (m_bodies.size() && (m_solverInfo->m_reportSolverAnalytics & 1))
		{
			m_deformableSolver->m_analyticsData.m_islandId = islandId;
			m_islandAnalyticsData.push_back(m_solver->m_analyticsData);
		}
		m_bodies.resize(0);
		m_softBodies.resize(0);
		m_manifolds.resize(0);
		m_constraints.resize(0);
		m_multiBodyConstraints.resize(0);
	}
};

#endif /* DeformableBodyInplaceSolverIslandCallback_h */
// poly34.cpp : solution of cubic and quartic equation
// (c) Khashin S.I. http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
// khash2 (at) gmail.com
// Thanks to Alexandr Rakhmanin <rakhmanin (at) gmail.com>
// public domain
//
#include <math.h>

#include "poly34.h"  // solution of cubic and quartic equation
#define TwoPi 6.28318530717958648
const btScalar eps = SIMD_EPSILON;

//=============================================================================
// _root3, root3 from http://prografix.narod.ru
//=============================================================================
static SIMD_FORCE_INLINE btScalar _root3(btScalar x)
{
	btScalar s = 1.;
	while (x < 1.)
	{
		x *= 8.;
		s *= 0.5;
	}
	while (x > 8.)
	{
		x *= 0.125;
		s *= 2.;
	}
	btScalar r = 1.5;
	r -= 1. / 3. * (r - x / (r * r));
	r -= 1. / 3. * (r - x / (r * r));
	r -= 1. / 3. * (r - x / (r * r));
	r -= 1. / 3. * (r - x / (r * r));
	r -= 1. / 3. * (r - x / (r * r));
	r -= 1. / 3. * (r - x / (r * r));
	return r * s;
}

btScalar SIMD_FORCE_INLINE root3(btScalar x)
{
	if (x > 0)
		return _root3(x);
	else if (x < 0)
		return -_root3(-x);
	else
		return 0.;
}

// x - array of size 2
// return 2: 2 real roots x[0], x[1]
// return 0: pair of complex roots: x[0]i*x[1]
int SolveP2(btScalar* x, btScalar a, btScalar b)
{  // solve equation x^2 + a*x + b = 0
	btScalar D = 0.25 * a * a - b;
	if (D >= 0)
	{
		D = sqrt(D);
		x[0] = -0.5 * a + D;
		x[1] = -0.5 * a - D;
		return 2;
	}
	x[0] = -0.5 * a;
	x[1] = sqrt(-D);
	return 0;
}
//---------------------------------------------------------------------------
// x - array of size 3
// In case 3 real roots: => x[0], x[1], x[2], return 3
//         2 real roots: x[0], x[1],          return 2
//         1 real root : x[0], x[1]  i*x[2], return 1
int SolveP3(btScalar* x, btScalar a, btScalar b, btScalar c)
{  // solve cubic equation x^3 + a*x^2 + b*x + c = 0
	btScalar a2 = a * a;
	btScalar q = (a2 - 3 * b) / 9;
	if (q < 0)
		q = eps;
	btScalar r = (a * (2 * a2 - 9 * b) + 27 * c) / 54;
	// equation x^3 + q*x + r = 0
	btScalar r2 = r * r;
	btScalar q3 = q * q * q;
	btScalar A, B;
	if (r2 <= (q3 + eps))
	{  //<<-- FIXED!
		btScalar t = r / sqrt(q3);
		if (t < -1)
			t = -1;
		if (t > 1)
			t = 1;
		t = acos(t);
		a /= 3;
		q = -2 * sqrt(q);
		x[0] = q * cos(t / 3) - a;
		x[1] = q * cos((t + TwoPi) / 3) - a;
		x[2] = q * cos((t - TwoPi) / 3) - a;
		return (3);
	}
	else
	{
		//A =-pow(fabs(r)+sqrt(r2-q3),1./3);
		A = -root3(fabs(r) + sqrt(r2 - q3));
		if (r < 0)
			A = -A;
		B = (A == 0 ? 0 : q / A);

		a /= 3;
		x[0] = (A + B) - a;
		x[1] = -0.5 * (A + B) - a;
		x[2] = 0.5 * sqrt(3.) * (A - B);
		if (fabs(x[2]) < eps)
		{
			x[2] = x[1];
			return (2);
		}
		return (1);
	}
}  // SolveP3(btScalar *x,btScalar a,btScalar b,btScalar c) {
//---------------------------------------------------------------------------
// a>=0!
void CSqrt(btScalar x, btScalar y, btScalar& a, btScalar& b)  // returns:  a+i*s = sqrt(x+i*y)
{
	btScalar r = sqrt(x * x + y * y);
	if (y == 0)
	{
		r = sqrt(r);
		if (x >= 0)
		{
			a = r;
			b = 0;
		}
		else
		{
			a = 0;
			b = r;
		}
	}
	else
	{  // y != 0
		a = sqrt(0.5 * (x + r));
		b = 0.5 * y / a;
	}
}
//---------------------------------------------------------------------------
int SolveP4Bi(btScalar* x, btScalar b, btScalar d)  // solve equation x^4 + b*x^2 + d = 0
{
	btScalar D = b * b - 4 * d;
	if (D >= 0)
	{
		btScalar sD = sqrt(D);
		btScalar x1 = (-b + sD) / 2;
		btScalar x2 = (-b - sD) / 2;  // x2 <= x1
		if (x2 >= 0)                  // 0 <= x2 <= x1, 4 real roots
		{
			btScalar sx1 = sqrt(x1);
			btScalar sx2 = sqrt(x2);
			x[0] = -sx1;
			x[1] = sx1;
			x[2] = -sx2;
			x[3] = sx2;
			return 4;
		}
		if (x1 < 0)  // x2 <= x1 < 0, two pair of imaginary roots
		{
			btScalar sx1 = sqrt(-x1);
			btScalar sx2 = sqrt(-x2);
			x[0] = 0;
			x[1] = sx1;
			x[2] = 0;
			x[3] = sx2;
			return 0;
		}
		// now x2 < 0 <= x1 , two real roots and one pair of imginary root
		btScalar sx1 = sqrt(x1);
		btScalar sx2 = sqrt(-x2);
		x[0] = -sx1;
		x[1] = sx1;
		x[2] = 0;
		x[3] = sx2;
		return 2;
	}
	else
	{  // if( D < 0 ), two pair of compex roots
		btScalar sD2 = 0.5 * sqrt(-D);
		CSqrt(-0.5 * b, sD2, x[0], x[1]);
		CSqrt(-0.5 * b, -sD2, x[2], x[3]);
		return 0;
	}  // if( D>=0 )
}  // SolveP4Bi(btScalar *x, btScalar b, btScalar d)    // solve equation x^4 + b*x^2 d
//---------------------------------------------------------------------------
#define SWAP(a, b) \
	{              \
		t = b;     \
		b = a;     \
		a = t;     \
	}
static void dblSort3(btScalar& a, btScalar& b, btScalar& c)  // make: a <= b <= c
{
	btScalar t;
	if (a > b)
		SWAP(a, b);  // now a<=b
	if (c < b)
	{
		SWAP(b, c);  // now a<=b, b<=c
		if (a > b)
			SWAP(a, b);  // now a<=b
	}
}
//---------------------------------------------------------------------------
int SolveP4De(btScalar* x, btScalar b, btScalar c, btScalar d)  // solve equation x^4 + b*x^2 + c*x + d
{
	//if( c==0 ) return SolveP4Bi(x,b,d); // After that, c!=0
	if (fabs(c) < 1e-14 * (fabs(b) + fabs(d)))
		return SolveP4Bi(x, b, d);  // After that, c!=0

	int res3 = SolveP3(x, 2 * b, b * b - 4 * d, -c * c);  // solve resolvent
	// by Viet theorem:  x1*x2*x3=-c*c not equals to 0, so x1!=0, x2!=0, x3!=0
	if (res3 > 1)  // 3 real roots,
	{
		dblSort3(x[0], x[1], x[2]);  // sort roots to x[0] <= x[1] <= x[2]
		// Note: x[0]*x[1]*x[2]= c*c > 0
		if (x[0] > 0)  // all roots are positive
		{
			btScalar sz1 = sqrt(x[0]);
			btScalar sz2 = sqrt(x[1]);
			btScalar sz3 = sqrt(x[2]);
			// Note: sz1*sz2*sz3= -c (and not equal to 0)
			if (c > 0)
			{
				x[0] = (-sz1 - sz2 - sz3) / 2;
				x[1] = (-sz1 + sz2 + sz3) / 2;
				x[2] = (+sz1 - sz2 + sz3) / 2;
				x[3] = (+sz1 + sz2 - sz3) / 2;
				return 4;
			}
			// now: c<0
			x[0] = (-sz1 - sz2 + sz3) / 2;
			x[1] = (-sz1 + sz2 - sz3) / 2;
			x[2] = (+sz1 - sz2 - sz3) / 2;
			x[3] = (+sz1 + sz2 + sz3) / 2;
			return 4;
		}  // if( x[0] > 0) // all roots are positive
		// now x[0] <= x[1] < 0, x[2] > 0
		// two pair of comlex roots
		btScalar sz1 = sqrt(-x[0]);
		btScalar sz2 = sqrt(-x[1]);
		btScalar sz3 = sqrt(x[2]);

		if (c > 0)  // sign = -1
		{
			x[0] = -sz3 / 2;
			x[1] = (sz1 - sz2) / 2;  // x[0]i*x[1]
			x[2] = sz3 / 2;
			x[3] = (-sz1 - sz2) / 2;  // x[2]i*x[3]
			return 0;
		}
		// now: c<0 , sign = +1
		x[0] = sz3 / 2;
		x[1] = (-sz1 + sz2) / 2;
		x[2] = -sz3 / 2;
		x[3] = (sz1 + sz2) / 2;
		return 0;
	}  // if( res3>1 )    // 3 real roots,
	// now resoventa have 1 real and pair of compex roots
	// x[0] - real root, and x[0]>0,
	// x[1]i*x[2] - complex roots,
	// x[0] must be >=0. But one times x[0]=~ 1e-17, so:
	if (x[0] < 0)
		x[0] = 0;
	btScalar sz1 = sqrt(x[0]);
	btScalar szr, szi;
	CSqrt(x[1], x[2], szr, szi);  // (szr+i*szi)^2 = x[1]+i*x[2]
	if (c > 0)                    // sign = -1
	{
		x[0] = -sz1 / 2 - szr;  // 1st real root
		x[1] = -sz1 / 2 + szr;  // 2nd real root
		x[2] = sz1 / 2;
		x[3] = szi;
		return 2;
	}
	// now: c<0 , sign = +1
	x[0] = sz1 / 2 - szr;  // 1st real root
	x[1] = sz1 / 2 + szr;  // 2nd real root
	x[2] = -sz1 / 2;
	x[3] = szi;
	return 2;
}  // SolveP4De(btScalar *x, btScalar b, btScalar c, btScalar d)    // solve equation x^4 + b*x^2 + c*x + d
//-----------------------------------------------------------------------------
btScalar N4Step(btScalar x, btScalar a, btScalar b, btScalar c, btScalar d)  // one Newton step for x^4 + a*x^3 + b*x^2 + c*x + d
{
	btScalar fxs = ((4 * x + 3 * a) * x + 2 * b) * x + c;  // f'(x)
	if (fxs == 0)
		return x;                                       //return 1e99; <<-- FIXED!
	btScalar fx = (((x + a) * x + b) * x + c) * x + d;  // f(x)
	return x - fx / fxs;
}
//-----------------------------------------------------------------------------
// x - array of size 4
// return 4: 4 real roots x[0], x[1], x[2], x[3], possible multiple roots
// return 2: 2 real roots x[0], x[1] and complex x[2]i*x[3],
// return 0: two pair of complex roots: x[0]i*x[1],  x[2]i*x[3],
int SolveP4(btScalar* x, btScalar a, btScalar b, btScalar c, btScalar d)
{  // solve equation x^4 + a*x^3 + b*x^2 + c*x + d by Dekart-Euler method
	// move to a=0:
	btScalar d1 = d + 0.25 * a * (0.25 * b * a - 3. / 64 * a * a * a - c);
	btScalar c1 = c + 0.5 * a * (0.25 * a * a - b);
	btScalar b1 = b - 0.375 * a * a;
	int res = SolveP4De(x, b1, c1, d1);
	if (res == 4)
	{
		x[0] -= a / 4;
		x[1] -= a / 4;
		x[2] -= a / 4;
		x[3] -= a / 4;
	}
	else if (res == 2)
	{
		x[0] -= a / 4;
		x[1] -= a / 4;
		x[2] -= a / 4;
	}
	else
	{
		x[0] -= a / 4;
		x[2] -= a / 4;
	}
	// one Newton step for each real root:
	if (res > 0)
	{
		x[0] = N4Step(x[0], a, b, c, d);
		x[1] = N4Step(x[1], a, b, c, d);
	}
	if (res > 2)
	{
		x[2] = N4Step(x[2], a, b, c, d);
		x[3] = N4Step(x[3], a, b, c, d);
	}
	return res;
}
//-----------------------------------------------------------------------------
#define F5(t) (((((t + a) * t + b) * t + c) * t + d) * t + e)
//-----------------------------------------------------------------------------
btScalar SolveP5_1(btScalar a, btScalar b, btScalar c, btScalar d, btScalar e)  // return real root of x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
{
	int cnt;
	if (fabs(e) < eps)
		return 0;

	btScalar brd = fabs(a);  // brd - border of real roots
	if (fabs(b) > brd)
		brd = fabs(b);
	if (fabs(c) > brd)
		brd = fabs(c);
	if (fabs(d) > brd)
		brd = fabs(d);
	if (fabs(e) > brd)
		brd = fabs(e);
	brd++;  // brd - border of real roots

	btScalar x0, f0;       // less than root
	btScalar x1, f1;       // greater than root
	btScalar x2, f2, f2s;  // next values, f(x2), f'(x2)
	btScalar dx = 0;

	if (e < 0)
	{
		x0 = 0;
		x1 = brd;
		f0 = e;
		f1 = F5(x1);
		x2 = 0.01 * brd;
	}  // positive root
	else
	{
		x0 = -brd;
		x1 = 0;
		f0 = F5(x0);
		f1 = e;
		x2 = -0.01 * brd;
	}  // negative root

	if (fabs(f0) < eps)
		return x0;
	if (fabs(f1) < eps)
		return x1;

	// now x0<x1, f(x0)<0, f(x1)>0
	// Firstly 10 bisections
	for (cnt = 0; cnt < 10; cnt++)
	{
		x2 = (x0 + x1) / 2;  // next point
		//x2 = x0 - f0*(x1 - x0) / (f1 - f0);        // next point
		f2 = F5(x2);  // f(x2)
		if (fabs(f2) < eps)
			return x2;
		if (f2 > 0)
		{
			x1 = x2;
			f1 = f2;
		}
		else
		{
			x0 = x2;
			f0 = f2;
		}
	}

	// At each step:
	// x0<x1, f(x0)<0, f(x1)>0.
	// x2 - next value
	// we hope that x0 < x2 < x1, but not necessarily
	do
	{
		if (cnt++ > 50)
			break;
		if (x2 <= x0 || x2 >= x1)
			x2 = (x0 + x1) / 2;  // now  x0 < x2 < x1
		f2 = F5(x2);             // f(x2)
		if (fabs(f2) < eps)
			return x2;
		if (f2 > 0)
		{
			x1 = x2;
			f1 = f2;
		}
		else
		{
			x0 = x2;
			f0 = f2;
		}
		f2s = (((5 * x2 + 4 * a) * x2 + 3 * b) * x2 + 2 * c) * x2 + d;  // f'(x2)
		if (fabs(f2s) < eps)
		{
			x2 = 1e99;
			continue;
		}
		dx = f2 / f2s;
		x2 -= dx;
	} while (fabs(dx) > eps);
	return x2;
}  // SolveP5_1(btScalar a,btScalar b,btScalar c,btScalar d,btScalar e)    // return real root of x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
//-----------------------------------------------------------------------------
int SolveP5(btScalar* x, btScalar a, btScalar b, btScalar c, btScalar d, btScalar e)  // solve equation x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
{
	btScalar r = x[0] = SolveP5_1(a, b, c, d, e);
	btScalar a1 = a + r, b1 = b + r * a1, c1 = c + r * b1, d1 = d + r * c1;
	return 1 + SolveP4(x + 1, a1, b1, c1, d1);
}  // SolveP5(btScalar *x,btScalar a,btScalar b,btScalar c,btScalar d,btScalar e)    // solve equation x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
//-----------------------------------------------------------------------------
// poly34.h : solution of cubic and quartic equation
// (c) Khashin S.I. http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
// khash2 (at) gmail.com

#ifndef POLY_34
#define POLY_34
#include "LinearMath/btScalar.h"
// x - array of size 2
// return 2: 2 real roots x[0], x[1]
// return 0: pair of complex roots: x[0]i*x[1]
int SolveP2(btScalar* x, btScalar a, btScalar b);  // solve equation x^2 + a*x + b = 0

// x - array of size 3
// return 3: 3 real roots x[0], x[1], x[2]
// return 1: 1 real root x[0] and pair of complex roots: x[1]i*x[2]
int SolveP3(btScalar* x, btScalar a, btScalar b, btScalar c);  // solve cubic equation x^3 + a*x^2 + b*x + c = 0

// x - array of size 4
// return 4: 4 real roots x[0], x[1], x[2], x[3], possible multiple roots
// return 2: 2 real roots x[0], x[1] and complex x[2]i*x[3],
// return 0: two pair of complex roots: x[0]i*x[1],  x[2]i*x[3],
int SolveP4(btScalar* x, btScalar a, btScalar b, btScalar c, btScalar d);  // solve equation x^4 + a*x^3 + b*x^2 + c*x + d = 0 by Dekart-Euler method

// x - array of size 5
// return 5: 5 real roots x[0], x[1], x[2], x[3], x[4], possible multiple roots
// return 3: 3 real roots x[0], x[1], x[2] and complex x[3]i*x[4],
// return 1: 1 real root x[0] and two pair of complex roots: x[1]i*x[2],  x[3]i*x[4],
int SolveP5(btScalar* x, btScalar a, btScalar b, btScalar c, btScalar d, btScalar e);  // solve equation x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0

//-----------------------------------------------------------------------------
// And some additional functions for internal use.
// Your may remove this definitions from here
int SolveP4Bi(btScalar* x, btScalar b, btScalar d);                              // solve equation x^4 + b*x^2 + d = 0
int SolveP4De(btScalar* x, btScalar b, btScalar c, btScalar d);                  // solve equation x^4 + b*x^2 + c*x + d = 0
void CSqrt(btScalar x, btScalar y, btScalar& a, btScalar& b);                    // returns as a+i*s,  sqrt(x+i*y)
btScalar N4Step(btScalar x, btScalar a, btScalar b, btScalar c, btScalar d);     // one Newton step for x^4 + a*x^3 + b*x^2 + c*x + d
btScalar SolveP5_1(btScalar a, btScalar b, btScalar c, btScalar d, btScalar e);  // return real root of x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
#endif
