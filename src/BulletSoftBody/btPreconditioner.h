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
#include <cmath>

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
	btAlignedObjectArray<btMatrix3x3> m_invBlocks;
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
		if (m_implicit)
		{
			buildImplicitBlocks();
		}
		else
		{
			buildDiagonalA(m_inv_A);
			for (int i = 0; i < m_inv_A.size(); ++i)
			{
				//            printf("A[%d] = %f, %f, %f \n", i, m_inv_A[i][0], m_inv_A[i][1], m_inv_A[i][2]);
				for (int d = 0; d < 3; ++d)
				{
					m_inv_A[i][d] = (m_inv_A[i][d] == 0) ? 0.0 : 1.0 / m_inv_A[i][d];
				}
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



	btVector3 applyInverseNodeBlock(int index, const btVector3& x) const
	{
		return m_implicit ? m_invBlocks[index] * x : x * m_inv_A[index];
	}

	// Scale before Cholesky to avoid a determinant-based inverse on ill-scaled blocks.
	static bool invertPositiveBlock(const btMatrix3x3& block, btMatrix3x3& inverse, bool& regularized)
	{
		btScalar scale = 0;
		for (int r = 0; r < 3; ++r)
			for (int c = 0; c < 3; ++c)
			{
				if (!std::isfinite((double)block[r][c]))
					return false;
				scale = btMax(scale, btFabs(block[r][c]));
			}
		if (!(scale > 0))
			return false;
		const btMatrix3x3 normalized = (block * (btScalar(1) / scale) + block.transpose() * (btScalar(1) / scale)) * btScalar(0.5);
		btScalar shift = 0;
		for (int attempt = 0; attempt < 6; ++attempt)
		{
			btScalar l[3][3] = {};
			bool positive = true;
			for (int r = 0; r < 3 && positive; ++r)
				for (int c = 0; c <= r; ++c)
				{
					btScalar value = normalized[r][c] + (r == c ? shift : btScalar(0));
					for (int k = 0; k < c; ++k)
						value -= l[r][k] * l[c][k];
					if (r == c)
					{
						if (!(value > btScalar(0))) { positive = false; break; }
						l[r][c] = btSqrt(value);
					}
					else
						l[r][c] = value / l[c][c];
				}
			if (positive)
			{
				for (int column = 0; column < 3; ++column)
				{
					btScalar y[3] = {}, x[3] = {};
					for (int r = 0; r < 3; ++r)
					{
						btScalar value = r == column ? btScalar(1) : btScalar(0);
						for (int k = 0; k < r; ++k) value -= l[r][k] * y[k];
						y[r] = value / l[r][r];
					}
					for (int r = 2; r >= 0; --r)
					{
						btScalar value = y[r];
						for (int k = r + 1; k < 3; ++k) value -= l[k][r] * x[k];
						x[r] = value / l[r][r];
						inverse[r][column] = x[r] / scale;
						positive = positive && std::isfinite((double)inverse[r][column]);
					}
				}
				if (positive)
				{
					inverse = (inverse + inverse.transpose()) * btScalar(0.5);
					regularized = attempt > 0;
					return true;
				}
			}
			shift = attempt == 0 ? btMax(btScalar(1e-8), btScalar(64) * SIMD_EPSILON) : shift * btScalar(10);
		}
		return false;
	}

	void buildImplicitBlocks()
	{
		BT_PROFILE("buildImplicitBlocks");
		m_invBlocks.resize(m_inv_A.size());
		const btMatrix3x3 identity = btMatrix3x3::getIdentity();
		int index = 0;
		for (int b = 0; b < m_softBodies.size(); ++b)
			for (int n = 0; n < m_softBodies[b]->m_nodes.size(); ++n, ++index)
			{
				const btSoftBody::Node& node = m_softBodies[b]->m_nodes[n];
				m_invBlocks[index] = identity * ((node.m_frozen <= 0 && node.m_im > 0) ? btScalar(1) / node.m_im : btScalar(0));
			}
		TVStack diagonal;
		for (int f = 0; f < m_lf.size(); ++f)
		{
			if (m_lf[f]->addImplicitForceDifferentialBlocks(m_dt, m_invBlocks))
				continue;
			// Forces without a block implementation retain their existing approximation.
			diagonal.resize(m_inv_A.size());
			for (int n = 0; n < diagonal.size(); ++n) diagonal[n].setZero();
			m_lf[f]->buildDampingForceDifferentialDiagonal(-m_dt, diagonal);
			for (int n = 0; n < diagonal.size(); ++n)
				for (int d = 0; d < 3; ++d) m_invBlocks[n][d][d] += diagonal[n][d];
		}
		index = 0;
		for (int b = 0; b < m_softBodies.size(); ++b)
			for (int n = 0; n < m_softBodies[b]->m_nodes.size(); ++n, ++index)
			{
				const btSoftBody::Node& node = m_softBodies[b]->m_nodes[n];
				btMatrix3x3 inverse;
				bool regularized = false;
				if (node.m_frozen > 0 || node.m_im <= 0)
					inverse = identity * btScalar(0);
				else if (!invertPositiveBlock(m_invBlocks[index], inverse, regularized))
				{
					inverse = identity * node.m_im;
				}
				m_invBlocks[index] = inverse;
				m_inv_A[index] = btVector3(inverse[0][0], inverse[1][1], inverse[2][2]);
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
				diagA[counter] = (node.m_frozen > 0) ? btVector3(0, 0, 0) : btVector3(1.0 / node.m_im, 1.0 / node.m_im, 1.0 / node.m_im);
				++counter;
			}
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
					if (m_implicit)
					{
						t[j] += lm.m_dirs[j].dot(applyInverseNodeBlock(lm.m_indices[i], lm.m_dirs[j])) * lm.m_weights[i] * lm.m_weights[i];
					}
					else
					{
						for (int d = 0; d < 3; ++d)
						{
							t[j] += inv_A[lm.m_indices[i]][d] * lm.m_dirs[j][d] * lm.m_dirs[j][d] * lm.m_weights[i] * lm.m_weights[i];
						}
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
			b[i] = applyInverseNodeBlock(i, x[i]);
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
			b[i] = applyInverseNodeBlock(i, x[i]);
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
			b[i] = applyInverseNodeBlock(i, x[i] - b[i]);
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
			b[i] += applyInverseNodeBlock(i, t[i]);
		}

		for (int i = 0; i < m_inv_S.size(); ++i)
		{
			b[i + offset] -= x[i + offset] * m_inv_S[i];
		}
	}
#endif
};

#endif /* BT_PRECONDITIONER_H */
