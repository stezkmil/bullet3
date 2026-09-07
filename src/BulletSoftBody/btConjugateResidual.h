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
	btScalar m_initialResidual;
	btScalar m_finalResidual;
	btScalar m_targetResidual;
	btScalar m_relativeTolerance;
	bool m_stagnated;

public:
	btConjugateResidual(const int max_it_in)
		: Base(max_it_in, 1e-8),
		  m_initialResidual(0),
		  m_finalResidual(0),
		  m_targetResidual(0),
		  m_relativeTolerance(1e-4),
		  m_stagnated(false)
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
		m_initialResidual = residual_norm;
		m_finalResidual = residual_norm;
		m_targetResidual = btMax(Base::m_tolerance, m_relativeTolerance * residual_norm);
		m_stagnated = false;
		best_x = x;
		best_r = residual_norm;
		if (residual_norm <= m_targetResidual)
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
		const int stagnationWindow = 25;
		const btScalar minimumWindowReduction = btScalar(1e-3);
		btScalar windowInitialBest = best_r;
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
				m_finalResidual = norm_r;
				if (norm_r <= m_targetResidual)
				{
					return k;
				}
			}
			// Returning the best iterate after a flat residual window avoids spending
			// the remainder of the iteration budget on roundoff-level progress.  The
			// test uses the best residual, so temporary CR oscillations do not trigger it.
			if (k % stagnationWindow == 0)
			{
				const btScalar requiredBest = windowInitialBest * (btScalar(1) - minimumWindowReduction);
				if (!(best_r < requiredBest))
				{
					x = best_x;
					m_finalResidual = best_r;
					m_stagnated = true;
					if (verbose)
					{
						std::cout << "ConjugateResidual stagnated at iteration " << k << ", residual = " << best_r << std::endl;
					}
					return k;
				}
				windowInitialBest = best_r;
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
		m_finalResidual = best_r;
		return Base::m_maxIterations;
	}

	btScalar getInitialResidual() const
	{
		return m_initialResidual;
	}

	btScalar getFinalResidual() const
	{
		return m_finalResidual;
	}

	btScalar getTargetResidual() const
	{
		return m_targetResidual;
	}

	btScalar getResidualRatio() const
	{
		return m_initialResidual > btScalar(0) ? m_finalResidual / m_initialResidual : btScalar(0);
	}

	bool getStagnated() const
	{
		return m_stagnated;
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
