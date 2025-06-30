/*! \file btGImpactShape.h
\author Francisco Leon Najera
*/
/*
This source file is part of GIMPACT Library.

For the latest info, see http://gimpact.sourceforge.net/

Copyright (c) 2007 Francisco Leon Najera. C.C. 80087371.
email: projectileman@yahoo.com


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

#ifndef BT_GIMPACT_BVH_CONCAVE_COLLISION_ALGORITHM_EVAL_H
#define BT_GIMPACT_BVH_CONCAVE_COLLISION_ALGORITHM_EVAL_H

#include "LinearMath/btVector3.h"
#include "LinearMath/btScalar.h"
#include "LinearMath/btTransform.h"
#include "gim_pair.h"

#include <tbb/tbb.h>

enum class btFindOnlyFirstPairEnum : uint8_t
{
	DISABLED = 0,
	PENETRATING = 1,
	TOUCHING = 2
};

struct btGImpactIntermediateResult
{
	btVector3 point;
	btVector3 normal;
	btScalar depth;
	btScalar unmodified_depth;
	int index0, index1;
};

// Experimentally switching to vectors, because they give better results (when destroying the vectors) when the number of contacts is huge (contact count collecting mode)
typedef tbb::enumerable_thread_specific<std::vector<btGImpactIntermediateResult>> ThreadLocalGImpactResult;

class btGImpactMeshShapePart;

struct btGimpactVsGimpactGroupedParams
{
	const btGImpactMeshShapePart* shape0;
	const btGImpactMeshShapePart* shape1;
	btTransform orgtrans0;
	btTransform orgtrans1;
	btTransform lastSafeTrans0;
	btTransform lastSafeTrans1;
	bool doUnstuck;
	int previouslyConsumedTime;
	btGimpactVsGimpactGroupedParams() {}
};

struct btGImpactPairEval
{
	static bool EvalPair(const GIM_PAIR& pair,
						 btGimpactVsGimpactGroupedParams& grpParams, btFindOnlyFirstPairEnum findOnlyFirstTriPair, bool isSelfCollision,
						 ThreadLocalGImpactResult* perThreadIntermediateResults,
						 std::vector<btGImpactIntermediateResult>* intermediateResults);
};

#endif  //BT_GIMPACT_BVH_CONCAVE_COLLISION_ALGORITHM_H
