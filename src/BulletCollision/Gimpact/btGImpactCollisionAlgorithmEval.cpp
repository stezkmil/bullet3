#include "btGImpactCollisionAlgorithmEval.h"
#include "btGImpactShape.h"

bool btGImpactPairEval::EvalPair(const GIM_PAIR& pair,
								 btGimpactVsGimpactGroupedParams& grpParams, bool findOnlyFirstPenetratingPair,
						 ThreadLocalGImpactResult* perThreadIntermediateResults,
						 std::list<btGImpactIntermediateResult>* intermediateResults)
{
	btPrimitiveTriangle ptri0;
	btPrimitiveTriangle ptri1;
	GIM_TRIANGLE_CONTACT contact_data;

	grpParams.triface0 = pair.m_index1;
	grpParams.triface1 = pair.m_index2;

	grpParams.shape0->getPrimitiveTriangle(grpParams.triface0, ptri0);
	grpParams.shape1->getPrimitiveTriangle(grpParams.triface1, ptri1);

	btPrimitiveTriangle ptri0Backup = ptri0;
	btPrimitiveTriangle ptri1Backup = ptri1;

	ptri0.applyTransform(grpParams.orgtrans0);
	ptri1.applyTransform(grpParams.orgtrans1);

	if (ptri0.validity_test() && ptri1.validity_test())
	{
		//build planes
		ptri0.buildTriPlane();
		ptri1.buildTriPlane();

		if (ptri0.overlap_test(ptri1))
		{
			if (ptri0.find_triangle_collision_alt_method_outer(ptri1, contact_data, gMarginZoneRecoveryStrengthFactor, grpParams.lastSafeTrans0,
															   grpParams.lastSafeTrans1, ptri0Backup, ptri1Backup, grpParams.doUnstuck))
			{
				if (contact_data.m_point_count >= 1)
				{
					bool insert = true;
					if (findOnlyFirstPenetratingPair && contact_data.m_penetration_depth >= 0.0) // Ever since there is marginEpsilon in btPrimitiveTriangle::find_triangle_collision_alt_method_outer, the equal sign in the ">= 0.0" comparison has to be there, because 0.0 is now a very common non-penetrating result even in the margin zone
						insert = false;
					if (insert)
					{
						if (perThreadIntermediateResults)
							perThreadIntermediateResults->local().push_back({contact_data.m_points[0], contact_data.m_separating_normal, -contact_data.m_penetration_depth});
						if (intermediateResults)
							intermediateResults->push_back({contact_data.m_points[0], contact_data.m_separating_normal, -contact_data.m_penetration_depth});
					}
					return contact_data.m_penetration_depth < 0.0;
				}
			}
		}
	}
	return false;
}