#include "btGImpactCollisionAlgorithmEval.h"
#include "btGImpactShape.h"

bool btGImpactPairEval::EvalPair(const GIM_PAIR& pair,
								 btGimpactVsGimpactGroupedParams& grpParams, btFindOnlyFirstPairEnum findOnlyFirstTriPair,
								 bool isSelfCollision,
								 std::vector<int> node0Prev, std::vector<int> node1Prev,
								 ThreadLocalGImpactResult* perThreadIntermediateResults,
								 std::vector<btGImpactIntermediateResult>* intermediateResults)
{
	btPrimitiveTriangle ptri0;
	btPrimitiveTriangle ptri1;
	GIM_TRIANGLE_CONTACT contact_data;

	//fprintf(stderr, "pair.m_index1 %d pair.m_index2 %d\n", pair.m_index1, pair.m_index2);
	grpParams.shape0->getPrimitiveTriangle(pair.m_index1, ptri0, false);
	grpParams.shape1->getPrimitiveTriangle(pair.m_index2, ptri1, false);

	btPrimitiveTriangle ptri0Backup;
	btPrimitiveTriangle ptri1Backup;
	if (!grpParams.shape0->getPrimitiveTriangleSafe(pair.m_index1, ptri0Backup))
		ptri0Backup = ptri0;
	if (!grpParams.shape1->getPrimitiveTriangleSafe(pair.m_index2, ptri1Backup))
		ptri1Backup = ptri1;

	// Not sure whether softs (hasSafeX == true) should be subjected to this applyTransform. TODO verify.
	ptri0.applyTransform(grpParams.orgtrans0);
	ptri1.applyTransform(grpParams.orgtrans1);

	if (ptri0.validity_test() && ptri1.validity_test())
	{
		if (isSelfCollision)
		{
			// This discards triangle pairs which are too close to each other for self collisions. Otherwise all triangles
			// would collide with their topological neighbours.

			// This is not exactly fast. It would be better to have the triangle distances (on the original undeformed mesh) cached somehow,
			// but I currently do not know how to do it efficiently. This also does not work right for meshes which are
			// already tightly coiled in their original undeformed shapes (such triangles will always pass through), but I
			// have not met such shapes in practice yet.

			// Before, this was done by discarding immediately neighbouring triangles, but that was not enough. I still got too
			// many close-neighbour-collisions which seriously destabilised the simulation.

			// I have also tried to base this on bvh ancestry. Triangles with common bvh parents were assumed to lie close to each other.
			// Triangles with common bvh grand-parents were assumed to lie further from each other and would be tested. This approach also
			// was not good enough. There were still plenty of triangles which hand common only distant ancestry, but still lied geometriaclly
			// close.

			btPrimitiveTriangle ptri_orig0;
			btPrimitiveTriangle ptri_orig1;

			//fprintf(stderr, "pair.m_index1 %d pair.m_index2 %d\n", pair.m_index1, pair.m_index2);
			grpParams.shape0->getPrimitiveTriangle(pair.m_index1, ptri_orig0, true);
			grpParams.shape1->getPrimitiveTriangle(pair.m_index2, ptri_orig1, true);

			ptri_orig0.buildTriPlane();
			ptri_orig1.buildTriPlane();

			btScalar dist_sq_out;
			btVector3 a_closest_out, b_closest_out;
			auto ret = ptri_orig0.triangle_triangle_distance(ptri_orig1, dist_sq_out, a_closest_out, b_closest_out);
			auto dist = sqrtf(dist_sq_out);

			const btScalar marginInflationFactor = 5.0;
			btScalar margin = std::max((ptri_orig0.m_margin + ptri_orig1.m_margin) * marginInflationFactor, 1.0);

			if (ret && dist_sq_out != 0.0 && dist < margin)
				return false;
			else if (ret && dist_sq_out == 0.0)
				return false;

			// It would seem natural just to compare the triangle indices to check if the triangles are adjacent, but we are doing mesh subdivisions which
			// do not do create shared indices (to make the subdivisions less costly), so all that we can do now is to compare the vertex coordinates. This
			// is not completely bulletproof, but probably the best I can do now.
			//constexpr btScalar vertIsEqualEpsilon = 0.001;
			//for (auto i = 0; i < 3; ++i)
			//	for (auto j = 0; j < 3; ++j)
			//		if ((ptri0.m_vertices[i] - ptri1.m_vertices[j]).length2() < vertIsEqualEpsilon)
			//			return false;  // It is an adjacent triangle. We ignore adjacent triangles in self collisions now.
		}

		//build planes
		ptri0.buildTriPlane();
		ptri1.buildTriPlane();

		if (ptri0.overlap_test(ptri1))
		{
			if (ptri0.find_triangle_collision_alt_method_outer(ptri1, contact_data, gMarginZoneRecoveryStrengthFactor, grpParams.lastSafeTrans0,
															   grpParams.lastSafeTrans1, ptri0Backup, ptri1Backup, grpParams.doUnstuck, isSelfCollision))
			{
				if (contact_data.m_point_count >= 1)
				{
					bool insert = true;
					if (findOnlyFirstTriPair == btFindOnlyFirstPairEnum::PENETRATING && contact_data.m_penetration_depth >= 0.0)  // Ever since there is marginEpsilon in btPrimitiveTriangle::find_triangle_collision_alt_method_outer, the equal sign in the ">= 0.0" comparison has to be there, because 0.0 is now a very common non-penetrating result even in the margin zone
						insert = false;
					if (insert)
					{
						if (perThreadIntermediateResults)
							perThreadIntermediateResults->local().push_back({contact_data.m_points[0], contact_data.m_separating_normal, -contact_data.m_penetration_depth, contact_data.m_unmodified_depth, pair.m_index1, pair.m_index2});
						if (intermediateResults)
							intermediateResults->push_back({contact_data.m_points[0], contact_data.m_separating_normal, -contact_data.m_penetration_depth, contact_data.m_unmodified_depth, pair.m_index1, pair.m_index2});
					}
					if (findOnlyFirstTriPair == btFindOnlyFirstPairEnum::PENETRATING)
						return contact_data.m_penetration_depth < 0.0;
					else
						return true;
				}
			}
		}
	}
	return false;
}
