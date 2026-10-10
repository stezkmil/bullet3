// SPDX-FileCopyrightText: Copyright (c) 2025 The Newton Developers
// SPDX-License-Identifier: Apache-2.0
// Experimental C++ adaptation of Newton 009158e62b86; triangle discovery is Bullet-specific.
#ifndef BT_DEFORMABLE_VBD_SOLVER_H
#define BT_DEFORMABLE_VBD_SOLVER_H
#include "btDeformableVbdKernels.h"
#include "btDeformableVbdSettings.h"
#include "btDeformableVbdCollision.h"
#include "btDeformableVbdParallel.h"
#include "BulletCollision/CollisionShapes/btTriangleShape.h"
#include "BulletCollision/Gimpact/btGImpactShape.h"
#include "BulletCollision/NarrowPhaseCollision/btGjkPairDetector.h"
#include "BulletCollision/NarrowPhaseCollision/btVoronoiSimplexSolver.h"
#include "BulletCollision/NarrowPhaseCollision/btGjkEpaPenetrationDepthSolver.h"
#include "BulletCollision/NarrowPhaseCollision/btPointCollector.h"
#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <vector>
#include <limits>

// Use the same dot, translation, and unit-conversion order as nativeTriangle.
static btVector3 btVbdNativeExtractPosition(const btTransform &transform, const btVector3 &vertex, btScalar units)
{
	const btScalar inverse = btScalar(1) / units;
	const auto &basis = transform.getBasis();
	const auto &origin = transform.getOrigin();
	btVector3 result;
	for (int axis = 0; axis < 3; ++axis)
	{
		const auto &row = basis[axis];
		const btScalar first = vertex.x() * row.x() + vertex.y() * row.y();
		const btScalar rotated = first + vertex.z() * row.z();
		result[axis] = (rotated + origin[axis]) * inverse;
	}
	result[3] = 0;
	return result;
}

// All quantities are SI. This backend deliberately does not invoke Bullet's coupled Newton solve.
class btDeformableVbdSolver
{
  public:
	using Settings = btDeformableVbdSettings;
	explicit btDeformableVbdSolver(btDeformableVbdCollisionCache *cache = nullptr)
		: ownedCollisionCache(cache ? nullptr : new btDeformableVbdCollisionCache()),
		  softMesh((cache ? cache : ownedCollisionCache.get())->soft), barrierMesh((cache ? cache : ownedCollisionCache.get())->barrier),
		  gpu((cache ? cache : ownedCollisionCache.get())->gpu),
		  guardPositions((cache ? cache : ownedCollisionCache.get())->guardPositions),
		  guardStamp((cache ? cache : ownedCollisionCache.get())->guardStamp),
		  guardReferenceStamp((cache ? cache : ownedCollisionCache.get())->guardReferenceStamp),
		  pairDistances((cache ? cache : ownedCollisionCache.get())->pairDistances),
		  pairSortScratch((cache ? cache : ownedCollisionCache.get())->pairSortScratch)
	{
	}
	struct Tet
	{
		std::array<int, 4> nodes;
		btMatrix3x3 inverseRest;
		btScalar volume, mu, lambda, damping;
	};
	struct Triangle
	{
		btVector3 x[3];
		btScalar friction;
		int owner = -1;
		int rigid = -1;
		btVector3 local[3];
	};
	struct Support
	{
		int node;
		btMatrix3x3 jacobian;
	};
	struct SurfaceVertex
	{
		std::vector<Support> support;
		btVector3 offset = btVector3(0, 0, 0);
		btVector3 position(const std::vector<btVector3> &positions) const
		{
			btVector3 p = offset;
			for (const auto &s : support)
				p += s.jacobian * positions[s.node];
			return p;
		}
	};
	struct Contact
	{
		SurfaceVertex vertex;
		btVector3 point, normal;
		btScalar friction;
		int rigid = -1, rigidA = -1;
		int convexA = -1, convexB = -1, triangleA = -1, triangleB = -1;
		btScalar referenceGap = 0;
		btVector3 localA = btVector3(0, 0, 0);
		btVector3 localPoint = btVector3(0, 0, 0);
	};
	struct NodeSpring
	{
		int node;
		btVector3 target;
		btScalar stiffness, damping, maxForce;
	};
	std::vector<NodeSpring> springs;
	static void springBlock(const NodeSpring &spring, const btVector3 &position, const btVector3 &previous, btScalar dt, btVector3 &force,
							btMatrix3x3 &hessian)
	{
		const btVector3 delta = position - spring.target;
		const btScalar length = delta.length();
		const btMatrix3x3 identity = btMatrix3x3::getIdentity();
		if (spring.stiffness * length <= spring.maxForce || length == 0)
		{
			force -= delta * spring.stiffness;
			hessian += identity * spring.stiffness;
		}
		else
		{
			const btVector3 direction = delta / length;
			force -= direction * spring.maxForce;
			hessian += (identity - btVbd::outer(direction, direction)) * (spring.maxForce / length);
		}
		// Freeze the damping axis for this step to keep its block symmetric positive semidefinite.
		const btVector3 oldDelta = previous - spring.target;
		const btMatrix3x3 damping = oldDelta.length2() > SIMD_EPSILON * SIMD_EPSILON
										? btVbd::outer(oldDelta.normalized(), oldDelta.normalized()) * spring.damping
										: identity * spring.damping;
		force -= damping * ((position - previous) / dt);
		hessian += damping * (1 / dt);
	}
	struct Rigid
	{
		btTransform pose = btTransform::getIdentity(), previous = btTransform::getIdentity(), reference = btTransform::getIdentity(),
					predicted = btTransform::getIdentity();
		btVector3 velocity = btVector3(0, 0, 0), angularVelocity = btVector3(0, 0, 0), force = btVector3(0, 0, 0),
				  torque = btVector3(0, 0, 0);
		btVector3 inertia = btVector3(1, 1, 1);
		btVector3 linearFactor = btVector3(1, 1, 1), angularFactor = btVector3(1, 1, 1);
		btScalar mass = 0, radius = 0, linearDamping = 0, angularDamping = 0;
	};
	std::vector<Rigid> rigids;
	struct JointRow
	{
		int a = -1, b = -1;
		btVector3 linearA = btVector3(0, 0, 0), angularA = btVector3(0, 0, 0), linearB = btVector3(0, 0, 0), angularB = btVector3(0, 0, 0);
		btScalar rhs = 0, cfm = 0, lower = -SIMD_INFINITY, upper = SIMD_INFINITY, impulse = 0;
	};
	std::vector<JointRow> joints;
	struct Attachment
	{
		int node, rigid;
		btVector3 localPoint;
		btScalar stiffness = 1e7, damping = 10;
	};
	std::vector<Attachment> attachments;
	Settings settings;
	int gpuGuardCalls = 0, gpuGuardFallbacks = 0;
	int gpuSurfaceCalls = 0, gpuSurfaceFallbacks = 0;
	std::string gpuError;
	std::vector<btVector3> x, velocity, external;
	std::vector<btVector3> recoveryPositions;
	std::vector<btScalar> mass, massDamping;
	std::vector<Tet> tets;
	std::vector<std::array<int, 3>> surface;
	std::vector<SurfaceVertex> surfaceVertices;
	bool mappedSurface = false;
	std::vector<Triangle> barriers;
	struct NativeBarrier
	{
		const btGImpactMeshShapePart *shape;
		btTransform transform; // Existing shape coordinates to Bullet world units.
		btScalar friction;
		int owner;
		std::vector<int> triangles;
	};
	std::vector<NativeBarrier> nativeBarriers;
	struct ConvexBarrier
	{
		const btConvexShape *shape;
		btTransform transform; // Bullet collision units, matching collisionUnitsPerMeter.
		btScalar friction;
		int owner = -1;
		int rigid = -1;
		btTransform localTransform = btTransform::getIdentity();
	};
	std::vector<ConvexBarrier> convexBarriers;
	std::vector<int> surfaceOwners;
	std::vector<btScalar> surfaceFrictions;
	std::vector<bool> selfContact;
	std::vector<std::vector<bool>> collisionAllowed;
	struct MovingPlane
	{
		int a, b;
		btVector3 normal;
	};
	std::vector<MovingPlane> movingPlanes;
	std::vector<Contact> contacts;
	struct Plane
	{
		std::array<int, 3> nodes;
		btVector3 point, normal;
		int rigid = -1, triangle = -1, convex = -1;
		btScalar referenceGap = 0;
	};
	std::vector<Plane> planes;
	int broadphasePairs = 0;
	int recoveredIntersections = 0;
	int colorCount = 0;
	btScalar minimumJ = 1;
	const char *error = nullptr;

	static btVector3 closestWeights(const btVector3 &p, const btVector3 *t)
	{
		const btVector3 ab = t[1] - t[0], ac = t[2] - t[0], ap = p - t[0];
		const btScalar d1 = ab.dot(ap), d2 = ac.dot(ap);
		if (d1 <= 0 && d2 <= 0)
			return btVector3(1, 0, 0);
		const btVector3 bp = p - t[1];
		const btScalar d3 = ab.dot(bp), d4 = ac.dot(bp);
		if (d3 >= 0 && d4 <= d3)
			return btVector3(0, 1, 0);
		const btScalar vc = d1 * d4 - d3 * d2;
		if (vc <= 0 && d1 >= 0 && d3 <= 0)
		{
			const btScalar v = d1 / (d1 - d3);
			return btVector3(1 - v, v, 0);
		}
		const btVector3 cp = p - t[2];
		const btScalar d5 = ab.dot(cp), d6 = ac.dot(cp);
		if (d6 >= 0 && d5 <= d6)
			return btVector3(0, 0, 1);
		const btScalar vb = d5 * d2 - d1 * d6;
		if (vb <= 0 && d2 >= 0 && d6 <= 0)
		{
			const btScalar w = d2 / (d2 - d6);
			return btVector3(1 - w, 0, w);
		}
		const btScalar va = d3 * d6 - d5 * d4;
		if (va <= 0 && d4 >= d3 && d5 >= d6)
		{
			const btScalar w = (d4 - d3) / (d4 - d3 + d5 - d6);
			return btVector3(0, 1 - w, w);
		}
		const btScalar sum = va + vb + vc;
		return btVector3(va / sum, vb / sum, vc / sum);
	}
	static btVector3 rotationVector(btQuaternion q)
	{
		q.normalize();
		if (q.w() < 0)
			q = btQuaternion(-q.x(), -q.y(), -q.z(), -q.w());
		const btScalar length = btSqrt(q.x() * q.x() + q.y() * q.y() + q.z() * q.z());
		return length > 1e-12 ? btVector3(q.x(), q.y(), q.z()) * (2 * btAtan2(length, q.w()) / length) : btVector3(0, 0, 0);
	}
	static btQuaternion rotationIncrement(const btVector3 &v)
	{
		const btScalar angle = v.length();
		return angle > 1e-12 ? btQuaternion(v / angle, angle) : btQuaternion::getIdentity();
	}
	static btMatrix3x3 pointRotationJacobian(const btVector3 &r)
	{
		return btMatrix3x3(0, r.z(), -r.y(), -r.z(), 0, r.x(), r.y(), -r.x(), 0);
	}
	void contactBlock(const Contact &c, btScalar dt, btVector3 &force, btMatrix3x3 &hessian) const
	{
		const btVector3 p = c.rigidA < 0 ? c.vertex.position(x) : rigids[c.rigidA].pose * c.localA;
		const btVector3 p0 = c.rigidA < 0 ? c.vertex.position(previous) : rigids[c.rigidA].previous * c.localA;
		const btVector3 target = c.rigid < 0 ? c.point : rigids[c.rigid].pose * c.localPoint;
		const btVector3 oldTarget = c.rigid < 0 ? c.point : rigids[c.rigid].previous * c.localPoint;
		btVbd::contactBlock(c.normal.dot(p - target), settings.radius, c.normal, (p - p0) - (target - oldTarget), settings.ke, settings.kd,
							c.friction, settings.frictionEpsilon, dt, force, hessian);
	}
	struct ScalarMappingData
	{
		std::vector<std::array<int, 2>> ranges;
		std::vector<btVbd::ScalarMappingSupport> supports;
		bool allScalar = true;
	};
	struct SurfaceTopology
	{
		std::shared_ptr<const std::vector<std::set<int>>> neighbors;
		btScalar amplification = 0;
		std::shared_ptr<const ScalarMappingData> scalarMappings;
	};
	SurfaceTopology surfaceTopology() const
	{
		return {baseNeighbors, mappingAmplification, scalarMappings};
	}

	void includeSurfaceBounds(btScalar scale, btVector3 &lo, btVector3 &hi) const
	{
		if (settings.workers <= 1 || surfaceVertices.size() < 4096)
		{
			for (int i = 0; i < int(surfaceVertices.size()); ++i)
			{
				const btVector3 point = surfacePosition(i, x) / scale;
				lo.setMin(point);
				hi.setMax(point);
			}
			return;
		}
		const int count = int(surfaceVertices.size()), blockSize = 256;
		std::vector<std::pair<btVector3, btVector3>> bounds((count + blockSize - 1) / blockSize, {lo, hi});
		btVbdParallelFor(
			int(bounds.size()), settings.workers,
			[&](int block)
			{
				auto &range = bounds[block];
				const int end = btMin(count, (block + 1) * blockSize);
				for (int i = block * blockSize; i < end; ++i)
				{
					const btVector3 point = surfacePosition(i, x) / scale;
					range.first.setMin(point);
					range.second.setMax(point);
				}
			},
			2, 1);
		for (const auto &range : bounds)
		{
			lo.setMin(range.first);
			hi.setMax(range.second);
		}
	}
	// Reuse requires an exact match of mapping, connectivity, transforms, and node indices.
	bool initialize(const SurfaceTopology *cachedTopology = nullptr)
	{
		error = nullptr;
		guardSurfaceVertices.clear();
		guardPlaneNodes.clear();
		if (!cachedTopology)
			gpu.invalidateMapping();
		gpuPrepared = cachedTopology && gpu.hasMapping();
		gpuReferenceValid = false;
		gpuReferenceStamp = 0;
		gpuGuardCalls = gpuGuardFallbacks = 0;
		gpuSurfaceCalls = gpuSurfaceFallbacks = 0;
		gpuError.clear();
		const int n = int(x.size());
		if (!n || tets.empty() || velocity.size() != x.size() || external.size() != x.size() || mass.size() != x.size() ||
			massDamping.size() != x.size())
		{
			error = "invalid_node_arrays";
			return false;
		}
		attachmentAdjacency.assign(n, {});
		for (int k = 0; k < int(attachments.size()); ++k)
		{
			const auto &a = attachments[k];
			if (a.node < 0 || a.node >= n || a.rigid < 0 || a.rigid >= int(rigids.size()) || !(a.stiffness > 0) || a.damping < 0)
			{
				error = "invalid_attachment";
				return false;
			}
			attachmentAdjacency[a.node].push_back(k);
		}
		springAdjacency.assign(n, {});
		for (int k = 0; k < int(springs.size()); ++k)
		{
			const auto &spring = springs[k];
			if (spring.node < 0 || spring.node >= n || !std::isfinite(double(spring.target.length2())) ||
				!std::isfinite(double(spring.stiffness)) || !std::isfinite(double(spring.damping)) ||
				!std::isfinite(double(spring.maxForce)) || spring.stiffness < 0 || spring.damping < 0 || spring.maxForce < 0)
			{
				error = "invalid_grab_spring";
				return false;
			}
			springAdjacency[spring.node].push_back(k);
		}
		adjacency.assign(n, {});
		colors.assign(n, -1);
		groups.clear();
		if (!mappedSurface)
		{
			surface.clear();
			surfaceVertices.clear();
		}
		std::map<std::array<int, 3>, int> faces;
		for (int k = 0; k < int(tets.size()); ++k)
		{
			const Tet &tet = tets[k];
			if (!(tet.volume > 0) || !(tet.mu > 0) || !(tet.mu + tet.lambda > 0) || !(tet.damping >= 0))
			{
				error = "invalid_rest_tet";
				return false;
			}
			for (int j = 0; j < 4; ++j)
			{
				if (tet.nodes[j] < 0 || tet.nodes[j] >= n)
				{
					error = "invalid_tet_node";
					return false;
				}
				adjacency[tet.nodes[j]].push_back(std::make_pair(k, j));
				if (!mappedSurface)
				{
					std::array<int, 3> face;
					int c = 0;
					for (int i = 0; i < 4; ++i)
						if (i != j)
							face[c++] = tet.nodes[i];
					std::sort(face.begin(), face.end());
					++faces[face];
				}
			}
		}
		if (!mappedSurface)
		{
			for (int i = 0; i < n; ++i)
			{
				SurfaceVertex v;
				v.support.push_back({i, btMatrix3x3::getIdentity()});
				surfaceVertices.push_back(v);
			}
			for (const auto &face : faces)
				if (face.second == 1)
					surface.push_back(face.first);
		}
		if (surface.empty() || surfaceVertices.empty())
		{
			error = "empty_collision_surface";
			return false;
		}
		if (cachedTopology && mappedSurface && cachedTopology->neighbors && cachedTopology->neighbors->size() == x.size() &&
			cachedTopology->scalarMappings && cachedTopology->scalarMappings->ranges.size() == surfaceVertices.size())
		{
			baseNeighbors = cachedTopology->neighbors;
			mappingAmplification = cachedTopology->amplification;
			scalarMappings = cachedTopology->scalarMappings;
		}
		else
		{
			mappingAmplification = 0;
			auto scalarData = std::make_shared<ScalarMappingData>();
			scalarData->ranges.reserve(surfaceVertices.size());
			for (const auto &v : surfaceVertices)
			{
				btScalar amplification = 0;
				bool scalar = true;
				if (v.support.empty() || !std::isfinite(double(v.offset.length2())))
				{
					error = "invalid_surface_mapping";
					return false;
				}
				for (const auto &a : v.support)
				{
					if (a.node < 0 || a.node >= n)
					{
						error = "invalid_surface_node";
						return false;
					}
					for (int d = 0; d < 3; ++d)
						if (!std::isfinite(double(a.jacobian[d].length2())))
						{
							error = "invalid_surface_weight";
							return false;
						}
				}
				for (const auto &a : v.support)
				{
					const auto &j = a.jacobian;
					scalar = scalar && j[0][1] == 0 && j[0][2] == 0 && j[1][0] == 0 && j[1][2] == 0 && j[2][0] == 0 && j[2][1] == 0 &&
							 j[0][0] == j[1][1] && j[0][0] == j[2][2];
					// A diagonal mapping's exact operator norm avoids unnecessary detailed guards.
					if (j[0][1] == 0 && j[0][2] == 0 && j[1][0] == 0 && j[1][2] == 0 && j[2][0] == 0 && j[2][1] == 0)
						amplification += btMax(btFabs(j[0][0]), btMax(btFabs(j[1][1]), btFabs(j[2][2])));
					else
						amplification += btSqrt(j[0].length2() + j[1].length2() + j[2].length2());
				}
				mappingAmplification = btMax(mappingAmplification, amplification);
				const int begin = int(scalarData->supports.size());
				if (scalar)
					for (const auto &support : v.support)
						scalarData->supports.push_back({support.node, support.jacobian[0][0]});
				scalarData->ranges.push_back({{scalar ? begin : -1, int(scalarData->supports.size())}});
				scalarData->allScalar = scalarData->allScalar && scalar;
			}
			scalarMappings = scalarData;
			auto mutableNeighbors = std::make_shared<std::vector<std::set<int>>>(n);
			baseNeighbors = mutableNeighbors;
			auto &surfaceNeighbors = *mutableNeighbors;
			for (const auto &face : surface)
			{
				std::set<int> nodes;
				for (int v : face)
				{
					if (v < 0 || v >= int(surfaceVertices.size()))
					{
						error = "invalid_surface_triangle";
						return false;
					}
					for (const auto &a : surfaceVertices[v].support)
						nodes.insert(a.node);
				}
				// A mapped triangle can couple vertices in different tetrahedra.
				for (int i : nodes)
					surfaceNeighbors[i].insert(nodes.begin(), nodes.end());
			}
		}
		if (surfaceOwners.empty())
			surfaceOwners.assign(surface.size(), 0);
		if (surfaceFrictions.empty())
			surfaceFrictions.assign(surface.size(), .5);
		if (surfaceOwners.size() != surface.size() || surfaceFrictions.size() != surface.size())
		{
			error = "invalid_surface_properties";
			return false;
		}
		colorVertices(false);
		return true;
	}
	bool step(btScalar dt)
	{
		error = nullptr;
		bool finite = std::isfinite(double(dt));
		for (btScalar value : {settings.radius, settings.gap, settings.ke, settings.kd, settings.frictionEpsilon, settings.volumeFloor,
							   settings.recoveryDistance, settings.collisionUnitsPerMeter})
			finite = finite && std::isfinite(double(value));
		if (!finite || !(dt > 0) || settings.iterations < 1 || settings.iterations > 1000 || settings.workers < 1 ||
			settings.workers > 256 || !(settings.radius > 0) || !(settings.gap >= settings.radius) || !(settings.frictionEpsilon > 0) ||
			settings.ke < 0 || settings.kd < 0 || !(settings.volumeFloor > 0 && settings.volumeFloor < 1) ||
			settings.recoveryDistance < 0 || !(settings.collisionUnitsPerMeter > 0))
		{
			error = "invalid_settings";
			return false;
		}
		previous = x;
		for (auto &r : rigids)
		{
			r.previous = r.pose;
			r.predicted = r.pose;
			if (!(r.mass > 0))
			{
				r.previous.setOrigin(r.pose.getOrigin() - r.velocity * dt);
				r.previous.setRotation((rotationIncrement(-r.angularVelocity * dt) * r.pose.getRotation()).normalized());
			}
			if (r.mass > 0)
			{
				const auto basis = r.pose.getBasis();
				const btMatrix3x3 inverseInertia =
					basis * btMatrix3x3(1 / r.inertia.x(), 0, 0, 0, 1 / r.inertia.y(), 0, 0, 0, 1 / r.inertia.z()) * basis.transpose();
				const btVector3 v = (r.velocity * btPow(btScalar(1) - r.linearDamping, dt) + r.force * (dt / r.mass)) * r.linearFactor;
				const btVector3 w =
					(r.angularVelocity * btPow(btScalar(1) - r.angularDamping, dt) + inverseInertia * r.torque * dt) * r.angularFactor;
				r.predicted.setOrigin(r.pose.getOrigin() + v * dt);
				r.predicted.setRotation((rotationIncrement(w * dt) * r.pose.getRotation()).normalized());
			}
		}
		recoveredIntersections = 0;
		// Validation and initial contacts see identical geometry; share this one snapshot.
		btVbdPairSet initialPairs = collisionPairs();
		const bool reuseInitialPairs = validSurface(&initialPairs);
		if (!reuseInitialPairs)
		{
			if (!recoverAcceptedPositions() && !recoverConvexIntersections())
			{
				error = "initial_surface_intersection";
				return false;
			}
			previous = x;
		}
		if (!detect(reuseInitialPairs ? &initialPairs : nullptr))
			return false;
		initialPairs.clear();
		inertia = x;
		for (int i = 0; i < int(x.size()); ++i)
			if (mass[i] > 0)
				inertia[i] += velocity[i] * dt + (external[i] / mass[i] - velocity[i] * massDamping[i]) * (dt * dt);
		std::vector<btVector3> proposed = inertia;
		truncate(proposed);
		x = proposed;
		for (int r = 0; r < int(rigids.size()); ++r)
			if (rigids[r].mass > 0)
				moveRigid(r, rigids[r].predicted);
		if (!detect())
			return false;
		for (int iteration = 0; iteration < settings.iterations; ++iteration)
		{
			for (const auto &group : groups)
			{
				proposed = x;
				btVbdParallelFor(int(group.size()), settings.workers,
								 [&](int item)
								 {
									 const int i = group[item];
									 if (mass[i] > 0)
									 {
										 const btScalar inertiaScale = mass[i] / (dt * dt);
										 btVector3 force = (inertia[i] - x[i]) * inertiaScale;
										 btMatrix3x3 hessian = btMatrix3x3::getIdentity() * inertiaScale;
										 for (const auto &a : adjacency[i])
										 {
											 const Tet &tet = tets[a.first];
											 btVector3 p[4], p0[4];
											 for (int j = 0; j < 4; ++j)
											 {
												 p[j] = x[tet.nodes[j]];
												 p0[j] = previous[tet.nodes[j]];
											 }
											 btVbd::tetraBlock(p, p0, tet.inverseRest, tet.volume, a.second, tet.mu, tet.lambda,
															   tet.damping, dt, force, hessian);
										 }
										 for (int k : springAdjacency[i])
											 springBlock(springs[k], x[i], previous[i], dt, force, hessian);
										 for (int k : attachmentAdjacency[i])
										 {
											 const auto &a = attachments[k];
											 const auto &r = rigids[a.rigid];
											 const btVector3 target = r.pose * a.localPoint, oldTarget = r.previous * a.localPoint;
											 force -= (x[i] - target) * a.stiffness +
													  ((x[i] - previous[i]) - (target - oldTarget)) * (a.damping / dt);
											 hessian += btMatrix3x3::getIdentity() * (a.stiffness + a.damping / dt);
										 }
										 for (int k : contactAdjacency[i])
										 {
											 const Contact &c = contacts[k];
											 btMatrix3x3 jacobian = btMatrix3x3::getIdentity() * btScalar(0);
											 for (const auto &a : c.vertex.support)
												 if (a.node == i)
													 jacobian += a.jacobian;
											 btVector3 cf(0, 0, 0);
											 btMatrix3x3 ch = btMatrix3x3::getIdentity() * btScalar(0);
											 contactBlock(c, dt, cf, ch);
											 force += jacobian.transpose() * cf;
											 hessian += jacobian.transpose() * ch * jacobian;
										 }
										 if (btFabs(hessian.determinant()) > 1e-8)
											 proposed[i] += hessian.inverse() * force;
									 }
								 });
				truncate(proposed);
				x = proposed;
			}
			solveRigids(dt);
			solveJoints(dt);
		}
		minimumJ = SIMD_INFINITY;
		if (!validSurface())
		{
			error = "intersecting_surface";
			return false;
		}
		for (const Tet &tet : tets)
		{
			btVector3 p[4];
			for (int j = 0; j < 4; ++j)
				p[j] = x[tet.nodes[j]];
			const btScalar J = btVbd::deformation(p, tet.inverseRest).determinant();
			if (!(J > 0) || !std::isfinite(double(J)))
			{
				error = "invalid_final_tet";
				return false;
			}
			minimumJ = btMin(minimumJ, J);
		}
		for (int i = 0; i < int(x.size()); ++i)
		{
			velocity[i] = (x[i] - previous[i]) / dt;
			if (!std::isfinite(double(velocity[i].length2())))
			{
				error = "nonfinite_velocity";
				return false;
			}
		}
		for (auto &r : rigids)
			if (r.mass > 0)
			{
				r.velocity = (r.pose.getOrigin() - r.previous.getOrigin()) / dt;
				r.angularVelocity = rotationVector(r.pose.getRotation() * r.previous.getRotation().inverse()) / dt;
			}
		return true;
	}

  private:
	std::unique_ptr<btDeformableVbdCollisionCache> ownedCollisionCache;
	btDeformableVbdCollisionMesh &softMesh, &barrierMesh;
	btDeformableVbdGpu &gpu;
	std::vector<std::vector<std::pair<int, int>>> adjacency;
	std::vector<std::vector<int>> groups, contactAdjacency, springAdjacency, attachmentAdjacency;
	std::vector<int> colors;
	std::vector<btVector3> previous, reference, inertia;
	using GuardPosition = btDeformableVbdGuardPosition;
	std::vector<GuardPosition> &guardPositions;
	std::vector<unsigned char> guardChangedNodes;
	std::vector<std::vector<int>> guardSurfaceVertices;
	std::vector<int> guardAffectedVertices;
	std::vector<std::vector<int>> guardPlaneNodes;
	std::vector<unsigned int> guardPlaneStamps;
	std::vector<int> guardAffectedPlanes;
	void affectedGuardPlanes()
	{
		if (guardPlaneNodes.empty())
		{
			guardPlaneNodes.resize(x.size());
			guardPlaneStamps.assign(planes.size(), 0);
			std::vector<int> seen(x.size(), -1);
			for (int p = 0; p < int(planes.size()); ++p)
				for (int v : planes[p].nodes)
					for (const auto &support : surfaceVertices[v].support)
						if (seen[support.node] != p)
						{
							seen[support.node] = p;
							guardPlaneNodes[support.node].push_back(p);
						}
		}
		guardAffectedPlanes.clear();
		for (int node = 0; node < int(x.size()); ++node)
			if (guardChangedNodes[node])
				for (int p : guardPlaneNodes[node])
					if (guardPlaneStamps[p] != guardStamp)
					{
						guardPlaneStamps[p] = guardStamp;
						guardAffectedPlanes.push_back(p);
					}
		// Each feature only tightens the common minimum; selection order does not change its bound.
	}

	std::vector<int> guardTets;
	std::vector<btScalar> guardCommonLimits;
	std::vector<std::uint32_t> guardTetMasks;
	bool selectGuardTets(const std::vector<btVector3> &proposed)
	{
		if (tets.size() < 512)
			return false;
		// Bounded coordinates guarantee finite determinants for unchanged tetrahedra.
		static const btScalar safeCoordinate = btScalar(std::cbrt(double(std::numeric_limits<btScalar>::max())) / 8);
		for (const auto &p : x)
			for (int d = 0; d < 3; ++d)
				if (!std::isfinite(double(p[d])) || btFabs(p[d]) > safeCoordinate)
					return false;
		guardTetMasks.assign((tets.size() + 31) / 32, 0);
		for (int node = 0; node < int(x.size()); ++node)
			if (proposed[node] != x[node])
				for (const auto &entry : adjacency[node])
					guardTetMasks[entry.first / 32] |= std::uint32_t(1) << (entry.first % 32);
		guardTets.clear();
		// Enumerate marked tetrahedra in original order without sorting the affected set.
		static const unsigned char bitIndex[32] = {0,  1,  28, 2,  29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4,  8,
												   31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6,  11, 5,  10, 9};
		for (int word = 0; word < int(guardTetMasks.size()); ++word)
		{
			std::uint32_t mask = guardTetMasks[word];
			while (mask)
			{
				const std::uint32_t lowest = mask & (std::uint32_t(0) - mask);
				guardTets.push_back(word * 32 + bitIndex[(lowest * std::uint32_t(0x077CB531)) >> 27]);
				mask &= mask - 1;
			}
		}
		return true;
	}
	std::vector<btScalar> guardVertexLimits;
	unsigned int &guardStamp, &guardReferenceStamp;
	btVector3 surfacePosition(int vertex, const std::vector<btVector3> &positions) const
	{
		const auto &mapping = surfaceVertices[vertex];
		if (!scalarMappings || scalarMappings->ranges[vertex][0] < 0)
			return mapping.position(positions);
		const auto &range = scalarMappings->ranges[vertex];
		const auto *data = scalarMappings->supports.data();
		const btVector3 result = btVbd::scalarMappingPosition(mapping.offset, positions.data(), data + range[0], data + range[1]);
		// Preserve the matrix path's invalid-input and overflow propagation.
		if (!std::isfinite(double(result.x())) || !std::isfinite(double(result.y())) || !std::isfinite(double(result.z())))
			return mapping.position(positions);
		return result;
	}
	GuardPosition &guardPosition(int vertex, const std::vector<btVector3> &proposed)
	{
		auto &entry = guardPositions[vertex];
		if (entry.stamp != guardStamp)
		{
			const auto &mapping = surfaceVertices[vertex];
			entry.changed = false;
			for (const auto &support : mapping.support)
				entry.changed = entry.changed || guardChangedNodes[support.node];
			entry.current = surfacePosition(vertex, x);
			entry.proposed = entry.changed ? surfacePosition(vertex, proposed) : entry.current;
			if (entry.referenceStamp != guardReferenceStamp)
			{
				entry.reference = surfacePosition(vertex, reference);
				entry.referenceStamp = guardReferenceStamp;
			}
			entry.stamp = guardStamp;
		}
		return entry;
	}
	struct ContactKeyHash
	{
		size_t operator()(const std::vector<long long> &key) const
		{
			size_t hash = 0;
			for (auto value : key)
				hash ^= std::hash<long long>{}(value) + size_t(0x9e3779b9) + (hash << 6) + (hash >> 2);
			return hash;
		}
	};
	std::unordered_set<std::vector<long long>, ContactKeyHash> contactKeys;
	std::vector<Support> contactSupportScratch;
	using PairDistance = btDeformableVbdPairDistance;
	std::vector<PairDistance> &pairDistances;
	std::vector<GIM_PAIR> &pairSortScratch;
	std::vector<unsigned char> validationBlockFlags, nativeCollectedNodes;
	std::vector<std::vector<unsigned char>> nativeSourceCollected;
	std::vector<std::vector<int>> nativeSourceCandidates;
	std::vector<int> nativeRigidTriangles;

	bool gpuPrepared = false, gpuReferenceValid = false;
	unsigned int gpuReferenceStamp = 0;
	std::vector<btVbdGpuVec> gpuCurrent, gpuProposed, gpuReference;
	bool prepareGpuMapping()
	{
		if (sizeof(btScalar) != sizeof(double))
		{
			gpuError = "GPU guard prototype requires double-precision Bullet";
			return false;
		}
		if (!gpu.ready())
		{
			gpuError = gpu.error();
			return false;
		}
		if (!gpuPrepared)
		{
			std::vector<btVbdGpuMapping> maps(surfaceVertices.size());
			std::vector<btVbdGpuSupport> supports;
			supports.reserve(surfaceVertices.size() * 4);
			for (int v = 0; v < int(surfaceVertices.size()); ++v)
			{
				const auto &source = surfaceVertices[v];
				auto &target = maps[v];
				target.begin = int(supports.size());
				target.offset = {double(source.offset.x()), double(source.offset.y()), double(source.offset.z())};
				for (const auto &support : source.support)
				{
					btVbdGpuSupport item{};
					item.node = support.node;
					for (int row = 0; row < 3; ++row)
						for (int col = 0; col < 3; ++col)
							item.j[3 * row + col] = double(support.jacobian[row][col]);
					supports.push_back(item);
				}
				target.end = int(supports.size());
			}
			if (!gpu.mapping(maps, supports))
			{
				gpuError = gpu.error();
				return false;
			}
			gpuPrepared = true;
		}
		return true;
	}
	bool gpuMappedBound(const std::vector<btVector3> &proposed, btScalar bound, btScalar &limit)
	{
		if (!prepareGpuMapping())
		{
			++gpuGuardFallbacks;
			return false;
		}
		const bool refresh = !gpuReferenceValid || gpuReferenceStamp != guardReferenceStamp;
		gpuCurrent.resize(x.size());
		gpuProposed.resize(x.size());
		gpuReference.resize(x.size());
		for (int i = 0; i < int(x.size()); ++i)
		{
			gpuCurrent[i] = {double(x[i].x()), double(x[i].y()), double(x[i].z())};
			gpuProposed[i] = {double(proposed[i].x()), double(proposed[i].y()), double(proposed[i].z())};
			if (refresh)
				gpuReference[i] = {double(reference[i].x()), double(reference[i].y()), double(reference[i].z())};
		}
		double result = 1;
		if (!gpu.guard(gpuCurrent, gpuProposed, gpuReference, refresh, double(bound), result))
		{
			gpuError = gpu.error();
			++gpuGuardFallbacks;
			return false;
		}
		gpuReferenceStamp = guardReferenceStamp;
		gpuReferenceValid = true;
		// The CUDA kernel uses double precision without fused multiply-add contractions.
		limit = btMin(limit, btScalar(result));
		++gpuGuardCalls;
		return true;
	}
	btScalar mappingAmplification = 0;
	std::shared_ptr<const std::vector<std::set<int>>> baseNeighbors;
	std::shared_ptr<const ScalarMappingData> scalarMappings;
	void colorVertices(bool includeContacts)
	{
		const int n = int(x.size());
		std::vector<std::set<int>> contactNeighbors;
		const auto *neighbors = baseNeighbors.get();
		if (includeContacts)
		{
			contactNeighbors = *baseNeighbors;
			neighbors = &contactNeighbors;
		}
		if (includeContacts)
			for (const auto &contact : contacts)
				for (const auto &a : contact.vertex.support)
					for (const auto &b : contact.vertex.support)
						contactNeighbors[a.node].insert(b.node);
		colors.assign(n, -1);
		groups.clear();
		std::vector<int> used(n, -1);
		for (int i = 0; i < n; ++i)
		{
			for (int j : (*neighbors)[i])
				if (colors[j] >= 0)
					used[colors[j]] = i;
			for (const auto &a : adjacency[i])
				for (int j : tets[a.first].nodes)
					if (colors[j] >= 0)
						used[colors[j]] = i;
			int color = 0;
			while (color < n && used[color] == i)
				++color;
			colors[i] = color;
			if (int(groups.size()) <= color)
				groups.resize(color + 1);
			groups[color].push_back(i);
		}
		colorCount = int(groups.size());
	}
	bool allowed(int a, int b) const
	{
		return collisionAllowed.empty() ||
			   (a >= 0 && b >= 0 && a < int(collisionAllowed.size()) && b < int(collisionAllowed[a].size()) && collisionAllowed[a][b]);
	}
	bool movingPairAllowed(int a, int b) const
	{
		if (a >= b || !allowed(surfaceOwners[a], surfaceOwners[b]))
			return false;
		if (surfaceOwners[a] != surfaceOwners[b])
			return true;
		const int owner = surfaceOwners[a];
		if (owner < 0 || owner >= int(selfContact.size()) || !selfContact[owner])
			return false;
		// Neighbouring embedded triangles share material support; they are not self-contact candidates.
		for (int va : surface[a])
			for (int vb : surface[b])
				for (const auto &na : surfaceVertices[va].support)
					for (const auto &nb : surfaceVertices[vb].support)
						if (na.node == nb.node)
							return false;
		return true;
	}
	void movingPairsRecursive(int a, int b, const std::vector<int> &owners, btVbdPairSet &pairs)
	{
		const int mixed = (-2147483647 - 1), ownerA = owners[a], ownerB = owners[b];
		if (ownerA != mixed && ownerB != mixed)
		{
			if (!allowed(ownerA, ownerB) && !allowed(ownerB, ownerA))
				return;
			if (ownerA == ownerB && (ownerA < 0 || ownerA >= int(selfContact.size()) || !selfContact[ownerA]))
				return;
		}
		auto &tree = softMesh.tree;
		btAABB boundsA, boundsB;
		tree.getNodeBound(a, boundsA);
		tree.getNodeBound(b, boundsB);
		if (!boundsA.has_collision(boundsB))
			return;
		const bool leafA = tree.isLeafNode(a), leafB = tree.isLeafNode(b);
		if (leafA && leafB)
		{
			const int i = btMin(tree.getNodeData(a), tree.getNodeData(b));
			const int j = btMax(tree.getNodeData(a), tree.getNodeData(b));
			if (movingPairAllowed(i, j))
				pairs.push_back({i, j});
		}
		else if (a == b)
		{
			const int left = tree.getLeftNode(a), right = tree.getRightNode(a);
			movingPairsRecursive(left, left, owners, pairs);
			movingPairsRecursive(left, right, owners, pairs);
			movingPairsRecursive(right, right, owners, pairs);
		}
		else if (!leafA)
		{
			movingPairsRecursive(tree.getLeftNode(a), b, owners, pairs);
			movingPairsRecursive(tree.getRightNode(a), b, owners, pairs);
		}
		else
		{
			movingPairsRecursive(a, tree.getLeftNode(b), owners, pairs);
			movingPairsRecursive(a, tree.getRightNode(b), owners, pairs);
		}
	}
	btVbdPairSet movingPairs()
	{
		btVbdPairSet pairs;
		if (surface.empty())
			return pairs;
		// Reject disabled same-body subtrees before generating dense render-mesh self pairs.
		auto &tree = softMesh.tree;
		const auto &owners = softMesh.subtreeOwners(surfaceOwners);
		if (!owners.empty())
			movingPairsRecursive(0, 0, owners, pairs);
		return pairs;
	}
	SurfaceVertex interpolate(const std::array<int, 3> &face, const btVector3 &weights) const
	{
		SurfaceVertex result;
		std::map<int, btMatrix3x3> support;
		for (int j = 0; j < 3; ++j)
		{
			const auto &v = surfaceVertices[face[j]];
			result.offset += v.offset * weights[j];
			for (const auto &a : v.support)
			{
				auto found = support.find(a.node);
				if (found == support.end())
					support.emplace(a.node, a.jacobian * weights[j]);
				else
					found->second += a.jacobian * weights[j];
			}
		}
		for (const auto &a : support)
			if (a.second[0].length2() + a.second[1].length2() + a.second[2].length2() > 0)
				result.support.push_back({a.first, a.second});
		return result;
	}

	btVector3 planePoint(const Plane &plane) const
	{
		if (plane.rigid < 0)
			return plane.point;
		const auto &rigid = rigids[plane.rigid];
		if (plane.triangle >= 0)
		{
			const auto &triangle = barriers[plane.triangle];
			btVector3 point = rigid.pose * triangle.local[0];
			for (int j = 1; j < 3; ++j)
			{
				const btVector3 p = rigid.pose * triangle.local[j];
				if (plane.normal.dot(p) > plane.normal.dot(point))
					point = p;
			}
			return point;
		}
		const auto &shape = convexBarriers[plane.convex];
		btTransform pose = rigid.pose;
		pose.setOrigin(pose.getOrigin() * settings.collisionUnitsPerMeter);
		pose *= shape.localTransform;
		return (pose * shape.shape->localGetSupportingVertex(pose.getBasis().transpose() * plane.normal)) / settings.collisionUnitsPerMeter;
	}
	void setPlaneGap(Plane &plane)
	{
		plane.referenceGap = SIMD_INFINITY;
		const auto point = planePoint(plane);
		for (int v : plane.nodes)
			plane.referenceGap = btMin(plane.referenceGap, plane.normal.dot(surfaceVertices[v].position(reference) - point));
	}
	btScalar planeGap(const Plane &plane) const
	{
		btScalar gap = SIMD_INFINITY;
		const auto point = planePoint(plane);
		for (int v : plane.nodes)
			gap = btMin(gap, plane.normal.dot(surfaceVertices[v].position(x) - point));
		return gap;
	}
	void moveRigid(int index, const btTransform &target)
	{
		auto &r = rigids[index];
		const btTransform old = r.pose;
		std::vector<std::pair<int, btScalar>> required;
		for (int p = 0; p < int(planes.size()); ++p)
			if (planes[p].rigid == index)
				required.push_back(
					{p, btMin(planeGap(planes[p]), btMax(btScalar(0), planes[p].referenceGap) * btScalar(.05) + btScalar(1e-6))});
		std::vector<std::pair<int, btScalar>> contactRequired;
		for (int k = 0; k < int(contacts.size()); ++k)
			if (contacts[k].vertex.support.empty() && (contacts[k].rigid == index || contacts[k].rigidA == index))
				contactRequired.push_back(
					{k, btMin(contactGap(contacts[k]), btMax(btScalar(0), contacts[k].referenceGap) * btScalar(.05) + btScalar(1e-6))});
		btScalar fraction = 1;
		for (int attempt = 0; attempt < 24; ++attempt)
		{
			r.pose.setOrigin(old.getOrigin() + (target.getOrigin() - old.getOrigin()) * fraction);
			r.pose.setRotation(old.getRotation().slerp(target.getRotation(), fraction).normalized());
			bool valid = (r.pose.getOrigin() - r.reference.getOrigin()).length() +
							 rotationVector(r.pose.getRotation() * r.reference.getRotation().inverse()).length() * r.radius <=
						 settings.gap * btScalar(.425);
			for (const auto &p : required)
				valid = valid && planeGap(planes[p.first]) + btScalar(1e-12) >= p.second;
			for (const auto &c : contactRequired)
				valid = valid && contactGap(contacts[c.first]) + btScalar(1e-12) >= c.second;
			if (valid)
				return;
			fraction *= btScalar(.5);
		}
		r.pose = old;
	}
	btMatrix3x3 inverseRigidInertia(int index) const
	{
		const auto &r = rigids[index];
		const auto basis = r.pose.getBasis();
		return basis * btMatrix3x3(1 / r.inertia.x(), 0, 0, 0, 1 / r.inertia.y(), 0, 0, 0, 1 / r.inertia.z()) * basis.transpose();
	}
	void solveJoints(btScalar dt)
	{
		// Bullet's constraint rows retain their limits and motors; project impulses in SI units.
		for (auto &row : joints)
		{
			btScalar response = row.cfm, speed = 0;
			for (int side = 0; side < 2; ++side)
			{
				const int index = side ? row.b : row.a;
				if (index < 0)
					continue;
				const auto &r = rigids[index];
				const auto &linear = side ? row.linearB : row.linearA;
				const auto &angular = side ? row.angularB : row.angularA;
				speed += (linear.dot(r.pose.getOrigin() - r.previous.getOrigin()) +
						  angular.dot(rotationVector(r.pose.getRotation() * r.previous.getRotation().inverse()))) /
						 dt;
				if (r.mass > 0)
					response += linear.dot(linear * r.linearFactor) / r.mass +
								angular.dot((inverseRigidInertia(index) * angular) * r.angularFactor);
			}
			if (!(response > 1e-18))
				continue;
			const btScalar impulse = btClamped(row.impulse + (row.rhs - speed - row.cfm * row.impulse) / response, row.lower, row.upper);
			const btScalar change = impulse - row.impulse;
			row.impulse = impulse;
			for (int side = 0; side < 2; ++side)
			{
				const int index = side ? row.b : row.a;
				if (index < 0 || !(rigids[index].mass > 0))
					continue;
				const auto &r = rigids[index];
				const auto &linear = side ? row.linearB : row.linearA;
				const auto &angular = side ? row.angularB : row.angularA;
				btTransform target = r.pose;
				target.setOrigin(target.getOrigin() + linear * r.linearFactor * (change * dt / r.mass));
				target.setRotation(
					(rotationIncrement((inverseRigidInertia(index) * angular) * r.angularFactor * (change * dt)) * target.getRotation())
						.normalized());
				moveRigid(index, target);
			}
		}
	}
	void solveRigids(btScalar dt)
	{
		for (int index = 0; index < int(rigids.size()); ++index)
		{
			auto &r = rigids[index];
			if (!(r.mass > 0))
				continue;
			const auto basis = r.pose.getBasis();
			const btMatrix3x3 inertia =
				basis * btMatrix3x3(r.inertia.x(), 0, 0, 0, r.inertia.y(), 0, 0, 0, r.inertia.z()) * basis.transpose();
			btVector3 force = (r.predicted.getOrigin() - r.pose.getOrigin()) * (r.mass / (dt * dt));
			btVector3 torque = inertia * rotationVector(r.predicted.getRotation() * r.pose.getRotation().inverse()) / (dt * dt);
			btMatrix3x3 linear = btMatrix3x3::getIdentity() * (r.mass / (dt * dt)), angular = inertia * (1 / (dt * dt));
			for (const auto &c : contacts)
				if (c.rigid == index || c.rigidA == index)
				{
					btVector3 f(0, 0, 0);
					btMatrix3x3 h = btMatrix3x3::getIdentity() * btScalar(0);
					contactBlock(c, dt, f, h);
					const bool sideA = c.rigidA == index;
					const btMatrix3x3 jacobian = pointRotationJacobian(r.pose.getBasis() * (sideA ? c.localA : c.localPoint));
					const btScalar sign = sideA ? 1 : -1;
					force += f * sign;
					torque += jacobian.transpose() * f * sign;
					linear += h;
					angular += jacobian.transpose() * h * jacobian;
				}
			for (const auto &a : attachments)
				if (a.rigid == index)
				{
					const btVector3 point = r.pose * a.localPoint, oldPoint = r.previous * a.localPoint;
					const btVector3 f =
						(x[a.node] - point) * a.stiffness + ((x[a.node] - previous[a.node]) - (point - oldPoint)) * (a.damping / dt);
					const btMatrix3x3 jacobian = pointRotationJacobian(r.pose.getBasis() * a.localPoint);
					const btScalar h = a.stiffness + a.damping / dt;
					force += f;
					torque += jacobian.transpose() * f;
					linear += btMatrix3x3::getIdentity() * h;
					angular += jacobian.transpose() * jacobian * h;
				}
			for (const auto &row : joints)
			{
				if (row.a == index)
				{
					force += row.linearA * (row.impulse / dt);
					torque += row.angularA * (row.impulse / dt);
				}
				if (row.b == index)
				{
					force += row.linearB * (row.impulse / dt);
					torque += row.angularB * (row.impulse / dt);
				}
			}
			btTransform target = r.pose;
			if (btFabs(linear.determinant()) > 1e-12)
				target.setOrigin(target.getOrigin() + (linear.inverse() * force) * r.linearFactor);
			if (btFabs(angular.determinant()) > 1e-18)
				target.setRotation((rotationIncrement((angular.inverse() * torque) * r.angularFactor) * target.getRotation()).normalized());
			moveRigid(index, target);
		}
	}
	void addContact(std::array<int, 3> vertices, btVector3 weights, const btVector3 &point, btScalar friction, int rigid = -1)
	{
		SurfaceVertex vertex;
		// Construction is serial; reuse the small support buffer instead of allocating tree nodes.
		auto &support = contactSupportScratch;
		support.clear();
		for (int j = 0; j < 3; ++j)
		{
			const auto &v = surfaceVertices[vertices[j]];
			vertex.offset += v.offset * weights[j];
			for (const auto &a : v.support)
			{
				auto it = std::find_if(support.begin(), support.end(), [&](const Support &entry) { return entry.node == a.node; });
				if (it == support.end())
					support.push_back({a.node, a.jacobian * weights[j]});
				else
					it->jacobian += a.jacobian * weights[j];
			}
		}
		std::sort(support.begin(), support.end(), [](const Support &a, const Support &b) { return a.node < b.node; });
		vertex.support.reserve(support.size());
		for (const auto &a : support)
			if (a.jacobian[0].length2() + a.jacobian[1].length2() + a.jacobian[2].length2() > 0)
				vertex.support.push_back(a);
		const btVector3 delta = vertex.position(x) - point;
		const btScalar length = delta.length();
		if (length >= settings.gap)
			return;
		if (length < 1e-12)
		{
			error = "initial_touch_requires_contact_normal";
			return;
		}
		const btVector3 normal = delta / length;
		std::vector<long long> key;
		key.reserve(10 + 5 * vertex.support.size());
		key.push_back(rigid);
		for (int d = 0; d < 3; ++d)
		{
			key.push_back(std::llround(double(point[d]) * 1e8));
			key.push_back(std::llround(double(normal[d]) * 1e6));
			key.push_back(std::llround(double(vertex.offset[d]) * 1e8));
		}
		for (const auto &a : vertex.support)
		{
			key.push_back(a.node);
			// Encode the same quantized matrix without storing its zero entries.
			const size_t maskIndex = key.size();
			key.push_back(0);
			long long mask = 0;
			for (int row = 0; row < 3; ++row)
				for (int col = 0; col < 3; ++col)
				{
					const btScalar component = a.jacobian[row][col];
					const long long value = component == 0 ? 0 : std::llround(double(component) * 1e8);
					if (value)
					{
						mask |= 1LL << (3 * row + col);
						key.push_back(value);
					}
				}
			key[maskIndex] = mask;
		}
		if (contactKeys.insert(std::move(key)).second)
		{
			Contact contact{std::move(vertex), point, normal, friction};
			contact.rigid = rigid;
			if (rigid >= 0)
				contact.localPoint = rigids[rigid].pose.invXform(point);
			contacts.push_back(std::move(contact));
		}
	}
	int nativeTriangle(NativeBarrier &source, int primitive)
	{
		// Native primitive IDs are dense; allocate only when this source contributes.
		if (primitive >= int(source.triangles.size()))
			source.triangles.resize(source.shape->getPrimitiveManager()->get_primitive_count(), -1);
		if (source.triangles[primitive] >= 0)
			return source.triangles[primitive];
		btPrimitiveTriangle input;
		source.shape->getPrimitiveManager()->get_primitive_triangle(primitive, input, false);
		Triangle output;
		for (int j = 0; j < 3; ++j)
			output.x[j] = (source.transform * input.m_vertices[j]) / settings.collisionUnitsPerMeter;
		if ((output.x[1] - output.x[0]).cross(output.x[2] - output.x[0]).length2() < 1e-30)
			return -1;
		output.friction = source.friction;
		output.owner = source.owner;
		const int index = int(barriers.size());
		barriers.push_back(output);
		source.triangles[primitive] = index;
		return index;
	}
	std::vector<int> nativeExtractOffsets;
	std::vector<std::array<btVector3, 3>> nativeExtractVertices;
	std::vector<unsigned char> nativeExtractValid;
	bool prepareNativeTriangles()
	{
		nativeExtractOffsets.resize(nativeBarriers.size() + 1);
		int count = 0;
		for (int i = 0; i < int(nativeBarriers.size()); ++i)
		{
			nativeExtractOffsets[i] = count;
			count += int(nativeSourceCandidates[i].size());
		}
		nativeExtractOffsets.back() = count;
		if (settings.workers <= 1 || count < 8192)
			return false;
		nativeExtractVertices.resize(count);
		nativeExtractValid.resize(count);
		// Prepare immutable primitive data before dispatching readers across source boundaries.
		for (int i = 0; i < int(nativeBarriers.size()); ++i)
			if (nativeExtractOffsets[i] != nativeExtractOffsets[i + 1])
			{
				auto &source = nativeBarriers[i];
				source.shape->lockChildShapes();
				const auto *manager = source.shape->getPrimitiveManager();
				if (source.triangles.size() < size_t(manager->get_primitive_count()))
					source.triangles.resize(manager->get_primitive_count(), -1);
				manager->begin_geometry_query();
			}
		btVbdParallelFor((count + 255) / 256, settings.workers,
						 [&](int block)
						 {
							 const int begin = block * 256, end = btMin(count, begin + 256);
							 int sourceIndex = int(std::upper_bound(nativeExtractOffsets.begin(), nativeExtractOffsets.end(), begin) -
												   nativeExtractOffsets.begin()) -
											   1;
							 for (int i = begin; i < end; ++i)
							 {
								 while (i >= nativeExtractOffsets[sourceIndex + 1])
									 ++sourceIndex;
								 const auto &source = nativeBarriers[sourceIndex];
								 const int primitive = nativeSourceCandidates[sourceIndex][i - nativeExtractOffsets[sourceIndex]];
								 nativeExtractValid[i] = 0;
								 if (source.triangles[primitive] >= 0)
									 continue;
								 btPrimitiveTriangle input;
								 source.shape->getPrimitiveManager()->get_primitive_triangle(primitive, input, false);
								 Triangle output;
								 for (int j = 0; j < 3; ++j)
									 output.x[j] =
										 btVbdNativeExtractPosition(source.transform, input.m_vertices[j], settings.collisionUnitsPerMeter);
								 nativeExtractValid[i] = !((output.x[1] - output.x[0]).cross(output.x[2] - output.x[0]).length2() < 1e-30);
								 for (int j = 0; j < 3; ++j)
									 nativeExtractVertices[i][j] = output.x[j];
							 }
						 },
						 2, 1);
		for (int i = int(nativeBarriers.size()) - 1; i >= 0; --i)
			if (nativeExtractOffsets[i] != nativeExtractOffsets[i + 1])
			{
				const auto &source = nativeBarriers[i];
				source.shape->getPrimitiveManager()->end_geometry_query();
				source.shape->unlockChildShapes();
			}
		barriers.reserve(barriers.size() + count);
		return true;
	}
	void appendNativeTriangles(int sourceIndex)
	{
		auto &source = nativeBarriers[sourceIndex];
		for (int i = nativeExtractOffsets[sourceIndex]; i < nativeExtractOffsets[sourceIndex + 1]; ++i)
		{
			const int primitive = nativeSourceCandidates[sourceIndex][i - nativeExtractOffsets[sourceIndex]];
			if (!nativeExtractValid[i] || source.triangles[primitive] >= 0)
				continue;
			Triangle output;
			for (int j = 0; j < 3; ++j)
				output.x[j] = nativeExtractVertices[i][j];
			output.friction = source.friction;
			output.owner = source.owner;
			const int index = int(barriers.size());
			barriers.push_back(output);
			source.triangles[primitive] = index;
		}
	}

	template <class Collect>
	void nativeSoftPairs(const NativeBarrier &source, int softNode, int nativeNode, const std::vector<int> &owners,
						 std::vector<unsigned char> &collected, const Collect &collect)
	{
		if (owners[softNode] != (-2147483647 - 1) && !allowed(owners[softNode], source.owner))
			return;
		const auto &native = *source.shape->getBoxSet();
		if (collected[nativeNode])
			return;
		// Collection needs a triangle only once; later passes still test all soft/rigid pairs.
		if (native.isLeafNode(nativeNode) && native.getNodeData(nativeNode) < int(source.triangles.size()) &&
			source.triangles[native.getNodeData(nativeNode)] >= 0)
		{
			collected[nativeNode] = 1;
			return;
		}
		btAABB a, b;
		softMesh.tree.getNodeBound(softNode, a);
		native.getNodeBound(nativeNode, b);
		b.appy_transform(source.transform);
		const btVector3 padding(settings.gap * settings.collisionUnitsPerMeter, settings.gap * settings.collisionUnitsPerMeter,
								settings.gap * settings.collisionUnitsPerMeter);
		b.m_min -= padding;
		b.m_max += padding;
		if (!a.has_collision(b))
			return;
		const bool softLeaf = softMesh.tree.isLeafNode(softNode), nativeLeaf = native.isLeafNode(nativeNode);
		if (softLeaf && nativeLeaf)
			collected[nativeNode] = collect(native.getNodeData(nativeNode));
		else if (!softLeaf && (nativeLeaf || (a.m_max - a.m_min).length2() > (b.m_max - b.m_min).length2()))
		{
			nativeSoftPairs(source, softMesh.tree.getLeftNode(softNode), nativeNode, owners, collected, collect);
			nativeSoftPairs(source, softMesh.tree.getRightNode(softNode), nativeNode, owners, collected, collect);
		}
		else
		{
			nativeSoftPairs(source, softNode, native.getLeftNode(nativeNode), owners, collected, collect);
			nativeSoftPairs(source, softNode, native.getRightNode(nativeNode), owners, collected, collect);
			collected[nativeNode] = collected[native.getLeftNode(nativeNode)] && collected[native.getRightNode(nativeNode)];
		}
	}
	void collectNativeBarriers()
	{
		if (nativeBarriers.empty())
			return;
		// Native extraction only appends static triangles; collect moving query sources once per pass.
		nativeRigidTriangles.clear();
		for (int k = 0; k < int(barriers.size()); ++k)
			if (barriers[k].rigid >= 0)
				nativeRigidTriangles.push_back(k);
		const auto &owners = softMesh.subtreeOwners(surfaceOwners);
		size_t nativeNodes = 0;
		if (settings.workers > 1 && nativeBarriers.size() > 1 && !owners.empty())
			for (const auto &source : nativeBarriers)
				nativeNodes += source.shape->getBoxSet()->getNodeCount();
		const bool parallelSources = nativeNodes >= 4096;
		if (parallelSources)
		{
			nativeSourceCollected.resize(nativeBarriers.size());
			nativeSourceCandidates.resize(nativeBarriers.size());
			// Queries read immutable trees; extraction below preserves source and triangle order.
			btVbdParallelFor(
				int(nativeBarriers.size()), settings.workers,
				[&](int i)
				{
					const auto &source = nativeBarriers[i];
					const int count = source.shape->getBoxSet()->getNodeCount();
					auto &collected = nativeSourceCollected[i];
					auto &candidates = nativeSourceCandidates[i];
					collected.assign(count, 0);
					candidates.clear();
					if (count)
						nativeSoftPairs(source, 0, 0, owners, collected,
										[&](int primitive)
										{
											candidates.push_back(primitive);
											return true;
										});
				},
				2, 1);
		}
		const bool preparedTriangles = parallelSources && prepareNativeTriangles();
		for (int sourceIndex = 0; sourceIndex < int(nativeBarriers.size()); ++sourceIndex)
		{
			auto &source = nativeBarriers[sourceIndex];
			if (!source.shape->getBoxSet()->getNodeCount())
				continue;
			source.shape->lockChildShapes();
			if (preparedTriangles)
				appendNativeTriangles(sourceIndex);
			else if (parallelSources)
				for (int primitive : nativeSourceCandidates[sourceIndex])
					nativeTriangle(source, primitive);
			else if (!owners.empty())
			{
				// Only collected leaves justify pruning; reset for each source and discovery pass.
				nativeCollectedNodes.assign(source.shape->getBoxSet()->getNodeCount(), 0);
				nativeSoftPairs(source, 0, 0, owners, nativeCollectedNodes,
								[&](int primitive) { return nativeTriangle(source, primitive) >= 0; });
			}
			auto collectBox = [&](btAABB box, int owner)
			{
				if (!allowed(owner, source.owner))
					return;
				const btScalar gap = settings.gap * settings.collisionUnitsPerMeter;
				const btVector3 padding(gap, gap, gap);
				box.m_min -= padding;
				box.m_max += padding;
				box.appy_transform(source.transform.inverse());
				btAlignedObjectArray<int> candidates;
				source.shape->getBoxSet()->boxQuery(box, candidates);
				for (int i = 0; i < candidates.size(); ++i)
					nativeTriangle(source, candidates[i]);
			};
			for (const auto &shape : convexBarriers)
				if (shape.rigid >= 0)
				{
					btAABB box;
					shape.shape->getAabb(shape.transform, box.m_min, box.m_max);
					collectBox(box, shape.owner);
				}
			for (int k : nativeRigidTriangles)
			{
				btAABB box;
				box.calc_from_triangle(barriers[k].x[0] * settings.collisionUnitsPerMeter,
									   barriers[k].x[1] * settings.collisionUnitsPerMeter,
									   barriers[k].x[2] * settings.collisionUnitsPerMeter);
				const int owner = barriers[k].owner;
				collectBox(box, owner);
			}
			source.shape->unlockChildShapes();
		}
	}
	const btVector3 &collisionPosition(int vertex) const
	{
		return *softMesh.candidatePositions.current(vertex);
	}
	bool gpuSurfacePositions()
	{
		if (!prepareGpuMapping())
		{
			++gpuSurfaceFallbacks;
			return false;
		}
		gpuCurrent.resize(x.size());
		for (int i = 0; i < int(x.size()); ++i)
			gpuCurrent[i] = {double(x[i].x()), double(x[i].y()), double(x[i].z())};
		if (!gpu.evaluate(gpuCurrent, int(surfaceVertices.size())))
		{
			gpuError = gpu.error();
			++gpuSurfaceFallbacks;
			return false;
		}
		++gpuSurfaceCalls;
		return true;
	}
	void updateCollisionMeshes()
	{
		// A collision pass owns a snapshot; solver iterations keep evaluating their changing state directly.
		softMesh.candidatePositions.end();
		// Small surfaces do not amortize the CUDA launch and download overhead.
		// Packed tetrahedral interpolation can avoid the GPU round trip with several CPU workers.
		const bool packedCpuSurface = settings.workers >= 4 && scalarMappings && scalarMappings->allScalar &&
									  scalarMappings->supports.size() <= 4 * surfaceVertices.size();
		const bool gpuSurface = settings.gpuGuards && !packedCpuSurface && surfaceVertices.size() >= 32768 && gpuSurfacePositions();
		softMesh.candidatePositions.beginCurrent(
			int(surfaceVertices.size()),
			[&](int i, btVector3 &p)
			{
				if (gpuSurface)
				{
					const auto &v = gpu.positions()[i];
					p = btVector3(v.x, v.y, v.z);
				}
				else
					p = surfacePosition(i, x);
			},
			[&](int count, const auto &operation) { btVbdParallelGeometry(count, settings.workers, operation); });
		for (auto &triangle : barriers)
			if (triangle.rigid >= 0)
				for (int j = 0; j < 3; ++j)
					triangle.x[j] = rigids[triangle.rigid].pose * triangle.local[j];
		for (auto &shape : convexBarriers)
			if (shape.rigid >= 0)
			{
				btTransform pose = rigids[shape.rigid].pose;
				pose.setOrigin(pose.getOrigin() * settings.collisionUnitsPerMeter);
				shape.transform = pose * shape.localTransform;
			}

		const btScalar units = settings.collisionUnitsPerMeter;
		for (int side = 0; side < 2; ++side)
		{
			if (side)
				collectNativeBarriers();
			auto &mesh = side ? barrierMesh : softMesh;
			const size_t count = side ? barriers.size() : surface.size();
			const btScalar margin = settings.radius * units / 2;
			const btScalar padding = (settings.gap - settings.radius) * units / 2;
			const bool changed = mesh.updateTriangles(int(count), margin, padding, settings.workers, [&](int k, int j)
													  { return (side ? barriers[k].x[j] : collisionPosition(surface[k][j])) * units; });
			// Static geometry retains valid bounds, including when only body ownership changes.
			if (changed)
				mesh.update(settings.workers);
		}
	}
	btVbdPairSet collisionPairs()
	{
		updateCollisionMeshes();
		btVbdPairSet pairs;
		if (!surface.empty() && !barriers.empty())
		{
			if (settings.workers > 1 && surface.size() + barriers.size() >= 4096)
				btGImpactBvh::find_collision_parallel(&softMesh.tree, btTransform::getIdentity(), &barrierMesh.tree,
													  btTransform::getIdentity(), pairs, settings.workers * 4,
													  [&](int n, const auto &fn) { btVbdParallelFor(n, settings.workers, fn, 2, 1); });
			else
				btGImpactBvh::find_collision(&softMesh.tree, btTransform::getIdentity(), &barrierMesh.tree, btTransform::getIdentity(),
											 pairs);
		}
		// Contact reduction must not depend on BVH layout or unrelated distant triangles.
		btVbdSortCollisionPairsParallel(pairs, pairSortScratch, settings.workers);
		return pairs;
	}
	btVector3 contactSupport(const Contact &c, bool sideA) const
	{
		const int convex = sideA ? c.convexA : c.convexB, triangle = sideA ? c.triangleA : c.triangleB;
		const btVector3 direction = c.normal * (sideA ? btScalar(-1) : btScalar(1));
		if (convex >= 0)
		{
			const auto &shape = convexBarriers[convex];
			btTransform pose = shape.transform;
			if (shape.rigid >= 0)
			{
				pose = rigids[shape.rigid].pose;
				pose.setOrigin(pose.getOrigin() * settings.collisionUnitsPerMeter);
				pose *= shape.localTransform;
			}
			return (pose * shape.shape->localGetSupportingVertex(pose.getBasis().transpose() * direction)) /
				   settings.collisionUnitsPerMeter;
		}
		const auto &face = barriers[triangle];
		btVector3 point = face.rigid < 0 ? face.x[0] : rigids[face.rigid].pose * face.local[0];
		for (int j = 1; j < 3; ++j)
		{
			const auto p = face.rigid < 0 ? face.x[j] : rigids[face.rigid].pose * face.local[j];
			if (p.dot(direction) > point.dot(direction))
				point = p;
		}
		return point;
	}
	btScalar contactGap(const Contact &c) const
	{
		const btVector3 a = (c.convexA >= 0 || c.triangleA >= 0) ? contactSupport(c, true)
																 : (c.rigidA < 0 ? c.vertex.position(x) : rigids[c.rigidA].pose * c.localA);
		const btVector3 b =
			(c.convexB >= 0 || c.triangleB >= 0) ? contactSupport(c, false) : (c.rigid < 0 ? c.point : rigids[c.rigid].pose * c.localPoint);
		return c.normal.dot(a - b);
	}
	bool addRigidContact(int a, int b, const btVector3 &pa, const btVector3 &pb, const btVector3 &normal, btScalar distance,
						 btScalar friction, bool collect)
	{
		if (!std::isfinite(double(distance)) || !(distance > 0))
		{
			error = "rigid_surface_intersection";
			return false;
		}
		if (!collect || distance >= settings.gap)
			return true;
		Contact contact;
		contact.vertex.offset = pa;
		contact.point = pb;
		contact.normal = normal;
		contact.friction = friction;
		contact.rigidA = a;
		contact.rigid = b;
		contact.referenceGap = distance;
		if (a >= 0)
			contact.localA = rigids[a].pose.invXform(pa);
		if (b >= 0)
			contact.localPoint = rigids[b].pose.invXform(pb);
		contacts.push_back(contact);
		return true;
	}
	bool rigidContacts(bool collect)
	{
		const btScalar units = settings.collisionUnitsPerMeter;
		auto query = [&](const btConvexShape *a, const btTransform &ta, int ra, const btConvexShape *b, const btTransform &tb, int rb,
						 btScalar friction, int convexA, int convexB, int triangleB)
		{
			btVector3 alo, ahi, blo, bhi;
			a->getAabb(ta, alo, ahi);
			b->getAabb(tb, blo, bhi);
			const btVector3 padding(settings.gap * units, settings.gap * units, settings.gap * units);
			alo -= padding;
			ahi += padding;
			if (!TestAabbAgainstAabb2(alo, ahi, blo, bhi))
				return true;
			btVoronoiSimplexSolver simplex;
			btGjkEpaPenetrationDepthSolver penetration;
			btGjkPairDetector detector(a, b, &simplex, &penetration);
			btDiscreteCollisionDetectorInterface::ClosestPointInput input;
			input.m_transformA = ta;
			input.m_transformB = tb;
			btPointCollector result;
			detector.getClosestPoints(input, result, nullptr);
			if (!result.m_hasResult)
			{
				error = "rigid_convex_distance_failed";
				return false;
			}
			const btVector3 pb = result.m_pointInWorld / units, pa = pb + result.m_normalOnBInWorld * (result.m_distance / units);
			const auto count = contacts.size();
			if (!addRigidContact(ra, rb, pa, pb, result.m_normalOnBInWorld, result.m_distance / units, friction, collect))
				return false;
			if (contacts.size() > count)
			{
				auto &c = contacts.back();
				c.convexA = convexA;
				c.convexB = convexB;
				c.triangleB = triangleB;
			}
			return true;
		};
		for (int a = 0; a < int(convexBarriers.size()); ++a)
		{
			const auto &ca = convexBarriers[a];
			for (int b = a + 1; b < int(convexBarriers.size()); ++b)
			{
				const auto &cb = convexBarriers[b];
				if ((ca.owner >= 0 && ca.owner == cb.owner) || (ca.rigid >= 0 && ca.rigid == cb.rigid) || (ca.rigid < 0 && cb.rigid < 0) ||
					!allowed(ca.owner, cb.owner))
					continue;
				if (!query(ca.shape, ca.transform, ca.rigid, cb.shape, cb.transform, cb.rigid, ca.friction * cb.friction, a, b, -1))
					return false;
			}
			btAABB box;
			ca.shape->getAabb(ca.transform, box.m_min, box.m_max);
			const btVector3 padding(settings.gap * units, settings.gap * units, settings.gap * units);
			box.m_min -= padding;
			box.m_max += padding;
			btAlignedObjectArray<int> candidates;
			if (!barriers.empty())
				barrierMesh.tree.boxQuery(box, candidates);
			for (int item = 0; item < candidates.size(); ++item)
			{
				const int b = candidates[item];
				const auto &cb = barriers[b];
				if ((ca.owner >= 0 && ca.owner == cb.owner) || (ca.rigid >= 0 && ca.rigid == cb.rigid) || (ca.rigid < 0 && cb.rigid < 0) ||
					!allowed(ca.owner, cb.owner))
					continue;
				const auto &tri = barrierMesh.triangles[b];
				btTriangleShape shape(tri.m_vertices[0], tri.m_vertices[1], tri.m_vertices[2]);
				shape.setMargin(0);
				if (!query(ca.shape, ca.transform, ca.rigid, &shape, btTransform::getIdentity(), cb.rigid, ca.friction * cb.friction, a, -1,
						   b))
					return false;
			}
		}
		for (int a = 0; a < int(barriers.size()); ++a)
			if (barriers[a].rigid >= 0)
			{
				btAABB box;
				barrierMesh.get_primitive_box(a, box);
				btAlignedObjectArray<int> candidates;
				barrierMesh.tree.boxQuery(box, candidates);
				for (int item = 0; item < candidates.size(); ++item)
				{
					const int b = candidates[item];
					if ((barriers[a].owner >= 0 && barriers[a].owner == barriers[b].owner) || (barriers[a].rigid == barriers[b].rigid) ||
						(barriers[b].rigid >= 0 && a >= b) || !allowed(barriers[a].owner, barriers[b].owner))
						continue;
					auto &ta = barrierMesh.triangles[a];
					auto &tb = barrierMesh.triangles[b];
					if (!ta.overlap_test(tb))
						continue;
					btScalar distance2;
					btVector3 pa, pb;
					if (!ta.triangle_triangle_distance(tb, distance2, pa, pb) || !(distance2 > 0))
					{
						error = "rigid_mesh_intersection";
						return false;
					}
					const auto count = contacts.size();
					if (!addRigidContact(barriers[a].rigid, barriers[b].rigid, pa / units, pb / units, (pa - pb).normalized(),
										 btSqrt(distance2) / units, barriers[a].friction * barriers[b].friction, collect))
						return false;
					if (contacts.size() > count)
					{
						contacts.back().triangleA = a;
						contacts.back().triangleB = b;
					}
				}
			}
		return true;
	}
	bool convexContacts(bool collect)
	{
		const btScalar units = settings.collisionUnitsPerMeter;
		for (const auto &barrier : convexBarriers)
		{
			btAABB box;
			barrier.shape->getAabb(barrier.transform, box.m_min, box.m_max);
			const btVector3 padding(settings.gap * units, settings.gap * units, settings.gap * units);
			box.m_min -= padding;
			box.m_max += padding;
			btAlignedObjectArray<int> candidates;
			softMesh.tree.boxQuery(box, candidates);
			for (int c = 0; c < candidates.size(); ++c)
			{
				const int a = candidates[c];
				if (!allowed(surfaceOwners[a], barrier.owner))
					continue;
				const auto &triangle = softMesh.triangles[a];
				btTriangleShape shape(triangle.m_vertices[0], triangle.m_vertices[1], triangle.m_vertices[2]);
				shape.setMargin(0);
				btVoronoiSimplexSolver simplex;
				btGjkEpaPenetrationDepthSolver penetration;
				btGjkPairDetector query(&shape, barrier.shape, &simplex, &penetration);
				btDiscreteCollisionDetectorInterface::ClosestPointInput input;
				input.m_transformA.setIdentity();
				input.m_transformB = barrier.transform;
				btPointCollector result;
				query.getClosestPoints(input, result, nullptr);
				if (!result.m_hasResult || !std::isfinite(double(result.m_distance)))
				{
					error = "convex_distance_failed";
					return false;
				}
				if (!(result.m_distance > 0))
				{
					error = "convex_surface_intersection";
					return false;
				}
				if (!collect || result.m_distance >= settings.gap * units)
					continue;
				const btVector3 normal = result.m_normalOnBInWorld, point = result.m_pointInWorld / units;
				btVector3 vertices[3];
				for (int j = 0; j < 3; ++j)
					vertices[j] = collisionPosition(surface[a][j]);
				const btVector3 onSoft = point + normal * (result.m_distance / units);
				addContact(surface[a], closestWeights(onSoft, vertices), point, barrier.friction * surfaceFrictions[a], barrier.rigid);
				planes.push_back({surface[a], point, normal, barrier.rigid, -1, int(&barrier - convexBarriers.data())});
				setPlaneGap(planes.back());
			}
		}
		return true;
	}
	bool detect(const btVbdPairSet *checkedPairs = nullptr)
	{
		contacts.clear();
		planes.clear();
		guardPlaneNodes.clear();
		movingPlanes.clear();
		contactKeys.clear();
		reference = x;
		if (++guardReferenceStamp == 0)
		{
			for (auto &entry : guardPositions)
				entry.referenceStamp = 0;
			++guardReferenceStamp;
		}
		for (auto &rigid : rigids)
			rigid.reference = rigid.pose;
		const auto freshPairs = checkedPairs ? btVbdPairSet() : collisionPairs();
		const auto &pairs = checkedPairs ? *checkedPairs : freshPairs;

		broadphasePairs = pairs.size();
		auto evaluatePair = [&](int index, PairDistance &result)
		{
			result.status = 0;
			const auto &pair = pairs[index];
			const int a = pair.m_index1, b = pair.m_index2;
			if (!allowed(surfaceOwners[a], barriers[b].owner))
				return;
			auto ta = softMesh.triangles[a];
			const auto &tb = barrierMesh.triangles[b];
			if (!ta.overlap_test(tb))
				return;
			if (!ta.triangle_triangle_distance(tb, result.distance2, result.closestSoft, result.closestRigid))
			{
				result.status = 1;
				return;
			}
			// No last-safe history exists here; ambiguous touches remain errors.
			if (!(result.distance2 > 0) || !std::isfinite(double(result.distance2)))
			{
				result.status = 2;
				return;
			}
			const btScalar search = ta.m_margin + tb.m_margin + ta.m_discoveryPadding + tb.m_discoveryPadding;
			if (result.distance2 >= search * search)
				return;
			result.status = 3;
		};
		const bool parallelDistances = settings.workers > 1 && pairs.size() >= 1024;
		if (parallelDistances)
		{
			pairDistances.resize(pairs.size());
			btVbdParallelGeometry(int(pairs.size()), settings.workers, [&](int i) { evaluatePair(i, pairDistances[i]); });
		}
		// Read-only distances run independently; errors and contacts retain pair order.
		for (int pairIndex = 0; pairIndex < int(pairs.size()); ++pairIndex)
		{
			PairDistance local;
			if (!parallelDistances)
				evaluatePair(pairIndex, local);
			const auto &result = parallelDistances ? pairDistances[pairIndex] : local;
			if (!result.status)
				continue;
			if (result.status != 3)
			{
				error = result.status == 1 ? "gimpact_distance_failed" : "gimpact_intersection_or_ambiguous_touch";
				return false;
			}
			const auto &pair = pairs[pairIndex];
			const int a = pair.m_index1, b = pair.m_index2;
			const auto &closestSoft = result.closestSoft;
			const auto &closestRigid = result.closestRigid;
			const btVector3 normal = (closestSoft - closestRigid).normalized();
			const btVector3 point = closestSoft / settings.collisionUnitsPerMeter;
			const btVector3 other = closestRigid / settings.collisionUnitsPerMeter;
			const auto &face = surface[a];
			const btVector3 vertices[] = {collisionPosition(face[0]), collisionPosition(face[1]), collisionPosition(face[2])};
			addContact(face, closestWeights(point, vertices), other,
					   barriers[b].friction * (barriers[b].owner < 0 ? btScalar(1) : surfaceFrictions[a]), barriers[b].rigid);
			// DAT must protect the entire triangle, even if the closest point lies on one corner.
			planes.push_back({face, other, normal, barriers[b].rigid, b, -1});
			setPlaneGap(planes.back());
		}
		if (!convexContacts(true) || !rigidContacts(true))
			return false;
		const auto dynamicPairs = movingPairs();
		for (const auto &pair : dynamicPairs)
		{
			const int a = pair.m_index1, b = pair.m_index2;
			if (!movingPairAllowed(a, b))
				continue;
			++broadphasePairs;
			auto &ta = softMesh.triangles[a];
			auto &tb = softMesh.triangles[b];
			if (!ta.overlap_test(tb))
				continue;
			btScalar distance2;
			btVector3 pa, pb;
			if (!ta.triangle_triangle_distance(tb, distance2, pa, pb) || !(distance2 > 0) || !std::isfinite(double(distance2)))
			{
				error = "moving_surface_intersection_or_ambiguous_touch";
				return false;
			}
			const btScalar units = settings.collisionUnitsPerMeter;
			if (distance2 >= settings.gap * settings.gap * units * units)
				continue;
			btVector3 av[3], bv[3];
			for (int j = 0; j < 3; ++j)
			{
				av[j] = collisionPosition(surface[a][j]);
				bv[j] = collisionPosition(surface[b][j]);
			}
			SurfaceVertex relative = interpolate(surface[a], closestWeights(pa / units, av));
			const SurfaceVertex other = interpolate(surface[b], closestWeights(pb / units, bv));
			relative.offset -= other.offset;
			for (const auto &entry : other.support)
			{
				bool merged = false;
				for (auto &current : relative.support)
					if (current.node == entry.node)
					{
						current.jacobian -= entry.jacobian;
						merged = true;
						break;
					}
				if (!merged)
					relative.support.push_back({entry.node, entry.jacobian * btScalar(-1)});
			}
			const btVector3 normal = (pa - pb).normalized();
			contacts.push_back({relative, btVector3(0, 0, 0), normal, btMax(btScalar(0), surfaceFrictions[a] * surfaceFrictions[b])});
			movingPlanes.push_back({a, b, normal});
		}
		if (!movingPlanes.empty())
			colorVertices(true);
		contactAdjacency.assign(x.size(), {});
		for (int k = 0; k < int(contacts.size()); ++k)
			for (const auto &a : contacts[k].vertex.support)
				contactAdjacency[a.node].push_back(k);
		return error == nullptr;
	}
	static bool segmentCrossesTriangle(const btVector3 &a, const btVector3 &b, const btVector3 *t)
	{
		const btVector3 d = b - a, e1 = t[1] - t[0], e2 = t[2] - t[0], h = d.cross(e2);
		const btScalar det = e1.dot(h);
		if (btFabs(det) < 1e-12 * btSqrt(d.length2() * e1.length2() * e2.length2()))
			return false;
		const btVector3 s = a - t[0], q = s.cross(e1);
		const btScalar u = s.dot(h) / det, v = d.dot(q) / det, r = e2.dot(q) / det;
		return u >= 0 && v >= 0 && u + v <= 1 && r >= 0 && r <= 1;
	}
	bool recoverAcceptedPositions()
	{
		// Revalidate accepted geometry against today's obstacles; never infer a side from a penetrated triangle.
		if (!(settings.recoveryDistance > 0) || recoveryPositions.size() != x.size())
			return false;
		for (int n = 0; n < int(x.size()); ++n)
			if ((x[n] - recoveryPositions[n]).length() > settings.recoveryDistance || (mass[n] <= 0 && x[n] != recoveryPositions[n]))
				return false;
		const auto original = x;
		x = recoveryPositions;
		for (const auto &tet : tets)
		{
			btVector3 p[4];
			for (int j = 0; j < 4; ++j)
				p[j] = x[tet.nodes[j]];
			if (!(btVbd::deformation(p, tet.inverseRest).determinant() > 0))
			{
				x = original;
				return false;
			}
		}
		error = nullptr;
		if (!validSurface())
		{
			x = original;
			return false;
		}
		recoveredIntersections = 1;
		return true;
	}
	bool recoverConvexIntersections()
	{
		// EPA provides signed convex penetration directions; ambiguous mesh intersections remain rejected.
		if (!(settings.recoveryDistance > 0))
			return false;
		const auto original = x;
		std::vector<btScalar> minimumVolume;
		for (const auto &tet : tets)
		{
			btVector3 p[4];
			for (int j = 0; j < 4; ++j)
				p[j] = x[tet.nodes[j]];
			minimumVolume.push_back(btMin(settings.volumeFloor, btVbd::deformation(p, tet.inverseRest).determinant()));
		}
		for (int pass = 0; pass < 8; ++pass)
		{
			collisionPairs();
			bool changed = false;
			for (const auto &barrier : convexBarriers)
			{
				btAABB box;
				barrier.shape->getAabb(barrier.transform, box.m_min, box.m_max);
				btAlignedObjectArray<int> candidates;
				softMesh.tree.boxQuery(box, candidates);
				for (int item = 0; item < candidates.size(); ++item)
				{
					const int face = candidates[item];
					if (!allowed(surfaceOwners[face], barrier.owner))
						continue;
					btVector3 vertices[3];
					for (int j = 0; j < 3; ++j)
						vertices[j] = surfaceVertices[surface[face][j]].position(x) * settings.collisionUnitsPerMeter;
					btTriangleShape triangle(vertices[0], vertices[1], vertices[2]);
					triangle.setMargin(0);
					btVoronoiSimplexSolver simplex;
					btGjkEpaPenetrationDepthSolver penetration;
					btGjkPairDetector query(&triangle, barrier.shape, &simplex, &penetration);
					btDiscreteCollisionDetectorInterface::ClosestPointInput input;
					input.m_transformA.setIdentity();
					input.m_transformB = barrier.transform;
					btPointCollector result;
					query.getClosestPoints(input, result, nullptr);
					if (!result.m_hasResult || !std::isfinite(double(result.m_distance)))
					{
						x = original;
						return false;
					}
					if (result.m_distance > 0)
						continue;
					if (-result.m_distance / settings.collisionUnitsPerMeter > settings.recoveryDistance)
					{
						x = original;
						return false;
					}
					const btVector3 normal = result.m_normalOnBInWorld;
					if (!std::isfinite(double(normal.length2())) || normal.length2() < .9)
					{
						x = original;
						return false;
					}
					const btScalar plane = normal.dot(result.m_pointInWorld) / settings.collisionUnitsPerMeter +
										   btMin(settings.radius, settings.recoveryDistance) * btScalar(.1);
					for (int corner : surface[face])
					{
						const auto &vertex = surfaceVertices[corner];
						const btScalar depth = plane - normal.dot(vertex.position(x));
						if (depth <= 0)
							continue;
						btScalar response = 0;
						for (const auto &support : vertex.support)
							if (mass[support.node] > 0)
								response += (support.jacobian.transpose() * normal).length2() / mass[support.node];
						if (!(response > 1e-18))
						{
							x = original;
							return false;
						}
						for (const auto &support : vertex.support)
							if (mass[support.node] > 0)
								x[support.node] += support.jacobian.transpose() * normal * (depth / (response * mass[support.node]));
					}
					changed = true;
				}
			}
			bool bounded = true;
			for (int n = 0; n < int(x.size()); ++n)
				bounded = bounded && (x[n] - original[n]).length() <= settings.recoveryDistance;
			for (int k = 0; k < int(tets.size()); ++k)
			{
				const auto &tet = tets[k];
				btVector3 p[4];
				for (int j = 0; j < 4; ++j)
					p[j] = x[tet.nodes[j]];
				const btScalar j = btVbd::deformation(p, tet.inverseRest).determinant();
				bounded = bounded && j > 0 && j + btScalar(1e-12) >= minimumVolume[k];
			}
			if (!bounded || !changed)
			{
				x = original;
				return false;
			}
			error = nullptr;
			if (validSurface())
			{
				recoveredIntersections = 1;
				return true;
			}
		}
		x = original;
		return false;
	}
	template <class Valid> bool validPairs(const btVbdPairSet &pairs, const Valid &valid)
	{
		if (settings.workers <= 1 || pairs.size() < 4096)
		{
			for (const auto &pair : pairs)
				if (!valid(pair))
					return false;
			return true;
		}
		const int batch = 64, count = int(pairs.size()), blocks = (count + batch - 1) / batch;
		validationBlockFlags.resize(blocks);
		btVbdParallelFor(blocks, settings.workers,
						 [&](int block)
						 {
							 const int end = btMin(count, (block + 1) * batch);
							 bool accepted = true;
							 for (int i = block * batch; i < end; ++i)
								 if (!valid(pairs[i]))
								 {
									 accepted = false;
									 break;
								 }
							 validationBlockFlags[block] = accepted;
						 });
		return std::find(validationBlockFlags.begin(), validationBlockFlags.end(), 0) == validationBlockFlags.end();
	}
	bool validSurface(const btVbdPairSet *checkedPairs = nullptr)
	{
		const auto freshPairs = checkedPairs ? btVbdPairSet() : collisionPairs();
		const auto &pairs = checkedPairs ? *checkedPairs : freshPairs;
		// Pair checks only read geometry; error-producing convex/rigid queries stay serial.
		if (!validPairs(pairs,
						[&](const GIM_PAIR &pair)
						{
							if (!allowed(surfaceOwners[pair.m_index1], barriers[pair.m_index2].owner))
								return true;
							const auto &face = surface[pair.m_index1];
							const auto &tri = barriers[pair.m_index2];
							const btVector3 p[] = {collisionPosition(face[0]), collisionPosition(face[1]), collisionPosition(face[2])};
							for (int j = 0; j < 3; ++j)
								if (segmentCrossesTriangle(p[j], p[(j + 1) % 3], tri.x) ||
									segmentCrossesTriangle(tri.x[j], tri.x[(j + 1) % 3], p))
									return false;
							return true;
						}))
			return false;
		if (!convexContacts(false) || !rigidContacts(false))
			return false;
		return validPairs(movingPairs(),
						  [&](const GIM_PAIR &pair)
						  {
							  const int a = pair.m_index1, b = pair.m_index2;
							  if (!movingPairAllowed(a, b))
								  return true;
							  btVector3 av[3], bv[3];
							  for (int j = 0; j < 3; ++j)
							  {
								  av[j] = collisionPosition(surface[a][j]);
								  bv[j] = collisionPosition(surface[b][j]);
							  }
							  for (int j = 0; j < 3; ++j)
								  if (segmentCrossesTriangle(av[j], av[(j + 1) % 3], bv) ||
									  segmentCrossesTriangle(bv[j], bv[(j + 1) % 3], av))
									  return false;
							  return true;
						  });
	}

	void truncate(std::vector<btVector3> &proposed)
	{
		std::vector<btScalar> factors(x.size(), 1);
		for (int i = 0; i < int(x.size()); ++i)
		{
			const btScalar length = (proposed[i] - reference[i]).length(), bound = btScalar(.5 * .85) * settings.gap;
			if (length > bound)
				factors[i] = bound / length;
		}
		if (!mappedSurface && rigids.empty())
			for (const Plane &planeContact : planes)
			{
				const auto &c = planeContact;
				btScalar approach = 0, scale = 1, gap = SIMD_INFINITY;
				for (int i : c.nodes)
				{
					gap = btMin(gap, c.normal.dot(reference[i] - c.point));
					approach = btMax(approach, -c.normal.dot(proposed[i] - reference[i]));
					for (int d = 0; d < 3; ++d)
						scale = btMax(scale, btFabs(reference[i][d]));
				}
				for (int d = 0; d < 3; ++d)
					scale = btMax(scale, btFabs(c.point[d]));
				const btScalar epsilon = btMax(btScalar(1e-6), btScalar(4 * 1.1920928955078125e-7) * scale);
				gap = btMax(btScalar(0), gap);
				btScalar fraction = .5;
				if (gap >= 2 * epsilon && approach > 0)
					fraction = btMax(btScalar(.05), epsilon / gap);
				const btVector3 plane = c.point + c.normal * (fraction * gap);
				for (int i : c.nodes)
					factors[i] = btMin(factors[i], btVbd::planarBound(reference[i], proposed[i] - reference[i], c.normal, plane, epsilon));
			}
		for (int i = 0; i < int(x.size()); ++i)
			proposed[i] = reference[i] + (proposed[i] - reference[i]) * factors[i];
		// Local bounds avoid slowing unrelated vertices; the common bound covers prediction.
		factors.assign(x.size(), 1);
		const bool indexedTets = selectGuardTets(proposed);
		const int guardTetCount = int(indexedTets ? guardTets.size() : tets.size());
		for (int guardTet = 0; guardTet < guardTetCount; ++guardTet)
		{
			const Tet &tet = tets[indexedTets ? guardTets[guardTet] : guardTet];
			btVector3 old[4], next[4];
			int changed = -1, count = 0;
			for (int j = 0; j < 4; ++j)
			{
				int i = tet.nodes[j];
				old[j] = x[i];
				next[j] = proposed[i];
				if ((next[j] - old[j]).length2() > 0)
				{
					changed = i;
					++count;
				}
			}
			if (count == 1)
				factors[changed] = btMin(factors[changed], btVbd::volumeBound(old, next, tet.volume * 6, settings.volumeFloor));
		}
		for (int i = 0; i < int(x.size()); ++i)
			proposed[i] = x[i] + (proposed[i] - x[i]) * factors[i];
		btScalar common = 1;
		if (settings.workers > 1 && guardTetCount >= 4096)
		{
			guardCommonLimits.resize(guardTetCount);
			// Volume work is heavier than geometry batches; dispatch explicit independent blocks.
			btVbdParallelFor((guardTetCount + 127) / 128, settings.workers,
							 [&](int block)
							 {
								 const int end = btMin(guardTetCount, (block + 1) * 128);
								 for (int guardTet = block * 128; guardTet < end; ++guardTet)
								 {
									 const Tet &tet = tets[indexedTets ? guardTets[guardTet] : guardTet];
									 btVector3 old[4], next[4];
									 for (int j = 0; j < 4; ++j)
									 {
										 old[j] = x[tet.nodes[j]];
										 next[j] = proposed[tet.nodes[j]];
									 }
									 guardCommonLimits[guardTet] = btVbd::volumeBound(old, next, tet.volume * 6, settings.volumeFloor);
								 }
							 },
							 2, 1);
			for (btScalar limit : guardCommonLimits)
				common = btMin(common, limit);
		}
		else
		{
			for (int guardTet = 0; guardTet < guardTetCount; ++guardTet)
			{
				const Tet &tet = tets[indexedTets ? guardTets[guardTet] : guardTet];
				btVector3 old[4], next[4];
				for (int j = 0; j < 4; ++j)
				{
					old[j] = x[tet.nodes[j]];
					next[j] = proposed[tet.nodes[j]];
				}
				common = btMin(common, btVbd::volumeBound(old, next, tet.volume * 6, settings.volumeFloor));
			}
		}
		if (common < 1)
			for (int i = 0; i < int(x.size()); ++i)
				proposed[i] = x[i] + (proposed[i] - x[i]) * common;
		if (mappedSurface || !rigids.empty())
		{
			// Geometry changes between color groups; only the contact reference survives a guard pass.
			guardChangedNodes.resize(x.size());
			for (int i = 0; i < int(x.size()); ++i)
				guardChangedNodes[i] = proposed[i] != x[i];
			if (++guardStamp == 0)
			{
				for (auto &entry : guardPositions)
					entry.stamp = entry.affectedStamp = 0;
				std::fill(guardPlaneStamps.begin(), guardPlaneStamps.end(), 0);
				++guardStamp;
			}
			// One common factor preserves signed mapping weights and every accepted volume path.
			btScalar limit = 1;
			const btScalar bound = btScalar(.5 * .85) * settings.gap;
			// Mapping norm bounds certify small displacements, including signed weights.
			btScalar maximumMove2 = 0;
			for (int i = 0; i < int(proposed.size()); ++i)
				maximumMove2 = btMax(maximumMove2, (proposed[i] - reference[i]).length2());
			const bool detailed = !(mappingAmplification * btSqrt(maximumMove2) < btScalar(.99) * bound);
			if (detailed || !planes.empty() || !movingPlanes.empty())
				guardPositions.resize(surfaceVertices.size());
			// Small mapped surfaces do not amortize GPU guard launch and synchronization costs.
			if (detailed && !(settings.gpuGuards && surfaceVertices.size() >= 8192 && gpuMappedBound(proposed, bound, limit)))
			{
				// Build only when needed; mappings stay fixed between initialize calls.
				if (guardSurfaceVertices.empty())
				{
					guardSurfaceVertices.resize(x.size());
					for (int v = 0; v < int(surfaceVertices.size()); ++v)
						for (const auto &support : surfaceVertices[v].support)
							guardSurfaceVertices[support.node].push_back(v);
				}
				// Unchanged mapped vertices already satisfy the accepted bound and cannot tighten this update.
				guardAffectedVertices.clear();
				for (int i = 0; i < int(x.size()); ++i)
					if (guardChangedNodes[i])
						for (int v : guardSurfaceVertices[i])
							if (guardPositions[v].affectedStamp != guardStamp)
							{
								guardPositions[v].affectedStamp = guardStamp;
								guardAffectedVertices.push_back(v);
							}
				guardVertexLimits.resize(guardAffectedVertices.size());
				btVbdParallelGeometry(int(guardAffectedVertices.size()), settings.workers,
									  [&](int index)
									  {
										  btScalar vertexLimit = 1;
										  const auto &p = guardPosition(guardAffectedVertices[index], proposed);
										  const btVector3 start = p.current, d = p.proposed - start, r = start - p.reference;
										  if ((r + d).length2() > bound * bound && d.length2() > 0)
										  {
											  const btScalar a = d.length2(), b = r.dot(d), c = r.length2() - bound * bound;
											  vertexLimit =
												  btClamped((-b + btSqrt(btMax(btScalar(0), b * b - a * c))) / a, btScalar(0), btScalar(1));
										  }
										  guardVertexLimits[index] = vertexLimit;
									  });
				for (btScalar vertexLimit : guardVertexLimits)
					limit = btMin(limit, vertexLimit);
			}
			const bool indexedPlanes = planes.size() >= 512;
			if (indexedPlanes)
				affectedGuardPlanes();
			const int planeCount = int(indexedPlanes ? guardAffectedPlanes.size() : planes.size());
			for (int planeIndex = 0; planeIndex < planeCount; ++planeIndex)
			{
				const auto &c = planes[indexedPlanes ? guardAffectedPlanes[planeIndex] : planeIndex];
				bool changed = indexedPlanes;
				if (!changed)
					for (int v : c.nodes)
						changed = guardPosition(v, proposed).changed || changed;
				// An unchanged feature contributes exactly one to the displacement bound.
				if (!changed)
					continue;
				const btVector3 contactPoint = planePoint(c);
				btScalar gap = SIMD_INFINITY, approach = 0, scale = 1;
				for (int v : c.nodes)
				{
					const auto &cached = guardPosition(v, proposed);
					const btVector3 r = cached.reference, p = cached.proposed;
					gap = btMin(gap, c.normal.dot(r - contactPoint));
					approach = btMax(approach, -c.normal.dot(p - r));
					for (int d = 0; d < 3; ++d)
						scale = btMax(scale, btFabs(r[d]));
				}
				for (int d = 0; d < 3; ++d)
					scale = btMax(scale, btFabs(contactPoint[d]));
				const btScalar epsilon = btMax(btScalar(1e-6), btScalar(4 * 1.1920928955078125e-7) * scale);
				gap = btMax(btScalar(0), c.rigid >= 0 ? c.referenceGap : gap);
				const btScalar fraction = c.rigid >= 0						   ? btScalar(.05)
										  : gap >= 2 * epsilon && approach > 0 ? btMax(btScalar(.05), epsilon / gap)
																			   : btScalar(.5);
				const btVector3 plane = contactPoint + c.normal * (fraction * gap);
				for (int v : c.nodes)
				{
					const auto &cached = guardPosition(v, proposed);
					const btVector3 start = cached.current;
					limit = btMin(limit, btVbd::planarBound(start, cached.proposed - start, c.normal, plane, epsilon));
				}
			}
			for (const auto &plane : movingPlanes)
			{
				bool changed = false;
				for (int a : surface[plane.a])
					changed = guardPosition(a, proposed).changed || changed;
				for (int b : surface[plane.b])
					changed = guardPosition(b, proposed).changed || changed;
				if (!changed)
					continue;
				btScalar gap = SIMD_INFINITY;
				for (int a : surface[plane.a])
					for (int b : surface[plane.b])
						gap = btMin(gap, plane.normal.dot(guardPosition(a, proposed).reference - guardPosition(b, proposed).reference));
				const btScalar epsilon = btScalar(1e-6);
				const btVector3 boundary = plane.normal * (btMax(btScalar(0), gap) * btScalar(.05));
				for (int a : surface[plane.a])
					for (int b : surface[plane.b])
					{
						const btVector3 start = guardPosition(a, proposed).current - guardPosition(b, proposed).current;
						const btVector3 end = guardPosition(a, proposed).proposed - guardPosition(b, proposed).proposed;
						limit = btMin(limit, btVbd::planarBound(start, end - start, plane.normal, boundary, epsilon));
					}
			}
			for (int i = 0; i < int(x.size()); ++i)
				proposed[i] = x[i] + (proposed[i] - x[i]) * limit;
		}
	}
};
#endif
