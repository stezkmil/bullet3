// Numerical fixtures evaluated using unmodified Newton 009158e62b86 CPU kernels.
#include "BulletSoftBody/btDeformableVbdSolver.h"
#include "BulletCollision/CollisionShapes/btTriangleMesh.h"

TEST(DeformableVbd, GImpactPaddingDoesNotChangeClearance)
{
	btPrimitiveTriangle a, b;
	a.m_vertices[0] = btVector3(0, 0, 0);
	a.m_vertices[1] = btVector3(2, 0, 0);
	a.m_vertices[2] = btVector3(0, 2, 0);
	b = a;
	for (int j = 0; j < 3; ++j)
		b.m_vertices[j].setZ(1);
	a.m_margin = b.m_margin = .05;
	a.buildTriPlane();
	b.buildTriPlane();
	EXPECT_FALSE(a.overlap_test(b));
	a.m_discoveryPadding = 5;
	EXPECT_TRUE(a.overlap_test(b));
	btScalar distance2;
	btVector3 pa, pb;
	ASSERT_TRUE(a.triangle_triangle_distance(b, distance2, pa, pb));
	EXPECT_NEAR(distance2, 1, 1e-10);
	EXPECT_NEAR(a.m_margin + b.m_margin, .1, 1e-10);
	btTriangleMesh mesh;
	mesh.addTriangle(a.m_vertices[0], a.m_vertices[1], a.m_vertices[2]);
	btGImpactMeshShape shape(&mesh);
	shape.setMargin(.05);
	shape.updateBound();
	btVector3 lo, hi;
	shape.getAabb(btTransform::getIdentity(), lo, hi);
	const btVector3 original = lo;
	shape.getMeshPart(0)->setDiscoveryPadding(5);
	shape.postUpdate();
	shape.updateBound();
	shape.getAabb(btTransform::getIdentity(), lo, hi);
	EXPECT_LT(lo.z(), original.z() - 4.9);
	shape.getMeshPart(0)->setDiscoveryPadding(0);
	shape.postUpdate();
	shape.updateBound();
	shape.getAabb(btTransform::getIdentity(), lo, hi);
	EXPECT_NEAR(lo.z(), original.z(), 1e-8);
}

static void setupVbdGroundTest(btDeformableVbdSolver &s, btScalar height)
{
	s.x = {btVector3(0, 0, height), btVector3(.1, 0, height), btVector3(0, .1, height), btVector3(0, 0, height + .1)};
	s.velocity.assign(4, btVector3(0, 0, 0));
	s.external.assign(4, btVector3(0, 0, 0));
	s.mass.assign(4, .25);
	s.massDamping.assign(4, 0);
	btDeformableVbdSolver::Tet t;
	t.nodes = {{0, 1, 2, 3}};
	t.inverseRest = btMatrix3x3::getIdentity() * 10;
	t.volume = .001 / 6;
	t.mu = 1785.714285714;
	t.lambda = 7142.857142857;
	t.damping = 4.464285714;
	s.tets.push_back(t);
	btDeformableVbdSolver::Triangle ground;
	ground.x[0] = btVector3(-1, -1, 0);
	ground.x[1] = btVector3(3, -1, 0);
	ground.x[2] = btVector3(-1, 3, 0);
	ground.friction = .25;
	s.barriers.push_back(ground);
}

TEST(DeformableVbd, DiscoveryOnlyContactsExertNoForce)
{
	for (btScalar gap : {btScalar(.0005), btScalar(.005)})
	{
		btDeformableVbdSolver s;
		setupVbdGroundTest(s, .003);
		s.settings.gap = gap;
		const auto original = s.x;
		ASSERT_TRUE(s.initialize());
		ASSERT_TRUE(s.step(.002));
		EXPECT_EQ(s.contacts.empty(), gap < .003);
		for (int i = 0; i < 4; ++i)
			EXPECT_NEAR((s.x[i] - original[i]).length(), 0, 1e-12);
	}
}

TEST(DeformableVbd, RejectsPenetrationWithoutInventingRecoveryNormal)
{
	for (btScalar height : {btScalar(-.02), btScalar(0)})
	{
		btDeformableVbdSolver s;
		setupVbdGroundTest(s, height);
		const auto original = s.x;
		ASSERT_TRUE(s.initialize());
		EXPECT_FALSE(s.step(.002));
		EXPECT_TRUE(s.contacts.empty());
		for (int i = 0; i < 4; ++i)
			EXPECT_EQ(s.x[i], original[i]);
	}
}
TEST(DeformableVbd, MatchesNewtonTetBlocks)
{
	const double previous[4][3] = {{0.0, 0.0, 0.0},
								   {0.10000000149011612, 0.0, 0.0},
								   {0.02500000037252903, 0.125, 0.0},
								   {-0.0062500000931322575, 0.010416666977107525, 0.0833333358168602}};
	const double positions[4][3] = {{0.0010000000474974513, -0.003000000026077032, 0.0020000000949949026},
									{0.09400000423192978, 0.00800000037997961, -0.0010000000474974513},
									{0.029000001028180122, 0.11800000071525574, 0.008999999612569809},
									{-0.0042500002309679985, 0.01341666653752327, 0.07233333587646484}};
	const double force[4][3] = {{-3.0368757247924805, 1.2570722103118896, -6.005086898803711},
								{4.110369682312012, -2.275933265686035, 0.9823492765426636},
								{-1.452612042427063, 2.7874488830566406, -1.8292895555496216},
								{0.3791177272796631, -1.7685877084732056, 6.852027893066406}};
	const double hessian[4][3][3] = {{{437.2457275390625, 131.248291015625, 189.12075805664062},
									  {131.248291015625, 353.1363525390625, 145.3734893798828},
									  {189.12075805664062, 145.3734893798828, 473.69573974609375}},
									 {{269.7172546386719, -17.83623504638672, 14.1205472946167},
									  {-17.83623695373535, 113.11551666259766, 7.047466278076172},
									  {14.1205472946167, 7.047466278076172, 91.38557434082031}},
									 {{64.83447265625, 1.0534424781799316, 0.580785870552063},
									  {1.0534429550170898, 170.42938232421875, -10.286266326904297},
									  {0.580785870552063, -10.28626537322998, 59.10634231567383}},
									 {{141.61050415039062, 15.006302833557129, 4.401141166687012},
									  {15.006302833557129, 148.3335723876953, 10.881412506103516},
									  {4.401141166687012, 10.881412506103516, 371.5632629394531}}};
	btVector3 p[4], p0[4];
	for (int i = 0; i < 4; ++i)
	{
		p[i] = btVector3(positions[i][0], positions[i][1], positions[i][2]);
		p0[i] = btVector3(previous[i][0], previous[i][1], previous[i][2]);
	}
	const btMatrix3x3 inv(10, -2, 1, 0, 8, -1, 0, 0, 12);
	for (int i = 0; i < 4; ++i)
	{
		btVector3 f(0, 0, 0);
		btMatrix3x3 h = btMatrix3x3::getIdentity() * btScalar(0);
		btVbd::tetraBlock(p, p0, inv, 1.0 / (960 * 6), i, 1785.714285714, 7142.857142857, 4.464285714, .002, f, h);
		for (int a = 0; a < 3; ++a)
		{
			EXPECT_NEAR(f[a], force[i][a], 2e-5);
			for (int b = 0; b < 3; ++b)
				EXPECT_NEAR(h[a][b], hessian[i][a][b], 2e-4);
		}
	}
}
TEST(DeformableVbd, MatchesNewtonRegularizedAndSlidingFriction)
{
	const double force[2][3] = {{-6.124999046325684, -4.593749523162842, 70.09999084472656},
								{-13.999996185302734, -10.499998092651367, 70.09999084472656}};
	const double hessian[2][3][3] = {{{1531249.75, 0.0, 0.0}, {0.0, 1531249.75, 0.0}, {0.0, 0.0, 1005000.0}},
									 {{34999.9921875, 0.0, 0.0}, {0.0, 34999.9921875, 0.0}, {0.0, 0.0, 1005000.0}}};
	for (int i = 0; i < 2; ++i)
	{
		btVector3 f(0, 0, 0);
		btMatrix3x3 h = btMatrix3x3::getIdentity() * btScalar(0);
		btVbd::contactBlock(.00003, .0001, btVector3(0, 0, 1), i ? btVector3(.0004, .0003, -.00002) : btVector3(.000004, .000003, -.00002),
							1e6, 10, .25, .01, .002, f, h);
		for (int a = 0; a < 3; ++a)
		{
			EXPECT_NEAR(f[a], force[i][a], 2e-5);
			for (int b = 0; b < 3; ++b)
				EXPECT_NEAR(h[a][b], hessian[i][a][b], .5);
		}
	}
}
TEST(DeformableVbd, VolumeGuardDetectsCrossingEvenWhenEndpointIsPositive)
{
	// (1-2t)(1-3t) is positive at both endpoints but crosses zero twice.
	const btScalar bound = btVbd::firstVolumeBound(0, 6, -5, 1);
	EXPECT_NEAR(bound, .3, 1e-10);
}
TEST(DeformableVbd, RigidTranslationPreservesTetShape)
{
	btDeformableVbdSolver s;
	s.x = {btVector3(0, 0, 0), btVector3(.1, 0, 0), btVector3(0, .1, 0), btVector3(0, 0, .1)};
	s.velocity.assign(4, btVector3(.1, .2, .3));
	s.external.assign(4, btVector3(0, 0, -2.5));
	s.mass.assign(4, .25);
	s.massDamping.assign(4, 0);
	btDeformableVbdSolver::Tet t;
	t.nodes = {{0, 1, 2, 3}};
	t.inverseRest = btMatrix3x3::getIdentity() * 10;
	t.volume = .001 / 6;
	t.mu = 1785.714285714;
	t.lambda = 7142.857142857;
	t.damping = 4.464285714;
	s.tets.push_back(t);
	ASSERT_TRUE(s.initialize());
	EXPECT_EQ(s.colorCount, 4);
	ASSERT_TRUE(s.step(.002));
	EXPECT_NEAR((s.x[1] - s.x[0] - btVector3(.1, 0, 0)).length(), 0, 1e-12);
	EXPECT_NEAR(s.x[0].z(), .00056, 1e-12);
	EXPECT_NEAR(s.minimumJ, 1, 1e-12);
}

TEST(DeformableVbd, CubeFallsIntoGroundAndWallWithoutInversion)
{
	btDeformableVbdSolver solver;
	for (int i = 0; i < 8; ++i)
		solver.x.push_back(btVector3((i & 1) * .1, .05 + ((i >> 1) & 1) * .1, .05 + ((i >> 2) & 1) * .1));
	solver.velocity.assign(8, btVector3(0, 0, 0));
	solver.external.assign(8, btVector3(0, -1.25, -1.25));
	solver.mass.assign(8, .125);
	solver.massDamping.assign(8, .1);
	const int indices[5][4] = {{0, 1, 2, 4}, {1, 2, 3, 7}, {1, 4, 5, 7}, {2, 4, 6, 7}, {1, 2, 4, 7}};
	for (const auto &ids : indices)
	{
		btDeformableVbdSolver::Tet t;
		for (int j = 0; j < 4; ++j)
			t.nodes[j] = ids[j];
		btMatrix3x3 rest =
			btVbd::columns(solver.x[ids[1]] - solver.x[ids[0]], solver.x[ids[2]] - solver.x[ids[0]], solver.x[ids[3]] - solver.x[ids[0]]);
		if (rest.determinant() < 0)
		{
			std::swap(t.nodes[1], t.nodes[2]);
			rest = btVbd::columns(rest.getColumn(1), rest.getColumn(0), rest.getColumn(2));
		}
		t.inverseRest = rest.inverse();
		t.volume = rest.determinant() / 6;
		t.mu = 1785.714285714;
		t.lambda = 7142.857142857;
		t.damping = 4.464285714;
		solver.tets.push_back(t);
	}
	const btVector3 ground[4] = {btVector3(-1, -1, 0), btVector3(1, -1, 0), btVector3(1, 1, 0), btVector3(-1, 1, 0)};
	for (int axis = 0; axis < 2; ++axis)
		for (int k = 0; k < 2; ++k)
		{
			btDeformableVbdSolver::Triangle tri;
			const int ids[3] = {0, k + 1, k + 2};
			tri.friction = .25;
			for (int j = 0; j < 3; ++j)
			{
				tri.x[j] = ground[ids[j]];
				if (axis)
					std::swap(tri.x[j][1], tri.x[j][2]);
			}
			solver.barriers.push_back(tri);
		}
	ASSERT_TRUE(solver.initialize());
	int bothContacts = 0;
	for (int step = 0; step < 1000; ++step)
	{
		ASSERT_TRUE(solver.step(.002)) << (solver.error ? solver.error : "");
		for (const auto &p : solver.x)
		{
			ASSERT_GE(p.y(), 0);
			ASSERT_GE(p.z(), 0);
		}
		ASSERT_GE(solver.minimumJ, .0499);
		bool floor = false, wall = false;
		for (const auto &c : solver.contacts)
		{
			floor = floor || c.normal.z() > .99;
			wall = wall || c.normal.y() > .99;
		}
		if (floor && wall)
			++bothContacts;
	}
	EXPECT_GT(bothContacts, 500);
	btVector3 center(0, 0, 0);
	for (const auto &p : solver.x)
		center += p / 8;
	EXPECT_LT(center.y(), .07);
	EXPECT_LT(center.z(), .07);
}

TEST(DeformableVbd, MappedSurfaceUsesSignedWeightsAndProtectsRenderGeometry)
{
	btDeformableVbdSolver s;
	setupVbdGroundTest(s, .01005);
	s.mappedSurface = true;
	// The render face is 10 mm below the tetrahedral boundary, using extrapolation.
	for (int i = 0; i < 3; ++i)
	{
		btDeformableVbdSolver::SurfaceVertex v;
		v.support.push_back({i, btMatrix3x3::getIdentity()});
		v.support.push_back({0, btMatrix3x3::getIdentity() * btScalar(.1)});
		v.support.push_back({3, btMatrix3x3::getIdentity() * btScalar(-.1)});
		s.surfaceVertices.push_back(v);
	}
	s.surface.push_back({{0, 1, 2}});
	s.settings.gap = .0002;
	ASSERT_TRUE(s.initialize());
	ASSERT_TRUE(s.step(.0002));
	ASSERT_FALSE(s.contacts.empty());
	bool negative = false;
	for (const auto &c : s.contacts)
		for (const auto &a : c.vertex.support)
			if (a.node == 3 && a.jacobian[0][0] < 0)
				negative = true;
	EXPECT_TRUE(negative);
	for (int step = 0; step < 50; ++step)
	{
		// Force a predictor which would cross the ground without mapped DAT.
		s.velocity.assign(4, btVector3(0, 0, -1));
		ASSERT_TRUE(s.step(.0002));
		for (const auto &v : s.surfaceVertices)
			EXPECT_GT(v.position(s.x).z(), 0);
	}
}

TEST(DeformableVbd, MappedContactCombinesTwelveNodesAndTransfersVirtualWork)
{
	btDeformableVbdSolver s;
	setupVbdGroundTest(s, .03);
	const auto initial = s.x;
	for (int k = 1; k < 3; ++k)
	{
		for (auto p : initial)
			s.x.push_back(p);
		auto t = s.tets[0];
		for (int &n : t.nodes)
			n += k * 4;
		s.tets.push_back(t);
	}
	s.velocity.assign(12, btVector3(0, 0, 0));
	s.external.assign(12, btVector3(0, 0, 0));
	s.mass.assign(12, .25);
	s.massDamping.assign(12, 0);
	s.mappedSurface = true;
	const btVector3 corners[] = {btVector3(-.01, -.01, .00005), btVector3(.01, -.01, .00005), btVector3(0, .01, .00005)};
	for (int k = 0; k < 3; ++k)
	{
		btDeformableVbdSolver::SurfaceVertex v;
		for (int j = 0; j < 4; ++j)
			v.support.push_back({k * 4 + j, btMatrix3x3(0, -2, 0, 1, 0, 0, 0, 0, .5) * btScalar(.25)});
		v.offset = corners[k] - v.position(s.x);
		s.surfaceVertices.push_back(v);
	}
	s.surface.push_back({{0, 1, 2}});
	// A small rigid triangle above the interior produces an interior soft closest point.
	s.barriers[0].x[0] = btVector3(-.001, -.001, 0);
	s.barriers[0].x[1] = btVector3(.001, -.001, 0);
	s.barriers[0].x[2] = btVector3(0, .001, 0);
	ASSERT_TRUE(s.initialize());
	ASSERT_TRUE(s.step(.0002));
	ASSERT_FALSE(s.contacts.empty());
	const auto &c = s.contacts[0];
	ASSERT_EQ(c.vertex.support.size(), 12u);
	const btVector3 force(.7, -.3, 1.2);
	const btVector3 direction(.2, .4, -.1);
	const btScalar epsilon = 1e-7;
	for (const auto &a : c.vertex.support)
	{
		auto plus = s.x, minus = s.x;
		plus[a.node] += direction * epsilon;
		minus[a.node] -= direction * epsilon;
		const btScalar work = force.dot(c.vertex.position(plus) - c.vertex.position(minus)) / (2 * epsilon);
		EXPECT_NEAR(work, (a.jacobian.transpose() * force).dot(direction), 1e-8);
	}
}

TEST(DeformableVbd, InvalidMappedSurfaceDoesNotFallBackToTetBoundary)
{
	btDeformableVbdSolver s;
	setupVbdGroundTest(s, .01);
	s.mappedSurface = true;
	EXPECT_FALSE(s.initialize());
	s.surfaceVertices.resize(1);
	s.surfaceVertices[0].support.push_back({4, btMatrix3x3::getIdentity()});
	s.surface.push_back({{0, 0, 0}});
	EXPECT_FALSE(s.initialize());
	s.surfaceVertices[0].support[0].node = 0;
	s.surface[0][2] = 1;
	EXPECT_FALSE(s.initialize());
}

#include "BulletSoftBody/btDeformableMousePickingForce.h"
TEST(DeformableVbd, GrabSpringJacobianMatchesFiniteDifferenceIncludingForceCap)
{
	for (btScalar cap : {btScalar(.3), btScalar(300)})
	{
		btDeformableVbdSolver::NodeSpring spring{0, btVector3(.1, .2, .3), 100, .7, cap};
		const btVector3 p(.12, .23, .31), old(.11, .21, .32);
		btVector3 force(0, 0, 0);
		btMatrix3x3 h = btMatrix3x3::getIdentity() * btScalar(0);
		btDeformableVbdSolver::springBlock(spring, p, old, .002, force, h);
		for (int d = 0; d < 3; ++d)
		{
			btVector3 plus = p, minus = p, fp(0, 0, 0), fm(0, 0, 0);
			plus[d] += 1e-7;
			minus[d] -= 1e-7;
			btMatrix3x3 scratch = btMatrix3x3::getIdentity() * btScalar(0);
			btDeformableVbdSolver::springBlock(spring, plus, old, .002, fp, scratch);
			btDeformableVbdSolver::springBlock(spring, minus, old, .002, fm, scratch);
			for (int r = 0; r < 3; ++r)
				EXPECT_NEAR(-(fp[r] - fm[r]) / 2e-7, h[r][d], 1e-5);
		}
	}
}

TEST(DeformableVbd, GrabWorldSupportsRotationReleaseAndMillimeterUnits)
{
	for (btScalar unit : {btScalar(1), btScalar(.001)})
	{
		btSoftBodyRigidBodyCollisionConfiguration config;
		btCollisionDispatcher dispatcher(&config);
		btDbvtBroadphase broadphase;
		btDeformableBodySolver solver;
		btDeformableMultiBodyConstraintSolver constraints;
		constraints.setDeformableSolver(&solver);
		btDeformableMultiBodyDynamicsWorld world(&dispatcher, &broadphase, &constraints, &config, &solver);
		world.setVbdSolver(true, unit, 10);
		world.setGravity(btVector3(0, 0, 0));
		btVector3 p[] = {btVector3(0, 0, 0), btVector3(.1, 0, 0), btVector3(0, .1, 0), btVector3(0, 0, .1)};
		for (auto &v : p)
			v /= unit;
		const btScalar masses[] = {.25, .25, .25, .25};
		btSoftBody body(&world.getWorldInfo(), 4, p, masses);
		body.appendTetra(0, 1, 2, 3);
		body.initializeDmInverse();
		body.m_tetraScratches.resize(1);
		body.m_tetraScratchesTn.resize(1);
		body.m_cfg.collisions = 0;
		body.m_cfg.drag = 0;
		btDeformableLinearElasticityForce elastic(1000 * unit, 4000 * unit, 0, 0);
		elastic.addSoftBody(&body);
		world.addSoftBody(&body);
		world.addForce(&elastic);
		btDeformableMousePickingForce grab(500, 2, nullptr, &body.m_tetras[0], nullptr, btTransform::getIdentity(), 20 / unit);
		grab.addSoftBody(&body);
		world.addForce(&grab);
		btTransform target(btQuaternion(btVector3(0, 0, 1), .15), btVector3(.02 / unit, 0, .01 / unit));
		grab.setMouseTransform(target);
		std::vector<btDeformableMousePickingForce::NodeSpring> springs;
		grab.appendNodeSprings(springs);
		ASSERT_EQ(springs.size(), 4u);
		for (int n = 0; n < 4; ++n)
			EXPECT_NEAR((springs[n].target - target * p[n]).length() * unit, 0, 1e-10);
		for (int step = 0; step < 80; ++step)
		{
			ASSERT_EQ(world.stepSimulation(.002, 0), 1);
			ASSERT_FALSE(world.hasCoupledStepFailed());
		}
		EXPECT_GT(body.m_nodes[0].m_x.x() * unit, .003);
		EXPECT_GT(body.m_nodes[0].m_x.z() * unit, .001);
		world.removeForce(&grab);
		for (int step = 0; step < 20; ++step)
		{
			ASSERT_EQ(world.stepSimulation(.002, 0), 1);
			ASSERT_FALSE(world.hasCoupledStepFailed());
		}
		world.removeForce(&elastic);
		world.removeSoftBody(&body);
	}
}

TEST(DeformableVbd, GrabAgainstGroundRetainsContactProtection)
{
	btDeformableVbdSolver s;
	setupVbdGroundTest(s, .00015);
	s.mappedSurface = true;
	for (int i = 0; i < 4; ++i)
	{
		btDeformableVbdSolver::SurfaceVertex v;
		v.support.push_back({i, btMatrix3x3::getIdentity()});
		s.surfaceVertices.push_back(v);
	}
	s.surface = {{{0, 1, 2}}, {{0, 1, 3}}, {{0, 2, 3}}, {{1, 2, 3}}};
	for (int i = 0; i < 4; ++i)
		s.springs.push_back({i, s.x[i] + btVector3(.03, 0, -.03), 10000, 5, 30});
	ASSERT_TRUE(s.initialize());
	for (int step = 0; step < 150; ++step)
	{
		ASSERT_TRUE(s.step(.002));
		for (const auto &v : s.surfaceVertices)
			EXPECT_GT(v.position(s.x).z(), 0);
	}
	EXPECT_GT(s.x[0].x(), .001);
}

static void setupVbdTwoBodies(btDeformableVbdSolver &s, bool sameOwner = false)
{
	s.x = {btVector3(0, 0, 0),		btVector3(0, .1, 0),	  btVector3(.1, 0, 0),		btVector3(0, 0, -.1),
		   btVector3(0, 0, .00015), btVector3(.1, 0, .00015), btVector3(0, .1, .00015), btVector3(0, 0, .10015)};
	s.velocity.assign(8, btVector3(0, 0, 0));
	s.external.assign(8, btVector3(0, 0, 0));
	s.mass.assign(8, .25);
	s.massDamping.assign(8, 0);
	s.mappedSurface = true;
	for (int b = 0; b < 2; ++b)
	{
		btVector3 p[4];
		for (int j = 0; j < 4; ++j)
			p[j] = s.x[b * 4 + j];
		btMatrix3x3 rest;
		rest.setValue(p[1].x() - p[0].x(), p[2].x() - p[0].x(), p[3].x() - p[0].x(), p[1].y() - p[0].y(), p[2].y() - p[0].y(),
					  p[3].y() - p[0].y(), p[1].z() - p[0].z(), p[2].z() - p[0].z(), p[3].z() - p[0].z());
		btDeformableVbdSolver::Tet tet{{{b * 4, b * 4 + 1, b * 4 + 2, b * 4 + 3}}, rest.inverse(), rest.determinant() / 6, 1000, 4000, 0};
		s.tets.push_back(tet);
		for (int j = 0; j < 4; ++j)
		{
			btDeformableVbdSolver::SurfaceVertex v;
			v.support.push_back({b * 4 + j, btMatrix3x3::getIdentity()});
			s.surfaceVertices.push_back(v);
		}
		for (int omit = 0; omit < 4; ++omit)
		{
			std::array<int, 3> face;
			int k = 0;
			for (int j = 0; j < 4; ++j)
				if (j != omit)
					face[k++] = b * 4 + j;
			s.surface.push_back(face);
			s.surfaceOwners.push_back(sameOwner ? 0 : b);
			s.surfaceFrictions.push_back(.5);
		}
	}
	s.selfContact = {sameOwner};
}

TEST(DeformableVbd, TwoMovingSurfacesAndSelfContactUseRelativeResponse)
{
	for (bool sameOwner : {false, true})
	{
		btDeformableVbdSolver s;
		setupVbdTwoBodies(s, sameOwner);
		s.settings.iterations = 20;
		for (int i = 0; i < 8; ++i)
			s.velocity[i] = btVector3(0, 0, i < 4 ? .01 : -.01);
		ASSERT_TRUE(s.initialize());
		for (int step = 0; step < 100; ++step)
		{
			ASSERT_TRUE(s.step(.002));
			ASSERT_FALSE(s.movingPlanes.empty());
			EXPECT_GT(s.x[4].z() - s.x[0].z(), 0);
			for (const auto &c : s.contacts)
			{
				btMatrix3x3 sum = btMatrix3x3::getIdentity() * btScalar(0);
				for (const auto &a : c.vertex.support)
					sum += a.jacobian;
				for (int row = 0; row < 3; ++row)
					EXPECT_NEAR(sum[row].length(), 0, 1e-10);
			}
		}
	}
}

TEST(DeformableVbd, MovingPairFiltersAndSelfAdjacency)
{
	btDeformableVbdSolver s;
	setupVbdTwoBodies(s);
	s.collisionAllowed = {{true, false}, {false, true}};
	ASSERT_TRUE(s.initialize());
	ASSERT_TRUE(s.step(.002));
	EXPECT_TRUE(s.contacts.empty());
	btDeformableVbdSolver single;
	setupVbdGroundTest(single, .03);
	single.barriers.clear();
	single.selfContact = {true};
	ASSERT_TRUE(single.initialize());
	ASSERT_TRUE(single.step(.002));
	EXPECT_TRUE(single.contacts.empty());
}

TEST(DeformableVbd, WorldStepsMultipleBodiesAndRejectsTogether)
{
	btDeformableVbdSolver source;
	setupVbdTwoBodies(source);
	btSoftBodyRigidBodyCollisionConfiguration config;
	btCollisionDispatcher dispatcher(&config);
	btDbvtBroadphase broadphase;
	btDeformableBodySolver solver;
	btDeformableMultiBodyConstraintSolver constraints;
	constraints.setDeformableSolver(&solver);
	btDeformableMultiBodyDynamicsWorld world(&dispatcher, &broadphase, &constraints, &config, &solver);
	world.setVbdSolver(true);
	world.setGravity(btVector3(0, 0, 0));
	const btScalar masses[] = {.25, .25, .25, .25};
	btSoftBody a(&world.getWorldInfo(), 4, source.x.data(), masses), b(&world.getWorldInfo(), 4, source.x.data() + 4, masses);
	btDeformableLinearElasticityForce material(1000, 4000, 0, 0);
	for (auto *body : {&a, &b})
	{
		body->appendTetra(0, 1, 2, 3);
		body->initializeDmInverse();
		body->m_tetraScratches.resize(1);
		body->m_tetraScratchesTn.resize(1);
		body->m_cfg.drag = 0;
		world.addSoftBody(body);
		material.addSoftBody(body);
	}
	world.addForce(&material);
	for (int i = 0; i < 4; ++i)
	{
		a.m_nodes[i].m_v = btVector3(0, 0, .01);
		b.m_nodes[i].m_v = btVector3(0, 0, -.01);
	}
	for (int step = 0; step < 100; ++step)
	{
		ASSERT_EQ(world.stepSimulation(.002, 0), 1);
		ASSERT_FALSE(world.hasCoupledStepFailed());
	}
	for (int i = 0; i < 4; ++i)
		b.m_nodes[i].m_x.setZ(b.m_nodes[i].m_x.z() - .005);
	const btVector3 oldA = a.m_nodes[0].m_x, oldB = b.m_nodes[0].m_x;
	EXPECT_EQ(world.stepSimulation(.002, 0), 0);
	EXPECT_TRUE(world.hasCoupledStepFailed());
	EXPECT_EQ(a.m_nodes[0].m_x, oldA);
	EXPECT_EQ(b.m_nodes[0].m_x, oldB);
	world.removeForce(&material);
	world.removeSoftBody(&a);
	world.removeSoftBody(&b);
}

#include "BulletCollision/CollisionShapes/btBoxShape.h"
#include "BulletCollision/CollisionShapes/btSphereShape.h"
#include "BulletCollision/CollisionShapes/btCapsuleShape.h"
TEST(DeformableVbd, ConvexBoxesSpheresAndCapsulesUseNativeShapeDistance)
{
	btBoxShape box(btVector3(200, 200, 100));
	btSphereShape sphere(100);
	btCapsuleShape capsule(100, 200);
	for (const btConvexShape *shape : {static_cast<const btConvexShape *>(&box), static_cast<const btConvexShape *>(&sphere),
									   static_cast<const btConvexShape *>(&capsule)})
	{
		btDeformableVbdSolver s;
		setupVbdGroundTest(s, .00015);
		s.barriers.clear();
		btTransform transform;
		transform.setIdentity();
		transform.setOrigin(btVector3(0, 0, -100));
		s.convexBarriers.push_back({shape, transform, .5, -1});
		s.external.assign(4, btVector3(0, 0, -.25));
		ASSERT_TRUE(s.initialize());
		for (int step = 0; step < 50; ++step)
		{
			ASSERT_TRUE(s.step(.002));
			ASSERT_FALSE(s.planes.empty());
		}
	}
}

TEST(DeformableVbd, ConvexIntersectionRejectsBeforeChangingState)
{
	btSphereShape sphere(100);
	btTransform transform;
	transform.setIdentity();
	transform.setOrigin(btVector3(0, 0, -100));
	btDeformableVbdSolver s;
	setupVbdGroundTest(s, -.01);
	s.barriers.clear();
	s.convexBarriers.push_back({&sphere, transform, .5, -1});
	ASSERT_TRUE(s.initialize());
	const auto original = s.x;
	EXPECT_FALSE(s.step(.002));
	for (int n = 0; n < 4; ++n)
		EXPECT_EQ(s.x[n], original[n]);
}

TEST(DeformableVbd, ParallelVertexBlocksMatchSerialAcrossWorkerCounts)
{
	btDeformableVbdSolver serial;
	setupVbdGroundTest(serial, 1);
	serial.barriers.clear();
	const auto points = serial.x;
	const auto tet = serial.tets[0];
	for (int body = 1; body < 100; ++body)
	{
		for (const auto &p : points)
			serial.x.push_back(p + btVector3(body, 0, 0));
		auto next = tet;
		for (int &node : next.nodes)
			node += 4 * body;
		serial.tets.push_back(next);
	}
	serial.velocity.assign(400, btVector3(0, 0, 0));
	serial.mass.assign(400, .25);
	serial.massDamping.assign(400, 0);
	serial.external.assign(400, btVector3(0, 0, -2.5));
	ASSERT_TRUE(serial.initialize());
	for (int workers : {2, 4, int(btMin(256u, btMax(2u, std::thread::hardware_concurrency() - 1)))})
	{
		btDeformableVbdSolver parallel;
		parallel.x = serial.x;
		parallel.velocity = serial.velocity;
		parallel.mass = serial.mass;
		parallel.massDamping = serial.massDamping;
		parallel.external = serial.external;
		parallel.tets = serial.tets;
		parallel.settings.workers = workers;
		ASSERT_TRUE(parallel.initialize());
		const auto start = serial.x, velocity = serial.velocity;
		for (int step = 0; step < 10; ++step)
		{
			ASSERT_TRUE(serial.step(.002));
			ASSERT_TRUE(parallel.step(.002));
		}
		for (int i = 0; i < 400; ++i)
		{
			EXPECT_EQ(serial.x[i], parallel.x[i]);
			EXPECT_EQ(serial.velocity[i], parallel.velocity[i]);
		}
		serial.x = start;
		serial.velocity = velocity;
	}
}

TEST(DeformableVbd, DynamicRigidReceivesContactResponseWithinSolve)
{
	btDeformableVbdSolver s;
	setupVbdGroundTest(s, .00015);
	s.barriers.clear();
	s.settings.iterations = 40;
	btBoxShape box(btVector3(200, 200, 100));
	btTransform pose;
	pose.setIdentity();
	pose.setOrigin(btVector3(0, 0, -.1));
	btDeformableVbdSolver::Rigid rigid;
	rigid.pose = pose;
	rigid.mass = 1;
	rigid.inertia = btVector3(.02, .02, .03);
	rigid.radius = .3;
	s.rigids.push_back(rigid);
	btTransform collisionPose = pose;
	collisionPose.setOrigin(pose.getOrigin() * 1000);
	s.convexBarriers.push_back({&box, collisionPose, .5, -1, 0, btTransform::getIdentity()});
	s.velocity.assign(4, btVector3(0, 0, -.1));
	ASSERT_TRUE(s.initialize());
	for (int step = 0; step < 100; ++step)
	{
		ASSERT_TRUE(s.step(.002));
	}
	EXPECT_LT(s.rigids[0].pose.getOrigin().z(), -.1001);
	EXPECT_LT(s.rigids[0].velocity.z(), -.001);
	EXPECT_GT(s.velocity[0].z(), -.1);
	btVector3 momentum = s.rigids[0].velocity;
	for (int i = 0; i < 4; ++i)
		momentum += s.velocity[i] * s.mass[i];
	EXPECT_NEAR(momentum.z(), -.1, .02);
}

TEST(DeformableVbd, AttachmentPullsBothSoftAndRigidAndTransfersTorque)
{
	btDeformableVbdSolver s;
	setupVbdGroundTest(s, 1);
	s.barriers.clear();
	s.settings.iterations = 40;
	btDeformableVbdSolver::Rigid rigid;
	rigid.pose.setIdentity();
	rigid.pose.setOrigin(btVector3(0, 0, 1));
	rigid.mass = .5;
	rigid.inertia = btVector3(.01, .01, .01);
	rigid.radius = .2;
	s.rigids.push_back(rigid);
	s.attachments.push_back({1, 0, btVector3(.1, 0, 0), 10000, 10});
	s.external[1] = btVector3(0, 2, 0);
	ASSERT_TRUE(s.initialize());
	for (int step = 0; step < 100; ++step)
		ASSERT_TRUE(s.step(.002));
	EXPECT_GT(s.rigids[0].pose.getOrigin().y(), 1e-5);
	EXPECT_GT(btDeformableVbdSolver::rotationVector(s.rigids[0].pose.getRotation()).z(), 1e-5);
	EXPECT_LT((s.x[1] - s.rigids[0].pose * btVector3(.1, 0, 0)).length(), .001);
}

TEST(DeformableVbd, WorldAcceptsDynamicRigidAndDeformableAnchor)
{
	btSoftBodyRigidBodyCollisionConfiguration config;
	btCollisionDispatcher dispatcher(&config);
	btDbvtBroadphase broadphase;
	btDeformableBodySolver solver;
	btDeformableMultiBodyConstraintSolver constraints;
	constraints.setDeformableSolver(&solver);
	btDeformableMultiBodyDynamicsWorld world(&dispatcher, &broadphase, &constraints, &config, &solver);
	world.setVbdSolver(true);
	world.setGravity(btVector3(0, 0, 0));
	const btVector3 p[] = {btVector3(0, 0, 0), btVector3(.1, 0, 0), btVector3(0, .1, 0), btVector3(0, 0, .1)};
	const btScalar masses[] = {.25, .25, .25, .25};
	btSoftBody body(&world.getWorldInfo(), 4, p, masses);
	body.appendTetra(0, 1, 2, 3);
	body.initializeDmInverse();
	body.m_tetraScratches.resize(1);
	body.m_tetraScratchesTn.resize(1);
	body.m_cfg.drag = 0;
	btDeformableLinearElasticityForce material(1000, 4000, 0, 0);
	material.addSoftBody(&body);
	world.addSoftBody(&body);
	world.addForce(&material);
	btBoxShape box(btVector3(.02, .02, .02));
	btVector3 inertia;
	box.calculateLocalInertia(.5, inertia);
	btRigidBody::btRigidBodyConstructionInfo info(.5, nullptr, &box, inertia);
	btRigidBody rigid(info);
	btTransform pose;
	pose.setIdentity();
	pose.setOrigin(btVector3(1, 0, 0));
	rigid.setWorldTransform(pose);
	world.addRigidBody(&rigid);
	body.appendDeformableAnchor(0, &rigid, 1);
	for (int step = 0; step < 100; ++step)
	{
		rigid.applyCentralForce(btVector3(0, 1, 0));
		ASSERT_EQ(world.stepSimulation(.002, 0), 1);
		ASSERT_FALSE(world.hasCoupledStepFailed());
	}
	EXPECT_GT(rigid.getWorldTransform().getOrigin().y(), 0);
	EXPECT_GT(body.m_nodes[0].m_x.y(), 0);
	EXPECT_LT((body.m_nodes[0].m_x - rigid.getWorldTransform() * body.m_deformableAnchors[0].m_local).length(), .002);
	body.removeDeformableAnchor(0);
	world.removeForce(&material);
	world.removeSoftBody(&body);
	world.removeRigidBody(&rigid);
}

TEST(DeformableVbd, DynamicRigidFallsOntoConvexGround)
{
	btDeformableVbdSolver s;
	setupVbdGroundTest(s, 2);
	s.barriers.clear();
	s.settings.iterations = 30;
	btBoxShape ground(btVector3(1000, 1000, 100));
	btSphereShape ball(50);
	btTransform groundPose = btTransform::getIdentity();
	groundPose.setOrigin(btVector3(0, 0, -100));
	s.convexBarriers.push_back({&ground, groundPose, .5, 0});
	btDeformableVbdSolver::Rigid r;
	r.pose.setOrigin(btVector3(0, 0, .052));
	r.mass = 1;
	r.inertia = btVector3(.001, .001, .001);
	r.radius = .05;
	r.force = btVector3(0, 0, -9.81);
	s.rigids.push_back(r);
	btTransform ballPose = r.pose;
	ballPose.setOrigin(ballPose.getOrigin() * 1000);
	s.convexBarriers.push_back({&ball, ballPose, .5, 1, 0, btTransform::getIdentity()});
	ASSERT_TRUE(s.initialize());
	for (int step = 0; step < 400; ++step)
	{
		ASSERT_TRUE(s.step(.002));
		EXPECT_GT(s.rigids[0].pose.getOrigin().z(), .05);
	}
	EXPECT_LT(s.rigids[0].pose.getOrigin().z(), .0502);
	EXPECT_NEAR(s.rigids[0].velocity.z(), 0, .002);
}

#include "BulletDynamics/ConstraintSolver/btPoint2PointConstraint.h"
#include "BulletDynamics/ConstraintSolver/btHingeConstraint.h"
TEST(DeformableVbd, WorldPreservesPointJointInMetresAndMillimetres)
{
	for (btScalar unit : {btScalar(1), btScalar(.001)})
	{
		btSoftBodyRigidBodyCollisionConfiguration config;
		btCollisionDispatcher dispatcher(&config);
		btDbvtBroadphase broadphase;
		btDeformableBodySolver solver;
		btDeformableMultiBodyConstraintSolver constraints;
		constraints.setDeformableSolver(&solver);
		btDeformableMultiBodyDynamicsWorld world(&dispatcher, &broadphase, &constraints, &config, &solver);
		world.setVbdSolver(true, unit, 30);
		world.setGravity(btVector3(0, 0, -9.81 / unit));
		const btVector3 p[] = {btVector3(5, 0, 0) / unit, btVector3(5.1, 0, 0) / unit, btVector3(5, .1, 0) / unit,
							   btVector3(5, 0, .1) / unit};
		const btScalar masses[] = {.25, .25, .25, .25};
		btSoftBody body(&world.getWorldInfo(), 4, p, masses);
		body.appendTetra(0, 1, 2, 3);
		body.initializeDmInverse();
		body.m_tetraScratches.resize(1);
		body.m_tetraScratchesTn.resize(1);
		body.m_cfg.drag = 0;
		btDeformableLinearElasticityForce material(1000, 4000, 0, 0);
		material.addSoftBody(&body);
		world.addSoftBody(&body);
		world.addForce(&material);
		btSphereShape shape(.05 / unit);
		btVector3 inertia;
		shape.calculateLocalInertia(1, inertia);
		btRigidBody::btRigidBodyConstructionInfo info(1, nullptr, &shape, inertia);
		btRigidBody rigid(info);
		world.addRigidBody(&rigid);
		btPoint2PointConstraint joint(rigid, btVector3(0, 0, .1 / unit));
		btJointFeedback feedback{};
		joint.setJointFeedback(&feedback);
		world.addConstraint(&joint);
		for (int step = 0; step < 100; ++step)
		{
			ASSERT_EQ(world.stepSimulation(.002, 0), 1);
			ASSERT_FALSE(world.hasCoupledStepFailed());
		}
		EXPECT_LT(rigid.getWorldTransform().getOrigin().length() * unit, 1e-5);
		EXPECT_GT(joint.getAppliedImpulse(), 0);
		EXPECT_NEAR(feedback.m_appliedForceBodyA.z() * unit, 9.81, .01);
		world.removeConstraint(&joint);
		world.removeRigidBody(&rigid);
		world.removeForce(&material);
		world.removeSoftBody(&body);
	}
}
TEST(DeformableVbd, WorldHingeMotorUsesAngularImpulseUnits)
{
	for (btScalar unit : {btScalar(1), btScalar(.001)})
	{
		btSoftBodyRigidBodyCollisionConfiguration config;
		btCollisionDispatcher dispatcher(&config);
		btDbvtBroadphase broadphase;
		btDeformableBodySolver solver;
		btDeformableMultiBodyConstraintSolver constraints;
		constraints.setDeformableSolver(&solver);
		btDeformableMultiBodyDynamicsWorld world(&dispatcher, &broadphase, &constraints, &config, &solver);
		world.setVbdSolver(true, unit, 30);
		world.setGravity(btVector3(0, 0, -9.81 / unit));
		const btVector3 p[] = {btVector3(5, 0, 0) / unit, btVector3(5.1, 0, 0) / unit, btVector3(5, .1, 0) / unit,
							   btVector3(5, 0, .1) / unit};
		const btScalar masses[] = {.25, .25, .25, .25};
		btSoftBody body(&world.getWorldInfo(), 4, p, masses);
		body.appendTetra(0, 1, 2, 3);
		body.initializeDmInverse();
		body.m_tetraScratches.resize(1);
		body.m_tetraScratchesTn.resize(1);
		body.m_cfg.drag = 0;
		btDeformableLinearElasticityForce material(1000, 4000, 0, 0);
		material.addSoftBody(&body);
		world.addSoftBody(&body);
		world.addForce(&material);
		btSphereShape shape(.05 / unit);
		btVector3 inertia;
		shape.calculateLocalInertia(1, inertia);
		btRigidBody::btRigidBodyConstructionInfo info(1, nullptr, &shape, inertia);
		btRigidBody rigid(info);
		world.addRigidBody(&rigid);
		btHingeConstraint joint(rigid, btVector3(0, 0, 0), btVector3(0, 0, 1));
		joint.enableAngularMotor(true, 1, .01 / (unit * unit));
		world.addConstraint(&joint);
		for (int step = 0; step < 100; ++step)
		{
			ASSERT_EQ(world.stepSimulation(.002, 0), 1);
			ASSERT_FALSE(world.hasCoupledStepFailed());
		}
		EXPECT_LT(rigid.getWorldTransform().getOrigin().length() * unit, 1e-5);
		EXPECT_NEAR(btFabs(rigid.getAngularVelocity().z()), 1, .01);
		EXPECT_GT(btFabs(btDeformableVbdSolver::rotationVector(rigid.getWorldTransform().getRotation()).z()), .1);
		world.removeConstraint(&joint);
		world.removeRigidBody(&rigid);
		world.removeForce(&material);
		world.removeSoftBody(&body);
	}
}

TEST(DeformableVbd, RigidMotionFactorsRetainLockedAxes)
{
	btDeformableVbdSolver s;
	setupVbdGroundTest(s, 2);
	s.barriers.clear();
	btDeformableVbdSolver::Rigid rigid;
	rigid.mass = 1;
	rigid.linearFactor = btVector3(0, 1, 1);
	rigid.angularFactor = btVector3(1, 1, 0);
	rigid.force = btVector3(10, 10, 0);
	rigid.torque = btVector3(0, 0, 10);
	s.rigids.push_back(rigid);
	ASSERT_TRUE(s.initialize());
	for (int step = 0; step < 100; ++step)
		ASSERT_TRUE(s.step(.002));
	EXPECT_EQ(s.rigids[0].pose.getOrigin().x(), 0);
	EXPECT_GT(s.rigids[0].pose.getOrigin().y(), 0);
	EXPECT_EQ(btDeformableVbdSolver::rotationVector(s.rigids[0].pose.getRotation()).z(), 0);
}

TEST(DeformableVbd, ShallowConvexOverlapRecoversWithoutLaunchingBody)
{
	btBoxShape box(btVector3(200, 200, 100));
	btTransform pose = btTransform::getIdentity();
	pose.setOrigin(btVector3(0, 0, -100));
	for (bool recover : {false, true})
	{
		btDeformableVbdSolver s;
		setupVbdGroundTest(s, -.00002);
		s.barriers.clear();
		s.convexBarriers.push_back({&box, pose, .5, -1});
		s.settings.recoveryDistance = recover ? .0001 : 0;
		ASSERT_TRUE(s.initialize());
		const auto original = s.x;
		if (!recover)
		{
			EXPECT_FALSE(s.step(.002));
			EXPECT_EQ(s.x, original);
			continue;
		}
		// A short physical step isolates recovery displacement from subsequent contact acceleration.
		ASSERT_TRUE(s.step(1e-8));
		EXPECT_EQ(s.recoveredIntersections, 1);
		for (int n = 0; n < 4; ++n)
		{
			EXPECT_GT(s.x[n].z(), 0);
			EXPECT_LT((s.x[n] - original[n]).length(), s.settings.recoveryDistance);
			EXPECT_LT(s.velocity[n].length(), .001);
		}
		for (int step = 0; step < 100; ++step)
			ASSERT_TRUE(s.step(.002));
	}
}

TEST(DeformableVbd, MeshOverlapUsesRevalidatedAcceptedPositions)
{
	for (bool history : {false, true})
	{
		btDeformableVbdSolver s;
		setupVbdGroundTest(s, .00002);
		s.external.assign(4, btVector3(0, 0, 0));
		ASSERT_TRUE(s.initialize());
		if (history)
			s.recoveryPositions = s.x;
		for (int n = 0; n < 3; ++n)
			s.x[n].setZ(s.x[n].z() - .00004);
		const auto intersecting = s.x;
		if (!history)
		{
			EXPECT_FALSE(s.step(.002));
			EXPECT_EQ(s.x, intersecting);
			continue;
		}
		ASSERT_TRUE(s.step(.002));
		EXPECT_EQ(s.recoveredIntersections, 1);
		for (const auto &p : s.x)
			EXPECT_GT(p.z(), 0);
	}
}

#include "BulletCollision/CollisionShapes/btCompoundShape.h"
TEST(DeformableVbd, WorldCompoundChildrenRespectLocalTransforms)
{
	btDeformableVbdSolver source;
	setupVbdGroundTest(source, .00015);
	btSoftBodyRigidBodyCollisionConfiguration config;
	btCollisionDispatcher dispatcher(&config);
	btDbvtBroadphase broadphase;
	btDeformableBodySolver solver;
	btDeformableMultiBodyConstraintSolver constraints;
	constraints.setDeformableSolver(&solver);
	btDeformableMultiBodyDynamicsWorld world(&dispatcher, &broadphase, &constraints, &config, &solver);
	world.setVbdSolver(true);
	world.setGravity(btVector3(0, 0, 0));
	const btScalar masses[] = {.25, .25, .25, .25};
	btSoftBody body(&world.getWorldInfo(), 4, source.x.data(), masses);
	body.appendTetra(0, 1, 2, 3);
	body.initializeDmInverse();
	body.m_tetraScratches.resize(1);
	body.m_tetraScratchesTn.resize(1);
	body.m_cfg.drag = 0;
	btDeformableLinearElasticityForce material(1000, 4000, 0, 0);
	material.addSoftBody(&body);
	world.addSoftBody(&body);
	world.addForce(&material);
	btBoxShape tile(btVector3(.2, .2, .1));
	tile.setMargin(.001);
	btCompoundShape compound;
	for (btScalar dx : {btScalar(-.2), btScalar(.2)})
	{
		btTransform child = btTransform::getIdentity();
		child.setOrigin(btVector3(dx, 0, -.1));
		compound.addChildShape(child, &tile);
	}
	btRigidBody::btRigidBodyConstructionInfo info(0, nullptr, &compound);
	btRigidBody ground(info);
	world.addRigidBody(&ground);
	for (int step = 0; step < 100; ++step)
	{
		for (int n = 0; n < 4; ++n)
			body.m_nodes[n].m_f = btVector3(0, 0, -.25);
		ASSERT_EQ(world.stepSimulation(.002, 0), 1);
		ASSERT_FALSE(world.hasCoupledStepFailed());
		for (int n = 0; n < 4; ++n)
			EXPECT_GT(body.m_nodes[n].m_x.z(), 0);
	}
	world.removeRigidBody(&ground);
	world.removeForce(&material);
	world.removeSoftBody(&body);
}

TEST(DeformableVbd, RigidOnlySceneRetainsOriginalBulletPath)
{
	btSoftBodyRigidBodyCollisionConfiguration config;
	btCollisionDispatcher dispatcher(&config);
	btDbvtBroadphase broadphase;
	btDeformableBodySolver solver;
	btDeformableMultiBodyConstraintSolver constraints;
	constraints.setDeformableSolver(&solver);
	btDeformableMultiBodyDynamicsWorld world(&dispatcher, &broadphase, &constraints, &config, &solver);
	world.setVbdSolver(true);
	world.setGravity(btVector3(0, 0, -9.81));
	btSphereShape shape(.05);
	btVector3 inertia;
	shape.calculateLocalInertia(1, inertia);
	btRigidBody::btRigidBodyConstructionInfo info(1, nullptr, &shape, inertia);
	btRigidBody rigid(info);
	world.addRigidBody(&rigid);
	for (int step = 0; step < 50; ++step)
	{
		ASSERT_EQ(world.stepSimulation(.002, 0), 1);
		ASSERT_FALSE(world.hasCoupledStepFailed());
	}
	EXPECT_NEAR(rigid.getLinearVelocity().z(), -.981, .001);
	EXPECT_NEAR(rigid.getWorldTransform().getOrigin().z(), -.050031, .001);
	world.removeRigidBody(&rigid);
}

TEST(DeformableVbd, InvalidSettingsRejectWithoutMutation)
{
	btDeformableVbdSolver s;
	setupVbdGroundTest(s, 1);
	ASSERT_TRUE(s.initialize());
	const auto initial = s.x;
	s.settings.workers = 10000;
	EXPECT_FALSE(s.step(.002));
	EXPECT_EQ(s.x, initial);
	s.settings.workers = 1;
	s.settings.ke = SIMD_INFINITY * SIMD_INFINITY;
	EXPECT_FALSE(s.step(.002));
	EXPECT_EQ(s.x, initial);
}

TEST(DeformableVbd, WorldHingeLimitsBoundMotorMotion)
{
	for (btScalar target : {btScalar(-1), btScalar(1)})
		for (btScalar unit : {btScalar(1), btScalar(.001)})
		{
			btSoftBodyRigidBodyCollisionConfiguration config;
			btCollisionDispatcher dispatcher(&config);
			btDbvtBroadphase broadphase;
			btDeformableBodySolver solver;
			btDeformableMultiBodyConstraintSolver constraints;
			constraints.setDeformableSolver(&solver);
			btDeformableMultiBodyDynamicsWorld world(&dispatcher, &broadphase, &constraints, &config, &solver);
			world.setVbdSolver(true, unit, 30);
			world.setGravity(btVector3(0, 0, -9.81 / unit));
			const btVector3 p[] = {btVector3(5, 0, 0) / unit, btVector3(5.1, 0, 0) / unit, btVector3(5, .1, 0) / unit,
								   btVector3(5, 0, .1) / unit};
			const btScalar masses[] = {.25, .25, .25, .25};
			btSoftBody body(&world.getWorldInfo(), 4, p, masses);
			body.appendTetra(0, 1, 2, 3);
			body.initializeDmInverse();
			body.m_tetraScratches.resize(1);
			body.m_tetraScratchesTn.resize(1);
			body.m_cfg.drag = 0;
			btDeformableLinearElasticityForce material(1000, 4000, 0, 0);
			material.addSoftBody(&body);
			world.addSoftBody(&body);
			world.addForce(&material);
			btSphereShape shape(.05 / unit);
			btVector3 inertia;
			shape.calculateLocalInertia(1, inertia);
			btRigidBody::btRigidBodyConstructionInfo info(1, nullptr, &shape, inertia);
			btRigidBody rigid(info);
			world.addRigidBody(&rigid);
			btHingeConstraint joint(rigid, btVector3(0, 0, 0), btVector3(0, 0, 1));
			joint.enableAngularMotor(true, target, .01 / (unit * unit));
			joint.setLimit(-.05, .05);
			world.addConstraint(&joint);
			for (int step = 0; step < 100; ++step)
			{
				ASSERT_EQ(world.stepSimulation(.002, 0), 1);
				ASSERT_FALSE(world.hasCoupledStepFailed());
			}
			EXPECT_LT(rigid.getWorldTransform().getOrigin().length() * unit, 1e-5);
			EXPECT_NEAR(rigid.getAngularVelocity().z(), 0, .02);
			EXPECT_LE(btFabs(btDeformableVbdSolver::rotationVector(rigid.getWorldTransform().getRotation()).z()), .055);
			world.removeConstraint(&joint);
			world.removeRigidBody(&rigid);
			world.removeForce(&material);
			world.removeSoftBody(&body);
		}
}

TEST(DeformableVbd, CachedCollisionHierarchyMatchesFreshGeometry)
{
	btDeformableVbdCollisionCache cache;
	for (int pass = 0; pass < 8; ++pass)
	{
		btDeformableVbdSolver cached(&cache), fresh;
		for (auto *solver : {&cached, &fresh})
		{
			setupVbdGroundTest(*solver, .003);
			solver->mappedSurface = true;
			for (int i = 0; i < 4; ++i)
			{
				btDeformableVbdSolver::SurfaceVertex vertex;
				vertex.support.push_back({i, btMatrix3x3::getIdentity()});
				solver->surfaceVertices.push_back(vertex);
			}
			solver->surface = {{{0, 2, 1}}, {{0, 1, 3}}, {{0, 3, 2}}, {{1, 2, 3}}};
			if (pass == 2 || pass == 5)
				solver->barriers.clear();
			else
			{
				for (auto &point : solver->barriers[0].x)
					point.setZ(pass % 2 ? -.1 : 0);
				if (pass > 3)
					solver->barriers.push_back(solver->barriers[0]);
			}
			for (auto &velocity : solver->velocity)
				velocity.setX(pass % 2 ? .01 : 3);
			ASSERT_TRUE(solver->initialize());
			ASSERT_TRUE(solver->step(.002)) << solver->error;
		}
		EXPECT_EQ(cached.contacts.size(), fresh.contacts.size());
		for (int i = 0; i < 4; ++i)
			EXPECT_EQ(cached.x[i], fresh.x[i]);
	}
}

TEST(DeformableVbd, MovingPairPruningPreservesTriangleOrderFilter)
{
	btDeformableVbdSolver s;
	setupVbdTwoBodies(s);
	s.collisionAllowed = {{false, true}, {false, false}};
	ASSERT_TRUE(s.initialize());
	ASSERT_TRUE(s.step(.002));
	EXPECT_FALSE(s.movingPlanes.empty());
}

TEST(DeformableVbd, NativeBarrierMatchesFlattenedMeshAtCandidatePositions)
{
	btTriangleMesh geometry;
	geometry.addTriangle(btVector3(-1000, -1000, 0), btVector3(3000, -1000, 0), btVector3(-1000, 3000, 0));
	for (int i = 0; i < 100; ++i)
	{
		const btScalar farCoordinate = 10000 + i * 100;
		geometry.addTriangle(btVector3(farCoordinate, 0, 0), btVector3(farCoordinate + 50, 0, 0), btVector3(farCoordinate, 50, 0));
	}
	btGImpactMeshShape shape(&geometry);
	shape.setLocalScaling(btVector3(2, 1, 1));
	shape.updateBound();
	btTransform transform(btQuaternion(0, 0, 1, 0), btVector3(30, 20, 0));
	btDeformableVbdSolver native, flat;
	setupVbdGroundTest(native, .003);
	setupVbdGroundTest(flat, .003);
	// Keep the fixture geometry inside its explicit discovery range.
	native.settings.gap = flat.settings.gap = .005;
	native.barriers.clear();
	flat.barriers.clear();
	native.settings.collisionUnitsPerMeter = flat.settings.collisionUnitsPerMeter = 1024;
	native.nativeBarriers.push_back({shape.getMeshPart(0), transform, .25, 1, {}});
	shape.getMeshPart(0)->lockChildShapes();
	for (int i = 0; i < 101; ++i)
	{
		btPrimitiveTriangle triangle;
		shape.getMeshPart(0)->getPrimitiveManager()->get_primitive_triangle(i, triangle, false);
		btDeformableVbdSolver::Triangle output;
		for (int j = 0; j < 3; ++j)
			output.x[j] = (transform * triangle.m_vertices[j]) / 1024;
		output.friction = .25;
		output.owner = 1;
		flat.barriers.push_back(output);
	}
	shape.getMeshPart(0)->unlockChildShapes();
	ASSERT_TRUE(native.initialize());
	ASSERT_TRUE(flat.initialize());
	for (int step = 0; step < 30; ++step)
	{
		native.velocity.assign(4, btVector3(.03, 0, -.1));
		flat.velocity = native.velocity;
		ASSERT_TRUE(native.step(.002));
		ASSERT_TRUE(flat.step(.002));
		ASSERT_EQ(native.barriers.size(), 1u);
		for (int j = 0; j < 3; ++j)
			EXPECT_EQ(native.barriers[0].x[j], flat.barriers[0].x[j]);
		for (int i = 0; i < 4; ++i)
			EXPECT_NEAR((native.x[i] - flat.x[i]).length(), 0, 1e-8) << "step " << step;
	}
	ASSERT_EQ(native.barriers.size(), 1u);
	for (int j = 0; j < 3; ++j)
		EXPECT_EQ(native.barriers[0].x[j], flat.barriers[0].x[j]);
}
TEST(DeformableVbd, NativeBarrierHonorsFiltersBeforeExtractingTriangles)
{
	btTriangleMesh geometry;
	geometry.addTriangle(btVector3(-1000, -1000, 0), btVector3(3000, -1000, 0), btVector3(-1000, 3000, 0));
	btGImpactMeshShape shape(&geometry);
	shape.updateBound();
	btDeformableVbdSolver s;
	setupVbdGroundTest(s, .003);
	// Keep the fixture geometry inside its explicit discovery range.
	s.settings.gap = .005;
	s.barriers.clear();
	s.nativeBarriers.push_back({shape.getMeshPart(0), btTransform::getIdentity(), .25, 1, {}});
	s.collisionAllowed = {{true, false}, {false, true}};
	ASSERT_TRUE(s.initialize());
	ASSERT_TRUE(s.step(.002));
	EXPECT_TRUE(s.barriers.empty());
	EXPECT_TRUE(s.contacts.empty());
}

TEST(DeformableVbd, WorldMappingCacheChecksInPlaceEditsAndTransformChanges)
{
	btSoftBodyRigidBodyCollisionConfiguration config;
	btCollisionDispatcher dispatcher(&config);
	btDbvtBroadphase broadphase;
	btDeformableBodySolver solver;
	btDeformableMultiBodyConstraintSolver constraints;
	constraints.setDeformableSolver(&solver);
	btDeformableMultiBodyDynamicsWorld world(&dispatcher, &broadphase, &constraints, &config, &solver);
	world.setVbdSolver(true);
	world.setGravity(btVector3(0, 0, 0));
	btVector3 positions[] = {btVector3(0, 0, 0), btVector3(.1, 0, 0), btVector3(0, .1, 0), btVector3(0, 0, .1)};
	btScalar masses[] = {.25, .25, .25, .25};
	struct MappedBody : btSoftBody
	{
		using btSoftBody::btSoftBody;
		std::vector<btVertexToTetraMapping> mapping;
		const std::vector<btVertexToTetraMapping> *getCollisionShapeVertexToSimTetra() const override
		{
			return &mapping;
		}
	} body(&world.getWorldInfo(), 4, positions, masses);
	body.appendTetra(0, 1, 2, 3);
	body.initializeDmInverse();
	body.m_tetraScratches.resize(1);
	body.m_tetraScratchesTn.resize(1);
	for (int i = 0; i < 4; ++i)
		body.m_nodes[i].m_frozen = 1;
	btTriangleMesh geometry;
	geometry.addTriangle(positions[0], positions[1], positions[2]);
	btGImpactMeshShape shape(&geometry);
	shape.updateBound();
	body.setCollisionShape(&shape);
	body.mapping.resize(3);
	for (int i = 0; i < 3; ++i)
	{
		body.mapping[i].vertexToTetra = 0;
		body.mapping[i].baryCoordInTetra = btVector4(0, 0, 0, 0);
		body.mapping[i].baryCoordInTetra[i] = 1;
	}
	btDeformableLinearElasticityForce material(1000, 4000, 0, 0);
	material.addSoftBody(&body);
	world.addSoftBody(&body);
	world.addForce(&material);
	for (int step = 0; step < 3; ++step)
	{
		EXPECT_EQ(world.stepSimulation(.002, 0), 1);
		EXPECT_FALSE(world.hasCoupledStepFailed());
	}
	btTransform changed = btTransform::getIdentity();
	changed.setOrigin(btVector3(.2, 0, 0));
	body.setWorldTransform(changed);
	EXPECT_EQ(world.stepSimulation(.002, 0), 1);
	EXPECT_FALSE(world.hasCoupledStepFailed());
	body.mapping[0].baryCoordInTetra = btVector4(0, 1, 0, 0);
	EXPECT_EQ(world.stepSimulation(.002, 0), 0);
	EXPECT_TRUE(world.hasCoupledStepFailed());
	world.removeForce(&material);
	world.removeSoftBody(&body);
}
TEST(DeformableVbd, NativeBarrierQueriesMovingRigidAwayFromSoftSurface)
{
	btTriangleMesh geometry;
	geometry.addTriangle(btVector3(-1000, -1000, 0), btVector3(3000, -1000, 0), btVector3(-1000, 3000, 0));
	btGImpactMeshShape shape(&geometry);
	shape.updateBound();
	btDeformableVbdSolver s;
	setupVbdGroundTest(s, 5);
	// Keep the fixture geometry inside its explicit discovery range.
	s.settings.gap = .005;
	s.barriers.clear();
	s.nativeBarriers.push_back({shape.getMeshPart(0), btTransform::getIdentity(), .5, 1, {}});
	btSphereShape sphere(100);
	btDeformableVbdSolver::Rigid rigid;
	rigid.mass = 1;
	rigid.inertia = btVector3(.01, .01, .01);
	rigid.radius = .1;
	rigid.pose = btTransform::getIdentity();
	rigid.pose.setOrigin(btVector3(0, 0, .101));
	s.rigids.push_back(rigid);
	btTransform collisionPose = rigid.pose;
	collisionPose.setOrigin(rigid.pose.getOrigin() * 1000);
	s.convexBarriers.push_back({&sphere, collisionPose, .5, 2, 0, btTransform::getIdentity()});
	ASSERT_TRUE(s.initialize());
	ASSERT_TRUE(s.step(.002));
	EXPECT_EQ(s.barriers.size(), 1u);
	bool found = false;
	for (const auto &contact : s.contacts)
		found = found || contact.rigidA == 0;
	EXPECT_TRUE(found);
}

TEST(GImpactVertexCache, CurrentOnlyRefreshDoesNotReconstructSafePositions)
{
	btGImpactVertexCache cache;
	int calls = 0;
	btScalar position = 1;
	auto reconstruct = [&](int i, btVector3 &p)
	{
		++calls;
		p = btVector3(position + i, 0, 0);
	};
	cache.beginCurrent(3, reconstruct);
	EXPECT_EQ(calls, 3);
	EXPECT_TRUE(cache.safe(0) == nullptr);
	cache.beginCurrent(3, reconstruct);
	EXPECT_EQ(calls, 3);
	cache.end();
	position = 2;
	EXPECT_EQ(cache.current(0)->x(), 1);
	cache.end();
	cache.beginCurrent(3, reconstruct);
	EXPECT_EQ(calls, 6);
	EXPECT_EQ(cache.current(0)->x(), 2);
	EXPECT_TRUE(cache.safe(0) == nullptr);
	cache.end();
	cache.begin(3,
				[&](int i, btVector3 &p, btVector3 &safe)
				{
					p = btVector3(i, 0, 0);
					safe = btVector3(-i, 0, 0);
				});
	EXPECT_TRUE(cache.safe(0) != nullptr);
	cache.end();
	cache.beginCurrent(3, reconstruct);
	EXPECT_TRUE(cache.safe(0) == nullptr);
	cache.end();
}
TEST(DeformableVbd, OwnershipCacheTracksRelabelingAndTreeRebuilds)
{
	btDeformableVbdCollisionMesh mesh;
	mesh.triangles.resize(4);
	for (int i = 0; i < 4; ++i)
	{
		auto &t = mesh.triangles[i];
		t.m_margin = 0;
		t.m_discoveryPadding = 0;
		t.m_vertices[0] = btVector3(i * 2, 0, 0);
		t.m_vertices[1] = btVector3(i * 2 + 1, 0, 0);
		t.m_vertices[2] = btVector3(i * 2, 1, 0);
	}
	mesh.update();
	std::vector<int> owners = {0, 0, 1, 1};
	auto verify = [&]()
	{
		const auto &labels = mesh.subtreeOwners(owners);
		ASSERT_EQ(labels.size(), size_t(mesh.tree.getNodeCount()));
		for (int i = 0; i < mesh.tree.getNodeCount(); ++i)
			if (mesh.tree.isLeafNode(i))
				EXPECT_EQ(labels[i], owners[mesh.tree.getNodeData(i)]);
			else
			{
				const int a = labels[mesh.tree.getLeftNode(i)], b = labels[mesh.tree.getRightNode(i)];
				EXPECT_EQ(labels[i], a == b ? a : (-2147483647 - 1));
			}
	};
	verify();
	const auto *storage = mesh.subtreeOwners(owners).data();
	mesh.triangles[0].m_vertices[0].setX(-100);
	mesh.update();
	verify();
	EXPECT_EQ(storage, mesh.subtreeOwners(owners).data());
	owners = {1, 0, 1, 0};
	verify();
	mesh.triangles.push_back(mesh.triangles[0]);
	owners.push_back(2);
	mesh.update();
	verify();
}
TEST(DeformableVbd, ParallelCollisionRefitMatchesEverySerialBound)
{
	btDeformableVbdCollisionMesh serial, parallel;
	for (int count : {32771, 17001, 33, 0, 40001})
	{
		serial.triangles.resize(count);
		for (int i = 0; i < count; ++i)
		{
			auto &t = serial.triangles[i];
			t.m_margin = .1;
			t.m_discoveryPadding = .4;
			t.m_vertices[0] = btVector3(i % 137, i / 137, (i % 7) * .1);
			t.m_vertices[1] = t.m_vertices[0] + btVector3(.8, 0, .2);
			t.m_vertices[2] = t.m_vertices[0] + btVector3(0, .7, -.1);
		}
		parallel.triangles = serial.triangles;
		serial.update();
		parallel.update(4);
		for (int workers : {2, 4, int(btMax(2u, std::thread::hardware_concurrency()))})
		{
			for (auto &t : serial.triangles)
				for (auto &p : t.m_vertices)
					p += btVector3(.001 * p.y(), -.003 * p.x(), .01);
			parallel.triangles = serial.triangles;
			serial.update();
			parallel.update(workers);
			if (!count)
				continue;
			ASSERT_EQ(serial.tree.getNodeCount(), parallel.tree.getNodeCount());
			for (int node = 0; node < serial.tree.getNodeCount(); ++node)
			{
				btAABB a, b;
				serial.tree.getNodeBound(node, a);
				parallel.tree.getNodeBound(node, b);
				EXPECT_EQ(a.m_min, b.m_min);
				EXPECT_EQ(a.m_max, b.m_max);
			}
		}
	}
}
TEST(GImpactVertexCache, ParallelCurrentSnapshotMatchesSerialAcrossRefreshes)
{
	btGImpactVertexCache serial, parallel;
	for (int count : {65539, 31, 0, 33001})
		for (int workers : {1, 2, 4, int(btMax(2u, std::thread::hardware_concurrency()))})
		{
			std::vector<btVector3> input(count);
			for (int i = 0; i < count; ++i)
				input[i] = btVector3(i * .1, i * -.3, i * .7);
			auto evaluate = [&](int i, btVector3 &p) { p = input[i] * .3 + btVector3(.1, .2, .3); };
			serial.beginCurrent(count, evaluate);
			parallel.beginCurrent(count, evaluate, [&](int n, const auto &fn) { btVbdParallelGeometry(n, workers, fn); });
			for (int i = 0; i < count; ++i)
			{
				EXPECT_EQ(*serial.current(i), *parallel.current(i));
				EXPECT_TRUE(parallel.safe(i) == nullptr);
			}
			serial.end();
			parallel.end();
		}
}
#include <chrono>
TEST(DeformableVbd, DISABLED_DenseMovingSurfaceBenchmark)
{
	btDeformableVbdCollisionMesh mesh;
	const int vertices = 301797, triangles = 193802;
	std::vector<btDeformableVbdSolver::SurfaceVertex> mapped(vertices);
	std::vector<btVector3> nodes(732);
	for (int i = 0; i < int(nodes.size()); ++i)
		nodes[i] = btVector3(i % 17, i / 17, i % 3);
	for (int i = 0; i < vertices; ++i)
		for (int j = 0; j < 4; ++j)
			mapped[i].support.push_back({(i + j) % 732, btMatrix3x3::getIdentity() * btScalar(.25)});
	mesh.triangles.resize(triangles);
	for (int i = 0; i < triangles; ++i)
	{
		auto &t = mesh.triangles[i];
		t.m_margin = .1;
		t.m_discoveryPadding = .4;
		for (int j = 0; j < 3; ++j)
			t.m_vertices[j] = btVector3(i % 500, i / 500, 0) + btVector3(j == 1, j == 2, 0);
	}
	mesh.update();
	const int logical = int(btMax(2u, std::thread::hardware_concurrency()));
	for (int workers : {1, 2, 4, 8, logical - 1, logical})
	{
		double surface = 0, refit = 0;
		for (int step = 0; step < 55; ++step)
		{
			for (auto &p : nodes)
				p += btVector3(.01, .02, -.03);
			for (auto &t : mesh.triangles)
				for (auto &p : t.m_vertices)
					p += btVector3(.0001 * p.y(), -.0001 * p.x(), .001);
			const auto a = std::chrono::steady_clock::now();
			mesh.candidatePositions.beginCurrent(
				vertices, [&](int i, btVector3 &p) { p = mapped[i].position(nodes); },
				[&](int count, const auto &fn) { btVbdParallelGeometry(count, workers, fn); });
			const auto b = std::chrono::steady_clock::now();
			mesh.update(workers);
			const auto c = std::chrono::steady_clock::now();
			mesh.candidatePositions.end();
			if (step >= 5)
			{
				surface += std::chrono::duration<double, std::milli>(b - a).count();
				refit += std::chrono::duration<double, std::milli>(c - b).count();
			}
		}
		printf("GEOMETRY_BENCH workers=%d surface_ms=%.6f refit_ms=%.6f\n", workers, surface / 50, refit / 50);
	}
}

TEST(DeformableVbd, MappedGuardRefreshesAcrossChangingLoadAndContactReferences)
{
	btDeformableVbdSolver s;
	setupVbdGroundTest(s, .01005);
	s.mappedSurface = true;
	for (int i = 0; i < 3; ++i)
	{
		btDeformableVbdSolver::SurfaceVertex v;
		v.support.push_back({i, btMatrix3x3::getIdentity()});
		v.support.push_back({0, btMatrix3x3::getIdentity() * btScalar(.1)});
		v.support.push_back({3, btMatrix3x3::getIdentity() * btScalar(-.1)});
		s.surfaceVertices.push_back(v);
	}
	s.surface.push_back({{0, 1, 2}});
	s.settings.gap = .0002;
	ASSERT_TRUE(s.initialize());
	for (int step = 0; step < 80; ++step)
	{
		s.velocity.assign(4, btVector3(0, 0, 0));
		s.velocity[step % 4] = btVector3(.2 * ((step % 3) - 1), .1, -1);
		ASSERT_TRUE(s.step(.0002)) << "step=" << step;
		for (const auto &v : s.surfaceVertices)
			EXPECT_GT(v.position(s.x).z(), 0);
		EXPECT_GT(s.minimumJ, 0);
	}
}

TEST(DeformableVbd, NativeCollectionRetainsEarlierTrianglesAndFindsNewRegions)
{
	btTriangleMesh geometry;
	for (int i = 0; i < 2; ++i)
	{
		const btVector3 offset(i * 10000, 0, 0);
		geometry.addTriangle(offset + btVector3(-1000, -1000, 0), offset + btVector3(3000, -1000, 0), offset + btVector3(-1000, 3000, 0));
	}
	btGImpactMeshShape shape(&geometry);
	shape.updateBound();
	btDeformableVbdSolver s;
	setupVbdGroundTest(s, .003);
	// Keep the fixture geometry inside its explicit discovery range.
	s.settings.gap = .005;
	s.external.assign(4, btVector3(0, 0, 0));
	s.barriers.clear();
	s.nativeBarriers.push_back({shape.getMeshPart(0), btTransform::getIdentity(), .5, 1, {}});
	ASSERT_TRUE(s.initialize());
	ASSERT_TRUE(s.step(.002));
	ASSERT_EQ(s.barriers.size(), 1u);
	const auto first = s.barriers[0];
	for (auto &p : s.x)
		p.setX(p.x() + 10);
	s.velocity.assign(4, btVector3(0, 0, 0));
	ASSERT_TRUE(s.initialize());
	ASSERT_TRUE(s.step(.002));
	ASSERT_EQ(s.barriers.size(), 2u);
	for (int i = 0; i < 3; ++i)
		EXPECT_EQ(first.x[i], s.barriers[0].x[i]);
	ASSERT_TRUE(s.step(.002));
	EXPECT_EQ(s.barriers.size(), 2u);
}

TEST(DeformableVbd, DenseMappedMotionGuardMatchesSerialAcrossMappingChanges)
{
	btDeformableVbdSolver serial, parallel;
	for (auto *s : {&serial, &parallel})
	{
		setupVbdGroundTest(*s, .01);
		s->barriers.clear();
		s->external.assign(4, btVector3(0, 0, 0));
		s->mappedSurface = true;
		s->settings.gap = .0002;
		for (int i = 0; i < 32000; ++i)
		{
			btDeformableVbdSolver::SurfaceVertex v;
			v.support.push_back({i % 4, btMatrix3x3::getIdentity() * btScalar(1.1)});
			v.support.push_back({(i + 1) % 4, btMatrix3x3::getIdentity() * btScalar(-.1)});
			s->surfaceVertices.push_back(v);
		}
		s->surface.push_back({{0, 1, 2}});
		ASSERT_TRUE(s->initialize());
	}
	parallel.settings.workers = 8;
	for (int step = 0; step < 8; ++step)
	{
		for (auto *s : {&serial, &parallel})
		{
			if (step == 4)
			{
				for (auto &v : s->surfaceVertices)
					for (auto &support : v.support)
						support.node = (support.node + 1) % 4;
				ASSERT_TRUE(s->initialize());
			}
			s->velocity.assign(4, btVector3(.8, -.1, .2));
			ASSERT_TRUE(s->step(.002));
			EXPECT_GT(s->minimumJ, 0);
		}
		for (int i = 0; i < 4; ++i)
			EXPECT_EQ(serial.x[i], parallel.x[i]);
	}
}
