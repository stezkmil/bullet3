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
TEST(DeformableVbd, WorldMappingCacheStillValidatesSharedVerticesAfterGeometryEdits)
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
	geometry.addTriangle(positions[0], positions[1], positions[2], true);
	geometry.addTriangle(positions[0], positions[2], positions[3], true);
	btGImpactMeshShape shape(&geometry);
	shape.updateBound();
	body.setCollisionShape(&shape);
	body.mapping.resize(4);
	for (int i = 0; i < 4; ++i)
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
	// Geometry positions are deliberately absent from the topology key.
	// Change a vertex used only by the second triangle while reusing that key.
	unsigned char *vertexBase, *indexBase;
	int vertexCount, vertexStride, indexStride, faceCount;
	PHY_ScalarType vertexType, indexType;
	geometry.getLockedVertexIndexBase(&vertexBase, vertexCount, vertexType, vertexStride, &indexBase, indexStride, faceCount, indexType);
	ASSERT_EQ(vertexCount, 4);
	if (vertexType == PHY_DOUBLE)
		reinterpret_cast<double *>(vertexBase + 3 * vertexStride)[0] += .01;
	else
		reinterpret_cast<float *>(vertexBase + 3 * vertexStride)[0] += .01f;
	geometry.unLockVertexBase(0);
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

TEST(DeformableVbd, DISABLED_GpuMappedMotionGuardMatchesCpuAcrossMappingChanges)
{
	btDeformableVbdSolver serial, parallel;
	for (auto *s : {&serial, &parallel})
	{
		setupVbdGroundTest(*s, .01);
		s->barriers.clear();
		s->external.assign(4, btVector3(0, 0, 0));
		s->mappedSurface = true;
		s->settings.gap = .0002;
		for (int i = 0; i < 40000; ++i)
		{
			btDeformableVbdSolver::SurfaceVertex v;
			const btMatrix3x3 matrix(1.1, .02, -.03, -.015, 1.08, .04, .02, .03, 1.12);
			v.support.push_back({i % 4, matrix});
			v.support.push_back({(i + 1) % 4, btMatrix3x3::getIdentity() - matrix});
			v.offset = btVector3(.0001 * (i % 3), -.0002 * (i % 5), .0001 * (i % 7));
			s->surfaceVertices.push_back(v);
		}
		s->surface.push_back({{0, 1, 2}});
		ASSERT_TRUE(s->initialize());
	}
	parallel.settings.workers = 8;
	parallel.settings.gpuGuards = true;
	for (int step = 0; step < 8; ++step)
	{
		for (auto *s : {&serial, &parallel})
		{
			if (step == 4)
			{
				for (auto &v : s->surfaceVertices)
					for (auto &support : v.support)
					{
						support.node = (support.node + 1) % 4;
						support.jacobian = btMatrix3x3::getIdentity() * (support.jacobian[0][0] > 0 ? btScalar(1.1) : btScalar(-.1));
					}
				ASSERT_TRUE(s->initialize());
			}
			else
			{
				const auto topology = s->surfaceTopology();
				ASSERT_TRUE(s->initialize(&topology));
			}
			s->velocity.assign(4, btVector3(.8, -.1, .2));
			ASSERT_TRUE(s->step(.002));
			EXPECT_GT(s->minimumJ, 0);
		}
		for (int i = 0; i < 4; ++i)
			EXPECT_NEAR((serial.x[i] - parallel.x[i]).length(), 0, 1e-10);
		EXPECT_GT(parallel.gpuGuardCalls, 0);
		if (step < 4)
			EXPECT_GT(parallel.gpuSurfaceCalls, 0);
		else
			EXPECT_EQ(parallel.gpuSurfaceCalls, 0);
		EXPECT_EQ(parallel.gpuSurfaceFallbacks, 0);
		EXPECT_EQ(parallel.gpuGuardFallbacks, 0) << parallel.gpuError;
	}
}

TEST(DeformableVbd, DISABLED_GpuMappedGuardRetainsContactAndVolumeProtection)
{
	btDeformableVbdSolver s;
	setupVbdGroundTest(s, .01005);
	s.mappedSurface = true;
	s.settings.gpuGuards = true;
	for (int i = 0; i < 8193; ++i)
	{
		btDeformableVbdSolver::SurfaceVertex v;
		v.support.push_back({i % 3, btMatrix3x3::getIdentity()});
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
	EXPECT_GT(s.gpuGuardCalls, 0);
	EXPECT_EQ(s.gpuGuardFallbacks, 0) << s.gpuError;
}

TEST(DeformableVbd, DISABLED_GpuMappedGuardCoversSharedAndDuplicateSupports)
{
	btDeformableVbdGpu gpu;
	ASSERT_TRUE(gpu.ready()) << gpu.error();
	std::vector<btVbdGpuMapping> maps;
	std::vector<btVbdGpuSupport> supports;
	std::vector<btVbdGpuVec> reference(129), current(129), proposed(129);
	for (int node = 0; node < 129; ++node)
		reference[node] = {.001 * node, -.002 * node, .003 * node};
	for (int v = 0; v < 8192; ++v)
	{
		btVbdGpuMapping map{};
		map.begin = int(supports.size());
		map.offset = {.0001 * (v % 3), -.0002 * (v % 7), .0003 * (v % 5)};
		for (int k = 0; k < 3; ++k)
		{
			btVbdGpuSupport support{};
			support.node = (v + (k == 1 ? 1 : 0)) % 128;
			const double weight = k == 0 ? 1.2 : (k == 1 ? -.3 : .1);
			support.j[0] = support.j[4] = support.j[8] = weight;
			support.j[1] = .02 * weight;
			support.j[5] = -.03 * weight;
			supports.push_back(support);
		}
		map.end = int(supports.size());
		maps.push_back(map);
	}
	ASSERT_TRUE(gpu.mapping(maps, supports)) << gpu.error();
	const double bound = .001;
	auto position = [&](const btVbdGpuMapping &map, const std::vector<btVbdGpuVec> &nodes)
	{
		btVector3 result(map.offset.x, map.offset.y, map.offset.z);
		for (int k = map.begin; k < map.end; ++k)
		{
			const auto &s = supports[k];
			const auto &x = nodes[s.node];
			for (int row = 0; row < 3; ++row)
				result[row] += (s.j[3 * row] * x.x + s.j[3 * row + 1] * x.y) + s.j[3 * row + 2] * x.z;
		}
		return result;
	};
	for (int trial = 0; trial < 9; ++trial)
	{
		if (trial == 2 || trial == 4 || trial == 7)
		{
			for (auto &support : supports)
			{
				support.j[1] = trial == 4 ? .02 * support.j[0] : 0;
				support.j[5] = trial == 4 ? -.03 * support.j[0] : 0;
			}
			ASSERT_TRUE(gpu.mapping(maps, supports)) << gpu.error();
		}
		ASSERT_TRUE(gpu.evaluate(reference, int(maps.size()))) << gpu.error();
		for (size_t v = 0; v < maps.size(); ++v)
		{
			const auto expectedPosition = position(maps[v], reference);
			const auto actualPosition = gpu.positions()[v];
			EXPECT_NEAR(expectedPosition.x(), actualPosition.x, 1e-14);
			EXPECT_NEAR(expectedPosition.y(), actualPosition.y, 1e-14);
			EXPECT_NEAR(expectedPosition.z(), actualPosition.z, 1e-14);
		}

		if (trial == 5)
			for (auto &x : reference)
				x.z += .01;
		if (trial == 6)
		{
			for (auto &s : supports)
				s.node = (s.node + 17) % 128;
			ASSERT_TRUE(gpu.mapping(maps, supports)) << gpu.error();
		}
		current = reference;
		for (auto &x : current)
			x.x += .0001;
		proposed = current;
		for (int node = 0; node < 129; ++node)
			if (trial == 0 || trial == 4 || (trial != 3 && trial != 8 && (node == 13 || (trial == 2 && node == 14))) ||
				(trial == 3 && node == 128))
				proposed[node].y += .004;
		double expected = 1;
		for (const auto &map : maps)
		{
			const btVector3 start = position(map, current), d = position(map, proposed) - start;
			const btVector3 a0 = start - position(map, reference);
			if ((a0 + d).length2() > bound * bound && d.length2() > 0)
			{
				const double a = d.length2(), b = a0.dot(d), c = a0.length2() - bound * bound;
				expected = btMin(expected, btMax(0., (-b + btSqrt(btMax(0., b * b - a * c))) / a));
			}
		}
		double actual = -1;
		ASSERT_TRUE(gpu.guard(current, proposed, reference, trial == 0 || trial == 5, bound, actual)) << gpu.error();
		EXPECT_NEAR(actual, expected, 1e-12) << "trial=" << trial;
	}
}

TEST(DeformableVbd, DISABLED_GpuUnavailableFallsBackToCpu)
{
	btDeformableVbdSolver serial, parallel;
	for (auto *s : {&serial, &parallel})
	{
		setupVbdGroundTest(*s, .01);
		s->barriers.clear();
		s->external.assign(4, btVector3(0, 0, 0));
		s->mappedSurface = true;
		s->settings.gap = .0002;
		for (int i = 0; i < 40000; ++i)
		{
			btDeformableVbdSolver::SurfaceVertex v;
			btMatrix3x3 matrix = btMatrix3x3::getIdentity() * btScalar(1.1);
			matrix[0][1] = .01;
			v.support.push_back({i % 4, matrix});
			v.support.push_back({(i + 1) % 4, btMatrix3x3::getIdentity() * btScalar(-.1)});
			s->surfaceVertices.push_back(v);
		}
		s->surface.push_back({{0, 1, 2}});
		ASSERT_TRUE(s->initialize());
	}
	parallel.settings.workers = 8;
	parallel.settings.gpuGuards = true;
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
		EXPECT_EQ(parallel.gpuGuardCalls, 0);
		EXPECT_GT(parallel.gpuGuardFallbacks, 0);
		EXPECT_GT(parallel.gpuSurfaceFallbacks, 0);
		EXPECT_EQ(parallel.gpuSurfaceCalls, 0);
		EXPECT_FALSE(parallel.gpuError.empty());
	}
}

TEST(DeformableVbd, DISABLED_GpuSurfaceEvaluationMatchesCpuAndBenchmarksDownloads)
{
	btDeformableVbdGpu gpu;
	for (int count : {0, 1, 8193, 301797})
	{
		std::vector<btDeformableVbdSolver::SurfaceVertex> vertices(count);
		std::vector<btVbdGpuMapping> maps(count);
		std::vector<btVbdGpuSupport> supports;
		std::vector<btVector3> x(732), cpu(count), converted(count);
		std::vector<btVbdGpuVec> nodes(732);
		for (int i = 0; i < count; ++i)
		{
			auto &v = vertices[i];
			v.offset = btVector3(.01 * (i % 7), -.02 * (i % 11), .03);
			maps[i].begin = int(supports.size());
			maps[i].offset = {v.offset.x(), v.offset.y(), v.offset.z()};
			for (int j = 0; j < 4; ++j)
			{
				int node = (i + j * 7) % 732;
				btMatrix3x3 matrix(.2, -.01, .03, -.04, .3, .02, .01, -.02, j == 0 ? -.1 : .4);
				v.support.push_back({node, matrix});
				btVbdGpuSupport item{};
				item.node = node;
				for (int r = 0; r < 3; ++r)
					for (int c = 0; c < 3; ++c)
						item.j[3 * r + c] = matrix[r][c];
				supports.push_back(item);
			}
			maps[i].end = int(supports.size());
		}
		ASSERT_TRUE(gpu.mapping(maps, supports)) << gpu.error();
		double cpuMs = 0, gpuMs = 0;
		for (int step = 0; step < 55; ++step)
		{
			for (int i = 0; i < 732; ++i)
				x[i] = btVector3(.01 * (i % 17) + .001 * step, .03 * (i % 11) - .002 * step, .07 * (i % 3) + .003 * step);
			auto a = std::chrono::steady_clock::now();
			btVbdParallelGeometry(count, 8, [&](int i) { cpu[i] = vertices[i].position(x); });
			auto b = std::chrono::steady_clock::now();
			for (int i = 0; i < 732; ++i)
				nodes[i] = {x[i].x(), x[i].y(), x[i].z()};
			ASSERT_TRUE(gpu.evaluate(nodes, count)) << gpu.error();
			btVbdParallelGeometry(count, 8,
								  [&](int i)
								  {
									  const auto &p = gpu.positions()[i];
									  converted[i] = btVector3(p.x, p.y, p.z);
								  });
			auto c = std::chrono::steady_clock::now();
			if (step >= 5)
			{
				cpuMs += std::chrono::duration<double, std::milli>(b - a).count();
				gpuMs += std::chrono::duration<double, std::milli>(c - b).count();
			}
			if (step == 0 || step == 54)
			{
				double difference = 0;
				for (int i = 0; i < count; ++i)
					difference = btMax(difference, double((cpu[i] - converted[i]).length()));
				EXPECT_LT(difference, 1e-14);
			}
		}
		printf("SURFACE_EVAL count=%d cpu8_ms=%.6f gpu_download_convert_ms=%.6f\n", count, cpuMs / 50, gpuMs / 50);
	}
}

TEST(DeformableVbd, ParallelTriangleUpdatesPreservePlanesAndChangeDetection)
{
	btDeformableVbdCollisionMesh serial, parallel;
	for (int count : {32771, 17, 0, 40001})
		for (int change = 0; change < 4; ++change)
		{
			const btScalar margin = change >= 2 ? .2 : .1;
			const btScalar padding = change >= 2 ? .7 : .4;
			auto position = [&](int i, int j)
			{
				return btVector3(i % 137, i / 137, (i % 7) * .1) + btVector3(j == 1 ? .8 : 0, j == 2 ? .7 : 0, j * .2) +
					   (change == 3 ? btVector3(.01, .03, -.05) : btVector3(0, 0, 0));
			};
			const bool a = serial.updateTriangles(count, margin, padding, 1, position);
			const bool b = parallel.updateTriangles(count, margin, padding, 8, position);
			EXPECT_EQ(a, b);
			EXPECT_EQ(a, change == 0 || (count > 0 && change >= 2));
			ASSERT_EQ(serial.triangles.size(), parallel.triangles.size());
			for (int i = 0; i < count; ++i)
			{
				const auto &x = serial.triangles[i], &y = parallel.triangles[i];
				for (int j = 0; j < 3; ++j)
					EXPECT_EQ(x.m_vertices[j], y.m_vertices[j]);
				for (int j = 0; j < 4; ++j)
					EXPECT_EQ(x.m_plane[j], y.m_plane[j]);
				EXPECT_EQ(x.m_margin, y.m_margin);
				EXPECT_EQ(x.m_discoveryPadding, y.m_discoveryPadding);
			}
		}
}

TEST(DeformableVbd, ContiguousPairSortPreservesLegacyOrderAndDuplicates)
{
	std::vector<GIM_PAIR> scratch;
	for (int count : {0, 1, 511, 512, 17000, 50000, 17})
	{
		btPairSet reference;
		btVbdPairSet actual;
		for (int i = count - 1; i >= 0; --i)
			reference.push_back(GIM_PAIR(((i / 2) * 977) % 137 - 61, ((i / 2) * 71) % 509 - 207));
		actual.assign(reference.begin(), reference.end());
		reference.sort([](const GIM_PAIR &a, const GIM_PAIR &b)
					   { return a.m_index1 != b.m_index1 ? a.m_index1 < b.m_index1 : a.m_index2 < b.m_index2; });
		btVbdSortCollisionPairs(actual, scratch);
		ASSERT_EQ(reference.size(), actual.size());
		auto a = actual.begin();
		for (const auto &r : reference)
		{
			EXPECT_EQ(r.m_index1, a->m_index1);
			EXPECT_EQ(r.m_index2, a->m_index2);
			++a;
		}
	}
}
TEST(DeformableVbd, DISABLED_CandidatePairSortBenchmark)
{
	std::vector<GIM_PAIR> scratch;
	for (int count : {17000, 50000})
	{
		btPairSet source;
		for (int i = count - 1; i >= 0; --i)
			source.push_back(GIM_PAIR((i * 977) % 137, (i * 71) % 509));
		double listMs = 0, contiguousMs = 0;
		for (int run = 0; run < 25; ++run)
		{
			auto a = source;
			btVbdPairSet b(source.begin(), source.end());
			auto start = std::chrono::steady_clock::now();
			a.sort([](const GIM_PAIR &x, const GIM_PAIR &y)
				   { return x.m_index1 != y.m_index1 ? x.m_index1 < y.m_index1 : x.m_index2 < y.m_index2; });
			auto middle = std::chrono::steady_clock::now();
			btVbdSortCollisionPairs(b, scratch);
			auto end = std::chrono::steady_clock::now();
			if (run >= 5)
			{
				listMs += std::chrono::duration<double, std::milli>(middle - start).count();
				contiguousMs += std::chrono::duration<double, std::milli>(end - middle).count();
			}
		}
		printf("PAIR_SORT count=%d list_ms=%.6f contiguous_ms=%.6f\n", count, listMs / 20, contiguousMs / 20);
	}
}

#include <limits>
TEST(DeformableVbd, DISABLED_GpuScalarMappingsRejectNonFinitePositions)
{
	for (double invalid : {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
		for (double weight : {0., 1., -.3})
			for (bool evaluate : {false, true})
			{
				btDeformableVbdGpu gpu;
				std::vector<btVbdGpuMapping> maps(1);
				maps[0].end = 1;
				std::vector<btVbdGpuSupport> supports(1);
				supports[0].j[0] = supports[0].j[4] = supports[0].j[8] = weight;
				ASSERT_TRUE(gpu.mapping(maps, supports)) << gpu.error();
				std::vector<btVbdGpuVec> current(1), proposed(1), reference(1);
				proposed[0].y = invalid;
				double limit = 1;
				if (evaluate)
					EXPECT_FALSE(gpu.evaluate(proposed, 1));
				else
					EXPECT_FALSE(gpu.guard(current, proposed, reference, true, .001, limit));
				EXPECT_FALSE(gpu.error().empty());
			}
}

TEST(DeformableVbd, DenseParallelMappingValidationRejectsLateGeometryEdits)
{
	for (int workers : {1, 8})
	{
		btSoftBodyRigidBodyCollisionConfiguration config;
		btCollisionDispatcher dispatcher(&config);
		btDbvtBroadphase broadphase;
		btDeformableBodySolver solver;
		btDeformableMultiBodyConstraintSolver constraints;
		constraints.setDeformableSolver(&solver);
		btDeformableMultiBodyDynamicsWorld world(&dispatcher, &broadphase, &constraints, &config, &solver);
		world.setVbdSolver(true);
		btDeformableVbdSettings settings;
		settings.workers = workers;
		world.setVbdSettings(settings);
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
		for (int i = 0; i < 6001; ++i)
			geometry.addTriangle(positions[0], positions[1], positions[2], false);
		btGImpactMeshShape shape(&geometry);
		shape.updateBound();
		body.setCollisionShape(&shape);
		body.mapping.resize(18003);
		for (int i = 0; i < 18003; ++i)
		{
			body.mapping[i].vertexToTetra = 0;
			body.mapping[i].baryCoordInTetra = btVector4(0, 0, 0, 0);
			body.mapping[i].baryCoordInTetra[i % 3] = 1;
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
		// Geometry positions are deliberately absent from the topology key.
		// Change a vertex in the last parallel block while reusing that key.
		unsigned char *vertexBase, *indexBase;
		int vertexCount, vertexStride, indexStride, faceCount;
		PHY_ScalarType vertexType, indexType;
		geometry.getLockedVertexIndexBase(&vertexBase, vertexCount, vertexType, vertexStride, &indexBase, indexStride, faceCount,
										  indexType);
		ASSERT_EQ(vertexCount, 18003);
		if (vertexType == PHY_DOUBLE)
			reinterpret_cast<double *>(vertexBase + (vertexCount - 1) * vertexStride)[0] += .01;
		else
			reinterpret_cast<float *>(vertexBase + (vertexCount - 1) * vertexStride)[0] += .01f;
		geometry.unLockVertexBase(0);
		EXPECT_EQ(world.stepSimulation(.002, 0), 0);
		EXPECT_TRUE(world.hasCoupledStepFailed());
		world.removeForce(&material);
		world.removeSoftBody(&body);
	}
}

TEST(DeformableVbd, ParallelBvhBuildPreservesEveryNodeAndPrimitiveOrder)
{
	for (int count : {0, 1, 2, 7, 4097, 32771})
		for (int pattern = 0; pattern < 3; ++pattern)
		{
			GIM_BVH_DATA_ARRAY input;
			input.resize(count);
			for (int i = 0; i < count; ++i)
			{
				const btVector3 center =
					pattern == 0 ? btVector3(1, 2, 3)
								 : (pattern == 1 ? btVector3(i * .01, 0, 0) : btVector3((i * 7919) % 101, (i * 3571) % 127, (i * 37) % 61));
				input[i].m_bound.m_min = center - btVector3(.2, .3, .4);
				input[i].m_bound.m_max = center + btVector3(.2, .3, .4);
				input[i].m_data = count - 1 - i;
			}
			auto serialBoxes = input;
			btBvhTree serial;
			if (count)
				serial.build_tree(serialBoxes);
			for (int workers : {1, 8, 31})
			{
				auto boxes = input;
				btBvhTree parallel;
				parallel.build_tree_parallel(boxes, workers * 4, [&](int n, const auto &fn) { btVbdParallelFor(n, workers, fn, 2, 1); });
				ASSERT_EQ(serial.getNodeCount(), parallel.getNodeCount());
				for (int i = 0; i < serial.getNodeCount(); ++i)
				{
					ASSERT_EQ(serial.isLeafNode(i), parallel.isLeafNode(i));
					btAABB a, b;
					serial.getNodeBound(i, a);
					parallel.getNodeBound(i, b);
					EXPECT_EQ(a.m_min, b.m_min);
					EXPECT_EQ(a.m_max, b.m_max);
					if (serial.isLeafNode(i))
						EXPECT_EQ(serial.getNodeData(i), parallel.getNodeData(i));
					else
						EXPECT_EQ(serial.getEscapeNodeIndex(i), parallel.getEscapeNodeIndex(i));
				}
				for (int i = 0; i < count; ++i)
					EXPECT_EQ(serialBoxes[i].m_data, boxes[i].m_data);
			}
		}
}

TEST(DeformableVbd, NativePrimitiveLookupHandlesGrowingSourceAndKeepsEarlierIndices)
{
	btTriangleMesh geometry;
	for (int i = 0; i < 2; ++i)
	{
		const btVector3 offset(i * 10000, 0, 0);
		geometry.addTriangle(offset + btVector3(-1000, -1000, 0), offset + btVector3(3000, -1000, 0), offset + btVector3(-1000, 3000, 0));
	}
	btGImpactMeshShape shape(&geometry);
	shape.updateBound();
	btTriangleMesh firstGeometry;
	firstGeometry.addTriangle(btVector3(-1000, -1000, 0), btVector3(3000, -1000, 0), btVector3(-1000, 3000, 0));
	btGImpactMeshShape firstShape(&firstGeometry);
	firstShape.updateBound();
	btDeformableVbdSolver s;
	setupVbdGroundTest(s, .003);
	// Keep the fixture geometry inside its explicit discovery range.
	s.settings.gap = .005;
	s.external.assign(4, btVector3(0, 0, 0));
	s.barriers.clear();
	s.nativeBarriers.push_back({firstShape.getMeshPart(0), btTransform::getIdentity(), .5, 1, {}});
	ASSERT_TRUE(s.initialize());
	ASSERT_TRUE(s.step(.002));
	ASSERT_EQ(s.barriers.size(), 1u);
	const auto first = s.barriers[0];
	s.nativeBarriers[0].shape = shape.getMeshPart(0);
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

TEST(DeformableVbd, ParallelSurfaceValidationRejectsLateCrossingAmongManySafePairs)
{
	for (int workers : {1, 8, 31})
		for (bool crossing : {false, true})
		{
			btDeformableVbdSolver s;
			setupVbdGroundTest(s, .002);
			s.settings.workers = workers;
			s.settings.gap = .005;
			s.settings.recoveryDistance = 0;
			const auto ground = s.barriers[0];
			s.barriers.assign(8193, ground);
			if (crossing)
			{
				auto &tri = s.barriers.back();
				tri.x[0] = btVector3(.05, -1, -1);
				tri.x[1] = btVector3(.05, 3, -1);
				tri.x[2] = btVector3(.05, -1, 3);
			}
			const auto original = s.x;
			ASSERT_TRUE(s.initialize());
			const bool accepted = s.step(.002);
			EXPECT_EQ(accepted, !crossing) << "workers=" << workers;
			if (crossing)
			{
				EXPECT_STREQ(s.error, "initial_surface_intersection");
				for (int i = 0; i < 4; ++i)
					EXPECT_EQ(original[i], s.x[i]);
			}
			else
				EXPECT_GT(s.minimumJ, .99);
		}
}

TEST(DeformableVbd, CollectedNativeBranchesPreserveEveryTriangleAcrossMotion)
{
	btTriangleMesh geometry;
	for (int group = 0; group < 2; ++group)
		for (int i = 0; i < 64; ++i)
		{
			const btVector3 offset(group * 10000 + i, 0, 0);
			geometry.addTriangle(offset + btVector3(-1000, -1000, 0), offset + btVector3(3000, -1000, 0),
								 offset + btVector3(-1000, 3000, 0));
		}
	btGImpactMeshShape shape(&geometry);
	shape.updateBound();
	btDeformableVbdSolver s;
	setupVbdGroundTest(s, .003);
	s.settings.gap = .005;
	s.barriers.clear();
	s.nativeBarriers.push_back({shape.getMeshPart(0), btTransform::getIdentity(), .5, 1, {}});
	for (int visit = 0; visit < 3; ++visit)
	{
		if (visit)
			for (auto &p : s.x)
				p.setX(p.x() + (visit == 1 ? 10 : -10));
		s.velocity.assign(4, btVector3(0, 0, 0));
		ASSERT_TRUE(s.initialize());
		ASSERT_TRUE(s.step(.002));
		EXPECT_EQ(s.barriers.size(), visit ? 128u : 64u);
		const auto &indices = s.nativeBarriers[0].triangles;
		ASSERT_EQ(indices.size(), 128u);
		std::unordered_set<int> unique;
		for (int i = 0; i < (visit ? 128 : 64); ++i)
		{
			ASSERT_GE(indices[i], 0);
			ASSERT_LT(indices[i], int(s.barriers.size()));
			unique.insert(indices[i]);
			const btScalar expected = btScalar((i / 64) * 10000 + i % 64 - 1000) / s.settings.collisionUnitsPerMeter;
			EXPECT_EQ(s.barriers[indices[i]].x[0].x(), expected);
		}
		EXPECT_EQ(unique.size(), visit ? 128u : 64u);
		ASSERT_TRUE(s.step(.002));
		EXPECT_EQ(s.barriers.size(), visit ? 128u : 64u);
	}
}

TEST(DeformableVbd, ContiguousGImpactTraversalPreservesListCoverageAndOrder)
{
	for (int count : {0, 1, 7, 129, 513})
	{
		btDeformableVbdCollisionMesh a, b;
		for (int side = 0; side < 2; ++side)
		{
			auto &mesh = side ? b : a;
			mesh.updateTriangles(count, .01, .1, 1,
								 [&](int i, int j)
								 {
									 const btVector3 offset((i * 31) % 13, (i * 17) % 11, (i * 7) % 5);
									 return offset + (j == 0 ? btVector3(0, 0, 0) : j == 1 ? btVector3(2, 0, 0) : btVector3(0, 2, 0));
								 });
			mesh.update(1);
		}
		for (int transform = 0; transform < 4; ++transform)
		{
			btTransform ta = btTransform::getIdentity(), tb = btTransform::getIdentity();
			tb.setOrigin(btVector3(transform * .37, -.2 * transform, transform == 3 ? 100 : .05));
			tb.setRotation(btQuaternion(btVector3(0, 0, 1), transform * .19));
			btPairSet reference;
			btVbdPairSet actual;
			// Both APIs append to the caller's existing results.
			reference.push_back({-1, -2});
			actual.push_back({-1, -2});
			btGImpactBvh::find_collision(&a.tree, ta, &b.tree, tb, reference);
			btGImpactBvh::find_collision(&a.tree, ta, &b.tree, tb, actual);
			for (int workers : {1, 8, 31})
			{
				btVbdPairSet parallel{{-1, -2}};
				btGImpactBvh::find_collision_parallel(&a.tree, ta, &b.tree, tb, parallel, workers * 4,
													  [&](int n, const auto &fn) { btVbdParallelFor(n, workers, fn, 2, 1); });
				ASSERT_EQ(actual.size(), parallel.size());
				for (size_t i = 0; i < actual.size(); ++i)
				{
					EXPECT_EQ(actual[i].m_index1, parallel[i].m_index1);
					EXPECT_EQ(actual[i].m_index2, parallel[i].m_index2);
				}
			}
			ASSERT_EQ(reference.size(), actual.size());
			int i = 0;
			for (const auto &pair : reference)
			{
				EXPECT_EQ(pair.m_index1, actual[i].m_index1);
				EXPECT_EQ(pair.m_index2, actual[i].m_index2);
				++i;
			}
		}
	}
}

TEST(DeformableVbd, IndexedPlaneGuardsMatchEquivalentSmallContactSet)
{
	btDeformableVbdSolver sparseGuard, dense;
	for (auto *s : {&sparseGuard, &dense})
	{
		setupVbdGroundTest(*s, .00015);
		s->mappedSurface = true;
		s->settings.gap = .002;
		for (int i = 0; i < 3; ++i)
		{
			btDeformableVbdSolver::SurfaceVertex v;
			v.support.push_back({i, btMatrix3x3::getIdentity() * btScalar(1.1)});
			v.support.push_back({(i + 1) % 3, btMatrix3x3::getIdentity() * btScalar(-.1)});
			v.support.push_back({i, btMatrix3x3::getIdentity() * btScalar(0)});
			s->surfaceVertices.push_back(v);
		}
		s->surface.push_back({{0, 1, 2}});
	}
	const auto ground = dense.barriers[0];
	dense.barriers.assign(700, ground);
	ASSERT_TRUE(sparseGuard.initialize());
	ASSERT_TRUE(dense.initialize());
	for (int step = 0; step < 40; ++step)
	{
		if (step == 20)
		{
			dense.barriers.assign(530, ground);
			ASSERT_TRUE(sparseGuard.initialize());
			ASSERT_TRUE(dense.initialize());
		}
		for (auto *s : {&sparseGuard, &dense})
		{
			s->velocity.assign(4, btVector3(0, 0, 0));
			s->velocity[step % 4] = btVector3(.01 * ((step % 3) - 1), .003, -.02);
			ASSERT_TRUE(s->step(.0002));
			EXPECT_GT(s->minimumJ, 0);
		}
		ASSERT_GE(dense.planes.size(), 512u);
		ASSERT_LT(sparseGuard.planes.size(), 512u);
		ASSERT_EQ(sparseGuard.contacts.size(), dense.contacts.size());
		for (int i = 0; i < 4; ++i)
			EXPECT_EQ(sparseGuard.x[i], dense.x[i]) << "step=" << step << " node=" << i;
	}
}

TEST(DeformableVbd, ParallelDistancesPreserveContactOrderAndMotion)
{
	for (int workers : {8, 31})
	{
		btDeformableVbdSolver serial, parallel;
		for (auto *s : {&serial, &parallel})
		{
			setupVbdGroundTest(*s, .00015);
			s->settings.gap = .002;
			const auto ground = s->barriers[0];
			s->barriers.assign(5000, ground);
			for (int i = 0; i < 5000; ++i)
				for (auto &p : s->barriers[i].x)
					p.setZ(-.0000001 * (i % 11));
			ASSERT_TRUE(s->initialize());
		}
		parallel.settings.workers = workers;
		for (int step = 0; step < 8; ++step)
		{
			for (auto *s : {&serial, &parallel})
			{
				s->velocity.assign(4, btVector3(0, 0, 0));
				s->velocity[step % 4] = btVector3(.01, .003, -.02);
				ASSERT_TRUE(s->step(.0002));
				ASSERT_GE(s->broadphasePairs, 64 * 256);
			}
			ASSERT_EQ(serial.contacts.size(), parallel.contacts.size());
			ASSERT_EQ(serial.planes.size(), parallel.planes.size());
			for (size_t i = 0; i < serial.contacts.size(); ++i)
			{
				const auto &a = serial.contacts[i];
				const auto &b = parallel.contacts[i];
				EXPECT_EQ(a.point, b.point);
				EXPECT_EQ(a.normal, b.normal);
				EXPECT_EQ(a.vertex.offset, b.vertex.offset);
				ASSERT_EQ(a.vertex.support.size(), b.vertex.support.size());
				for (size_t j = 0; j < a.vertex.support.size(); ++j)
				{
					EXPECT_EQ(a.vertex.support[j].node, b.vertex.support[j].node);
					for (int row = 0; row < 3; ++row)
						EXPECT_EQ(a.vertex.support[j].jacobian[row], b.vertex.support[j].jacobian[row]);
				}
			}
			for (size_t i = 0; i < serial.planes.size(); ++i)
			{
				EXPECT_EQ(serial.planes[i].nodes, parallel.planes[i].nodes);
				EXPECT_EQ(serial.planes[i].point, parallel.planes[i].point);
				EXPECT_EQ(serial.planes[i].normal, parallel.planes[i].normal);
			}
			for (int i = 0; i < 4; ++i)
			{
				EXPECT_EQ(serial.x[i], parallel.x[i]);
				EXPECT_EQ(serial.velocity[i], parallel.velocity[i]);
			}
		}
	}
}

TEST(DeformableVbd, ParallelMappedBoundsMatchSerialIncludingSignedWeights)
{
	btDeformableVbdSolver s;
	for (int i = 0; i < 97; ++i)
		s.x.push_back(btVector3((i % 11) - 5, (i % 7) - 3, (i % 13) - 6) * btScalar(.001));
	for (int count : {0, 1, 4095, 4096, 10003})
	{
		s.surfaceVertices.clear();
		for (int i = 0; i < count; ++i)
		{
			btDeformableVbdSolver::SurfaceVertex v;
			v.offset = btVector3(.00001 * (i % 3), -.00003, .00004);
			v.support.push_back({i % 97, btMatrix3x3::getIdentity() * btScalar(1.2)});
			v.support.push_back({(i * 7 + 3) % 97, btMatrix3x3::getIdentity() * btScalar(-.2)});
			s.surfaceVertices.push_back(v);
		}
		for (btScalar scale : {btScalar(.001), btScalar(1), btScalar(2.5)})
		{
			btVector3 expectedLo(-.0001, -.0002, -.0003), expectedHi(.0002, .0001, .0003);
			for (const auto &v : s.surfaceVertices)
			{
				const auto p = v.position(s.x) / scale;
				expectedLo.setMin(p);
				expectedHi.setMax(p);
			}
			for (int workers : {1, 8, 31})
			{
				s.settings.workers = workers;
				btVector3 lo(-.0001, -.0002, -.0003), hi(.0002, .0001, .0003);
				s.includeSurfaceBounds(scale, lo, hi);
				// Optimized division can differ by one ULP between the serial and worker loops.
				for (int d = 0; d < 3; ++d)
				{
					const btScalar tolerance = 8 * std::numeric_limits<btScalar>::epsilon() *
											   btMax(btScalar(1), btMax(btFabs(expectedLo[d]), btFabs(expectedHi[d])));
					EXPECT_NEAR(expectedLo[d], lo[d], tolerance);
					EXPECT_NEAR(expectedHi[d], hi[d], tolerance);
				}
			}
		}
	}
}

TEST(DeformableVbd, SmallSurfaceUsesCpuGuardsWhenGpuAccelerationIsEnabled)
{
	btDeformableVbdSolver cpu, automatic;
	for (auto *s : {&cpu, &automatic})
	{
		setupVbdGroundTest(*s, .01005);
		s->mappedSurface = true;
		s->settings.gap = .0002;
		for (int i = 0; i < 3; ++i)
		{
			btDeformableVbdSolver::SurfaceVertex v;
			v.support.push_back({i, btMatrix3x3::getIdentity()});
			v.support.push_back({0, btMatrix3x3::getIdentity() * btScalar(.1)});
			v.support.push_back({3, btMatrix3x3::getIdentity() * btScalar(-.1)});
			s->surfaceVertices.push_back(v);
		}
		s->surface.push_back({{0, 1, 2}});
		ASSERT_TRUE(s->initialize());
	}
	automatic.settings.gpuGuards = true;
	for (int step = 0; step < 40; ++step)
	{
		for (auto *s : {&cpu, &automatic})
		{
			s->velocity.assign(4, btVector3(0, 0, 0));
			s->velocity[step % 4] = btVector3(.2 * ((step % 3) - 1), .1, -1);
			ASSERT_TRUE(s->step(.0002));
			EXPECT_GT(s->minimumJ, 0);
		}
		for (int i = 0; i < 4; ++i)
			EXPECT_EQ(cpu.x[i], automatic.x[i]);
	}
	EXPECT_EQ(automatic.gpuGuardCalls, 0);
	EXPECT_EQ(automatic.gpuGuardFallbacks, 0);
	EXPECT_TRUE(automatic.gpuError.empty());
}

TEST(DeformableVbd, UnchangedTetrahedraDoNotAlterActiveVolumeGuards)
{
	btDeformableVbdSolver fullScan, dense;
	for (auto *s : {&fullScan, &dense})
	{
		setupVbdGroundTest(*s, 1);
		s->barriers.clear();
		s->mappedSurface = true;
		s->settings.gap = .2;
		s->settings.iterations = 2;
		for (int i = 0; i < 4; ++i)
		{
			btDeformableVbdSolver::SurfaceVertex v;
			v.support.push_back({i, btMatrix3x3::getIdentity()});
			s->surfaceVertices.push_back(v);
		}
		s->surface.push_back({{0, 1, 2}});
	}
	for (int i = 0; i < 4; ++i)
	{
		dense.x.push_back(dense.x[i] + btVector3(2, 0, 0));
		dense.velocity.push_back(btVector3(0, 0, 0));
		dense.external.push_back(btVector3(0, 0, 0));
		dense.mass.push_back(0);
		dense.massDamping.push_back(0);
	}
	auto stationary = dense.tets[0];
	stationary.nodes = {{4, 5, 6, 7}};
	dense.tets.resize(512, stationary);
	ASSERT_TRUE(fullScan.initialize());
	ASSERT_TRUE(dense.initialize());
	for (int step = 0; step < 30; ++step)
	{
		if (step == 10 || step == 20)
		{
			dense.tets.resize(step == 10 ? 511 : 700, stationary);
			ASSERT_TRUE(fullScan.initialize());
			ASSERT_TRUE(dense.initialize());
		}
		for (auto *s : {&fullScan, &dense})
		{
			s->velocity.assign(s->x.size(), btVector3(0, 0, 0));
			s->velocity[step % 4] = btVector3(.2, -.1, -50);
			ASSERT_TRUE(s->step(.002));
			EXPECT_GT(s->minimumJ, 0);
		}
		for (int i = 0; i < 4; ++i)
		{
			EXPECT_EQ(fullScan.x[i], dense.x[i]);
			EXPECT_EQ(fullScan.velocity[i], dense.velocity[i]);
		}
	}
}

TEST(DeformableVbd, DiagonalMappingBoundsCoverSignedAnisotropicDisplacements)
{
	for (const btMatrix3x3 &matrix :
		 {btMatrix3x3::getIdentity(), btMatrix3x3(-2, 0, 0, 0, .5, 0, 0, 0, 3), btMatrix3x3(1, .2, 0, -.3, 2, .4, 0, .1, 1)})
	{
		btDeformableVbdSolver s;
		setupVbdGroundTest(s, 1);
		s.mappedSurface = true;
		for (int i = 0; i < 3; ++i)
		{
			btDeformableVbdSolver::SurfaceVertex v;
			v.support.push_back({i, matrix});
			v.support.push_back({3, matrix * btScalar(-.25)});
			s.surfaceVertices.push_back(v);
		}
		s.surface.push_back({{0, 1, 2}});
		ASSERT_TRUE(s.initialize());
		const auto topology = s.surfaceTopology();
		const bool diagonal = matrix[0][1] == 0;
		const btScalar expected = diagonal ? btMax(btFabs(matrix[0][0]), btMax(btFabs(matrix[1][1]), btFabs(matrix[2][2])))
										   : btSqrt(matrix[0].length2() + matrix[1].length2() + matrix[2].length2());
		EXPECT_NEAR(topology.amplification, 1.25 * expected, 1e-6);
		for (int pass = 0; pass < 200; ++pass)
		{
			std::vector<btVector3> moved = s.x;
			btScalar maximum = 0;
			for (int i = 0; i < 4; ++i)
			{
				btVector3 d(btSin(btScalar(pass + i)), btCos(btScalar(pass * 3 + i)), btSin(btScalar(pass * 7 - i)));
				moved[i] += d;
				maximum = btMax(maximum, d.length());
			}
			for (const auto &v : s.surfaceVertices)
				EXPECT_LE((v.position(moved) - v.position(s.x)).length(), topology.amplification * maximum + btScalar(1e-5));
		}
		ASSERT_TRUE(s.initialize(&topology));
		EXPECT_EQ(s.surfaceTopology().amplification, topology.amplification);
	}
}

TEST(DeformableVbd, ParallelNativeSourcesPreserveCollectionOrderAndMotion)
{
	// Three sources exercise serial extraction; six exceed the aggregate parallel threshold.
	for (int sourceCount : {3, 6})
	{
		btTriangleMesh geometry;
		for (int i = 0; i < 2048; ++i)
		{
			const btVector3 offset((i % 16) * .01, 0, 0);
			geometry.addTriangle(offset + btVector3(-1000, -1000, 0), offset + btVector3(3000, -1000, 0),
								 offset + btVector3(-1000, 3000, 0));
		}
		geometry.addTriangle(btVector3(0, 0, 0), btVector3(0, 0, 0), btVector3(0, 0, 0));
		btGImpactMeshShape shape(&geometry);
		shape.updateBound();
		for (int workers : {8, 31})
		{
			btDeformableVbdSolver serial, parallel;
			for (auto *s : {&serial, &parallel})
			{
				setupVbdGroundTest(*s, .003);
				s->settings.gap = .005;
				s->barriers.clear();
				s->collisionAllowed.assign(sourceCount + 1, std::vector<bool>(sourceCount + 1, true));
				s->collisionAllowed[0][2] = false;
				for (int source = 0; source < sourceCount; ++source)
				{
					auto transform = btTransform::getIdentity();
					transform.setOrigin(btVector3(source * 10, 0, -source * .1));
					s->nativeBarriers.push_back({shape.getMeshPart(0), transform, .25, source + 1, {}});
				}
				ASSERT_TRUE(s->initialize());
			}
			parallel.settings.workers = workers;
			for (int step = 0; step < 6; ++step)
			{
				for (auto *s : {&serial, &parallel})
				{
					s->velocity.assign(4, btVector3(.02, 0, -.1));
					ASSERT_TRUE(s->step(.002));
				}
				ASSERT_EQ(serial.barriers.size(), parallel.barriers.size());
				for (size_t k = 0; k < serial.barriers.size(); ++k)
				{
					EXPECT_EQ(serial.barriers[k].owner, parallel.barriers[k].owner);
					for (int j = 0; j < 3; ++j)
						EXPECT_EQ(serial.barriers[k].x[j], parallel.barriers[k].x[j]);
				}
				for (int source = 0; source < sourceCount; ++source)
					EXPECT_EQ(serial.nativeBarriers[source].triangles, parallel.nativeBarriers[source].triangles);
				EXPECT_TRUE(parallel.nativeBarriers[1].triangles.empty());
				EXPECT_EQ(serial.contacts.size(), parallel.contacts.size());
				EXPECT_EQ(serial.planes.size(), parallel.planes.size());
				for (int i = 0; i < 4; ++i)
				{
					EXPECT_EQ(serial.x[i], parallel.x[i]);
					EXPECT_EQ(serial.velocity[i], parallel.velocity[i]);
				}
			}
		}
	}
}

TEST(DeformableVbd, PersistentGuardStorageInvalidatesAcrossSolverRecreation)
{
	btDeformableVbdCollisionCache cache;
	std::vector<btVector3> state;
	for (int pass = 0; pass < 30; ++pass)
	{
		if (pass == 15)
		{
			cache.guardStamp = ~0u;
			cache.guardReferenceStamp = ~0u;
		}
		btDeformableVbdSolver persistent(&cache), fresh;
		for (auto *s : {&persistent, &fresh})
		{
			setupVbdGroundTest(*s, .00015);
			s->mappedSurface = true;
			s->settings.gap = .0002;
			if (!state.empty())
				s->x = state;
			const int count = pass % 4 == 0 ? 39 : 3;
			for (int i = 0; i < count; ++i)
			{
				btDeformableVbdSolver::SurfaceVertex v;
				v.support.push_back({i % 3, btMatrix3x3::getIdentity() * btScalar(1.1)});
				v.support.push_back({3, btMatrix3x3::getIdentity() * btScalar(-.1)});
				v.offset = btVector3(0, 0, .01 + .000001 * pass);
				s->surfaceVertices.push_back(v);
			}
			s->surface.push_back({{0, 1, 2}});
			ASSERT_TRUE(s->initialize());
			s->velocity.assign(4, btVector3(0, 0, 0));
			s->velocity[pass % 4] = btVector3(.2, -.1, -2);
			ASSERT_TRUE(s->step(.0002));
			EXPECT_GT(s->minimumJ, 0);
		}
		for (int i = 0; i < 4; ++i)
		{
			EXPECT_EQ(persistent.x[i], fresh.x[i]);
			EXPECT_EQ(persistent.velocity[i], fresh.velocity[i]);
		}
		state = fresh.x;
		EXPECT_GT(cache.guardStamp, 0u);
		EXPECT_GT(cache.guardReferenceStamp, 0u);
	}
	EXPECT_FALSE(cache.guardPositions.empty());
}

TEST(DeformableVbd, ContactDeduplicationUsesQuantizedMappingCoefficients)
{
	for (btScalar perturbation : {btScalar(0), btScalar(1e-10), btScalar(1e-5)})
	{
		btDeformableVbdSolver baseline, duplicate;
		for (auto *s : {&baseline, &duplicate})
		{
			setupVbdGroundTest(*s, .00015);
			s->mappedSurface = true;
			s->settings.gap = .002;
			for (int i = 0; i < 3; ++i)
			{
				btDeformableVbdSolver::SurfaceVertex v;
				v.support.push_back({i, btMatrix3x3::getIdentity()});
				s->surfaceVertices.push_back(v);
			}
			s->surface.push_back({{0, 1, 2}});
		}
		for (int i = 0; i < 3; ++i)
		{
			auto v = duplicate.surfaceVertices[i];
			// The perturbed column multiplies a zero coordinate, leaving collision geometry identical.
			v.support[0].jacobian[2][i == 2 ? 0 : 1] = perturbation;
			duplicate.surfaceVertices.push_back(v);
		}
		duplicate.surface.push_back({{3, 4, 5}});
		for (auto *s : {&baseline, &duplicate})
		{
			ASSERT_TRUE(s->initialize());
			ASSERT_TRUE(s->step(.0002));
		}
		ASSERT_FALSE(baseline.contacts.empty());
		if (perturbation < btScalar(1e-8))
			EXPECT_EQ(baseline.contacts.size(), duplicate.contacts.size());
		else
			EXPECT_GT(duplicate.contacts.size(), baseline.contacts.size());
	}
}

TEST(DeformableVbd, ParallelCommonVolumeGuardsMatchSerialForSharedNodes)
{
	for (int workers : {8, 31})
	{
		btDeformableVbdSolver serial, parallel;
		for (auto *s : {&serial, &parallel})
		{
			setupVbdGroundTest(*s, 1);
			s->barriers.clear();
			s->mappedSurface = true;
			s->settings.gap = 1;
			s->settings.iterations = 2;
			const auto tet = s->tets[0];
			s->tets.assign(4200, tet);
			for (int i = 0; i < 4; ++i)
			{
				btDeformableVbdSolver::SurfaceVertex v;
				v.support.push_back({i, btMatrix3x3::getIdentity()});
				s->surfaceVertices.push_back(v);
			}
			s->surface.push_back({{0, 1, 2}});
			ASSERT_TRUE(s->initialize());
		}
		parallel.settings.workers = workers;
		for (int step = 0; step < 6; ++step)
		{
			for (auto *s : {&serial, &parallel})
			{
				s->velocity.assign(4, btVector3(0, 0, 0));
				s->velocity[1] = btVector3(-100, 0, 0);
				s->velocity[2] = btVector3(0, -100, 0);
				ASSERT_TRUE(s->step(.002));
				EXPECT_GT(s->minimumJ, 0);
			}
			EXPECT_EQ(serial.minimumJ, parallel.minimumJ);
			for (int i = 0; i < 4; ++i)
			{
				EXPECT_EQ(serial.x[i], parallel.x[i]);
				EXPECT_EQ(serial.velocity[i], parallel.velocity[i]);
			}
		}
	}
}

TEST(DeformableVbd, SharedSurfaceTopologySurvivesItsOriginalSolver)
{
	auto setup = [](btDeformableVbdSolver &s)
	{
		setupVbdGroundTest(s, .003);
		s.mappedSurface = true;
		s.settings.gap = .005;
		for (int i = 0; i < 3; ++i)
		{
			btDeformableVbdSolver::SurfaceVertex v;
			v.support.push_back({i, btMatrix3x3::getIdentity()});
			s.surfaceVertices.push_back(v);
		}
		s.surface.push_back({{0, 1, 2}});
	};
	btDeformableVbdSolver cached, fresh;
	setup(cached);
	setup(fresh);
	std::weak_ptr<const std::vector<std::set<int>>> lifetime;
	{
		btDeformableVbdSolver original;
		setup(original);
		ASSERT_TRUE(original.initialize());
		const auto topology = original.surfaceTopology();
		lifetime = topology.neighbors;
		ASSERT_TRUE(cached.initialize(&topology));
		EXPECT_EQ(cached.surfaceTopology().neighbors.get(), topology.neighbors.get());
	}
	EXPECT_FALSE(lifetime.expired());
	ASSERT_TRUE(fresh.initialize());
	for (int step = 0; step < 20; ++step)
	{
		for (auto *s : {&cached, &fresh})
		{
			s->velocity.assign(4, btVector3(.03, 0, -.1));
			ASSERT_TRUE(s->step(.002));
		}
		for (int i = 0; i < 4; ++i)
		{
			EXPECT_EQ(cached.x[i], fresh.x[i]);
			EXPECT_EQ(cached.velocity[i], fresh.velocity[i]);
		}
	}
	const auto oldTopology = cached.surfaceTopology();
	ASSERT_TRUE((*oldTopology.neighbors)[3].empty());
	cached.surfaceVertices[0].support.push_back({3, btMatrix3x3::getIdentity() * btScalar(.1)});
	ASSERT_TRUE(cached.initialize());
	const auto newTopology = cached.surfaceTopology();
	EXPECT_NE(oldTopology.neighbors.get(), newTopology.neighbors.get());
	EXPECT_TRUE((*oldTopology.neighbors)[3].empty());
	EXPECT_FALSE((*newTopology.neighbors)[3].empty());
}

TEST(DeformableVbd, PersistentPairWorkspaceDoesNotReuseOldContactResults)
{
	btDeformableVbdCollisionCache cache;
	const int counts[] = {5000, 16, 0, 5200, 400, 5000, 1, 5000};
	for (int pass = 0; pass < 8; ++pass)
	{
		for (auto &entry : cache.pairDistances)
		{
			entry.status = 3;
			entry.distance2 = -1;
			entry.closestSoft = entry.closestRigid = btVector3(1000, -1000, 1000);
		}
		for (auto &entry : cache.pairSortScratch)
			entry.m_index1 = entry.m_index2 = -123;
		btDeformableVbdSolver persistent(&cache), fresh;
		for (auto *s : {&persistent, &fresh})
		{
			setupVbdGroundTest(*s, .00015);
			s->settings.gap = .002;
			s->settings.workers = pass % 3 == 1 ? 1 : (pass % 2 ? 31 : 8);
			const auto ground = s->barriers[0];
			s->barriers.assign(counts[pass], ground);
			for (int i = 0; i < counts[pass]; ++i)
				for (auto &p : s->barriers[i].x)
					p.setZ(-.0000001 * ((i + pass) % 11));
			ASSERT_TRUE(s->initialize());
			s->velocity.assign(4, btVector3(0, 0, 0));
			s->velocity[pass % 4] = btVector3(.01, .003, -.02);
			ASSERT_TRUE(s->step(.0002));
		}
		EXPECT_EQ(persistent.broadphasePairs, fresh.broadphasePairs);
		ASSERT_EQ(persistent.contacts.size(), fresh.contacts.size());
		ASSERT_EQ(persistent.planes.size(), fresh.planes.size());
		for (size_t i = 0; i < persistent.contacts.size(); ++i)
		{
			EXPECT_EQ(persistent.contacts[i].point, fresh.contacts[i].point);
			EXPECT_EQ(persistent.contacts[i].normal, fresh.contacts[i].normal);
		}
		for (size_t i = 0; i < persistent.planes.size(); ++i)
		{
			EXPECT_EQ(persistent.planes[i].nodes, fresh.planes[i].nodes);
			EXPECT_EQ(persistent.planes[i].point, fresh.planes[i].point);
			EXPECT_EQ(persistent.planes[i].normal, fresh.planes[i].normal);
		}
		for (int i = 0; i < 4; ++i)
		{
			EXPECT_EQ(persistent.x[i], fresh.x[i]);
			EXPECT_EQ(persistent.velocity[i], fresh.velocity[i]);
		}
		if (pass == 0)
		{
			EXPECT_GE(cache.pairDistances.size(), size_t(64 * 256));
			EXPECT_FALSE(cache.pairSortScratch.empty());
		}
	}
}

TEST(DeformableVbd, ParallelPairSortPreservesSignedOrderAndDuplicates)
{
	std::vector<GIM_PAIR> scratch;
	for (int count : {0, 511, 512, 131071, 131072, 200000, 600000, 17})
	{
		btVbdPairSet source;
		unsigned random = 12345;
		for (int i = 0; i < count; ++i)
		{
			random = random * 1664525u + 1013904223u;
			const int a = int(random);
			random = random * 1664525u + 1013904223u;
			source.push_back(i % 5 ? GIM_PAIR(a, int(random)) : GIM_PAIR(-1, 2147483647));
		}
		auto expected = source;
		std::sort(expected.begin(), expected.end(), [](const GIM_PAIR &a, const GIM_PAIR &b)
				  { return a.m_index1 != b.m_index1 ? a.m_index1 < b.m_index1 : a.m_index2 < b.m_index2; });
		for (int workers : {1, 8, 16, 31})
		{
			auto actual = source;
			btVbdSortCollisionPairsParallel(actual, scratch, workers);
			ASSERT_EQ(actual.size(), expected.size());
			for (size_t i = 0; i < actual.size(); ++i)
			{
				ASSERT_EQ(actual[i].m_index1, expected[i].m_index1);
				ASSERT_EQ(actual[i].m_index2, expected[i].m_index2);
			}
		}
	}
}

TEST(DeformableVbd, ActiveVolumeGuardsCoverBitmapBoundaries)
{
	for (int activeIndex : {0, 1, 15, 31, 32, 63, 127, 511, 512, 6047})
	{
		btDeformableVbdSolver fullScan, dense;
		for (auto *s : {&fullScan, &dense})
		{
			setupVbdGroundTest(*s, 1);
			s->barriers.clear();
			s->mappedSurface = true;
			s->settings.gap = 1;
			s->settings.iterations = 2;
			for (int i = 0; i < 4; ++i)
			{
				btDeformableVbdSolver::SurfaceVertex v;
				v.support.push_back({i, btMatrix3x3::getIdentity()});
				s->surfaceVertices.push_back(v);
			}
			s->surface.push_back({{0, 1, 2}});
		}
		for (int i = 0; i < 4; ++i)
		{
			dense.x.push_back(dense.x[i] + btVector3(2, 0, 0));
			dense.velocity.push_back(btVector3(0, 0, 0));
			dense.external.push_back(btVector3(0, 0, 0));
			dense.mass.push_back(0);
			dense.massDamping.push_back(0);
		}
		const auto active = dense.tets[0];
		auto stationary = active;
		stationary.nodes = {{4, 5, 6, 7}};
		dense.tets.assign(6048, stationary);
		dense.tets[activeIndex] = active;
		ASSERT_TRUE(fullScan.initialize());
		ASSERT_TRUE(dense.initialize());
		for (int step = 0; step < 4; ++step)
		{
			for (auto *s : {&fullScan, &dense})
			{
				s->velocity.assign(s->x.size(), btVector3(0, 0, 0));
				s->velocity[step] = btVector3(.2, -.1, -500);
				ASSERT_TRUE(s->step(.002));
				EXPECT_GT(s->minimumJ, 0);
			}
			for (int i = 0; i < 4; ++i)
			{
				EXPECT_EQ(fullScan.x[i], dense.x[i]);
				EXPECT_EQ(fullScan.velocity[i], dense.velocity[i]);
			}
		}
	}
}

TEST(DeformableVbd, PackedSurfaceMappingsFollowTopologyChanges)
{
	btDeformableVbdSolver fresh, cached;
	for (auto *s : {&fresh, &cached})
	{
		setupVbdGroundTest(*s, .1);
		s->mappedSurface = true;
		for (int i = 0; i < 4; ++i)
		{
			btDeformableVbdSolver::SurfaceVertex v;
			v.support.push_back({i, btMatrix3x3::getIdentity() * btScalar(1.2)});
			v.support.push_back({(i + 1) % 4, btMatrix3x3::getIdentity() * btScalar(-.2)});
			if (i == 1)
				v.support[0].jacobian[0][1] = .1;
			if (i == 2)
				v.support[0].jacobian[1][1] = 2;
			s->surfaceVertices.push_back(v);
		}
		s->surface.push_back({{0, 1, 2}});
	}
	ASSERT_TRUE(fresh.initialize());
	const auto topology = fresh.surfaceTopology();
	ASSERT_TRUE(bool(topology.scalarMappings));
	ASSERT_EQ(topology.scalarMappings->ranges.size(), 4u);
	EXPECT_EQ(topology.scalarMappings->supports.size(), 4u);
	EXPECT_EQ(topology.scalarMappings->ranges[0][0], 0);
	EXPECT_EQ(topology.scalarMappings->ranges[1][0], -1);
	EXPECT_EQ(topology.scalarMappings->ranges[2][0], -1);
	EXPECT_EQ(topology.scalarMappings->ranges[3][0], 2);
	ASSERT_TRUE(cached.initialize(&topology));
	EXPECT_EQ(cached.surfaceTopology().scalarMappings.get(), topology.scalarMappings.get());
	for (int pass = 0; pass < 8; ++pass)
	{
		for (auto *s : {&fresh, &cached})
		{
			s->velocity.assign(4, btVector3(.01, -.02, -.01));
			ASSERT_TRUE(s->step(.0002));
		}
		for (int i = 0; i < 4; ++i)
			EXPECT_EQ(fresh.x[i], cached.x[i]);
	}
	cached.surfaceVertices[0].support[0].jacobian[1][2] = .3;
	ASSERT_TRUE(cached.initialize());
	EXPECT_EQ(cached.surfaceTopology().scalarMappings->ranges[0][0], -1);
	EXPECT_EQ(topology.scalarMappings->ranges[0][0], 0);
}

TEST(DeformableVbd, PackedScalarPositionsMatchMatrixAccumulation)
{
	unsigned random = 29751;
	auto value = [&]()
	{
		random = random * 1664525u + 1013904223u;
		return btScalar(int(random % 20001) - 10000) / 10000;
	};
	std::vector<btVector3> positions;
	for (int i = 0; i < 128; ++i)
		positions.push_back(btVector3(value(), value(), value()));
	for (int i = 0; i < 4096; ++i)
	{
		btDeformableVbdSolver::SurfaceVertex mapping;
		mapping.offset = btVector3(value(), value(), value());
		std::vector<btVbd::ScalarMappingSupport> packed;
		for (int j = 0; j < 4; ++j)
		{
			const btScalar weight = i % 7 == 0 ? btScalar(0) : value() * 2;
			const int node = (i * 977 + j * 31) % int(positions.size());
			mapping.support.push_back({node, btMatrix3x3::getIdentity() * weight});
			packed.push_back({node, weight});
		}
		const auto expected = mapping.position(positions);
		const auto actual = btVbd::scalarMappingPosition(mapping.offset, positions.data(), packed.data(), packed.data() + packed.size());
		for (int d = 0; d < 4; ++d)
			ASSERT_EQ(expected[d], actual[d]);
	}
}

TEST(DeformableVbd, NativeCollectionFindsMovingTrianglesAmongStaticEntries)
{
	btTriangleMesh geometry;
	for (int region = 0; region < 2; ++region)
	{
		const btVector3 offset(region * 5000, 0, 0);
		geometry.addTriangle(offset + btVector3(-1000, -1000, 0), offset + btVector3(3000, -1000, 0), offset + btVector3(-1000, 3000, 0));
	}
	btGImpactMeshShape shape(&geometry);
	shape.updateBound();
	for (int ignored : {0, 16})
	{
		btDeformableVbdSolver s;
		setupVbdGroundTest(s, 5);
		s.settings.gap = .005;
		s.barriers.clear();
		s.nativeBarriers.push_back({shape.getMeshPart(0), btTransform::getIdentity(), .5, 1, {}});
		btDeformableVbdSolver::Rigid rigid;
		rigid.mass = 1;
		rigid.inertia = btVector3(.01, .01, .01);
		rigid.radius = .1;
		rigid.pose = btTransform::getIdentity();
		s.rigids.push_back(rigid);
		for (int i = 0; i < ignored; ++i)
		{
			btDeformableVbdSolver::Triangle t;
			t.owner = 3;
			t.friction = .5;
			t.x[0] = btVector3(10 + i, 10, 10);
			t.x[1] = t.x[0] + btVector3(.1, 0, 0);
			t.x[2] = t.x[0] + btVector3(0, .1, 0);
			s.barriers.push_back(t);
		}
		btDeformableVbdSolver::Triangle moving;
		moving.owner = 2;
		moving.rigid = 0;
		moving.friction = .5;
		moving.local[0] = btVector3(0, 0, .003);
		moving.local[1] = btVector3(.01, 0, .003);
		moving.local[2] = btVector3(0, .01, .003);
		for (int j = 0; j < 3; ++j)
			moving.x[j] = moving.local[j];
		s.barriers.push_back(moving);
		ASSERT_TRUE(s.initialize());
		for (int region = 0; region < 2; ++region)
		{
			s.rigids[0].pose.setOrigin(btVector3(region * 5, 0, 0));
			ASSERT_TRUE(s.step(.0002));
			EXPECT_EQ(s.barriers.size(), size_t(ignored + 2 + region));
			EXPECT_EQ(s.barriers[ignored].rigid, 0);
			EXPECT_GE(s.nativeBarriers[0].triangles[region], ignored + 1);
		}
	}
}

TEST(DeformableVbd, PackedScalarSurfaceAvoidsGpuRoundTripWithCpuWorkers)
{
	btDeformableVbdSolver s;
	setupVbdGroundTest(s, .01);
	s.barriers.clear();
	s.mappedSurface = true;
	s.settings.workers = 8;
	s.settings.gpuGuards = true;
	s.settings.gap = .002;
	for (int i = 0; i < 40000; ++i)
	{
		btDeformableVbdSolver::SurfaceVertex v;
		v.support.push_back({i % 4, btMatrix3x3::getIdentity()});
		s.surfaceVertices.push_back(v);
	}
	s.surface.push_back({{0, 1, 2}});
	ASSERT_TRUE(s.initialize());
	const auto topology = s.surfaceTopology();
	ASSERT_TRUE(topology.scalarMappings->allScalar);
	ASSERT_TRUE(s.initialize(&topology));
	ASSERT_TRUE(s.step(.0002));
	EXPECT_EQ(s.gpuSurfaceCalls, 0);
	EXPECT_EQ(s.gpuSurfaceFallbacks, 0);
	EXPECT_TRUE(s.gpuError.empty());
	s.surfaceVertices.back().support[0].jacobian[0][1] = .01;
	ASSERT_TRUE(s.initialize());
	EXPECT_FALSE(s.surfaceTopology().scalarMappings->allScalar);
	EXPECT_TRUE(topology.scalarMappings->allScalar);
}

TEST(GImpactVertexCache, ParallelFullSnapshotPreservesCurrentSafeAndNestedQueries)
{
	btGImpactVertexCache serial, parallel;
	for (int pass = 0; pass < 5; ++pass)
	{
		const int count = pass % 2 ? 32777 : 19001;
		std::vector<unsigned char> calls(count, 0);
		auto reconstruct = [&](int i, btVector3 &current, btVector3 &safe)
		{
			current = btVector3(i * .001 + pass, i * .003 - pass, i * .007 + 2 * pass);
			safe = btVector3(current.x() - .1, current.y() + .2, current.z() - .3);
		};
		serial.begin(count, reconstruct);
		parallel.begin(
			count,
			[&](int i, btVector3 &current, btVector3 &safe)
			{
				++calls[i];
				reconstruct(i, current, safe);
			},
			[&](int n, const auto &operation) { btVbdParallelGeometry(n, 8, operation); });
		for (int i = 0; i < count; ++i)
		{
			ASSERT_EQ(calls[i], 1);
			ASSERT_TRUE(parallel.current(i) != nullptr);
			ASSERT_TRUE(parallel.safe(i) != nullptr);
			EXPECT_EQ(*serial.current(i), *parallel.current(i));
			EXPECT_EQ(*serial.safe(i), *parallel.safe(i));
		}
		int nestedDispatches = 0;
		parallel.begin(count, reconstruct, [&](int, const auto &) { ++nestedDispatches; });
		EXPECT_EQ(nestedDispatches, 0);
		parallel.end();
		ASSERT_TRUE(parallel.current(0) != nullptr);
		parallel.end();
		serial.end();
		EXPECT_TRUE(parallel.current(0) == nullptr);
		EXPECT_TRUE(parallel.safe(0) == nullptr);
		parallel.beginCurrent(count,
							  [&](int i, btVector3 &current)
							  {
								  btVector3 safe;
								  reconstruct(i, current, safe);
							  });
		EXPECT_TRUE(parallel.safe(0) == nullptr);
		parallel.begin(count, reconstruct, [&](int n, const auto &operation) { btVbdParallelGeometry(n, 8, operation); });
		ASSERT_TRUE(parallel.safe(count - 1) != nullptr);
		btVector3 current, safe;
		reconstruct(count - 1, current, safe);
		EXPECT_EQ(*parallel.current(count - 1), current);
		EXPECT_EQ(*parallel.safe(count - 1), safe);
		parallel.end();
		parallel.end();
		EXPECT_TRUE(parallel.current(0) == nullptr);
	}
}
