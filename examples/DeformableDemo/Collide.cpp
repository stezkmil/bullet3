/*
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

#include "Collide.h"
///btBulletDynamicsCommon.h is the main Bullet include file, contains most common include files.
#include "btBulletDynamicsCommon.h"
#include "BulletSoftBody/btDeformableMultiBodyDynamicsWorld.h"
#include "BulletSoftBody/btSoftBody.h"
#include "BulletSoftBody/btDeformableGravityForce.h"
#include "BulletSoftBody/btSoftBodyHelpers.h"
#include "BulletSoftBody/btDeformableBodySolver.h"
#include "BulletSoftBody/btSoftBodyRigidBodyCollisionConfiguration.h"
#include "BulletDynamics/Featherstone/btMultiBodyConstraintSolver.h"
#include "../CommonInterfaces/CommonParameterInterface.h"
#include "../OpenGLWindow/SimpleCamera.h"
#include <stdio.h>  //printf debugging
#include <stdlib.h>
#include <sstream>
#include <string>

#include "../CommonInterfaces/CommonDeformableBodyBase.h"
#include "../Utils/b3ResourcePath.h"

///The Collide shows the contact between volumetric deformable objects and rigid objects.
static btScalar E = 68000000;
static btScalar nu = 0.3;
static btScalar damping_alpha = 0.1;
static btScalar damping_beta = 0.01;
static btScalar total_mass = btScalar(5.655);
static btScalar gravity_magnitude = btScalar(10000);
static btScalar gravity_ramp_duration = btScalar(10);
static btScalar internal_time_step = btScalar(0.0002);
static int max_physics_substeps_per_update = 8;
static int use_generated_cantilever_mesh = 0;
static int alternate_tet_split = 0;
static int mesh_cells_x = 76;
static int mesh_cells_y = 2;
static int mesh_cells_z = 2;
static const int kReportInterval = 30;
static const float CAMERA_FRUSTUM_NEAR = 1.f;
static const float CAMERA_FRUSTUM_FAR = 1000000.f;

namespace FooSpace
{
struct OriginalCantileverMesh
{
#include "../SoftDemo/cantilever.inl"
};

static const btScalar kMinX = btScalar(-1010.0009999500051);
static const btScalar kMaxX = btScalar(10.000999950004996);
static const btScalar kMinY = btScalar(-16.580920953957236);
static const btScalar kMaxY = btScalar(10.261236938148294);
static const btScalar kMinZ = btScalar(-6.5809209539572358);
static const btScalar kMaxZ = btScalar(20.261236938148294);

struct GeneratedMesh
{
	std::string m_nodes;
	std::string m_elements;
	int m_numNodes;
	int m_numTets;
};

static int nodeIndex(int ix, int iy, int iz, int nodesY, int nodesZ)
{
	return (ix * nodesY + iy) * nodesZ + iz;
}

struct LocalVertex
{
	int x;
	int y;
	int z;
};

static btScalar signedTetraVolumeTimesSix(const LocalVertex& a,
										  const LocalVertex& b,
										  const LocalVertex& c,
										  const LocalVertex& d)
{
	const btVector3 pa(btScalar(a.x), btScalar(a.y), btScalar(a.z));
	const btVector3 pb(btScalar(b.x), btScalar(b.y), btScalar(b.z));
	const btVector3 pc(btScalar(c.x), btScalar(c.y), btScalar(c.z));
	const btVector3 pd(btScalar(d.x), btScalar(d.y), btScalar(d.z));
	return (pb - pa).cross(pc - pa).dot(pd - pa);
}

static GeneratedMesh createCantileverPrismMesh(int cellsX, int cellsY, int cellsZ, bool alternateSplit)
{
	GeneratedMesh mesh;
	cellsX = btMax(1, cellsX);
	cellsY = btMax(1, cellsY);
	cellsZ = btMax(1, cellsZ);

	const int nodesX = cellsX + 1;
	const int nodesY = cellsY + 1;
	const int nodesZ = cellsZ + 1;
	mesh.m_numNodes = nodesX * nodesY * nodesZ;
	mesh.m_numTets = cellsX * cellsY * cellsZ * 5;

	std::ostringstream nodeStream;
	nodeStream.setf(std::ios::fixed);
	nodeStream.precision(15);
	nodeStream << mesh.m_numNodes << " 3 0 0\n";

	for (int ix = 0; ix < nodesX; ++ix)
	{
		const btScalar fx = btScalar(ix) / cellsX;
		const btScalar x = kMinX + (kMaxX - kMinX) * fx;
		for (int iy = 0; iy < nodesY; ++iy)
		{
			const btScalar fy = btScalar(iy) / cellsY;
			const btScalar y = kMinY + (kMaxY - kMinY) * fy;
			for (int iz = 0; iz < nodesZ; ++iz)
			{
				const btScalar fz = btScalar(iz) / cellsZ;
				const btScalar z = kMinZ + (kMaxZ - kMinZ) * fz;
				const int index = nodeIndex(ix, iy, iz, nodesY, nodesZ);
				nodeStream << index << " " << x << " " << y << " " << z << "\n";
			}
		}
	}

	static const LocalVertex kBaseTets[5][4] = {
		{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
		{{1, 1, 0}, {1, 1, 1}, {0, 1, 0}, {1, 0, 0}},
		{{1, 1, 1}, {0, 0, 1}, {0, 1, 1}, {0, 1, 0}},
		{{1, 0, 1}, {0, 0, 1}, {1, 1, 1}, {1, 0, 0}},
		{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}},
	};

	std::ostringstream elementStream;
	elementStream << mesh.m_numTets << " 4 0\n";
	int tetIndex = 0;
	for (int ix = 0; ix < cellsX; ++ix)
	{
		for (int iy = 0; iy < cellsY; ++iy)
		{
			for (int iz = 0; iz < cellsZ; ++iz)
			{
				const bool reflectX = alternateSplit && ((ix & 1) != 0);
				const bool reflectY = alternateSplit && ((iy & 1) != 0);
				const bool reflectZ = alternateSplit && ((iz & 1) != 0);

				for (int t = 0; t < 5; ++t)
				{
					LocalVertex tetVertices[4];
					int tetNodeIndices[4];
					for (int v = 0; v < 4; ++v)
					{
						LocalVertex vertex = kBaseTets[t][v];
						vertex.x = reflectX ? 1 - vertex.x : vertex.x;
						vertex.y = reflectY ? 1 - vertex.y : vertex.y;
						vertex.z = reflectZ ? 1 - vertex.z : vertex.z;
						tetVertices[v] = vertex;
						tetNodeIndices[v] = nodeIndex(ix + vertex.x, iy + vertex.y, iz + vertex.z, nodesY, nodesZ);
					}

					if (signedTetraVolumeTimesSix(tetVertices[0], tetVertices[1], tetVertices[2], tetVertices[3]) < 0)
					{
						const int tmpNode = tetNodeIndices[1];
						tetNodeIndices[1] = tetNodeIndices[2];
						tetNodeIndices[2] = tmpNode;
					}

					elementStream << tetIndex++ << " "
								  << tetNodeIndices[0] << " "
								  << tetNodeIndices[1] << " "
								  << tetNodeIndices[2] << " "
								  << tetNodeIndices[3] << "\n";
				}
			}
		}
	}

	mesh.m_nodes = nodeStream.str();
	mesh.m_elements = elementStream.str();
	return mesh;
}
}  // namespace FooSpace

class Collide : public CommonDeformableBodyBase
{
	btDeformableLinearElasticityForce* m_linearElasticity;
	btDeformableGravityForce* m_gravityForce;
	btSoftBody* m_softBody;
	FooSpace::GeneratedMesh m_generatedMesh;
	btAlignedObjectArray<int> m_clampedNodes;
	btAlignedObjectArray<int> m_freeEndNodes;
	btScalar m_referenceFreeEndAvgZ;
	btScalar m_elapsedSimulationTime;
	btScalar m_frameTimeAccumulator;
	btScalar m_appliedGravityMagnitude;
	btScalar m_lastIncomingDeltaTime;
	btScalar m_incomingDeltaTimeSum;
	btScalar m_minIncomingDeltaTime;
	btScalar m_maxIncomingDeltaTime;
	int m_incomingDeltaTimeSamples;
	int m_lastAvailableSubsteps;
	int m_lastExecutedSubsteps;
	int m_reportCounter;

public:
	Collide(struct GUIHelperInterface* helper)
		: CommonDeformableBodyBase(helper)
	{
		m_linearElasticity = 0;
		m_gravityForce = 0;
		m_softBody = 0;
		m_referenceFreeEndAvgZ = 0;
		m_elapsedSimulationTime = 0;
		m_frameTimeAccumulator = 0;
		m_appliedGravityMagnitude = 0;
		m_lastIncomingDeltaTime = 0;
		m_incomingDeltaTimeSum = 0;
		m_minIncomingDeltaTime = SIMD_INFINITY;
		m_maxIncomingDeltaTime = 0;
		m_incomingDeltaTimeSamples = 0;
		m_lastAvailableSubsteps = 0;
		m_lastExecutedSubsteps = 0;
		m_reportCounter = 0;
	}
	virtual ~Collide()
	{
	}
	void initPhysics();

	void exitPhysics();
	void processCommandLineArgs(int argc, char* argv[]);
	void initializeCantileverBoundaryConditions();
	void reportCantileverState();

	void resetCamera()
	{
		float dist = 1000;
		float pitch = 0;
		float yaw = 0;
		float targetPos[3] = {0, 0, 0};
		m_guiHelper->resetCamera(dist, yaw, pitch, targetPos[0], targetPos[1], targetPos[2]);

		if (m_guiHelper->getRenderInterface() && m_guiHelper->getRenderInterface()->getActiveCamera())
		{
			SimpleCamera* camera = static_cast<SimpleCamera*>(m_guiHelper->getRenderInterface()->getActiveCamera());
			camera->setCameraFrustumNear(CAMERA_FRUSTUM_NEAR);
			camera->setCameraFrustumFar(CAMERA_FRUSTUM_FAR);
		}
	}

	void stepSimulation(float deltaTime)
	{
		m_linearElasticity->setPoissonRatio(nu);
		m_linearElasticity->setYoungsModulus(E);
		m_linearElasticity->setDamping(damping_alpha, damping_beta);
		m_lastIncomingDeltaTime = deltaTime;
		m_incomingDeltaTimeSum += deltaTime;
		m_minIncomingDeltaTime = btMin(m_minIncomingDeltaTime, btScalar(deltaTime));
		m_maxIncomingDeltaTime = btMax(m_maxIncomingDeltaTime, btScalar(deltaTime));
		++m_incomingDeltaTimeSamples;
		m_frameTimeAccumulator += deltaTime;
		const int availableSubsteps = int((m_frameTimeAccumulator + SIMD_EPSILON) / internal_time_step);
		m_lastAvailableSubsteps = availableSubsteps;
		m_lastExecutedSubsteps = 0;
		if (availableSubsteps <= 0)
		{
			return;
		}

		const int substepsToRun = btMin(max_physics_substeps_per_update, availableSubsteps);
		m_lastExecutedSubsteps = substepsToRun;
		for (int step = 0; step < substepsToRun; ++step)
		{
			const btScalar nextSimulationTime = m_elapsedSimulationTime + internal_time_step;
			btScalar rampFraction = btScalar(1);
			if (gravity_ramp_duration > 0)
			{
				rampFraction = btMin(btScalar(1), nextSimulationTime / gravity_ramp_duration);
			}
			m_appliedGravityMagnitude = gravity_magnitude * rampFraction;
			if (m_gravityForce)
			{
				m_gravityForce->m_gravity = btVector3(0, 0, -m_appliedGravityMagnitude);
			}

			m_dynamicsWorld->stepSimulation(internal_time_step, 0);
			m_elapsedSimulationTime = nextSimulationTime;
			m_frameTimeAccumulator = btMax(btScalar(0), m_frameTimeAccumulator - internal_time_step);
			reportCantileverState();
		}
	}

	virtual void renderScene()
	{
		CommonDeformableBodyBase::renderScene();
		btDeformableMultiBodyDynamicsWorld* deformableWorld = getDeformableDynamicsWorld();

		for (int i = 0; i < deformableWorld->getSoftBodyArray().size(); i++)
		{
			btSoftBody* psb = (btSoftBody*)deformableWorld->getSoftBodyArray()[i];
			{
				btSoftBodyHelpers::DrawFrame(psb, deformableWorld->getDebugDrawer());
				btSoftBodyHelpers::Draw(psb, deformableWorld->getDebugDrawer(), deformableWorld->getDrawFlags());
			}
		}
	}
};

void Collide::initializeCantileverBoundaryConditions()
{
	m_clampedNodes.clear();
	m_freeEndNodes.clear();
	if (!m_softBody)
	{
		return;
	}

	btScalar minX = SIMD_INFINITY;
	btScalar maxX = -SIMD_INFINITY;
	for (int i = 0; i < m_softBody->m_nodes.size(); ++i)
	{
		const btScalar x = m_softBody->m_nodes[i].m_x.x();
		minX = btMin(minX, x);
		maxX = btMax(maxX, x);
	}

	const btScalar faceTolerance = btMax(btScalar(1e-4), (maxX - minX) * btScalar(1e-4));
	m_referenceFreeEndAvgZ = 0;
	for (int i = 0; i < m_softBody->m_nodes.size(); ++i)
	{
		btSoftBody::Node& node = m_softBody->m_nodes[i];
		const btScalar x = node.m_x.x();
		if (btFabs(x - maxX) <= faceTolerance)
		{
			node.m_frozen = 1;
			m_clampedNodes.push_back(i);
		}
		else
		{
			node.m_frozen = 0;
		}

		if (btFabs(x - minX) <= faceTolerance)
		{
			m_freeEndNodes.push_back(i);
			m_referenceFreeEndAvgZ += node.m_x.z();
		}
	}

	if (m_freeEndNodes.size() > 0)
	{
		m_referenceFreeEndAvgZ /= btScalar(m_freeEndNodes.size());
	}
}

void Collide::reportCantileverState()
{
	if (!m_softBody)
	{
		return;
	}

	++m_reportCounter;
	if ((m_reportCounter % kReportInterval) != 0)
	{
		return;
	}

	btScalar freeEndAvgZ = 0;
	btScalar avgSpeed = 0;
	for (int i = 0; i < m_freeEndNodes.size(); ++i)
	{
		const btSoftBody::Node& node = m_softBody->m_nodes[m_freeEndNodes[i]];
		freeEndAvgZ += node.m_x.z();
		avgSpeed += node.m_v.length();
	}

	if (m_freeEndNodes.size() > 0)
	{
		const btScalar invCount = btScalar(1) / m_freeEndNodes.size();
		freeEndAvgZ *= invCount;
		avgSpeed *= invCount;
	}

	const btScalar tipDeflection = m_referenceFreeEndAvgZ - freeEndAvgZ;
	const btScalar avgIncomingDeltaTime = m_incomingDeltaTimeSamples > 0 ? m_incomingDeltaTimeSum / m_incomingDeltaTimeSamples : btScalar(0);
	const btScalar minIncomingDeltaTime = m_incomingDeltaTimeSamples > 0 ? m_minIncomingDeltaTime : btScalar(0);
	fprintf(stdout,
			"cantilever report: generated_mesh=%d alt_split=%d sim_time=%.9f applied_gravity=%.9f tip_deflection=%.9f free_end_avg_z=%.9f avg_speed=%.9f mass=%.9f dt=%.9f frame_dt=%.9f frame_dt_avg=%.9f frame_dt_min=%.9f frame_dt_max=%.9f accum=%.9f avail=%d ran=%d\n",
			use_generated_cantilever_mesh, alternate_tet_split, m_elapsedSimulationTime, m_appliedGravityMagnitude, tipDeflection, freeEndAvgZ, avgSpeed, total_mass, internal_time_step, m_lastIncomingDeltaTime, avgIncomingDeltaTime, minIncomingDeltaTime, m_maxIncomingDeltaTime, m_frameTimeAccumulator, m_lastAvailableSubsteps, m_lastExecutedSubsteps);
	fflush(stdout);
}

void Collide::initPhysics()
{
	m_guiHelper->setUpAxis(2);

	///collision configuration contains default setup for memory, collision setup
	m_collisionConfiguration = new btSoftBodyRigidBodyCollisionConfiguration();

	///use the default collision dispatcher. For parallel processing you can use a diffent dispatcher (see Extras/BulletMultiThreaded)
	m_dispatcher = new btCollisionDispatcher(m_collisionConfiguration);

	m_broadphase = new btDbvtBroadphase();
	btDeformableBodySolver* deformableBodySolver = new btDeformableBodySolver();

	btDeformableMultiBodyConstraintSolver* sol = new btDeformableMultiBodyConstraintSolver();
	sol->setDeformableSolver(deformableBodySolver);
	m_solver = sol;

	m_dynamicsWorld = new btDeformableMultiBodyDynamicsWorld(m_dispatcher, m_broadphase, sol, m_collisionConfiguration, deformableBodySolver);
	const btVector3 zeroGravity(0, 0, 0);
	m_dynamicsWorld->setGravity(zeroGravity);
	getDeformableDynamicsWorld()->getWorldInfo().m_gravity = zeroGravity;
	m_guiHelper->createPhysicsDebugDrawer(m_dynamicsWorld);

	if (use_generated_cantilever_mesh)
	{
		m_generatedMesh = FooSpace::createCantileverPrismMesh(mesh_cells_x, mesh_cells_y, mesh_cells_z, alternate_tet_split != 0);
		m_softBody = btSoftBodyHelpers::CreateFromTetGenData(getDeformableDynamicsWorld()->getWorldInfo(),
															 m_generatedMesh.m_elements.c_str(),
															 0,
															 m_generatedMesh.m_nodes.c_str(),
															 false, true, true);
	}
	else
	{
		m_softBody = btSoftBodyHelpers::CreateFromTetGenData(getDeformableDynamicsWorld()->getWorldInfo(),
															 FooSpace::OriginalCantileverMesh::getElements(),
															 0,
															 FooSpace::OriginalCantileverMesh::getNodes(),
															 false, true, true);
	}
	getDeformableDynamicsWorld()->addSoftBody(m_softBody);
	m_softBody->getCollisionShape()->setMargin(0.1);
	m_softBody->setTotalMass(total_mass);
	m_softBody->m_cfg.kKHR = 1;
	m_softBody->m_cfg.kCHR = 1;
	m_softBody->m_cfg.kDF = 0;
	m_softBody->m_sleepingThreshold = 0;
	btSoftBodyHelpers::generateBoundaryFaces(m_softBody);
	initializeCantileverBoundaryConditions();

	btDeformableLinearElasticityForce* linearElasticity = new btDeformableLinearElasticityForce(100, 100, 0.01);
	m_linearElasticity = linearElasticity;
	getDeformableDynamicsWorld()->addForce(m_softBody, linearElasticity);
	m_forces.push_back(linearElasticity);

	m_appliedGravityMagnitude = (gravity_ramp_duration > 0) ? btScalar(0) : gravity_magnitude;
	const btVector3 explicitGravity(0, 0, -m_appliedGravityMagnitude);
	m_gravityForce = new btDeformableGravityForce(explicitGravity);
	getDeformableDynamicsWorld()->addForce(m_softBody, m_gravityForce);
	m_forces.push_back(m_gravityForce);
	//getDeformableDynamicsWorld()->setImplicit(true);
	//getDeformableDynamicsWorld()->setLineSearch(false);
	//getDeformableDynamicsWorld()->setUseProjection(true);
	//getDeformableDynamicsWorld()->getSolverInfo().m_deformable_erp = 0.3;
	//getDeformableDynamicsWorld()->getSolverInfo().m_deformable_maxErrorReduction = btScalar(200);
	getDeformableDynamicsWorld()->getSolverInfo().m_leastSquaresResidualThreshold = 1e-3;
	//getDeformableDynamicsWorld()->getSolverInfo().m_splitImpulse = true;
	getDeformableDynamicsWorld()->getSolverInfo().m_numIterations = 500;
	m_elapsedSimulationTime = 0;
	m_frameTimeAccumulator = 0;

	fprintf(stdout,
			"cantilever setup: generated_mesh=%d alt_split=%d cells=(%d,%d,%d) nodes=%d tets=%d clamped_nodes=%d free_end_nodes=%d mass=%.9f dt=%.9f max_substeps=%d gravity=%.9f ramp=%.9f damping_alpha=%.9f damping_beta=%.9f E=%.9f nu=%.9f\n",
			use_generated_cantilever_mesh,
			alternate_tet_split,
			mesh_cells_x,
			mesh_cells_y,
			mesh_cells_z,
			m_softBody->m_nodes.size(),
			m_softBody->m_tetras.size(),
			m_clampedNodes.size(),
			m_freeEndNodes.size(),
			total_mass,
			internal_time_step,
			max_physics_substeps_per_update,
			gravity_magnitude,
			gravity_ramp_duration,
			damping_alpha,
			damping_beta,
			E,
			nu);
	fflush(stdout);

	m_guiHelper->autogenerateGraphicsObjects(m_dynamicsWorld);

	//	{
	//		SliderParams slider("Young's Modulus", &E);
	//		slider.m_minVal = 0;
	//		slider.m_maxVal = 2000;
	//		if (m_guiHelper->getParameterInterface())
	//			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(slider);
	//	}
	//	{
	//		SliderParams slider("Poisson Ratio", &nu);
	//		slider.m_minVal = 0.05;
	//		slider.m_maxVal = 0.49;
	//		if (m_guiHelper->getParameterInterface())
	//			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(slider);
	//	}
	//	{
	//		SliderParams slider("Mass Damping", &damping_alpha);
	//		slider.m_minVal = 0.001;
	//		slider.m_maxVal = 0.01;
	//		if (m_guiHelper->getParameterInterface())
	//			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(slider);
	//	}
	//    {
	//        SliderParams slider("Stiffness Damping", &damping_beta);
	//        slider.m_minVal = 0.001;
	//        slider.m_maxVal = 0.01;
	//        if (m_guiHelper->getParameterInterface())
	//            m_guiHelper->getParameterInterface()->registerSliderFloatParameter(slider);
	//    }
	{
		SliderParams slider("Young's Modulus", &E);
		slider.m_minVal = 0;
		slider.m_maxVal = 100'000'000;
		if (m_guiHelper->getParameterInterface())
			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(slider);
	}
	{
		SliderParams slider("Poisson Ratio", &nu);
		slider.m_minVal = 0.05;
		slider.m_maxVal = 0.49;
		if (m_guiHelper->getParameterInterface())
			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(slider);
	}
	{
		SliderParams slider("Mass Damping", &damping_alpha);
		slider.m_minVal = 0;
		slider.m_maxVal = 1;
		if (m_guiHelper->getParameterInterface())
			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(slider);
	}
	{
		SliderParams slider("Stiffness Damping", &damping_beta);
		slider.m_minVal = 0;
		slider.m_maxVal = 0.1;
		if (m_guiHelper->getParameterInterface())
			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(slider);
	}
	{
		SliderParams slider("Gravity Ramp", &gravity_ramp_duration);
		slider.m_minVal = 0;
		slider.m_maxVal = 120;
		if (m_guiHelper->getParameterInterface())
			m_guiHelper->getParameterInterface()->registerSliderFloatParameter(slider);
	}
}

void Collide::processCommandLineArgs(int argc, char* argv[])
{
	for (int i = 0; i < argc; ++i)
	{
		if (!argv[i])
		{
			continue;
		}

		const std::string arg(argv[i]);
		if (arg.find("--use_generated_mesh=") == 0)
		{
			use_generated_cantilever_mesh = atoi(arg.substr(21).c_str()) != 0;
		}
		else if (arg.find("--alternate_tet_split=") == 0)
		{
			alternate_tet_split = atoi(arg.substr(22).c_str()) != 0;
		}
		else if (arg.find("--mesh_cells_x=") == 0)
		{
			mesh_cells_x = btMax(1, atoi(arg.substr(15).c_str()));
		}
		else if (arg.find("--mesh_cells_y=") == 0)
		{
			mesh_cells_y = btMax(1, atoi(arg.substr(15).c_str()));
		}
		else if (arg.find("--mesh_cells_z=") == 0)
		{
			mesh_cells_z = btMax(1, atoi(arg.substr(15).c_str()));
		}
		else if (arg.find("--youngs_modulus=") == 0)
		{
			E = btScalar(atof(arg.substr(17).c_str()));
		}
		else if (arg.find("--poisson_ratio=") == 0)
		{
			nu = btScalar(atof(arg.substr(16).c_str()));
		}
		else if (arg.find("--mass=") == 0)
		{
			total_mass = btScalar(atof(arg.substr(7).c_str()));
		}
		else if (arg.find("--internal_time_step=") == 0)
		{
			internal_time_step = btMax(btScalar(1e-6), btScalar(atof(arg.substr(21).c_str())));
		}
		else if (arg.find("--gravity=") == 0)
		{
			gravity_magnitude = btScalar(atof(arg.substr(10).c_str()));
		}
		else if (arg.find("--gravity_ramp_duration=") == 0)
		{
			gravity_ramp_duration = btMax(btScalar(0), btScalar(atof(arg.substr(24).c_str())));
		}
		else if (arg.find("--max_physics_substeps=") == 0)
		{
			max_physics_substeps_per_update = btMax(1, atoi(arg.substr(23).c_str()));
		}
		else if (arg.find("--damping_alpha=") == 0)
		{
			damping_alpha = btMax(btScalar(0), btScalar(atof(arg.substr(16).c_str())));
		}
		else if (arg.find("--damping_beta=") == 0)
		{
			damping_beta = btMax(btScalar(0), btScalar(atof(arg.substr(15).c_str())));
		}
	}
}

void Collide::exitPhysics()
{
	//cleanup in the reverse order of creation/initialization
	removePickingConstraint();
	//remove the rigidbodies from the dynamics world and delete them
	int i;
	for (i = m_dynamicsWorld->getNumCollisionObjects() - 1; i >= 0; i--)
	{
		btCollisionObject* obj = m_dynamicsWorld->getCollisionObjectArray()[i];
		btRigidBody* body = btRigidBody::upcast(obj);
		if (body && body->getMotionState())
		{
			delete body->getMotionState();
		}
		m_dynamicsWorld->removeCollisionObject(obj);
		delete obj;
	}
	// delete forces
	for (int j = 0; j < m_forces.size(); j++)
	{
		btDeformableLagrangianForce* force = m_forces[j];
		delete force;
	}
	m_forces.clear();

	//delete collision shapes
	for (int j = 0; j < m_collisionShapes.size(); j++)
	{
		btCollisionShape* shape = m_collisionShapes[j];
		delete shape;
	}
	m_collisionShapes.clear();

	delete m_dynamicsWorld;

	delete m_solver;

	delete m_broadphase;

	delete m_dispatcher;

	delete m_collisionConfiguration;
}

class CommonExampleInterface* CollideCreateFunc(struct CommonExampleOptions& options)
{
	return new Collide(options.m_guiHelper);
}
