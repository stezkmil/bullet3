#include "BulletSoftBody/btDeformableMultiBodyDynamicsWorld.h"
#include "BulletSoftBody/btDeformableMultiBodyConstraintSolver.h"
#include "BulletSoftBody/btDeformableBodySolver.h"
#include "BulletSoftBody/btDeformableLinearElasticityForce.h"
#include "BulletSoftBody/btDeformableNodalForce.h"
#include "BulletSoftBody/btDeformableMousePickingForce.h"
#include "BulletSoftBody/btSoftBodyRigidBodyCollisionConfiguration.h"
#include "BulletCollision/BroadphaseCollision/btDbvtBroadphase.h"
#include "BulletCollision/CollisionDispatch/btCollisionDispatcher.h"
#include "BulletCollision/CollisionShapes/btTriangleMesh.h"
#include "BulletCollision/CollisionShapes/btBvhTriangleMeshShape.h"
#include "BulletCollision/Gimpact/btGImpactShape.h"
#include <cstdio>
#include <cmath>
static int failures=0;
static void require(bool value,const char* message)
{
    if(!value) { std::fprintf(stderr,"FAIL: %s\n",message);++failures; }
}

struct MappedBody:btSoftBody
{
    std::vector<btVertexToTetraMapping> mapping;
    using btSoftBody::btSoftBody;
    const std::vector<btVertexToTetraMapping>* getCollisionShapeVertexToSimTetra() const override
    { return mapping.empty()?nullptr:&mapping; }
};

static void testWorld(btScalar metersPerUnit,bool mapped=false,bool disparateMass=false,btScalar maxMassRatio=50,bool translated=false)
{
    btSoftBodyRigidBodyCollisionConfiguration configuration;
    btCollisionDispatcher dispatcher(&configuration);
    btDbvtBroadphase broadphase;
    btDeformableMultiBodyConstraintSolver constraints;
    btDeformableBodySolver solver;
    btDeformableMultiBodyDynamicsWorld world(&dispatcher,&broadphase,&constraints,&configuration,&solver);
    world.setGpuXpbdTgsSolver(true,metersPerUnit,2,8,maxMassRatio);
    world.setGravity(btVector3(0,0,0));
    const btScalar size=btScalar(.01)/metersPerUnit;
    btVector3 positions[]={btVector3(0,0,0),btVector3(size,0,0),btVector3(0,size,0),btVector3(0,0,size)};
    const btVector3 translation=translated?btVector3(.4,-.776,.711)/metersPerUnit:btVector3(0,0,0);
    for(auto& p:positions) p+=translation;
    btScalar masses[]={1,1,1,1};
    if(disparateMass) masses[3]=.001;
    btTriangleMesh collisionMesh;
    MappedBody body(&world.getWorldInfo(),4,positions,masses);
    body.appendTetra(0,1,2,3);
    body.initializeDmInverse();body.setFriction(0);body.m_cfg.collisions=0;
    const btScalar surfaceOffset=mapped?btScalar(.00002)/metersPerUnit:0;
    if(mapped)
    {
        const int faces[][3]={{0,2,1},{0,1,3},{0,3,2},{1,2,3}};
        const btVector3 offset(0,0,surfaceOffset);
        for(const auto& face:faces) collisionMesh.addTriangle(positions[face[0]]-offset,positions[face[1]]-offset,positions[face[2]]-offset);
        auto* collisionShape=new btGImpactMeshShape(&collisionMesh);collisionShape->updateBound();
        auto* part=collisionShape->getMeshPart(0);part->lockChildShapes();
        const auto* manager=part->getTrimeshPrimitiveManager();
        for(unsigned i=0;i<manager->get_vertex_count();++i)
        {
            btVector3 point;manager->get_vertex(i,point,false);
            const auto p=(point-translation)/size;
            body.mapping.push_back({0,btVector4(1-p.x()-p.y()-p.z(),p.x(),p.y(),p.z()),btVector4(0,0,0,0)});
        }
        part->unlockChildShapes();
        delete body.getCollisionShape();body.setCollisionShape(collisionShape);
    }
    world.addSoftBody(&body);
    btDeformableLinearElasticityForce material(1000,1000,0,0);
    world.addForce(&body,&material);
    btAlignedObjectArray<int> selected;for(int i=0;i<4;++i) selected.push_back(i);
    btDeformableNodalForce force(&body,selected,btVector3(btScalar(4)/metersPerUnit,0,0));
    world.addForce(&force);
    const btScalar dt=.001;
    require(world.stepSimulation(dt,0)==1,"world completes a forced GPU step");
    if(world.hasDeformableStepFailed()) std::fprintf(stderr,"world error: %s\n",world.getGpuXpbdTgsError().c_str());
    require(!world.hasDeformableStepFailed(),"GPU step accepted");
    btVector3 mean(0,0,0);btScalar totalMass=0;
    for(int i=0;i<4;++i)
    {
        const btScalar effectiveMass=maxMassRatio>0?btMax(masses[i],btScalar(1)/maxMassRatio):masses[i];
        totalMass+=effectiveMass;mean+=body.m_nodes[i].m_v*metersPerUnit*effectiveMass;
    }
    mean/=totalMass;
    require(std::fabs(double(mean.x()-4*dt/totalMass))<2e-6,"nodal force transfers expected conditioned momentum in either unit system");
    require(std::fabs(double(body.m_nodes[3].m_im-1/masses[3]))<1e-6,"conditioning preserves the source nodal mass");
    require(std::fabs(double(mean.y()))<2e-6&&std::fabs(double(mean.z()))<2e-6,"material preserves transverse momentum");
    const btVector3 before=body.m_nodes[0].m_x;
    body.setFriction(-.5);
    require(world.stepSimulation(dt,0)==0,"unsupported input rejects step");
    require(world.hasDeformableStepFailed(),"rejection is visible to caller");
    require(body.m_nodes[0].m_x==before,"rejection preserves accepted position");
    body.setFriction(.5);
    require(world.stepSimulation(dt,0)==1&&!world.hasDeformableStepFailed(),"corrected input can resume");
    world.removeForce(&force);
    btTriangleMesh mesh;
    const btScalar width=10*size, plane=-surfaceOffset-btScalar(.00001)/metersPerUnit;
    mesh.addTriangle(btVector3(-width,-width,plane),btVector3(width,-width,plane),btVector3(width,width,plane));
    mesh.addTriangle(btVector3(-width,-width,plane),btVector3(width,width,plane),btVector3(-width,width,plane));
    btBvhTriangleMeshShape shape(&mesh,true);
    btCollisionObject ground;ground.setCollisionShape(&shape);ground.setCollisionFlags(btCollisionObject::CF_STATIC_OBJECT);
    ground.setWorldTransform(btTransform(btQuaternion::getIdentity(),translation));
    world.addCollisionObject(&ground);
    body.setCheckCollideWith(true);
    for(int pass=0;pass<5;++pass)
    {
        const bool disabled=pass==1||pass==3;
        ground.setCollisionFlags(btCollisionObject::CF_STATIC_OBJECT|(pass==3?btCollisionObject::CF_NO_CONTACT_RESPONSE:0));
        body.setIgnoreCollisionCheck(&ground,pass==1);
        for(int i=0;i<4;++i)
        {
            body.m_nodes[i].m_x=positions[i];body.m_nodes[i].m_q=positions[i];
            body.m_nodes[i].m_v=btVector3(0,0,-btScalar(.1)/metersPerUnit);
        }
        require(world.stepSimulation(dt,0)==1&&!world.hasDeformableStepFailed(),"static mesh contact step accepted");
        if(world.hasDeformableStepFailed()) std::fprintf(stderr,"contact error: %s\n",world.getGpuXpbdTgsError().c_str());
        btScalar minimum=body.m_nodes[0].m_x.z();
        for(int i=1;i<4;++i) minimum=btMin(minimum,body.m_nodes[i].m_x.z());
        const double depth=double((minimum-translation.z()-surfaceOffset-plane)*metersPerUnit);
        require(disabled?depth<-.00005:depth>-.000002,disabled?"disabled pair permits crossing":"enabled/re-enabled pair stops at static mesh");
    }
    world.removeCollisionObject(&ground);
    body.setFriction(0);
    for(int pass=0;pass<3;++pass)
    {
        std::vector<btGpuXpbdTgsWorldTarget> targets;
        for(int i=0;i<4;++i) targets.push_back({&body,i,body.m_nodes[i].m_x+btVector3(size*.01,0,0)});
        require(world.setGpuXpbdTgsKinematicTargets(targets),"world accepts kinematic patch");
        if(pass==1)
        {
            auto invalid=targets;invalid.push_back(targets[0]);
            require(!world.setGpuXpbdTgsKinematicTargets(invalid),"duplicate target rejected without replacing patch");
        }
        require(world.stepSimulation(dt,0)==1&&!world.hasDeformableStepFailed(),"persistent kinematic patch step accepted");
        for(int i=0;i<4;++i)
        {
            require((body.m_nodes[i].m_x-targets[i].position).length()*metersPerUnit<1e-7,"GPU patch reaches requested endpoint");
            require(std::fabs(double(body.m_nodes[i].m_v.x()*metersPerUnit-.1))<2e-5,"GPU patch reports prescribed velocity");
        }
    }
    world.clearGpuXpbdTgsKinematicTargets(&body);
    require(world.stepSimulation(dt,0)==1&&!world.hasDeformableStepFailed(),"kinematic patch releases back to dynamics");
    mean.setZero();
    for(int i=0;i<4;++i)
    {
        const btScalar effectiveMass=maxMassRatio>0?btMax(masses[i],btScalar(1)/maxMassRatio):masses[i];
        mean+=body.m_nodes[i].m_v*metersPerUnit*effectiveMass;
    }
    require(std::fabs(double(mean.x()/totalMass-.1))<2e-5,"released patch retains physical momentum");
    for(int i=0;i<4;++i)
    {
        body.m_nodes[i].m_x=positions[i];body.m_nodes[i].m_q=positions[i];
        body.m_nodes[i].m_v=btVector3(btScalar(.01)/metersPerUnit,0,0);
    }
    require(world.setGpuXpbdTgsSettling(&body,btScalar(1)/metersPerUnit,1000),"world accepts body settling policy");
    require(world.stepSimulation(dt,0)==1&&!world.hasDeformableStepFailed(),"GPU settling step accepted");
    mean.setZero();
    for(int i=0;i<4;++i)
    {
        const btScalar effectiveMass=maxMassRatio>0?btMax(masses[i],btScalar(1)/maxMassRatio):masses[i];
        mean+=body.m_nodes[i].m_v*metersPerUnit*effectiveMass;
    }
    require(std::fabs(double(mean.x()/totalMass-.0025))<2e-6,"GPU settling damps motion with correct units and substep count");
    require(world.setGpuXpbdTgsSettling(nullptr,0,0),"settling can be disabled");
    world.setGpuXpbdTgsSolver(false,metersPerUnit);
    require(body.m_tetraScratches.size()==body.m_tetras.size(),"disabling GPU restores CPU deformation storage");
    const auto& tetra=body.m_tetras[0];
    const btVector3 edge=body.m_nodes[1].m_x-body.m_nodes[0].m_x;
    require((tetra.m_F*btVector3(size,0,0)-edge).length()*metersPerUnit<1e-7,"CPU deformation cache matches accepted GPU pose on disable");
    world.setGpuXpbdTgsSolver(true,metersPerUnit,2,8,maxMassRatio);
    require(world.stepSimulation(dt,0)==1&&!world.hasDeformableStepFailed(),"GPU can restart from accepted pose after switching solvers");
    world.removeForce(&body,&material);world.removeSoftBody(&body);
}

static btVector3 firstForcedDisplacement(const btVector3& translation)
{
    btSoftBodyRigidBodyCollisionConfiguration configuration;
    btCollisionDispatcher dispatcher(&configuration);
    btDbvtBroadphase broadphase;
    btDeformableMultiBodyConstraintSolver constraints;
    btDeformableBodySolver solver;
    btDeformableMultiBodyDynamicsWorld world(&dispatcher,&broadphase,&constraints,&configuration,&solver);
    world.setGpuXpbdTgsSolver(true,.001,16,50);
    world.setGravity(btVector3(0,0,0));
    btVector3 positions[]={translation,translation+btVector3(1,0,0),translation+btVector3(0,1,0),translation+btVector3(0,0,1)};
    btScalar masses[]={.00025,.00025,.00025,.00025};
    btSoftBody body(&world.getWorldInfo(),4,positions,masses);
    body.appendTetra(0,1,2,3);body.initializeDmInverse();body.setFriction(0);body.m_cfg.collisions=0;
    world.addSoftBody(&body);
    btDeformableLinearElasticityForce material(1000,1000,0,0);
    world.addForce(&body,&material);
    btAlignedObjectArray<int> selected;for(int i=0;i<4;++i) selected.push_back(i);
    btDeformableNodalForce force(&body,selected,btVector3(0,10,0));
    world.addForce(&force);
    require(world.stepSimulation(.0002,0)==1&&!world.hasDeformableStepFailed(),"first small-force step accepted");
    btVector3 displacement(0,0,0);
    for(int i=0;i<4;++i) displacement+=(body.m_nodes[i].m_x-positions[i])/4;
    require(displacement.y()>.0001,"small force moves body on its first step at scene coordinates");
    world.removeForce(&force);world.removeForce(&body,&material);world.removeSoftBody(&body);
    return displacement;
}

static void testMousePicking(btScalar metersPerUnit,bool facePicking)
{
    btSoftBodyRigidBodyCollisionConfiguration configuration;
    btCollisionDispatcher dispatcher(&configuration);btDbvtBroadphase broadphase;
    btDeformableMultiBodyConstraintSolver constraints;btDeformableBodySolver solver;
    btDeformableMultiBodyDynamicsWorld world(&dispatcher,&broadphase,&constraints,&configuration,&solver);
    world.setGpuXpbdTgsSolver(true,metersPerUnit,4,8);world.setGravity(btVector3(0,0,0));
    btVector3 origin=btVector3(.4,-.776,.711)/metersPerUnit;
    btScalar size=.01/metersPerUnit;
    btVector3 positions[]={origin,origin+btVector3(size,0,0),origin+btVector3(0,size,0),origin+btVector3(0,0,size)};
    btScalar masses[]={.00025,.00025,.00025,.00025};
    btSoftBody body(&world.getWorldInfo(),4,positions,masses);
    body.appendTetra(0,1,2,3);body.appendFace(0,1,2);body.initializeDmInverse();body.setFriction(0);body.m_cfg.collisions=0;
    world.addSoftBody(&body);
    btDeformableLinearElasticityForce material(1000,1000,0,0);world.addForce(&body,&material);
    btTransform target(btQuaternion::getIdentity(),origin);
    btDeformableMousePickingForce mouse(10,0,facePicking?&body.m_faces[0]:nullptr,
        facePicking?nullptr:&body.m_tetras[0],nullptr,target,.001/metersPerUnit);
    mouse.addSoftBody(&body);world.addForce(&mouse);
    target.setOrigin(origin+btVector3(.1/metersPerUnit,0,0));mouse.setMouseTransform(target);
    for(int i=0;i<10;++i)
        require(world.stepSimulation(.001,0)==1&&!world.hasDeformableStepFailed(),"mouse picking step accepted");
    btVector3 mean(0,0,0);
    for(int i=0;i<4;++i) mean+=body.m_nodes[i].m_v*metersPerUnit/4;
    const btScalar expected=facePicking?.03:.04;
    require(std::fabs(double(mean.x()-expected))<2.e-4,"face/tetra picking respects force cap and scene units");
    world.removeForce(&mouse);
    require(world.stepSimulation(.001,0)==1&&!world.hasDeformableStepFailed(),"mouse release step accepted");
    btVector3 released(0,0,0);
    for(int i=0;i<4;++i) released+=body.m_nodes[i].m_v*metersPerUnit/4;
    require((released-mean).length()<2.e-5,"removing mouse force preserves body momentum");
    world.removeForce(&body,&material);world.removeSoftBody(&body);
}

static void testRefinedBoundary()
{
    btSoftBodyRigidBodyCollisionConfiguration configuration;btCollisionDispatcher dispatcher(&configuration);
    btDbvtBroadphase broadphase;btDeformableMultiBodyConstraintSolver constraints;btDeformableBodySolver solver;
    btDeformableMultiBodyDynamicsWorld world(&dispatcher,&broadphase,&constraints,&configuration,&solver);
    world.setGpuXpbdTgsSolver(true,.001,2,8,50);world.setGravity(btVector3(0,0,0));
    btVector3 positions[]={{-10,-10,.01},{10,-10,.01},{0,0,10.01},{0,10,.01},{0,0,.01}};
    btScalar masses[]={1,1,1,1,1};btTriangleMesh mesh;
    MappedBody body(&world.getWorldInfo(),5,positions,masses);
    int tets[][4]={{0,1,4,2},{1,3,4,2},{3,0,4,2}};
    for(auto& t:tets)
    {
        if((positions[t[1]]-positions[t[0]]).cross(positions[t[2]]-positions[t[0]]).dot(positions[t[3]]-positions[t[0]])<0)
            std::swap(t[0],t[1]);
        body.appendTetra(t[0],t[1],t[2],t[3]);
    }
    body.initializeDmInverse();body.m_cfg.collisions=0;body.setFriction(0);
    mesh.addTriangle(positions[0],positions[3],positions[1]);
    auto* shape=new btGImpactMeshShape(&mesh);shape->updateBound();
    auto* part=shape->getMeshPart(0);part->lockChildShapes();const auto* manager=part->getTrimeshPrimitiveManager();
    for(unsigned i=0;i<manager->get_vertex_count();++i)
    {
        btVector3 p;manager->get_vertex(i,p,false);int id=0;while(id<5&&(positions[id]-p).length2()>1.e-10) ++id;
        bool found=false;
        for(int t=0;t<3&&!found;++t) for(int j=0;j<4;++j) if(tets[t][j]==id)
        { btVector4 weights(0,0,0,0);weights[j]=1;body.mapping.push_back({unsigned(t),weights,btVector4(0,0,0,0)});found=true;break; }
    }
    part->unlockChildShapes();delete body.getCollisionShape();body.setCollisionShape(shape);world.addSoftBody(&body);
    btDeformableLinearElasticityForce material(1,1,0,0);world.addForce(&body,&material);
    btTriangleMesh floor;floor.addTriangle(btVector3(-100,-100,0),btVector3(100,-100,0),btVector3(0,100,0));
    btBvhTriangleMeshShape floorShape(&floor,true);btCollisionObject obstacle;obstacle.setCollisionShape(&floorShape);
    obstacle.setCollisionFlags(btCollisionObject::CF_STATIC_OBJECT);world.addCollisionObject(&obstacle);
    body.m_nodes[4].m_v=btVector3(0,0,-100);
    require(world.stepSimulation(.001,0)==1&&!world.hasDeformableStepFailed(),"refined boundary step accepted");
    require(body.m_nodes[4].m_x.z()>-.002,"a boundary node inside a coarse collision face cannot fall through the floor");
    world.removeCollisionObject(&obstacle);world.removeForce(&body,&material);world.removeSoftBody(&body);
}
int main()
{
    testRefinedBoundary();
    testWorld(1);testWorld(.001);testWorld(1,true);testWorld(.001,true);
    testWorld(1,false,true);testWorld(.001,false,true);testWorld(.001,false,true,0);
    testWorld(1,true,false,50,true);testWorld(.001,true,false,50,true);
    const auto originDisplacement=firstForcedDisplacement(btVector3(0,0,0));
    const auto sceneDisplacement=firstForcedDisplacement(btVector3(400,-776,711));
    require((originDisplacement-sceneDisplacement).length()<1.e-6,"force response is invariant to the scene translation");
    testMousePicking(1,false);testMousePicking(.001,false);
    testMousePicking(1,true);testMousePicking(.001,true);
    std::printf("GPU XPBD world tests: %d failures\n",failures);
    return failures?1:0;
}
