#include "BulletDynamics/Dynamics/btActionInterface.h"
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
#include "BulletCollision/CollisionShapes/btBoxShape.h"
#include "BulletDynamics/ConstraintSolver/btGeneric6DofSpring2Constraint.h"
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
    world.setGpuXpbdTgsImpulseReporting(true);
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
        require(world.getGpuXpbdTgsImpulseError().empty(),"world impulse readback succeeds");
        const auto& reports=world.getGpuXpbdTgsContactImpulses();
        require(disabled?reports.empty():!reports.empty(),"world reports only solved response contacts");
        btScalar impulseZ=0;
        for(const auto& report:reports)
        {
            require((report.bodyA==&body&&report.bodyB==&ground)||(report.bodyB==&body&&report.bodyA==&ground),"reported world body identity");
            impulseZ+=(report.bodyA==&body?1:-1)*report.normal.z()*report.impulse;
        }
        if(!disparateMass)
        {
            btScalar momentumChange=0;
            for(int i=0;i<4;++i) momentumChange+=body.m_nodes[i].m_v.z()+btScalar(.1)/metersPerUnit;
            require(btFabs((impulseZ-momentumChange)*metersPerUnit)<1.e-4,"world impulse converted to scene units matches momentum change");
        }
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
static void testRigidWorld(btScalar metersPerUnit)
{
    btSoftBodyRigidBodyCollisionConfiguration configuration;btCollisionDispatcher dispatcher(&configuration);btDbvtBroadphase broadphase;
    btDeformableMultiBodyConstraintSolver constraints;btDeformableBodySolver solver;
    btDeformableMultiBodyDynamicsWorld world(&dispatcher,&broadphase,&constraints,&configuration,&solver);
    world.setGpuXpbdTgsSolver(true,metersPerUnit,2,20,0);world.setGravity(btVector3(0,0,0));
    const btScalar unit=btScalar(.001)/metersPerUnit;
    btBoxShape shape(btVector3(1,1,1)*unit);btVector3 inertia;shape.calculateLocalInertia(2,inertia);
    btRigidBody body(2,nullptr,&shape,inertia);world.addRigidBody(&body);
    body.applyCentralForce(btVector3(4,0,0)*unit);body.applyTorque(btVector3(0,0,2)*unit*unit);
    world.stepSimulation(.1,0);
    require(!world.hasDeformableStepFailed(),world.getGpuXpbdTgsError().c_str());
    require(btFabs(body.getLinearVelocity().x()/unit-.2)<1.e-7,"world rigid force integration and length units");
    require(btFabs(body.getAngularVelocity().z()-btScalar(.2)*unit*unit/inertia.z())<1.e-7,"world rigid torque and inertia units");
    const auto velocity=body.getLinearVelocity();world.stepSimulation(.1,0);
    require(!world.hasDeformableStepFailed(),"world rigid cached second step");
    require((body.getLinearVelocity()-velocity).length()<1.e-7,"world clears transient rigid force");
    world.setGravity(btVector3(0,-10,0)*unit);world.stepSimulation(.01,0);
    require(!world.hasDeformableStepFailed()&&btFabs(body.getLinearVelocity().y()/unit+.1)<1.e-7,"world applies rigid gravity exactly once");
    world.setGravity(btVector3(0,0,0));
    world.removeRigidBody(&body);
    // Convex rigid geometry must actually collide, rather than merely integrate.
    btBoxShape groundShape(btVector3(10,.5,10)*unit);btRigidBody ground(0,nullptr,&groundShape);
    btTransform groundPose;groundPose.setIdentity();groundPose.setOrigin(btVector3(0,-.5,0)*unit);ground.setWorldTransform(groundPose);
    world.addRigidBody(&ground);body.setWorldTransform(btTransform(btQuaternion::getIdentity(),btVector3(0,1.001,0)*unit));
    body.setLinearVelocity(btVector3(0,-1,0)*unit);body.setAngularVelocity(btVector3(0,0,0));world.addRigidBody(&body);
    for(int i=0;i<10;++i) {world.stepSimulation(.005,0);require(!world.hasDeformableStepFailed(),"world rigid ground step");}
    require(body.getWorldTransform().getOrigin().y()/unit>.99,"world convex rigid ground penetration");
    ground.setCollisionFlags(btCollisionObject::CF_KINEMATIC_OBJECT);ground.setActivationState(DISABLE_DEACTIVATION);
    ground.setInterpolationWorldTransform(ground.getWorldTransform());
    auto movingGround=ground.getWorldTransform();movingGround.getOrigin().setY(movingGround.getOrigin().y()+btScalar(.01)*unit);
    ground.setWorldTransform(movingGround);world.stepSimulation(.01,0);
    require(!world.hasDeformableStepFailed(),"kinematic rigid contact step");
    require(body.getLinearVelocity().y()/unit>.5,"kinematic rigid contact transfers normal velocity");
    world.removeRigidBody(&body);world.removeRigidBody(&ground);
}

static void testAnchoredWorld(btScalar metersPerUnit)
{
    btSoftBodyRigidBodyCollisionConfiguration configuration;btCollisionDispatcher dispatcher(&configuration);btDbvtBroadphase broadphase;
    btDeformableMultiBodyConstraintSolver constraints;btDeformableBodySolver solver;
    btDeformableMultiBodyDynamicsWorld world(&dispatcher,&broadphase,&constraints,&configuration,&solver);
    world.setGpuXpbdTgsSolver(true,metersPerUnit,2,20,0);world.setGravity(btVector3(0,0,0));
    const btScalar unit=btScalar(.001)/metersPerUnit;
    btVector3 positions[]={btVector3(0,0,0),btVector3(unit,0,0),btVector3(0,unit,0),btVector3(0,0,unit)};
    btScalar masses[]={1,1,1,1};MappedBody soft(&world.getWorldInfo(),4,positions,masses);
    soft.appendTetra(0,1,2,3);soft.initializeDmInverse();soft.m_cfg.collisions=0;world.addSoftBody(&soft);
    btDeformableLinearElasticityForce material(1000,1000,0,0);world.addForce(&soft,&material);
    btBoxShape shape(btVector3(unit,unit,unit));btVector3 inertia;shape.calculateLocalInertia(2,inertia);
    btRigidBody rigid(2,nullptr,&shape,inertia);world.addRigidBody(&rigid);
    soft.setIgnoreCollisionCheck(&rigid,true);rigid.setIgnoreCollisionCheck(&soft,true);
    soft.appendDeformableAnchor(0,&rigid,1);
    for(int i=0;i<5;++i)
    {
        rigid.applyCentralForce(btVector3(4*unit,0,0));world.stepSimulation(.001,0);
        require(!world.hasDeformableStepFailed(),"world coupled anchor step");
        btVector3 momentum=rigid.getLinearVelocity()*2;
        for(int n=0;n<4;++n) momentum+=soft.m_nodes[n].m_v;
        require((momentum/unit-btVector3(.004*(i+1),0,0)).length()<1.e-7,"world anchor conserves coupled momentum in both units");
        require((soft.m_nodes[0].m_x-rigid.getWorldTransform().getOrigin()).length()/unit<1.e-7,"world anchor remains attached");
    }
    rigid.setCollisionFlags(btCollisionObject::CF_KINEMATIC_OBJECT);rigid.setActivationState(DISABLE_DEACTIVATION);
    auto target=rigid.getWorldTransform();target.getOrigin().setX(target.getOrigin().x()+btScalar(.01)*unit);rigid.setWorldTransform(target);
    world.stepSimulation(.01,0);
    require(!world.hasDeformableStepFailed(),"kinematic soft anchor step");
    require((soft.m_nodes[0].m_x-target.getOrigin()).length()/unit<1.e-6,"kinematic soft anchor follows prescribed target");
    soft.removeDeformableAnchorByUserIndex(1);
    world.stepSimulation(.001,0);require(!world.hasDeformableStepFailed(),"world anchor removal while GPU state is cached");
    world.removeRigidBody(&rigid);
    world.stepSimulation(.001,0);require(!world.hasDeformableStepFailed(),"world rigid removal while soft body remains");
    auto held=soft.m_nodes[0].m_x;soft.setCollisionFlags(soft.getCollisionFlags()|btCollisionObject::CF_STATIC_OBJECT);
    world.stepSimulation(.001,0);
    require(!world.hasDeformableStepFailed()&&(soft.m_nodes[0].m_x-held).length()/unit<1.e-10,"static soft flag from host callbacks freezes GPU nodes");
    world.removeForce(&soft,&material);world.removeSoftBody(&soft);
}

static void testKinematicJoint(btScalar metersPerUnit)
{
    btSoftBodyRigidBodyCollisionConfiguration configuration;btCollisionDispatcher dispatcher(&configuration);btDbvtBroadphase broadphase;
    btDeformableMultiBodyConstraintSolver constraints;btDeformableBodySolver solver;
    btDeformableMultiBodyDynamicsWorld world(&dispatcher,&broadphase,&constraints,&configuration,&solver);
    world.setGpuXpbdTgsSolver(true,metersPerUnit,2,20,0);world.setGravity(btVector3(0,0,0));
    const btScalar unit=btScalar(.001)/metersPerUnit;
    btBoxShape shape(btVector3(unit,unit,unit));btVector3 inertia;shape.calculateLocalInertia(2,inertia);
    btRigidBody body(2,nullptr,&shape,inertia),kinematic(0,nullptr,&shape);
    kinematic.setCollisionFlags(btCollisionObject::CF_KINEMATIC_OBJECT);kinematic.setActivationState(DISABLE_DEACTIVATION);
    world.addRigidBody(&body);world.addRigidBody(&kinematic);
    btTransform frame;frame.setIdentity();btGeneric6DofSpring2Constraint joint(kinematic,body,frame,frame);
    for(int i=0;i<6;++i){joint.setLimit(i,0,0);joint.setParam(BT_CONSTRAINT_STOP_ERP,1,i);}
    world.addConstraint(&joint,true);
    const btTransform target(btQuaternion(btVector3(0,0,1),.01),btVector3(.01*unit,0,0));kinematic.setWorldTransform(target);
    require(world.stepSimulation(.01,2,.005)==2&&!world.hasDeformableStepFailed(),"kinematic joint advances multiple world substeps");
    if(world.hasDeformableStepFailed()) std::fprintf(stderr,"kinematic: %s\n",world.getGpuXpbdTgsError().c_str());
    require((body.getWorldTransform().getOrigin()-target.getOrigin()).length()/unit<1.e-6,"joint follows interpolated kinematic translation");
    require(body.getWorldTransform().getRotation().angleShortestPath(target.getRotation())<1.e-5,"joint follows kinematic rotation");
    require((joint.getCalculatedTransformB().getOrigin()-body.getWorldTransform().getOrigin()).length()/unit<1.e-8,"joint frame getters are refreshed after GPU solve");
    require(kinematic.getWorldTransform()==target,"kinematic target is not changed by reaction forces");
    world.stepSimulation(.01,0);require(!world.hasDeformableStepFailed()&&body.getLinearVelocity().length()/unit<1.e-6,"stopped kinematic target stops joint motion");
    world.removeConstraint(&joint);world.removeRigidBody(&kinematic);
    for(int order=0;order<6;++order)
    {
        body.setWorldTransform(btTransform(btQuaternion(.1,.2,.3),btVector3(0,0,0)));body.setLinearVelocity(btVector3(0,0,0));body.setAngularVelocity(btVector3(0,0,0));
        btGeneric6DofSpring2Constraint fixed(btTypedConstraint::getFixedBody(),body,frame,frame,static_cast<RotateOrder>(order));
        for(int i=0;i<6;++i){fixed.setLimit(i,0,0);fixed.setParam(BT_CONSTRAINT_STOP_ERP,.5,i);}
        world.addConstraint(&fixed);world.stepSimulation(.01,0);
        require(!world.hasDeformableStepFailed(),"all joint rotation orders supported");
        require(body.getWorldTransform().getRotation().angleShortestPath(btQuaternion::getIdentity())<.001,"joint angular projection converges for every rotation order");
        world.removeConstraint(&fixed);
    }
    world.removeRigidBody(&body);
}

struct CallbackProbe : btActionInterface
{
    static int cpuDispatches;
    static bool nearCallback(btBroadphasePair& pair,btCollisionDispatcher& dispatcher,const btDispatcherInfo& info)
    {
        ++cpuDispatches;
        return btCollisionDispatcher::defaultNearCallback(pair,dispatcher,info);
    }
    btRigidBody* body=nullptr;std::vector<int> order;int solverCalls=0;bool contacts=false;btScalar solverTime=0;
    bool actionContacts=false;
    void updateAction(btCollisionWorld* w,btScalar) override
    {
        order.push_back(5);actionContacts=false;
        for(int i=0;i<w->getDispatcher()->getNumManifolds();++i)
            actionContacts|=w->getDispatcher()->getManifoldByIndexInternal(i)->getNumContacts()>0;
    }
    void debugDraw(btIDebugDraw*) override {}
    static CallbackProbe& get(btDynamicsWorld* w) {return *static_cast<CallbackProbe*>(w->getWorldUserInfo());}
    static void pre(btDynamicsWorld* w,btScalar) {auto& p=get(w);p.order.push_back(1);p.body->applyCentralForce(btVector3(2,0,0));}
    static void collision(btDynamicsWorld* w,btScalar)
    {
        auto& p=get(w);p.order.push_back(2);p.contacts=false;
        for(int i=0;i<w->getDispatcher()->getNumManifolds();++i)
            p.contacts|=w->getDispatcher()->getManifoldByIndexInternal(i)->getNumContacts()>0;
    }
    static void solver(btScalar time,btDeformableMultiBodyDynamicsWorld* w)
    {
        auto& p=get(w);p.order.push_back(++p.solverCalls%2?3:4);p.solverTime=time;
        if(p.solverCalls%2)p.body->applyCentralForce(btVector3(3,0,0));
    }
    static void post(btDynamicsWorld* w,btScalar) {get(w).order.push_back(6);}
};
int CallbackProbe::cpuDispatches=0;
static void testCallbacks()
{
    btSoftBodyRigidBodyCollisionConfiguration configuration;btCollisionDispatcher dispatcher(&configuration);btDbvtBroadphase broadphase;
    btDeformableMultiBodyConstraintSolver constraints;btDeformableBodySolver solver;
    btDeformableMultiBodyDynamicsWorld world(&dispatcher,&broadphase,&constraints,&configuration,&solver);
    CallbackProbe::cpuDispatches=0;dispatcher.setNearCallback(CallbackProbe::nearCallback);
    world.setGpuXpbdTgsSolver(true,1,2,20,0);world.setGravity(btVector3(0,0,0));
    btBoxShape shape(btVector3(1,1,1));btVector3 inertia;shape.calculateLocalInertia(1,inertia);
    btRigidBody body(1,nullptr,&shape,inertia),sensor(0,nullptr,&shape);
    sensor.setCollisionFlags(sensor.getCollisionFlags()|btCollisionObject::CF_NO_CONTACT_RESPONSE);
    world.addRigidBody(&body);world.addRigidBody(&sensor);
    CallbackProbe probe;probe.body=&body;
    world.setInternalTickCallback(CallbackProbe::pre,&probe,true);
    world.setInternalTickCallback(CallbackProbe::collision,&probe,false,true);
    world.setInternalTickCallback(CallbackProbe::post,&probe,false);
    world.setSolverCallback(CallbackProbe::solver);world.addAction(&probe);
    require(world.stepSimulation(.01,0)==1&&!world.hasDeformableStepFailed(),"GPU callback step");
    require(probe.order==std::vector<int>({1,2,3,4,5,6}),"GPU callbacks and actions preserve Bullet order");
    require(probe.contacts,"notification-only contacts are available to callbacks");
    require(btFabs(body.getLinearVelocity().x()-.05)<1.e-8,"callback forces reach the same GPU step without sensor response");
    require(btFabs(probe.solverTime-.01)<1.e-8,"solver callback receives completed physical time");
    const auto queryPose=body.getWorldTransform();const auto queryVelocity=body.getLinearVelocity();const auto queryOrder=probe.order;
    world.performDiscreteCollisionDetection();
    require(dispatcher.getNumManifolds()>0&&body.getWorldTransform()==queryPose&&body.getLinearVelocity()==queryVelocity&&probe.order==queryOrder,
        "explicit GPU collision query retains reports without advancing physics or firing tick callbacks");
    require(CallbackProbe::cpuDispatches==0,"GPU steps and explicit queries do not run CPU narrowphase");
    sensor.setCollisionFlags(sensor.getCollisionFlags()|btCollisionObject::CF_FIND_ONLY_FIRST_TOUCH);
    world.performDiscreteCollisionDetection();
    require(CallbackProbe::cpuDispatches>0&&body.getWorldTransform()==queryPose&&probe.order==queryOrder,
        "explicit touch inspection preserves CPU margin semantics without advancing physics");
    sensor.setCollisionFlags(sensor.getCollisionFlags()&~btCollisionObject::CF_FIND_ONLY_FIRST_TOUCH);
    CallbackProbe::cpuDispatches=0;
    sensor.setWorldTransform(btTransform(btQuaternion::getIdentity(),btVector3(10,0,0)));probe.order.clear();
    require(world.stepSimulation(.01,0)==1&&!world.hasDeformableStepFailed(),"second GPU callback step");
    require(!probe.contacts,"separated notification contacts disappear on next callback");
    require(CallbackProbe::cpuDispatches==0,"GPU stepping after inspection does not retain a CPU collision pass");
    require(probe.order==std::vector<int>({1,2,3,4,5,6}),"callbacks fire once per outer step on cached GPU state");
    const auto beforeCountRequest=body.getWorldTransform();
    for(int mode=0;mode<2;++mode)
    {
        world.setGpuXpbdTgsSolver(true,1,2,20,0,mode);
        sensor.setCollisionFlags(sensor.getCollisionFlags()|btCollisionObject::CF_ONLY_GATHER_CONTACT_COUNTS);probe.order.clear();
        require(world.stepSimulation(.01,0)==0&&world.hasDeformableStepFailed(),"GPU qualities reject contact-count collection");
        require(world.getGpuXpbdTgsError().find("Contact-count collection")!=std::string::npos,"count rejection explains unsupported operation");
        require(body.getWorldTransform()==beforeCountRequest&&probe.order==std::vector<int>({1}),"count rejection precedes collision traversal and preserves poses");
    }
    sensor.setCollisionFlags(sensor.getCollisionFlags()&~btCollisionObject::CF_ONLY_GATHER_CONTACT_COUNTS);
    require(world.stepSimulation(.01,0)==1&&!world.hasDeformableStepFailed(),"clearing count mode allows GPU simulation to resume");
    world.setInternalTickCallback(nullptr,&probe,true);
    world.setInternalTickCallback(nullptr,&probe,false,true);
    world.setInternalTickCallback(nullptr,&probe,false);
    world.setSolverCallback(nullptr);
    world.removeRigidBody(&sensor);sensor.setWorldTransform(body.getWorldTransform());world.addRigidBody(&sensor);
    probe.order.clear();
    require(world.stepSimulation(.01,0)==1&&!world.hasDeformableStepFailed(),"action-only GPU step");
    require(probe.order==std::vector<int>({5})&&probe.actionContacts,"actions receive GPU contacts without tick callbacks");
    world.removeAction(&probe);world.removeRigidBody(&sensor);world.addRigidBody(&sensor);
    world.setInternalTickCallback(CallbackProbe::pre,&probe,true);probe.order.clear();
    require(world.stepSimulation(.01,0)==1&&!world.hasDeformableStepFailed(),"pre-tick-only GPU step");
    require(probe.order==std::vector<int>({1})&&dispatcher.getNumManifolds()>0,"pre-tick-only worlds retain contacts for their next callback");
    world.removeRigidBody(&sensor);world.removeRigidBody(&body);
    require(dispatcher.getNumManifolds()==0,"removing bodies releases GPU callback manifolds");
}

int main()
{
    testCallbacks();
    testRigidWorld(1);testRigidWorld(.001);
    testAnchoredWorld(1);testAnchoredWorld(.001);
    testKinematicJoint(1);testKinematicJoint(.001);

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
