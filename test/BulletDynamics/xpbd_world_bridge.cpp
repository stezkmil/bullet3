// Test-only ABI adapter: run the existing scenario driver through the actual Bullet world.
#include "BulletSoftBody/btDeformableMultiBodyDynamicsWorld.h"
#include "BulletSoftBody/btDeformableBodySolver.h"
#include "BulletSoftBody/btDeformableLinearElasticityForce.h"
#include "BulletSoftBody/btSoftBodyRigidBodyCollisionConfiguration.h"
#include "BulletSoftBody/btGpuXpbdTgsBackend.h"
#include "BulletCollision/BroadphaseCollision/btDbvtBroadphase.h"
#include "BulletCollision/CollisionDispatch/btCollisionDispatcher.h"
#include "BulletCollision/CollisionShapes/btTriangleIndexVertexArray.h"
#include "BulletCollision/Gimpact/btGImpactShape.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <map>
#include <string>
#include <vector>

namespace
{
btVector3 vector(btGpuXpbdVec p){return btVector3(p.x,p.y,p.z);}
btGpuXpbdVec vector(const btVector3& p){return {float(p.x()),float(p.y()),float(p.z())};}
struct Body:btSoftBody
{
    using btSoftBody::btSoftBody;
    std::vector<btVertexToTetraMapping> maps;
    const std::vector<btVertexToTetraMapping>* getCollisionShapeVertexToSimTetra() const override { return &maps; }
};
struct Range { int begin,end; };
struct Mesh
{
    std::vector<btScalar> points;
    std::vector<int> indices,sourceVertices;
    std::unique_ptr<btTriangleIndexVertexArray> meshInterface;
    std::unique_ptr<btGImpactMeshShape> staticShape;
    std::unique_ptr<btCollisionObject> obstacle;
};
struct State
{
    btSoftBodyRigidBodyCollisionConfiguration configuration;
    btCollisionDispatcher dispatcher{&configuration};
    btDbvtBroadphase broadphase;
    btDeformableMultiBodyConstraintSolver constraints;
    btDeformableBodySolver solver;
    btDeformableMultiBodyDynamicsWorld world{&dispatcher,&broadphase,&constraints,&configuration,&solver};
    std::unique_ptr<btGpuXpbdTgsBackend> cooking;
    std::vector<btGpuXpbdNode> nodes;
    std::vector<btGpuXpbdTet> tets;
    std::vector<btGpuXpbdSurfaceVertex> maps;
    std::vector<btGpuXpbdVec> forces;
    std::vector<btGpuXpbdTarget> targets;
    std::map<int,Range> ranges;
    std::map<int,std::unique_ptr<Body>> bodies;
    std::map<int,std::unique_ptr<btDeformableLinearElasticityForce>> materials;
    std::map<int,std::unique_ptr<Mesh>> meshes;
    std::map<std::array<int,4>,std::pair<int,int>> tetLookup;
    std::string error;
    bool built=false;
    ~State(){clear();}
    int fail(const char* message){error=message;return 0;}
    void clear()
    {
        for(auto& entry:materials) world.removeForce(entry.second.get());
        for(auto& entry:bodies) world.removeSoftBody(entry.second.get());
        for(auto& entry:meshes) if(entry.second->obstacle) world.removeCollisionObject(entry.second->obstacle.get());
        bodies.clear();materials.clear();meshes.clear();tetLookup.clear();world.invalidateGpuXpbdTgs();built=false;
    }
    int owner(int node) const
    {
        for(const auto& entry:ranges) if(node>=entry.second.begin&&node<entry.second.end) return entry.first;
        return -1;
    }
    int updateTargets()
    {
        std::vector<btGpuXpbdTgsWorldTarget> converted;
        for(const auto& target:targets)
        {
            int id=owner(target.node);if(id<0||!bodies.count(id)) return fail("Target outside body ranges");
            converted.push_back({bodies.at(id).get(),target.node-ranges.at(id).begin,vector(target.position)});
        }
        if(!world.setGpuXpbdTgsKinematicTargets(converted)) return fail(world.getGpuXpbdTgsError().c_str());
        return 1;
    }
};
void* create(){return new State;}
void destroy(void* state){delete static_cast<State*>(state);}
const char* error(void* state){return state?static_cast<State*>(state)->error.c_str():"World bridge creation failed";}
int upload(void* value,int count,const btGpuXpbdNode* nodes,int tetCount,const btGpuXpbdTet* tets)
{
    auto& s=*static_cast<State*>(value);
    s.clear();s.ranges.clear();s.targets.clear();s.maps.clear();
    s.cooking=std::make_unique<btGpuXpbdTgsBackend>();
    if(!s.cooking->ready()) return s.fail(s.cooking->error().c_str());
    if(!s.cooking->api()->upload(s.cooking->state(),count,nodes,tetCount,tets)) return s.fail(s.cooking->error().c_str());
    s.nodes.assign(nodes,nodes+count);s.tets.assign(tets,tets+tetCount);s.forces.assign(count,{0,0,0});return 1;
}
int embed(void* value,int count,const btGpuXpbdVec* points,int begin,int end,int body,btGpuXpbdSurfaceVertex* maps)
{
    auto& s=*static_cast<State*>(value);
    if(!s.cooking) return s.fail("Embedding must precede surface construction");
    if(!s.cooking->api()->embed(s.cooking->state(),count,points,begin,end,body,maps)) return s.fail(s.cooking->error().c_str());
    s.ranges[body]={begin,end};return 1;
}
int surface(void* value,int count,const btGpuXpbdSurfaceVertex* maps,int faceCount,const btGpuXpbdSurfaceTriangle* faces,double margin,double clearance)
{
    auto& s=*static_cast<State*>(value);
    if(count<=0||faceCount<=0||!maps||!faces||clearance!=0) return s.fail("Unsupported bridge surface input");
    s.clear();s.maps.assign(maps,maps+count);
    for(const auto& entry:s.ranges)
    {
        const int id=entry.first,begin=entry.second.begin,end=entry.second.end;
        std::vector<btVector3> positions;std::vector<btScalar> masses;
        for(int n=begin;n<end;++n){positions.push_back(vector(s.nodes[n].position));masses.push_back(s.nodes[n].inverseMass>0?1/btScalar(s.nodes[n].inverseMass):0);}
        auto body=std::make_unique<Body>(&s.world.getWorldInfo(),end-begin,positions.data(),masses.data());
        body->m_cfg.collisions=0;body->setFriction(0);
        const btGpuXpbdTet* first=nullptr;
        std::vector<const btGpuXpbdTet*> selected;
        for(const auto& tet:s.tets) if(s.owner(tet.nodes[0])==id)
        {
            for(int j=1;j<4;++j) if(s.owner(tet.nodes[j])!=id) return s.fail("Tetrahedron crosses body ranges");
            if(!first) first=&tet;
            if(tet.young!=first->young||tet.poisson!=first->poisson||tet.damping!=first->damping) return s.fail("Bridge requires uniform body material");
            std::array<int,4> key;std::copy(tet.nodes,tet.nodes+4,key.begin());std::sort(key.begin(),key.end());
            s.tetLookup[key]={id,body->m_tetras.size()};
            body->appendTetra(tet.nodes[0]-begin,tet.nodes[1]-begin,tet.nodes[2]-begin,tet.nodes[3]-begin);selected.push_back(&tet);
        }
        if(!first) return s.fail("Body has no tetrahedra");
        body->initializeDmInverse();
        for(int n=begin;n<end;++n){body->m_nodes[n-begin].m_im=s.nodes[n].inverseMass;body->m_nodes[n-begin].m_v=vector(s.nodes[n].velocity);}
        for(int t=0;t<int(selected.size());++t)
        {
            auto& output=body->m_tetras[t];const auto& input=*selected[t];output.m_element_measure=input.volume;
            for(int c=0;c<3;++c)for(int r=0;r<3;++r)output.m_Dm_inverse[r][c]=input.inverseRest[3*c+r];
        }
        auto material=std::make_unique<btDeformableLinearElasticityForce>();
        material->setYoungsModulus(first->young);material->setPoissonRatio(first->poisson);
        material->setDamping(s.nodes[begin].damping,first->damping);material->addSoftBody(body.get());
        s.world.addSoftBody(body.get());s.world.addForce(material.get());
        s.bodies.emplace(id,std::move(body));s.materials.emplace(id,std::move(material));
    }
    std::vector<int> local(count);
    for(int i=0;i<count;++i)
    {
        const auto& map=maps[i];auto& mesh=s.meshes[map.body];if(!mesh)mesh=std::make_unique<Mesh>();
        local[i]=int(mesh->sourceVertices.size());mesh->sourceVertices.push_back(i);
        btVector3 point=vector(map.offset);
        for(int j=0;j<4;++j)if(map.nodes[j]>=0)point+=vector(s.nodes[map.nodes[j]].position)*map.weights[j];
        for(int j=0;j<3;++j)mesh->points.push_back(point[j]);
        if(s.bodies.count(map.body))
        {
            if(map.offset.x||map.offset.y||map.offset.z) return s.fail("Bridge expects zero dynamic surface offsets");
            std::array<int,4> key;std::copy(map.nodes,map.nodes+4,key.begin());std::sort(key.begin(),key.end());
            auto found=s.tetLookup.find(key);
            if(found==s.tetLookup.end()||found->second.first!=map.body) return s.fail("Surface support does not identify a tetrahedron");
            auto& body=*s.bodies.at(map.body);const int begin=s.ranges.at(map.body).begin,index=found->second.second;
            btSoftBody::btVertexToTetraMapping converted{};converted.vertexToTetra=index;
            for(int j=0;j<4;++j)
            {
                int global=begin+int(body.m_tetras[index].m_n[j]-&body.m_nodes[0]);
                for(int k=0;k<4;++k)if(map.nodes[k]==global)converted.baryCoordInTetra[j]=map.weights[k];
            }
            body.maps.push_back(converted);
        }
        else for(int j=0;j<4;++j)if(map.nodes[j]>=0)return s.fail("Static surface has dynamic support");
    }
    for(int f=0;f<faceCount;++f)
    {
        if(faces[f].vertices[0]<0||faces[f].vertices[0]>=count)return s.fail("Invalid triangle vertex");
        const int id=maps[faces[f].vertices[0]].body;
        for(int j=0;j<3;++j)
        {
            int v=faces[f].vertices[j];if(v<0||v>=count||maps[v].body!=id)return s.fail("Triangle crosses body ranges");
            s.meshes.at(id)->indices.push_back(local[v]);
        }
    }
    for(auto& entry:s.meshes)
    {
        const int id=entry.first;auto& mesh=*entry.second;
        mesh.meshInterface=std::make_unique<btTriangleIndexVertexArray>(int(mesh.indices.size()/3),mesh.indices.data(),3*sizeof(int),int(mesh.points.size()/3),mesh.points.data(),3*sizeof(btScalar));
        auto* shape=new btGImpactMeshShape(mesh.meshInterface.get());shape->updateBound();
        if(s.bodies.count(id)){auto& body=*s.bodies.at(id);delete body.getCollisionShape();body.setCollisionShape(shape);}
        else
        {
            mesh.staticShape.reset(shape);mesh.obstacle=std::make_unique<btCollisionObject>();
            mesh.obstacle->setCollisionShape(shape);mesh.obstacle->setCollisionFlags(btCollisionObject::CF_STATIC_OBJECT);
            s.world.addCollisionObject(mesh.obstacle.get());
        }
    }
    s.world.setContactDiscoveryPadding(margin);s.world.setGravity(btVector3(0,0,0));s.cooking.reset();s.built=true;
    return s.updateTargets();
}
int forces(void* value,int count,const btGpuXpbdVec* forces)
{
    auto& s=*static_cast<State*>(value);if(count!=int(s.nodes.size())||!forces)return s.fail("Force count mismatch");
    s.forces.assign(forces,forces+count);return 1;
}
int targets(void* value,int count,const btGpuXpbdTarget* targets)
{
    auto& s=*static_cast<State*>(value);if(count<0||(count&&!targets))return s.fail("Invalid target input");
    if(count)s.targets.assign(targets,targets+count);else s.targets.clear();
    return s.built?s.updateTargets():1;
}
int step(void* value,const btGpuXpbdStep* config)
{
    auto& s=*static_cast<State*>(value);if(!s.built||!config||config->velocityIterations!=1)return s.fail("Invalid world bridge step");
    s.world.setGpuXpbdTgsSolver(true,.001,config->substeps,config->iterations,0);
    s.world.setGravity(vector(config->gravity));
    if(!s.world.setGpuXpbdTgsSettling(nullptr,0,0))return s.fail(s.world.getGpuXpbdTgsError().c_str());
    if(config->settlingSpeed||config->settlingDamping)
    {
        const int id=s.owner(config->settlingBegin);
        if(id<0||s.ranges.at(id).begin!=config->settlingBegin||s.ranges.at(id).end!=config->settlingEnd)return s.fail("Settling range must identify one body");
        if(!s.world.setGpuXpbdTgsSettling(s.bodies.at(id).get(),config->settlingSpeed,config->settlingDamping))return s.fail(s.world.getGpuXpbdTgsError().c_str());
    }
    for(const auto& entry:s.bodies)
    {
        const int begin=s.ranges.at(entry.first).begin;auto& body=*entry.second;
        for(int n=0;n<body.m_nodes.size();++n)body.m_nodes[n].m_f=vector(s.forces[begin+n]);
    }
    if(s.world.stepSimulation(config->dt,0)!=1||s.world.hasDeformableStepFailed())return s.fail(s.world.getGpuXpbdTgsError().c_str());
    for(const auto& entry:s.bodies)
    {
        const int begin=s.ranges.at(entry.first).begin;const auto& body=*entry.second;
        for(int n=0;n<body.m_nodes.size();++n){s.nodes[begin+n].position=vector(body.m_nodes[n].m_x);s.nodes[begin+n].velocity=vector(body.m_nodes[n].m_v);}
    }
    return 1;
}
int download(void* value,int count,btGpuXpbdNode* nodes)
{
    auto& s=*static_cast<State*>(value);if(count!=int(s.nodes.size())||!nodes)return s.fail("Download count mismatch");
    std::copy(s.nodes.begin(),s.nodes.end(),nodes);return 1;
}
int surfacePositions(void* value,int count,btGpuXpbdVec* output)
{
    // Diagnostic trajectory export only; collision skinning and all physics run on the GPU world.
    auto& s=*static_cast<State*>(value);if(count!=int(s.maps.size())||!output)return s.fail("Surface download count mismatch");
    for(int i=0;i<count;++i)
    {
        const auto& map=s.maps[i];btVector3 p=vector(map.offset);
        for(int j=0;j<4;++j)if(map.nodes[j]>=0)p+=vector(s.nodes[map.nodes[j]].position)*map.weights[j];
        output[i]=vector(p);
    }
    return 1;
}
int disabledPairs(void* value,int count,const btGpuXpbdDisabledPair*)
{
    auto& s=*static_cast<State*>(value);
    return count==0?1:s.fail("Test bridge does not expose pair filtering");
}
int updateNodes(void* value,int count,const btGpuXpbdNode* nodes)
{
    auto& s=*static_cast<State*>(value);
    if(!s.built||count!=int(s.nodes.size())||!nodes||!s.targets.empty())return s.fail("Invalid bridge node update or active targets");
    for(int i=0;i<count;++i)
    {
        const auto& n=nodes[i];
        if(!std::isfinite(n.position.x)||!std::isfinite(n.position.y)||!std::isfinite(n.position.z)||
           !std::isfinite(n.velocity.x)||!std::isfinite(n.velocity.y)||!std::isfinite(n.velocity.z)||
           n.inverseMass!=s.nodes[i].inverseMass||n.damping!=s.nodes[i].damping)
            return s.fail("Bridge checkpoint requires finite nodes with unchanged masses and damping");
    }
    s.nodes.assign(nodes,nodes+count);
    for(auto& entry:s.bodies)
    {
        const int begin=s.ranges.at(entry.first).begin;auto& body=*entry.second;
        for(int n=0;n<body.m_nodes.size();++n)
        {
            auto& node=body.m_nodes[n];node.m_x=vector(nodes[begin+n].position);node.m_q=node.m_x;
            node.m_v=vector(nodes[begin+n].velocity);node.m_vn=node.m_v;
        }
    }
    return 1;
}
const btGpuXpbdTgsApi api={create,destroy,upload,forces,step,download,error,surface,surfacePositions,embed,targets,disabledPairs,updateNodes};
}
extern "C" __declspec(dllexport) const btGpuXpbdTgsApi* btGetGpuXpbdTgsApi(int version){return version>=1&&version<=3?&api:nullptr;}
