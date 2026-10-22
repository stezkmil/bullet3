// GPU XPBD TGS world adapter. All iterative mechanics and contact work stays in the backend.
#include "btDeformableMultiBodyDynamicsWorld.h"
#include "btDeformableBodySolver.h"
#include "btDeformableLinearElasticityForce.h"
#include "btDeformableNodalForce.h"
#include "btDeformableMousePickingForce.h"
#include "btDeformableDiagnostics.h"
#include "btGpuXpbdTgsBackend.h"
#include "BulletCollision/Gimpact/btGImpactShape.h"
#include "BulletCollision/CollisionShapes/btCompoundShape.h"
#include "BulletCollision/CollisionShapes/btConcaveShape.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <map>
#include <queue>
#include <cfloat>
#include "BulletCollision/BroadphaseCollision/btDbvt.h"
#include "BulletCollision/NarrowPhaseCollision/btVoronoiSimplexSolver.h"

namespace
{
btGpuXpbdVec gpuVector(const btVector3& p, btScalar scale)
{
    return {double(p.x()*scale),double(p.y()*scale),double(p.z()*scale)};
}
btVector3 bulletVector(btGpuXpbdVec p, btScalar scale)
{
    return btVector3(p.x,p.y,p.z)/scale;
}
template<class T> void appendKey(std::vector<unsigned char>& key, const T& value)
{
    const auto* bytes=reinterpret_cast<const unsigned char*>(&value);
    key.insert(key.end(),bytes,bytes+sizeof(value));
}
void appendTransform(std::vector<unsigned char>& key,const btTransform& transform)
{
    for(int i=0;i<3;++i)
    {
        appendKey(key,transform.getOrigin()[i]);
        for(int j=0;j<3;++j) appendKey(key,transform.getBasis()[i][j]);
    }
}
bool identityBasis(const btMatrix3x3& basis)
{
    for(int i=0;i<3;++i) for(int j=0;j<3;++j)
        if(btFabs(basis[i][j]-btScalar(i==j))>btScalar(1.e-7)) return false;
    return true;
}
void shapeKey(std::vector<unsigned char>& key,const btCollisionShape* shape)
{
    appendKey(key,shape);
    for(int i=0;i<3;++i) appendKey(key,shape->getLocalScaling()[i]);
    if(shape->isCompound())
    {
        const auto* compound=static_cast<const btCompoundShape*>(shape);
        appendKey(key,compound->getNumChildShapes());
        for(int i=0;i<compound->getNumChildShapes();++i)
        {
            appendTransform(key,compound->getChildTransform(i));
            shapeKey(key,compound->getChildShape(i));
        }
    }
}
}

struct btDeformableGpuXpbdTgsWorld
{
    btGpuXpbdTgsBackend backend;
    std::vector<unsigned char> key;
    std::vector<btGpuXpbdNode> accepted;
    std::vector<btGpuXpbdDisabledPair> disabled;
    std::vector<btGpuXpbdSurfaceVertex> vertices;
    std::vector<btGpuXpbdSurfaceTriangle> triangles;
    std::string error;
    bool hasTargets=false;
    btVector3 origin=btVector3(0,0,0);
    bool hasOrigin=false;

    bool fail(const char* message) { error=message;return false; }
    bool checked(int result) { if(!result) error=backend.error();return result!=0; }

    void staticTriangle(const btVector3* points,const btTransform& transform,btScalar scale,int body)
    {
        btGpuXpbdSurfaceTriangle triangle{};
        for(int i=0;i<3;++i)
        {
            triangle.vertices[i]=int(vertices.size());
            vertices.push_back({{-1,-1,-1,-1},{0,0,0,0},gpuVector(transform*points[i]-origin,scale),body});
        }
        triangles.push_back(triangle);
    }

    bool mesh(btGImpactMeshShape* shape,btSoftBody* soft,int offset,const btTransform& transform,btScalar scale,int body)
    {
        const auto* mapping=soft?soft->getCollisionShapeVertexToSimTetra():nullptr;
        if(soft&&!mapping) return fail("GPU XPBD mapped mesh has no tetrahedron embedding");
        if(soft&&!identityBasis(transform.getBasis())) return fail("GPU XPBD requires world-space simulation nodes and collision embedding");
        for(int part=0;part<shape->getMeshPartCount();++part)
        {
            auto* piece=shape->getMeshPart(part);
            if(soft&&(piece->getLocalScaling()-btVector3(1,1,1)).length2()>btScalar(1.e-14))
                return fail("GPU XPBD requires collision scaling baked into simulation nodes");
            piece->lockChildShapes();
            struct Unlock { btGImpactMeshShapePart* piece;~Unlock(){piece->unlockChildShapes();} } unlock{piece};
            const auto* manager=piece->getTrimeshPrimitiveManager();
            std::map<unsigned int,int> indices;
            for(int face=0;face<manager->get_primitive_count();++face)
            {
                unsigned int ids[3];manager->get_primitive_indices(face,ids[0],ids[1],ids[2]);
                btGpuXpbdSurfaceTriangle triangle{};
                for(int c=0;c<3;++c)
                {
                    auto found=indices.find(ids[c]);
                    if(found==indices.end())
                    {
                        btGpuXpbdSurfaceVertex vertex={{-1,-1,-1,-1},{0,0,0,0},{0,0,0},body};
                        if(soft)
                        {
                            if(ids[c]>=mapping->size()) return fail("GPU XPBD collision embedding index is invalid");
                            const auto& entry=(*mapping)[ids[c]];
                            if(entry.vertexToTetra>=unsigned(soft->m_tetras.size())) return fail("GPU XPBD collision tetrahedron is invalid");
                            const auto& tetra=soft->m_tetras[entry.vertexToTetra];
                            btScalar weightSum=0;
                            for(int j=0;j<4;++j)
                            {
                                vertex.nodes[j]=offset+int(tetra.m_n[j]-&soft->m_nodes[0]);
                                vertex.weights[j]=double(entry.baryCoordInTetra[j]);
                                weightSum+=btScalar(vertex.weights[j]);
                            }
                            vertex.offset=gpuVector(transform.getOrigin()+origin*(weightSum-1),scale);
                        }
                        else
                        {
                            btVector3 point;manager->get_vertex(ids[c],point,false);
                            vertex.offset=gpuVector(transform*point-origin,scale);
                        }
                        const int id=int(vertices.size());vertices.push_back(vertex);
                        found=indices.emplace(ids[c],id).first;
                    }
                    triangle.vertices[c]=found->second;
                }
                triangles.push_back(triangle);
            }
        }
        return true;
    }

    bool obstacle(btCollisionShape* shape,const btTransform& transform,btScalar scale,int body)
    {
        if(shape->isCompound())
        {
            auto* compound=static_cast<btCompoundShape*>(shape);
            for(int i=0;i<compound->getNumChildShapes();++i)
                if(!obstacle(compound->getChildShape(i),transform*compound->getChildTransform(i),scale,body)) return false;
            return true;
        }
        if(shape->getShapeType()==GIMPACT_SHAPE_PROXYTYPE)
        {
            auto* impact=static_cast<btGImpactShapeInterface*>(shape);
            if(impact->getGImpactShapeType()!=CONST_GIMPACT_TRIMESH_SHAPE) return fail("GPU XPBD requires triangle mesh obstacles");
            return mesh(static_cast<btGImpactMeshShape*>(shape),nullptr,0,transform,scale,body);
        }
        if(!shape->isConcave()) return fail("GPU XPBD requires triangle mesh static obstacles");
        struct Collector:btTriangleCallback
        {
            btDeformableGpuXpbdTgsWorld& output;const btTransform& transform;btScalar scale;int body;
            Collector(btDeformableGpuXpbdTgsWorld& o,const btTransform& t,btScalar s,int b):output(o),transform(t),scale(s),body(b){}
            void processTriangle(btVector3* p,int,int) override { output.staticTriangle(p,transform,scale,body); }
        } collector(*this,transform,scale,body);
        static_cast<btConcaveShape*>(shape)->processAllTriangles(&collector,btVector3(-BT_LARGE_FLOAT,-BT_LARGE_FLOAT,-BT_LARGE_FLOAT),btVector3(BT_LARGE_FLOAT,BT_LARGE_FLOAT,BT_LARGE_FLOAT));
        return true;
    }

    void boundarySamples(btSoftBody* soft,int offset,int body,const std::vector<btGpuXpbdNode>& rest,int firstTriangle)
    {
        if(firstTriangle>=int(triangles.size())) return;
        std::map<std::array<int,3>,int> faces;
        for(int t=0;t<soft->m_tetras.size();++t) for(int opposite=0;opposite<4;++opposite)
        {
            std::array<int,3> ids;int k=0;
            for(int j=0;j<4;++j) if(j!=opposite) ids[k++]=offset+int(soft->m_tetras[t].m_n[j]-&soft->m_nodes[0]);
            std::sort(ids.begin(),ids.end());++faces[ids];
        }
        std::vector<bool> boundaryNodes(soft->m_nodes.size(),false),represented(soft->m_nodes.size(),false);
        for(const auto& face:faces) if(face.second==1) for(int n:face.first) boundaryNodes[n-offset]=true;
        btVector3 low=bulletVector(rest[offset].position,1),high=low;
        for(int n=0;n<soft->m_nodes.size();++n) { auto p=bulletVector(rest[offset+n].position,1);low.setMin(p);high.setMax(p); }
        const btScalar tolerance=32*FLT_EPSILON*(high-low).length();
        std::vector<std::array<btVector3,3>> geometry(triangles.size()-firstTriangle);
        btDbvt tree;
        for(int i=firstTriangle;i<int(triangles.size());++i)
        {
            auto& points=geometry[i-firstTriangle];
            for(int c=0;c<3;++c)
            {
                const auto& v=vertices[triangles[i].vertices[c]];points[c]=bulletVector(v.offset,1);
                for(int j=0;j<4;++j) if(v.nodes[j]>=0) points[c]+=bulletVector(rest[v.nodes[j]].position,1)*v.weights[j];
                for(int j=0;j<4;++j) if(v.nodes[j]>=offset&&v.nodes[j]<offset+soft->m_nodes.size()&&
                    (points[c]-bulletVector(rest[v.nodes[j]].position,1)).length2()<=tolerance*tolerance&&v.weights[j]>1-8*FLT_EPSILON)
                    represented[v.nodes[j]-offset]=true;
            }
            btVector3 a=points[0],b=a;for(int j=1;j<3;++j) { a.setMin(points[j]);b.setMax(points[j]); }
            tree.insert(btDbvtVolume::FromMM(a,b),&points);
        }
        struct Query:btDbvt::ICollide
        {
            btVector3 point;btScalar tolerance2;bool found=false;btVoronoiSimplexSolver solver;
            Query(const btVector3& p,btScalar t):point(p),tolerance2(t*t){}
            void Process(const btDbvtNode* node)
            {
                if(found) return;
                const auto& t=*static_cast<const std::array<btVector3,3>*>(node->data);
                if((t[1]-t[0]).cross(t[2]-t[0]).length2()==0) return;
                btSubSimplexClosestResult result;
                solver.closestPtPointTriangle(point,t[0],t[1],t[2],result);
                found=(point-result.m_closestPointOnSimplex).length2()<=tolerance2;
            }
        };
        const btVector3 padding(tolerance,tolerance,tolerance);
        for(int n=0;n<soft->m_nodes.size();++n) if(boundaryNodes[n]&&!represented[n])
        {
            auto p=bulletVector(rest[offset+n].position,1);Query query(p,tolerance);
            tree.collideTV(tree.m_root,btDbvtVolume::FromMM(p-padding,p+padding),query);
            // Sample only coincident material nodes; the collision mesh still defines the shape.
            if(query.found) vertices.push_back({{offset+n,-1,-1,-1},{1,0,0,0},{0,0,0},body});
        }
    }
    void boundary(btSoftBody* soft,int offset,int body)
    {
        const int first=int(vertices.size());
        for(int n=0;n<soft->m_nodes.size();++n) vertices.push_back({{offset+n,-1,-1,-1},{1,0,0,0},{0,0,0},body});
        struct Face { std::array<int,3> oriented;int count=0; };
        std::map<std::array<int,3>,Face> faces;
        for(int t=0;t<soft->m_tetras.size();++t) for(int opposite=0;opposite<4;++opposite)
        {
            const auto& tetra=soft->m_tetras[t];std::array<int,3> ids,local;int k=0;
            for(int j=0;j<4;++j) if(j!=opposite) { local[k]=j;ids[k++]=int(tetra.m_n[j]-&soft->m_nodes[0]); }
            // Preserve rest winding when the solver is enabled in a deformed pose.
            const auto dm=tetra.m_Dm_inverse.inverse();
            const btVector3 rest[]={btVector3(0,0,0),dm.getColumn(0),dm.getColumn(1),dm.getColumn(2)};
            const auto& a=rest[local[0]];const auto& b=rest[local[1]];const auto& c=rest[local[2]];
            if((b-a).cross(c-a).dot(rest[opposite]-a)>0) std::swap(ids[1],ids[2]);
            auto key=ids;std::sort(key.begin(),key.end());auto& face=faces[key];face.oriented=ids;++face.count;
        }
        for(const auto& entry:faces) if(entry.second.count==1)
        {
            const auto& ids=entry.second.oriented;triangles.push_back({{first+ids[0],first+ids[1],first+ids[2]}});
        }
    }
};

bool btDeformableMultiBodyDynamicsWorld::setGpuXpbdTgsKinematicTargets(const std::vector<btGpuXpbdTgsWorldTarget>& targets)
{
    std::map<std::pair<int,int>,bool> used;
    for(const auto& target:targets)
    {
        const int index=m_softBodies.findLinearSearch(target.body);
        if(index>=m_softBodies.size()||target.node<0||target.node>=target.body->m_nodes.size()||target.body->isStaticObject()||
           !std::isfinite(double(target.position.x()))||!std::isfinite(double(target.position.y()))||!std::isfinite(double(target.position.z()))||
           !used.emplace(std::make_pair(index,target.node),true).second)
        {
            m_gpuXpbdTgsError="Invalid or duplicate GPU XPBD kinematic target";return false;
        }
    }
    m_gpuXpbdTgsTargets=targets;return true;
}

void btDeformableMultiBodyDynamicsWorld::clearGpuXpbdTgsKinematicTargets(btSoftBody* body)
{
    auto& targets=m_gpuXpbdTgsTargets;
    targets.erase(std::remove_if(targets.begin(),targets.end(),[&](const auto& target){return !body||target.body==body;}),targets.end());
}

bool btDeformableMultiBodyDynamicsWorld::setGpuXpbdTgsSettling(btSoftBody* body,btScalar speed,btScalar damping)
{
    if(!std::isfinite(double(speed))||!std::isfinite(double(damping))||speed<0||damping<0||
       (body&&(m_softBodies.findLinearSearch(body)>=m_softBodies.size()||body->isStaticObject()))||(!body&&(speed!=0||damping!=0)))
    {
        m_gpuXpbdTgsError="Invalid GPU XPBD settling body or parameters";return false;
    }
    m_gpuXpbdTgsSettlingBody=body;m_gpuXpbdTgsSettlingSpeed=speed;m_gpuXpbdTgsSettlingDamping=damping;return true;
}

void btDeformableMultiBodyDynamicsWorld::setGpuXpbdTgsSolver(bool enabled,btScalar metersPerUnit,int substeps,int iterations,btScalar maxMassRatio,int materialMode)
{
    if(m_gpuXpbdTgsEnabled!=enabled||m_gpuXpbdTgsMetersPerUnit!=metersPerUnit)
    {
        if(m_gpuXpbdTgsEnabled&&!enabled)
            for(int b=0;b<m_softBodies.size();++b)
            {
                auto* body=m_softBodies[b];
                body->m_tetraScratches.resize(body->m_tetras.size());
                body->m_tetraScratchesTn.resize(body->m_tetras.size());
                body->advanceDeformation();
                body->updateLastSafeWorldTransform();
            }
        m_gpuXpbdTgsWorld.reset();m_gpuXpbdTgsError.clear();m_deformableStepFailed=false;
    }
    m_gpuXpbdTgsEnabled=enabled;m_gpuXpbdTgsMetersPerUnit=metersPerUnit;
    m_gpuXpbdTgsSubsteps=substeps;m_gpuXpbdTgsIterations=iterations;
    m_gpuXpbdTgsMaxMassRatio=maxMassRatio;
    m_gpuXpbdTgsMaterialMode=materialMode;
}

void btDeformableMultiBodyDynamicsWorld::gpuXpbdTgsSingleStepSimulation(btScalar timeStep)
{
    auto fail=[&](const std::string& reason)
    {
        m_deformableStepFailed=true;m_gpuXpbdTgsError=reason;
        btDeformableDiagnostics::write("STEP_FAILED","solver=gpu_xpbd_tgs reason=%s",reason.c_str());
    };
    const btScalar units=m_gpuXpbdTgsMetersPerUnit*1000;
    if(!std::isfinite(double(units))||units<=0||!std::isfinite(double(timeStep))||timeStep<=0)
        { fail("Invalid GPU XPBD units or time step");return; }
    if(!std::isfinite(double(m_gpuXpbdTgsMaxMassRatio))||m_gpuXpbdTgsMaxMassRatio<0||
       (m_gpuXpbdTgsMaxMassRatio>0&&m_gpuXpbdTgsMaxMassRatio<1))
        { fail("GPU XPBD mass ratio must be zero or at least one");return; }
    if(m_multiBodies.size()||m_constraints.size()||m_solverCallback)
        { fail("GPU XPBD articulated bodies, rigid joints and solver callbacks are not implemented");return; }
    if(m_internalPreTickCallback) (*m_internalPreTickCallback)(this,timeStep);
    if(!m_gpuXpbdTgsWorld) m_gpuXpbdTgsWorld=std::make_shared<btDeformableGpuXpbdTgsWorld>();
    auto& cache=*m_gpuXpbdTgsWorld;
    if(!cache.backend.ready()) { fail(cache.backend.error());return; }
    const auto* api=cache.backend.api();void* state=cache.backend.state();
    std::vector<btSoftBody*> bodies;std::vector<int> offsets,owners;
    std::vector<btGpuXpbdNode> nodes;std::vector<btGpuXpbdTet> tets;std::vector<btGpuXpbdVec> forces;
    std::vector<btGpuXpbdSpring> springs;
    std::vector<unsigned char> key;appendKey(key,units);
    for(int b=0;b<m_softBodies.size();++b)
    {
        auto* body=m_softBodies[b];if(body->isStaticObject()) continue;
        // Keep solver coordinates local by subtracting the origin in host precision.
        if(!cache.hasOrigin&&body->m_nodes.size())
        {
            cache.origin=body->m_nodes[0].m_x;cache.hasOrigin=true;
        }
        if(!body->m_tetras.size()||body->m_anchors.size()||body->m_deformableAnchors.size()||
           (body->m_cfg.collisions&btSoftBody::fCollision::CL_SELF))
            { fail("GPU XPBD currently requires tetrahedral bodies without anchors or self-collision");return; }
        int owner=m_collisionObjects.findLinearSearch(body);
        if(owner>=m_collisionObjects.size()) { fail("GPU XPBD body is absent from the world");return; }
        const int offset=int(nodes.size());bodies.push_back(body);offsets.push_back(offset);owners.push_back(owner);
        btDeformableLinearElasticityForce* material=nullptr;
        std::vector<btVector3> external(body->m_nodes.size(),btVector3(0,0,0));
        auto& objective=*m_deformableBodySolver->m_objective;
        for(int f=0;f<objective.m_lf.size();++f)
        {
            auto* force=objective.m_lf[f];bool attached=false;
            for(int j=0;j<force->m_softBodies.size();++j) attached|=force->m_softBodies[j]==body;
            if(!attached) continue;
            if(force->getForceType()==BT_LINEAR_ELASTICITY_FORCE)
            {
                if(material) { fail("Multiple GPU XPBD elastic materials on one body");return; }
                material=static_cast<btDeformableLinearElasticityForce*>(force);
            }
            else if(force->getForceType()==BT_NODAL_FORCE)
            {
                auto* nodal=static_cast<btDeformableNodalForce*>(force);int count=0;
                for(int j=0;j<nodal->nodeIndices().size();++j)
                {
                    int n=nodal->nodeIndices()[j];
                    if(n>=0&&n<body->m_nodes.size()&&body->m_nodes[n].m_im>0&&body->m_nodes[n].m_frozen<=0) ++count;
                }
                if(count&&body->isActive()) for(int j=0;j<nodal->nodeIndices().size();++j)
                {
                    int n=nodal->nodeIndices()[j];
                    if(n>=0&&n<body->m_nodes.size()&&body->m_nodes[n].m_im>0&&body->m_nodes[n].m_frozen<=0)
                        external[n]+=nodal->totalForce()/count;
                }
            }
            else if(force->getForceType()==BT_MOUSE_PICKING_FORCE)
            {
                std::vector<btDeformableMousePickingForce::NodeSpring> picked;
                static_cast<btDeformableMousePickingForce*>(force)->appendNodeSprings(picked);
                for(const auto& spring:picked)
                {
                    int n=0;while(n<body->m_nodes.size()&&spring.node!=&body->m_nodes[n]) ++n;
                    if(n==body->m_nodes.size()) continue;
                    springs.push_back({offset+n,gpuVector(spring.target-cache.origin,units),
                        double(spring.stiffness),double(spring.damping),double(spring.maxForce*units)});
                }
            }
            else if(force->getForceType()!=BT_GRAVITY_FORCE) { fail("Unsupported GPU XPBD force type");return; }
        }
        if(!material||material->m_mu<=0||material->m_lambda<0) { fail("GPU XPBD requires a valid linear elastic material");return; }
        const btScalar mu=material->m_mu,lambda=material->m_lambda;
        btScalar minInverseMass=SIMD_INFINITY;
        for(int n=0;n<body->m_nodes.size();++n)
            if(body->m_nodes[n].m_im>0) minInverseMass=btMin(minInverseMass,body->m_nodes[n].m_im);
        // PhysX setMass raises the lightest nodal masses without renormalizing the total.
        const btScalar maxInverseMass=m_gpuXpbdTgsMaxMassRatio>0?minInverseMass*m_gpuXpbdTgsMaxMassRatio:SIMD_INFINITY;
        for(int n=0;n<body->m_nodes.size();++n)
        {
            const auto& node=body->m_nodes[n];const btScalar inv=node.m_frozen>0?0:btMin(node.m_im,maxInverseMass);
            nodes.push_back({gpuVector(node.m_x-cache.origin,units),double(inv),gpuVector(inv>0?node.m_v:btVector3(0,0,0),units),double(material->m_damping_alpha)});
            const auto gravity=inv>0?getGravity()*(body->m_gravityFactor/inv):btVector3(0,0,0);
            forces.push_back(gpuVector(external[n]+node.m_f+gravity,units));
        }
        appendKey(key,body);appendKey(key,body->m_nodes.size());appendKey(key,body->getCollisionShapeVertexToSimTetra());
        for(int i=0;i<body->m_tetras.size();++i)
        {
            const auto& source=body->m_tetras[i];btGpuXpbdTet tet{};
            for(int j=0;j<4;++j) tet.nodes[j]=offset+int(source.m_n[j]-&body->m_nodes[0]);
            for(int c=0;c<3;++c) for(int r=0;r<3;++r) tet.inverseRest[3*c+r]=double(source.m_Dm_inverse[r][c]/units);
            tet.volume=double(source.m_element_measure*units*units*units);
            tet.young=double(mu*(3*lambda+2*mu)/((lambda+mu)*units));tet.poisson=double(lambda/(2*(lambda+mu)));
            tet.damping=double(material->m_damping_beta);tets.push_back(tet);appendKey(key,tet);
        }
    }
    if(nodes.empty()) { fail("GPU XPBD requires at least one dynamic tetrahedral body");return; }
    std::vector<btGpuXpbdTarget> targets;
    for(const auto& target:m_gpuXpbdTgsTargets)
    {
        auto found=std::find(bodies.begin(),bodies.end(),target.body);
        if(found==bodies.end()||target.node<0||target.node>=target.body->m_nodes.size())
            { fail("GPU XPBD kinematic target body or node is no longer available");return; }
        targets.push_back({offsets[found-bodies.begin()]+target.node,gpuVector(target.position-cache.origin,units)});
    }
    std::sort(targets.begin(),targets.end(),[](const auto& a,const auto& b){return a.node<b.node;});
    std::vector<double> inverseMasses;inverseMasses.reserve(nodes.size());
    for(const auto& node:nodes) inverseMasses.push_back(node.inverseMass);
    std::vector<double> friction;
    for(int i=0;i<m_collisionObjects.size();++i)
    {
        const auto* object=m_collisionObjects[i];appendKey(key,object);appendTransform(key,object->getWorldTransform());
        const btScalar coefficient=object->getFriction();
        if(!std::isfinite(double(coefficient))||coefficient<0) { fail("GPU XPBD friction must be finite and nonnegative");return; }
        friction.push_back(double(coefficient));appendKey(key,friction.back());
        shapeKey(key,object->getCollisionShape());
        if(!object->isStaticObject()&&std::find(bodies.begin(),bodies.end(),object)==bodies.end())
            { fail("GPU XPBD moving rigid obstacles are not implemented");return; }
    }
    const double margin=double(m_contactDiscoveryPadding>=0?btMax(btScalar(.0001),m_contactDiscoveryPadding*units):btScalar(.1));
    appendKey(key,margin);
    if(cache.key!=key)
    {
        cache.key.clear();cache.accepted.clear();cache.disabled.clear();cache.vertices.clear();cache.triangles.clear();
        // Reconstruct the undeformed geometry from element rest edges, even after a reset to a deformed pose.
        struct Edge { int node;btVector3 delta; };
        std::vector<std::vector<Edge>> edges(nodes.size());
        for(int b=0;b<int(bodies.size());++b) for(int t=0;t<bodies[b]->m_tetras.size();++t)
        {
            const auto& tetra=bodies[b]->m_tetras[t];auto rest=tetra.m_Dm_inverse.inverse();
            int a=offsets[b]+int(tetra.m_n[0]-&bodies[b]->m_nodes[0]);
            for(int j=1;j<4;++j)
            {
                int c=offsets[b]+int(tetra.m_n[j]-&bodies[b]->m_nodes[0]);auto delta=rest.getColumn(j-1)*units;
                edges[a].push_back({c,delta});edges[c].push_back({a,-delta});
            }
        }
        auto rest=nodes;std::vector<bool> visited(nodes.size(),false);std::queue<int> pending;
        for(int seed=0;seed<int(nodes.size());++seed) if(!visited[seed])
        {
            visited[seed]=true;pending.push(seed);
            while(!pending.empty())
            {
                int a=pending.front();pending.pop();
                for(const auto& edge:edges[a]) if(!visited[edge.node])
                {
                    rest[edge.node].position=gpuVector(bulletVector(rest[a].position,1)+edge.delta,1);
                    visited[edge.node]=true;pending.push(edge.node);
                }
            }
        }
        if(!cache.checked(api->upload(state,int(rest.size()),rest.data(),int(tets.size()),tets.data()))||
           !cache.checked(api->updateNodes(state,int(nodes.size()),nodes.data()))) { fail(cache.error);return; }
        cache.hasTargets=false;
        for(int b=0;b<int(bodies.size());++b)
        {
            auto* body=bodies[b];auto* shape=body->getCollisionShape();
            if(body->getCollisionShapeVertexToSimTetra())
            {
                if(shape->getShapeType()!=GIMPACT_SHAPE_PROXYTYPE||static_cast<btGImpactShapeInterface*>(shape)->getGImpactShapeType()!=CONST_GIMPACT_TRIMESH_SHAPE)
                    { fail("GPU XPBD mapped collision surface requires a GImpact triangle mesh");return; }
                const int firstTriangle=int(cache.triangles.size());
                if(!cache.mesh(static_cast<btGImpactMeshShape*>(shape),body,offsets[b],body->getWorldTransform(),units,owners[b])) { fail(cache.error);return; }
                cache.boundarySamples(body,offsets[b],owners[b],rest,firstTriangle);
            }
            else cache.boundary(body,offsets[b],owners[b]);
        }
        for(int i=0;i<m_collisionObjects.size();++i) if(m_collisionObjects[i]->isStaticObject())
            if(!cache.obstacle(m_collisionObjects[i]->getCollisionShape(),m_collisionObjects[i]->getWorldTransform(),units,i)) { fail(cache.error);return; }
        if(!cache.checked(api->surface(state,int(cache.vertices.size()),cache.vertices.data(),int(cache.triangles.size()),cache.triangles.data(),margin,0)))
            { fail(cache.error);return; }
        if(!cache.checked(api->friction(state,int(friction.size()),friction.data()))) { fail(cache.error);return; }
        cache.key=key;cache.accepted=nodes;
    }
    else if(cache.accepted.size()!=nodes.size()||std::memcmp(nodes.data(),cache.accepted.data(),nodes.size()*sizeof(nodes[0]))!=0)
    {
        if(cache.hasTargets)
        {
            if(!cache.checked(api->targets(state,0,nullptr))) { fail(cache.error);return; }
            cache.hasTargets=false;
        }
        if(!cache.checked(api->updateNodes(state,int(nodes.size()),nodes.data()))) { fail(cache.error);return; }
    }
    std::vector<btGpuXpbdDisabledPair> disabled;
    for(int a=0;a<m_collisionObjects.size();++a) for(int b=a+1;b<m_collisionObjects.size();++b)
    {
        auto* first=m_collisionObjects[a];auto* second=m_collisionObjects[b];
        const auto* pa=first->getBroadphaseHandle();const auto* pb=second->getBroadphaseHandle();
        if(!first->hasContactResponse()||!second->hasContactResponse()||!pa||!pb||!(pa->m_collisionFilterGroup&pb->m_collisionFilterMask)||!(pb->m_collisionFilterGroup&pa->m_collisionFilterMask)||
           !first->checkCollideWith(second)||!second->checkCollideWith(first)) disabled.push_back({a,b});
    }
    if(disabled.size()!=cache.disabled.size()||(!disabled.empty()&&std::memcmp(disabled.data(),cache.disabled.data(),disabled.size()*sizeof(disabled[0]))!=0))
    {
        if(!cache.checked(api->disabledPairs(state,int(disabled.size()),disabled.data()))) { fail(cache.error);return; }
        cache.disabled=disabled;
    }
    btGpuXpbdStep config={double(timeStep),m_gpuXpbdTgsSubsteps,{0,0,0},0,0,m_gpuXpbdTgsIterations,0,0,1,m_gpuXpbdTgsMaterialMode};
    if(m_gpuXpbdTgsSettlingBody)
    {
        auto found=std::find(bodies.begin(),bodies.end(),m_gpuXpbdTgsSettlingBody);
        if(found==bodies.end()) { fail("GPU XPBD settling body is no longer available");return; }
        config.settlingBegin=offsets[found-bodies.begin()];
        config.settlingEnd=config.settlingBegin+m_gpuXpbdTgsSettlingBody->m_nodes.size();
        config.settlingSpeed=double(m_gpuXpbdTgsSettlingSpeed*units);config.settlingDamping=double(m_gpuXpbdTgsSettlingDamping);
    }
    if(!cache.checked(api->targets(state,int(targets.size()),targets.data()))) { fail(cache.error);return; }
    cache.hasTargets=!targets.empty();
    if(!cache.checked(api->springs(state,int(springs.size()),springs.data()))) { fail(cache.error);return; }
    if(!cache.checked(api->forces(state,int(forces.size()),forces.data()))||!cache.checked(api->step(state,&config))||
       !cache.checked(api->download(state,int(nodes.size()),nodes.data()))) { fail(cache.error);return; }
    cache.accepted=nodes;
    // Target kernels pin their nodes, but subsequent host comparisons use the original effective masses.
    for(int n=0;n<int(nodes.size());++n)
    {
        cache.accepted[n].inverseMass=inverseMasses[n];
        if(inverseMasses[n]<=0) cache.accepted[n].velocity={0,0,0};
    }
    for(int b=0;b<int(bodies.size());++b)
    {
        auto* body=bodies[b];
        for(int n=0;n<body->m_nodes.size();++n)
        {
            auto& node=body->m_nodes[n];const auto& result=nodes[offsets[b]+n];
            node.m_x=bulletVector(result.position,units)+cache.origin;node.m_q=node.m_x;
            node.m_v=bulletVector(result.velocity,units);node.m_vn=node.m_v;node.m_splitv.setZero();node.m_f.setZero();
            auto& accepted=cache.accepted[offsets[b]+n];
            accepted.position=gpuVector(node.m_x-cache.origin,units);
            accepted.velocity=gpuVector(accepted.inverseMass>0?node.m_v:btVector3(0,0,0),units);
        }
        // Only the host display and broadphase mirror need refreshing here.
        body->updateNormals();body->updateBounds();updateSingleAabb(body);
    }
    m_gpuXpbdTgsError.clear();m_internalTime+=timeStep;m_solverInfo.m_timeStep=timeStep;m_dispatchInfo.m_timeStep=timeStep;
    if(m_internalTickCallback) (*m_internalTickCallback)(this,timeStep);
    ++m_executed_step_counter;m_dispatchInfo.m_stepCounter=m_executed_step_counter;
}
