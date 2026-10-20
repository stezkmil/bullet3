#include "btGpuXpbdTgsApi.h"
#include "xpbd_ptx.h"
#include <cuda.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <functional>
#include <numeric>
#include <string>
#include <vector>
#include <array>
#include <map>

namespace
{
struct State
{
    CUcontext context=nullptr;
    CUmodule module=nullptr;
    CUstream stream=nullptr;
    CUgraph graph=nullptr;
    CUgraphExec graphExec=nullptr;
    btGpuXpbdStep graphConfig={};
    CUfunction predict=nullptr,solve=nullptr,settle=nullptr;
    CUfunction initializeCopies=nullptr,averageCopies=nullptr;
    CUdeviceptr copies=0,copyActive=0,solverTets=0;
    int copyLanes=0;
    CUfunction skin=nullptr,refit=nullptr,discoverContacts=nullptr,countContacts=nullptr,solveContacts=nullptr,applyContacts=nullptr,embed=nullptr;
    CUfunction startDrive=nullptr,drive=nullptr,releaseDrive=nullptr;
    CUfunction applyContactVelocities=nullptr;
    CUfunction discoverTriangleContacts=nullptr;
    CUfunction prepareContacts=nullptr;
    CUfunction refineTriangleContacts=nullptr;
    CUfunction discoverVolumeContacts=nullptr;
    CUdeviceptr volumeBodies=0;
    int volumeCount=0;
    CUdeviceptr tetBodies=0;
    CUdeviceptr remapOffsets=0,remapTets=0;
    std::vector<btGpuXpbdNode> restNodes;
    CUdeviceptr incident=0;
    CUdeviceptr targets=0,driveStart=0,inverseMasses=0;
    std::vector<int> targetIds;
    int targetCount=0;
    CUdeviceptr springs=0;
    int springCount=0;
    CUdeviceptr nodes=0,tets=0,rotations=0,forces=0,failure=0;
    CUdeviceptr acceptedNodes=0,acceptedRotations=0;
    CUdeviceptr disabledPairs=0;
    int disabledPairCount=0;
    CUdeviceptr friction=0;
    int frictionCount=0;
    CUdeviceptr accumulated=0;
    CUdeviceptr maps=0,positions=0,triangles=0,tree=0,leaves=0,boxes=0,ready=0,corrections=0;
    CUdeviceptr contactCache=0,contactCount=0,references=0;
    int contactCapacity=1048576;
    int count=0;
    int vertexCount=0,triangleCount=0,treeCount=0;
    float margin=0,clearance=0;
    std::vector<int> batches;
    std::string error;
    bool check(CUresult result)
    {
        if(result==CUDA_SUCCESS) return true;
        const char* message=nullptr;cuGetErrorString(result,&message);
        error=message?message:"CUDA error";return false;
    }
    void clear()
    {
        clearSurface();
        for(auto p:{nodes,tets,rotations,forces,failure,targets,driveStart,inverseMasses,copies,copyActive,solverTets,accumulated,acceptedNodes,acceptedRotations,springs}) if(p) cuMemFree(p);
        springs=0;springCount=0;
        acceptedNodes=acceptedRotations=0;
        accumulated=0;
        copies=copyActive=solverTets=0;copyLanes=0;
        targets=driveStart=inverseMasses=0;targetCount=0;targetIds.clear();
        nodes=tets=rotations=forces=failure=0;count=0;batches.clear();
        restNodes.clear();
    }
    void clearSurface()
    {
        clearGraph();
        if(friction) cuMemFree(friction);friction=0;frictionCount=0;
        if(disabledPairs) cuMemFree(disabledPairs);disabledPairs=0;disabledPairCount=0;
        if(incident) cuMemFree(incident);incident=0;
        if(tetBodies) cuMemFree(tetBodies);tetBodies=0;
        if(remapOffsets) cuMemFree(remapOffsets);remapOffsets=0;
        if(remapTets) cuMemFree(remapTets);remapTets=0;
        if(volumeBodies) cuMemFree(volumeBodies);volumeBodies=0;volumeCount=0;
        for(auto p:{maps,positions,triangles,tree,leaves,boxes,ready,corrections,contactCache,contactCount,references}) if(p) cuMemFree(p);
        contactCache=contactCount=references=0;
        maps=positions=triangles=tree=leaves=boxes=ready=corrections=0;vertexCount=triangleCount=treeCount=0;
    }
    void clearGraph()
    {
        if(graphExec) cuGraphExecDestroy(graphExec);
        if(graph) cuGraphDestroy(graph);
        graphExec=nullptr;graph=nullptr;
    }
    bool copy(CUdeviceptr& p,const void* source,size_t size)
    {
        return !size || (check(cuMemAlloc(&p,size)) && check(cuMemcpyHtoD(p,source,size)));
    }
    bool launch(CUfunction f,int n,void** args)
    {
        return !n || check(cuLaunchKernel(f,(n+127)/128,1,1,128,1,1,0,stream,args,nullptr));
    }
};
thread_local std::string creationError;
struct Current
{
    State& s;bool valid;
    explicit Current(State& state):s(state),valid(s.check(cuCtxPushCurrent(s.context))) {}
    ~Current() { if(valid) { CUcontext previous;cuCtxPopCurrent(&previous); } }
};
void destroy(void* value)
{
    if(!value) return;auto& s=*static_cast<State*>(value);
    if(s.context)
    {
        { Current current(s);if(current.valid) { s.clear();if(s.module) cuModuleUnload(s.module);if(s.stream) cuStreamDestroy(s.stream); } }
        cuDevicePrimaryCtxRelease(0);
    }
    delete &s;
}
void* create()
{
    auto* s=new State;
    if(s->check(cuInit(0)) && s->check(cuDevicePrimaryCtxRetain(&s->context,0)))
    {
        Current current(*s);
        const char* ptx=xpbdPtx;
        if(current.valid && s->check(cuStreamCreate(&s->stream,CU_STREAM_DEFAULT)) && s->check(cuModuleLoadData(&s->module,ptx)) &&
           s->check(cuModuleGetFunction(&s->predict,s->module,"predict")) &&
           s->check(cuModuleGetFunction(&s->solve,s->module,"solveTets")) &&
           s->check(cuModuleGetFunction(&s->initializeCopies,s->module,"initializeCopies")) &&
           s->check(cuModuleGetFunction(&s->averageCopies,s->module,"averageCopies")) &&
           s->check(cuModuleGetFunction(&s->settle,s->module,"settle")) &&
           s->check(cuModuleGetFunction(&s->skin,s->module,"skin")) &&
           s->check(cuModuleGetFunction(&s->refit,s->module,"refit")) &&
           s->check(cuModuleGetFunction(&s->discoverContacts,s->module,"discoverContacts")) &&
           s->check(cuModuleGetFunction(&s->discoverTriangleContacts,s->module,"discoverTriangleContacts")) &&
           s->check(cuModuleGetFunction(&s->prepareContacts,s->module,"prepareContacts")) &&
           s->check(cuModuleGetFunction(&s->refineTriangleContacts,s->module,"refineTriangleContacts")) &&
           s->check(cuModuleGetFunction(&s->discoverVolumeContacts,s->module,"discoverVolumeContacts")) &&
           s->check(cuModuleGetFunction(&s->countContacts,s->module,"countContacts")) &&
           s->check(cuModuleGetFunction(&s->solveContacts,s->module,"solveContacts")) &&
           s->check(cuModuleGetFunction(&s->applyContacts,s->module,"applyContacts")) &&
           s->check(cuModuleGetFunction(&s->applyContactVelocities,s->module,"applyContactVelocities")) &&
           s->check(cuModuleGetFunction(&s->embed,s->module,"embed")) &&
           s->check(cuModuleGetFunction(&s->startDrive,s->module,"startDrive")) &&
           s->check(cuModuleGetFunction(&s->drive,s->module,"drive")) &&
           s->check(cuModuleGetFunction(&s->releaseDrive,s->module,"releaseDrive"))) return s;
    }
    creationError=s->error;destroy(s);return nullptr;
}
bool finite(btGpuXpbdVec v) { return std::isfinite(v.x)&&std::isfinite(v.y)&&std::isfinite(v.z); }
int upload(void* value,int count,const btGpuXpbdNode* nodes,int tetCount,const btGpuXpbdTet* tets)
{
    auto& s=*static_cast<State*>(value);Current current(s);if(!current.valid) return 0;
    if(count<=0 || tetCount<0 || !nodes || (tetCount && !tets)) { s.error="Invalid mesh counts";return 0; }
    for(int i=0;i<count;++i)
        if(!finite(nodes[i].position)||!finite(nodes[i].velocity)||!std::isfinite(nodes[i].inverseMass)||nodes[i].inverseMass<0||
           !std::isfinite(nodes[i].damping)||nodes[i].damping<0) { s.error="Invalid node";return 0; }
    std::vector<std::vector<int>> nodeColors(count),colors;
    for(int i=0;i<tetCount;++i)
    {
        const auto& t=tets[i];
        for(int j=0;j<4;++j)
        {
            if(t.nodes[j]<0||t.nodes[j]>=count) { s.error="Invalid tetrahedron index";return 0; }
            for(int k=0;k<j;++k) if(t.nodes[j]==t.nodes[k]) { s.error="Repeated tetrahedron node";return 0; }
        }
        for(float v:t.inverseRest) if(!std::isfinite(v)) { s.error="Invalid rest matrix";return 0; }
        if(!std::isfinite(t.volume)||t.volume<=0||!std::isfinite(t.young)||t.young<=0||
           !std::isfinite(t.poisson)||t.poisson<0||t.poisson>=.5f||!std::isfinite(t.damping)||t.damping<0)
        { s.error="Invalid tetrahedron material";return 0; }
        std::vector<bool> used(colors.size(),false);
        for(int node:t.nodes) for(int c:nodeColors[node]) used[c]=true;
        int c=0;while(c<int(used.size())&&used[c]) ++c;
        if(c==int(colors.size())) colors.emplace_back();
        colors[c].push_back(i);for(int node:t.nodes) nodeColors[node].push_back(c);
    }
    s.clear();s.batches.push_back(0);
    std::vector<btGpuXpbdTet> ordered;ordered.reserve(tetCount);
    s.copyLanes=int((colors.size()+7)/8);
    std::vector<int> active(size_t(count)*s.copyLanes,0);
    std::vector<btGpuXpbdTet> work;
    const int passes=8;
    for(int pass=0;pass<passes;++pass)
    {
        for(int color=pass;color<int(colors.size());color+=passes)
        {
            int lane=color/8;
            for(int i:colors[color])
            {
                ordered.push_back(tets[i]);auto tet=tets[i];
                for(int& node:tet.nodes) { active[lane*count+node]=1;node+=lane*count; }
                work.push_back(tet);
            }
        }
        s.batches.push_back(int(ordered.size()));
    }
    struct Quat { float x,y,z,w; };
    std::vector<Quat> rotations(tetCount,{0,0,0,1});std::vector<btGpuXpbdVec> zero(count,{0,0,0});int failure=0;
    std::vector<float> masses(count);for(int i=0;i<count;++i) masses[i]=nodes[i].inverseMass;
    if(!s.copy(s.nodes,nodes,count*sizeof(*nodes))||!s.copy(s.tets,ordered.data(),ordered.size()*sizeof(*tets))||
       !s.copy(s.acceptedNodes,nodes,count*sizeof(*nodes))||!s.copy(s.acceptedRotations,rotations.data(),rotations.size()*sizeof(Quat))||
       !s.copy(s.rotations,rotations.data(),rotations.size()*sizeof(Quat))||!s.copy(s.forces,zero.data(),zero.size()*sizeof(zero[0]))||!s.copy(s.accumulated,zero.data(),zero.size()*sizeof(zero[0]))||
       !s.copy(s.failure,&failure,sizeof(failure))||!s.copy(s.inverseMasses,masses.data(),masses.size()*sizeof(float))||
       !s.copy(s.solverTets,work.data(),work.size()*sizeof(*tets))||!s.copy(s.copyActive,active.data(),active.size()*sizeof(int))||
       (s.copyLanes && !s.check(cuMemAlloc(&s.copies,size_t(count)*s.copyLanes*sizeof(*nodes))))) { s.clear();return 0; }
    s.restNodes.assign(nodes,nodes+count);s.count=count;return 1;
}
int forces(void* value,int count,const btGpuXpbdVec* values)
{
    auto& s=*static_cast<State*>(value);Current current(s);if(!current.valid) return 0;
    if(count!=s.count||!values) { s.error="Force count mismatch";return 0; }
    for(int i=0;i<count;++i) if(!finite(values[i])) { s.error="Nonfinite force";return 0; }
    return s.check(cuMemcpyHtoD(s.forces,values,count*sizeof(*values)));
}
int updateNodes(void* value,int count,const btGpuXpbdNode* nodes)
{
    auto& s=*static_cast<State*>(value);Current current(s);if(!current.valid) return 0;
    if(count!=s.count||count<=0||!nodes) { s.error="Node state count mismatch";return 0; }
    if(s.targetCount) { s.error="Release kinematic targets before replacing node state";return 0; }
    std::vector<float> masses(count);
    for(int i=0;i<count;++i)
    {
        const auto& n=nodes[i];
        if(!finite(n.position)||!finite(n.velocity)||!std::isfinite(n.inverseMass)||n.inverseMass<0||
           !std::isfinite(n.damping)||n.damping<0) { s.error="Invalid node state";return 0; }
        masses[i]=n.inverseMass;
    }
    return s.check(cuMemcpyHtoD(s.nodes,nodes,size_t(count)*sizeof(*nodes))) &&
        s.check(cuMemcpyHtoD(s.inverseMasses,masses.data(),masses.size()*sizeof(float)));
}
bool updateSurface(State& s)
{
    void* args[]={&s.nodes,&s.maps,&s.positions,&s.vertexCount};
    return s.launch(s.skin,s.vertexCount,args);
}
int disabledPairs(void* value,int count,const btGpuXpbdDisabledPair* pairs)
{
    auto& s=*static_cast<State*>(value);Current current(s);if(!current.valid) return 0;
    if(count<0 || (count && !pairs)) { s.error="Invalid disabled collision pairs";return 0; }
    std::vector<unsigned long long> keys;keys.reserve(count);
    for(int i=0;i<count;++i)
    {
        int a=pairs[i].bodyA,b=pairs[i].bodyB;
        if(a<0||b<0) { s.error="Invalid collision body ID";return 0; }
        if(a>b) std::swap(a,b);
        keys.push_back((static_cast<unsigned long long>(a)<<32)|static_cast<unsigned int>(b));
    }
    std::sort(keys.begin(),keys.end());keys.erase(std::unique(keys.begin(),keys.end()),keys.end());
    CUdeviceptr replacement=0;
    if(!s.copy(replacement,keys.data(),keys.size()*sizeof(keys[0])))
    {
        if(replacement) cuMemFree(replacement);
        return 0;
    }
    s.clearGraph();
    if(s.disabledPairs) cuMemFree(s.disabledPairs);
    s.disabledPairs=replacement;s.disabledPairCount=int(keys.size());return 1;
}
int friction(void* value,int count,const float* coefficients)
{
    auto& s=*static_cast<State*>(value);Current current(s);if(!current.valid) return 0;
    if(count<0||(count&&!coefficients)) { s.error="Invalid friction table";return 0; }
    for(int i=0;i<count;++i) if(!std::isfinite(coefficients[i])||coefficients[i]<0)
        { s.error="Friction must be finite and nonnegative";return 0; }
    CUdeviceptr replacement=0;
    if(!s.copy(replacement,coefficients,count*sizeof(float)))
    {
        if(replacement) cuMemFree(replacement);return 0;
    }
    s.clearGraph();
    if(s.friction) cuMemFree(s.friction);
    s.friction=replacement;s.frictionCount=count;return 1;
}
int surface(void* value,int vertexCount,const btGpuXpbdSurfaceVertex* vertices,int triangleCount,
    const btGpuXpbdSurfaceTriangle* triangles,float margin,float clearance)
{
    auto& s=*static_cast<State*>(value);Current current(s);if(!current.valid) return 0;
    if(vertexCount<=0||triangleCount<=0||!vertices||!triangles||!std::isfinite(margin)||!std::isfinite(clearance)||margin<=0||clearance<0||clearance>margin)
    { s.error="Invalid collision surface";return 0; }
    std::vector<btGpuXpbdNode> nodes(s.count);
    if(!s.check(cuMemcpyDtoH(nodes.data(),s.nodes,nodes.size()*sizeof(nodes[0])))) return 0;
    std::vector<btGpuXpbdVec> points(vertexCount),centers(triangleCount);
    for(int i=0;i<vertexCount;++i)
    {
        const auto& v=vertices[i];auto p=v.offset;float sum=0;bool dynamic=false;
        if(!finite(p)||v.body<0) { s.error="Invalid surface vertex";return 0; }
        for(int j=0;j<4;++j)
        {
            if(v.nodes[j]<-1||v.nodes[j]>=s.count||!std::isfinite(v.weights[j])) { s.error="Invalid surface support";return 0; }
            if(v.nodes[j]>=0)
            {
                dynamic=true;sum+=v.weights[j];auto n=nodes[v.nodes[j]].position;
                p.x+=n.x*v.weights[j];p.y+=n.y*v.weights[j];p.z+=n.z*v.weights[j];
            }
        }
        if(dynamic&&std::fabs(sum-1)>1.e-4f) { s.error="Surface weights must sum to one";return 0; }
        points[i]=p;
    }
    std::vector<int> incident(vertexCount,0);
    for(int i=0;i<triangleCount;++i)
    {
        auto& c=centers[i];c={0,0,0};int body=-1;
        for(int id:triangles[i].vertices)
        {
            if(id<0||id>=vertexCount) { s.error="Invalid collision triangle index";return 0; }
            incident[id]=1;
            if(body>=0&&vertices[id].body!=body) { s.error="Triangle spans bodies";return 0; }body=vertices[id].body;
            c.x+=points[id].x/3;c.y+=points[id].y/3;c.z+=points[id].z/3;
        }
    }
    struct Bvh { int left,right,parent,triangle,body; };
    struct Box { btGpuXpbdVec low,high; };
    std::vector<Bvh> tree;std::vector<int> leaves,order(triangleCount);std::iota(order.begin(),order.end(),0);
    std::function<int(int,int,int)> build=[&](int begin,int end,int parent)
    {
        int index=int(tree.size());tree.push_back({-1,-1,parent,-1,-1});
        if(end-begin==1)
        {
            tree[index].triangle=order[begin];tree[index].body=vertices[triangles[order[begin]].vertices[0]].body;
            leaves.push_back(index);return index;
        }
        float lo[3]={INFINITY,INFINITY,INFINITY},hi[3]={-INFINITY,-INFINITY,-INFINITY};
        for(int k=begin;k<end;++k) { auto p=centers[order[k]];float v[]={p.x,p.y,p.z};for(int j=0;j<3;++j) {lo[j]=std::min(lo[j],v[j]);hi[j]=std::max(hi[j],v[j]);} }
        int axis=0;for(int j=1;j<3;++j) if(hi[j]-lo[j]>hi[axis]-lo[axis]) axis=j;
        int mid=(begin+end)/2;
        std::nth_element(order.begin()+begin,order.begin()+mid,order.begin()+end,[&](int a,int b)
        { auto x=centers[a],y=centers[b];return axis==0?x.x<y.x:(axis==1?x.y<y.y:x.z<y.z); });
        int left=build(begin,mid,index),right=build(mid,end,index);tree[index].left=left;tree[index].right=right;
        tree[index].body=tree[left].body==tree[right].body?tree[left].body:-1;return index;
    };
    build(0,triangleCount,-1);
    int tetCount=s.batches.back();std::vector<btGpuXpbdTet> tets(tetCount);
    if(tetCount&&!s.check(cuMemcpyDtoH(tets.data(),s.tets,tets.size()*sizeof(tets[0])))) return 0;
    std::vector<int> parent(s.count),componentBody(s.count,-1),tetBodies(tetCount,-1);std::iota(parent.begin(),parent.end(),0);
    auto root=[&](int n) { while(parent[n]!=n) { parent[n]=parent[parent[n]];n=parent[n]; }return n; };
    for(const auto& t:tets) for(int j=1;j<4;++j) parent[root(t.nodes[j])]=root(t.nodes[0]);
    for(int i=0;i<vertexCount;++i) for(int n:vertices[i].nodes) if(n>=0)
    {
        int component=root(n);
        if(componentBody[component]>=0&&componentBody[component]!=vertices[i].body) { s.error="Simulation component spans collision bodies";return 0; }
        componentBody[component]=vertices[i].body;
    }
    for(int i=0;i<tetCount;++i) tetBodies[i]=componentBody[root(tets[i].nodes[0])];
    std::vector<int> volumeBodies=tetBodies;std::sort(volumeBodies.begin(),volumeBodies.end());
    volumeBodies.erase(std::unique(volumeBodies.begin(),volumeBodies.end()),volumeBodies.end());
    volumeBodies.erase(std::remove(volumeBodies.begin(),volumeBodies.end(),-1),volumeBodies.end());
    // Match PhysX exact-surface cooking: intersect rest AABBs, then include
    // the three vertex embedding tetrahedra, including extrapolated vertices.
    std::vector<Box> tetBounds(tetCount);
    std::map<std::array<int,4>,int> embeddingTets;
    auto extend=[](Box& box,btGpuXpbdVec p)
    {
        box.low={std::min(box.low.x,p.x),std::min(box.low.y,p.y),std::min(box.low.z,p.z)};
        box.high={std::max(box.high.x,p.x),std::max(box.high.y,p.y),std::max(box.high.z,p.z)};
    };
    for(int i=0;i<tetCount;++i)
    {
        auto p=s.restNodes[tets[i].nodes[0]].position;tetBounds[i]={p,p};
        std::array<int,4> key;
        for(int j=0;j<4;++j) { key[j]=tets[i].nodes[j];extend(tetBounds[i],s.restNodes[key[j]].position); }
        std::sort(key.begin(),key.end());embeddingTets[key]=i;
    }
    std::vector<int> offsets(triangleCount+1,0),candidates;
    for(int i=0;i<triangleCount;++i)
    {
        Box bound={{INFINITY,INFINITY,INFINITY},{-INFINITY,-INFINITY,-INFINITY}};
        bool dynamic=false;std::vector<int> embedded;
        for(int id:triangles[i].vertices)
        {
            const auto& v=vertices[id];auto p=v.offset;std::array<int,4> key;
            for(int j=0;j<4;++j)
            {
                key[j]=v.nodes[j];if(key[j]<0) continue;dynamic=true;
                auto x=s.restNodes[key[j]].position;
                p.x+=x.x*v.weights[j];p.y+=x.y*v.weights[j];p.z+=x.z*v.weights[j];
            }
            extend(bound,p);std::sort(key.begin(),key.end());auto entry=embeddingTets.find(key);
            if(entry!=embeddingTets.end()) embedded.push_back(entry->second);
        }
        if(dynamic)
        {
            int body=vertices[triangles[i].vertices[0]].body;
            for(int j=0;j<tetCount;++j)
            {
                const auto& b=tetBounds[j];
                if(tetBodies[j]==body&&bound.low.x<=b.high.x&&bound.high.x>=b.low.x&&
                    bound.low.y<=b.high.y&&bound.high.y>=b.low.y&&bound.low.z<=b.high.z&&bound.high.z>=b.low.z) candidates.push_back(j);
            }
            for(int j:embedded) if(std::find(candidates.begin()+offsets[i],candidates.end(),j)==candidates.end()) candidates.push_back(j);
        }
        offsets[i+1]=int(candidates.size());
    }
    s.clearSurface();std::vector<Box> boxes(tree.size());std::vector<int> ready(tree.size(),0);
    s.volumeCount=int(volumeBodies.size());
    std::vector<btGpuXpbdVec> corrections(s.count,{0,0,0});
    if(!s.copy(s.volumeBodies,volumeBodies.data(),volumeBodies.size()*sizeof(int))||
       !s.copy(s.remapOffsets,offsets.data(),offsets.size()*sizeof(int))||!s.copy(s.remapTets,candidates.data(),candidates.size()*sizeof(int))||
       !s.copy(s.tetBodies,tetBodies.data(),tetBodies.size()*sizeof(int))||!s.copy(s.maps,vertices,vertexCount*sizeof(*vertices))||!s.copy(s.positions,points.data(),points.size()*sizeof(points[0]))||
       !s.copy(s.triangles,triangles,triangleCount*sizeof(*triangles))||!s.copy(s.incident,incident.data(),incident.size()*sizeof(int))||!s.copy(s.tree,tree.data(),tree.size()*sizeof(Bvh))||
       !s.copy(s.leaves,leaves.data(),leaves.size()*sizeof(int))||!s.copy(s.boxes,boxes.data(),boxes.size()*sizeof(Box))||
       !s.copy(s.ready,ready.data(),ready.size()*sizeof(int))||!s.copy(s.corrections,corrections.data(),corrections.size()*sizeof(corrections[0])))
    { s.clearSurface();return 0; }
    if(!s.check(cuMemAlloc(&s.contactCache,s.contactCapacity*176))||!s.check(cuMemAlloc(&s.contactCount,sizeof(int)))||!s.check(cuMemAlloc(&s.references,s.count*sizeof(int))))
    { s.clearSurface();return 0; }
    s.vertexCount=vertexCount;s.triangleCount=triangleCount;s.treeCount=int(tree.size());s.margin=margin;s.clearance=clearance;return 1;
}
int surfacePositions(void* value,int count,btGpuXpbdVec* positions)
{
    auto& s=*static_cast<State*>(value);Current current(s);if(!current.valid) return 0;
    if(count!=s.vertexCount||!positions) { s.error="Surface download count mismatch";return 0; }
    return updateSurface(s)&&s.check(cuMemcpyDtoH(positions,s.positions,count*sizeof(*positions)));
}
int embed(void* value,int count,const btGpuXpbdVec* points,int begin,int end,int body,btGpuXpbdSurfaceVertex* maps)
{
    auto& s=*static_cast<State*>(value);Current current(s);if(!current.valid) return 0;
    if(count<=0||!points||!maps||begin<0||end>s.count||begin>=end||body<0||s.batches.empty()) { s.error="Invalid embedding request";return 0; }
    for(int i=0;i<count;++i) if(!finite(points[i])) { s.error="Nonfinite embedding point";return 0; }
    CUdeviceptr input=0,output=0;int failure=0,tetCount=s.batches.back();
    bool ok=s.copy(input,points,count*sizeof(*points))&&s.check(cuMemAlloc(&output,count*sizeof(*maps)))&&
        s.check(cuMemcpyHtoD(s.failure,&failure,sizeof(failure)));
    if(ok)
    {
        for(int start=0;start<count&&ok;start+=8192)
        {
            int batch=std::min(8192,count-start);CUdeviceptr in=input+start*sizeof(*points),out=output+start*sizeof(*maps);
            void* args[]={&s.nodes,&s.tets,&tetCount,&in,&out,&batch,&begin,&end,&body,&s.failure};
            ok=s.launch(s.embed,batch,args)&&s.check(cuCtxSynchronize());
        }
        ok=ok&&s.check(cuMemcpyDtoH(maps,output,count*sizeof(*maps)))&&s.check(cuMemcpyDtoH(&failure,s.failure,sizeof(failure)));
    }
    if(input) cuMemFree(input);if(output) cuMemFree(output);
    if(failure) { s.error="No eligible embedding tetrahedron";return 0; }return ok;
}
int targets(void* value,int count,const btGpuXpbdTarget* values)
{
    auto& s=*static_cast<State*>(value);Current current(s);if(!current.valid) return 0;
    if(count<0||count>s.count||(count&&!values)) { s.error="Invalid target count";return 0; }
    std::vector<int> ids;std::vector<bool> used(s.count,false);
    for(int i=0;i<count;++i)
    {
        if(values[i].node<0||values[i].node>=s.count||!finite(values[i].position)||used[values[i].node]) { s.error="Invalid kinematic target";return 0; }
        ids.push_back(values[i].node);used[values[i].node]=true;
    }
    if(ids!=s.targetIds)
    {
        s.clearGraph();
        void* args[]={&s.nodes,&s.targets,&s.inverseMasses,&s.targetCount};
        if(!s.launch(s.releaseDrive,s.targetCount,args)||!s.check(cuStreamSynchronize(s.stream))) return 0;
        if(s.targets) cuMemFree(s.targets);if(s.driveStart) cuMemFree(s.driveStart);
        s.targets=s.driveStart=0;s.targetCount=0;s.targetIds.clear();
        if(count&&(!s.copy(s.targets,values,count*sizeof(*values))||!s.check(cuMemAlloc(&s.driveStart,count*sizeof(btGpuXpbdVec))))) return 0;
        s.targetCount=count;s.targetIds=ids;return 1;
    }
    return !count||s.check(cuMemcpyHtoD(s.targets,values,count*sizeof(*values)));
}
int springs(void* value,int count,const btGpuXpbdSpring* values)
{
    auto& s=*static_cast<State*>(value);Current current(s);if(!current.valid) return 0;
    if(count<0||(count&&!values)) { s.error="Invalid mouse spring array";return 0; }
    for(int i=0;i<count;++i)
    {
        const auto& v=values[i];
        if(v.node<0||v.node>=s.count||!finite(v.target)||!std::isfinite(v.stiffness)||v.stiffness<0||
           !std::isfinite(v.damping)||v.damping<0||!std::isfinite(v.maxForce)||v.maxForce<0)
            { s.error="Invalid mouse spring";return 0; }
    }
    if(count!=s.springCount)
    {
        CUdeviceptr replacement=0;
        if(!s.copy(replacement,values,count*sizeof(*values)))
        {
            if(replacement) cuMemFree(replacement);return 0;
        }
        s.clearGraph();if(s.springs) cuMemFree(s.springs);
        s.springs=replacement;s.springCount=count;return 1;
    }
    return !count||s.check(cuMemcpyHtoD(s.springs,values,count*sizeof(*values)));
}
bool enqueueStep(State& s,const btGpuXpbdStep* config)
{
    float externalDt=config->dt/config->substeps,dt=externalDt/config->iterations;auto gravity=config->gravity;
    float speed=config->settlingSpeed,damping=config->settlingDamping;
    int settlingBegin=config->settlingBegin,settlingEnd=config->settlingEnd?config->settlingEnd:s.count;
    float fullDt=config->dt;
    void* startArgs[]={&s.nodes,&s.targets,&s.driveStart,&s.targetCount,&fullDt};
    if(!s.launch(s.startDrive,s.targetCount,startArgs)) return 0;
    for(int k=0;k<config->substeps*config->iterations;++k)
    {
        float forceDt=k%config->iterations==0?externalDt:0;
        if(k%config->iterations==0 && !s.check(cuMemsetD32Async(s.accumulated,0,s.count*3,s.stream))) return 0;
        // PhysX TGS prepares contacts before integrating the first increment.
        if(s.vertexCount&&k%config->iterations==0)
        {
            if(!updateSurface(s)||!s.check(cuMemsetD32Async(s.ready,0,s.treeCount,s.stream))||!s.check(cuMemsetD32Async(s.contactCount,0,1,s.stream))) return 0;
            void* refitArgs[]={&s.positions,&s.triangles,&s.tree,&s.leaves,&s.boxes,&s.ready,&s.triangleCount};
            if(!s.launch(s.refit,s.triangleCount,refitArgs)) return 0;
            void* contactArgs[]={&s.nodes,&s.maps,&s.positions,&s.triangles,&s.tree,&s.boxes,&s.contactCache,&s.contactCount,&s.contactCapacity,&s.vertexCount,&s.margin,&s.failure,&s.incident,&s.disabledPairs,&s.disabledPairCount};
            if(!s.launch(s.discoverContacts,s.vertexCount,contactArgs)) return 0;
            void* triangleArgs[]={&s.nodes,&s.maps,&s.positions,&s.triangles,&s.tree,&s.boxes,&s.contactCache,&s.contactCount,&s.contactCapacity,&s.triangleCount,&s.margin,&s.failure,&s.volumeBodies,&s.volumeCount,&s.disabledPairs,&s.disabledPairCount};
            if(!s.launch(s.discoverTriangleContacts,s.triangleCount,triangleArgs)) return 0;
            void* volumeArgs[]={&s.maps,&s.positions,&s.triangles,&s.tree,&s.boxes,&s.incident,&s.volumeBodies,&s.volumeCount,&s.vertexCount,&s.contactCache,&s.contactCount,&s.contactCapacity,&s.margin,&s.failure,&s.disabledPairs,&s.disabledPairCount};
            if(!s.launch(s.discoverVolumeContacts,s.vertexCount*s.volumeCount,volumeArgs)) return 0;
            void* refineArgs[]={&s.nodes,&s.maps,&s.positions,&s.triangles,&s.contactCache,&s.contactCount,&s.contactCapacity,&s.margin,&s.failure};
            if(!s.check(cuLaunchKernel(s.refineTriangleContacts,1024,1,1,32,4,1,0,s.stream,refineArgs,nullptr))) return 0;
            int tetCount=s.batches.back();
            void* prepArgs[]={&s.nodes,&s.tets,&s.tetBodies,&tetCount,&s.remapOffsets,&s.remapTets,&s.accumulated,&s.maps,&s.triangles,&s.contactCache,&s.contactCount,&s.contactCapacity,&s.friction,&s.frictionCount};
            if(!s.launch(s.prepareContacts,16384,prepArgs)) return 0;
        }
        float fraction=float(k+1)/float(config->substeps*config->iterations);
        float increment=1.f/float(config->substeps*config->iterations);
        void* driveArgs[]={&s.nodes,&s.targets,&s.driveStart,&s.accumulated,&s.targetCount,&fraction,&increment};
        if(!s.launch(s.drive,s.targetCount,driveArgs)) return 0;
        void* predictArgs[]={&s.nodes,&s.forces,&s.accumulated,&s.count,&dt,&gravity,&forceDt,&s.springs,&s.springCount};
        if(!s.launch(s.predict,s.count,predictArgs)) return 0;
        void* copyArgs[]={&s.nodes,&s.copies,&s.count,&s.copyLanes};
        if(!s.launch(s.initializeCopies,s.count*s.copyLanes,copyArgs)) return 0;
        for(size_t c=1;c<s.batches.size();++c)
        {
            int begin=s.batches[c-1],count=s.batches[c]-begin,first=k%config->iterations==0;
            void* args[]={&s.copies,&s.solverTets,&s.rotations,&begin,&count,&dt,&first};
            if(!s.launch(s.solve,count,args)) return 0;
        }
        void* averageArgs[]={&s.nodes,&s.copies,&s.copyActive,&s.accumulated,&s.count,&s.copyLanes,&dt};
        if(!s.launch(s.averageCopies,s.count,averageArgs)) return 0;
        if(s.vertexCount)
        {
            if(!s.check(cuMemsetD32Async(s.corrections,0,s.count*3,s.stream))||!s.check(cuMemsetD32Async(s.references,0,s.count,s.stream))) return 0;
            int velocityOnly=0;
            void* countArgs[]={&s.nodes,&s.accumulated,&s.maps,&s.triangles,&s.contactCache,&s.contactCount,&s.contactCapacity,&s.references,&s.clearance,&velocityOnly,&dt};
            if(!s.launch(s.countContacts,16384,countArgs)) return 0;
            void* solveArgs[]={&s.nodes,&s.accumulated,&s.maps,&s.triangles,&s.contactCache,&s.contactCount,&s.contactCapacity,&s.references,&s.corrections,&s.clearance,&velocityOnly,&dt};
            if(!s.launch(s.solveContacts,16384,solveArgs)) return 0;
            void* applyArgs[]={&s.nodes,&s.corrections,&s.accumulated,&s.count,&dt};
            if(!s.launch(s.applyContacts,s.count,applyArgs)) return 0;
        }
        if((k+1)%config->iterations==0)
        {
            for(int v=0;v<config->velocityIterations && s.vertexCount;++v)
            {
                if(!s.check(cuMemsetD32Async(s.corrections,0,s.count*3,s.stream))||!s.check(cuMemsetD32Async(s.references,0,s.count,s.stream))) return 0;
                int velocityOnly=1;
                void* countArgs[]={&s.nodes,&s.accumulated,&s.maps,&s.triangles,&s.contactCache,&s.contactCount,&s.contactCapacity,&s.references,&s.clearance,&velocityOnly,&dt};
                void* solveArgs[]={&s.nodes,&s.accumulated,&s.maps,&s.triangles,&s.contactCache,&s.contactCount,&s.contactCapacity,&s.references,&s.corrections,&s.clearance,&velocityOnly,&dt};
                void* applyArgs[]={&s.nodes,&s.corrections,&s.count};
                if(!s.launch(s.countContacts,16384,countArgs)||!s.launch(s.solveContacts,16384,solveArgs)||!s.launch(s.applyContactVelocities,s.count,applyArgs)) return 0;
            }
            void* settleArgs[]={&s.nodes,&s.accumulated,&s.count,&externalDt,&speed,&damping,&settlingBegin,&settlingEnd,&s.failure};
            if(!s.launch(s.settle,s.count,settleArgs)) return 0;
        }
    }
    return true;
}
int step(void* value,const btGpuXpbdStep* config)
{
    auto& s=*static_cast<State*>(value);Current current(s);if(!current.valid) return 0;
    if(!config||s.count<=0||!std::isfinite(config->dt)||config->dt<=0||config->substeps<=0||config->substeps>4096||config->iterations<=0||config->iterations>255||config->velocityIterations<0||config->velocityIterations>255||
       !finite(config->gravity)||!std::isfinite(config->settlingSpeed)||config->settlingSpeed<0||
       !std::isfinite(config->settlingDamping)||config->settlingDamping<0||config->settlingBegin<0||config->settlingEnd<0||config->settlingEnd>s.count) { s.error="Invalid step configuration";return 0; }
    int failure=0;
    if(!s.check(cuMemcpyHtoD(s.failure,&failure,sizeof(failure)))) return 0;
    if(!s.graphExec || std::memcmp(config,&s.graphConfig,sizeof(*config))!=0)
    {
        s.clearGraph();
        if(!s.check(cuStreamBeginCapture(s.stream,CU_STREAM_CAPTURE_MODE_THREAD_LOCAL))) return 0;
        bool captured=enqueueStep(s,config);
        CUresult end=cuStreamEndCapture(s.stream,&s.graph);
        if(!captured || !s.check(end)) { s.clearGraph();return 0; }
        if(!s.check(cuGraphInstantiateWithFlags(&s.graphExec,s.graph,0))) { s.clearGraph();return 0; }
        s.graphConfig=*config;
    }
    // Keep the transaction on device, including the warm-start rotations.
    const size_t nodeBytes=size_t(s.count)*sizeof(btGpuXpbdNode);
    const size_t rotationBytes=size_t(s.batches.back())*sizeof(float)*4;
    if(!s.check(cuMemcpyDtoDAsync(s.acceptedNodes,s.nodes,nodeBytes,s.stream)) ||
       (rotationBytes && !s.check(cuMemcpyDtoDAsync(s.acceptedRotations,s.rotations,rotationBytes,s.stream)))) return 0;
    if(!s.check(cuGraphLaunch(s.graphExec,s.stream)) || !s.check(cuStreamSynchronize(s.stream))) return 0;
    if(!s.check(cuMemcpyDtoH(&failure,s.failure,sizeof(failure)))) return 0;
    if(failure)
    {
        if(!s.check(cuMemcpyDtoDAsync(s.nodes,s.acceptedNodes,nodeBytes,s.stream)) ||
           (rotationBytes && !s.check(cuMemcpyDtoDAsync(s.rotations,s.acceptedRotations,rotationBytes,s.stream))) ||
           !s.check(cuStreamSynchronize(s.stream))) return 0;
        s.error=failure==4?"Contact capacity exceeded":(failure==2?"Collision traversal stack overflow":"Nonfinite GPU state");
        return 0;
    }
    return 1;
}
int download(void* value,int count,btGpuXpbdNode* nodes)
{
    auto& s=*static_cast<State*>(value);Current current(s);if(!current.valid) return 0;
    if(count!=s.count||!nodes) { s.error="Download count mismatch";return 0; }
    return s.check(cuMemcpyDtoH(nodes,s.nodes,count*sizeof(*nodes)));
}
const char* error(void* value) { return value?static_cast<State*>(value)->error.c_str():creationError.c_str(); }
const btGpuXpbdTgsApi api={create,destroy,upload,forces,step,download,error,surface,surfacePositions,embed,targets,disabledPairs,updateNodes,friction,springs};
}
extern "C" __declspec(dllexport) const btGpuXpbdTgsApi* btGetGpuXpbdTgsApi(int version)
{
    return version>=1&&version<=5?&api:nullptr;
}
