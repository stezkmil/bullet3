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
    CUfunction predict=nullptr,solve=nullptr,settle=nullptr,damp=nullptr;
    bool hasMaterialDamping=false;
    CUfunction predictRigids=nullptr,solveAttachments=nullptr,checkRigids=nullptr,bindRigids=nullptr,applyRigidContacts=nullptr;
    CUdeviceptr rigids=0,acceptedRigids=0,attachments=0,startRigids=0,rigidImpulse=0,rigidTorqueImpulse=0,rigidReferences=0;
    int rigidCount=0;
    std::vector<int> attachmentBatches;
    CUfunction solveJoints=nullptr;
    CUdeviceptr joints=0;
    std::vector<int> jointBatches;
    CUfunction dampFast=nullptr,solveFast=nullptr,initializeCopies=nullptr,averageCopies=nullptr;
    CUdeviceptr copies=0,copyActive=0,fastTets=0,rotations=0,acceptedRotations=0;
    int copyLanes=0;
    std::vector<int> fastBatches;
    CUfunction skin=nullptr,refit=nullptr,discoverContacts=nullptr,countContacts=nullptr,solveContacts=nullptr,applyContacts=nullptr,embed=nullptr;
    CUfunction startDrive=nullptr,drive=nullptr,releaseDrive=nullptr;
    CUfunction applyContactVelocities=nullptr;
    CUfunction discoverTriangleContacts=nullptr;
    CUfunction prepareContacts=nullptr;
    CUfunction refineTriangleContacts=nullptr,completeTriangleContacts=nullptr;
    CUfunction discoverVolumeContacts=nullptr;
    CUdeviceptr volumeBodies=0;
    int volumeCount=0;
    CUdeviceptr tetBodies=0;
    CUdeviceptr remapOffsets=0,remapTets=0;
    std::vector<btGpuXpbdNode> restNodes;
    CUdeviceptr incident=0,vertexTriangles=0;
    CUfunction reportContacts=nullptr,reportNotifications=nullptr,finishNotifications=nullptr;
    CUdeviceptr reports=0,reportCount=0,notificationPairs=0,notificationWinners=0;
    int notificationCapacity=0,lastReportCount=0;
    CUdeviceptr targets=0,driveStart=0,inverseMasses=0;
    std::vector<int> targetIds;
    int targetCount=0;
    CUdeviceptr springs=0;
    int springCount=0;
    CUdeviceptr nodes=0,tets=0,forces=0,failure=0;
    CUdeviceptr acceptedNodes=0;
    CUdeviceptr disabledPairs=0;
    int disabledPairCount=0;
    CUdeviceptr selfBodies=0;
    int selfBodyCount=0;
    CUdeviceptr friction=0;
    int frictionCount=0;
    CUdeviceptr accumulated=0;
    CUdeviceptr maps=0,positions=0,triangles=0,tree=0,leaves=0,boxes=0,ready=0,corrections=0;
    CUdeviceptr contactCache=0,contactCount=0,references=0;
    int contactCapacity=1048576;
    int count=0;
    int vertexCount=0,triangleCount=0,treeCount=0;
    double margin=0,clearance=0;
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
        if(joints) cuMemFree(joints);joints=0;jointBatches.clear();
        for(auto p:{rigids,acceptedRigids,attachments,startRigids,rigidImpulse,rigidTorqueImpulse,rigidReferences}) if(p) cuMemFree(p);
        rigids=acceptedRigids=attachments=startRigids=rigidImpulse=rigidTorqueImpulse=rigidReferences=0;rigidCount=0;attachmentBatches.clear();
        for(auto p:{copies,copyActive,fastTets,rotations,acceptedRotations}) if(p) cuMemFree(p);
        copies=copyActive=fastTets=rotations=acceptedRotations=0;copyLanes=0;fastBatches.clear();
        for(auto p:{nodes,tets,forces,failure,targets,driveStart,inverseMasses,accumulated,acceptedNodes,springs}) if(p) cuMemFree(p);
        springs=0;springCount=0;
        acceptedNodes=0;
        accumulated=0;
        targets=driveStart=inverseMasses=0;targetCount=0;targetIds.clear();
        nodes=tets=forces=failure=0;count=0;batches.clear();
        restNodes.clear();
    }
    void clearSurface()
    {
        clearGraph();
        if(friction) cuMemFree(friction);friction=0;frictionCount=0;
        if(disabledPairs) cuMemFree(disabledPairs);disabledPairs=0;disabledPairCount=0;
        if(selfBodies) cuMemFree(selfBodies);selfBodies=0;selfBodyCount=0;
        if(incident) cuMemFree(incident);incident=0;
        for(auto p:{vertexTriangles,reports,reportCount,notificationPairs,notificationWinners}) if(p) cuMemFree(p);
        vertexTriangles=reports=reportCount=notificationPairs=notificationWinners=0;
        notificationCapacity=lastReportCount=0;
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
           s->check(cuModuleGetFunction(&s->predictRigids,s->module,"predictRigids")) &&
           s->check(cuModuleGetFunction(&s->solveAttachments,s->module,"solveAttachments")) &&
           s->check(cuModuleGetFunction(&s->solveJoints,s->module,"solveJoints")) &&
           s->check(cuModuleGetFunction(&s->checkRigids,s->module,"checkRigids")) &&
           s->check(cuModuleGetFunction(&s->bindRigids,s->module,"bindRigids")) &&
           s->check(cuModuleGetFunction(&s->applyRigidContacts,s->module,"applyRigidContacts")) &&
           s->check(cuModuleGetFunction(&s->dampFast,s->module,"dampFastTets")) &&
           s->check(cuModuleGetFunction(&s->solveFast,s->module,"solveFastTets")) &&
           s->check(cuModuleGetFunction(&s->initializeCopies,s->module,"initializeCopies")) &&
           s->check(cuModuleGetFunction(&s->averageCopies,s->module,"averageCopies")) &&
           s->check(cuModuleGetFunction(&s->damp,s->module,"dampTets")) &&
           s->check(cuModuleGetFunction(&s->settle,s->module,"settle")) &&
           s->check(cuModuleGetFunction(&s->skin,s->module,"skin")) &&
           s->check(cuModuleGetFunction(&s->refit,s->module,"refit")) &&
           s->check(cuModuleGetFunction(&s->discoverContacts,s->module,"discoverContacts")) &&
           s->check(cuModuleGetFunction(&s->reportContacts,s->module,"reportContacts")) &&
           s->check(cuModuleGetFunction(&s->reportNotifications,s->module,"reportNotifications")) &&
           s->check(cuModuleGetFunction(&s->finishNotifications,s->module,"finishNotifications")) &&
           s->check(cuModuleGetFunction(&s->discoverTriangleContacts,s->module,"discoverTriangleContacts")) &&
           s->check(cuModuleGetFunction(&s->prepareContacts,s->module,"prepareContacts")) &&
           s->check(cuModuleGetFunction(&s->refineTriangleContacts,s->module,"refineTriangleContacts")) &&
           s->check(cuModuleGetFunction(&s->completeTriangleContacts,s->module,"completeTriangleContacts")) &&
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
        for(double v:t.inverseRest) if(!std::isfinite(v)) { s.error="Invalid rest matrix";return 0; }
        if(!std::isfinite(t.volume)||t.volume<=0||!std::isfinite(t.young)||t.young<=0||
           !std::isfinite(t.poisson)||t.poisson<0||t.poisson>=.5||!std::isfinite(t.damping)||t.damping<0)
        { s.error="Invalid tetrahedron material";return 0; }
        std::vector<bool> used(colors.size(),false);
        for(int node:t.nodes) for(int c:nodeColors[node]) used[c]=true;
        int c=0;while(c<int(used.size())&&used[c]) ++c;
        if(c==int(colors.size())) colors.emplace_back();
        colors[c].push_back(i);for(int node:t.nodes) nodeColors[node].push_back(c);
    }
    s.clear();s.batches.push_back(0);
    s.hasMaterialDamping=false;for(int i=0;i<tetCount;++i) s.hasMaterialDamping|=tets[i].damping>0;
    std::vector<btGpuXpbdTet> ordered;ordered.reserve(tetCount);
    for(const auto& color:colors)
    {
        for(int i:color) ordered.push_back(tets[i]);
        s.batches.push_back(int(ordered.size()));
    }
    std::vector<btGpuXpbdVec> zero(count,{0,0,0});int failure=0;
    std::vector<double> masses(count);for(int i=0;i<count;++i) masses[i]=nodes[i].inverseMass;
    if(!s.copy(s.nodes,nodes,count*sizeof(*nodes))||!s.copy(s.tets,ordered.data(),ordered.size()*sizeof(*tets))||
       !s.copy(s.acceptedNodes,nodes,count*sizeof(*nodes))||
       !s.copy(s.forces,zero.data(),zero.size()*sizeof(zero[0]))||!s.copy(s.accumulated,zero.data(),zero.size()*sizeof(zero[0]))||
       !s.copy(s.failure,&failure,sizeof(failure))||!s.copy(s.inverseMasses,masses.data(),masses.size()*sizeof(double))) { s.clear();return 0; }
    // Eight serial passes, with disjoint node copies across groups of eight colors.
    constexpr int passes=8;
    s.copyLanes=(int(colors.size())+passes-1)/passes;
    std::vector<int> active(size_t(s.copyLanes)*count,0);
    std::vector<btGpuXpbdTet> work;
    s.fastBatches.push_back(0);
    for(int pass=0;pass<passes;++pass)
    {
        for(int color=pass;color<int(colors.size());color+=passes)
            for(int i:colors[color])
            {
                auto tet=tets[i];int lane=color/passes;
                for(int& node:tet.nodes) { active[size_t(lane)*count+node]=1;node+=lane*count; }
                work.push_back(tet);
            }
        s.fastBatches.push_back(int(work.size()));
    }
    struct Quat { float x,y,z,w; };
    std::vector<Quat> rotations(tetCount,{0,0,0,1});
    if(!s.copy(s.fastTets,work.data(),work.size()*sizeof(work[0]))||
       !s.copy(s.copyActive,active.data(),active.size()*sizeof(int))||
       !s.copy(s.rotations,rotations.data(),rotations.size()*sizeof(Quat))||
       !s.copy(s.acceptedRotations,rotations.data(),rotations.size()*sizeof(Quat))||
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
    std::vector<double> masses(count);
    for(int i=0;i<count;++i)
    {
        const auto& n=nodes[i];
        if(!finite(n.position)||!finite(n.velocity)||!std::isfinite(n.inverseMass)||n.inverseMass<0||
           !std::isfinite(n.damping)||n.damping<0) { s.error="Invalid node state";return 0; }
        masses[i]=n.inverseMass;
    }
    return s.check(cuMemcpyHtoD(s.nodes,nodes,size_t(count)*sizeof(*nodes))) &&
        s.check(cuMemcpyHtoD(s.inverseMasses,masses.data(),masses.size()*sizeof(double)));
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
int friction(void* value,int count,const double* coefficients)
{
    auto& s=*static_cast<State*>(value);Current current(s);if(!current.valid) return 0;
    if(count<0||(count&&!coefficients)) { s.error="Invalid friction table";return 0; }
    for(int i=0;i<count;++i) if(!std::isfinite(coefficients[i])||coefficients[i]<0)
        { s.error="Friction must be finite and nonnegative";return 0; }
    CUdeviceptr replacement=0;
    if(!s.copy(replacement,coefficients,count*sizeof(double)))
    {
        if(replacement) cuMemFree(replacement);return 0;
    }
    s.clearGraph();
    if(s.friction) cuMemFree(s.friction);
    s.friction=replacement;s.frictionCount=count;return 1;
}
int surface(void* value,int vertexCount,const btGpuXpbdSurfaceVertex* vertices,int triangleCount,
    const btGpuXpbdSurfaceTriangle* triangles,double margin,double clearance)
{
    auto& s=*static_cast<State*>(value);Current current(s);if(!current.valid) return 0;
    if(vertexCount<=0||triangleCount<=0||!vertices||!triangles||!std::isfinite(margin)||!std::isfinite(clearance)||margin<=0||clearance<0||clearance>margin)
    { s.error="Invalid collision surface";return 0; }
    std::vector<btGpuXpbdNode> nodes(s.count);
    if(!s.check(cuMemcpyDtoH(nodes.data(),s.nodes,nodes.size()*sizeof(nodes[0])))) return 0;
    std::vector<btGpuXpbdRigid> rigidStates(s.rigidCount);
    if(s.rigidCount&&!s.check(cuMemcpyDtoH(rigidStates.data(),s.rigids,rigidStates.size()*sizeof(rigidStates[0])))) return 0;
    std::vector<btGpuXpbdVec> points(vertexCount),centers(triangleCount);
    for(int i=0;i<vertexCount;++i)
    {
        const auto& v=vertices[i];auto p=v.offset;double sum=0;bool dynamic=false;
        if(!finite(p)||v.body<0) { s.error="Invalid surface vertex";return 0; }
        if(v.nodes[0]<-1)
            for(int j=0;j<4;++j)
                if(v.weights[j]!=0||(j&&v.nodes[j]!=-1)) {s.error="Rigid surface support must contain only a rigid-local point";return 0;}
        for(int j=0;j<4;++j)
        {
            if((v.nodes[j]<-1&&(j!=0||-2-v.nodes[j]>=s.rigidCount))||v.nodes[j]>=s.count||!std::isfinite(v.weights[j])) { s.error="Invalid surface support";return 0; }
            if(v.nodes[j]>=0)
            {
                dynamic=true;sum+=v.weights[j];auto n=nodes[v.nodes[j]].position;
                p.x+=n.x*v.weights[j];p.y+=n.y*v.weights[j];p.z+=n.z*v.weights[j];
            }
        }
        if(dynamic&&std::fabs(sum-1)>1.e-4) { s.error="Surface weights must sum to one";return 0; }
        if(v.nodes[0]<-1)
        {
            const auto& r=rigidStates[-2-v.nodes[0]];double x=r.orientation[0],y=r.orientation[1],z=r.orientation[2],w=r.orientation[3];
            double tx=2*(y*p.z-z*p.y),ty=2*(z*p.x-x*p.z),tz=2*(x*p.y-y*p.x);
            p={p.x+w*tx+y*tz-z*ty+r.position.x,p.y+w*ty+z*tx-x*tz+r.position.y,p.z+w*tz+x*ty-y*tx+r.position.z};
        }
        points[i]=p;
    }
    std::vector<int> incident(vertexCount,0),vertexTriangles(vertexCount,-1);
    for(int i=0;i<triangleCount;++i)
    {
        auto& c=centers[i];c={0,0,0};int body=-1,rigid=-1;
        for(int id:triangles[i].vertices)
        {
            if(id<0||id>=vertexCount) { s.error="Invalid collision triangle index";return 0; }
            incident[id]=1;if(vertexTriangles[id]<0) vertexTriangles[id]=i;
            if(body>=0&&vertices[id].body!=body) { s.error="Triangle spans bodies";return 0; }body=vertices[id].body;
            const int support=vertices[id].nodes[0]<-1?-2-vertices[id].nodes[0]:-1;
            if(id==triangles[i].vertices[0]) rigid=support;
            else if(support!=rigid) {s.error="Triangle spans rigid supports";return 0;}
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
        double lo[3]={INFINITY,INFINITY,INFINITY},hi[3]={-INFINITY,-INFINITY,-INFINITY};
        for(int k=begin;k<end;++k) { auto p=centers[order[k]];double v[]={p.x,p.y,p.z};for(int j=0;j<3;++j) {lo[j]=std::min(lo[j],v[j]);hi[j]=std::max(hi[j],v[j]);} }
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
    std::vector<bool> volumeNode(s.count,false);
    for(const auto& tet:tets) for(int n:tet.nodes) volumeNode[n]=true;
    for(int i=0;i<vertexCount;++i) if(!incident[i]&&vertices[i].nodes[0]>=0&&
        vertices[i].nodes[1]<0&&vertices[i].nodes[2]<0&&vertices[i].nodes[3]<0&&vertices[i].weights[0]==1&&
        vertices[i].offset.x==0&&vertices[i].offset.y==0&&vertices[i].offset.z==0&&volumeNode[vertices[i].nodes[0]])
        incident[i]=2; // Material samples use the same face-only contact rule as mesh vertices.
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
       !s.copy(s.vertexTriangles,vertexTriangles.data(),vertexTriangles.size()*sizeof(int))||!s.copy(s.triangles,triangles,triangleCount*sizeof(*triangles))||!s.copy(s.incident,incident.data(),incident.size()*sizeof(int))||!s.copy(s.tree,tree.data(),tree.size()*sizeof(Bvh))||
       !s.copy(s.leaves,leaves.data(),leaves.size()*sizeof(int))||!s.copy(s.boxes,boxes.data(),boxes.size()*sizeof(Box))||
       !s.copy(s.ready,ready.data(),ready.size()*sizeof(int))||!s.copy(s.corrections,corrections.data(),corrections.size()*sizeof(corrections[0])))
    { s.clearSurface();return 0; }
    if(!s.check(cuMemAlloc(&s.contactCache,s.contactCapacity*312))||!s.check(cuMemAlloc(&s.contactCount,sizeof(int)))||!s.check(cuMemAlloc(&s.references,s.count*sizeof(int))))
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
int setRigids(void* value,int count,const btGpuXpbdRigid* bodies,bool update)
{
    auto& s=*static_cast<State*>(value);Current current(s);if(!current.valid) return 0;
    if(count<0||(count&&!bodies)) { s.error="Invalid rigid count";return 0; }
    for(int i=0;i<count;++i)
    {
        const auto& r=bodies[i];double q2=0;for(double v:r.orientation) q2+=v*v;
        if(!finite(r.position)||!finite(r.linearVelocity)||!finite(r.angularVelocity)||!finite(r.force)||!finite(r.torque)||
           !finite(r.inverseInertia)||!std::isfinite(r.inverseMass)||r.inverseMass<0||
           !std::isfinite(q2)||std::fabs(q2-1)>1.e-8||!std::isfinite(r.linearDamping)||r.linearDamping<0||
           !std::isfinite(r.angularDamping)||r.angularDamping<0||r.inverseInertia.x<0||r.inverseInertia.y<0||r.inverseInertia.z<0||
           (r.inverseMass==0&&(r.inverseInertia.x!=0||r.inverseInertia.y!=0||r.inverseInertia.z!=0))||
           (r.inverseMass>0&&(r.inverseInertia.x<=0||r.inverseInertia.y<=0||r.inverseInertia.z<=0)))
        { s.error="Invalid rigid state or inertia";return 0; }
    }
    if(update)
    {
        if(count!=s.rigidCount) {s.error="Rigid update count mismatch";return 0;}
        return !count||s.check(cuMemcpyHtoD(s.rigids,bodies,size_t(count)*sizeof(*bodies)));
    }
    CUdeviceptr next=0,accepted=0,start=0,impulse=0,torque=0,references=0;
    if(!s.copy(next,bodies,size_t(count)*sizeof(*bodies))||!s.copy(accepted,bodies,size_t(count)*sizeof(*bodies))||
       (count&&(!s.check(cuMemAlloc(&start,size_t(count)*sizeof(*bodies)))||
       !s.check(cuMemAlloc(&impulse,size_t(count)*sizeof(btGpuXpbdVec)))||
       !s.check(cuMemAlloc(&torque,size_t(count)*sizeof(btGpuXpbdVec)))||
       !s.check(cuMemAlloc(&references,size_t(count)*sizeof(int))))))
    {for(auto p:{next,accepted,start,impulse,torque,references}) if(p) cuMemFree(p);return 0;}
    // Replacing body IDs invalidates every surface and attachment referring to them.
    s.clearSurface();for(auto p:{s.rigids,s.acceptedRigids,s.attachments,s.startRigids,s.rigidImpulse,s.rigidTorqueImpulse,s.rigidReferences}) if(p) cuMemFree(p);
    if(s.joints) cuMemFree(s.joints);s.joints=0;s.jointBatches.clear();
    s.startRigids=start;s.rigidImpulse=impulse;s.rigidTorqueImpulse=torque;s.rigidReferences=references;
    s.rigids=next;s.acceptedRigids=accepted;s.rigidCount=count;s.attachments=0;s.attachmentBatches.clear();
    void* args[]={&s.rigids,&s.startRigids,&s.rigidImpulse,&s.rigidTorqueImpulse,&s.rigidReferences};
    return s.launch(s.bindRigids,1,args)&&s.check(cuStreamSynchronize(s.stream));
}
int rigids(void* value,int count,const btGpuXpbdRigid* bodies) { return setRigids(value,count,bodies,false); }
int updateRigids(void* value,int count,const btGpuXpbdRigid* bodies) { return setRigids(value,count,bodies,true); }
int downloadRigids(void* value,int count,btGpuXpbdRigid* bodies)
{
    auto& s=*static_cast<State*>(value);Current current(s);if(!current.valid) return 0;
    if(count!=s.rigidCount||(count&&!bodies)) { s.error="Rigid download count mismatch";return 0; }
    return !count||s.check(cuMemcpyDtoH(bodies,s.rigids,size_t(count)*sizeof(*bodies)));
}
int attachments(void* value,int count,const btGpuXpbdAttachment* input)
{
    auto& s=*static_cast<State*>(value);Current current(s);if(!current.valid) return 0;
    if(count<0||(count&&!input)) { s.error="Invalid attachment count";return 0; }
    std::vector<std::vector<int>> colors,entityColors(size_t(s.count)+s.rigidCount);
    for(int i=0;i<count;++i)
    {
        const auto& a=input[i];std::vector<int> entities;
        if(!std::isfinite(a.compliance)||a.compliance<0) {s.error="Invalid attachment compliance";return 0;}
        for(const auto& e:{a.a,a.b})
        {
            if(e.rigid< -1||e.rigid>=s.rigidCount||!finite(e.point)) {s.error="Invalid attachment endpoint";return 0;}
            if(e.rigid>=0) entities.push_back(s.count+e.rigid);
            else
            {
                double sum=0;bool supported=false;
                for(int j=0;j<4;++j)
                {
                    if(e.nodes[j]< -1||e.nodes[j]>=s.count||!std::isfinite(e.weights[j])||e.weights[j]<0||
                       (e.nodes[j]<0&&e.weights[j]!=0)) {s.error="Invalid attachment weights";return 0;}
                    if(e.nodes[j]>=0&&e.weights[j]>0) {entities.push_back(e.nodes[j]);sum+=e.weights[j];supported=true;}
                }
                if(supported&&std::fabs(sum-1)>1.e-8) {s.error="Attachment weights must sum to one";return 0;}
            }
        }
        std::sort(entities.begin(),entities.end());
        if(std::adjacent_find(entities.begin(),entities.end())!=entities.end()) {s.error="Attachment endpoints share a degree of freedom";return 0;}
        std::vector<bool> used(colors.size(),false);for(int entity:entities) for(int c:entityColors[entity]) used[c]=true;
        int c=0;while(c<int(used.size())&&used[c]) ++c;if(c==int(colors.size())) colors.emplace_back();
        colors[c].push_back(i);for(int entity:entities) entityColors[entity].push_back(c);
    }
    std::vector<btGpuXpbdAttachment> ordered;std::vector<int> batches{0};
    for(const auto& color:colors) {for(int i:color) ordered.push_back(input[i]);batches.push_back(int(ordered.size()));}
    CUdeviceptr next=0;if(!s.copy(next,ordered.data(),ordered.size()*sizeof(ordered[0]))) {if(next) cuMemFree(next);return 0;}
    s.clearGraph();if(s.attachments) cuMemFree(s.attachments);s.attachments=next;s.attachmentBatches=std::move(batches);return 1;
}
int selfCollision(void* value,int count,const int* bodies)
{
    auto& s=*static_cast<State*>(value);Current current(s);if(!current.valid) return 0;
    if(count<0||(count&&!bodies)) {s.error="Invalid self-collision body count";return 0;}
    std::vector<int> sorted;if(count) sorted.assign(bodies,bodies+count);
    for(int b:sorted) if(b<0) {s.error="Invalid self-collision body ID";return 0;}
    std::sort(sorted.begin(),sorted.end());sorted.erase(std::unique(sorted.begin(),sorted.end()),sorted.end());
    CUdeviceptr next=0;if(!s.copy(next,sorted.data(),sorted.size()*sizeof(int))) {if(next) cuMemFree(next);return 0;}
    s.clearGraph();if(s.selfBodies) cuMemFree(s.selfBodies);s.selfBodies=next;s.selfBodyCount=int(sorted.size());return 1;
}
int joints(void* value,int count,const btGpuXpbdJoint* input)
{
    auto& s=*static_cast<State*>(value);Current current(s);if(!current.valid) return 0;
    if(count<0||(count&&!input)) {s.error="Invalid joint count";return 0;}
    std::vector<std::vector<int>> colors,entityColors(s.rigidCount);
    for(int i=0;i<count;++i)
    {
        const auto& j=input[i];
        if(j.rotationOrder<0||j.rotationOrder>5||(j.bodies[0]>=0&&j.bodies[0]==j.bodies[1])) {s.error="Invalid joint bodies or rotation order";return 0;}
        for(int side=0;side<2;++side)
        {
            double norm=0;for(double q:j.orientations[side]) norm+=q*q;
            if(j.bodies[side]<-1||j.bodies[side]>=s.rigidCount||!finite(j.points[side])||!std::isfinite(norm)||std::fabs(norm-1)>1.e-8)
            {s.error="Invalid joint frame";return 0;}
        }
        for(const auto& a:j.axes)
        {
            for(double v:{a.lower,a.upper,a.stiffness,a.damping,a.equilibrium,a.velocity,a.maxForce,a.servoTarget,a.bounce,a.stopERP,a.stopCFM,a.motorERP,a.motorCFM})
                if(!std::isfinite(v)) {s.error="Nonfinite joint parameter";return 0;}
            if(a.stiffness<0||a.damping<0||a.maxForce<0||a.bounce<0||a.bounce>1||a.stopERP<0||a.stopERP>1||a.motorERP<0||a.motorERP>1||a.stopCFM<0||a.motorCFM<0||(a.flags&~31))
            {s.error="Invalid joint parameter range";return 0;}
        }
        std::vector<bool> used(colors.size(),false);
        for(int body:j.bodies) if(body>=0) for(int c:entityColors[body]) used[c]=true;
        int color=0;while(color<int(used.size())&&used[color]) ++color;if(color==int(colors.size())) colors.emplace_back();
        colors[color].push_back(i);for(int body:j.bodies) if(body>=0) entityColors[body].push_back(color);
    }
    std::vector<btGpuXpbdJoint> ordered;std::vector<int> batches{0};
    for(const auto& color:colors){for(int i:color) ordered.push_back(input[i]);batches.push_back(int(ordered.size()));}
    CUdeviceptr next=0;if(!s.copy(next,ordered.data(),ordered.size()*sizeof(ordered[0]))) {if(next) cuMemFree(next);return 0;}
    s.clearGraph();if(s.joints) cuMemFree(s.joints);s.joints=next;s.jointBatches=std::move(batches);return 1;
}
bool resetRigidCorrections(State& s)
{
    return !s.rigidCount || (s.check(cuMemsetD32Async(s.rigidImpulse,0,s.rigidCount*sizeof(btGpuXpbdVec)/sizeof(unsigned),s.stream))&&
        s.check(cuMemsetD32Async(s.rigidTorqueImpulse,0,s.rigidCount*sizeof(btGpuXpbdVec)/sizeof(unsigned),s.stream))&&
        s.check(cuMemsetD32Async(s.rigidReferences,0,s.rigidCount,s.stream)));
}
bool discoverSurfaceContacts(State& s)
{
    if(!s.vertexCount) return true;
            if(!updateSurface(s)||!s.check(cuMemsetD32Async(s.ready,0,s.treeCount,s.stream))||!s.check(cuMemsetD32Async(s.contactCount,0,1,s.stream))) return 0;
            void* refitArgs[]={&s.positions,&s.triangles,&s.tree,&s.leaves,&s.boxes,&s.ready,&s.triangleCount};
            if(!s.launch(s.refit,s.triangleCount,refitArgs)) return 0;
            void* contactArgs[]={&s.nodes,&s.maps,&s.positions,&s.triangles,&s.tree,&s.boxes,&s.contactCache,&s.contactCount,&s.contactCapacity,&s.vertexCount,&s.margin,&s.failure,&s.incident,&s.disabledPairs,&s.disabledPairCount,&s.selfBodies,&s.selfBodyCount};
            if(!s.launch(s.discoverContacts,s.vertexCount,contactArgs)) return 0;
            void* triangleArgs[]={&s.nodes,&s.maps,&s.positions,&s.triangles,&s.tree,&s.boxes,&s.contactCache,&s.contactCount,&s.contactCapacity,&s.triangleCount,&s.margin,&s.failure,&s.disabledPairs,&s.disabledPairCount,&s.selfBodies,&s.selfBodyCount};
            if(!s.launch(s.discoverTriangleContacts,s.triangleCount,triangleArgs)) return 0;
            void* volumeArgs[]={&s.maps,&s.positions,&s.triangles,&s.tree,&s.boxes,&s.incident,&s.volumeBodies,&s.volumeCount,&s.vertexCount,&s.contactCache,&s.contactCount,&s.contactCapacity,&s.margin,&s.failure,&s.disabledPairs,&s.disabledPairCount,&s.selfBodies,&s.selfBodyCount};
            if(!s.launch(s.discoverVolumeContacts,s.vertexCount*s.volumeCount,volumeArgs)) return 0;
            void* refineArgs[]={&s.nodes,&s.maps,&s.positions,&s.triangles,&s.contactCache,&s.contactCount,&s.contactCapacity,&s.margin,&s.failure};
            if(!s.check(cuLaunchKernel(s.refineTriangleContacts,1024,1,1,32,4,1,0,s.stream,refineArgs,nullptr))) return 0;
            void* completeArgs[]={&s.positions,&s.triangles,&s.contactCache,&s.contactCount,&s.contactCapacity,&s.margin,&s.failure};
            if(!s.launch(s.completeTriangleContacts,16384,completeArgs)) return 0;
            int tetCount=s.batches.back();
            void* prepArgs[]={&s.nodes,&s.tets,&s.tetBodies,&tetCount,&s.remapOffsets,&s.remapTets,&s.accumulated,&s.maps,&s.triangles,&s.contactCache,&s.contactCount,&s.contactCapacity,&s.friction,&s.frictionCount};
            if(!s.launch(s.prepareContacts,16384,prepArgs)) return 0;
    return true;
}
bool enqueueStep(State& s,const btGpuXpbdStep* config)
{
    double externalDt=config->dt/config->substeps,dt=externalDt/config->iterations;auto gravity=config->gravity;
    double speed=config->settlingSpeed,damping=config->settlingDamping;
    int settlingBegin=config->settlingBegin,settlingEnd=config->settlingEnd?config->settlingEnd:s.count;
    double fullDt=config->dt;
    void* startArgs[]={&s.nodes,&s.targets,&s.driveStart,&s.targetCount,&fullDt};
    if(!s.launch(s.startDrive,s.targetCount,startArgs)) return 0;
    for(int k=0;k<config->substeps*config->iterations;++k)
    {
        double forceDt=dt;
        if(k%config->iterations==0 && !s.check(cuMemsetD32Async(s.accumulated,0,s.count*sizeof(btGpuXpbdVec)/sizeof(unsigned),s.stream))) return 0;
        if(s.rigidCount&&k%config->iterations==0&&!s.check(cuMemcpyDtoDAsync(s.startRigids,s.rigids,size_t(s.rigidCount)*sizeof(btGpuXpbdRigid),s.stream))) return 0;
        // PhysX TGS prepares contacts before integrating the first increment.
        if(s.vertexCount&&k%config->iterations==0)
        {
            if(!discoverSurfaceContacts(s)) return false;
        }
        double fraction=double(k+1)/double(config->substeps*config->iterations);
        double increment=1./double(config->substeps*config->iterations);
        void* driveArgs[]={&s.nodes,&s.targets,&s.driveStart,&s.accumulated,&s.targetCount,&fraction,&increment};
        if(!s.launch(s.drive,s.targetCount,driveArgs)) return 0;
        const bool fast=config->materialMode==1;
        void* copyArgs[]={&s.nodes,&s.copies,&s.copyActive,&s.count,&s.copyLanes};
        int positions=0;
        void* averageArgs[]={&s.nodes,&s.copies,&s.copyActive,&s.accumulated,&s.count,&s.copyLanes,&dt,&positions};
        auto& batches=fast?s.fastBatches:s.batches;
        CUdeviceptr materialNodes=fast?s.copies:s.nodes,materialTets=fast?s.fastTets:s.tets;
        if(s.hasMaterialDamping)
        {
            if(fast&&!s.launch(s.initializeCopies,s.count*s.copyLanes,copyArgs)) return 0;
            for(size_t c=1;c<batches.size();++c)
            {
                int begin=batches[c-1],count=batches[c]-begin;
                void* args[]={&materialNodes,&materialTets,&begin,&count,&dt};
                if(!s.launch(fast?s.dampFast:s.damp,count,args)) return 0;
            }
            if(fast&&!s.launch(s.averageCopies,s.count,averageArgs)) return 0;
        }
        void* predictArgs[]={&s.nodes,&s.forces,&s.accumulated,&s.count,&dt,&gravity,&forceDt,&s.springs,&s.springCount};
        if(!s.launch(s.predict,s.count,predictArgs)) return 0;
        void* rigidArgs[]={&s.rigids,&s.rigidCount,&dt,&gravity};
        if(!s.launch(s.predictRigids,s.rigidCount,rigidArgs)) return 0;
        if(fast&&!s.launch(s.initializeCopies,s.count*s.copyLanes,copyArgs)) return 0;
        for(size_t c=1;c<batches.size();++c)
        {
            int begin=batches[c-1],count=batches[c]-begin,first=k%config->iterations==0;
            void* args[]={&s.nodes,&s.tets,&s.accumulated,&begin,&count,&dt};
            void* fastArgs[]={&s.copies,&s.fastTets,&s.rotations,&begin,&count,&dt,&first};
            if(!s.launch(fast?s.solveFast:s.solve,count,fast?fastArgs:args)) return 0;
        }
        positions=1;
        if(fast&&!s.launch(s.averageCopies,s.count,averageArgs)) return 0;
        if(s.vertexCount)
        {
            if(!s.check(cuMemsetD32Async(s.corrections,0,s.count*sizeof(btGpuXpbdVec)/sizeof(unsigned),s.stream))||!s.check(cuMemsetD32Async(s.references,0,s.count,s.stream))) return 0;
            if(!resetRigidCorrections(s)) return 0;
            int velocityOnly=0;
            void* countArgs[]={&s.nodes,&s.accumulated,&s.maps,&s.triangles,&s.contactCache,&s.contactCount,&s.contactCapacity,&s.references,&s.clearance,&velocityOnly,&dt};
            if(!s.launch(s.countContacts,16384,countArgs)) return 0;
            void* solveArgs[]={&s.nodes,&s.accumulated,&s.maps,&s.triangles,&s.contactCache,&s.contactCount,&s.contactCapacity,&s.references,&s.corrections,&s.clearance,&velocityOnly,&dt};
            if(!s.launch(s.solveContacts,16384,solveArgs)) return 0;
            void* applyArgs[]={&s.nodes,&s.corrections,&s.accumulated,&s.count,&dt};
            if(!s.launch(s.applyContacts,s.count,applyArgs)) return 0;
            void* rigidApply[]={&s.rigids,&s.rigidCount,&dt,&velocityOnly};
            if(!s.launch(s.applyRigidContacts,s.rigidCount,rigidApply)) return 0;
        }
        for(size_t c=1;c<s.attachmentBatches.size();++c)
        {
            int begin=s.attachmentBatches[c-1],count=s.attachmentBatches[c]-begin;
            void* args[]={&s.nodes,&s.rigids,&s.accumulated,&s.attachments,&begin,&count,&dt};
            if(!s.launch(s.solveAttachments,count,args)) return 0;
        }
        for(size_t c=1;c<s.jointBatches.size();++c)
        {
            int begin=s.jointBatches[c-1],count=s.jointBatches[c]-begin;
            void* args[]={&s.rigids,&s.joints,&begin,&count,&dt};
            if(!s.launch(s.solveJoints,count,args)) return 0;
        }
        if((k+1)%config->iterations==0)
        {
            for(int v=0;v<config->velocityIterations && s.vertexCount;++v)
            {
                if(!s.check(cuMemsetD32Async(s.corrections,0,s.count*sizeof(btGpuXpbdVec)/sizeof(unsigned),s.stream))||!s.check(cuMemsetD32Async(s.references,0,s.count,s.stream))) return 0;
                if(!resetRigidCorrections(s)) return 0;
                int velocityOnly=1;
                void* countArgs[]={&s.nodes,&s.accumulated,&s.maps,&s.triangles,&s.contactCache,&s.contactCount,&s.contactCapacity,&s.references,&s.clearance,&velocityOnly,&dt};
                void* solveArgs[]={&s.nodes,&s.accumulated,&s.maps,&s.triangles,&s.contactCache,&s.contactCount,&s.contactCapacity,&s.references,&s.corrections,&s.clearance,&velocityOnly,&dt};
                void* applyArgs[]={&s.nodes,&s.corrections,&s.count};
                if(!s.launch(s.countContacts,16384,countArgs)||!s.launch(s.solveContacts,16384,solveArgs)||!s.launch(s.applyContactVelocities,s.count,applyArgs)) return 0;
                void* rigidApply[]={&s.rigids,&s.rigidCount,&dt,&velocityOnly};
                if(!s.launch(s.applyRigidContacts,s.rigidCount,rigidApply)) return 0;
            }
            void* settleArgs[]={&s.nodes,&s.accumulated,&s.count,&externalDt,&speed,&damping,&settlingBegin,&settlingEnd,&s.failure};
            if(!s.launch(s.settle,s.count,settleArgs)) return 0;
        }
    }
    void* rigidCheck[]={&s.rigids,&s.rigidCount,&s.failure};
    return s.launch(s.checkRigids,s.rigidCount,rigidCheck);
}
bool sameStep(const btGpuXpbdStep& a,const btGpuXpbdStep& b)
{
    // The step structure has alignment padding; comparing its bytes can spuriously rebuild the CUDA graph.
    return a.dt==b.dt&&a.substeps==b.substeps&&a.gravity.x==b.gravity.x&&a.gravity.y==b.gravity.y&&a.gravity.z==b.gravity.z&&
        a.settlingSpeed==b.settlingSpeed&&a.settlingDamping==b.settlingDamping&&a.iterations==b.iterations&&
        a.settlingBegin==b.settlingBegin&&a.settlingEnd==b.settlingEnd&&a.velocityIterations==b.velocityIterations&&a.materialMode==b.materialMode;
}
int step(void* value,const btGpuXpbdStep* config)
{
    auto& s=*static_cast<State*>(value);Current current(s);if(!current.valid) return 0;
    if(!config||config->materialMode<0||config->materialMode>1||s.count<=0||!std::isfinite(config->dt)||config->dt<=0||config->substeps<=0||config->substeps>4096||config->iterations<=0||config->iterations>255||config->velocityIterations<0||config->velocityIterations>255||
       !finite(config->gravity)||!std::isfinite(config->settlingSpeed)||config->settlingSpeed<0||
       !std::isfinite(config->settlingDamping)||config->settlingDamping<0||config->settlingBegin<0||config->settlingEnd<0||config->settlingEnd>s.count) { s.error="Invalid step configuration";return 0; }
    int failure=0;
    if(!s.check(cuMemcpyHtoD(s.failure,&failure,sizeof(failure)))) return 0;
    if(!s.graphExec || !sameStep(*config,s.graphConfig))
    {
        s.clearGraph();
        if(!s.check(cuStreamBeginCapture(s.stream,CU_STREAM_CAPTURE_MODE_THREAD_LOCAL))) return 0;
        bool captured=enqueueStep(s,config);
        CUresult end=cuStreamEndCapture(s.stream,&s.graph);
        if(!captured || !s.check(end)) { s.clearGraph();return 0; }
        if(!s.check(cuGraphInstantiateWithFlags(&s.graphExec,s.graph,0))) { s.clearGraph();return 0; }
        s.graphConfig=*config;
    }
    // Preserve the last accepted state on numerical or contact-capacity rejection.
    const size_t rigidBytes=size_t(s.rigidCount)*sizeof(btGpuXpbdRigid);
    if(rigidBytes&&!s.check(cuMemcpyDtoDAsync(s.acceptedRigids,s.rigids,rigidBytes,s.stream))) return 0;
    const size_t nodeBytes=size_t(s.count)*sizeof(btGpuXpbdNode);
    const size_t rotationBytes=size_t(s.batches.back())*4*sizeof(float);
    if(config->materialMode==1 && rotationBytes && !s.check(cuMemcpyDtoDAsync(s.acceptedRotations,s.rotations,rotationBytes,s.stream))) return 0;
    if(!s.check(cuMemcpyDtoDAsync(s.acceptedNodes,s.nodes,nodeBytes,s.stream))) return 0;
    if(!s.check(cuGraphLaunch(s.graphExec,s.stream)) || !s.check(cuStreamSynchronize(s.stream))) return 0;
    if(!s.check(cuMemcpyDtoH(&failure,s.failure,sizeof(failure)))) return 0;
    if(failure)
    {
        if(rigidBytes&&!s.check(cuMemcpyDtoDAsync(s.rigids,s.acceptedRigids,rigidBytes,s.stream))) return 0;
        if(config->materialMode==1 && rotationBytes && !s.check(cuMemcpyDtoDAsync(s.rotations,s.acceptedRotations,rotationBytes,s.stream))) return 0;
        if(!s.check(cuMemcpyDtoDAsync(s.nodes,s.acceptedNodes,nodeBytes,s.stream)) ||
           !s.check(cuStreamSynchronize(s.stream))) return 0;
        s.error=(failure==4?"Contact capacity exceeded":(failure==2?"Collision traversal stack overflow":"Nonfinite GPU state"));
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
int contacts(void* value,int pairCount,const btGpuXpbdNotificationPair* pairs,int* outputCount)
{
    auto& s=*static_cast<State*>(value);Current current(s);if(!current.valid) return 0;
    if(pairCount<0||(pairCount&&!pairs)||!outputCount|| (s.triangleCount && pairCount>2147483647/s.triangleCount))
        {s.error="Invalid GPU notification query";return 0;}
    s.lastReportCount=0;*outputCount=0;if(!s.vertexCount) return 1;
    if(!s.reports&&!s.check(cuMemAlloc(&s.reports,size_t(s.contactCapacity)*sizeof(btGpuXpbdContactReport)))) return 0;
    if(!s.reportCount&&!s.check(cuMemAlloc(&s.reportCount,sizeof(int)))) return 0;
    for(int i=0;i<pairCount;++i) if(pairs[i].bodyA<0||pairs[i].bodyB<0||pairs[i].bodyA==pairs[i].bodyB)
        {s.error="Invalid notification body IDs";return 0;}
    if(pairCount>s.notificationCapacity)
    {
        if(s.notificationPairs) cuMemFree(s.notificationPairs);s.notificationPairs=0;
        if(s.notificationWinners) cuMemFree(s.notificationWinners);s.notificationWinners=0;
        s.notificationCapacity=0;
        if(!s.check(cuMemAlloc(&s.notificationPairs,size_t(pairCount)*sizeof(*pairs)))||
           !s.check(cuMemAlloc(&s.notificationWinners,size_t(pairCount)*sizeof(unsigned long long)))) return 0;
        s.notificationCapacity=pairCount;
    }
    if(pairCount&&(!s.check(cuMemcpyHtoD(s.notificationPairs,pairs,size_t(pairCount)*sizeof(*pairs)))||
        !s.check(cuMemsetD32Async(s.notificationWinners,0xffffffff,2*pairCount,s.stream)))) return 0;
    if(!s.check(cuMemsetD32Async(s.failure,0,1,s.stream))||!s.check(cuMemsetD32Async(s.reportCount,0,1,s.stream))||
       !s.check(cuMemsetD32Async(s.accumulated,0,s.count*sizeof(btGpuXpbdVec)/sizeof(unsigned),s.stream))||
       !discoverSurfaceContacts(s)) return 0;
    void* reportArgs[]={&s.positions,&s.maps,&s.triangles,&s.vertexTriangles,&s.contactCache,&s.contactCount,&s.reports,&s.reportCount,&s.contactCapacity,&s.failure};
    if(!s.launch(s.reportContacts,s.contactCapacity,reportArgs)) return 0;
    void* notificationArgs[]={&s.maps,&s.positions,&s.triangles,&s.tree,&s.boxes,&s.notificationPairs,&pairCount,&s.triangleCount,&s.notificationWinners,&s.reports,&s.reportCount,&s.contactCapacity,&s.failure};
    if(!s.launch(s.reportNotifications,pairCount*s.triangleCount,notificationArgs)) return 0;
    void* finishArgs[]={&s.positions,&s.triangles,&s.notificationPairs,&pairCount,&s.notificationWinners,&s.reports,&s.reportCount,&s.contactCapacity,&s.failure};
    if(!s.launch(s.finishNotifications,pairCount,finishArgs)||!s.check(cuStreamSynchronize(s.stream))) return 0;
    int failure=0;
    if(!s.check(cuMemcpyDtoH(&failure,s.failure,sizeof(failure)))||!s.check(cuMemcpyDtoH(&s.lastReportCount,s.reportCount,sizeof(int)))) return 0;
    if(failure||s.lastReportCount>s.contactCapacity) {s.error=failure==4?"GPU callback contact capacity exceeded":failure==2?"GPU callback traversal stack overflow":"Nonfinite GPU callback contact witness";return 0;}
    *outputCount=s.lastReportCount;return 1;
}
int downloadContacts(void* value,int count,btGpuXpbdContactReport* output)
{
    auto& s=*static_cast<State*>(value);Current current(s);if(!current.valid) return 0;
    if(count!=s.lastReportCount||(count&&!output)) {s.error="Invalid GPU callback report download";return 0;}
    return !count||s.check(cuMemcpyDtoH(output,s.reports,size_t(count)*sizeof(*output)));
}
const btGpuXpbdTgsApi api={create,destroy,upload,forces,step,download,error,surface,surfacePositions,embed,targets,disabledPairs,updateNodes,friction,springs,rigids,downloadRigids,attachments,updateRigids,joints,selfCollision,contacts,downloadContacts};
}
extern "C" __declspec(dllexport) const btGpuXpbdTgsApi* btGetGpuXpbdTgsApi(int version)
{
    return version==BT_GPU_XPBD_ABI_VERSION?&api:nullptr;
}
