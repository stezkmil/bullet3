#include "btGpuXpbdTgsApi.h"
#include "device_types.h"
#include <tuple>
#include <type_traits>
#include <vector>

namespace
{
template<class T> struct Fields;
struct FieldsbtGpuXpbdVec { template<class T> static auto get(T& v) { return std::tie(v.x,v.y,v.z); } };
template<> struct Fields<btGpuXpbdVec> : FieldsbtGpuXpbdVec {};
template<> struct Fields<gpu::btGpuXpbdVec> : FieldsbtGpuXpbdVec {};
struct FieldsbtGpuXpbdNode { template<class T> static auto get(T& v) { return std::tie(v.position,v.inverseMass,v.velocity,v.damping); } };
template<> struct Fields<btGpuXpbdNode> : FieldsbtGpuXpbdNode {};
template<> struct Fields<gpu::btGpuXpbdNode> : FieldsbtGpuXpbdNode {};
struct FieldsbtGpuXpbdTet { template<class T> static auto get(T& v) { return std::tie(v.nodes,v.inverseRest,v.volume,v.young,v.poisson,v.damping); } };
template<> struct Fields<btGpuXpbdTet> : FieldsbtGpuXpbdTet {};
template<> struct Fields<gpu::btGpuXpbdTet> : FieldsbtGpuXpbdTet {};
struct FieldsbtGpuXpbdStep { template<class T> static auto get(T& v) { return std::tie(v.dt,v.substeps,v.gravity,v.settlingSpeed,v.settlingDamping,v.iterations,v.settlingBegin,v.settlingEnd,v.velocityIterations,v.materialMode); } };
template<> struct Fields<btGpuXpbdStep> : FieldsbtGpuXpbdStep {};
template<> struct Fields<gpu::btGpuXpbdStep> : FieldsbtGpuXpbdStep {};
struct FieldsbtGpuXpbdSurfaceVertex { template<class T> static auto get(T& v) { return std::tie(v.nodes,v.weights,v.offset,v.body); } };
template<> struct Fields<btGpuXpbdSurfaceVertex> : FieldsbtGpuXpbdSurfaceVertex {};
template<> struct Fields<gpu::btGpuXpbdSurfaceVertex> : FieldsbtGpuXpbdSurfaceVertex {};
struct FieldsbtGpuXpbdSurfaceTriangle { template<class T> static auto get(T& v) { return std::tie(v.vertices,v.part,v.index); } };
template<> struct Fields<btGpuXpbdSurfaceTriangle> : FieldsbtGpuXpbdSurfaceTriangle {};
template<> struct Fields<gpu::btGpuXpbdSurfaceTriangle> : FieldsbtGpuXpbdSurfaceTriangle {};
struct FieldsbtGpuXpbdContactReport { template<class T> static auto get(T& v) { return std::tie(v.bodyA,v.bodyB,v.partA,v.indexA,v.partB,v.indexB,v.pointA,v.pointB,v.normal,v.distance); } };
template<> struct Fields<btGpuXpbdContactReport> : FieldsbtGpuXpbdContactReport {};
template<> struct Fields<gpu::btGpuXpbdContactReport> : FieldsbtGpuXpbdContactReport {};
struct FieldsbtGpuXpbdContactImpulse { template<class T> static auto get(T& v) { return std::tie(v.contact,v.impulse); } };
template<> struct Fields<btGpuXpbdContactImpulse> : FieldsbtGpuXpbdContactImpulse {};
template<> struct Fields<gpu::btGpuXpbdContactImpulse> : FieldsbtGpuXpbdContactImpulse {};
struct FieldsbtGpuXpbdTarget { template<class T> static auto get(T& v) { return std::tie(v.node,v.position); } };
template<> struct Fields<btGpuXpbdTarget> : FieldsbtGpuXpbdTarget {};
template<> struct Fields<gpu::btGpuXpbdTarget> : FieldsbtGpuXpbdTarget {};
struct FieldsbtGpuXpbdDisabledPair { template<class T> static auto get(T& v) { return std::tie(v.bodyA,v.bodyB); } };
template<> struct Fields<btGpuXpbdDisabledPair> : FieldsbtGpuXpbdDisabledPair {};
template<> struct Fields<gpu::btGpuXpbdDisabledPair> : FieldsbtGpuXpbdDisabledPair {};
struct FieldsbtGpuXpbdNotificationPair { template<class T> static auto get(T& v) { return std::tie(v.bodyA,v.bodyB,v.touching); } };
template<> struct Fields<btGpuXpbdNotificationPair> : FieldsbtGpuXpbdNotificationPair {};
template<> struct Fields<gpu::btGpuXpbdNotificationPair> : FieldsbtGpuXpbdNotificationPair {};
struct FieldsbtGpuXpbdSpring { template<class T> static auto get(T& v) { return std::tie(v.node,v.target,v.stiffness,v.damping,v.maxForce); } };
template<> struct Fields<btGpuXpbdSpring> : FieldsbtGpuXpbdSpring {};
template<> struct Fields<gpu::btGpuXpbdSpring> : FieldsbtGpuXpbdSpring {};
struct FieldsbtGpuXpbdRigid { template<class T> static auto get(T& v) { return std::tie(v.position,v.orientation,v.linearVelocity,v.angularVelocity,v.inverseMass,v.inverseInertia,v.force,v.torque,v.linearDamping,v.angularDamping); } };
template<> struct Fields<btGpuXpbdRigid> : FieldsbtGpuXpbdRigid {};
template<> struct Fields<gpu::btGpuXpbdRigid> : FieldsbtGpuXpbdRigid {};
struct FieldsbtGpuXpbdEndpoint { template<class T> static auto get(T& v) { return std::tie(v.rigid,v.nodes,v.weights,v.point); } };
template<> struct Fields<btGpuXpbdEndpoint> : FieldsbtGpuXpbdEndpoint {};
template<> struct Fields<gpu::btGpuXpbdEndpoint> : FieldsbtGpuXpbdEndpoint {};
struct FieldsbtGpuXpbdAttachment { template<class T> static auto get(T& v) { return std::tie(v.a,v.b,v.compliance); } };
template<> struct Fields<btGpuXpbdAttachment> : FieldsbtGpuXpbdAttachment {};
template<> struct Fields<gpu::btGpuXpbdAttachment> : FieldsbtGpuXpbdAttachment {};
struct FieldsbtGpuXpbdJointAxis { template<class T> static auto get(T& v) { return std::tie(v.lower,v.upper,v.stiffness,v.damping,v.equilibrium,v.velocity,v.maxForce,v.servoTarget,v.bounce,v.stopERP,v.stopCFM,v.motorERP,v.motorCFM,v.flags); } };
template<> struct Fields<btGpuXpbdJointAxis> : FieldsbtGpuXpbdJointAxis {};
template<> struct Fields<gpu::btGpuXpbdJointAxis> : FieldsbtGpuXpbdJointAxis {};
struct FieldsbtGpuXpbdJoint { template<class T> static auto get(T& v) { return std::tie(v.bodies,v.points,v.orientations,v.axes,v.rotationOrder); } };
template<> struct Fields<btGpuXpbdJoint> : FieldsbtGpuXpbdJoint {};
template<> struct Fields<gpu::btGpuXpbdJoint> : FieldsbtGpuXpbdJoint {};
template<class D,class S> void convert(D& d,const S& s);
template<class D,class S,size_t... I> void convertFields(D d,S s,std::index_sequence<I...>)
{ (convert(std::get<I>(d),std::get<I>(s)),...); }
template<class D,class S> void convert(D& d,const S& s)
{
    if constexpr(std::is_arithmetic_v<D>) d=static_cast<D>(s);
    else if constexpr(std::is_array_v<D>) {for(size_t i=0;i<std::extent_v<D>;++i) convert(d[i],s[i]);}
    else
    {
        auto dst=Fields<D>::get(d);auto src=Fields<S>::get(s);
        convertFields(dst,src,std::make_index_sequence<std::tuple_size_v<decltype(dst)>>{});
        // FP32 accuracy: split/join the public FP64 node state without discarding GPU residuals.
        if constexpr(std::is_same_v<D,gpu::btGpuXpbdNode>)
        {
            d.positionResidual={float(s.position.x-d.position.x),float(s.position.y-d.position.y),float(s.position.z-d.position.z)};
            d.velocityResidual={float(s.velocity.x-d.velocity.x),float(s.velocity.y-d.velocity.y),float(s.velocity.z-d.velocity.z)};
        }
        else if constexpr(std::is_same_v<D,gpu::btGpuXpbdTet>)
        {
            // FP32 accuracy: derive reference edges from the rest matrix, never the possibly deformed input nodes.
            const auto cross=[](btGpuXpbdVec a,btGpuXpbdVec b)
            {return btGpuXpbdVec{a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};};
            btGpuXpbdVec a={s.inverseRest[0],s.inverseRest[1],s.inverseRest[2]},
                b={s.inverseRest[3],s.inverseRest[4],s.inverseRest[5]},c={s.inverseRest[6],s.inverseRest[7],s.inverseRest[8]};
            auto x=cross(b,c),y=cross(c,a),z=cross(a,b);double det=a.x*x.x+a.y*x.y+a.z*x.z;
            const btGpuXpbdVec edges[]={{x.x/det,y.x/det,z.x/det},{x.y/det,y.y/det,z.y/det},{x.z/det,y.z/det,z.z/det}};
            for(int j=0;j<3;++j)
            {
                convert(d.restEdges[j],edges[j]);
                d.restEdgeResiduals[j]={float(edges[j].x-d.restEdges[j].x),float(edges[j].y-d.restEdges[j].y),float(edges[j].z-d.restEdges[j].z)};
            }
        }
        else if constexpr(std::is_same_v<S,gpu::btGpuXpbdNode>)
        {
            d.position.x+=s.positionResidual.x;d.position.y+=s.positionResidual.y;d.position.z+=s.positionResidual.z;
            d.velocity.x+=s.velocityResidual.x;d.velocity.y+=s.velocityResidual.y;d.velocity.z+=s.velocityResidual.z;
        }
    }
}
template<class D,class S> std::vector<D> input(int count,const S* values)
{
    std::vector<D> result(values&&count>0?count:0);
    for(size_t i=0;i<result.size();++i) convert(result[i],values[i]);
    return result;
}
template<class D,class S> void output(D* values,const std::vector<S>& source)
{ if(values) for(size_t i=0;i<source.size();++i) convert(values[i],source[i]); }
#define INPUT(name,type) \
int name(void* state,int count,const type* values) \
{ auto data=input<gpu::type>(count,values);return gpuApi()->name(state,count,values?data.data():nullptr); }
#define OUTPUT(name,type) \
int name(void* state,int count,type* values) \
{ std::vector<gpu::type> data(values&&count>0?count:0);int ok=gpuApi()->name(state,count,values?data.data():nullptr);if(ok) output(values,data);return ok; }
INPUT(forces,btGpuXpbdVec)
INPUT(targets,btGpuXpbdTarget)
INPUT(disabledPairs,btGpuXpbdDisabledPair)
INPUT(updateNodes,btGpuXpbdNode)
INPUT(springs,btGpuXpbdSpring)
INPUT(rigids,btGpuXpbdRigid)
INPUT(attachments,btGpuXpbdAttachment)
INPUT(updateRigids,btGpuXpbdRigid)
INPUT(joints,btGpuXpbdJoint)
OUTPUT(download,btGpuXpbdNode)
OUTPUT(surfacePositions,btGpuXpbdVec)
OUTPUT(downloadRigids,btGpuXpbdRigid)
OUTPUT(downloadContacts,btGpuXpbdContactReport)
#undef INPUT
#undef OUTPUT
int upload(void* state,int count,const btGpuXpbdNode* nodes,int tetCount,const btGpuXpbdTet* tets)
{
    auto n=input<gpu::btGpuXpbdNode>(count,nodes);auto t=input<gpu::btGpuXpbdTet>(tetCount,tets);
    return gpuApi()->upload(state,count,nodes?n.data():nullptr,tetCount,tets?t.data():nullptr);
}
int step(void* state,const btGpuXpbdStep* config)
{
    gpu::btGpuXpbdStep value={};if(config) convert(value,*config);
    return gpuApi()->step(state,config?&value:nullptr);
}
int surface(void* state,int vertexCount,const btGpuXpbdSurfaceVertex* vertices,int triangleCount,
    const btGpuXpbdSurfaceTriangle* triangles,double margin,double clearance)
{
    auto v=input<gpu::btGpuXpbdSurfaceVertex>(vertexCount,vertices);auto t=input<gpu::btGpuXpbdSurfaceTriangle>(triangleCount,triangles);
    return gpuApi()->surface(state,vertexCount,vertices?v.data():nullptr,triangleCount,triangles?t.data():nullptr,
        static_cast<float>(margin),static_cast<float>(clearance));
}
int embed(void* state,int count,const btGpuXpbdVec* points,int begin,int end,int body,btGpuXpbdSurfaceVertex* maps)
{
    auto p=input<gpu::btGpuXpbdVec>(count,points);std::vector<gpu::btGpuXpbdSurfaceVertex> m(maps&&count>0?count:0);
    int ok=gpuApi()->embed(state,count,points?p.data():nullptr,begin,end,body,maps?m.data():nullptr);
    if(ok) output(maps,m);return ok;
}
int friction(void* state,int count,const double* coefficients)
{
    auto data=input<float>(count,coefficients);return gpuApi()->friction(state,count,coefficients?data.data():nullptr);
}
int contacts(void* state,int count,const btGpuXpbdNotificationPair* pairs,int* result)
{
    auto data=input<gpu::btGpuXpbdNotificationPair>(count,pairs);return gpuApi()->contacts(state,count,pairs?data.data():nullptr,result);
}
int contactImpulses(void* state,int* count,btGpuXpbdContactImpulse* values,int capacity)
{
    std::vector<gpu::btGpuXpbdContactImpulse> data(values&&capacity>0?capacity:0);
    int ok=gpuApi()->contactImpulses(state,count,values?data.data():nullptr,capacity);
    if(ok&&values&&count) {data.resize(*count);output(values,data);}return ok;
}
const btGpuXpbdTgsApi api={gpuApi()->create,gpuApi()->destroy,upload,forces,step,download,gpuApi()->error,
    surface,surfacePositions,embed,targets,disabledPairs,updateNodes,friction,springs,rigids,downloadRigids,attachments,
    updateRigids,joints,gpuApi()->selfCollision,contacts,downloadContacts,gpuApi()->enableContactImpulses,contactImpulses};
}
extern "C" __declspec(dllexport) const btGpuXpbdTgsApi* btGetGpuXpbdTgsApi(int version)
{ return version==BT_GPU_XPBD_ABI_VERSION?&api:nullptr; }
