#include "btGpuXpbdTgsApi.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>
#ifdef XPBD_TEST_RUNTIME_LOADER
#include "btGpuXpbdTgsBackend.h"
#else
extern "C" __declspec(dllimport) const btGpuXpbdTgsApi* btGetGpuXpbdTgsApi(int);
#endif
namespace
{
int failures=0;
#ifdef XPBD_TEST_FP32
constexpr double scalarRoundoff=2.e-6;
using SolverScalar=float;
#else
constexpr double scalarRoundoff=0;
using SolverScalar=double;
#endif
double tolerance(double strict,double scale=1) { return strict>scalarRoundoff*scale?strict:scalarRoundoff*scale; }

int (*accurateStep)(void*,const btGpuXpbdStep*)=nullptr;
int mediumStep(void* world,const btGpuXpbdStep* config)
{
    auto step=*config;step.materialMode=1;return accurateStep(world,&step);
}

void require(bool value,const char* label) { if(!value) { std::fprintf(stderr,"FAIL: %s\n",label);++failures; } }
double distance(btGpuXpbdVec a,btGpuXpbdVec b)
{ return std::sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z)); }
btGpuXpbdRigid rigidState(double inverseMass=1)
{
    btGpuXpbdRigid r={};r.orientation[3]=1;r.inverseMass=inverseMass;
    r.inverseInertia={inverseMass,inverseMass,inverseMass};return r;
}
btGpuXpbdEndpoint softEndpoint(int node)
{
    btGpuXpbdEndpoint e={};e.rigid=-1;for(int& n:e.nodes) n=-1;e.nodes[0]=node;e.weights[0]=node>=0?1:0;return e;
}
btGpuXpbdEndpoint rigidEndpoint(int rigid,btGpuXpbdVec point={0,0,0})
{
    auto e=softEndpoint(-1);e.rigid=rigid;e.point=point;return e;
}
void testCompensatedState(const btGpuXpbdTgsApi* api,void* world)
{
    // FP32 accuracy: both motion and impulses are smaller than the high part's ULP per increment.
    btGpuXpbdNode node={{1.e6,0,0},1,{.001,0,0},0};
    btGpuXpbdVec force={1.e-5,0,0};
    btGpuXpbdStep step={.01,1,{0,0,0},0,0,50,0,0,0};
    require(api->upload(world,1,&node,0,nullptr)&&api->forces(world,1,&force),"compensated state fixture");
    for(int i=0;i<100;++i)
        require(api->step(world,&step)&&api->download(world,1,&node)&&api->updateNodes(world,1,&node),
            "compensated state survives application repacking");
    const double h=step.dt/step.iterations,increments=100*step.iterations;
    const double expected=1.e6+increments*h*.001+force.x*h*h*increments*(increments+1)/2;
    require(std::fabs(node.position.x-expected)<5.e-9,"sub-ULP displacement survives large-coordinate API roundtrips");
    require(std::fabs(node.velocity.x-(.001+force.x))<1.e-10,"sub-ULP impulses survive API roundtrips");
}
void testContactReports(const btGpuXpbdTgsApi* api,void* world)
{
    btGpuXpbdNode node={};
    require(api->upload(world,1,&node,0,nullptr)!=0,"report fixture upload");
    const btGpuXpbdVec points[]={{-1,-1,0},{1,-1,0},{0,1,0},{0,-.5,-1},{0,-.5,1},{0,.5,0}};
    btGpuXpbdSurfaceVertex vertices[6]={};
    for(int i=0;i<6;++i)
    {
        for(int& n:vertices[i].nodes) n=-1;
        vertices[i].offset=points[i];vertices[i].body=i/3;
    }
    btGpuXpbdSurfaceTriangle triangles[]={{{0,1,2},7,41},{{3,4,5},9,83},{{3,4,5},9,84}};
    btGpuXpbdDisabledPair pair={0,1};btGpuXpbdNotificationPair notification={0,1,0};int count=-1;
    require(api->surface(world,6,vertices,3,triangles,.1,0)&&api->disabledPairs(world,1,&pair),"notification report fixture");
    require(api->contacts(world,1,&notification,&count)&&count==1,"GPU notification first-hit query");
    btGpuXpbdContactReport report={};
    if(count==1)
    {
        require(api->downloadContacts(world,1,&report)!=0,"GPU report download");
        require(report.bodyA==0&&report.bodyB==1&&report.partA==7&&report.indexA==41&&report.partB==9&&report.indexB==83,"original mesh triangle provenance survives GPU query");
        require(report.distance<0&&std::isfinite(report.normal.x),"penetrating report has finite signed geometry");
        for(int i=0;i<4;++i)
        {
            btGpuXpbdContactReport repeat={};
            require(api->contacts(world,1,&notification,&count)&&count==1&&api->downloadContacts(world,1,&repeat)&&
                std::memcmp(&repeat,&report,sizeof(report))==0,"first-hit report is independent of CUDA scheduling");
        }
    }
    require(api->contacts(world,0,nullptr,&count)&&count==0,"notification reports do not persist without a requested pair");
    for(int i=3;i<6;++i) vertices[i].offset.x+=10;
    require(api->surface(world,6,vertices,2,triangles,.1,0)&&api->contacts(world,1,&notification,&count)&&count==0,"separated notification pair loses contact");
    require(api->download(world,1,&node)&&distance(node.position,{0,0,0})==0,"report query never advances simulation");
    node.position={0,0,-.01};node.inverseMass=1;
    require(api->updateNodes(world,1,&node)!=0,"response report node update");
    const btGpuXpbdVec responsePoints[]={{-.1,-.1,0},{.1,-.1,0},{0,.1,0},{-1,-1,0},{1,-1,0},{0,1,0}};
    for(int i=0;i<6;++i)
    {
        vertices[i].offset=responsePoints[i];
        if(i<3) {vertices[i].nodes[0]=0;vertices[i].weights[0]=1;}
    }
    require(api->surface(world,6,vertices,2,triangles,.1,0)&&api->contacts(world,0,nullptr,&count)&&count>0,"response contacts reuse GPU narrowphase data");
    std::vector<btGpuXpbdContactReport> response(count);
    require(api->downloadContacts(world,count,response.data())!=0,"response report download");
    for(const auto& contact:response)
        require(contact.partA==7&&contact.indexA==41&&contact.partB==9&&contact.indexB==83&&contact.distance<0,
            "vertex response reports retain original incident triangle IDs and penetration");
    std::vector<btGpuXpbdSurfaceVertex> overlapping(vertices,vertices+6);
    for(int i=3;i<6;++i) {auto vertex=vertices[i];vertex.body=2;overlapping.push_back(vertex);}
    btGpuXpbdSurfaceTriangle interleaved[]={{{3,4,5},9,83},{{0,1,2},7,41},{{6,7,8},11,91},{{3,4,5},9,84}};
    require(api->surface(world,9,overlapping.data(),4,interleaved,.1,0)&&api->disabledPairs(world,1,&pair)&&
        api->contacts(world,0,nullptr,&count)&&count>0,"body BVH retains allowed overlapping pair");
    response.resize(count);
    require(api->downloadContacts(world,count,response.data())!=0,"body BVH filtered reports");
    for(const auto& contact:response)
        require(contact.bodyA==0&&contact.bodyB==2,"ignored body subtree produces no contacts while overlapping allowed body does");
}
void testCallbackContactReuse(const btGpuXpbdTgsApi* api,void* world)
{
    for(bool rigid: {false,true}) for(int substeps: {1,2}) for(int edit=0;edit<9;++edit)
    {
        btGpuXpbdVec expectedPosition={},expectedVelocity={};
        for(bool query: {false,true})
        {
            btGpuXpbdNode node={{0,.001,0},1,{0,-1,0},0};
            auto body=rigidState();body.position=node.position;body.linearVelocity=node.velocity;
            btGpuXpbdSurfaceVertex maps[]={{{rigid?-2:0,-1,-1,-1},{rigid?0.:1.,0,0,0},{0,0,0},0},
                {{-1,-1,-1,-1},{0,0,0,0},{-2,0,-2},1},{{-1,-1,-1,-1},{0,0,0,0},{0,0,2},1},{{-1,-1,-1,-1},{0,0,0,0},{2,0,-2},1}};
            btGpuXpbdSurfaceTriangle triangle={{1,2,3}};
            if(edit==8) maps[0].offset.x=.5;
            btGpuXpbdStep step={.01,substeps,{0,-10,0},0,0,10,0,0,1};
            require(api->upload(world,1,&node,0,nullptr)&&(!rigid||api->rigids(world,1,&body))&&
                api->surface(world,4,maps,1,&triangle,.1,0),"callback reuse fixture");
            btGpuXpbdDisabledPair ignored={0,1};
            if(edit==5) require(api->disabledPairs(world,1,&ignored)!=0,"initial callback exclusion");
            int count=0;
            if(query) require(api->contacts(world,0,nullptr,&count)!=0,"callback discovery before edits");
            if(edit==1) {node.velocity.x=1;body.linearVelocity.x=1;body.force.x=2;}
            if(edit==2) {node.position.y=1;body.position.y=1;}
            if(edit==3) {node.inverseMass=0;body=rigidState(0);body.position=node.position;}
            if(edit==8&&rigid) {body.orientation[2]=std::sqrt(.5);body.orientation[3]=std::sqrt(.5);}
            require(api->updateNodes(world,1,&node)&&(!rigid||api->updateRigids(world,1,&body)),"callback state repack");
            if(edit==4) require(api->disabledPairs(world,1,&ignored)!=0,"callback disables pair");
            if(edit==5) require(api->disabledPairs(world,0,nullptr)!=0,"callback enables pair");
            if(edit==6)
            {
                double coefficients[]={.5,.5};
                require(api->friction(world,2,coefficients)!=0,"callback changes friction");
            }
            if(edit==7&&!rigid)
            {
                btGpuXpbdTarget target={0,{0,.1,0}};
                require(api->targets(world,1,&target)!=0,"callback sets target");
            }
            // A second step must refresh contacts even when the first reused a callback query.
            require(api->step(world,&step)&&api->step(world,&step)&&api->download(world,1,&node)&&
                (!rigid||api->downloadRigids(world,1,&body)),"callback reuse and following step");
            auto position=rigid?body.position:node.position,velocity=rigid?body.linearVelocity:node.velocity;
            if(edit==0) require(position.y>=-1.e-8&&velocity.y>=-1.e-8,"reused contacts stop the falling body at the plane");
            if(!query) {expectedPosition=position;expectedVelocity=velocity;}
            else require(distance(position,expectedPosition)<tolerance(1.e-9)&&distance(velocity,expectedVelocity)<tolerance(1.e-9),
                "callback discovery reuse or invalidation matches fresh discovery");
        }
    }
}
void testContactImpulses(const btGpuXpbdTgsApi* api,void* world)
{
    for(int substeps: {1,4})
    {
        btGpuXpbdNode initial={{0,0,0},1,{0,0,0},0},baseline={};
        btGpuXpbdSurfaceVertex maps[]={{{0,-1,-1,-1},{1,0,0,0},{0,0,0},0},
            {{-1,-1,-1,-1},{0,0,0,0},{-2,0,-2},1},{{-1,-1,-1,-1},{0,0,0,0},{0,0,2},1},{{-1,-1,-1,-1},{0,0,0,0},{2,0,-2},1}};
        btGpuXpbdSurfaceTriangle triangle={{1,2,3},7,42};
        btGpuXpbdStep step={.01,substeps,{0,-10,0},0,0,10,0,0,1};
        require(api->upload(world,1,&initial,0,nullptr)&&api->surface(world,4,maps,1,&triangle,.1,0),"impulse fixture");
        require(api->step(world,&step)&&api->download(world,1,&baseline),"reporting disabled baseline");
        int count=-1;
        require(api->contactImpulses(world,&count,nullptr,0)&&count==0,"reporting disabled by default");
        require(api->updateNodes(world,1,&initial)&&api->enableContactImpulses(world,1)&&api->step(world,&step),"reported step");
        btGpuXpbdNode reported={};require(api->download(world,1,&reported)!=0,"reported state");
        require(std::memcmp(&reported,&baseline,sizeof(reported))==0,"reporting does not change physical state");
        require(api->contactImpulses(world,&count,nullptr,0)&&count==1,"substeps aggregate one contact");
        btGpuXpbdContactImpulse impulse={};
        if(count==1)
        {
            require(api->contactImpulses(world,&count,&impulse,1)!=0,"impulse readback");
            require(std::fabs(impulse.impulse-.1)<tolerance(1.e-8)&&impulse.contact.normal.y>.999,"reported impulse balances gravity across all substeps");
            require(impulse.contact.bodyA==0&&impulse.contact.bodyB==1&&impulse.contact.indexB==42,"impulse body and triangle identity");
        }
        auto invalid=step;invalid.dt=0;
        require(!api->step(world,&invalid)&&api->contactImpulses(world,&count,nullptr,0)&&count==0,"failed step clears stale reports");
        require(api->enableContactImpulses(world,0)&&api->step(world,&step)&&api->contactImpulses(world,&count,nullptr,0)&&count==0,"disable clears reports and captured reporting work");
    }
}
void testRigidCoupling(const btGpuXpbdTgsApi* api,void* world)
{
    btGpuXpbdNode node={};btGpuXpbdStep step={.1,2,{0,0,0},0,0,10,0,0,0,0};
    require(api->upload(world,1,&node,0,nullptr)!=0,"rigid fixture upload");
    auto r=rigidState(.5);r.inverseInertia={1,2,3};r.force={4,0,0};r.torque={0,0,2};
    require(api->rigids(world,1,&r)&&api->step(world,&step)&&api->downloadRigids(world,1,&r),"rigid force and torque integration");
    require(std::fabs(r.linearVelocity.x-.2)<tolerance(1.e-12),"rigid force impulse");
    require(std::fabs(r.angularVelocity.z-.6)<tolerance(1.e-11),"rigid torque impulse");
    require(std::fabs(r.orientation[0]*r.orientation[0]+r.orientation[1]*r.orientation[1]+r.orientation[2]*r.orientation[2]+r.orientation[3]*r.orientation[3]-1)<tolerance(1.e-12),"rigid normalized orientation");
    // An off-centre anchor must partition a node impulse into rigid translation and spin.
    node={ {0,1,0},1,{1,0,0},0 };r=rigidState();step.dt=1.e-5;step.substeps=step.iterations=1;
    require(api->upload(world,1,&node,0,nullptr)&&api->rigids(world,1,&r),"off-centre setup");
    btGpuXpbdAttachment anchor={softEndpoint(0),rigidEndpoint(0,{0,1,0}),0};
    require(api->attachments(world,1,&anchor)&&api->step(world,&step)&&api->download(world,1,&node)&&api->downloadRigids(world,1,&r),"off-centre GPU attachment");
    require(std::fabs(node.velocity.x-2./3)<tolerance(1.e-9)&&std::fabs(r.linearVelocity.x-1./3)<tolerance(1.e-9)&&std::fabs(r.angularVelocity.z+1./3)<tolerance(1.e-8),"off-centre effective mass includes inertia");
    require(std::fabs(node.velocity.x+r.linearVelocity.x-1)<tolerance(1.e-12),"rigid-soft linear momentum");
    double angular=node.position.x*node.velocity.y-node.position.y*node.velocity.x+r.position.x*r.linearVelocity.y-r.position.y*r.linearVelocity.x+r.angularVelocity.z;
    require(std::fabs(angular+1)<tolerance(1.e-8),"rigid-soft angular momentum");
    // A fixed endpoint and a soft-soft attachment use the same generalized solve.
    btGpuXpbdNode pair[2]={{{0,0,0},1,{1,0,0},0},{{0,0,0},.5,{0,0,0},0}};
    anchor={softEndpoint(0),softEndpoint(1),0};
    require(api->upload(world,2,pair,0,nullptr)&&api->attachments(world,1,&anchor)&&api->step(world,&step)&&api->download(world,2,pair),"soft-soft attachment");
    require(std::fabs(pair[0].velocity.x-1./3)<tolerance(1.e-12)&&std::fabs(pair[1].velocity.x-1./3)<tolerance(1.e-12),"soft-soft attachment impulse balance");
    anchor={softEndpoint(0),softEndpoint(-1),0};node={ {0,0,0},1,{1,0,0},0 };
    require(api->upload(world,1,&node,0,nullptr)&&api->attachments(world,1,&anchor)&&api->step(world,&step)&&api->download(world,1,&node),"world attachment");
    require(std::fabs(node.position.x)<tolerance(1.e-12)&&std::fabs(node.velocity.x)<tolerance(1.e-12),"world attachment holds node");
    for(int i=1;i<=3;++i)
    {
        anchor.b.point.x=.01*i;
        require(api->attachments(world,1,&anchor)&&api->step(world,&step)&&api->download(world,1,&node),"update captured attachment target");
        require(std::fabs(node.position.x-anchor.b.point.x)<tolerance(1.e-10),"captured attachment follows updated target");
    }
    // Rejected steps must roll back both sides of a coupled island.
    r=rigidState();r.force={std::numeric_limits<SolverScalar>::max(),0,0};r.inverseMass=100;
    require(api->rigids(world,1,&r)!=0,"finite rigid overflow input");require(api->downloadRigids(world,1,&r)!=0,"quantized rigid input");auto accepted=r;
    step.dt=1;require(!api->step(world,&step),"rigid overflow rejection");
    require(api->downloadRigids(world,1,&r)&&std::memcmp(&r,&accepted,sizeof(r))==0,"rigid rollback exact state");
    r=rigidState();r.orientation[3]=2;require(!api->rigids(world,1,&r),"invalid rigid orientation rejected");
}
void testRigidContacts(const btGpuXpbdTgsApi* api,void* world)
{
    btGpuXpbdNode node={{0,.001,0},1,{0,-1,0},0};auto r=rigidState();
    btGpuXpbdSurfaceVertex maps[]={{{0,-1,-1,-1},{1,0,0,0},{0,0,0},0},
        {{-2,-1,-1,-1},{0,0,0,0},{-2,0,-2},1},{{-2,-1,-1,-1},{0,0,0,0},{0,0,2},1},{{-2,-1,-1,-1},{0,0,0,0},{2,0,-2},1}};
    btGpuXpbdSurfaceTriangle triangle={{1,2,3}};
    require(api->upload(world,1,&node,0,nullptr)&&api->rigids(world,1,&r)&&api->surface(world,4,maps,1,&triangle,.1,0),"rigid-soft contact upload");
    btGpuXpbdStep step={.01,2,{0,0,0},0,0,20,0,0,1};
    require(api->step(world,&step)&&api->download(world,1,&node)&&api->downloadRigids(world,1,&r),"rigid-soft contact step");
    require(r.linearVelocity.y<-.1,"soft contact pushes dynamic rigid");
    require(std::fabs(node.velocity.y+r.linearVelocity.y+1)<tolerance(1.e-9),"rigid-soft contact conserves momentum");
    require(node.position.y-r.position.y> -1.e-5,"rigid-soft contact limits penetration");
    // A rigid point against a fixed plane exercises the rigid-only contact path.
    node={};r=rigidState();r.position={0,.001,0};r.linearVelocity={0,-1,0};maps[0]={{-2,-1,-1,-1},{0,0,0,0},{0,0,0},0};
    for(int i=1;i<4;++i) maps[i].nodes[0]=-1;
    require(api->upload(world,1,&node,0,nullptr)&&api->rigids(world,1,&r)&&api->surface(world,4,maps,1,&triangle,.1,0),"rigid ground upload");
    require(api->step(world,&step)&&api->downloadRigids(world,1,&r),"rigid ground step");
    require(r.position.y> -1.e-6&&r.linearVelocity.y> -1.e-5,"rigid ground contact stops fall");
    // Tangential inertia differs from normal inertia for a contact below the COM.
    r=rigidState();r.position={0,1,0};r.linearVelocity={1,0,0};maps[0].offset={0,-1,0};
    require(api->rigids(world,1,&r)&&api->surface(world,4,maps,1,&triangle,.1,0),"rigid friction upload");
    double friction[2]={10,10};require(api->friction(world,2,friction)!=0,"rigid friction coefficient");
    step={.001,1,{0,-100,0},0,0,1,0,0,0};
    require(api->step(world,&step)&&api->downloadRigids(world,1,&r),"rigid friction step");
    std::printf("rigid friction vx=%g vy=%g omegaZ=%g\n",r.linearVelocity.x,r.linearVelocity.y,r.angularVelocity.z);
    require(std::fabs(r.linearVelocity.x-.5)<1.e-6&&std::fabs(r.angularVelocity.z+.5)<1.e-6,"rigid friction uses tangential effective inertia");
    require(std::fabs(r.linearVelocity.x+r.angularVelocity.z)<1.e-6,"rigid contact slip removed without reversing it");
    // Replacing rigid IDs must drop old surfaces; state-only updates must retain them.
    auto malformed=maps[0];maps[0].nodes[1]=0;maps[0].weights[1]=1;
    require(!api->surface(world,4,maps,1,&triangle,.1,0),"mixed rigid-soft surface support rejected");maps[0]=malformed;
    r=rigidState();r.position={0,1,0};r.linearVelocity={0,-1,0};
    require(api->updateRigids(world,1,&r)&&api->step(world,&step)&&api->downloadRigids(world,1,&r),"rigid state update retains contacts");
    require(r.linearVelocity.y>-.001,"updated rigid is stopped by retained surface");
    require(api->rigids(world,0,nullptr)&&api->step(world,&step),"removing rigids clears dependent surface references");
}
btGpuXpbdJoint freeJoint(int a,int b)
{
    btGpuXpbdJoint j={};j.bodies[0]=a;j.bodies[1]=b;j.orientations[0][3]=j.orientations[1][3]=1;
    for(auto& axis:j.axes){axis.lower=1;axis.upper=-1;axis.stopERP=1;axis.motorERP=1;}return j;
}
void testRigidJoints(const btGpuXpbdTgsApi* api,void* world)
{
    btGpuXpbdNode node={};require(api->upload(world,1,&node,0,nullptr)!=0,"joint fixture");
    auto r=rigidState();r.force={2,3,4};r.torque={1,2,3};auto j=freeJoint(-1,0);
    for(auto& axis:j.axes) axis.lower=axis.upper=0;
    btGpuXpbdStep step={.01,1,{0,0,0},0,0,20,0,0,0};
    require(api->rigids(world,1,&r)&&api->joints(world,1,&j)&&api->step(world,&step)&&api->downloadRigids(world,1,&r),"fixed GPU joint");
    require(distance(r.position,{0,0,0})<tolerance(1.e-10)&&distance(r.linearVelocity,{0,0,0})<tolerance(1.e-8),"fixed joint holds position under force");
    std::printf("fixed joint residual angular speed: %.9g\n",distance(r.angularVelocity,{0,0,0}));
    require(distance(r.angularVelocity,{0,0,0})<1.e-7,"fixed joint holds rotation under torque");
    for(int i=1;i<=3;++i)
    {
        j.points[0].x=.01*i;
        require(api->joints(world,1,&j)&&api->step(world,&step)&&api->downloadRigids(world,1,&r),"update captured mouse joint frame");
        require(std::fabs(r.position.x-j.points[0].x)<tolerance(1.e-9),"captured joint follows updated frame");
    }
    require(api->joints(world,0,nullptr)&&api->step(world,&step)&&api->downloadRigids(world,1,&r),"remove captured joint");
    require(r.position.x>j.points[0].x,"removed joint releases body");
    require(api->joints(world,1,&j)&&api->step(world,&step)&&api->downloadRigids(world,1,&r),"restore captured joint");
    require(std::fabs(r.position.x-j.points[0].x)<tolerance(1.e-9),"restored joint holds body");
    r=rigidState(.5);j=freeJoint(-1,0);j.axes[0].flags=2;j.axes[0].velocity=1;j.axes[0].maxForce=.4;step.dt=.1;
    require(api->rigids(world,1,&r)&&api->joints(world,1,&j)&&api->step(world,&step)&&api->downloadRigids(world,1,&r),"force limited joint motor");
    require(std::fabs(r.linearVelocity.x-.02)<tolerance(1.e-10),"motor impulse equals force times time");
    r=rigidState();r.force={2,0,0};j=freeJoint(-1,0);j.axes[0].flags=1;j.axes[0].stiffness=100;j.axes[0].damping=20;step.dt=.02;
    require(api->rigids(world,1,&r)&&api->joints(world,1,&j),"spring joint upload");
    for(int i=0;i<100;++i) require(api->step(world,&step)!=0,"spring joint step");
    require(api->downloadRigids(world,1,&r)&&std::fabs(r.position.x-.02)<1.e-6,"joint spring reaches force over stiffness equilibrium");
    r=rigidState();r.linearVelocity={1,0,0};j=freeJoint(-1,0);j.axes[0].lower=-.01;j.axes[0].upper=.01;
    require(api->rigids(world,1,&r)&&api->joints(world,1,&j)&&api->step(world,&step)&&api->downloadRigids(world,1,&r),"joint upper stop");
    std::printf("joint limit position error %.9g velocity %.9g\n",r.position.x-.01,r.linearVelocity.x);
    // FP32 accuracy: the residual velocity is two float ULPs at unit speed.
    require(std::fabs(r.position.x-.01)<tolerance(1.e-8)&&std::fabs(r.linearVelocity.x)<tolerance(1.e-7),"joint limit stops at its bound");
    r=rigidState();j=freeJoint(-1,0);j.axes[0].flags=6;j.axes[0].velocity=1;j.axes[0].maxForce=100;j.axes[0].servoTarget=.05;
    require(api->rigids(world,1,&r)&&api->joints(world,1,&j),"joint servo setup");
    for(int i=0;i<20;++i) require(api->step(world,&step)!=0,"joint servo step");
    require(api->downloadRigids(world,1,&r)&&std::fabs(r.position.x-.05)<1.e-7,"joint servo reaches target without drift");
    btGpuXpbdRigid pair[2]={rigidState(),rigidState(.5)};pair[0].linearVelocity={1,0,0};j=freeJoint(0,1);j.axes[0].lower=j.axes[0].upper=0;
    require(api->rigids(world,2,pair)&&api->joints(world,1,&j)&&api->step(world,&step)&&api->downloadRigids(world,2,pair),"two body joint");
    require(std::fabs(pair[0].linearVelocity.x-1./3)<tolerance(1.e-9)&&std::fabs(pair[1].linearVelocity.x-1./3)<tolerance(1.e-9),"joint transfers momentum by physical masses");
    pair[0]=rigidState(0);pair[0].linearVelocity={1,0,0};pair[0].angularVelocity={0,0,1};pair[1]=rigidState(.5);
    j=freeJoint(0,1);for(auto& a:j.axes)a.lower=a.upper=0;
    require(api->rigids(world,2,pair)&&api->joints(world,1,&j)&&api->step(world,&step)&&api->downloadRigids(world,2,pair),"kinematic core joint");
    require(std::fabs(pair[0].position.x-step.dt)<tolerance(1.e-10)&&std::fabs(pair[1].position.x-step.dt)<1.e-7,"kinematic core joint follows translation");
}
void testSelfCollision(const btGpuXpbdTgsApi* api,void* world)
{
    btGpuXpbdNode nodes[6]={{{-2,.001,-2},1,{0,-1,0},0},{{2,.001,-2},1,{0,-1,0},0},{{0,.001,2},1,{0,-1,0},0},
        {{-2,0,-2},1,{0,0,0},0},{{0,0,2},1,{0,0,0},0},{{2,0,-2},1,{0,0,0},0}};
    btGpuXpbdSurfaceVertex maps[6]={};for(int i=0;i<6;++i) {for(int& n:maps[i].nodes)n=-1;maps[i].nodes[0]=i;maps[i].weights[0]=1;}
    btGpuXpbdSurfaceTriangle triangles[2]={{{0,1,2}},{{3,4,5}}};
    btGpuXpbdStep step={.01,4,{0,0,0},0,0,20,0,0,2};int self=0;
    for(int enabled=0;enabled<2;++enabled)
    {
        require(api->upload(world,6,nodes,0,nullptr)&&api->surface(world,6,maps,2,triangles,.1,0)&&api->selfCollision(world,enabled,&self),"self collision surface setup");
        btGpuXpbdNode result[6];require(api->step(world,&step)&&api->download(world,6,result),"self collision step");
        double momentum=0,lower=0;for(int i=0;i<6;++i){momentum+=result[i].velocity.y;if(i>=3)lower+=result[i].velocity.y;}
        require(std::fabs(momentum+3)<tolerance(1.e-8),"self collision preserves total momentum");
        require(enabled?lower<-.1:std::fabs(lower)<tolerance(1.e-12),"self collision enabled only for requested bodies");
    }
    triangles[1]={{0,2,1}};
    require(api->upload(world,3,nodes,0,nullptr)&&api->surface(world,3,maps,2,triangles,.1,0)&&api->selfCollision(world,1,&self),"adjacent self collision setup");
    btGpuXpbdNode result[3];require(api->step(world,&step)&&api->download(world,3,result),"adjacent self collision step");
    require(std::fabs(result[0].velocity.y+1)<tolerance(1.e-12),"shared mesh support does not self collide");
}
void testMaterialModeSwitch(const btGpuXpbdTgsApi* api,void* world)
{
    btGpuXpbdNode nodes[4]={{{1.e6,0,0},4,{0,0,0},0},{{1.e6+1,0,0},4,{0,0,0},0},
        {{1.e6,1,0},4,{0,0,0},0},{{1.e6,0,1},4,{0,0,0},0}};
    btGpuXpbdTet tet={{0,1,2,3},{1,0,0,0,1,0,0,0,1},1./6,7.e7,.4,.01};
    btGpuXpbdVec force[4]={{1,0,0},{1,0,0},{1,0,0},{1,0,0}};
    require(api->upload(world,4,nodes,1,&tet)&&api->forces(world,4,force),"mode switch upload");
    btGpuXpbdStep step={.0002,1,{0,0,0},0,0,50,0,0,0,0};
    for(int mode:{0,1,0,1}) { step.materialMode=mode;require(api->step(world,&step)!=0,"switch material mode on live state"); }
    require(api->download(world,4,nodes)!=0,"mode switch download");
    double momentum=0;for(const auto& n:nodes) momentum+=n.velocity.x/n.inverseMass;
    require(std::fabs(momentum-.0032)<tolerance(1.e-8),"mode switching preserves applied impulse");
    require(nodes[0].position.x>1.e6+1.e-7,"Medium retains small translation at large coordinates");
    auto accepted=std::vector<btGpuXpbdNode>(nodes,nodes+4);
    force[0].x=std::numeric_limits<SolverScalar>::max();
    require(api->forces(world,4,force)!=0,"Medium finite overflow force accepted");
    require(api->step(world,&step)==0,"Medium rejects nonfinite projection");
    require(api->download(world,4,nodes)!=0,"Medium rollback download");
    require(std::memcmp(nodes,accepted.data(),sizeof(nodes))==0,"Medium rollback restores accepted nodes");
    step.materialMode=2;require(!api->step(world,&step),"invalid material mode rejected");
}
void testDampedEquilibrium(const btGpuXpbdTgsApi* api,void* world)
{
    // Uniaxial StVK equilibrium with nu=0: F = E V stretch (stretch^2-1)/2.
    // Damping must dissipate motion without changing this static load-extension relation.
    double stretch=1;
    for(int i=0;i<10;++i) stretch-=(10000./12*stretch*(stretch*stretch-1)-1)/(10000./12*(3*stretch*stretch-1));
    for(double beta:{0.,.01,.1})
    {
        btGpuXpbdNode nodes[4]={{{0,0,0},0,{0,0,0},10},{{1,0,0},1,{0,0,0},10},
            {{0,1,0},0,{0,0,0},10},{{0,0,1},0,{0,0,0},10}};
        btGpuXpbdTet tet={{0,1,2,3},{1,0,0,0,1,0,0,0,1},1./6,10000,0,beta};
        btGpuXpbdVec force[4]={{0,0,0},{1,0,0},{0,0,0},{0,0,0}};
        btGpuXpbdStep cfg={.001,1,{0,0,0},0,0,20,0,0,0};
        require(api->upload(world,4,nodes,1,&tet)&&api->forces(world,4,force),"damped equilibrium upload");
        for(int i=0;i<2000;++i) if(!api->step(world,&cfg)) { require(false,api->error(world));break; }
        require(api->download(world,4,nodes)!=0,"damped equilibrium download");
        std::printf("equilibrium beta=%g: actual %.9g expected %.9g extension error %.3f%%\n",beta,nodes[1].position.x,stretch,100*std::fabs(nodes[1].position.x-stretch)/(stretch-1));
        require(std::fabs(nodes[1].position.x-stretch)<(stretch-1)*.002,"damping preserves static StVK stiffness");
    }
}
void testMaterialScale(const btGpuXpbdTgsApi* api,void* world)
{
    // Analytic StVK energy/gradient for F=diag(1.01,1,1), in three unit systems.
    // Metre-sized small elements used to lose all stiffness to an absolute epsilon.
    for(double scale:{.001,1.,1000.}) for(double nu:{0.,.3,.49})
    {
        const double length=.01*scale,young=7.e10/scale,h=2.e-7;
        btGpuXpbdNode nodes[4]={{{0,0,0},4000,{0,0,0},0},{{length*1.01,0,0},4000,{0,0,0},0},
            {{0,length,0},4000,{0,0,0},0},{{0,0,length},4000,{0,0,0},0}};
        auto before=std::vector<btGpuXpbdNode>(nodes,nodes+4);
        const double volume=length*length*length/6;
        btGpuXpbdTet tet={{0,1,2,3},{1/length,0,0,0,1/length,0,0,0,1/length},volume,young,nu,0};
        const double strain=(1.01*1.01-1)*.5,mu=1/(2*(1+nu)),lambda=nu/((1+nu)*(1-2*nu));
        const double c2=volume*(2*mu+lambda)*strain*strain;
        const double gx=volume*1.01*(2*mu+lambda)*strain/length,gy=volume*lambda*strain/length;
        const double multiplier=-c2/(8000*(gx*gx+2*gy*gy)+c2/(young*h*h));
        btGpuXpbdVec gradient[4]={{-gx,-gy,-gy},{gx,0,0},{0,gy,0},{0,0,gy}};
        btGpuXpbdStep cfg={h,1,{0,0,0},0,0,1,0,0,0};
        require(api->upload(world,4,nodes,1,&tet)!=0&&api->step(world,&cfg)!=0&&api->download(world,4,nodes)!=0,"scaled material step");
        for(int j=0;j<4;++j)
        {
            auto expected=before[j].position;const double factor=4000*multiplier;
            expected.x+=factor*gradient[j].x;expected.y+=factor*gradient[j].y;expected.z+=factor*gradient[j].z;
            require(distance(nodes[j].position,expected)<tolerance(length*1.e-13,length),"StVK correction agrees with analytic XPBD at every scale");
        }
    }
    // Regrouping identical advancing increments cannot change constant-force integration.
    btGpuXpbdNode result[2];
    for(int grouping=0;grouping<2;++grouping)
    {
        btGpuXpbdNode node={{0,0,0},2,{0,0,0},0};btGpuXpbdVec force={3,0,0};
        btGpuXpbdStep cfg={.01,grouping?10:1,{0,0,0},0,0,grouping?1:10,0,0,0};
        require(api->upload(world,1,&node,0,nullptr)&&api->forces(world,1,&force)&&api->step(world,&cfg)&&api->download(world,1,&result[grouping]),"force grouping step");
    }
    require(distance(result[0].position,result[1].position)<tolerance(1.e-15),"force is integrated at every advancing increment");
    require(std::fabs(result[0].position.x-.00033)<tolerance(1.e-15)&&std::fabs(result[0].velocity.x-.06)<tolerance(1.e-15),"force displacement and impulse are correct");
}
void testPartitionMomentum(const btGpuXpbdTgsApi* api,void* world)
{
    // Nine elements sharing one node exercise multiple graph colors and unequal valence.
    for(double young:{5.,1.e6,7.e7}) for(double damping:{0.,.01})
    {
        std::vector<btGpuXpbdNode> nodes(28);
        nodes[0]={{0,0,0},2,{1,0,0},0};
        std::vector<btGpuXpbdTet> tets;
        for(int i=0;i<9;++i)
        {
            int a=1+3*i;
            nodes[a]={{1,0,0},3,{0,0,0},0};
            nodes[a+1]={{0,1,0},4,{0,0,0},0};
            nodes[a+2]={{0,0,1},5,{0,0,0},0};
            tets.push_back({{0,a,a+1,a+2},{1,0,0,0,1,0,0,0,1},1./6,young,.4,damping});
        }
        require(api->upload(world,int(nodes.size()),nodes.data(),int(tets.size()),tets.data())!=0,"partition momentum upload");
        btGpuXpbdStep step={.001,2,{0,0,0},0,0,20};
        for(int i=0;i<10;++i) require(api->step(world,&step)!=0,"partition momentum step");
        require(api->download(world,int(nodes.size()),nodes.data())!=0,"partition momentum download");
        btGpuXpbdVec momentum={0,0,0};
        for(const auto& n:nodes)
        {
            momentum.x+=n.velocity.x/n.inverseMass;
            momentum.y+=n.velocity.y/n.inverseMass;
            momentum.z+=n.velocity.z/n.inverseMass;
        }
        std::printf("partition momentum error young=%g damping=%g: %.9g\n",young,damping,distance(momentum,{.5,0,0}));
        require(distance(momentum,{.5,0,0})<2.e-5,"unequal partition copies conserve momentum");
    }
}
void testFriction(const btGpuXpbdTgsApi* api,void* world)
{
    auto run=[&](double mu,double speed,double height,double push,int velocityIterations,bool movingGround)
    {
        btGpuXpbdNode nodes[]={{{0,height,0},1,{speed,0,0},0},
            {{-20,0,-20},movingGround?1.:0.,{0,0,0},0},
            {{0,0,20},movingGround?1.:0.,{0,0,0},0},
            {{20,0,-20},movingGround?1.:0.,{0,0,0},0}};
        btGpuXpbdSurfaceVertex maps[4];
        for(int i=0;i<4;++i)
            maps[i]=i&&!movingGround?btGpuXpbdSurfaceVertex{{-1,-1,-1,-1},{0,0,0,0},nodes[i].position,1}:
                btGpuXpbdSurfaceVertex{{i,-1,-1,-1},{1,0,0,0},{0,0,0},i?1:0};
        btGpuXpbdSurfaceTriangle face={{1,2,3}};
        require(api->upload(world,4,nodes,0,nullptr)!=0,"friction nodes");
        require(api->surface(world,4,maps,1,&face,.1,0)!=0,"friction surface");
        double coefficients[]={mu,1};
        require(api->friction(world,2,coefficients)!=0,"friction coefficients");
        btGpuXpbdVec forces[]={{push,-10,0},{0,0,0},{0,0,0},{0,0,0}};
        require(api->forces(world,4,forces)!=0,"friction load");
        btGpuXpbdStep step={.001,1,{0,0,0},0,0,4,0,0,velocityIterations};
        for(int i=0;i<100;++i) require(api->step(world,&step)!=0,"friction step");
        require(api->download(world,4,nodes)!=0,"friction result");
        if(movingGround)
        {
            btGpuXpbdVec momentum={0,0,0};
            for(auto& n:nodes) { momentum.x+=n.velocity.x;momentum.y+=n.velocity.y;momentum.z+=n.velocity.z; }
            require(distance(momentum,{speed,-1,0})<2.e-4,"friction preserves two-body momentum");
            require(nodes[1].velocity.x>0,"friction transfers tangent momentum to second body");
        }
        return nodes[0];
    };
    auto free=run(0,1,0,0,1,false);
    std::printf("friction free vx=%g x=%g\n",free.velocity.x,free.position.x);
    require(std::fabs(free.velocity.x-1)<1.e-5,"zero friction retains tangent velocity");
    auto sliding=run(.5,1,0,0,1,false);
    std::printf("friction sliding vx=%g x=%g\n",sliding.velocity.x,sliding.position.x);
    require(std::fabs(sliding.velocity.x-.5)<.025,"sliding deceleration matches Coulomb friction");
    auto manyPasses=run(.5,1,0,0,4,false);
    require(std::fabs(manyPasses.velocity.x-sliding.velocity.x)<.01,"velocity iterations share one friction budget");
    auto sticking=run(.5,0,0,2,1,false);
    require(std::fabs(sticking.velocity.x)<.002&&std::fabs(sticking.position.x)<.0002,"sub-limit lateral load sticks");
    double noFriction[]={0,1};
    require(api->friction(world,2,noFriction)!=0,"friction changes without mesh upload");
    btGpuXpbdStep changedStep={.001,1,{0,0,0},0,0,4,0,0,1};
    require(api->step(world,&changedStep)!=0,"step after friction change");
    btGpuXpbdNode changed[4];api->download(world,4,changed);
    require(changed[0].velocity.x>.0015,"disabling friction immediately releases sticking");
    auto airborne=run(.5,1,1,0,1,false);
    require(std::fabs(airborne.velocity.x-1)<1.e-5,"friction does not act across open gaps");
    run(.5,1,0,0,1,true);
    double invalid[]={-1};
    require(!api->friction(world,1,invalid),"negative friction rejected");
    invalid[0]=std::numeric_limits<double>::quiet_NaN();
    require(!api->friction(world,1,invalid),"nonfinite friction rejected");
}
void testBroadGroundContact(const btGpuXpbdTgsApi* api,void* world)
{
    btGpuXpbdNode nodes[]={{{-.05,.02,-.05},1,{0,0,0},0},
        {{0,.02,.05},1,{0,0,0},0},{{.05,.02,-.05},1,{0,0,0},0}};
    btGpuXpbdSurfaceVertex maps[7];
    for(int i=0;i<3;++i) maps[i]={{i,-1,-1,-1},{1,0,0,0},{0,0,0},0};
    const btGpuXpbdVec corners[]={{-5,0,-5},{-5,0,5},{5,0,5},{5,0,-5}};
    for(int i=0;i<4;++i) maps[i+3]={{-1,-1,-1,-1},{0,0,0,0},corners[i],1};
    btGpuXpbdSurfaceTriangle faces[]={{{0,1,2}},{{3,4,5}},{{3,5,6}}};
    require(api->upload(world,3,nodes,0,nullptr)!=0,"broad ground nodes");
    require(api->surface(world,7,maps,3,faces,.1,0)!=0,"broad ground surface");
    btGpuXpbdStep step={.001,16,{0,-9.81,0},0,0,50,0,0,1};
    double minimum=1;
    for(int i=0;i<150;++i)
    {
        require(api->step(world,&step)!=0,"broad ground falling step");
        require(api->download(world,3,nodes)!=0,"broad ground result");
        for(auto& node:nodes) minimum=std::fmin(minimum,node.position.y);
    }
    std::printf("broad ground minimum=%g\n",minimum);
    require(minimum>-.00001,"every incident surface vertex remains above the large ground face");
    for(auto& node:nodes) require(std::fabs(node.position.y)<.00001,"all corners settle on the floor");
}
void testOutsideStaticFace(const btGpuXpbdTgsApi* api,void* world)
{
    btGpuXpbdNode nodes[]={{{-.02,-.01,.2},1,{0,0,0},0},
        {{-.03,-.01,.2},1,{0,0,0},0},{{-.02,-.01,.3},1,{0,0,0},0}};
    btGpuXpbdSurfaceVertex maps[6];
    for(int i=0;i<3;++i) maps[i]={{i,-1,-1,-1},{1,0,0,0},{0,0,0},0};
    const btGpuXpbdVec corners[]={{0,0,0},{0,0,1},{1,0,0}};
    for(int i=0;i<3;++i) maps[i+3]={{-1,-1,-1,-1},{0,0,0,0},corners[i],1};
    btGpuXpbdSurfaceTriangle faces[]={{{0,1,2}},{{3,4,5}}};
    require(api->upload(world,3,nodes,0,nullptr)!=0,"outside face nodes");
    require(api->surface(world,6,maps,2,faces,.1,0)!=0,"outside face surface");
    btGpuXpbdStep step={.001,1,{0,0,0},0,0,4,0,0,1};
    require(api->step(world,&step)!=0,"outside face step");
    api->download(world,3,nodes);
    for(auto& node:nodes) require(distance(node.velocity,{0,0,0})<1.e-5&&std::fabs(node.position.y+.01)<1.e-6,
        "a point beyond a static face's edge is not pushed by its extended plane");
}
void testMouseSprings(const btGpuXpbdTgsApi* api,void* world)
{
    btGpuXpbdNode node={{0,0,0},1,{0,0,0},0};
    require(api->upload(world,1,&node,0,nullptr)!=0,"mouse spring node");
    btGpuXpbdSpring spring={0,{1,0,0},100,10,1000};
    btGpuXpbdStep step={.1,1,{0,0,0},0,0,1};
    require(api->springs(world,1,&spring)!=0,"mouse spring upload");
    require(api->step(world,&step)!=0,"mouse spring step");
    api->download(world,1,&node);
    require(std::fabs(node.velocity.x-10./3)<1.e-5,"implicit spring matches analytic damped step");
    require(std::fabs(node.position.x-1./3)<1.e-5,"spring moves picked node");
    spring.target={-1,0,0};
    require(api->springs(world,1,&spring)!=0,"mouse target changes without graph rebuild");
    require(api->step(world,&step)!=0,"moved target step");
    api->download(world,1,&node);
    require(node.velocity.x<0,"picked node follows moved target");
    require(api->springs(world,0,nullptr)!=0,"mouse release");
    double releaseVelocity=node.velocity.x;
    require(api->step(world,&step)!=0,"released node step");api->download(world,1,&node);
    require(std::fabs(node.velocity.x-releaseVelocity)<1.e-6,"release preserves velocity");
    node={{0,0,0},1,{0,0,0},0};
    api->updateNodes(world,1,&node);
    spring={0,{100,0,0},1.e12,0,.5};
    require(api->springs(world,1,&spring)!=0,"force-limited spring");
    require(api->step(world,&step)!=0,"stiff capped spring remains finite");api->download(world,1,&node);
    require(std::fabs(node.velocity.x-.05)<1.e-6,"elastic mouse force respects per-node cap");
    node.inverseMass=0;api->updateNodes(world,1,&node);
    auto pinned=node.position;
    api->step(world,&step);api->download(world,1,&node);
    require(distance(node.position,pinned)==0,"mouse spring does not move a fixed node");
    spring.node=1;require(!api->springs(world,1,&spring),"invalid picked node rejected");
    spring.node=0;spring.damping=-1;require(!api->springs(world,1,&spring),"invalid picking damping rejected");
}
}
// Collision geometry can overlap another element without changing its material ownership.
static void testSoftContactMaterial(const btGpuXpbdTgsApi* api,void* world)
{
    for(int mode : {0,1})
    {
        btGpuXpbdNode nodes[12]={};
        const btGpuXpbdVec rest[4]={{0,0,0},{1,0,0},{0,1,0},{0,0,1}};
        btGpuXpbdTet tets[3]={};
        for(int body=0;body<3;++body)
        {
            tets[body]={{body*4,body*4+1,body*4+2,body*4+3},{1,0,0,0,1,0,0,0,1},1./6,10000,.4,0};
            for(int j=0;j<4;++j)
            {
                auto p=rest[j];
                if(body==1) p.y+=2;
                if(body==2) {p.x+=.2;p.y+=2.00005;p.z+=.2;}
                nodes[body*4+j]={p,1,{0,body==2?-1.:0.,0},0};
            }
        }
        btGpuXpbdSurfaceVertex maps[5]={
            {{0,1,2,3},{1,0,0,0},{0,2,0},0},
            {{0,1,2,3},{0,0,0,1},{0,2,0},0},
            {{0,1,2,3},{0,1,0,0},{0,2,0},0},
            {{4,5,6,7},{1,0,0,0},{10,0,0},0},
            {{8,9,10,11},{1,0,0,0},{0,0,0},1}};
        btGpuXpbdSurfaceTriangle face={{0,1,2}};
        btGpuXpbdStep step={.0001,1,{0,0,0},0,0,1,0,0,1};step.materialMode=mode;
        require(api->upload(world,12,nodes,3,tets)&&api->surface(world,5,maps,1,&face,.01,0)&&
            api->step(world,&step)&&api->download(world,12,nodes),"soft contact material setup and solve");
        double supportedMomentum=0,totalMomentum=0;
        for(int i=0;i<12;++i) {totalMomentum+=nodes[i].velocity.y;if(i<4)supportedMomentum+=nodes[i].velocity.y;}
        require(supportedMomentum<-.01,"soft contact impulse reaches the surface embedding");
        for(int i=4;i<8;++i) require(distance(nodes[i].velocity,{0,0,0})<tolerance(1.e-8),"nearby tetrahedron does not steal soft contact");
        require(std::fabs(totalMomentum+4)<1.e-6,"embedded soft contact conserves momentum");
    }
}

void testInversionRecovery(const btGpuXpbdTgsApi* api,void* world)
{
    for(double scale : { .001,1.,1000. })
    {
        btGpuXpbdTet tet={{0,1,2,3},{1/scale,0,0,0,1/scale,0,0,0,1/scale},scale*scale*scale/6,50,.4,0};
        btGpuXpbdStep step={.01*scale,1,{0,0,0},0,0,1,0,0,1};
        for(double height : {-1.,0.,1.})
        {
            btGpuXpbdNode nodes[]={{{0,0,0},0},{{scale,0,0},0},{{0,scale,0},0},{{0,0,height*scale},1/(scale*scale*scale),{0,0,0},4/scale}};
            require(api->upload(world,4,nodes,1,&tet)!=0,"inversion recovery upload");
            btGpuXpbdNode uploaded[4];require(api->download(world,4,uploaded)!=0,"fixed node upload baseline");
            for(int k=0;k<500;++k) require(api->step(world,&step)!=0,"inversion recovery step");
            require(api->download(world,4,nodes)!=0,"inversion recovery download");
            require(std::fabs(nodes[3].position.z/scale-1)<.01,"inverted and flat tetrahedra recover rest shape");
            // FP32 accuracy: require exact immobility of the stored state, independently of upload rounding.
            for(int j=0;j<3;++j) require(distance(nodes[j].position,uploaded[j].position)==0,"recovery preserves fixed nodes");
            if(height==1) require(distance(nodes[3].position,{0,0,scale})<tolerance(scale*1.e-10,scale),"recovery leaves rest equilibrium unchanged");
        }
        btGpuXpbdNode collapsed[4]={};
        for(auto& node:collapsed) node.inverseMass=1;
        require(api->upload(world,4,collapsed,1,&tet)&&api->step(world,&step)&&api->download(world,4,collapsed),"fully collapsed element remains finite");
    }
}

int main()
{
#ifdef XPBD_TEST_RUNTIME_LOADER
    btGpuXpbdTgsBackend backend;
    if(!backend.ready()) { std::fprintf(stderr,"%s\n",backend.error().c_str());return 1; }
    const auto* api=backend.api();void* world=backend.state();
#else
    require(btGetGpuXpbdTgsApi(5)==nullptr,"reject previous float ABI");
    const auto* api=btGetGpuXpbdTgsApi(BT_GPU_XPBD_ABI_VERSION);void* world=api->create();
#endif
    if(!world) { std::fprintf(stderr,"%s\n",api->error(nullptr));return 1; }
    testSoftContactMaterial(api,world);
    testOutsideStaticFace(api,world);
    testBroadGroundContact(api,world);
    testFriction(api,world);
    testMouseSprings(api,world);
    testPartitionMomentum(api,world);
    testMaterialModeSwitch(api,world);
    testContactReports(api,world);
    testCompensatedState(api,world);
    testCallbackContactReuse(api,world);
    testContactImpulses(api,world);
    testRigidCoupling(api,world);
    testRigidContacts(api,world);
    testRigidJoints(api,world);
    testSelfCollision(api,world);
    testInversionRecovery(api,world);
    auto mediumApi=*api;accurateStep=api->step;mediumApi.step=mediumStep;
    testContactImpulses(&mediumApi,world);
    testPartitionMomentum(&mediumApi,world);
    testOutsideStaticFace(&mediumApi,world);
    testBroadGroundContact(&mediumApi,world);
    testFriction(&mediumApi,world);
    testMouseSprings(&mediumApi,world);

    btGpuXpbdTet tet={{0,1,2,3},{1,0,0,0,1,0,0,0,1},1./6,10000,.4,.01};
    testMaterialScale(api,world);
    testDampedEquilibrium(api,world);
    std::vector<btGpuXpbdNode> rest={{{0,0,0},4,{0,0,0},0},{{1,0,0},4,{0,0,0},0},{{0,1,0},4,{0,0,0},0},{{0,0,1},4,{0,0,0},0}};
    btGpuXpbdStep cfg={.001,10,{0,0,0},0,0,1};
    auto run=[&](std::vector<btGpuXpbdNode> input,int steps)
    {
        require(api->upload(world,4,input.data(),1,&tet)!=0,"upload");
        for(int k=0;k<steps;++k) if(!api->step(world,&cfg)) { require(false,api->error(world));break; }
        require(api->download(world,4,input.data())!=0,"download");return input;
    };
    auto state=run(rest,100);
    for(int i=0;i<4;++i) require(distance(state[i].position,rest[i].position)<1.e-6,"rest equilibrium");
    btGpuXpbdVec samples[2]={{.2,.3,.1},{1.2,.1,.1}};
    btGpuXpbdSurfaceVertex embedded[2];
    require(api->embed(world,2,samples,0,4,7,embedded)!=0,"GPU embedding");
    for(int i=0;i<2;++i)
    {
        btGpuXpbdVec p={0,0,0};double sum=0;
        for(int j=0;j<4;++j) { auto v=rest[embedded[i].nodes[j]].position;double w=embedded[i].weights[j];sum+=w;p.x+=v.x*w;p.y+=v.y*w;p.z+=v.z*w; }
        require(distance(p,samples[i])<1.e-6&&std::fabs(sum-1)<1.e-6,"embedding reproduces interior and exterior points");
    }
    auto translated=rest;for(auto& n:translated) n.velocity={1,2,-3};
    state=run(translated,100);
    for(int i=0;i<4;++i)
    {
        auto p=rest[i].position;p.x+=.1;p.y+=.2;p.z-=.3;
        require(distance(state[i].position,p)<1.e-4,"uniform translation");
        require(distance(state[i].velocity,{1,2,-3})<1.e-3,"uniform velocity");
    }
    auto rotated=rest;for(auto& n:rotated) { auto p=n.position;n.position={p.x*std::cos(.7)-p.y*std::sin(.7),p.x*std::sin(.7)+p.y*std::cos(.7),p.z}; }
    state=run(rotated,100);
    for(int i=0;i<4;++i) require(distance(state[i].position,rotated[i].position)<1.e-4,"rigid rotation equilibrium");
    auto strained=rest;strained[1].position.x=1.2;
    state=run(strained,100);
    btGpuXpbdVec momentum={0,0,0},com={0,0,0};
    for(auto& n:state) { momentum.x+=n.velocity.x/4;momentum.y+=n.velocity.y/4;momentum.z+=n.velocity.z/4;com.x+=n.position.x/4;com.y+=n.position.y/4;com.z+=n.position.z/4; }
    require(distance(momentum,{0,0,0})<1.e-4,"internal forces preserve linear momentum");
    require(distance(com,{.3,.25,.25})<1.e-4,"internal forces preserve COM");
    require(std::fabs(distance(state[0].position,state[1].position)-1)<.05,"elastic recovery");
    auto pinned=strained;pinned[0].inverseMass=0;state=run(pinned,100);
    require(distance(state[0].position,pinned[0].position)==0,"pinned node");
    auto local=rest;local[1].position.x=1.25;
    auto savedCfg=cfg;cfg.dt=.0001;cfg.substeps=1;
    auto nearOrigin=run(local,1);
    for(auto& n:local) { n.position.x+=4096;n.position.y+=4096;n.position.z+=4096; }
    auto farOrigin=run(local,1);
    for(int i=0;i<4;++i)
    {
        double error=distance(nearOrigin[i].velocity,farOrigin[i].velocity);
        if(error>1.e-5) std::fprintf(stderr,"translation material velocity error node %d: %.9g\n",i,error);
        require(error<1.e-5,"material corrections independent of world translation");
    }
    cfg=savedCfg;
    std::vector<btGpuXpbdNode> contactNodes={{{0,.01,0},1,{0,-1,0},0},{{-2,0,-2},1,{0,0,0},0},{{0,0,2},1,{0,0,0},0},{{2,0,-2},1,{0,0,0},0}};
    cfg.velocityIterations=1;
    std::vector<btGpuXpbdSurfaceVertex> maps(4);
    for(int i=0;i<4;++i) maps[i]={{i,-1,-1,-1},{1,0,0,0},{0,0,0},i==0?0:1};
    btGpuXpbdSurfaceTriangle face={{1,2,3}};
    require(api->upload(world,4,contactNodes.data(),0,nullptr)!=0,"contact nodes");
    require(api->surface(world,4,maps.data(),1,&face,.1,0)!=0,"dynamic surface upload");
    for(int k=0;k<100;++k) require(api->step(world,&cfg)!=0,"dynamic contact step");
    state=contactNodes;require(api->download(world,4,state.data())!=0,"contact download");
    momentum={0,0,0};for(auto& n:state) { momentum.x+=n.velocity.x;momentum.y+=n.velocity.y;momentum.z+=n.velocity.z; }
    require(distance(momentum,{0,-1,0})<1.e-4,"dynamic contact transfers equal and opposite momentum");
    require(state[0].velocity.y>-.99 && state[1].velocity.y<0,"both bodies respond to contact");
    for(int i=1;i<4;++i) maps[i]={{-1,-1,-1,-1},{0,0,0,0},contactNodes[i].position,1};
    require(api->upload(world,4,contactNodes.data(),0,nullptr)!=0,"reset contact nodes");
    require(api->surface(world,4,maps.data(),1,&face,.1,0)!=0,"static surface upload");
    for(int k=0;k<100;++k) require(api->step(world,&cfg)!=0,"static contact step");
    require(api->download(world,4,state.data())!=0,"static contact download");
    require(state[0].position.y>=-.001,"static triangle prevents penetration");
    btGpuXpbdVec surfacePositions[4];require(api->surfacePositions(world,4,surfacePositions)!=0,"GPU skin download");
    require(distance(surfacePositions[0],state[0].position)<1.e-7,"dynamic GPU skinning");
    require(distance(surfacePositions[1],contactNodes[1].position)<1.e-7,"static GPU skinning");
    savedCfg=cfg;cfg.dt=.0001;cfg.substeps=1;cfg.iterations=1;
    contactNodes[0].position={0,.00015,0};contactNodes[0].velocity={0,-1,0};
    require(api->upload(world,4,contactNodes.data(),0,nullptr)!=0,"look-ahead contact nodes");
    require(api->surface(world,4,maps.data(),1,&face,.1,0)!=0,"look-ahead contact surface");
    require(api->step(world,&cfg)!=0,"look-ahead contact step");api->download(world,4,state.data());
    require(std::fabs(state[0].position.y-.00005)<1.e-7,"velocity contact projection preserves position");
    require(std::fabs(state[0].velocity.y+.5)<1.e-4,"velocity pass limits predicted penetration across an open gap");
    contactNodes[0].position={0,.0004,0};
    require(api->upload(world,4,contactNodes.data(),0,nullptr)!=0,"open gap nodes");
    require(api->surface(world,4,maps.data(),1,&face,.1,0)!=0,"open gap surface");
    require(api->step(world,&cfg)!=0,"open gap step");api->download(world,4,state.data());
    require(std::fabs(state[0].velocity.y+1)<1.e-5,"velocity pass leaves a non-closing gap unchanged");
    contactNodes[0].position={0,.00005,0};contactNodes[0].velocity={0,1,0};
    require(api->upload(world,4,contactNodes.data(),0,nullptr)!=0,"separating nodes");
    require(api->surface(world,4,maps.data(),1,&face,.1,0)!=0,"separating surface");
    require(api->step(world,&cfg)!=0,"separating step");api->download(world,4,state.data());
    require(std::fabs(state[0].velocity.y-1)<1.e-5,"velocity pass does not attract separating bodies");
    cfg=savedCfg;
    contactNodes[0].position={2.01,-.001,-2};contactNodes[0].velocity={0,0,0};
    require(api->upload(world,4,contactNodes.data(),0,nullptr)!=0,"edge contact nodes");
    require(api->surface(world,4,maps.data(),1,&face,.1,0)!=0,"edge contact surface");
    require(api->step(world,&cfg)!=0,"edge contact step");api->download(world,4,state.data());
    require(std::fabs(state[0].velocity.x)<1.e-5&&std::fabs(state[0].velocity.z)<1.e-5,"frictionless face contact does not pull sideways behind an edge");
    require(api->upload(world,4,rest.data(),0,nullptr)!=0,"drive nodes");
    btGpuXpbdTarget target={0,{1,2,3}};
    require(api->targets(world,1,&target)!=0,"drive target");
    cfg.iterations=3;
    require(api->step(world,&cfg)!=0,"drive step");
    require(api->download(world,4,state.data())!=0,"drive download");
    require(distance(state[0].position,target.position)<1.e-6&&state[0].inverseMass==0,"drive reaches target");
    target.position={2,3,4};require(api->targets(world,1,&target)!=0,"update target");
    require(api->step(world,&cfg)!=0,"replayed drive graph");api->download(world,4,state.data());
    require(distance(state[0].position,target.position)<1.e-6,"graph uses updated target buffer");
    require(api->targets(world,0,nullptr)!=0,"release target");api->download(world,4,state.data());
    require(state[0].inverseMass==rest[0].inverseMass,"release restores mass");
    std::vector<btGpuXpbdNode> crossed={{{-2,.00015,-.1},1,{0,-1,0},0},{{2,.00015,-.1},1,{0,-1,0},0},{{2,.00015,.1},1,{0,-1,0},0},
        {{-.1,0,-2},1,{0,0,0},0},{{0,0,2},1,{0,0,0},0},{{.1,0,-2},1,{0,0,0},0}};
    btGpuXpbdSurfaceTriangle crossingFaces[2]={{{0,1,2}},{{3,4,5}}};
    std::vector<btGpuXpbdSurfaceVertex> crossingMaps(6);
    for(int i=0;i<6;++i) crossingMaps[i]={{i,-1,-1,-1},{1,0,0,0},{0,0,0},i<3?0:1};
    cfg={.0001,1,{0,0,0},0,0,1,0,0,1};
    require(api->upload(world,6,crossed.data(),0,nullptr)!=0,"crossing edge nodes");
    require(api->surface(world,6,crossingMaps.data(),2,crossingFaces,.01,0)!=0,"crossing edge surfaces");
    require(api->step(world,&cfg)!=0,"crossing edge step");
    auto crossedResult=crossed;api->download(world,6,crossedResult.data());
    double sumVelocity=0,otherVelocity=0;
    for(int i=0;i<6;++i) { sumVelocity+=crossedResult[i].velocity.y;if(i>=3) otherVelocity+=crossedResult[i].velocity.y; }
    require(otherVelocity<-.01,"crossing edges transfer momentum without vertex-face proximity");
    require(std::fabs(sumVelocity+3)<1.e-4,"triangle contacts preserve equal and opposite momentum");
    // Crossed volume faces must collide even when no vertex is near the other face.
    for(int materialMode : {0,1})
    {
        cfg.materialMode=materialMode;
        auto edgeVolumes=crossed;
        edgeVolumes.push_back({{-2,1.00015,-.1},1,{0,-1,0},0});
        edgeVolumes.push_back({{-.1,-1,-2},1,{0,0,0},0});
        btGpuXpbdTet edgeTets[2]={
            {{0,1,2,6},{.25,0,0,0,0,1,-5,5,0},.8/6,10000,.4,0},
            {{3,4,5,7},{0,5,0,0,0,-1,.25,-.125,0},.8/6,10000,.4,0}};
        require(api->upload(world,8,edgeVolumes.data(),2,edgeTets)&&api->surface(world,6,crossingMaps.data(),2,crossingFaces,.01,0),"crossing volume edge setup");
        require(api->step(world,&cfg)&&api->download(world,8,edgeVolumes.data()),"crossing volume edge step");
        double edgeMomentum=0,lowerEdgeMomentum=0;
        for(int i=0;i<8;++i) {edgeMomentum+=edgeVolumes[i].velocity.y;if((i>=3&&i<6)||i==7)lowerEdgeMomentum+=edgeVolumes[i].velocity.y;}
        require(lowerEdgeMomentum<-.01,"crossing volume edges transfer momentum without nearby vertices");
        require(std::fabs(edgeMomentum+4)<1.e-4,"crossing volume edges conserve momentum");
    }
    cfg.materialMode=0;
    for(int i=3;i<6;++i) crossingMaps[i]={{-1,-1,-1,-1},{0,0,0,0},crossed[i].position,1};
    require(api->upload(world,6,crossed.data(),0,nullptr)!=0,"static crossing edge nodes");
    require(api->surface(world,6,crossingMaps.data(),2,crossingFaces,.01,0)!=0,"static crossing edge surfaces");
    require(api->step(world,&cfg)!=0,"static crossing edge step");api->download(world,6,crossedResult.data());
    sumVelocity=0;for(int i=0;i<3;++i) sumVelocity+=crossedResult[i].velocity.y;
    require(sumVelocity>-2.99,"static triangle catches crossing edges");
    auto crossingStep=crossed;
    for(int i=0;i<3;++i) crossingStep[i].position.y=.00005;
    require(api->upload(world,6,crossingStep.data(),0,nullptr)!=0,"first increment crossing nodes");
    require(api->surface(world,6,crossingMaps.data(),2,crossingFaces,.01,0)!=0,"first increment crossing plane");
    require(api->step(world,&cfg)!=0,"first increment crossing step");api->download(world,6,crossedResult.data());
    sumVelocity=0;for(int i=0;i<3;++i) sumVelocity+=crossedResult[i].velocity.y;
    require(sumVelocity>-2.99,"contact prepared before integration survives crossing static face");
    crossed[0].position.y=.0033;crossed[1].position.y=crossed[2].position.y=-.003;
    require(api->upload(world,6,crossed.data(),0,nullptr)!=0,"backface triangle nodes");
    require(api->surface(world,6,crossingMaps.data(),2,crossingFaces,.01,0)!=0,"backface triangle surfaces");
    require(api->step(world,&cfg)!=0,"backface triangle step");api->download(world,6,crossedResult.data());
    for(int i=0;i<3;++i) require(std::fabs(crossedResult[i].velocity.y+1)<1.e-5,"static mesh culls triangle centroid behind face");
    auto remapNodes=rest;for(auto& n:remapNodes) n.velocity={0,-1,0};
    btGpuXpbdSurfaceVertex remapMaps[4]={{{0,1,2,3},{-.2,1.2,0,0},{0,.00015,0},0},
        {{-1,-1,-1,-1},{0,0,0,0},{-4,0,-4},1},{{-1,-1,-1,-1},{0,0,0,0},{0,0,4},1},{{-1,-1,-1,-1},{0,0,0,0},{4,0,-4},1}};
    require(api->upload(world,4,remapNodes.data(),1,&tet)!=0,"contact remap tetrahedron");
    require(api->surface(world,4,remapMaps,1,&face,.01,0)!=0,"contact remap surface");
    require(api->step(world,&cfg)!=0,"contact remap step");api->download(world,4,remapNodes.data());
    require(std::fabs(remapNodes[1].velocity.y+.5)<1.e-3,"outside contact maps to nearest tetrahedron point");
    require(std::fabs(remapNodes[0].velocity.y+1)<1.e-3,"contact remap does not extrapolate negative force weights");
    // The entire collision triangle lies outside the tet AABB. Cooking must
    // retain its embedding tet, rather than falling back to negative weights.
    btGpuXpbdSurfaceVertex outsideMaps[6]={remapMaps[0],remapMaps[0],remapMaps[0],remapMaps[1],remapMaps[2],remapMaps[3]};
    outsideMaps[0].offset.z=-.001;outsideMaps[1].offset.x=.001;outsideMaps[2].offset.z=.001;
    remapNodes=rest;for(auto& n:remapNodes) n.velocity={0,-1,0};
    require(api->upload(world,4,remapNodes.data(),1,&tet)!=0,"outside triangle tet");
    require(api->surface(world,6,outsideMaps,2,crossingFaces,.01,0)!=0,"outside triangle candidates");
    require(api->step(world,&cfg)!=0,"outside triangle contact step");api->download(world,4,remapNodes.data());
    require(std::fabs(remapNodes[0].velocity.y+1)<1.e-3,"outside triangle retains clamped embedding tet");
    require(remapNodes[1].velocity.y>-.99,"outside triangle transfers contact through cooked candidate");
    // A nearby part of the same body must not steal an embedded surface's contact.
    auto ownershipNodes=rest;ownershipNodes.insert(ownershipNodes.end(),rest.begin(),rest.end());
    btGpuXpbdTet ownershipTets[2]={tet,tet};for(int& n:ownershipTets[1].nodes) n+=4;
    for(int i=0;i<4;++i)
    {
        ownershipNodes[i].position.y+=.00005;ownershipNodes[i].velocity={0,-1,0};
        ownershipNodes[i+4].position.y-=.1;
    }
    btGpuXpbdSurfaceVertex ownershipMaps[7]={};
    const int ownerIds[3]={0,1,3};
    for(int i=0;i<3;++i)
    {
        ownershipMaps[i]={{0,1,2,3},{0,0,0,0},{0,0,0},0};ownershipMaps[i].weights[ownerIds[i]]=1;
        ownershipMaps[i+3]=remapMaps[i+1];
    }
    ownershipMaps[6]={{4,5,6,7},{1,0,0,0},{0,0,0},0};
    require(api->upload(world,8,ownershipNodes.data(),2,ownershipTets)!=0,"contact material ownership upload");
    require(api->surface(world,7,ownershipMaps,2,crossingFaces,.01,0)!=0,"contact material ownership surface");
    require(api->step(world,&cfg)!=0,"contact material ownership step");
    require(api->download(world,8,ownershipNodes.data())!=0,"contact material ownership download");
    double ownerY=0;for(int i=0;i<4;++i) ownerY+=ownershipNodes[i].velocity.y;
    require(ownerY>-3.99,"contact acts on the embedded surface material");
    for(int i=4;i<8;++i) require(distance(ownershipNodes[i].velocity,{0,0,0})<1.e-5,"nearby material does not steal contact");
    contactNodes[0].position={0,-.001,0};contactNodes[0].velocity={0,0,0};cfg={.002,1,{0,0,0},0,0,2,0,0,1};
    require(api->upload(world,4,contactNodes.data(),0,nullptr)!=0,"penetration recovery nodes");
    require(api->surface(world,4,maps.data(),1,&face,.1,0)!=0,"penetration recovery plane");
    require(api->step(world,&cfg)!=0,"penetration recovery step");api->download(world,4,state.data());
    require(std::fabs(state[0].velocity.y)<1.e-5,"static contact releases excess accumulated impulse");
    require(std::fabs(state[0].position.y)<1.e-6,"penetration recovery does not launch a resting node");
    // Exercise the EPA/degenerate-simplex paths, not only separated GJK pairs.
    for(int mode=0;mode<3;++mode)
    {
        auto pair=crossed;
        for(int i=0;i<6;++i)
        {
            pair[i].velocity={0,0,0};pair[i].position.y=0;
            crossingMaps[i]={{i,-1,-1,-1},{1,0,0,0},{0,0,0},i<3?0:1};
        }
        if(mode==1) { pair[0].position.y=-.01;pair[1].position.y=pair[2].position.y=.01; }
        if(mode==2) { pair[1].position=pair[0].position; } // collapsed face is ignored
        cfg={.0001,1,{0,0,0},0,0,2,0,0,1};
        require(api->upload(world,6,pair.data(),0,nullptr)!=0,"EPA test upload");
        require(api->surface(world,6,crossingMaps.data(),2,crossingFaces,.01,0)!=0,"EPA test surfaces");
        require(api->step(world,&cfg)!=0,"EPA and coplanar contact step");
        require(api->download(world,6,crossedResult.data())!=0,"EPA test download");
        btGpuXpbdVec momentum={0,0,0};double speed=0;
        for(const auto& node:crossedResult)
        {
            require(std::isfinite(node.position.x)&&std::isfinite(node.position.y)&&std::isfinite(node.position.z)&&
                std::isfinite(node.velocity.x)&&std::isfinite(node.velocity.y)&&std::isfinite(node.velocity.z),"EPA finite state");
            momentum.x+=node.velocity.x;momentum.y+=node.velocity.y;momentum.z+=node.velocity.z;
            speed+=distance(node.velocity,{0,0,0});
        }
        require(distance(momentum,{0,0,0})<1.e-3,"EPA contacts preserve total momentum");
        if(mode!=1) require(speed<1.e-3,"coplanar or collapsed contact does not launch resting nodes");
    }
    std::vector<btGpuXpbdNode> tinySupport(2,{{0,.00015,0},1,{0,-1,0},0});
    btGpuXpbdSurfaceVertex tinyMaps[4]={{{0,1,-1,-1},{.0005,.9995,0,0},{0,0,0},0},remapMaps[1],remapMaps[2],remapMaps[3]};
    cfg={.0001,1,{0,0,0},0,0,1,0,0,1};
    require(api->upload(world,2,tinySupport.data(),0,nullptr)!=0,"tiny support nodes");
    require(api->surface(world,4,tinyMaps,1,&face,.01,0)!=0,"tiny support surface");
    require(api->step(world,&cfg)!=0,"tiny support contact step");api->download(world,2,tinySupport.data());
    require(std::fabs(tinySupport[0].velocity.y+1)<1.e-6,"tiny weight does not create a contact reference");
    require(std::fabs(tinySupport[1].velocity.y-(-1.+.5/.9995))<1.e-4,"contact denominator uses counted support");
    auto preciseNodes=rest;
    for(auto& n:preciseNodes) { n.position.y+=1024.0001220703125;n.velocity={0,-1,0}; }
    btGpuXpbdSurfaceVertex preciseMaps[4]={{{0,1,2,3},{1,0,0,0},{0,0,0},0},remapMaps[1],remapMaps[2],remapMaps[3]};
    for(int i=1;i<4;++i) preciseMaps[i].offset.y=1024;
    cfg={.0002,1,{0,0,0},0,0,200,0,0,1};
    require(api->upload(world,4,preciseNodes.data(),1,&tet)!=0,"sub-ULP contact nodes");
    require(api->surface(world,4,preciseMaps,1,&face,.01,0)!=0,"sub-ULP contact plane");
    require(api->step(world,&cfg)!=0,"sub-ULP accumulated contact step");api->download(world,4,preciseNodes.data());
    require(preciseNodes[0].velocity.y>-.5,"contact resolves motion smaller than world-position precision");
    auto volumes=rest;volumes.insert(volumes.end(),rest.begin(),rest.end());
    for(int i=0;i<4;++i) { volumes[i].position.y+=.00015;volumes[i].velocity={0,-1,0};volumes[i+4].position.y=-volumes[i+4].position.y; }
    btGpuXpbdTet volumeTets[2]={tet,tet};
    for(int& n:volumeTets[1].nodes) n+=4;volumeTets[1].inverseRest[4]=-1;
    btGpuXpbdSurfaceVertex volumeMaps[6];int volumeIds[6]={0,1,3,4,7,5};
    for(int i=0;i<6;++i)
    {
        int body=i<3?0:1;volumeMaps[i]={{body*4,body*4+1,body*4+2,body*4+3},{0,0,0,0},{0,0,0},body};
        volumeMaps[i].weights[volumeIds[i]-body*4]=1;
    }
    cfg={.0001,1,{0,0,0},0,0,1,0,0,1};
    require(api->upload(world,8,volumes.data(),2,volumeTets)!=0,"volume vertex contact upload");
    require(api->surface(world,6,volumeMaps,2,crossingFaces,.01,0)!=0,"volume vertex contact surfaces");
    require(api->step(world,&cfg)!=0,"volume vertex contact step");api->download(world,8,volumes.data());
    double totalY=0,lowerY=0;
    for(int i=0;i<8;++i) { totalY+=volumes[i].velocity.y;if(i>=4) lowerY+=volumes[i].velocity.y; }
    require(lowerY<-.01,"signed deformable vertex contacts transfer momentum");
    require(std::fabs(totalY+4)<1.e-3,"signed deformable vertex contacts conserve momentum");
    // Separate collision vertices retain their tetrahedral support for self contact too.
    auto selfVolumes=rest;selfVolumes.insert(selfVolumes.end(),rest.begin(),rest.end());
    for(int i=0;i<4;++i) {selfVolumes[i].position.y+=.00015;selfVolumes[i].velocity={0,-1,0};selfVolumes[i+4].position.y=-selfVolumes[i+4].position.y;}
    btGpuXpbdSurfaceVertex selfMaps[6];std::memcpy(selfMaps,volumeMaps,sizeof(selfMaps));
    for(auto& map:selfMaps)map.body=0;int selfId=0;
    require(api->upload(world,8,selfVolumes.data(),2,volumeTets)&&api->surface(world,6,selfMaps,2,crossingFaces,.01,0)&&api->selfCollision(world,1,&selfId),"mapped volumetric self contact setup");
    require(api->step(world,&cfg)&&api->download(world,8,selfVolumes.data()),"mapped volumetric self contact solve");
    totalY=0;lowerY=0;for(int i=0;i<8;++i){totalY+=selfVolumes[i].velocity.y;if(i>=4)lowerY+=selfVolumes[i].velocity.y;}
    require(lowerY<-.01&&std::fabs(totalY+4)<1.e-3,"mapped volumetric self contact transfers and conserves momentum");
    require(api->selfCollision(world,0,nullptr)&&api->step(world,&cfg),"disable self collision on a captured graph");
    // A collision surface outside its simulation tetrahedron has signed support.
    // Dropping its negative contributions would create net contact momentum.
    auto extrapolatedVolumes=rest;extrapolatedVolumes.insert(extrapolatedVolumes.end(),rest.begin(),rest.end());
    btGpuXpbdSurfaceVertex extrapolatedMaps[6];
    std::memcpy(extrapolatedMaps,volumeMaps,sizeof(volumeMaps));
    for(int i=0;i<4;++i)
    {
        extrapolatedVolumes[i].position.y+=.04995;
        extrapolatedVolumes[i].velocity={0,-1,0};
        extrapolatedVolumes[i+4].position.y=-extrapolatedVolumes[i+4].position.y;
    }
    for(int i=0;i<3;++i) { extrapolatedMaps[i].weights[0]+=.05;extrapolatedMaps[i].weights[2]-=.05; }
    require(api->upload(world,8,extrapolatedVolumes.data(),2,volumeTets)!=0,"extrapolated volume contact upload");
    require(api->surface(world,6,extrapolatedMaps,2,crossingFaces,.01,0)!=0,"extrapolated volume contact surface");
    require(api->step(world,&cfg)!=0,"extrapolated volume contact step");
    require(api->download(world,8,extrapolatedVolumes.data())!=0,"extrapolated volume contact download");
    totalY=0;lowerY=0;
    for(int i=0;i<8;++i) { totalY+=extrapolatedVolumes[i].velocity.y;if(i>=4) lowerY+=extrapolatedVolumes[i].velocity.y; }
    require(lowerY<-.01,"extrapolated volume contacts transfer momentum");
    require(std::fabs(totalY+4)<1.e-3,"extrapolated volume contacts conserve momentum");
    auto checkFilter=[&](const std::vector<btGpuXpbdNode>& input,int tetCount,const btGpuXpbdTet* tets,
        int vertexCount,const btGpuXpbdSurfaceVertex* vertices,int triangleCount,const btGpuXpbdSurfaceTriangle* triangles)
    {
        const int count=int(input.size());auto free=input,filtered=input;
        cfg={.0001,1,{0,0,0},0,0,1,0,0,1};
        require(api->upload(world,count,input.data(),tetCount,tets)!=0,"pair filter free mesh");
        require(api->step(world,&cfg)!=0,"pair filter free step");
        require(api->download(world,count,free.data())!=0,"pair filter free state");
        require(api->upload(world,count,input.data(),tetCount,tets)!=0,"pair filter mesh");
        require(api->surface(world,vertexCount,vertices,triangleCount,triangles,.01,0)!=0,"pair filter surface");
        btGpuXpbdDisabledPair disabled[]={{1,0},{0,1},{4,3},{2,1}};
        require(api->disabledPairs(world,4,disabled)!=0,"pair filter normalizes unordered duplicates");
        require(api->step(world,&cfg)!=0,"pair filter disabled step");
        require(api->download(world,count,filtered.data())!=0,"pair filter disabled state");
        for(int i=0;i<count;++i)
            require(distance(free[i].position,filtered[i].position)<1.e-7 && distance(free[i].velocity,filtered[i].velocity)<1.e-6,
                "disabled GPU contacts match free motion");
        btGpuXpbdDisabledPair invalid={-1,0};
        require(api->disabledPairs(world,1,&invalid)==0,"pair filter rejects negative IDs");
        require(api->disabledPairs(world,0,nullptr)!=0,"pair filter reenable");
        require(api->step(world,&cfg)!=0,"pair filter reenabled step");
        require(api->download(world,count,filtered.data())!=0,"pair filter reenabled state");
        double change=0;for(int i=0;i<count;++i) change+=distance(free[i].velocity,filtered[i].velocity);
        require(change>.01,"reenabling pair invalidates graph and restores contacts");
    };
    contactNodes[0].position={0,.00025,0};contactNodes[0].velocity={0,-1,0};
    checkFilter(contactNodes,0,nullptr,4,maps.data(),1,&face);
    crossed[0].position.y=crossed[1].position.y=crossed[2].position.y=.00025;
    for(int i=0;i<6;++i)
    {
        crossed[i].velocity=i<3?btGpuXpbdVec{0,-1,0}:btGpuXpbdVec{0,0,0};
        crossingMaps[i]={{i,-1,-1,-1},{1,0,0,0},{0,0,0},i<3?0:1};
    }
    checkFilter(crossed,0,nullptr,6,crossingMaps.data(),2,crossingFaces);
    for(int i=3;i<6;++i) crossingMaps[i]={{-1,-1,-1,-1},{0,0,0,0},crossed[i].position,1};
    checkFilter(crossed,0,nullptr,6,crossingMaps.data(),2,crossingFaces);
    volumes=rest;volumes.insert(volumes.end(),rest.begin(),rest.end());
    for(int i=0;i<4;++i) { volumes[i].position.y+=.00025;volumes[i].velocity={0,-1,0};volumes[i+4].position.y=-volumes[i+4].position.y; }
    checkFilter(volumes,2,volumeTets,6,volumeMaps,2,crossingFaces);
    auto movingPlaneNodes=contactNodes;
    movingPlaneNodes[0].position={0,.00015,0};movingPlaneNodes[0].velocity={0,0,0};
    btGpuXpbdSurfaceVertex movingPlaneMaps[4];btGpuXpbdTarget planeTargets[3];
    for(int i=0;i<4;++i)
    {
        movingPlaneMaps[i]={{i,-1,-1,-1},{1,0,0,0},{0,0,0},i==0?0:1};
        if(i) { planeTargets[i-1]={i,movingPlaneNodes[i].position};planeTargets[i-1].position.y=.0002; }
    }
    cfg={.0002,1,{0,0,0},0,0,2,0,0,1};
    require(api->upload(world,4,movingPlaneNodes.data(),0,nullptr)!=0,"kinematic contact nodes");
    require(api->surface(world,4,movingPlaneMaps,1,&face,.01,0)!=0,"kinematic contact surface");
    require(api->targets(world,3,planeTargets)!=0,"kinematic plane targets");
    require(api->step(world,&cfg)!=0,"kinematic plane contact step");
    require(api->download(world,4,movingPlaneNodes.data())!=0,"kinematic contact result");
    require(movingPlaneNodes[0].position.y>.00019,"contact separation includes kinematic displacement");
    require(std::fabs(movingPlaneNodes[0].velocity.y-1.)<1.e-4,"kinematic plane transfers its normal velocity");
    require(api->updateNodes(world,4,movingPlaneNodes.data())==0,"state replacement rejects active targets");
    auto edited=rest;for(auto& n:edited) n.position.y+=.001;
    btGpuXpbdSurfaceVertex editMaps[4]={{{0,1,2,3},{1,0,0,0},{0,0,0},0},remapMaps[1],remapMaps[2],remapMaps[3]};
    cfg={.0001,1,{0,0,0},0,0,1,0,0,1};
    require(api->upload(world,4,edited.data(),1,&tet)!=0,"state update mesh");
    require(api->surface(world,4,editMaps,1,&face,.01,0)!=0,"state update collision surface");
    require(api->step(world,&cfg)!=0,"state update captures graph");
    for(auto& n:edited) { n.position.y+=.1;n.velocity={0,1,0};n.inverseMass=8; }
    require(api->updateNodes(world,4,edited.data())!=0,"replace state without recooking");
    btGpuXpbdVec editedSurface[4];
    require(api->surfacePositions(world,4,editedSurface)!=0,"state update retains collision mapping");
    require(distance(editedSurface[0],edited[0].position)<1.e-6,"skinning uses replaced node state");
    require(api->download(world,4,edited.data())!=0,"quantized edited state");
    auto invalidState=edited;invalidState.back().position.x=std::numeric_limits<double>::quiet_NaN();
    require(api->updateNodes(world,4,invalidState.data())==0,"reject invalid edited state");
    state.resize(4);require(api->download(world,4,state.data())!=0,"download after invalid edit");
    require(std::memcmp(state.data(),edited.data(),4*sizeof(btGpuXpbdNode))==0,"invalid edit preserves prior GPU state");
    require(api->step(world,&cfg)!=0,"replay graph after state update");
    require(api->download(world,4,state.data())!=0,"download after edited step");
    for(int i=0;i<4;++i)
    {
        auto expected=edited[i].position;expected.y+=cfg.dt;
        require(distance(state[i].position,expected)<1.e-6&&state[i].inverseMass==8,"edited graph uses updated positions velocities and masses");
    }
    edited=state;for(auto& n:edited) n.velocity={0,0,0};edited[1].position.x+=.1;
    require(api->updateNodes(world,4,edited.data())!=0,"deform edited state");cfg.dt=.001;
    require(api->step(world,&cfg)!=0,"solve deformation after state edit");
    require(api->download(world,4,state.data())!=0,"read edited deformation");
    require(state[1].position.x<edited[1].position.x-1.e-5,"state replacement preserves original elastic rest shape");
    btGpuXpbdTarget editTarget={0,state[0].position};
    require(api->targets(world,1,&editTarget)!=0&&api->targets(world,0,nullptr)!=0,"edited state target roundtrip");
    require(api->download(world,4,state.data())!=0&&state[0].inverseMass==8,"target release restores edited inverse mass");
    auto accepted=rest;accepted[1].position.x=1.1;
    cfg={.001,1,{0,0,0},0,0,2,0,0,0};
    require(api->upload(world,4,accepted.data(),1,&tet)!=0,"rollback mesh");
    require(api->step(world,&cfg)!=0,"rollback accepted step");
    require(api->download(world,4,accepted.data())!=0,"rollback accepted state");
    btGpuXpbdVec extreme[4]={{std::numeric_limits<SolverScalar>::max(),0,0},{0,0,0},{0,0,0},{0,0,0}};
    require(api->forces(world,4,extreme)!=0,"finite overflow force accepted");
    require(api->step(world,&cfg)==0,"nonfinite computed state rejects step");
    state.resize(4);
    require(api->download(world,4,state.data())!=0,"download rejected step");
    require(std::memcmp(state.data(),accepted.data(),4*sizeof(btGpuXpbdNode))==0,"rejected step restores exact accepted nodes");
    btGpuXpbdVec zeroForces[4]={};
    require(api->forces(world,4,zeroForces)!=0,"clear overflow force");
    require(api->step(world,&cfg)!=0,"retry after rejected step");
    require(api->download(world,4,state.data())!=0,"retry state");
    auto retry=state;
    require(api->upload(world,4,accepted.data(),1,&tet)!=0,"clean recovery comparison");
    require(api->step(world,&cfg)!=0,"clean recovery step");
    require(api->download(world,4,state.data())!=0,"clean recovery state");
    for(int i=0;i<4;++i)
        require(distance(state[i].position,retry[i].position)<1.e-6 && distance(state[i].velocity,retry[i].velocity)<1.e-4,
            "retry matches uncorrupted simulation");
    require(api->upload(world,-1,nullptr,0,nullptr)==0,"reject invalid input");
#ifndef XPBD_TEST_RUNTIME_LOADER
    api->destroy(world);
#endif
    std::printf("GPU XPBD core: %d failures\n",failures);return failures?1:0;
}
