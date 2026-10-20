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
void require(bool value,const char* label) { if(!value) { std::fprintf(stderr,"FAIL: %s\n",label);++failures; } }
float distance(btGpuXpbdVec a,btGpuXpbdVec b)
{ return std::sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z)); }
void testFriction(const btGpuXpbdTgsApi* api,void* world)
{
    auto run=[&](float mu,float speed,float height,float push,int velocityIterations,bool movingGround)
    {
        btGpuXpbdNode nodes[]={{{0,height,0},1,{speed,0,0},0},
            {{-20,0,-20},movingGround?1.f:0.f,{0,0,0},0},
            {{0,0,20},movingGround?1.f:0.f,{0,0,0},0},
            {{20,0,-20},movingGround?1.f:0.f,{0,0,0},0}};
        btGpuXpbdSurfaceVertex maps[4];
        for(int i=0;i<4;++i)
            maps[i]=i&&!movingGround?btGpuXpbdSurfaceVertex{{-1,-1,-1,-1},{0,0,0,0},nodes[i].position,1}:
                btGpuXpbdSurfaceVertex{{i,-1,-1,-1},{1,0,0,0},{0,0,0},i?1:0};
        btGpuXpbdSurfaceTriangle face={{1,2,3}};
        require(api->upload(world,4,nodes,0,nullptr)!=0,"friction nodes");
        require(api->surface(world,4,maps,1,&face,.1f,0)!=0,"friction surface");
        float coefficients[]={mu,1};
        require(api->friction(world,2,coefficients)!=0,"friction coefficients");
        btGpuXpbdVec forces[]={{push,-10,0},{0,0,0},{0,0,0},{0,0,0}};
        require(api->forces(world,4,forces)!=0,"friction load");
        btGpuXpbdStep step={.001f,1,{0,0,0},0,0,4,0,0,velocityIterations};
        for(int i=0;i<100;++i) require(api->step(world,&step)!=0,"friction step");
        require(api->download(world,4,nodes)!=0,"friction result");
        if(movingGround)
        {
            btGpuXpbdVec momentum={0,0,0};
            for(auto& n:nodes) { momentum.x+=n.velocity.x;momentum.y+=n.velocity.y;momentum.z+=n.velocity.z; }
            require(distance(momentum,{speed,-1,0})<2.e-4f,"friction preserves two-body momentum");
            require(nodes[1].velocity.x>0,"friction transfers tangent momentum to second body");
        }
        return nodes[0];
    };
    auto free=run(0,1,0,0,1,false);
    std::printf("friction free vx=%g x=%g\n",free.velocity.x,free.position.x);
    require(std::fabs(free.velocity.x-1)<1.e-5f,"zero friction retains tangent velocity");
    auto sliding=run(.5f,1,0,0,1,false);
    std::printf("friction sliding vx=%g x=%g\n",sliding.velocity.x,sliding.position.x);
    require(std::fabs(sliding.velocity.x-.5f)<.025f,"sliding deceleration matches Coulomb friction");
    auto manyPasses=run(.5f,1,0,0,4,false);
    require(std::fabs(manyPasses.velocity.x-sliding.velocity.x)<.01f,"velocity iterations share one friction budget");
    auto sticking=run(.5f,0,0,2,1,false);
    require(std::fabs(sticking.velocity.x)<.002f&&std::fabs(sticking.position.x)<.0002f,"sub-limit lateral load sticks");
    float noFriction[]={0,1};
    require(api->friction(world,2,noFriction)!=0,"friction changes without mesh upload");
    btGpuXpbdStep changedStep={.001f,1,{0,0,0},0,0,4,0,0,1};
    require(api->step(world,&changedStep)!=0,"step after friction change");
    btGpuXpbdNode changed[4];api->download(world,4,changed);
    require(changed[0].velocity.x>.0015f,"disabling friction immediately releases sticking");
    auto airborne=run(.5f,1,1,0,1,false);
    require(std::fabs(airborne.velocity.x-1)<1.e-5f,"friction does not act across open gaps");
    run(.5f,1,0,0,1,true);
    float invalid[]={-1};
    require(!api->friction(world,1,invalid),"negative friction rejected");
    invalid[0]=std::numeric_limits<float>::quiet_NaN();
    require(!api->friction(world,1,invalid),"nonfinite friction rejected");
}
void testMouseSprings(const btGpuXpbdTgsApi* api,void* world)
{
    btGpuXpbdNode node={{0,0,0},1,{0,0,0},0};
    require(api->upload(world,1,&node,0,nullptr)!=0,"mouse spring node");
    btGpuXpbdSpring spring={0,{1,0,0},100,10,1000};
    btGpuXpbdStep step={.1f,1,{0,0,0},0,0,1};
    require(api->springs(world,1,&spring)!=0,"mouse spring upload");
    require(api->step(world,&step)!=0,"mouse spring step");
    api->download(world,1,&node);
    require(std::fabs(node.velocity.x-10.f/3)<1.e-5f,"implicit spring matches analytic damped step");
    require(std::fabs(node.position.x-1.f/3)<1.e-5f,"spring moves picked node");
    spring.target={-1,0,0};
    require(api->springs(world,1,&spring)!=0,"mouse target changes without graph rebuild");
    require(api->step(world,&step)!=0,"moved target step");
    api->download(world,1,&node);
    require(node.velocity.x<0,"picked node follows moved target");
    require(api->springs(world,0,nullptr)!=0,"mouse release");
    float releaseVelocity=node.velocity.x;
    require(api->step(world,&step)!=0,"released node step");api->download(world,1,&node);
    require(std::fabs(node.velocity.x-releaseVelocity)<1.e-6f,"release preserves velocity");
    node={{0,0,0},1,{0,0,0},0};
    api->updateNodes(world,1,&node);
    spring={0,{100,0,0},1.e12f,0,.5f};
    require(api->springs(world,1,&spring)!=0,"force-limited spring");
    require(api->step(world,&step)!=0,"stiff capped spring remains finite");api->download(world,1,&node);
    require(std::fabs(node.velocity.x-.05f)<1.e-6f,"elastic mouse force respects per-node cap");
    node.inverseMass=0;api->updateNodes(world,1,&node);
    auto pinned=node.position;
    api->step(world,&step);api->download(world,1,&node);
    require(distance(node.position,pinned)==0,"mouse spring does not move a fixed node");
    spring.node=1;require(!api->springs(world,1,&spring),"invalid picked node rejected");
    spring.node=0;spring.damping=-1;require(!api->springs(world,1,&spring),"invalid picking damping rejected");
}
}
int main()
{
#ifdef XPBD_TEST_RUNTIME_LOADER
    btGpuXpbdTgsBackend backend;
    if(!backend.ready()) { std::fprintf(stderr,"%s\n",backend.error().c_str());return 1; }
    const auto* api=backend.api();void* world=backend.state();
#else
    const auto* api=btGetGpuXpbdTgsApi(5);void* world=api->create();
#endif
    if(!world) { std::fprintf(stderr,"%s\n",api->error(nullptr));return 1; }
    testFriction(api,world);
    testMouseSprings(api,world);
    btGpuXpbdTet tet={{0,1,2,3},{1,0,0,0,1,0,0,0,1},1.f/6,10000,.4f,.01f};
    std::vector<btGpuXpbdNode> rest={{{0,0,0},4,{0,0,0},0},{{1,0,0},4,{0,0,0},0},{{0,1,0},4,{0,0,0},0},{{0,0,1},4,{0,0,0},0}};
    btGpuXpbdStep cfg={.001f,10,{0,0,0},0,0,1};
    auto run=[&](std::vector<btGpuXpbdNode> input,int steps)
    {
        require(api->upload(world,4,input.data(),1,&tet)!=0,"upload");
        for(int k=0;k<steps;++k) if(!api->step(world,&cfg)) { require(false,api->error(world));break; }
        require(api->download(world,4,input.data())!=0,"download");return input;
    };
    auto state=run(rest,100);
    for(int i=0;i<4;++i) require(distance(state[i].position,rest[i].position)<1.e-6f,"rest equilibrium");
    btGpuXpbdVec samples[2]={{.2f,.3f,.1f},{1.2f,.1f,.1f}};
    btGpuXpbdSurfaceVertex embedded[2];
    require(api->embed(world,2,samples,0,4,7,embedded)!=0,"GPU embedding");
    for(int i=0;i<2;++i)
    {
        btGpuXpbdVec p={0,0,0};float sum=0;
        for(int j=0;j<4;++j) { auto v=rest[embedded[i].nodes[j]].position;float w=embedded[i].weights[j];sum+=w;p.x+=v.x*w;p.y+=v.y*w;p.z+=v.z*w; }
        require(distance(p,samples[i])<1.e-6f&&std::fabs(sum-1)<1.e-6f,"embedding reproduces interior and exterior points");
    }
    auto translated=rest;for(auto& n:translated) n.velocity={1,2,-3};
    state=run(translated,100);
    for(int i=0;i<4;++i)
    {
        auto p=rest[i].position;p.x+=.1f;p.y+=.2f;p.z-=.3f;
        require(distance(state[i].position,p)<1.e-4f,"uniform translation");
        require(distance(state[i].velocity,{1,2,-3})<1.e-3f,"uniform velocity");
    }
    auto rotated=rest;for(auto& n:rotated) { auto p=n.position;n.position={p.x*std::cos(.7f)-p.y*std::sin(.7f),p.x*std::sin(.7f)+p.y*std::cos(.7f),p.z}; }
    state=run(rotated,100);
    for(int i=0;i<4;++i) require(distance(state[i].position,rotated[i].position)<1.e-4f,"rigid rotation equilibrium");
    auto strained=rest;strained[1].position.x=1.2f;
    state=run(strained,100);
    btGpuXpbdVec momentum={0,0,0},com={0,0,0};
    for(auto& n:state) { momentum.x+=n.velocity.x/4;momentum.y+=n.velocity.y/4;momentum.z+=n.velocity.z/4;com.x+=n.position.x/4;com.y+=n.position.y/4;com.z+=n.position.z/4; }
    require(distance(momentum,{0,0,0})<1.e-4f,"internal forces preserve linear momentum");
    require(distance(com,{.3f,.25f,.25f})<1.e-4f,"internal forces preserve COM");
    require(std::fabs(distance(state[0].position,state[1].position)-1)<.05f,"elastic recovery");
    auto pinned=strained;pinned[0].inverseMass=0;state=run(pinned,100);
    require(distance(state[0].position,pinned[0].position)==0,"pinned node");
    auto local=rest;local[1].position.x=1.25f;
    auto savedCfg=cfg;cfg.dt=.0001f;cfg.substeps=1;
    auto nearOrigin=run(local,1);
    for(auto& n:local) { n.position.x+=4096;n.position.y+=4096;n.position.z+=4096; }
    auto farOrigin=run(local,1);
    for(int i=0;i<4;++i)
    {
        float error=distance(nearOrigin[i].velocity,farOrigin[i].velocity);
        if(error>1.e-5f) std::fprintf(stderr,"translation material velocity error node %d: %.9g\n",i,error);
        require(error<1.e-5f,"material corrections independent of world translation");
    }
    cfg=savedCfg;
    std::vector<btGpuXpbdNode> contactNodes={{{0,.01f,0},1,{0,-1,0},0},{{-2,0,-2},1,{0,0,0},0},{{0,0,2},1,{0,0,0},0},{{2,0,-2},1,{0,0,0},0}};
    cfg.velocityIterations=1;
    std::vector<btGpuXpbdSurfaceVertex> maps(4);
    for(int i=0;i<4;++i) maps[i]={{i,-1,-1,-1},{1,0,0,0},{0,0,0},i==0?0:1};
    btGpuXpbdSurfaceTriangle face={{1,2,3}};
    require(api->upload(world,4,contactNodes.data(),0,nullptr)!=0,"contact nodes");
    require(api->surface(world,4,maps.data(),1,&face,.1f,0)!=0,"dynamic surface upload");
    for(int k=0;k<100;++k) require(api->step(world,&cfg)!=0,"dynamic contact step");
    state=contactNodes;require(api->download(world,4,state.data())!=0,"contact download");
    momentum={0,0,0};for(auto& n:state) { momentum.x+=n.velocity.x;momentum.y+=n.velocity.y;momentum.z+=n.velocity.z; }
    require(distance(momentum,{0,-1,0})<1.e-4f,"dynamic contact transfers equal and opposite momentum");
    require(state[0].velocity.y>-.99f && state[1].velocity.y<0,"both bodies respond to contact");
    for(int i=1;i<4;++i) maps[i]={{-1,-1,-1,-1},{0,0,0,0},contactNodes[i].position,1};
    require(api->upload(world,4,contactNodes.data(),0,nullptr)!=0,"reset contact nodes");
    require(api->surface(world,4,maps.data(),1,&face,.1f,0)!=0,"static surface upload");
    for(int k=0;k<100;++k) require(api->step(world,&cfg)!=0,"static contact step");
    require(api->download(world,4,state.data())!=0,"static contact download");
    require(state[0].position.y>=-.001f,"static triangle prevents penetration");
    btGpuXpbdVec surfacePositions[4];require(api->surfacePositions(world,4,surfacePositions)!=0,"GPU skin download");
    require(distance(surfacePositions[0],state[0].position)<1.e-7f,"dynamic GPU skinning");
    require(distance(surfacePositions[1],contactNodes[1].position)<1.e-7f,"static GPU skinning");
    savedCfg=cfg;cfg.dt=.0001f;cfg.substeps=1;cfg.iterations=1;
    contactNodes[0].position={0,.00015f,0};contactNodes[0].velocity={0,-1,0};
    require(api->upload(world,4,contactNodes.data(),0,nullptr)!=0,"look-ahead contact nodes");
    require(api->surface(world,4,maps.data(),1,&face,.1f,0)!=0,"look-ahead contact surface");
    require(api->step(world,&cfg)!=0,"look-ahead contact step");api->download(world,4,state.data());
    require(std::fabs(state[0].position.y-.00005f)<1.e-7f,"velocity contact projection preserves position");
    require(std::fabs(state[0].velocity.y+.5f)<1.e-4f,"velocity pass limits predicted penetration across an open gap");
    contactNodes[0].position={0,.0004f,0};
    require(api->upload(world,4,contactNodes.data(),0,nullptr)!=0,"open gap nodes");
    require(api->surface(world,4,maps.data(),1,&face,.1f,0)!=0,"open gap surface");
    require(api->step(world,&cfg)!=0,"open gap step");api->download(world,4,state.data());
    require(std::fabs(state[0].velocity.y+1)<1.e-5f,"velocity pass leaves a non-closing gap unchanged");
    contactNodes[0].position={0,.00005f,0};contactNodes[0].velocity={0,1,0};
    require(api->upload(world,4,contactNodes.data(),0,nullptr)!=0,"separating nodes");
    require(api->surface(world,4,maps.data(),1,&face,.1f,0)!=0,"separating surface");
    require(api->step(world,&cfg)!=0,"separating step");api->download(world,4,state.data());
    require(std::fabs(state[0].velocity.y-1)<1.e-5f,"velocity pass does not attract separating bodies");
    cfg=savedCfg;
    contactNodes[0].position={2.01f,-.001f,-2};contactNodes[0].velocity={0,0,0};
    require(api->upload(world,4,contactNodes.data(),0,nullptr)!=0,"edge contact nodes");
    require(api->surface(world,4,maps.data(),1,&face,.1f,0)!=0,"edge contact surface");
    require(api->step(world,&cfg)!=0,"edge contact step");api->download(world,4,state.data());
    require(std::fabs(state[0].velocity.x)<1.e-5f&&std::fabs(state[0].velocity.z)<1.e-5f,"frictionless face contact does not pull sideways behind an edge");
    require(api->upload(world,4,rest.data(),0,nullptr)!=0,"drive nodes");
    btGpuXpbdTarget target={0,{1,2,3}};
    require(api->targets(world,1,&target)!=0,"drive target");
    cfg.iterations=3;
    require(api->step(world,&cfg)!=0,"drive step");
    require(api->download(world,4,state.data())!=0,"drive download");
    require(distance(state[0].position,target.position)<1.e-6f&&state[0].inverseMass==0,"drive reaches target");
    target.position={2,3,4};require(api->targets(world,1,&target)!=0,"update target");
    require(api->step(world,&cfg)!=0,"replayed drive graph");api->download(world,4,state.data());
    require(distance(state[0].position,target.position)<1.e-6f,"graph uses updated target buffer");
    require(api->targets(world,0,nullptr)!=0,"release target");api->download(world,4,state.data());
    require(state[0].inverseMass==rest[0].inverseMass,"release restores mass");
    std::vector<btGpuXpbdNode> crossed={{{-2,.00015f,-.1f},1,{0,-1,0},0},{{2,.00015f,-.1f},1,{0,-1,0},0},{{2,.00015f,.1f},1,{0,-1,0},0},
        {{-.1f,0,-2},1,{0,0,0},0},{{0,0,2},1,{0,0,0},0},{{.1f,0,-2},1,{0,0,0},0}};
    btGpuXpbdSurfaceTriangle crossingFaces[2]={{{0,1,2}},{{3,4,5}}};
    std::vector<btGpuXpbdSurfaceVertex> crossingMaps(6);
    for(int i=0;i<6;++i) crossingMaps[i]={{i,-1,-1,-1},{1,0,0,0},{0,0,0},i<3?0:1};
    cfg={.0001f,1,{0,0,0},0,0,1,0,0,1};
    require(api->upload(world,6,crossed.data(),0,nullptr)!=0,"crossing edge nodes");
    require(api->surface(world,6,crossingMaps.data(),2,crossingFaces,.01f,0)!=0,"crossing edge surfaces");
    require(api->step(world,&cfg)!=0,"crossing edge step");
    auto crossedResult=crossed;api->download(world,6,crossedResult.data());
    float sumVelocity=0,otherVelocity=0;
    for(int i=0;i<6;++i) { sumVelocity+=crossedResult[i].velocity.y;if(i>=3) otherVelocity+=crossedResult[i].velocity.y; }
    require(otherVelocity<-.01f,"crossing edges transfer momentum without vertex-face proximity");
    require(std::fabs(sumVelocity+3)<1.e-4f,"triangle contacts preserve equal and opposite momentum");
    for(int i=3;i<6;++i) crossingMaps[i]={{-1,-1,-1,-1},{0,0,0,0},crossed[i].position,1};
    require(api->upload(world,6,crossed.data(),0,nullptr)!=0,"static crossing edge nodes");
    require(api->surface(world,6,crossingMaps.data(),2,crossingFaces,.01f,0)!=0,"static crossing edge surfaces");
    require(api->step(world,&cfg)!=0,"static crossing edge step");api->download(world,6,crossedResult.data());
    sumVelocity=0;for(int i=0;i<3;++i) sumVelocity+=crossedResult[i].velocity.y;
    require(sumVelocity>-2.99f,"static triangle catches crossing edges");
    auto crossingStep=crossed;
    for(int i=0;i<3;++i) crossingStep[i].position.y=.00005f;
    require(api->upload(world,6,crossingStep.data(),0,nullptr)!=0,"first increment crossing nodes");
    require(api->surface(world,6,crossingMaps.data(),2,crossingFaces,.01f,0)!=0,"first increment crossing plane");
    require(api->step(world,&cfg)!=0,"first increment crossing step");api->download(world,6,crossedResult.data());
    sumVelocity=0;for(int i=0;i<3;++i) sumVelocity+=crossedResult[i].velocity.y;
    require(sumVelocity>-2.99f,"contact prepared before integration survives crossing static face");
    crossed[0].position.y=.0033f;crossed[1].position.y=crossed[2].position.y=-.003f;
    require(api->upload(world,6,crossed.data(),0,nullptr)!=0,"backface triangle nodes");
    require(api->surface(world,6,crossingMaps.data(),2,crossingFaces,.01f,0)!=0,"backface triangle surfaces");
    require(api->step(world,&cfg)!=0,"backface triangle step");api->download(world,6,crossedResult.data());
    for(int i=0;i<3;++i) require(std::fabs(crossedResult[i].velocity.y+1)<1.e-5f,"static mesh culls triangle centroid behind face");
    auto remapNodes=rest;for(auto& n:remapNodes) n.velocity={0,-1,0};
    btGpuXpbdSurfaceVertex remapMaps[4]={{{0,1,2,3},{-.2f,1.2f,0,0},{0,.00015f,0},0},
        {{-1,-1,-1,-1},{0,0,0,0},{-4,0,-4},1},{{-1,-1,-1,-1},{0,0,0,0},{0,0,4},1},{{-1,-1,-1,-1},{0,0,0,0},{4,0,-4},1}};
    require(api->upload(world,4,remapNodes.data(),1,&tet)!=0,"contact remap tetrahedron");
    require(api->surface(world,4,remapMaps,1,&face,.01f,0)!=0,"contact remap surface");
    require(api->step(world,&cfg)!=0,"contact remap step");api->download(world,4,remapNodes.data());
    require(std::fabs(remapNodes[1].velocity.y+.5f)<1.e-3f,"outside contact maps to nearest tetrahedron point");
    require(std::fabs(remapNodes[0].velocity.y+1)<1.e-3f,"contact remap does not extrapolate negative force weights");
    // The entire collision triangle lies outside the tet AABB. Cooking must
    // retain its embedding tet, rather than falling back to negative weights.
    btGpuXpbdSurfaceVertex outsideMaps[6]={remapMaps[0],remapMaps[0],remapMaps[0],remapMaps[1],remapMaps[2],remapMaps[3]};
    outsideMaps[0].offset.z=-.001f;outsideMaps[1].offset.x=.001f;outsideMaps[2].offset.z=.001f;
    remapNodes=rest;for(auto& n:remapNodes) n.velocity={0,-1,0};
    require(api->upload(world,4,remapNodes.data(),1,&tet)!=0,"outside triangle tet");
    require(api->surface(world,6,outsideMaps,2,crossingFaces,.01f,0)!=0,"outside triangle candidates");
    require(api->step(world,&cfg)!=0,"outside triangle contact step");api->download(world,4,remapNodes.data());
    require(std::fabs(remapNodes[0].velocity.y+1)<1.e-3f,"outside triangle retains clamped embedding tet");
    require(remapNodes[1].velocity.y>-.99f,"outside triangle transfers contact through cooked candidate");
    // A nearby part of the same body must not steal an embedded surface's contact.
    auto ownershipNodes=rest;ownershipNodes.insert(ownershipNodes.end(),rest.begin(),rest.end());
    btGpuXpbdTet ownershipTets[2]={tet,tet};for(int& n:ownershipTets[1].nodes) n+=4;
    for(int i=0;i<4;++i)
    {
        ownershipNodes[i].position.y+=.00005f;ownershipNodes[i].velocity={0,-1,0};
        ownershipNodes[i+4].position.y-=.1f;
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
    require(api->surface(world,7,ownershipMaps,2,crossingFaces,.01f,0)!=0,"contact material ownership surface");
    require(api->step(world,&cfg)!=0,"contact material ownership step");
    require(api->download(world,8,ownershipNodes.data())!=0,"contact material ownership download");
    float ownerY=0;for(int i=0;i<4;++i) ownerY+=ownershipNodes[i].velocity.y;
    require(ownerY>-3.99f,"contact acts on the embedded surface material");
    for(int i=4;i<8;++i) require(distance(ownershipNodes[i].velocity,{0,0,0})<1.e-5f,"nearby material does not steal contact");
    contactNodes[0].position={0,-.001f,0};contactNodes[0].velocity={0,0,0};cfg={.002f,1,{0,0,0},0,0,2,0,0,1};
    require(api->upload(world,4,contactNodes.data(),0,nullptr)!=0,"penetration recovery nodes");
    require(api->surface(world,4,maps.data(),1,&face,.1f,0)!=0,"penetration recovery plane");
    require(api->step(world,&cfg)!=0,"penetration recovery step");api->download(world,4,state.data());
    require(std::fabs(state[0].velocity.y)<1.e-5f,"static contact releases excess accumulated impulse");
    require(std::fabs(state[0].position.y)<1.e-6f,"penetration recovery does not launch a resting node");
    // Exercise the EPA/degenerate-simplex paths, not only separated GJK pairs.
    for(int mode=0;mode<3;++mode)
    {
        auto pair=crossed;
        for(int i=0;i<6;++i)
        {
            pair[i].velocity={0,0,0};pair[i].position.y=0;
            crossingMaps[i]={{i,-1,-1,-1},{1,0,0,0},{0,0,0},i<3?0:1};
        }
        if(mode==1) { pair[0].position.y=-.01f;pair[1].position.y=pair[2].position.y=.01f; }
        if(mode==2) { pair[1].position=pair[0].position; } // collapsed face is ignored
        cfg={.0001f,1,{0,0,0},0,0,2,0,0,1};
        require(api->upload(world,6,pair.data(),0,nullptr)!=0,"EPA test upload");
        require(api->surface(world,6,crossingMaps.data(),2,crossingFaces,.01f,0)!=0,"EPA test surfaces");
        require(api->step(world,&cfg)!=0,"EPA and coplanar contact step");
        require(api->download(world,6,crossedResult.data())!=0,"EPA test download");
        btGpuXpbdVec momentum={0,0,0};float speed=0;
        for(const auto& node:crossedResult)
        {
            require(std::isfinite(node.position.x)&&std::isfinite(node.position.y)&&std::isfinite(node.position.z)&&
                std::isfinite(node.velocity.x)&&std::isfinite(node.velocity.y)&&std::isfinite(node.velocity.z),"EPA finite state");
            momentum.x+=node.velocity.x;momentum.y+=node.velocity.y;momentum.z+=node.velocity.z;
            speed+=distance(node.velocity,{0,0,0});
        }
        require(distance(momentum,{0,0,0})<1.e-3f,"EPA contacts preserve total momentum");
        if(mode!=1) require(speed<1.e-3f,"coplanar or collapsed contact does not launch resting nodes");
    }
    std::vector<btGpuXpbdNode> tinySupport(2,{{0,.00015f,0},1,{0,-1,0},0});
    btGpuXpbdSurfaceVertex tinyMaps[4]={{{0,1,-1,-1},{.0005f,.9995f,0,0},{0,0,0},0},remapMaps[1],remapMaps[2],remapMaps[3]};
    cfg={.0001f,1,{0,0,0},0,0,1,0,0,1};
    require(api->upload(world,2,tinySupport.data(),0,nullptr)!=0,"tiny support nodes");
    require(api->surface(world,4,tinyMaps,1,&face,.01f,0)!=0,"tiny support surface");
    require(api->step(world,&cfg)!=0,"tiny support contact step");api->download(world,2,tinySupport.data());
    require(std::fabs(tinySupport[0].velocity.y+1)<1.e-6f,"tiny weight does not create a contact reference");
    require(std::fabs(tinySupport[1].velocity.y-(-1.f+.5f/.9995f))<1.e-4f,"contact denominator uses counted support");
    auto preciseNodes=rest;
    for(auto& n:preciseNodes) { n.position.y+=1024.0001220703125f;n.velocity={0,-1,0}; }
    btGpuXpbdSurfaceVertex preciseMaps[4]={{{0,1,2,3},{1,0,0,0},{0,0,0},0},remapMaps[1],remapMaps[2],remapMaps[3]};
    for(int i=1;i<4;++i) preciseMaps[i].offset.y=1024;
    cfg={.0002f,1,{0,0,0},0,0,200,0,0,1};
    require(api->upload(world,4,preciseNodes.data(),1,&tet)!=0,"sub-ULP contact nodes");
    require(api->surface(world,4,preciseMaps,1,&face,.01f,0)!=0,"sub-ULP contact plane");
    require(api->step(world,&cfg)!=0,"sub-ULP accumulated contact step");api->download(world,4,preciseNodes.data());
    require(preciseNodes[0].velocity.y>-.5f,"contact resolves motion smaller than world-position precision");
    auto volumes=rest;volumes.insert(volumes.end(),rest.begin(),rest.end());
    for(int i=0;i<4;++i) { volumes[i].position.y+=.00015f;volumes[i].velocity={0,-1,0};volumes[i+4].position.y=-volumes[i+4].position.y; }
    btGpuXpbdTet volumeTets[2]={tet,tet};
    for(int& n:volumeTets[1].nodes) n+=4;volumeTets[1].inverseRest[4]=-1;
    btGpuXpbdSurfaceVertex volumeMaps[6];int volumeIds[6]={0,1,3,4,7,5};
    for(int i=0;i<6;++i)
    {
        int body=i<3?0:1;volumeMaps[i]={{body*4,body*4+1,body*4+2,body*4+3},{0,0,0,0},{0,0,0},body};
        volumeMaps[i].weights[volumeIds[i]-body*4]=1;
    }
    cfg={.0001f,1,{0,0,0},0,0,1,0,0,1};
    require(api->upload(world,8,volumes.data(),2,volumeTets)!=0,"volume vertex contact upload");
    require(api->surface(world,6,volumeMaps,2,crossingFaces,.01f,0)!=0,"volume vertex contact surfaces");
    require(api->step(world,&cfg)!=0,"volume vertex contact step");api->download(world,8,volumes.data());
    float totalY=0,lowerY=0;
    for(int i=0;i<8;++i) { totalY+=volumes[i].velocity.y;if(i>=4) lowerY+=volumes[i].velocity.y; }
    require(lowerY<-.01f,"signed deformable vertex contacts transfer momentum");
    require(std::fabs(totalY+4)<1.e-3f,"signed deformable vertex contacts conserve momentum");
    // A collision surface outside its simulation tetrahedron has signed support.
    // Dropping its negative contributions would create net contact momentum.
    auto extrapolatedVolumes=rest;extrapolatedVolumes.insert(extrapolatedVolumes.end(),rest.begin(),rest.end());
    btGpuXpbdSurfaceVertex extrapolatedMaps[6];
    std::memcpy(extrapolatedMaps,volumeMaps,sizeof(volumeMaps));
    for(int i=0;i<4;++i)
    {
        extrapolatedVolumes[i].position.y+=.04995f;
        extrapolatedVolumes[i].velocity={0,-1,0};
        extrapolatedVolumes[i+4].position.y=-extrapolatedVolumes[i+4].position.y;
    }
    for(int i=0;i<3;++i) { extrapolatedMaps[i].weights[0]+=.05f;extrapolatedMaps[i].weights[2]-=.05f; }
    require(api->upload(world,8,extrapolatedVolumes.data(),2,volumeTets)!=0,"extrapolated volume contact upload");
    require(api->surface(world,6,extrapolatedMaps,2,crossingFaces,.01f,0)!=0,"extrapolated volume contact surface");
    require(api->step(world,&cfg)!=0,"extrapolated volume contact step");
    require(api->download(world,8,extrapolatedVolumes.data())!=0,"extrapolated volume contact download");
    totalY=0;lowerY=0;
    for(int i=0;i<8;++i) { totalY+=extrapolatedVolumes[i].velocity.y;if(i>=4) lowerY+=extrapolatedVolumes[i].velocity.y; }
    require(lowerY<-.01f,"extrapolated volume contacts transfer momentum");
    require(std::fabs(totalY+4)<1.e-3f,"extrapolated volume contacts conserve momentum");
    auto checkFilter=[&](const std::vector<btGpuXpbdNode>& input,int tetCount,const btGpuXpbdTet* tets,
        int vertexCount,const btGpuXpbdSurfaceVertex* vertices,int triangleCount,const btGpuXpbdSurfaceTriangle* triangles)
    {
        const int count=int(input.size());auto free=input,filtered=input;
        cfg={.0001f,1,{0,0,0},0,0,1,0,0,1};
        require(api->upload(world,count,input.data(),tetCount,tets)!=0,"pair filter free mesh");
        require(api->step(world,&cfg)!=0,"pair filter free step");
        require(api->download(world,count,free.data())!=0,"pair filter free state");
        require(api->upload(world,count,input.data(),tetCount,tets)!=0,"pair filter mesh");
        require(api->surface(world,vertexCount,vertices,triangleCount,triangles,.01f,0)!=0,"pair filter surface");
        btGpuXpbdDisabledPair disabled[]={{1,0},{0,1},{4,3},{2,1}};
        require(api->disabledPairs(world,4,disabled)!=0,"pair filter normalizes unordered duplicates");
        require(api->step(world,&cfg)!=0,"pair filter disabled step");
        require(api->download(world,count,filtered.data())!=0,"pair filter disabled state");
        for(int i=0;i<count;++i)
            require(distance(free[i].position,filtered[i].position)<1.e-7f && distance(free[i].velocity,filtered[i].velocity)<1.e-6f,
                "disabled GPU contacts match free motion");
        btGpuXpbdDisabledPair invalid={-1,0};
        require(api->disabledPairs(world,1,&invalid)==0,"pair filter rejects negative IDs");
        require(api->disabledPairs(world,0,nullptr)!=0,"pair filter reenable");
        require(api->step(world,&cfg)!=0,"pair filter reenabled step");
        require(api->download(world,count,filtered.data())!=0,"pair filter reenabled state");
        float change=0;for(int i=0;i<count;++i) change+=distance(free[i].velocity,filtered[i].velocity);
        require(change>.01f,"reenabling pair invalidates graph and restores contacts");
    };
    contactNodes[0].position={0,.00025f,0};contactNodes[0].velocity={0,-1,0};
    checkFilter(contactNodes,0,nullptr,4,maps.data(),1,&face);
    crossed[0].position.y=crossed[1].position.y=crossed[2].position.y=.00025f;
    for(int i=0;i<6;++i)
    {
        crossed[i].velocity=i<3?btGpuXpbdVec{0,-1,0}:btGpuXpbdVec{0,0,0};
        crossingMaps[i]={{i,-1,-1,-1},{1,0,0,0},{0,0,0},i<3?0:1};
    }
    checkFilter(crossed,0,nullptr,6,crossingMaps.data(),2,crossingFaces);
    for(int i=3;i<6;++i) crossingMaps[i]={{-1,-1,-1,-1},{0,0,0,0},crossed[i].position,1};
    checkFilter(crossed,0,nullptr,6,crossingMaps.data(),2,crossingFaces);
    volumes=rest;volumes.insert(volumes.end(),rest.begin(),rest.end());
    for(int i=0;i<4;++i) { volumes[i].position.y+=.00025f;volumes[i].velocity={0,-1,0};volumes[i+4].position.y=-volumes[i+4].position.y; }
    checkFilter(volumes,2,volumeTets,6,volumeMaps,2,crossingFaces);
    auto movingPlaneNodes=contactNodes;
    movingPlaneNodes[0].position={0,.00015f,0};movingPlaneNodes[0].velocity={0,0,0};
    btGpuXpbdSurfaceVertex movingPlaneMaps[4];btGpuXpbdTarget planeTargets[3];
    for(int i=0;i<4;++i)
    {
        movingPlaneMaps[i]={{i,-1,-1,-1},{1,0,0,0},{0,0,0},i==0?0:1};
        if(i) { planeTargets[i-1]={i,movingPlaneNodes[i].position};planeTargets[i-1].position.y=.0002f; }
    }
    cfg={.0002f,1,{0,0,0},0,0,2,0,0,1};
    require(api->upload(world,4,movingPlaneNodes.data(),0,nullptr)!=0,"kinematic contact nodes");
    require(api->surface(world,4,movingPlaneMaps,1,&face,.01f,0)!=0,"kinematic contact surface");
    require(api->targets(world,3,planeTargets)!=0,"kinematic plane targets");
    require(api->step(world,&cfg)!=0,"kinematic plane contact step");
    require(api->download(world,4,movingPlaneNodes.data())!=0,"kinematic contact result");
    require(movingPlaneNodes[0].position.y>.00019f,"contact separation includes kinematic displacement");
    require(std::fabs(movingPlaneNodes[0].velocity.y-1.f)<1.e-4f,"kinematic plane transfers its normal velocity");
    require(api->updateNodes(world,4,movingPlaneNodes.data())==0,"state replacement rejects active targets");
    auto edited=rest;for(auto& n:edited) n.position.y+=.001f;
    btGpuXpbdSurfaceVertex editMaps[4]={{{0,1,2,3},{1,0,0,0},{0,0,0},0},remapMaps[1],remapMaps[2],remapMaps[3]};
    cfg={.0001f,1,{0,0,0},0,0,1,0,0,1};
    require(api->upload(world,4,edited.data(),1,&tet)!=0,"state update mesh");
    require(api->surface(world,4,editMaps,1,&face,.01f,0)!=0,"state update collision surface");
    require(api->step(world,&cfg)!=0,"state update captures graph");
    for(auto& n:edited) { n.position.y+=.1f;n.velocity={0,1,0};n.inverseMass=8; }
    require(api->updateNodes(world,4,edited.data())!=0,"replace state without recooking");
    btGpuXpbdVec editedSurface[4];
    require(api->surfacePositions(world,4,editedSurface)!=0,"state update retains collision mapping");
    require(distance(editedSurface[0],edited[0].position)<1.e-6f,"skinning uses replaced node state");
    auto invalidState=edited;invalidState.back().position.x=std::numeric_limits<float>::quiet_NaN();
    require(api->updateNodes(world,4,invalidState.data())==0,"reject invalid edited state");
    state.resize(4);require(api->download(world,4,state.data())!=0,"download after invalid edit");
    require(std::memcmp(state.data(),edited.data(),4*sizeof(btGpuXpbdNode))==0,"invalid edit preserves prior GPU state");
    require(api->step(world,&cfg)!=0,"replay graph after state update");
    require(api->download(world,4,state.data())!=0,"download after edited step");
    for(int i=0;i<4;++i)
    {
        auto expected=edited[i].position;expected.y+=cfg.dt;
        require(distance(state[i].position,expected)<1.e-6f&&state[i].inverseMass==8,"edited graph uses updated positions velocities and masses");
    }
    edited=state;for(auto& n:edited) n.velocity={0,0,0};edited[1].position.x+=.1f;
    require(api->updateNodes(world,4,edited.data())!=0,"deform edited state");cfg.dt=.001f;
    require(api->step(world,&cfg)!=0,"solve deformation after state edit");
    require(api->download(world,4,state.data())!=0,"read edited deformation");
    require(state[1].position.x<edited[1].position.x-1.e-5f,"state replacement preserves original elastic rest shape");
    btGpuXpbdTarget editTarget={0,state[0].position};
    require(api->targets(world,1,&editTarget)!=0&&api->targets(world,0,nullptr)!=0,"edited state target roundtrip");
    require(api->download(world,4,state.data())!=0&&state[0].inverseMass==8,"target release restores edited inverse mass");
    auto accepted=rest;accepted[1].position.x=1.1f;
    cfg={.001f,1,{0,0,0},0,0,2,0,0,0};
    require(api->upload(world,4,accepted.data(),1,&tet)!=0,"rollback mesh");
    require(api->step(world,&cfg)!=0,"rollback accepted step");
    require(api->download(world,4,accepted.data())!=0,"rollback accepted state");
    btGpuXpbdVec extreme[4]={{std::numeric_limits<float>::max(),0,0},{0,0,0},{0,0,0},{0,0,0}};
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
        require(distance(state[i].position,retry[i].position)<1.e-6f && distance(state[i].velocity,retry[i].velocity)<1.e-4f,
            "retry matches uncorrupted simulation");
    require(api->upload(world,-1,nullptr,0,nullptr)==0,"reject invalid input");
#ifndef XPBD_TEST_RUNTIME_LOADER
    api->destroy(world);
#endif
    std::printf("GPU XPBD core: %d failures\n",failures);return failures?1:0;
}
