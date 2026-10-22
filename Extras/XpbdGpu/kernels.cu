// Contact/TGS integration adapted from PhysX softBodyGM.cu/deformableUtils.cuh.
// Copyright (c) 2001-2004 NovodeX AG. All rights reserved.
// Copyright (c) 2004-2008 AGEIA Technologies, Inc. All rights reserved.
// Copyright (c) 2008-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0
// Modified for Bullet: standalone POD storage, colored batches, CUDA kernels without PhysX dependencies.
// The StVK material projection below is derived from PositionBasedDynamics (MIT); see its separate license.

struct V { double x,y,z; };
struct Node { V p; double w; V v; double damping; };
struct Tet { int n[4]; double rest[9]; double volume,young,poisson,damping; };
struct M { V c[3]; };
__device__ V operator+(V a,V b) { return {a.x+b.x,a.y+b.y,a.z+b.z}; }
__device__ V operator-(V a,V b) { return {a.x-b.x,a.y-b.y,a.z-b.z}; }
__device__ V operator*(V a,double b) { return {a.x*b,a.y*b,a.z*b}; }
__device__ double dot(V a,V b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
__device__ V cross(V a,V b) { return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x}; }
__device__ V mul(M a,V b) { return a.c[0]*b.x+a.c[1]*b.y+a.c[2]*b.z; }
__device__ M mul(M a,M b) { return {mul(a,b.c[0]),mul(a,b.c[1]),mul(a,b.c[2])}; }
__device__ M gradient(M tensor,M rest)
{
    return {{mul(tensor,V{rest.c[0].x,rest.c[1].x,rest.c[2].x}),
             mul(tensor,V{rest.c[0].y,rest.c[1].y,rest.c[2].y}),
             mul(tensor,V{rest.c[0].z,rest.c[1].z,rest.c[2].z})}};
}
struct Spring { int node;V target;double stiffness,damping,maxForce; };
extern "C" __global__ void predict(Node* nodes,const V* force,V* accumulated,int count,double dt,V gravity,double forceDt,
    const Spring* springs,int springCount)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=count) return;
    Node p=nodes[i];if(p.w==0) return;
    p.v=p.v*exp(-p.damping*dt)+force[i]*(p.w*forceDt)+gravity*dt;
    for(int j=0;j<springCount;++j) if(springs[j].node==i)
    {
        Spring s=springs[j];
        // Compliant attachment: backward-Euler spring and isotropic damping, capped elastic force.
        double a=1.+dt*p.w*s.damping;
        V elastic=((s.target-p.p)*a-p.v*dt)*(s.stiffness/(a+dt*dt*p.w*s.stiffness));
        double length=sqrt(dot(elastic,elastic));
        if(length>s.maxForce) elastic=elastic*(s.maxForce/length);
        p.v=(p.v+elastic*(dt*p.w))*(1./a);
    }
    V dx=p.v*dt;p.p=p.p+dx;accumulated[i]=accumulated[i]+dx;nodes[i]=p;
}
// Split Rayleigh viscosity from elasticity. Damping acts on the accepted velocity,
// before adding this increment's external impulse; it therefore vanishes at static equilibrium.
// Six orthonormal symmetric strain modes diagonalize the isotropic constitutive tensor.
extern "C" __global__ void dampTets(Node* nodes,const Tet* tets,int begin,int count,double dt)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=count) return;
    Tet t=tets[begin+i];if(t.damping<=0) return;
    Node p[4];for(int j=0;j<4;++j) p[j]=nodes[t.n[j]];
    M rest={{{t.rest[0],t.rest[1],t.rest[2]},{t.rest[3],t.rest[4],t.rest[5]},{t.rest[6],t.rest[7],t.rest[8]}}};
    M f=mul(M{{p[1].p-p[0].p,p[2].p-p[0].p,p[3].p-p[0].p}},rest);
    double mu=t.young/(2*(1+t.poisson)),lambda=t.young*t.poisson/((1+t.poisson)*(1-2*t.poisson));
    const double a=0.57735026918962576451,b=0.70710678118654752440,c=0.40824829046386301637;
    M modes[6]={{{{a,0,0},{0,a,0},{0,0,a}}},{{{b,0,0},{0,-b,0},{0,0,0}}},
        {{{c,0,0},{0,c,0},{0,0,-2*c}}},{{{0,b,0},{b,0,0},{0,0,0}}},
        {{{0,0,b},{0,0,0},{b,0,0}}},{{{0,0,0},{0,0,b},{0,b,0}}}};
    for(int mode=0;mode<6;++mode)
    {
        M g=gradient(mul(f,modes[mode]),rest);
        V grad[4]={(g.c[0]+g.c[1]+g.c[2])*-1,g.c[0],g.c[1],g.c[2]};
        double weight=0,rate=0;
        for(int j=0;j<4;++j) { weight+=p[j].w*dot(grad[j],grad[j]);rate+=dot(grad[j],p[j].v); }
        double viscosity=t.damping*t.volume*(2*mu+(mode==0?3*lambda:0))*dt;
        double impulse=-viscosity*rate/(1+viscosity*weight);
        for(int j=0;j<4;++j) p[j].v=p[j].v+grad[j]*(p[j].w*impulse);
    }
    for(int j=0;j<4;++j) nodes[t.n[j]].v=p[j].v;
}
// StVK energy XPBD, following PositionBasedDynamics XPBD::solve_FEMTetraConstraint.
// See LICENSE-PositionBasedDynamics.txt. Each pass advances time, so its multiplier starts at zero.
extern "C" __global__ void solveTets(Node* nodes,const Tet* tets,V* accumulated,int begin,int count,double dt)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=count) return;
    Tet t=tets[begin+i];Node p[4];for(int j=0;j<4;++j) p[j]=nodes[t.n[j]];
    M rest={{{t.rest[0],t.rest[1],t.rest[2]},{t.rest[3],t.rest[4],t.rest[5]},{t.rest[6],t.rest[7],t.rest[8]}}};
    M f=mul(M{{p[1].p-p[0].p,p[2].p-p[0].p,p[3].p-p[0].p}},rest);
    M strain;
    for(int j=0;j<3;++j)
        strain.c[j]={(dot(f.c[0],f.c[j])-(j==0))*.5,(dot(f.c[1],f.c[j])-(j==1))*.5,(dot(f.c[2],f.c[j])-(j==2))*.5};
    double mu=1/(2*(1+t.poisson)),lambda=t.poisson/((1+t.poisson)*(1-2*t.poisson));
    double trace=strain.c[0].x+strain.c[1].y+strain.c[2].z;
    double norm=0;for(int j=0;j<3;++j) norm+=dot(strain.c[j],strain.c[j]);
    double c2=2*t.volume*(mu*norm+.5*lambda*trace*trace);
    if(c2<=0) return;
    M stress;for(int j=0;j<3;++j) stress.c[j]=strain.c[j]*(2*mu);
    stress.c[0].x+=lambda*trace;stress.c[1].y+=lambda*trace;stress.c[2].z+=lambda*trace;
    M g=gradient(mul(f,stress),rest);for(int j=0;j<3;++j) g.c[j]=g.c[j]*t.volume;
    V grad[4]={(g.c[0]+g.c[1]+g.c[2])*-1,g.c[0],g.c[1],g.c[2]};
    double weight=0;
    for(int j=0;j<4;++j) weight+=p[j].w*dot(grad[j],grad[j]);
    // Algebraically eliminate sqrt(2U) from the XPBD constraint and its gradient.
    // No dimensional epsilon: valid small tetrahedra must retain their stiffness.
    double denominator=weight+c2/(t.young*dt*dt);
    if(denominator<=0) return;
    double multiplier=-c2/denominator;
    for(int j=0;j<4;++j)
    {
        V dx=grad[j]*(p[j].w*multiplier);
        p[j].p=p[j].p+dx;p[j].v=p[j].v+dx*(1/dt);
        nodes[t.n[j]]=p[j];accumulated[t.n[j]]=accumulated[t.n[j]]+dx;
    }
}
// Medium uses float material arithmetic in a local frame; accepted positions,
// integration, contact response and correction accumulation retain double precision.
namespace FastMaterial {
struct V { float x,y,z; };
struct Node { V p; float w; V v; float damping; };
struct Tet { int n[4]; float rest[9]; float volume,young,poisson,damping; };
struct M { V c[3]; };
struct Q { float x,y,z,w; };
__device__ V operator+(V a,V b) { return {a.x+b.x,a.y+b.y,a.z+b.z}; }
__device__ V operator-(V a,V b) { return {a.x-b.x,a.y-b.y,a.z-b.z}; }
__device__ V operator*(V a,float b) { return {a.x*b,a.y*b,a.z*b}; }
__device__ float dot(V a,V b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
__device__ V cross(V a,V b) { return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x}; }
__device__ V mul(M a,V b) { return a.c[0]*b.x+a.c[1]*b.y+a.c[2]*b.z; }
__device__ M mul(M a,M b) { return {mul(a,b.c[0]),mul(a,b.c[1]),mul(a,b.c[2])}; }
__device__ M rotation(Q q)
{
    return {{{1-2*(q.y*q.y+q.z*q.z),2*(q.x*q.y+q.z*q.w),2*(q.x*q.z-q.y*q.w)},
             {2*(q.x*q.y-q.z*q.w),1-2*(q.x*q.x+q.z*q.z),2*(q.y*q.z+q.x*q.w)},
             {2*(q.x*q.z+q.y*q.w),2*(q.y*q.z-q.x*q.w),1-2*(q.x*q.x+q.y*q.y)}}};
}
__device__ Q extract(M f,Q q,int iterations)
{
    for(int k=0;k<iterations;++k)
    {
        M r=rotation(q);
        V omega=(cross(r.c[0],f.c[0])+cross(r.c[1],f.c[1])+cross(r.c[2],f.c[2]))*
            (1.f/(fabsf(dot(r.c[0],f.c[0])+dot(r.c[1],f.c[1])+dot(r.c[2],f.c[2]))+1.e-6f));
        float angle=sqrtf(dot(omega,omega));
        if(angle<1.e-6f) break;
        V a=omega*(sinf(angle*.5f)/angle), b={q.x,q.y,q.z};
        float w=cosf(angle*.5f);
        V c=a*q.w+b*w+cross(a,b);
        q={c.x,c.y,c.z,w*q.w-dot(a,b)};
        float s=rsqrtf(q.x*q.x+q.y*q.y+q.z*q.z+q.w*q.w);
        q={q.x*s,q.y*s,q.z*s,q.w*s};
    }
    return q;
}
__device__ M gradient(M tensor,M rest)
{
    return {{mul(tensor,V{rest.c[0].x,rest.c[1].x,rest.c[2].x}),
             mul(tensor,V{rest.c[0].y,rest.c[1].y,rest.c[2].y}),
             mul(tensor,V{rest.c[0].z,rest.c[1].z,rest.c[2].z})}};
}
__device__ void constrain(const Node* p,M tensor,M rest,float C,float compliance,float damping,float invDt,V* correction)
{
    M g=gradient(tensor,rest);
    V grad[4]={(g.c[0]+g.c[1]+g.c[2])*-1.f,g.c[0],g.c[1],g.c[2]};
    float weight=0,rate=0;
    for(int j=0;j<4;++j) { weight+=p[j].w*dot(grad[j],grad[j]);rate+=dot(grad[j],p[j].v); }
    float denominator=(1+damping*invDt)*weight+compliance;
    if(denominator<1.e-15f) return;
    float dl=(-C-damping*rate)/denominator;
    for(int j=0;j<4;++j) correction[j]=correction[j]+grad[j]*(dl*p[j].w);
}

__device__ Tet loadTet(const ::Tet& t)
{
    Tet result;for(int i=0;i<4;++i)result.n[i]=t.n[i];
    for(int i=0;i<9;++i)result.rest[i]=float(t.rest[i]);
    result.volume=float(t.volume);result.young=float(t.young);result.poisson=float(t.poisson);result.damping=float(t.damping);return result;
}
__device__ void loadNodes(const ::Node* nodes,const Tet& t,Node* p)
{
    const ::V origin=nodes[t.n[0]].p;
    for(int j=0;j<4;++j)
    {
        const ::Node n=nodes[t.n[j]];
        p[j]={{float(n.p.x-origin.x),float(n.p.y-origin.y),float(n.p.z-origin.z)},float(n.w),
              {float(n.v.x),float(n.v.y),float(n.v.z)},float(n.damping)};
    }
}
__device__ void solveFastTet(::Node* nodes,const ::Tet* tets,Q* rotations,int i,double stepDt,int first)
{
    float dt=float(stepDt);Tet t=loadTet(tets[i]);Node p[4];loadNodes(nodes,t,p);
    M rest={{{t.rest[0],t.rest[1],t.rest[2]},{t.rest[3],t.rest[4],t.rest[5]},{t.rest[6],t.rest[7],t.rest[8]}}};
    M edges={{p[1].p-p[0].p,p[2].p-p[0].p,p[3].p-p[0].p}};
    M f=mul(edges,rest);
    Q q=extract(f,rotations[i],first?100:4);rotations[i]=q;M r=rotation(q);
    M tensor;float c2=1.e-18f;
    for(int j=0;j<3;++j) { tensor.c[j]=f.c[j]-r.c[j];c2+=dot(tensor.c[j],tensor.c[j]); }
    float c=sqrtf(c2);for(int j=0;j<3;++j) tensor.c[j]=tensor.c[j]*(1.f/c);
    float mu=t.young/(2*(1+t.poisson));float invDt=1.f/dt;
    V correction[4]={{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
    constrain(p,tensor,rest,c,invDt*invDt/(2*mu*t.volume),0,invDt,correction);
    if(t.poisson>0)
    {
        // Preserve sub-ULP shape corrections when constructing the volume constraint.
        for(int j=0;j<3;++j) edges.c[j]=edges.c[j]+(correction[j+1]-correction[0]);
        f=mul(edges,rest);
        tensor={{cross(f.c[1],f.c[2]),cross(f.c[2],f.c[0]),cross(f.c[0],f.c[1])}};
        float lambda=t.young*t.poisson/((1+t.poisson)*(1-2*t.poisson));
        constrain(p,tensor,rest,dot(f.c[0],tensor.c[0])-1,invDt*invDt/(lambda*t.volume),0,invDt,correction);
    }
    for(int j=0;j<4;++j)
    {
        ::V dx={correction[j].x,correction[j].y,correction[j].z};
        nodes[t.n[j]].p=::operator+(nodes[t.n[j]].p,dx);
        nodes[t.n[j]].v=::operator+(nodes[t.n[j]].v,::operator*(dx,1./stepDt));
    }
}
extern "C" __global__ void solveFastTets(::Node* nodes,const ::Tet* tets,Q* rotations,int begin,int count,double stepDt,int first)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<count) solveFastTet(nodes,tets,rotations,begin+i,stepDt,first);
}
extern "C" __global__ void dampFastTets(::Node* nodes,const ::Tet* tets,int begin,int count,double stepDt)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=count) return;
    float dt=float(stepDt);Tet t=loadTet(tets[begin+i]);if(t.damping<=0) return;
    Node p[4];loadNodes(nodes,t,p);V correction[4]={};
    M rest={{{t.rest[0],t.rest[1],t.rest[2]},{t.rest[3],t.rest[4],t.rest[5]},{t.rest[6],t.rest[7],t.rest[8]}}};
    M f=mul(M{{p[1].p-p[0].p,p[2].p-p[0].p,p[3].p-p[0].p}},rest);
    float mu=t.young/(2*(1+t.poisson)),lambda=t.young*t.poisson/((1+t.poisson)*(1-2*t.poisson));
    const float a=0.57735026918962576451,b=0.70710678118654752440,c=0.40824829046386301637;
    M modes[6]={{{{a,0,0},{0,a,0},{0,0,a}}},{{{b,0,0},{0,-b,0},{0,0,0}}},
        {{{c,0,0},{0,c,0},{0,0,-2*c}}},{{{0,b,0},{b,0,0},{0,0,0}}},
        {{{0,0,b},{0,0,0},{b,0,0}}},{{{0,0,0},{0,0,b},{0,b,0}}}};
    for(int mode=0;mode<6;++mode)
    {
        M g=gradient(mul(f,modes[mode]),rest);
        V grad[4]={(g.c[0]+g.c[1]+g.c[2])*-1,g.c[0],g.c[1],g.c[2]};
        float weight=0,rate=0;
        for(int j=0;j<4;++j) { weight+=p[j].w*dot(grad[j],grad[j]);rate+=dot(grad[j],p[j].v); }
        float viscosity=t.damping*t.volume*(2*mu+(mode==0?3*lambda:0))*dt;
        float impulse=-viscosity*rate/(1+viscosity*weight);
        for(int j=0;j<4;++j) { V dv=grad[j]*(p[j].w*impulse);p[j].v=p[j].v+dv;correction[j]=correction[j]+dv; }
    }
    for(int j=0;j<4;++j) nodes[t.n[j]].v=::operator+(nodes[t.n[j]].v,::V{correction[j].x,correction[j].y,correction[j].z});
}
}
extern "C" __global__ void initializeCopies(const Node* nodes,Node* copies,const int* active,int count,int lanes)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=count*lanes) return;
    int node=i%count,used=0;for(int lane=0;lane<lanes;++lane) used+=active[lane*count+node];
    // Split physical mass between active copies so their averaged impulses conserve momentum.
    Node copy=nodes[node];copy.w*=used;copies[i]=copy;
}
extern "C" __global__ void averageCopies(Node* nodes,const Node* copies,const int* active,V* accumulated,int count,int lanes,double dt,int positions)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=count||nodes[i].w==0) return;
    V delta={0,0,0};int used=0;
    for(int lane=0;lane<lanes;++lane) if(active[lane*count+i]) { delta=delta+(copies[lane*count+i].v-nodes[i].v);++used; }
    if(used) { delta=delta*(1.f/used);nodes[i].v=nodes[i].v+delta;if(positions) { V dx=delta*dt;nodes[i].p=nodes[i].p+dx;accumulated[i]=accumulated[i]+dx; } }
}
extern "C" __global__ void settle(Node* nodes,const V* accumulated,int count,double dt,double speed,double damping,int begin,int end,int* failure)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=count) return;
    Node p=nodes[i];
    if(!isfinite(p.p.x)||!isfinite(p.p.y)||!isfinite(p.p.z)||!isfinite(p.v.x)||!isfinite(p.v.y)||!isfinite(p.v.z)) atomicExch(failure,1);
    if(p.w>0)
    {
        V average=accumulated[i]*(1./dt);
        // PhysX's default per-iteration force mode disables final velocity averaging.
        if(i>=begin&&i<end&&dot(average,average)<speed*speed) p.v=p.v*fmax(0.,1-damping*dt);
    }
    nodes[i]=p;
}
struct Target { int node;V position; };
extern "C" __global__ void startDrive(Node* nodes,const Target* targets,V* start,int count,double dt)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=count) return;
    Target t=targets[i];start[i]=nodes[t.node].p;nodes[t.node].w=0;nodes[t.node].v=(t.position-start[i])*(1./dt);
}
extern "C" __global__ void drive(Node* nodes,const Target* targets,const V* start,V* accumulated,int count,double fraction,double increment)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=count) return;
    Target t=targets[i];V travel=t.position-start[i];
    nodes[t.node].p=start[i]+travel*fraction;
    accumulated[t.node]=accumulated[t.node]+travel*increment;
}
extern "C" __global__ void releaseDrive(Node* nodes,const Target* targets,const double* masses,int count)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=count) return;
    int node=targets[i].node;nodes[node].w=masses[node];
}

struct Map { int nodes[4];double weights[4];V offset;int body; };
struct Tri { int v[3]; };
struct Bvh { int left,right,parent,triangle,body; };
struct Box { V low,high; };
__device__ V minimum(V a,V b) { return {fmin(a.x,b.x),fmin(a.y,b.y),fmin(a.z,b.z)}; }
__device__ V maximum(V a,V b) { return {fmax(a.x,b.x),fmax(a.y,b.y),fmax(a.z,b.z)}; }
extern "C" __global__ void skin(const Node* nodes,const Map* maps,V* positions,int count)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=count) return;
    Map m=maps[i];V p=m.offset;
    for(int j=0;j<4;++j) if(m.nodes[j]>=0) p=p+nodes[m.nodes[j]].p*m.weights[j];
    positions[i]=p;
}
extern "C" __global__ void refit(const V* positions,const Tri* tris,const Bvh* tree,const int* leaves,Box* boxes,int* ready,int count)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=count) return;
    int index=leaves[i];Tri t=tris[tree[index].triangle];
    V a=positions[t.v[0]],b=positions[t.v[1]],c=positions[t.v[2]];
    boxes[index]={minimum(a,minimum(b,c)),maximum(a,maximum(b,c))};
    __threadfence();
    int parent=tree[index].parent;
    while(parent>=0)
    {
        if(atomicAdd(ready+parent,1)==0) return;
        Box a=boxes[tree[parent].left],b=boxes[tree[parent].right];
        boxes[parent]={minimum(a.low,b.low),maximum(a.high,b.high)};
        __threadfence();parent=tree[parent].parent;
    }
}
__device__ V closestWeights(V p,V a,V b,V c)
{
    V ab=b-a,ac=c-a,ap=p-a;double d1=dot(ab,ap),d2=dot(ac,ap);
    if(d1<=0&&d2<=0) return {1,0,0};
    V bp=p-b;double d3=dot(ab,bp),d4=dot(ac,bp);
    if(d3>=0&&d4<=d3) return {0,1,0};
    double vc=d1*d4-d3*d2;
    if(vc<=0&&d1>=0&&d3<=0) { double v=d1/(d1-d3);return {1-v,v,0}; }
    V cp=p-c;double d5=dot(ab,cp),d6=dot(ac,cp);
    if(d6>=0&&d5<=d6) return {0,0,1};
    double vb=d5*d2-d1*d6;
    if(vb<=0&&d2>=0&&d6<=0) { double w=d2/(d2-d6);return {1-w,0,w}; }
    double va=d3*d6-d5*d4;
    if(va<=0&&(d4-d3)>=0&&(d5-d6)>=0) { double w=(d4-d3)/((d4-d3)+(d5-d6));return {0,1-w,w}; }
    double inv=1./(va+vb+vc);return {1-(vb+vc)*inv,vb*inv,vc*inv};
}
__device__ void add(V* p,V v) { atomicAdd(&p->x,v.x);atomicAdd(&p->y,v.y);atomicAdd(&p->z,v.z); }
struct Contact { int vertex,triangle;V bary,normal;int triangleA;V baryA;Map simA,simB;double bias;int prepared;double impulse;int state;V tangentImpulse;double friction; };
static_assert(sizeof(Contact)==312,"Contact allocation size");
__device__ V mappedPosition(const Node* nodes,Map m)
{
    V p=m.offset;for(int j=0;j<4;++j) if(m.nodes[j]>=0) p=p+nodes[m.nodes[j]].p*m.weights[j];return p;
}
__device__ int contactSupport(Contact contact,const Map* maps,const Tri* tris,int* ids,double* weights)
{
    if(contact.prepared==1)
    {
        Map support[2]={contact.simA,contact.simB};int used=0;
        for(int j=0;j<2;++j) for(int k=0;k<4;++k)
        {
            int id=support[j].nodes[k];if(id<0) continue;double w=support[j].weights[k]*(j?-1.:1.);
            int slot=0;while(slot<used&&ids[slot]!=id) ++slot;
            if(slot==used) { ids[used]=id;weights[used++]=0; }weights[slot]+=w;
        }
        return used;
    }
    Tri t=tris[contact.triangle];int used=0;
    int a[3]={contact.vertex,contact.vertex,contact.vertex};V baryA={1,0,0};
    if(contact.triangleA>=0) { Tri ta=tris[contact.triangleA];for(int j=0;j<3;++j) a[j]=ta.v[j];baryA=contact.baryA; }
    Map support[6]={maps[a[0]],maps[a[1]],maps[a[2]],maps[t.v[0]],maps[t.v[1]],maps[t.v[2]]};
    double factor[6]={baryA.x,baryA.y,baryA.z,-contact.bary.x,-contact.bary.y,-contact.bary.z};
    for(int j=0;j<6;++j) for(int k=0;k<4;++k)
    {
        int id=support[j].nodes[k];if(id<0) continue;
        double w=support[j].weights[k]*factor[j];int slot=0;while(slot<used&&ids[slot]!=id) ++slot;
        if(slot==used) { ids[used]=id;weights[used++]=0; }weights[slot]+=w;
    }
    return used;
}
__device__ V mappedDelta(const V* accumulated,Map m)
{
    V delta={0,0,0};for(int j=0;j<4;++j) if(m.nodes[j]>=0) delta=delta+accumulated[m.nodes[j]]*m.weights[j];return delta;
}
__device__ V sourceTriangleDelta(Contact contact,const V* accumulated,const Map* maps,const Tri* tris)
{
    Tri t=tris[contact.triangleA];
    return mappedDelta(accumulated,maps[t.v[0]])*contact.baryA.x+
        mappedDelta(accumulated,maps[t.v[1]])*contact.baryA.y+mappedDelta(accumulated,maps[t.v[2]])*contact.baryA.z;
}
__device__ double separation(Contact contact,const Node* nodes,const V* accumulated,const Map* maps,const Tri* tris)
{
    // PhysX tracks contact motion using accumulated displacement, avoiding the
    // loss of small corrections when subtracting rounded world positions.
    if(contact.prepared) return dot((contact.prepared==2?sourceTriangleDelta(contact,accumulated,maps,tris):mappedDelta(accumulated,contact.simA))-
        mappedDelta(accumulated,contact.simB),contact.normal)+contact.bias;
    Tri t=tris[contact.triangle];
    V point;
    if(contact.triangleA<0) point=mappedPosition(nodes,maps[contact.vertex]);
    else { Tri a=tris[contact.triangleA];point=mappedPosition(nodes,maps[a.v[0]])*contact.baryA.x+mappedPosition(nodes,maps[a.v[1]])*contact.baryA.y+mappedPosition(nodes,maps[a.v[2]])*contact.baryA.z; }
    V delta=point-mappedPosition(nodes,maps[t.v[0]])*contact.bary.x-
        mappedPosition(nodes,maps[t.v[1]])*contact.bary.y-mappedPosition(nodes,maps[t.v[2]])*contact.bary.z;
    return dot(delta,contact.normal);
}
__device__ bool pairDisabled(int a,int b,const unsigned long long* pairs,int count)
{
    if(a<0||b<0) return false; // A mixed BVH branch must still be traversed.
    if(a==b) return true;
    if(a>b) { int t=a;a=b;b=t; }
    unsigned long long key=(static_cast<unsigned long long>(a)<<32)|static_cast<unsigned int>(b);
    int lo=0,hi=count;
    while(lo<hi) { int mid=lo+(hi-lo)/2;if(pairs[mid]<key) lo=mid+1;else hi=mid; }
    return lo<count&&pairs[lo]==key;
}
extern "C" __global__ void discoverContacts(const Node* nodes,const Map* maps,const V* positions,const Tri* tris,
    const Bvh* tree,const Box* boxes,Contact* cache,int* contactCount,int capacity,int count,double margin,int* failure,const int* incident,
    const unsigned long long* disabledPairs,int disabledPairCount)
{
    int vertex=blockIdx.x*blockDim.x+threadIdx.x;if(vertex>=count) return;
    Map point=maps[vertex];V p=positions[vertex];
    bool movingPoint=false;for(int j=0;j<4;++j) movingPoint|=point.nodes[j]>=0;
    if(incident[vertex]&&!movingPoint) return;
    // Static vertices are queried too: narrow pegs must also collide against dynamic faces.
    int stack[64],top=0;stack[top++]=0;
    while(top)
    {
        int index=stack[--top];Bvh node=tree[index];Box box=boxes[index];
        if(pairDisabled(node.body,point.body,disabledPairs,disabledPairCount) || p.x<box.low.x-margin||p.y<box.low.y-margin||p.z<box.low.z-margin||
           p.x>box.high.x+margin||p.y>box.high.y+margin||p.z>box.high.z+margin) continue;
        if(node.triangle<0)
        {
            if(top>61) { atomicExch(failure,2);return; }
            stack[top++]=node.left;stack[top++]=node.right;continue;
        }
        Tri t=tris[node.triangle];
        // Triangle witnesses alone do not constrain every corner of a broad contact patch.
        if(incident[vertex])
        {
            bool moving=false;for(int k=0;k<3;++k) for(int j=0;j<4;++j) moving|=maps[t.v[k]].nodes[j]>=0;
            if(moving) continue;
        }
        V a=positions[t.v[0]],b=positions[t.v[1]],c=positions[t.v[2]];
        V normal=cross(b-a,c-a);double length2=dot(normal,normal);if(length2<1.e-20) continue;
        normal=normal*rsqrt(length2);
        if(incident[vertex])
        {
            // Supplemental face contacts must not extend the face past a hole's rim.
            // Edge contacts are already handled by the triangle-pair narrowphase.
            V u=b-a,v=c-a,q=p-a;double uu=dot(u,u),uv=dot(u,v),vv=dot(v,v);
            double denominator=uu*vv-uv*uv;
            if(denominator<=0) continue;
            double w1=(dot(q,u)*vv-dot(q,v)*uv)/denominator;
            double w2=(dot(q,v)*uu-dot(q,u)*uv)/denominator;
            if(w1<-1.e-6||w2<-1.e-6||w1+w2>1.+1.e-6) continue;
        }
        V bary=closestWeights(p,a,b,c);
        V delta=p-(a*bary.x+b*bary.y+c*bary.z);double d2=dot(delta,delta);if(d2>margin*margin) continue;
        double signedDistance=sqrt(d2)*(dot(delta,normal)>=0?1.:-1.);
        bool dynamic=false;
        for(int j=0;j<4;++j) if(point.nodes[j]>=0 && nodes[point.nodes[j]].w>0) dynamic=true;
        for(int k=0;k<3;++k) for(int j=0;j<4;++j) if(maps[t.v[k]].nodes[j]>=0 && nodes[maps[t.v[k]].nodes[j]].w>0) dynamic=true;
        if(!dynamic) continue;
        // Behind a one-sided face, edge distance must not introduce a tangential attraction.
        // Interior witnesses use the face normal; rounded barycentrics otherwise create false friction.
        if(d2>1.e-16 && signedDistance>0 && fmin(bary.x,fmin(bary.y,bary.z))<=1.e-6)
            normal=delta*(1./signedDistance);
        int slot=atomicAdd(contactCount,1);
        if(slot>=capacity) { atomicExch(failure,4);continue; }
        cache[slot]={vertex,node.triangle,bary,normal,-1,{1,0,0}};
        bool material=incident[vertex]!=0;for(int j=0;j<4;++j) if(point.weights[j]<0) material=false;
        if(material) cache[slot].state=16; // Retain a non-extrapolated vertex's material embedding.
    }
}
__device__ double clamp01(double x) { return fmin(1.,fmax(0.,x)); }
__device__ V corner(int i) { return {i==0?1.:0.,i==1?1.:0.,i==2?1.:0.}; }
// Adapted from PhysX distanceSegmentSegment.cuh; guard division for parallel edges.
__device__ void edgeParameters(V p,V d,V q,V e,double& s,double& t)
{
    V r=p-q;double a=dot(d,d),b=dot(d,e),c=dot(d,r),f=dot(e,r),ee=dot(e,e);
    double denominator=a*ee-b*b;
    double seed=denominator>1.e-6?clamp01((b*f-c*ee)/denominator):.5;
    t=clamp01((b*seed+f)*(ee>1.e-6?1./ee:0.));
    s=clamp01((b*t-c)*(a>1.e-6?1./a:0.));
}
__device__ V trianglePoint(const V* a,V w) { return a[0]*w.x+a[1]*w.y+a[2]*w.z; }
__device__ void closestCandidate(const V* a,const V* b,V wa,V wb,double& distance,V& bestA,V& bestB)
{
    V delta=trianglePoint(a,wa)-trianglePoint(b,wb);double d=dot(delta,delta);
    if(d<distance) { distance=d;bestA=wa;bestB=wb; }
}
__device__ void edgeFaceCandidate(const V* a,const V* b,int i,bool reverse,double& distance,V& bestA,V& bestB)
{
    int j=(i+1)%3;V normal=cross(b[1]-b[0],b[2]-b[0]),edge=a[j]-a[i];
    double denominator=dot(edge,normal);if(fabs(denominator)<1.e-16) return;
    double t=dot(b[0]-a[i],normal)/denominator;if(t<0||t>1) return;
    V p=a[i]+edge*t,wb=closestWeights(p,b[0],b[1],b[2]);
    V d=p-trianglePoint(b,wb);if(dot(d,d)>1.e-10) return;
    V wa=corner(i)*(1-t)+corner(j)*t;
    if(reverse) closestCandidate(b,a,wb,wa,distance,bestA,bestB);
    else closestCandidate(a,b,wa,wb,distance,bestA,bestB);
}
__device__ bool movingTriangle(Tri t,const Map* maps,const Node* nodes)
{
    for(int i=0;i<3;++i) for(int j=0;j<4;++j)
        if(maps[t.v[i]].nodes[j]>=0&&nodes[maps[t.v[i]].nodes[j]].w>0) return true;
    return false;
}
// Candidate metric adapted from PhysX triangletriangle.cuh (Apache-2.0).
// Keep signed SAT separation separate from the narrowphase witness distance.
__device__ double triangleSeparation(const V* a,const V* b,V an,V bn,double tolerance)
{
    double sep=fmin(dot(an,b[0]-a[0]),fmin(dot(an,b[1]-a[0]),dot(an,b[2]-a[0])));
    sep=fmax(sep,fmin(dot(bn,a[0]-b[0]),fmin(dot(bn,a[1]-b[0]),dot(bn,a[2]-b[0]))));
    if(sep>tolerance) return sep;
    V ae[3]={a[1]-a[0],a[2]-a[0],a[2]-a[1]};
    V be[3]={b[1]-b[0],b[2]-b[0],b[2]-b[1]};
    for(int i=0;i<3;++i) for(int j=0;j<3;++j)
    {
        V axis=cross(ae[i],be[j]);double n2=dot(axis,axis);
        if(n2<=1.e-5) continue;
        axis=axis*rsqrt(n2);
        double alo=dot(axis,a[0]),ahi=alo,blo=dot(axis,b[0]),bhi=blo;
        for(int k=1;k<3;++k)
        {
            double ap=dot(axis,a[k]),bp=dot(axis,b[k]);
            alo=fmin(alo,ap);ahi=fmax(ahi,ap);blo=fmin(blo,bp);bhi=fmax(bhi,bp);
        }
        sep=fmax(sep,fmax(blo-ahi,alo-bhi));
        if(sep>tolerance) return sep;
    }
    return sep;
}
extern "C" __global__ void discoverTriangleContacts(const Node* nodes,const Map* maps,const V* positions,const Tri* tris,
    const Bvh* tree,const Box* boxes,Contact* cache,int* contactCount,int capacity,int count,double margin,int* failure,
    const int* volumeBodies,int volumeCount,const unsigned long long* disabledPairs,int disabledPairCount)
{
    int indexA=blockIdx.x*blockDim.x+threadIdx.x;if(indexA>=count) return;
    Tri ta=tris[indexA];if(!movingTriangle(ta,maps,nodes)) return;
    V a[3]={positions[ta.v[0]],positions[ta.v[1]],positions[ta.v[2]]};
    V an=cross(a[1]-a[0],a[2]-a[0]);if(dot(an,an)<1.e-20) return;
    an=an*rsqrt(dot(an,an));
    Box bound={minimum(a[0],minimum(a[1],a[2])),maximum(a[0],maximum(a[1],a[2]))};
    int body=maps[ta.v[0]].body,stack[64],top=0;stack[top++]=0;
    Contact best[4];double distances[4]={margin,margin,margin,margin};int found=0;
    while(top)
    {
        int index=stack[--top];Bvh node=tree[index];Box box=boxes[index];
        if(pairDisabled(node.body,body,disabledPairs,disabledPairCount)||bound.high.x<box.low.x-margin||bound.low.x>box.high.x+margin||
            bound.high.y<box.low.y-margin||bound.low.y>box.high.y+margin||bound.high.z<box.low.z-margin||bound.low.z>box.high.z+margin) continue;
        if(node.triangle<0)
        {
            if(top>61) { atomicExch(failure,2);return; }
            stack[top++]=node.left;stack[top++]=node.right;continue;
        }
        Tri tb=tris[node.triangle];bool dynamicB=movingTriangle(tb,maps,nodes);
        bool volumeA=false,volumeB=false;
        for(int j=0;j<volumeCount;++j) { volumeA|=volumeBodies[j]==body;volumeB|=volumeBodies[j]==maps[tb.v[0]].body; }
        if(volumeA&&volumeB) continue; // separate signed vertex/surface path below
        if(dynamicB&&node.triangle<indexA) continue;
        V b[3]={positions[tb.v[0]],positions[tb.v[1]],positions[tb.v[2]]};
        V bn=cross(b[1]-b[0],b[2]-b[0]);double n2=dot(bn,bn);if(n2<1.e-20) continue;
        if(!dynamicB&&dot((a[0]+a[1]+a[2])*(1./3.)-b[0],bn)<0) continue;
        bn=bn*rsqrt(n2);
        double separation=triangleSeparation(a,b,an,bn,margin);
        if(separation>=margin) continue;
        V wa=corner(0),wb=corner(0),normal=bn;
        int slot=found<4?found++:0;
        if(found==4) { for(int j=1;j<4;++j) if(distances[j]>distances[slot]) slot=j;if(separation>=distances[slot]) continue; }
        best[slot]={ta.v[0],node.triangle,wb,normal,indexA,wa};distances[slot]=separation;
    }
    for(int i=0;i<found;++i)
    {
        int slot=atomicAdd(contactCount,1);if(slot>=capacity) { atomicExch(failure,4);continue; }cache[slot]=best[i];
    }
}
// Adapted from PhysX SbSbMeshTreeTraverser exact-surface path: one nearest
// signed contact per source collision vertex and target deformable volume.
extern "C" __global__ void discoverVolumeContacts(const Map* maps,const V* positions,const Tri* tris,const Bvh* tree,
    const Box* boxes,const int* incident,const int* volumeBodies,int volumeCount,int vertexCount,
    Contact* cache,int* contactCount,int capacity,double margin,int* failure,const unsigned long long* disabledPairs,int disabledPairCount)
{
    int job=blockIdx.x*blockDim.x+threadIdx.x;
    if(job>=vertexCount*volumeCount) return;
    int vertex=job%vertexCount,target=volumeBodies[job/vertexCount],source=maps[vertex].body;
    if(!incident[vertex]||pairDisabled(source,target,disabledPairs,disabledPairCount)) return;
    bool sourceVolume=false;for(int j=0;j<volumeCount;++j) sourceVolume|=volumeBodies[j]==source;
    if(!sourceVolume) return;
    V p=positions[vertex];int stack[64],top=0;stack[top++]=0;
    double best=margin*margin;Contact result={};result.triangle=-1;
    while(top)
    {
        int index=stack[--top];Bvh node=tree[index];Box box=boxes[index];
        if((node.body>=0&&node.body!=target)||p.x<box.low.x-margin||p.x>box.high.x+margin||
            p.y<box.low.y-margin||p.y>box.high.y+margin||p.z<box.low.z-margin||p.z>box.high.z+margin) continue;
        if(node.triangle<0)
        {
            if(top>61) { atomicExch(failure,2);return; }
            stack[top++]=node.left;stack[top++]=node.right;continue;
        }
        Tri t=tris[node.triangle];V a=positions[t.v[0]],b=positions[t.v[1]],c=positions[t.v[2]];
        V n=cross(b-a,c-a);double area=dot(n,n);if(area<=0||!isfinite(area)) continue;
        V bary=closestWeights(p,a,b,c),delta=p-(a*bary.x+b*bary.y+c*bary.z);
        double d2=dot(delta,delta);if(d2>=best) continue;
        double length=sqrt(d2),side=dot(delta,n)<0?-1.:1.;
        V normal=length>1.e-10?delta*(side/length):n*rsqrt(area);
        result={vertex,node.triangle,bary,normal,-1,{1,0,0}};result.state=16;best=d2;
    }
    if(result.triangle>=0)
    {
        int slot=atomicAdd(contactCount,1);if(slot>=capacity) { atomicExch(failure,4);return; }cache[slot]=result;
    }
}
// One warp per pair, preserving the cooperative layout of PhysX's GJK/EPA.
extern "C" __global__ void refineTriangleContacts(const Node* nodes,const Map* maps,const V* positions,const Tri* tris,
    Contact* cache,const int* contactCount,int capacity,double margin,int* failure)
{
    struct Scratch { physx::PxVec3 a[4],b[4];schlock::GjkCachedData cached;schlock::GjkOutput output;physx::squawk::EpaScratch epa; };
    __shared__ __align__(16) unsigned char storage[4*sizeof(Scratch)];
    Scratch& s=reinterpret_cast<Scratch*>(storage)[threadIdx.y];
    int count=min(*contactCount,capacity);
    for(int i=blockIdx.x*blockDim.y+threadIdx.y;i<count;i+=gridDim.x*blockDim.y)
    {
        Contact c=cache[i];if(c.triangleA<0) continue;
        Tri ta=tris[c.triangleA],tb=tris[c.triangle];
        if(threadIdx.x<3)
        {
            V a=positions[ta.v[threadIdx.x]],b=positions[tb.v[threadIdx.x]];
            s.a[threadIdx.x]={float(a.x),float(a.y),float(a.z)};s.b[threadIdx.x]={float(b.x),float(b.y),float(b.z)};
        }
        if(threadIdx.x==0) { s.cached.size=0;s.a[3]=s.b[3]=physx::PxVec3(0); }
        __syncwarp();
        auto dir=(s.a[0]+s.a[1]+s.a[2]-s.b[0]-s.b[1]-s.b[2])*(1./3.);
        if(dir.magnitudeSquared()==0) dir={1,0,0};
        const double ratio=1.-.000225;
        auto result=physx::squawk::gjk(s.a,3,s.b,3,dir,margin,ratio,s.output,s.cached);
        __syncwarp();
        if(result!=schlock::GjkResult::eDISTANT && (result!=schlock::GjkResult::eCLOSE||s.output.closestPointDir.dot(s.output.direction)<.999))
        {
            result=physx::squawk::epa(s.epa,s.a,3,s.b,3,&s.cached,ratio,1.e-6,s.output);
            __syncwarp();
            if(s.output.degenerate)
            {
                if(threadIdx.x==0) s.cached.size=0;__syncwarp();
                physx::squawk::gjk(s.a,3,s.b,3,physx::PxVec3(0,1,0),margin,ratio,s.output,s.cached);
                __syncwarp();
                result=physx::squawk::epa(s.epa,s.a,3,s.b,3,&s.cached,ratio,1.e-6,s.output);
            }
        }
        __syncwarp();
        if(threadIdx.x==0)
        {
            if(result==schlock::GjkResult::eDISTANT) c.state=4;
            else
            {
                auto n=s.output.direction.getNormalized();V normal={-n.x,-n.y,-n.z};
                auto pa=s.output.closestPointA,pb=s.output.closestPointB;
                if(!isfinite(n.x)||!isfinite(n.y)||!isfinite(n.z)||!isfinite(pa.x)||!isfinite(pa.y)||!isfinite(pa.z)||
                    !isfinite(pb.x)||!isfinite(pb.y)||!isfinite(pb.z)) { atomicExch(failure,1);c.state=4; }
                else if(!movingTriangle(tb,maps,nodes)&&dot(normal,cross(positions[tb.v[1]]-positions[tb.v[0]],positions[tb.v[2]]-positions[tb.v[0]]))<0) c.state=4;
                else
                {
                    c.normal=normal;
                    // GJK's second witness is reconstructed from its distance bound;
                    // retain it exactly for the PhysX contact-to-tet remapping.
                    c.simA.offset={pa.x,pa.y,pa.z};c.simB.offset={pb.x,pb.y,pb.z};c.state|=8;
                    c.baryA=closestWeights({pa.x,pa.y,pa.z},positions[ta.v[0]],positions[ta.v[1]],positions[ta.v[2]]);
                    c.bary=closestWeights({pb.x,pb.y,pb.z},positions[tb.v[0]],positions[tb.v[1]],positions[tb.v[2]]);
                }
            }
            cache[i]=c;
        }
        __syncwarp();
    }
}
__device__ bool remapContactPoint(V point,int body,const Node* nodes,const Tet* tets,const int* tetBodies,int tetCount,
    int triangle,const int* offsets,const int* candidates,Map& result)
{
    double best=3.402823466e38;bool found=false;
    const int faces[4][3]={{1,2,3},{0,3,2},{0,1,3},{0,2,1}};
    int begin=triangle<0?0:offsets[triangle],end=triangle<0?tetCount:offsets[triangle+1];
    for(int candidate=begin;candidate<end;++candidate)
    {
        int i=triangle<0?candidate:candidates[candidate];
        if(tetBodies[i]!=body) continue;
        Tet tet=tets[i];V x[4];for(int j=0;j<4;++j) x[j]=nodes[tet.n[j]].p;
        V u=x[1]-x[0],v=x[2]-x[0],w=x[3]-x[0],q=point-x[0];double det=dot(u,cross(v,w));
        if(fabs(det)>1.e-16)
        {
            double b[4];b[1]=dot(q,cross(v,w))/det;b[2]=dot(u,cross(q,w))/det;b[3]=dot(u,cross(v,q))/det;b[0]=1-b[1]-b[2]-b[3];
            if(b[0]>=0&&b[1]>=0&&b[2]>=0&&b[3]>=0)
            {
                result={{tet.n[0],tet.n[1],tet.n[2],tet.n[3]},{b[0],b[1],b[2],b[3]},{0,0,0},body};return true;
            }
        }
        for(int face=0;face<4;++face)
        {
            int a=faces[face][0],b=faces[face][1],c=faces[face][2];V n=cross(x[b]-x[a],x[c]-x[a]);if(dot(n,n)<1.e-20) continue;
            V bc=closestWeights(point,x[a],x[b],x[c]);V delta=point-(x[a]*bc.x+x[b]*bc.y+x[c]*bc.z);double d=dot(delta,delta);
            if(d<best)
            {
                best=d;found=true;result={{tet.n[0],tet.n[1],tet.n[2],tet.n[3]},{0,0,0,0},{0,0,0},body};
                result.weights[a]=bc.x;result.weights[b]=bc.y;result.weights[c]=bc.z;
            }
        }
    }
    return found;
}
extern "C" __global__ void prepareContacts(const Node* nodes,const Tet* tets,const int* tetBodies,int tetCount,
    const int* offsets,const int* candidates,const V* accumulated,const Map* maps,const Tri* tris,Contact* cache,const int* contactCount,int capacity,
    const double* friction,int frictionCount)
{
    int count=min(*contactCount,capacity);
    for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<count;i+=blockDim.x*gridDim.x)
    {
        Contact c=cache[i];if(c.state&4) continue;Tri tb=tris[c.triangle];
        V pb=mappedPosition(nodes,maps[tb.v[0]])*c.bary.x+mappedPosition(nodes,maps[tb.v[1]])*c.bary.y+mappedPosition(nodes,maps[tb.v[2]])*c.bary.z;
        V pa;bool dynamicA=false,dynamicB=false;int bodyA=maps[c.vertex].body,bodyB=maps[tb.v[0]].body;
        if(c.triangleA<0) { pa=mappedPosition(nodes,maps[c.vertex]);for(int j=0;j<4;++j) dynamicA|=maps[c.vertex].nodes[j]>=0; }
        else
        {
            Tri ta=tris[c.triangleA];pa=mappedPosition(nodes,maps[ta.v[0]])*c.baryA.x+mappedPosition(nodes,maps[ta.v[1]])*c.baryA.y+mappedPosition(nodes,maps[ta.v[2]])*c.baryA.z;
            for(int j=0;j<3;++j) for(int k=0;k<4;++k) dynamicA|=maps[ta.v[j]].nodes[k]>=0;
        }
        for(int j=0;j<3;++j) for(int k=0;k<4;++k) dynamicB|=maps[tb.v[j]].nodes[k]>=0;
        if(c.state&8) { pa=c.simA.offset;pb=c.simB.offset; }
        c.simA={{-1,-1,-1,-1},{0,0,0,0},pa,bodyA};c.simB={{-1,-1,-1,-1},{0,0,0,0},pb,bodyB};
        // PhysX's rigid-volume writer remaps the rigid witness point onto the volume.
        V remapA=dynamicA&&!dynamicB?pb:pa;
        bool materialSupport=dynamicA&&!dynamicB&&c.triangleA>=0;
        if(materialSupport)
        {
            Tri ta=tris[c.triangleA];
            for(int j=0;j<3;++j) for(int k=0;k<4;++k) if(maps[ta.v[j]].weights[k]<0) materialSupport=false;
        }
        bool a=!dynamicA||materialSupport;
        if(c.state&16) { c.simA=maps[c.vertex];a=true; }
        else if(dynamicA&&!materialSupport) a=remapContactPoint(remapA,bodyA,nodes,tets,tetBodies,tetCount,c.triangleA,offsets,candidates,c.simA);
        bool b=!dynamicB||remapContactPoint(pb,bodyB,nodes,tets,tetBodies,tetCount,c.triangle,offsets,candidates,c.simB);
        // Keep the force on the collision triangle's material when nearby elements fold over it.
        c.prepared=a&&b?(materialSupport?2:1):0;
        c.impulse=0;c.state=(!dynamicA||!dynamicB)?1:0;
        c.tangentImpulse={0,0,0};
        c.friction=bodyA>=0&&bodyA<frictionCount&&bodyB>=0&&bodyB<frictionCount?
            fmin(10.,friction[bodyA]*friction[bodyB]):0.;
        c.bias=dot(pa-pb,c.normal)-dot((materialSupport?sourceTriangleDelta(c,accumulated,maps,tris):mappedDelta(accumulated,c.simA))-
            mappedDelta(accumulated,c.simB),c.normal);
        cache[i]=c;
    }
}
extern "C" __global__ void countContacts(const Node* nodes,const V* accumulated,const Map* maps,const Tri* tris,Contact* cache,
    const int* contactCount,int capacity,int* references,double clearance,int velocityOnly,double dt)
{
    int count=min(*contactCount,capacity);
    for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<count;i+=blockDim.x*gridDim.x)
    {
        Contact c=cache[i];if(c.state&4) continue;double C=separation(c,nodes,accumulated,maps,tris)-clearance;
        int ids[24];double weights[24];int used=contactSupport(c,maps,tris,ids,weights);
        if(velocityOnly)
        {
            V velocity={0,0,0};for(int j=0;j<used;++j) velocity=velocity+nodes[ids[j]].v*weights[j];
            C=fmax(0.,C)/dt+dot(velocity,c.normal);
        }
        if(C>=0&&!(c.state&2)&&c.impulse<=0) continue;
        if(c.state&1) cache[i].state|=2;
        // PhysX bumpDbRefCountTet ignores barycentric contributions <= 1e-3.
        for(int j=0;j<used;++j) if(fabs(weights[j])>1.e-3 && nodes[ids[j]].w>0) atomicAdd(references+ids[j],1);
    }
}
extern "C" __global__ void solveContacts(const Node* nodes,const V* accumulated,const Map* maps,const Tri* tris,Contact* cache,
    const int* contactCount,int capacity,const int* references,V* corrections,double clearance,int velocityOnly,double dt)
{
    int count=min(*contactCount,capacity);
    for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<count;i+=blockDim.x*gridDim.x)
    {
        Contact c=cache[i];if(c.state&4) continue;double C=separation(c,nodes,accumulated,maps,tris)-clearance;
        int ids[24];double weights[24];int used=contactSupport(c,maps,tris,ids,weights);double denominator=0;
        if(velocityOnly)
        {
            V velocity={0,0,0};for(int j=0;j<used;++j) velocity=velocity+nodes[ids[j]].v*weights[j];
            // Match PhysX: remove only velocity that closes the remaining gap within dt.
            C=fmax(0.,C)/dt+dot(velocity,c.normal);
        }
        if(C>=0&&!(c.state&2)&&c.impulse<=0) continue;
        for(int j=0;j<used;++j) denominator+=weights[j]*weights[j]*nodes[ids[j]].w*references[ids[j]];
        if(denominator<1.e-15) continue;
        double velocityError=velocityOnly?C:C/dt;
        double impulse=fmax(-velocityError/denominator,-c.impulse);
        cache[i].impulse=c.impulse+impulse;
        V deltaImpulse=c.normal*impulse;
        if(c.friction>0)
        {
            // PhysX TGS: oppose contact drift, with one Coulomb budget across all iterations.
            V relative={0,0,0};
            for(int j=0;j<used;++j)
                relative=relative+(velocityOnly?nodes[ids[j]].v:accumulated[ids[j]]*(1./dt))*weights[j];
            V tangent=relative-c.normal*dot(relative,c.normal);
            V target=c.tangentImpulse-tangent*(1./denominator);
            double limit=c.friction*cache[i].impulse,length=sqrt(dot(target,target));
            if(length>limit) target=target*(limit/length);
            deltaImpulse=deltaImpulse+target-c.tangentImpulse;
            cache[i].tangentImpulse=target;
        }
        V correction=deltaImpulse*(velocityOnly?1.:dt);
        // Refcount-inflated scatter divided by max(refcount,1) at application:
        // positive counts cancel, but an unreferenced node receives no impulse.
        for(int j=0;j<used;++j) if(references[ids[j]]>0) add(corrections+ids[j],correction*(weights[j]*nodes[ids[j]].w));
    }
}
extern "C" __global__ void applyContactVelocities(Node* nodes,const V* corrections,int count)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<count) nodes[i].v=nodes[i].v+corrections[i];
}
extern "C" __global__ void applyContacts(Node* nodes,const V* corrections,V* accumulated,int count,double dt)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=count) return;
    nodes[i].p=nodes[i].p+corrections[i];nodes[i].v=nodes[i].v+corrections[i]*(1./dt);accumulated[i]=accumulated[i]+corrections[i];
}
extern "C" __global__ void embed(const Node* nodes,const Tet* tets,int tetCount,const V* points,Map* maps,
    int count,int begin,int end,int body,int* failure)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=count) return;
    V p=points[i];double best=3.402823466e38;Map result={{-1,-1,-1,-1},{0,0,0,0},{0,0,0},body};
    const int faces[4][3]={{1,2,3},{0,3,2},{0,1,3},{0,2,1}};
    for(int t=0;t<tetCount;++t)
    {
        Tet tet=tets[t];bool eligible=true;for(int j=0;j<4;++j) if(tet.n[j]<begin||tet.n[j]>=end) eligible=false;
        if(!eligible) continue;
        V x[4];for(int j=0;j<4;++j) x[j]=nodes[tet.n[j]].p;
        M rest={{{tet.rest[0],tet.rest[1],tet.rest[2]},{tet.rest[3],tet.rest[4],tet.rest[5]},{tet.rest[6],tet.rest[7],tet.rest[8]}}};
        V b=mul(rest,p-x[0]);double w[4]={1-b.x-b.y-b.z,b.x,b.y,b.z};
        bool inside=true;for(int j=0;j<4;++j) if(w[j]<0) inside=false;
        double distance=0;
        if(!inside)
        {
            distance=3.402823466e38;
            for(int f=0;f<4;++f)
            {
                V a=x[faces[f][0]],b=x[faces[f][1]],c=x[faces[f][2]],bary=closestWeights(p,a,b,c);
                V d=p-(a*bary.x+b*bary.y+c*bary.z);distance=fmin(distance,dot(d,d));
            }
        }
        if(distance<best)
        {
            best=distance;for(int j=0;j<4;++j) {result.nodes[j]=tet.n[j];result.weights[j]=w[j];}
            if(inside) break;
        }
    }
    if(result.nodes[0]<0) atomicExch(failure,3);
    maps[i]=result;
}
