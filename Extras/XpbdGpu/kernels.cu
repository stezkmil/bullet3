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
// Euler conventions from Bullet btGeneric6DofSpring2Constraint.cpp (zlib license).
__device__ double btGetMatrixElem(const M& m,int i) { V c=m.c[i/3];return i%3==0?c.x:(i%3==1?c.y:c.z); }
__device__ bool matrixToEulerXYZ(const M& mat, double* xyz)
{
	// rot =  cy*cz          -cy*sz           sy
	//        cz*sx*sy+cx*sz  cx*cz-sx*sy*sz -cy*sx
	//       -cx*cz*sy+sx*sz  cz*sx+cx*sy*sz  cx*cy

	double fi = btGetMatrixElem(mat, 2);
	if (fi < double(1.0f))
	{
		if (fi > double(-1.0f))
		{
			xyz[0] = atan2(-btGetMatrixElem(mat, 5), btGetMatrixElem(mat, 8));
			xyz[1] = asin(btGetMatrixElem(mat, 2));
			xyz[2] = atan2(-btGetMatrixElem(mat, 1), btGetMatrixElem(mat, 0));
			return true;
		}
		else
		{
			// WARNING.  Not unique.  XA - ZA = -atan2(r10,r11)
			xyz[0] = -atan2(btGetMatrixElem(mat, 3), btGetMatrixElem(mat, 4));
			xyz[1] = -1.57079632679489661923;
			xyz[2] = double(0.0);
			return false;
		}
	}
	else
	{
		// WARNING.  Not unique.  XAngle + ZAngle = atan2(r10,r11)
		xyz[0] = atan2(btGetMatrixElem(mat, 3), btGetMatrixElem(mat, 4));
		xyz[1] = 1.57079632679489661923;
		xyz[2] = 0.0;
	}
	return false;
}

__device__ bool matrixToEulerXZY(const M& mat, double* xyz)
{
	// rot =  cy*cz          -sz           sy*cz
	//        cy*cx*sz+sx*sy  cx*cz        sy*cx*sz-cy*sx
	//        cy*sx*sz-cx*sy  sx*cz        sy*sx*sz+cx*cy

	double fi = btGetMatrixElem(mat, 1);
	if (fi < double(1.0f))
	{
		if (fi > double(-1.0f))
		{
			xyz[0] = atan2(btGetMatrixElem(mat, 7), btGetMatrixElem(mat, 4));
			xyz[1] = atan2(btGetMatrixElem(mat, 2), btGetMatrixElem(mat, 0));
			xyz[2] = asin(-btGetMatrixElem(mat, 1));
			return true;
		}
		else
		{
			xyz[0] = -atan2(-btGetMatrixElem(mat, 6), btGetMatrixElem(mat, 8));
			xyz[1] = double(0.0);
			xyz[2] = 1.57079632679489661923;
			return false;
		}
	}
	else
	{
		xyz[0] = atan2(-btGetMatrixElem(mat, 6), btGetMatrixElem(mat, 8));
		xyz[1] = 0.0;
		xyz[2] = -1.57079632679489661923;
	}
	return false;
}

__device__ bool matrixToEulerYXZ(const M& mat, double* xyz)
{
	// rot =  cy*cz+sy*sx*sz  cz*sy*sx-cy*sz  cx*sy
	//        cx*sz           cx*cz           -sx
	//        cy*sx*sz-cz*sy  sy*sz+cy*cz*sx  cy*cx

	double fi = btGetMatrixElem(mat, 5);
	if (fi < double(1.0f))
	{
		if (fi > double(-1.0f))
		{
			xyz[0] = asin(-btGetMatrixElem(mat, 5));
			xyz[1] = atan2(btGetMatrixElem(mat, 2), btGetMatrixElem(mat, 8));
			xyz[2] = atan2(btGetMatrixElem(mat, 3), btGetMatrixElem(mat, 4));
			return true;
		}
		else
		{
			xyz[0] = 1.57079632679489661923;
			xyz[1] = -atan2(-btGetMatrixElem(mat, 1), btGetMatrixElem(mat, 0));
			xyz[2] = double(0.0);
			return false;
		}
	}
	else
	{
		xyz[0] = -1.57079632679489661923;
		xyz[1] = atan2(-btGetMatrixElem(mat, 1), btGetMatrixElem(mat, 0));
		xyz[2] = 0.0;
	}
	return false;
}

__device__ bool matrixToEulerYZX(const M& mat, double* xyz)
{
	// rot =  cy*cz   sy*sx-cy*cx*sz   cx*sy+cy*sz*sx
	//        sz           cz*cx           -cz*sx
	//        -cz*sy  cy*sx+cx*sy*sz   cy*cx-sy*sz*sx

	double fi = btGetMatrixElem(mat, 3);
	if (fi < double(1.0f))
	{
		if (fi > double(-1.0f))
		{
			xyz[0] = atan2(-btGetMatrixElem(mat, 5), btGetMatrixElem(mat, 4));
			xyz[1] = atan2(-btGetMatrixElem(mat, 6), btGetMatrixElem(mat, 0));
			xyz[2] = asin(btGetMatrixElem(mat, 3));
			return true;
		}
		else
		{
			xyz[0] = double(0.0);
			xyz[1] = -atan2(btGetMatrixElem(mat, 7), btGetMatrixElem(mat, 8));
			xyz[2] = -1.57079632679489661923;
			return false;
		}
	}
	else
	{
		xyz[0] = double(0.0);
		xyz[1] = atan2(btGetMatrixElem(mat, 7), btGetMatrixElem(mat, 8));
		xyz[2] = 1.57079632679489661923;
	}
	return false;
}

__device__ bool matrixToEulerZXY(const M& mat, double* xyz)
{
	// rot =  cz*cy-sz*sx*sy    -cx*sz   cz*sy+cy*sz*sx
	//        cy*sz+cz*sx*sy     cz*cx   sz*sy-cz*xy*sx
	//        -cx*sy              sx     cx*cy

	double fi = btGetMatrixElem(mat, 7);
	if (fi < double(1.0f))
	{
		if (fi > double(-1.0f))
		{
			xyz[0] = asin(btGetMatrixElem(mat, 7));
			xyz[1] = atan2(-btGetMatrixElem(mat, 6), btGetMatrixElem(mat, 8));
			xyz[2] = atan2(-btGetMatrixElem(mat, 1), btGetMatrixElem(mat, 4));
			return true;
		}
		else
		{
			xyz[0] = -1.57079632679489661923;
			xyz[1] = double(0.0);
			xyz[2] = -atan2(btGetMatrixElem(mat, 2), btGetMatrixElem(mat, 0));
			return false;
		}
	}
	else
	{
		xyz[0] = 1.57079632679489661923;
		xyz[1] = double(0.0);
		xyz[2] = atan2(btGetMatrixElem(mat, 2), btGetMatrixElem(mat, 0));
	}
	return false;
}

__device__ bool matrixToEulerZYX(const M& mat, double* xyz)
{
	// rot =  cz*cy   cz*sy*sx-cx*sz   sz*sx+cz*cx*sy
	//        cy*sz   cz*cx+sz*sy*sx   cx*sz*sy-cz*sx
	//        -sy          cy*sx         cy*cx

	double fi = btGetMatrixElem(mat, 6);
	if (fi < double(1.0f))
	{
		if (fi > double(-1.0f))
		{
			xyz[0] = atan2(btGetMatrixElem(mat, 7), btGetMatrixElem(mat, 8));
			xyz[1] = asin(-btGetMatrixElem(mat, 6));
			xyz[2] = atan2(btGetMatrixElem(mat, 3), btGetMatrixElem(mat, 0));
			return true;
		}
		else
		{
			xyz[0] = double(0.0);
			xyz[1] = 1.57079632679489661923;
			xyz[2] = -atan2(btGetMatrixElem(mat, 1), btGetMatrixElem(mat, 2));
			return false;
		}
	}
	else
	{
		xyz[0] = double(0.0);
		xyz[1] = -1.57079632679489661923;
		xyz[2] = atan2(-btGetMatrixElem(mat, 1), -btGetMatrixElem(mat, 2));
	}
	return false;
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

// Generalized rigid mass and point Jacobians follow PhysX's rigid/deformable
// attachment formulation. State and arithmetic use double precision here.
struct Q { double x,y,z,w; };
struct Rigid { V p;Q q;V v,omega;double w;V invI,force,torque;double damping,angularDamping; };
struct Endpoint { int rigid;int nodes[4];double weights[4];V point; };
struct Attachment { Endpoint a,b;double compliance; };
struct JointAxis { double lower,upper,stiffness,damping,equilibrium,velocity,maxForce,servoTarget,bounce,stopERP,stopCFM,motorERP,motorCFM;int flags; };
struct Joint { int bodies[2];V points[2];Q orientations[2];JointAxis axes[6];int rotationOrder; };
__device__ Q qmul(Q a,Q b)
{
    V u={a.x,a.y,a.z},v={b.x,b.y,b.z};V c=u*b.w+v*a.w+cross(u,v);
    return {c.x,c.y,c.z,a.w*b.w-dot(u,v)};
}
__device__ Q qconj(Q q) { return {-q.x,-q.y,-q.z,q.w}; }
__device__ V rotate(Q q,V v) { V u={q.x,q.y,q.z};return v+cross(u,cross(u,v)+v*q.w)*2.; }
__device__ Q advance(Q q,V delta)
{
    double angle=sqrt(dot(delta,delta));double scale=angle>1.e-12?sin(angle*.5)/angle:.5;
    Q r=qmul(Q{delta.x*scale,delta.y*scale,delta.z*scale,cos(angle*.5)},q);
    double norm=sqrt(r.x*r.x+r.y*r.y+r.z*r.z+r.w*r.w);
    return {r.x/norm,r.y/norm,r.z/norm,r.w/norm};
}
__device__ V inverseInertia(Rigid r,V v)
{
    V local=rotate(qconj(r.q),v);return rotate(r.q,V{local.x*r.invI.x,local.y*r.invI.y,local.z*r.invI.z});
}
__device__ V angularMomentum(Rigid r)
{
    V v=rotate(qconj(r.q),r.omega);
    return rotate(r.q,V{r.invI.x>0?v.x/r.invI.x:0,r.invI.y>0?v.y/r.invI.y:0,r.invI.z>0?v.z/r.invI.z:0});
}
extern "C" __global__ void predictRigids(Rigid* bodies,int count,double dt,V gravity)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=count) return;
    Rigid r=bodies[i];if(r.w==0) {r.p=r.p+r.v*dt;r.q=advance(r.q,r.omega*dt);bodies[i]=r;return;}
    r.v=r.v*exp(-r.damping*dt)+(r.force*r.w+gravity)*dt;
    V L=angularMomentum(r)*exp(-r.angularDamping*dt)+r.torque*dt;
    // Midpoint orientation keeps torque-free motion from explicitly amplifying spin.
    Rigid mid=r;V omega=inverseInertia(r,L);
    for(int j=0;j<3;++j) {mid.q=advance(r.q,omega*(dt*.5));omega=inverseInertia(mid,L);}
    r.p=r.p+r.v*dt;r.q=advance(r.q,omega*dt);r.omega=inverseInertia(r,L);bodies[i]=r;
}

__device__ Rigid* activeRigids;
__device__ Rigid* startRigids;
__device__ V* rigidImpulse;
__device__ V* rigidTorqueImpulse;
__device__ int* rigidReferences;
extern "C" __global__ void bindRigids(Rigid* bodies,Rigid* start,V* impulses,V* torques,int* references)
{
    if(blockIdx.x||threadIdx.x) return;
    activeRigids=bodies;startRigids=start;rigidImpulse=impulses;rigidTorqueImpulse=torques;rigidReferences=references;
}
struct Map { int nodes[4];double weights[4];V offset;int body; };
__device__ int rigidIndex(Map m) { return m.nodes[0]<-1?-2-m.nodes[0]:-1; }
__device__ V rigidPosition(Map m) { Rigid r=activeRigids[rigidIndex(m)];return r.p+rotate(r.q,m.offset); }
__device__ V rigidDelta(Map m)
{
    Rigid a=activeRigids[rigidIndex(m)],b=startRigids[rigidIndex(m)];
    return a.p-b.p+rotate(a.q,m.offset)-rotate(b.q,m.offset);
}
__device__ V rigidVelocity(Map m)
{
    Rigid r=activeRigids[rigidIndex(m)];return r.v+cross(r.omega,rotate(r.q,m.offset));
}
__device__ double rigidWeight(Map m,V axis)
{
    if(rigidIndex(m)<0) return 0;Rigid r=activeRigids[rigidIndex(m)];V j=cross(rotate(r.q,m.offset),axis);
    return r.w+dot(j,inverseInertia(r,j));
}
__device__ V rigidResponse(Map m,V impulse)
{
    if(rigidIndex(m)<0) return {0,0,0};
    Rigid r=activeRigids[rigidIndex(m)];V arm=rotate(r.q,m.offset);
    return (impulse*r.w+cross(inverseInertia(r,cross(arm,impulse)),arm))*rigidReferences[rigidIndex(m)];
}
struct Tri { int v[3];int part,index; };
__device__ bool sharedSupport(Map a,Map b)
{
    for(int i=0;i<4;++i) if(a.nodes[i]>=0) for(int j=0;j<4;++j) if(a.nodes[i]==b.nodes[j]) return true;
    return false;
}
__device__ bool adjacentPoint(Map a,Tri b,const Map* maps)
{
    for(int i=0;i<3;++i) if(sharedSupport(a,maps[b.v[i]])) return true;return false;
}
struct Bvh { int left,right,parent,triangle,body; };
struct Box { V low,high; };
__device__ V minimum(V a,V b) { return {fmin(a.x,b.x),fmin(a.y,b.y),fmin(a.z,b.z)}; }
__device__ V maximum(V a,V b) { return {fmax(a.x,b.x),fmax(a.y,b.y),fmax(a.z,b.z)}; }
extern "C" __global__ void skin(const Node* nodes,const Map* maps,V* positions,int count)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=count) return;
    Map m=maps[i];V p=m.offset;
    for(int j=0;j<4;++j) if(m.nodes[j]>=0) p=p+nodes[m.nodes[j]].p*m.weights[j];
    positions[i]=rigidIndex(m)>=0?rigidPosition(m):p;
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
    if(rigidIndex(m)>=0) return rigidPosition(m);
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
    if(rigidIndex(m)>=0) return rigidDelta(m);
    V delta={0,0,0};for(int j=0;j<4;++j) if(m.nodes[j]>=0) delta=delta+accumulated[m.nodes[j]]*m.weights[j];return delta;
}
__device__ V sourceTriangleDelta(Contact contact,const V* accumulated,const Map* maps,const Tri* tris)
{
    if(contact.triangleA<0) return mappedDelta(accumulated,maps[contact.vertex]);
    Tri t=tris[contact.triangleA];
    return mappedDelta(accumulated,maps[t.v[0]])*contact.baryA.x+
        mappedDelta(accumulated,maps[t.v[1]])*contact.baryA.y+mappedDelta(accumulated,maps[t.v[2]])*contact.baryA.z;
}
__device__ V targetTriangleDelta(Contact contact,const V* accumulated,const Map* maps,const Tri* tris)
{
    Tri t=tris[contact.triangle];
    return mappedDelta(accumulated,maps[t.v[0]])*contact.bary.x+
        mappedDelta(accumulated,maps[t.v[1]])*contact.bary.y+mappedDelta(accumulated,maps[t.v[2]])*contact.bary.z;
}
__device__ double separation(Contact contact,const Node* nodes,const V* accumulated,const Map* maps,const Tri* tris)
{
    // PhysX tracks contact motion using accumulated displacement, avoiding the
    // loss of small corrections when subtracting rounded world positions.
    if(contact.prepared) return dot((contact.prepared>=2?sourceTriangleDelta(contact,accumulated,maps,tris):mappedDelta(accumulated,contact.simA))-
        (contact.prepared==3?targetTriangleDelta(contact,accumulated,maps,tris):mappedDelta(accumulated,contact.simB)),contact.normal)+contact.bias;
    Tri t=tris[contact.triangle];
    V point;
    if(contact.triangleA<0) point=mappedPosition(nodes,maps[contact.vertex]);
    else { Tri a=tris[contact.triangleA];point=mappedPosition(nodes,maps[a.v[0]])*contact.baryA.x+mappedPosition(nodes,maps[a.v[1]])*contact.baryA.y+mappedPosition(nodes,maps[a.v[2]])*contact.baryA.z; }
    V delta=point-mappedPosition(nodes,maps[t.v[0]])*contact.bary.x-
        mappedPosition(nodes,maps[t.v[1]])*contact.bary.y-mappedPosition(nodes,maps[t.v[2]])*contact.bary.z;
    return dot(delta,contact.normal);
}
__device__ bool pairDisabled(int a,int b,const unsigned long long* pairs,int count,const int* selfBodies,int selfBodyCount)
{
    if(a<0||b<0) return false; // A mixed BVH branch must still be traversed.
    if(a==b) {int lo=0,hi=selfBodyCount;while(lo<hi){int m=lo+(hi-lo)/2;if(selfBodies[m]<a)lo=m+1;else hi=m;}return lo==selfBodyCount||selfBodies[lo]!=a;}
    if(a>b) { int t=a;a=b;b=t; }
    unsigned long long key=(static_cast<unsigned long long>(a)<<32)|static_cast<unsigned int>(b);
    int lo=0,hi=count;
    while(lo<hi) { int mid=lo+(hi-lo)/2;if(pairs[mid]<key) lo=mid+1;else hi=mid; }
    return lo<count&&pairs[lo]==key;
}
extern "C" __global__ void discoverContacts(const Node* nodes,const Map* maps,const V* positions,const Tri* tris,
    const Bvh* tree,const Box* boxes,Contact* cache,int* contactCount,int capacity,int count,double margin,int* failure,const int* incident,
    const unsigned long long* disabledPairs,int disabledPairCount,const int* selfBodies,int selfBodyCount)
{
    int vertex=blockIdx.x*blockDim.x+threadIdx.x;if(vertex>=count) return;
    Map point=maps[vertex];V p=positions[vertex];
    bool movingPoint=rigidIndex(point)>=0;for(int j=0;j<4;++j) movingPoint|=point.nodes[j]>=0;
    if(incident[vertex]&&!movingPoint) return;
    // Static vertices are queried too: narrow pegs must also collide against dynamic faces.
    int stack[64],top=0;stack[top++]=0;
    while(top)
    {
        int index=stack[--top];Bvh node=tree[index];Box box=boxes[index];
        if(pairDisabled(node.body,point.body,disabledPairs,disabledPairCount,selfBodies,selfBodyCount) || p.x<box.low.x-margin||p.y<box.low.y-margin||p.z<box.low.z-margin||
           p.x>box.high.x+margin||p.y>box.high.y+margin||p.z>box.high.z+margin) continue;
        if(node.triangle<0)
        {
            if(top>61) { atomicExch(failure,2);return; }
            stack[top++]=node.left;stack[top++]=node.right;continue;
        }
        Tri t=tris[node.triangle];
        if(point.body==node.body&&adjacentPoint(point,t,maps)) continue;
        // Triangle witnesses alone do not constrain every corner of a broad contact patch.
        if(incident[vertex])
        {
            bool moving=rigidIndex(maps[t.v[0]])>=0;for(int k=0;k<3;++k) for(int j=0;j<4;++j) moving|=maps[t.v[k]].nodes[j]>=0;
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
        bool dynamic=rigidIndex(point)>=0||rigidIndex(maps[t.v[0]])>=0;
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
__device__ double closestTriangleWeights(const V* a,const V* b,V& wa,V& wb)
{
    double best=1.e100;wa=wb=corner(0);
    for(int k=0;k<3;++k)
    {
        closestCandidate(a,b,corner(k),closestWeights(a[k],b[0],b[1],b[2]),best,wa,wb);
        closestCandidate(a,b,closestWeights(b[k],a[0],a[1],a[2]),corner(k),best,wa,wb);
        edgeFaceCandidate(a,b,k,false,best,wa,wb);edgeFaceCandidate(b,a,k,true,best,wa,wb);
        for(int j=0;j<3;++j)
        {
            double u,v;edgeParameters(a[k],a[(k+1)%3]-a[k],b[j],b[(j+1)%3]-b[j],u,v);
            closestCandidate(a,b,corner(k)*(1-u)+corner((k+1)%3)*u,corner(j)*(1-v)+corner((j+1)%3)*v,best,wa,wb);
        }
    }
    return best;
}
__device__ bool movingTriangle(Tri t,const Map* maps,const Node* nodes)
{
    if(rigidIndex(maps[t.v[0]])>=0) return activeRigids[rigidIndex(maps[t.v[0]])].w>0;
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
    const unsigned long long* disabledPairs,int disabledPairCount,const int* selfBodies,int selfBodyCount)
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
        if(pairDisabled(node.body,body,disabledPairs,disabledPairCount,selfBodies,selfBodyCount)||bound.high.x<box.low.x-margin||bound.low.x>box.high.x+margin||
            bound.high.y<box.low.y-margin||bound.low.y>box.high.y+margin||bound.high.z<box.low.z-margin||bound.low.z>box.high.z+margin) continue;
        if(node.triangle<0)
        {
            if(top>61) { atomicExch(failure,2);return; }
            stack[top++]=node.left;stack[top++]=node.right;continue;
        }
        Tri tb=tris[node.triangle];bool dynamicB=movingTriangle(tb,maps,nodes);
        if(body==node.body)
        {
            bool adjacent=false;for(int k=0;k<3;++k) adjacent|=adjacentPoint(maps[ta.v[k]],tb,maps);
            if(adjacent) continue;
        }
        // Volume vertex contacts cannot cover crossings between triangle edges.
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
    Contact* cache,int* contactCount,int capacity,double margin,int* failure,const unsigned long long* disabledPairs,int disabledPairCount,const int* selfBodies,int selfBodyCount)
{
    int job=blockIdx.x*blockDim.x+threadIdx.x;
    if(job>=vertexCount*volumeCount) return;
    int vertex=job%vertexCount,target=volumeBodies[job/vertexCount],source=maps[vertex].body;
    if(!incident[vertex]||pairDisabled(source,target,disabledPairs,disabledPairCount,selfBodies,selfBodyCount)) return;
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
        Tri t=tris[node.triangle];
        if(source==target&&adjacentPoint(maps[vertex],t,maps)) continue;
        V a=positions[t.v[0]],b=positions[t.v[1]],c=positions[t.v[2]];
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
                    !isfinite(pb.x)||!isfinite(pb.y)||!isfinite(pb.z))
                {
                    c.state|=32; // Complete singular triangle witnesses in a separate kernel.
                }
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
extern "C" __global__ void completeTriangleContacts(const V* positions,const Tri* tris,Contact* cache,const int* contactCount,int capacity,double margin,int* failure)
{
    int count=min(*contactCount,capacity);
    for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<count;i+=blockDim.x*gridDim.x)
    {
        if(!(cache[i].state&32)) continue;
        Contact c=cache[i];c.state&=~32;Tri ta=tris[c.triangleA],tb=tris[c.triangle];
        // EPA can have a singular witness for intersecting zero-thickness triangles.
        // Keeping this rare geometric fallback separate preserves the GJK kernel's register budget.
        V a[3]={positions[ta.v[0]],positions[ta.v[1]],positions[ta.v[2]]};
        V b[3]={positions[tb.v[0]],positions[tb.v[1]],positions[tb.v[2]]};
        double d2=closestTriangleWeights(a,b,c.baryA,c.bary);
        if(!isfinite(d2)||!isfinite(dot(c.normal,c.normal))) {atomicExch(failure,1);c.state=4;}
        else
        {
            c.simA.offset=trianglePoint(a,c.baryA);c.simB.offset=trianglePoint(b,c.bary);c.state|=8;
            if(d2>margin*margin) c.state|=4;
        }
        cache[i]=c;
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
        // Contact gradients must follow the embeddings that move both soft surfaces.
        bool surfaceSupport=dynamicA&&dynamicB;
        if((c.state&8)&&!surfaceSupport) { pa=c.simA.offset;pb=c.simB.offset; }
        c.simA={{-1,-1,-1,-1},{0,0,0,0},pa,bodyA};c.simB={{-1,-1,-1,-1},{0,0,0,0},pb,bodyB};
        // PhysX's rigid-volume writer remaps the rigid witness point onto the volume.
        V remapA=dynamicA&&!dynamicB?pb:pa;
        bool materialSupport=dynamicA&&!dynamicB&&c.triangleA>=0;
        if(materialSupport)
        {
            Tri ta=tris[c.triangleA];
            for(int j=0;j<3;++j) for(int k=0;k<4;++k) if(maps[ta.v[j]].weights[k]<0) materialSupport=false;
        }
        bool a=!dynamicA||materialSupport||surfaceSupport;
        if(c.state&16) { c.simA=maps[c.vertex];a=true; }
        else if(dynamicA&&!materialSupport&&!surfaceSupport) a=remapContactPoint(remapA,bodyA,nodes,tets,tetBodies,tetCount,c.triangleA,offsets,candidates,c.simA);
        bool b=!dynamicB||surfaceSupport||remapContactPoint(pb,bodyB,nodes,tets,tetBodies,tetCount,c.triangle,offsets,candidates,c.simB);
        // Keep the force on the collision triangle's material when nearby elements fold over it.
        c.prepared=surfaceSupport?3:(a&&b?(materialSupport?2:1):0);
        int rigidA=rigidIndex(maps[c.vertex]),rigidB=rigidIndex(maps[tb.v[0]]);
        if(rigidA>=0) {c.simA.nodes[0]=-2-rigidA;c.simA.offset=rotate(qconj(activeRigids[rigidA].q),pa-activeRigids[rigidA].p);}
        if(rigidB>=0) {c.simB.nodes[0]=-2-rigidB;c.simB.offset=rotate(qconj(activeRigids[rigidB].q),pb-activeRigids[rigidB].p);}
        c.impulse=0;c.state=(!(dynamicA||(rigidA>=0&&activeRigids[rigidA].w>0))||!(dynamicB||(rigidB>=0&&activeRigids[rigidB].w>0)))?1:0;
        c.tangentImpulse={0,0,0};
        c.friction=bodyA>=0&&bodyA<frictionCount&&bodyB>=0&&bodyB<frictionCount?
            fmin(10.,friction[bodyA]*friction[bodyB]):0.;
        c.bias=dot(pa-pb,c.normal)-dot(((materialSupport||surfaceSupport)?sourceTriangleDelta(c,accumulated,maps,tris):mappedDelta(accumulated,c.simA))-
            (surfaceSupport?targetTriangleDelta(c,accumulated,maps,tris):mappedDelta(accumulated,c.simB)),c.normal);
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
            if(rigidIndex(c.simA)>=0) velocity=velocity+rigidVelocity(c.simA);
            if(rigidIndex(c.simB)>=0) velocity=velocity-rigidVelocity(c.simB);
            C=fmax(0.,C)/dt+dot(velocity,c.normal);
        }
        if(C>=0&&!(c.state&2)&&c.impulse<=0) continue;
        if(c.state&1) cache[i].state|=2;
        if(rigidIndex(c.simA)>=0) atomicAdd(rigidReferences+rigidIndex(c.simA),1);
        if(rigidIndex(c.simB)>=0) atomicAdd(rigidReferences+rigidIndex(c.simB),1);
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
            if(rigidIndex(c.simA)>=0) velocity=velocity+rigidVelocity(c.simA);
            if(rigidIndex(c.simB)>=0) velocity=velocity-rigidVelocity(c.simB);
            C=fmax(0.,C)/dt+dot(velocity,c.normal);
        }
        if(C>=0&&!(c.state&2)&&c.impulse<=0) continue;
        for(int j=0;j<used;++j) denominator+=weights[j]*weights[j]*nodes[ids[j]].w*references[ids[j]];
        double softWeight=denominator;
        if(rigidIndex(c.simA)>=0) denominator+=rigidWeight(c.simA,c.normal)*rigidReferences[rigidIndex(c.simA)];
        if(rigidIndex(c.simB)>=0) denominator+=rigidWeight(c.simB,c.normal)*rigidReferences[rigidIndex(c.simB)];
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
            if(rigidIndex(c.simA)>=0) relative=relative+(velocityOnly?rigidVelocity(c.simA):rigidDelta(c.simA)*(1./dt));
            if(rigidIndex(c.simB)>=0) relative=relative-(velocityOnly?rigidVelocity(c.simB):rigidDelta(c.simB)*(1./dt));
            V tangent=relative-c.normal*dot(relative,c.normal);
            V target=c.tangentImpulse;double tangentLength=sqrt(dot(tangent,tangent));
            if(tangentLength>0)
            {
                V u=tangent*(1./tangentLength),v=cross(c.normal,u);
                V responseU=u*softWeight+rigidResponse(c.simA,u)+rigidResponse(c.simB,u);
                V responseV=v*softWeight+rigidResponse(c.simA,v)+rigidResponse(c.simB,v);
                double a=dot(u,responseU),b=dot(v,responseU),d=dot(v,responseV),det=a*d-b*b;
                if(det>0) target=target-u*(tangentLength*d/det)+v*(tangentLength*b/det);
            }
            double limit=c.friction*cache[i].impulse,length=sqrt(dot(target,target));
            if(length>limit) target=target*(limit/length);
            deltaImpulse=deltaImpulse+target-c.tangentImpulse;
            cache[i].tangentImpulse=target;
        }
        Map ends[2]={c.simA,c.simB};
        for(int side=0;side<2;++side) if(rigidIndex(ends[side])>=0)
        {
            int id=rigidIndex(ends[side]);V impulse=deltaImpulse*(side?-1.:1.);
            add(rigidImpulse+id,impulse);add(rigidTorqueImpulse+id,cross(rotate(activeRigids[id].q,ends[side].offset),impulse));
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

__device__ V endpointPosition(Endpoint e,const Node* nodes,const Rigid* bodies)
{
    if(e.rigid>=0) {Rigid r=bodies[e.rigid];return r.p+rotate(r.q,e.point);}
    V p=e.point;for(int j=0;j<4;++j) if(e.nodes[j]>=0) p=p+nodes[e.nodes[j]].p*e.weights[j];return p;
}
__device__ V endpointResponse(Endpoint e,const Node* nodes,const Rigid* bodies,V impulse)
{
    if(e.rigid>=0) {Rigid r=bodies[e.rigid];V arm=rotate(r.q,e.point);return impulse*r.w+cross(inverseInertia(r,cross(arm,impulse)),arm);}
    double w=0;for(int j=0;j<4;++j) if(e.nodes[j]>=0) w+=nodes[e.nodes[j]].w*e.weights[j]*e.weights[j];return impulse*w;
}
__device__ void correctEndpoint(Endpoint e,Node* nodes,Rigid* bodies,V* accumulated,V impulse,double dt)
{
    if(e.rigid>=0)
    {
        Rigid r=bodies[e.rigid];if(r.w==0) return;
        V L=angularMomentum(r)+cross(rotate(r.q,e.point),impulse);
        V dx=impulse*(r.w*dt),angle=inverseInertia(r,cross(rotate(r.q,e.point),impulse))*dt;
        r.v=r.v+impulse*r.w;r.p=r.p+dx;r.q=advance(r.q,angle);r.omega=inverseInertia(r,L);bodies[e.rigid]=r;
    }
    else for(int j=0;j<4;++j) if(e.nodes[j]>=0)
    {
        int n=e.nodes[j];V dv=impulse*(nodes[n].w*e.weights[j]);nodes[n].v=nodes[n].v+dv;
        V dx=dv*dt;nodes[n].p=nodes[n].p+dx;accumulated[n]=accumulated[n]+dx;
    }
}
extern "C" __global__ void solveAttachments(Node* nodes,Rigid* bodies,V* accumulated,const Attachment* attachments,int begin,int count,double dt)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=count) return;Attachment a=attachments[begin+i];
    V error=endpointPosition(a.a,nodes,bodies)-endpointPosition(a.b,nodes,bodies);
    V axis[3]={{1,0,0},{0,1,0},{0,0,1}},columns[3];double alpha=a.compliance/(dt*dt);
    for(int j=0;j<3;++j) columns[j]=endpointResponse(a.a,nodes,bodies,axis[j])+endpointResponse(a.b,nodes,bodies,axis[j])+axis[j]*alpha;
    double det=dot(columns[0],cross(columns[1],columns[2]));if(det<=0) return;
    V impulse={-dot(error,cross(columns[1],columns[2]))/(det*dt),-dot(columns[0],cross(error,columns[2]))/(det*dt),-dot(columns[0],cross(columns[1],error))/(det*dt)};
    correctEndpoint(a.a,nodes,bodies,accumulated,impulse,dt);correctEndpoint(a.b,nodes,bodies,accumulated,impulse*-1,dt);
}
__device__ Rigid jointBody(const Rigid* bodies,int id)
{
    if(id>=0) return bodies[id];Rigid r={};r.q.w=1;return r;
}
__device__ void jointImpulse(Rigid& r,V linear,V angular,double impulse,double dt)
{
    if(r.w==0) return;
    V L=angularMomentum(r)+angular*impulse;
    V dv=linear*(r.w*impulse);r.v=r.v+dv;r.p=r.p+dv*dt;
    r.q=advance(r.q,inverseInertia(r,angular)*(impulse*dt));r.omega=inverseInertia(r,L);
}
extern "C" __global__ void solveJoints(Rigid* bodies,const Joint* joints,int begin,int count,double dt)
{
    int index=blockIdx.x*blockDim.x+threadIdx.x;if(index>=count) return;
    Joint j=joints[begin+index];Rigid a=jointBody(bodies,j.bodies[0]),b=jointBody(bodies,j.bodies[1]);
    const V basis[3]={{1,0,0},{0,1,0},{0,0,1}};
    const int orders[6][3]={{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    for(int row=0;row<6;++row)
    {
        JointAxis c=j.axes[row];if(c.lower>c.upper&&!(c.flags&3)) continue;
        Q qa=qmul(a.q,j.orientations[0]),qb=qmul(b.q,j.orientations[1]);
        V pa=a.p+rotate(a.q,j.points[0]),pb=b.p+rotate(b.q,j.points[1]);
        V la={0,0,0},lb={0,0,0},aa,ab;double position;
        if(row<3)
        {
            V axis=rotate(qa,basis[row]);position=dot(pb-pa,axis);la=axis*-1;lb=axis;
            // The axis rotates with frame A: both reaction forces act at B's pivot.
            aa=cross(pb-a.p,la);ab=cross(pb-b.p,lb);
        }
        else
        {
            Q relative=qmul(qconj(qa),qb);M m={{rotate(relative,basis[0]),rotate(relative,basis[1]),rotate(relative,basis[2])}};
            double angles[3];
            switch(j.rotationOrder) {
                case 0:matrixToEulerXYZ(m,angles);break;case 1:matrixToEulerXZY(m,angles);break;
                case 2:matrixToEulerYXZ(m,angles);break;case 3:matrixToEulerYZX(m,angles);break;
                case 4:matrixToEulerZXY(m,angles);break;default:matrixToEulerZYX(m,angles);break;
            }
            position=angles[row-3];int first=orders[j.rotationOrder][0],middle=orders[j.rotationOrder][1],last=orders[j.rotationOrder][2];
            V x=rotate(qb,basis[first]),z=rotate(qa,basis[last]);double cosine=dot(x,z),denominator=fmax(1.e-12,1-cosine*cosine);
            if(row-3==first) aa=(x-z*cosine)*(1./denominator);
            else if(row-3==last) aa=(z-x*cosine)*(1./denominator);
            else {double parity=(first==0&&middle==1)||(first==1&&middle==2)||(first==2&&middle==0)?1.:-1.;aa=cross(z,x)*(parity/sqrt(denominator));}
            ab=aa*-1;
        }
        double weight=a.w*dot(la,la)+b.w*dot(lb,lb)+dot(aa,inverseInertia(a,aa))+dot(ab,inverseInertia(b,ab));
        if(weight<=0) continue;
        double rate=dot(la,a.v)+dot(lb,b.v)+dot(aa,a.omega)+dot(ab,b.omega),impulse=0;
        if(c.flags&1)
        {
            // Bullet Spring2 optional sampling caps, evaluated at the TGS microstep.
            if(c.flags&8) c.stiffness=fmin(c.stiffness,1/(16*dt*dt*weight));
            if(c.flags&16) c.damping=fmin(c.damping,1/(dt*weight));
            double spring=-(c.stiffness*dt*(position-c.equilibrium)+c.damping*dt*rate)/(1+(c.stiffness*dt*dt+c.damping*dt)*weight);
            impulse+=spring;rate+=spring*weight;position+=spring*weight*dt;
        }
        if(c.flags&2)
        {
            double desired=c.velocity;
            if(c.flags&4) {double speed=fabs(c.velocity);desired=fmax(-speed,fmin(speed,rate+(c.servoTarget-position)*c.motorERP/dt));}
            double motor=(desired-rate)/(weight+c.motorCFM/dt),limit=c.maxForce*dt;
            motor=fmax(-limit,fmin(limit,motor));impulse+=motor;rate+=motor*weight;position+=motor*weight*dt;
        }
        if(c.lower<=c.upper)
        {
            double error=position-fmax(c.lower,fmin(c.upper,position));
            if(error!=0||c.lower==c.upper)
            {
                double correction=-c.stopERP*error/dt;
                if(c.lower<c.upper&&error*rate>0) correction+=-c.bounce*rate;
                impulse+=correction/(weight+c.stopCFM/dt);
            }
        }
        jointImpulse(a,la,aa,impulse,dt);jointImpulse(b,lb,ab,impulse,dt);
    }
    if(j.bodies[0]>=0) bodies[j.bodies[0]]=a;if(j.bodies[1]>=0) bodies[j.bodies[1]]=b;
}
extern "C" __global__ void checkRigids(const Rigid* bodies,int count,int* failure)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=count) return;Rigid r=bodies[i];
    if(!isfinite(r.p.x)||!isfinite(r.p.y)||!isfinite(r.p.z)||!isfinite(r.v.x)||!isfinite(r.v.y)||!isfinite(r.v.z)||
       !isfinite(r.omega.x)||!isfinite(r.omega.y)||!isfinite(r.omega.z)||!isfinite(r.q.x)||!isfinite(r.q.y)||!isfinite(r.q.z)||!isfinite(r.q.w)) atomicExch(failure,1);
}

extern "C" __global__ void applyRigidContacts(Rigid* bodies,int count,double dt,int velocityOnly)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=count) return;Rigid r=bodies[i];if(r.w==0) return;
    V dv=rigidImpulse[i]*r.w,L=angularMomentum(r)+rigidTorqueImpulse[i];r.v=r.v+dv;
    if(!velocityOnly) {r.p=r.p+dv*dt;r.q=advance(r.q,inverseInertia(r,rigidTorqueImpulse[i])*dt);}
    r.omega=inverseInertia(r,L);bodies[i]=r;
}

struct ContactReport
{
    int bodyA,bodyB,partA,indexA,partB,indexB;
    V pointA,pointB,normal;
    double distance;
};
static_assert(sizeof(Tri)==20 && sizeof(ContactReport)==104,"Contact report GPU ABI mismatch");
__device__ void appendReport(ContactReport report,ContactReport* output,int* count,int capacity,int* failure)
{
    int slot=atomicAdd(count,1);
    if(slot>=capacity) { atomicExch(failure,4);return; }
    output[slot]=report;
}
extern "C" __global__ void reportContacts(const V* positions,const Map* maps,const Tri* tris,const int* vertexTriangles,
    const Contact* contacts,const int* count,ContactReport* output,int* outputCount,int capacity,int* failure)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=min(*count,capacity)) return;
    Contact c=contacts[i];if(c.state&4) return;
    int source=c.triangleA>=0?c.triangleA:vertexTriangles[c.vertex];if(source<0) return;
    Tri a=tris[source],b=tris[c.triangle];
    V pa=c.triangleA>=0?positions[a.v[0]]*c.baryA.x+positions[a.v[1]]*c.baryA.y+positions[a.v[2]]*c.baryA.z:positions[c.vertex];
    V pb=positions[b.v[0]]*c.bary.x+positions[b.v[1]]*c.bary.y+positions[b.v[2]]*c.bary.z;
    double distance=dot(pa-pb,c.normal);
    if(distance>0||!isfinite(distance)) return;
    appendReport({maps[a.v[0]].body,maps[b.v[0]].body,a.part,a.index,b.part,b.index,pa,pb,c.normal,distance},output,outputCount,capacity,failure);
}
// Notification-only pairs share the response BVH, but never enter the constraint list.
// The first triangle pair in mesh order wins, independently of CUDA scheduling.
struct NotificationPair { int x,y,touching; };
extern "C" __global__ void reportNotifications(const Map* maps,const V* positions,const Tri* tris,const Bvh* tree,const Box* boxes,
    const NotificationPair* pairs,int pairCount,int triangleCount,unsigned long long* winners,ContactReport* output,int* outputCount,int capacity,int* failure)
{
    int job=blockIdx.x*blockDim.x+threadIdx.x;if(job>=pairCount*triangleCount) return;
    int pair=job/triangleCount,source=job%triangleCount;NotificationPair bodies=pairs[pair];
    Tri ta=tris[source];if(maps[ta.v[0]].body!=bodies.x||(atomicAdd(winners+pair,0ULL)>>32)<static_cast<unsigned>(source)) return;
    V a[3]={positions[ta.v[0]],positions[ta.v[1]],positions[ta.v[2]]};
    V an=cross(a[1]-a[0],a[2]-a[0]);double a2=dot(an,an);if(a2<1.e-20) return;an=an*rsqrt(a2);
    Box bound={minimum(a[0],minimum(a[1],a[2])),maximum(a[0],maximum(a[1],a[2]))};
    int stack[64],top=0;stack[top++]=0;
    while(top&&(atomicAdd(winners+pair,0ULL)>>32)>=static_cast<unsigned>(source))
    {
        int index=stack[--top];Bvh node=tree[index];Box box=boxes[index];
        if((node.body>=0&&node.body!=bodies.y)||bound.high.x<box.low.x||bound.low.x>box.high.x||
            bound.high.y<box.low.y||bound.low.y>box.high.y||bound.high.z<box.low.z||bound.low.z>box.high.z) continue;
        if(node.triangle<0)
        {
            if(top>61) {atomicExch(failure,2);return;}
            stack[top++]=node.left;stack[top++]=node.right;continue;
        }
        Tri tb=tris[node.triangle];V b[3]={positions[tb.v[0]],positions[tb.v[1]],positions[tb.v[2]]};
        V bn=cross(b[1]-b[0],b[2]-b[0]);double b2=dot(bn,bn);if(b2<1.e-20) continue;bn=bn*rsqrt(b2);
        if(triangleSeparation(a,b,an,bn,0)>0) continue;
        V wa,wb;double best=closestTriangleWeights(a,b,wa,wb);
        if(best>1.e-16) continue;
        int deepest=0;
        for(int k=1;k<3;++k) if(dot(a[k]-b[0],bn)<dot(a[deepest]-b[0],bn)) deepest=k;
        double depth=dot(a[deepest]-b[0],bn);
        if(depth>0||(depth==0&&!bodies.touching)) continue;
        atomicMin(winners+pair,(static_cast<unsigned long long>(source)<<32)|static_cast<unsigned>(node.triangle));
    }
}

extern "C" __global__ void finishNotifications(const V* positions,const Tri* tris,const NotificationPair* pairs,int pairCount,
    const unsigned long long* winners,ContactReport* output,int* outputCount,int capacity,int* failure)
{
    int pair=blockIdx.x*blockDim.x+threadIdx.x;if(pair>=pairCount||winners[pair]==~0ULL) return;
    Tri ta=tris[winners[pair]>>32],tb=tris[static_cast<unsigned>(winners[pair])];
    V a[3]={positions[ta.v[0]],positions[ta.v[1]],positions[ta.v[2]]};
    V b[3]={positions[tb.v[0]],positions[tb.v[1]],positions[tb.v[2]]};
    V bn=cross(b[1]-b[0],b[2]-b[0]);bn=bn*rsqrt(dot(bn,bn));
    int deepest=0;for(int k=1;k<3;++k) if(dot(a[k]-b[0],bn)<dot(a[deepest]-b[0],bn)) deepest=k;
    V pa=a[deepest],pb=trianglePoint(b,closestWeights(pa,b[0],b[1],b[2]));
    appendReport({pairs[pair].x,pairs[pair].y,ta.part,ta.index,tb.part,tb.index,pa,pb,bn,dot(pa-pb,bn)},output,outputCount,capacity,failure);
}
