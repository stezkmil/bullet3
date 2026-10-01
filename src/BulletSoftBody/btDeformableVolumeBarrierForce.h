#ifndef BT_DEFORMABLE_VOLUME_BARRIER_FORCE_H
#define BT_DEFORMABLE_VOLUME_BARRIER_FORCE_H
#include "btDeformableLagrangianForce.h"
#include <cmath>
#include <limits>

// Additional compression resistance for the experimental coupled solver.
// Zero above activation; divergent at the same floor used by step validation.
class btDeformableVolumeBarrierForce : public btDeformableLagrangianForce
{
public:
    struct Material { btSoftBody* body; btScalar bulk; };
    btAlignedObjectArray<Material> materials;
    struct CachedElement
    {
        const btSoftBody::Tetra* tet;
        btMatrix3x3 f;
        btVector3 gradient[4];
        btScalar first, second, weight;
    };
    btAlignedObjectArray<CachedElement> cached;
    bool cacheActive = false;
    void prepareImplicitForceDifferential(btScalar) override
    {
        cached.clear();
        for(int b=0;b<materials.size();++b) if(active(materials[b].body))
            for(int e=0;e<materials[b].body->m_tetras.size();++e)
            {
                const auto& t=materials[b].body->m_tetras[e];const auto f=deformation(t);
                if(f.determinant()>=activation())continue;
                CachedElement entry;entry.tet=&t;entry.f=f;entry.weight=materials[b].bulk*t.m_element_measure;
                btScalar energy;density(f.determinant(),energy,entry.first,entry.second);gradients(t,f,entry.gradient);cached.push_back(entry);
            }
        cacheActive=true;
    }
    void finishImplicitForceDifferential() override { cacheActive=false;cached.clear(); }

    static btScalar floor() { return btScalar(.05); }
    static btScalar activation() { return btScalar(.5); }
    btDeformableLagrangianForceType getForceType() override { return BT_VOLUME_BARRIER_FORCE; }
    static bool movable(const btSoftBody::Node* n) { return n->m_im > 0 && n->m_frozen <= 0; }
    static bool active(const btSoftBody* b) { return b->isActive() && !b->isStaticObject(); }
    static btMatrix3x3 deformation(const btSoftBody::Tetra& t)
    {
        return btMatrix3x3(t.m_n[1]->m_q-t.m_n[0]->m_q,t.m_n[2]->m_q-t.m_n[0]->m_q,t.m_n[3]->m_q-t.m_n[0]->m_q).transpose()*t.m_Dm_inverse;
    }
    static btMatrix3x3 cofactor(const btMatrix3x3& f)
    {
        return btMatrix3x3(f.getColumn(1).cross(f.getColumn(2)),f.getColumn(2).cross(f.getColumn(0)),f.getColumn(0).cross(f.getColumn(1))).transpose();
    }
    static btMatrix3x3 cofactorDifferential(const btMatrix3x3& f,const btMatrix3x3& df)
    {
        return btMatrix3x3(df.getColumn(1).cross(f.getColumn(2))+f.getColumn(1).cross(df.getColumn(2)),
            df.getColumn(2).cross(f.getColumn(0))+f.getColumn(2).cross(df.getColumn(0)),
            df.getColumn(0).cross(f.getColumn(1))+f.getColumn(0).cross(df.getColumn(1))).transpose();
    }
    static void density(btScalar j, btScalar& energy, btScalar& first, btScalar& second)
    {
        energy=first=second=0;
        if (j >= activation()) return;
        if (!(j > floor())) { energy=std::numeric_limits<btScalar>::infinity(); return; }
        const btScalar x=j-floor(), width=activation()-floor(), z=j-activation();
        const btScalar logarithm=btLog(x/width);
        energy=-z*z*logarithm;
        first=-2*z*logarithm-z*z/x;
        second=-2*logarithm-4*z/x+z*z/(x*x);
    }
    static void gradients(const btSoftBody::Tetra& t,const btMatrix3x3& f,btVector3 g[4])
    {
        const btMatrix3x3 a=cofactor(f)*t.m_Dm_inverse.transpose();
        for(int k=1;k<4;++k)g[k]=a.getColumn(k-1);
        g[0]=-g[1]-g[2]-g[3];
    }
    void addScaledForces(btScalar scale,TVStack& force) override
    {
        for(int b=0;b<materials.size();++b) if(active(materials[b].body))
            for(int e=0;e<materials[b].body->m_tetras.size();++e)
            {
                const auto& t=materials[b].body->m_tetras[e];const auto f=deformation(t);
                btScalar energy,first,second;density(f.determinant(),energy,first,second);
                btVector3 g[4];gradients(t,f,g);
                for(int k=0;k<4;++k)force[t.m_n[k]->index]-=scale*materials[b].bulk*t.m_element_measure*first*g[k];
            }
    }
    void addScaledElasticForceDifferential(btScalar scale,const TVStack& dx,TVStack& out) override
    {
        const bool temporary = !cacheActive;
        if(temporary)prepareImplicitForceDifferential(0);
        for(int e=0;e<cached.size();++e)
        {
            const auto& entry=cached[e];const auto& t=*entry.tet;btScalar dj=0;
            for(int k=0;k<4;++k)dj+=entry.gradient[k].dot(dx[t.m_n[k]->index]);
            const auto df=Ds(t.m_n[0]->index,t.m_n[1]->index,t.m_n[2]->index,t.m_n[3]->index,dx)*t.m_Dm_inverse;
            const auto dg123=cofactorDifferential(entry.f,df)*t.m_Dm_inverse.transpose();
            btVector3 dg[4];for(int k=1;k<4;++k)dg[k]=dg123.getColumn(k-1);dg[0]=-dg[1]-dg[2]-dg[3];
            for(int k=0;k<4;++k)out[t.m_n[k]->index]-=scale*entry.weight*(entry.second*dj*entry.gradient[k]+entry.first*dg[k]);
        }
        if(temporary)finishImplicitForceDifferential();
    }

    bool addImplicitForceDifferentialBlocks(btScalar dt,btAlignedObjectArray<btMatrix3x3>& blocks) override
    {
        for(int b=0;b<materials.size();++b) if(active(materials[b].body))
            for(int e=0;e<materials[b].body->m_tetras.size();++e)
            {
                const auto& t=materials[b].body->m_tetras[e];const auto f=deformation(t);
                btScalar energy,first,second;density(f.determinant(),energy,first,second);
                btVector3 g[4];gradients(t,f,g);
                // Determinant is affine in one vertex, so its diagonal Hessian is zero.
                for(int k=0;k<4;++k)if(movable(t.m_n[k]))
                    blocks[t.m_n[k]->index]+=btMatrix3x3(g[k]*g[k].x(),g[k]*g[k].y(),g[k]*g[k].z())*(dt*dt*materials[b].bulk*t.m_element_measure*second);
            }
        return true;
    }
    double totalEnergy(btScalar) override
    {
        double sum=0;
        for(int b=0;b<materials.size();++b)if(active(materials[b].body))
            for(int e=0;e<materials[b].body->m_tetras.size();++e)
            {
                const auto& t=materials[b].body->m_tetras[e];btScalar energy,first,second;
                density(deformation(t).determinant(),energy,first,second);sum+=materials[b].bulk*t.m_element_measure*energy;
            }
        return sum;
    }
    bool admissible() const
    {
        for(int b=0;b<materials.size();++b)if(active(materials[b].body))
            for(int e=0;e<materials[b].body->m_tetras.size();++e)
            { const btScalar j=deformation(materials[b].body->m_tetras[e]).determinant();if(!std::isfinite(double(j))||j<=floor())return false; }
        return true;
    }
    btScalar safeStep(btScalar dt,const TVStack& dv) const
    {
        btScalar step=1;
        for(int b=0;b<materials.size();++b)if(active(materials[b].body))
            for(int e=0;e<materials[b].body->m_tetras.size();++e)
            {
                const auto& t=materials[b].body->m_tetras[e];const auto f=deformation(t);const btScalar j=f.determinant();
                if(!std::isfinite(double(j))||j<=floor())return 0;
                btVector3 delta[4];for(int k=0;k<4;++k)delta[k]=movable(t.m_n[k])?dv[t.m_n[k]->index]*dt:btVector3(0,0,0);
                const auto df=btMatrix3x3(delta[1]-delta[0],delta[2]-delta[0],delta[3]-delta[0]).transpose()*t.m_Dm_inverse;
                const auto relative=f.inverse()*df;
                const btScalar norm=btSqrt(relative[0].length2()+relative[1].length2()+relative[2].length2());
                // Bound all singular values throughout the step, including its interior.
                if(norm>0)step=btMin(step,btScalar(.9)*(1-btPow(floor()/j,btScalar(1.0/3.0)))/norm);
            }
        return step;
    }
    void addScaledExplicitForce(btScalar,TVStack&) override {}
    void addScaledDampingForce(btScalar,TVStack&) override {}
    void addScaledDampingForceDifferential(btScalar,const TVStack&,TVStack&) override {}
    void buildDampingForceDifferentialDiagonal(btScalar,TVStack&) override {}
};
#endif
