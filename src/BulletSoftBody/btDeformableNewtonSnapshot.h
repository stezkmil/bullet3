#ifndef BT_DEFORMABLE_NEWTON_SNAPSHOT_H
#define BT_DEFORMABLE_NEWTON_SNAPSHOT_H

#include "btDeformableBodySolver.h"
#include "btDeformableContactForce.h"
#include "btDeformableVolumeBarrierForce.h"
#include "btDeformableLinearElasticityForce.h"
#include "btDeformableNodalForce.h"
#include "btDeformableGravityForce.h"
#include <cstdio>
#include <memory>
#include <vector>

// Diagnostic format for the implicit unconstrained KKT path. No pointers or
// collision objects are persisted; contacts already contain their surface maps.
class btDeformableNewtonSnapshot
{
	typedef btAlignedObjectArray<btVector3> Vectors;
	btSoftBodyWorldInfo info;
	std::vector<std::unique_ptr<btSoftBody> > bodies;
	std::vector<std::unique_ptr<btDeformableLagrangianForce> > forces;
	FILE* file = nullptr;
	bool reading = false, good = true;
	template<class T> void value(T& v)
	{
		if (!good) return;
		good = (reading ? std::fread(&v,sizeof(T),1,file) : std::fwrite(&v,sizeof(T),1,file)) == 1;
	}
	void vector(btVector3& v) { for(int i=0;i<3;++i) value(v[i]); }
	void matrix(btMatrix3x3& m) { for(int i=0;i<3;++i) vector(m[i]); }
	void scratch(btSoftBody::TetraScratch& s)
	{ matrix(s.m_F); value(s.m_trace); value(s.m_J); matrix(s.m_cofF); matrix(s.m_corotation); }
	int count(int existing, int maximum = 1000000)
	{ int n=existing; value(n); if(n<0 || n>maximum) { good=false; return 0; } return n; }
	void vectors(Vectors& v)
	{ int n=count(v.size()); if(reading)v.resize(n); for(int i=0;i<n && good;++i)vector(v[i]); }
	int bodyIndex(btDeformableBodySolver& s, btSoftBody* b)
	{ for(int i=0;i<s.m_softBodies.size();++i)if(s.m_softBodies[i]==b)return i; good=false; return -1; }
	void bodyReference(btDeformableBodySolver& s, btSoftBody*& b)
	{
		int index=reading?0:bodyIndex(s,b); value(index);
		if(index<0 || index>=s.m_softBodies.size()){good=false;return;}
		if(reading)b=s.m_softBodies[index];
	}
	void nodeReference(btDeformableBodySolver& s, btSoftBody::Node*& n)
	{
		int index=reading?0:n->index; value(index);
		if(index<0 || index>=s.m_objective->m_nodes.size()){good=false;return;}
		if(reading)n=s.m_objective->m_nodes[index];
	}
	bool supported(btDeformableBodySolver& s)
	{
		if(!s.m_implicit || s.m_useProjection || s.m_objective->m_projection.m_lagrangeMultipliers.size() ||
			s.m_objective->m_preconditioner!=s.m_objective->m_KKTPreconditioner)return false;
		for(int i=0;i<s.m_objective->m_lf.size();++i)
		{
			const auto type=s.m_objective->m_lf[i]->getForceType();
			if(type!=BT_LINEAR_ELASTICITY_FORCE && type!=BT_GRAVITY_FORCE && type!=BT_NODAL_FORCE &&
				type!=BT_CONTACT_FORCE && type!=BT_VOLUME_BARRIER_FORCE)return false;
		}
		return true;
	}
	void transfer(btDeformableBodySolver& s)
	{
		int magic=0x4e575431, scalarBytes=sizeof(btScalar);
		value(magic); value(scalarBytes);
		if(magic!=0x4e575431 || scalarBytes!=sizeof(btScalar)){good=false;return;}
		value(dt); value(step); value(baselineEnergy); value(slope); value(initialScale);
		value(assembledElastic);value(rotationCorrection);
		value(linearRecoveryBefore);
		int trials=count(scales.size(),128); if(reading){scales.resize(trials);energies.resize(trials);}
		for(int i=0;i<trials;++i){value(scales[i]);value(energies[i]);}
		int nb=count(s.m_softBodies.size(),10000);
		std::vector<int> activations;
		for(int b=0;b<nb && good;++b)
		{
			int nn=count(reading?0:s.m_softBodies[b]->m_nodes.size());
			if(reading)
			{
				Vectors positions; positions.resize(nn,btVector3(0,0,0));
				btAlignedObjectArray<btScalar> masses; masses.resize(nn,btScalar(1));
				bodies.emplace_back(new btSoftBody(&info,nn,positions.size()?&positions[0]:nullptr,masses.size()?&masses[0]:nullptr));
				s.m_softBodies.push_back(bodies.back().get());
			}
			auto& body=*s.m_softBodies[b];
			int id=body.getUserIndex(), activation=body.getActivationState(), flags=body.getCollisionFlags();
			value(id);value(activation);value(flags);value(body.m_gravityFactor);
			activations.push_back(activation);
			if(reading){body.setUserIndex(id);body.forceActivationState(activation);body.setCollisionFlags(flags);}
			for(int n=0;n<nn;++n)
			{
				auto& node=body.m_nodes[n];
				vector(node.m_x);vector(node.m_q);vector(node.m_v);vector(node.m_vn);vector(node.m_splitv);
				value(node.m_im);value(node.m_frozen);
			}
			int nt=count(body.m_tetras.size());
			if(reading){body.m_tetraScratches.resize(nt);body.m_tetraScratchesTn.resize(nt);}
			bool completeTn = true;
			for(int t=0;t<nt && good;++t)
			{
				int nodes[4]={};
				for(int n=0;n<4;++n)
				{
					if(!reading)nodes[n]=int(body.m_tetras[t].m_n[n]-&body.m_nodes[0]);
					value(nodes[n]);if(nodes[n]<0 || nodes[n]>=nn)good=false;
				}
				if(!good)return;
				if(reading)body.appendTetra(nodes[0],nodes[1],nodes[2],nodes[3]);
				auto& tet=body.m_tetras[t];
				matrix(tet.m_Dm_inverse);matrix(tet.m_F);value(tet.m_element_measure);value(tet.m_rv);
				for(int j=0;j<3;++j)for(int k=0;k<4;++k)value(tet.m_P_inv[j][k]);
				scratch(body.m_tetraScratches[t]);
				// Preserve the frozen damping state; direct callers may have no predictor.
				int hasTn=reading?0:int(body.m_tetraScratchesTn.size()>t);value(hasTn);
				if (hasTn)
					scratch(body.m_tetraScratchesTn[t]);
				else
					completeTn = false;
			}
			if (reading && !completeTn)
				body.m_tetraScratchesTn.clear();
		}
		if(!good)return;
		if(reading)
		{
			const auto all=s.m_softBodies; s.setImplicit(true);s.m_useProjection=false;s.reinitialize(all,dt);
			s.setPreconditioner(btDeformableBackwardEulerObjective::KKT_preconditioner);
		}
		value(s.m_newtonIteration);value(s.m_maxNewtonIterations);value(s.m_newtonTolerance);
		value(s.m_contactWeightedTarget);value(s.m_implicitRecoveryUsed);value(s.m_lineSearch);
		vectors(s.m_dv);vectors(s.m_backupVelocity);vectors(s.m_ddv);vectors(s.m_residual);
		if(s.m_dv.size()!=s.m_numNodes || s.m_backupVelocity.size()!=s.m_numNodes || s.m_ddv.size()!=s.m_numNodes || s.m_residual.size()!=s.m_numNodes){good=false;return;}
		int nf=count(s.m_objective->m_lf.size(),10000);
		for(int f=0;f<nf && good;++f)
		{
			int type=reading?0:int(s.m_objective->m_lf[f]->getForceType());value(type);
			btDeformableLagrangianForce* force=reading?nullptr:s.m_objective->m_lf[f];
			if(reading)
			{
				switch(type)
				{
				case BT_LINEAR_ELASTICITY_FORCE: force=new btDeformableLinearElasticityForce();break;
				case BT_GRAVITY_FORCE: force=new btDeformableGravityForce(btVector3(0,0,0));break;
				case BT_NODAL_FORCE: force=new btDeformableNodalForce(nullptr,btAlignedObjectArray<int>(),btVector3(0,0,0));break;
				case BT_CONTACT_FORCE: force=new btDeformableContactForce(dt);break;
				case BT_VOLUME_BARRIER_FORCE: force=new btDeformableVolumeBarrierForce();break;
				default:good=false;return;
				}
				forces.emplace_back(force);s.m_objective->m_lf.push_back(force);
			}
			int countBodies=count(force->m_softBodies.size(),nb);
			for(int b=0;b<countBodies;++b){btSoftBody* body=reading?nullptr:force->m_softBodies[b];bodyReference(s,body);if(reading && good)force->m_softBodies.push_back(body);}
			if(type==BT_LINEAR_ELASTICITY_FORCE)
			{
				auto& e=*static_cast<btDeformableLinearElasticityForce*>(force);
				value(e.m_mu);value(e.m_lambda);value(e.m_damping_alpha);value(e.m_damping_beta);
				if(reading)e.updateYoungsModulusAndPoissonRatio();
			}
			else if(type==BT_GRAVITY_FORCE)vector(static_cast<btDeformableGravityForce*>(force)->m_gravity);
			else if(type==BT_NODAL_FORCE)
			{
				auto& n=*static_cast<btDeformableNodalForce*>(force);btVector3 total=n.totalForce();vector(total);
				btAlignedObjectArray<int> indices=n.nodeIndices();int ni=count(indices.size());if(reading)indices.resize(ni);
				for(int i=0;i<ni;++i)value(indices[i]);
				if(reading){n.setNodeIndices(indices);n.setForce(total);}
			}
			else if(type==BT_VOLUME_BARRIER_FORCE)
			{
				auto& m=static_cast<btDeformableVolumeBarrierForce*>(force)->materials;int nm=count(m.size(),nb);if(reading)m.resize(nm);
				for(int i=0;i<nm;++i){bodyReference(s,m[i].body);value(m[i].bulk);}
			}
			else if(type==BT_CONTACT_FORCE)
			{
				auto& c=*static_cast<btDeformableContactForce*>(force);value(c.dt);int nc=count(c.contacts.size());if(reading)c.contacts.resize(nc);
				for(int i=0;i<nc && good;++i)
				{
					auto& contact=c.contacts[i];vector(contact.normal);vector(contact.tangentImpulse);
					value(contact.gap);value(contact.friction);value(contact.rho);value(contact.tangentRho);value(contact.normalImpulse);
					int nn=count(contact.nodes.size(),s.m_numNodes);if(reading)contact.nodes.resize(nn);
					for(int n=0;n<nn;++n){nodeReference(s,contact.nodes[n].node);matrix(contact.nodes[n].jacobian);}
				}
			}
		}
		if(reading)
		{
			s.m_objective->m_implicitConstraintDv=s.m_dv;
			for(int b=0;b<nb;++b)s.m_softBodies[b]->forceActivationState(activations[b]);
		}
	}
public:
	// Audit the numerical operator separately from the long Krylov recurrence.
	Vectors operatorProbe, preconditionerProbe;
	static void probe(btDeformableBodySolver& s, Vectors& product, Vectors& conditioned)
	{
		auto& objective=*s.m_objective;
		for(int f=0;f<objective.m_lf.size();++f)objective.m_lf[f]->prepareImplicitForceDifferential(objective.m_dt);
		objective.m_KKTPreconditioner->reinitialize(true);
		product.resize(s.m_numNodes);conditioned.resize(s.m_numNodes);
		objective.multiply(s.m_residual,product);
		for(int n=0;n<s.m_numNodes;++n)conditioned[n]=objective.m_KKTPreconditioner->applyInverseNodeBlock(n,s.m_residual[n]);
		for(int f=0;f<objective.m_lf.size();++f)objective.m_lf[f]->finishImplicitForceDifferential();
	}
	btScalar dt=0, baselineEnergy=0, slope=0, initialScale=1;
	long long step=-1;
	int assembledElastic=1, rotationCorrection=1;
	bool linearRecoveryBefore=false;
	btAlignedObjectArray<btScalar> scales, energies;
	btDeformableBodySolver solver;
	bool save(const char* path, btDeformableBodySolver& source)
	{
		if(!supported(source))return false;
		// TODO: Remove all environment-variable lookups across this feature before merging the feature branch.
		// The final implementation must not read environment variables.
		const char* assembly=std::getenv("BULLET_DEFORMABLE_ASSEMBLED_ELASTIC");
		const char* rotation=std::getenv("BULLET_DEFORMABLE_ROTATION_CORRECTION");
		assembledElastic=!(assembly && assembly[0]=='0');rotationCorrection=!(rotation && rotation[0]=='0');
		reading=false;good=true;file=std::fopen(path,"wb");if(!file)return false;
		transfer(source);
		if(good) { probe(source,operatorProbe,preconditionerProbe); vectors(operatorProbe);vectors(preconditionerProbe);value(source.m_objective->m_contactCoarseEnabled); }
		if(std::fclose(file)!=0)good=false;file=nullptr;
		return good;
	}
	bool load(const char* path)
	{
		if(!bodies.empty() || !forces.empty())return false;
		reading=true;good=true;file=std::fopen(path,"rb");if(!file)return false;
		transfer(solver);
		if(good)
		{
			const int next=std::fgetc(file);
			if(next!=EOF)
			{
				std::ungetc(next,file);vectors(operatorProbe);vectors(preconditionerProbe);
				if(operatorProbe.size()!=solver.m_numNodes || preconditionerProbe.size()!=solver.m_numNodes)good=false;
				const int flag=std::fgetc(file);
				if(flag!=EOF) { std::ungetc(flag,file);value(solver.m_objective->m_contactCoarseEnabled);if(std::fgetc(file)!=EOF)good=false; }
			}
		}
		std::fclose(file);file=nullptr;
		return good;
	}
	btScalar energy() { return solver.m_objective->totalEnergy(dt)+solver.kineticEnergy(); }
	btScalar weightedResidual()
	{
		Vectors r;r.resize(solver.m_numNodes,btVector3(0,0,0));
		for(int n=0;n<r.size();++n)
		{
			const auto* node=solver.m_objective->m_nodes[n];
			if(node->m_im>0 && node->m_frozen<=0)r[n]=-solver.m_dv[n]/node->m_im;
		}
		solver.m_objective->computeResidual(dt,r);
		btScalar squared=0;
		for(int n=0;n<r.size();++n)squared+=r[n].dot(solver.m_objective->m_KKTPreconditioner->applyInverseNodeBlock(n,r[n]));
		return btSqrt(squared);
	}
	const Vectors& direction() const { return solver.m_ddv; }
	const Vectors& residual() const { return solver.m_residual; }
	void recomputeDirection() { solver.m_implicitRecoveryUsed=linearRecoveryBefore;for(int n=0;n<solver.m_ddv.size();++n)solver.m_ddv[n].setZero();solver.computeStep(solver.m_ddv,solver.m_residual); }
};
#endif
