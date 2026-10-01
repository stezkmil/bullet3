// Optional, read-only diagnostics for the custom deformable contact path.
#ifndef BT_DEFORMABLE_DIAGNOSTICS_H
#define BT_DEFORMABLE_DIAGNOSTICS_H

#include "btSoftBody.h"
#include <cmath>
#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <mutex>

namespace btDeformableDiagnostics
{
struct Output
{
	FILE* file = nullptr;
	std::mutex mutex;
	Output()
	{
		const char* path = std::getenv("BULLET_DEFORMABLE_DIAGNOSTICS");
		if (!path || !*path) return;
		file = std::fopen(path, "a");
		if (file)
		{
			std::setvbuf(file, nullptr, _IOFBF, 65536);
			std::fprintf(file, "SESSION version=1 scalar_bytes=%zu units=scene strain=dimensionless\n", sizeof(btScalar));
			std::fflush(file);
		}
		else std::fprintf(stderr, "Bullet deformable diagnostics: cannot open %s\n", path);
	}
	~Output() { if (file) std::fclose(file); }
};

inline Output& output()
{
	static Output value;
	return value;
}

struct Context
{
	const void* world;
	long long step;
	double dt;
	int rigidSurfaceSamples = 0;
};

inline Context*& current()
{
	static thread_local Context* value = nullptr;
	return value;
}

inline bool enabled() { return current() != nullptr; }

inline void write(const char* kind, const char* format, ...)
{
	if (!enabled()) return;
	auto& out = output();
	std::lock_guard<std::mutex> lock(out.mutex);
	std::fprintf(out.file, "%s world=%p step=%lld dt=%.9g ", kind, current()->world, current()->step, current()->dt);
	va_list args;
	va_start(args, format);
	std::vfprintf(out.file, format, args);
	va_end(args);
	std::fputc('\n', out.file);
}

class StepScope
{
	Context m_context;
	Context* m_previous;
	std::chrono::steady_clock::time_point m_start;
public:
	StepScope(const void* world, long long step, btScalar dt) : m_context{world, step, double(dt)}, m_previous(current()), m_start(std::chrono::steady_clock::now())
	{
		if (output().file) current() = &m_context;
	}
	~StepScope()
	{
		if (enabled())
		{
			const double elapsed = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - m_start).count();
			write("STEP_TIMING", "elapsed_ms=%.6f", elapsed);
			auto& out = output();
			std::lock_guard<std::mutex> lock(out.mutex);
			std::fflush(out.file);
		}
		current() = m_previous;
	}
};

inline double kinetic(const btSoftBody::Node& node, const btVector3& velocity)
{
	return node.m_im > 0 ? 0.5 * double(velocity.length2()) / double(node.m_im) : 0;
}

inline void bodies(const char* stage, const btAlignedObjectArray<btSoftBody*>& softBodies, bool measureStrain = false)
{
	if (!enabled()) return;
	for (int b = 0; b < softBodies.size(); ++b)
	{
		const btSoftBody& body = *softBodies[b];
		double mass = 0, ke = 0, maxSpeed = 0, maxSplit = 0;
		btVector3 momentum(0, 0, 0);
		int frozen = 0, fixed = 0, constrained = 0, invalid = 0;
		for (int i = 0; i < body.m_nodes.size(); ++i)
		{
			const auto& n = body.m_nodes[i];
			frozen += n.m_frozen > 0;
			fixed += n.m_im <= 0;
			constrained += n.m_constrained != 0;
			const double speed = double(n.m_v.length());
			const double split = double(n.m_splitv.length());
			if (!std::isfinite(speed) || !std::isfinite(split) || !std::isfinite(double(n.m_x.length2()))) ++invalid;
			maxSpeed = btMax(maxSpeed, speed);
			maxSplit = btMax(maxSplit, split);
			if (n.m_im > 0)
			{
				mass += 1.0 / double(n.m_im);
				ke += kinetic(n, n.m_v);
				momentum += n.m_v / n.m_im;
			}
		}
		const btVector3 comVelocity = mass > 0 ? momentum / btScalar(mass) : btVector3(0, 0, 0);
		write("BODY", "stage=%s id=%d ptr=%p active=%d activation=%d nodes=%d frozen=%d fixed=%d constrained=%d invalid=%d mass=%.9g ke=%.9g com_vx=%.9g com_vy=%.9g com_vz=%.9g max_v=%.9g max_split_v=%.9g nn_contacts=%d rigid_contacts=%d anchors=%d friction=%.9g contact_stiffness=%.9g drag=%.9g",
			stage, body.getUserIndex(), static_cast<const void*>(&body), int(body.isActive()), body.getActivationState(), body.m_nodes.size(), frozen, fixed, constrained, invalid,
			mass, ke, double(comVelocity.x()), double(comVelocity.y()), double(comVelocity.z()), maxSpeed, maxSplit,
			body.m_nodeNodeContacts.size(), body.m_nodeRigidContacts.size() + body.m_faceRigidContacts.size(), body.m_deformableAnchors.size(), double(body.m_cfg.kDF), double(body.m_softVsSoftContactStiffness), double(body.m_cfg.drag));
		if (!measureStrain || body.m_tetras.size() == 0) continue;
		double minJ = std::numeric_limits<double>::infinity(), maxStrain = 0, maxTrace = 0, sumStrain = 0, weight = 0;
		int inverted = 0, invalidTets = 0, worstTet = -1;
		for (int i = 0; i < body.m_tetras.size(); ++i)
		{
			const auto& t = body.m_tetras[i];
			const btVector3 a = t.m_n[1]->m_x - t.m_n[0]->m_x;
			const btVector3 b = t.m_n[2]->m_x - t.m_n[0]->m_x;
			const btVector3 c = t.m_n[3]->m_x - t.m_n[0]->m_x;
			const btMatrix3x3 F = btMatrix3x3(a.x(), b.x(), c.x(), a.y(), b.y(), c.y(), a.z(), b.z(), c.z()) * t.m_Dm_inverse;
			const btMatrix3x3 E = (F.transpose() * F - btMatrix3x3::getIdentity()) * btScalar(0.5);
			const double norm = std::sqrt(double(E[0].length2() + E[1].length2() + E[2].length2()));
			const double J = double(F.determinant());
			if (!std::isfinite(norm) || !std::isfinite(J)) { ++invalidTets; continue; }
			minJ = btMin(minJ, J);
			inverted += J <= 0;
			if (norm > maxStrain) { maxStrain = norm; worstTet = i; }
			maxTrace = btMax(maxTrace, std::fabs(double(E[0][0] + E[1][1] + E[2][2])));
			const double volume = std::fabs(double(t.m_element_measure));
			sumStrain += norm * volume;
			weight += volume;
		}
		write("STRAIN", "stage=%s id=%d tets=%d min_J=%.9g inverted=%d invalid=%d max_norm=%.9g mean_norm=%.9g max_abs_trace=%.9g worst_tet=%d",
			stage, body.getUserIndex(), body.m_tetras.size(), minJ, inverted, invalidTets, maxStrain, weight > 0 ? sumStrain / weight : 0, maxTrace, worstTet);
	}
}

struct ContactSummary
{
	int applied = 0, boosted = 0, immobile = 0;
	double energyDelta = 0, positiveEnergy = 0, maxDv = 0, maxDepth = 0;
	int other = -1, node0 = -1, node1 = -1, constrained0 = 0, constrained1 = 0;
	double worstVnBefore = 0, worstVnAfter = 0, worstImpulse = 0, worstBias = 0;

	void add(const btSoftBody::DeformableNodeNodeContact& c, const btVector3& v0, const btVector3& v1, btScalar impulse, btScalar bias)
	{
		++applied;
		const auto& n0 = *c.m_node0;
		const auto& n1 = *c.m_node1;
		boosted += (n0.m_constrained && n0.m_frozen <= 0 && n0.m_im > 0) || (n1.m_constrained && n1.m_frozen <= 0 && n1.m_im > 0);
		immobile += n0.m_frozen > 0 || n1.m_frozen > 0 || n0.m_im <= 0 || n1.m_im <= 0;
		const double delta = kinetic(n0, n0.m_v) + kinetic(n1, n1.m_v) - kinetic(n0, v0) - kinetic(n1, v1);
		energyDelta += delta;
		positiveEnergy += btMax(0.0, delta);
		maxDepth = btMax(maxDepth, -double(c.m_offset));
		const double dv = btMax(double((n0.m_v - v0).length()), double((n1.m_v - v1).length()));
		if (dv > maxDv || applied == 1)
		{
			maxDv = dv;
			other = c.m_colObj->getUserIndex();
			node0 = n0.local_index;
			node1 = n1.local_index;
			constrained0 = n0.m_constrained;
			constrained1 = n1.m_constrained;
			worstVnBefore = double((v0 - v1).dot(c.m_normal));
			worstVnAfter = double((n0.m_v - n1.m_v).dot(c.m_normal));
			worstImpulse = double(impulse);
			worstBias = double(bias);
		}
	}

	void report(const btSoftBody& body) const
	{
		if (body.m_nodeNodeContacts.size() == 0) return;
		write("SOFT_CONTACT", "id=%d contacts=%d applied=%d boosted=%d immobile=%d delta_ke=%.9g positive_delta_ke=%.9g max_dv=%.9g max_depth=%.9g worst_other=%d node0=%d node1=%d constrained0=%d constrained1=%d vn_before=%.9g vn_after=%.9g impulse=%.9g bias=%.9g",
			body.getUserIndex(), body.m_nodeNodeContacts.size(), applied, boosted, immobile, energyDelta, positiveEnergy, maxDv, maxDepth, other, node0, node1, constrained0, constrained1, worstVnBefore, worstVnAfter, worstImpulse, worstBias);
	}
};
}
#endif
