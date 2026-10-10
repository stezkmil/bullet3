// Bullet's GImpact traversal and triangle narrow phase for the VBD collision surface.
#ifndef BT_DEFORMABLE_VBD_COLLISION_H
#define BT_DEFORMABLE_VBD_COLLISION_H
#include "BulletCollision/Gimpact/btGImpactBvh.h"
#include <vector>
#include <array>
#include <algorithm>
#include "btDeformableVbdGpu.h"
#include "btDeformableVbdParallel.h"
#include "BulletCollision/Gimpact/btGImpactVertexCache.h"

using btVbdPairSet = std::vector<GIM_PAIR>;

// Preserve deterministic contact order while avoiding pointer-chasing during large sorts.
inline void btVbdSortCollisionPairs(btVbdPairSet &pairs, std::vector<GIM_PAIR> &scratch)
{
	const auto less = [](const GIM_PAIR &a, const GIM_PAIR &b)
	{ return a.m_index1 != b.m_index1 ? a.m_index1 < b.m_index1 : a.m_index2 < b.m_index2; };
	if (pairs.size() < 512)
	{
		std::sort(pairs.begin(), pairs.end(), less);
		return;
	}
	scratch.resize(pairs.size());
	unsigned int varying[2] = {0, 0};
	for (const auto &p : pairs)
	{
		varying[0] |= unsigned(p.m_index2) ^ unsigned(pairs[0].m_index2);
		varying[1] |= unsigned(p.m_index1) ^ unsigned(pairs[0].m_index1);
	}
	for (int field = 0; field < 2; ++field)
		for (int shift = 0; shift < 32; shift += 8)
		{
			if (!(varying[field] & (255u << shift)))
				continue;
			size_t offsets[256] = {};
			auto digit = [&](const GIM_PAIR &p) { return (((unsigned(field ? p.m_index1 : p.m_index2) ^ 0x80000000u) >> shift) & 255u); };
			for (const auto &p : pairs)
				++offsets[digit(p)];
			size_t total = 0;
			for (auto &offset : offsets)
			{
				size_t count = offset;
				offset = total;
				total += count;
			}
			for (const auto &p : pairs)
				scratch[offsets[digit(p)]++] = p;
			pairs.swap(scratch);
		}
}

// Stable per-chunk bucket offsets preserve the serial pair order without atomic writes.
inline void btVbdSortCollisionPairsParallel(btVbdPairSet &pairs, std::vector<GIM_PAIR> &scratch, int workers)
{
	if (workers <= 1 || pairs.size() < 131072)
	{
		btVbdSortCollisionPairs(pairs, scratch);
		return;
	}
	scratch.resize(pairs.size());
	unsigned varying[2] = {0, 0};
	for (const auto &p : pairs)
	{
		varying[0] |= unsigned(p.m_index2) ^ unsigned(pairs[0].m_index2);
		varying[1] |= unsigned(p.m_index1) ^ unsigned(pairs[0].m_index1);
	}
	constexpr int Bits = 8;
	constexpr unsigned Bins = 1u << Bits, Mask = Bins - 1;
	struct alignas(64) Counts
	{
		std::array<size_t, Bins> values;
	};
	const int jobs = workers <= 1 ? 1 : std::max(1, std::min(workers * 2, int((pairs.size() + 8191) / 8192)));
	std::vector<Counts> counts(jobs);
	auto range = [&](int j) { return std::pair<size_t, size_t>(pairs.size() * size_t(j) / jobs, pairs.size() * size_t(j + 1) / jobs); };
	for (int field = 0; field < 2; ++field)
		for (int shift = 0; shift < 32; shift += Bits)
		{
			if (!(varying[field] & (Mask << shift)))
				continue;
			auto digit = [&](const GIM_PAIR &p) { return ((unsigned(field ? p.m_index1 : p.m_index2) ^ 0x80000000u) >> shift) & Mask; };
			auto count = [&](int job)
			{
				auto &c = counts[job].values;
				c.fill(0);
				auto r = range(job);
				for (size_t i = r.first; i < r.second; ++i)
					++c[digit(pairs[i])];
			};
			btVbdParallelFor(jobs, workers, count, 2, 1);
			size_t total = 0;
			for (unsigned digit = 0; digit < Bins; ++digit)
				for (int job = 0; job < jobs; ++job)
				{
					auto &entry = counts[job].values[digit];
					const auto n = entry;
					entry = total;
					total += n;
				}
			auto scatter = [&](int job)
			{
				auto &c = counts[job].values;
				auto r = range(job);
				for (size_t i = r.first; i < r.second; ++i)
					scratch[c[digit(pairs[i])]++] = pairs[i];
			};
			btVbdParallelFor(jobs, workers, scatter, 2, 1);
			pairs.swap(scratch);
		}
}

class btDeformableVbdCollisionMesh : public btPrimitiveManagerBase
{
	int m_builtCount = 0;
	bool m_ownersValid = false;
	std::vector<int> m_triangleOwners, m_subtreeOwners;
	std::vector<int> m_refitRoots, m_refitTop;
	std::vector<unsigned char> m_triangleChanges;
	void buildRefitSchedule()
	{
		m_refitRoots.clear();
		m_refitTop.clear();
		// Subtrees occupy contiguous ranges; keep the parallel walks cache-friendly.
		for (int node = 0; node < tree.getNodeCount();)
		{
			const int size = tree.isLeafNode(node) ? 1 : tree.getEscapeNodeIndex(node);
			if (size <= 512)
			{
				m_refitRoots.push_back(node);
				node += size;
			}
			else
				m_refitTop.push_back(node++);
		}
	}
	void refitNode(int node)
	{
		btAABB bound;
		if (tree.isLeafNode(node))
			get_primitive_box(tree.getNodeData(node), bound);
		else
		{
			bound.invalidate();
			btAABB child;
			tree.getNodeBound(tree.getLeftNode(node), child);
			bound.merge(child);
			tree.getNodeBound(tree.getRightNode(node), child);
			bound.merge(child);
		}
		tree.setNodeBound(node, bound);
	}

  public:
	std::vector<btPrimitiveTriangle> triangles;
	btGImpactBvh tree;
	btGImpactVertexCache candidatePositions;
	btDeformableVbdCollisionMesh() : tree(this) {}
	btDeformableVbdCollisionMesh(const btDeformableVbdCollisionMesh &) = delete;
	btDeformableVbdCollisionMesh &operator=(const btDeformableVbdCollisionMesh &) = delete;
	btPrimitiveManagerBase *clone() const override
	{
		auto *copy = new btDeformableVbdCollisionMesh();
		copy->triangles = triangles;
		copy->update();
		return copy;
	}
	bool is_trimesh() const override
	{
		return true;
	}
	int get_primitive_count() const override
	{
		return int(triangles.size());
	}
	void get_primitive_triangle(int index, btPrimitiveTriangle &triangle, bool) const override
	{
		triangle = triangles[index];
	}
	bool get_primitive_triangle_safe(int index, btPrimitiveTriangle &triangle) const override
	{
		triangle = triangles[index];
		return true;
	}
	void get_primitive_indices(int index, unsigned int &a, unsigned int &b, unsigned int &c) const override
	{
		a = 3 * index;
		b = a + 1;
		c = a + 2;
	}
	void get_primitive_box(int index, btAABB &box) const override
	{
		const auto &t = triangles[index];
		box.calc_from_triangle_margin(t.m_vertices[0], t.m_vertices[1], t.m_vertices[2], t.m_margin + t.m_discoveryPadding);
	}
	const std::vector<int> &subtreeOwners(const std::vector<int> &owners)
	{
		if (m_ownersValid && m_triangleOwners == owners)
			return m_subtreeOwners;
		m_triangleOwners = owners;
		m_subtreeOwners.assign(triangles.empty() ? 0 : tree.getNodeCount(), (-2147483647 - 1));
		for (int node = int(m_subtreeOwners.size()) - 1; node >= 0; --node)
		{
			if (tree.isLeafNode(node))
				m_subtreeOwners[node] = owners[tree.getNodeData(node)];
			else
			{
				const int left = m_subtreeOwners[tree.getLeftNode(node)], right = m_subtreeOwners[tree.getRightNode(node)];
				if (left == right)
					m_subtreeOwners[node] = left;
			}
		}
		m_ownersValid = true;
		return m_subtreeOwners;
	}
	template <class Position> bool updateTriangles(int count, btScalar margin, btScalar padding, int workers, const Position &position)
	{
		const bool resized = triangles.size() != size_t(count);
		triangles.resize(count);
		const int batch = 256, blocks = (count + batch - 1) / batch;
		m_triangleChanges.resize(blocks);
		btVbdParallelFor(blocks, workers,
						 [&](int block)
						 {
							 bool changed = false;
							 const int end = btMin(count, (block + 1) * batch);
							 for (int k = block * batch; k < end; ++k)
							 {
								 auto &t = triangles[k];
								 bool triangleChanged = resized || t.m_margin != margin || t.m_discoveryPadding != padding;
								 for (int j = 0; j < 3; ++j)
								 {
									 const btVector3 point = position(k, j);
									 triangleChanged = triangleChanged || t.m_vertices[j] != point;
									 t.m_vertices[j] = point;
								 }
								 t.m_margin = margin;
								 t.m_discoveryPadding = padding;
								 if (triangleChanged)
									 t.buildTriPlane();
								 changed = changed || triangleChanged;
							 }
							 m_triangleChanges[block] = changed;
						 });
		for (auto changed : m_triangleChanges)
			if (changed)
				return true;
		return resized;
	}
	void update(int workers = 1)
	{
		if (triangles.empty())
			return;
		if (m_builtCount != int(triangles.size()))
		{
			if (workers > 1 && triangles.size() >= 4096)
				tree.buildSetParallel(workers * 4,
									  [&](int count, const auto &operation) { btVbdParallelFor(count, workers, operation, 2, 1); });
			else
				tree.buildSet();
			buildRefitSchedule();
			m_ownersValid = false;
			m_builtCount = int(triangles.size());
		}
		else if (workers <= 1)
			tree.update();
		else
		{
			btVbdParallelFor(int(m_refitRoots.size()), workers,
							 [&](int i)
							 {
								 const int root = m_refitRoots[i];
								 const int size = tree.isLeafNode(root) ? 1 : tree.getEscapeNodeIndex(root);
								 for (int node = root + size - 1; node >= root; --node)
									 refitNode(node);
							 });
			// All independent subtrees are complete before merging their ancestors.
			for (auto node = m_refitTop.rbegin(); node != m_refitTop.rend(); ++node)
				refitNode(*node);
		}
	}
};
struct btDeformableVbdPairDistance
{
	btVector3 closestSoft, closestRigid;
	btScalar distance2;
	int status = 0;
};
struct btDeformableVbdGuardPosition
{
	btVector3 reference, current, proposed;
	unsigned int stamp = 0, referenceStamp = 0, affectedStamp = 0;
	bool changed = false;
};
// Geometry is refreshed each step; retaining the hierarchy avoids rebuilding static meshes.
struct btDeformableVbdCollisionCache
{
	btDeformableVbdCollisionMesh soft, barrier;
	btDeformableVbdGpu gpu;
	// Only storage persists; each distance and sorting pass overwrites its results.
	std::vector<btDeformableVbdPairDistance> pairDistances;
	std::vector<GIM_PAIR> pairSortScratch;
	// Retain scratch storage; generation stamps invalidate geometry between guard passes.
	std::vector<btDeformableVbdGuardPosition> guardPositions;
	unsigned int guardStamp = 0, guardReferenceStamp = 0;
};
#endif
