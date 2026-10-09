// Bullet's GImpact traversal and triangle narrow phase for the VBD collision surface.
#ifndef BT_DEFORMABLE_VBD_COLLISION_H
#define BT_DEFORMABLE_VBD_COLLISION_H
#include "BulletCollision/Gimpact/btGImpactBvh.h"
#include <vector>
#include "btDeformableVbdParallel.h"
#include "BulletCollision/Gimpact/btGImpactVertexCache.h"

class btDeformableVbdCollisionMesh : public btPrimitiveManagerBase
{
	int m_builtCount = 0;
	bool m_ownersValid = false;
	std::vector<int> m_triangleOwners, m_subtreeOwners;
	std::vector<int> m_refitRoots, m_refitTop;
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
	void update(int workers = 1)
	{
		if (triangles.empty())
			return;
		if (m_builtCount != int(triangles.size()))
		{
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
// Geometry is refreshed each step; retaining the hierarchy avoids rebuilding static meshes.
struct btDeformableVbdCollisionCache
{
	btDeformableVbdCollisionMesh soft, barrier;
};
#endif
