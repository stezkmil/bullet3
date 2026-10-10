#ifndef BT_GIMPACT_VERTEX_CACHE_H
#define BT_GIMPACT_VERTEX_CACHE_H

#include "LinearMath/btAlignedObjectArray.h"
#include "LinearMath/btVector3.h"

// Prepare on the calling thread; readers and geometry must remain unchanged until end().
class btGImpactVertexCache
{
	btAlignedObjectArray<btVector3> m_current, m_safe;
	int m_depth;
	bool m_hasSafe = false;

public:
	btGImpactVertexCache() : m_depth(0) {}
	// A cloned manager must never inherit an active query or its geometry.
	btGImpactVertexCache(const btGImpactVertexCache&) : m_depth(0) {}
	btGImpactVertexCache& operator=(const btGImpactVertexCache&)
	{
		m_depth = 0;
		return *this;
	}

	template <class Reconstruct>
	void begin(int count, const Reconstruct& reconstruct)
	{
		begin(count, reconstruct, [](int n, const auto& operation) {
			for (int i = 0; i < n; ++i) operation(i);
		});
	}
	// The dispatcher completes all reconstruction before the snapshot becomes readable.
	template <class Reconstruct, class Execute>
	void begin(int count, const Reconstruct& reconstruct, const Execute& execute)
	{
		if (m_depth == 0 || !m_hasSafe)
		{
			m_hasSafe = true;
			m_current.resize(count);
			m_safe.resize(count);
			execute(count, [&](int i) { reconstruct(i, m_current[i], m_safe[i]); });
		}
		++m_depth;
	}

	template <class Reconstruct>
	void beginCurrent(int count, const Reconstruct& reconstruct)
	{
		beginCurrent(count, reconstruct, [](int n, const auto& operation) {
			for (int i = 0; i < n; ++i) operation(i);
		});
	}
	// Execute must finish all writes before returning; allocation stays on the caller.
	template <class Reconstruct, class Execute>
	void beginCurrent(int count, const Reconstruct& reconstruct, const Execute& execute)
	{
		if (m_depth == 0)
		{
			m_hasSafe = false;
			m_current.resize(count);
			execute(count, [&](int i) { reconstruct(i, m_current[i]); });
		}
		++m_depth;
	}

	void end() { if (m_depth > 0) --m_depth; }
	const btVector3* current(unsigned int index) const { return m_depth ? &m_current[index] : 0; }
	const btVector3* safe(unsigned int index) const { return m_depth && m_hasSafe ? &m_safe[index] : 0; }
};

#endif
