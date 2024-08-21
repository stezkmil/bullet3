#include "btSoftBody.h"
#include "../BulletCollision/CollisionShapes/btTriangleIndexVertexArray.h"

template <typename ShapeType> class btSoftBodyWithCollisionShape : public btSoftBody
{
public:
	struct btVertexToTetraMapping
	{
		unsigned vertexToTetra;               // Index of first tetra index to which the collision mesh vertex is attached
		btVector4 baryCoordInTetra;           // Barycentric coordiante of collision shape vertex in that tetra
		btVector4 baryCoordNormalEndInTetra;  // Barycentric coordiante of collision shape vertex normal end in that tetra
	};

private:
	// Uses technique described here https://nvidia-omniverse.github.io/PhysX/physx/5.4.0/docs/SoftBodies.html?highlight=soft%20bodies#
	// (both simulation shape and collision shape are used in a soft body)
	bool hasCollisionShape = false;
	std::vector<btVector3> verticesCollisionShape;    // TODO shared pointers to avoid a vector copy
	std::vector<btVector3> normalsCollisionShape;     // TODO shared pointers to avoid a vector copy
	std::vector<unsigned int> indicesCollisionShape;  // TODO shared pointers to avoid a vector copy
	std::vector<btVertexToTetraMapping> collisionShapeVertexToSimTetra;

	btTriangleIndexVertexArray triangleIndexVertexArray;
	btIndexedMesh indexedMesh;
	std::unique_ptr<ShapeType> collisionMeshShape;

public:
	btSoftBodyWithCollisionShape(btSoftBodyWorldInfo* worldInfo, int node_count, const btVector3* simulationShapeCoordinates, const btScalar* m,
								 const std::vector<std::array<int, 4>>& simulationShapeTetras, bool createLinks, bool hasCollisionShape, const std::vector<unsigned int>& indicesCollisionShape,
								 const std::vector<btVector3>& verticesCollisionShape, const std::vector<btVector3>& normalsCollisionShape,
								 const std::vector<btVertexToTetraMapping>& collisionShapeVertexToSimTetra) : hasCollisionShape(hasCollisionShape), indicesCollisionShape(indicesCollisionShape), verticesCollisionShape(verticesCollisionShape), normalsCollisionShape(normalsCollisionShape), collisionShapeVertexToSimTetra(collisionShapeVertexToSimTetra), btSoftBody(worldInfo, node_count, simulationShapeCoordinates, m)
	{
		if (hasCollisionShape)
		{
			indexedMesh.m_vertexBase = reinterpret_cast<const unsigned char*>(verticesCollisionShape.data());
			indexedMesh.m_vertexStride = sizeof(btScalar) * 3;
			indexedMesh.m_numVertices = static_cast<int>(verticesCollisionShape.size());
			indexedMesh.m_triangleIndexBase = reinterpret_cast<const unsigned char*>(indicesCollisionShape.data());
			indexedMesh.m_triangleIndexStride = sizeof(unsigned int) * 3;
			indexedMesh.m_numTriangles = static_cast<int>(indicesCollisionShape.size() / 3);
			indexedMesh.m_indexType = PHY_INTEGER;
			indexedMesh.m_vertexType = PHY_FLOAT;

			triangleIndexVertexArray.addIndexedMesh(indexedMesh, PHY_INTEGER);

			collisionMeshShape = std::make_unique<ShapeType>(&triangleIndexVertexArray);

			if (collisionMeshShape->getShapeType() != INVALID_SHAPE_PROXYTYPE)
			{
				collisionMeshShape->setMargin(getCollisionShape()->getMargin());
			}
		}
	}

	virtual ~btSoftBodyWithCollisionShape() {}

	bool HasCollisionShape() const
	{
		return hasCollisionShape;
	}

	const std::vector<unsigned int>& GetCollisionShapeIndices() const
	{
		return indicesCollisionShape;
	}

	const std::vector<btVector3>& GetCollisionShapeVertices() const
	{
		return verticesCollisionShape;
	}

	const std::vector<btVector3>& GetCollisionShapeNormals() const
	{
		return normalsCollisionShape;
	}

	const std::vector<btVertexToTetraMapping>& GetCollisionShapeVertexToSimTetra() const
	{
		return collisionShapeVertexToSimTetra;
	}
};