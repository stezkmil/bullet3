#include <btBulletDynamicsCommon.h>
#include <BulletSoftBody/btSoftBody.h>
#include <BulletSoftBody/btSoftBodyHelpers.h>
#include <BulletSoftBody/btDeformableMultiBodyDynamicsWorld.h>
#include <BulletSoftBody/btSoftBodyRigidBodyCollisionConfiguration.h>
#include <BulletCollision/Gimpact/btGImpactShape.h>
#include <BulletCollision/Gimpact/btGImpactCollisionAlgorithm.h>

#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <array>
#include <algorithm>
#include <memory>
#include <cmath>

// Define Vector3 alias to match ported code structure
using Vector3 = btVector3;
using Index = unsigned int;

const float transformationCompareEpsilon = 1e-6f;

struct btVector3Hash
{
    std::size_t operator()(const btVector3& k) const
    {
        size_t h1 = std::hash<float>{}(k.getX());
        size_t h2 = std::hash<float>{}(k.getY());
        size_t h3 = std::hash<float>{}(k.getZ());
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

// --- Helper Functions Ported/Adapted from Bullet3Utils and Bullet3DynamicsForCDModule ---

btVector3 ApplyBarycentricCoordinatesInTetra(
    const btVector4& baryCoord, const btVector3& a, const btVector3& b, const btVector3& c, const btVector3& d
)
{
    return baryCoord.x() * a + baryCoord.y() * b + baryCoord.z() * c + baryCoord.w() * d;
}

btVector4 GetBarycentricCoordinatesInTetra(
    const btVector3& p, const btVector3& a, const btVector3& b, const btVector3& c, const btVector3& d
)
{
    btVector3 vap = p - d;
    btVector3 vad = a - d;
    btVector3 vbd = b - d;
    btVector3 vcd = c - d;
    
    btMatrix3x3 m(vad.x(), vbd.x(), vcd.x(),
                  vad.y(), vbd.y(), vcd.y(),
                  vad.z(), vbd.z(), vcd.z());
    
    // We want to solve M * |u v w|^T = vap
    // So |u v w|^T = M_inv * vap
    
    // Invert. Note: Safe inverse might be needed if tetra is degenerate.
    btMatrix3x3 mInv = m.inverse();
    btVector3 uvw = mInv * vap;
    
    float u = uvw.x();
    float v = uvw.y();
    float w = uvw.z();
    float t = 1.0f - u - v - w;
    
    return btVector4(u, v, w, t);
}

int findOrAddVoxelVertex(
    const btVector3& vertex, std::vector<btVector3>& vertices, std::unordered_map<btVector3, int, btVector3Hash>& vertexMap
)
{
    auto it = vertexMap.find(vertex);
    if (it != vertexMap.end())
    {
        return it->second;
    }
    else
    {
        int newIndex = static_cast<int>(vertices.size());
        vertices.emplace_back(vertex);
        vertexMap[vertex] = newIndex;
        return newIndex;
    }
}

bool isPointInsideTetrahedron(
    const btVector3& point, const btVector3& v0, const btVector3& v1, const btVector3& v2, const btVector3& v3
)
{
    auto inside = [](const btVector3& point, const btVector3& v0, const btVector3& v1, const btVector3& v2, const btVector3& v3)
    {
        auto normal = (v1 - v0).cross(v2 - v0);
        auto planed = -normal.dot(v0);
        if (normal.dot(v3) + planed > 0.0)
        {
            normal = -normal;
            planed = -planed;
        }
        auto dist = normal.dot(point) + planed;
        return dist < transformationCompareEpsilon;
    };
    return inside(point, v1, v0, v2, v3) && inside(point, v1, v2, v3, v0) && inside(point, v3, v2, v0, v1) && inside(point, v3, v0, v1, v2);
}

void addTetrahedraForCube(
    const btVector3& v000, const btVector3& v100, const btVector3& v010, const btVector3& v001,
    const btVector3& v110, const btVector3& v101, const btVector3& v011, const btVector3& v111,
    int vertexIndex,
    const std::vector<int>& triIndices,
    const btVector3& vertex,
    std::vector<btVector3>& simVertices,
    std::vector<std::array<int, 4>>& simTetraIndices,
    std::vector<int>& vertexToTet,
    std::vector<std::vector<int>>& triToTet,
    std::unordered_map<btVector3, int, btVector3Hash>& voxelVertexMap,
    bool& cubeAdded
)
{
    if (!cubeAdded)
    {
        int index000 = findOrAddVoxelVertex(v000, simVertices, voxelVertexMap);
        int index100 = findOrAddVoxelVertex(v100, simVertices, voxelVertexMap);
        int index010 = findOrAddVoxelVertex(v010, simVertices, voxelVertexMap);
        int index001 = findOrAddVoxelVertex(v001, simVertices, voxelVertexMap);
        int index110 = findOrAddVoxelVertex(v110, simVertices, voxelVertexMap);
        int index101 = findOrAddVoxelVertex(v101, simVertices, voxelVertexMap);
        int index011 = findOrAddVoxelVertex(v011, simVertices, voxelVertexMap);
        int index111 = findOrAddVoxelVertex(v111, simVertices, voxelVertexMap);

        simTetraIndices.push_back({ index000, index100, index010, index001 });
        simTetraIndices.push_back({ index110, index111, index010, index100 });
        simTetraIndices.push_back({ index111, index001, index011, index010 });
        simTetraIndices.push_back({ index101, index001, index111, index100 });
        simTetraIndices.push_back({ index100, index010, index001, index111 });

        cubeAdded = true;
    }

    int tetraIndex1 = static_cast<int>(simTetraIndices.size() - 5);
    int tetraIndex2 = static_cast<int>(simTetraIndices.size() - 4);
    int tetraIndex3 = static_cast<int>(simTetraIndices.size() - 3);
    int tetraIndex4 = static_cast<int>(simTetraIndices.size() - 2);
    int tetraIndex5 = static_cast<int>(simTetraIndices.size() - 1);

    if (vertexIndex >= 0)
    {
        if (isPointInsideTetrahedron(vertex, v000, v100, v010, v001)) vertexToTet[vertexIndex] = tetraIndex1;
        else if (isPointInsideTetrahedron(vertex, v110, v111, v010, v100)) vertexToTet[vertexIndex] = tetraIndex2;
        else if (isPointInsideTetrahedron(vertex, v111, v001, v011, v010)) vertexToTet[vertexIndex] = tetraIndex3;
        else if (isPointInsideTetrahedron(vertex, v101, v001, v111, v100)) vertexToTet[vertexIndex] = tetraIndex4;
        else if (isPointInsideTetrahedron(vertex, v100, v010, v001, v111)) vertexToTet[vertexIndex] = tetraIndex5;
        else {
             // Fallback
             vertexToTet[vertexIndex] = tetraIndex1; 
        }
    }
    else
    {
        for (auto triIndex : triIndices)
        {
            triToTet[triIndex].push_back(tetraIndex1);
            triToTet[triIndex].push_back(tetraIndex2);
            triToTet[triIndex].push_back(tetraIndex3);
            triToTet[triIndex].push_back(tetraIndex4);
            triToTet[triIndex].push_back(tetraIndex5);
        }
    }
}

// Helper AABB class since we don't have VRUT::AABB
struct SimpleAABB {
    btVector3 m_min, m_max;
    SimpleAABB(const btVector3& _min, const btVector3& _max) : m_min(_min), m_max(_max) {}
    bool Intersects(const btVector3& v0, const btVector3& v1, const btVector3& v2) {
        btVector3 triMin = v0, triMax = v0;
        triMin.setMin(v1); triMin.setMin(v2);
        triMax.setMax(v1); triMax.setMax(v2);
        if (triMax.x() < m_min.x() || triMin.x() > m_max.x()) return false;
        if (triMax.y() < m_min.y() || triMin.y() > m_max.y()) return false;
        if (triMax.z() < m_min.z() || triMin.z() > m_max.z()) return false;
        return true; // Approximation
    }
};

auto querySoftCollisionMeshBVH(
    const btBvhTriangleMeshShape& bvhTriMeshShape, const btVector3& aabbMin, const btVector3& aabbMax
)
{
    class TriangleCallback: public btTriangleCallback
    {
        const btVector3& aabbMin;
        const btVector3& aabbMax;
        std::map<int, btVector3> overlappingVerticesMap;
        std::vector<int> overlappingTriangles;
        SimpleAABB vrutAABB;

    public:
        TriangleCallback(const btVector3& aabbMin, const btVector3& aabbMax) :
            btTriangleCallback(true), aabbMin(aabbMin), aabbMax(aabbMax), vrutAABB(aabbMin, aabbMax)
        {}

        virtual void processTriangle(btVector3* triangle, int partId, int triangleIndex) override {}

        virtual void processTriangleEx(btVector3* triangle, int* triangleVertexIndex, int partId, int triangleIndex) override
        {
            for (int i = 0; i < 3; ++i)
            {
                auto vert = triangle[i];
                if (vert.x() >= aabbMin.x() && vert.x() <= aabbMax.x() && vert.y() >= aabbMin.y() && vert.y() <= aabbMax.y()
                    && vert.z() >= aabbMin.z() && vert.z() <= aabbMax.z())
                {
                    overlappingVerticesMap.insert({ triangleVertexIndex[i], vert });
                }
            }
            if (vrutAABB.Intersects(triangle[0], triangle[1], triangle[2]))
                overlappingTriangles.push_back(triangleIndex);
        }

        const std::map<int, btVector3>& GetOverlappingVerticesMap() const { return overlappingVerticesMap; }
        std::vector<int> GetOverlappingTriangles() const { return overlappingTriangles; }
    };

    TriangleCallback triangleCallback(aabbMin, aabbMax);
    bvhTriMeshShape.processAllTriangles(&triangleCallback, aabbMin, aabbMax);

    return std::make_pair(triangleCallback.GetOverlappingVerticesMap(), triangleCallback.GetOverlappingTriangles());
}

struct VoxelGridCollector
{
    std::vector<btVector3>& vertices;
    std::vector<std::array<int, 4>>& tetraIndices;
    std::vector<int>& vertexToTet;
    std::vector<std::vector<int>>& triToTet;
};

void createVoxelGrid(
    const btBvhTriangleMeshShape& collisionShape,
    const btVector3& collisionShapeAabbMin,
    const btVector3& collisionShapeAabbMax,
    const btScalar maxVoxelSize,
    const btVector3& scale,
    uint32_t& voxelsCreated,
    VoxelGridCollector* collector // Changed optional to pointer for simplicity/compatibility
)
{
    btVector3 range = collisionShapeAabbMax - collisionShapeAabbMin;
    int gridX = static_cast<int>(ceil(range.x() / maxVoxelSize));
    int gridY = static_cast<int>(ceil(range.y() / maxVoxelSize));
    int gridZ = static_cast<int>(ceil(range.z() / maxVoxelSize));

    std::unordered_map<btVector3, int, btVector3Hash> voxelVertexMap;

    for (int x = 0; x < gridX; ++x)
    {
        for (int y = 0; y < gridY; ++y)
        {
            for (int z = 0; z < gridZ; ++z)
            {
                btVector3 aabbMin = collisionShapeAabbMin + btVector3(x * maxVoxelSize, y * maxVoxelSize, z * maxVoxelSize);
                btVector3 aabbMax = aabbMin + btVector3(maxVoxelSize, maxVoxelSize, maxVoxelSize);

                auto overlap = querySoftCollisionMeshBVH(collisionShape, aabbMin, aabbMax);
                auto& vertices = overlap.first;
                auto& triangles = overlap.second;

                if (!vertices.empty() || !triangles.empty())
                {
                    voxelsCreated++;
                    if (collector)
                    {
                        bool cubeAdded = false;
                        btVector3 v000 = aabbMin;
                        btVector3 v100 = aabbMin + btVector3(maxVoxelSize, 0, 0);
                        btVector3 v010 = aabbMin + btVector3(0, maxVoxelSize, 0);
                        btVector3 v001 = aabbMin + btVector3(0, 0, maxVoxelSize);
                        btVector3 v110 = aabbMin + btVector3(maxVoxelSize, maxVoxelSize, 0);
                        btVector3 v101 = aabbMin + btVector3(maxVoxelSize, 0, maxVoxelSize);
                        btVector3 v011 = aabbMin + btVector3(0, maxVoxelSize, maxVoxelSize);
                        btVector3 v111 = aabbMin + btVector3(maxVoxelSize, maxVoxelSize, maxVoxelSize);

                        std::vector<int> emptyTriIndices; // Needed for vertex call
                        for (auto& vertPair : vertices)
                        {
                            addTetrahedraForCube(v000, v100, v010, v001, v110, v101, v011, v111, vertPair.first, emptyTriIndices, vertPair.second,
                                collector->vertices, collector->tetraIndices, collector->vertexToTet, collector->triToTet, voxelVertexMap, cubeAdded);
                        }
                        
                         if (!triangles.empty()) {
                            addTetrahedraForCube(v000, v100, v010, v001, v110, v101, v011, v111, -1, triangles, btVector3(0,0,0),
                                collector->vertices, collector->tetraIndices, collector->vertexToTet, collector->triToTet, voxelVertexMap, cubeAdded);
                        }
                    }
                }
            }
        }
    }
}

btScalar findVoxelSizeForVoxelGridWithMaxVoxels(
    const btBvhTriangleMeshShape& collisionShape,
    const btVector3& collisionShapeAabbMin,
    const btVector3& collisionShapeAabbMax,
    const uint32_t maxVoxelCount,
    const btVector3& scale
)
{
    btVector3 range = collisionShapeAabbMax - collisionShapeAabbMin;
    btScalar volume = range.x() * range.y() * range.z();
    btScalar voxelSize = pow(volume / maxVoxelCount, 1.0f / 3.0f);
    
    // Iterative adjustment
    for (int i = 0; i < 10; ++i)
    {
        uint32_t voxelsCreated = 0;
        createVoxelGrid(collisionShape, collisionShapeAabbMin, collisionShapeAabbMax, voxelSize, scale, voxelsCreated, nullptr);
        if (voxelsCreated == 0) return voxelSize; // Should not happen
        float ratio = (float)maxVoxelCount / voxelsCreated;
        voxelSize *= pow(ratio, 1.0f / 3.0f);
    }
    return voxelSize;
}

void createVoxelTetrahedronMesh(
    std::shared_ptr<std::vector<Vector3>>& verticesCollisionMesh,
    std::shared_ptr<std::vector<Vector3>>& normalsCollisionMesh,
    std::shared_ptr<std::vector<Index>>& indicesCollisionMesh,
    std::vector<btVector3>& simulationVerticesForSoftBody,
    std::vector<std::array<int, 4>>& simulationTetrasForSoftBody,
    bool hasCollisionShape,
    std::shared_ptr<std::vector<btSoftBody::btVertexToTetraMapping>>& collisionMeshVertexToTetra,
    std::shared_ptr<std::vector<std::vector<int>>>& collisionMeshTriToTetra,
    const uint32_t maxVoxelCount,
    const btVector3& scale
)
{
    btTriangleIndexVertexArray* m_indexVertexArrays = new btTriangleIndexVertexArray(
        (int)indicesCollisionMesh->size() / 3, (int*)indicesCollisionMesh->data(), 3 * sizeof(unsigned int),
        (int)verticesCollisionMesh->size(), (btScalar*)verticesCollisionMesh->data(), sizeof(btVector3)
    );
    
    btBvhTriangleMeshShape* trimeshShape = new btBvhTriangleMeshShape(m_indexVertexArrays, true);
    btVector3 aabbMin, aabbMax;
    trimeshShape->getAabb(btTransform::getIdentity(), aabbMin, aabbMax);
    
    float voxelSize = findVoxelSizeForVoxelGridWithMaxVoxels(*trimeshShape, aabbMin, aabbMax, maxVoxelCount, scale);
    
    std::vector<int> vertexToTet(verticesCollisionMesh->size(), -1);
    collisionMeshTriToTetra = std::make_shared<std::vector<std::vector<int>>>(indicesCollisionMesh->size() / 3);
    
    uint32_t voxelsCreated = 0;
    VoxelGridCollector collector = { simulationVerticesForSoftBody, simulationTetrasForSoftBody, vertexToTet, *collisionMeshTriToTetra };
    
    createVoxelGrid(*trimeshShape, aabbMin, aabbMax, voxelSize, scale, voxelsCreated, &collector);
    
    collisionMeshVertexToTetra = std::make_shared<std::vector<btSoftBody::btVertexToTetraMapping>>(verticesCollisionMesh->size());
    for(size_t i=0; i<vertexToTet.size(); ++i) {
        if(vertexToTet[i] != -1) {
            (*collisionMeshVertexToTetra)[i].vertexToTetra = vertexToTet[i];
            const auto& tetra = simulationTetrasForSoftBody[vertexToTet[i]];
             (*collisionMeshVertexToTetra)[i].baryCoordInTetra = GetBarycentricCoordinatesInTetra(
                (*verticesCollisionMesh)[i], simulationVerticesForSoftBody[tetra[0]], simulationVerticesForSoftBody[tetra[1]],
                simulationVerticesForSoftBody[tetra[2]], simulationVerticesForSoftBody[tetra[3]]
             );
             
             // Simple normal end mapping (offset by normal)
             btVector3 n = (*normalsCollisionMesh)[i];
             btVector3 p_n = (*verticesCollisionMesh)[i] + n * 0.1f; // Short normal segment
             (*collisionMeshVertexToTetra)[i].baryCoordNormalEndInTetra = GetBarycentricCoordinatesInTetra(
                 p_n, simulationVerticesForSoftBody[tetra[0]], simulationVerticesForSoftBody[tetra[1]],
                 simulationVerticesForSoftBody[tetra[2]], simulationVerticesForSoftBody[tetra[3]]
             );
        }
    }
    
    delete trimeshShape;
    delete m_indexVertexArrays; 
}


// --- Custom GImpact Primitive Manager ---

class TrimeshDeformedPrimitiveManager : public btGImpactMeshShapePart::TrimeshPrimitiveManager
{
    std::shared_ptr<std::vector<btSoftBody::btVertexToTetraMapping>> collisionShapeVertexToSimTetra;
    std::shared_ptr<std::vector<std::vector<int>>> collisionShapeTriToSimTetra;
    btSoftBody* softBody;
    std::shared_ptr<std::vector<Vector3>> originalVertices;
    std::shared_ptr<std::vector<Index>> indices;

public:
    TrimeshDeformedPrimitiveManager(
        std::shared_ptr<std::vector<btSoftBody::btVertexToTetraMapping>> collisionShapeVertexToSimTetra,
        std::shared_ptr<std::vector<std::vector<int>>> collisionShapeTriToSimTetra,
        btSoftBody* sb,
        std::shared_ptr<std::vector<Vector3>> verts,
        std::shared_ptr<std::vector<Index>> inds
    ) : btGImpactMeshShapePart::TrimeshPrimitiveManager(),
        collisionShapeVertexToSimTetra(collisionShapeVertexToSimTetra),
        collisionShapeTriToSimTetra(collisionShapeTriToSimTetra),
        softBody(sb),
        originalVertices(verts),
        indices(inds)
    {
         // Setup Base Properties for non-virtual function access in base class
         m_margin = 0.01f;
         numverts = (int)verts->size();
         numfaces = (int)inds->size() / 3;
         
         #ifdef BT_USE_DOUBLE_PRECISION
         type = PHY_DOUBLE;
         #else
         type = PHY_FLOAT;
         #endif
         stride = sizeof(btVector3);
         vertexbase = (const unsigned char*)originalVertices->data();
         
         indicestype = PHY_INTEGER;
         indexstride = 3 * sizeof(unsigned int);
         indexbase = (const unsigned char*)indices->data();
    }
    
    virtual btPrimitiveManagerBase* clone() const override { return new TrimeshDeformedPrimitiveManager(*this); }
    
    virtual bool is_trimesh() const override { return true; }
    virtual int get_primitive_count() const override { return numfaces; }
    virtual int get_vertex_count() const { return numverts; }
    
    virtual void get_vertex(unsigned int vertex_index, btVector3& vertex, bool get_original) const override {
        if(get_original) {
            vertex = (*originalVertices)[vertex_index];
        } else {
             const auto& mapping = (*collisionShapeVertexToSimTetra)[vertex_index];
             int tIdx = mapping.vertexToTetra;
             if(tIdx >= 0 && softBody) {
                  const auto& n0 = softBody->m_tetras[tIdx].m_n[0];
                  const auto& n1 = softBody->m_tetras[tIdx].m_n[1];
                  const auto& n2 = softBody->m_tetras[tIdx].m_n[2];
                  const auto& n3 = softBody->m_tetras[tIdx].m_n[3];
                  
                  vertex = ApplyBarycentricCoordinatesInTetra(mapping.baryCoordInTetra, n0->m_x, n1->m_x, n2->m_x, n3->m_x);
             } else {
                  vertex = (*originalVertices)[vertex_index];
             }
        }
    }
    
    virtual void get_primitive_triangle(int prim_index, btPrimitiveTriangle& triangle, bool get_original) const override {
        unsigned int i0, i1, i2;
        // Call base get_indices (non-virtual but we set up the pointers for it to work)
        btGImpactMeshShapePart::TrimeshPrimitiveManager::get_indices(prim_index, i0, i1, i2);
        get_vertex(i0, triangle.m_vertices[0], get_original);
        get_vertex(i1, triangle.m_vertices[1], get_original);
        get_vertex(i2, triangle.m_vertices[2], get_original);
        triangle.m_margin = m_margin;
    }
    
    virtual bool get_primitive_triangle_safe(int prim_index, btPrimitiveTriangle& triangle) const override {
        return false;
    }
};

// --- Debug output ---

void framestart() { std::cout << "framestart()" << std::endl; }
void frameend() { std::cout << "frameend()" << std::endl; }
void drawline(const btVector3& from, const btVector3& to, const btVector3& color) {
    std::cout << "drawline \"line\" [" << from.x() << "," << from.y() << "," << from.z() << "][" 
              << to.x() << "," << to.y() << "," << to.z() << "][" 
              << color.x() << "," << color.y() << "," << color.z() << ",1]" << std::endl;
}


// --- Main ---

int main()
{
    std::cout << "Starting main..." << std::endl;
    // 1. Setup World
    btSoftBodyRigidBodyCollisionConfiguration* config = new btSoftBodyRigidBodyCollisionConfiguration();
    btCollisionDispatcher* dispatcher = new btCollisionDispatcher(config);
    btBroadphaseInterface* broadphase = new btDbvtBroadphase();
    btDeformableMultiBodyConstraintSolver* solver = new btDeformableMultiBodyConstraintSolver();
    solver->setDeformableSolver(new btDeformableBodySolver());
    
    btDeformableMultiBodyDynamicsWorld* world = new btDeformableMultiBodyDynamicsWorld(dispatcher, broadphase, solver, config);
    
    world->setGravity(btVector3(0, -9.8, 0));
    btGImpactCollisionAlgorithm::registerAlgorithm(dispatcher);
    std::cout << "World setup done." << std::endl;

    // 2. Rigid Body (Static Box)
    {
        btTriangleMesh* mesh = new btTriangleMesh();
        
        // Simple cube mesh
        btVector3 v0(-2,-2,-2), v1(2,-2,-2), v2(2,-2,2), v3(-2,-2,2);
        btVector3 v4(-2,2,-2), v5(2,2,-2), v6(2,2,2), v7(-2,2,2);
        
        // Bottom
        mesh->addTriangle(v0, v1, v2); mesh->addTriangle(v0, v2, v3);
        // Top
        mesh->addTriangle(v4, v6, v5); mesh->addTriangle(v4, v7, v6);
        // Front
        mesh->addTriangle(v3, v2, v6); mesh->addTriangle(v3, v6, v7);
        // Back
        mesh->addTriangle(v0, v5, v1); mesh->addTriangle(v0, v4, v5);
        // Left
        mesh->addTriangle(v0, v3, v7); mesh->addTriangle(v0, v7, v4);
        // Right
        mesh->addTriangle(v1, v5, v6); mesh->addTriangle(v1, v6, v2);
        
        btGImpactMeshShape* gimpactShape = new btGImpactMeshShape(mesh);
        gimpactShape->setLocalScaling(btVector3(1,1,1));
        gimpactShape->updateBound();
        
        btTransform startTransform;
        startTransform.setIdentity();
        startTransform.setOrigin(btVector3(0, -5, 0));
        
        btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
        btRigidBody::btRigidBodyConstructionInfo rbInfo(0, myMotionState, gimpactShape, btVector3(0,0,0)); // Mass 0 for static
        btRigidBody* body = new btRigidBody(rbInfo);
        
        world->addRigidBody(body);
        std::cout << "Rigid body added." << std::endl;
    }
    
    // 3. Soft Body
    btSoftBody* softBody = nullptr;
    {
        // 3a. Create detailed collision mesh (Subdivided Cube)
        // For simplicity, generate a regular grid of triangles
        auto vertices = std::make_shared<std::vector<btVector3>>();
        auto normals = std::make_shared<std::vector<btVector3>>();
        auto indices = std::make_shared<std::vector<unsigned int>>();
        
        int res = 4; // 4x4x4 subdivisions
        float size = 4.0f;
        float step = size / res;
        btVector3 offset(-size/2, -size/2, -size/2);
        
        // Helper to add quad
        auto addQuad = [&](btVector3 p0, btVector3 p1, btVector3 p2, btVector3 p3) {
            int base = (int)vertices->size();
            vertices->push_back(p0); vertices->push_back(p1); vertices->push_back(p2); vertices->push_back(p3);
            btVector3 n = (p1-p0).cross(p2-p0).normalized();
            normals->push_back(n); normals->push_back(n); normals->push_back(n); normals->push_back(n);
            indices->push_back(base); indices->push_back(base+1); indices->push_back(base+2);
            indices->push_back(base); indices->push_back(base+2); indices->push_back(base+3);
        };
        
        // Generate box shell
        for(int i=0; i<res; ++i) {
            for(int j=0; j<res; ++j) {
                float u = i*step; float v = j*step;
                float u1 = (i+1)*step; float v1 = (j+1)*step;
                // Top
                addQuad(offset+btVector3(u,size,v), offset+btVector3(u1,size,v), offset+btVector3(u1,size,v1), offset+btVector3(u,size,v1));
                // Bottom
                addQuad(offset+btVector3(u,0,v1), offset+btVector3(u1,0,v1), offset+btVector3(u1,0,v), offset+btVector3(u,0,v));
                // Front
                addQuad(offset+btVector3(u,v,size), offset+btVector3(u1,v,size), offset+btVector3(u1,v1,size), offset+btVector3(u,v1,size));
                 // Back
                addQuad(offset+btVector3(u1,v,0), offset+btVector3(u,v,0), offset+btVector3(u,v1,0), offset+btVector3(u1,v1,0));
                 // Left
                addQuad(offset+btVector3(0,v,u1), offset+btVector3(0,v,u), offset+btVector3(0,v1,u), offset+btVector3(0,v1,u1));
                // Right
                addQuad(offset+btVector3(size,v,u), offset+btVector3(size,v,u1), offset+btVector3(size,v1,u1), offset+btVector3(size,v1,u));
            }
        }
        std::cout << "Collision mesh generated." << std::endl;
        
        // 3b. Voxelize
        std::vector<btVector3> simVertices;
        std::vector<std::array<int, 4>> simTetras;
        auto vertexToTetra = std::make_shared<std::vector<btSoftBody::btVertexToTetraMapping>>();
        std::shared_ptr<std::vector<std::vector<int>>> triToTetra;
        
        std::cout << "Voxelizing..." << std::endl;
        createVoxelTetrahedronMesh(vertices, normals, indices, simVertices, simTetras, true, vertexToTetra, triToTetra, 100 /*maxVoxel*/, btVector3(1,1,1));
        std::cout << "Voxelization done. Verts: " << simVertices.size() << " Tetras: " << simTetras.size() << std::endl;
        
        // 3c. Create SoftBody
        softBody = new btSoftBody(&world->getWorldInfo(), (int)simVertices.size(), simVertices.data(), 0);
        
        for(auto& t : simTetras) {
            softBody->appendTetra(t[0], t[1], t[2], t[3]);
        }
        std::cout << "Soft body created." << std::endl;
        
        softBody->generateClusters(0); 
        softBody->m_cfg.collisions = btSoftBody::fCollision::CL_SS | btSoftBody::fCollision::CL_RS;
        
        softBody->m_materials[0]->m_kLST = 0.5;
        softBody->m_materials[0]->m_kAST = 0.5;
        softBody->m_materials[0]->m_kVST = 0.5;
        
        softBody->generateBendingConstraints(2);
        softBody->randomizeConstraints();
        
        TrimeshDeformedPrimitiveManager* primitiveManager = new TrimeshDeformedPrimitiveManager(vertexToTetra, triToTetra, softBody, vertices, indices);
        
        btTriangleIndexVertexArray* meshInterface = new btTriangleIndexVertexArray(
            (int)indices->size()/3, (int*)indices->data(), 3*sizeof(int),
            (int)vertices->size(), (btScalar*)vertices->data(), sizeof(btVector3)
        );
        
        btGImpactMeshShape* gimpactSoftShape = new btGImpactMeshShape(meshInterface, primitiveManager);
        gimpactSoftShape->setMargin(0.01f);
        gimpactSoftShape->updateBound();
        
        softBody->setCollisionShape(gimpactSoftShape);
        
        btTransform trans; trans.setIdentity(); trans.setOrigin(btVector3(0, 10, 0));
        softBody->transform(trans);
        
        world->addSoftBody(softBody);
        std::cout << "Soft body added to world." << std::endl;
    }
    
    // 4. Simulation Loop
    for (int i = 0; i < 200; i++)
    {
        framestart();
        
        world->stepSimulation(1.0f / 60.0f);
        
        // Print soft body nodes
        for (int n = 0; n < softBody->m_nodes.size(); ++n)
        {
            const auto& node = softBody->m_nodes[n];
            // std::cout << "v [" << node.m_x.x() << "," << node.m_x.y() << "," << node.m_x.z() << "]" << std::endl;
        }

        // Draw wireframe for soft body tetras (optional debug)
        const btVector3 color(0,1,0); 
        for(int t=0; t<softBody->m_tetras.size(); ++t) {
            const auto& tetra = softBody->m_tetras[t];
             drawline(tetra.m_n[0]->m_x, tetra.m_n[1]->m_x, color);
             drawline(tetra.m_n[0]->m_x, tetra.m_n[2]->m_x, color);
             drawline(tetra.m_n[0]->m_x, tetra.m_n[3]->m_x, color);
             drawline(tetra.m_n[1]->m_x, tetra.m_n[2]->m_x, color);
             drawline(tetra.m_n[1]->m_x, tetra.m_n[3]->m_x, color);
             drawline(tetra.m_n[2]->m_x, tetra.m_n[3]->m_x, color);
        }

        frameend();
    }

    // Cleanup (Simplified)
    delete world;
    delete solver;
    delete broadphase;
    delete dispatcher;
    delete config;
    
    return 0;
}
