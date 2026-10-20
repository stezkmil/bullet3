#ifndef BT_GPU_XPBD_TGS_API_H
#define BT_GPU_XPBD_TGS_API_H

// POD boundary independent of Bullet's float/double scalar configuration.
struct btGpuXpbdVec { float x, y, z; };
struct btGpuXpbdNode
{
    btGpuXpbdVec position;
    float inverseMass;
    btGpuXpbdVec velocity;
    float damping;
};
struct btGpuXpbdTet
{
    int nodes[4];
    float inverseRest[9]; // column major
    float volume, young, poisson, damping;
};
struct btGpuXpbdStep
{
    float dt;
    int substeps;
    btGpuXpbdVec gravity;
    float settlingSpeed, settlingDamping;
    int iterations;
    int settlingBegin,settlingEnd; // end=0 selects all nodes.
    int velocityIterations;
};
struct btGpuXpbdSurfaceVertex
{
    int nodes[4]; // -1 denotes an unused support; all -1 denotes a static vertex.
    float weights[4];
    btGpuXpbdVec offset;
    int body;
};
struct btGpuXpbdSurfaceTriangle { int vertices[3]; };
struct btGpuXpbdTarget { int node;btGpuXpbdVec position; };
struct btGpuXpbdDisabledPair { int bodyA, bodyB; };
struct btGpuXpbdSpring
{
    int node;
    btGpuXpbdVec target;
    float stiffness,damping,maxForce;
};
struct btGpuXpbdTgsApi
{
    void* (*create)();
    void (*destroy)(void*);
    int (*upload)(void*, int, const btGpuXpbdNode*, int, const btGpuXpbdTet*);
    int (*forces)(void*, int, const btGpuXpbdVec*);
    // Numerical/contact-capacity rejection preserves the accepted simulation state.
    // A CUDA execution/context error requires recreating the backend.
    int (*step)(void*, const btGpuXpbdStep*);
    int (*download)(void*, int, btGpuXpbdNode*);
    const char* (*error)(void*);
    int (*surface)(void*, int, const btGpuXpbdSurfaceVertex*, int, const btGpuXpbdSurfaceTriangle*, float, float);
    int (*surfacePositions)(void*, int, btGpuXpbdVec*);
    int (*embed)(void*, int, const btGpuXpbdVec*, int, int, int, btGpuXpbdSurfaceVertex*);
    int (*targets)(void*, int, const btGpuXpbdTarget*);
    // ABI 2: replace the disabled body-pair set; count=0 enables all distinct bodies.
    // Uploading a mesh or collision surface clears this set.
    int (*disabledPairs)(void*, int, const btGpuXpbdDisabledPair*);
    // ABI 3: replace node state without changing rest geometry, topology or filters.
    // Release kinematic targets before replacing their node state.
    int (*updateNodes)(void*, int, const btGpuXpbdNode*);
    // ABI 4: friction indexed by surface body ID; missing entries are frictionless.
    // Surface/mesh uploads clear the table. Pair coefficients use Bullet's product rule.
    int (*friction)(void*, int, const float*);
    // ABI 5: compliant mouse attachments; count=0 releases all springs.
    // Target positions use the same coordinates as nodes; stiffness/damping are kg/s^2 and kg/s.
    int (*springs)(void*, int, const btGpuXpbdSpring*);
};
using btGpuXpbdTgsGetApi = const btGpuXpbdTgsApi* (*)(int);
static_assert(sizeof(btGpuXpbdNode)==32 && sizeof(btGpuXpbdTet)==68 && sizeof(btGpuXpbdStep)==44, "GPU XPBD ABI mismatch");
#endif
