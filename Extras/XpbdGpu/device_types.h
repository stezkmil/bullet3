#pragma once

// Internal FP32 storage; the public Bullet API remains ABI 13 with FP64 PODs.
namespace gpu
{
struct btGpuXpbdVec { float x, y, z; };
struct btGpuXpbdNode
{
    btGpuXpbdVec position;
    float inverseMass;
    btGpuXpbdVec velocity;
    float damping;
    // FP32 accuracy: retain sub-ULP position and velocity updates.
    btGpuXpbdVec positionResidual,velocityResidual;
};
struct btGpuXpbdTet
{
    int nodes[4];
    float inverseRest[9];
    float volume, young, poisson, damping;
    // FP32 accuracy: rest edges support cancellation-resistant displacement gradients.
    btGpuXpbdVec restEdges[3],restEdgeResiduals[3];
};
struct btGpuXpbdStep
{
    float dt;
    int substeps;
    btGpuXpbdVec gravity;
    float settlingSpeed, settlingDamping;
    int iterations;
    int settlingBegin,settlingEnd;
    int velocityIterations;
    int materialMode;
};
struct btGpuXpbdSurfaceVertex
{
    int nodes[4];
    float weights[4];
    btGpuXpbdVec offset;
    int body;
};
struct btGpuXpbdSurfaceTriangle { int vertices[3];int part=-1,index=-1; };
struct btGpuXpbdContactReport
{
    int bodyA,bodyB,partA,indexA,partB,indexB;
    btGpuXpbdVec pointA,pointB,normal;
    float distance;
};

struct btGpuXpbdContactImpulse { btGpuXpbdContactReport contact; float impulse; };

struct btGpuXpbdTarget { int node;btGpuXpbdVec position; };
struct btGpuXpbdDisabledPair { int bodyA, bodyB; };

struct btGpuXpbdNotificationPair { int bodyA,bodyB;int touching; };

struct btGpuXpbdSpring
{
    int node;
    btGpuXpbdVec target;
    float stiffness,damping,maxForce;
};
struct btGpuXpbdRigid
{
    btGpuXpbdVec position;
    float orientation[4];
    btGpuXpbdVec linearVelocity, angularVelocity;
    float inverseMass;
    btGpuXpbdVec inverseInertia;
    btGpuXpbdVec force, torque;
    float linearDamping, angularDamping;
};
struct btGpuXpbdEndpoint
{
    int rigid;
    int nodes[4];
    float weights[4];
    btGpuXpbdVec point;
};
struct btGpuXpbdAttachment
{
    btGpuXpbdEndpoint a,b;
    float compliance;
};
struct btGpuXpbdJointAxis
{
    float lower,upper,stiffness,damping,equilibrium;
    float velocity,maxForce,servoTarget,bounce,stopERP,stopCFM,motorERP,motorCFM;
    int flags;
};
struct btGpuXpbdJoint
{
    int bodies[2];
    btGpuXpbdVec points[2];
    float orientations[2][4];
    btGpuXpbdJointAxis axes[6];
    int rotationOrder;
};
struct btGpuXpbdTgsApi
{
    void* (*create)();
    void (*destroy)(void*);
    int (*upload)(void*, int, const btGpuXpbdNode*, int, const btGpuXpbdTet*);
    int (*forces)(void*, int, const btGpuXpbdVec*);

    int (*step)(void*, const btGpuXpbdStep*);
    int (*download)(void*, int, btGpuXpbdNode*);
    const char* (*error)(void*);
    int (*surface)(void*, int, const btGpuXpbdSurfaceVertex*, int, const btGpuXpbdSurfaceTriangle*, float, float);
    int (*surfacePositions)(void*, int, btGpuXpbdVec*);
    int (*embed)(void*, int, const btGpuXpbdVec*, int, int, int, btGpuXpbdSurfaceVertex*);
    int (*targets)(void*, int, const btGpuXpbdTarget*);

    int (*disabledPairs)(void*, int, const btGpuXpbdDisabledPair*);

    int (*updateNodes)(void*, int, const btGpuXpbdNode*);

    int (*friction)(void*, int, const float*);

    int (*springs)(void*, int, const btGpuXpbdSpring*);

    int (*rigids)(void*, int, const btGpuXpbdRigid*);
    int (*downloadRigids)(void*, int, btGpuXpbdRigid*);
    int (*attachments)(void*, int, const btGpuXpbdAttachment*);
    int (*updateRigids)(void*, int, const btGpuXpbdRigid*);
    int (*joints)(void*, int, const btGpuXpbdJoint*);
    int (*selfCollision)(void*, int, const int*);

    int (*contacts)(void*, int, const btGpuXpbdNotificationPair*, int*);
    int (*downloadContacts)(void*, int, btGpuXpbdContactReport*);
    int (*enableContactImpulses)(void*, int);

    int (*contactImpulses)(void*, int*, btGpuXpbdContactImpulse*, int);
};

}
const gpu::btGpuXpbdTgsApi* gpuApi();

static_assert(sizeof(gpu::btGpuXpbdNode)==56 && sizeof(gpu::btGpuXpbdTet)==140 && sizeof(gpu::btGpuXpbdStep)==48,"FP32 simulation layout");
static_assert(sizeof(gpu::btGpuXpbdSurfaceVertex)==48 && sizeof(gpu::btGpuXpbdRigid)==100 && sizeof(gpu::btGpuXpbdEndpoint)==48,"FP32 coupling layout");
static_assert(sizeof(gpu::btGpuXpbdContactReport)==64 && sizeof(gpu::btGpuXpbdContactImpulse)==68,"FP32 report layout");
