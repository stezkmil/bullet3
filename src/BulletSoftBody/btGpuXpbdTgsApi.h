#ifndef BT_GPU_XPBD_TGS_API_H
#define BT_GPU_XPBD_TGS_API_H

// ABI 12 adds optional solved contact impulse reporting.
// POD boundary independent of Bullet's float/double scalar configuration.
constexpr int BT_GPU_XPBD_ABI_VERSION = 12;
struct btGpuXpbdVec { double x, y, z; };
struct btGpuXpbdNode
{
    btGpuXpbdVec position;
    double inverseMass;
    btGpuXpbdVec velocity;
    double damping;
};
struct btGpuXpbdTet
{
    int nodes[4];
    double inverseRest[9]; // column major
    double volume, young, poisson, damping;
};
struct btGpuXpbdStep
{
    double dt;
    int substeps;
    btGpuXpbdVec gravity;
    double settlingSpeed, settlingDamping;
    int iterations;
    int settlingBegin,settlingEnd; // end=0 selects all nodes.
    int velocityIterations;
    int materialMode; // 0: accurate direct StVK, 1: fast partitioned corotated.
};
struct btGpuXpbdSurfaceVertex
{
    int nodes[4]; // -1: unused; nodes[0] = -2-rigidIndex: rigid-local offset; all -1: static.
    double weights[4];
    btGpuXpbdVec offset;
    int body;
};
struct btGpuXpbdSurfaceTriangle { int vertices[3];int part=-1,index=-1; };
struct btGpuXpbdContactReport
{
    int bodyA,bodyB,partA,indexA,partB,indexB;
    btGpuXpbdVec pointA,pointB,normal;
    double distance;
};
// Normal impulse integrated over one outer step, aggregated by surface contact identity.
struct btGpuXpbdContactImpulse { btGpuXpbdContactReport contact; double impulse; };
static_assert(sizeof(btGpuXpbdContactImpulse)==112, "Contact impulse GPU ABI mismatch");
struct btGpuXpbdTarget { int node;btGpuXpbdVec position; };
struct btGpuXpbdDisabledPair { int bodyA, bodyB; };
struct btGpuXpbdNotificationPair { int bodyA,bodyB;int touching; };
static_assert(sizeof(btGpuXpbdSurfaceTriangle)==20 && sizeof(btGpuXpbdContactReport)==104 && sizeof(btGpuXpbdNotificationPair)==12, "Contact report GPU ABI mismatch");
struct btGpuXpbdSpring
{
    int node;
    btGpuXpbdVec target;
    double stiffness,damping,maxForce;
};
struct btGpuXpbdRigid
{
    btGpuXpbdVec position;
    double orientation[4]; // x,y,z,w; principal inertia frame.
    btGpuXpbdVec linearVelocity, angularVelocity;
    double inverseMass;
    btGpuXpbdVec inverseInertia;
    btGpuXpbdVec force, torque;
    double linearDamping, angularDamping;
};
struct btGpuXpbdEndpoint
{
    int rigid; // >=0: rigid index; -1: weighted soft nodes or fixed world point.
    int nodes[4];
    double weights[4];
    btGpuXpbdVec point; // rigid-local offset, or world offset for soft/fixed points.
};
struct btGpuXpbdAttachment
{
    btGpuXpbdEndpoint a,b;
    double compliance; // m/N in the backend length unit; zero gives a hard attachment.
};
struct btGpuXpbdJointAxis
{
    double lower,upper,stiffness,damping,equilibrium;
    double velocity,maxForce,servoTarget,bounce,stopERP,stopCFM,motorERP,motorCFM;
    int flags; // Spring: 1, motor: 2, servo: 4, stiffness cap: 8, damping cap: 16.
};
struct btGpuXpbdJoint
{
    int bodies[2]; // -1: fixed world frame.
    btGpuXpbdVec points[2];
    double orientations[2][4];
    btGpuXpbdJointAxis axes[6];
    int rotationOrder; // Bullet RotateOrder.
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
    int (*surface)(void*, int, const btGpuXpbdSurfaceVertex*, int, const btGpuXpbdSurfaceTriangle*, double, double);
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
    int (*friction)(void*, int, const double*);
    // ABI 5: compliant mouse attachments; count=0 releases all springs.
    // Target positions use the same coordinates as nodes; stiffness/damping are kg/s^2 and kg/s.
    int (*springs)(void*, int, const btGpuXpbdSpring*);
    // Mesh upload clears rigid bodies and attachments. Rigid upload clears surfaces and attachments.
    int (*rigids)(void*, int, const btGpuXpbdRigid*);
    int (*downloadRigids)(void*, int, btGpuXpbdRigid*);
    int (*attachments)(void*, int, const btGpuXpbdAttachment*);
    int (*updateRigids)(void*, int, const btGpuXpbdRigid*);
    int (*joints)(void*, int, const btGpuXpbdJoint*);
    int (*selfCollision)(void*, int, const int*); // Replace enabled surface body IDs; surface upload clears this set.
    // Query current geometry without advancing time. Notification pairs stop at their first hit.
    int (*contacts)(void*, int, const btGpuXpbdNotificationPair*, int*);
    int (*downloadContacts)(void*, int, btGpuXpbdContactReport*);
    int (*enableContactImpulses)(void*, int);
    // Null output queries the count; otherwise capacity must cover all reports. Cleared on each step.
    int (*contactImpulses)(void*, int*, btGpuXpbdContactImpulse*, int);
};
static_assert(sizeof(btGpuXpbdRigid)==200 && sizeof(btGpuXpbdEndpoint)==80 && sizeof(btGpuXpbdAttachment)==168, "Rigid GPU ABI mismatch");
using btGpuXpbdTgsGetApi = const btGpuXpbdTgsApi* (*)(int);
static_assert(sizeof(btGpuXpbdNode)==64 && sizeof(btGpuXpbdTet)==120 && sizeof(btGpuXpbdStep)==80, "GPU XPBD ABI mismatch");
#endif
