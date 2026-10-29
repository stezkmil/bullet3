# GPU XPBD/TGS accuracy baseline

Elasticity uses the St. Venant–Kirchhoff energy constraint from
PositionBasedDynamics `XPBD::solve_FEMTetraConstraint` (MIT), evaluated without
inversion repair. Young’s modulus and Poisson’s ratio determine the Lamé
coefficients. Rest volume is included in both energy and gradient.

Every advancing increment integrates force and applies one material pass.
Tetrahedra are greedily graph-colored; colors execute in order and each kernel
updates the shared physical nodes directly. No copied-node mass splitting or
averaging is used. A fresh XPBD multiplier is appropriate because each pass
advances time. This is still a local XPBD/TGS solver, not a global FEM solve.

For normalized energy U (Young’s modulus factored out), the constraint is
C = sqrt(2U). The implementation eliminates this square root algebraically:

    dx_i = w_i grad(U)_i * (-2U) / (sum(w_j |grad(U)_j|^2) + 2U/(E h^2))

Stiffness-proportional damping is a separate implicit velocity projection over
six orthonormal symmetric strain modes (bulk, two diagonal deviatoric and three
shear modes). It runs before force integration, so it vanishes at static rest
and does not dilute the elastic projection. Each local mode dissipates kinetic
energy and conserves momentum. Mass-proportional damping is exponential and
also precedes force integration. This local operator splitting does not certify
transient damping accuracy at arbitrary increments; those require convergence
tests too.
No absolute, dimension-dependent gradient cutoff is used.

ABI 7 stores positions, velocities, material data and contact state in double
precision. Older runtimes are rejected; rebuild the adapter and deploy its
matching DLL together. The borrowed PhysX triangle GJK/EPA witness calculation
still uses its original float precision.

The increment h = outer dt / (substeps * iterations) and mesh resolution must
be convergence-tested. Finite increments and linear tetrahedral discretization
can change measured stiffness. Inversion robustness and automatic timestep
selection are outside this change.

Tests:
- XpbdCoreTest / XpbdLoaderTest: analytic material projection over three unit
  scales and Poisson ratios, force grouping, momentum, contacts and friction.
- Test_btGpuXpbdWorld: adapter units, mappings, mouse springs, collision filters
  and notification-only obstacles.
- VRUT `gpuMaterialAccuracy.js`: existing uniformly loaded cantilever geometry,
  selectable Young’s modulus, substeps, iterations, damping and tetra mesher.
  Construction runs at zero load with notification detectors registered first;
  the gravity load starts only after the fixture and output geometry exist.
