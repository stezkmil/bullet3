# XPBD/TGS cantilever validation — 2026-10-20

Production Bullet GPU backend, exercised through VRUT’s existing uniformly loaded
cantilever geometry and both tetrahedral mesh generators. Double-precision
material/state arithmetic; outer timestep 0.0002 s; 50 advancing iterations per
substep. Mass damping alpha=60 s^-1 accelerates settling.

| E (GPa) | Mesh | Substeps | Material beta (s) | Tip deflection (mm) | Test target (mm) | Result |
|---:|---|---:|---:|---:|---:|---|
| 15 | voxel | 1 | 0 | 13.203191 | 11.500000 | pass |
| 15 | voxel | 2 | 0 | 11.028232 | 11.500000 | pass |
| 15 | voxel (unloaded) | 16 | 0 | 0.000000 | 0.000000 | pass |
| 15 | conforming | 2 | 0 | 10.294547 | 11.500000 | pass |
| 70 | voxel | 4 | 0 | 2.357316 | 2.464286 | pass |
| 15 | voxel | 2 | 0.01 | 10.896727 | 11.500000 | pass |

The loaded cases pass the existing 15% deflection tolerance. The unloaded case
has zero measured drift (0.001 mm tolerance), with no rejected solver step.
The final test script requires two successive stable samples and at least 750
loaded updates. The initial 1- and 2-substep undamped voxel runs used the older
single-stable-sample stopping rule; their values should not be read as tighter
equilibrium estimates than that rule supports.

The coarse 1-substep result is substantially softer than the 2-substep result.
This is finite-increment error, not a change in material constants. These checks
establish static stiffness within the beam test tolerance; they do not establish
arbitrary-step accuracy or converged transient damping. Mesh discretization
error remains visible. Inversion handling was deliberately not changed.

All three native suites passed: XpbdCoreTest, XpbdLoaderTest and
Test_btGpuXpbdWorld. Their coverage includes analytic StVK projection over three
coordinate scales and Poisson ratios; unchanged static stiffness under damping;
constant-force increment grouping; momentum; friction; mouse springs; mapped
contacts; and notification-only collision objects.

The test setup registers the notification detectors before the beam and keeps
gravity zero until the fixture exists. This avoids accidentally measuring a
body already distorted by setup-time contacts.

Reproduce from the VRUT bin directory with separate quoted command arguments:

```powershell
$testArgs = @(
    "-noautoexec", "-G0", "-RJSScripting", "-Dkernel.logToConsole=1",
    "-DJSScripting.scriptPath=C:/Devel/vrut/base_modules/JSScripting/scripts/jsScripting/gpuMaterialAccuracy.js",
    "-DJSScripting.runScript=1", "-DTest.BeamYoung=15000000", "-DTest.BeamSubsteps=2"
)
& ./vrut.exe @testArgs
```

Set Test.BeamBeta=0.01 for damping, Test.BeamConforming=1 for conforming
tetrahedra, Test.BeamYoung=70000000 for 70 GPa, or Test.BeamUnloaded=1 for
the unloaded control. Reject any run whose console log reports a solver step
rejection, even if its geometry stops moving.
