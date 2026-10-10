Experimental CUDA mapped-surface operations

Scope: detailed mapped-surface displacement limits inside VBD color updates,
and dense mapped-position evaluation for collision snapshots (32,768+ vertices).
CPU contact-plane guards, tetrahedral volume guards, elastic solves, and GImpact
collision discovery/validation remain active. This is not a complete GPU solver.

Build the standalone CMake project in this directory with an installed CUDA toolkit
and Python 3 on Windows. For example from the Bullet root:
  cmake -S Extras/VbdGpu -B ../build-vbd-gpu -G "Visual Studio 18 2026" -A x64
  cmake --build ../build-vbd-gpu --config Release
Copy Release/BulletVbdGpu.dll beside vrut.exe (and beside the Bullet test executable
when running GPU tests). NVRTC is used at build time only; the resulting DLL embeds
PTX and depends on the installed NVIDIA driver. Target: compute capability 7.5+.

Enable "GPU surface operations" under Experimental VBD, or set the VRUT parameter
CollisionDetection.DeformableVbdGpuGuards=1 after loading the scene. Default is off.
This prototype requires double-precision Bullet. Unavailable/failed backends revert
to CPU evaluation; VBD_STEP records gpu_guards/gpu_fallbacks and VBD_GPU_FALLBACK
records the reason. gpu_surfaces/gpu_surface_fallbacks count collision snapshots. A failed backend remains disabled until the world is recreated.

The opt-in GPU tests are DISABLED_GpuMappedMotionGuardMatchesCpuAcrossMappingChanges,
DISABLED_GpuMappedGuardRetainsContactAndVolumeProtection, and
DISABLED_GpuMappedGuardCoversSharedAndDuplicateSupports,
DISABLED_GpuScalarMappingsRejectNonFinitePositions, and
DISABLED_GpuSurfaceEvaluationMatchesCpuAndBenchmarksDownloads. Run them explicitly
with --gtest_also_run_disabled_tests and a matching --gtest_filter. The separate
DISABLED_GpuUnavailableFallsBackToCpu test requires the DLL to be absent.

Mappings/device allocations persist in the world's collision cache. Node positions
are uploaded for each detailed guard call and one reduced scalar is downloaded.
Reference mapped positions are refreshed when contact discovery changes its reference.
Dense collision snapshots download all mapped positions for CPU GImpact processing.
Transfer and synchronization costs must be included in performance comparisons.

Mappings consisting entirely of scalar-weight identity Jacobians use compact
node/weight records and specialized kernels. Signed and zero weights are supported.
General Jacobians retain the full-matrix path; classification is renewed whenever
the exact mapping changes. Both paths use double precision with FMA disabled.
