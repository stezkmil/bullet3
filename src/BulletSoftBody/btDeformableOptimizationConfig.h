#ifndef BT_DEFORMABLE_OPTIMIZATION_CONFIG_H
#define BT_DEFORMABLE_OPTIMIZATION_CONFIG_H
// Optimized kernels: 0 selects original arithmetic; 63 enables all paths.
// Override the mask when comparing subsets in full-scene stability tests.
// Bits: 1 scaleAndAdd, 2 elastic, 4 translation input, 8 coupling,
//       16 mass, 32 multAndAddTo.
#ifndef BT_DEFORMABLE_OPTIMIZATION_MASK
#define BT_DEFORMABLE_OPTIMIZATION_MASK 63
#endif
#endif
