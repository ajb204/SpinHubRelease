ApplyIter generalization benchmark
===================================

Starting from the current 2D best architecture.

Changes:
- slice1D::ApplyIter uses sparseDB when SPARSE is true, dense fallback otherwise.
- slice3D::ApplyIter uses sparseDB when SPARSE is true, dense fallback otherwise.
- 3D bore mode explicitly rebuilds each synchronized 1D sparseDB before update.
- 4D mode explicitly rebuilds each 2D sparseDB after the 4D DS is copied back.
- No convolution code, cutoff, kernel, or convergence criterion is changed.

Coverage:
1D; 3D mode 1; 3D mode 3 with bore; 3D full; 4D mode 4; 4D mode 2.
