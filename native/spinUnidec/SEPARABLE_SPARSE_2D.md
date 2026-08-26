# Sparse separable 2D convolution benchmark

This build starts from the current best `ApplyIter()` optimisation (the 4.7 s baseline).

The only algorithmic change is the sparse 2D convolution in `src/slice2D.cpp`.

## What changed

The original sparse convolution evaluates every active output against every active source:

    O(n^2)

The experimental implementation exploits the separable kernel

    K(i,j) = pki[i] * pkj[j]

and performs two 1D convolutions.

Because `BuildSparseDB()` is ordered by `j`, sources are grouped into occupied source rows. The first pass convolves along `i` for each occupied source row and stores a compact intermediate of:

    unique_source_j * si

rather than the full `si * sj` spectrum.

The second pass evaluates only the currently active output points, combining the occupied intermediate rows along `j`.

For the profiling case previously measured, this is approximately 281 * 2048 = 575,488 intermediate doubles rather than 2048 * 1023 = 2,095,104.

## Benchmark output

Builds with `PROFILE_SEPARABLE_2D` enabled print:

    [PROFILE_SEPARABLE2D] sparse_n=... unique_j=... intermediate=... first_i_ms=... second_j_ms=... total_ms=...

The normal computation remains unchanged outside the sparse convolution and the existing active-index `ApplyIter()` optimisation.

## A/B switch

The source defaults to the separable path:

    #define SEPARABLE_SPARSE_2D 1

To benchmark the old O(n^2) sparse convolution from the same source tree, compile with:

    -DSEPARABLE_SPARSE_2D=0

No cutoff or other approximation is introduced.

## Numerical note

The separable formulation is algebraically equivalent to the original separable kernel, but floating-point summation order changes. Compare the final spectra and convergence metrics against the 4.7 s baseline before treating it as numerically interchangeable.
