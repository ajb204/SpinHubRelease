# Optimized ApplyIter benchmark

This build starts from the original `spinUnidec_v3(1).zip` source.

Only `slice2D::ApplyIter()` was changed.

## Change

`decon::ApplyIter()` calls `calcspec()` immediately before `slice2D::ApplyIter()`.
For the sparse 2D route, `slice2D::CalcSpec()` rebuilds `sparseDB` from the current non-zero `DB` entries.

The original `ApplyIter()` then scanned the entire `si * sj` dense array and tested each value for zero.

The benchmark version instead loops directly over `sparseDB`, using the already-known linear index `ii` for each active source:

    DB[ii] = DB[ii] * fabs(DI[ii] / DS[ii]);
    tack += fabs(DB[ii]);

No change was made to the numerical update formula, `CalcSpec()`, the sparse convolution, or the construction of `sparseDB`.

## Expected effect

For `si * sj` much larger than the number of active sparse sources, this removes the dense `si * sj` traversal from `ApplyIter()`.

The sparse list is rebuilt by `CalcSpec()` immediately before the update, so no additional support/cache logic is introduced.

## Validation

`src/slice2D.cpp` passes a syntax-only C++ compilation check when the repository's existing missing standard-library compatibility includes are supplied on the compiler command line. Those unrelated include issues are unchanged from the original source.
