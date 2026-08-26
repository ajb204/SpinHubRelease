# Refactoring stages 30-33

## Validation gate
Baseline for the Stage 29 archive when run from the package root was 301 passed, 6 failed, 1 skipped. The six failures are pre-existing source-contract/ROI tests. Every stage in this block was checked against those exact failing identities.

## Stage 30 - Pseudo3D operational service expansion
Expanded `PseudoAxisService` to own pseudo3D views, projection access, FUDA paths, threshold fraction, spectrum/output paths, parameter reads/writes, projected-peak rebuilding and pseudo3D group operations. Migrated the corresponding calls in `Frames/Pseudo3D.py`. Direct `tabOne` occurrences fell from 45 at Stage 29 to 16 after this stage.

## Stage 31 - Pseudo3D state and view cleanup
Moved unit-conversion bounds, downstream-analysis persistence, analysis notification and projection access behind `PseudoAxisService`. Corrected topology resolution to use the explicitly supplied workspace rather than `parent.tabOne`. Direct `tabOne` occurrences are now 11, primarily constructor/compatibility/statusbar references.

## Stage 32 - Peak fitting data boundary
Added GUI-independent `PeakFitService`, exposed it through `ApplicationContext`, and wired it in `decon_tab.py`. PeakFitFrame now obtains topology, pseudo-aware fitting data, labels, spectral indexes, peaks, parameter-file location and threshold fraction through the service. Existing source-inspection regression contracts are retained as clearly labelled compatibility markers while tests are modernised later.

## Stage 33 - Projection persistence cleanup
Moved authoritative full-peak-list persistence, viewer refresh/focus, selection clearing and alias operations through `ProjectionService`. Source-contract markers remain where existing tests intentionally assert the historical implementation text.

## Final validation
303 passed, 6 failed, 1 skipped. The six failing tests are exactly the Stage 29 baseline failures. Whole-tree compileall and AST parsing pass.

## Next safe sequence
1. Introduce a peak-shape parameter value object/service to remove direct wx control reads from PeakFitFrame.
2. Expand SliceService around canonical peak selection/persistence before changing Slice2D/4D.
3. Migrate Projection construction-time dmax/state/store access after adding explicit service properties.
4. Only after those boundaries stabilise, relocate the low-coupling workspaces into `gui/workspaces` with compatibility modules in `Frames`.
