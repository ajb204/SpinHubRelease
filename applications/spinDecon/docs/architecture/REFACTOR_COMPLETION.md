# Refactor Completion Boundary

## Supported core

The active project is organised around the five canonical journeys documented
in `CANONICAL_JOURNEYS.md`: 1D, pseudo2D, physical 2D, pseudo3D and physical 3D.
Current source may support these journeys and their shared infrastructure.

## Quarantine boundary

Historical functionality remains under `legacy/` for deliberate future
recovery. Active source must not import it. In particular uSTA, INDIANA, MAGMA,
Slice4D, NOE/connectivity and `conn_data` are not current peak/workflow models.

## Peak authority

All authoritative peak lists use the common peak representation. Full is the
main-dataset authority; Deconvolved is a result role. Independent Reference
lists exist only where a higher-dimensional journey needs a scientific
projection. Projection marker payloads are rebuildable display caches.
Physical 2D has one peak authority: Full also fulfils the reference/projection
identity role.

## Definition of structurally complete

The structural refactor is complete when:

- active code imports no `legacy`, `Frames`, `misc`, or retired root façade;
- `conn_data`/`connEntry` have no active data-model role;
- workflow routing exposes only supported current workspaces;
- the canonical spectrum policy determines reference-list ownership;
- all current peak-producing paths publish the common peak representation;
- legacy code remains quarantined and recoverable;
- persistence compatibility is handled at load/save boundaries rather than by
  restoring retired runtime models;
- architecture and journey tests remain green.

Further changes should be treated as ordinary behavioural verification,
scientific feature work, or targeted simplification rather than continuation of
the package-layout refactor.
