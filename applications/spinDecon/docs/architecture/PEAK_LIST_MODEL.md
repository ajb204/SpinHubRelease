# Canonical Peak-List Model

## Architectural rule

All current peak lists use one common, dimension-independent representation.
Dimensionality changes the number and meaning of coordinates; it does not create
a different peak-list type.

A dataset or workflow may legitimately retain several peak lists at the same
time. Their **role**, not their storage type, distinguishes them. Typical roles
are:

- `reference` / projection peak list - peaks used to define or inspect a lower-dimensional projection/reference;
- `full` - the authoritative complete peak list for the main dataset;
- `decon` - the calculated/deconvolved peak list;
- derived projected lists - views derived from an authoritative list for display.

Canonical nD records contain a peak identity/name, a coordinate tuple,
axis-labelled values where available, optional intensity, source-row metadata,
and an extensible `analysis` mapping. The same schema is used for 1D, pseudo2D,
2D, pseudo3D and 3D journeys and is intended to extend to future dimensionalities.

## Authority and ownership

The Full Peak List is authoritative for the complete main-dataset peak
collection. A reference/projection list may coexist where the scientific journey
requires it. A deconvolved list is a separate result collection and must not be
encoded as connectivity data or overwrite the Full Peak List merely because it
shares coordinates with it.

`DataStore.peak_lists` is the current shared storage boundary for these named
roles. `projected_peak_lists` contains derived display projections and is not a
second source of peak authority.

## Legacy connectivity

`conn_data` and `connEntry` are deprecated connectivity/NOE structures and are
not part of the current peak model. Their historical implementation is retained
under `legacy/` for later recovery of useful NOE/connectivity behaviour.

If connectivity/NOE functionality is restored, it must use a dedicated
connection model that references canonical peak identities. It must not restore
`conn_data` as a peak-list representation or source of peak authority.

## Journey-specific identity rules

The detailed journey policy is specified in `CANONICAL_JOURNEYS.md` and encoded
by `domain.spectrum_policy`.

- **1D:** Full is authoritative; no independent reference list.
- **pseudo2D:** Full 1D is authoritative and also supplies reference frequencies.
- **physical 2D:** Full, reference and projection concepts identify the same 2D
  peak collection. Full remains the sole authority; projected marker payloads
  are derived caches only.
- **pseudo3D:** an independent 2D reference projection may coexist with the main
  list because it defines spectral locations followed through the pseudo axis.
- **physical 3D:** an independent 2D reference projection coexists with the
  authoritative Full 3D list.

A **projection peak list** used on a projected/protected view is not synonymous
with an independent **reference peak list**. Projection payloads exist to render
canonical peaks in a different view and must be rebuildable from their source
peak list.
