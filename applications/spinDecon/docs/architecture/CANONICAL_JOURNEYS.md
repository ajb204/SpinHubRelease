# Canonical Analysis Journeys

## Governing model

The current application supports five canonical journeys: **1D, pseudo2D, 2D,
pseudo3D and 3D**.  Dataset topology decides the journey.  The journey then
decides which spectra/views are useful and which peak-list roles are
scientifically independent.

All authoritative peak collections use the common peak-list representation.
Dimensionality is data carried by that representation, not a reason to create a
new peak-list type.

Three peak concepts must be kept distinct:

1. **Full Peak List** — authoritative peak identities for the main spectral
   dataset.
2. **Reference Peak List** — a scientifically independent projection used to
   locate/constrain work on a higher-dimensional dataset.  It exists only where
   the journey actually needs such a projection.
3. **Projection peak list** — derived marker coordinates used to show known peak
   locations on projected/protected views.  It is display data derived from an
   authoritative list, not a second source of peak identity.

The **Deconvolved Peak List** is a result role using the same canonical peak
schema.  It never uses `conn_data`.

## 1D

A physical 1D spectrum has one spectral axis and one authoritative Full Peak
List.  There is no independent reference spectrum/list.  Projected marker data,
if a view requires it, is derived from Full.

Journey: prepare spectrum -> determine peak shape -> pick Full peaks -> review ->
optional fit/deconvolution -> Deconvolved Peak List.

## pseudo2D

Pseudo2D has one spectral axis plus one pseudo axis.  The Full 1D Peak List is
both the authoritative peak collection and the frequency reference used to
extract/follow each peak through the pseudo-axis series.  It must not be copied
into an independent reference authority merely because a reference-style view
is useful.

Journey: prepare series -> determine peak shape -> establish/review Full 1D peaks
-> extract intensity series across pseudo axis -> review series -> downstream
experiment analysis where applicable.  Deconvolution results use the common
peak schema.

## physical 2D

Physical 2D is the important singularity.  The **Full, reference and projection
peak concepts refer to the same physical 2D peak collection**.  Full is the
source of authority.  A projected/protected view may cache transformed marker
coordinates, but it must not create a second authoritative reference list.

Journey: prepare 2D spectrum -> determine peak shape -> pick Full 2D peaks ->
review/check peaks -> restrained fit/deconvolution -> review fitting results ->
Deconvolved Peak List.

## pseudo3D

Pseudo3D has two spectral axes plus one pseudo axis.  Here an independent 2D
reference/projection is scientifically meaningful: it identifies the spectral
locations to follow through the pseudo dimension.  The main/full collection and
that reference projection may therefore coexist, while projected marker lists
remain derived display data.

Journey: prepare series -> determine peak shape -> establish 2D reference peaks
-> extract/follow peaks through pseudo axis -> review intensity series ->
downstream analysis (for example fitting/relaxation/diffusion as appropriate).

## physical 3D

Physical 3D has a genuinely independent lower-dimensional 2D reference
projection.  The reference list constrains/organises work in the 3D spectrum;
the Full 3D Peak List remains authoritative for peaks in the main dataset.
Protected 2D views display derived projections of these canonical lists.

Journey: prepare 3D spectrum -> determine peak shape -> establish 2D reference
projection/list -> pick Full 3D peaks -> inspect/review in 2D slices ->
Deconvolved Peak List/results as produced by the analysis.

## Spectrum/view decision rule

A spectrum or view belongs in the active architecture only when it serves one of
these journeys: the main spectrum, a scientifically required reference
projection, a pseudo-axis series, a deconvolved result spectrum, or a derived
view needed to inspect those objects.  A view does not gain authority over peak
identity merely because it displays peaks.

This rule is also the quarantine test.  uSTA, INDIANA, MAGMA, Slice4D,
NOE/`conn_data` and historical compatibility surfaces remain under `legacy/`
until deliberately recovered into this model.
