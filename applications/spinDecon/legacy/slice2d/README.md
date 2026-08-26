# Legacy Slice2D connectivity/NOE code

This directory preserves the pre-Stage-148C Slice2D implementation containing
`conn_data`, NOE, and MAGMA behaviours.  It is quarantine/reference code only
and must not be imported by the active application.

The modern Slice2D journey uses the authoritative Full Peak List.  If NOE or
connectivity functionality is restored later, recover the useful algorithms
behind a dedicated connection/NOE model that references Full Peak List peak
identities; do not restore `conn_data` as canonical peak state.
