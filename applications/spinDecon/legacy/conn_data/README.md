# conn_data quarantine

Stage 148D quarantines residual GUI functionality whose only purpose was the
historical `conn_data` connectivity library.  The Full Peak List is the
authoritative peak collection in the modern 1D/pseudo2D/2D/pseudo3D/3D
journeys.

Snapshots preserve pre-removal implementations for future NOE/connectivity
recovery.  Any future recovery should use an explicit relationship model keyed
to Full Peak List identities rather than restoring `conn_data` as a competing
peak store.

`conn_data` still exists inside the NMR controller/DataStore in Stage 148D
because deconvolution result loading for current 2D/3D journeys still passes
through that representation.  That remaining use is migration debt and was
not removed without a dedicated replacement.
