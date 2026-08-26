# Stages 77-81

## Stage 77 - Remove canonical deconFrame helper imports

Moved peak-list parsing to the GUI-independent `misc/peak_io.py` helper and changed
canonical Projection/Slice/Peak Review code to import peak records, connection records,
parameter parsing, and local numeric helpers from their actual owners rather than
`Frames/deconFrame.py`. Fixed stale relative imports created by earlier physical moves in
the processing, uSTA, and UniDec canonical modules.

## Stage 78 - Complete processing GUI physical migration

Migrated source-contract tests to the canonical `gui/dialogs/processing/` modules, then
replaced the five historical processing modules in `Frames/` with compatibility imports.
The tests still inspect the full implementations; they now inspect their authoritative
locations rather than duplicated historical sources.

## Stage 79 - NMRWorkspace becomes canonical

Physically migrated `Frames/deconFrame.py` to `gui/workspaces/nmr.py` and renamed the
implementation class to `NMRWorkspace`. `deconFrame` remains as a compatibility alias and
the old Frames module is now a two-line wrapper. `decon_tab.py` and UniDec use the
canonical class. Source-contract tests were migrated to the canonical implementation.

## Stage 80 - Downstream pseudo-axis workspaces

Physically migrated CPMG and decay-analysis implementations to `gui/workspaces/cpmg.py`
and `gui/workspaces/decay.py`. Pseudo3D and project-summary code now use those canonical
locations; historical Frames modules are wrappers.

## Stage 81 - Legacy quarantine and architecture guards

Quarantined the disconnected `SetDataStoreFrame`, old shared `frameFeatures`, and the old
NOE workspace under `legacy/`, retaining compatibility wrappers. NOE/connectivity code is
preserved because useful functionality may be restored later, but it is explicitly not a
peak authority. Added architecture tests preventing canonical GUI/integration modules from
reintroducing `deconFrame` imports and protecting the physical migrations.

## Authority rule

The Full Peak List remains authoritative. Legacy `conn_data` and preserved NOE code are
compatibility/history only. Any restored connectivity functionality should use a dedicated
connection model/service referencing canonical peak identities.
