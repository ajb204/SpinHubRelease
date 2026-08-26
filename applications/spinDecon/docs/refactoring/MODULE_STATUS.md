# Module status ledger

| Module | Status | Rationale / action |
|---|---|---|
| Frames/Pseudo2Dold.py | LEGACY-HIGH-CONFIDENCE | Explicit old implementation; quarantine after dependency check. |
| catia_tab.py | INVESTIGATE-PROBABLE-LEGACY | Parallel application shell. |
| Frames/catiaApp.py | INVESTIGATE-PROBABLE-LEGACY | CATIA-specific shell/support. |
| Frames/STD_frame.py | INVESTIGATE-PROBABLE-DUPLICATE | Compare with Frames/uSTA/STD_frame.py. |
| Frames/NOEframe.py | INVESTIGATE | Confirm callbacks/dynamic use. |
| LoadFilePopUp.py | INVESTIGATE | Confirm dynamic use. |
| PeakShapeOptimizer.py | INVESTIGATE-ACTIVE | Scientific utility; determine current callers. |
| notepad.py | INVESTIGATE | Utility/standalone candidate. |
| multi_tab.py | INVESTIGATE | Possible historical shell. |
| shiftXPostFilter.py | INVESTIGATE | Scientific utility; determine current callers. |
| Frames/SetDataStoreFrame.py | INVESTIGATE | Possible transitional state GUI. |
| Frames/uindecNMRFrame.py | ACTIVE-OPTIONAL | UniDec integration candidate. |
| Frames/frameFeatures.py | INVESTIGATE-INFRASTRUCTURE | Generic GUI helper but may be hierarchy-coupled. |

## Stage 17 legacy quarantine

- `Frames/Pseudo2Dold.py` — **LEGACY (high confidence)**. The implementation
  is now quarantined at `legacy/pseudo2d_old.py`; the old path is a compatibility
  import only. No active internal importer was found.
- `Frames/STD_frame.py` — **PROBABLE LEGACY/DUPLICATE**. The active application
  imports `Frames/uSTA/STD_frame.py`; retained pending external-use confirmation.
- `Frames/NOEframe.py` — **LEGACY-CANDIDATE**. Referenced only by the old
  `multi_tab.py` shell in static internal analysis; retained pending quarantine
  of that shell.
- `catia_tab.py`, `Frames/catiaApp.py`, `multi_tab.py` — **LEGACY-CANDIDATE
  APPLICATION SHELLS**. Kept in place because their standalone entry-point use
  cannot be disproved by the regression suite.
