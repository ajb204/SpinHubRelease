# Process label-state hydration fix

The Process session now hydrates persisted `label1`..`label4` values into
`ProjectState.gui_settings` before `ProcessFrame.GetLabs()` runs.  This makes
shared live state authoritative before Conversion or Projection is opened.

`ProcessFrame.get_spectral_labels()` is the Process-family semantic accessor.
Projection asks this accessor first and only uses its older disk/state lookup as
a compatibility fallback for isolated callers.

The session-start hydration does not mark the project dirty.  Conversion child
hydration still uses `seed_gui_settings`, so an already-edited live label is
never overwritten merely by opening/reopening Conversion.

Regression result: 55 focused GUI-state tests pass. Full suite: 254 passed,
1 skipped, 6 pre-existing unrelated failures.
