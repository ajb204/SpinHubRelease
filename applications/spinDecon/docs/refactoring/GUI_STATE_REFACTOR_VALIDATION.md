# GUI State / Save / Script Refactor - Final Validation

## Final state contract

1. Open GUI controls are the most current representation of user edits.
2. `ProjectState.gui_settings` is the authoritative live cross-window state.
3. The parameter/system file is a persistence boundary, not an inter-window message bus.
4. Process **Save** is the single parameter commit action.
5. Before processing or conversion execution, current GUI state is centrally committed.
6. Processing scripts receive an immutable `ProcessingScriptState` captured from current widgets with live-state fallback.
7. All 1D, 2D, pseudo-3D, 3D and 4D NMRPipe builders consume that explicit state.
8. Projection phase slider values remain transient preview state until **Re-process** promotes them into processing state.
9. Script-editor Save actions remain separate because they save script text, not parameter state.

## Adversarial validation

`tests/test_gui_state_end_to_end_journeys.py` deliberately gives disk, shared state and widgets conflicting values. It verifies that the widget value wins, is persistable, and is the value supplied to every dimensional processing route. It also verifies projection promotion, close/reopen live-state overlay, immutable processing snapshots and conversion snapshot ordering.

Focused state/script regression suite: 66 passed.
Full test suite at Step 16: 243 passed, 1 skipped, 6 pre-existing unrelated failures.
