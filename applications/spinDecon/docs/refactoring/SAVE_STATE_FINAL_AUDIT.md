# Save-state refactor final audit

## Final persistence contract

- Process-family widgets publish current values to `ProjectState.gui_settings`.
- File > Save / Ctrl-S performs the normal project save and, when focus belongs to the Process family, calls `ProcessFrame.save_current_gui_state()`.
- `ProcessFrame.save_current_gui_state()` is the sole Process-family parameter-file commit path.
- Processing and conversion execution commit current GUI state before script generation.
- Processing scripts consume an immutable `ProcessingScriptState` captured from current widgets/live state.
- Projection phase slider values remain transient until Re-process promotes them into processing state.
- Closing Process compares persistable values with the accepted baseline; Yes saves, No restores the live baseline without writing disk.

## Cleanup performed

- Removed obsolete `ProcessFrame.fileNotSaved` state.
- Removed the obsolete `ProcessFrame.OnButtonSave` compatibility handler.
- Removed obsolete parameter-save handlers from Processing, Conversion and ProcessProjection.
- Removed now-unused `update_parameter_file` imports from those child frames.
- Retained script-editor Save handlers: they save script text and are a separate operation.
- Retained `ProcessingFrame.reload_from_file()` and the optional loader support in `_ensure_processing_frame()` as a conservative compatibility API; active auto-processing never requests a stale reload.
- Updated regression tests so they assert the obsolete parameter-save handlers are absent rather than requiring compatibility delegates.

## Audit results

Within the Process-family frames, direct `update_parameter_file()` persistence is confined to `ProcessFrame.save_current_gui_state()`. Other application subsystems retain their own independent persistence paths where appropriate.

Focused save/state/script regression tests: 54 passed.

Complete suite: 251 passed, 1 skipped, 6 failed. The six failures are the same unrelated pre-existing failures present before this cleanup and are outside the Process save-state refactor.
