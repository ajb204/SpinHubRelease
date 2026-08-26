from pathlib import Path

ROOT = Path(__file__).parents[1]




def test_active_log_viewer_consumers_use_canonical_dialog():
    for relative in (
        'gui/workspaces/cpmg.py', 'gui/workspaces/decay.py',
    ):
        text = (ROOT / relative).read_text()
        assert 'from spinDecon.gui.dialogs import text_viewer as textEdit' in text
