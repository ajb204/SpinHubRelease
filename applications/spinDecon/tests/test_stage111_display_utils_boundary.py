from pathlib import Path

ROOT = Path(__file__).parents[1]


def test_active_gui_uses_canonical_display_helpers():
    active = [
        ROOT / 'gui' / 'dialogs' / 'processing' / 'projections.py',
        ROOT / 'gui' / 'workspaces' / 'projection.py',
        ROOT / 'gui' / 'workspaces' / 'peak_review.py',
    ]
    for path in active:
        source = path.read_text()
        assert 'decon.misc.display_utils' not in source
