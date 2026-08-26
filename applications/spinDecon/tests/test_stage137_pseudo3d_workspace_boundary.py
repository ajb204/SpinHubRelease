from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

def test_pseudo3d_uses_canonical_workspace_for_gui_shell_access():
    source = (ROOT / 'gui/workspaces/pseudo3d.py').read_text()
    executable = [line for line in source.splitlines() if 'self.tabOne' in line and not line.lstrip().startswith('#')]
    assert executable == []
    assert 'self.app_context.nmr_workspace' in source
    assert "getattr(self.nmr_workspace, 'statusbar', None)" in source
