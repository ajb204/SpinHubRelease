from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_pseudo2d_always_has_service_and_usta_is_not_active_analysis():
    source = (ROOT / 'gui/workspaces/pseudo2d.py').read_text()
    assert 'self.pseudo_service = PseudoAxisService(tabOne)' in source
    assert "return ['Diffusion', 'Decay']" in source
    assert "return ['uSTA', 'Diffusion', 'Decay']" not in source


def test_pseudo2d_normal_path_does_not_branch_on_missing_service():
    source = (ROOT / 'gui/workspaces/pseudo2d.py').read_text()
    assert 'self.pseudo_service is not None' not in source
    assert 'self.pseudo_service is None' in source  # constructor compatibility bridge only
