from pathlib import Path

ROOT = Path(__file__).parents[1]


def test_notebook_uses_canonical_nmr_workspace_internally():
    source = (ROOT / "app" / "notebook.py").read_text()
    assert "self.nmr_workspace = NMRWorkspace(" in source
    assert "attach_analysis_services(self.app_context, self.nmr_workspace)" in source
    assert "self.tabOne = self.nmr_workspace" in source
    # tabOne is retained only as the backwards-compatible public alias in the notebook shell.
    executable = [line.strip() for line in source.splitlines() if "self.tabOne" in line and not line.lstrip().startswith("#")]
    assert executable == ["self.tabOne = self.nmr_workspace"]
