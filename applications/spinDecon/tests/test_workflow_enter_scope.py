"""Regression guard for Workflow's Enter-key scope."""
from pathlib import Path


def test_workflow_enter_is_local_to_visible_page():
    source = (Path(__file__).parents[1] / "gui" / "workspaces" / "workflow.py").read_text()
    assert "button.SetDefault()" not in source
    assert "self.Bind(wx.EVT_CHAR_HOOK, self._on_char_hook)" in source
    assert "self.notebook.GetPage(selected) is self" in source
    assert "self._enter_action_key = stage.key" in source
