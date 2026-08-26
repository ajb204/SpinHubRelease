"""Compatibility entry point for the quarantined legacy CATIA application."""
from spinDecon.legacy.catia.catia_tab import NotebookDemo, MyApp

__all__ = ["NotebookDemo", "MyApp"]

if __name__ == "__main__":
    app = MyApp()
    app.MainLoop()
