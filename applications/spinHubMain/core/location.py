"""Explicit browser location state, independent of the process working directory."""
from pathlib import Path


class BrowserLocation:
    """Track the local area displayed by SpinHub and its navigation history."""

    def __init__(self, root=None):
        self._root = self._normalise(root if root is not None else Path.cwd())
        self._history = []

    @staticmethod
    def _normalise(path):
        path = Path(path).expanduser().resolve(strict=False)
        if not path.is_dir():
            raise ValueError(f'Browser location is not a directory: {path}')
        return path

    @property
    def root(self):
        return self._root

    @property
    def can_back(self):
        return bool(self._history)

    @property
    def can_up(self):
        return self._root.parent != self._root

    def go(self, path):
        target = self._normalise(path)
        if target != self._root:
            self._history.append(self._root)
            self._root = target
        return self._root

    def up(self):
        if not self.can_up:
            return self._root
        return self.go(self._root.parent)

    def back(self):
        if not self._history:
            return self._root
        self._root = self._history.pop()
        return self._root
