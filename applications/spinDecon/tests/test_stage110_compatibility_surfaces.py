from pathlib import Path

ROOT = Path(__file__).parents[1]


def test_frames_is_compatibility_surface_not_implementation_home():
    offenders = []
    for path in ROOT.joinpath('Frames').rglob('*.py'):
        if path.name == '__init__.py':
            continue
        # Compatibility wrappers are intentionally tiny.  A large file here
        # means implementation has leaked back into the historical package.
        if path.stat().st_size > 1024:
            offenders.append(str(path.relative_to(ROOT)))
    assert offenders == []
