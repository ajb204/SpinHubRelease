#!/usr/bin/env python3
"""Create the project venv and install Python dependencies."""
from pathlib import Path
import os
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
VENV = ROOT / ".venv"

def run(*args):
    print("+", *map(str, args))
    subprocess.run([str(a) for a in args], check=True)

def bootstrap_python():
    requested = os.environ.get("PYTHON_BOOTSTRAP")
    if requested:
        exe = Path(shutil.which(requested) or requested).resolve()
        if not exe.exists():
            raise SystemExit(f"PYTHON_BOOTSTRAP not found: {requested}")
        return exe

    # Do not accidentally seed this project from an activated venv belonging
    # to another checkout. Search PATH again with that venv removed.
    active = os.environ.get("VIRTUAL_ENV")
    paths = os.environ.get("PATH", "").split(os.pathsep)
    if active:
        active_bin = str((Path(active) / "bin").resolve())
        paths = [x for x in paths if str(Path(x).resolve()) != active_bin]
    exe = shutil.which("python3", path=os.pathsep.join(paths))
    if not exe:
        raise SystemExit("No system Python 3 found. Set PYTHON_BOOTSTRAP=/path/to/python3")
    exe = Path(exe).resolve()
    if VENV in exe.parents:
        raise SystemExit("Refusing to bootstrap .venv from itself")
    return exe

def main():
    if not VENV.exists():
        seed = bootstrap_python()
        print(f"Bootstrap Python: {seed}")
        run(seed, "-m", "venv", VENV)
    py = VENV / "bin" / "python"
    run(py, "-m", "pip", "install", "--upgrade", "pip", "setuptools", "wheel")
    run(py, "-m", "pip", "install", "-r", ROOT / "requirements.txt")

if __name__ == "__main__":
    main()
