#!/usr/bin/env python3
"""Create a distributable macOS DMG containing spinDecon.app.

The disk image presents spinDecon.app alongside an Applications symlink so
users can install by dragging the app into Applications.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import plistlib
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_APP = ROOT / "dist" / "spinDecon.app"
DEFAULT_DIST = ROOT / "dist"


def require_tool(name: str) -> str:
    path = shutil.which(name)
    if not path:
        raise SystemExit(f"Required macOS tool not found: {name}")
    return path


def app_version(app: Path) -> str:
    plist = app / "Contents" / "Info.plist"
    try:
        with plist.open("rb") as fh:
            info = plistlib.load(fh)
        return str(info.get("CFBundleShortVersionString") or info.get("CFBundleVersion") or "unknown")
    except (OSError, plistlib.InvalidFileException):
        return "unknown"


def main() -> int:
    parser = argparse.ArgumentParser(description="Create the spinDecon macOS installer DMG")
    parser.add_argument("--app", type=Path, default=DEFAULT_APP)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--volume-name", default="spinDecon Installer")
    args = parser.parse_args()

    if sys.platform != "darwin":
        raise SystemExit("DMG creation must be run on macOS")

    hdiutil = require_tool("hdiutil")
    codesign = require_tool("codesign")
    app = args.app.expanduser().resolve()
    if not app.is_dir():
        raise SystemExit(f"Application not found: {app}\nRun 'make app' first.")

    # Refuse to package a damaged/invalid app.
    subprocess.run([codesign, "--verify", "--deep", "--strict", "--verbose=2", str(app)], check=True)

    version = app_version(app)
    output = args.output or (DEFAULT_DIST / f"spinDecon-{version}-macOS.dmg")
    output = output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        output.unlink()

    with tempfile.TemporaryDirectory(prefix="spindecon-dmg-") as tmp:
        stage = Path(tmp) / args.volume_name
        stage.mkdir()

        # Preserve the app bundle exactly as signed. copytree(symlinks=True)
        # avoids dereferencing Python.framework aliases.
        shutil.copytree(app, stage / app.name, symlinks=True, copy_function=shutil.copy2)
        os.symlink("/Applications", stage / "Applications")

        subprocess.run([
            hdiutil, "create",
            "-volname", args.volume_name,
            "-srcfolder", str(stage),
            "-ov",
            "-format", "UDZO",
            "-imagekey", "zlib-level=9",
            str(output),
        ], check=True)

    subprocess.run([hdiutil, "verify", str(output)], check=True)
    print(f"Created installer: {output}")
    print("Install by opening the DMG and dragging spinDecon.app to Applications.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
