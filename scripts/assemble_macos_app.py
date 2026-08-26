#!/usr/bin/env python3
"""Assemble a self-contained macOS spinDecon.app from a completed spinHub build.

The build venv supplies third-party packages.  Its base CPython framework is
copied into Contents/Frameworks so the resulting app does not depend on the
Homebrew Python installation used to build it.  External NMR software such as
NMRPipe/MDDNMR is deliberately not bundled.
"""
from __future__ import annotations

import argparse
import json
import os
import plistlib
import shutil
import stat
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_APP = ROOT / "dist" / "spinDecon.app"
BUNDLE_ID = "uk.ac.ox.spinDecon"
EXCLUDE_NAMES = {".DS_Store", ".git", ".pytest_cache", "__pycache__", ".build"}


def ignore(_directory: str, names: list[str]) -> set[str]:
    ignored = {n for n in names if n in EXCLUDE_NAMES or n.startswith("._")}
    ignored.update(n for n in names if n.endswith((".pyc", ".pyo")))
    return ignored


def copytree(src: Path, dst: Path, *, symlinks: bool = True) -> None:
    shutil.copytree(src, dst, symlinks=symlinks, ignore=ignore)


def require(path: Path, description: str) -> None:
    if not path.exists():
        raise SystemExit(f"Cannot build app: missing {description}: {path}")


def run_text(args: list[str]) -> str:
    return subprocess.check_output(args, text=True).strip()


def python_build_info(venv_python: Path) -> dict[str, str]:
    code = r'''import json, pathlib, sys, sysconfig
print(json.dumps({
 "base_prefix": sys.base_prefix,
 "executable": sys.executable,
 "version": f"{sys.version_info.major}.{sys.version_info.minor}",
 "site_packages": sysconfig.get_paths()["purelib"],
}))'''
    return json.loads(run_text([str(venv_python), "-c", code]))


def find_python_framework(base_prefix: Path) -> Path:
    candidates = [base_prefix, *base_prefix.parents]
    for p in candidates:
        if p.name == "Python.framework":
            return p
        q = p / "Python.framework"
        if q.is_dir():
            return q
    raise SystemExit(
        "Cannot build self-contained app: the build Python is not a macOS "
        f"framework build (base prefix: {base_prefix})."
    )


def copy_python_runtime(venv: Path, frameworks: Path, payload: Path) -> tuple[Path, str]:
    """Copy CPython framework and venv site-packages into the app.

    We intentionally do not copy .venv/bin: console-script shebangs and Python
    symlinks there point at the build machine.  spinDecon invokes the bundled
    interpreter directly and gets packages from the private runtime directory.
    """
    build_python = venv / "bin" / "python"
    info = python_build_info(build_python)
    version = info["version"]
    source_framework = find_python_framework(Path(info["base_prefix"]).resolve())
    dest_framework = frameworks / "Python.framework"
    print(f"Bundling Python {version} framework from {source_framework}")
    copytree(source_framework, dest_framework, symlinks=True)

    # Homebrew's framework contains a relative site-packages symlink whose
    # target lives outside Python.framework.  Once the framework is relocated
    # into the app that link is dangling.  Third-party packages are bundled
    # separately below under Resources/spinHub/python-runtime, so the framework
    # link is neither needed nor valid in the distributable app.
    framework_site = (
        dest_framework / "Versions" / version / "lib"
        / f"python{version}" / "site-packages"
    )
    if framework_site.is_symlink():
        print(f"Removing non-relocatable Python site-packages symlink: {framework_site}")
        framework_site.unlink()

    source_site = Path(info["site_packages"])
    runtime = payload / "python-runtime"
    dest_site = runtime / "lib" / f"python{version}" / "site-packages"
    dest_site.parent.mkdir(parents=True, exist_ok=True)
    copytree(source_site, dest_site, symlinks=True)

    bundled_python = dest_framework / "Versions" / version / "bin" / f"python{version}"
    require(bundled_python, "bundled Python executable")
    return bundled_python, version


def rewrite_python_framework_links(framework: Path) -> None:
    """Make the copied Python.framework independent of its build-machine path.

    Homebrew's framework executables may reference Python via an absolute path
    or via @rpath while carrying an absolute Homebrew LC_RPATH.  For every
    Mach-O inside the copied framework, point Python-framework dependencies
    directly at the bundled framework library using @loader_path.  This avoids
    relying on LC_RPATH altogether.  Absolute Homebrew/MacPorts-style rpaths
    are removed as well.
    """
    version_dirs = [p for p in (framework / "Versions").iterdir()
                    if p.is_dir() and not p.is_symlink() and p.name != "Current"]
    python_libs = {p.name: p / "Python" for p in version_dirs if (p / "Python").is_file()}

    for path in framework.rglob("*"):
        if path.is_symlink() or not path.is_file():
            continue
        try:
            deps = run_text(["otool", "-L", str(path)])
        except (subprocess.CalledProcessError, OSError):
            continue

        for line in deps.splitlines()[1:]:
            dep = line.strip().split(" (", 1)[0]
            marker = "Python.framework/Versions/"
            if marker not in dep:
                continue
            suffix = dep.split(marker, 1)[1]
            parts = suffix.split("/", 1)
            if len(parts) != 2 or parts[1] != "Python":
                continue
            version = parts[0]
            target = python_libs.get(version)
            if target is None:
                continue
            rel = os.path.relpath(target, start=path.parent)
            new = "@loader_path/" + rel
            if dep != new:
                subprocess.run(
                    ["install_name_tool", "-change", dep, new, str(path)],
                    check=True,
                )

        # Remove absolute package-manager rpaths inherited from the build Python.
        # They are neither necessary nor desirable in a relocatable .app.
        try:
            load_commands = run_text(["otool", "-l", str(path)])
        except (subprocess.CalledProcessError, OSError):
            continue
        lines = load_commands.splitlines()
        rpaths: list[str] = []
        for i, line in enumerate(lines):
            if line.strip() == "cmd LC_RPATH":
                for following in lines[i + 1:i + 5]:
                    stripped = following.strip()
                    if stripped.startswith("path "):
                        rpaths.append(stripped[5:].split(" (offset", 1)[0])
                        break
        for rpath in rpaths:
            if rpath.startswith(("/opt/homebrew/", "/usr/local/", "/opt/local/")):
                subprocess.run(
                    ["install_name_tool", "-delete_rpath", rpath, str(path)],
                    check=True,
                )

    # Give the framework library a relocatable install id. Consumers above use
    # @loader_path, but a relocatable id also keeps tooling output portable.
    for version, lib in python_libs.items():
        subprocess.run(
            ["install_name_tool", "-id",
             f"@rpath/Python.framework/Versions/{version}/Python", str(lib)],
            check=True,
        )


def reject_external_symlinks(app: Path) -> None:
    """Reject symlinks that would make the staged app non-portable."""
    bad: list[str] = []
    for p in app.rglob("*"):
        if not p.is_symlink():
            continue
        target = os.readlink(p)
        if os.path.isabs(target):
            bad.append(f"{p.relative_to(app)} -> {target} [absolute]")
            continue
        if not p.exists():
            bad.append(f"{p.relative_to(app)} -> {target} [dangling]")
            continue
        try:
            p.resolve(strict=True).relative_to(app)
        except (ValueError, FileNotFoundError):
            bad.append(f"{p.relative_to(app)} -> {target} [external]")
    if bad:
        raise SystemExit(
            "Cannot sign app: non-portable symlink(s) remain:\n  "
            + "\n  ".join(bad[:30])
        )


def write_launcher(path: Path, version: str) -> None:
    text = f'''#!/bin/sh
set -u
HERE="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
CONTENTS="$(CDPATH= cd -- "$HERE/.." && pwd)"
RESOURCES="$CONTENTS/Resources"
ROOT="$RESOURCES/spinHub"
PYTHON="$CONTENTS/Frameworks/Python.framework/Versions/{version}/bin/python{version}"
SITE="$ROOT/python-runtime/lib/python{version}/site-packages"
APP="$ROOT/applications/spinDecon/spinDecon"

export PATH="$ROOT/bin:/usr/bin:/bin:/usr/sbin:/sbin:${{PATH:-}}"
export PYTHONPATH="$SITE${{PYTHONPATH:+:$PYTHONPATH}}"
export DYLD_FRAMEWORK_PATH="$CONTENTS/Frameworks${{DYLD_FRAMEWORK_PATH:+:$DYLD_FRAMEWORK_PATH}}"
export SPINDECON_APP_BUNDLE="$CONTENTS"
export SPINDECON_ROOT="$ROOT"
export SPINDECON_BUNDLED_PYTHON=1

if [ ! -x "$PYTHON" ]; then
    osascript -e 'display alert "spinDecon cannot start" message "The bundled Python runtime is missing or invalid. Please reinstall spinDecon." as critical' >/dev/null 2>&1 || true
    exit 1
fi
exec "$PYTHON" "$APP" "$@"
'''
    path.write_text(text, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def compile_staged_python(build_python: Path, root: Path, version: str) -> None:
    """Precompile staged sources with the build interpreter.

    Do not execute the copied Python.framework here. install_name_tool has just
    modified its Mach-O files, invalidating Homebrew's signature; on macOS the
    kernel may SIGKILL that executable until the finished app is re-signed.
    The build venv uses the same Python major/minor ABI, so its checked-hash
    bytecode is valid for the staged runtime.
    """
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "python-runtime" / "lib" / f"python{version}" / "site-packages")
    subprocess.run(
        [str(build_python), "-m", "compileall", "-q", "-f", "-j", "0",
         "--invalidation-mode", "checked-hash",
         str(root / "applications")], check=True, env=env,
    )


def smoke_test(python: Path, root: Path, version: str) -> None:
    """Exercise the signed embedded runtime without modifying the app bundle.

    The application has already been code-signed when this test runs.  Python
    normally writes ``__pycache__/*.pyc`` files while importing modules; doing
    that inside a signed bundle changes its sealed resources and makes the
    subsequent DMG verification fail.  ``-B`` plus PYTHONDONTWRITEBYTECODE make
    this a read-only relocation/import test.
    """
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "python-runtime" / "lib" / f"python{version}" / "site-packages")
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    subprocess.run(
        [str(python), "-B", "-c", "import wx,numpy,scipy,matplotlib,nmrglue,fuda,fudaIO,fudalib,dataIO; print('Bundled Python runtime: OK')"],
        check=True, env=env,
    )



def sign_macos_app(app: Path) -> None:
    """Ad-hoc sign the staged bundle without traversing framework symlink aliases.

    ``codesign --deep`` is deliberately avoided here. Python.framework contains
    the normal ``Versions/Current`` symlink tree; recursive deep signing can
    encounter those aliases as separate bundle objects. We sign real Mach-O
    files first, then nested bundles/frameworks, and the outer app last.
    """
    contents = app / "Contents"
    framework = contents / "Frameworks" / "Python.framework"

    def is_macho(path: Path) -> bool:
        if path.is_symlink() or not path.is_file():
            return False
        try:
            out = subprocess.check_output(["file", "-b", str(path)], stderr=subprocess.DEVNULL)
        except (OSError, subprocess.CalledProcessError):
            return False
        return b"Mach-O" in out

    # os.walk(..., followlinks=False) is important: never descend through
    # Python.framework/Versions/Current or the framework's top-level aliases.
    macho_files: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(contents, followlinks=False):
        base = Path(dirpath)
        dirnames[:] = [d for d in dirnames if not (base / d).is_symlink()]
        for name in filenames:
            candidate = base / name
            if is_macho(candidate):
                macho_files.append(candidate)

    # Deepest files first. Re-signing is harmless and ensures every extension
    # module/dylib modified or copied into the bundle has a valid signature.
    for path in sorted(macho_files, key=lambda x: len(x.parts), reverse=True):
        subprocess.run(["codesign", "--force", "--sign", "-", str(path)], check=True)

    nested_python_app = framework / "Versions"
    real_versions = [p for p in nested_python_app.iterdir() if p.is_dir() and not p.is_symlink()] if nested_python_app.is_dir() else []
    for version_dir in real_versions:
        python_app = version_dir / "Resources" / "Python.app"
        if python_app.is_dir():
            subprocess.run(["codesign", "--force", "--sign", "-", str(python_app)], check=True)

    if framework.is_dir():
        subprocess.run(["codesign", "--force", "--sign", "-", str(framework)], check=True)

    subprocess.run(["codesign", "--force", "--sign", "-", str(app)], check=True)
    subprocess.run(["codesign", "--verify", "--deep", "--strict", "--verbose=2", str(app)], check=True)

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_APP)
    parser.add_argument("--version", default=os.environ.get("SPINDECON_VERSION", "3.0"))
    args = parser.parse_args()
    if sys.platform != "darwin":
        raise SystemExit("spinDecon.app assembly must be run on macOS")

    venv = ROOT / ".venv"
    require(venv, "project virtual environment (run make setup)")
    require(ROOT / "applications" / "spinDecon" / "spinDecon", "spinDecon entry point")
    require(ROOT / "bin", "built bin directory")

    app = args.output.expanduser().resolve()
    if app.exists(): shutil.rmtree(app)
    contents = app / "Contents"
    macos, resources, frameworks = contents / "MacOS", contents / "Resources", contents / "Frameworks"
    payload = resources / "spinHub"
    macos.mkdir(parents=True); payload.mkdir(parents=True); frameworks.mkdir(parents=True)

    copytree(ROOT / "applications", payload / "applications")
    copytree(ROOT / "bin", payload / "bin")
    if (ROOT / "requirements.txt").exists(): shutil.copy2(ROOT / "requirements.txt", payload / "requirements.txt")

    bundled_python, pyver = copy_python_runtime(venv, frameworks, payload)
    rewrite_python_framework_links(frameworks / "Python.framework")
    write_launcher(macos / "spinDecon", pyver)

    # Install the application icon before signing so it is included in the
    # bundle resource seal.  The checked-in ICNS is generated from the
    # canonical PNG in assets/macos/.
    icon_source = ROOT / "assets" / "macos" / "spinDecon.icns"
    require(icon_source, "spinDecon macOS application icon")
    shutil.copy2(icon_source, resources / "spinDecon.icns")

    plist = {
        "CFBundleDevelopmentRegion": "English", "CFBundleDisplayName": "spinDecon",
        "CFBundleExecutable": "spinDecon", "CFBundleIdentifier": BUNDLE_ID,
        "CFBundleIconFile": "spinDecon.icns",
        "CFBundleInfoDictionaryVersion": "6.0", "CFBundleName": "spinDecon",
        "CFBundlePackageType": "APPL", "CFBundleShortVersionString": args.version,
        "CFBundleVersion": args.version, "LSMinimumSystemVersion": "11.0",
        "NSHighResolutionCapable": True, "NSPrincipalClass": "NSApplication",
    }
    with (contents / "Info.plist").open("wb") as fh: plistlib.dump(plist, fh, sort_keys=True)

    # Precompile with the still-valid build interpreter. The copied framework
    # has had install names rewritten and must not be executed until signed.
    compile_staged_python(venv / "bin" / "python", payload, pyver)
    reject_external_symlinks(app)

    sign_macos_app(app)

    # Only execute the embedded interpreter after codesign has repaired the
    # signatures invalidated by install_name_tool. This is also the relocation
    # smoke test for the packaged Python and extension modules.
    smoke_test(bundled_python, payload, pyver)
    print(f"Created {app}")
    print(f"Bundled Python: {pyver} ({bundled_python})")
    print(f"Open with: open {shlex_quote(str(app))}")
    return 0


def shlex_quote(value: str) -> str:
    import shlex
    return shlex.quote(value)


if __name__ == "__main__": raise SystemExit(main())
