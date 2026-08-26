#!/usr/bin/env python3
from __future__ import annotations
import shutil, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / 'bin'
LIBDIR = BIN / 'lib'
SYSTEM_PREFIXES = ('/System/Library/', '/usr/lib/')
FORBIDDEN_PREFIXES = ('/opt/homebrew/', '/usr/local/', '/opt/local/', '/sw/')


def run(*args: str, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(args, text=True, capture_output=True, check=check)


def deps(path: Path) -> list[str]:
    out = run('otool', '-L', str(path)).stdout.splitlines()[1:]
    return [line.strip().split(' (compatibility version', 1)[0]
            for line in out if line.strip()]


def rpaths(path: Path) -> list[str]:
    lines = run('otool', '-l', str(path)).stdout.splitlines()
    result: list[str] = []
    for i, line in enumerate(lines):
        if line.strip() == 'cmd LC_RPATH':
            for later in lines[i + 1:i + 5]:
                s = later.strip()
                if s.startswith('path '):
                    result.append(s[5:].split(' (offset', 1)[0])
                    break
    return result


def is_system(dep: str) -> bool:
    return dep.startswith(SYSTEM_PREFIXES) or dep.startswith('/System/Library/')


def expand_token(path: str, owner: Path) -> Path:
    path = path.replace('@loader_path', str(owner.parent))
    path = path.replace('@executable_path', str(BIN))
    return Path(path).expanduser()


def resolve_dep(dep: str, owner: Path) -> Path | None:
    if dep.startswith('/'):
        p = Path(dep)
        return p if p.exists() else None
    if dep.startswith('@loader_path/'):
        p = expand_token(dep, owner).resolve()
        return p if p.exists() else None
    if dep.startswith('@executable_path/'):
        p = expand_token(dep, owner).resolve()
        return p if p.exists() else None
    if dep.startswith('@rpath/'):
        suffix = dep[len('@rpath/'):]
        # First honour the Mach-O LC_RPATH entries.  Homebrew GCC uses this for
        # libgfortran -> libquadmath/libgcc_s, so these dependencies must be
        # followed before the package can be considered portable.
        for rp in rpaths(owner):
            if rp.startswith('@rpath'):
                continue
            candidate = expand_token(rp, owner) / suffix
            if candidate.exists():
                return candidate.resolve()
        # A previously bundled dependency may already be alongside its owner.
        for base in (owner.parent, LIBDIR):
            candidate = base / suffix
            if candidate.exists():
                return candidate.resolve()
    return None


def change(owner: Path, old: str, new: str) -> None:
    subprocess.run(['install_name_tool', '-change', old, new, str(owner)], check=True)


def set_id(lib: Path) -> None:
    subprocess.run(['install_name_tool', '-id', f'@loader_path/{lib.name}', str(lib)], check=True)


def sign(path: Path) -> None:
    subprocess.run(['codesign', '--force', '--sign', '-', str(path)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def bundle_for(exe: Path) -> None:
    # Queue entries carry both the file being modified and, where applicable,
    # its original source path.  LC_RPATH resolution is done against the source
    # copy because install_name_tool/copying can change the useful context.
    queue: list[tuple[Path, Path]] = [(exe, exe)]
    seen: set[Path] = set()
    while queue:
        owner, resolve_owner = queue.pop(0)
        owner = owner.resolve()
        if owner in seen:
            continue
        seen.add(owner)
        owner_is_bundled_lib = owner.parent == LIBDIR.resolve()
        for dep in deps(owner):
            if is_system(dep):
                continue
            src = resolve_dep(dep, resolve_owner)
            if src is None:
                # A dependency already rewritten into bin/lib is valid only if
                # the referenced file really exists.
                rewritten = resolve_dep(dep, owner)
                if rewritten is not None:
                    continue
                raise RuntimeError(f'Cannot resolve dependency {dep!r} required by {owner}')
            if is_system(str(src)):
                continue
            LIBDIR.mkdir(parents=True, exist_ok=True)
            dst = LIBDIR / src.name
            if not dst.exists():
                shutil.copy2(src, dst)
                subprocess.run(['chmod', 'u+w', str(dst)], check=True)
                set_id(dst)
            new = f'@loader_path/{dst.name}' if owner_is_bundled_lib else f'@loader_path/lib/{dst.name}'
            if dep != new:
                change(owner, dep, new)
            queue.append((dst, src))


def macho_files() -> list[Path]:
    result = []
    for p in BIN.iterdir():
        if not p.is_file() or p.is_symlink():
            continue
        if 'Mach-O' in run('file', str(p), check=False).stdout:
            result.append(p)
    return result


def validate(paths: list[Path]) -> None:
    bad: list[tuple[Path, str]] = []
    all_paths = paths + ([p for p in LIBDIR.iterdir() if p.is_file()] if LIBDIR.exists() else [])
    for p in all_paths:
        for dep in deps(p):
            if is_system(dep):
                continue
            if dep.startswith(FORBIDDEN_PREFIXES):
                bad.append((p, dep))
                continue
            if dep.startswith('@loader_path/') or dep.startswith('@executable_path/'):
                if resolve_dep(dep, p) is None:
                    bad.append((p, f'{dep} (missing bundled file)'))
            elif dep.startswith('@rpath/'):
                # Release artifacts should not retain unresolved @rpath runtime
                # dependencies; rewrite them to @loader_path during bundling.
                bad.append((p, f'{dep} (unresolved @rpath)'))
            elif dep.startswith('/'):
                bad.append((p, dep))
    if bad:
        for p, dep in bad:
            print(f'ERROR: non-portable dependency: {p}: {dep}', file=sys.stderr)
        raise SystemExit(1)


def main() -> int:
    if sys.platform != 'darwin':
        print('Runtime bundling: not macOS; nothing to do')
        return 0
    LIBDIR.mkdir(parents=True, exist_ok=True)
    paths = macho_files()
    for p in paths:
        bundle_for(p)
    for p in LIBDIR.iterdir():
        if p.is_file():
            sign(p)
    for p in paths:
        sign(p)
    validate(paths)
    print(f'Portable macOS runtime: bundled {len(list(LIBDIR.glob("*.dylib")))} dylib(s) in {LIBDIR}')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
