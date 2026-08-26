#!/usr/bin/env python3
from __future__ import annotations
import os, platform, shutil, subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VENV = ROOT / '.venv'
BIN = ROOT / 'bin'
LIBS = ROOT / 'extern' / 'libs'
SYSTEM_PREFIXES = ('/System/Library/', '/usr/lib/')
FORBIDDEN_PREFIXES = ('/opt/homebrew/', '/usr/local/', '/opt/local/', '/sw/')


def mark(ok): return 'OK' if ok else 'MISSING'
def check_file(label, path, executable=False):
    ok = path.exists() and (not executable or os.access(path, os.X_OK))
    print(f'[{mark(ok):7}] {label}: {path}')
    return ok

def macho(path: Path) -> bool:
    return 'Mach-O' in subprocess.run(['file', str(path)], capture_output=True, text=True).stdout

def deps(path: Path) -> list[str]:
    out = subprocess.run(['otool', '-L', str(path)], capture_output=True, text=True, check=True).stdout.splitlines()[1:]
    return [x.strip().split(' (compatibility version', 1)[0] for x in out if x.strip()]

def local_dep_exists(owner: Path, dep: str) -> bool:
    if dep.startswith('@loader_path/'):
        return (owner.parent / dep[len('@loader_path/'):]).resolve().exists()
    if dep.startswith('@executable_path/'):
        return (BIN / dep[len('@executable_path/'):]).resolve().exists()
    return False

def portable_runtime_check() -> tuple[bool, list[str]]:
    problems: list[str] = []
    candidates = [p for p in BIN.iterdir() if p.is_file() and not p.is_symlink() and macho(p)]
    libdir = BIN / 'lib'
    if libdir.exists():
        candidates.extend(p for p in libdir.iterdir() if p.is_file() and macho(p))
    for p in candidates:
        for dep in deps(p):
            if dep.startswith(SYSTEM_PREFIXES) or dep.startswith('/System/Library/'):
                continue
            if dep.startswith(FORBIDDEN_PREFIXES):
                problems.append(f'{p.relative_to(BIN)}: forbidden absolute dependency {dep}')
            elif dep.startswith('@rpath/'):
                problems.append(f'{p.relative_to(BIN)}: unresolved runtime dependency {dep}')
            elif dep.startswith('@loader_path/') or dep.startswith('@executable_path/'):
                if not local_dep_exists(p, dep):
                    problems.append(f'{p.relative_to(BIN)}: bundled dependency missing: {dep}')
            elif dep.startswith('/'):
                problems.append(f'{p.relative_to(BIN)}: non-system absolute dependency {dep}')
    return not problems, problems

def main():
    print(f'spin suite doctor - {platform.system()} {platform.machine()}')
    checks=[]
    py = VENV / 'bin' / 'python'
    checks.append(check_file('venv Python', py, True))
    for tool in ('make','cc','c++','gfortran'):
        ok=shutil.which(tool) is not None; checks.append(ok); print(f'[{mark(ok):7}] tool: {tool}')
    if py.exists():
        code='import importlib.util; mods='+repr(['numpy','scipy','matplotlib','wx','nmrglue','texttable','pymupdf'])+'; print("\\n".join(f"{m}:{bool(importlib.util.find_spec(m))}" for m in mods))'
        out=subprocess.run([str(py),'-c',code],text=True,capture_output=True)
        for line in out.stdout.splitlines():
            name,val=line.rsplit(':',1); ok=val=='True'; checks.append(ok); print(f'[{mark(ok):7}] Python module: {name}')
    for lib in ('libfftw3.a','libfftw3_threads.a','libfftw3f.a','libfftw3f_threads.a','liblbfgs.a'):
        checks.append(check_file(lib, LIBS/lib))
    checks.append(check_file('spinDecon launcher', BIN/'spinDecon'))
    checks.append(check_file('spinHub launcher', BIN/'spinHub'))
    if platform.system() == 'Darwin':
        ok, problems = portable_runtime_check()
        checks.append(ok); print(f'[{mark(ok):7}] portable macOS dylib dependencies')
        for problem in problems: print(f'          {problem}')
    if py.exists():
        env=os.environ.copy(); env['PYTHONPATH']=str(ROOT/'applications')
        p=subprocess.run([str(py),'-c','import spinDecon; import spinHubMain'],env=env,capture_output=True,text=True)
        ok=p.returncode==0; checks.append(ok); print(f'[{mark(ok):7}] import spinDecon + spinHubMain')
        if not ok and p.stderr: print('          '+p.stderr.splitlines()[-1])
    print(f'\nResult: {sum(checks)}/{len(checks)} checks passed')
    return 0 if all(checks) else 1

if __name__ == '__main__': raise SystemExit(main())
