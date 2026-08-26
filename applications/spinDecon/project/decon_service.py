"""Shared launcher for the decon external binary.

This keeps init-file writing and process launch in one place so the GUI tabs
only assemble run specifications.
"""
from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass
from typing import Iterable, Mapping, Optional, Sequence, Union


@dataclass
class DeconLaunchResult:
    returncode: int
    command: str
    cwd: str


class DeconService:
    def __init__(self, default_cwd: Optional[str] = None):
        self.default_cwd = default_cwd or os.getcwd()

    @staticmethod
    def init_path(filename: str = "decon.init", cwd: Optional[str] = None) -> str:
        return os.path.join(cwd or os.getcwd(), filename)

    def write_init_dict(
        self,
        values: Mapping[str, object],
        filename: str = "decon.init",
        cwd: Optional[str] = None,
    ) -> str:
        path = self.init_path(filename, cwd=cwd)
        with open(path, "w", encoding="utf-8") as handle:
            for key, value in values.items():
                handle.write(f"{key}\t{value}\n")
        return path

    def write_init_lines(
        self,
        lines: Iterable[str],
        filename: str = "decon.init",
        cwd: Optional[str] = None,
    ) -> str:
        path = self.init_path(filename, cwd=cwd)
        with open(path, "w", encoding="utf-8") as handle:
            for line in lines:
                handle.write(f"{line.rstrip()}\n")
        return path

    def launch(
        self,
        command: Union[str, Sequence[str]],
        cwd: Optional[str] = None,
        background: bool = False,
    ):
        workdir = cwd or self.default_cwd
        if isinstance(command, str):
            if background:
                return subprocess.Popen(
                    command,
                    cwd=workdir,
                    shell=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                )
            completed = subprocess.run(command, cwd=workdir, shell=True)
            return DeconLaunchResult(
                returncode=completed.returncode,
                command=command,
                cwd=workdir,
            )

        if background:
            return subprocess.Popen(command, cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        completed = subprocess.run(command, cwd=workdir)
        return DeconLaunchResult(
            returncode=completed.returncode,
            command=" ".join(command),
            cwd=workdir,
        )

    def run(
        self,
        command: Union[str, Sequence[str]],
        values: Optional[Mapping[str, object]] = None,
        lines: Optional[Iterable[str]] = None,
        filename: str = "decon.init",
        cwd: Optional[str] = None,
        background: bool = False,
    ):
        if values is not None:
            self.write_init_dict(values, filename=filename, cwd=cwd)
        elif lines is not None:
            self.write_init_lines(lines, filename=filename, cwd=cwd)
        return self.launch(command, cwd=cwd, background=background)
