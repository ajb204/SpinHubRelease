"""Helpers for reading and updating legacy parameter files.

The project stores settings as loose text lines of the form:

    key = value optional trailing tokens

These helpers preserve the original file format while centralising the logic
used by several GUI windows.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping


@dataclass(frozen=True)
class ParameterMatch:
    key: str
    value: str
    tail: str = ""


def _read_text_lines(path: str | Path) -> list[str]:
    file_path = Path(path)
    if not file_path.exists():
        return []
    return file_path.read_text().splitlines(keepends=True)


def _split_legacy_line(line: str) -> tuple[list[str], str]:
    """Return whitespace tokens plus the original newline suffix."""
    newline = ""
    if line.endswith("\r\n"):
        newline = "\r\n"
        line = line[:-2]
    elif line.endswith("\n"):
        newline = "\n"
        line = line[:-1]
    return line.split(), newline


def _extract_match(line: str, key: str) -> ParameterMatch | None:
    tokens, _ = _split_legacy_line(line)
    if len(tokens) < 2:
        return None
    if tokens[0] != key or tokens[1] != "=":
        return None
    value = tokens[2] if len(tokens) >= 3 else ""
    tail = " ".join(tokens[3:]) if len(tokens) > 3 else ""
    return ParameterMatch(key=key, value=value, tail=tail)

def parse_value(path: str | Path, key: str, default: Any = "") -> Any:
    for line in _read_text_lines(path):
        match = _extract_match(line, key)
        if match is not None:
            return match.value if match.tail == "" else f"{match.value} {match.tail}".strip()
    return default


def parse_float(path: str | Path, key: str, default: float = 0.0) -> float:
    for line in _read_text_lines(path):
        match = _extract_match(line, key)
        if match is not None:
            try:
                return float(match.value)
            except Exception:
                return default
    return default


def parse_int(path: str | Path, key: str, default: int = 0) -> int:
    for line in _read_text_lines(path):
        match = _extract_match(line, key)
        if match is not None:
            try:
                return int(float(match.value))
            except Exception:
                return default
    return default


def parse_all_strings(path: str | Path, key: str) -> str | int:
    values: list[str] = []
    for line in _read_text_lines(path):
        tokens, _ = _split_legacy_line(line)
        if len(tokens) > 2 and tokens[0] == key and tokens[1] == "=":
            values.extend(tokens[2:])
    if not values:
        return 0
    return " ".join(values) + " "


def _format_line(key: str, value: Any, tail: str) -> str:
    if tail:
        return f"{key} = {value} {tail}\n"
    return f"{key} = {value}\n"


def update_parameter_file(
    target_path: str | Path,
    updates: Mapping[str, Any],
    *,
    source_path: str | Path | None = None,
) -> Path:
    """Update a legacy parameter file in place.

    Existing lines are preserved. Matching keys keep their trailing tokens,
    comments, and any legacy formatting after the value.
    """

    target = Path(target_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    source = Path(source_path) if source_path is not None else target
    existing_lines = _read_text_lines(source)
    remaining = dict(updates)
    updated_keys = set(updates)
    output_lines: list[str] = []
    replaced_keys: set[str] = set()

    for line in existing_lines:
        tokens, newline = _split_legacy_line(line)
        if len(tokens) >= 2 and tokens[0] in updated_keys and tokens[1] == "=":
            key = tokens[0]
            if key in replaced_keys:
                continue
            if key in remaining:
                value = remaining.pop(key)
                # Preserve the original line when the stored value is already
                # identical.  Besides avoiding needless file churn, this makes
                # save operations genuinely change only settings that changed.
                if len(tokens) >= 3 and tokens[2] == str(value):
                    output_lines.append(line)
                    replaced_keys.add(key)
                    continue
                tail = " ".join(tokens[3:]) if len(tokens) > 3 else ""
                output_lines.append(_format_line(key, value, tail).rstrip("\n") + newline)
                replaced_keys.add(key)
            continue
        output_lines.append(line)

    for key, value in remaining.items():
        output_lines.append(f"{key} = {value}\n")

    new_text = "".join(output_lines)
    if target.exists():
        try:
            if target.read_text() == new_text:
                return target
        except Exception:
            pass
    target.write_text(new_text)
    return target

def remove_parameter_keys(target_path: str | Path, keys) -> Path:
    """Remove named scalar parameters from a legacy parameter file.

    Used for settings whose current value equals the application default: the
    absence of the key records that the default, rather than an override, is in
    effect.
    """
    target = Path(target_path)
    if not target.exists():
        return target
    remove = set(keys)
    output_lines = []
    for line in _read_text_lines(target):
        tokens, _newline = _split_legacy_line(line)
        if len(tokens) >= 2 and tokens[0] in remove and tokens[1] == "=":
            continue
        output_lines.append(line)
    new_text = "".join(output_lines)
    if target.read_text() != new_text:
        target.write_text(new_text)
    return target


def write_structured_parameter_file(target_path: str | Path, data: Mapping[str, Any]) -> Path:
    """Write the structured legacy session files used by CATIA/CPMG/Decay views."""

    target = Path(target_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []

    for key, vals in data.items():
        if isinstance(vals, Mapping):
            if vals and all(isinstance(v, Mapping) for v in vals.values()):
                for outer_key, inner_map in vals.items():
                    for inner_key, inner_val in inner_map.items():
                        lines.append(f"{key} {outer_key} {inner_key} {inner_val}\n")
            else:
                for outer_key, value in vals.items():
                    if isinstance(value, (tuple, list)):
                        joined = " ".join(str(v) for v in value)
                        lines.append(f"{key} {outer_key} {joined}\n")
                    else:
                        lines.append(f"{key} {outer_key} {value}\n")
        else:
            lines.append(f"{key} {vals}\n")

    target.write_text("".join(lines))
    return target

def read_structured_parameter_file(path: str | Path) -> dict[str, Any]:
    """Parse the structured legacy session files used by CATIA/CPMG/Decay views."""

    parsed: dict[str, Any] = {}
    for line in _read_text_lines(path):
        tokens, _ = _split_legacy_line(line)
        if not tokens:
            continue

        key = tokens[0]
        if len(tokens) == 2:
            parsed[key] = tokens[1]
            continue

        if key not in parsed or not isinstance(parsed.get(key), dict):
            parsed[key] = {}
        section = parsed[key]

        if len(tokens) == 3:
            section[tokens[1]] = tokens[2]
            continue

        if key == 'par' and len(tokens) >= 4:
            entry = section.setdefault(tokens[1], {})
            entry[tokens[2]] = ' '.join(tokens[3:])
        elif key == 'dataset' and len(tokens) >= 4:
            section[tokens[1]] = tuple(tokens[2:])
        else:
            section[tokens[1]] = ' '.join(tokens[2:])

    return parsed

