"""Shared helpers for Process-family NMR dimension labels.

The ProcessFrame owns live label state.  This module contains only non-GUI
normalisation/discovery helpers so child windows never need to initialise the
label store themselves.
"""
from __future__ import annotations

import os
import re


def clean_dimension_label(value) -> str:
    text = '' if value is None else str(value)
    text = text.strip().replace(' ', '').replace('<', '').replace('>', '')
    return '' if text in ('', '0', 'None') else text


def canonical_spectral_labels(labels: list[str]) -> list[str]:
    """Disambiguate repeated raw labels using the historical _1/_2 scheme."""
    out: list[str] = []
    for raw in labels:
        label = clean_dimension_label(raw)
        if not label:
            label = 'H1'
        matches = [i for i, existing in enumerate(out)
                   if existing == label or re.fullmatch(re.escape(label) + r'_\d+', existing)]
        if matches:
            # Convert the first unsuffixed occurrence and number this occurrence.
            first = matches[0]
            if out[first] == label:
                out[first] = label + '_1'
            label = label + '_' + str(len(matches) + 1)
        out.append(label)
    return out


def discover_bruker_labels(raw_dir: str, count: int) -> list[str]:
    """Read NUC1 from acqus/acquNs without importing a GUI frame."""
    names = ('acqus', 'acqu2s', 'acqu3s', 'acqu4s')
    labels: list[str] = []
    for idx in range(min(max(int(count), 0), 4)):
        path = os.path.join(raw_dir or '', names[idx])
        value = ''
        try:
            with open(path, 'r', errors='replace') as handle:
                for line in handle:
                    if line.startswith('##$NUC1='):
                        value = clean_dimension_label(line.split('=', 1)[1])
                        break
        except OSError:
            pass
        labels.append(value)
    return labels
