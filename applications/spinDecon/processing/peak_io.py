"""Peak-list file I/O shared by GUI workspaces.

This module intentionally has no wx dependency. The Full Peak List remains the
application authority; this helper only parses an explicit peak-list file for
import/review workflows.
"""
from __future__ import annotations

from spinDecon.domain.peaks import peakEntry


def read_peak_list(infile):
    peaks = []
    with open(infile, "r") as peakfile:
        for line in peakfile:
            fields = line.split()
            if not fields:
                continue
            try:
                float(fields[1])
            except (IndexError, TypeError, ValueError):
                continue
            peaks.append(peakEntry(fields))
    return peaks


readpeaklist = read_peak_list
