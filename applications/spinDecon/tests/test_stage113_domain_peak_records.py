from pathlib import Path

from spinDecon.domain.peaks import peakEntry


def test_canonical_peak_records_are_available():
    assert peakEntry is not None


def test_active_code_uses_domain_peak_records():
    root = Path(__file__).resolve().parents[1]
    active = [root / p for p in ("gui", "integrations", "analysis", "project", "workflow", "processing")]
    offenders = []
    for base in active:
        for path in base.rglob("*.py"):
            text = path.read_text(errors="ignore")
            if "decon.misc.peaks" in text:
                offenders.append(str(path.relative_to(root)))
    assert offenders == []
