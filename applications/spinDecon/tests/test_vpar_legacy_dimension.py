import ast
from pathlib import Path


def test_vpar_has_no_direct_ordering_of_legacy_dim_against_numbers():
    source_path = Path(__file__).parents[1] / 'processing' / 'vpar_decon.py'
    tree = ast.parse(source_path.read_text())
    offenders = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Compare):
            continue
        text = ast.get_source_segment(source_path.read_text(), node) or ''
        if 'self.dim' not in text:
            continue
        if any(isinstance(op, (ast.Gt, ast.GtE, ast.Lt, ast.LtE)) for op in node.ops):
            offenders.append((node.lineno, text))
    assert offenders == []
