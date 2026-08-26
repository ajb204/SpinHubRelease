from pathlib import Path


def test_navigation_svg_assets_live_with_plotting_toolbar():
    root = Path(__file__).resolve().parents[1]
    assets = root / "gui" / "plotting" / "assets"
    expected = {
        "fid_to_spectrum.svg", "contours.svg", "redraw_pencil.svg",
        "vertical_trace.svg", "rotate_90.svg", "arrow_up.svg", "one_d.svg",
        "arrow_down.svg", "decon_overlay.svg", "orthogonal.svg", "peaks.svg",
        "horizontal_trace.svg", "pickaxe.svg", "sliders.svg",
    }
    assert assets.is_dir()
    assert expected <= {p.name for p in assets.glob("*.svg")}
    assert not (root / "Frames" / "assets").exists()


def test_toolbar_resolves_assets_beside_canonical_plotting_module():
    root = Path(__file__).resolve().parents[1]
    source = (root / "gui" / "plotting" / "toolbar.py").read_text()
    assert 'Path(__file__).with_name("assets")' in source
