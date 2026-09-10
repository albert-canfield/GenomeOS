"""Task 4.2: segment count and period follow the clock-and-wavefront relation."""

from pathlib import Path

from genomeos.runtime.segmentation import run_segmentation


def test_segments_form_at_the_clock_period():
    r = run_segmentation(length=120, hours=100, dt=0.1)
    assert r.period_h == 5.0
    assert r.frozen_cells > 30
    # segments in the frozen region ≈ cells / (speed × period), allowing one edge effect
    assert abs(r.count - r.expected) <= 1.5, (r.count, r.expected, r.segments)
    widths = [e - s for s, e in r.segments[1:-1]]
    assert widths and max(widths) - min(widths) <= 2  # regular, like somites
    unc = r.uncertainty().to_dict()
    assert unc["molecular"]["label"] in ("high", "medium")
    assert "wavefront_speed" in [p for p in r.module.parameters]


def test_slower_clock_gives_fewer_segments(tmp_path):
    src = Path("data/demo/segmentation_clock.bio").read_text()
    slow_src = src.replace("param period_h = 5.0 h", "param period_h = 10.0 h")
    assert slow_src != src
    p = tmp_path / "slow.bio"
    p.write_text(slow_src)
    slow = run_segmentation(p, length=120, hours=100, dt=0.1)
    fast = run_segmentation(length=120, hours=100, dt=0.1)
    assert slow.period_h == 10.0
    assert fast.count > slow.count * 1.5, (fast.count, slow.count)
