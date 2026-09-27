# SPDX-License-Identifier: AGPL-3.0-or-later
"""The reader family's normalised readings, registered 2026-09-27 before they were computed.

Three per-biosample readings did not mean what their names say between biosamples: `nodes_open`
is a median split of each biosample against itself, `genes_poised` tracks the H3K27me3 peak call,
and `genes_read` tracks DNase depth. Each now has a normalised counterpart shipped beside it by
`genomeos.genome.reader.normalise_family`, and `node_open_threshold` stays the within-biosample
rank it is. See docs/NODES-READER-WRITER.md, "The reader family normalised, to close
milestone 1.1".
"""

from __future__ import annotations

import gzip
import json
import math
from pathlib import Path

import pytest

from genomeos.genome import reader
from genomeos.genome.domains import Domain

RESULTS = Path(__file__).resolve().parent.parent / "data" / "results"


def test_log_fit_recovers_the_line_and_the_residual_is_standardised():
    xs = [1e4, 2e4, 5e4, 1e5, 3e5]
    ys = [3.0 + 2.0 * math.log(x) for x in xs]
    a, b, s = reader.log_fit(xs, ys)
    assert round(a, 6) == 3.0 and round(b, 6) == 2.0 and s < 1e-9
    # one biosample reads more than its depth predicts: it, and only it, sits clearly above zero
    ys[2] += 5.0
    z = reader.depth_residuals(xs, ys)
    assert max(range(5), key=lambda i: z[i]) == 2 and z[2] > 1.0
    assert abs(sum(z)) < 1e-9


def test_leave_one_out_is_not_orthogonal_to_depth_by_construction():
    """The in-sample residual is uncorrelated with ln(depth) whatever the data, which is why the
    registered depth bar binds on the leave-one-out form."""
    xs = [1e4, 2e4, 4e4, 8e4, 1.6e5, 3.2e5]
    ys = [10.0, 30.0, 45.0, 55.0, 60.0, 62.0]  # saturating: a line in ln(x) is the wrong form
    z_in = reader.depth_residuals(xs, ys)
    lx = [math.log(x) for x in xs]
    mx, mz = sum(lx) / 6, sum(z_in) / 6
    assert abs(sum((u - mx) * (v - mz) for u, v in zip(lx, z_in, strict=True))) < 1e-9
    z_loo = reader.loo_residuals(xs, ys)
    # the deepest point, extrapolated from the other five, is where a wrong form shows most
    assert abs(z_loo[-1]) > abs(z_in[-1])


def _panel(offsets: dict[str, float]) -> dict[str, dict]:
    depth = {"a": 80_000, "b": 150_000, "c": 200_000, "d": 260_000, "e": 330_000, "f": 520_000}
    out = {}
    for c, x in depth.items():
        base = 1000.0 * math.log(x)
        out[c] = {
            "peaks": x,
            "h3k27me3_peaks": x // 10,
            "nodes_open_at_reference": base + offsets.get(c, 0.0),
            "genes_poised": base / 10,
            "genes_read": base + offsets.get(c, 0.0),
            "genes_read_open": base,
        }
    return out


def test_normalise_family_ships_every_registered_reading_beside_an_edge_flag():
    out = reader.normalise_family(_panel({"c": 400.0, "e": -300.0}))
    for name in reader.NORMALISED_READINGS:
        assert all(name in row and name + "_at_edge" in row for row in out.values())
    # depth alone produces no difference between biosamples
    assert all(abs(row["genes_read_open_depth_residual"]) < 1e-6 for row in out.values())
    # what depth does not predict is what the reading keeps, in the right order
    z = {c: row["genes_read_depth_residual"] for c, row in out.items()}
    assert max(z, key=z.get) == "c" and min(z, key=z.get) == "e"
    # the deepest and shallowest biosamples are flagged as extrapolations, nobody else
    edge = {c for c, row in out.items() if row["genes_read_depth_residual_at_edge"]}
    assert edge == {"a", "f"}


def test_a_panel_with_a_hole_normalises_nothing_rather_than_a_different_panel():
    p = _panel({})
    p["b"]["h3k27me3_peaks"] = None
    out = reader.normalise_family(p)
    assert all(row["genes_poised_mark_residual"] is None for row in out.values())
    assert all(row["genes_read_depth_residual"] is not None for row in out.values())


def test_reference_density_is_frozen_and_the_median_rank_is_untouched():
    assert reader.NODE_OPEN_REFERENCE_DENSITY == 4.79
    # node_open_threshold stays the within-biosample median rank with its 1.0 floor
    assert reader.node_open_threshold([1.0, 3.0, 5.0, 7.0]) == 5.0
    assert reader.node_open_threshold([0.1, 0.2, 0.3]) == reader.NODE_OPEN_MIN_DENSITY
    assert reader.NORMALISED_READINGS["genes_poised_mark_residual"] == ("genes_poised", "h3k27me3_peaks")


def _read(tmp_path, monkeypatch, per_node: list[int]) -> dict:
    monkeypatch.setattr(reader, "RESULTS", tmp_path)
    with gzip.open(tmp_path / "dnase_S_c.bed.gz", "wt") as fh:
        fh.write("# test\n")
        for i, n in enumerate(per_node):
            for k in range(n):
                s = i * 100_000 + 1_000 + k * 5_000
                fh.write(f"{s}\t{s + 150}\t1.0\n")
    ann = type("A", (), {"genes": {}})()
    domains = [Domain(f"c:D{i}", "c", i * 100_000, (i + 1) * 100_000) for i in range(len(per_node))]
    return reader.read_chromosome("S", "c", ann, domains, [], marks=False)


def test_deeper_file_leaves_the_median_split_alone_and_moves_the_absolute_count(tmp_path, monkeypatch):
    """The defect and the repair side by side: doubling every node's peaks keeps `nodes_open` at
    half the nodes, while `nodes_open_at_reference` rises, which is why the reading shipped
    between biosamples is its depth residual and not the count."""
    shallow = _read(tmp_path, monkeypatch, [1, 3, 5, 7])
    deep = _read(tmp_path, monkeypatch, [2, 6, 10, 14])
    assert shallow["nodes_open"] == deep["nodes_open"] == 2
    assert shallow["nodes_open_at_reference"] == 2 and deep["nodes_open_at_reference"] == 3
    assert shallow["h3k27me3_peaks"] is None  # no mark read, no covariate claimed
    assert "nodes_open_depth_residual" in shallow["evidence"]["depth"]


_NORM = RESULTS / "reader_normalised.json"
_RAR = RESULTS / "reader_rarefied.json"


@pytest.mark.skipif(not _NORM.exists() or not _RAR.exists(), reason="result files not on disk")
def test_the_published_share_reverses_on_the_residual_and_survives_at_matched_depth():
    """Both results are kept visible: the registered residual reverses K562 against HepG2 on read
    share, and cutting both to the same number of peaks keeps K562 above with the gap shrunk. K562
    is alone at the top of DNase depth, so its residual is an extrapolation."""
    norm = json.loads(_NORM.read_text())
    assert all(norm["kept"].values())
    for f in ("genes_read_open_depth_residual", "genes_read_depth_residual"):
        assert norm["k562_vs_hepg2_first"][f]["verdict"] == "reverses"
    assert norm["per_biosample"]["K562"]["genes_read_depth_residual_at_edge"] is True
    runs = json.loads(_RAR.read_text())["runs"]
    at_hepg2 = runs["176634"]["agreement"]
    assert at_hepg2["genes_read_open"]["k562_vs_hepg2"] == "K562 > HepG2"
    assert at_hepg2["genes_read"]["k562_vs_hepg2"] == "K562 > HepG2"
    # the node reading's reversal is the one matched depth confirms
    assert at_hepg2["nodes_open_at_reference"]["k562_vs_hepg2"] == "HepG2 >= K562"
    gap_published = (14_828 - 13_015) / 20_094
    gap_matched = (at_hepg2["genes_read_open"]["k562"] - at_hepg2["genes_read_open"]["hepg2"]) / 20_094
    assert 0 < gap_matched < gap_published / 2
