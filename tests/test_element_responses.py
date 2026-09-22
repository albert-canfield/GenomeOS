# SPDX-License-Identifier: AGPL-3.0-or-later
"""The cache-backed reader: a signed value or a named silence, never a silent zero."""

import gzip
import json

from genomeos.attribution.targets import (
    NOT_CACHED,
    NOT_IN_WINDOW,
    NOT_ON_TRACK,
    SCORED,
    ElementResponses,
)


def _gene(name, cells, drop=None):
    g = {"gene": name, "by_cell": dict(cells)}
    if drop is not None:
        g["max_drop_log2fc"] = drop
    return g


def _root(tmp_path):
    archive = {
        "E1": {
            "id": "E1",
            "chrom": "chr1",
            "genes": [
                _gene("TOP", {"K562": -0.40, "GM12878": -0.10}, drop=-0.55),
                _gene("OTHER", {"K562": -0.12}, drop=-0.30),
                _gene("RISES", {"K562": 0.08, "GM12878": 0.01}, drop=-0.02),
            ],
        }
    }
    with gzip.open(tmp_path / "chr1.json.gz", "wt") as fh:
        json.dump(archive, fh)
    loose = tmp_path / "chr2"
    loose.mkdir()
    (loose / "E9.json").write_text(json.dumps({"id": "E9", "genes": [_gene("L", {"K562": -0.2}, -0.2)]}))
    return ElementResponses(tmp_path)


def test_a_gene_the_compact_table_would_not_name_still_has_a_value(tmp_path):
    r = _root(tmp_path)
    got = r.response("chr1", "E1", "OTHER", "K562")
    assert got.value == -0.12 and got.reason == SCORED and got.answered and got
    assert r.value("chr1", "E1", "OTHER", "K562") == -0.12


def test_each_silence_is_named_and_none_of_them_is_zero(tmp_path):
    r = _root(tmp_path)
    outside = r.response("chr1", "E1", "NOWHERE", "K562")
    assert outside.value is None and outside.reason == NOT_IN_WINDOW and not outside.answered
    no_track = r.response("chr1", "E1", "OTHER", "IMR-90")
    assert no_track.value is None and no_track.reason == NOT_ON_TRACK
    missing = r.response("chr1", "E404", "TOP", "K562")
    assert missing.value is None and missing.reason == NOT_CACHED
    uncached_chrom = r.response("chrZ", "E1", "TOP", "K562")
    assert uncached_chrom.value is None and uncached_chrom.reason == NOT_CACHED
    assert r.genes("chr1", "E404") is None


def test_a_scored_zero_is_not_a_silence(tmp_path):
    (tmp_path / "chr3").mkdir()
    (tmp_path / "chr3" / "E0.json").write_text(json.dumps({"genes": [_gene("Z", {"K562": 0.0}, 0.0)]}))
    r = ElementResponses(tmp_path)
    zero = r.response("chr3", "E0", "Z", "K562")
    assert zero.value == 0.0 and zero.reason == SCORED and zero.answered
    assert zero  # truthy: the bool asks "did the sweep answer", never "is the number non-zero"
    assert not r.response("chr3", "E0", "GONE", "K562")


def test_the_ranking_the_compact_table_is_the_head_of(tmp_path):
    r = _root(tmp_path)
    assert r.ranked("chr1", "E1") == [("TOP", -0.55), ("OTHER", -0.30), ("RISES", -0.02)]
    assert [g for g, _ in r.ranked("chr1", "E1", "K562")] == ["TOP", "OTHER", "RISES"]
    assert r.ranked("chr1", "E404") == []
    assert r.response("chr1", "E1", "TOP").value == -0.55  # no cell: the strongest fall anywhere


def test_loose_files_answer_and_one_chromosome_is_held_at_a_time(tmp_path):
    r = _root(tmp_path)
    assert r.value("chr2", "E9", "L", "K562") == -0.2
    assert r.value("chr1", "E1", "TOP", "K562") == -0.40
    assert r.value("chr2", "E9", "L", "K562") == -0.2
    assert r._loaded == ["chr2", "chr1", "chr2"]  # a caller that hops pays a decompression each time
    assert set(r.genes("chr1", "E1")) == {"TOP", "OTHER", "RISES"}
