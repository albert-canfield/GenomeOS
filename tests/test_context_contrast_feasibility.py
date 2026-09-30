# SPDX-License-Identifier: AGPL-3.0-or-later
"""Audit B, cached context contrast, checkpoint 1 (lane-contrast, 2026-09-29): the feasibility rules.
Synthetic answers and pairs only: no benchmark file, no cache, no request."""

from __future__ import annotations

import gzip
import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location(
    "context_contrast_feasibility", ROOT / "scripts/context_contrast_feasibility.py"
)
ccf = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = ccf
_SPEC.loader.exec_module(ccf)
ho = ccf.ho
FULL = ccf.FULL_TRACKS


def row(gene="G", n=FULL, **cells):
    return {"gene": gene, "n_tracks": n, "by_cell": dict(cells)}


def answer(eid="EH1", start=100, end=400, genes=None, small=True):
    genes = list(genes or [])
    if small:  # an uncensored answer stores values below 0.05
        genes.append(row("FILLER", K562=0.001, HepG2=0.0, GM12878=-0.002, **{"IMR-90": 0.0}))
    return {"id": eid, "chrom": "chr1", "start": start, "end": end, "tracks": FULL, "genes": genes}


def unit(start=150, end=350, gene="G", gene_id="ENSG1", tss=10_000):
    return ho.Unit(
        source="crispri:Gasperini2019",
        chrom="chr1",
        start=start,
        end=end,
        outcome="not_significant_well_powered",
        gene=gene,
        gene_id=gene_id,
        tss=tss,
        cell="K562",
        split="training",
        study="Gasperini2019",
    )


# --- negatives first -------------------------------------------------------------------------------
def test_only_lane_prior_admissible_endpoints_are_counted():
    for src, part, chroms, _role in ccf.ENDPOINT_SETS:
        assert src in ccf.pot.ADMISSIBLE_SOURCES
        assert ccf.pot.STUDIES[src]["class"] != ccf.pot.INADMISSIBLE
        assert set(chroms) == set(ccf.pot.FRESH if part == "fresh" else ccf.pot.SEEN)
    assert ccf.PRIMARY[:2] == ("crispri:Gasperini2019", "fresh")
    assert [(s, p) for s, p, _, _ in ccf.ENDPOINT_SETS if "Schraivogel" in s] == [
        ("crispri:Schraivogel2020", "seen")
    ]


def test_a_missing_cell_is_never_zero():
    r = row(K562=-0.3, HepG2=-0.1)  # GM12878 and IMR-90 missing
    assert ccf.contrast(r, 1) == pytest.approx(-0.2)
    assert ccf.contrast(r, 2) is None and ccf.contrast(r, 3) is None
    assert ccf.contrast(row(HepG2=-0.1, GM12878=0.0, **{"IMR-90": 0.2}), 1) is None  # no K562: no contrast


def test_an_answer_written_at_the_default_threshold_is_censored():
    censored = answer(genes=[row(K562=-0.3, HepG2=0.2, GM12878=-0.06, **{"IMR-90": 0.07})], small=False)
    assert not ccf.uncensored(censored)
    assert ccf.row_status(censored, censored["genes"][0]) == "answer_censored"
    assert ccf.uncensored(answer(genes=[row(K562=-0.3)]))


def test_rows_short_of_or_above_the_full_track_count_are_not_reliable():
    a = answer()
    base = {"K562": -0.2, "HepG2": 0.0, "GM12878": 0.1, "IMR-90": 0.3}
    assert ccf.row_status(a, row(n=FULL, **base)) == "reliable"
    assert ccf.row_status(a, row(n=2 * FULL, **base)) == "row_merges_two_genes"
    assert ccf.row_status(a, row(n=FULL - 1, **base)) == "row_dropped_exact_zero_tracks"
    no_k562 = {k: v for k, v in base.items() if k != "K562"}
    assert ccf.row_status(a, row(n=FULL, **no_k562)) == "same_cell_missing"


def test_a_pair_with_no_answer_or_no_listed_gene_abstains_with_its_reason():
    assert ccf.pair_record(unit(), [])["reason"] == "no_cached_answer_overlaps"
    far = answer(genes=[row("OTHER", K562=-0.3, HepG2=0.0, GM12878=0.0, **{"IMR-90": 0.0})])
    assert ccf.pair_record(unit(tss=10_000_000), [far])["reason"] == "gene_not_listed_out_of_reach"
    assert ccf.pair_record(unit(tss=1_000), [far])["reason"] == "gene_not_listed_in_reach"
    assert ccf.pair_record(unit(tss=None), [far])["reason"] == "gene_not_listed_no_tss"
    rec = ccf.pair_record(unit(), [far])
    assert rec["any_k"] == -1 and "contrast" not in rec


def test_the_duplicate_rule_fires_when_the_contrast_is_the_same_cell_value():
    recs = []
    for i in range(50):
        v = -0.01 * i
        recs.append(
            {
                "reason": "reliable",
                "contrast": {3: v},
                "same": v,
                "other_median": 0.0,
                "compiled": {"named": False},
            }
        )
    d = ccf.duplication(recs)
    assert d["d2_fires"] and d["duplicate"]


def test_the_duplicate_rule_fires_when_the_compiled_target_determines_the_contrast():
    recs = []
    for i in range(60):
        named = i % 2 == 0
        s = 0.1 + 0.01 * i
        c = -2.0 * s if named else 0.0
        comp = (
            {"named": True, "action": "activates", "strength": s, "cell_is_same": True}
            if named
            else {"named": False}
        )
        recs.append(
            {
                "reason": "reliable",
                "contrast": {3: c},
                "same": (-1) ** i * 0.01 * i,
                "other_median": 0.0,
                "compiled": comp,
            }
        )
    d = ccf.duplication(recs)
    assert d["d1_fires"] and d["d1_r2_on_compiled_target"] == pytest.approx(1.0)


def test_anything_short_of_a_clear_go_is_named():
    sem = {"matched_answers_uncensored_share": 1.0}
    dup = {"share_within_rounding": 0.1, "duplicate": False}
    good = {
        f"{ccf.PRIMARY[0]} {ccf.PRIMARY[1]}": {
            "reliable_k3_with_activity_and_distance": {"clears_floors": True},
            "coverage_reliable_k3_full": 0.9,
            "coverage_reliable_k3_full_contact": 0.5,
        }
    }
    v = ccf.assess(sem, good, dup)
    assert v["clear_go"] and v["reading"] == "go" and not v["contact_on_the_same_pairs"]
    low = {k: {**x, "coverage_reliable_k3_full": 0.6} for k, x in good.items()}
    assert ccf.assess(sem, low, dup)["reading"] == "no_go_completeness"
    assert ccf.assess(sem, good, {**dup, "duplicate": True})["reading"] == "no_go_duplicate"
    assert ccf.assess({"matched_answers_uncensored_share": 0.5}, good, dup)["reading"] == "no_go_semantics"
    assert ccf.assess(sem, good, {**dup, "share_within_rounding": 0.7})["reading"] == "no_go_semantics"


# --- the positive cases ----------------------------------------------------------------------------
def test_the_contrast_is_the_same_cell_value_minus_the_median_of_the_others():
    r = row(K562=-0.5, HepG2=-0.1, GM12878=0.2, **{"IMR-90": -0.4})
    assert ccf.contrast(r, 3) == pytest.approx(-0.5 - (-0.1))
    assert ccf.contrast(row(K562=0.1, HepG2=0.3, GM12878=-0.1), 2) == pytest.approx(0.1 - 0.1)


def test_the_pair_takes_the_reliable_row_of_the_largest_overlap_and_falls_back_to_the_ensembl_id():
    cells = {"HepG2": 0.0, "GM12878": 0.1, "IMR-90": 0.2}
    small = answer("EH_small", 300, 360, [row("G", K562=-0.9, **cells)])
    big = answer("EH_big", 120, 380, [row("ENSG1", K562=-0.3, **cells)])
    unreliable = answer("EH_bad", 140, 360, [row("G", n=FULL - 3, K562=-0.7, **cells)])
    rec = ccf.pair_record(unit(), [small, big, unreliable])
    assert rec["reason"] == "reliable" and rec["element"] == "EH_big"
    assert rec["matched_by_chosen"] == "ensembl_id"
    assert rec["contrast"][3] == pytest.approx(-0.3 - 0.1)
    assert rec["any_k"] == 3


def test_stream_answers_gives_what_json_load_gives(tmp_path):
    data = {f"EH{i}": answer(f"EH{i}", i, i + 10, [row(f"G{i}", K562=-0.01 * i)]) for i in range(40)}
    p = tmp_path / "chr1.json.gz"
    with gzip.open(p, "wt") as fh:
        json.dump(data, fh)
    for chunk in (7, 64, 1 << 20):
        assert dict(ccf.stream_answers(p, chunk=chunk)) == data
    assert [k for k, _ in ccf.stream_answers(p, chunk=5)] == list(data)


def test_stream_answers_refuses_a_truncated_archive(tmp_path):
    p = tmp_path / "chr1.json.gz"
    text = json.dumps({"EH1": answer("EH1"), "EH2": answer("EH2")})
    with gzip.open(p, "wt") as fh:
        fh.write(text[: len(text) - 40])
    with pytest.raises(ValueError):
        list(ccf.stream_answers(p, chunk=16))


def test_a_loose_answer_replaces_the_archive_entry_and_only_overlapping_answers_are_kept(tmp_path):
    arch = {"EH1": answer("EH1", 100, 200), "EH2": answer("EH2", 5_000, 5_100)}
    with gzip.open(tmp_path / "chr1.json.gz", "wt") as fh:
        json.dump(arch, fh)
    (tmp_path / "chr1").mkdir()
    loose = answer("EH1", 100, 200, [row("NEW", K562=-0.2)])
    (tmp_path / "chr1" / "EH1.json").write_text(json.dumps(loose))
    kept, st = ccf.cached_answers("chr1", ccf.overlap_filter([(150, 160)]), cache=tmp_path)
    assert set(kept) == {"EH1"} and kept["EH1"]["_from"] == "loose"
    assert kept["EH1"]["genes"][0]["gene"] == "NEW"
    assert st["archive_answers"] == 2 and st["loose_answers"] == 1


def test_contact_is_read_from_the_cache_only(tmp_path):
    src = ccf.hc.ContactSource(knowledge=tmp_path, readers={"K562": None})
    assert src.contact("K562", "chr1", 1_000, 90_000) is None
    assert src.misses == 1 and src._open == {"K562": None}
    assert not list(tmp_path.iterdir())


# --- added with the reviewer's reading (2026-09-29): wording and locus accounting, the gate unchanged ------
def test_loci_are_counted_on_the_parent_universe_not_recomputed_after_filtering():
    a = unit(start=100, end=200, gene="X", gene_id="EX")
    b = unit(start=5_000_100, end=5_000_200, gene="Y", gene_id="EY")
    c = unit(start=300, end=400, gene="Y", gene_id="EY")  # joins a (same bin) and b (same gene)
    parent = dict(zip([a, b, c], ho.loci([a, b, c]), strict=True))
    fl = ccf.parent_loci(ccf.floors([a, b]), [a, b], parent)
    assert fl["independent_loci_parent_universe"] == 1
    assert fl["loci_recomputed_after_filtering_descriptive"] == 2
    assert "loci" not in fl and fl["clears_floors"] is False


def test_the_notes_never_say_the_feature_cannot_help():
    assert "does not say that a context contrast cannot help" in ccf.READING_NOTES["scope"]
    assert "not registered, not launched" in ccf.READING_NOTES["possible_future_design"]
    assert "consistency check" in ccf.READING_NOTES["reliable_rule_reading"]
    assert "not independent proof" in ccf.READING_NOTES["reliable_rule_reading"]


def test_the_loss_steps_and_the_coverage_wording_follow_the_records():
    us = [unit(start=100 + 10 * i, end=105 + 10 * i, gene=f"G{i}", gene_id=f"E{i}") for i in range(4)]
    base = {
        "activity": True,
        "distance": True,
        "contact": True,
        "contact_zero": False,
        "compiled": {"named": False},
    }
    recs = [
        {**base, "unit": us[0], "reason": "no_cached_answer_overlaps", "any_k": -1, "answers_overlapping": 0},
        {
            **base,
            "unit": us[1],
            "reason": "gene_not_listed_out_of_reach",
            "any_k": -1,
            "answers_overlapping": 1,
        },
        {
            **base,
            "unit": us[2],
            "reason": "row_dropped_exact_zero_tracks",
            "any_k": 3,
            "answers_overlapping": 1,
            "answers_matched": 1,
        },
        {
            **base,
            "unit": us[3],
            "reason": "reliable",
            "any_k": 3,
            "answers_overlapping": 1,
            "answers_matched": 1,
            "contrast": {1: -0.1, 2: -0.1, 3: -0.1},
        },
    ]
    src, part, chroms, role = ccf.PRIMARY
    summary = ccf.with_parent_loci(ccf.set_summary(src, part, chroms, role, us, recs), us, recs)
    assert [x["pairs"] for x in summary["where_coverage_is_lost"]] == [4, 3, 2, 1]
    full = summary["reliable_k3_with_activity_and_distance"]
    assert full["independent_loci_parent_universe"] == 1 and "loci" not in full
    sets = {f"{src} {part}": summary}
    dup = {
        "d1_r2_on_compiled_target": 0.5,
        "named_pairs": 0,
        "pairs": 1,
        "d2_spearman_contrast_same_cell": 0.1,
        "k562_activating_links_with_negative_contrast": "none",
    }
    verdict = {"clear_go": False, "failed": ["no_go_completeness"]}
    notes = ccf.reading_notes(sets, dup, verdict, {"answers_by_model_version": {"unrequested": 3}})
    assert notes["coverage_wording"].startswith(
        "insufficient coverage for this registered complete-case design"
    )
    assert "1 of the 2 pairs that hold a same-cell value" in notes["coverage_wording"]
    assert notes["where_coverage_is_lost_primary"] == [4, 3, 2, 1]
