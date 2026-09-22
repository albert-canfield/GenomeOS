# SPDX-License-Identifier: AGPL-3.0-or-later
"""The non-coding frame: membership is the fourth frame's own branch, and n is 2, not 14."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomeos.benchmark import loci_noncoding as nc

RESULTS = Path("data/results")
FOURTH = RESULTS / "loci_fourth.json"
BRANCH_REASON = "a regulated target is not protein coding in GENCODE"

needs_knowledge = pytest.mark.skipif(
    not (nc.GENE_TABLE.exists() and (nc.crispri.KNOWLEDGE / nc.SOURCE).exists()),
    reason="the GENCODE table or the held-out CRISPR benchmark file is not fetched",
)
needs_fourth = pytest.mark.skipif(not FOURTH.exists(), reason="the fourth frame has not been run")


@pytest.fixture(scope="module")
def branch():
    return nc.dropped_for_non_coding()


@needs_knowledge
@needs_fourth
def test_membership_is_the_fourth_frames_own_branch(branch):
    """The re-derivation must reproduce the committed verdicts exactly, element for element.

    This is the whole claim that no curation happened here: if the cached GENCODE table and the
    per-chromosome GFF3 `loci_fourth` reads ever disagree, this fails rather than quietly drawing a
    different frame.
    """
    committed = {
        (v["chrom"], v["element"][0], v["element"][1], v["cell"], tuple(v["targets"]))
        for v in json.loads(FOURTH.read_text())["draw"]["verdicts"]
        if (v["reason"] or "").startswith(BRANCH_REASON)
    }
    mine = {(r["chrom"], r["start"], r["end"], r["cell"], tuple(r["targets"])) for r in branch}
    assert mine == committed
    assert len(mine) == 14


@needs_knowledge
def test_the_ensembl_join_finds_twelve_stale_symbols(branch):
    split = nc.classify(nc.resolve(branch))
    assert len(split["stale_symbol"]) == 12
    assert len(split["non_coding"]) == 2
    stale = {x["published_symbol"]: x["current_symbol"] for r in split["stale_symbol"] for x in r["resolved"]}
    assert stale == {"SSFA2": "ITPRID2", "SARS": "SARS1", "WDR61": "SKIC8"}
    kept = {t for r in split["non_coding"] for t in r["non_coding_targets"]}
    assert kept == {"LINC00885", "CCDC26"}


@needs_knowledge
def test_every_resolved_target_is_resolvable(branch):
    """An unresolvable Ensembl id would make the split a guess. None of the fourteen is."""
    for r in nc.resolve(branch):
        for x in r["resolved"]:
            assert x["ensembl_id"], x
            assert x["gene_type"], x


@needs_knowledge
def test_the_non_coding_two_are_lncrna(branch):
    split = nc.classify(nc.resolve(branch))
    types = {x["published_symbol"]: x["gene_type"] for r in split["non_coding"] for x in r["resolved"]}
    assert types == {"LINC00885": "lncRNA", "CCDC26": "lncRNA"}


@needs_knowledge
def test_the_registration_states_n_twice_and_a_zero_budget():
    p = nc.PREREGISTRATION
    assert p["n_as_the_branch_reports_it"] == 14
    assert p["n_after_the_ensembl_join"] == 2
    assert p["request_budget"] == 0
    assert nc.REQUEST_BUDGET == 0
    assert "NOT the nearest coding TSS" in p["the_right_baseline"]
    assert "known_before_the_registration" in p
    assert "failure" in p["outcomes"]


def test_plan_refuses_to_exceed_the_registered_budget(monkeypatch):
    """The budget is zero, so a locus that needed a request must abort the run, not fund it."""
    monkeypatch.setattr(
        nc.loci,
        "read_reach",
        lambda ch, e: {"askable": True, "targets_in_reach": [], "targets_out_of_reach": []},
    )
    monkeypatch.setattr(nc.loci, "_deletion_rows", lambda *a, **k: [])
    monkeypatch.setattr(nc.loci, "Chromosome", lambda *a, **k: type("C", (), {"close": lambda self: None})())
    e = nc.Expect(
        locus="X",
        chrom="chr1",
        element=(10, 20),
        window=(0, 100),
        classes=("program",),
        targets=("X",),
        direction="activates",
        tissues=(),
    )
    with pytest.raises(RuntimeError, match="registered a budget of 0"):
        nc.plan((e,))


def test_counterfactual_drops_the_coding_first_break(monkeypatch):
    """The registered baseline must rank by effect across BOTH keys, not stop at the coding one."""
    rows = [
        {
            "predicted": {"gene": "LNC", "log2_fold_change": -1.7},
            "predicted_coding": {"gene": "COD", "log2_fold_change": -0.14},
        }
    ]
    monkeypatch.setattr(nc.loci, "_deletion_rows", lambda *a, **k: rows)
    e = nc.Expect(
        locus="L",
        chrom="chr8",
        element=(10, 20),
        window=(0, 100),
        classes=("program",),
        targets=("LNC",),
        direction="activates",
        tissues=(),
    )
    out = nc.counterfactual((e,))
    assert out["k"] == 1 and out["n"] == 1
    assert out["per_locus"][0]["would_name"] == "LNC"


def test_read_deletion_no_longer_hides_a_non_coding_answer_behind_a_coding_one(monkeypatch):
    """The defect, now fixed, and the history kept because the history is what licenses the fix.

    Until 2026-09-21 this test asserted the OPPOSITE - ``out["target"] == "COD"`` and LNC absent from
    the ranking - and it was written that way on purpose: `loci.read_deletion` walked
    ``("predicted_coding", "predicted")`` and broke at the first key that held a gene, so a deletion
    whose strongest effect was on a non-coding gene was reported as whatever coding gene came second.
    Section 18 of docs/LOCI-BENCHMARK.md found that at the H19 ICR as a caveat that changed no
    verdict. Section 23 caught it with a positive - the chr8 CCDC26 element, where the row below is
    the real one: the published lncRNA at -1.691 against a coding gene 204 kb away at -0.1433, a
    twelfth of the effect - and pinned the defective behaviour here so that whoever repaired it would
    be told that the section had to be re-read rather than left standing.

    It has been repaired, and section 23 was re-read: section 24 carries the before and after of
    every frame, and `loci_reread.PREREGISTRATION` was committed before any of those numbers existed.
    Both keys are now ranked together, so the reader returns LNC; `coding_first_target` keeps the old
    answer beside it so the rates published before that date stay computable.
    """
    rows = [
        {
            "predicted": {"gene": "LNC", "log2_fold_change": -1.7, "action": "activates", "tissue": "t"},
            "predicted_coding": {
                "gene": "COD",
                "log2_fold_change": -0.14,
                "action": "activates",
                "tissue": "t",
            },
            "run": "sweep",
        }
    ]
    monkeypatch.setattr(nc.loci, "_deletion_rows", lambda *a, **k: rows)
    ch = type("C", (), {"chrom": "chr8"})()
    out = nc.loci.read_deletion(ch, 10, 20)
    assert out["target"] == "LNC"
    assert [g["gene"] for g in out["targets"]] == ["LNC", "COD"]
    # the element counts once towards each gene, not twice towards the one both keys could name
    assert [g["elements"] for g in out["targets"]] == [1, 1]
    # and the reading the benchmark published before 2026-09-21 is kept beside the ranking
    assert out["coding_first_target"] == "COD"
    assert [g["gene"] for g in out["coding_first_targets"]] == ["COD"]
