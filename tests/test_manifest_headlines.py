"""Review item R9, second pass: the results behind the README's headline figures carry a complete
manifest, were rebuilt in a clean checkout, and still say what the docs quote."""

import json
import sys
from pathlib import Path

import pytest

from genomeos import manifest as mf

RESULTS = Path("data/results")
RECORD = json.loads((RESULTS / "manifest_headlines.json").read_text())
REBUILT = {r["result"]: r for r in RECORD["rebuilt"]}


def load(name: str) -> dict:
    return json.loads((RESULTS / f"{name}.json").read_text())


@pytest.mark.parametrize("name", sorted(REBUILT))
def test_a_headline_result_carries_a_complete_manifest_from_a_clean_checkout(name):
    m = load(name)[mf.KEY]
    assert m["complete"] and mf.validate({k: v for k, v in m.items() if k in mf.REQUIRED}) == []
    assert m["code"]["dirty"] is False and m["code"]["git_sha"] == REBUILT[name]["rebuilt_at"]
    assert all(len(i["sha256"]) == 64 and "partition" in i for i in m["inputs"])


@pytest.mark.parametrize("name", sorted(REBUILT))
def test_a_headline_is_recorded_as_reproduced(name):
    r = REBUILT[name]
    assert r["status"] in ("reproduced", "reproduced_figures") and r["alphagenome_requests"] == 0
    assert r["status"] == "reproduced" or r["difference_explained"]
    assert r["second_rerun_differences"] == [] or r["status"] == "reproduced_figures"


def test_the_quoted_figures_are_the_ones_on_disk():
    assert load("node_containment_audit")["modelled"]["controls"]["uniform"]["excess_points"] == 2.898
    assert load("node_containment_measured")["measured"]["controls"]["uniform"]["excess_points"] == 5.893
    ctl = load("constrained_unknown_targets")["matched_random_control"]["real_unknown"]["moves_a_gene"]
    assert (ctl["block_rate"], ctl["random_window_rate"]) == (0.6234, 0.8604)
    tb = load("therapeutic_benchmark")
    assert (tb["targets_recovered"], tb["verdicts_correct"], tb["top_mechanism_defensible"]) == (9, 9, 9)
    ru = load("unknown_coverage")["real_unknown"]
    assert (ru["measured_bp"], ru["bp"]) == (160_447, 30_602_182)


def test_crispri_inputs_carry_the_benchmark_split():
    from genomeos.attribution.measured import CRISPRI_SPLIT_OF

    for name in ("node_containment_measured", "unknown_coverage"):
        got = {Path(i["path"]).name: i["partition"] for i in load(name)[mf.KEY]["inputs"]}
        assert {n: got[n] for n in CRISPRI_SPLIT_OF} == CRISPRI_SPLIT_OF


def test_the_first_run_check_reads_history_not_the_file_on_disk():
    sys.path.insert(0, "scripts")
    import constrained_unknown_targets as cut

    first = cut.first_run()
    if first is None:
        pytest.skip("no git history here")
    assert first["date"] == "2026-09-16" and "matched_random_control" not in first
    assert load("constrained_unknown_targets")["reproduced_from_2026_09_16"]["all"] is True


def test_crispri_published_is_named_as_pending_until_it_is_rebuilt():
    pending = {p["result"] for p in RECORD["pending"]}
    assert "crispri_published" in pending or "crispri_published" in REBUILT
