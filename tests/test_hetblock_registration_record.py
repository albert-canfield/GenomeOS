# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-hetblock's registration record, guarded by the test that ships WITH it and never before it.

What is guarded here: that data/results/hetblock_registration.json was written only AFTER the harness
read HARNESS_REPRODUCES; that every BOUND field in it is the committed module's field VERBATIM, so the
run cannot differ from what was fixed before any value of this lane was read; that the one change from
lane-scalehet's test is the exchangeable unit, with the statistic, the secondary statistic, the draw
count, the seed, MIN_GENES, MIN_CHROMOSOMES, the 98.33% interval and the tolerance all scalehet's or
magcal's; that ALPHA is the LAXER 0.05 with 0.0167 only reported; that the exposure is stated twice
and names the nine readings and the post-hoc spread already tested and rejected once; that BOTH
outcomes are written out here, before either is seen, and the withdrawal is not softened; that no
reading proposes a regional model; and that NO STATISTIC OF THIS LANE is in the file.

This file reads data/results/hetblock_registration.json and is committed in the same commit as it
(`22e12de`). No pytest.skip is keyed on that record (`4f2652c`).
"""

from __future__ import annotations

import json
from pathlib import Path

import tracked_paths as tp

from genomeos.attribution import hetblock as hb
from genomeos.attribution import magcal as mc
from genomeos.attribution import scalehet as sh

ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / "data/results/hetblock_registration.json"
HARNESS_RECORD = ROOT / "data/results/hetblock_harness.json"


def as_json(value: object) -> object:
    """tuples become lists through json.dump: normalise, so VERBATIM means verbatim and not 'close'."""
    if isinstance(value, (list, tuple)):
        return [as_json(v) for v in value]
    if isinstance(value, dict):
        return {k: as_json(v) for k, v in value.items()}
    return value


def record() -> dict:
    tp.must_be_present(RECORD)
    return json.loads(RECORD.read_text())


def test_the_registration_is_tracked_and_is_step_three_of_this_lane() -> None:
    d = record()
    assert d["result"] == "hetblock_registration"
    assert d["lane"] == hb.LANE
    assert d["step"].startswith("3")


def test_it_could_not_be_written_until_scalehets_figures_had_come_back() -> None:
    d = record()
    tp.must_be_present(HARNESS_RECORD)
    assert d["harness_outcome"] == "HARNESS_REPRODUCES"
    reproduced = d["scalehet_figures_reproduce"]
    assert "REPRODUCE" in reproduced
    assert "bit-identically" in reproduced
    # the committed figures this lane re-tests, carried into the registration and not retyped
    assert "166.04792148316798" in reproduced
    assert "0.003996003996003996" in reproduced
    assert "0.013986013986013986" in reproduced


def test_every_registered_field_is_the_committed_modules_field_VERBATIM() -> None:
    d = record()
    assert len(hb.BOUND) == 57
    for key, name in hb.BOUND.items():
        assert key in d, f"BOUND field {name} is not in the registration"
        want = as_json(getattr(hb, name))
        assert d[key] == want, f"the registration's {key} is not the module's {name}"


def test_the_only_change_from_scalehet_is_the_unit_everything_else_is_its_object() -> None:
    d = record()
    assert d["min_genes"] == sh.MIN_GENES == mc.MIN_GENES == 100
    assert d["min_chromosomes"] == sh.MIN_CHROMOSOMES
    assert d["control_alpha"] == sh.CONTROL_ALPHA
    assert d["form_interval"] == list(sh.FORM_INTERVAL)
    assert d["tolerance"] == sh.TOLERANCE
    assert d["permutations"] == {"count": sh.PERMUTATIONS["count"], "seed": sh.PERMUTATIONS["seed"]}
    assert d["form"] == sh.FORM
    assert d["source"] == sh.SOURCE
    assert d["cells"] == list(sh.CELLS)
    assert d["scalehet_committed"]["admitted_count"] == 23


def test_alpha_is_the_LAXER_0_05_and_the_stricter_number_is_only_reported() -> None:
    d = record()
    assert d["control_alpha"] == 0.05
    assert d["stricter_alpha"] == 0.0167
    assert d["control_alpha"] > d["stricter_alpha"]
    assert "0.05" in d["alpha"]


def test_the_unit_is_the_1mb_locus_block_and_the_reason_the_element_unit_was_not_used_is_stated() -> None:
    d = record()
    assert d["governing_unit"] == "locus_block_1mb"
    assert d["block_span"] == hb.BLOCK_SPAN == 1_000_000
    assert d["robust_span"] == hb.ROBUST_SPAN == 5_000_000
    assert d["min_shared_elements"] == hb.MIN_SHARED_ELEMENTS
    assert d["element_sharing_graph"]["elements_driving_two_or_more_genes_of_this_table"] == 0
    why = d["block_unit_why_not_elements"]
    assert "HAS NO EDGES" in why
    assert "LOCUS PROXIMITY" in why
    assert "freshrna" in d["block_unit"]


def test_the_null_holds_every_chromosomes_gene_count_exactly_and_the_singleton_case_is_scalehets() -> None:
    d = record()
    assert "gene count" in d["null"].lower() or "OBSERVED cumulative" in d["null"]
    assert d["singleton_block_equivalence"]["bit_identical"] is True


def test_the_exposure_is_stated_twice_and_names_the_rejection_already_made() -> None:
    d = record()
    twice = " ".join(d["prior_exposure"]) + " " + d["prior_exposure_clause"]
    for reading in ("7d", "7h", "capped-sum", "M2h", "M2i", "scalesplit", "recovered", "scalehet"):
        assert reading in twice, f"the exposure does not name {reading}"
    assert "post hoc" in twice.lower()
    assert "LOOKED AT AND ACTED ON" in twice or "looked at and acted on" in twice.lower()


def test_BOTH_outcomes_are_written_out_here_before_either_is_seen() -> None:
    d = record()
    assert d["both_outcomes_are_written_here"]
    outcomes = d["outcomes"]
    assert set(outcomes) == {"SURVIVES_BLOCK_NULL", "WITHDRAWN_NOT_ESTABLISHED", "MIXED", "INCONCLUSIVE"}
    assert outcomes == hb.OUTCOMES


def test_the_withdrawal_outcome_is_written_plainly_and_is_not_softened() -> None:
    w = record()["outcomes"]["WITHDRAWN_NOT_ESTABLISHED"]
    assert "SCALEHET'S REJECTION IS WITHDRAWN" in w
    assert "NOT ESTABLISHED" in w
    assert "should no longer carry" in w
    assert "NOT SHOWN THAT ONE CONSTANT IS THE WHOLE CORRECTION EITHER" in w
    for softener in ("may suggest", "weak evidence", "borderline", "trend towards", "nearly"):
        assert softener not in w.lower()


def test_no_registered_outcome_proposes_a_regional_model() -> None:
    for name, text in record()["outcomes"].items():
        low = text.lower()
        assert "no regional model is proposed" in low or "untouched" in low, name
        assert "we propose" not in low
        assert "+0.2146" in text or "INCONCLUSIVE" in name or "2%" in text


def test_the_global_scale_and_its_variance_share_are_never_strengthened() -> None:
    d = record()
    text = json.dumps(d["outcomes"])
    assert "+0.2146" in text
    assert "about 2% of the variance" in text
    assert "explains" not in text.lower().replace("explained", "")


def test_NO_STATISTIC_OF_THIS_LANE_IS_IN_THE_REGISTRATION() -> None:
    d = record()
    assert d["nothing_of_this_lanes_statistic_is_in_this_file"]
    for key in ("q_block", "p_q_block", "block_null_q", "outcome", "reading"):
        assert key not in d, f"the registration must not carry {key}"


def test_every_input_the_run_will_read_is_pinned_by_sha256_and_declared_relatively() -> None:
    d = record()
    pinned = d["pinned_sha256"]
    assert pinned
    for path in pinned:
        assert not path.startswith("/"), f"input paths are RELATIVE (d72f959): {path}"
    for entry in d["result_manifest"]["inputs"]:
        assert not entry["path"].startswith("/")
        assert entry["sha256"]
    assert d["readers"], "every reader is declared"
    assert d["alphagenome_requests"] == 0
