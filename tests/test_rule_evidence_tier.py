# SPDX-License-Identifier: AGPL-3.0-or-later
"""The rule-evidence-tier PROPOSAL, tested against its six registered predictions.

Registered at `709aebc` (artefact) over the design committed at `0d4247b`, both before any of this
existed. Each prediction is tested here rather than asserted in the proposal, including the two that
could have found the proposal redundant.

NOTHING HERE WRITES A `.bio` FILE OR CHANGES A NUMBER IN ONE. The one test that needs a mutated
program mutates the PARSED IR in memory, and a separate test digests every `.bio` file in the demo
and organism corpus before and after the whole module runs and requires every digest to be
unchanged.
"""

from __future__ import annotations

import copy
import hashlib
import itertools
import json
import subprocess
import sys
from pathlib import Path

import pytest

from genomeos import manifest as mf
from genomeos.ir import Action, EvidenceKind
from genomeos.lang.parser import parse_file
from genomeos.provenance import rule_evidence_tier as ret
from genomeos.provenance import rule_number_sources as rns
from genomeos.runtime.grn import NetworkRuntime
from genomeos.runtime.uncertainty import UncertaintyReport, report_for_network

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import rule_evidence_tier_cost as cost_script  # noqa: E402

PROGRAM = ROOT / "data" / "demo" / "gastrulation.bio"
REGISTRATION = ROOT / "data" / "results" / "rule_evidence_tier_proposal_registration.json"


# --- 1. the axis is the census's, as the same objects ----------------------------------------------


def test_the_tiers_are_the_censuss_own_objects_not_a_second_copy():
    """Asserted with `is`, so the language cannot grow a divergent spelling of the census's classes."""
    assert ret.TIERS is rns.CASCADE
    assert ret.TIER_DEFINITIONS is rns.CLASSES
    assert ret.TIERED_FIELDS is rns.RULE_FIELDS
    assert ret.INTERACTION_STANDING is rns.ATTRIBUTION


def test_the_axis_has_five_values_and_only_four_are_tiers():
    assert len(ret.TIERS) == 4
    assert len(ret.AXIS_VALUES) == 5
    assert ret.NOT_ASSESSED in ret.AXIS_VALUES
    assert ret.NOT_ASSESSED not in ret.TIERS
    assert not ret.is_tier(ret.NOT_ASSESSED)
    assert all(ret.is_tier(t) for t in ret.TIERS)
    # and `not_assessed` is not a renaming of the census's own classes
    assert ret.NOT_ASSESSED not in ret.TIER_DEFINITIONS


def test_not_assessed_has_no_rank_so_it_cannot_be_ordered_against_a_tier():
    """The guard that keeps 'nobody looked' from being compared with 'looked and found nothing'."""
    with pytest.raises(ValueError, match="absence of a tier"):
        ret.rank(ret.NOT_ASSESSED)
    with pytest.raises(ValueError, match="not on the axis"):
        ret.rank("measured")
    assert [ret.rank(t) for t in ret.TIERS] == [0, 1, 2, 3]


# --- 2. P5: the floor, and what it refuses to claim ------------------------------------------------


def test_P5_an_empty_set_of_marks_floors_to_not_assessed_and_never_to_a_tier():
    """The vacuous-truth bug this guards: an unmarked result reading as measured."""
    assert ret.floor([]) == ret.NOT_ASSESSED


@pytest.mark.parametrize("n", [1, 2, 3])
def test_P5_any_unassessed_mark_makes_the_floor_undetermined(n):
    for combo in itertools.product(ret.AXIS_VALUES, repeat=n):
        got = ret.floor(combo)
        if ret.NOT_ASSESSED in combo:
            assert got == ret.NOT_ASSESSED, combo
        else:
            assert got == max(combo, key=ret.rank), combo
            assert ret.rank(got) == max(ret.rank(c) for c in combo), combo


def test_the_floor_refuses_a_value_off_the_axis():
    with pytest.raises(ValueError, match="not on the axis"):
        ret.floor(["M_measured_quoted", "probably_fine"])


# --- 3. P1: what a gate would cost ----------------------------------------------------------------


def _census_marks():
    return ret.census_slot_marks()


def test_P1_a_gate_at_measured_leaves_no_rule_of_the_program_integrable():
    """Registered prediction P1. Held: a gate is not a stricter corpus, it is an unrunnable one."""
    slots, _ = _census_marks()
    ids = sorted({s.rule for s in slots})
    assert len(ids) == 7
    assert ret.gate_survivors(ids, slots, "M_measured_quoted") == []
    assert ret.gate_survivors(ids, slots, "F_fitted_or_modelled") == []
    # and the gate only stops being fatal two tiers down, where it admits almost everything
    assert len(ret.gate_survivors(ids, slots, "E_existence_only")) == 6
    assert len(ret.gate_survivors(ids, slots, "U_unsourced")) == 7


def test_the_floor_over_the_census_marks_is_the_lowest_tier_it_assigned():
    slots, _ = _census_marks()
    ids = sorted({s.rule for s in slots})
    fl = ret.result_floor(ids, slots)
    assert fl["slots"] == 21
    assert fl["lowest_tier_it_rests_on"] == "U_unsourced"
    assert fl["determined"] is True
    assert fl["counts"][ret.NOT_ASSESSED] == 0


def test_a_result_over_the_program_as_committed_is_undetermined_not_unsourced():
    """The day-one state: the program carries no mark, so the floor is UNDETERMINED. Claiming
    `U_unsourced` here would assert a search of this program that the program does not record."""
    slots, _ = _census_marks()
    ids = sorted({s.rule for s in slots})
    fl = ret.result_floor(ids, [])
    assert fl["lowest_tier_it_rests_on"] == ret.NOT_ASSESSED
    assert fl["determined"] is False
    assert fl["counts"][ret.NOT_ASSESSED] == 21
    assert "never been assessed" in fl["statement"]


# --- 4. P2 and P3: whether the existing instrument already answers this -----------------------------


def test_P2_the_existing_uncertainty_report_does_not_already_answer_this():
    """Registered prediction P2, computed after it was registered and reported either way.

    HELD. The committed report labels this program's molecular level `medium` at 0.473, while the
    census classes none of its 21 rule numbers as measured. The label is not wrong on its own terms -
    it aggregates each DECLARATION's evidence kind and confidence, which are about the INTERACTION -
    but a reader cannot get the tier of the NUMBERS out of it, and `medium` is easy to mistake for it.
    """
    module = parse_file(PROGRAM)
    rep = report_for_network(module, module.rules).to_dict()["molecular"]
    assert rep["label"] == "medium"
    assert rep["label"] != "low"
    assert 0.47 <= rep["confidence"] <= 0.48
    slots, _ = _census_marks()
    assert ret.floor(s.tier for s in slots) == "U_unsourced"
    assert sum(1 for s in slots if s.tier == "M_measured_quoted") == 0


def test_P3a_the_existing_aggregation_is_a_mean_so_one_weak_item_is_diluted():
    """A mean cannot state a floor: six grounded items hide a seventh that rests on nothing."""
    rep = UncertaintyReport()
    for i in range(6):
        rep.add("molecular", 0.9, EvidenceKind.EXPERIMENTAL, f"strong{i}")
    rep.add("molecular", 0.0, EvidenceKind.NONE, "rests_on_nothing")
    lr = rep.levels["molecular"]
    assert lr.label == "high"  # with an item scoring 0.0 in the set
    assert ret.floor(["M_measured_quoted"] * 6 + ["U_unsourced"]) == "U_unsourced"


def test_P3b_the_existing_weakest_field_is_not_the_minimum():
    """Registered prediction P3(b). HELD, and it is a defect of an existing module, REPORTED AND NOT
    FIXED here: `genomeos/runtime/uncertainty.py` guards `weakest` with `score < lr.confidence`
    where `lr.confidence` is the RUNNING MEAN, so `weakest` ends up holding the last item that fell
    below a moving average rather than the lowest-scoring one. Its owner decides what to do."""
    scores = {"a": 1.0, "b": 0.5, "c": 0.6}
    rep = UncertaintyReport()
    for name, s in scores.items():
        rep.add("molecular", s, EvidenceKind.EXPERIMENTAL, name)
    reported = rep.levels["molecular"].weakest
    true_argmin = min(scores, key=lambda k: scores[k])
    assert true_argmin == "b"
    assert reported == "c"
    assert reported != true_argmin


# --- 5. P4: the tier cannot gate, as a property of the import closure -------------------------------


@pytest.mark.parametrize(
    "entry",
    ["genomeos/runtime/grn.py", "genomeos/runtime/gastrulation.py", "genomeos/runtime/uncertainty.py"],
)
def test_P4_no_runtime_module_can_even_read_the_tier(entry):
    """Stronger than a passing simulation: it holds for inputs nobody ran. If the tier is ever made
    to gate, this test is what fails."""
    closure = mf.counting_path(entry, ROOT)
    assert "genomeos/provenance/rule_evidence_tier.py" not in closure
    # and the one-way direction is the point: the tier module does not import the runtime either
    tier_path = "genomeos/provenance/rule_evidence_tier.py"
    assert "genomeos/runtime/grn.py" not in mf.counting_path(tier_path, ROOT)


# --- 6. P6: the coherence rule, validated against the only adjudication that exists ----------------


def test_P6_the_coherence_rule_holds_on_the_censuss_own_21_adjudications():
    """Registered prediction P6. HELD: a rule derived from the census's class definitions must hold
    on the census's own output, or the two-level design is unsound."""
    slots, rules = _census_marks()
    assert len(slots) == 21
    assert len(rules) == 7
    assert ret.coherence_problems(slots, rules) == []


def test_the_coherence_check_fires_on_a_planted_E_with_no_source_for_the_interaction():
    slots = [ret.SlotMark("R", "strength", "E_existence_only", "planted")]
    rules = [ret.RuleMark("R", "no_source_cited", "planted")]
    problems = ret.coherence_problems(slots, rules)
    assert [p["rule_broken"] for p in problems] == ["E_needs_a_source_for_the_interaction"]


def test_the_coherence_check_fires_on_a_planted_U_whose_citation_supports_the_rule():
    slots = [ret.SlotMark("R", "hill", "U_unsourced", "planted")]
    rules = [ret.RuleMark("R", "supports_declaration", "planted")]
    problems = ret.coherence_problems(slots, rules)
    assert [p["rule_broken"] for p in problems] == ["U_needs_the_citation_not_to_establish_the_interaction"]


def test_the_coherence_check_fires_on_a_bare_tier_that_names_no_adjudication():
    """A tier asserts a class someone decided, so a mark with no adjudication is a claim no artefact
    supports. `not_assessed` needs none, because it asserts nothing."""
    rules = [ret.RuleMark("R", "supports_declaration", "planted")]
    bare = [ret.SlotMark("R", "threshold", "E_existence_only", "")]
    assert [p["rule_broken"] for p in ret.coherence_problems(bare, rules)] == [
        "a_tier_names_its_adjudication"
    ]
    unmarked = [ret.SlotMark("R", "threshold", ret.NOT_ASSESSED, "")]
    assert ret.coherence_problems(unmarked, rules) == []


def test_a_tier_off_the_axis_is_caught_rather_than_counted():
    rules = [ret.RuleMark("R", "supports_declaration", "planted")]
    slots = [ret.SlotMark("R", "strength", "measured", "planted")]
    assert [p["rule_broken"] for p in ret.coherence_problems(slots, rules)] == ["off_the_axis"]


# --- 7. the population: what the dynamics actually rest on -----------------------------------------


def test_the_parser_synthesises_rules_the_program_does_not_write():
    """The finding that narrowed the population after registration, checked rather than asserted."""
    module = parse_file(PROGRAM)
    assert len(module.rules) == 11
    written = [r for r in module.rules if cost_script.SYNTHESISED_ID not in r.id]
    synthesised = [r for r in module.rules if cost_script.SYNTHESISED_ID in r.id]
    assert len(written) == 7
    assert len(synthesised) == 4
    # every synthesised rule carries all three numbers at the language default, stated in no program
    for r in synthesised:
        assert r.action is Action.PRODUCE
        assert (r.strength, r.threshold, r.hill) == (1.0, 1.0, 2.0)
    assert len(ret.tiered_rules(module.rules)) == 7


def test_the_synthesised_rules_three_numbers_are_inert_in_the_runtime():
    """Why those 12 numbers are outside the tier population, established from BEHAVIOUR and not from
    a reading of `grn.py`: perturbing them does not move the trajectory, and perturbing a number
    that IS in the population does. No `.bio` file is touched; the parsed IR is mutated in memory."""
    base = parse_file(PROGRAM)
    kw = dict(hours=6.0, dt=0.05, initial={"Sox2": 1.0}, record_every=10**9, clamp={"Nodal": 3.0})
    before = NetworkRuntime(base).run(**kw).final()

    inert = copy.deepcopy(base)
    touched = 0
    for r in inert.rules:
        if r.action is Action.PRODUCE:
            r.strength, r.threshold, r.hill = 0.01, 99.0, 11.0
            touched += 1
    assert touched == 4
    assert NetworkRuntime(inert).run(**kw).final() == before

    live = copy.deepcopy(base)
    for r in live.rules:
        if r.action is Action.INHIBIT:
            r.strength = 0.01
    assert NetworkRuntime(live).run(**kw).final() != before


# --- 8. nothing in this proposal writes a program ---------------------------------------------------


def _digests() -> dict[str, str]:
    out = {}
    for root in ("data/demo", "data/organisms/celegans", "genomeos/std"):
        for p in sorted((ROOT / root).rglob("*.bio")):
            out[str(p.relative_to(ROOT))] = hashlib.sha256(p.read_bytes()).hexdigest()
    return out


def test_no_bio_file_is_written_by_anything_in_this_proposal():
    before = _digests()
    assert before, "the corpus under test must not be empty or this test proves nothing"
    slots, rules = _census_marks()
    ids = sorted({s.rule for s in slots})
    ret.result_floor(ids, slots)
    ret.coherence_problems(slots, rules)
    ret.gate_survivors(ids, slots, "E_existence_only")
    ret.registration()
    cost_script.cost(scan_generated=False)
    assert _digests() == before


def test_the_proposal_has_no_code_path_that_writes_a_program():
    """Belt and braces: no `open(..., "w")`, no `write_text` and no `.bio` write anywhere in the
    proposal's own three files."""
    for rel in (
        "genomeos/provenance/rule_evidence_tier.py",
        "scripts/rule_evidence_tier_cost.py",
        "scripts/rule_evidence_tier_register.py",
    ):
        src = (ROOT / rel).read_text()
        assert "write_text" not in src, rel
        assert "open(" not in src or '"w"' not in src, rel
        assert ".unlink(" not in src, rel


# --- 9. the registration still says what it was registered as ---------------------------------------


def test_the_registration_artefact_matches_the_module_and_is_still_empty():
    payload = json.loads(REGISTRATION.read_text())
    assert payload["axis_values"] == list(ret.AXIS_VALUES)
    assert payload["tiers"] == list(ret.TIERS)
    assert payload["not_assessed"] == ret.NOT_ASSESSED
    assert payload["not_assessed_is_the_default_and_is_not_a_tier"] is True
    assert payload["tiered_fields"] == list(ret.TIERED_FIELDS)
    assert len(payload["predictions"]) == len(ret.PREDICTIONS) == 6
    # a registration may not acquire a finding afterwards
    assert payload["adoption_cost"] == {}
    assert payload["adoption_cost_entries"] == 0
    assert payload["prediction_outcomes"] == {}
    assert payload["prediction_outcomes_entries"] == 0
    assert payload["recommendation"] == ""
    assert payload["recommendation_length"] == 0
    assert payload["result_manifest"]["complete"] is True


def test_the_registration_does_not_restate_the_censuss_counts():
    """The census's figures are cited to 4d4003d, not copied into this proposal's artefact."""
    text = REGISTRATION.read_text()
    assert "4d4003d" in text
    assert "rule_number_sources_census.json" in text
    for forbidden in ("0 of 21", "17 have", "21 rule numbers of"):
        assert forbidden not in text, forbidden


# --- 10. the adoption cost, re-derived ---------------------------------------------------------------


#: Where the 24 git-ignored compiled chromosomes live. Named as a store rather than counted, so the
#: marker below skips BY NAME in a checkout that does not have them.
GENERATED_STORE = "data/knowledge/compiled"


def _tracked_bio() -> set[str]:
    """Every `.bio` path git tracks. A count over these is the same in every checkout."""
    out = subprocess.run(
        ["git", "-C", str(ROOT), "ls-files", "--", "*.bio"], capture_output=True, text=True, check=True
    )
    return {line for line in out.stdout.splitlines() if line}


def test_the_adoption_cost_of_the_hand_authored_corpus_is_what_is_reported():
    """The figures this proposal reports, re-derived from the corpus every run. The generated half is
    a 399 MB line scan and is not re-derived here; `--scan` does it.

    SPLIT 2026-10-03, and the split is the finding. `cost()` counts the FILESYSTEM, so
    `programs_total` was 67 on the laptop the figure was pinned on and 43 in CI, where
    `data/knowledge/compiled`'s 24 git-ignored chromosomes do not exist: `ef99f3c` was green here and
    red in both CI jobs from the moment it landed. Neither number was wrong - the ASSERTION was, for
    claiming a machine. So the tracked counts are asserted here, unconditionally and identically
    everywhere, and the local total is asserted in the test below under `needs_local_data`.

    Every hand-authored figure stays unconditional, and the reason is checked rather than assumed:
    all 33 hand-authored programs are TRACKED, which the first assertion measures. If one ever is
    not, that assertion fails and says so instead of letting a hand-authored figure drift with the
    machine - which is the whole class of defect this split is undoing.
    """
    c = cost_script.cost(scan_generated=False)
    tracked = _tracked_bio()
    hand = [f for f in c["files"] if f["kind"] == "hand_authored"]
    untracked_hand = sorted(f["path"] for f in hand if f["path"] not in tracked)
    assert untracked_hand == [], (
        "a hand-authored program git does not track makes every hand-authored figure below a fact "
        f"about this machine: {untracked_hand}"
    )
    tracked_files = [f for f in c["files"] if f["path"] in tracked]
    tracked_generated = [f for f in tracked_files if f["kind"] == "generated"]
    # the repository's own corpus: 43 = 33 hand-authored + 10 generated-and-committed
    assert len(tracked_files) == 43
    assert len(tracked_generated) == 10
    assert len(tracked_files) == len(hand) + len(tracked_generated)
    # the relationship is asserted and not only the totals, so a program added on either side cannot
    # leave these disagreeing again without saying which side it joined
    assert c["programs_total"] == c["programs_hand_authored"] + c["programs_generated"]
    assert c["programs_hand_authored"] == 33
    assert c["rules_hand_authored_written"] == 22
    assert c["rules_hand_authored_compiled"] == 41
    assert c["rules_hand_authored_synthesised_by_the_parser"] == 19
    assert c["rules_hand_authored_with_tiered_numbers"] == 20
    assert c["slots_hand_authored"] == 60
    assert c["programs_hand_authored_with_at_least_one_rule"] == 7
    assert c["default"] == ret.NOT_ASSESSED
    assert c["rules_generated_scanned"] is None


@pytest.mark.needs_local_data(
    GENERATED_STORE, how="compile the chromosomes: genomeos budget --chrom chrN --bio"
)
def test_the_local_generated_corpus_adds_the_twenty_four_compiled_chromosomes():
    """The other half of the split, kept rather than lost: 67 is the right total WHERE the store is.

    The 2026-10-03 figures stand unchanged in their own place. `programs_total` went 66 -> 67 and
    `programs_generated` 33 -> 34 when lane-headerfix added
    `data/organisms/human/noncoding_chr21_v2.bio` - v1 regenerated with the corrected cell-sentence
    header and v1's pinned bytes left intact - and not one hand-authored figure moved.
    """
    c = cost_script.cost(scan_generated=False)
    tracked = _tracked_bio()
    paths = {f["path"] for f in c["files"]}
    untracked = sorted(paths - tracked)
    assert len(untracked) == 24
    assert all(p.startswith(GENERATED_STORE + "/") for p in untracked), untracked
    assert all(f["kind"] == "generated" for f in c["files"] if f["path"] in set(untracked))
    assert c["programs_total"] == 67
    assert c["programs_generated"] == 34
    assert c["programs_total"] == len(paths & tracked) + len(untracked)


@pytest.mark.needs_local_data(
    GENERATED_STORE, how="compile the chromosomes: genomeos budget --chrom chrN --bio"
)
def test_the_hand_authored_figures_are_the_same_in_a_checkout_without_the_local_store(monkeypatch):
    """The CI case, PROVED here instead of inferred. The corpus enumeration is injected with the
    git-ignored store taken out - which is what a fresh checkout is - and every hand-authored figure
    above is unchanged, which is what makes those assertions safe to keep unconditional. The total
    comes out 43, the figure both CI jobs reported.

    Injection and not deletion: nothing is removed from this machine, and `cost()` is the real one.

    AND IT NEEDS THE STORE IT IS NAMED FOR, which is the joke the marker above is undoing. The
    injection subtracts the untracked programs from the real enumeration, so it only simulates a
    fresh checkout where untracked programs EXIST: its own first assertion says so -- "this machine
    holds no untracked program, so this proves nothing". In a genuine store-free checkout tracked ==
    total == 43, `len(fresh) < len(real)` is `43 < 43`, and the test FAILS in exactly the situation
    it is named for. CI run 37124753522 on 5c4589d reported that failure word for word. So it carries
    the same marker as the local-totals test above and SKIPS BY NAME there, saying "not run here"
    rather than claiming a red.

    Nothing it asserts is weakened, and this is not a substitute for the real check either: the
    hand-authored figures stay UNCONDITIONAL in
    `test_the_adoption_cost_of_the_hand_authored_corpus_is_what_is_reported`, which runs everywhere
    and measures that all 33 are tracked. This test adds the counterfactual on top of that, where the
    store is present to subtract.
    """
    tracked = _tracked_bio()
    real = cost_script.bio_files()
    fresh = [p for p in real if str(p.relative_to(cost_script.ROOT)) in tracked]
    assert len(fresh) < len(real), "this machine holds no untracked program, so this proves nothing"
    monkeypatch.setattr(cost_script, "bio_files", lambda: fresh)
    c = cost_script.cost(scan_generated=False)
    assert c["programs_total"] == 43
    assert c["programs_hand_authored"] == 33
    assert c["programs_generated"] == 10
    assert c["rules_hand_authored_written"] == 22
    assert c["rules_hand_authored_compiled"] == 41
    assert c["rules_hand_authored_synthesised_by_the_parser"] == 19
    assert c["rules_hand_authored_with_tiered_numbers"] == 20
    assert c["slots_hand_authored"] == 60
    assert c["programs_hand_authored_with_at_least_one_rule"] == 7
    assert c["scan_calibration"]["calibrated"] is True


def test_the_store_the_marker_names_is_git_ignored_so_the_marker_skips_and_never_raises():
    """`needs_local_data` RAISES on a path git tracks, by design, so the name it is given has to be
    the machine-local store and not a committed path. Asked of git, as that module asks it."""
    sys.path.insert(0, str(ROOT / "tests"))
    import local_data

    ignored, answered_about = local_data.ignored_by_name(GENERATED_STORE)
    assert ignored, f"git does not ignore {answered_about}, so the marker would raise rather than skip"
    # and NOT whether the store is on this machine: that is the very kind of assertion being removed
    # from this file today, and it would be red in CI for being true here


def test_the_line_scan_is_calibrated_against_the_parse_before_it_is_believed():
    """The scan is the only instrument that fits the generated corpus, so the figure it produces
    there is worth exactly what this calibration says."""
    c = cost_script.cost(scan_generated=False)
    assert c["scan_calibration"]["disagreements"] == []
    assert c["scan_calibration"]["calibrated"] is True
    assert c["scan_calibration"]["files_compared"] == 33
