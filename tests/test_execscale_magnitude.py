# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-execscale2: the calibrated magnitude as an opt-in runtime parameter, and what it must not move.

`tests/test_execscale.py` holds the census instrument that was registered BEFORE this parameter
existed (`38f1389`) and the census record itself (`81d3949`). This file holds the change the census
was taken for: `grn.MEASURED_MAGNITUDE_CORRECTION` as an opt-in parameter of the multiplicative
composition with default 1.0, and `compile.MAGNITUDE_CALIBRATED` as the compile option that emits a
programme declaring it.

WHAT IS ASSERTED, and each item is here because it is a way the change could go wrong:

1. At the DEFAULT the composition is the one that was there before, bit for bit, and no other rule is
   touched by the parameter existing. That is the whole safety property of the change.
2. The parameter is NOT INERT, and what it does is M2h's form: it scales the CENTRED log drive by
   exactly the correction, which is the quantity `ac5580b` measured (a_F * P_F(G, X), the gated run's
   log2 level in a cell minus its mean over the thirteen cells).
3. A correction other than 1.0 is REFUSED under a form M2h did not fit it on, because M2h fitted one
   scale per form and 0.2146 is the multiplicative form's alone.
4. Neither default moved: the runtime composes at 1.0 and the compiler emits nothing, so no committed
   programme and no default run moves.
5. The record's byte claim, read BLOB AGAINST RECORD at the commit the census was taken on, and the
   lane property that no programme's bytes differ from that record except by a commit after it.

WHY ITS OWN FILE. Adding these to `tests/test_execscale.py` was tried and abandoned: the honest
version of (5) replaces a working-tree comparison that a peer's programme commit had already turned
red, and replacing it DELETES lines an earlier commit of this lane put there, which this checkout's
rules do not allow the lane that is adding to do. So nothing committed is removed and the additions
live here. The inherited red is reported to the coordinator, not papered over, and (5) is the
property this lane is actually answerable for.
"""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path

import pytest

from genomeos.attribution import compile as compile_mod
from genomeos.lang import parse
from genomeos.runtime import grn
from tests.tracked_paths import must_be_present

ROOT = Path(__file__).resolve().parents[1]
RECORD = "data/results/execscale_census.json"
#: the commit the census record was taken on, named so its byte claim is a claim about a tree that
#: exists for ever rather than about whatever the working tree holds today
CENSUS_COMMIT = "81d3949"

TINY_MAG = (
    "module tiny.mag\n"
    "gene g { max: 10; basal: 0.0 }\n"
    "protein S {}\n"
    "rule S activates g { strength: 1.0; threshold: 1.0; hill: 2 }\n"
)


def _record() -> dict:
    must_be_present(RECORD)
    return json.loads((ROOT / RECORD).read_text())


def _blob(rev: str, path: str) -> bytes:
    """One file's bytes at one commit, read from git and not from the working tree."""
    return subprocess.run(["git", "show", f"{rev}:{path}"], cwd=ROOT, capture_output=True, check=True).stdout


def _tiny_module(extra: str = ""):
    m = parse(TINY_MAG + extra)
    m.entities["g"].attrs[grn.MULTIPLICATIVE_GAIN] = 1.0
    return m


# ---- the constant and the defaults ---------------------------------------------------------


def test_the_constant_is_M2hs_fixed_number_and_the_DEFAULT_IS_ONE():
    assert grn.MEASURED_MAGNITUDE_CORRECTION == 0.2146
    assert grn.DEFAULT_MAGNITUDE_CORRECTION == 1.0
    assert grn.MAGNITUDE_CORRECTION_PARAM == "magnitude_correction"
    src = (ROOT / "genomeos" / "runtime" / "grn.py").read_text()
    # the name says it is a measured correction, M2h's sha is in the docstring, and the global average
    # is called one, because a reader of this parameter must not read it as a tuning knob
    assert "ac5580b" in src and "not a tuning knob" in src
    assert "GLOBAL AVERAGE" in src and "d9e3ddd" in src
    # the engine file keeps its Apache-2.0 header (tests/test_engine_licence_headers.py owns the rule)
    assert src.startswith("# SPDX-License-Identifier: Apache-2.0\n")


def test_at_the_DEFAULT_the_composition_is_the_one_that_was_there_before():
    """Default 1.0 is the whole safety property: the drive is exp(gain * W), unmultiplied."""
    terms = [0.3, 0.7, 0.2]
    want = math.exp(1.3 * sum(terms))
    assert grn.combine_activators(terms, "multiplicative", gain=1.3) == want
    assert grn.combine_activators(terms, "multiplicative", gain=1.3, magnitude_correction=1.0) == want
    # and every other rule is untouched by the parameter existing
    for rule in grn.ACTIVATOR_RULES:
        assert grn.combine_activators(terms, rule) == grn.combine_activators(
            terms, rule, magnitude_correction=grn.DEFAULT_MAGNITUDE_CORRECTION
        )
    # a gene with no activator takes A = 1 under every magnitude, as it always did
    assert grn.combine_activators([], "multiplicative", gain=1.3, magnitude_correction=0.2146) == 1.0


def test_PLANT_the_parameter_is_NOT_INERT_and_scales_the_centred_log_drive_by_exactly_c():
    """The correction must do something, and the something must be M2h's form.

    M2h scaled the CENTRED prediction: the run's log2 level in a cell minus its mean over the cells.
    Under this rule log A = c * gain * W, so the difference in log drive between two cells is scaled
    by exactly c -- which is what is checked here, because that is the quantity the measurement is
    about. The base of the log does not enter: scaling commutes with a change of base.
    """
    c = grn.MEASURED_MAGNITUDE_CORRECTION
    hi, lo = [0.8, 0.6], [0.2, 0.1]
    a_hi = grn.combine_activators(hi, "multiplicative", gain=1.0)
    a_lo = grn.combine_activators(lo, "multiplicative", gain=1.0)
    c_hi = grn.combine_activators(hi, "multiplicative", gain=1.0, magnitude_correction=c)
    c_lo = grn.combine_activators(lo, "multiplicative", gain=1.0, magnitude_correction=c)
    assert c_hi != a_hi and c_lo != a_lo  # not inert
    assert math.log(c_hi) - math.log(c_lo) == pytest.approx(c * (math.log(a_hi) - math.log(a_lo)))
    # and it is the same number whichever way it is applied, because the exponent is linear in it
    assert c_hi == pytest.approx(grn.combine_activators(hi, "multiplicative", gain=c))


def test_the_correction_is_REFUSED_on_a_form_it_was_not_fitted_on_and_on_a_bad_number():
    for rule in grn.ACTIVATOR_RULES + ("competition",):
        with pytest.raises(ValueError, match="multiplicative"):
            grn.combine_activators([0.5], rule, magnitude_correction=0.2146)
    for bad in (0.0, -0.5, float("nan"), float("inf")):
        with pytest.raises(ValueError, match="magnitude_correction"):
            grn.combine_activators([0.5], "multiplicative", magnitude_correction=bad)


# ---- the runtime ---------------------------------------------------------------------------


def test_the_runtime_takes_the_correction_and_a_RUN_at_the_default_is_the_run_that_always_was():
    m = _tiny_module()
    assert grn.NetworkRuntime(m).magnitude_correction == 1.0
    plain = grn.NetworkRuntime(m, activator_combination="multiplicative")
    assert plain.magnitude_correction == 1.0
    asked = grn.NetworkRuntime(
        m,
        activator_combination="multiplicative",
        magnitude_correction=grn.MEASURED_MAGNITUDE_CORRECTION,
    )
    assert asked.magnitude_correction == grn.MEASURED_MAGNITUDE_CORRECTION
    held = {"S": 2.0}
    a = plain.run(2.0, dt=0.05, clamp=held).final()["g.mRNA"]
    b = asked.run(2.0, dt=0.05, clamp=held).final()["g.mRNA"]
    explicit = (
        grn.NetworkRuntime(m, activator_combination="multiplicative", magnitude_correction=1.0)
        .run(2.0, dt=0.05, clamp=held)
        .final()["g.mRNA"]
    )
    assert a == explicit  # asking for the default is the same run, bit for bit
    assert b != a and b < a  # and asking for the measured magnitude is a different, smaller drive


def test_a_PROGRAMME_may_declare_the_correction_and_one_declaring_none_runs_at_the_default():
    declared = _tiny_module(
        f"param magnitude_correction = {grn.MEASURED_MAGNITUDE_CORRECTION} "
        '{ evidence: inferred "M2h ac5580b"; confidence: 0.5 }\n'
    )
    assert (
        grn.NetworkRuntime(declared, activator_combination="multiplicative").magnitude_correction
        == grn.MEASURED_MAGNITUDE_CORRECTION
    )
    # the caller's argument wins over the programme's declaration, and no declaration means 1.0
    assert (
        grn.NetworkRuntime(
            declared, activator_combination="multiplicative", magnitude_correction=1.0
        ).magnitude_correction
        == 1.0
    )
    assert grn.NetworkRuntime(_tiny_module()).magnitude_correction == 1.0
    # a declared correction under a form it was not fitted on is refused, not silently applied
    with pytest.raises(ValueError, match="MULTIPLICATIVE"):
        grn.NetworkRuntime(declared)


def test_no_COMMITTED_programme_declares_the_correction_so_no_committed_run_moves():
    for row in _record()["programs"]:
        must_be_present(row["path"])
        assert grn.MAGNITUDE_CORRECTION_PARAM not in (ROOT / row["path"]).read_text()


# ---- the compile option --------------------------------------------------------------------


def test_the_compile_OPTION_is_opt_in_and_its_DEFAULT_IS_UNCHANGED():
    assert compile_mod.MAGNITUDE_DEFAULT == compile_mod.MAGNITUDE_COMPILED == "compiled"
    assert compile_mod.MAGNITUDES == ("compiled", "calibrated")
    assert compile_mod.magnitude_lines() == []
    assert compile_mod.magnitude_lines(compile_mod.MAGNITUDE_COMPILED) == []
    with pytest.raises(ValueError, match="no magnitude"):
        compile_mod.magnitude_lines("something_else")
    # and the emission default this lane was told not to move is where it was
    assert compile_mod.EMISSION_DEFAULT == compile_mod.ONE_TARGET


def test_the_compile_OPTION_emits_a_param_the_runtime_READS_and_describes_it_honestly():
    lines = compile_mod.magnitude_lines(compile_mod.MAGNITUDE_CALIBRATED)
    body = "\n".join(lines)
    assert f"param magnitude_correction = {grn.MEASURED_MAGNITUDE_CORRECTION}" in body
    assert "ac5580b" in body and "d9e3ddd" in body and "GLOBAL AVERAGE" in body
    assert "measured correction and not a tuning" in body
    # what it emits is a programme the runtime reads back at the measured magnitude
    rt = grn.NetworkRuntime(parse(TINY_MAG + body + "\n"), activator_combination="multiplicative")
    assert rt.magnitude_correction == grn.MEASURED_MAGNITUDE_CORRECTION
    # and it adds no rule and no entity, so no `# test:` count of any programme can move
    assert not [ln for ln in lines if ln.startswith(("rule ", "gene ", "element ", "region "))]
    assert len(parse(TINY_MAG + body + "\n").rules) == len(parse(TINY_MAG).rules)


# ---- what this lane must not have moved ----------------------------------------------------


def test_the_record_is_byte_identical_to_the_TREE_IT_WAS_TAKEN_ON():
    """The claim's own check, blob against record: every programme's digest and size at `81d3949`.

    BLOB TO RECORD, NOT BLOB TO WORKING TREE, for the reason
    `tests/test_licence_headers_are_insertions.py` gives for its own sides: a working-tree side
    "would have turned a one-off proof into a standing ban on ever editing" those files. The census
    is a record of ONE tree, so that is the tree it is checked against, and it stays true for ever.
    """
    d = _record()
    assert len(d["programs"]) == 58
    for row in d["programs"]:
        raw = _blob(CENSUS_COMMIT, row["path"])
        assert (hashlib.sha256(raw).hexdigest(), len(raw)) == (row["sha256"], row["bytes"]), row["path"]


def test_no_programme_MOVED_except_by_a_COMMIT_after_the_census():
    """What this lane's change must not do, as a property that survives the peers' own commits.

    A working-tree comparison against the record cannot say this by itself, because other lanes
    commit programmes of their own: `data/organisms/human/haematopoiesis.bio` moved 25 bytes at
    19:20 on 2026-10-12 in `71f9d47`, which restored `fraction:` where `share:` had been, after the
    census was taken at 16:52. So the property asserted is the one this lane is answerable for: every
    programme's bytes are its own HEAD blob's bytes -- no uncommitted edit of anyone's has moved one
    -- and a programme differing from the record is one that a commit after the census changed. A
    byte this lane moved would be an uncommitted difference with no commit behind it and fails here,
    in whichever file it is, with none exempted by name.
    """
    d = _record()
    moved = []
    for row in d["programs"]:
        path = row["path"]
        must_be_present(path)
        now = (ROOT / path).read_bytes()
        assert hashlib.sha256(now).hexdigest() == hashlib.sha256(_blob("HEAD", path)).hexdigest(), path
        if hashlib.sha256(now).hexdigest() != row["sha256"]:
            log = subprocess.run(
                ["git", "log", "--oneline", f"{CENSUS_COMMIT}..HEAD", "--", path],
                cwd=ROOT,
                capture_output=True,
                text=True,
                check=True,
            ).stdout.strip()
            assert log, f"{path} differs from the census and no commit after {CENSUS_COMMIT} touched it"
            moved.append(path)
    # a moved byte must not have moved a COUNT: the census is held to counts as well as to bytes
    for path in moved:
        if path in d["compiled_rule_totals"]:
            assert len(parse((ROOT / path).read_text()).rules) == d["compiled_rule_totals"][path], path


# ---- the chr21 demo record -----------------------------------------------------------------
# These read data/results/execscale_demo.json, which is committed in the same commit as they are.

DEMO = "data/results/execscale_demo.json"


def _demo() -> dict:
    must_be_present(DEMO)
    return json.loads((ROOT / DEMO).read_text())


def test_the_demo_record_says_it_claims_NOTHING_ABOUT_RNA_and_names_the_average():
    d = _demo()
    assert d["result_manifest"]["complete"] is True
    assert "NOTHING ABOUT RNA" in d["claim"] or "NOTHING about RNA" in d["does_not_claim"]
    assert "nothing is fitted" in d["does_not_claim"]
    assert d["magnitude"]["value"] == grn.MEASURED_MAGNITUDE_CORRECTION
    assert d["magnitude"]["source_commit"] == "ac5580b"
    assert "global AVERAGE" in d["magnitude"]["it_is_a_global_average"]
    assert "d9e3ddd" in d["magnitude"]["it_is_a_global_average"]
    assert d["magnitude"]["default_in_the_runtime"] == 1.0
    # the direction of every change is arithmetic under this rule, and the record says so itself
    assert "not a finding" in " ".join(d["difference"])


def test_the_demo_programme_differs_from_the_compiled_one_by_the_BLOCK_AND_NOTHING_ELSE():
    p = _demo()["programmes"]
    assert p["chrom"] == "chr21"
    assert p["lines_added"] == p["magnitude_block_lines"] == 11
    assert p["added_lines_are_the_magnitude_block"] is True
    assert p["declared_tests_identical"] is True
    # the counts the census is held to are the counts the calibrated programme still declares
    assert p["declared_tests"] == [
        "# test: entities >= 6172",
        "# test: unknowns == 20",
        "# test: rules == 5176",
    ]
    assert p["compiled_sha256"] != p["calibrated_sha256"]


def test_the_demo_RAN_at_the_measured_magnitude_read_off_the_programme_and_was_not_vacuous():
    d = _demo()
    plain, cal = d["run_at_the_compiled_magnitude"], d["run_at_the_calibrated_magnitude"]
    assert plain["magnitude_correction"] == 1.0
    assert cal["magnitude_correction"] == grn.MEASURED_MAGNITUDE_CORRECTION
    assert cal["read_from"] == "the programme's own param"
    # a run with no active rule composes nothing whatever the magnitude is, so the count is the
    # evidence that this run was not that run
    assert plain["active_rules"] == cal["active_rules"] == 372
    assert cal["genes_with_at_least_one_activator"] == 51
    assert cal["context"] == {"cell_type": "K562"}
    # and the two runs are the same run but for the magnitude
    for key in ("genes", "active_rules", "clamped_sources", "species", "mrna_species"):
        assert plain[key] == cal[key], key


def test_the_demo_difference_is_the_51_genes_that_have_an_activator_in_that_cell():
    f = _demo()["difference"]
    assert f["species_compared"] == 193
    assert f["species_whose_level_moved"] == 51
    assert f["species_unmoved"] == 142
    assert f["calibrated_is_lower_for"] == 51
    assert f["calibrated_is_higher_for"] == 0
    r = f["moved_level_ratio_calibrated_over_compiled"]
    assert r["n"] == 51
    assert 0.16 < r["min"] < 0.17
    assert 0.76 < r["median"] < 0.77
