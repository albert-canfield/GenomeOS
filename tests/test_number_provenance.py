# SPDX-License-Identifier: AGPL-3.0-or-later
"""The number-provenance census: the classification partitions, and its premises still hold."""

from __future__ import annotations

import dataclasses
from pathlib import Path

import pytest

from genomeos.ir import model as ir_model
from genomeos.lang import number_provenance as np

GASTRULATION = "data/demo/gastrulation.bio"


# --- the defaults transcribed into SCOPE must still be the language's --------------------------


def _default(cls: type, name: str) -> object:
    for f in dataclasses.fields(cls):
        if f.name == name:
            return f.default
    raise AssertionError(f"{cls.__name__} has no field {name}")


@pytest.mark.parametrize(
    ("kind", "field", "cls", "attr"),
    [
        ("rule", "strength", ir_model.Rule, "strength"),
        ("rule", "threshold", ir_model.Rule, "threshold"),
        ("rule", "hill", ir_model.Rule, "hill"),
        ("gene", "basal", ir_model.Gene, "basal_rate"),
        ("protein", "initial", ir_model.Protein, "initial"),
        ("field", "diffusion", ir_model.Field, "diffusion"),
        ("field", "decay", ir_model.Field, "decay"),
        ("signal", "threshold", ir_model.Signal, "threshold"),
        ("event", "rate", ir_model.Event, "rate"),
        ("regime", "threshold", ir_model.Regime, "threshold"),
        ("stage", "from", ir_model.Stage, "start"),
        ("timer", "lengthening", ir_model.Timer, "lengthening"),
        ("decision", "priority", ir_model.Decision, "priority"),
    ],
)
def test_scope_transcribes_the_language_default(kind: str, field: str, cls: type, attr: str) -> None:
    """A census that says `omitted; the language supplies 2.0` must be reading the real default."""
    slot = next(s for s in np.SCOPE if s.kind == kind and s.field == field)
    assert slot.language_default == pytest.approx(float(_default(cls, attr)))


def test_the_rule_defaults_are_the_ones_the_census_names() -> None:
    r = ir_model.Rule(id="r", source="A", action="activates", target="B")
    assert (r.strength, r.threshold, r.hill) == (1.0, 1.0, 2.0)


# --- the cascade partitions ----------------------------------------------------------------------


def test_every_slot_of_every_program_gets_exactly_one_known_class() -> None:
    for p in np.programs():
        rows = np.census_file(p)
        for r in rows:
            assert r["cls"] in np.CLASSES, (p, r)


def test_written_and_omitted_do_not_overlap() -> None:
    omitted = {"D_language_default", "D_no_number"}
    for p in np.programs():
        for r in np.census_file(p):
            assert (r["cls"] in omitted) == (not r["written"]), (p, r)


def test_the_denominator_is_written_plus_defaults() -> None:
    """Every scoped slot is either in the denominator or in one of the two outside classes."""
    for p in np.programs():
        rows = np.census_file(p)
        inside = [r for r in rows if r["cls"] not in ("D_no_number", "N_non_numeric")]
        outside = [r for r in rows if r["cls"] in ("D_no_number", "N_non_numeric")]
        assert len(inside) + len(outside) == len(rows)
        for r in inside:
            assert r["written"] or r["cls"] == "D_language_default", r


# --- the citation-furniture strip ----------------------------------------------------------------


def test_a_figure_label_is_not_a_quoted_number() -> None:
    """`hill: 3` against `Fig. 3C` must not be read as sourced. This is the lane's whole point."""
    rows = np.census_file(GASTRULATION)
    hill = next(
        r
        for r in rows
        if r["kind"] == "rule" and r["declaration"] == "Sox17 inhibits TBXT" and r["field"] == "hill"
    )
    assert "Fig. 3C" in hill["evidence_text"]
    assert hill["value"] == 3.0
    assert hill["cls"] == "E_existence_only"
    assert 3.0 not in hill["numbers_quoted_after_strip"]


def test_the_strip_removes_years_volumes_pages_and_accessions() -> None:
    s = np.strip_furniture("Lolas et al. 2014, PNAS 111:4478, Fig. 3C; UniProt Q9H6I2; GENCODE v50")
    assert np.numbers_in(s) == []


def test_the_strip_keeps_a_quantity() -> None:
    s = np.strip_furniture("Collier et al. 1996: f(d) = d^2 / (0.01 + d^2)")
    assert 0.01 in np.numbers_in(s)
    assert 1996 not in np.numbers_in(s)


# --- the premises this census was built on -------------------------------------------------------


def test_nothing_in_gastrulation_represses_sox17() -> None:
    """`d340e5c`: SOX17's repression factor is exactly 1.00000 because no rule inhibits it.

    If a lane adds one, this test fails and the census's framing has to be rewritten rather than
    quietly carried forward.
    """
    decls = np.read_program(GASTRULATION)
    inhibitors = [d.name for d in decls if d.kind == "rule" and d.name.endswith("inhibits SOX17")]
    assert inhibitors == []


def test_the_gastrulation_header_claims_mutual_repression_the_file_does_not_contain() -> None:
    """Reported, never edited: data/demo is shared. The claim is false and the test pins it."""
    text = Path(GASTRULATION).read_text()
    assert "Mutual repression makes the choice sharp" in text
    decls = np.read_program(GASTRULATION)
    rules = [d.name for d in decls if d.kind == "rule"]
    assert "Sox17 inhibits TBXT" in rules and "Sox17 inhibits SOX2" in rules
    assert not any(r.endswith("inhibits SOX17") for r in rules)


def test_the_inhibiting_convention_is_a_repeated_value_not_a_language_default() -> None:
    """`53d9936`: four inhibiting rules write the same `threshold: 2.0; hill: 3`.

    The language's defaults are 1.0 and 2.0, so these numbers are the module's own convention. The
    census must tag them `repeated_in_file` and not `equals_language_default`.
    """
    rows = np.census_file(GASTRULATION)
    thr = [r for r in rows if r["kind"] == "rule" and r["field"] == "threshold" and r["value"] == 2.0]
    hill = [r for r in rows if r["kind"] == "rule" and r["field"] == "hill" and r["value"] == 3.0]
    # the brief this lane was given said "the module\'s own `threshold 2.0; hill 3`" as though the
    # pair covered the four inhibiting rules. It does not: `Sox2 inhibits TBXT` writes threshold
    # 4.0, so threshold 2.0 is on three rules, and hill 3 is on five, the activating
    # `Tbxt activates SOX17` included. The convention is wider than the inhibiting rules.
    assert len(thr) == 3 and len(hill) == 5
    for r in thr + hill:
        assert "repeated_in_file" in r["convention_tags"]
        assert "equals_language_default" not in r["convention_tags"]


def test_every_gastrulation_strength_is_one_and_that_is_the_language_default() -> None:
    rows = np.census_file(GASTRULATION)
    strengths = [r for r in rows if r["kind"] == "rule" and r["field"] == "strength"]
    assert len(strengths) == 7
    ones = [r for r in strengths if r["value"] == 1.0]
    assert len(ones) == 6, "six rules write strength 1.0; the seventh writes 0.6"
    for r in ones:
        assert "equals_language_default" in r["convention_tags"]


# --- the census writes nothing ------------------------------------------------------------------


def test_the_census_does_not_write_to_any_program() -> None:
    before = {p: p.read_bytes() for p in np.programs()}
    for p in np.programs():
        np.census_file(p)
    for p, b in before.items():
        assert p.read_bytes() == b, f"{p} was modified"


def test_every_derivation_is_a_candidate_the_matcher_exposed() -> None:
    """An S_derived promotion may only reach a slot that classified E and quotes some number."""
    candidates = set()
    for p in np.programs():
        for r in np.census_file(p):
            if r.get("derivation_candidate") or r["cls"] == "S_derived":
                candidates.add((r["file"], r["line"], r["field"]))
    for key in np.DERIVATIONS:
        assert key in candidates, f"{key} was promoted but the matcher never exposed it"


def test_the_registration_states_the_population_the_code_reads() -> None:
    """The registration's population must be the one `programs()` enumerates.

    This asserts the MECHANISM and not a file count. An earlier version of this test asserted
    `len(np.programs()) == 19`, the figure in this repository, and that is a property of this
    repository's data rather than of the code: `scripts/package_engine.py` ships `genomeos/lang`,
    where this module lives, together with only six of the nineteen programs, so the packaged
    engine's suite qualified the test and could not satisfy it - 6 against 19, red for every lane.
    The population figure is a registered datum and belongs in the result, where
    `scripts/number_provenance_register.py` already digests every program by name; what belongs
    here is that the code reads exactly the declared directories, whatever a tree holds of them.
    """
    reg = np.registration()
    assert reg["population"]["directories"] == list(np.POPULATION_DIRS)
    assert len(reg["scoped_fields"]) == len(np.SCOPE)
    assert set(reg["classes"]) == set(np.CLASSES)
    found = np.programs()
    assert found, "the population is empty: no .bio program under any declared directory"
    for p in found:
        assert p.suffix == ".bio", p
        assert p.parent.as_posix() in np.POPULATION_DIRS, p
    assert Path(GASTRULATION) in found, "the module this census is about must be in the population"
    # nothing is filtered out behind the declared directories: the population is every .bio in them
    independently = {q for d in np.POPULATION_DIRS for q in Path(d).glob("*.bio")}
    assert set(found) == independently


def test_the_corpus_is_nineteen_programs_in_the_source_repository() -> None:
    """The registered population figure, checked only where the whole corpus can be present.

    The condition is `.git`, not a program's path. A witness path would defeat itself:
    `scripts/package_engine.py` scans the tests it selects for path literals under `data/` and
    copies what it finds, so naming a program the engine does not ship would make the engine ship
    it. `.git` is under neither `data/` nor `docs/`, so no scan can satisfy it, and it is present in
    the checkout and in the verdict worktree - the two places the figure is a real claim.

    It skips by name rather than passing quietly, so a reader can tell `not run here` from `passed`.
    `needs_local_data` is the marker for this shape of thing, but it raises on a path git does not
    ignore and every program here is committed, so the condition is stated directly instead.
    """
    found = np.programs()
    if not Path(".git").exists():
        pytest.skip(
            f"not the source repository (no .git): this tree holds {len(found)} programs. The "
            "packaged engine ships genomeos/lang, where this module lives, with only part of the "
            "corpus, so the repository's population figure is not a claim about it. The figure is "
            "pinned in data/results/number_provenance_registration.json, which digests every "
            "program by name."
        )
    assert len(found) == 19
