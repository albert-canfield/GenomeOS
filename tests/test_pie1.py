# SPDX-License-Identifier: AGPL-3.0-or-later
"""Area E, the pie-1 discrepancy (scripts/pie1_registration.py): the closed list registered before any
rule changed, held to what was registered.

The registration is the REGISTRATION dict in committed code, and its digest is pinned here, so a
candidate, a source, a prediction or the falsifier cannot change without this test failing and a commit
saying why. Every registered edit is checked to apply exactly once to the programs at the base commit,
alone and in every registered pair, so no arm can be void for a reason the registration could have
caught.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from genomeos.lang import parse_file
from genomeos.organism.reference import ReferenceLineage
from genomeos.runtime.body import Body
from scripts.celegans_nolookup import counts, decided, founder_types, terminal_cells
from scripts.pie1_registration import (
    BASE_COMMIT,
    CANDIDATES,
    EXCLUDED,
    ORG,
    PAIRS,
    PREDICTIONS,
    REGISTRATION,
    SOURCES,
    apply_edits,
    base_files,
    registration_digest,
)
from tests.committed_data import committed

#: registered 2026-10-10 by lane-pie1, before any candidate was grown
REGISTRATION_DIGEST = "dca737065c4d9e9b59f35a2aed579401bed65c33bb3276db29d45ef2fd651052"
#: the three programs the candidates edit, as they were at the base commit (1c2ae67)
BASE_SHA256 = {
    "founders_contacts.bio": "35f52d00f904f2d2cdb923d373983c01bdd672160cab6640d88eafe7b0a58f83",
    "embryo_nolookup.bio": "72d2e4be67ea63456c53efa3e3a0180656d6e018983454b69cf56def500d4bac",
    "embryo_contacts.bio": "aee021970ebafa59c7c8372e3ca38c0a4987e24bc8da9cde6925914bf2134169",
}


@pytest.fixture(scope="module")
def files() -> dict[str, str]:
    return base_files()


def test_the_registration_is_the_one_committed_before_any_run():
    assert registration_digest() == REGISTRATION_DIGEST


def test_the_list_is_closed():
    assert set(CANDIDATES) == {"C1", "C2", "C3"}
    assert set(EXCLUDED) == {"X1", "X2", "X3", "X4", "X5"}
    assert PAIRS == (("C1", "C2"), ("C1", "C3"), ("C2", "C3"))
    assert set(PREDICTIONS) == {"baseline", "C1", "C2", "C3", "C1+C2", "C1+C3", "C2+C3"}
    assert REGISTRATION["base"]["commit"] == BASE_COMMIT


def test_every_candidate_cites_only_fetched_sources_and_every_exclusion_says_why():
    for cid, c in CANDIDATES.items():
        assert c["sources"] and set(c["sources"]) <= set(SOURCES), cid
        assert c["mechanism"] and c["edits"], cid
    for s in SOURCES.values():
        assert s["cite"] and s["fetched"] and s["quote"]
    assert all(e["what"] and e["why"] for e in EXCLUDED.values())


def test_the_base_programs_are_the_registered_ones(files):
    for name, digest in BASE_SHA256.items():
        assert hashlib.sha256(files[name].encode()).hexdigest() == digest, name


@pytest.mark.parametrize("arm", [*CANDIDATES, *("+".join(p) for p in PAIRS)])
def test_every_registered_edit_applies_exactly_once(files, arm):
    edited, void = apply_edits(files, arm)
    assert void == []
    assert edited != files


def test_an_unregistered_candidate_is_refused(files):
    with pytest.raises(SystemExit):
        apply_edits(files, "C4")


# ---- the run (01c64ba) and the one candidate it applied (2026-10-10) ----------------------------------

RESULT = "data/results/pie1_corrections.json"


@committed(RESULT)
def test_the_committed_result_is_the_registered_run_and_c3_alone_survives():
    d = json.loads(Path(RESULT).read_text())
    assert d["registration_digest"] == REGISTRATION_DIGEST and d["base_commit"] == BASE_COMMIT
    assert d["result_manifest"]["complete"] is True
    assert all(a["void"] == [] for a in d["arms"].values())
    v = d["verdicts"]
    assert v["survivors"] == ["C3"] and v["surviving_pairs"] == [] and v["apply"] == "C3"
    assert {k: a["survives"] for k, a in v["arms"].items()} == {
        "C1": False,
        "C2": False,
        "C3": True,
        "C1+C2": False,
        "C1+C3": False,
        "C2+C3": False,
    }
    # each rejection is the one registered in advance
    assert any("PAR-2" in b for b in v["arms"]["C1"]["breaks"])
    assert not any("MOM-2" in b for b in v["arms"]["C1"]["breaks"])
    assert any("MOM-2 terminal" in b for b in v["arms"]["C2"]["breaks"])
    assert any("MOM-5" in b for b in v["arms"]["C2"]["breaks"])
    for arm in d["arms"].values():  # the wild type is untouched in every arm
        w = arm["wild_type"]
        assert (w["decided"], w["right"], w["wrong"], w["seeds_all_56_56_0"]) == (56, 56, 0, True)
    pie1 = {k: _row(a["knockouts"]["PIE-1"]["terminal"]) for k, a in d["arms"].items()}
    assert pie1["baseline"] == (0, 20, 0, 20) and pie1["C3"] == (34, 20, 0, 54)
    held = v["arms"]["C3"]["pie1_readings"]
    assert held["R3 E is not MS-like"] and not held["R5 extra intestine (Intestine above the wild type's 34)"]
    # disclosed in 01c64ba: the baseline's Du et al. 2014 count on the no-lookup program is 10 / 1, not the
    # registered 11 / 0, because tissue_E retypes pop-1's E-like MS to Intestine before 300 min
    du = d["arms"]["baseline"]["du2014"]
    assert (du["embryo_contacts"]["reproduced"], du["embryo_contacts"]["missed"]) == (11, 0)
    assert du["embryo_nolookup"]["observed"]["pop-1"]["MS"] == "missed"


def _row(c: dict) -> tuple[int, int, int, int]:
    return (c["Intestine"], c["Muscle"], c["Reproductive"], c["decided"])


def test_the_programs_carry_the_registered_c3_edits_verbatim_and_nothing_else(files):
    edited, void = apply_edits(files, "C3")
    assert void == []
    for name in {f for f, _, _ in CANDIDATES["C3"]["edits"]}:
        assert (Path(ORG) / name).read_text() == edited[name], name


def test_the_applied_program_keeps_the_wild_type_and_moves_only_pie1():
    ref = ReferenceLineage.load()
    cells = terminal_cells(ref)
    module = parse_file(Path(ORG) / "embryo_nolookup.bio")
    targets = {d.id: d.to for d in module.decisions if d.action == "differentiate"}
    wt = decided(Body(module, means=True).run(until=800), cells, ref, targets)
    assert (wt["decided"], wt["right"], wt["wrong"]) == (56, 56, 0)
    pie1 = Body(module, means=True, knockouts={"PIE-1"}).run(until=800)
    f = founder_types(pie1, targets, ("P2", "ABp", "E"))
    assert f == {"P2": "EMSPrecursor", "ABp": "ABaPrecursor", "E": "EPrecursor"}
    assert _row(counts(decided(pie1, cells, ref, targets))) == (34, 20, 0, 54)
    # par-2 keeps Du et al. 2014's E-to-MS: P2 holds neither PIE-1 nor POS-1 (the registered premise)
    par2 = Body(module, means=True, knockouts={"PAR-2"}).run(until=800)
    assert founder_types(par2, targets, ("E",)) == {"E": "MSPrecursor"}
    mom2 = Body(module, means=True, knockouts={"MOM-2"}).run(until=800)
    assert founder_types(mom2, targets, ("E",)) == {"E": "MSPrecursor"}
