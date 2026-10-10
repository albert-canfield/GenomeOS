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

import pytest

from scripts.pie1_registration import (
    BASE_COMMIT,
    CANDIDATES,
    EXCLUDED,
    PAIRS,
    PREDICTIONS,
    REGISTRATION,
    SOURCES,
    apply_edits,
    base_files,
    registration_digest,
)

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
