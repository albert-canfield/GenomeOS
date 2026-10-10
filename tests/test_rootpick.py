# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-rootpick: the root rule, its ground, and the one thing that would make it worthless.

The guard that matters most is `test_the_ground_is_quoted_from_the_bridge_verbatim`. The entire rule
rests on what the bridge says the strength convention does not license, so a paraphrase that drifted
would leave a rule with no ground while still reading as though it had one.
"""

from __future__ import annotations

import ast
import json
from pathlib import Path

import pytest

from genomeos.attribution import rootpick as rp

ROOT = Path(__file__).resolve().parent.parent
CENSUS = ROOT / "data/results/capsumrna_census.json"


def _unwrapped(text: str) -> str:
    """bridge.py's statements live in wrapped `#:` comments, so a sentence spans lines.

    Normalising the comment marker and the line breaks tests the real claim -- that the sentence is in
    that file -- without the check being defeated by where the author happened to wrap. It does NOT
    normalise words: a changed word still fails, which is the point. The first version of this test
    compared raw bytes, failed on exactly this wrapping, and caught my own quote drifting from the
    source it claims to quote; that is why the normalisation is named here rather than done quietly.
    """
    lines = [ln.lstrip().removeprefix("#:").removeprefix("#").strip() for ln in text.splitlines()]
    return " ".join(" ".join(lines).split())


def test_the_ground_is_quoted_from_the_bridge_verbatim() -> None:
    """Each fragment of GROUND must appear in bridge.py's own text, modulo comment wrapping."""
    bridge = _unwrapped((ROOT / "genomeos/attribution/bridge.py").read_text())
    for fragment in (
        "V and s are NOT separately identified: only their product is",
        "DECLARED_STRENGTH fixes the split by convention and by nothing else",
        "it is not a measurement that the element saturates its",
        "not the ceiling, not the strength",
    ):
        assert fragment in bridge, f"GROUND quotes {fragment!r}, which is not in bridge.py"
        assert fragment in _unwrapped(rp.GROUND), f"{fragment!r} left GROUND"


def test_the_capped_branch_is_the_one_the_bridge_names_as_unlicensed() -> None:
    """The rule's logic in one assertion: the capped root needs the ceiling, and that is unlicensed."""
    from genomeos.attribution import bridge

    assert "not the ceiling" in bridge.DECLARED_STRENGTH_LICENSES
    assert rp.licensed_root("fit_not_unique") == "uncapped"
    assert rp.licensed_root("capped_root_only") == "unknown"


def test_every_census_status_has_a_verdict_and_an_unknown_one_raises() -> None:
    for status in rp.VERDICT_BY_STATUS:
        assert rp.licensed_root(status) in ("uncapped", "unknown", "no_fit")
    with pytest.raises(ValueError, match="unknown capped-sum status"):
        rp.licensed_root("a_status_nobody_registered")


def test_no_status_differs_from_the_mean_under_this_rule() -> None:
    """The finding, as arithmetic from the rule rather than a measurement."""
    assert [s for s in rp.VERDICT_BY_STATUS if rp.differs_from_the_mean(s)] == []


def test_only_sum_capped_is_ambiguous_and_the_others_say_why_not() -> None:
    assert "AMBIGUOUS" in rp.FORMS["sum_capped"]
    for form in ("mean", "competition", "multiplicative"):
        assert "single-valued" in rp.FORMS[form], form


def test_the_rule_is_applied_to_the_committed_census_and_admits_what_it_says() -> None:
    counts = json.loads(CENSUS.read_text())["reader_selection"]["capped_sum_status"]
    out = rp.reclassify(counts)
    assert out["admitted"] == counts["runnable"] + counts["fit_not_unique"]
    assert out["unknown_under_this_rule"] == counts["capped_root_only"]
    assert out["differs_from_the_mean"] == 0


def test_PLANTED_a_rule_that_took_the_capped_root_would_cross_the_gate() -> None:
    """Why FORECLOSES is in the registration: the rejected choice is the one that buys a result.

    Not a counterfactual anyone may adopt -- it is the disclosure that the chosen rule is the costly
    one. If this ever fails, the asymmetry the registration rests on has gone and the disclosure is
    stale.
    """
    census = json.loads(CENSUS.read_text())["reader_selection"]
    would_be = (
        census["two_roots"]["informative_on_the_capped_root"]
        + census["capped_sum_status"]["capped_root_only"]
    )
    assert would_be == 128, would_be
    assert would_be > 100, "the gate lane-capsumrna registered"
    assert "128" in rp.FORECLOSES and "100" in rp.FORECLOSES


def test_the_module_reads_no_rna_and_no_outcome() -> None:
    """Static: it must not import a measured-RNA reader, so a blind registration stays blind."""
    tree = ast.parse((ROOT / "genomeos/attribution/rootpick.py").read_text())
    imported = {n.module or "" for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)} | {
        a.name for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names
    }
    forbidden = [m for m in imported if any(x in m for x in ("rna", "combine2", "measured", "capsumrna"))]
    assert forbidden == [], f"a blind registration must not import an outcome source: {forbidden}"


def test_the_registration_states_its_falsifier_exposure_and_limits() -> None:
    for name in ("ROOT_RULE", "GROUND", "FALSIFIER", "INCONCLUSIVE", "PRIOR_EXPOSURE", "DOES_NOT_CLAIM"):
        assert getattr(rp, name).strip(), name
    assert "NO RNA value" in rp.PRIOR_EXPOSURE
    assert "would NOT refute it" in rp.FALSIFIER
