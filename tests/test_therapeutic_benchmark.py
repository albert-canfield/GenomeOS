# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The therapeutic pipeline against targets whose answer is already known.

The benchmark itself needs the network and several minutes
(`scripts/therapeutic_benchmark.py`); these tests read its committed result,
so a regression that loses a known target fails in CI within a second.

Three questions are kept apart on purpose. Recovering the target and calling
the route correctly are asserted, because we believe them. Whether the
mechanism ranked first is defensible is *recorded* rather than asserted,
because the pipeline currently gets it wrong three times out of six and
pretending otherwise would make the benchmark decorative. The count is
pinned: it may fall, never rise.
"""

from __future__ import annotations

import pytest

from genomeos.results import load_result

RESULT = load_result("therapeutic_benchmark")
pytestmark = pytest.mark.skipif(not RESULT, reason="benchmark result not present")

#: Known defects, counted today. Lower this when one is fixed; a rise is a regression.
KNOWN_MECHANISM_DEFECTS = 2  # was 3; the agonist-against-a-driver case was fixed 2026-09-11


def rows():
    return RESULT["rows"]


def test_every_known_target_is_recovered():
    lost = [r["gene"] for r in rows() if not r["recovered"]]
    assert not lost, f"the pipeline lost known targets: {lost}"


def test_the_approved_antibody_targets_are_found_as_surface_targets():
    """EGFR and ERBB2 have approved antibodies. If these are not surface targets, nothing is."""
    surface = [r for r in rows() if r["expected"] == "surface"]
    assert len(surface) >= 2
    for r in surface:
        assert r["pass"], f"{r['gene']} has an approved antibody but was not called a surface target"
        assert (r["surface_accessibility"] or 0) > 0.5


def test_no_surface_route_is_claimed_for_an_intracellular_target():
    """The four small-molecule cases pass by *not* claiming a route that does not exist."""
    out = [r for r in rows() if r["expected"] == "out_of_scope_expected"]
    assert len(out) >= 4
    for r in out:
        assert r["pass"], f"{r['gene']}: {r['verdict']}"
        assert (r["surface_accessibility"] or 0) <= 0.5


def test_an_intracellular_driver_still_gets_the_peptide_route():
    """GenomeOS's distinctive claim: a nuclear or cytoplasmic mutation is never 'no target'."""
    missense = [r for r in rows() if r["expected"] == "out_of_scope_expected"]
    assert all(r.get("peptide_route") for r in missense), (
        "an intracellular missense driver must still offer the peptide/HLA route"
    )


def test_the_known_mechanism_defects_do_not_grow():
    flagged = [r["gene"] for r in rows() if not r.get("mechanism_sane", True)]
    assert len(flagged) <= KNOWN_MECHANISM_DEFECTS, (
        f"the top-ranked mechanism became indefensible for more targets: {flagged}"
    )
    for r in rows():
        assert "mechanism_sane" in r, "the benchmark stopped scoring the mechanism question"


def test_the_benchmark_states_what_it_does_not_cover():
    assert "small molecule" in RESULT["note"].lower()
    assert RESULT["cases"] == len(rows()) >= 6
