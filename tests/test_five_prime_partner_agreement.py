# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The class and the reason agree with the surface score on a fusion's 5' partner.

The structural-variant case (lane-svcase, `2b5a967`,
data/results/therapeutic_sv_case_result.json) is the first benchmark row to reach
`product_extent_unknown`, and its FGFR3 row published two defects:

1. `target_class` read `direct_surface` while `surface_accessibility` was `None`.
   The class was read from the curated annotation and the score from the
   product, so one row named a surface target and carried no surface figure.
2. The score's reason said "no junction position is reported" on a record that
   reports both breakpoints (chr4 1,806,934 to 1,739,702, GRCh37), carried
   verbatim into `fusion_junction`. What is missing is the reading of the
   junction, not the junction: the case result's own manifest says "no
   breakpoint is read into exon content, so the product's extent stays
   unknown".

The class keeps the module's existing vocabulary: `unknown` is the class that
says accessibility cannot be established, so it is the honest value wherever
the score refuses a surface reading for a stated reason. No class is added.

Built from constructed candidates so they run without caches, plus one on the
committed fixture that skips when the caches are absent.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import genomeos.therapeutics.mechanisms as mech
from genomeos.therapeutics.model import Localisation, Region, TherapeuticTargetCandidate, VariantOrigin
from genomeos.therapeutics.pipeline import classify, score_candidate

ROOT = Path(__file__).resolve().parent.parent
FIXTURES = ROOT / "data/demo/alterations"
VEP_CACHE = ROOT / "data/knowledge/vep/vep_cache.jsonl"
PROTEINS = ROOT / "data/knowledge/proteins"

#: The case record's junction, as `read_sv_table` carries it into the alteration.
JUNCTION = "4:1806934 to 4:1739702 (GRCh37)"
DENIAL = "no junction position is reported"


def fgfr3_like() -> Localisation:
    """FGFR3's curated arrangement as the case reads it: type I, signal peptide 1-22."""
    return Localisation(
        compartments={"plasma_membrane": 0.95, "cell_surface": 0.9},
        plasma_membrane=True,
        extracellular_regions=[Region("topo_dom", 23, 375, "Extracellular")],
        transmembrane_regions=[Region("transmem", 376, 396, "Helical")],
        intracellular_regions=[Region("topo_dom", 397, 806, "Cytoplasmic")],
        signal_peptide=Region("signal_peptide", 1, 22, "Signal peptide"),
        topology="single-pass type I",
        orientation="N-terminus outside",
    )


def fusion(orientation: str, junction: str = "") -> VariantOrigin:
    o = VariantOrigin(gene="FGFR3", alteration_kind="fusion")
    o.fusion_orientation = orientation
    o.fusion_partner = "TACC3"
    o.fusion_junction = junction
    return o


def candidate(loc: Localisation, origins: list[VariantOrigin]) -> TherapeuticTargetCandidate:
    return TherapeuticTargetCandidate(gene="FGFR3", localization=loc, origins=origins)


def surface(c: TherapeuticTargetCandidate) -> tuple[float | None, str]:
    comp = score_candidate(c, precedent_available=False).components["surface_accessibility"]
    # A dimension answered None keeps its sentence as `unknown_reason`, not `basis`.
    return comp.value, comp.basis or comp.unknown_reason


# ---- defect 1: the class ------------------------------------------------------------


def _refused_states() -> dict[str, TherapeuticTargetCandidate]:
    """Every fusion state in which the score answers None for a stated reason."""
    uncurated = fgfr3_like()
    uncurated.signal_peptide = None
    return {
        "5' partner, junction reported": candidate(fgfr3_like(), [fusion("5'", JUNCTION)]),
        "5' partner, no junction": candidate(fgfr3_like(), [fusion("5'")]),
        "orientation unreported": candidate(fgfr3_like(), [fusion("")]),
        "5' of one fusion, 3' of another": candidate(fgfr3_like(), [fusion("5'"), fusion("3'")]),
        "3' partner, arrangement uncurated": candidate(uncurated, [fusion("3'")]),
    }


@pytest.mark.parametrize("state", sorted(_refused_states()))
def test_no_surface_class_where_the_score_refuses_a_surface_reading(state):
    """The FGFR3 row's first defect, and every state of the same shape.

    Where the score answers None for a reason, the class may not say
    `direct_surface`: the class is `unknown`, and its reason carries the
    score's own sentence so the two fields of one row say the same thing.
    """
    c = _refused_states()[state]
    value, basis = surface(c)
    assert value is None, f"{state}: the score is expected to refuse here, got {value}"
    klass, reason = classify(c, "altered_gene")
    assert klass != "direct_surface", f"{state}: classed a surface target with no surface figure"
    assert klass == "unknown", f"{state}: {klass}"
    assert basis in reason, f"{state}: the class reason does not carry the score's reason"


def test_the_five_prime_branch_is_the_one_the_case_reaches():
    c = _refused_states()["5' partner, junction reported"]
    assert mech.product_extent_unknown(c)
    assert not mech.ectodomain_lost(c) and not mech.ectodomain_uncertain(c)


def test_a_surface_receptor_the_score_answers_is_still_a_surface_target():
    """The control: the repair is the width of the refusal and no wider."""
    c = candidate(fgfr3_like(), [])
    value, _ = surface(c)
    assert value is not None and value > 0.5
    assert classify(c, "altered_gene")[0] == "direct_surface"


def test_a_three_prime_partner_with_a_curated_arrangement_is_still_intracellular():
    """`ectodomain_lost` keeps its answer of 0.0 and its class."""
    c = candidate(fgfr3_like(), [fusion("3'")])
    assert surface(c)[0] == 0.0
    assert classify(c, "altered_gene")[0] == "intracellular_only"


# ---- defect 2: the reason -----------------------------------------------------------


def test_a_reported_junction_is_named_and_never_denied():
    """The FGFR3 row's second defect: the reason must say what is actually missing."""
    value, basis = surface(candidate(fgfr3_like(), [fusion("5'", JUNCTION)]))
    assert value is None
    assert DENIAL not in basis, "the record reports the junction; the reason says it does not"
    assert JUNCTION in basis, "the junction the record states is quoted as stated"
    assert "exon content" in basis, "what is missing is the reading of the junction into exons"


def test_an_unreported_junction_is_still_said_to_be_unreported():
    value, basis = surface(candidate(fgfr3_like(), [fusion("5'")]))
    assert value is None
    assert DENIAL in basis


def test_a_junction_reported_for_some_fusions_says_how_many():
    value, basis = surface(candidate(fgfr3_like(), [fusion("5'", JUNCTION), fusion("")]))
    assert value is None
    assert DENIAL not in basis
    assert JUNCTION in basis and "1 of 2" in basis


# ---- the committed case -------------------------------------------------------------

needs_caches = pytest.mark.skipif(
    not (FIXTURES.exists() and VEP_CACHE.exists() and (PROTEINS / "FGFR3.json").exists()),
    reason="needs the cached VEP annotations and the compiled FGFR3 protein",
)


@needs_caches
def test_the_committed_case_row_agrees_with_itself(tmp_path):
    """The FGFR3-TACC3 fixture, offline, as the case runner scores the ten rows."""
    from genomeos.therapeutics import analyse_vcf

    a = analyse_vcf(
        str(FIXTURES / "fgfr3_tacc3_urothelial.vcf"),
        sv=str(FIXTURES / "fgfr3_tacc3_urothelial.sv"),
        top_genes=8,
        net=False,
        indirect=True,
        log=None,
        vep_append_to=tmp_path,
    )
    fgfr3 = next(c for c in a["candidates"] if c.gene == "FGFR3")
    assert [o.fusion_junction for o in fgfr3.origins] == [JUNCTION]
    comp = fgfr3.scores.components["surface_accessibility"]
    assert comp.value is None
    reason = comp.unknown_reason
    assert DENIAL not in reason and JUNCTION in reason
    assert fgfr3.target_class == "unknown"
    assert reason in fgfr3.class_reason
