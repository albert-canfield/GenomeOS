# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The third route into the candidate list: a target the tumour never altered.

CD19 and BCMA are not mutated, amplified or rearranged. A pipeline that can
only reach a gene through an alteration cannot propose either, which is a hole
in the middle of the target space. These tests check that the hole is closed
and that closing it did not cost the pipeline its refusals.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from genomeos.therapeutics.atlas import PACKAGED
from genomeos.therapeutics.providers import PatientRnaProvider
from genomeos.therapeutics.scan import LINEAGE_CAVEAT, MIN_RATIO, limitations, scan

ROOT = Path(__file__).resolve().parent.parent
DEMO = ROOT / "data/demo/expression"
VCF = DEMO / "bcell_lymphoma.vcf"
RNA = DEMO / "bcell_lymphoma_rna.tsv"
PROTEINS = ROOT / "data/knowledge/proteins"
VEP_CACHE = ROOT / "data/knowledge/vep/vep_cache.jsonl"

pytestmark = pytest.mark.skipif(
    not (PACKAGED.exists() and RNA.exists()), reason="run scripts/normal_tissue.py"
)


def hits():
    return scan(PatientRnaProvider.from_file(RNA), limit=6)


def test_a_gene_the_tumour_never_altered_becomes_a_target():
    passing = [h for h in hits() if h.passes]
    assert "CD19" in {h.gene for h in passing}
    cd19 = next(h for h in passing if h.gene == "CD19")
    assert cd19.ratio > MIN_RATIO
    assert cd19.normal_tissue == "spleen", "the healthy home of the lineage, and it is named"
    assert cd19.cd_marker


def test_nothing_is_proposed_without_the_patients_own_rna():
    """A cohort describes a cancer type; it cannot say what this tumour displays."""
    assert scan(PatientRnaProvider()) == []


def test_the_classic_b_cell_targets_are_turned_down_and_the_reason_is_kept():
    """The instructive half of the screen, and a real limit of the measure.

    CD20, CD79A, CD79B and CD22 are the targets of the most successful
    antibodies in oncology and every one of them fails a selectivity screen,
    because a lineage antigen has a healthy home full of the same lineage. The
    screen reports them as turned down with the number, rather than dropping
    them silently, because the threshold is a sort order and not a verdict.
    """
    turned_down = {h.gene: h for h in hits() if not h.passes}
    assert {"MS4A1", "CD79A", "CD79B", "CD22"} <= set(turned_down)
    ms4a1 = turned_down["MS4A1"]
    assert ms4a1.normal_max > 200  # healthy spleen carries more than this tumour's RNA does
    assert "not about usefulness" in ms4a1.rejected_because
    assert "247" in ms4a1.rejected_because or f"{ms4a1.normal_max:g}" in ms4a1.rejected_because


def test_every_expression_candidate_carries_what_the_cd19_story_costs():
    cd19 = next(h for h in hits() if h.gene == "CD19")
    notes = limitations(cd19)
    assert LINEAGE_CAVEAT in notes
    assert any("does not show that the protein is made" in n for n in notes)
    assert any("different normalisations" in n for n in notes)


@pytest.mark.skipif(
    not (VCF.exists() and VEP_CACHE.exists() and (PROTEINS / "CD19.json").exists()),
    reason="needs the cached VEP annotations and compiled proteins",
)
def test_the_scan_reaches_the_ranking_and_says_it_is_unaltered():
    from genomeos.therapeutics import analyse_vcf

    a = analyse_vcf(str(VCF), rna=str(RNA), top_genes=4, net=False, indirect=False, scan_expression=6)
    genes = {c.gene: c for c in a["candidates"]}
    assert "CD19" in genes, "unreachable before the scan: no alteration touches it"
    cd19 = genes["CD19"]
    assert cd19.target_class == "direct_surface"
    assert cd19.origins == [], "it has no tumour-DNA origin, and the report must not invent one"
    assert "not altered in this tumour" in cd19.class_reason
    assert LINEAGE_CAVEAT in cd19.limitations
    assert any(m.what == "surface protein in this tumour" for m in cd19.ledger.missing)
    assert [h["gene"] for h in a["expression_hits"] if not h["passes"]], "rejections are reported"


def test_the_scan_is_off_unless_asked_for():
    from genomeos.therapeutics import analyse_vcf

    if not (VCF.exists() and VEP_CACHE.exists()):
        pytest.skip("needs the cached VEP annotations")
    a = analyse_vcf(str(VCF), rna=str(RNA), top_genes=4, net=False, indirect=False)
    assert a["expression_hits"] == []
    assert "CD19" not in {c.gene for c in a["candidates"]}
