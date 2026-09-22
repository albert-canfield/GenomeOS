# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Copy-number and structural variants: the knowledge, the reading, the effect.

The distillation itself needs the network (`genomeos cancer alterations`);
these tests read its committed result, so a change that loses the events or
mis-grades them fails in CI within a second. The pipeline tests at the bottom
are the controls that say what the profiles actually add: each fixture's VCF
carries a variant in a different gene, so the gene under test can reach the
candidate list only through its copy-number or structural call.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from genomeos.cancer.alterations import (
    GeneAlteration,
    detect_cna_format,
    grade,
    read_cna_table,
    read_sv_table,
)
from genomeos.results import load_result

K = load_result("cancer_alterations_msk_impact_2017")
MUT = load_result("cancer_msk_impact_2017")

ROOT = Path(__file__).resolve().parent.parent
DEMO = ROOT / "data/demo/alterations"
VEP_CACHE = ROOT / "data/knowledge/vep/vep_cache.jsonl"
PROTEINS = ROOT / "data/knowledge/proteins"

needs_knowledge = pytest.mark.skipif(not K, reason="run genomeos cancer alterations")
needs_caches = pytest.mark.skipif(
    not (DEMO.exists() and VEP_CACHE.exists() and (PROTEINS / "ERBB2.json").exists()),
    reason="needs the cached VEP annotations and compiled proteins",
)


def freq(gene: str, kind: str) -> float:
    return ((K["genes"].get(gene) or {}).get(kind) or {}).get("frequency", 0.0)


# --- the distilled cohort -------------------------------------------------------------


@needs_knowledge
def test_the_study_has_both_profiles():
    assert K["cna_profile"] == "msk_impact_2017_cna"
    assert K["sv_profile"] == "msk_impact_2017_structural_variants"
    assert K["samples"] > 10_000


@needs_knowledge
def test_the_canonical_deep_deletion_and_amplification_are_found():
    """CDKN2A is deleted and ERBB2 amplified, or the table is not describing cancer."""
    assert freq("CDKN2A", "deep_deletion") > 0.05
    assert freq("CDKN2A", "deep_deletion") > freq("CDKN2A", "amplification") * 20
    assert freq("ERBB2", "amplification") > 0.03
    assert freq("PTEN", "deep_deletion") > 0.02
    assert freq("RB1", "deep_deletion") > 0.01


@needs_knowledge
def test_amplification_finds_genes_the_mutation_table_calls_passengers():
    """The reason this table exists: frequency by mutation alone under-counts.

    CCND1 and MYC are mutated in well under 1% of these tumours and amplified
    in several percent. A pipeline reading only the mutation table sees two
    passengers.
    """
    if not MUT:
        pytest.skip("run genomeos cancer distil")
    for gene in ("CCND1", "MYC"):
        mutated = MUT["genes"][gene]["frequency"]
        assert mutated < 0.01
        assert freq(gene, "amplification") > 4 * mutated


@needs_knowledge
def test_the_cancer_type_is_where_the_frequency_lives():
    """ERBB2 amplification is 4% of all tumours and 14% of breast tumours."""
    breast = K["genes"]["ERBB2"]["amplification"]["by_cancer_type"]["Breast Cancer"]
    assert breast > 0.10 > freq("ERBB2", "amplification")
    glioma = K["genes"]["CDKN2A"]["deep_deletion"]["by_cancer_type"]["Glioma"]
    assert glioma > 0.25


@needs_knowledge
def test_the_recurrent_fusions_have_their_known_partners():
    partners = dict(K["genes"]["ALK"]["fusion_partners"])
    assert partners.get("EML4", 0) > 20, "EML4-ALK is the fusion of lung adenocarcinoma"
    ret = dict(K["genes"]["RET"]["fusion_partners"])
    assert {"KIF5B", "CCDC6"} <= set(ret)


@needs_knowledge
def test_a_panel_studys_silence_is_recorded_as_silence():
    """MSK-IMPACT calls only deep deletions and amplifications, not shallow ones."""
    kinds = {k for entry in K["genes"].values() for k in entry}
    assert "gain" not in kinds and "shallow_deletion" not in kinds
    assert "absence of a call" in K["note"]


# --- reading a patient's table --------------------------------------------------------


def test_the_two_copy_number_conventions_are_told_apart_and_never_guessed():
    """A GISTIC 2 is an amplification; 2 copies is an untouched gene."""
    assert detect_cna_format([12.0, 2.0])[0] == "copies"
    assert detect_cna_format([-2.0, 2.0])[0] == "gistic"
    assert detect_cna_format([1.4, 2.0])[0] == "copies"
    fmt, reason = detect_cna_format([0.0, 1.0, 2.0])
    assert fmt == "copies", "an ambiguous table must not be read as a table of amplifications"
    assert "invents no amplification" in reason


def test_copies_and_discrete_calls_produce_the_same_kinds(tmp_path):
    copies = tmp_path / "copies.tsv"
    copies.write_text("gene\tcopies\nERBB2\t12\nPTEN\t0\nCDK4\t3\nFLAT1\t2\n")
    kinds = {a.gene: a.kind for a in read_cna_table(str(copies))}
    assert kinds == {"ERBB2": "amplification", "PTEN": "deep_deletion", "CDK4": "gain"}
    assert "FLAT1" not in kinds, "an untouched gene is not an alteration"

    discrete = tmp_path / "gistic.tsv"
    discrete.write_text("ERBB2\t2\nPTEN\t-2\nCDK4\t1\nFLAT1\t0\n")
    kinds = {a.gene: a.kind for a in read_cna_table(str(discrete))}
    assert kinds == {"ERBB2": "amplification", "PTEN": "deep_deletion", "CDK4": "gain"}
    assert all(a.copies is None for a in read_cna_table(str(discrete))), (
        "a discrete call is a call, not a copy count, and must never be turned into one"
    )


@needs_knowledge
def test_grading_attaches_the_cohort_frequency_and_scores_a_recurrent_fusion_higher():
    amp, fusion, rare = grade(
        [
            GeneAlteration("ERBB2", "amplification", copies=12.0),
            GeneAlteration("ALK", "fusion", partner="EML4"),
            GeneAlteration("ALK", "fusion", partner="NOTAREALGENE"),
        ]
    )[0:3]
    by_gene = {(a.gene, a.partner): a for a in (amp, fusion, rare)}
    eml4 = by_gene[("ALK", "EML4")]
    other = by_gene[("ALK", "NOTAREALGENE")]
    assert eml4.recurrent_partner and not other.recurrent_partner
    assert eml4.score > other.score
    erbb2 = by_gene[("ERBB2", "")]
    assert erbb2.cohort_frequency and erbb2.cohort_frequency > 0.03
    assert any("Breast Cancer" in e or "Esophagogastric" in e for e in erbb2.evidence)


def test_a_deep_deletion_says_it_removes_the_product():
    (a,) = grade([GeneAlteration("CDKN2A", "deep_deletion", gistic=-2)], knowledge={})
    assert a.removes_product
    assert any("never a target for a binder" in e for e in a.evidence)
    assert not GeneAlteration("ERBB2", "amplification").removes_product


def test_a_structural_table_reads_the_partner():
    (a,) = read_sv_table(str(DEMO / "alk_eml4_fusion.sv"))
    assert (a.gene, a.partner, a.kind) == ("ALK", "EML4", "fusion")
    assert a.label() == "ALK-EML4 fusion"


# --- the controls: what the profiles add to the pipeline ------------------------------


def _run(vcf: str, **kw):
    from genomeos.therapeutics import analyse_vcf

    kw.setdefault("indirect", False)
    return analyse_vcf(str(DEMO / vcf), top_genes=8, net=False, **kw)


@needs_caches
@needs_knowledge
def test_an_amplified_oncogene_reaches_the_ranking_with_no_variant_of_its_own():
    """The control. ERBB2 is in this VCF nowhere; the tumour has 12 copies of it.

    Before copy number reached the ranking, the same inputs produced one
    candidate and it was not ERBB2 — the canonical surface target in oncology,
    with an approved antibody, absent because no point mutation touched it.
    """
    a = _run("erbb2_amplification_only.vcf", cnv=str(DEMO / "erbb2_amplification_only.cnv"))
    ranked = [c.gene for c in a["candidates"]]
    assert ranked[0] == "ERBB2"
    erbb2 = a["candidates"][0]
    assert erbb2.target_class == "direct_surface"
    assert not any(o.chromosome for o in erbb2.origins), "no variant of its own"
    assert [o.alteration_kind for o in erbb2.origins] == ["amplification"]
    assert erbb2.origins[0].copy_number == 12.0
    assert erbb2.origins[0].driver_frequency, "the cohort frequency of the amplification is attached"


@needs_caches
@needs_knowledge
def test_a_gene_reached_by_its_own_alteration_is_not_proposed_again_as_a_hypothesis():
    """Found by scoring the copy-number route end to end, which needs `indirect`.

    The control above passes `indirect=False`, so for as long as it was the
    only test of these inputs the defect could not appear. With the pathway
    route on, PIK3CA is mutated and intracellular, ERBB2 is one of its
    physical-evidence partners, and the exclusion set held only the disrupted
    drivers — so ERBB2, already a candidate carrying its 12 copies, was added a
    second time as a hypothesis about PIK3CA. Both entries scored 0.494: the
    amplification bought nothing, which is why the duplicate was invisible in
    the ranking and visible only as a repeated name.

    The weaker entry is the one that has to go. A gene reached by its own
    alteration is not a guess about a neighbour.
    """
    a = _run(
        "erbb2_amplification_only.vcf",
        cnv=str(DEMO / "erbb2_amplification_only.cnv"),
        indirect=True,
    )
    ranked = [c.gene for c in a["candidates"]]
    assert len(ranked) == len(set(ranked)), f"a gene is proposed twice: {ranked}"
    erbb2 = next(c for c in a["candidates"] if c.gene == "ERBB2")
    assert erbb2.target_class == "direct_surface", "the surviving entry is the one with the evidence"
    assert [o.alteration_kind for o in erbb2.origins] == ["amplification"]


@needs_caches
def test_a_deleted_tumour_suppressor_is_never_offered_as_a_target():
    """A homozygous deletion is actionable and it is not a target.

    The gene has a plasma-membrane-adjacent annotation and a 7.6% deletion
    frequency, both of which would have been read as reasons to aim at it. The
    tumour makes none of the protein, so every mechanism that has to recognise
    a product fails, and it fails with that sentence rather than a low score.
    """
    a = _run("cdkn2a_deleted.vcf", cnv=str(DEMO / "cdkn2a_deleted.cnv"))
    genes = {c.gene: c for c in a["candidates"]}
    assert "CDKN2A" in genes, "the deletion still reaches the analysis: it is a finding"
    cdkn2a = genes["CDKN2A"]
    assert cdkn2a.target_class == "unsuitable"
    assert cdkn2a.best_mechanism is None
    failed = [m for m in cdkn2a.therapeutic_mechanisms if m.gates_failed]
    assert any("gene_product_present" in g for m in failed for g in m.gates_failed)
    assert all(m.compatibility == 0.0 for m in cdkn2a.therapeutic_mechanisms)
    # and the amplified gene in the same GISTIC table came through as a target
    assert "EGFR" in genes and genes["EGFR"].target_class == "direct_surface"


@needs_caches
@needs_knowledge
def test_a_fusion_reaches_the_ranking_and_its_junction_is_not_invented():
    a = _run("alk_eml4_fusion.vcf", sv=str(DEMO / "alk_eml4_fusion.sv"))
    alk = next(c for c in a["candidates"] if c.gene == "ALK")
    assert alk.origins[0].alteration_kind == "fusion"
    assert alk.origins[0].fusion_partner == "EML4"
    assert alk.origins[0].recurrent_partner
    assert any("not reconstructed" in limitation for limitation in alk.limitations), (
        "the junction sequence is not reconstructed, and the report has to say so"
    )


@needs_caches
@needs_knowledge
def test_a_fusion_product_is_not_given_the_whole_genes_outside():
    """The same defect as BRAF's, one layer down, found by the benchmark work.

    ALK is a single-pass receptor and every curated compartment says so, so the
    pipeline offered a blocking antibody at 0.75 as EML4-ALK's preferred
    mechanism. In the tumour EML4-ALK is a cytoplasmic kinase: ALK is the 3'
    partner and its extracellular domain is not in the product. There is no
    approved antibody against it and there could not be; crizotinib, alectinib
    and lorlatinib are small molecules that work inside the cell.

    The curated compartment describes the full-length protein, and a fusion
    keeps one side of a junction. When this test was written GenomeOS held no
    orientation, so the surface requirement was unanswered and every surface
    mechanism was provisional. The orientation arrived on 2026-09-21, and an
    unanswered question became an answer: ALK is the 3' partner, so the product
    carries neither its signal peptide nor its ectodomain, and a binder that
    needs an extracellular epitope is refused rather than left open.
    """
    a = _run("alk_eml4_fusion.vcf", sv=str(DEMO / "alk_eml4_fusion.sv"))
    alk = next(c for c in a["candidates"] if c.gene == "ALK")
    assert alk.best_mechanism is None, "an antibody against a cytoplasmic kinase headed this list"
    surface = [m for m in alk.therapeutic_mechanisms if m.mechanism == "blocking_antibody"]
    assert surface and not surface[0].viable, (
        "the surface requirement is answered no, so the mechanism is refused and not merely provisional"
    )
    assert not any(m.viable for m in alk.therapeutic_mechanisms), (
        "every route GenomeOS models needs the outside of the protein; the approved ALK drugs are "
        "small molecules, which is the honest answer here and not a gap"
    )


@needs_caches
@needs_knowledge
def test_the_class_and_the_score_describe_the_product_and_not_the_gene():
    """The half this test used to pin as unreachable, closed by reading one more field.

    It asserted today's wrong answer on purpose — `direct_surface` and an
    accessibility of 1.0, both read off the full-length gene — and said the
    measurement that would close it was in no table GenomeOS reads: the
    junction, the orientation, or transcript evidence for the retained
    domains.

    THAT WAS WRONG, AND IT WAS WRONG BY ONE REQUEST PARAMETER. The cBioPortal
    structural-variant endpoint carries all of it; the client asked for
    `projection=SUMMARY`, which drops it. Under DETAILED the same row gives
    `site1HugoSymbol=EML4` and `site2HugoSymbol=ALK` — site1 is the 5' partner,
    so ALK is the 3' — the two breakpoint positions, and an annotation reading
    "EML4 exons 1-20 with ALK exons 20-29", ALK's ectodomain being exons 1-19.

    A test that records a defect has to name what would close it, and this one
    named something it had not checked was absent.
    """
    a = _run("alk_eml4_fusion.vcf", sv=str(DEMO / "alk_eml4_fusion.sv"))
    alk = next(c for c in a["candidates"] if c.gene == "ALK")
    assert alk.origins[0].fusion_partner == "EML4"
    assert alk.origins[0].fusion_orientation == "3'", "the measurement that closed this"
    assert alk.target_class == "intracellular_only", "classed from the product, not the curated gene"
    assert alk.scores.value("surface_accessibility") == 0.0
    assert "signal peptide" in alk.class_reason and "3'" in alk.class_reason, (
        "the class has to say why, because 'intracellular' about a curated surface receptor is "
        "the surprising answer and the reason is the evidence for it"
    )


def test_a_five_prime_partner_keeps_its_own_ectodomain():
    """The rule is about which end the gene contributes, not about fusions.

    Turning the same fixture round must not produce the same answer: a gene
    contributed as the 5' partner keeps its N-terminus, so its signal peptide
    and ectodomain are in the product and the surface question is open again
    rather than refused. Without this the rule would read 'a fusion is never a
    surface target', which is false — and it is how the previous fix went
    wrong, by answering a question it had only stopped asking.
    """
    import genomeos.therapeutics.mechanisms as mech

    a = _run("alk_eml4_fusion.vcf", sv=str(DEMO / "alk_eml4_fusion.sv"))
    alk = next(c for c in a["candidates"] if c.gene == "ALK")
    assert mech.ectodomain_lost(alk)
    for o in alk.origins:
        o.fusion_orientation = "5'"
    assert not mech.ectodomain_lost(alk), "a 5' partner keeps the N-terminus it contributes"
    for o in alk.origins:
        o.fusion_orientation = ""
    assert not mech.ectodomain_lost(alk), "no orientation reported is unknown, not 'no ectodomain'"


@needs_caches
def test_an_amplification_is_never_read_as_a_measurement_of_protein():
    a = _run("erbb2_amplification_only.vcf", cnv=str(DEMO / "erbb2_amplification_only.cnv"))
    erbb2 = a["candidates"][0]
    assert not erbb2.tumour.expression.known, "copy number is not expression"
    assert erbb2.scores.value("tumour_expression") <= 0.5
    assert "never shows that it does" in erbb2.scores.components["tumour_expression"].basis
    assert any("never shows that it makes it" in x for x in erbb2.limitations)


@needs_caches
def test_the_supplied_alterations_are_reported_as_an_input_level():
    a = _run("erbb2_amplification_only.vcf", cnv=str(DEMO / "erbb2_amplification_only.cnv"))
    assert "copy_number" in a["data_level"]["inputs_present"]
    assert [x["label"] for x in a["alterations"]] == ["ERBB2 amplification (12 copies)"]
    wanted = {m["input"] for m in a["missing_data"]}
    assert "copy number" not in wanted
    assert "structural variants" in wanted


# --- the two open mechanism-ranking defects -------------------------------------------

BENCH = ROOT / "data/demo/benchmark"
needs_bench = pytest.mark.skipif(
    not (BENCH.exists() and VEP_CACHE.exists() and (PROTEINS / "BRAF.json").exists()),
    reason="needs the benchmark fixtures and the caches",
)


def _bench(vcf: str, **kw):
    from genomeos.therapeutics import analyse_vcf

    return analyse_vcf(str(BENCH / vcf), top_genes=4, net=False, indirect=False, **kw)


@needs_bench
def test_an_unanswered_requirement_heads_no_list():
    """The benchmark's two remaining mechanism defects, in one sentence.

    BRAF is cytoplasmic. Its membrane compartment is curated at 0.45, below the
    0.6 reachability threshold, so the surface requirement is neither met nor
    refused. Every surface mechanism was therefore scored provisionally, capped
    at 0.25 — and still headed the list, because everything else scored zero. A
    mechanism that was merely not ruled out is a question, not the best option.
    """
    braf = next(c for c in _bench("braf_v600e_melanoma.vcf")["candidates"] if c.gene == "BRAF")
    assert braf.best_mechanism is None
    nearest = braf.best_provisional_mechanism
    assert nearest is not None, "it is still scored and still listed, with its open requirement"
    assert nearest.provisional_requirements == ["surface_accessible"]
    assert nearest.viable and not nearest.established


@needs_bench
def test_the_answer_appears_as_soon_as_the_data_does():
    """The control that says this is a fix and not a silencing.

    The same cytoplasmic driver, with the patient's HLA genotype supplied,
    gets a preferred mechanism: the peptide/HLA route, which is the correct
    answer for a protein no binder reaches.
    """
    braf = next(
        c for c in _bench("braf_v600e_melanoma.vcf", hla=["HLA-A*02:01"])["candidates"] if c.gene == "BRAF"
    )
    assert braf.best_mechanism is not None
    assert braf.best_mechanism.mechanism in ("tcr_based", "tcr_mimic")
    assert braf.best_mechanism.established


@needs_bench
def test_a_reachable_target_keeps_its_preferred_mechanism():
    """The other side of the control: nothing was taken from the cases that work."""
    egfr = next(c for c in _bench("egfr_l858r_lung.vcf")["candidates"] if c.gene == "EGFR")
    assert egfr.best_mechanism is not None
    assert egfr.best_mechanism.mechanism == "blocking_antibody"
    assert egfr.best_mechanism.established
    assert egfr.best_mechanism.compatibility > 0.5
