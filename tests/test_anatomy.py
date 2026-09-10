"""Counting the blocks and elements of a genome."""

from pathlib import Path

import pytest

from genomeos.genome import Annotation, Genome, anatomy_of, default_gencode, design_lessons
from genomeos.genome.anatomy import sequence_elements
from genomeos.results import load_result

CHRM = Path("data/reference/chrM.fa.gz")


def test_sequence_elements_on_synthetic():
    s = "ACGT" * 50 + "A" * 15 + "AC" * 12 + "TTAGGG" * 4 + "N" * 120 + "CG" * 150 + "GC" * 150
    seq, el = sequence_elements(s)
    assert el["homopolymer_runs_ge12"] == 1 and el["microsatellites_dinuc_ge20bp"] >= 1
    assert el["telomeric_blocks"] == 1 and el["telomeric_bp"] == 24
    assert seq["gaps"] == 1 and seq["gap_bp"] == 120
    assert el["cpg_islands"] >= 1  # the CG/GC tail is a CpG island


@pytest.mark.skipif(not (CHRM.exists() and default_gencode({"chrM"})), reason="fetch chrM first")
def test_mitochondrial_anatomy_is_compact():
    ann = Annotation.from_gff3(default_gencode({"chrM"}), {"chrM"})
    seq = Genome.from_fasta(CHRM).chromosomes["chrM"].sequence
    a = anatomy_of("mt", seq, ann, "chrM")
    f = a.composition_fraction()
    assert a.genes["total"] == 37 and a.genes["protein_coding"] == 13
    assert f["cds"] > 0.6 and f["intergenic"] < 0.1 and a.structure["introns_total"] == 0
    assert a.structure["single_exon_fraction"] == 1.0
    assert sum(a.composition_bp.values()) == a.length
    assert "no introns" in design_lessons([a])[0]


def test_saved_comparison_tells_the_organisation_story():
    r = load_result("anatomy_comparison")
    if r is None:
        pytest.skip("run the comparison first")
    by = {a["name"]: a for a in r["anatomies"]}
    mt, human, worm = by["human mtDNA (chrM)"], by["human chr21"], by["C. elegans III"]
    assert (
        mt["composition_fraction"]["cds"]
        > worm["composition_fraction"]["cds"]
        > human["composition_fraction"]["cds"]
    )
    assert human["structure"]["intron_length_median"] > 10 * worm["structure"]["intron_length_median"]
    assert human["composition_fraction"]["intron"] > 0.4
    assert worm["layout"]["coding_genes_per_mb"] > 20 * human["layout"]["coding_genes_per_mb"]
