"""Tasks 1.1 and 1.4: GENCODE annotation compiles to BioIR and every complete
protein-coding transcript translates cleanly. Needs the downloaded data."""

from pathlib import Path

import pytest

from genomeos.genome import Annotation, Genome
from genomeos.runtime import (
    STANDARD_CODE,
    VERTEBRATE_MITOCHONDRIAL_CODE,
    coding_sequence,
    find_orfs,
    translate,
    translate_transcript,
)

GFF = Path("data/reference/gencode.v50.annotation.gff3.gz")
CHRM = Path("data/reference/chrM.fa.gz")
CHR21 = Path("data/reference/chr21.fa.gz")
needs_data = pytest.mark.skipif(
    not (GFF.exists() and CHRM.exists() and CHR21.exists()), reason="fetch data first"
)


@pytest.fixture(scope="module")
def ann():
    return Annotation.from_gff3(GFF, {"chr21", "chrM"})


@pytest.fixture(scope="module")
def genomes():
    return {"chrM": Genome.from_fasta(CHRM), "chr21": Genome.from_fasta(CHR21)}


@needs_data
def test_gencode_parses_chr21_and_chrm(ann):
    assert 1100 < len(ann.genes) < 1300
    assert ann.summary()["protein_coding"] > 200
    assert ann.gene("APP").locus.chrom == "chr21"  # amyloid precursor protein lives on chr21
    m = ann.to_module("gencode.test")
    assert m.entities[ann.gene("APP").id].kind == "gene"
    assert len(m.entities) > 9000


@needs_data
def test_all_thirteen_mitochondrial_proteins(ann, genomes):
    m = ann.to_module("mt")
    mt = [g for g in ann.protein_coding() if g.locus.chrom == "chrM"]
    assert len(mt) == 13
    expected = {
        "MT-ND1": 318,
        "MT-ND2": 347,
        "MT-CO1": 513,
        "MT-CO2": 227,
        "MT-ATP8": 68,
        "MT-ATP6": 226,
        "MT-CO3": 261,
        "MT-ND3": 115,
        "MT-ND4L": 98,
        "MT-ND4": 459,
        "MT-ND5": 603,
        "MT-ND6": 174,
        "MT-CYB": 380,
    }
    for g in mt:
        tx = m.entities[g.id].transcripts[0]
        prot = translate_transcript(genomes["chrM"], tx, VERTEBRATE_MITOCHONDRIAL_CODE)
        assert prot.startswith("M"), g.symbol
        assert len(prot) == expected[g.symbol], (g.symbol, len(prot))


@needs_data
def test_orf_finder_with_mitochondrial_starts_covers_every_mt_gene(ann, genomes):
    seq = genomes["chrM"].chromosomes["chrM"].sequence
    orfs = list(find_orfs(seq, "chrM", min_aa=50, table=VERTEBRATE_MITOCHONDRIAL_CODE))
    m = ann.to_module("mt")
    for g in (g for g in ann.protein_coding() if g.locus.chrom == "chrM"):
        cds = m.entities[g.id].transcripts[0].cds
        hit = [
            o
            for o in orfs
            if o.locus.strand == cds.strand
            and o.locus.start <= cds.start
            and o.locus.end >= cds.end
            and (o.locus.start - cds.start) % 3 == 0
        ]
        assert hit, g.symbol


@needs_data
def test_complete_protein_coding_transcripts_translate_without_internal_stops(ann, genomes):
    m = ann.to_module("pc")
    ok = bad_start = internal = 0
    for g in ann.protein_coding():
        genome = genomes[g.locus.chrom]
        table = VERTEBRATE_MITOCHONDRIAL_CODE if g.locus.chrom == "chrM" else STANDARD_CODE
        for tx in m.entities[g.id].transcripts:
            if not tx.cds_segments or tx.attrs["transcript_type"] != "protein_coding":
                continue
            if "cds_start_NF" in tx.tags or "cds_end_NF" in tx.tags:
                continue
            full = translate(coding_sequence(genome, tx), table=table, to_stop=False, initiator=True)
            if "*" in full[:-1]:
                internal += 1
            elif not full.startswith("M"):
                bad_start += 1
            else:
                ok += 1
    assert internal == 0
    assert ok > 2500
    assert bad_start / (ok + bad_start) < 0.005  # a handful of GENCODE models start off-ATG
