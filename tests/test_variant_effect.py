"""Task 1.6: variant effect classification, checked against ClinVar's labels on chr21."""

from collections import Counter
from pathlib import Path

import pytest

from genomeos.genome import (
    Annotation,
    Chromosome,
    Genome,
    IndexedGenome,
    Locus,
    Sequence,
    Variant,
    default_gencode,
    iter_vcf,
)
from genomeos.ir import Transcript
from genomeos.runtime import classify, classify_all

GFF = default_gencode({"chr21", "chrM"}) or Path("missing")
CHR21 = Path("data/reference/chr21.fa.gz")
CLINVAR = Path("data/reference/clinvar_chr21_MT.vcf")


def _toy():
    #            0         1         2         3
    #            0123456789012345678901234567890123
    chrom = "GGATGAAACCCGTAAGTTTTAGGGGTTTTAAGGCC"  # exon1 2-13 (ATGAAACCC GTA) intron exon2 ...
    g = Genome("t")
    g.chromosomes["c"] = Chromosome("c", Sequence(chrom))
    # CDS: ATG AAA CCC | (intron GTAAGTTTTAG) | GGG TTT TAA  -> M K P G F *
    tx = Transcript(
        id="tx1",
        kind="transcript",
        exons=[Locus("c", 2, 11), Locus("c", 22, 31)],
        cds_segments=[Locus("c", 2, 11), Locus("c", 22, 31)],
        cds=Locus("c", 2, 31),
    )
    return g, tx


def test_toy_consequences():
    g, tx = _toy()
    assert (
        classify(g, tx, Variant("c", 5, "A", ("G",), gt=(1, 1))).consequence == "missense_variant"
    )  # AAA->GAA K->E
    assert (
        classify(g, tx, Variant("c", 7, "A", ("G",), gt=(1, 1))).consequence == "synonymous_variant"
    )  # AAA->AAG K
    e = classify(g, tx, Variant("c", 5, "A", ("T",), gt=(1, 1)))  # AAA->TAA stop
    assert e.consequence == "nonsense" and e.protein_change == "K2*"
    assert classify(g, tx, Variant("c", 6, "AA", ("A",), gt=(1, 1))).consequence == "frameshift_variant"
    assert classify(g, tx, Variant("c", 5, "AAAC", ("A",), gt=(1, 1))).consequence == "inframe_deletion"
    assert classify(g, tx, Variant("c", 3, "T", ("C",), gt=(1, 1))).consequence == "start_lost"
    assert classify(g, tx, Variant("c", 15, "G", ("C",), gt=(1, 1))).consequence == "intron_variant"
    assert classify(g, tx, Variant("c", 11, "G", ("C",), gt=(1, 1))).consequence == "splice_site"
    # second exon: GGG -> GAG is G->E (missense); third codon position GGG -> GGA stays G
    assert classify(g, tx, Variant("c", 23, "G", ("A",), gt=(1, 1))).consequence == "missense_variant"
    assert classify(g, tx, Variant("c", 24, "G", ("A",), gt=(1, 1))).consequence == "synonymous_variant"


@pytest.mark.skipif(not (GFF.exists() and CHR21.exists() and CLINVAR.exists()), reason="fetch data first")
def test_clinvar_chr21_consequences_agree():
    ann = Annotation.from_gff3(GFF, {"chr21"})
    m = ann.to_module("chr21")
    genome = IndexedGenome(CHR21)
    coding: dict[str, list[Transcript]] = {}
    for g in ann.protein_coding():
        txs = [
            t
            for t in m.entities[g.id].transcripts
            if t.cds_segments and t.attrs["transcript_type"] == "protein_coding"
        ]
        if txs:
            coding[g.symbol] = txs
    wanted = {"nonsense", "missense_variant", "synonymous_variant", "frameshift_variant"}
    # ClinVar calls any premature termination "nonsense", including frameshifts
    compatible = {
        "nonsense": {"nonsense", "frameshift_variant"},
        "frameshift_variant": {"frameshift_variant"},
        "missense_variant": {"missense_variant"},
        "synonymous_variant": {"synonymous_variant"},
    }
    agree: Counter[str] = Counter()
    total: Counter[str] = Counter()
    confusion: Counter[tuple[str, str]] = Counter()
    for n, v in enumerate(iter_vcf(CLINVAR, {"21"}, pass_only=False)):
        if n % 4:  # a deterministic quarter of the ~53k records keeps the test near ten seconds
            continue
        info = dict(kv.split("=", 1) for kv in v.info.split(";") if "=" in kv)
        labels = {x.split("|")[1] for x in info.get("MC", "").split(",") if "|" in x}
        if len(labels) != 1 or not labels & wanted or len(v.alts) != 1:
            continue
        label = labels.pop()
        genes = [x.split(":")[0] for x in info.get("GENEINFO", "").split("|")]
        txs = [t for g in genes for t in coding.get(g, [])]
        if not txs:
            continue
        v21 = Variant("chr21", v.pos, v.ref, v.alts, gt=(1, 1))
        got = classify_all(genome, txs, v21)[0].consequence
        total[label] += 1
        agree[label] += got in compatible[label]
        if got not in compatible[label]:
            confusion[(label, got)] += 1
    genome.close()
    assert sum(total.values()) > 2000
    for label in wanted:
        assert total[label] > 30, label
        assert agree[label] / total[label] > 0.9, (
            label,
            agree[label],
            total[label],
            confusion.most_common(8),
        )
