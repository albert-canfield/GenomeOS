"""Phase 5: a second organism compiles through the same pipeline (C. elegans, Ensembl WBcel235)."""

from pathlib import Path

import pytest

from genomeos.genome import Annotation, Genome
from genomeos.runtime import coding_sequence, translate

GFF = Path("data/reference/celegans/WBcel235.63.gff3.gz")
FA = Path("data/reference/celegans/WBcel235.dna.toplevel.fa.gz")


@pytest.mark.skipif(not (GFF.exists() and FA.exists()), reason="download the C. elegans files first")
def test_celegans_chromosome_iii_compiles_and_translates():
    ann = Annotation.from_gff3(GFF, {"III"})
    assert ann.summary()["protein_coding"] > 2000
    genome = Genome.from_fasta(FA)
    assert "III" in genome.chromosomes and 13_000_000 < genome.chromosomes["III"].length < 14_500_000
    m = ann.to_module("celegans.III")
    ok = internal = no_start = 0
    for g in ann.protein_coding():
        for tx in m.entities[g.id].transcripts:
            if not tx.cds_segments or tx.attrs["transcript_type"] != "protein_coding":
                continue
            full = translate(coding_sequence(genome, tx), to_stop=False, initiator=True)
            if "*" in full[:-1]:
                internal += 1
            elif not full.startswith("M"):
                no_start += 1
            else:
                ok += 1
    total = ok + internal + no_start
    assert total > 2000
    assert internal / total < 0.01, (ok, internal, no_start)
    assert ok / total > 0.97
    # a famous worm gene: unc-22 (twitchin) is on IV, but III carries lin-12 (Notch) and egl-1 among others
    assert ann.gene("lin-12").type == "protein_coding"
