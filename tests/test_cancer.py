"""Healthy vs tumour on chromosome 21 with the distilled cBioPortal knowledge."""

from pathlib import Path

import pytest

from genomeos.cancer import agent_packet, annotate, somatic, suggest_cancer_type, surface_targets
from genomeos.cancer.cbioportal import DRIVER_PANEL
from genomeos.genome import Annotation, IndexedGenome, default_gencode
from genomeos.results import load_result

CHR21 = Path("data/reference/chr21.fa.gz")
K = load_result("cancer_msk_impact_2017")


def test_panel_and_distilled_knowledge():
    assert "TP53" in DRIVER_PANEL and "RUNX1" in DRIVER_PANEL
    if K is None:
        pytest.skip("run genomeos cancer distil")
    assert K["samples"] > 10_000 and K["genes"]["TP53"]["frequency"] > 0.3
    assert K["genes"]["KRAS"]["hotspots"][0][0] in ("G12D", "G12V", "G12C")
    assert K["genes"]["KRAS"]["by_cancer_type"]["Pancreatic Cancer"] > 0.5  # KRAS in most pancreatic tumours


@pytest.mark.skipif(
    not (CHR21.exists() and default_gencode({"chr21"}) and K), reason="needs chr21 + distilled knowledge"
)
def test_somatic_comparison_ranks_a_runx1_truncation(tmp_path):
    ann = Annotation.from_gff3(default_gencode({"chr21"}), {"chr21"})
    g = IndexedGenome(CHR21)
    runx1 = ann.gene("RUNX1")
    m = ann.to_module("t")
    tx = next(t for t in m.entities[runx1.id].transcripts if "Ensembl_canonical" in t.tags)
    # build a nonsense variant in RUNX1 (minus strand): find a codon we can turn into a stop
    from genomeos.coords import Locus
    from genomeos.genome import Variant
    from genomeos.runtime import classify

    hit = None
    for seg in sorted(tx.cds_segments, key=lambda l: l.start):
        for pos in range(seg.start, seg.end):
            ref = str(g.fetch(Locus("chr21", pos, pos + 1)))
            for alt in "ACGT":
                if alt == ref:
                    continue
                e = classify(g, tx, Variant("chr21", pos, ref, (alt,), gt=(1, 1)))
                if e.consequence == "nonsense":
                    hit = (pos, ref, alt)
                    break
            if hit:
                break
        if hit:
            break
    assert hit
    pos, ref, alt = hit
    header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    germline = f"chr21\t{runx1.locus.start - 5000}\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n"
    normal = tmp_path / "normal.vcf"
    tumour = tmp_path / "tumour.vcf"
    normal.write_text(header + germline)
    tumour.write_text(header + germline + f"chr21\t{pos + 1}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t0/1\n")
    som = somatic(str(normal), str(tumour), {"chr21"})
    assert len(som) == 1  # the shared germline variant is removed
    ranked = annotate(som, ann, g, K)
    g.close()
    top = ranked[0]
    assert top.gene == "RUNX1" and top.consequence == "nonsense" and top.driver_frequency
    assert top.score > 1.0 and any("cBioPortal" in e for e in top.evidence)
    types = suggest_cancer_type({"RUNX1"}, K)
    assert types and all("inferred" in t["evidence"] for t in types)
    targets = surface_targets({"EGFR", "TP53", "RUNX1"})
    assert [t["gene"] for t in targets] == ["EGFR"]  # only EGFR is on the cell surface
    packet = agent_packet("tumour", ranked, types, targets)
    assert packet["inputs"]["somatic_variants_ranked"][0]["gene"] == "RUNX1"
    assert "not clinical advice" in " ".join(packet["constraints"])
