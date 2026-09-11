# SPDX-License-Identifier: AGPL-3.0-or-later
"""Re-verify GenomeOS translation against UniProt on every local chromosome (no network), then
aggregate the genome-wide summary and triage the disagreements into categories."""

from __future__ import annotations

import glob
import json
import time
from pathlib import Path

from genomeos.genome import Annotation, IndexedGenome, default_gencode
from genomeos.molecules.verify import verify_chromosome
from genomeos.results import load_result, save_result

ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


def triage(row: dict) -> str:
    """Why the reference does not give the curated protein, from machine-checkable signatures:
    an untagged CDS that is not a whole number of codons is a frameshifting reference allele; a
    translation that stops far short of the annotated CDS is a nonsense allele in the reference;
    the rest is isoform or exon-boundary choice, or a different product altogether."""
    sim, ours, uni = row["similarity_best_isoform"], row["ours_aa"], row["uniprot_aa"]
    frac = row.get("translated_fraction_of_cds")
    if row.get("chrom") == "chrM":
        # MT-CO3, MT-CYB, MT-ND1..ND4 end mid-codon by design (the stop is completed by polyadenylation)
        return "mitochondrial incomplete stop codon, completed by polyadenylation: by design"
    if row.get("cds_out_of_frame"):
        return "frameshift allele in hg38: the canonical CDS is not a whole number of codons (untagged)"
    if frac is not None and frac < 0.7 and ours < 0.7 * uni:
        # two-sided: our translation stops early AND the curated protein is substantially longer
        return "nonsense allele in hg38: translation stops before 70% of the CDS and of the curated protein"
    if sim < 0.35:
        return "different protein (UniProt entry or gene model does not describe the same product)"
    if abs(ours - uni) <= 0.3 * max(ours, uni):
        return "same length, different sequence: frame or isoform choice"
    return "different length, partial similarity: exon boundary or isoform choice"


def main() -> None:
    t0 = time.time()
    for chrom in ORDER:
        fa = Path(f"data/reference/{chrom}.fa.gz")
        gff = default_gencode({chrom})
        if not fa.exists() or not gff:
            continue
        ann = Annotation.from_gff3(gff, {chrom})
        g = IndexedGenome(fa)
        try:
            r = verify_chromosome(chrom, ann, g)
        finally:
            g.close()
        save_result(f"translation_vs_uniprot_{chrom}", r)
        print(
            f"{chrom}: {r['genes_checked']} genes, {r['exact_some_isoform_fraction']:.1%} exact for some "
            f"isoform, {r['disagreement_count']} disagreements",
            flush=True,
        )
    rows = []
    for f in glob.glob("data/results/translation_vs_uniprot_chr*.json"):
        d = json.loads(Path(f).read_text())
        for x in d["disagreements"]:
            rows.append({**x, "chrom": d["chrom"], "category": triage(x)})
    cats: dict[str, list[str]] = {}
    for x in rows:
        cats.setdefault(x["category"], []).append(f"{x['gene']} ({x['ours_aa']}/{x['uniprot_aa']} aa)")
    per = {}
    for f in glob.glob("data/results/translation_vs_uniprot_chr*.json"):
        v = json.loads(Path(f).read_text())
        per[v["chrom"]] = {
            "genes": v["genes_checked"],
            "exact_canonical": v["exact_canonical_fraction"],
            "exact_some_isoform": v["exact_some_isoform_fraction"],
            "disagreements": v["disagreement_count"],
        }
    vg = sum(v["genes"] for v in per.values())
    flagged = []
    for f in glob.glob("data/results/translation_vs_uniprot_chr*.json"):
        d = json.loads(Path(f).read_text())
        flagged += [{**x, "chrom": d["chrom"]} for x in d.get("out_of_frame_canonical", [])]
    exact_flagged = [x["gene"] for x in flagged if x["exact_some_isoform"]]
    summary = {
        "out_of_frame_signature": {
            "genes_flagged": len(flagged),
            "translate_exactly_anyway": sorted(exact_flagged),
            "note": "an untagged canonical CDS that is not a whole number of codons; the genes that still "
            "translate exactly are the mitochondrial genes whose CDS ends mid-codon by design",
        },
        "genes_checked": vg,
        "exact_canonical": round(sum(v["exact_canonical"] * v["genes"] for v in per.values()) / vg, 4),
        "exact_some_isoform": round(sum(v["exact_some_isoform"] * v["genes"] for v in per.values()) / vg, 4),
        "disagreements": len(rows),
        "categories": {
            k: {"count": len(v), "genes": sorted(v)}
            for k, v in sorted(cats.items(), key=lambda kv: -len(kv[1]))
        },
    }
    g = load_result("proteome_genome_wide") or {}
    g["translation_vs_uniprot"] = {"summary": summary, "per_chromosome": per}
    save_result("proteome_genome_wide", g)
    print(
        json.dumps({k: v for k, v in summary.items() if k not in ("categories", "out_of_frame_signature")}),
        flush=True,
    )
    o = summary["out_of_frame_signature"]
    print(
        f"  out-of-frame signature: {o['genes_flagged']} genes flagged, "
        f"{len(o['translate_exactly_anyway'])} translate exactly anyway: "
        f"{', '.join(o['translate_exactly_anyway'])}",
        flush=True,
    )
    for k, v in summary["categories"].items():
        print(
            f"  {v['count']:>3}  {k}: {', '.join(v['genes'][:8])}{' …' if v['count'] > 8 else ''}", flush=True
        )
    print(f"done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
