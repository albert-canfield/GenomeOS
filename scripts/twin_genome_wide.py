"""The test human on every chromosome: apply HG002's benchmark variants to each local
chromosome, keep the statistics (variants applied per haplotype, SNVs, indels, length
change, reference mismatches), and keep no haplotype FASTA (stream, distil, discard).
Result: data/results/hg002_twin_by_chromosome.json."""

from __future__ import annotations

import time
from pathlib import Path

from genomeos.coords import Locus
from genomeos.genome import IndexedGenome, apply_variants, iter_vcf
from genomeos.genome.fetch import individual_vcf_path
from genomeos.results import load_result, save_result

ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def main() -> None:
    out = load_result("hg002_twin_by_chromosome") or {"sample": "HG002", "chromosomes": {}}
    t0 = time.time()
    for chrom in ORDER:
        if chrom in out["chromosomes"]:
            continue
        fa = Path(f"data/reference/{chrom}.fa.gz")
        vcf = individual_vcf_path(chrom)
        if not fa.exists() or not vcf.exists():
            print(f"{chrom}: skipped (needs the sequence and HG002's rows)", flush=True)
            continue
        t1 = time.time()
        variants = list(iter_vcf(vcf, {chrom}))
        genome = IndexedGenome(fa)
        ref = genome.fetch(Locus(chrom, 0, genome.lengths[chrom]))
        genome.close()
        row = {
            "variants": len(variants),
            "snv": sum(1 for v in variants if v.is_snv),
            "phased": sum(1 for v in variants if v.phased),
        }
        for h in (0, 1):
            _, st = apply_variants(ref, variants, h)
            row[f"hap{h + 1}"] = st
        row["seconds"] = round(time.time() - t1)
        out["chromosomes"][chrom] = row
        save_result("hg002_twin_by_chromosome", out)
        print(
            f"{chrom}: {row['variants']:,} PASS variants; hap1 {row['hap1']}; hap2 {row['hap2']} "
            f"({row['seconds']} s)",
            flush=True,
        )
    ch = out["chromosomes"]
    out["totals"] = {
        "chromosomes": len(ch),
        "variants": sum(r["variants"] for r in ch.values()),
        "snv": sum(r["snv"] for r in ch.values()),
        "reference_mismatches": sum(
            r[h].get("ref_mismatch", 0) for r in ch.values() for h in ("hap1", "hap2")
        ),
        "evidence": "measured: GIAB HG002 v4.2.1 benchmark applied to hg38; haplotype FASTA not kept",
        "seconds": round(time.time() - t0),
    }
    save_result("hg002_twin_by_chromosome", out)
    print(f"done: {out['totals']}", flush=True)


if __name__ == "__main__":
    main()
