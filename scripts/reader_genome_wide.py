"""Reader v1 on every chromosome for a few cell types: fetch each cell type's ENCODE DNase
peaks once (all chromosomes kept), then read every local chromosome. Keeps one summary
(data/results/reader_genome_wide.json) with per-chromosome and genome-wide counts and the
genes read in one cell type but not another."""

from __future__ import annotations

import sys
import time
from pathlib import Path

from genomeos.genome import Annotation, IndexedGenome, default_gencode
from genomeos.genome.domains import infer_domains
from genomeos.genome.reader import compare, fetch_peaks, load_peaks, read_chromosome, slug
from genomeos.genome.regulatory import load_ccres
from genomeos.results import load_result, save_result

ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
CELLS = sys.argv[1:] or ["K562", "HepG2"]


def main() -> None:
    out = load_result("reader_genome_wide") or {"cell_types": CELLS, "chromosomes": {}}
    t0 = time.time()
    local = [c for c in ORDER if Path(f"data/reference/{c}.fa.gz").exists() and default_gencode({c})]
    for cell in CELLS:
        missing = {c for c in local if not load_peaks(cell, c)}
        if missing:
            info = fetch_peaks(cell, missing)
            print(f"{cell}: {info['accession']}, peaks kept for {len(info['kept'])} chromosomes", flush=True)
    for chrom in local:
        if chrom in out["chromosomes"]:
            continue
        ann = Annotation.from_gff3(default_gencode({chrom}), {chrom})
        g = IndexedGenome(f"data/reference/{chrom}.fa.gz")
        length = g.lengths[chrom]
        g.close()
        ccres = load_ccres(chrom)
        domains = infer_domains(chrom, length, ccres, ann) if ccres else []
        rows = {}
        reads = []
        for cell in CELLS:
            r = read_chromosome(cell, chrom, ann, domains, ccres)
            save_result(f"reader_{slug(cell)}_{chrom}", {k: v for k, v in r.items() if k != "_read_all"})
            rows[cell] = {
                k: r[k]
                for k in (
                    "peaks",
                    "coding_genes",
                    "genes_read",
                    "read_fraction",
                    "enhancers_active",
                    "enhancers_active_fraction",
                    "nodes",
                    "nodes_open",
                    "nodes_silent",
                )
            }
            reads.append(r)
        if len(reads) > 1:
            c = compare(reads[0], reads[1])
            rows["only_" + CELLS[0]] = c["read_in_a_only"][:30]
            rows["only_" + CELLS[1]] = c["read_in_b_only"][:30]
            rows["read_in_both"] = c["read_in_both"]
        out["chromosomes"][chrom] = rows
        save_result("reader_genome_wide", out)
        print(
            f"{chrom}: "
            + "; ".join(
                f"{cell} reads {rows[cell]['genes_read']}/{rows[cell]['coding_genes']}" for cell in CELLS
            ),
            flush=True,
        )
    ch = out["chromosomes"]
    out["totals"] = {
        cell: {
            "coding_genes": sum(r[cell]["coding_genes"] for r in ch.values()),
            "genes_read": sum(r[cell]["genes_read"] for r in ch.values()),
            "enhancers_active": sum(r[cell]["enhancers_active"] for r in ch.values()),
            "nodes_silent": sum(r[cell]["nodes_silent"] for r in ch.values()),
        }
        for cell in CELLS
    }
    out["totals"]["seconds"] = round(time.time() - t0)
    out["evidence"] = "experimental: ENCODE DNase-seq peaks; inferred: read = promoter open"
    save_result("reader_genome_wide", out)
    print(f"done: {out['totals']}", flush=True)


if __name__ == "__main__":
    main()
