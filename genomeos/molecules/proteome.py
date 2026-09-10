"""Compile the proteins of a chromosome and keep one coverage summary.

Every protein-coding gene of the chromosome goes through the compiler; the
compiled definitions stay in the local knowledge cache and the coverage
table (how many proteins have a sequence, a structure, a function, ...) is
saved as a result. This is the measurable version of "how much of the
proteome do we know".
"""

from __future__ import annotations

import sys
import time
from typing import Any

from genomeos.molecules.compiler import compile_protein, coverage

QUESTIONS = [
    "genomic_origin",
    "sequence",
    "name",
    "function",
    "domains",
    "pathways",
    "interactions",
    "expression",
    "structure_experimental",
    "structure_predicted",
    "disease",
]


def compile_chromosome(
    chrom: str, gff3: str | None = None, limit: int | None = None, log=sys.stdout
) -> dict[str, Any]:
    from genomeos.genome import Annotation, default_gencode

    ann = Annotation.from_gff3(gff3 or default_gencode({chrom}), {chrom})
    symbols = sorted(
        g.symbol for g in ann.genes.values() if g.type == "protein_coding" and not g.symbol.startswith("ENSG")
    )
    if limit:
        symbols = symbols[:limit]
    counts = dict.fromkeys(QUESTIONS, 0)
    per_gene: dict[str, dict[str, bool]] = {}
    no_entry: list[str] = []
    t0 = time.time()
    for i, sym in enumerate(symbols, 1):
        d = compile_protein(sym)
        cov = d.get("coverage") or coverage(d)
        per_gene[sym] = cov
        if not cov["sequence"]:
            no_entry.append(sym)
        for q in QUESTIONS:
            counts[q] += cov[q]
        if i % 10 == 0 or i == len(symbols):
            print(f"{chrom}: {i}/{len(symbols)} compiled ({time.time() - t0:.0f} s)", file=log, flush=True)
    n = len(symbols)
    return {
        "chrom": chrom,
        "coding_genes": n,
        "coverage_counts": counts,
        "coverage_fraction": {q: round(counts[q] / n, 3) if n else None for q in QUESTIONS},
        "no_reviewed_entry": no_entry,
        "per_gene": per_gene,
        "seconds": round(time.time() - t0),
        "note": "definitions cached in data/knowledge/proteins/ (not committed); this summary is the result",
    }
