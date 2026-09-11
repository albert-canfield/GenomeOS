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


CHROM_LENGTHS = {  # hg38, for windowing the Ensembl gene listing
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
    "chrM": 16569,
}


def coding_symbols_from_ensembl(chrom: str, window: int = 5_000_000) -> list[str]:
    """Protein-coding gene symbols of a chromosome from Ensembl REST (5 Mb windows, with the
    compiler's retry and back-off), for chromosomes without local gene models."""
    from genomeos.molecules.compiler import _get

    name = chrom.removeprefix("chr")
    if name == "M":
        name = "MT"
    length = CHROM_LENGTHS.get(chrom)
    if not length:
        raise ValueError(f"unknown chromosome {chrom}")
    out: set[str] = set()
    for start in range(1, length + 1, window):
        end = min(start + window - 1, length)
        url = (
            f"https://rest.ensembl.org/overlap/region/human/{name}:{start}-{end}"
            "?feature=gene;biotype=protein_coding;content-type=application/json"
        )
        for g in _get(url, timeout=120, retries=5):
            sym = g.get("external_name")
            if sym and not sym.startswith("ENSG"):
                out.add(sym)
    return sorted(out)


def origins_from_annotation(annotation, chrom: str) -> dict[str, dict[str, Any]]:
    """The genomic-origin section from local gene models: no Ensembl call needed (the slowest source)."""
    out: dict[str, dict[str, Any]] = {}
    for g in annotation.genes.values():
        if g.locus.chrom != chrom or g.type != "protein_coding":
            continue
        txs = []
        for t in g.transcripts.values():
            cds_len = sum(loc.length for loc, _phase in t.cds)
            txs.append(
                {
                    "transcript": t.id.split(".")[0],
                    "name": t.name,
                    "biotype": t.type,
                    "canonical": "Ensembl_canonical" in t.tags,
                    "protein": t.id.split(".")[0] if cds_len else None,
                    "protein_length": (cds_len // 3 - 1) if cds_len >= 6 else None,
                }
            )
        txs.sort(key=lambda x: (not x["canonical"], x["name"] or ""))
        out[g.symbol] = {
            "gene_id": g.id.split(".")[0],
            "biotype": g.type,
            "locus": f"{chrom}:{g.locus.start + 1}-{g.locus.end}({g.locus.strand.value})",
            "description": None,
            "transcripts": txs,
            "protein_products": sum(1 for x in txs if x["protein"]),
        }
    return out


def compile_chromosome(
    chrom: str, gff3: str | None = None, limit: int | None = None, log=sys.stdout, workers: int = 8
) -> dict[str, Any]:
    from genomeos.genome import Annotation, default_gencode

    gff = gff3 or default_gencode({chrom})
    origins: dict[str, dict[str, Any]] = {}
    if gff:
        ann = Annotation.from_gff3(gff, {chrom})
        symbols = sorted(
            g.symbol
            for g in ann.genes.values()
            if g.type == "protein_coding" and not g.symbol.startswith("ENSG")
        )
        symbol_source = "GENCODE (local models)"
        origins = origins_from_annotation(ann, chrom)
    else:
        symbols = coding_symbols_from_ensembl(chrom)
        symbol_source = "Ensembl REST gene list"
    print(f"{chrom}: {len(symbols)} protein-coding genes from {symbol_source}", file=log, flush=True)
    if limit:
        symbols = symbols[:limit]
    counts = dict.fromkeys(QUESTIONS, 0)
    per_gene: dict[str, dict[str, bool]] = {}
    no_entry: list[str] = []
    t0 = time.time()
    # the work is waiting on public databases, so several genes compile at once (each writes its
    # own cache file; a gene that fails is retried once, then recorded without an entry)
    from concurrent.futures import ThreadPoolExecutor

    def one(sym: str) -> tuple[str, dict[str, bool]]:
        for attempt in range(2):
            try:
                d = compile_protein(sym, origin=origins.get(sym))
                return sym, d.get("coverage") or coverage(d)
            except Exception:  # noqa: BLE001
                if attempt:
                    return sym, dict.fromkeys(QUESTIONS, False)
                time.sleep(5)
        return sym, dict.fromkeys(QUESTIONS, False)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        for i, (sym, cov) in enumerate(pool.map(one, symbols), 1):
            per_gene[sym] = cov
            if not cov["sequence"]:
                no_entry.append(sym)
            for q in QUESTIONS:
                counts[q] += cov[q]
            if i % 10 == 0 or i == len(symbols):
                print(
                    f"{chrom}: {i}/{len(symbols)} compiled ({time.time() - t0:.0f} s)", file=log, flush=True
                )
                try:
                    from genomeos.jobs import heartbeat

                    heartbeat("proteome_genome_wide")
                except OSError:
                    pass
    n = len(symbols)
    return {
        "chrom": chrom,
        "coding_genes": n,
        "coverage_counts": counts,
        "coverage_fraction": {q: round(counts[q] / n, 3) if n else None for q in QUESTIONS},
        "no_reviewed_entry": no_entry,
        "per_gene": per_gene,
        "seconds": round(time.time() - t0),
        "symbol_source": symbol_source,
        "note": "definitions cached in data/knowledge/proteins/ (not committed); this summary is the result",
    }
