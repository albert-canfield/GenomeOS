# SPDX-License-Identifier: AGPL-3.0-or-later
"""Reader v1 on every chromosome for a few cell types: fetch each cell type's ENCODE DNase
peaks once (all chromosomes kept), then read every local chromosome. Keeps one summary
(data/results/reader_genome_wide.json) with per-chromosome and genome-wide counts and the
genes read in one cell type but not another."""

from __future__ import annotations

import gzip
import sys
import time
from pathlib import Path

from genomeos import manifest as mf
from genomeos.genome import Annotation, IndexedGenome, default_gencode
from genomeos.genome import epigenome as ep
from genomeos.genome import reader as rd
from genomeos.genome.domains import infer_domains
from genomeos.genome.reader import compare, fetch_peaks, load_peaks, normalise_family, read_chromosome, slug
from genomeos.genome.regulatory import load_ccres
from genomeos.jobs import heartbeat
from genomeos.results import load_result, save_result

ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
# the two lines of the first run, then a spread of tissues and lineages with released GRCh38 DNase peaks
DEFAULT_CELLS = [
    "K562",
    "HepG2",
    "GM12878",
    "H1",
    "IMR-90",
    "SK-N-SH",
    "cardiac muscle cell",
    "keratinocyte",
    "hepatocyte",
    "astrocyte",
    "CD14-positive monocyte",
]
RERUN = "--rerun" in sys.argv  # start the summary afresh (after a change to what "read" means)
CELLS = [a for a in sys.argv[1:] if a != "--rerun"] or DEFAULT_CELLS

ENTRY = "scripts/reader_genome_wide.py"
#: the code this roll-up is answerable for: itself, the reader and the layer the marks come from
OWN_CODE = (ENTRY, "genomeos/genome/reader.py", "genomeos/genome/epigenome.py")


def chromosome_manifest(cell: str, chrom: str) -> dict:
    """The provenance of one per-chromosome reading (2026-10-10, lane-readerneural).

    A per-chromosome name that is not on the legacy allowlist (any biosample added after the result
    contract, e.g. spinal cord) enters the registry only with a complete manifest, so every reading
    carries one: the files `read_chromosome` and its setup open for this chromosome, each hashed. The
    reading itself is unchanged; this records what it read."""
    gff = default_gencode({chrom})
    ref = Path(f"data/reference/{chrom}.fa.gz")
    files = [gff, ref, Path(f"{ref}.fai"), Path(f"{ref}.gzi"), Path(f"data/results/ccres_{chrom}.bed.gz")]
    files += [rd.peaks_path(cell, chrom)] + [ep.peaks_path(cell, m, chrom) for m in rd.POISED_MARKS]
    # the accession from this chromosome's own peak file, which is declared below; the layer's
    # dnase_accession opens another chromosome's file, which this result did not read
    with gzip.open(rd.peaks_path(cell, chrom), "rt") as fh:
        head = fh.readline()
    accession = head[1:].split()[0] if head.startswith("#") else "unknown: no header line"
    return {
        "sources": [
            {
                "accession": f"ENCODE DNase-seq narrowPeak for {cell}: {accession}",
                "version": "released GRCh38, the largest peak file for the biosample term "
                "(reader.find_dnase_file)",
                "licence": "ENCODE data use policy",
            },
            {
                "accession": f"ENCODE Histone ChIP-seq peaks for {cell} ({', '.join(rd.POISED_MARKS)}), "
                "as pinned in data/results/epigenome_manifest.json",
                "version": "default analysis, replicated before pseudoreplicated narrowPeak",
                "licence": "ENCODE data use policy",
            },
            {"accession": "GENCODE v50 comprehensive annotation (GFF3)", "version": "v50"},
            {"accession": "GRCh38 reference sequence", "version": "GRCh38"},
            {
                "accession": "ENCODE cCREs (registry)",
                "version": "as cached in data/results/ccres_<chrom>.bed.gz",
            },
        ],
        "inputs": [
            mf.files_entry(f"reader inputs for {cell} on {chrom}", [f for f in files if Path(f).exists()])
        ],
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "code": {},
        "parameters": {
            "promoter_window_bp": rd.PROMOTER_WINDOW,
            "node_open_min_density_per_100kb": rd.NODE_OPEN_MIN_DENSITY,
            "node_open_reference_density_per_100kb": rd.NODE_OPEN_REFERENCE_DENSITY,
            "marks": list(rd.POISED_MARKS),
        },
        "exclusions": ["genes that are not protein_coding in GENCODE v50 are not called"],
        "partitions": "n/a: a reading of one biosample on one chromosome, no evaluation split",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE),
    }


def panel_or_refuse(cells: list[str], published: list[str], what: str) -> None:
    """Refuse a roll-up whose panel is narrower than the published one, naming what would go.

    The normalised readings are standardised residuals of a fit over the panel of biosamples
    (`normalise_family` in `genome/reader.py`), so the panel is their denominator: drop one
    biosample and a, b and s all move and every residual in the file is a different number. The
    module already refuses a panel with a missing *field*, on the stated ground that "a panel with
    a hole is a different panel"; a panel with a missing *biosample* is the same hole and nothing
    refused it. Two routes reach it, and both run without a word today: `--rerun` starts the
    summary afresh from the biosamples named on the command line, and DEFAULT_CELLS is narrower
    than the panel on disk; and the completeness filter below drops any biosample absent from one
    chromosome's row. So this is checked, not assumed, and a narrower panel is a decision that has
    to be made out loud.
    """
    dropped = [c for c in published if c not in cells]
    if dropped:
        raise ValueError(
            f"{what} would normalise over {len(cells)} biosamples and drop {len(dropped)} of the "
            f"{len(published)} the published panel was fitted over: {', '.join(dropped)}. Every "
            "residual in `normalised` is standardised by a fit over the panel, so this moves all "
            "of them and says nothing. Name every biosample of the published panel to keep it; a "
            "narrower panel is a different panel and belongs in a result of its own, not on top "
            "of this one."
        )


def _sum_or_none(ch: dict, cell: str, field: str) -> int | None:
    vals = [r[cell].get(field) for r in ch.values()]
    return None if any(v is None for v in vals) else sum(vals)


def main() -> None:
    previous = load_result("reader_genome_wide") or {}
    published = list(previous.get("cell_types") or [])
    # a --rerun overwrites the summary as each chromosome finishes, so its panel is checked before
    # the first write rather than after the last
    if RERUN:
        panel_or_refuse(CELLS, published, "--rerun")
    out = (None if RERUN else previous) or {"cell_types": [], "chromosomes": {}}
    out["cell_types"] = list(dict.fromkeys([*published, *CELLS]))
    t0 = time.time()
    local = [c for c in ORDER if Path(f"data/reference/{c}.fa.gz").exists() and default_gencode({c})]
    for cell in CELLS:
        missing = {c for c in local if not load_peaks(cell, c)}
        if missing:
            info = fetch_peaks(cell, missing)
            print(f"{cell}: {info['accession']}, peaks kept for {len(info['kept'])} chromosomes", flush=True)
    for chrom in local:
        rows = out["chromosomes"].get(chrom, {})
        todo = [cell for cell in CELLS if cell not in rows]
        if not todo:
            continue
        # the window each per-chromosome result is reconciled against opens with its own reading
        mf.trace_begin()
        ann = Annotation.from_gff3(default_gencode({chrom}), {chrom})
        g = IndexedGenome(f"data/reference/{chrom}.fa.gz")
        length = g.lengths[chrom]
        g.close()
        ccres = load_ccres(chrom)
        domains = infer_domains(chrom, length, ccres, ann) if ccres else []
        reads = []
        for cell in todo:
            heartbeat("reader_genome_wide")
            r = read_chromosome(cell, chrom, ann, domains, ccres)
            save_result(
                f"reader_{slug(cell)}_{chrom}",
                {k: v for k, v in r.items() if k != "_read_all"},
                manifest=chromosome_manifest(cell, chrom),
            )
            rows[cell] = {
                k: r[k]
                for k in (
                    "peaks",
                    "coding_genes",
                    "genes_read",
                    "genes_read_open",
                    "genes_poised",
                    "genes_read_by_marks",
                    "read_fraction",
                    "enhancers_active",
                    "enhancers_active_fraction",
                    "enhancers_active_per_100k_peaks",
                    "nodes",
                    "nodes_open",
                    "nodes_open_at_reference",
                    "h3k27me3_peaks",
                    "nodes_silent",
                )
            }
            reads.append(r)
        if len(reads) > 1 and "read_in_both" not in rows:
            c = compare(reads[0], reads[1])
            rows["only_" + todo[0]] = c["read_in_a_only"][:30]
            rows["only_" + todo[1]] = c["read_in_b_only"][:30]
            rows["read_in_both"] = c["read_in_both"]
        out["chromosomes"][chrom] = rows
        save_result("reader_genome_wide", out)
        print(
            f"{chrom}: "
            + "; ".join(
                f"{cell} reads {rows[cell]['genes_read']}/{rows[cell]['coding_genes']}" for cell in todo
            ),
            flush=True,
        )
    ch = out["chromosomes"]
    cells = [c for c in out["cell_types"] if all(c in r for r in ch.values())]
    panel_or_refuse(cells, published, "the roll-up")
    out["cell_types"] = cells
    out["totals"] = {
        cell: {
            "coding_genes": sum(r[cell]["coding_genes"] for r in ch.values()),
            "genes_read": sum(r[cell]["genes_read"] for r in ch.values()),
            "genes_read_open": sum(
                r[cell].get("genes_read_open", r[cell]["genes_read"]) for r in ch.values()
            ),
            "genes_poised": sum(r[cell].get("genes_poised") or 0 for r in ch.values()),
            "genes_read_by_marks": sum(r[cell].get("genes_read_by_marks") or 0 for r in ch.values()),
            "marks_used_on": sum(1 for r in ch.values() if r[cell].get("genes_poised") is not None),
            "enhancers_active": sum(r[cell]["enhancers_active"] for r in ch.values()),
            # the count is substantially DNase depth, so the rate travels with it; `peaks` is
            # carried into the totals so any consumer can re-derive or re-normalise it
            "peaks": sum(r[cell]["peaks"] for r in ch.values()),
            "enhancers_active_per_100k_peaks": round(
                sum(r[cell]["enhancers_active"] for r in ch.values())
                / max(1, sum(r[cell]["peaks"] for r in ch.values()))
                * 100_000,
                1,
            ),
            "nodes_silent": sum(r[cell]["nodes_silent"] for r in ch.values()),
            # the raw counts and covariates the normalised readings are fitted on (2026-09-27);
            # None where a row predates them, and normalise_family then leaves that reading None
            "nodes_open_at_reference": _sum_or_none(ch, cell, "nodes_open_at_reference"),
            "h3k27me3_peaks": _sum_or_none(ch, cell, "h3k27me3_peaks"),
        }
        for cell in cells
    }
    # beside the raw totals, never instead of them: a residual belongs to this panel of biosamples
    out["normalised"] = normalise_family(out["totals"])
    out["totals"]["seconds"] = round(time.time() - t0)
    out["evidence"] = (
        "experimental: ENCODE DNase-seq peaks and Histone ChIP-seq peaks; inferred: read = promoter open "
        "and not poised (H3K27me3 without H3K27ac), or closed with H3K4me3 and H3K27ac, where the marks "
        "were read"
    )
    save_result("reader_genome_wide", out)
    print(f"done: {out['totals']}", flush=True)


if __name__ == "__main__":
    main()
