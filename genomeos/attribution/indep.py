# SPDX-License-Identifier: AGPL-3.0-or-later
"""The attribution layer held against evidence no script in this repository had read: one benchmark.

Review item 11 (2026-09-28) asked for "one independent biological benchmark before expanding
simulation scope" and for generalisation "on genuinely independent evidence before claiming discovery
performance". lane-split's audit found the ENCODE CRISPRi held-out file reused by at least ten scored
results, so it can no longer be that benchmark. The census of what the attribution layer has read
(docs/ATTRIBUTION.md, "Pre-registration: one independent benchmark", 2026-09-28) lists the ENCODE
EPCrisprBenchmark training and held-out files (Nasser 2021, Gasperini 2019, Schraivogel 2020, Xie
2019, Morris 2023, Klann 2021, Reilly 2021, Ray 2025 DC-TAP, Guckelberger 2024), the ENCODE4
lentiMPRA (ENCSR106SZM), Kircher 2019 saturation mutagenesis, VISTA, GTEx v8 eQTLs with DAP-G
fine-mapping, the GWAS Catalog, ClinVar, the reporter-allele sets behind E1, the ENCODE cCRE
registry V3 and the ENCODE epigenome and 4DN Hi-C tracks used as features.

The benchmark chosen is outside all of it: the IGVF K562 CRISPRi Perturb-seq of the MHC locus
(Gersbach lab, Duke; analysis set IGVFDS7132YVKO, released 2026-06-08, CC BY 4.0). More than 500
non-promoter DHSs across the 3.5 Mb MHC were silenced with dCas9-KRAB and every gene of a targeted
panel read by single-cell RNA; the element-level FRACTEL table (IGVFFI4093WUVB) gives, for each
element and gene, a signed effect and a Benjamini-Hochberg FDR. It is a new study, in a new locus,
released after every table the project's features, thresholds and calibrations were built from.

What is scored is what the project already holds: the all-element AlphaGenome deletion sweep
(data/knowledge/alphagenome/all_elements, and its per-element response cache
data/knowledge/alphagenome/elements, every gene in the scorer's 1 Mb window on K562's own track).
No model request is made; a pair the sweep did not score is missing, never imputed.

`PREREGISTERED` was committed before any pair was scored. Only the table's shape was read before it
(columns, rows, elements, genes, chromosome, and the count of rows under FDR 0.05 in either
direction), plus the label-blind coverage of the pairs by the sweep.
"""

from __future__ import annotations

import bisect
import csv
import gzip
import math
import random
from collections import defaultdict
from pathlib import Path
from typing import Any

from genomeos.attribution import crispri

CACHE = Path("data/cache/indep")
SCREEN_FILE = "IGVFFI4093WUVB.tsv.gz"
SCE2G_FILE = "scE2G_K562_chr6_MHC.tsv.gz"  # distilled from IGVFFI1706PNVV: chr6 rows inside 27-36 Mb
REFERENCE = Path("data/reference")
CELL = "K562"
CHROM = "chr6"
FDR = 0.05
PROMOTER_WINDOW = 1_000  # an element within this of any GENCODE v50 gene TSS is a promoter, excluded
REACH = 500_000  # the scorer's half-window: a gene farther from the element is outside what the sweep saw
MIN_POSITIVES = 20  # fewer covered positives than this and the test is not read
MIN_COVERAGE = 0.5  # a smaller covered share of eligible pairs and the test is not read
BOOTSTRAPS = 2_000
SEED = 20260928

SOURCE: dict[str, Any] = {
    "accession": "IGVFFI4093WUVB",
    "analysis_set": "IGVFDS7132YVKO",
    "doi": "10.65695/IGVFDS7132YVKO",
    "url": "https://api.data.igvf.org/tabular-files/IGVFFI4093WUVB/@@download/IGVFFI4093WUVB.tsv.gz",
    "md5": "c6d0c44d5db44888743feb6b60c33a79",
    "released": "2026-06-08",
    "lab": "Charles Gersbach, Duke (IGVF)",
    "assay": "CRISPRi (dCas9-KRAB) Perturb-seq, targeted panel, FRACTEL element-level test",
    "cell": "K562",
    "assembly": "GRCh38",
    "annotation": "GENCODE 32",
    "coordinates": "BED-like, 0-based half-open (IGVF column specification)",
    "licence": "CC BY 4.0 (IGVF data, Registry of Open Data on AWS)",
}
PUBLISHED_MODEL: dict[str, Any] = {
    "accession": "IGVFFI1706PNVV",
    "prediction_set": "IGVFDS5428HHMB",
    "model": "scE2G v1.2 (Multiome), K562-CRISPRi, unfiltered element-gene scores, column `Score`",
    "url": "https://api.data.igvf.org/tabular-files/IGVFFI1706PNVV/@@download/IGVFFI1706PNVV.tsv.gz",
    "md5": "96232c244a50d7e522d8830cd0665aed",
    "released": "2025-06-25",
    "licence": "CC BY 4.0 (IGVF data)",
    "note": "trained on the ENCODE CRISPRi benchmark, which does not contain this screen",
}

PREREGISTERED = (
    "Registered 2026-09-28 by lane-indep, before any pair of IGVFFI4093WUVB was scored. Seen before "
    "registration: the column names, 34,279 rows, 581 elements, 59 genes, all on chr6, 385 rows at "
    "FDR < 0.05 in either direction, and the label-blind coverage of pairs by the sweep: of 20,228 "
    "eligible pairs 4,651 are covered (23%); within 500 kb of the TSS 4,472 of 6,544 (68%), beyond it "
    "almost none, because the sweep scored only the genes inside its 1 Mb window. The reach below was "
    "fixed from that label-blind count.\n"
    "\n"
    "PAIRS. Each row is one element (DHS) and one gene. Excluded before scoring, each counted: an "
    "element within 1 kb of any GENCODE v50 gene TSS on chr6 (a promoter or a TSS control); a gene "
    "GENCODE v50 does not carry on chr6 (matched by Ensembl id, then symbol); a pair whose gene and "
    "overlapping element also appear in either ENCODE EPCrisprBenchmark file (reused evidence). "
    "REACH: the primary population is the eligible pairs whose element midpoint lies within 500 kb of "
    "the gene's TSS, the half-width of the window the sweep scored; pairs beyond it are counted as out "
    "of the sweep's reach and not scored, since the sweep never saw those genes.\n"
    "\n"
    "LABELS (R2's outcome classes, kept apart). decrease: FDR < 0.05 and effect < 0 (the positive). "
    "increase: FDR < 0.05 and effect > 0, excluded from the primary test and reported alone, never "
    "read as 'no effect'. not detected: FDR >= 0.05 (the negative). The screen publishes no power "
    "column, so 'not detected' is not 'no effect', and the negatives carry that caveat.\n"
    "\n"
    "COVERAGE. A pair is covered when a sweep element overlaps the DHS and the per-element cache "
    "holds a K562 value for the pair's gene. Uncovered pairs are missing: counted, never imputed, "
    "never scored as zero. The test is not read (outcome 'not readable') if fewer than 20 covered "
    "positives remain or the covered share of in-reach decrease-or-not-detected pairs is below 0.5. "
    "Selection check, reported: the distance baseline's AP on the in-reach pairs the sweep did not "
    "cover against the pairs it did, so an easier covered subset is visible.\n"
    "\n"
    "PREDICTOR. The deletion: the largest predicted fall of the pair's gene on K562's own track over "
    "the overlapping sweep elements (crispri.deletion_drop, floored at zero).\n"
    "\n"
    "BASELINES. (1) distance: minus the log distance from the element's midpoint to the gene's GENCODE "
    "v50 TSS (floored at 1 kb). (2) the current annotation: 1 when the gene is an overlapping "
    "element's compiled target (the compact table's top predicted gene, coding or any), else 0; and "
    "the node rule, 1 when the gene is the nearest TSS in the element's CTCF node. (3) the published "
    "model: scE2G v1.2 in K562 (IGVFFI1706PNVV), the largest `Score` over its elements overlapping the "
    "DHS for the same Ensembl gene; pairs it does not score are missing for it and it is compared on "
    "the pairs both cover.\n"
    "\n"
    "METRIC. Average precision (AUPRC, ties as one threshold) on covered pairs, with a 95% interval "
    "from 2,000 bootstrap resamples of ELEMENTS (all pairs of a drawn element move together; seed "
    "20260928), and the paired difference on the same resamples. Resampling genes (59) is reported as "
    "a sensitivity reading. The base rate (positives over covered pairs) is reported beside every AP.\n"
    "\n"
    "PASS. The 95% interval of AP(deletion) - AP(distance) lies above 0.\n"
    "FALSIFIER. The deletion fails on independent evidence if that interval lies below 0, or if the "
    "upper bound of AP(deletion) does not exceed the base rate.\n"
    "Otherwise NOT ESTABLISHED: the interval straddles 0.\n"
    "Secondary, reported and not part of the verdict: the deletion against the annotation, the node "
    "rule and scE2G (paired, on shared pairs); the gene-resampled interval; the increases alone.\n"
    "\n"
    "READINGS. PASS: on one new K562 screen in one locus, the deletion ranks measured targets above "
    "distance; a statement about K562 and the MHC, not about other cells or loci. FAIL: the ranking "
    "that the ENCODE held-out file supported does not transfer to evidence the project never read, "
    "and every claim built on that file (target bands, the calibration) must be restated as "
    "unvalidated outside it. NOT ESTABLISHED: this screen cannot tell the deletion from distance, "
    "which is not support. NOT READABLE: the sweep does not cover the screen, and the lane's result is "
    "the coverage. In every case one locus bounds the claim: the element bootstrap does not sample "
    "locus-to-locus variation, and the MHC is unusual (HLA polymorphism, K562's low class I "
    "expression). AlphaGenome was trained on ENCODE K562 tracks, not on this screen; the objection "
    "that it has seen K562 chromatin stays open and is stated, not closed."
)


# ------------------------------------------------------------------------------------------
# Loading
# ------------------------------------------------------------------------------------------


def load_screen(path: Path | None = None) -> list[dict[str, Any]]:
    """The FRACTEL element-gene rows, typed; nothing is dropped here."""
    path = path or CACHE / SCREEN_FILE
    out = []
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out.append(
                {
                    "chrom": r["intended_target_chr"],
                    "start": int(r["intended_target_start"]),
                    "end": int(r["intended_target_end"]),
                    "dhs": r["dhs"],
                    "ensembl": r["ensembl_id"].split(".")[0],
                    "symbol": r["gene_symbol"],
                    "effect": float(r["FRACTEL_effect_size"]),
                    "fdr": float(r["FRACTEL_pval_fdr_corr"]),
                }
            )
    return out


def outcome(row: dict[str, Any]) -> str:
    """R2's classes: a significant increase is its own outcome, never 'no effect'."""
    if row["fdr"] < FDR:
        return "decrease" if row["effect"] < 0 else "increase"
    return "not_detected"


def gencode(chrom: str = CHROM, reference: Path = REFERENCE) -> dict[str, Any]:
    """GENCODE v50 genes on one chromosome: TSS by sweep name, Ensembl id to sweep name, every TSS."""
    from genomeos.genome.annotation import iter_gff3

    tss: dict[str, int] = {}
    by_ensembl: dict[str, str] = {}
    by_symbol: dict[str, str] = {}
    for row in iter_gff3(reference / f"gencode_v50_{chrom}.gff3.gz", types={"gene"}):
        ens = row.attrs.get("gene_id", "").split(".")[0]
        name = row.attrs.get("gene_name") or ens
        tss.setdefault(name, row.end - 1 if row.strand == "-" else row.start - 1)
        by_ensembl.setdefault(ens, name)
        if row.attrs.get("gene_name"):
            by_symbol.setdefault(row.attrs["gene_name"], name)
    return {"tss": tss, "by_ensembl": by_ensembl, "by_symbol": by_symbol, "all_tss": sorted(tss.values())}


def benchmark_pairs(knowledge: Path = crispri.KNOWLEDGE) -> dict[str, list[tuple[int, int]]]:
    """Every ENCODE EPCrisprBenchmark element (training and held-out) per chr6 gene, for the reuse rule."""
    out: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for name in (crispri.TRAINING, crispri.HELDOUT):
        with gzip.open(knowledge / name, "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["chrom"] != CHROM:
                    continue
                span = (int(r["chromStart"]), int(r["chromEnd"]))
                out[r["measuredGeneSymbol"]].append(span)
                ens = (r.get("measuredGeneEnsemblId") or "").split(".")[0]
                if ens:
                    out[ens].append(span)
    return out


def load_sce2g(path: Path | None = None) -> dict[str, list[tuple[int, int, float]]]:
    """scE2G's K562 scores near the MHC, per Ensembl gene: (start, end, Score)."""
    path = path or CACHE / SCE2G_FILE
    out: dict[str, list[tuple[int, int, float]]] = defaultdict(list)
    if not path.exists():
        return out
    with gzip.open(path, "rt") as fh:
        rows = (line for line in fh if not line.startswith("#"))
        for r in csv.DictReader(rows, delimiter="\t"):
            out[r["GeneEnsemblID"].split(".")[0]].append(
                (int(r["ElementStart"]), int(r["ElementEnd"]), float(r["Score"]))
            )
    return out


# ------------------------------------------------------------------------------------------
# Eligibility, coverage and features (label-blind)
# ------------------------------------------------------------------------------------------


def near_tss(all_tss: list[int], start: int, end: int, window: int = PROMOTER_WINDOW) -> bool:
    i = bisect.bisect_left(all_tss, start - window)
    return i < len(all_tss) and all_tss[i] < end + window


def prepare(
    rows: list[dict[str, Any]],
    genes: dict[str, Any],
    reused: dict[str, list[tuple[int, int]]],
    table: crispri.DeletionTable,
    cache: crispri.ElementCache,
    sce2g: dict[str, list[tuple[int, int, float]]] | None = None,
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    """Eligible pairs with their features; the exclusions counted. Reads no label."""
    excluded = {"promoter": 0, "gene_not_in_gencode_v50": 0, "reused_in_encode_benchmark": 0}
    out = []
    for r in rows:
        if near_tss(genes["all_tss"], r["start"], r["end"]):
            excluded["promoter"] += 1
            continue
        name = genes["by_ensembl"].get(r["ensembl"]) or genes["by_symbol"].get(r["symbol"])
        if name is None or name not in genes["tss"]:
            excluded["gene_not_in_gencode_v50"] += 1
            continue
        spans = [s for key in (name, r["ensembl"]) for s in reused.get(key, ())]
        if any(a < r["end"] and b > r["start"] for a, b in spans):
            excluded["reused_in_encode_benchmark"] += 1
            continue
        els = table.overlapping(r["chrom"], r["start"], r["end"])
        top, values = crispri.deletion_values(els, name, CELL, cache, r["chrom"])
        mid = (r["start"] + r["end"]) // 2
        raw = abs(mid - genes["tss"][name])
        d = max(crispri.MIN_DISTANCE, raw)
        sc = None
        if sce2g is not None:
            hits = [s for a, b, s in sce2g.get(r["ensembl"], ()) if a < r["end"] and b > r["start"]]
            sc = max(hits) if hits else None
        out.append(
            {
                **r,
                "gene": name,
                "element": r["dhs"],
                "covered": bool(els) and bool(values),
                "overlaps_sweep": bool(els),
                "deletion": crispri.deletion_drop(values) if values else None,
                "tss_distance": raw,
                "in_reach": raw <= REACH,
                "distance": -math.log(d),
                "annotation": top,
                "node": float(any((e.get("inferred") or {}).get("gene") == name for e in els)),
                "sce2g": sc,
            }
        )
    return out, excluded


# ------------------------------------------------------------------------------------------
# The metric
# ------------------------------------------------------------------------------------------


def average_precision(scores: list[float], labels: list[bool]) -> float | None:
    """AUPRC as a step sum over distinct thresholds (tied scores enter together)."""
    n_pos = sum(labels)
    if n_pos == 0:
        return None
    order = sorted(range(len(scores)), key=lambda i: -scores[i])
    ap, tp, fp, prev_recall, k = 0.0, 0, 0, 0.0, 0
    while k < len(order):
        s = scores[order[k]]
        while k < len(order) and scores[order[k]] == s:
            if labels[order[k]]:
                tp += 1
            else:
                fp += 1
            k += 1
        recall = tp / n_pos
        ap += (recall - prev_recall) * (tp / (tp + fp))
        prev_recall = recall
    return ap


def ap_of(pairs: list[dict[str, Any]], key: str) -> float | None:
    return average_precision([p[key] for p in pairs], [p["label"] == "decrease" for p in pairs])


def bootstrap(
    pairs: list[dict[str, Any]], keys: list[str], unit: str = "element", n: int = BOOTSTRAPS, seed: int = SEED
) -> dict[str, list[float]]:
    """AP per predictor, and keys[0] minus each other predictor, over resampled clusters of pairs."""
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for p in pairs:
        groups[p[unit]].append(p)
    ids = sorted(groups)
    rng = random.Random(seed)
    draws: dict[str, list[float]] = {k: [] for k in keys}
    for k in keys[1:]:
        draws[f"{keys[0]}-{k}"] = []
    for _ in range(n):
        sample = [p for i in (rng.choice(ids) for _ in ids) for p in groups[i]]
        aps = {k: ap_of(sample, k) for k in keys}
        if any(v is None for v in aps.values()):
            continue
        for k in keys:
            draws[k].append(aps[k])
        for k in keys[1:]:
            draws[f"{keys[0]}-{k}"].append(aps[keys[0]] - aps[k])
    return draws


def interval(xs: list[float]) -> list[float] | None:
    if not xs:
        return None
    s = sorted(xs)
    return [round(s[int(0.025 * (len(s) - 1))], 4), round(s[int(0.975 * (len(s) - 1))], 4)]


def verdict(delta_ci: list[float] | None, ap_ci: list[float] | None, base: float, readable: bool) -> str:
    """The registered reading: PASS, FAIL, NOT ESTABLISHED or NOT READABLE."""
    if not readable or delta_ci is None or ap_ci is None:
        return "NOT READABLE"
    if delta_ci[1] < 0 or ap_ci[1] <= base:
        return "FAIL"
    if delta_ci[0] > 0:
        return "PASS"
    return "NOT ESTABLISHED"
