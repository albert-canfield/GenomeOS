# SPDX-License-Identifier: AGPL-3.0-or-later
"""M5b step 1: a fine-mapped eQTL arm calibrated on sequence outside the 882 real-unknown blocks, sized.

Writes `data/results/eqtl_arm_calibration.json`. The order is the method (genomeos/attribution/eqtl_arm.py):
the blocks are read and gated first, the exclusion is built from them and checked against every block,
and only then is the GTEx DAP-G track opened, behind that exclusion. The CRISPRi benchmark's measured
pairs, their matched windows and their flanks are drawn outside it; the track is read only there.

Inputs, all public and held locally (git-ignored): the GTEx v8 DAP-G bigBed and the dbSNP 155 Common
bigBed, each fetched whole from UCSC into data/knowledge/eqtl_arm; the CRISPRi benchmark tables in
data/knowledge/crispri; the element archive and the organiser inputs the population is built from;
GENCODE v50. 0 AlphaGenome requests, no money.
"""

from __future__ import annotations

import argparse
import json
import random
import statistics
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri, organise  # noqa: E402
from genomeos.attribution import eqtl_arm as ea  # noqa: E402
from genomeos.attribution import measured_unknown as mu  # noqa: E402
from genomeos.results import save_result  # noqa: E402

NAME = "eqtl_arm_calibration"
ENTRY = "scripts/eqtl_arm_calibrate.py"
OWN_CODE = (ENTRY, "genomeos/attribution/eqtl_arm.py")
SCALES = ("element", "block")


def say(*a: Any) -> None:
    print(*a, flush=True)


def gencode(chrom: str) -> tuple[dict[str, str], dict[str, list[tuple[str, int]]]]:
    """Unversioned Ensembl id -> symbol, and symbol -> [(chrom, TSS)], from the project's GENCODE."""
    from genomeos.genome import Annotation
    from genomeos.genome.annotation import default_gencode

    path = default_gencode({chrom})
    if path is None:
        return {}, {}
    ann = Annotation.from_gff3(path, {chrom})
    sym: dict[str, str] = {}
    tss: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for gid, g in ann.genes.items():
        sym[gid.split(".")[0]] = g.symbol
        if g.locus.chrom != chrom:
            continue
        plus = str(getattr(g.locus.strand, "value", g.locus.strand)) == "+"
        tss[g.symbol].append((chrom, g.locus.start if plus else g.locus.end - 1))
    return sym, tss


def exclusion_check(rows: list[dict[str, Any]], ex: ea.Exclusion) -> dict[str, Any]:
    """Every block's first and last base, and its midpoint, lie inside the exclusion; and the bases a
    margin outside it lie outside. Refuses to continue otherwise."""
    inside = all(
        ex.holds(r["chrom"], r["start"])
        and ex.holds(r["chrom"], r["end"] - 1)
        and ex.holds(r["chrom"], (r["start"] + r["end"]) // 2)
        and ex.hits(r["chrom"], r["start"], r["end"])
        for r in rows
    )
    widened = all(
        ex.holds(r["chrom"], max(0, r["start"] - ex.margin))
        and ex.holds(r["chrom"], r["end"] + ex.margin - 1)
        for r in rows
    )
    if not (inside and widened):
        raise SystemExit("the exclusion does not hold every block: nothing is read")
    return {
        "blocks": ex.blocks,
        "margin_each_side": ex.margin,
        "merged_spans": sum(len(v) for v in ex.spans.values()),
        "bases": ex.bases(),
        "every_block_inside": inside,
        "every_widened_edge_inside": widened,
        "digest": ex.digest(),
    }


def run(experiments: int, bootstraps: int) -> dict[str, Any]:
    t0 = time.time()
    rows: list[dict[str, Any]] = []
    attribution: list[dict[str, Any]] = []
    for c in mu.CHROMS:
        blocks = organise.blocks(c)
        els = mu._elements(c)
        rows += mu.population_of(c, blocks, els)
        attribution += ea.block_attribution(c, blocks, els, mu.MIN_LOG2)
    gate = mu.gate(rows)
    say("population", gate)
    if not gate["ok"]:
        raise SystemExit("the population is not the clause-2 population; nothing is read")
    ex = ea.Exclusion(rows)
    excl = exclusion_check(rows, ex)
    say("exclusion", excl)
    block_length = int(round(statistics.median(r["length"] for r in rows)))
    lengths = {"element": None, "block": block_length}

    pairs = ea.element_pairs(crispri.load(crispri.TRAINING) + crispri.load(crispri.HELDOUT))
    say("benchmark pairs", len(pairs), "positive", sum(p["regulated"] for p in pairs))
    units: dict[str, list[ea.Unit]] = {}
    dropped: dict[str, dict[str, int]] = {}
    for sc in SCALES:
        units[sc], dropped[sc] = ea.build_units(pairs, ex, lengths[sc], sc)
        say(sc, "units", len(units[sc]), "dropped", dropped[sc])
    avoid = ea.positive_spans(pairs)

    track = ea.DapgTrack(ex)
    common = ea.CommonVariants()
    tss_of: dict[str, list[tuple[str, int]]] = defaultdict(list)
    scored: dict[str, list[dict[str, Any]]] = {sc: [] for sc in SCALES}
    unmatched: dict[str, dict[str, int]] = {sc: {"positive": 0, "negative": 0} for sc in SCALES}
    far: dict[str, dict[str, list[int]]] = {sc: defaultdict(list) for sc in SCALES}
    track_rows = 0
    for c in mu.CHROMS:
        sym, tss = gencode(c)
        for k, v in tss.items():
            tss_of[k] += v
        fine = track.chromosome(c)
        track_rows += len(fine.pos)
        codes = ea.gene_codes(fine, sym)
        size = track.bb.chroms[c][1] if c in track.bb.chroms else 0
        for sc in SCALES:
            for u in units[sc]:
                if u.chrom != c:
                    continue
                u.common = common.count(c, u.start, u.end)
                rng = random.Random(f"{ea.SEED}:{u.key}")
                got = ea.draw_windows(
                    u, ea.WINDOWS + 1, common.count, ex, avoid.get((c, u.gene), []), size, rng
                )
                if got is None:
                    unmatched[sc]["positive" if u.regulated else "negative"] += 1
                    far[sc]["unmatched" if u.regulated else "unmatched_negative"].append(u.distance)
                    continue
                far[sc]["matched" if u.regulated else "matched_negative"].append(u.distance)
                u.windows = got
                u.flanks = ea.flanks_of(u, ex, size)
                scored[sc].append(ea.scored_unit(u, fine, codes.get(u.gene, set())))
        say(c, "track rows outside the exclusion", len(fine.pos), {sc: len(v) for sc, v in scored.items()})

    avail = ea.available(attribution, tss_of)
    say("available", {k: v for k, v in avail.items() if k != "distance_quartiles"})
    rng = np.random.default_rng(ea.SEED)
    out_scales: dict[str, Any] = {}
    for sc in SCALES:
        pos = [r for r in scored[sc] if r["regulated"]]
        neg = [r for r in scored[sc] if not r["regulated"]]
        n_avail = avail[f"usable_at_the_{sc}_scale"]
        s_pos, s_neg = ea.summarise(pos, rng, bootstraps), ea.summarise(neg, rng, bootstraps)
        say(sc, "positives", {k: s_pos.get(k) for k in ("units", "rate_at_the_units", "rate_at_the_windows")})
        say(sc, "negatives", {k: s_neg.get(k) for k in ("units", "rate_at_the_units", "rate_at_the_windows")})
        sizes = ea.sizing(pos, neg, n_avail, rng, experiments=experiments)
        decision = ea.fair_test(sizes, n_avail)
        say(sc, "decision", decision)
        out_scales[sc] = {
            "unit": "the CRISPRi element with 500 bases each side"
            if sc == "element"
            else f"a {block_length}-base window (the 882 blocks' median length) centred on the element",
            "length": lengths[sc],
            "dropped_before_drawing": dropped[sc],
            "unmatched_after_drawing": unmatched[sc],
            "median_distance_to_the_tss": {
                k: int(statistics.median(v)) for k, v in sorted(far[sc].items()) if v
            },
            "positives": s_pos,
            "tested_negatives": s_neg,
            "sizing": sizes,
            "fair_test": decision,
        }
    return {
        "population": {"gate": gate, "blocks_naming_any_gene": sum(1 for r in rows if r["named_genes"])},
        "exclusion": excl,
        "benchmark": {
            "pairs": len(pairs),
            "positive": sum(p["regulated"] for p in pairs),
            "tested_negative": sum(not p["regulated"] for p in pairs),
            "unit": "one (element, measured gene) pair; a pair tested in several screens counts once and "
            "is positive if any screen called it",
        },
        "track_rows_decoded_outside_the_exclusion": track_rows,
        "available": avail,
        "scales": out_scales,
        "seconds": round(time.time() - t0, 1),
    }


READ_GROUPS = (
    ("data/knowledge/eqtl_arm/", "GTEx v8 DAP-G and dbSNP 155 Common, fetched whole from UCSC"),
    ("data/knowledge/crispri/", "the CRISPRi benchmark (EPCrisprBenchmark), as held"),
    ("data/knowledge/alphagenome/all_elements/", "the element archive the population is built from"),
    ("data/results/", "the organiser's per-chromosome inputs, as committed"),
    ("data/reference/", "GENCODE v50 per chromosome"),
    ("data/", "other held inputs the organiser reads"),
)


def result_manifest(reads: list[str], experiments: int, bootstraps: int) -> dict[str, Any]:
    groups: dict[str, list[str]] = {}
    for path in sorted(set(reads)):
        label = next(lab for pre, lab in READ_GROUPS if path.startswith(pre))
        groups.setdefault(label, []).append(path)
    inputs = [mf.files_entry(label, paths, partition="read by the run") for label, paths in groups.items()]
    return {
        "sources": [
            {
                "accession": "UCSC gtexDapg (GTEx v8 cis-eQTL fine-mapping, DAP-G, 49 tissues)",
                "version": ea.DAPG_URL,
            },
            {"accession": "UCSC dbSnp155Common (dbSNP build 155, common variants)", "version": ea.COMMON_URL},
            {
                "accession": "EPCrisprBenchmark (ENCODE-rE2G) training K562 and held-out 5 cell types",
                "version": crispri.BASE_URL,
            },
            {"accession": "GENCODE v50", "version": "data/reference/gencode_v50_chr*.gff3.gz as held"},
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "pip": ea.PIP,
            "margin": ea.MARGIN,
            "exclusion_margin": ea.EXCLUSION_MARGIN,
            "cis_window": ea.CIS_WINDOW,
            "tss_guard": ea.TSS_GUARD,
            "windows_per_unit": ea.WINDOWS,
            "pseudo_block_windows": 1,
            "distance_band": ea.DISTANCE_BAND,
            "min_band_in_unit_lengths": ea.MIN_BAND_LENGTHS,
            "count_band": ea.COUNT_BAND,
            "draw_tries_per_window": ea.DRAW_TRIES,
            "seed": ea.SEED,
            "tissue_of_cell": ea.TISSUE_OF_CELL,
            "alpha_one_sided": ea.ALPHA,
            "power": ea.POWER,
            "shares_right": list(ea.SHARES),
            "n_grid": list(ea.N_GRID),
            "experiments": experiments,
            "bootstraps": bootstraps,
            "fooled_bound": ea.FOOLED_BOUND,
            "min_log2": mu.MIN_LOG2,
        },
        "exclusions": [
            "every DAP-G record whose variant lies inside the 882 blocks widened by 1 kb: dropped on its "
            "position before any other field is split; not counted",
            "every calibration unit, window and flank touching the widened blocks",
            "units beyond the 1 Mb cis window or holding the gene's TSS with its 1 kb guard",
            "tested negatives overlapping a measured positive element of the same gene",
            "units for which the draw found fewer than 5 matched windows (counted in the result)",
        ],
        "partitions": "n/a: a calibration on measured pairs outside the blocks; nothing is trained",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument("--no-save", action="store_true")
    ap.add_argument("--experiments", type=int, default=ea.EXPERIMENTS)
    ap.add_argument("--bootstraps", type=int, default=ea.BOOTSTRAPS)
    args = ap.parse_args(argv)
    mf.trace_begin()
    body = run(args.experiments, args.bootstraps)
    record: dict[str, Any] = {
        "result": NAME,
        "date": time.strftime("%Y-%m-%d"),
        "question": "on sequence outside the 882 real-unknown blocks, how often does a GTEx v8 fine-mapped "
        "cis-eQTL name the gene a CRISPRi screen measured an element to regulate, how often at matched "
        "windows, how clustered are its hits, and how many blocks would a test of the blocks' attribution "
        "need?",
        "what_an_eqtl_is": "an association between a variant and a gene's expression across people; not a "
        "perturbation and not mechanism. Fine-mapping localises the association to a variant with a "
        "posterior; it does not say which element acts",
        "fair_test_rule": ea.FAIR_TEST_RULE,
        "window_design": "four control windows and one pseudo-block per unit, within a quarter of the unit's "
        "distance to the TSS (at least three unit lengths) and a quarter of its common-variant count. A "
        "design probe that read coordinates and common-variant counts only, before any eQTL value, found "
        "that ten disjoint windows in a tenth-wide band fit around under half of the positives at the "
        "element scale and about a tenth at the block scale, the losses near the TSS",
        "eqtl_values_read_inside_the_blocks": 0,
        **body,
        "alphagenome_requests": 0,
        "money": "none",
    }
    record["result_manifest"] = result_manifest(mf.traced_reads(), args.experiments, args.bootstraps)
    if args.no_save:
        json.dump({k: record[k] for k in ("available", "scales")}, sys.stdout, indent=1, default=str)
        return 0
    path = save_result(NAME, record)
    say("wrote", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
