# SPDX-License-Identifier: AGPL-3.0-or-later
"""The node-containment claim audited, then held against measured CRISPRi pairs.

    uv run python scripts/node_containment_audit.py [--workers 6]

The claim under audit, quoted seven times in the docs: the default `ctcf_only` node caller keeps an
element's most-moved coding gene inside the element's own node for 75.4% of 440,377 archive
elements, against 72.5% for "as many boundaries placed at random" - +2.9 points.

Three things the claim has never carried, and this script produces all three:

1. **The baseline, defined.** `scripts/oriented_domains.py` draws its random boundary set as
   `len(edges)` uniform positions on the chromosome, with no other matching, and - the part that
   matters - without the 50 kb minimum-size merge that every real node passes through
   (`genome/domains.py _domains_from`). Real nodes are therefore floored at 50 kb and random nodes
   are not, which lowers random containment by construction. This script scores the published
   control and three stricter ones on the same pairs:

   - `uniform`        the published control, reproduced (as many uniform positions, no merge)
   - `uniform_merged` the same positions passed through the same 50 kb merge, so the node length
                      floor is matched and only the placement differs
   - `circular`       the real boundary set rotated by a random offset modulo the chromosome, which
                      preserves the inter-boundary spacing distribution exactly - hence the node
                      count and the node length distribution - and destroys only the relation
                      between a boundary and the genes near it
   - `count_matched`  uniform positions drawn so that the post-merge edge count equals the real one

2. **An interval.** The 440,377 pairs are not independent: elements in one node share one boundary
   pair, so the containment outcome of every element in a node is decided by the same two numbers,
   and many elements share a target gene. Two resamplings are reported for that reason:

   - **over chromosomes** (the headline): 24 units resampled with replacement, pooled shares
     recomputed inside each resample. This is the interval to quote, because the chromosome is the
     largest unit over which boundary placement is plausibly independent. Its weakness is stated in
     the output: 24 units, and chrY contributes 0.1% of the pairs.
   - **over nodes**: elements resampled in node-sized clusters, which keeps the element count and
     respects the dependence the pooled share actually has.

   A paired per-chromosome sign test is reported beside them, because "ahead on 18 of 24
   chromosomes" is the only dispersion the claim has ever been quoted with and it appears in no
   result file.

3. **The denominator, identified.** `coding_named` is counted here field by field so that the three
   published numbers - 961,227 scored, 593,765 banded, 440,377 coding-target - can be told apart.

Then stage 2: the same containment statistic on a measured element-gene set. The ENCODE CRISPRi
benchmark's `Regulated=TRUE` pairs (data/knowledge/crispri) give (element midpoint, measured gene
TSS) in place of (element midpoint, TSS of the gene the deletion moves most). The registration for
that arm, with its power arithmetic, is in docs/NODES-READER-WRITER.md and was committed before this
script was run.

No model is called and nothing is fetched.
"""

from __future__ import annotations

import bisect
import json
import random
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri  # noqa: E402
from genomeos.attribution.measured import CRISPRI_SPLIT_OF  # noqa: E402
from genomeos.genome import Annotation, default_gencode  # noqa: E402
from genomeos.genome.domains import MIN_DOMAIN, infer_domains  # noqa: E402
from genomeos.genome.regulatory import load_ccres  # noqa: E402
from genomeos.results import save_result  # noqa: E402

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
ARCHIVE = Path("data/knowledge/alphagenome/all_elements")
SHUFFLES = 20  # the draw count scripts/oriented_domains.py uses, kept so the control reproduces
SEED = 7  # the same per-chromosome seed, so `uniform` reproduces the published 72.5%
BOOTSTRAPS = 2000
CONTROLS = ("uniform", "uniform_merged", "circular", "count_matched")
CRISPRI_TABLES = (
    "EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz",
    "EPCrisprBenchmark_combined_data.heldout_5_cell_types.GRCh38.tsv.gz",
)


def chrom_length(chrom: str) -> int | None:
    """From the faidx beside the reference, so no bgzip index is opened for a number."""
    for p in (Path(f"data/reference/{chrom}.fa.gz.fai"), Path(f"data/reference/{chrom}.fa.fai")):
        if p.exists():
            for line in p.read_text().splitlines():
                f = line.split("\t")
                if f[0] == chrom:
                    return int(f[1])
    return None


def _tss(g) -> int:
    return g.locus.end - 1 if g.locus.strand.value == "-" else g.locus.start


def _inside(starts: list[int], pairs: list[tuple[int, int]]) -> int:
    """How many (position, TSS) pairs fall in one node, nodes starting at `starts`."""
    return sum(1 for m, t in pairs if bisect.bisect_right(starts, m) == bisect.bisect_right(starts, t))


def _merge_starts(boundaries: list[int], length: int, min_size: int = MIN_DOMAIN) -> list[int]:
    """The node starts `genome/domains.py _domains_from` produces from these boundaries.

    Copied in form, not called, because `_domains_from` also needs cCREs and an annotation to
    assign contents, and only the edges matter here. Pinned against the real caller in
    tests/test_node_containment_audit.py."""
    bounds = [0, *sorted(boundaries), length]
    edges = [bounds[0]]
    for b in bounds[1:]:
        if b - edges[-1] < min_size and len(edges) > 1:
            continue
        edges.append(b)
    if edges[-1] != length:
        edges[-1] = length
    return edges[:-1]  # the node starts: every edge but the terminal one


def controls(
    real_starts: list[int], n_edges: int, length: int, pairs: list[tuple[int, int]], rng: random.Random
) -> dict[str, dict]:
    """Every control's per-draw inside counts, and the node geometry each draw produced."""
    out: dict[str, dict] = {}
    real_boundaries = real_starts[1:]

    def record(name: str, draws: list[tuple[list[int], int]]) -> None:
        ins = [_inside(st, pairs) for st, _ in draws]
        out[name] = {
            "draws": ins,
            "mean_inside": sum(ins) / len(ins),
            "mean_nodes": sum(len(st) for st, _ in draws) / len(draws),
            "mean_median_node_length": sum(m for _, m in draws) / len(draws),
        }

    def geom(starts: list[int]) -> tuple[list[int], int]:
        lens = sorted(b - a for a, b in zip(starts, [*starts[1:], length], strict=False))
        return starts, lens[len(lens) // 2] if lens else 0

    # 1. the published control: as many uniform positions, no minimum-size merge
    record(
        "uniform",
        [geom([0, *sorted(rng.randint(1, length - 1) for _ in range(n_edges))]) for _ in range(SHUFFLES)],
    )
    # 2. the same positions through the same 50 kb merge: the node length floor matched
    record(
        "uniform_merged",
        [
            geom(_merge_starts([rng.randint(1, length - 1) for _ in range(n_edges)], length))
            for _ in range(SHUFFLES)
        ],
    )
    # 3. the real spacing distribution, rotated: node count and node lengths matched exactly
    rot = []
    for _ in range(SHUFFLES):
        off = rng.randint(1, length - 1)
        rot.append(geom([0, *sorted(((b + off) % length) or 1 for b in real_boundaries)]))
    record("circular", rot)
    # 4. uniform positions drawn until the post-merge edge count matches the real one
    cm = []
    target = len(real_starts)
    cap = max(n_edges, 4 * (length // MIN_DOMAIN))  # the merge saturates near length/MIN_DOMAIN edges
    for _ in range(SHUFFLES):
        k = n_edges
        st = _merge_starts([rng.randint(1, length - 1) for _ in range(k)], length)
        while len(st) < target and k < cap:
            k = min(cap, int(k * 1.4) + 1)
            st = _merge_starts([rng.randint(1, length - 1) for _ in range(k)], length)
        if len(st) > target:  # the merge overshoots; trim at random so the node count is exact
            st = [0, *sorted(rng.sample(st[1:], target - 1))]
        cm.append(geom(st))
    record("count_matched", cm)
    return out


def chrom_row(chrom: str) -> dict:
    length = chrom_length(chrom)
    gff = default_gencode({chrom})
    ccres = load_ccres(chrom)
    if not length or not gff or not ccres:
        return {"chrom": chrom, "skipped": True}
    ann = Annotation.from_gff3(gff, {chrom})
    doms = infer_domains(chrom, length, ccres, ann)
    starts = [d.start for d in doms]
    n_edges = len(doms) - 1
    coding = {
        x.symbol: x for x in ann.genes.values() if x.locus.chrom == chrom and x.type == "protein_coding"
    }
    anyg = {x.symbol: x for x in ann.genes.values() if x.locus.chrom == chrom}

    # the denominator, counted field by field
    rows = json.loads((ARCHIVE / f"{chrom}.json").read_text()) if (ARCHIVE / f"{chrom}.json").exists() else []
    census = {
        "archive_rows": len(rows),
        "has_predicted": 0,
        "has_predicted_coding": 0,
        "predicted_coding_gene_in_annotation": 0,
        "predicted_coding_gene_protein_coding": 0,
        "has_confidence": 0,
        "predicted_coding_has_confidence": 0,
        # R4: a link written since states no confidence and carries a certainty record instead
        "has_certainty": 0,
        "predicted_coding_has_certainty": 0,
    }
    pairs: list[tuple[int, int]] = []
    pairs_any: list[tuple[int, int]] = []
    for e in rows:
        mid = (e["start"] + e["end"]) // 2
        p = e.get("predicted") or {}
        pc = e.get("predicted_coding") or {}
        if p.get("gene"):
            census["has_predicted"] += 1
        if p.get("confidence") is not None:
            census["has_confidence"] += 1
        if p.get("certainty") is not None:
            census["has_certainty"] += 1
        if pc.get("gene"):
            census["has_predicted_coding"] += 1
            if pc.get("confidence") is not None:
                census["predicted_coding_has_confidence"] += 1
            if pc.get("certainty") is not None:
                census["predicted_coding_has_certainty"] += 1
            if pc["gene"] in anyg:
                census["predicted_coding_gene_in_annotation"] += 1
            if pc["gene"] in coding:
                census["predicted_coding_gene_protein_coding"] += 1
                pairs.append((mid, _tss(coding[pc["gene"]])))
        if p.get("gene") in anyg:
            pairs_any.append((mid, _tss(anyg[p["gene"]])))

    lens = sorted(d.length for d in doms)
    row: dict = {
        "chrom": chrom,
        "length": length,
        "nodes": len(doms),
        "edges": n_edges,
        "median_node_length": lens[len(lens) // 2],
        "min_node_length": lens[0],
        "census": census,
        "coding_named": len(pairs),
        "coding_inside": _inside(starts, pairs),
        "any_named": len(pairs_any),
        "any_inside": _inside(starts, pairs_any),
    }
    if pairs:
        rng = random.Random(SEED)
        row["controls"] = controls(starts, n_edges, length, pairs, rng)
    # the node each element's containment outcome belongs to, for the clustered resample
    clusters: dict[int, list[int]] = {}
    for m, t in pairs:
        i = bisect.bisect_right(starts, m)
        clusters.setdefault(i, []).append(1 if i == bisect.bisect_right(starts, t) else 0)
    row["node_clusters"] = [len(v) for v in clusters.values()]
    row["node_cluster_inside"] = [sum(v) for v in clusters.values()]
    return row


# ---------------------------------------------------------------- intervals


def pooled(rows: list[dict], control: str) -> tuple[float, float, float]:
    named = sum(r["coding_named"] for r in rows)
    inside = sum(r["coding_inside"] for r in rows)
    ctrl = sum(r["controls"][control]["mean_inside"] for r in rows if "controls" in r)
    return inside / named, ctrl / named, inside / named - ctrl / named


def bootstrap_chromosomes(rows: list[dict], control: str, rng: random.Random) -> dict:
    """Resample the 24 chromosomes with replacement; recompute the pooled excess in each resample."""
    ex = []
    for _ in range(BOOTSTRAPS):
        s = [rows[rng.randrange(len(rows))] for _ in rows]
        named = sum(r["coding_named"] for r in s)
        if not named:
            continue
        inside = sum(r["coding_inside"] for r in s)
        ctrl = sum(r["controls"][control]["mean_inside"] for r in s if "controls" in r)
        ex.append((inside - ctrl) / named)
    ex.sort()
    lo, hi = ex[int(0.025 * len(ex))], ex[int(0.975 * len(ex)) - 1]
    return {
        "point": round(pooled(rows, control)[2] * 100, 3),
        "ci95": [round(lo * 100, 3), round(hi * 100, 3)],
        "crosses_zero": lo <= 0 <= hi,
        "resamples": len(ex),
        "unit": "chromosome",
    }


def bootstrap_nodes(rows: list[dict], control: str, rng: random.Random) -> dict:
    """Resample elements in node-sized clusters: the element count is kept and so is the dependence
    (every element in a node shares one boundary pair, so its outcome is not independent)."""
    sizes = [n for r in rows for n in r["node_clusters"]]
    ins = [n for r in rows for n in r["node_cluster_inside"]]
    total = sum(sizes)
    base_ctrl = sum(r["controls"][control]["mean_inside"] for r in rows if "controls" in r) / total
    ex = []
    for _ in range(BOOTSTRAPS):
        n = i = 0
        while n < total:
            j = rng.randrange(len(sizes))
            n += sizes[j]
            i += ins[j]
        ex.append(i / n - base_ctrl)
    ex.sort()
    lo, hi = ex[int(0.025 * len(ex))], ex[int(0.975 * len(ex)) - 1]
    return {
        "point": round(pooled(rows, control)[2] * 100, 3),
        "ci95": [round(lo * 100, 3), round(hi * 100, 3)],
        "crosses_zero": lo <= 0 <= hi,
        "note": "the control share is held fixed; only the measured share is resampled",
        "unit": "node cluster",
    }


def sign_test(rows: list[dict], control: str) -> dict:
    """Ahead on how many chromosomes, and the exact two-sided binomial p at 0.5."""
    from math import comb

    per = []
    for r in rows:
        if "controls" not in r:
            continue
        m = r["coding_inside"] / r["coding_named"]
        c = r["controls"][control]["mean_inside"] / r["coding_named"]
        per.append(
            {"chrom": r["chrom"], "excess_points": round((m - c) * 100, 2), "pairs": r["coding_named"]}
        )
    ahead = sum(1 for p in per if p["excess_points"] > 0)
    n = len(per)
    p_one = sum(comb(n, k) for k in range(ahead, n + 1)) / 2**n
    return {
        "ahead": ahead,
        "of": n,
        "p_one_sided": round(p_one, 4),
        "p_two_sided": round(min(1.0, 2 * p_one), 4),
        "per_chromosome": sorted(per, key=lambda x: x["excess_points"]),
    }


# ---------------------------------------------------------------- stage 2


def measured_arm(rows: list[dict]) -> dict:
    """The same containment statistic on CRISPRi `Regulated=TRUE` pairs.

    (element midpoint, measured gene TSS) replaces (element midpoint, TSS of the gene the deletion
    moves most). Everything downstream - the node sets, the four controls, the resampling - is the
    code above, unchanged."""
    pairs = []
    for t in CRISPRI_TABLES:
        if (crispri.KNOWLEDGE / t).exists():
            pairs += crispri.load(t)
    reg = [p for p in pairs if p.regulated and p.tss is not None]
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    seen: set[tuple] = set()
    distinct: list = []
    for p in reg:
        by_chrom.setdefault(p.chrom, []).append((p.midpoint, p.tss))
        if p.element not in seen:
            seen.add(p.element)
            distinct.append(p)
    cells: dict[str, int] = {}
    for p in reg:
        cells[p.cell] = cells.get(p.cell, 0) + 1
    dists = sorted(abs(p.distance) for p in reg)
    inventory = {
        "regulated_pairs": len(reg),
        "distinct_element_cell": len(seen),
        "cells": cells,
        "median_distance_to_tss": int(dists[len(dists) // 2]) if dists else None,
        "beyond_10kb": sum(1 for d in dists if d > 10_000),
        "chromosomes": sorted(by_chrom),
    }
    out: dict = {"inventory": inventory, "rows": []}
    for r in rows:
        mp = by_chrom.get(r["chrom"])
        if not mp:
            continue
        length, chrom = r["length"], r["chrom"]
        ccres = load_ccres(chrom)
        gff = default_gencode({chrom})
        ann = Annotation.from_gff3(gff, {chrom})
        doms = infer_domains(chrom, length, ccres, ann)
        starts = [d.start for d in doms]
        rng = random.Random(SEED)
        cs = controls(starts, len(doms) - 1, length, mp, rng)
        clusters: dict[int, list[int]] = {}
        for m, t in mp:
            i = bisect.bisect_right(starts, m)
            clusters.setdefault(i, []).append(1 if i == bisect.bisect_right(starts, t) else 0)
        out["rows"].append(
            {
                "chrom": chrom,
                "length": length,
                "coding_named": len(mp),
                "coding_inside": _inside(starts, mp),
                "controls": cs,
                "node_clusters": [len(v) for v in clusters.values()],
                "node_cluster_inside": [sum(v) for v in clusters.values()],
            }
        )
    return out


def readings(rows: list[dict], label: str) -> dict:
    rng = random.Random(11)
    out: dict = {"label": label, "pairs": sum(r["coding_named"] for r in rows), "controls": {}}
    out["measured_share"] = round(sum(r["coding_inside"] for r in rows) / max(1, out["pairs"]), 4)
    for c in CONTROLS:
        share, ctrl, ex = pooled(rows, c)
        out["controls"][c] = {
            "control_share": round(ctrl, 4),
            "excess_points": round(ex * 100, 3),
            "mean_median_control_node_length": int(
                sum(r["controls"][c]["mean_median_node_length"] for r in rows if "controls" in r)
                / max(1, sum(1 for r in rows if "controls" in r))
            ),
            "mean_control_nodes": int(sum(r["controls"][c]["mean_nodes"] for r in rows if "controls" in r)),
            "bootstrap_chromosomes": bootstrap_chromosomes(rows, c, rng),
            "bootstrap_nodes": bootstrap_nodes(rows, c, rng),
            "sign_test": sign_test(rows, c),
        }
    return out


@mf.depends_on_models("alphagenome")  # the sweep's model, as far as the disk says (R9)
def manifest(stage: int, rows: list[dict]) -> dict:
    """The provenance contract (review item R9). Stage 1 reads the model's archive; stage 2 reads the
    CRISPRi benchmark in its place. Both read the same node caller over the same annotation."""
    chroms = [r["chrom"] for r in rows]
    inputs = []
    for c in chroms:
        for p in (Path(f"data/reference/{c}.fa.gz.fai"), Path(f"data/reference/{c}.fa.fai")):
            if p.exists():
                inputs.append(mf.input_entry(p, partition=None))
                break
        inputs.append(mf.input_entry(default_gencode({c}), partition=None))
        inputs.append(mf.input_entry(Path(f"data/results/ccres_{c}.bed.gz"), partition=None))
        if (ARCHIVE / f"{c}.json").exists():  # stage 2 reads it too, for the rows it reuses
            inputs.append(mf.input_entry(ARCHIVE / f"{c}.json", partition=None))
    sources = [
        {
            "accession": "UCSC hg38 chromosome FASTA (lengths from its faidx)",
            "version": "GRCh38, as served at hgdownload goldenPath/hg38/chromosomes; pinned here by sha256",
        },
        {"accession": "GENCODE human gene annotation", "version": "release 50 (GRCh38.p14)"},
        {
            "accession": "ENCODE SCREEN registry of cCREs (CTCF-only elements as node boundaries)",
            "version": "V3, downloads.wenglab.org/V3/GRCh38-cCREs.bed",
        },
        {
            "accession": "ENCODE SCREEN cCREs scored by AlphaGenome deletion (all-element archive)",
            "version": "AlphaGenome as served during the 2026-09 all-element sweep (unpinned); "
            "pinned here by sha256",
            "path": str(ARCHIVE),
        },
    ]
    if stage == 2:
        sources.append(
            {
                "accession": "EngreitzLab/CRISPR_comparison resources/crispr_data (Gschwind et al. 2025)",
                "version": "main branch, unpinned upstream; fetched 2026-09-16; pinned here by sha256",
                "url": crispri.BASE_URL,
            }
        )
        for t in CRISPRI_TABLES:
            if (crispri.KNOWLEDGE / t).exists():
                inputs.append(mf.input_entry(crispri.KNOWLEDGE / t, partition=CRISPRI_SPLIT_OF.get(t)))
    exclusions: list = [
        {"chromosomes_skipped": [c for c in CHROMS if c not in chroms], "why": "no length, GENCODE or cCREs"},
    ]
    if stage == 1:
        exclusions.append(
            "archive elements with no predicted coding gene, or whose gene is not protein_coding in "
            "GENCODE 50, are not pairs (each step counted in denominator_census)"
        )
    else:
        exclusions.append("CRISPRi pairs not Regulated=TRUE, or whose gene has no TSS in the table")
    return {
        "sources": sources,
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "stage": stage,
            "node_caller": "genome.domains.infer_domains, the default ctcf_only caller",
            "min_node_bp": MIN_DOMAIN,
            "shuffles": SHUFFLES,
            "seed_per_chromosome": SEED,
            "bootstraps": BOOTSTRAPS,
            "bootstrap_seed": 11,
            "controls": list(CONTROLS),
            "workers_do_not_change_numbers": True,
        },
        "exclusions": exclusions,
        "partitions": (
            "n/a: model answers over the whole genome, nothing fitted and nothing held out"
            if stage == 1
            else {
                "training": "EPCrisprBenchmark training_K562 Regulated=TRUE pairs",
                "heldout": "EPCrisprBenchmark heldout_5_cell_types Regulated=TRUE pairs",
                "pooled": "the containment statistic pools both: it fits nothing, so no partition is a test",
            }
        ),
    }


def main() -> None:
    args = sys.argv[1:]
    workers = int(args[args.index("--workers") + 1]) if "--workers" in args else 6
    stage = int(args[args.index("--stage") + 1]) if "--stage" in args else 1
    with ProcessPoolExecutor(workers) as pool:
        rows = [r for r in pool.map(chrom_row, CHROMS) if not r.get("skipped")]
    print("chromosomes:", [r["chrom"] for r in rows], flush=True)

    base = {
        "baseline_definitions": {
            "uniform": "the published control: as many uniform positions as the caller has edges, no "
            "50 kb minimum-size merge, 20 draws, seed 7 per chromosome "
            "(scripts/oriented_domains.py human_chrom)",
            "uniform_merged": "the same positions through the same 50 kb merge as every real node, so "
            "the node length floor is matched and only placement differs",
            "circular": "the real boundary set rotated by a uniform offset modulo the chromosome: node "
            "count and node length distribution matched (the one gap containing position 0 is split by "
            "the wrap), only the relation between a boundary and the genes near it destroyed",
            "count_matched": "uniform positions drawn until the post-merge edge count reaches the real one",
        },
        "evidence": {
            "nodes": "inferred: CTCF-only ENCODE elements as boundary proxies, no Hi-C",
            "modelled": "predicted: AlphaGenome deletion archive, most-moved coding gene",
            "measured": "experimental: ENCODE CRISPRi enhancer-gene benchmark, Regulated=TRUE",
        },
    }
    thin = [{k: v for k, v in r.items() if k not in ("node_clusters", "node_cluster_inside")} for r in rows]

    if stage == 1:
        census: dict[str, int] = {}
        for r in rows:
            for k, v in r["census"].items():
                census[k] = census.get(k, 0) + v
        census["coding_named"] = sum(r["coding_named"] for r in rows)
        census["any_named"] = sum(r["any_named"] for r in rows)
        print("denominator census", json.dumps(census), flush=True)
        modelled = readings(rows, "modelled: the deletion's most-moved coding gene")
        print("modelled", json.dumps({k: v for k, v in modelled.items() if k != "controls"}), flush=True)
        for c, v in modelled["controls"].items():
            print(
                f"  {c}: excess {v['excess_points']} ci_chrom {v['bootstrap_chromosomes']['ci95']} "
                f"ci_node {v['bootstrap_nodes']['ci95']} "
                f"sign {v['sign_test']['ahead']}/{v['sign_test']['of']} "
                f"p={v['sign_test']['p_two_sided']} "
                f"ctrl_med_len {v['mean_median_control_node_length']} ctrl_nodes {v['mean_control_nodes']}",
                flush=True,
            )
        save_result(
            "node_containment_audit",
            {"denominator_census": census, "per_chromosome": thin, "modelled": modelled} | base,
            manifest=manifest(1, rows),
        )
        print("saved stage 1", flush=True)
        return

    meas = measured_arm(rows)
    print("crispri inventory", json.dumps(meas["inventory"]), flush=True)
    measured = readings(meas["rows"], "measured: CRISPRi Regulated=TRUE element-gene pairs")
    print("measured", json.dumps({k: v for k, v in measured.items() if k != "controls"}), flush=True)
    for c, v in measured["controls"].items():
        print(
            f"  {c}: excess {v['excess_points']} ci_chrom {v['bootstrap_chromosomes']['ci95']} "
            f"ci_node {v['bootstrap_nodes']['ci95']} "
            f"sign {v['sign_test']['ahead']}/{v['sign_test']['of']} p={v['sign_test']['p_two_sided']}",
            flush=True,
        )
    save_result(
        "node_containment_measured",
        {
            "crispri_inventory": meas["inventory"],
            "measured": measured,
            "per_chromosome": [
                {k: v for k, v in r.items() if k not in ("node_clusters", "node_cluster_inside")}
                for r in meas["rows"]
            ],
        }
        | base,
        manifest=manifest(2, rows),
    )
    print("saved stage 2", flush=True)


if __name__ == "__main__":
    main()
