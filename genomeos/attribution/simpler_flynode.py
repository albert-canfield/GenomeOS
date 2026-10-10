# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Simpler genomes, the fly's boundary node: a measured boundary set against the worm's rule and nearest TSS.

Area I, simpler genomes, after the fly (genomeos/attribution/simpler_fly.py, docs/SIMPLER-GENOMES.md
sections 8 and 9). On Calderon et al. 2022 embryos the worm's gene-structure rule beat nearest TSS,
+2.70 points [+1.32, +4.10], but between genes it did not separate (+0.98 [-0.42, +2.31], 917
elements), so that run could not tell regulation from the host gene's own transcription. The human
node is a CTCF-bounded domain. The fly has measured boundaries (dCTCF, BEAF-32, CP190; Hi-C domain
boundaries) and the worm has none, so the fly can ask the human form of the question: does the gene
chosen inside an element's boundary-delimited domain match the element's tissue more often than the
gene with the nearest 5' end, where the chosen gene's own transcription through the element cannot
be the reason?

THE BOUNDARY SET. Stadler, Haines and Eisen 2017 (eLife 6:e29550, CC BY 4.0), Hi-C on hand-sorted
nuclear-cycle-14 embryos (2-3 h) at 500-bp resolution, supplementary file 1, which gives every list in
dm6 coordinates as the authors mapped them, so nothing is lifted. Primary: the authors' manual call
(Table S4, made by eye from the heatmaps genome-wide), the most complete of the three lists. Beside,
deciding nothing: the computational list (Table S5) and the merged representative set (Table S7).
Each boundary is a 500-bp bin; its position here is the bin's start plus 250.

THE NODE (fixed before any RNA value was read). An element's domain runs from the boundary at or left
of its midpoint to the next boundary right of it (an arm end when there is none); a gene is in the
domain when its 5' end lies there or inside either boundary's bin (boundary bins are shared, since
fly boundaries sit mostly at promoters). The rule names the gene with the nearest 5' end among the
domain's genes; when the arm-wide nearest 5' end is in the domain this is nearest TSS itself, and when
the domain holds no 5' end the rule falls back to nearest TSS. So the rule differs from nearest TSS
exactly where the nearest 5' end lies across a measured boundary and the domain holds another.

THE POPULATION (fixed with the position-only census in view, before any RNA value was read). The
lane's brief named the 917 between-gene elements as primary. There the boundary rule changes the
nearest-TSS gene for 90 elements, under the stop rule's 200, so between genes this boundary set cannot
decide. The primary population keeps the brief's reason and passes the stop rule: the distal elements
that lie outside the body of the gene the boundary rule names (every between-gene element, and the
elements inside a gene body when the rule names another gene). On it the rule's gene never contains
the element. The between-gene elements and all distal elements are reported beside.

The endpoint, the matched control, the bootstrap and the positive control are simpler_fly's, imported
and unchanged; so is the element population the subsets are drawn from.

WHAT THIS CANNOT ESTABLISH. Accessibility is not function, and a boundary is not a regulatory unit by
itself. The boundaries come from 2-3 h embryos while the marker calls come mostly from 6-16 h; that
fly boundaries persist after zygotic genome activation is the sources' report (Hug et al. 2017,
Cavalheiro et al. 2023), not checked here.
"""

from __future__ import annotations

import bisect
import gzip
import math
import random
from pathlib import Path
from typing import Any

from genomeos.attribution import simpler
from genomeos.attribution import simpler_fly as sf
from genomeos.attribution.simpler import Chromosome, Gene

STORE = Path("data/knowledge/simpler_flynode")
BOUNDARY_FILE = STORE / "elife-29550-supp1-v2.xlsx"
FETCH_RECORD = STORE / "fetch_record.json"
BOUNDARY_URL = "https://cdn.elifesciences.org/articles/29550/elife-29550-supp1-v2.xlsx"
BOUNDARY_SOURCE = (
    "Stadler MR, Haines JE, Eisen MB (2017) Convergence of topological domain boundaries, insulators, "
    "and polytene interbands revealed by high-resolution mapping of chromatin contacts in the early "
    "Drosophila melanogaster embryo. eLife 6:e29550, doi:10.7554/eLife.29550; supplementary file 1 "
    "(GEO GSE100370 holds the reads)"
)
BOUNDARY_LICENCE = "CC BY 4.0 (eLife: the article and its supplementary files)"

#: The three lists in supplementary file 1. The first is the primary; the others are beside.
SHEETS = {
    "manual": "Manual Boundaries",
    "computational": "Computational Boundaries",
    "representative": "Representative Boundaries",
}
PRIMARY_SET = "manual"
#: The Hi-C bin of every list; a boundary's position is its bin's start plus half of it.
BIN_BP = 500
HALF_BIN = BIN_BP // 2

#: Draws of the rotated-boundary comparison (beside, deciding nothing).
ROTATION_DRAWS = 1_000
#: Tissues are shown beside when they hold at least this many primary elements.
MIN_TISSUE = 50

#: The gene-structure rule's "between genes": its case is anything but a host gene.
HOST = "host"
#: The boundary rule's cases.
IN_DOMAIN = "nearest_in_domain"
MOVES = "boundary_moves"
EMPTY = "empty_domain"
CASES = (IN_DOMAIN, MOVES, EMPTY)

# simpler_fly's statistics and thresholds, unchanged.
SEED = sf.SEED
ALPHA = sf.ALPHA
MIN_DISCORDANT = sf.MIN_DISCORDANT
CLUSTER_BP = sf.CLUSTER_BP

#: lane-fly's registered result this run follows from (docs/SIMPLER-GENOMES.md section 8).
FLY = {
    "registration": "fd81e37",
    "all_distal_delta_points": 2.70,
    "all_distal_ci95_points": [1.32, 4.10],
    "between_genes_delta_points": 0.98,
    "between_genes_ci95_points": [-0.42, 2.31],
    "between_genes_elements": 917,
}


# --------------------------------------------------------------------------------------------- reading


def read_boundaries(path: Path, sheet: str) -> tuple[dict[str, list[int]], dict[str, Any]]:
    """One list of supplementary file 1: dm6 chromosome -> sorted boundary positions (bin start + 250).

    Every list gives the dm6 chromosome in column 0 and the bin start in column 1 (the computational
    and representative lists also give the bin end in column 2, always start + 499). Rows that are not
    a dm6 chromosome and an integer are headers and are counted, not read.
    """
    import openpyxl

    wb = openpyxl.load_workbook(path, read_only=True)
    per: dict[str, set[int]] = {}
    rows = data = 0
    widths: dict[int, int] = {}
    for r in wb[SHEETS[sheet]].iter_rows(values_only=True):
        rows += 1
        if not r or not isinstance(r[0], str) or not r[0].startswith("chr") or not isinstance(r[1], int):
            continue
        if r[0].removeprefix("chr") not in sf.CHROMOSOMES:
            continue
        data += 1
        if len(r) > 2 and isinstance(r[2], int):
            widths[r[2] - r[1]] = widths.get(r[2] - r[1], 0) + 1
        per.setdefault(r[0], set()).add(r[1] + HALF_BIN)
    bounds = {c: sorted(v) for c, v in sorted(per.items())}
    summary = {
        "sheet": SHEETS[sheet],
        "rows_in_sheet": rows,
        "boundary_rows": data,
        "boundaries": sum(len(v) for v in bounds.values()),
        "boundaries_by_chromosome": {c: len(v) for c, v in bounds.items()},
        "bin_end_minus_start": {str(k): v for k, v in sorted(widths.items())},
    }
    return bounds, summary


def read_arm_lengths(path: Path = sf.GFF3) -> dict[str, int]:
    """dm6 arm lengths from the GFF3's sequence-region lines, named the way the peaks name them."""
    out: dict[str, int] = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if not line.startswith("#"):
                break
            if line.startswith("##sequence-region"):
                _, name, _start, end = line.split()
                if name in sf.CHROMOSOMES:
                    out["chr" + name] = int(end)
    return out


# --------------------------------------------------------------------------------------------- the node


def domain_of(bounds: list[int], pos: int) -> tuple[float, float]:
    """The closed span whose 5' ends count as in the element's domain, boundary bins included."""
    i = bisect.bisect_right(bounds, pos)
    lo = bounds[i - 1] - HALF_BIN if i > 0 else -math.inf
    hi = bounds[i] + HALF_BIN if i < len(bounds) else math.inf
    return lo, hi


def boundary_rule(chrom: Chromosome, bounds: list[int], pos: int) -> tuple[Gene | None, str]:
    """The boundary node: the gene with the nearest 5' end inside the element's domain."""
    n = chrom.nearest_tss(pos)
    if n is None:
        return None, "none"
    lo, hi = domain_of(bounds, pos)
    if lo <= n.tss <= hi:
        return n, IN_DOMAIN
    a = bisect.bisect_left(chrom.tss, lo)
    b = bisect.bisect_right(chrom.tss, hi)
    cands = chrom.by_tss[a:b]
    if not cands:
        return n, EMPTY
    return min(cands, key=lambda g: (abs(g.tss - pos), g.gid)), MOVES


def attribute(
    rows: list[dict[str, Any]], chroms: dict[str, Chromosome], bounds: dict[str, list[int]]
) -> list[dict[str, Any]]:
    """`simpler.attribute` rows with the boundary rule's gene, case and 5'-end distance added.

    Also records whether the boundary gene and the nearest-TSS gene contain the element's midpoint.
    Reads gene positions and strands and the boundary positions only.
    """
    out = []
    for r in rows:
        ch = chroms[r["chrom"]]
        g, case = boundary_rule(ch, bounds.get(r["chrom"], []), r["mid"])
        if g is None:
            raise ValueError(f"no gene on {r['chrom']}")
        hosts = {h.gid for h in ch.hosts(r["mid"])}
        out.append(
            {
                **r,
                "boundary_gene": g.gid,
                "boundary_case": case,
                "boundary_dist": abs(g.tss - r["mid"]),
                "boundary_gene_contains_element": g.gid in hosts,
                "nearest_gene_contains_element": r["nearest_gene"] in hosts,
            }
        )
    return out


def distal(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [r for r in rows if r["kind"] == sf.DISTAL]


def primary(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """The primary population: distal elements outside the body of the gene the boundary rule names."""
    return [r for r in distal(rows) if not r["boundary_gene_contains_element"]]


def between_genes(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """The brief's population: distal elements whose midpoint lies in no protein-coding gene body."""
    return [r for r in distal(rows) if r["domain_case"] != HOST]


# --------------------------------------------------------------------------------------------- census


def _mde(k: int, n: int) -> float | None:
    """The 80%-power, two-sided 5% difference in points before clustering, as simpler.census bounds it."""
    return round(100 * 2.8 * math.sqrt(k) / n, 3) if n else None


def _population(rows: list[dict[str, Any]]) -> dict[str, Any]:
    n = len(rows)
    moves = [r for r in rows if r["boundary_gene"] != r["nearest_gene"]]
    k_gs = sum(1 for r in rows if r["boundary_gene"] != r["domain_gene"])
    return {
        "elements": n,
        "elements_by_tissue": {t: sum(1 for r in rows if r["tissue"] == t) for t in sf.TISSUES},
        "boundary_cases": {c: sum(1 for r in rows if r["boundary_case"] == c) for c in CASES},
        "differ_from_nearest_tss": len(moves),
        "differ_from_nearest_tss_inside_a_gene_body": sum(1 for r in moves if r["domain_case"] == HOST),
        "differ_from_nearest_tss_where_the_nearest_gene_contains_the_element": sum(
            1 for r in moves if r["nearest_gene_contains_element"]
        ),
        "differ_from_nearest_tss_where_the_boundary_gene_contains_the_element": sum(
            1 for r in moves if r["boundary_gene_contains_element"]
        ),
        "differ_from_gene_structure": k_gs,
        "gene_structure_names_a_host": sum(1 for r in rows if r["domain_case"] == HOST),
        "gene_structure_differs_from_nearest_tss": sum(
            1 for r in rows if r["domain_gene"] != r["nearest_gene"]
        ),
        "clusters_1mb": len({(r["chrom"], r["mid"] // CLUSTER_BP) for r in rows}),
        "mde_points_unclustered_vs_nearest_tss": _mde(len(moves), n),
        "mde_points_unclustered_vs_gene_structure": _mde(k_gs, n),
    }


def _domain_sizes(rows: list[dict[str, Any]], bounds: dict[str, list[int]], arm: dict[str, int]) -> dict:
    sizes = []
    for r in rows:
        lo, hi = domain_of(bounds.get(r["chrom"], []), r["mid"])
        a = max(lo + HALF_BIN, 0) if lo != -math.inf else 0
        b = min(hi - HALF_BIN, arm[r["chrom"]]) if hi != math.inf else arm[r["chrom"]]
        sizes.append(int(b - a))
    sizes.sort()
    if not sizes:
        return {}
    return {
        "median_bp": sizes[len(sizes) // 2],
        "quartiles_bp": [sizes[len(sizes) // 4], sizes[(3 * len(sizes)) // 4]],
    }


def _near(bounds: list[int], pos: int) -> list[int]:
    i = bisect.bisect_left(bounds, pos)
    return [bounds[j] for j in (i - 1, i) if 0 <= j < len(bounds)]


def census(
    rows: list[dict[str, Any]],
    alternatives: dict[str, list[dict[str, Any]]],
    bounds: dict[str, list[int]],
    arm: dict[str, int],
) -> dict[str, Any]:
    """Population counts and power in advance, from positions and the element-side tissue call only.

    `rows` carry the primary list's boundary gene; `alternatives` maps each other list to rows
    carrying its own.
    """
    prim = primary(rows)
    out = {
        "primary_outside_the_boundary_gene": _population(prim),
        "between_genes": _population(between_genes(rows)),
        "all_distal": _population(distal(rows)),
        "domain_size_of_primary_elements": _domain_sizes(prim, bounds, arm),
        "distal_elements_within_a_boundary_bin": sum(
            1
            for r in distal(rows)
            if any(abs(b - r["mid"]) <= HALF_BIN for b in _near(bounds.get(r["chrom"], []), r["mid"]))
        ),
        "alternative_lists_primary": {k: _population(primary(v)) for k, v in alternatives.items()},
        "positive_control_promoters": sum(1 for r in rows if r["kind"] == sf.PROMOTER),
    }
    out["stop_rule_met"] = (
        out["primary_outside_the_boundary_gene"]["differ_from_nearest_tss"] < MIN_DISCORDANT
    )
    return out


# --------------------------------------------------------------------------------------------- scoring


def compare(rows: list[dict[str, Any]], expr: dict[str, tuple[float, ...]], a: str, b: str) -> dict[str, Any]:
    """Concordance of the genes under keys `a` and `b` on the same elements, and their difference a - b."""
    ca = [int(sf.concordant(expr.get(r[a]), r["tissue"])) for r in rows]
    cb = [int(sf.concordant(expr.get(r[b]), r["tissue"])) for r in rows]
    diffs = [{"chrom": r["chrom"], "mid": r["mid"], "d": x - y} for r, x, y in zip(rows, ca, cb, strict=True)]
    bs = simpler.cluster_bootstrap(diffs, "d") if rows else None
    n = len(rows)
    return {
        "n": n,
        "a": a,
        "b": b,
        "differ": sum(1 for r in rows if r[a] != r[b]),
        "a_concordant": sum(ca),
        "b_concordant": sum(cb),
        "a_rate_pct": simpler._rate(sum(ca), n),
        "b_rate_pct": simpler._rate(sum(cb), n),
        "delta_points": round(100 * bs["point"], 3) if bs else None,
        "delta_ci95_points": [round(100 * bs["lo"], 3), round(100 * bs["hi"], 3)] if bs else None,
        "clusters": bs["clusters"] if bs else 0,
    }


def label(ci: list[float] | None, other: str) -> str:
    """BOUNDARY_RULE_AHEAD, `other`, or NO_DIFFERENCE_DETECTED, from the 95% interval of boundary - other."""
    if ci is None:
        return "NO_DIFFERENCE_DETECTED"
    lab = simpler.comparison_label(*ci)
    return {"DOMAIN_RULE_AHEAD": "BOUNDARY_RULE_AHEAD", "NEAREST_AHEAD": other}.get(lab, lab)


def rotation_null(
    rows: list[dict[str, Any]],
    chroms: dict[str, Chromosome],
    bounds: dict[str, list[int]],
    arm: dict[str, int],
    expr: dict[str, tuple[float, ...]],
    observed: int,
    draws: int = ROTATION_DRAWS,
    seed: int = SEED,
) -> dict[str, Any]:
    """The rule on the same elements with each arm's boundaries rotated by one uniform offset.

    The number and spacing of boundaries are kept; only where they lie changes.
    """
    rng = random.Random(seed)
    sims = []
    for _ in range(draws):
        rot = {}
        for c in sorted(bounds):
            off = rng.randrange(arm[c])
            rot[c] = sorted((b + off) % arm[c] for b in bounds[c])
        k = 0
        for r in rows:
            g, _ = boundary_rule(chroms[r["chrom"]], rot.get(r["chrom"], []), r["mid"])
            k += sf.concordant(expr.get(g.gid) if g else None, r["tissue"])
        sims.append(k)
    ge = sum(1 for s in sims if s >= observed)
    n = len(rows)
    return {
        "draws": draws,
        "observed": observed,
        "mean_rotated": round(sum(sims) / draws, 3),
        "mean_rotated_rate_pct": round(100 * sum(sims) / draws / n, 3) if n else None,
        "draws_at_or_above_observed": ge,
        "p_one_sided": round((1 + ge) / (1 + draws), 5),
    }


def _pts(x: float) -> str:
    return f"{x:+.2f}"


def _fill(res: dict[str, Any]) -> dict[str, str]:
    p = res["primary"]
    g = res["primary_vs_gene_structure"]
    bg = res["between_genes"]
    s = res["all_distal"]
    sg = res["all_distal_vs_gene_structure"]
    ctl = res["control"]
    rot = res["rotation_beside"]
    cen = res["census"]["primary_outside_the_boundary_gene"]
    return {
        "n": f"{p['n']:,}",
        "kb": f"{p['a_concordant']:,}",
        "kn": f"{p['b_concordant']:,}",
        "cb": f"{p['a_rate_pct']:.1f}",
        "cn": f"{p['b_rate_pct']:.1f}",
        "k": f"{p['differ']:,}",
        "d": _pts(p["delta_points"]),
        "lo": _pts(p["delta_ci95_points"][0]),
        "hi": _pts(p["delta_ci95_points"][1]),
        "knh": f"{cen['differ_from_nearest_tss_where_the_nearest_gene_contains_the_element']:,}",
        "mde": f"{cen['mde_points_unclustered_vs_nearest_tss']:.1f}",
        "kg": f"{g['b_concordant']:,}",
        "cg": f"{g['b_rate_pct']:.1f}",
        "k2": f"{g['differ']:,}",
        "gh": f"{cen['gene_structure_names_a_host']:,}",
        "d2": _pts(g["delta_points"]),
        "lo2": _pts(g["delta_ci95_points"][0]),
        "hi2": _pts(g["delta_ci95_points"][1]),
        "nbg": f"{bg['n']:,}",
        "kbg": f"{bg['differ']:,}",
        "abg": f"{bg['a_concordant']:,}",
        "cabg": f"{bg['a_rate_pct']:.1f}",
        "bbg": f"{bg['b_concordant']:,}",
        "cbbg": f"{bg['b_rate_pct']:.1f}",
        "dbg": _pts(bg["delta_points"]),
        "lobg": _pts(bg["delta_ci95_points"][0]),
        "hibg": _pts(bg["delta_ci95_points"][1]),
        "N": f"{s['n']:,}",
        "kA": f"{s['differ']:,}",
        "kbA": f"{s['a_concordant']:,}",
        "cbA": f"{s['a_rate_pct']:.1f}",
        "knA": f"{s['b_concordant']:,}",
        "cnA": f"{s['b_rate_pct']:.1f}",
        "dA": _pts(s["delta_points"]),
        "loA": _pts(s["delta_ci95_points"][0]),
        "hiA": _pts(s["delta_ci95_points"][1]),
        "kgA": f"{sg['b_concordant']:,}",
        "cgA": f"{sg['b_rate_pct']:.1f}",
        "dgA": _pts(sg["delta_points"]),
        "logA": _pts(sg["delta_ci95_points"][0]),
        "higA": _pts(sg["delta_ci95_points"][1]),
        "cc": f"{100 * ctl['expected_rate']:.1f}",
        "p": simpler._fmt_p(ctl["p_one_sided"]),
        "rm": f"{rot['mean_rotated_rate_pct']:.1f}",
        "pr": simpler._fmt_p(rot["p_one_sided"]),
        "t": f"{len(sf.TISSUES)}",
        "nb": f"{res['boundaries']:,}",
    }


def reading(res: dict[str, Any]) -> str:
    """The registered outcome wording, filled in. The templates are REGISTRATION['outcome_wording']."""
    w = REGISTRATION["outcome_wording"]
    fill = _fill(res)
    lab = res["comparison"]
    lab2 = res["comparison_vs_gene_structure"]
    falsified = res["falsifier"] == "FALSIFIED"
    rotated = res["rotation_beside"]["p_one_sided"] < ALPHA
    few = res["between_genes"]["differ"] < MIN_DISCORDANT
    parts = [
        w["head"],
        w[lab],
        w["meaning_" + lab],
        w["nearest_host"],
        w["gene_structure"],
        w["gene_structure_" + lab2],
        w["FALSIFIED" if falsified else "NOT_FALSIFIED"],
        w["rotation_" + ("ahead" if rotated else "not_ahead")],
        w["between_genes_" + ("few" if few else "enough")],
        w["all_distal"],
        w["tail"],
    ]
    return " ".join(x.format(**fill) for x in parts)


#: Written and committed before `score` read any RNA value. Not to be changed afterwards.
REGISTRATION: dict[str, Any] = {
    "claim": (
        "On D. melanogaster embryo distal elements accessible in clusters of exactly one tissue and lying "
        "outside the body of the gene the boundary rule names, the gene with the nearest 5' end inside "
        "the element's domain, delimited by measured Hi-C boundaries, is expressed highest in the "
        "element's tissue more often than the gene with the nearest 5' end on the arm."
    ),
    "why_this_test": (
        "lane-fly (registration fd81e37) found the worm's gene-structure rule ahead of nearest TSS on all "
        "2,518 distal elements (+2.70 [+1.32, +4.10]) but not between genes (+0.98 [-0.42, +2.31], 917 "
        "elements), so it could not separate regulation from the host gene's transcription. The human "
        "node is a boundary-delimited domain; the fly, unlike the worm, has measured boundaries."
    ),
    "boundary_set": {
        "source": BOUNDARY_SOURCE,
        "licence": BOUNDARY_LICENCE,
        "url": BOUNDARY_URL,
        "primary": (
            "Table S4, the authors' manual boundary call made from the Hi-C heatmaps at 500-bp resolution "
            "over the whole genome, in the dm6 coordinates the file gives; the most complete of the "
            "three lists (the authors estimate 4,000 to 5,500 boundaries in all)."
        ),
        "beside": (
            "Table S5 (computational, the authors' directionality score, 500-bp bins) and Table S7 (the "
            "representative merge of the two), each run through the primary comparison on its own "
            "primary population, deciding nothing."
        ),
        "why_this_set": (
            "embryos; dm6 as published, so nothing is lifted; 500-bp resolution, finer than the spacing "
            "of fly genes (a 10-kb boundary call, as in Cavalheiro et al. 2023, is coarser than many fly "
            "intergenic gaps); boundaries called by the authors, so no caller or threshold is chosen "
            "here; a stated licence (CC BY 4.0); 386,317 bytes."
        ),
        "stage": (
            "nuclear cycle 14, 2-3 h. The marker calls come mostly from 6-16 h and the RNA from 0-20 h. "
            "That boundaries formed at zygotic genome activation persist later is the sources' report "
            "(Hug et al. 2017, Cavalheiro et al. 2023), not checked here."
        ),
        "position": "each boundary is a 500-bp bin; its position is the bin start + 250",
    },
    "node": (
        "Boundary rule: the element's domain runs from the boundary at or left of its midpoint to the "
        "next boundary right of it, or to the arm end; a protein-coding gene is in the domain when its "
        "5' end lies within [left boundary - 250, right boundary + 250] (the boundary bins are shared by "
        "the two domains they separate). The rule names the in-domain gene with the nearest 5' end "
        "(ties by gene ID). If the arm-wide nearest 5' end is in the domain this is nearest TSS itself "
        "(case nearest_in_domain); if the domain holds no 5' end the rule falls back to nearest TSS "
        "(case empty_domain); otherwise it names a different gene (case boundary_moves). Reads Ensembl "
        "104 gene positions and strands and the boundary positions only."
    ),
    "comparators": {
        "nearest_tss": (
            "simpler_fly's baseline: nearest annotated 5' end among the arm's protein-coding genes"
        ),
        "gene_structure": (
            "simpler_fly's node, the worm's rule (simpler.domain_rule): host gene, else the flanking gene "
            "whose 5' end faces the element"
        ),
    },
    "population": (
        "Drawn from simpler_fly's 2,518 single-tissue distal elements of registration fd81e37, unchanged. "
        "Primary: the 1,634 whose midpoint lies outside the body of the gene the boundary rule names (all "
        "917 between-gene elements and 717 elements inside a gene body for which the rule names another "
        "gene), so the named gene never contains the element. Beside: the 917 between-gene elements "
        "(midpoint in no protein-coding gene body) and all 2,518 distal elements."
    ),
    "census_decision": (
        "The lane's brief named the 917 between-gene elements as the primary population. The census, "
        "computed from positions and the element-side tissue call only, before any RNA value was read, "
        "shows the boundary rule names a different gene from nearest TSS for 90 of them, under the stop "
        "rule's 200; between genes this boundary set cannot decide. The primary population keeps the "
        "brief's reason (the attributed gene's own transcription through the element is excluded) and "
        "holds 234 elements where the rules differ. The choice was made with that census in view; no "
        "RNA value informed it."
    ),
    "endpoint": "simpler_fly's, unchanged (concordant: >= 5 CPM in the element's tissue and highest of 14)",
    "primary_comparison": (
        "Boundary rule against nearest TSS on the primary population. delta = concordance(boundary rule) "
        "- concordance(nearest TSS), in points, with simpler_fly's percentile cluster bootstrap over 1-Mb "
        "windows (2,000 resamples, seed 20261010). BOUNDARY_RULE_AHEAD if the lower bound > 0; "
        "NEAREST_AHEAD if the upper bound < 0; NO_DIFFERENCE_DETECTED otherwise. This label is the "
        "outcome. Where the rules differ the nearest-TSS gene sometimes contains the element (115 of 234 "
        "in the census), which can favour nearest TSS; the comparison is conservative for the boundary "
        "rule there."
    ),
    "second_registered_comparison": (
        "Boundary rule against the gene-structure rule on the same elements, the same statistic; "
        "BOUNDARY_RULE_AHEAD, GENE_STRUCTURE_AHEAD or NO_DIFFERENCE_DETECTED. Reported with its own label; "
        "it does not change the outcome. On 717 of the primary elements the gene-structure rule names the "
        "host, so its own transcription can favour it. Two comparisons are made on one population; each "
        "is read at 95% with no correction for the pair."
    ),
    "beside_populations": (
        "Both comparisons on the between-gene elements and on all distal elements, with numbers, deciding "
        "nothing. Between genes the census holds 90 differing elements, under the stop rule's 200, and "
        "the wording says so."
    ),
    "matched_control": (
        "simpler_fly's, unchanged, applied to the boundary rule's gene: per element a protein-coding gene "
        "of the same arm whose 5' end lies in the same distance bin, drawn uniformly, the attributed gene "
        "excluded; 1,000 draws, seed 20261010; one-sided p = (1 + draws >= observed) / 1,001."
    ),
    "positive_control": (
        "simpler_fly's, unchanged: the single-tissue promoter elements, nearest-TSS gene against its "
        "matched control. If its p >= 0.05 the readout is declared unresponsive and the reading withheld."
    ),
    "falsifier": (
        "On the primary population, the boundary rule is falsified as a carrier of tissue information on "
        "this endpoint if its concordance is not above its matched-distance control (one-sided p >= 0.05)."
    ),
    "beside": (
        "Deciding nothing: (1) the rule on the primary elements with each arm's boundaries rotated by one "
        "uniform offset (number and spacing kept), 1,000 draws, seed 20261010, one-sided p = (1 + draws "
        "with concordance >= observed) / 1,001; the population stays fixed, so a rotated rule may name a "
        "gene that contains the element; (2) the primary comparison with the computational and the "
        "representative lists, each on its own primary population; (3) the primary elements where the "
        "rules differ, alone; (4) each tissue with 50 or more primary elements; (5) the gene-structure "
        "rule against nearest TSS on the primary elements."
    ),
    "power_in_advance": (
        "In the census before any RNA value was read: per population, the elements that differ and the "
        "smallest difference with 80% power at two-sided 5% before any inflation for clustering, "
        "2.8 sqrt(k) / n points."
    ),
    "stop_rule": (
        "One run. If fewer than 200 primary elements have different genes under the boundary rule and "
        "nearest TSS, stop unscored (UNDERPOWERED). No boundary list, rule, population, threshold or "
        "endpoint is changed after an RNA value is read; another definition would be a new registration, "
        "committed before it is run."
    ),
    "never_read": (
        "for the attribution: any accessibility or RNA value; at all: simpler_fly.NEVER_READ, and no "
        "sheet of the boundary file but Tables S4, S5 and S7"
    ),
    "cannot_establish": (
        "Accessibility is not function, and a boundary is not a regulatory unit by itself. Concordance "
        "says the attributed gene is expressed where the element is open, not that the element acts on "
        "it. The boundaries come from 2-3 h embryos, the marker calls mostly from 6-16 h, and "
        "accessibility and RNA from different nuclei. Fly boundaries sit mostly at active promoters, so "
        "the boundary node and the gene-structure node are not independent. The named gene's own "
        "transcription through the element is excluded; a neighbouring gene's transcription is not."
    ),
    "outcome_wording": {
        "head": (
            "On {n} D. melanogaster embryo distal elements accessible in clusters of exactly one of {t} "
            "tissues (Calderon et al. 2022, 0-20 h, the source's own marker call) and lying outside the "
            "body of the gene the boundary rule names, the gene with the nearest 5' end inside the "
            "element's domain, delimited by {nb} Hi-C boundaries of nuclear-cycle-14 embryos (Stadler et "
            "al. 2017, the authors' manual call), is expressed at 5 counts per million or more and "
            "highest in the element's own tissue, in RNA pooled over that tissue's cells, for {kb} of {n} "
            "({cb}%), against {kn} of {n} ({cn}%) for the gene with the nearest 5' end on the arm; the "
            "rules differ on {k} of {n} elements; difference {d} points, 95% {lo} to {hi} over 1-Mb "
            "clusters."
        ),
        "BOUNDARY_RULE_AHEAD": "On this endpoint the boundary rule is ahead of nearest TSS.",
        "NEAREST_AHEAD": "On this endpoint nearest TSS is ahead of the boundary rule.",
        "NO_DIFFERENCE_DETECTED": "No difference detected between the boundary rule and nearest TSS.",
        "meaning_BOUNDARY_RULE_AHEAD": (
            "The gene the boundary rule names never contains these elements, so its lead does not come "
            "from the element sitting inside that gene's transcribed body."
        ),
        "meaning_NEAREST_AHEAD": (
            "Where the nearest 5' end lies across a measured boundary, the gene across it matches the "
            "element's tissue more often than the gene inside the domain."
        ),
        "meaning_NO_DIFFERENCE_DETECTED": (
            "With {k} elements differing, a difference smaller than about {mde} points could go undetected."
        ),
        "nearest_host": (
            "For {knh} of the {k} elements where the rules differ, the nearest-TSS gene contains the "
            "element, which can favour nearest TSS."
        ),
        "gene_structure": (
            "Against the gene-structure rule on the same elements ({kg} of {n}, {cg}%; the two differ on "
            "{k2}, and the gene-structure rule names the host for {gh}), the difference is {d2} points, "
            "95% {lo2} to {hi2}:"
        ),
        "gene_structure_BOUNDARY_RULE_AHEAD": "the boundary rule is ahead of the gene-structure rule.",
        "gene_structure_GENE_STRUCTURE_AHEAD": "the gene-structure rule is ahead of the boundary rule.",
        "gene_structure_NO_DIFFERENCE_DETECTED": "no difference detected between the two nodes.",
        "NOT_FALSIFIED": (
            "The boundary rule's genes match the element's tissue more often than genes at a matched "
            "distance ({cc}%, one-sided p {p}), so position carries tissue information here."
        ),
        "FALSIFIED": (
            "The boundary rule's genes match the element's tissue no more often than genes at a matched "
            "distance ({cc}%, one-sided p {p}): on this endpoint the boundary rule carries no tissue "
            "information beyond distance, and it is falsified as a carrier of it."
        ),
        "rotation_ahead": (
            "With the same boundaries rotated along each arm the rule averages {rm}% on these elements, "
            "below the measured boundaries' (one-sided p {pr}), so where the boundaries lie, not only how "
            "many there are, carries information here."
        ),
        "rotation_not_ahead": (
            "With the same boundaries rotated along each arm the rule averages {rm}% on these elements "
            "(one-sided p {pr}), so the measured positions are not shown to carry more than the "
            "boundaries' number and spacing."
        ),
        "between_genes_few": (
            "Between genes (beside, deciding nothing) the boundary rule names a different gene from "
            "nearest TSS for only {kbg} of {nbg} elements, fewer than the stop rule's 200, so this "
            "boundary set cannot decide the between-gene question: there it is concordant for {abg} "
            "({cabg}%) against {bbg} ({cbbg}%) for nearest TSS, {dbg} points [{lobg}, {hibg}]."
        ),
        "between_genes_enough": (
            "Between genes (beside, deciding nothing) the boundary rule names a different gene from "
            "nearest TSS for {kbg} of {nbg} elements and is concordant for {abg} ({cabg}%) against {bbg} "
            "({cbbg}%) for nearest TSS, {dbg} points [{lobg}, {hibg}]."
        ),
        "all_distal": (
            "On all {N} distal elements (beside, deciding nothing; the rules differ on {kA}) the boundary "
            "rule is concordant for {kbA} ({cbA}%) against {knA} ({cnA}%) for nearest TSS, {dA} points "
            "[{loA}, {hiA}], and against {kgA} ({cgA}%) for the gene-structure rule, {dgA} points [{logA}, "
            "{higA}]."
        ),
        "tail": (
            "Accessibility is not function, and a boundary is not a regulatory unit by itself: this says "
            "the attributed gene is expressed where the element is open, not that the element acts on "
            "it. The boundaries come from 2-3 h embryos, the marker calls mostly from 6-16 h, and "
            "accessibility and RNA from different nuclei."
        ),
    },
}


def positive_control(
    rows: list[dict[str, Any]], idx: sf.ControlIndex, expr: dict[str, tuple[float, ...]]
) -> dict[str, Any]:
    """simpler_fly's positive control, as `simpler_fly.score` computes it."""
    prom = [r for r in rows if r["kind"] == sf.PROMOTER]
    kp = sum(int(sf.concordant(expr.get(r["nearest_gene"]), r["tissue"])) for r in prom)
    probs = [
        idx.probability(r["chrom"], r["tissue"], r["mid"], r["nearest_gene"], r["nearest_dist"]) for r in prom
    ]
    pos: dict[str, Any] = {
        "n": len(prom),
        "nearest_concordant": kp,
        "nearest_rate_pct": simpler._rate(kp, len(prom)),
    }
    pos["control"] = simpler.control_test(probs, kp)
    pos["responsive"] = pos["control"]["p_one_sided"] < ALPHA
    return pos


def score(
    chroms: dict[str, Chromosome],
    rows: list[dict[str, Any]],
    alternatives: dict[str, list[dict[str, Any]]],
    bounds: dict[str, list[int]],
    arm: dict[str, int],
    expr: dict[str, tuple[float, ...]],
) -> dict[str, Any]:
    """The registered run. `rows` and `alternatives` come from `attribute`; `expr` from simpler_fly."""
    cen = census(rows, alternatives, bounds, arm)
    if cen["stop_rule_met"]:
        return {"outcome": "UNDERPOWERED", "census": cen}
    prim = primary(rows)
    b_, n_, g_ = "boundary_gene", "nearest_gene", "domain_gene"
    idx = sf.ControlIndex(chroms, expr)
    first = compare(prim, expr, b_, n_)
    probs = [idx.probability(r["chrom"], r["tissue"], r["mid"], r[b_], r["boundary_dist"]) for r in prim]
    control = simpler.control_test(probs, first["a_concordant"])
    beside: dict[str, Any] = {
        "primary_where_the_rules_differ": compare([r for r in prim if r[b_] != r[n_]], expr, b_, n_),
        "primary_gene_structure_vs_nearest_tss": compare(prim, expr, g_, n_),
    }
    for name, alt in alternatives.items():
        beside[f"{name}_list_vs_nearest_tss"] = compare(primary(alt), expr, b_, n_)
    for t in sf.TISSUES:
        sub = [r for r in prim if r["tissue"] == t]
        if len(sub) >= MIN_TISSUE:
            beside["tissue_" + t] = compare(sub, expr, b_, n_)
    res: dict[str, Any] = {
        "census": cen,
        "boundaries": sum(len(v) for v in bounds.values()),
        "primary": first,
        "primary_vs_gene_structure": compare(prim, expr, b_, g_),
        "between_genes": compare(between_genes(rows), expr, b_, n_),
        "between_genes_vs_gene_structure": compare(between_genes(rows), expr, b_, g_),
        "all_distal": compare(distal(rows), expr, b_, n_),
        "all_distal_vs_gene_structure": compare(distal(rows), expr, b_, g_),
        "control": control,
        "positive_control": positive_control(rows, idx, expr),
        "rotation_beside": rotation_null(prim, chroms, bounds, arm, expr, first["a_concordant"]),
        "beside": beside,
        "attributed_genes_missing_from_rna_rows": sum(
            1 for r in distal(rows) for k in (b_, n_, g_) if r[k] not in expr
        ),
        "fly": FLY,
    }
    res["comparison"] = label(first["delta_ci95_points"], "NEAREST_AHEAD")
    res["comparison_vs_gene_structure"] = label(
        res["primary_vs_gene_structure"]["delta_ci95_points"], "GENE_STRUCTURE_AHEAD"
    )
    res["falsifier"] = "FALSIFIED" if control["p_one_sided"] >= ALPHA else "NOT_FALSIFIED"
    if not res["positive_control"]["responsive"]:
        res["outcome"] = "READOUT_UNRESPONSIVE"
        res["reading"] = (
            "The positive control did not separate from its matched control, so the readout is "
            "declared unresponsive and the primary reading is withheld."
        )
    else:
        res["outcome"] = res["comparison"]
        res["reading"] = reading(res)
    return res
