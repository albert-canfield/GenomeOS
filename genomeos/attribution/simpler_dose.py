# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Simpler genomes, stage 3: the worm's constraint-dose view, registered on a second source.

WHY. Stage 2 (`simpler_tiers`, registration 752f7fc) found no difference for the registered tier on
Serizay et al. 2020 (+2.63 [-1.35, +7.03], FALSIFIED). Beside it, deciding nothing, the dose view (the
top third of elements by constrained fraction against the bottom third, phyloP135way, cut 2.08) gave
+6.69 [+2.11, +11.68]. A post-hoc view is not a result. This module registers that view as a
hypothesis on an independent source before any of its RNA is read.

THE SOURCE. Jänes et al. 2018 (eLife 7:e37344, CC BY 4.0): ATAC-seq of whole wild-type animals at six
stages (embryo, L1, L2, L3, L4, young adult) and strand-specific long nuclear RNA (chromatin-associated,
over 200 nt, "long cap") at the same six stages, both by the same study. Elements come from the
article's source data (Figure 1, 2 and 4 source data 1); the RNA from the study's GEO tracks
(GSE114483, the linear, DESeq2-size-factor-normalised coverage averaged over two replicates).

ADAPTATIONS, each fixed before the RNA was read (REGISTRATION states them):
- the element's group is a stage, not a tissue: the source gives no per-stage accessibility call, so
  among the putative enhancers the source calls developmentally dynamic, the element's stage is the
  stage of its highest ATAC height of the six (ties to the earlier stage);
- the gene's RNA is coverage-TPM: the mean long-cap coverage over the union of the gene's exons on its
  own strand, scaled per stage to sum to one million over the protein-coding genes; this is TPM's
  construction (reads per base, normalised to a million) on coverage instead of counts;
- concordant: at least 5 coverage-TPM at the element's stage and strictly higher than at each of the
  other five stages (stage 1's 5 TPM and "highest", on six stages instead of five tissues).

Nothing here reads the source's own element-to-gene columns (`NEVER_READ`), and nothing sends a
model request.
"""

from __future__ import annotations

import bisect
import gzip
import math
from pathlib import Path
from typing import Any

import numpy as np

from genomeos.attribution import simpler
from genomeos.attribution import simpler_tiers as st
from genomeos.attribution.bigwig import BigWig, coalesce

STORE = Path("data/knowledge/simpler_dose")
ELIFE = "https://cdn.elifesciences.org/articles/37344"
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE114nnn/GSE114483/suppl"
FIG1 = STORE / "elife-37344-fig1-data1-v2.txt"
FIG2 = STORE / "elife-37344-fig2-data1-v2.txt"
FIG4 = STORE / "elife-37344-fig4-data1-v2.txt"
FETCH_RECORD = STORE / "fetch_record.json"
SOURCE = (
    "Jänes J et al. (2018) Chromatin accessibility dynamics across C. elegans development and ageing. "
    "eLife 7:e37344; source data (Figures 1, 2 and 4) and GEO GSE114483 (long cap RNA, wild type)"
)
LICENCE = (
    "article and source data CC BY 4.0 (eLife); GEO GSE114483: NCBI places no restriction on use, the "
    "submitters state none"
)

#: The six wild-type developmental stages, in time order (the source's labels).
STAGES = ("emb", "l1", "l2", "l3", "l4", "ya")
STRANDS = ("fwd", "rev")
LCAP = {(s, d): STORE / f"GSE114483_lcap_wt_{s}_linear_{d}.bw" for s in STAGES for d in STRANDS}
ELIFE_FILES = (FIG1, FIG2, FIG4)

ELEMENT_TYPE = simpler.ELEMENT_TYPE
PROMOTER_TYPE = simpler.PROMOTER_TYPE
DETECTION_TPM = simpler.DETECTION_TPM
#: The source's own element-to-gene columns (Figure 2 source data). Never read.
NEVER_READ = (
    "promoter_gene_id_fwd",
    "promoter_locus_id_fwd",
    "promoter_gene_biotype_fwd",
    "promoter_gene_id_rev",
    "promoter_locus_id_rev",
    "promoter_gene_biotype_rev",
    "associated_gene_id",
    "associated_locus_id",
)
#: Below this many elements in the top or the bottom third (strata holding both), stop unscored.
MIN_PER_ARM = 200
#: The stage-2 view's tertile edges of the constrained fraction (simpler_worm_tiers_registration census).
STAGE2_EDGES = (0.08, 0.190476)
#: Stage 2's primary-track cut (simpler_worm_tiers.json, phyloP135way).
STAGE2_CUT = 2.08
PRIMARY = st.PRIMARY
PRIMARY_URL = st.TRACKS[PRIMARY]["url"]


# --------------------------------------------------------------------------------------------- elements


def _table(path: Path):
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) == len(head):
                yield dict(zip(head, p, strict=True))


def stage_of(heights: list[float]) -> str:
    """The stage of the highest ATAC height; a tie goes to the earlier stage."""
    best = max(heights)
    return STAGES[heights.index(best)]


def read_elements(
    fig1: Path = FIG1,
    fig2: Path = FIG2,
    fig4: Path = FIG4,
    kinds: tuple[str, ...] = (ELEMENT_TYPE, PROMOTER_TYPE),
) -> tuple[list[simpler.Element], dict[str, Any]]:
    """Developmentally dynamic elements of the given classes, each with its stage of highest
    accessibility. Reads coordinates, the class, the dynamic call and the six wild-type ATAC heights;
    the RNA-derived flag columns and the gene columns are never read."""
    dynamic: dict[tuple[str, str, str], bool] = {}
    for r in _table(fig4):
        dynamic[(r["chrom_ce10"], r["start_ce10"], r["end_ce10"])] = r["devel_is_dynamic"] == "True"
    heights: dict[tuple[str, str, str], list[float]] = {}
    for r in _table(fig1):
        heights[(r["chrom"], r["start"], r["end"])] = [float(r[f"atac_wt_{s}_height"]) for s in STAGES]
    out: list[simpler.Element] = []
    funnel: dict[str, Any] = {"rows": 0, "by_class": {}, "dynamic_by_class": {}, "no_ce11": 0, "ties": 0}
    for r in _table(fig2):
        funnel["rows"] += 1
        kind = r["annot"]
        funnel["by_class"][kind] = funnel["by_class"].get(kind, 0) + 1
        key = (r["chrom_ce10"], r["start_ce10"], r["end_ce10"])
        if kind not in kinds or not dynamic.get(key, False):
            continue
        funnel["dynamic_by_class"][kind] = funnel["dynamic_by_class"].get(kind, 0) + 1
        if not r["chrom_ce11"] or r["start_ce11"] in ("", "NA") or r["end_ce11"] in ("", "NA"):
            funnel["no_ce11"] += 1
            continue
        h = heights[key]
        if h.count(max(h)) > 1:
            funnel["ties"] += 1
        out.append(
            simpler.Element(r["chrom_ce11"], int(r["start_ce11"]), int(r["end_ce11"]), kind, stage_of(h))
        )
    return out, funnel


def attribute(elements: list[simpler.Element], chroms: dict[str, simpler.Chromosome]) -> list[dict[str, Any]]:
    """Each element with the domain rule's gene and the nearest-TSS gene (positions and strands only)."""
    rows = []
    for e in elements:
        ch = chroms.get(e.chrom)
        if ch is None:
            continue
        d, case = simpler.domain_rule(ch, e.mid)
        n = ch.nearest_tss(e.mid)
        if d is None or n is None:
            continue
        rows.append(
            {
                "chrom": e.chrom,
                "start": e.start,
                "end": e.end,
                "mid": e.mid,
                "kind": e.kind,
                "tissue": e.tissue,  # the stage; simpler_tiers' strata read this key
                "domain_gene": d.gid,
                "domain_case": case,
                "domain_dist": abs(d.tss - e.mid),
                "nearest_gene": n.gid,
                "nearest_dist": abs(n.tss - e.mid),
            }
        )
    return rows


def overlaps_any(rows: list[dict[str, Any]], others: list[tuple[str, int, int]]) -> list[bool]:
    """Per row, whether [start, end) overlaps any of `others` (chrom, start, end)."""
    by: dict[str, list[tuple[int, int]]] = {}
    for c, s, e in others:
        by.setdefault(c, []).append((s, e))
    merged = {}
    for c, v in by.items():
        s, e = st.merge(v)
        merged[c] = (s.tolist(), e.tolist())
    out = []
    for r in rows:
        m = merged.get(r["chrom"])
        if m is None:
            out.append(False)
            continue
        starts, ends = m
        i = bisect.bisect_right(starts, r["end"] - 1)
        out.append(i > 0 and ends[i - 1] > r["start"])
    return out


# --------------------------------------------------------------------------------------------- RNA


def read_gene_exons(path: Path = simpler.GFF3) -> dict[str, tuple[str, str, list[tuple[int, int]]]]:
    """Per protein-coding gene: chromosome, strand and the merged exons of all its transcripts,
    0-based half-open (the GFF3 is 1-based closed)."""
    genes: dict[str, tuple[str, str]] = {}
    tx: dict[str, str] = {}
    exons: dict[str, list[tuple[int, int]]] = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[0] == "MtDNA":
                continue
            attrs = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
            if p[2] == "gene":
                if "biotype=protein_coding" in p[8]:
                    gid = attrs.get("gene_id") or attrs["ID"].split(":")[-1]
                    genes[gid] = ("chr" + p[0], p[6])
            elif p[2] == "exon":
                parent = attrs.get("Parent", "")
                if parent.startswith("transcript:"):
                    exons.setdefault(parent.split(":", 1)[1], []).append((int(p[3]) - 1, int(p[4])))
            elif attrs.get("ID", "").startswith("transcript:") and attrs.get("Parent", "").startswith(
                "gene:"
            ):
                tx[attrs["ID"].split(":", 1)[1]] = attrs["Parent"].split(":", 1)[1]
    per: dict[str, list[tuple[int, int]]] = {}
    for t, ivs in exons.items():
        g = tx.get(t)
        if g in genes:
            per.setdefault(g, []).extend(ivs)
    out = {}
    for g, ivs in per.items():
        s, e = st.merge(ivs)
        out[g] = (genes[g][0], genes[g][1], list(zip(s.tolist(), e.tolist(), strict=True)))
    return out


def chrom_values(bw: BigWig, chrom: str) -> np.ndarray:
    """The whole chromosome as per-base float32, zero where the track has no value."""
    size = bw.chroms[chrom][1]
    vals = np.zeros(size, dtype=np.float32)
    cid = bw.chroms[chrom][0]
    items = bw.leaf_items(chrom, [(0, size)])
    for g in coalesce(items):
        base = g[0].offset
        blob = bw.src.read(base, g[-1].offset + g[-1].size - base)
        for sec_cid, raw in st._sections(bw, g, blob):
            if sec_cid != cid:
                continue
            for start, a in st._runs(raw):
                end = min(start + len(a), size)
                if end > start:
                    vals[start:end] = a[: end - start]
    return vals


def exon_sums(vals: np.ndarray, exons: list[tuple[int, int]]) -> tuple[float, int]:
    """Summed coverage over the exons, and their length."""
    total = 0.0
    length = 0
    for s, e in exons:
        e = min(e, len(vals))
        if e > s:
            total += float(vals[s:e].sum(dtype=np.float64))
            length += e - s
    return total, length


def coverage(
    paths: dict[tuple[str, str], Path], gene_exons: dict[str, tuple[str, str, list[tuple[int, int]]]]
) -> dict[tuple[str, str], dict[str, float]]:
    """Per (stage, track strand) and gene: mean coverage over the gene's exon union."""
    by_chrom: dict[str, list[str]] = {}
    for g, (c, _, _) in gene_exons.items():
        by_chrom.setdefault(c, []).append(g)
    out: dict[tuple[str, str], dict[str, float]] = {}
    for key, path in paths.items():
        bw = BigWig(path)
        cov: dict[str, float] = {}
        for c in sorted(by_chrom):
            if c not in bw.chroms:
                continue
            vals = chrom_values(bw, c)
            for g in by_chrom[c]:
                tot, n = exon_sums(vals, gene_exons[g][2])
                cov[g] = tot / n if n else 0.0
        bw.close()
        out[key] = cov
    return out


def orientation(
    cov: dict[tuple[str, str], dict[str, float]],
    gene_exons: dict[str, tuple[str, str, list[tuple[int, int]]]],
) -> dict[str, Any]:
    """The registered strand check: summed mean coverage on each gene's own strand against the other,
    per stage. 'as_labelled' if sense wins at all six stages, 'swapped' if antisense wins at all six,
    'failed' otherwise."""
    per = {}
    for s in STAGES:
        sense = anti = 0.0
        for g, (_, strand, _) in gene_exons.items():
            f, r = cov[(s, "fwd")].get(g, 0.0), cov[(s, "rev")].get(g, 0.0)
            sense += f if strand == "+" else r
            anti += r if strand == "+" else f
        per[s] = {"sense": round(sense, 3), "antisense": round(anti, 3)}
    if all(v["sense"] > v["antisense"] for v in per.values()):
        verdict = "as_labelled"
    elif all(v["antisense"] > v["sense"] for v in per.values()):
        verdict = "swapped"
    else:
        verdict = "failed"
    return {"per_stage": per, "verdict": verdict}


def expression(
    cov: dict[tuple[str, str], dict[str, float]],
    gene_exons: dict[str, tuple[str, str, list[tuple[int, int]]]],
    verdict: str,
    genes: set[str],
) -> dict[str, tuple[float, ...]]:
    """Coverage-TPM per gene and stage: the gene's own-strand mean coverage, scaled per stage to sum to
    one million over `genes` (the protein-coding genes of the six nuclear chromosomes)."""
    flip = verdict == "swapped"
    sense: dict[str, list[float]] = {}
    for g in genes:
        if g not in gene_exons:
            continue
        strand = gene_exons[g][1]
        own = "fwd" if (strand == "+") != flip else "rev"
        sense[g] = [cov[(s, own)].get(g, 0.0) for s in STAGES]
    totals = [sum(v[i] for v in sense.values()) for i in range(len(STAGES))]
    return {
        g: tuple(1e6 * v[i] / totals[i] if totals[i] > 0 else 0.0 for i in range(len(STAGES)))
        for g, v in sense.items()
    }


def concordant(tpm: tuple[float, ...] | None, stage: str) -> bool:
    """At least DETECTION_TPM at the stage and strictly higher there than at each of the other five."""
    if tpm is None:
        return False
    i = STAGES.index(stage)
    v = tpm[i]
    return v >= DETECTION_TPM and all(v > tpm[j] for j in range(len(STAGES)) if j != i)


class ControlIndex(simpler.ControlIndex):
    """simpler.ControlIndex with the six stages in place of the five tissues: a matched-control gene
    is drawn from the element's chromosome at the same 5'-end distance bin as the attributed gene."""

    def __init__(self, chroms: dict[str, simpler.Chromosome], expr: dict[str, tuple[float, ...]]):
        self.chroms = chroms
        self.prefix = {}
        for c, ch in chroms.items():
            for t in STAGES:
                acc, p = 0, [0]
                for g in ch.by_tss:
                    acc += concordant(expr.get(g.gid), t)
                    p.append(acc)
                self.prefix[(c, t)] = p
        self.expr = expr

    def probability(self, chrom: str, tissue: str, mid: int, gene: str, dist: int) -> float:
        b = simpler.distance_bin(dist)
        lo = 0 if b == 0 else simpler.DISTANCE_BINS[b - 1]
        hi = simpler.DISTANCE_BINS[b] if b < len(simpler.DISTANCE_BINS) else math.inf
        n1, c1 = self._count(chrom, tissue, mid - hi + 1, mid - lo + 1)
        n2, c2 = self._count(chrom, tissue, mid + lo, mid + hi)
        if lo == 0:
            n0, c0 = self._count(chrom, tissue, mid, mid + 1)
            n1, c1 = n1 - n0, c1 - c0
        total, conc = n1 + n2, c1 + c2
        own = concordant(self.expr.get(gene), tissue)
        total, conc = total - 1, conc - int(own)
        if total > 0:
            return conc / total
        others = [g for g in self.chroms[chrom].genes if g.gid != gene]
        if not others:
            return 0.0
        sub = min(others, key=lambda g: (abs(abs(g.tss - mid) - dist), g.gid))
        return float(concordant(self.expr.get(sub.gid), tissue))


# --------------------------------------------------------------------------------------------- the dose


def attach_counts(rows: list[dict[str, Any]], counts: list[list[Any]]) -> list[dict[str, Any]]:
    """The population with each element's exon-free scored bases and bases at or above the cut, joined
    by (chrom, start, end), its stage-2 tier and the measured length used for the length strata."""
    by_key = {(c[0], c[1], c[2]): (c[3], c[4]) for c in counts}
    out = []
    for r in rows:
        m, a = by_key[(r["chrom"], r["start"], r["end"])]
        out.append({**r, "measured": m, "above": a, "tier": st.tier(m, a), "length_measured": m})
    return out


def by_edges(rows: list[dict[str, Any]], edges: tuple[float, float] | list[float]) -> list[dict[str, Any]]:
    """Top third ('constrained') against bottom third ('not_constrained') at fixed edges; the middle
    and the unmeasured are left out. simpler_tiers.by_fraction_tertile, the stage-2 view's own rule."""
    return st.by_fraction_tertile(rows, list(edges))


def middle_vs_bottom(rows: list[dict[str, Any]], edges: list[float]) -> list[dict[str, Any]]:
    """The middle third relabelled 'constrained' against the bottom third: the gradient, beside."""
    out = []
    for r in rows:
        if r["measured"] <= 0:
            continue
        f = r["above"] / r["measured"]
        if f < edges[0]:
            out.append({**r, "tier": "not_constrained"})
        elif not (f >= edges[1] and f > edges[0]):
            out.append({**r, "tier": "constrained"})
    return out


def _support(rows: list[dict[str, Any]], ledges: list[float]) -> tuple[int, int, int, int]:
    groups: dict[tuple, list[int]] = {}
    for r in rows:
        g = groups.setdefault(st.stratum(r, ledges), [0, 0])
        g[0 if r["tier"] == "constrained" else 1] += 1
    sup = [k for k, (a, b) in groups.items() if a and b]
    return len(groups), len(sup), sum(groups[k][0] for k in sup), sum(groups[k][1] for k in sup)


def _mde(nc: int, nu: int) -> float | None:
    return round(100 * 2.8 * math.sqrt(0.25 * (1 / nc + 1 / nu)), 3) if nc and nu else None


def census(
    rows: list[dict[str, Any]], promoters: list[dict[str, Any]], overlap: list[bool]
) -> dict[str, Any]:
    """Population counts, edges, common support and power in advance, from positions, the ATAC calls
    and the conservation track only. No RNA value is read."""
    measured = [r for r in rows if r["measured"] > 0]
    fedges = st.fraction_edges(rows)
    ledges = st.length_edges([r["measured"] for r in rows])
    dose = by_edges(rows, fedges)
    n_s, n_sup, n_top, n_bot = _support(dose, ledges)
    stage2 = by_edges(rows, STAGE2_EDGES)
    _, _, s2_top, s2_bot = _support(stage2, ledges)
    fresh = [r for r, o in zip(rows, overlap, strict=True) if not o]
    fresh_dose = by_edges(fresh, fedges)
    _, _, f_top, f_bot = _support(fresh_dose, ledges)
    disc = [r for r in rows if r["domain_gene"] != r["nearest_gene"]]
    k, n = len(disc), len(rows)

    def arms(sub: list[dict[str, Any]]) -> dict[str, int]:
        return {
            "top": sum(1 for r in sub if r["tier"] == "constrained"),
            "bottom": sum(1 for r in sub if r["tier"] == "not_constrained"),
        }

    return {
        "elements": n,
        "elements_by_stage": {s: sum(1 for r in rows if r["tissue"] == s) for s in STAGES},
        "measured": len(measured),
        "not_measured": n - len(measured),
        "overlap_stage2_population": sum(overlap),
        "overlap_stage2_population_measured": sum(
            1 for r, o in zip(rows, overlap, strict=True) if o and r["measured"] > 0
        ),
        "fraction_tertile_edges": [round(x, 6) for x in fedges],
        "length_edges_measured_bases": [round(x, 3) for x in ledges],
        "dose_arms": arms(dose),
        "dose_arms_by_stage": {s: arms([r for r in dose if r["tissue"] == s]) for s in STAGES},
        "dose_arms_inside_gene_body": arms([r for r in dose if r["domain_case"] == "host"]),
        "dose_arms_between_genes": arms([r for r in dose if r["domain_case"] != "host"]),
        "strata": n_s,
        "strata_with_both_arms": n_sup,
        "top_in_support": n_top,
        "bottom_in_support": n_bot,
        "mde_points_unclustered": _mde(n_top, n_bot),
        "stage2_edges_arms_in_support": {"top": s2_top, "bottom": s2_bot},
        "not_overlapping_stage2_arms_in_support": {"top": f_top, "bottom": f_bot},
        "tier_5pct": {
            "constrained": sum(1 for r in measured if r["tier"] == "constrained"),
            "not_constrained": sum(1 for r in measured if r["tier"] == "not_constrained"),
        },
        "secondary": {
            "discordant": k,
            "discordant_fraction": round(k / n, 4) if n else None,
            "discordant_inside_gene_body": sum(1 for r in disc if r["domain_case"] == "host"),
            "mde_points_unclustered": round(100 * 2.8 * math.sqrt(k) / n, 3) if n else None,
            "stop_rule_met": k < simpler.MIN_DISCORDANT,
        },
        "positive_control_promoters": len(promoters),
        "positive_control_promoters_by_stage": {
            s: sum(1 for r in promoters if r["tissue"] == s) for s in STAGES
        },
        "clusters_1mb": len({(r["chrom"], r["mid"] // simpler.CLUSTER_BP) for r in rows}),
        "stop_rule_met": n_top < MIN_PER_ARM or n_bot < MIN_PER_ARM,
    }


# --------------------------------------------------------------------------------------------- scoring


def dose_label(lo: float, hi: float) -> str:
    if lo > 0:
        return "DOSE_ADDS"
    if hi < 0:
        return "DOSE_BELOW"
    return "NO_DIFFERENCE_DETECTED"


def compare(rows: list[dict[str, Any]], expr: dict[str, tuple[float, ...]]) -> dict[str, Any]:
    """simpler.compare on the stage endpoint: domain rule against nearest TSS, cluster bootstrap."""
    for r in rows:
        r["c_domain"] = int(concordant(expr.get(r["domain_gene"]), r["tissue"]))
        r["c_nearest"] = int(concordant(expr.get(r["nearest_gene"]), r["tissue"]))
        r["d"] = r["c_domain"] - r["c_nearest"]
    n = len(rows)
    kd = sum(r["c_domain"] for r in rows)
    kn = sum(r["c_nearest"] for r in rows)
    bs = simpler.cluster_bootstrap(rows, "d") if rows else None
    return {
        "n": n,
        "differ": sum(1 for r in rows if r["domain_gene"] != r["nearest_gene"]),
        "domain_concordant": kd,
        "nearest_concordant": kn,
        "domain_rate_pct": round(100 * kd / n, 3) if n else None,
        "nearest_rate_pct": round(100 * kn / n, 3) if n else None,
        "delta_points": round(100 * bs["point"], 3) if bs else None,
        "delta_ci95_points": [round(100 * bs["lo"], 3), round(100 * bs["hi"], 3)] if bs else None,
        "clusters": bs["clusters"] if bs else 0,
    }


def control(rows: list[dict[str, Any]], ctl: ControlIndex, gene_key: str, dist_key: str, hit_key: str):
    probs = [ctl.probability(r["chrom"], r["tissue"], r["mid"], r[gene_key], r[dist_key]) for r in rows]
    return simpler.control_test(probs, sum(r[hit_key] for r in rows))


def score(
    rows: list[dict[str, Any]],
    promoters: list[dict[str, Any]],
    overlap: list[bool],
    chroms: dict[str, simpler.Chromosome],
    expr: dict[str, tuple[float, ...]],
    cen: dict[str, Any],
    strand: dict[str, Any],
) -> dict[str, Any]:
    """The registered run. `expr` (the endpoint) is read only here."""
    if strand["verdict"] == "failed":
        res: dict[str, Any] = {"outcome": "INSTRUMENT_FAILED", "strand_check": strand}
        res["reading"] = REGISTRATION["outcome_wording"]["INSTRUMENT_FAILED"]
        return res
    ctl = ControlIndex(chroms, expr)
    for r in promoters:
        r["c_nearest"] = int(concordant(expr.get(r["nearest_gene"]), r["tissue"]))
    pc = control(promoters, ctl, "nearest_gene", "nearest_dist", "c_nearest")
    pos = {
        "n": len(promoters),
        "concordant": sum(r["c_nearest"] for r in promoters),
        "rate_pct": round(100 * sum(r["c_nearest"] for r in promoters) / len(promoters), 3)
        if promoters
        else None,
        "matched_control": pc,
        "separates": pc["p_one_sided"] < simpler.ALPHA,
    }
    base: dict[str, Any] = {"strand_check": strand, "positive_control": pos}
    if not pos["separates"]:
        res = {**base, "outcome": "UNRESPONSIVE"}
        res["reading"] = reading(res)
        return res
    if cen["stop_rule_met"]:
        res = {**base, "outcome": "UNDERPOWERED"}
        res["reading"] = reading(res)
        return res
    for r in rows:
        r["c_domain"] = int(concordant(expr.get(r["domain_gene"]), r["tissue"]))
    fedges = cen["fraction_tertile_edges"]
    ledges = cen["length_edges_measured_bases"]
    primary = st.matched_difference(by_edges(rows, fedges), ledges)
    label = dose_label(*primary["delta_ci95_points"])
    falsifier = "FALSIFIED" if primary["permutation"]["p_one_sided"] >= simpler.ALPHA else "NOT_FALSIFIED"
    fresh = [r for r, o in zip(rows, overlap, strict=True) if not o]
    beside = {
        "stage2_edges": st.matched_difference(by_edges(rows, STAGE2_EDGES), ledges),
        "not_overlapping_stage2": st.matched_difference(by_edges(fresh, fedges), ledges),
        "inside_gene_body": st.matched_difference(
            by_edges([r for r in rows if r["domain_case"] == "host"], fedges), ledges
        ),
        "between_genes": st.matched_difference(
            by_edges([r for r in rows if r["domain_case"] != "host"], fedges), ledges
        ),
        "middle_vs_bottom": st.matched_difference(middle_vs_bottom(rows, fedges), ledges),
        "tier_5pct": st.matched_difference(rows, ledges),
    }
    for s in STAGES:
        beside[f"stage_{s}"] = st.matched_difference(
            by_edges([r for r in rows if r["tissue"] == s], fedges), ledges
        )
    res = {
        **base,
        "outcome": label,
        "falsifier": falsifier,
        "primary": primary,
        "beside": beside,
        "domain_genes_absent_from_rna": sum(1 for r in rows if r["domain_gene"] not in expr),
    }
    # The secondary: stage 1's comparison on this source, with its own label.
    if cen["secondary"]["stop_rule_met"]:
        res["secondary"] = {"outcome": "UNDERPOWERED"}
    else:
        cmp_all = compare(rows, expr)
        ctl_d = control(rows, ctl, "domain_gene", "domain_dist", "c_domain")
        sec = {
            **cmp_all,
            "outcome": simpler.comparison_label(*cmp_all["delta_ci95_points"]),
            "matched_control": ctl_d,
            "falsifier": "FALSIFIED" if ctl_d["p_one_sided"] >= simpler.ALPHA else "NOT_FALSIFIED",
            "inside_gene_body": compare([r for r in rows if r["domain_case"] == "host"], expr),
            "between_genes": compare([r for r in rows if r["domain_case"] != "host"], expr),
        }
        res["secondary"] = sec
    res["reading"] = reading(res)
    return res


def _ci(x: dict[str, Any]) -> str:
    lo, hi = x["delta_ci95_points"]
    return f"{x['delta_points']:+.2f} [{lo:+.2f}, {hi:+.2f}]"


def reading(res: dict[str, Any]) -> str:
    """The registered wording, filled in. Templates: REGISTRATION['outcome_wording']."""
    w = REGISTRATION["outcome_wording"]
    out = res["outcome"]
    if out == "INSTRUMENT_FAILED":
        return w[out]
    pos = res["positive_control"]
    pvals = {
        "pn": f"{pos['n']:,}",
        "pk": f"{pos['concordant']:,}",
        "pr": f"{pos['rate_pct']:.1f}%",
        "pc": f"{100 * pos['matched_control']['expected_rate']:.1f}%",
        "pp": f"{pos['matched_control']['p_one_sided']:.3g}",
    }
    if out == "UNRESPONSIVE":
        return w[out].format(**pvals)
    if out == "UNDERPOWERED":
        return w[out]
    p = res["primary"]
    lo, hi = p["delta_ci95_points"]
    vals = {
        "nc": f"{p['constrained_in_support']:,}",
        "kc": f"{p['constrained_concordant_in_support']:,}",
        "nu": f"{p['not_constrained_in_support']:,}",
        "rc": f"{p['constrained_rate_pct']:.1f}%",
        "ru": f"{p['matched_not_constrained_rate_pct']:.1f}%",
        "d": f"{p['delta_points']:+.2f}",
        "lo": f"{lo:+.2f}",
        "hi": f"{hi:+.2f}",
        "perm": f"{p['permutation']['p_one_sided']:.3g}",
        "falsifier": res["falsifier"],
    }
    text = w["common"].format(**vals) + " " + w[out]
    fresh = res["beside"]["not_overlapping_stage2"]
    if out == "DOSE_ADDS" and (not fresh.get("delta_ci95_points") or fresh["delta_ci95_points"][0] <= 0):
        text += " " + w["overlap_caveat"].format(
            fresh=_ci(fresh) if fresh.get("delta_ci95_points") else "not estimable"
        )
    text += " " + w["positive_control"].format(**pvals)
    sec = res["secondary"]
    if sec["outcome"] == "UNDERPOWERED":
        text += " " + w["secondary_underpowered"]
    else:
        slo, shi = sec["delta_ci95_points"]
        svals = {
            "n": f"{sec['n']:,}",
            "kd": f"{sec['domain_concordant']:,}",
            "rd": f"{sec['domain_rate_pct']:.1f}%",
            "kn": f"{sec['nearest_concordant']:,}",
            "rn": f"{sec['nearest_rate_pct']:.1f}%",
            "k": f"{sec['differ']:,}",
            "d": f"{sec['delta_points']:+.2f}",
            "lo": f"{slo:+.2f}",
            "hi": f"{shi:+.2f}",
            "cr": f"{100 * sec['matched_control']['expected_rate']:.1f}%",
            "cp": f"{sec['matched_control']['p_one_sided']:.3g}",
        }
        text += " " + w["secondary_common"].format(**svals) + " " + w[f"secondary_{sec['outcome']}"]
        key = (
            "secondary_control_passed" if sec["falsifier"] == "NOT_FALSIFIED" else "secondary_control_failed"
        )
        text += " " + w[key].format(**svals)
    return text + " " + w["always"]


REGISTRATION: dict[str, Any] = {
    "hypothesis": (
        "Stage 2's dose view, as a registered claim on a second source: among C. elegans putative "
        "enhancers, the domain-rule gene is concordant with the element's own group more often for the "
        "top third of elements by constrained fraction than for the bottom third, at the same group, "
        "position class, distance to that gene and measured length."
    ),
    "why_registered_here": (
        "On Serizay et al. 2020 the view gave +6.69 [+2.11, +11.68] beside a FALSIFIED registered tier "
        "(stage 2, registration 752f7fc). It was seen after that registration's run, so it cannot be "
        "re-tested on that source; it is registered here before any RNA value of this source is read."
    ),
    "source": (
        "Jänes et al. 2018 (eLife 7:e37344, CC BY 4.0): ATAC-seq of whole wild-type animals at six "
        "stages (embryo, L1, L2, L3, L4, young adult), 42,245 typed accessible elements (Figure 2 source "
        "data 1, ce11 coordinates), the source's developmental-dynamics call per element (Figure 4 "
        "source data 1, devel_is_dynamic) and per-stage ATAC heights (Figure 1 source data 1); "
        "strand-specific long nuclear RNA (chromatin-associated, over 200 nt) at the same six stages, "
        "GEO GSE114483, linear coverage normalised by DESeq2 size factors and averaged over two "
        "biological replicates. Both measured by the same study at the same stages; the source does not "
        "say whether ATAC and RNA came from the same nuclei preparations."
    ),
    "independence": (
        "The RNA (stage, whole animal) and the element selection (developmentally dynamic, not "
        "tissue-specific) are new; the elements are not: Serizay 2020 extends the same lab's annotation, "
        "so many elements here are elements of stage 2's population and carry the same constrained "
        "fraction. The census counts the overlap; the elements that overlap none of stage 2's 5,814 are "
        "reported beside, and a DOSE_ADDS reading carries a registered caveat if that subset's interval "
        "does not exclude zero."
    ),
    "population": (
        "Putative enhancers (the source's 'putative_enhancer' class) that the source calls "
        "developmentally dynamic (devel_is_dynamic True), on a chromosome with protein-coding genes, "
        "for which the domain rule names a gene. Ageing (glp-1) samples are not used."
    ),
    "adaptations": {
        "group": (
            "Stage 1's group was a tissue, called by the source ('accessible in exactly one of five "
            "tissues'). This source gives no per-stage call, so the element's group is the stage of its "
            "highest ATAC height among the six wild-type stages (atac_wt_<stage>_height, SPMR-normalised "
            "peak signal), among the elements the source calls dynamic; a tie goes to the earlier stage."
        ),
        "rna": (
            "Stage 1 read the source's TPM. This source publishes no gene table, so the gene's value at a "
            "stage is coverage-TPM: the mean of the long-cap linear coverage over the union of the "
            "exons of all the gene's WBcel235 transcripts, on the gene's own strand, scaled per stage "
            "so that the protein-coding genes of the six nuclear chromosomes sum to one million. "
            "Coverage over a gene that overlaps another on the same strand is counted for both."
        ),
        "endpoint": (
            "Concordant: at least 5 coverage-TPM at the element's stage and strictly higher than at "
            "each of the other five stages; a gene with no value is not concordant (stage 1's rule with "
            "six stages in place of five tissues)."
        ),
        "strand_check": (
            "Before scoring, per stage: summed own-strand mean coverage over the protein-coding genes "
            "against summed opposite-strand. Sense above antisense at all six stages: the tracks are "
            "read as labelled (fwd = plus strand). Antisense above sense at all six: fwd and rev are "
            "swapped. Anything else: INSTRUMENT_FAILED, nothing scored."
        ),
    },
    "constraint": (
        "Stage 2's instrument, unchanged: UCSC ce11 phyloP135way read by HTTP range request, stage 2's "
        "cut 2.08 (re-derived there at 5% FDR by Zoonomia's procedure; its CDS-above-intron instrument "
        "check passed), exon bases of any WBcel235 transcript removed from each element first; per "
        "element the scored bases left (measured) and those at or above the cut (above). Constrained "
        "fraction = above / measured; elements with no scored base left are not measured and excluded."
    ),
    "dose": (
        "Top third against bottom third of the constrained fraction, at tertile edges computed over "
        "this population's measured elements in the census, before any RNA value is read: top = "
        "fraction at or above the upper edge and above the lower edge, bottom = below the lower edge, "
        "the middle left out (simpler_tiers.by_fraction_tertile, the stage-2 view's own rule)."
    ),
    "matching": (
        "Strata: stage (6) x inside a gene body or between genes (2) x distance bin of the domain-rule "
        "gene's 5' end (stage 1's 8 bins) x measured-length tertile (edges fixed in the census). Only "
        "strata holding both arms count. The matched control of the dose test is the bottom third "
        "within the same strata."
    ),
    "statistic": (
        "delta = concordance(top third) - concordance(bottom third), standardised to the top third's "
        "strata, in points; 95% percentile interval over 1-Mb windows, 2,000 resamples, seed 20261010 "
        "(simpler_tiers.matched_difference, unchanged)."
    ),
    "falsifier": (
        "The dose claim is falsified if delta is not above its within-stratum permutation null: the "
        "arm label shuffled within each stratum 1,000 times (seed 20261011), one-sided p = (1 + "
        "permutations with delta >= observed) / 1,001; FALSIFIED if p >= 0.05."
    ),
    "decision": (
        "DOSE_ADDS if the interval's lower bound > 0; DOSE_BELOW if its upper bound < 0; "
        "NO_DIFFERENCE_DETECTED otherwise."
    ),
    "positive_control": (
        "Developmentally dynamic coding promoters (the source's 'coding_promoter' class), stage by the "
        "same rule, nearest-TSS gene: concordance against a gene drawn from the same chromosome at the "
        "same 5'-end distance bin (stage 1's matched control, 1,000 draws). If one-sided p >= 0.05 the "
        "readout is declared unresponsive (UNRESPONSIVE) and neither reading is given."
    ),
    "secondary": (
        "Stage 1's comparison on this source, with its own label, deciding nothing about the dose: on "
        "the whole population (measured or not), the domain rule against nearest TSS, stage 1's "
        "statistic (difference in concordance, 95% cluster bootstrap over 1-Mb windows, 2,000 "
        "resamples); DOMAIN_RULE_AHEAD, NEAREST_AHEAD or NO_DIFFERENCE_DETECTED; matched control and "
        "falsifier stage 1's (the domain rule's gene against a gene at the same distance bin; falsified "
        "as a carrier of stage information if one-sided p >= 0.05). Fewer than 200 elements where the "
        "rules differ: the secondary is not scored."
    ),
    "beside": (
        "Deciding nothing: stage 2's own numeric edges (fraction >= 0.190476 against < 0.08); the "
        "elements overlapping none of stage 2's population; inside gene bodies and between genes; each "
        "stage; the middle third against the bottom third; stage 2's registered 5% tier on this "
        "source; the secondary's inside and between strata."
    ),
    "power": (
        "Stated in the census before the run: 2.8 x sqrt(0.25 x (1/n_top + 1/n_bottom)) in support, "
        "before clustering."
    ),
    "confound": (
        "Stages 1 and 2 named transcription: an element inside a gene body can be open when its host "
        "is transcribed, an element upstream of an annotated 5' end can lie in a transcribed outron, and "
        "transcribed sequence can be constrained for splicing or RNA-level reasons. Whole-animal RNA "
        "adds one more: a stage concordance can come from any tissue. A difference therefore does not "
        "separate regulatory constraint from transcription-linked constraint."
    ),
    "stop_rule": (
        "One run. Fewer than 200 elements in the top or the bottom third in strata holding both: stop "
        "unscored (UNDERPOWERED). Nothing (source, population, stage rule, RNA construction, cut, edges, "
        "strata, endpoint) changes after any RNA value is read; another definition is a new "
        "registration committed before it is run."
    ),
    "outcome_wording": {
        "common": (
            "On C. elegans putative enhancers whose accessibility changes across development (Jänes et "
            "al. 2018, whole animals, six stages from embryo to young adult), each assigned to its stage "
            "of highest accessibility, with constraint measured on their exon-free bases (UCSC ce11 "
            "phyloP135way, stage 2's cut 2.08), the domain-rule gene is expressed at 5 coverage-TPM or "
            "more in long nuclear RNA and highest at the element's own stage for {kc} of {nc} elements in "
            "the top third by constrained fraction ({rc}), against {ru} for the {nu} in the bottom third "
            "at the same stage, position class, distance and measured length; difference {d} points, "
            "95% {lo} to {hi} over 1-Mb clusters (within-stratum permutation p {perm}; {falsifier})."
        ),
        "DOSE_ADDS": (
            "The dose view seen beside stage 2 replicates on a second source: on this endpoint the "
            "elements with the largest constrained fraction carry stage information beyond the domain "
            "rule's position that those with the smallest do not."
        ),
        "NO_DIFFERENCE_DETECTED": (
            "No difference was detected: the dose view seen beside stage 2 does not replicate on this "
            "source. This is not evidence that constrained elements lack function."
        ),
        "DOSE_BELOW": (
            "The elements with the largest constrained fraction have genes less often stage-concordant "
            "than matched elements with the smallest, the opposite of the dose view seen beside stage 2; "
            "one reading named in advance is that constrained elements sit near genes expressed across "
            "stages, which an endpoint asking for the highest stage counts against."
        ),
        "overlap_caveat": (
            "On the elements that overlap none of stage 2's population the difference is {fresh}, an "
            "interval that does not exclude zero, so the replication rests partly on elements stage 2 "
            "already held."
        ),
        "positive_control": (
            "The positive control separates: of {pn} developmentally dynamic coding promoters, the "
            "nearest-TSS gene is stage-concordant for {pk} ({pr}) against {pc} at a matched distance "
            "(one-sided p {pp})."
        ),
        "secondary_common": (
            "Secondary, stage 1's comparison on this source: on {n} elements the domain-rule gene is "
            "stage-concordant for {kd} ({rd}) against {kn} ({rn}) for the gene with the nearest 5' end; "
            "the rules differ on {k}; difference {d} points, 95% {lo} to {hi} over 1-Mb clusters."
        ),
        "secondary_DOMAIN_RULE_AHEAD": (
            "On this endpoint the domain rule is ahead of nearest TSS, as in stage 1."
        ),
        "secondary_NO_DIFFERENCE_DETECTED": (
            "No difference was detected between the domain rule and nearest TSS."
        ),
        "secondary_NEAREST_AHEAD": (
            "On this endpoint nearest TSS is ahead of the domain rule, the opposite of stage 1."
        ),
        "secondary_control_passed": (
            "The domain rule's genes match the element's stage more often than genes at a matched "
            "distance ({cr}, one-sided p {cp})."
        ),
        "secondary_control_failed": (
            "The domain rule's genes do not match the element's stage more often than genes at a matched "
            "distance ({cr}, one-sided p {cp}): falsified as a carrier of stage information here."
        ),
        "secondary_underpowered": (
            "Secondary: fewer than 200 elements where the domain rule and nearest TSS differ; not scored."
        ),
        "always": (
            "Conservation and accessibility are not function: this says the attributed gene is "
            "expressed when the element is open, not that the element acts on it, and it does not "
            "separate regulatory constraint from transcription-linked constraint. The RNA is from whole "
            "animals, so a stage concordance can come from any tissue."
        ),
        "INSTRUMENT_FAILED": (
            "INSTRUMENT_FAILED: the long-cap RNA tracks show neither own-strand coverage above "
            "opposite-strand at all six stages nor the reverse at all six, so the RNA's strand cannot be "
            "read; nothing is scored."
        ),
        "UNRESPONSIVE": (
            "UNRESPONSIVE: of {pn} developmentally dynamic coding promoters, the nearest-TSS gene is "
            "stage-concordant for {pk} ({pr}) against {pc} at a matched distance (one-sided p {pp}); the "
            "readout does not separate, so neither reading is given."
        ),
        "UNDERPOWERED": (
            "UNDERPOWERED: fewer than 200 elements in the top or the bottom third in strata holding "
            "both; stopped unscored."
        ),
    },
}
