# SPDX-License-Identifier: AGPL-3.0-or-later
"""Measured grammar: which bases of a regulatory element change its activity, and what marks them.

Area J's panel (knowledge/across.py) showed that motif sites held across species do not separate VISTA
enhancers from conserved negatives, and could not say why: nothing in it measured which bases matter.
Kircher et al. 2019 (Nature Communications 10:3583) did, for 21 regulatory elements: saturation
mutagenesis read out by a massively parallel reporter assay, nearly every single-base substitution of
each element with its effect on activity (a log2 coefficient) and a p-value, 44,658 measurements on
GRCh38 (data: github.com/kircherlab/MPRA_SaturationMutagenesis, data/elements.tsv.gz; GEO GSE126550).

With that ground truth this module asks, per element and pooled:
- conservation: does mammalian phyloP (Zoonomia, 241 species) rank the bases that matter above the rest?
- motif sites: how many of the bases that matter fall inside a JASPAR site, against the share of all bases
  sites cover, across a sweep of match thresholds (calibration against measurement, not a search for a
  passing test);
- both: are conserved bases inside a site more often functional than conserved bases outside one;
- families: which TFClass families' sites contain a functional base more often than other sites.

A base is functional when some substitution at it is significant (p below 1e-5 with at least 10
barcodes, the data portal's defaults), and strongly functional when that substitution also changes
activity by |log2| >= 0.25 (about 19%): at deep coverage a p-value alone flags small effects, so the share of
"functional" bases swings with depth (1% of FOXE1, 77% of IRF4), and both definitions are reported.
Deletions are left out, being a different kind of change. Besides the threshold sweep, each base gets a
continuous motif score, the best relative score of any JASPAR hit covering it, whose AUC needs no cut. Several
loci were measured more than once (other cells or time points); pooled statistics use one primary
experiment per locus, the one with the most measurements, and repeats are used only to report agreement.
Evidence: the measurement is `experimental`; motif hits are `predicted`; the enrichments are `inferred`.
"""

from __future__ import annotations

import csv
import gzip
import time
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Any

from genomeos.results import RESULTS_DIR, save_result

DATA_URL = (
    "https://raw.githubusercontent.com/kircherlab/MPRA_SaturationMutagenesis/master/data/elements.tsv.gz"
)
DATA_PATH = Path("data/knowledge/satmut/elements.tsv.gz")
MIN_BARCODES = 10
P_SIGNIFICANT = 1e-5
THRESHOLDS = (0.80, 0.85, 0.90, 0.95)
SITE_THRESHOLD = (
    0.95  # below it JASPAR sites cover 97% or more of an element's bases (measured: see the sweep)
)
STRONG_EFFECT = 0.25  # |log2 coefficient| for the "strong" definition, about a 19% change in activity
PHYLOP_CONSTRAINED = 2.27
BASES = set("ACGT")
EVIDENCE = {
    "measurement": "experimental: saturation mutagenesis MPRA, Kircher et al. 2019 (GSE126550), GRCh38",
    "functional": f"a base with a substitution at p < {P_SIGNIFICANT} and >= {MIN_BARCODES} barcodes",
    "motifs": "predicted: JASPAR 2026 CORE vertebrates, every hit above the threshold, TFClass families",
    "conservation": "curated: Zoonomia phyloP over 241 placental mammals, per base",
    "tests": "inferred: one-sided Fisher and Mann–Whitney, Benjamini–Hochberg over families",
}


def load(path: Path = DATA_PATH, release: str = "GRCh38") -> list[dict]:
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(DATA_URL, headers={"User-Agent": "GenomeOS/0.9 (satmut)"})
        with urllib.request.urlopen(req, timeout=300) as r:  # noqa: S310
            path.write_bytes(r.read())
    with gzip.open(path, "rt") as fh:
        rows = []
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["Release"] != release or r["Pos"] in ("", "NA") or r["Alt"] in ("", "NA"):
                continue
            rows.append(
                {
                    "element": r["Element"],
                    "chrom": f"chr{r['Chrom']}",
                    "pos": int(r["Pos"]),
                    "ref": r["Ref"],
                    "alt": r["Alt"],
                    "barcodes": int(r["Barcodes"]),
                    "coef": float(r["Coefficient"]),
                    "p": float(r["pValue"]),
                }
            )
    return rows


def by_element(rows: list[dict]) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        out[r["element"]].append(r)
    return dict(out)


def loci(experiments: dict[str, list[dict]]) -> dict[str, list[str]]:
    """Experiments grouped by the locus they measured (overlapping spans on one chromosome).

    The first experiment of each group is the primary one: the most measurements, then the name.
    """
    spans = {
        e: (rs[0]["chrom"], min(r["pos"] for r in rs), max(r["pos"] for r in rs))
        for e, rs in experiments.items()
    }
    groups: list[list[str]] = []
    for e in sorted(spans, key=lambda k: (spans[k][0], spans[k][1])):
        c, lo, hi = spans[e]
        for g in groups:
            gc, glo, ghi = spans[g[0]]
            if gc == c and lo <= ghi and glo <= hi:
                g.append(e)
                break
        else:
            groups.append([e])
    out = {}
    for g in groups:
        g.sort(key=lambda e: (-len(experiments[e]), e))
        out[g[0]] = g
    return out


def base_table(
    rows: list[dict], min_barcodes: int = MIN_BARCODES, p_cut: float = P_SIGNIFICANT
) -> dict[int, dict]:
    """Per position: how many substitutions were measured, whether one is significant, the largest effect."""
    out: dict[int, dict] = {}
    for r in rows:
        if r["alt"] not in BASES:
            continue
        b = out.setdefault(
            r["pos"], {"ref": r["ref"], "measured": 0, "functional": False, "strong": False, "effect": 0.0}
        )
        if r["barcodes"] < min_barcodes:
            continue
        b["measured"] += 1
        if r["p"] < p_cut and abs(r["coef"]) > abs(b["effect"]):
            b["functional"] = True
            b["effect"] = round(r["coef"], 4)
            b["strong"] = abs(r["coef"]) >= STRONG_EFFECT
    return {p: b for p, b in out.items() if b["measured"]}


def auc(positive: list[float], negative: list[float]) -> float | None:
    """Probability that a positive scores above a negative (ties count half): the area under the ROC curve."""
    if not positive or not negative:
        return None
    allv = sorted([(v, 1) for v in positive] + [(v, 0) for v in negative])
    rank_sum, i = 0.0, 0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1][0] == allv[i][0]:
            j += 1
        mean_rank = (i + j) / 2 + 1
        rank_sum += mean_rank * sum(1 for k in range(i, j + 1) if allv[k][1] == 1)
        i = j + 1
    n1, n0 = len(positive), len(negative)
    return round((rank_sum - n1 * (n1 + 1) / 2) / (n1 * n0), 4)


def site_coverage(
    hits: list[tuple[str, int, int, int, float]], offset: int, threshold: float
) -> dict[int, set[str]]:
    """Genome position (1-based) to the set of factors whose site at or above `threshold` covers it."""
    cov: dict[int, set[str]] = defaultdict(set)
    for f, st, en, _, rel in hits:
        if rel < threshold:
            continue
        for k in range(st, en):
            cov[offset + k + 1].add(f)
    return cov


def motif_score(hits: list[tuple[str, int, int, int, float]], offset: int) -> dict[int, float]:
    """Genome position (1-based) to the best relative score of any hit covering it."""
    best: dict[int, float] = {}
    for _f, st, en, _sd, rel in hits:
        for k in range(st, en):
            pos = offset + k + 1
            if rel > best.get(pos, 0.0):
                best[pos] = rel
    return best


def sweep(
    bases: dict[int, dict], coverage: dict[float, dict[int, set[str]]], key: str = "functional"
) -> list[dict]:
    """For each threshold: share of bases in a site, recall and precision of sites for functional bases."""
    from genomeos.knowledge.across import fisher_greater

    n = len(bases)
    fn = sum(1 for b in bases.values() if b[key])
    rows = []
    for t, cov in sorted(coverage.items()):
        inside = [p for p in bases if p in cov]
        a = sum(1 for p in inside if bases[p][key])  # functional, in a site
        b = len(inside) - a
        c = fn - a
        d = n - len(inside) - c
        rows.append(
            {
                "threshold": t,
                "bases_in_sites": round(len(inside) / n, 4) if n else None,
                "recall": round(a / fn, 4) if fn else None,
                "precision": round(a / len(inside), 4) if inside else None,
                "baseline": round(fn / n, 4) if n else None,
                "enrichment": round((a / len(inside)) / (fn / n), 3) if inside and fn and n else None,
                "p": fisher_greater(a, b, c, d) if n else None,
                "table": [a, b, c, d],
            }
        )
    return rows


def analyse_locus(
    name: str, rows: list[dict], motifs, families: dict, phylop: dict[int, float], sequence: str, start0: int
) -> dict[str, Any]:
    """One experiment: functional bases (both definitions) against conservation, motif sites and families."""
    from genomeos.genome.motifs import all_hits, family_unit

    bases = base_table(rows)
    hits = all_hits(motifs, sequence, min_relative=min(THRESHOLDS))
    coverage = {t: site_coverage(hits, start0, t) for t in THRESHOLDS}
    mscore = motif_score(hits, start0)
    cov = coverage[SITE_THRESHOLD]
    out: dict[str, Any] = {
        "experiment": name,
        "chrom": rows[0]["chrom"],
        "start": min(bases),
        "end": max(bases),
        "bases_measured": len(bases),
        "hits_at_lowest_threshold": len(hits),
    }
    for key in ("functional", "strong"):
        fun = [p for p, b in bases.items() if b[key]]
        non = [p for p, b in bases.items() if not b[key]]
        conserved = [p for p in bases if phylop.get(p, float("-inf")) >= PHYLOP_CONSTRAINED]
        cons_in = [p for p in conserved if p in cov]
        cons_out = [p for p in conserved if p not in cov]
        out[key] = {
            "bases": len(fun),
            "share": round(len(fun) / len(bases), 4) if bases else None,
            "phylop_auc": auc(
                [phylop[p] for p in fun if p in phylop], [phylop[p] for p in non if p in phylop]
            ),
            "motif_auc": auc([mscore.get(p, 0.0) for p in fun], [mscore.get(p, 0.0) for p in non]),
            "sweep": sweep(bases, coverage, key),
            "conserved_counts": {
                "in_site": [sum(bases[p][key] for p in cons_in), len(cons_in)],
                "outside": [sum(bases[p][key] for p in cons_out), len(cons_out)],
            },
        }
    sites = []
    for f, st, en, _sd, rel in hits:
        if rel < SITE_THRESHOLD:
            continue
        span = [start0 + k + 1 for k in range(st, en)]
        measured = [p for p in span if p in bases]
        if not measured:
            continue
        sites.append(
            {
                "factor": f,
                "family": family_unit(f, families),
                "start": span[0],
                "end": span[-1],
                "functional": any(bases[p]["functional"] for p in measured),
                "strong": any(bases[p]["strong"] for p in measured),
            }
        )
    out["sites"] = sites
    return out


def replicate_agreement(experiments: dict[str, list[dict]], group: list[str]) -> list[dict]:
    """Pearson correlation of substitution effects measured in two experiments of the same locus."""
    from genomeos.knowledge.profiling import pearson

    out = []
    first = {(r["pos"], r["alt"]): r["coef"] for r in experiments[group[0]] if r["alt"] in BASES}
    for other in group[1:]:
        second = {(r["pos"], r["alt"]): r["coef"] for r in experiments[other] if r["alt"] in BASES}
        shared = sorted(set(first) & set(second))
        out.append(
            {
                "experiments": [group[0], other],
                "shared_substitutions": len(shared),
                "r": pearson([first[k] for k in shared], [second[k] for k in shared])
                if len(shared) > 2
                else None,
            }
        )
    return out


def pool(results: list[dict], key: str = "functional") -> dict[str, Any]:
    """Pooled over primary experiments for one definition: sweep, AUCs, the conserved contrast, families."""
    from genomeos.knowledge.across import bh, fisher_greater

    pooled_sweep = []
    for i, t in enumerate(THRESHOLDS):
        a = b = c = d = 0
        for r in results:
            ta, tb, tc, td = r[key]["sweep"][i]["table"]
            a, b, c, d = a + ta, b + tb, c + tc, d + td
        n, fn, inside = a + b + c + d, a + c, a + b
        pooled_sweep.append(
            {
                "threshold": t,
                "bases_in_sites": round(inside / n, 4) if n else None,
                "recall": round(a / fn, 4) if fn else None,
                "precision": round(a / inside, 4) if inside else None,
                "baseline": round(fn / n, 4) if n else None,
                "enrichment": round((a / inside) / (fn / n), 3) if inside and fn and n else None,
                "p": fisher_greater(a, b, c, d),
            }
        )
    ci = [sum(r[key]["conserved_counts"]["in_site"][k] for r in results) for k in (0, 1)]
    co = [sum(r[key]["conserved_counts"]["outside"][k] for r in results) for k in (0, 1)]
    conserved = {
        "in_site": {"functional": ci[0], "bases": ci[1], "share": round(ci[0] / ci[1], 4) if ci[1] else None},
        "outside": {"functional": co[0], "bases": co[1], "share": round(co[0] / co[1], 4) if co[1] else None},
        "p": fisher_greater(ci[0], ci[1] - ci[0], co[0], co[1] - co[0]) if ci[1] and co[1] else None,
    }
    fam_sites: dict[str, list[bool]] = defaultdict(list)
    for r in results:
        seen: set[tuple[str, int]] = set()
        for s_ in r["sites"]:
            k = (s_["family"], s_["start"])  # one family's overlapping matrices at one start count once
            if k in seen:
                continue
            seen.add(k)
            fam_sites[s_["family"]].append(s_[key])
    total = sum(len(v) for v in fam_sites.values())
    total_f = sum(sum(v) for v in fam_sites.values())
    fam_rows = []
    for fam, v in fam_sites.items():
        a, n = sum(v), len(v)
        if n < 5:
            continue
        c, m = total_f - a, total - n
        fam_rows.append(
            {
                "family": fam,
                "sites": n,
                "functional_sites": a,
                "share": round(a / n, 3),
                "other_share": round(c / m, 3) if m else None,
                "p": fisher_greater(a, n - a, c, m - c),
            }
        )
    for row, q in zip(fam_rows, bh([r["p"] for r in fam_rows]), strict=True):
        row["q"] = round(q, 6)
    fam_rows.sort(key=lambda r: r["p"])

    def med(vals: list[float]) -> float | None:
        vals = sorted(v for v in vals if v is not None)
        return vals[len(vals) // 2] if vals else None

    ph = [r[key]["phylop_auc"] for r in results]
    mo = [r[key]["motif_auc"] for r in results]
    return {
        "definition": key,
        "loci": len(results),
        "bases_measured": sum(r["bases_measured"] for r in results),
        "bases_functional": sum(r[key]["bases"] for r in results),
        "phylop_auc_median": med(ph),
        "phylop_auc_above_half": sum(1 for x in ph if x is not None and x > 0.5),
        "motif_auc_median": med(mo),
        "motif_auc_above_half": sum(1 for x in mo if x is not None and x > 0.5),
        "sweep": pooled_sweep,
        "conserved": conserved,
        "families_tested": len(fam_rows),
        "families_q05": [r for r in fam_rows if r["q"] <= 0.05],
        "families": fam_rows[:30],
    }


def run(results_dir: Path = RESULTS_DIR, progress=None) -> dict[str, Any]:
    from genomeos.attribution.bigwig import BigWig
    from genomeos.attribution.constraint import PHYLOP_241_URL
    from genomeos.coords import Locus
    from genomeos.genome.fetch import REFERENCE
    from genomeos.genome.genome import Genome
    from genomeos.genome.motifs import load_families, load_motifs

    t0 = time.time()
    experiments = by_element(load())
    groups = loci(experiments)
    motifs = load_motifs(relative=min(THRESHOLDS))
    families = load_families()
    genomes: dict[str, Genome] = {}
    per_experiment: list[dict] = []
    primary: list[dict] = []
    agreement: list[dict] = []
    bw = BigWig(PHYLOP_241_URL)
    try:
        for lead, group in groups.items():
            for name in group:
                rows = experiments[name]
                chrom = rows[0]["chrom"]
                lo, hi = min(r["pos"] for r in rows), max(r["pos"] for r in rows)
                if chrom not in genomes:
                    genomes[chrom] = Genome.from_fasta(REFERENCE / f"{chrom}.fa.gz")
                seq = str(genomes[chrom].fetch(Locus(chrom, lo - 1, hi)))
                ivs = [(p - 1, p) for p in range(lo, hi + 1)]
                stats = bw.summarise(chrom, ivs, PHYLOP_CONSTRAINED)
                phylop = {p: st.maximum for p, st in zip(range(lo, hi + 1), stats, strict=True) if st.bases}
                res = analyse_locus(name, rows, motifs, families, phylop, seq, lo - 1)
                res["locus"] = lead
                res["primary"] = name == lead
                per_experiment.append(res)
                if name == lead:
                    primary.append(res)
                if progress:
                    fx, st_ = res["functional"], res["strong"]
                    progress(
                        f"{name}: {fx['bases']}/{res['bases_measured']} functional ({st_['bases']} strong); "
                        f"AUC phyloP {fx['phylop_auc']} motif {fx['motif_auc']}"
                    )
            agreement += replicate_agreement(experiments, group)
    finally:
        bw.close()
    return {
        "experiments": len(experiments),
        "loci_count": len(groups),
        "loci": {k: v for k, v in groups.items()},
        "pooled": {"functional": pool(primary, "functional"), "strong": pool(primary, "strong")},
        "replicate_agreement": agreement,
        "per_experiment": per_experiment,
        "thresholds": {
            "min_barcodes": MIN_BARCODES,
            "p_significant": P_SIGNIFICANT,
            "motif_thresholds": list(THRESHOLDS),
            "site_threshold": SITE_THRESHOLD,
            "strong_effect_log2": STRONG_EFFECT,
            "phylop_constrained": PHYLOP_CONSTRAINED,
        },
        "evidence": EVIDENCE,
        "source": DATA_URL,
        "cost": {"seconds": round(time.time() - t0, 1)},
    }


def run_and_save(results_dir: Path = RESULTS_DIR, progress=None) -> dict[str, Any]:
    out = run(results_dir, progress)
    save_result("satmut_grammar", out, results_dir)
    return out
