# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Simpler genomes, C. elegans first: one positional attribution held against measured tissue RNA.

Area I's direction of 2026-09-11 named three routes to function for the non-coding genome:
attribution, closure tests and simpler genomes. This module is the first step of the third. The
worm has about 100 Mb, a gene every 5 kb, and one public table in which the same five tissues were
measured for both chromatin accessibility and RNA: Serizay et al. 2020, Genome Research 30:1752
(young adults; germline, neurons, muscle, hypodermis, intestine; Supplemental Table S2, CC BY 4.0).

WHAT TRANSFERS FROM THE HUMAN PIPELINE WITHOUT A MODEL. The human element-to-gene step is a
deletion in AlphaGenome, which this module never calls (no model request, and AlphaGenome is human
and mouse only). What is left of the human attribution without the model is the node model: an
element acts inside its node and on the nearest gene there. The human node is a CTCF domain, and
C. elegans has no CTCF, so the node cannot be carried over as it is. Its worm instance here is
positional and reads only the genome annotation (WBcel235, Ensembl gene models):

    domain rule   an element inside a protein-coding gene is attributed to that gene (the host;
                  the one with the nearest 5' end if genes overlap). An element between genes is
                  attributed to the flanking gene whose 5' end faces it, i.e. the gene it sits
                  upstream of. If both flanking genes face it (a divergent pair) the nearer 5'
                  end wins; if neither does (a convergent pair, between two 3' ends) the nearer
                  5' end of the two wins, and that case is counted apart.

The baseline is the nearest annotated 5' end (nearest TSS) among all protein-coding genes. The
two rules differ only where an element is nearer some other gene's 5' end than the 5' end of the
gene it sits in or upstream of; every comparison below is decided on those elements.

THE ENDPOINT THE ATTRIBUTION NEVER USED. An element accessible in exactly one of the five tissues
(the source's own tissue call) is concordant with its attributed gene when that gene's measured
RNA in the same tissue is at least 5 TPM (the source's detection threshold) and higher than in
each of the other four. The attribution reads positions and strands only: no accessibility, no
RNA, and none of the source's own element-to-gene columns (`associated_gene_id`, the promoter
gene columns, `nb_associated_enhancers`).

WHAT THIS CANNOT ESTABLISH. Accessibility is not function: the elements are defined by measured
accessibility, and concordance says that the attributed gene is expressed where the element is
open, not that the element acts on it. An element inside a gene body can be open where its host
is transcribed for that reason alone, so the host case favours the domain rule by construction of
the assay; the intergenic elements are reported apart for that reason. The annotated 5' end is
not always the transcription start: most worm transcripts are trans-spliced, and the true start
can sit kilobases upstream in an outron.
"""

from __future__ import annotations

import bisect
import gzip
import hashlib
import math
import random
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

#: The source table and the genome annotation, both git-ignored local copies.
STORE = Path("data/knowledge/simpler_worm")
TABLE_S2 = STORE / "supp_gr.265934.120_Supplemental_Table_S2.xlsx"
GFF3 = Path("data/reference/celegans/WBcel235.63.gff3.gz")
#: Where the table was fetched: Europe PMC's bundle of the article's supplementary files. The
#: journal's own link redirects to a login, and PMC answers a scripted request with a challenge.
TABLE_S2_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7706728/supplementaryFiles"
TABLE_S2_MEMBER = "supp_gr.265934.120_Supplemental_Table_S2.xlsx"
SOURCE = (
    "Serizay J, Dong Y, Janes J, Chesney M, Cerrato C, Ahringer J (2020) Distinctive regulatory "
    "architectures of germline-active and somatic genes in C. elegans. Genome Res 30:1752-1765, "
    "Supplemental Table S2 (sheets ATAC_metrics and RNA_metrics); CC BY 4.0"
)

#: The five tissues, in the column order both sheets use.
TISSUES = ("Germline", "Neurons", "Muscle", "Hypod.", "Intest.")
#: The element population: the source's distal class (initiation signal, no elongation).
ELEMENT_TYPE = "putative_enhancer"
#: The positive control's population: the source's coding promoters.
PROMOTER_TYPE = "coding_promoter"
#: The source's detection threshold for a gene, in TPM.
DETECTION_TPM = 5.0
#: Distance bins for the matched control (bp, upper edges; the last bin is open).
DISTANCE_BINS = (1_000, 2_000, 5_000, 10_000, 20_000, 50_000, 100_000)
CONTROL_DRAWS = 1_000
BOOTSTRAP = 2_000
#: Elements are resampled in 1-Mb windows, so neighbouring elements that share genes move together.
CLUSTER_BP = 1_000_000
SEED = 20261010
#: Below this many elements on which the two rules name different genes, the run stops unscored.
MIN_DISCORDANT = 200
ALPHA = 0.05

#: Sheet columns (0-based) read from Table S2.
ATAC_COLS = {"chrom": 0, "start": 1, "end": 2, "type": 3, "rpm": (17, 22), "annotation": 31}
RNA_COLS = {"wbid": 4, "tpm": (8, 13)}
#: Columns of Table S2 that carry the source's own element-to-gene links. Never read.
NEVER_READ = (
    "promoter_gene_id_fwd",
    "promoter_locus_id_fwd",
    "promoter_gene_id_rev",
    "promoter_locus_id_rev",
    "associated_gene_id",
    "associated_locus_id",
    "nb_associated_promoters",
    "nb_associated_enhancers",
)


@dataclass(frozen=True)
class Gene:
    gid: str
    name: str
    chrom: str
    start: int
    end: int
    strand: str

    @property
    def tss(self) -> int:
        return self.start if self.strand == "+" else self.end


@dataclass(frozen=True)
class Element:
    chrom: str
    start: int
    end: int
    kind: str
    tissue: str

    @property
    def mid(self) -> int:
        return (self.start + self.end) // 2


class Chromosome:
    """The protein-coding genes of one chromosome, indexed by start and by 5' end."""

    def __init__(self, genes: list[Gene]):
        self.genes = sorted(genes, key=lambda g: (g.start, g.end))
        self.starts = [g.start for g in self.genes]
        self.max_len = max((g.end - g.start for g in self.genes), default=0)
        self.by_end = sorted(self.genes, key=lambda g: g.end)
        self.ends = [g.end for g in self.by_end]
        self.by_tss = sorted(self.genes, key=lambda g: g.tss)
        self.tss = [g.tss for g in self.by_tss]

    def hosts(self, pos: int) -> list[Gene]:
        lo = bisect.bisect_left(self.starts, pos - self.max_len)
        hi = bisect.bisect_right(self.starts, pos)
        return [g for g in self.genes[lo:hi] if g.start <= pos <= g.end]

    def left(self, pos: int) -> Gene | None:
        i = bisect.bisect_left(self.ends, pos)
        return self.by_end[i - 1] if i > 0 else None

    def right(self, pos: int) -> Gene | None:
        i = bisect.bisect_right(self.starts, pos)
        return self.genes[i] if i < len(self.genes) else None

    def nearest_tss(self, pos: int) -> Gene | None:
        if not self.tss:
            return None
        i = bisect.bisect_left(self.tss, pos)
        cands = [self.by_tss[j] for j in (i - 1, i) if 0 <= j < len(self.tss)]
        return min(cands, key=lambda g: (abs(g.tss - pos), g.gid))


def _nearer(pos: int, genes: list[Gene]) -> Gene:
    return min(genes, key=lambda g: (abs(g.tss - pos), g.gid))


def domain_rule(chrom: Chromosome, pos: int) -> tuple[Gene | None, str]:
    """The worm node: host gene, else the flanking gene the element sits upstream of."""
    hosts = chrom.hosts(pos)
    if hosts:
        return _nearer(pos, hosts), "host"
    left, right = chrom.left(pos), chrom.right(pos)
    flank = [g for g in (left, right) if g is not None]
    if not flank:
        return None, "none"
    if len(flank) == 1:
        return flank[0], "one_flank"
    facing = []
    if left is not None and left.strand == "-":
        facing.append(left)
    if right is not None and right.strand == "+":
        facing.append(right)
    if len(facing) == 1:
        return facing[0], "faces_one"
    if len(facing) == 2:
        return _nearer(pos, facing), "divergent"
    return _nearer(pos, flank), "convergent_fallback"


def distance_bin(d: int) -> int:
    return bisect.bisect_right(DISTANCE_BINS, abs(d))


# --------------------------------------------------------------------------------------------- reading


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_genes(path: Path = GFF3) -> dict[str, Chromosome]:
    """Protein-coding genes of the six nuclear chromosomes, named the way Table S2 names them."""
    per: dict[str, list[Gene]] = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "gene" or "biotype=protein_coding" not in p[8]:
                continue
            if p[0] == "MtDNA":
                continue
            attrs = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
            gid = attrs.get("gene_id") or attrs["ID"].split(":")[-1]
            chrom = "chr" + p[0]
            per.setdefault(chrom, []).append(
                Gene(gid, attrs.get("Name", gid), chrom, int(p[3]), int(p[4]), p[6])
            )
    return {c: Chromosome(g) for c, g in per.items()}


def _sheet_rows(path: Path, sheet: str):
    import openpyxl

    wb = openpyxl.load_workbook(path, read_only=True)
    rows = wb[sheet].iter_rows(values_only=True)
    for _ in range(3):  # two comment rows and the header
        next(rows)
    for r in rows:
        if r and r[0] is not None:
            yield r


def read_elements(
    path: Path = TABLE_S2, kinds: tuple[str, ...] = (ELEMENT_TYPE, PROMOTER_TYPE)
) -> list[Element]:
    """Elements of the given classes accessible in exactly one tissue (the source's own call)."""
    out = []
    for r in _sheet_rows(path, "ATAC_metrics"):
        kind, tissue = r[ATAC_COLS["type"]], r[ATAC_COLS["annotation"]]
        if kind in kinds and tissue in TISSUES:
            out.append(Element(r[0], int(r[1]), int(r[2]), kind, tissue))
    return out


def read_expression(path: Path = TABLE_S2) -> dict[str, tuple[float, ...]]:
    """The endpoint: per gene, TPM in the five tissues. Read only by `score`."""
    a, b = RNA_COLS["tpm"]
    out = {}
    for r in _sheet_rows(path, "RNA_metrics"):
        tpm = tuple(float(x) if x is not None else 0.0 for x in r[a:b])
        out[str(r[RNA_COLS["wbid"]])] = tpm
    return out


def concordant(tpm: tuple[float, ...] | None, tissue: str) -> bool:
    """Expressed at >= DETECTION_TPM in the tissue and higher there than in each of the other four."""
    if tpm is None:
        return False
    i = TISSUES.index(tissue)
    v = tpm[i]
    return v >= DETECTION_TPM and all(v > tpm[j] for j in range(len(TISSUES)) if j != i)


# --------------------------------------------------------------------------------------------- pairing


def attribute(elements: list[Element], chroms: dict[str, Chromosome]) -> list[dict[str, Any]]:
    """Each element with the domain rule's gene, its case, and the nearest-TSS gene. No endpoint read."""
    rows = []
    for e in elements:
        ch = chroms.get(e.chrom)
        if ch is None:
            continue
        d, case = domain_rule(ch, e.mid)
        n = ch.nearest_tss(e.mid)
        if d is None or n is None:
            continue
        rows.append(
            {
                "chrom": e.chrom,
                "mid": e.mid,
                "kind": e.kind,
                "tissue": e.tissue,
                "domain_gene": d.gid,
                "domain_case": case,
                "domain_dist": abs(d.tss - e.mid),
                "nearest_gene": n.gid,
                "nearest_dist": abs(n.tss - e.mid),
            }
        )
    return rows


def census(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Population counts and power in advance, from positions and the element-side tissue call only."""
    enh = [r for r in rows if r["kind"] == ELEMENT_TYPE]
    prom = [r for r in rows if r["kind"] == PROMOTER_TYPE]
    disc = [r for r in enh if r["domain_gene"] != r["nearest_gene"]]
    cases: dict[str, int] = {}
    disc_cases: dict[str, int] = {}
    for r in enh:
        cases[r["domain_case"]] = cases.get(r["domain_case"], 0) + 1
    for r in disc:
        disc_cases[r["domain_case"]] = disc_cases.get(r["domain_case"], 0) + 1
    by_tissue = {t: sum(1 for r in enh if r["tissue"] == t) for t in TISSUES}
    disc_by_tissue = {t: sum(1 for r in disc if r["tissue"] == t) for t in TISSUES}
    n, k = len(enh), len(disc)
    clusters = len({(r["chrom"], r["mid"] // CLUSTER_BP) for r in enh})
    # Per-element difference d is in {-1, 0, +1} and is non-zero only where the rules differ, so
    # Var(d) <= k/n and SE(delta) <= sqrt(k)/n; 2.8 SE is the 80%-power, two-sided 5% difference,
    # before any inflation for clustering.
    mde = 2.8 * math.sqrt(k) / n if n else None
    return {
        "elements": n,
        "elements_by_tissue": by_tissue,
        "domain_cases": cases,
        "discordant": k,
        "discordant_fraction": round(k / n, 4) if n else None,
        "discordant_by_case": disc_cases,
        "discordant_by_tissue": disc_by_tissue,
        "host_elements": cases.get("host", 0),
        "intergenic_elements": n - cases.get("host", 0),
        "intergenic_discordant": k - disc_cases.get("host", 0),
        "clusters_1mb": clusters,
        "mde_points_unclustered": round(100 * mde, 3) if mde is not None else None,
        "positive_control_promoters": len(prom),
        "stop_rule_met": k < MIN_DISCORDANT,
    }


# --------------------------------------------------------------------------------------------- scoring


class ControlIndex:
    """Per chromosome and tissue, prefix counts of concordant genes in 5'-end order.

    A matched-control gene is drawn uniformly from the protein-coding genes of the element's
    chromosome whose 5' end lies in the same distance bin as the attributed gene's, the attributed
    gene excluded. Drawing one and asking whether it is concordant is a Bernoulli trial with
    probability f = (concordant candidates) / (candidates), so f is computed exactly from prefix
    counts and the draws are simulated from it.
    """

    def __init__(self, chroms: dict[str, Chromosome], expr: dict[str, tuple[float, ...]]):
        self.chroms = chroms
        self.prefix: dict[tuple[str, str], list[int]] = {}
        for c, ch in chroms.items():
            for t in TISSUES:
                acc, p = 0, [0]
                for g in ch.by_tss:
                    acc += concordant(expr.get(g.gid), t)
                    p.append(acc)
                self.prefix[(c, t)] = p
        self.expr = expr

    def _count(self, chrom: str, tissue: str, lo_pos: float, hi_pos: float) -> tuple[int, int]:
        """Genes with 5' end in [lo_pos, hi_pos): (total, concordant)."""
        tss = self.chroms[chrom].tss
        i, j = bisect.bisect_left(tss, lo_pos), bisect.bisect_left(tss, hi_pos)
        p = self.prefix[(chrom, tissue)]
        return j - i, p[j] - p[i]

    def probability(self, chrom: str, tissue: str, mid: int, gene: str, dist: int) -> float:
        b = distance_bin(dist)
        lo = 0 if b == 0 else DISTANCE_BINS[b - 1]
        hi = DISTANCE_BINS[b] if b < len(DISTANCE_BINS) else math.inf
        # |tss - mid| in [lo, hi): two half-open windows, one each side of the element.
        n1, c1 = self._count(chrom, tissue, mid - hi + 1, mid - lo + 1)
        n2, c2 = self._count(chrom, tissue, mid + lo, mid + hi)
        if lo == 0:  # the element's own position was counted in both windows
            n0, c0 = self._count(chrom, tissue, mid, mid + 1)
            n1, c1 = n1 - n0, c1 - c0
        total, conc = n1 + n2, c1 + c2
        own = concordant(self.expr.get(gene), tissue)
        total, conc = total - 1, conc - int(own)  # the attributed gene is in its own bin
        if total > 0:
            return conc / total
        ch = self.chroms[chrom]
        others = [g for g in ch.genes if g.gid != gene]
        if not others:
            return 0.0
        sub = min(others, key=lambda g: (abs(abs(g.tss - mid) - dist), g.gid))
        return float(concordant(self.expr.get(sub.gid), tissue))


def control_test(
    probs: list[float], observed: int, draws: int = CONTROL_DRAWS, seed: int = SEED
) -> dict[str, Any]:
    rng = random.Random(seed)
    sims = []
    for _ in range(draws):
        sims.append(sum(1 for p in probs if rng.random() < p))
    ge = sum(1 for s in sims if s >= observed)
    n = len(probs)
    return {
        "expected": round(sum(probs), 3),
        "expected_rate": round(sum(probs) / n, 5) if n else None,
        "draws": draws,
        "draws_at_or_above_observed": ge,
        "p_one_sided": round((1 + ge) / (1 + draws), 5),
    }


def cluster_bootstrap(
    rows: list[dict[str, Any]], diff_key: str, reps: int = BOOTSTRAP, seed: int = SEED
) -> dict[str, Any]:
    clusters: dict[tuple[str, int], list[int]] = {}
    for r in rows:
        clusters.setdefault((r["chrom"], r["mid"] // CLUSTER_BP), []).append(r[diff_key])
    keys = sorted(clusters)
    sums = [sum(clusters[k]) for k in keys]
    sizes = [len(clusters[k]) for k in keys]
    n = sum(sizes)
    point = sum(sums) / n if n else float("nan")
    rng = random.Random(seed)
    stats = []
    m = len(keys)
    for _ in range(reps):
        s = z = 0
        for _ in range(m):
            i = rng.randrange(m)
            s += sums[i]
            z += sizes[i]
        stats.append(s / z)
    stats.sort()
    lo = stats[int(0.025 * reps)]
    hi = stats[int(0.975 * reps) - 1]
    return {"point": point, "lo": lo, "hi": hi, "clusters": m, "n": n, "reps": reps}


def _rate(k: int, n: int) -> float | None:
    return round(100 * k / n, 3) if n else None


def compare(rows: list[dict[str, Any]], expr: dict[str, tuple[float, ...]]) -> dict[str, Any]:
    for r in rows:
        r["c_domain"] = int(concordant(expr.get(r["domain_gene"]), r["tissue"]))
        r["c_nearest"] = int(concordant(expr.get(r["nearest_gene"]), r["tissue"]))
        r["d"] = r["c_domain"] - r["c_nearest"]
    n = len(rows)
    kd = sum(r["c_domain"] for r in rows)
    kn = sum(r["c_nearest"] for r in rows)
    bs = cluster_bootstrap(rows, "d") if rows else None
    return {
        "n": n,
        "domain_concordant": kd,
        "nearest_concordant": kn,
        "domain_rate_pct": _rate(kd, n),
        "nearest_rate_pct": _rate(kn, n),
        "delta_points": round(100 * bs["point"], 3) if bs else None,
        "delta_ci95_points": [round(100 * bs["lo"], 3), round(100 * bs["hi"], 3)] if bs else None,
        "clusters": bs["clusters"] if bs else 0,
    }


def comparison_label(lo: float, hi: float) -> str:
    if lo > 0:
        return "DOMAIN_RULE_AHEAD"
    if hi < 0:
        return "NEAREST_AHEAD"
    return "NO_DIFFERENCE_DETECTED"


def _fmt_p(p: float) -> str:
    return f"{p:.3g}"


def reading(res: dict[str, Any]) -> str:
    """The registered outcome wording, filled in. The templates are REGISTRATION['outcome_wording']."""
    p = res["primary"]
    ctl = res["control"]
    falsified = ctl["p_one_sided"] >= ALPHA
    lab = comparison_label(*p["delta_ci95_points"])
    w = REGISTRATION["outcome_wording"]
    fill = {
        "n": f"{p['n']:,}",
        "kd": f"{p['domain_concordant']:,}",
        "kn": f"{p['nearest_concordant']:,}",
        "cd": f"{p['domain_rate_pct']:.1f}",
        "cn": f"{p['nearest_rate_pct']:.1f}",
        "d": f"{p['delta_points']:+.2f}",
        "lo": f"{p['delta_ci95_points'][0]:+.2f}",
        "hi": f"{p['delta_ci95_points'][1]:+.2f}",
        "cc": f"{100 * ctl['expected_rate']:.1f}",
        "p": _fmt_p(ctl["p_one_sided"]),
        "k": f"{res['census']['discordant']:,}",
    }
    parts = [w["head"].format(**fill), w[lab].format(**fill)]
    parts.append((w["FALSIFIED"] if falsified else w["NOT_FALSIFIED"]).format(**fill))
    if lab == "DOMAIN_RULE_AHEAD":
        ig = res["strata"]["intergenic"]
        if ig["delta_ci95_points"] is None or ig["delta_ci95_points"][0] <= 0:
            parts.append(w["host_caveat"])
    parts.append(w["tail"])
    return " ".join(parts)


#: Written and committed before `score` read the RNA sheet. Not to be changed afterwards.
REGISTRATION: dict[str, Any] = {
    "claim": (
        "On C. elegans distal elements accessible in exactly one of five tissues, the gene the worm "
        "domain rule attributes is expressed highest in the element's tissue more often than the gene "
        "with the nearest annotated 5' end."
    ),
    "attribution": (
        "Domain rule (the worm instance of the human node model; C. elegans has no CTCF): host "
        "protein-coding gene if the element's midpoint lies in one (nearest 5' end among overlapping "
        "hosts); otherwise the flanking protein-coding gene whose 5' end faces the element (the element "
        "is upstream of it); divergent pair: nearer 5' end; convergent pair: nearer 5' end of the two, "
        "counted apart. Reads WBcel235 Ensembl gene models (positions, strands) only."
    ),
    "baseline": "Nearest annotated 5' end among all protein-coding genes of the chromosome (nearest TSS).",
    "population": (
        "Serizay et al. 2020 Table S2 ATAC_metrics rows with annot_dev_age_tissues == putative_enhancer "
        "and Annotation in {Germline, Neurons, Muscle, Hypod., Intest.} (accessible in exactly one "
        "tissue by the source's own call); midpoint = (start + end) // 2; chromosomes I-V and X."
    ),
    "endpoint": (
        "Concordant: the attributed gene's TPM (RNA_metrics, young adults) in the element's tissue is "
        ">= 5 and strictly higher than in each of the other four tissues. A gene absent from RNA_metrics "
        "is not concordant (kept in the denominator)."
    ),
    "never_read_for_attribution": list(NEVER_READ),
    "statistic": (
        "delta = concordance(domain rule) - concordance(nearest TSS), in points, on the same elements; "
        "95% interval: percentile cluster bootstrap over 1-Mb windows, 2,000 resamples, seed 20261010."
    ),
    "matched_control": (
        "Per element, a protein-coding gene of the same chromosome drawn uniformly among those whose "
        "5' end lies in the same distance bin as the domain-rule gene's (bins 0-1, 1-2, 2-5, 5-10, 10-20, "
        "20-50, 50-100, >100 kb), the domain-rule gene excluded; if none, the gene whose 5'-end distance "
        "is closest to the domain-rule gene's. 1,000 draws, seed 20261010; one-sided p = (1 + draws "
        "with concordance >= observed) / 1,001."
    ),
    "positive_control": (
        "Single-tissue coding_promoter elements, nearest-TSS gene, against its own matched-distance "
        "control built the same way. If its p >= 0.05 the readout is declared unresponsive and the "
        "primary reading is withheld."
    ),
    "falsifier": (
        "The domain rule is falsified as a carrier of tissue information on this endpoint if its "
        "concordance is not above its matched-distance control (one-sided p >= 0.05)."
    ),
    "decision": (
        "DOMAIN_RULE_AHEAD if the interval's lower bound > 0; NEAREST_AHEAD if its upper bound < 0; "
        "NO_DIFFERENCE_DETECTED otherwise."
    ),
    "strata_beside": (
        "Reported beside the primary, deciding nothing: host vs intergenic elements, each tissue, the "
        "discordant elements alone, and the nearest-TSS gene against its own matched control. If the "
        "primary is DOMAIN_RULE_AHEAD and the intergenic stratum's lower bound is <= 0, the host caveat "
        "is appended to the reading."
    ),
    "stop_rule": (
        "One run. If fewer than 200 elements have different genes under the two rules, stop unscored "
        "(UNDERPOWERED). No rule, population, threshold or endpoint is changed after the RNA sheet is "
        "read; a second definition would be a new registration, committed before it is run."
    ),
    "outcome_wording": {
        "head": (
            "On {n} C. elegans putative enhancers accessible in exactly one of five tissues (Serizay et "
            "al. 2020, young adults), the gene the domain rule attributes is expressed at 5 TPM or more "
            "and highest in the element's own tissue for {kd} of {n} ({cd}%), against {kn} of {n} ({cn}%) "
            "for the gene with the nearest 5' end; the rules differ on {k} elements; difference {d} "
            "points, 95% {lo} to {hi} over 1-Mb clusters."
        ),
        "DOMAIN_RULE_AHEAD": "On this endpoint the domain rule is ahead of nearest TSS.",
        "NEAREST_AHEAD": "On this endpoint nearest TSS is ahead of the domain rule.",
        "NO_DIFFERENCE_DETECTED": "No difference detected between the domain rule and nearest TSS.",
        "NOT_FALSIFIED": (
            "The domain rule's genes match the element's tissue more often than genes at a matched "
            "distance ({cc}%, one-sided p {p}), so position carries tissue information here."
        ),
        "FALSIFIED": (
            "The domain rule's genes match the element's tissue no more often than genes at a matched "
            "distance ({cc}%, one-sided p {p}): on this endpoint the attribution carries no tissue "
            "information beyond distance, and it is falsified as a carrier of it."
        ),
        "host_caveat": (
            "The difference is not carried by the intergenic elements; inside gene bodies accessibility "
            "can follow the host gene's own transcription in the same nuclei, so this run does not "
            "separate regulation from transcription."
        ),
        "tail": (
            "Accessibility is not function: this says the attributed gene is expressed where the "
            "element is open, not that the element acts on it."
        ),
    },
}


def score(
    chroms: dict[str, Chromosome], rows: list[dict[str, Any]], expr: dict[str, tuple[float, ...]]
) -> dict[str, Any]:
    """The registered run. `rows` come from `attribute`; `expr` from `read_expression`."""
    enh = [r for r in rows if r["kind"] == ELEMENT_TYPE]
    prom = [r for r in rows if r["kind"] == PROMOTER_TYPE]
    cen = census(rows)
    if cen["stop_rule_met"]:
        return {"outcome": "UNDERPOWERED", "census": cen}
    idx = ControlIndex(chroms, expr)
    primary = compare(enh, expr)
    probs = [
        idx.probability(r["chrom"], r["tissue"], r["mid"], r["domain_gene"], r["domain_dist"]) for r in enh
    ]
    control = control_test(probs, primary["domain_concordant"])
    # positive control: promoters, nearest-TSS gene against its own matched control
    for r in prom:
        r["c_nearest"] = int(concordant(expr.get(r["nearest_gene"]), r["tissue"]))
    kp = sum(r["c_nearest"] for r in prom)
    pprobs = [
        idx.probability(r["chrom"], r["tissue"], r["mid"], r["nearest_gene"], r["nearest_dist"]) for r in prom
    ]
    pos = {"n": len(prom), "nearest_concordant": kp, "nearest_rate_pct": _rate(kp, len(prom))}
    pos["control"] = control_test(pprobs, kp)
    pos["responsive"] = pos["control"]["p_one_sided"] < ALPHA
    # beside, deciding nothing
    nprobs = [
        idx.probability(r["chrom"], r["tissue"], r["mid"], r["nearest_gene"], r["nearest_dist"]) for r in enh
    ]
    strata = {
        "host": compare([r for r in enh if r["domain_case"] == "host"], expr),
        "intergenic": compare([r for r in enh if r["domain_case"] != "host"], expr),
        "discordant_only": compare([r for r in enh if r["domain_gene"] != r["nearest_gene"]], expr),
    }
    for t in TISSUES:
        strata["tissue_" + t] = compare([r for r in enh if r["tissue"] == t], expr)
    missing = sum(1 for r in enh for g in (r["domain_gene"], r["nearest_gene"]) if g not in expr)
    res = {
        "census": cen,
        "primary": primary,
        "control": control,
        "positive_control": pos,
        "nearest_against_its_control": control_test(nprobs, primary["nearest_concordant"]),
        "strata_beside": strata,
        "strata": strata,
        "attributed_genes_missing_from_rna_sheet": missing,
        "comparison": comparison_label(*primary["delta_ci95_points"]),
        "falsifier": "FALSIFIED" if control["p_one_sided"] >= ALPHA else "NOT_FALSIFIED",
    }
    if not pos["responsive"]:
        res["outcome"] = "READOUT_UNRESPONSIVE"
        res["reading"] = (
            "The positive control did not separate from its matched control, so the readout is "
            "declared unresponsive and the primary reading is withheld."
        )
    else:
        res["outcome"] = res["comparison"]
        res["reading"] = reading(res)
    del res["strata"]
    return res


def fetch_table(dest: Path = TABLE_S2) -> dict[str, Any]:
    """Fetch the article's supplementary bundle, keep Table S2 only, discard the rest."""
    import urllib.request

    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.parent / "serizay2020_supp.zip"
    with urllib.request.urlopen(TABLE_S2_URL, timeout=300) as r, open(tmp, "wb") as f:
        while chunk := r.read(1 << 20):
            f.write(chunk)
    bundle = tmp.stat().st_size
    with zipfile.ZipFile(tmp) as z, open(dest, "wb") as f:
        f.write(z.read(TABLE_S2_MEMBER))
    tmp.unlink()
    return {"bundle_bytes": bundle, "kept_bytes": dest.stat().st_size, "sha256": sha256(dest)}
