# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Simpler genomes, the fly: the worm's domain-rule test repeated on Drosophila embryos.

Area I, simpler genomes, after the worm (genomeos/attribution/simpler.py, docs/SIMPLER-GENOMES.md
sections 3, 4 and 8). The worm test found that a gene-structure node (host gene, else the gene the
element sits upstream of) names a gene expressed highest in the element's tissue more often than
nearest TSS, on one source (Serizay et al. 2020, young adults). This module asks the same question,
with the same rule, baseline, control, falsifier and stop rule, on a second species and a second
source: Calderon et al. 2022 (Science 377:eabn5800), single-cell accessibility (sci-ATAC) and
single-cell RNA (sci-RNA) over 0-20 h Drosophila melanogaster embryos, each modality annotated with
tissue labels by the source.

THE NODE. Drosophila has insulator proteins (dCTCF, BEAF-32, CP190, Su(Hw)) and TADs, so unlike
the worm a boundary-delimited node could be built. It is not, and the decision is written in
docs/SIMPLER-GENOMES.md section 8 before anything was scored: this is a replication, so the node is
the worm's rule unchanged (`simpler.domain_rule`, imported, not copied). A boundary node would be a
separate registration.

THE POPULATION. The source's own call: Table S4 lists, per time window and ATAC cluster, the peaks
more accessible in that cluster than in the window's other cells (Seurat markers, adjusted p < 0.05).
Table S2 names each cluster. An element is accessible in exactly one tissue when every positive
marker call it has, over all windows, falls in clusters of one tissue of `TISSUES` and none falls in a
cluster with a label outside them. The labels of both modalities are grouped into the same tissues by
name in `ATAC_TISSUE` and `RNA_TISSUE`; nothing about expression informs the grouping.

THE ENDPOINT. Per tissue, the UMI counts of every RNA cell the source annotates with one of that
tissue's labels, summed over all windows (as the source's own pseudobulk script sums them) and
scaled to counts per million. Concordant: the attributed gene's CPM in the element's tissue is at
least 5 and higher than in each other tissue. The attribution reads gene positions and strands only.

WHAT THIS CANNOT ESTABLISH. Accessibility is not function. The two assays come from different nuclei
of comparably staged embryos, and the RNA is pooled over 0-20 h while an element's marker call may
come from one window. An element inside a gene body can be open where its host is transcribed for
that reason alone.
"""

from __future__ import annotations

import bisect
import gzip
import math
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from genomeos.attribution import simpler
from genomeos.attribution.simpler import Chromosome, Element, Gene

STORE = Path("data/knowledge/simpler_fly")
SUPP = STORE / "SupplementaryTables.xlsx"
PEAKS = STORE / "peaks.bed"
GFF3 = STORE / "Drosophila_melanogaster.BDGP6.32.104.gff3.gz"
RNA = STORE / "rna_umi_by_manual_annot.tsv.gz"
FETCH_RECORD = STORE / "fetch_record.json"

SOURCE = (
    "Calderon D, Blecher-Gonen R, Huang X, Secchia S, Kentro J, Daza RM, Martin B, Dulja A, Schaub C, "
    "Trapnell C, Larschan E, O'Connor-Giles KM, Furlong EEM, Shendure J (2022) The continuum of "
    "Drosophila embryonic development at single-cell resolution. Science 377:eabn5800, "
    "doi:10.1126/science.abn5800; supplementary tables, peak set and cell annotations from the "
    "authors' site (DEAP), RNA counts from GEO GSE190147"
)
GENE_MODELS = "Ensembl release 104, Drosophila melanogaster BDGP6.32 (dm6) GFF3"

#: The tissues both modalities name, in a fixed order. Each is a group of the source's labels.
TISSUES = (
    "neural",
    "ventral midline",
    "epidermis",
    "foregut and hindgut",
    "tracheal system",
    "salivary gland",
    "Malpighian tubule",
    "muscle",
    "fat body",
    "plasmatocytes",
    "midgut",
    "amnioserosa",
    "germ cell",
    "yolk",
)

#: ATAC: Table S2 `refined_broad_annotation` -> tissue. Every label in the table is here or excluded.
ATAC_TISSUE = {
    "Neural": "neural",
    "Neural prim.": "neural",
    "Brain": "neural",
    "Brain prim.": "neural",
    "Ventral midline": "ventral midline",
    "Epidermis": "epidermis",
    "Epidermis prim.": "epidermis",
    "Hindgut": "foregut and hindgut",
    "Hindgut prim.": "foregut and hindgut",
    "Foregut prim.": "foregut and hindgut",
    "Pharnyx": "foregut and hindgut",  # the source's spelling
    "Tracheal system": "tracheal system",
    "Tracheal system prim.": "tracheal system",
    "Salivary gland": "salivary gland",
    "Malpighian tubule": "Malpighian tubule",
    "Muscle": "muscle",
    "Muscle prim.": "muscle",
    "Fat body": "fat body",
    "Plasmatocytes": "plasmatocytes",
    "Midgut": "midgut",  # includes the source's one Proventriculus cluster
    "Midgut prim.": "midgut",
    "Amnioserosa": "amnioserosa",
    "Amnioserosa anlage": "amnioserosa",
    "Germ cell": "germ cell",
    "Yolk": "yolk",
}
#: ATAC labels that name no tissue both modalities resolve: a germ layer's or a region's anlage,
#: the blastoderm, head ectoderm (no RNA label), and unknown.
ATAC_EXCLUDED = frozenset(
    {
        "Blastoderm",
        "Ectoderm anlage",
        "Endoderm anlage",
        "Mesoderm anlage",
        "Mesectoderm anlage",
        "Head ectoderm",
        "Head ectoderm prim.",
        "Unknown",
    }
)
#: RNA: rna_meta.rds `manual_annot` -> tissue. Every label in the table is here or excluded.
RNA_TISSUE = {
    "brain": "neural",
    "brain primordium": "neural",
    "ventral nerve cord": "neural",
    "ventral nerve cord primordium": "neural",
    "peripheral nervous system": "neural",
    "peripheral nervous system primordium": "neural",
    "sensory nervous system": "neural",
    "sensory nervous system primordium": "neural",
    "glia": "neural",
    "midline primordium": "ventral midline",
    "epidermis": "epidermis",
    "epidermis primordium": "epidermis",
    "foregut / hindgut / anal pad": "foregut and hindgut",
    "foregut / hindgut primordium": "foregut and hindgut",
    "hindgut / anal pad": "foregut and hindgut",
    "hindgut primordium": "foregut and hindgut",
    "tracheal system": "tracheal system",
    "tracheal primordium": "tracheal system",
    "spiracle (tracheal system)": "tracheal system",
    "spiracles (tracheal system)": "tracheal system",
    "salivary gland": "salivary gland",
    "salivary gland primordium": "salivary gland",
    "salivary gland body primordium": "salivary gland",
    "Malpighian tubule": "Malpighian tubule",
    "somatic muscle": "muscle",
    "somatic muscle primordium": "muscle",
    "muscle": "muscle",
    "visceral muscle": "muscle",
    "visceral muscle primordium": "muscle",
    "fat body": "fat body",
    "plasmatocytes": "plasmatocytes",
    "plasmatocytes anlage": "plasmatocytes",
    "midgut": "midgut",
    "midgut primordium": "midgut",
    "proventriculus": "midgut",  # as the source's ATAC annotation places it
    "amnioserosa": "amnioserosa",
    "germ cell": "germ cell",
    "yolk nuclei": "yolk",
}
#: RNA labels that name no tissue both modalities resolve (crystal cells have no ATAC cluster).
RNA_EXCLUDED = frozenset(
    {
        "maternal",
        "ubiquitous",
        "unknown",
        "cellular blastoderm / pole cell",
        "ectoderm primordium",
        "ectoderm anlage",
        "ectoderm anlage in statu nascendi",
        "anlage in statu nascendi",
        "mesoderm primordium",
        "trunk mesoderm primordium",
        "endoderm primordium",
        "crystal cell",
        "procrystal cell",
    }
)

#: The element classes. The primary population is distal; the positive control is promoters.
DISTAL = "distal"
PROMOTER = "promoter"
#: A promoter element's midpoint lies within this many bp of a protein-coding gene's 5' end.
PROMOTER_BP = 500
#: A distal element's midpoint lies farther than this from every annotated transcript 5' end.
DISTAL_BP = 1_000
#: Transcript feature types whose 5' ends count against "distal" (every biotype in the GFF3).
TRANSCRIPT_TYPES = frozenset(
    {"mRNA", "ncRNA", "pseudogenic_transcript", "tRNA", "snoRNA", "pre_miRNA", "rRNA", "snRNA"}
)
CHROMOSOMES = ("2L", "2R", "3L", "3R", "4", "X", "Y")
#: Counts per million in the tissue pseudobulk; set at the worm's 5 TPM before any value was read.
DETECTION_CPM = 5.0
MARKER_P_ADJ = 0.05

#: Table S4 and Table S2 columns (0-based).
S4_COLS = {"log_fc": 1, "p_adj": 4, "cluster": 5, "peak": 6, "sample": 7}
S2_COLS = {"time": 0, "cluster": 1, "broad": 9}
#: Sheets of the supplementary workbook that carry RNA values or RNA-derived links. Never read.
NEVER_READ = ("TableS3", "TableS5", "TableS6", "TableS8", "TableS9", "TableS10")

# The worm's statistics, unchanged.
DISTANCE_BINS = simpler.DISTANCE_BINS
CONTROL_DRAWS = simpler.CONTROL_DRAWS
BOOTSTRAP = simpler.BOOTSTRAP
CLUSTER_BP = simpler.CLUSTER_BP
SEED = simpler.SEED
MIN_DISCORDANT = simpler.MIN_DISCORDANT
ALPHA = simpler.ALPHA

#: The worm's registered result this run replicates or not (docs/SIMPLER-GENOMES.md section 4).
WORM = {"delta_points": 2.89, "ci95_points": [1.83, 3.91], "registration": "595d7c2"}


# --------------------------------------------------------------------------------------------- reading


def _attrs(field: str) -> dict[str, str]:
    return dict(kv.split("=", 1) for kv in field.strip().split(";") if "=" in kv)


def read_annotation(path: Path = GFF3) -> tuple[dict[str, Chromosome], dict[str, list[int]], dict[str, str]]:
    """Protein-coding genes per chromosome, every transcript 5' end per chromosome, and name -> FBgn.

    Genes are `simpler.Gene` with the FlyBase gene ID as `gid` and the symbol as `name`; chromosomes
    are named the way the peak set names them (`chr2L`, ...).
    """
    per: dict[str, list[Gene]] = {}
    tss: dict[str, list[int]] = {}
    name_to_id: dict[str, str] = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[0] not in CHROMOSOMES:
                continue
            chrom = "chr" + p[0]
            if p[2] in TRANSCRIPT_TYPES:
                tss.setdefault(chrom, []).append(int(p[3]) if p[6] == "+" else int(p[4]))
            elif p[2] == "gene":
                a = _attrs(p[8])
                if a.get("biotype") != "protein_coding":
                    continue
                gid = a["gene_id"]
                name = a.get("Name", gid)
                per.setdefault(chrom, []).append(Gene(gid, name, chrom, int(p[3]), int(p[4]), p[6]))
                name_to_id[name] = gid
    return {c: Chromosome(g) for c, g in per.items()}, {c: sorted(v) for c, v in tss.items()}, name_to_id


def read_peaks(path: Path = PEAKS) -> set[str]:
    with open(path) as f:
        return {"-".join(line.split()[:3]) for line in f if line.strip()}


def _sheet(path: Path, name: str) -> Iterable[tuple]:
    if name in NEVER_READ:
        raise ValueError(f"{name} carries RNA values or RNA-derived links and is never read")
    import openpyxl

    wb = openpyxl.load_workbook(path, read_only=True)
    rows = wb[name].iter_rows(values_only=True)
    next(rows)  # header
    for r in rows:
        if r and r[0] is not None:
            yield r


def read_cluster_labels(path: Path = SUPP) -> dict[tuple[str, int], str]:
    """(time window, ATAC cluster) -> the source's refined broad annotation (Table S2)."""
    c = S2_COLS
    return {(str(r[c["time"]]), int(r[c["cluster"]])): str(r[c["broad"]]) for r in _sheet(path, "TableS2")}


def read_markers(path: Path = SUPP) -> list[tuple[str, str, int, float, float]]:
    """Table S4 rows: (peak, window, cluster, avg_logFC, p_val_adj)."""
    c = S4_COLS
    out = []
    for r in _sheet(path, "TableS4"):
        window = str(r[c["sample"]]).replace("fly.new.", "")
        out.append(
            (str(r[c["peak"]]), window, int(r[c["cluster"]]), float(r[c["log_fc"]]), float(r[c["p_adj"]]))
        )
    return out


def tissue_calls(
    markers: list[tuple[str, str, int, float, float]], labels: dict[tuple[str, int], str]
) -> tuple[dict[str, str], dict[str, Any]]:
    """Peak -> its one tissue, for peaks whose positive marker calls all fall in one tissue's clusters."""
    seen = set(labels.values())
    unknown = seen - set(ATAC_TISSUE) - ATAC_EXCLUDED
    if unknown:
        raise ValueError(f"ATAC labels neither grouped nor excluded: {sorted(unknown)}")
    groups: dict[str, set[str]] = {}
    pos = neg = 0
    for peak, window, cluster, lfc, padj in markers:
        if padj >= MARKER_P_ADJ or lfc <= 0:
            neg += 1
            continue
        pos += 1
        lab = labels[(window, cluster)]
        groups.setdefault(peak, set()).add(ATAC_TISSUE.get(lab, "excluded:" + lab))
    calls = {}
    multi = excluded = 0
    for peak, ts in groups.items():
        if any(t.startswith("excluded:") for t in ts):
            excluded += 1
        elif len(ts) > 1:
            multi += 1
        else:
            calls[peak] = next(iter(ts))
    funnel = {
        "marker_rows": len(markers),
        "marker_rows_not_positive_or_not_significant": neg,
        "positive_marker_rows": pos,
        "peaks_with_a_positive_marker_call": len(groups),
        "peaks_with_a_call_in_an_excluded_label": excluded,
        "peaks_called_in_two_or_more_tissues": multi,
        "peaks_called_in_exactly_one_tissue": len(calls),
    }
    return calls, funnel


def classify(chrom: str, mid: int, chroms: dict[str, Chromosome], tss: dict[str, list[int]]) -> str | None:
    """Promoter, distal, or neither (None). Reads annotated positions only."""
    ch = chroms.get(chrom)
    if ch is None:
        return None
    n = ch.nearest_tss(mid)
    if n is not None and abs(n.tss - mid) <= PROMOTER_BP:
        return PROMOTER
    t = tss.get(chrom, [])
    i = bisect.bisect_left(t, mid)
    near = min((abs(t[j] - mid) for j in (i - 1, i) if 0 <= j < len(t)), default=math.inf)
    if near > DISTAL_BP:
        return DISTAL
    return None


def elements(
    calls: dict[str, str], peaks: set[str], chroms: dict[str, Chromosome], tss: dict[str, list[int]]
) -> tuple[list[Element], dict[str, int]]:
    out = []
    counts = {"not_in_peak_set": 0, "chromosome_not_used": 0, "neither_class": 0, DISTAL: 0, PROMOTER: 0}
    for peak, tissue in sorted(calls.items()):
        if peak not in peaks:
            counts["not_in_peak_set"] += 1
            continue
        chrom, s, e = peak.rsplit("-", 2)
        if chrom.removeprefix("chr") not in CHROMOSOMES:
            counts["chromosome_not_used"] += 1
            continue
        el = Element(chrom, int(s), int(e), "", tissue)
        kind = classify(chrom, el.mid, chroms, tss)
        if kind is None:
            counts["neither_class"] += 1
            continue
        counts[kind] += 1
        out.append(Element(chrom, int(s), int(e), kind, tissue))
    return out, counts


def read_expression(
    path: Path, name_to_id: dict[str, str]
) -> tuple[dict[str, tuple[float, ...]], dict[str, Any]]:
    """The endpoint: per protein-coding gene (FBgn), CPM in each of `TISSUES`. Read only by `score`."""
    with gzip.open(path, "rt") as f:
        head = f.readline().rstrip("\n").split("\t")[1:]
        labels = [h.strip().replace("\xa0", "") for h in head]
        unknown = set(labels) - set(RNA_TISSUE) - RNA_EXCLUDED
        if unknown:
            raise ValueError(f"RNA labels neither grouped nor excluded: {sorted(unknown)}")
        col_t = [TISSUES.index(RNA_TISSUE[lab]) if lab in RNA_TISSUE else -1 for lab in labels]
        sums: dict[str, list[float]] = {}
        totals = [0.0] * len(TISSUES)
        for line in f:
            p = line.rstrip("\n").split("\t")
            v = [0.0] * len(TISSUES)
            for j, x in zip(col_t, p[1:], strict=True):
                if j >= 0:
                    v[j] += float(x)
            for j in range(len(TISSUES)):
                totals[j] += v[j]
            if p[0] in name_to_id:
                sums[name_to_id[p[0]]] = v
    expr = {
        g: tuple(1e6 * v[j] / totals[j] if totals[j] else 0.0 for j in range(len(TISSUES)))
        for g, v in sums.items()
    }
    return expr, {
        "protein_coding_genes_in_rna": len(expr),
        "tissue_umi_totals": dict(zip(TISSUES, totals, strict=True)),
    }


def concordant(cpm: tuple[float, ...] | None, tissue: str) -> bool:
    """At least DETECTION_CPM in the tissue and higher there than in each other tissue."""
    if cpm is None:
        return False
    i = TISSUES.index(tissue)
    v = cpm[i]
    return v >= DETECTION_CPM and all(v > cpm[j] for j in range(len(TISSUES)) if j != i)


# --------------------------------------------------------------------------------------------- census


def census(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Population counts and power in advance, from positions and the element-side tissue call only."""
    enh = [r for r in rows if r["kind"] == DISTAL]
    prom = [r for r in rows if r["kind"] == PROMOTER]
    disc = [r for r in enh if r["domain_gene"] != r["nearest_gene"]]
    cases: dict[str, int] = {}
    disc_cases: dict[str, int] = {}
    for r in enh:
        cases[r["domain_case"]] = cases.get(r["domain_case"], 0) + 1
    for r in disc:
        disc_cases[r["domain_case"]] = disc_cases.get(r["domain_case"], 0) + 1
    n, k = len(enh), len(disc)
    mde = 2.8 * math.sqrt(k) / n if n else None
    return {
        "elements": n,
        "elements_by_tissue": {t: sum(1 for r in enh if r["tissue"] == t) for t in TISSUES},
        "domain_cases": dict(sorted(cases.items())),
        "discordant": k,
        "discordant_fraction": round(k / n, 4) if n else None,
        "discordant_by_case": dict(sorted(disc_cases.items())),
        "discordant_by_tissue": {t: sum(1 for r in disc if r["tissue"] == t) for t in TISSUES},
        "host_elements": cases.get("host", 0),
        "intergenic_elements": n - cases.get("host", 0),
        "intergenic_discordant": k - disc_cases.get("host", 0),
        "clusters_1mb": len({(r["chrom"], r["mid"] // CLUSTER_BP) for r in enh}),
        "mde_points_unclustered": round(100 * mde, 3) if mde is not None else None,
        "positive_control_promoters": len(prom),
        "positive_control_promoters_by_tissue": {
            t: sum(1 for r in prom if r["tissue"] == t) for t in TISSUES
        },
        "stop_rule_met": k < MIN_DISCORDANT,
    }


# --------------------------------------------------------------------------------------------- scoring


class ControlIndex:
    """`simpler.ControlIndex` over the fly's tissues: the same matched-distance draw, exactly.

    A matched-control gene is drawn uniformly from the protein-coding genes of the element's
    chromosome whose 5' end lies in the same distance bin as the attributed gene's, the attributed
    gene excluded; the chance it is concordant is computed exactly from prefix counts.
    """

    def __init__(self, chroms: dict[str, Chromosome], expr: dict[str, tuple[float, ...]]):
        self.chroms = chroms
        self.expr = expr
        self.prefix: dict[tuple[str, str], list[int]] = {}
        for c, ch in chroms.items():
            for t in TISSUES:
                acc, p = 0, [0]
                for g in ch.by_tss:
                    acc += concordant(expr.get(g.gid), t)
                    p.append(acc)
                self.prefix[(c, t)] = p

    def _count(self, chrom: str, tissue: str, lo_pos: float, hi_pos: float) -> tuple[int, int]:
        tss = self.chroms[chrom].tss
        i, j = bisect.bisect_left(tss, lo_pos), bisect.bisect_left(tss, hi_pos)
        p = self.prefix[(chrom, tissue)]
        return j - i, p[j] - p[i]

    def probability(self, chrom: str, tissue: str, mid: int, gene: str, dist: int) -> float:
        b = simpler.distance_bin(dist)
        lo = 0 if b == 0 else DISTANCE_BINS[b - 1]
        hi = DISTANCE_BINS[b] if b < len(DISTANCE_BINS) else math.inf
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
        ch = self.chroms[chrom]
        others = [g for g in ch.genes if g.gid != gene]
        if not others:
            return 0.0
        sub = min(others, key=lambda g: (abs(abs(g.tss - mid) - dist), g.gid))
        return float(concordant(self.expr.get(sub.gid), tissue))


def compare(rows: list[dict[str, Any]], expr: dict[str, tuple[float, ...]]) -> dict[str, Any]:
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
        "domain_concordant": kd,
        "nearest_concordant": kn,
        "domain_rate_pct": simpler._rate(kd, n),
        "nearest_rate_pct": simpler._rate(kn, n),
        "delta_points": round(100 * bs["point"], 3) if bs else None,
        "delta_ci95_points": [round(100 * bs["lo"], 3), round(100 * bs["hi"], 3)] if bs else None,
        "clusters": bs["clusters"] if bs else 0,
    }


def reading(res: dict[str, Any]) -> str:
    """The registered outcome wording, filled in. The templates are REGISTRATION['outcome_wording']."""
    p = res["primary"]
    ctl = res["control"]
    falsified = ctl["p_one_sided"] >= ALPHA
    lab = simpler.comparison_label(*p["delta_ci95_points"])
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
        "p": simpler._fmt_p(ctl["p_one_sided"]),
        "k": f"{res['census']['discordant']:,}",
        "t": f"{len(TISSUES)}",
    }
    parts = [w["head"].format(**fill), w[lab].format(**fill), w["replication_" + lab].format(**fill)]
    parts.append((w["FALSIFIED"] if falsified else w["NOT_FALSIFIED"]).format(**fill))
    if lab == "DOMAIN_RULE_AHEAD":
        ig = res["strata"]["intergenic"]
        if ig["delta_ci95_points"] is None or ig["delta_ci95_points"][0] <= 0:
            parts.append(w["host_caveat"])
    parts.append(w["tail"])
    return " ".join(parts)


#: Written and committed before `score` read any RNA value. Not to be changed afterwards.
REGISTRATION: dict[str, Any] = {
    "claim": (
        "On Drosophila melanogaster embryo distal elements accessible in clusters of exactly one tissue, "
        "the gene the worm's domain rule attributes is expressed highest in the element's tissue more "
        "often than the gene with the nearest annotated 5' end: the worm result (Serizay 2020, +2.89 "
        "points, registration 595d7c2) repeated on a second species and a second source."
    ),
    "node": (
        "The worm's gene-structure rule, unchanged (simpler.domain_rule, imported). Drosophila has "
        "insulators and TADs, but no measured boundary set is used: this is a replication, and changing "
        "the node together with the species would leave a disagreement with the worm unattributable. "
        "Decision written in docs/SIMPLER-GENOMES.md section 8 before any RNA value was read."
    ),
    "attribution": (
        "Domain rule: host protein-coding gene if the element's midpoint lies in one (nearest 5' end "
        "among overlapping hosts); otherwise the flanking protein-coding gene whose 5' end faces the "
        "element; divergent pair: nearer 5' end; convergent pair: nearer 5' end of the two, counted "
        "apart. Reads Ensembl 104 BDGP6.32 gene positions and strands only."
    ),
    "baseline": (
        "Nearest annotated 5' end among all protein-coding genes of the chromosome arm (nearest TSS)."
    ),
    "tissues": list(TISSUES),
    "tissue_grouping": {"atac_refined_broad_annotation": ATAC_TISSUE, "rna_manual_annot": RNA_TISSUE},
    "labels_excluded": {"atac": sorted(ATAC_EXCLUDED), "rna": sorted(RNA_EXCLUDED)},
    "population": (
        "Calderon 2022 Table S4 marker rows with p_val_adj < 0.05 and avg_logFC > 0, joined to Table S2 "
        "on (time window, cluster); a peak is single-tissue when every such call falls in clusters of "
        "one tissue and none in an excluded label. Distal: midpoint ((start + end) // 2 of the peak) "
        "more than 1,000 bp from every annotated transcript 5' end of any biotype; chromosome arms 2L, "
        "2R, 3L, 3R, 4, X, Y."
    ),
    "endpoint": (
        "Per tissue, UMI counts summed over every RNA cell whose manual_annot is grouped into it, all "
        "time windows pooled, scaled to counts per million of that tissue's total. Concordant: the "
        "attributed gene's CPM in the element's tissue is >= 5 and strictly higher than in each of the "
        "other 13 tissues. A protein-coding gene absent from the RNA rows is not concordant (kept in "
        "the denominator)."
    ),
    "never_read": (
        "for the attribution: any accessibility or RNA value; at all: supplementary sheets "
        + ", ".join(NEVER_READ)
        + " (RNA markers, RNA-derived gene clusters and motif-expression links, NNLS RNA-ATAC links)"
    ),
    "statistic": (
        "delta = concordance(domain rule) - concordance(nearest TSS), in points, on the same elements; "
        "95% interval: percentile cluster bootstrap over 1-Mb windows, 2,000 resamples, seed 20261010 "
        "(simpler.cluster_bootstrap)."
    ),
    "matched_control": (
        "Per element, a protein-coding gene of the same chromosome arm drawn uniformly among those whose "
        "5' end lies in the same distance bin as the domain-rule gene's (bins 0-1, 1-2, 2-5, 5-10, 10-20, "
        "20-50, 50-100, >100 kb), the domain-rule gene excluded; if none, the gene whose 5'-end distance "
        "is closest. 1,000 draws, seed 20261010; one-sided p = (1 + draws with concordance >= observed) "
        "/ 1,001 (simpler.control_test)."
    ),
    "positive_control": (
        "Single-tissue promoter elements (midpoint within 500 bp of a protein-coding gene's annotated 5' "
        "end), nearest-TSS gene, against its own matched-distance control built the same way. If its "
        "p >= 0.05 the readout is declared unresponsive and the primary reading is withheld."
    ),
    "falsifier": (
        "The domain rule is falsified as a carrier of tissue information on this endpoint if its "
        "concordance is not above its matched-distance control (one-sided p >= 0.05)."
    ),
    "decision": (
        "DOMAIN_RULE_AHEAD if the interval's lower bound > 0; NEAREST_AHEAD if its upper bound < 0; "
        "NO_DIFFERENCE_DETECTED otherwise. Replication of the worm is read from the same label: "
        "DOMAIN_RULE_AHEAD replicates its direction, NO_DIFFERENCE_DETECTED does not replicate it, "
        "NEAREST_AHEAD reverses it."
    ),
    "strata_beside": (
        "Reported beside the primary, deciding nothing: host vs intergenic elements, each tissue, the "
        "discordant elements alone, and the nearest-TSS gene against its own matched control. If the "
        "primary is DOMAIN_RULE_AHEAD and the intergenic stratum's lower bound is <= 0, the host caveat "
        "is appended to the reading."
    ),
    "stop_rule": (
        "One run. If fewer than 200 elements have different genes under the two rules, stop unscored "
        "(UNDERPOWERED). No rule, grouping, population, threshold or endpoint is changed after an RNA "
        "value is read; a second definition would be a new registration, committed before it is run."
    ),
    "cannot_establish": (
        "Accessibility is not function. Concordance says the attributed gene is expressed where the "
        "element is open, not that the element acts on it. The assays come from different nuclei; the "
        "RNA is pooled over 0-20 h while a marker call may come from one window; tissue labels were "
        "assigned by the source separately in each modality and grouped here by name. Inside gene "
        "bodies accessibility can follow the host's own transcription."
    ),
    "outcome_wording": {
        "head": (
            "On {n} D. melanogaster embryo distal elements accessible in clusters of exactly one of {t} "
            "tissues (Calderon et al. 2022, 0-20 h, the source's own marker call), the gene the domain "
            "rule attributes is expressed at 5 counts per million or more and highest in the element's "
            "own tissue, in RNA pooled over that tissue's cells, for {kd} of {n} ({cd}%), against {kn} of "
            "{n} ({cn}%) for the gene with the nearest 5' end; the rules differ on {k} elements; "
            "difference {d} points, 95% {lo} to {hi} over 1-Mb clusters."
        ),
        "DOMAIN_RULE_AHEAD": "On this endpoint the domain rule is ahead of nearest TSS.",
        "NEAREST_AHEAD": "On this endpoint nearest TSS is ahead of the domain rule.",
        "NO_DIFFERENCE_DETECTED": "No difference detected between the domain rule and nearest TSS.",
        "replication_DOMAIN_RULE_AHEAD": (
            "The worm's result (+2.89 points on Serizay et al. 2020) is replicated in direction on a "
            "second species and source."
        ),
        "replication_NO_DIFFERENCE_DETECTED": (
            "The worm's result (+2.89 points on Serizay et al. 2020) is not replicated on this source."
        ),
        "replication_NEAREST_AHEAD": (
            "The worm's result (+2.89 points on Serizay et al. 2020) is reversed on this source."
        ),
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
            "can follow the host gene's own transcription, so this run does not separate regulation from "
            "transcription."
        ),
        "tail": (
            "Accessibility is not function: this says the attributed gene is expressed where the "
            "element is open, not that the element acts on it, and the two assays come from different "
            "nuclei."
        ),
    },
}


def score(
    chroms: dict[str, Chromosome], rows: list[dict[str, Any]], expr: dict[str, tuple[float, ...]]
) -> dict[str, Any]:
    """The registered run. `rows` come from `simpler.attribute`; `expr` from `read_expression`."""
    enh = [r for r in rows if r["kind"] == DISTAL]
    prom = [r for r in rows if r["kind"] == PROMOTER]
    cen = census(rows)
    if cen["stop_rule_met"]:
        return {"outcome": "UNDERPOWERED", "census": cen}
    idx = ControlIndex(chroms, expr)
    primary = compare(enh, expr)
    probs = [
        idx.probability(r["chrom"], r["tissue"], r["mid"], r["domain_gene"], r["domain_dist"]) for r in enh
    ]
    control = simpler.control_test(probs, primary["domain_concordant"])
    for r in prom:
        r["c_nearest"] = int(concordant(expr.get(r["nearest_gene"]), r["tissue"]))
    kp = sum(r["c_nearest"] for r in prom)
    pprobs = [
        idx.probability(r["chrom"], r["tissue"], r["mid"], r["nearest_gene"], r["nearest_dist"]) for r in prom
    ]
    pos: dict[str, Any] = {
        "n": len(prom),
        "nearest_concordant": kp,
        "nearest_rate_pct": simpler._rate(kp, len(prom)),
    }
    pos["control"] = simpler.control_test(pprobs, kp)
    pos["responsive"] = pos["control"]["p_one_sided"] < ALPHA
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
    res: dict[str, Any] = {
        "census": cen,
        "primary": primary,
        "control": control,
        "positive_control": pos,
        "nearest_against_its_control": simpler.control_test(nprobs, primary["nearest_concordant"]),
        "strata_beside": strata,
        "strata": strata,
        "attributed_genes_missing_from_rna_rows": missing,
        "comparison": simpler.comparison_label(*primary["delta_ci95_points"]),
        "falsifier": "FALSIFIED" if control["p_one_sided"] >= ALPHA else "NOT_FALSIFIED",
        "worm": WORM,
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
