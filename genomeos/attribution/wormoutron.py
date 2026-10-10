# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Simpler genomes, area I: does the worm's between-gene lead survive excluding outrons and operons.

Stage 1 (`genomeos/attribution/simpler.py`, docs/SIMPLER-GENOMES.md sections 3 and 4) measured the
worm domain rule against nearest TSS on 5,814 single-tissue putative enhancers and read
DOMAIN_RULE_AHEAD, +2.89 points [+1.83, +3.91]. Its between-gene stratum, where no host gene can be
named, was +3.73 [+2.79, +4.74] on 2,335 elements (docs/SIMPLER-GENOMES.md line 124).

That stratum is what tonight's fly finding makes worth re-examining. lane-flysource measured in
Drosophila that the domain rule's lead over nearest TSS appears only where a host gene can be named
and is absent between genes (+0.03 [-0.50, +0.60] on 789 elements, docs/SIMPLER-GENOMES.md section
11), consistent with the fly's lead coming from host-gene transcription rather than regulation. The
worm is NOT qualified by that: its between-gene lead held. But stage 1 registered its own caveat for
exactly this stratum and never tested it:

    "Most worm transcripts are trans-spliced, and the outron upstream of the annotated 5' end is
    transcribed with the gene, so an element between genes that the rule attributes to the gene it
    sits upstream of can also lie in that gene's transcribed sequence."

THIS MODULE IS THAT TEST. It asks whether the between-gene lead survives removing the elements that
the caveat names: those in a gene's outron, and those whose attribution is confounded by operon
transcription. Nothing else about stage 1 changes -- the attribution, the comparator (nearest
annotated 5' end), the endpoint, the 1-Mb cluster bootstrap and the reading vocabulary are imported
from `simpler.py`, not restated.

PRIOR EXPOSURE, STATED PLAINLY. Serizay et al. 2020's RNA was read in stage 1. This is a re-analysis
of a measurement this project has already read, not a fresh test. The exclusion masks are built from
positions and annotations only -- WormBase WS298 operon, trans-splice-acceptor and TSS-region
features, and the Ensembl WBcel235 gene models -- and the population they leave is fixed and
committed before `score` opens the RNA sheet again. That removes the chance of tuning the mask to
the answer; it does not make the RNA unseen.

The map (fetched 2026-10-10, no model request, no key, no money):

    c_elegans.PRJNA13758.WS298.annotations.gff3.gz, WormBase release WS298, streamed from the EBI
    mirror and filtered in one pass to the six feature classes this lane needs. The streamed file's
    md5 equals the md5 WormBase publishes for it in WS298/CHECKSUMS, and two independent streams gave
    the same sha256, so the filtered extract is of the authentic complete release file.

    operon|operon              1,385  WormBase operons, with their member WBGene ids
    dicistronic_mRNA|operon       46  dicistronic messages, treated as operons here
    deprecated_operon|operon     111  NOT USED
    SL1|SL1_acceptor_site     43,207  measured SL1 trans-splice acceptor sites
    SL2|SL2_acceptor_site     13,444  measured SL2 trans-splice acceptor sites
    TSS_region|TSS_region     73,499  WormBase's clustered transcription start site regions

    The TSS regions are WormBase's own annotation of the four published C. elegans TSS data sets --
    the WS298 release letter names papers WBPaper00042246, WBPaper00042529 and WBPaper00042354 for
    the TSS wiggles, and the wiggle archive's members are TSS_Chen, TSS_Kruesi, TSS_Saito and TSS_Gu,
    i.e. Chen et al. 2013, Kruesi et al. 2013, Saito et al. 2013 and Gu et al. Saito et al. 2013 is
    the outron paper: it enriched for transcripts that still carried their outrons in order to place
    TSSs for trans-spliced genes, which is the measurement this lane needs.

    Assembly: WBcel235 (ce11). WormBase has used this assembly since WS235, so the WS298 features and
    the Ensembl WBcel235 gene models stage 1 read are on the same build and NOTHING IS LIFTED OVER.

    Licence: WormBase data are distributed for free reuse; this lane records the release, the URL and
    the published checksum rather than asserting a licence identifier it did not read, because the
    WormBase policy page refused an unauthenticated fetch (HTTP 403, Cloudflare).

What this module cannot establish is unchanged from stage 1. Accessibility is not function, and
conservation is not function. A concordance says the attributed gene is expressed where the element
is open, never that the element regulates the gene. Removing outron-affected elements removes one
named confound on one map; it does not turn the remainder into measured regulation.
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from genomeos.attribution import simpler
from genomeos.attribution.simpler import (
    ALPHA,
    CLUSTER_BP,
    ELEMENT_TYPE,
    TISSUES,
    Chromosome,
    ControlIndex,
    cluster_bootstrap,
    compare,
    control_test,
    distance_bin,
)

STORE = Path("data/knowledge/wormoutron")

#: The filtered WS298 extract: the six feature classes above, GFF3 columns unchanged.
MAP_FILE = STORE / "ws298_operon_sl_tss.gff3"
#: The whole-file source|type census written in the same single pass, so the filter can be audited.
CENSUS_FILE = STORE / "ws298_source_type_census.txt"
#: WormBase's published checksums for WS298, as served; the gff3's md5 line is the provenance.
CHECKSUMS_FILE = STORE / "ws298_CHECKSUMS.txt"

MAP_SOURCE = (
    "WormBase release WS298, c_elegans.PRJNA13758.WS298.annotations.gff3.gz, "
    "https://ftp.ebi.ac.uk/pub/databases/wormbase/releases/WS298/species/c_elegans/PRJNA13758/, "
    "published md5 f41b0e0afb7bdcf15d7b35c9a1f7f3a9 (WS298/CHECKSUMS), "
    "streamed sha256 4c84b835c5eba84d2128947ae5646329c0172e25c29b40b416b9bf7807146030"
)
REMOTE_MD5 = "f41b0e0afb7bdcf15d7b35c9a1f7f3a9"
REMOTE_SHA256 = "4c84b835c5eba84d2128947ae5646329c0172e25c29b40b416b9bf7807146030"

OPERON_SOURCES = ("operon", "dicistronic_mRNA")
OPERON_SOURCES_NOT_USED = ("deprecated_operon",)
SL_TYPES = ("SL1_acceptor_site", "SL2_acceptor_site")
TSS_TYPE = "TSS_region"

#: Beyond this, no TSS region upstream of a gene's annotated 5' end is taken to be that gene's, and
#: the gene is given no outron rather than an outron of the cap's length. Fixed before the census ran.
OUTRON_MAX_BP = 20_000
#: An SL1 or SL2 acceptor site this close to a gene's annotated 5' end is evidence that the gene is
#: trans-spliced there. Fixed before the census ran.
SL_WINDOW_BP = 100

#: Stage 1's between-gene stratum, as published in docs/SIMPLER-GENOMES.md line 124. The comparison
#: this lane reads against, not a quantity it recomputes loosely.
STAGE1_BETWEEN_GENES = {
    "elements": 2335,
    "domain_rate_pct": 53.9,
    "nearest_rate_pct": 50.2,
    "delta_points": 3.73,
    "delta_ci95_points": [2.79, 4.74],
    "registration": "595d7c2 / eede7d9, docs/SIMPLER-GENOMES.md sections 3 and 4",
}

#: Fixed BEFORE the census was computed and not changed afterwards. Below either floor the run is
#: INCONCLUSIVE and unscored. With n surviving elements and k of them where the rules differ, the
#: smallest difference at 80% power and two-sided 5% is about 2.8*sqrt(k)/n before any inflation for
#: clustering; n = 800 with k = 100 puts that at about 3.5 points, which is the smallest floor at
#: which stage 1's +3.73 could be seen at all.
MIN_SURVIVING = 800
MIN_SURVIVING_DISCORDANT = 100


@dataclass(frozen=True)
class Interval:
    chrom: str
    start: int
    end: int
    strand: str
    name: str


@dataclass(frozen=True)
class Operon:
    chrom: str
    start: int
    end: int
    strand: str
    name: str
    members: tuple[str, ...]
    source: str


def _attrs(field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in field.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


def read_map(path: Path = MAP_FILE) -> dict[str, Any]:
    """The WS298 extract, parsed. Positions and strands only; no RNA, no accessibility."""
    operons: list[Operon] = []
    deprecated = 0
    sl: list[Interval] = []
    tss: list[Interval] = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        # The WS298 GFF3 names the sequences I..V, X, MtDNA; simpler.read_genes and Table S2 both
        # name them chrI..chrX, so the map is renamed here and nowhere else.
        chrom = f[0] if f[0].startswith("chr") else "chr" + f[0]
        source, kind = f[1], f[2]
        start, end, strand = int(f[3]), int(f[4]), f[6]
        a = _attrs(f[8])
        name = a.get("Name", "")
        if kind == "operon":
            if source in OPERON_SOURCES_NOT_USED:
                deprecated += 1
                continue
            if source not in OPERON_SOURCES:
                continue
            members = tuple(g for g in a.get("genes", "").split(",") if g)
            operons.append(Operon(chrom, start, end, strand, name, members, source))
        elif kind in SL_TYPES:
            sl.append(Interval(chrom, start, end, strand, name))
        elif kind == TSS_TYPE:
            tss.append(Interval(chrom, start, end, strand, name))
    return {
        "operons": operons,
        "deprecated_operons_skipped": deprecated,
        "sl_sites": sl,
        "tss_regions": tss,
    }


class TssIndex:
    """TSS regions of one chromosome and strand, in coordinate order."""

    def __init__(self, regions: list[Interval]):
        self.by_start = sorted(regions, key=lambda r: (r.start, r.end))
        self.starts = [r.start for r in self.by_start]
        self.by_end = sorted(regions, key=lambda r: (r.end, r.start))
        self.ends = [r.end for r in self.by_end]

    def last_at_or_before(self, pos: int) -> Interval | None:
        i = bisect.bisect_right(self.starts, pos)
        return self.by_start[i - 1] if i > 0 else None

    def first_at_or_after(self, pos: int) -> Interval | None:
        i = bisect.bisect_left(self.ends, pos)
        return self.by_end[i] if i < len(self.by_end) else None


def _strand_index(items: list[Interval]) -> dict[tuple[str, str], TssIndex]:
    groups: dict[tuple[str, str], list[Interval]] = {}
    for it in items:
        groups.setdefault((it.chrom, it.strand), []).append(it)
    return {k: TssIndex(v) for k, v in groups.items()}


def sl_spliced_genes(
    chroms: dict[str, Chromosome], sl_sites: list[Interval], window: int = SL_WINDOW_BP
) -> set[str]:
    """Genes with an SL1 or SL2 acceptor site within `window` bp of the annotated 5' end, same strand.

    Measured trans-splicing at the gene's mature 5' end. This is what makes an upstream outron
    possible at all, so it qualifies the narrower of the two outron definitions below.
    """
    idx: dict[tuple[str, str], list[int]] = {}
    for s in sl_sites:
        idx.setdefault((s.chrom, s.strand), []).append(s.start)
    for key in idx:
        idx[key].sort()
    out: set[str] = set()
    for c, ch in chroms.items():
        for g in ch.genes:
            pos = idx.get((c, g.strand))
            if not pos:
                continue
            i = bisect.bisect_left(pos, g.tss - window)
            if i < len(pos) and pos[i] <= g.tss + window:
                out.add(g.gid)
    return out


def outrons(
    chroms: dict[str, Chromosome],
    tss_regions: list[Interval],
    max_bp: int = OUTRON_MAX_BP,
) -> dict[str, tuple[int, int]]:
    """Per gene, the interval between its annotated 5' end and the TSS region upstream of it.

    The outron is the transcribed sequence that trans-splicing removes from the message, so an
    element inside it lies in the gene's transcript without lying in the gene as annotated. The
    definition is structural and was fixed before the census: the nearest TSS region upstream of the
    gene's annotated 5' end on the SAME strand; truncated so the interval does not reach past another
    protein-coding gene's annotated 5' end on that strand; and no outron at all if the nearest such
    TSS region is more than `max_bp` away, rather than an outron of the cap's length.
    """
    idx = _strand_index(tss_regions)
    out: dict[str, tuple[int, int]] = {}
    for c, ch in chroms.items():
        same = {
            "+": sorted(g.tss for g in ch.genes if g.strand == "+"),
            "-": sorted(g.tss for g in ch.genes if g.strand == "-"),
        }
        for g in ch.genes:
            ti = idx.get((c, g.strand))
            if ti is None:
                continue
            t = g.tss
            if g.strand == "+":
                r = ti.last_at_or_before(t)
                if r is None or r.start >= t:
                    continue
                lo, hi = r.start, t
                if t - lo > max_bp:
                    continue
                others = same["+"]
                k = bisect.bisect_left(others, t) - 1
                while k >= 0 and others[k] == t:
                    k -= 1
                if k >= 0 and others[k] > lo:
                    lo = others[k] + 1
            else:
                r = ti.first_at_or_after(t)
                if r is None or r.end <= t:
                    continue
                lo, hi = t, r.end
                if hi - t > max_bp:
                    continue
                others = same["-"]
                j = bisect.bisect_right(others, t)
                while j < len(others) and others[j] == t:
                    j += 1
                if j < len(others) and others[j] < hi:
                    hi = others[j] - 1
            if lo <= hi:
                out[g.gid] = (lo, hi)
    return out


def operon_members(operons: list[Operon], chroms: dict[str, Chromosome]) -> dict[str, dict[str, Any]]:
    """Per operon member we can place, whether it is the 5'-most member of its operon.

    A non-5'-most member has no promoter of its own: it is transcribed from the operon's single
    upstream promoter and SL2-trans-spliced out of the polycistron. An element attributed to such a
    gene because it sits upstream of it is therefore attributed across an intercistronic region, which
    is the second half of the caveat. Where fewer than two members can be placed in the gene models
    the order is unknown, and the member is marked downstream, which excludes it -- the conservative
    choice, fixed before the census.
    """
    pos: dict[str, Any] = {}
    for ch in chroms.values():
        for g in ch.genes:
            pos[g.gid] = g
    out: dict[str, dict[str, Any]] = {}
    for op in operons:
        placed = [pos[m] for m in op.members if m in pos]
        if len(placed) >= 2:
            first = (
                min(placed, key=lambda g: g.start) if op.strand == "+" else max(placed, key=lambda g: g.end)
            )
            order_known = True
        else:
            first, order_known = None, False
        for g in placed:
            out[g.gid] = {
                "operon": op.name,
                "source": op.source,
                "order_known": order_known,
                "is_first": bool(first is not None and g.gid == first.gid),
                "members_total": len(op.members),
                "members_placed": len(placed),
                "agrees_with_list_order": bool(
                    first is not None and op.members and first.gid == op.members[0]
                ),
            }
    return out


def operon_spans(operons: list[Operon]) -> dict[str, list[tuple[int, int]]]:
    spans: dict[str, list[tuple[int, int]]] = {}
    for op in operons:
        spans.setdefault(op.chrom, []).append((op.start, op.end))
    for c in spans:
        spans[c].sort()
    return spans


class SpanIndex:
    """Closed intervals of one chromosome, with a containment test that scans no more than it must."""

    def __init__(self, intervals: list[tuple[int, int]]):
        self.ivs = sorted(intervals)
        self.starts = [a for a, _ in self.ivs]
        self.max_end: list[int] = []
        running = -1
        for _, b in self.ivs:
            running = max(running, b)
            self.max_end.append(running)

    def covers(self, pos: int) -> bool:
        i = bisect.bisect_right(self.starts, pos) - 1
        while i >= 0:
            if self.max_end[i] < pos:
                return False
            if self.ivs[i][1] >= pos:
                return True
            i -= 1
        return False


def _in_outron(outron: dict[str, tuple[int, int]], gid: str, pos: int) -> bool:
    iv = outron.get(gid)
    return iv is not None and iv[0] <= pos <= iv[1]


def flag(
    rows: list[dict[str, Any]],
    chroms: dict[str, Chromosome],
    mapdata: dict[str, Any],
) -> dict[str, Any]:
    """Mark each element with every mask. No RNA value and no endpoint is read here."""
    outron = outrons(chroms, mapdata["tss_regions"])
    spliced = sl_spliced_genes(chroms, mapdata["sl_sites"])
    members = operon_members(mapdata["operons"], chroms)
    spans = {c: SpanIndex(v) for c, v in operon_spans(mapdata["operons"]).items()}
    gene_chrom = {g.gid: c for c, ch in chroms.items() for g in ch.genes}
    per_chrom: dict[str, list[tuple[int, int]]] = {}
    for gid, iv in outron.items():
        per_chrom.setdefault(gene_chrom[gid], []).append(iv)
    any_idx = {c: SpanIndex(v) for c, v in per_chrom.items()}

    def downstream(gid: str) -> bool:
        m = members.get(gid)
        return m is not None and not m["is_first"]

    for r in rows:
        mid, c = r["mid"], r["chrom"]
        dom, near = r["domain_gene"], r["nearest_gene"]
        r["outron_domain"] = int(_in_outron(outron, dom, mid))
        r["outron_nearest"] = int(_in_outron(outron, near, mid))
        r["outron_either"] = int(bool(r["outron_domain"] or r["outron_nearest"]))
        r["outron_either_sl"] = int(
            bool((r["outron_domain"] and dom in spliced) or (r["outron_nearest"] and near in spliced))
        )
        oi = any_idx.get(c)
        r["outron_any_gene"] = int(oi is not None and oi.covers(mid))
        si = spans.get(c)
        r["operon_span"] = int(si is not None and si.covers(mid))
        r["operon_downstream_domain"] = int(downstream(dom))
        r["operon_downstream_nearest"] = int(downstream(near))
        r["operon_mask"] = int(
            bool(r["operon_span"] or r["operon_downstream_domain"] or r["operon_downstream_nearest"])
        )
        r["excluded"] = int(bool(r["outron_either"] or r["operon_mask"]))
        r["excluded_strict"] = int(bool(r["outron_any_gene"] or r["operon_mask"]))
        r["excluded_sl"] = int(bool(r["outron_either_sl"] or r["operon_mask"]))
    return {
        "outron_genes": len(outron),
        "sl_spliced_genes": len(spliced),
        "operon_member_genes": len(members),
        "operon_members_order_unknown": sum(1 for m in members.values() if not m["order_known"]),
        "operon_members_downstream": sum(1 for m in members.values() if not m["is_first"]),
        "operon_first_agrees_with_list_order": sum(
            1 for m in members.values() if m["agrees_with_list_order"]
        ),
        "outron_length_bp": _spread([hi - lo + 1 for lo, hi in outron.values()]),
    }


def _spread(values: list[int]) -> dict[str, Any]:
    if not values:
        return {"n": 0}
    v = sorted(values)
    return {
        "n": len(v),
        "min": v[0],
        "p25": v[len(v) // 4],
        "median": v[len(v) // 2],
        "p75": v[(3 * len(v)) // 4],
        "p90": v[int(len(v) * 0.9)],
        "max": v[-1],
    }


def between_genes(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Stage 1's between-gene stratum: single-tissue putative enhancers with no host gene."""
    return [r for r in rows if r["kind"] == ELEMENT_TYPE and r["domain_case"] != "host"]


def census(rows: list[dict[str, Any]], mapstats: dict[str, Any]) -> dict[str, Any]:
    """The population the masks leave, from positions and annotations only. Fixes what `score` reads."""
    bg = between_genes(rows)
    kept = [r for r in bg if not r["excluded"]]
    disc = [r for r in kept if r["domain_gene"] != r["nearest_gene"]]
    n, k = len(kept), len(disc)
    mde = 2.8 * (k**0.5) / n if n else None
    return {
        "stage1_between_genes": STAGE1_BETWEEN_GENES,
        "between_genes_recomputed": len(bg),
        "matches_stage1": len(bg) == STAGE1_BETWEEN_GENES["elements"],
        "map": mapstats,
        "excluded_by": {
            "outron_of_either_rule_gene": sum(r["outron_either"] for r in bg),
            "outron_of_domain_rule_gene": sum(r["outron_domain"] for r in bg),
            "outron_of_nearest_tss_gene": sum(r["outron_nearest"] for r in bg),
            "outron_with_sl_evidence": sum(r["outron_either_sl"] for r in bg),
            "outron_of_any_gene": sum(r["outron_any_gene"] for r in bg),
            "inside_an_operon_span": sum(r["operon_span"] for r in bg),
            "domain_gene_is_a_downstream_operon_member": sum(r["operon_downstream_domain"] for r in bg),
            "nearest_gene_is_a_downstream_operon_member": sum(r["operon_downstream_nearest"] for r in bg),
            "operon_mask": sum(r["operon_mask"] for r in bg),
            "primary_union": sum(r["excluded"] for r in bg),
            "strict_union": sum(r["excluded_strict"] for r in bg),
            "sl_union": sum(r["excluded_sl"] for r in bg),
        },
        "surviving": n,
        "surviving_fraction": round(n / len(bg), 4) if bg else None,
        "surviving_discordant": k,
        "surviving_clusters_1mb": len({(r["chrom"], r["mid"] // CLUSTER_BP) for r in kept}),
        "surviving_by_tissue": {t: sum(1 for r in kept if r["tissue"] == t) for t in TISSUES},
        "surviving_by_case": _by_case(kept),
        "minimum_surviving": MIN_SURVIVING,
        "minimum_surviving_discordant": MIN_SURVIVING_DISCORDANT,
        "mde_points_unclustered": round(100 * mde, 3) if mde is not None else None,
        "inconclusive_unscored": n < MIN_SURVIVING or k < MIN_SURVIVING_DISCORDANT,
        "strict_surviving": sum(1 for r in bg if not r["excluded_strict"]),
        "sl_surviving": sum(1 for r in bg if not r["excluded_sl"]),
    }


def _by_case(rows: list[dict[str, Any]]) -> dict[str, int]:
    out: dict[str, int] = {}
    for r in rows:
        out[r["domain_case"]] = out.get(r["domain_case"], 0) + 1
    return out


def survival_label(point: float, lo: float, hi: float) -> str:
    """The registered decision. `point` and the bounds are in points, the bootstrap's own units."""
    if lo > 0 and point >= STAGE1_BETWEEN_GENES["delta_ci95_points"][0]:
        return "LEAD_SURVIVES"
    if lo > 0:
        return "LEAD_ATTENUATED"
    if hi < 0:
        return "NEAREST_AHEAD_AFTER_EXCLUSION"
    return "LEAD_DOES_NOT_SURVIVE"


#: Written and committed before `score` re-read the RNA sheet. Not to be changed afterwards.
REGISTRATION: dict[str, Any] = {
    "question": (
        "Does the worm's between-gene lead of the domain rule over nearest TSS survive EXCLUDING the "
        "elements that fall in outrons or in trans-splicing-affected 5' regions by the WormBase WS298 "
        "map -- the confound stage 1 registered as a caveat for this very stratum and never tested?"
    ),
    "prior_exposure": (
        "Serizay et al. 2020's RNA was read in stage 1 (registration 595d7c2, result eede7d9, "
        "docs/SIMPLER-GENOMES.md sections 3 and 4), so THIS IS A RE-ANALYSIS OF A MEASUREMENT THIS "
        "PROJECT HAS ALREADY READ, NOT A FRESH TEST. What is new is the WormBase WS298 outron and "
        "operon map, which no earlier run read. The masks are built from positions and annotations "
        "alone, their definitions were fixed before the census was computed, and the surviving "
        "population is counted and committed before the RNA sheet is opened again; that bars tuning "
        "the mask to the answer, and it does not make the RNA unseen. The reading must be read as a "
        "re-analysis throughout."
    ),
    "map": MAP_SOURCE,
    "map_defines": (
        "operon|operon and dicistronic_mRNA|operon: operon spans with member WBGene ids "
        "(deprecated_operon is not used); SL1_acceptor_site and SL2_acceptor_site: measured "
        "trans-splice acceptor positions; TSS_region: WormBase's clustered transcription start site "
        "regions from the four published C. elegans TSS data sets (TSS_Chen, TSS_Kruesi, TSS_Saito, "
        "TSS_Gu in the WS298 TSS wiggle archive; WS298's release letter names WBPaper00042246, "
        "WBPaper00042529 and WBPaper00042354). Assembly WBcel235 (ce11), the same build as the "
        "Ensembl gene models stage 1 read; nothing is lifted over."
    ),
    "population": (
        "Stage 1's between-gene stratum: the Serizay et al. 2020 Table S2 putative enhancers "
        "accessible in exactly one of the five tissues whose midpoint lies in no protein-coding gene "
        "(domain_case != 'host'), recomputed by simpler.attribute from the same inputs and required to "
        "come to stage 1's 2,335 elements. EXCLUSION, applied to that stratum and nothing else: an "
        "element is removed if (a) its midpoint lies in the outron of the gene the domain rule "
        "attributes OR of the nearest-5'-end gene, or (b) its midpoint lies inside a WormBase operon "
        "span, or (c) either rule's gene is an operon member that is not the 5'-most member of its "
        "operon (order unknown counts as not 5'-most). A gene's outron is the interval between its "
        "annotated 5' end and the nearest TSS region upstream of it on the same strand, truncated so "
        "it does not reach past another protein-coding gene's annotated 5' end on that strand, and "
        "absent altogether if that TSS region is more than 20,000 bp away. The population the "
        "exclusion leaves is fixed by the census in this registration, before any RNA value is read."
    ),
    "attribution": (
        "Unchanged and imported, not restated: simpler.domain_rule. " + simpler.REGISTRATION["attribution"]
    ),
    "baseline": simpler.REGISTRATION["baseline"],
    "endpoint": simpler.REGISTRATION["endpoint"],
    "statistic": simpler.REGISTRATION["statistic"],
    "matched_control": simpler.REGISTRATION["matched_control"],
    "decision": (
        "On the surviving population: LEAD_SURVIVES if the 95% interval's lower bound > 0 AND the "
        "point estimate is at or above stage 1's between-gene lower bound of +2.79 points; "
        "LEAD_ATTENUATED if the lower bound > 0 but the point estimate is below +2.79; "
        "NEAREST_AHEAD_AFTER_EXCLUSION if the upper bound < 0; LEAD_DOES_NOT_SURVIVE otherwise, i.e. "
        "when the interval includes 0."
    ),
    "falsifier": (
        "The claim under test -- that the worm's between-gene lead is not an artefact of outron and "
        "operon transcription -- is falsified if the excluded population's interval crosses 0, or if "
        "its point estimate falls below stage 1's lower bound of +2.79 points. The first reads "
        "LEAD_DOES_NOT_SURVIVE, the second LEAD_ATTENUATED; both are stated first and plainly, and "
        "the worm would then join the fly."
    ),
    "inconclusive": (
        "INCONCLUSIVE means the run cannot decide either way and is reported unscored, with no "
        "comparison read. It occurs when too few elements survive the exclusion: fewer than 800 "
        "surviving elements, or fewer than 100 of the survivors on which the two rules name different "
        "genes. Both floors were fixed before the census was computed. The reason for 800 and 100: "
        "the smallest difference at 80% power and two-sided 5% is about 2.8*sqrt(k)/n points before "
        "any inflation for clustering, which at n = 800 and k = 100 is about 3.5 points -- the "
        "smallest population in which stage 1's +3.73 could be seen at all. A run whose interval is "
        "wide but whose floors are met is NOT inconclusive: it reads LEAD_DOES_NOT_SURVIVE, because a "
        "lead that an adequate population cannot confirm is not a lead that survived."
    ),
    "strata_beside": (
        "Reported beside the primary, deciding nothing: the excluded elements themselves; the "
        "strict exclusion (the outron of ANY protein-coding gene, union the operon mask); the "
        "SL-restricted exclusion (outrons only of genes with a measured SL1 or SL2 acceptor site "
        "within 100 bp of the annotated 5' end, union the operon mask); the outron mask and the "
        "operon mask each alone; each tissue among the survivors; and the survivors' domain-rule genes "
        "against their own matched-distance control."
    ),
    "stop_rule": (
        "One run. No mask, cap, floor, population, threshold or endpoint is changed after the RNA "
        "sheet is read; a second definition would be a new registration, committed before it is run. "
        "If the run is re-executed in a worktree, this registration is committed there before the run "
        "too, and the two are never merged."
    ),
    "what_it_cannot_establish": (
        "Accessibility is not function and conservation is not function. A concordance says the "
        "attributed gene is expressed where the element is open, never that the element regulates the "
        "gene. Survival of this exclusion says only 'survives exclusion of outron-affected elements on "
        "this map' -- one map, one stage (young adults), five tissues, one source, and a re-analysis of "
        "RNA this project has already read. It does not separate regulation from transcription in "
        "general, and it says nothing about the host-gene stratum, which is not the population here."
    ),
    "outcome_wording": {
        "head": (
            "On the {n} of stage 1's 2,335 C. elegans between-gene putative enhancers that survive "
            "excluding outron- and operon-affected elements by the WormBase WS298 map, the gene the "
            "domain rule attributes is expressed at 5 TPM or more and highest in the element's own "
            "tissue for {kd} of {n} ({cd}%), against {kn} of {n} ({cn}%) for the gene with the nearest "
            "5' end; the rules differ on {k} of the survivors; difference {d} points, 95% {lo} to {hi} "
            "over 1-Mb clusters, against stage 1's +3.73 [+2.79, +4.74] on the unexcluded 2,335."
        ),
        "LEAD_SURVIVES": (
            "The worm's between-gene lead survives exclusion of outron-affected elements on this map."
        ),
        "LEAD_ATTENUATED": (
            "The worm's between-gene lead is attenuated by the exclusion: it is still above nearest "
            "TSS, but the point estimate falls below stage 1's lower bound, so the registered "
            "falsifier is met and the lead is NOT read as having survived."
        ),
        "LEAD_DOES_NOT_SURVIVE": (
            "The worm's between-gene lead does NOT survive the exclusion: with outron- and "
            "operon-affected elements removed, the interval includes 0 and no difference between the "
            "domain rule and nearest TSS is detected. On this map the worm joins the fly, where the "
            "lead also appears only where a gene's own transcription can explain it."
        ),
        "NEAREST_AHEAD_AFTER_EXCLUSION": (
            "After the exclusion nearest TSS is ahead of the domain rule between genes, so the worm's "
            "between-gene lead does not survive and reverses on this map."
        ),
        "INCONCLUSIVE": (
            "The exclusion leaves too few elements to decide: {n} survivors and {k} on which the rules "
            "differ, against the registered floors of 800 and 100. The run is unscored and no "
            "comparison is read."
        ),
        "control_informative": (
            "Among the survivors the domain rule's genes match the element's tissue more often than "
            "genes at a matched distance ({cc}%, one-sided p {p}), so position still carries tissue "
            "information in the population the exclusion leaves."
        ),
        "control_uninformative": (
            "Among the survivors the domain rule's genes match the element's tissue no more often than "
            "genes at a matched distance ({cc}%, one-sided p {p}): in the population the exclusion "
            "leaves, the attribution carries no tissue information beyond distance."
        ),
        "exposure": (
            "This is a re-analysis of RNA this project has already read in stage 1, not a fresh test."
        ),
        "tail": (
            "Accessibility is not function: this says the attributed gene is expressed where the "
            "element is open, not that the element acts on it."
        ),
    },
}


def _fmt_p(p: float) -> str:
    return f"{p:.3g}"


def reading(res: dict[str, Any]) -> str:
    """The registered outcome wording, filled in. The templates are REGISTRATION['outcome_wording']."""
    w = REGISTRATION["outcome_wording"]
    cen = res["census"]
    if res["label"] == "INCONCLUSIVE":
        return " ".join(
            [
                w["INCONCLUSIVE"].format(n=f"{cen['surviving']:,}", k=f"{cen['surviving_discordant']:,}"),
                w["exposure"],
            ]
        )
    p = res["primary"]
    ctl = res["control"]
    fill = {
        "n": f"{p['n']:,}",
        "kd": f"{p['domain_concordant']:,}",
        "kn": f"{p['nearest_concordant']:,}",
        "cd": f"{p['domain_rate_pct']:.1f}",
        "cn": f"{p['nearest_rate_pct']:.1f}",
        "d": f"{p['delta_points']:+.2f}",
        "lo": f"{p['delta_ci95_points'][0]:+.2f}",
        "hi": f"{p['delta_ci95_points'][1]:+.2f}",
        "k": f"{cen['surviving_discordant']:,}",
        "cc": f"{100 * ctl['expected_rate']:.1f}",
        "p": _fmt_p(ctl["p_one_sided"]),
    }
    parts = [w["head"].format(**fill), w[res["label"]].format(**fill)]
    key = "control_informative" if ctl["p_one_sided"] < ALPHA else "control_uninformative"
    parts.append(w[key].format(**fill))
    parts.append(w["exposure"])
    parts.append(w["tail"])
    return " ".join(parts)


def score(
    rows: list[dict[str, Any]],
    chroms: dict[str, Chromosome],
    expr: dict[str, tuple[float, ...]],
    cen: dict[str, Any],
) -> dict[str, Any]:
    """One run on the committed census's surviving population. The first RNA read of this lane."""
    bg = between_genes(rows)
    kept = [r for r in bg if not r["excluded"]]
    out: dict[str, Any] = {"census": cen}
    if cen["inconclusive_unscored"]:
        out["label"] = "INCONCLUSIVE"
        out["primary"] = None
        out["control"] = None
        out["strata"] = {}
        out["reading"] = reading(out)
        return out
    primary = compare(list(kept), expr)
    out["primary"] = primary
    label = survival_label(
        primary["delta_points"], primary["delta_ci95_points"][0], primary["delta_ci95_points"][1]
    )
    out["label"] = label
    if label == "LEAD_SURVIVES":
        out["falsifier"] = "NOT_FALSIFIED"
    elif label == "LEAD_ATTENUATED":
        out["falsifier"] = "FALSIFIED: the point estimate falls below stage 1's lower bound of +2.79"
    else:
        out["falsifier"] = "FALSIFIED: the interval does not keep the lead above 0"
    idx = ControlIndex(chroms, expr)
    probs, observed = [], 0
    for r in kept:
        probs.append(idx.probability(r["chrom"], r["tissue"], r["mid"], r["domain_gene"], r["domain_dist"]))
        observed += r["c_domain"]
    out["control"] = control_test(probs, observed)
    out["strata"] = _strata(bg, expr)
    out["comparison_to_stage1"] = {
        "stage1": STAGE1_BETWEEN_GENES,
        "point_minus_stage1_point": round(primary["delta_points"] - STAGE1_BETWEEN_GENES["delta_points"], 3),
        "at_or_above_stage1_lower_bound": primary["delta_points"]
        >= STAGE1_BETWEEN_GENES["delta_ci95_points"][0],
    }
    out["reading"] = reading(out)
    return out


def _subset(rows: list[dict[str, Any]], expr: dict[str, tuple[float, ...]]) -> dict[str, Any]:
    if not rows:
        return {"n": 0}
    res = compare(list(rows), expr)
    if res["delta_ci95_points"] is not None:
        res["label"] = simpler.comparison_label(*res["delta_ci95_points"])
    return res


def _strata(bg: list[dict[str, Any]], expr: dict[str, tuple[float, ...]]) -> dict[str, Any]:
    kept = [r for r in bg if not r["excluded"]]
    out: dict[str, Any] = {
        "all_between_genes_unexcluded": _subset(bg, expr),
        "excluded_elements": _subset([r for r in bg if r["excluded"]], expr),
        "strict_exclusion_survivors": _subset([r for r in bg if not r["excluded_strict"]], expr),
        "sl_restricted_exclusion_survivors": _subset([r for r in bg if not r["excluded_sl"]], expr),
        "outron_mask_only_survivors": _subset([r for r in bg if not r["outron_either"]], expr),
        "operon_mask_only_survivors": _subset([r for r in bg if not r["operon_mask"]], expr),
        "surviving_discordant_only": _subset(
            [r for r in kept if r["domain_gene"] != r["nearest_gene"]], expr
        ),
    }
    for t in TISSUES:
        out[f"surviving_{t}"] = _subset([r for r in kept if r["tissue"] == t], expr)
    return out


__all__ = [
    "MAP_FILE",
    "MAP_SOURCE",
    "MIN_SURVIVING",
    "MIN_SURVIVING_DISCORDANT",
    "OUTRON_MAX_BP",
    "REGISTRATION",
    "STAGE1_BETWEEN_GENES",
    "SpanIndex",
    "between_genes",
    "census",
    "cluster_bootstrap",
    "distance_bin",
    "flag",
    "operon_members",
    "operon_spans",
    "outrons",
    "read_map",
    "reading",
    "score",
    "sl_spliced_genes",
    "survival_label",
]
