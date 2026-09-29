# SPDX-License-Identifier: AGPL-3.0-or-later
"""Hide evidence, predict it: the leave-one-source-out harness (item 13 C4, registered 2026-09-28).

Every later approach of the coherence programme (C1 the debugger, C2 the boundaries, C3 the families)
is scored by this module and by nothing else, so its rules are fixed here, before any score was read,
and changing one is a new registration rather than an edit.

**What it asks.** One measured source is held out; a labelling built without it predicts that
source's own endpoint; the prediction is scored with an interval over loci. Rotating the held-out
source asks whether coherence transfers across measurements: whether sequence, chromatin and the
other assays predict what a CRISPRi screen, a reporter, a transgenic embryo or a population measured.

**What it is not.** An internal development benchmark only (`STATUS`). Every source here has been read
by this project before, several of them repeatedly: the ENCODE CRISPRi held-out file alone is read by
at least ten scored results (the split audit of 2026-09-28, lane-split), so no score this module
returns is a fresh or external validation, and every score says so (`MAY_NOT_BE_CALLED`).

**The pieces, each a registered constant.**

- Sources (`SOURCE_KINDS`): each CRISPRi study (the benchmark's `Dataset` field, R5), each lentiMPRA
  cell, VISTA, saturation mutagenesis and GTEx eQTL, each a separately holdable unit.
- Endpoints (`ENDPOINTS`): each assay's own, never pooled. CRISPRi keeps R2's outcomes apart: a
  significant decrease and a significant increase are two endpoints, each against the well-powered
  nulls, and underpowered nulls and missing effects are counted, never scored. lentiMPRA is read per
  tile (R6), never aggregated. Activity is not regulation, association is not perturbation, and no
  endpoint substitutes for another.
- The split (`evidence`): by study, and by locus. Holding out S removes S and its provenance siblings;
  every record of S's endpoint family whose interval shares a base with an S unit is masked as well
  (lane-split's `split_overlap` relation, "related intervals"), so an element measured by two studies
  cannot leak across the split; and the CRISPRi benchmark's held-out file is never evidence at all
  (R5, `measured.development_only`).
- The metric (`PRIMARY`, `SECONDARY`) and its interval: a percentile bootstrap over loci, never over
  pairs (`LOCUS_RULE`), with coverage and abstentions reported beside every value.
- Baselines (`LABELLINGS`): distance to TSS, and the unchanged labels (the compiled programs'
  predicted layer as it stands). The reference labelling `rest` is the prediction from the remaining
  sources, an unfitted rule (`REST_RULE`).

`score(labels, held_out_source)` is the call the pilot makes. It is deterministic (fixed seed, fixed
resample count), cached (in memory and under data/cache/holdout), and cheap: numpy over at most a few
hundred thousand units, no model request, no network.
"""

from __future__ import annotations

import bisect
import csv
import gzip
import hashlib
import json
import math
import re
from collections import defaultdict
from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any

from genomeos.attribution import measured as ms
from genomeos.attribution import mpra, vista
from genomeos.attribution.crispri import KNOWLEDGE as CRISPRI_KNOWLEDGE

# ==================================================================================================
# The registration (2026-09-28, item 13 C4, lane-c4). Written before any score of any labelling on any
# source was computed; docs/ATTRIBUTION.md carries the same text under its dated heading.
# ==================================================================================================
REGISTERED = "2026-09-28"
STATUS = "internal development benchmark"
REUSE = (
    "every source is evidence this project has already read, several repeatedly; the ENCODE CRISPRi "
    "held-out file is read by at least ten scored results (lane-split, 2026-09-28) and is a reused "
    "benchmark. No score from this harness is a fresh or external validation"
)
MAY_BE_CALLED = (
    "an internal development benchmark reading: on evidence this project has already read, labelling L "
    "predicts held-out source S's own endpoint with metric M [95% locus-bootstrap interval] at coverage "
    "C, beside distance to TSS and the unchanged labels"
)
MAY_NOT_BE_CALLED = (
    "an external, independent or fresh validation",
    "evidence of biological accuracy beyond the sources held here",
    "a test of regulation when the endpoint is reporter activity, in-vivo activity, base sensitivity "
    "or association (activity is not regulation; association is not perturbation)",
    "a comparison between sources or endpoints (different units, prevalences and assays)",
    "out of sample for anything tuned on the CRISPRi benchmark, whose held-out file is a reused benchmark",
    "a verdict on what the evidence can predict when `rest` scores low: `rest` is an unfitted rule, so "
    "a low score says that rule transfers poorly",
)

#: the endpoint family of each assay: records of the held-out source's family at related intervals
#: are masked; records of another family are evidence in full (a reporter predicting a screen is the
#: question, a second screen of the same pair is a replicate)
FAMILY = {
    "crispri": "endogenous_regulation",
    "lentimpra": "reporter_activity",
    "vista": "in_vivo_activity",
    "satmut": "base_sensitivity",
    "gtex": "natural_variation_association",
}
#: what a unit of each assay is, which decides how a labelling is asked about it
SOURCE_KINDS = {
    "crispri": "pair",
    "lentimpra": "element",
    "vista": "element",
    "satmut": "base",
    "gtex": "pair",
}
#: each assay's endpoints, modelled separately. Binary endpoints map a unit's outcome to 1, 0 or
#: excluded (None); `activity` is continuous
DECREASE_EP, INCREASE_EP = "decrease", "increase"
ENDPOINTS = {
    "crispri": (DECREASE_EP, INCREASE_EP),
    "lentimpra": ("activity", "active"),
    "vista": ("positive",),
    "satmut": ("functional",),
    "gtex": ("associated",),
}
CONTINUOUS = {"activity"}
ENDPOINT_TEXT = {
    DECREASE_EP: "CRISPRi significant decrease of the measured gene (1) against a well-powered null (0); "
    "significant increases, underpowered nulls and missing effects excluded and counted (R2)",
    INCREASE_EP: "CRISPRi significant increase of the measured gene (1) against a well-powered null (0); "
    "significant decreases, underpowered nulls and missing effects excluded and counted (R2). An "
    "increase is an effect, not evidence of a silencer",
    "activity": "one lentiMPRA tile's log2(RNA/DNA) in that cell, per tile, never aggregated (R6)",
    "active": f"one lentiMPRA tile at or above {mpra.ACTIVE} log2(RNA/DNA) in that cell (1) or below (0)",
    "positive": "VISTA transgenic mouse e11.5: positive in any tissue (1) or negative (0)",
    "functional": "saturation mutagenesis: a measured base of a locus's primary experiment where some "
    "substitution is significant (1) or none is (0); repeats listed by measured.py, never pooled",
    "associated": "GTEx v8: a significant cis-eQTL variant inside the element for this protein-coding gene "
    "in any of 49 tissues (1), or a protein-coding gene with a TSS within CIS_WINDOW and no such variant "
    "(0: not significant and not tested are not told apart, weaker than a well-powered null)",
}
PRIMARY = {
    DECREASE_EP: "average_precision",
    INCREASE_EP: "average_precision",
    "activity": "spearman",
    "active": "average_precision",
    "positive": "average_precision",
    "functional": "average_precision",
    "associated": "average_precision",
}
SECONDARY = "auroc"  # for every endpoint; for `activity` it is the AUROC of `active`
METRIC_TEXT = {
    "average_precision": "average precision (area under the precision-recall step curve over distinct "
    "score thresholds, ties grouped), with the prevalence beside it",
    "spearman": "Spearman rank correlation (average ranks); under resampling the full-sample ranks are "
    "kept and weighted, not recomputed",
    "auroc": "area under the ROC curve (Mann-Whitney, ties counted one half)",
}
#: the interval: a percentile bootstrap over loci, the same resamples for every labelling of a source
RESAMPLES = 1000
SEED = 20260928
LEVEL = 0.95
LOCUS_BIN = 1_000_000
LOCUS_RULE = (
    "a locus is a connected component of a source's units, joined when their interval midpoints fall "
    f"in the same {LOCUS_BIN:,} bp bin of a chromosome or, for pair endpoints, when they name the same "
    "gene; each resample draws as many loci as the source has, with replacement"
)
#: floors below which an endpoint is described and not scored
MIN_POSITIVES = ms.MIN_FOR_A_COMPARISON
MIN_NEGATIVES = ms.MIN_FOR_A_COMPARISON
MIN_UNITS = ms.MIN_FOR_A_COMPARISON
MIN_LOCI = 10
#: abstention: a labelling returns None for a unit it says nothing about; None ranks below every
#: stated score, and coverage is reported beside the metric
ABSTENTION_RULE = "None ranks below every stated score (tied among themselves); coverage reported beside"
#: GTEx candidates: protein-coding genes with a TSS this close to the element's midpoint (GTEx's cis
#: window), plus every protein-coding gene GTEx ties to the element
CIS_WINDOW = 1_000_000
SPLIT_RULE = (
    "holding out S removes S and every provenance sibling of S (CRISPRi studies sharing a value of the "
    "benchmark's Reference column; the three lentiMPRA cells, one library ENCSR106SZM); masks every "
    "record of S's endpoint family whose interval shares at least one base with an S unit (the related "
    "intervals of measured.split_overlap); and never admits a pair from the CRISPRi benchmark's "
    "held-out file as evidence (R5)"
)
#: the labellings the first run scores for every source
REST, DISTANCE, UNCHANGED = "rest", "distance", "unchanged"
LABELLINGS = (REST, DISTANCE, UNCHANGED)
DISTANCE_RULE = (
    "minus the distance from the unit's midpoint to the TSS of the unit's gene (CRISPRi: the benchmark's "
    "own TSS columns; GTEx: GENCODE v50), or for element and base units to the nearest GENCODE v50 "
    "protein-coding TSS"
)
UNCHANGED_RULE = (
    "the compiled programs' predicted layer as it stands (data/knowledge/compiled/noncoding_*.bio, every "
    "`rule` with `evidence: predicted`), matched to a unit by at least one shared base. Pair endpoints: "
    "the largest strength among matched links naming the unit's gene in the endpoint's direction "
    "(decrease: activates; increase: inhibits; associated: either), 0 when the matched links name "
    "another gene or direction. Element and base endpoints: the largest strength of any matched link. "
    "None when no predicted element shares a base. The measured twins are never read"
)
REST_GENE_WEIGHT = 10
REST_RULE = (
    "unfitted: score = 10 x G + E + closeness, where G counts the sources in the view that support this "
    "gene at this interval (a CRISPRi training pair sharing a base, same gene, significant in the "
    "endpoint's direction; a GTEx association of an element sharing a base, same gene), E counts the "
    "sources that say the interval does something (lentiMPRA active by R6's rule in some cell over the "
    "tiles sharing a base; a VISTA positive sharing a base; a functional saturation-mutagenesis base "
    "inside; for non-regulation endpoints also a significant CRISPRi training pair on any gene and a "
    "GTEx association on any gene), and closeness = 1 / (1 + distance / 100 kb) with the distance of "
    "DISTANCE_RULE. G is at most 2 and E at most 5, so the order is G, then E, then distance"
)
CACHE_DIR = Path("data/cache/holdout")
CACHE_VERSION = "c4-1"
CHROMS = tuple(f"chr{i}" for i in range(1, 23)) + ("chrX", "chrY")
COMPILED_DIR = Path("data/knowledge/compiled")
GTEX_KNOWLEDGE = Path("data/knowledge/gtex")
#: a pseudo-source name: the CRISPRi benchmark's held-out file as evidence. Never allowed in `reads`
CRISPRI_HELDOUT_FILE = "crispri_heldout_file"
MODEL = "alphagenome"


class LeakError(ValueError):
    """A labelling read the held-out source, a sibling, a same-endpoint source outside the masked view,
    or the CRISPRi benchmark's test set."""


# --- units ----------------------------------------------------------------------------------------
@dataclass(frozen=True, slots=True)
class Unit:
    """One measured observation of one source: a pair, a tile or element, or a base."""

    source: str
    chrom: str
    start: int  # 0-based, half-open
    end: int
    outcome: str = ""
    value: float | None = None
    gene: str = ""  # pair endpoints: the gene's symbol
    gene_id: str = ""  # pair endpoints: Ensembl id without version
    tss: int | None = None  # pair endpoints: the gene's TSS, 0-based
    cell: str = ""
    split: str = ""  # CRISPRi: the benchmark partition of the pair
    study: str = ""

    @property
    def assay(self) -> str:
        return self.source.split(":", 1)[0]

    @property
    def mid(self) -> int:
        return (self.start + self.end) // 2


def assay_of(source: str) -> str:
    return source.split(":", 1)[0]


def target(unit: Unit, endpoint: str) -> float | None:
    """The endpoint's value for a unit: 1/0 for binary endpoints, the reading for `activity`, None when
    the unit is excluded from that endpoint."""
    a = unit.assay
    if a == "crispri":
        if unit.outcome == ms.NULL_INFORMATIVE:
            return 0.0
        if endpoint == DECREASE_EP:
            return 1.0 if unit.outcome == ms.DECREASE else None
        if endpoint == INCREASE_EP:
            return 1.0 if unit.outcome == ms.INCREASE else None
        raise ValueError(endpoint)
    if a == "lentimpra":
        if unit.value is None:
            return None
        return unit.value if endpoint == "activity" else float(unit.value >= mpra.ACTIVE)
    if a == "vista":
        return float(unit.outcome == "positive")
    if a == "satmut":
        return float(unit.outcome == "functional")
    if a == "gtex":
        return float(unit.outcome == "associated")
    raise ValueError(a)


# --- loaders (local caches only; nothing is fetched) ----------------------------------------------
def _crispri_extra(path: Path) -> list[dict[str, str]]:
    """The benchmark's valid rows in file order, for the columns measured.CrispriPair does not carry."""
    with gzip.open(path, "rt") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t") if r.get("ValidConnection") == "TRUE"]


def crispri_units(knowledge: Path = CRISPRI_KNOWLEDGE) -> dict[str, list[Unit]]:
    """One source per CRISPRi study (the benchmark's `Dataset`), both files, every valid pair."""
    out: dict[str, list[Unit]] = defaultdict(list)
    for name in ms.CRISPRI_FILES:
        p = knowledge / name
        if not p.exists():
            continue
        with gzip.open(p, "rt") as fh:
            pairs, _ = ms.parse_crispri(fh, None, ms.CRISPRI_SPLIT_OF[name], name)
        for pair, row in zip(pairs, _crispri_extra(p), strict=True):
            tss = row.get("startTSS")
            out[f"crispri:{pair.study}"].append(
                Unit(
                    source=f"crispri:{pair.study}",
                    chrom=pair.chrom,
                    start=pair.start,
                    end=pair.end,
                    outcome=pair.outcome,
                    value=None if pair.effect_size != pair.effect_size else pair.effect_size,
                    gene=pair.gene,
                    gene_id=(row.get("measuredGeneEnsemblId") or "").split(".")[0],
                    tss=int(tss) if tss not in (None, "", "NA") else None,
                    cell=pair.cell,
                    split=pair.split,
                    study=pair.study,
                )
            )
    return dict(out)


@lru_cache(maxsize=1)
def _references_cached() -> dict[str, frozenset[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    for name in ms.CRISPRI_FILES:
        p = CRISPRI_KNOWLEDGE / name
        if p.exists():
            for r in _crispri_extra(p):
                out[f"crispri:{r.get('Dataset', '')}"].add(r.get("Reference", ""))
    return {k: frozenset(v) for k, v in out.items()}


def crispri_references() -> dict[str, frozenset[str]]:
    """The benchmark's Reference values per study, which decide provenance siblings."""
    return _references_cached()


def lentimpra_units(knowledge: Path = mpra.KNOWLEDGE) -> dict[str, list[Unit]]:
    """One source per lentiMPRA cell; one unit per tile (ENCODE element interval) in that cell (R6)."""
    out: dict[str, list[Unit]] = {}
    for cell, acc in mpra.FILES.items():
        p = knowledge / f"{acc}.bed.gz"
        if not p.exists():
            continue
        into: dict[str, mpra.Element] = {}
        with gzip.open(p, "rt") as fh:
            mpra.parse(fh, cell, None, into)
        src = f"lentimpra:{cell}"
        out[src] = [
            Unit(
                source=src,
                chrom=e.chrom,
                start=e.start,
                end=e.end,
                outcome="active" if e.activity[cell] >= mpra.ACTIVE else "inactive",
                value=e.activity[cell],
                cell=cell,
                study=mpra.LIBRARY,
            )
            for e in into.values()
            if cell in e.activity
        ]
    return out


def vista_units(knowledge: Path = vista.KNOWLEDGE) -> dict[str, list[Unit]]:
    p = vista.locus_path(knowledge)
    if not p.exists():
        return {}
    with gzip.open(p, "rt") as fh:
        elements = vista.parse_loci(fh, None)
    return {
        "vista": [
            Unit(source="vista", chrom=e.chrom, start=e.start, end=e.end, outcome=e.status, study=e.id)
            for e in elements
        ]
    }


def satmut_units() -> dict[str, list[Unit]]:
    """One unit per measured base of each locus's primary experiment (repeats are never pooled)."""
    units = []
    for chrom in CHROMS:
        for e in ms.load_satmut(chrom):
            for pos, b in sorted(e.bases.items()):
                units.append(
                    Unit(
                        source="satmut",
                        chrom=chrom,
                        start=pos - 1,
                        end=pos,
                        outcome="functional" if b["functional"] else "inert",
                        study=e.experiment,
                    )
                )
    return {"satmut": units} if units else {}


@dataclass(frozen=True, slots=True)
class Gene:
    chrom: str
    tss: int  # 0-based
    gene_id: str  # without version
    symbol: str
    gene_type: str


def gencode_paths() -> list[Path]:
    from genomeos.genome.annotation import GENCODE_FULL

    if GENCODE_FULL.exists():
        return [GENCODE_FULL]
    return [p for c in CHROMS if (p := Path("data/reference") / f"gencode_v50_{c}.gff3.gz").exists()]


_ATTR = re.compile(r"(gene_id|gene_type|gene_name)=([^;]+)")


@lru_cache(maxsize=1)
def genes() -> tuple[Gene, ...]:
    """Every GENCODE v50 gene with its TSS, from the local reference (cached under CACHE_DIR)."""
    paths = gencode_paths()
    key = _files_key(paths)
    cached = _cache_read("genes", key)
    if cached is not None:
        return tuple(Gene(*g) for g in cached)
    out: list[Gene] = []
    for p in paths:
        with gzip.open(p, "rt") as fh:
            for line in fh:
                if "\tgene\t" not in line or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if f[2] != "gene":
                    continue
                attrs = dict(_ATTR.findall(f[8]))
                tss = int(f[3]) - 1 if f[6] == "+" else int(f[4]) - 1
                out.append(
                    Gene(
                        f[0],
                        tss,
                        attrs.get("gene_id", "").split(".")[0],
                        attrs.get("gene_name", ""),
                        attrs.get("gene_type", ""),
                    )
                )
    out.sort(key=lambda g: (g.chrom, g.tss))
    _cache_write("genes", key, [[g.chrom, g.tss, g.gene_id, g.symbol, g.gene_type] for g in out])
    return tuple(out)


@lru_cache(maxsize=1)
def coding_tss() -> dict[str, list[int]]:
    """Sorted protein-coding TSS positions per chromosome."""
    by: dict[str, list[int]] = defaultdict(list)
    for g in genes():
        if g.gene_type == "protein_coding":
            by[g.chrom].append(g.tss)
    return {c: sorted(v) for c, v in by.items()}


@lru_cache(maxsize=1)
def gene_by_id() -> dict[str, Gene]:
    return {g.gene_id: g for g in genes()}


def nearest_coding_distance(chrom: str, pos: int) -> int | None:
    starts = coding_tss().get(chrom)
    if not starts:
        return None
    i = bisect.bisect_left(starts, pos)
    return min(abs(starts[j] - pos) for j in (i - 1, i) if 0 <= j < len(starts))


def _gtex_hits(knowledge: Path) -> dict[tuple[str, int, int], set[str]]:
    hits: dict[tuple[str, int, int], set[str]] = defaultdict(set)
    for p in sorted(knowledge.glob("hits_*.tsv")):
        with p.open() as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                for eid in f[8].split(","):
                    chrom, _, span = eid.partition(":")
                    s, _, e = span.partition("-")
                    hits[(chrom, int(s), int(e))].add(f[5].split(".")[0])
    return hits


def gtex_units(knowledge: Path = GTEX_KNOWLEDGE) -> dict[str, list[Unit]]:
    """Element x protein-coding gene candidates at every element GTEx ties to a gene.

    The elements are the ones the eQTL distillation indexed (genomeos/attribution/eqtl.py keeps only
    variants inside elements an earlier target run predicted for), so the unit set was chosen by those
    runs; the registration says so."""
    hits = _gtex_hits(knowledge)
    if not hits:
        return {}
    by_id = gene_by_id()
    coding_by_chrom: dict[str, list[Gene]] = defaultdict(list)
    for g in genes():
        if g.gene_type == "protein_coding":
            coding_by_chrom[g.chrom].append(g)
    starts = {c: [g.tss for g in v] for c, v in coding_by_chrom.items()}
    units: list[Unit] = []
    for (chrom, s, e), gids in sorted(hits.items()):
        mid = (s + e) // 2
        cand = coding_by_chrom.get(chrom, [])
        lo = bisect.bisect_left(starts.get(chrom, []), mid - CIS_WINDOW)
        hi = bisect.bisect_right(starts.get(chrom, []), mid + CIS_WINDOW)
        chosen = {g.gene_id: g for g in cand[lo:hi]}
        for gid in gids:
            g = by_id.get(gid)
            if g is not None and g.gene_type == "protein_coding" and g.chrom == chrom:
                chosen.setdefault(gid, g)
        for gid, g in sorted(chosen.items()):
            units.append(
                Unit(
                    source="gtex",
                    chrom=chrom,
                    start=s,
                    end=e,
                    outcome="associated" if gid in gids else "not_associated",
                    gene=g.symbol,
                    gene_id=gid,
                    tss=g.tss,
                    study="GTEx v8",
                )
            )
    return {"gtex": units}


def gtex_excluded(knowledge: Path = GTEX_KNOWLEDGE) -> dict[str, int]:
    """Element-gene associations left out of the GTEx units: the gene is not protein-coding in GENCODE
    v50, or not in it at all. Counted, never scored."""
    by_id = gene_by_id()
    pairs = {(k, gid) for k, gids in _gtex_hits(knowledge).items() for gid in gids}
    out = {"associations": len(pairs), "gene_not_in_gencode_v50": 0, "gene_not_protein_coding": 0}
    for _, gid in pairs:
        g = by_id.get(gid)
        if g is None:
            out["gene_not_in_gencode_v50"] += 1
        elif g.gene_type != "protein_coding":
            out["gene_not_protein_coding"] += 1
    return out


@lru_cache(maxsize=1)
def all_units() -> dict[str, tuple[Unit, ...]]:
    """Every holdable source, by name, from the local caches."""
    out: dict[str, list[Unit]] = {}
    for part in (crispri_units(), lentimpra_units(), vista_units(), satmut_units(), gtex_units()):
        out.update(part)
    return {k: tuple(sorted(v, key=_unit_order)) for k, v in out.items()}


def _unit_order(u: Unit) -> tuple:
    return (u.chrom, u.start, u.end, u.gene, u.cell, u.outcome)


def siblings(source: str, references: dict[str, frozenset[str]] | None = None) -> tuple[str, ...]:
    """Provenance siblings under SPLIT_RULE (the source itself excluded)."""
    a = assay_of(source)
    if a == "lentimpra":
        return tuple(f"lentimpra:{c}" for c in mpra.FILES if f"lentimpra:{c}" != source)
    if a == "crispri":
        refs = references if references is not None else crispri_references()
        mine = refs.get(source, frozenset())
        return tuple(sorted(s for s, r in refs.items() if s != source and r & mine))
    return ()


# --- loci -----------------------------------------------------------------------------------------
def loci(units: Iterable[Unit]) -> list[int]:
    """Locus codes 0..L-1 under LOCUS_RULE, in the order of `units`."""
    parent: dict[Any, Any] = {}

    def find(x: Any) -> Any:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    keys = []
    for u in units:
        b = ("bin", u.chrom, u.mid // LOCUS_BIN)
        rb = find(b)
        if SOURCE_KINDS[u.assay] == "pair" and (u.gene_id or u.gene):
            rg = find(("gene", u.gene_id or u.gene))
            if rg != rb:
                parent[rg] = rb
        keys.append(b)
    codes: dict[Any, int] = {}
    return [codes.setdefault(find(k), len(codes)) for k in keys]


# --- the evidence view: everything but the held-out source, masked by locus -----------------------
class _Index:
    """Units (or anything with chrom, start, end) searchable by shared base."""

    def __init__(self, items: Iterable[Any]) -> None:
        by: dict[str, list[Any]] = defaultdict(list)
        for u in items:
            by[u.chrom].append(u)
        self.by = {c: sorted(v, key=lambda u: u.start) for c, v in by.items()}
        self.starts = {c: [u.start for u in v] for c, v in self.by.items()}
        self.reach = max((u.end - u.start for v in self.by.values() for u in v), default=0)

    def near(self, chrom: str, start: int, end: int) -> list[Any]:
        items = self.by.get(chrom)
        if not items:
            return []
        st = self.starts[chrom]
        lo = bisect.bisect_left(st, start - self.reach)
        hi = bisect.bisect_left(st, end)
        return [u for u in items[lo:hi] if u.end > start and u.start < end]

    def __len__(self) -> int:
        return sum(len(v) for v in self.by.values())


@dataclass
class EvidenceView:
    """What a labelling may read while `without` is held out (SPLIT_RULE)."""

    without: str
    sources: dict[str, tuple[Unit, ...]]
    masked: dict[str, int] = field(default_factory=dict)
    removed: tuple[str, ...] = ()
    heldout_file_pairs_dropped: int = 0

    def of_assay(self, assay: str) -> list[str]:
        return sorted(s for s in self.sources if assay_of(s) == assay)

    def units_of(self, assay: str) -> list[Unit]:
        return [u for s in self.of_assay(assay) for u in self.sources[s]]


def evidence(
    without: str,
    units: dict[str, tuple[Unit, ...]] | None = None,
    references: dict[str, frozenset[str]] | None = None,
) -> EvidenceView:
    """The view for holding out `without`: the source and its siblings removed, same-family records at
    related intervals masked, and no pair of the CRISPRi held-out file anywhere."""
    units = units if units is not None else all_units()
    if without not in units:
        raise KeyError(f"unknown source {without!r}")
    if references is None and assay_of(without) == "crispri":
        references = crispri_references()
    gone = {without, *siblings(without, references)}
    held = _Index(units[without])
    fam = FAMILY[assay_of(without)]
    out: dict[str, tuple[Unit, ...]] = {}
    masked: dict[str, int] = {}
    dropped = 0
    for name, us in units.items():
        if name in gone:
            continue
        keep = []
        n_mask = 0
        for u in us:
            if u.assay == "crispri" and u.split != ms.TRAINING:
                dropped += 1
                continue  # the benchmark's test set is never evidence (R5)
            if FAMILY[u.assay] == fam and held.near(u.chrom, u.start, u.end):
                n_mask += 1
                continue
            keep.append(u)
        if n_mask:
            masked[name] = n_mask
        if keep:
            out[name] = tuple(keep)
    return EvidenceView(
        without=without,
        sources=out,
        masked=masked,
        removed=tuple(sorted(gone)),
        heldout_file_pairs_dropped=dropped,
    )


def _as_pair(u: Unit, split: str) -> ms.CrispriPair:
    return ms.CrispriPair(
        chrom=u.chrom,
        start=u.start,
        end=u.end,
        gene=u.gene,
        cell=u.cell,
        dataset=u.study,
        reference="",
        regulated=u.outcome == ms.DECREASE,
        significant=u.outcome in (ms.DECREASE, ms.INCREASE),
        effect_size=u.value if u.value is not None else float("nan"),
        p_adjusted=1.0,
        split=split,
    )


def relation_counts(held_out: str, units: dict[str, tuple[Unit, ...]] | None = None) -> dict[str, Any]:
    """How the held-out source's units relate to the rest of its endpoint family before masking.

    CRISPRi sources reuse lane-split's `measured.split_overlap` itself: the held-out study's pairs play
    the held-out split and every other study's pairs the training split. Other sources are counted in
    the same strictest-first order, without the gene categories."""
    units = units if units is not None else all_units()
    fam = FAMILY[assay_of(held_out)]
    others = [u for n, us in units.items() if n != held_out and FAMILY[assay_of(n)] == fam for u in us]
    if assay_of(held_out) == "crispri":
        pairs = [_as_pair(u, ms.HELDOUT) for u in units[held_out]] + [
            _as_pair(u, ms.TRAINING) for u in others
        ]
        r = ms.split_overlap(pairs)
        return {
            "reused": "measured.split_overlap",
            "counts": r["counts"],
            "related": r["heldout_related_to_training"],
        }
    idx = _Index(others)
    counts = dict.fromkeys(("identical", "near_identical", "overlapping", "independent"), 0)
    for u in units[held_out]:
        near = idx.near(u.chrom, u.start, u.end)
        if any(v.start == u.start and v.end == u.end for v in near):
            counts["identical"] += 1
        elif any(ms.reciprocal_overlap(u.start, u.end, v.start, v.end) >= ms.NEAR_IDENTICAL for v in near):
            counts["near_identical"] += 1
        elif near:
            counts["overlapping"] += 1
        else:
            counts["independent"] += 1
    return {
        "reused": "measured.split_overlap's order, without the gene categories",
        "counts": counts,
        "related": len(units[held_out]) - counts["independent"],
    }


# --- labellings -----------------------------------------------------------------------------------
@dataclass(frozen=True)
class Labels:
    """A labelling as the harness scores it.

    `predict(unit, endpoint)` returns a score (higher = more of the endpoint) or None to abstain.
    `reads` names every source and input the labelling read (assay sources by their source name,
    `MODEL` for AlphaGenome output, `CRISPRI_HELDOUT_FILE` for the benchmark's test set). A labelling
    built from `evidence(without=S)` sets `built_without=S`, the only way it may read S's family."""

    name: str
    predict: Callable[[Unit, str], float | None]
    reads: frozenset[str] = frozenset()
    built_without: str | None = None
    note: str = ""

    @property
    def model_dependent(self) -> bool:
        return MODEL in self.reads


def check_provenance(
    labels: Labels, held_out: str, references: dict[str, frozenset[str]] | None = None
) -> None:
    """Raise LeakError when `labels` read what holding out `held_out` forbids (SPLIT_RULE)."""
    if references is None and assay_of(held_out) == "crispri":
        references = crispri_references()
    gone = {held_out, *siblings(held_out, references)}
    bad = sorted(gone & labels.reads)
    if bad:
        raise LeakError(f"{labels.name} read {bad} while {held_out} is held out")
    if CRISPRI_HELDOUT_FILE in labels.reads:
        raise LeakError(f"{labels.name} read the CRISPRi benchmark's held-out file (R5)")
    fam = FAMILY[assay_of(held_out)]
    same = sorted(r for r in labels.reads if assay_of(r) in FAMILY and FAMILY[assay_of(r)] == fam)
    if same and labels.built_without != held_out:
        raise LeakError(
            f"{labels.name} read {same}, of {held_out}'s endpoint family, outside "
            f"evidence(without={held_out!r})"
        )


def distance_of(unit: Unit) -> int | None:
    if SOURCE_KINDS[unit.assay] == "pair":
        return abs(unit.mid - unit.tss) if unit.tss is not None else None
    return nearest_coding_distance(unit.chrom, unit.mid)


def distance_labels() -> Labels:
    def predict(u: Unit, endpoint: str) -> float | None:
        d = distance_of(u)
        return None if d is None else -float(d)

    return Labels(DISTANCE, predict, frozenset({"gencode_v50"}), note=DISTANCE_RULE)


@dataclass(frozen=True, slots=True)
class Link:
    """One predicted link of a compiled program: element, gene, action, strength, cell."""

    chrom: str
    start: int
    end: int
    element: str
    gene: str
    action: str  # activates | inhibits
    strength: float
    cell: str


_RULE = re.compile(
    r"^rule (\S+) (activates|inhibits) (\S+) \{ strength: ([0-9.eE+-]+); when: cell_type = ([^;]+); "
    r"evidence: (\w+)"
)
LOCUS_RE = re.compile(r"^(chr[0-9XYM]+):(\d+)-(\d+)$")


def program_blocks(path: Path) -> Iterator[tuple[str, str, dict[str, str]]]:
    """(kind, name, fields) for every `element` and `region` block and every `rule` line of a compiled
    program, streamed. A rule's fields hold action, gene, strength, cell, evidence and `heldout`."""
    kind: str | None = None
    name = ""
    fields: dict[str, str] = {}
    with path.open() as fh:
        for line in fh:
            if kind is not None:
                if line.startswith("}"):
                    yield kind, name, fields
                    kind = None
                    continue
                k, sep, v = line.strip().partition(": ")
                if sep:
                    fields[k] = v
                continue
            if line.startswith("rule "):
                m = _RULE.match(line)
                if m:
                    yield (
                        "rule",
                        m.group(1),
                        {
                            "action": m.group(2),
                            "gene": m.group(3),
                            "strength": m.group(4),
                            "cell": m.group(5).strip(),
                            "evidence": m.group(6),
                            "heldout": "1" if ms.HELDOUT_MARK in line else "",
                        },
                    )
                continue
            if line.startswith(("element ", "region ")) and line.rstrip().endswith("{"):
                parts = line.split()
                kind, name = parts[0], parts[1]
                fields = {}


def compiled_programs(compiled: Path = COMPILED_DIR) -> list[Path]:
    return [p for c in CHROMS if (p := compiled / f"noncoding_{c}.bio").exists()]


def compiled_links(compiled: Path = COMPILED_DIR) -> dict[str, list[Link]]:
    """Every predicted link of the compiled programs, per chromosome, sorted (cached under CACHE_DIR)."""
    paths = compiled_programs(compiled)
    key = _files_key(paths)
    cached = _cache_read("compiled_links", key)
    if cached is not None:
        return {c: [Link(*x) for x in v] for c, v in cached.items()}
    out: dict[str, list[Link]] = defaultdict(list)
    for p in paths:
        loci_of: dict[str, tuple[str, int, int]] = {}
        for kind, name, f in program_blocks(p):
            if kind == "element" and not name.endswith("_measured"):
                m = LOCUS_RE.match(f.get("locus", ""))
                if m:
                    loci_of[name] = (m.group(1), int(m.group(2)), int(m.group(3)))
            elif kind == "rule" and f["evidence"] == "predicted" and name in loci_of:
                c, s, e = loci_of[name]
                out[c].append(Link(c, s, e, name, f["gene"], f["action"], float(f["strength"]), f["cell"]))
    res = {c: sorted(v, key=lambda x: (x.start, x.end, x.gene)) for c, v in out.items()}
    _cache_write(
        "compiled_links",
        key,
        {
            c: [[x.chrom, x.start, x.end, x.element, x.gene, x.action, x.strength, x.cell] for x in v]
            for c, v in res.items()
        },
    )
    return res


def unchanged_labels(links: dict[str, list[Link]] | None = None) -> Labels:
    idx = _Index(x for v in (links if links is not None else compiled_links()).values() for x in v)

    def predict(u: Unit, endpoint: str) -> float | None:
        near = idx.near(u.chrom, u.start, u.end)
        if not near:
            return None
        if SOURCE_KINDS[u.assay] != "pair":
            return max(x.strength for x in near)
        want = {DECREASE_EP: ("activates",), INCREASE_EP: ("inhibits",)}.get(
            endpoint, ("activates", "inhibits")
        )
        hit = [x.strength for x in near if x.gene == u.gene and x.action in want]
        return max(hit) if hit else 0.0

    return Labels(UNCHANGED, predict, frozenset({MODEL, "encode_ccre_registry"}), note=UNCHANGED_RULE)


def rest_labels(view: EvidenceView) -> Labels:
    """The prediction of the held-out endpoint from the remaining sources (REST_RULE)."""
    crispri_idx = _Index(view.units_of("crispri"))
    gtex_idx = _Index(u for u in view.sources.get("gtex", ()) if u.outcome == "associated")
    mpra_idx = [_Index(view.sources[s]) for s in view.of_assay("lentimpra")]
    vista_idx = _Index(view.sources.get("vista", ()))
    satmut_idx = _Index(view.sources.get("satmut", ()))

    def predict(u: Unit, endpoint: str) -> float | None:
        s, e = u.start, u.end
        g = el = 0
        if SOURCE_KINDS[u.assay] == "pair":
            want = {DECREASE_EP: (ms.DECREASE,), INCREASE_EP: (ms.INCREASE,)}.get(
                endpoint, (ms.DECREASE, ms.INCREASE)
            )
            if any(p.gene == u.gene and p.outcome in want for p in crispri_idx.near(u.chrom, s, e)):
                g += 1
            if u.gene_id and any(h.gene_id == u.gene_id for h in gtex_idx.near(u.chrom, s, e)):
                g += 1
        else:
            if any(p.outcome in (ms.DECREASE, ms.INCREASE) for p in crispri_idx.near(u.chrom, s, e)):
                el += 1
            if gtex_idx.near(u.chrom, s, e):
                el += 1
        for idx in mpra_idx:
            tiles = [t.value for t in idx.near(u.chrom, s, e) if t.value is not None]
            if tiles and ms.reporter_label(tiles) == ms.LABEL_ACTIVE:
                el += 1
                break
        if any(v.outcome == "positive" for v in vista_idx.near(u.chrom, s, e)):
            el += 1
        if any(b.outcome == "functional" for b in satmut_idx.near(u.chrom, s, e)):
            el += 1
        d = distance_of(u)
        closeness = 0.0 if d is None else 1.0 / (1.0 + d / 100_000)
        return REST_GENE_WEIGHT * g + el + closeness

    return Labels(
        REST, predict, frozenset(view.sources) | {"gencode_v50"}, built_without=view.without, note=REST_RULE
    )


# --- metrics --------------------------------------------------------------------------------------
def _order(scores: list[float | None]) -> tuple[Any, Any, Any]:
    """Scores as floats with abstentions at -inf, the descending order, and each score group's end."""
    import numpy as np

    x = np.array([-math.inf if v is None else float(v) for v in scores], dtype=float)
    order = np.argsort(-x, kind="stable")
    xs = x[order]
    ends = np.flatnonzero(np.r_[xs[1:] != xs[:-1], True])
    return x, order, ends


def binary_metrics(scores: list[float | None], y: list[float], weights: Any) -> dict[str, Any]:
    """Weighted average precision and AUROC for each row of `weights` (R x N)."""
    import numpy as np

    _, order, ends = _order(scores)
    yv = np.asarray(y, dtype=float)[order]
    w = np.atleast_2d(np.asarray(weights, dtype=float))[:, order]
    starts = np.r_[0, ends[:-1] + 1]
    pos_g = np.add.reduceat(w * yv, starts, axis=1)
    neg_g = np.add.reduceat(w * (1.0 - yv), starts, axis=1)
    tp = np.cumsum(pos_g, axis=1)
    fp = np.cumsum(neg_g, axis=1)
    n_pos = tp[:, -1]
    n_neg = fp[:, -1]
    with np.errstate(invalid="ignore", divide="ignore"):
        precision = np.where(tp + fp > 0, tp / np.where(tp + fp > 0, tp + fp, 1.0), 0.0)
        ap = np.sum(pos_g * precision, axis=1) / n_pos
        auc = np.sum(pos_g * ((n_neg[:, None] - fp) + 0.5 * neg_g), axis=1) / (n_pos * n_neg)
    bad = (n_pos <= 0) | (n_neg <= 0)
    ap = np.where(bad, np.nan, ap)
    auc = np.where(bad, np.nan, auc)
    return {"average_precision": ap, "auroc": auc}


def _ranks(x: Any) -> Any:
    """Average ranks (1-based), ties shared."""
    import numpy as np

    order = np.argsort(x, kind="stable")
    xs = x[order]
    n = len(x)
    ranks = np.empty(n, dtype=float)
    ends = np.flatnonzero(np.r_[xs[1:] != xs[:-1], True]) if n else np.array([], dtype=int)
    start = 0
    for end in ends:
        ranks[order[start : end + 1]] = (start + end) / 2.0 + 1.0
        start = end + 1
    return ranks


def spearman(scores: list[float | None], values: list[float], weights: Any) -> Any:
    """Weighted Pearson correlation of the full-sample average ranks, per row of `weights`."""
    import numpy as np

    x, _, _ = _order(scores)
    rx = _ranks(x)
    ry = _ranks(np.asarray(values, dtype=float))
    w = np.atleast_2d(np.asarray(weights, dtype=float))
    sw = w.sum(axis=1)
    mx = (w @ rx) / sw
    my = (w @ ry) / sw
    dx = rx[None, :] - mx[:, None]
    dy = ry[None, :] - my[:, None]
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.sum(w * dx * dy, axis=1) / np.sqrt(
            np.sum(w * dx * dx, axis=1) * np.sum(w * dy * dy, axis=1)
        )


def metrics(endpoint: str, scores: list[float | None], y: list[float], weights: Any) -> dict[str, Any]:
    if endpoint in CONTINUOUS:
        b = binary_metrics(scores, [float(v >= mpra.ACTIVE) for v in y], weights)
        return {
            "spearman": spearman(scores, y, weights),
            "auroc": b["auroc"],
            "average_precision": b["average_precision"],
        }
    return binary_metrics(scores, y, weights)


def resample_counts(n_loci: int, resamples: int = RESAMPLES, seed: int = SEED) -> Any:
    """R x L draw counts: each resample draws as many loci as there are, with replacement. A unit's
    weight in a resample is its locus's count; the draws depend only on the number of loci, so every
    labelling of a source is scored on the same resamples."""
    import numpy as np

    rng = np.random.default_rng(seed)
    counts = np.zeros((resamples, n_loci), dtype=float)
    for r in range(resamples):
        counts[r] = np.bincount(rng.integers(0, n_loci, size=n_loci), minlength=n_loci)
    return counts


# --- the call the pilot makes ---------------------------------------------------------------------
_MEMO: dict[str, dict[str, Any]] = {}
BLOCK_CELLS = 4_000_000  # resamples x units per numpy block, so memory stays bounded on the largest source


def score(
    labels: Labels,
    held_out_source: str,
    endpoint: str | None = None,
    *,
    units: dict[str, tuple[Unit, ...]] | None = None,
    chroms: Iterable[str] | None = None,
    resamples: int = RESAMPLES,
    seed: int = SEED,
    cache: bool = True,
    references: dict[str, frozenset[str]] | None = None,
) -> dict[str, Any]:
    """The labelling's prediction of the held-out source's endpoint, with a 95% locus-bootstrap interval.

    `endpoint` defaults to the source's first (primary) endpoint. `chroms` restricts evaluation to a
    chromosome partition. Deterministic for fixed `resamples` and `seed`; cached by the predictions
    themselves, so two labellings that predict the same numbers share an entry."""
    return evaluate(
        [labels],
        held_out_source,
        endpoint,
        units=units,
        chroms=chroms,
        resamples=resamples,
        seed=seed,
        cache=cache,
        references=references,
    )[0]


def difference(ra: dict[str, Any], rb: dict[str, Any]) -> dict[str, Any]:
    """`ra` minus `rb` from two results of one `evaluate(..., keep_draws=True)` call: the same units and
    the same resamples, so the interval is paired."""
    import numpy as np

    out = {"a": ra["labels"], "b": rb["labels"], "source": ra["source"], "endpoint": ra["endpoint"]}
    if ra["status"] != "scored" or rb["status"] != "scored":
        return {**out, "status": "described_not_scored", "why": ra.get("why") or rb.get("why")}
    d = np.asarray(ra["_draws"], dtype=float) - np.asarray(rb["_draws"], dtype=float)
    d = d[~np.isnan(d)]
    q = [(1 - LEVEL) / 2, 1 - (1 - LEVEL) / 2]
    return {
        **out,
        "status": "scored",
        "metric": ra["metric"],
        "difference": round(ra["value"] - rb["value"], 4),
        "interval": [round(float(x), 4) for x in np.quantile(d, q)] if len(d) else None,
        "level": LEVEL,
        "a_value": ra["value"],
        "b_value": rb["value"],
    }


def compare(
    a: Labels, b: Labels, held_out_source: str, endpoint: str | None = None, **kw: Any
) -> dict[str, Any]:
    """`a` minus `b` on the same held-out units and the same resamples, with a paired interval."""
    ra, rb = evaluate([a, b], held_out_source, endpoint, keep_draws=True, **kw)
    return difference(ra, rb)


def evaluate(
    labellings: list[Labels],
    held_out_source: str,
    endpoint: str | None = None,
    *,
    units: dict[str, tuple[Unit, ...]] | None = None,
    chroms: Iterable[str] | None = None,
    resamples: int = RESAMPLES,
    seed: int = SEED,
    cache: bool = True,
    references: dict[str, frozenset[str]] | None = None,
    keep_draws: bool = False,
) -> list[dict[str, Any]]:
    """`score` for several labellings of one held-out source, on the same units and resamples."""
    import numpy as np

    units = units if units is not None else all_units()
    if held_out_source not in units:
        raise KeyError(f"unknown source {held_out_source!r}")
    assay = assay_of(held_out_source)
    endpoint = endpoint or ENDPOINTS[assay][0]
    if endpoint not in ENDPOINTS[assay]:
        raise ValueError(f"{held_out_source} has endpoints {ENDPOINTS[assay]}, not {endpoint!r}")
    for lab in labellings:
        check_provenance(lab, held_out_source, references)
    only = set(chroms) if chroms is not None else None
    excluded: dict[str, int] = defaultdict(int)
    kept: list[Unit] = []
    ys: list[float] = []
    for u in units[held_out_source]:
        if only is not None and u.chrom not in only:
            continue
        t = target(u, endpoint)
        if t is None:
            excluded[u.outcome or "no_value"] += 1
            continue
        kept.append(u)
        ys.append(t)
    codes = loci(kept)
    n_loci = max(codes) + 1 if codes else 0
    base: dict[str, Any] = {
        "source": held_out_source,
        "endpoint": endpoint,
        "metric": PRIMARY[endpoint],
        "units": len(kept),
        "loci": n_loci,
        "excluded": dict(sorted(excluded.items())),
        "chromosomes": sorted(only) if only is not None else "all",
        "benchmark_status": STATUS,
    }
    if endpoint in CONTINUOUS:
        base["active_units"] = int(sum(1 for v in ys if v >= mpra.ACTIVE))
        short = len(kept) < MIN_UNITS
        why = f"{len(kept)} units, below the floor of {MIN_UNITS}" if short else ""
    else:
        pos = int(sum(ys))
        base["positives"], base["negatives"] = pos, len(ys) - pos
        base["prevalence"] = round(pos / len(ys), 4) if ys else None
        short = pos < MIN_POSITIVES or len(ys) - pos < MIN_NEGATIVES
        why = (
            f"{pos} positives and {len(ys) - pos} negatives; the floor is {MIN_POSITIVES} of each"
            if short
            else ""
        )
    if not short and n_loci < MIN_LOCI:
        short, why = True, f"{n_loci} loci, below the floor of {MIN_LOCI}"
    results = []
    counts = code_arr = None
    prim = PRIMARY[endpoint]
    q = [(1 - LEVEL) / 2, 1 - (1 - LEVEL) / 2]
    for lab in labellings:
        preds = [lab.predict(u, endpoint) for u in kept]
        stated = sum(1 for p in preds if p is not None)
        r = {
            **base,
            "labels": lab.name,
            "reads": sorted(lab.reads),
            "model_dependent": lab.model_dependent,
            "coverage": {
                "stated": stated,
                "abstained": len(kept) - stated,
                "share": round(stated / len(kept), 4) if kept else None,
            },
        }
        if short:
            results.append({**r, "status": "described_not_scored", "why": why})
            continue
        key = _score_key(held_out_source, endpoint, kept, ys, preds, resamples, seed, only)
        hit = _MEMO.get(key) or (_cache_read("score", key) if cache else None)
        if hit is not None and not (keep_draws and "_draws" not in hit):
            _MEMO[key] = hit
            out = {**r, **{k: v for k, v in hit.items() if k != "_draws"}}
            if keep_draws:
                out["_draws"] = hit["_draws"]
            results.append(out)
            continue
        if counts is None:
            counts = resample_counts(n_loci, resamples, seed)
            code_arr = np.asarray(codes, dtype=int)
        full = metrics(endpoint, preds, ys, np.ones((1, len(kept))))
        prim_draws, sec_draws = [], []
        step = max(1, BLOCK_CELLS // max(1, len(kept)))
        for i in range(0, resamples, step):
            m = metrics(endpoint, preds, ys, counts[i : i + step][:, code_arr])
            prim_draws.append(m[prim])
            sec_draws.append(m[SECONDARY])
        vals = np.concatenate(prim_draws)
        sec = np.concatenate(sec_draws)
        ok, sok = vals[~np.isnan(vals)], sec[~np.isnan(sec)]
        body: dict[str, Any] = {
            "status": "scored",
            "value": round(float(full[prim][0]), 4),
            "interval": [round(float(x), 4) for x in np.quantile(ok, q)] if len(ok) else None,
            "level": LEVEL,
            "secondary": {
                "metric": "auroc_of_active" if endpoint in CONTINUOUS else SECONDARY,
                "value": round(float(full[SECONDARY][0]), 4),
                "interval": [round(float(x), 4) for x in np.quantile(sok, q)] if len(sok) else None,
            },
            "bootstrap": {
                "unit": "locus",
                "resamples": resamples,
                "seed": seed,
                "failed_resamples": int(np.isnan(vals).sum()),
            },
            "_draws": [None if math.isnan(v) else round(float(v), 6) for v in vals],
        }
        _MEMO[key] = body
        if cache:
            _cache_write("score", key, body)
        out = {**r, **{k: v for k, v in body.items() if k != "_draws"}}
        if keep_draws:
            out["_draws"] = body["_draws"]
        results.append(out)
    for out in results:
        if "_draws" in out:
            out["_draws"] = [math.nan if v is None else v for v in out["_draws"]]
    return results


def _score_key(
    source: str,
    endpoint: str,
    kept: list[Unit],
    ys: list[float],
    preds: list,
    resamples: int,
    seed: int,
    only: set[str] | None,
) -> str:
    h = hashlib.sha256()
    h.update(
        json.dumps(
            [CACHE_VERSION, source, endpoint, resamples, seed, sorted(only) if only else None]
        ).encode()
    )
    for u, y, p in zip(kept, ys, preds, strict=True):
        h.update(f"{u.chrom}:{u.start}-{u.end}|{u.gene}|{u.cell}|{y}|{p!r}\n".encode())
    return h.hexdigest()


# --- the small disk cache -------------------------------------------------------------------------
def _files_key(paths: list[Path]) -> str:
    """A cache key from each file's name, size and modification time. It only decides whether a parse
    can be reused; the result manifest records each input's sha256."""
    h = hashlib.sha256(CACHE_VERSION.encode())
    for p in paths:
        st = p.stat()
        h.update(f"{p}|{st.st_size}|{st.st_mtime_ns}\n".encode())
    return h.hexdigest()


def _cache_path(kind: str, key: str) -> Path:
    return CACHE_DIR / kind / f"{key}.json.gz"


def _cache_read(kind: str, key: str) -> Any:
    p = _cache_path(kind, key)
    if not p.exists():
        return None
    try:
        with gzip.open(p, "rt") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def _cache_write(kind: str, key: str, value: Any) -> None:
    p = _cache_path(kind, key)
    try:
        p.parent.mkdir(parents=True, exist_ok=True)
        tmp = p.with_name(p.name + ".tmp")
        with gzip.open(tmp, "wt") as fh:
            json.dump(value, fh)
        tmp.replace(p)
    except OSError:
        pass


def sources(units: dict[str, tuple[Unit, ...]] | None = None) -> list[str]:
    """Every holdable source present, in a fixed order."""
    order = {a: i for i, a in enumerate(FAMILY)}
    return sorted(units if units is not None else all_units(), key=lambda s: (order[assay_of(s)], s))


def registration() -> dict[str, Any]:
    """The registered constants, as a result file states them."""
    return {
        "registered": REGISTERED,
        "status": STATUS,
        "reuse": REUSE,
        "may_be_called": MAY_BE_CALLED,
        "may_not_be_called": list(MAY_NOT_BE_CALLED),
        "families": FAMILY,
        "endpoints": {k: list(v) for k, v in ENDPOINTS.items()},
        "endpoint_text": ENDPOINT_TEXT,
        "primary_metric": PRIMARY,
        "secondary_metric": SECONDARY,
        "metric_text": METRIC_TEXT,
        "split_rule": SPLIT_RULE,
        "locus_rule": LOCUS_RULE,
        "resamples": RESAMPLES,
        "seed": SEED,
        "level": LEVEL,
        "floors": {
            "positives": MIN_POSITIVES,
            "negatives": MIN_NEGATIVES,
            "units": MIN_UNITS,
            "loci": MIN_LOCI,
        },
        "abstention": ABSTENTION_RULE,
        "labellings": {REST: REST_RULE, DISTANCE: DISTANCE_RULE, UNCHANGED: UNCHANGED_RULE},
        "cis_window": CIS_WINDOW,
    }
