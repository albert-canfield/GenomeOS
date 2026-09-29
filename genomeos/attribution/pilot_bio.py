# SPDX-License-Identifier: AGPL-3.0-or-later
"""The coherence pilot on real evidence: gate 2 (item 12 S3; item 13 C1-C3 and the metric).

The debugger of `genomeos/attribution/pilot.py`, unchanged, is run on ENCODE cCREs of the validation
chromosomes with the evidence C4's harness allows when one source is held out
(`holdout.evidence(without=S)`), and its revised labels are scored on that held-out source with
`holdout.evaluate` beside three baselines: the unchanged labels, the independent-block variant and
distance to TSS. The registration (chromosomes, endpoints, pass rule, falsifier, readings, metric and
compute budget) is in this module's constants and in docs/ATTRIBUTION.md, written before any score of
any pilot labelling was read. An internal development result on withheld sources, never a validation.

Only the query's coordinates, gene, TSS and cell are read from a held-out unit; its outcome and effect
are never read by a labelling. The held-out source's own outcomes are read afterwards, by the harness
and by `validate`, which counts how many committed corrections improved the prediction of the units
they touch.
"""

from __future__ import annotations

import bisect
import gzip
import hashlib
import math
import time
from collections import Counter, defaultdict
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from genomeos.attribution import holdout as ho
from genomeos.attribution import measured as ms
from genomeos.attribution import pilot as pl

# ==================================================================================================
# Gate 2's registration (2026-09-29, lane-pilot), fixed before any pilot labelling was scored.
# ==================================================================================================
REGISTERED = "2026-09-29"
STATUS = pl.STATUS
#: the validation partition: the chromosomes the pilot is run and scored on. Chosen from unit counts
#: alone (no score of any labelling was read on any chromosome), as the smallest set found that keeps
#: all four CRISPRi decrease endpoints above the harness's floors; the other 17 are the development
#: partition, used only to time and debug the pilot, never scored
VALIDATION = ("chr1", "chr6", "chr8", "chr9", "chr10", "chr11", "chr19")
DEVELOPMENT = tuple(c for c in ho.CHROMS if c not in VALIDATION)
CHROMOSOME_CHOICE = (
    "chromosomes are the validation partitions. The seven are the smallest set found by unit counts "
    "alone that keeps all four CRISPRi decrease endpoints above the harness's floors (20 positives, 20 "
    "negatives, 10 loci): Gasperini2019 150/2,637 in 201 loci, Morris 21/71 in 33, Schraivogel2020 "
    "23/1,276 in 12 (chr8 and chr11 hold all of its units), Xie 26/148 in 22. Odd or even chromosomes "
    "lose Schraivogel2020 and one of Morris or Xie. Of the 15 endpoints C4 scored, 14 stay scorable; "
    "WTC11_DC_TAP increase (20 positives genome-wide) cannot be scored on any subset. The other 17 "
    "chromosomes are the development partition: the pilot was timed and debugged there without any "
    "held-out score being read, and nothing is scored there"
)
PRIMARY = (
    ("crispri:Gasperini2019", "decrease"),
    ("crispri:Morris", "decrease"),
    ("crispri:Schraivogel2020", "decrease"),
    ("crispri:Xie", "decrease"),
)
SECONDARY = (
    ("crispri:Gasperini2019", "increase"),
    ("lentimpra:K562", "activity"),
    ("lentimpra:K562", "active"),
    ("lentimpra:HepG2", "activity"),
    ("lentimpra:HepG2", "active"),
    ("lentimpra:WTC11", "activity"),
    ("lentimpra:WTC11", "active"),
    ("vista", "positive"),
    ("satmut", "functional"),
    ("gtex", "associated"),
)
HELD_OUT = tuple(dict.fromkeys(s for s, _ in (*PRIMARY, *SECONDARY)))
#: the binary endpoint each held-out source validates corrections against
VALIDATION_ENDPOINT = {
    "crispri": ho.DECREASE_EP,
    "lentimpra": "active",
    "vista": "positive",
    "satmut": "functional",
    "gtex": "associated",
}
PILOT, INDEPENDENT, PRIOR, PILOT_COVERAGE = (
    "pilot",
    "independent_block",
    "prior_only",
    "pilot_at_unchanged_coverage",
)
BASELINES = (ho.UNCHANGED, INDEPENDENT, ho.DISTANCE)
DESCRIPTIVE = (PRIOR, PILOT_COVERAGE)
COMPARISONS = (
    (PILOT, ho.UNCHANGED),
    (PILOT, INDEPENDENT),
    (PILOT, ho.DISTANCE),
    (PILOT, PRIOR),
    (PILOT_COVERAGE, ho.UNCHANGED),
    (INDEPENDENT, ho.DISTANCE),
    (PRIOR, ho.DISTANCE),
    (ho.UNCHANGED, ho.DISTANCE),
)
USEFUL_COVERAGE = 0.80
PASS_ENDPOINTS = 2
PASS_RULE = (
    "the pilot passes if, on at least 2 of the 4 primary endpoints (CRISPRi significant decrease held "
    "out: Gasperini2019, Morris, Schraivogel2020, Xie; validation chromosomes only), its paired "
    "difference in average precision has a 95% locus-bootstrap interval above zero against each of the "
    "three baselines (the unchanged labels, the independent-block variant and distance to TSS) at a "
    "coverage of at least 0.80 of the endpoint's units, and on no primary endpoint is its interval "
    "against any baseline entirely below zero"
)
FALSIFIER = (
    "a pass is void (i) if the prior-only labelling, which reads no holdable source, passes the same "
    "rule against the same three baselines on the same endpoints and the pilot's interval against it "
    "includes zero on them: the gain is then the priors', not the debugger's; (ii) if the pass against "
    "the unchanged labels is coverage alone: the pilot restricted to the unchanged labels' coverage does "
    "not beat them on those endpoints (then the pass against the unchanged labels is reported as "
    "coverage). Any LeakError is fatal and never caught"
)
READINGS = {
    "pass": "an internal development result on withheld sources: joint revision of the labels by the "
    "debugger predicts the held-out CRISPRi decreases better than the unchanged labels, the same "
    "revision one block at a time, and distance to TSS, on the endpoints named, at the coverage given. "
    "Phase C may follow (item 13), and genome-wide deployment waits for a fresh set",
    "beats_unchanged_and_distance_not_independent": "a failure of the coherence pilot: revision helps "
    "but coupling adds nothing measurable; the joint debugger stops as a discontinued investigation and "
    "the per-block revision is a separate finding, to be registered on its own before any use",
    "beats_unchanged_only": "a failure: the revision improves on the compiled labels but not on distance "
    "to TSS; the pilot stops as a discontinued investigation",
    "beats_none": "a failure: the pilot stops and is recorded as a discontinued investigation, apart "
    "from implemented capabilities",
    "void": "a failure of the kind the falsifier names, reported as such; the pilot stops",
}
STOP_RULE = (
    "item 13, binding: if the pilot does not improve prediction on withheld evidence at useful coverage "
    "it stops and is recorded as a discontinued investigation, separately from implemented capabilities. "
    "A higher internal coherence score is never success"
)
METRIC = (
    "validated corrections per compute-hour: a committed correction (a target, activity, split or merge "
    "change of the best alternative that every surviving alternative shares) is validated when, over "
    "the held-out source's units it touches on the validation chromosomes, sum((y - prevalence) x "
    "delta) > 0, where delta is +1 for a unit the correction newly predicts (its new target's pairs, a "
    "context switched on) and -1 for one it stops predicting; an error when that sum is below zero; "
    "neutral at zero; untested when it touches no unit. Compute-hours are the pilot's own CPU time "
    "(building neighbourhoods and searching), measured with time.process_time; the harness's scoring "
    "time is reported beside. Abstentions (neighbourhoods whose best change is not shared by every "
    "survivor, and changes not committed) and observation groups marked inadequate are reported beside"
)
BUDGET_CPU_HOURS = 2.0
BUDGET_RULE = (
    "at most 2 CPU-hours for the pilot's own work over all held-out sources; per neighbourhood the "
    "search caps of pilot.SEARCH apply and a capped neighbourhood is counted. One process, cached data, "
    "no per-element response cache, no model request; peak memory is measured"
)
S4_BESIDE = (
    "descriptive, beside the scores and not in the pass rule: the pilot's committed labels on the "
    "solved blocks, as target and context claims, and the unchanged labels' claims on the same blocks, "
    "each judged by correctness.judge with sources=[S]; target accuracy, role accuracy and coverage "
    "are reported apart, with coverage beside every accuracy"
)
CCRE_TEMPLATE = "data/results/ccres_{chrom}.bed.gz"
PEAK_TEMPLATE = "data/knowledge/epigenome/peaks/{sample}_H3K27ac_{chrom}.bed.gz"
READS_ALWAYS = frozenset({"gencode_v50", "encode_ccre_registry", ho.MODEL, "encode_h3k27ac"})


def registration() -> dict[str, Any]:
    return {
        "registered": REGISTERED,
        "status": STATUS,
        "validation_chromosomes": list(VALIDATION),
        "development_chromosomes": list(DEVELOPMENT),
        "chromosome_choice": CHROMOSOME_CHOICE,
        "primary_endpoints": [list(x) for x in PRIMARY],
        "secondary_endpoints": [list(x) for x in SECONDARY],
        "baselines": list(BASELINES),
        "descriptive_labellings": list(DESCRIPTIVE),
        "comparisons": [f"{a} - {b}" for a, b in COMPARISONS],
        "useful_coverage": USEFUL_COVERAGE,
        "pass_endpoints": PASS_ENDPOINTS,
        "pass_rule": PASS_RULE,
        "falsifier": FALSIFIER,
        "readings": READINGS,
        "stop_rule": STOP_RULE,
        "metric": METRIC,
        "budget_cpu_hours": BUDGET_CPU_HOURS,
        "budget_rule": BUDGET_RULE,
        "s4_beside": S4_BESIDE,
        "validation_endpoint": VALIDATION_ENDPOINT,
        "model": pl.registration(),
    }


# --- the fixed inputs: blocks, compiled links, chromatin, genes -------------------------------------
@dataclass
class Blocks:
    """ENCODE cCREs of one chromosome, searchable by overlap."""

    chrom: str
    ids: list[str]
    starts: list[int]
    ends: list[int]
    reach: int

    @classmethod
    def load(cls, chrom: str, template: str = CCRE_TEMPLATE) -> Blocks:
        rows = []
        with gzip.open(template.format(chrom=chrom), "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split("\t", 4)
                rows.append((int(f[1]), int(f[2]), f[3]))
        rows.sort()
        return cls(
            chrom,
            [r[2] for r in rows],
            [r[0] for r in rows],
            [r[1] for r in rows],
            max((e - s for s, e, _ in rows), default=0),
        )

    def near(self, start: int, end: int) -> list[int]:
        lo = bisect.bisect_left(self.starts, start - self.reach)
        hi = bisect.bisect_left(self.starts, end)
        return [i for i in range(lo, hi) if self.ends[i] > start]

    def root(self, i: int) -> tuple[str, int, int]:
        return (self.ids[i], self.starts[i], self.ends[i])


def load_peaks(chrom: str, template: str = PEAK_TEMPLATE) -> pl.Peaks:
    by: dict[str, list[tuple[int, int]]] = {}
    for s in pl.BREADTH_BIOSAMPLES:
        p = Path(template.format(sample=s, chrom=chrom))
        if not p.exists():
            continue
        rows = []
        with gzip.open(p, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split("\t", 4)
                k = 1 if f[0].startswith("chr") else 0  # the cache keeps start, end, score per chromosome
                rows.append((int(f[k]), int(f[k + 1])))
        by[s] = rows
    return pl.Peaks(by)


def gene_tss(chrom: str) -> dict[str, int]:
    """Symbol -> GENCODE v50 TSS on `chrom`, a protein-coding gene preferred when a symbol repeats."""
    out: dict[str, tuple[int, int]] = {}
    for g in ho.genes():
        if g.chrom != chrom or not g.symbol:
            continue
        rank = 0 if g.gene_type == "protein_coding" else 1
        if g.symbol not in out or rank < out[g.symbol][0]:
            out[g.symbol] = (rank, g.tss)
    return {k: v[1] for k, v in out.items()}


def compiled_targets(links: dict[str, list[ho.Link]], chrom: str) -> dict[str, tuple[str, float, str]]:
    return {x.element: (x.gene, x.strength, x.cell) for x in links.get(chrom, ())}


@dataclass
class Fixed:
    """Everything the pilot reads that no held-out source changes, per validation chromosome."""

    blocks: dict[str, Blocks]
    peaks: dict[str, pl.Peaks]
    tss: dict[str, dict[str, int]]
    links: dict[str, dict[str, tuple[str, float, str]]]

    @classmethod
    def load(cls, chroms: Iterable[str] = VALIDATION) -> Fixed:
        links = ho.compiled_links()
        cs = list(chroms)
        return cls(
            {c: Blocks.load(c) for c in cs},
            {c: load_peaks(c) for c in cs},
            {c: gene_tss(c) for c in cs},
            {c: compiled_targets(links, c) for c in cs},
        )


# --- the view as observations ----------------------------------------------------------------------
def view_obs(view: ho.EvidenceView, fixed: Fixed, chrom: str) -> tuple[list[pl.Obs], Counter]:
    """The held-out view's observations on one chromosome, as the pilot reads them. Counted: what is
    not used (CRISPRi increases and underpowered nulls, GTEx non-associations, saturation-mutagenesis
    bases outside every block)."""
    out: list[pl.Obs] = []
    skipped: Counter = Counter()

    def add(kind: str, s: int, e: int, pos: bool, gene: str = "", ctx: str = "", source: str = "") -> None:
        out.append(pl.make_obs(len(out), kind, chrom, s, e, pos, gene=gene, ctx=ctx, source=source))

    for u in view.units_of("crispri"):
        if u.chrom != chrom:
            continue
        if u.outcome == ms.DECREASE:
            add("crispri", u.start, u.end, True, u.gene, u.cell, u.source)
        elif u.outcome == ms.NULL_INFORMATIVE:
            add("crispri", u.start, u.end, False, u.gene, u.cell, u.source)
        else:
            skipped[f"crispri_{u.outcome or 'no_outcome'}"] += 1
    for u in view.sources.get("gtex", ()):
        if u.chrom != chrom:
            continue
        if u.outcome == "associated":
            add("gtex", u.start, u.end, True, u.gene, "", "gtex")
        else:
            skipped["gtex_not_associated"] += 1
    for src in view.of_assay("lentimpra"):
        for u in view.sources[src]:
            if u.chrom == chrom:
                add("reporter", u.start, u.end, u.outcome == "active", "", u.cell, src)
    for u in view.sources.get("vista", ()):
        if u.chrom == chrom:
            add("vista", u.start, u.end, u.outcome == "positive", "", "invivo", "vista")
    blocks = fixed.blocks[chrom]
    groups: dict[tuple[str, int], list[ho.Unit]] = defaultdict(list)
    for u in view.sources.get("satmut", ()):
        if u.chrom != chrom:
            continue
        hits = blocks.near(u.start, u.end)
        if not hits:
            skipped["satmut_base_outside_every_block"] += 1
        for i in hits:
            groups[(u.study, i)].append(u)
    for _key, us in sorted(groups.items()):
        s = min(u.start for u in us)
        e = max(u.end for u in us)
        add("satmut", s, e, any(u.outcome == "functional" for u in us), "", "other", "satmut")
    return out, skipped


def responsiveness(view: ho.EvidenceView) -> pl.Responsiveness:
    dec: dict[str, list[tuple[int, int]]] = defaultdict(list)
    null: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for u in view.units_of("crispri"):
        if u.outcome == ms.DECREASE:
            dec[u.gene].append((u.start, u.end))
        elif u.outcome == ms.NULL_INFORMATIVE:
            null[u.gene].append((u.start, u.end))
    return pl.Responsiveness(dict(dec), dict(null))


class _ObsIndex:
    def __init__(self, obs: list[pl.Obs]) -> None:
        self.obs = sorted(obs, key=lambda o: (o.start, o.end, o.idx))
        self.starts = [o.start for o in self.obs]
        self.reach = max((o.end - o.start for o in self.obs), default=0)

    def near(self, start: int, end: int) -> list[pl.Obs]:
        lo = bisect.bisect_left(self.starts, start - self.reach)
        hi = bisect.bisect_left(self.starts, end)
        return [o for o in self.obs[lo:hi] if o.end > start]


# --- neighbourhoods and their solutions -------------------------------------------------------------
@dataclass
class Solved:
    """One labelling's solved neighbourhoods on one chromosome, and which one each block belongs to."""

    mode: str
    chrom: str
    hood_of: dict[str, tuple[pl.Hood, pl.Outcome]] = field(default_factory=dict)
    hoods: list[tuple[pl.Hood, pl.Outcome]] = field(default_factory=list)
    cpu: float = 0.0
    reused: int = 0
    solved: int = 0
    capped: int = 0


def relevant_blocks(units: Iterable[ho.Unit], blocks: Blocks) -> list[int]:
    """The blocks a held-out source's units overlap: where its predictions will be asked. Coordinates
    only."""
    out: set[int] = set()
    for u in units:
        out.update(blocks.near(u.start, u.end))
    return sorted(out)


def _components(relevant: list[int], blocks: Blocks, idx: _ObsIndex) -> list[list[int]]:
    """Blocks joined by an observation over two or more of them, or by adjacency within the merge gap,
    reached from the relevant blocks."""
    seen: dict[int, int] = {}
    comps: list[list[int]] = []
    gap = int(pl.SEARCH["merge_gap"])
    for r in relevant:
        if r in seen:
            continue
        comp = [r]
        seen[r] = len(comps)
        k = 0
        while k < len(comp):
            i = comp[k]
            k += 1
            nbrs: set[int] = set()
            for o in idx.near(blocks.starts[i], blocks.ends[i]):
                nbrs.update(blocks.near(o.start, o.end))
            for j in (i - 1, i + 1):
                if 0 <= j < len(blocks.ids):
                    a, b = min(i, j), max(i, j)
                    if 0 <= blocks.starts[b] - blocks.ends[a] <= gap:
                        nbrs.add(j)
            for j in nbrs:
                if j not in seen:
                    seen[j] = len(comps)
                    comp.append(j)
        comps.append(sorted(comp))
    return comps


def _hood_key(
    mode: str, roots: list[tuple[str, int, int]], fixed: set[str], obs: list[pl.Obs], resp: dict[str, float]
) -> str:
    h = hashlib.sha256()
    h.update(repr((mode, roots, sorted(fixed), sorted(resp.items()))).encode())
    for o in obs:
        h.update(repr((o.kind, o.start, o.end, o.positive, o.gene, o.ctx)).encode())
    return h.hexdigest()


def solve(
    mode: str,
    chrom: str,
    relevant: list[int],
    obs: list[pl.Obs],
    fixed: Fixed,
    resp: pl.Responsiveness | None,
    cache: dict[str, pl.Outcome] | None = None,
) -> Solved:
    """Build and solve the neighbourhoods of `mode` that hold the relevant blocks."""
    c0 = time.process_time()
    blocks = fixed.blocks[chrom]
    priors = pl.Priors(
        chrom, fixed.tss[chrom], fixed.links[chrom], fixed.peaks[chrom], resp if mode == "joint" else None
    )
    idx = _ObsIndex(obs)
    out = Solved(mode, chrom)
    cache = cache if cache is not None else {}
    work: list[tuple[list[int], set[int]]] = []
    if mode == "joint":
        cap = int(pl.SEARCH["max_roots"])
        want = set(relevant)
        for comp in _components(relevant, blocks, idx):
            for w in range(0, len(comp), cap):
                window = comp[w : w + cap]
                if not want & set(window):
                    continue
                inside = set(window)
                outside: set[int] = set()
                for i in window:
                    for o in idx.near(blocks.starts[i], blocks.ends[i]):
                        outside.update(j for j in blocks.near(o.start, o.end) if j not in inside)
                work.append((window, outside))
    else:
        work = [([i], set()) for i in relevant]
    for window, outside in work:
        members = sorted(set(window) | outside)
        roots = [blocks.root(i) for i in members]
        seen_o: dict[int, pl.Obs] = {}
        for i in window:
            for o in idx.near(blocks.starts[i], blocks.ends[i]):
                seen_o.setdefault(o.idx, o)
        hood_obs = (
            []
            if mode == "prior"
            else [
                pl.make_obs(
                    k, o.kind, o.chrom, o.start, o.end, o.positive, gene=o.gene, ctx=o.ctx, source=o.source
                )
                for k, o in enumerate(sorted(seen_o.values(), key=lambda o: o.idx))
            ]
        )
        fixed_ids = {blocks.ids[i] for i in outside}
        hood = pl.Hood(chrom, roots, hood_obs, priors, fixed=fixed_ids, mode=mode)
        resp_vals = (
            {f"{g}@{r}": v for g in hood.genes for r, s, e in roots if (v := resp.bonus(g, s, e))}
            if (resp is not None and mode == "joint")
            else {}
        )
        key = _hood_key(mode, roots, fixed_ids, hood_obs, resp_vals)
        res = cache.get(key)
        if res is None:
            res = pl.search(hood)
            cache[key] = res
            out.solved += 1
        else:
            out.reused += 1
            hood = res_hood(res, hood)
        out.capped += int(res.capped)
        out.hoods.append((hood, res))
        for i in window:
            out.hood_of[blocks.ids[i]] = (hood, res)
    out.cpu = time.process_time() - c0
    return out


def res_hood(res: pl.Outcome, hood: pl.Hood) -> pl.Hood:
    """A reused outcome's states name parts the new Hood has not registered yet: register them."""
    states = [res.start, res.best] + [st for st, _ in res.surviving] + [st for st, _ in res.ranked]
    for st in states:
        for r, parts in st.parts.items():
            for p in parts:
                if p not in hood.iv:
                    if "+" in p:
                        a, b = p.split("+", 1)
                        hood.iv[p] = (hood.iv[a][0], max(hood.iv[a][1], hood.iv[b][1]))
                        hood.part_roots[p] = (a, b)
                    else:
                        s, e = (int(x) for x in p.rsplit("|", 1)[1].split("-"))
                        hood.iv[p] = (s, e)
                        hood.part_roots[p] = (r,)
    return hood


# --- labellings for the harness ---------------------------------------------------------------------
def unit_query(u: ho.Unit, endpoint: str) -> dict[str, Any]:
    """What a labelling may read of a held-out unit: its coordinates, gene, TSS and cell."""
    a = u.assay
    if a == "crispri":
        return {"gene": u.gene, "ctx": u.cell, "tss": u.tss}
    if a == "gtex":
        return {"gene": u.gene, "ctx": None, "tss": u.tss}
    if a == "lentimpra":
        return {"gene": None, "ctx": u.cell, "tss": None}
    if a == "vista":
        return {"gene": None, "ctx": "invivo", "tss": None}
    return {"gene": None, "ctx": "other", "tss": None}


def labels_from(
    name: str,
    solved: dict[str, Solved],
    fixed: Fixed,
    reads: frozenset[str],
    built_without: str,
    restrict: ho.Labels | None = None,
    note: str = "",
) -> ho.Labels:
    """A harness labelling from solved neighbourhoods: the largest block score over the blocks a unit
    overlaps, None where it overlaps no solved block. `restrict` abstains wherever that labelling does."""
    memo: dict[tuple, float | None] = {}

    def predict(u: ho.Unit, endpoint: str) -> float | None:
        if restrict is not None and restrict.predict(u, endpoint) is None:
            return None
        sv = solved.get(u.chrom)
        if sv is None:
            return None
        q = unit_query(u, endpoint)
        key = (u.chrom, u.start, u.end, q["gene"], q["ctx"], q["tss"])
        if key in memo:
            return memo[key]
        blocks = fixed.blocks[u.chrom]
        best = None
        for i in blocks.near(u.start, u.end):
            got = sv.hood_of.get(blocks.ids[i])
            if got is None:
                continue
            hood, res = got
            v = pl.root_score(
                hood, res, blocks.ids[i], u.start, u.end, gene=q["gene"], ctx=q["ctx"], tss=q["tss"]
            )
            if v is not None and (best is None or v > best):
                best = v
        memo[key] = best
        return best

    return ho.Labels(name, predict, reads, built_without=built_without, note=note)


# --- the metric: validated corrections -------------------------------------------------------------
def _delta_units(
    change: dict[str, Any],
    hood: pl.Hood,
    res: pl.Outcome,
    units: list[tuple[ho.Unit, float]],
    idx: Any,
) -> list[tuple[ho.Unit, float, int]]:
    """The held-out units a committed change touches, each with +1 (newly predicted) or -1 (no longer)."""
    out: list[tuple[ho.Unit, float, int]] = []
    kind = change["kind"]
    best = res.best
    if kind in ("target", "activity"):
        pid = change["part"]
        s, e = hood.iv[pid]
        for u, y in idx(s, e):
            if kind == "target":
                if u.gene and u.gene == change["now"]:
                    out.append((u, y, 1))
                elif u.gene and u.gene == change["was"]:
                    out.append((u, y, -1))
            else:
                ctx = change["context"]
                sign = 1 if change["now"] else -1
                q = unit_query(u, "")
                element_unit = q["gene"] is None
                pair_on_target = bool(q["gene"]) and u.gene == best.labels[pid].target
                if q["ctx"] == ctx and (element_unit or pair_on_target):
                    out.append((u, y, sign))
        return out
    # a split or a merge: each resulting part, over each block it covers, against that block's own
    # starting label. Correction, 2026-09-29 (lane-pilot), after gate 2's registered run: the run's
    # counter compared a merged part with its first block's starting label only, so the second block's
    # change was never counted; `MERGE_COUNT_AS_RUN` keeps that behaviour for the run's own figures
    start = res.start
    root = change["root"]
    for p in best.parts[root]:
        lab = best.labels[p]
        portions = hood.portions(p)
        if MERGE_COUNT_AS_RUN[0]:
            portions = [(root, hood.iv[p], 1.0)]
        for r, (s, e), _ in portions:
            was = start.labels[r]
            for u, y in idx(s, e):
                q = unit_query(u, "")
                if q["gene"]:
                    if lab.target != was.target:
                        if u.gene == lab.target:
                            out.append((u, y, 1))
                        elif u.gene == was.target:
                            out.append((u, y, -1))
                elif q["ctx"] is not None:
                    now_on, was_on = q["ctx"] in lab.active, q["ctx"] in was.active
                    if now_on != was_on:
                        out.append((u, y, 1 if now_on else -1))
    return out


#: True reproduces gate 2's registered run exactly (a merge counted against its first block only)
MERGE_COUNT_AS_RUN = [False]


def validate(solved: dict[str, Solved], held: list[ho.Unit], endpoint: str) -> dict[str, Any]:
    """METRIC's counts for one held-out source: every committed correction of the pilot's solved
    neighbourhoods, validated against the source's units it touches."""
    ys = [(u, ho.target(u, endpoint)) for u in held]
    ys = [(u, y) for u, y in ys if y is not None]
    prevalence = sum(y for _, y in ys) / len(ys) if ys else 0.0
    by_chrom: dict[str, list[tuple[ho.Unit, float]]] = defaultdict(list)
    for u, y in ys:
        by_chrom[u.chrom].append((u, y))
    tallies: Counter = Counter()
    by_kind: dict[str, Counter] = defaultdict(Counter)
    status: Counter = Counter()
    abstained_changes = 0
    excused = 0
    examples: list[dict[str, Any]] = []
    for chrom, sv in solved.items():
        rows = sorted(by_chrom.get(chrom, []), key=lambda t: t[0].start)
        starts = [u.start for u, _ in rows]
        reach = max((u.end - u.start for u, _ in rows), default=0)

        def idx(s: int, e: int, rows=rows, starts=starts, reach=reach) -> list[tuple[ho.Unit, float]]:
            lo = bisect.bisect_left(starts, s - reach)
            hi = bisect.bisect_left(starts, e)
            return [(u, y) for u, y in rows[lo:hi] if u.end > s]

        seen: set[int] = set()
        for hood, res in sv.hoods:
            if id(res) in seen:
                continue
            seen.add(id(res))
            status[res.status] += 1
            abstained_changes += len(res.abstained)
            excused += len(res.excused)
            for c in res.committed:
                touched = _delta_units(c, hood, res, rows, idx)
                if not touched:
                    verdict = "untested"
                else:
                    score = sum((y - prevalence) * d for _, y, d in touched)
                    verdict = "validated" if score > 1e-12 else ("error" if score < -1e-12 else "neutral")
                tallies[verdict] += 1
                by_kind[c["kind"]][verdict] += 1
                if verdict in ("validated", "error") and len(examples) < 20:
                    examples.append(
                        {
                            "verdict": verdict,
                            "change": c,
                            "touched": [
                                {
                                    "unit": f"{u.chrom}:{u.start}-{u.end}",
                                    "gene": u.gene or None,
                                    "cell": u.cell or None,
                                    "y": y,
                                    "delta": d,
                                }
                                for u, y, d in touched[:6]
                            ],
                        }
                    )
    return {
        "endpoint": endpoint,
        "prevalence": round(prevalence, 4),
        "committed": sum(tallies.values()),
        "validated": tallies["validated"],
        "errors": tallies["error"],
        "neutral": tallies["neutral"],
        "untested": tallies["untested"],
        "by_kind": {k: dict(v) for k, v in sorted(by_kind.items())},
        "neighbourhoods_by_status": dict(status),
        "abstained_changes": abstained_changes,
        "observation_groups_marked_inadequate": excused,
        "examples": examples,
    }


# --- what a disputed neighbourhood returns -----------------------------------------------------------
def neighbourhood_report(hood: pl.Hood, res: pl.Outcome, max_alternatives: int = 5) -> dict[str, Any]:
    """Ranked alternatives, the supporting and conflicting observations of each, the families and the
    next measurement; and the smallest conflicting sets of the heaviest violated observations."""
    st0 = res.start
    v0 = set(hood.violated(st0))
    alts = []
    for st, e in res.ranked[:max_alternatives]:
        v = set(hood.violated(st))
        alts.append(
            {
                "energy": round(e, 4),
                "changes": pl.changes(hood, st0, st),
                "supporting": [hood.obs[i].describe() for i in sorted(v0 - v)],
                "conflicting": [hood.obs[i].describe() for i in sorted(v)],
                "survives": e <= res.best_energy + pl.SEARCH["margin"],
            }
        )
    heavy = sorted(v0, key=lambda i: (-hood.obs[i].weight, i))[: int(pl.SEARCH["mus_per_hood"])]
    return {
        "chrom": hood.chrom,
        "blocks": [r for r in hood.order if r not in hood.fixed],
        "status": res.status,
        "start_energy": round(res.start_energy, 4),
        "best_energy": round(res.best_energy, 4),
        "alternatives": alts,
        "committed": res.committed,
        "abstained": res.abstained,
        "families": [{"size": f["size"], "explains": f["explains"], "of": f["of"]} for f in res.families],
        "next_measurement": res.next_measurement,
        "smallest_conflicting_sets": [pl.mus(hood, st0, i) for i in heavy],
    }


def claims_of(solved: dict[str, Solved]) -> list[Any]:
    """The pilot's committed labels on its solved blocks as S4 claims: a target claim per active cell
    context (or unstated), and a context claim per active cell."""
    from genomeos.attribution import correctness as cx

    out = []
    cells = ("K562", "HepG2", "WTC11")
    for chrom, sv in solved.items():
        seen: set[int] = set()
        for hood, res in sv.hoods:
            if id(res) in seen:
                continue
            seen.add(id(res))
            st = res.best
            for p in hood.all_parts(st):
                if not hood.movable(p):
                    continue
                lab = st.labels[p]
                if not lab.target:
                    continue
                s, e = hood.iv[p]
                on = [c for c in cells if c in lab.active]
                for c in on or [""]:
                    out.append(cx.Claim(p, chrom, s, e, cx.TARGET, lab.target, lab.target, c))
                for c in on:
                    out.append(cx.Claim(p, chrom, s, e, cx.CONTEXT, c, lab.target, c))
    return out


def unchanged_claims_on(solved: dict[str, Solved], claims: list[Any]) -> list[Any]:
    """The unchanged labels' claims on the blocks the pilot solved (their target and context claims)."""
    ids = {r for sv in solved.values() for hood, _ in sv.hoods for r in hood.order if r not in hood.fixed}
    return [c for c in claims if c.element in ids]


def cpu_hours(seconds: float) -> float:
    return seconds / 3600.0


def finite(x: float | None) -> float | None:
    return None if x is None or (isinstance(x, float) and math.isnan(x)) else x


METRIC_KEYS = (
    "committed",
    "validated",
    "errors",
    "neutral",
    "untested",
    "abstained_changes",
    "observation_groups_marked_inadequate",
)


def assess(comparisons: list[dict[str, Any]], scores: list[dict[str, Any]]) -> dict[str, Any]:
    """PASS_RULE, FALSIFIER and READINGS applied to the primary endpoints, as registered."""
    by = {(c["source"], c["endpoint"], c["a"], c["b"]): c for c in comparisons}
    cov = {(s["source"], s["endpoint"], s["labels"]): s["coverage"]["share"] for s in scores}

    def side(src: str, ep: str, a: str, b: str) -> int:
        """+1 when the paired interval is above zero, -1 when below, 0 otherwise or unscored."""
        c = by.get((src, ep, a, b))
        if not c or c.get("status") != "scored" or not c.get("interval"):
            return 0
        lo, hi = c["interval"]
        return 1 if lo > 0 else (-1 if hi < 0 else 0)

    rows = []
    for src, ep in PRIMARY:
        useful = (cov.get((src, ep, PILOT)) or 0) >= USEFUL_COVERAGE
        row = {
            "source": src,
            "endpoint": ep,
            "pilot_coverage": cov.get((src, ep, PILOT)),
            "useful_coverage": useful,
            "above": {b: side(src, ep, PILOT, b) > 0 for b in BASELINES},
            "below": {b: side(src, ep, PILOT, b) < 0 for b in BASELINES},
            "pilot_above_prior": side(src, ep, PILOT, PRIOR) > 0,
            "prior_above": {b: side(src, ep, PRIOR, b) > 0 for b in (ho.UNCHANGED, ho.DISTANCE)},
            "pilot_at_unchanged_coverage_above_unchanged": side(src, ep, PILOT_COVERAGE, ho.UNCHANGED) > 0,
        }
        row["passes"] = bool(useful and all(row["above"].values()))
        rows.append(row)
    n_pass = sum(r["passes"] for r in rows)
    worse = [r["source"] for r in rows if any(r["below"].values())]
    passed = n_pass >= PASS_ENDPOINTS and not worse
    passing = [r for r in rows if r["passes"]]
    prior_explains = passed and all(
        not r["pilot_above_prior"] and all(r["prior_above"].values()) for r in passing
    )
    coverage_only = passed and all(not r["pilot_at_unchanged_coverage_above_unchanged"] for r in passing)
    if passed and prior_explains:
        reading = "void"
    elif passed:
        reading = "pass"
    else:
        count = {b: sum(1 for r in rows if r["useful_coverage"] and r["above"][b]) for b in BASELINES}
        if count[ho.UNCHANGED] >= PASS_ENDPOINTS and count[ho.DISTANCE] >= PASS_ENDPOINTS:
            reading = "beats_unchanged_and_distance_not_independent"
        elif count[ho.UNCHANGED] >= PASS_ENDPOINTS:
            reading = "beats_unchanged_only"
        else:
            reading = "beats_none"
    return {
        "endpoints": rows,
        "endpoints_passing": n_pass,
        "endpoints_worse_than_a_baseline": worse,
        "passed": reading == "pass",
        "reading": reading,
        "reading_text": READINGS[reading],
        "falsifier": {
            "prior_only_explains_the_pass": prior_explains,
            "pass_against_unchanged_is_coverage": coverage_only,
        },
        "stop_rule_fires": reading != "pass",
    }
