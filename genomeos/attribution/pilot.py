# SPDX-License-Identifier: AGPL-3.0-or-later
"""The coherence pilot: a debugger, boundaries and families over local neighbourhoods (item 13 C1, C2
and C3; item 12 S3).

The question it asks is not "is this block X or Y" one block at a time, but which smallest set of
assumptions cannot explain the observations of a neighbourhood together, which repairs remove the
contradiction, and which repairs the evidence cannot tell apart. It is scored on evidence it did not
see (genomeos/attribution/holdout.py, C4) and by nothing else: a higher internal coherence score is
never success.

**The model.** A neighbourhood holds blocks (ENCODE cCREs; a block can be split or merged), genes with
their TSS, and observations. A block's label is one target gene or none, and the set of contexts it is
active in. Observations are weighted soft constraints over the blocks their interval overlaps, read as
an OR: a CRISPRi decrease of gene g in cell c is explained when some overlapping block targets g and is
active in c; a reporter tile, a transgenic embryo or a saturation-mutagenesis locus is explained when
some overlapping block is active in its context; a GTEx association when some overlapping block
targets its gene. A well-powered null is violated when that holds. Priors: distance to each candidate
gene's TSS (a power law, not fitted), the compiled predicted target (the unchanged labels), and the
annotation's activity (H3K27ac peaks). Hard constraints are only the mandatory ones: a part lies inside
its block, parts do not overlap, no part is shorter than `SEARCH["min_part"]`. All weights are
registered constants (`W`); nothing is fitted.

**C1, the debugger.** From the starting labels, the violated observations name the conflicting
blocks; `mus` reduces a violated observation to a smallest conflicting set of observations and label
assumptions (deletion-based, exact over the involved blocks); repairs are generated only on the
variables of the conflict (change the target, change the activity context, split the block, merge it
with its neighbour, or mark an observation group's model inadequate) and evaluated singly and in pairs,
the second move drawn from the conflict the first leaves. Only the affected blocks' energy is
recomputed. The best repair is applied while it lowers the energy; up to `SEARCH["rounds"]` rounds.

**C2, boundaries.** A split is proposed only where evidence changes inside a block (an observation's
edge or an H3K27ac peak's edge), costs `W["split"]` (fragmentation), and is counted as a correction
that must improve prediction of withheld evidence like any other. A merge joins adjacent blocks under
one label at `W["merge"]`.

**C3, families.** Alternatives within `SEARCH["margin"]` of the best survive. Surviving alternatives
that predict the same value for every observation of the neighbourhood form one family: the evidence
supports their shared relationship but cannot say which of them supplies it. A change of the best
alternative is committed only if every surviving alternative shares it; otherwise the pilot abstains
and names the cheapest measurement (`MEASUREMENT_COST`) whose predicted outcome differs among the
survivors, by expected information per unit cost.

**What it never uses.** R8's two retired score transformations (element competition and gene budget,
`joint_pretest.RETIRED`): no share of a score over genes or over elements is computed anywhere here; a
pair's score is its link odds against "no target", never a fraction of other genes' or blocks' scores.
No model request; the per-element response cache is not opened; the only model output read is the
compiled predicted layer, which is what the unchanged labels are made of.
"""

from __future__ import annotations

import bisect
import hashlib
import math
import random
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from typing import Any

# ==================================================================================================
# The registration of the model and the search (2026-09-29, lane-pilot), fixed before the synthetic gate
# was run and before any score of any pilot labelling was read. docs/ATTRIBUTION.md carries the text.
# ==================================================================================================
REGISTERED = "2026-09-29"
STATUS = "internal development result on withheld sources"
NONE = ""  # no target
#: the contexts every block states an activity for; an observation may name another cell, which then
#: becomes a context of that neighbourhood
CONTEXTS = ("K562", "HepG2", "WTC11", "invivo", "other")
#: cells with their own experimental chromatin in the local cache (ENCODE H3K27ac replicated peaks)
OWN_CHROMATIN = {"K562": "K562", "HepG2": "HepG2"}
#: every H3K27ac biosample in the local cache; the share with a peak is a block's activity breadth
BREADTH_BIOSAMPLES = (
    "CD14_positive_monocyte",
    "GM12878",
    "H1",
    "HepG2",
    "IMR_90",
    "K562",
    "SK_N_SH",
    "astrocyte",
    "cardiac_muscle_cell",
    "hepatocyte",
    "keratinocyte",
    "ovary",
    "testis",
)
#: the compiled model cells that name one of the contexts
MODEL_CELL_CONTEXT = {"K562": "K562", "HepG2": "HepG2", "WTC11": "WTC11"}
#: weights in nats (log-odds), registered and not fitted
W: dict[str, float] = {
    # observations: the cost of leaving one unexplained (positive) or contradicted (negative)
    "crispri_decrease": 3.0,  # a significant decrease against a well-powered screen
    "crispri_null": 1.5,  # a well-powered null: weaker, knockdown can be incomplete
    "gtex_association": 0.5,  # association is not perturbation; linkage; any tissue
    "reporter_active": 1.0,  # episomal activity (R6: one tile)
    "reporter_inactive": 0.3,  # an inactive episomal tile says little about the endogenous locus
    "vista_positive": 1.0,
    "vista_negative": 0.5,
    "satmut_functional": 1.0,
    "satmut_inert": 0.3,
    # gene-level coupling (joint only): a gene the evidence saw respond to CRISPRi at another element
    "responsive": 1.0,
    "unresponsive": 0.5,  # at least `unresponsive_min_nulls` well-powered nulls elsewhere, no decrease
    "unresponsive_min_nulls": 3,
    # activity prior, logit per context
    "activity_prior": -1.0,
    "own_peak": 2.0,
    "own_no_peak": -1.0,
    "breadth": 2.0,  # times the share of the 13 biosamples with a peak; contexts without own chromatin
    "model_activity": 0.5,  # the compiled link's model cell is this context
    # target prior, log-weight
    "model_target_base": 1.0,  # the compiled target: base + min(strength, cap)
    "model_target_strength_cap": 1.0,
    "distance_gamma": 1.0,  # -gamma * ln(1 + d / d0): a contact power law, not fitted
    "distance_d0": 5000.0,
    "no_target": -3.0,  # the log-weight of "no target": that of a gene at about 100 kb
    # structural repairs
    "split": 2.0,  # fragmentation
    "merge": 1.0,
    "inadequate": 2.0,  # per (block, observation kind) excused
}
OBS_WEIGHT = {
    ("crispri", True): "crispri_decrease",
    ("crispri", False): "crispri_null",
    ("gtex", True): "gtex_association",
    ("reporter", True): "reporter_active",
    ("reporter", False): "reporter_inactive",
    ("vista", True): "vista_positive",
    ("vista", False): "vista_negative",
    ("satmut", True): "satmut_functional",
    ("satmut", False): "satmut_inert",
}
#: the observation kinds whose model can be marked inadequate at a block. CRISPRi is the endpoint of
#: regulation and is never excused; its weaker null weight already carries incomplete knockdown
INADEQUATE_KINDS = ("reporter", "vista", "satmut", "gtex")
SEARCH: dict[str, float] = {
    "rounds": 6,
    "margin": 1.0,  # alternatives within this many nats of the best survive
    "improve": 0.05,  # a repair is applied only if it lowers the energy by more than this
    "min_part": 50,
    "merge_gap": 500,
    "max_roots": 60,  # a larger component is cut into windows, the outside blocks held fixed
    "max_singles": 400,
    "max_pairs": 20_000,
    "mus_per_hood": 3,  # smallest conflicting sets computed for the heaviest violated observations
    "mus_max_assignments": 50_000,
}
MEASUREMENT_COST = {"reporter": 1.0, "crispri": 4.0}


def registration() -> dict[str, Any]:
    return {
        "registered": REGISTERED,
        "status": STATUS,
        "contexts": list(CONTEXTS),
        "own_chromatin": OWN_CHROMATIN,
        "breadth_biosamples": list(BREADTH_BIOSAMPLES),
        "model_cell_context": MODEL_CELL_CONTEXT,
        "weights": W,
        "inadequate_kinds": list(INADEQUATE_KINDS),
        "search": SEARCH,
        "measurement_cost": MEASUREMENT_COST,
        "retired_never_used": "element competition and gene budget (joint_pretest.RETIRED); no share of "
        "a score over genes or elements is computed",
    }


# --- the pieces -----------------------------------------------------------------------------------
@dataclass(frozen=True, slots=True)
class Obs:
    """One observation as a soft constraint over the blocks its interval overlaps."""

    idx: int
    kind: str  # crispri | gtex | reporter | vista | satmut
    chrom: str
    start: int
    end: int
    positive: bool
    gene: str = ""
    ctx: str = ""
    source: str = ""
    weight: float = 0.0

    def describe(self) -> dict[str, Any]:
        return {
            "kind": self.kind,
            "interval": f"{self.chrom}:{self.start}-{self.end}",
            "gene": self.gene or None,
            "context": self.ctx or None,
            "outcome": "positive" if self.positive else "negative",
            "source": self.source,
            "weight": self.weight,
        }


def make_obs(idx: int, kind: str, chrom: str, start: int, end: int, positive: bool, **kw: Any) -> Obs:
    return Obs(idx, kind, chrom, start, end, positive, weight=W[OBS_WEIGHT[(kind, positive)]], **kw)


@dataclass(frozen=True, slots=True)
class Label:
    target: str
    active: frozenset[str]


@dataclass
class State:
    """One assignment of a neighbourhood: its parts per block, their labels, and excused groups."""

    parts: dict[str, tuple[str, ...]]  # root block -> its current parts (a merge maps two roots to one)
    labels: dict[str, Label]
    excused: frozenset[tuple[str, str]] = frozenset()  # (root, observation kind)
    nsplit: int = 0
    nmerge: int = 0
    moves: tuple[tuple, ...] = ()

    def key(self) -> tuple:
        return (
            tuple(sorted(self.parts.items())),
            tuple(sorted((p, lab.target, tuple(sorted(lab.active))) for p, lab in self.labels.items())),
            tuple(sorted(self.excused)),
        )


def _overlaps(s1: int, e1: int, s2: int, e2: int) -> bool:
    return s1 < e2 and s2 < e1


class Peaks:
    """Peak intervals per biosample on one chromosome, searchable by overlap."""

    def __init__(self, by_sample: dict[str, list[tuple[int, int]]]) -> None:
        self.by = {k: sorted(v) for k, v in by_sample.items()}
        self.starts = {k: [s for s, _ in v] for k, v in self.by.items()}
        self.reach = {k: max((e - s for s, e in v), default=0) for k, v in self.by.items()}

    def hits(self, sample: str, start: int, end: int) -> list[tuple[int, int]]:
        v = self.by.get(sample)
        if not v:
            return []
        st = self.starts[sample]
        lo = bisect.bisect_left(st, start - self.reach[sample])
        hi = bisect.bisect_left(st, end)
        return [(s, e) for s, e in v[lo:hi] if e > start and s < end]

    def any(self, sample: str, start: int, end: int) -> bool:
        return bool(self.hits(sample, start, end))


class Responsiveness:
    """Gene-level coupling from measured CRISPRi (joint mode only): whether the evidence saw a gene
    respond at another element. The element asked about is excluded, so direct evidence at a block is
    never counted twice."""

    def __init__(
        self, decreases: dict[str, list[tuple[int, int]]], nulls: dict[str, list[tuple[int, int]]]
    ) -> None:
        self.dec = decreases
        self.null = nulls

    def bonus(self, gene: str, start: int, end: int) -> float:
        if not gene:
            return 0.0
        if any(not _overlaps(s, e, start, end) for s, e in self.dec.get(gene, ())):
            return W["responsive"]
        n = sum(1 for s, e in self.null.get(gene, ()) if not _overlaps(s, e, start, end))
        if n >= W["unresponsive_min_nulls"]:
            return -W["unresponsive"]
        return 0.0


@dataclass
class Priors:
    """What a block's label costs before any holdable observation: distance, the compiled target and
    the annotation's activity."""

    chrom: str
    tss: dict[str, int]  # gene symbol -> TSS on this chromosome (0-based)
    links: dict[str, tuple[str, float, str]]  # root -> (gene, strength, model cell)
    peaks: Peaks
    responsiveness: Responsiveness | None = None

    @staticmethod
    def distance_logw(mid: int, tss: int) -> float:
        return -W["distance_gamma"] * math.log1p(abs(mid - tss) / W["distance_d0"])

    def link_of(self, roots: tuple[str, ...]) -> tuple[str, float, str] | None:
        best = None
        for r in roots:
            x = self.links.get(r)
            if x is not None and (best is None or x[1] > best[1]):
                best = x
        return best

    def target_logw(
        self, iv: tuple[int, int], roots: tuple[str, ...], gene: str, tss: int | None = None
    ) -> float:
        if gene == NONE:
            return W["no_target"]
        t = tss if tss is not None else self.tss.get(gene)
        mid = (iv[0] + iv[1]) // 2
        lw = self.distance_logw(mid, t) if t is not None else self.distance_logw(0, 10**9)
        link = self.link_of(roots)
        if link is not None and link[0] == gene:
            lw += W["model_target_base"] + min(link[1], W["model_target_strength_cap"])
        if self.responsiveness is not None:
            lw += self.responsiveness.bonus(gene, iv[0], iv[1])
        return lw

    def act_logit(self, iv: tuple[int, int], roots: tuple[str, ...], ctx: str) -> float:
        z = W["activity_prior"]
        own = OWN_CHROMATIN.get(ctx)
        if own is not None and own in self.peaks.by:
            z += W["own_peak"] if self.peaks.any(own, iv[0], iv[1]) else W["own_no_peak"]
        else:
            n = sum(1 for s in BREADTH_BIOSAMPLES if self.peaks.any(s, iv[0], iv[1]))
            z += W["breadth"] * n / len(BREADTH_BIOSAMPLES)
        link = self.link_of(roots)
        if link is not None and MODEL_CELL_CONTEXT.get(link[2]) == ctx:
            z += W["model_activity"]
        return z

    def peak_edges(self, iv: tuple[int, int]) -> set[int]:
        out: set[int] = set()
        for s in OWN_CHROMATIN.values():
            for a, b in self.peaks.hits(s, iv[0], iv[1]):
                out.update((a, b))
        return out


# --- one neighbourhood ----------------------------------------------------------------------------
class Hood:
    """A neighbourhood: blocks (roots), observations over them, priors, and the energy of a state.

    `mode` is "joint" (the pilot), "independent" (one block, its own copy of every observation over
    it, no coupling: one Hood per block) or "prior" (no observations). `fixed` roots are held at their
    starting labels: they are the outside blocks a window cut of a large component keeps, so an
    observation that crosses the cut is still evaluated against both sides."""

    def __init__(
        self,
        chrom: str,
        roots: list[tuple[str, int, int]],
        obs: list[Obs],
        priors: Priors,
        *,
        fixed: Iterable[str] = (),
        mode: str = "joint",
    ) -> None:
        self.chrom = chrom
        self.mode = mode
        self.priors = priors
        self.order = [r for r, _, _ in sorted(roots, key=lambda x: (x[1], x[2], x[0]))]
        self.iv: dict[str, tuple[int, int]] = {r: (s, e) for r, s, e in roots}
        self.part_roots: dict[str, tuple[str, ...]] = {r: (r,) for r, _, _ in roots}
        self.fixed = set(fixed)
        self.obs = [] if mode == "prior" else list(obs)
        self.obs_roots: dict[int, tuple[str, ...]] = {}
        self.root_obs: dict[str, list[int]] = {r: [] for r in self.iv}
        starts = sorted((s, e, r) for r, (s, e) in self.iv.items())
        st_list = [x[0] for x in starts]
        reach = max((e - s for s, e, _ in starts), default=0)
        for i, o in enumerate(self.obs):
            lo = bisect.bisect_left(st_list, o.start - reach)
            hi = bisect.bisect_left(st_list, o.end)
            rs = tuple(r for s, e, r in starts[lo:hi] if _overlaps(s, e, o.start, o.end))
            self.obs_roots[i] = rs
            for r in rs:
                self.root_obs[r].append(i)
        linked = {x[0] for r in self.iv if (x := priors.links.get(r)) is not None}
        self.genes = sorted({o.gene for o in self.obs if o.gene} | linked)
        self.ctxs = sorted(set(CONTEXTS) | {o.ctx for o in self.obs if o.ctx})
        self._tlw: dict[tuple, float] = {}
        self._alg: dict[tuple, float] = {}
        self._mass: dict[str, list[tuple[str, tuple[int, int], float]]] = {}

    # priors, memoised per part
    def portions(self, pid: str) -> list[tuple[str, tuple[int, int], float]]:
        """A part's share of each block it covers: (block, the covered interval, the covered share).
        A part's prior is the share-weighted sum of each block's own prior over its covered interval,
        so a split conserves a block's prior and a merge neither adds nor drops one, and no block
        lends its compiled target to its neighbour: structure is decided by evidence and the
        structural costs."""
        v = self._mass.get(pid)
        if v is None:
            s, e = self.iv[pid]
            v = []
            for r in self.part_roots[pid]:
                rs, re_ = self.iv[r]
                a, b = max(s, rs), min(e, re_)
                if b > a:
                    v.append((r, (a, b), (b - a) / max(1, re_ - rs)))
            self._mass[pid] = v
        return v

    def tlogw(self, pid: str, gene: str, tss: int | None = None) -> float:
        k = (pid, gene, tss)
        v = self._tlw.get(k)
        if v is None:
            v = sum(f * self.priors.target_logw(iv, (r,), gene, tss) for r, iv, f in self.portions(pid))
            self._tlw[k] = v
        return v

    def alogit(self, pid: str, ctx: str) -> float:
        k = (pid, ctx)
        v = self._alg.get(k)
        if v is None:
            v = sum(f * self.priors.act_logit(iv, (r,), ctx) for r, iv, f in self.portions(pid))
            self._alg[k] = v
        return v

    def prior_energy(self, pid: str, lab: Label, tss: int | None = None) -> float:
        e = -self.tlogw(pid, lab.target, tss if lab.target else None)
        for c in lab.active:
            e -= self.alogit(pid, c)
        return e

    def s0(self) -> State:
        """The starting labels: the compiled target and the activity the annotation implies."""
        labels = {}
        for r in self.order:
            link = self.priors.links.get(r)
            target = link[0] if link is not None else NONE
            active = frozenset(c for c in self.ctxs if self.alogit(r, c) > 0)
            labels[r] = Label(target, active)
        return State(parts={r: (r,) for r in self.order}, labels=labels)

    # observations
    def parts_of_obs(self, st: State, i: int) -> list[str]:
        o = self.obs[i]
        out: list[str] = []
        for r in self.obs_roots[i]:
            for p in st.parts[r]:
                if p not in out and _overlaps(*self.iv[p], o.start, o.end):
                    out.append(p)
        return out

    def predicts(self, st: State, i: int) -> bool:
        o = self.obs[i]
        return any(_explains(o, st.labels[p]) for p in self.parts_of_obs(st, i))

    def excused(self, st: State, i: int) -> bool:
        if not st.excused:
            return False
        k = self.obs[i].kind
        return any((r, k) in st.excused for r in self.obs_roots[i])

    def obs_cost(self, st: State, i: int) -> float:
        if self.excused(st, i):
            return 0.0
        o = self.obs[i]
        p = self.predicts(st, i)
        if o.positive:
            return 0.0 if p else o.weight
        return o.weight if p else 0.0

    def obs_touching(self, pids: Iterable[str]) -> list[int]:
        out: set[int] = set()
        for p in pids:
            s, e = self.iv[p]
            for r in self.part_roots[p]:
                for i in self.root_obs[r]:
                    o = self.obs[i]
                    if _overlaps(s, e, o.start, o.end):
                        out.add(i)
        return sorted(out)

    def all_parts(self, st: State) -> list[str]:
        seen: dict[str, None] = {}
        for r in self.order:
            for p in st.parts[r]:
                seen.setdefault(p, None)
        return list(seen)

    @staticmethod
    def structural(st: State) -> float:
        return W["split"] * st.nsplit + W["merge"] * st.nmerge + W["inadequate"] * len(st.excused)

    def energy(self, st: State) -> float:
        e = sum(self.prior_energy(p, st.labels[p]) for p in self.all_parts(st))
        e += sum(self.obs_cost(st, i) for i in range(len(self.obs)))
        return e + self.structural(st)

    def local(self, st: State, pids: Iterable[str]) -> float:
        pids = list(dict.fromkeys(pids))
        e = sum(self.prior_energy(p, st.labels[p]) for p in pids)
        return e + sum(self.obs_cost(st, i) for i in self.obs_touching(pids))

    def local_with(self, st: State, pid: str, lab: Label, tss: int | None = None) -> float:
        """The local energy of `pid` under label `lab`, everything else as in `st` (restored after)."""
        old = st.labels[pid]
        st.labels[pid] = lab
        try:
            e = self.prior_energy(pid, lab, tss)
            return e + sum(self.obs_cost(st, i) for i in self.obs_touching([pid]))
        finally:
            st.labels[pid] = old

    def violated(self, st: State) -> list[int]:
        return [i for i in range(len(self.obs)) if self.obs_cost(st, i) > 0]

    def fits(self, st: State) -> bool:
        return not st.excused and not self.violated(st)

    def signature(self, st: State) -> tuple:
        return tuple((self.predicts(st, i), self.excused(st, i)) for i in range(len(self.obs)))

    def movable(self, pid: str) -> bool:
        return not any(r in self.fixed for r in self.part_roots[pid])

    def conflict_parts(self, st: State, viol: Iterable[int]) -> list[str]:
        out: dict[str, None] = {}
        for i in viol:
            for p in self.parts_of_obs(st, i):
                if self.movable(p):
                    out.setdefault(p, None)
        return list(out)

    # repairs
    def breakpoints(self, pid: str, touching: list[int]) -> list[int]:
        """Where evidence changes inside the part: an observation's edge or an H3K27ac peak's edge."""
        s, e = self.iv[pid]
        lo, hi = s + SEARCH["min_part"], e - SEARCH["min_part"]
        cand: set[int] = set()
        for i in touching:
            o = self.obs[i]
            cand.update(x for x in (o.start, o.end) if lo <= x <= hi)
        cand.update(x for x in self.priors.peak_edges((s, e)) if lo <= x <= hi)
        return sorted(cand)

    def adjacent(self, st: State, pid: str) -> list[str]:
        if pid not in self.order or st.parts[pid] != (pid,):
            return []
        i = self.order.index(pid)
        out = []
        for j in (i - 1, i + 1):
            if 0 <= j < len(self.order):
                q = self.order[j]
                if st.parts[q] == (q,) and self.movable(q):
                    a, b = sorted((pid, q), key=lambda r: self.iv[r][0])
                    if 0 <= self.iv[b][0] - self.iv[a][1] <= SEARCH["merge_gap"]:
                        out.append(q)
        return out

    def moves(self, st: State, pids: Iterable[str], full: bool = False) -> list[tuple]:
        """Repairs on the given parts. The pilot proposes only the genes and contexts the parts'
        observations name (and the compiled target, and no target); `full` proposes every gene and
        context of the neighbourhood (the brute-force checker of the synthetic gate)."""
        out: list[tuple] = []
        for p in pids:
            if not self.movable(p) or p not in st.labels:
                continue
            lab = st.labels[p]
            touching = self.obs_touching([p])
            if full:
                genes = set(self.genes) | {NONE}
                ctxs = set(self.ctxs)
            else:
                genes = {self.obs[i].gene for i in touching if self.obs[i].gene} | {NONE}
                link = self.priors.link_of(self.part_roots[p])
                if link is not None:
                    genes.add(link[0])
                ctxs = {self.obs[i].ctx for i in touching if self.obs[i].ctx}
            for g in sorted(genes - {lab.target}):
                out.append(("target", p, g))
            for c in sorted(ctxs):
                out.append(("act", p, c))
            if len(self.part_roots[p]) == 1:
                for x in self.breakpoints(p, touching):
                    out.append(("split", p, x))
            for q in self.adjacent(st, p):
                out.append(("merge", p, q))
            kinds = {
                self.obs[i].kind
                for i in touching
                if self.obs[i].kind in INADEQUATE_KINDS and (full or self.obs_cost(st, i) > 0)
            }
            for k in sorted(kinds):
                for r in self.part_roots[p]:
                    if (r, k) not in st.excused:
                        out.append(("excuse", r, k))
        return list(dict.fromkeys(out))

    def apply(self, st: State, m: tuple) -> tuple[State, list[str], list[str]]:
        """The state after move `m`, with the parts it touched before and after."""
        kind, a, b = m
        parts = dict(st.parts)
        labels = dict(st.labels)
        excused, nsplit, nmerge = st.excused, st.nsplit, st.nmerge
        if kind == "target":
            labels[a] = Label(b, labels[a].active)
            before = after = [a]
        elif kind == "act":
            lab = labels[a]
            labels[a] = Label(lab.target, lab.active ^ {b})
            before = after = [a]
        elif kind == "split":
            (root,) = self.part_roots[a]
            s, e = self.iv[a]
            left, right = f"{root}|{s}-{b}", f"{root}|{b}-{e}"
            for pid, iv in ((left, (s, b)), (right, (b, e))):
                self.iv.setdefault(pid, iv)
                self.part_roots.setdefault(pid, (root,))
            lab = labels.pop(a)
            labels[left] = labels[right] = lab
            ps = list(parts[root])
            k = ps.index(a)
            parts[root] = tuple(ps[:k] + [left, right] + ps[k + 1 :])
            nsplit += 1
            before, after = [a], [left, right]
        elif kind == "merge":
            p, q = sorted((a, b), key=lambda r: self.iv[r][0])
            mid = f"{p}+{q}"
            self.iv.setdefault(mid, (self.iv[p][0], max(self.iv[p][1], self.iv[q][1])))
            self.part_roots.setdefault(mid, (p, q))
            lab = labels[a]
            labels.pop(p)
            labels.pop(q)
            labels[mid] = lab
            parts[p] = parts[q] = (mid,)
            nmerge += 1
            before, after = [p, q], [mid]
        elif kind == "excuse":
            excused = excused | {(a, b)}
            before = after = list(st.parts[a])
        else:
            raise ValueError(f"unknown move {m!r}")
        return State(parts, labels, excused, nsplit, nmerge, (*st.moves, m)), before, after

    def step(self, st: State, e: float, m: tuple) -> tuple[State, float]:
        st2, before, after = self.apply(st, m)
        de = self.local(st2, after) - self.local(st, before) + self.structural(st2) - self.structural(st)
        return st2, e + de


def _explains(o: Obs, lab: Label) -> bool:
    if o.kind == "crispri":
        return lab.target == o.gene and o.ctx in lab.active
    if o.kind == "gtex":
        return lab.target == o.gene
    return o.ctx in lab.active


# --- C1: the smallest conflicting set -------------------------------------------------------------
def mus(hood: Hood, st: State, i: int) -> dict[str, Any]:
    """A smallest set of observations and label assumptions, containing the violated observation `i`,
    that no assignment of the involved blocks can satisfy together. Deletion-based and exact by
    enumeration over the involved parts, their breakpoints, the genes and contexts the observations
    name; blocks outside are held as in `st`."""
    parts = [p for p in hood.parts_of_obs(st, i) if hood.movable(p)]
    touching = hood.obs_touching(parts)
    others = [j for j in touching if j != i and not hood.excused(st, j)]
    named = [i, *others]
    ctxs = sorted({hood.obs[j].ctx for j in named if hood.obs[j].ctx})
    genes = sorted(
        {hood.obs[j].gene for j in named if hood.obs[j].gene} | {st.labels[p].target for p in parts} | {NONE}
    )
    asm: list[tuple] = []
    for p in parts:
        asm.append(("target", p, st.labels[p].target))
        for c in ctxs:
            asm.append(("act", p, c, c in st.labels[p].active))
        if len(hood.part_roots[p]) == 1 and hood.breakpoints(p, touching):
            asm.append(("boundary", p))
    labs = [
        Label(g, frozenset(c for c, on in zip(ctxs, bits, strict=True) if on))
        for g in genes
        for bits in _bits(len(ctxs))
    ]
    budget = [int(SEARCH["mus_max_assignments"])]

    def structures(p: str, keep: bool) -> list[list[tuple[str, tuple[int, int]]]]:
        out = [[(p, hood.iv[p])]]
        if not keep and len(hood.part_roots[p]) == 1:
            s, e = hood.iv[p]
            for x in hood.breakpoints(p, touching):
                out.append([(f"{p}|{s}-{x}", (s, x)), (f"{p}|{x}-{e}", (x, e))])
        return out

    def consistent(obs_set: list[int], asm_set: list[tuple]) -> bool | None:
        fixed_target = {a[1]: a[2] for a in asm_set if a[0] == "target"}
        fixed_act = {(a[1], a[2]): a[3] for a in asm_set if a[0] == "act"}
        keep = {a[1] for a in asm_set if a[0] == "boundary"}
        for combo in _product([structures(p, p in keep) for p in parts]):
            flat = [(p, iv) for p, subs in zip(parts, combo, strict=True) for _, iv in subs]
            choices = [
                [
                    lab
                    for lab in labs
                    if fixed_target.get(p, lab.target) == lab.target
                    and all(fixed_act.get((p, c), c in lab.active) == (c in lab.active) for c in ctxs)
                ]
                for p, _ in flat
            ]
            for pick in _product(choices):
                budget[0] -= 1
                if budget[0] < 0:
                    return None
                assign = [(iv, lab) for (_, iv), lab in zip(flat, pick, strict=True)]
                if all(_satisfied(hood, st, j, parts, assign) for j in obs_set):
                    return True
        return False

    first = consistent(list(named), list(asm))
    if first is None:
        return {"minimal": False, "why": "enumeration cap reached", "observation": hood.obs[i].describe()}
    if first:
        return {
            "minimal": False,
            "why": "no conflict among the involved blocks",
            "observation": hood.obs[i].describe(),
        }

    def shrink(order: str) -> tuple[list[int], list[tuple]] | None:
        core_obs, core_asm = list(named), list(asm)
        steps = [("o", j) for j in others] + [("a", a) for a in asm]
        if order == "assumptions_first":
            steps = [("a", a) for a in asm] + [("o", j) for j in others]
        for what, x in steps:
            t_obs = [y for y in core_obs if y != x] if what == "o" else core_obs
            t_asm = [y for y in core_asm if y != x] if what == "a" else core_asm
            r = consistent(t_obs, t_asm)
            if r is None:
                return None
            if not r:
                core_obs, core_asm = t_obs, t_asm
        return core_obs, core_asm

    def describe(core: tuple[list[int], list[tuple]]) -> dict[str, Any]:
        core_obs, core_asm = core
        pos_genes = {hood.obs[j].gene for j in core_obs if hood.obs[j].positive and hood.obs[j].gene}
        rules = ["a block has one target"] if len(pos_genes) > 1 else []
        return {
            "observations": [hood.obs[j].describe() for j in core_obs],
            "assumptions": [_asm_text(a) for a in core_asm],
            "rules": rules,
            "size": len(core_obs) + len(core_asm) + len(rules),
        }

    small = shrink("observations_first")
    among = shrink("assumptions_first")
    if small is None or among is None:
        return {"minimal": False, "why": "enumeration cap reached", "observation": hood.obs[i].describe()}
    out = {"minimal": True, "observation": hood.obs[i].describe(), **describe(small)}
    b = describe(among)
    if b["size"] < out["size"]:
        out, b = {"minimal": True, "observation": hood.obs[i].describe(), **b}, describe(small)
    out["conflict_among_observations"] = b
    return out


def _asm_text(a: tuple) -> str:
    if a[0] == "target":
        return f"{a[1]} targets {a[2] or 'no gene'}"
    if a[0] == "act":
        return f"{a[1]} is {'active' if a[3] else 'inactive'} in {a[2]}"
    return f"{a[1]} is one block (its boundary is kept)"


def _bits(n: int) -> list[tuple[bool, ...]]:
    return [tuple(bool(k >> j & 1) for j in range(n)) for k in range(2**n)]


def _product(lists: list[list[Any]]) -> Iterator[tuple]:
    if not lists:
        yield ()
        return
    head, *rest = lists
    for x in head:
        for tail in _product(rest):
            yield (x, *tail)


def _satisfied(hood: Hood, st: State, j: int, parts: list[str], assign: list[tuple]) -> bool:
    o = hood.obs[j]
    hit = False
    for r in hood.obs_roots[j]:
        for p in st.parts[r]:
            if p not in parts and _overlaps(*hood.iv[p], o.start, o.end) and _explains(o, st.labels[p]):
                hit = True
    for iv, lab in assign:
        if _overlaps(iv[0], iv[1], o.start, o.end) and _explains(o, lab):
            hit = True
    return hit if o.positive else not hit


# --- the search, the families and the next measurement --------------------------------------------
@dataclass
class Outcome:
    """What the pilot returns for one neighbourhood."""

    disputed: bool
    start: State
    start_energy: float
    best: State
    best_energy: float
    surviving: list[tuple[State, float]] = field(default_factory=list)
    ranked: list[tuple[State, float]] = field(default_factory=list)
    status: str = "undisputed"  # undisputed | kept | resolved | abstained | excused
    committed: list[dict[str, Any]] = field(default_factory=list)
    abstained: list[dict[str, Any]] = field(default_factory=list)
    excused: list[tuple[str, str]] = field(default_factory=list)
    families: list[dict[str, Any]] = field(default_factory=list)
    next_measurement: dict[str, Any] | None = None
    evaluations: int = 0
    capped: bool = False

    def weights(self) -> list[tuple[State, float]]:
        if not self.surviving:
            return [(self.best, 1.0)]
        z = [math.exp(-(e - self.best_energy)) for _, e in self.surviving]
        t = sum(z)
        return [(s, w / t) for (s, _), w in zip(self.surviving, z, strict=True)]


def search(hood: Hood, pairs: bool = True) -> Outcome:
    """C1's repair search from the starting labels: singles and, if `pairs`, pairs of repairs on the
    conflict, the best applied while it lowers the energy by more than SEARCH["improve"]."""
    st0 = hood.s0()
    e0 = hood.energy(st0)
    if not hood.violated(st0):
        return Outcome(False, st0, e0, st0, e0, surviving=[(st0, e0)], ranked=[(st0, e0)])
    seen: dict[tuple, tuple[State, float]] = {st0.key(): (st0, e0)}
    st, e = st0, e0
    evals = 0
    capped = False
    for _ in range(int(SEARCH["rounds"])):
        viol = hood.violated(st)
        if not viol:
            break
        singles = hood.moves(st, hood.conflict_parts(st, viol))
        if len(singles) > SEARCH["max_singles"]:
            singles, capped = singles[: int(SEARCH["max_singles"])], True
        cands: list[tuple[State, float]] = []
        firsts = []
        for m in singles:
            s1, e1 = hood.step(st, e, m)
            evals += 1
            cands.append((s1, e1))
            firsts.append((m, s1, e1))
        if pairs:
            n_pairs = 0
            for m, s1, e1 in firsts:
                if n_pairs >= SEARCH["max_pairs"]:
                    capped = True
                    break
                c1 = hood.conflict_parts(s1, hood.violated(s1))
                c1 = list(dict.fromkeys(c1 + [p for p in _touched(hood, s1, m) if hood.movable(p)]))
                for m2 in hood.moves(s1, c1):
                    if _same_variable(m, m2):
                        continue
                    s2, e2 = hood.step(s1, e1, m2)
                    evals += 1
                    n_pairs += 1
                    cands.append((s2, e2))
                    if n_pairs >= SEARCH["max_pairs"]:
                        break
        if not cands:
            break
        for s, x in cands:
            k = s.key()
            if k not in seen or x < seen[k][1]:
                seen[k] = (s, x)
        best = min(cands, key=lambda t: (t[1], len(t[0].moves)))
        if best[1] < e - SEARCH["improve"]:
            st, e = best
        else:
            break
    return summarise(hood, st0, e0, st, e, list(seen.values()), evals, capped)


def _touched(hood: Hood, st: State, m: tuple) -> list[str]:
    kind, a, _ = m
    if kind == "split":
        (root,) = hood.part_roots[a]
        return list(st.parts[root])
    if kind in ("merge", "excuse"):
        return list(st.parts[a])
    return [a]


def _same_variable(m1: tuple, m2: tuple) -> bool:
    if m1[0] != m2[0]:
        return False
    if m1[0] == "act":
        return m1[1] == m2[1] and m1[2] == m2[2]
    if m1[0] == "target":
        return m1[1] == m2[1]
    return m1 == m2


def root_view(hood: Hood, st: State, root: str) -> tuple:
    return tuple(
        (hood.iv[p], st.labels[p].target, tuple(sorted(st.labels[p].active))) for p in st.parts[root]
    )


def changes(hood: Hood, a: State, b: State) -> list[dict[str, Any]]:
    """What differs from `a` to `b`, block by block: structure, target, activity; and excused groups."""
    out: list[dict[str, Any]] = []
    for r in hood.order:
        if r in hood.fixed:
            continue
        pa, pb = a.parts[r], b.parts[r]
        if pa != pb:
            kind = "merge" if any("+" in p for p in pb) else "split"
            if kind == "merge" and pb[0].split("+")[0] != r:
                continue  # a merge is reported once, under its first block
            out.append(
                {
                    "kind": kind,
                    "root": r,
                    "parts": [
                        {"interval": list(iv), "target": t or None, "active": list(act)}
                        for iv, t, act in root_view(hood, b, r)
                    ],
                    "was": [
                        {"target": t or None, "active": list(act)} for _, t, act in root_view(hood, a, r)
                    ],
                }
            )
            continue
        for p in pb:
            la, lb = a.labels[p], b.labels[p]
            if la.target != lb.target:
                out.append(
                    {
                        "kind": "target",
                        "root": r,
                        "part": p,
                        "was": la.target or None,
                        "now": lb.target or None,
                    }
                )
            for c in sorted(la.active ^ lb.active):
                out.append({"kind": "activity", "root": r, "part": p, "context": c, "now": c in lb.active})
    for r, k in sorted(b.excused - a.excused):
        out.append({"kind": "inadequate", "root": r, "observations": k})
    return out


def change_key(c: dict[str, Any]) -> tuple:
    if c["kind"] in ("split", "merge"):
        return (
            c["kind"],
            c["root"],
            tuple((tuple(x["interval"]), x["target"], tuple(x["active"])) for x in c["parts"]),
        )
    if c["kind"] == "target":
        return ("target", c["part"], c["now"])
    if c["kind"] == "activity":
        return ("activity", c["part"], c["context"], c["now"])
    return ("inadequate", c["root"], c["observations"])


def summarise(
    hood: Hood,
    st0: State,
    e0: float,
    best: State,
    eb: float,
    explored: list[tuple[State, float]],
    evals: int,
    capped: bool,
) -> Outcome:
    """Rank what was explored, keep the survivors, split the best's changes into committed and
    abstained, and name the families and the next measurement."""
    ranked = sorted(explored, key=lambda t: (t[1], len(t[0].moves), t[0].key()))
    surviving = [(s, x) for s, x in ranked if x <= eb + SEARCH["margin"]]
    best_changes = changes(hood, st0, best)
    others = [{change_key(c) for c in changes(hood, st0, s)} for s, _ in surviving if s.key() != best.key()]
    committed, abstained = [], []
    for c in best_changes:
        if c["kind"] == "inadequate":
            continue
        k = change_key(c)
        (committed if all(k in oc for oc in others) else abstained).append(c)
    alternatives_change = any(oc - {k for k in oc if k[0] == "inadequate"} for oc in others)
    if not best_changes:
        status = "abstained" if alternatives_change else "kept"
    elif abstained:
        status = "abstained"
    elif committed:
        status = "resolved"
    else:
        status = "excused"
    out = Outcome(
        True,
        st0,
        e0,
        best,
        eb,
        surviving=surviving,
        ranked=ranked[:10],
        status=status,
        committed=committed,
        abstained=abstained,
        excused=sorted(best.excused),
        evaluations=evals,
        capped=capped,
    )
    out.families = families(hood, st0, surviving)
    out.next_measurement = next_measurement(hood, out)
    return out


def families(hood: Hood, st0: State, surviving: list[tuple[State, float]]) -> list[dict[str, Any]]:
    """Surviving alternatives grouped by what they predict for every observation of the neighbourhood:
    members of one family are the ones the evidence cannot tell apart."""
    groups: dict[tuple, list[tuple[State, float]]] = defaultdict(list)
    for s, x in surviving:
        groups[hood.signature(s)].append((s, x))
    out = []
    for sig, members in sorted(groups.items(), key=lambda kv: min(x for _, x in kv[1])):
        out.append(
            {
                "size": len(members),
                "members": [{"energy": round(x, 4), "changes": changes(hood, st0, s)} for s, x in members],
                "explains": sum(1 for i, (p, exc) in enumerate(sig) if not exc and p == hood.obs[i].positive),
                "of": len(sig),
            }
        )
    return out


def next_measurement(hood: Hood, out: Outcome) -> dict[str, Any] | None:
    """The cheapest candidate measurement (a reporter tile or a CRISPRi pair over a block or a part of
    it) whose predicted outcome differs between the best alternative and its strongest surviving
    rival, ties broken by expected information over all survivors; when no candidate separates those
    two, the one with the most expected information per unit cost over all survivors. None when one
    alternative survives or no candidate separates any."""
    ws = out.weights()
    if len(ws) < 2:
        return None
    first, rival = ws[0][0], max(ws[1:], key=lambda t: t[1])[0]
    roots: dict[str, None] = {}
    for s, _ in ws:
        for c in changes(hood, out.start, s):
            for r in hood.part_roots.get(c.get("part", c["root"]), (c["root"],)):
                roots.setdefault(r, None)
    intervals: dict[tuple[str, int, int], None] = {}
    for s, _ in ws:
        for r in roots:
            for p in s.parts[r]:
                ps, pe = hood.iv[p]
                intervals.setdefault((r, ps, pe), None)
            rs, re_ = hood.iv[r]
            intervals.setdefault((r, rs, re_), None)
    genes = sorted(
        {s.labels[p].target for s, _ in ws for r in roots for p in s.parts[r] if s.labels[p].target}
        | {hood.obs[i].gene for r in roots for i in hood.root_obs[r] if hood.obs[i].gene}
    )
    ctxs = sorted({hood.obs[i].ctx for r in roots for i in hood.root_obs[r] if hood.obs[i].ctx} or {"K562"})
    best = None
    for r, s_, e_ in intervals:
        for c in ctxs:
            for kind, g in [("reporter", "")] + [("crispri", g) for g in genes]:
                probe = Obs(-1, kind, hood.chrom, s_, e_, True, gene=g, ctx=c)
                p = sum(w for s, w in ws if _predict_probe(hood, s, probe))
                if p <= 1e-9 or p >= 1 - 1e-9:
                    continue
                h = -(p * math.log2(p) + (1 - p) * math.log2(1 - p))
                cost = MEASUREMENT_COST[kind]
                separates = _predict_probe(hood, first, probe) != _predict_probe(hood, rival, probe)
                key = (0 if separates else 1, cost if separates else -h / cost, -h, kind, s_, e_, g, c)
                if best is None or key < best[0]:
                    best = (
                        key,
                        {
                            "measure": kind,
                            "block": r,
                            "interval": f"{hood.chrom}:{s_}-{e_}",
                            "gene": g or None,
                            "context": c,
                            "cost": cost,
                            "bits": round(h, 4),
                            "separates_best_from_rival": separates,
                            "p_positive_under_survivors": round(p, 4),
                        },
                    )
    return best[1] if best else None


def _predict_probe(hood: Hood, st: State, o: Obs) -> bool:
    return any(
        _overlaps(*hood.iv[p], o.start, o.end) and _explains(o, st.labels[p]) for p in hood.all_parts(st)
    )


# --- scores for the harness: link odds and activity, mixed over the survivors ---------------------
def link_probability(hood: Hood, st: State, pid: str, gene: str, ctx: str | None, tss: int | None) -> float:
    """P(pid targets `gene`, and is active in `ctx` when one is named), against "no target" (and
    inactive), everything else as in `st`: a contrast with no target, never a share over other genes."""
    lab = st.labels[pid]
    acts: tuple[bool | None, ...] = (True, False) if ctx else (None,)
    es = []
    for t in (gene, NONE):
        for a in acts:
            active = lab.active if a is None else (lab.active | {ctx} if a else lab.active - {ctx})
            es.append(hood.local_with(st, pid, Label(t, active), tss if t == gene else None))
    m = min(es)
    z = [math.exp(-(x - m)) for x in es]
    return z[0] / sum(z)


def activity_probability(hood: Hood, st: State, pid: str, ctx: str) -> float:
    lab = st.labels[pid]
    e1 = hood.local_with(st, pid, Label(lab.target, lab.active | {ctx}))
    e0 = hood.local_with(st, pid, Label(lab.target, lab.active - {ctx}))
    return 1.0 / (1.0 + math.exp(max(-50.0, min(50.0, e1 - e0))))


def root_score(
    hood: Hood,
    out: Outcome,
    root: str,
    start: int,
    end: int,
    *,
    gene: str | None = None,
    ctx: str | None = None,
    tss: int | None = None,
) -> float | None:
    """A block's score for a unit over [start, end): a pair (gene, and ctx when the endpoint is
    cell-specific) or an activity (ctx), mixed over the surviving alternatives; the largest over the
    block's parts the unit touches."""
    total = 0.0
    hit = False
    for st, w in out.weights():
        best = None
        for p in st.parts[root]:
            if not _overlaps(*hood.iv[p], start, end):
                continue
            v = (
                link_probability(hood, st, p, gene, ctx, tss)
                if gene
                else activity_probability(hood, st, p, ctx or "")
            )
            best = v if best is None else max(best, v)
        if best is not None:
            total += w * best
            hit = True
    return total if hit else None


# ==================================================================================================
# Gate 1: the synthetic gate (item 12 S3). Registered 2026-09-29 before it was run.
# ==================================================================================================
SYNTH_SEED = 20260929
SYNTH_INSTANCES = 100  # valid generated instances per motif, beside instance 0, which is hand-built
SYNTH_MAX_DRAWS = 5000  # draws per motif before the generator gives up
SYNTH_PASS_SHARE = 0.95
SYNTH_MOTIFS = ("split_retarget", "swap", "context_handoff")
SYNTH_FAMILY_MOTIF = "indistinguishable"
SYNTH_CHROM = "chrS"
SYNTH_PASS_RULE = (
    "on each of the three hand-built cases (which must also pass the validity check) and on at least 95% "
    "of the 100 generated instances of each of the three two-correction motifs, the pilot's returned "
    "labels equal the planted truth and both planted corrections are committed (no abstention on them); "
    "and on the hand-built case and at least 95% of the 100 instances of the family motif, the two "
    "indistinguishable alternatives both survive in one family, the pilot abstains on which block "
    "supplies the link, and its next measurement's predicted outcome differs between them. If the gate "
    "fails, the pilot stops and nothing biological is scored"
)
SYNTH_VALIDITY = (
    "a two-correction instance is kept only if, by brute force over every single repair and every pair "
    "of repairs in the full move space (every gene of the neighbourhood and no target, every context, "
    "every breakpoint, merge and excuse), (1) no single repair fits the evidence, (2) no single repair "
    "lowers the energy by more than SEARCH['improve'], so a one-repair-at-a-time search is stuck at the "
    "starting labels, (3) the planted pair fits, and (4) the planted pair's energy is below every other "
    "alternative of at most two repairs by more than SEARCH['margin']. A family instance is kept only "
    "if both single alternatives fit, predict the same for every observation, and lie within "
    "SEARCH['margin'] of each other and of the best alternative of at most two repairs"
)
SYNTH_CONTROLS = (
    "descriptive, not part of the pass rule: the same search without pairs (greedy, one repair at a "
    "time); the independent-block variant (each block its own neighbourhood with its own copy of every "
    "observation over it, no coupling); and a noisy copy of each instance with one extra random "
    "observation"
)
SYNTH_FALSIFIER = (
    "the gate is void if the greedy control recovers the planted truth on any kept instance: the "
    "validity check would then not be doing its job"
)


@dataclass
class Instance:
    motif: str
    roots: list[tuple[str, int, int]]
    obs: list[Obs]
    priors: Priors
    truth_moves: tuple[tuple, ...]
    hand_built: bool
    params: dict[str, Any]

    def hood(self, mode: str = "joint") -> Hood:
        return Hood(SYNTH_CHROM, list(self.roots), list(self.obs), self.priors, mode=mode)


def _add(obs: list[Obs], kind: str, s: int, e: int, positive: bool, gene: str = "", ctx: str = "") -> None:
    obs.append(make_obs(len(obs), kind, SYNTH_CHROM, s, e, positive, gene=gene, ctx=ctx, source="synthetic"))


def build_instance(motif: str, rng: random.Random | None) -> Instance:
    """One instance of a motif; `rng=None` gives the hand-built instance 0."""

    def u(lo: float, hi: float, fixed: float) -> float:
        return fixed if rng is None else rng.uniform(lo, hi)

    def ri(lo: int, hi: int, fixed: int) -> int:
        return fixed if rng is None else rng.randint(lo, hi)

    base = 1_000_000
    obs: list[Obs] = []
    peaks: dict[str, list[tuple[int, int]]] = {"K562": [], "HepG2": []}
    roots: list[tuple[str, int, int]] = []
    links: dict[str, tuple[str, float, str]] = {}
    tss: dict[str, int] = {}
    if motif == "split_retarget":
        # one block whose two halves regulate different genes: split it and retarget one half
        length = ri(1200, 3000, 2000)
        cut = base + int(length * u(0.35, 0.65, 0.5))
        roots.append(("B", base, base + length))
        d_a, d_b = ri(15_000, 120_000, 40_000), ri(15_000, 120_000, 40_000)
        tss.update({"GA": base - d_a, "GB": base + length + d_b})
        links["B"] = ("GA", u(0.2, 0.8, 0.5), "K562")
        peaks["K562"].append((base - 100, base + length + 100))
        _add(obs, "crispri", base, cut, True, "GB", "K562")
        _add(obs, "crispri", base, cut, False, "GA", "K562")
        _add(obs, "gtex", base, cut, True, "GB")
        _add(obs, "crispri", cut, base + length, True, "GA", "K562")
        _add(obs, "crispri", cut, base + length, False, "GB", "K562")
        _add(obs, "reporter", base, cut, True, ctx="K562")
        _add(obs, "reporter", cut, base + length, True, ctx="K562")
        _decoy(obs, roots, links, tss, peaks, base + length + ri(3_000, 20_000, 8_000), "GC", rng)
        truth = (("split", "B", cut), ("target", f"B|{base}-{cut}", "GB"))
        params = {"length": length, "cut": cut - base, "d_a": d_a, "d_b": d_b, "strength": links["B"][1]}
    elif motif == "swap":
        # two blocks under two screens that each span both: exchange their targets together
        l1, gap, l2 = ri(200, 600, 350), ri(100, 800, 300), ri(200, 600, 350)
        b1 = (base, base + l1)
        b2 = (base + l1 + gap, base + l1 + gap + l2)
        roots += [("B1", *b1), ("B2", *b2)]
        d_a, d_b = ri(10_000, 150_000, 50_000), ri(10_000, 150_000, 50_000)
        tss.update({"GA": base - d_a, "GB": b2[1] + d_b})
        links["B1"] = ("GA", u(0.6, 1.5, 1.0), "K562")
        links["B2"] = ("GB", u(0.6, 1.5, 1.0), "K562")
        peaks["K562"] += [(b1[0] - 50, b1[1] + 50), (b2[0] - 50, b2[1] + 50)]
        span = (b1[0] - ri(50, 400, 200), b2[1] + ri(50, 400, 200))
        _add(obs, "crispri", *span, True, "GA", "K562")
        _add(obs, "crispri", *span, True, "GB", "K562")
        _add(obs, "crispri", *b1, False, "GA", "K562")
        _add(obs, "crispri", *b1, True, "GB", "K562")
        _add(obs, "crispri", *b2, True, "GA", "K562")
        _add(obs, "crispri", *b2, False, "GB", "K562")
        _decoy(obs, roots, links, tss, peaks, b2[1] + ri(3_000, 20_000, 8_000), "GC", rng)
        truth = (("target", "B1", "GB"), ("target", "B2", "GA"))
        params = {"l1": l1, "gap": gap, "l2": l2, "d_a": d_a, "d_b": d_b}
    elif motif == "context_handoff":
        # a block that is not active in K562 after all, and its neighbour that carries the link instead
        l1, gap, l2 = ri(200, 600, 350), ri(100, 800, 300), ri(200, 600, 350)
        b1 = (base, base + l1)
        b2 = (base + l1 + gap, base + l1 + gap + l2)
        roots += [("B1", *b1), ("B2", *b2)]
        d_a = ri(165_000, 215_000, 190_000)
        tss.update({"GA": b2[1] + d_a, "GB": base - ri(20_000, 200_000, 80_000)})
        links["B1"] = ("GA", u(0.5, 1.2, 0.8), "HepG2")
        peaks["K562"] += [(b1[0] - 50, b1[1] + 50), (b2[0] - 50, b2[1] + 50)]
        peaks["HepG2"].append((b1[0] - 50, b1[1] + 50))
        span = (b1[0] - ri(50, 400, 200), b2[1] + ri(50, 400, 200))
        third = (b1[1] - b1[0]) // 3
        _add(obs, "crispri", *span, True, "GA", "K562")
        _add(obs, "crispri", *b1, False, "GA", "K562")
        _add(obs, "gtex", *b1, True, "GA")
        _add(obs, "reporter", b1[0], b1[0] + third, False, ctx="K562")
        _add(obs, "reporter", b1[0] + third, b1[0] + 2 * third, False, ctx="K562")
        _add(obs, "reporter", b1[0] + 2 * third, b1[1], False, ctx="K562")
        _add(obs, "reporter", *b1, True, ctx="HepG2")
        _add(obs, "reporter", *b2, True, ctx="K562")
        _add(obs, "gtex", *b2, True, "GA")
        truth = (("act", "B1", "K562"), ("target", "B2", "GA"))
        params = {"l1": l1, "gap": gap, "l2": l2, "d_a": d_a}
    elif motif == SYNTH_FAMILY_MOTIF:
        # one screen spanning two blocks, nothing else to tell which of them regulates the gene
        l1, gap, l2 = ri(200, 600, 350), ri(100, 800, 300), ri(200, 600, 350)
        b1 = (base, base + l1)
        b2 = (base + l1 + gap, base + l1 + gap + l2)
        roots += [("B1", *b1), ("B2", *b2)]
        d_a = ri(20_000, 200_000, 60_000)
        tss.update({"GA": (b1[0] + b2[1]) // 2 - d_a, "GC": b2[1] + ri(300_000, 600_000, 400_000)})
        links["B1"] = ("GC", u(0.05, 0.3, 0.1), "placenta")
        peaks["K562"] += [(b1[0] - 50, b1[1] + 50), (b2[0] - 50, b2[1] + 50)]
        span = (b1[0] - ri(50, 400, 200), b2[1] + ri(50, 400, 200))
        _add(obs, "crispri", *span, True, "GA", "K562")
        truth = (("target", "B1", "GA"), ("target", "B2", "GA"))  # the two single alternatives
        params = {"l1": l1, "gap": gap, "l2": l2, "d_a": d_a}
    else:
        raise ValueError(motif)
    priors = Priors(SYNTH_CHROM, tss, links, Peaks(peaks), None)
    return Instance(motif, roots, obs, priors, truth, rng is None, params)


def _decoy(
    obs: list[Obs],
    roots: list,
    links: dict,
    tss: dict,
    peaks: dict,
    at: int,
    gene: str,
    rng: random.Random | None,
) -> None:
    """A block beside the motif whose own evidence agrees with its labels: nothing to repair there."""
    length = 300 if rng is None else rng.randint(200, 500)
    roots.append(("D", at, at + length))
    tss.setdefault(gene, at + (20_000 if rng is None else rng.randint(5_000, 40_000)))
    links["D"] = (gene, 1.0 if rng is None else rng.uniform(0.6, 1.5), "K562")
    peaks["K562"].append((at - 50, at + length + 50))
    _add(obs, "crispri", at, at + length, True, gene, "K562")


def all_alternatives(hood: Hood, st0: State, e0: float) -> list[tuple[State, float]]:
    """Every alternative of at most two repairs in the full move space, by brute force."""
    out = [(st0, e0)]
    for m in hood.moves(st0, hood.all_parts(st0), full=True):
        s1, e1 = hood.step(st0, e0, m)
        out.append((s1, e1))
        for m2 in hood.moves(s1, hood.all_parts(s1), full=True):
            if not _same_variable(m, m2):
                out.append(hood.step(s1, e1, m2))
    return out


def truth_state(hood: Hood, inst: Instance) -> tuple[State, float]:
    st = hood.s0()
    e = hood.energy(st)
    for m in inst.truth_moves:
        st, e = hood.step(st, e, m)
    return st, e


def valid(inst: Instance) -> tuple[bool, str]:
    """SYNTH_VALIDITY, by brute force."""
    hood = inst.hood()
    st0 = hood.s0()
    e0 = hood.energy(st0)
    alts = all_alternatives(hood, st0, e0)
    if inst.motif == SYNTH_FAMILY_MOTIF:
        a = [hood.step(st0, e0, m) for m in inst.truth_moves]
        if not all(hood.fits(s) for s, _ in a):
            return False, "a planted single does not fit"
        if hood.signature(a[0][0]) != hood.signature(a[1][0]):
            return False, "the two alternatives predict differently"
        best = min(x for _, x in alts)
        if abs(a[0][1] - a[1][1]) > SEARCH["margin"] or max(a[0][1], a[1][1]) > best + SEARCH["margin"]:
            return False, "the two alternatives are not both within the margin of the best"
        return True, "ok"
    singles = [(s, x) for s, x in alts if len(s.moves) == 1]
    if any(hood.fits(s) for s, _ in singles):
        return False, "a single repair fits"
    if any(x < e0 - SEARCH["improve"] for _, x in singles):
        return False, "a single repair lowers the energy"
    ts, te = truth_state(hood, inst)
    if not hood.fits(ts):
        return False, "the planted pair does not fit"
    tk = ts.key()
    rival = min((x for s, x in alts if s.key() != tk), default=math.inf)
    if not te < rival - SEARCH["margin"]:
        return False, "the planted pair is not the best alternative by the margin"
    return True, "ok"


def generate(motif: str, n: int = SYNTH_INSTANCES, seed: int = SYNTH_SEED) -> tuple[list[Instance], Counter]:
    """Instance 0 (hand-built, kept whatever its validity, which is reported) and `n` valid random
    instances; rejected draws are counted by reason."""
    rng = random.Random(f"{seed}:{motif}")
    out = [build_instance(motif, None)]
    rejected: Counter = Counter()
    draws = 0
    while len(out) < n + 1 and draws < SYNTH_MAX_DRAWS:
        draws += 1
        inst = build_instance(motif, rng)
        ok, why = valid(inst)
        if ok:
            out.append(inst)
        else:
            rejected[why] += 1
    rejected["draws"] = draws
    return out, rejected


def solve_independent(inst: Instance) -> tuple[Hood, State]:
    """The independent-block control: each block its own neighbourhood with its own copy of every
    observation over it; the final labels are the union of the per-block answers."""
    joint = inst.hood()
    st = joint.s0()
    parts = dict(st.parts)
    labels = dict(st.labels)
    nsplit = 0
    for rid, s, e in inst.roots:
        own = [o for o in inst.obs if _overlaps(s, e, o.start, o.end)]
        own = [
            make_obs(i, o.kind, o.chrom, o.start, o.end, o.positive, gene=o.gene, ctx=o.ctx, source=o.source)
            for i, o in enumerate(own)
        ]
        h = Hood(SYNTH_CHROM, [(rid, s, e)], own, inst.priors, mode="independent")
        out = search(h)
        for p in parts[rid]:
            labels.pop(p, None)
        parts[rid] = out.best.parts[rid]
        for p in parts[rid]:
            labels[p] = out.best.labels[p]
            joint.iv.setdefault(p, h.iv[p])
            joint.part_roots.setdefault(p, h.part_roots[p])
        nsplit += out.best.nsplit
    return joint, State(parts, labels, frozenset(), nsplit)


def _labels_equal(hood: Hood, a: State, b: State) -> bool:
    for r in hood.order:
        if tuple(a.parts[r]) != tuple(b.parts[r]):
            return False
        if any(a.labels[p] != b.labels[p] for p in a.parts[r]):
            return False
    return True


def _family_check(hood: Hood, inst: Instance, out: Outcome) -> dict[str, Any]:
    st0 = hood.s0()
    e0 = hood.energy(st0)
    alts = [hood.step(st0, e0, m)[0] for m in inst.truth_moves]
    surv = {s.key(): s for s, _ in out.surviving}
    both = all(a.key() in surv for a in alts)
    one_family = False
    want = [tuple(sorted(change_key(c) for c in changes(hood, st0, a))) for a in alts]
    for fam in out.families:
        keys = {tuple(sorted(change_key(c) for c in mem["changes"])) for mem in fam["members"]}
        if all(w in keys for w in want):
            one_family = True
    abstains = out.status == "abstained" and not any(c["kind"] == "target" for c in out.committed)
    nm = out.next_measurement
    separates = False
    if nm is not None and both:
        s_, e_ = (int(x) for x in nm["interval"].split(":")[1].split("-"))
        probe = Obs(-1, nm["measure"], hood.chrom, s_, e_, True, gene=nm["gene"] or "", ctx=nm["context"])
        separates = len({_predict_probe(hood, a, probe) for a in alts}) == 2
    return {
        "both_survive": both,
        "one_family": one_family,
        "abstains": abstains,
        "next_measurement": nm,
        "separates": separates,
        "passed": bool(both and one_family and abstains and separates),
    }


def _noisy(inst: Instance, k: int, seed: int) -> Instance:
    """A copy with one extra random observation over one of the blocks."""
    rng = random.Random(f"{seed}:noise:{inst.motif}:{k}")
    obs = list(inst.obs)
    s, e = rng.choice([(x[1], x[2]) for x in inst.roots])
    kind = rng.choice(("crispri", "reporter", "gtex"))
    gene = rng.choice(sorted(inst.priors.tss)) if kind != "reporter" else ""
    positive = True if kind == "gtex" else rng.random() < 0.5
    ctx = "" if kind == "gtex" else "K562"
    obs.append(make_obs(len(obs), kind, SYNTH_CHROM, s, e, positive, gene=gene, ctx=ctx, source="noise"))
    return Instance(inst.motif, inst.roots, obs, inst.priors, inst.truth_moves, inst.hand_built, inst.params)


def run_synthetic_gate(n: int = SYNTH_INSTANCES, seed: int = SYNTH_SEED) -> dict[str, Any]:
    """Gate 1, as registered: the pilot on every kept instance, the controls beside it."""
    per: dict[str, Any] = {}
    digests: dict[str, str] = {}
    for motif in (*SYNTH_MOTIFS, SYNTH_FAMILY_MOTIF):
        insts, rejected = generate(motif, n, seed)
        digests[motif] = instance_digest(insts)
        rows = []
        for k, inst in enumerate(insts):
            hood = inst.hood()
            ok_valid, why = valid(inst) if k == 0 else (True, "ok")
            out = search(hood)
            row: dict[str, Any] = {
                "instance": k,
                "hand_built": inst.hand_built,
                "valid": ok_valid,
                "why": why,
            }
            if motif == SYNTH_FAMILY_MOTIF:
                row.update(_family_check(hood, inst, out))
            else:
                ts, _ = truth_state(hood, inst)
                planted = {change_key(c) for c in changes(hood, hood.s0(), ts)}
                row["found"] = out.best.key() == ts.key()
                row["committed_planted"] = planted <= {change_key(c) for c in out.committed}
                row["passed"] = bool(row["found"] and row["committed_planted"])
                g = search(inst.hood(), pairs=False)
                row["greedy_found"] = g.best.key() == ts.key()
                jh, ind = solve_independent(inst)
                row["independent_found"] = _labels_equal(jh, ind, ts)
                nout = search(_noisy(inst, k, seed).hood())
                row["noisy_found"] = nout.best.key() == ts.key()
            row["passed"] = bool(row["passed"] and row["valid"])
            row["evaluations"] = out.evaluations
            row["status"] = out.status
            rows.append(row)
        per[motif] = _summarise_motif(motif, rows, rejected, n)
    passed = all(per[m]["passed"] for m in (*SYNTH_MOTIFS, SYNTH_FAMILY_MOTIF))
    void = any(per[m]["greedy_found"] > 0 for m in SYNTH_MOTIFS)
    return {
        "motifs": per,
        "instance_digests": digests,
        "passed": passed and not void,
        "void": void,
        "verdict": "void" if void else ("pass" if passed else "fail"),
    }


def _summarise_motif(motif: str, rows: list[dict[str, Any]], rejected: Counter, n: int) -> dict[str, Any]:
    hand, gen = rows[0], rows[1:]
    passed_gen = sum(1 for r in gen if r["passed"])
    share = passed_gen / len(gen) if gen else 0.0
    out: dict[str, Any] = {
        "hand_built": hand,
        "generated": len(gen),
        "generated_passed": passed_gen,
        "generated_share": round(share, 4),
        "rejected_draws": dict(rejected),
        "passed": bool(hand["passed"] and len(gen) >= n and share >= SYNTH_PASS_SHARE),
        "evaluations_median": sorted(r["evaluations"] for r in rows)[len(rows) // 2],
        "status": dict(Counter(r["status"] for r in rows)),
        "failures": [r for r in gen if not r["passed"]][:5],
    }
    if motif != SYNTH_FAMILY_MOTIF:
        out["of"] = len(rows)
        out["greedy_found"] = sum(1 for r in rows if r["greedy_found"])
        out["independent_found"] = sum(1 for r in rows if r["independent_found"])
        out["noisy_found"] = sum(1 for r in rows if r["noisy_found"])
    return out


def instance_digest(insts: list[Instance]) -> str:
    """A digest of the instances, so a rerun can show it scored the same ones."""
    h = hashlib.sha256()
    for i in insts:
        h.update(repr((i.motif, sorted(i.params.items()), i.truth_moves)).encode())
    return h.hexdigest()
