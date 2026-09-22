# SPDX-License-Identifier: Apache-2.0
# Part of the BioLang engine (language, IR, VM, standard library); see LICENSING.md.
"""CpG methylation as a state machine the VM *runs*, division by division.

The epigenome layer of GenomeOS reads methylation as an annotation: a fraction per
200 bp bin. This engine instead runs the published mechanism that produces that
fraction, so a program can predict how a methylome changes as a lineage divides.

Per CpG site (a dyad: the C on each strand), three states:

    U   unmethylated       neither strand carries 5mC
    H   hemimethylated     one strand carries it
    M   methylated         both strands carry it

One cell division, in the order the enzymes act:

1. **Replication.** Each daughter keeps one parental strand and gets one new,
   unmethylated one, so M -> H. An H dyad's methylated strand goes to one
   daughter only: following a single lineage, H -> H or U with equal
   probability. U -> U. (Watson & Crick semiconservative copying; the
   consequence for methylation is Holliday & Pugh 1975, Science 187:226.)
2. **Maintenance.** UHRF1 reads the hemimethylated dyad and recruits DNMT1,
   which restores it with fidelity `f` below 1: H -> M with probability f
   (Bostick et al. 2007, Science 317:1760; Sharif et al. 2007, Nature 450:908).
3. **De novo methylation.** DNMT3A/B methylate an unmethylated strand with
   probability `d` per strand per division (Okano et al. 1999, Cell 99:247).
4. **Active erasure.** TET1-3 oxidise 5mC and it is lost with probability `e`
   per methylated strand per division (Tahiliani et al. 2009, Science 324:930;
   Ito et al. 2011, Science 333:1300).

**Neighbour dependence.** Maintenance and de novo efficiency rise with the local
CpG density: a lone CpG is a poor substrate, a cluster is a good one, which is
the published mechanism for why isolated CpGs erode while dense ones do not
(collaborative models: Haerter et al. 2014, NAR 42:2235; Lovkvist et al. 2016,
NAR 44:5123; the genomic observation is Zhou et al. 2018, Nat Genet 50:591 —
solo-WCGW CpGs, isolated and flanked by A or T, lose methylation with mitotic
age). Here, with `n` other CpGs inside the window and a density weight
`h(n) = 1 - exp(-n / neighbour_scale)`:

    maintenance(n, flank) = f_solo(flank) + (f_dense - f_solo(flank)) * h(n)
    de_novo(n)            = d_solo        + (d_dense - d_solo)        * h(n)

`h(0) = 0` (a solo CpG gets the solo values) and `h -> 1` in a dense cluster.
The flank enters only through the solo fidelity, because a flanking base cannot
matter once a neighbour is close enough to recruit the enzyme.

Two ways to run it, both exact about what they are:

- `expected_trajectory` / `steady_state` propagate the distribution over
  (U, H, M) through the per-division transition matrix. No sampling: this is
  the mean over an infinite population of identically-parameterised CpGs.
- `simulate` runs individual CpGs from a seeded `random.Random`, which is what a
  lineage of one cell does and what a small number of sites really looks like.

Every parameter is a BioIR `Parameter` with evidence and confidence, and the
ones the literature does not establish per CpG per division are `UNKNOWN`
(`PARAMETERS` below). `MethylationParams.from_module` refuses to run a program
that has not bound them: a rate nobody has measured is the caller's assumption
to declare, not the engine's to invent.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field, replace
from typing import Any

from genomeos.ir import UNKNOWN, Evidence, EvidenceKind, Module, Parameter

# ---- evidence ------------------------------------------------------------------------

E_LAIRD = Evidence(
    EvidenceKind.EXPERIMENTAL,
    "Laird et al. 2004, PNAS 101:204; Genereux et al. 2005, PNAS 102:5802",
    note="hairpin-bisulfite double-strand reads: maintenance fidelity roughly 0.95-0.99 per CpG "
    "per division, de novo 0.01-0.2, both locus- and cell-type-dependent",
)
E_USHIJIMA = Evidence(
    EvidenceKind.EXPERIMENTAL,
    "Ushijima et al. 2003, PNAS 100:2709",
    note="clonal populations put per-site fidelity above 0.99, the top of the published range",
)
E_ZHOU = Evidence(
    EvidenceKind.EXPERIMENTAL,
    "Zhou et al. 2018, Nat Genet 50:591",
    note="solo-WCGW CpGs (no other CpG within 35 bp, A/T on both sides) lose methylation with "
    "mitotic age in partially methylated domains; no per-division fidelity is reported",
)
E_ENDICOTT = Evidence(
    EvidenceKind.EXPERIMENTAL,
    "Endicott et al. 2022, Nat Commun 13:6659",
    note="cell division drives the loss in primary cells in culture; the ordering is by divisions, "
    "the per-CpG rate is not resolved",
)
E_UHRF1 = Evidence(
    EvidenceKind.EXPERIMENTAL,
    "Bostick et al. 2007, Science 317:1760; Sharif et al. 2007, Nature 450:908",
    note="UHRF1 binds the hemimethylated dyad and recruits DNMT1; maintenance is a read-write step",
)
E_DNMT3 = Evidence(
    EvidenceKind.EXPERIMENTAL,
    "Okano et al. 1999, Cell 99:247",
    note="DNMT3A/B are the de novo methyltransferases; the per-division rate per context is not "
    "established outside a few loci",
)
E_TET = Evidence(
    EvidenceKind.EXPERIMENTAL,
    "Tahiliani et al. 2009, Science 324:930; Ito et al. 2011, Science 333:1300",
    note="TET1-3 oxidise 5mC; genome-average active turnover per division in a somatic cell is not "
    "established",
)
E_COLLABORATIVE = Evidence(
    EvidenceKind.CURATED,
    "Haerter et al. 2014, NAR 42:2235; Lovkvist et al. 2016, NAR 44:5123",
    note="collaborative (neighbour-dependent) methylation: recruitment by nearby methylated CpGs. "
    "The models are published; the CpG count at which the effect saturates is not",
)

#: What the engine needs, with evidence and confidence; UNKNOWN where the literature does not
#: establish a per-CpG-per-division number. A program binds the UNKNOWN ones as its own
#: assumptions (see `genomeos/std/methylation.bio` and `data/demo/methylation_erosion.bio`).
PARAMETERS: dict[str, Parameter] = {
    "maintenance_fidelity_dense": Parameter(
        "maintenance_fidelity_dense", 0.97, "per CpG per division", E_LAIRD, 0.5
    ),
    "maintenance_fidelity_solo_wcgw": Parameter(
        "maintenance_fidelity_solo_wcgw", UNKNOWN, "per CpG per division", E_ZHOU, 0.0
    ),
    "maintenance_fidelity_solo_scgs": Parameter(
        "maintenance_fidelity_solo_scgs", UNKNOWN, "per CpG per division", E_ZHOU, 0.0
    ),
    "de_novo_dense": Parameter("de_novo_dense", UNKNOWN, "per strand per division", E_DNMT3, 0.0),
    "de_novo_solo": Parameter("de_novo_solo", UNKNOWN, "per strand per division", E_DNMT3, 0.0),
    "tet_erasure": Parameter("tet_erasure", UNKNOWN, "per strand per division", E_TET, 0.0),
    "neighbour_scale": Parameter("neighbour_scale", UNKNOWN, "CpGs", E_COLLABORATIVE, 0.0),
    "neighbour_window_bp": Parameter("neighbour_window_bp", 35.0, "bp", E_ZHOU, 0.9),
}

#: The prefix a BioLang program uses: `param methylation.tet_erasure = 0.0 { ... }`.
PARAM_PREFIX = "methylation."

# ---- state ---------------------------------------------------------------------------

U, H, M = 0, 1, 2
STATE_NAMES = ("unmethylated", "hemimethylated", "methylated")


@dataclass(frozen=True, slots=True)
class CpGContext:
    """A class of CpG site: how many other CpGs sit inside the window, and the flanking bases.

    `flank` is "W" for A/T on both sides (the W of solo-WCGW) and "S" for a C/G on either side.
    """

    name: str
    neighbours: int
    flank: str = "W"

    def __post_init__(self) -> None:
        if self.neighbours < 0:
            raise ValueError(f"{self.name}: neighbours must be 0 or more")
        if self.flank not in ("W", "S"):
            raise ValueError(f"{self.name}: flank must be W (A/T on both sides) or S")


#: The contexts the published observation is about: an isolated CpG in A/T sequence, an isolated
#: one with a G/C flank, a CpG with a couple of neighbours, and one inside a cluster.
DEFAULT_CONTEXTS: tuple[CpGContext, ...] = (
    CpGContext("solo_wcgw", 0, "W"),
    CpGContext("solo_scgs", 0, "S"),
    CpGContext("sparse", 2, "W"),
    CpGContext("dense", 8, "S"),
)


class UnknownParametersError(Exception):
    """A run was asked for with parameters the literature does not establish and nobody bound."""

    def __init__(self, names: list[str]) -> None:
        self.names = names
        super().__init__(
            "methylation parameters are UNKNOWN and were not bound by the program: "
            + ", ".join(names)
            + ". Declare them (with evidence and a confidence) or pass them to the engine; "
            "the engine will not invent a rate."
        )


@dataclass(frozen=True, slots=True)
class MethylationParams:
    """The mechanism's rates. Every field is per CpG (or per strand) per cell division."""

    maintenance_fidelity_dense: float
    maintenance_fidelity_solo_wcgw: float
    maintenance_fidelity_solo_scgs: float
    de_novo_dense: float
    de_novo_solo: float
    tet_erasure: float
    neighbour_scale: float
    neighbour_window_bp: float = 35.0
    assumptions: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        for name in (
            "maintenance_fidelity_dense",
            "maintenance_fidelity_solo_wcgw",
            "maintenance_fidelity_solo_scgs",
            "de_novo_dense",
            "de_novo_solo",
            "tet_erasure",
        ):
            v = getattr(self, name)
            if not 0.0 <= float(v) <= 1.0:
                raise ValueError(f"{name} is a probability per division: {v} is outside 0..1")
        if float(self.neighbour_scale) <= 0.0:
            raise ValueError("neighbour_scale must be positive (CpGs)")

    # -- the neighbour dependence ------------------------------------------------------

    def density_weight(self, neighbours: int) -> float:
        """0 for a solo CpG, rising to 1 in a cluster."""
        return 1.0 - math.exp(-neighbours / self.neighbour_scale)

    def solo_fidelity(self, flank: str) -> float:
        return self.maintenance_fidelity_solo_wcgw if flank == "W" else self.maintenance_fidelity_solo_scgs

    def maintenance(self, ctx: CpGContext) -> float:
        solo = self.solo_fidelity(ctx.flank)
        return solo + (self.maintenance_fidelity_dense - solo) * self.density_weight(ctx.neighbours)

    def de_novo(self, ctx: CpGContext) -> float:
        return self.de_novo_solo + (self.de_novo_dense - self.de_novo_solo) * self.density_weight(
            ctx.neighbours
        )

    # -- provenance --------------------------------------------------------------------

    @classmethod
    def from_module(cls, module: Module, assume: dict[str, float] | None = None) -> MethylationParams:
        """Read `param methylation.*` from a compiled program; refuse what is still UNKNOWN."""
        given = dict(assume or {})
        values: dict[str, float] = {}
        unknown: list[str] = []
        notes: list[str] = []
        for name, default in PARAMETERS.items():
            p = module.parameters.get(PARAM_PREFIX + name)
            if name in given:
                values[name] = float(given[name])
                notes.append(f"{name} = {values[name]:g} passed to the engine by the caller")
                continue
            if p is None or p.value is UNKNOWN:
                if default.value is UNKNOWN:
                    unknown.append(name)
                else:
                    values[name] = float(default.value)
                    notes.append(
                        f"{name} = {values[name]:g} from the engine's literature value "
                        f"({default.evidence.source})"
                    )
                continue
            values[name] = float(p.value)
            if p.confidence < 0.4 or p.evidence.kind is EvidenceKind.INFERRED:
                notes.append(
                    f"{name} = {values[name]:g} is the program's assumption "
                    f"(confidence {p.confidence:g}: {p.evidence.source or 'no source'})"
                )
        if unknown:
            raise UnknownParametersError(sorted(unknown))
        return cls(**values, assumptions=tuple(notes))

    def parameters(self) -> list[Parameter]:
        """The values actually used, each carrying the evidence of the fact it stands for."""
        out = []
        for name, default in PARAMETERS.items():
            out.append(
                Parameter(
                    name,
                    float(getattr(self, name)),
                    default.unit,
                    default.evidence,
                    default.confidence if default.value is not UNKNOWN else 0.2,
                )
            )
        return out


# ---- one division, as a distribution over (U, H, M) ----------------------------------


def _replication(dist: tuple[float, float, float]) -> tuple[float, float, float]:
    u, h, m = dist
    return (u + 0.5 * h, 0.5 * h + m, 0.0)


def _maintenance(dist: tuple[float, float, float], f: float) -> tuple[float, float, float]:
    u, h, m = dist
    return (u, h * (1.0 - f), m + h * f)


def _de_novo(dist: tuple[float, float, float], d: float) -> tuple[float, float, float]:
    u, h, m = dist
    return (u * (1.0 - d) ** 2, u * 2.0 * d * (1.0 - d) + h * (1.0 - d), u * d * d + h * d + m)


def _erasure(dist: tuple[float, float, float], e: float) -> tuple[float, float, float]:
    u, h, m = dist
    return (u + h * e + m * e * e, h * (1.0 - e) + m * 2.0 * e * (1.0 - e), m * (1.0 - e) ** 2)


def division(
    dist: tuple[float, float, float], ctx: CpGContext, params: MethylationParams
) -> tuple[float, float, float]:
    """One cell division: replication, maintenance, de novo, erasure, in that order."""
    d = _replication(dist)
    d = _maintenance(d, params.maintenance(ctx))
    d = _de_novo(d, params.de_novo(ctx))
    return _erasure(d, params.tet_erasure)


def level(dist: tuple[float, float, float]) -> float:
    """The methylated share of strands, which is what bisulfite sequencing measures."""
    return dist[M] + 0.5 * dist[H]


def initial_distribution(initial_level: float, hemi: float = 0.0) -> tuple[float, float, float]:
    """A starting distribution with this methylation level: M and U, plus an optional H share."""
    if not 0.0 <= initial_level <= 1.0:
        raise ValueError(f"initial_level must be within 0..1, not {initial_level}")
    if not 0.0 <= hemi <= 1.0:
        raise ValueError(f"hemi must be within 0..1, not {hemi}")
    m = initial_level - 0.5 * hemi
    if m < -1e-12 or m > 1.0 - hemi + 1e-12:
        raise ValueError(f"level {initial_level} is not reachable with a hemimethylated share of {hemi}")
    m = min(max(m, 0.0), 1.0 - hemi)
    return (1.0 - hemi - m, hemi, m)


def expected_trajectory(
    ctx: CpGContext,
    params: MethylationParams,
    divisions: int,
    initial_level: float = 1.0,
    hemi: float = 0.0,
) -> list[float]:
    """The exact mean methylation level after 0, 1, ... `divisions` divisions (no sampling)."""
    if divisions < 0:
        raise ValueError("divisions must be 0 or more")
    dist = initial_distribution(initial_level, hemi)
    out = [level(dist)]
    for _ in range(divisions):
        dist = division(dist, ctx, params)
        out.append(level(dist))
    return out


def steady_state(
    ctx: CpGContext,
    params: MethylationParams,
    tol: float = 1e-12,
    max_divisions: int = 200_000,
) -> tuple[float, int]:
    """(level, divisions to reach it) of the fixed point this context converges to."""
    dist = initial_distribution(1.0)
    prev = level(dist)
    for i in range(1, max_divisions + 1):
        dist = division(dist, ctx, params)
        cur = level(dist)
        if abs(cur - prev) < tol:
            return cur, i
        prev = cur
    return prev, max_divisions


def fidelity_for_level(
    target_level: float,
    divisions: int,
    params: MethylationParams,
    ctx: CpGContext,
    knob: str = "maintenance_fidelity_dense",
    initial_level: float = 1.0,
    tol: float = 1e-10,
) -> float | None:
    """The maintenance fidelity that reproduces `target_level` after `divisions` divisions.

    `knob` is the field to solve for (dense, solo-WCGW or solo-SCGS fidelity). The level after a
    fixed number of divisions rises monotonically with fidelity, so this is a bisection. Returns
    None when the target is outside what any fidelity in 0..1 can reach with these other rates:
    an unreachable observation is a refutation of the parameter set, not a number to report.
    """
    knob = {"solo": "maintenance_fidelity_solo_wcgw", "dense": "maintenance_fidelity_dense"}.get(knob, knob)
    if knob not in (
        "maintenance_fidelity_dense",
        "maintenance_fidelity_solo_wcgw",
        "maintenance_fidelity_solo_scgs",
    ):
        raise ValueError(f"cannot solve for {knob!r}: only a maintenance fidelity")

    def reached(f: float) -> float:
        return expected_trajectory(ctx, replace(params, **{knob: f}), divisions, initial_level=initial_level)[
            -1
        ]

    lo, hi = reached(0.0), reached(1.0)
    if not lo - tol <= target_level <= hi + tol:
        return None
    a, b = 0.0, 1.0
    for _ in range(200):
        mid = 0.5 * (a + b)
        if reached(mid) < target_level:
            a = mid
        else:
            b = mid
        if b - a < tol:
            break
    return 0.5 * (a + b)


def divisions_for_level(
    target_level: float,
    params: MethylationParams,
    ctx: CpGContext,
    initial_level: float = 1.0,
    max_divisions: int = 100_000,
) -> int | None:
    """The first division count at which the expected level reaches `target_level` from above.

    None when the trajectory never gets there (its steady state is above the target).
    """
    dist = initial_distribution(initial_level)
    if level(dist) <= target_level:
        return 0
    for i in range(1, max_divisions + 1):
        dist = division(dist, ctx, params)
        if level(dist) <= target_level:
            return i
    return None


# ---- the stochastic run --------------------------------------------------------------


@dataclass(slots=True)
class MethylationRun:
    """What the engine returns: a trajectory per CpG context, and the steady state of each.

    `levels` are the sampled means over `cpgs` individual sites; `expected` are the exact means
    of the same Markov chain. `levels` keys are prefixed `methylation.` so the toolchain can
    name them in a `# test:` line, and `methylation.gap` is the contrast the prediction is about.
    """

    contexts: tuple[CpGContext, ...]
    divisions: int
    cpgs: int
    seed: int
    params: MethylationParams
    times: list[float] = field(default_factory=list)
    levels: dict[str, list[float]] = field(default_factory=dict)
    expected: dict[str, list[float]] = field(default_factory=dict)
    steady: dict[str, float] = field(default_factory=dict)
    steady_divisions: dict[str, int] = field(default_factory=dict)
    final_states: dict[str, dict[str, int]] = field(default_factory=dict)
    gap_contexts: tuple[str, str] = ("dense", "solo_wcgw")

    @property
    def species(self) -> list[str]:
        return list(self.levels)

    def peaks(self, name: str) -> int:
        xs = self.levels.get(name) or []
        return sum(1 for i in range(1, len(xs) - 1) if xs[i] > xs[i - 1] and xs[i] >= xs[i + 1])

    def summary(self) -> dict[str, Any]:
        return {
            "divisions": self.divisions,
            "cpgs_per_context": self.cpgs,
            "seed": self.seed,
            "levels": {k: round(v[-1], 4) for k, v in self.levels.items()},
            "expected": {k: round(v[-1], 4) for k, v in self.expected.items()},
            "steady_state": {k: round(v, 4) for k, v in self.steady.items()},
            "divisions_to_steady_state": dict(self.steady_divisions),
            "maintenance_fidelity": {c.name: round(self.params.maintenance(c), 4) for c in self.contexts},
            "de_novo": {c.name: round(self.params.de_novo(c), 4) for c in self.contexts},
            "assumptions": list(self.params.assumptions),
        }


def simulate(
    params: MethylationParams,
    divisions: int = 100,
    contexts: tuple[CpGContext, ...] = DEFAULT_CONTEXTS,
    cpgs: int = 500,
    seed: int = 0,
    initial_level: float = 1.0,
) -> MethylationRun:
    """Run `cpgs` sites of each context through `divisions` divisions of one lineage.

    Deterministic under `seed`: the same seed gives the same trajectory, site for site.
    """
    if divisions < 0:
        raise ValueError("divisions must be 0 or more")
    if cpgs < 1:
        raise ValueError("cpgs must be 1 or more")
    run = MethylationRun(
        contexts=tuple(contexts),
        divisions=divisions,
        cpgs=cpgs,
        seed=seed,
        params=params,
        times=[float(i) for i in range(divisions + 1)],
    )
    n_methylated = int(round(initial_level * cpgs))
    for ctx in contexts:
        rng = random.Random(f"{seed}:{ctx.name}")
        f, d, e = params.maintenance(ctx), params.de_novo(ctx), params.tet_erasure
        states = [M] * n_methylated + [U] * (cpgs - n_methylated)
        series = [sum(2 if s == M else s == H for s in states) / (2.0 * cpgs)]
        for _ in range(divisions):
            for i, s in enumerate(states):
                if s == M:  # replication leaves one methylated strand
                    s = H
                elif s == H:  # the methylated strand goes to one daughter only
                    s = H if rng.random() < 0.5 else U
                if s == H and rng.random() < f:  # UHRF1 + DNMT1
                    s = M
                if s == U:  # DNMT3A/B on either strand
                    hits = (rng.random() < d) + (rng.random() < d)
                    s = (U, H, M)[hits]
                elif s == H and rng.random() < d:
                    s = M
                if e > 0.0:  # TET
                    if s == M:
                        hits = (rng.random() < e) + (rng.random() < e)
                        s = (M, H, U)[hits]
                    elif s == H and rng.random() < e:
                        s = U
                states[i] = s
            series.append(sum(2 if s == M else s == H for s in states) / (2.0 * cpgs))
        key = PARAM_PREFIX + ctx.name
        run.levels[key] = series
        run.expected[key] = expected_trajectory(ctx, params, divisions, initial_level=initial_level)
        st, n = steady_state(ctx, params)
        run.steady[key] = st
        run.steady_divisions[key] = n
        run.final_states[key] = {
            STATE_NAMES[U]: states.count(U),
            STATE_NAMES[H]: states.count(H),
            STATE_NAMES[M]: states.count(M),
        }
    names = {c.name for c in contexts}
    a, b = run.gap_contexts
    if a in names and b in names:
        run.levels[PARAM_PREFIX + "gap"] = [
            x - y for x, y in zip(run.levels[PARAM_PREFIX + a], run.levels[PARAM_PREFIX + b], strict=True)
        ]
        run.expected[PARAM_PREFIX + "gap"] = [
            x - y for x, y in zip(run.expected[PARAM_PREFIX + a], run.expected[PARAM_PREFIX + b], strict=True)
        ]
        run.steady[PARAM_PREFIX + "gap"] = run.steady[PARAM_PREFIX + a] - run.steady[PARAM_PREFIX + b]
        run.steady_divisions[PARAM_PREFIX + "gap"] = max(
            run.steady_divisions[PARAM_PREFIX + a], run.steady_divisions[PARAM_PREFIX + b]
        )
    return run


def run_module(
    module: Module,
    divisions: int | None = None,
    cpgs: int | None = None,
    seed: int | None = None,
    assume: dict[str, float] | None = None,
) -> MethylationRun:
    """Run the methylation engine on a compiled program's `param methylation.*` block.

    A program declares the rates and, optionally, `methylation.divisions`,
    `methylation.cpgs`, `methylation.seed` and `methylation.initial_level`.
    """
    params = MethylationParams.from_module(module, assume=assume)

    def setting(name: str, default: float) -> float:
        p = module.parameters.get(PARAM_PREFIX + name)
        return default if p is None or p.value is UNKNOWN else float(p.value)

    return simulate(
        params,
        divisions=int(divisions if divisions is not None else setting("divisions", 100)),
        cpgs=int(cpgs if cpgs is not None else setting("cpgs", 500)),
        seed=int(seed if seed is not None else setting("seed", 0)),
        initial_level=setting("initial_level", 1.0),
    )


def declares_methylation(module: Module) -> bool:
    """Whether a program asks for this engine (it declares any `methylation.*` parameter)."""
    return any(name.startswith(PARAM_PREFIX) for name in module.parameters)
