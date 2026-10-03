# SPDX-License-Identifier: Apache-2.0
# Part of the BioLang engine (language, IR, VM, standard library); see LICENSING.md.
"""Gene regulatory network runtime.

Executes a BioIR Module as a continuous-time system: mRNA and protein levels
change according to PRODUCE / ACTIVATE / INHIBIT rules. Regulation uses Hill
functions, the standard phenomenological form for transcription-factor binding.

    transcription(gene) = basal + max_rate * A(gene) * R(gene)

    A = mean over activators of  s_i * H(x_i)         (1 if no activators)
    R = product over inhibitors of (1 - s_j * H(x_j))
    H(x) = x^n / (K^n + x^n)

    d[mRNA]/dt    = transcription - delta_m * [mRNA]
    d[protein]/dt = k_tl * [mRNA] - delta_p * [protein]

Integration: RK4 when noise == 0, Euler-Maruyama otherwise. Everything is in
arbitrary concentration units and hours; parameters in the module override the
defaults below and carry their own evidence.

The combination rules are model assumptions, named in `ACTIVATOR_COMBINATION` and
`INHIBITOR_COMBINATION`. Under the mean, adding an activator that is low where
another is high lowers that gene's drive there (lane-sign, 2026-09-28).
`combine_activators` computes the mean and the three alternatives the census in
docs/DESIGN-MINIMAL-CELL.md compares; the runtime uses the one `ACTIVATOR_COMBINATION` names.
The registered choice between them is `COMBINATION_TEST` (2026-09-28).

Unresolved inputs (review R3, 2026-09-28). A regulator the runtime holds no state
for, and a regulated gene with no `max_rate`, are still integrated as zero, as they
always were, but every `run` now reports them by name (`NetworkRuntime.diagnose`):
on the trajectory's `unresolved` list and as an `UnresolvedModelWarning`, or as an
`UnresolvedModel` error when the runtime is strict. A source whose name ends in
`EXPLICIT_ZERO_SUFFIX` is a zero the caller declared (network_experiment's edge
knockout) and is not unresolved. A state is given by `clamp`; `initial` does not
hold a species the runtime does not integrate.

Units (review R3, lane-rates, 2026-09-28). `DEFAULTS` below are arbitrary units per hour and are
not measurements of anything. A module whose parameters carry measured values puts the whole run
into real units: `mrna_half_life` and a protein's `half_life_h` in HOURS give delta_m and delta_p
in 1/h, `translation_rate` in molecules per mRNA per hour, and a gene's `basal_rate` and `max_rate`
in molecules per cell per hour, so that an mRNA level reads as molecules per cell. A measured
half-life is a degradation constant and never a transcription rate; a measured transcription rate is
the gene's TOTAL rate at the state it was measured in and never `max_rate`, which is a ceiling no
steady-state measurement observes. `attribution/bridge.py` holds the full mapping, the sources and
what each of them does not license (`MEASURED_RATE_UNITS`, `split_measured_rate`). Mixing a measured
parameter with a defaulted one mixes units: the defaults are left in place only so that a module
that declares nothing still runs, and such a run is in a.u.

Provenance on the output (item 12 S5, 2026-09-28). A number measured somewhere other than the cell
being simulated travels on the trajectory, never only in a header: a gene carrying the attributes in
`RATE_PROVENANCE_ATTRS` (the application's bridge writes them when it puts a measured or borrowed
rate on a gene) is listed on `Trajectory.transferred` with its source, the species and cell it was
measured in, its tier and the cell context of this run, and every kinetic parameter the run took
from `DEFAULTS` because the module declared none is listed on `Trajectory.defaulted`. A level read
from a run with a transferred rate and a defaulted half-life is in mixed units, and the two lists
together say so. The runtime labels; it never judges whether a number is validated.
"""

from __future__ import annotations

import math
import random
import warnings
from dataclasses import dataclass, field

from genomeos.ir import UNKNOWN, Action, Gene, Module, Protein, Rule

LN2 = math.log(2)

DEFAULTS = {
    "translation_rate": 5.0,  # protein per mRNA per hour
    "mrna_half_life": LN2,  # hours  (delta_m = 1 / hour)
    "protein_half_life": LN2 / 5,  # hours (delta_p = 5 / hour) unless protein declares one
    "noise": 0.0,  # multiplicative noise amplitude
}

#: how a gene's activators and inhibitors combine: declared model assumptions, not findings
ACTIVATOR_COMBINATION = "mean"  # A = mean of s_i * H(x_i): a low extra activator dilutes the drive
INHIBITOR_COMBINATION = "product"  # R = product of (1 - s_j * H(x_j))
#: the activator rules the combination registration (lane-combine, 2026-09-28) compares; each is
#: over the terms t_i = s_i * H(x_i) in [0, 1], and all four agree when a gene has one activator
ACTIVATOR_RULES = ("mean", "sum_capped", "max", "or")
#: The choice test between the rules, registered 2026-09-28 before its one run (lane-combine;
#: docs/DESIGN-MINIMAL-CELL.md, "The activator combination rule: a choice test"). Each case is a
#: published or project measurement whose outcome a rule either can or cannot produce whatever its
#: numbers are; a rule is contradicted by a case only in that parameter-free sense. No strength,
#: threshold, Hill coefficient or basal rate is re-tuned for any rule, and the gastrulation census is
#: reported under each rule but is not a criterion (choosing a rule by that fit is a tuning).
COMBINATION_TEST = {
    "cases": {
        # Bothma et al. 2015 eLife 4:e07956 (PMC4532966): one reporter, primary / shadow / both
        # enhancers, the pair P against the singles S1, S2 above basal. Predictions: mean P = (S1+S2)/2;
        # max P = max(S1, S2); sum_capped max <= P <= S1+S2; or max <= P <= S1+S2.
        "bothma2015_kni_early": {"observed": "super-additive, P > S1 + S2", "contradicts": ACTIVATOR_RULES},
        "bothma2015_kni_late": {"observed": "additive, P = S1 + S2", "contradicts": ("mean", "max")},
        "bothma2015_hb_anterior": {
            "observed": "sub-additive; removing the shadow enhancer has no effect",
            "contradicts": (),
        },
        "bothma2015_hb_central": {"observed": "additive, P = S1 + S2", "contradicts": ("mean", "max")},
        "bothma2015_sna": {
            "observed": "P below the shadow enhancer alone, P < max(S1, S2)",
            "contradicts": ("sum_capped", "max", "or"),
        },
        # the project's CRISPRi training split: single removals only; a gene-cell with two significant
        # decreases is impossible under max (removing a non-maximal activator changes nothing)
        "crispri_training_two_decreases": {"observed": "computed at the run", "contradicts": "max if >= 10"},
    },
    "not_judged": {
        "zhou2024_gasperini_doubles": "compares additive and multiplicative GLM links for enhancer pairs; "
        "each rule's double-removal prediction depends on the gene's other activators and saturation",
        "crispri_training_mean": "under the mean the single-removal changes of a gene's activators sum to "
        "zero, but untested or undetected activators can carry the balancing increases",
        "gastrulation_tbxt_knockdown": "Lolas 2014 Fig. 3A (Brachyury depletion lowers Sox17): whether "
        "the model reproduces it depends on the module's unsourced strengths and thresholds",
        "gastrulation_census": "a fit criterion, not a mechanism test",
    },
    # a case that contradicts all four rules says the family is inadequate and discriminates nothing
    "exclude_cases_contradicting_every_rule": True,
    # switch from the mean to R only if the mean is contradicted by >= 1 judged case and R by none; if
    # no R or several survive, the mean is retained and the choice is recorded as undetermined
    "pass_rule": "switch to the single surviving rule, else mean retained, undetermined",
    "crispri_min_gene_cells_for_max": 10,
}
#: a rule source ending in this reads zero by the caller's declaration, never by omission
EXPLICIT_ZERO_SUFFIX = "@zero"


#: gene attributes naming where a gene's rates were measured (item 12 S5): source key, species and
#: cell, and tier (measured | borrowed_median). Read back onto every trajectory; a missing one is
#: reported as `unstated`, never filled in with a species
RATE_PROVENANCE_ATTRS = ("rate_source", "rate_species_cell", "rate_tier")
#: the quantity a gene's `basal_rate` and `max_rate` are split from when they carry that provenance
TRANSFERRED_RATE_QUANTITY = "total transcription rate T (split into basal_rate and max_rate)"
UNSTATED_PROVENANCE = "unstated"


@dataclass(frozen=True, slots=True)
class Transferred:
    """One number a run used that was measured somewhere other than the cell it simulates (S5)."""

    quantity: str
    subject: str  # the gene it belongs to
    source: str
    species_cell: str  # where it was measured
    tier: str  # measured | borrowed_median | unstated
    used_in: str  # the run's cell_type context, "" when the run names none


@dataclass(frozen=True, slots=True)
class Unresolved:
    """One input the model needs and was not given: why, for what, and how to give it."""

    reason: str  # regulator_state_missing | gene_parameter_missing | (bridge reasons)
    subject: str
    detail: str = ""


class UnresolvedModelError(ValueError):
    """A strict run refused because the model is not resolved; `issues` names every missing input."""

    def __init__(self, issues: list[Unresolved]) -> None:
        self.issues = list(issues)
        super().__init__(
            "unresolved model: "
            + "; ".join(
                f"{i.reason} {i.subject}" + (f" ({i.detail})" if i.detail else "") for i in self.issues
            )
        )


#: the name the R3 registration used
UnresolvedModel = UnresolvedModelError


class UnresolvedModelWarning(UserWarning):
    """A non-strict run integrated missing inputs as zero; the trajectory lists them."""


@dataclass(slots=True)
class Trajectory:
    times: list[float]
    species: list[str]
    levels: dict[str, list[float]] = field(default_factory=dict)
    unresolved: list[Unresolved] = field(default_factory=list)  # what the run read as zero undeclared
    transferred: list[Transferred] = field(default_factory=list)  # numbers measured elsewhere (S5)
    defaulted: list[str] = field(default_factory=list)  # kinetic parameters left at a.u. DEFAULTS

    def provenance(self) -> dict:
        """The run's parameter provenance as plain data, for any summary that quotes a level from it."""
        return {
            "transferred": [
                {
                    "quantity": t.quantity,
                    "subject": t.subject,
                    "source": t.source,
                    "species_cell": t.species_cell,
                    "tier": t.tier,
                    "used_in": t.used_in,
                }
                for t in self.transferred
            ],
            "defaulted": list(self.defaulted),
        }

    def final(self) -> dict[str, float]:
        return {s: self.levels[s][-1] for s in self.species}

    def peaks(self, species: str) -> int:
        """Count local maxima; a crude oscillation detector."""
        xs = self.levels[species]
        return sum(1 for i in range(1, len(xs) - 1) if xs[i - 1] < xs[i] >= xs[i + 1])


def _hill(x: float, k: float, n: float) -> float:
    if x <= 0:
        return 0.0
    xn = x**n
    return xn / (k**n + xn)


def combine_activators(terms: list[float], rule: str | None = None) -> float:
    """A gene's activator drive A from its terms t_i = s_i * H(x_i), under `rule`.

    mean: sum(t) / n (the runtime's rule since the first commit); sum_capped: min(1, sum(t));
    max: max(t); or: 1 - prod(1 - t_i), independent binding of any one activator sufficing.
    `rule` defaults to `ACTIVATOR_COMBINATION`, read at call time.
    """
    rule = ACTIVATOR_COMBINATION if rule is None else rule
    if not terms:
        return 1.0
    if rule == "mean":
        return sum(terms) / len(terms)
    if rule == "sum_capped":
        return min(1.0, sum(terms))
    if rule == "max":
        return max(terms)
    if rule == "or":
        miss = 1.0
        for t in terms:
            miss *= 1.0 - t
        return 1.0 - miss
    raise ValueError(f"unknown activator combination {rule!r}; one of {', '.join(ACTIVATOR_RULES)}")


class NetworkRuntime:
    def __init__(
        self,
        module: Module,
        context: dict[str, str] | None = None,
        seed: int | None = None,
        strict: bool = False,
    ) -> None:
        self.module = module
        self.strict = strict  # refuse to run an unresolved model instead of warning (R3)
        self.context = context or {}
        self.rng = random.Random(seed)
        self.params = dict(DEFAULTS)
        for name, p in module.parameters.items():
            if name in self.params:
                self.params[name] = p.value

        self.genes: list[Gene] = module.genes()
        self.proteins: list[Protein] = module.proteins()
        # `Module.active_rules` rather than a bare `applies` filter: a rule whose source or target is
        # not expressed in this cell type is not a rule this cell runs, and docs/BIOLANG-v0.2.md
        # documents the contextual gate as both the `when` clause AND the `expresses` list. Until
        # 2026-09-21 this used only the first, so `genomeos run --context cell_type=X` overstated what
        # a cell expresses.
        self.active_rules: list[Rule] = module.active_rules(self.context)
        # and dropping the rules is NOT enough, which is why the silenced set is kept as well. A gene
        # with no activator takes `a = 1.0` below and transcribes at its FULL max_rate, so removing a
        # silenced gene's activator would have driven it to maximum instead of to zero — the opposite
        # of silencing, and a plausible-looking curve either way. Measured on cell_context.bio:
        # SYN1 read 80.05 before, 80.05 with the rules dropped, and 0 only when the gene itself is
        # held off.
        self.silenced: set[str] = module.silenced_genes(self.context)
        # this runtime has no places and therefore no volumes; a concentration threshold is only a
        # number once a compartment divides it, so it belongs to a located program (v0.4 §4.1)
        molar = [r.id for r in self.active_rules if r.threshold_unit]
        if molar:
            raise ValueError(
                f"rules state their thresholds as concentrations ({', '.join(molar[:3])}), which needs a "
                "compartment with an absolute_volume: run a located program, or state amounts"
            )

        # index regulators per gene, translation per protein
        self.activators: dict[str, list[Rule]] = {g.id: [] for g in self.genes}
        self.inhibitors: dict[str, list[Rule]] = {g.id: [] for g in self.genes}
        self.produces: dict[str, list[str]] = {p.id: [] for p in self.proteins}  # protein <- genes
        for r in self.active_rules:
            if r.action is Action.ACTIVATE and r.target in self.activators:
                self.activators[r.target].append(r)
            elif r.action is Action.INHIBIT and r.target in self.inhibitors:
                self.inhibitors[r.target].append(r)
            elif r.action is Action.PRODUCE and r.target in self.produces:
                self.produces[r.target].append(r.source)

        self.species = [f"{g.id}.mRNA" for g in self.genes] + [p.id for p in self.proteins]
        self.transferred, self.defaulted = self._provenance()

    def _provenance(self) -> tuple[list[Transferred], list[str]]:
        """Which numbers this run takes from elsewhere, and which from DEFAULTS (item 12 S5)."""
        cell = str(self.context.get("cell_type", ""))
        transferred = []
        for g in self.genes:
            if g.id in self.silenced or not any(k in g.attrs for k in RATE_PROVENANCE_ATTRS):
                continue
            source, species_cell, tier = (
                str(g.attrs.get(k, UNSTATED_PROVENANCE)) for k in RATE_PROVENANCE_ATTRS
            )
            transferred.append(Transferred(TRANSFERRED_RATE_QUANTITY, g.id, source, species_cell, tier, cell))
        declared = set(self.module.parameters)
        used = []
        if self.genes:
            used.append("mrna_half_life")
        if any(self.produces.values()):
            used.append("translation_rate")
        if any(p.half_life_h is UNKNOWN for p in self.proteins):
            used.append("protein_half_life")
        return transferred, [name for name in used if name not in declared]

    def diagnose(self, clamp: dict[str, float] | None = None) -> list[Unresolved]:
        """Every input a run with this clamp would read as zero without anyone having said so (R3).

        A regulator of an expressed gene must be a species, be clamped, or carry
        `EXPLICIT_ZERO_SUFFIX`; a gene an active rule regulates must declare `max_rate`.
        """
        held = set(clamp or {})
        species = set(self.species)
        genes = {g.id: g for g in self.genes}
        out: list[Unresolved] = []
        seen: set[tuple[str, str]] = set()
        for r in self.active_rules:
            if r.action not in (Action.ACTIVATE, Action.INHIBIT) or r.target not in genes:
                continue
            if r.target in self.silenced:
                continue
            src = r.source
            missing = src not in species and src not in held and not src.endswith(EXPLICIT_ZERO_SUFFIX)
            if missing and ("regulator_state_missing", src) not in seen:
                seen.add(("regulator_state_missing", src))
                out.append(
                    Unresolved(
                        "regulator_state_missing",
                        src,
                        f"regulates {r.target} (rule {r.id}); not a species and not clamped, read as 0",
                    )
                )
            if "max_rate" not in genes[r.target].attrs and ("gene_parameter_missing", r.target) not in seen:
                seen.add(("gene_parameter_missing", r.target))
                out.append(
                    Unresolved(
                        "gene_parameter_missing", r.target, "regulated but declares no max_rate, read as 0"
                    )
                )
        return out

    # ---- dynamics --------------------------------------------------------

    def _derivatives(self, state: dict[str, float]) -> dict[str, float]:
        d: dict[str, float] = {}
        delta_m = LN2 / self.params["mrna_half_life"]
        k_tl = self.params["translation_rate"]

        for g in self.genes:
            acts = self.activators[g.id]
            inhs = self.inhibitors[g.id]
            a = 1.0
            if acts:
                a = combine_activators(
                    [r.strength * _hill(state.get(r.source, 0.0), r.threshold, r.hill) for r in acts]
                )
            rep = 1.0
            for r in inhs:
                rep *= 1.0 - r.strength * _hill(state.get(r.source, 0.0), r.threshold, r.hill)
            max_rate = float(g.attrs.get("max_rate", 0.0))
            key = f"{g.id}.mRNA"
            if g.id in self.silenced:
                d[key] = -delta_m * state[key]  # not expressed here: it decays and is not made
                continue
            d[key] = g.basal_rate + max_rate * a * rep - delta_m * state[key]

        for p in self.proteins:
            hl = p.half_life_h if p.half_life_h is not UNKNOWN else self.params["protein_half_life"]
            delta_p = LN2 / float(hl)
            production = sum(k_tl * state[f"{gid}.mRNA"] for gid in self.produces[p.id])
            d[p.id] = production - delta_p * state[p.id]
        return d

    def run(
        self,
        hours: float,
        dt: float = 0.01,
        initial: dict[str, float] | None = None,
        record_every: int = 10,
        clamp: dict[str, float] | None = None,
    ) -> Trajectory:
        """Integrate for `hours`. `clamp` holds species at fixed levels (external
        signals such as a morphogen the module does not itself produce, or an mRNA
        knocked out at zero).

        A clamped species is held inside every Runge-Kutta stage, not only between
        steps: its derivative is zero and every stage state carries the clamp value.
        Until 2026-09-28 the clamp was re-applied only after a whole step, so a
        clamped mRNA was transcribed inside the stages and its protein translated off
        them (the three-gene ring's protein A read 1.87 with its mRNA clamped at zero,
        a tenth of the unperturbed 19.2), and a clamped signal decayed inside them.
        Census and values: data/results/grn_clamp_census.json (before), _after.json."""
        held = dict(clamp or {})
        unresolved = self.diagnose(held)
        if unresolved and self.strict:
            raise UnresolvedModel(unresolved)
        if unresolved:
            warnings.warn(str(UnresolvedModel(unresolved)), UnresolvedModelWarning, stacklevel=2)

        def derivatives(st: dict[str, float]) -> dict[str, float]:
            st.update(held)  # every stage reads the clamp, including a signal outside `species`
            d = self._derivatives(st)
            for s in held:
                if s in d:
                    d[s] = 0.0
            return d

        state = dict.fromkeys(self.species, 0.0)
        if initial:
            state.update(initial)
        if clamp:
            state.update(clamp)
        noise = self.params["noise"]
        steps = int(round(hours / dt))
        traj = Trajectory(
            times=[],
            species=list(self.species),
            levels={s: [] for s in self.species},
            unresolved=unresolved,
            transferred=list(self.transferred),
            defaulted=list(self.defaulted),
        )

        def record(t: float) -> None:
            traj.times.append(t)
            for s in self.species:
                traj.levels[s].append(state[s])

        record(0.0)
        for step in range(1, steps + 1):
            if noise > 0:
                d = derivatives(state)
                for s in self.species:  # a clamped species draws its noise too, then is restored below,
                    state[s] = max(  # so the random stream is the one a seed gave before 2026-09-28
                        0.0,
                        state[s]
                        + d[s] * dt
                        + noise * math.sqrt(dt) * self.rng.gauss(0, 1) * math.sqrt(max(state[s], 1e-9)),
                    )
            else:
                k1 = derivatives(state)
                s2 = {s: state[s] + 0.5 * dt * k1[s] for s in self.species}
                k2 = derivatives(s2)
                s3 = {s: state[s] + 0.5 * dt * k2[s] for s in self.species}
                k3 = derivatives(s3)
                s4 = {s: state[s] + dt * k3[s] for s in self.species}
                k4 = derivatives(s4)
                for s in self.species:
                    state[s] = max(0.0, state[s] + dt / 6 * (k1[s] + 2 * k2[s] + 2 * k3[s] + k4[s]))
            if clamp:
                state.update(clamp)
            if step % record_every == 0:
                record(step * dt)
        if steps % record_every != 0:  # always record the final state
            record(steps * dt)
        return traj
