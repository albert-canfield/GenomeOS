# SPDX-License-Identifier: AGPL-3.0-or-later
"""Review R3 (2026-09-28): from executable annotation to a parameterised simulation, and no further.

A compiled program (`attribution/compile.py`) is **executable annotation**: it parses, its rules are
gated on their cell, `bio test` checks its counts, and a runtime will integrate it. It is not a
simulation of anything. Its rule `strength` is an observation's magnitude in that observation's own
unit (|log2 fold change| clipped at 1 for an AlphaGenome deletion, |fractional change| for a CRISPRi
screen), its target genes are stubs with no transcription parameters, and the element it names as a
source is not a species the runtime holds. `runtime/grn.py` read all three silently: a missing source
state as 0.0 and a missing `max_rate` as 0.0, so a compiled rule "ran" and moved nothing.

This module is the bridge, registered here before any of it was built. It says what one observation
may be turned into, what the model must be given before it can be turned into anything, and what is
reported when it cannot.

What it will NOT claim. A deletion log2 fold change is an observed effect of removing a piece of DNA
on one gene's steady-state expression in one cell or track. It is not a rate constant, not a binding
affinity, not a dose response and not a statement about any perturbation other than that removal.
The bridge fits one dimensionless rule strength so that the model reproduces that one observed
response, given transcription parameters and Hill assumptions the caller declares; everything else
the simulation then says is the model's, not the observation's.
"""

from __future__ import annotations

import copy
import re
from dataclasses import dataclass, field

from genomeos.ir import Action, Module, Rule
from genomeos.runtime.grn import Unresolved

REGISTERED = "2026-09-28"

# ---- the observation ------------------------------------------------------------------------------
#: what one compiled rule observed: removing element E changes gene G's steady-state expression in
#: cell C by the fold RHO = expression(E removed) / expression(E intact)
OBSERVABLE = "steady-state expression fold on removing the element: expression(removed) / expression(intact)"
#: how each evidence source's number becomes RHO, and the unit it arrives in
OBSERVATION_UNITS = {
    "predicted": (
        "log2 fold change on deletion (AlphaGenome), read from the evidence note 'effect <x> log2 fold"
        " change'; RHO = 2**x. The rule strength is |x| clipped at 1 and is used only when no note"
        " exists and it is below 1; a strength of exactly 1 is censored and unresolved"
    ),
    "experimental": (
        "fractional change in expression on CRISPRi silencing (ENCODE benchmark EffectSize), carried as"
        " the rule strength |f| with the sign from the action (activates: f < 0); RHO = 1 + f"
        ", used only when the link rests on one regulated pair (EXPERIMENTAL_SUMMARY)"
    ),
}
#: amendment of 2026-09-28, before the build (lane-assay's census A8 of measured.rule_links): a compiled
#: experimental rule's strength is the LARGEST |EffectSize| among the regulated pairs of one (element,
#: gene, cell) in the partition used, a maximum and not an observation whenever there is more than one
#: pair. The compiler writes each link's pair count on the measured element's evidence note
#: ("<gene> in <cell> from <n> pairs"); the bridge maps a link that rests on one pair and reports one
#: that rests on several as `observation_summarised`, never fitting to the maximum. No summary rule
#: (median, mean, meta-analysis) is chosen here; choosing one is its own registration.
EXPERIMENTAL_SUMMARY = "one pair maps; several pairs are unresolved (the compiled strength is their maximum)"

# ---- the model it maps into ----------------------------------------------------------------------
#: the element is a regulator whose state is its presence: dimensionless, 1 intact, 0 removed
SOURCE_STATE = "element presence, dimensionless: 1.0 intact, 0.0 removed; held by clamp, never defaulted"
INTACT, REMOVED = 1.0, 0.0
#: the target gene's transcription parameters, in the runtime's arbitrary units per hour
REQUIRED_GENE_PARAMETERS = {
    "basal_rate": "a.u./h, > 0: transcription with the element's regulation absent",
    "max_rate": "a.u./h, > 0: the regulated term's ceiling",
}
#: Hill assumptions, declared, not measured: threshold K and coefficient n of H(x) = x^n / (K^n + x^n)
HILL_THRESHOLD, HILL_COEFFICIENT = 1.0, 2.0
#: H at the intact state under those assumptions
H_INTACT = INTACT**HILL_COEFFICIENT / (HILL_THRESHOLD**HILL_COEFFICIENT + INTACT**HILL_COEFFICIENT)
#: the runtime's combination rules, named as model assumptions (lane-sign, b1f3405 and 5cbce26): a
#: gene's activators combine as a MEAN of s_i * H(x_i), so a second activator that is low where the
#: first is high HALVES the drive there; adding an activator can lower expression under this rule.
#: Inhibitors multiply. The bridge does not change either rule (one rule change per registration).
ACTIVATOR_COMBINATION = "mean"
INHIBITOR_COMBINATION = "product"
#: the mapping, one mechanism alone on its gene in its cell (b = basal_rate, V = max_rate, h = H_INTACT):
#:   activates (RHO < 1):  RHO = b / (b + V*s*h)            =>  s = b * (1/RHO - 1) / (V*h)
#:   inhibits  (RHO > 1):  RHO = (b + V) / (b + V*(1 - s*h)) =>  s = (b + V) * (1 - 1/RHO) / (V*h)
#: a fitted s outside (0, 1] means the declared parameters cannot produce the observed response
MAPPING = {
    "activates": "s = b * (1/RHO - 1) / (V * h)",
    "inhibits": "s = (b + V) * (1 - 1/RHO) / (V * h)",
}

# ---- one mechanism, one parameter -----------------------------------------------------------------
#: a mechanism is (element, target gene, cell context); `<id>` and `<id>_measured` are one element
MECHANISM_KEY = ("element id without the _measured suffix", "target gene", "when clauses")
#: which representation a mechanism is simulated from when several exist; the rest are recorded as
#: superseded, never added. Another citation of the same kind with the same observation is the same
#: parameter; with a different observation it is a conflict, reported and not averaged
PRECEDENCE = ("experimental", "curated", "predicted", "inferred", "none")
CITATION_TOLERANCE = 1e-6  # |RHO_a - RHO_b| within this is the same observation

# ---- the diagnostic -------------------------------------------------------------------------------
#: every reason a compiled mechanism does not become a parameter; each is reported by name
UNRESOLVED = {
    "regulator_state_missing": "the source is not a species, not clamped and not a declared zero",
    "gene_parameter_missing": "the target declares no max_rate or no basal_rate > 0",
    "observation_missing": "no effect with a unit could be read from the rule",
    "observation_censored": "the only number is a strength clipped at 1",
    "observation_summarised": "an experimental strength is the maximum of several pairs, not one reading",
    "not_identifiable": "more than one mechanism regulates the gene in this context",
    "conflicting_observations": "two citations of one mechanism report different responses",
    "response_out_of_range": "the fitted strength falls outside (0, 1] for the declared parameters",
}

# ---- the double-counting audit, asked before it was answered ---------------------------------------
AUDIT_QUESTION = (
    "In the 24 compiled programs, how many (element, gene, cell) mechanisms carry both a predicted rule"
    " and a measured rule that one run in that cell would integrate side by side, and how many"
    " (gene, cell) pairs have more than one active regulatory rule, so that no single-deletion"
    " observation identifies a strength under the mean rule?"
)

# ---- acceptance, fixed before the build -----------------------------------------------------------
ACCEPTANCE = (
    "end-to-end: a program compiled from one element whose deletion is observed at log2 fold change -1"
    " in K562, parameterised with declared basal_rate and max_rate, reproduces expression(removed) /"
    " expression(intact) = 0.5 within 1% when the runtime is run with the element clamped at 1 and at 0",
    "diagnostic: running a program whose regulator has no state raises UnresolvedModel naming"
    " regulator_state_missing in strict mode and records it on the trajectory otherwise; a source"
    " renamed '<id>@zero' (network_experiment's edge knockout) is a declared zero and raises nothing",
    "citations: adding a second citation of the same mechanism (a duplicate rule citing another source"
    " with the same observation, or a predicted rule beside a measured one) leaves the simulated"
    " strength unchanged; two representations are never summed or averaged",
)


# ==== the build (after the registration in 06b7557 and its amendment) ==============================

_EFFECT = re.compile(r"effect ([+-]?\d+(?:\.\d*)?(?:e[+-]?\d+)?) log2 fold change")
_LINK = re.compile(r"(\S+) in (\S+) from (\d+) pairs")
MEASURED_SUFFIX = "_measured"


@dataclass(frozen=True, slots=True)
class Observation:
    """One observed removal response: RHO = expression(removed) / expression(intact)."""

    rho: float
    kind: str  # the evidence kind it came from
    how: str  # how RHO was read


@dataclass(slots=True)
class Bridged:
    """A simulation built from annotation, and everything that did not make it in, by name."""

    module: Module  # one fitted rule per resolved mechanism; use it only when `unresolved` is empty
    context: dict[str, str]
    strengths: dict[tuple, float] = field(default_factory=dict)  # mechanism -> fitted s
    observations: dict[tuple, Observation] = field(default_factory=dict)
    required_state: dict[str, float] = field(default_factory=dict)  # clamp these, intact = 1
    superseded: list[Rule] = field(default_factory=list)  # lower-precedence citations, never added
    unresolved: list[Unresolved] = field(default_factory=list)


def element_of(source: str) -> str:
    """`<id>` and `<id>_measured` are one element."""
    return source[: -len(MEASURED_SUFFIX)] if source.endswith(MEASURED_SUFFIX) else source


def mechanism_key(rule: Rule) -> tuple:
    return (element_of(rule.source), rule.target, tuple(sorted(rule.when.items())))


def _kind(rule: Rule) -> str:
    return str(getattr(rule.evidence.kind, "value", rule.evidence.kind))


def observation_of(rule: Rule, module: Module) -> Observation | Unresolved:
    """The removal response a compiled rule states, or why it states none (OBSERVATION_UNITS)."""
    kind = _kind(rule)
    if kind == "predicted":
        m = _EFFECT.search(rule.evidence.note)
        if m:
            return Observation(2.0 ** float(m.group(1)), kind, "evidence note, log2 fold change")
        if 0.0 < rule.strength < 1.0:
            x = -rule.strength if rule.action is Action.ACTIVATE else rule.strength
            return Observation(2.0**x, kind, "rule strength |log2 fold change|, below the clip")
        return Unresolved("observation_censored", rule.source, f"strength {rule.strength} and no effect note")
    if kind == "experimental":
        el = module.entities.get(rule.source)
        cell = rule.when.get("cell_type", "")
        n = None
        for g, c, k in _LINK.findall(el.evidence.note if el is not None else ""):
            if g == rule.target and c == cell:
                n = int(k)
        if n != 1:
            said = "no pair count stated" if n is None else f"the largest of {n} pairs"
            return Unresolved("observation_summarised", rule.source, f"{rule.target} in {cell}: {said}")
        f = -rule.strength if rule.action is Action.ACTIVATE else rule.strength
        return Observation(1.0 + f, kind, "rule strength |EffectSize| of one CRISPRi pair")
    return Unresolved("observation_missing", rule.source, f"{kind} evidence states no removal response")


def fit(rho: float, basal: float, max_rate: float) -> float:
    """The strength that makes one mechanism alone reproduce RHO (MAPPING); may fall outside (0, 1]."""
    if rho < 1.0:
        return basal * (1.0 / rho - 1.0) / (max_rate * H_INTACT)
    return (basal + max_rate) * (1.0 - 1.0 / rho) / (max_rate * H_INTACT)


def parameterize(
    module: Module,
    context: dict[str, str],
    genes: dict[str, dict[str, float]] | None = None,
    build: bool = True,
) -> Bridged:
    """One fitted rule per resolved mechanism active in `context`; every other one reported by name.

    `build=False` classifies without copying the module or writing any rule (the audit's mode).
    """
    genes = genes or {}
    out = copy.deepcopy(module) if build else module
    regulatory = (Action.ACTIVATE, Action.INHIBIT)
    active = [r for r in module.active_rules(context) if r.action in regulatory]
    if build:
        out.rules = [r for r in out.rules if r.action not in regulatory]
    b = Bridged(module=out, context=dict(context))

    by_key: dict[tuple, list[Rule]] = {}
    for r in active:
        by_key.setdefault(mechanism_key(r), []).append(r)
    by_gene: dict[str, list[tuple]] = {}
    for key in by_key:
        by_gene.setdefault(key[1], []).append(key)

    for gene, keys in sorted(by_gene.items()):
        if len(keys) > 1:
            b.unresolved.append(
                Unresolved(
                    "not_identifiable",
                    gene,
                    f"{len(keys)} mechanisms ({', '.join(k[0] for k in keys)}) under the"
                    f" {ACTIVATOR_COMBINATION}/{INHIBITOR_COMBINATION} rules",
                )
            )
            continue
        (key,) = keys
        reps = sorted(by_key[key], key=lambda r: PRECEDENCE.index(_kind(r)) if _kind(r) in PRECEDENCE else 99)
        top = _kind(reps[0])
        chosen = [r for r in reps if _kind(r) == top]
        b.superseded += [r for r in reps if _kind(r) != top]
        obs = [observation_of(r, module) for r in chosen]
        missing = [o for o in obs if isinstance(o, Unresolved)]
        if missing:
            b.unresolved.append(missing[0])
            continue
        rhos = [o.rho for o in obs if isinstance(o, Observation)]
        if max(rhos) - min(rhos) > CITATION_TOLERANCE:
            b.unresolved.append(
                Unresolved(
                    "conflicting_observations",
                    key[0],
                    f"{gene}: {', '.join(f'{x:.4g}' for x in rhos)}",
                )
            )
            continue
        g = out.entities.get(gene)
        basal = genes.get(gene, {}).get("basal_rate", getattr(g, "basal_rate", 0.0))
        vmax = genes.get(gene, {}).get("max_rate", (g.attrs.get("max_rate") if g is not None else None))
        if g is None or not basal or basal <= 0 or not vmax or vmax <= 0:
            b.unresolved.append(
                Unresolved("gene_parameter_missing", gene, "needs basal_rate > 0 and max_rate > 0 (a.u./h)")
            )
            continue
        ob = obs[0]
        assert isinstance(ob, Observation)
        s = fit(ob.rho, basal, vmax) if ob.rho != 1.0 else 0.0
        if not 0.0 < s <= 1.0:
            b.unresolved.append(
                Unresolved("response_out_of_range", key[0], f"{gene}: RHO {ob.rho:.4g} needs s = {s:.4g}")
            )
            continue
        b.strengths[key] = s
        b.observations[key] = ob
        b.required_state[key[0]] = INTACT
        if not build:
            continue
        g.basal_rate = float(basal)
        g.attrs["max_rate"] = float(vmax)
        src = chosen[0]
        out.rules.append(
            Rule(
                id=f"{key[0]}__{gene}",
                source=key[0],
                action=Action.ACTIVATE if ob.rho < 1.0 else Action.INHIBIT,
                target=gene,
                strength=s,
                threshold=HILL_THRESHOLD,
                hill=HILL_COEFFICIENT,
                when=dict(src.when),
                evidence=src.evidence,
            )
        )
    return b
