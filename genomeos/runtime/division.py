# SPDX-License-Identifier: Apache-2.0
# Part of the BioLang engine (language, IR, VM, standard library); see LICENSING.md.
"""Stage 4's arithmetic: what a cell's contents do between divisions, and what a division does to them.

Before this module a `divide` event gave each daughter what the parent had, which is matter made out
of nothing — the failure §9.1 item 4 already refused once for stoichiometry. Here a division with
`partition: contents` gives the two daughters the parent's contents *between* them, so the totals are
equal before and after, and the per-cell amount halves. That is the whole of stage 4's mechanism, and
it is deliberately small: the claim it makes is about conservation, not about a cell cycle.

**What it does not decide.** Nothing here supplies a half-life. A protein's `half_life` is the
program's fact, with the program's evidence, and a protein that declares none is **not** given one:
it is not degraded at all and is reported under `dilution_only`, because inventing a half-life and
inventing infinity are both inventions, and only the second one can be labelled. That labelling is
the point — a protein whose only removal is dilution is the one case where a division changes an
answer that turnover could not reach (`runtime/division_gate.py` registers why).

**The three modes, and why the default is the old behaviour.**

- `duplicate` — each daughter gets what the parent had. Not conservative, and the default, because it
  is what every committed program has been run under. Changing that silently would move results that
  nobody re-ran; a program opts in to stage 4 by saying so.
- `contents` — partition. A species is split binomially when the program is counting molecules and
  there are few enough of them for one molecule to matter (a whole number, strictly below the
  regime's `threshold`, under `units: copies`); otherwise by exact halves. Both branches conserve,
  and each division records which branch each species took, so the choice is in the result rather
  than in this docstring only.
- `binomial` — always stochastic. It refuses a species whose amount is not a whole number, since a
  molecule goes to one daughter or the other and half a molecule goes nowhere.

**The checkpoint is the `when` clause the event already had.** §8 asks for "a size or resource
checkpoint"; a size checkpoint is `when: size >= 2`, read through `ir.model.matches` like every other
guard in the language, and a resource checkpoint is the event's `cost` clauses against stage 2's
`Economy`: a cell that can pay in full divides, and one that cannot does **not** divide partially. A
half-paid division is not a biological state, so a factor below 1.0 postpones the division rather
than shrinking it.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

from genomeos.ir import UNKNOWN, Event, Module, Protein, Regime

if TYPE_CHECKING:  # pragma: no cover - typing only, and the engine must not need the economy
    from genomeos.runtime.economy import Economy


class DivisionError(ValueError):
    """A division that cannot be carried out as the program wrote it."""


@dataclass
class Split:
    """One division: what each daughter got, and by which rule each species was split."""

    a: dict[str, float]
    b: dict[str, float]
    rule: dict[str, str]  # species -> "binomial" | "halves" | "duplicated"

    def conserves(self, parent: dict[str, float]) -> bool:
        """Exactly, not approximately: a partition that loses a molecule is a bug, not a rounding."""
        return all(self.a[s] + self.b[s] == amount for s, amount in parent.items())


@dataclass
class Division:
    """A `divide` event, the half-lives the program declared, and the regime that sizes a molecule.

    `half_lives` maps a species id to hours, or to UNKNOWN for one the program does not give. It is
    built from the module's proteins and is never filled in from anywhere else.
    """

    event: Event
    half_lives: dict[str, Any] = field(default_factory=dict)
    threshold: float = 50.0
    units: str = "au"

    @classmethod
    def from_module(cls, module: Module, event_id: str | None = None) -> Division:
        events = [e for e in (module.events or []) if event_id is None or e.id == event_id]
        if not events:
            raise DivisionError(
                f"no event {event_id!r} in this module"
                if event_id
                else "this module declares no event, so nothing divides"
            )
        if len(events) > 1:
            raise DivisionError(
                "this module declares more than one event and none was named: "
                f"{', '.join(sorted(e.id for e in events))}"
            )
        regime = module.regime or Regime()
        half_lives = {e.id: e.half_life_h for e in module.entities.values() if isinstance(e, Protein)}
        return cls(
            event=events[0],
            half_lives=half_lives,
            threshold=float(regime.threshold),
            units=regime.units,
        )

    @property
    def mode(self) -> str:
        return self.event.partition

    @property
    def dilution_only(self) -> tuple[str, ...]:
        """Species the program gives no half-life: not degraded, and labelled rather than assumed."""
        return tuple(sorted(s for s, t in self.half_lives.items() if t is UNKNOWN))

    # ---- between divisions -----------------------------------------------------------------
    def decay(self, amounts: dict[str, float], hours: float) -> dict[str, float]:
        """Degradation over `hours`, by each species' own declared half-life and by nothing else."""
        if hours < 0:
            raise DivisionError(f"cannot run a cell backwards: {hours} h")
        out = {}
        for species, amount in amounts.items():
            t = self.half_lives.get(species, UNKNOWN)
            out[species] = amount if t is UNKNOWN else amount * 2.0 ** (-hours / float(t))
        return out

    # ---- the checkpoint --------------------------------------------------------------------
    def ready(
        self,
        context: dict[str, str] | None = None,
        economy: Economy | None = None,
    ) -> tuple[bool, str]:
        """Whether the division may fire, and the reason when it may not.

        The size checkpoint is the event's own `when`, through the one matcher every block shares.
        The resource checkpoint is the event's `cost` clauses: a cell that cannot pay in full waits.
        """
        if not self.event.applies(context or {}):
            return False, f"checkpoint unmet: {self.event.when} against {context or {}}"
        if economy is not None and self.event.costs:
            factor = economy.scale({self.event.id: 1.0}).get(self.event.id, 1.0)
            if factor < 1.0:
                return False, f"cannot pay for the division: factor {factor:g} of what it costs"
        return True, "ready"

    # ---- the division ----------------------------------------------------------------------
    def split(self, amounts: dict[str, float], rng: random.Random | None = None) -> Split:
        """The two daughters. Every branch but `duplicate` conserves each species exactly."""
        if self.mode not in ("duplicate", "contents", "binomial"):
            raise DivisionError(f"unknown partition mode {self.mode!r}")
        rng = rng or random.Random(0)
        a: dict[str, float] = {}
        b: dict[str, float] = {}
        rule: dict[str, str] = {}
        for species in sorted(amounts):
            amount = amounts[species]
            if amount < 0:
                raise DivisionError(f"{species} holds {amount}, and a negative amount cannot be split")
            if self.mode == "duplicate":
                a[species] = b[species] = amount
                rule[species] = "duplicated"
                continue
            if self._stochastic(amount):
                drawn = float(rng.binomialvariate(int(amount), 0.5))
                a[species], b[species] = drawn, amount - drawn
                rule[species] = "binomial"
            else:
                a[species] = amount / 2.0
                b[species] = amount - a[species]
                rule[species] = "halves"
        return Split(a=a, b=b, rule=rule)

    def _stochastic(self, amount: float) -> bool:
        """Whether one molecule matters here. Stated rather than inferred inside the loop."""
        whole = float(amount).is_integer()
        if self.mode == "binomial":
            if not whole:
                raise DivisionError(
                    f"partition: binomial needs whole molecules and got {amount}; a molecule goes to "
                    "one daughter or the other, so use `partition: contents` for a continuous amount"
                )
            return True
        return whole and amount < self.threshold and self.units == "copies"

    # ---- several divisions -----------------------------------------------------------------
    def run(
        self,
        amounts: dict[str, float],
        divisions: int,
        interval_h: float,
        context: dict[str, str] | None = None,
        economy: Economy | None = None,
        seed: int = 0,
    ) -> dict[str, Any]:
        """`divisions` cycles of: degrade for `interval_h`, test the checkpoint, divide if it holds.

        One daughter is followed, which is what a per-cell amount means. The trace carries every
        cycle, including the ones where the checkpoint refused, so a run with no divisions is
        distinguishable from a run that never started.
        """
        rng = random.Random(seed)
        current = dict(amounts)
        trace: list[dict[str, Any]] = [{"division": 0, "hours": 0.0, "amounts": dict(current)}]
        fired = 0
        conserved = True
        for n in range(1, divisions + 1):
            grown = self.decay(current, interval_h)
            may, why = self.ready(context, economy)
            if may:
                got = self.split(grown, rng)
                conserved = conserved and (self.mode == "duplicate" or got.conserves(grown))
                current = got.a
                fired += 1
                row = {"division": n, "divided": True, "rule": got.rule}
            else:
                current = grown
                row = {"division": n, "divided": False, "why": why}
            trace.append({**row, "hours": n * interval_h, "amounts": dict(current)})
        return {
            "mode": self.mode,
            "event": self.event.id,
            "divisions_asked": divisions,
            "divisions_fired": fired,
            "interval_h": interval_h,
            "seed": seed,
            "dilution_only": list(self.dilution_only),
            "conserved": conserved,
            "trace": trace,
            "final": dict(current),
        }

    # ---- what the arithmetic predicts, written out rather than fitted ----------------------
    def expected(
        self,
        start: float,
        half_life_h: Any,
        divisions: int,
        interval_h: float,
        cycles: int | None = None,
    ) -> float:
        """The closed form a run must reproduce: `A0 · 2^-n · 2^(-mT/t)`, or without the `2^-n`
        when the mode duplicates, or without the turnover term when no half-life is declared.

        `divisions` is how many times the cell actually divided and `cycles` how many intervals it
        lived through, which are the same number only when the checkpoint never refused. Keeping
        them apart is what makes the checkpoint arm a test rather than a tautology: a cell whose
        checkpoint is never met still ages.
        """
        hours = (divisions if cycles is None else cycles) * interval_h
        dilution = 1.0 if self.mode == "duplicate" else 2.0**-divisions
        turnover = 1.0 if half_life_h is UNKNOWN else 2.0 ** (-hours / float(half_life_h))
        return start * dilution * turnover

    def steady_state(self, synthesis: float, half_life_h: Any, interval_h: float | None) -> float:
        """`s / (k_degradation + k_division)`. `interval_h` None is a cell that never divides; a
        species with no half-life and no division has no removal at all, which is stated rather than
        returned as an infinity for somebody else to trip over."""
        k_deg = 0.0 if half_life_h is UNKNOWN else math.log(2.0) / float(half_life_h)
        k_div = 0.0 if interval_h is None else math.log(2.0) / float(interval_h)
        if k_deg + k_div == 0.0:
            raise DivisionError(
                "a species with no declared half-life in a cell that does not divide has no removal, "
                "so it has no steady state: say which of the two the program means"
            )
        return synthesis / (k_deg + k_div)
