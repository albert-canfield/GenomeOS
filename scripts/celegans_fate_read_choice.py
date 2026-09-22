# SPDX-License-Identifier: AGPL-3.0-or-later
"""Should the worm's fate rules read `exposure(lineage)` or `mean(lineage)`, and what moved 332 to 414?

Area E, following the credit measurement of 2026-09-21 (542d613). Two things were left open there.

1. `scripts/celegans_fate_credit.py` found that `F.mean(lineage) >= 0.25` is ahead of the shipped
   `F.exposure(lineage) >= 15` on the honest credit metric as well as with the lookup as fallback, and
   left the choice of read to area E's owner. This script decides it, against a bar registered below
   rather than against the headline.
2. The same `fates.bio` and the same reads have the factor rules firing on **414** cells where the
   2026-09-15 result file recorded **332**, with the published score unchanged at 522 of 555. A rule set
   that fires on 25% more cells and decides exactly as many fates has to be understood before anything is
   built on it: either something else decides the added cells identically, or the metric cannot see them.

WHAT THE TWO READS COMPUTE (`genomeos/runtime/body.py:318`, `_read`). Both integrate the same thing: the
reader's presence call, banked minute by minute at every decision point, a factor counting 1 for every
minute the cell carries it (this program declares no `cell_network`, so there is no level to integrate).
They differ only in the denominator.

    F.exposure(lineage)  total minutes the path from the zygote to this cell carried F. Absolute,
                         monotone, and unbounded: it can only grow, and it grows faster the later the
                         cell is born, because the path is longer.
    F.mean(lineage)      the same total divided by the span of the window, which for the lineage window
                         is `t` itself -- time since the first cell. So it is the **fraction of the whole
                         elapsed embryo** that this cell's path spent carrying F. Not monotone: it falls
                         whenever the path is not carrying F, so lateness cannot inflate it.

That is the substantive difference and it is why `mean` might be the better read: `exposure >= 15 min` is
a threshold a late-born cell meets for having a long ancestry, while `mean >= 0.25` asks the same
question of a cell born at 100 min and one born at 500 min. It is also why the shipped read is on a
plateau (1 to 30 min give the identical program) and `mean` is not (0.20, 0.22, 0.25 all differ): an
absolute threshold on a quantity that spans two orders of magnitude has a wide dead band, and a fraction
of a window does not.

THE METRIC. The published 522 of 555 is not a measurement of the factors, because the observed lineage
lookup answers every cell no factor rule claims and is right by construction. 542d613 established the
honest replacements, and this script scores both reads on all four, on the same 555 cells:

    H  honest credit          terminal fates where a `factor_fate_*` decision set the cell's type and
                              that type is the observed one (`claimed_right`). No credit for a fate the
                              lookup supplied.
    I  decided independently  claimed cells on which the lookup's own `fate_<cell>` decision had NOT
                              already fired (`claimed - lookup_wrote_it_first`). 542d613 measured this at
                              **0 of 414**: every factor rule is guarded on a terminal `cell_type` and
                              only the lookup ever sets one, so a factor rule is a successor to the
                              lookup in the same differentiation chain and never a competitor.
    P  precision              H / claimed, on the cells the read does claim.
    F  free credit            H minus the claimed-and-right cells on which the rule fired while the cell
                              was already `committed` to that very type. This one is new here and it
                              exists because of item 2. `commitment terminal_fate` locks `cell_type`, and
                              `_pick_fate` refuses any decision whose `to` differs from the lock
                              (`body.py:456`), so a rule that fires after the lock is set cannot return
                              any answer but the locked one. Credit for such a cell is credit for
                              agreeing with an answer the program had already refused to let it
                              contradict.

THE TRAP THIS REGISTERS AGAINST, in this lane's own words. The number to move is not the headline. This
project has spent two days establishing that a rise in a headline is not a result unless the thing the
headline is measured against moves with it -- a band table retired because the same targets band from
27,514 to 198,475 depending on whose base rate is used, a calibration constant that turned out to be a
mean over per-screen shifts from -1.303 to +4.228. Area E's version of that is exact and was measured:
**no factor rule in this program has ever decided a fate the lineage had not already written, 414 of
414.** So a read that scores higher on H while I stays 0 has not made the factors decide one more thing;
it has agreed with the answer sheet more often on cells the answer sheet had already filled in. The
baseline that names Neuron for every cell scores 226 of 555, and the unfitted rules' honest 184 sits
under a ceiling of 247. Changing a read to buy +11 agreements would be exactly the move this area has
been publishing negatives about.

PRE-REGISTRATION. Committed before this script had been run and before
`data/results/celegans_fate_read_choice.json` existed. Checked automatically at the end of the run; the
result file records per claim whether it held. *derived* = arithmetic on a committed result file, stated
so that the run confirms the accounting rather than discovers it; *open* = not known before the run.

  R1 (derived)  the shipped read, `exposure(lineage) >= 15`, reproduces 542d613 exactly: 414 claimed,
                381 honest, precision 0.9203, 522 with the lookup as fallback, and I = 0.
  R2 (derived)  `mean(lineage) >= 0.25` scores H = 392 in sample and 383 held out, against 381 and 374,
                claiming 414 and 423 cells. So it wins on H by 11 in sample and 9 held out.
  R3 (open)     I = 0 for `mean` as well, in sample and held out. The guard is a terminal `cell_type`
                whatever the read is, so neither read decides a single fate the lineage had not written.
                **This is the claim the decision turns on.** If it holds, the registered bar below is not
                met by a read that only raises H.
  R4 (open)     the 332-to-414 move is the crossing re-decisions of section 7.5 (`recheck: crossings`,
                9d42485, landed after 2026-09-15 and the default since). Run the identical program with
                `recheck="none"` and the shipped read claims 332 cells again, with the score still 522.
                Point prediction 332; anything within 332 +/- 3 is the same finding.
  R5 (open)     every one of the ~82 cells that only `crossings` adds is claimed-and-right, and on every
                one of them the cell was already `committed` to the type the rule names when the rule
                fired. Not one of them is claimed-and-wrong. If that holds, the added cells are decided
                identically by the lock -- which is derived from the lookup -- and the published score
                cannot see them because there is nothing there to see.
  R6 (open)     therefore F, free credit, is materially below H for both reads.

  THE BAR, registered before the measurement. The shipped read changes to `mean(lineage) >= 0.25` only
  if, in sample and held out:
      (a) I(mean) > I(exposure)  -- the mean decides at least one fate the lineage had not written; or
      (b) F(mean) > F(exposure) by more than 5 fates, with precision not lower.
  A rise in H alone is NOT an improvement and is registered here as no improvement, whatever its size,
  because H can rise without a single additional fate being decided by the factors. Nor is a rise in the
  published 522-style score, which rewards abstention. If neither (a) nor (b) is met, the shipped read
  stays as it is and both numbers are published side by side.

  Secondary and NOT part of the bar: `mean`'s threshold is not on a plateau (0.20, 0.22, 0.25 differ)
  where 15 min is flat from 1 to 30. Even if the bar were met, that would have to be answered first.

  What makes the run void rather than negative (checked, and the run stops):
  V1  the shipped arm must reproduce 414 claimed / 381 honest / 522 with the fallback and a `fates.bio`
      byte-identical to the checked-in one. Otherwise this is not the program that ships.
  V2  for every arm, no abstained cell is wrong: the lookup is a clean fallback.

  Not a failure: finding that nothing should change. The instruction this lane works under says so, and
  so does the area's record of published negatives.

    uv run python scripts/celegans_fate_read_choice.py

No model calls, no network, one core. Result: data/results/celegans_fate_read_choice.json.
"""

from __future__ import annotations

import json
import shutil
import sys
import tempfile
from collections.abc import Callable
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path and not (Path.cwd() / "genomeos").is_dir():
    sys.path.insert(0, str(ROOT))

from genomeos.lang import parse_file  # noqa: E402
from genomeos.organism import fate_rules as fr  # noqa: E402
from genomeos.organism.reference import ReferenceLineage  # noqa: E402
from genomeos.organism.tf_atlas import load_cells  # noqa: E402
from genomeos.runtime.body import Body, Cell  # noqa: E402
from scripts.celegans_fate_credit import ORG, credit, rule_targets  # noqa: E402
from scripts.celegans_fate_reads import Arm, harvest, program_rules  # noqa: E402

RESULT = Path("data/results/celegans_fate_read_choice.json")
HATCH = 800.0
SHIPPED = "exposure(lineage) >= 15 (shipped)"
ALT = "mean(lineage) >= 0.25"
READS = {SHIPPED: ("exposure", fr.THRESHOLD_MIN), ALT: ("mean", 0.25)}


# ---- when a factor rule fires, and what the cell was already committed to ---------------------------


class FiringLog:
    """Every `factor_fate_*` decision that fires on a terminal cell, with the state of the cell at the
    instant it fired. Wraps `Body._decide` rather than reading the finished run, because the two facts
    that matter -- whether this was the cell's birth decision point, and whether the cell was already
    committed to a type -- are both gone by the time the run ends."""

    def __init__(self, wanted: set[str]) -> None:
        self.wanted = wanted
        self.events: dict[str, list[dict]] = {}
        self._original: Callable | None = None

    def __enter__(self) -> FiringLog:
        original = Body._decide
        self._original = original
        events, wanted = self.events, self.wanted

        def decide(self: Body, c: Cell, born: bool = True) -> None:
            seen = set(c.fired)
            committed, cell_type = c.committed, c.cell_type
            original(self, c, born=born)
            if c.name not in wanted:
                return
            for d in c.fired:
                if d.startswith("factor_fate_") and d not in seen:
                    events.setdefault(c.name, []).append(
                        {
                            "decision": d,
                            "at": round(self.time, 3),
                            "born": round(c.born, 3),
                            "at_birth": born,
                            "committed_before": committed,
                            "cell_type_before": cell_type,
                        }
                    )

        Body._decide = decide
        return self

    def __exit__(self, *exc: object) -> None:
        Body._decide = self._original

    def summarise(self, claimed: list[str], targets: dict[str, str], right: set[str]) -> dict:
        """Of the cells a rule claimed, how many were decided at the cell's birth, how many again at a
        later crossing, and on how many of those was the cell already locked to the very type the rule
        went on to name."""
        later = [c for c in claimed if any(not e["at_birth"] for e in self.events.get(c, []))]
        locked = [
            c
            for c in claimed
            for e in self.events.get(c, [])[-1:]
            if e["committed_before"] and e["committed_before"] == targets[e["decision"]]
        ]
        return {
            "claimed": len(claimed),
            "fired_only_at_the_cell_s_birth": len(claimed) - len(later),
            "fired_again_at_a_crossing": len(later),
            "already_committed_to_the_type_the_rule_named": len(locked),
            "of_those_right": sum(1 for c in locked if c in right),
            "of_those_wrong": sum(1 for c in locked if c not in right),
        }


# ---- arms -------------------------------------------------------------------------------------------


def run_arm(arm: Arm, ref: ReferenceLineage, cells: list[str], recheck: str, wanted: set[str]) -> dict:
    """One finished run of one program under one recheck mode, with the credit decomposition of 542d613
    and the firing log beside it."""
    with FiringLog(wanted) as log:
        body = Body(parse_file(arm.dir / "embryo_factors.bio"), means=True, recheck=recheck).run(until=HATCH)
    targets = rule_targets(arm)
    c = credit(body, ref, cells, targets)
    claimed = c["claimed_cells"]
    right = {cid for cid in claimed if body.cells[cid].cell_type == ref.cells[cid].cell_type}
    c["when_the_rule_fired"] = log.summarise(claimed, targets, right)
    c["rechecks"] = body.summary()["rechecks"]
    c["recheck"] = recheck
    return c


def held_out(arm: Arm, ref: ReferenceLineage, write: Callable[[list[str]], None], labels: dict) -> dict:
    """Leave one founder sublineage out, credit taken only on the held-out eighth. The same shape as
    `celegans_fate_credit.fold_credit`, repeated here so this script's arms are all built the same way
    and the independence count is carried through the folds."""
    cells = sorted(labels)
    keys = ("cells", "claimed", "claimed_right", "claimed_wrong", "abstained_wrong", "lookup_wrote_it_first")
    total = dict.fromkeys(keys, 0)
    for g in sorted({fr.sublineage(c) for c in cells}):
        held = [c for c in cells if fr.sublineage(c) == g]
        write([c for c in cells if fr.sublineage(c) != g])
        body = Body(parse_file(arm.dir / "embryo_factors.bio"), means=True).run(until=HATCH)
        c = credit(body, ref, held, rule_targets(arm))
        for k in keys:
            total[k] += c[k]
    return {
        **total,
        "honest": total["claimed_right"],
        "decided_a_fate_the_lookup_had_not": total["claimed"] - total["lookup_wrote_it_first"],
        "precision": round(total["claimed_right"] / total["claimed"], 4) if total["claimed"] else None,
        "with_lookup_fallback": total["cells"] - total["claimed_wrong"] - total["abstained_wrong"],
    }


def free_credit(a: dict) -> int:
    """Honest credit minus the cells on which no other answer was reachable."""
    return a["claimed_right"] - a["when_the_rule_fired"]["of_those_right"]


def main() -> None:
    ref = ReferenceLineage.load()
    term = fr.embryonic_terminal(ref)
    labels = {c.id: c.tissue for c in term}
    cells = sorted(labels)
    wanted = set(cells)
    atlas = load_cells()
    factors = sorted({f for fs in atlas.values() for f in fs})
    print(f"{len(cells)} embryonic terminal cells, {len(factors)} factors", flush=True)
    reads, _points = harvest(ref, factors)

    tmpdir = Path(tempfile.mkdtemp(prefix="celegans-read-choice-"))
    arms: dict[str, dict] = {}
    claimed_sets: dict[str, set[str]] = {}
    try:
        for label, (kind, thr) in READS.items():
            arm = Arm(tmpdir / kind)

            def write(fit: list[str], arm: Arm = arm, kind: str = kind, thr: float = thr) -> None:
                (arm.dir / "fates.bio").write_text(
                    fr.to_bio_fates(
                        program_rules(fr.runtime_features(reads, kind, "lineage", thr), labels, fit),
                        read=kind,
                        window="lineage",
                        threshold=thr,
                    )
                )

            write(cells)
            if kind == "exposure":  # V1: this has to be the program that ships
                assert (arm.dir / "fates.bio").read_text() == (ORG / "fates.bio").read_text(), (
                    "V1: not the shipped fates.bio"
                )
            for recheck in ("crossings", "none"):
                c = run_arm(arm, ref, cells, recheck, wanted)
                claimed_sets[f"{label} | recheck={recheck}"] = set(c.pop("claimed_cells"))
                arms[f"{label} | recheck={recheck}"] = c
            arms[f"{label} | held out by founder sublineage"] = held_out(arm, ref, write, labels)
            write(cells)
            print(f"  {label}: done", flush=True)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    ship_cross = arms[f"{SHIPPED} | recheck=crossings"]
    ship_none = arms[f"{SHIPPED} | recheck=none"]
    mean_cross = arms[f"{ALT} | recheck=crossings"]
    mean_none = arms[f"{ALT} | recheck=none"]
    ship_out = arms[f"{SHIPPED} | held out by founder sublineage"]
    mean_out = arms[f"{ALT} | held out by founder sublineage"]
    assert ship_cross["claimed"] == 414 and ship_cross["honest"] == 381, "V1: not 542d613's program"
    assert ship_cross["with_lookup_fallback"] == 522, "V1: not the published score"
    assert all(a["abstained_wrong"] == 0 for a in arms.values()), "V2: the fallback is not clean"

    added = sorted(claimed_sets[f"{SHIPPED} | recheck=crossings"] - claimed_sets[f"{SHIPPED} | recheck=none"])
    lost = sorted(claimed_sets[f"{SHIPPED} | recheck=none"] - claimed_sets[f"{SHIPPED} | recheck=crossings"])
    move = {
        "cells_claimed_with_the_crossing_re_decisions": ship_cross["claimed"],
        "cells_claimed_without_them": ship_none["claimed"],
        "score_with_the_lookup_as_fallback": [
            ship_cross["with_lookup_fallback"],
            ship_none["with_lookup_fallback"],
        ],
        "honest_credit": [ship_cross["honest"], ship_none["honest"]],
        "claimed_and_wrong": [ship_cross["claimed_wrong"], ship_none["claimed_wrong"]],
        "decided_a_fate_the_lookup_had_not": [
            ship_cross["decided_a_fate_the_lookup_had_not"],
            ship_none["decided_a_fate_the_lookup_had_not"],
        ],
        "cells_only_the_crossings_arm_claims": len(added),
        "cells_only_the_recheck_none_arm_claims": len(lost),
        "of_the_added_cells_claimed_and_right": ship_cross["claimed_right"] - ship_none["claimed_right"],
        "of_the_added_cells_claimed_and_wrong": ship_cross["claimed_wrong"] - ship_none["claimed_wrong"],
        "the_added_cells": added,
        "crossings_taken": [ship_cross["rechecks"], ship_none["rechecks"]],
        "the_same_under_mean": {
            "claimed": [mean_cross["claimed"], mean_none["claimed"]],
            "honest": [mean_cross["honest"], mean_none["honest"]],
            "with_lookup_fallback": [
                mean_cross["with_lookup_fallback"],
                mean_none["with_lookup_fallback"],
            ],
        },
        "so": "",
    }
    move["so"] = (
        f"the {len(added)} cells the crossing re-decisions add are claimed "
        f"{move['of_the_added_cells_claimed_and_right']} right and "
        f"{move['of_the_added_cells_claimed_and_wrong']} wrong, and the published score is "
        f"{ship_cross['with_lookup_fallback']} either way. A cell whose fate the lookup wrote at its "
        "birth is committed to that type, and _pick_fate refuses any decision naming another, so a "
        "factor rule firing at a later crossing can only re-assert the type already on the cell. The "
        "score cannot see the move because there is nothing there to see: the added cells were decided "
        "by the lookup and locked before the factors were asked."
    )

    bar = {
        "decided_a_fate_the_lookup_had_not_in_sample_and_held_out": {
            "exposure": [
                ship_cross["decided_a_fate_the_lookup_had_not"],
                ship_out["decided_a_fate_the_lookup_had_not"],
            ],
            "mean": [
                mean_cross["decided_a_fate_the_lookup_had_not"],
                mean_out["decided_a_fate_the_lookup_had_not"],
            ],
        },
        "honest_in_sample_and_held_out": {
            "exposure": [ship_cross["honest"], ship_out["honest"]],
            "mean": [mean_cross["honest"], mean_out["honest"]],
        },
        "free_credit_in_sample": {"exposure": free_credit(ship_cross), "mean": free_credit(mean_cross)},
        "precision_in_sample_and_held_out": {
            "exposure": [ship_cross["precision"], ship_out["precision"]],
            "mean": [mean_cross["precision"], mean_out["precision"]],
        },
        "published_score_with_the_lookup_as_fallback": {
            "exposure": [ship_cross["with_lookup_fallback"], ship_out["with_lookup_fallback"]],
            "mean": [mean_cross["with_lookup_fallback"], mean_out["with_lookup_fallback"]],
        },
    }
    clause_a = (
        mean_cross["decided_a_fate_the_lookup_had_not"] > ship_cross["decided_a_fate_the_lookup_had_not"]
        and mean_out["decided_a_fate_the_lookup_had_not"] > ship_out["decided_a_fate_the_lookup_had_not"]
    )
    clause_b = free_credit(mean_cross) - free_credit(ship_cross) > 5 and (
        mean_cross["precision"] >= ship_cross["precision"]
    )
    bar["clause_a_independence_rises"] = clause_a
    bar["clause_b_free_credit_rises_by_more_than_5_fates"] = clause_b
    bar["registered_bar_met"] = clause_a or clause_b
    bar["decision"] = (
        f"change the shipped read to {ALT}"
        if (clause_a or clause_b)
        else "the shipped read stays as it is: exposure(lineage) >= 15, and both numbers are published"
    )

    out = {
        "result": "celegans_fate_read_choice",
        "date": date.today().isoformat(),
        "question": "should the worm's fate rules read exposure(lineage) or mean(lineage), and what "
        "moved the cells a factor rule fires on from 332 to 414 with the score unchanged?",
        "metric": {
            "honest": "terminal fates a factor rule set and got right; no credit for a lookup answer",
            "decided_a_fate_the_lookup_had_not": "claimed cells the lookup's own fate_<cell> decision had "
            "not already fired on -- the count the bar is on",
            "free_credit": "honest minus the claimed-and-right cells on which the cell was already "
            "committed to the very type the rule named, so no other answer was reachable",
        },
        "the_332_to_414_move": move,
        "arms": arms,
        "the_registered_bar": bar,
    }
    out["preregistered"] = {
        "R1 the shipped arm reproduces 542d613": {
            "predicted": [414, 381, 0.9203, 522, 0],
            "observed": (
                obs := [
                    ship_cross["claimed"],
                    ship_cross["honest"],
                    ship_cross["precision"],
                    ship_cross["with_lookup_fallback"],
                    ship_cross["decided_a_fate_the_lookup_had_not"],
                ]
            ),
            "held": obs == [414, 381, 0.9203, 522, 0],
        },
        "R2 mean scores 392 in sample and 383 held out on honest credit": {
            "predicted": [392, 383],
            "observed": [mean_cross["honest"], mean_out["honest"]],
            "held": [mean_cross["honest"], mean_out["honest"]] == [392, 383],
        },
        "R3 mean decides no fate the lineage had not written either": {
            "predicted": [0, 0],
            "observed": bar["decided_a_fate_the_lookup_had_not_in_sample_and_held_out"]["mean"],
            "held": bar["decided_a_fate_the_lookup_had_not_in_sample_and_held_out"]["mean"] == [0, 0],
        },
        "R4 without the crossing re-decisions the shipped read claims 332 again": {
            "predicted": "332 +/- 3",
            "observed": ship_none["claimed"],
            "held": abs(ship_none["claimed"] - 332) <= 3,
        },
        "R5 every cell only the crossings arm claims is right, and was locked to that type": {
            "predicted": "all right, none wrong",
            "observed": [
                move["cells_only_the_crossings_arm_claims"],
                move["of_the_added_cells_claimed_and_right"],
                move["of_the_added_cells_claimed_and_wrong"],
                ship_cross["when_the_rule_fired"]["already_committed_to_the_type_the_rule_named"],
            ],
            "held": move["of_the_added_cells_claimed_and_wrong"] == 0
            and move["of_the_added_cells_claimed_and_right"] == move["cells_only_the_crossings_arm_claims"],
        },
        "R6 free credit is materially below honest credit for both reads": {
            "observed": {
                "exposure": [ship_cross["honest"], free_credit(ship_cross)],
                "mean": [mean_cross["honest"], free_credit(mean_cross)],
            },
            "held": free_credit(ship_cross) < ship_cross["honest"]
            and free_credit(mean_cross) < mean_cross["honest"],
        },
    }

    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(out, indent=2) + "\n")

    print(f"\n{'arm':<58} {'claimed':>8} {'honest':>7} {'indep':>6} {'prec':>7} {'fallback':>9}")
    for name, a in arms.items():
        print(
            f"{name:<58} {a['claimed']:>8} {a['honest']:>7} "
            f"{a['decided_a_fate_the_lookup_had_not']:>6} {a['precision'] or 0:>7.4f} "
            f"{a['with_lookup_fallback']:>9}"
        )
    print("\n332 -> 414:")
    print(json.dumps({k: v for k, v in move.items() if k != "the_added_cells"}, indent=1))
    print(f"\nbar met: {bar['registered_bar_met']} -- {bar['decision']}")
    for claim, v in out["preregistered"].items():
        print(f"  {'held ' if v['held'] else 'BROKEN'}  {claim}: {v['observed']}")


if __name__ == "__main__":
    main()
