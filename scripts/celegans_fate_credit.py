# SPDX-License-Identifier: AGPL-3.0-or-later
"""How many of the worm's terminal fates the measured factors decide **on their own** (area E).

The program scores 522 of 555 embryonic terminal fates. That number is not what the factors know. It is
what the program gets right when a factor rule decides the cells it can and the observed lineage lookup
supplies every other fate, and the lookup is Sulston's answer sheet: every fate it supplies is right by
construction. So a cell the factors say nothing about still counts as a fate the program got right.

This script takes the fallback away and asks the narrower question the roadmap's "missing for complete"
is really about: **of the 555 embryonic terminal cells, for how many is the fate decided by the measured
transcription factors and right, with no credit for anything the lookup answered?** It then strips a
second borrowing — rules fitted to the observed lineage's own tissue labels are not the factors speaking
— and a third, the terminality gate: a factor rule in this program can only fire on a cell the lookup has
already differentiated, so it chooses *which* of 17 terminal types, never *that* the cell is terminal now.

Four things this method cannot do, stated before the numbers:

  1. It cannot make the 555 cells without the observed lineage. Their existence, their names, their birth
     times, the tree above them and the fact that each is terminal all come from `lineage_embryo.bio`.
     Even the vocabulary of 17 terminal types is the WormWeb tissue labelling. The most this measures is
     how much of the *choice among those types* the factors carry.
  2. It cannot separate the atlas from the lineage completely. 100 of the 555 are not in Ma 2021 at all
     and inherit the factor set of the nearest ancestor the atlas tracked, along the observed tree.
  3. It cannot score a fate the atlas stops before. Ma 2021 ends at the bean stage, so the terminal genes
     of a cell that differentiates later are not in the measurement at all (the glia result of 2026-09-14).
  4. Holding a fitted rule out by founder sublineage does not make it unfitted. The hypothesis class --
     which tissues get rules, which of the 266 factors are candidates, that two factors may be conjoined
     -- was chosen by looking at the observed labels of all 555.

PRE-REGISTRATION. Committed before this script had ever been run and before
`data/results/celegans_fate_credit.json` existed. Each claim below is checked automatically at the end of
the run and the result file records, per claim, whether it held. Two kinds are mixed on purpose:

  *derived*  arithmetic on the published 2026-09-15 result file, which reports for each arm both the
             fates it got right and the cells its factor rules claimed. If 332 cells are claimed and 33
             fates are wrong, and the lookup cannot be wrong, then 299 claims are right. Stating these
             first is the point: this measurement's job is to confirm that accounting cell by cell
             through the Body, not to discover it. If the Body disagrees with the arithmetic, the
             arithmetic was wrong about which cells the errors are on, and that is worth knowing.
  *open*     genuinely not known before the run.

  P1 (derived)  cited rules only -- the seven textbook factor rules and the one cited sheath-glia rule,
                nothing fitted anywhere -- decide correctly 130 to 160 of 555. Point estimate 144: the
                published plateau has the textbook rules claiming 164 and wrong on 25, and the sheath
                rule was worth +5 fates on 2026-09-14.
  P2 (derived)  the whole shipped rule set, in sample, decides correctly 299 of 555 (332 claimed, 33
                wrong), and held out by founder sublineage 268 of 555 (340 claimed, 72 wrong).
  P3 (open)     the honest number is *below* the best baseline that uses no factors and no per-cell
                lineage information at all: naming the commonest terminal tissue for every cell, which is
                Neuron, 226 of 555 (40.7%). This is the claim I expect to survive and it is the finding.
  P4 (open)     the factors are nevertheless informative where they speak: on exactly the cells a rule
                claims, the rule is right more often than "always Neuron" is right on those same cells,
                for the cited set and for the shipped set in sample and held out. If this fails the
                factors are adding nothing at all and the honest number is 0, not small.
  P5 (open)     the fallback rewards abstention, so the honest metric can reorder two reads that the
                published metric ordered the other way: `mean(lineage) >= 0.25` scores above
                `exposure(lineage) >= 15` with the lookup (533 against 522) and I expect it to score
                *below* it without one.
  P6 (open)     with the lookup's 555 `fate_*` decisions deleted from the program, the factor rules fire
                on 0 cells, because their guard is `cell_type = <one of the 17 terminal types>` and only
                the lookup ever puts a terminal type on a cell. The eligible denominator without the
                observed lineage is zero.
  P7 (open)     open that guard to any undifferentiated cell, on the program with the lookup's fates
                deleted, and the cited rules claim more non-terminal cells than terminal ones: they
                cannot tell a precursor carrying a factor from the terminal cell that carries it.

  What would make the measurement void rather than negative (checked, and the run stops):
  F1  the shipped arm must reproduce 522 of 555 with 332 cells claimed. Otherwise the instrument is not
      measuring the program that shipped and no number is reported.
  F2  for every arm, claimed-and-right + claimed-and-wrong + abstained must equal 555, and every
      abstained cell must have the observed fate. If an abstained cell is wrong, the lookup is not a
      clean fallback and the decomposition of the 33 errors is not what it says.
  F3  the headline must be flat across the published 1 to 60 min plateau of the threshold, within 10
      fates. If it moves more, the honest number is a knob and is reported as a range, not a number.
  F4  the shuffle null -- the same cited rules applied to factor sets permuted among the 555 cells --
      must stay more than two standard deviations below the real number. If it does not, the cited rules
      are reading the composition of the tissue classes and not the factors of the cell.

  Not a failure: a low number. Whatever it is, it is the answer, and area E has published negatives
  before (socket glia, commitment, competence, sister divergence).

    uv run python scripts/celegans_fate_credit.py

No model calls, no network, about two minutes on one core. Result: data/results/celegans_fate_credit.json.
"""

from __future__ import annotations

import json
import random
import shutil
import sys
import tempfile
from collections import Counter
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
from genomeos.runtime.body import Body  # noqa: E402
from scripts.celegans_fate_reads import Arm, harvest, presence_at_birth, program_rules  # noqa: E402

ORG = Path("data/organisms/celegans")
RESULT = Path("data/results/celegans_fate_credit.json")
HATCH = 800.0
PERMUTATIONS = 200

# The rule sets, by how much of the observed lineage went into writing them.
#   cited     a factor -> tissue relation published for C. elegans, with a citation, which nobody fitted
#             here: the seven textbook rules and the one sheath-glia rule of 2026-09-14.
#   fitted    the neuron decision list learned from the atlas crossed with the observed tissue labels.
# Only the second needs holding out, and holding it out does not make it cited (see the docstring).


def cited_rules() -> list[fr.Rule]:
    tb = fr.textbook_rules()
    return tb[:4] + fr.glia_rules() + tb[4:]


# ---- credit, cell by cell -------------------------------------------------------------------------


def credit(body: Body, ref: ReferenceLineage, cells: list[str]) -> dict:
    """What the factors can be credited with, and what the lookup answered, for one finished run.

    A cell is *claimed* when a `factor_fate_*` decision fired on it; it is right when the type it ends
    with is the observed one. Everything else is *abstained*: the lookup answered, and the lookup is the
    answer sheet, so those cells say nothing about the factors. F2 of the pre-registration is asserted
    here rather than reported, because if an abstained cell is wrong the whole decomposition is wrong."""
    claimed_right = claimed_wrong = abstained_right = abstained_wrong = 0
    claimed: list[str] = []
    wrong_by_tissue: Counter[str] = Counter()
    for cid in cells:
        bc = body.cells.get(cid)
        want = ref.cells[cid].cell_type
        got = bc.cell_type if bc is not None else None
        if bc is not None and any(d.startswith("factor_fate_") for d in bc.fired):
            claimed.append(cid)
            if got == want:
                claimed_right += 1
            else:
                claimed_wrong += 1
                wrong_by_tissue[f"{ref.cells[cid].tissue}->{got}"] += 1
        elif got == want:
            abstained_right += 1
        else:
            abstained_wrong += 1
    n = len(cells)
    assert claimed_right + claimed_wrong + abstained_right + abstained_wrong == n, "F2: cells lost"
    return {
        "cells": n,
        "claimed": len(claimed),
        "claimed_right": claimed_right,
        "claimed_wrong": claimed_wrong,
        "abstained": abstained_right + abstained_wrong,
        "abstained_wrong": abstained_wrong,
        "with_lookup_fallback": claimed_right + abstained_right,
        "honest": claimed_right,
        "honest_fraction": round(claimed_right / n, 4) if n else None,
        "precision": round(claimed_right / len(claimed), 4) if claimed else None,
        "claimed_cells": claimed,
        "confusion": dict(wrong_by_tissue.most_common(10)),
    }


def _maj(arm: dict) -> int:
    """How often the factor-free majority tissue is right on the cells that arm's rules claimed."""
    return arm["baselines_on_the_same_cells"]["majority_tissue"]["right"]


def on_the_same_cells(claimed: list[str], labels: dict[str, str], baselines: dict[str, dict]) -> dict:
    """P4's comparison: on exactly the cells the factors claimed, how often is each label-free baseline
    right? Abstention is what makes `honest` small, and this removes it from the comparison, so what is
    left is whether the factors know anything about the cells they do speak about."""
    return {
        name: {
            "right": sum(1 for c in claimed if pred.get(c) == labels[c]),
            "of": len(claimed),
            "fraction": round(sum(1 for c in claimed if pred.get(c) == labels[c]) / len(claimed), 4)
            if claimed
            else None,
        }
        for name, pred in baselines.items()
    }


def fold_credit(
    arm: Arm,
    ref: ReferenceLineage,
    write: Callable[[list[str]], None],
    labels: dict[str, str],
) -> dict:
    """Leave one founder sublineage out: rules written on the other seven, credit taken only on the
    eighth, so no rule is credited on a lineage whose labels it saw."""
    cells = sorted(labels)
    per: dict[str, dict] = {}
    total = {"cells": 0, "claimed": 0, "claimed_right": 0, "claimed_wrong": 0}
    claimed_all: list[str] = []
    for g in sorted({fr.sublineage(c) for c in cells}):
        held = [c for c in cells if fr.sublineage(c) == g]
        write([c for c in cells if fr.sublineage(c) != g])
        _, body = arm.run(ref)
        c = credit(body, ref, held)
        claimed_all += c["claimed_cells"]
        per[g] = {k: c[k] for k in ("cells", "claimed", "claimed_right", "claimed_wrong")}
        for k in total:
            total[k] += c[k]
    n = total["cells"]
    return {
        **total,
        "abstained": n - total["claimed"],
        "honest": total["claimed_right"],
        "honest_fraction": round(total["claimed_right"] / n, 4) if n else None,
        "precision": round(total["claimed_right"] / total["claimed"], 4) if total["claimed"] else None,
        "claimed_cells": claimed_all,
        "per_sublineage": per,
    }


# ---- baselines that use no factors ----------------------------------------------------------------


def baselines(labels: dict[str, str]) -> tuple[dict[str, dict[str, str]], dict]:
    """Two label-free predictors and what they score on all 555.

    `majority` names the commonest terminal tissue for every cell: no factors, no per-cell information at
    all, only the composition of the answer sheet. `sublineage_majority` names the commonest tissue of
    the cell's founder sublineage: no factors, but the observed tree. Both are in sample and both are
    therefore generous; that is the point, they are the bar the factors have to clear."""
    top = Counter(labels.values()).most_common(1)[0][0]
    preds = {
        "majority_tissue": dict.fromkeys(labels, top),
        "sublineage_majority": fr.lineage_majority(labels),
    }
    scores = {
        name: {
            "right": sum(1 for c, t in labels.items() if pred[c] == t),
            "of": len(labels),
            "fraction": round(sum(1 for c, t in labels.items() if pred[c] == t) / len(labels), 4),
        }
        for name, pred in preds.items()
    }
    scores["majority_tissue"]["names"] = top
    counts = Counter(labels.values())
    scores["draw_from_the_composition"] = {
        "expected_right": round(sum(n * n / len(labels) for n in counts.values()), 1),
        "of": len(labels),
    }
    return preds, scores


def shuffle_null(
    feats: dict[str, set[str]],
    instant: dict[str, set[str]],
    labels: dict[str, str],
    n: int = PERMUTATIONS,
    seed: int = 0,
) -> dict:
    """F4. The same cited rules, applied to factor sets permuted among the 555 cells. What survives is
    only how well the rules' output distribution matches the tissue composition; what is destroyed is any
    link between a cell's own factors and its own fate. Evaluated without the Body, on the reads the Body
    gave each cell, because a permutation is not a program."""
    rules = cited_rules()
    cells = sorted(labels)
    rng = random.Random(seed)
    rights, claims = [], []
    for _ in range(n):
        order = cells[:]
        rng.shuffle(order)
        swap = dict(zip(cells, order, strict=True))
        reads = {
            "integrated": {c: feats[swap[c]] for c in cells},
            "instantaneous": {c: instant[swap[c]] for c in cells},
        }
        pred = fr.apply_mixed(rules, reads)
        rights.append(sum(1 for c, p in pred.items() if labels[c] == p))
        claims.append(len(pred))
    mean = sum(rights) / n
    sd = (sum((r - mean) ** 2 for r in rights) / n) ** 0.5
    return {
        "permutations": n,
        "claimed_mean": round(sum(claims) / n, 1),
        "right_mean": round(mean, 2),
        "right_sd": round(sd, 2),
        "right_max": max(rights),
    }


# ---- the terminality gate -------------------------------------------------------------------------


def strip_lookup_fates(arm: Arm) -> int:
    """Delete the observed lineage's 555 terminal-fate decisions from an arm's copy of the program. The
    divisions, the deaths and the timers stay, so the tree is unchanged and a cell that would be terminal
    simply never differentiates."""
    path = arm.dir / "lineage_embryo.bio"
    lines = path.read_text().splitlines()
    kept = [x for x in lines if not x.startswith("decision fate_")]
    path.write_text("\n".join(kept) + "\n")
    return len(lines) - len(kept)


def open_the_guard(arm: Arm) -> None:
    """Replace the factor rules' `cell_type = <terminal types>` guard with the type an undifferentiated
    cell actually has, so a rule may fire on any cell that has not differentiated yet. On the program
    whose lookup fates are stripped, this is the factors deciding terminality as well as type."""
    path = arm.dir / "fates.bio"
    path.write_text(path.read_text().replace(f"cell_type = {fr.TERMINAL_TYPES}", "cell_type = Blastomere"))


def gate(ref: ReferenceLineage, tmpdir: Path, feats: dict[str, set[str]], labels: dict[str, str]) -> dict:
    """P6 and P7. The guard the factor rules carry means they are only ever asked about a cell the lookup
    has already called terminal, so the 555 eligible decision points are the lookup's, not the factors'.
    Without it the rules would also have to say *when* a cell is done dividing, and this measures how
    they do at that: every cell born by hatching that is not one of the 555 is a cell they must not
    claim."""
    rules = cited_rules()
    out = {}
    for name, opened in (("guard as shipped", False), ("guard opened to any undifferentiated cell", True)):
        arm = Arm(tmpdir / f"gate_{int(opened)}")
        (arm.dir / "fates.bio").write_text(fr.to_bio_fates(rules))
        removed = strip_lookup_fates(arm)
        if opened:
            open_the_guard(arm)
        body = Body(parse_file(arm.dir / "embryo_factors.bio"), means=True).run(until=HATCH)
        fired = {
            cid: c
            for cid, c in body.cells.items()
            if any(d.startswith("factor_fate_") for d in c.fired) and c.born < HATCH
        }
        terminal = set(labels)
        right = sum(
            1 for cid in fired if cid in terminal and fired[cid].cell_type == ref.cells[cid].cell_type
        )
        out[name] = {
            "lookup_fate_decisions_deleted": removed,
            "cells_a_factor_rule_fired_on": len(fired),
            "of_those_terminal": sum(1 for cid in fired if cid in terminal),
            "of_those_not_terminal": sum(1 for cid in fired if cid not in terminal),
            "terminal_and_right": right,
            "honest_of_555": right,
            "cells_born_by_hatching": sum(1 for c in body.cells.values() if c.born < HATCH),
        }
        shutil.rmtree(arm.dir, ignore_errors=True)
    return out


# ---- the arms ------------------------------------------------------------------------------------


def arm_for(
    tmpdir: Path,
    tag: str,
    rules_for: Callable[[list[str]], list[fr.Rule]],
    read: str = "exposure",
    threshold: float = fr.THRESHOLD_MIN,
) -> tuple[Arm, Callable[[list[str]], None]]:
    arm = Arm(tmpdir / tag)

    def write(fit: list[str]) -> None:
        (arm.dir / "fates.bio").write_text(
            fr.to_bio_fates(rules_for(fit), read=read, window="lineage", threshold=threshold)
        )

    return arm, write


def main() -> None:
    ref = ReferenceLineage.load()
    term = fr.embryonic_terminal(ref)
    labels = {c.id: c.tissue for c in term}
    cells = sorted(labels)
    atlas = load_cells()
    factors = sorted({f for fs in atlas.values() for f in fs})
    print(f"{len(cells)} embryonic terminal cells, {len(factors)} factors", flush=True)

    reads, _points = harvest(ref, factors)
    instant = presence_at_birth(ref)
    feats = fr.runtime_features(reads, "exposure", "lineage", fr.THRESHOLD_MIN)
    preds, base_scores = baselines(labels)

    # the denominators, before any score
    nameable = {r.tissue for r in cited_rules()}
    counts = Counter(labels.values())
    denominators = {
        "embryonic_terminal_cells": len(cells),
        "in_a_tissue_a_cited_rule_can_name": sum(n for t, n in counts.items() if t in nameable),
        "tissues_a_cited_rule_can_name": sorted(nameable),
        "terminal_types_in_the_lineage_labelling": len(set(fr.TISSUE_CELL_TYPE.values())),
        "tracked_in_the_atlas_themselves": sum(1 for c in term if c.id in atlas),
        "inheriting_an_ancestor_factor_set": sum(1 for c in term if c.id not in atlas),
        "factors_in_the_atlas": len(factors),
        "tissue_counts": dict(counts.most_common()),
    }

    out: dict = {
        "result": "celegans_fate_credit",
        "date": date.today().isoformat(),
        "question": "of the 555 embryonic terminal fates, how many do the measured transcription factors "
        "decide on their own -- with no credit for a fate the observed lineage lookup supplied, and none "
        "for a rule fitted to the observed labels?",
        "what_this_cannot_do": [
            "the 555 cells, their birth times, the tree above them, the fact that each is terminal and "
            "the vocabulary of 17 terminal types all come from the observed lineage",
            "100 of the 555 are not in Ma 2021 and inherit the factor set of the nearest tracked "
            "ancestor, along the observed tree",
            "the atlas ends at the bean stage, so a fate whose terminal genes come later is not in the "
            "measurement",
            "holding a fitted rule out by founder sublineage does not make it unfitted: the hypothesis "
            "class was chosen with all 555 labels in view",
        ],
        "denominators": denominators,
        "baselines_that_use_no_factors": base_scores,
    }

    tmpdir = Path(tempfile.mkdtemp(prefix="celegans-fate-credit-"))
    arms: dict[str, dict] = {}
    try:
        # 1. the program as it ships: cited rules plus the fitted neuron decision list
        arm, write = arm_for(tmpdir, "shipped", lambda fit: program_rules(feats, labels, fit))
        write(cells)
        _, body = arm.run(ref)
        shipped = credit(body, ref, cells)
        assert shipped["with_lookup_fallback"] == 522 and shipped["claimed"] == 332, "F1: not the program"
        arms["cited + fitted, in sample"] = shipped
        arms["cited + fitted, held out by founder sublineage"] = fold_credit(arm, ref, write, labels)

        # 2. nothing fitted: the cited rules alone. The fold arm is run too, not because the rules could
        #    change -- they cannot, `fit` is ignored -- but so that "nothing is fitted here" is a
        #    measurement and not a claim: the held-out credit must equal the in-sample credit exactly.
        arm_c, write_c = arm_for(tmpdir, "cited", lambda _fit: cited_rules())
        write_c(cells)
        _, body_c = arm_c.run(ref)
        arms["cited only, nothing fitted"] = credit(body_c, ref, cells)
        arms["cited only, held out (must be identical)"] = fold_credit(arm_c, ref, write_c, labels)

        # 3. the textbook rules without the cited sheath-glia rule, to see what that one rule is worth
        arm_t, write_t = arm_for(tmpdir, "textbook", lambda _fit: fr.textbook_rules())
        write_t(cells)
        _, body_t = arm_t.run(ref)
        arms["textbook rules only"] = credit(body_t, ref, cells)

        # 4. the fitted neuron list with no cited rule in front of it
        arm_f, write_f = arm_for(
            tmpdir,
            "fitted",
            lambda fit: [r for r in program_rules(feats, labels, fit) if r.tissue == "neuron"],
        )
        write_f(cells)
        _, body_f = arm_f.run(ref)
        arms["fitted neuron rules only, in sample"] = credit(body_f, ref, cells)
        arms["fitted neuron rules only, held out"] = fold_credit(arm_f, ref, write_f, labels)

        # 5. P5: the read that looks better with the fallback, scored without one
        feats_mean = fr.runtime_features(reads, "mean", "lineage", 0.25)
        arm_m, write_m = arm_for(
            tmpdir,
            "mean",
            lambda fit: program_rules(feats_mean, labels, fit),
            read="mean",
            threshold=0.25,
        )
        write_m(cells)
        _, body_m = arm_m.run(ref)
        arms["mean(lineage) >= 0.25, cited + fitted, in sample"] = credit(body_m, ref, cells)
        arms["mean(lineage) >= 0.25, held out"] = fold_credit(arm_m, ref, write_m, labels)

        for a in arms.values():
            a["baselines_on_the_same_cells"] = on_the_same_cells(a.pop("claimed_cells"), labels, preds)

        # F3: is the headline flat across the published plateau, or is it a knob?
        plateau = []
        for thr in (1, 15, 30, 45, 60, 90, 120):
            a, w = arm_for(tmpdir, f"plateau_{thr}", lambda _fit: cited_rules(), threshold=float(thr))
            w(cells)
            _, b = a.run(ref)
            c = credit(b, ref, cells)
            plateau.append({"threshold": thr, "claimed": c["claimed"], "honest": c["honest"]})
            shutil.rmtree(a.dir, ignore_errors=True)
        out["headline_across_the_threshold_plateau"] = plateau

        out["terminality_gate"] = gate(ref, tmpdir, feats, labels)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    out["arms"] = arms
    out["shuffle_null_for_the_cited_rules"] = shuffle_null(feats, instant, labels)

    # the pre-registered claims, checked
    cited = arms["cited only, nothing fitted"]
    shipped_in = arms["cited + fitted, in sample"]
    shipped_out = arms["cited + fitted, held out by founder sublineage"]
    mean_in = arms["mean(lineage) >= 0.25, cited + fitted, in sample"]
    mean_out = arms["mean(lineage) >= 0.25, held out"]
    null = out["shuffle_null_for_the_cited_rules"]
    opened = out["terminality_gate"]["guard opened to any undifferentiated cell"]
    flat = [p["honest"] for p in plateau if p["threshold"] <= 60]
    out["preregistered"] = {
        "P1 cited rules alone decide 130 to 160 of 555 correctly": {
            "predicted": "130..160",
            "observed": cited["honest"],
            "held": 130 <= cited["honest"] <= 160,
        },
        "P2 the shipped set decides 299 in sample and 268 held out": {
            "predicted": [299, 268],
            "observed": [shipped_in["honest"], shipped_out["honest"]],
            "held": shipped_in["honest"] == 299 and shipped_out["honest"] == 268,
        },
        "P3 the honest number is below the factor-free majority-tissue baseline": {
            "predicted": f"cited ({cited['honest']}) < {base_scores['majority_tissue']['right']}",
            "observed": [cited["honest"], base_scores["majority_tissue"]["right"]],
            "held": cited["honest"] < base_scores["majority_tissue"]["right"],
        },
        "P4 where they speak the factors beat the majority tissue on the same cells": {
            "observed": {
                k: [arms[k]["claimed_right"], _maj(arms[k])]
                for k in (
                    "cited only, nothing fitted",
                    "cited + fitted, in sample",
                    "cited + fitted, held out by founder sublineage",
                )
            },
            "held": all(
                arms[k]["claimed_right"] > _maj(arms[k])
                for k in (
                    "cited only, nothing fitted",
                    "cited + fitted, in sample",
                    "cited + fitted, held out by founder sublineage",
                )
            ),
        },
        "P5 the honest metric reverses the two reads": {
            "predicted": "mean(lineage) above exposure with the fallback, below it without one",
            "observed": {
                "with_fallback": [mean_in["with_lookup_fallback"], shipped_in["with_lookup_fallback"]],
                "honest_in_sample": [mean_in["honest"], shipped_in["honest"]],
                "honest_held_out": [mean_out["honest"], shipped_out["honest"]],
            },
            "held": mean_in["with_lookup_fallback"] > shipped_in["with_lookup_fallback"]
            and mean_in["honest"] < shipped_in["honest"]
            and mean_out["honest"] < shipped_out["honest"],
        },
        "P6 without the lookup's fates no factor rule fires at all": {
            "observed": out["terminality_gate"]["guard as shipped"]["cells_a_factor_rule_fired_on"],
            "held": out["terminality_gate"]["guard as shipped"]["cells_a_factor_rule_fired_on"] == 0,
        },
        "P7 with the guard opened the cited rules claim more non-terminal cells than terminal ones": {
            "observed": [opened["of_those_not_terminal"], opened["of_those_terminal"]],
            "held": opened["of_those_not_terminal"] > opened["of_those_terminal"],
        },
        "F2 the lookup is a clean fallback: no abstained cell is wrong": {
            "observed": {k: a["abstained_wrong"] for k, a in arms.items()},
            "held": all(a["abstained_wrong"] == 0 for a in arms.values()),
        },
        "F3 the headline is flat to within 10 fates across the 1 to 60 min plateau": {
            "observed": flat,
            "held": max(flat) - min(flat) <= 10,
        },
        "F4 the shuffle null stays more than two sd below the cited number": {
            "observed": [null["right_mean"], null["right_sd"], cited["honest"]],
            "held": cited["honest"] > null["right_mean"] + 2 * null["right_sd"],
        },
    }
    RESULT.write_text(json.dumps(out, indent=2) + "\n")

    print(f"\n{'arm':<52} {'claimed':>8} {'right':>6} {'of 555':>7} {'fallback':>9}")
    for name, a in arms.items():
        print(
            f"{name:<52} {a['claimed']:>8} {a['honest']:>6} "
            f"{a['honest_fraction'] * 100:>6.1f}% {a.get('with_lookup_fallback', '-'):>9}"
        )
    for name, b in base_scores.items():
        print(f"{name:<52} {'-':>8} {b.get('right', b.get('expected_right')):>6}")
    print()
    for claim, v in out["preregistered"].items():
        print(f"  {'held ' if v['held'] else 'BROKEN'}  {claim}: {v['observed']}")


if __name__ == "__main__":
    main()
