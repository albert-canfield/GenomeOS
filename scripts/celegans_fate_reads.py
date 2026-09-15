# SPDX-License-Identifier: AGPL-3.0-or-later
"""The worm's terminal fates on the runtime's own integrated reads (area E, Next step 1).

Until 2026-09-15 the factor rules in `fates.bio` read `ELT-2_integrated = present`, a per-cell lookup
generated into `exposure.bio` by `scripts/celegans_fates.py`. That lookup called a factor present at 20%
of its largest path exposure *over the finished run*, which is a statistic no deciding cell can hold. The
Body now computes `F.exposure(window)` and `F.mean(window)` itself (BIOLANG-v0.4-ECONOMY.md §7.2a), in
absolute units, so a rule can state a threshold a cell could actually read: minutes carried along a named
window, or a fraction of it.

This script rewrites the rules on that read and scores the result honestly:

  1. it runs the lookup-only worm (`embryo.bio`) through the Body once with every read declared, and takes
     the value the runtime put in each embryonic terminal cell's context **at the decision point where its
     fate is taken** — not at the end of the run, and not from the atlas;
  2. it sweeps threshold x read x window on the cited textbook rules, so the threshold is shown to be a
     plateau rather than a fitted knob;
  3. it learns the neuron rules, writes `fates.bio`, and scores the program through the Body and the
     existing LineageDiff three ways: in sample, held out one founder sublineage at a time (the rules
     never see the lineage they are scored on), and against the precomputed baseline regenerated and run
     through the same Body so the comparison is like for like.

    uv run python scripts/celegans_fate_reads.py

No model calls, no network, about three minutes on one core. Result: data/results/celegans_fate_reads.json.
"""

from __future__ import annotations

import json
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
from genomeos.organism import atlas_levels  # noqa: E402
from genomeos.organism import fate_rules as fr  # noqa: E402
from genomeos.organism.diff import compare  # noqa: E402
from genomeos.organism.reference import ReferenceLineage  # noqa: E402
from genomeos.organism.tf_atlas import load_cells  # noqa: E402
from genomeos.runtime.body import Body  # noqa: E402

ORG = Path("data/organisms/celegans")
RESULT = Path("data/results/celegans_fate_reads.json")
HATCH = 800.0
ADULT = 6000.0

# Every factor any rule in the program reads, plus the textbook and glia factors, so the harvest covers
# the hypothesis class the learner searches. Collected from the atlas rather than guessed.
WINDOWS = ("cell", "lineage")
KINDS = ("exposure", "mean")


# ---- 1. what the runtime read at each terminal cell's decision point ----------------------------


def harvest(ref: ReferenceLineage, factors: list[str]) -> tuple[dict[str, dict[str, float]], Counter]:
    """Run the lookup-only worm and take every declared read from the context of each embryonic terminal
    cell at the decision point where its fate is settled. `_reads` is set directly because the reads a
    program pays for are the ones its own rules name, and here the rules do not exist yet. The counter
    says how many decision points each terminal cell gets, which is the whole story for the `cell`
    window: a cell that decides once, at birth, has lived for nothing when it decides."""
    wanted = {c.id for c in fr.embryonic_terminal(ref)}
    body = Body(parse_file(ORG / "embryo.bio"), means=True)
    body._reads = sorted(  # noqa: SLF001 (measuring the runtime's own read, before any rule names it)
        (f"{f}.{k}({w})", f, k, w) for f in factors for k in KINDS for w in WINDOWS
    )
    out: dict[str, dict[str, float]] = {}
    points: Counter[str] = Counter()
    original = Body._resolve

    def resolve(self: Body, c, born: bool = True) -> None:  # noqa: ANN001
        original(self, c, born=born)
        if c.name in wanted:
            points[c.name] += 1
            ctx = self.context(c)
            out.setdefault(c.name, {k: float(ctx[k]) for k, *_ in self._reads if k in ctx})  # noqa: SLF001

    Body._resolve = resolve
    try:
        body.run(until=HATCH)
    finally:
        Body._resolve = original
    return out, Counter(points.values())


def presence_at_birth(ref: ReferenceLineage) -> dict[str, set[str]]:
    """The instantaneous read as the Body gives it: the factors a terminal cell carries when it decides."""
    wanted = {c.id for c in fr.embryonic_terminal(ref)}
    body = Body(parse_file(ORG / "embryo.bio"), means=True)
    out: dict[str, set[str]] = {}
    original = Body._resolve

    def resolve(self: Body, c, born: bool = True) -> None:  # noqa: ANN001
        original(self, c, born=born)
        if c.name in wanted:
            out.setdefault(c.name, set(c.factors))

    Body._resolve = resolve
    try:
        body.run(until=HATCH)
    finally:
        Body._resolve = original
    return out


# ---- 2. the rule list ---------------------------------------------------------------------------


def program_rules(feats: dict[str, set[str]], labels: dict[str, str], cells: list[str]) -> list[fr.Rule]:
    """The same shape as the program has carried since 2026-09-14: the cited textbook rules for intestine
    and body-wall muscle, the cited sheath-glia rule (instantaneous), the cited hypodermis rules, then the
    neuron rules learned from the atlas. Only the neuron rules are fitted, so only they need holding out."""
    textbook = fr.textbook_rules()
    train = {c: feats[c] for c in cells}
    neuron = [r for r in fr.learn(train, {c: labels[c] for c in cells}) if r.tissue == "neuron"]
    return textbook[:4] + fr.glia_rules() + textbook[4:] + neuron


# ---- 3. scoring through the Body and the existing diff -------------------------------------------


class Arm:
    """A worm program in a scratch copy of data/organisms/celegans, so `fates.bio` can be rewritten per
    fold without touching the checkout."""

    def __init__(self, tmp: Path, baseline: bool = False) -> None:
        self.dir = tmp
        shutil.copytree(ORG, self.dir, dirs_exist_ok=True)
        text = (ORG / "embryo_factors.bio").read_text()
        if baseline and "import exposure.bio" not in text:
            text = text.replace("import fates.bio", "import exposure.bio\nimport fates.bio")
        if not baseline:
            text = text.replace("import exposure.bio\n", "")
            (self.dir / "exposure.bio").unlink(missing_ok=True)
        (self.dir / "embryo_factors.bio").write_text(text)

    def run(self, ref: ReferenceLineage, until: float = HATCH) -> tuple:
        body = Body(parse_file(self.dir / "embryo_factors.bio"), means=True).run(until=until)
        return compare(body, ref, until=until), body


def decided(body: Body, cells: list[str]) -> int:
    return sum(
        1
        for c in cells
        if (x := body.cells.get(c)) is not None and any(d.startswith("factor_fate_") for d in x.fired)
    )


def score_arm(
    arm: Arm,
    ref: ReferenceLineage,
    write: Callable[[list[str]], None],
    labels: dict[str, str],
    folds: bool = True,
) -> dict:
    """`write(cells)` puts a fates.bio (and, for the baseline, an exposure.bio) fitted on `cells` into the
    arm's directory. In sample it is fitted on all 555; held out, eight times on seven founder sublineages
    and scored only on the eighth, so no rule is credited on a lineage it saw."""
    cells = sorted(labels)
    write(cells)
    d, body = arm.run(ref)
    out = {
        "in_sample": {
            "fates_correct": d.fates_correct,
            "fates_checked": d.fates_checked,
            "decided_by_factors": decided(body, cells),
            "cells_born": d.matched,
            "deaths_matched": d.deaths_matched,
        }
    }
    d_adult, _ = arm.run(ref, until=ADULT)
    out["in_sample"]["to_adult"] = [d_adult.fates_correct, d_adult.fates_checked]
    out["in_sample"]["per_tissue"] = per_tissue(body, ref, cells)
    out["in_sample"]["confusion"] = confusion(body, ref, cells)
    if not folds:
        return out
    right = checked = dec = 0
    for g in sorted({fr.sublineage(c) for c in cells}):
        held = [c for c in cells if fr.sublineage(c) == g]
        write([c for c in cells if fr.sublineage(c) != g])
        _, b = arm.run(ref)
        for cid in held:
            checked += 1
            bc = b.cells.get(cid)
            right += bc is not None and bc.cell_type == ref.cells[cid].cell_type
            dec += bc is not None and any(x.startswith("factor_fate_") for x in bc.fired)
    out["held_out"] = {"fates_correct": right, "fates_checked": checked, "decided_by_factors": dec}
    return out


def per_tissue(body: Body, ref: ReferenceLineage, cells: list[str]) -> dict:
    per: dict[str, list[int]] = {}
    for cid in cells:
        rc = ref.cells[cid]
        bc = body.cells.get(cid)
        row = per.setdefault(rc.tissue, [0, 0])
        row[1] += 1
        row[0] += bc is not None and bc.cell_type == rc.cell_type
    return {t: {"correct": k, "cells": n} for t, (k, n) in sorted(per.items(), key=lambda kv: -kv[1][1])}


def confusion(body: Body, ref: ReferenceLineage, cells: list[str]) -> dict:
    c: Counter[str] = Counter()
    for cid in cells:
        rc, bc = ref.cells[cid], body.cells.get(cid)
        if bc is not None and bc.cell_type != rc.cell_type:
            c[f"{rc.tissue}->{bc.cell_type}"] += 1
    return dict(c.most_common(10))


# ---- the threshold plateau, on the cited rules only, without the Body ----------------------------


def plateau(reads: dict, labels: dict[str, str], instant: dict[str, set[str]]) -> list[dict]:
    """The cited textbook rules alone (no fitting anywhere), applied to the runtime read at a grid of
    thresholds, with the lineage lookup as the fallback: `wrong` is what the rules get wrong and
    `decided` is how many cells they claim. The instantaneous row has no threshold of its own — the 20%
    presence call is inside the reader."""
    rules = fr.textbook_rules()
    rows = []
    grids = {
        ("exposure", "lineage"): (1, 15, 30, 45, 60, 75, 90, 105, 120, 150, 180),
        ("mean", "lineage"): (0.01, 0.05, 0.1, 0.15, 0.2, 0.22, 0.25, 0.3, 0.4, 0.5),
        ("exposure", "cell"): (1, 15, 30),
        ("mean", "cell"): (0.01, 0.1, 0.2),
    }
    for (kind, window), grid in grids.items():
        for thr in grid:
            feats = fr.runtime_features(reads, kind, window, thr)
            pred = fr.apply_rules(rules, feats)
            wrong = sum(1 for c, p in pred.items() if labels[c] != p)
            rows.append(
                {
                    "read": f"{kind}({window})",
                    "threshold": thr,
                    "decided": len(pred),
                    "right": len(pred) - wrong,
                    "wrong": wrong,
                    "with_lookup_fallback": len(labels) - wrong,
                }
            )
    pred = fr.apply_rules(rules, instant)
    wrong = sum(1 for c, p in pred.items() if labels[c] != p)
    rows.append(
        {
            "read": "instantaneous (F = present)",
            "threshold": None,
            "decided": len(pred),
            "right": len(pred) - wrong,
            "wrong": wrong,
            "with_lookup_fallback": len(labels) - wrong,
        }
    )
    return rows


def cadence_stability(ref: ReferenceLineage, tmpdir: Path, reads: dict, labels: dict[str, str]) -> dict:
    """A fate rule on an integrated read is not stable when the cell decides again, and this measures how
    far. `F.exposure(lineage)` grows for as long as the cell carries the factor, so a guard that was false
    when the fate was settled becomes true later for no reason but the clock. The Body refuses a
    lower-precedence decision at a later decision point only when it *could already have applied* then
    (§7.3), and a threshold on a growing integral is never in that set. With a `cell_network` cadence the
    worm re-decides every live cell every step, and the score falls. The lookup this rewrite replaces was
    a time-invariant guard and did not move, which is why the engine's cadence invariant held for it."""
    cells = sorted(labels)
    rows = {}
    for kind, thr in (("exposure", fr.THRESHOLD_MIN), ("mean", 0.25)):
        arm = Arm(tmpdir / f"cadence_{kind}")
        feats = fr.runtime_features(reads, kind, "lineage", thr)
        (arm.dir / "fates.bio").write_text(
            fr.to_bio_fates(program_rules(feats, labels, cells), read=kind, window="lineage", threshold=thr)
        )
        for cadence in (0, 6):
            module = parse_file(arm.dir / "embryo_factors.bio")
            module.organism.cell_network = cadence
            body = Body(module, seed=None).run(until=HATCH)
            d = compare(body, ref, until=HATCH)
            s = body.summary()
            rows[f"{kind}(lineage) >= {thr:g}, cadence {cadence} min"] = {
                "fates_correct": d.fates_correct,
                "fates_checked": d.fates_checked,
                "overruled_fates": s["overruled_fates"],
                "revised_fates": s["revised_fates"],
                "network_steps": s["network_steps"],
            }
        shutil.rmtree(arm.dir, ignore_errors=True)
    return rows


# ---- main -----------------------------------------------------------------------------------------


def main() -> None:
    ref = ReferenceLineage.load()
    term = fr.embryonic_terminal(ref)
    labels = {c.id: c.tissue for c in term}
    cells = sorted(labels)
    atlas = load_cells()
    factors = sorted({f for fs in atlas.values() for f in fs})
    print(f"{len(cells)} embryonic terminal cells, {len(factors)} factors in the atlas", flush=True)

    reads, points = harvest(ref, factors)
    instant = presence_at_birth(ref)
    print(f"runtime reads harvested; decision points per terminal cell: {dict(points)}", flush=True)
    empty_cell_window = sum(
        1 for r in reads.values() if not any(v > 0 for k, v in r.items() if "(cell)" in k)
    )

    out: dict = {
        "result": "celegans_fate_reads",
        "date": date.today().isoformat(),
        "question": "with the Body computing F.exposure(window) and F.mean(window) itself, can the worm's "
        "terminal fate rules be written on a read a deciding cell could hold, and does that score better "
        "or worse than the precomputed lookup it replaces?",
        "reads_are": "absolute (minutes, or a fraction of the named window) and integrate carriage, not "
        "concentration: this program declares no cell_network, so a factor counts 1 while carried "
        "(BIOLANG-v0.4-ECONOMY.md §7.2a)",
        "decision_points_per_terminal_cell": {str(k): v for k, v in sorted(points.items())},
        "terminal_cells_whose_cell_window_is_empty_when_they_decide": empty_cell_window,
        "threshold_plateau": plateau(reads, labels, instant),
    }

    tmpdir = Path(tempfile.mkdtemp(prefix="celegans-fate-reads-"))
    try:
        arms: dict[str, dict] = {}

        # the read the program adopts, and the other three windows for the comparison
        for kind, window, thr in (
            ("exposure", "lineage", fr.THRESHOLD_MIN),
            ("mean", "lineage", 0.25),
            ("exposure", "cell", 1.0),
            ("mean", "cell", 0.01),
        ):
            arm = Arm(tmpdir / f"{kind}_{window}_{thr:g}")
            feats = fr.runtime_features(reads, kind, window, thr)

            def write(  # noqa: E501
                fit: list[str],
                arm: Arm = arm,
                feats: dict = feats,
                k: str = kind,
                w: str = window,
                t: float = thr,
            ) -> None:
                (arm.dir / "fates.bio").write_text(
                    fr.to_bio_fates(program_rules(feats, labels, fit), read=k, window=w, threshold=t)
                )

            arms[f"{kind}({window}) >= {thr:g}"] = score_arm(arm, ref, write, labels)
            print(f"  {kind}({window}) >= {thr:g}: {arms[f'{kind}({window}) >= {thr:g}']}", flush=True)

        # the instantaneous read, through the same Body and the same rule shape
        arm = Arm(tmpdir / "instantaneous")

        def write_instant(fit: list[str], arm: Arm = arm) -> None:
            rules = [
                fr.Rule(r.tissue, r.factors, r.source, r.precision, r.support, read="instantaneous")
                for r in program_rules(instant, labels, fit)
            ]
            (arm.dir / "fates.bio").write_text(fr.to_bio_fates(rules))

        arms["instantaneous (F = present)"] = score_arm(arm, ref, write_instant, labels)

        # the baseline: the precomputed lookup this rewrite replaces, regenerated and run the same way
        levels = atlas_levels.load_levels()
        exp = fr.exposures(ref, levels, term)
        integ = fr.integrated(exp)
        inst_atlas = fr.instantaneous(ref, atlas, term)
        base = Arm(tmpdir / "baseline", baseline=True)

        def write_baseline(fit: list[str]) -> None:
            rules = program_rules(integ, labels, fit)
            (base.dir / "fates.bio").write_text(_legacy_fates(rules))
            wanted = {f for r in rules if r.read == "integrated" for f in r.factors}
            (base.dir / "exposure.bio").write_text(fr.to_bio_exposure(integ, inst_atlas, atlas, wanted))

        arms["baseline: precomputed _integrated lookup"] = score_arm(base, ref, write_baseline, labels)
        out["arms"] = arms

        # why the threshold is 15 min and not the value that scores best: the whole program, through the
        # Body, at a grid of thresholds. Below 60 min nothing moves at all; above it the score rises only
        # because the rules claim fewer cells, which is buying a number with abstention.
        curve = []
        for kind, grid in (("exposure", (1, 15, 30, 60, 90, 120, 180)), ("mean", (0.1, 0.2, 0.25, 0.4))):
            for thr in grid:
                a = Arm(tmpdir / f"curve_{kind}_{thr:g}")
                f2 = fr.runtime_features(reads, kind, "lineage", thr)

                def w2(fit: list[str], a: Arm = a, f2: dict = f2, k: str = kind, t: float = thr) -> None:
                    (a.dir / "fates.bio").write_text(
                        fr.to_bio_fates(program_rules(f2, labels, fit), read=k, window="lineage", threshold=t)
                    )

                s = score_arm(a, ref, w2, labels)
                curve.append(
                    {
                        "read": f"{kind}(lineage)",
                        "threshold": thr,
                        "in_sample": s["in_sample"]["fates_correct"],
                        "in_sample_decided": s["in_sample"]["decided_by_factors"],
                        "held_out": s["held_out"]["fates_correct"],
                        "held_out_decided": s["held_out"]["decided_by_factors"],
                    }
                )
                shutil.rmtree(a.dir, ignore_errors=True)
                print(f"  curve {kind}({thr:g}): {curve[-1]}", flush=True)
        out["threshold_curve_through_the_body"] = curve
        out["stability_under_a_re_decision_cadence"] = cadence_stability(ref, tmpdir, reads, labels)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    # what the program ships
    feats = fr.runtime_features(reads, "exposure", "lineage", fr.THRESHOLD_MIN)
    rules = program_rules(feats, labels, cells)
    (ORG / "fates.bio").write_text(fr.to_bio_fates(rules, read="exposure", window="lineage"))
    out["program_rules"] = [
        {
            "tissue": r.tissue,
            "factors": list(r.factors),
            "read": r.read,
            "precision": r.precision,
            "support": r.support,
            "source": r.source,
        }
        for r in rules
    ]
    out["threshold"] = {"minutes": fr.THRESHOLD_MIN, "basis": fr.THRESHOLD_BASIS}
    out["generated"] = [str(ORG / "fates.bio")]
    RESULT.write_text(json.dumps(out, indent=2) + "\n")
    for name, a in out["arms"].items():
        i, h = a["in_sample"], a.get("held_out", {})
        print(
            f"{name:<44} in sample {i['fates_correct']}/{i['fates_checked']} "
            f"(decided {i['decided_by_factors']}), to adult {i['to_adult'][0]}/{i['to_adult'][1]}, "
            f"held out {h.get('fates_correct')}/{h.get('fates_checked')} "
            f"(decided {h.get('decided_by_factors')})"
        )


def _legacy_fates(rules: list[fr.Rule]) -> str:
    """The pre-2026-09-15 emitter, inlined here because only the baseline arm needs it: every integrated
    rule reads a `_integrated` factor set by the generated exposure.bio lookup."""
    lines = [
        "# The precomputed lookup this rewrite replaces; regenerated only to score the old read.",
        "module organism.celegans.fates",
        "import cell_types.bio",
        "",
    ]
    for i, r in enumerate(rules):
        suffix = fr.INTEGRATED_SUFFIX if r.read == "integrated" else ""
        cond = ", ".join(f"{f}{suffix} = present" for f in r.factors)
        lines.append(
            f"decision factor_fate_{i:02d} {{ action: differentiate; "
            f"when: cell_type = {fr.TERMINAL_TYPES}, stage = {fr.EMBRYONIC_STAGES}, {cond}; "
            f"to: {fr.TISSUE_CELL_TYPE[r.tissue]}; priority: {len(rules) - i}; "
            'evidence: inferred "the precomputed read, kept for the baseline arm"; confidence: 0.7 }'
        )
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    main()
