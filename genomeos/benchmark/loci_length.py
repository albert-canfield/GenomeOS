# SPDX-License-Identifier: AGPL-3.0-or-later
"""Is "a value sits on constrained sequence" an artefact of element length?

The known-locus benchmark reads one claim as separating its panel from its matched controls:
**a value sits on constrained sequence**, 9 of 17 at the panel against 15 of 85 at the controls.
Every other claim it scores has been withdrawn as an artefact of coverage or of matching. The
suspicion this module tests, registered before anything was computed, is that the survivor is an
artefact of **element length**: a longer window is more likely to contain a variable position that
happens to sit on constrained sequence, for reasons that have nothing to do with being a regulatory
element.

The module builds a third control arm, prints the length distribution of every arm so a reader can
see whether the matching worked, and runs the comparison through `genomeos.compare.standardised`
both ways: holding the covariates fixed, and holding coverage fixed beside them.

`PREREGISTRATION` below is the whole of the registration and was committed before the first number
existed; `verdict()` applies its decision rule mechanically, so the rule cannot be re-read after the
number arrives.
"""

from __future__ import annotations

import bisect
import random
import time
from pathlib import Path
from statistics import median
from typing import Any

from genomeos.benchmark.loci import (
    NEGATIVE_KEEP_OUT,
    PANEL,
    Chromosome,
    Expect,
    Window,
    read_constraint,
    read_deletion,
    read_syntax_values_many,
)
from genomeos.compare import Strata, standardised
from genomeos.results import RESULTS_DIR, load_result, save_result

# ----------------------------------------------------------------- the registration, in code

#: How a control window is drawn for the new arm. Exact, because `candidate_windows` in `loci.py`
#: already draws a window of exactly the positive element's length and this arm keeps that.
LENGTH_TOLERANCE_BP = 0
CONTROLS_PER_LOCUS = 5
DRAW_SEED = 4703  # not loci.py's 11: the new arm must not be the old arm under another name
DRAW_ATTEMPTS = 4_000  # per locus, before giving up on a locus with a crowded chromosome
VISTA_FLANK = 5_000
GWAS_FLANK = 5_000

#: The claim. A window scores a hit when some position where one of the people we hold differs from
#: the reference sits above the Zoonomia phyloP threshold (`loci.read_syntax_values_many`).
HIT = "storage"

#: Registered strata. `length` is the covariate under test; `deletion_scored` is coverage as a
#: stratum; `has_values` separates a window that was looked at and holds no variable position from
#: one that holds several and none of them is constrained.
STRATA: dict[str, Strata] = {
    # primary, covariates only
    "P1_length": Strata(length=(1_000, 5_000)),
    # primary, coverage held fixed beside it
    "P2_length_coverage": Strata(length=(1_000, 5_000), deletion_scored=(0.5,)),
    # primary, the claim's own denominator held fixed: does the window hold a value at all
    "P3_length_has_values": Strata(length=(1_000, 5_000), has_values=(0.5,)),
    # secondary, the whole covariate set the existing controls were matched on
    "S1_covariates": Strata(length=(1_000, 5_000), gc=(0.45,), distance=(25_000,)),
    # secondary, that set with coverage added
    "S2_covariates_coverage": Strata(
        length=(1_000, 5_000), gc=(0.45,), distance=(25_000,), deletion_scored=(0.5,)
    ),
}
#: Which stratification decides the verdict, per control arm. Arm B is matched on GC and distance
#: already, so length is the covariate left to hold fixed. Arm C is matched on NOTHING but length, so
#: a low rate there could be its distance to a coding TSS rather than its length, and the covariate
#: set has to be standardised away before the verdict reads it. Revision of 2026-09-17, before
#: scoring; see THE REVISION in `PREREGISTRATION`.
PRIMARY_BY_ARM = {"B_existing": "P1_length", "C_length_only": "S1_covariates"}
PRIMARY_WITH_COVERAGE_BY_ARM = {
    "B_existing": "P2_length_coverage",
    "C_length_only": "S2_covariates_coverage",
}
PRIMARY = "P1_length"  # the default reading, and arm B's
PRIMARY_WITH_COVERAGE = "P2_length_coverage"

SURVIVES_DIFFERENCE = 0.20
SURVIVES_P = 0.05
SURVIVES_WITH_COVERAGE = 0.15
ARTEFACT_DIFFERENCE = 0.10
ARTEFACT_RAW_FLOOR = 0.25
MIN_TARGETS_MATCHED = 14
MAX_DROPPED = 3

PREREGISTRATION = """\
Registered 2026-09-17, in this file, before any number in it existed. Revised 2026-09-17 after an
outside read and BEFORE scoring; the revision is written out under THE REVISION below rather than
folded silently into the text. Read together with the limitations at the end, which are part of the
registration and not a caveat added afterwards.

THIS LANE IS NOT BLINDED, and the reader meets that here rather than after the tables. Before the
blind was widened to the whole of docs/LOCI-BENCHMARK.md, this lane had read a paragraph of section 9
that quotes the later section's results, so it knows a second lane's figures for this claim and that
lane's own remark about element length. The blind failed on timing, not on scope, and no widening
could un-read it. What survives is worth stating and is not nothing: this registration was committed
before any number in it existed, its decision rule is mechanical and symmetric between the outcomes,
and its control selection is seeded and executes rather than being described. That is a
pre-registered unblinded test -- weaker than a blinded one, much stronger than an unregistered read.

THE QUESTION. The known-locus benchmark reports "a value sits on constrained sequence" at 9 of 17
published loci against 15 of 85 matched control windows. Every other claim it scores has been
withdrawn. This lane asks whether the survivor is an artefact of element length.

WHAT IS FIXED AND BY WHOM. The seventeen loci are the existing panel and nobody is choosing them
now; they were fixed before this lane existed, in `loci.PANEL`, with their citations. This lane
chooses CONTROL WINDOWS only. The selection rule for them is the code below, so it executes rather
than being described.

THE ARMS.
  A  panel        the 17 published elements, read over `Expect.element`.
  B  existing     the 85 control windows already in `data/results/loci_benchmark.json`, 5 per locus,
                  drawn by `loci.candidate_windows` and narrowed by `loci.pick_negatives`: a window
                  of exactly the element's length, GC within 0.04, distance to the nearest coding
                  TSS within 35%, constrained fraction the closest of the candidates.
  C  length-only  NEW. 5 windows per locus, drawn here, of exactly the element's length (tolerance
                  0 bp), on the same chromosome, outside every panel locus with 200 kb of flank,
                  outside VISTA elements and GWAS Catalog hits with 5 kb of flank -- and matched on
                  NOTHING ELSE. No GC, no distance, no constraint. Seed 4703, fixed here.

WHY C IS DRAWN THAT WAY, AND THE CONSTRUCTION FACT THAT ANSWERED THE QUESTION BEFORE IT WAS ASKED.
The brief asked for a third arm matched on length IN ADDITION to what the existing controls match
on. Reading `candidate_windows` before registering shows that arm B is ALREADY matched on length
exactly -- `length = element[1] - element[0]` and the drawn window is `s .. s + length` -- so "B plus
length" is B, and a third arm built that way would be the second arm under another name. THE
EXISTING COMPARISON THEREFORE ALREADY HOLDS ELEMENT LENGTH FIXED, WINDOW BY WINDOW, AT ZERO
TOLERANCE, and 9 of 17 against 15 of 85 cannot be a naive length artefact -- not because a test
looked and found none, but because the design makes one impossible. That retires the question as the
brief posed it, and `length_match_check` verifies it from the committed coordinates rather than from
a reading of the code.

WHAT ARM C IS FOR, THEN. Not "does length explain the excess", which is answered above, but the
COMPLEMENT question: what does matching on length alone buy, with everything else left unmatched? It
bounds how much of the panel's rate is reachable by length and nothing else. Nobody should read it as
the original test.

THE STATISTIC. `genomeos.compare.standardised(targets, controls, strata, hit="storage")`, direct
standardisation, the control arm carrying the matched n rather than the number of control rows, with
`dropped_for_want_of_a_control` and both groups' covariate medians reported. Not reimplemented here.

THE DENOMINATORS. Arm A 17, arm B 85 (5 per locus), arm C 85 (5 per locus, fewer only if a locus's
chromosome cannot supply five windows outside the keep-outs, which is reported per locus). The
standardised comparison's n is `targets_matched`, at most 17, and that is the n that sets its
precision -- not 85.

THE STRATA, both ways, side by side:
  P1_length             length (cuts 1 kb, 5 kb)                              covariates only
  P2_length_coverage    P1 plus deletion_scored (cut 0.5)                     coverage held fixed
  P3_length_has_values  P1 plus has_values (cut 0.5)                          the claim's denominator
  S1_covariates         length, GC (0.45), distance (25 kb)                   covariates only
  S2_covariates_coverage  S1 plus deletion_scored                             coverage held fixed
All five are computed and printed for both control arms. Which one DECIDES differs by arm, because
the arms are matched on different things:
  arm B decides on P1_length and P2_length_coverage      (already matched on GC and distance)
  arm C decides on S1_covariates and S2_covariates_coverage (matched on nothing but length)

THE OUTCOMES. Write prim(X) for the deciding covariates-only stratification of arm X and cov(X) for
the same with coverage added, per the two lines above; raw(X) is the unstandardised difference.
  "the claim survives"          d(prim(B)) >= 0.20 at one-sided p <= 0.05 AND d(prim(C)) >= 0.20 at
                                one-sided p <= 0.05 AND d(cov(B)) >= 0.15 AND d(cov(C)) >= 0.15, with
                                targets_matched >= 14 and dropped <= 3 in each.
  "the claim is an artefact of length"
                                for either control arm X: raw(X) >= 0.25 while d(prim(X)) <= 0.10 --
                                the excess exists unstandardised and disappears when the covariates
                                are held fixed. A point-estimate statement; see the limitation below.
  "undecided at this n"         everything else, and in particular any reading with
                                targets_matched < 14 or dropped > 3.

THE TWO BRANCHES ARE NOT EQUALLY HARD, ON PURPOSE, AND THE VERDICT PRINTS THE REASON. Survives needs
both arms, two stratifications, a p-value and a matched-n floor; artefact needs one arm and two point
estimates. The strong form of the negative needs about 100 loci and the panel has 17, so a weak
negative registered honestly beats a strong one that cannot be delivered. The asymmetry is in
difficulty, not in direction, and it is printed beside the verdict rather than defended afterwards.

THE EXPECTED OUTCOME IS "UNDECIDED", and that is registered rather than merely allowed. Seventeen
positives is what the panel has; the standardised comparison's standard error is set by them. At a
panel rate near 0.53 the standard error of the difference is about 0.13, so the smallest difference
separable from zero at one-sided p 0.05 is about 0.21. The claim's unstandardised excess is about
0.35, so "survives" is reachable only if standardisation leaves the effect nearly untouched; any
attenuation into 0.10--0.21 lands in undecided and no arrangement of the control arms can rescue it.
To bound the difference BELOW 0.10 with 95% confidence -- the strong form of "artefact" -- needs a
standard error near 0.05, which is about 100 published loci. The panel has 17 and this lane is not
choosing loci. So the strong form of the negative is NOT attainable here, and the "artefact" branch
above is deliberately written as a collapse of a point estimate, which is weaker, and is to be
reported as weaker.
  n needed: "survives" needs the effect to stay near 0.35 at n = 17; a durable "survives" at an
  attenuated effect of 0.20 needs about 35 loci; "artefact" in its strong form needs about 100.

THE COMBINATION RULE, fixed now because it will be under most pressure when the number arrives.
When this lane finishes there are three figures about this claim: the panel's own 9/17, a second
lane's figure on loci this lane has not scored, and this lane's arms. THEY ARE NEVER POOLED, and no
"two of the three agree" reading is to be taken. Three reasons, all structural: this lane's arm A IS
the panel's 17 positives, so pooling it with the panel's 9/17 counts the same seventeen observations
twice; the second lane's loci come from a different sampling frame chosen by a different rule; and
arms B and C are control arms, which are not estimates of the claim's rate at published loci at all.
They are reported side by side with their frames named, and a disagreement between them stays a
disagreement.

COVERAGE, AS NAMED COUNTED CATEGORIES. A zero that means measured-and-absent and a zero that means
never-looked are separate rows in every arm's table:
  never_looked_syntax     no syntax reading came back for the window (a failed range read)
  no_variable_position    looked at, and nobody in the trio differs from the reference here
  values_none_constrained looked at, values present, none of them above the phyloP threshold
  values_in_syntax        the hit
and for the deletion layer, `deletion_scored` against `no_element_scored`, the latter being
never-looked. Coverage enters the comparison as the stratum `deletion_scored`, because a finding
from this project today is that standardising on GC, distance and constraint does not remove a
coverage artefact: none of those covariates says whether anybody ever spent a measurement here.

MODEL REQUESTS. Zero are registered and zero are to be spent. Every reading this lane needs is
already computed or is a range read over a public track.

WHAT THE AUTHOR KNEW, stated because a leak through selection leaves no trace.
  1. The loci were not chosen by this lane. The control windows were, by the seeded rule above.
  2. Before writing this, the author had read the panel arm's per-locus storage hits against its
     per-locus element lengths out of `data/results/loci_benchmark.json` -- outcome data on arm A --
     and had run a structural dry run that produced arm A against arm B under four of the five
     stratifications. Arm B IS the existing comparison, and neither reading can have steered arm C's
     draw, which is seeded and mechanical. Named here rather than left implicit.
  3. The author had also seen, in a paragraph of `docs/LOCI-BENCHMARK.md` section 9 that quotes a
     later section, a second lane's figures for this claim and that lane's own remark about element
     length. That is why this lane is reported as unblinded, at the top of this text rather than
     here. Nothing in the rule above was chosen to agree or disagree with it.

THE REVISION, 2026-09-17, after an outside read by the coordinating session and BEFORE any figure in
this module was scored. The order of operations is the whole point, so it is written down rather
than folded in.
  a. WHICH STRATIFICATION DECIDES, PER ARM. As first registered, P1_length decided for both control
     arms. Arm C is unmatched on distance to a coding TSS by construction, so a low rate there could
     be distance rather than length and the artefact branch could have fired for the wrong reason.
     Arm C now decides on S1_covariates and S2_covariates_coverage; arm B still decides on P1_length
     and P2_length_coverage. All five stratifications are still computed and printed for both arms,
     so nothing is hidden by the change -- only which one carries the verdict.
  b. WHAT ARM C IS FOR. Reframed from "does length explain the excess" to the complement question,
     because the construction fact above answers the first one without an experiment.
  c. THE BLINDING. Restated at the top as failed, on timing rather than on scope, rather than
     described as a limitation at the foot.
  Nothing else moved: not the thresholds, not the tolerance, not the denominators, not the
  combination rule, not the expected outcome.
"""


# ------------------------------------------------------------------ the new arm, drawn mechanically
def keep_out(ch: Chromosome, results_dir: Path = RESULTS_DIR) -> list[tuple[int, int]]:
    """Every interval a control window must avoid: the panel with its flank, VISTA, GWAS hits.

    The same three exclusions `loci.candidate_windows` applies, rebuilt here because this arm has to
    drop that function's GC, distance and constraint matching and the two are one expression there.
    """
    out = [
        (e.window[0] - NEGATIVE_KEEP_OUT, e.window[1] + NEGATIVE_KEEP_OUT)
        for e in PANEL
        if e.chrom == ch.chrom
    ]
    vista = load_result(f"vista_{ch.chrom}", results_dir) or {}
    out += [(x["start"] - VISTA_FLANK, x["end"] + VISTA_FLANK) for x in vista.get("rows", [])]
    gwas = Path("data/knowledge/gwas/hits.tsv")
    if gwas.exists():
        with gwas.open() as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0] == ch.chrom:
                    out.append((int(f[1]) - GWAS_FLANK, int(f[1]) + GWAS_FLANK))
    out.sort()
    return out


def blocked(intervals: list[tuple[int, int]], starts: list[int], s: int, e: int) -> bool:
    i = bisect.bisect_right(starts, e)
    return any(a < e and s < b for a, b in intervals[max(0, i - 400) : i])


def length_only_windows(
    ch: Chromosome,
    expect: Expect,
    keep: int = CONTROLS_PER_LOCUS,
    seed: int = DRAW_SEED,
    results_dir: Path = RESULTS_DIR,
) -> list[Window]:
    """`keep` windows of exactly the element's length, matched on length and on nothing else."""
    rng = random.Random(f"{seed}:{expect.locus}")
    length = expect.element[1] - expect.element[0]
    ivs = keep_out(ch, results_dir)
    starts = [k[0] for k in ivs]
    out: list[Window] = []
    for _ in range(DRAW_ATTEMPTS):
        if len(out) >= keep:
            break
        s = rng.randrange(1_000_000, max(1_000_001, ch.length - length - 1_000_000))
        e = s + length
        if blocked(ivs, starts, s, e):
            continue
        gc = ch.gc(s, e)
        if gc is None:  # an all-N window is not a control, it is a gap in the assembly
            continue
        near = ch.nearest_coding((s + e) // 2)
        if near is None:
            continue
        if any(w.start < e and s < w.end for w in out):
            continue
        out.append(Window(ch.chrom, s, e, gc=gc, distance=near[1], of_locus=expect.locus))
    return out


# ------------------------------------------------------------------------------------- the rows
def row_of(
    locus: str,
    chrom: str,
    start: int,
    end: int,
    gc: float | None,
    distance: int | None,
    constrained: float | None,
    syntax: dict[str, Any] | None,
    deletion: dict[str, Any] | None,
) -> dict[str, Any]:
    """One window as `compare.standardised` wants it, with the two zeros kept apart.

    `values` is None when nothing looked; 0 when the window was read and holds no variable position.
    `storage` is the claim. `has_values` and `deletion_scored` are the two coverage strata.
    """
    values = None if syntax is None else syntax.get("values")
    in_syntax = None if syntax is None else syntax.get("values_in_syntax")
    return {
        "locus": locus,
        "chrom": chrom,
        "start": start,
        "end": end,
        "length": end - start,
        "gc": gc,
        "distance": distance,
        "constrained": constrained,
        "values": values,
        "values_in_syntax": in_syntax,
        "storage": bool(in_syntax),
        "has_values": 1 if values else 0,
        "deletion_scored": 1 if (deletion or {}).get("elements_scored") else 0,
        "coverage": coverage_of(values, in_syntax),
    }


def coverage_of(values: int | None, in_syntax: int | None) -> str:
    """The registered four categories: never-looked is not measured-and-absent."""
    if values is None:
        return "never_looked_syntax"
    if not values:
        return "no_variable_position"
    return "values_in_syntax" if in_syntax else "values_none_constrained"


COVERAGE_CATEGORIES = (
    "never_looked_syntax",
    "no_variable_position",
    "values_none_constrained",
    "values_in_syntax",
)


def quantiles(xs: list[float]) -> dict[str, Any]:
    s = sorted(xs)
    if not s:
        return {"n": 0}
    return {
        "n": len(s),
        "min": s[0],
        "q1": s[len(s) // 4],
        "median": median(s),
        "q3": s[(3 * len(s)) // 4],
        "max": s[-1],
    }


def summarise_arm(name: str, rows: list[dict[str, Any]], what: str) -> dict[str, Any]:
    """Length, GC, median distance to a coding TSS, and coverage, for one arm. Every arm prints all
    four, so a reader never has to go and find one of them."""
    hits = sum(1 for r in rows if r["storage"])
    with_values = [r for r in rows if r["has_values"]]
    return {
        "arm": name,
        "what": what,
        "n": len(rows),
        "length": quantiles([r["length"] for r in rows]),
        "gc": quantiles([r["gc"] for r in rows if r["gc"] is not None]),
        "distance_to_coding_tss": quantiles([r["distance"] for r in rows if r["distance"] is not None]),
        "constrained_fraction": quantiles([r["constrained"] for r in rows if r["constrained"] is not None]),
        "coverage": {c: sum(1 for r in rows if r["coverage"] == c) for c in COVERAGE_CATEGORIES},
        "deletion_scored": sum(r["deletion_scored"] for r in rows),
        "no_element_scored": sum(1 for r in rows if not r["deletion_scored"]),
        "storage": {
            "k": hits,
            "n": len(rows),
            "rate": round(hits / len(rows), 4) if rows else None,
        },
        "storage_given_a_value_exists": {
            "k": sum(1 for r in with_values if r["storage"]),
            "n": len(with_values),
            "rate": round(sum(1 for r in with_values if r["storage"]) / len(with_values), 4)
            if with_values
            else None,
        },
    }


def length_match_check(panel: list[dict[str, Any]], arm: list[dict[str, Any]]) -> dict[str, Any]:
    """Whether an arm really is length-matched to the panel, per locus, from the coordinates."""
    by_locus = {r["locus"]: r["length"] for r in panel}
    worst = 0
    bad = []
    for r in arm:
        want = by_locus.get(r["locus"])
        if want is None:
            continue
        d = abs(r["length"] - want)
        worst = max(worst, d)
        if d > LENGTH_TOLERANCE_BP:
            bad.append({"locus": r["locus"], "want": want, "got": r["length"]})
    return {
        "tolerance_bp": LENGTH_TOLERANCE_BP,
        "worst_difference_bp": worst,
        "windows_outside_tolerance": len(bad),
        "examples": bad[:5],
        "matched": not bad,
    }


def compare_arms(panel: list[dict[str, Any]], controls: list[dict[str, Any]]) -> dict[str, Any]:
    """Every registered stratification of one comparison, covariates only and coverage held fixed."""
    return {name: standardised(panel, controls, strata, hit=HIT) for name, strata in STRATA.items()}


# ------------------------------------------------------------------------------- the verdict
def _ok(out: dict[str, Any]) -> bool:
    return (
        out["targets_matched"] >= MIN_TARGETS_MATCHED and out["dropped_for_want_of_a_control"] <= MAX_DROPPED
    )


def verdict(comparisons: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """The registered decision rule, applied mechanically to the two control arms' comparisons."""
    reasons: list[str] = []
    survives = True
    artefact = []
    for arm in ("B_existing", "C_length_only"):
        prim, cov = PRIMARY_BY_ARM[arm], PRIMARY_WITH_COVERAGE_BY_ARM[arm]
        a, b = comparisons[arm][prim], comparisons[arm][cov]
        d1, p_1 = a["matched"]["difference"], a["matched"]["p_one_sided"]
        d2 = b["matched"]["difference"]
        good = (
            _ok(a)
            and _ok(b)
            and d1 is not None
            and d2 is not None
            and d1 >= SURVIVES_DIFFERENCE
            and p_1 is not None
            and p_1 <= SURVIVES_P
            and d2 >= SURVIVES_WITH_COVERAGE
        )
        survives = survives and good
        reasons.append(
            f"{arm}: decides on {prim} -> difference {d1} at p {p_1}"
            f" (matched {a['targets_matched']}, dropped {a['dropped_for_want_of_a_control']}),"
            f" with coverage held fixed ({cov}) difference {d2}"
        )
        raw = a["raw"]["difference"]
        if raw is not None and d1 is not None and raw >= ARTEFACT_RAW_FLOOR and d1 <= ARTEFACT_DIFFERENCE:
            artefact.append(arm)
    if artefact:
        call = "the claim is an artefact of length"
    elif survives:
        call = "the claim survives"
    else:
        call = "undecided at this n"
    return {
        "verdict": call,
        "registered_words": [
            "the claim survives",
            "the claim is an artefact of length",
            "undecided at this n",
        ],
        "artefact_arms": artefact,
        "reasons": reasons,
        "decides_on": dict(PRIMARY_BY_ARM),
        "decides_on_with_coverage": dict(PRIMARY_WITH_COVERAGE_BY_ARM),
        "rule": (
            f"survives: the deciding difference >= {SURVIVES_DIFFERENCE} at one-sided p <="
            f" {SURVIVES_P} at BOTH control arms, and >= {SURVIVES_WITH_COVERAGE} at both with"
            f" coverage held fixed, with targets_matched >= {MIN_TARGETS_MATCHED} and dropped <="
            f" {MAX_DROPPED}; artefact: raw >= {ARTEFACT_RAW_FLOOR} while the deciding difference <="
            f" {ARTEFACT_DIFFERENCE} at either arm; otherwise undecided"
        ),
        "the_branches_are_not_equally_hard": (
            "on purpose: survives needs both arms, two stratifications, a p-value and a matched-n"
            " floor, while artefact needs one arm and two point estimates. The strong form of the"
            " negative needs about 100 loci and the panel has 17, so a weak negative registered"
            " honestly beats a strong one that cannot be delivered. The asymmetry is in difficulty,"
            " not in direction"
        ),
        "expected_before_the_run": "undecided at this n",
    }


# ------------------------------------------------------------------------------------- the build
def panel_rows(source: dict[str, Any]) -> list[dict[str, Any]]:
    """Arm A out of the committed benchmark: the same readings the 9/17 was counted from."""
    out = []
    for r in source["loci"]:
        e = r["expected"]
        out.append(
            row_of(
                e["locus"],
                e["chrom"],
                e["element"][0],
                e["element"][1],
                r.get("gc"),
                r.get("distance_to_coding_tss"),
                r.get("constrained_fraction"),
                r["readings"].get("syntax_values"),
                r["readings"].get("deletion"),
            )
        )
    return out


def existing_control_rows(source: dict[str, Any]) -> list[dict[str, Any]]:
    """Arm B out of the committed benchmark. The committed row keeps `values_in_syntax` but not
    `values`, so the two zeros cannot be told apart from the file: those rows are recomputed by
    `build` and this function leaves `values` unknown until they are."""
    out = []
    for w in source["negatives"]:
        c = w["claims"]
        out.append(
            row_of(
                w["of_locus"],
                w["chrom"],
                w["start"],
                w["end"],
                w.get("gc"),
                w.get("distance_to_coding_tss"),
                w.get("constrained_fraction"),
                {"values": None, "values_in_syntax": w.get("values_in_syntax")},
                {"elements_scored": 1 if c.get("deletion_scored") else 0},
            )
        )
    return out


def _fill_syntax(
    ch: Chromosome, rows: list[dict[str, Any]], progress=None
) -> None:  # pragma: no cover - network
    """Read syntax against values over one chromosome's windows, in one bigWig session."""
    if not rows:
        return
    got = read_syntax_values_many(ch, [(r["start"], r["end"]) for r in rows], None)
    for r, sv in zip(rows, got, strict=True):
        r["values"] = sv.get("values")
        r["values_in_syntax"] = sv.get("values_in_syntax")
        r["storage"] = bool(sv.get("values_in_syntax"))
        r["has_values"] = 1 if sv.get("values") else 0
        r["coverage"] = coverage_of(r["values"], r["values_in_syntax"])
    if progress:
        progress(f"{ch.chrom}: syntax against values over {len(rows)} windows")


def build(results_dir: Path = RESULTS_DIR, network: bool = True, progress=None) -> dict[str, Any]:
    """The three arms, their covariates, their coverage, the comparison both ways, and the verdict."""
    t0 = time.time()
    say = progress or (lambda _m: None)
    source = load_result("loci_benchmark", results_dir)
    if source is None:
        raise FileNotFoundError("data/results/loci_benchmark.json is not here; run scripts/loci_benchmark.py")
    arm_a = panel_rows(source)
    arm_b = existing_control_rows(source)

    arm_c: list[dict[str, Any]] = []
    per_locus: dict[str, int] = {}
    chroms: dict[str, Chromosome] = {}
    try:
        for e in PANEL:
            ch = chroms.get(e.chrom) or chroms.setdefault(e.chrom, Chromosome(e.chrom, results_dir))
            wins = length_only_windows(ch, e, results_dir=results_dir)
            per_locus[e.locus] = len(wins)
            for w in wins:
                arm_c.append(
                    row_of(
                        e.locus,
                        w.chrom,
                        w.start,
                        w.end,
                        w.gc,
                        w.distance,
                        None,
                        None,
                        read_deletion(ch, w.start, w.end, results_dir),
                    )
                )
            say(f"{e.locus}: {len(wins)} length-only controls drawn")

        if network:
            for chrom, ch in chroms.items():
                mine = [r for r in arm_c if r["chrom"] == chrom]
                if mine:
                    read = read_constraint(chrom, [(r["start"], r["end"]) for r in mine])
                    for r, c in zip(mine, read, strict=True):
                        r["constrained"] = (c.get("mammals") or {}).get("fraction_above")
                _fill_syntax(ch, mine, say)
                _fill_syntax(ch, [r for r in arm_b if r["chrom"] == chrom], say)
    finally:
        for ch in chroms.values():
            ch.close()

    arms = {
        "A_panel": summarise_arm("A_panel", arm_a, "the 17 published elements"),
        "B_existing": summarise_arm(
            "B_existing", arm_b, "5 per locus, matched on length, GC, distance and constraint"
        ),
        "C_length_only": summarise_arm(
            "C_length_only", arm_c, "5 per locus, matched on length and on nothing else"
        ),
    }
    comparisons = {
        "B_existing": compare_arms(arm_a, arm_b),
        "C_length_only": compare_arms(arm_a, arm_c),
    }
    return {
        "preregistration": PREREGISTRATION,
        "model_requests": 0,
        "claim": "a value sits on constrained sequence",
        "published_reading": "9/17 at the panel against 15/85 at the matched controls",
        "arms": arms,
        "length_matching": {
            "B_existing": length_match_check(arm_a, arm_b),
            "C_length_only": length_match_check(arm_a, arm_c),
        },
        "controls_drawn_per_locus": per_locus,
        "comparisons": comparisons,
        "verdict": verdict(comparisons),
        "combination_rule": (
            "never pooled with the panel's own 9/17 or with any other lane's figure: arm A IS the"
            " panel's seventeen positives, the other lane's loci come from a different sampling"
            " frame, and arms B and C are control arms rather than estimates of the claim's rate"
        ),
        "rows": {"A_panel": arm_a, "B_existing": arm_b, "C_length_only": arm_c},
        "network": network,
        "cost": {"seconds": round(time.time() - t0, 1), "model_requests": 0},
    }


def run_and_save(results_dir: Path = RESULTS_DIR, network: bool = True, progress=None) -> dict[str, Any]:
    out = build(results_dir=results_dir, network=network, progress=progress)
    save_result("loci_length_controls", out, results_dir)
    return out
