# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""guardcensus (2026-10-12): which guards have no planted violation that turns them red?

THE QUESTION. A guard nobody has shown can fail is a guard we do not know works. This census asks, of
the suite's guards, which ones have a planted violation that turns them RED -- and which have none.

THE RULE IS REGISTERED HERE BEFORE ANYTHING IS COUNTED, and that is the whole value of it: a census whose
classification is fixed after the list is seen can be tuned to the list, and this one cannot.

WHY A FOURTH CONDITION EXISTS AT ALL, and it is this lane's inheritance rather than its invention.
globcensus stage 2 classified six committed results as realising a silent shrink, from their manifests.
Its test had three conditions and they were all satisfied -- short member count, member on disk today, no
exclusion naming it -- and the finding was still WRONG: lane-shrinkimpact recomputed on the complete
stores and the shrink was realised in none of the six. One was a LABELLING ARTEFACT (the member was read
and filed under another group's label) and the rest were EMPTY CHROMOSOMES (the member exists and
contributes nothing eligible). All three cases look identical from the manifest. The supervisor's
correction added the condition that settles them: the missing member must CARRY DATA the result's own
scope would have used -- which can only be known by opening it.

SO THE CENSUS OVER DESCRIPTIONS COULD ESTABLISH A MECHANISM AND NOT THAT IT FIRED. That is the lesson
this lane is built on, and its direct analogue is the whole of CANNOT_FIRE_TEST below: reading a guard's
SOURCE can establish that a violation looks unconstructible, and cannot establish that the guard cannot
fire. Only planting the violation and watching the guard stay green establishes that.
"""

from __future__ import annotations

import ast
import random

#: The four conditions. Registered before counting; none is sufficient alone.
CANNOT_FIRE_TEST = (
    "a guard is counted CANNOT-FIRE only when all four hold: (1) it is a guard -- an assertion, raise"
    " or refusal whose job is to reject something; (2) no test in the suite plants a violation of it;"
    " (3) reading its source shows no input that would reject -- the condition is unreachable, vacuous,"
    " or compares values that can never differ; and (4) A VIOLATION WAS ACTUALLY PLANTED AND THE GUARD"
    " STAYED GREEN. Conditions (1) to (3) are description; (4) is measurement, and without it the"
    " census can only say a violation was not FOUND, never that none EXISTS"
)

#: The condition this lane exists because of, named so the provenance of the rule is on the record.
WHY_CONDITION_4 = (
    "globcensus stage 2 met its own three conditions on six results and was still wrong, because a"
    " labelling artefact and an empty chromosome are indistinguishable from a genuine shrink in the"
    " manifest. The supervisor's fourth condition there -- the missing member must carry data the"
    " result's own scope would have used -- is the same move: stop reading the description and open"
    " the thing. Here the thing to open is the guard, by planting a violation and running it"
)

#: What a failure to plant means, written down so it cannot be read as a result.
IF_NO_VIOLATION_CAN_BE_PLANTED = (
    "NOT_ESTABLISHED, and counted apart. A violation this lane could not construct is a gap in this"
    " lane, not a verdict about the guard. Reporting it as cannot-fire would repeat stage 2's error"
    " exactly: concluding from the absence of evidence that was never sought in the right place"
)

#: The generating mechanism, which the supervisor asked to be named as a mechanism rather than as more
#: cases. Eleven guards have gone inert this way since Saturday.
GENERATING_MECHANISM = (
    "COMMITTING A RESULT DISARMS THE SKIP GUARD KEYED ON ITS ABSENCE. A test written while a result does"
    " not yet exist guards itself with a skip on that absence; the day the result is committed the skip"
    " can never fire again, and nothing notices, because a skip that stops firing looks exactly like a"
    " test that passes. Eleven guards went inert this way since 2026-10-10. The prevention is in"
    " tests/test_inert_skip_guards.py's docstring at 4f2652c: write tracked_paths.must_be_present in"
    " the FIRST version, never a skip you will later have to convert -- because converting removes a"
    " line and therefore needs a force the harness denies"
)

#: The cases measured in the last three days, each with what was wrong and how it surfaced.
CASES = (
    {
        "case": "tests/test_closure_cell2.py:237",
        "kind": "inert skip",
        "what": "skipif keyed on a path git TRACKS, so the skip could never fire",
        "surfaced_by": "tests/test_inert_skip_guards.py, the standing invariant",
    },
    {
        "case": "score's boundary guard",
        "kind": "could never PASS",
        "what": "compared int dict keys against their JSON string copies, so no input satisfied it",
        "surfaced_by": "a lane reading the comparison rather than its verdict",
        "note": "the mirror image of the others: a guard that always fires is as uninformative as one"
        " that never can",
    },
    {
        "case": "lane-rebuild0's leaf comparison script (mine, UNDIAGNOSED)",
        "kind": "instrument that could not detect",
        "what": "returned '0 differences' where a byte-compare showed 151,927 against 90,107 bytes and"
        " a declared sha256 differing at char 49595",
        "surfaced_by": "reading the byte counts behind the zero; the defect itself was never diagnosed",
        "note": "included because an undiagnosed defect in a MEASURING instrument is exactly what this"
        " census should be able to find, and because it is mine",
    },
    {
        "case": "scripts/manifest_rebuild.py:309",
        "kind": "vacuous field",
        "what": "'reconciles': total == compared + (total - compared) -- the right-hand side IS total,"
        " so it is True for every input, including compared greater than total and negative counts",
        "surfaced_by": "the coordinator, on random values; verified independently",
        "note": "it sat inside the tool used to check the only positive result in its series, and had"
        " been cited as evidence",
    },
    {
        "case": "the eleven inert skips of 2026-10-10 onward",
        "kind": "the generating mechanism, not eleven separate cases",
        "what": "see GENERATING_MECHANISM",
        "surfaced_by": "tests/test_inert_skip_guards.py and its conversions at e07e5ce and 65bdee3",
    },
)

#: Scope, fixed now so it cannot widen to fit what is found.
SCOPE = (
    "the suite under tests/, plus guards in genomeos/ and scripts/ that a test is responsible for"
    " exercising. A LIST AND A COUNT, no repairs: nothing is fixed by this lane, and a guard found"
    " unexercised is reported rather than given a plant, because writing the plant is the owning lane's"
    " work and doing it here would hide how many were unexercised"
)

#: What the census will not claim, whatever the count turns out to be.
DOES_NOT_CLAIM = (
    "not that an unexercised guard is wrong -- most will be right and merely unproven; not that the"
    " count is the number of defects, since condition (4) is expensive and NOT_ESTABLISHED will be the"
    " largest class; not that a guard with a plant is therefore correct, because a plant can be"
    " satisfied by the wrong cause, which is the lesson of test_PLANTED_an_EFFECTIVE_line of"
    " 2026-10-03; and nothing about guards outside SCOPE"
)


# ------------------------------------------------------------- added 2026-10-12, after the rule landed
# Three further cases, appended rather than folded into CASES so the order they arrived in is visible.
# The first is the most interesting thing in this census, because nothing is broken in it.

#: The ordering case. Named as a mechanism, as the coordinator asked, because the instance is incidental.
ORDERING_MECHANISM = (
    "A PREVENTION RULE THAT MAKES A TEST RAISE EARLIER CAN MOVE THE DEFECT FROM 'SILENTLY INERT' TO"
    " 'RED FOR TWO MINUTES OF HISTORY', AND A PUSH THAT LANDS IN THOSE TWO MINUTES PAYS FOR IT."
    " GENERATING_MECHANISM's cure is tracked_paths.must_be_present from the first version. That is"
    " right. But must_be_present correctly RAISES NotTrackedError while the path is untracked, so a"
    " test committed BEFORE its record turns every sha in between red. lane-scalehet:"
    " tests/test_scalehet.py:570/581/589 read data/results/scalehet_result.json at 036eec6; the record"
    " landed at 40b0a8d two minutes later; a push captured a sha inside the gap and got 3 failed,"
    " 7,258 passed, while the same file passes 57 at HEAD. THE RULE: a test that reads a result record"
    " is committed IN THE SAME COMMIT as that record, or after it, never before"
)

#: Why it belongs in a census of guards that cannot fire, which is not obvious.
WHY_ORDERING_IS_THIS_FAMILY = (
    "because nothing in it is broken. The guard is correct, the prevention rule is correct, and the"
    " failure is in the ORDERING BETWEEN THEM: obeying either one alone produces a red that neither"
    " mentions. The other cases here are instruments that cannot detect; this one is two instruments"
    " that are each right and jointly wrong, which is a defect no single-guard audit can see and the"
    " reason this census counts mechanisms and not only sites"
)

#: Two further cases from lane-cired2, verified here at the source rather than carried.
CASES_ADDED = (
    {
        "case": "tests/extras_lock.py:401, read_module",
        "kind": "guard blind to its own subject",
        "what": "an ImportFrom is recorded by its MODULE only -- `names = [node.module] if"
        " node.level == 0 and node.module else []` -- so every `from package import submodule` edge is"
        " INVISIBLE to test_extra_only_imports.py, the guard written for exactly that defect",
        "surfaced_by": "lane-cired2; verified here by reading extras_lock.py:401 directly",
        "note": "the sharpest kind in this census: not a guard that cannot fire, but one that cannot"
        " SEE the class of thing it exists to catch, so its green is silence rather than evidence",
    },
    {
        "case": "requires_extra's AND semantics",
        "kind": "correct rule, inexpressible case",
        "what": "each named extra must be installed, so 'scipy by whichever extra supplies it' cannot"
        " be written. Costing 36 real tests on this machine, where scipy is installed and"
        " process_bigraph is not",
        "surfaced_by": "lane-cired2",
        "note": "the AND is DELIBERATE and its reason is good -- conftest.py records that seven tests"
        " passed on pandas because some other extra dragged it in. So this is a cost of a correct"
        " rule, not a defect in it, and it is counted as such rather than as a fault",
    },
)

#: The inherited condition, carried here verbatim because the supervisor ruled it belongs in this lane.
FOURTH_CONDITION_INHERITED = (
    "a missing member must CARRY DATA the result's own scope would have used (supervisor, 2026-10-12,"
    " after globcensus stage 2 was refuted). Registered here before counting, and generalised into"
    " CANNOT_FIRE_TEST's condition (4): description cannot settle what only opening the thing can"
)


# ------------------------------------------------ condition (3), MEASURED: the vacuity leg, 2026-10-12
# CANNOT_FIRE_TEST's condition (3) has a clause that can be measured instead of read: "compares
# values that can never differ". For a comparison built only of names, constants and arithmetic,
# values can be SUBSTITUTED and the comparison EVALUATED, thousands of times. That is condition
# (4)'s discipline -- plant the violation, do not reason about it -- applied to arithmetic.

#: What the instrument is, and what it is not.
VACUITY_INSTRUMENT = (
    "every ==/!= comparison in genomeos/, scripts/ and tests/ whose leaves are only names,"
    " constants and arithmetic is compiled and evaluated under random substitutions drawn from a"
    " pool containing THE EXPRESSION'S OWN CONSTANTS. A comparison whose verdict never changes over"
    " at least 150 usable evaluations is one whose two sides cannot differ. It is not a proof: it is"
    " a search that plants many violations and reports whether any was found"
)

#: The four defects found in this instrument, each of which changed the count. Kept because the
#: sequence IS the finding: a census is only as good as the instrument, and mine was wrong 4 times.
VACUITY_INSTRUMENT_WAS_WRONG_FOUR_TIMES = (
    {
        "version": 1,
        "count": 1600,
        "defect": "the pool held random INTEGERS ONLY, so `action == 'activates'` could not come"
        " out True and was reported always-False. The probe was measuring its own pool",
    },
    {
        "version": 2,
        "count": 20,
        "defect": "collecting ast.Constant alone missed -32768 (a UnaryOp) and 512 * 1024 * 1024"
        " and 490.0+50.0+7.0+3.0 (BinOps), so reachable comparisons were reported unreachable."
        " Whole sides are folded now",
    },
    {
        "version": 3,
        "count": 7,
        "defect": "a NaN in the pool DELETED THE ONE TRUE POSITIVE, because nan == nan is False."
        " The discriminator is not the pool but the SHAPE: `v == v` has syntactically identical"
        " sides and is the deliberate NaN idiom; a tautology has different sides that are"
        " algebraically equal. Self-comparisons are set aside by shape and counted separately."
        " Also: trials that RAISED were counted as agreeing, so `ne - ns == e - s` was called"
        " invariant on about 99 usable trials out of 400 -- a floor on usable trials fixed it",
    },
    {
        "version": 4,
        "count": 1,
        "defect": "the pool builder sorted constants into ints and strs, and BYTES fell through"
        " the gap, so six bytes comparisons -- k == b's', magic == b'\\x1f\\x8b', data == b'hello'"
        " -- were called tautologies because the only value satisfying them had been removed."
        " Every constant stays in the pool now, whatever its type",
    },
)

#: The reading. A NEGATIVE result, and the honest one.
VACUITY_RESULT = (
    "ONE tautology in genomeos/, scripts/ and tests/ together: scripts/manifest_rebuild.py:309,"
    " `total == compared + (total - compared)`, always True over 2,555 usable plants -- and it was"
    " ALREADY IN CASES before this leg ran, so the leg found nothing new. 16 self-comparisons were"
    " set aside as the NaN idiom, which is correct code. 0 sites went unmeasured. So the 'compares"
    " values that can never differ' clause of condition (3) is very nearly EMPTY in this codebase,"
    " and the clause's value is as a rule for reading new code rather than as a yield of old code"
)

#: Why the sequence of counts matters more than the count.
VACUITY_LESSON = (
    "four versions of one instrument gave 1,600, 20, 7 and 1, and only the last is right. Every"
    " earlier number was the POOL measuring itself: a value the code requires had been left out, so"
    " the comparison could not be satisfied and was reported as unsatisfiable. This is the census's"
    " own subject turned on the census -- an instrument that cannot detect reports absence, and"
    " absence reads like a finding. The validation that caught each one was a KNOWN POSITIVE held"
    " aside: manifest_rebuild.py:309 had to appear in the output, and version 3 was rejected"
    " precisely because it stopped appearing"
)

#: The 16 set aside, by path, so the exclusion can be audited rather than trusted.
NAN_IDIOM_SET_ASIDE = (
    "genomeos/attribution/runrules.py:627",
    "genomeos/cli.py:614",
    "genomeos/genome/epigenome.py:554",
    "genomeos/predict/alphagenome_adapter.py:407",
    "genomeos/predict/contact_maps.py:63",
    "genomeos/runtime/segmentation.py:44",
    "genomeos/web/server.py:1874",
    "scripts/context_contrast_feasibility.py:585",
    "scripts/cptac_phospho.py:348",
    "scripts/epigenome.py:368",
    "scripts/epigenome.py:441",
    "scripts/epigenome.py:925",
    "scripts/placement_audit.py:381",
    "scripts/placement_audit.py:561",
    "scripts/placement_audit.py:830",
    "scripts/residue_census.py:210",
)


def _fold(node: ast.expr) -> object | None:
    """The value of a sub-expression that names nothing, or None.

    A negative literal is a UnaryOp and ``512 * 1024 * 1024`` a BinOp, so collecting
    ``ast.Constant`` alone misses the value the name is actually compared against.
    """
    if any(isinstance(n, ast.Name) for n in ast.walk(node)):
        return None
    try:
        return eval(  # noqa: S307 -- a constant-only sub-expression of this repo's own source
            compile(ast.Expression(body=node), "<fold>", "eval"), {"__builtins__": {}}, {}
        )
    except Exception:
        return None


def pool_for(node: ast.Compare) -> list[object]:
    """The substitution pool: the expression's OWN constants, plus a few generic values.

    Every collected constant stays in, whatever its type -- see version 4 in
    :data:`VACUITY_INSTRUMENT_WAS_WRONG_FOUR_TIMES`.
    """
    vals = [n.value for n in ast.walk(node) if isinstance(n, ast.Constant)]
    for side in [node.left, *node.comparators]:
        folded = _fold(side)
        if folded is not None:
            vals.append(folded)
    return [*vals, 0, 1, 2, -1, 7, 50, 3, 11, "", "x"]


def classify_comparison(
    node: ast.Compare, *, trials: int = 4000, min_usable: int = 150
) -> tuple[str, bool | None, int]:
    """``("tautology" | "nan_idiom" | "varies" | "not_measured", verdict, usable_trials)``.

    A trial that RAISES is not a trial that agreed: a verdict needs ``min_usable`` successful
    evaluations or the site is reported as not measured rather than as invariant.
    """
    safe = (
        ast.Name,
        ast.Constant,
        ast.BinOp,
        ast.UnaryOp,
        ast.Add,
        ast.Sub,
        ast.Mult,
        ast.FloorDiv,
        ast.Mod,
        ast.USub,
        ast.UAdd,
        ast.Compare,
        ast.Eq,
        ast.NotEq,
        ast.Load,
    )
    if len(node.ops) != 1 or not isinstance(node.ops[0], (ast.Eq, ast.NotEq)):
        return "varies", None, 0
    if not all(isinstance(n, safe) for n in ast.walk(node)):
        return "varies", None, 0
    names = sorted({x.id for x in ast.walk(node) if isinstance(x, ast.Name)})
    if not names:
        return "varies", None, 0

    pool = pool_for(node)
    code = compile(ast.Expression(body=node), "<probe>", "eval")
    seen: set[bool] = set()
    usable = 0
    for _ in range(trials):
        env = {k: random.choice(pool) for k in names}
        try:
            got = bool(eval(code, {"__builtins__": {}}, env))  # noqa: S307 -- this repo's source
        except Exception:
            continue
        usable += 1
        seen.add(got)
        if len(seen) > 1:
            return "varies", None, usable
    if usable < min_usable:
        return "not_measured", None, usable
    verdict = seen.pop()
    # BY SHAPE, not by pool: identical sides are the deliberate NaN idiom, not a tautology.
    if ast.unparse(node.left) == ast.unparse(node.comparators[0]):
        return "nan_idiom", verdict, usable
    return "tautology", verdict, usable
