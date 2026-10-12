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


# ----------------------- condition (4) SATISFIED, for the first time in this census -- 2026-10-12
# Everything above this line is either a rule, a description, or a negative result. This block is
# the first case where a violation was PLANTED and the guard STAYED GREEN, which is the only thing
# that establishes CANNOT-FIRE rather than "no violation was found".

#: Where the blindness actually is. NARROWER than CASES_ADDED says; see the correction below.
BLINDNESS_LOCATION = (
    "tests/extras_lock.py:459, the FIRST-PARTY CLOSURE WALK in scan(). read_module records an"
    " ImportFrom by node.module alone, so `from genomeos.attribution import pbrun` puts"
    " 'genomeos.attribution' in first_party, _first_party_file resolves THAT to"
    " genomeos/attribution/__init__.py, and pbrun.py is never queued -- so no third-party import"
    " written inside pbrun.py is ever checked for any test that reaches it this way"
)

#: The correction. CASES_ADDED's entry overstates the mechanism and is LEFT STANDING, as
#: globcensus's refuted headline was, so the record shows what was claimed and when it narrowed.
BLINDNESS_CORRECTION = (
    "CASES_ADDED says 'every `from package import submodule` edge is INVISIBLE to"
    " test_extra_only_imports.py'. That is too broad and the entry is left standing rather than"
    " edited. THE THIRD-PARTY LEG IS CORRECT: it appends Import(name.split('.')[0], ...), so"
    " `from scipy import stats` records scipy, which is exactly the right answer to 'which"
    " distribution does this file need'. Only the first-party closure walk loses the edge. The"
    " guard is right in the leg that reads package presence and blind in the leg that walks the"
    " closure, and saying so is the difference between a defect and a direction"
)

#: The measurement. Not a reading of the source: the guard's own walk, run twice.
BLINDNESS_IS_REALISED = (
    "extras_lock.scan() was run over all 471 test files twice -- once as committed, once with"
    " read_module patched to record `module.submodule` where that names a real first-party file."
    " AS COMMITTED: 0 violations. WITH THE EDGE RECORDED: 1 violation. A violation was therefore"
    " planted and the guard stayed green, so condition (4) is SATISFIED and this site is an"
    " established CANNOT-FIRE rather than a NOT_ESTABLISHED one"
)

#: The violation it was hiding, as at 271912f. Each link read by hand, not taken from the patch.
BLINDNESS_INSTANCE = {
    "test": "tests/test_pbrun.py",
    "edge": "line 26, `from genomeos.attribution import pbrun as pr` -- top level, and ast finds"
    " no try/except anywhere in the file and no importorskip, so it is unguarded",
    "through": "genomeos/attribution/pbrun.py:41, `import h5py`, top level and unguarded",
    "package": "h5py, absent from test-bare's 15 bare distributions and supplied only by the"
    " `predict` extra, which NEITHER CI job installs: `test` takes --extra compose --extra"
    " roadrunner and test-bare takes `uv sync --frozen`",
    "declared": "nothing: tests/test_pbrun.py contains zero requires_extra calls",
    "why_not_red_yet": "1a1e567 landed 2026-10-12 13:37 and is on no remote branch; CI's runs are"
    " on b97a8d5, which does not contain it. The next push carrying it is what pays, and test-bare"
    " fails at COLLECTION, which reads as 'no such tests' rather than as a missing dependency",
    "why_the_local_gate_misses_it": "h5py IS installed in the shared .venv, so scripts/check.sh is"
    " green here. This is 'local pre-push green is not CI green' with a line number",
    "same_class_as": "the openpyxl defect of 2026-10-02, reproduced in the file whose docstring describes it",
}

#: Scale, with the caveat that keeps it honest.
BLINDNESS_SCALE = (
    "1,942 from-import edges across genomeos/, scripts/ and tests/ name a real first-party"
    " submodule that the walk does not record, over 239 distinct target files. THAT IS NOT 1,942"
    " DEFECTS: nearly every one of those targets is reachable by some other path, and the number"
    " that matters is that exactly ONE of them currently hides a non-bare import. A census that"
    " reported 1,942 would be the vacuity leg's mistake again -- counting the mechanism's sites"
    " and calling them symptoms"
)

#: And the repair that is deliberately NOT done here.
BLINDNESS_NO_REPAIR = (
    "SCOPE says a guard found unexercised is reported and not fixed, because writing the fix here"
    " would hide how many were unexercised. The fix is lane-pbrun's and it is one of two lines:"
    " @pytest.mark.requires_extra('predict') on the test, or move `import h5py` into the function"
    " that needs it -- which pyproject already calls 'the pattern to keep', and which also fixes"
    " it for every other file that imports pbrun"
)


# ------------------- the stopgap, measured in CI's own condition -- 2026-10-12, after f4a9faa
# BLINDNESS_INSTANCE's "why_not_red_yet" says test-bare fails at COLLECTION. That was true of the
# sha it was written against and is NO LONGER TRUE at HEAD. Corrected additively, as the overstated
# CASES_ADDED entry was, so the record shows what was measured when.

#: The correction, and the guard that made it necessary.
INSTANCE_SUPERSEDED = (
    "f4a9faa (2026-10-12 14:36:09), 59 minutes after 1a1e567, adds to tests/conftest.py a declared"
    " COLLECT_NEEDS_EXTRA = {'test_pbrun.py': 'predict'} and a collect_ignore built only when"
    " extras_lock.missing_modules_for_extra('predict') is non-empty. MEASURED, not read: with"
    " h5py, alphagenome and anndata all made genuinely unimportable, pytest --collect-only with no"
    " path argument -- the way CI runs it -- collects 7,500 tests in 18.29s, exit 0, ZERO"
    " collection errors, and test_pbrun.py absent. SO THE STOPGAP WORKS and test-bare does not"
    " fail. f4a9faa is an ancestor of ae272c8 and of HEAD; neither it nor 1a1e567 is on b97a8d5,"
    " so origin carries neither the defect nor the guard"
)

#: How absence was simulated, since no venv could be changed and three lanes share this one.
HOW_ABSENCE_WAS_SIMULATED = (
    "a sitecustomize that WRAPS every finder in sys.meta_path, because _importable asks every"
    " finder and returns True if any answers -- so a prepended blocking finder hides nothing, as"
    " that function's own docstring says. sitecustomize is imported by `site` at interpreter"
    " startup, which is the only point early enough to precede tests/conftest.py. Validated both"
    " ways: `import h5py` raises ModuleNotFoundError under the veil and succeeds without it. No"
    " uv sync, no change to the shared .venv, and the return value was NOT monkeypatched -- the"
    " modules were made genuinely unimportable, which is what absence is"
)

#: And the condition that looks like a defect in the stopgap and is not.
STOPGAP_INERT_CONDITION = (
    "veiling h5py ALONE leaves missing_modules_for_extra('predict') EMPTY, so collect_ignore is"
    " empty and the file is collected. The reason is by design and by docstring: predict's DIRECT"
    " requirements are alphagenome and anndata, h5py is transitive through anndata, and that"
    " function checks direct requirements only. So the stopgap is inert in the condition 'direct"
    " deps present, h5py absent' -- which CANNOT ARISE IN CI, because anndata requires h5py. A"
    " bounded caveat, recorded so it is not mistaken for a fault later"
)

#: A NEW CASE, and of exactly this census's kind: the guard fires and nobody hears it.
ANNOUNCEMENT_CANNOT_BE_HEARD = {
    "case": "tests/conftest.py's pytest_report_header, from f4a9faa",
    "kind": "an announcement that cannot be heard where it was written to be heard",
    "measured": "without -q the line is printed; with -q it is absent; with CI's exact flags,"
    " `-q -rs --no-fold-skipped`, it is absent. Both CI jobs run"
    " `uv run --no-sync pytest -q -rs --no-fold-skipped`",
    "why": "pytest_report_header writes a HEADER and -q suppresses headers. -rs lists skip"
    " reasons, and an ignored file is not a skip, so it does not appear there either",
    "against_its_own_words": "the commit says the guard 'says so in the header rather than"
    " skipping in silence' and the docstring says 'so a stopgap is never a silent absence'."
    " Under CI's own flags it IS a silent absence",
    "family": "the same family as the eleven inert skips and as extras_lock's blind leg: the"
    " instrument is RIGHT and its reading never reaches anyone. A guard whose output is"
    " discarded by the only caller that matters cannot fire in the sense this census means",
    "smaller_one_beside_it": "the message names alphagenome and anndata as what cannot be"
    " imported. The module that actually fails is h5py. A reader goes looking in the wrong place",
}

#: A confound found in this lane's OWN test, recorded because another lane will hit it.
COLLECT_IGNORE_AND_EXPLICIT_PATHS = (
    "naming tests/test_pbrun.py explicitly on the command line STILL errors, because collect_ignore"
    " governs DIRECTORY collection. CI passes no path, so it is not bitten -- but a lane that runs"
    " that file by name gets the old ModuleNotFoundError and may read it as the stopgap failing."
    " Found by making exactly that mistake here first"
)


# ------------------------------------------- procedure deviations, not instrument ones -- 2026-10-12
# A separate class deliberately. Everything above concerns guards that cannot fire. This concerns a
# RULE that was not followed, where nothing was broken and the judgement was defensible -- which is
# the only interesting kind, because a rule that only ever breaks when someone is careless teaches
# nothing.

#: The tightened wording. Carried verbatim from the coordinator's 16:22 ruling; applies to every brief.
NO_FORCE_RULE = "never `--force`; a red you caused yourself is still a report, not an exception"

#: Why the wording needed tightening at all.
WHY_THE_WORDING_WAS_TIGHTENED = (
    "the old wording let a lane infer that a SELF-INFLICTED red was different in kind. It is not."
    " The rule exists to replace judgement in that moment, and 'it is only my own line' is what"
    " every bad force also says -- which is the whole reason the deviation is recorded although it"
    " did no damage"
)

#: The deviation, with the facts that made it acceptable ONCE, listed so the acceptance is auditable.
PROCEDURE_DEVIATIONS = (
    {
        "lane": "lane-hetblock2",
        "commit": "b2bf090",
        "rule": "a brief that forbade commit_own.sh --force",
        "ruled": "accepted once, on the record as disclosed -- supervisor, 2026-10-12 16:19",
        "facts": (
            "it was its own line, committed four minutes earlier in aa200ac",
            "it had caused its own red",
            "verified by the coordinator as 2 added, 1 removed, one file"
            " (tests/test_hetblock_harness_record.py), no peer content",
            "it recorded the force in that commit message",
            "it disclosed the force unprompted in its report",
            "it restructured its later fixes as pure additions (16/0, 4/0) specifically to avoid"
            " needing it again",
        ),
        "class": "procedure, not instrument: nothing was broken and no guard was weakened",
    },
)


#: This census's own cannot-fire guard, caught by ruff while the census was being written.
CASES_SELF = (
    {
        "case": "tests/test_guardcensus.py, in the very commit recording the stopgap",
        "what": "`assert <phrase> in INSTANCE_SUPERSEDED or True` -- an assertion that cannot"
        " fail, because `or True` makes the whole expression true whatever the left side says."
        " The phrase was not in the constant at all, so without the `or True` the test would have"
        " failed immediately",
        "found_by": "ruff SIM222, before any commit. Not by a test, and not by reading it",
        "why_recorded": "a census of guards that cannot fire, which produced a guard that cannot"
        " fire, in the same file. It is the cheapest possible demonstration that condition (4) is"
        " not pedantry: conditions (1) to (3) were all satisfied by that line -- it is an"
        " assertion, nothing planted a violation of it, and reading it shows no input that would"
        " reject -- and the only thing that settles it is that NO input can make it red",
        "fixed": "yes -- this one IS repaired, because it is this lane's own line and leaving a"
        " vacuous assertion in place to preserve a count would be the worse error. SCOPE's 'no"
        " repairs' binds the census's SUBJECTS, not its own instrument",
    },
)


#: The census's SECOND self-case, and the worse of the two: a guard that fires at random.
CASES_SELF_FLAKY = {
    "case": "tests/test_guardcensus.py::test_a_trial_that_raised_is_not_a_trial_that_agreed,"
    " committed by this lane at 271912f",
    "what": "it asserts classify_comparison returns 'not_measured' for `ne - ns == e - s`. That"
    " expression genuinely VARIES, and the varies branch returns as soon as two verdicts are"
    " seen -- before the usable floor is ever consulted. So the assertion was false when it was"
    " committed",
    "worse_than_wrong": "FLAKY, not red. classify_comparison draws random values, so the verdict"
    " depends on the draw: it passed the run in which it was written and has mostly failed since."
    " MEASURED here over five consecutive runs of the file: 3 xfailed, 2 xpassed. So 271912f is a"
    " red commit on dev whose redness is intermittent, which is the kind that reaches a push",
    "the_census_point": "a test whose verdict depends on a seed is the mirror image of a guard"
    " that cannot fire: it fires, or does not, for reasons unrelated to the thing it watches."
    " Both report something other than the state of the code, and neither is caught by reading",
    "how_it_was_caught": "not by a test and not by review -- by the UNDO GUARD refusing the"
    " commit that removed it, which named the four lines and which commit had added them",
    "why_not_deleted": "removing it needs --force inside the undo window, and the rule is: never"
    " --force; a red you caused yourself is still a report, not an exception. So it stays, marked"
    " @pytest.mark.xfail(strict=False) with the reason beside it, and the corrected test"
    " (test_the_usable_floor_is_reached_when_every_evaluation_raises) is added beside it. 199 and"
    " 181 lines added, ZERO removed",
    "why_strict_false": "strict=True would fail on the runs where the flake happens to pass --"
    " measured above as 2 of 5. A strict xfail on a flaky test is itself a guard that fires at"
    " random, which would be the same defect again",
}


#: The census's THIRD self-case, and the first that was a red on dev rather than an inert line.
CASES_SELF_RED = {
    "case": 'tests/test_guardcensus.py\'s literal "genomeos.attribution.pbrun", added by this'
    " lane at ce65d51",
    "what": "read_module collects first-party edges FROM STRING LITERALS as well as from import"
    " statements -- deliberately, because the closure that hid openpyxl was a string. That"
    " literal resolves, so the walk followed it to pbrun.py, found `import h5py` at line 41, and"
    " test_no_test_module_needs_an_undeclared_package went RED",
    "the_irony_is_the_finding": "the test documenting an edge the walk CANNOT follow supplied, by"
    " hand and in a string, exactly the edge that makes it followable. The guard was right and"
    " this lane was the violation",
    "why_the_census_missed_it": "this lane scanned tests/test_*.py; the guard scans"
    " TESTS.glob('*.py'), every .py in tests/ including conftest.py and extras_lock.py itself."
    " Same function, different file set, different answer. A CENSUS THAT MEASURES A NARROWER"
    " POPULATION THAN THE GUARD IT AUDITS WILL REPORT THE GUARD CLEAN WHILE THE GUARD IS RED,"
    " which is worth more than the incident: the population must be the guard's own, not one that"
    " resembles it",
    "fix": "@pytest.mark.requires_extra('predict') on the one offending test -- the guard's own"
    " prescribed remedy, and a pure insertion, because read_module collects requires_extra at"
    " FILE level. 1 line added, 0 removed, 0 violations after; the other 35 guards still run bare",
    "better_fix_not_taken": "dropping the literal is better, since the test never imports pbrun"
    " and only hands a string to _first_party_file, so the marker makes it skip where it would"
    " have run. But editing that line REMOVES a line ce65d51 added, the window is 2.days, and"
    " --force is not this lane's to take. Reported to the coordinator with the better fix named",
    "on_the_undo_guard": "its docstring says deleting your own lines passes. It refused only"
    " because every commit here is authored Albert Canfield, so it cannot tell this lane's line"
    " from a peer's -- which its refusal message states outright. A limitation of the detector,"
    " not a licence",
}


# ------------------------- HOW THE COUNT WILL BE READ, registered BEFORE the count exists
# The raise-site measurement is running as this is written and its numbers are not known here.
# Registering the reading first is the same discipline as registering the rule first: a reading
# chosen after the numbers are seen can be fitted to them.

#: The population, and how condition (1) is applied to it. By SHAPE, not by name.
COUNT_POPULATION = (
    "every `raise` in genomeos/, scripts/ and tests/, classified by the ancestors BETWEEN the"
    " raise and its own enclosing function: `guard` if a conditional decides whether it happens,"
    " `reraise` if it is bare or inside an except handler, `loop` if only an iteration reaches"
    " it, `unconditional` if it sits on a straight path. Only `guard` is in the census, because"
    " nothing can plant a violation of a raise that reaching the call already causes. Measured"
    " at 25db524: 1,736 raise sites, of which 1,550 guards, 88 reraise, 94 unconditional, 4 loop."
    " A first version of the classifier counted any conditional anywhere in the FILE and gave"
    " 1,556 guards; restricting it to the raise's own function moved it to 1,550, so the idiom"
    " `if not valid: raise ...` really is dominant and the definition is not doing the work"
)

#: What NEVER FIRED will and will not mean. Fixed before the number is known.
COUNT_READING = (
    "a guard the whole suite never fires is a guard NO TEST PLANTS A VIOLATION OF -- condition"
    " (2), measured rather than grepped. It is NOT a defect and NOT a cannot-fire: per"
    " DOES_NOT_CLAIM, most will be right and merely unproven, and only condition (4) can tell"
    " the two apart. THE THREE DIRECTORIES MUST BE READ SEPARATELY AND THIS IS THE POINT OF"
    " REGISTERING IT NOW: scripts/ holds entry points the suite mostly does not import, so a"
    " high never-fired count there is EXPECTED and nearly uninformative; genomeos/ is the tested"
    " library, so a never-fired guard there is the interesting class; tests/ is the smallest and"
    " most surprising, because a guard in the suite that the suite never fires is one nobody"
    " wrote a planting test for although they had the file open"
)

#: And the trap this census already walked into once, so the count does not repeat it.
COUNT_POPULATION_WARNING = (
    "CASES_SELF_RED is the reason this is written down: this lane reported '0 violations' from"
    " tests/test_*.py while the guard it was auditing scanned tests/*.py and was RED. A census"
    " that measures a narrower population than the guard it audits will report that guard clean"
    " while the guard is red, and that 0 passed through two readers before the guard disproved"
    " it. So the population above is taken from the TREE, not from a plausible glob, and the"
    " directory split is declared rather than chosen once the totals are in"
)


# ========================== THE SUITE-SCALE LEG IS ABANDONED, and why -- 2026-10-12 19:20
# Condition (2) was to be measured, not grepped: sys.monitoring's RAISE event recording which of
# the guards the whole suite ever fires. That leg is abandoned. It is recorded in full because a
# census that quietly drops its central measurement is worth less than one that says it could not
# afford it.

#: What it cost, measured twice and wrong twice.
SUITE_LEG_ABANDONED = (
    "the traced suite reached 44% in 160 minutes -- about 15x the untraced suite -- while holding"
    " the MACHINE-WIDE suite lock, which every lane's push gate queues behind. The remainder"
    " would have been some four and a half hours more. Killed at 19:19 and the lock released."
    " Nothing was salvaged: atexit does not run on SIGTERM, so no partial result survived"
)

#: The cause, which is not the callback and so could not be optimised away.
WHY_IT_WAS_SLOW = (
    "enabling RAISE monitoring GLOBALLY deoptimises the interpreter for all code, so the cost is"
    " not in the callback and tuning the callback cannot remove it. A first version rebuilt"
    " code.co_lines() on every exception and ran 15x slow; each (code, offset) pair was then"
    " resolved once and remembered, a three-file benchmark showed 8.84s against 8.81s -- inside"
    " noise -- and the full suite ran AT THE SAME 15x. The dedupe was real and irrelevant"
)

#: And the reason the benchmark lied, which is this lane's third instance of one mistake.
WHY_THE_BENCHMARK_LIED = (
    "the three files it benchmarked raise few exceptions, so it measured a population that does"
    " not resemble the suite. THIS IS THE THIRD TIME THIS LANE MADE THE SAME MISTAKE IN ONE DAY:"
    " the vacuity pool that could not express the values the code requires; the violation scan"
    " over tests/test_*.py while the guard scanned tests/*.py; and now a benchmark over files"
    " that do not raise. Each time the instrument measured something narrower than its subject"
    " and each time the answer looked clean. That is one finding, not three, and it is the most"
    " transferable thing this census produced: AN INSTRUMENT MUST BE EXERCISED ON THE POPULATION"
    " IT WILL JUDGE, never on a convenient sample of it"
)

#: What a future lane should do instead, with the reading that a sample forces.
SUITE_LEG_ALTERNATIVE = (
    "a bounded sample inverts the reading and must be declared that way BEFORE running: in a"
    " partial run only the POSITIVE is sound. A guard the sample fires certainly has a planting"
    " test; a guard it does not fire is NOT MEASURED -- not even 'no violation found' -- because"
    " another test outside the sample may plant it. So a sample can establish which guards ARE"
    " exercised and can never establish the count of guards that are not, which was the"
    " question. Measuring the question as asked needs the untraced suite plus per-test"
    " instrumentation cheaper than RAISE monitoring, or a tolerated 5-hour exclusive run"
)

#: THE COUNT THAT IS SOUND, and it is static. A list and a count, as the brief asked.
STATIC_COUNT = {
    "measured_at": "HEAD a2f0fa9, by reading the working tree, which three lanes are editing",
    "tree_moves_under_the_count": "1,736 raise sites at 25db524 and 1,742 an hour later. The"
    " difference is peers committing, NOT instrument drift -- which is the practical reason a"
    " census must name its tree, and this one names a sha it read rather than a clean worktree",
    "raise_sites": 1742,
    "guards": 1556,
    "reraise": 88,
    "unconditional": 94,
    "loop": 4,
    "guards_by_directory": {"genomeos/": 781, "scripts/": 750, "tests/": 25},
    "top_guard_types": {
        "SystemExit": 617,
        "ValueError": 377,
        "BioLangError": 115,
        "RefusedError": 79,
        "FileNotFoundError": 62,
        "ApiError": 59,
    },
    "what_it_establishes": "the POPULATION of the census and nothing about any guard's state."
    " Condition (2) is unmeasured for all 1,556, so every one of them is NOT_ESTABLISHED except"
    " the single site where condition (4) was satisfied -- which is exactly what DOES_NOT_CLAIM"
    " predicted: 'condition (4) is expensive and NOT_ESTABLISHED will be the largest class'",
}


#: Standing rule, every session, from the supervisor at 19:25 on the case below. Carried verbatim.
INSTRUMENTED_RUNS_RULE = (
    "an INSTRUMENTED or traced run never takes the machine-wide suite lock. It runs in a clone,"
    " scoped to the files it needs, never the whole suite. Any run that holds the lock writes its"
    " expected duration on its board entry, and a lock held past 45 minutes is reported to the"
    " coordinator"
)

#: The case that earned it, which is this lane's.
WHY_THAT_RULE_EXISTS = (
    "this lane's traced suite held the lock from 16:25 to 19:19 at 15x speed and reached 44%."
    " Every gate on the machine queues behind that lock, so the project recorded NO COMMIT"
    " between 16:54 and 19:19 and this lane was the cause rather than a bystander. The lock did"
    " exactly what it was built to do -- serialise -- and the fault was entirely in what was"
    " allowed to hold it: a measurement, not a verdict, with no duration estimate and no"
    " threshold at which anyone would look"
)

#: And the honesty clause the first estimate failed.
DURATION_ESTIMATES_ARE_MEASURED = (
    "the estimate must be MEASURED over a sample long enough to mean something. This lane posted"
    " '~20 min left' from a 90-second window early in the run; the truth was four and a half"
    " hours. The first 3% of a pytest run is collection and startup, so a rate taken there"
    " describes the wrong thing. On the replacement run a 75-second sample showed NO progress at"
    " all and was reported as 'rate unknown' rather than as an estimate; a four-minute sample"
    " then gave 2.50%/min and an expected end of 20:03"
)


#: A correction to this lane's own wording, from the supervisor at 19:33. Additive, as the others.
VEIL_IS_INSTRUMENTATION = (
    "this lane called the bare run 'UNINSTRUMENTED, normal speed'. That is WRONG and the wording"
    " is the dangerous part. The import veil wraps every finder in sys.meta_path, which is"
    " precisely why the run is about 70% slower than the ~25-minute baseline, so it IS"
    " instrumentation and it falls under INSTRUMENTED_RUNS_RULE. Calling it uninstrumented is how"
    " a rule gets sidestepped by vocabulary rather than by argument -- the slowdown was visible,"
    " this lane even named it in the same sentence, and still filed the run under the exemption"
)

#: And the method that replaces the veil, which tests the condition instead of simulating it.
BARE_CONDITION_METHOD = (
    "from the next one on, the bare condition is tested THE WAY CI TESTS IT: a"
    " `git clone --no-local` with a fresh venv of the 13 base packages, no veil, OUTSIDE the"
    " suite lock -- as lane-cired2 did. That needs no lock because it is a separate clone with"
    " its own environment, and it needs no veil because the packages are genuinely absent rather"
    " than hidden. The veil was a way to avoid touching the shared .venv, which was the right"
    " constraint and the wrong solution: a real bare environment answers the same question"
    " without instrumenting anything, and without holding anything other lanes need"
)
