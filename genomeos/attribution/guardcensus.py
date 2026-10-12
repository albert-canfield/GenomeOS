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
