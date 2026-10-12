# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""globscope (2026-10-12): the CLASSIFICATION RULE for the 18 `exists()`-filtered input groups.

WRITTEN AND COMMITTED BEFORE A SINGLE CALL SITE IS READ. Nothing in this module was decided after
looking at a site; that is the whole discipline of the lane and the only reason its answer is worth
anything. The 18 paths and line numbers are inherited from `data/results/globcensus.json` and no
source line of any of them has been opened at the time this file is committed.

THE DEFECT, inherited verbatim in substance from `globcensus.py` and not restated as this lane's own:
a declared input group built as a filesystem glob filtered by `.exists()` makes the declaration a
function of the MACHINE rather than of the COMMIT, so it silently shrinks when a store is partial
instead of refusing.

WHAT IS ALREADY ESTABLISHED and is carried, not re-measured:

* 82 `files_entry` calls, **18** of them `exists()`-filtered -- not 17, which is what a one-line grep
  sees; `scripts/trio_q100_phase.py:72` spans four lines with its filter on line 74. This lane uses an
  AST scan for the same reason and never a grep.
* 130 group declarations over 88 names: 82 constant, 6 varying.
* The 2026-10-11 finding was a NEGATIVE: the reader v1 DNase narrowPeak group declares 312 in ALL 13
  results carrying it, no shrink in any artefact. So **the defect is real as a mechanism in 18 places
  and so far unrealised where it was found**.
* The limit, held by a test and not by a caveat: **a group CONSTANT at a shrunken value is invisible to
  that method, so the 82 constant groups -- including the 13 that all say 312 -- are NOT shown safe.**
* Stage 2 of `globcensus` claimed "THE DEFECT IS REALISED" over six declarations and was REFUTED by
  lane-shrinkimpact (`8d620bc`, folded `9a68bf2`): **realised in 0 of 6**, and the wrong headline is
  left standing beside the correction. `STAGE_2_SUMMARY` there is "the mechanism is real in 18 places
  and realised in 0 of the 6 checked". Those six declarations fall in 2 of the 6 varying groups
  (`per_element_response_cache` x4, GENCODE 50 x2); the other 4 varying groups are untouched by it.

WHAT THIS LANE DOES: classifies each of the 18 call sites into one of five classes, and each of the 6
varying groups into one of three. It REPAIRS NOTHING. A classification that implies a fix is a finding
for the report and not an edit. `globcensus.py`, its test and `globcensus.json` are not touched.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent

#: The inherited count, pinned as the same object it is read from, never retyped as a literal here.
INHERITED_COUNT_SOURCE = "data/results/globcensus.json"
INHERITED_COUNT_FIELD = ("call_sites", "exists_filtered")


def inherited_site_count(root: Path | None = None) -> int:
    """The 18, read from `globcensus.json` itself so a later edit there cannot silently disagree."""
    base = root if root is not None else ROOT
    d = json.loads((base / INHERITED_COUNT_SOURCE).read_text())
    for key in INHERITED_COUNT_FIELD:
        d = d[key]
    return int(d)


# ------------------------------------------------------------------ the classes, fixed before looking

#: The five classes a call site may fall into. MUTUALLY EXCLUSIVE and decided in the order of
#: `DECISION_ORDER` below, so no site can be argued into a second one after the fact.
CLASSES = {
    "REFUSES": (
        "the site raises, exits or otherwise returns an error on a short set rather than filtering it"
        " away. DECIDABLE TEST, from the code alone: the function enclosing the `files_entry` call"
        " contains a `raise` or a `SystemExit`/`sys.exit` whose reachability depends on a test over the"
        " globbed collection -- its length, its emptiness, a set difference against the expected"
        " members, or a per-member `.exists()` assertion. A refusal makes the remaining questions moot,"
        " because a short set never reaches the manifest"
    ),
    "CANNOT_SHRINK": (
        "every file the UNFILTERED glob can match is tracked in git, so a missing file is a"
        " tracked-file absence that the repository's existing guards already catch, and the `.exists()`"
        " filter can only ever drop a file whose absence is already loud. DECIDABLE TEST: the glob's"
        " directory prefix resolves to a literal path from the code, the unfiltered pattern matches at"
        " least one file on disk, and `git ls-files` reports EVERY matched file as tracked. THIS MUST BE"
        " PROVEN AND NOT ASSUMED: a directory under `data/knowledge`, `data/cache` or `data/reference`"
        " is git-ignored and is never this class, however tracked-looking its name"
    ),
    "SCOPE_IS_DECLARED_ELSEWHERE": (
        "the result's own artefact states the expected scope from a source INDEPENDENT of the glob --"
        " an expected count, a chromosome list or a sample list -- so a shrink is detectable from the"
        " artefact by comparing the two. DECIDABLE TEST: the same manifest-writing function records a"
        " second entry (a parameter, a value entry or another input) whose value derives from the same"
        " constant or argument that feeds the glob and NOT from the filtered list, so that a short glob"
        " leaves the artefact internally inconsistent. A scope derived from the filtered list itself"
        " does NOT qualify, because it shrinks with it"
    ),
    "CAN_SHRINK_SILENTLY": (
        "the site filters on `.exists()` and records whatever is present, with no refusal, no tracked"
        " directory and no independent record of the expected scope. DECIDABLE TEST: the glob's"
        " directory resolves, at least one matched file is NOT tracked in git, and neither the REFUSES"
        " nor the SCOPE_IS_DECLARED_ELSEWHERE test is met. These are the live sites of the defect"
    ),
    "UNDECIDABLE_FROM_THE_CODE": (
        "the facts the four tests above need cannot be resolved statically -- the glob's directory or"
        " its expected member set is built from a runtime argument, a configuration value, or a"
        " function of another module whose return is data-dependent -- so no class can be decided"
        " without running the program or opening a store. DECLARED IN ADVANCE AS A LEGITIMATE OUTCOME:"
        " a site in this class is reported as undecidable with the fact that is missing named, and is"
        " never guessed into one of the other four"
    ),
}

#: The order in which the tests are applied. The FIRST class whose test is met is the site's class.
DECISION_ORDER = (
    "REFUSES",
    "UNDECIDABLE_FROM_THE_CODE",
    "CANNOT_SHRINK",
    "SCOPE_IS_DECLARED_ELSEWHERE",
    "CAN_SHRINK_SILENTLY",
)

#: Why UNDECIDABLE sits SECOND and not last, which is the one ordering decision that could be bent.
WHY_UNDECIDABLE_IS_SECOND = (
    "a refusal is readable from the enclosing function whatever the glob resolves to, so REFUSES is"
    " decidable first. Every one of the three remaining classes needs the glob's member set resolved:"
    " CANNOT_SHRINK needs every matchable file shown tracked, and CAN_SHRINK_SILENTLY needs at least"
    " one shown untracked. If the member set does not resolve, NONE of the three is decidable, and"
    " placing UNDECIDABLE last would let the residual class CAN_SHRINK_SILENTLY collect every site"
    " whose directory simply could not be read -- which would inflate the finding this lane exists to"
    " measure. Placing it second makes an unresolved glob say so instead of accusing the site"
)


def classify_site(
    *,
    refusal_in_enclosing_scope: bool,
    glob_root_resolved: bool,
    matched_files_on_disk: int,
    all_matched_files_tracked: bool | None,
    independent_scope_recorded: bool,
) -> str:
    """`DECISION_ORDER` as a function. Pure: it takes the facts and returns the class, nothing else.

    `all_matched_files_tracked` is `None` when the question could not be put -- the glob root did not
    resolve, or it resolved and matched nothing, in which case tracked-ness is unmeasured rather than
    false.
    """
    if refusal_in_enclosing_scope:
        return "REFUSES"
    if not glob_root_resolved or matched_files_on_disk == 0 or all_matched_files_tracked is None:
        return "UNDECIDABLE_FROM_THE_CODE"
    if all_matched_files_tracked:
        return "CANNOT_SHRINK"
    if independent_scope_recorded:
        return "SCOPE_IS_DECLARED_ELSEWHERE"
    return "CAN_SHRINK_SILENTLY"


# ------------------------------------------- the six varying groups, the rule fixed before looking

#: The three classes a VARYING group's declaration may fall into.
VARIATION_CLASSES = {
    "LEGITIMATE_SCOPE": (
        "the result's OWN declared scope accounts for the count it declares: either the result's"
        " parameters state a member list or a count equal to the declared member count, or an exclusion"
        " in that same result names every missing member. A chr21-only result declaring 1 where a"
        " genome-wide one declares 24, with the narrower scope STATED IN THE ARTEFACT, is this class"
    ),
    "SILENT_SHRINK": (
        "all four conditions of `globcensus.is_shrink_v2` hold, carried whole and not weakened: the"
        " declared count is below the result's OWN declared scope; the missing members are present in"
        " the store today; no exclusion in that result names them; AND the missing members CARRY DATA"
        " the result's own scope would have used. The fourth condition is lane-shrinkimpact's and is"
        " the one that refuted stage 2 -- without it a labelling artefact and an empty chromosome look"
        " exactly like a shrink from the manifest alone"
    ),
    "UNDECIDABLE": (
        "the artefact does not state a scope the count can be read against, or it states one but the"
        " fourth condition cannot be settled without opening the members and recomputing, which this"
        " lane does not do. IF THE ARTEFACT DOES NOT STATE ITS SCOPE THE HONEST CLASS IS UNDECIDABLE"
        " AND NOT 'PROBABLY FINE'. This class is declared in advance as a legitimate outcome and is the"
        " expected one wherever a recomputation would be needed"
    ),
}

#: The fourth condition is not this lane's to re-derive, so it names its source.
FOURTH_CONDITION_SOURCE = (
    "genomeos/attribution/globscope.py defers to globcensus.FOURTH_CONDITION and STAGE_2_REFUTED"
    " (lane-shrinkimpact, 8d620bc, folded 9a68bf2): realised in 0 of 6, the wrong headline left"
    " standing beside the correction"
)


def classify_variation(
    *,
    artefact_states_a_scope: bool,
    declared_count: int | None,
    scope_count: int | None,
    missing_members_present_in_store_today: bool | None,
    missing_members_named_in_exclusions: bool | None,
    missing_members_carry_eligible_data: bool | None,
) -> str:
    """`VARIATION_CLASSES` as a function, applied per (group, result) declaration."""
    if not artefact_states_a_scope or declared_count is None or scope_count is None:
        return "UNDECIDABLE"
    if declared_count >= scope_count:
        return "LEGITIMATE_SCOPE"
    if missing_members_named_in_exclusions:
        return "LEGITIMATE_SCOPE"
    if (
        missing_members_present_in_store_today is None
        or missing_members_named_in_exclusions is None
        or missing_members_carry_eligible_data is None
    ):
        return "UNDECIDABLE"
    if missing_members_present_in_store_today and missing_members_carry_eligible_data:
        return "SILENT_SHRINK"
    return "LEGITIMATE_SCOPE"


# -------------------------------------------------------- the outcomes, written before the scan runs

#: What will be reported in each case. Fixed now so no finding can be dressed up after the fact.
PREWRITTEN_OUTCOMES = {
    "all_18_cannot_shrink": (
        "the mechanism is INERT at every site: every glob reads a git-tracked directory, so the defect"
        " of 2026-10-11 has no live site in this repository and the existing tracked-file guards"
        " already cover it. The report would then say the mechanism is real in 18 places, realised in 0"
        " of the 6 checked, and ineffective at all 18 by construction -- and the named fix would be"
        " documentation only"
    ),
    "some_can_shrink_silently": (
        "those sites are named with path:line as the LIVE sites of the defect, with the untracked store"
        " each one globs. The report names a fix and does not perform it: the shape would be a refusal"
        " or an expected-member list recorded beside the group. No count of live sites is a count of"
        " realised shrinks -- lane-shrinkimpact found 0 of 6 realised -- and the report will not let"
        " the two be read as one number"
    ),
    "most_undecidable": (
        "reported as the answer, with the missing fact named per site. ALL UNDECIDABLE FROM THE CODE IS"
        " A COMPLETE ANSWER and not a failure of the lane: it would establish that the defect cannot be"
        " triaged statically, which is itself the finding, and it would name what a decidable version"
        " of the question requires"
    ),
    "the_six_varying": (
        "each gets one of LEGITIMATE_SCOPE, SILENT_SHRINK or UNDECIDABLE with the ARTEFACT FIELD that"
        " decided it quoted. UNDECIDABLE is expected for any group whose fourth condition needs a"
        " recomputation, and is reported as such rather than softened"
    ),
}

#: What this lane does not claim, written before the scan so a later reader cannot widen it.
DOES_NOT_CLAIM = (
    "not that any committed result is wrong; not that a CAN_SHRINK_SILENTLY site has ever shrunk, which"
    " is a different question and was answered 0 of 6 by lane-shrinkimpact; not that a CANNOT_SHRINK"
    " site is safe in any sense beyond the one its test establishes, namely that its absence would be a"
    " tracked-file absence; not anything about the 82 CONSTANT groups, which remain NOT SHOWN SAFE"
    " because a group constant at a shrunken value is invisible to the census method; and not anything"
    " about the 64 `files_entry` calls that carry no `exists()` filter, which are outside the 18"
)

#: The limit inherited from globcensus, carried word for word and never strengthened.
INHERITED_LIMIT = (
    "a group CONSTANT at a shrunken value is invisible to that method, so the 82 constant groups --"
    " including the 13 that all say 312 -- are NOT shown safe"
)

#: The inherited reading, carried word for word.
INHERITED_READING = "the defect is real as a mechanism in 18 places and so far unrealised where it was found"

#: The boundary every negative in this lane's report must carry in the same sentence.
NEGATIVE_BOUNDARY = (
    "18 sites, from the AST scan of genomeos/ and scripts/ at this lane's commit, over the 6 varying"
    " groups of the 130 group declarations in data/results/*.json -- and nothing wider"
)

#: The 2026-09-27 stop list, read before this lane was briefed. Its six items are the node model,
#: no more biosamples, nothing built in area J, no AlphaGenome requests on the 882 blocks, area E's
#: fate ceiling, and BioForge's confidence calibration. THIS LANE TOUCHES NONE OF THEM: it is a static
#: classification of call sites in committed code, makes 0 model requests, spends no money, opens no
#: new data source and adds no biosample.
STOP_LIST_READ = "2026-09-27, six items, none touched: static classification only, 0 model requests"
