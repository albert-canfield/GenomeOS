# SPDX-License-Identifier: AGPL-3.0-or-later
"""The four repressions S4 judged wrong, traced to where their direction came from (item 12, S4 follow-up).

    uv run python scripts/repression_trace.py

S4 (`data/results/attribution_correctness.json`, `genomeos/attribution/correctness.py`) judged every
`represses_target` claim any held evidence reached. There were four, and all four came out `incorrect`
with the detail "the gene moved the other way". Each is judged in a cell other than the one the claim
names, and each carries `refutable: false`. This script traces each claim back to the compiled rule, the
run row and the cached model answer that produced it, and traces the deciding observation back to its raw
benchmark row. It does the same for the one `activates_target` claim S4 judged wrong (CCND1, HCT116), as a
comparison.

An internal development trace of existing claims. It is not a validation, it builds no attribution
machinery and no scorer, and it changes no label and no verdict. It makes no model request. It reads the
per-element response cache (the model's answers already on disk) and two rows of the CRISPRi benchmark's
held-out file, only to trace verdicts S4 already decided with them. Nothing is fitted, selected or
relabelled from either. Writes data/results/repression_trace.json.
"""

from __future__ import annotations

from typing import Any

NAME = "repression_trace"

# ==================================================================================================
# The registration (2026-09-29, lane-repress). The cause classes and the rule that picks the primary
# class are fixed here, before any trace code was written and before the result was computed.
# ==================================================================================================
REGISTERED = "2026-09-29"
STATUS = "an internal development trace of existing claims, not a validation"
WRITTEN_AFTER = (
    "The brief asked for the classes before any tracing. This lane had already read the five compiled "
    "rules, their rows in the whole-chromosome run, the five cached model answers and the matching raw "
    "CRISPRi rows when it wrote them down, so this registration is not blind to the data. Classes (a) to "
    "(e) are the brief's five, with its definitions. Class (f) and the primary-class rule were "
    "written after that first reading, and are marked as such."
)
QUESTION = (
    "for each represses_target claim S4 judged incorrect (four) and the one activates_target claim it judged "
    "incorrect (CCND1), where did the claimed direction come from, and why does it disagree with the "
    "observation that decided the verdict"
)

#: the claims traced, exactly as S4's result states them (element, locus, gene, the claim's cell, and the
#: deciding observation's source, gene, cell and split)
CASES = (
    {
        "element": "EH38E2447945",
        "locus": "chr6:14095446-14095736",
        "value": "represses_target",
        "gene": "CD83",
        "cell": "K562",
        "deciding": {"source": "crispri:Nasser2021", "gene": "CD83", "cell": "GM12878", "split": "heldout"},
    },
    {
        "element": "EH38E3896696",
        "locus": "chr9:97960870-97961216",
        "value": "represses_target",
        "gene": "HEMGN",
        "cell": "Brain_Cerebellar_Hemisphere",
        "deciding": {"source": "crispri:Gasperini2019", "gene": "HEMGN", "cell": "K562", "split": "training"},
    },
    {
        "element": "EH38E3426791",
        "locus": "chr20:31608511-31608860",
        "value": "represses_target",
        "gene": "ID1",
        "cell": "Whole_Blood",
        "deciding": {"source": "crispri:Gasperini2019", "gene": "ID1", "cell": "K562", "split": "training"},
    },
    {
        "element": "EH38E3940111",
        "locus": "chrX:103233290-103233585",
        "value": "represses_target",
        "gene": "BEX4",
        "cell": "hair_follicular_keratinocyte",
        "deciding": {"source": "crispri:Gasperini2019", "gene": "BEX4", "cell": "K562", "split": "training"},
    },
)
#: the comparison case: not one of the four
COMPARISON = {
    "element": "EH38E1548926",
    "locus": "chr11:69588295-69588640",
    "value": "activates_target",
    "gene": "CCND1",
    "cell": "CD8_positive__alpha_beta_memory_T_cell",
    "deciding": {"source": "crispri:HCT116", "gene": "CCND1", "cell": "HCT116", "split": "heldout"},
}

SIGN, CELL, CONTRADICTION, INDIRECT, JUDGE, SPLIT = "a", "b", "c", "d", "e", "f"
PRESENT, ABSENT, UNKNOWN = "present", "absent", "unknown"
CLASSES = {
    SIGN: {
        "name": "sign or convention error",
        "definition": "a sign flipped, or a convention misread, somewhere between the source and the label, "
        "or between the raw measurement and the observation. On the label side the steps are: the model's "
        "log2 fold change (alternate over reference, where the alternate allele is the deletion), the "
        "choice between the gene's largest drop and its largest rise (enhancer_target.predict_target), the "
        "rule verb the compiler writes (activates or inhibits), and the claim value S4 reads from the verb. "
        "On the measurement side the steps are: the benchmark's EffectSize and Significant columns, the "
        "pair's outcome (measured.CrispriPair.outcome), the observation kind (correctness.CRISPRI_KIND) "
        "and the S4 table's reading of that kind on the claim's value",
        "decided_by": "step by step: present when a step's stored sign or direction is not the one the "
        "previous step's value implies under that step's own stated convention. The label side is "
        "recomputed from the cached answer with predict_target and the compiler's verb rule, and compared "
        "with the stored run row, the compiled rule line and the S4 claim. The measurement side is "
        "recomputed from the raw benchmark row and compared with the observation S4 recorded",
        "added_after_first_reading": False,
    },
    CELL: {
        "name": "cell mismatch",
        "definition": "the claim's stated cell and the deciding observation's cell differ, so by the S4 "
        "table the observation may not bear on the claim at all. The table's CRISPRi activity and context "
        "cells speak of 'the cell screened'",
        "decided_by": "present when no deciding observation's cell equals the claim's cell after S4's own "
        "normalisation (correctness.norm_cell)",
        "added_after_first_reading": False,
    },
    CONTRADICTION: {
        "name": "genuine contradiction",
        "definition": "the same reading in the same cell has the opposite direction: the model's own value "
        "for the named gene on the measured cell's track has the opposite sign to the measured effect",
        "decided_by": "the cached answer's value for the gene on the measured cell's track: `by_cell`, and "
        "the largest drop or rise where the track it came from is named for that cell. Present when a value "
        "for that cell is opposite in sign to the measurement. Absent when every such value has the "
        "measurement's sign. Unknown when the cache holds no value for that cell. The cache keeps one "
        "`by_cell` value per cell name, so where a cell has several tracks the others are not seen, and "
        "this is stated beside each case",
        "added_after_first_reading": False,
    },
    INDIRECT: {
        "name": "indirect or ambiguous",
        "definition": "another gene may carry the effect: another gene in the window, a shared promoter or "
        "a multi-target element",
        "decided_by": "present when any of: the model's any-gene prediction names a different gene from its "
        "coding prediction; the run's compact verdict says 'another gene in the same domain'; or the same "
        "screen changed another gene significantly on silencing this element (a multi-target element). "
        "Each of the three is reported separately, and the annotated genes nearest the element are listed "
        "as description only",
        "added_after_first_reading": False,
    },
    JUDGE: {
        "name": "judge inconsistency",
        "definition": "the verdict `incorrect` alongside `refutable: false`: the judge gave a verdict its "
        "own registered rules (correctness.TABLE, correctness.AXIS_RULE) do not give, or its rules "
        "contradict each other on this claim",
        "decided_by": "correctness.verdict_of re-run on the claim with the observations S4 held for the "
        "element must give the committed verdict and the committed deciding observations. What `refutable` "
        "means is read from correctness.py and stated. Separately, the committed S4 result is counted per "
        "axis by verdict and by `refutable`, to show where `incorrect` can occur with `refutable` false",
        "added_after_first_reading": False,
    },
    SPLIT: {
        "name": "direction split across tracks",
        "definition": "the model moves the named gene both ways on different tracks, and the label takes "
        "the direction and the cell of the single most extreme track (enhancer_target.predict_target: the "
        "larger of the largest drop and the largest rise over every track wins)",
        "decided_by": "present when the cached answer records both a drop and a rise for the named gene "
        "(max_drop_log2fc < 0 < max_rise_log2fc). The live scorer records only values beyond 0.05, so both "
        "being non-zero means both passed that recording threshold. The winning margin is reported",
        "added_after_first_reading": True,
    },
}
PRIMARY_RULE = (
    "written after the first reading. The primary class is the first that applies, in this order: (a) "
    "when present, because a flipped sign explains any disagreement; (e) when present; (c) when present, "
    "because the disagreement then survives aligning the cells; (b) when present and (c) is absent or "
    "unknown, because the claim and the observation are then in different cells and the model's own "
    "reading of the measured cell does not contradict the measurement, or cannot be read. (d) and (f) are "
    "recorded beside and are never primary: S4 matches the deciding observation to the named gene, so "
    "another gene cannot be what disagrees, and (f) says where a direction came from, not why a "
    "measurement disagrees with it"
)
NO_CHANGE = (
    "no label and no verdict is changed by this trace. A bug, if one is found, is fixed with a test that "
    "fails before the fix. When the fix would move a committed result's figures, the before-and-after "
    "counts are written and the result is not regenerated: that is the coordinator's call"
)
MAY_BE_CALLED = (
    "an internal development trace of five existing claims: where each claimed direction came from and "
    "which registered class explains its disagreement with the observation that decided it"
)
MAY_NOT_BE_CALLED = (
    "a validation of the labels, the model or the judge",
    "a rate: five claims are traced, chosen because S4 judged them wrong, and no share of anything is "
    "estimated from them",
    "evidence about claims S4 did not judge",
)


def registration() -> dict[str, Any]:
    """The registered constants, as the result file states them."""
    return {
        "registered": REGISTERED,
        "status": STATUS,
        "written_after": WRITTEN_AFTER,
        "question": QUESTION,
        "cases": list(CASES),
        "comparison": COMPARISON,
        "classes": CLASSES,
        "primary_rule": PRIMARY_RULE,
        "no_change": NO_CHANGE,
        "may_be_called": MAY_BE_CALLED,
        "may_not_be_called": list(MAY_NOT_BE_CALLED),
    }


def primary(classes: dict[str, str]) -> str:
    """The primary class under PRIMARY_RULE, from each class's reading (present, absent or unknown)."""
    for k in (SIGN, JUDGE, CONTRADICTION):
        if classes.get(k) == PRESENT:
            return k
    if classes.get(CELL) == PRESENT:
        return CELL
    raise ValueError(f"no registered class is primary for {classes}")


def main() -> None:
    raise SystemExit("registration only; the trace follows in its own commit")


if __name__ == "__main__":
    main()
