# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""wormsource: is the worm's between-gene lead a property of the WORM, or of Serizay et al. 2020?

THE WHOLE THREE-SPECIES OUTLIER RESTS ON ONE SOURCE. The worm is the outlier of three on these
sources, and every between-gene figure it contributes comes from a single study. Nothing has tested
whether the lead belongs to the organism or to that study. This lane does, on an independent
source in the same species.

**Durham et al. 2021** (Genome Research 31:1952, gr.271791.120) is sci-ATAC-seq in L2 *C. elegans*
built deliberately TO MATCH available sci-RNA-seq, tissue by tissue. It is a second, independent
C. elegans source beside Serizay et al. 2020, which is the worm's registered one.

**IT CHANGES NEITHER THE MECHANISM NOR THE LINEAGE AND CANNOT ADDRESS THE FOURTH-GENOME QUESTION
AT ALL.** Durham is the same species, so it has the same operons, the same SL trans-splicing and
the same ancestry. It cannot separate "the worm" from "operons and trans-splicing"; that needs a
genome that keeps one and drops the other, and `simpler_genome4` records that no public source for
such a genome was found by its six queries. This lane answers a DIFFERENT question and a result
here must not be read as answering that one.

**THE CONFOUNDS, AND THERE ARE THREE, NOT ONE.** The tissue set is the first: if Durham's tissues
differ from Serizay's, a difference in the reading is explained by TISSUE SET rather than by source.
The rule below fixes that -- compute on the INTERSECTION of tissues present in both sources and
RE-READ SERIZAY ON THAT SAME INTERSECTION.

**But that is NECESSARY AND NOT SUFFICIENT, and an earlier version of this registration wrongly said
the two readings would then differ in source alone.** They do not. **Durham 2021 is L2 LARVAE and
SINGLE-CELL** (sci-ATAC / sci-RNA); **Serizay 2020 is YOUNG ADULTS and BULK TISSUE** -- the worm's
own registered endpoint reads "the attributed gene's TPM (RNA_metrics, young adults)". The tissue
intersection controls the tissue set; it controls NEITHER STAGE NOR ASSAY. So what this lane
compares is **two sources that also differ in stage and in assay**, and that compound is what it
tests. Every outcome below says "source, stage and assay together", and none of them may attribute
a result to source alone.

**ORDER, and it is this project's own recorded trap:** code, then the registration, then the
result, then the tests that read the result. A record and the tests that read it are never one
commit, because `must_be_present` correctly raises while a path is untracked and a test committed
first reddens every sha in between. Nothing is fetched before the registration is committed; a
registration written after a fetch is worthless however honest its author, and git is what makes
it honest.

**THE 2026-09-27 STOP LIST was read and this lane touches none of its six items**: nothing more on
the node model until wave 0's measured test lands (met since 2026-09-27); no more biosamples;
nothing built in area J; no AlphaGenome requests on the 882 blocks; area E's fate ceiling is not
chased; BioForge's confidence is not calibrated. The biosample clause's subject is the human
reader's 13 biosamples, `READER_TERMS` in `genomeos/genome/context_evidence.py`, and **this lane
adds none, names none and widens none**. A C. elegans tissue in a comparative-genomics population
is not a human reader biosample. No model request, no money, no person's genome.
"""

from __future__ import annotations

import re

from genomeos.attribution import simpler

# ------------------------------------------------------------------- 1. the question, and its limit

#: The question, as asked.
QUESTION = (
    "is the worm's between-gene lead a property of THE WORM, or of SERIZAY ET AL. 2020? The"
    " entire three-species outlier rests on one source and nothing has tested that dependency"
)

#: The key every result of this lane must carry. Same mechanism as simpler_genome4's.
LIMIT_KEY = "what_this_cannot_reach"

#: The limit, fixed here so a result cannot restate it more weakly.
NOT_THE_FOURTH_GENOME_QUESTION = (
    "DURHAM IS THE SAME SPECIES, SO IT CHANGES NEITHER THE MECHANISM NOR THE LINEAGE AND CANNOT"
    " ADDRESS THE FOURTH-GENOME QUESTION AT ALL. It has the same operons, the same SL"
    " trans-splicing and the same ancestry as Serizay's worms. Separating 'the worm' from"
    " 'operons and trans-splicing' needs a genome that keeps one and drops the other, and no"
    " public source for one was found. A result here answers only whether the lead survives a"
    " change of SOURCE within C. elegans"
)


def result_carries_its_limit(result: dict) -> bool:
    """Whether a result carries :data:`LIMIT_KEY` with the limit's own wording.

    Reused from `simpler_genome4` on the coordinator's instruction, and for its reason: the test
    that matters is not the one proving a correct result passes but the one proving a WEAKENED
    RESTATEMENT FAILS, because a softened sentence is what a summary actually produces.
    """
    text = result.get(LIMIT_KEY)
    if not isinstance(text, str):
        return False
    return (
        "CANNOT ADDRESS THE FOURTH-GENOME QUESTION AT ALL" in text
        and "only whether the lead survives a change of SOURCE" in text
    )


# -------------------------------------------------- 2. the replication target, both figures, pinned

#: The record both figures are read from, so they are pinned to it rather than retyped.
SOURCE_OF_RECORD = "docs/SIMPLER-GENOMES.md"

#: THE PRIMARY, named primary. Carried word for word.
PRIMARY_TARGET = "+3.73 [+2.79, +4.74]"

#: Its population, from the same table.
PRIMARY_ELEMENTS = 2335

#: The figure surviving the worm's own outron and trans-splicing caveat. NOT the primary.
OUTRON_SURVIVING_TARGET = "+4.47 [+3.16, +5.78]"

#: Its population.
OUTRON_SURVIVING_ELEMENTS = 1745

#: Which is which, stated so a later reader cannot promote the larger number.
WHICH_IS_PRIMARY = (
    "the PRIMARY is +3.73 [+2.79, +4.74] on 2,335 between-gene elements. +4.47 [+3.16, +5.78] on"
    " 1,745 is the figure that SURVIVES the outron and trans-splicing exclusion and is secondary;"
    " the worm's own lane did not call it a strengthening, because the intervals overlap, and"
    " this lane does not either"
)


def targets_in_the_record(text: str) -> bool:
    """Whether both carried figures are still in the record's text, with their populations."""
    return all(token in text for token in (PRIMARY_TARGET, OUTRON_SURVIVING_TARGET, "2,335", "1,745"))


# ------------------------------------------------------- 3. the same endpoint, not a similar one

#: Carried from the shared module AS THE SAME OBJECT, so it cannot drift from the worm's own.
MIN_DISCORDANT = simpler.MIN_DISCORDANT

#: And the rule if Durham's input will not fit the endpoint.
ENDPOINT_IS_NOT_ADJUSTED = (
    "the same code path, the same discordance definition and the same MIN_DISCORDANT. IF THE"
    " ENDPOINT HAS TO CHANGE AT ALL TO ACCEPT DURHAM'S INPUT, THAT IS A FINDING TO REPORT BEFORE"
    " RUNNING, NOT A CHANGE TO MAKE. An endpoint adjusted to fit a second source cannot test"
    " whether the first source's reading survives it"
)

# ----------------------------------------- 4. the tissue-set confound, fixed before any value

#: The rule, fixed before the data is looked at.
TISSUE_INTERSECTION_RULE = (
    "compute on the INTERSECTION of tissues present in BOTH sources; report the intersection's"
    " size, and name which tissues each source contributes that the other does not; and RE-READ"
    " SERIZAY ON THAT SAME INTERSECTION, so the two readings differ in SOURCE ALONE. A comparison"
    " of Durham's tissues against Serizay's full five would confound source with tissue set, and"
    " the lane would answer nothing while appearing to answer something. THIS IS NECESSARY AND"
    " NOT SUFFICIENT: it controls the tissue set and controls neither stage nor assay, so even"
    " with it the comparison remains source-plus-stage-plus-assay"
)

#: And the complete answer declared in advance for the case where it cannot be done.
COMPARISON_CANNOT_BE_MADE_IS_COMPLETE = (
    "if the tissue intersection yields fewer than MIN_DISCORDANT elements getting different genes"
    " under the two rules, the honest answer is that THE COMPARISON CANNOT BE MADE. That is"
    " declared here, before any value is read, as a COMPLETE ANSWER -- not a failure, and not a"
    " licence to widen the tissue set until the count is reached"
)


#: The three differences between the two sources, AS FIELDS, so a later reader sees all three
#: without reading a sentence. The coordinator's instruction, after the supervisor caught that an
#: earlier version of this registration claimed the tissue intersection isolated source.
THREE_DIFFERENCES = {
    "stage": {"serizay_2020": "young adults", "durham_2021": "L2 larvae"},
    "assay": {"serizay_2020": "bulk tissue", "durham_2021": "single-cell (sci-ATAC / sci-RNA)"},
    "source": {"serizay_2020": "Serizay et al. 2020", "durham_2021": "Durham et al. 2021"},
}

#: Stated plainly, because it is what the whole lane is allowed to conclude.
WHAT_IS_COMPARED = (
    "two sources THAT ALSO DIFFER IN STAGE AND IN ASSAY. The comparison tests that compound and"
    " cannot decompose it: no outcome here attributes anything to source alone, and the tissue"
    " intersection does not change that -- it removes a FOURTH difference, the tissue set"
)

# -------------------------------------- 5. the outcomes, all pre-written, before any RNA value

#: Every outcome, with its criterion, and whether it falsifies or concludes anything.
OUTCOMES = {
    "REPLICATED": {
        "criterion": "Durham's delta interval has its LOWER BOUND > 0 on the intersection, AND"
        " it OVERLAPS Serizay's interval re-read on the same intersection",
        "concludes": "the lead holds across SOURCE, STAGE AND ASSAY TOGETHER. That is a MORE"
        " ROBUST agreement than a same-stage, same-assay replication would be, because it"
        " survived three changes at once -- AND IT STILL DOES NOT ISOLATE WHICH FACTOR IT"
        " SURVIVED, so the positive must not be over-read either",
        "falsifies": "nothing; it is the registered expectation",
    },
    "NOT_REPLICATED": {
        "criterion": "Durham's interval has its UPPER BOUND <= 0, or the two intervals DO NOT OVERLAP",
        "concludes": "the lead does not hold on an L2 single-cell source; whether that is"
        " source, stage or assay is not separable here",
        "falsifies": "the reading that the lead is a property of the worm AS MEASURED BY SERIZAY"
        " 2020. Every downstream statement about the worm must then carry this outcome. It is"
        " written here before any value is read so that it cannot be softened afterwards -- and"
        " ALSO so that it cannot be STRENGTHENED into an attribution this comparison cannot"
        " make: stage and assay differ too, so no outcome here names source on its own",
    },
    "INCONCLUSIVE": {
        "criterion": "Durham's interval CONTAINS 0 and also overlaps Serizay's",
        "concludes": "NOTHING. The two readings are compatible with each other and with no"
        " difference at all, which is a real outcome and not a weak version of REPLICATED",
        "falsifies": "nothing",
    },
    "COMPARISON_CANNOT_BE_MADE": {
        "criterion": "the tissue intersection's discordant count is below MIN_DISCORDANT",
        "concludes": "nothing about the worm; only that this source cannot test it",
        "falsifies": "nothing",
    },
}

#: The registered NOT_REPLICATED sentence, carried word for word from the coordinator so it cannot
#: be paraphrased into an attribution.
NOT_REPLICATED_WORDING = (
    "the lead does not hold on an L2 single-cell source; whether that is source, stage or assay"
    " is not separable here"
)


def result_attributes_only_to_source(text: str) -> bool:
    """Whether a NOT_REPLICATED reading wrongly blames SOURCE alone.

    The refusal that matters, in the shape of :func:`result_carries_its_limit`: the plausible
    weakening a summary produces here is not a softening but a SHARPENING -- "the lead is
    source-dependent" is shorter, more quotable and says more than the comparison can support,
    because stage and assay differ too. True means the text must be rejected.

    **ITS LIMIT, stated rather than discovered later: a string match cannot tell a CLAIM from a
    PROHIBITION of that claim.** This function flagged an early draft of
    ``OUTCOMES["NOT_REPLICATED"]["falsifies"]`` that said the result must not be strengthened
    *into* "source-dependent" -- text forbidding the attribution, not making it. The draft was
    reworded rather than the check loosened, because a checker that learned to ignore the phrase
    in some contexts would stop catching it in the one context that matters. So this is a
    heuristic over wording, it is deliberately crude, and a result that discusses the phrase
    rather than asserting it must name stage and assay in the same text to pass.
    """
    lowered = text.lower()
    names_source = "source-dependent" in lowered or "depends on the source" in lowered
    names_the_others = "stage" in lowered and "assay" in lowered
    return names_source and not names_the_others


WHICH_OUTCOMES_CONCLUDE_NOTHING = ("INCONCLUSIVE", "COMPARISON_CANNOT_BE_MADE")

# ------------------------------------------------------------------------------ 6. multiplicity

MULTIPLICITY = (
    "ONE pre-specified comparison on ONE pre-specified endpoint: the between-gene stratum, delta"
    " = concordance(domain rule) - concordance(nearest TSS), on the tissue intersection. No other"
    " interval in any output is a result of this lane, and nobody may later read a chosen one"
    " into it"
)

# ------------------------------------------------------------------------- 7. order, and 8. data

ORDER = (
    "code, then the registration, then the result, then the tests that read the result. A record"
    " and its tests are never one commit: must_be_present correctly raises while a path is"
    " untracked, so a test committed first reddens every sha in between -- measured on"
    " lane-scalehet, 3 failed and 7,258 passed on a sha inside a two-minute gap. Nothing is"
    " fetched before the registration is committed"
)

DATA_TERMS = (
    "free only, streamed, every input declared with a RELATIVE input_entry path because an"
    " absolute one makes the result unrebuildable elsewhere. No purchase, no key, 0 model"
    " requests. IF THE DATA NEEDS A LOGIN OR A TRANSFER THAT COSTS ANYTHING, STOP AND REPORT --"
    " that is Albert's decision and neither this lane's nor the coordinator's"
)


def stated_min_discordant() -> int:
    """The MIN_DISCORDANT named in :data:`COMPARISON_CANNOT_BE_MADE_IS_COMPLETE`'s own sentence.

    The boundary-cannot-drift-from-its-evidence pattern, which the coordinator adopted after this
    lane's sibling wrote NINE where its records summed to SIX. Here the sentence names the
    constant by name rather than by value, which is the strongest form: there is no number in the
    prose to drift.
    """
    if not re.search(r"\bMIN_DISCORDANT\b", COMPARISON_CANNOT_BE_MADE_IS_COMPLETE):
        raise ValueError("the sentence no longer names MIN_DISCORDANT")
    return MIN_DISCORDANT
