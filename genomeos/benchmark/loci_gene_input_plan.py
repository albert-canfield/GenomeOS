"""A REGISTRATION ONLY: what a re-score under the repaired gene-input reader would be (2026-10-03).

This module holds no code that scores anything. It is committed on its own, before any locus is
re-read under `loci_gene_input.GENE_INPUT_V2`, so that the population, the denominators, the
comparison and the reading are all on the record ahead of the numbers. **It carries NO COUNT.** Not
one locus is classified by it, no frame is re-read, no rate is computed, and the fields below that
would hold results are absent rather than empty - a field that would be filled by a run does not
exist here, so nothing can be read as a preview.

Why a registration at all. The repair this plans for is the same kind of change section 24 made to
`read_deletion`, and that change was registered first for a reason that is on the record in
`loci_reread.PREREGISTRATION` and in section 25's own trap clause: "a benchmark that widens what
counts as a hit after seeing its rate is a benchmark tuning itself". The reader change committed
alongside this file is strictly narrower than that - it moves no rate at all, because v1 is the
default - but the moment anyone asks what v2 would read, they are proposing to move a published
number, and the honest order is population first, reading second, numbers last.

The reading this plans for is NOT YET AUTHORISED. The decision to run it belongs to the benchmark's
owner, and `AUTHORISATION` says so in those terms.
"""

from __future__ import annotations

from typing import Any

from genomeos.benchmark.loci_gene_input import GENE_INPUT_V1, GENE_INPUT_V2
from genomeos.benchmark.loci_miss import GRADED_FRAMES, PUBLISHED_STRICT, TOLERATED
from genomeos.benchmark.loci_miss import PREREGISTRATION as MISS_PREREGISTRATION

#: the governing registration this plan is bound by, by sha. Its rules are IMPORTED above - the
#: frames and their published strict rates, the two classes a tolerant reading may ever count, and
#: the graded frames' denominator - and are not restated with different words anywhere below.
GOVERNING_REGISTRATION = "bf33233"

PLAN: dict[str, Any] = {
    "written": (
        "2026-10-03, committed ALONE and before any locus was re-read under"
        f" {GENE_INPUT_V2!r}. This module carries no count, classifies no locus, reads no frame and"
        " computes no rate. It is a plan for a re-score, not a re-score"
    ),
    "the_question": (
        "`loci.read_gene_input` read `predicted_coding` alone until 2026-10-03 and now has a named,"
        " opt-in second reading. The question a re-score would answer is: how many loci does the"
        " gene-input layer's rank-1 gene change at, in which direction, and does the benchmark's"
        " headline - a union over the derived layers - move at all"
    ),
    "authorisation": (
        "NOT AUTHORISED BY THIS FILE. The decision to run it is the benchmark's owner's, because"
        " the output would be a movement in published numbers. Nothing in this repository calls"
        f" {GENE_INPUT_V2!r}, and the reader's default stays {GENE_INPUT_V1!r} until the owner says"
        " otherwise. A session that runs the plan below without that decision has moved a rate on"
        " its own authority"
    ),
    "the_population": (
        "every locus of the six frames `loci_miss.FRAMES` names that carries a `gene_input` reading"
        " in its committed result - all of them, hits and misses alike, with no selection on what"
        " the layer said. 147 readings exist, which is a property of the committed frames and is"
        " stated here as the population's size so that a later run cannot quietly narrow it; it is"
        " not a count of anything measured. Each frame keeps its OWN published denominator, imported"
        " from `loci_miss`: the two graded frames exclude controls and unaskable loci, every curated"
        " panel reports over every locus it chose"
    ),
    "what_would_be_re_read_and_what_would_not": (
        "ONLY the gene-input layer, and only through `read_gene_input(..., reading=v2)`. The"
        " deletion, eQTL, coding, node, reader, satmut and window readings are read from the"
        " committed frames unchanged, because re-reading them would put two movements behind one"
        " number - the objection `loci_reread.PREREGISTRATION` raised against this very change, and"
        " the one the versioning is for. The hit rule stays symbol equality. The element rows are"
        " the same rows: `loci._deletion_rows` is called with the same window and the same results"
        " directory, and a locus whose rows cannot be reconstructed identically - `elements_scored`"
        " not equal to the committed value - is REFUSED and named, never re-scored on a different"
        " row set"
    ),
    "the_rate_and_its_three_denominators": (
        "three numbers per frame, each on its own denominator and never pooled with another."
        " (a) the GENE-INPUT LAYER'S OWN rate: loci where that layer's rank-1 gene is a published"
        " target symbol, over the frame's published denominator, under v1 and under v2 side by side."
        " (b) the HEADLINE derived rate, the union over the derived layers, under v1 and v2 - and v1"
        " MUST reproduce `loci_miss.PUBLISHED_STRICT` for that frame exactly or the re-score is"
        " reading a different scoring and nothing in it may be read. (c) the count of loci whose"
        " gene-input rank-1 gene changed at all, in both directions, which is the measurement that"
        " does not depend on any hit rule"
    ),
    "the_comparison_against_the_baseline": (
        "the control is the nearest-coding-TSS-in-node rule, as it already is in `loci_miss`, and it"
        " is read from the committed frames UNCHANGED under both arms - it has no gene-input layer,"
        " so no reading of this reader can move it. That is what makes the comparison the point: the"
        " baseline is held fixed by construction, so any movement in (a) or (b) is the model's and"
        " not the measurement's. `loci_miss._headline`'s own rule travels with it and binds this"
        " plan: 'a tolerance applied to the model and not to the baseline it is measured against"
        " would flatter the model by construction, so it is applied to both or the comparison says"
        " nothing'"
    ),
    "the_prohibition_this_plan_is_bound_by": (
        f"imported verbatim from {GOVERNING_REGISTRATION}, and it binds every arm above:"
        f" {MISS_PREREGISTRATION['only_overlap_can_be_tolerated']}"
    ),
    "how_that_prohibition_constrains_this_plan_specifically": (
        "no arm of the re-score may report a rate in which a distance class counts. The only"
        f" tolerant reading that may ever appear is over {TOLERATED}, it is reported BESIDE the"
        " strict rate and labelled a description, and it is applied to the baseline in the same"
        " breath or not at all. A plan that scored the gene-input layer as a hit because its named"
        " gene is near the published target would score the nearest-gene baseline as a hit too, and"
        " that baseline is the control the whole fourth frame was drawn to beat. Said here because a"
        " lane was briefed against this prohibition on 2026-10-02 and a plan is exactly where it"
        " gets broken next"
    ),
    "expected_direction_written_before_the_run": (
        "a loss or nothing, in every frame but one, for a construction reason and not a hope."
        " abs(predicted) >= abs(predicted_coding) holds per element per head, so under v2 a gene the"
        " coding head named can be displaced and can never be promoted: a published CODING target"
        " that ranked first can only lose rank. The fourth frame and its corrected draw REQUIRE the"
        " published target to be protein coding, and three of the four other frames carry coding"
        " targets only, so in five of the six frames v2 can lower the gene-input layer's rate and"
        " cannot raise it. `loci_noncoding` (n = 2) is the one frame where a gain is possible, and"
        " section 23 already records that its gene-input layer scores 0 of 2 BY CONSTRUCTION under"
        " v1 - 'the layer reads `predicted_coding` only' - so that frame is where the repair should"
        " show first and where two loci cannot support a rate"
    ),
    "what_a_change_would_mean": (
        "four outcomes, with what each licenses, written before the numbers so that none of them can"
        " be chosen afterwards."
        " ONE, the gene-input layer's rate falls and the headline does not move: the expected"
        " outcome. It means the layer was being credited for coding-restricted answers at loci where"
        " some other derived layer carried the hit anyway, and the published headline was never"
        " resting on the defect. Section 24's and 25's rates stand and the repair is a correction to"
        " a component, which is what should be reported."
        " TWO, the headline falls: the published derived rate was resting in part on a reader that"
        " answered a narrower question than the one it was scored on. Every section quoting that"
        " frame's rate must then carry the corrected number, the owner decides whether the frame is"
        " re-reported or superseded, and the fall is reported as a fall - not averaged with anything"
        " and not offset against the two loci of `loci_noncoding`."
        " THREE, the headline rises: possible only where a published NON-CODING target gains, which"
        " is `loci_noncoding` and nothing else, and n = 2 cannot support a rate. It would be"
        " reported as two loci named, with no rate attached."
        " FOUR, nothing moves anywhere: the strongest outcome for the benchmark and the one to"
        " report just as plainly. It would mean the gene-input layer's summed window is dominated by"
        " elements whose two heads name the same gene, so the restriction was never binding on this"
        " layer - and it would close section 24's and section 25's item rather than leaving it open"
    ),
    "falsifiers_that_stop_the_run": [
        "the v1 arm does not reproduce `loci_miss.PUBLISHED_STRICT` for some frame - the re-score is"
        " reading a different scoring and nothing in it may be read",
        "a locus's element rows cannot be reconstructed with the committed `elements_scored` - it is"
        " refused and named, not re-scored on a different row set",
        "the baseline's k moves in either arm - the nearest-TSS rule has no gene-input layer and"
        " cannot move, so a movement means the join or the denominator is wrong",
        "any AlphaGenome request is spent - every reading is already on disk",
        "a reported rate counts a distance class - the prohibition above, broken",
        "a committed result, registration or figure is rewritten by the run rather than a new result"
        " written beside it",
    ],
    "what_this_plan_does_not_do": (
        "it does not score a locus, name a count, state a rate, or read a frame. It does not change"
        " the reader's default. It does not touch `loci_miss`, `loci_reread` or any committed"
        " result. It does not propose a new hit rule, and it does not reopen the question"
        f" {GOVERNING_REGISTRATION} settled about tolerating distance"
    ),
    "cost_of_running_it": (
        "0 AlphaGenome requests, no network, no money. Every element row and every committed reading"
        " is already on disk. The one real cost is memory: the element rows of most of these loci"
        " live in the git-ignored per-chromosome archives of the all-element sweep, about 30 MB of"
        " JSON each, and a run must read them ONE CHROMOSOME AT A TIME under a registered resident"
        " ceiling that raises. The reader change committed with this plan opened none of them, which"
        " is why its byte-identity evidence covers 5 loci end to end and 147 by key set and"
        " invariant rather than 147 end to end"
    ),
    "the_frames_and_their_published_strict_rates": {
        "imported": "loci_miss.PUBLISHED_STRICT, by sha " + GOVERNING_REGISTRATION,
        "rates": {k: list(v) for k, v in PUBLISHED_STRICT.items()},
        "graded": list(GRADED_FRAMES),
        "reading": (
            "these are the numbers the v1 arm must reproduce, quoted from the governing registration"
            " rather than recomputed here. Nothing in this module computed them"
        ),
    },
}
