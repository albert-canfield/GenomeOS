# SPDX-License-Identifier: AGPL-3.0-or-later
"""What a miss means: overlap and distance classes over every scored frame, registered first.

The loci benchmark scores a hit when the gene a derived layer ranks first is one of the published
target's symbols, and a miss otherwise. That rule has never distinguished two things that are not
the same finding:

- the model named a gene **somewhere else**, which is the failure the benchmark exists to count;
- the model named **the right place under another name** - an antisense, host, readthrough or
  microRNA transcript whose body lies inside the published target's - which is a fact about what
  GENCODE calls a locus and not about what the model predicted.

Section 24 of docs/LOCI-BENCHMARK.md found the second case by accident. Of the 21 loci where the
repaired `loci.read_deletion` changed the deletion layer's answer, 19 are locatable and 5 promote a
gene whose body **overlaps** the published target's - IGF2-AS over IGF2, ENSG00000259006 over MC1R,
ENSG00000288879 over HMGA1, ENSG00000240739 over SLC2A3, CCDC26 over CCDC26 - with 7 more inside
100 kb. Every one is a miss today. The same phenomenon broke the fourth frame's positive control,
and the coordinator's decision of 2026-09-22 (docs/ROADMAP.md item 4) ruled that section 22 stays
interpretable, that **the hit rule stays symbol equality everywhere so that no rate moves**, and
that the general question comes here with its own registration rather than a ruling there.

**This module describes existing misses. It does not change the hit rule and it computes no rate
that replaces one.** Whether the benchmark should ever report a second, overlap-tolerant rate beside
the strict one is a question with a number attached, and the number is fixed in `PREREGISTRATION`
below before a single miss is classified.

**The trap, stated plainly: a benchmark that widens what counts as a hit after seeing its rate is a
benchmark tuning itself.** Every rate in this benchmark has fallen under scrutiny - 0.200 at n = 50,
then 0.180, then 0.172 - and a session that has just watched its headline fall and then discovers a
rule under which several misses become hits is in exactly the position where a benchmark stops
measuring anything. Three things are done about it, all before any counting: the classes are defined
in GENCODE coordinates and committed first; the thresholds that would make the case for a second
rate, and the thresholds that would close the question, are committed with them; and any tolerant
rate is reported BESIDE the strict one, labelled a description, and never instead of it.

**Cost: zero AlphaGenome requests, by construction.** Every reading classified here is already on
disk - the frames' saved results and `loci_reread`'s re-read of them - and GENCODE coordinates are a
local annotation. Nothing in this module calls a scorer, opens a client or reads an interval the
sweep has not already bought.
"""

from __future__ import annotations

from typing import Any

#: the frames whose misses are classified, with the doc section each is reported in. The first five
#: are `loci_reread`'s frames, re-read through the repaired reader; the sixth is the corrected draw,
#: which was scored through the repaired reader in the first place.
FRAMES: dict[str, str] = {
    "loci_benchmark": "9 - the seventeen panel",
    "loci_candidates": "19 - the nine candidates",
    "loci_third": "21 - the third set",
    "loci_fourth": "22 - the fourth frame, drawn by a rule",
    "loci_noncoding": "23 - the non-coding frame",
    "loci_fourth_rejoined": "24 - the fourth frame, the corrected draw",
}

#: The classes, in the order they are tested. The first that matches is the gene's class, so the
#: order is part of the definition and is fixed here with the classes themselves.
CLASSES: list[dict[str, str]] = [
    {
        "name": "exact",
        "definition": (
            "the gene the layer ranks first IS one of the published target symbols. This is a hit"
            " under the benchmark's rule, and it is in the list as the classifier's control: every"
            " strict hit must land here and a classifier that puts one anywhere else is broken"
        ),
    },
    {
        "name": "overlaps_the_target_body",
        "definition": (
            "a different symbol whose GENCODE gene body intersects a published target's gene body"
            " on the same chromosome: start < target.end and end > target.start, the same test"
            " `loci_reread.neighbourhood` uses. This is the class - and the ONLY class - that can"
            " feed the tolerant rate, because it is the one where the model named the published"
            " target's own span under a name the panel did not write down"
        ),
    },
    {
        "name": "within_10_kb",
        "definition": "no overlap; the smallest body-to-body gap to a published target is 1 to 10,000 bp",
    },
    {
        "name": "within_100_kb",
        "definition": "no overlap; the smallest body-to-body gap is 10,001 to 100,000 bp",
    },
    {
        "name": "same_ctcf_node",
        "definition": (
            "further than 100 kb, but the named gene and a published target are both inside the"
            " element's CTCF node, read off the `node` layer's own `genes_in_node` list as that"
            " layer stored it. It is tested after the distance classes so that a near gene is"
            " reported as near rather than as a node member"
        ),
    },
    {
        "name": "elsewhere",
        "definition": (
            "on the same chromosome, further than 100 kb, and not in the element's CTCF node. This"
            " is the class that means the benchmark's miss is a real miss"
        ),
    },
    {
        "name": "not_locatable",
        "definition": (
            "the named symbol has no GENCODE record on the locus's chromosome, or GENCODE uses it"
            " at more than one locus there - Y_RNA is used at 57 loci on chr2 - so a distance to"
            " whichever copy comes first is arithmetic about the wrong gene. Counted and listed"
            " apart, never folded into any other class and never into a tolerant rate"
        ),
    },
    {
        "name": "named_nothing",
        "definition": "the layer ranked no gene first, so there is no gene to place",
    },
]

PREREGISTRATION: dict[str, Any] = {
    "written": (
        "2026-09-22, and committed on its own before a single miss was classified. The counts in"
        " section 24 that motivate it - 5 overlaps and 7 more inside 100 kb among 19 locatable"
        " movements - are over the 21 loci where the REPAIRED READER CHANGED THE ANSWER, which is a"
        " different and much smaller set than the misses classified here; they are the reason to"
        " measure and not a preview of the result"
    ),
    "the_question": (
        "the benchmark scores a miss whenever the gene a derived layer ranks first is not one of"
        " the published target's symbols, and it has never measured whether a miss means `named a"
        " gene somewhere else` or `named the right place under another name`"
    ),
    "what_this_is_not": (
        "it is a DESCRIPTION of misses that already exist and NOT a new hit rule. `loci.score_locus`"
        " is not touched, no saved result is rewritten, and no rate published in sections 9, 18, 19,"
        " 21, 22, 23 or 24 moves by this work. The coordinator's decision of 2026-09-22 already"
        " fixed that the hit rule stays symbol equality everywhere; this module is bound by it"
    ),
    "what_is_classified": (
        "every locus in every frame scored through the repaired reader that is a MISS under the"
        " strict rule - `score.target_hit_derived` false - in the six frames of `FRAMES`. The"
        " strict HITS are classified too and are the control: each must come out `exact`"
    ),
    "which_gene_is_classified": (
        "the benchmark's headline hit is a union over the derived layers, so the locus's class is"
        " the BEST class, in `CLASSES` order, over the rank-1 gene of each derived layer present"
        " (deletion, eqtl, gene_input, coding). Taking the best is what makes the locus-level class"
        " comparable to a union hit; the per-layer classes are reported beside it, and the deletion"
        " layer's own class is reported separately because that layer is the model's answer about"
        " the element and is where section 24's movements happened"
    ),
    "the_classes": [c["name"] for c in CLASSES],
    "the_classes_are_fixed_first": (
        "the eight classes and the order they are tested in are defined in `CLASSES` above, in"
        " GENCODE coordinates, and committed before any locus is classified. Nothing is added,"
        " merged or re-cut after the counts are seen"
    ),
    "only_overlap_can_be_tolerated": (
        "the tolerant rate, if it is ever reported, counts a locus when a derived layer's rank-1"
        " gene is `exact` OR `overlaps_the_target_body`. The distance classes are DESCRIPTION ONLY"
        " and feed no rate, ever. The reason is not taste: a rule that admitted `within_10_kb`"
        " would score the nearest-gene baseline as a hit, and that baseline is the control the"
        " whole fourth frame was drawn to beat. A benchmark whose hit rule tolerates distance"
        " cannot then report that the model follows proximity"
    ),
    "what_would_make_the_case_for_a_second_rate": (
        "BOTH of: (a) pooled over the six frames, at least 15% of strict misses classify as"
        " `overlaps_the_target_body`; and (b) the overlap-tolerant rate exceeds the strict rate by"
        " at least 0.05 in at least one frame of n >= 50 - the fourth frame or the corrected draw -"
        " which is at least three loci there. (a) alone would be a handful of curated loci; (b)"
        " alone could be three loci out of a pooled distribution that is otherwise clean. A second"
        " permanently reported rate has to earn both"
    ),
    "what_would_close_the_question": (
        "BOTH of: (a) pooled, fewer than 5% of strict misses classify as"
        " `overlaps_the_target_body`; and (b) no frame's overlap-tolerant rate exceeds its strict"
        " rate by more than one locus. That outcome says the misses are genuinely elsewhere, the"
        " strict rule is right, and section 24's five overlaps were a property of the movement set"
        " rather than of the benchmark. It is to be reported exactly as plainly as the other"
        " outcome, and it STRENGTHENS the benchmark rather than costing it anything"
    ),
    "in_between": (
        "between 5% and 15%, or either threshold met without the other: the distribution is"
        " published, NO second rate is adopted, and the question is recorded as measured but not"
        " settled. It would then be re-asked only on a frame drawn for it, not on these"
    ),
    "the_trap": (
        "a benchmark that widens what counts as a hit after seeing its rate is a benchmark tuning"
        " itself. This session can see, before it starts, that a tolerant rule would raise numbers"
        " that have fallen three times in two days. That is why the classes are fixed before the"
        " counting, why the thresholds are fixed with them, why only overlap - never distance - can"
        " ever be tolerated, and why any tolerant rate is reported BESIDE the strict rate, labelled"
        " a description, and never in place of it. If the thresholds are met, the output of this"
        " work is a proposal to the benchmark's owner carrying its own registration, not a rate"
        " change made here"
    ),
    "the_classifier_control": (
        "every locus that is a strict hit must classify as `exact`, in the locus-level class and in"
        " the class of the layer that scored the hit. A hit that classifies as anything else means"
        " the symbol lookup, the chromosome join or the class order is wrong and NOTHING else in"
        " the output may be read. `classify_all` reports the count and the harness raises"
    ),
    "falsifiers": [
        "a strict hit classifies as anything but `exact` - the classifier is broken, stop",
        (
            "a locus classifies as `overlaps_the_target_body` at a layer whose named gene is on a"
            " different chromosome from the published target - the join is wrong"
        ),
        (
            "the counts per frame do not sum to that frame's loci - misses plus hits - so some"
            " locus was silently dropped"
        ),
        (
            "the strict rate recomputed from the classified rows differs from the rate section 24"
            " publishes for that frame - this module is reading a different scoring than the"
            " benchmark published, and the classification is about something else"
        ),
        (
            "any AlphaGenome request is spent. Every reading is already on disk and GENCODE is"
            " local; a request means this module is scoring rather than describing"
        ),
    ],
    "cost": "0 AlphaGenome requests, checked before anything was written rather than estimated after",
}
