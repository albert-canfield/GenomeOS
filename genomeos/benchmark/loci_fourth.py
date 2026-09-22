# SPDX-License-Identifier: AGPL-3.0-or-later
"""A fourth set of published enhancer-gene loci, drawn by a stated rule rather than chosen by eye.

docs/LOCI-BENCHMARK.md now holds three sampling frames, and section 21 measured the thing they have
in common. The derived target rate reads 15/17 (0.882), 8/9 (0.889) and 8/9 (0.889) - a spread of
0.007 across three sets curated on different days by different people for different reasons. The
nearest-gene rule over the same three reads 0.647, 0.333 and 0.778, a spread of 0.445. The headline
is stable and the thing it is supposed to be measured against is not.

Section 21 also said why the baseline moves: at the third set seven of nine published targets are the
element's own nearest coding TSS, because an element sitting in its target's intron or promoter is
what a clean perturbation experiment usually looks like. So the benchmark's three frames differ
mostly in *geometry*, and nobody has asked the question that geometry raises:

  **is the derived rate a property of the model, or of the fact that most curated elements sit next
  to their targets?**

A fourth frame chosen by hand cannot answer that, because whoever picks the loci picks the geometry.
This module therefore does not contain a list of loci. It contains the **rule that draws them**, so
the frame is fixed before anybody has seen which loci it returns: every published pair is taken from
one external file in a fixed order, the rule is arithmetic over GENCODE, and there is no step at
which a locus is kept or dropped on how it is expected to score.

The rule, in one sentence: **every element in the ENCODE CRISPR benchmark's held-out arm whose
regulated target is NOT its own nearest coding TSS**. That is the geometry on which the shortcut is
wrong by construction, which is what makes it worth a day of requests - and it is also why the
nearest-coding-TSS-anywhere baseline is reported here as arithmetic rather than as a measurement.
The benchmark's own heuristic layer (`loci.read_node`, the nearest coding TSS *inside the CTCF node*)
is a different rule and is a genuine reading at these loci; the two are kept apart below and in the
section this run writes.

Nothing here is a new scorer. Every reading, hit rule, denominator and matched control comes from
`genomeos.benchmark.loci`, called with a different panel, exactly as `loci_candidates` and
`loci_third` do. The reach filter is `loci.read_reach` - a gene BODY test against the scorer's
1,048,576 bp input, never a TSS-distance test.

`PREREGISTRATION` and `assess` were written and committed before the file was read, before any locus
was known and before any request was spent. A negative is the outcome, not a reason to re-cut.
"""

from __future__ import annotations

import contextlib
from pathlib import Path
from typing import Any

from genomeos.attribution import crispri
from genomeos.benchmark import loci, loci_candidates, loci_third
from genomeos.benchmark.loci import Expect
from genomeos.results import RESULTS_DIR, load_result, save_result

NAME = "loci_fourth"
#: deletions of intervals this frame states because no registry drew one over the perturbed element.
INTERVALS = "loci_fourth_intervals"
#: this frame's own distilled GTEx hits. `loci.stream_gtex` clears a directory whose window manifest
#: changed, so two panels sharing one directory would delete each other's cache.
GTEX_DIR = Path("data/knowledge/loci_fourth")
#: the file the frame is drawn from: the held-out arm of the ENCODE CRISPR benchmark. NOT the K562
#: training arm, on which this project has already fitted logistic weights (data/results/
#: crispri_benchmark.json, 2026-09-16). Using the arm those weights were fitted on would put a fitted
#: threshold inside the loci that test it.
SOURCE = crispri.HELDOUT
#: the most loci the frame may contain, applied in genome order after the rule, so the cap can never
#: act as a quality filter. Set at 24 before the file was read; raised to 60 on 2026-09-21, after the
#: first 24 had been scored, to take in the 36 the first cap left undrawn. See `CAP_RAISE`: the rule
#: that decides which elements pass is untouched, and 60 is every element that passes it, so the cap
#: no longer binds at all and cannot be re-cut later on anything that has been seen.
MAX_LOCI = 60
#: how many the first draw took, in genome order: chr1 to chr10. Kept so the run can report the two
#: halves apart, because genome order is not random with respect to gene density.
FIRST_DRAW = 24
#: the most model requests this run may spend. A locus whose element the finished sweep has already
#: deleted costs nothing; the rest cost one each, and loci past the budget are dropped from the END of
#: genome order, with the drop reported.
REQUEST_BUDGET = 30
#: how far a candidate element must stay from every locus in the first three frames, so the fourth
#: frame is independent of them and the matched control draw stays clean. `loci.NEGATIVE_KEEP_OUT`.
KEEP_OUT = loci.NEGATIVE_KEEP_OUT
#: half the locus window a drawn locus gets, as in the other three frames.
WINDOW_HALF = 500_000
#: what each CRISPRi cell line counts as in this benchmark's vocabulary, fixed before the draw.
#: Only K562 and GM12878 have a counterpart among `loci.READER_CELLS`; the other three lines have
#: none, and no GTEx tissue is claimed for a carcinoma or an iPSC line because that would hand the
#: cell axis a hit the reader never earned. Where the counterpart is empty the cell axis is not
#: judgeable at that locus, and that is registered rather than discovered.
CELL_COUNTERPART: dict[str, dict[str, tuple[str, ...]]] = {
    "K562": {
        "cells": ("K562",),
        "gtex_tissues": (),
        "tissues": ("K562 chronic myelogenous leukaemia", "erythroid"),
    },
    "GM12878": {
        "cells": ("GM12878",),
        "gtex_tissues": ("Cells_EBV-transformed_lymphocytes",),
        "tissues": ("GM12878 lymphoblastoid", "B lymphocyte"),
    },
    "HCT116": {"cells": (), "gtex_tissues": (), "tissues": ("HCT116 colorectal carcinoma",)},
    "Jurkat": {"cells": (), "gtex_tissues": (), "tissues": ("Jurkat T-lymphoblast",)},
    "WTC11": {"cells": (), "gtex_tissues": (), "tissues": ("WTC11 induced pluripotent stem cell",)},
}

# ------------------------------------------------- the cap raise, registered before the second run
#: The second registration, written and committed on 2026-09-21 BEFORE the 36 undrawn loci were
#: scored and before the three requests they cost were spent. It is small, because raising a cap is
#: a small thing next to stating a rule - but it is a registration, and it is here rather than in a
#: message because the first 24 have been seen and this run has to be checkable against what was
#: expected of it before it ran.
CAP_RAISE: dict[str, Any] = {
    "written": (
        "2026-09-21, after commit 188f4ce scored the first 24 and before the remaining 36 were drawn,"
        " read or scored. `assess` and `select` are byte-for-byte the rule that drew the first 24:"
        " the only change is MAX_LOCI, 24 to 60"
    ),
    "what_changes": (
        f"MAX_LOCI {FIRST_DRAW} -> {MAX_LOCI}. 60 is the number of held-out elements that PASS the"
        " rule, so after this the cap does not bind and there is nothing left for a later hand to"
        " choose. n goes from 20 graded to 50 graded: 60 drawn, 10 of which the free reach filter"
        " kills before a request exists"
    ),
    "what_does_not_change": (
        "the rule, the hit rules, the readers, the denominators, the keep-out, the window, the cell"
        " counterparts and the positive control. Nothing that decides which elements are in the frame"
        " is touched, and nothing that decides whether a reading is a hit is touched. Changing a"
        " selection rule after seeing its first result is the artefact this benchmark has been caught"
        " by before; raising a cap that was registered as a cap is not that, and this file would be"
        " the wrong place to find out"
    ),
    "expected_request_cost": (
        "3 requests, counted by `plan` before any was spent and before any score of the new 36 was"
        " read. `plan` is free arithmetic - reach against the scorer's 1,048,576 bp input, and whether"
        " the finished all-chromosome sweep has already deleted an element inside the perturbed"
        " interval - and it reads no model answer, so costing the run in advance cannot feed back"
        " into what the run selects. 47 of the 50 graded intervals are already covered by the sweep;"
        " the three that are not are SEPHS2_Jurkat_chr16_30472k, TRAPPC2L_K562_chr16_88496k and"
        " VAPA_K562_chr18_9889k. The registered budget of 30 is unchanged and is not approached"
    ),
    "what_the_next_36_would_have_to_read": (
        "THE POINT OF REGISTERING THIS. The first 24 read 6/24 over every drawn locus and 6/20 where"
        " the model could answer (0.300), against a chance floor of 0.137. For the frame's combined"
        " rate to be called STABLE rather than a number that happened at n = 20, the 30 graded loci"
        " of the second half have to land in the same place: a second-half derived rate inside about"
        " 0.15 and 0.45 where the model could answer, which is the interval the first half's 6 of 20"
        " would be an unremarkable draw from, and a combined rate that therefore stays well below the"
        " 0.867 to 0.889 the three curated frames read and above the chance floor. If that is what"
        " comes back, the negative of 188f4ce is the frame's result at n = 50 and not an n = 20"
        " accident, and the deletion layer's 1 of 20 is the number to watch: it should stay near 0.05"
        " and it should keep naming the nearest coding TSS instead at something near 14 of 20"
    ),
    "what_a_divergence_would_mean": (
        "the two halves are NOT interchangeable and the registration says so before the numbers"
        " exist. The cap was applied in genome order, so the first 24 are chr1 to chr10 and the next"
        " 36 are chr11 to chrX. Genome order is not random with respect to gene density: the later"
        " chromosomes include the most gene-dense stretches in the genome (chr17, chr19, chr22),"
        " where the nearest coding TSS is closer, more coding genes sit inside the scorer's window,"
        " and the chance floor per locus is LOWER. So a second half that reads BELOW the first is the"
        " expected direction of a density effect rather than a new finding, and must be reported"
        " against its own chance floor, which this run computes per frame half rather than assuming."
        " A second half that reads materially ABOVE the first would be the interesting divergence,"
        " because it would run against that gradient, and it would mean the first half's 0.300 is a"
        " property of chr1-chr10 and not of the geometry the rule selects. Either way the two halves"
        " are reported separately, with their floors, and neither is dropped in favour of the other."
        " The combined rate over all 50 is the frame's rate, because all 60 were drawn by one rule"
        " and no locus was kept or dropped on anything anybody saw"
    ),
    "what_is_not_claimed": (
        "the two halves are not independent replications. They are one rule applied to one file, split"
        " by an arbitrary cap, so a difference between them is a statement about the genome's"
        " arrangement and about n, not about two experiments agreeing"
    ),
}

# ------------------------------------------------------------------ the registration, before the draw
PREREGISTRATION: dict[str, Any] = {
    "written": (
        "2026-09-21, before data/knowledge/crispri/EPCrisprBenchmark_combined_data."
        "heldout_5_cell_types.GRCh38.tsv.gz was read, before any locus in this frame was known, and"
        " before any request was spent. The rule below is `assess` and `select`; the frame is whatever"
        " they return"
    ),
    "why_a_fourth_frame": (
        "docs/LOCI-BENCHMARK.md section 21: across three frames the derived target rate moves 0.007"
        " (15/17, 8/9, 8/9) and the nearest-gene rule moves 0.445 (0.647, 0.333, 0.778). Section 21"
        " named the cause - at the third set 7 of 9 published targets are the element's own nearest"
        " coding TSS - and so raised the question no existing frame can answer: is the derived rate a"
        " property of the model or of the geometry the curators happened to pick. A fourth frame is"
        " worth requests only if its geometry is fixed by a rule written down before the loci are"
        " known, which is what this file is"
    ),
    "frame": (
        "the FOURTH sampling frame. Not merged into loci.PANEL, loci_candidates.CANDIDATES or"
        " loci_third.THIRD; its rates are reported beside 15/17, 8/9 and 8/9 with each frame named,"
        " and are never pooled with any of them"
    ),
    "the_rule": {
        "stated": "before the loci are picked. The frame is a function of the file and of GENCODE,"
        " and no locus is kept or dropped on anything the model says",
        "source": (
            "EPCrisprBenchmark_combined_data.heldout_5_cell_types.GRCh38.tsv.gz - the HELD-OUT arm of"
            " the ENCODE CRISPR benchmark (Gschwind et al. 2025, EngreitzLab/CRISPR_comparison), 190"
            " CRISPR positives over K562, GM12878, HCT116, Jurkat and WTC11. The K562 TRAINING arm is"
            " excluded because data/results/crispri_benchmark.json fitted logistic weights on it"
        ),
        "steps_in_order": [
            "1. `crispri.parse` keeps only rows with ValidConnection TRUE, which drops pairs that"
            " overlap a promoter or a target exon. Called, not reimplemented",
            "2. keep rows with Regulated TRUE: a significant NEGATIVE effect on the target when the"
            " element is silenced, which is the file's own definition of a positive",
            "3. group by perturbed element and cell line; an element's targets are ALL of its"
            " regulated genes in that cell",
            "4. every target must be protein coding in GENCODE. An element with a non-coding"
            " regulated target is dropped and counted apart, because the deletion layer ranks"
            " `predicted_coding` first, so a non-coding target tests the H19 defect rather than the"
            " geometry question, and mixing the two would confound both",
            "5. THE GEOMETRY RULE: no target of the element may be the nearest coding TSS to the"
            " element's midpoint (GENCODE canonical TSS, protein coding, whole chromosome). An"
            " element where the shortcut is right does not belong in this frame",
            f"6. the element midpoint must be at least {KEEP_OUT} bp from every locus in the first"
            " three frames, so the fourth frame is independent of them and the matched control draw"
            " stays clean",
            f"7. genome order (chr1..chr22, chrX, chrY, then start), capped at {MAX_LOCI}. The cap is"
            " applied to the order, never to the scores, so it cannot act as a quality filter. It was"
            f" {FIRST_DRAW} at the first run and is {MAX_LOCI} now, which is every element the rule"
            " passes: `CAP_RAISE` registers that change and was committed before the second run",
            "8. `loci.read_reach` last, as in every other frame: a target whose gene BODY lies"
            " outside the scorer's 1,048,576 bp input is unaskable, is never sent a request, and is"
            " reported as a reach fatality rather than as a miss",
        ],
        "first_matching_reason_is_recorded": "an element is rejected for the first step it fails in"
        " the order above, so the rejection counts are a partition and can be added up",
    },
    "why_this_is_not_loci_this_model_will_do_well_on": (
        "the rule selects the geometry on which the cheap answer is WRONG. At every locus in this"
        " frame the nearest coding TSS is a gene the publication says is not the target, so a model"
        " that has learned `name the nearest gene` scores zero here by construction. Nothing in the"
        " rule reads a model output, a deletion result, an effect size, a chromatin score or a"
        " distance threshold: it reads the file's own positive flag and GENCODE arithmetic. The frame"
        " is adversarial to the benchmark's own headline, not flattering to it, and the registered"
        " failure outcome below is the one this lane expects to have to write if the headline was"
        " geometry all along"
    ),
    "loci": (
        f"however many the rule returns, capped at {MAX_LOCI} in genome order (see `CAP_RAISE`; the"
        f" first run's cap was {FIRST_DRAW}). The count is not chosen:"
        " it is a property of the file. If the rule returns fewer than 5, the frame is reported as"
        " unassemblable at that n and no rate is quoted from it - that is a statement about how rare"
        " a long-range CRISPRi positive is, which is itself worth the section"
    ),
    "requests": {
        "budget": REQUEST_BUDGET,
        "how": "zero where the finished all-chromosome sweep has already deleted an element inside the"
        " perturbed interval; one stated-interval deletion for each element it has not. Counted before"
        " anything is spent and reported per locus",
        "expected": "low. data/results/crispri_benchmark.json records 3,849 of 4,378 held-out pairs"
        " already sitting on a deleted element, so most of this frame is expected to cost nothing. The"
        " budget exists for the rest",
        "over_budget": "loci are dropped from the END of genome order until the count fits, and the"
        " drop is reported. Never dropped by score, by cell type or by distance",
    },
    "predictions": {
        "target_derived": (
            "THE REGISTERED QUESTION. If the three frames' 0.88 is a property of the model, the rate"
            " here should stay well clear of the chance floor even though every locus is a long-range"
            " skip. If it is a property of the geometry the curators picked, it should fall towards"
            " the floor. No point prediction is offered because the two hypotheses are what is being"
            " separated, and a number written down now would only be a preference"
        ),
        "target_heuristic": (
            "LOW, and it is NOT zero by construction. The benchmark's heuristic layer is"
            " `loci.read_node` - the nearest coding TSS INSIDE the element's CTCF node - which is a"
            " different rule from the one that drew this frame. It can still be right where the node"
            " boundary excludes the nearer gene, so a non-zero reading here is a real measurement of"
            " how much the node adds to raw proximity, and is registered as the interesting part of"
            " the baseline line rather than as a beaten baseline"
        ),
        "nearest_coding_tss_anywhere": (
            "0 of n, BY CONSTRUCTION, because that is the rule that drew the frame. It is reported so"
            " the construction is visible and is never quoted as a result. Any section written from"
            " this run that compares the derived rate to this number rather than to the chance floor"
            " is misreading its own frame"
        ),
        "reach_fatalities": (
            "1 to 4, and registered as HIGHER than the third set's 0 of 9 for a stated reason: the"
            " geometry rule selects elements that skip over a nearer gene, and skipping correlates"
            " with distance, so this frame should meet the 524 kb reach limit more often than a"
            " hand-curated one. If it comes back 0 that is a surprise worth a line; if it comes back"
            " above 4 the frame is mostly unaskable and that is the result"
        ),
        "direction": (
            "NOT ASKABLE IN THIS FRAME, and registered as such before the draw. Every CRISPR positive"
            " in this file is an element whose silencing LOWERS the target, so every locus here is an"
            " activator and the repressor arm has n = 0. The activator arm is reported as a reading"
            " beside the other frames' and no claim about the deletion layer's sign is made from it."
            " Section 21's repressor question stays open at n = 2 and nothing here moves it"
        ),
        "cell": (
            "judgeable only at the K562 and GM12878 loci, which are the two CRISPRi lines with a"
            " counterpart among the eleven reader cells. HCT116, Jurkat and WTC11 have none and are"
            " registered in advance as unjudgeable on the cell axis rather than as misses"
        ),
        "value_on_syntax": (
            "expected not to separate, as at the candidates and at the third set. CRISPRi elements are"
            " a few hundred bp, so this frame is all short elements and section 20's length gradient"
            " says a set of short elements cannot carry the panel's 9/17. n is declared too small to"
            " decide it before the number exists"
        ),
        "naming_some_target": "near 90% at the matched controls, as everywhere else in this"
        " benchmark, so that claim will not separate and is not expected to",
        "non_coding_targets_dropped": (
            "counted, and the count is itself an answer to section 21's open item. Three frames have"
            " failed to find a published non-coding target clean enough to state; this rule will say"
            " how many of 190 held-out CRISPR positives have one, which is a fact about the literature"
            " that no amount of further curation was going to produce"
        ),
    },
    "what_each_outcome_would_mean": {
        "derived_rate_holds_near_the_other_three": (
            "the strongest statement this benchmark has made. The headline survives the geometry that"
            " kills the shortcut, so the 0.88 the three frames agree on is not the curators' choice of"
            " nearby elements. This would also mean the nearest-gene baseline's 0.445 spread is purely"
            " a property of the loci and the derived rate's 0.007 is not"
        ),
        "derived_rate_falls_to_the_chance_floor": (
            "THE FAILURE, and it is reported as plainly as the positive would be. It would mean the"
            " three frames' agreement was carried by their geometry: the benchmark would have been"
            " measuring proximity in a more expensive way, and every derived rate in sections 9, 19"
            " and 21 would have to be read with the element-to-target distance printed beside it. This"
            " lane writes that section if that is the number"
        ),
        "in_between": (
            "reported as the number with no verdict attached, with the per-layer split - deletion,"
            " eQTL, summed window - so a reader can see which layer survives the long range and which"
            " one was the proximity. No verdict is invented to make the section land"
        ),
        "too_few_loci": (
            "fewer than 5 after the rule: no rate is quoted, and the section reports how rare a"
            " long-range CRISPRi positive is in the held-out arm. That is a real finding about what"
            " the field has measured and about what any enhancer-gene benchmark can test"
        ),
    },
    "what_would_count_as_a_failure_of_the_RUN_rather_than_a_result": (
        "the positive control missing its target. The control is drawn from the same file by a stated"
        " rule - among the elements this frame REJECTED for having the shortcut right, in reach and"
        " clear of the earlier frames, the one with the smallest published element-to-TSS distance -"
        " so it is the easiest case the same data can offer. If the derived layers cannot name the"
        " target when the target is both the nearest gene and the closest one, the reading is broken"
        " and nothing else in the run is interpretable. The control is reported and is excluded from"
        " every rate"
    ),
    "hit_rules": {
        "source": "loci.score_target, loci.score_cell, loci.score_direction, unchanged and unwrapped",
        "target": "strict: a derived layer's FIRST-ranked gene is one of the published targets. The"
        " derived layers are the deletion, the GTEx eQTL and the summed-window readings; `among` is"
        " the lenient reading and is never counted in any rate",
        "provenance": "derived, heuristic and looked_up on separate lines, as in all three other sets",
    },
    "denominators": [
        "all loci the rule returns",
        "loci where the model could answer (loci.read_reach askable), the same second denominator the"
        " panel reports as target_derived_where_the_model_could_answer",
        "both quoted beside 15/17, 8/9 and 8/9. No fourth denominator is invented and no rate pools"
        " the four frames",
    ],
    "negative_controls": (
        "loci.pick_negatives, five matched windows per locus, the same four covariates - length"
        " exactly, GC within 0.04, distance to the nearest coding TSS within 35%, Zoonomia constrained"
        " fraction closest - drawn away from VISTA elements and GWAS hits and away from every locus in"
        " ALL FOUR frames with 200 kb of flank. `genomeos.compare.standardised` does the standardising"
        " and `input_presence` runs before any coverage stratum, because sections 8, 19 and 21 all"
        " found that GC, distance and constraint say nothing about whether anybody ever spent a"
        " measurement on the window"
    ),
}

#: the hunk `genomeos/benchmark/loci.py` needed, LANDED 2026-09-21 with the cap raise. It is kept as
#: the record of what the runtime patch was doing, and `register_stated_intervals` is kept because it
#: is idempotent and because a test asserts the two say the same thing.
NEEDED_LOCI_HUNK = {
    "landed": "2026-09-21, in the commit that raised the cap. loci.STATED_INTERVAL_RESULTS now lists"
    " this frame, so the deletion this frame paid for is read back from the file rather than from a"
    " global the run patched at import time",
    "file": "genomeos/benchmark/loci.py",
    "constant": "STATED_INTERVAL_RESULTS",
    "from": ("loci_stated_intervals", "loci_candidate_intervals", "loci_third_intervals"),
    "to": ("loci_stated_intervals", "loci_candidate_intervals", "loci_third_intervals", INTERVALS),
    "why": "loci._deletion_rows reads the panels' self-stated intervals from that tuple and labels"
    " them `stated_interval`. Without the fourth name, a deletion this frame paid for is written, is"
    " never read back, and the locus grades as though no request had been made",
}

ENSEMBL_JOIN: dict[str, Any] = {
    "landed": "2026-09-21, as the second of the two stages registered in `loci_reread.PREREGISTRATION`",
    "the_defect": (
        "step 4 of the rule asks whether a regulated target is protein coding, and it asked it of the"
        " benchmark file's `measuredGeneSymbol` column, which is as old as the screens that filled"
        " it. The same file carries `measuredGeneEnsemblId`, which does not go stale. Section 23"
        " joined it and found SSFA2 is ITPRID2 (10 elements), SARS is SARS1 and WDR61 is SKIC8, so"
        " twelve of the fourteen elements this frame rejected as non-coding were stale symbols and"
        " eight of them pass every step of this rule under the current symbol"
    ),
    "what_changed": (
        "`select` resolves every published symbol through the file's own Ensembl id and GENCODE's"
        " gene table before `assess` sees it. `assess`, the rule itself, is untouched - it is the"
        " names it is given that were wrong, not the steps it applies to them"
    ),
    "it_cuts_both_ways": (
        "a stale symbol can also have been DRAWN wrongly: if a drawn locus's target resolves to the"
        " element's nearest coding TSS, the rule's own step 5 rejects it and the frame loses a locus"
        " it had scored. The re-draw is reported as drawn-as-implemented against drawn-as-written,"
        " with both denominators, precisely because the correction is not one-directional"
    ),
    "why_it_is_a_flag_and_not_a_replacement": (
        "`join_ensembl_ids=False` reproduces the draw section 22 reported. A frame whose result can"
        " no longer be regenerated is a frame whose published numbers cannot be checked"
    ),
    "the_cap_binds_again_and_what_is_done_about_it": (
        "WRITTEN AFTER THE FREE DRAW WAS TAKEN AND BEFORE ANYTHING WAS SCORED, and it is disclosed"
        " that way because the membership of both draws was known when this was decided; no score"
        " was. The corrected rule passes 68 held-out elements rather than 60, so MAX_LOCI binds"
        " again, and applying it in genome order would ADD the 8 recovered elements on chr2 and"
        " chr15 and DROP 8 already-scored loci off the end - VAPA (two), MRPL4, ADGRE2, PDCD5,"
        " EIF3K, CEBPB and MSN, on chr18 to chrX. That draw would confound the correction with the"
        " loss of eight loci the correction says nothing about."
        " `CAP_RAISE` already fixed what the cap means: it set 60 because 60 was every element the"
        " rule then passed, 'so the cap no longer binds at all and cannot be re-cut later on"
        " anything that has been seen'. Applied to the rule as written, that principle gives 68, and"
        " 68 is a SUPERSET of the 60 - nothing section 22 scored is dropped and no locus is chosen -"
        " which is what makes it the reading that cannot be a cherry-pick. The corrected frame is"
        " therefore scored uncapped, and the cap-60 membership is reported in the section beside it"
        " so the alternative draw is visible rather than merely rejected"
    ),
}

CHROM_ORDER = {f"chr{c}": i for i, c in enumerate([*range(1, 23), "X", "Y"])}


def current_symbols(source: str = SOURCE) -> dict[str, str]:
    """Published symbol to the symbol GENCODE uses now, for the ones that have been renamed.

    The join is the benchmark file's own `measuredGeneEnsemblId` against GENCODE's gene table. No
    alias list is written here and no symbol is corrected from memory, which is what makes this a
    join rather than a curation. `loci_noncoding` is imported inside the function because it imports
    this module at the top of its own.
    """
    from genomeos.benchmark import loci_noncoding

    table = loci_noncoding.gene_table()
    out: dict[str, str] = {}
    for published, gid in loci_noncoding.ensembl_ids(source).items():
        now = (table.get(gid) or {}).get("symbol")
        if now and now != published:
            out[published] = now
    return out


def with_current_symbols(rows: list[dict[str, Any]], renamed: dict[str, str]) -> list[dict[str, Any]]:
    """Rewrite each element's regulated targets to the symbols GENCODE uses now, keeping the old ones."""
    out = []
    for r in rows:
        published = list(r["targets"])
        targets = sorted({renamed.get(t, t) for t in published})
        distances = {renamed.get(t, t): d for t, d in r["distances"].items()}
        out.append({**r, "targets": targets, "distances": distances, "published_targets": published})
    return out


def register_stated_intervals() -> None:
    """Make `loci._deletion_rows` read this frame's stated intervals too. Idempotent, append only."""
    if INTERVALS not in loci.STATED_INTERVAL_RESULTS:
        loci.STATED_INTERVAL_RESULTS = (*loci.STATED_INTERVAL_RESULTS, INTERVALS)


def earlier_frames() -> tuple[Expect, ...]:
    """Every locus of the first three frames, which this one may not repeat or sit beside."""
    return (*loci.PANEL, *loci_candidates.CANDIDATES, *loci_third.THIRD)


@contextlib.contextmanager
def keep_out_all_four_sets(panel: tuple[Expect, ...]):
    """Draw negative controls away from every locus in all four frames.

    `loci.candidate_windows` builds its keep-out list from the module-level `loci.PANEL`, so the only
    way to extend it without editing that file is to extend the global for the draw and restore it
    after, including on an exception. `loci_third.keep_out_all_three_sets` does the same for three.
    """
    original = loci.PANEL
    extra = tuple(e for e in (*earlier_frames(), *panel) if e not in original)
    loci.PANEL = (*original, *extra)
    try:
        yield loci.PANEL
    finally:
        loci.PANEL = original


def genome_key(row: dict[str, Any]) -> tuple[int, int, int, str]:
    """chr1..chr22, chrX, chrY, then start, then end, then cell: the frame's fixed order."""
    return (CHROM_ORDER.get(row["chrom"], 99), row["start"], row["end"], row.get("cell", ""))


def heldout(source: str = SOURCE, knowledge: Path = crispri.KNOWLEDGE) -> list[crispri.Pair]:
    """The held-out arm, through `crispri.parse`, which already drops invalid connections."""
    return crispri.load(source, knowledge)


def group_elements(pairs: list[crispri.Pair]) -> list[dict[str, Any]]:
    """One row per perturbed element and cell line, carrying every gene it regulates.

    Step 3 of the rule. The regulated flag is the file's own: a significant negative effect on the
    target when the element is silenced. A tested-but-not-regulated gene is kept only as a count.
    """
    by_element: dict[tuple[str, int, int, str], dict[str, Any]] = {}
    for p in pairs:
        key = (p.chrom, p.start, p.end, p.cell)
        row = by_element.setdefault(
            key,
            {
                "chrom": p.chrom,
                "start": p.start,
                "end": p.end,
                "cell": p.cell,
                "mid": p.midpoint,
                "targets": [],
                "distances": {},
                "datasets": set(),
                "tested_genes": 0,
            },
        )
        row["tested_genes"] += 1
        row["datasets"].add(p.dataset)
        if p.regulated:
            row["targets"].append(p.gene)
            row["distances"][p.gene] = int(p.distance)
    out = []
    for row in by_element.values():
        if not row["targets"]:
            continue
        row["targets"] = sorted(set(row["targets"]))
        row["datasets"] = sorted(row["datasets"])
        row["distance"] = min(row["distances"][g] for g in row["targets"])
        out.append(row)
    return sorted(out, key=genome_key)


def assess(
    row: dict[str, Any],
    nearest: tuple[str, int] | None,
    coding: set[str],
    taken: list[tuple[str, int]],
    keep_out: int = KEEP_OUT,
) -> dict[str, Any]:
    """Apply steps 4 to 6 of the rule to one element. The first failing step is the reason.

    Pure arithmetic over the file's own row, GENCODE symbols and the earlier frames' midpoints. It
    reads no model output, no effect size and no chromatin score, which is what makes the frame
    statable in advance.
    """
    verdict = {
        "locus": None,
        "chrom": row["chrom"],
        "element": [row["start"], row["end"]],
        "cell": row["cell"],
        "targets": list(row["targets"]),
        "distance": row["distance"],
        "nearest_coding": nearest[0] if nearest else None,
        "nearest_coding_distance": nearest[1] if nearest else None,
        "reason": None,
    }
    non_coding = [t for t in row["targets"] if t not in coding]
    if non_coding:
        verdict["reason"] = "a regulated target is not protein coding in GENCODE"
        verdict["non_coding_targets"] = non_coding
        return verdict
    if nearest is None:
        verdict["reason"] = "no coding TSS on this chromosome"
        return verdict
    if nearest[0] in row["targets"]:
        verdict["reason"] = "the nearest coding TSS IS a published target - the shortcut is right here"
        return verdict
    near = [c for c, m in taken if c == row["chrom"] and abs(m - row["mid"]) < keep_out]
    if near:
        verdict["reason"] = f"within {keep_out} bp of a locus in an earlier frame"
        return verdict
    verdict["locus"] = locus_name(row)
    return verdict


def locus_name(row: dict[str, Any]) -> str:
    """A stable, readable name: target, cell line and the element's position in kb."""
    return f"{row['targets'][0]}_{row['cell']}_{row['chrom']}_{row['start'] // 1000}k"


def as_expect(row: dict[str, Any], verdict: dict[str, Any], length: int) -> Expect:
    """One drawn locus as an `Expect`, with the trap the rule guarantees it has."""
    counterpart = CELL_COUNTERPART.get(row["cell"], {"cells": (), "gtex_tissues": (), "tissues": ()})
    mid = row["mid"]
    targets = tuple(row["targets"])
    # the symbol the screen published, kept in the citation whenever GENCODE has since renamed it:
    # the locus is graded against the current symbol and cited under the one the paper used.
    published = tuple(row.get("published_targets") or targets)
    renamed = (
        f". The screen published {', '.join(published)}; GENCODE now calls"
        f" {'/'.join(targets)} the same gene, joined on the file's own Ensembl id"
        if published != targets
        else ""
    )
    return Expect(
        locus=verdict["locus"],
        chrom=row["chrom"],
        element=(row["start"], row["end"]),
        window=(max(0, mid - WINDOW_HALF), min(length, mid + WINDOW_HALF)),
        classes=("program",),
        targets=targets,
        direction="activates",
        tissues=tuple(counterpart["tissues"]),
        gtex_tissues=tuple(counterpart["gtex_tissues"]),
        cells=tuple(counterpart["cells"]),
        distance=row["distance"],
        nearest_gene_trap=verdict["nearest_coding"],
        element_source="annotated",
        element_citation=(
            f"ENCODE CRISPR benchmark (Gschwind et al. 2025), held-out arm, {'/'.join(row['datasets'])}"
            f", {row['cell']}: silencing {row['chrom']}:{row['start']}-{row['end']} significantly lowers"
            f" {', '.join(targets)} ({row['tested_genes']} gene(s) tested at this element). The interval"
            " is the screen's own perturbed element, not one this panel drew" + renamed
        ),
        citations=(
            "Gschwind et al. 2025, the ENCODE CRISPR benchmark combined dataset"
            " (EngreitzLab/CRISPR_comparison), held-out arm over K562, GM12878, HCT116, Jurkat, WTC11",
            f"the original screen: {'/'.join(row['datasets'])}",
        ),
        answer_from=("CRISPRi", *row["datasets"]),
        note=(
            f"drawn by rule, not chosen: the nearest coding TSS to this element is"
            f" {verdict['nearest_coding']} at {verdict['nearest_coding_distance']} bp, and the published"
            f" target {targets[0]} is at {row['distance']} bp. A model that names the nearest gene is"
            " wrong here by construction, which is the whole reason the locus is in this frame"
        ),
    )


def select(
    results_dir: Path = RESULTS_DIR,
    source: str = SOURCE,
    cap: int = MAX_LOCI,
    progress=None,
    join_ensembl_ids: bool = True,
) -> dict[str, Any]:
    """Draw the frame. Free: the file, GENCODE and arithmetic. No request, no network.

    Returns the drawn panel, the control locus, and every rejection with its reason, so the frame can
    be checked against the rule by anyone who reads the result file.

    `join_ensembl_ids` resolves each published symbol through the file's own `measuredGeneEnsemblId`
    before the rule is applied, which is the rule as written; `False` reproduces the draw section 22
    reported, which is the rule as it was implemented. `ENSEMBL_JOIN` has the reason for both.
    """
    say = progress or (lambda _m: None)
    pairs = heldout(source)
    elements = group_elements(pairs)
    renamed = current_symbols(source) if join_ensembl_ids else {}
    if renamed:
        elements = with_current_symbols(elements, renamed)
    taken = [(e.chrom, (e.element[0] + e.element[1]) // 2) for e in earlier_frames()]
    verdicts: list[dict[str, Any]] = []
    drawn: list[Expect] = []
    rejected_shortcut: list[tuple[dict[str, Any], dict[str, Any], int]] = []
    for chrom in sorted({r["chrom"] for r in elements}, key=lambda c: CHROM_ORDER.get(c, 99)):
        mine = [r for r in elements if r["chrom"] == chrom]
        ch = loci.Chromosome(chrom, results_dir)
        try:
            coding = {sym for _t, sym in ch.coding}
            for row in mine:
                v = assess(row, ch.nearest_coding(row["mid"]), coding, taken)
                verdicts.append(v)
                if v["reason"] is None:
                    drawn.append(as_expect(row, v, ch.length))
                elif v["reason"].startswith("the nearest coding TSS IS"):
                    rejected_shortcut.append((row, v, ch.length))
        finally:
            ch.close()
        say(f"{chrom}: {len(mine)} perturbed element(s), {sum(1 for v in verdicts if not v['reason'])} kept")

    over_cap = [e.locus for e in drawn[cap:]]
    panel = tuple(drawn[:cap])
    control = None
    if rejected_shortcut:
        row, v, length = min(rejected_shortcut, key=lambda x: (x[0]["distance"], genome_key(x[0])))
        control = as_expect(row, {**v, "locus": f"CONTROL_{locus_name(row)}"}, length)
    reasons: dict[str, int] = {}
    for v in verdicts:
        reasons[v["reason"] or "kept"] = reasons.get(v["reason"] or "kept", 0) + 1
    return {
        "source": source,
        "perturbed_elements_with_a_regulated_gene": len(elements),
        "kept": len(panel),
        "over_cap": over_cap,
        "rejected_by_reason": reasons,
        "ensembl_join": {
            "applied": bool(renamed),
            "renamed_symbols": dict(sorted(renamed.items())) if renamed else {},
            "elements_with_a_renamed_target": sum(
                1 for r in elements if r.get("published_targets", r["targets"]) != r["targets"]
            ),
            "registration": ENSEMBL_JOIN,
        },
        "panel": panel,
        "control": control,
        "verdicts": verdicts,
        "rule": PREREGISTRATION["the_rule"],
    }


def plan(draw: dict[str, Any], results_dir: Path = RESULTS_DIR) -> list[dict[str, Any]]:
    """Which drawn loci the model could answer and what each would cost. No request, no network."""
    register_stated_intervals()
    panel = [*draw["panel"]] + ([draw["control"]] if draw["control"] else [])
    rows: list[dict[str, Any]] = []
    for chrom in sorted({e.chrom for e in panel}, key=lambda c: CHROM_ORDER.get(c, 99)):
        ch = loci.Chromosome(chrom, results_dir)
        try:
            for e in (x for x in panel if x.chrom == chrom):
                reach = loci.read_reach(ch, e)
                s, en = e.element
                scored = loci._deletion_rows(chrom, s, en, results_dir)  # noqa: SLF001
                rows.append(
                    {
                        "locus": e.locus,
                        "chrom": chrom,
                        "element": [s, en],
                        "length": en - s,
                        "targets": list(e.targets),
                        "cell": e.cells[0] if e.cells else (e.tissues[0] if e.tissues else None),
                        "trap": e.nearest_gene_trap,
                        "distance": e.distance,
                        "control": e.locus.startswith("CONTROL_"),
                        "askable": reach["askable"],
                        "targets_in_reach": reach["targets_in_reach"],
                        "targets_out_of_reach": reach["targets_out_of_reach"],
                        "graded": reach["askable"],
                        "ccres_over_element": len(ch.ccres_in(s, en)),
                        "already_scored_over_element": len(scored),
                        "requests": 0 if (not reach["askable"] or scored) else 1,
                    }
                )
        finally:
            ch.close()
    rows.sort(key=lambda r: (CHROM_ORDER.get(r["chrom"], 99), r["element"][0]))
    return within_budget(rows)


def within_budget(rows: list[dict[str, Any]], budget: int = REQUEST_BUDGET) -> list[dict[str, Any]]:
    """Drop loci from the END of genome order until the request count fits the registered budget."""
    spend = sum(r["requests"] for r in rows)
    for r in reversed(rows):
        if spend <= budget:
            break
        if r["requests"] and not r["control"]:
            r["graded"] = False
            r["dropped_for_budget"] = True
            r["requests"] = 0
            spend -= 1
    return rows


def graded(rows: list[dict[str, Any]], draw: dict[str, Any]) -> tuple[Expect, ...]:
    """The loci the run grades: every drawn one the reach filter and the budget keep, plus control."""
    keep = {r["locus"] for r in rows if r["graded"]}
    panel = [*draw["panel"]] + ([draw["control"]] if draw["control"] else [])
    return tuple(e for e in panel if e.locus in keep)


def reach_fatalities(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """How many drawn loci the free filter killed - a number in its own right, whatever it is."""
    real = [r for r in rows if not r["control"]]
    dead = [r for r in real if not r["askable"]]
    return {
        "loci": len(real),
        "died_at_the_reach_filter": len(dead),
        "registered_expectation": PREREGISTRATION["predictions"]["reach_fatalities"],
        "which": [
            {
                "locus": r["locus"],
                "target": r["targets_out_of_reach"][0]["target"] if r["targets_out_of_reach"] else None,
                "distance_kb": (r["targets_out_of_reach"][0].get("distance_bp", 0) // 1000)
                if r["targets_out_of_reach"]
                else None,
            }
            for r in dead
        ],
        "requests_saved": len(dead),
        "beside_the_other_frames": {
            "the_seventeen": "3 of 17 unaskable (SHH_ZRS, SOX9_PierreRobin, FTO_IRX3's IRX5)",
            "the_nine": "2 of 11 died, both exhibits placed on a bare offset",
            "the_third_set": "0 of 9, which was the registered prediction there",
            "the_fourth_frame": f"{len(dead)} of {len(real)}",
        },
        "reading": (
            "this frame was registered to lose MORE loci here than the third set did, because the"
            " geometry rule selects elements that skip a nearer gene and skipping correlates with"
            " distance. Whatever the number, it is a measurement of how much of the long-range"
            " regulatory literature a 1 Mb model can be asked about at all"
        ),
    }


def is_control(row: dict[str, Any]) -> bool:
    return str(row.get("locus", "")).startswith("CONTROL_")


def drawn_rows(result: dict[str, Any]) -> list[dict[str, Any]]:
    """The scored loci without the positive control, which is what every rate is computed over."""
    return [r for r in result["loci"] if not is_control(r)]


def without_the_control(result: dict[str, Any]) -> dict[str, Any]:
    """Re-aggregate over the drawn loci alone.

    `loci.build` aggregates over every row it was handed, and this frame hands it the positive
    control as well, because the control has to go through the identical readers to be worth
    anything. The registration says the control is excluded from every rate, so the rate the frame
    reports is this one and `aggregate_with_the_control` is kept beside it rather than thrown away.
    `loci.trim_locus` keeps each row whole, so this is the same function over a smaller list.
    """
    return loci.aggregate(drawn_rows(result))


def baselines(result: dict[str, Any]) -> dict[str, Any]:
    """The two nearest-gene rules kept apart: the one that drew the frame, and the node one.

    The first is 0 by construction and is printed so the construction is visible. The second is the
    benchmark's own heuristic layer and is a real reading at these loci.
    """
    anywhere_hits, node_hits, n = 0, 0, 0
    per_locus = []
    for r in drawn_rows(result):
        node = (r.get("readings") or {}).get("node") or {}
        want = set(r["expected"]["targets"])
        anywhere = node.get("nearest_coding_anywhere")
        in_node = node.get("target")
        n += 1
        anywhere_hits += int(anywhere in want)
        node_hits += int(in_node in want)
        per_locus.append(
            {
                "locus": r["locus"],
                "targets": sorted(want),
                "nearest_coding_anywhere": anywhere,
                "nearest_in_node": in_node,
                "node_named_nothing": in_node is None,
            }
        )
    return {
        "nearest_coding_anywhere": {
            "k": anywhere_hits,
            "n": n,
            "expected": "0 by construction: it is the rule that drew the frame",
            "reading": "printed so the construction is visible. Never quoted as a beaten baseline",
        },
        "nearest_tss_in_node": {
            "k": node_hits,
            "n": n,
            "rate": round(node_hits / n, 4) if n else None,
            "reading": "the benchmark's own heuristic layer, a different rule from the one that drew"
            " the frame, and a genuine reading here: it is right only where the CTCF node excludes the"
            " nearer coding gene. This is the number the section compares with 0.647, 0.333 and 0.778",
        },
        "silences": sum(1 for p in per_locus if p["node_named_nothing"]),
        "silences_reading": (
            "section 21 found the node rule's misses at long range are SILENCES rather than wrong"
            " answers, so a baseline that abstains at long range is not the baseline its rate"
            " describes. This frame is all long range, so the silence count is the direct measurement"
            " of that, and it is reported beside the rate rather than folded into it"
        ),
        "per_locus": per_locus,
    }


def what_the_layers_named(result: dict[str, Any]) -> dict[str, Any]:
    """Per derived layer: how often it named the published target, and what it named instead.

    This is the reading the frame exists to produce. Every locus here has a nearest coding TSS that
    the publication says is NOT the target, so `named the trap` and `named the target` are the two
    answers a layer can give, and the count between them says whether the layer is reading regulation
    or reading proximity. The three earlier frames could not ask it, because at most of their loci
    the two answers are the same gene.
    """
    rows = drawn_rows(result)
    layers: dict[str, dict[str, Any]] = {}
    for r in rows:
        for name, v in (
            ((r.get("score") or {}).get("scored") or {}).get("target", {}).get("by_layer", {}).items()
        ):
            layers.setdefault(name, {"provenance": v["provenance"], "hit": 0, "n": 0})
            layers[name]["n"] += 1
            layers[name]["hit"] += int(v["hit"])
    instead: dict[str, dict[str, int]] = {}
    for layer, key in (("deletion", "target"), ("node", "target")):
        tally = {"the published target": 0, "the nearest coding TSS": 0, "another gene": 0, "nothing": 0}
        for r in rows:
            named = ((r.get("readings") or {}).get(layer) or {}).get(key)
            trap = r["expected"]["nearest_gene_trap"]
            if named is None:
                tally["nothing"] += 1
            elif named in r["expected"]["targets"]:
                tally["the published target"] += 1
            elif named == trap:
                tally["the nearest coding TSS"] += 1
            else:
                tally["another gene"] += 1
        instead[layer] = tally
    return {
        "by_layer": layers,
        "what_it_named_instead": instead,
        "reading": (
            "the derived rate is a union over the deletion, eQTL and summed-window layers, so the"
            " per-layer line says which of them survives a locus where the target is not the nearest"
            " gene. `what_it_named_instead` is the sharper reading: a layer that names the nearest"
            " coding TSS at most of these loci is reporting proximity, whatever it scores at a frame"
            " where proximity and the published answer agree"
        ),
    }


def direction_readings(result: dict[str, Any]) -> dict[str, Any]:
    """The activator arm only, because the rule can draw nothing else. Registered before the draw."""
    rows = []
    for r in drawn_rows(result):
        d = ((r.get("score") or {}).get("scored") or {}).get("direction") or {}
        rows.append(
            {
                "locus": r["locus"],
                "published": r["expected"]["direction"],
                "judged": bool(d.get("judged")),
                "hit": bool(d.get("hit_derived")),
                "model_action": d.get("model_action"),
            }
        )
    judged = [r for r in rows if r["judged"]]
    return {
        "per_locus": rows,
        "activators": {"k": sum(1 for r in judged if r["hit"]), "n": len(judged)},
        "repressors": {"k": 0, "n": 0},
        "registered": PREREGISTRATION["predictions"]["direction"],
        "not_askable": (
            "every CRISPR positive in this file is an element whose silencing lowers its target, so"
            " the repressor arm is empty by construction and this frame cannot move the sign question"
            " section 21 left open at n = 2. The activator reading is reported and nothing is claimed"
            " about the deletion layer's sign from it"
        ),
    }


def positive_control(result: dict[str, Any]) -> dict[str, Any]:
    """The easiest case the same file can offer. A no here invalidates the rest of the run."""
    row = next((r for r in result["loci"] if is_control(r)), None)
    if row is None:
        return {"ran": False, "passed": None, "registered": PREREGISTRATION["the_rule"]}
    target = ((row.get("score") or {}).get("scored") or {}).get("target") or {}
    deletion = (row.get("readings") or {}).get("deletion") or {}
    return {
        "locus": row["locus"],
        "ran": True,
        "passed": bool(target.get("hit_derived")),
        "hit_by": target.get("hit_derived_by") or [],
        "targets": row["expected"]["targets"],
        "distance": row["expected"]["distance"],
        "deletion_named_first": deletion.get("target"),
        "deletion_log2_fold_change": (deletion.get("targets") or [{}])[0].get("best"),
        "excluded_from_every_rate": True,
        "registered": PREREGISTRATION["what_would_count_as_a_failure_of_the_RUN_rather_than_a_result"],
    }


def compare_with_controls(result: dict[str, Any]) -> dict[str, Any]:
    """The four claims at this frame against its matched windows, standardised.

    `genomeos.compare` does the standardising; this only shapes the rows, and `input_presence` runs
    before any coverage stratum for the reason sections 19 and 21 both record.
    """
    from genomeos.compare import Strata, imbalance, input_presence, standardised

    def row(gc, dist, con, claims) -> dict[str, Any] | None:
        if gc is None or dist is None:
            return None
        return {"gc": gc, "log_distance": (dist or 1) ** 0.5, "constrained": con or 0.0, **claims}

    def presence(con, claims) -> dict[str, Any]:
        return {
            "deletion_scored": True if claims.get("deletion_scored") else None,
            "eqtl_present": True if claims.get("eqtl_present") else None,
            "constrained": con,
        }

    claimed = [(r, loci.claims(r["readings"])) for r in drawn_rows(result)]
    targets = [
        row(r.get("gc"), r.get("distance_to_coding_tss"), r.get("constrained_fraction"), c)
        for r, c in claimed
    ]
    controls = [
        row(w.get("gc"), w.get("distance_to_coding_tss"), w.get("constrained_fraction"), w["claims"])
        for w in result["negatives"]
    ]
    t_presence = [presence(r.get("constrained_fraction"), c) for r, c in claimed]
    c_presence = [presence(w.get("constrained_fraction"), w["claims"]) for w in result["negatives"]]
    targets = [r for r in targets if r]
    controls = [r for r in controls if r]
    strata = Strata(gc=(0.40, 0.48, 0.56), log_distance=(100.0, 250.0, 500.0), constrained=(0.02, 0.10))
    covered = Strata(
        gc=(0.40, 0.48, 0.56),
        log_distance=(100.0, 250.0, 500.0),
        constrained=(0.02, 0.10),
        deletion_scored=(0.5,),
    )
    for r in targets + controls:
        r["deletion_scored"] = float(bool(r["deletion_scored"]))
    out: dict[str, Any] = {
        "reading": "direct standardisation on the covariates the windows were matched on; the"
        " `*_given_coverage` lines add whether an element inside the window has actually been deleted",
        "input_presence": input_presence(
            t_presence,
            c_presence,
            {
                "a deletion spent inside the window": "deletion_scored",
                "a GTEx eQTL distilled for the window": "eqtl_present",
                "the mammalian constraint track over the window": "constrained",
            },
        ),
        "imbalance": imbalance(targets, controls, strata),
    }
    for claim in ("target", "cell", "direction", "storage"):
        out[claim] = standardised(targets, controls, strata, hit=claim)
        out[f"{claim}_given_coverage"] = standardised(targets, controls, covered, hit=claim)
    return out


def beside_the_other_three(result: dict[str, Any], results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """This frame's headline rates beside 15/17, 8/9 and 8/9, all four frames named."""
    panel = (load_result("loci_benchmark", results_dir) or {}).get("aggregate") or {}
    cands = (load_result(loci_candidates.NAME, results_dir) or {}).get("aggregate") or {}
    third = (load_result(loci_third.NAME, results_dir) or {}).get("aggregate") or {}
    mine = result["aggregate"]

    def four(key: str) -> dict[str, Any]:
        return {
            "the_fourth_frame": {k: (mine.get(key) or {}).get(k) for k in ("k", "n", "rate")},
            "the_third_set": {k: (third.get(key) or {}).get(k) for k in ("k", "n", "rate")},
            "the_nine_candidates": {k: (cands.get(key) or {}).get(k) for k in ("k", "n", "rate")},
            "the_seventeen_panel": {k: (panel.get(key) or {}).get(k) for k in ("k", "n", "rate")},
        }

    return {
        "frames": {
            "the_seventeen_panel": "loci.PANEL, hand-curated famous loci, 2026-09-13",
            "the_nine_candidates": "loci_candidates.CANDIDATES, perturbation-sourced, 2026-09-17",
            "the_third_set": "loci_third.THIRD, perturbation-sourced, direction and length, 2026-09-21",
            "the_fourth_frame": "loci_fourth.select, DRAWN BY RULE from the ENCODE CRISPR benchmark's"
            " held-out arm: every element whose regulated target is not its own nearest coding TSS."
            " The first frame in this benchmark that nobody chose locus by locus, 2026-09-21",
        },
        "target_derived": four("target_derived"),
        "target_derived_where_the_model_could_answer": four("target_derived_where_the_model_could_answer"),
        "target_heuristic": four("target_heuristic"),
        "target_looked_up": four("target_looked_up"),
        "chance_floor": {
            "the_fourth_frame": mine.get("target_by_chance"),
            "the_third_set": third.get("target_by_chance"),
            "the_nine_candidates": cands.get("target_by_chance"),
            "the_seventeen_panel": panel.get("target_by_chance"),
        },
        "reading": (
            "four frames, four rates, never one. The first three were chosen locus by locus and agree"
            " to 0.007 while the nearest-gene rule over them moves 0.445. This one was drawn by a rule"
            " that fixes the geometry the first three left free, so its disagreement with them - in"
            " either direction - is the measurement, and its agreement with them would be the"
            " strongest statement the benchmark has made"
        ),
    }


def halves(result: dict[str, Any]) -> dict[str, Any]:
    """The frame split at the first cap: the 24 scored on 2026-09-21 and the 36 added after it.

    `CAP_RAISE` registered this split before the second half was scored, and it registered why the
    two halves are not interchangeable: the cap was applied in genome order, so the first half is
    chr1 to chr10 and the second is chr11 to chrX, and the later chromosomes carry the gene-dense
    stretches. Each half therefore gets its own chance floor rather than sharing the frame's.
    """
    order = [e["locus"] for e in result.get("drawn_panel", [])]
    first, second = set(order[:FIRST_DRAW]), set(order[FIRST_DRAW:])
    rows = drawn_rows(result)

    def half(name: str, want: set[str], drawn: int) -> dict[str, Any]:
        mine = [r for r in rows if r["locus"] in want]
        agg = loci.aggregate(mine) if mine else {}
        askable = agg.get("target_derived_where_the_model_could_answer", {})
        floor = agg.get("target_by_chance", {})
        return {
            "half": name,
            # the chromosome lives on the expectation, not on the row: `loci.trim_locus` keeps the
            # `expected` block whole and does not copy the coordinate up to the top level.
            "chromosomes": sorted(
                {r["expected"]["chrom"] for r in mine}, key=lambda c: CHROM_ORDER.get(c, 99)
            ),
            "drawn": drawn,
            "graded": len(mine),
            "died_at_the_reach_filter": drawn - len(mine),
            "target_derived_over_every_drawn_locus": {
                "k": agg.get("target_derived", {}).get("k"),
                "n": drawn,
                "rate": round(agg["target_derived"]["k"] / drawn, 3) if mine and drawn else None,
            },
            "target_derived_where_the_model_could_answer": {k: askable.get(k) for k in ("k", "n", "rate")},
            "target_heuristic": {k: agg.get("target_heuristic", {}).get(k) for k in ("k", "n", "rate")},
            "target_looked_up": {k: agg.get("target_looked_up", {}).get(k) for k in ("k", "n", "rate")},
            "chance_floor": {
                "expected": floor.get("expected"),
                "of": floor.get("of"),
                "rate": round(floor["expected"] / floor["of"], 3) if floor.get("of") else None,
            },
            "layers": what_the_layers_named({"loci": mine}),
            "loci": [r["locus"] for r in mine],
        }

    return {
        "why": CAP_RAISE["what_a_divergence_would_mean"],
        "first_24": half(
            f"the first {FIRST_DRAW}, scored 2026-09-21 at the registered cap", first, FIRST_DRAW
        ),
        "next_36": half(
            f"the {len(order) - FIRST_DRAW} the first cap left undrawn, same rule",
            second,
            len(order) - FIRST_DRAW,
        ),
        "registered_before_the_second_half_was_scored": CAP_RAISE["what_the_next_36_would_have_to_read"],
        "reading": (
            "the two halves are one rule applied to one file and split by an arbitrary cap, so a"
            " difference between them is a statement about the genome's arrangement and about n, not"
            " about two experiments agreeing. The frame's rate is the combined one over every locus"
            " the rule returned; these are printed so a reader can see whether it is carried by half"
            " of it"
        ),
    }


def readings(out: dict[str, Any], results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """Everything this frame computes from the scored rows. Free, and re-runnable on a saved result.

    Split out of `run` so a result can be re-aggregated without re-reading a single track: the rows
    are fixed once `loci.build` has written them, and every number below is a function of those rows.
    """
    out["aggregate_with_the_control"] = out["aggregate"]
    out["aggregate"] = without_the_control(out)
    out["aggregate_reading"] = (
        "the headline rates are over the DRAWN loci only. The positive control went through the same"
        " readers, so loci.build aggregated it in; the registration says it is excluded from every"
        " rate, and `aggregate_with_the_control` is kept beside this one rather than discarded"
    )
    out["layers"] = what_the_layers_named(out)
    out["halves"] = halves(out)
    out["cap_raise"] = CAP_RAISE
    out["baselines"] = baselines(out)
    out["directions"] = direction_readings(out)
    out["positive_control"] = positive_control(out)
    out["needed_loci_hunk"] = NEEDED_LOCI_HUNK
    out["against_controls"] = compare_with_controls(out)
    out["beside_the_other_three"] = beside_the_other_three(out, results_dir)
    return out


def reaggregate(results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """Recompute every derived block of a saved result. No request, no network, no re-reading."""
    out = load_result(NAME, results_dir)
    if not out:
        raise FileNotFoundError(f"no {NAME} result to re-aggregate")
    if "aggregate_with_the_control" in out:  # idempotent: always aggregate from build's own figure
        out["aggregate"] = out["aggregate_with_the_control"]
    readings(out, results_dir)
    save_result(NAME, out, results_dir)
    return out


def run(
    results_dir: Path = RESULTS_DIR,
    network: bool = True,
    negatives: bool = True,
    gtex_dir: Path = GTEX_DIR,
    progress=None,
    draw: dict[str, Any] | None = None,
    name: str = NAME,
) -> dict[str, Any]:
    """Draw the frame, filter by reach and budget, then read and score through `loci.build`.

    `name` is the result the run is saved under, and it exists so a re-drawn frame can be written
    beside the one that has been published rather than over it. The 2026-09-21 Ensembl join is
    scored as `loci_fourth_rejoined` for exactly that reason: section 22's numbers stay where they
    are and the corrected draw is reported next to them.
    """
    register_stated_intervals()
    draw = draw if draw is not None else select(results_dir, progress=progress)
    rows = plan(draw, results_dir)
    panel = graded(rows, draw)
    with keep_out_all_four_sets(panel):
        out = loci.build(
            panel=panel,
            results_dir=results_dir,
            gtex_dir=gtex_dir,
            network=network,
            negatives=negatives,
            progress=progress,
        )
    out["result"] = name
    out["preregistration"] = PREREGISTRATION
    out["draw"] = {k: v for k, v in draw.items() if k not in ("panel", "control")}
    out["drawn_panel"] = [e.as_dict() for e in draw["panel"]]
    out["control_locus"] = draw["control"].as_dict() if draw["control"] else None
    out["plan"] = rows
    out["reach"] = reach_fatalities(rows)
    # the rows are what the run PAID for - range reads over public tracks, and the requests before
    # them - and every block below is a pure function of them. They are saved before the readings
    # are computed so that a defect in a reading costs a re-aggregation (free, `--reaggregate`)
    # rather than the whole read. 2026-09-21: a KeyError in `halves` threw away 25 minutes of reads
    # that were already finished and correct.
    save_result(name, out, results_dir)
    readings(out, results_dir)
    out["note"] = (
        "A FOURTH, separately registered frame of published enhancer-gene loci, DRAWN BY A RULE"
        " written before the file was read rather than curated locus by locus, and graded by"
        " genomeos.benchmark.loci with no change to any scorer or hit rule. Reported beside the"
        " seventeen, the nine and the third set with each frame named, never pooled. " + loci.NOTE
    )
    save_result(name, out, results_dir)
    return out
