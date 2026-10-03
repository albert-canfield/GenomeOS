# SPDX-License-Identifier: AGPL-3.0-or-later
"""Does a compiled rule's argmax cell carry cell-type information?

`lane-clause1` established (registration `ba8c41f`, result `b4622ad`) that on its own population of
166 predicted rows the rule's cell is an argmax on every row - 150 by `max_drop`, 16 by `max_rise`,
0 assigned any other way. `direction_v2.ASSIGNED_CELL` is `predicted_coding['tissue']`, and
`predict.enhancer_target.predict_target` sets that field to `max_drop_tissue` or `max_rise_tissue`:
the tissue of the larger of the biggest fall and the biggest rise over every track the scorer
returned, chosen jointly with the gene. So a rule's "cell" is the place an extreme happened to fall
and not a context anyone chose before the values were read.

That leaves the question this module asks and does not answer by assertion: **is the argmax cell
nevertheless informative about where the element acts?** It is answerable because some element-gene
pairs were MEASURED in a named cell. For those pairs the rate at which the compiled cell IS the
measured cell can be set beside the rate at which that same cell label appears over every compiled
rule genome-wide - the base rate. A rate at or near the base rate means the cell label is not
evidence of where the rule acts; a rate well above it means the label carries something.

The asymmetry in `CONFOUND` is registered before any count and is the reason the two readings are
not symmetric: the CRISPRi benchmark chose which elements to test, and it tested elements in regions
active in the cell it was testing, so a positive difference is explained by that selection as
readily as by the model. A difference at or near zero is NOT explained by it, because the selection
pushes the other way. This module can therefore refute the proposition cleanly and can only fail to
refute it, never confirm it.

Nothing here edits a rule, a threshold or a `.bio` file: `direction_v2`, `compile` and
`enhancer_target` are imported and read.

READ DISCIPLINE. `enhancer_target.load_cached` falls back to the per-chromosome archive and caches
every archive it opens in a module-level dict that never evicts; a peer's inline loop over 93
elements reached 23 GB RSS on 2026-10-02 and drove the machine's free disk under the floor. No
function here opens an archive or calls that loader. The compact per-chromosome element table is
read by `stream_table`, which decodes ONE element at a time with `json.JSONDecoder.raw_decode` and
retains four scalars per element, so the largest chromosome's 59 MB table is never held as a parsed
object graph. `check_rss` raises at `RSS_CEILING_BYTES` after every chromosome rather than warning.
"""

from __future__ import annotations

import json
import math
import random
import resource
import sys
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from genomeos.attribution import compile as cp
from genomeos.attribution import crispri as cr
from genomeos.attribution import direction_v2 as dv2
from genomeos.attribution import targets as tg
from genomeos.predict import enhancer_target as et

# ---- the files this lane reads, and only these -------------------------------------------------

#: The committed genome-wide census of every compiled rule's `when: cell_type` label. The base rate
#: comes from here and from nowhere else, so it is a published figure this lane did not compute.
CENSUS = Path("data/results/context_evidence.json")
CENSUS_RELATIVE = "data/results/context_evidence.json"

#: The two cached ENCODE CRISPRi benchmark tables. `crispri.load` parses them; this lane reads a
#: pair's cell, locus and gene and reads no outcome column, no DHS value and no power column.
TRAINING = cr.TRAINING
HELDOUT = cr.HELDOUT

#: The compact per-chromosome element tables the all-enhancer sweep wrote, pointed at by each
#: committed `enhancer_targets_all_<chrom>` summary's `elements_where`. Plain JSON, never an archive.
COMPACT = cr.ELEMENTS

CHROMS = tuple(f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"])

# ---- the question, the population and the denominator ------------------------------------------

QUESTION = (
    "A compiled rule's cell is an argmax over the scorer's tracks: `predicted_coding['tissue']` is "
    "`max_drop_tissue` or `max_rise_tissue`, the place the largest movement happened to fall. Does "
    "that argmax carry cell-type information at all? For element-gene pairs MEASURED in a named "
    "cell, at what rate is the compiled cell label equal to the measured cell, against the rate at "
    "which that label appears over every compiled rule genome-wide?"
)

POPULATION = (
    "One arm per measured cell. An arm's population is every DISTINCT attributed element that "
    "overlaps at least one valid pair of the two cached ENCODE CRISPRi benchmark tables whose "
    "`CellType` is that cell. `attributed` is `genomeos.attribution.targets.attributed`'s own "
    "population - the elements of the three deletion runs that name a coding target gene, "
    "deduplicated by element id with the whole-chromosome run first - which is exactly the "
    "population `compile.compile_chromosome` emits one rule for, and therefore exactly the "
    "population the base rate is taken over. A pair is valid by the benchmark's own "
    "`ValidConnection` column, which `crispri.parse` applies; promoter and exon overlaps are "
    "already dropped there. Overlap is `OVERLAP_RULE`."
)

DENOMINATOR = (
    "DISTINCT ELEMENTS, not pairs. The compiled cell is a property of the element's rule and not of "
    "the pair, so counting an element once per tested gene would weight elements by how many genes "
    "the screen happened to assay near them. The pair-level figure is reported beside it as a "
    "secondary count and is never the headline. Both are stated with the arm's cell name."
)

DENOMINATOR_IS_SEPARATE = (
    "This lane's denominators are its own and are added to nothing. They are NOT to be added to or "
    "compared with the 112, 93, 73, 67, 26, 19, 18, 14 or 10 of tonight's other lanes: several of "
    "those overlap, several are different populations, and two were wrongly added together earlier "
    "today. Every figure below is reported with the arm and the denominator it was taken over."
)

WILL_REPORT = (
    "per arm: the measured cell; the number of valid pairs in that cell; how many of them overlap "
    "an attributed element; the number of DISTINCT attributed elements they overlap (the "
    "denominator); how many of those elements carry the arm's cell as their compiled label (the "
    "numerator); the observed rate; its Wilson 95% interval; the registered base rate; the "
    "difference; and which of the three registered readings the difference selects. Plus, over the "
    "same elements, the ten most frequent compiled labels with each one's own genome-wide base rate "
    "(the CROSS_ARM check), and the cross-arm matrix of every arm's rate on every other arm's cell."
)

# ---- the base rate, and how it is computed -----------------------------------------------------

BASE_RATE_RULE = (
    "base_rate(cell) = census['per_cell'][label]['rules'] / census['rules'], where census is the "
    "COMMITTED data/results/context_evidence.json and label is `compile.context(cell)` - the same "
    "function the compiler used to write the `when: cell_type` value the census counted. The census "
    "denominator is 440,589: every rule `compile.compile_chromosome` emits over the 24 chromosomes, "
    "counted as the compiler emits it, and its own `per_cell` entries sum to exactly that. No part "
    "of this base rate is computed by this lane and none of it is derived from the measured pairs, "
    "so it cannot have been chosen after seeing an observed figure."
)

BASE_RATE_IS_NOT_UNIFORM = (
    "The base rate is NOT 1/371. K562 is the single most frequent compiled label genome-wide, so "
    "its base rate is about twenty-three times what a uniform draw over the scorer's tracks would "
    "give. That is precisely why the base rate and not the track count is the control: a high rate "
    "of `label == K562` on K562-measured elements is what the genome-wide label distribution "
    "produces on its own, before any cell-type information is supposed."
)

SECOND_BASE_RATE = (
    "A second base rate is computed in the SAME streaming pass as the observed rate, from this "
    "lane's own reading of the element tables: the share of all attributed elements of all 24 "
    "chromosomes whose compiled label is the arm's cell. It exists to check the committed census "
    "rather than to replace it. The registered reading uses the COMMITTED base rate; the second is "
    "reported beside it with the difference, and a difference between the two is reported as a "
    "difference and not reconciled away."
)

#: The arms, in the order they are reported, named by the benchmark tables' own `CellType` string so
#: a pair is found by the name its own table gives it. The five cell types of the two cached tables:
#: K562 is the training table and one held-out arm, the other four are held out.
ARMS = ("K562", "GM12878", "HCT116", "Jurkat", "WTC11")

#: The compiled labels that name each arm's cell line, and the whole of the reason for the map.
#:
#: MEASURED, and it is a correction to a premise this lane held while writing the first draft of this
#: module: the two sides of the comparison do not spell every cell line the same way. The benchmark's
#: `CellType` for the Jurkat arm is `Jurkat`; the model's own track name for that line is recorded in
#: the committed census under TWO labels, `Jurkat__Clone_E6_1` (1,127 rules, CLO:0007045) and
#: `Jurkat` (5 rules, EFO:0002796), and `crispri.PREREGISTERED_PUBLISHED` already names CLO:0007045
#: as the term the model's output_metadata carries for Jurkat. A single string equality would
#: therefore have tested the arm against 5 of its 1,132 compiled rules and would have measured the
#: spelling rather than the question. Both labels are in the Jurkat set and the base rate is their
#: sum, so numerator and denominator share the set.
#:
#: Every other arm has exactly one label in the census. GM12891 and GM12892 are NOT in the GM12878
#: set: they are different individuals of the same trio, carry their own ontology terms
#: (EFO:0002785, EFO:0002786) and merging them would be merging three people.
LABELS_FOR: dict[str, tuple[str, ...]] = {
    "K562": ("K562",),
    "GM12878": ("GM12878",),
    "HCT116": ("HCT116",),
    "Jurkat": ("Jurkat", "Jurkat__Clone_E6_1"),
    "WTC11": ("WTC11",),
}

ALIAS_RULE = (
    "An arm's cell is matched against the SET of compiled labels that name that cell line, not "
    "against one string, and the same set fixes the base rate and the numerator so the two stay "
    "comparable. The set is LABELS_FOR, and it is one label for every arm but Jurkat. Jurkat has "
    "two, `Jurkat` (5 rules, EFO:0002796) and `Jurkat__Clone_E6_1` (1,127 rules, CLO:0007045), "
    "because the benchmark spells the line `Jurkat` while the model's own track name carries the "
    "clone; `crispri.PREREGISTERED_PUBLISHED` already names CLO:0007045 as the term the model's "
    "output_metadata gives Jurkat. Testing one spelling would have measured the spelling. "
    "GM12891 and GM12892 are deliberately NOT in the GM12878 set: they are other individuals of "
    "the same trio with their own ontology terms, and merging them would merge three people. "
    "A label the census does not carry contributes 0 rules and is reported as 0, never as absent."
)

# ---- the falsifier, stated before any count ----------------------------------------------------

TOLERANCE = 0.02  # absolute, on the difference between the observed rate and the committed base rate
USABLE = 0.25  # the share of an arm's elements that must carry the arm's cell for the label to be usable
MIN_ELEMENTS = 30  # fewer than this and the arm reports counts and NO rate reading

#: AMENDMENT 1's constants. The three above are UNCHANGED by it.
MIN_CLUSTERS = cr.MIN_CLUSTERS_FOR_AN_INTERVAL  # 10, the project's own floor, imported not restated
BOOTSTRAPS = cr.BOOTSTRAPS  # 2000, the project's own draw count
MIN_RESAMPLES = cr.MIN_RESAMPLES  # 1000, below which a percentile bound carries less than it says
SEED = 20261002  # fixed here, so the interval is the same number on every run of this committed code
LOCUS_SPAN = 1_000_000  # the span of a locus cluster when chromosomes do not give MIN_CLUSTERS
Z95 = 1.959964

AMENDMENT_1 = (
    "AMENDED on the coordinator's direction after the registration at 67e14d7 and BEFORE any count "
    "of this lane is committed. Additive: TOLERANCE 0.02, USABLE 0.25 and MIN_ELEMENTS 30 are "
    "unchanged, and no reading is made easier to reach.\n"
    "THE DEFECT IT FIXES, in the coordinator's own terms: the registered readings were selected by "
    "the POINT estimate of d against a +/-0.02 band. At 30 elements and a rate near 0.06 the "
    "standard error is about 0.043, more than twice the band, so on a small arm reading (2) - `NO "
    "cell-type information is detected, stated plainly` - would fire from noise. An underpowered "
    "null presented as a finding is the one outcome this lane must not produce.\n"
    "(a) A reading is now decided by the INTERVAL and never by the point. (3) detected only if the "
    "interval's LOWER bound on d exceeds +0.02; (1) anti-correlated only if its UPPER bound is "
    "below -0.02; (2) no information only if the WHOLE interval lies inside +/-0.02, which makes it "
    "an equivalence result and not a failure to reject.\n"
    "(b) A NEW reading (4) INCONCLUSIVE when the interval crosses a tolerance bound. Its words are "
    "`the data cannot tell`, never `no information`: those are different claims and eliding the "
    "difference is what produced several of tonight's overstatements.\n"
    "(c) The DECIDING interval is CLUSTERED - by chromosome, or by locus when chromosomes do not "
    "give at least 10 clusters. Under 10 clusters no interval is printed and the arm is inconclusive "
    "BY RULE and not by judgement. The Wilson interval over elements is kept beside it as SECONDARY "
    "and stays labelled as the binomial, unclustered interval it is.\n"
    "(d) POWER IS STATED PER ARM: the elements an arm needs for a +/-0.02 half-width, and whether it "
    "has them. An arm that does not is declared underpowered as a property of the arm.\n"
    "DISCLOSED, because the record must not read better than it was: when this amendment was "
    "written the un-amended run had ALREADY completed in a worktree at 67e14d7 and its figures had "
    "been seen. They are named in the amendment's commit message. The amendment is therefore NOT "
    "blind, and the one thing that can be checked rather than trusted is its direction: every "
    "change here can only turn a reading into INCONCLUSIVE or leave it standing. None of them can "
    "create a detection, widen a detection, or convert an inconclusive arm into a finding."
)

#: AMENDMENT 2's one new constant. A design effect is never allowed to NARROW an interval.
MIN_DESIGN_EFFECT = 1.0

AMENDMENT_2 = (
    "AMENDED a second time on the coordinator's direction, after the result at 7a44d59, and this "
    "amendment is POST-HOC and says so. It changes WHICH INSTRUMENT decides a reading on one shape "
    "of arm. No threshold moves: TOLERANCE 0.02, USABLE 0.25, MIN_ELEMENTS 30 and MIN_CLUSTERS 10 "
    "are unchanged.\n"
    "THE DEFECT. At 0 successes, or at n successes, a cluster bootstrap is DEGENERATE: every "
    "resample returns the same rate, so the interval has zero width and is estimating nothing. "
    "7a44d59 published WTC11's deciding interval as [0.0, 0.0] from 12 clusters and 2,000 "
    "resamples with `clustered_interval_decides_the_reading: true`, so reading (2) for that arm "
    "rested on a zero-width interval. The file did not say the interval was degenerate, and a "
    "reader could not have told it from a narrow one.\n"
    "THE RULE, stated generally rather than for the arm that exposed it: at 0 or n successes a "
    "cluster bootstrap decides nothing. The reading is then decided by a WILSON bound computed on "
    "an EFFECTIVE sample size, n_eff = n / design_effect, and the degenerate case is NAMED in the "
    "file (`bootstrap_degenerate`, `deciding_interval`) rather than silently passed through. An arm "
    "with 0 < k < n still routes to the CLUSTERED interval, unchanged. Both routings are planted by "
    "test, in both directions.\n"
    "THE DESIGN EFFECT IS MEASURED, NOT ASSUMED, because a degenerate arm cannot estimate its own "
    "intracluster correlation and an assumption would decide the reading. It is taken from the arms "
    "of this same result that CAN estimate one - 0 < k < n and a clustered interval - as "
    "(clustered width / Wilson width)^2, floored at MIN_DESIGN_EFFECT so it can never narrow an "
    "interval, and the LARGEST such value is applied: a 0-success arm offers no evidence that its "
    "own clustering is weaker than its siblings'. When no arm can estimate one the fallback is the "
    "most conservative there is, one observation per cluster (n_eff = clusters), which is "
    "intracluster correlation of 1. Every value, the arm it came from and the one taken are "
    "reported.\n"
    "DISCLOSED, and this is the part that is a correction rather than a confirmation. Before "
    "committing this rule I computed what three candidate choices imply for WTC11's 0 of 397. "
    "Independence (design effect 1) keeps the Wilson upper bound at 0.009583, below the tolerance, "
    "and reading (2) survives - which is what the coordinator's direction expected. The measured "
    "design effect from this result's own estimable arms is about 10.2 from K562, giving n_eff near "
    "39 and a Wilson upper bound near 0.067, so reading (2) does NOT survive and WTC11 becomes (4), "
    "THE DATA CANNOT TELL. Intracluster correlation of 1 gives n_eff = 12 and the same answer. The "
    "rule above was chosen because measuring the design effect is the right method and assuming "
    "independence is the thing clustering exists to doubt - not because of which answer it gives - "
    "and the outcome it gives is reported as a contradiction of the expectation that prompted it. "
    "It is also the conservative direction amendment 1 established: towards INCONCLUSIVE, never "
    "towards a finding.\n"
    "ALSO ADDED, and it moves nothing: `bootstrap_identical_share`, the share of resamples that "
    "returned the modal rate. A fully degenerate arm scores 1.0; HCT116's 1 of 248 is near-"
    "degenerate and will score high without triggering the rule, which is reported so a reader can "
    "see the near case rather than have this lane invent a second threshold for it."
)

FALSIFIER = (
    "Stated before any count was taken (67e14d7) and AMENDED before any count was committed "
    "(AMENDMENT_1, which only tightens). d = observed_rate - committed_base_rate, both over the "
    "arm's own denominator, and every branch below is decided by the DECIDING interval on d and "
    "never by the point estimate. The deciding interval is the CLUSTERED bootstrap, except on an "
    "arm at 0 or n successes, where that bootstrap is degenerate and AMENDMENT_2 routes the reading "
    "to a Wilson bound on an effective sample size instead.\n"
    "(1) upper bound of d < -0.02: the argmax lands on the measured cell LESS often than the "
    "label's own genome-wide frequency. The label is not evidence of where the rule acts, and on "
    "this arm it points away from it.\n"
    "(2) the WHOLE interval on d inside +/-0.02: NO cell-type information. This is an equivalence "
    "result and not a failure to reject: the arm has excluded any difference larger than the "
    "tolerance. On this arm the compiled cell is the label's background frequency and nothing more, "
    "so a rule's `cell` is NOT evidence of where it acts - and that is the finding, stated plainly.\n"
    "(3) lower bound of d > +0.02: a difference IS detected. Read with CONFOUND, which says in "
    "advance that this branch does not establish that the argmax carries cell-type information, "
    "because the benchmark's choice of which elements to test produces the same sign.\n"
    "(4) the interval crosses a tolerance bound: INCONCLUSIVE. THE DATA CANNOT TELL. This is not "
    "`no information` and may never be reported as one, and it is not a detection either.\n"
    "A SECOND, separate threshold, on usability rather than on detection: observed_rate >= 0.25 - "
    "the compiled cell names the measured cell on at least a quarter of the arm's elements, so a "
    "reader could use it; observed_rate < 0.25 - whatever the interval says, the label names the "
    "measured cell on a minority of elements and may not be read as the place the rule acts.\n"
    "An arm with fewer than 30 elements reports its counts and NO rate reading at all. An arm whose "
    "clusters number fewer than 10 prints NO interval and is INCONCLUSIVE BY RULE. The Wilson "
    "interval over elements is reported beside the clustered one as secondary, labelled binomial "
    "and unclustered wherever it appears, and decides nothing."
)

AMENDMENT_TIMELINE = (
    "WHAT WAS KNOWN WHEN, with the shas, because a claim about when something was known is worth "
    "nothing without them. `amendment_1` says it was made before any count of this lane was "
    "COMMITTED. That is true and it is misleading by omission, so the omission is repaired here.\n"
    "THE BLIND HALF, and it is the half the answer rests on. 67e14d7 (2026-10-02 22:20) fixed, "
    "before a single CRISPRi pair had been joined to a single element and before any element table "
    "had been opened: the population; the element denominator; the five base rates as numbers (K562 "
    "0.062292, GM12878 0.014540, Jurkat 0.002569, HCT116 0.001348, WTC11 0.000218) and the formula "
    "that produced them from the committed census; the usability threshold USABLE = 0.25; the "
    "tolerance 0.02; the element floor 30; and the SELECTION CONFOUND, which is the reason a "
    "positive difference may never be read as the argmax carrying cell-type information. The "
    "headline finding - that no arm reaches 0.25, so a rule's cell may not be read as the place it "
    "acts - rests on that threshold and on those base rates, and both were fixed at 67e14d7.\n"
    "THE PART THAT IS NOT BLIND. The first count RAN before amendment 1 was written. Its result was "
    "never committed and is not in the history, but the figures HAD BEEN SEEN when amendment 1 was "
    "written at 21b80b1 (22:30) and committed into the registration at 18224ef (22:31). So "
    "amendment 1 is NOT blind and must not be read as a pre-registration of the uncertainty rule. "
    "Amendment 2, at d3d0686, came later still - after the committed result at 7a44d59 (22:35) - "
    "and is post-hoc in its own first sentence.\n"
    "WHAT CAN BE CHECKED INSTEAD OF TRUSTED: the DIRECTION of both amendments. Every change in "
    "either can only turn a reading INCONCLUSIVE or leave it standing. None can create a detection, "
    "widen one, or turn an inconclusive arm into a finding. Amendment 1 removed two detections and "
    "no nulls; amendment 2 removed one null and no detections.\n"
    "AND WHAT THEY COST, named so the gain is not left implicit. Under the ORIGINAL point rule of "
    "67e14d7, GM12878 at 4 of 45 and Jurkat at 1 of 41 would both have read (3) DETECTION - two "
    "detections published on 45 and 41 elements. Under amendment 1 both are (4), THE DATA CANNOT "
    "TELL, inconclusive by rule at 9 and 5 resampling clusters. Under amendment 2 WTC11's 0 of 397 "
    "stops being an equivalence null decided by a zero-width interval."
)

REPORTING_CODE_CHANGED_AFTER_A_RUN = (
    "DISCLOSED because a reader who finds a code change between a registration and a result should "
    "not have to guess at it. 6c9a30e (2026-10-02 22:33) changed sentence-GENERATING code in "
    "scripts/argmaxcell_count.py after the amended registration at 18224ef and before the "
    "committed result at 7a44d59. The answer to the only question that matters about it is YES: it "
    "was written AFTER this lane had seen the OUTPUT of a run - the first run under amendment 1, "
    "whose printed verdict is where the fault was found.\n"
    "WHAT IT CHANGED. The verdict's split sentence said the reading-(2) arms were `at or below the "
    "base rate`. That was FALSE for HCT116, whose 0.004032 is above its base rate of 0.001348; the "
    "phrase was the point rule's wording left standing after amendment 1 made the reading an "
    "interval equivalence. It now names each arm by the reading its own interval selected. The "
    "same commit added, per arm, the tolerance and the observed rate as MULTIPLES of that arm's own "
    "base rate, because reading (2) is an equivalence against an ABSOLUTE 0.02 and that is 14.8 "
    "times HCT116's base rate and 91.7 times WTC11's.\n"
    "WHAT IT DID NOT CHANGE: no threshold, no branch rule, no population, no denominator, no base "
    "rate and no reading. It is a change to what the result SAYS about figures it had already "
    "computed, and `argmaxcell.reading` was not touched by it."
)

RELATIVE_BANDS_FROM_NOW_ON = (
    "ADOPTED on the supervisor's direction, and recorded here as a rule for the NEXT registration "
    "rather than as a change to this one: an equivalence or tolerance band is stated RELATIVE to "
    "its base rate, or it carries a stated reason why an absolute band is right. This lane's band "
    "is absolute (0.02) and carries no such reason, which is the defect; the patch kept for THIS "
    "result is the printed multiple of each arm's own base rate, so a reader can see that the band "
    "is 0.3 times K562's base rate and 91.7 times WTC11's and therefore which arms the equivalence "
    "reading is strong on. Nothing is regenerated for this beyond keeping those multiples."
)

CONFOUND = (
    "Registered before any count, and it makes the two readings asymmetric. The CRISPRi benchmark "
    "chose which elements to test, and it tested candidate elements in regions active in the cell "
    "it was testing. So the elements of the K562 arm are not a random draw from the compiled "
    "population: they are enriched for K562-active sequence before any model touched them. That "
    "enrichment pushes the observed rate UP. Therefore:\n"
    "- a difference at or near zero is NOT explained by the selection, because the selection pushes "
    "the other way, and it is clean evidence that the argmax carries no cell-type information;\n"
    "- a positive difference is explained by the selection at least as readily as by the model, and "
    "this lane cannot separate the two with the data it reads. A positive difference is therefore "
    "reported as a difference and MUST NOT be reported as the argmax carrying cell-type "
    "information.\n"
    "This lane can refute the proposition and can fail to refute it. It cannot confirm it."
)

CROSS_ARM = (
    "The one check that bears on the confound without new data, reported for every pair of arms: "
    "the rate of `label == C` over the elements of arm A, for every arm cell C. If the argmax "
    "tracked the cell a measurement was taken in, arm A's own cell would stand out among the "
    "columns of row A. If instead each column is near that column's genome-wide base rate whatever "
    "the row, the compiled label is tracking the genome-wide label distribution and not the "
    "measured cell. This is a comparison of measured rates and sets no threshold of its own."
)

# ---- what a name match is, and is not ----------------------------------------------------------

LABEL_RULE = (
    "An element's compiled label is `compile.context(predicted_coding['tissue'])`, which is "
    "`compile.ident` on a non-empty name and `unknown` on an absent one. The arm's cell is passed "
    "through the same function before comparison, so both sides of the equality are the string the "
    "compiler would write into `when: cell_type`. Two raw tissue names that differ only in "
    "characters `ident` replaces collapse to one label; that is a property of the census this base "
    "rate comes from, so the observed side must share it or the two would not be comparable."
)

NAME_MATCH_IS_NOT_A_CONFIRMATION = (
    "Equality of two labels is equality of two strings. It does not establish that the scorer's "
    "track and the screen's cells are the same material, the same passage or the same laboratory, "
    "and no file this lane reads could establish any of that."
)

METADATA_COPY_LIMITATION = (
    "Carried word for word from the correction a peer established, because it bounds what any cell "
    "statement of this lane may say: the saved copy of the model client's track metadata carries NO "
    "biosample accession and NO experiment accession. `one biosample` is therefore unsupportable "
    "from it; `one biosample NAME` is the most that can be said, and a two-track cell's two tracks "
    "are polyA plus RNA-seq and total RNA-seq of one biosample name, `not known to be biological "
    "replicates` in `_cell_summary`'s own words. This lane reads that metadata copy not at all and "
    "makes no claim that needs it."
)

CELL2_GROUP_CAVEAT = "an operational grouping, NOT established biological independence"

GROUPING = (
    "No grouping is applied. Elements are counted individually and the interval is binomial over "
    f"elements. If a later run groups them, `cell2.group`'s caveat travels with it: {CELL2_GROUP_CAVEAT}."
)

# ---- what this lane may not claim --------------------------------------------------------------

VALIDATES_NOTHING = (
    "This lane validates nothing about the direction rule in either direction. A peer's registered "
    "reading stands and is carried word for word: v2 can withhold a direction and never reverse "
    "one, and an unresolved class is NOT absence of regulation. Nothing measured here is evidence "
    "for or against v1's directions, and no reading below may be cited as such. The question here "
    "is about the CELL a rule names, not about the direction it calls."
)

NO_RECOMMENDATION = (
    "No rule, threshold, clause or `.bio` file is changed by this lane, and none is proposed for "
    "change. This is a measurement. If the finding is that a rule's cell is not evidence of where "
    "it acts, that is a finding for Albert and the supervisor; changing the cell assignment would "
    "invalidate every count already taken under it."
)

CANNOT_ESTABLISH = (
    "whether the scorer's track for a cell name and the screen's cells of that name are the same "
    "material: no file read here carries an accession on either side",
    "whether a compiled label that is NOT the measured cell is wrong: the element may act in both, "
    "and a screen in one cell measures nothing about the other",
    "whether the argmax carries cell-type information when the observed rate is ABOVE the base "
    "rate: the benchmark's selection of tested elements produces the same sign (CONFOUND)",
    "any direction, magnitude or effect: no outcome column of the benchmark tables is read, and no "
    "predicted value is read beyond the tissue field of the compact table",
    "anything about the 98% non-coding space outside the attributed elements: an element with no "
    "predicted coding target emits no rule and is in neither the base rate nor any arm",
)

COUNTS_NAMED = (
    "pairs_in_cell",
    "pairs_on_an_attributed_element",
    "elements",
    "elements_whose_label_is_the_measured_cell",
    "observed_rate",
    "deciding_interval",
    "deciding_ci95_on_the_rate",
    "deciding_ci95_on_the_difference",
    "bootstrap_degenerate",
    "bootstrap_identical_share",
    "effective_sample_size",
    "design_effect_applied",
    "clustered_ci95_on_the_rate",
    "clustered_ci95_on_the_difference",
    "clusters",
    "power",
    "reading_branch",
    "wilson95_secondary",
    "committed_base_rate",
    "difference_point_estimate",
    "second_base_rate_this_run",
    "label_distribution_top10",
    "cross_arm_matrix",
    "attributed_elements_streamed",
    "chromosomes_streamed",
    "max_rss_bytes",
)

REFUSALS = (
    "the committed census is absent or untracked: raise. A committed artefact's absence is a "
    "defect, not a difference, and a guard that skips on it would be measuring the machine",
    "the census `per_cell` entries do not sum to its own `rules`: raise, because the base rate's "
    "denominator would then not be the population it claims",
    "a compact element table named by a committed summary is not on this machine: raise. "
    "`targets.run_elements` already raises for this reason and the same rule holds here",
    "resident memory passes RSS_CEILING_BYTES after any chromosome: raise, naming the chromosome "
    "and the figure. Never widen the ceiling to make a run pass",
    "any call that would open a per-element archive or the cached-element loader: there is no such "
    "call, and the test suite asserts the module's source contains neither name",
    "an arm with fewer than 30 elements: counts are reported and NO rate reading is",
    "an arm with fewer than 10 resampling clusters: NO interval is printed and the arm is "
    "INCONCLUSIVE BY RULE, which is reading (4) and is never reported as `no information`",
    "an arm at 0 or n successes: the cluster bootstrap is DEGENERATE and decides nothing. The "
    "reading is decided by a Wilson bound on an effective sample size and the degeneracy is NAMED "
    "in the file; it is never passed through as though it were a narrow interval (AMENDMENT_2)",
)

# ---- the read discipline -----------------------------------------------------------------------

#: 2 GiB. Measured against a gate sample before the full run; the full run's own figure is reported.
RSS_CEILING_BYTES = 2 * 1024**3

READ_DISCIPLINE = (
    "No archive is opened and the cached-element loader is never called. The compact "
    "per-chromosome element table is decoded ONE element at a time by `stream_table` with "
    "`json.JSONDecoder.raw_decode`, and four scalars per element are retained, so the parsed object "
    "graph of a 59 MB table is never held. One chromosome is in hand at a time and its text is "
    "released before the next is opened. `check_rss` raises at RSS_CEILING_BYTES after every "
    f"chromosome. The ceiling is {RSS_CEILING_BYTES} bytes and the run reports its own peak."
)

#: The overlap rule, copied from `crispri.DeletionTable.overlapping` so the arm's elements are the
#: ones that module would have found. `test_argmaxcell` asserts the two agree element for element.
OVERLAP_RULE = (
    "`crispri.DeletionTable.overlapping`'s own rule, reimplemented over four scalars per element "
    "and asserted equal to it by test: an element overlaps [start, end) when element.end > start "
    f"and element.start > start - {cr.REACH}, searched from bisect_left(starts, end) - 1 downwards. "
    f"{cr.REACH} is `crispri.REACH`, the longest element in the table, and is imported rather than "
    "restated."
)

ARGMAX_RULE = getattr(dv2, "ASSIGNED_CELL", "")

ARGMAX_IS_A_SELECTION = (
    "`direction_v2.ASSIGNED_CELL` is `predicted_coding['tissue']`, and "
    "`enhancer_target.predict_target` sets that field to `max_drop_tissue` or `max_rise_tissue` - "
    "the tissue of the larger of the biggest fall and the biggest rise over every track, chosen "
    "jointly with the gene. An argmax is a selection: the cell is where an extreme fell, not a "
    "context chosen before the values were read. `lane-clause1` measured this on 166 of 166 rows of "
    "its own population (150 by max_drop, 16 by max_rise, 0 otherwise), registration ba8c41f, "
    "result b4622ad. That count is carried here, not recomputed, and this lane's denominators are "
    "separate from its 166."
)

EXPLORATION_PRECEDED_THIS = (
    "Stated plainly rather than implied. Before this registration was written this lane read: "
    "`genomeos/predict/enhancer_target.py` for CELLS, the cached-element loader and predict_target; "
    "`genomeos/attribution/direction_v2.py` for the assigned-cell field; "
    "`genomeos/attribution/targets.py` and `compile.py` for the attributed population and the label "
    "function; `genomeos/attribution/crispri.py` for the pair parser, the overlap rule and REACH; "
    "and the committed `data/results/context_evidence.json` for its `per_cell` block. From that "
    "last file it had SEEN the base rates before writing them down: K562 27,445/440,589 = 0.062292, "
    "HepG2 14,779 = 0.033544, GM12878 6,406 = 0.014540, IMR_90 2,942 = 0.006677, HCT116 594 = "
    "0.001348, WTC11 96 = 0.000218, Jurkat 5 and Jurkat__Clone_E6_1 1,127 = 0.002569 together, "
    "and it had checked for a second spelling of every arm's line in the census's own label list "
    "before fixing LABELS_FOR. It had seen NO observed rate, NO arm "
    "denominator and NO element of any arm: not one CRISPRi pair had been joined to an element and "
    "not one element table had been opened when this was written. What this registration binds is "
    "the population, the denominator, the base rate and its formula, the falsifier, the confound "
    "and the count list, all of them before the observed side exists."
)


# ---- the base rate ------------------------------------------------------------------------------


def census(path: Path = CENSUS) -> dict[str, Any]:
    """The committed genome-wide label census, with the refusals that make it a base rate applied.

    Absence RAISES. A committed artefact's absence is a defect, not a difference: a guard that
    skipped here would report a pass while measuring nothing, which is the defect class
    `must_be_committed()` was written against at 4f44dbf.
    """
    if not path.exists():
        raise FileNotFoundError(
            f"{path} is absent from this checkout. It is a committed result and the base rate of "
            "this lane is read from it; its absence is a defect and never a reason to skip"
        )
    with path.open() as fh:
        d = json.load(fh)
    per_cell = d.get("per_cell") or {}
    total = d.get("rules")
    if not per_cell or not isinstance(total, int) or total <= 0:
        raise ValueError(f"{path} carries no usable per_cell census (REFUSALS)")
    summed = sum(int(v["rules"]) for v in per_cell.values())
    if summed != total:
        raise ValueError(
            f"{path}: per_cell rules sum to {summed} but the census reports {total} rules, so the "
            "base rate's denominator is not the population it claims (REFUSALS)"
        )
    return d


def labels_for(cell: str) -> tuple[str, ...]:
    """`ALIAS_RULE`: the compiled labels that name one arm's cell line, the compiler's own spelling.

    An arm not in LABELS_FOR falls back to `compile.context(cell)` alone, so the function is total
    and a caller cannot get a silent empty set.
    """
    return LABELS_FOR.get(cell) or (cp.context(cell),)


def base_rate(cell: str, d: dict[str, Any]) -> dict[str, Any]:
    """`BASE_RATE_RULE` over `ALIAS_RULE`'s label set. The same set fixes the numerator."""
    per_cell = d.get("per_cell") or {}
    labels = labels_for(cell)
    rows = {label: (per_cell.get(label) or {}) for label in labels}
    rules = sum(int(r.get("rules") or 0) for r in rows.values())
    total = int(d["rules"])
    return {
        "cell": cell,
        "labels": list(labels),
        "rules_per_label": {label: int(r.get("rules") or 0) for label, r in rows.items()},
        "rules_with_this_label": rules,
        "rules_total": total,
        "rate": round(rules / total, 6),
        "how": BASE_RATE_RULE,
        "alias_rule": ALIAS_RULE,
        "source": CENSUS_RELATIVE,
        "ontology_terms": {label: r.get("ontology_term") for label, r in rows.items()},
    }


# ---- the streaming element reader ---------------------------------------------------------------

_DECODER = json.JSONDecoder()
_WS = " \t\r\n,"


def _results_dir() -> Path:
    from genomeos.results import RESULTS_DIR

    return RESULTS_DIR


def table_path(chrom: str, results_dir: Path | None = None) -> Path:
    """The compact element table a committed `enhancer_targets_all_<chrom>` summary points at.

    A summary whose table is not on this machine RAISES, which is `targets.run_elements`'s own rule.
    """
    from genomeos.results import load_result

    rd = results_dir or _results_dir()
    summary = load_result(f"enhancer_targets_all_{chrom}", rd) or {}
    where = summary.get("elements_where")
    if not where:
        raise FileNotFoundError(
            f"enhancer_targets_all_{chrom} names no elements_where, so the compact table cannot be "
            "located; this lane will not fall back to an archive (READ_DISCIPLINE)"
        )
    p = Path(where)
    if not p.is_absolute() and not p.exists():
        p = rd.parent.parent / where
    if not p.exists():
        raise FileNotFoundError(
            f"enhancer_targets_all_{chrom} points at {where}, which is not on this machine; rerun "
            "the job. This lane will not read the chromosome archive instead (READ_DISCIPLINE)"
        )
    return p


def stream_table(path: Path) -> Iterator[tuple[str, int, int, str]]:
    """(id, start, end, compiled label) for every element of one compact table that names a gene.

    One element is decoded at a time and four scalars are retained; the element's own dict is
    released before the next is decoded. No archive is opened and no cached-element loader is used.
    """
    text = path.read_text()
    i = text.find("[")
    if i < 0:
        raise ValueError(f"{path} is not a JSON list of elements")
    i += 1
    n = len(text)
    while True:
        while i < n and text[i] in _WS:
            i += 1
        if i >= n or text[i] == "]":
            return
        e, i = _DECODER.raw_decode(text, i)
        pc = e.get("predicted_coding") or {}
        if not pc.get("gene"):
            continue
        yield (str(e["id"]), int(e["start"]), int(e["end"]), cp.context(pc.get("tissue")))


def small_run_elements(chrom: str, results_dir: Path | None = None) -> Iterator[tuple[str, int, int, str]]:
    """The two sampled runs' attributed elements, which are inline in their committed summaries.

    `targets.RUNS` order is kept: the whole-chromosome run first (streamed by `stream_table`), then
    `constrained_targets`, then `enhancer_targets`. Deduplication by id is the caller's, exactly as
    `targets.attributed` does it, so the population here is that function's population.
    """
    for name in tg.RUNS[1:]:
        for e in tg.run_elements(name, chrom, results_dir or _results_dir()):
            pc = e.get("predicted_coding") or {}
            if not pc.get("gene"):
                continue
            yield (str(e["id"]), int(e["start"]), int(e["end"]), cp.context(pc.get("tissue")))


def attributed_table(chrom: str, results_dir: Path | None = None) -> list[tuple[int, int, str, str]]:
    """`targets.attributed(chrom)`'s population as (start, end, id, label), sorted by start.

    Four scalars per element and nothing else, so a chromosome of 45,000 elements costs a few
    megabytes instead of the table's 59.
    """
    seen: set[str] = set()
    out: list[tuple[int, int, str, str]] = []
    for eid, start, end, label in stream_table(table_path(chrom, results_dir)):
        if eid in seen:
            continue
        seen.add(eid)
        out.append((start, end, eid, label))
    for eid, start, end, label in small_run_elements(chrom, results_dir):
        if eid in seen:
            continue
        seen.add(eid)
        out.append((start, end, eid, label))
    out.sort(key=lambda t: t[0])
    return out


def overlapping(
    table: list[tuple[int, int, str, str]], starts: list[int], start: int, end: int
) -> list[tuple[int, int, str, str]]:
    """`OVERLAP_RULE`: `crispri.DeletionTable.overlapping`'s rule over the lean table."""
    import bisect

    out = []
    j = bisect.bisect_left(starts, end) - 1
    while j >= 0 and table[j][0] > start - cr.REACH:
        if table[j][1] > start:
            out.append(table[j])
        j -= 1
    return out


# ---- resident memory ----------------------------------------------------------------------------


def max_rss_bytes() -> int:
    """Peak resident set size of this PROCESS, in BYTES on every platform.

    macOS reports `ru_maxrss` in bytes and Linux in kilobytes; the unit is normalised here rather
    than left to the caller, because a ceiling compared against the wrong unit is a ceiling that
    never fires.

    THIS CHECK IS VALID ONLY IN A PROCESS THAT DOES NOTHING SUBSTANTIAL BEFORE THE WORK IT GUARDS -
    which is true of `scripts/argmaxcell_count.py`, a fresh interpreter that checks after each
    chromosome, and FALSE of a pytest session. `ru_maxrss` is a MONOTONIC HIGH-WATER: measured on
    this machine, 15.6 MB before allocating 300 MB, 330.2 MB after, and still 330.2 MB after
    freeing it. So it cannot attribute a byte to the work in front of it, and a test reading it
    live inherits every earlier test's peak - a full suite reached 4,406 MB, over this 2 GiB
    ceiling, and this file's test survived only by running early in the alphabet. A test passes
    `peak=` instead (`check_rss`); do not reuse this live reading inside one. A before-and-after
    DELTA is not the alternative: if an earlier test peaked higher the high-water does not move,
    the delta reads 0, and the guard becomes one that cannot fire.
    """
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(raw) if sys.platform == "darwin" else int(raw) * 1024


def check_rss(where: str, ceiling: int = RSS_CEILING_BYTES, peak: int | None = None) -> int:
    """The peak, or a RuntimeError naming it. Raises; never warns and never widens the ceiling.

    `peak` is the figure in bytes to judge. None means read this process live, which is what the
    count script does and the only context `max_rss_bytes` is valid in; a test injects the figure
    it means to test, because the live reading there is the whole session's high-water.
    """
    peak = max_rss_bytes() if peak is None else int(peak)
    if peak > ceiling:
        raise RuntimeError(
            f"resident memory reached {peak} bytes after {where}, over the registered ceiling of "
            f"{ceiling}. The run stops. Widening the ceiling to make it pass is not a fix "
            "(REFUSALS, READ_DISCIPLINE)"
        )
    return peak


# ---- the arms -----------------------------------------------------------------------------------


def wilson(k: int, n: int, z: float = 1.959964) -> list[float] | None:
    """A Wilson 95% interval on k/n, or None when n is 0. Binomial over elements, NOT clustered."""
    if n <= 0:
        return None
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(max(0.0, centre - half), 6), round(min(1.0, centre + half), 6)]


INTERVAL_IS_BINOMIAL = (
    "a Wilson 95% interval over ELEMENTS, assuming independent elements. It is NOT a clustered "
    "interval: elements of one locus are not independent, so it is narrower than a chromosome- or "
    "locus-clustered interval would be. It is reported because the readings turn on a difference "
    "from a base rate and a reader is entitled to the sampling width; it is labelled here and "
    "wherever it is printed, and no reading of this lane rests on it alone."
)


@dataclass
class Arm:
    """One measured cell's arm: its pairs, its elements and the counts registered in COUNTS_NAMED.

    `element_label` is the arm's whole population: one entry per DISTINCT attributed element that
    any of the arm's pairs overlaps, holding that element's compiled label. Every count below is
    derived from it, so the element denominator and the label tally cannot disagree.
    """

    cell: str
    labels: tuple[str, ...] = ()
    pairs: int = 0
    pairs_on_an_element: int = 0
    pair_level_matching: int = 0
    #: element id -> (chromosome, start, compiled label). The chromosome and start are carried for
    #: AMENDMENT_1's clustering and for nothing else; no coordinate of an element is reported.
    element_rows: dict[str, tuple[str, int, str]] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.labels = tuple(self.labels) or labels_for(self.cell)

    def add(self, chrom: str, elements: list[tuple[int, int, str, str]]) -> bool:
        """Record one pair's overlapping elements; True when any carries one of the arm's labels."""
        self.pairs += 1
        if not elements:
            return False
        self.pairs_on_an_element += 1
        hit = False
        for start, _end, eid, label in elements:
            self.element_rows[eid] = (chrom, start, label)
            hit = hit or label in self.labels
        if hit:
            self.pair_level_matching += 1
        return hit

    @property
    def elements(self) -> int:
        return len(self.element_rows)

    @property
    def elements_matching(self) -> int:
        return sum(1 for _c, _s, label in self.element_rows.values() if label in self.labels)

    def clusters(self, kind: str) -> dict[Any, list[bool]]:
        """The arm's elements grouped into resampling clusters, each as its match flags.

        `chromosome` is the first choice. `locus` is the fallback AMENDMENT_1 (c) names: a cluster is
        one chromosome's elements within a LOCUS_SPAN window of each other, keyed by the window. It
        IS a grouping, so `CELL2_GROUP_CAVEAT` travels with every figure taken from it.
        """
        out: dict[Any, list[bool]] = {}
        for chrom, start, label in self.element_rows.values():
            key = chrom if kind == "chromosome" else (chrom, start // LOCUS_SPAN)
            out.setdefault(key, []).append(label in self.labels)
        return out

    def cluster_choice(self) -> tuple[str, dict[Any, list[bool]]]:
        """`chromosome` when it reaches MIN_CLUSTERS, else `locus`; the chosen kind is reported."""
        by_chrom = self.clusters("chromosome")
        if len(by_chrom) >= MIN_CLUSTERS:
            return "chromosome", by_chrom
        return "locus", self.clusters("locus")

    @property
    def label_counts(self) -> dict[str, int]:
        """How many of this arm's elements carry each compiled label, the arm's whole distribution."""
        out: dict[str, int] = {}
        for _c, _s, label in self.element_rows.values():
            out[label] = out.get(label, 0) + 1
        return out

    @property
    def rate(self) -> float | None:
        return round(self.elements_matching / self.elements, 6) if self.elements else None


def cluster_interval(
    clusters: dict[Any, list[bool]], draws: int = BOOTSTRAPS, seed: int = SEED
) -> dict[str, Any]:
    """A percentile interval on the arm's rate, resampling whole CLUSTERS with replacement.

    `AMENDMENT_1` (c): this is the interval that DECIDES a reading. Under MIN_CLUSTERS clusters it
    returns no `ci95` at all and says why, because a percentile bound read from a handful of
    resampling units carries no more information than those units do - `crispri`'s own rule, whose
    constant is imported rather than restated.
    """
    keys = sorted(clusters, key=repr)
    made = {
        "clusters": len(keys),
        "cluster_minimum": MIN_CLUSTERS,
        "resamples_requested": draws,
        "resamples_minimum": MIN_RESAMPLES,
        "seed": seed,
        "method": "percentile bootstrap over whole clusters, resampled with replacement",
        "met_minimum": draws >= MIN_RESAMPLES,
    }
    if len(keys) < MIN_CLUSTERS:
        return {
            "ci95": None,
            "made": made,
            "why_no_interval": (
                f"{len(keys)} resampling clusters is below the floor of {MIN_CLUSTERS}, so no "
                "interval is printed and this arm is INCONCLUSIVE BY RULE rather than by judgement"
            ),
        }
    rng = random.Random(seed)
    rates = []
    for _ in range(draws):
        k = n = 0
        for _ in keys:
            flags = clusters[keys[rng.randrange(len(keys))]]
            n += len(flags)
            k += sum(flags)
        rates.append(k / n if n else 0.0)
    rates.sort()
    lo = rates[int(0.025 * (len(rates) - 1))]
    hi = rates[int(0.975 * (len(rates) - 1))]
    modal = max(set(rates), key=rates.count)
    made["identical_share"] = round(rates.count(modal) / len(rates), 4)
    made["degenerate"] = made["identical_share"] >= 1.0
    return {"ci95": [round(lo, 6), round(hi, 6)], "made": made, "why_no_interval": None}


def elements_needed(base: float, half_width: float = TOLERANCE, z: float = Z95) -> int:
    """`AMENDMENT_1` (d): the elements an arm needs for a +/-half_width interval at a given rate.

    z^2 p (1-p) / half_width^2, the normal-approximation sample size. It is reported at TWO rates
    and the reason is stated rather than left for a reader to notice: at the BASE rate - which is
    what the coordinator asked for, and which is the width of the null the arm is testing against -
    and at the arm's OBSERVED rate, which is the width the arm's own interval actually has. For a
    cell whose base rate is near zero the first number is small and the second is not, and quoting
    only the first would make an arm look powered that is not.
    """
    return max(1, math.ceil(z * z * base * (1 - base) / (half_width * half_width)))


def power(arm: Arm, base: float) -> dict[str, Any]:
    """Whether this arm can carry a reading at all, as a property of the arm."""
    observed = arm.rate if arm.rate is not None else 0.0
    at_base = elements_needed(base)
    at_observed = elements_needed(observed)
    return {
        "elements": arm.elements,
        "committed_base_rate": round(base, 6),
        "elements_needed_at_the_base_rate": at_base,
        "has_them_at_the_base_rate": arm.elements >= at_base,
        "observed_rate": arm.rate,
        "elements_needed_at_the_observed_rate": at_observed,
        "has_them_at_the_observed_rate": arm.elements >= at_observed,
        "underpowered": arm.elements < max(at_base, at_observed),
        "how": elements_needed.__doc__,
    }


def bootstrap_degenerate(k: int, n: int) -> bool:
    """`AMENDMENT_2`: at 0 or n successes every cluster resample returns the same rate."""
    return n > 0 and (k == 0 or k == n)


def design_effect(clustered: list[float] | None, wil: list[float] | None) -> float | None:
    """(clustered width / Wilson width)^2, floored at MIN_DESIGN_EFFECT. None when either is absent.

    Floored because a design effect below one would NARROW an interval, and clustering is never a
    reason to claim more precision than the independent case.
    """
    if not clustered or not wil:
        return None
    w = wil[1] - wil[0]
    if w <= 0:
        return None
    return max(MIN_DESIGN_EFFECT, ((clustered[1] - clustered[0]) / w) ** 2)


def measured_design_effect(arms: dict[str, Arm], draws: int = BOOTSTRAPS, seed: int = SEED) -> dict[str, Any]:
    """`AMENDMENT_2`'s design effect: the largest one estimable from this result's own arms.

    An arm can estimate one when 0 < k < n and its clusters reach MIN_CLUSTERS. Every value is
    reported with the arm it came from, so a reader can see what the number rests on.
    """
    per_arm: dict[str, float | None] = {}
    for name, arm in arms.items():
        k, n = arm.elements_matching, arm.elements
        if n < MIN_ELEMENTS or bootstrap_degenerate(k, n):
            per_arm[name] = None
            continue
        ci = cluster_interval(arm.cluster_choice()[1], draws=draws, seed=seed)
        per_arm[name] = design_effect(ci["ci95"], wilson(k, n))
    usable = {a: v for a, v in per_arm.items() if v is not None}
    taken = max(usable.values()) if usable else None
    return {
        "per_arm": {a: (round(v, 3) if v is not None else None) for a, v in per_arm.items()},
        "estimable_arms": sorted(usable),
        "taken": None if taken is None else round(taken, 3),
        "taken_from": (None if taken is None else sorted(a for a, v in usable.items() if v == taken)[0]),
        "rule": (
            "(clustered width / Wilson width)^2 per estimable arm, floored at "
            f"{MIN_DESIGN_EFFECT}, largest taken. An arm can estimate one when 0 < k < n and its "
            f"clusters reach {MIN_CLUSTERS}. With none estimable the fallback is one observation "
            "per cluster, which is intracluster correlation of 1 (AMENDMENT_2)"
        ),
    }


def reading(
    arm: Arm,
    base: float,
    draws: int = BOOTSTRAPS,
    seed: int = SEED,
    deff: float | None = None,
) -> dict[str, Any]:
    """`FALSIFIER` under `AMENDMENT_1`, applied. The branch is chosen by the CLUSTERED interval."""
    n, k = arm.elements, arm.elements_matching
    if n < MIN_ELEMENTS:
        return {
            "rate_reported": False,
            "reading": "counts only",
            "why": (
                f"{n} elements is below the registered floor of {MIN_ELEMENTS}, so this arm reports "
                "its counts and no rate reading. The counts stand; the reading does not exist"
            ),
            "counts_only": {"elements": n, "elements_whose_label_is_the_measured_cell": k},
            "power": power(arm, base),
        }
    rate = k / n
    d = rate - base
    kind, clusters = arm.cluster_choice()
    ci = cluster_interval(clusters, draws=draws, seed=seed)
    band = ci["ci95"]
    wil = wilson(k, n)
    # AMENDMENT_2: a degenerate bootstrap decides nothing, so the Wilson bound on an EFFECTIVE
    # sample size decides instead, and the degeneracy is named rather than passed through.
    degenerate = bootstrap_degenerate(k, n)
    deciding = "clustered_bootstrap"
    n_eff = k_eff = None
    used_deff = None
    if degenerate:
        deciding = "wilson_on_an_effective_sample_size"
        used_deff = deff if deff is not None else (n / len(clusters) if clusters else 1.0)
        used_deff = max(MIN_DESIGN_EFFECT, float(used_deff))
        n_eff = max(1, round(n / used_deff))
        k_eff = min(n_eff, round(k * n_eff / n))
        band = wilson(k_eff, n_eff)
        ci = {
            **ci,
            "why_no_interval": (
                f"the cluster bootstrap is DEGENERATE at {k} of {n} successes: every resample "
                "returns the same rate, so it has zero width and estimates nothing. It decides "
                f"nothing here. The reading is decided by a Wilson bound on an effective sample "
                f"size of {n_eff} = {n} / design effect {used_deff:.3f} (AMENDMENT_2)"
            ),
        }
    d_lo = None if band is None else round(band[0] - base, 6)
    d_hi = None if band is None else round(band[1] - base, 6)
    if band is None:
        branch, carries = 4, None
        detected = (
            "INCONCLUSIVE BY RULE: THE DATA CANNOT TELL. "
            + str(ci["why_no_interval"])
            + ". This is NOT `no cell-type information` and may never be reported as one; it is also "
            "not a detection. The point estimate and the counts are reported and decide nothing"
        )
    elif d_lo is not None and d_lo > TOLERANCE:
        branch, carries = 3, None
        detected = (
            "a difference IS detected: the deciding interval's lower bound on d is above the "
            "registered tolerance. It does NOT establish that the argmax carries cell-type "
            "information: the benchmark chose elements active in the cell it was testing, which "
            "produces this sign on its own, and this lane cannot separate the two (CONFOUND)"
        )
    elif d_hi is not None and d_hi < -TOLERANCE:
        branch, carries = 1, False
        detected = (
            "the argmax lands on the measured cell LESS often than that label's own genome-wide "
            "frequency, and the deciding interval's upper bound is below the tolerance. On this "
            "arm the compiled cell is not evidence of where the rule acts, and it points away from it"
        )
    elif d_lo is not None and d_hi is not None and d_lo >= -TOLERANCE and d_hi <= TOLERANCE:
        branch, carries = 2, False
        detected = (
            "NO cell-type information, as an EQUIVALENCE result and not a failure to reject: the "
            "whole deciding interval on d lies inside the registered tolerance, so this arm has "
            "excluded any difference larger than it. The rate at which the compiled cell is the "
            "measured cell is that label's background frequency and nothing more, so a rule's "
            "`cell` is NOT evidence of where it acts. The confound pushes the other way, so this "
            "reading is not explained by the benchmark's selection of tested elements"
        )
    else:
        branch, carries = 4, None
        detected = (
            "INCONCLUSIVE: THE DATA CANNOT TELL. The deciding interval on d crosses a tolerance "
            "bound, so this arm neither detects a difference larger than the tolerance nor excludes "
            "one. This is NOT `no cell-type information` and may never be reported as one"
        )
    return {
        "rate_reported": True,
        "reading": f"({branch})",
        "observed_rate": round(rate, 6),
        "committed_base_rate": round(base, 6),
        "difference_point_estimate": round(d, 6),
        "point_estimate_decides_nothing": (
            "carried because a reader will want it; under AMENDMENT_1 the branch is chosen by the "
            "clustered interval alone"
        ),
        "tolerance": TOLERANCE,
        "deciding_interval": deciding,
        "deciding_ci95_on_the_rate": band,
        "deciding_ci95_on_the_difference": None if band is None else [d_lo, d_hi],
        "bootstrap_degenerate": degenerate,
        "bootstrap_identical_share": ci["made"].get("identical_share"),
        "effective_sample_size": n_eff,
        "effective_successes": k_eff,
        "design_effect_applied": None if used_deff is None else round(used_deff, 3),
        "clustered_ci95_on_the_rate": None if degenerate else band,
        "clustered_ci95_on_the_difference": (None if (degenerate or band is None) else [d_lo, d_hi]),
        "clustered_by": kind,
        "cluster_grouping_caveat": (CELL2_GROUP_CAVEAT if kind == "locus" else None),
        "clustered_interval_provenance": ci["made"],
        "why_no_interval": ci["why_no_interval"],
        "clustered_interval_decides_the_reading": not degenerate,
        "amendment_2": AMENDMENT_2,
        "wilson95_secondary": wil,
        "wilson95_is": INTERVAL_IS_BINOMIAL,
        "wilson95_excludes_the_base_rate": None if wil is None else not (wil[0] <= base <= wil[1]),
        "detection": detected,
        "argmax_carries_cell_type_information": carries,
        "power": power(arm, base),
        "usable": rate >= USABLE,
        "usability": (
            f"observed rate {rate:.4f} is at or above the registered {USABLE}: the compiled cell "
            "names the measured cell on at least a quarter of this arm's elements"
            if rate >= USABLE
            else f"observed rate {rate:.4f} is below the registered {USABLE}: whatever the interval "
            "says, the compiled cell names the measured cell on a MINORITY of this arm's elements "
            "and may not be read as the place the rule acts"
        ),
    }


def top_labels(arm: Arm, d: dict[str, Any], n: int = 10) -> list[dict[str, Any]]:
    """The arm's most frequent compiled labels, each beside its own genome-wide base rate."""
    total = arm.elements or 1
    out = []
    for label, count in sorted(arm.label_counts.items(), key=lambda kv: (-kv[1], kv[0]))[:n]:
        row = (d.get("per_cell") or {}).get(label) or {}
        bt = int(row.get("rules") or 0) / int(d["rules"])
        out.append(
            {
                "label": label,
                "elements": count,
                "share_of_arm": round(count / total, 6),
                "genome_wide_base_rate": round(bt, 6),
                "ratio_to_base_rate": round((count / total) / bt, 3) if bt else None,
            }
        )
    return out


#: Names that would mean an archive had been opened: the cached-element loader, the archive helper
#: it calls, and the decompression module either would need. `names_used` asserts their absence from
#: this module's own CODE rather than from its text, so the sentence in the docstring above that has
#: to name the loader in order to say it is not used does not make the check pass or fail.
ARCHIVE_NAMES = ("load" + "_cached", "_archive", "gzip", "GzipFile")

_SOURCE = Path(__file__)


def names_used(source: str | None = None) -> set[str]:
    """Every identifier this module's code refers to: imports, names and attribute names.

    Read with `ast`, so a string literal and a docstring are not identifiers. A text scan would have
    to decide whether the docstring sentence that names the loader in order to say it is never
    called counts as a use, and whichever way it decided it would be measuring prose.
    """
    import ast

    tree = ast.parse(source if source is not None else _SOURCE.read_text())
    out: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute):
            out.add(node.attr)
        elif isinstance(node, ast.Name):
            out.add(node.id)
        elif isinstance(node, ast.Import):
            for a in node.names:
                out.add(a.name.split(".")[0])
                out.add((a.asname or a.name).split(".")[0])
        elif isinstance(node, ast.ImportFrom):
            out.add((node.module or "").split(".")[0])
            for a in node.names:
                out.add(a.name)
                out.add(a.asname or a.name)
    return out


def source_opens_no_archive() -> bool:
    """True when this module's code refers to no archive loader. Asserted by test, not claimed."""
    return not (names_used() & set(ARCHIVE_NAMES))


# Named so a reader can see the module knows what it is not using. et is imported for CELLS only.
RETAINED_CELLS_TODAY = tuple(et.CELLS)
RETENTION_IS_NOT_THIS_QUESTION = (
    "`enhancer_target.CELLS` is the four cell lines whose own track value is kept per gene. It "
    "bounds which cells a `by_cell` value exists for and so which rows clause (1) can reach - a "
    "peer's registered reading. It does NOT bound this question: the argmax tissue is over every "
    "track the scorer returned, so a compiled label may be any of the census's 317 and HCT116's 594 "
    "rules are assessable here although HCT116 carries no `by_cell` column at all."
)
