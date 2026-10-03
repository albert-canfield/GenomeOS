# SPDX-License-Identifier: AGPL-3.0-or-later
"""Item (h): the astrocyte-argmax question, REGISTERED and nothing else.

THE QUESTION. On the 1,232 paid astrocyte-element answers of AstroREG-2, at what rate is the
argmax label -- the track carrying the largest movement in a gene row's whole track vector -- a
brain/astrocyte track, against (a) that label set's base rate over the compiled rules and (b) the
K562 arm's rate for the same label set?

NO ANSWER EXISTS AND NONE CAN. The astrocyte answers have not been bought. This module computes a
CONTROL, a PANEL SHARE and a POWER figure, all of them properties of things that already exist, and
it computes no observed rate of any arm. The fields a run would fill are ABSENT from the payload
rather than present and empty, so nothing here can be read as a preview.

THIS FILE AUTHORISES NOTHING. Albert's second approval names four other conditions beside the one
this registration satisfies, and `astrorun` enforces each on its own. Committing this file makes
`check_astroargmax_registration_is_committed` stop refusing; it does not make a run permitted, and
it moves no cap, ledger, CI check or sign-off.

WHY THE ORDER OF THIS FILE IS THE POINT. The sections are: BANDS, then the label-set RULE in words,
then the ENUMERATION, then the CONTROL, then the DISCLOSURE. Until the bands are committed a control
is an input a band could be drawn around -- a band width chosen after seeing a control can be chosen
to clear it -- so the bands precede the control for the same reason a pre-registration precedes a
run, and not for tidiness. While the enumeration was being built, no compiled-rule count of any term
was consulted: reading the panel's CURIEs and labels is reading the axis the answer comes back on,
but reading which of them carry many compiled rules would let the set lean toward or away from its
own control.

0 AlphaGenome requests. No money. No network to the service. No key is read.
"""

from __future__ import annotations

import collections
import csv
import math
from pathlib import Path
from typing import Any

#: The committed census of every compiled rule's `when: cell_type` label. The CONTROL is read from
#: its per_cell block and from nowhere else; no part of it is computed by this lane.
CENSUS_RELATIVE = "data/results/context_evidence.json"
CENSUS = Path(CENSUS_RELATIVE)

#: The committed dry run whose `requests` list IS the 1,232. Read for the cluster sizes only.
PLAN_RELATIVE = "data/results/astroreg2_request_plan.json"
PLAN = Path(PLAN_RELATIVE)

#: The assay panel the answer's track axis is drawn from. GIT-IGNORED (`.gitignore:20 data/cache`),
#: so it is declared by sha256 and never committed; see THE_PANEL_IS_DECLARED_BY_HASH.
PANEL_RELATIVE = "data/cache/entex/alphagenome_track_metadata_copy.csv"
PANEL = Path(PANEL_RELATIVE)
PANEL_SHA256 = "5157dcca52644952bee2402b03cfa315c576d9f21de229ea9bf9b37e3068abf7"
PANEL_BYTES = 403465

#: The panel rows the gene-axis argmax is taken over.
PANEL_OUTPUT = "rna_seq"

QUESTION = (
    "On the 1,232 paid astrocyte-element answers of AstroREG-2, at what rate is the ARGMAX LABEL -- "
    "the track carrying the largest movement in a gene row's whole track vector -- a CNS track, "
    "against (a) that label set's share of the 440,589 compiled rules and (b) the K562 arm's rate "
    "for the same label set? No astrocyte answer exists, so no rate of either arm exists here"
)

NOT_AN_AUTHORISATION = (
    "this file does not authorise the run and does not spend anything. Albert's second approval has "
    "five clauses -- the cap of 1,232, the separate ledger, today's CI green, THIS registration "
    'committed, and the supervisor\'s words "dry run reviewed" -- and astrorun enforces each one '
    "separately and refuses on each separately. Committing this file satisfies exactly one of them. "
    "No cap, ledger, CI check, sign-off, threshold, test, pin or registration is moved, weakened, "
    "exempted or allowlisted by it"
)

ORDERING_IS_LOAD_BEARING = (
    "the bands come FIRST in this file and were decided before the control was computed. A control "
    "is an input a band could be drawn around: a band width chosen after seeing a control can be "
    "chosen to clear it, and nothing in the artefact would show that it had been. So the order is "
    "bands, then the label-set rule in words, then the enumeration under it, then the control, then "
    "the disclosure. Blindness depends on that order, which is why it is stated rather than implied"
)

# ----------------------------------------------------------------- 1. THE BANDS, before the control

#: The relative floor the enrichment reading must clear: the lower end of the two-sided interval on
#: the RATIO observed-rate / control.
RATIO_ENRICHED_FLOOR = 2.0

#: The absolute level the enrichment reading must ALSO clear, as a share of the 1,232 answers.
ABSOLUTE_LEVEL_FLOOR = 0.10

#: The equivalence band, stated RELATIVE: the two-sided interval on the ratio must lie wholly inside
#: [1 / (1 + EQUIVALENCE_BAND), 1 + EQUIVALENCE_BAND].
EQUIVALENCE_BAND = 0.25

#: The depletion reading's ceiling on the ratio's UPPER end.
RATIO_DEPLETED_CEILING = 0.8

BANDS_BEFORE_THE_CONTROL = (
    "every number in this section was fixed before `control()` was called, and none of them is a "
    "function of the control's value: the two deciding bands are RATIOS to the control and a share "
    "of the answers, so they are the same numbers whatever the control turns out to be. That is the "
    "property an after-the-fact band cannot have"
)

WHY_RELATIVE_AND_NOT_ABSOLUTE = (
    "the band is RELATIVE because an absolute band at a low control is close to a vacuous "
    "equivalence claim. argmaxcell's band was an absolute 0.02 and its own patch records what that "
    "cost: 0.02 is 0.3 times K562's base rate and 91.7 times WTC11's, so one band meant a different "
    "claim on every arm, and on the low arms 'the rate is within 0.02 of its base rate' was "
    "satisfied by anything from zero to many times the base rate. A control near a few per cent has "
    "the same problem here. A ratio band of +/-25% says the same thing at any control value, which "
    "is what makes the equivalence reading worth printing"
)

WHY_A_RATIO_FLOOR_IS_NOT_ENOUGH = (
    "the enrichment reading needs BOTH the ratio floor and the absolute level floor, and either "
    "failing fails it. On 2026-10-02 this project published a 2.9-4.1x relative excess whose "
    "absolute level was 6.0%, so 94% of pairs disagreed with the property the ratio suggested; the "
    "ratio was right and the sentence it supported was not. A level floor of 0.10 states in advance "
    "how much of the data a positive reading has to describe: at least one answer in ten. It is a "
    "declared convention and not a measured threshold, and it is stated as one"
)

WHY_TWO_POINT_ZERO = (
    "2.0 is the smallest ratio floor that is not reachable by this registration's own boundary "
    "decisions. The borderline calls below move the label set by 6 rows (spinal cord), 8 rows (the "
    "CNS-derived cell lines) and 8 rows (the neural-crest-derived lines) out of 667, so a defensible "
    "alternative enumeration moves the set's panel share by up to about a third of itself. A floor "
    "of 1.5 could therefore be cleared or missed by a decision about the pituitary; a floor of 2.0 "
    "cannot. The three alternative enumerations are themselves registered below and reported, so a "
    "reader does not have to take that argument on trust"
)

#: Clusters: whole chromosomes. Below this many, no bootstrap interval is printed and the comparison
#: is INCONCLUSIVE BY RULE.
MIN_CLUSTERS = 10

#: 1 Mb loci, used as clusters only if chromosomes do not reach MIN_CLUSTERS.
LOCUS_SPAN = 1_000_000

BOOTSTRAP_DRAWS = 2000
BOOTSTRAP_SEED = 20261003

#: Below this many DISTINCT resamples the percentile interval is reported as degenerate and decides
#: nothing. The identical-resample share is reported for every interval whatever its value.
MIN_DISTINCT_RESAMPLES = 100

READINGS = {
    "enriched": (
        "the lower end of the two-sided clustered interval on observed_rate / control is at or above "
        f"{RATIO_ENRICHED_FLOOR} AND the observed rate itself is at or above {ABSOLUTE_LEVEL_FLOOR}. "
        "Both, because a ratio without a level misleads. Reading: the argmax carries CNS information "
        "beyond what the compiled-rule label distribution produces on its own -- subject in full to "
        "TRAINING_EXPOSURE, which this design cannot see past"
    ),
    "depleted": (
        f"the upper end of that interval is at or below {RATIO_DEPLETED_CEILING}. Reading: the argmax "
        "lands on a CNS track LESS often than the compiled-rule distribution alone would give. That "
        "is a result and not a failure"
    ),
    "equivalent": (
        "the two-sided interval on the ratio lies wholly inside "
        f"[{1 / (1 + EQUIVALENCE_BAND):.4f}, {1 + EQUIVALENCE_BAND:.4f}]. Reading: on these answers "
        "the argmax label is NOT evidence of the cell the element was measured in; its rate is what "
        "the label distribution produces by itself. That is the finding and it will be written in "
        "those words"
    ),
    "inconclusive": (
        "anything else, including fewer than "
        f"{MIN_CLUSTERS} clusters and a degenerate bootstrap. Worded: THE DATA CANNOT TELL. Not 'no "
        "information' -- an interval that spans both a doubling and a halving carries information "
        "about neither, and saying 'no information' would misreport an instrument that worked"
    ),
}

INSTRUMENT_ASSUMPTIONS = {
    "clusters": (
        "whole chromosomes. The astrocyte arm's 1,232 requests fall on 23 chromosomes, sized 5 to "
        "137 requests, which is in the committed dry run and is read here"
    ),
    "assumed_independent": "chromosomes, and nothing finer",
    "known_not_independent": [
        "requests within a chromosome: the dry run's list is a TILING, 1,232 registry elements over "
        "882 distinct screen-element sets, up to 5 adjacent registry elements serving one screen "
        "element, so neighbouring requests see overlapping sequence",
        "the gene rows within one answer: one deletion, many genes",
        "the 667 track columns within one row: the panel's own rows are 371 distinct names over 2 "
        "strands plus repeats, so a two-strand name supplies two correlated columns to the argmax",
    ],
    "degenerate_bootstrap": (
        "at 0 or n successes a cluster bootstrap is DEGENERATE: every resample returns the same "
        "value, the percentile interval collapses to a point and it decides nothing. The rule is to "
        "route to a Clopper-Pearson exact bound (or Wilson, named which) on an EFFECTIVE sample size "
        "n_eff = n / design_effect, and to print beside it that the bound is NOT the cluster "
        "bootstrap's and that the reading is INCONCLUSIVE unless the band is cleared by the exact "
        "bound alone"
    ),
    "identical_resample_share": (
        "reported for EVERY interval, degenerate or not: the share of the "
        f"{BOOTSTRAP_DRAWS} resamples whose statistic equals the most common resampled value. An "
        "interval whose resamples are mostly identical has fewer effective draws than its draw count "
        f"says, and below {MIN_DISTINCT_RESAMPLES} distinct resamples it is reported as degenerate"
    ),
    "design_effect": (
        "Kish's, from the cluster sizes: (sum n_c)^2 / (C * sum n_c^2) inverted, reported beside "
        "every interval. It is computed from sizes alone and assumes nothing about the outcome"
    ),
}

EVERY_RATIO_CARRIES_ITS_LEVEL = (
    "no ratio is printed anywhere in the result without, in the same record: the absolute level it "
    "is a ratio OF, the control it is a ratio TO, the numerator count and the denominator. The "
    "2026-10-02 failure is the reason -- a 2.9-4.1x excess was published whose absolute level was "
    "6.0%, so 94% of the pairs disagreed with what the ratio was taken to mean. A ratio is a "
    "comparison and never a description of the typical answer"
)

# -------------------------------------------- 2. THE LABEL-SET RULE, in words before any CURIE

LABEL_SET_RULE = (
    "WHAT COUNTS as a CNS track, decided and written before any CURIE of the panel was read for "
    "this purpose. (1) CNS TISSUE: a term naming the brain or any region, lobe, nucleus or part of "
    "it, and the spinal cord or any segment of it. (2) CNS-DERIVED CELLS: a term naming a cell type "
    "of the central nervous system or a progenitor committed to it -- astrocytes, neurons, glia, "
    "neural progenitors, and the neural stem cell and neurosphere terms that name aggregates of "
    "them. WHAT DOES NOT COUNT: peripheral nervous tissue; a term whose stated referent mixes neural "
    "and non-neural tissue and that the panel gives no way to split; a lineage term broad enough to "
    "contain non-neural fates; and a cell line, whose panel label is an opaque code that says "
    "nothing about its tissue of origin on its face. The reading is taken on the set this rule "
    "produces; the borderline calls below are decided individually with their reasons, and three "
    "alternative enumerations are registered so a reader can see what each decision is worth"
)

BORDERLINE_CALLS = {
    "peripheral_nerve_tibial_and_sciatic": {
        "terms": ["UBERON:0001323 tibial nerve", "UBERON:0001322 sciatic nerve"],
        "decision": "OUT",
        "reason": "both are peripheral nerve. The question asks about a brain/astrocyte track and "
        "the rule says CNS; a peripheral nerve is neither, and its resident cells (Schwann cells, "
        "not astrocytes) are a different lineage. Admitting them would make 'CNS' mean 'nervous "
        "system anywhere', which is a different question and a larger label set. They are the "
        "clearest case where a neural-sounding name is not a CNS track",
    },
    "pituitary_gland": {
        "terms": ["UBERON:0000007 pituitary gland"],
        "decision": "OUT",
        "reason": "the gland is two organs in one term. The adenohypophysis is oral ectoderm and not "
        "neural tissue; the neurohypophysis is a down-growth of the diencephalon and is. The panel "
        "carries ONE pituitary term and no lobe terms, so there is no way to take the neural half "
        "and leave the other, and a term that is part non-neural cannot be counted as a CNS track "
        "without overstating what it is. Out, and named here so the decision is visible rather than "
        "an omission",
    },
    "retina": {
        "terms": ["(absent) retina", "UBERON:0000019 camera-type eye"],
        "decision": "OUT, and the first term DOES NOT EXIST on this panel",
        "reason": "retina is developmentally CNS and would have been IN under the rule. The panel "
        "carries no retina term at all -- a MEASURED absence, checked against the 285 terms -- so "
        "the decision has no effect on any count. The eye term that is present is not a substitute: "
        "'camera-type eye' includes cornea, lens and sclera, which are not neural, and it is the "
        "same mixed-referent problem as the pituitary. Recorded because a reader who knows retina is "
        "CNS should be able to see that it was looked for and not found",
    },
    "neuroblastoma_lines": {
        "terms": [
            "EFO:0003072 SK-N-SH",
            "EFO:0005721 SK-N-DZ",
            "EFO:0005725 BE2C",
            "CL:0011012 neural crest cell",
        ],
        "decision": "OUT of the registered set; IN the registered alternative CNS_AND_ALL_NEURAL_LINES",
        "reason": "neuroblastoma arises from the sympathoadrenal lineage of the NEURAL CREST, which "
        "is peripheral: these lines are neural in origin and not central, and the same argument that "
        "excludes the tibial nerve excludes them. The neural crest term itself is excluded for the "
        "stronger reason that the crest also gives melanocytes, craniofacial cartilage and bone, so "
        "it is a lineage term containing non-neural fates. They are 8 panel rows and the alternative "
        "that admits them is enumerated and reported, so the decision is priced rather than hidden",
    },
    "cns_derived_cell_lines": {
        "terms": [
            "EFO:0002101 A172",
            "EFO:0002184 H4",
            "EFO:0005237 U-87 MG",
            "EFO:0005234 PFSK-1",
            "EFO:0005698 Daoy",
        ],
        "decision": "OUT of the registered set; IN the registered alternative CNS_AND_CNS_LINES",
        "reason": "these five are glioma, neuroglioma, glioblastoma, CNS primitive neuroectodermal "
        "and medulloblastoma lines, so by tissue of origin they belong under the rule. What keeps "
        "them out of the REGISTERED set is the rule's last clause: the panel label is 'A172', an "
        "opaque code, and classifying it needs knowledge from outside the committed artefact, which "
        "is precisely the kind of input this registration is trying not to depend on. Excluding them "
        "also biases nothing in the direction of a positive finding. They are 8 rows, enumerated as "
        "an alternative and reported beside the reading",
    },
    "spinal_cord": {
        "terms": ["UBERON:0002240 spinal cord", "UBERON:0006469 C1 segment of cervical spinal cord"],
        "decision": "IN the registered set; OUT of the registered alternative BRAIN_ONLY",
        "reason": "the spinal cord is CNS by the definition of CNS, and the astrocytes the bought "
        "answers are about are CNS astrocytes wherever they sit. The counter-argument is that the "
        "question says 'brain', and 6 of the registered set's 50 rows are cord. Both readings are "
        "defensible, so the decision is IN and the brain-only alternative is enumerated and reported "
        "with its own control and its own bands",
    },
    "broad_lineage_terms": {
        "terms": [
            "CL:0000221 ectodermal cell",
            "CL:0000222 mesodermal cell",
            "CL:0000223 endodermal cell",
            "EFO:0003042 H1",
            "EFO:0003045 H9",
            "EFO:0005904 H7",
            "EFO:0007089 HUES64",
            "EFO:0009747 WTC11",
            "NTR:0000856 mesendoderm",
        ],
        "decision": "OUT",
        "reason": "ectoderm contains skin as well as neural plate, and a pluripotent or germ-layer "
        "term contains every fate. A rule that admitted them would count a label that is CNS among "
        "other things as a CNS label, which is the mixed-referent problem again and at its widest",
    },
    "clo_prefix": {
        "terms": ["CLO:0007045 Jurkat, Clone E6-1"],
        "decision": "OUT",
        "reason": "the panel's ONLY CLO term, and a T-cell leukaemia line. Recorded because the "
        "enumeration is required to cover all five prefixes present on the panel: CLO is covered and "
        "contributes nothing, so the registered set draws on four prefixes. That is a measured fact "
        "about this panel and not a rule about prefixes",
    },
}

NOT_BY_DESCENDANT_RESOLUTION = (
    "the set is an ENUMERATED LIST of CURIEs and not a descendant query, and the argument for that "
    "is specific to this machine and this panel rather than a preference.\n"
    "(1) IT CANNOT BE RUN HERE. There is no UBERON file on this machine: data/knowledge holds "
    "cl-basic.obo and go-basic.obo and nothing else, checked. 283 of the panel's 667 rna_seq rows "
    "carry UBERON terms, so 'descendants of UBERON:0000955' is not a rule that can be evaluated at "
    "all, and a rule that cannot be evaluated is not a rule.\n"
    "(2) IT COULD NOT BE COMMITTED EITHER. data/knowledge is git-ignored (.gitignore:17) and has 0 "
    "tracked files, so the ontology a descendant query resolved against could not travel with the "
    "registration. A later reader would resolve against whatever version their machine had.\n"
    "(3) IT MOVES. 'Descendants of UBERON:0000955' is a different set in different UBERON releases, "
    "so the registered quantity would change without the registration changing -- exactly what a "
    "pre-registration exists to prevent. An enumerated list of 27 CURIEs is fixed by the commit and "
    "is checkable with no ontology file at all.\n"
    "(4) AND IT WOULD DECIDE 22% OF THE PANEL BY ACCIDENT. A UBERON-plus-CL rule leaves every EFO, "
    "NTR and CLO row unclassifiable: 130 + 15 + 1 = 146 of 667 rows, 21.89%, which includes the five "
    "CNS-derived cell lines and the neurosphere term. Those rows would fall out of the label set "
    "because of which ontology the panel happened to spell them in, not because anybody decided "
    "anything about them. The enumeration decides them, in BORDERLINE_CALLS, with reasons that can "
    "be disagreed with."
)

# ------------------------------------------- 3. THE ENUMERATION, fixed before the control was read

ENUMERATION_PRECEDED_THE_CONTROL = (
    "every CURIE below was chosen from the panel's own terms and labels. No compiled-rule count of "
    "any term, and no compiled-rule count of anything, was consulted while the list was being built: "
    "reading the axis is reading what was bought, but reading which terms carry many compiled rules "
    "would let the set lean toward or away from the control it is about to be compared with. The one "
    "exposure that did happen is disclosed in DISCLOSURE, by term and by figure"
)

#: Brain and its parts. 17 terms.
CNS_BRAIN: tuple[tuple[str, str], ...] = (
    ("UBERON:0000955", "brain"),
    ("UBERON:0001870", "frontal cortex"),
    ("UBERON:0001871", "temporal lobe"),
    ("UBERON:0001872", "parietal lobe"),
    ("UBERON:0001873", "caudate nucleus"),
    ("UBERON:0001874", "putamen"),
    ("UBERON:0001876", "amygdala"),
    ("UBERON:0001882", "nucleus accumbens"),
    ("UBERON:0001894", "diencephalon"),
    ("UBERON:0001898", "hypothalamus"),
    ("UBERON:0001954", "Ammon's horn"),
    ("UBERON:0002021", "occipital lobe"),
    ("UBERON:0002037", "cerebellum"),
    ("UBERON:0002038", "substantia nigra"),
    ("UBERON:0002245", "cerebellar hemisphere"),
    ("UBERON:0009834", "dorsolateral prefrontal cortex"),
    ("UBERON:0009835", "anterior cingulate cortex"),
)

#: Spinal cord. 2 terms. IN the registered set, OUT of BRAIN_ONLY.
CNS_SPINAL_CORD: tuple[tuple[str, str], ...] = (
    ("UBERON:0002240", "spinal cord"),
    ("UBERON:0006469", "C1 segment of cervical spinal cord"),
)

#: CNS-derived cells and the committed-progenitor terms. 7 CL terms and 1 NTR term.
CNS_CELLS: tuple[tuple[str, str], ...] = (
    ("CL:0000047", "neuronal stem cell"),
    ("CL:0000100", "motor neuron"),
    ("CL:0000121", "Purkinje cell"),
    ("CL:0000127", "astrocyte"),
    ("CL:0000679", "glutamatergic neuron"),
    ("CL:0002319", "neural cell"),
    ("CL:0011020", "neural progenitor cell"),
    ("NTR:0000427", "neurosphere"),
)

#: CNS-derived cell LINES. OUT of the registered set; the CNS_AND_CNS_LINES alternative adds them.
CNS_CELL_LINES: tuple[tuple[str, str], ...] = (
    ("EFO:0002101", "A172"),
    ("EFO:0002184", "H4"),
    ("EFO:0005234", "PFSK-1"),
    ("EFO:0005237", "U-87 MG"),
    ("EFO:0005698", "Daoy"),
)

#: Neural-crest-derived, i.e. peripheral in origin. OUT; the widest alternative adds them.
NEURAL_CREST_TERMS: tuple[tuple[str, str], ...] = (
    ("CL:0011012", "neural crest cell"),
    ("EFO:0003072", "SK-N-SH"),
    ("EFO:0005721", "SK-N-DZ"),
    ("EFO:0005725", "BE2C"),
)

#: The REGISTERED label set: the reading is taken on this one. 27 terms.
REGISTERED_SET = CNS_BRAIN + CNS_SPINAL_CORD + CNS_CELLS

#: The three alternatives, each enumerated here and each reported with its own control and bands.
#: Registered as alternatives BEFORE the control, so none of them can be promoted afterwards.
ALTERNATIVE_SETS: dict[str, tuple[tuple[str, str], ...]] = {
    "BRAIN_ONLY": CNS_BRAIN + CNS_CELLS,
    "CNS_AND_CNS_LINES": REGISTERED_SET + CNS_CELL_LINES,
    "CNS_AND_ALL_NEURAL_LINES": REGISTERED_SET + CNS_CELL_LINES + NEURAL_CREST_TERMS,
}

WHICH_SET_DECIDES = (
    "CNS_PRIMARY, the 27 terms of REGISTERED_SET, decides the reading. The three alternatives are "
    "reported with their own rates, controls and band verdicts; none of them may be promoted to the "
    "reading after the fact, and a disagreement between them is reported as a disagreement and not "
    "reconciled by choosing one"
)

# ---------------------------------------------------------------------- the panel: axis, not population

AXIS_NOT_POPULATION = (
    "`tracks_total` and `track_names` on an answer are the RESPONSE'S AXIS -- what was bought, in "
    "column order. They are NOT a population and NOT a denominator. The panel this axis is drawn "
    "from carries 667 rna_seq rows under 371 DISTINCT NAMES: 100 names on one row, 246 on two "
    "(strand + and strand -), 25 on more. So the axis has 667 entries and 371 names, and a rate "
    "'per track' and a rate 'per name' are different quantities on the same answer. Neither is a "
    "denominator for the question: the question's denominator is the 1,232 ANSWERS"
)

MEASURED_CORRECTION_TO_A_CODE_CONSTANT = (
    "astrorun.NAMES_ARE_NOT_UNIQUE states that the real RNA-seq metadata carries '667 rows under "
    "316 distinct names', collapsing 351. MEASURED on the file declared by hash here, the figure is "
    "667 rows under 371 distinct names, collapsing 296. 316 is a different quantity: "
    "astroreg_registration.json records it as the number of distinct WINNING names in the chr21 "
    "sweep, '316 distinct winning names and 371 tracks scored per window'. The two were conflated. "
    "Recorded as a measurement because the collapse is what makes comparison (b) an unmatched "
    "comparison (see K562_ARM). astrorun.py is another lane's file and is modified in the working "
    "tree, so nothing of it is touched here: this is the measurement, not an edit"
)

THE_PANEL_IS_DECLARED_BY_HASH = (
    "data/cache/entex/alphagenome_track_metadata_copy.csv is GIT-IGNORED (.gitignore:20, "
    "data/cache), so a FRESH CHECKOUT CANNOT RECOMPUTE the panel share or the enumeration's row "
    "counts from the commit alone. Albert has ruled that data/cache stays local and is never "
    "committed, so declaring it by sha256 is the route -- the one 55 committed results already "
    "take. It is declared here by path, sha256 and byte count, and every figure derived from it is "
    "written into this registration so that what cannot be recomputed can still be CHECKED against "
    "the hash where the file exists, and read where it does not"
)

#: The cross-check's substring patterns and the fields they are matched over, case-folded.
PANEL_SHARE_PATTERNS = ("brain", "astro", "cortex", "cerebellum", "neuron", "glia")
PANEL_SHARE_FIELDS = ("name", "biosample_name")

PANEL_SHARE_IS_A_PANEL_SHARE = (
    "THIS IS A PROPERTY OF THE ASSAY PANEL. It is not a property of the genome, it is not a property "
    "of anything the compiled rules say, and it is NEVER the rule and never a denominator. It is "
    "printed beside the control and labelled, because an unlabelled percentage sitting next to a "
    "compiled-rule rate is how a reader mistakes one for the other -- the label does the work, not "
    "the number. It is reported for one reason: a substring sweep and the enumeration should land in "
    "the same neighbourhood, and where they do not, the disagreement is itself the argument for the "
    "enumeration. They do not: the substring rule misses 20 of the registered set's rows and adds 10 "
    "rows that are not CNS at all, which is reported term by term"
)

# ------------------------------------------------------------- 4. THE CONTROL: COMPILED RULES ONLY

CONTROL_RULE = (
    "control = sum over census per_cell entries whose `ontology_term` is in the label set, of that "
    "entry's `rules`, divided by the census's own `rules` (440,589: every rule compile.compile_"
    "chromosome emits over the 24 chromosomes). The join is by CURIE, which the census carries per "
    "label, so no name matching is involved. No part of the control is computed by this lane from "
    "the astrocyte side, and none of it can have been chosen after seeing an observed figure, "
    "because no observed figure exists"
)

NEVER_A_TRACK_SHARE = (
    "the control is the label set's share of COMPILED RULES. It is NEVER its share of TRACK NAMES. "
    "argmaxcell's registration says why, and it is quoted rather than paraphrased: 'That is "
    "precisely why the base rate and not the track count is the control: a high rate of `label == "
    "K562` on K562-measured elements is what the genome-wide label distribution produces on its own, "
    "before any cell-type information is supposed.' A label's genome-wide frequency already reflects "
    "properties of its tracks -- scale, dynamic range, how often that track wins an argmax anywhere "
    "-- and a track-name share throws all of that away and replaces it with a count of columns. This "
    "is a LIVE error and not a hypothetical: a coordinator used a track share as a control on "
    "2026-10-02, after that registration had already forbidden it in writing"
)

WHAT_A_RATE_AT_THE_CONTROL_MEANS = (
    "stated here so it cannot be softened later: if the observed rate sits inside the equivalence "
    "band around the control, then on these answers the ARGMAX LABEL IS NOT EVIDENCE of the cell "
    "the element was measured in. That is the finding, it will be written in those words, it is not "
    "a failure of the measurement, and it is not a claim that astrocyte regulation is absent "
    "anywhere -- it is a fact about what an argmax over this track panel carries"
)

TRAINING_EXPOSURE = (
    "REGISTERED IN ADVANCE AND NEVER SOFTENED, and it is the strongest limitation on the whole "
    "question. AlphaGenome's tracks were trained on ENCODE and GTEx data from these same biosamples "
    "-- the panel's own `data_source` column says encode on 613 of the 667 rows and gtex on 54 -- so "
    "an agreement between the argmax and a brain label may be RECALL OF SEEN DATA rather than "
    "biology, and THIS DESIGN CANNOT DISTINGUISH THE TWO. There is no held-out brain assay here, no "
    "biosample outside the training distribution, and no second instrument. An enriched reading "
    "therefore supports 'the argmax lands on CNS tracks more often than the compiled-rule "
    "distribution alone would give' and does NOT support 'the model knows this element acts in "
    "astrocytes'. This sentence travels with every quote of any reading of this question, in this "
    "form of words, and may not be shortened to 'a caveat about training data'"
)

K562_ARM = (
    "comparison (b) is the K562 arm's rate for the SAME label set, and it is registered here as a "
    "quantity to be computed at reading time and NOT computed now. Definition: the argmax label of "
    "the all-element deletion sweep's committed per-chromosome element tables, restricted to "
    "elements overlapping the ENCODE K562 CRISPRi screen's valid pairs, classified against the same "
    "27 CURIEs by the same join. THE ASYMMETRY, registered before either rate exists: the sweep's "
    "cached argmax is over NAME-COLLAPSED columns -- 371 names, which is why its own record speaks "
    "of 371 tracks scored per window -- while a paid astrocyte answer keeps all 667 columns, which "
    "is the whole point of item (g). A two-strand name therefore gets two chances to win on the "
    "astrocyte arm and one on the K562 arm, and the collapse is not label-neutral. So the astrocyte "
    "arm's rate is reported BOTH ways: over the 667 columns, and over the axis collapsed to 371 "
    "names exactly as the sweep collapses it, and comparison (b) is read ONLY on the collapsed "
    "pair. The 667-column rate is reported beside it and is not compared to K562"
)

# --------------------------------------------------------------------------------- 5. DISCLOSURE

DISCLOSURE = (
    "WHAT WAS SEEN BEFORE THIS REGISTRATION WAS WRITTEN, AND WHEN. All of it on 2026-10-03, in this "
    "order.\n"
    "(1) BEFORE the bands: astrorun.py's approval clauses, its track-count constants (371, 667) and "
    "its item-(g) vector checks; data/results/astroreg2_registration.json and "
    "astroreg2_request_plan.json in full, including the 1,232/1,322/90 split and the per-chromosome "
    "request counts; data/results/argmaxcell_registration.json in full, which is the model for this "
    "file and from which the base rates of its five arms were seen (K562 27,445/440,589 among them); "
    "and compile.context.\n"
    "(2) THEN the bands were decided and written. They are ratios and a share of the answers, so no "
    "figure seen above fixes any of them.\n"
    "(3) THEN the panel: the 667 rna_seq rows of the file declared by hash, their 285 CURIEs with "
    "labels and row counts, the 371-name collapse and the substring cross-check. This is the axis.\n"
    "(4) THEN the enumeration was fixed, and the four sets' term lists and panel row counts were "
    "computed: registered set 27 terms / 50 rows, BRAIN_ONLY 25 / 44, CNS_AND_CNS_LINES 32 / 58, "
    "CNS_AND_ALL_NEURAL_LINES 36 / 66.\n"
    "(5) ONE EXPOSURE, DISCLOSED RATHER THAN MINIMISED. While reading the census's SHAPE to find out "
    "whether per_cell entries carry a CURIE at all, the first eight entries printed, and one of them "
    "is in the registered set: neural_cell (CL:0002319) at 9,955 of 440,589 rules. The seven others "
    "-- K562, CD14-positive monocyte, HepG2, placenta, Whole_Blood, psoas_muscle, small_intestine -- "
    "are not in it. The enumeration was ALREADY FIXED when that printed, and it was not revised "
    "afterwards: no term was added or removed after step (4). The exposure is recorded because a "
    "reader cannot check that from the file and has to be told.\n"
    "(6) THEN the control, the alternatives' controls and the power figures were computed, which is "
    "this file.\n"
    "WHAT WAS NOT SEEN, AND DOES NOT EXIST: no astrocyte answer, no astrocyte argmax, no astrocyte "
    "rate, no K562 argmax distribution of any kind, and no winning-name distribution from the sweep. "
    "Not one of the 1,232 requests has been sent."
)

# ------------------------------------------------------------------------- what this cannot claim

CANNOT_ESTABLISH = (
    "whether an agreement between the argmax and a CNS label is biology or recall of training data: "
    "see TRAINING_EXPOSURE. No outcome of this question settles it",
    "whether the element acts in astrocytes. An argmax is the place the largest movement happened to "
    "fall, over a panel this lane did not choose; it is a SELECTION and not a measurement of a cell",
    "whether the label set is the right one. Three alternatives are registered because it is a "
    "judgement, and the reading is reported against all four",
    "anything about the 2,174 non-rna_seq rows of the panel. The question is about the gene-axis "
    "argmax, which is taken over rna_seq columns",
    "any statement needing a biosample or experiment accession: the declared metadata copy carries "
    "neither, so 'one biosample NAME' is the most that can be said of a track and this lane says "
    "less than that",
)

REFUSALS = (
    "this registration is committed ALONE. No count, no result, no answer, no rate, no interval",
    "the fields a run would fill are ABSENT from the payload rather than present and empty, so "
    "nothing in it can be read as a preview of an answer that does not exist",
    "no threshold, floor, cap, pin, test or registration of any other lane is moved, weakened, "
    "exempted or allowlisted, including astrorun's five clauses and its 1,232 cap",
    "no AlphaGenome request, no network to the service, no key read, no money",
)


# ------------------------------------------------------------------------------------- computation


def panel_rows() -> list[dict[str, str]]:
    """The panel's rna_seq rows. Reads the axis and nothing else; refuses on a hash mismatch."""
    from genomeos import manifest as mf

    digest, size, _ = mf.sha256_of(PANEL)
    if digest != PANEL_SHA256 or size != PANEL_BYTES:
        raise ValueError(
            f"{PANEL_RELATIVE} is {size} bytes with sha256 {digest}; this registration declares "
            f"{PANEL_BYTES} bytes and {PANEL_SHA256}. The panel is git-ignored, so the hash is the "
            "only thing that identifies it, and a different file is a different axis"
        )
    with PANEL.open(newline="", encoding="utf-8") as fh:
        return [r for r in csv.DictReader(fh) if r.get("output") == PANEL_OUTPUT]


def panel_summary(rows: list[dict[str, str]]) -> dict[str, Any]:
    """The axis's shape: rows, distinct names, how the names collapse, and where they come from."""
    per_name = collections.Counter(r["name"] for r in rows)
    return {
        "rows": len(rows),
        "distinct_names": len(per_name),
        "rows_that_collapse_under_name": len(rows) - len(per_name),
        "names_with_one_row": sum(1 for v in per_name.values() if v == 1),
        "names_with_two_rows": sum(1 for v in per_name.values() if v == 2),
        "names_with_more_than_two_rows": sum(1 for v in per_name.values() if v > 2),
        "distinct_ontology_terms": len(set(r["ontology_curie"] for r in rows)),
        "rows_by_prefix": dict(
            sorted(collections.Counter(r["ontology_curie"].split(":")[0] for r in rows).items())
        ),
        "rows_by_data_source": dict(sorted(collections.Counter(r["data_source"] for r in rows).items())),
        "axis_not_population": AXIS_NOT_POPULATION,
        "measured_correction": MEASURED_CORRECTION_TO_A_CODE_CONSTANT,
    }


def enumerate_set(terms: tuple[tuple[str, str], ...], rows: list[dict[str, str]]) -> dict[str, Any]:
    """One label set, term by term, with each term's row count and name count ON THE PANEL."""
    by_term = collections.Counter(r["ontology_curie"] for r in rows)
    names: dict[str, set[str]] = collections.defaultdict(set)
    for r in rows:
        names[r["ontology_curie"]].add(r["name"])
    listed = [
        {
            "curie": curie,
            "label": label,
            "panel_rows": by_term.get(curie, 0),
            "panel_names": len(names.get(curie, ())),
        }
        for curie, label in terms
    ]
    total = sum(e["panel_rows"] for e in listed)
    return {
        "terms": len(listed),
        "included": listed,
        "absent_from_the_panel": [e["curie"] for e in listed if e["panel_rows"] == 0],
        "panel_rows": total,
        "panel_names": sum(e["panel_names"] for e in listed),
        "panel_row_share": round(total / len(rows), 6) if rows else None,
        "panel_row_share_is_a_panel_share": PANEL_SHARE_IS_A_PANEL_SHARE,
        "prefixes_represented": sorted(set(e["curie"].split(":")[0] for e in listed)),
    }


def panel_substring_crosscheck(
    terms: tuple[tuple[str, str], ...], rows: list[dict[str, str]]
) -> dict[str, Any]:
    """The substring sweep, and every row it disagrees with the enumeration about.

    Reported as a CROSS-CHECK and labelled a panel share. The disagreement is the point.
    """
    members = set(c for c, _ in terms)

    def matched(r: dict[str, str]) -> bool:
        text = " ".join((r.get(f) or "") for f in PANEL_SHARE_FIELDS).lower()
        return any(p in text for p in PANEL_SHARE_PATTERNS)

    hit = [r for r in rows if matched(r)]
    fp = collections.Counter(
        (r["ontology_curie"], r["biosample_name"]) for r in hit if r["ontology_curie"] not in members
    )
    miss = collections.Counter(
        (r["ontology_curie"], r["biosample_name"])
        for r in rows
        if r["ontology_curie"] in members and not matched(r)
    )
    return {
        "label": "PANEL SHARE -- a property of the assay panel, not of the genome and not of the rules",
        "patterns": list(PANEL_SHARE_PATTERNS),
        "fields": list(PANEL_SHARE_FIELDS),
        "matching_rows": len(hit),
        "panel_rows": len(rows),
        "panel_share": round(len(hit) / len(rows), 6) if rows else None,
        "never_the_rule": PANEL_SHARE_IS_A_PANEL_SHARE,
        "disagrees_with_the_enumeration_on": sum(fp.values()) + sum(miss.values()),
        "matched_but_not_cns": [
            {"curie": c, "label": lab, "panel_rows": n, "why": "a substring, not a referent"}
            for (c, lab), n in sorted(fp.items())
        ],
        "cns_but_not_matched": [
            {"curie": c, "label": lab, "panel_rows": n} for (c, lab), n in sorted(miss.items())
        ],
        "what_the_disagreement_shows": (
            "'astro' is a substring of 'gastrocnemius' and of 'gastroesophageal', and 'cortex' is a "
            "substring of 'cortex of kidney' and of three renal cortex interstitium terms: that is "
            "10 rows a substring rule calls brain and no reader would. In the other direction the "
            "rule misses spinal cord, four lobes of the brain, the diencephalon, Purkinje cells, "
            "neural cells, neural progenitor cells and neurospheres: 20 rows. A rule that is wrong "
            "in both directions on 30 of 667 rows is a cross-check and cannot be the label set"
        ),
    }


def census() -> dict[str, Any]:
    """The committed census, with its per_cell block's sum checked against its own `rules`."""
    from genomeos.results import load_result

    cen = load_result("context_evidence")
    if not cen:
        raise ValueError(f"{CENSUS_RELATIVE} is absent; the control has no source")
    per = cen.get("per_cell") or {}
    total = int(cen["rules"])
    summed = sum(int(v["rules"]) for v in per.values())
    if summed != total:
        raise ValueError(
            f"the census's per_cell entries sum to {summed} and it reports {total} rules. The "
            "control's denominator and its numerator would then come from different populations, so "
            "this registration refuses to be written"
        )
    return {"rules": total, "per_cell": per, "distinct_labels": len(per)}


def control(terms: tuple[tuple[str, str], ...], cen: dict[str, Any]) -> dict[str, Any]:
    """The label set's share of COMPILED RULES. Never of track names: see NEVER_A_TRACK_SHARE."""
    members = set(c for c, _ in terms)
    per = cen["per_cell"]
    matched = {
        label: int(entry["rules"]) for label, entry in per.items() if entry.get("ontology_term") in members
    }
    found = set(per[label].get("ontology_term") for label in matched)
    rules = sum(matched.values())
    return {
        "rate": round(rules / cen["rules"], 8),
        "rules_with_a_label_in_the_set": rules,
        "rules_total": cen["rules"],
        "census_labels_matched": dict(sorted(matched.items(), key=lambda kv: -kv[1])),
        "census_labels_matched_count": len(matched),
        "curies_in_the_set_with_no_census_label": sorted(members - found),
        "census_labels_carrying_no_ontology_term": sum(1 for e in per.values() if not e.get("ontology_term")),
        "how": CONTROL_RULE,
        "never_a_track_share": NEVER_A_TRACK_SHARE,
        "source": CENSUS_RELATIVE,
    }


def cluster_sizes() -> dict[str, int]:
    """The 1,232 requests per chromosome, from the committed dry run's own list. Axis, not outcome."""
    from genomeos.results import load_result

    plan = load_result("astroreg2_request_plan")
    if not plan:
        raise ValueError(f"{PLAN_RELATIVE} is absent; the cluster sizes have no source")
    reqs = plan.get("requests") or []
    return dict(sorted(collections.Counter(str(r["chrom"]) for r in reqs).items()))


def design(sizes: dict[str, int]) -> dict[str, Any]:
    """Kish's effective cluster count and the design effect, from sizes alone."""
    vals = list(sizes.values())
    n = sum(vals)
    sq = sum(v * v for v in vals)
    n_eff_clusters = (n * n) / sq if sq else 0.0
    return {
        "clusters": len(vals),
        "requests": n,
        "smallest_cluster": min(vals) if vals else 0,
        "largest_cluster": max(vals) if vals else 0,
        "effective_clusters_kish": round(n_eff_clusters, 4),
        "design_effect_from_sizes": round(len(vals) / n_eff_clusters, 4) if n_eff_clusters else None,
        "how": "n_eff = (sum n_c)^2 / sum n_c^2 over the per-chromosome request counts; the design "
        "effect is C / n_eff. Both are computed from cluster SIZES only and assume nothing about any "
        "outcome, so they are stated in advance",
        "clusters_are": INSTRUMENT_ASSUMPTIONS["clusters"],
    }


def clusters_needed(p_alt: float, p_null: float, rho: float, mean_size: float) -> int:
    """Chromosome-clusters needed for the lower 2.5% bound on p_alt to clear p_null, at icc `rho`.

    A normal approximation on the log rate with the standard cluster inflation
    1 + (m - 1) * rho. Stated in CLUSTER units because the clusters are what the interval resamples:
    an element count would overstate the information in 23 correlated chromosomes.
    """
    if not (0 < p_null < p_alt < 1):
        raise ValueError(f"needs 0 < p_null {p_null} < p_alt {p_alt} < 1")
    z = 1.959963984540054
    inflation = 1.0 + (mean_size - 1.0) * rho
    gap = math.log(p_alt) - math.log(p_null)
    var_per_element = (1.0 - p_alt) / p_alt
    n = (z * z) * var_per_element * inflation / (gap * gap)
    return max(1, math.ceil(n / mean_size))


def power_in_cluster_units(control_rate: float, sizes: dict[str, int]) -> dict[str, Any]:
    """Stated in advance, per comparison, in the units of the clusters.

    Placed AFTER the control in the payload and computed from it: the BANDS are fixed independently
    of the control, and these figures are what those fixed bands cost in chromosomes.
    """
    d = design(sizes)
    mean_size = d["requests"] / d["clusters"]
    alt = min(0.999, control_rate * RATIO_ENRICHED_FLOOR * 1.25)
    null = control_rate * RATIO_ENRICHED_FLOOR
    rows = {}
    for rho in (0.0, 0.05, 0.2):
        rows[f"icc_{rho}"] = {
            "clusters_needed": clusters_needed(alt, null, rho, mean_size),
            "clusters_available": d["clusters"],
            "effective_clusters_kish": d["effective_clusters_kish"],
        }
    return {
        "comparison_a_vs_the_control": {
            "alternative_registered_in_advance": round(alt, 8),
            "alternative_is": f"{RATIO_ENRICHED_FLOOR} x {1.25} x the control, i.e. the ratio floor "
            "cleared with 25% to spare. Declared here so the power figure describes a stated effect "
            "and not whatever is found",
            "null_the_lower_bound_must_clear": round(null, 8),
            "absolute_level_floor_must_also_be_cleared": ABSOLUTE_LEVEL_FLOOR,
            "mean_cluster_size": round(mean_size, 4),
            "per_intracluster_correlation": rows,
            "units": "whole chromosomes, which are the clusters the interval resamples",
            "how": "normal approximation on the log rate with cluster inflation 1 + (m - 1) * rho, "
            "z = 1.96 two-sided. It is an approximation and is reported as one; the interval that "
            "DECIDES is the cluster bootstrap, and where that is degenerate, the exact bound",
            "and_the_honest_part": (
                "the astrocyte arm has 23 chromosomes and an effective "
                f"{d['effective_clusters_kish']} of them by Kish, because the sizes run from "
                f"{d['smallest_cluster']} to {d['largest_cluster']} requests. Where the table above "
                "asks for more clusters than 23, the comparison is UNDERPOWERED BY REGISTRATION and "
                "its reading will be INCONCLUSIVE -- the data cannot tell -- whatever the point "
                "estimate turns out to be. That is registered now rather than discovered later"
            ),
        },
        "comparison_b_vs_the_K562_arm": {
            "what_is_compared": "two rates, on the NAME-COLLAPSED axis only: see K562_ARM",
            "clusters": "whole chromosomes on both arms; the astrocyte arm's 23 are known and sized "
            "here, the K562 arm's count and size distribution are NOT known to this registration "
            "and are reported with its rate",
            "power_cannot_be_stated_as_one_number_here": (
                "a two-sample power figure needs both arms' cluster counts and sizes, and the K562 "
                "arm's are not known without reading its element set, which this registration does "
                "not read. What IS registered in advance: the SAME bands apply, the comparison is "
                f"INCONCLUSIVE BY RULE if either arm has fewer than {MIN_CLUSTERS} clusters, and the "
                "clusters_needed table above is reported for both arms at their own mean cluster "
                "sizes once those are known. Stating a number here that depends on a quantity this "
                "lane has not read would be a figure with nothing behind it"
            ),
        },
        "degenerate_bootstrap": INSTRUMENT_ASSUMPTIONS["degenerate_bootstrap"],
        "identical_resample_share": INSTRUMENT_ASSUMPTIONS["identical_resample_share"],
    }
