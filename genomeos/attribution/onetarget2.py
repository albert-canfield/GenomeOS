# SPDX-License-Identifier: AGPL-3.0-or-later
"""The remaining one-target consumers: the census, and one registration with a falsifier per module.

Every deletion run this project has made writes a COMPACT table: one target gene per element, the
strongest mover in the scorer's window. `attribution/targets.py` holds both readers -- `run_elements`
(the compact table) and `ElementResponses` (the sweep's own per-element cache, where every gene in
the window is kept with its signed change) -- and `predict/enhancer_target.cached_prediction`
projects that same window down to one gene for the layers that never call the API.

Nine consumers have been moved off the compact head one at a time, each with its own registration.
This module registers the REST OF THEM TOGETHER, under one registration, with a falsifier that can
fail on its own for each module. That is the point of doing it once: a single falsifier for a set
would let most of the set pass on the strength of one member.

What the projection does and does not lose is already measured and is not re-opened here. The head
of the window's at-the-bar list IS the compact table's `predicted` gene wherever the table names
one, and the list is empty exactly when the table names none -- 849,469 elements compared, 0 head
disagreements and 0 coding-head disagreements (`data/results/onetarget_consumers.json`, 2026-09-27,
commit d717b28). So a consumer that asks only WHETHER an element moved something, or only HOW LARGE
the strongest move was, cannot change. What a one-target reader loses is WHICH genes and HOW MANY,
and only a consumer whose own output depends on that is a one-target consumer. Two modules in the
census are invariant by construction for exactly this reason and are registered as NOT to be moved:
moving them would change nothing, or would move a number with no reason to move it.

This file is a REGISTRATION. It holds no measurement of its own and no field a run would fill.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

REGISTERED = "2026-10-03"
LANE = "lane-onetarget"
WAVE = "section 5 item 10, wave 2"

#: The row this registration answers, quoted from docs/ROADMAP.md's wave table.
ROW = (
    "the remaining one-target consumers under **one** registration with a falsifier per module | "
    "five modules, one shared reader; one registration instead of five | 3"
)

QUESTION = (
    "Which modules still read ONE target gene per element, what does each of them claim on the "
    "strength of that one gene, and what single figure of each one's own output would have to move "
    "for its reading to be wrong? One registration, one falsifier per module, no module's claim "
    "widened or dropped on the way in."
)

#: The count the row states, and the count the tree holds. The row's figure is the coordinator's
#: estimate; this census is the measurement of the scope, taken before any consumer is moved.
#: The scope, declared as numbers so a guard can re-derive every one of them from CENSUS and refuse
#: a figure that has drifted. These are SCOPE, not a result: they count modules, not elements.
SCOPE_COUNTS = {
    "row_says": 5,
    "census_modules": 20,
    "identity": 18,
    "invariant": 2,
    # 5 as registered 2026-10-03; 6 since the coordinator's ruling the same day on
    # attribution/silenceragree.py, recorded in RECORD_ONLY below and in the prose amendment.
    "record_only": 6,
    "already_moved": 9,
    "excluded_by_name": 9,
    "identity_free_today": 11,
    "identity_peer_held": 5,
    "identity_frozen_closure": 2,
    "closure_files_reading_a_head": 3,
    "scripts_out_of_scope": 24,
}

ROW_SAYS_FIVE_THE_TREE_SAYS_OTHERWISE = (
    "The row says five modules. The tree does not. Under `genomeos/` there are 18 modules whose own "
    "output depends on WHICH or HOW MANY genes an element targets (class `identity`) and 2 that "
    "read a head and cannot move (class `invariant`), making a census of 20; beside them 5 modules "
    "read a head without being a live consumer of one (class `record`) and 9 have already been "
    "moved, each under its own registration. **Amended 2026-10-03, the sentence above kept as "
    "written: the record count is 6, not 5. `attribution/silenceragree.py` was ruled record-only "
    "that day, on the same ground as `respmap_v2.py` -- it audits the published head rather than "
    "consuming it, and re-reading the window would measure a quantity the project does not "
    "publish. The ruling does not move any figure of this registration's own wave, which is still "
    "18.** `scripts/` holds 24 further head-reading entry points, "
    "named below and OUT of this registration's scope. So the real figure for this wave is 18, not "
    "5. Of those 18, two are frozen inside the paid study's 52-file import closure "
    "(`genome/motifs.py`, `genome/regulation.py`) and a third closure file, "
    "`predict/enhancer_target.py`, is the projection itself; five are held on the work board by "
    "live peer lanes, which the commit gate refuses rather than a decision refusing; leaving 11 "
    "free to move today. The roadmap row is not edited by this lane: the count is reported here and "
    "goes to the coordinator as a row."
)

EXPLORATION_PRECEDED_THIS = (
    "This lane read the tree, the committed 2026-09-27 result, the four earlier one-target "
    "registrations and the roadmap's own account of them BEFORE writing this file. This "
    "registration is therefore not a prediction of the census figures above: the census is the "
    "scope, taken by `head_read_sites()` over the tree, and `check_census_matches_tree()` refuses "
    "if the declared lines are not the lines the tree holds. What this file binds and cannot know "
    "is each module's falsifier: no figure of any consumer's output has been read through the "
    "window, and no consumer has been re-run."
)

# ---------------------------------------------------------------------------
# the one shared reader
# ---------------------------------------------------------------------------

THE_SHARED_READER = {
    "module": "genomeos/attribution/targets.py",
    "compact": "run_elements(name, chrom, results_dir) at targets.py:40 -- one run's elements from "
    "the compact table, or from the local table a summary's `elements_where` points at; a summary "
    "whose table is missing raises rather than reading as a run with no elements",
    "deduplicated": "attributed(chrom, results_dir, coding=True) at targets.py:59 -- every element "
    "with a named (coding) head across the three runs, by id, with its run under `origin`. This is "
    "the function most of the census reaches the head through",
    "window": "ElementResponses at targets.py:127 -- the sweep's per-element response cache: every "
    "gene it scored at an element, with a NAMED SILENCE (not in the window, not on that track, not "
    "cached) where it has no number, never a zero standing in for a missing answer",
    "at_the_bar": "genes_at_bar(record, min_effect) at targets.py:79 -- every gene of one cached "
    "element reaching the bar, strongest first, by `predict_target`'s own size rule",
    "projection_for_layers_that_never_call_the_api": "cached_prediction(chrom, element_id) at "
    "genomeos/predict/enhancer_target.py:589 -- it loads the WHOLE window and returns "
    "predict_target(hit['genes']), one gene. Two consumers of it already hold the window in hand "
    "and discard it, so their fix needs no new reader and no request",
    "one_registration_not_five": "every module below reads one of these four entry points, so one "
    "reader is shared and one registration covers them; the falsifiers are not shared",
}

THE_INVARIANT = (
    "the window's at-the-bar head is the compact table's `predicted` gene wherever the table names "
    "one, and the at-the-bar list is empty exactly when the table names none: 849,469 elements, 0 "
    "head disagreements, 0 coding-head disagreements (data/results/onetarget_consumers.json, "
    "2026-09-27). Carried here word for word as that result registered it and NOT strengthened: it "
    "is a measurement over the elements that run's arm compared, not a theorem about every element "
    "on every chromosome, and the same result records 10 elements of 4,794 with no "
    "`predicted_coding` field whose cache DOES name a coding gene, so the coding head is not "
    "invariant in the way the any-gene head is"
)

WHAT_THE_PROJECTION_LOSES = (
    "the identity and the count of the other genes at the bar, and the per-gene silence. It does "
    "not lose whether the element moved something, and it does not lose the magnitude of the "
    "strongest move, which is the window maximum by `predict_target`'s own rule"
)

# ---------------------------------------------------------------------------
# the census
# ---------------------------------------------------------------------------

#: The two spellings a compact head is read under. A site is a `.get("<key>")` or a `["<key>"]`.
HEAD_KEYS = ("predicted", "predicted" + "_coding")

#: Files of the paid study's 52-file import closure. `sender_closure()` hashes the WORKING TREE, so
#: even an uncommitted edit to one of these refuses the frozen send. They are named here so that the
#: wave cannot reach them by accident, and `check_no_closure_file_is_moving()` refuses if it tries.
CLOSURE_FROZEN = (
    "genomeos/genome/motifs.py",
    "genomeos/genome/regulation.py",
    "genomeos/predict/enhancer_target.py",
)

#: Classes. `identity`: the module's own output depends on which or how many genes the element
#: targets, so it is a one-target consumer. `invariant`: it reads a head and cannot move, because it
#: reads only the yes/no or only the head's magnitude. `record`: it reads a head without being a
#: live consumer of one -- it is the projection itself, or it audits a quote, or it exists to
#: document the defect.
CLASSES = ("identity", "invariant", "record")

#: Statuses. `free`: this wave may move it. `frozen_closure`: a closure file, Albert's call, not
#: this lane's. `peer_held`: a live peer lane holds the path on the work board; the commit gate
#: refuses a staged path a peer holds, so these wait for the holder, not for a decision.
STATUSES = ("free", "frozen_closure", "peer_held")

CENSUS: tuple[dict[str, Any], ...] = (
    {
        "module": "genomeos/attribution/closure.py",
        "lines": (73,),
        "reads": "targets.attributed -> e['predicted_coding']",
        "claims_today": "the elements of a gene's regulatory closure, grouped BY the head gene's "
        "symbol, each with its per-cell signed value",
        "class": "identity",
        "status": "free",
        "falsifier": "a gene whose closure gains or loses an element when the elements are grouped "
        "by every gene at the bar instead of by the head alone. If no gene's element set changes, "
        "the closure was never one-target in effect and this module is withdrawn from the wave.",
    },
    {
        "module": "genomeos/attribution/compile.py",
        "lines": (332, 838, 880),
        "reads": "cp._attributed -> e['predicted_coding'] at three sites: the rule's action (332), "
        "the compiled gene list (838) and the rule's basis sentence (880)",
        "claims_today": "one `predicted_deletion_target` rule per element naming one gene, the "
        "chromosome's gene list as the distinct head genes plus the measured genes, and an entity "
        "count and three `# test:` lines computed from those",
        "class": "identity",
        "status": "peer_held",
        "falsifier": "a compiled chromosome whose `# test: entities >= N` or `# test: rules == N` "
        "line changes under the window reading. Those two lines are the program's own self-test, so "
        "a change there is a change to what the committed programs assert about themselves, and it "
        "is the single figure that decides whether this module may be moved in place at all.",
    },
    {
        "module": "genomeos/attribution/organise.py",
        "lines": (164,),
        "reads": "targets.attributed -> e['predicted_coding']['gene']",
        "claims_today": "one `target` per attributed element in each organised block, which is what "
        "the block organiser's `attributed_elements` and `targets` counts are built from",
        "class": "identity",
        "status": "peer_held",
        "falsifier": "the per-block `targets` count on any of the 24 organised chromosomes moving "
        "by more than the number of elements that gain a second at-the-bar gene. A larger move "
        "means the window changed which elements are attributed at all, not only how many genes "
        "each names, and that would contradict THE_INVARIANT rather than extend it.",
    },
    {
        "module": "genomeos/attribution/variation.py",
        "lines": (197,),
        "reads": "the two sampled runs' inline `elements` -> e['predicted_coding']",
        "claims_today": "one `target` and one `log2_fold_change` per element on the human variation "
        "axis, from the constrained and uniform samples only",
        "class": "identity",
        "status": "free",
        "falsifier": "any row's `log2_fold_change` changing. The head's magnitude is the window "
        "maximum by construction, so this number must NOT move; if it does, `predict_target`'s size "
        "rule and `genes_at_bar`'s are not the same rule and every magnitude in the census is in "
        "doubt.",
    },
    {
        "module": "genomeos/attribution/measured.py",
        "lines": (1092,),
        "reads": "targets.attributed -> e['predicted_coding']",
        "claims_today": "a `predicted_gene` set beside each element's MEASURED response, which is "
        "the field any later comparison of model against measurement reads",
        "class": "identity",
        "status": "free",
        "falsifier": "an element whose measured regulated gene is at the bar in the window but is "
        "not the head. Every such element is one where this module today says the model named a "
        "different gene than the measurement, and the window says it named this one too. One such "
        "element falsifies the reading; zero of them, over the elements the layer covers, says the "
        "one-target projection costs this comparison nothing and the module is withdrawn.",
    },
    {
        "module": "genomeos/attribution/not_open_profile.py",
        "lines": (254,),
        "reads": "cp._attributed -> e['predicted_coding']",
        "claims_today": "one predicted `Rule` per element, carrying the head gene and the head's "
        "`log2_fold_change` as the rule's effect",
        "class": "identity",
        "status": "free",
        "falsifier": "the number of predicted rules changing. One element must yield one rule per "
        "gene it names, so the rule COUNT is expected to rise; a count that does not rise at all "
        "means no element in this profile's population names a second gene, and the module is "
        "withdrawn from the wave on its own evidence.",
    },
    {
        "module": "genomeos/attribution/gwas.py",
        "lines": (153, 154, 155, 164, 166),
        "reads": "the caller's `elements` -> e['predicted'] (any-gene head, not the coding head)",
        "claims_today": "the share of GWAS-hit elements whose predicted target AGREES with the "
        "catalogue's mapped gene, computed against the head gene only",
        "class": "identity",
        "status": "free",
        "falsifier": "the agreement share not rising. This is the one module in the census whose "
        "headline figure can only move one way under the window -- an element agrees if ANY gene at "
        "the bar is the mapped gene, which is a superset of agreeing on the head -- so a share that "
        "falls or an element that loses its agreement is a defect in the move, not a finding, and "
        "stops the wave.",
    },
    {
        "module": "genomeos/attribution/argmaxcell.py",
        "lines": (677, 692),
        "reads": "a streaming decode of the compact table, and small_run_elements -> "
        "e['predicted_coding'] for its `tissue` and as the gate on a named coding gene",
        "claims_today": "one CONTEXT per element, from the head gene's strongest track, and an "
        "element with no coding head contributes none",
        "class": "identity",
        "status": "peer_held",
        "falsifier": "an element whose at-the-bar genes do not all resolve to the one context the "
        "head resolves to. If every gene at an element agrees on the context, then context is a "
        "property of the element and not of the chosen gene, and this module is invariant after all "
        "-- which is a finding and moves it out of `identity`.",
    },
    {
        "module": "genomeos/attribution/context_evidence.py",
        "lines": (441,),
        "reads": "cp._attributed -> e['predicted_coding'] for its `tissue` only",
        "claims_today": "one (context, start, end) per attributed element, from the head gene's "
        "strongest track, pooled with the measured rows' own cells",
        "class": "identity",
        "status": "peer_held",
        "falsifier": "the per-context element counts changing. This module's figure is a base rate "
        "over contexts, so a count that moves says the context mix of the attributed elements is a "
        "property of the projection; a count that does not move says it is not, and the module "
        "joins `invariant`.",
    },
    {
        "module": "genomeos/attribution/vista.py",
        "lines": (260, 277, 278, 280, 312, 326),
        "reads": "enhancer_target.score_element's output -> row['predicted'] for `predicted_group` "
        "and `tissue_agrees` (277, 278, 280); the over-length skip writes a None head (260); and "
        "the summary counts `with_predicted_target` and `strong` off the head (312, 326)",
        "claims_today": "whether the predicted track GROUP of one gene agrees with the VISTA "
        "element's observed expression groups, and the per-status share of elements with a "
        "predicted target at all",
        "class": "identity",
        "status": "free",
        "falsifier": "a VISTA row whose `tissue_agrees` is False on the head and True for some "
        "other gene at the bar. One such row falsifies the stored agreement rate as an agreement "
        "rate of the ELEMENT; zero of them leaves the stored rate standing unchanged.",
    },
    {
        "module": "genomeos/decompile.py",
        "lines": (124,),
        "reads": "the two sampled runs' inline `elements`, SELECTING on "
        "e['predicted_coding']['gene'] == the requested symbol",
        "claims_today": "the elements the decompiler shows for a gene: exactly those whose head IS "
        "that gene, so an element that moves the gene as a second mover is invisible",
        "class": "identity",
        "status": "free",
        "falsifier": "a gene whose decompiled element set is empty today and non-empty under the "
        "window. That is the strongest form of this defect -- a gene the decompiler says has no "
        "scored element when the sweep scored one at the bar for it -- and one instance establishes "
        "it; none means the selection never excluded anything.",
    },
    {
        "module": "genomeos/genome/motifs.py",
        "lines": (522,),
        "reads": "the two sampled runs' inline `elements` -> e['predicted_coding']['gene']",
        "claims_today": "one `target` per element in the motif population, which the motif "
        "enrichment is grouped by",
        "class": "identity",
        "status": "frozen_closure",
        "falsifier": "an enriched motif whose target set changes under the window. NOT RUN THIS "
        "WAVE: this file is in the 52-file import closure and even an uncommitted edit to it "
        "refuses the frozen paid send, so the falsifier is registered and left unfired. It is "
        "Albert's call, not this lane's.",
    },
    {
        "module": "genomeos/genome/regulation.py",
        "lines": (163,),
        "reads": "enhancer_target.cached_prediction(chrom, e.id) -> pred['gene'], then "
        "`predicted_this_gene = pred['gene'] == g.symbol`",
        "claims_today": "whether the model predicts THIS gene at this enhancer -- decided on one "
        "gene, although `cached_prediction` loaded the whole window one line earlier and threw it "
        "away. The fix here costs 0 requests and needs no new reader",
        "class": "identity",
        "status": "frozen_closure",
        "falsifier": "an enhancer whose `predicted_this_gene` is False and whose window has the "
        "gene at the bar. NOT RUN THIS WAVE: closure file, Albert's call. Registered here because "
        "this is the clearest instance in the census of the window being in hand and discarded, and "
        "it must not be lost by being blocked.",
    },
    {
        "module": "genomeos/genome/regdiff.py",
        "lines": (101, 176, 177, 194, 195),
        "reads": "enhancer_target.cached_prediction at regdiff.py:74 writes one "
        "{gene, log2_fold_change, tissue} per variant's element at regdiff.py:85 (a dict KEY, so "
        "not a site the scan matches); the head is then read as a rank term (101), as the two "
        "`with_predicted_target_only_*` counts (176, 177) and as the printed line (194, 195)",
        "claims_today": "`with_predicted_target_only_a` / `_only_b`: how many of two individuals' "
        "differing regulatory variants sit in an element with a predicted target, and which gene",
        "class": "identity",
        "status": "free",
        "falsifier": "either `with_predicted_target_only_a` or `with_predicted_target_only_b` "
        "changing. These two counts are conditioned on `predicted` being truthy, which is the "
        "yes/no, so by THE_INVARIANT they must NOT move; a move means `cached_prediction` and the "
        "window disagree about whether anything was predicted at all. The GENE may change and that "
        "is the expected gain.",
    },
    {
        "module": "genomeos/report.py",
        "lines": (98, 99, 100, 104),
        "reads": "regulation.py's rows -> e['predicted'] for the gene report's "
        "`strongest_predicted`, over the elements flagged `predicted_this_gene`",
        "claims_today": "the strongest predicted element FOR THIS GENE, selected from a set that "
        "regulation.py built on one gene per element",
        "class": "identity",
        "status": "free",
        "falsifier": "a gene whose `strongest_predicted` element changes identity. This module "
        "cannot be moved before regulation.py is, and regulation.py is frozen, so the registered "
        "expectation is that this falsifier CANNOT FIRE this wave: it is named so that a change "
        "appearing here without regulation.py moving would be a defect and not a result.",
    },
    {
        "module": "genomeos/benchmark/rearrangements.py",
        "lines": (427,),
        "reads": "loci._deletion_rows -> e['predicted_coding'] or e['predicted']",
        "claims_today": "`names_the_recipient`: whether the rearrangement's recipient gene is among "
        "the HEAD genes of the scored elements within 100 kb, plus `genes_named` truncated to 8",
        "class": "identity",
        "status": "free",
        "falsifier": "a rearrangement whose `names_the_recipient` goes False -> True. This is a "
        "yes/no about a NAMED gene rather than about any gene, so THE_INVARIANT does not protect "
        "it, and one flip establishes that the benchmark's negative was a property of the "
        "projection. `genes_named`'s truncation to 8 must be widened or declared before the count "
        "is read, because a longer list silently truncated would hide the gain.",
    },
    {
        "module": "genomeos/knowledge/across.py",
        "lines": (326,),
        "reads": "the committed vista_chr*.json rows -> x['predicted']['gene'], falling back to "
        "x['inferred']['gene']",
        "claims_today": "one `gene` per VISTA element in the across-species view, with a predicted "
        "head preferred over the inferred nearest TSS",
        "class": "identity",
        "status": "free",
        "falsifier": "an element whose across-species gene changes. This module reads COMMITTED "
        "vista results, so it cannot move until vista.py is re-run; the registered expectation is "
        "that this falsifier cannot fire on the stored files, and it is named so that the "
        "dependency is on the record rather than discovered later.",
    },
    {
        "module": "genomeos/cli.py",
        "lines": (3840, 3848, 4109),
        "not_head_lines": {
            3363: "a MISSENSE variant's predicted protein change",
            3365: "the same missense variant's predicted score",
            3366: "the same missense variant's predicted class",
            4229: "a contact-map comparison's count of MEASURED Hi-C boundaries",
            4276: "the same count, for insulation minima from the predicted contact map",
        },
        "reads": "a stored element row -> r['predicted'], r['predicted_coding'], e['predicted']",
        "claims_today": "the printed lines `predicted target: <gene>`, `strongest coding gene: "
        "<gene>` and the regulation table's predicted cell; `predicted target: none (no gene moves "
        "by <min> log2)` where the head is empty, which is correct under THE_INVARIANT",
        "class": "identity",
        "status": "peer_held",
        "falsifier": "a printed line that states one gene where the window has several, without "
        "saying how many. The falsifier is the wording, not a number: the `none` line is already "
        "right and must not be touched, and any new line must say the count or it is a wider claim "
        "than the reader has.",
    },
    {
        "module": "genomeos/attribution/confidence_calibration.py",
        "lines": (245, 891),
        "reads": "an element row -> e['predicted_coding'] for `gene` as a gate and for "
        "`confidence` / abs(`log2_fold_change`) as the number",
        "claims_today": "the confidence the compiler writes on this element's block and rule, "
        "recomputed identically; an element with no predicted target states nothing and is None",
        "class": "invariant",
        "status": "free",
        "falsifier": "any element's stated confidence changing under the window. Registered as "
        "EXPECTED NOT TO MOVE and as a module NOT to be moved: the gate is the yes/no and the "
        "number is the head's magnitude, which is the window maximum by `predict_target`'s own "
        "rule. If it moves, THE_INVARIANT is wrong and the whole wave stops; if it does not, this "
        "module is left exactly as it is and that is its result.",
    },
    {
        "module": "genomeos/response_map2.py",
        "lines": (413,),
        "reads": "an element row -> e['predicted_coding'] for `basis` only",
        "claims_today": "the interval's `basis` string; no gene symbol is read and no count is taken",
        "class": "invariant",
        "status": "free",
        "falsifier": "a `basis` string changing. Registered as EXPECTED NOT TO MOVE and as a module "
        "NOT to be moved: `basis` is a property of the element's record, not of which gene was "
        "chosen. A change here would mean the per-gene records carry different bases, which would "
        "be a finding about the sweep's own output and not about this consumer.",
    },
)

#: Modules that read a compact head without being a live consumer of one. Named so the census is
#: neither padded by them nor silently missing them.
RECORD_ONLY: tuple[dict[str, str], ...] = (
    {
        "module": "genomeos/predict/enhancer_target.py",
        "why": "it IS the projection: predict_target writes `predicted` and `predicted_coding` from "
        "the window it already holds. Also a closure file",
    },
    {
        "module": "genomeos/attribution/respmap_v2.py",
        "why": "it reads the head out of a PUBLISHED assertion's own quote, to audit what was "
        "published. Re-reading the window there would change the audit's subject",
    },
    {
        "module": "genomeos/attribution/silenceragree.py",
        "why": "it AUDITS the published head rather than consuming it: the first external check of "
        "`silencer_like`, whose two registered readings are named `reading_1_as_shipped` and "
        "`reading_2_cell_matched` and are both defined ON the shipped head. Re-reading the window "
        "would measure a quantity the project does not publish, so the agreement figure would "
        "answer no question that was asked -- the same reason as respmap_v2.py above. Of its three "
        "matching lines, 333 and 695 read the head's value to recompute the shipped call and its "
        'strata base rate; 244 reads only `r.get("predicted") is not None`, the PRESENCE of a '
        "head and never its content, as part of a row fingerprint",
    },
    {
        "module": "genomeos/benchmark/loci_noncoding.py",
        "why": "a census module that exists to document the one-key-per-element defect in "
        "loci.read_deletion; its head reads are the evidence, not a consumer",
    },
    {
        "module": "genomeos/benchmark/loci_reread.py",
        "why": "the same, for the re-read of the locus benchmark",
    },
    {
        "module": "genomeos/benchmark/loci_gene_input.py",
        "why": "the gene-input layer's own registered reader, written 2026-10-03 by another lane and "
        "held by it on the board; its head tuple is the declared v1 reading",
    },
)

#: Already moved, each under its own registration. This registration replaces none of them and
#: re-opens none of their figures.
ALREADY_MOVED: tuple[dict[str, str], ...] = (
    {"module": "genomeos/attribution/crispri.py", "when": "2026-09-22, lane-cache, bda8a83"},
    {"module": "genomeos/attribution/eqtl.py", "when": "2026-09-22, lane-reader2, 350c4c3"},
    {"module": "genomeos/attribution/target_calibration.py", "when": "2026-09-22, lane-calib, 5881d45"},
    {"module": "genomeos/attribution/motif_transfer.py", "when": "2026-09-27, lane-onetarget, d717b28"},
    {"module": "genomeos/attribution/syntax_tiling.py", "when": "2026-09-27, lane-onetarget, d717b28"},
    {"module": "genomeos/attribution/candidates.py", "when": "2026-09-27, lane-onetarget, d717b28"},
    {
        "module": "genomeos/benchmark/loci.py",
        "when": "2026-09-27, lane-onetarget, d717b28 (read_gene_input; its window rule was FALSIFIED "
        "and the one-target sum was KEPT, with window_reading carried beside it)",
    },
    {"module": "genomeos/attribution/unknown_scoring.py", "when": "wave 0, in progress"},
    {"module": "genomeos/attribution/joint_pretest.py", "when": "reads ElementResponses by default"},
)

#: Modules whose `predicted` spelling is NOT a compact element head. Excluded by name so the census
#: cannot be padded and so a reader can check the exclusion rather than take it on trust.
EXCLUDED_BY_NAME: tuple[dict[str, str], ...] = (
    {"module": "genomeos/attribution/n1_perturb_response.py", "why": "a pipeline STEP named 'predicted'"},
    {"module": "genomeos/genome/missense.py", "why": "a missense variant's predicted class and score"},
    {"module": "genomeos/genome/diff.py", "why": "the same missense prediction, in the variant diff"},
    {"module": "genomeos/genome/individuals.py", "why": "the same missense prediction, per individual"},
    {"module": "genomeos/therapeutics/structure.py", "why": "PREDICTED protein STRUCTURES (AlphaFold)"},
    {"module": "genomeos/attribution/targets.py", "why": "the reader itself"},
    {
        "module": "genomeos/attribution/onetarget2.py",
        "why": "this registration. Its matching lines are PROSE quoting each consumer's own read "
        "inside `reads` and `falsifier`; it reads no element row and opens no table at run time, "
        "and the drift guard found it on the first run, which is the guard working",
    },
    {"module": "genomeos/attribution/crispri_direction.py", "why": "its own stored rows, moved family"},
    {"module": "genomeos/attribution/crispri_direction_both.py", "why": "its own stored rows, moved family"},
)

#: `scripts/` entry points that read a compact head. OUT OF SCOPE for this registration: each writes
#: a committed result, and re-running one is its own decision with its own cost and its own
#: registration. Declared by path and counted so the figure cannot drift unnoticed.
SCRIPTS_OUT_OF_SCOPE = (
    "scripts/c4_ablation.py",
    "scripts/c5_paired_variation_probe.py",
    "scripts/consequence_targets.py",
    "scripts/constrained_targets.py",
    "scripts/direction_v2.py",
    "scripts/enhancer_targets.py",
    "scripts/enhancer_targets_all.py",
    "scripts/enhancer_targets_all_genome_wide.py",
    "scripts/loci_candidates.py",
    "scripts/loci_fourth.py",
    "scripts/loci_score.py",
    "scripts/loci_stated_intervals.py",
    "scripts/loci_third.py",
    "scripts/measured_layer.py",
    "scripts/node_containment_audit.py",
    "scripts/oriented_domains.py",
    "scripts/placement_audit.py",
    "scripts/repress_population.py",
    "scripts/repression_trace.py",
    "scripts/response_map_increment4_count.py",
    "scripts/satmut_vista.py",
    "scripts/shortlist_in_real_unknown.py",
    "scripts/variant_confidence_census.py",
    "scripts/vista_score.py",
)

# ---------------------------------------------------------------------------
# what this registration authorises, and what it refuses
# ---------------------------------------------------------------------------

AUTHORISES_NO_CONCLUSION = (
    "This file and the result file it writes state a census, a plan and a falsifier per module. "
    "They hold no measurement of any consumer's output, no verdict and no figure this wave will "
    "produce, and the fields a run would fill DO NOT EXIST in them rather than sitting empty, so "
    "that nothing here can read as a preview. The only figures present are the SCOPE -- how many "
    "modules read a head and at which lines, which `check_census_matches_tree()` re-derives from "
    "the tree -- and PRIOR EVIDENCE quoted from a committed result with its date and commit. It "
    "authorises no conclusion of any kind, and in particular it does not authorise calling any "
    "consumer's current figure wrong."
)

CANNOT_ESTABLISH = (
    "that the window's extra genes are CORRECT: they are the same model's output read without the "
    "projection, not a measurement, and nothing here makes them evidence",
    "that a consumer whose figure moves was wrong: a figure computed on the head is a true figure "
    "about the head, and the wave's finding is about what it is a figure OF",
    "that a consumer whose figure does not move is free of the defect: it may carry no element with "
    "a second gene at the bar in its own population, which is a fact about that population",
    "anything about the node model, area J, area E's fate ceiling, BioForge's confidence or any "
    "biosample: this wave touches none of them and the stop list is not approached",
    "anything requiring an AlphaGenome request: the window is on disk and 0 requests are spent",
)

REFUSALS = (
    "check_census_matches_tree(): refuses if a declared module's declared lines are not the lines "
    "the tree holds, if the tree holds a head-reading module under genomeos/ that is in no list, or "
    "if a declared module is not in the tree at all",
    "check_one_falsifier_per_module(): refuses if any census module has no falsifier, or if two "
    "modules share one verbatim -- a falsifier shared between modules is the defect this row names",
    "check_no_closure_file_is_moving(): refuses if a closure file is given any status but "
    "frozen_closure, so the wave cannot reach the frozen paid study's import closure",
    "check_no_preview(payload): refuses if the registration payload holds a field a run would fill, "
    "or a float, or any of the forbidden result spellings outside the registry's own name stamp",
    "check_scope_counts(): refuses if any figure quoted in this registration's own prose is not "
    "what CENSUS holds, so the registration cannot misstate its own scope",
    "check_the_quoted_invariant(): refuses unless data/results/onetarget_consumers.json really "
    "holds the 849,469 compared and 0 disagreements THE_INVARIANT quotes, and unless THE_INVARIANT "
    "states that figure. The one piece of prior evidence this wave rests on is checked against the "
    "file, not recalled",
)

#: Refusals the MOVES added, kept apart from REFUSALS so that the registered tuple stays exactly as
#: 24adf33 registered it and the test asserting its length keeps passing unedited. Superseding is
#: additive: a registration's own list is not rewritten because later work needed more refusals.
REFUSALS_ADDED_BY_THE_MOVES = (
    "check_registered_lines_are_not_edited(): refuses any change to the lines 24adf33 registered; "
    "a moved module's new lines go in MOVED_LINES additively, so the pre-move claim is kept",
    "check_head_invariant(census, where): refuses on ANY any-gene head disagreement in any "
    "module's move, because that refutes the ground the whole wave stands on. One refusal for one "
    "invariant, never a substitute for a module's own falsifier; the CODING head is not checked "
    "here, by the registered limit",
)

#: Field names a run would fill. None of them may exist in the registration payload.
FORBIDDEN_FIELDS = (
    "verdict",
    "verdicts",
    "moved",
    "unmoved",
    "disagreements",
    "head_disagreements",
    "compared",
    "gained",
    "lost",
    "figures",
    "outcome",
    "conclusion",
    "fired",
    "passed",
    "failed",
)

WILL_WRITE = "data/results/onetarget2.json"

#: The committed result THE_INVARIANT is quoted from, and the two figures the quote states. The
#: quote is checked against the file rather than taken on trust: a registration that misquotes its
#: own prior evidence is worse than one that cites none.
PRIOR_EVIDENCE = "data/results/onetarget_consumers.json"
PRIOR_EVIDENCE_ARM = "candidates"
PRIOR_EVIDENCE_COMPARED = 849469
PRIOR_EVIDENCE_DISAGREEMENTS = 0

OWN_CODE = (
    "genomeos/attribution/onetarget2.py",
    "scripts/onetarget2_register.py",
    "tests/test_onetarget2.py",
)

_ALT = "|".join(HEAD_KEYS)
_SITE = re.compile(rf"""(?:get\(\s*["'](?:{_ALT})["']|\[\s*["'](?:{_ALT})["']\s*\])""")


#: The lines the registration declared AS COMMITTED at 24adf33, frozen. A move changes the tree, so
#: MOVED_LINES below says where each moved module's head reads sit now; this dict says where they
#: sat when the claim was registered, and `check_registered_lines_are_not_edited()` refuses any
#: change to it. Superseding is additive: the pre-move record is kept, never overwritten.
REGISTERED_LINES: dict[str, tuple[int, ...]] = {
    "genomeos/attribution/closure.py": (73,),
    "genomeos/attribution/compile.py": (332, 838, 880),
    "genomeos/attribution/organise.py": (164,),
    "genomeos/attribution/variation.py": (197,),
    "genomeos/attribution/measured.py": (1092,),
    "genomeos/attribution/not_open_profile.py": (254,),
    "genomeos/attribution/gwas.py": (153, 154, 155, 164, 166),
    "genomeos/attribution/argmaxcell.py": (677, 692),
    "genomeos/attribution/context_evidence.py": (441,),
    "genomeos/attribution/vista.py": (260, 277, 278, 280, 312, 326),
    "genomeos/decompile.py": (124,),
    "genomeos/genome/motifs.py": (522,),
    "genomeos/genome/regulation.py": (163,),
    "genomeos/genome/regdiff.py": (101, 176, 177, 194, 195),
    "genomeos/report.py": (98, 99, 100, 104),
    "genomeos/benchmark/rearrangements.py": (427,),
    "genomeos/knowledge/across.py": (326,),
    "genomeos/cli.py": (3840, 3848, 4109),
    "genomeos/attribution/confidence_calibration.py": (245, 891),
    "genomeos/response_map2.py": (413,),
}

#: Per module moved in this wave: where its compact-head reads sit AFTER the move, and the commit.
#: A moved module keeps its compact reading -- that is the point of the move -- so its head reads do
#: not disappear; they move, and usually one is ADDED where the window reading keeps the compact
#: fallback for an element the cache does not hold.
MOVED_LINES: dict[str, dict[str, Any]] = {
    "genomeos/attribution/closure.py": {
        "now": (92, 150),
        "was": (73,),
        "what": "attributed_elements is unchanged and still the reading every committed closure "
        "figure rests on (its head read moved 73 -> 92 on the docstring this move added). The "
        "second read, 150, is window_elements' COMPACT FALLBACK: the head it keeps for an element "
        "the response cache does not hold, and for one it holds without a coding window",
    },
    "genomeos/attribution/measured.py": {
        "now": (1141,),
        "was": (1092,),
        "what": "rows() keeps its one `predicted_gene` per row and every committed field; its head "
        "read moved 1092 -> 1141 behind the new window_agreement(). With a reader each row gains "
        "ONE new key, `window`, and nothing else moves",
    },
    "genomeos/attribution/not_open_profile.py": {
        "now": (262, 434),
        "was": (254,),
        "what": "rules() is NOT CHANGED AT ALL -- the diff against its committed form removes zero "
        "lines -- which is what keeps tests/test_not_open_profile.py's element-for-element pin to "
        "rule_loci true by construction; its head read moved 254 -> 262 only because the new source "
        "constants sit above it. The second read, 434, is the new window_rules(), which CALLS "
        "rules() and appends SOURCE_PREDICTED_WINDOW rules: a population in no committed program, "
        "kept out of SOURCES so no caller counting compiled rules moves",
    },
    "genomeos/attribution/variation.py": {
        "now": (298,),
        "was": (197,),
        "what": "elements_of() is NOT CHANGED AT ALL -- the diff against its committed form removes "
        "zero lines -- so every committed variation_chr*.json figure is produced by the unchanged "
        "function; its one read moved 197 -> 298 only because window_elements_of() and the two "
        "constants sit above it. window_elements_of() CALLS it and does not read a head itself: it "
        "reads the ROW elements_of already built, compares the stored magnitude with the window's "
        "value for the same gene, and RAISES if they differ, because the registration says a move "
        "there refutes the invariant and stops the wave rather than becoming a finding",
    },
    "genomeos/attribution/vista.py": {
        "now": (260, 277, 278, 280, 358, 438, 452),
        "was": (260, 277, 278, 280, 312, 326),
        "what": "score() is NOT TOUCHED: 260 (the over-length skip writing a None head), 277, 278 "
        "and 280 (predicted_group and tissue_agrees) are all at their registered lines, because the "
        "window arm sits BELOW score() and above summarise(). summarise()'s two reads moved 312, "
        "326 -> 430, 444 and are unchanged, so with_predicted_target and strong do not move. The "
        "NEW read, 350, is inside window_tissue_agreement(), which returns {} without a reader; it "
        "reads the same head in order to find the rows whose head DISAGREED, which is the "
        "population the registered falsifier asks about. tissue_agrees is never recomputed and "
        "never overwritten",
    },
    "genomeos/benchmark/rearrangements.py": {
        "now": (450,),
        "was": (427,),
        "what": "the ONE registered read is unchanged and still reads the same two heads in the same "
        "order; it moved 427 -> 450 only because the two constants and the widened signature sit "
        "above it. names_the_recipient is NOT touched, so no committed benchmark figure moves: the "
        "falsifier is read from names_the_recipient_in_the_window beside it, and the registration's "
        "demand about genes_named's truncation to 8 is met by DECLARING it (genes_named_total, "
        "genes_named_truncated_to) rather than widening a committed field",
    },
    "genomeos/attribution/gwas.py": {
        "now": (203, 232, 282, 283, 284, 293, 295),
        "was": (153, 154, 155, 164, 166),
        "what": "summarise() keeps every committed figure and every registered read: the three "
        "partition lines moved 153, 154, 155 -> 274, 275, 276 and the agreement selection 164, 166 "
        "-> 285, 287, pushed down by the window arm added above them and not rewritten. The two NEW "
        "reads, 203 and 224, are inside window_agreement(), which returns {} without a reader so "
        "summarise contributes no new key at all; 203 is the head half of the comparison the "
        "registered falsifier names and 224 names the head beside the gene that was gained. The "
        "catalogue's mapped-gene rule is read through a NEW helper, _mapped_genes, so the two "
        "committed expressions that produce the stored shares are left exactly as they are",
    },
    "genomeos/genome/regdiff.py": {
        "now": (128, 205, 220, 295, 296, 314, 315),
        "was": (101, 176, 177, 194, 195),
        "what": "every registered read is still there and still reads the same head: the rank term "
        "moved 101 -> 128, the two with_predicted_target counts 176, 177 -> 295, 296 and the printed "
        "line 194, 195 -> 314, 315, all of them pushed down by the constants and the window arm "
        "added above them and none of them rewritten. The two NEW reads, 205 and 220, are inside "
        "that arm and read the SAME head to compare it with the window: 205 is the head half of the "
        "count the registered falsifier names, and 220 names the head beside the genes the "
        "projection dropped. `variants_in_elements` still gets its one `predicted` gene from "
        "`cached_prediction` unchanged, and the window is a new key beside it that exists only when "
        "a reader was asked for",
    },
}


#: THE TWO FREE MODULES THIS WAVE DELIBERATELY DID NOT MOVE, with the dependency that stops them
#: CHECKED against the record rather than recalled. Both are class `identity` and status `free` in
#: the census -- nothing holds them on the work board and neither is in the frozen closure -- and
#: both are registered with a falsifier that CANNOT FIRE this wave. A change appearing in either
#: without its upstream moving would therefore be a DEFECT and not a result, which is exactly why
#: the registration named them: the dependency is on the record, not discovered later.
NOT_MOVED_THIS_WAVE: tuple[dict[str, Any], ...] = (
    {
        "module": "genomeos/report.py",
        "falsifier": "a gene whose `strongest_predicted` element changes identity",
        "cannot_fire_because": "its rows come from genomeos/genome/regulation.py, which IS one of "
        "the 52 files of the frozen import closure and is not this wave's to touch. report.py "
        "itself is NOT in the closure, so the block is the upstream and not the module.",
        "checked": "report.py imports genomeos.genome.regulation; "
        "genomeos/genome/regulation.py in sender_closure() is True; "
        "genomeos/report.py in sender_closure() is False.",
        "what_would_be_a_defect": "any change to strongest_predicted while regulation.py is "
        "unchanged, because the set it selects from is built one gene per element upstream",
    },
    {
        "module": "genomeos/knowledge/across.py",
        "falsifier": "an element whose across-species gene changes",
        "cannot_fire_because": "it reads the COMMITTED data/results/vista_chr*.json files, so it "
        "cannot move until vista.py is RE-RUN and those files are rewritten. The registration says "
        "this in those words, and it is NOT the same dependency as report.py's: across.py is not "
        "downstream of regulation.py at all.",
        "checked": "across.py reads data/results/vista_chr*.json; 23 such files on disk, all "
        "tracked by git and none modified in the working tree, so the gene it reads today is the "
        "committed one.",
        "what_would_be_a_defect": "any change to the across-species gene without a new vista run, "
        "because the stored head is the only thing this module reads",
        "what_this_wave_DID_establish_about_it": "vista.py's own window arm fired on 13 rows of "
        "2,223, so the stored vista heads ARE narrower than the window at those elements. That is a "
        "measurement about vista.py's rows and NOT a change to across.py's output, which still "
        "reads the committed head and is untouched by this wave.",
    },
)


#: AMENDMENTS TO THE REGISTRATION, additive and each with the module that forced it. The rule
#: 24adf33 set is that a module which cannot move without widening its claim says so HERE rather
#: than widening it in silence. Nothing here weakens a falsifier or a refusal; each entry adds a
#: quantity the shared reader did not provide and names what asserts it did not become a second rule.
AMENDMENTS: tuple[dict[str, Any], ...] = (
    {
        "amendment": "the shared reader gains the TRACK beside the gene at the bar",
        "added": "attribution/targets.tracks_at_bar",
        "forced_by": "genomeos/attribution/vista.py",
        "why": "vista's registered falsifier is about tissue_agrees, which compares the predicted "
        "track GROUP with the element's observed expression groups. WindowReading carries the gene "
        "and the signed change and drops the track, so a second gene at the bar has no group to "
        "compare and the falsifier cannot be asked at all without it.",
        "what_is_not_widened": "no falsifier, no refusal and no census line changes; the compact "
        "head stays the head and tissue_agrees itself is untouched. The answer read off the window "
        "is a NEW key beside it.",
        "how_it_is_kept_from_becoming_a_second_size_rule": "tracks_at_bar applies genes_at_bar's own "
        "rule to the same fields in the same order, and tests/test_vista_window.py asserts that its "
        "(gene, signed) projection EQUALS genes_at_bar's output on fixtures and on real chromosome "
        "data, so the two cannot drift apart unnoticed.",
    },
)


class CensusDriftError(RuntimeError):
    """The declared census is not what the tree holds."""


class SharedFalsifierError(RuntimeError):
    """Two modules were given one falsifier, or a module was given none."""


class ClosureFileMovingError(RuntimeError):
    """A file of the paid study's import closure was marked movable."""


class PreviewInRegistrationError(RuntimeError):
    """The registration holds something a run would fill."""


def _root() -> Path:
    return Path(__file__).resolve().parents[2]


def head_read_sites(root: Path | None = None, package: str = "genomeos") -> dict[str, list[int]]:
    """Every file under `package` with a compact-head read, and the 1-based lines it reads one on.

    Mechanical, so the census cannot drift from the tree: the declared lines are checked against
    this and nothing else.
    """
    base = (root or _root()) / package
    out: dict[str, list[int]] = {}
    for p in sorted(base.rglob("*.py")):
        lines = [i + 1 for i, line in enumerate(p.read_text().splitlines()) if _SITE.search(line)]
        if lines:
            out[p.relative_to(root or _root()).as_posix()] = lines
    return out


def declared() -> dict[str, tuple[int, ...]]:
    return {c["module"]: tuple(c["lines"]) for c in CENSUS}


def effective_lines() -> dict[str, tuple[int, ...]]:
    """Per module, the lines the TREE should match today: `MOVED_LINES['now']` once it has moved,
    the registered lines until then. One definition, shared by the guard and its tests, so the two
    cannot drift and a passing plant cannot be a plant against a stale rule."""
    return {
        c["module"]: tuple(MOVED_LINES[c["module"]]["now"])
        if c["module"] in MOVED_LINES
        else tuple(c["lines"])
        for c in CENSUS
    }


def not_head() -> dict[str, dict[int, str]]:
    """Per module, the matching lines that are NOT a compact element head, each with its reason."""
    return {c["module"]: dict(c.get("not_head_lines") or {}) for c in CENSUS}


def check_census_matches_tree(root: Path | None = None) -> dict[str, list[int]]:
    """Refuse unless the declared census is exactly what the tree holds."""
    found = head_read_sites(root)
    dec = declared()
    accounted = (
        set(dec)
        | {m["module"] for m in RECORD_ONLY}
        | {m["module"] for m in ALREADY_MOVED}
        | {m["module"] for m in EXCLUDED_BY_NAME}
    )
    unknown = sorted(set(found) - accounted)
    if unknown:
        raise CensusDriftError(
            f"{len(unknown)} head-reading module(s) under genomeos/ are in no list: {unknown}. "
            "Classify each one or name it in EXCLUDED_BY_NAME with its reason."
        )
    eff = effective_lines()
    for mod, lines in sorted(dec.items()):
        lines = eff[mod]
        if mod not in found:
            raise CensusDriftError(f"{mod} is declared in the census and reads no compact head in the tree")
        other = set(not_head().get(mod, {}))
        if set(lines) & other:
            raise CensusDriftError(
                f"{mod}: line(s) {sorted(set(lines) & other)} are declared both as a head read and as not one"
            )
        if set(found[mod]) != set(lines) | other:
            raise CensusDriftError(
                f"{mod}: the census accounts for lines {sorted(set(lines) | other)} and the tree "
                f"matches on {found[mod]}. Re-read the module before touching it: a line that "
                "moved may be a peer's edit, and a line that vanished may be a claim that was "
                "dropped. A line that is not a compact head goes in `not_head_lines` with its "
                "reason, never left out."
            )
    return found


def check_registered_lines_are_not_edited() -> None:
    """Refuse if CENSUS's declared lines have been changed from what 24adf33 registered.

    A module that moves changes the tree, and the temptation is to update the census to match. That
    would quietly replace the claim this wave was registered against. MOVED_LINES carries the new
    lines additively instead, and this refusal keeps the registered ones where they are.
    """
    if declared() != REGISTERED_LINES:
        bad = {
            m: (REGISTERED_LINES.get(m), declared().get(m))
            for m in set(REGISTERED_LINES) | set(declared())
            if REGISTERED_LINES.get(m) != declared().get(m)
        }
        raise CensusDriftError(
            f"CENSUS's declared lines no longer match what 24adf33 registered: {bad} (registered, "
            "now). A moved module's new lines go in MOVED_LINES; the registered ones do not change."
        )
    for mod, rec in MOVED_LINES.items():
        if tuple(rec.get("was") or ()) != REGISTERED_LINES.get(mod):
            raise CensusDriftError(
                f"MOVED_LINES[{mod!r}]['was'] is {rec.get('was')} and the registration said "
                f"{REGISTERED_LINES.get(mod)}. A move record may not restate the claim."
            )


def check_one_falsifier_per_module() -> None:
    """Refuse unless every census module has a falsifier of its own, no two alike."""
    seen: dict[str, str] = {}
    for c in CENSUS:
        f = (c.get("falsifier") or "").strip()
        if not f:
            raise SharedFalsifierError(f"{c['module']} has no falsifier; the row asks for one per module")
        if f in seen:
            raise SharedFalsifierError(
                f"{c['module']} and {seen[f]} were given one falsifier verbatim. A falsifier shared "
                "between modules lets one of them pass on the other's evidence, which is the "
                "defect this row exists to avoid."
            )
        seen[f] = c["module"]


def check_no_closure_file_is_moving() -> None:
    """Refuse unless every closure file in the census is held at `frozen_closure`."""
    for c in CENSUS:
        if c["module"] in CLOSURE_FROZEN and c["status"] != "frozen_closure":
            raise ClosureFileMovingError(
                f"{c['module']} is in the paid study's 52-file import closure and is marked "
                f"{c['status']!r}. sender_closure() hashes the WORKING TREE, so even an uncommitted "
                "edit here refuses the frozen send. This is Albert's call, not a lane's."
            )


def check_no_preview(payload: dict[str, Any]) -> None:
    """Refuse unless the payload holds no field a run would fill and no measurement of its own."""

    def walk(node: Any, path: str) -> None:
        if isinstance(node, dict):
            for k, v in node.items():
                if k in FORBIDDEN_FIELDS and path:
                    raise PreviewInRegistrationError(
                        f"{path}.{k} is a field a run would fill; a registration may not hold it"
                    )
                walk(v, f"{path}.{k}" if path else k)
        elif isinstance(node, list | tuple):
            for i, v in enumerate(node):
                walk(v, f"{path}[{i}]")
        elif isinstance(node, bool):
            return
        elif isinstance(node, float):
            raise PreviewInRegistrationError(f"{path} is a float; no measurement may be in a registration")

    walk(payload, "")


def check_scope_counts() -> None:
    """Refuse unless every declared scope figure is what CENSUS and the lists actually hold."""
    got = {
        "row_says": 5,
        "census_modules": len(CENSUS),
        "identity": len(by_class("identity")),
        "invariant": len(by_class("invariant")),
        "record_only": len(RECORD_ONLY),
        "already_moved": len(ALREADY_MOVED),
        "excluded_by_name": len(EXCLUDED_BY_NAME),
        "identity_free_today": len(free_to_move()),
        "identity_peer_held": len(
            [c for c in CENSUS if c["class"] == "identity" and c["status"] == "peer_held"]
        ),
        "identity_frozen_closure": len(
            [c for c in CENSUS if c["class"] == "identity" and c["status"] == "frozen_closure"]
        ),
        "closure_files_reading_a_head": len(CLOSURE_FROZEN),
        "scripts_out_of_scope": len(SCRIPTS_OUT_OF_SCOPE),
    }
    if got != SCOPE_COUNTS:
        bad = {k: (SCOPE_COUNTS.get(k), v) for k, v in got.items() if SCOPE_COUNTS.get(k) != v}
        raise CensusDriftError(
            f"declared scope counts do not match the census: {bad} (declared, actual). The prose in "
            "ROW_SAYS_FIVE_THE_TREE_SAYS_OTHERWISE quotes these figures, so a drift here is a "
            "registration that misstates its own scope."
        )


class MisquotedPriorEvidenceError(RuntimeError):
    """THE_INVARIANT does not say what the committed result it cites says."""


def check_the_quoted_invariant(path: Path | None = None) -> dict[str, Any]:
    """Refuse unless the committed result really holds the figures THE_INVARIANT quotes."""
    import json

    p = path or (_root() / PRIOR_EVIDENCE)
    if not p.exists():
        raise MisquotedPriorEvidenceError(
            f"{PRIOR_EVIDENCE} is not on this machine; THE_INVARIANT may not be quoted without it"
        )
    hc = (json.loads(p.read_text()).get(PRIOR_EVIDENCE_ARM) or {}).get("head_control") or {}
    want = {
        "compared": PRIOR_EVIDENCE_COMPARED,
        "head_disagreements": PRIOR_EVIDENCE_DISAGREEMENTS,
        "coding_head_disagreements": PRIOR_EVIDENCE_DISAGREEMENTS,
    }
    got = {k: hc.get(k) for k in want}
    if got != want:
        raise MisquotedPriorEvidenceError(
            f"{PRIOR_EVIDENCE}[{PRIOR_EVIDENCE_ARM}].head_control holds {got} and THE_INVARIANT "
            f"quotes {want}. Re-read the committed result before quoting it."
        )
    if f"{PRIOR_EVIDENCE_COMPARED:,}" not in THE_INVARIANT:
        raise MisquotedPriorEvidenceError(
            f"THE_INVARIANT does not state the {PRIOR_EVIDENCE_COMPARED:,} elements it rests on"
        )
    return dict(hc)


class InvariantRefutedError(RuntimeError):
    """The window's at-the-bar head is not the compact table's gene. The wave stops here."""


def check_head_invariant(census: dict[str, Any], where: str) -> None:
    """Refuse on any ANY-GENE head disagreement. ONE refusal, because it is ONE invariant.

    This is deliberately not a per-module falsifier and is not one of the eighteen. The eighteen
    falsifiers each ask what a module's own figure does; this asks whether the ground the whole wave
    stands on is still there. A single disagreement means `predict_target` and `genes_at_bar` no
    longer apply the same size rule with the same tie order, every module's move is void, and the
    wave stops rather than reporting a result.

    The CODING head is explicitly NOT checked here, and that is the registered limit doing work:
    the 2026-09-27 result found 10 elements of 4,794 with no `predicted_coding` field although the
    cache named a coding gene, so a coding-head disagreement is expected, is counted by name, and
    is not a refutation.
    """
    n = int(census.get("head_disagreements") or census.get("head_disagrees") or 0)
    if n:
        raise InvariantRefutedError(
            f"{where}: {n} element(s) where the window's at-the-bar head is not the compact "
            f"table's `predicted` gene. THE WAVE STOPS. The registered invariant "
            f"({PRIOR_EVIDENCE_COMPARED:,} elements, {PRIOR_EVIDENCE_DISAGREEMENTS} disagreements, "
            "2026-09-27) is refuted on this population, so no module's move may be read as a "
            "result. Report it before anything else."
        )


def check_all(root: Path | None = None, payload: dict[str, Any] | None = None) -> None:
    check_registered_lines_are_not_edited()
    check_census_matches_tree(root)
    check_scope_counts()
    check_one_falsifier_per_module()
    check_no_closure_file_is_moving()
    check_the_quoted_invariant()
    if payload is not None:
        check_no_preview(payload)


def by_class(klass: str) -> tuple[str, ...]:
    return tuple(c["module"] for c in CENSUS if c["class"] == klass)


def by_status(status: str) -> tuple[str, ...]:
    return tuple(c["module"] for c in CENSUS if c["status"] == status)


def free_to_move() -> tuple[str, ...]:
    """The identity-class consumers this wave may actually touch today."""
    return tuple(c["module"] for c in CENSUS if c["class"] == "identity" and c["status"] == "free")
