# SPDX-License-Identifier: AGPL-3.0-or-later
"""Does a block execute the value stored in it? The shortlist, the controls and the pre-registered test.

Albert's third class: a program is a unit whose read-out changes with the value held in the storage
beside it. The human panel (attribution/human_panel.py) catalogued storage units, 200-bp windows
where the 89 assemblies carry a few recurring values; some of those values are measured to change
a read-out, as a fine-mapped GTEx eQTL (the allele goes with a gene's expression across donors) or
as an MPRA allele pair (the allele changes a reporter's activity). This module turns that
correlation into a test that can fail:

1. **The shortlist** (no model call): every storage unit on the read chromosomes whose recurring
   value *is* a measured variant, with the alleles and their panel counts, the gnomAD frequency and
   allele number, the GTEx slope and p-value per tissue for the gene (streamed from GTEx v8 by
   `attribution/eqtl.py`'s reader into this lane's own directory), the MPRA log2FC and FDR per cell,
   the block and tier it sits in, and its GC and replication-timing strata.
2. **Matched controls** (no model call): for each shortlist unit, storage units nearby (same
   chromosome, within 250 kb, so the same genes are in reach of a 1 Mb window) with the same GC
   stratum, the same timing tertile and the same value structure (two, three, or four or more
   recurring values; substitutions), no significant GTEx pair and no MPRA row at any of their
   events, and a panel r-squared under 0.2 with the tested variant. A control borrows the unit's
   gene, tissue or cell and measured sign: it asks how often the model's allele effect on that
   gene has that sign when the variant is not a measured regulator.
3. **The criterion** (`CRITERION`), written before any request: an allele-dependent predicted effect
   that agrees in sign with the measured direction more often than in the matched controls, with
   the thresholds, the alpha spent at each interim look, the order of the requests and the stopping
   rule all fixed in code and in the committed result.
4. **The run** (`run_pairs`), which only spends requests once the quota is handed over: the existing
   variant path (`predict/alphagenome_adapter.py`'s RNA-seq gene scorer, both alleles in one request
   over a 1 Mb window, cached through `predict/individual_effects.score_variant`), one request per
   unit and one per control, in pairs, so a stop at any point leaves balanced arms.

The model was trained on ENCODE tracks, and GTEx expression is close kin to them, so agreement with
an eQTL is partly the model reading back what it learned; the MPRA allele pairs are the cleaner
test and are scored and reported apart.

5. **The widened E1** (`E1_WIDE`, `wide_pairs`, `--wide`), pre-registered 2026-09-14 after E2 and E3:
   the storage-unit framing is dropped and every significant MPRAVarDB allele pair genome-wide is
   tested, in the cell lines the model has RNA-seq tracks for, against a control the same study
   measured in the same cell and found not significant. The read-out gene is chosen without the
   model (DAP-G at the variant, else the nearest protein-coding TSS), because this is the endpoint
   whose outcome is an external measurement rather than expression data of the kind the model was
   trained on.
"""

from __future__ import annotations

import functools
import gzip
import json
import math
import random
import re
import time
from bisect import bisect_left
from collections import Counter
from collections.abc import Callable
from pathlib import Path
from typing import Any

from genomeos.attribution.human_panel import (
    CACHE as PANEL_CACHE,
)
from genomeos.attribution.human_panel import (
    MIN_RECURRING,
    Overlaps,
    Panel,
    events,
    gc_stratum,
    track_rows,
)

LANE = PANEL_CACHE / "executor"
GTEX_DIR = LANE / "gtex"
PREDICTIONS = LANE / "predictions"
CONTROL_WINDOW = 250_000  # controls sit within this distance of the unit's tested variant
LD_MAX_R2 = 0.2  # a control must not be in panel linkage with the tested variant
CONTROLS_PER_UNIT = 2  # raised from 1 on 2026-09-14 after the first 308 requests: see AMENDMENTS
MATCH_LEVELS = (
    ("gc_stratum", "rt_stratum", "value_bucket"),  # 0: every stratum
    ("gc_stratum", "rt_stratum"),  # 1: the number of values relaxed
    ("gc_stratum",),  # 2: timing relaxed as well; never GC
)
MPRA_FDR = 0.05
MPRA_P = 0.01  # a row without an FDR counts as measured at this nominal p
MPRA_EXCLUDE = "3'UTR"  # those libraries measure mRNA stability, not transcription: reported apart
EQTL_P = 1e-5
DAPG_PIP = 0.5  # fine-mapped: the primary eQTL endpoint, where the value is probably the cause
ENDPOINTS = ("E1_mpra", "E2_eqtl", "E3_eqtl_linked")
SEED = 20260914

# -- the widened E1: every significant MPRAVarDB allele pair, genome-wide ------------------------
MPRA_URL_KEY = "mpravardb"  # human_panel.BIGBEDS: gbdb/hg38/mpra/mpravardb/mpravardb.bb
MPRA_CELLS = ("GM12878", "HepG2", "K562", "Jurkat", "HeLa")  # MPRA lines the model has tracks for
MPRA_NULL_P = 0.2  # a control's own measurement must be this far from significant
WIDE_MIN_DISTANCE = 200  # a control outside the tested variant's own 200-bp window
WIDE_MAX_DISTANCE = 250_000  # and inside the 1 Mb prediction window, as in E2 and E3
SATURATION = "Saturation mutagenesis"  # the one study that mutates whole elements: a stratum apart
WIDE_ENDPOINTS = ("E1W_mpra_wide", "E1S_mpra_saturation")
MAIN_CHROMS = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y"])

AMENDMENTS: tuple[dict[str, str], ...] = (
    {
        "when": "2026-09-14, after a two-request smoke test and before any endpoint was scored",
        "what": "a no-call grid of 0, 0.001 and 0.01 reported beside the primary 0.001",
        "why": "single-variant gene effects came out near 0.001 log2, so the threshold choice is shown",
    },
    {
        "when": "2026-09-14, after the first 308 requests (155 pairs, E1 and E2)",
        "what": (
            "E1's controls are drawn within 250 kb like E2's rather than chromosome-wide; the cache filter "
            "keeps every gene's tracks in the read-out cell; two controls per unit instead of one"
        ),
        "why": (
            "19 of E1's 20 controls returned no value at all: a chromosome-wide control's borrowed gene was "
            "outside its window and the cache filter had dropped the cell's other genes, so the arm could "
            "not be scored. The unit arm is unchanged and its answers come from cache; only controls are "
            "added. Both readings are reported, before and after"
        ),
    },
    {
        "when": "2026-09-14, after E1 and E2 were scored and before any E3 row was read",
        "what": (
            "a unit with two controls is counted once in the unit arm and its controls once each, "
            "instead of once per pair"
        ),
        "why": (
            "the second amendment gave each unit two controls, and the first tally counted the unit's own "
            "answer twice: E2's unit arm read 68 of 100 where 56 units were answered. The difference barely "
            "moves (+0.269 to +0.268) but the p-value it deserves does, from 0.00036 to 0.0021, still a "
            "success by the criterion; every number from here counts each unit once and each control once"
        ),
    },
)

HOLD_OUT = {
    "when": "2026-09-14, written after E2's amended reading and before any E3 request",
    "question": (
        "is E2's agreement about the value being the cause, or about the model's behaviour in eQTL-rich "
        "neighbourhoods? E3's variants are significant eQTLs that DAP-G does not fine-map, so most are "
        "linked to a cause rather than causal"
    ),
    "prediction": (
        "if a unit executes its value, E3's agreement difference should be clearly smaller than E2's "
        "(dilution by linkage); if E3 matches E2, the difference is not about causality and E2's reading "
        "is a property of the neighbourhood, which would be the more important result"
    ),
    "budget": "2,000 requests, 1,000 pairs in the pre-registered order, stopping rule unchanged",
    "read": (
        "written 2026-09-14 by the session that ran E3, before any E3 row was read, because 'clearly "
        "smaller' needed a number: compare_endpoints takes E2's agreement difference minus E3's with a "
        "normal-approximation one-sided 95% bound. 'diluted by linkage' if that bound puts the gap above "
        "zero; 'E3 matches E2' if the gap does not clear zero while E3 on its own meets the success "
        "criterion (difference at or above 0.10 at p 0.01), which is the reading that says the effect is "
        "the neighbourhood and not the value; 'neither separated' otherwise"
    ),
}

CRITERION: dict[str, Any] = {
    "written": "2026-09-14, before any model request for this test",
    "claim_tested": (
        "a storage unit executes its value: swapping the stored allele changes the predicted read-out of "
        "the measured gene or cell in the measured direction, more often than the same swap at a matched "
        "unit that stores a value nobody has measured"
    ),
    "unit_of_test": (
        "one variant per shortlist unit: the recurring value event that is the measured variant (rsID "
        "match to gnomAD at the event, same alternative base); units whose measured variant is not one "
        "of their recurring values are listed but not tested"
    ),
    "prediction": (
        "AlphaGenome RNA-seq gene scorer through predict/alphagenome_adapter.py: both alleles in one "
        "request over the 1 Mb window centred on the variant; the read-out is the log2 fold change of "
        "alternative over reference on the implicated gene"
    ),
    "readout_eqtl": (
        "the gene GTEx names; the track whose GTEx tissue name matches the eQTL's most significant tissue, "
        "else the mean over the gene's tracks (declared per unit)"
    ),
    "readout_mpra": (
        "the MPRA's cell line tracks on the gene the variant's eQTL names; without an eQTL, the gene with "
        "the largest absolute predicted effect in that cell's tracks, chosen identically for the control"
    ),
    "measured_sign_eqtl": "sign of GTEx v8 slope for the alternative allele (hg38 ref/alt of the variant id)",
    "measured_sign_mpra": (
        "sign of MPRAVarDB log2FC, which UCSC's track description defines as log2(alt RNA/DNA) - "
        "log2(ref RNA/DNA): alternative over reference"
    ),
    "no_call": "a predicted |log2FC| below 0.001 is a no-call in either arm and is counted, not scored",
    "endpoints": {
        "E1_mpra": (
            "sign agreement among units whose value has an MPRA allele pair at FDR <= 0.05 (or nominal "
            "p <= 0.01 where no FDR is given), 3'UTR stability libraries excluded, against matched controls "
            "drawn chromosome-wide; the cleaner test, reported apart"
        ),
        "E2_eqtl": (
            "sign agreement among units whose value is DAP-G fine-mapped (PIP >= 0.5) for a gene with "
            "GTEx p <= 1e-5 and one sign across tissues, against controls within 250 kb; partly the model "
            "reading back ENCODE-like training data"
        ),
        "E3_eqtl_linked": (
            "the other significant eQTL values, which may only be linked to the cause; exploratory, run only "
            "after E1 and E2 are decided and only if the quota holder agrees, never a claim"
        ),
    },
    "controls": (
        "one per unit, never reused: same chromosome and GC stratum, and where the neighbourhood allows "
        "the same timing tertile and number of recurring values (2, 3, 4+), relaxed in that order and "
        "recorded as match_level (0 all, 1 values relaxed, 2 timing relaxed); a substitution value with "
        "no GTEx pair and no MPRA row anywhere in the unit; panel r-squared under 0.2 with the tested "
        "variant; within 250 kb for E2 and E3, anywhere on the chromosome for E1; the control borrows "
        "the unit's gene, tissue or cell and measured sign"
    ),
    "power": (
        "a difference of 0.10 in agreement needs about 495 pairs (80% power, one-sided alpha 0.01); with "
        "the 20 MPRA pairs chr21 and chr22 offer, only a difference of 0.5 or more is detectable, so a "
        "negative E1 here is a lack of power and not evidence of absence; E2's 135 pairs detect about "
        "0.19. Both primary endpoints together are 155 pairs, so the first look (150 pairs) falls near "
        "their end and the later looks matter only if E3 is run"
    ),
    "success": (
        "for an endpoint: agreement in units minus agreement in controls >= 0.10, one-sided Fisher exact "
        "p at or below the alpha of the look, and the same direction of difference on chr21 and chr22 "
        "separately; anything else is a negative"
    ),
    "secondary": [
        "absolute predicted effect larger in units than in controls (Mann-Whitney, reported only)",
        "rank correlation of predicted log2FC with GTEx slope among units (reported only)",
        "agreement against 0.5 without controls (reported only; the controls decide)",
        "a no-call grid of 0, 0.001 and 0.01 on |predicted log2FC|, added 2026-09-14 after a two-request "
        "smoke test showed single-variant gene effects of about 0.001, before any endpoint was scored; "
        "the primary keeps 0.001",
    ],
    "alpha_spending": {
        "look_1_pairs": 150,
        "alpha_1": 0.0005,
        "look_2_pairs": 400,
        "alpha_2": 0.002,
        "final_alpha": 0.01,
    },
    "stopping": (
        "requests go out in unit-control pairs in the order below; at each look an endpoint stops for "
        "success "
        "at its alpha, or for futility if the upper 95% bound of the agreement difference is below 0.05; "
        "the run stops when both endpoints have stopped, when the ordered list is exhausted, or when the "
        "quota holder asks for it back; a stopped run reports requests spent against units answered"
    ),
    "order": [
        "1: E1 MPRA units, strongest |log2FC| first",
        "2: E2 fine-mapped eQTL units, early-replicating with at most three values, smallest p first",
        "3: the other E2 units, smallest p first",
        "4: E3 linked eQTL units, only on request",
    ],
    "caveat": (
        "AlphaGenome learned from ENCODE RNA-seq, DNase and histone tracks, and GTEx expression is the same "
        "kind of measurement; agreement with an eQTL is partly the model reading back its inputs. MPRA "
        "reporter activity is not a training track in that way, so E1 is the cleaner test"
    ),
}

E2_REPLICATION: dict[str, Any] = {
    "written": (
        "2026-09-15, before any model request of the replication and before a pair of it was scored; the "
        "human panel had just reached 19 chromosomes and the shortlist 123,667 units with a measured value"
    ),
    "why": (
        "E2 is the strongest result this lane has (+0.268, units 38 of 56 against controls 30 of 73, p "
        "0.0021) and E3 showed its agreement follows causality rather than the neighbourhood. Both rest on "
        "161 pairs from chr21 and chr22, the two chromosomes that happened to be scored first, and chr22 "
        "carried the difference while chr21 did not reach significance alone. A claim that holds on those "
        "two and nowhere else has caught this project four times already, so the question is asked again "
        "where the answer can be different"
    ),
    "population": (
        "storage units on the catalogued chromosomes other than chr21 and chr22, which are the discovery "
        "set and are excluded from the replication claim entirely. E2R_eqtl_replication: the unit's "
        "recurring value is a variant DAP-G fine-maps (PIP at or above 0.5) for a gene with GTEx p at or "
        "below 1e-5 and one sign across tissues. E3R_eqtl_linked_replication: the other significant eQTL "
        "values, which are usually linked to a cause rather than causal"
    ),
    "unchanged": (
        "everything the first run fixed stays fixed: the read-out (the gene GTEx names, the track matching "
        "the most significant tissue else the mean over the gene's tracks), the measured sign (GTEx slope "
        "for the alternative allele), the no-call at 0.001 with the 0, 0.001 and 0.01 grid, two matched "
        "controls per unit from the same chromosome within 250 kb with the same GC stratum and, where the "
        "neighbourhood allows, the same timing tertile and value count, panel r-squared under 0.2, no "
        "measured value of their own, borrowing the unit's gene, tissue and measured sign, each unit "
        "counted once and each control once"
    ),
    "order": (
        "round robin over chromosomes, strongest evidence first inside each: the first requests answer "
        "whether the effect appears on chromosomes the original never saw, not which chromosome has the "
        "most units. A unit's two pairs stay together, so a stop leaves no half-answered unit. E2R is "
        "spent before E3R"
    ),
    "alpha_spending": {"look_1_pairs": 150, "alpha_1": 0.0005, "look_2_pairs": 400, "alpha_2": 0.002},
    "final_alpha": 0.01,
    "stopping": (
        "unchanged: at 150 and 400 pairs an endpoint stops for success at its alpha or for futility if the "
        "upper 95% bound of the difference is below 0.05; otherwise the run ends at the budget and is read "
        "at 0.01. E3R never stops the run and is read against E2R by the hold-out's own rule "
        "(compare_endpoints)"
    ),
    "replicates": (
        "the pooled difference on the replication chromosomes is at or above 0.10 with one-sided Fisher p "
        "at or below the alpha of the look, AND no single chromosome's removal takes the pooled difference "
        "below 0.05 (leave_one_out). The second half is the point: it is what would have caught a claim "
        "carried by one chromosome, and E2 itself would have passed it only because chr21 and chr22 both "
        "leaned the same way"
    ),
    "fails_to_replicate": (
        "the pooled difference is at or below 0 with the budget spent, or the run stops for futility. Then "
        "E2's +0.268 belongs to chr21 and chr22 or to 161 pairs, the executor claim loses the result it "
        "rested on, and E3's dilution becomes a statement about two chromosomes"
    ),
    "resolves_nothing": (
        "a pooled difference between 0 and 0.10, or p above the alpha with the upper 95% bound above 0.05, "
        "or a pooled difference that only clears 0.10 while one chromosome carries it (leave-one-out below "
        "0.05). The last case is reported as 'not replicated as a genome-wide claim', never as a success"
    ),
    "dilution_again": (
        "E3R is the hold-out asked at power: if a unit executes its value, E3R's difference should again be "
        "clearly smaller than E2R's, with the gap's one-sided 95% lower bound above zero. If the two match "
        "this time, the reading is the neighbourhood and not causality, whatever the first run found"
    ),
    "what_a_success_would_still_not_mean": (
        "this endpoint is one model agreeing with fine-mapping. The widened E1, the endpoint whose outcome "
        "is an external measurement, came back at +0.014 with an upper bound of 0.073 and resolved nothing, "
        "and the benchmark lane has since shown that a model producing a direction is evidence of having "
        "been scored rather than of an element doing anything. The only part of this work that does not "
        "depend on AlphaGenome is the measured-against-measured reading: MPRA direction against GTEx slope, "
        "0.670 where DAP-G fine-maps the variant and 0.518 over all of them, with no model in between. A "
        "replication here strengthens the pattern; it does not move the evidence outside the model"
    ),
    "caveat": (
        "the same ENCODE caveat as E2: the model learned from RNA-seq tracks and GTEx expression is the "
        "same kind of measurement. Five chromosomes (chr1 to chr5) have no panel yet, so 'genome-wide' "
        "means 17 replication chromosomes, and chrY contributes almost nothing"
    ),
}

TWO_INSTRUMENTS: dict[str, Any] = {
    "written": (
        "2026-09-15, before the genome-wide table was built and before any of its numbers were seen; the "
        "first reading of it, on 179 fine-mapped variants in one cell type, is in the result it replaces"
    ),
    "question": (
        "do two instruments agree about what an allele does, and does their agreement follow causality? "
        "A massively parallel reporter assay measures the alternative allele's effect on a plasmid "
        "fragment's transcription; a GTEx eQTL measures the same allele against a gene's expression in "
        "hundreds of donors' own chromosomes. DAP-G says whether the allele is probably the cause of the "
        "eQTL or a passenger of one. No model is involved at any point, which is why this is worth doing"
    ),
    "population": (
        "every variant with a significant MPRAVarDB row (FDR at or below 0.05, or nominal p at or below "
        "0.01 where the study reports none; 3'UTR stability libraries excluded) and a significant GTEx v8 "
        "pair at the same position with the same reference and alternative base. Every cell line, "
        "including those the model has no track for: this comparison does not need one"
    ),
    "the_eqtl_used": (
        "the pair DAP-G fine-maps at that very variant (highest PIP) where there is one, else the most "
        "significant pair at the variant. One row per (variant, cell line, study): the MPRA direction is a "
        "property of the cell line and the study that measured it, so rows are not pooled across them "
        "before they are reported"
    ),
    "split": (
        "fine-mapped, meaning DAP-G gives the variant a posterior at or above 0.5 for the gene used, "
        "against the rest, which are usually passengers. Declared before the table was built"
    ),
    "control": (
        "the same comparison on variants whose MPRA row is not significant (FDR above 0.2, or p above 0.2 "
        "where no FDR): two instruments, one of which measured nothing, should agree at chance in both "
        "strata. If the fine-mapped stratum of this control also reads near 0.67, the split is an artefact "
        "of the data's construction and not a property of the alleles"
    ),
    "reported": (
        "per cell line and per study with the sample size beside every figure, and pooled only after "
        "those, because the studies differ in library, cell and depth and a pooled figure hides which of "
        "them carries it"
    ),
    "holds_at_scale": (
        "the fine-mapped stratum agrees at 0.60 or better, it beats the rest by at least 0.08 with a "
        "one-sided Fisher p at or below 0.01, the same direction appears in the two largest cell lines "
        "separately, and the non-significant control sits near 0.5 in both strata"
    ),
    "thins_out": (
        "the difference between the strata is 0.03 or less, or the fine-mapped stratum agrees at 0.55 or "
        "less. Then the one piece of evidence in this area that owes nothing to the model is no better "
        "than the model-dependent ones, and that is the more important result, to be said as plainly as "
        "the E2 replication was"
    ),
    "resolves_nothing": "anything between those, reported as such and not as support",
    "what_it_cannot_mean": [
        "MPRA libraries are built from eQTL and GWAS variants, so the overlap between the two "
        "instruments is not a random sample of the genome and no statement about the genome follows from "
        "it; it is a statement about alleles that were already suspected of doing something",
        "a plasmid fragment is not a gene in its chromosome: agreement says the allele does something in "
        "both assays and in the same direction, not that the eQTL's gene is the reporter element's target",
        "fine-mapping is itself a statistical instrument, and the split is between what DAP-G calls causal "
        "and what it does not, not between truth and falsehood",
    ],
}

E1_WIDE: dict[str, Any] = {
    "written": "2026-09-14, before any model request of the widened endpoint and before any pair was scored",
    "why": (
        "E1 was the endpoint that answers the standing objection to E2 and E3 — that the model's agreement "
        "with a GTEx eQTL is partly the model reading back training data of the same kind — and it had 12 "
        "pairs, because the shortlist was built from the human panel's storage units and the panel has run "
        "only on chr21 and chr22. The measurement is not so limited: MPRAVarDB is genome-wide. This "
        "endpoint drops the storage-unit framing and tests every significant allele pair the database "
        "holds, in the cell lines the model has RNA-seq tracks for"
    ),
    "claim_tested": (
        "swapping an allele that a reporter assay measured as changing activity moves the model's predicted "
        "read-out of the local gene in the measured direction more often than swapping an allele the same "
        "assay measured in the same cell and found not to change activity"
    ),
    "population": (
        "MPRAVarDB rows on the main assembly with an FDR at or below 0.05, or a nominal p at or below 0.01 "
        "where the study reports no FDR; 3'UTR stability libraries excluded; single-base ref and alt; the "
        "cell line one of GM12878, HepG2, K562, Jurkat or HeLa, the MPRA lines the model has RNA-seq tracks "
        "for. One test per variant: the significant row with the largest absolute log2FC"
    ),
    "strata": (
        "E1W_mpra_wide, the headline, is every such variant outside the saturation-mutagenesis study; "
        "E1S_mpra_saturation is that study alone, which mutates a handful of whole elements and so is a "
        "within-element test at five loci, reported apart and never the headline. Both are declared here, "
        "before either was assembled or scored"
    ),
    "readout": (
        "the model's RNA-seq gene tracks in the MPRA's own cell line, both alleles in one request over the "
        "1 Mb window; the gene is chosen without the model: the gene DAP-G names at this very variant "
        "(highest PIP) if there is one, else the nearest protein-coding TSS in GENCODE. The control is read "
        "on the same gene in the same cell, so the no-call threshold cannot filter the two arms differently"
    ),
    "readout_changed_from_e1": (
        "the original E1 chose the read-out gene as the one with the largest absolute predicted effect in "
        "the cell when no eQTL named a gene, which is the model choosing its own read-out and favours the "
        "unit arm at the no-call threshold. The widened endpoint chooses the gene from DAP-G or the "
        "annotation alone, and refuses the pair rather than falling back on a model-chosen gene"
    ),
    "measured_sign": (
        "sign of MPRAVarDB log2FC, which UCSC's track description defines as log2(alt RNA/DNA) - "
        "log2(ref RNA/DNA): alternative over reference. The sign convention is checked once against GTEx "
        "lymphocyte eQTLs before the run and reported with the result"
    ),
    "controls": (
        "one per unit, never reused: a variant the same study measured in the same cell line and found not "
        "significant (FDR above 0.2, or nominal p above 0.2 where no FDR is given), with no significant row "
        "in any study or cell anywhere, single-base ref and alt, between 200 bp and 250 kb from the tested "
        "variant so that it is outside the tested variant's own 200-bp window and inside the prediction "
        "window; the nearest such variant is taken. It borrows the unit's gene, cell and measured sign. "
        "This is a measured null from the same library rather than the unmeasured storage unit E2 and E3 "
        "used, because off chr21 and chr22 there is no human panel to draw one from and because the "
        "same-library null is the stronger match: same assay, same cell, same neighbourhood"
    ),
    "no_call": "a predicted |log2FC| below 0.001 is a no-call in either arm, with the 0, 0.001 and 0.01 grid",
    "order": "strongest measured |log2FC| first, which is the subset where the measurement is most certain",
    "budget": (
        "1,000 pairs and at most 2,000 requests for E1W_mpra_wide; if the quota holder grants it, a further "
        "300 pairs and 600 requests for E1S_mpra_saturation, reported apart"
    ),
    "alpha_spending": {"look_1_pairs": 150, "alpha_1": 0.0005, "look_2_pairs": 400, "alpha_2": 0.002},
    "final_alpha": 0.01,
    "stopping": (
        "pairs in the pre-registered order; at 150 and 400 pairs an endpoint stops for success at its alpha "
        "or for futility if the upper 95% bound of the difference is below 0.05; otherwise the run ends at "
        "the budget and is read at alpha 0.01"
    ),
    "success": (
        "agreement in units minus agreement in matched nulls at or above 0.10, one-sided Fisher exact p at "
        "or below the alpha of the look, and the same direction of difference in GM12878 and in Jurkat "
        "separately (the two cells that carry the headline stratum)"
    ),
    "confirms_e2_e3": (
        "a success here is the external answer to the ENCODE objection: the model tracks a measured "
        "direction that is not one of its training tracks, and E2's +0.268 with E3's dilution then reads as "
        "a property of the sequence and not of the model's memory of GTEx"
    ),
    "contradicts_e2_e3": (
        "a difference at or below 0 with the budget spent, or a futility stop, says the model does not "
        "track an external measurement of direction at all, and E2 and E3 are then best explained by the "
        "model reproducing expression data of the kind it was trained on. The executor claim would be "
        "confined to eQTL-shaped read-outs and would not be a statement about sequence"
    ),
    "fails_to_resolve": (
        "a difference between 0 and 0.10, or p above 0.01 with the upper 95% bound above 0.05, decides "
        "nothing: the reporter measures an episomal fragment and the model predicts a chromosomal gene, so "
        "a small positive difference is what a weak but real correspondence would also look like"
    ),
    "secondary_fine_mapped_split": (
        "declared 2026-09-14 before any request, after the convention check and from it: among the tested "
        "GM12878 variants that also have a significant eQTL in EBV-transformed lymphocytes, the reporter's "
        "direction and the eQTL's direction agree 120 of 179 times (0.670) where DAP-G fine-maps the "
        "variant and 946 of 1,826 (0.518) over all of them. Two independent measurements therefore show "
        "the same split E2 and E3 found in the model. E1W is read again inside those strata: if the model "
        "is tracking the measurement rather than its own training data, its agreement should be higher on "
        "the fine-mapped variants than on the rest, by about the same margin. This is a prediction, not a "
        "criterion: the endpoint's success or failure is decided by the difference against the matched "
        "nulls as stated above"
    ),
    "convention": (
        "the control borrows the unit's measured sign, so mirroring the log2FC convention would flip both "
        "arms and change the sign of the difference: this endpoint stands or falls with the convention. It "
        "is therefore checked before any request, genome-wide, against GTEx slopes in EBV-transformed "
        "lymphocytes for the GM12878 variants (wide_sign_check), and the check is committed with the plan. "
        "A strongly negative difference at the end would be read first as evidence about the convention and "
        "only then as evidence about the claim"
    ),
    "caveat": (
        "an MPRA measures a short fragment on a plasmid, not the gene in its chromosome; the assay's own "
        "reproducibility across studies is modest; and the read-out gene is the nearest TSS for most "
        "variants, which is wrong for a third of distal elements by this project's own deletion work. All "
        "three blunt the endpoint rather than bias it towards success"
    ),
}


# ------------------------------------------------------------------------------------------------
# GTEx slopes for every storage unit, through eqtl.py's reader
# ------------------------------------------------------------------------------------------------
def unit_id(chrom: str, start: int, end: int) -> str:
    return f"{chrom}:{start}-{end}"


def load_catalogue(chrom: str, cache: Path = PANEL_CACHE) -> list[dict[str, Any]]:
    with gzip.open(cache / chrom / "storage_catalogue.json.gz", "rt") as fh:
        return json.load(fh)


def distil_gtex(chroms: list[str], progress=None, directory: Path = GTEX_DIR) -> dict[str, Any]:
    """Every significant GTEx v8 pair whose variant lies in a storage unit, one file per tissue."""
    from genomeos.attribution.eqtl import Intervals, distil

    iv = Intervals()
    for chrom in chroms:
        for u in load_catalogue(chrom):
            iv.add(chrom, u["start"], u["end"], unit_id(chrom, u["start"], u["end"]))
    return distil(iv.freeze(), knowledge=directory, progress=progress)


def gtex_by_variant(directory: Path = GTEX_DIR) -> dict[tuple[str, int], list[dict[str, Any]]]:
    """(chrom, 1-based position) -> significant pairs, every tissue."""
    out: dict[tuple[str, int], list[dict[str, Any]]] = {}
    for p in sorted(directory.glob("hits_*.tsv")):
        with p.open() as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                out.setdefault((f[1], int(f[2])), []).append(
                    {
                        "tissue": f[0],
                        "ref": f[3],
                        "alt": f[4],
                        "gene_id": f[5],
                        "slope": float(f[6]),
                        "pval": float(f[7]),
                    }
                )
    return out


def gene_ids(chroms: list[str]) -> dict[str, str]:
    """Ensembl gene id (unversioned) -> symbol, from GENCODE."""
    from genomeos.genome.annotation import Annotation, default_gencode

    out: dict[str, str] = {}
    for chrom in chroms:
        gff = default_gencode({chrom})
        if gff is None:
            continue
        for gid, g in Annotation.from_gff3(gff, {chrom}).genes.items():
            out[gid.split(".")[0]] = g.symbol
    return out


# ------------------------------------------------------------------------------------------------
# The shortlist
# ------------------------------------------------------------------------------------------------
def _float(x: Any) -> float | None:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(v) else v


def value_bucket(recurring_values: int) -> str:
    return "2" if recurring_values <= 2 else "3" if recurring_values == 3 else "4+"


def build_shortlist(
    chroms: list[str],
    gtex: dict[tuple[str, int], list[dict[str, Any]]],
    symbols: dict[str, str],
    mpra_rows: dict[str, list[dict]] | None = None,
    dapg_rows: dict[str, list[dict]] | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """(shortlist, pool): units with a measured value, and the storage units with none (the control pool)."""
    shortlist: list[dict[str, Any]] = []
    pool: list[dict[str, Any]] = []
    for chrom in chroms:
        panel = Panel(PANEL_CACHE / chrom)
        mp = Overlaps((mpra_rows or {}).get(chrom, []))
        dp = Overlaps((dapg_rows or {}).get(chrom, []))
        for u in load_catalogue(chrom):
            uid = unit_id(chrom, u["start"], u["end"])
            gc = gc_stratum(panel.gc_fraction([(u["start"], u["end"])]))
            base = {
                "unit": uid,
                "chrom": chrom,
                "start": u["start"],
                "end": u["end"],
                "tier": u["tier"],
                "block": u["block"],
                "gc_stratum": gc,
                "rt_stratum": u["rt_stratum"],
                "recurring_values": u["domain"]["recurring_values"],
                "value_bucket": value_bucket(u["domain"]["recurring_values"]),
                "kinds": u["kinds"],
                "domain": u["domain"],
                "values": u["values"][:4],
            }
            events_measured = []
            any_gtex = False
            for f in u["frequencies"]:
                if f["panel"]["minor"] < MIN_RECURRING or f["kind"] != "snv":
                    continue
                g = f.get("gnomad") or {}
                pairs = [
                    p
                    for p in gtex.get((chrom, f["pos"]), [])
                    if p["ref"] == g.get("ref") and p["alt"] == g.get("alt")
                ]
                any_gtex = any_gtex or bool(gtex.get((chrom, f["pos"])))
                mpra = [
                    {
                        "cell": r.get("cellLine"),
                        "log2fc": _float(r.get("log2FC")),
                        "fdr": _float(r.get("fdr")),
                        "pvalue": _float(r.get("pvalue")),
                        "study": (r.get("mpraStudy") or "")[:90],
                    }
                    for r in mp.over(f["pos"] - 1, f["pos"])
                    if r["chromStart"] == f["pos"] - 1 and r.get("alt") == g.get("alt")
                ]
                dapg = sorted(
                    {
                        (r.get("geneName"), r.get("tissue"), _float(r.get("pip")))
                        for r in dp.over(f["pos"] - 1, f["pos"])
                    },
                    key=lambda x: -(x[2] or 0),
                )
                if not pairs and not mpra:
                    continue
                by_gene: dict[str, dict[str, Any]] = {}
                for p in pairs:
                    sym = symbols.get(p["gene_id"], p["gene_id"])
                    cur = by_gene.get(sym)
                    if cur is None or p["pval"] < cur["pval"]:
                        by_gene[sym] = {**p, "gene": sym, "tissues": 0}
                    by_gene[sym]["tissues"] += 1
                for sym, row in by_gene.items():
                    signs = {
                        math.copysign(1, p["slope"])
                        for p in pairs
                        if symbols.get(p["gene_id"], p["gene_id"]) == sym
                    }
                    row["consistent_sign_across_tissues"] = len(signs) == 1
                events_measured.append(
                    {
                        "pos": f["pos"],
                        "rsid": g.get("rsid"),
                        "ref": g.get("ref"),
                        "alt": g.get("alt"),
                        "panel_minor": f["panel"]["minor"],
                        "panel_haplotypes": f["panel"]["haplotypes"],
                        "gnomad_af": g.get("af"),
                        "gnomad_an": g.get("an"),
                        "eqtl": sorted(by_gene.values(), key=lambda r: r["pval"])[:5],
                        "dapg": [{"gene": a, "tissue": b, "pip": c} for a, b, c in dapg[:5]],
                        "mpra": mpra,
                    }
                )
            if events_measured:
                shortlist.append({**base, "measured": events_measured})
            elif not any_gtex and not u["eqtl"] and not u["mpra"] and "snv" in u["kinds"]:
                pool.append(
                    {
                        **base,
                        "snvs": [
                            f
                            for f in u["frequencies"]
                            if f["kind"] == "snv" and f["panel"]["minor"] >= MIN_RECURRING and f.get("gnomad")
                        ],
                    }
                )
    return shortlist, pool


def mpra_significant(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """MPRA rows measured at FDR 0.05, or at nominal p 0.01 where the database gives no FDR."""
    out = []
    for m in rows:
        if not m["log2fc"] or MPRA_EXCLUDE in (m["study"] or ""):
            continue
        if (m["fdr"] is not None and m["fdr"] <= MPRA_FDR) or (
            m["fdr"] is None and (m["pvalue"] or 1) <= MPRA_P
        ):
            out.append(m)
    return out


def tested_variant(unit: dict[str, Any]) -> dict[str, Any] | None:
    """The unit's test: its best-measured value event, the endpoint it serves, gene, tissue or cell, sign.

    E1 takes a significant MPRA row. E2 takes an eQTL whose gene DAP-G fine-maps to this variant
    (PIP >= 0.5) with one sign across GTEx tissues. E3 takes any other eQTL at p <= 1e-5, where the
    variant may only be linked to the cause; it is exploratory and runs last.
    """
    rank = {"E1_mpra": 0, "E2_eqtl": 1, "E3_eqtl_linked": 2}
    best = None
    for ev in unit["measured"]:
        base = {"pos": ev["pos"], "ref": ev["ref"], "alt": ev["alt"], "rsid": ev["rsid"]}
        eq = [e for e in ev["eqtl"] if e["pval"] <= EQTL_P]
        fine = {g["gene"]: g for g in ev.get("dapg", []) if (g["pip"] or 0) >= DAPG_PIP}
        sig = mpra_significant(ev["mpra"])
        if sig:
            m = max(sig, key=lambda m: abs(m["log2fc"]))
            cand = {
                **base,
                "endpoint": "E1_mpra",
                "cell": m["cell"],
                "gene": eq[0]["gene"] if eq else None,
                "measured_sign": 1 if m["log2fc"] > 0 else -1,
                "measured_effect": m["log2fc"],
                "measured_fdr": m["fdr"],
                "study": m["study"],
                "rank_key": -abs(m["log2fc"]),
            }
        elif eq:
            primary = [e for e in eq if e["gene"] in fine and e.get("consistent_sign_across_tissues", True)]
            e = primary[0] if primary else eq[0]
            cand = {
                **base,
                "endpoint": "E2_eqtl" if primary else "E3_eqtl_linked",
                "tissue": fine[e["gene"]]["tissue"] if primary else e["tissue"],
                "gene": e["gene"],
                "measured_sign": 1 if e["slope"] > 0 else -1,
                "measured_effect": e["slope"],
                "pval": e["pval"],
                "pip": fine[e["gene"]]["pip"] if primary else None,
                "rank_key": e["pval"],
            }
        else:
            continue
        if best is None or (rank[cand["endpoint"]], cand["rank_key"]) < (
            rank[best["endpoint"]],
            best["rank_key"],
        ):
            best = cand
    return best


# ------------------------------------------------------------------------------------------------
# The replication: E2 and E3 asked again on every catalogued chromosome but the discovery pair
# ------------------------------------------------------------------------------------------------
DISCOVERY_CHROMS = ("chr21", "chr22")  # where E2 and E3 were found: excluded from the replication
GTEX_WIDE = LANE / "gtex_wide"  # the genome-wide distillation, kept apart from the first run's
REPLICATION_ENDPOINTS = ("E2R_eqtl_replication", "E3R_eqtl_linked_replication")
REPLICATION_OF = {"E2_eqtl": REPLICATION_ENDPOINTS[0], "E3_eqtl_linked": REPLICATION_ENDPOINTS[1]}


def catalogued_chroms(cache: Path = PANEL_CACHE) -> list[str]:
    """The chromosomes the human panel has catalogued, in genome order."""
    have = {p.parent.name for p in cache.glob("chr*/storage_catalogue.json.gz")}
    return [c for c in MAIN_CHROMS if c in have]


def shard_gtex(directory: Path = GTEX_WIDE) -> Path:
    """The distilled hits split into one file per chromosome, written once.

    Without this every chromosome's assembly would scan all 49 tissue files again.
    """
    out = directory / "by_chrom"
    if out.exists() and any(out.glob("*.tsv")):
        return out
    out.mkdir(parents=True, exist_ok=True)
    handles: dict[str, Any] = {}
    for p in sorted(directory.glob("hits_*.tsv")):
        with p.open() as fh:
            next(fh, None)
            for line in fh:
                chrom = line.split("\t", 3)[1]
                h = handles.get(chrom)
                if h is None:
                    h = handles[chrom] = (out / f"{chrom}.tsv").open("w")
                h.write(line)
    for h in handles.values():
        h.close()
    return out


def gtex_of_chrom(chrom: str, shards: Path) -> dict[tuple[str, int], list[dict[str, Any]]]:
    """One chromosome's significant pairs, keyed the way gtex_by_variant keys them."""
    out: dict[tuple[str, int], list[dict[str, Any]]] = {}
    p = shards / f"{chrom}.tsv"
    if not p.exists():
        return out
    with p.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out.setdefault((f[1], int(f[2])), []).append(
                {
                    "tissue": f[0],
                    "ref": f[3],
                    "alt": f[4],
                    "gene_id": f[5],
                    "slope": float(f[6]),
                    "pval": float(f[7]),
                }
            )
    return out


def _units_in_order(pairs: list[dict[str, Any]]) -> list[list[dict[str, Any]]]:
    """A chromosome's pairs grouped by unit, each group kept together and in its own order."""
    groups: dict[str, list[dict[str, Any]]] = {}
    for p in pairs:
        groups.setdefault(p["unit"], []).append(p)
    return list(groups.values())


def round_robin(pairs: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """One unit from each chromosome in turn, endpoint by endpoint: breadth before depth.

    The first requests then answer whether the effect appears away from the discovery chromosomes,
    which is the question, instead of filling up on whichever chromosome has the most units.
    """
    out: list[dict[str, Any]] = []
    for endpoint in REPLICATION_ENDPOINTS:
        queues: dict[str, list[list[dict[str, Any]]]] = {}
        for p in pairs:
            if p["endpoint"] == endpoint:
                queues.setdefault(p["unit"].split(":")[0], []).append(p)
        ordered = {c: _units_in_order(v) for c, v in queues.items()}
        while any(ordered.values()):
            for c in sorted(ordered, key=lambda c: MAIN_CHROMS.index(c)):
                if ordered[c]:
                    out.extend(ordered[c].pop(0))
    for i, p in enumerate(out, 1):
        p["order"] = i
    return out


def replication_pairs(
    chroms: list[str] | None = None, shards: Path | None = None, progress=None
) -> dict[str, Any]:
    """Assemble E2R and E3R chromosome by chromosome, keeping only the pairs. No model request."""
    chroms = chroms or [c for c in catalogued_chroms() if c not in DISCOVERY_CHROMS]
    shards = shards or shard_gtex()
    say = progress or (lambda *a: None)
    pairs: list[dict[str, Any]] = []
    counts: dict[str, Any] = {}
    unmatched: Counter = Counter()
    for chrom in chroms:
        gtex = gtex_of_chrom(chrom, shards)
        if not gtex:
            say(f"{chrom}: no distilled GTEx pairs, skipped")
            continue
        symbols = gene_ids([chrom])
        mpra = {chrom: track_over_units("mpravardb", chrom)}
        dapg = {chrom: track_over_units("gtex_dapg", chrom)}
        shortlist, pool = build_shortlist([chrom], gtex, symbols, mpra, dapg)
        tests = []
        for u in shortlist:
            t = tested_variant(u)
            if t and t["endpoint"] in REPLICATION_OF:
                tests.append({"unit": u, "test": t, "order_key": order_key(u, t)})
        matched = match_controls(tests, pool, random.Random(SEED))
        for p in matched["pairs"]:
            p["endpoint"] = REPLICATION_OF[p["endpoint"]]
            p.pop("order", None)
        pairs.extend(matched["pairs"])
        for k, v in matched["unmatched"].items():
            unmatched[REPLICATION_OF.get(k, k)] += v
        counts[chrom] = {
            "units_with_a_measured_value": len(shortlist),
            "tests": dict(Counter(REPLICATION_OF[t["test"]["endpoint"]] for t in tests)),
            "pairs": dict(Counter(p["endpoint"] for p in matched["pairs"])),
            "control_pool": len(pool),
        }
        say(f"{chrom}: {counts[chrom]['tests']} tests, {counts[chrom]['pairs']} pairs")
        del gtex, shortlist, pool, tests, matched
    return {"pairs": round_robin(pairs), "by_chromosome": counts, "unmatched": dict(unmatched)}


def leave_one_out(rows: list[dict[str, Any]], endpoint: str) -> dict[str, Any]:
    """The pooled difference with each chromosome removed in turn: does one chromosome carry it?"""
    chroms = sorted({d["chrom"] for d in rows if d["endpoint"] == endpoint})
    out = {c: tally([d for d in rows if d["chrom"] != c], endpoint) for c in chroms}
    diffs = [(c, t["difference"]) for c, t in out.items() if t["difference"] is not None]
    worst = min(diffs, key=lambda x: x[1]) if diffs else None
    return {
        "by_chromosome_removed": out,
        "lowest": {"chromosome_removed": worst[0], "difference": worst[1]} if worst else None,
        "holds_above_0.05": bool(worst and worst[1] >= 0.05),
    }


# ------------------------------------------------------------------------------------------------
# The widened E1: MPRAVarDB read genome-wide, its nulls as controls
# ------------------------------------------------------------------------------------------------
def mpra_genome_wide(cache: Path | None = None, progress=None) -> list[dict[str, Any]]:
    """Every MPRAVarDB row on the main assembly, through the panel's own bigBed reader, cached once."""
    from genomeos.attribution.human_panel import open_bigbed

    path = cache or LANE / "mpravardb_rows.json.gz"
    if path.exists():
        with gzip.open(path, "rt") as fh:
            return json.load(fh)
    bb = open_bigbed(MPRA_URL_KEY)
    rows: list[dict[str, Any]] = []
    for chrom in MAIN_CHROMS:
        if chrom not in bb.chroms:
            continue
        got = bb.query(chrom, [(0, bb.chroms[chrom][1])])
        rows.extend(got)
        if progress:
            progress(f"{chrom}: {len(got)} MPRA rows, {len(rows)} in all")
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        json.dump(rows, fh)
    return rows


def mpra_call(row: dict[str, Any]) -> dict[str, Any]:
    """A raw MPRAVarDB row as the fields this test reads."""
    return {
        "chrom": row["chrom"],
        "pos": int(row["chromStart"]) + 1,
        "ref": (row.get("ref") or "").strip(),
        "alt": (row.get("alt") or "").strip(),
        "rsid": row.get("rsid") or None,
        "cell": row.get("cellLine"),
        "log2fc": _float(row.get("log2FC")),
        "fdr": _float(row.get("fdr")),
        "pvalue": _float(row.get("pvalue")),
        "study": (row.get("mpraStudy") or "")[:90],
    }


def wide_usable(call: dict[str, Any]) -> bool:
    """A row this endpoint can use at all: a single-base swap, in a cell the model has tracks for."""
    return (
        call["cell"] in MPRA_CELLS
        and len(call["ref"]) == 1
        and len(call["alt"]) == 1
        and call["ref"] in "ACGT"
        and call["alt"] in "ACGT"
        and MPRA_EXCLUDE not in (call["study"] or "")
    )


def wide_split(rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """(tests, nulls): one significant row per variant, and the measured-null rows, by the population rule."""
    calls = [c for c in (mpra_call(r) for r in rows) if wide_usable(c)]
    significant: dict[tuple[str, int, str], dict[str, Any]] = {}
    ever_significant: set[tuple[str, int, str]] = set()
    for c in calls:
        key = (c["chrom"], c["pos"], c["alt"])
        if not mpra_significant([c]):
            continue
        ever_significant.add(key)
        cur = significant.get(key)
        if cur is None or abs(c["log2fc"]) > abs(cur["log2fc"]):
            significant[key] = c
    nulls = [
        c
        for c in calls
        if (c["chrom"], c["pos"], c["alt"]) not in ever_significant
        and (c["fdr"] if c["fdr"] is not None else c["pvalue"] or 1) > MPRA_NULL_P
    ]
    tests = []
    for c in significant.values():
        tests.append(
            {
                **c,
                "endpoint": "E1S_mpra_saturation" if c["study"].startswith(SATURATION) else "E1W_mpra_wide",
                "measured_sign": 1 if c["log2fc"] > 0 else -1,
                "rank_key": -abs(c["log2fc"]),
            }
        )
    return tests, nulls


def protein_coding_tss(chrom: str) -> tuple[list[int], list[str]]:
    """The chromosome's protein-coding transcription starts, sorted, with their symbols."""
    from genomeos.genome.annotation import Annotation, Strand, default_gencode

    gff = default_gencode({chrom})
    if gff is None:
        return [], []
    out = []
    for g in Annotation.from_gff3(gff, {chrom}).genes.values():
        if g.type != "protein_coding":
            continue
        loc = g.locus
        out.append(((loc.start + 1) if loc.strand is Strand.PLUS else loc.end, g.symbol))
    out.sort()
    return [p for p, _ in out], [s for _, s in out]


def nearest_gene(tss: tuple[list[int], list[str]], pos: int) -> tuple[str | None, int | None]:
    """The nearest protein-coding TSS and its distance; the model has no say in this choice."""
    starts, symbols = tss
    if not starts:
        return None, None
    i = bisect_left(starts, pos)
    best = min(
        (j for j in (i - 1, i) if 0 <= j < len(starts)), key=lambda j: abs(starts[j] - pos), default=None
    )
    return (symbols[best], abs(starts[best] - pos)) if best is not None else (None, None)


def dapg_at(chrom: str, positions: list[int]) -> dict[int, dict[str, Any]]:
    """The DAP-G gene with the highest PIP at each 1-based position, from the genome-wide track."""
    from genomeos.attribution.human_panel import _eqtl_pos, track_rows

    if not positions:
        return {}
    ivs = [(p - 1, p) for p in sorted(set(positions))]
    rows, _ = track_rows("gtex_dapg", chrom, ivs)
    best: dict[int, dict[str, Any]] = {}
    for r in rows:
        p = _eqtl_pos(r)
        if p is None:
            continue
        pip = _float(r.get("pip")) or 0.0
        pos1 = p + 1
        cur = best.get(pos1)
        if cur is None or pip > cur["pip"]:
            best[pos1] = {"gene": r.get("geneName"), "tissue": r.get("tissue"), "pip": pip}
    return best


def wide_pairs(rows: list[dict[str, Any]], progress=None) -> dict[str, Any]:
    """Tests, their read-out genes and their matched measured-null controls, chromosome by chromosome."""
    tests, nulls = wide_split(rows)
    by_chrom_null: dict[str, list[dict[str, Any]]] = {}
    for n in nulls:
        by_chrom_null.setdefault(n["chrom"], []).append(n)
    for v in by_chrom_null.values():
        v.sort(key=lambda c: c["pos"])
    pairs: list[dict[str, Any]] = []
    unmatched: Counter = Counter()
    gene_source: Counter = Counter()
    for chrom in sorted({t["chrom"] for t in tests}, key=lambda c: MAIN_CHROMS.index(c)):
        mine = [t for t in tests if t["chrom"] == chrom]
        tss = protein_coding_tss(chrom)
        dapg = dapg_at(chrom, [t["pos"] for t in mine])
        pool = by_chrom_null.get(chrom, [])
        positions = [c["pos"] for c in pool]
        used: set[tuple[int, str]] = set()
        for t in sorted(mine, key=lambda t: t["rank_key"]):
            d = dapg.get(t["pos"])
            if d and d["gene"]:
                gene, how = d["gene"], f"DAP-G at the variant (PIP {d['pip']:.3f})"
            else:
                gene, dist = nearest_gene(tss, t["pos"])
                how = f"nearest protein-coding TSS ({dist} bp)" if gene else "no gene"
            if not gene:
                unmatched[t["endpoint"] + ":no_gene"] += 1
                continue
            gene_source["dapg" if d and d["gene"] else "nearest_tss"] += 1
            i = bisect_left(positions, t["pos"])
            best = None
            for j in range(max(0, i - 400), min(len(pool), i + 400)):
                c = pool[j]
                if (c["pos"], c["alt"]) in used or c["study"] != t["study"] or c["cell"] != t["cell"]:
                    continue
                dist = abs(c["pos"] - t["pos"])
                if dist < WIDE_MIN_DISTANCE or dist > WIDE_MAX_DISTANCE:
                    continue
                if best is None or dist < best[1]:
                    best = (c, dist)
            if best is None:
                unmatched[t["endpoint"] + ":no_null"] += 1
                continue
            c, dist = best
            used.add((c["pos"], c["alt"]))
            pairs.append(
                {
                    "endpoint": t["endpoint"],
                    "unit": f"{chrom}:{t['pos']}",
                    "test": {
                        k: t[k] for k in ("pos", "ref", "alt", "rsid", "log2fc", "fdr", "pvalue", "study")
                    },
                    "control": {
                        "unit": f"{chrom}:{c['pos']}",
                        "pos": c["pos"],
                        "ref": c["ref"],
                        "alt": c["alt"],
                        "rsid": c["rsid"],
                        "log2fc": c["log2fc"],
                        "fdr": c["fdr"],
                        "pvalue": c["pvalue"],
                        "distance": c["pos"] - t["pos"],
                    },
                    "borrowed": {
                        "gene": gene,
                        "tissue": None,
                        "cell": t["cell"],
                        "measured_sign": t["measured_sign"],
                        "strict": True,
                    },
                    "gene_how": how,
                    "rank_key": t["rank_key"],
                }
            )
        if progress:
            progress(f"{chrom}: {len(mine)} tests, {len(pairs)} pairs so far")
    pairs.sort(key=lambda p: (WIDE_ENDPOINTS.index(p["endpoint"]), p["rank_key"]))
    for i, p in enumerate(pairs, 1):
        p["order"] = i
    return {
        "pairs": pairs,
        "unmatched": dict(unmatched),
        "gene_source": dict(gene_source),
        "tests": len(tests),
    }


# ------------------------------------------------------------------------------------------------
# Matched controls
# ------------------------------------------------------------------------------------------------
def r2(a: list[Any], b: list[Any]) -> float | None:
    """Panel r-squared between two biallelic columns of tokens (hg38 state against the rest)."""
    pairs = [
        (x == "=", y == "=")
        for x, y in zip(a, b, strict=True)
        if x is not None and y is not None and x != "*" and y != "*"
    ]
    n = len(pairs)
    if n < 10:
        return None
    pa = sum(x for x, _ in pairs) / n
    pb = sum(y for _, y in pairs) / n
    if pa in (0, 1) or pb in (0, 1):
        return None
    pab = sum(x and y for x, y in pairs) / n
    d = pab - pa * pb
    return d * d / (pa * (1 - pa) * pb * (1 - pb))


def snv_tokens(panel: Panel, pos1: int) -> list[Any] | None:
    for ev in events(panel, [(pos1 - 1, pos1)]):
        if ev.kind == "snv" and ev.start == pos1 - 1:
            return ev.tokens
    return None


def match_controls(
    tests: list[dict[str, Any]], pool: list[dict[str, Any]], rng: random.Random, k: int = CONTROLS_PER_UNIT
) -> dict[str, Any]:
    """One matched control per tested unit, never reused; unmatched units are counted."""
    panels: dict[str, Panel] = {}
    by_chrom: dict[str, list[dict[str, Any]]] = {}
    for c in pool:
        by_chrom.setdefault(c["chrom"], []).append(c)
    used: set[str] = set()
    pairs = []
    unmatched = Counter()
    for t in tests:
        u, v = t["unit"], t["test"]
        window = CONTROL_WINDOW  # E1 too, from 2026-09-14: see AMENDMENTS
        if u["chrom"] not in panels:
            panels[u["chrom"]] = Panel(PANEL_CACHE / u["chrom"])
        panel = panels[u["chrom"]]
        tv = snv_tokens(panel, v["pos"])
        near = [
            c
            for c in by_chrom.get(u["chrom"], [])
            if c["unit"] not in used
            and c["snvs"]
            and (window is None or abs(c["start"] - v["pos"]) <= window)
        ]
        rng.shuffle(near)
        chosen = []
        level = None
        for lvl, keys in enumerate(MATCH_LEVELS):
            for c in near:
                if any(c[key] != u[key] for key in keys) or c["unit"] in {x[0]["unit"] for x in chosen}:
                    continue
                snv = max(c["snvs"], key=lambda f: f["panel"]["minor"])
                ct = snv_tokens(panel, snv["pos"])
                if tv is None or ct is None:
                    continue
                rr = r2(tv, ct)
                if rr is not None and rr >= LD_MAX_R2:
                    continue
                chosen.append((c, snv, rr))
                level = lvl
                if len(chosen) == k:
                    break
            if len(chosen) == k:
                break
        if not chosen:
            unmatched[v["endpoint"]] += 1
            continue
        for c, snv, rr in chosen:
            used.add(c["unit"])
            g = snv["gnomad"]
            pairs.append(
                {
                    "endpoint": v["endpoint"],
                    "unit": u["unit"],
                    "test": {k2: v[k2] for k2 in v if k2 != "rank_key"},
                    "control": {
                        "unit": c["unit"],
                        "pos": snv["pos"],
                        "ref": g["ref"],
                        "alt": g["alt"],
                        "rsid": g.get("rsid"),
                        "panel_minor": snv["panel"]["minor"],
                        "gnomad_af": g.get("af"),
                        "distance": snv["pos"] - v["pos"],
                        "r2_with_test": round(rr, 3) if rr is not None else None,
                    },
                    "borrowed": {
                        "gene": v.get("gene"),
                        "tissue": v.get("tissue"),
                        "cell": v.get("cell"),
                        "measured_sign": v["measured_sign"],
                    },
                    "strata": {
                        "gc": u["gc_stratum"],
                        "rt": u["rt_stratum"],
                        "values": u["value_bucket"],
                    },
                    "match_level": level,
                    "order_key": t["order_key"],
                }
            )
    pairs.sort(key=lambda p: p["order_key"])
    for i, p in enumerate(pairs, 1):
        p["order"] = i
    return {"pairs": pairs, "unmatched": dict(unmatched)}


def order_key(unit: dict[str, Any], test: dict[str, Any]) -> tuple:
    """The pre-registered order: MPRA, then fine-mapped eQTLs early and small, the other fine-mapped,
    then the linked eQTLs last."""
    if test["endpoint"] == "E1_mpra":
        return (1, test["rank_key"])
    if test["endpoint"] == "E2_eqtl":
        early_small = unit["rt_stratum"] == 2 and unit["recurring_values"] <= 3
        return (2 if early_small else 3, test["rank_key"])
    return (4, test["rank_key"])


# ------------------------------------------------------------------------------------------------
# Before any request: is MPRAVarDB's log2FC alternative over reference?
# ------------------------------------------------------------------------------------------------
LCL_TISSUE = "Cells_EBV-transformed_lymphocytes"


def sign_convention_check(shortlist: list[dict[str, Any]]) -> dict[str, Any]:
    """Among variants with a significant MPRA row in a lymphoblastoid line and a GTEx eQTL in EBV
    lymphocytes, how often the two signs agree. Well above half supports reading log2FC as
    alternative over reference; well below half says the opposite; near half settles nothing."""
    agree = n = 0
    for u in shortlist:
        for ev in u["measured"]:
            lcl = [e for e in ev["eqtl"] if e["tissue"] == LCL_TISSUE]
            mp = [
                m
                for m in ev["mpra"]
                if (m["cell"] or "").upper().startswith("GM")
                and m["log2fc"]
                and m["fdr"] is not None
                and m["fdr"] <= MPRA_FDR
            ]
            if lcl and mp:
                n += 1
                agree += (lcl[0]["slope"] > 0) == (mp[0]["log2fc"] > 0)
    p = binom_two_sided(agree, n) if n else None
    return {
        "variants": n,
        "signs_agree": agree,
        "share": round(agree / n, 3) if n else None,
        "binomial_p_two_sided": p,
        "reading": (
            "log2FC read as alternative over reference"
            if n and agree / n > 0.5 and p is not None and p <= 0.05
            else "convention reversed relative to GTEx"
            if n and agree / n < 0.5 and p is not None and p <= 0.05
            else "not settled by these variants; E1 is reported under the stated convention and its mirror"
        ),
    }


BINOM_EXACT_MAX = 1000  # above this the exact sum overflows a float; the normal approximation is used


def binom_two_sided(k: int, n: int, p: float = 0.5) -> float:
    if n > BINOM_EXACT_MAX:
        z = abs(k - n * p) / math.sqrt(n * p * (1 - p))
        return round(math.erfc(z / math.sqrt(2)), 8)
    probs = [math.comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(n + 1)]
    obs = probs[k]
    return round(min(1.0, sum(q for q in probs if q <= obs * (1 + 1e-9))), 6)


def fisher_greater(a: int, n1: int, c: int, n2: int) -> float:
    """One-sided Fisher exact p that the first proportion (a of n1) exceeds the second (c of n2)."""
    total = a + c
    n = n1 + n2

    def hyper(x: int) -> float:
        return math.comb(n1, x) * math.comb(n2, total - x) / math.comb(n, total)

    lo = max(0, total - n2)
    hi = min(n1, total)
    return round(min(1.0, sum(hyper(x) for x in range(max(a, lo), hi + 1))), 8)


def difference_upper(a: int, n1: int, c: int, n2: int, z: float = 1.645) -> float | None:
    """Upper one-sided 95% bound of p1 - p2 (normal approximation)."""
    if not n1 or not n2:
        return None
    p1, p2 = a / n1, c / n2
    se = math.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
    return p1 - p2 + z * se


# ------------------------------------------------------------------------------------------------
# Reading a prediction, scoring the pairs
# ------------------------------------------------------------------------------------------------
NO_CALL = 0.001


def _norm(name: str | None) -> str:
    return "".join(ch for ch in (name or "").lower() if ch.isalnum())


def _pick_gene(tracks: list[tuple[str, str, float]]) -> str | None:
    best: dict[str, float] = {}
    for g, _, v in tracks:
        best[g] = max(best.get(g, 0.0), abs(v))
    return max(best, key=best.get) if best else None


def readout(
    effects: list[list[Any]],
    gene: str | None,
    tissue: str | None,
    cell: str | None,
    fallback: bool = True,
) -> dict[str, Any]:
    """The predicted log2FC for the read-out the criterion names, and how it was chosen.

    With fallback off (the widened E1) the named gene is the only read-out: the model never chooses
    its own, so the no-call threshold cannot filter the two arms differently.
    """
    rows = [(g, t, float(v)) for g, t, v in effects]
    if cell:
        rows = [r for r in rows if _norm(cell) in _norm(r[1])]
        if not rows:
            return {"value": None, "gene": gene, "how": "no track for the cell"}
        gene = gene or (_pick_gene(rows) if fallback else None)
        if gene is None:
            return {"value": None, "gene": None, "how": "no gene named for a strict read-out"}
    mine = [r for r in rows if r[0] == gene]
    if not mine and cell and fallback:
        gene = _pick_gene(rows)
        mine = [r for r in rows if r[0] == gene]
        if mine:
            return {
                "value": sum(r[2] for r in mine) / len(mine),
                "gene": gene,
                "how": f"the named gene is outside this window; strongest gene in {cell}",
            }
    if not mine:
        return {"value": None, "gene": gene, "how": "gene not in the window's output"}
    if tissue:
        match = [r for r in mine if _norm(r[1]) == _norm(tissue)]
        if match:
            return {"value": match[0][2], "gene": gene, "how": f"track {match[0][1]}"}
    if cell:
        v = sum(r[2] for r in mine) / len(mine)
        return {"value": v, "gene": gene, "how": f"mean of {len(mine)} {cell} tracks"}
    v = sum(r[2] for r in mine) / len(mine)
    return {"value": v, "gene": gene, "how": f"mean of the gene's {len(mine)} tracks (no tissue match)"}


Scorer = Callable[[str, int, str, str], list[tuple[str, str, float]]]


def live_scorer() -> Scorer:
    """The existing variant path, unchanged: the adapter's RNA-seq gene scorer with no threshold.

    Building the scorer does not call the model; each call of the scorer is one request.
    """
    from genomeos.predict.alphagenome_adapter import AlphaGenomeAdapter

    return AlphaGenomeAdapter()._live_scorer(threshold=0.0)


def kept(scorer: Scorer, gene: str | None, cell: str | None = None, keep_min: float = 0.05) -> Scorer:
    """Every track of the gene under test, every track of the cell under test, and anything that moves.

    The filter only shrinks the cache; what it keeps must cover the read-out of both arms, which is why
    a cell read-out keeps that cell's tracks for every gene in the window.
    """

    def score(chrom: str, pos: int, ref: str, alt: str) -> list[tuple[str, str, float]]:
        return [
            (g, t, v)
            for g, t, v in scorer(chrom, pos, ref, alt)
            if g == gene or abs(v) >= keep_min or (cell and _norm(cell) in _norm(t))
        ]

    return score


QUOTA_WAIT = 30.0  # seconds to wait when the service reports its per-minute token quota exhausted
QUOTA_TRIES = 8


def with_quota_waits(fn, tries: int = QUOTA_TRIES, wait: float = QUOTA_WAIT, progress=None):
    """Call fn, waiting out the service's per-minute quota; anything else is raised at once."""
    for attempt in range(tries):
        try:
            return fn()
        except Exception as e:  # noqa: BLE001 - the client raises its own error types
            if "RESOURCE_EXHAUSTED" not in str(e) or attempt == tries - 1:
                raise
            if progress:
                progress(f"quota exhausted, waiting {wait:.0f} s ({attempt + 1}/{tries})")
            time.sleep(wait)
    raise RuntimeError("unreachable")  # pragma: no cover


def predict_side(
    scorer: Scorer, chrom: str, side: dict[str, Any], borrowed: dict[str, Any], cache: Path
) -> dict[str, Any]:
    """One request (or a cached answer) for one variant, read the way the criterion says."""
    from genomeos.predict.individual_effects import cache_path, score_variant

    cell = borrowed.get("cell")
    folder = cache / (f"cell_{cell}" if cell else borrowed["gene"])
    hit = cache_path(chrom, side["pos"], side["ref"], side["alt"], folder).exists()
    effects = score_variant(
        kept(scorer, borrowed.get("gene"), cell), chrom, side["pos"], side["ref"], side["alt"], folder
    )
    r = readout(
        effects,
        borrowed.get("gene"),
        borrowed.get("tissue"),
        borrowed.get("cell"),
        fallback=not borrowed.get("strict"),
    )
    v = r["value"]
    call = None if v is None or abs(v) < NO_CALL else (1 if v > 0 else -1)
    return {
        **r,
        "request": not hit,
        "call": call,
        "agrees": None if call is None else call == borrowed["measured_sign"],
    }


def recall(done: list[dict[str, Any]], endpoint: str, no_call: float) -> dict[str, Any]:
    """The endpoint scored again at another no-call threshold, from the values already paid for."""
    rows = []
    for d in done:
        if d["endpoint"] != endpoint:
            continue
        out = dict(d)
        for side in ("test_result", "control_result"):
            r = dict(d[side])
            v = r.get("value")
            call = None if v is None or abs(v) < no_call else (1 if v > 0 else -1)
            r["call"] = call
            r["agrees"] = None if call is None else call == d["borrowed"]["measured_sign"]
            out[side] = r
        rows.append(out)
    return {"no_call": no_call, **tally(rows, endpoint)}


def tally(done: list[dict[str, Any]], endpoint: str, chrom: str | None = None) -> dict[str, Any]:
    """Sign agreement in units and in their matched controls, with the one-sided test."""
    rows = [d for d in done if d["endpoint"] == endpoint and (chrom is None or d["chrom"] == chrom)]
    # a unit with two controls appears in two pairs: it is counted once, each control once
    units = {d["unit"]: d["test_result"]["agrees"] for d in rows}
    controls = {(d["unit"], d["control"]["unit"]): d["control_result"]["agrees"] for d in rows}
    u = [x for x in units.values() if x is not None]
    c = [x for x in controls.values() if x is not None]
    a, n1, k, n2 = sum(u), len(u), sum(c), len(c)
    d = (a / n1 - k / n2) if n1 and n2 else None
    return {
        "pairs": len(rows),
        "units_called": n1,
        "units_agree": a,
        "controls_called": n2,
        "controls_agree": k,
        "difference": round(d, 3) if d is not None else None,
        "p_one_sided": fisher_greater(a, n1, k, n2) if n1 and n2 else None,
        "difference_upper_95": round(difference_upper(a, n1, k, n2) or 0.0, 3) if n1 and n2 else None,
    }


def compare_endpoints(
    rows_a: list[dict[str, Any]],
    rows_b: list[dict[str, Any]],
    a: str = "E2_eqtl",
    b: str = "E3_eqtl_linked",
    z: float = 1.645,
) -> dict[str, Any]:
    """The hold-out read: one endpoint's agreement difference minus the other's, with a one-sided bound.

    The rule is HOLD_OUT['read'], fixed before any row of the second endpoint was read. A gap whose
    lower bound clears zero is the dilution the executor claim predicts; no gap, with the second
    endpoint meeting the success criterion on its own, says the difference belongs to the
    neighbourhood and not to the value being the cause.
    """
    ta, tb = tally(rows_a, a), tally(rows_b, b)
    if ta["difference"] is None or tb["difference"] is None:
        return {a: ta, b: tb, "gap": None, "reading": "one endpoint has no answered pair"}
    var = 0.0
    for t in (ta, tb):
        for agree, n in (
            (t["units_agree"], t["units_called"]),
            (t["controls_agree"], t["controls_called"]),
        ):
            p = agree / n
            var += p * (1 - p) / n
    gap = ta["difference"] - tb["difference"]
    lower = gap - z * math.sqrt(var)
    holds = tb["difference"] >= 0.10 and (tb["p_one_sided"] or 1) <= 0.01
    return {
        a: ta,
        b: tb,
        "gap": round(gap, 3),
        "gap_lower_95": round(lower, 3),
        "second_endpoint_meets_success": holds,
        "reading": (
            f"{b} is diluted against {a} as the executor claim predicts"
            if lower > 0
            else f"{b} matches {a}: the difference is a property of the neighbourhood, not of the value "
            "being the cause"
            if holds
            else "neither separated: the two endpoints are not distinguished and the second does not "
            "stand on its own"
        ),
    }


def group_of(row: dict[str, Any], key: str) -> Any:
    """The value a split names: a column of the pair, or one the control borrowed (the cell line)."""
    return row["borrowed"].get(key) if key in row.get("borrowed", {}) else row.get(key)


def by_group(done: list[dict[str, Any]], endpoint: str, key: str) -> dict[str, dict[str, Any]]:
    """The endpoint tallied inside each value of a split (chromosome, or the MPRA's cell line)."""
    rows = [d for d in done if d["endpoint"] == endpoint]
    out = {}
    for v in sorted({str(group_of(d, key)) for d in rows}):
        out[v] = tally([d for d in rows if str(group_of(d, key)) == v], endpoint)
    return out


def decide(
    done: list[dict[str, Any]], endpoint: str, alpha: float, final: bool = False, split: str = "chrom"
) -> str | None:
    """'success', 'futility', 'negative' (at the final look) or None (carry on), by the criterion."""
    t = tally(done, endpoint)
    if not t["units_called"] or not t["controls_called"]:
        return "negative" if final else None
    per_group = list(by_group(done, endpoint, split).values())
    same_way = all(
        (x["difference"] or 0) > 0 for x in per_group if x["units_called"] and x["controls_called"]
    )
    if t["difference"] >= 0.10 and t["p_one_sided"] <= alpha and same_way:
        return "success"
    if t["difference_upper_95"] < 0.05:
        return "futility"
    return "negative" if final else None


def run_pairs(
    pairs: list[dict[str, Any]],
    scorer: Scorer,
    cache: Path = PREDICTIONS,
    max_requests: int | None = None,
    progress=None,
    max_pairs: int | None = None,
    looks_at: tuple[str, ...] = ("E1_mpra", "E2_eqtl"),
    endpoints: tuple[str, ...] = ENDPOINTS,
    split: str = "chrom",
) -> dict[str, Any]:
    """Spend requests in the pre-registered order, look at the pre-registered points, stop by the rule.

    max_pairs stops at a pair count as well as at a request count: the budget is written in both, and a
    rerun over the cache answers the same pairs again for nothing.
    """
    looks = CRITERION["alpha_spending"]
    done: list[dict[str, Any]] = []
    requests = 0
    stopped: dict[str, str] = {}
    for p in pairs:
        if (looks_at and len(stopped) >= len(looks_at)) or (
            max_requests is not None and requests >= max_requests
        ):
            break
        if max_pairs is not None and len(done) >= max_pairs:
            break
        if p["endpoint"] in stopped:
            continue
        chrom = p["unit"].split(":")[0]
        t = with_quota_waits(
            functools.partial(predict_side, scorer, chrom, p["test"], p["borrowed"], cache), progress=progress
        )
        c = with_quota_waits(
            functools.partial(predict_side, scorer, chrom, p["control"], p["borrowed"], cache),
            progress=progress,
        )
        requests += int(t["request"]) + int(c["request"])
        done.append({**p, "chrom": chrom, "test_result": t, "control_result": c})
        n = len(done)
        for key, alpha in (("look_1_pairs", looks["alpha_1"]), ("look_2_pairs", looks["alpha_2"])):
            if n == looks[key]:
                for e in looks_at:  # an exploratory endpoint is not looked at and never stops the run
                    call = decide(done, e, alpha, split=split)
                    if call and e not in stopped:
                        stopped[e] = f"{call} at {n} pairs"
        if progress and n % 25 == 0:
            progress(f"{n} pairs, {requests} requests, {stopped or 'running'}")
    out = {"pairs_done": len(done), "requests_spent": requests, "stopped": stopped}
    for e in endpoints:
        out[e] = tally(done, e)
        out[e]["by_chromosome"] = {c: tally(done, e, c) for c in sorted({d["chrom"] for d in done})}
        if split != "chrom":
            out[e][f"by_{split}"] = by_group(done, e, split)
        out[e]["verdict"] = stopped.get(e) or (
            decide(done, e, looks["final_alpha"], final=True, split=split) or "negative"
        )
    out["rows"] = done
    return out


def wide_sign_check(pairs: list[dict[str, Any]], progress=None) -> dict[str, Any]:
    """The direction convention, checked genome-wide before any request and without one.

    MPRAVarDB's log2FC is read as alternative over reference. GTEx's slope is signed the same way,
    so among the tested GM12878 variants that also have a significant eQTL in EBV-transformed
    lymphocytes the two signs should agree more often than not. This settles nothing about the
    executor claim; it settles whether the endpoint's measured sign is the right way up.
    """
    from genomeos.attribution.eqtl import Intervals, distil, hits_path

    tested = [p for p in pairs if p["borrowed"]["cell"] == "GM12878"]
    directory = LANE / "gtex_lcl"
    if not hits_path(LCL_TISSUE, directory).exists():
        iv = Intervals()
        for p in tested:
            chrom, pos = p["unit"].split(":")[0], int(p["unit"].split(":")[1])
            iv.add(chrom, pos - 1, pos, p["unit"])
        distil(iv.freeze(), knowledge=directory, tissues=[LCL_TISSUE], progress=progress)
    slopes: dict[tuple[str, int], float] = {}
    with hits_path(LCL_TISSUE, directory).open() as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            slopes[(f[1], int(f[2]))] = float(f[6])
    agree = n = 0
    fine_agree = fine_n = 0
    for p in tested:
        chrom, pos = p["unit"].split(":")[0], int(p["unit"].split(":")[1])
        s = slopes.get((chrom, pos))
        if s is None:
            continue
        same = (s > 0) == (p["test"]["log2fc"] > 0)
        n += 1
        agree += same
        m = re.search(r"PIP ([0-9.]+)", p.get("gene_how") or "")
        if m and float(m.group(1)) >= DAPG_PIP:
            fine_n += 1
            fine_agree += same
    share = round(agree / n, 3) if n else None
    return {
        "tissue": LCL_TISSUE,
        "variants": n,
        "signs_agree": agree,
        "share": share,
        "binomial_p_two_sided": binom_two_sided(agree, n) if n else None,
        "fine_mapped_only": {
            "variants": fine_n,
            "signs_agree": fine_agree,
            "share": round(fine_agree / fine_n, 3) if fine_n else None,
            "binomial_p_two_sided": binom_two_sided(fine_agree, fine_n) if fine_n else None,
            "note": (
                "the variants DAP-G fine-maps at this very position: the population E2 tested, where "
                "linkage does not stand between the allele and the eQTL's direction"
            ),
        },
        "reading": (
            "log2FC read as alternative over reference, as the endpoint assumes"
            if share and share > 0.5
            else "the convention may be mirrored; read the endpoint's sign with that in mind"
            if share
            else "no overlap: the convention is not checked by these variants"
        ),
    }


def replication_summary(built: dict[str, Any], examples: int = 8) -> dict[str, Any]:
    """The committed plan of the replication: what is available, what the first requests buy, the budget."""
    pairs = built["pairs"]
    per = {e: [p for p in pairs if p["endpoint"] == e] for e in REPLICATION_ENDPOINTS}
    units = {e: len({p["unit"] for p in v}) for e, v in per.items()}
    first = {e: [p for p in per[e][:400]] for e in REPLICATION_ENDPOINTS}
    budget = {}
    # the ask, fixed here before the key arrives: E2R buys the replication, E3R the dilution contrast
    for e, want in (("E2R_eqtl_replication", 600), ("E3R_eqtl_linked_replication", 400)):
        u = min(want, units[e])
        budget[e] = {
            "units": u,
            "pairs": sum(1 for p in per[e][: 2 * u]),
            "requests": 3 * u,  # one for the unit, one for each of its two controls
        }
    return {
        "chromosomes": sorted({p["unit"].split(":")[0] for p in pairs}, key=MAIN_CHROMS.index),
        "discovery_chromosomes_excluded": list(DISCOVERY_CHROMS),
        "by_chromosome": built["by_chromosome"],
        "unmatched": built["unmatched"],
        "pairs": {e: len(v) for e, v in per.items()},
        "units": units,
        "chromosomes_in_the_first_400_pairs": {
            e: dict(Counter(p["unit"].split(":")[0] for p in v)) for e, v in first.items()
        },
        "match_levels": {e: dict(Counter(str(p.get("match_level")) for p in v)) for e, v in per.items()},
        "measured_positive_share": {
            e: _share([p["borrowed"]["measured_sign"] > 0 for p in v]) for e, v in per.items()
        },
        "control_distance_median_kb": {
            e: _median_abs([p["control"]["distance"] for p in v], 1000) for e, v in per.items()
        },
        "budget": budget,
        "pre_registration": E2_REPLICATION,
        "first_pairs": [
            {k: p[k] for k in ("order", "endpoint", "unit", "test", "control", "borrowed", "match_level")}
            for p in pairs[:examples]
        ],
        "model_requests_spent": 0,
    }


GTEX_AT_MPRA = LANE / "gtex_at_mpra"  # GTEx distilled over the MPRA-tested positions, not over units
MATCHED_TISSUE = {  # the GTEx tissue closest to each MPRA cell line, for the declared secondary
    "GM12878": "Cells_EBV-transformed_lymphocytes",
    "Jurkat": "Cells_EBV-transformed_lymphocytes",
    "HepG2": "Liver",
    "K562": "Whole_Blood",
    "HEL92.1.7": "Whole_Blood",
}


def distil_gtex_at(positions: list[tuple[str, int]], directory: Path = GTEX_AT_MPRA, progress=None):
    """Every significant GTEx pair at the given 1-based positions, one file per tissue. No model call."""
    from genomeos.attribution.eqtl import Intervals, distil

    iv = Intervals()
    for chrom, pos in positions:
        iv.add(chrom, pos - 1, pos, f"{chrom}:{pos}")
    return distil(iv.freeze(), knowledge=directory, progress=progress)


def gtex_at_positions(directory: Path = GTEX_AT_MPRA) -> dict[tuple[str, int], list[dict[str, Any]]]:
    """The distilled pairs, keyed by (chromosome, 1-based position)."""
    out: dict[tuple[str, int], list[dict[str, Any]]] = {}
    for p in sorted(directory.glob("hits_*.tsv")):
        with p.open() as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                out.setdefault((f[1], int(f[2])), []).append(
                    {
                        "tissue": f[0],
                        "ref": f[3],
                        "alt": f[4],
                        "gene_id": f[5],
                        "slope": float(f[6]),
                        "pval": float(f[7]),
                    }
                )
    return out


def two_instrument_rows(
    mpra: list[dict[str, Any]], gtex: dict[tuple[str, int], list[dict[str, Any]]], progress=None
) -> list[dict[str, Any]]:
    """One row per (variant, cell line, study) that both instruments measured, with the two directions.

    The eQTL used is the one DAP-G fine-maps at that very variant, else the most significant pair
    there. Rows carry the stratum, the cell, the study and both signs; nothing is pooled here.
    """
    calls = [mpra_call(r) for r in mpra]
    wanted: dict[str, list[int]] = {}
    for c in calls:
        if (c["chrom"], c["pos"]) in gtex:
            wanted.setdefault(c["chrom"], []).append(c["pos"])
    dapg: dict[tuple[str, int], dict[str, Any]] = {}
    for chrom, positions in sorted(wanted.items()):
        for pos, best in dapg_at(chrom, positions).items():
            dapg[(chrom, pos)] = best
        if progress:
            progress(f"{chrom}: {len(positions)} positions with both instruments")
    rows = []
    for c in calls:
        if c["log2fc"] is None or MPRA_EXCLUDE in (c["study"] or ""):
            continue
        pairs = [
            p for p in gtex.get((c["chrom"], c["pos"]), []) if p["ref"] == c["ref"] and p["alt"] == c["alt"]
        ]
        if not pairs:
            continue
        fine = dapg.get((c["chrom"], c["pos"]))
        used = None
        if fine and (fine["pip"] or 0) >= DAPG_PIP:
            named = [p for p in pairs if p["tissue"] == fine["tissue"]]
            used = min(named or pairs, key=lambda p: p["pval"])
        else:
            used = min(pairs, key=lambda p: p["pval"])
        significant = bool(mpra_significant([c]))
        null = (c["fdr"] if c["fdr"] is not None else (c["pvalue"] or 1)) > MPRA_NULL_P
        if not significant and not null:
            continue
        rows.append(
            {
                "variant": f"{c['chrom']}:{c['pos']}",
                "rsid": c["rsid"],
                "cell": c["cell"],
                "study": c["study"][:60],
                "mpra_log2fc": c["log2fc"],
                "mpra_significant": significant,
                "gtex_tissue": used["tissue"],
                "gtex_gene": used["gene_id"],
                "gtex_slope": used["slope"],
                "gtex_pval": used["pval"],
                "pip": (fine or {}).get("pip"),
                "fine_mapped": bool(fine and (fine["pip"] or 0) >= DAPG_PIP),
                "tissue_matches_cell": MATCHED_TISSUE.get(c["cell"]) == used["tissue"],
                "agree": (c["log2fc"] > 0) == (used["slope"] > 0),
            }
        )
    return rows


def _agreement(rows: list[dict[str, Any]]) -> dict[str, Any]:
    n = len(rows)
    k = sum(r["agree"] for r in rows)
    return {
        "rows": n,
        "agree": k,
        "share": round(k / n, 3) if n else None,
        "binomial_p_two_sided": binom_two_sided(k, n) if n else None,
    }


def two_instruments(rows: list[dict[str, Any]], significant: bool = True) -> dict[str, Any]:
    """The reading TWO_INSTRUMENTS asks for: split by fine-mapping, reported per cell and per study."""
    mine = [r for r in rows if r["mpra_significant"] is significant]
    fine = [r for r in mine if r["fine_mapped"]]
    rest = [r for r in mine if not r["fine_mapped"]]
    out: dict[str, Any] = {
        "fine_mapped": _agreement(fine),
        "the_rest": _agreement(rest),
        "matched_tissue_only": {
            "fine_mapped": _agreement([r for r in fine if r["tissue_matches_cell"]]),
            "the_rest": _agreement([r for r in rest if r["tissue_matches_cell"]]),
        },
    }
    a, b = out["fine_mapped"], out["the_rest"]
    if a["rows"] and b["rows"]:
        out["difference"] = round(a["share"] - b["share"], 3)
        out["p_one_sided"] = fisher_greater(a["agree"], a["rows"], b["agree"], b["rows"])
    for key, field in (("by_cell", "cell"), ("by_study", "study")):
        groups = sorted({r[field] for r in mine})
        out[key] = {
            g: {
                "fine_mapped": _agreement([r for r in fine if r[field] == g]),
                "the_rest": _agreement([r for r in rest if r[field] == g]),
            }
            for g in groups
        }
    return out


def fine_mapped_split(rows: list[dict[str, Any]], endpoint: str) -> dict[str, Any]:
    """E1_WIDE's declared secondary: the endpoint inside the fine-mapped variants and outside them.

    The two measurements agree 0.670 where DAP-G fine-maps the variant and 0.518 elsewhere; if the
    model tracks the measurement rather than its own training data, its agreement should split the
    same way. The read-out gene records the posterior it was chosen with, which is the split here.
    """

    def pip(row: dict[str, Any]) -> float:
        m = re.search(r"PIP ([0-9.]+)", row.get("gene_how") or "")
        return float(m.group(1)) if m else 0.0

    fine = [r for r in rows if pip(r) >= DAPG_PIP]
    rest = [r for r in rows if pip(r) < DAPG_PIP]
    return {
        "fine_mapped": tally(fine, endpoint),
        "the_rest": tally(rest, endpoint),
        "prediction": E1_WIDE["secondary_fine_mapped_split"],
    }


def wide_summary(
    built: dict[str, Any],
    rows: list[dict[str, Any]],
    sign_check: dict[str, Any] | None = None,
    examples: int = 10,
) -> dict[str, Any]:
    """The committed plan of the widened E1: what the database holds, what is testable, and the budget."""
    calls = [mpra_call(r) for r in rows]
    usable = [c for c in calls if wide_usable(c)]
    sig = [c for c in usable if mpra_significant([c])]
    pairs = built["pairs"]
    per = {e: [p for p in pairs if p["endpoint"] == e] for e in WIDE_ENDPOINTS}
    return {
        "database": {
            "rows_on_the_main_assembly": len(calls),
            "distinct_variants": len({(c["chrom"], c["pos"], c["alt"]) for c in calls}),
            "significant_rows_any_cell": sum(
                bool(mpra_significant([c])) for c in calls if MPRA_EXCLUDE not in (c["study"] or "")
            ),
            "rows_in_a_cell_the_model_has": len(usable),
            "significant_rows_in_those_cells": len(sig),
            "significant_variants_in_those_cells": len({(c["chrom"], c["pos"], c["alt"]) for c in sig}),
            "cells": dict(Counter(c["cell"] for c in sig)),
            "studies": dict(Counter(c["study"][:60] for c in sig).most_common()),
            "chromosomes": dict(Counter(c["chrom"] for c in sig).most_common()),
        },
        "tests_built": built["tests"],
        "pairs": {e: len(v) for e, v in per.items()},
        "unmatched": built["unmatched"],
        "gene_source": built["gene_source"],
        "by_cell": {e: dict(Counter(p["borrowed"]["cell"] for p in v)) for e, v in per.items()},
        "by_chromosome": {
            e: dict(Counter(p["unit"].split(":")[0] for p in v).most_common()) for e, v in per.items()
        },
        "loci_20kb": {e: _loci(v) for e, v in per.items()},
        "control_distance_median_bp": {
            e: _median_abs([p["control"]["distance"] for p in v], 1) for e, v in per.items()
        },
        "control_distance_bands": {
            e: dict(Counter(_band(abs(p["control"]["distance"])) for p in v)) for e, v in per.items()
        },
        "measured_positive_share": {
            e: _share([p["borrowed"]["measured_sign"] > 0 for p in v]) for e, v in per.items()
        },
        "sign_convention": sign_check,
        "pre_registration": E1_WIDE,
        "budget": {
            "requests_per_pair": 2,
            "headline_pairs": min(1000, len(per["E1W_mpra_wide"])),
            "headline_requests": 2 * min(1000, len(per["E1W_mpra_wide"])),
            "saturation_pairs": min(300, len(per["E1S_mpra_saturation"])),
            "saturation_requests": 2 * min(300, len(per["E1S_mpra_saturation"])),
            "seconds_per_request_measured_elsewhere": 3.9,
        },
        "first_pairs": [
            {k: p[k] for k in ("order", "endpoint", "unit", "test", "control", "borrowed", "gene_how")}
            for p in pairs[:examples]
        ],
        "model_requests_spent": 0,
    }


def _band(d: int) -> str:
    return "under 5 kb" if d < 5_000 else "5 to 50 kb" if d < 50_000 else "50 to 250 kb"


def _loci(pairs: list[dict[str, Any]], window: int = 20_000) -> int:
    """How many separate neighbourhoods the pairs sit in: tested variants merged within 20 kb."""
    pos = sorted((p["unit"].split(":")[0], int(p["unit"].split(":")[1])) for p in pairs)
    n = 0
    last: tuple[str, int] | None = None
    for c, x in pos:
        if last is None or last[0] != c or x - last[1] > window:
            n += 1
        last = (c, x)
    return n


def plan(pairs: list[dict[str, Any]]) -> dict[str, Any]:
    """The request budget and the order, so the run only spends when the quota arrives."""
    by_endpoint = Counter(p["endpoint"] for p in pairs)
    looks = CRITERION["alpha_spending"]
    return {
        "requests_per_pair": 2,
        "pairs": len(pairs),
        "requests_if_all_run": 2 * len(pairs),
        "pairs_by_endpoint": dict(by_endpoint),
        "first_look_after_requests": 2 * looks["look_1_pairs"],
        "second_look_after_requests": 2 * looks["look_2_pairs"],
        "order": CRITERION["order"],
        "seconds_per_request_measured_elsewhere": 3.9,
        "hours_if_all_run_sequentially": round(2 * len(pairs) * 3.9 / 3600, 1),
        "cache": str(PREDICTIONS),
        "note": "a cached variant costs no request; the run reports requests spent against units answered",
    }


def track_over_units(key: str, chrom: str) -> list[dict[str, Any]]:
    """A bigBed track's rows over the chromosome's storage units; eQTL rows placed on their variant."""
    from genomeos.attribution.human_panel import _eqtl_pos

    rows, _ = track_rows(key, chrom, [(u["start"], u["end"]) for u in load_catalogue(chrom)])
    if key == "gtex_dapg":
        rows = [dict(r, chromStart=p, chromEnd=p + 1) for r in rows for p in [_eqtl_pos(r)] if p is not None]
    return rows


def assemble(chroms: list[str], progress=None) -> dict[str, Any]:
    """Part one, no model call: shortlist, tests, matched controls, the convention check and the plan."""
    say = progress or (lambda *a: None)
    if not list(GTEX_DIR.glob("hits_*.tsv")):
        say("streaming GTEx v8 significant pairs over the storage units")
        distil_gtex(chroms, progress=say)
    gtex = gtex_by_variant()
    symbols = gene_ids(chroms)
    mpra = {c: track_over_units("mpravardb", c) for c in chroms}
    dapg = {c: track_over_units("gtex_dapg", c) for c in chroms}
    say("shortlist")
    shortlist, pool = build_shortlist(chroms, gtex, symbols, mpra, dapg)
    tests = []
    for u in shortlist:
        t = tested_variant(u)
        if t:
            tests.append({"unit": u, "test": t, "order_key": order_key(u, t)})
    say("matched controls")
    matched = match_controls(tests, pool, random.Random(SEED))
    return {"shortlist": shortlist, "pool": pool, "tests": tests, "matched": matched}


def _share(xs: list[bool]) -> float | None:
    return round(sum(xs) / len(xs), 3) if xs else None


def summarise(assembled: dict[str, Any], examples: int = 12) -> dict[str, Any]:
    """The committed summary: counts, the pre-registration, the plan and the first pairs in order."""
    sl, tests, matched = assembled["shortlist"], assembled["tests"], assembled["matched"]
    pairs = matched["pairs"]
    events_ = [e for u in sl for e in u["measured"]]
    by_chrom: dict[str, Any] = {}
    for c in sorted({u["chrom"] for u in sl}):
        us = [u for u in sl if u["chrom"] == c]
        ts = [t for t in tests if t["unit"]["chrom"] == c]
        by_chrom[c] = {
            "units_with_a_measured_value": len(us),
            "tests": dict(Counter(t["test"]["endpoint"] for t in ts)),
            "pairs": dict(Counter(p["endpoint"] for p in pairs if p["unit"].startswith(c + ":"))),
        }
    return {
        "chromosomes": sorted(by_chrom),
        "by_chromosome": by_chrom,
        "measured_value_events": len(events_),
        "events_with_gtex_pairs": sum(bool(e["eqtl"]) for e in events_),
        "events_with_dapg_pip_05": sum(any((g["pip"] or 0) >= DAPG_PIP for g in e["dapg"]) for e in events_),
        "events_with_mpra_rows": sum(bool(e["mpra"]) for e in events_),
        "events_with_significant_mpra": sum(bool(mpra_significant(e["mpra"])) for e in events_),
        "tests": dict(Counter(t["test"]["endpoint"] for t in tests)),
        "pairs": dict(Counter(p["endpoint"] for p in pairs)),
        "unmatched": matched["unmatched"],
        "control_pool": len(assembled["pool"]),
        "match_levels": {
            e: dict(Counter(str(p.get("match_level")) for p in pairs if p["endpoint"] == e))
            for e in ENDPOINTS
        },
        "measured_positive_share": {
            e: _share([p["borrowed"]["measured_sign"] > 0 for p in pairs if p["endpoint"] == e])
            for e in ENDPOINTS
        },
        "control_distance_median_kb": _median_abs([p["control"]["distance"] for p in pairs], 1000),
        "control_r2_max": max((p["control"]["r2_with_test"] or 0 for p in pairs), default=None),
        "criterion": CRITERION,
        "amendments": list(AMENDMENTS),
        "hold_out": HOLD_OUT,
        "plan": plan(pairs),
        "first_pairs": [
            {k: p[k] for k in ("order", "endpoint", "unit", "test", "control", "borrowed")}
            for p in pairs[:examples]
        ],
        "sample_sizes": {
            "panel": "89 haplotype assemblies per value",
            "gnomad": "gnomAD v4.1.1 genomes, allele number per variant (about 152,000)",
            "gtex": "GTEx v8 significant pairs per tissue (73 to 706 donors by tissue)",
            "mpra": "MPRAVarDB as reported per study and cell line",
        },
        "model_requests_spent": 0,
    }


def _median_abs(xs: list[int], scale: int) -> float | None:
    v = sorted(abs(x) / scale for x in xs)
    return round(v[len(v) // 2], 1) if v else None
