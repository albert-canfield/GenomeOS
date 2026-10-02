# SPDX-License-Identifier: AGPL-3.0-or-later
"""What a fourth response-map increment could honestly add, counted before anything was built.

    uv run python scripts/response_map_increment4_count.py

The answer is NOTHING, and the numbers are the deliverable. Increment 3's committed count names
four candidates for a further kind of evidence. This counts each one at its own denominator, asks
of each the only question that decides it -- whether it can attach to a NATIVE LOCUS under the rule
increment 2 registered and increment 3 carried word for word -- and refuses all five arms they
divide into. No increment is built, nothing is added to the map, nothing is scored or fitted, and
no cut-off is introduced: every rule, span and call is imported from where it was registered.

Why a refusal is committed rather than left in a note. A refusal that lives only on a work board
gets re-litigated from scratch by someone who does not know the candidates were already
enumerated. The numbers below are what stops that: 227 / 35 / 3 further loci with ZERO native-locus
observations, 1,444 of 1,505 compiled elements UNASSESSED for eQTLs, 19 addable loci of which 6 are
already built, and 264 widening loci at which no three-input chain can exist.

It also carries, in a field of its own rather than buried in prose, THE NUMBER THIS COUNT WITHDREW
BEFORE QUOTING IT: the grouping offers 92 addable loci "reached" by an eQTL, and 86 of those 92 are
joined to one only by the 1 Mb proximity arm -- a DIFFERENT element within 1 Mb. Under increment
3's own attachment rule the figure is 6. See WITHDRAWN.

A NOTE FOR THE NEXT WRITER ON THIS PATH, and it is a trap the project has already paid for twice:
`response_map2.candidates` opens the 24 per-chromosome element tables through a result file's
`elements_where` POINTER field, and `data/results/eqtl_targets.json` carries `cache_target_effect`
values whose `window_source` is the per-element response cache. Both must be declared or
`save_result` refuses the result. 17 published results read a table they never declared.
"""

from __future__ import annotations

import collections
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import response_map_increment3_count as c3  # noqa: E402

from genomeos import manifest as mf  # noqa: E402
from genomeos import response_map2 as rm2  # noqa: E402
from genomeos.attribution import cell2, eqtl  # noqa: E402
from genomeos.attribution import context_evidence as ce  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.genome import reader  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

RESULT = "response_map_increment4_count"
CHROMS = rm2.CHROMS
ENTRY = "scripts/response_map_increment4_count.py"
OWN_CODE = (
    "scripts/response_map_increment4_count.py",
    "scripts/eqtl_crispri_frame.py",
    "tests/test_response_map_increment4.py",
)

#: The count this one starts from, and every helper taken from it rather than written again, so one
#: implementation answers "is there a retained hit here" and "could one have been seen here".
INCREMENT_3_COUNT = "data/results/response_map_increment3_count.json"

#: Imported, never chosen here. A count that invents its own rule cannot refuse anything.
IMPORTED = (
    "genomeos.response_map2.candidates - the three-input status of every CRISPRi pair, in the "
    "counting order scripts/response_map_coverage.py grouped in",
    "genomeos.response_map2.map1_keys, .STATUS_ORDER and .STATUS_ADDABLE - the best-status-per-locus rule",
    "genomeos.attribution.cell2.group - the independent-locus grouping, imported where registered",
    "scripts/response_map_increment3_count.py eqtl_index, eqtl_hits_by_chrom, covered_by and hits_in "
    "- increment 3's own frame and attachment tests, imported from its file rather than rewritten",
    "scripts/response_map_increment3_count.EQTL_MARGIN - scripts/eqtl_targets.MARGIN = 500, the "
    "margin the retained hits were distilled under",
    "genomeos.attribution.eqtl.KNOWLEDGE and load_hits' own column order - the retained GTEx hits",
)

#: Increment 2's rule, carried word for word through increment 3 and carried here unchanged. It is
#: the rule that refuses three of the four candidates, and it is quoted rather than paraphrased
#: because paraphrasing it is how it would be relaxed.
NATIVE_LOCUS_RULE = c3.NATIVE_LOCUS_RULE
REPORTER_WORDING = (
    "they read a sequence, or the bases inside it, not the native locus after a perturbation of it"
)

#: What increment 3 did with lentiMPRA, so that admitting it here could not be mistaken for a
#: precedent. It is the distinction the bigger number would have collapsed.
INCREMENT_3_PRECEDENT = (
    "increment 3 admitted lentiMPRA as `tile-in-cell` observations, labelled as such, AT LOCI THAT "
    "ALREADY CARRIED a perturbation and the reader's own reading, and still reported ZERO "
    "native-locus observations from it. A locus whose only evidence is a reporter reading carries no "
    "native-locus observation at all, so building such loci is not a labelled sequence arm: it is "
    "the map's population redefined to the loci the native-locus rule excludes"
)

#: Reported in its own field because catching it is the point. A count of loci joined by a grouping
#: is not a count of loci carrying the evidence.
WITHDRAWN = {
    "the_number": 92,
    "what_it_would_have_been_quoted_as": (
        "addable loci carrying two native-locus kinds - a CRISPRi perturbation and a GTEx cis-eQTL"
    ),
    "the_honest_figure": 6,
    "overstatement": "fifteen-fold",
    "why_it_is_wrong": (
        "`cell2.group` unions two keys on a shared gene or on 1 Mb proximity on one chromosome. An "
        "eQTL element 1 Mb from an addable locus therefore lands in that locus's group without the "
        "variant being anywhere inside the element the chain stands on. 86 of the 92 are joined only "
        "by that proximity arm, so they are a DIFFERENT element within 1 Mb and not a second "
        "observation of the same locus"
    ),
    "the_rule_that_gives_6": (
        "increment 3's own attachment rule, in its registered wording: 'a retained hit whose "
        "variant position falls inside the compiled element widened by 500 bases each side, which "
        "is the margin the retained set was distilled under (scripts/eqtl_targets.MARGIN), imported "
        "and not chosen here'. The element it must fall in is one the chain stands on. Applied to "
        "the 166 shown elements it gives 6 loci, which is what increment 3 built"
    ),
    "the_window_is_not_what_separates_92_from_6": (
        "both figures use the same 500-base window, because `hits_in` widens the element by the "
        "distil margin itself. What separates them is WHICH ELEMENT the hit must fall in - one the "
        "chain stands on, or any element at the locus - so the overstatement is about attachment "
        "and not about a tolerance"
    ),
    "the_shape_it_is": (
        "the same shape as the 198 loci quoted as 198 places the model is wrong when 66 held a "
        "measured decrease, and as the 48 measured links that cleared a floor while only 21 were "
        "answerable in their own cell. A count is not a measurement of the thing wanted"
    ),
    "caught": "in this lane's own working, before the figure was quoted anywhere",
}

#: The trap, carried into the result so the next writer on this path reads it from the file.
POINTER_FIELD_NOTE = (
    "`response_map2.candidates` opens the 24 per-chromosome element tables through a result file's "
    "`elements_where` POINTER field - about 1.4 GB that no reading of the source can name - and "
    "`data/results/eqtl_targets.json` carries `cache_target_effect` values whose `window_source` is "
    "the per-element response cache. Both must be declared or `save_result` refuses the result. 17 "
    "published results read a table they never declared, which is why inputs are TRACED and not "
    "believed"
)

#: Files opened that no glob of data/results reaches. Carried from increment 3's own declaration,
#: because this count runs the same two readers over the same chromosomes.
OPENED_BESIDE_THE_RESULTS = c3.OPENED_BESIDE_THE_RESULTS
GENCODE = c3.GENCODE


class RefusedError(ValueError):
    """Raised instead of returning a payload whose reconciliation does not hold."""


def arms(root: Path = ROOT) -> dict[str, Any]:
    """Count every candidate a fourth increment could be built on, and what each one can attach to."""
    cands, counted = rm2.candidates(list(CHROMS), RESULTS_DIR, ce.Readers(), ce.mapping())
    covered = rm2.map1_keys(root)
    crispri_keys = [c.key for c in cands]
    groups = cell2.group([*crispri_keys, *covered])
    mine = groups[: len(cands)]
    best: dict[int, str] = {}
    for g, c in zip(mine, cands, strict=True):
        if g not in best or rm2.STATUS_ORDER[c.status] < rm2.STATUS_ORDER[best[g]]:
            best[g] = c.status
    addable = {g for g, s in best.items() if s == rm2.STATUS_ADDABLE}

    # the elements the built chains stand on (increment 3's `shown`), and every element any CRISPRi
    # pair attaches to, kept APART: the difference between them is what the withdrawn number was
    shown: dict[int, set[str]] = collections.defaultdict(set)
    shown_iv: dict[str, tuple[str, int, int]] = {}
    att_iv: dict[str, tuple[str, int, int]] = {}
    att_groups: dict[str, set[int]] = collections.defaultdict(set)
    group_pairs: dict[int, list[Any]] = collections.defaultdict(list)
    for g, c in zip(mine, cands, strict=True):
        group_pairs[g].append(c)
        if c.element is None:
            continue
        eid = c.element["id"]
        att_iv[eid] = (c.pair.chrom, c.element["start"], c.element["end"])
        att_groups[eid].add(g)
        if g in addable and c.status == rm2.STATUS_ADDABLE:
            shown[g].add(eid)
            shown_iv[eid] = (c.pair.chrom, c.element["start"], c.element["end"])

    # the measured layer's other three assays, as their own locus keys, the way increment 3 keyed them
    other_keys: dict[str, list[cell2.LocusKey]] = {a: [] for a in ms.ASSAYS}
    for chrom in CHROMS:
        elements = c3.cp._attributed(chrom, RESULTS_DIR)
        _layer, rows = c3.cp._measured_rows(chrom, elements, RESULTS_DIR)
        for r in rows:
            for a in r["assays"]:
                if a == "crispri":
                    continue
                other_keys[a].append(cell2.LocusKey(a, chrom, r["start"], r["end"], f"{a}|{r['id']}"))

    # the eQTL arm's keys, built exactly as increment 3 built them, from the committed result
    with open(root / c3.EQTL_RESULT) as fh:
        eq = json.load(fh)
    eqtl_keys: list[cell2.LocusKey] = []
    key_elem: list[str] = []
    elem_rows: dict[str, dict[str, Any]] = {}
    for setname, blk in eq["sets"].items():
        for r in blk["elements_with_eqtl"]:
            s, e = r["id"].split(":")[1].split("-")
            elem_rows[r["id"]] = {**r, "set": setname, "s": int(s), "e": int(e)}
            for gene in r["egenes"]:
                eqtl_keys.append(cell2.LocusKey("gtex_donors", r["chrom"], int(s), int(e), gene))
                key_elem.append(r["id"])

    # --- one grouping per candidate, against every CRISPRi key ---------------------------------
    further: dict[str, dict[str, Any]] = {}
    for kind, keys in [*other_keys.items(), ("eqtl", eqtl_keys)]:
        if kind == "crispri":
            continue
        g2 = cell2.group([*crispri_keys, *keys])
        theirs = g2[len(crispri_keys) :]
        reached = set(g2[: len(crispri_keys)])
        further[kind] = {
            "observations": len(keys),
            "independent_loci_of_this_kind": len(set(theirs)),
            "of_those_the_crispri_arm_already_reaches": len(set(theirs) & reached),
            "further_independent_loci": len(set(theirs) - reached),
        }
        if kind == "eqtl":
            eqtl_group_of = list(theirs)
            eqtl_reached = reached
            eqtl_cr_best: dict[int, str] = {}
            for g, c in zip(g2[: len(crispri_keys)], cands, strict=True):
                if g not in eqtl_cr_best or (rm2.STATUS_ORDER[c.status] < rm2.STATUS_ORDER[eqtl_cr_best[g]]):
                    eqtl_cr_best[g] = c.status

    # --- the eQTL arm, under increment 3's own frame and attachment tests ----------------------
    index, frame_n = c3.eqtl_index(root)
    hits = c3.eqtl_hits_by_chrom(eqtl.KNOWLEDGE)
    margin = c3.EQTL_MARGIN

    def test(iv: dict[str, tuple[str, int, int]]) -> tuple[set[str], dict[str, list[dict[str, Any]]]]:
        assessable: set[str] = set()
        hit: dict[str, list[dict[str, Any]]] = {}
        for eid, (chrom, s, e) in iv.items():
            if c3.covered_by(index.get(chrom, []), s - margin, e + margin):
                assessable.add(eid)
            got = c3.hits_in(hits.get(chrom, []), s, e)
            if got:
                hit[eid] = got
        return assessable, hit

    sh_ass, sh_hit = test(shown_iv)
    at_ass, at_hit = test(att_iv)
    hit_groups = {g for eid in at_hit for g in att_groups[eid]}
    ass_groups = {g for eid in at_ass for g in att_groups[eid]}
    shown_hit_loci = {g for g, eids in shown.items() if eids & set(sh_hit)}
    addable_hit_any = {g for g in hit_groups if best[g] == rm2.STATUS_ADDABLE}

    def outcomes(sel: set[int]) -> dict[str, Any]:
        pairs: collections.Counter[str] = collections.Counter()
        loci: collections.Counter[str] = collections.Counter()
        for g in sel:
            seen = set()
            for c in group_pairs[g]:
                pairs[c.pair.outcome] += 1
                seen.add(c.pair.outcome)
            for o in seen:
                loci[o] += 1
        return {
            "pairs": dict(pairs),
            "loci_with_at_least_one": dict(loci),
            "loci_rows_do_not_sum": (
                "a locus can carry pairs of more than one outcome, so these rows are counted "
                "independently and do not sum to the locus count"
            ),
        }

    slope_neg = sum(1 for v in at_hit.values() for h in v if h["slope"] < 0)
    slope_pos = sum(1 for v in at_hit.values() for h in v if h["slope"] >= 0)

    deepening = {
        "what_it_would_be": (
            "a SECOND native-locus kind at the loci the map is already built over, which is the only "
            "thing that makes the existing map deeper rather than wider"
        ),
        "crispri_attached_compiled_elements": len(att_iv),
        "inside_the_distillation_frame_at_all": len(at_ass),
        "outside_the_frame_unassessed": len(att_iv) - len(at_ass),
        "unassessed_is_not_an_absence": (
            "the retained set was distilled against the intervals of three sampled element sets, so "
            "at an element outside that frame no eQTL could have been seen whatever is really there. "
            "Its absence is UNASSESSED and is not an absence of eQTLs"
        ),
        "frame_intervals": frame_n,
        "elements_with_at_least_one_retained_hit": len(at_hit),
        "attachment_rule": (
            "increment 3's own, imported and not chosen: 'a retained hit whose variant position "
            "falls inside the compiled element widened by 500 bases each side, which is the margin "
            "the retained set was distilled under (scripts/eqtl_targets.MARGIN)'. `hits_in` applies "
            "that widening itself, so every figure on this arm is at that window"
        ),
        "retained_records_on_them": sum(len(v) for v in at_hit.values()),
        "slope_negative": slope_neg,
        "slope_positive": slope_pos,
        "tissues": len({h["tissue"] for v in at_hit.values() for h in v}),
        "crispri_loci_total": len(set(mine)),
        "loci_the_frame_could_have_seen_at_all": len(ass_groups),
        "loci_with_an_eqtl_inside_an_attached_element": len(hit_groups),
        "loci_with_a_hit_by_best_status": dict(collections.Counter(best[g] for g in hit_groups)),
        "loci_assessable_by_best_status": dict(collections.Counter(best[g] for g in ass_groups)),
        "addable_loci_with_a_hit_on_ANY_attached_element": len(addable_hit_any),
        "addable_loci_with_a_hit_on_a_SHOWN_element": len(shown_hit_loci),
        "already_built_by_increment_3": len(shown_hit_loci),
        "loci_this_arm_would_ADD": 0,
        "why_zero": (
            f"{len(addable_hit_any)} addable loci carry an eQTL inside SOME element a CRISPRi pair "
            f"attaches to, but only {len(shown_hit_loci)} carry it inside an element a built chain "
            "stands on - and those are the ones increment 3 built. The difference of "
            f"{len(addable_hit_any) - len(shown_hit_loci)} requires the attachment relaxed from the "
            "chain's element to any element at the locus. Increment 3 chose the shown elements "
            "deliberately: they are 'the elements the built chains actually stand on, which is where "
            "another kind of evidence has to attach to be evidence about the same chain'. Relaxing "
            "that is relaxing an imported rule, so this arm adds nothing"
        ),
        "what_the_crispri_side_says_at_the_hit_loci": outcomes(hit_groups),
        "what_the_crispri_side_says_at_the_hit_loci_that_are_not_addable": outcomes(
            {g for g in hit_groups if best[g] != rm2.STATUS_ADDABLE}
        ),
        "the_only_way_to_grow_it": (
            "re-distil GTEx v8 against the CRISPRi element set, which turns the 1,444 unassessed "
            "elements into an answer. Registered and run separately as `eqtl_crispri_frame`; it is "
            "not assumed here and no number here depends on it"
        ),
        "refused_under_cause": "adds_no_locus_under_the_rules_as_registered",
    }

    further_groups = set(eqtl_group_of) - eqtl_reached
    further_elems: set[str] = set()
    for g, eid in zip(eqtl_group_of, key_elem, strict=True):
        if g in further_groups:
            further_elems.add(eid)
    w_neg = w_pos = w_rec = 0
    w_tis: set[str] = set()
    w_genes: set[str] = set()
    for eid in further_elems:
        r = elem_rows[eid]
        for h in c3.hits_in(hits.get(r["chrom"], []), r["s"], r["e"]):
            w_rec += 1
            w_tis.add(h["tissue"])
            w_genes.add(h["gene_id"])
            if h["slope"] < 0:
                w_neg += 1
            else:
                w_pos += 1
    widening = {
        "what_it_would_be": (
            "the loci of the eQTL arm the CRISPRi arm does not reach, built as the map's own loci"
        ),
        "further_independent_loci": len(further_groups),
        "elements": len(further_elems),
        "retained_records": w_rec,
        "element_egene_keys": sum(len(elem_rows[e]["egenes"]) for e in further_elems),
        "distinct_egenes": len(w_genes),
        "tissues": len(w_tis),
        "slope_negative": w_neg,
        "slope_positive": w_pos,
        "by_sampled_set": dict(collections.Counter(elem_rows[e]["set"] for e in further_elems)),
        "the_model_named_a_predicted_target": sum(1 for e in further_elems if elem_rows[e].get("predicted")),
        "the_model_named_none": sum(1 for e in further_elems if not elem_rows[e].get("predicted")),
        "an_observation_of_the_native_locus": True,
        "native_locus_kinds_per_locus": 1,
        "perturbation_observations": 0,
        "why_no_perturbation": (
            "zero by construction, not by result: 'further' means the CRISPRi arm does not reach "
            "these loci, so the map's central relation - the locus perturbed and the consequence "
            "read out - is absent at every one of them"
        ),
        "reader_state_observations": 0,
        "why_no_reader_state": (
            "the eQTL locus key names no cell: its context is `gtex_donors`. Increment 3 took the "
            "reader's biosample from the CRISPRi pair's own cell, and with no cell named a reader "
            "state would mean CHOOSING a biosample - a new choice rather than an imported call, and "
            "GTEx tissues are not the reader's biosamples. So no three-input chain can exist here"
        ),
        "three_input_chains_possible": 0,
        "population": (
            "three sampled element sets (enhancer_targets_chr*, constrained_targets_chr*, "
            "vista_chr*), not the measured perturbation layer the map has been built over. This "
            "widens the frame rather than deepening the evidence"
        ),
        "refused_under_cause": "one_native_locus_kind_and_no_possible_chain",
    }

    return {
        "cands": cands,
        "counted": counted,
        "mine": mine,
        "best": best,
        "addable": addable,
        "shown_elements": len(shown_iv),
        "shown_elements_with_a_retained_hit": len(sh_hit),
        "shown_elements_inside_the_frame": len(sh_ass),
        "shown_hit_loci": shown_hit_loci,
        "shown_assessable_loci": {g for g, eids in shown.items() if eids & sh_ass},
        "further": further,
        "eqtl_groups": set(eqtl_group_of),
        "eqtl_reached": {g for g in set(eqtl_group_of) if g in eqtl_reached},
        "eqtl_further": further_groups,
        "eqtl_keys": len(eqtl_keys),
        "eqtl_elements": len(elem_rows),
        "eqtl_cr_best": eqtl_cr_best,
        "deepening": deepening,
        "widening": widening,
        "other_keys": other_keys,
    }


def candidate_blocks(a: dict[str, Any]) -> list[dict[str, Any]]:
    """One block per candidate increment 3's count names, with the question that decides it."""
    reporters = {
        "lentimpra": (
            "reporter activity: a 200 bp copy driving an integrated lentiviral reporter, outside "
            "its native locus (Agarwal et al. 2025, Nature)"
        ),
        "vista": "developmental activity: a transgenic reporter in the mouse embryo at e11.5",
        "satmut": "base-level sensitivity: single substitutions of the element read out in a reporter",
    }
    out = []
    for kind, what in reporters.items():
        f = a["further"][kind]
        out.append(
            {
                "candidate": kind,
                "what_the_observation_is": what,
                "further_independent_loci": f["further_independent_loci"],
                "observations": f["observations"],
                "independent_loci_of_this_kind": f["independent_loci_of_this_kind"],
                "of_those_the_crispri_arm_already_reaches": f["of_those_the_crispri_arm_already_reaches"],
                "an_observation_of_the_native_locus": False,
                "native_locus_observations": 0,
                "why_not": REPORTER_WORDING,
                "increment_3_precedent": INCREMENT_3_PRECEDENT,
                "verdict": "REFUSED",
                "refused_under_cause": "zero_native_locus_observations",
                "what_admitting_it_would_cost": (
                    "the native-locus rule relaxed. That is the only way these loci enter the map, "
                    "and it is why they are refused rather than labelled"
                ),
                "upper_bound": (
                    "the further-locus figure is an UPPER BOUND besides: a reporter observation names "
                    "no measured gene, so only the 1 Mb proximity arm of the grouping can fire for it "
                    "and with gene identity the groups could only be fewer"
                ),
            }
        )
    f = a["further"]["eqtl"]
    out.append(
        {
            "candidate": "eqtl",
            "what_the_observation_is": (
                "an association, not a perturbation: a significant cis-eQTL is a variant inside the "
                "element whose alleles go with a gene's expression across GTEx donors, in one tissue"
            ),
            "source": eqtl.EVIDENCE,
            "an_observation_of_the_native_locus": True,
            "why": (
                "the variant sits on the chromosome in place and the expression is of the donor's "
                "own gene, so the reading is of the native locus. It is natural variation and an "
                "association, never an experimental perturbation, and association is not mechanism"
            ),
            "further_independent_loci": f["further_independent_loci"],
            "observations": f["observations"],
            "observations_unit": (
                "element-eGene keys, not variant-gene-tissue records and emphatically not places in "
                "the genome"
            ),
            "independent_loci_of_this_kind": f["independent_loci_of_this_kind"],
            "of_those_the_crispri_arm_already_reaches": f["of_those_the_crispri_arm_already_reaches"],
            "it_splits_into_two_arms_and_both_are_refused": {
                "deepening": a["deepening"],
                "widening": a["widening"],
            },
            "verdict": "REFUSED on both arms",
        }
    )
    return out


def reconcile(a: dict[str, Any], blocks: list[dict[str, Any]]) -> dict[str, Any]:
    """Every candidate accounted for, and every locus of the two populations. Refuses otherwise."""
    causes: collections.Counter[str] = collections.Counter()
    for b in blocks:
        if b["candidate"] == "eqtl":
            causes[a["deepening"]["refused_under_cause"]] += 1
            causes[a["widening"]["refused_under_cause"]] += 1
        else:
            causes[b["refused_under_cause"]] += 1
    arms_named = 5  # three reporter candidates, and the eQTL candidate's two arms
    rec = {
        "candidates_increment_3_names": len(blocks),
        "arms_they_divide_into": arms_named,
        "arms_built": 0,
        "arms_refused_by_cause": dict(causes),
        "arms_accounted_for": sum(causes.values()),
        "rule": (
            "arms_built + sum(arms_refused_by_cause) == arms_they_divide_into, and no arm may leave "
            "the count without a named cause. `main` refuses to save a payload where this does not "
            "hold"
        ),
        "crispri_loci": len(set(a["mine"])),
        "crispri_loci_by_best_status": dict(collections.Counter(a["best"].values())),
        "crispri_loci_accounted_for": sum(collections.Counter(a["best"].values()).values()),
        "eqtl_loci_of_this_kind": len(a["eqtl_groups"]),
        "eqtl_loci_reached_by_the_crispri_arm": len(a["eqtl_reached"]),
        "eqtl_loci_further": len(a["eqtl_further"]),
        "eqtl_loci_accounted_for": len(a["eqtl_reached"]) + len(a["eqtl_further"]),
    }
    rec["holds"] = (
        rec["arms_built"] + rec["arms_accounted_for"] == arms_named
        and rec["crispri_loci_accounted_for"] == rec["crispri_loci"]
        and rec["eqtl_loci_accounted_for"] == rec["eqtl_loci_of_this_kind"]
    )
    return rec


def anchors(a: dict[str, Any], root: Path = ROOT) -> dict[str, Any]:
    """Every anchor of the committed count, re-derived here, against what that file says.

    This is what makes the candidate figures above trustworthy rather than merely stated: the same
    code path reproduces thirteen committed quantities before it is asked for a new one.
    """
    with open(root / INCREMENT_3_COUNT) as fh:
        c = json.load(fh)
    k = {x["kind"]: x for x in c["kinds_of_evidence"]}
    d = c["denominator"]
    mine = {
        "pairs": a["counted"]["pairs"],
        "crispri_loci": len(set(a["mine"])),
        "addable_loci": len(a["addable"]),
        "shown_elements": a["shown_elements"],
        "shown_elements_with_a_retained_hit": a["shown_elements_with_a_retained_hit"],
        "shown_elements_inside_the_frame": a["shown_elements_inside_the_frame"],
        "addable_loci_with_an_eqtl_on_a_shown_element": len(a["shown_hit_loci"]),
        "addable_loci_the_frame_could_see": len(a["shown_assessable_loci"]),
        "eqtl_loci_of_this_kind": a["further"]["eqtl"]["independent_loci_of_this_kind"],
        "eqtl_loci_reached": a["further"]["eqtl"]["of_those_the_crispri_arm_already_reaches"],
        "eqtl_loci_further": a["further"]["eqtl"]["further_independent_loci"],
        "eqtl_element_egene_keys": a["eqtl_keys"],
        "eqtl_elements": a["eqtl_elements"],
    }
    theirs = {
        "pairs": d["measured_perturbation_pairs"],
        "crispri_loci": d["independent_loci_of_the_measured_layer"],
        "addable_loci": d["addable_loci_re_derived_here"],
        "shown_elements": d["shown_elements_of_those_loci"],
        "shown_elements_with_a_retained_hit": k["eqtl"]["shown_elements_with_a_retained_hit"],
        "shown_elements_inside_the_frame": k["eqtl"][
            "shown_elements_the_retained_set_could_have_seen_at_all"
        ],
        "addable_loci_with_an_eqtl_on_a_shown_element": k["eqtl"]["loci_on_the_existing_110"],
        "addable_loci_the_frame_could_see": k["eqtl"]["loci_the_retained_set_could_have_seen_at_all"],
        "eqtl_loci_of_this_kind": k["eqtl"]["further"]["independent_loci_of_this_kind"],
        "eqtl_loci_reached": k["eqtl"]["further"]["of_those_the_crispri_arm_already_reaches"],
        "eqtl_loci_further": k["eqtl"]["further"]["further_independent_loci"],
        "eqtl_element_egene_keys": k["eqtl"]["further"]["observations"],
        "eqtl_elements": k["eqtl"]["further_elements_with_an_eqtl_in_the_committed_result"],
    }
    differ = {q: [mine[q], theirs[q]] for q in mine if mine[q] != theirs[q]}
    return {
        "counted_at": INCREMENT_3_COUNT,
        "counted_sha256": mf.input_entry(root / INCREMENT_3_COUNT)["sha256"],
        "re_derived_here": mine,
        "as_the_committed_count_states_them": theirs,
        "anchors_compared": len(mine),
        "anchors_differing": differ,
        "all_reproduce": not differ,
        "why_this_is_here": (
            "a candidate figure is trustworthy when the code that produced it first reproduces the "
            "committed quantities of the same population. Thirteen are compared and named, so a "
            "reader can see which were checked rather than take the agreement on trust"
        ),
    }


def main() -> int:
    a = arms()
    blocks = candidate_blocks(a)
    # the reconciliation is checked BEFORE anything else is computed, so a count that does not
    # account for every arm cannot reach the anchor check, the manifest or the writer
    rec = reconcile(a, blocks)
    if not rec["holds"]:
        raise RefusedError(f"reconciliation does not hold, no payload is served: {rec}")
    anc = anchors(a)
    if not anc["all_reproduce"]:
        raise RefusedError(
            "the committed count's anchors do not reproduce, so no figure here is comparable "
            f"with it: {anc['anchors_differing']}"
        )

    payload: dict[str, Any] = {
        "result": RESULT,
        "date": "2026-10-02",
        "question": (
            "what a fourth response-map increment could honestly add, counted over every candidate "
            "increment 3's committed count names, before anything is built"
        ),
        "verdict": "NOTHING. Every candidate is refused, and the numbers are the deliverable",
        "counting_only": (
            "no increment is built, nothing is added to the map, nothing is scored and nothing is "
            "fitted. No cut-off is introduced: every rule, span and call is imported from where it "
            "was registered. 0 model requests, no downloads, no money: every input was on disk"
        ),
        "why_a_refusal_is_committed": (
            "a refusal of this kind that lives only on a work board is re-litigated from scratch by "
            "someone who does not know the candidates were already enumerated. The numbers are what "
            "stops that, and they belong in a file with a manifest"
        ),
        "native_locus_rule": NATIVE_LOCUS_RULE,
        "reporter_wording_carried_word_for_word": REPORTER_WORDING,
        "increment_3_precedent_on_lentimpra": INCREMENT_3_PRECEDENT,
        "imported_not_chosen": list(IMPORTED),
        "candidates": blocks,
        "the_number_withdrawn_before_it_was_quoted": WITHDRAWN,
        "note_for_the_next_writer_on_this_path": POINTER_FIELD_NOTE,
        "anchors": anc,
        "reconciliation": rec,
        "cannot_establish": [
            "nothing here is a verdict about biology. It is a count of what evidence exists and of "
            "where it stops, and a refusal to build on it",
            "that the reporter assays' loci hold no regulatory function. They hold no NATIVE-LOCUS "
            "observation, which is a statement about what their assays read and not about the "
            "sequence",
            "that the 1,444 unassessed elements carry no eQTL. No eQTL could have been seen at them "
            "at all, and unassessed is not absent",
            "that the 25 loci carrying both a perturbation and an eQTL are 25 agreements. A second "
            "kind of evidence is not confirmation of the first, and GTEx tissues are not the CRISPRi "
            "cell lines",
            "that the widening arm's 264 loci are 264 places in the genome where anything is "
            "established. They are loci of three sampled element sets carrying one association each",
            "the independent-locus grouping is an operational grouping that pools cell types, NOT "
            "established biological independence, and nothing may cite it as one",
        ],
        "alphagenome_requests": 0,
        "money": "none: every input was already on disk",
    }

    inputs = [mf.input_entry(ce.TRACK_METADATA, partition=None)]
    for f in ms.CRISPRI_FILES:
        p = ms.CRISPRI_KNOWLEDGE / f
        if p.exists():
            inputs.append(mf.input_entry(p, partition=ms.CRISPRI_SPLIT_OF[f]))
    seen = {str(i["path"]) for i in inputs}
    paths = [
        p
        for cell in ce.READER_TERMS
        for chrom in CHROMS
        if (p := RESULTS_DIR / reader.peaks_path(cell, chrom).name).exists()
    ]
    for g in (
        "budget_chr*.json",
        "budget_axes_chr*.json",
        "variation_chr*.json",
        "duplication_chr*.json",
        "domains_chr*.json",
        "unknown_chr*.json",
        "reader_*_chr*.json",
        "ccres_chr*.bed.gz",
        "enhancer_targets_chr*.json",
        "enhancer_targets_all_chr*.json",
        "constrained_targets_chr*.json",
        "vista_chr*.json",
    ):
        for chrom in CHROMS:
            paths += sorted(RESULTS_DIR.glob(g.replace("chr*", chrom)))
    paths.append(Path(c3.EQTL_RESULT))
    paths.append(Path(INCREMENT_3_COUNT))
    paths += sorted(eqtl.KNOWLEDGE.glob("hits_*.tsv"))
    p = eqtl.KNOWLEDGE / "distil_summary.json"
    if p.exists():
        paths.append(p)
    for group in OPENED_BESIDE_THE_RESULTS.values():
        paths += [Path(x) for x in group]
    for chrom in CHROMS:
        q = Path(GENCODE.format(chrom=chrom))
        if q.exists():
            paths.append(q)
    for q in paths:
        if q.as_posix() not in seen:
            seen.add(q.as_posix())
            inputs.append(mf.input_entry(q, partition=None))

    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "EngreitzLab/CRISPR_comparison (Gschwind et al. 2025)",
                "version": "the two local benchmark tables under data/knowledge/crispri",
            },
            {
                "accession": "GTEx v8 single-tissue cis-eQTLs, significant variant-gene pairs",
                "version": (
                    "the local distillation under data/knowledge/gtex, kept only inside the three "
                    "sampled element sets; see data/knowledge/gtex/distil_summary.json"
                ),
            },
            {
                "accession": "the response map's third increment and its count",
                "version": "data/results/response_map_increment3_count.json, by sha256",
            },
        ],
        "inputs": inputs,
        "inputs_opened_beside_the_declared_results": {
            "why_they_are_listed": (
                "an input is what the writer opened, not what it said it read. The audit hook in "
                "`genomeos.manifest` records every open under `data/` and `save_result` refuses a "
                "result whose traced reads exceed its declared inputs"
            ),
            "groups": {k: len(v) for k, v in OPENED_BESIDE_THE_RESULTS.items()},
            "the_pointer_case": POINTER_FIELD_NOTE,
        },
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "independent_locus_span": cell2.INDEPENDENT_LOCUS_SPAN,
            "reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
            "openness_call": ce.OPENNESS_CALL,
            "eqtl_margin_bp": c3.EQTL_MARGIN,
            "chromosomes": list(CHROMS),
        },
        "exclusions": [
            "no pair is excluded: the denominator is every CRISPRi pair the measured layer parses as "
            "valid on these chromosomes, both arms and every cell type",
            "lentiMPRA, VISTA and saturation mutagenesis are refused as candidate populations under "
            "a named cause - they read a sequence, not the native locus after a perturbation - and "
            "their further-locus counts are reported rather than dropped",
            "an eQTL outside the distillation's sampled frame is excluded as unassessed, by name, "
            "and is not counted as an absence of eQTLs",
        ],
        "partitions": {
            "candidates": "one block per candidate increment 3's count names",
            "arms": "the eQTL candidate's deepening and widening arms, refused separately",
        },
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }

    path = save_result(RESULT, payload)
    print(f"{RESULT}: {path}")
    print(f"  verdict: {payload['verdict']}")
    for b in blocks:
        print(
            f"  {b['candidate']}: further loci {b['further_independent_loci']}, "
            f"native-locus {b['an_observation_of_the_native_locus']} -> {b['verdict']}"
        )
    d = a["deepening"]
    print(
        f"  deepening: {d['outside_the_frame_unassessed']} of "
        f"{d['crispri_attached_compiled_elements']} elements unassessed; "
        f"{d['addable_loci_with_a_hit_on_ANY_attached_element']} addable loci with a hit on any "
        f"element, {d['addable_loci_with_a_hit_on_a_SHOWN_element']} on a shown one; "
        f"adds {d['loci_this_arm_would_ADD']}"
    )
    print(
        f"  withdrawn: {WITHDRAWN['the_number']} -> {WITHDRAWN['the_honest_figure']} "
        f"({WITHDRAWN['overstatement']})"
    )
    print(f"  anchors: {anc['anchors_compared']} compared, all reproduce {anc['all_reproduce']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
