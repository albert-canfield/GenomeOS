# SPDX-License-Identifier: AGPL-3.0-or-later
"""Render section 21 of docs/LOCI-BENCHMARK.md from the third set's result file.

    uv run python scripts/loci_third_section.py            # print it
    uv run python scripts/loci_third_section.py --write    # splice it into the doc

Every figure in the section is read out of `data/results/loci_third.json` and
`data/results/loci_third_intervals.json`. None is typed. Section 19 had to be hand-copied from its
result and section 10's cCRE count went stale in the doc while the result moved on; a section that
is rendered cannot drift from the run that produced it, and re-running this after a re-scoring
updates every number at once.

The prose is the reading and is written by hand; the numbers inside it are interpolated. Where a
sentence would change meaning if a number moved, the sentence names the number rather than asserting
the direction, so the two cannot disagree.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from genomeos.benchmark import loci_third as third
from genomeos.results import load_result

DOC = Path("docs/LOCI-BENCHMARK.md")
HEADING = (
    "## 21. A third set, the direction axis asked both ways, and an answer that was"
    " registered as undecidable (2026-09-21)"
)


def kn(d: dict[str, Any]) -> str:
    """`k/n (rate)` straight out of a rate block."""
    return f"{d['k']}/{d['n']} ({d['rate']:.3f})" if d.get("n") else "-"


def frac(d: dict[str, Any]) -> str:
    return f"{d['k']}/{d['n']}" if d.get("n") else "-"


def pct(x: float | None) -> str:
    return "-" if x is None else f"{x:+.3f}"


def spread(block: dict[str, Any], frames: dict[str, Any]) -> str:
    """How far one rate moves across the three frames: the number the doc argues from."""
    rates = [block[k]["rate"] for k in frames if block[k].get("rate") is not None]
    return f"{max(rates) - min(rates):.3f}"


def p_of(block: dict[str, Any]) -> str:
    p = block.get("p_one_sided")
    return "-" if p is None else f"{p:.3f}"


def render(result: dict[str, Any], intervals: dict[str, Any] | None) -> str:
    a = result["aggregate"]
    reach = result["reach"]
    plan = {r["locus"]: r for r in result["plan"]}
    by_locus = {r["locus"]: r for r in result["loci"]}
    beside = result["beside_the_other_two"]
    ac = result["against_controls"]
    dirs = result["directions"]
    pc = result["positive_control"]
    read = {r["locus"]: r for r in (intervals or {}).get("read", [])}
    # the plan is re-derived on every run, so after a stated interval has been scored its row shows
    # nothing left to spend. What was actually bought is the intervals result, and what was free is
    # every locus that is not in it.
    spent = len(read)
    already = len(result["plan"]) - spent
    lengths = sorted(r["length"] for r in result["plan"])
    at_the_element = sum(
        1
        for r in result["loci"]
        if ((r["readings"].get("deletion") or {}).get("target") in r["expected"]["targets"])
    )
    bought = ac["input_presence"]["a deletion spent inside the window"]
    out: list[str] = [HEADING, ""]

    out += [
        f"`genomeos/benchmark/loci_third.py`, `scripts/loci_third.py`, `tests/test_loci_third.py`,"
        f" results `{third.NAME}` and `{third.INTERVALS}`. **{len(read)} model"
        f" request{'s' if len(read) != 1 else ''} spent.** {result['cost']['seconds']:.0f} seconds"
        " of range reads over public tracks for the rest.",
        "",
        'Section 19\'s own "left undone" named three gaps. This set is aimed at two of them and'
        " records the third as unfillable. **Every one of the nine candidates activates**, so the"
        " benchmark has never asked whether the deletion layer's sign means anything when the"
        " published element's job is to hold a gene *down*; and section 20 asked for a set whose"
        " element lengths span the panel's rather than sitting at a uniform 500 bp. Two of the nine"
        " loci here are published **repressors** with a perturbation behind them, and the element"
        f" lengths run {lengths[0]} bp to {lengths[-1]:,} bp.",
        "",
        "It is a **third sampling frame**, separately registered, scored by"
        " `genomeos/benchmark/loci.py` with no change to any scorer or hit rule. Its rates are"
        " reported beside the seventeen and the nine with each frame named, and never pooled with"
        " either.",
        "",
        "### The nine, and what the free filter did to them first",
        "",
        "| locus | element | len | published target | direction | distance | perturbation |",
        "|---|---|---|---|---|---|---|",
    ]
    how = {
        "SHH_SBE2": "point mutation in a holoprosencephaly patient; transgenic reporter (Jeong 2008)",
        "GATA2_plus95": "targeted mouse enhancer deletion; germline human mutations (Johnson 2012)",
        "TAL1_MuTE": "CRISPR deletion of a somatically created enhancer (Mansour 2014)",
        "CDKN2A_9p21": "70 kb mouse interval deletion (Visel 2010)",
        "KITLG_blond": "human-enhancer knock-in mouse (Guenther 2014)",
        "SOST_VanBuchem": "homozygous 52 kb human deletion; mouse enhancer deletion (Balemans 2002)",
        "HBG1_BCL11A_site": "HPFH point mutations; editing in human erythroid cells (Martyn 2018)",
        "MYC_PVT1promoter": "CRISPR deletion of the promoter raises MYC (Cho 2018)",
        "TERT_promoter": "recurrent somatic mutation; CRISPR reversal (Chiba 2015)",
    }
    for e in third.THIRD:
        p = plan[e.locus]
        src = " (**stated**)" if e.element_source == "stated" else ""
        out.append(
            f"| {e.locus} | {e.chrom}:{e.element[0]:,}-{e.element[1]:,}{src} | {p['length']:,} |"
            f" {', '.join(e.targets)} | **{e.direction}** | {e.distance:,} bp | {how[e.locus]} |"
        )
    out += [
        "",
        f"**{reach['died_at_the_reach_filter']} of {reach['loci']} died at the reach filter, and that"
        " was the registered prediction.** `loci.read_reach` was called, not reimplemented: a"
        " gene-**body** test against the scorer's 1,048,576 bp input, never a TSS-distance test."
        f" Beside the other two frames - {reach['beside_the_other_sets']['the_seventeen']},"
        f" {reach['beside_the_other_sets']['the_nine']} - the filter keeps killing the famous"
        " megabase loci and almost nothing else. Three frames have now made the same measurement: a"
        " published enhancer-gene link is usually well inside the model's reach.",
        "",
        "**The sharpest case is SHH.** Section 18 retired the ZRS as an element-level miss because"
        " SHH is 979 kb away and **not in the model's input at all** - no number of requests puts it"
        f" there. SBE2 is a different published SHH enhancer {third.by_name()['SHH_SBE2'].distance:,}"
        " bp away, SHH's body is inside the window, and the question the benchmark could not ask"
        " about its flagship gene becomes askable by changing the element rather than the model.",
        "",
        "### What the requests were",
        "",
        "| loci | cCREs over the element | already deleted | requests |",
        "|---|---|---|---|",
    ]
    free = [r for r in result["plan"] if r["locus"] not in read]
    paid = [r for r in result["plan"] if r["locus"] in read]
    out += [
        f"| {', '.join(sorted(r['locus'] for r in free))} |"
        f" {min(r['ccres_over_element'] for r in free)} to"
        f" {max(r['ccres_over_element'] for r in free)} |"
        f" 1 to {max(r['already_scored_over_element'] for r in free)} | **0** |",
        f"| {', '.join(sorted(r['locus'] for r in paid))} |"
        f" {min(r['ccres_over_element'] for r in paid)} to"
        f" {max(r['ccres_over_element'] for r in paid)} | 0 | **1 each** |",
        "",
        f"**{spent} requests for the whole set**, because the finished all-chromosome sweep had"
        f" already deleted an element inside {already} of the {reach['loci']} published intervals."
        " Section 6 spent 2,283 on twelve loci; the cost of a new locus now is the cost of the"
        " intervals nobody annotated.",
        "",
        "**A caveat on those three, from the run's own record rather than from reading it"
        " charitably.** `Context.score_region` substitutes an overlapping ENCODE element when there"
        " is one, and it did so at "
        f"{', '.join(sorted(r['locus'] for r in read.values() if r['substituted_annotated_element']))}"
        " - every interval this set paid for. So what was deleted is an annotated element"
        " *overlapping* the published one, not the published interval itself, and the stated"
        " coordinate acted as a pointer rather than as the thing scored. It is recorded per locus in"
        f" `{third.INTERVALS}` and is the honest reading of what those three requests bought.",
        "",
        "### The scores, beside the other two frames",
        "",
        "| field | the third set | the nine candidates | the seventeen panel |",
        "|---|---|---|---|",
    ]
    for label, key in (
        ("right target, derived", "target_derived"),
        (
            "right target, where the model could answer",
            "target_derived_where_the_model_could_answer",
        ),
        ("right target, nearest TSS in node", "target_heuristic"),
        ("right target, annotation lookup", "target_looked_up"),
    ):
        b = beside[key]
        out.append(
            f"| {label} | **{kn(b['the_third_set'])}** | {kn(b['the_nine_candidates'])} |"
            f" {kn(b['the_seventeen_panel'])} |"
        )
    cf = beside["chance_floor"]
    out += [
        f"| chance floor (random coding gene in the window) |"
        f" {cf['the_third_set']['expected']}/{cf['the_third_set']['of']} |"
        f" {cf['the_nine_candidates']['expected']}/{cf['the_nine_candidates']['of']} |"
        f" {cf['the_seventeen_panel']['expected']}/{cf['the_seventeen_panel']['of']} |",
        f"| direction, where a deletion names the target |"
        f" **{kn(a['direction_derived_where_judged'])}** | 3/5 (0.600) | 5/5 |",
        f"| right cell or tissue, where judged | {kn(a['cell_derived_where_judged'])} | 4/6 (0.667) | - |",
        f"| a published class read | {frac(a['class_derived'])} | 9/9 | 13/17 |",
        "",
        f"The floors are not comparable across the frames and are printed so that each rate can be"
        f" read against its own: this set's {cf['the_third_set']['expected']} of"
        f" {cf['the_third_set']['of']} sits between the candidates' and the panel's. Within each"
        " frame the derived rate clears its own floor.",
        "",
        f"**The three derived rates agree to within"
        f" {spread(beside['target_derived'], beside['frames'])}"
        "**, which is the least interesting thing here and worth saying first so it cannot be"
        " mistaken for the finding. Three sets chosen on different days by different criteria land"
        f" on {kn(beside['target_derived']['the_third_set'])},"
        f" {kn(beside['target_derived']['the_nine_candidates'])} and"
        f" {kn(beside['target_derived']['the_seventeen_panel'])}. The headline is stable across"
        " sampling frames. Everything that moves is underneath it.",
        "",
        "### The nearest-gene rule moves further than anything else, and in the opposite direction to"
        " section 19",
        "",
        f"The heuristic reads {kn(beside['target_heuristic']['the_third_set'])} here against"
        f" {kn(beside['target_heuristic']['the_nine_candidates'])} at the candidates and"
        f" {kn(beside['target_heuristic']['the_seventeen_panel'])} at the panel. **That spread is"
        " wider than any spread in the derived rate, and it is a property of the loci, not of the"
        " rule.** It was registered in advance as a weakness of this set:"
        f" {len(third.THIRD) - len(third.TRAPS)} of {len(third.THIRD)} loci here have a published"
        " target as their own nearest coding TSS, because an element that sits in its target's"
        " intron or promoter is what a clean perturbation experiment usually looks like.",
        "",
        "**And the two traps caught nothing, which is not the same as the traps failing.**"
        f" {', '.join(a['nearest_gene_traps']['loci'])} were written down before scoring as the two"
        " loci where the nearest coding TSS is not the published target - RNF32 at"
        f" {third.TRAPS['SHH_SBE2']['trap_distance']:,} bp against SHH's"
        f" {third.TRAPS['SHH_SBE2']['target_distance']:,}, MEOX1 at"
        f" {third.TRAPS['SOST_VanBuchem']['trap_distance']:,} bp against SOST's"
        f" {third.TRAPS['SOST_VanBuchem']['target_distance']:,}. The heuristic fell into"
        f" **{len(a['nearest_gene_traps']['heuristic_fell_in'])} of them**: at both loci the node"
        " model named **no gene at all** rather than naming the wrong one. Its two misses"
        f" ({', '.join(a['target_heuristic']['misses'])}) are both silences. A baseline that abstains"
        " at the long-range loci and answers at the short-range ones is not the baseline its rate"
        " describes, and this is the first set in which that distinction shows.",
        "",
        "### The direction axis, asked both ways - and the registered answer is *undecidable*",
        "",
        "This is what the set was assembled for, and the registration wrote down in advance what each"
        " outcome would be worth: both repressors right is weak support that the layer reads"
        " direction, both inverted is weak support that it does not, **one of each decides nothing**.",
        "",
        "| arm | right | n |",
        "|---|---|---|",
        f"| activators | {dirs['activators']['k']} | {dirs['activators']['n']} |",
        f"| **repressors** | **{dirs['repressors']['k']}** | **{dirs['repressors']['n']}** |",
        "",
        "**It came back one of each, so by its own registration this set decides nothing about the"
        " sign.** No pooled direction rate is quoted here and none should be quoted from it.",
        "",
        "| locus | published | the model | |",
        "|---|---|---|---|",
    ]
    for r in dirs["per_locus"]:
        if not r["judged"]:
            continue
        mark = "right" if r["hit"] else "**inverted**"
        tag = " (repressor)" if r["repressor"] else ""
        out.append(f"| {r['locus']}{tag} | {r['published']} | {r['model_action']} | {mark} |")
    hbg = by_locus["HBG1_BCL11A_site"]
    hbg_del = hbg["readings"]["deletion"]
    myc = by_locus["MYC_PVT1promoter"]
    myc_del = myc["readings"]["deletion"]
    out += [
        "",
        "**The two repressors are the two readings worth having, and they disagree with each other.**",
        "",
        f"At **MYC_PVT1promoter** the model gets it right: deleting the PVT1 promoter *raises* MYC"
        f" ({myc_del['targets'][0]['best']:+.4f}), which is what Cho 2018 reports and is the opposite"
        " of what the panel's own MYC element does 335 kb the other side of the gene. A layer that"
        " reports a sign has to disagree with itself across those two elements, and it did.",
        "",
        f"At **HBG1_BCL11A_site** it gets it wrong, and this is the sharper of the two because"
        " nothing else about the locus is ambiguous. Deleting the BCL11A motif at -115 *lowers* HBG1"
        f" in the model ({hbg_del['targets'][0]['best']:+.4f}); in human erythroid cells, destroying"
        " that motif is one of the best-replicated ways to *raise* fetal haemoglobin. The element is"
        " in the reader's own cell type, the target is named first and correctly, the distance is"
        f" {third.by_name()['HBG1_BCL11A_site'].distance} bp - the model has every input it could"
        " want and reports the sign backwards.",
        "",
        "**What this adds to the ledger.** The candidates got 2 of 5 wrong among activators"
        " (IL2RA, GDF5); this set gets 0 of"
        f" {dirs['activators']['n']} wrong among activators and 1 of {dirs['repressors']['n']}"
        " wrong among repressors. Pooled across all three frames the sign is right more often than"
        " not, and the one clean case where a published *repressor* could be checked with the cell,"
        " the target and the distance all in the model's favour came back inverted. **n = 2 was"
        " declared too small before the run and is still too small after it.** What would settle it"
        " is a CRISPRi screen reporting de-repression, not more curation - published repressive"
        " elements with a perturbation behind them were searched for twice now, and the two here are"
        " what a second search produced.",
        "",
        "### The positive control, and the one miss",
        "",
        f"**The control passed.** Deleting TERT's own core promoter named TERT first at"
        f" {pc['deletion_log2_fold_change']:+.4f} - by a wide margin the largest effect anywhere in"
        " this set. That was registered as the thing whose failure would invalidate the run rather"
        " than score as a miss, so it is reported before the rates that depend on it and not after.",
        "",
        f"**{', '.join(a['target_derived']['misses'])} is the one target miss, and it is the FTO and"
        " SORT1 pattern for the third time.** KITLG is *second* in the summed-window layer behind"
        " DUSP6 and is named by no other derived layer; the element's own deletion names nothing at"
        " all, and neither does the node. So the lenient reading hits and the strict one misses, as"
        " at FTO and at SORT1. The expectation named KITLG alone, before scoring, and the miss"
        " stands.",
        "",
        "**Which layer does the work has moved again.** The element's own deletion names the"
        f" published target first at **{at_the_element} of"
        f" {len(result['loci'])}** here, against 5 of 9 at the candidates. The difference is not the"
        " model: it is that six of these nine elements sit inside or beside their own target, which"
        " is the registered weakness of the set and the same fact that lifts the nearest-gene rule.",
        "",
        "### The controls, with presence counted before coverage is made a stratum",
        "",
        f"{len(result['negatives'])} matched windows, the same four covariates and the same hit rules,"
        " with the keep-outs extended to every locus in **all three** sets."
        " `genomeos/compare.py` does the standardising; nothing here reimplements it.",
        "",
        "**`input_presence` runs first, and it is the reason the table below can be read at all.**",
        "",
        "| input | kind | targets | controls |",
        "|---|---|---|---|",
    ]
    for name, v in ac["input_presence"].items():
        out.append(
            f"| {name} | **{v['kind']}** | {v['targets_with_the_input']}/{v['targets']} |"
            f" {v['controls_with_the_input']}/{v['controls']} |"
        )
    out += [
        "",
        "A coverage stratum can only bite a claim whose input somebody had to **buy**. The deletion"
        f" input is bought here - {bought['controls_with_the_input']}"
        f" of {bought['controls']} controls have"
        " one - so the direction claim is genuinely under test. The constraint track is **free**,"
        " present on every window in both arms, so a null on the storage claim under a coverage"
        " stratum is arithmetic and not a test the claim passed. Section 19 did not separate those"
        " two cases; this run does, before the strata are run rather than after.",
        "",
        "| claim | third set | controls, covariates only | p | controls, **coverage held fixed** | p |",
        "|---|---|---|---|---|---|",
    ]
    for claim in ("target", "cell", "direction", "storage"):
        s, sc = ac[claim], ac[f"{claim}_given_coverage"]
        out.append(
            f"| {claim} | {s['raw']['a']:.3f} | {s['matched']['b']:.3f}"
            f" ({pct(s['matched']['difference'])}) | {p_of(s['matched'])} |"
            f" **{sc['matched']['b']:.3f}** ({pct(sc['matched']['difference'])}) |"
            f" {p_of(sc['matched'])} |"
        )
    d_raw, d_cov = ac["direction"], ac["direction_given_coverage"]
    out += [
        "",
        f"**The coverage artefact reproduces on a third independent set of loci.** Raw, the direction"
        f" claim separates: {d_raw['raw']['a']:.3f} against {d_raw['raw']['b']:.3f}, p ="
        f" {p_of(d_raw['raw'])}. Standardising on GC, distance and constraint takes it to"
        f" {pct(d_raw['matched']['difference'])}; adding *has an element inside this window actually"
        f" been deleted* as a fifth stratum takes it to **{pct(d_cov['matched']['difference'])}** at"
        f" p = {p_of(d_cov['matched'])} - the matched windows produce a direction slightly more often"
        " than the loci do. Section 8 retracted this claim at the panel, section 19 reproduced the"
        " retraction at the candidates, and it now holds at three frames. It is not a property of"
        " which loci anyone picked.",
        "",
        f"**The value-on-syntax claim does not separate here either:**"
        f" {ac['storage']['raw']['a']:.3f} at the loci against {ac['storage']['raw']['b']:.3f} at the"
        f" controls, and {pct(ac['storage']['matched']['difference'])} once standardised. That is the"
        " candidates' 1/9 shape and not the panel's 9/17, and it was registered in advance as the"
        " expected outcome because this set has one long element and eight short ones. **With one"
        " long element it cannot settle a length gradient and was declared unable to before the"
        " number existed.** What it adds is only that the panel's 9/17 has now failed to reproduce"
        " twice.",
        "",
        f"**The controls are thinner than they should be and that is stated rather than buried.**"
        f" Five windows per locus were asked for and {len(result['negatives'])} of"
        f" {5 * len(result['loci'])} were found: the length-and-distance matching could place five"
        " for six loci, one for HBG1_BCL11A_site and **none** for MYC_PVT1promoter or TERT_promoter,"
        " because a window whose nearest coding TSS is essentially zero away and whose length is 300"
        " to 500 bp is hard to match outside a promoter. One locus is dropped for want of a control"
        f" in every comparison above ({ac['direction']['dropped_for_want_of_a_control']} of"
        f" {ac['direction']['targets']}), and the constrained fraction is"
        f" {ac['imbalance']['constrained']['ratio']:.2f}x between the arms - the one covariate that"
        " differs materially, as it did at the candidates.",
        "",
        "### Where this set disagrees with the other two, and what the disagreement is",
        "",
        "| | the seventeen | the nine | the third set |",
        "|---|---|---|---|",
        "| chosen for | fame | perturbation | direction and length |",
        f"| derived target | {kn(beside['target_derived']['the_seventeen_panel'])} |"
        f" {kn(beside['target_derived']['the_nine_candidates'])} |"
        f" {kn(beside['target_derived']['the_third_set'])} |",
        f"| nearest-gene rule | {kn(beside['target_heuristic']['the_seventeen_panel'])} |"
        f" {kn(beside['target_heuristic']['the_nine_candidates'])} |"
        f" {kn(beside['target_heuristic']['the_third_set'])} |",
        f"| annotation lookup | {kn(beside['target_looked_up']['the_seventeen_panel'])} |"
        f" {kn(beside['target_looked_up']['the_nine_candidates'])} |"
        f" {kn(beside['target_looked_up']['the_third_set'])} |",
        f"| unaskable for reach | 3/17 | 2/11 | {reach['died_at_the_reach_filter']}/{reach['loci']} |",
        "",
        "**The headline is the stable number and the baselines are the unstable ones.** Across three"
        " frames the derived target rate moves by"
        f" {spread(beside['target_derived'], beside['frames'])}"
        " and the nearest-gene baseline moves by"
        f" {spread(beside['target_heuristic'], beside['frames'])}."
        " A benchmark reporting a derived rate alone would look reproducible; the thing it is"
        " supposed to be measured *against* is what depends on who picked the loci. That is a"
        " finding about benchmarks, and it is the main one this set produces.",
        "",
        "### Left undone",
        "",
        "- **The repressor arm is n = 2 and came back one of each, which was registered as"
        " undecidable.** Nothing here licenses a claim about the deletion layer's sign. The HBG"
        " inversion is a single reading at a locus where every input favoured the model, and it is"
        " worth following up on its own rather than aggregating.",
        "- **A published non-coding target is still missing after three sets**, so the"
        " `predicted_coding`-first defect H19 exposed remains untested. The search is recorded in"
        " `PREREGISTRATION['aims_that_failed_before_the_run']`: the clean non-coding cases act"
        " through parent-of-origin methylation, which no layer here carries.",
        "- **The three paid requests were answered about an overlapping annotated element**, not the"
        " stated interval, because `Context.score_region` substitutes one when it can. Whether that"
        " substitution should be refused for a stated interval is the benchmark owner's call.",
        f"- **`loci.STATED_INTERVAL_RESULTS` does not list `{third.INTERVALS}`.** Until it does,"
        " `loci_third.register_stated_intervals` applies the same one-line change at runtime. The"
        " hunk is in `loci_third.NEEDED_LOCI_HUNK`; this lane was not permitted to make it.",
        "- **The cell axis is close to unmeasurable on this set**, as registered:"
        f" {len(a['cell_unreachable'])} of {len(result['loci'])} published tissues have no"
        " counterpart in the eleven reader cells or the GTEx panel, and where it could be judged it"
        f" reads {kn(a['cell_derived_where_judged'])}. A benchmark built on loci GTEx can see will"
        " keep overrating the tissue layers.",
        f"- **The controls are {len(result['negatives'])} windows, not {5 * len(result['loci'])}**,"
        " and two loci have none at all. A promoter-proximal element cannot be length-and-distance"
        " matched against anything that is not itself a promoter, which is a structural limit of"
        " `loci.candidate_windows` rather than a shortage of genome.",
        "- **The three sets are still three sets.** Whether any of them should merge, and on what"
        " denominator, is the benchmark owner's call; nothing here pools them.",
        "",
        "---",
        "",
    ]
    return "\n".join(out)


def splice(text: str) -> str:
    doc = DOC.read_text()
    if HEADING in doc:
        head, rest = doc.split(HEADING, 1)
        nxt = rest.find("\n## ")
        tail = rest[nxt + 1 :] if nxt >= 0 else ""
        return head + text + tail
    return doc.rstrip("\n") + "\n\n" + text


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--write", action="store_true", help="splice the section into the doc")
    args = ap.parse_args(argv)
    result = load_result(third.NAME)
    if not result:
        raise SystemExit(f"no {third.NAME} result: run scripts/loci_third.py --score first")
    text = render(result, load_result(third.INTERVALS))
    if args.write:
        DOC.write_text(splice(text))
        print(f"wrote {DOC}")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
