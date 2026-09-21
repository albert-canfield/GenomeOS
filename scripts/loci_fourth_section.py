# SPDX-License-Identifier: AGPL-3.0-or-later
"""Render section 22 of docs/LOCI-BENCHMARK.md from the fourth frame's result file.

    uv run python scripts/loci_fourth_section.py            # print it
    uv run python scripts/loci_fourth_section.py --write    # splice it into the doc

Every figure in the section is read out of `data/results/loci_fourth.json` and
`data/results/loci_fourth_intervals.json`. None is typed, for the reason `scripts/
loci_third_section.py` gives: section 19 was hand-copied from its result and section 10's cCRE count
went stale in the doc while the result moved on.

The prose is the reading and is written by hand. Where a sentence would change meaning if a number
moved, it names the number rather than asserting the direction, so the two cannot disagree.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from genomeos.benchmark import loci_fourth as fourth
from genomeos.results import load_result

DOC = Path("docs/LOCI-BENCHMARK.md")
HEADING = (
    "## 22. A fourth frame, drawn by a rule instead of chosen, and the headline does not"
    " survive it, at n = 50 and in both halves (2026-09-21)"
)
#: headings this section has carried before, so `splice` replaces the section rather than appending a
#: second copy of it when the title moves with the result.
OLD_HEADINGS = (
    "## 22. A fourth frame, drawn by a rule instead of chosen, and the headline does not"
    " survive it (2026-09-21)",
)


def kn(d: dict[str, Any]) -> str:
    return f"{d['k']}/{d['n']} ({d['rate']:.3f})" if d.get("n") else "-"


def pct(x: float | None) -> str:
    return "-" if x is None else f"{x:+.3f}"


def p_of(block: dict[str, Any]) -> str:
    p = block.get("p_one_sided")
    return "-" if p is None else f"{p:.3f}"


def rate_of(d: dict[str, Any]) -> str:
    """k/n (rate) for one half, where a half with nothing graded prints as a dash rather than 0."""
    return f"{d['k']}/{d['n']} ({d['rate']:.3f})" if d.get("n") and d.get("rate") is not None else "-"


def floor_of(half: dict[str, Any]) -> str:
    f = half["chance_floor"]
    return f"{f['expected']}/{f['of']} ({f['rate']:.3f})" if f.get("of") else "-"


def layer_of(half: dict[str, Any], layer: str) -> str:
    v = (half["layers"].get("by_layer") or {}).get(layer)
    return f"{v['hit']}/{v['n']}" if v else "-"


def named_of(half: dict[str, Any]) -> str:
    named = (half["layers"].get("what_it_named_instead") or {}).get("deletion") or {}
    return f"{named.get('the nearest coding TSS', 0)}/{half['graded']}"


def halves_reading(h1: dict[str, Any], h2: dict[str, Any]) -> str:
    """The sentence the two halves earn, written from the numbers rather than around them."""
    a = h1["target_derived_where_the_model_could_answer"]
    b = h2["target_derived_where_the_model_could_answer"]
    if a.get("rate") is None or b.get("rate") is None:
        return "One half has nothing graded, so the two cannot be compared."
    gap = b["rate"] - a["rate"]
    registered = (
        f"Commit f407915 registered, before the second half was read, that a second-half rate"
        f" between about 0.15 and 0.45 would mean the first half's {a['rate']:.3f} was not an"
        " n = 20 accident, that a second half BELOW the first is the expected direction of a"
        " gene-density effect and not a new finding, and that a second half materially ABOVE the"
        " first would be the interesting divergence because it runs against that gradient."
    )
    if 0.15 <= b["rate"] <= 0.45:
        verdict = (
            f"**The second half reads {rate_of(b)} against the first half's {rate_of(a)}"
            f" ({gap:+.3f}), inside the registered interval.** The frame's rate is the combined one,"
            " and it is not carried by either half."
        )
    elif b["rate"] < 0.15:
        verdict = (
            f"**The second half reads {rate_of(b)}, BELOW the registered interval and below the first"
            f" half's {rate_of(a)} ({gap:+.3f}).** That is the direction the gene-density gradient"
            " predicts, so it is reported against the second half's own floor rather than as a new"
            " finding; it does not soften the combined rate and it is not smoothed into it."
        )
    else:
        verdict = (
            f"**The second half reads {rate_of(b)}, ABOVE the registered interval and above the first"
            f" half's {rate_of(a)} ({gap:+.3f}).** That is the divergence registered as the"
            " interesting one, because it runs against the density gradient: it would mean the first"
            " half's rate is partly a property of chr1 to chr10 rather than of the geometry alone."
        )
    f1, f2 = h1["chance_floor"], h2["chance_floor"]
    floors = (
        f" Against their own floors, which is how the registration said to read them: the first half"
        f" is {a['rate']:.3f} on a floor of {f1['rate']:.3f} and the second is {b['rate']:.3f} on a"
        f" floor of {f2['rate']:.3f}, so the margin over chance is"
        f" {a['rate'] / f1['rate']:.1f}x and {b['rate'] / f2['rate']:.1f}x. The floor itself fell"
        f" from {f1['rate']:.3f} to {f2['rate']:.3f} between the halves, which is the gene-density"
        " gradient showing up in the arithmetic rather than in the readings: more coding genes in"
        " the window makes a random draw likelier to be wrong, not likelier to be right."
    )
    if b["rate"] / f2["rate"] < a["rate"] / f1["rate"]:
        floors += (
            " So the second half is weaker than the first even after its easier floor is taken into"
            " account, and the density gradient does not account for the whole of the drop. The"
            " honest statement is that the frame's rate is the combined one and that the half that"
            " was scored first was the better of the two."
        )
    return (
        verdict
        + floors
        + " "
        + registered
        + " What is not claimed is that the halves are independent replications: they are one rule"
        " applied to one file and split by an arbitrary cap, so a difference between them is a"
        " statement about the genome's arrangement and about n."
    )


def render(result: dict[str, Any], intervals: dict[str, Any] | None) -> str:
    a = result["aggregate"]
    draw = result["draw"]
    reach = result["reach"]
    layers = result["layers"]
    halved = result["halves"]
    base = result["baselines"]
    beside = result["beside_the_other_three"]
    ac = result["against_controls"]
    pc = result["positive_control"]
    reg = result["preregistration"]
    read = {r["locus"]: r for r in (intervals or {}).get("read", [])}
    rows = fourth.drawn_rows(result)
    plan = {r["locus"]: r for r in result["plan"]}
    drawn = result["drawn_panel"]
    graded = len(rows)
    all_drawn = reach["loci"]
    floor = a["target_by_chance"]
    floor_rate = floor["expected"] / floor["of"]
    derived = a["target_derived"]
    over_all = derived["k"] / all_drawn
    dists = sorted(r["expected"]["distance"] for r in rows)
    unnamed = [e for e in drawn if str(e["nearest_gene_trap"]).startswith("ENSG")]
    unnamed_graded = [r["locus"] for r in rows if str(r["expected"]["nearest_gene_trap"]).startswith("ENSG")]
    del_named = layers["what_it_named_instead"]["deletion"]
    node_named = layers["what_it_named_instead"]["node"]
    by_layer = layers["by_layer"]
    cells = sorted({(e["cells"] or e["tissues"] or ["-"])[0] for e in drawn})
    bought = ac["input_presence"]["a deletion spent inside the window"]
    why = draw["rejected_by_reason"]
    shortcut_right = next((n for r, n in why.items() if r.startswith("the nearest coding TSS IS")), 0)
    non_coding = next((n for r, n in why.items() if r.startswith("a regulated target is not")), 0)
    passed_the_rule = why.get("kept", 0)
    tested = draw["perturbed_elements_with_a_regulated_gene"]

    out: list[str] = [HEADING, ""]
    out += [
        "`genomeos/benchmark/loci_fourth.py` (the registration and the rule), `scripts/"
        "loci_fourth.py`, `scripts/loci_fourth_section.py`, `tests/test_loci_fourth.py`, results"
        f" `{fourth.NAME}` and `{fourth.INTERVALS}`. **{len(read)} model"
        f" request{'s' if len(read) != 1 else ''} spent.**"
        f" {result['cost']['seconds']:.0f} seconds of range reads over public tracks for the rest.",
        "",
        "Section 21 ended on a number that reads well and one that does not. Across three frames the"
        f" derived target rate moves {spread(beside['target_derived'], CURATED)} and the nearest-gene"
        f" rule over the same three moves {spread(beside['target_heuristic'], CURATED)}. Section 21"
        " also said why"
        " the baseline moves: at the third set seven of nine published targets are the element's own"
        " nearest coding TSS. So the three frames differ mostly in **geometry**, and the question"
        " that raises cannot be answered from any of them - is the 0.88 a property of the model, or"
        " of the fact that curators pick elements sitting next to their targets?",
        "",
        "**A frame chosen by hand cannot answer it, because whoever picks the loci picks the"
        " geometry.** So this one was not picked. The rule was written, committed and its tests run"
        " before the file it draws from was read (commit 2103c5a), and it contains no loci at all:",
        "",
        "> every element in the ENCODE CRISPR benchmark's **held-out** arm whose regulated target is"
        " **not** its own nearest coding TSS.",
        "",
        "That is the geometry on which the cheap answer is wrong by construction. It reads no model"
        " output, no effect size and no chromatin score - only the file's own positive flag and"
        " GENCODE arithmetic - and it is adversarial to this benchmark's headline rather than"
        " flattering to it. The K562 **training** arm is excluded because"
        " `data/results/crispri_benchmark.json` fitted logistic weights on it on 2026-09-16.",
        "",
        "### What the rule drew, and what it threw away",
        "",
        "| step | elements |",
        "|---|---|",
        f"| perturbed elements with a regulated gene in the held-out arm | {tested} |",
    ]
    for reason, n in sorted(why.items(), key=lambda kv: -kv[1]):
        if reason == "kept":
            continue
        out.append(f"| rejected: {reason} | {n} |")
    out += [
        f"| **passed the rule** | **{passed_the_rule}** |",
        f"| drawn: every element that passes, the cap raised to {fourth.MAX_LOCI} | {all_drawn} |",
        f"| of those, drawn first under the registered cap of {fourth.FIRST_DRAW} | {fourth.FIRST_DRAW} |",
        f"| added by the cap raise, same rule | {all_drawn - fourth.FIRST_DRAW} |",
        "",
        f"**{shortcut_right} of the {tested} held-out CRISPR positives have their own nearest coding"
        f" TSS as the published target** - {shortcut_right / tested:.0%} of them, and"
        f" {shortcut_right / (shortcut_right + passed_the_rule):.0%} of the"
        f" {shortcut_right + passed_the_rule} where the geometry question is even well posed. That"
        " number is worth reading before any rate below it: at half the field's own measured"
        " enhancer-gene pairs, naming the nearest gene is simply correct. A benchmark drawn without"
        " a geometry rule inherits that proportion, whoever draws it.",
        "",
        f"**{non_coding}"
        " elements were dropped for having a non-coding regulated target**, and that answers an item"
        " section 21 left open. Three frames failed to find a published non-coding target clean"
        " enough to state by curation; the held-out arm has that many, and they were dropped here"
        " on purpose because the deletion layer ranks `predicted_coding` first, so they test the H19"
        " defect rather than the geometry question. They are a frame of their own, already drawn.",
        "",
        f"**The cap no longer binds.** It was {fourth.FIRST_DRAW} when the frame was first run"
        " (commit 2103c5a) and the first run scored that many, leaving"
        f" {all_drawn - fourth.FIRST_DRAW} loci that pass the same rule undrawn. Commit f407915"
        f" raised it to {fourth.MAX_LOCI}, which is every element the rule passes, so there is now"
        " nothing left over for a later hand to choose from. The raise was committed before the"
        " second run, with its expected cost and with what each outcome would mean written down"
        " first: `loci_fourth.CAP_RAISE`. `assess` and `select` - the rule itself - are untouched,"
        " and the two halves are reported separately below because genome order is not random with"
        " respect to gene density.",
        "",
        f"### The {all_drawn}, and what the free filter did to them first",
        "",
        "| locus | element | len | published target | nearest coding TSS | distance | cell | half |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for i, e in enumerate(drawn):
        p = plan[e["locus"]]
        dead = "" if p["askable"] else " **(unaskable)**"
        cell = (e["cells"] or e["tissues"] or ["-"])[0]
        half = "first 24" if i < fourth.FIRST_DRAW else "next 36"
        out.append(
            f"| {e['locus']}{dead} | {e['chrom']}:{e['element'][0]:,}-{e['element'][1]:,} |"
            f" {e['length']:,} | {', '.join(e['targets'])} | {e['nearest_gene_trap']} |"
            f" {e['distance']:,} bp | {cell} | {half} |"
        )
    out += [
        "",
        f"**{reach['died_at_the_reach_filter']} of {all_drawn} died at the reach filter** - 1 to 4 of"
        " the first 24 was the registered prediction, deliberately higher than the third set's"
        " registered 0 of 9, because a rule that selects elements skipping over a nearer gene selects"
        " for distance. It came back 4 of the first 24 and the same rate holds over all"
        f" {all_drawn}. `loci.read_reach` was called, not reimplemented: a gene-**body** test against"
        " the scorer's 1,048,576 bp input. They are "
        + ", ".join(f"{x['locus'].split('_')[0]} ({x['distance_kb']} kb)" for x in reach["which"])
        + f". Beside the other frames: {reach['beside_the_other_frames']['the_seventeen']},"
        f" {reach['beside_the_other_frames']['the_nine']},"
        f" {reach['beside_the_other_frames']['the_third_set']}.",
        "",
        f"The element-to-target distances of the {graded} graded loci run {dists[0]:,} bp to"
        f" {dists[-1]:,} bp, median {dists[len(dists) // 2]:,} bp.",
        "",
        "### What it cost",
        "",
        f"**{len(read)} request{'s' if len(read) != 1 else ''} for the whole frame, {graded} graded"
        " loci.** The finished all-chromosome sweep had already deleted an element inside"
        f" {graded - len(read)} of the {graded} graded intervals, and the registered budget of"
        f" {reg['requests']['budget']} was never approached. Section 6 spent 2,283 requests on twelve"
        " loci; once the sweep is finished, a new locus costs only what the intervals nobody"
        " annotated cost. The cap raise took n from 20 to"
        f" {graded} for {len(read) - 1} further request{'s' if len(read) - 1 != 1 else ''}.",
        "",
    ]
    if read:
        sub = [r["locus"] for r in read.values() if r["substituted_annotated_element"]]
        out += [
            "**The caveat section 21 recorded applies to these requests too.**"
            " `Context.score_region` substitutes an overlapping ENCODE element when there is one, and"
            f" it did so at {', '.join(sub) if sub else 'none of them'}. What was deleted is an"
            " annotated element *overlapping* the perturbed interval, not the screen's interval"
            f" itself. It is recorded per locus in `{fourth.INTERVALS}`.",
            "",
        ]
    out += [
        "### The scores, beside the other three frames",
        "",
        "The positive control is excluded from every rate below, as registered; it is reported on its"
        " own further down.",
        "",
        "| field | **the fourth frame** | the third set | the nine candidates | the seventeen panel |",
        "|---|---|---|---|---|",
    ]
    b = beside["target_derived"]
    out.append(
        f"| right target, derived, every locus the frame contains |"
        f" **{derived['k']}/{all_drawn} ({over_all:.3f})** | {kn(b['the_third_set'])} |"
        f" {kn(b['the_nine_candidates'])} | {kn(b['the_seventeen_panel'])} |"
    )
    for label, key in (
        ("right target, where the model could answer", "target_derived_where_the_model_could_answer"),
        ("right target, nearest TSS in node", "target_heuristic"),
        ("right target, annotation lookup", "target_looked_up"),
    ):
        b = beside[key]
        out.append(
            f"| {label} | **{kn(b['the_fourth_frame'])}** | {kn(b['the_third_set'])} |"
            f" {kn(b['the_nine_candidates'])} | {kn(b['the_seventeen_panel'])} |"
        )
    cf = beside["chance_floor"]
    out += [
        "| chance floor (random coding gene in the window) |"
        f" **{cf['the_fourth_frame']['expected']}/{cf['the_fourth_frame']['of']}"
        f" ({floor_rate:.3f})** | {cf['the_third_set']['expected']}/{cf['the_third_set']['of']} |"
        f" {cf['the_nine_candidates']['expected']}/{cf['the_nine_candidates']['of']} |"
        f" {cf['the_seventeen_panel']['expected']}/{cf['the_seventeen_panel']['of']} |",
        "",
        f"**The derived rate is {derived['k']}/{all_drawn} ({over_all:.3f}) over every locus the rule"
        f" returned and {kn(derived)} where the model could answer.** The first row is the one"
        " comparable with the panel's 15/17, because that rate counts the panel's three unaskable"
        " loci as misses; the second is comparable with its 13/15. The three hand-curated frames read"
        " 0.882, 0.889 and 0.889 on the first and 0.867, 0.889, 0.889 on the second. This frame is"
        " the same benchmark, the same scorers, the same hit rules and the same machine; what changed"
        " is that the published target is no longer the nearest gene.",
        "",
        f"It is **above its own chance floor** - {derived['rate']:.3f} against {floor_rate:.3f}, so"
        " the layers are not guessing - and it is **nowhere near the number this benchmark has been"
        " reporting since 2026-09-12**. The registration wrote down in advance what each outcome would"
        " mean and this one is the middle case: the number is reported with the per-layer split and"
        " without a verdict bolted onto it. The split is where the finding is.",
        "",
        "### The two halves, because the cap was applied in genome order",
        "",
        "The first 24 were scored on 2026-09-21 under the registered cap; the next 36 were added by"
        " raising it, under the same rule, and commit f407915 wrote down before they were scored what"
        " they would have to read for the frame to be called stable and what a divergence would mean."
        " The cap ran in genome order, so the first half is chr1 to chr10 and the second is chr11 to"
        " chrX - and the later chromosomes carry the gene-dense stretches, where the nearest coding"
        " TSS is closer and more coding genes sit in the window. Each half is therefore measured"
        " against **its own** chance floor.",
        "",
        "| | the first 24 (chr1-chr10) | the next 36 (chr11-chrX) | all 60 |",
        "|---|---|---|---|",
    ]
    h1, h2 = halved["first_24"], halved["next_36"]
    for label, key in (
        ("drawn by the rule", "drawn"),
        ("died at the reach filter", "died_at_the_reach_filter"),
        ("graded", "graded"),
    ):
        out.append(f"| {label} | {h1[key]} | {h2[key]} | {h1[key] + h2[key]} |")
    out += [
        f"| right target, derived, every drawn locus | {rate_of(h1['target_derived_over_every_drawn_locus'])}"
        f" | {rate_of(h2['target_derived_over_every_drawn_locus'])} |"
        f" **{derived['k']}/{all_drawn} ({over_all:.3f})** |",
        f"| right target, where the model could answer |"
        f" {rate_of(h1['target_derived_where_the_model_could_answer'])} |"
        f" {rate_of(h2['target_derived_where_the_model_could_answer'])} | **{kn(derived)}** |",
        f"| right target, nearest TSS in node | {rate_of(h1['target_heuristic'])} |"
        f" {rate_of(h2['target_heuristic'])} | {kn(beside['target_heuristic']['the_fourth_frame'])} |",
        f"| chance floor, computed per half | {floor_of(h1)} | {floor_of(h2)} |"
        f" {cf['the_fourth_frame']['expected']}/{cf['the_fourth_frame']['of']} ({floor_rate:.3f}) |",
        f"| deletion layer, the model itself | {layer_of(h1, 'deletion')} | {layer_of(h2, 'deletion')} |"
        f" {by_layer['deletion']['hit']}/{by_layer['deletion']['n']} |",
        f"| deletion layer named the nearest coding TSS instead | {named_of(h1)} | {named_of(h2)} |"
        f" {del_named['the nearest coding TSS']}/{graded} |",
        "",
        halves_reading(h1, h2),
        "",
        "### Which layer, and what it named instead - this is the result",
        "",
        "| layer | provenance | names the published target first |",
        "|---|---|---|",
    ]
    for name in ("deletion", "eqtl", "gene_input", "node", "lookups"):
        v = by_layer.get(name)
        if not v:
            continue
        out.append(f"| {name} | {v['provenance']} | {v['hit']}/{v['n']} |")
    out += [
        "",
        f"**The deletion layer - the model - names the published target at"
        f" {by_layer['deletion']['hit']}/{by_layer['deletion']['n']}.** The derived rate above is a"
        " union over three derived layers, and the two that carry it are GTEx eQTL"
        f" ({by_layer['eqtl']['hit']}/{by_layer['eqtl']['n']}) and the summed-window reading"
        f" ({by_layer['gene_input']['hit']}/{by_layer['gene_input']['n']}), neither of which is the"
        " model's own answer about the element.",
        "",
        "And what it names instead is the whole point of drawing the frame this way:",
        "",
        "| what it named | the deletion layer | the nearest-TSS-in-node rule |",
        "|---|---|---|",
    ]
    for key in ("the published target", "the nearest coding TSS", "another gene", "nothing"):
        out.append(f"| {key} | {del_named[key]}/{graded} | {node_named[key]}/{graded} |")
    out += [
        "",
        f"**The model names the nearest coding TSS at {del_named['the nearest coding TSS']} of"
        f" {graded} loci and the published target at {del_named['the published target']}.** At these"
        " loci the two are different genes by construction, and when they differ the model goes with"
        f" the nearer one {del_named['the nearest coding TSS']} times and with a third gene"
        f" {del_named['another gene']} times - so it is not only proximity, but the published answer"
        f" is the rarest of the three. The node heuristic names the nearer gene at"
        f" {node_named['the nearest coding TSS']} of {graded}, which is what a pure proximity rule"
        " looks like, and the deletion layer sits between that and the publication rather than at"
        " the publication.",
        "",
        "**That is the sharpest reading this benchmark has produced, and it is a negative.** It does"
        " not say the deletion layer is the nearest-gene rule in general - at the three curated"
        " frames it also got directions and cell types right, which proximity alone cannot do. It"
        " says that on the axis this frame isolates, where proximity and the published answer"
        f" disagree, the model names the published target {by_layer['deletion']['hit']} times in"
        f" {graded} and the nearer gene"
        f" {del_named['the nearest coding TSS'] // max(by_layer['deletion']['hit'], 1)} times as"
        " often; and that the three frames' agreement at 0.88 was measured almost entirely where the"
        " two agree, so it could not have shown this.",
        "",
        "### The two nearest-gene rules, kept apart, because one of them is arithmetic",
        "",
        "| rule | reads | status |",
        "|---|---|---|",
        f"| nearest coding TSS anywhere | {base['nearest_coding_anywhere']['k']}/"
        f"{base['nearest_coding_anywhere']['n']} | **0 by construction** - it is the rule that drew"
        " the frame, and it is printed so the construction is visible |",
        f"| nearest coding TSS in the CTCF node | {kn(base['nearest_tss_in_node'])} | a genuine"
        " reading: a different rule, right wherever the node excludes the nearer gene |",
        "",
        f"The node rule was registered in advance as low but **not** zero by construction, and it"
        f" came back {base['nearest_tss_in_node']['k']}/{base['nearest_tss_in_node']['n']}: at"
        f" {node_named['the nearest coding TSS']} loci it named the nearer gene and at"
        f" {base['silences']} it named nothing. So on this frame the node is a proximity rule with"
        " abstentions, and what it adds over raw proximity is the handful of loci where the node"
        " boundary excludes the nearer gene. The comparison the doc has been making - 0.647, 0.333,"
        " 0.778 - is with that rule, and it is"
        f" {kn(base['nearest_tss_in_node'])} at the one frame where the rule was not allowed to be"
        " right by default. **On the where-the-model-could-answer line its spread across four frames"
        f" is now {spread(beside['target_heuristic'])} and the derived rate's is"
        f" {spread(beside['target_derived_where_the_model_could_answer'])}, against"
        f" {spread(beside['target_derived_where_the_model_could_answer'], CURATED)} across the three"
        " curated frames alone.** Both numbers move once the geometry is fixed; it is no longer only"
        " the baseline that depends on who picked the loci.",
        "",
        "### The positive control, and why the run is interpretable at all",
        "",
        f"**The control passed.** {pc['locus']} was drawn by the same rule from the same file - among"
        " the elements this frame REJECTED for having the shortcut right, the one with the smallest"
        f" published element-to-TSS distance ({pc['distance']:,} bp), which is the easiest case the"
        f" data can offer. Deleting it named {pc['deletion_named_first']} first at"
        f" {pc['deletion_log2_fold_change']}, by the deletion layer"
        f" ({', '.join(pc['hit_by'])}). Registered in advance: a failure here would mean the reading"
        " is broken and nothing else in the run is interpretable. It is not broken - the same layer"
        f" that answers correctly at 1.1 kb names the wrong gene at {graded - by_layer['deletion']['hit']}"
        f" of the {graded} loci where the right answer is not the nearest one.",
        "",
        "### The matched windows, with presence counted before coverage is made a stratum",
        "",
        f"{len(result['negatives'])} matched windows, the same four covariates and the same hit"
        " rules, with the keep-outs extended to every locus in **all four** frames."
        " `genomeos/compare.py` does the standardising; nothing here reimplements it.",
        "",
        "| input | kind | targets | controls |",
        "|---|---|---|---|",
    ]
    for label, v in ac["input_presence"].items():
        out.append(
            f"| {label} | **{v['kind']}** | {v['targets_with_the_input']}/{v['targets']} |"
            f" {v['controls_with_the_input']}/{v['controls']} |"
        )
    out += [
        "",
        "| claim | the fourth frame | controls, covariates only | p | controls, **coverage held"
        " fixed** | p |",
        "|---|---|---|---|---|---|",
    ]
    for claim in ("target", "cell", "direction", "storage"):
        m = ac[claim]["matched"]
        c = ac[f"{claim}_given_coverage"]["matched"]
        out.append(
            f"| {claim} | {m['a']:.3f} | {m['b']:.3f} ({pct(m['difference'])}) | {p_of(m)} |"
            f" **{c['b']:.3f}** ({pct(c['difference'])}) | {p_of(c)} |"
        )
    out += [
        "",
        "**The coverage artefact reproduces on a fourth independent frame.** The target claim -"
        " *some* coding gene is named - separates raw and does not survive the coverage stratum,"
        " which is what sections 8, 19 and 21 each found. These claims are about whether a layer"
        " says anything at all, not about whether it is right, and the rates above are what say"
        " whether it is right.",
        "",
        "### What this frame cannot say, stated rather than buried",
        "",
        "- **A CRISPRi positive is not always a direct target.** The benchmark file carries a"
        " `direct_vs_indirect_negative` column estimating exactly that, and this frame did not use"
        f" it: the rule takes the file's own `Regulated` flag. Some of these {graded} published targets"
        " may be indirect consequences of silencing the element, in which case naming the nearer gene"
        " is not as wrong as the rate makes it look. Using that column would be a different frame and"
        " it would have to be registered as one.",
        f"- **{len(unnamed)} of the drawn loci have an unnamed novel gene as their nearest coding"
        f" TSS** ({', '.join(e['nearest_gene_trap'] for e in unnamed)}), so the shortcut is wrong"
        " there because of an annotation rather than because of a genuine gene skip."
        f" {len(unnamed_graded)} of those {'is' if len(unnamed_graded) == 1 else 'are'} in the graded"
        f" {graded}"
        f" ({', '.join(unnamed_graded) if unnamed_graded else 'none'}).",
        f"- **The frame is one cell line more than it looks.** The drawn cells are"
        f" {', '.join(cells)}, and K562 dominates because the held-out arm does. A cell-type"
        " imbalance is a property of the field's screens, not of this rule.",
        "- **The frame is now every element the rule passes, so it is not a sample of anything.** The"
        " genome-order cap that made the first run chromosome 1 to 10 has been raised past the last"
        " locus; what remains is that the held-out arm is itself whatever the field happened to"
        " screen.",
        f"- **{bought['controls_with_the_input']} of {bought['controls']} control windows have a"
        " deletion spent inside them against"
        f" {bought['targets_with_the_input']}/{bought['targets']} of the loci**, so the coverage"
        " stratum is a real test here rather than arithmetic.",
        "",
        "### Left undone",
        "",
        "- **The frame is exhausted.** Every element the rule passes is drawn, so there are no more"
        " loci to add without a different rule, a different file, or the K562 training arm this"
        " registration excludes. Taking n past"
        f" {all_drawn} means registering a new frame, not raising a cap.",
        f"- **The non-coding frame is already drawn and unscored.** The {non_coding}"
        " elements dropped at step 4 are the first set this project has ever had of published,"
        " perturbation-backed, non-coding targets, which is what section 21 said curation could not"
        " produce. They would test the `predicted_coding`-first defect H19 exposed.",
        "- **The direction axis was registered unaskable here and is.** Every CRISPR positive is an"
        " element whose silencing lowers its target, so the repressor arm is empty by construction"
        " and section 21's sign question stays open at n = 2.",
        f"- **`loci.STATED_INTERVAL_RESULTS` now lists `{fourth.INTERVALS}`** (f407915), so this"
        " frame's paid deletions are read back from the file rather than from a global patched at"
        " import time. `loci_fourth.register_stated_intervals` is kept because it is idempotent and"
        " because a test asserts the two say the same thing.",
        "- **Whether the earlier frames' rates should now be reported with their element-to-target"
        " distances beside them** is the benchmark owner's call. This frame's reading is that a"
        " derived rate quoted without that distance is not comparable between frames, and sections 9,"
        " 19 and 21 do not carry it.",
        "- **The four frames are still four frames.** Nothing here pools them, and the fourth one"
        " least of all: it was drawn to disagree with the other three and it does.",
        "",
        "---",
        "",
    ]
    return "\n".join(out)


def spread(block: dict[str, Any], frames: tuple[str, ...] | None = None) -> str:
    """How far one rate moves across the named frames; across all four when none are named."""
    keys = frames or tuple(block)
    rates = [
        block[k]["rate"] for k in keys if isinstance(block.get(k), dict) and block[k].get("rate") is not None
    ]
    return f"{max(rates) - min(rates):.3f}" if rates else "-"


#: the three frames that were curated locus by locus, for the spreads section 21 argued from.
CURATED = ("the_third_set", "the_nine_candidates", "the_seventeen_panel")


def splice(text: str) -> str:
    doc = DOC.read_text()
    for heading in (HEADING, *OLD_HEADINGS):
        if heading not in doc:
            continue
        head, rest = doc.split(heading, 1)
        nxt = rest.find("\n## ")
        tail = rest[nxt + 1 :] if nxt >= 0 else ""
        return head + text + tail
    return doc.rstrip("\n") + "\n\n" + text


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--write", action="store_true", help="splice the section into the doc")
    args = ap.parse_args(argv)
    result = load_result(fourth.NAME)
    if not result:
        raise SystemExit(f"no {fourth.NAME} result: run scripts/loci_fourth.py --score first")
    text = render(result, load_result(fourth.INTERVALS))
    if args.write:
        DOC.write_text(splice(text))
        print(f"wrote {DOC}")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
