# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""What `compile.py` emits per element now, and what "every coding target at its bar" would emit.

ROADMAP §5 item 7g's next step: "the compiler emits every coding target at its bar (a `compile.py`
change, registered against the compiled programs' census tests)". This module is the census half, and
it is written and committed BEFORE anything in `compile.py` moves, because the counts it reports are
what any later change has to be read against.

**What compile.py emits now.** `compile.compile_chromosome` walks `_attributed(chrom)` -- the elements
`targets.attributed` gives it, which are the elements of the whole-chromosome sweep that carry a
`predicted_coding` record -- and for each one writes a single `targets:` property and a single `rule`
line, both naming `predicted_coding["gene"]`. One element, one target, one rule. The other genes the
sweep moved at that element are in the response cache (`targets.ElementResponses`, every gene of the
scorer's window) and are not written: the single target is a property of the compact projection, not
of the sweep (`targetchoice`'s opening paragraph, 79aa1e2).

**What "every coding target at its bar" would emit.** The bar is the compiler's own,
`predict.enhancer_target.MIN_EFFECT` = 0.1 |log2 fold change|, the same object, imported here and not
copied; the size rule is `targets.genes_at_bar`'s, which is `predict_target`'s own (the larger of the
predicted fall and the predicted rise over every track, ties to the fall, strongest first). The
population of elements does not change: only the elements that already carry a compiled rule are
reached, so **no element without a compiled target gains one** -- a `region` block (the budget's
unknown blocks, `U_<chrom>_<start>`) is not a scored element and has no cache record, and an
`element <id>_measured` block of the experimental layer carries no predicted rule either way.

The candidate list is `targetchoice.candidate_list`, the list registered at b746a23 and read on the
held-out CRISPRi pairs at 715d830: the compiled target first, then every protein-coding gene of the
window at the bar, strongest first, deduplicated. An element the cache does not hold keeps its
compiled target alone.

**What the series says this is worth.** 7/7d/7h read composed forms against measured RNA and every one
was worse than no rule. M2h (ac5580b) then found the compiled rules carry real cell-pattern
information at about a fifth of their compiled magnitude (held-out Delta MSE +0.04339 [98.33% CI
+0.0379, +0.04894], fitted scale +0.2146, about 2% of the variance), and M2i (f6e71c6) replicated that
on RNA fresh to this project with nothing fitted. Its ceiling: independence from AlphaGenome's training
domain is unknown for 13 of 13 reader cells (18de375). An emitted target is a model's attribution --
AlphaGenome's predicted change on deleting the element -- and never a measured regulatory
relationship, whichever count this census reports.

No model request, no key, no money: every input is already on this machine.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

from genomeos.attribution import readergate as rgate
from genomeos.attribution import targetchoice as tc
from genomeos.predict.enhancer_target import MIN_EFFECT

LANE = "lane-emit7g"
ITEM = "ROADMAP §5 item 7g, the next step: compile.py emits every coding target at its bar"

#: the bar, as the object the compiler already applies. Asserted by `is` in the tests, never copied.
BAR = MIN_EFFECT
BAR_NOTE = (
    "the bar is predict.enhancer_target.MIN_EFFECT, 0.1 |log2 fold change|, the same object the"
    " compiler names its one target with; the size rule is targets.genes_at_bar's, the larger of the"
    " predicted fall and the predicted rise over every track, ties to the fall, strongest first, which"
    " is the rule predict_target names the compiled target with. Nothing here is tuned or chosen"
)
#: the 24 generated programs every census of this series reads
COMPILED_DIR = Path("data/knowledge/compiled")
CHROMS = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y"])
POPULATION = (
    "the 24 compiled non-coding programs as they stand (data/knowledge/compiled/noncoding_<chrom>.bio,"
    " git-ignored, pinned by sha256 in this census), which are the files M1 ran and every census of"
    " this series reads. Per program: every `element` block, every `region` block, and every `rule`"
    " line split into predicted and experimental by readergate.rule_targets"
)
EMISSION_NOW = (
    "one `targets:` property and one `rule` line per attributed element, both naming"
    " predicted_coding['gene'] (compile.compile_chromosome's attributed-elements section). An element"
    " of the experimental layer (`<id>_measured`) carries no predicted rule and a `region` block"
    " carries none either"
)
EMISSION_AFTER = (
    "one `rule` line per candidate target of the same elements: targetchoice.candidate_list(heads,"
    " responses.at_bar(chrom, element, MIN_EFFECT), coding) -- the compiled target first, then every"
    " protein-coding gene of the response cache's window at the bar, strongest first, deduplicated. The"
    " population of elements is unchanged, so no element without a compiled target gains one"
)
#: the fields a per-candidate rule line needs, and where each comes from in the cache record
RULE_FIELDS = (
    "gene: the candidate's symbol",
    "action: activates when the larger side is the predicted fall, inhibits when it is the predicted"
    " rise (compile.compiled_action's own reading of the sign, unchanged)",
    "log2: the signed change targets.genes_at_bar returns, which is max_drop_log2fc for a fall and"
    " max_rise_log2fc for a rise",
    "tissue: max_drop_tissue for a fall and max_rise_tissue for a rise, the cache record's own field",
    "strength: min(1.0, |log2|) rounded to 3 places, as the compiled rule already states it",
)

# ---- what a program states as it stands --------------------------------------------------------------

BLOCKS = ("module", "domain", "gene", "region", "element")


def program_counts(text: str) -> dict[str, int]:
    """The blocks and rule lines of one compiled program as it stands, by kind.

    `rules_predicted` is what the emission would replace; `rules_experimental` is what an assay
    measured and is never touched. `elements_with_a_predicted_rule` is the population the emission
    reaches, and `elements_without_one` is counted beside it so a zero cannot come from not looking.
    """
    out: Counter[str] = Counter()
    elements: set[str] = set()
    for ln in text.splitlines():
        head = ln.split(" ", 1)[0] if ln and not ln.startswith((" ", "#")) else ""
        if head in BLOCKS:
            out[f"{head}_blocks"] += 1
            if head == "element":
                elements.add(ln.split(" ")[1])
    with_rule: set[str] = set()
    for element, _gene, kind in rgate.rule_targets(text):
        out[f"rules_{kind}"] += 1
        if kind == "predicted":
            with_rule.add(element)
    out["elements_with_a_predicted_rule"] = len(with_rule)
    out["elements_without_a_predicted_rule"] = len(elements - with_rule)
    out["rules_total"] = out["rules_predicted"] + out["rules_experimental"]
    return {k: out.get(k, 0) for k in sorted({*out, *(f"{b}_blocks" for b in BLOCKS)})}


#: every count `targetchoice.chromosome_choice` can raise, stated even at 0: a zero that comes from not
#: looking and a zero that was looked for read the same in a Counter, and this project has been bitten
#: twice by the difference
CHOICE_KEYS = (
    "elements",
    "experimental_rule_lines_removed",
    "compiled_targets_1",
    "compiled_targets_2_or_more",
    "not_cached",
    "no_coding_gene_at_the_bar",
    "compiled_target_not_first_at_the_bar",
    "compiled_target_not_at_the_bar",
)

# ---- what the emission would write -------------------------------------------------------------------


def candidate_rules(
    responses: Any, chrom: str, element_id: str, heads: Sequence[str], coding: set[str]
) -> list[dict[str, Any]]:
    """One row per rule the emission would write at this element, in the order it would write them.

    The head rows carry no cache reading of their own: `compile.py` already states the compiled
    target's own number from `predicted_coding`, and this census neither re-reads nor overwrites it.
    A row beyond the heads carries the cache's signed change and the track it came from, so a later
    `compile.py` change can write a rule line from this row alone.
    """
    bar = responses.at_bar(chrom, element_id, BAR)
    names = tc.candidate_list(heads, bar, coding)
    record = responses.element(chrom, element_id) if bar is not None else None
    by_gene = {g["gene"]: g for g in (record or {}).get("genes") or [] if g.get("gene")}
    rows: list[dict[str, Any]] = []
    signed = dict(bar or ())
    for i, gene in enumerate(names):
        row: dict[str, Any] = {"gene": gene, "rank": i, "head": gene in set(heads)}
        value = signed.get(gene)
        if value is not None:
            g = by_gene.get(gene, {})
            fall = value < 0
            row |= {
                "log2": value,
                "action": "activates" if fall else "inhibits",
                "tissue": g.get("max_drop_tissue") if fall else g.get("max_rise_tissue"),
                "strength": round(min(1.0, abs(float(value))), 3),
            }
        rows.append(row)
    return rows


def chromosome_census(chrom: str, text: str, coding: set[str], responses: Any) -> dict[str, Any]:
    """One chromosome's census: what the program states now and what the emission would state."""
    now = program_counts(text)
    look = tc.chromosome_choice(chrom, text, coding, responses)
    sizes = [len(c) for c in look.candidates.values()]
    heads = {e: len({lk.symbol for lk in look.look.links[e]}) for e in look.candidates}
    after = sum(sizes)
    gained = sum(1 for e, c in look.candidates.items() if len(c) > heads[e])
    return {
        "now": now,
        "choice_counts": {k: look.counts.get(k, 0) for k in sorted({*look.counts, *CHOICE_KEYS})},
        "elements_reached": len(look.candidates),
        "rules_predicted_now": sum(heads.values()),
        "rules_predicted_after": after,
        "rules_added": after - sum(heads.values()),
        "elements_gaining_a_target": gained,
        "elements_keeping_one_target": len(sizes) - gained,
        "elements_gaining_a_target_that_had_none": 0,
        "targets_per_element_after": tc.size_histogram(sizes),
        "largest_candidate_list": max(sizes, default=0),
    }


SUM_KEYS = (
    "elements_reached",
    "rules_predicted_now",
    "rules_predicted_after",
    "rules_added",
    "elements_gaining_a_target",
    "elements_keeping_one_target",
    "elements_gaining_a_target_that_had_none",
)


def combine(per_chrom: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """The genome-wide totals of `chromosome_census` rows, every count added and none recomputed."""
    rows = list(per_chrom.values())
    out: dict[str, Any] = {k: sum(r[k] for r in rows) for k in SUM_KEYS}
    out["now"] = _add_dicts(r["now"] for r in rows)
    out["choice_counts"] = _add_dicts(r["choice_counts"] for r in rows)
    out["targets_per_element_after"] = _add_dicts(r["targets_per_element_after"] for r in rows)
    out["largest_candidate_list"] = max((r["largest_candidate_list"] for r in rows), default=0)
    out["chromosomes"] = len(rows)
    return out


def _add_dicts(ds: Iterable[dict[str, int]]) -> dict[str, int]:
    total: Counter[str] = Counter()
    for d in ds:
        total.update(d)
    return dict(sorted(total.items()))
