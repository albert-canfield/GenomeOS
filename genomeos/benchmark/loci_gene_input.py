"""The gene-input layer's reader, versioned: v1 as it has always read, v2 repaired (2026-10-03).

`loci.read_gene_input` is the gene-input layer: every already-scored element inside a locus's
window votes, and the votes are summed per gene. It has read ``predicted_coding`` and nothing
else since it was written, and that is the one deletion-derived reader the repair of 2026-09-21
deliberately did not touch. Both places that left it undone say so in those words:

* docs/LOCI-BENCHMARK.md section 24, "Left undone": "**`read_gene_input` still reads
  `predicted_coding` and nothing else.** Section 23 registered that as a property of a separate
  layer - a summed-window reading rather than an element answer - and it was deliberately left out
  of this change so that two movements did not sit behind one number. It is the next thing anyone
  looking at this should measure."
* docs/LOCI-BENCHMARK.md section 25, "Left undone", which is the item this module answers:
  "**`read_gene_input` still reads `predicted_coding` alone**, as section 24 left it. Four of the
  seven overlaps in the fourth frame come from the eQTL layer, so the layer that has not been
  repaired is also the one contributing least to this reading."

Nothing here recomputes a rate. `REPAIR` states what v2 does and why; `credit_rows` is the one
function the two readings differ in; and `loci.read_gene_input` keeps `GENE_INPUT_V1` as its
default so that every committed frame reads exactly as it did. The re-scoring that would turn v2
into a number is registered separately in `loci_gene_input_plan` and is not run by this module.

The shape is the one this file already uses for a second reading, and it is taken rather than
invented. `read_gene_input`'s own docstring says a `window_reading` "is added beside the layer and
does not replace it", and `read_deletion` carries its pre-repair answer beside the repaired one as
`coding_first_target` and `coding_first_targets` so that "every rate this benchmark published
before that date stays computable from the shipped reader rather than only from the git history".
v2 does both: the repaired ranking is the layer's answer, and the v1 ranking travels beside it
under `coding_head_*`.
"""

from __future__ import annotations

from typing import Any

#: The reading as the layer has always had it: each element votes once, for its ``predicted_coding``
#: head, and an element whose coding head names no gene does not vote at all. This is the DEFAULT in
#: `loci.read_gene_input` and its output is byte-identical to the pre-2026-10-03 reader, which
#: `tests/test_loci_gene_input.py` checks against the committed frames rather than against the code.
GENE_INPUT_V1 = "coding_head_only"

#: The repaired reading, opt-in and named. Each element enters BOTH stored heads into one ranking -
#: ``predicted_coding`` first, then ``predicted`` - and counts once towards a gene however many of
#: its heads name that gene. This is `read_deletion`'s repair of 2026-09-21 applied to the summed
#: window, and the walk order is its walk order so that an exact tie keeps the v1 answer.
GENE_INPUT_V2 = "both_keys_ranked"

READINGS = (GENE_INPUT_V1, GENE_INPUT_V2)

#: the two stored heads of one deletion, in the order v2 walks them. Coding first, so that on an
#: exact tie the stable sort keeps the gene v1 would have named.
HEADS = ("predicted_coding", "predicted")

REPAIR: dict[str, Any] = {
    "written": (
        "2026-10-03, before `loci.read_gene_input` was touched. This is a registration of a READER"
        " CHANGE and not of a result: no rate is recomputed by this module, no committed result,"
        " registration or figure is rewritten, and no locus is scored"
    ),
    "the_item": (
        "docs/LOCI-BENCHMARK.md section 25, Left undone, verbatim: '`read_gene_input` still reads"
        " `predicted_coding` alone, as section 24 left it. Four of the seven overlaps in the fourth"
        " frame come from the eQTL layer, so the layer that has not been repaired is also the one"
        " contributing least to this reading.' Section 24's Left undone adds: 'It is the next thing"
        " anyone looking at this should measure'"
    ),
    "what_v1_is": (
        "each already-scored element inside the window contributes abs(log2_fold_change) to the one"
        " gene its `predicted_coding` head names, and contributes nothing when that head names no"
        " gene. v1 is the DEFAULT and is byte-identical to the reader as it stood before this"
        " module existed. A v1 reading carries NO version key, by requirement: stamping one would"
        " change the dict and the committed frames would no longer read identically. The absence of"
        " the `reading` key IS v1"
    ),
    "what_v2_does": (
        "for every already-scored element inside the window both stored heads are entered into the"
        " same summed ranking - `predicted_coding` first, then `predicted` - and an element counts"
        " once towards a gene however many of its heads name that gene, so an element whose"
        " strongest effect happens to be on a coding gene does not tally double. The layer's answer"
        " is the gene with the largest summed abs(log2_fold_change) whatever its biotype. The v1"
        " ranking is computed alongside and returned in the same dict as `coding_head_genes`,"
        " `coding_head_target` and `coding_head_rank_of_first_published_target`, so every rate this"
        " benchmark published under v1 stays computable from the shipped reader"
    ),
    "why_this_reading_and_not_the_favourable_one": (
        "the argument is `loci_reread.PREREGISTRATION`'s and is not restated here with different"
        " words: the layer's question does not depend on the moved gene's biotype, and"
        " `predicted_coding` is 'a restriction of the same computation rather than a better"
        " estimate of it' (`loci.read_deletion`'s docstring, 2026-09-21). What that passage ALSO"
        " says is the reason this change was held back, and it is answered rather than ignored:"
        " 'folding a second layer into this change would put two movements behind one number'."
        " That was an argument about COMBINING the two repairs in one re-read, not an argument that"
        " this layer should stay restricted. It is answered by versioning instead of combining -"
        " the deletion layer's repair is already committed and already measured in section 24, so"
        " v2 here is the SECOND movement arriving on its own, separable by construction because v1"
        " remains the default and both readings are returned side by side from one call. Nothing"
        " requires the two to be read together, and this module recomputes neither"
    ),
    "direction_this_can_move_the_rate": (
        "the same asymmetry `loci_reread` registered applies, for the same construction reason."
        " abs(predicted) >= abs(predicted_coding) holds per element per head by construction, so a"
        " coding gene that ranked first under v1 can be DISPLACED under v2 and can never be"
        " promoted; only a published non-coding target can gain. The fourth frame's drawing rule"
        " requires the published target to be protein coding, so in that frame v2 can lower the"
        " gene-input layer's rate and cannot raise it. Written here before any re-score, as the"
        " expected direction and not as a measurement"
    ),
    "what_this_does_not_touch": (
        "`loci.score_locus`'s hit rule stays symbol equality. `read_deletion`, `read_eqtl` and"
        " every other layer are untouched. `summed_by_cell` and `predicted_coding_by_cell` are"
        " untouched, being coding-specific by name. `_gene_input_window` - the `window_reading`"
        " added beside the layer - keeps reading the compact coding head for elements the response"
        " cache does not hold, because that is a different reading with its own registration"
    ),
    "the_eqtl_sentence_is_about_a_layer_this_repair_cannot_reach": (
        "section 25's item puts the four eQTL overlaps beside the unrepaired reader, and the two"
        " are not connected by this change. `loci.read_eqtl` ranks GTEx v8 significant"
        " single-tissue cis-eQTL pairs by p-value and never reads `predicted` or"
        " `predicted_coding` at all, so NO version of `read_gene_input` can move an eQTL-layer"
        " answer. Checked against the committed `loci_miss` result rather than asserted: of the"
        " fourth frame's seven `overlapping_misses`, four are layer `eqtl` and three are layer"
        " `deletion`, and NONE is layer `gene_input`. The sentence's own reading - that the"
        " unrepaired layer is the one contributing least - is what the result shows, and it"
        " contributes nothing at all to those seven"
    ),
    "the_governing_rules_are_imported_not_restated": (
        "the overlap and nearest-gene rules belong to `loci_miss.PREREGISTRATION` and"
        " `loci_miss.CLASSES`, committed in bf33233. Anything here that needs a class calls"
        " `loci_miss.classify_gene`. Its prohibition travels with it and binds this module and any"
        " plan of this module's: the distance classes are DESCRIPTION ONLY and feed no rate, ever,"
        " because 'a rule that admitted `within_10_kb` would score the nearest-gene baseline as a"
        " hit, and that baseline is the control the whole fourth frame was drawn to beat'"
    ),
    "cost": "0 AlphaGenome requests, no network, no money. Every input is already on disk",
}


def credit_rows(rows: list[dict[str, Any]], reading: str = GENE_INPUT_V1) -> dict[str, dict[str, Any]]:
    """The per-gene summed credit over one window's already-scored elements, under one reading.

    The ONE function the two readings differ in, so that v1's identity is a property of a single
    code path and not of two copies that have to be kept in step. Returns the accumulator keyed by
    gene symbol, in insertion order, exactly as `loci.read_gene_input` built it inline before this
    module existed; the caller sorts it.

    v1 walks `predicted_coding` alone. v2 walks `HEADS` in order with a per-element set of genes
    already credited, which is `read_deletion`'s repaired walk applied to the summed window.
    """
    if reading not in READINGS:
        raise ValueError(f"unknown gene-input reading {reading!r}; one of {READINGS}")
    by_gene: dict[str, dict[str, Any]] = {}
    heads = (HEADS[0],) if reading == GENE_INPUT_V1 else HEADS
    for e in rows:
        counted: set[str] = set()
        for key in heads:
            p = e.get(key) or {}
            if not p.get("gene"):
                continue
            # an element counts once towards a gene however many of its heads name it, or the
            # element tallies double wherever the strongest effect happens to be on a coding gene.
            if p["gene"] in counted:
                continue
            counted.add(p["gene"])
            g = by_gene.setdefault(
                p["gene"], {"gene": p["gene"], "elements": 0, "activating": 0, "summed": 0.0}
            )
            g["elements"] += 1
            g["activating"] += int(p["action"] == "activates")
            g["summed"] = round(g["summed"] + abs(p["log2_fold_change"]), 3)
    return by_gene


def rank(by_gene: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    """The accumulator sorted by summed credit, descending. Stable, so `credit_rows`'s insertion
    order breaks a tie - and under v2 the coding head of an element is inserted first, which makes
    the tie-break the one that preserves the v1 answer."""
    return sorted(by_gene.values(), key=lambda g: -g["summed"])


def rank_of_first_target(ranked: list[dict[str, Any]], targets: Any) -> int | None:
    """The 1-based rank of the first published target in a ranking, or None. `loci`'s expression,
    moved here so v1 and v2 cannot drift apart in it."""
    return next((i + 1 for i, g in enumerate(ranked) if g["gene"] in targets), None)
