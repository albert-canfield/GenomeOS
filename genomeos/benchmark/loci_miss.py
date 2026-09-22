# SPDX-License-Identifier: AGPL-3.0-or-later
"""What a miss means: overlap and distance classes over every scored frame, registered first.

The loci benchmark scores a hit when the gene a derived layer ranks first is one of the published
target's symbols, and a miss otherwise. That rule has never distinguished two things that are not
the same finding:

- the model named a gene **somewhere else**, which is the failure the benchmark exists to count;
- the model named **the right place under another name** - an antisense, host, readthrough or
  microRNA transcript whose body lies inside the published target's - which is a fact about what
  GENCODE calls a locus and not about what the model predicted.

Section 24 of docs/LOCI-BENCHMARK.md found the second case by accident. Of the 21 loci where the
repaired `loci.read_deletion` changed the deletion layer's answer, 19 are locatable and 5 promote a
gene whose body **overlaps** the published target's - IGF2-AS over IGF2, ENSG00000259006 over MC1R,
ENSG00000288879 over HMGA1, ENSG00000240739 over SLC2A3, CCDC26 over CCDC26 - with 7 more inside
100 kb. Every one is a miss today. The same phenomenon broke the fourth frame's positive control,
and the coordinator's decision of 2026-09-22 (docs/ROADMAP.md item 4) ruled that section 22 stays
interpretable, that **the hit rule stays symbol equality everywhere so that no rate moves**, and
that the general question comes here with its own registration rather than a ruling there.

**This module describes existing misses. It does not change the hit rule and it computes no rate
that replaces one.** Whether the benchmark should ever report a second, overlap-tolerant rate beside
the strict one is a question with a number attached, and the number is fixed in `PREREGISTRATION`
below before a single miss is classified.

**The trap, stated plainly: a benchmark that widens what counts as a hit after seeing its rate is a
benchmark tuning itself.** Every rate in this benchmark has fallen under scrutiny - 0.200 at n = 50,
then 0.180, then 0.172 - and a session that has just watched its headline fall and then discovers a
rule under which several misses become hits is in exactly the position where a benchmark stops
measuring anything. Three things are done about it, all before any counting: the classes are defined
in GENCODE coordinates and committed first; the thresholds that would make the case for a second
rate, and the thresholds that would close the question, are committed with them; and any tolerant
rate is reported BESIDE the strict one, labelled a description, and never instead of it.

**Cost: zero AlphaGenome requests, by construction.** Every reading classified here is already on
disk - the frames' saved results and `loci_reread`'s re-read of them - and GENCODE coordinates are a
local annotation. Nothing in this module calls a scorer, opens a client or reads an interval the
sweep has not already bought.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from genomeos.benchmark import loci_reread
from genomeos.benchmark.loci import RESULTS_DIR
from genomeos.results import load_result

#: the frames whose misses are classified, with the doc section each is reported in. The first five
#: are `loci_reread`'s frames, re-read through the repaired reader; the sixth is the corrected draw,
#: which was scored through the repaired reader in the first place.
FRAMES: dict[str, str] = {
    "loci_benchmark": "9 - the seventeen panel",
    "loci_candidates": "19 - the nine candidates",
    "loci_third": "21 - the third set",
    "loci_fourth": "22 - the fourth frame, drawn by a rule",
    "loci_noncoding": "23 - the non-coding frame",
    "loci_fourth_rejoined": "24 - the fourth frame, the corrected draw",
}

#: The classes, in the order they are tested. The first that matches is the gene's class, so the
#: order is part of the definition and is fixed here with the classes themselves.
CLASSES: list[dict[str, str]] = [
    {
        "name": "exact",
        "definition": (
            "the gene the layer ranks first IS one of the published target symbols. This is a hit"
            " under the benchmark's rule, and it is in the list as the classifier's control: every"
            " strict hit must land here and a classifier that puts one anywhere else is broken"
        ),
    },
    {
        "name": "overlaps_the_target_body",
        "definition": (
            "a different symbol whose GENCODE gene body intersects a published target's gene body"
            " on the same chromosome: start < target.end and end > target.start, the same test"
            " `loci_reread.neighbourhood` uses. This is the class - and the ONLY class - that can"
            " feed the tolerant rate, because it is the one where the model named the published"
            " target's own span under a name the panel did not write down"
        ),
    },
    {
        "name": "within_10_kb",
        "definition": "no overlap; the smallest body-to-body gap to a published target is 1 to 10,000 bp",
    },
    {
        "name": "within_100_kb",
        "definition": "no overlap; the smallest body-to-body gap is 10,001 to 100,000 bp",
    },
    {
        "name": "same_ctcf_node",
        "definition": (
            "further than 100 kb, but the named gene and a published target are both inside the"
            " element's CTCF node, read off the `node` layer's own `genes_in_node` list as that"
            " layer stored it. It is tested after the distance classes so that a near gene is"
            " reported as near rather than as a node member"
        ),
    },
    {
        "name": "elsewhere",
        "definition": (
            "on the same chromosome, further than 100 kb, and not in the element's CTCF node. This"
            " is the class that means the benchmark's miss is a real miss"
        ),
    },
    {
        "name": "not_locatable",
        "definition": (
            "the named symbol has no GENCODE record on the locus's chromosome, or GENCODE uses it"
            " at more than one locus there - Y_RNA is used at 57 loci on chr2 - so a distance to"
            " whichever copy comes first is arithmetic about the wrong gene. Counted and listed"
            " apart, never folded into any other class and never into a tolerant rate"
        ),
    },
    {
        "name": "named_nothing",
        "definition": "the layer ranked no gene first, so there is no gene to place",
    },
]

PREREGISTRATION: dict[str, Any] = {
    "written": (
        "2026-09-22, and committed on its own before a single miss was classified. The counts in"
        " section 24 that motivate it - 5 overlaps and 7 more inside 100 kb among 19 locatable"
        " movements - are over the 21 loci where the REPAIRED READER CHANGED THE ANSWER, which is a"
        " different and much smaller set than the misses classified here; they are the reason to"
        " measure and not a preview of the result"
    ),
    "the_question": (
        "the benchmark scores a miss whenever the gene a derived layer ranks first is not one of"
        " the published target's symbols, and it has never measured whether a miss means `named a"
        " gene somewhere else` or `named the right place under another name`"
    ),
    "what_this_is_not": (
        "it is a DESCRIPTION of misses that already exist and NOT a new hit rule. `loci.score_locus`"
        " is not touched, no saved result is rewritten, and no rate published in sections 9, 18, 19,"
        " 21, 22, 23 or 24 moves by this work. The coordinator's decision of 2026-09-22 already"
        " fixed that the hit rule stays symbol equality everywhere; this module is bound by it"
    ),
    "what_is_classified": (
        "every locus in every frame scored through the repaired reader that is a MISS under the"
        " strict rule - `score.target_hit_derived` false - in the six frames of `FRAMES`. The"
        " strict HITS are classified too and are the control: each must come out `exact`"
    ),
    "which_gene_is_classified": (
        "the benchmark's headline hit is a union over the derived layers, so the locus's class is"
        " the BEST class, in `CLASSES` order, over the rank-1 gene of each derived layer present"
        " (deletion, eqtl, gene_input, coding). Taking the best is what makes the locus-level class"
        " comparable to a union hit; the per-layer classes are reported beside it, and the deletion"
        " layer's own class is reported separately because that layer is the model's answer about"
        " the element and is where section 24's movements happened"
    ),
    "the_classes": [c["name"] for c in CLASSES],
    "the_classes_are_fixed_first": (
        "the eight classes and the order they are tested in are defined in `CLASSES` above, in"
        " GENCODE coordinates, and committed before any locus is classified. Nothing is added,"
        " merged or re-cut after the counts are seen"
    ),
    "only_overlap_can_be_tolerated": (
        "the tolerant rate, if it is ever reported, counts a locus when a derived layer's rank-1"
        " gene is `exact` OR `overlaps_the_target_body`. The distance classes are DESCRIPTION ONLY"
        " and feed no rate, ever. The reason is not taste: a rule that admitted `within_10_kb`"
        " would score the nearest-gene baseline as a hit, and that baseline is the control the"
        " whole fourth frame was drawn to beat. A benchmark whose hit rule tolerates distance"
        " cannot then report that the model follows proximity"
    ),
    "what_would_make_the_case_for_a_second_rate": (
        "BOTH of: (a) pooled over the six frames, at least 15% of strict misses classify as"
        " `overlaps_the_target_body`; and (b) the overlap-tolerant rate exceeds the strict rate by"
        " at least 0.05 in at least one frame of n >= 50 - the fourth frame or the corrected draw -"
        " which is at least three loci there. (a) alone would be a handful of curated loci; (b)"
        " alone could be three loci out of a pooled distribution that is otherwise clean. A second"
        " permanently reported rate has to earn both"
    ),
    "what_would_close_the_question": (
        "BOTH of: (a) pooled, fewer than 5% of strict misses classify as"
        " `overlaps_the_target_body`; and (b) no frame's overlap-tolerant rate exceeds its strict"
        " rate by more than one locus. That outcome says the misses are genuinely elsewhere, the"
        " strict rule is right, and section 24's five overlaps were a property of the movement set"
        " rather than of the benchmark. It is to be reported exactly as plainly as the other"
        " outcome, and it STRENGTHENS the benchmark rather than costing it anything"
    ),
    "in_between": (
        "between 5% and 15%, or either threshold met without the other: the distribution is"
        " published, NO second rate is adopted, and the question is recorded as measured but not"
        " settled. It would then be re-asked only on a frame drawn for it, not on these"
    ),
    "the_trap": (
        "a benchmark that widens what counts as a hit after seeing its rate is a benchmark tuning"
        " itself. This session can see, before it starts, that a tolerant rule would raise numbers"
        " that have fallen three times in two days. That is why the classes are fixed before the"
        " counting, why the thresholds are fixed with them, why only overlap - never distance - can"
        " ever be tolerated, and why any tolerant rate is reported BESIDE the strict rate, labelled"
        " a description, and never in place of it. If the thresholds are met, the output of this"
        " work is a proposal to the benchmark's owner carrying its own registration, not a rate"
        " change made here"
    ),
    "the_classifier_control": (
        "every locus that is a strict hit must classify as `exact`, in the locus-level class and in"
        " the class of the layer that scored the hit. A hit that classifies as anything else means"
        " the symbol lookup, the chromosome join or the class order is wrong and NOTHING else in"
        " the output may be read. `classify_all` reports the count and the harness raises"
    ),
    "falsifiers": [
        "a strict hit classifies as anything but `exact` - the classifier is broken, stop",
        (
            "a locus classifies as `overlaps_the_target_body` at a layer whose named gene is on a"
            " different chromosome from the published target - the join is wrong"
        ),
        (
            "the counts per frame do not sum to that frame's loci - misses plus hits - so some"
            " locus was silently dropped"
        ),
        (
            "the strict rate recomputed from the classified rows differs from the rate section 24"
            " publishes for that frame - this module is reading a different scoring than the"
            " benchmark published, and the classification is about something else"
        ),
        (
            "any AlphaGenome request is spent. Every reading is already on disk and GENCODE is"
            " local; a request means this module is scoring rather than describing"
        ),
    ],
    "cost": "0 AlphaGenome requests, checked before anything was written rather than estimated after",
}

#: what section 24 publishes for each frame under the repaired reader, as (k, n) on the frame's own
#: headline denominator. Falsifier 4: the strict rate recomputed here must reproduce it exactly.
PUBLISHED_STRICT: dict[str, tuple[int, int]] = {
    "loci_benchmark": (15, 17),
    "loci_candidates": (7, 9),
    "loci_third": (7, 9),
    "loci_fourth": (9, 50),
    "loci_noncoding": (1, 2),
    "loci_fourth_rejoined": (10, 58),
}

#: the frames whose headline is the graded rate - drawn loci only, and only where the published
#: target lies inside the model's input. The curated panels report over every locus they chose.
GRADED_FRAMES = ("loci_fourth", "loci_fourth_rejoined")

_ORDER = [c["name"] for c in CLASSES]

# the two classes that an overlap-tolerant rate would count, and nothing else, ever.
TOLERATED = ("exact", "overlaps_the_target_body")

DERIVED_LAYERS = ("deletion", "eqtl", "gene_input", "coding")


def _annotation(chrom: str, cache: dict[str, Any]) -> Any:
    """GENCODE gene models for one chromosome. Local, free, and cached across frames."""
    if chrom not in cache:
        from genomeos.genome import Annotation, default_gencode

        gff = default_gencode({chrom})
        if gff is None:
            raise FileNotFoundError(f"{chrom} has no GENCODE annotation fetched")
        cache[chrom] = Annotation.from_gff3(gff, {chrom})
    return cache[chrom]


def _record(ann: Any, name: str) -> tuple[Any, int]:
    """The GENCODE record for a named gene, and how many loci carry that symbol.

    The lookup is `loci_reread.neighbourhood`'s, deliberately: id key first, then a scan by symbol,
    so the two agree on the five overlaps section 24 published. A symbol GENCODE uses more than once
    on the chromosome has no single body, which is what the second return value is for.
    """
    copies = [g for g in ann.genes.values() if g.symbol == name]
    rec = ann.genes.get(name) or (copies[0] if copies else None)
    return rec, len(copies)


def classify_gene(ann: Any, named: str | None, targets: list[str], node_genes: list[str]) -> dict[str, Any]:
    """Where one named gene sits relative to the published targets, in the registered classes."""
    if not named:
        return {"named": None, "class": "named_nothing", "gap_bp": None, "nearest_target": None}
    if named in targets:
        return {"named": named, "class": "exact", "gap_bp": 0, "nearest_target": named}
    rec, copies = _record(ann, named)
    if rec is None or copies > 1:
        return {
            "named": named,
            "class": "not_locatable",
            "gap_bp": None,
            "nearest_target": None,
            "why": (
                "no GENCODE record on this chromosome" if rec is None else f"{copies} copies of the symbol"
            ),
        }
    gap: int | None = None
    nearest: str | None = None
    overlapping: str | None = None
    for t in targets:
        trec, tcopies = _record(ann, t)
        if trec is None or tcopies > 1:
            continue
        if rec.locus.start < trec.locus.end and rec.locus.end > trec.locus.start:
            overlapping, gap, nearest = t, 0, t
            break
        d = min(abs(rec.locus.start - trec.locus.end), abs(trec.locus.start - rec.locus.end))
        if gap is None or d < gap:
            gap, nearest = d, t
    if overlapping is not None:
        return {
            "named": named,
            "class": "overlaps_the_target_body",
            "gap_bp": 0,
            "nearest_target": overlapping,
        }
    if gap is None:
        return {
            "named": named,
            "class": "not_locatable",
            "gap_bp": None,
            "nearest_target": None,
            "why": "no published target has a single GENCODE body on this chromosome",
        }
    if gap <= 10_000:
        cls = "within_10_kb"
    elif gap <= 100_000:
        cls = "within_100_kb"
    elif named in node_genes and any(t in node_genes for t in targets):
        cls = "same_ctcf_node"
    else:
        cls = "elsewhere"
    return {"named": named, "class": cls, "gap_bp": gap, "nearest_target": nearest}


def classify_locus(row: dict[str, Any], cache: dict[str, Any]) -> dict[str, Any]:
    """One scored locus: the class of each derived layer's rank-1 gene, and the locus's own class.

    The locus class is the best of the layer classes in `CLASSES` order, because the benchmark's
    headline hit is a union over the derived layers.
    """
    chrom = row["expected"]["chrom"]
    ann = _annotation(chrom, cache)
    targets = list(row["expected"]["targets"])
    node_genes = list((row["readings"].get("node") or {}).get("genes_in_node") or [])
    by_layer = row["score"]["scored"]["target"]["by_layer"]
    layers: dict[str, Any] = {}
    for layer in DERIVED_LAYERS:
        v = by_layer.get(layer)
        if v is None or v.get("provenance") != "derived":
            continue
        named = (v.get("named") or [None])[0]
        layers[layer] = {**classify_gene(ann, named, targets, node_genes), "hit": bool(v.get("hit"))}
    ranked = [(_ORDER.index(v["class"]), layer) for layer, v in layers.items()]
    best = min(ranked)[1] if ranked else None
    # The nearest-coding-TSS-in-node rule gets the SAME description applied to it. It is not a
    # derived layer and takes no part in the locus class; it is here because a tolerance applied to
    # the model and not to the baseline it is measured against would flatter the model by
    # construction, and this is the one addition that can only make the tolerant reading look worse.
    node = by_layer.get("node") or {}
    node_named = (node.get("named") or [None])[0]
    node_class = classify_gene(ann, node_named, targets, node_genes)
    return {
        "locus": row["locus"],
        "chrom": chrom,
        "targets": targets,
        "node_class": node_class["class"],
        "node_named": node_named,
        "node_hit": bool(node.get("hit")),
        "strict_hit": bool(row["score"]["target_hit_derived"]),
        "class": layers[best]["class"] if best else "named_nothing",
        "class_from_layer": best,
        "deletion_class": (layers.get("deletion") or {}).get("class"),
        "deletion_named": (layers.get("deletion") or {}).get("named"),
        "deletion_gap_bp": (layers.get("deletion") or {}).get("gap_bp"),
        "gap_bp": layers[best]["gap_bp"] if best else None,
        "unaskable": bool(row["score"].get("deletion_unaskable")),
        "layers": layers,
    }


def _headline(frame: str, classed: list[dict[str, Any]]) -> dict[str, Any]:
    """The frame's own published denominator, the strict rate over it, and the tolerant description."""
    keep = classed
    if frame in GRADED_FRAMES:
        keep = [c for c in classed if not c["locus"].startswith("CONTROL_") and not c["unaskable"]]
    strict = sum(1 for c in keep if c["strict_hit"])
    tolerant = sum(1 for c in keep if c["strict_hit"] or c["class"] in TOLERATED)
    n = len(keep)
    return {
        "denominator": (
            "the drawn loci where the published target is inside the model's input"
            if frame in GRADED_FRAMES
            else "every locus the panel chose"
        ),
        "n": n,
        "strict_k": strict,
        "strict_rate": round(strict / n, 3) if n else None,
        "tolerant_k": tolerant,
        "tolerant_rate": round(tolerant / n, 3) if n else None,
        "loci_it_would_add": [c["locus"] for c in keep if not c["strict_hit"] and c["class"] in TOLERATED],
        "baseline_strict_k": sum(1 for c in keep if c["node_hit"]),
        "baseline_tolerant_k": sum(1 for c in keep if c["node_hit"] or c["node_class"] in TOLERATED),
        "baseline_reading": (
            "the nearest-coding-TSS-in-node rule under the same two readings. A tolerance applied"
            " to the model and not to the baseline it is measured against would flatter the model"
            " by construction, so it is applied to both or the comparison says nothing"
        ),
        "reading": (
            "the tolerant column is a DESCRIPTION of what the strict rate would read if a gene body"
            " overlapping the published target's counted. It is not this benchmark's rate and is"
            " never reported in place of the strict one"
        ),
    }


def classify_frame(
    frame: str, rows: list[dict[str, Any]], cache: dict[str, Any] | None = None
) -> dict[str, Any]:
    """Every locus of one frame, classified, with the misses tallied and the hits used as a control."""
    cache = {} if cache is None else cache
    classed = [classify_locus(r, cache) for r in rows]
    misses = [c for c in classed if not c["strict_hit"]]
    hits = [c for c in classed if c["strict_hit"]]
    tally = {name: 0 for name in _ORDER}
    for c in misses:
        tally[c["class"]] += 1
    del_tally = {name: 0 for name in _ORDER}
    for c in misses:
        if c["deletion_class"]:
            del_tally[c["deletion_class"]] += 1
    broken = [c["locus"] for c in hits if c["class"] != "exact"]
    return {
        "frame": frame,
        "section": FRAMES.get(frame, "?"),
        "loci": len(classed),
        "hits": len(hits),
        "misses": len(misses),
        "miss_classes": tally,
        "miss_classes_deletion_layer_only": del_tally,
        "headline": _headline(frame, classed),
        "control_hits_that_are_not_exact": broken,
        "overlapping_misses": [
            {
                "locus": c["locus"],
                "named": (c["layers"].get(c["class_from_layer"]) or {}).get("named"),
                "layer": c["class_from_layer"],
                "targets": c["targets"],
            }
            for c in misses
            if c["class"] == "overlaps_the_target_body"
        ],
        "per_locus": classed,
    }


def rows_for(frame: str, results_dir: Path) -> list[dict[str, Any]]:
    """The frame's loci AS THE REPAIRED READER READS THEM.

    The five frames of `loci_reread.FRAMES` were scored before the fix, so they are re-read through
    `loci_reread.reread_frame`, which refuses to return anything unless the frozen pre-fix loop
    reproduces every stored reading field for field. The corrected draw was scored through the
    repaired reader in the first place and is loaded as it stands. No request either way.
    """
    if frame in loci_reread.FRAMES:
        return loci_reread.reread_frame(frame, results_dir)["loci_after"]
    out = load_result(frame, results_dir)
    if not out:
        raise FileNotFoundError(f"no {frame} result to classify")
    return list(out["loci"])


def _share_without(frames: dict[str, Any], drop: str) -> float | None:
    """The pooled overlap share with one frame left out, so the pool's dependence is visible."""
    misses = sum(v["misses"] for f, v in frames.items() if f != drop)
    over = sum(v["miss_classes"]["overlaps_the_target_body"] for f, v in frames.items() if f != drop)
    return round(over / misses, 3) if misses else None


def classify_all(results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """Every miss in every frame scored through the repaired reader, sorted into the fixed classes."""
    cache: dict[str, Any] = {}
    frames: dict[str, Any] = {}
    for frame in FRAMES:
        frames[frame] = classify_frame(frame, rows_for(frame, results_dir), cache)
    broken = {
        f: v["control_hits_that_are_not_exact"]
        for f, v in frames.items()
        if v["control_hits_that_are_not_exact"]
    }
    if broken:
        raise AssertionError(
            f"the classifier control failed: strict hits that did not classify as an exact symbol"
            f" match: {broken}. Nothing else in this output may be read"
        )
    off = {
        f: (v["headline"]["strict_k"], v["headline"]["n"])
        for f, v in frames.items()
        if f in PUBLISHED_STRICT and (v["headline"]["strict_k"], v["headline"]["n"]) != PUBLISHED_STRICT[f]
    }
    if off:
        raise AssertionError(
            f"the strict rate recomputed here does not reproduce what section 24 published: {off}."
            f" This is describing a different scoring than the benchmark's and may not be read"
        )
    pooled = {name: 0 for name in _ORDER}
    pooled_del = {name: 0 for name in _ORDER}
    misses = 0
    for v in frames.values():
        misses += v["misses"]
        for name in _ORDER:
            pooled[name] += v["miss_classes"][name]
            pooled_del[name] += v["miss_classes_deletion_layer_only"][name]
    share = round(pooled["overlaps_the_target_body"] / misses, 3) if misses else None
    moved = {f: v["headline"]["tolerant_k"] - v["headline"]["strict_k"] for f, v in frames.items()}
    big = [f for f in ("loci_fourth", "loci_fourth_rejoined") if frames[f]["headline"]["n"] >= 50]
    rate_gap = max(
        (frames[f]["headline"]["tolerant_rate"] - frames[f]["headline"]["strict_rate"] for f in big),
        default=0.0,
    )
    case_for = bool(share is not None and share >= 0.15 and rate_gap >= 0.05)
    closes = bool(share is not None and share < 0.05 and max(moved.values(), default=0) <= 1)
    return {
        "result": "loci_miss",
        "preregistration": PREREGISTRATION,
        "classes": CLASSES,
        "frames": frames,
        "pooled": {
            "misses": misses,
            "classes": pooled,
            "classes_deletion_layer_only": pooled_del,
            "overlap_share_of_misses": share,
            "largest_tolerant_rate_gap_at_n_over_50": round(rate_gap, 3),
            "tolerant_minus_strict_per_frame": moved,
            # the two fourth-frame entries are two readings of one draw sharing most of their loci,
            # so the pool is not six independent frames and the share is reported with each of them
            # dropped in turn. Disclosure, not a re-cut: the verdict above stands as registered.
            "overlap_share_dropping_one_fourth_frame_reading": {
                f: _share_without(frames, f) for f in ("loci_fourth", "loci_fourth_rejoined") if f in frames
            },
            "one_fewer_overlap_would_read": (
                round((pooled["overlaps_the_target_body"] - 1) / misses, 3) if misses else None
            ),
        },
        "verdict": {
            "case_for_a_second_rate": case_for,
            "question_closed": closes,
            "reading": (
                "measured against the thresholds committed in `PREREGISTRATION` before any locus was"
                " classified: a second rate needs at least a 15% overlap share AND a gap of at least"
                " 0.05 at n >= 50; the question closes on an overlap share under 5% AND no frame"
                " moving by more than one locus; anything else publishes the distribution and adopts"
                " nothing"
            ),
        },
        "requests_spent": 0,
        "note": (
            "a description of misses that already exist, not a hit rule. `loci.score_locus` is"
            " untouched, no saved result is rewritten, and no published rate moves. The tolerant"
            " column is what the strict rate WOULD read if a gene body overlapping the published"
            " target's counted, and it is reported beside the strict rate and never instead of it"
        ),
    }
