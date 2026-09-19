# SPDX-License-Identifier: AGPL-3.0-or-later
"""The 2,006 Mb the budget never cut into a block: what it is, and whether the offset lives there.

`panel_background.py` measured what the panel's GC- and replication-timing-matched background holds
and found that composition explains at most half of why every non-coding tier reads above it: exons,
conserved elements and promoter proximity together account for 41% (fossil), 56% (regulatory) and 25%
(neutral), leaving 44% to 75% standing, with segmental duplication pulling the other way. It then
named its own leftover: **2,006 Mb of the 2,808 Mb background -- 71% -- lies outside every unknown
block**, reading 6.39 recurring events per kilobase against the tiers' 7.07 to 7.40, and nobody had
measured what that sequence is. This module measures it.

The split is not invented. `genome/unknown.unknown_blocks` cuts a chromosome into blocks with exactly
one rule::

    cursor = 0
    for every gene in the annotation, by start:
        if gene.start - cursor >= min_size:  a block is cut from cursor to gene.start
        cursor = max(cursor, gene.end)
    if length - cursor >= min_size:          a last block is cut

So a base is never cut into a block for exactly one of three reasons, and they are exhaustive:

- it lies inside the span of an annotated gene -- **any** gene, since `Annotation` loads `gene`,
  `ncRNA_gene` and `pseudogene` alike, so introns, UTRs, lncRNA bodies and pseudogene bodies are all
  skipped along with the coding exons the rule is aimed at;
- it lies in an inter-gene gap shorter than `min_size` (1,000 bases), ordinary intergenic sequence
  that the rule declines only because of its length;
- it lies after the last gene in a tail shorter than `min_size`.

`replay_the_cut` walks that loop with the bookkeeping the original discards, and the blocks it
recovers are checked against the committed budget block by block before anything is read from them.
The gene spans are then split by what they are -- CDS, UTR and intron of coding genes, non-coding
gene bodies, pseudogene bodies -- so that a biological class and a rule of the instrument can be told
apart.

The question the leftover poses is then asked the one way that is not a story, the way the previous
lane asked it: the background is rebuilt with each class taken out of it through the panel's own
`build_background`, and the tier ratios are recomputed. One arm of that decomposition is degenerate
and says so in the result rather than in prose: taking out *every* never-cut base leaves a background
made of blocks, which is a tier average and not a control, and a tier average is not evidence about a
tier. The comparison that carries the weight is `compare.standardised` over 100-base windows, where
both arms are real groups of sequence and coverage is held apart from the covariates.
"""

from __future__ import annotations

from array import array
from bisect import bisect_left
from collections import Counter
from pathlib import Path
from typing import Any

from genomeos.attribution.human_panel import (
    BACKGROUND_BIN,
    CACHE,
    GC_STRATA,
    MIN_RECURRING,
    TIERS,
    merge_intervals,
    rt_over,
    rt_stratum,
)
from genomeos.attribution.panel_background import (
    CELL,
    CLAIM_WITH_COVERAGE,
    CLAIM_WITHOUT_COVERAGE,
    COVERAGE_EDGES,
    COVERAGE_MEASURE,
    TSS_EDGES,
    WINDOW_CAP,
    Rebuilt,
    background_bins,
    composition,
    cover_cells,
    coverage_over,
    explained,
    informing,
    nearest,
    rebuild,
    span_bp,
    tier_ratios,
)
from genomeos.compare import Strata, input_presence, standardised

MIN_SIZE = 1_000  # genome.unknown.unknown_blocks' own default: the length the budget's blocks are cut with
DOMINANT_MIN = BACKGROUND_BIN // 2  # a kilobase belongs to the class holding at least half of it
DEGENERATE_MIN = 0.95  # a background this nearly made of blocks is a tier average, not a control
CARRIES_MOST = 0.5  # a class that moves this share of the offset carries it
LIKE_THE_REST = 0.5  # the length-rule class reading this much of the gene bodies' gap: the rule, not biology
# Below this share of never-cut bases declined for their LENGTH rather than for being inside a gene,
# the cutting rule and the biological class are the same predicate and no measurement over these
# bases can tell them apart. It is a property of the rule, not of the sample, so more chromosomes
# will not rescue it.
SEPARABLE_MIN = 0.05

# The classes, in the priority a base is assigned by. Priority matters only where annotations nest --
# a lncRNA inside a coding intron, a pseudogene inside a lncRNA -- and how many bases it decides is
# reported per class as `bases_given_to_a_higher_priority_class` rather than left to the reader.
CLASSES = (
    "coding_gene_cds",
    "coding_gene_utr",
    "coding_gene_intron",
    "noncoding_gene_body",
    "pseudogene_body",
    "intergenic_gap_under_min_size",
    "intergenic_tail_under_min_size",
)
GENE_BODY = CLASSES[:5]
INTERGENIC = CLASSES[5:]
# the one class that is ordinary intergenic sequence declined for its length alone: the discriminator
# between "a biological class is quiet" and "the cutting rule selects quiet sequence"
LENGTH_RULE_ONLY = "intergenic_gap_under_min_size"
UNACCOUNTED = "outside_every_block_and_unaccounted"
IN_A_BLOCK = "in_a_block"  # the window's arm when the budget did cut it: not a class, the other side
COVERAGE_FULL = 0.999  # at or above this the whole panel informed the window
RARE = INTERGENIC  # the classes taken whole rather than subsampled: they are the discriminator

WHY_NOT_CUT = {
    "coding_gene_cds": "inside an annotated gene span: the cursor skips the whole gene, coding exon",
    "coding_gene_utr": "inside an annotated gene span: the cursor skips the whole gene, untranslated exon",
    "coding_gene_intron": "inside an annotated gene span: the cursor skips the whole gene, intron",
    "noncoding_gene_body": "inside an annotated gene span: Annotation loads ncRNA_gene as a gene",
    "pseudogene_body": "inside an annotated gene span: Annotation loads pseudogene as a gene",
    "intergenic_gap_under_min_size": f"inter-gene gap shorter than min_size ({MIN_SIZE} bases)",
    "intergenic_tail_under_min_size": f"tail after the last gene, shorter than min_size ({MIN_SIZE} bases)",
    UNACCOUNTED: "not in a block and in none of the seven classes: the cut and the annotation disagree",
}

SOURCES = {
    "blocks": "the committed budget_<chrom> blocks, cut by genome.unknown.unknown_blocks",
    "gene_model": "GENCODE 50 (GRCh38.p14): gene, ncRNA_gene and pseudogene spans, exons and CDS",
    "background": "human_panel.build_background, rebuilt from the local store, nothing re-derived",
    "coverage": COVERAGE_MEASURE,
}

# The verdicts this lane is allowed to reach, as data rather than as prose. A fifth does not exist,
# and none may be paraphrased into "nothing is quiet" or "everything is explained".
VERDICT_SAME_PARTITION = (
    "the cutting rule and the biological class are the same predicate on all but a fraction of the "
    "never-cut bases, because the budget declines a base for being inside an annotated gene span, so "
    "these data cannot assign the offset to one rather than the other; what they do say is that the "
    "panel's background is mostly gene-body sequence and that gene-body sequence carries fewer "
    "recurring events than block sequence at matched GC, timing, distance to a coding TSS and coverage"
)
VERDICT_INSTRUMENT = (
    "the offset tracks sequence the budget declined to cut into a block, not a biological class: "
    "ordinary intergenic sequence declined for its length alone is as quiet as the gene bodies"
)
VERDICT_CLASS = "the offset tracks a biological class of the never-cut sequence"
VERDICT_NEITHER = "neither a class of the never-cut sequence nor the cutting rule accounts for it"
VERDICT_UNASSESSED = "the discriminating class was not assessed, so no reading is available"

# The verdict is about MASS: which reading the 2,006 Mb can be assigned to. The mechanism is a
# different question on a different sample -- the handful of megabases where the two readings come
# apart -- and it gets its own sentence so that a small, clear result cannot be quoted as a large one.
MECHANISM_GENE_BODY = (
    "on the bases where the two readings come apart -- ordinary intergenic sequence the budget "
    "declined for its length alone -- the sequence is NOT quiet, so the quietness travels with the "
    "gene body rather than with the decision not to cut; the sample is a fraction of the never-cut "
    "sequence and gene-proximal by construction, and the sentence is about the mechanism, not the mass"
)
MECHANISM_THE_CUT = (
    "on the bases where the two readings come apart -- ordinary intergenic sequence the budget "
    "declined for its length alone -- the sequence IS quiet, so the quietness travels with the "
    "decision not to cut rather than with the gene body; the sample is a fraction of the never-cut "
    "sequence and gene-proximal by construction, and the sentence is about the mechanism, not the mass"
)
MECHANISM_UNASSESSED = "too few windows are declined for their length alone to read the mechanism"
MECHANISM_MIN_WINDOWS = 500  # below this the length-rule arm is reported without a mechanism sentence
# A threshold picked by hand can decide a sentence by a hair, and a reader cannot see that from the
# sentence. Within this many points of the line the mechanism reading is ON THE LINE, and the result
# says so beside the sentence rather than leaving the margin to be recomputed by whoever doubts it.
MECHANISM_ON_THE_LINE = 0.005
ON_THE_LINE = (
    "the mechanism reading is decided by a margin of {margin} points and would flip if the threshold "
    "moved: it is the direction the numbers point, not a result that survives its own threshold"
)

NO_COVERAGE = "no kilobase of the group carries a coverage value"
NOT_ASSESSED = "not assessed"


def verdict(
    separable: bool,
    blocks_above_the_never_cut: float | None,
    length_rule_against_the_rest: float | None,
    length_rule_windows: int = 0,
) -> dict[str, Any]:
    """Which sentences the numbers allow -- one about the mass, one about the mechanism.

    `blocks_above_the_never_cut` is how much more often a block window carries a recurring event than
    a never-cut window with the same GC, timing, distance to a coding TSS and coverage: the gap the
    cut opens. `length_rule_against_the_rest` is the same quantity for ordinary intergenic sequence
    declined for its length alone, against the rest of the never-cut. Sitting near zero means that
    sequence behaves like a gene body although it is not one, which is the cut and not the class;
    sitting near the first means it behaves like a block, which is the class and not the cut.

    `separable` decides the first sentence, and it is a property of the cutting rule rather than of
    the sample: where the budget declines a base for being inside a gene span, "never cut" and "gene
    body" are the same predicate and no sample size separates them.
    """
    gap = blocks_above_the_never_cut
    line = None if gap is None else LIKE_THE_REST * gap
    facts: dict[str, Any] = {
        "separable": separable,
        "blocks_above_the_never_cut": gap,
        "length_rule_against_the_rest_of_the_never_cut": length_rule_against_the_rest,
        "length_rule_windows": length_rule_windows,
        "like_the_rest_threshold": LIKE_THE_REST,
        "mechanism_min_windows": MECHANISM_MIN_WINDOWS,
        # where the length-rule sequence sits on the line from the never-cut sequence (0) to the
        # blocks (1); the threshold above is a cut on exactly this number
        "length_rule_share_of_the_block_gap": (
            round(length_rule_against_the_rest / gap, 4)
            if gap and length_rule_against_the_rest is not None
            else None
        ),
        "mechanism_margin": (
            round(length_rule_against_the_rest - line, 5)
            if line is not None and length_rule_against_the_rest is not None
            else None
        ),
        "mechanism_on_the_line_threshold": MECHANISM_ON_THE_LINE,
    }
    if line is None or length_rule_against_the_rest is None or length_rule_windows < MECHANISM_MIN_WINDOWS:
        facts["mechanism"] = MECHANISM_UNASSESSED
        facts["mechanism_on_the_line"] = False
    else:
        facts["mechanism"] = (
            MECHANISM_GENE_BODY if length_rule_against_the_rest >= line else MECHANISM_THE_CUT
        )
        margin = facts["mechanism_margin"]
        facts["mechanism_on_the_line"] = abs(margin) < MECHANISM_ON_THE_LINE
        if facts["mechanism_on_the_line"]:
            facts["mechanism_caveat"] = ON_THE_LINE.format(margin=abs(margin))
    if blocks_above_the_never_cut is None:
        return {"verdict": VERDICT_UNASSESSED, **facts}
    if blocks_above_the_never_cut <= 0:
        return {"verdict": VERDICT_NEITHER, **facts}
    if not separable:
        return {"verdict": VERDICT_SAME_PARTITION, **facts}
    if length_rule_against_the_rest is None:
        return {"verdict": VERDICT_UNASSESSED, **facts}
    if length_rule_against_the_rest <= line:
        return {"verdict": VERDICT_INSTRUMENT, **facts}
    return {"verdict": VERDICT_CLASS, **facts}


# ================================================================================================
# interval arithmetic, so that the classes are disjoint by construction rather than by hope
# ================================================================================================
def subtract(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Every base of `a` that is not in `b`. Both sides are merged first, so no base counts twice."""
    a = merge_intervals(list(a))
    b = merge_intervals(list(b))
    ends = [e for _, e in b]
    out: list[tuple[int, int]] = []
    for s, e in a:
        cur = s
        i = bisect_left(ends, s + 1)
        while i < len(b) and b[i][0] < e:
            if b[i][0] > cur:
                out.append((cur, b[i][0]))
            cur = max(cur, b[i][1])
            if cur >= e:
                break
            i += 1
        if cur < e:
            out.append((cur, e))
    return out


def bases(intervals: list[tuple[int, int]]) -> int:
    return sum(e - s for s, e in merge_intervals(list(intervals)))


# ================================================================================================
# the cut, replayed with the bookkeeping the original discards
# ================================================================================================
def replay_the_cut(
    spans: list[tuple[int, int]], length: int, min_size: int = MIN_SIZE
) -> dict[str, list[tuple[int, int]]]:
    """`unknown_blocks`' loop, keeping what it throws away: the blocks, and why the rest was declined.

    `spans` are the annotated gene loci sorted by start, exactly as the original sorts them. The four
    outputs partition [0, length) by construction: a base is in a block, inside a gene span, in a gap
    under `min_size`, or in a tail under `min_size`.
    """
    blocks: list[tuple[int, int]] = []
    gaps: list[tuple[int, int]] = []
    cursor = 0
    for s, e in spans:
        if s - cursor >= min_size:
            blocks.append((cursor, s))
        elif s > cursor:
            gaps.append((cursor, s))
        cursor = max(cursor, e)
    tail: list[tuple[int, int]] = []
    if length - cursor >= min_size:
        blocks.append((cursor, length))
    elif length > cursor:
        tail.append((cursor, length))
    return {
        "blocks": sorted(blocks),
        "gene_spans": merge_intervals(spans),
        "intergenic_gap_under_min_size": merge_intervals(gaps),
        "intergenic_tail_under_min_size": merge_intervals(tail),
    }


def gene_model(chrom: str, gff: Path | None = None) -> dict[str, Any]:
    """Gene spans by kind, the exons and CDS of coding genes, and the coding TSS, from GENCODE.

    `gff` is the annotation to read; the default is the one `genome.annotation` would pick for the
    chromosome, and a caller passes its own only in a test with no GENCODE on disk.
    """
    from genomeos.genome.annotation import Annotation, default_gencode

    gff = gff or default_gencode({chrom})
    if gff is None:
        return {}
    ann = Annotation.from_gff3(gff, {chrom})
    spans: list[tuple[int, int]] = []
    coding: list[tuple[int, int]] = []
    noncoding: list[tuple[int, int]] = []
    pseudo: list[tuple[int, int]] = []
    exons: list[tuple[int, int]] = []
    cds: list[tuple[int, int]] = []
    tss: list[int] = []
    for g in sorted(ann.genes.values(), key=lambda g: g.locus.start):
        if g.locus.chrom != chrom:
            continue
        iv = (g.locus.start, g.locus.end)
        spans.append(iv)
        if g.type == "protein_coding":
            coding.append(iv)
            for t in g.transcripts.values():
                exons += [(x.start, x.end) for x in t.exons]
                cds += [(c.start, c.end) for c, _ in t.cds]
            tss.append(g.locus.start if "MINUS" not in str(g.locus.strand) else g.locus.end)
        elif "pseudogene" in g.type:
            pseudo.append(iv)
        else:
            noncoding.append(iv)
    return {
        "genes": len(spans),
        "spans_in_start_order": spans,
        "coding_spans": merge_intervals(coding),
        "noncoding_spans": merge_intervals(noncoding),
        "pseudogene_spans": merge_intervals(pseudo),
        "coding_exons": merge_intervals(exons),
        "coding_cds": merge_intervals(cds),
        "tss": sorted(tss),
    }


def classify(model: dict[str, Any], cut: dict[str, list[tuple[int, int]]], length: int) -> dict[str, Any]:
    """The seven disjoint classes of never-cut sequence, plus whatever falls outside all of them.

    Disjointness is enforced by subtraction in the order `CLASSES` declares, and the bases each
    subtraction took from a lower-priority class are reported so that a nested annotation is a number
    rather than a silent choice.
    """
    coding = model["coding_spans"]
    exons = model["coding_exons"]
    cds = model["coding_cds"]
    noncoding = model["noncoding_spans"]
    pseudo = model["pseudogene_spans"]
    ivs: dict[str, list[tuple[int, int]]] = {
        "coding_gene_cds": cds,
        # exonic on some transcript of the gene and coding on none of them: UTR, and the exons of a
        # coding gene's non-coding transcripts, which are untranslated for the same reason
        "coding_gene_utr": subtract(exons, cds),
        "coding_gene_intron": subtract(coding, exons),
        "noncoding_gene_body": subtract(noncoding, coding),
        "pseudogene_body": subtract(pseudo, coding + noncoding),
        "intergenic_gap_under_min_size": cut["intergenic_gap_under_min_size"],
        "intergenic_tail_under_min_size": cut["intergenic_tail_under_min_size"],
    }
    taken = {
        "noncoding_gene_body": bases(noncoding) - bases(ivs["noncoding_gene_body"]),
        "pseudogene_body": bases(pseudo) - bases(ivs["pseudogene_body"]),
    }
    covered = merge_intervals([iv for v in ivs.values() for iv in v] + cut["blocks"])
    ivs[UNACCOUNTED] = subtract([(0, length)], covered)
    return {
        "intervals": ivs,
        "bases_given_to_a_higher_priority_class": taken,
        "bases_in_a_block_and_in_a_class": bases(
            subtract([iv for v in ivs.values() for iv in v], subtract([(0, length)], cut["blocks"]))
        ),
    }


def reproduces_the_budget(cut: dict[str, list[tuple[int, int]]], committed: list[dict]) -> dict[str, Any]:
    """Does the replayed loop give back the committed budget's blocks, one for one?

    Nothing below means anything if it does not: the classes are defined as what the cut declined, and
    a cut that is not the budget's declines something else.
    """
    mine = sorted((s, e) for s, e in cut["blocks"])
    theirs = sorted((b["start"], b["end"]) for b in committed)
    a, b = set(mine), set(theirs)
    return {
        "blocks_replayed": len(mine),
        "blocks_committed": len(theirs),
        "identical": mine == theirs,
        "only_in_the_replay": len(a - b),
        "only_in_the_committed_budget": len(b - a),
        "bases_replayed": bases(mine),
        "bases_committed": bases(theirs),
        "min_size": MIN_SIZE,
    }


# ================================================================================================
# the background's kilobases, split by the class that holds them
# ================================================================================================
def class_arrays(rb: Rebuilt, ivs: dict[str, list[tuple[int, int]]]) -> dict[str, array]:
    """A coverage array per class, plus the canonical CDS cut and the tiers `background_bins` needs."""
    cells = rb.length // CELL + 2
    sets: dict[str, list[tuple[int, int]]] = {"cds_canonical": rb.cds_all}
    sets.update(ivs)
    for t in TIERS:
        sets[f"in_{t}"] = rb.tier_intervals(t)
    return {k: cover_cells(v, cells) for k, v in sets.items()}


def dominant(row: dict[str, Any], names: list[str]) -> str:
    """The class holding at least half of a kilobase, or `mixed` when none does."""
    best = max(names, key=lambda n: row[n])
    return best if row[best] >= DOMINANT_MIN else "mixed"


def group(rows: list[dict[str, Any]], names: list[str]) -> dict[str, Any]:
    """`panel_background.composition` with the covariates every group here has to print beside it.

    A group whose kilobases carry no coverage value at all is named rather than given a median of
    nothing: measured-and-absent and never-looked stay apart.
    """
    if not rows:
        return {"kilobases": 0, "bases": 0, "why": "the class has no kilobase in the background"}
    if not any(r["coverage"] is not None for r in rows):
        return {
            "kilobases": len(rows),
            "bases": sum(r["kept"] for r in rows),
            "coverage_median": None,
            "coverage_assessed": 0,
            "why": NO_COVERAGE,
        }
    out = composition(rows, names)
    out["coverage_measure"] = COVERAGE_MEASURE
    return out


def split_the_leftover(rows: list[dict[str, Any]], names: list[str]) -> tuple[dict[str, Any], dict[str, Any]]:
    """The leftover's bases by class, and its kilobases grouped by the class that dominates them.

    Two denominators, both printed, because they are not the same number: `bases_in_the_kilobases` is
    the whole of the kilobases the background counts as outside every block, and `bases_kept` is what
    survives the canonical-CDS cut the background applies -- the 2,006 Mb figure. Class bases are
    counted against the first, since a class is a property of the sequence and not of the cut.
    """
    span = BACKGROUND_BIN * len(rows)
    kept = sum(r["kept"] for r in rows)
    by_class: dict[str, Any] = {}
    for n in names:
        b = sum(r[n] for r in rows)
        by_class[n] = {
            "bases": b,
            "share_of_the_kilobases": round(b / span, 4) if span else None,
            "why_the_budget_did_not_cut_it": WHY_NOT_CUT[n],
        }
    accounted = sum(v["bases"] for v in by_class.values())
    split = {
        "kilobases": len(rows),
        "bases_in_the_kilobases": span,
        "bases_kept": kept,
        "bases_accounted_by_a_class": accounted,
        "bases_in_no_class": span - accounted,
        "by_class": by_class,
    }
    labels = Counter(dominant(r, names) for r in rows)
    groups = {
        label: group([r for r in rows if dominant(r, names) == label], names)
        for label in sorted(labels, key=lambda x: -labels[x])
    }
    return split, groups


# ================================================================================================
# the decomposition: the background rebuilt with each class taken out of it
# ================================================================================================
def separability(split: dict[str, Any]) -> dict[str, Any]:
    """Can the cutting rule and the biological class be told apart on these bases at all?

    The budget declines a base for one of two reasons: it is inside an annotated gene span, or the
    inter-gene gap it sits in is shorter than `min_size`. On the first kind the two candidate
    explanations -- "the budget did not cut it" and "it is a gene body" -- are the SAME predicate, and
    no comparison over those bases can separate them. Only the second kind is ordinary intergenic
    sequence declined for a reason with no biology in it, and its share is what decides whether the
    question is answerable rather than what the answer is.
    """
    by = split["by_class"]
    total = sum(by[n]["bases"] for n in by)
    length_rule = sum(by[n]["bases"] for n in INTERGENIC if n in by)
    share = length_rule / total if total else None
    return {
        "never_cut_bases": total,
        "bases_declined_for_being_inside_a_gene_span": total - length_rule,
        "bases_declined_for_the_length_rule_alone": length_rule,
        "share_that_can_separate_the_two_readings": round(share, 6) if share is not None else None,
        "threshold": SEPARABLE_MIN,
        "separable": bool(share is not None and share >= SEPARABLE_MIN),
        "reads": (
            "the two readings are the same predicate on this chromosome's never-cut bases"
            if share is not None and share < SEPARABLE_MIN
            else "enough sequence is declined for its length alone to tell the two readings apart"
        ),
    }


def arms(ivs: dict[str, list[tuple[int, int]]]) -> dict[str, list[tuple[int, int]]]:
    """Every exclusion this lane makes, each named for what it takes out of the background."""
    out = {f"without_{n}": ivs[n] for n in ivs}
    out["without_every_gene_body"] = merge_intervals([iv for n in GENE_BODY for iv in ivs[n]])
    out["without_the_intergenic_never_cut"] = merge_intervals([iv for n in INTERGENIC for iv in ivs[n]])
    out["without_every_never_cut_base"] = merge_intervals([iv for v in ivs.values() for iv in v])
    return out


def block_share_left(rows: list[dict[str, Any]], names: list[str], taken: list[str]) -> float | None:
    """How much of the background an exclusion leaves is block sequence, at kilobase resolution.

    An exclusion that leaves a background made of blocks has turned the comparison into a tier
    average against a tier, which this project does not accept as a control. The share is the number
    that says so, and `DEGENERATE_MIN` is the line. `rows` is the WHOLE background, not the leftover,
    and the denominator is the kilobase, since a class is counted over whole kilobases.
    """
    span = BACKGROUND_BIN * len(rows)
    if not span:
        return None
    never_cut = sum(r[n] for r in rows for n in names)
    gone = sum(r[n] for r in rows for n in taken)
    left = span - gone
    return round(1 - max(0, never_cut - gone) / left, 4) if left > 0 else None


def decompose(
    rb: Rebuilt, ivs: dict[str, list[tuple[int, int]]], rows: list[dict[str, Any]], names: list[str]
) -> dict[str, Any]:
    """Every tier's ratio with each class taken out of the background, through the panel's own code."""
    taken_by: dict[str, list[str]] = {f"without_{n}": [n] for n in ivs}
    taken_by["without_every_gene_body"] = list(GENE_BODY)
    taken_by["without_the_intergenic_never_cut"] = list(INTERGENIC)
    taken_by["without_every_never_cut_base"] = list(ivs)
    out: dict[str, Any] = {"as_the_panel_builds_it": tier_ratios(rb, rb.cds_all)}
    for name, excluded in arms(ivs).items():
        if not excluded:
            out[name] = {"assessed": False, "why": "the class has no bases on this chromosome"}
            continue
        share = block_share_left(rows, names, taken_by[name])
        degenerate = share is not None and share >= DEGENERATE_MIN
        entry: dict[str, Any] = {
            "assessed": not degenerate,
            "block_share_of_the_remaining_background": share,
            "degenerate": degenerate,
            **tier_ratios(rb, merge_intervals(rb.cds_all + excluded)),
        }
        if degenerate:
            entry["why"] = (
                "the remaining background is almost entirely block sequence, so this arm is a tier "
                "average against a tier and is reported without being counted as evidence"
            )
        out[name] = entry
    return out


# ================================================================================================
# the two-group comparison, where both arms are real groups of sequence
# ================================================================================================
def class_windows(
    rb: Rebuilt, feats: dict[str, array], tss: list[int], names: list[str], cap: int = WINDOW_CAP
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """100-base windows carrying their class, whether they sit in a block, and their coverage.

    Windows with no coverage value are KEPT here, with `coverage: None`, so that `input_presence` can
    count what is present before any coverage stratum is applied to it. They are dropped only when the
    stratified comparison is run, and the number dropped is reported.

    Long chromosomes are subsampled by a fixed step, as `panel_background.window_rows` subsamples
    them, with one exception that is counted rather than hidden: every window of a class in `RARE` is
    taken whole. Those classes are the discriminator this lane turns on -- ordinary intergenic
    sequence the budget declined for its length alone -- and they are a thousandth of the sequence, so
    a fixed step would leave them with no power at all while adding nothing to the other arms.
    """
    panel = rb.panel
    per_block = informing(panel)
    ev_cell: Counter = Counter()
    for ev in rb.evs:
        if ev.minor >= MIN_RECURRING:
            ev_cell[ev.start // CELL] += 1
    cells = sorted(
        {
            c
            for i in range(len(panel.block_start))
            for c in range(panel.block_start[i] // CELL, (panel.block_end[i] - 1) // CELL + 1)
        }
    )
    step = max(1, -(-len(cells) // cap))
    taken = set(cells[::step])
    rare = {c for c in cells if c not in taken and any(feats[n][c] >= CELL / 2 for n in RARE if n in feats)}
    rows: list[dict[str, Any]] = []
    thin = 0
    for c in sorted(taken | rare):
        s, e = c * CELL, (c + 1) * CELL
        cov, aligned = coverage_over(panel, per_block, s, e)
        gc = panel.gc_fraction([(s, e)])
        if aligned < CELL / 2 or gc is None:
            thin += 1
            continue
        held = {n: span_bp(feats[n], s, e) for n in names}
        best = max(held, key=lambda n: held[n])
        # blocks are disjoint, so their bases add; the window is placed where half of it lies
        block_bp = sum(span_bp(feats[f"in_{t}"], s, e) for t in TIERS)
        if block_bp >= CELL / 2:
            placed: str | None = IN_A_BLOCK
        elif held[best] >= CELL / 2:
            placed = best
        else:
            placed = None
        rows.append(
            {
                "start": s,
                "length": aligned,
                "gc": round(gc, 4),
                "rt": rt_stratum(rt_over(rb.rt, [(s, e)]), rb.edges) or 0,
                "tss": nearest(tss, s + CELL // 2),
                "coverage": round(cov, 4) if cov is not None else None,
                "recurring": bool(ev_cell.get(c, 0)),
                # the arm this window belongs to: a block, one class of never-cut sequence, or
                # neither because it is split between them with nothing holding half
                "annotation_class": placed,
                "coding": span_bp(feats["cds_canonical"], s, e) > 0,
            }
        )
    return rows, {
        "step": step,
        "cells_in_panel_blocks": len(cells),
        "cells_taken_by_the_step": len(taken),
        "cells_added_for_a_rare_class": len(rare),
        "rare_classes_taken_whole": list(RARE),
        "windows": len(rows),
        "windows_too_thin_or_without_gc": thin,
        "windows_without_a_coverage_value": sum(1 for r in rows if r["coverage"] is None),
        "windows_in_a_block": sum(1 for r in rows if r["annotation_class"] == IN_A_BLOCK),
        # measured, and placed in neither arm: split between a block and a class, or between classes
        "windows_placed_in_neither_arm": sum(1 for r in rows if r["annotation_class"] is None),
        "coverage_measure": COVERAGE_MEASURE,
    }


LADDER = {
    "gc_and_timing": Strata(length=(CELL - 1,), gc=GC_STRATA, rt=(0.5, 1.5)),
    "gc_timing_and_tss": Strata(length=(CELL - 1,), gc=GC_STRATA, rt=(0.5, 1.5), tss=TSS_EDGES),
    "gc_timing_tss_and_coverage": Strata(
        length=(CELL - 1,), gc=GC_STRATA, rt=(0.5, 1.5), tss=TSS_EDGES, coverage=COVERAGE_EDGES
    ),
}
INPUTS = {
    "the panel's alignment over the window": "coverage",
    "the local annotation of the window": "annotation_class",
}


def compare_groups(targets: list[dict[str, Any]], controls: list[dict[str, Any]]) -> dict[str, Any]:
    """One two-group comparison: presence of the inputs first, then the three stratifications.

    Presence before the stratum, always. A coverage stratum only bites a claim whose input had to be
    bought, and saying which BEFORE running it costs nothing and stops "conditioning changed nothing"
    from reading as a test passed when no test was administered.
    """
    out: dict[str, Any] = {
        "targets": len(targets),
        "controls": len(controls),
        "input_presence": input_presence(targets, controls, INPUTS),
    }
    t = [r for r in targets if r["coverage"] is not None]
    c = [r for r in controls if r["coverage"] is not None]
    out["targets_without_a_coverage_value"] = len(targets) - len(t)
    out["controls_without_a_coverage_value"] = len(controls) - len(c)
    # presence of the input and the level of it are different things, and only the second is what a
    # coverage stratum bites on: a lane where presence is free can still have a level that varies
    out["targets_below_full_coverage"] = sum(1 for r in t if r["coverage"] < COVERAGE_FULL)
    out["controls_below_full_coverage"] = sum(1 for r in c if r["coverage"] < COVERAGE_FULL)
    if not t or not c:
        out["assessed"] = False
        out["why"] = "one of the two arms has no window with a coverage value"
        return out
    out["assessed"] = True
    for name, strata in LADDER.items():
        out[name] = standardised(t, c, strata, hit="recurring")
    return out


def comparisons(rows: list[dict[str, Any]], names: list[str]) -> dict[str, Any]:
    """The instrument's question and the classes' question, both against a control that is sequence.

    `every_block_against_the_never_cut` is the instrument statement: block sequence against the
    sequence the cutting rule declined, holding GC, timing, distance to a coding TSS and coverage.
    Neither arm is a tier average. Each class is then set against the REST of the never-cut sequence,
    so that a class carrying the quietness shows up without a tier appearing in either arm.
    """
    pool = [r for r in rows if not r["coding"]]
    never_cut = [r for r in pool if r["annotation_class"] not in (IN_A_BLOCK, None)]
    out: dict[str, Any] = {
        "coding_windows_excluded": sum(1 for r in rows if r["coding"]),
        "windows_in_a_block": sum(1 for r in pool if r["annotation_class"] == IN_A_BLOCK),
        "windows_never_cut": len(never_cut),
        "windows_placed_in_neither_arm": sum(1 for r in pool if r["annotation_class"] is None),
        "every_block_against_the_never_cut": compare_groups(
            [r for r in pool if r["annotation_class"] == IN_A_BLOCK], never_cut
        ),
    }
    for n in names:
        t = [r for r in never_cut if r["annotation_class"] == n]
        c = [r for r in never_cut if r["annotation_class"] != n]
        out[f"{n}_against_the_rest_of_the_never_cut"] = (
            compare_groups(t, c) if t else {"assessed": False, "why": "no window of this class", "targets": 0}
        )
    return out


# ================================================================================================
# one chromosome
# ================================================================================================
def analyse(
    chrom: str, cache: Path = CACHE, results_dir: Path | None = None, gff: Path | None = None
) -> dict[str, Any]:
    """Split one chromosome's never-cut background by what it is, and ask where the offset lives."""
    rb = rebuild(chrom, cache, results_dir)
    model = gene_model(chrom, gff)
    if not model:
        return {"chrom": chrom, "assessed": False, "why": "GENCODE not read for this chromosome"}
    cut = replay_the_cut(model["spans_in_start_order"], rb.length)
    reproduces = reproduces_the_budget(cut, rb.committed.get("blocks", []) or rb.blocks)
    cls = classify(model, cut, rb.length)
    ivs = cls["intervals"]
    names = [*CLASSES, UNACCOUNTED]
    feats = class_arrays(rb, ivs)
    tss = model["tss"]
    rows, dropped = background_bins(rb, feats, tss)
    leftover = [r for r in rows if not any(r[f"in_{t}"] for t in TIERS)]
    split, groups = split_the_leftover(leftover, names)
    split["separability"] = separability(split)
    win, win_cost = class_windows(rb, feats, tss, names)
    comp = comparisons(win, names)
    held = [
        v["gc_timing_tss_and_coverage"]["targets_matched"]
        for v in comp.values()
        if isinstance(v, dict) and v.get("assessed")
    ]
    offset = decompose(rb, ivs, rows, names)
    return {
        "chrom": chrom,
        "assessed": True,
        "assemblies": rb.panel.n,
        "evidence": SOURCES,
        "claim_available": CLAIM_WITH_COVERAGE if held and all(held) else CLAIM_WITHOUT_COVERAGE,
        "reproduces_the_budget": reproduces,
        "chromosome_bases": rb.length,
        # what could have been measured, before what was: the denominators of everything below
        "eligible": {
            "chromosome_bases": rb.length,
            "kilobases_in_the_panel_alignment": len(rows) + sum(dropped.values()),
            "classes_declared": list(names),
            "classes_with_bases_on_this_chromosome": [n for n in names if ivs[n]],
            "tiers_declared": list(TIERS),
            "tiers_with_blocks_on_this_chromosome": [t for t in TIERS if rb.tier_intervals(t)],
            "windows_in_the_panel_alignment": win_cost["cells_in_panel_blocks"],
        },
        "cut": {
            "min_size": MIN_SIZE,
            "genes_in_the_annotation": model["genes"],
            "block_bases": bases(cut["blocks"]),
            "gene_span_bases": bases(cut["gene_spans"]),
            "gap_bases_under_min_size": bases(cut["intergenic_gap_under_min_size"]),
            "tail_bases_under_min_size": bases(cut["intergenic_tail_under_min_size"]),
            "bases_given_to_a_higher_priority_class": cls["bases_given_to_a_higher_priority_class"],
            "bases_in_a_block_and_in_a_class": cls["bases_in_a_block_and_in_a_class"],
        },
        "background": group(rows, names),
        "leftover": group(leftover, names),
        "split": split,
        "by_dominant_class": groups,
        "offset": offset,
        "explained": {t: explained(offset, t) for t in TIERS},
        "comparisons": comp,
        "coverage": {
            "kilobases_in_panel_blocks": len(rows) + sum(dropped.values()),
            "kilobases_measured": len(rows),
            "kilobases_not_measured": dropped,
            "kilobases_outside_every_block": len(leftover),
            "classes_assessed": [n for n in names if ivs[n]],
            "classes_not_assessed": {n: "read, and empty on this chromosome" for n in names if not ivs[n]},
            "windows": win_cost,
        },
    }
