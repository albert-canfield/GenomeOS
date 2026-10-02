# SPDX-License-Identifier: AGPL-3.0-or-later
"""Which gene does a compiled rule name, and which gene did the deletion answer measure?

A compiled rule names its target by **symbol**: ``rule EH38E2130399 activates TPTE``. A gene symbol
is not unique over a chromosome - GENCODE v50 carries twelve genes called ``Y_RNA`` on chr21 alone -
and the AlphaGenome deletion answer the rule is built from measured **one** gene, the one in the
element's own 1 Mb window. So a symbol that resolves to another locus of the same symbol makes the
rule assert a relation between an element and a gene the model never scored.

This module resolves both sides to an Ensembl gene id and says whether they are the same gene. It
counts and reports. **It establishes about no rule that the rule is wrong, it moves no verdict, and
it deletes, rewrites and relabels no compiled rule.** A count of rules whose two sides resolve
differently is a count of symbol resolutions, not a measure of how often the model's predictions are
wrong; see :func:`limitations`.

What each side is, and why
--------------------------

**The compiled side** is the gene the project's own symbol resolution returns for the rule as
written. That resolver is :func:`genomeos.attribution.pilot_bio.gene_tss`: "Symbol -> GENCODE v50 TSS
on ``chrom``, a protein-coding gene preferred when a symbol repeats", over
:func:`genomeos.attribution.holdout.genes` which is sorted by ``(chrom, tss)``. It is the resolver
behind the distance that lane-notopen's 1,678 came from, so :func:`compiled_locus` reproduces its
choice exactly and carries the ``gene_id`` of the row it picks, which ``gene_tss`` throws away.
A target token that is itself a versionless Ensembl id - the deletion answer names a gene with no
symbol by its id, as GENCODE does - is resolved as that id, which ``gene_tss`` cannot do because it
skips symbol-less genes. That one extension beyond the project's resolver is named here, registered,
and counted on its own.

**The measured side** is the gene the deletion answer measured. **The answer does not record an
Ensembl gene id.** ``data/knowledge/alphagenome/all_elements/<chrom>.json`` records
``predicted_coding`` as ``{"gene": "<name>", "action": ..., "log2_fold_change": ..., "tissue": ...}``
and nothing else about the gene's identity, where ``<name>`` is a symbol, or a bare versionless
Ensembl id when the gene has no symbol. The id is therefore **recovered, not read**, from two things
that *are* recorded:

1. ``data/knowledge/alphagenome/elements/<chrom>.json.gz`` records, per element, every gene the
   scorer read in that element's window, by the same naming convention (``genes_in_window`` and the
   ``genes`` rows). How many times the rule's target name appears in that list is a recorded fact
   about the answer, and when it is not exactly once no id is recoverable from the record.
2. the window the method records - "delete the element as a variant in its 1 Mb window" - taken as
   the element midpoint plus or minus :data:`SCORER_HALF_WINDOW`, the half-window the project already
   fixed at ``not_open_profile.SCORER_HALF_WINDOW``. A gene is a candidate when its **gene body**
   overlaps that interval, because that is what puts a gene in front of the scorer.

Both sides read the same annotation, GENCODE v50 under ``data/reference/gencode_v50_<chrom>.gff3.gz``,
recorded by sha256 in the result's manifest.

The outcome of one rule
-----------------------

Exactly one of :data:`AGREE`, :data:`DIFFER` or an unresolvable cause from :data:`CAUSES`, decided in
the fixed precedence of :func:`classify` and never re-ordered after a count. The three sum to the
population with nothing in a residue.
"""

from __future__ import annotations

import gzip
import re
from collections import Counter
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from genomeos.attribution.not_open_profile import SCORER_HALF_WINDOW
from genomeos.predict.enhancer_target import STRONG_EFFECT

# ---- where everything is read from --------------------------------------------------------------

REFERENCE = Path("data/reference")
#: The compiled programs themselves. Their rule lines are the population: 24 files whose `rule`
#: lines number 440,589, the figure context_evidence_registration.json fixed as the denominator.
COMPILED = Path("data/knowledge/compiled")
#: The per-element record of what the deletion scorer read in each element's window.
WINDOW_CACHE = Path("data/knowledge/alphagenome/elements")

CHROMS = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y"])

#: Half the 1 Mb window the deletion scorer reads, imported from where the project fixed it rather
#: than restated. A gene whose body overlaps the element midpoint plus or minus this could be in
#: front of the scorer; one whose body does not, could not.
HALF_WINDOW = SCORER_HALF_WINDOW
#: The sensitivity variant registered before the counts: AlphaGenome's 1 Mb sequence interval is
#: 2**20 = 1,048,576 bp, so half of it is 524,288 rather than the round 500,000 the project fixed.
#: Every count is reported at HALF_WINDOW and the difference at this half-window is reported beside
#: it, so the figure cannot be read as depending on the rounder number.
SENSITIVITY_HALF_WINDOW = 524_288

# ---- the two outcomes that are a comparison, and every cause that is not ------------------------

AGREE = "agree"
DIFFER = "differ"

#: Why a rule's two sides cannot be compared at all. Each is a statement about the record, never
#: about the rule: a rule counted here is not a rule in doubt, it is a rule this test cannot reach.
CAUSES: dict[str, str] = {
    "target_absent_from_the_chromosome_annotation": (
        "the rule's target token is neither a gene_name nor a versionless gene_id of any gene on the "
        "element's own chromosome in GENCODE v50, so neither side resolves to an id. The annotation "
        "AlphaGenome scored with and GENCODE v50 are not guaranteed to be the same release"
    ),
    "rule_of_the_experimental_layer": (
        "the rule's evidence is a CRISPRi screen, not an AlphaGenome deletion answer, so there is no "
        "deletion answer to compare against. The screen records its own Ensembl gene id "
        "(holdout.Unit.gene_id, `measuredGeneEnsemblId`) and that is a different comparison, reported "
        "separately and never added into these counts"
    ),
    "element_absent_from_the_window_cache": (
        "the element has no per-element record in data/knowledge/alphagenome/elements/<chrom>.json.gz, "
        "so what the scorer read in its window is not on disk and the measured side cannot be "
        "recovered without a request, which this lane does not make"
    ),
    "target_not_in_the_recorded_window_list": (
        "the element's recorded window gene list does not contain the rule's target name at all, so "
        "the record does not say which gene the answer named and no id is recoverable from it"
    ),
    "target_named_more_than_once_in_the_recorded_window_list": (
        "the scorer's own record lists more than one gene under the rule's target name in this "
        "element's window, so the answer's `gene` field does not say which of them it measured and no "
        "id is recoverable from the record"
    ),
    "no_annotated_locus_of_that_name_overlaps_the_scorer_window": (
        "the name is in the recorded window list, but no GENCODE v50 gene of that name on this "
        "chromosome has a gene body overlapping the scorer window, so the measured gene cannot be "
        "placed in the annotation both sides are resolved against"
    ),
    "two_or_more_annotated_loci_of_that_name_in_the_scorer_window": (
        "two or more GENCODE v50 genes of that name have a gene body inside the scorer window, so the "
        "window does not pick one out and no id is recoverable"
    ),
}

OUTCOMES = (AGREE, DIFFER, *CAUSES)

MIS_RESOLUTION = (
    "a rule whose target token resolves, by the project's own symbol resolution, to a different "
    "Ensembl gene id from the gene the deletion answer measured. It says the rule names a gene the "
    "model never scored; it does NOT say the predicted relation is wrong, and nothing here "
    "establishes that any one rule is wrong"
)

# ---- the annotation ----------------------------------------------------------------------------


@dataclass(frozen=True)
class GeneRow:
    """One GENCODE v50 gene: its id without a version, its name, its type and its span."""

    gene_id: str
    symbol: str
    gene_type: str
    start: int  # 0-based, inclusive
    end: int  # 0-based, exclusive
    tss: int  # 0-based

    @property
    def coding(self) -> bool:
        return self.gene_type == "protein_coding"


def annotation_path(chrom: str, reference: Path = REFERENCE) -> Path:
    return reference / f"gencode_v50_{chrom}.gff3.gz"


def annotation(chrom: str, reference: Path = REFERENCE) -> tuple[GeneRow, ...]:
    """Every GENCODE v50 gene of `chrom`, sorted by TSS as `holdout.genes` sorts them.

    The TSS is the convention the project already uses everywhere: the start of a `+` gene and the
    end of a `-` gene, both 0-based (`holdout.genes`, `element_types.gene_starts`).
    """
    from genomeos.genome.annotation import iter_gff3

    path = annotation_path(chrom, reference)
    if not path.exists():
        return ()
    rows = []
    for row in iter_gff3(path, types={"gene"}):
        gid = row.attrs.get("gene_id", "").split(".")[0]
        rows.append(
            GeneRow(
                gene_id=gid,
                symbol=row.attrs.get("gene_name", ""),
                gene_type=row.attrs.get("gene_type", ""),
                start=row.start - 1,
                end=row.end,
                tss=row.end - 1 if row.strand == "-" else row.start - 1,
            )
        )
    rows.sort(key=lambda g: g.tss)
    return tuple(rows)


def by_symbol(rows: Iterable[GeneRow]) -> dict[str, tuple[GeneRow, ...]]:
    out: dict[str, list[GeneRow]] = {}
    for g in rows:
        if g.symbol:
            out.setdefault(g.symbol, []).append(g)
    return {k: tuple(v) for k, v in out.items()}


def by_id(rows: Iterable[GeneRow]) -> dict[str, GeneRow]:
    return {g.gene_id: g for g in rows}


# ---- the compiled side: the project's own symbol resolution, with the id carried ----------------

COMPILED_RESOLUTION = (
    "genomeos.attribution.pilot_bio.gene_tss's own choice, with the gene_id it discards carried: "
    "among the GENCODE v50 genes of the element's chromosome whose gene_name is the rule's target "
    "token, a protein_coding gene is preferred, and among equal types the lowest TSS wins, because "
    "holdout.genes is sorted by (chrom, tss) and gene_tss keeps the first of each symbol. A token "
    "that is itself a versionless gene_id of that chromosome resolves to that gene directly, which "
    "gene_tss cannot do because it skips genes with no symbol; those rules are counted on their own."
)


def compiled_locus(
    token: str, symbols: Mapping[str, tuple[GeneRow, ...]], ids: Mapping[str, GeneRow]
) -> GeneRow | None:
    """The gene the project's own symbol resolution returns for a rule's target token, or None."""
    direct = ids.get(token)
    if direct is not None:
        return direct
    cand = symbols.get(token)
    if not cand:
        return None
    return min(cand, key=lambda g: (0 if g.coding else 1, g.tss))


# ---- the measured side: the gene the deletion answer could have measured ------------------------

MEASURED_RESOLUTION = (
    "the deletion answer records its target as a name and NO Ensembl gene id, so the id is recovered "
    "from two recorded things and never read: the element's own recorded window gene list "
    "(data/knowledge/alphagenome/elements/<chrom>.json.gz), which must name the target exactly once, "
    "and the window the method records, taken as the element midpoint plus or minus "
    f"{HALF_WINDOW} (not_open_profile.SCORER_HALF_WINDOW), inside which the gene's BODY must lie. "
    "Exactly one GENCODE v50 gene of that name with a body overlapping the window is the measured "
    "gene; none or several is an unresolvable cause and never a guess. A target token that is a "
    "versionless gene_id is read directly, because there the answer does record the id."
)


def window_candidates(
    token: str,
    midpoint: int,
    symbols: Mapping[str, tuple[GeneRow, ...]],
    half_window: int = HALF_WINDOW,
) -> tuple[GeneRow, ...]:
    """The genes of that name whose body overlaps the scorer's window around the element."""
    lo, hi = midpoint - half_window, midpoint + half_window
    return tuple(g for g in symbols.get(token, ()) if g.end > lo and g.start < hi)


# ---- a compiled rule, as the compiled program states it -----------------------------------------

SOURCE_PREDICTED = "predicted_alphagenome_deletion"
SOURCE_MEASURED = "measured_crispri"

#: `strength` on a compiled rule line is the observation's magnitude in its own unit, |log2 fold
#: change| clipped at 1 for a predicted rule. Clipping at 1 cannot move a value across 0.3, so
#: banding the rule's own `strength` at STRONG_EFFECT gives the same two bands as the cached
#: prediction's own `strength` field, which is what lane-notopen banded by.
EFFECT_BAND_CALL = (
    f"the rule line's own `strength:`, banded at predict.enhancer_target.STRONG_EFFECT = "
    f"{STRONG_EFFECT}: strong at or above, weak below. Clipping at 1 cannot move a value across "
    f"0.3, so this is the same band the cached prediction's `strength` field carries."
)

_RULE = re.compile(r"^rule (\S+) (\S+) (\S+) \{(.*)\}\s*$")
_ELEMENT = re.compile(r"^element (\S+) \{$")
_LOCUS = re.compile(r"^  locus: (\S+):(\d+)-(\d+)$")
_ACTIVITY = re.compile(r"^  activity: (.*)$")
_STRENGTH = re.compile(r"strength: ([0-9.]+)")
_WHEN = re.compile(r"when: cell_type = ([^;}\s]+)")


@dataclass(frozen=True)
class Rule:
    """One compiled rule line, with the element block's locus and activity axis beside it."""

    element: str
    chrom: str
    start: int
    end: int
    verb: str
    target: str
    cell: str
    strength: float | None
    activity: str

    @property
    def source(self) -> str:
        return SOURCE_MEASURED if self.element.endswith("_measured") else SOURCE_PREDICTED

    @property
    def midpoint(self) -> int:
        return (self.start + self.end) // 2

    @property
    def effect_band(self) -> str | None:
        if self.strength is None:
            return None
        return "strong" if self.strength >= STRONG_EFFECT else "weak"


def program_path(chrom: str, compiled: Path = COMPILED) -> Path:
    return compiled / f"noncoding_{chrom}.bio"


def rules(chrom: str, compiled: Path = COMPILED) -> Iterator[Rule]:
    """Every `rule` line of a compiled program, in the order the program states them.

    The element block always precedes its rules, so the locus and the `activity:` axis are taken
    from the block of the id the rule names. A rule naming an element with no block is an error and
    raises, because a rule this function could not locate would be a rule silently dropped from the
    population.
    """
    path = program_path(chrom, compiled)
    blocks: dict[str, tuple[str, int, int, str]] = {}
    current: str | None = None
    locus: tuple[str, int, int] | None = None
    activity = ""
    with path.open() as fh:
        for line in fh:
            if m := _ELEMENT.match(line.rstrip("\n")):
                current, locus, activity = m.group(1), None, ""
                continue
            if current is not None:
                if m := _LOCUS.match(line.rstrip("\n")):
                    locus = (m.group(1), int(m.group(2)), int(m.group(3)))
                    continue
                if m := _ACTIVITY.match(line.rstrip("\n")):
                    activity = m.group(1).strip()
                    continue
                if line.startswith("}"):
                    if locus is not None:
                        blocks[current] = (*locus, activity)
                    current, locus, activity = None, None, ""
                    continue
            if m := _RULE.match(line.rstrip("\n")):
                eid, verb, target, body = m.group(1), m.group(2), m.group(3), m.group(4)
                if eid not in blocks:
                    raise ValueError(f"{path}: rule on {eid!r} with no element block before it")
                chrm, start, end, act = blocks[eid]
                s = _STRENGTH.search(body)
                w = _WHEN.search(body)
                yield Rule(
                    element=eid,
                    chrom=chrm,
                    start=start,
                    end=end,
                    verb=verb,
                    target=target,
                    cell=w.group(1) if w else "unknown",
                    strength=float(s.group(1)) if s else None,
                    activity=act,
                )


# ---- what the scorer read in each element's window, as the cache records it ---------------------

_CACHE_TOKEN = re.compile(rb'"id": "(EH38E\d+)"|"gene": "([^"]+)"')


def window_gene_counts(
    chrom: str,
    wanted: Mapping[str, frozenset[str]],
    cache: Path = WINDOW_CACHE,
    chunk: int = 1 << 22,
) -> tuple[dict[str, Counter[str]], int, set[str]]:
    """Per element, how many times each wanted name appears in its recorded window gene list.

    The cache is one JSON object per chromosome and is streamed, never loaded: the file is read in
    chunks and scanned for the two keys that matter, so memory is the size of `wanted` and not the
    size of the cache. Returns the counts, how many element records the file holds, and which of the
    wanted elements it holds a record for - a count of zero for a name is not the same thing as no
    record at all, and the two are separate causes.
    """
    path = cache / f"{chrom}.json.gz"
    out: dict[str, Counter[str]] = {}
    present: set[str] = set()
    seen = 0
    if not path.exists():
        return out, seen, present
    want_for: frozenset[str] = frozenset()
    current = ""
    with gzip.open(path, "rb") as fh:
        tail = b""
        while True:
            block = fh.read(chunk)
            if not block:
                break
            buf = tail + block
            last = 0
            for m in _CACHE_TOKEN.finditer(buf):
                last = m.end()
                if m.group(1) is not None:
                    current = m.group(1).decode()
                    seen += 1
                    want_for = wanted.get(current, frozenset())
                    if want_for:
                        present.add(current)
                    continue
                if not want_for:
                    continue
                name = m.group(2).decode()
                if name in want_for:
                    out.setdefault(current, Counter())[name] += 1
            tail = buf[last:]
            if len(tail) > (1 << 20):  # no token in a megabyte: nothing can still be mid-match
                tail = tail[-512:]
    return out, seen, present


# ---- the outcome of one rule --------------------------------------------------------------------


@dataclass(frozen=True)
class Outcome:
    """What this test could say about one rule, and the two ids it could resolve."""

    outcome: str  # AGREE, DIFFER, or a key of CAUSES
    compiled_id: str | None
    measured_id: str | None
    token_was_an_id: bool

    @property
    def resolved(self) -> bool:
        return self.outcome in (AGREE, DIFFER)


def classify(
    rule: Rule,
    symbols: Mapping[str, tuple[GeneRow, ...]],
    ids: Mapping[str, GeneRow],
    occurrences: int | None,
    half_window: int = HALF_WINDOW,
) -> Outcome:
    """One rule's outcome, in a precedence fixed before any count was read and never re-ordered.

    `occurrences` is how many times the rule's target name appears in the element's recorded window
    gene list, or None when the element has no record in the cache at all.
    """
    compiled = compiled_locus(rule.target, symbols, ids)
    token_is_id = rule.target in ids
    if compiled is None:
        return Outcome("target_absent_from_the_chromosome_annotation", None, None, token_is_id)
    if rule.source == SOURCE_MEASURED:
        return Outcome("rule_of_the_experimental_layer", compiled.gene_id, None, token_is_id)
    if occurrences is None:
        return Outcome("element_absent_from_the_window_cache", compiled.gene_id, None, token_is_id)
    if occurrences == 0:
        return Outcome("target_not_in_the_recorded_window_list", compiled.gene_id, None, token_is_id)
    if occurrences > 1:
        return Outcome(
            "target_named_more_than_once_in_the_recorded_window_list",
            compiled.gene_id,
            None,
            token_is_id,
        )
    if token_is_id:
        measured: GeneRow | None = ids[rule.target]
    else:
        cand = window_candidates(rule.target, rule.midpoint, symbols, half_window)
        if not cand:
            return Outcome(
                "no_annotated_locus_of_that_name_overlaps_the_scorer_window",
                compiled.gene_id,
                None,
                token_is_id,
            )
        if len(cand) > 1:
            return Outcome(
                "two_or_more_annotated_loci_of_that_name_in_the_scorer_window",
                compiled.gene_id,
                None,
                token_is_id,
            )
        measured = cand[0]
    same = compiled.gene_id == measured.gene_id
    return Outcome(AGREE if same else DIFFER, compiled.gene_id, measured.gene_id, token_is_id)


# ---- the axes the correlation is reported on, all imported ---------------------------------------


def distance_band(distance: int | None) -> str:
    """The project's own bands, imported from executor, with its residual row above where it stops."""
    from genomeos.attribution import executor as ex

    if distance is None:
        return "no distance"
    if distance > ex.WIDE_MAX_DISTANCE:
        return f"beyond {ex.WIDE_MAX_DISTANCE // 1000} kb"
    return ex._band(distance)


def compiled_distance(rule: Rule, compiled: GeneRow | None) -> int | None:
    """The element-midpoint-to-TSS distance of the compiled side, as lane-notopen measured it."""
    from genomeos.attribution.target_calibration import tss_distance

    if compiled is None:
        return None
    d = tss_distance(rule.midpoint, compiled.tss)
    return None if d is None else int(d)


# ---- what this cannot establish -----------------------------------------------------------------


def limitations() -> dict[str, str]:
    """Carried with every result this module writes, unchanged in force."""
    return {
        "a_count_is_not_a_contribution": (
            "this is a count of symbol resolutions. It says how many rules name a target whose id "
            "differs from the id of the gene the deletion answer measured. It does NOT say the "
            "model's predictions are wrong by that proportion, it is not a measure of any feature's "
            "contribution, and it establishes about no single rule that the rule is wrong"
        ),
        "the_measured_id_is_recovered_not_read": (
            "the deletion answer records a gene NAME and no Ensembl gene id. The measured side is "
            "recovered from the element's recorded window gene list and the window the method "
            "records. Where that recovery does not pick out exactly one gene the rule is counted as "
            "unresolvable with its cause named, never resolved on a preference"
        ),
        "two_annotations_are_not_one": (
            "AlphaGenome scored with its own gene annotation and both sides here are resolved "
            "against GENCODE v50. Where the two releases disagree about a name, the rule lands in "
            "`target_absent_from_the_chromosome_annotation` or in "
            "`no_annotated_locus_of_that_name_overlaps_the_scorer_window`, which are statements about "
            "the annotations and not about the rule"
        ),
        "a_differing_id_is_not_a_repair": (
            "nothing here deletes, rewrites or relabels a compiled rule. If a repair is warranted it "
            "is proposed and left to a ruling"
        ),
        "the_window_is_a_model_of_the_scorer": (
            "the scorer's window is taken from the method it records, as the element midpoint plus "
            "or minus not_open_profile.SCORER_HALF_WINDOW, with gene-body overlap. Every count is "
            f"reported again at a half-window of {SENSITIVITY_HALF_WINDOW} - half AlphaGenome's "
            "2**20 bp interval - so no figure can be read as depending on the rounder number"
        ),
        "correlation_is_not_cause": (
            "where the differing share is reported by activity, by effect band or by distance it is "
            "reported descriptively. No cause is inferred and no axis is claimed to explain another"
        ),
    }


def imported_readings() -> dict[str, Any]:
    """The readings this lane carries from where they were registered, word for word."""
    return {
        "lane_notopen_1678": (
            "1,678 rules carry a TSS further from the element than half the 1 Mb deletion window, so "
            "a repeated gene symbol resolved to another locus - a fault in the distance, not the rule"
        ),
        "lane_notopen_1678_where": "data/results/not_open_profile.json, docs/ATTRIBUTION.md",
        "lane_notopen_1678_population": (
            "the 55,084 rules in state not_open_in_reader, out of the 81,635 assessable rules, out of "
            "the 440,589 rules the compiler emits"
        ),
        "lane_notopen_activity_axis": (
            "0.8511 of the 30,480 assessable repressions against 0.5702 of the 51,103 assessable activations"
        ),
        "lane_notopen_refused_to_claim": (
            "Nothing establishes that any one of the 8,184 rules is wrong; the list says where to "
            "look and what kind of place each one is."
        ),
    }
