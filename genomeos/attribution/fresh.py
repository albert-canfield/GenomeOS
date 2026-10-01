# SPDX-License-Identifier: AGPL-3.0-or-later
"""Is there a CRISPRi enhancer-gene set this project has never read, outside K562, big enough to test?

Every CRISPRi comparison this project holds is K562, or a cell type it has already scored. The
second-cell-type lane closed eligible by locus count but with its power not computable, and its
never-scored stratum held 10 loci and 16 positives. So the prior question is whether a genuinely
unread set exists at all.

This module decides that on labels and coordinates only. Three floors are fixed here, before any
candidate is opened, and none may be moved after a count:

* independent loci at least `cell2.POOLED_LOCUS_FLOOR`, under `cell2.INDEPENDENT_LOCUS_RULE`;
* measured positives at least `POSITIVE_FLOOR`;
* the cell context is not K562, and no outcome from the set has been read, scored or quoted here.

The grouping is `cell2`'s, imported rather than restated, so a locus count here means exactly what it
means there: an operational grouping for deciding whether a stratum has enough distinct places in the
genome to be worth measuring, *not* established biological independence. Nothing downstream may cite
a locus count from this module as an independence guarantee.

The blinding is enforced, not promised. `PERMITTED_COLUMNS` is the whole set of fields a candidate
table may be read for, and `permit` raises on anything else, so an effect size, a fold change, a
p-value, an FDR or a model score cannot enter this module even by accident. Counting positives from a
published hit-call column is the point of the lane; reading the magnitude of an effect is forbidden.

Freshness is a property of the *record*, not an assumption: a set is fresh only if no outcome from it
has been scored, read or quoted here, and `Candidate.read_here` carries how that was established.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

from genomeos.attribution import cell2

#: Measured positives below this is a no-go. Set before any candidate table was opened, alongside the
#: locus floor this lane takes from `cell2`; recorded here so that no later run can move it after a
#: count. 30 is the smallest positive count at which the registered locus convention can distribute
#: positives over the 20-locus floor without a locus carrying the result on its own.
POSITIVE_FLOOR = 30

#: The locus floor, taken from the second-cell-type lane rather than restated, so that a count here
#: and a count there mean the same thing and can be compared.
LOCUS_FLOOR = cell2.POOLED_LOCUS_FLOOR

#: The cell context the project's CRISPRi result already covers, so not a fresh context.
EXCLUDED_CELL = cell2.PRIMARY_CELL

#: Every field a candidate's table may be read for. A column outside this set is refused by `permit`.
#: `Hit` is a published hit call, a label; the lane counts labels and never reads a magnitude.
PERMITTED_COLUMNS = frozenset(
    {
        "Pair",
        "Enhancer",
        "Enh",
        "GeneSymbol",
        "Gene",
        "GeneID",
        "Hit",
        "Tested",
        "Wellpowered",
        "WellPowered_at_FC_0.15",
        "WellPowered_at_FC_0.25",
        "EnhancerCoord",
        "Coord",
        "LinkedGene",
        "GeneTSS",
    }
)

#: Column names that name an outcome's magnitude. Listed so the refusal can say why, and so a reader
#: of this module can see that the lane knew what it was declining to read.
FORBIDDEN_SUBSTRINGS = (
    "log2fc",
    "logfoldchange",
    "foldchange",
    "effectsize",
    "effect_size",
    "zscore",
    "z_",
    "_p",
    "p_",
    "pvalue",
    "fdr",
    "score",
    "auprc",
    "expression",
    "sensitivity",
)

_COORD = re.compile(r"^(chr[0-9A-Za-z_]+):([0-9]+)-([0-9]+)$")


class BlindingError(ValueError):
    """Raised when a column outside `PERMITTED_COLUMNS` is asked for."""


def permit(columns: list[str]) -> list[str]:
    """`columns`, unchanged, if every one of them is permitted; otherwise refuse and say which.

    This is the lane's blinding condition in code. It runs before a candidate workbook is opened, so
    the forbidden column is never loaded, let alone read.
    """
    bad = [c for c in columns if c not in PERMITTED_COLUMNS]
    if bad:
        lowered = {c: [s for s in FORBIDDEN_SUBSTRINGS if s in c.lower()] for c in bad}
        named = ", ".join(
            f"{c} (names an outcome magnitude: {'; '.join(lowered[c])})" if lowered[c] else c for c in bad
        )
        raise BlindingError(
            f"this lane reads labels, coordinates and power flags only; refused: {named}. "
            "Counting positives from a published hit-call column is permitted; reading the "
            "magnitude of an effect, its p-value, its FDR or a model score is not"
        )
    return list(columns)


def parse_coord(text: str) -> tuple[str, int, int]:
    """`chr1:112390877-112391209` as (chrom, start, end). Raises on anything it cannot spell out,
    because a silently dropped element would change a count without changing the count's claim."""
    match = _COORD.match(str(text).strip())
    if not match:
        raise ValueError(f"not a chrom:start-end interval: {text!r}")
    start, end = int(match.group(2)), int(match.group(3))
    if end < start:
        raise ValueError(f"interval ends before it starts: {text!r}")
    return match.group(1), start, end


class _Element:
    """The blinded view of one measured pair, shaped for `cell2.locus_key`."""

    __slots__ = ("cell", "chrom", "start", "end", "gene")

    def __init__(self, cell: str, chrom: str, start: int, end: int, gene: str) -> None:
        self.cell, self.chrom, self.start, self.end, self.gene = cell, chrom, start, end, gene


def locus_keys(cell: str, positives: list[dict[str, Any]], coord_field: str, gene_field: str):
    """The `cell2` locus keys of a candidate's positives. Every row must parse; none is skipped."""
    keys = []
    for row in positives:
        chrom, start, end = parse_coord(row[coord_field])
        keys.append(cell2.locus_key(_Element(cell, chrom, start, end, str(row[gene_field]))))
    return keys


@dataclass
class Candidate:
    """One candidate set, and what was established about it rather than assumed."""

    name: str
    cell_context: str
    assay: str
    accession: str
    pairs_tested: int | None
    positives: int | None
    independent_loci: int | None
    read_here: str
    access: str
    notes: str = ""
    distinct_genes: int | None = None
    chromosomes: int | None = None

    @property
    def is_k562(self) -> bool:
        return EXCLUDED_CELL.lower() in self.cell_context.lower()

    #: The only phrase that asserts freshness. A bare "no", or a "not established", must not pass:
    #: the claim is that *no outcome* from the set has been read, and it has to be said in those
    #: words, with how it was established after it.
    FRESH_PREFIX = "no outcome"

    @property
    def fresh(self) -> bool:
        """Fresh only on a positive statement that no outcome has been read; absence of a claim is
        not freshness, so an unestablished record reads as not fresh."""
        return self.read_here.strip().lower().startswith(self.FRESH_PREFIX)

    def misses(self) -> list[str]:
        """Every floor this candidate fails, with the shortfall. An unknown count is a miss, because a
        floor that cannot be checked has not been cleared."""
        out: list[str] = []
        if self.is_k562:
            out.append(f"cell context is {EXCLUDED_CELL}, which the result already covers")
        if not self.fresh:
            out.append(f"not fresh: {self.read_here}")
        if self.independent_loci is None:
            out.append(f"independent loci not established, so the floor of {LOCUS_FLOOR} is unchecked")
        elif self.independent_loci < LOCUS_FLOOR:
            out.append(
                f"independent loci {self.independent_loci}, short of {LOCUS_FLOOR} "
                f"by {LOCUS_FLOOR - self.independent_loci}"
            )
        if self.positives is None:
            out.append(f"positives not established, so the floor of {POSITIVE_FLOOR} is unchecked")
        elif self.positives < POSITIVE_FLOOR:
            out.append(
                f"measured positives {self.positives}, short of {POSITIVE_FLOOR} "
                f"by {POSITIVE_FLOOR - self.positives}"
            )
        return out

    @property
    def eligible(self) -> bool:
        return not self.misses()

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "cell_context": self.cell_context,
            "assay": self.assay,
            "accession": self.accession,
            "pairs_tested": self.pairs_tested,
            "measured_positives": self.positives,
            "independent_loci": self.independent_loci,
            "distinct_measured_genes": self.distinct_genes,
            "chromosomes": self.chromosomes,
            "read_here": self.read_here,
            "fresh": self.fresh,
            "access": self.access,
            "misses": self.misses(),
            "eligible": self.eligible,
            "notes": self.notes,
        }


@dataclass
class Verdict:
    """The lane's answer over all candidates. A no-go is a finished result, not a failure."""

    candidates: list[Candidate] = field(default_factory=list)

    @property
    def eligible(self) -> list[Candidate]:
        return [c for c in self.candidates if c.eligible]

    def to_dict(self) -> dict:
        passing = self.eligible
        return {
            "floors": {
                "independent_loci_at_least": LOCUS_FLOOR,
                "measured_positives_at_least": POSITIVE_FLOOR,
                "cell_context_not": EXCLUDED_CELL,
                "never_read_here": True,
                "locus_rule": cell2.INDEPENDENT_LOCUS_RULE,
                "locus_span_bp": cell2.INDEPENDENT_LOCUS_SPAN,
                "fixed_before_any_candidate_was_opened": True,
                "note": (
                    "the locus rule and span are imported from genomeos/attribution/cell2.py, not "
                    "restated, so a locus count here means what it means there. It is an operational "
                    "grouping, not established biological independence"
                ),
            },
            "candidates_examined": len(self.candidates),
            "candidates": [c.to_dict() for c in self.candidates],
            "eligible_count": len(passing),
            "verdict": "go" if passing else "no-go",
            "eligible_sets": [c.name for c in passing],
        }
