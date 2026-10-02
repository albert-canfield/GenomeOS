# SPDX-License-Identifier: AGPL-3.0-or-later
"""What it would cost to give direction rule v2 a second value for the cells it has only one of.

`data/results/respmap_direction_v2.json` (lane-propb, registration 230509a) leaves 112 of the
response map's 166 published predicted assertions UNRESOLVED under direction rule v2, under exactly
two reasons: `one_value_only` on 93 rows and `one_track_seen_twice` on 19. This module prices the
93 and treats the 19 apart. It buys nothing and sends nothing: every figure it reports is read from
a committed result, from the deletion sweep's own cache on disk, and from the copy of the model
client's track metadata already on this machine.

The premise this lane was given, and what the code says instead
---------------------------------------------------------------
The brief framed the 93 as rules "assigned to a cell the sweep never requested". The code refuses
that reading, and `PREMISE_CORRECTED` states the refusal in full. The sweep requests no cell at
all: `predict.alphagenome_adapter.AlphaGenomeAdapter._live_scorer` calls
`client.score_variant(interval=..., variant=..., variant_scorers=[RECOMMENDED_VARIANT_SCORERS
["RNA_SEQ"]])` and passes no `ontology_terms`, no biosample and no cell argument of any kind, so
one response carries the whole RNA-seq track axis. Every one of the 93 cached gene rows records
`n_tracks` 371 and names the assigned cell as the tissue of the retained extreme itself -- measured
by the read `AD_HOC_READ_INCIDENT` describes, kept as that read's output and never re-run, because
the committed program reads no element answer of the population at all.
The cell was in the response and its value is the value v2 reads. What is missing is a
SECOND value for that cell, and the reason is `predict.enhancer_target.aggregate`, which keeps a
per-cell value only `if tissue in cells` for `cells=CELLS` -- four names -- and folds every other
track into `min`/`max` and discards it.

So this is not a coverage gap in what was bought. It is a retention gap in what was kept.

What a request can and cannot buy
---------------------------------
`REQUEST_UNIT` records, against the code, that a request is per ELEMENT and not per cell. It
follows that the price of the 93 is the number of distinct elements among the rows whose cell could
yield a second value at all -- and that number is decided by the track axis, not by money:

* a cell carried by two or more distinct RNA-seq track names can yield two values, so a re-run of
  that element with the cell retained produces the input v2's clause (1) asks for;
* a cell carried by exactly one distinct RNA-seq track name can never yield two values. No request
  count resolves it, because the second value does not exist in any response. `UNPURCHASABLE`
  states that as a property of the track axis.

Nothing here says a purchased second value resolves anything. `WHAT_IT_DOES_NOT_BUY` is explicit:
clause (2) unanimity, amendment 1 distinctness and clause (3) the magnitude floor all still have to
hold of the pair, and none of them is predictable from the axis.

The readings this module carries, and may not strengthen
-------------------------------------------------------
`direction_v2.UNRESOLVED_MEANS`, `direction_v2.FLIPS_ARE_STRUCTURALLY_IMPOSSIBLE` and
`direction_v2.ASSIGNED_CELL` are imported as objects and re-exported unchanged, never retyped. v2
can withhold a direction and can never reverse one; an unresolved class is not absence of
regulation; and the cell on a row is the compiled rule's assigned tissue and not a measured cell.

0 model requests, no network, no money. This module imports no client and holds no key.
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from genomeos.attribution.direction_v2 import (
    ASSIGNED_CELL,
    CACHE_IS_LEGACY,
    FLIPS_ARE_STRUCTURALLY_IMPOSSIBLE,
    UNRESOLVED_MEANS,
)
from genomeos.predict.enhancer_target import CELLS

QUESTION = (
    "For the 93 published predicted assertions direction rule v2 leaves unresolved under "
    "`one_value_only`, and the 19 it leaves unresolved under `one_track_seen_twice`: which of their "
    "assigned cells is carried by an AlphaGenome RNA-seq track at all, how many tracks each is "
    "carried by, how many requests a second value for them would take, and how many of them cannot "
    "be given a second value at any price"
)

POPULATION = (
    "the rows of the committed data/results/respmap_direction_v2.json whose v2_unresolved_reason is "
    "`one_value_only` (93 rows) or `one_track_seen_twice` (19 rows), read from that file's own "
    "assertion rows and from no tally. The other 54 rows of its 166 are resolved and are not priced. "
    "These two counts are lane-propb's, carried unchanged; they are not re-derived here and they are "
    "added to no other count"
)

# ---- the premise this lane was given, refused against the code ------------------------------------

PREMISE_CORRECTED = (
    "the brief described the 93 as rules assigned to `a cell the sweep never requested`, and the "
    "code refuses it. (1) The sweep requests no cell: the only live request path, "
    "predict.alphagenome_adapter.AlphaGenomeAdapter._live_scorer, calls client.score_variant with an "
    "interval, a variant and one variant scorer, and passes no ontology_terms, no biosample and no "
    "cell argument, so the response carries the whole RNA-seq track axis whatever cell is of "
    "interest. (2) Every one of the 93 cached gene rows carries n_tracks 371, the full axis. (3) On "
    "each of the 93 the assigned cell is the tissue of the single retained value itself - "
    "max_drop_tissue on 79 rows and max_rise_tissue on 14 - so that cell's own track was in the "
    "response and was the extreme over all 371. The cell was requested, measured and read. What is "
    "absent is a second value for it, and the cause is retention, not coverage: "
    "predict.enhancer_target.aggregate writes a per-cell value only `if tissue in cells`, with "
    f"cells defaulting to CELLS {list(CELLS)}, and discards every other track's value after folding "
    "it into min and max. The reframing the brief drew from its premise - `we never asked about that "
    "cell` - therefore does not hold, and the figures below are reported as a retention gap"
)

REQUEST_UNIT = (
    "one request buys one element, and cells are not a request dimension. Read from the code: "
    "predict.enhancer_target.score_element returns a cached answer without a request and otherwise "
    "calls `scorer(chrom, start, seq, seq[0])` exactly once, then `aggregate(effects)`; the scorer "
    "is one client.score_variant call carrying no cell argument; and `aggregate` is the step that "
    "selects cells, after the response has arrived. So the request count for a set of rows is the "
    "number of DISTINCT ELEMENTS those rows name, and adding a cell to the retained set costs 0 "
    "requests by itself. It is a code change to aggregate's `cells`, not a purchase"
)

PURCHASABLE_RULE = (
    "a row of the 93 is PURCHASABLE if and only if two or more distinct RNA-seq track names in the "
    "model client's own track metadata carry its assigned cell's label under TRACK_LABEL_RULE. Two "
    "names mean two columns of the response, so re-running that element with the cell retained "
    "yields two values and satisfies v2's clause (1). One name means one column: the response has "
    "one value for that cell and will have one value on every future response too, so no request "
    "count produces a second. The rule is about the axis and never about the values, which are not "
    "read here"
)

UNPURCHASABLE = (
    "a cell carried by exactly one RNA-seq track name is unpurchasable at any price for v2's "
    "purposes. This is not a statement that the cell is unmeasured, that its track is poor, or that "
    "its rule is wrong: the cell is measured, by one track, and v2's clause (1) asks for two values "
    "from one cell. A second value does not exist to be bought. The only things that would change "
    "the count are a change to the model's track axis, which this project does not control, or a "
    "change to v2's clause (1), which is a rule change and not a purchase"
)

WHAT_IT_DOES_NOT_BUY = (
    "a purchased second value is an input to the sign rule and not a resolution. Having two "
    "numerically present values, v2 still requires: they are numerically DISTINCT (amendment 1, "
    "which refuses two identical values whatever their tracks), their signs are unanimous with no "
    "zero (clause 2), and the largest |value| clears MAGNITUDE_FLOOR (clause 3). None of the three "
    "is predictable from the track axis, so the request count below is a count of rows that would "
    "GET a second value and is not a count of rows that would resolve. It also buys nothing for the "
    "unpurchasable rows, nothing for the 19, nothing about any cell's biology, and no direction for "
    "any rule: a re-run can only move a row from unresolved to resolved or leave it unresolved, by "
    "FLIPS_ARE_STRUCTURALLY_IMPOSSIBLE"
)

HOW_ESTABLISHED = (
    "every per-cell row says how its track count was established, and the two are not the same "
    "claim. `client_track_metadata` means the cell's label is equal, character for character, to "
    "the label a track of the client's own metadata carries under TRACK_LABEL_RULE, and the count is "
    "the number of distinct track names that do so. `name_only` would mean a label matched by "
    "spelling, inflection or synonym rather than by equality, and `absent` that no track carries it. "
    "No row is resolved by judgement, by an ontology lookup or by a synonym table: a label either "
    "equals a track's own label or it does not, and whichever it is, is recorded"
)

TRACK_LABEL_RULE = (
    "a track's label is predict.alphagenome_adapter.tissue_names' own rule, applied to the metadata "
    "rows rather than to a response: the gtex_tissue when the row carries one that is not empty and "
    "not the string `nan`, else the biosample_name. This is the same rule that decides the `tissue` "
    "a value arrives under in the effects, so a cell label on a compiled rule is a track label by "
    "construction. tests/test_cellcover.py holds this function against the adapter's own "
    "tissue_names on a stub response so the two cannot drift"
)

AXIS_IDENTIFIED = (
    "the metadata copy's RNA-seq rows are identified with the response's track axis by an equality "
    "neither side was fitted to: the copy holds 667 RNA-seq rows over exactly 371 DISTINCT track "
    "names (the surplus is the strand column, 271 `+` with 271 `-` and 125 `.`), and a real "
    "response's own record on disk gives its axis width as 371. 371 == 371, so one response column "
    "is one distinct RNA-seq track name, with the two strands of a stranded track scored as one "
    "column. The program asserts the equality and refuses rather than reporting a per-cell count it "
    "cannot ground. The response side is read from ONE loose per-element answer "
    "(AXIS_WITNESS), whose `model.tracks` the adapter sets to the length of the response's own track "
    "axis and beside which it records `model.tracks_sha256`, a digest of that axis's names"
)

AXIS_WITNESS_IS_NOT_THE_POPULATION = (
    "stated so the identification is not read as stronger than it is. The loose answer that "
    "witnesses the axis width is not one of the 93: the 93's elements exist only inside the "
    "compressed per-chromosome archives, and reading one of those is the heaviest read in this "
    "project (see AD_HOC_READ_INCIDENT), so the program does not do it. That every cached gene row "
    "of the 93 records n_tracks 371 was MEASURED, by an ad hoc read this project has since "
    "forbidden, and is reported under `measured_by_a_read_now_forbidden` rather than recomputed. The "
    "witness answer was written by the same adapter against the same scorer, which is why its axis "
    "width is evidence about the axis; that it was written on the same model revision as the 93 is "
    "NOT established, and if the axis changed between them the per-cell counts would be mispriced "
    "and nothing here would detect it"
)

AD_HOC_READ_INCIDENT = (
    "recorded against this lane. While exploring, this lane ran two inline `python3 -c` loops that "
    "called predict.enhancer_target.load_cached once per row over the 93 and then the 19. "
    "`load_cached` falls back to the chromosome archive, so each loop decompressed whole archives: "
    "one reached 23 GB RSS, took swap to 13.1 of 14.3 GB, grew the swap file to 14 GB and took free "
    "disk to 7.8 GiB, under the project's 10 GB floor. It had exited on its own before a peer "
    "looked. The coordinator has since made the rule general: no inline exploration over the data "
    "stores, because a one-liner has no gate sample, no RSS ceiling and no name on the board, so "
    "nobody can tell what is paging the machine or who owns it. The figures that read produced are "
    "kept, named as its output, and are not re-run; the committed program reads no element answer of "
    "the population at all, which is why it needs the witness above"
)

METADATA_COPY_LIMITATION = (
    "the track metadata is read from the copy a peer session saved to disk "
    "(data/cache/entex/alphagenome_track_metadata_copy.csv, the same authority "
    "attribution.context_evidence uses for its label-to-term mapping). Re-fetching it is a request "
    "this lane may not make. That the copy is the client's CURRENT output_metadata, and not a "
    "superseded one, is therefore NOT ESTABLISHED here; the 371 == 371 agreement with the sweep's "
    "own cached n_tracks is the only check this lane can run on it, and it is reported as that. A "
    "cell whose track count changed after the copy was saved would be mispriced, and no figure here "
    "detects that"
)

THE_NINETEEN = (
    "the 19 `one_track_seen_twice` rows are not priced with the 93 and their count is not added to "
    "it. They have two retained values that are numerically identical, and amendment 1 refuses them "
    "because the cache cannot say whether that is one track counted twice or two tracks agreeing. "
    "The question is whether track identity can be recorded for a CACHED answer. It cannot: the "
    "response is not on disk. What the cache keeps per gene row is the extremes with their tissue "
    "labels, `by_cell` (the last emitted value per cell name), `n_tracks` and nothing else; the "
    "effects the writer saw are (gene, tissue, value) triples carrying no column, so the column was "
    "already gone before `aggregate` ran, and enhancer_target._cell_summary keeps only a sorted "
    "MULTISET per cell. For a cached answer the answer is therefore `cannot be recorded, at any "
    "price` - a finding, not a gap. For a FRESH response it is already recorded and needs no new "
    "code: alphagenome_adapter.recorded_axis keeps each cell's values in COLUMN ORDER, one entry per "
    "column of that cell, for its `cells` default CELL_TRACKS; every cell the 19 carry is in that "
    "default, so a new response for those elements would show two entries and settle whether one "
    "track was counted twice or two agreed. That is positional identity within one response and not "
    "a track name, which is weaker than a name but is exactly what amendment 1 asks. It would still "
    "not resolve the 19, because amendment 1 as implemented refuses two numerically identical values "
    "whatever their columns: consulting the identity is a change to the rule, which is not a "
    "purchase and is not proposed here. The element count is reported so the request side is priced, "
    "with no claim that it resolves anything"
)

NO_RECOMMENDATION = (
    "this module makes no recommendation about whether to spend. It reports a count, what the count "
    "covers and what it does not buy. Nothing here argues for the purchase, no figure is framed as a "
    "return, and the decision is not this lane's"
)

CELL2_GROUP_NOT_USED = (
    "no cell is grouped with another. Each of the 59 assigned cells is counted on its own label, so "
    "attribution.cell2's operational grouping - which is an operational grouping and NOT established "
    "biological independence - is neither used nor needed, and no T-cell subset is merged with "
    "another"
)

CANNOT_ESTABLISH = (
    "whether a purchased second value would resolve any row: that needs the values, which need the "
    "requests this lane may not make",
    "whether the saved track metadata copy is the client's current output_metadata: re-fetching it "
    "is a request (METADATA_COPY_LIMITATION)",
    "the money price of one request: no committed artefact of this project states a per-request "
    "price, so the cost is reported in requests and in nothing else",
    "whether the two tracks of a two-track cell are biological replicates: in every case found they "
    "are one biosample under two assay titles, which _cell_summary already says in those words is "
    "not known to be replication",
    "whether a cell absent from the metadata is absent from the model: no cell of this population "
    "was absent, so the branch was never taken and its reading is untested",
)

READINGS_CARRIED = {
    "unresolved_means": UNRESOLVED_MEANS,
    "flips_are_structurally_impossible": FLIPS_ARE_STRUCTURALLY_IMPOSSIBLE,
    "assigned_cell": ASSIGNED_CELL,
    "cache_is_legacy": CACHE_IS_LEGACY,
}

# ---- the inputs ------------------------------------------------------------------------------------

#: The copy of the model client's track metadata already on this machine. The same file
#: `attribution.context_evidence.TRACK_METADATA` names, read for its RNA-seq rows only.
TRACK_METADATA = Path("data/cache/entex/alphagenome_track_metadata_copy.csv")

#: The committed result whose rows are the population.
V2_RESULT = Path("data/results/respmap_direction_v2.json")

#: The output the RNA-seq gene scorer reads, as the metadata names it.
RNA_OUTPUT = "rna_seq"

#: The two reasons priced, apart. Both are members of `direction_v2.UNRESOLVED_REASONS`.
ONE_VALUE_ONLY = "one_value_only"
ONE_TRACK_SEEN_TWICE = "one_track_seen_twice"

#: v2's clause (1): a direction may not rest on one track.
VALUES_NEEDED = 2

#: The track-axis width, asserted rather than assumed (AXIS_IDENTIFIED).
EXPECTED_AXIS = 371

#: One loose per-element answer whose `model.tracks` is the response's own axis width. A single small
#: JSON file read by path: never `load_cached`, which falls back to a compressed chromosome archive.
AXIS_WITNESS = Path("data/knowledge/alphagenome/elements_hct116/chr12/EH38E1630471.json")

#: The population's elements are archive-only, measured by the read AD_HOC_READ_INCIDENT describes.
MEASURED_BY_A_READ_NOW_FORBIDDEN = {
    "what": "each cached gene row's n_tracks, and which extreme field names the assigned cell, over "
    "the 93 and the 19",
    "one_value_only_n_tracks": {"371": 93},
    "one_value_only_assigned_cell_is": {"max_drop_tissue": 79, "max_rise_tissue": 14},
    "one_track_seen_twice_n_tracks": {"371": 19},
    "one_track_seen_twice_cells": {"K562": 17, "GM12878": 2},
    "by_cell_entries_per_row": 4,
    "schema": "no cached answer of either set carries `model`, `by_cell_summary` or a gene-axis "
    "record: the row keys are gene, n_tracks, mean_log2fc, max_drop_log2fc, max_drop_tissue, "
    "max_rise_log2fc, max_rise_tissue and by_cell, and nothing else",
    "how": "two inline loops over load_cached, which paged the machine (AD_HOC_READ_INCIDENT). The "
    "committed program does not repeat them and cannot reproduce these five figures; they are "
    "reported as that read's output and carry its authority and no more",
    "reproducible": False,
}


def track_label(row: dict[str, Any]) -> str:
    """One metadata row's label, by `alphagenome_adapter.tissue_names`' rule (TRACK_LABEL_RULE)."""
    gtex = row.get("gtex_tissue")
    if gtex and str(gtex) not in ("nan", ""):
        return str(gtex)
    return str(row.get("biosample_name") or "")


@dataclass(frozen=True)
class Axis:
    """The RNA-seq track axis of the client's metadata, one entry per distinct track name."""

    tracks: dict[str, dict[str, Any]]  # track name -> its metadata row
    rows_read: int
    source: str

    @property
    def width(self) -> int:
        return len(self.tracks)

    def by_label(self) -> dict[str, list[dict[str, Any]]]:
        out: dict[str, list[dict[str, Any]]] = {}
        for r in self.tracks.values():
            out.setdefault(track_label(r), []).append(r)
        return out


def rna_axis(path: Path = TRACK_METADATA) -> Axis:
    """The RNA-seq track axis, read from the saved metadata copy and never from the model."""
    if not path.exists():
        raise FileNotFoundError(
            f"the model client's track metadata is not on this machine: {path}. The per-cell track "
            "counts are read from it and are never guessed from a cell's name."
        )
    tracks: dict[str, dict[str, Any]] = {}
    rows_read = 0
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("output") != RNA_OUTPUT:
                continue
            rows_read += 1
            tracks.setdefault(str(row.get("name") or ""), row)
    return Axis(tracks=tracks, rows_read=rows_read, source=str(path))


@dataclass(frozen=True)
class CellCoverage:
    """One assigned cell: how many RNA-seq tracks carry it, and how that was established."""

    cell: str
    rows: int
    rna_tracks: int
    established_how: str
    ontology_terms: tuple[str, ...] = ()
    assay_titles: tuple[str, ...] = ()
    data_sources: tuple[str, ...] = ()
    biosample_types: tuple[str, ...] = ()
    label_via: str = ""

    @property
    def purchasable(self) -> bool:
        """Whether a second value for this cell exists in any response (PURCHASABLE_RULE)."""
        return self.rna_tracks >= VALUES_NEEDED

    def to_dict(self) -> dict[str, Any]:
        return {
            "cell": self.cell,
            "rows": self.rows,
            "rna_seq_tracks": self.rna_tracks,
            "track_exists": self.rna_tracks > 0,
            "established_how": self.established_how,
            "second_value_exists_in_any_response": self.purchasable,
            "ontology_terms": list(self.ontology_terms),
            "assay_titles": list(self.assay_titles),
            "data_sources": list(self.data_sources),
            "biosample_types": list(self.biosample_types),
            "label_via": self.label_via,
        }


def coverage(cell_rows: Counter[str], axis: Axis) -> list[CellCoverage]:
    """One `CellCoverage` per assigned cell, ordered by row count then name."""
    by_label = axis.by_label()
    out: list[CellCoverage] = []
    for cell, n in cell_rows.items():
        hits = by_label.get(cell, [])
        out.append(
            CellCoverage(
                cell=cell,
                rows=n,
                rna_tracks=len(hits),
                established_how="client_track_metadata" if hits else "absent",
                ontology_terms=tuple(sorted({str(r.get("ontology_curie") or "") for r in hits})),
                assay_titles=tuple(sorted(str(r.get("Assay title") or "") for r in hits)),
                data_sources=tuple(sorted({str(r.get("data_source") or "") for r in hits})),
                biosample_types=tuple(sorted({str(r.get("biosample_type") or "") for r in hits})),
                label_via=(
                    "gtex_tissue"
                    if any(r.get("gtex_tissue") not in (None, "", "nan") for r in hits)
                    else ("biosample_name" if hits else "")
                ),
            )
        )
    out.sort(key=lambda c: (-c.rows, c.cell))
    return out


@dataclass
class Price:
    """The request count for one set of rows, with what it covers."""

    rows: int
    elements: tuple[str, ...]
    cells: tuple[str, ...]

    @property
    def requests(self) -> int:
        """One request per distinct element (REQUEST_UNIT)."""
        return len(self.elements)

    def to_dict(self) -> dict[str, Any]:
        return {
            "rows": self.rows,
            "requests": self.requests,
            "request_unit": "one request per distinct element",
            "elements": list(self.elements),
            "cells_covered": list(self.cells),
        }


def price(rows: list[dict[str, Any]]) -> Price:
    """The request count for `rows`, counted over distinct elements and never over rows or cells."""
    return Price(
        rows=len(rows),
        elements=tuple(sorted({str(r["element"]) for r in rows})),
        cells=tuple(sorted({str(r["cell"]) for r in rows})),
    )


def split(rows: list[dict[str, Any]], cover: list[CellCoverage]) -> dict[str, list[dict[str, Any]]]:
    """`rows` split by whether a second value for the assigned cell exists at all."""
    buyable = {c.cell for c in cover if c.purchasable}
    out: dict[str, list[dict[str, Any]]] = {"purchasable": [], "unpurchasable_at_any_price": []}
    for r in rows:
        key = "purchasable" if str(r["cell"]) in buyable else "unpurchasable_at_any_price"
        out[key].append(r)
    return out


# ---- the counts, named before the run --------------------------------------------------------------

COUNTS_NAMED = (
    "population_rows: the rows read from the committed v2 result under each of the two reasons, "
    "separately and never pooled",
    "distinct_assigned_cells: how many distinct cell labels those rows carry",
    "axis_width: the distinct RNA-seq track names in the saved metadata copy beside one real "
    "response record's own axis width, reported together so the identification is visible; the "
    "population's own n_tracks is carried under measured_by_a_read_now_forbidden and not recomputed",
    "cells_with_no_track: assigned cells no RNA-seq track carries, which are unpurchasable because "
    "the cell is absent and not because it is thin",
    "cells_with_one_track: assigned cells exactly one RNA-seq track carries, the unpurchasable count",
    "cells_with_two_or_more_tracks: the purchasable cells",
    "rows_unpurchasable_at_any_price and rows_purchasable: the same split counted over rows",
    "requests: distinct elements among the purchasable rows, one request each",
    "nineteen: the same figures for the 19, reported apart, with whether track identity can be "
    "recorded for a cached answer",
    "per_cell: one row per assigned cell with its track count and how it was established",
)


def axis_witness(path: Path = AXIS_WITNESS) -> dict[str, Any]:
    """One response's own axis width, from a single loose per-element answer read by path.

    Deliberately not `enhancer_target.load_cached`: that falls back to the compressed chromosome
    archive, and the population's elements are archive-only (AD_HOC_READ_INCIDENT). This opens one
    small JSON file, reads two fields of its `model` record and keeps nothing else.
    """
    import json

    if not path.exists():
        raise FileNotFoundError(
            f"the axis witness is not on this machine: {path}. AXIS_IDENTIFIED rests on a real "
            "response's own axis width, and the per-cell counts are not reported without one."
        )
    if path.suffix != ".json":
        raise ValueError(f"the axis witness must be a loose per-element answer, not {path}")
    answer = json.loads(path.read_text())
    model = answer.get("model") or {}
    return {
        "path": str(path),
        "tracks": int(model.get("tracks") or answer.get("tracks") or 0),
        "tracks_sha256": str(model.get("tracks_sha256") or ""),
        "scorer": str(model.get("scorer") or ""),
        "model_version": str(model.get("model_version") or ""),
        "note": AXIS_WITNESS_IS_NOT_THE_POPULATION,
    }


@dataclass
class Report:
    """Everything the run reports, assembled from committed rows and from the metadata copy."""

    reason: str
    rows: int
    cells: list[CellCoverage]
    purchasable: Price
    unpurchasable: Price
    axis_width: int
    witness_width: int = 0

    def counts(self) -> dict[str, Any]:
        no_track = [c for c in self.cells if c.rna_tracks == 0]
        one = [c for c in self.cells if c.rna_tracks == 1]
        many = [c for c in self.cells if c.rna_tracks >= VALUES_NEEDED]
        return {
            "population_rows": self.rows,
            "distinct_assigned_cells": len(self.cells),
            "axis_width_from_metadata": self.axis_width,
            "axis_width_from_a_response_record": self.witness_width,
            "cells_with_no_track": len(no_track),
            "cells_with_one_track": len(one),
            "cells_with_two_or_more_tracks": len(many),
            "rows_on_cells_with_no_track": sum(c.rows for c in no_track),
            "rows_unpurchasable_at_any_price": self.unpurchasable.rows,
            "rows_purchasable": self.purchasable.rows,
            "requests_to_give_every_purchasable_row_a_second_value": self.purchasable.requests,
            "requests_that_would_buy_anything_for_the_unpurchasable_rows": 0,
        }
