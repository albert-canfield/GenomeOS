# SPDX-License-Identifier: AGPL-3.0-or-later
"""The attributed elements of a chromosome, from whichever deletion runs exist.

Three runs name a target gene per element: the whole-chromosome scoring
(`enhancer_targets_all_<chrom>`, every element in a node; its committed result is
a summary whose `elements_where` points at the local table), the constrained
sample and the uniform sample. The whole-chromosome run is read first; the samples
add only elements it does not hold. The closure, the organiser and the compiler
all read through here, so none of them silently falls back to a sample when the
whole table exists.

Every table those runs write is compact: one target gene per element. That is a
property of the projection, not of the sweep, and `ElementResponses` below reads
the sweep's own per-element response cache instead, where every gene in the
scorer's window is kept with its signed change. A caller that asks the compact
table about a gene it does not name gets nothing back and has, until 2026-09-22,
written a zero — see `crispri.ElementCache`, the first consumer moved off it.
"""

from __future__ import annotations

import gzip
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from genomeos.results import RESULTS_DIR, load_result

RUNS = ("enhancer_targets_all", "constrained_targets", "enhancer_targets")
ORIGIN = {"enhancer_targets_all": "all", "constrained_targets": "constrained", "enhancer_targets": "uniform"}
ELEMENT_CACHE = Path("data/knowledge/alphagenome/elements")

SCORED = "scored"
NOT_IN_WINDOW = "the gene is not in the scorer's window at this element"
NOT_ON_TRACK = "the gene was scored, but not on this cell's own track"
NOT_CACHED = "this element is not in the response cache"


def run_elements(name: str, chrom: str, results_dir: Path = RESULTS_DIR) -> list[dict[str, Any]]:
    """One run's elements: inline, or from the local table a summary points at. A summary whose
    table is missing raises, rather than reading as a run with no elements."""
    r = load_result(f"{name}_{chrom}", results_dir) or {}
    if "elements" in r:
        return r["elements"]
    where = r.get("elements_where")
    if not where:
        return []
    p = Path(where)
    if not p.is_absolute() and not p.exists():
        p = results_dir.parent.parent / where  # a relative path is relative to the project root
    if not p.exists():
        raise FileNotFoundError(
            f"{name}_{chrom} points at {where}, which is not on this machine; rerun the job"
        )
    t = json.loads(p.read_text())
    return t if isinstance(t, list) else t.get("elements", [])


def attributed(chrom: str, results_dir: Path = RESULTS_DIR, coding: bool = True) -> list[dict[str, Any]]:
    """Every element with a named (coding) target across the runs, deduplicated by id, each with the
    run it came from under `origin`."""
    seen: set[str] = set()
    out: list[dict[str, Any]] = []
    for name in RUNS:
        for e in run_elements(name, chrom, results_dir):
            pc = e.get("predicted_coding" if coding else "predicted") or {}
            if not pc.get("gene") or e["id"] in seen:
                continue
            seen.add(e["id"])
            out.append({**e, "origin": ORIGIN[name]})
    return out


def runs_present(chrom: str, results_dir: Path = RESULTS_DIR) -> dict[str, int]:
    return {name: len(run_elements(name, chrom, results_dir)) for name in RUNS}


def genes_at_bar(record: dict[str, Any] | None, min_effect: float = 0.1) -> list[tuple[str, float]] | None:
    """Every gene of one cached element that reaches `min_effect`, strongest first, as (gene, signed).

    The size is `predict.enhancer_target.predict_target`'s own: the larger of the predicted fall and
    the predicted rise over every track, ties to the fall. So the head of this list is the element's
    compact `predicted` gene whenever it has one, and the list is empty exactly when the compact table
    records none -- the window changes which and how many genes move, never whether one does. None
    means the element was not cached, which is not the same as an empty list.
    """
    if record is None:
        return None
    out = []
    for g in record.get("genes") or []:
        drop, rise = g.get("max_drop_log2fc"), g.get("max_rise_log2fc")
        if not g.get("gene") or drop is None:
            continue
        fall = -float(drop)
        up = float(rise) if rise is not None else float("-inf")
        size, signed = (fall, float(drop)) if fall >= up else (up, up)
        if size >= min_effect:
            out.append((g["gene"], round(signed, 4), size))
    out.sort(key=lambda t: -t[2])
    return [(name, signed) for name, signed, _ in out]


@dataclass(frozen=True)
class Response:
    """What the sweep predicted for one gene at one element, or why it did not say.

    `value` is a signed log2 fold change on deleting the element: negative is a predicted fall.
    It is None whenever the sweep has no number, and `reason` then says which silence it was —
    the gene was outside the scorer's window, the gene was scored but not on that cell's own
    track, or the element was never cached. A caller must not read None as zero; the whole point
    of this type is that "no prediction" and "a predicted change of zero" are different answers,
    and the compact tables could not tell them apart.
    """

    value: float | None
    reason: str

    @property
    def answered(self) -> bool:
        return self.value is not None

    def __bool__(self) -> bool:  # `if response:` means "the sweep answered", never "non-zero"
        return self.value is not None


class ElementResponses:
    """The sweep's per-element response cache: every gene it scored at an element, not one.

    The compact tables this module loads (`run_elements`, `attributed`) keep a single target gene
    per element, so a question about any other gene has no answer in them and has historically
    been recorded as a zero. The same sweep wrote this cache with `threshold=0.0`
    (`scripts/enhancer_targets_all.py`, `worker_scorer`): every gene in the scorer's 1 Mb window,
    with a signed log2 fold change on each of K562, HepG2, GM12878 and IMR-90's own track and the
    strongest fall across all 371 tracks. Where both sources carry a value they agree exactly
    (0 disagreements over 24 chromosomes, 2026-09-22), because the cache is the same run with the
    censoring removed rather than a later re-score.

    What it costs, measured 2026-09-22. The tree is 775 MB gzipped over 24 chromosomes and must
    never be read whole. One chromosome's archive is held at a time and switching drops the
    previous one, so the bill is the largest chromosome, not the tree: chr21 (9.8 MB gzipped,
    12,158 elements) takes 0.6 s and about 0.6 GB resident, chr1 (77 MB, 88,302 elements) takes
    5.7 s and peaks near 2.9 GB. That is heavy enough to matter — a caller walks its elements
    chromosome by chromosome, as `crispri.annotate` does, one reader per process, and a caller
    that hops between chromosomes pays a full decompression every time (`_loaded` records them).
    Elements the archive does not hold are looked up as loose per-element files under
    `<root>/<chrom>/<id>.json`, which cost nothing to keep.

    `crispri.ElementCache` is the same reader, written first for the CRISPRi benchmark; this is
    the general one, and its `value()` returns what that one returns for the same arguments.
    """

    def __init__(self, root: Path = ELEMENT_CACHE) -> None:
        self.root = Path(root)
        self._chrom: str | None = None
        self._archive: dict[str, Any] = {}
        self._loaded: list[str] = []  # the chromosomes decompressed, in order, for cost reporting

    def _load(self, chrom: str) -> None:
        if self._chrom == chrom:
            return
        p = self.root / f"{chrom}.json.gz"
        archive: dict[str, Any] = {}
        if p.exists():
            with gzip.open(p, "rt") as fh:
                archive = json.load(fh)
        self._chrom, self._archive = chrom, archive
        self._loaded.append(chrom)

    def element(self, chrom: str, element_id: str) -> dict[str, Any] | None:
        """One element's cached response, or None when the sweep never scored that element."""
        self._load(chrom)
        hit = self._archive.get(element_id)
        if hit is None:
            p = self.root / chrom / f"{element_id}.json"
            if p.exists():
                hit = json.loads(p.read_text())
        return hit

    def genes(self, chrom: str, element_id: str) -> dict[str, dict[str, Any]] | None:
        """Every gene in the scorer's window at this element, by symbol; None when not cached.

        An empty dict would mean "scored, and no gene in the window", which is not the same as
        None, "never scored" — so the two are kept apart here as well.
        """
        hit = self.element(chrom, element_id)
        if hit is None:
            return None
        return {g["gene"]: g for g in hit.get("genes") or [] if g.get("gene")}

    def response(self, chrom: str, element_id: str, gene: str, cell: str | None = None) -> Response:
        """What the sweep predicted for this gene at this element, in this cell, or why it did not.

        `cell` names one of the four lines the sweep scored on its own track. With `cell=None` the
        value is the gene's strongest predicted fall across every track (`max_drop_log2fc`), which
        is what the compact table ranks its single target on.
        """
        gs = self.genes(chrom, element_id)
        if gs is None:
            return Response(None, NOT_CACHED)
        g = gs.get(gene)
        if g is None:
            return Response(None, NOT_IN_WINDOW)
        if cell is None:
            v = g.get("max_drop_log2fc")
            return Response(float(v), SCORED) if v is not None else Response(None, NOT_ON_TRACK)
        v = (g.get("by_cell") or {}).get(cell)
        return Response(float(v), SCORED) if v is not None else Response(None, NOT_ON_TRACK)

    def value(self, chrom: str, element_id: str, gene: str, cell: str | None = None) -> float | None:
        """The signed change, or None. `response()` says which silence a None is."""
        return self.response(chrom, element_id, gene, cell).value

    def at_bar(
        self, chrom: str, element_id: str, min_effect: float = 0.1, genes: set[str] | None = None
    ) -> list[tuple[str, float]] | None:
        """The window's genes at the bar, strongest first (`genes_at_bar`), optionally only those in
        `genes` (a coding set, say); None when the element is not cached."""
        got = genes_at_bar(self.element(chrom, element_id), min_effect)
        if got is None or genes is None:
            return got
        return [(g, v) for g, v in got if g in genes]

    def ranked(
        self, chrom: str, element_id: str, cell: str | None = None, by: str = "drop"
    ) -> list[tuple[str, float]]:
        """The window's genes ordered strongest first, as (gene, signed change); [] when not cached.

        This is the ranking whose head the compact table keeps, so a caller can ask where a
        measured gene sits in it rather than only whether it is the head. `by="drop"` orders on
        the predicted fall alone, largest fall first. `by="effect"` orders on the larger of the
        fall and the rise, which is the rule `predict.enhancer_target.predict_target` names a
        target with, so with `cell=None` its head is the element's `predicted` gene — an element
        whose head clears `MIN_EFFECT` is the one the compact table records and every other gene
        in the list is what the compact table dropped.
        """
        gs = self.genes(chrom, element_id) or {}
        out = []
        for name, g in gs.items():
            if cell is not None:
                v = (g.get("by_cell") or {}).get(cell)
                if v is None:
                    continue
                out.append((name, float(v)))
                continue
            drop, rise = g.get("max_drop_log2fc"), g.get("max_rise_log2fc")
            if drop is None:
                continue
            if by == "effect" and rise is not None and float(rise) > -float(drop):
                out.append((name, float(rise)))
            else:
                out.append((name, float(drop)))
        key = (lambda kv: -abs(kv[1])) if by == "effect" else (lambda kv: kv[1])
        return sorted(out, key=key)


# ---------------------------------------------------------------------------
# the one window reading every consumer gains, so none of them invents its own
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class WindowReading:
    """One element read BOTH ways: the compact head the table kept, and the window it came from.

    Registered in `attribution/onetarget2.py` (2026-10-03). Eighteen consumers under `genomeos/`
    read one target gene per element; this is the single record each of them gains, so that the
    window arrives in one shape rather than eighteen. Nothing here replaces a head: `head` is the
    compact table's own gene, carried through unchanged, and the window is carried BESIDE it.

    `at_bar` is None when the element is not in the response cache -- which is NOT an empty window
    and NOT a zero. A consumer must keep its compact reading for such an element and count it by
    name; `not_cached` says so, and `extra` is empty rather than misleading.

    The head and the window agree by construction and that agreement is measured, not assumed:
    `predict_target` and `genes_at_bar` apply the same size rule with the same tie order, and
    849,469 elements were compared with 0 head disagreements (2026-09-27, `d717b28`). The LIMIT of
    that measurement travels with it: it is one arm's elements and not a theorem, and 10 of 4,794
    elements had no `predicted_coding` field although the cache named a coding gene -- so
    `coding_head_agrees` can be False where `head_agrees` is True, and a consumer that depends on
    the CODING head specifically must read `coding_head_agrees` rather than assume it.
    """

    head: str | None
    coding_head: str | None
    at_bar: tuple[tuple[str, float], ...] | None
    coding_at_bar: tuple[tuple[str, float], ...] | None

    @property
    def not_cached(self) -> bool:
        return self.at_bar is None

    @property
    def genes(self) -> tuple[str, ...]:
        """Every gene at the bar, strongest first; empty when not cached or nothing is at the bar."""
        return tuple(g for g, _ in (self.at_bar or ()))

    @property
    def coding_genes(self) -> tuple[str, ...]:
        return tuple(g for g, _ in (self.coding_at_bar or ()))

    @property
    def extra(self) -> tuple[str, ...]:
        """At the bar and NOT the compact head: exactly what the one-target projection dropped."""
        return tuple(g for g in self.genes if g != self.head)

    @property
    def coding_extra(self) -> tuple[str, ...]:
        return tuple(g for g in self.coding_genes if g != self.coding_head)

    @property
    def head_agrees(self) -> bool | None:
        """Is the window's at-the-bar head the compact table's `predicted` gene? None: not cached.

        False refutes the registered invariant and stops the wave; it is never a result.
        """
        if self.at_bar is None:
            return None
        return (self.at_bar[0][0] if self.at_bar else None) == self.head

    @property
    def coding_head_agrees(self) -> bool | None:
        """The same for the coding head, which the 2026-09-27 result found is NOT equally safe."""
        if self.coding_at_bar is None:
            return None
        return (self.coding_at_bar[0][0] if self.coding_at_bar else None) == self.coding_head


def head_gene(element: dict[str, Any], coding: bool = True) -> str | None:
    """The compact table's own head gene for one element, read and not recomputed."""
    return ((element.get("predicted_coding" if coding else "predicted")) or {}).get("gene")


#: A SECOND function of this name exists: `attribution/unknown_scoring.py:129`'s
#: `window_reading(responses, chrom, element_id, coding)`, written by wave 0's lane, which takes an
#: element ID rather than the element and returns a plain dict rather than a `WindowReading`. The two
#: do not collide -- neither module imports the other -- but they are two shapes of one reading, which
#: is the problem this wave exists to stop, reappearing one level down. Unifying them is NOT done
#: here: unknown_scoring is wave 0's module and is in progress, and consolidating a live lane's reader
#: underneath it would be the kind of silent widening this registration refuses. Reported instead.
def window_reading(
    responses: ElementResponses | None,
    chrom: str,
    element: dict[str, Any],
    coding: set[str] | None = None,
    min_effect: float = 0.1,
) -> WindowReading:
    """One element's `WindowReading`. With `responses=None` nothing is opened and nothing is read.

    `responses=None` is the control every consumer keeps: the record then carries the compact heads
    alone and `not_cached` is True, so a consumer written against this helper reproduces its own
    committed figures exactly when it is called the old way.

    `coding` is the chromosome's coding symbols. Without it the coding window is not formed, because
    a coding head beside an uncoded window would invite a comparison between two different
    populations -- which is the defect this whole wave is about.
    """
    head, chead = head_gene(element, coding=False), head_gene(element, coding=True)
    if responses is None:
        return WindowReading(head, chead, None, None)
    eid = str(element["id"])
    bar = responses.at_bar(chrom, eid, min_effect)
    cbar = None if (bar is None or coding is None) else [t for t in bar if t[0] in coding]
    return WindowReading(
        head,
        chead,
        None if bar is None else tuple(bar),
        None if cbar is None else tuple(cbar),
    )


def prediction_window_reading(
    record: dict[str, Any] | None,
    head: str | None,
    coding_head: str | None = None,
    coding: set[str] | None = None,
    min_effect: float = 0.1,
) -> WindowReading:
    """The same `WindowReading`, for a consumer that reads the sweep's cache RECORD, not a table row.

    `window_reading` above takes a compact table element, whose head is a `predicted` field. Some
    consumers never see that table: they call `predict.enhancer_target.cached_prediction`, which
    opens the element's cached record, collapses the window with `predict_target` and returns one
    gene. `genome/regdiff.py` and `genome/regulation.py` are the two, named in the registration
    (`attribution/onetarget2.py`, THE_SHARED_READER) as the pair whose fix needs no new reader and
    no request, because the record is already in hand one call down.

    So this takes the record and the head the caller already read, and returns the one shape. It
    does NOT recompute the head: `head` is `cached_prediction`'s own gene, carried through, exactly
    as `window_reading` carries the table's. `record=None` is the control -- nothing was opened, or
    the element is not cached -- and the reading then carries the heads alone with `not_cached` True,
    which is what lets a consumer reproduce its committed figures when it is called the old way.

    `coding` is the chromosome's coding symbols, and without it the coding window is not formed,
    for the reason `window_reading` gives: a coding head beside an uncoded window would compare two
    different populations.
    """
    bar = genes_at_bar(record, min_effect)
    cbar = None if (bar is None or coding is None) else [t for t in bar if t[0] in coding]
    return WindowReading(
        head,
        coding_head,
        None if bar is None else tuple(bar),
        None if cbar is None else tuple(cbar),
    )


def tracks_at_bar(
    record: dict[str, Any] | None, min_effect: float = 0.1
) -> list[tuple[str, float, str | None]] | None:
    """`genes_at_bar` WITH the track each gene's size was taken on, as (gene, signed, tissue).

    AN AMENDMENT TO THE REGISTRATION, recorded in `onetarget2.AMENDMENTS`, because it widens what the
    shared reader provides. `WindowReading` carries the gene and the signed change and drops the
    track; `attribution/vista.py` cannot move without it, since its whole claim is whether the
    predicted TRACK GROUP agrees with the element's observed expression groups, and a gene with no
    track has no group to compare. Rather than let vista reimplement the size rule privately, the
    track is returned here beside the same gene list.

    The rule is `genes_at_bar`'s own, applied to the same fields in the same order: the larger of the
    predicted fall and the predicted rise, ties to the fall, and the tissue of whichever side won.
    That the two agree is not assumed: `tests/test_vista_window.py` asserts that the (gene, signed)
    projection of this list is EQUAL to `genes_at_bar`'s output, element for element, on fixtures and
    on real chromosome data, so the two cannot drift into two size rules.

    None means the element was not cached, exactly as in `genes_at_bar`.
    """
    if record is None:
        return None
    out = []
    for g in record.get("genes") or []:
        drop, rise = g.get("max_drop_log2fc"), g.get("max_rise_log2fc")
        if not g.get("gene") or drop is None:
            continue
        fall = -float(drop)
        up = float(rise) if rise is not None else float("-inf")
        if fall >= up:
            size, signed, tissue = fall, float(drop), g.get("max_drop_tissue")
        else:
            size, signed, tissue = up, up, g.get("max_rise_tissue")
        if size >= min_effect:
            out.append((g["gene"], round(signed, 4), size, tissue))
    out.sort(key=lambda t: -t[2])
    return [(name, signed, tissue) for name, signed, _, tissue in out]
