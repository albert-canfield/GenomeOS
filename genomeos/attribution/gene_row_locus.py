# SPDX-License-Identifier: AGPL-3.0-or-later
"""Can a cached scorer answer's gene row be tied to a locus?

lane-identity resolved every one of the 440,589 compiled rules' targets to an Ensembl gene id and
compared each against the gene its deletion answer measured. 0 differ -- but the test did not fire:
where a target name has one locus on its chromosome both sides resolve to it whatever the window
says, and 440,035 of the 440,589 are of that kind. The rules a differing id was possible for are
the 547 whose name has two annotated loci inside the scorer window, and it recorded all 547 as
unresolvable, because both loci lie inside the window and the window does not choose between them.
Its stated repair was for the scoring chain to record the Ensembl gene id AlphaGenome returns per
gene row rather than only its name.

This module asks the prior question, from what is on disk and nothing else: does a cached gene row
already carry anything that ties it to a locus?

**What a gene row contains.** Every gene row in every cached answer carries exactly eight keys --
``gene``, ``n_tracks``, ``mean_log2fc``, ``max_drop_log2fc``, ``max_drop_tissue``,
``max_rise_log2fc``, ``max_rise_tissue``, ``by_cell`` -- and no ninth. There is no gene id, no
coordinate, no strand, no transcript and no gene type. :func:`row_key_vocabulary` establishes that
by streaming, never by loading a chromosome whole.

**Where it is dropped.** Twice, both in committed code:

``genomeos/predict/alphagenome_adapter.py:227`` reads one column of the response's gene axis,
``genes = list(adata.obs.get("gene_name", []))``. Whatever else that axis carries is never read, and
the row's own position in the axis is not kept either, so by the time the answer reaches
:func:`genomeos.predict.enhancer_target.aggregate` a gene is a bare string.

``genomeos/predict/enhancer_target.py:198`` then keys the accumulator by that string
(``by_gene.setdefault(gene, ...)``), so two response rows carrying one ``gene_name`` merge into one
cached row: ``n_tracks`` is summed, ``mean_log2fc`` is averaged over both, ``max_drop_log2fc`` is the
extreme of either with nothing saying which, and ``by_cell`` keeps "the last value read wins"
(line 215) -- one locus's number standing for both.

**That merge is measurable, and it is what proves the response distinguished the loci.** One response
row can contribute at most one value per column, and the response's own column count is recorded:
``model.tracks = 371`` under one ``tracks_sha256`` over every answer that carries a run record. So a
row whose ``n_tracks`` exceeds 371 cannot be one response row. :data:`TRACKS` fixes that bound.

**What that does and does not settle.** :func:`classify_element` takes every cached element, every
gene name the response returned, and the GENCODE v50 loci of that name whose gene BODY lies inside
the scorer window -- lane-identity's own candidate rule, at its own ``SCORER_HALF_WINDOW``. Where two
or more such loci exist the pair is of the kind the 547 are, and the classes are:

``merged``
    the symbol row's ``n_tracks`` exceeds :data:`TRACKS`. The response held both loci and the chain
    pooled them. Nothing in the pooled numbers says which contributed what, so this is not settled.

``pinned``
    the element ALSO carries rows named by the bare versionless Ensembl ids of exactly (loci - 1) of
    those loci, and the symbol row is not merged. AlphaGenome's ``gene_name`` for those loci IS their
    id, because its annotation gives them no symbol, so the symbol row is the one locus left and is
    named by elimination. **This is settled, from bytes already on disk.**

``partly_pinned``
    some but not all of the other loci appear by id: narrows, does not settle.

``nothing``
    no locus-identifying signal at all.

``pinned`` is locus-level and not symbol-level, which is the test that matters here: the two loci
share one GENCODE symbol and the cache gives them DIFFERENT names, and no property of a symbol can
differ between two loci that have the same symbol. Contrast ``predicted_coding``, which restricts the
chain's target choice to symbols that are protein_coding somewhere on the chromosome -- a property of
the SYMBOL, which is why it identifies no locus and resolved none of the 547.

**And the condition on it.** Pinning by elimination assumes AlphaGenome's gene set holds the same
loci of that name in that window as GENCODE v50 does. It does not always: MATR3, NOX5 and ZNF724 have
one GENCODE v50 locus each and yet their cached rows are merges of two response rows, so AlphaGenome's
annotation carries loci GENCODE v50 does not. That is the same caveat ``gene_identity`` already
records -- "the annotation AlphaGenome scored with and GENCODE v50 are not the same" -- and it bounds
any repair built on ids the service returns.

Nothing here reads the compiled rules, resolves any rule, or touches ``gene_identity.json``. It
counts at the (element, name) denominator, which is **not** the 547: see :func:`what_this_is_not`.
"""

from __future__ import annotations

import collections
import gzip
import re
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any

from genomeos.attribution.gene_identity import HALF_WINDOW, GeneRow, annotation

#: The response's own RNA-seq column count, as the chain records it per answer
#: (``model.tracks``, under a single ``tracks_sha256``). One response row contributes at most one
#: value per column, so a cached gene row above this cannot be one response row. The chain asked
#: with ``threshold=0.0`` (``scripts/enhancer_targets_all.py``, ``worker_scorer``), which drops a
#: value of exactly zero, so a single row's ``n_tracks`` is at most this and often just below.
TRACKS = 371

#: The eight keys every cached gene row carries, and the whole of what one says. None of them names
#: a locus: ``gene`` is AlphaGenome's ``gene_name``, the rest are magnitudes and track labels.
GENE_ROW_KEYS = (
    "by_cell",
    "gene",
    "max_drop_log2fc",
    "max_drop_tissue",
    "max_rise_log2fc",
    "max_rise_tissue",
    "mean_log2fc",
    "n_tracks",
)

#: Where the locus-identifying information is dropped, by file and line, as read 2026-10-02.
DROPPED_AT = (
    {
        "file": "genomeos/predict/alphagenome_adapter.py",
        "line": 227,
        "code": 'genes = list(adata.obs.get("gene_name", []))',
        "what_is_dropped": (
            "every column of the response's gene axis but gene_name, and the row's own position in "
            "that axis. After this line a gene is a bare string and no later code can recover which "
            "row of the response it was"
        ),
    },
    {
        "file": "genomeos/predict/enhancer_target.py",
        "line": 198,
        "code": "g = by_gene.setdefault(gene, {...})",
        "what_is_dropped": (
            "the separation between two response rows that carry the same gene_name: they merge into "
            "one cached row, n_tracks summed, mean_log2fc averaged over both, max_drop_log2fc the "
            "extreme of either with nothing saying which, and by_cell 'the last value read wins' "
            "(line 215), one locus's number standing for both"
        ),
    },
)

#: A bare versionless Ensembl gene id, which is what AlphaGenome's ``gene_name`` holds for a gene its
#: annotation gives no symbol. The one locus-identifying string that reaches the cache.
ENSG = re.compile(r"^ENSG\d+$")

#: Matches a cached gene row's first two fields in the raw bytes. The row is written by
#: ``json.dumps`` of a dict built in insertion order, so ``gene`` and ``n_tracks`` are adjacent and
#: can be read without parsing -- which is how a 59 MB chromosome is read without loading it whole.
ROW = re.compile(rb'"gene":\s*"([^"]{1,64})",\s*"n_tracks":\s*(\d+)')

#: Any JSON object key in the raw bytes, for the vocabulary scan.
KEY = re.compile(rb'"([^"\\]{1,64})"\s*:')

CHUNK = 1 << 22

MERGED = "the cached row is a merge of more than one response row, so the numbers are pooled"
PINNED = "every other locus of that name in the window is named by its bare Ensembl id in this element"
PARTLY = "some but not all of the other loci are named by their bare Ensembl id"
NOTHING = "no locus-identifying signal in the element's own rows"

CLASSES = {"merged": MERGED, "nothing": NOTHING, "partly_pinned": PARTLY, "pinned": PINNED}


def _chunks(path: Path, overlap: int) -> Iterator[tuple[int, bytes]]:
    """The file's bytes in chunks, each prefixed with the last `overlap` bytes of the one before, and
    the absolute offset that prefix starts at. A match is keyed by its absolute offset so a match
    inside the overlap is counted once, and `overlap` must exceed the longest match."""
    fn = gzip.open if path.suffix == ".gz" else open
    tail = b""
    off = 0
    with fn(path, "rb") as fh:  # type: ignore[operator]
        while True:
            buf = fh.read(CHUNK)
            if not buf:
                return
            yield off - len(tail), tail + buf
            off += len(buf)
            tail = (tail + buf)[-overlap:]


def row_key_vocabulary(paths: Iterable[Path]) -> dict[str, int]:
    """Every distinct JSON object key in those files, with how often it occurs, by streaming.

    This is the test for a ninth gene-row field. A key that carried locus-identifying information
    per gene row would have to appear here, however few rows held it. The element ids are keys too
    (an archive is a dict keyed by element id), so the result mixes schema keys with ids; a caller
    separates them by frequency or by shape. Counts are approximate where a key falls inside the
    overlap between chunks; the key SET is complete.
    """
    out: collections.Counter[str] = collections.Counter()
    for p in paths:
        seen: set[int] = set()
        for base, data in _chunks(Path(p), 80):
            for m in KEY.finditer(data):
                a = base + m.start()
                if a in seen:
                    continue
                seen.add(a)
                out[m.group(1).decode()] += 1
    return dict(out)


def row_scan(paths: Iterable[Path]) -> dict[str, Any]:
    """Every cached gene row's name and ``n_tracks``, summarised, by streaming.

    ``merged_rows`` counts the rows above :data:`TRACKS`, which no single response row can reach.
    ``merged_at_least`` bands them by ``ceil(n_tracks / TRACKS)``, which is a LOWER BOUND on how many
    response rows merged and not the number: the sweep's ``threshold=0.0`` drops a value of exactly
    zero, so a constituent row can contribute fewer than :data:`TRACKS` values and two rows can land
    anywhere up to twice it. ``ensembl_named_rows`` counts the rows whose name is a bare Ensembl id,
    the only locus-identifying string the cache holds.
    """
    rows = 0
    merged = 0
    ensembl = 0
    band: collections.Counter[int] = collections.Counter()
    names: collections.Counter[str] = collections.Counter()
    for p in paths:
        seen: set[int] = set()
        for base, data in _chunks(Path(p), 120):
            for m in ROW.finditer(data):
                a = base + m.start()
                if a in seen:
                    continue
                seen.add(a)
                rows += 1
                name = m.group(1).decode()
                nt = int(m.group(2))
                if ENSG.match(name):
                    ensembl += 1
                if nt > TRACKS:
                    merged += 1
                    band[-(-nt // TRACKS)] += 1
                    names[name] += 1
    return {
        "gene_rows": rows,
        "merged_rows": merged,
        "merged_share": round(merged / rows, 6) if rows else None,
        "merged_at_least": dict(sorted(band.items())),
        "merged_names": dict(names.most_common()),
        "ensembl_named_rows": ensembl,
        "what_merged_rows_is_not": (
            "it is not a count of wrong rows, of wrong predictions or of wrong rules. It counts the "
            "cached rows that pooled two or more response rows under one gene_name, which is a loss "
            "of a distinction the response made, not an error in a number"
        ),
    }


def repeated_names(chrom: str) -> dict[str, list[GeneRow]]:
    """The GENCODE v50 gene names of `chrom` that more than one gene carries, with those genes."""
    by: dict[str, list[GeneRow]] = collections.defaultdict(list)
    for g in annotation(chrom):
        if g.symbol:
            by[g.symbol].append(g)
    return {n: gs for n, gs in by.items() if len(gs) > 1}


def coding_names(chrom: str) -> set[str]:
    """Every GENCODE v50 name of `chrom` with at least one protein_coding locus.

    This is a property of the NAME and not of a locus, and it is exactly the property
    ``predicted_coding`` restricts the chain's target choice by -- which is why it picks out no
    locus. It is used here only as the generous filter on which names a compiled target could be.
    """
    return {g.symbol for g in annotation(chrom) if g.symbol and g.coding}


def loci_in_window(loci: list[GeneRow], start: int, end: int, half: int = HALF_WINDOW) -> list[GeneRow]:
    """Those loci whose gene BODY overlaps the scorer window: the element midpoint plus or minus
    `half`. The rule and the half-window are lane-identity's, imported rather than restated."""
    mid = (start + end) // 2
    lo, hi = mid - half, mid + half
    return [g for g in loci if g.start < hi and g.end > lo]


def classify_element(
    record: dict[str, Any],
    repeated: dict[str, list[GeneRow]],
    coding: set[str],
    half: int = HALF_WINDOW,
) -> list[dict[str, Any]]:
    """One cached element's same-named-locus ambiguities, each classified. See the module docstring."""
    start, end = record.get("start"), record.get("end")
    if start is None or end is None:
        return []
    rows = record.get("genes") or []
    present = {r["gene"] for r in rows if r.get("gene")}
    out = []
    for r in rows:
        name = r.get("gene")
        if not name or name not in repeated:
            continue
        inwin = loci_in_window(repeated[name], start, end, half)
        if len(inwin) < 2:
            continue
        by_id = {g.gene_id for g in inwin} & present
        if (r.get("n_tracks") or 0) > TRACKS:
            verdict, settles = "merged", None
        elif len(by_id) == len(inwin) - 1:
            verdict = "pinned"
            settles = next(g.gene_id for g in inwin if g.gene_id not in by_id)
        elif by_id:
            verdict, settles = "partly_pinned", None
        else:
            verdict, settles = "nothing", None
        out.append(
            {
                "element": record.get("id"),
                "chrom": record.get("chrom"),
                "start": start,
                "end": end,
                "name": name,
                "loci_in_window": len(inwin),
                "gene_ids_in_window": sorted(g.gene_id for g in inwin),
                "gene_types_in_window": sorted(g.gene_type for g in inwin),
                "named_by_bare_ensembl_id_in_this_element": sorted(by_id),
                "n_tracks": r.get("n_tracks"),
                "verdict": verdict,
                "the_symbol_row_is_then": settles,
                "name_is_protein_coding_somewhere_on_the_chromosome": name in coding,
            }
        )
    return out


def tally(pairs: Iterable[dict[str, Any]]) -> dict[str, Any]:
    """The classes at every denominator, and which field settles the ones that are settled."""
    c: collections.Counter[str] = collections.Counter()
    pinned_names: collections.Counter[str] = collections.Counter()
    loci: collections.Counter[int] = collections.Counter()
    for p in pairs:
        v = p["verdict"]
        c["ambiguous"] += 1
        c[v] += 1
        loci[p["loci_in_window"]] += 1
        if p["name_is_protein_coding_somewhere_on_the_chromosome"]:
            c["ambiguous_coding_name"] += 1
            c[f"{v}_coding_name"] += 1
        if v == "pinned":
            pinned_names[p["name"]] += 1
    return {
        "ambiguous_element_name_pairs": c["ambiguous"],
        "loci_in_window": dict(sorted(loci.items())),
        "by_verdict": {k: c[k] for k in sorted(CLASSES)},
        "verdict_meanings": dict(CLASSES),
        "settled_by_what_is_on_disk": c["pinned"],
        "settled_by_which_field": (
            "the gene row's own `gene` field, where AlphaGenome's gene_name for a locus IS its bare "
            "versionless Ensembl id because its annotation gives that locus no symbol. Nothing else "
            "in the eight-key row carries a locus"
        ),
        "restricted_to_names_with_a_protein_coding_locus_on_the_chromosome": {
            "ambiguous": c["ambiguous_coding_name"],
            **{k: c[f"{k}_coding_name"] for k in sorted(CLASSES)},
        },
        "distinct_pinned_names": len(pinned_names),
        "pinned_names": dict(pinned_names.most_common()),
    }


def what_this_is_not() -> dict[str, str]:
    """Every figure this module reports, and what it is not. A count is not a measurement of the
    thing one wants, and the denominators here are not lane-identity's."""
    return {
        "the_denominator": (
            "every figure is counted over (cached element, gene name) pairs -- one per gene the "
            "response returned in each cached element's window. It is NOT the 440,589 compiled rules "
            "and it is NOT the 547. One rule names one element and one target, so one rule is at most "
            "one such pair, but this module does not read the compiled rules and therefore "
            "establishes no figure at the rule denominator"
        ),
        "the_547": (
            "this module does not recount, re-resolve or revise the 547. What it establishes about "
            "them is a property of the class they belong to: a compiled target comes from a "
            "`predicted_coding` block, so its name has a protein_coding locus on the chromosome, and "
            "no pair whose name has one is settled by anything on disk"
        ),
        "merged_rows": (
            "a merge is a lost distinction, not a wrong number. It does not say the prediction is "
            "wrong, that any rule is wrong, or what share of predictions is affected"
        ),
        "pinned": (
            "pinning names the locus the symbol row stands for by elimination, under the assumption "
            "that AlphaGenome's gene set holds the same loci of that name in that window as GENCODE "
            "v50 does. MATR3, NOX5 and ZNF724 show that assumption failing: one GENCODE v50 locus "
            "each, and a cached row that merges two response rows. So `pinned` is a reading under a "
            "stated condition, not an established identity"
        ),
        "the_eight_keys": (
            "that a gene row has eight keys is read from the bytes of every cached answer, not from "
            "a docstring. It does not establish what the AlphaGenome response contained, only what "
            "the chain kept of it"
        ),
        "no_interval": (
            "these are exhaustive counts over the files on disk, not an estimate from a sample, so "
            "no confidence interval is computed and none would mean anything"
        ),
    }
