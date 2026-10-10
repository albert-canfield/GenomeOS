# SPDX-License-Identifier: AGPL-3.0-or-later
"""The reader's call on one gene in every biosample, with a call made by the marks told apart.

Reader v1 (`genome/reader.py`) calls each coding gene's promoter (TSS ± 1 kb) per biosample and
chromosome. Four outcomes come from its result files, and one more is the absence of a call:

- ``open``: a DNase peak over the promoter, not poised. Read, on the DNase file.
- ``read_by_marks``: no DNase peak, but H3K4me3 and H3K27ac peaks over the promoter. Read, on the
  histone marks alone: a shallow DNase experiment misses promoters, so this count is inverse DNase
  depth by construction (reader.py's depth note).
- ``poised``: a DNase peak with H3K27me3 and no H3K27ac. Not read; held, not shut.
- ``closed``: no DNase peak, and the marks do not carry both H3K4me3 and H3K27ac. Not read.
- ``unknown``: the reader makes no call, and the reason says why. It is never drawn as closed.

The Blocks lane's ``<cell>_read`` attribute folds the first two into "read", so the call made by the
marks was shown nowhere. This module reads the same result files and keeps them apart. A result file
lists the closed, poised and read-by-marks genes in full and the open ones only as its top 25 by
signal, so "open" for any other gene is read from absence. Absence means open only when the reader
called the same coding genes the annotation holds; the coding-gene count is checked, and a gene on
none of the lists of a reader whose count differs is UNKNOWN, not open.

Read is the promoter accessible or marked: necessary for transcription, not proof of it.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from genomeos.genome.reader import slug

#: the calls, in the order the Blocks card lists them
CALLS = ("open", "read_by_marks", "poised", "closed", "unknown")

#: whether a call counts as read (None: no call)
READ = {"open": True, "read_by_marks": True, "poised": False, "closed": False, "unknown": None}

#: what each call rests on
BASIS = {
    "open": "DNase",
    "read_by_marks": "marks",
    "poised": "DNase and marks",
    "closed": "DNase and marks",
    "unknown": None,
}

#: plain words for each call; none says the gene is transcribed
GLOSS = {
    "open": "read: a DNase peak over the promoter (TSS ± 1 kb), and not poised",
    "read_by_marks": "read by marks: no DNase peak over the promoter, but H3K4me3 and H3K27ac peaks there",
    "poised": "not read: a DNase peak with H3K27me3 and no H3K27ac over the promoter (held, not shut)",
    "closed": "not read: no DNase peak over the promoter, and the marks do not read it",
    "unknown": "no call: the reason says why; not closed",
}

#: every reason an UNKNOWN can carry, in plain words
REASONS = {
    "not_in_annotation": "the gene is not in the annotation of this chromosome",
    "not_protein_coding": "the reader calls protein-coding promoters only",
    "no_reader_result": "no reader result for this biosample on this chromosome",
    "coding_set_differs": (
        "the gene is on none of the reader's lists, and the reader called a different number of coding "
        "genes than the annotation holds, so absence does not show the promoter open"
    ),
}

EVIDENCE = (
    "inferred, from reader v1 result files (data/results/reader_<biosample>_<chrom>.json): ENCODE DNase-seq "
    "peaks and, where read, H3K4me3, H3K27ac and H3K27me3 Histone ChIP-seq peaks over the promoter "
    "(TSS ± 1 kb). Read is the promoter accessible or marked: necessary for transcription, not proof of "
    "it. A call by the marks counts promoters the DNase file missed, so it rises as DNase depth falls."
)


def _panel(results_dir: Path, chrom: str) -> list[tuple[str, Path]]:
    """(biosample, file) for every biosample the reader holds: the genome-wide panel's order first,
    then any other reading on this chromosome. The same-chromosome comparison files are not readings."""
    out: list[tuple[str, Path]] = []
    seen: set[str] = set()
    gw = results_dir / "reader_genome_wide.json"
    if gw.exists():
        for name in json.loads(gw.read_text()).get("cell_types") or []:
            out.append((name, results_dir / f"reader_{slug(name)}_{chrom}.json"))
            seen.add(slug(name))
    for f in sorted(results_dir.glob(f"reader_*_{chrom}.json")):
        s = f.name[len("reader_") : -len(f"_{chrom}.json")]
        if s in seen or "_vs_" in s:
            continue
        r = json.loads(f.read_text())
        if "silent_genes" not in r:
            continue
        out.append((r.get("cell_type") or s, f))
        seen.add(s)
    return out


def _unknown(reason: str, detail: str) -> dict[str, Any]:
    return {"call": "unknown", "read": None, "basis": None, "reason": reason, "reason_text": detail}


def call_one(
    symbol: str, gene_type: str | None, chrom: str, coding_on_chrom: int | None, r: dict | None
) -> dict[str, Any]:
    """The reader's call on one gene from one result file (None: no file)."""
    if gene_type is None:
        return _unknown("not_in_annotation", f"{symbol} is not in the annotation of {chrom}")
    if gene_type != "protein_coding":
        return _unknown(
            "not_protein_coding", f"the reader calls protein-coding promoters only; {symbol} is {gene_type}"
        )
    if r is None or "silent_genes" not in r:
        return _unknown("no_reader_result", f"no reader result on {chrom}")
    marks = r.get("marks_used")
    note = None if marks else "marks not read here: poised and read by marks were not tested"
    out: dict[str, Any] = {"marks_used": marks, "signal": None, "by_absence": False, "note": note}
    if symbol in set(r.get("read_by_marks") or []):
        call = "read_by_marks"
    elif symbol in set(r.get("poised_genes") or []):
        call = "poised"
    elif symbol in set(r["silent_genes"]):
        call = "closed"
    else:
        top = {t["gene"]: t.get("signal") for t in r.get("top_read") or []}
        if symbol in top:
            call, out["signal"] = "open", top[symbol]
        elif coding_on_chrom is None or r.get("coding_genes") != coding_on_chrom:
            return _unknown(
                "coding_set_differs",
                f"{symbol} is on none of the reader's lists, and the reader called "
                f"{r.get('coding_genes')} coding genes on {chrom} against the annotation's {coding_on_chrom}",
            )
        else:
            call, out["by_absence"] = "open", True
    basis = BASIS[call] if marks else "DNase"
    return {"call": call, "read": READ[call], "basis": basis, "reason": None, "reason_text": None, **out}


def read_calls(
    symbol: str,
    chrom: str,
    gene_type: str | None,
    coding_on_chrom: int | None,
    results_dir: Path,
    reader_terms: set[str] | None = None,
) -> dict[str, Any]:
    """The reader's call on `symbol` in every biosample it holds on `chrom`.

    `gene_type` is the annotation's type for the gene (None: not annotated) and `coding_on_chrom` the
    annotation's protein-coding gene count on the chromosome, against which absence is checked.
    `reader_terms` names the biosamples compile and readergate gate on; a reading outside it is shown
    and flagged, not dropped."""
    if reader_terms is None:
        from genomeos.attribution.context_evidence import READER_TERMS

        reader_terms = set(READER_TERMS)
    calls = []
    for name, f in _panel(results_dir, chrom):
        r = json.loads(f.read_text()) if f.exists() else None
        calls.append(
            {
                "biosample": name,
                "result": f.stem,
                "in_reader_terms": name in reader_terms,
                **call_one(symbol, gene_type, chrom, coding_on_chrom, r),
            }
        )
    counts = {c: sum(1 for x in calls if x["call"] == c) for c in CALLS}
    return {
        "gene": symbol,
        "chrom": chrom,
        "gene_type": gene_type,
        "biosamples": len(calls),
        "read": counts["open"] + counts["read_by_marks"],
        "counts": counts,
        "calls": calls,
        "gloss": GLOSS,
        "reasons": REASONS,
        "reader_terms": sorted(reader_terms),
        "evidence": EVIDENCE,
    }
