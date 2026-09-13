# SPDX-License-Identifier: AGPL-3.0-or-later
"""Enhancer activity measured in cell lines (ENCODE4 lentiMPRA) against what GenomeOS says, per cell.

VISTA asks whether a sequence is an enhancer in an embryo; GTEx asks which gene; a massively
parallel reporter assay asks how much a sequence drives transcription in a given cell line, for
tens of thousands of sequences at once. The Ahituv lab's ENCODE4 joint lentiMPRA library tested
the same 53,990 elements (200 bp, chosen from DNase peaks of each line) in K562, HepG2 and WTC11,
so activity can be compared across cells for the same sequence. The three element files (1.1 MB
each) are streamed once into data/knowledge/mpra, local, where the per-element rows also stay; the
committed result per chromosome is the summary.

Every element is read blind, per cell: which ENCODE cCRE class sits there; whether the cell's
own measured DNase (the reader) is open over it; how constrained it is; and, where the model is
on, AlphaGenome's predicted DNase for that cell line over the element. The score holds each of
these against the measured activity in the same cell (K562 against K562) and, separately, across
cells (HepG2's prediction against K562's activity), so the result says whether the layer knows
the cell and not only the element. Activity is log2(RNA/DNA); an element counts as active at or
above ACTIVE. Measured layers are `experimental`; the model's is `predicted`.
"""

from __future__ import annotations

import gzip
import json
import math
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

LIBRARY = "ENCSR106SZM"
FILES = {"K562": "ENCFF802FUV", "HepG2": "ENCFF475FKV", "WTC11": "ENCFF769REH"}
EXPERIMENTS = {"K562": "ENCSR203UFY", "HepG2": "ENCSR405QCT", "WTC11": "ENCSR336MKI"}
ENCODE_FILE = "https://www.encodeproject.org/files/{acc}/@@download/{acc}.bed.gz"
KNOWLEDGE = Path("data/knowledge/mpra")
ACTIVE = 1.0  # log2(RNA/DNA) at or above which an element counts as active
CONSTRAINED = 0.2
MODELLED = ("K562", "HepG2")  # cell lines with an AlphaGenome DNase track; WTC11 is measured only
EVIDENCE = (
    "experimental: ENCODE4 lentiMPRA (Ahituv lab, UCSF), joint library ENCSR106SZM in K562, HepG2 and WTC11, "
    "log2(RNA/DNA) per 200 bp element"
)
CLASS_ORDER = ("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3", "none")


@dataclass
class Element:
    chrom: str
    start: int
    end: int
    name: str
    activity: dict[str, float] = field(default_factory=dict)

    @property
    def key(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"


def fetch(knowledge: Path = KNOWLEDGE) -> dict[str, Path]:
    """The three element files, streamed once."""
    knowledge.mkdir(parents=True, exist_ok=True)
    out = {}
    for cell, acc in FILES.items():
        dest = knowledge / f"{acc}.bed.gz"
        if not dest.exists():
            req = urllib.request.Request(
                ENCODE_FILE.format(acc=acc), headers={"User-Agent": "GenomeOS/0.1 (stream)"}
            )
            with urllib.request.urlopen(req, timeout=300) as resp:  # noqa: S310
                dest.write_bytes(resp.read())
        out[cell] = dest
    return out


def parse(lines, cell: str, chrom: str | None, into: dict[str, Element]) -> None:
    """Rows of one cell's file folded into elements; forward and reverse copies are averaged."""
    sums: dict[str, list[float]] = {}
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 7 or (chrom and f[0] != chrom):
            continue
        key = f"{f[0]}:{f[1]}-{f[2]}"
        e = into.get(key)
        if e is None:
            e = into[key] = Element(f[0], int(f[1]), int(f[2]), f[3].removesuffix("_Reversed:"))
        sums.setdefault(key, []).append(float(f[6]))
    for key, vals in sums.items():
        into[key].activity[cell] = round(sum(vals) / len(vals), 4)


def load(chrom: str | None = None, knowledge: Path = KNOWLEDGE) -> list[Element]:
    paths = fetch(knowledge)
    into: dict[str, Element] = {}
    for cell, p in paths.items():
        with gzip.open(p, "rt") as fh:
            parse(fh, cell, chrom, into)
    out = [e for e in into.values() if e.activity]
    out.sort(key=lambda e: (e.chrom, e.start))
    return out


def annotate(
    elements: list[Element],
    ccres: list,
    peaks: dict[str, Any],
    stats=None,
) -> list[dict[str, Any]]:
    """Registry class, the reader's openness per cell and constraint for each element, blind to the assay."""
    import bisect

    cc = sorted(ccres, key=lambda c: c.start)
    starts = [c.start for c in cc]
    rows = []
    for k, e in enumerate(elements):
        i = bisect.bisect_left(starts, e.start - 5_000)
        classes = set()
        while i < len(cc) and cc[i].start < e.end:
            if cc[i].end > e.start:
                classes.add(cc[i].cls)
            i += 1
        cls = next((c for c in CLASS_ORDER if c in classes), "none")
        s = stats[k] if stats else None
        rows.append(
            {
                "key": e.key,
                "name": e.name,
                "start": e.start,
                "end": e.end,
                "activity": dict(e.activity),
                "active": {c: v >= ACTIVE for c, v in e.activity.items()},
                "ccre_class": cls,
                "open": {c: idx.covered_bp(e.start, e.end) > 0 for c, idx in peaks.items()},
                "constrained_fraction": round(s.fraction_above or 0.0, 3)
                if s is not None and s.bases
                else None,
            }
        )
    return rows


def save_rows(chrom: str, rows: list[dict[str, Any]], knowledge: Path = KNOWLEDGE) -> Path:
    """The per-element rows stay local (a few MB per chromosome); the committed result keeps the summary."""
    p = knowledge / f"rows_{chrom}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(rows))
    return p


def load_rows(chrom: str, knowledge: Path = KNOWLEDGE) -> list[dict[str, Any]]:
    p = knowledge / f"rows_{chrom}.json"
    return json.loads(p.read_text()) if p.exists() else []


def attach_predicted(rows: list[dict[str, Any]], predicted: dict[str, dict[str, float]]) -> None:
    for r in rows:
        r["predicted_dnase"] = predicted.get(r["key"])


def spearman(xs: list[float], ys: list[float]) -> float | None:
    """Rank correlation without scipy; ties get their average rank."""
    n = len(xs)
    if n < 3:
        return None

    def ranks(v: list[float]) -> list[float]:
        order = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            for k in range(i, j + 1):
                r[order[k]] = (i + j) / 2 + 1
            i = j + 1
        return r

    rx, ry = ranks(xs), ranks(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(rx, ry, strict=True))
    vx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    vy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return round(cov / (vx * vy), 3) if vx and vy else None


def _frac(rows: list[dict[str, Any]], pred) -> float | None:
    return round(sum(1 for r in rows if pred(r)) / len(rows), 3) if rows else None


def _precision_recall(rows: list[dict[str, Any]], flag, truth) -> dict[str, Any]:
    flagged = [r for r in rows if flag(r)]
    true = [r for r in rows if truth(r)]
    both = [r for r in flagged if truth(r)]
    return {
        "flagged": len(flagged),
        "precision": round(len(both) / len(flagged), 3) if flagged else None,
        "recall": round(len(both) / len(true), 3) if true else None,
        "base_rate": _frac(rows, truth),
    }


def summarise(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Per cell: activity against the registry, the reader, constraint and the predicted track, same cell
    and across cells."""
    cells = sorted({c for r in rows for c in r["activity"]})
    out: dict[str, Any] = {"elements": len(rows), "active_threshold_log2": ACTIVE, "cells": {}}
    for cell in cells:
        rs = [r for r in rows if cell in r["activity"]]
        active = lambda r, c=cell: r["active"].get(c, False)  # noqa: E731
        by_class = {}
        for cls in CLASS_ORDER:
            sub = [r for r in rs if r["ccre_class"] == cls]
            if sub:
                by_class[cls] = {"elements": len(sub), "active": _frac(sub, active)}
        entry: dict[str, Any] = {
            "elements": len(rs),
            "active": _frac(rs, active),
            "median_activity": sorted(r["activity"][cell] for r in rs)[len(rs) // 2] if rs else None,
            "active_by_ccre_class": by_class,
        }
        with_c = [r for r in rs if r["constrained_fraction"] is not None]
        if with_c:
            entry["active_when_constrained"] = _frac(
                [r for r in with_c if r["constrained_fraction"] >= CONSTRAINED], active
            )
            entry["active_when_not_constrained"] = _frac(
                [r for r in with_c if r["constrained_fraction"] < CONSTRAINED], active
            )
        reader = {}
        for other in sorted({c for r in rs for c in r["open"]}):
            reader[other] = _precision_recall(rs, lambda r, o=other: r["open"].get(o, False), active)
        if reader:
            entry["reader_open_predicts_active"] = {
                "same_cell": reader.get(cell),
                "other_cells": {k: v for k, v in reader.items() if k != cell},
            }
        pred = {}
        for other in MODELLED:
            sub = [r for r in rs if (r.get("predicted_dnase") or {}).get(other) is not None]
            if len(sub) >= 3:
                xs = [r["predicted_dnase"][other] for r in sub]
                ys = [r["activity"][cell] for r in sub]
                q = sorted(xs)
                lo, hi = q[len(q) // 4], q[(3 * len(q)) // 4]
                pred[other] = {
                    "elements": len(sub),
                    "spearman": spearman(xs, ys),
                    "active_in_top_quartile": _frac(
                        [r for r in sub if r["predicted_dnase"][other] >= hi], active
                    ),
                    "active_in_bottom_quartile": _frac(
                        [r for r in sub if r["predicted_dnase"][other] <= lo], active
                    ),
                }
        if pred:
            entry["predicted_dnase_vs_activity"] = {
                "same_cell": pred.get(cell),
                "other_cells": {k: v for k, v in pred.items() if k != cell},
            }
        out["cells"][cell] = entry
    # specificity between the two modelled lines: does the layer know which cell?
    a, b = MODELLED
    spec = [r for r in rows if a in r["active"] and b in r["active"] and r["active"][a] != r["active"][b]]
    if spec:
        on = lambda r: a if r["active"][a] else b  # noqa: E731
        off = lambda r: b if r["active"][a] else a  # noqa: E731
        rd = [r for r in spec if a in r["open"] and b in r["open"]]
        pd_ = [
            r
            for r in spec
            if (r.get("predicted_dnase") or {}).get(a) is not None and r["predicted_dnase"].get(b) is not None
        ]
        out["specific_elements"] = {
            "elements": len(spec),
            "active_in": {
                a: sum(1 for r in spec if r["active"][a]),
                b: sum(1 for r in spec if r["active"][b]),
            },
            "reader_open_in_the_active_cell_only": _frac(
                rd, lambda r: r["open"][on(r)] and not r["open"][off(r)]
            ),
            "reader_open_in_the_inactive_cell_only": _frac(
                rd, lambda r: r["open"][off(r)] and not r["open"][on(r)]
            ),
            "predicted_dnase_higher_in_the_active_cell": _frac(
                pd_, lambda r: r["predicted_dnase"][on(r)] > r["predicted_dnase"][off(r)]
            ),
            "predicted_judged": len(pd_),
        }
    out["evidence"] = EVIDENCE
    return out
