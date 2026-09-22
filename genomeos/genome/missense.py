# SPDX-License-Identifier: AGPL-3.0-or-later
"""The missense effect predicted rather than ranked: AlphaMissense scores for a person's variants.

The coding inventory ranks a person's missense variants by what UniProt says about the residue (a
site, a domain, nothing), which is annotation, not effect. AlphaMissense (DeepMind, 2023) predicts a
pathogenicity score for every possible missense variant of the human proteome, per transcript. Its
hg38 table is 643 MB compressed, 71 million rows; it is streamed once per person from DeepMind's public
bucket, only the rows for that person's missense variants are kept, and the kept rows stay under the
person's own git-ignored directory (`alphamissense.json`), never under data/results.

The scores enter as `predicted` evidence (confidence the score itself, capped at 0.7), with the model's
own classes (likely_pathogenic ≥ 0.564, likely_benign ≤ 0.34, ambiguous between); a variant absent from
the table (an indel, a non-canonical transcript, a stop) stays unscored and says so. AlphaMissense is
released under CC BY-NC-SA 4.0: research and personal use, not commercial.
"""

from __future__ import annotations

import gzip
import io
import json
import time
import urllib.request
from pathlib import Path
from typing import Any

URL = "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
EVIDENCE = "predicted: AlphaMissense (DeepMind 2023, CC BY-NC-SA 4.0) pathogenicity per missense variant"
CONFIDENCE_CAP = 0.7
CLASS_ORDER = {"likely_pathogenic": 0, "ambiguous": 1, "likely_benign": 2, None: 3}
Key = tuple[str, int, str, str]


def cache_path(name: str, root: Path) -> Path:
    return root / name / "alphamissense.json"


def stream_scores(
    wanted: set[Key], url: str = URL, progress=None, rows=None
) -> tuple[dict[str, list[dict[str, Any]]], dict[str, Any]]:
    """One pass over the table; the rows whose variant is in `wanted`, keyed "chrom:pos:ref:alt"."""
    t0 = time.time()
    if rows is None:
        req = urllib.request.Request(url, headers={"User-Agent": "GenomeOS/0.1 (stream)"})
        resp = urllib.request.urlopen(req, timeout=3600)  # noqa: S310
        rows = io.TextIOWrapper(gzip.GzipFile(fileobj=io.BufferedReader(resp, 1 << 20)), encoding="utf-8")
    out: dict[str, list[dict[str, Any]]] = {}
    seen = 0
    by_chrom_pos = {(c, p) for c, p, _, _ in wanted}
    for line in rows:
        if line.startswith("#"):
            continue
        seen += 1
        if progress and seen % 5_000_000 == 0:
            progress(f"AlphaMissense: {seen:,} rows read, {len(out):,} of {len(wanted):,} variants found")
        f = line.rstrip("\n").split("\t")
        if len(f) < 10:
            continue
        pos = int(f[1])
        if (f[0], pos) not in by_chrom_pos or (f[0], pos, f[2], f[3]) not in wanted:
            continue
        out.setdefault(f"{f[0]}:{pos}:{f[2]}:{f[3]}", []).append(
            {
                "transcript": f[6],
                "protein_variant": f[7],
                "score": float(f[8]),
                "class": f[9],
                "uniprot": f[5],
            }
        )
    cost = {"rows_read": seen, "seconds": round(time.time() - t0), "variants_found": len(out)}
    return out, cost


def scores_for(
    name: str, variants: list[dict[str, Any]], root: Path, progress=None, url: str = URL, rows=None
) -> dict[str, list[dict[str, Any]]]:
    """The person's scores, from the cache under their directory or one streaming pass."""
    p = cache_path(name, root)
    cached = json.loads(p.read_text()) if p.exists() else {"scores": {}, "asked": []}
    keys = {f"{v['chrom']}:{v['pos']}:{v['ref']}:{v['alt']}" for v in variants}
    missing = keys - set(cached["asked"])
    if missing:
        wanted = {(k.split(":")[0], int(k.split(":")[1]), k.split(":")[2], k.split(":")[3]) for k in missing}
        found, cost = stream_scores(wanted, url, progress, rows)
        cached["scores"].update(found)
        cached["asked"] = sorted(set(cached["asked"]) | missing)
        cached["last_pass"] = {**cost, "date": time.strftime("%Y-%m-%d"), "source": url}
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(cached))
    return cached["scores"]


def pick(entries: list[dict[str, Any]], transcript: str | None) -> dict[str, Any] | None:
    """The score on the person's transcript when the table has it, else the highest across transcripts."""
    if not entries:
        return None
    base = (transcript or "").split(".")[0]
    for e in entries:
        if base and e["transcript"].split(".")[0] == base:
            return {**e, "on_transcript": "the person's"}
    best = max(entries, key=lambda e: e["score"])
    return {**best, "on_transcript": "highest of the table's"}


def annotate(missense: list[dict[str, Any]], scores: dict[str, list[dict[str, Any]]]) -> dict[str, Any]:
    """Each missense row gains `predicted` (score, class, transcript); the rows are re-ranked by effect."""
    for m in missense:
        key = f"{m['chrom']}:{m['pos']}:{m['ref']}:{m['alt']}"
        chosen = pick(scores.get(key, []), m.get("transcript"))
        m["predicted"] = (
            {
                "score": chosen["score"],
                "class": chosen["class"],
                "confidence": round(min(CONFIDENCE_CAP, chosen["score"]), 3),
                "transcript": chosen["transcript"],
                "protein_variant": chosen["protein_variant"],
                "on_transcript": chosen["on_transcript"],
                "basis": EVIDENCE,
            }
            if chosen
            else None
        )
    ranked = sorted(
        missense,
        key=lambda m: (
            CLASS_ORDER[(m["predicted"] or {}).get("class")],
            -((m["predicted"] or {}).get("score") or 0.0),
            m["zygosity"] != "homozygous",
            m["gene"],
        ),
    )
    classes = {"likely_pathogenic": 0, "ambiguous": 0, "likely_benign": 0, "unscored": 0}
    for m in missense:
        classes[(m["predicted"] or {}).get("class") or "unscored"] += 1
    return {
        "missense_predicted": classes,
        "missense_likely_pathogenic_homozygous": sum(
            1
            for m in missense
            if (m["predicted"] or {}).get("class") == "likely_pathogenic" and m["zygosity"] == "homozygous"
        ),
        "missense_by_effect": ranked[:60],
        "missense_evidence": EVIDENCE,
    }
