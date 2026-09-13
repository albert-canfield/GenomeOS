# SPDX-License-Identifier: AGPL-3.0-or-later
"""The predicted target gene held against GTEx eQTLs, the measured ground truth for "which gene".

VISTA says whether a sequence is an enhancer and where; it cannot say which gene. GTEx can: a
significant cis-eQTL is a variant whose alleles go with the expression of a gene in a tissue, in
hundreds of donors. When such a variant sits inside an element, the gene it moves is a measured
target of that element (with the usual caveat that association is not mechanism and linkage
carries the signal a little way). This module streams GTEx v8's single-tissue eQTL archive (1.56
GB, 49 tissues, about 71 million significant pairs) from its public bucket, member by member over
HTTP ranges, keeps only the pairs whose variant falls inside an element GenomeOS has a prediction
for, and discards the rest. What stays is a few megabytes under data/knowledge/gtex, local.

The score then asks, element by element: of the genes GTEx ties to this element, is the one the
deletion model named among them, and is the one the node model named (the nearest coding TSS)?
Evidence: the eQTLs are `experimental`, the comparison is a test of `predicted` and `inferred`.
"""

from __future__ import annotations

import bisect
import gzip
import io
import json
import time
import urllib.request
from pathlib import Path
from typing import Any

GTEX_EQTL_TAR = (
    "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL.tar"
)
KNOWLEDGE = Path("data/knowledge/gtex")
EVIDENCE = "experimental: GTEx v8 single-tissue cis-eQTLs (significant variant-gene pairs), hg38"
USER_AGENT = "GenomeOS/0.1 (stream)"
MAX_ELEMENT = 8_000  # elements are never longer, so the interval search looks back this far
HIT_COLUMNS = ("tissue", "chrom", "pos", "ref", "alt", "gene_id", "slope", "pval_nominal", "elements")


class Intervals:
    """Elements per chromosome, searched by position."""

    def __init__(self) -> None:
        self._by_chrom: dict[str, list[tuple[int, int, str]]] = {}
        self._starts: dict[str, list[int]] = {}
        self.n = 0

    def add(self, chrom: str, start: int, end: int, element_id: str) -> None:
        self._by_chrom.setdefault(chrom, []).append((start, end, element_id))
        self.n += 1

    def freeze(self) -> Intervals:
        for c, rows in self._by_chrom.items():
            rows.sort()
            self._starts[c] = [r[0] for r in rows]
        return self

    def at(self, chrom: str, pos: int) -> list[str]:
        """Element ids covering the 1-based position."""
        rows = self._by_chrom.get(chrom)
        if not rows:
            return []
        i = bisect.bisect_right(self._starts[chrom], pos - 1) - 1
        out = []
        while i >= 0 and rows[i][0] > pos - 1 - MAX_ELEMENT:
            s, e, eid = rows[i]
            if s < pos <= e:
                out.append(eid)
            i -= 1
        return out


class _Slice:
    """A bounded read() over an open stream, so gzip stops at the member's end."""

    def __init__(self, fh, size: int) -> None:
        self.fh = fh
        self.left = size

    def read(self, n: int = -1) -> bytes:
        if self.left <= 0:
            return b""
        n = self.left if n < 0 else min(n, self.left)
        chunk = self.fh.read(n)
        self.left -= len(chunk)
        return chunk


def _open_range(source: str, offset: int, size: int, timeout: int = 600):
    if source.startswith(("http://", "https://")):
        req = urllib.request.Request(
            source, headers={"Range": f"bytes={offset}-{offset + size - 1}", "User-Agent": USER_AGENT}
        )
        return _Slice(urllib.request.urlopen(req, timeout=timeout), size)  # noqa: S310
    fh = open(source, "rb")  # noqa: SIM115
    fh.seek(offset)
    return _Slice(fh, size)


def tar_members(source: str = GTEX_EQTL_TAR) -> list[tuple[str, int, int]]:
    """(name, data offset, size) of every member, read from the 512-byte headers alone."""
    out = []
    off = 0
    while True:
        h = _open_range(source, off, 512).read(512)
        if len(h) < 512 or h[0] == 0:
            break
        name = h[0:100].split(b"\0", 1)[0].decode()
        size = int((h[124:136].split(b"\0", 1)[0].strip() or b"0"), 8)
        prefix = h[345:500].split(b"\0", 1)[0].decode() if h[257:262] == b"ustar" else ""
        if h[156:157] in (b"0", b"\0"):
            out.append((f"{prefix}/{name}" if prefix else name, off + 512, size))
        off += 512 + ((size + 511) // 512) * 512
    return out


def stream_pairs(source: str, offset: int, size: int):
    """Rows of one significant-pairs member: (chrom, pos, ref, alt, gene_id, slope, pval)."""
    raw = _open_range(source, offset, size)
    with gzip.GzipFile(fileobj=raw) as gz, io.TextIOWrapper(gz, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {k: i for i, k in enumerate(header)}
        iv, ig, isl, ip = col["variant_id"], col["gene_id"], col["slope"], col["pval_nominal"]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            v = f[iv].split("_")  # chr1_13550_G_A_b38
            if len(v) < 4:
                continue
            yield v[0], int(v[1]), v[2], v[3], f[ig].split(".")[0], f[isl], f[ip]


def tissue_of(member: str) -> str:
    return member.rsplit("/", 1)[-1].split(".v8.")[0]


def hits_path(tissue: str, knowledge: Path = KNOWLEDGE) -> Path:
    return knowledge / f"hits_{tissue}.tsv"


def distil(
    intervals: Intervals,
    knowledge: Path = KNOWLEDGE,
    source: str = GTEX_EQTL_TAR,
    tissues: list[str] | None = None,
    progress=None,
) -> dict[str, Any]:
    """Stream every tissue's significant pairs once; keep the pairs inside the elements; a file per tissue."""
    knowledge.mkdir(parents=True, exist_ok=True)
    members = [m for m in tar_members(source) if "signif_variant_gene_pairs" in m[0]]
    if tissues:
        members = [m for m in members if tissue_of(m[0]) in tissues]
    t0 = time.time()
    scanned = kept = 0
    mb = 0.0
    done = []
    for name, off, size in members:
        tissue = tissue_of(name)
        dest = hits_path(tissue, knowledge)
        if dest.exists():
            done.append(tissue)
            continue
        n = k = 0
        tmp = dest.with_suffix(".part")
        with tmp.open("w") as out:
            out.write("\t".join(HIT_COLUMNS) + "\n")
            for chrom, pos, ref, alt, gene, slope, pval in stream_pairs(source, off, size):
                n += 1
                ids = intervals.at(chrom, pos)
                if ids:
                    k += 1
                    out.write(
                        f"{tissue}\t{chrom}\t{pos}\t{ref}\t{alt}\t{gene}\t{slope}\t{pval}\t{','.join(ids)}\n"
                    )
        tmp.rename(dest)
        scanned += n
        kept += k
        mb += size / 1e6
        done.append(tissue)
        if progress:
            progress(f"{tissue}: {n:,} pairs, {k:,} inside elements ({len(done)}/{len(members)} tissues)")
    summary = {
        "source": source,
        "tissues": len(done),
        "pairs_scanned_this_run": scanned,
        "hits_this_run": kept,
        "mb_streamed_this_run": round(mb, 1),
        "seconds": round(time.time() - t0),
        "elements_indexed": intervals.n,
        "kept_where": str(knowledge),
        "evidence": EVIDENCE,
    }
    (knowledge / "distil_summary.json").write_text(json.dumps(summary, indent=1))
    return summary


def load_hits(knowledge: Path = KNOWLEDGE, chrom: str | None = None) -> dict[str, list[dict[str, Any]]]:
    """Hits per element id, over every tissue distilled so far."""
    out: dict[str, list[dict[str, Any]]] = {}
    for p in sorted(knowledge.glob("hits_*.tsv")):
        with p.open() as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if chrom and f[1] != chrom:
                    continue
                row = {
                    "tissue": f[0],
                    "chrom": f[1],
                    "pos": int(f[2]),
                    "gene_id": f[5],
                    "slope": float(f[6]),
                    "pval": float(f[7]),
                }
                for eid in f[8].split(","):
                    out.setdefault(eid, []).append(row)
    return out


def _norm(t: str | None) -> str:
    return (t or "").lower().replace("_", " ").replace("-", " ").strip()


def tissue_matches(predicted: str | None, gtex: str) -> bool | None:
    """Whether a predicted track name is the GTEx tissue; None when the track is not a GTEx-style name."""
    p, g = _norm(predicted), _norm(gtex)
    if not p:
        return None
    if "_" not in (predicted or ""):
        return True if p == g else None
    if p == g:
        return True
    head = g.split(" ")[0]  # brain, heart, liver, lung, ...
    if head in ("brain", "heart", "liver", "lung", "kidney", "spleen", "testis", "pancreas", "stomach"):
        return head in p
    return False


def score_element(
    element: dict[str, Any], hits: list[dict[str, Any]], symbol_of: dict[str, str]
) -> dict[str, Any]:
    """One element: the eGenes GTEx ties to it against the predicted and the inferred target."""
    genes: dict[str, set[str]] = {}
    for h in hits:
        sym = symbol_of.get(h["gene_id"], h["gene_id"])
        genes.setdefault(sym, set()).add(h["tissue"])
    pred = (element.get("predicted") or {}).get("gene")
    pred_coding = (element.get("predicted_coding") or {}).get("gene")
    inf = (element.get("inferred") or {}).get("gene")
    pred_tissue = (element.get("predicted") or {}).get("tissue")
    tissue_agree = None
    if pred in genes and pred_tissue:
        checks = [tissue_matches(pred_tissue, t) for t in genes[pred]]
        tissue_agree = any(c for c in checks) if any(c is not None for c in checks) else None
    return {
        "id": element["id"],
        "egenes": sorted(genes),
        "n_egenes": len(genes),
        "eqtls": len(hits),
        "predicted": pred,
        "predicted_coding": pred_coding,
        "predicted_in_egenes": (pred in genes) if pred else None,
        "predicted_coding_in_egenes": (pred_coding in genes) if pred_coding else None,
        "inferred": inf,
        "inferred_in_egenes": (inf in genes) if inf else None,
        "predicted_tissue_is_an_eqtl_tissue": tissue_agree,
    }


def summarise(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Over the elements with at least one eQTL: how often each attribution names a measured eGene."""
    with_e = [r for r in rows if r["n_egenes"]]

    def frac(key: str) -> tuple[float | None, int]:
        judged = [r for r in with_e if r[key] is not None]
        return (round(sum(r[key] for r in judged) / len(judged), 3) if judged else None), len(judged)

    p, np_ = frac("predicted_in_egenes")
    pc, npc = frac("predicted_coding_in_egenes")
    i, ni = frac("inferred_in_egenes")
    t, nt = frac("predicted_tissue_is_an_eqtl_tissue")
    both = [r for r in with_e if r["predicted_in_egenes"] is not None and r["inferred_in_egenes"] is not None]
    disagree = [r for r in both if r["predicted"] != r["inferred"]]
    coding = [
        r
        for r in with_e
        if r["predicted_coding_in_egenes"] is not None and r["inferred_in_egenes"] is not None
    ]
    cdis = [r for r in coding if r["predicted_coding"] != r["inferred"]]
    either = (
        round(sum(r["predicted_in_egenes"] or r["inferred_in_egenes"] for r in both) / len(both), 3)
        if both
        else None
    )
    return {
        "elements": len(rows),
        "with_eqtl": len(with_e),
        "fraction_with_eqtl": round(len(with_e) / len(rows), 3) if rows else None,
        "mean_egenes_when_present": round(sum(r["n_egenes"] for r in with_e) / len(with_e), 2)
        if with_e
        else None,
        "predicted_target_is_an_egene": p,
        "predicted_judged": np_,
        "predicted_coding_target_is_an_egene": pc,
        "predicted_coding_judged": npc,
        "inferred_target_is_an_egene": i,
        "inferred_judged": ni,
        "coding_judged": len(coding),
        "when_coding_disagrees": {
            "elements": len(cdis),
            "predicted_right": sum(r["predicted_coding_in_egenes"] for r in cdis),
            "inferred_right": sum(r["inferred_in_egenes"] for r in cdis),
        },
        "both_judged": len(both),
        "either_target_is_an_egene": either,
        "when_they_disagree": {
            "elements": len(disagree),
            "predicted_right": sum(r["predicted_in_egenes"] for r in disagree),
            "inferred_right": sum(r["inferred_in_egenes"] for r in disagree),
            "both_right": sum(r["predicted_in_egenes"] and r["inferred_in_egenes"] for r in disagree),
            "neither_right": sum(not (r["predicted_in_egenes"] or r["inferred_in_egenes"]) for r in disagree),
        },
        "predicted_tissue_is_an_eqtl_tissue": t,
        "tissue_judged": nt,
        "evidence": EVIDENCE,
    }


__all__ = [
    "Intervals",
    "distil",
    "load_hits",
    "score_element",
    "stream_pairs",
    "summarise",
    "tar_members",
    "tissue_matches",
    "tissue_of",
]
