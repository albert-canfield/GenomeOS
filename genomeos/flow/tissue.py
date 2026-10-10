# SPDX-License-Identifier: AGPL-3.0-or-later
"""The trace on the transcript a tissue makes, beside the canonical one.

`trace_gene` follows a gene's canonical transcript. A tissue may make another isoform: APP695, not
APP770, in the cerebellum. A variant's consequence predicted on the canonical transcript can then
describe a protein that tissue does not make: an exon the tissue skips, another reading frame, another
residue number. `trace_gene_in_tissue` traces the tissue's dominant isoform (GTEx v8 median transcript
TPM, read from the local cache only) and `TissueTrace.substitute` reports the consequence on it with
the canonical consequence beside it.

The tissue's dominant isoform is ESTABLISHED only when all of these hold; otherwise the call is
UNKNOWN with the reason, and the trace never falls back to the canonical transcript:

- the gene has GTEx isoform data in the local cache (`no_isoform_data`),
- the tissue is one GTEx measured for the gene (`tissue_not_measured`; names match exactly, or after
  case and punctuation are folded, never by a partial match),
- the top isoform's median TPM in the tissue is at least MIN_TPM (`low_expression`),
- the next isoform's median TPM is below TIE_RATIO times the top's (`tie`),
- the top isoform is in the local GENCODE models (`dominant_not_in_local_models`) and codes a protein
  there (`dominant_non_coding`: the tissue's most abundant transcript of the gene makes no protein in
  the models, so which protein the tissue makes is not established by this rule).

What a call is and is not: GTEx's median across donors, quantified on GENCODE v26 models and traced on
the local GENCODE v50 model of the same transcript id. It is not any one person's own expression, and a
consequence on either transcript is a prediction of the local trace, not a measurement.
"""

from __future__ import annotations

import json
import re
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from genomeos.flow.trace import CentralDogmaTrace, canonical_transcript, trace
from genomeos.ir.model import Gene, Transcript

#: The dominant isoform's own GTEx median TPM in the tissue must reach this. 1 TPM is the threshold
#: the project's GTEx readers already call "expressed" (molecules/rna.py, summarise_tissues).
MIN_TPM = 1.0
#: The next isoform at or above this fraction of the top's median TPM makes the call a tie.
TIE_RATIO = 0.9

ESTABLISHED = "ESTABLISHED"
UNKNOWN = "UNKNOWN"
REASONS = (
    "no_isoform_data",
    "tissue_not_measured",
    "low_expression",
    "tie",
    "dominant_not_in_local_models",
    "dominant_non_coding",
)
#: How a consequence on the tissue's isoform compares with the canonical one. The first two after
#: `same` are the disagreement the census counts; `change_differs` (same class and residue number,
#: another amino acid or codon) is counted beside it.
AGREEMENTS = ("same", "class_changed", "residue_changed", "change_differs")

ISOFORM_CACHE = Path("data/knowledge/expression")
ISOFORM_EVIDENCE = "measured: GTEx v8 median transcript TPM per tissue (RSEM isoform quantification)"


def cached_isoforms(symbol: str, cache_dir: Path = ISOFORM_CACHE) -> dict[str, Any] | None:
    """A gene's GTEx isoform table from the local cache (`genomeos rna isoforms` writes it), or None.
    Never fetches."""
    p = Path(cache_dir) / f"gtex_isoforms_{symbol.upper()}.json"
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def _fold(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def measured_tissues(isoforms: dict[str, Any] | None) -> list[str]:
    """Every tissue the gene's isoform table measures, sorted."""
    iso = (isoforms or {}).get("isoforms") or {}
    return sorted({t for v in iso.values() for t in (v.get("tissues") or {})})


def match_tissue(tissue: str, available: Iterable[str]) -> str | None:
    """The measured tissue a name means: exact, else the one name equal after folding case and
    punctuation. None when nothing matches or the folded name is ambiguous."""
    names = list(available)
    if tissue in names:
        return tissue
    hits = [n for n in names if _fold(n) == _fold(tissue)]
    return hits[0] if len(hits) == 1 else None


def tissue_isoform(
    isoforms: dict[str, Any] | None, tissue: str, min_tpm: float = MIN_TPM, tie_ratio: float = TIE_RATIO
) -> dict[str, Any]:
    """The tissue's dominant isoform from a cached GTEx isoform table, ESTABLISHED or UNKNOWN with
    the reason. Ranks every measured transcript by its median TPM in the tissue (equal TPMs are
    ordered by transcript id so the call is deterministic, and are a tie by the rule anyway)."""
    out: dict[str, Any] = {
        "tissue": tissue,
        "status": UNKNOWN,
        "reason": None,
        "transcript": None,
        "tpm": None,
        "share": None,
        "next_transcript": None,
        "next_tpm": None,
        "rule": {"min_tpm": min_tpm, "tie_ratio": tie_ratio},
        "evidence": ISOFORM_EVIDENCE,
    }
    iso = (isoforms or {}).get("isoforms") or {}
    if not iso:
        out["reason"] = "no_isoform_data"
        return out
    name = match_tissue(tissue, measured_tissues(isoforms))
    if name is None:
        out["reason"] = "tissue_not_measured"
        return out
    out["tissue"] = name
    ranked = sorted(
        ((float((v.get("tissues") or {}).get(name, 0.0)), tid) for tid, v in iso.items()),
        key=lambda x: (-x[0], x[1]),
    )
    top_tpm, top = ranked[0]
    next_tpm, nxt = ranked[1] if len(ranked) > 1 else (0.0, None)
    total = sum(x[0] for x in ranked)
    out.update(
        {
            "transcript": top,
            "tpm": round(top_tpm, 3),
            "share": round(top_tpm / total, 3) if total > 0 else None,
            "next_transcript": nxt,
            "next_tpm": round(next_tpm, 3),
        }
    )
    if top_tpm < min_tpm:
        out["reason"] = "low_expression"
    elif nxt is not None and next_tpm >= tie_ratio * top_tpm:
        out["reason"] = "tie"
    else:
        out["status"] = ESTABLISHED
    return out


def _unversioned(tid: str) -> str:
    return tid.split(".")[0]


def resolve_in_models(call: dict[str, Any], transcripts: list[Transcript]) -> tuple[dict, Transcript | None]:
    """An ESTABLISHED call held against the local models: the transcript it names, or the call turned
    UNKNOWN because the models lack that transcript or it codes no protein there."""
    if call["status"] != ESTABLISHED:
        return call, None
    tx = {_unversioned(t.id): t for t in transcripts}.get(_unversioned(call["transcript"]))
    if tx is None:
        return {**call, "status": UNKNOWN, "reason": "dominant_not_in_local_models"}, None
    if not tx.cds_segments:
        return {**call, "status": UNKNOWN, "reason": "dominant_non_coding"}, None
    return call, tx


def agreement(canonical: dict[str, Any], isoform: dict[str, Any]) -> str:
    """How the consequence on the tissue's isoform compares with the canonical one (AGREEMENTS)."""
    if canonical.get("consequence") != isoform.get("consequence"):
        return "class_changed"
    if canonical.get("residue") != isoform.get("residue"):
        return "residue_changed"
    if canonical.get("hgvs_p") != isoform.get("hgvs_p") or canonical.get("hgvs_c") != isoform.get("hgvs_c"):
        return "change_differs"
    return "same"


def _brief(sub: dict[str, Any]) -> dict[str, Any]:
    return {k: sub.get(k) for k in ("consequence", "hgvs_p", "hgvs_c", "residue", "region")}


@dataclass(slots=True)
class TissueTrace:
    """A gene traced on the tissue's dominant isoform (when established) and on the canonical one."""

    gene: str
    call: dict[str, Any]
    canonical: CentralDogmaTrace
    canonical_id: str
    isoform: CentralDogmaTrace | None = None
    isoform_id: str | None = None

    @property
    def status(self) -> str:
        return ESTABLISHED if self.isoform is not None else UNKNOWN

    @property
    def same_as_canonical(self) -> bool:
        return self.isoform_id is not None and _unversioned(self.isoform_id) == _unversioned(
            self.canonical_id
        )

    def substitute(self, genomic_pos: int, ref: str, alt: str) -> dict[str, Any]:
        """A single-base change on the tissue's isoform, with the canonical consequence beside it.
        UNKNOWN (and the reason) when the tissue's isoform is not established; never the canonical
        consequence in its place."""
        canon = self.canonical.substitute(genomic_pos, ref, alt)
        out: dict[str, Any] = {
            "tissue": self.call["tissue"],
            "status": self.status,
            "reason": self.call["reason"],
            "transcript": self.isoform.transcript if self.isoform is not None else None,
            "transcript_id": self.isoform_id,
            "canonical_transcript": self.canonical.transcript,
            "canonical": _brief(canon),
        }
        if self.isoform is None:
            out.update({"consequence": UNKNOWN, "hgvs_p": None, "residue": None, "agreement": None})
            return out
        sub = self.isoform.substitute(genomic_pos, ref, alt)
        out.update(
            {
                "consequence": sub.get("consequence"),
                "hgvs_p": sub.get("hgvs_p"),
                "residue": sub.get("residue"),
                "isoform": sub,
                "agreement": agreement(canon, sub),
            }
        )
        return out

    def to_dict(self) -> dict[str, Any]:
        return {
            "gene": self.gene,
            "tissue": self.call["tissue"],
            "status": self.status,
            "reason": self.call["reason"],
            "call": self.call,
            "canonical_transcript": self.canonical.transcript,
            "canonical_transcript_id": self.canonical_id,
            "canonical_protein_length": len(self.canonical.protein),
            "transcript": self.isoform.transcript if self.isoform is not None else None,
            "transcript_id": self.isoform_id,
            "protein_length": len(self.isoform.protein) if self.isoform is not None else None,
            "same_as_canonical": self.same_as_canonical,
            "note": "the isoform is GTEx's population median for the tissue, not one person's own "
            "expression; consequences are predicted by the local trace",
        }


def trace_gene_in_tissue(
    genome,
    gene: Gene,
    transcripts: list[Transcript],
    tissue: str,
    isoforms: dict[str, Any] | None,
    min_tpm: float = MIN_TPM,
    tie_ratio: float = TIE_RATIO,
) -> TissueTrace | None:
    """Trace a gene on the tissue's dominant isoform beside its canonical transcript. `isoforms` is the
    gene's GTEx isoform table (`cached_isoforms`); None makes the tissue's call UNKNOWN. Returns None
    only when the gene has no transcript at all, as `trace_gene` does."""
    canon_tx = canonical_transcript(transcripts)
    if canon_tx is None:
        return None
    canon = trace(genome, canon_tx, gene.symbol)
    call, tx = resolve_in_models(tissue_isoform(isoforms, tissue, min_tpm, tie_ratio), transcripts)
    if tx is None:
        return TissueTrace(gene.symbol, call, canon, canon_tx.id)
    iso = canon if _unversioned(tx.id) == _unversioned(canon_tx.id) else trace(genome, tx, gene.symbol)
    return TissueTrace(gene.symbol, call, canon, canon_tx.id, iso, tx.id)


# ---- the disagreement census over one person's coding SNVs ---------------------------------------------


def gene_census(
    genome,
    gene: Gene,
    transcripts: list[Transcript],
    isoforms: dict[str, Any] | None,
    variants: Iterable[tuple[int, str, str]],
    tissues: list[str],
    min_tpm: float = MIN_TPM,
    tie_ratio: float = TIE_RATIO,
) -> dict[str, Any] | None:
    """One gene's part of the census. `variants` are (0-based position, ref, alt) single-base changes.

    The gene's universe is every variant in the CDS of the canonical transcript or of the dominant
    isoform of any tissue in `tissues` where that isoform is established. For every tissue each
    universe variant is UNKNOWN (with the reason) or compared (AGREEMENTS). Returns counts only, never
    a variant. None when the gene has no coding canonical transcript."""
    canon_tx = canonical_transcript(transcripts)
    if canon_tx is None or not canon_tx.cds_segments:
        return None
    canon_key = _unversioned(canon_tx.id)
    traces: dict[str, CentralDogmaTrace] = {canon_key: trace(genome, canon_tx, gene.symbol)}
    calls: dict[str, tuple[dict[str, Any], str | None]] = {}
    for t in tissues:
        call, tx = resolve_in_models(tissue_isoform(isoforms, t, min_tpm, tie_ratio), transcripts)
        key = _unversioned(tx.id) if tx is not None else None
        if tx is not None and key not in traces:
            traces[key] = trace(genome, tx, gene.symbol)
        calls[t] = (call, key)
    subs: list[dict[str, dict[str, Any]]] = []
    for pos, ref, alt in variants:
        per = {k: tr.substitute(pos, ref, alt) for k, tr in traces.items()}
        if any(s.get("region") == "CDS" for s in per.values()):
            subs.append(per)
    per_tissue: dict[str, dict[str, Any]] = {}
    for t, (call, key) in calls.items():
        c = _empty_counts()
        c["variants"] = len(subs)
        for per in subs:
            if key is None:
                c["unknown"][call["reason"]] += 1
                continue
            c["established"] += 1
            if key == canon_key:
                c["same_transcript"] += 1
                c["same"] += 1
                continue
            c["other_transcript"] += 1
            c[agreement(per[canon_key], per[key])] += 1
        c["disagree"] = c["class_changed"] + c["residue_changed"]
        c["call"] = {
            "status": ESTABLISHED if key else UNKNOWN,
            "reason": call["reason"],
            "dominant_is_canonical": (key == canon_key) if key else None,
        }
        per_tissue[t] = c
    return {
        "gene": gene.symbol,
        "variants": len(subs),
        "coding_on_canonical": sum(1 for per in subs if per[canon_key].get("region") == "CDS"),
        "transcripts_traced": len(traces),
        "per_tissue": per_tissue,
    }


def _empty_counts() -> dict[str, Any]:
    return {
        "variants": 0,
        "unknown": {r: 0 for r in REASONS},
        "established": 0,
        "same_transcript": 0,
        "other_transcript": 0,
        **{a: 0 for a in AGREEMENTS},
        "disagree": 0,
    }


def combine(genes: list[dict[str, Any]], tissues: list[str]) -> dict[str, dict[str, Any]]:
    """The census over genes: per tissue, the sums of `gene_census`'s counts, the genes whose
    established isoform there is not the canonical one, and the genes with a disagreement there."""
    per_tissue: dict[str, dict[str, Any]] = {}
    for t in tissues:
        c = _empty_counts()
        c["genes_established_other_transcript"] = 0
        c["genes_disagreeing"] = 0
        for g in genes:
            x = g["per_tissue"][t]
            for k, v in x.items():
                if k == "unknown":
                    for r, n in v.items():
                        c["unknown"][r] += n
                elif k != "call":
                    c[k] += v
            if x["call"]["dominant_is_canonical"] is False:
                c["genes_established_other_transcript"] += 1
            if x["disagree"]:
                c["genes_disagreeing"] += 1
        per_tissue[t] = c
    return per_tissue
