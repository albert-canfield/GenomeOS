"""Verify the central dogma engine at chromosome scale.

For every protein-coding gene with local models and a compiled UniProt
definition, translate the canonical transcript from the reference sequence
and compare it with the reviewed UniProt sequence. Identity 100% on the
same length means GenomeOS reads the gene exactly as the curators do; a
different length usually means UniProt's canonical isoform is another
transcript; anything else is worth a look. No network: models, sequence
and definitions are all local.
"""

from __future__ import annotations

import json
from typing import Any

from genomeos.molecules.compiler import CACHE
from genomeos.runtime.central_dogma import coding_sequence, table_for, translate_cds


def _identity(a: str, b: str) -> float:
    n = max(len(a), len(b))
    return sum(1 for x, y in zip(a, b, strict=False) if x == y) / n if n else 0.0


def _similarity(a: str, b: str) -> float:
    """Gapped similarity (difflib ratio): an N-terminal offset or one extra exon still scores high,
    a wrong protein does not."""
    from difflib import SequenceMatcher

    return SequenceMatcher(None, a, b, autojunk=False).ratio() if a and b else 0.0


def verify_chromosome(chrom: str, annotation, genome) -> dict[str, Any]:
    module = annotation.to_module("verify")
    rows: list[dict[str, Any]] = []
    for g in annotation.protein_coding():
        if g.locus.chrom != chrom or g.symbol.startswith("ENSG"):
            continue
        p = CACHE / f"{g.symbol}.json"
        if not p.exists():
            continue
        ident = json.loads(p.read_text())["sections"].get("identity", {}).get("items") or {}
        uni = ident.get("sequence")
        if not uni:
            continue
        txs = [t for t in module.entities[g.id].transcripts if t.cds_segments]
        if not txs:
            continue
        canon = next((t for t in txs if "Ensembl_canonical" in t.tags), txs[0])
        ours = translate_cds(coding_sequence(genome, canon), table_for(canon, chrom))
        # best over all coding isoforms, to tell "other isoform" from "real disagreement"
        best = max(
            (_identity(translate_cds(coding_sequence(genome, t), table_for(t, chrom)), uni) for t in txs),
            default=0.0,
        )
        sim = (
            max(
                (
                    _similarity(translate_cds(coding_sequence(genome, t), table_for(t, chrom)), uni)
                    for t in txs
                ),
                default=0.0,
            )
            if best < 0.9999
            else 1.0
        )
        cds_len = sum(seg.length for seg in canon.cds_segments) - canon.cds_phase
        incomplete = any(t.endswith("_NF") for t in canon.tags)
        rows.append(
            {
                "similarity_best_isoform": round(sim, 4),
                # mechanism signatures for a reference that does not encode the curated protein
                "cds_out_of_frame": cds_len % 3 != 0 and not incomplete,
                "translated_fraction_of_cds": round(len(ours) * 3 / cds_len, 3) if cds_len else None,
                "gene": g.symbol,
                "accession": ident.get("accession"),
                "transcript": canon.attrs.get("name", canon.id),
                "ours_aa": len(ours),
                "uniprot_aa": len(uni),
                "identity_canonical": round(_identity(ours, uni), 4),
                "identity_best_isoform": round(best, 4),
                "exact": ours == uni,
                "same_length": len(ours) == len(uni),
            }
        )
    n = len(rows)
    exact = sum(1 for r in rows if r["exact"])
    best_exact = sum(1 for r in rows if r["identity_best_isoform"] >= 0.9999)
    other_isoform = sum(1 for r in rows if not r["exact"] and r["identity_best_isoform"] >= 0.9999)
    offset = [r for r in rows if r["identity_best_isoform"] < 0.98 and r["similarity_best_isoform"] >= 0.9]
    disagree = [r for r in rows if r["similarity_best_isoform"] < 0.9]
    return {
        "chrom": chrom,
        "genes_checked": n,
        "exact_canonical": exact,
        "exact_canonical_fraction": round(exact / n, 4) if n else None,
        "exact_some_isoform": best_exact,
        "exact_some_isoform_fraction": round(best_exact / n, 4) if n else None,
        "canonical_differs_but_another_isoform_matches": other_isoform,
        "same_protein_different_boundaries": len(offset),
        "disagreements": sorted(disagree, key=lambda r: r["similarity_best_isoform"])[:40],
        "disagreement_count": len(disagree),
        # the frameshift signature over every verified gene, not only the disagreements: a gene that
        # fires it and still translates exactly (the mitochondrial genes ending mid-codon) shows the
        # signature detects a property of the annotation, not a fault of the engine
        "out_of_frame_canonical": [
            {"gene": r["gene"], "exact_some_isoform": r["identity_best_isoform"] >= 0.9999}
            for r in rows
            if r["cds_out_of_frame"]
        ],
        "evidence": "derived: GenomeOS translation of GENCODE models on the reference vs UniProt/Swiss-Prot",
        "note": "a length difference with a matching isoform is a canonical-choice difference; a high gapped "
        "similarity with low identity is the same protein with different start or exon boundaries",
    }
