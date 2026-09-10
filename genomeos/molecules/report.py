"""One gene, all its molecules: our transcripts and translated protein, the
UniProt record, and the AlphaFold structure, with agreement checks."""

from __future__ import annotations

from .alphafold import fetch_structure
from .uniprot import uniprot_entry


def _identity(a: str, b: str) -> float:
    n = min(len(a), len(b))
    if not n:
        return 0.0
    return sum(1 for i in range(n) if a[i] == b[i]) / max(len(a), len(b))


def protein_report(
    gene: str, annotation=None, genome=None, structure: bool = True, coords: bool = False
) -> dict:
    out: dict = {"gene": gene}
    ours = None
    if annotation is not None and genome is not None:
        from genomeos.runtime import translate_transcript

        try:
            g = annotation.gene(gene)
            m = annotation.to_module("p")
            txs = [t for t in m.entities[g.id].transcripts if t.cds_segments]
            canon = next((t for t in txs if "Ensembl_canonical" in t.tags), txs[0] if txs else None)
            if canon is not None:
                prot = translate_transcript(genome, canon)
                ours = {
                    "transcript": canon.attrs.get("name", canon.id),
                    "length": len(prot),
                    "sequence": prot,
                    "transcripts": len(g.transcripts),
                    "coding_transcripts": len(txs),
                    "evidence": "GENCODE + GenomeOS translation",
                }
        except KeyError:
            ours = None
    out["ours"] = ours
    entry = uniprot_entry(gene)
    out["uniprot"] = entry.to_dict() if entry else None
    if entry and ours:
        out["agreement"] = {
            "identity": round(_identity(ours["sequence"], entry.sequence), 4),
            "same_length": ours["length"] == entry.length,
        }
    if structure and entry:
        s = fetch_structure(entry.accession)
        out["structure"] = s.to_dict(coords=coords) if s else None
        if s and entry and s.residues:
            out["structure"]["sequence_matches_uniprot"] = "".join(s.residues) == entry.sequence
    return out
