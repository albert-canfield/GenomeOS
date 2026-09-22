# SPDX-License-Identifier: AGPL-3.0-or-later
"""The translation disagreements, read with the isoform the body actually makes.

For each gene where the reference's canonical transcript does not translate to the curated UniProt
protein (96 genes, triaged by mechanism in scripts/verify_genome_wide.py), ask GTEx which transcript
dominates across tissues, translate that one locally, and compare it with UniProt. Three answers are
possible and each means something different: the dominant isoform matches UniProt (the annotation's
canonical choice was the disagreement, not the genome); the dominant isoform is the canonical one and
still differs (the reference allele is the disagreement); or no measured isoform matches (a different
product, or a gene GTEx does not resolve).

    uv run python scripts/disagreements_isoforms.py
"""

from __future__ import annotations

import glob
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from verify_genome_wide import triage  # noqa: E402

from genomeos.genome import Annotation, Genome, default_gencode  # noqa: E402
from genomeos.molecules.rna import gtex_isoforms  # noqa: E402
from genomeos.molecules.verify import _identity  # noqa: E402
from genomeos.results import save_result  # noqa: E402
from genomeos.runtime.central_dogma import coding_sequence, table_for, translate_cds  # noqa: E402

CACHE = Path("data/knowledge/proteins")


def main() -> None:
    rows = []
    for f in sorted(glob.glob("data/results/translation_vs_uniprot_chr*.json")):
        with open(f) as fh:
            d = json.load(fh)
        for r in d["disagreements"]:
            rows.append({**r, "chrom": d["chrom"], "triage": triage({**r, "chrom": d["chrom"]})})
    by_chrom: dict[str, list[dict]] = {}
    for r in rows:
        by_chrom.setdefault(r["chrom"], []).append(r)
    out = []
    t0 = time.time()
    for chrom, rs in sorted(by_chrom.items()):
        gff = default_gencode({chrom})
        fa = Path(f"data/reference/{chrom}.fa.gz")
        if not gff or not fa.exists():
            continue
        ann = Annotation.from_gff3(gff, {chrom})
        genome = Genome.from_fasta(fa)
        module = ann.to_module("iso")
        for r in rs:
            try:
                g = ann.gene(r["gene"])
            except KeyError:
                continue
            p = CACHE / f"{r['gene']}.json"
            defn = json.loads(p.read_text()) if p.exists() else {}
            uni = ((defn.get("sections", {}).get("identity") or {}).get("items") or {}).get("sequence") or ""
            txs = [t for t in module.entities[g.id].transcripts if t.cds_segments]
            ident_by_tx = {}
            for t in txs:
                tid = t.id.split(".")[0]
                try:
                    ident_by_tx[tid] = round(
                        _identity(translate_cds(coding_sequence(genome, t), table_for(t, chrom)), uni), 4
                    )
                except Exception:  # noqa: BLE001
                    ident_by_tx[tid] = 0.0
            canonical = next((t.id.split(".")[0] for t in txs if "Ensembl_canonical" in t.tags), None)
            labels = ((defn.get("sections", {}).get("genomic_origin") or {}).get("items") or {}).get(
                "transcripts"
            ) or []
            try:
                iso = gtex_isoforms(r["gene"], transcripts=labels)
            except Exception as ex:  # noqa: BLE001
                iso = {"error": str(ex)[:80]}
            dominant = None
            if iso.get("isoforms"):
                dominant = max(iso["isoforms"], key=lambda k: iso["isoforms"][k]["tissues_dominant"])
                if iso["isoforms"][dominant]["tissues_dominant"] == 0:
                    dominant = None
            dom_ident = ident_by_tx.get(dominant) if dominant else None
            best_tx = max(ident_by_tx, key=ident_by_tx.get) if ident_by_tx else None
            if iso.get("error") or dominant is None:
                verdict = "no expressed isoform in GTEx (or gene not in GTEx)"
            elif dom_ident is not None and dom_ident >= 0.9999:
                verdict = (
                    "the isoform the body makes matches UniProt: the canonical choice was the disagreement"
                )
            elif dominant == canonical:
                verdict = (
                    "the body makes the canonical transcript and it differs: "
                    "the reference allele is the disagreement"
                )
            elif dom_ident is None:
                verdict = "the dominant transcript is non-coding in the local models"
            else:
                verdict = "the body makes another isoform and it differs too"
            out.append(
                {
                    "gene": r["gene"],
                    "chrom": chrom,
                    "accession": r["accession"],
                    "triage": r["triage"],
                    "canonical": canonical,
                    "identity_canonical": r["identity_canonical"],
                    "dominant_transcript": dominant,
                    "dominant_in_tissues": iso["isoforms"][dominant]["tissues_dominant"] if dominant else 0,
                    "identity_dominant": dom_ident,
                    "best_transcript": best_tx,
                    "identity_best": ident_by_tx.get(best_tx) if best_tx else None,
                    "verdict": verdict,
                }
            )
            print(
                f"{r['gene']:12} {r['triage'][:40]:40} dominant {dominant} ident {dom_ident} "
                f"→ {verdict[:60]}",
                flush=True,
            )
    counts: dict[str, int] = {}
    cross: dict[str, dict[str, int]] = {}
    for x in out:
        counts[x["verdict"]] = counts.get(x["verdict"], 0) + 1
        cross.setdefault(x["triage"].split(":")[0], {})[x["verdict"].split(":")[0]] = (
            cross.setdefault(x["triage"].split(":")[0], {}).get(x["verdict"].split(":")[0], 0) + 1
        )
    save_result(
        "translation_disagreements_isoforms",
        {
            "genes": len(out),
            "verdicts": dict(sorted(counts.items(), key=lambda kv: -kv[1])),
            "triage_by_verdict": cross,
            "rows": out,
            "seconds": round(time.time() - t0),
            "evidence": "measured: GTEx v8 transcript TPM; derived: local translation of each transcript "
            "against UniProt",
        },
    )
    print("verdicts:", json.dumps(counts, indent=1), flush=True)


if __name__ == "__main__":
    main()
