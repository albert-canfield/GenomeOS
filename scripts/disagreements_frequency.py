# SPDX-License-Identifier: AGPL-3.0-or-later
"""The residue-level disagreements settled by population allele frequency.

For the genes whose canonical translation has the curated protein's length but not its sequence, every
residue that differs names a codon; the single-base change of that codon that would give UniProt's
residue is a variant with a genomic position and an allele. If that allele is the common one in gnomAD,
the reference carries a minor allele and UniProt describes the population; if it is rare or unknown,
the gene model and the UniProt entry describe different things. One VEP request per site, cached by
the lookup layer.

    uv run python scripts/disagreements_frequency.py
"""

from __future__ import annotations

import glob
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from verify_genome_wide import triage  # noqa: E402

from genomeos.flow import trace_gene  # noqa: E402
from genomeos.flow.trace import COMPLEMENT  # noqa: E402
from genomeos.genome import Annotation, IndexedGenome, default_gencode  # noqa: E402
from genomeos.genome.lookup import _vep  # noqa: E402
from genomeos.results import save_result  # noqa: E402
from genomeos.runtime.central_dogma import STANDARD_CODE  # noqa: E402

CACHE = Path("data/knowledge/proteins")
MAX_SITES = 6


def curated_alleles(tr, uni: str) -> list[dict]:
    """For each differing residue: the single-base change giving UniProt's residue, in genomic terms."""
    out = []
    ours = tr.protein
    n = min(len(ours), len(uni))
    for i in range(n):
        if ours[i] == uni[i]:
            continue
        codon = tr.cds[i * 3 : i * 3 + 3]
        positions = tr.genomic_of_residue(i + 1)
        if len(codon) != 3 or len(positions) != 3:
            continue
        found = None
        for k in range(3):
            for base in "ACGU":
                if base == codon[k]:
                    continue
                trial = codon[:k] + base + codon[k + 1 :]
                if STANDARD_CODE.get(trial) == uni[i]:
                    found = (k, base)
                    break
            if found:
                break
        if not found:
            out.append({"residue": i + 1, "ours": ours[i], "uniprot": uni[i], "single_base": False})
            continue
        k, base = found
        ref_rna, alt_rna = codon[k], base
        ref_dna, alt_dna = ref_rna.replace("U", "T"), alt_rna.replace("U", "T")
        if tr.strand.value == "-":
            ref_dna, alt_dna = COMPLEMENT[ref_dna], COMPLEMENT[alt_dna]
        out.append(
            {
                "residue": i + 1,
                "ours": ours[i],
                "uniprot": uni[i],
                "single_base": True,
                "pos": positions[k] + 1,
                "ref": ref_dna,
                "alt": alt_dna,
            }
        )
        if len(out) >= MAX_SITES:
            break
    return out


def main() -> None:
    rows = []
    for f in sorted(glob.glob("data/results/translation_vs_uniprot_chr*.json")):
        with open(f) as fh:
            d = json.load(fh)
        for r in d["disagreements"]:
            r = {**r, "chrom": d["chrom"]}
            r["triage"] = triage(r)
            if r["triage"].startswith("same length"):
                rows.append(r)
    by_chrom: dict[str, list[dict]] = {}
    for r in rows:
        by_chrom.setdefault(r["chrom"], []).append(r)
    out = []
    t0 = time.time()
    for chrom, rs in sorted(by_chrom.items()):
        gff = default_gencode({chrom})
        fa = Path(f"data/reference/{chrom}.fa")
        if not gff or not fa.exists():
            continue
        ann = Annotation.from_gff3(gff, {chrom})
        genome = IndexedGenome(str(fa))
        module = ann.to_module("freq")
        try:
            for r in rs:
                try:
                    g = ann.gene(r["gene"])
                except KeyError:
                    continue
                p = CACHE / f"{r['gene']}.json"
                defn = json.loads(p.read_text()) if p.exists() else {}
                uni = ((defn.get("sections", {}).get("identity") or {}).get("items") or {}).get(
                    "sequence"
                ) or ""
                tr = trace_gene(genome, g, module.entities[g.id].transcripts)
                if tr is None or not uni:
                    continue
                sites = curated_alleles(tr, uni)
                for s in sites:
                    if not s["single_base"]:
                        continue
                    try:
                        v = _vep(chrom, s["pos"], s["ref"], s["alt"])
                    except Exception as ex:  # noqa: BLE001 - VEP hiccups; record and go on
                        v = {"consequence": f"error: {str(ex)[:40]}"}
                    s["gnomad_af"] = v.get("gnomad_af")
                    s["dbsnp"] = (v.get("dbsnp") or [None])[0]
                    s["vep"] = v.get("consequence")
                afs = [s["gnomad_af"] for s in sites if s.get("gnomad_af") is not None]
                known = sum(1 for s in sites if s.get("dbsnp"))
                if not sites:
                    verdict = "no residue-level difference found (length or alignment)"
                elif afs and min(afs) >= 0.5:
                    verdict = "the curated allele is the common one: the reference carries a minor allele"
                elif afs and max(afs) >= 0.01:
                    verdict = "the curated allele is a known polymorphism, not the majority"
                elif known:
                    verdict = "the curated allele is in dbSNP but rare or unmeasured"
                else:
                    verdict = (
                        "the curated allele is unknown to dbSNP and gnomAD: gene model or entry difference"
                    )
                out.append(
                    {
                        "gene": r["gene"],
                        "chrom": chrom,
                        "accession": r["accession"],
                        "differing_residues": sum(1 for a, b in zip(tr.protein, uni, strict=False) if a != b),
                        "sites": sites,
                        "verdict": verdict,
                    }
                )
                print(
                    f"{r['gene']:12} {len(sites)} sites, AF {[round(a, 3) for a in afs][:6]}, "
                    f"dbSNP {known}: {verdict[:70]}",
                    flush=True,
                )
        finally:
            genome.close()
    counts: dict[str, int] = {}
    for x in out:
        counts[x["verdict"]] = counts.get(x["verdict"], 0) + 1
    save_result(
        "translation_disagreements_frequency",
        {
            "genes": len(out),
            "verdicts": dict(sorted(counts.items(), key=lambda kv: -kv[1])),
            "rows": out,
            "seconds": round(time.time() - t0),
            "evidence": "derived: the single-base change giving the curated residue; measured: gnomAD allele "
            "frequency and dbSNP via Ensembl VEP",
        },
    )
    print("verdicts:", json.dumps(counts, indent=1), flush=True)


if __name__ == "__main__":
    main()
