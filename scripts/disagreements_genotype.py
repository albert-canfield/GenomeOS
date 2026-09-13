# SPDX-License-Identifier: AGPL-3.0-or-later
"""The translation disagreements held against real genotypes: does a person carry the curated protein?

For each gene whose reference canonical transcript does not translate to the UniProt protein, take the
SNVs a local individual carries inside that transcript's CDS, apply them (alt allele, every site) to the
coding sequence, translate, and compare with UniProt again. If the identity rises to 1.0 the reference
carries a minor allele and this person carries the curated protein; if it does not move, the
disagreement is not an allele the person has. Indels inside the CDS are counted, not applied (a
frameshift allele of the reference would need one to restore the frame).

    uv run python scripts/disagreements_genotype.py [--name HG002] [--name HG003 ...]
"""

from __future__ import annotations

import argparse
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
from genomeos.genome.individuals import rows_in, vcf_path  # noqa: E402
from genomeos.molecules.verify import _identity  # noqa: E402
from genomeos.results import save_result  # noqa: E402
from genomeos.runtime.central_dogma import table_for, translate_cds  # noqa: E402

CACHE = Path("data/knowledge/proteins")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", action="append", default=None, help="local individuals (default HG002)")
    args = ap.parse_args()
    names = args.name or ["HG002"]
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
        fa = Path(f"data/reference/{chrom}.fa")
        if not gff or not fa.exists():
            continue
        ann = Annotation.from_gff3(gff, {chrom})
        genome = IndexedGenome(str(fa))
        module = ann.to_module("geno")
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
                canon = next((t for t in module.entities[g.id].transcripts if t.id == tr.transcript), None)
                table = table_for(canon, chrom) if canon is not None else None
                before = _identity(translate_cds(tr.cds, table) if table else tr.protein, uni)
                for name in names:
                    vcf = vcf_path(name, chrom)
                    if vcf is None:
                        continue
                    cds = list(tr.cds)
                    snvs, indels, applied = 0, 0, []
                    for f in rows_in(vcf, g.locus.start + 1, g.locus.end):
                        ref, alt = f[3], f[4].split(",")[0]
                        where = tr.residue_of(int(f[1]) - 1) or {}
                        if where.get("region") != "CDS":
                            continue
                        if len(ref) != 1 or len(alt) != 1:
                            indels += 1
                            continue
                        snvs += 1
                        a = alt.upper()
                        if tr.strand.value == "-":
                            a = COMPLEMENT.get(a, "N")
                        cds[where["cds"]] = a.replace("T", "U")
                        applied.append(f"{f[1]}{ref}>{alt}({f[9].split(':')[0]})")
                    after = (
                        _identity(translate_cds("".join(cds), table) if table else "", uni)
                        if snvs
                        else before
                    )
                    if snvs == 0 and indels == 0:
                        verdict = "no variant of this person inside the CDS"
                    elif after >= 0.9999:
                        verdict = (
                            "this person's alleles give the curated protein: "
                            "the reference carries a minor allele"
                        )
                    elif after > before + 1e-6:
                        verdict = "this person's SNVs move the translation towards the curated protein"
                    elif indels and not snvs:
                        verdict = "only indels inside the CDS (not applied): a frame change is possible"
                    else:
                        verdict = "this person's variants leave the disagreement as it is"
                    out.append(
                        {
                            "gene": r["gene"],
                            "chrom": chrom,
                            "triage": r["triage"],
                            "individual": name,
                            "cds_snvs": snvs,
                            "cds_indels": indels,
                            "identity_reference": round(before, 4),
                            "identity_with_alleles": round(after, 4),
                            "applied": applied[:8],
                            "verdict": verdict,
                        }
                    )
                    print(
                        f"{name} {r['gene']:12} {r['triage'][:34]:34} snv {snvs:2} indel {indels:2} "
                        f"{before:.3f} → {after:.3f}  {verdict[:60]}",
                        flush=True,
                    )
        finally:
            genome.close()
    counts: dict[str, dict[str, int]] = {}
    for x in out:
        c = counts.setdefault(x["individual"], {})
        c[x["verdict"]] = c.get(x["verdict"], 0) + 1
    save_result(
        "translation_disagreements_genotype",
        {
            "individuals": names,
            "genes": len({x["gene"] for x in out}),
            "verdicts": counts,
            "rows": out,
            "seconds": round(time.time() - t0),
            "evidence": "measured genotypes; derived: local translation of the canonical CDS with the "
            "person's SNVs applied, against UniProt",
        },
    )
    print("verdicts:", json.dumps(counts, indent=1), flush=True)


if __name__ == "__main__":
    main()
