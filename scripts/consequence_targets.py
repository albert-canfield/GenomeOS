# SPDX-License-Identifier: AGPL-3.0-or-later
"""Consequence: GWAS Catalog lead variants and ClinVar pathogenic non-coding variants on the elements.

    uv run python scripts/consequence_targets.py [--distil]

Elements are the ones every committed result carries a prediction for (uniform, constrained, VISTA).
Saves consequence_targets.json. `--distil` re-streams the GWAS Catalog (or when no hits exist).
"""

from __future__ import annotations

import argparse
import glob
import json
import sys

from genomeos.attribution import gwas
from genomeos.attribution.eqtl import Intervals
from genomeos.genome import clinvar
from genomeos.jobs import heartbeat
from genomeos.results import save_result

SETS = {
    "uniform": "data/results/enhancer_targets_chr*.json",
    "constrained": "data/results/constrained_targets_chr*.json",
    "vista": "data/results/vista_chr*.json",
}
CODING = (
    "missense",
    "nonsense",
    "frameshift",
    "synonymous",
    "inframe",
    "stop_lost",
    "start_lost",
    "initiator",
)


def elements() -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    for name, pattern in SETS.items():
        for f in sorted(glob.glob(pattern)):
            with open(f) as fh:
                d = json.load(fh)
            rows = d.get("rows") if isinstance(d.get("rows"), list) else d.get("elements") or []
            for e in rows:
                if "start" in e:
                    out.setdefault(name, []).append(
                        {**e, "chrom": d["chrom"], "key": f"{d['chrom']}:{e['start']}-{e['end']}"}
                    )
    return out


def clinvar_noncoding(sets: dict[str, list[dict]]) -> dict[str, dict]:
    """ClinVar pathogenic variants inside the elements whose consequence is not a coding change."""
    iv = Intervals()
    for rows in sets.values():
        for e in rows:
            iv.add(e["chrom"], e["start"], e["end"], e["key"])
    iv.freeze()
    chroms = {e["chrom"] for rows in sets.values() for e in rows}
    hits: dict[str, list[dict]] = {}
    for chrom in sorted(chroms):
        for (pos, ref, alt), row in clinvar.load_chromosome(chrom).items():
            cons = row.get("consequence", "")
            if any(c in cons for c in CODING) or not cons:
                continue
            for eid in iv.at(chrom, pos):
                hits.setdefault(eid, []).append({"pos": pos, "ref": ref, "alt": alt, **row})
    out = {}
    for name, rows in sets.items():
        with_hit = [e for e in rows if e["key"] in hits]
        pred_ok = [
            e
            for e in with_hit
            if (e.get("predicted") or {}).get("gene") in {h["gene"] for h in hits[e["key"]]}
        ]
        inf_ok = [
            e
            for e in with_hit
            if (e.get("inferred") or {}).get("gene") in {h["gene"] for h in hits[e["key"]]}
        ]
        out[name] = {
            "elements": len(rows),
            "with_pathogenic_noncoding_variant": len(with_hit),
            "variants": sum(len(hits[e["key"]]) for e in with_hit),
            "clinvar_gene_is_the_predicted_target": len(pred_ok),
            "clinvar_gene_is_the_inferred_target": len(inf_ok),
            "examples": [
                {
                    "element": e.get("id"),
                    "key": e["key"],
                    "predicted": (e.get("predicted") or {}).get("gene"),
                    "inferred": (e.get("inferred") or {}).get("gene"),
                    "clinvar": [
                        {
                            k: h[k]
                            for k in ("pos", "gene", "significance", "conditions", "consequence", "stars")
                        }
                        for h in hits[e["key"]][:3]
                    ],
                }
                for e in with_hit[:40]
            ],
        }
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--distil", action="store_true")
    ap.add_argument("--source", default=gwas.GWAS_ZIP, help="URL or local path of the catalog zip")
    args = ap.parse_args()
    sets = elements()
    if args.distil or not (gwas.KNOWLEDGE / "hits.tsv").exists():
        iv = gwas.index([e for rows in sets.values() for e in rows])

        def progress(msg: str) -> None:
            heartbeat("consequence_targets")
            print(msg, flush=True)

        print(json.dumps(gwas.distil(iv, source=args.source, progress=progress), indent=1))
    heartbeat("consequence_targets")
    hits = gwas.load_hits()
    out: dict = {
        "gwas": {},
        "clinvar": clinvar_noncoding(sets),
        "evidence": {"gwas": gwas.EVIDENCE, "clinvar": clinvar.EVIDENCE},
    }
    for name, rows in sets.items():
        s = gwas.summarise(rows, hits)
        out["gwas"][name] = s
        bp = s["by_partition"]
        print(
            f"gwas {name}: {s['with_lead_variant']} of {s['elements']} hold a lead variant "
            f"({s['fraction_with_lead_variant']}, shifted control {s['shifted_control_fraction']}, "
            f"enrichment {s['enrichment_over_shifted']}); named {bp['named_target']['with_lead_variant']} "
            f"vs no target {bp['no_target']['with_lead_variant']}; "
            f"strong {bp['strong_effect']['with_lead_variant']}; "
            f"constrained {bp['constrained']['with_lead_variant']} "
            f"vs {bp['not_constrained']['with_lead_variant']}"
        )
    for name, s in out["clinvar"].items():
        print(
            f"clinvar {name}: {s['with_pathogenic_noncoding_variant']} elements hold a pathogenic non-coding "
            f"variant ({s['variants']}); gene = predicted {s['clinvar_gene_is_the_predicted_target']}, "
            f"= inferred {s['clinvar_gene_is_the_inferred_target']}"
        )
    save_result("consequence_targets", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
