# SPDX-License-Identifier: AGPL-3.0-or-later
"""Predicted and inferred enhancer targets against GTEx eQTLs.

    uv run python scripts/eqtl_targets.py --distil   # stream GTEx once, keep the pairs inside our elements
    uv run python scripts/eqtl_targets.py            # score every element with a prediction

Saves eqtl_targets.json. Elements carry a margin of MARGIN bases on each side when the pairs are kept,
since an eQTL a few hundred bases from an element is in linkage with it.

Elements are every element a committed result carries a prediction for: the uniform samples
(enhancer_targets_chr*), the constrained tier (constrained_targets_chr*) and VISTA's (vista_chr*).
"""

from __future__ import annotations

import argparse
import glob
import json
import sys

from genomeos.attribution import eqtl
from genomeos.jobs import heartbeat
from genomeos.results import save_result

MARGIN = 500
SETS = {
    "uniform": "data/results/enhancer_targets_chr*.json",
    "constrained": "data/results/constrained_targets_chr*.json",
    "vista": "data/results/vista_chr*.json",
}


def elements() -> dict[str, list[dict]]:
    """Every element with its chromosome, coordinates, prediction and inference, by set."""
    out: dict[str, list[dict]] = {}
    for name, pattern in SETS.items():
        for f in sorted(glob.glob(pattern)):
            with open(f) as fh:
                d = json.load(fh)
            rows = d.get("rows") if isinstance(d.get("rows"), list) else d.get("elements") or []
            for e in rows:
                if "start" not in e:
                    continue
                out.setdefault(name, []).append({**e, "chrom": d["chrom"]})
    return out


def symbol_map(chroms: set[str]) -> dict[str, str]:
    from genomeos.genome import Annotation, default_gencode

    out: dict[str, str] = {}
    for c in sorted(chroms):
        path = default_gencode({c})
        if path is None:
            continue
        ann = Annotation.from_gff3(path, {c})
        out.update({gid.split(".")[0]: g.symbol for gid, g in ann.genes.items()})
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--distil", action="store_true", help="stream the GTEx archive and keep the hits")
    ap.add_argument("--tissues", default="", help="comma list of GTEx tissues (default: all 49)")
    args = ap.parse_args()
    sets = elements()
    if args.distil or not any(eqtl.KNOWLEDGE.glob("hits_*.tsv")):
        iv = eqtl.Intervals()
        seen = set()
        for rows in sets.values():
            for e in rows:
                key = (e["chrom"], e["start"], e["end"])
                if key not in seen:
                    seen.add(key)
                    iv.add(
                        e["chrom"],
                        max(0, e["start"] - MARGIN),
                        e["end"] + MARGIN,
                        f"{e['chrom']}:{e['start']}-{e['end']}",
                    )
        iv.freeze()
        tissues = [t for t in args.tissues.split(",") if t]

        def progress(msg: str) -> None:
            heartbeat("eqtl_distil")
            print(msg, flush=True)

        s = eqtl.distil(iv, tissues=tissues or None, progress=progress)
        print(json.dumps(s, indent=1))
        if args.distil:
            return 0
    heartbeat("eqtl_targets")
    chroms = {e["chrom"] for rows in sets.values() for e in rows}
    symbols = symbol_map(chroms)
    hits = eqtl.load_hits()
    out: dict = {"sets": {}, "margin_bp": MARGIN, "evidence": eqtl.EVIDENCE}
    for name, rows in sets.items():
        scored = []
        for e in rows:
            key = f"{e['chrom']}:{e['start']}-{e['end']}"
            r = eqtl.score_element({**e, "id": key}, hits.get(key, []), symbols)
            r["element"] = e["id"]
            r["chrom"] = e["chrom"]
            scored.append(r)
        out["sets"][name] = {
            "summary": eqtl.summarise(scored),
            "elements_with_eqtl": [r for r in scored if r["n_egenes"]],
        }
        s = out["sets"][name]["summary"]
        print(
            f"{name}: {s['with_eqtl']} of {s['elements']} elements carry an eQTL; "
            f"predicted target is an eGene {s['predicted_target_is_an_egene']} (n={s['predicted_judged']}), "
            f"inferred "
            f"{s['inferred_target_is_an_egene']} (n={s['inferred_judged']}); when they disagree "
            f"{s['when_they_disagree']}"
        )
    save_result("eqtl_targets", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
