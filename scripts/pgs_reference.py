# SPDX-License-Identifier: AGPL-3.0-or-later
"""Reference distributions for the default polygenic scores, and what LD does to them on the HPRC panel.

    uv run python scripts/pgs_reference.py --register   # freeze RULES before the panel is read
    uv run python scripts/pgs_reference.py --fetch      # allele frequencies from Ensembl, cached
    uv run python scripts/pgs_reference.py              # data/results/pgs_reference.json

The result holds group-level numbers only: per score and population the coverage, mean and sd of
the modelled distribution, the panel's variance ratio per HPRC group (R1) and the panel's z summary
against gnomAD (R2). No person's genotype or score is written to it.
"""

from __future__ import annotations

import argparse
import gzip
import json
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomeos.genome import pgs_reference as ref  # noqa: E402
from genomeos.genome.polygenic import DEFAULT_SCORES, metadata, scoring_file, weight_table  # noqa: E402


def _head() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout.strip()


def _say(msg: str) -> None:
    print(msg, flush=True)


def register(path: Path) -> None:
    if path.exists():
        raise SystemExit(f"{path} exists: a registration is written once")
    out = {
        "what": "pre-registration of the readings of the polygenic reference distributions (ROADMAP §5 "
        "item 13, area D), written before the HPRC panel is read for any score",
        "date": time.strftime("%Y-%m-%d"),
        "code_commit": _head(),
        "module": "genomeos/genome/pgs_reference.py",
        "scores": list(DEFAULT_SCORES),
        "model": "mean = sum 2 p w; variance = sum 2 p (1 - p) w^2 (Hardy-Weinberg, independence); "
        "percentile = 100 * Phi((S - mean) / sd) over the variants with a group frequency and a "
        "person genotype",
        "frequencies": "gnomAD v4 genome groups (afr, ami, amr, asj, eas, fin, mid, nfe, sas) and 1000 "
        "Genomes phase 3 superpopulations, through Ensembl REST; scores above "
        f"{ref.MAX_FETCH_VARIANTS:,} variants are not looked up",
        "panel": "HPRC release 1 in the UCSC hprc90way Cactus alignment, as distilled under "
        "data/knowledge/human_panel (88 haplotype assemblies of 44 people); groups are the HPRC "
        "superpopulation labels",
        "rules": ref.RULES,
        "not_claimed": "no accuracy of any person's percentile is claimed: the percentile is "
        "descriptive, a position in a modelled distribution, and is not compared with any measured "
        "trait",
    }
    path.write_text(json.dumps(out, indent=1) + "\n")
    _say(f"wrote {path}")


def fetch() -> None:
    for s in DEFAULT_SCORES:
        t = ref.fetch_frequencies(s, progress=_say)
        _say(f"{s}: " + ("not fetched (too many variants)" if t is None else f"{len(t['variants']):,} rows"))


def build(out_path: Path, panel: bool = True) -> dict:
    t0 = time.time()
    people = None
    scores = {}
    for s in DEFAULT_SCORES:
        meta = metadata(s)
        with gzip.open(scoring_file(s), "rt") as fh:
            _, rows = weight_table(fh)
        table = ref.load_frequencies(s)
        entry: dict = {
            "name": meta.get("name"),
            "trait": meta.get("trait"),
            "publication": meta.get("publication"),
            "score_ancestry": meta.get("ancestry"),
            "variants_in_score": len(rows),
            "snv_rows": sum(1 for r in rows if len(r[2]) == 1 and len(r[3]) == 1),
        }
        if table is None:
            entry["frequencies"] = {
                "status": "not fetched",
                "reason": f"{len(rows):,} variants exceed the per-variant lookup cap of "
                f"{ref.MAX_FETCH_VARIANTS:,}; a dense score also breaks the independence variance",
            }
        elif not ref.aligned(table, rows):
            raise SystemExit(f"{s}: cached frequency table does not match the scoring file")
        else:
            entry["frequencies"] = {
                "status": "fetched",
                "date": table["date"],
                "requests": table["requests"],
                "bytes": table["bytes"],
                "ids_found_by_position": table["ids_found_by_position"],
                "lookup": dict(Counter(v[7] for v in table["variants"]).most_common()),
            }
            entry["groups"] = {p: ref.group_distribution(table, p) for p in ref.POPULATIONS}
            entry["author_check"] = ref.author_check(table)
        if panel:
            if people is None:
                assemblies = json.loads((ref.PANEL / "chr21" / "meta.json").read_text())["assemblies"]
                people = ref.hprc_people(assemblies)
            chk = ref.panel_check(rows, ref.panel_dosages(rows, people=people, progress=_say), people, table)
            entry["panel"] = {
                "snv_rows_read": chk["snv_rows"],
                "ld": {
                    "groups": {
                        g: {k: v for k, v in o.items() if k != "calibration"}
                        for g, o in chk["groups"].items()
                    },
                    "pooled": chk["pooled"],
                },
                "calibration": {g: o["calibration"] for g, o in chk["groups"].items() if "calibration" in o},
            }
        scores[s] = entry
        _say(f"{s} done ({time.time() - t0:.0f} s)")
    reg = None
    if ref.REGISTRATION.exists():
        reg = json.loads(ref.REGISTRATION.read_text()).get("code_commit")
    out = {
        "what": "reference distributions of the default polygenic scores from public allele "
        "frequencies, and the HPRC panel's check of the independence variance (R1) and of the "
        "gnomAD model (R2)",
        "date": time.strftime("%Y-%m-%d"),
        "code_commit": _head(),
        "registration": {"path": str(ref.REGISTRATION), "code_commit": reg},
        "sources": {
            "frequencies": {
                "service": "Ensembl REST, https://rest.ensembl.org (POST /variation/homo_sapiens?pops=1)",
                "gnomad": "gnomAD v4 genomes per genetic-ancestry group, CC0 1.0 "
                "(https://gnomad.broadinstitute.org/terms)",
                "thousand_genomes": "1000 Genomes phase 3 superpopulations, open access (IGSR)",
            },
            "panel": "HPRC release 1 (Liao et al. 2023) through UCSC hprc90way, distilled under "
            "data/knowledge/human_panel; sample names from NCBI Datasets, populations from "
            + ref.HPRC_METADATA,
        },
        "panel_people": None,
        "rules": ref.RULES,
        "scores": scores,
        "note": ref.NOTE,
        "seconds": round(time.time() - t0, 1),
    }
    if people:
        per_person = {p["sample"]: p["superpopulation"] or "unassigned" for p in people.values()}
        out["panel_people"] = dict(Counter(per_person.values()).most_common())
    out_path.write_text(json.dumps(out, indent=1) + "\n")
    _say(f"wrote {out_path}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--register", action="store_true")
    ap.add_argument("--fetch", action="store_true")
    ap.add_argument("--no-panel", action="store_true")
    ap.add_argument("--out", default=str(ref.RESULT))
    a = ap.parse_args()
    if a.register:
        register(ref.REGISTRATION)
        return 0
    if a.fetch:
        fetch()
        return 0
    build(Path(a.out), panel=not a.no_panel)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
