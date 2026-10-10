# SPDX-License-Identifier: AGPL-3.0-or-later
"""Reference distributions for the default polygenic scores, and what LD does to them on the HPRC panel.

    uv run python scripts/pgs_reference.py --register   # freeze RULES before the panel is read
    uv run python scripts/pgs_reference.py --fetch      # allele frequencies from Ensembl, cached
    uv run python scripts/pgs_reference.py              # data/results/pgs_reference.json

The result holds group-level numbers only: per score and population the coverage, mean and sd of
the modelled distribution, the panel's variance ratio per HPRC group (R1) and the panel's z summary
against gnomAD (R2). No person's genotype or score is written to it. Both files go through
`genomeos.results.save_result`, with every file read under data/ hashed in the manifest.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.genome import pgs_reference as ref  # noqa: E402
from genomeos.genome.polygenic import DEFAULT_SCORES, metadata, scoring_file, weight_table  # noqa: E402
from genomeos.results import save_result  # noqa: E402

OWN_CODE = (
    "genomeos/genome/pgs_reference.py",
    "genomeos/genome/polygenic.py",
    "scripts/pgs_reference.py",
    "tests/test_pgs_reference.py",
)
SOURCES = [
    {
        "accession": "PGS Catalog harmonised scoring files (hmPOS GRCh38) and score metadata for "
        + ", ".join(DEFAULT_SCORES),
        "version": "as fetched on 2026-09-12 into data/knowledge/pgs, digested in inputs; EBI terms of use",
    },
    {
        "accession": "Ensembl REST POST /variation/homo_sapiens?pops=1: gnomAD v4 genome groups "
        "(gnomADg:*, CC0 1.0) and 1000 Genomes phase 3 superpopulations",
        "version": "per-score tables fetched on the date each records, cached as "
        "data/knowledge/pgs/<id>_af.json.gz and digested in inputs",
    },
    {
        "accession": "HPRC release 1 in the UCSC hprc90way Cactus alignment (88 haplotype assemblies "
        "of 44 people, Liao et al. 2023), distilled by genomeos/attribution/human_panel.py",
        "version": "the distilled store under data/knowledge/human_panel, digested in inputs",
    },
    {
        "accession": "panel sample names from NCBI Datasets; populations from " + ref.HPRC_METADATA,
        "version": "cached as data/knowledge/pgs/hprc_samples.json, digested in inputs",
    },
]

PARTITIONS = {
    "panel groups": "groups of at least 10 people are pooled (R1) and judged (R2); smaller groups are "
    "reported only"
}


def _say(msg: str) -> None:
    print(msg, flush=True)


def _inputs() -> list[dict]:
    """Every file read under data/ so far, hashed: the panel store and the reference as two groups."""
    reads = sorted(set(mf.traced_reads()))
    panel = [r for r in reads if r.startswith("data/knowledge/human_panel/")]
    genome = [r for r in reads if r.startswith("data/reference/")]
    out = [
        mf.input_entry(r, partition=None)
        for r in reads
        if r not in panel and r not in genome and Path(r).is_file()
    ]
    if panel:
        out.append(
            mf.files_entry(
                "HPRC panel store (meta, sites, blocks per chromosome)", panel, partition="reference panel"
            )
        )
    if genome:
        out.append(mf.files_entry("hg38 reference sequence and its indexes", genome, partition=None))
    return out


def _manifest(parameters: dict, exclusions: list[str], partitions: dict) -> dict:
    entries = _inputs()
    return {
        "sources": SOURCES,
        "inputs": entries,
        "input_count": len(entries),
        "assembly": "GRCh38/hg38: PGS Catalog hmPOS positions, Ensembl GRCh38 mappings and the "
        "hg38-referenced HPRC alignment",
        "coordinates": {"base": 1, "interval": "closed"},
        "parameters": parameters,
        "exclusions": exclusions,
        "partitions": partitions,
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def _rows(s: str) -> list[tuple]:
    with gzip.open(scoring_file(s), "rt") as fh:
        return weight_table(fh)[1]


def _people() -> dict:
    assemblies = json.loads((ref.PANEL / "chr21" / "meta.json").read_text())["assemblies"]
    return ref.hprc_people(assemblies)


def _per_group(people: dict) -> dict[str, int]:
    per_person = {p["sample"]: p["superpopulation"] or "unassigned" for p in people.values()}
    return dict(Counter(per_person.values()).most_common())


def register() -> Path:
    if (ref.REGISTRATION).exists():
        raise SystemExit(f"{ref.REGISTRATION} exists: a registration is written once")
    variants = {s: len(_rows(s)) for s in DEFAULT_SCORES}
    groups = _per_group(_people())
    payload = {
        "what": "pre-registration of the readings of the polygenic reference distributions (ROADMAP §5 "
        "item 13, area D), written before the HPRC panel is read for any score",
        "module": "genomeos/genome/pgs_reference.py",
        "scores": variants,
        "panel_people": groups,
        "model": "mean = sum 2 p w; variance = sum 2 p (1 - p) w^2 (Hardy-Weinberg, independence); "
        "percentile = 100 * Phi((S - mean) / sd) over the variants with a group frequency and a "
        "person genotype",
        "frequencies": "gnomAD v4 genome groups (" + ", ".join(ref.GNOMAD) + ") and 1000 Genomes "
        f"phase 3 superpopulations, through Ensembl REST; scores above {ref.MAX_FETCH_VARIANTS:,} "
        "variants are not looked up",
        "rules": ref.RULES,
        "not_claimed": "no accuracy of any person's percentile is claimed: the percentile is "
        "descriptive, a position in a modelled distribution, and is not compared with any measured "
        "trait",
    }
    payload["result_manifest"] = _manifest(
        {
            "rules": ref.RULES,
            "max_fetch_variants": ref.MAX_FETCH_VARIANTS,
            "populations": list(ref.POPULATIONS),
        },
        ["the panel's hs1 (T2T-CHM13) row is not a person and is not read"],
        PARTITIONS,
    )
    p = save_result("pgs_reference_registration", payload)
    _say(f"registered: {p}")
    return p


def fetch() -> None:
    for s in DEFAULT_SCORES:
        t = ref.fetch_frequencies(s, progress=_say)
        _say(f"{s}: " + ("not fetched (too many variants)" if t is None else f"{len(t['variants']):,} rows"))


def build(panel: bool = True) -> Path:
    t0 = time.time()
    people = _people() if panel else None
    reg = json.loads(ref.REGISTRATION.read_text()) if ref.REGISTRATION.exists() else {}
    scores = {}
    exclusions = ["the panel's hs1 (T2T-CHM13) row is not a person and is not read"]
    for s in DEFAULT_SCORES:
        meta = metadata(s)
        rows = _rows(s)
        table = ref.load_frequencies(s)
        snv = sum(1 for r in rows if len(r[2]) == 1 and len(r[3]) == 1)
        entry: dict = {
            "name": meta.get("name"),
            "trait": meta.get("trait"),
            "publication": meta.get("publication"),
            "score_ancestry": meta.get("ancestry"),
            "variants_in_score": len(rows),
            "snv_rows": snv,
        }
        if snv < len(rows):
            exclusions.append(f"{s}: {len(rows) - snv:,} non-SNV rows are not read on the panel")
        if table is None:
            entry["frequencies"] = {
                "status": "not fetched",
                "reason": f"{len(rows):,} variants exceed the per-variant lookup cap of "
                f"{ref.MAX_FETCH_VARIANTS:,}; a dense score also breaks the independence variance",
            }
            exclusions.append(f"{s}: no allele frequencies ({len(rows):,} variants above the lookup cap)")
        elif not ref.aligned(table, rows):
            raise SystemExit(f"{s}: cached frequency table does not match the scoring file")
        else:
            lookup = dict(Counter(v[7] for v in table["variants"]).most_common())
            entry["frequencies"] = {
                "status": "fetched",
                "date": table["date"],
                "requests": table["requests"],
                "bytes": table["bytes"],
                "ids_found_by_position": table["ids_found_by_position"],
                "lookup": lookup,
            }
            entry["groups"] = {p: ref.group_distribution(table, p) for p in ref.POPULATIONS}
            entry["author_check"] = ref.author_check(table)
            unmatched = len(rows) - lookup.get("matched", 0)
            if unmatched:
                exclusions.append(f"{s}: {unmatched:,} rows carry no frequency (lookup: {lookup})")
        if panel:
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
    payload = {
        "what": "reference distributions of the default polygenic scores from public allele "
        "frequencies, and the HPRC panel's check of the independence variance (R1) and of the "
        "gnomAD model (R2)",
        "registration": {
            "result": "pgs_reference_registration",
            "code": (reg.get("result_manifest") or {}).get("code"),
        },
        "panel_people": _per_group(people) if people else None,
        "rules": ref.RULES,
        "scores": scores,
        "note": ref.NOTE,
        "seconds": round(time.time() - t0, 1),
    }
    payload["result_manifest"] = _manifest(
        {
            "rules": ref.RULES,
            "max_fetch_variants": ref.MAX_FETCH_VARIANTS,
            "populations": list(ref.POPULATIONS),
        },
        exclusions,
        PARTITIONS,
    )
    p = save_result("pgs_reference", payload)
    _say(f"wrote {p}")
    return p


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--register", action="store_true")
    ap.add_argument("--fetch", action="store_true")
    ap.add_argument("--no-panel", action="store_true")
    a = ap.parse_args()
    if a.register:
        register()
    elif a.fetch:
        fetch()
    else:
        build(panel=not a.no_panel)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
