#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The search for an individual-level panel of the person's declared group (ROADMAP area D).

What this opens: each candidate panel's OWN sample metadata (a panel file, a sample sheet, an EGA
dataset record) -- who is in the panel and which population the panel says it sampled.

What this does NOT open: not one genotype of any panel person, and no polygenic score on any panel.
Each candidate is read against the gate frozen in `genomeos/genome/pgs_panel.py` (committed
`167be45`) before these counts were taken. If no candidate meets the gate for the person's declared
group, this inventory is the lane's result: nothing is registered, no score is computed, and
HG002's UNKNOWN for PGS000297 (height, 3,290 SNPs) and PGS000018 (metaGRS) stands.

A percentile is a position among people of the stated group on one score; it is not a risk.
"""

from __future__ import annotations

import csv
import gzip
import json
import re
import sys
import time
import urllib.request
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.genome import pgs_panel as pp  # noqa: E402
from genomeos.genome.pgs_reference import GNOMAD, GNOMAD_LABEL  # noqa: E402
from genomeos.results import save_result  # noqa: E402

NAME = "pgs_panel_inventory"
ENTRY = "scripts/pgs_panel_inventory.py"
OWN_CODE = (
    "genomeos/genome/pgs_panel.py",
    "scripts/pgs_panel_inventory.py",
    "tests/test_pgs_panel.py",
)
CACHE = Path("data/knowledge/pgs_panel")
USER_AGENT = "GenomeOS/0.1"

# What each candidate's metadata is cached as, and the extra files the EBC releases need.
FILES: dict[str, dict[str, str]] = {
    "1kgp3": {"1kgp3.panel": "metadata"},
    "hgdp": {"hgdp_metadata.txt": "metadata"},
    "sgdp": {"sgdp_metadata.txt": "metadata"},
    "hprc2": {"hprc2_sample_metadata.csv": "metadata"},
    "gnomad_hgdp_1kg": {"gnomad_hgdp_1kg_sample_meta.tsv.bgz": "metadata"},
    "aadr_ho": {"aadr_v66.p1_HO.PUB.ind": "metadata"},
    "ebc_jew2010": {
        "ebc_jew2010_samples.xlsx": "https://evolbio.ut.ee/jew/Jew_paper_samples.xlsx",
        "ebc_jew2010.fam": "https://evolbio.ut.ee/jew/jew_paper_data_dbSNP-b131_pos-b37_1KG_strand.fam",
    },
    "ebc_khazar2013": {
        "ebc_khazar2013_full_panel.fam": "extracted from "
        "https://evolbio.ut.ee/khazar/full_data_panel_for_khazar_paper.tar.bz2",
        "ebc_khazar2013_samples.ind": "https://evolbio.ut.ee/khazar/samples.ind",
    },
    "tagc": {"tagc_ega_datasets.json": "metadata"},
    "pgp": {},
    "giab": {},
    "controlled_cohorts": {},
}


def fetch(url: str, dest: Path, log=None) -> Path:
    """Cache one metadata file. Nothing here is a genotype file."""
    if dest.exists():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as out:
        out.write(r.read())
    print(f"fetched {url} -> {dest}", file=log)
    return dest


def read_1kgp3(path: Path) -> list[dict[str, Any]]:
    """The phase 3 panel file: the 2,504 people IGSR publishes as the unrelated set."""
    out = []
    with open(path) as f:
        next(f)
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                out.append({"id": parts[0], "label": parts[1]})
    return out


def read_hgdp(path: Path) -> list[dict[str, Any]]:
    with open(path) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return [{"id": r["sample"], "label": r["population"]} for r in rows if r.get("sample")]


def read_sgdp(path: Path) -> tuple[list[dict[str, Any]], int]:
    """Only the fully public samples; the signed-letter ones fail the access gate by name."""
    with open(path, encoding="utf-8", errors="replace") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    cat = [c for c in (rows[0] if rows else {}) if c.startswith("SGDP-lite category")]
    key = cat[0] if cat else None
    public = [r for r in rows if key and (r.get(key) or "").strip() == "X"]
    people = [{"id": r["Sample_ID"], "label": r["Population_ID"]} for r in public]
    return people, len(rows) - len(public)


def read_hprc2(path: Path) -> list[dict[str, Any]]:
    """HPRC release 2: a child of a trio whose parents are also in the release counts as related."""
    with open(path) as f:
        rows = list(csv.DictReader(f))
    ids = {r["sample_id"] for r in rows}
    out = []
    for r in rows:
        related = bool(r.get("paternal_id") in ids or r.get("maternal_id") in ids)
        out.append(
            {
                "id": r["sample_id"],
                "label": r.get("population_descriptor") or r.get("population_abbreviation") or "",
                "related": related,
            }
        )
    return out


def read_gnomad_subset(path: Path) -> list[dict[str, Any]]:
    """gnomAD's HGDP+1KG subset, by the PROJECTS' own population label (hgdp_tgp_meta).

    gnomAD's `gnomad_population_inference` is an inferred posterior, never a label here."""
    out = []
    with gzip.open(path, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        i_s, i_meta = hdr.index("s"), hdr.index("hgdp_tgp_meta")
        i_rel = hdr.index("relatedness_inference")
        for line in f:
            v = line.rstrip("\n").split("\t")
            try:
                meta = json.loads(v[i_meta])
                rel = json.loads(v[i_rel])
            except (ValueError, IndexError):
                continue
            out.append(
                {
                    "id": v[i_s],
                    "label": meta.get("population") or "",
                    "related": bool(rel.get("related")),
                }
            )
    return out


AADR_DROP = re.compile(r"(QCremove|Ignore|_dup|-o$)", re.IGNORECASE)


def read_aadr(path: Path) -> tuple[list[dict[str, Any]], int]:
    """The .ind file's population column; AADR's own flagged labels are dropped, not re-labelled."""
    people, dropped = [], 0
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            label = parts[2]
            if AADR_DROP.search(label):
                dropped += 1
                continue
            people.append({"id": parts[0], "label": label})
    return people, dropped


def read_ebc_jew2010(xlsx: Path) -> list[dict[str, Any]]:
    """Behar et al. 2010's own sample sheet: `Sample ID` and `Population`."""
    import openpyxl  # a script-only dependency (pyproject's dev extra)

    wb = openpyxl.load_workbook(xlsx, read_only=True)
    rows = list(wb.worksheets[0].iter_rows(values_only=True))
    hdr = [str(h or "") for h in rows[0]]
    i_id, i_pop = hdr.index("Sample ID"), hdr.index("Population")
    return [{"id": str(r[i_id]), "label": str(r[i_pop] or "")} for r in rows[1:] if r and r[i_id] is not None]


def read_ebc_khazar2013(fam: Path, sheet_people: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Behar et al. 2013's released panel carries no population column: the .fam's family and
    individual columns are the sample id twice, and the paper keeps the labels in its supplement.
    The only machine-readable labels for its people are the 2010 sheet's, for the samples the two
    releases share; everyone else is recorded with an empty label and so counts for no group."""
    ids = []
    with open(fam) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                ids.append(parts[1])
    by_id = {p["id"]: p["label"] for p in sheet_people}
    return [{"id": i, "label": by_id.get(i, "")} for i in ids]


def candidate_rows(log=None) -> list[dict[str, Any]]:
    """Every candidate with the people its own metadata names, cached under CACHE."""
    by_key = {c["key"]: c for c in pp.CANDIDATES}
    rows: list[dict[str, Any]] = []

    def path_of(key: str, fname: str) -> Path:
        c, spec = by_key[key], FILES[key][fname]
        url = c["metadata"] if spec == "metadata" else spec
        dest = CACHE / fname
        if not dest.exists() and url and url.startswith("http"):
            fetch(url, dest, log=log)
        return dest

    p = path_of("1kgp3", "1kgp3.panel")
    rows.append(
        {
            "key": "1kgp3",
            "people": read_1kgp3(p),
            "files": [p],
            "note": "the 2,504 people of the phase 3 unrelated set, by the project's own population code",
        }
    )

    p = path_of("hgdp", "hgdp_metadata.txt")
    rows.append(
        {
            "key": "hgdp",
            "people": read_hgdp(p),
            "files": [p],
            "note": "this release publishes no relatedness flag, so per-group counts are upper "
            "bounds on the unrelated people",
        }
    )

    p = path_of("sgdp", "sgdp_metadata.txt")
    people, held = read_sgdp(p)
    rows.append(
        {
            "key": "sgdp",
            "people": people,
            "files": [p],
            "note": f"the fully public samples only; {held} rows are signed-letter or do-not-use "
            "and fail the access gate by name",
        }
    )

    p = path_of("hprc2", "hprc2_sample_metadata.csv")
    rows.append(
        {
            "key": "hprc2",
            "people": read_hprc2(p),
            "files": [p],
            "note": "a release-2 child whose parents are also in the release is counted as related",
        }
    )

    p = path_of("gnomad_hgdp_1kg", "gnomad_hgdp_1kg_sample_meta.tsv.bgz")
    rows.append(
        {
            "key": "gnomad_hgdp_1kg",
            "people": read_gnomad_subset(p),
            "files": [p],
            "note": "labels are the HGDP and 1000 Genomes projects' own; gnomAD's inferred group "
            "posterior is not read as a label",
        }
    )

    p = path_of("aadr_ho", "aadr_v66.p1_HO.PUB.ind")
    people, dropped = read_aadr(p)
    rows.append(
        {
            "key": "aadr_ho",
            "people": people,
            "files": [p],
            "note": f"{dropped} rows carry AADR's own QCremove, Ignore or outlier marks and are "
            "dropped; this file does not separate present-day from ancient people, so the count is "
            "an upper bound on present-day people",
        }
    )

    xlsx = path_of("ebc_jew2010", "ebc_jew2010_samples.xlsx")
    fam2010 = path_of("ebc_jew2010", "ebc_jew2010.fam")
    sheet = read_ebc_jew2010(xlsx)
    with open(fam2010) as f:
        genotyped = {line.split()[1] for line in f if line.split()}
    people2010 = [p for p in sheet if p["id"] in genotyped]
    rows.append(
        {
            "key": "ebc_jew2010",
            "people": people2010,
            "files": [xlsx, fam2010],
            "note": f"{len(sheet)} people in the sheet, {len(people2010)} of them in the released "
            "genotype set",
        }
    )

    fam = CACHE / "ebc_khazar2013_full_panel.fam"
    if not fam.exists():
        raise SystemExit(
            "data/knowledge/pgs_panel/ebc_khazar2013_full_panel.fam is missing: extract "
            "khazar3_1-22_maf0.01mind0.035geno0.005.fam from "
            "https://evolbio.ut.ee/khazar/full_data_panel_for_khazar_paper.tar.bz2"
        )
    ind = path_of("ebc_khazar2013", "ebc_khazar2013_samples.ind")
    people2013 = read_ebc_khazar2013(fam, sheet)
    rows.append(
        {
            "key": "ebc_khazar2013",
            "people": people2013,
            "files": [fam, ind],
            "note": "the release carries no population column; only the people it shares with the "
            "2010 release have a machine-readable label, and the rest count for no group",
        }
    )

    rows.append(
        {
            "key": "giab",
            "people": [
                {"id": "HG002", "label": "Ashkenazi Jewish"},
                {"id": "HG003", "label": "Ashkenazi Jewish"},
                {"id": "HG004", "label": "Ashkenazi Jewish"},
            ],
            "files": [],
            "note": "the person and the person's two parents: by the gate's unrelated rule the "
            "person's own panel counts none of them",
        }
    )

    p = path_of("tagc", "tagc_ega_datasets.json")
    ega = json.loads(p.read_text())
    rows.append(
        {
            "key": "tagc",
            "people": None,
            "declared_size": sum(d.get("num_samples") or 0 for d in ega),
            "files": [p],
            "note": "EGA's own record: "
            + "; ".join(
                f"{d['accession_id']} {d.get('num_samples')} samples, access {d.get('access_type')}, "
                f"policy {d.get('policy_accession_id')}"
                for d in ega
            ),
        }
    )

    rows.append(
        {
            "key": "pgp",
            "people": None,
            "files": [],
            "note": "the participant survey that carries the ancestry answers does not answer: "
            "https://my.pgp-hms.org/google_surveys/1 returns a 301 to itself (checked "
            "2026-10-10), so no count can be taken from it. The project also uploads consumer "
            "array and whole-genome files per participant with no one population label",
        }
    )

    rows.append(
        {
            "key": "controlled_cohorts",
            "people": None,
            "files": [],
            "note": "each needs an access application, which the lane may not make",
        }
    )
    return rows


def readings(group: str, evaluations: list[dict[str, Any]], outcome: str) -> dict[str, Any]:
    """The lane's reading, in the words the gate fixed."""
    best = max(
        (e for e in evaluations if e["group"] == group and e["people"] is not None),
        key=lambda e: e["people"],
        default=None,
    )
    open_best = max(
        (
            e
            for e in evaluations
            if e["group"] == group
            and e["people"] is not None
            and not any(f.startswith("access") for f in e["fails"])
        ),
        key=lambda e: e["people"],
        default=None,
    )
    if outcome == "PANEL_QUALIFIES":
        text = (
            f"A panel qualifies for {GNOMAD_LABEL[group]} ({group}): "
            + ", ".join(
                f"{e['candidate']} {e['people']} people"
                for e in evaluations
                if e["group"] == group and e["meets_size_access_label"]
            )
            + ". The registration of the per-sample score computation follows before any person's "
            "percentile is read."
        )
    else:
        text = (
            f"No public individual-level panel qualifies for {GNOMAD_LABEL[group]} ({group}): of "
            f"{len(pp.CANDIDATES)} candidates examined, the largest count of people carrying that "
            f"group's own label in an openly downloadable release is "
            f"{open_best['people'] if open_best else 0} "
            f"({open_best['candidate'] if open_best else 'none'}), against the gate's minimum of "
            f"{pp.GATE['min_people']} fixed before these counts were taken. The one panel of that "
            "group at the required size is controlled-access (tagc, 128 samples, EGA "
            "EGAD00001000781 under data-access committee EGAC00001000151), which the lane may not "
            "apply for. Nothing is registered, no score is computed on any panel, and the person's "
            "UNKNOWN for "
            + ", ".join(f"{k} ({v.split(':')[0]})" for k, v in pp.LEFT_UNKNOWN.items())
            + " stands. This is a statement about what is published, not about the person."
        )
    return {
        "outcome": outcome,
        "reading": text,
        "largest_count_any_access": None
        if best is None
        else {"candidate": best["candidate"], "people": best["people"]},
        "largest_count_open_access": None
        if open_best is None
        else {"candidate": open_best["candidate"], "people": open_best["people"]},
    }


def result_manifest(read_now: list[str], candidates: list[dict[str, Any]]) -> dict[str, Any]:
    inputs = [
        mf.input_entry(p, partition="not an evaluation partition: the search's inputs")
        for p in sorted(set(read_now))
        if (ROOT / p).exists()
    ]
    sources = [
        {"accession": f"{c['name']} -- sample metadata only", "version": c["terms"]} for c in pp.CANDIDATES
    ]
    sources.append(
        {
            "accession": "lane-pgs's committed polygenic reference result, for the scores left "
            "UNKNOWN and the independence ratios the panel size is set against",
            "version": pp.REFERENCE_RESULT + " (committed 0d19277, registered e9f0c95)",
        }
    )
    return {
        "sources": sources,
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "n/a: no genotype, variant or coordinate is read here; the inventory counts "
        "people and reads access terms",
        "coordinates": {"base": 0, "interval": "half-open", "note": "n/a: no interval is read"},
        "parameters": {
            "person": pp.PERSON,
            "declared_group": pp.declared_group(),
            "min_people": pp.GATE["min_people"],
            "access": pp.GATE["access"],
            "coverage_floor": pp.GATE["coverage_floor"],
            "candidates": [c["key"] for c in pp.CANDIDATES],
            "gate_committed_in": "genomeos/genome/pgs_panel.py (167be45), before these counts",
            "people_read": {
                r["key"]: (None if r["people"] is None else len(r["people"])) for r in candidates
            },
        },
        "exclusions": [
            "every genotype of every panel person: this search reads sample metadata only",
            "every polygenic score on every panel: none is computed here",
            "the person and the person's first-degree relatives, in the person's own panel",
            "a panel person a release's own metadata flags as related, or as a child of parents in "
            "the same release",
            "SGDP's signed-letter and do-not-use samples, and AADR's own QCremove, Ignore and outlier rows",
            "any group assignment by genetic similarity: only a panel's own population label",
        ],
        "partitions": "n/a: a search and a gate, with no outcome measured",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main(argv: list[str] | None = None) -> int:
    log = sys.stdout
    mf.trace_begin()
    started = time.time()
    group = pp.declared_group()
    rows = candidate_rows(log=log)
    by_key = {c["key"]: c for c in pp.CANDIDATES}

    records: list[dict[str, Any]] = []
    evaluations: list[dict[str, Any]] = []
    for r in rows:
        c = by_key[r["key"]]
        counts = None if r["people"] is None else pp.group_counts(r["key"], r["people"])
        n = None if counts is None else counts[group]
        e = pp.evaluate(c, group, n)
        evaluations.append(e)
        records.append(
            {
                "candidate": c["key"],
                "name": c["name"],
                "data": c["data"],
                "access": c["access"],
                "terms": c["terms"],
                "metadata": c["metadata"],
                "count_from": c["count_from"],
                "people_in_metadata": None if r["people"] is None else len(r["people"]),
                "people_counted": None if r["people"] is None else len(pp.counted(r["people"], pp.PERSON)),
                "declared_size": r.get("declared_size"),
                "group_counts": counts,
                "declared_group_people": n,
                "meets_size_access_label": e["meets_size_access_label"],
                "fails": e["fails"],
                "note": r["note"],
            }
        )

    outcome = pp.outcome(evaluations, group)
    other = {}
    for g in GNOMAD:
        if g == group:
            continue
        qualifying = [
            rec["candidate"]
            for rec in records
            if rec["group_counts"]
            and rec["access"] == pp.GATE["access"]
            and rec["group_counts"][g] >= pp.GATE["min_people"]
        ]
        counts = {
            rec["candidate"]: rec["group_counts"][g]
            for rec in records
            if rec["group_counts"] and rec["group_counts"][g]
        }
        if counts:
            other[g] = {"label": GNOMAD_LABEL[g], "counts": counts, "qualifying": qualifying}

    record: dict[str, Any] = {
        "question": "is there a public individual-level panel of the person's declared group for "
        "the polygenic scores lane-pgs left UNKNOWN?",
        "person": pp.PERSON,
        "declared_group": group,
        "declared_group_label": GNOMAD_LABEL[group],
        "declared_group_source": "the ancestry the GIAB donors declared (lane-pgs's DECLARED); "
        "never inferred, and never replaced by a group a panel does cover",
        "gate": pp.GATE,
        "scores_left_unknown": pp.LEFT_UNKNOWN,
        "scores_for_cross_check": pp.CROSS_CHECK,
        "candidates": records,
        "other_groups": other,
        "note": "a percentile is a position among people of the stated group on one score; it is "
        "not a risk. The other_groups counts are descriptive: no score is computed for anyone",
        "seconds": round(time.time() - started, 1),
    }
    record.update(readings(group, evaluations, outcome))
    record["result_manifest"] = result_manifest(mf.traced_reads(), rows)
    path = save_result(NAME, record)
    print(f"wrote {path}", file=log)
    print(f"outcome: {outcome}", file=log)
    print(record["reading"], file=log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
