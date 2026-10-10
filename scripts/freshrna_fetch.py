# SPDX-License-Identifier: AGPL-3.0-or-later
"""M2i data: ENCODE RNA-seq gene quantifications of the thirteen reader cell types whose VALUES this
project has never opened.

    uv run python scripts/freshrna_fetch.py [--dry-run]

Builds the READ set from the record (git-tracked files, the last 3,000 commit messages and the
untracked measured-RNA track records, minus the accessions that occur only as a REJECTED CANDIDATE of
item 7's own selection: `genomeos.attribution.freshrna.READ_SET_RULE`), then chooses one unread
qualifying experiment per reader biosample term by `freshrna.RNA_SELECTION`, from portal METADATA
only. It proves the non-overlap at accession level -- experiment, gene-quantification file and
REPLICATE LIBRARY -- downloads the chosen TSVs to data/cache/freshrna/rna/ (git-ignored), checks each
file's md5 against the portal, and records everything in data/results/freshrna_sources.json.

No value of any file is read: a file is hashed (its bytes) and its first line is read to check that
the columns the run will use (`gene_id`, `TPM`) are there. The values are opened only by
scripts/freshrna_run.py, after the registration is committed. No model request, no key, no money.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import freshrna as fr  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "freshrna_sources"
ENCODE = "https://www.encodeproject.org"
READERGATE_SOURCES = Path("data/results/readergate_rna_sources.json")
TRACK_RECORDS = "data/knowledge/rna/tracks_*.json"
ACC = re.compile(r"ENC(?:SR|FF|LB|BS)[0-9A-Z]{6}")
LICENCE = (
    "ENCODE Consortium data: free to download and use without restriction under the ENCODE data use"
    " policy; cite the ENCODE Consortium and the experiment accession"
    " (https://www.encodeproject.org/help/citing-encode/)"
)
BUDGET_BYTES = 500_000_000
OWN_CODE = (
    "genomeos/attribution/freshrna.py",
    "scripts/freshrna_fetch.py",
    "scripts/freshrna_register.py",
    "scripts/freshrna_run.py",
    "tests/test_freshrna.py",
)


def _get(url: str, timeout: int = 180) -> dict[str, Any]:
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "GenomeOS/0.1"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310
            return json.load(r)
    except urllib.error.HTTPError as e:
        if e.code == 404:  # the portal answers 404 for a search with no result
            return {"@graph": []}
        raise


# ---- the READ set, from the record -----------------------------------------------------------------
def read_set() -> dict[str, Any]:
    """READ_SET_RULE: every accession the record mentions, minus item 7's rejected candidates."""
    tracked = subprocess.run(  # noqa: S603
        ["git", "grep", "-ohIE", "ENC(SR|FF|LB|BS)[0-9A-Z]{6}"], capture_output=True, text=True, check=False
    ).stdout.split()
    messages = ACC.findall(
        subprocess.run(  # noqa: S603
            ["git", "log", "--format=%B", "-n", "3000"], capture_output=True, text=True, check=False
        ).stdout
    )
    tracks: list[str] = []
    for p in sorted(Path().glob(TRACK_RECORDS)):
        tracks += ACC.findall(p.read_text())
    src = json.loads(READERGATE_SOURCES.read_text())
    metadata_only: set[str] = set()
    for v in src["biosamples"].values():
        metadata_only |= {a for a in v["candidates"] if a != v["experiment"]}
        metadata_only |= set(v["files_excluded"])
    mentioned = set(tracked) | set(messages) | set(tracks)
    return {
        "accessions": sorted(mentioned - metadata_only),
        "mentioned_in_record": len(mentioned),
        "metadata_only_rejected_candidates": sorted(metadata_only),
        "sources_of_the_record": [
            "every ENCSR/ENCFF/ENCLB/ENCBS in any git-tracked file (git grep)",
            "every one in the last 3,000 commit messages (git log --format=%B -n 3000)",
            f"every one in the untracked measured-RNA track records {TRACK_RECORDS}",
        ],
        "readergate_read": {
            "experiments": sorted({v["experiment"] for v in src["biosamples"].values()}),
            "files": sorted({f["accession"] for v in src["biosamples"].values() for f in v["files"]}),
        },
    }


# ---- the portal -------------------------------------------------------------------------------------
def experiments(term: str, assay: str) -> dict[str, dict[str, Any]]:
    """accession -> metadata of the human, unperturbed, released experiments of one assay on one term."""
    q = (
        f"{ENCODE}/search/?type=Experiment&status=released&perturbed=false"
        f"&assay_title={urllib.parse.quote_plus(assay)}"
        f"&biosample_ontology.term_name={urllib.parse.quote(term)}"
        "&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens"
        "&format=json&limit=all&field=accession&field=date_released&field=lab.title"
        "&field=biosample_summary&field=audit.ERROR.category"
    )
    return {e["accession"]: e for e in _get(q).get("@graph", [])}


def files(term: str, assay: str) -> list[dict[str, Any]]:
    q = (
        f"{ENCODE}/search/?type=File&status=released&output_type=gene+quantifications&assembly=GRCh38"
        f"&genome_annotation=V29&file_format=tsv&assay_title={urllib.parse.quote_plus(assay)}"
        f"&biosample_ontology.term_name={urllib.parse.quote(term)}&format=json&limit=all"
        "&field=accession&field=dataset&field=href&field=md5sum&field=file_size"
        "&field=simple_biosample_summary&field=biological_replicates"
    )
    return _get(q).get("@graph", [])


def replicates(accessions: list[str]) -> dict[str, list[dict[str, Any]]]:
    """experiment accession -> its replicates' library and biosample accessions (the overlap proof)."""
    out: dict[str, list[dict[str, Any]]] = {a: [] for a in accessions}
    for i in range(0, len(accessions), 25):
        q = (
            f"{ENCODE}/search/?type=Replicate&format=json&limit=all"
            + "".join(f"&experiment.accession={a}" for a in accessions[i : i + 25])
            + "&field=experiment.accession&field=library.accession"
            "&field=library.biosample.accession&field=biological_replicate_number"
        )
        for r in _get(q).get("@graph", []):
            lib = r.get("library") or {}
            out.setdefault(r["experiment"]["accession"], []).append(
                {
                    "replicate": r.get("biological_replicate_number"),
                    "library": lib.get("accession"),
                    "biosample": (lib.get("biosample") or {}).get("accession"),
                }
            )
    return out


def _neg_date(d: str) -> str:
    """A key that sorts later dates first (YYYY-MM-DD, digit-complemented)."""
    return "".join(chr(ord("9") - ord(c) + ord("0")) if c.isdigit() else c for c in d)


def choose(term: str, read: set[str], read_biosamples: set[str]) -> dict[str, Any]:
    """RNA_SELECTION: the unread qualifying experiment this lane takes for one reader term."""
    rejected: list[dict[str, Any]] = []
    for assay in fr.RNA_ASSAYS:
        exps = experiments(term, assay)
        candidates: dict[str, list[dict[str, Any]]] = {}
        for f in files(term, assay):
            acc = f["dataset"].strip("/").split("/")[-1]
            summary = (f.get("simple_biosample_summary") or "").lower()
            if acc not in exps or exps[acc].get("audit", {}).get("ERROR"):
                continue
            if any(w in summary for w in fr.rg.RNA_EXCLUDE_WORDS):
                continue
            candidates.setdefault(acc, []).append(f)
        reps = replicates(sorted(candidates))
        unread: dict[str, list[dict[str, Any]]] = {}
        for acc, fl in candidates.items():
            seen = sorted(
                {acc}
                | {f["accession"] for f in fl}
                | {r["library"] for r in reps[acc] if r["library"]}
                | {r["biosample"] for r in reps[acc] if r["biosample"]}
            )
            hit = [a for a in seen if a in read]
            if hit:
                rejected.append({"experiment": acc, "assay": assay, "accessions_in_the_read_set": hit})
            else:
                unread[acc] = fl
        if not unread:
            continue

        def shares_biosample(acc: str, reps: dict[str, list[dict[str, Any]]] = reps) -> int:
            return 1 if {r["biosample"] for r in reps[acc] if r["biosample"]} & read_biosamples else 0

        best = sorted(
            unread,
            key=lambda a: (
                shares_biosample(a),
                -len(unread[a]),
                _neg_date(exps[a].get("date_released", "")),
                a,
            ),
        )[0]
        e = exps[best]
        return {
            "biosample": term,
            "assay": assay,
            "experiment": best,
            "experiment_url": f"{ENCODE}/experiments/{best}/",
            "date_released": e.get("date_released"),
            "lab": (e.get("lab") or {}).get("title"),
            "biosample_summary": e.get("biosample_summary"),
            "replicates": reps[best],
            "shares_a_read_biosample": bool(shares_biosample(best)),
            "files": sorted(unread[best], key=lambda f: f["accession"]),
            "other_unread_qualifying": sorted(set(unread) - {best}),
            "rejected_as_read": rejected,
        }
    raise LookupError(f"no unread qualifying RNA-seq gene quantification for {term!r}: {rejected}")


def download(f: dict[str, Any], out_dir: Path) -> dict[str, Any]:
    """One file to the cache, md5 checked against the portal; its sha256 and its header line."""
    url = ENCODE + f["href"]
    path = out_dir / f"{f['accession']}.tsv"
    if not path.exists() or path.stat().st_size != f.get("file_size"):
        tmp = path.with_suffix(".part")
        req = urllib.request.Request(url, headers={"User-Agent": "GenomeOS/0.1 (freshrna)"})
        with urllib.request.urlopen(req, timeout=600) as r, tmp.open("wb") as out:  # noqa: S310
            while chunk := r.read(1 << 20):
                out.write(chunk)
        tmp.rename(path)
    md5 = hashlib.md5()  # noqa: S324 - compared with the portal's published md5, not used for security
    sha = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(1 << 20):
            md5.update(chunk)
            sha.update(chunk)
    if md5.hexdigest() != f["md5sum"]:
        raise ValueError(f"{path}: md5 {md5.hexdigest()} is not the portal's {f['md5sum']}")
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
    missing = [c for c in ("gene_id", "TPM") if c not in header]
    if missing:
        raise ValueError(f"{path}: no {missing} column in the header line")
    return {
        "accession": f["accession"],
        "url": url,
        "path": path.as_posix(),
        "bytes": path.stat().st_size,
        "md5": f["md5sum"],
        "sha256": sha.hexdigest(),
        "biological_replicates": f.get("biological_replicates"),
        "header": header,
    }


def overlap(chosen: list[dict[str, Any]], read: dict[str, Any]) -> dict[str, Any]:
    """NON_OVERLAP, proved at accession level and committed, including the two biosample exceptions."""
    rg_read = read["readergate_read"]
    reps = replicates(sorted({c["experiment"] for c in chosen}))
    read_reps = replicates(rg_read["experiments"])
    mine = {
        "experiments": sorted({c["experiment"] for c in chosen}),
        "files": sorted({f["accession"] for c in chosen for f in c["files"]}),
        "libraries": sorted({r["library"] for v in reps.values() for r in v if r["library"]}),
        "biosamples": sorted({r["biosample"] for v in reps.values() for r in v if r["biosample"]}),
    }
    theirs = {
        "experiments": rg_read["experiments"],
        "files": rg_read["files"],
        "libraries": sorted({r["library"] for v in read_reps.values() for r in v if r["library"]}),
        "biosamples": sorted({r["biosample"] for v in read_reps.values() for r in v if r["biosample"]}),
    }
    inter = {k: sorted(set(mine[k]) & set(theirs[k])) for k in mine}
    return {
        "statement": fr.NON_OVERLAP,
        "read_set_rule": fr.READ_SET_RULE,
        "read_set_size": len(read["accessions"]),
        "this_lane": mine,
        "already_read": theirs,
        "intersection": inter,
        "disjoint_experiments": not inter["experiments"],
        "disjoint_files": not inter["files"],
        "disjoint_libraries": not inter["libraries"],
        "disjoint_biosamples": not inter["biosamples"],
        "cells_sharing_a_read_biosample": sorted(
            c["biosample"] for c in chosen if c["shares_a_read_biosample"]
        ),
        "intersection_with_the_whole_read_set": sorted(
            set(mine["experiments"] + mine["files"] + mine["libraries"]) & set(read["accessions"])
        ),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true", help="choose and print; download nothing")
    args = ap.parse_args()
    t0 = time.time()
    read = read_set()
    readset = set(read["accessions"])
    read_reps = replicates(read["readergate_read"]["experiments"])
    read_biosamples = {r["biosample"] for v in read_reps.values() for r in v if r["biosample"]}
    print(f"READ set: {len(readset)} accessions ({read['mentioned_in_record']} mentioned in the record)")
    chosen = [choose(term, readset, read_biosamples) for term in fr.READER_TERMS]
    total = sum(f.get("file_size", 0) for c in chosen for f in c["files"])
    for c in chosen:
        print(
            f"  {c['biosample']:24s} {c['assay']:19s} {c['experiment']} {len(c['files'])} files"
            f" {c['date_released']} shares_read_biosample={c['shares_a_read_biosample']} | {c['lab']}"
        )
    print(f"{sum(len(c['files']) for c in chosen)} files, {total / 1e6:.1f} MB")
    ov = overlap(chosen, read)
    print("intersection with what the project has read:", ov["intersection"])
    if ov["intersection"]["experiments"] or ov["intersection"]["files"] or ov["intersection"]["libraries"]:
        raise SystemExit(
            "refused: the chosen set shares an experiment, a file or a library with the read set"
        )
    if total > BUDGET_BYTES:
        raise SystemExit(f"refused: {total / 1e6:.0f} MB is over the {BUDGET_BYTES / 1e6:.0f} MB budget")
    if args.dry_run:
        return 0
    fr.RNA_DIR.mkdir(parents=True, exist_ok=True)
    for c in chosen:
        c["files"] = [download(f, fr.RNA_DIR) for f in c["files"]]
        print("fetched", c["biosample"], [f["accession"] for f in c["files"]], flush=True)
    cell_of = dict(zip(fr.READER_TERMS, fr.CELLS, strict=True))
    payload = {
        "lane": fr.LANE,
        "date": time.strftime("%Y-%m-%d"),
        "what": (
            "the measured RNA for M2i: one ENCODE RNA-seq experiment per reader biosample whose values no"
            " file, script, result or document of this project has ever opened"
        ),
        "selection_rule": fr.RNA_SELECTION,
        "what_is_fresh": fr.WHAT_IS_FRESH,
        "no_biosample_added": (
            "READER_TERMS is context_evidence's, unchanged: the thirteen reader cell types are exactly"
            " those of items 7, 7d, 7h and M2h, and no reader biosample is added (docs/ROADMAP.md §5's"
            " stop list, 'No more biosamples')"
        ),
        "values_read": "none: each file was hashed and its header line read; no TPM value was opened",
        "licence": LICENCE,
        "bytes_total": sum(f["bytes"] for c in chosen for f in c["files"]),
        "files_total": sum(len(c["files"]) for c in chosen),
        "assays": {cell_of[c["biosample"]]: c["assay"] for c in chosen},
        "biosamples": {cell_of[c["biosample"]]: c for c in chosen},
        "non_overlap": ov,
        "read_set": read,
        "alphagenome_requests": 0,
        "money": "none: public files, no key read",
        "seconds": round(time.time() - t0, 1),
    }
    inputs = [mf.input_entry(f["path"], partition=None) for c in chosen for f in c["files"]]
    inputs += [mf.input_entry(READERGATE_SOURCES.as_posix(), partition=None)]
    inputs += [mf.input_entry(p.as_posix(), partition=None) for p in sorted(Path().glob(TRACK_RECORDS))]
    manifest = {
        "sources": [
            {"accession": c["experiment"], "version": f"ENCODE portal, released {c['date_released']}"}
            for c in chosen
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38, GENCODE V29 gene models (ENCODE4 uniform RNA-seq pipeline)",
        "coordinates": "n/a: gene-level quantifications keyed by gene_id",
        "parameters": {"assays": list(fr.RNA_ASSAYS), "exclude_words": list(fr.rg.RNA_EXCLUDE_WORDS)},
        "exclusions": [
            "every experiment, file, library or biosample accession in the READ set (READ_SET_RULE)",
            "files of experiments that are perturbed, non-human, unreleased or carry an ERROR audit",
            "files whose biosample summary names a fraction, treatment, disease or modification",
            "every unread experiment but the one RNA_SELECTION takes per biosample",
        ],
        "partitions": "n/a: the measured outcome of a registered test; no value read here",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }
    print("saved", save_result(RESULT, payload, manifest=manifest))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
