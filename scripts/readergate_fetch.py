# SPDX-License-Identifier: AGPL-3.0-or-later
"""M2b data: the public ENCODE RNA-seq gene quantifications for the thirteen reader biosamples.

    uv run python scripts/readergate_fetch.py [--dry-run]

Chooses one polyA plus RNA-seq experiment per reader biosample by the rule written in
`genomeos.attribution.readergate.RNA_SELECTION`, from portal METADATA only, downloads that
experiment's GRCh38 V29 gene-quantification TSVs to data/cache/readergate/rna/ (git-ignored), checks
each file's md5 against the portal, and records URL, accession, size, md5, sha256 and licence in
data/results/readergate_rna_sources.json.

No value of any file is read: a file is hashed (its bytes), and its first line is read to check that
the columns the run will use (`gene_id`, `TPM`) are there. The values are opened only by
scripts/readergate_run.py, after the registration is committed. No model request, no key, no money.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import context_evidence as ce  # noqa: E402
from genomeos.attribution import readergate as rg  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "readergate_rna_sources"
ENCODE = "https://www.encodeproject.org"
LICENCE = (
    "ENCODE Consortium data: free to download and use without restriction under the ENCODE data use"
    " policy; cite the ENCODE Consortium and the experiment accession"
    " (https://www.encodeproject.org/help/citing-encode/)"
)
BUDGET_BYTES = 500_000_000
OWN_CODE = (
    "genomeos/attribution/compile.py",
    "genomeos/attribution/readergate.py",
    "scripts/readergate_fetch.py",
    "scripts/readergate_census.py",
    "scripts/readergate_register.py",
    "scripts/readergate_run.py",
    "tests/test_readergate.py",
)


def _get(url: str, timeout: int = 120) -> dict[str, Any]:
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "GenomeOS/0.1"})
    with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310
        return json.load(r)


def experiments(term: str) -> dict[str, dict[str, Any]]:
    """accession -> metadata of the human, unperturbed, released polyA plus experiments on one term."""
    q = (
        f"{ENCODE}/search/?type=Experiment&status=released&perturbed=false"
        f"&assay_title={urllib.parse.quote_plus(rg.RNA_ASSAY)}"
        f"&biosample_ontology.term_name={urllib.parse.quote(term)}"
        "&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens"
        "&format=json&limit=all&field=accession&field=date_released&field=lab.title"
        "&field=biosample_summary&field=audit.ERROR.category"
    )
    return {e["accession"]: e for e in _get(q).get("@graph", [])}


def files(term: str) -> list[dict[str, Any]]:
    q = (
        f"{ENCODE}/search/?type=File&status=released&output_type=gene+quantifications&assembly=GRCh38"
        f"&genome_annotation=V29&file_format=tsv&assay_title={urllib.parse.quote_plus(rg.RNA_ASSAY)}"
        f"&biosample_ontology.term_name={urllib.parse.quote(term)}&format=json&limit=all"
        "&field=accession&field=dataset&field=href&field=md5sum&field=file_size"
        "&field=simple_biosample_summary&field=biological_replicates&field=date_created"
    )
    return _get(q).get("@graph", [])


def _neg_date(d: str) -> str:
    """A key that sorts later dates first (YYYY-MM-DD, digit-complemented)."""
    return "".join(chr(ord("9") - ord(c) + ord("0")) if c.isdigit() else c for c in d)


def choose(term: str) -> dict[str, Any]:
    """The experiment RNA_SELECTION takes for one biosample, its files, and every candidate considered."""
    exps = experiments(term)
    candidates: dict[str, list[dict[str, Any]]] = {}
    excluded: dict[str, str] = {}
    for f in files(term):
        acc = f["dataset"].strip("/").split("/")[-1]
        summary = (f.get("simple_biosample_summary") or "").lower()
        if acc not in exps:
            excluded[f["accession"]] = f"experiment {acc} is not a human unperturbed released one"
        elif exps[acc].get("audit", {}).get("ERROR"):
            excluded[f["accession"]] = f"experiment {acc} carries an ERROR audit"
        elif any(w in summary for w in rg.RNA_EXCLUDE_WORDS):
            excluded[f["accession"]] = f"biosample summary {summary!r}"
        else:
            candidates.setdefault(acc, []).append(f)
    if not candidates:
        raise LookupError(f"no qualifying {rg.RNA_ASSAY} gene quantification for {term!r}")
    # most files, then the latest date_released, then the lowest accession
    best = sorted(
        candidates, key=lambda a: (-len(candidates[a]), _neg_date(exps[a].get("date_released", "")), a)
    )[0]
    e = exps[best]
    return {
        "biosample": term,
        "experiment": best,
        "experiment_url": f"{ENCODE}/experiments/{best}/",
        "date_released": e.get("date_released"),
        "lab": (e.get("lab") or {}).get("title"),
        "biosample_summary": e.get("biosample_summary"),
        "files": sorted(candidates[best], key=lambda f: f["accession"]),
        "candidates": {a: len(v) for a, v in sorted(candidates.items())},
        "files_excluded": excluded,
    }


def download(f: dict[str, Any], out_dir: Path) -> dict[str, Any]:
    """One file to the cache, md5 checked against the portal; its sha256 and its header line."""
    url = ENCODE + f["href"]
    path = out_dir / f"{f['accession']}.tsv"
    if not path.exists() or path.stat().st_size != f.get("file_size"):
        tmp = path.with_suffix(".part")
        req = urllib.request.Request(url, headers={"User-Agent": "GenomeOS/0.1 (readergate)"})
        with urllib.request.urlopen(req, timeout=300) as r, tmp.open("wb") as out:  # noqa: S310
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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true", help="choose and print; download nothing")
    args = ap.parse_args()
    t0 = time.time()
    chosen = [choose(term) for term in ce.READER_TERMS]
    total = sum(f.get("file_size", 0) for c in chosen for f in c["files"])
    for c in chosen:
        print(c["biosample"], c["experiment"], len(c["files"]), c["date_released"], c["candidates"])
    print(f"{sum(len(c['files']) for c in chosen)} files, {total / 1e6:.1f} MB")
    if total > BUDGET_BYTES:
        raise SystemExit(f"refused: {total / 1e6:.0f} MB is over the {BUDGET_BYTES / 1e6:.0f} MB budget")
    if args.dry_run:
        return 0
    rg.RNA_DIR.mkdir(parents=True, exist_ok=True)
    for c in chosen:
        c["files"] = [download(f, rg.RNA_DIR) for f in c["files"]]
        print("fetched", c["biosample"], [f["accession"] for f in c["files"]], flush=True)
    cell_of = dict(zip(ce.READER_TERMS, rg.CELLS, strict=True))
    payload = {
        "lane": rg.LANE,
        "date": time.strftime("%Y-%m-%d"),
        "what": "the measured RNA for M2b: one ENCODE polyA plus RNA-seq experiment per reader biosample",
        "selection_rule": rg.RNA_SELECTION,
        "values_read": "none: each file was hashed and its header line read; no TPM value was opened",
        "licence": LICENCE,
        "bytes_total": sum(f["bytes"] for c in chosen for f in c["files"]),
        "files_total": sum(len(c["files"]) for c in chosen),
        "biosamples": {cell_of[c["biosample"]]: c for c in chosen},
        "alphagenome_requests": 0,
        "money": "none: public files, no key read",
        "seconds": round(time.time() - t0, 1),
    }
    inputs = [mf.input_entry(f["path"], partition=None) for c in chosen for f in c["files"]]
    manifest = {
        "sources": [
            {"accession": c["experiment"], "version": f"ENCODE portal, released {c['date_released']}"}
            for c in chosen
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38, GENCODE V29 gene models (ENCODE4 uniform RNA-seq pipeline)",
        "coordinates": "n/a: gene-level quantifications keyed by gene_id",
        "parameters": {"assay": rg.RNA_ASSAY, "exclude_words": list(rg.RNA_EXCLUDE_WORDS)},
        "exclusions": [
            "files of experiments that are perturbed, non-human, unreleased or carry an ERROR audit",
            "files whose biosample summary names a fraction, treatment, disease or modification",
            "every experiment but the one RNA_SELECTION takes per biosample",
        ],
        "partitions": "n/a: the measured outcome of a registered test; no value read here",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }
    print("saved", save_result(RESULT, payload, manifest=manifest))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
