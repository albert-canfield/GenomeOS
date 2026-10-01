# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register amendment 2 to N1 before any measurement (data/results/n1_registration_amendment_2.json).

    uv run --frozen python scripts/n1_register_amendment_2.py

The response label becomes the producers' published Anderson-Darling result (Figshare+ 21632564, CC0),
BH-adjusted, at their documented level, and the question becomes one about assigned perturbations. This
reads structure and identities only, by HTTP range requests, cached in data/cache/n1/: the published
file's header line (its perturbation labels) and the first gene label after it, and the raw pseudobulk's
row and gene identities. No p-value, expression value or other measurement is read. The original
registration and amendment 1 stay as committed.
"""

from __future__ import annotations

import csv
import importlib.util
import io
import json
import sys
import urllib.request
import zlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import n1_perturb_response as n1  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

ORIGINAL = RESULTS_DIR / "n1_registration.json"
AMENDMENT_1 = RESULTS_DIR / f"{n1.AMENDMENT}.json"
CACHE = Path("data/cache/n1")
AD_HEADER = CACHE / "ad_header.json"
RAW_IDENTITY = CACHE / "raw_identity.json"
MODULE = Path("genomeos/attribution/n1_perturb_response.py")
RUNNER = Path("scripts/n1_run_v3.py")
PMC = "https://europepmc.org/article/PMC/PMC9380471"
DEPOSIT = "https://doi.org/10.25452/figshare.plus.21632564.v1"


def _register_module():
    spec = importlib.util.spec_from_file_location("n1_register", Path(__file__).with_name("n1_register.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def ad_header() -> dict:
    """The published file's header line and the first field of the next line, by range requests over the
    gzip prefix; decompression stops there and no value field is parsed. Cached."""
    if AD_HEADER.exists():
        return json.loads(AD_HEADER.read_text())
    d = zlib.decompressobj(16 + zlib.MAX_WBITS)
    buf, fetched, start, step = b"", 0, 0, 1 << 16
    while True:
        req = urllib.request.Request(
            n1.AD_SOURCE["url"],
            headers={"Range": f"bytes={start}-{start + step - 1}", "User-Agent": "GenomeOS/n1-structure"},
        )
        with urllib.request.urlopen(req, timeout=120) as r:  # noqa: S310
            if r.status != 206:
                raise RuntimeError(f"the server did not honour the range request (HTTP {r.status})")
            chunk = r.read()
        fetched += len(chunk)
        start += len(chunk)
        buf += d.decompress(chunk)
        nl = buf.find(b"\n")
        if nl >= 0 and buf.find(b",", nl + 1) >= 0:
            break
    header_line = buf[:nl].decode().rstrip("\r")
    rec = {
        "url": n1.AD_SOURCE["url"],
        "bytes_fetched": fetched,
        "columns": next(csv.reader([header_line]))[1:],
        "first_row_label": n1.first_field(buf[nl + 1 :].decode(errors="replace")),
    }
    CACHE.mkdir(parents=True, exist_ok=True)
    AD_HEADER.write_text(json.dumps(rec))
    return rec


def raw_identity() -> dict:
    """The raw pseudobulk's obs index and var gene_id by range requests, reduced to their digest. Cached."""
    if RAW_IDENTITY.exists():
        return json.loads(RAW_IDENTITY.read_text())
    import h5py

    reg = _register_module()
    f = n1.SOURCE["files"]["raw"]
    rf = reg.RangeFile(f["url"], f["bytes"])
    with h5py.File(io.BufferedReader(rf, buffer_size=reg.BLOCK), "r") as h:
        obs = [x.decode() for x in h["obs/gene_transcript"][()]]
        var = [x.decode() for x in h["var/gene_id"][()]]
    rec = {
        "url": f["url"],
        "bytes_fetched": rf.fetched,
        "datasets_read": ["obs/gene_transcript", "var/gene_id"],
        "identity_digest": n1.identity_digest(obs, var),
    }
    CACHE.mkdir(parents=True, exist_ok=True)
    RAW_IDENTITY.write_text(json.dumps(rec))
    return rec


def resource() -> dict:
    return {
        "chosen": "the producers' supplemental deposit, checked first: Figshare+ 21632564 (Replogle and "
        "Weissman, "
        "2022-11-29, CC0)",
        "link": DEPOSIT,
        "file": n1.AD_SOURCE,
        "deposit_quote": "Commonly requested supplemental files supporting the publication ... 1) All "
        "Anderson-Darling p-values for differential expression analysis",
        "statistic": [
            {
                "where": f"{PMC}, STAR Methods, 'Gene-level differential expression testing using the "
                "Anderson-Darling and Mann-Whitney tests'",
                "quote": "we z -normalize gene expression relative to control cells as described ... and "
                "for each "
                "gene test whether the distribution of normalized expression is identical between control "
                "cells "
                "bearing non-targeting sgRNAs and cells bearing each perturbation. We used two tests "
                "implemented "
                "in scipy: the Anderson-Darling test (scipy.stats.anderson_ksamp), which is broadly "
                "sensitive to "
                "changes in distribution",
            },
            {
                "where": "the same section",
                "quote": "For the Anderson-Darling test we extended the range of p values beyond those "
                "available "
                "in scipy's implementation by computing the p-value for many values of the test statistic "
                "using "
                "R's kSamples package and interpolating any intermediate values using "
                "scipy.interpolate.interp1d.",
            },
        ],
        "raw_or_adjusted": {
            "reading": "adjusted, Benjamini-Hochberg: the file is named 'BH-corrected' and the methods say "
            "so",
            "quote": "p -values in both cases were adjusted for multiple hypothesis testing using the "
            "Benjamini-Hochberg procedure to produce the final results.",
            "not_documented": "the family over which the adjustment ran (per perturbation, or the whole "
            "matrix)",
        },
        "documented_level": {
            "level": n1.AD_LEVEL,
            "quote": "at least 50 differentially expressed genes at a significance of p < 0.05 by "
            "Anderson-Darling "
            "test following Benjamini-Hochberg correction",
            "where": f"{PMC}, STAR Methods, 'Functional analyses of strong perturbations'",
        },
        "controls_as_published": "By the AD test, 2,935 of 9,608 genetic perturbations targeting a primary "
        "transcript (30.5%) compared with 12 of 585 controls (2.1%) caused >10 DEGs in K562 cells. (Results)",
        "principal_transcript": {
            "rule": "a factor with several rows is evaluated by its principal-transcript row (P1P2, else P1)",
            "quote": "targeted the principal ''P1'' transcript identified by the FANTOM consortium. (A "
            "handful of "
            "genes also had perturbations targeting the P2 transcript that did not generally have effects.)",
            "where": f"{PMC}, STAR Methods, 'Leverage scores'",
        },
        "missing_values": "not documented, and not inspected: the rule treats an empty, non-numeric or "
        "non-finite "
        "field, or one outside [0, 1], as missing, and a gene absent from the rows or listed twice as "
        "missing for "
        "every factor",
        "not_used": "Nadig et al. 2025 (Figshare 29498366) publishes unadjusted DESeq2 Wald p-values by gene "
        "symbol, with the combination of a gene's guides undocumented; the producers' own deposit was "
        "checked "
        "first and is used",
    }


def main() -> int:
    original = json.loads(ORIGINAL.read_text())
    plan = original["plan"]
    head = ad_header()
    rawid = raw_identity()
    cols = head["columns"]
    colset = set(cols)
    principal = {c["factor"]: n1.principal_row(c) for c in plan["candidates"]}
    ident = json.loads(Path("data/cache/n1/K562_gwps_normalized_bulk_01.identities.json").read_text())
    coverage = {
        "perturbation_columns": len(cols),
        "columns_equal_pseudobulk_rows_as_a_set": colset == set(ident["obs_index"]),
        "columns_in_pseudobulk_order": cols == ident["obs_index"],
        "non_targeting_columns": sum(1 for c in cols if n1.NON_TARGETING in c),
        "candidate_factors": len(plan["candidates"]),
        "candidate_principal_columns_present": sum(1 for v in principal.values() if v in colset),
        "candidate_rows_present": sum(
            1 for c in plan["candidates"] for r in c["rows"] if r["label"] in colset
        ),
        "multi_row_candidates_and_chosen_row": {
            c["factor"]: principal[c["factor"]] for c in plan["candidates"] if len(c["rows"]) > 1
        },
        "gene_axis": "Ensembl gene IDs as row labels",
        "first_row_label": head["first_row_label"],
        "first_row_label_index_in_pseudobulk_genes": ident["var_gene_id"].index(head["first_row_label"])
        if head["first_row_label"] in ident["var_gene_id"]
        else None,
        "universe_genes_in_file": "not countable before the authorised download: the gene labels sit inside "
        "the "
        "compressed values. The run's labels pass counts them first, parsing no value; an absent gene stays "
        "missing",
        "raw_pseudobulk_identity_matches_registration": rawid["identity_digest"] == plan["identity_digest"],
    }
    frozen = {
        "original_registration": {"path": str(ORIGINAL), "sha256": n1.sha256_file(ORIGINAL)},
        "amendment_1": {"path": str(AMENDMENT_1), "sha256": n1.sha256_file(AMENDMENT_1)},
        "module": {"path": str(MODULE), "sha256": n1.sha256_file(MODULE)},
        "runner": {"path": str(RUNNER), "sha256": n1.sha256_file(RUNNER)},
    }
    payload = {
        "status": "amendment registered before any measurement; not run",
        "lane": "lane-n1",
        "amends": {
            "original": "n1_registration (f9a9130)",
            "amendment_1": "n1_registration_amendment_1 (8ebaa9a)",
        },
        "question": "Among the 56 frozen candidate factors, do the element-attribution predictions rank the "
        "genes "
        "the producers' Anderson-Darling test calls differentially expressed (BH-adjusted p < 0.05) after "
        "the "
        "factor's assigned CRISPRi perturbation above the non-responding genes of the same "
        "control-expression "
        "decile, more than promoter proximity does?",
        "question_changes": [
            "assigned perturbations are evaluated: no claim that every knockdown was effective, and no "
            "knockdown-quality filter is used, so the undocumented knockdown fields are not read",
            "ineffective perturbations may weaken the signal",
            "an absent response does not establish that a factor has no regulatory role",
            "no factor is selected because it produced strong responses",
        ],
        "resource": resource(),
        "coverage_from_identities_and_structure": coverage,
        "rules": {
            "label": n1.CONSTANTS_V3["label"],
            "missing": n1.CONSTANTS_V3["missing"],
            "row_rule": n1.CONSTANTS_V3["row_rule"],
            "universe": "the frozen 860 genes present once in the file's rows and with a finite control "
            "expression; per factor, minus the perturbed gene, genes within 10 kb of its TSS, and genes "
            "whose "
            "p-value is missing for that factor",
            "strata": n1.CONSTANTS_V3["strata"] + ", over the analysis genes; the values are read only in "
            "the "
            "authorised run, after this rule is frozen; the unweighted mean avoids the inferred "
            "num_cells_filtered",
            "comparison": "amendment 1's: stratified auroc_v2 primary and unstratified beside it, the "
            "equal-weight "
            "mean difference, the cluster bootstrap, the floor of 30 and the two separate criteria",
            "replaced": "the operational label |T| >= 3, the calibration gate and the knockdown eligibility "
            "rule "
            "no longer apply; the normalized pseudobulk is not read",
        },
        "run_v3": {
            "needs": "the owner's separate authorisation, given to scripts/n1_run_v3.py --authorisation",
            "reads_in_order": [
                "the freeze check: no data file is opened",
                f"download the published file ({n1.AD_SOURCE['bytes']:,} bytes) and the raw pseudobulk "
                f"({n1.SOURCE['files']['raw']['bytes']:,} bytes); check Figshare's md5 "
                f"({n1.AD_SOURCE['md5']}, "
                f"{n1.SOURCE['files']['raw']['md5']}) and record each sha256",
                "the raw pseudobulk's identities (obs/gene_transcript, var/gene_id), checked against the "
                "digest",
                "raw X at the 585 non-targeting rows, decoded at the 860 universe columns only (strata)",
                "the published file's labels pass: the header, checked against the registered digest, and "
                "the first "
                "field of every line; the whole file is decompressed and no value is parsed",
                "the values pass: at the covered universe rows each line is split into text fields and only "
                "the 56 "
                "principal columns are converted",
            ],
            "never_read": "the normalized pseudobulk, any knockdown field, any value outside the 56 "
            "principal "
            "columns at the universe rows (decompressed and split as text at those rows, never converted)",
            "result": "data/results/n1_result_amendment_2.json, written once",
        },
        "superseded_runs": "scripts/n1_run.py refuses; scripts/n1_run_v2.py now refuses by its own freeze "
        "check, "
        "because the module gained amendment 2's functions",
        "structure_reads": {
            "published_header": {"bytes_fetched": head["bytes_fetched"], "url": head["url"]},
            "raw_identities": {"bytes_fetched": rawid["bytes_fetched"], "datasets": rawid["datasets_read"]},
        },
        "ad_header_digest": n1.header_digest(cols),
        "frozen": frozen,
        "constants": json.loads(json.dumps(n1.CONSTANTS_V3)),
    }
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "Replogle et al. 2022 supplemental files, Figshare+ 21632564 (header only)",
                "version": "v1",
            },
            {
                "accession": "Replogle et al. 2022 processed Perturb-seq, Figshare+ 20029387 (identities "
                "only)",
                "version": "v1",
            },
            {
                "accession": "Replogle et al. 2022, Cell, STAR Methods (documentation)",
                "version": "PMC9380471",
            },
            {
                "accession": "GenomeOS n1_registration and n1_registration_amendment_1",
                "version": "f9a9130, 8ebaa9a",
            },
        ],
        "inputs": [
            mf.input_entry(ORIGINAL, partition=None),
            mf.input_entry(AMENDMENT_1, partition=None),
            mf.input_entry(AD_HEADER, partition=None, url=n1.AD_SOURCE["url"]),
            mf.input_entry(RAW_IDENTITY, partition=None, url=n1.SOURCE["files"]["raw"]["url"]),
        ],
        "assembly": "n/a: an amendment; no coordinate is read",
        "coordinates": "n/a: an amendment; no coordinate is read",
        "parameters": dict(n1.CONSTANTS_V3),
        "exclusions": [],
        "partitions": "n/a: an amendment registered before any measurement; nothing is evaluated",
    }
    p = save_result(n1.AMENDMENT_2, payload)
    print(f"wrote {p}: {json.dumps({k: coverage[k] for k in list(coverage)[:8]})}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
