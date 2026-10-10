# SPDX-License-Identifier: AGPL-3.0-or-later
"""M2i registration: the fixed form, the FIXED SCALE, the population, the endpoint, the control, the
decision and the falsifier, written and committed BEFORE any value of the fresh RNA is read.

    uv run python scripts/freshrna_register.py

Writes data/results/freshrna_registration.json: `genomeos.attribution.freshrna.registered_fields()`,
the sha256 of every file the run may read, and the parts of the HARNESS that need no fresh value --
7h's committed genes table reproduces 7h's committed Delta MSE at scale 1, M2h's committed held-out
counts come back from the table and the split, the constant SCALE equals `ac5580b`'s committed scale,
and the fresh sources file's non-overlap holds. It also pins the primary gene list by the sha256 of
its sorted symbols, so the run can prove it read the registered population and no other.

NOTHING IS FITTED IN THIS LANE (`freshrna.NOTHING_IS_FITTED`): this script computes no parameter, and
`freshrna.fit_scale` raises if anything asks it to. No fresh TPM value is opened here. No model
request, no key, no money.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import freshrna as fr  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "freshrna_registration"
OWN_CODE = (
    "genomeos/attribution/freshrna.py",
    "scripts/freshrna_fetch.py",
    "scripts/freshrna_register.py",
    "scripts/freshrna_run.py",
    "tests/test_freshrna.py",
)
PINNED = (
    fr.P_SOURCE,
    fr.SCALE_SOURCE,
    "data/results/multrna_registration.json",
    "data/results/magcal_registration.json",
    fr.SOURCES.as_posix(),
)


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def seven_h_rows() -> list[dict[str, Any]]:
    """7h's committed genes table, as rows: symbol, chrom, the thirteen P and the thirteen old D."""
    d = json.loads(Path(fr.P_SOURCE).read_text())
    cols = d["genes_table"]["columns"]
    i = {c: cols.index(c) for c in ("symbol", "gene_id", "chrom", "P", "D")}
    return [
        {
            "symbol": r[i["symbol"]],
            "gene_id": r[i["gene_id"]],
            "chrom": r[i["chrom"]],
            "P": r[i["P"]],
            "D": r[i["D"]],
        }
        for r in d["genes_table"]["rows"]
    ]


def harness(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """HARNESS (2) to (5): everything that needs no fresh value, before the registration is written."""
    seven_h = json.loads(Path(fr.P_SOURCE).read_text())["multiplicative"]
    magcal = json.loads(Path(fr.SCALE_SOURCE).read_text())["per_form"]["multiplicative"]
    sources = json.loads(fr.SOURCES.read_text())
    pairs = sum(len(r["D"]) for r in rows)
    recomputed = fr.delta_mse(rows, scale=1.0)
    held = [r for r in rows if fr.is_held_out(r["symbol"])]
    train = [r for r in rows if not fr.is_held_out(r["symbol"])]
    out = {
        "statement": fr.HARNESS,
        "seven_h_genes": len(rows),
        "seven_h_genes_committed": seven_h["genes"],
        "seven_h_pairs": pairs,
        "seven_h_pairs_committed": seven_h["pairs"],
        "thirteen_everywhere": all(len(r["P"]) == 13 and len(r["D"]) == 13 for r in rows),
        "symbols_unique": len({r["symbol"] for r in rows}) == len(rows),
        "symbols_non_empty": all(r["symbol"] for r in rows),
        "delta_mse_at_scale_1_recomputed": recomputed,
        "delta_mse_at_scale_1_committed": seven_h["delta_mse"],
        "delta_mse_reproduced": abs(recomputed - seven_h["delta_mse"]) <= fr.TOLERANCE,
        "tolerance": fr.TOLERANCE,
        "scale_constant": fr.SCALE,
        "scale_committed_in_ac5580b": magcal["scale"],
        "scale_matches": magcal["scale"] == fr.SCALE,
        "m2h_held_out_genes": len(held),
        "m2h_held_out_genes_committed": magcal["split"]["held_out_genes"],
        "m2h_training_genes": len(train),
        "m2h_training_genes_committed": magcal["split"]["training_genes"],
        "m2h_held_out_pairs": sum(len(r["D"]) for r in held),
        "fresh_sources_disjoint": (
            sources["non_overlap"]["disjoint_experiments"]
            and sources["non_overlap"]["disjoint_files"]
            and sources["non_overlap"]["disjoint_libraries"]
        ),
        "fresh_files_on_disk": all(
            Path(f["path"]).exists() and sha256(f["path"]) == f["sha256"]
            for c in sources["biosamples"].values()
            for f in c["files"]
        ),
    }
    out["passed"] = all(
        (
            out["seven_h_genes"] == out["seven_h_genes_committed"],
            out["seven_h_pairs"] == out["seven_h_pairs_committed"],
            out["thirteen_everywhere"],
            out["symbols_unique"],
            out["symbols_non_empty"],
            out["delta_mse_reproduced"],
            out["scale_matches"],
            out["m2h_held_out_genes"] == out["m2h_held_out_genes_committed"],
            out["m2h_training_genes"] == out["m2h_training_genes_committed"],
            out["fresh_sources_disjoint"],
            out["fresh_files_on_disk"],
        )
    )
    return out


def main() -> int:
    t0 = time.time()
    rows = seven_h_rows()
    h = harness(rows)
    for k, v in h.items():
        if isinstance(v, bool) or k.endswith(("_genes", "_pairs", "_committed", "_recomputed")):
            print(f"  {k}: {v}")
    if not h["passed"]:
        raise SystemExit("refused: the harness does not pass; nothing is registered")
    held = sorted(r["symbol"] for r in rows if fr.is_held_out(r["symbol"]))
    train = sorted(r["symbol"] for r in rows if not fr.is_held_out(r["symbol"]))
    sources = json.loads(fr.SOURCES.read_text())
    payload = {
        "lane": fr.LANE,
        "date": time.strftime("%Y-%m-%d"),
        "registered": fr.REGISTERED,
        "what": (
            "M2i: M2h's fitted magnitude for the multiplicative composition, taken as a FIXED CONSTANT,"
            " tested against ENCODE RNA-seq of the same thirteen reader cell types whose values this"
            " project has never read. Registered before any fresh value was opened"
        ),
        "m2h": {
            "commit": "ac5580b",
            "registration": "fbce775",
            "code": "2b0b604",
            "reading": json.loads(Path(fr.SCALE_SOURCE).read_text())["per_form"]["multiplicative"]["reading"],
            "prior_exposure_clause": json.loads(Path(fr.SCALE_SOURCE).read_text())["prior_exposure_clause"],
        },
        **fr.registered_fields(),
        "primary_population": {
            "genes": len(held),
            "pairs": len(held) * len(fr.CELLS),
            "symbols_sha256": hashlib.sha256("\n".join(held).encode()).hexdigest(),
            "note": (
                "M2h's held-out genes of the multiplicative composition, before the fresh measurement's"
                " own resolution filter; the run reports the attrition and re-derives this digest"
            ),
        },
        "reported_beside": {
            "m2h_training_genes": len(train),
            "training_symbols_sha256": hashlib.sha256("\n".join(train).encode()).hexdigest(),
            "all_seven_h_genes": len(rows),
        },
        "harness_check": h,
        "fresh_rna": {
            "assays": sources["assays"],
            "experiments": {k: v["experiment"] for k, v in sources["biosamples"].items()},
            "files": {k: [f["accession"] for f in v["files"]] for k, v in sources["biosamples"].items()},
            "non_overlap": sources["non_overlap"],
            "values_read_so_far": sources["values_read"],
            "licence": sources["licence"],
            "bytes_total": sources["bytes_total"],
        },
        "pinned_sha256": {p: sha256(p) for p in PINNED},
        "alphagenome_requests": 0,
        "money": "none: public files only, no key read",
        "seconds": round(time.time() - t0, 1),
    }
    inputs = [mf.input_entry(p, partition=None) for p in PINNED]
    inputs += [
        mf.input_entry(f["path"], partition=None) for c in sources["biosamples"].values() for f in c["files"]
    ]
    manifest = {
        "sources": [
            {"accession": "ac5580b", "version": "M2h, the fitted scale of the multiplicative composition"},
            {"accession": "90e5c86", "version": "7h, the committed P and old D of the multiplicative form"},
        ]
        + [
            {"accession": c["experiment"], "version": f"ENCODE portal, released {c['date_released']}"}
            for c in sources["biosamples"].values()
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38, GENCODE V29 gene models",
        "coordinates": "n/a: gene-level quantifications keyed by gene_id",
        "parameters": {
            "scale": fr.SCALE,
            "nothing_is_fitted": fr.NOTHING_IS_FITTED,
            "interval": list(fr.FORM_INTERVAL),
            "resamples": fr.BOOTSTRAP["resamples"],
            "seed": fr.BOOTSTRAP["seed"],
            "permutations": fr.PERMUTATIONS,
        },
        "exclusions": ["no fresh TPM value is opened by this script"],
        "partitions": (
            "the primary population is M2h's held-out genes, pinned here by the sha256 of its sorted"
            " symbols; M2h's training genes and all of 7h's genes are reported beside, deciding nothing"
        ),
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }
    print("saved", save_result(RESULT, payload, manifest=manifest))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
