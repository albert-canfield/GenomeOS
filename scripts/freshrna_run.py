# SPDX-License-Identifier: AGPL-3.0-or-later
"""M2i, the one run: the registered reading of M2h's fixed magnitude against RNA this project has
never read, then -- separately, afterwards and labelled as such -- the post-hoc quantities.

    uv run python scripts/freshrna_run.py [--post-hoc]

Without --post-hoc it computes and writes ONLY the registered reading (data/results/freshrna_result.json).
With --post-hoc, run after that result is committed, it adds `freshrna.POST_HOC`'s quantities to
data/results/freshrna_posthoc.json: the 1-Mb locus-cluster interval, the scale the fresh RNA would
imply, the gain as a fraction of the fresh variance, the nine polyA-plus cells on their own, the
per-cell values, M2h's training genes and all of 7h's genes, and the per-cell correlation between the
fresh and the already-read measurement.

NOTHING IS FITTED IN THIS LANE. The scale is the constant `freshrna.SCALE` = +0.2145911381525326 from
`ac5580b`; `freshrna.fit_scale` raises if anything asks for a fit, and this script never calls it. No
model request, no key, no money.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import freshrna as fr  # noqa: E402
from genomeos.attribution import readergate as rg  # noqa: E402
from genomeos.results import save_result  # noqa: E402

REGISTRATION = Path("data/results/freshrna_registration.json")
GENE_TABLE = Path("data/cache/gencode_genes.tsv")
OWN_CODE = (
    "genomeos/attribution/freshrna.py",
    "scripts/freshrna_fetch.py",
    "scripts/freshrna_register.py",
    "scripts/freshrna_run.py",
    "tests/test_freshrna.py",
)


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def fmt(x: float, places: int = 5) -> str:
    return f"{x:+.{places}f}"


def gene_starts() -> dict[str, int]:
    """gene_id without version -> start, for POST_HOC (1)'s 1-Mb locus clusters."""
    out: dict[str, int] = {}
    with GENE_TABLE.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        gi, si = header.index("gene_id"), header.index("start")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[gi].endswith("_PAR_Y"):
                continue
            out.setdefault(f[gi].split(".")[0], int(f[si]))
    return out


def fresh_tpm(sources: dict[str, Any]) -> tuple[dict[str, dict[str, float]], dict[str, Any]]:
    """cell -> (gene_id -> mean TPM over the chosen experiment's files); and the per-table counts.

    This is where a fresh value is first opened. It happens only after the registration is committed.
    """
    per_cell: dict[str, dict[str, float]] = {}
    counts: dict[str, Any] = {}
    for cell in fr.CELLS:
        entry = sources["biosamples"][cell]
        tables = []
        for f in entry["files"]:
            if sha256(f["path"]) != f["sha256"]:
                raise ValueError(f"{f['path']}: sha256 is not the registration's")
            table, c = rg.load_quant(Path(f["path"]))
            tables.append(table)
            counts[f["accession"]] = c
        per_cell[cell] = rg.mean_tpm(tables)
    return per_cell, counts


def assemble(seven_h: list[dict[str, Any]], per_cell: dict[str, dict[str, float]]) -> tuple[list[dict], dict]:
    """POPULATION: 7h's rows whose gene_id the fresh measurement resolves in all thirteen cells, with
    D' recomputed from the fresh quantifications and P carried unchanged."""
    starts = gene_starts()
    rows: list[dict[str, Any]] = []
    dropped: list[str] = []
    for r in seven_h:
        gid = r["gene_id"]
        if any(gid not in per_cell[c] for c in fr.CELLS):
            dropped.append(r["symbol"])
            continue
        y = [rg._y(per_cell[c][gid]) for c in fr.CELLS]  # noqa: SLF001 - item 7's endpoint, carried
        rows.append(
            {
                "symbol": r["symbol"],
                "gene_id": gid,
                "chrom": r["chrom"],
                "start": starts.get(gid, 0),
                "P": r["P"],
                "D": rg._centre(y),  # noqa: SLF001 - item 7's within-gene centring, carried
                "D_old": r["D"],
                "y": y,
            }
        )
    return rows, {"resolved": len(rows), "dropped": len(dropped), "dropped_symbols": sorted(dropped)[:50]}


def side(rows: list[dict], label: str) -> dict[str, Any]:
    """The registered endpoint and its interval on one set of rows, at the constant SCALE."""
    dp, pp, dd, n = fr.sums(rows)
    d = fr.delta_mse(rows)
    lo, hi = fr.bootstrap_delta(rows, fr.BOOTSTRAP["resamples"], fr.BOOTSTRAP["seed"])
    return {
        "label": label,
        "genes": len(rows),
        "pairs": n,
        "delta_mse": d,
        "ci": [lo, hi],
        "mean_D2": dd / n if n else float("nan"),
        "mean_P2": pp / n if n else float("nan"),
        "mean_D_times_P": dp / n if n else float("nan"),
    }


def _seven_h_rows() -> list[dict[str, Any]]:
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


def registered_run() -> dict[str, Any]:
    t0 = time.time()
    reg = json.loads(REGISTRATION.read_text())
    sources = json.loads(fr.SOURCES.read_text())
    seven_h = _seven_h_rows()
    defects: list[str] = []
    for p, want in reg["pinned_sha256"].items():
        if sha256(p) != want:
            defects.append(f"{p}: sha256 moved since the registration")
    if reg["scale"] != fr.SCALE:
        defects.append("the constant SCALE is not the registered one")
    held = sorted(r["symbol"] for r in seven_h if fr.is_held_out(r["symbol"]))
    digest = hashlib.sha256("\n".join(held).encode()).hexdigest()
    if digest != reg["primary_population"]["symbols_sha256"]:
        defects.append("the primary gene list is not the registered one")
    if defects:
        return {
            "outcome": "HARNESS_DEFECT",
            "reading": fr.reading("HARNESS_DEFECT", what="; ".join(defects)),
            "harness_check": {"passed": False, "defects": defects},
        }

    per_cell, counts = fresh_tpm(sources)  # the first fresh value opened in this project
    all_rows, attrition = assemble(seven_h, per_cell)
    primary = [r for r in all_rows if fr.is_held_out(r["symbol"])]
    harness = {
        "statement": fr.HARNESS,
        "pinned_sha256_unmoved": True,
        "scale_is_the_registered_constant": True,
        "primary_gene_list_digest": digest,
        "fresh_table_counts": counts,
        "fresh_tables": len(counts),
        "passed": True,
    }
    if len(primary) < fr.MIN_GENES:
        return {
            "outcome": "NOT_TESTABLE",
            "reading": fr.full_reading(
                fr.reading("NOT_TESTABLE", form=fr.FORM_WORDS, g=f"{len(primary):,}", min=f"{fr.MIN_GENES:,}")
            ),
            "harness_check": harness,
            "attrition": attrition,
            "_rows": all_rows,
        }

    s = side(primary, "primary: M2h's held-out genes on the fresh measurement")
    control = fr.permuted_control(primary, fr.PERMUTATIONS["count"], fr.PERMUTATIONS["seed"])
    worse_or_equal = sum(1 for c in control if c >= s["delta_mse"])
    p = (1 + worse_or_equal) / (1 + fr.PERMUTATIONS["count"])
    outcome = fr.decide(s["ci"][0], s["ci"][1])
    sentence = fr.reading(
        outcome,
        form=fr.FORM_WORDS,
        d=fmt(s["delta_mse"]),
        ci=fr.FORM_INTERVAL_WORDS,
        lo=fmt(s["ci"][0]),
        hi=fmt(s["ci"][1]),
        n=f"{s['pairs']:,}",
        g=f"{s['genes']:,}",
        p=f"{p:.6f}",
    )
    return {
        "outcome": outcome,
        "reading": fr.full_reading(sentence),
        "verdict_sentence": sentence,
        "harness_check": harness,
        "attrition": attrition,
        "primary": s,
        "control_p": p,
        "control_delta_mse_median": sorted(control)[len(control) // 2],
        "control_delta_mse_max": max(control),
        "variance_fraction": s["delta_mse"] / s["mean_D2"] if s["mean_D2"] else float("nan"),
        "seconds": round(time.time() - t0, 1),
        "_rows": all_rows,
    }


def post_hoc(rows: list[dict], primary: list[dict]) -> dict[str, Any]:
    """POST_HOC: computed after the registered reading is written and committed. Decides nothing."""
    import numpy as np

    sources = json.loads(fr.SOURCES.read_text())
    polya = [i for i, c in enumerate(fr.CELLS) if sources["assays"][c] == "polyA plus RNA-seq"]
    nine = [
        {
            "symbol": r["symbol"],
            "chrom": r["chrom"],
            "start": r["start"],
            "P": [r["P"][i] for i in polya],
            "D": rg._centre([r["y"][i] for i in polya]),  # noqa: SLF001
        }
        for r in primary
    ]
    clo, chi, nclusters = fr.bootstrap_delta_clustered(
        primary, fr.BOOTSTRAP["resamples"], fr.BOOTSTRAP["seed"]
    )
    dp, pp, dd, n = fr.sums(primary)
    implied = fr.implied_scale(primary)
    per_cell = {}
    for i, cell in enumerate(fr.CELLS):
        num = sum(2 * fr.SCALE * r["P"][i] * r["D"][i] - fr.SCALE**2 * r["P"][i] ** 2 for r in primary)
        per_cell[cell] = num / len(primary)
    corr = {}
    for i, cell in enumerate(fr.CELLS):
        a = np.array([r["D"][i] for r in primary])
        b = np.array([r["D_old"][i] for r in primary])
        corr[cell] = float(np.corrcoef(a, b)[0, 1])
    train = [r for r in rows if not fr.is_held_out(r["symbol"])]
    return {
        "label": "POST HOC. Computed after the registered reading was written and committed. Decides nothing",
        "statement": fr.POST_HOC,
        "dependence_robust_interval": {
            "label": "POST HOC (1): the 1-Mb locus-cluster bootstrap, deciding nothing",
            "span_bp": 1_000_000,
            "clusters": nclusters,
            "delta_mse": fr.delta_mse(primary),
            "ci": [clo, chi],
            "registered_gene_level_ci": fr.bootstrap_delta(
                primary, fr.BOOTSTRAP["resamples"], fr.BOOTSTRAP["seed"]
            ),
            "note": (
                "genes that share elements or a locus are not independent, so the registered gene-level"
                " interval may be too narrow; this one resamples 1-Mb locus clusters instead. The"
                " registered interval remains the gene-level one and this figure changes no verdict"
            ),
        },
        "implied_scale": {
            "label": "POST HOC (2): the scale the FRESH RNA would imply, beside the constant used",
            "scale_used_fixed_constant": fr.SCALE,
            "scale_the_fresh_rna_would_imply": implied,
            "ratio_implied_over_used": (implied / fr.SCALE) if implied is not None else None,
            "m2h_training_scale_interval": json.loads(Path(fr.SCALE_SOURCE).read_text())["per_form"][
                "multiplicative"
            ]["scale_ci"],
            "note": (
                "DESCRIPTIVE AND POST HOC. It is sum(D'*P)/sum(P^2) on the primary population, computed"
                " only after the registered reading was committed; it was NOT fitted into the registered"
                " reading, which used the constant +0.2146 alone, and it is not a verdict of any kind"
            ),
        },
        "size_on_this_measurement": {
            "label": "POST HOC (3): SIZE's 2% figure recomputed on the fresh measurement",
            "mean_D2": dd / n,
            "mean_P2": pp / n,
            "mean_D_times_P": dp / n,
            "delta_mse": fr.delta_mse(primary),
            "fraction_of_variance": fr.delta_mse(primary) / (dd / n),
            "best_possible_gain": (dp / n) ** 2 / (pp / n),
            "implied_correlation": (dp / n) / ((dd / n) ** 0.5 * (pp / n) ** 0.5),
        },
        "nine_polya_cells": {
            "label": "POST HOC (4): the nine cells whose fresh experiment is polyA plus, D' recentred",
            "cells": [fr.CELLS[i] for i in polya],
            **{k: v for k, v in side(nine, "nine polyA plus cells").items() if k != "label"},
        },
        "per_cell_delta_mse": {"label": "POST HOC (5)", "values": per_cell},
        "other_populations": {
            "label": "POST HOC (6): M2h's training genes and all of 7h's genes",
            "m2h_training_genes": side(train, "M2h's training genes on the fresh measurement"),
            "all_seven_h_genes": side(rows, "all of 7h's genes on the fresh measurement"),
        },
        "fresh_versus_read_correlation": {
            "label": "POST HOC (7): per cell, the correlation of the fresh D' with item 7's D",
            "values": corr,
            "mean": sum(corr.values()) / len(corr),
            "note": "a descriptive check that the two measurements measure the same thing",
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--post-hoc", action="store_true", help="the POST_HOC stage, after the result is committed"
    )
    args = ap.parse_args()
    sources = json.loads(fr.SOURCES.read_text())
    pinned = (
        fr.P_SOURCE,
        fr.SCALE_SOURCE,
        REGISTRATION.as_posix(),
        fr.SOURCES.as_posix(),
        GENE_TABLE.as_posix(),
    )
    inputs = [mf.input_entry(p, partition=None) for p in pinned]
    inputs += [
        mf.input_entry(f["path"], partition=None) for c in sources["biosamples"].values() for f in c["files"]
    ]
    manifest = {
        "sources": [
            {"accession": "ac5580b", "version": "M2h, the fixed scale"},
            {"accession": "90e5c86", "version": "7h, the committed P rows"},
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
        "exclusions": [
            "genes 7h did not admit",
            "genes the fresh measurement does not resolve in all thirteen cells",
        ],
        "partitions": "the primary population is M2h's held-out genes, as the registration pins them",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }
    out = registered_run()
    rows = out.pop("_rows", [])
    if args.post_hoc:
        result = json.loads(Path("data/results/freshrna_result.json").read_text())
        if result["outcome"] != out["outcome"]:
            raise SystemExit(
                "refused: the registered reading does not reproduce; nothing post hoc is written"
            )
        primary = [r for r in rows if fr.is_held_out(r["symbol"])]
        payload = {
            "lane": fr.LANE,
            "date": time.strftime("%Y-%m-%d"),
            "what": "M2i POST HOC: computed after the registered reading was committed; decides nothing",
            "registered_reading": result["reading"],
            "registered_outcome": result["outcome"],
            "post_hoc": post_hoc(rows, primary),
            "alphagenome_requests": 0,
            "money": "none: public files only, no key read",
        }
        print("saved", save_result("freshrna_posthoc", payload, manifest=manifest))
        return 0
    payload = {
        "lane": fr.LANE,
        "date": time.strftime("%Y-%m-%d"),
        "registration": REGISTRATION.as_posix(),
        "registration_sha256": sha256(REGISTRATION),
        "what": (
            "M2i: the registered reading of M2h's fixed magnitude (+0.2146, a constant from ac5580b)"
            " against ENCODE RNA-seq of the same thirteen reader cell types that this project had never"
            " read. NOTHING IS FITTED IN THIS LANE"
        ),
        "question": fr.QUESTION,
        "nothing_is_fitted": fr.NOTHING_IS_FITTED,
        "scale": fr.SCALE,
        "scale_origin": fr.SCALE_ORIGIN,
        "scope": fr.SCOPE,
        "size": fr.SIZE,
        "what_is_fresh": fr.WHAT_IS_FRESH,
        "non_overlap": sources["non_overlap"],
        "assays": sources["assays"],
        "fresh_experiments": {k: v["experiment"] for k, v in sources["biosamples"].items()},
        "population": fr.POPULATION,
        "metric": fr.METRIC,
        "control": fr.CONTROL,
        "decision": fr.DECISION,
        "falsifier": fr.FALSIFIER,
        "stop_rule": fr.STOP_RULE,
        "limitations": list(fr.LIMITATIONS),
        "cells": list(fr.CELLS),
        "m2h_reading": json.loads(Path(fr.SCALE_SOURCE).read_text())["per_form"]["multiplicative"]["reading"],
        "alphagenome_requests": 0,
        "money": "none: public files only, no key read",
        **out,
    }
    print(out["outcome"])
    print(out["reading"][:1200])
    print("saved", save_result("freshrna_result", payload, manifest=manifest))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
