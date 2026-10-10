# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-maxform, the registered run: `or` on 7d's population, one scale fitted on training genes.

    uv run python scripts/maxform_run.py [--workers 4]

Runs only from the committed registration (data/results/maxform_registration.json); refuses if any
registered field differs from genomeos/attribution/maxform.py or any pinned file's sha256 differs from the
file on disk.

Order, as registered (maxform.DECISION):
  1. the 24 compiled programs are gated as 7d gated them, every gene runnable and informative under `or`
     is run by the runtime (activator_combination='or'), and the HARNESS is checked on all of them --
     each run against its closed form and every single-removal observation reproduced. Any miss stops the
     run here, before the RNA is opened (HARNESS_DEFECT), and the census counts must be reproduced too;
  2. only then are the RNA quantifications opened and the population assembled (POPULATION (f));
  3. the split is applied (maxform.SPLIT, lane-magcal's sha256 of the symbol);
  4. one scale is fitted on the TRAINING pairs alone and bootstrapped on training genes;
  5. only then is the held-out Delta MSE computed at that scale, with its 98.75% interval (Bonferroni over
     the family of four) and lane-magcal's within-gene cell permutation control, and magcal's four-way
     table applied as written;
  6. the head-to-heads against 7d's mean and 7h's multiplicative composition, and the post-hoc 1-Mb
     locus-cluster bootstrap, are reported beside it as descriptives that decide nothing.

No model request, no network, no money, 0 AlphaGenome requests.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from multiprocessing import get_context
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import combine2 as c2  # noqa: E402
from genomeos.attribution import compile as cp  # noqa: E402
from genomeos.attribution import magcal as mc  # noqa: E402
from genomeos.attribution import maxform as xf  # noqa: E402
from genomeos.attribution import readergate as rg  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import magcal_run as mcrun  # noqa: E402
from maxform_census import CHROMS, OWN_CODE, rate_of  # noqa: E402
from maxform_register import CENSUS, MC_RESULT  # noqa: E402
from maxform_register import inputs as registration_inputs  # noqa: E402
from readergate_census import peak_files  # noqa: E402
from readergate_run import load_rna  # noqa: E402
from runrules_census import load_rates  # noqa: E402

RESULT = "maxform_result"
REGISTRATION = RESULTS_DIR / "maxform_registration.json"
FORM = xf.FORM


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def check_registration() -> dict[str, Any]:
    """Refuse to run unless every BOUND field and every pinned sha256 is what was registered."""
    if not REGISTRATION.exists():
        raise SystemExit(f"refused: {REGISTRATION} does not exist; register before running")
    reg = json.loads(REGISTRATION.read_text())
    fields = xf.registered_fields()
    bad = [k for k in xf.BOUND if reg.get(k) != fields[k]]
    if bad:
        raise SystemExit(f"refused: the registration and maxform.py disagree on {', '.join(bad)}")
    group = mf.files_entry("reader v1 DNase narrowPeak files", [p for p in peak_files(CHROMS) if p.exists()])
    moved = [
        p
        for p, digest in reg["pinned_sha256"].items()
        if p != group["path"] and (not Path(p).is_file() or _sha(Path(p)) != digest)
    ]
    if reg["pinned_sha256"].get(group["path"]) != group["sha256"]:
        moved.append(group["path"])
    if moved:
        raise SystemExit(f"refused: pinned inputs changed or missing since the registration: {moved}")
    return reg


def _chromosome(chrom_and_genes: tuple[str, set[str]]) -> dict[str, Any]:
    """One program's genes runnable and informative under `or`, run in the thirteen reader cells."""
    from genomeos.lang.parser import parse

    chrom, genes = chrom_and_genes
    t0 = time.time()
    rates, borrowed, prov, lent = load_rates()
    text = (Path(COMPILED_DIR) / f"noncoding_{chrom}.bio").read_text()
    gated, info = cp.gate_program(text, chrom, selection=xf.SELECTION)
    module = parse(gated, f"noncoding_{chrom}_reader")
    out: dict[str, Any] = {
        "chrom": chrom,
        "program": {k: info[k] for k in ("base_sha256", "gated_sha256", "elements", "rules_gated")},
        "runs": [],
        "fitted": 0,
        "informative": 0,
    }
    for gr in c2.gene_rules(module, genes, chrom, rg.compiled_cells(text)):
        _fits, why = c2.classify(gr)
        if why:
            continue
        f = xf.fit(FORM, gr.rho)
        if not f.ok:
            continue
        out["fitted"] += 1
        s_closed = xf.predicted_log2(f, gr)
        if max(s_closed) - min(s_closed) <= xf.INFORMATIVE_LOG2:
            continue
        out["informative"] += 1
        rate, label = rate_of(gr.gene, rates, borrowed, prov, lent)
        got = xf.run_gene(rg.gated_gene_module(module, gr.gene), gr, f, rate, label)
        n_open = [sum(1 for el in gr.elements if gr.open_in(j)[el]) for j in range(len(xf.CELLS))]
        out["runs"].append(
            {
                "gene": gr.gene,
                "symbol": gr.symbol,
                "chrom": chrom,
                "n_rules": len(gr.elements),
                "n_inhibiting": sum(1 for r in gr.rho.values() if r > 1.0),
                "elements": list(gr.elements),
                "states": {el: list(gr.states[el]) for el in gr.elements},
                "open_share": [x / len(gr.elements) for x in n_open],
                "weights": dict(f.weights),
                "e": dict(f.e),
                "q": f.q,
                "z": f.z,
                "s": got["s"],
                "closed_s": got["closed_s"],
                "harness_miss": got["harness_miss"],
                "tier": label["tier"],
            }
        )
    out["seconds"] = round(time.time() - t0, 1)
    return out


def collect(workers: int) -> tuple[list[dict[str, Any]], dict[str, Any], dict[str, int]]:
    """Every gene of the 24 programs runnable and informative under `or`, run; the programs' digests."""
    texts = {c: (Path(COMPILED_DIR) / f"noncoding_{c}.bio").read_text() for c in CHROMS}
    multi, _why = c2.multi_rule_genes(texts)
    by_chrom: dict[str, set[str]] = {c: set() for c in CHROMS}
    for g, c in multi.items():
        by_chrom[c].add(g)
    order = sorted(CHROMS, key=lambda c: -(Path(COMPILED_DIR) / f"noncoding_{c}.bio").stat().st_size)
    res: dict[str, dict] = {}
    if workers > 1:
        with get_context("spawn").Pool(workers) as pool:
            for r in pool.imap_unordered(_chromosome, [(c, by_chrom[c]) for c in order]):
                res[r["chrom"]] = r
                print(r["chrom"], len(r["runs"]), r["seconds"], flush=True)
    else:
        for c in order:
            r = _chromosome((c, by_chrom[c]))
            res[r["chrom"]] = r
            print(r["chrom"], len(r["runs"]), r["seconds"], flush=True)
    runs = [x for c in CHROMS for x in res[c]["runs"]]
    programs = {c: res[c]["program"] for c in CHROMS}
    counts = {
        "fitted": sum(res[c]["fitted"] for c in CHROMS),
        "informative": sum(res[c]["informative"] for c in CHROMS),
        "run": len(runs),
    }
    return runs, programs, counts


def harness_check(runs: list[dict[str, Any]], counts: dict[str, int]) -> dict[str, Any]:
    """HARNESS, before the RNA: every run equals its closed form and the census counts are reproduced."""
    bad = [r["gene"] for r in runs if r["harness_miss"]]
    worst = max(
        (abs(a - b) for r in runs for a, b in zip(r["s"], r["closed_s"], strict=True)),
        default=None,
    )
    committed = json.loads(CENSUS.read_text())["reader_selection"]["this_lanes_forms"][FORM]
    problems = [f"{len(bad)} genes miss the harness: {bad[:20]}"] if bad else []
    for k in ("fitted", "informative", "run"):
        if counts[k] != committed[k]:
            problems.append(f"census {k}: run {counts[k]}, committed {committed[k]}")
    return {
        "gene_runs": len(runs),
        "failing": len(bad),
        "failing_examples": bad[:20],
        "max_abs_run_minus_closed_form_log2": worst,
        "census_counts_recomputed": counts,
        "problems": problems,
    }


def gene_starts() -> dict[str, int]:
    """symbol -> start, from the pinned gene table, for the POST HOC locus-cluster bootstrap only."""
    out: dict[str, int] = {}
    with rg.GENE_TABLE.open() as fh:
        head = fh.readline().rstrip("\n").split("\t")
        i_sym, i_start = head.index("symbol"), head.index("start")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            out.setdefault(parts[i_sym], int(parts[i_start]))
    return out


def _sides(rs: list[dict[str, Any]], scale: float) -> dict[str, Any]:
    dp, pp, dd, n = mc._sums(rs)
    return {
        "genes": len(rs),
        "pairs": n,
        "mean_D2": dd / n if n else None,
        "mean_P2": pp / n if n else None,
        "mean_scaled_P2": scale * scale * pp / n if n else None,
        "mean_D_times_P": dp / n if n else None,
    }


def score(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """The registered scale, the held-out comparison, the control and the verdict (maxform.DECISION 3-5)."""
    words = xf.FORM_WORDS[FORM]
    train, held = mc.split_rows(rows)
    scale = mc.fit_scale(train)
    base: dict[str, Any] = {
        "form": FORM,
        "form_words": words,
        "split": {"training_genes": len(train), "held_out_genes": len(held), "rule": xf.SPLIT},
    }
    if scale is None:
        clause = mc.scale_clause("not_fittable")
        return {
            **base,
            "outcome": "NOT_FITTABLE",
            "scale": None,
            "reading": mc.reading("NOT_FITTABLE", form=words, scale=clause),
        }
    alo, ahi = xf.bootstrap_scale(train, xf.BOOTSTRAP["resamples"], xf.BOOTSTRAP["seed"])
    tn = sum(len(r["D"]) for r in train)
    clause = mc.scale_clause(
        "fitted",
        a=f"{scale:+.4g}",
        ci=xf.FORM_INTERVAL_WORDS,
        alo=f"{alo:+.4g}",
        ahi=f"{ahi:+.4g}",
        tg=f"{len(train):,}",
        tn=f"{tn:,}",
    )
    negative = mc.scale_is_negative(alo, ahi)
    if negative:
        clause += "; " + mc.scale_clause(
            "negative",
            a=f"{scale:+.4g}",
            ci=xf.FORM_INTERVAL_WORDS,
            alo=f"{alo:+.4g}",
            ahi=f"{ahi:+.4g}",
            tg=f"{len(train):,}",
        )
    out: dict[str, Any] = {
        **base,
        "scale": scale,
        "scale_ci": [alo, ahi],
        "scale_interval_includes_zero": bool(alo <= 0.0 <= ahi),
        "scale_is_negative": bool(negative),
        "scale_clause": clause,
        "training_delta_mse_at_fitted_scale": mc.delta_mse(train, scale),
        "training_delta_mse_at_scale_1": mc.delta_mse(train, 1.0),
    }
    if len(held) < xf.MIN_GENES:
        out["outcome"] = "NOT_TESTABLE"
        out["reading"] = mc.reading(
            "NOT_TESTABLE", form=words, g=f"{len(held):,}", min=xf.MIN_GENES, scale=clause
        )
        return out
    dm = mc.delta_mse(held, scale)
    lo, hi = xf.bootstrap_delta(held, scale, xf.BOOTSTRAP["resamples"], xf.BOOTSTRAP["seed"])
    perm = mc.permuted_control(train, held, xf.PERMUTATIONS["count"], xf.PERMUTATIONS["seed"])
    p_ctrl = (1 + sum(1 for x in perm if x >= dm)) / (1 + len(perm))
    kind = mc.decide(lo, hi, alo, ahi)
    n = sum(len(r["D"]) for r in held)
    held_side = _sides(held, scale)
    out.update(
        {
            "outcome": kind,
            "held_out_delta_mse": dm,
            "held_out_ci": [lo, hi],
            "control_p": p_ctrl,
            "reading": mc.reading(
                kind,
                form=words,
                d=f"{dm:+.4g}",
                ci=xf.FORM_INTERVAL_WORDS,
                lo=f"{lo:+.4g}",
                hi=f"{hi:+.4g}",
                n=f"{n:,}",
                g=f"{len(held):,}",
                p=f"{p_ctrl:.3g}",
                scale=clause,
            ),
            # HONEST_SIZE: the gain against the variance it is read against, carried with any 'better'
            "honest_size": {
                "rule": xf.HONEST_SIZE,
                "mean_D2_held_out": held_side["mean_D2"],
                "delta_mse_as_a_fraction_of_mean_D2": (
                    dm / held_side["mean_D2"] if held_side["mean_D2"] else None
                ),
            },
            "post_hoc_cluster_bootstrap": xf.cluster_bootstrap_delta(
                held, scale, xf.BOOTSTRAP["resamples"], xf.BOOTSTRAP["seed"]
            ),
            "descriptives_not_verdicts": {
                "held_out_delta_mse_at_scale_1": mc.delta_mse(held, 1.0),
                "held_out_ci_at_scale_1": list(
                    xf.bootstrap_delta(held, 1.0, xf.BOOTSTRAP["resamples"], xf.BOOTSTRAP["seed"])
                ),
                "whole_population_delta_mse_at_scale_1": mc.delta_mse(rows, 1.0),
                "training_side": _sides(train, scale),
                "held_out_side": held_side,
                "held_out_delta_mse_by_cell": {
                    c: sum(2 * scale * r["P"][j] * r["D"][j] - scale * scale * r["P"][j] ** 2 for r in held)
                    / len(held)
                    for j, c in enumerate(xf.CELLS)
                },
                "held_out_delta_mse_by_rule_count": {
                    k: mc.delta_mse(rs, scale) if rs else None
                    for k, rs in {
                        "2-3": [r for r in held if r["n_rules"] <= 3],
                        "4-9": [r for r in held if 4 <= r["n_rules"] <= 9],
                        "10+": [r for r in held if r["n_rules"] >= 10],
                    }.items()
                },
                "held_out_delta_mse_by_inhibitor": {
                    k: mc.delta_mse(rs, scale) if rs else None
                    for k, rs in {
                        "has_inhibiting_rule": [r for r in held if r["n_inhibiting"] > 0],
                        "activating_only": [r for r in held if r["n_inhibiting"] == 0],
                    }.items()
                },
                "predicted_range_log2_quantiles_held_out": mc.quantiles(
                    [max(r["P"]) - min(r["P"]) for r in held]
                ),
                "control_delta_mse_median": sorted(perm)[len(perm) // 2],
                "control_delta_mse_95th_percentile": sorted(perm)[int(0.95 * len(perm)) - 1],
                "control_p_floor": 1 / (1 + xf.PERMUTATIONS["count"]),
            },
        }
    )
    return out


def head_to_heads(rows: list[dict[str, Any]], scale: float | None) -> dict[str, Any]:
    """DESCRIPTIVE (maxform.DECISION 6), deciding nothing: against 7d's mean and 7h's multiplicative form
    on the genes both run, from their committed rows, uncalibrated and at each form's own fitted scale."""
    own = {r["gene"]: r for r in rows}
    other = mcrun.load_rows()
    committed = json.loads(MC_RESULT.read_text())
    out: dict[str, Any] = {"decides": "nothing", "note": xf.DESCRIPTIVES[1]}
    for form in ("mean", "multiplicative"):
        rs = {r["gene"]: r for r in other[form]}
        a_other = mc.fit_scale(mc.split_rows(other[form])[0])
        out[form] = {
            "their_scale_recomputed": a_other,
            "their_scale_committed": (committed.get("per_form") or {}).get(form, {}).get("scale"),
            "uncalibrated": xf.head_to_head(
                own, rs, xf.BOOTSTRAP["resamples"], xf.BOOTSTRAP["seed"], (1.0, 1.0)
            ),
        }
        if scale is not None and a_other is not None:
            out[form]["at_each_forms_own_scale"] = xf.head_to_head(
                own, rs, xf.BOOTSTRAP["resamples"], xf.BOOTSTRAP["seed"], (scale, a_other)
            )
    return out


def _table(rows: list[dict[str, Any]]) -> dict[str, Any]:
    cols = ("gene", "symbol", "chrom", "start", "n_rules", "n_inhibiting", "P", "D")
    return {
        "columns": list(cols),
        "rows": [[r.get(c) for c in cols] for r in rows],
        "note": (
            "the 7d-shaped genes table this lane's reading is computed from: P is the centred log2 level"
            " under `or`, D the centred measured log2(TPM+1); `start` is carried for the post-hoc"
            " locus-cluster bootstrap only"
        ),
    }


def manifest(entries: list[dict[str, Any]], rna_read: bool) -> dict[str, Any]:
    return {
        "sources": [
            {"accession": "the 24 compiled non-coding programs 7d ran", "version": "pinned by sha256"},
            {"accession": "reader v1 ENCODE DNase-seq narrowPeak files, 13 biosamples", "version": "pinned"},
            {"accession": "ENCODE polyA plus RNA-seq gene quantifications, V29", "version": "pinned"},
            {"accession": "data/results/gene_rates.json (mouse NIH 3T3, e921fae)", "version": "pinned"},
        ],
        "inputs": entries,
        "input_count": len(entries),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "form": FORM,
            "split": xf.SPLIT,
            "min_genes": xf.MIN_GENES,
            "form_interval": list(xf.FORM_INTERVAL),
            "bootstrap": xf.BOOTSTRAP,
            "permutations": xf.PERMUTATIONS,
            "control_alpha": xf.CONTROL_ALPHA,
            "free_parameters": "one real scale for the form; no per-gene and no per-cell parameter",
        },
        "exclusions": ["as POPULATION and POPULATION_NOTE state"],
        "partitions": {
            "rule": xf.SPLIT,
            "training": "even first 32 bits of sha256(symbol): the scale is fitted here only",
            "held_out": "odd first 32 bits of sha256(symbol): the registered endpoint is read here only",
        }
        if rna_read
        else "n/a: stopped before the RNA was opened",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    t0 = time.time()
    reg = check_registration()
    runs, programs, counts = collect(args.workers)
    harness = harness_check(runs, counts)
    entries = [*registration_inputs(), mf.input_entry(REGISTRATION, partition=None)]
    head: dict[str, Any] = {
        "lane": xf.LANE,
        "date": time.strftime("%Y-%m-%d"),
        "roadmap_item": xf.ROADMAP_ITEM,
        "registration": "data/results/maxform_registration.json",
        "registered_date": reg["date"],
        "form": FORM,
        "harness": harness,
        "programs": programs,
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }
    if harness["problems"]:
        head.update(
            {
                "outcome": "HARNESS_DEFECT",
                "reading": mc.reading("HARNESS_DEFECT", what="; ".join(harness["problems"])),
                "rna_read": "none: the run stopped before the RNA was opened",
                "seconds": round(time.time() - t0, 1),
            }
        )
        print(head["reading"])
        print("saved", save_result(RESULT, head, manifest=manifest(entries, rna_read=False)))
        return 1
    tpm, rna_counts = load_rna()
    ids = rg.symbol_gene_ids()
    rows, why = c2.assemble(runs, tpm, ids)
    starts = gene_starts()
    for r in rows:
        r["start"] = starts.get(r["symbol"], 0)
    per = score(rows)
    h2h = head_to_heads(rows, per.get("scale"))
    reading = per["reading"] + ". " + xf.PRIOR_EXPOSURE_CLAUSE
    out = {
        **head,
        "population": {
            "genes": len(rows),
            "pairs": sum(len(r["D"]) for r in rows),
            "excluded_by_reason": dict(why),
            "rna_files": rna_counts,
        },
        "per_form": {FORM: per},
        "outcome": per["outcome"],
        "reading": reading,
        "head_to_heads_descriptive": h2h,
        "prior_exposure": xf.PRIOR_EXPOSURE,
        "falsifier": xf.FALSIFIER,
        "inconclusive": xf.INCONCLUSIVE,
        "limitations": list(xf.LIMITATIONS),
        "genes_table": _table(rows),
        "seconds": round(time.time() - t0, 1),
    }
    print(json.dumps({k: per[k] for k in ("outcome", "scale", "scale_ci") if k in per}, indent=1))
    print(reading)
    print("saved", save_result(RESULT, out, manifest=manifest(entries, rna_read=True), compact=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
