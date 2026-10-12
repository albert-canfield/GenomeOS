# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-hetblock, step 4: the run. Does scalehet's heterogeneity survive a locus-block null?

    uv run python scripts/hetblock_run.py

It refuses to start unless data/results/hetblock_registration.json exists, every BOUND field in it
equals the committed module's, and every pinned input's sha256 still matches the file on disk; and it
refuses to overwrite its own result. Then, exactly once: the observed Q and R on scalehet's own admitted
chromosomes (which must equal scalehet's committed values), the governing 1-Mb LOCUS-BLOCK null at the
registered count and seed, the registered coarser 5-Mb block null, the element-sharing component null,
the two one-sided p values, the registered outcome through `hetblock.decide_block`, and the registered
sentence for that outcome with its numbers filled and the registered 5-Mb clause chosen by the rule.

NO REGIONAL MODEL IS FITTED WHATEVER THE OUTCOME. No file of scalehet's, magcal's, scalesplit's,
freshrna's or multrna's is edited, no request is sent and no key is read.
"""

from __future__ import annotations

import json
import os
import sys
import time
from datetime import date
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import hetblock as hb  # noqa: E402
from genomeos.attribution import scalehet as sh  # noqa: E402
from genomeos.results import save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hetblock_harness import (  # noqa: E402
    MR_RESULT,
    OWN_CODE,
    READERS,
    SCALEHET_RESULT,
    SOURCE_FILES,
    admitted_parts,
    load_rows,
    with_gene_ids,
)

RESULT = "hetblock_result"
OUT = "data/results/hetblock_result.json"
REGISTRATION = "data/results/hetblock_registration.json"
HARNESS_RESULT = "data/results/hetblock_harness.json"


def fmt(x: float, places: int = 4) -> str:
    return f"{x:+.{places}f}"


def refuse_unless_registered() -> dict[str, Any]:
    """The registration must exist, bind this module word for word, and pin these inputs."""
    if Path(OUT).exists():
        raise SystemExit(f"refused: {OUT} already exists; the run is not re-run (STOP_RULE)")
    if not Path(REGISTRATION).exists():
        raise SystemExit(f"refused: {REGISTRATION} does not exist; register and commit it first")
    if not Path(HARNESS_RESULT).exists():
        raise SystemExit(f"refused: {HARNESS_RESULT} does not exist; run and commit the harness first")
    harness = json.loads(Path(HARNESS_RESULT).read_text())
    if harness["outcome"] != "HARNESS_REPRODUCES":
        raise SystemExit("refused: the committed harness did not reproduce lane-scalehet's figures")
    reg = json.loads(Path(REGISTRATION).read_text())
    drift = [f for f, v in hb.registered_fields().items() if reg.get(f) != v]
    if drift:
        raise SystemExit(f"refused: the registration differs from the module in {drift}")
    moved = [p for p, want in reg["pinned_sha256"].items() if mf.input_entry(p)["sha256"] != want]
    if moved:
        raise SystemExit(f"refused: a pinned input moved since the registration: {moved}")
    return reg


def placed_rows() -> tuple[list[dict[str, Any]], dict[str, Any], list[str], list[int], list[int]]:
    """scalehet's rows with each gene's start attached, and scalehet's own observed statistics."""
    rows = with_gene_ids(load_rows())
    starts = hb.gene_starts()
    missing = hb.with_starts(rows, starts)
    if missing:
        raise SystemExit(f"refused: {len(missing)} genes have no start, the first being {missing[0]}")
    rows = [{**r, "start": starts[str(r["gene_id"]).split(".")[0]]} for r in rows]
    obs, adm, sizes, keep = admitted_parts(rows)
    if obs.get("degenerate"):
        raise SystemExit(f"refused: the observed statistics are degenerate: {obs['degenerate']}")
    return rows, obs, adm, sizes, keep


def observed_matches_scalehet(obs: dict[str, Any]) -> list[str]:
    """The observed side is scalehet's, unchanged: this lane changes the NULL and nothing else."""
    problems = []
    for key, name in (("q", "q_observed"), ("r", "r_observed"), ("global_scale", "global_scale")):
        if abs(float(obs[key]) - float(hb.SCALEHET_COMMITTED[name])) > hb.TOLERANCE:
            problems.append(
                f"the observed {name} is {obs[key]!r}, not scalehet's committed"
                f" {hb.SCALEHET_COMMITTED[name]!r}: the observed side must be unchanged"
            )
    return problems


def summarise(null: dict[str, Any]) -> dict[str, Any]:
    """INTERVAL: a null's median, 98.33% percentile range, minimum and maximum."""
    q, r = sorted(null["q"]), sorted(null["r"])
    q_lo, q_hi = sh.percentile_range(null["q"])
    r_lo, r_hi = sh.percentile_range(null["r"])
    return {
        "q": {
            "median": q[len(q) // 2],
            "percentile_98_33": [q_lo, q_hi],
            "min": q[0],
            "max": q[-1],
        },
        "r": {
            "median": r[len(r) // 2],
            "percentile_98_33": [r_lo, r_hi],
            "min": r[0],
            "max": r[-1],
        },
        "blocks": null["blocks"],
        "split_blocks_per_draw": null["split_blocks_per_draw"],
        "count": null["count"],
        "seed": null["seed"],
    }


def one_null(
    dps: list[float], pps: list[float], blocks: list[list[int]], sizes: list[int], obs: dict[str, Any]
) -> tuple[dict[str, Any], dict[str, Any]]:
    """A block null, its summary and its two one-sided p values, through the registered functions."""
    null = hb.block_permutation_null(
        dps, pps, blocks, sizes, hb.PERMUTATIONS["count"], hb.PERMUTATIONS["seed"]
    )
    if null["degenerate"]:
        return null, {"degenerate": null["degenerate"]}
    s = summarise(null)
    s["p_q"] = sh.permutation_p(obs["q"], null["q"])
    s["p_r"] = sh.permutation_p(obs["r"], null["r"])
    s["p_q_below_alpha"] = bool(s["p_q"] <= hb.CONTROL_ALPHA)
    s["p_q_below_the_stricter_alpha"] = bool(s["p_q"] <= hb.STRICTER_ALPHA)
    s["p_r_below_alpha"] = bool(s["p_r"] <= hb.CONTROL_ALPHA)
    s["p_r_below_the_stricter_alpha"] = bool(s["p_r"] <= hb.STRICTER_ALPHA)
    s["degenerate"] = None
    return null, s


def per_chromosome_table(
    obs: dict[str, Any], adm: list[str], sizes: list[int], null: dict[str, Any]
) -> list[dict[str, Any]]:
    """DESCRIPTIVES: each admitted chromosome's fitted scale -- scalehet's values -- and its block band.

    The bands carry NO correction over the 23 chromosomes and are NOT a test of any one of them. No
    chromosome is read as a verdict about itself and no per-region scale is fitted.
    """
    out = []
    for i, c in enumerate(adm):
        lo, hi = sh.percentile_range(null["per_group"][i])
        a = obs["scales"][i]
        out.append(
            {
                "chrom": c,
                "genes": sizes[i],
                "scale": a,
                "deviation_from_global": a - obs["global_scale"],
                "weighted_squared_deviation": obs["pp"][i] * (a - obs["global_scale"]) ** 2,
                "block_null_scale_98_33_range": [lo, hi],
                "inside_its_own_block_null_band": bool(lo <= a <= hi),
                "not_a_verdict": "a descriptive band with no correction over the 23 chromosomes",
            }
        )
    return out


def context() -> dict[str, Any]:
    """The committed readings this lane re-tests, quoted verbatim and not re-decided."""
    s = json.loads(Path(SCALEHET_RESULT).read_text())
    m = json.loads(Path(MR_RESULT).read_text())
    return {
        "scalehet": {
            "source": SCALEHET_RESULT,
            "outcome": s["outcome"],
            "reading": s["reading"],
            "its_own_limitation_that_this_lane_prices": s["limitations"][2],
        },
        "seven_h": {"source": MR_RESULT, "outcome": m["outcome"], "reading": m["reading"]},
        "note": (
            "quoted and not re-decided. Nothing in this lane revises M2h's one scale +0.2146 or its about"
            " 2% of the variance, which stand as a GLOBAL AVERAGE either way, and a scale remains a"
            " MAGNITUDE CORRECTION and never an explanation of expression"
        ),
    }


def manifest(entries: list[dict[str, Any]], adm: list[str]) -> dict[str, Any]:
    return {
        "sources": [
            {"accession": "data/results/multrna_result.json (7h, lane-multrna)", "version": "pinned"},
            {
                "accession": "data/results/scalehet_result.json (lane-scalehet)",
                "version": "pinned; the reading this lane re-tests under a block null",
            },
            {
                "accession": "data/results/hetblock_registration.json (this lane, step 3)",
                "version": "pinned; committed before any statistic of this lane was read",
            },
            {
                "accession": "the 24 compiled non-coding programs and the GENCODE gene table",
                "version": "pinned by sha256; addresses only -- element sharing and gene starts",
            },
            {
                "accession": "ENCODE polyA plus RNA-seq gene quantifications, V29, 13 reader biosamples",
                "version": "read by 7d, 7h and M2h; this lane reads their committed rows, not the files",
            },
        ],
        "inputs": entries,
        "input_count": len(entries),
        "readers": READERS,
        "reader_count": len(READERS),
        "assembly": "GRCh38; genes as 7h matched them, by gene_id",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "form": hb.FORM,
            "grouping": "by chromosome, scalehet's grouping unchanged",
            "exchangeable_unit": "a 1-Mb locus block, (chromosome, gene start // 1,000,000)",
            "robust_unit": "a 5-Mb locus block, registered and deciding nothing",
            "statistic": "Q = sum over admitted chromosomes of pp_c * (a_c - a_g)^2; R = max a_c - min a_c",
            "null": (
                "the chromosome labels permuted among the genes BY WHOLE LOCUS BLOCKS, every"
                " chromosome's gene count held exactly at its observed value"
            ),
            "min_genes": hb.MIN_GENES,
            "min_chromosomes": hb.MIN_CHROMOSOMES,
            "min_blocks": hb.MIN_BLOCKS,
            "max_block_share": hb.MAX_BLOCK_SHARE,
            "admitted_chromosomes": adm,
            "permutations": hb.PERMUTATIONS,
            "control_alpha": hb.CONTROL_ALPHA,
            "stricter_alpha": hb.STRICTER_ALPHA,
            "form_interval": list(hb.FORM_INTERVAL),
            "tolerance": hb.TOLERANCE,
            "free_parameters": (
                "one closed-form scale per admitted chromosome and one global scale; the null refits"
                " every one of them on each draw, so nothing fitted here is unpriced"
            ),
        },
        "exclusions": [
            f"chromosomes below MIN_GENES = {hb.MIN_GENES} genes: chrY at 13 genes, scalehet's own"
            " exclusion, carried unchanged"
        ],
        "partitions": {
            "none": (
                "no training and held-out split: the quantity is a dispersion of fitted scales and the"
                " null prices the fitting exactly (NO_SPLIT, carried from scalehet)"
            )
        },
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    t0 = time.time()
    reg = refuse_unless_registered()
    rows, obs, adm, sizes, keep = placed_rows()
    problems = observed_matches_scalehet(obs)
    if problems:
        raise SystemExit("refused: " + "; ".join(problems))
    dps, pps = sh.gene_sums(rows)
    d, p = [dps[i] for i in keep], [pps[i] for i in keep]
    kept = [rows[i] for i in keep]

    blocks_1mb = hb.locus_blocks(kept, hb.BLOCK_SPAN)
    blocks_5mb = hb.locus_blocks(kept, hb.ROBUST_SPAN)
    element_components = hb.element_blocks(
        kept, {}
    )  # the harness committed that NO element drives two of these genes
    sizes_1mb = hb.block_sizes(blocks_1mb)
    sizes_5mb = hb.block_sizes(blocks_5mb)
    sizes_el = hb.block_sizes(element_components)

    degen = hb.degeneracy(len(adm), sizes_1mb["blocks"], sizes_1mb["largest_share"])
    null_1mb, s1 = one_null(d, p, blocks_1mb, sizes, obs)
    _null_5mb, s5 = one_null(d, p, blocks_5mb, sizes, obs)
    _null_el, sel = one_null(d, p, element_components, sizes, obs)
    degen = degen or s1.get("degenerate")

    out: dict[str, Any]
    if degen:
        numbers = {"why": degen}
        outcome = hb.decide_block(1.0, 1.0, degen)
        out = {
            "outcome": outcome,
            "degenerate": degen,
            "reading_numbers": numbers,
            "reading": hb.full_reading(outcome, **numbers),
        }
    else:
        outcome = hb.decide_block(s1["p_q"], s1["p_r"], None)
        robust = hb.robust_clause(
            s1["p_q"],
            s5["p_q"],
            q5=f"{obs['q']:.3f}",
            null_q5_median=f"{s5['q']['median']:.4f}",
            null_q5_lo=f"{s5['q']['percentile_98_33'][0]:.4f}",
            null_q5_hi=f"{s5['q']['percentile_98_33'][1]:.4f}",
            n_blocks5=f"{sizes_5mb['blocks']:,}",
        )
        numbers = {
            "unit": "1-Mb LOCUS BLOCK",
            "n_blocks": f"{sizes_1mb['blocks']:,}",
            "genes": f"{len(keep):,}",
            "chroms": len(adm),
            "global_scale": fmt(obs["global_scale"]),
            "q": f"{obs['q']:.3f}",
            "null_q_median": f"{s1['q']['median']:.4f}",
            "null_q_lo": f"{s1['q']['percentile_98_33'][0]:.4f}",
            "null_q_hi": f"{s1['q']['percentile_98_33'][1]:.4f}",
            "p_q": f"{s1['p_q']:.6f}",
            "alpha": hb.CONTROL_ALPHA,
            "stricter": hb.STRICTER_ALPHA,
            "r": fmt(obs["r"]),
            "null_r_median": fmt(s1["r"]["median"]),
            "null_r_lo": fmt(s1["r"]["percentile_98_33"][0]),
            "null_r_hi": fmt(s1["r"]["percentile_98_33"][1]),
            "p_r": f"{s1['p_r']:.6f}",
            "robust": robust,
        }
        out = {
            "outcome": outcome,
            "degenerate": None,
            "reading_numbers": numbers,
            "reading": hb.full_reading(outcome, **numbers),
            "per_chromosome": per_chromosome_table(obs, adm, sizes, null_1mb),
        }
    out.update(
        {
            "lane": hb.LANE,
            "step": "4 of 4: the run, once, after the registration was committed",
            "registration": REGISTRATION,
            "registration_date": reg["registered"],
            "harness": HARNESS_RESULT,
            "governing_unit": "locus_block_1mb",
            "observed": {
                "q": obs["q"],
                "r": obs["r"],
                "global_scale": obs["global_scale"],
                "note": (
                    "THE OBSERVED SIDE IS LANE-SCALEHET'S, UNCHANGED and checked equal to its committed"
                    " values within TOLERANCE. This lane changes the NULL and nothing else"
                ),
            },
            "scalehet_committed": dict(hb.SCALEHET_COMMITTED),
            "admitted": adm,
            "admitted_count": len(adm),
            "genes_admitted": len(keep),
            "pairs_admitted": len(keep) * len(hb.CELLS),
            "block_nulls": {
                "locus_block_1mb_GOVERNING": {**s1, "block_sizes": sizes_1mb},
                "locus_block_5mb_robustness_deciding_nothing": {**s5, "block_sizes": sizes_5mb},
                "element_sharing_components_which_price_nothing_here": {
                    **sel,
                    "block_sizes": sizes_el,
                    "note": (
                        "every component is a single gene, because no compiled element drives two of"
                        " these genes, so this null IS scalehet's gene-level null; it is reported to make"
                        " that a committed number rather than an assertion"
                    ),
                },
            },
            "null_q_1mb": null_1mb["q"],
            "null_r_1mb": null_1mb["r"],
            "how_much_the_block_unit_widens_the_null": {
                "gene_null_median_q": hb.SCALEHET_COMMITTED["null_q_median"],
                "block_null_median_q": s1["q"]["median"] if not degen else None,
                "ratio": (s1["q"]["median"] / hb.SCALEHET_COMMITTED["null_q_median"]) if not degen else None,
                "gene_level_p_q_scalehet_committed": hb.SCALEHET_COMMITTED["p_q"],
                "block_level_p_q": s1.get("p_q"),
                "note": (
                    "A DESCRIPTIVE AND NOT A VERDICT (DESCRIPTIVES, second paragraph): the quantity"
                    " scalehet could not price. It decides nothing"
                ),
            },
            "null_global_scale_is_invariant": (
                "by construction: permuting the labels by whole blocks redistributes the same genes, so"
                " the total sums and therefore the global fitted scale are the same in every draw,"
                " to within the order of the floating-point summation (TOLERANCE) and not bit-equal"
            ),
            "no_regional_model_was_fitted": (
                "none, whatever the outcome: no per-region scale is fitted, no cause is named, and no"
                " chromosome or block is read as a verdict about itself. That variation exists, or is not"
                " established, is the finding and this lane stops there"
            ),
            "question": hb.QUESTION,
            "what_is_not_asked": hb.WHAT_IS_NOT_ASKED,
            "block_unit": hb.BLOCK_UNIT,
            "block_unit_why_not_elements": hb.BLOCK_UNIT_WHY_NOT_ELEMENTS,
            "block_selection_rule": hb.BLOCK_SELECTION_RULE,
            "block_robustness": hb.BLOCK_ROBUSTNESS,
            "grouping": hb.GROUPING,
            "statistic": hb.STATISTIC,
            "null_description": hb.NULL,
            "permutations": hb.PERMUTATIONS,
            "alpha": hb.ALPHA,
            "interval": hb.INTERVAL,
            "which_governs": hb.WHICH_GOVERNS,
            "not_vacuous": hb.NOT_VACUOUS,
            "decision": hb.DECISION,
            "falsifier": hb.FALSIFIER,
            "prior_exposure": list(hb.PRIOR_EXPOSURE),
            "prior_exposure_clause": hb.PRIOR_EXPOSURE_CLAUSE,
            "stop_rule": hb.STOP_RULE,
            "limitations": list(hb.LIMITATIONS),
            "descriptives": list(hb.DESCRIPTIVES),
            "context_verbatim": context(),
            "cells": list(hb.CELLS),
            "readers": READERS,
            "alphagenome_requests": 0,
            "money": "none: no request is sent and no key is read",
            "result": RESULT,
            "date": date.today().isoformat(),
            "seconds": round(time.time() - t0, 1),
        }
    )
    entries = [mf.input_entry(p, partition="pinned by sha256") for p in (*SOURCE_FILES, REGISTRATION)]
    print("saved", save_result(RESULT, out, manifest=manifest(entries, adm)))
    print(out["outcome"])
    print(out["reading"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
