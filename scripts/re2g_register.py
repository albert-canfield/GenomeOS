# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register the paired comparison against ENCODE-rE2G before any comparator score is read
(data/results/re2g_registration.json).

    uv run --frozen python scripts/re2g_register.py

Reads the held-out benchmark table that is already in the cache, for the pair, positive and weight
counts that name each population, and the paper's Supplementary Tables zip, for the published figure the
gate is set against. It reads no predictor's score: the comparator scores do not exist at the time this
runs, which is the point of running it first.

Everything it writes comes from `genomeos.attribution.re2g`, so the registration and the code that would
apply it cannot drift apart. The pairs are read through `crispri.load`, so the label and the weight are
parsed by the same rule every other result in this module uses. docs/ATTRIBUTION.md, "The paired
comparison against ENCODE-rE2G, registered before any comparator score was read".
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri, re2g  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "re2g_registration"
HELDOUT = crispri.KNOWLEDGE / crispri.HELDOUT
SUPPLEMENT = Path("data/cache/re2g/41586_2026_10781_MOESM3_ESM.zip")
SUPPLEMENT_URL = (
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-026-10781-4/MediaObjects/"
    "41586_2026_10781_MOESM3_ESM.zip"
)


def strata(pairs: list[crispri.Pair]) -> dict[str, dict[str, Any]]:
    """Pairs, positives and weighted positives per cell type, by the module's own label and weight rule."""
    counts: dict[str, list[float]] = defaultdict(lambda: [0, 0, 0.0])
    for p in pairs:
        entry = counts[p.cell]
        entry[0] += 1
        if p.regulated:
            entry[1] += 1
            entry[2] += p.weight
    return {
        cell: {"pairs": int(v[0]), "positives": int(v[1]), "weighted_positives": round(v[2], 2)}
        for cell, v in sorted(counts.items(), key=lambda kv: -kv[1][0])
    }


def main() -> None:
    if not HELDOUT.exists():
        raise SystemExit(f"the held-out benchmark table is not in the cache: {HELDOUT}")
    if not SUPPLEMENT.exists():
        raise SystemExit(f"the paper's supplementary tables are not in the cache: {SUPPLEMENT}")

    pairs = crispri.load(crispri.HELDOUT)
    per_cell = strata(pairs)
    # Summed over the pairs, not over the per-cell figures: adding values already rounded to 2 dp loses
    # a hundredth of the published weighted total (157.38 against 157.39).
    pooled = {
        "pairs": len(pairs),
        "positives": sum(1 for p in pairs if p.regulated),
        "weighted_positives": round(sum(p.weight for p in pairs if p.regulated), 2),
    }

    payload: dict[str, Any] = {
        "status": "registered before any comparator score was read; the comparison is NOT run",
        "lane": "lane-re2g",
        "claim": re2g.CLAIM,
        "endpoint": re2g.ENDPOINT,
        "comparator": {
            "preference_order": list(re2g.COMPARATOR_PREFERENCE),
            "preference_1_per_pair_score_in_the_benchmark_table": re2g.COMPARATOR_IN_TABLE,
            "governs": "preference 2, the ENCODE portal's rE2G predictions per benchmark biosample",
        },
        "comparator_files": {k: dict(v) for k, v in re2g.COMPARATOR_FILES.items()},
        "comparator_not_obtained": {
            "gate_run": False,
            "reason": dict(re2g.BLOCKER),
            "bytes_required": re2g.COMPARATOR_BYTES_REQUIRED,
            "k562_only_attainable_range": {
                "pooled_weighted_auprc": list(re2g.K562_ONLY_RANGE),
                "population": "all 4,378 held-out pairs, pooled over the five cell types",
                "ends": "no skill inside K562; perfect separation inside K562",
                "brackets_the_target": True,
                "means": (
                    "the one prediction file that fits the download budget could land within the "
                    "tolerance of the published figure without carrying the published scores, so passing "
                    "the gate with it would establish nothing. Real labels and real weights, synthetic "
                    "scores; reproduced by scripts/re2g_k562_only_range.py"
                ),
            },
            "bytes_authorised": re2g.DOWNLOAD_BUDGET_BYTES,
            "sha256": (
                "not recorded: no prediction file was downloaded, so no model version or file hash is "
                "claimed here beyond the accessions above"
            ),
        },
        "join_rule": dict(re2g.JOIN_RULE),
        "join_rule_source": {
            "repository": "EngreitzLab/CRISPR_comparison",
            "tag": re2g.BENCHMARK_TAG,
            "commit": re2g.BENCHMARK_COMMIT,
            "files": (
                "workflow/scripts/crisprComparisonMergeFunctions.R (combineSingleExptPred steps 1 to 3, "
                "fillMissingPredictions) and workflow/scripts/createPredConfig.R (aggregate_function sum, "
                "fill_value 0)"
            ),
            "pinned_not_branch": (
                "read at the tag the paper's data availability names, not at main. The merge functions are "
                "byte-identical at both (sha256 "
                "fee907d41ff9a72740dd9cf57ace795536fa9adb67f8ff462c94777a6fea916c, 506 lines) and the "
                "defaults match, so the registered rule is the one the published comparison used"
            ),
            "when": "the rule is the benchmark's own, read from its source before any join was attempted",
        },
        "gate": re2g.GATE,
        "gate_target": {
            "weighted_auprc": re2g.GATE_TARGET,
            "interval": list(re2g.GATE_INTERVAL),
            "tolerance": re2g.GATE_TOLERANCE,
            "source": (
                "Supplementary Table 3, sheet 'Held-out benchmarks', row ENCODE-rE2G / Weighted AUPRC, "
                "read from the cached zip"
            ),
            "population": (
                "all 4,378 held-out pairs, pooled over the five cell types. The sheet carries six "
                "predictors times three metrics for one 'Held-out' dataset and no per-cell-type held-out "
                "AUPRC, so the pooled figure is the only gate value the source gives"
            ),
        },
        "populations": {k: dict(v) for k, v in re2g.POPULATIONS.items()},
        "populations_as_counted": {
            "per_cell_type": per_cell,
            "pooled": pooled,
            "counted_from": (
                "the held-out table's own CellType, Regulated and direct_vs_indirect_negative columns, "
                "through crispri.load, before any comparator score existed"
            ),
        },
        "feature_unavailable": {
            "strata": list(re2g.FEATURE_UNAVAILABLE),
            "reported_as": crispri.UNAVAILABLE_GAIN,
            "never": "a number, in either direction",
        },
        "statistic": re2g.STATISTIC,
        "statistic_parameters": {
            "draws_requested": re2g.DRAWS,
            "seed": re2g.SEED,
            "cluster": re2g.CLUSTER,
            "provenance_reported": ["clusters", "draws_requested", "draws_dropped", "resamples"],
        },
        "reading": {
            "lower_bound_above_zero": re2g.READS_BETTER,
            "interval_covers_zero": re2g.NO_DIFFERENCE,
            "upper_bound_below_zero": re2g.READS_WORSE,
            "no_other_wording": "these three strings are the only wording this comparison is reported in",
        },
        "independence_and_exposure": dict(re2g.INDEPENDENCE),
        "falsifier": re2g.FALSIFIER,
        "code": {
            "module": "genomeos/attribution/re2g.py",
            "tests": "tests/test_re2g.py (synthetic inputs only; no comparator score is read)",
            "register": "scripts/re2g_register.py",
            "k562_only_range": "scripts/re2g_k562_only_range.py",
        },
        "alphagenome_requests": 0,
    }

    inputs = [
        mf.input_entry(HELDOUT, partition="held-out", url=crispri.BASE_URL + crispri.HELDOUT),
        mf.input_entry(SUPPLEMENT, partition=None, url=SUPPLEMENT_URL),
    ]
    payload["frozen_inputs"] = [{k: e[k] for k in ("path", "sha256", "bytes") if k in e} for e in inputs]
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "Gschwind et al. 2026, Nature, doi:10.1038/s41586-026-10781-4, "
                "Supplementary Tables 3 and 12",
                "version": "version of record, supplementary zip 41586_2026_10781_MOESM3_ESM.zip",
            },
            {
                "accession": "EngreitzLab/CRISPR_comparison, "
                "EPCrisprBenchmark_combined_data.heldout_5_cell_types.GRCh38.tsv.gz",
                "version": "v1.0.0, the release the paper's data availability names",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "gate_target": re2g.GATE_TARGET,
            "gate_tolerance": re2g.GATE_TOLERANCE,
            "draws_requested": re2g.DRAWS,
            "seed": re2g.SEED,
            "cluster": re2g.CLUSTER,
            "aggregate_function": "sum",
            "fill_value": 0,
        },
        "exclusions": [
            "no pair is excluded: the benchmark's published pair set and labels are used whole",
            f"no number is reported for {', '.join(re2g.FEATURE_UNAVAILABLE)}: the deletion feature does "
            "not exist there",
        ],
        "partitions": {
            "primary": "held-out K562",
            "secondary": "all held-out pairs, pooled",
            "gm12878": "held-out GM12878, its own stratum",
        },
    }
    path = save_result(RESULT, payload)
    print(f"{RESULT}: {path}")
    print(
        f"  pooled {pooled['pairs']} pairs, {pooled['positives']} positives, "
        f"{pooled['weighted_positives']} weighted"
    )
    print(f"  gate {re2g.GATE_TARGET} +/- {re2g.GATE_TOLERANCE} on the pooled population")


if __name__ == "__main__":
    main()
