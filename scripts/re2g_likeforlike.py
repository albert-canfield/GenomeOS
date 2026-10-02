# SPDX-License-Identifier: AGPL-3.0-or-later
"""The like-for-like paired comparison: `dnase + distance + deletion` against ENCODE-rE2G
(data/results/re2g_likeforlike.json).

    uv run --frozen python scripts/re2g_likeforlike.py

Applies `re2g.SECOND_REGISTRATION`, committed before any figure of this comparison existed. Its reason
for existing: all 190 held-out positives carry H3K27ac at the tested element while 1,438 non-regulated
pairs do not, so an H3K27ac-reading model is handed a separation the comparator's DNase-only held-out
model cannot use. Here our activity term is DNase, matching the comparator's input.

Nothing of the comparator changes. The same five ENCODE portal files, the same join rule, the same fill
convention, the same per-pair scores that passed the gate. The like-for-like is achieved by changing our
features and never theirs; no refit, reweight or adjustment is applied to ENCODE-rE2G at any point.

Registered population: the 1,918 held-out K562 pairs. The pooled population is descriptive context.
HCT116, Jurkat and WTC11 report the unavailable-feature refusal, never a number. 0 model requests.
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri, re2g  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "re2g_likeforlike"
MODEL = "dnase + distance + deletion"
BASELINE = "dnase + distance"
#: The only population `re2g.SECOND_REGISTRATION` registers. Every other population in this result is
#: descriptive context and carries NO reading word: a registered reading on a descriptive arm gets lifted
#: as a finding, which is the overstatement this whole lane exists to have withdrawn.
REGISTERED_POPULATIONS = ("primary_k562",)
DESCRIPTIVE = (
    "no registered reading: descriptive context in re2g.SECOND_REGISTRATION, not a registered population"
)
#: Why the baseline arm's reading may stand: it is in the registration code, with its own reading, before
#: any figure of it existed. It is a pre-specified secondary arm, not a post-hoc one.
BASELINE_LABEL = (
    "pre-specified secondary arm, registered in code at bdba855 within the registered K562 population, "
    "before any figure of it existed"
)


def _paired_module() -> Any:
    spec = importlib.util.spec_from_file_location(
        "re2g_paired", Path(__file__).resolve().parent / "re2g_paired.py"
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def dnase_scores(
    training: list[crispri.Pair], heldout: list[crispri.Pair]
) -> tuple[dict[str, list[float]], dict[str, Any]]:
    """The DNase variant's per-pair scores, by the same recipe and the same frozen-weights discipline.

    The two DNase columns are set exactly as `crispri.dnase_only_diagnostic` sets them, so this is the
    same feature definition the committed diagnostic used and not a new one.
    """
    crispri.annotate(training + heldout, crispri.DeletionTable(), crispri.ElementCache())
    for p in training + heldout:
        p.features["covered"] = float(p.covered)
        d = max(crispri.MIN_DISTANCE, p.distance)
        p.features["log_dnase"] = math.log1p(max(p.dhs, 0.0))
        p.features["dnase_over_distance"] = p.features["log_dnase"] - math.log(d)
    covered = [p for p in training if p.covered]
    weights = {
        name: crispri.logistic_fit(crispri.matrix(covered, cols), [p.regulated for p in covered])
        for name, cols in crispri.DNASE_FEATURES.items()
    }
    scores = {
        name: crispri.logistic_score(weights[name], crispri.matrix(heldout, cols))
        for name, cols in crispri.DNASE_FEATURES.items()
    }
    return scores, {
        "features": {k: list(v) for k, v in crispri.DNASE_FEATURES.items()},
        "fitted_on": f"{len(covered)} covered training pairs of {len(training)}",
        "frozen_before_scoring_heldout": True,
        "feature_definition": (
            "log_dnase = log1p(DHS.RPM) and dnase_over_distance = log_dnase - log(max(1000, distance)), "
            "exactly as crispri.dnase_only_diagnostic sets them; H3K27ac is not read by any column"
        ),
    }


def main() -> None:
    paired = _paired_module()
    missing = [c for c in re2g.COMPARATOR_FILES if not paired.prediction_path(c).exists()]
    if missing:
        raise SystemExit(f"prediction files not in the cache for: {', '.join(missing)}")

    training, heldout = crispri.load(crispri.TRAINING), crispri.load(crispri.HELDOUT)
    ours_all, provenance = dnase_scores(training, heldout)
    ours, base = ours_all[MODEL], ours_all[BASELINE]

    theirs, join_report = paired.comparator_scores(heldout, "sum")
    labels, weights = [p.regulated for p in heldout], [p.weight for p in heldout]
    reproduced = crispri.benchmark_auprc(theirs, labels, weights)
    verdict = re2g.gate(reproduced)
    if not verdict["passed"]:
        raise SystemExit(f"the comparator scores no longer pass the gate: {verdict['reason']}")

    populations = {
        "primary_k562": [i for i, p in enumerate(heldout) if p.cell == "K562"],
        "pooled_descriptive": list(range(len(heldout))),
        "gm12878_descriptive": [i for i, p in enumerate(heldout) if p.cell == "GM12878"],
    }
    out: dict[str, Any] = {}
    for name, idx in populations.items():
        rows = [heldout[i] for i in idx]
        delta = re2g.paired_delta([ours[i] for i in idx], [theirs[i] for i in idx], rows)
        baseline = re2g.paired_delta([base[i] for i in idx], [theirs[i] for i in idx], rows)
        delta["population"] = {
            "name": name,
            "pairs": len(rows),
            "positives": sum(1 for p in rows if p.regulated),
            "weighted_positives": round(sum(p.weight for p in rows if p.regulated), 2),
        }
        delta["ours_auprc"] = crispri.bench_metrics([ours[i] for i in idx], rows, weighted=True)["auprc"]
        delta["re2g_auprc"] = crispri.bench_metrics([theirs[i] for i in idx], rows, weighted=True)["auprc"]
        registered = name in REGISTERED_POPULATIONS
        delta["dnase_plus_distance_alone"] = {
            "auprc": crispri.bench_metrics([base[i] for i in idx], rows, weighted=True)["auprc"],
            "delta_auprc": baseline["delta_auprc"],
            "ci95": baseline["ci95"],
            "reading": baseline["reading"] if registered else None,
            "label": BASELINE_LABEL if registered else DESCRIPTIVE,
        }
        # A descriptive population keeps its point estimate and its interval and loses the reading word.
        delta["registered_population"] = registered
        if not registered:
            delta["reading"] = None
            delta["descriptive"] = DESCRIPTIVE
        delta["deletion_features_share"] = round(
            (delta["delta_auprc"] or 0) - (baseline["delta_auprc"] or 0), 4
        )
        out[name] = delta
        print(
            f"  {name}: {MODEL} vs rE2G {delta['delta_auprc']} ci {delta['ci95']} "
            f"clusters {delta['clusters']} draws {delta['resamples']}/{delta['draws_requested']} "
            f"-> {delta['reading']}"
        )
        print(f"      {BASELINE} alone vs rE2G {baseline['delta_auprc']} ci {baseline['ci95']}")
    for cell in re2g.FEATURE_UNAVAILABLE:
        out[cell.lower()] = re2g.delta_where_available(cell, lambda: {})
        print(f"  {cell}: feature unavailable, no number reported")

    primary = out["primary_k562"]
    payload: dict[str, Any] = {
        # First key on purpose: this is the lane's binding outcome and it is a negative one.
        "read_this_first": (
            f"THE BINDING RESULT OF THIS LANE. On the registered primary population, the 1,918 held-out "
            f"K562 pairs, the frozen 'dnase + distance + deletion' model against ENCODE-rE2G reads "
            f'"{primary["reading"]}": {primary["delta_auprc"]:+}, interval {primary["ci95"]}, '
            f"{primary['clusters']} chromosome clusters, {primary['resamples']} of "
            f"{primary['draws_requested']} draws kept. This is the like-for-like comparison, so by its "
            "pre-registered falsifier README MAY NOT SAY THE DELETION MODEL RANKS BETTER THAN "
            "ENCODE-rE2G, and the H3K27ac comparison in re2g_paired.json may not be cited alone. The "
            "H3K27ac comparison reads 'ranks better than ENCODE-rE2G on these pairs' on the same pairs "
            "with a lower bound of +0.0022; the difference between the two is H3K27ac, which all 190 "
            "held-out positives carry and the comparator's held-out model does not read"
        ),
        "status": (
            "the like-for-like comparison: our feature set matched to the comparator's input, so DNase "
            "is shared and H3K27ac is read by neither model. This is the comparison that may be cited "
            "about ENCODE-rE2G; the H3K27ac comparison in re2g_paired.json is reported with its "
            "qualifiers and is never cited alone"
        ),
        "lane": "lane-re2g",
        "registration": "genomeos/attribution/re2g.py SECOND_REGISTRATION, committed before any figure "
        "of this comparison was computed",
        "registered_before_computing": dict(re2g.SECOND_REGISTRATION),
        "our_model": provenance,
        "comparator_unchanged": {
            "files": {c: v["file"] for c, v in sorted(re2g.COMPARATOR_FILES.items())},
            "join": join_report["_total"],
            "gate": verdict,
            "nothing_adjusted": (
                "the comparator's per-pair scores are the same ones that passed the gate; no refit, "
                "reweight or adjustment was applied to them here or anywhere in this lane"
            ),
        },
        "paired_delta": out,
        # The registered population only. A descriptive arm's figure lives in paired_delta and is not a
        # reading, so it cannot be quoted from here as one.
        "reading": {
            name: out[name]["reading"] for name in REGISTERED_POPULATIONS if out.get(name, {}).get("reading")
        },
        "registered_populations": list(REGISTERED_POPULATIONS),
        "descriptive_populations": [n for n in populations if n not in REGISTERED_POPULATIONS],
        "why_descriptive_arms_carry_no_reading": (
            "a registered reading word on a descriptive arm is lifted as a finding. The pooled arm of "
            "this comparison previously carried 'ranks better than ENCODE-rE2G on these pairs' three "
            "lines below a primary reading 'no difference detected', which is the overstatement this "
            "lane exists to have withdrawn. Only the output labels changed when this was fixed: no "
            "number moved, no threshold moved, and the analysis is the same one"
        ),
        "reading_must_travel_with": {
            "lower_bound": f"the K562 primary interval's lower bound is {primary['ci95'][0]}, unrounded",
            "the_comparator_is_a_reconstruction": (
                f"ENCODE-rE2G as reconstructed from the portal prediction files, reproducing its "
                f"published pooled weighted AUPRC {verdict['miss']} low ({verdict['reproduced']} "
                f"against {re2g.GATE_TARGET:.4f})"
            ),
            "dnase_is_still_shared": (
                "DNase is read by both models, ours from the benchmark table's DHS.RPM and theirs "
                "through the portal predictions, so this is a like-for-like comparison of inputs and "
                "not a comparison on independent evidence"
            ),
            "h3k27ac_is_read_by_neither": (
                "which is the point of this comparison: the benchmark's total positive selection on "
                "H3K27ac can no longer favour our side"
            ),
            "still_a_reused_benchmark": (
                "two frozen models on a benchmark this project has already scored; not fresh validation"
            ),
        },
        "prior_exposure_recomputed_not_quoted": {
            "note": (
                "the registration records the unpaired 200-resample figures that already existed "
                f"(0.4757 and 0.6393 pooled). They are exposure, not findings. Everything here is "
                f"recomputed at {re2g.DRAWS} draws"
            ),
            "pooled_dnase_plus_distance_now": out["pooled_descriptive"]["dnase_plus_distance_alone"]["auprc"],
            "pooled_full_model_now": out["pooled_descriptive"]["ours_auprc"],
        },
        "independence_and_exposure": dict(re2g.INDEPENDENCE),
        "falsifier": re2g.SECOND_REGISTRATION["falsifier"],
        "registration_unchanged": paired.registration_unchanged(),
        # This script's OWN import closure, not re2g_paired's: the lane's code list is shared
        # (paired.OWN_CODE names every file of it) but a counting path is a claim about
        # what could have entered THIS result, so it is computed from this entry.
        "code_cleanliness": mf.code_cleanliness("scripts/re2g_likeforlike.py", paired.OWN_CODE, paired.ROOT),
        "alphagenome_requests": 0,
    }
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "Gschwind et al. 2026, doi:10.1038/s41586-026-10781-4, Supplementary "
                "Tables 3 and 12",
                "version": "version of record",
            },
            {
                "accession": "EngreitzLab/CRISPR_comparison held-out benchmark and merge rule",
                "version": f"{re2g.BENCHMARK_TAG} ({re2g.BENCHMARK_COMMIT})",
            },
            {
                "accession": "ENCODE ENCODE-rE2G element gene links, "
                + ", ".join(f"{c}:{v['file']}" for c, v in sorted(re2g.COMPARATOR_FILES.items())),
                "version": "released 2024-06-06, software distal-regulation-encode_re2g",
            },
        ],
        "inputs": [
            mf.input_entry(
                crispri.KNOWLEDGE / crispri.HELDOUT,
                partition="held-out",
                url=crispri.BASE_URL + crispri.HELDOUT,
            ),
            mf.input_entry(
                crispri.KNOWLEDGE / crispri.TRAINING,
                partition="training",
                url=crispri.BASE_URL + crispri.TRAINING,
            ),
        ]
        + [
            mf.input_entry(
                paired.prediction_path(cell),
                partition=cell,
                url=paired.PORTAL.format(a=re2g.COMPARATOR_FILES[cell]["file"]),
                accession=re2g.COMPARATOR_FILES[cell]["file"],
            )
            for cell in sorted(re2g.COMPARATOR_FILES)
        ],
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "model": MODEL,
            "baseline": BASELINE,
            "draws_requested": re2g.DRAWS,
            "seed": re2g.SEED,
            "cluster": re2g.CLUSTER,
            "aggregate_function": "sum",
            "fill_value": 0,
            "gate_target": re2g.GATE_TARGET,
            "gate_tolerance": re2g.GATE_TOLERANCE,
        },
        "exclusions": [
            "no pair is excluded: the published pair set and labels are used whole",
            f"no number is reported for {', '.join(re2g.FEATURE_UNAVAILABLE)}: the deletion feature "
            "does not exist there",
        ],
        "partitions": {
            "primary_k562": "1918 held-out K562 pairs, the only registered population",
            "pooled_descriptive": f"all {len(heldout)} held-out pairs, descriptive",
            "gm12878_descriptive": "68 held-out GM12878 pairs, descriptive",
        },
    }
    path = save_result(RESULT, payload)
    print(f"{RESULT}: {path}")
    print("PRIMARY (K562, 1,918 pairs):", primary["delta_auprc"], primary["ci95"], "->", primary["reading"])
    print(json.dumps(payload["reading_must_travel_with"], indent=1))


if __name__ == "__main__":
    main()
