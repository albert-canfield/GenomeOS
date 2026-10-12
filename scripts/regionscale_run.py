# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-regionscale, step 3: THE ONE RUN. Does a per-region magnitude beat the single global scale?

    uv run python scripts/regionscale_run.py

It reads the committed registration, refuses to start unless that registration binds every BOUND
field of `genomeos/attribution/regionscale.py` word for word and every pinned input's sha256 is
unmoved, fits one scale per admitted chromosome on THAT chromosome's TRAINING genes alone, scores
M2h's own held-out genes under the per-region arm and under M2h's single global scale, and reports
DELTA with the registered interval, the registered permutation null and the registered sentence.

ONE RUN (regionscale.STOP_RULE). It refuses to overwrite its own result. No second seed, no second
split, no second interval width, no second regional unit and no variant is run here.

No model request, no network, no money, 0 AlphaGenome requests.
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
from genomeos.attribution import magcal as mc  # noqa: E402
from genomeos.attribution import regionscale as rs  # noqa: E402
from genomeos.results import save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from regionscale_harness import (  # noqa: E402
    OWN_CODE,
    READERS,
    SOURCE_FILES,
    load_rows,
)

RESULT = "regionscale_result"
OUT = "data/results/regionscale_result.json"
REGISTRATION = "data/results/regionscale_registration.json"
HARNESS_RESULT = "data/results/regionscale_harness.json"

#: DISCLOSED IN THE RESULT ITSELF, because a reader cannot check a commit message.
#: On 2026-10-12 at 19:38:34, before the registration was committed, lane-regionscale2 ran this script
#: ONCE as a smoke test, expecting it to refuse because no registration existed -- but the registration
#: had just been rewritten by the committed writer, so the three refusals passed and the script computed
#: DELTA, its bootstrap and its permutation null. NOTHING WAS WRITTEN INTO data/results: `save_result`
#: refused on the entry script, which git did not yet track, and quarantined the record at
#: data/quarantine/results/regionscale_result.json (sha256 37fd95a1..., 49,121 bytes, git-ignored). The
#: reading is printed only AFTER save_result returns, so it never printed; no value of it was read, and
#: the file was moved OUT of the repository into this session's scratchpad, unopened, so it cannot be
#: taken for the real result. WHAT SETTLES THE REGISTRATION'S HONESTY IS GIT AND NOT AN ASSURANCE, AND
#: NOT A FILE MTIME EITHER: every registered field -- the question, the held-out design, the admission
#: floors, the comparison, the statistic, the interval and its width, the alpha and its correction, the
#: null, the p-value form, the decision rule, the falsifier, BOTH pre-written outcome sentences, the stop
#: rule, the limitations and the descriptives -- lives in genomeos/attribution/regionscale.py, which is
#: UNCHANGED since 027d414 at 16:52:54, two hours and forty-five minutes before the trial, and the
#: registration is a verbatim copy of those 45 fields with a committed test asserting the copy word for
#: word. So no registered sentence could have been tuned to a number that did not exist when it was
#: fixed; the registration FILE was rewritten after the trial, and what it contains was not.
#: WHAT IT COST: the claim "no held-out value of this design had been computed before the registration
#: was committed" is NO LONGER TRUE and is made nowhere. The claim that is true, and is the one made, is
#: that no such value was READ, and that every registered sentence predates the computation in git.
TRIAL_RUN_DISCLOSURE = (
    "A TRIAL RUN HAPPENED BEFORE THE REGISTRATION WAS COMMITTED, AND IT IS DISCLOSED HERE RATHER THAN"
    " ONLY IN A COMMIT MESSAGE. lane-regionscale2 ran scripts/regionscale_run.py once at 19:38:34 on"
    " 2026-10-12 as a smoke test, expecting the registration-absent refusal; the registration had just"
    " been rewritten by the committed writer, so the refusals passed and the script computed DELTA, its"
    " paired bootstrap and its permutation null. NOTHING REACHED data/results: save_result refused"
    " because the entry script was not yet tracked by git and quarantined the record at"
    " data/quarantine/results/regionscale_result.json, sha256"
    " 37fd95a18592207ea70eb2e1d19ea4faaa67959ae6dd7c2f591e19da0fa2f1ea, 49,121 bytes. The reading prints"
    " only after save_result returns, so nothing printed; the file was moved out of the repository"
    " unopened and no value in it was read. WHAT SETTLES THE REGISTRATION'S HONESTY IS GIT, NOT AN"
    " ASSURANCE AND NOT A FILE MTIME: every registered field -- the question, the held-out design, the"
    " admission floors, the comparison, the statistic, the interval and its width, the alpha and its"
    " correction, the null, the p-value form, the decision rule, the falsifier, BOTH pre-written outcome"
    " sentences, the stop rule, the limitations and the descriptives -- lives in"
    " genomeos/attribution/regionscale.py, which is UNCHANGED since 027d414 at 16:52:54, two hours and"
    " forty-five minutes BEFORE the trial, and the registration is a verbatim copy of those 45 fields"
    " with a committed test asserting the copy word for word. No registered sentence could therefore"
    " have been tuned to a number that did not exist when it was fixed: the registration FILE was"
    " rewritten after the trial and what it CONTAINS was not. WHAT IS NO LONGER CLAIMABLE: that no"
    " held-out value of this design had been COMPUTED before the registration was committed. What is"
    " claimed, and is what the stop rule protects, is that none was READ, and that every registered"
    " sentence predates the computation in git"
)


# ---- the three refusals, before anything is fitted ------------------------------------------------
def refuse_unless_registered() -> dict[str, Any]:
    """STOP_RULE's three conditions, checked before a single scale is fitted.

    The registration must EXIST (it is committed with its tests before any held-out value of this
    design is read), it must bind every BOUND field WORD FOR WORD (a registration that drifted from
    the module is a registration of a different experiment), and every input it pinned must still
    have the sha256 it pinned (a moved input makes the registered numbers unreproducible). And the
    result is written ONCE: a second run of this script refuses rather than overwriting.
    """
    if Path(OUT).exists():
        raise SystemExit(f"refused: {OUT} already exists; ONE RUN (STOP_RULE)")
    if not Path(REGISTRATION).exists():
        raise SystemExit(f"refused: {REGISTRATION} does not exist; register and commit it first")
    reg = json.loads(Path(REGISTRATION).read_text())
    drift = [f for f, v in rs.registered_fields().items() if reg.get(f) != v]
    if drift:
        raise SystemExit(f"refused: the registration does not bind these fields word for word: {drift}")
    moved = []
    for path, want in reg["pinned_sha256"].items():
        got = mf.input_entry(path)["sha256"]
        if got != want:
            moved.append(f"{path}: pinned {want[:12]}, now {got[:12]}")
    if moved:
        raise SystemExit("refused: a pinned input moved since the registration: " + "; ".join(moved))
    h = json.loads(Path(HARNESS_RESULT).read_text())
    if h["outcome"] != "HARNESS_REPRODUCES":
        raise SystemExit("refused: the committed harness did not reproduce M2h's numbers")
    return reg


# ---- the two arms --------------------------------------------------------------------------------
def held_out_variance(held: list[dict[str, Any]]) -> tuple[float, int]:
    """The denominator of every "about X% of the variance": the no-rule baseline's held-out MSE.

    The harness computed this same quantity and confirmed M2h's committed "about 2% of the variance"
    against it; it is recomputed here because every percentage this run prints is taken against it.
    """
    dd = 0.0
    n = 0
    for r in held:
        for d in r["D"]:
            dd += d * d
            n += 1
    return (dd / n if n else float("nan")), n


def main() -> int:  # noqa: C901
    t0 = time.time()
    reg = refuse_unless_registered()

    rows = load_rows()
    train, held = mc.split_rows(rows)
    tr_dps, tr_pps, _tr_ns = rs.gene_sums(train)
    ho_dps, ho_pps, ho_ns = rs.gene_sums(held)
    tr_lab = [r["chrom"] for r in train]
    ho_lab = [r["chrom"] for r in held]

    # the global arm: ONE scale fitted on ALL training genes, which must be M2h's own committed number
    a_global = rs.scale_from_sums(sum(tr_dps), sum(tr_pps))
    if a_global is None:
        raise SystemExit("refused: the global fit is undefined on the training genes")
    committed = rs.GLOBAL_SCALE_COMMITTED
    global_check = {
        "fitted_here": a_global,
        "m2h_committed": committed,
        "abs_difference": abs(a_global - committed),
        "m2h_tolerance": rs.M2H_TOLERANCE,
        "agrees": abs(a_global - committed) <= rs.M2H_TOLERANCE,
        "what_it_is": (
            "a CHECK and not a value: the same fit on the same training genes must return M2h's"
            " committed scale, and it is reported so a reader can see that it did"
        ),
    }
    if not global_check["agrees"]:
        raise SystemExit(f"refused: the global scale fitted here is not M2h's committed one: {global_check}")

    counts = rs.split_counts(rows)
    admitted, fallback = rs.admitted_regions(counts)
    testable = rs.is_testable(admitted)
    design = {
        "regional_unit": "chromosome",
        "min_genes_each_half": rs.MIN_GENES,
        "min_regions": rs.MIN_REGIONS,
        "admitted": admitted,
        "admitted_count": len(admitted),
        "fallback_to_global_both_arms": fallback,
        "fallback_count": len(fallback),
        "genes_per_chromosome_per_half": {c: counts[c] for c in sorted(counts)},
        "held_out_genes_in_admitted": sum(counts[c]["held"] for c in admitted),
        "held_out_genes_in_fallback": sum(counts[c]["held"] for c in fallback),
        "testable": testable,
    }

    var, held_pairs = held_out_variance(held)
    nums: dict[str, Any] = {
        "held_genes": len(held),
        "held_pairs": held_pairs,
        "n_admitted": len(admitted),
        "n_fallback": len(fallback),
        "a_global": a_global,
        "min_regions": rs.MIN_REGIONS,
        "fw": rs.FORM_INTERVAL_WORDS,
        "nw": rs.NOMINAL_INTERVAL_WORDS,
        "n_perm": rs.PERMUTATIONS["count"],
    }

    # DECISION (2): not testable -> no DELTA is computed at all, and that is the finding
    if not testable:
        kind = rs.decide_region(0.0, 0.0, 0.0, 1.0, len(admitted))
        out = assemble(reg, kind, nums, design, global_check, var, None, t0)
        return publish(out)

    per_region, degenerate = rs.region_scales(tr_dps, tr_pps, tr_lab, admitted)
    # DECISION (3): an undefined fit -> no DELTA is reported
    if degenerate:
        nums["what"] = (
            "The training fit is undefined for "
            + ", ".join(degenerate)
            + ", which the registration admitted and which this run therefore reports instead of a"
            " comparison."
        )
        kind = rs.decide_region(0.0, 0.0, 0.0, 1.0, len(admitted), degenerate)
        out = assemble(reg, kind, nums, design, global_check, var, None, t0)
        return publish(out)

    scales_region = rs.scales_for(ho_lab, per_region, a_global)
    scales_global = [a_global] * len(ho_lab)
    dm_region = rs.delta_mse_rowwise(ho_dps, ho_pps, ho_ns, scales_region)
    dm_global = rs.delta_mse_rowwise(ho_dps, ho_pps, ho_ns, scales_global)
    delta = dm_region - dm_global

    boot = rs.bootstrap_paired(
        ho_dps,
        ho_pps,
        ho_ns,
        scales_region,
        scales_global,
        rs.BOOTSTRAP["resamples"],
        rs.BOOTSTRAP["seed"],
    )
    null = rs.permuted_delta(
        tr_dps,
        tr_pps,
        tr_lab,
        ho_dps,
        ho_pps,
        ho_ns,
        ho_lab,
        admitted,
        a_global,
        rs.PERMUTATIONS["count"],
        rs.PERMUTATIONS["seed"],
    )
    p_perm = rs.permutation_p(delta, null)
    lo, hi = boot["delta_ci"]
    nom_lo, nom_hi = boot["delta_ci_nominal"]
    kind = rs.decide_region(lo, hi, nom_lo, p_perm, len(admitted))

    nums.update(
        {
            "delta": delta,
            "delta_lo": lo,
            "delta_hi": hi,
            "nom_lo": nom_lo,
            "nom_hi": nom_hi,
            "p_perm": p_perm,
            "dm_region": dm_region,
            "dm_global": dm_global,
            "pct_region": 100.0 * dm_region / var,
            "pct_global": 100.0 * dm_global / var,
            "pct_delta": 100.0 * delta / var,
            "pct_delta_abs": abs(100.0 * delta / var),
        }
    )
    if kind == "MIXED":
        nums["why"] = rs.why_mixed(lo, p_perm)

    # ---- the descriptives, which decide nothing (DESCRIPTIVES) ----------------------------------
    per_chrom = []
    for c in admitted:
        idx = [i for i, lab in enumerate(ho_lab) if lab == c]
        d_ho = [ho_dps[i] for i in idx]
        p_ho = [ho_pps[i] for i in idx]
        n_ho = [ho_ns[i] for i in idx]
        a_c = per_region[c]
        own = rs.delta_mse_rowwise(d_ho, p_ho, n_ho, [a_c] * len(idx))
        glob = rs.delta_mse_rowwise(d_ho, p_ho, n_ho, [a_global] * len(idx))
        per_chrom.append(
            {
                "chrom": c,
                "training_genes": counts[c]["train"],
                "held_out_genes": counts[c]["held"],
                "scale": a_c,
                "distance_from_global": a_c - a_global,
                "held_out_delta_mse_own_scale": own,
                "held_out_delta_mse_global_scale": glob,
                "difference": own - glob,
            }
        )

    descriptives = {
        "per_admitted_chromosome": per_chrom,
        "per_admitted_chromosome_decides_nothing": rs.DESCRIPTIVES[0],
        "null_distribution_of_delta": {
            "median": sorted(null)[len(null) // 2],
            "percentile_range_98_33": list(rs.percentile_range(null)),
            "quartiles": rs.quantiles(null),
            "count": len(null),
            "seed": rs.PERMUTATIONS["seed"],
        },
        "arms": {
            "region_arm_delta_mse": dm_region,
            "region_arm_ci": boot["region_arm_ci"],
            "region_arm_percent_of_variance": nums["pct_region"],
            "global_arm_delta_mse": dm_global,
            "global_arm_ci": boot["global_arm_ci"],
            "global_arm_percent_of_variance": nums["pct_global"],
            "m2h_committed_held_out_delta_mse": rs.M2H_COMMITTED["held_out_delta_mse"],
            "m2h_committed_wording": "about 2% of the variance",
        },
        "global_scale_check": global_check,
        "p_against_the_stricter_divisor": {
            "p_perm": p_perm,
            "threshold": 0.0167,
            "clears": p_perm <= 0.0167,
            "what_it_is": (
                "a DESCRIPTIVE and never the verdict: 0.0167 is the 0.05/3 divisor M2h's 98.33% width"
                " came from and the threshold at which hetblock2's secondary statistic did NOT clear."
                " THE VERDICT USES 0.05 AND ONLY 0.05 (ALPHA)"
            ),
        },
        "held_out_variance": var,
        "held_out_variance_is": (
            "the no-rule baseline's own held-out MSE, the mean of D^2 over the held-out gene-cell"
            " pairs, and the denominator every 'about X% of the variance' here is taken against"
        ),
    }
    measured = {
        "delta": delta,
        "delta_ci": boot["delta_ci"],
        "delta_ci_nominal": boot["delta_ci_nominal"],
        "delta_draws_quantiles": boot["delta_draws_quantiles"],
        "resamples": boot["resamples"],
        "bootstrap_seed": boot["seed"],
        "p_perm": p_perm,
        "permutations": rs.PERMUTATIONS["count"],
        "permutation_resolution": 1.0 / (1.0 + rs.PERMUTATIONS["count"]),
        "region_arm_delta_mse": dm_region,
        "global_arm_delta_mse": dm_global,
        "per_region_scales": {c: per_region[c] for c in admitted},
    }
    out = assemble(reg, kind, nums, design, global_check, var, (measured, descriptives), t0)
    return publish(out)


def assemble(
    reg: dict[str, Any],
    kind: str,
    nums: dict[str, Any],
    design: dict[str, Any],
    global_check: dict[str, Any],
    var: float,
    computed: tuple[dict[str, Any], dict[str, Any]] | None,
    t0: float,
) -> dict[str, Any]:
    """The result record: the registered sentence with its numbers filled, and nothing else added."""
    out: dict[str, Any] = {
        "lane": "lane-regionscale2",
        "step": "3 of 4: the one run, after the registration and its tests were committed",
        "registration": REGISTRATION,
        "registered": reg["registered"],
        "date": date.today().isoformat(),
        "question": rs.QUESTION,
        "outcome": kind,
        "reading": rs.full_reading(kind, **nums),
        "which_governs": rs.WHICH_GOVERNS,
        "falsifier": rs.FALSIFIER,
        "decision": rs.DECISION,
        "prior_exposure": list(rs.PRIOR_EXPOSURE),
        "prior_exposure_clause": rs.PRIOR_EXPOSURE_CLAUSE,
        "held_out_design": rs.HELD_OUT_DESIGN,
        "design": design,
        "held_out_variance": var,
        "global_scale_check": global_check,
        "limitations": list(rs.LIMITATIONS),
        "no_further_model_is_proposed": (
            "WHATEVER THE OUTCOME, NO FURTHER MODEL IS PROPOSED BY THIS LANE: no per-region model, no"
            " shrunken estimator, no second regional unit and no next step. The per-chromosome table is"
            " a DESCRIPTIVE and is not a row of verdicts (STOP_RULE, DECISION step 7)"
        ),
        "trial_run_before_the_registration_was_committed": TRIAL_RUN_DISCLOSURE,
        "readers": READERS,
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }
    if computed is not None:
        measured, descriptives = computed
        out["measured"] = measured
        out["descriptives"] = descriptives
    else:
        out["nothing_was_fitted_or_scored"] = (
            "the registration forbids reading a held-out value under a design it has declared"
            " dishonest or a fit it has declared undefined, so NO DELTA, no per-region scale, no"
            " bootstrap and no permutation exists in this result"
        )
    out["seconds"] = round(time.time() - t0, 1)
    return out


def manifest(entries: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "sources": [
            {"accession": "data/results/multrna_result.json (7h, lane-multrna)", "version": "pinned"},
            {"accession": "data/results/magcal_result.json (M2h, lane-magcal)", "version": "pinned"},
            {
                "accession": "data/results/scalehet_result.json (lane-scalehet)",
                "version": "pinned; its rejection of a single global scale is why this lane exists",
            },
            {
                "accession": "data/results/hetblock_result.json (lane-hetblock2)",
                "version": "pinned; its locus-block re-test of that rejection, carried word for word",
            },
            {
                "accession": "data/results/regionscale_registration.json (this lane, step 2)",
                "version": "pinned; the run refuses to start unless it binds every BOUND field",
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
            "form": rs.FORM,
            "split_rule": "magcal.SPLIT, sha256 of the gene symbol, carried UNCHANGED",
            "regional_unit": "chromosome, inherited from scalehet and hetblock2",
            "min_genes_each_half": rs.MIN_GENES,
            "min_regions": rs.MIN_REGIONS,
            "form_interval": list(rs.FORM_INTERVAL),
            "nominal_interval": list(rs.NOMINAL_INTERVAL),
            "bootstrap": {
                "resamples": rs.BOOTSTRAP["resamples"],
                "seed": rs.BOOTSTRAP["seed"],
            },
            "permutations": rs.PERMUTATIONS,
            "control_alpha": rs.CONTROL_ALPHA,
            "tolerance": rs.TOLERANCE,
            "m2h_tolerance": rs.M2H_TOLERANCE,
        },
        "exclusions": [
            "none: the population is 7h's, unchanged. A chromosome that does not clear the ADMISSION"
            " floors is NOT excluded -- its held-out genes stay in the scored set and are scored by the"
            " global scale under BOTH arms"
        ],
        "partitions": {
            "this_lanes_own": (
                "M2h's sha256-of-symbol partition, taken WITHIN each chromosome: a region's scale is"
                " fitted on its training genes alone and scored on its held-out genes. Both arms are"
                " scored on the SAME held-out rows and neither arm saw them"
            ),
        },
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def publish(out: dict[str, Any]) -> int:
    entries = [
        mf.input_entry(
            p,
            partition=(
                "pinned by sha256; the source of the committed P and D"
                if p.endswith("multrna_result.json")
                else "pinned by sha256; committed numbers and readings, quoted and not re-decided"
            ),
        )
        for p in (*SOURCE_FILES, HARNESS_RESULT, REGISTRATION)
    ]
    print("saved", save_result(RESULT, out, manifest=manifest(entries)))
    print(out["outcome"], "|", out["reading"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
