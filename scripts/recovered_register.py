# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-recovered step 2: write the registration, BEFORE any P, D or Delta MSE of the population is read.

    uv run python scripts/recovered_register.py [--no-save]

This writes `data/results/recovered_registration.json` from `genomeos/attribution/recovered.py`'s BOUND
constants and nothing else: the question, the one form, the fixed scale and its source, the rows, the
endpoint side, the metric, the interval and the reason for its width, the bootstrap, the control and the
one way it differs from the endpoint, the harness, THE EXPOSURE, the supervisor's four-way table with the
scale interval it is evaluated on, the decision, the falsifier, the stop rule, the limitations and the
descriptives. It also pins the population committed at step 1 by its two sha256, so the run cannot read a
different set of genes than the one established before this registration existed.

NO RNA VALUE IS OPENED HERE EITHER. The only thing read from a result file is the population record's
symbol lists and its sha256, plus M2h's committed scale and scale interval, which this registration must
reproduce exactly.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import recovered as rc  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "recovered_registration"
POPULATION = Path("data/results/recovered_population.json")
MULTRNA = Path("data/results/multrna_result.json")
MAGCAL = Path("data/results/magcal_result.json")
FRESHRNA = Path("data/results/freshrna_result.json")
SCALEFIT_CENSUS = Path("data/results/scalefit_census.json")
OWN_CODE = (
    "genomeos/attribution/recovered.py",
    "scripts/recovered_population.py",
    "scripts/recovered_register.py",
    "scripts/recovered_run.py",
    "tests/test_recovered.py",
)


def population_pin() -> dict:
    """The population established at step 1, pinned by its two sha256 and re-derived from its lists."""
    d = json.loads(POPULATION.read_text())
    s = d["symbols"]
    if d["outcome"] != "POPULATION_ESTABLISHED" or not d["population_harness"]["agrees"]:
        raise SystemExit("the population record is not established; nothing is registered")
    for name, symbols, want in (
        ("population_sha256", s["all"], d["population_sha256"]),
        ("held_out_sha256", s["held_out"], d["held_out_sha256"]),
    ):
        got = rc.sha256_of_symbols(symbols)
        if got != want:
            raise SystemExit(f"{name}: the record's own symbol list does not reproduce it")
    return {
        "source": str(POPULATION),
        "genes": len(s["all"]),
        "held_out_genes": len(s["held_out"]),
        "training_genes": len(s["training"]),
        "population_sha256": d["population_sha256"],
        "held_out_sha256": d["held_out_sha256"],
        "re_derived_counts": d["re_derived_counts"],
        "overlap_with_m2h": d["overlap_with_m2h"],
    }


def scale_pin() -> dict:
    """M2h's committed scale and training interval, which this registration reproduces exactly."""
    d = json.loads(MAGCAL.read_text())
    if d["scales"]["multiplicative"] != rc.SCALE or d["scale_cis"]["multiplicative"] != list(rc.SCALE_CI):
        raise SystemExit("the registered scale or its interval differs from M2h's committed one")
    return {
        "source": f"{MAGCAL} (M2h, commit ac5580b)",
        "scale": rc.SCALE,
        "scale_rounded": rc.SCALE_ROUNDED,
        "scale_ci": list(rc.SCALE_CI),
        "scale_interval_excludes_zero": not (rc.SCALE_CI[0] <= 0.0 <= rc.SCALE_CI[1]),
        "is_refitted_here": False,
    }


def context_verbatim() -> dict:
    """M2h's and M2i's registered readings, quoted word for word, never re-decided."""
    m2h = json.loads(MAGCAL.read_text())
    m2i = json.loads(FRESHRNA.read_text())
    return {
        "m2h_multiplicative_reading": m2h["per_form"]["multiplicative"]["reading"],
        "m2h_variance_share": "about 2.0% of the held-out mean(D^2) of 2.1524",
        "m2i_verdict_sentence": m2i["verdict_sentence"],
        "m2i_variance_fraction": m2i["variance_fraction"],
        "scalefit_verdict": json.loads(SCALEFIT_CENSUS.read_text())["verdict"],
        "note": "quoted, not re-decided; summaries keep a lane's registered wording word for word",
    }


def registration() -> dict:
    return {
        "lane": rc.LANE,
        "date": time.strftime("%Y-%m-%d"),
        "registered": rc.REGISTERED,
        "step": "2 of 3: the registration, written BEFORE any P, D or Delta MSE of this population",
        "what": (
            "the RNA test lane-scalefit registered as separate: on the genes M2h's fixed 0.2146 makes"
            " fittable, is the scaled run closer to measured RNA than the no-rule baseline"
        ),
        **rc.registered_fields(),
        "population_pin": population_pin(),
        "scale_pin": scale_pin(),
        "context_verbatim": context_verbatim(),
        "bound_fields": sorted(rc.BOUND),
        "rna_read": "none: this step registers and reads no P, no D and no Delta MSE",
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }


def manifest() -> dict:
    inputs = [
        mf.input_entry(POPULATION, partition=None),
        mf.input_entry(MULTRNA, partition=None),
        mf.input_entry(MAGCAL, partition=None),
        mf.input_entry(FRESHRNA, partition=None),
        mf.input_entry(SCALEFIT_CENSUS, partition=None),
    ]
    return {
        "sources": [
            {
                "accession": "data/results/recovered_population.json (step 1)",
                "version": "pinned by sha256",
            },
            {"accession": "data/results/multrna_result.json (7h, the P and D rows)", "version": "pinned"},
            {
                "accession": "data/results/magcal_result.json (M2h, the scale and the split)",
                "version": "ac5580b",
            },
            {"accession": "data/results/freshrna_result.json (M2i, context)", "version": "pinned"},
            {"accession": "data/results/scalefit_census.json (lane-scalefit)", "version": "7f65544"},
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "form": rc.FORM,
            "population_form": rc.POPULATION_FORM,
            "scale": rc.SCALE,
            "scale_is_refitted": False,
            "interval": list(rc.INTERVAL),
            "bootstrap_resamples": rc.BOOTSTRAP["resamples"],
            "bootstrap_seed": rc.BOOTSTRAP["seed"],
            "permutations": rc.PERMUTATIONS,
            "min_genes": rc.MIN_GENES,
        },
        "exclusions": ["none in this step: a registration excludes nothing, it declares"],
        "partitions": "n/a: a registration; no outcome read",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args()
    out = registration()
    print("registered fields", len(rc.BOUND), flush=True)
    print("population", out["population_pin"]["genes"], "held out", out["population_pin"]["held_out_genes"])
    print("EXPOSURE:", out["exposure"][:400], flush=True)
    if not args.no_save:
        print("saved", save_result(RESULT, out, manifest=manifest()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
