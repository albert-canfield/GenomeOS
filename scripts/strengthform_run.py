# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-strengthform, the registered run: the strength-dependent form on held-out enhancer pairs.

    uv run python scripts/strengthform_run.py

Runs only from the committed registration (data/results/strengthform_registration.json), and refuses
if any registered field differs from `genomeos.attribution.strengthform`, if the table's sha256
differs from the pinned one, or if 7d's declared THETA has moved.

Order, so that a defect in the harness cannot be mistaken for a reading:
  1. the registration is checked field by field against the module, and the table against its digest;
  2. the HARNESS is checked before the outcome column is opened -- the form reproduces a single it is
     given, the census is recounted with no outcome read and must equal the committed census, and a
     planted pair of two equal maximal singles must read as each form's algebra says;
  3. only then are the enhancer-enhancer activities opened, scored under every form on the same
     combinations, and the reading rendered from the intervals alone.

No model request, no network, no money.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import combine2 as c2  # noqa: E402
from genomeos.attribution import strengthform as sf  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from strengthform_register import OWN_CODE, census, inputs  # noqa: E402

RESULT = "strengthform_result"
REGISTRATION = RESULTS_DIR / "strengthform_registration.json"
#: a partner that contributes nothing: the form must then return the other single's own boost
NULL_PARTNER = 1e-12


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check_registration() -> dict[str, Any]:
    """Refuse the run unless the committed registration still says what the module says."""
    if not REGISTRATION.exists():
        raise SystemExit(f"refused: {REGISTRATION} does not exist; register before running")
    reg = json.loads(REGISTRATION.read_text())
    refusals = []
    for key, value in sf.registered_fields().items():
        if reg.get(key) != value:
            refusals.append(f"registered {key} differs from the module")
    if reg.get("theta") != c2.THETA:
        refusals.append("7d's declared THETA has moved since the registration")
    got = _sha(sf.TABLE) if sf.TABLE.exists() else "absent"
    if got != sf.TABLE_SHA256:
        refusals.append(f"the table's sha256 is {got}, the pinned digest is {sf.TABLE_SHA256}")
    if refusals:
        raise SystemExit(f"refused: {refusals}")
    return reg


def check_harness(reg: dict[str, Any]) -> dict[str, Any]:
    """Every check that can be made before the outcome column is opened."""
    recounted = census()
    if recounted != reg["census"]:
        raise SystemExit("refused (HARNESS_DEFECT): the census does not reproduce the committed one")
    singles = []
    for bmax in (0.25, 1.0, 3.0):
        for b in (0.05, 0.5, bmax):
            if b > bmax:
                continue
            singles.append(
                {
                    "bmax": bmax,
                    "b": b,
                    "term": sf.terms(b, b, bmax)[0],
                    "with_a_null_partner": sf.predicted_boost(sf.FORM, b, NULL_PARTNER, bmax),
                }
            )
    worst_single = max(abs(row["with_a_null_partner"] - row["b"]) for row in singles)
    if worst_single > 1e-9:
        raise SystemExit(f"refused (HARNESS_DEFECT): a single is not reproduced, {worst_single}")
    planted: dict[str, Any] = {}
    b = 1.0
    for form in (sf.FORM, *sf.COMPARATORS):
        prediction = sf.predicted_boost(form, b, b, b)
        planted[form] = list(prediction) if isinstance(prediction, tuple) else prediction
    additive = math.log2(2.0 * 2.0**b - 1.0)
    if abs(planted["multiplicative"] - 2 * b) > 1e-12 or abs(planted["max"] - b) > 1e-12:
        raise SystemExit("refused (HARNESS_DEFECT): a planted pair's fixed-rule prediction is wrong")
    if not planted["max"] < planted[sf.FORM] < additive:
        raise SystemExit(
            "refused (HARNESS_DEFECT): the strength-dependent form must lie strictly between the"
            " stronger single and the additive expectation on two equal maximal singles"
        )
    return {
        "census_reproduced": True,
        "singles_reproduced": singles,
        "worst_single_reproduction_error_log2": worst_single,
        "planted_two_equal_maximal_singles": planted,
        "planted_additive_expectation": additive,
    }


def main() -> int:
    started = time.time()
    reg = check_registration()
    harness = check_harness(reg)
    table = sf.read_table(read_outcome=True)
    combos = table["combinations"]
    maxima = table["promoter_maxima"]
    errors = sf.errors_per_form(combos, maxima)
    clusters = [c.cluster for c in combos]
    intervals = sf.cluster_bootstrap(clusters, errors)
    own_mse = sf.mean_squared_error(errors[sf.FORM])
    per_promoter = Counter(c.promoter for c in combos)
    reading = sf.reading(
        intervals,
        own_mse=own_mse,
        n=len(combos),
        clusters=len(set(clusters)),
        promoters=len(per_promoter),
    )
    above_additive = 0
    below_stronger = 0
    for c in combos:
        observed = c.observed_boost
        lower, upper = sf.predicted_boost("additive_envelope", c.b1, c.b2, maxima[c.promoter])
        above_additive += observed > upper
        below_stronger += observed < lower
    entries = inputs()
    out: dict[str, Any] = {
        "lane": sf.LANE,
        "roadmap_item": sf.ROADMAP_ITEM,
        "registration": str(REGISTRATION),
        "reading": reading["sentence"],
        "words": reading["words"],
        "falsified": reading["falsified"],
        "mse": {name: sf.mean_squared_error(values) for name, values in errors.items()},
        "intervals": intervals,
        "population": {
            "informative_combinations": len(combos),
            "unordered_enhancer_pairs": len(set(clusters)),
            "combinations_per_promoter": dict(sorted(per_promoter.items())),
            "counts": table["counts"],
        },
        "secondary_deciding_nothing": {
            "what": sf.SECONDARY,
            "above_the_additive_expectation": above_additive,
            "below_the_stronger_single": below_stronger,
            "of": len(combos),
        },
        "harness": harness,
        "alphagenome_requests": 0,
        "money": "none",
        "seconds": round(time.time() - started, 1),
    }
    manifest = {
        "sources": [
            {"accession": sf.SOURCE["accession"], "version": "pinned by sha256"},
            {"accession": sf.SOURCE["publication"], "version": "2024"},
        ],
        "inputs": [*entries, mf.input_entry(REGISTRATION, partition=None)],
        "input_count": len(entries) + 1,
        "assembly": "n/a: a barcoded episomal reporter library, not genomic coordinates",
        "coordinates": {"base": "n/a", "interval": "n/a"},
        "parameters": {
            "form": sf.FORM,
            "theta": sf.THETA,
            "half_saturation": sf.HALF_SATURATION,
            "comparators": list(sf.COMPARATORS),
            "bootstrap_draws": sf.BOOTSTRAP_DRAWS,
            "bootstrap_seed": sf.BOOTSTRAP_SEED,
        },
        "exclusions": list(sf.INFORMATIVE),
        "partitions": "held out: this dataset has never been read by any earlier lane of the project",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }
    print(reading["sentence"])
    print("saved", save_result(RESULT, out, manifest=manifest))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
