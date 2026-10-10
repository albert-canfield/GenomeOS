# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-envelope: the additive envelope's definition, its two identities, and its census on 7d's design.

    uv run python scripts/envelope_census.py [--no-pairs] [--no-save]

Order, so that no step can be mistaken for a reading of measured RNA:

  1. the ENVELOPE'S DEFINITION is checked against the engine, numerically and with no data read
     (`envelope.identity_checks`): the interval's upper end, with the cap it names, is grn's
     `sum_capped`, and its lower end is grn's `max`;
  2. given that identity, the form's CENSUS on 7d's population is lane-capsumrna's committed census of
     the same engine rule at the same conventions (data/results/capsumrna_census.json), read field by
     field, and 7d's DECISION 2 (`combine2.MIN_GENES` = 100) is applied to the count of genes whose fit
     is unique and whose prediction is not the mean's. NO RNA VALUE IS OPENED anywhere in this script;
  3. descriptive and deciding nothing, the PAIR SIDE: on lane-strengthform's already-read,
     already-committed combinations, each end of the envelope as a point form beside the interval,
     with the share of combinations inside the interval. The bootstrap used is that lane's own, at its
     seed, and it must first reproduce that lane's five committed intervals exactly.

No model request, no network, no money, no new biosample.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import combine2 as c2  # noqa: E402
from genomeos.attribution import envelope as ev  # noqa: E402
from genomeos.attribution import strengthform as sf  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "envelope_census"
OWN_CODE = (
    "genomeos/attribution/envelope.py",
    "scripts/envelope_census.py",
    "tests/test_envelope.py",
)


def pair_side() -> dict[str, Any]:
    """Each end of the envelope as a point form, on lane-strengthform's committed combinations."""
    table = sf.read_table(read_outcome=True)
    combos = table["combinations"]
    maxima = table["promoter_maxima"]
    errors = sf.errors_per_form(combos, maxima)
    points = ev.point_errors(combos, maxima)
    inside = int(sum(points.pop("inside")))
    errors.update(points)
    clusters = [c.cluster for c in combos]
    intervals = ev.paired_bootstrap(clusters, errors)
    reproduced = ev.check_against_committed(intervals)
    mse = {name: sf.mean_squared_error(values) for name, values in errors.items()}
    rows = {
        "upper": intervals["envelope_upper_point"],
        "interval": intervals["additive_envelope"],
        "lower": intervals["envelope_lower_point"],
    }
    sentence = ev.PAIR_SENTENCE.format(
        n=len(combos),
        inside=inside,
        inside_pct=round(100.0 * inside / len(combos), 1),
        upper_mse=f"{mse['envelope_upper_point']:.4f}",
        interval_mse=f"{mse['additive_envelope']:.4f}",
        reference_mse=f"{mse[ev.REFERENCE_FORM]:.4f}",
        upper_delta=f"{rows['upper']['delta_mse']:+.4f}",
        upper_lo=f"{rows['upper']['ci_low']:+.4f}",
        upper_hi=f"{rows['upper']['ci_high']:+.4f}",
        upper_word=sf.word(rows["upper"]),
        interval_delta=f"{rows['interval']['delta_mse']:+.4f}",
        interval_lo=f"{rows['interval']['ci_low']:+.4f}",
        interval_hi=f"{rows['interval']['ci_high']:+.4f}",
        interval_word=sf.word(rows["interval"]),
        lower_delta=f"{rows['lower']['delta_mse']:+.4f}",
        lower_lo=f"{rows['lower']['ci_low']:+.4f}",
        lower_hi=f"{rows['lower']['ci_high']:+.4f}",
        lower_word=sf.word(rows["lower"]),
    )
    return {
        "what": ev.PAIR_SIDE,
        "source": str(ev.STRENGTHFORM_RESULT),
        "commits": list(ev.STRENGTHFORM_COMMITS),
        "combinations": len(combos),
        "unordered_enhancer_pairs": len(set(clusters)),
        "promoters": len(Counter(c.promoter for c in combos)),
        "inside_the_interval": inside,
        "inside_the_interval_share": round(inside / len(combos), 4),
        "mse": mse,
        "difference_from_the_strength_dependent_form": intervals,
        "words": {name: sf.word(row) for name, row in intervals.items()},
        "committed_intervals_reproduced": reproduced,
        "sentence": sentence,
    }


def census(*, pairs: bool = True) -> dict[str, Any]:
    identity = ev.identity_checks()
    counts = ev.capsum_counts()
    status = counts["capped_sum_status"]
    unique = counts["unique_fit"]
    fitted = status["runnable"] + status["fit_not_unique"]
    differ = unique["differs_from_the_mean"]["genes"]
    decision = ev.decide(identity, differ)
    sentence = ev.CENSUS_SENTENCE.format(
        identity_worst=f"{identity['worst']:.3g}",
        multi_element=f"{counts['eligible_for_a_fit']:,}",
        fitted=f"{fitted:,}",
        cells=len(ev.CELLS),
        largest=f"{unique['max_abs_difference_from_mean_log2']:.3g}",
        unique=f"{unique['genes']:,}",
        differ=differ,
        min_genes=ev.MIN_GENES,
        two_roots=counts["two_roots"]["genes"],
    )
    out: dict[str, Any] = {
        "result": RESULT,
        "date": ev.DATE,
        "lane": ev.LANE,
        "roadmap_item": ev.ROADMAP_ITEM,
        "form": ev.FORM,
        "lower_end_form": ev.LOWER_FORM,
        "definition": ev.DEFINITION,
        "identity_rule": ev.IDENTITY,
        "identity": identity,
        "rna_side": ev.RNA_SIDE,
        "capsumrna_census": counts,
        "uncapped_root_fits": fitted,
        "genes_whose_unique_fit_differs_from_the_mean": differ,
        "decision_rule": ev.DECISION,
        "decision": decision,
        "size_check": {
            "genes_whose_unique_fit_differs_from_the_mean": differ,
            "min_genes": ev.MIN_GENES,
            "under_7d_decision_2": "NOT_TESTABLE" if differ < ev.MIN_GENES else "testable",
            "registration_written": False,
            "rna_opened": False,
        },
        "two_root_ambiguity_is_lane_rootpicks": {
            "genes": counts["two_roots"]["genes"],
            "note": (
                "those genes have a second, capped root beside the uncapped one; the definition does"
                " not choose between them. That choice is lane-rootpick's registered decision and this"
                " lane picks no root"
            ),
        },
        "sentence": sentence,
        "alphagenome_requests": 0,
        "money": "none",
        "new_biosample": "none: the thirteen reader cells unchanged, READER_TERMS untouched",
    }
    if pairs:
        out["pair_side_descriptive"] = pair_side()
    return out


def manifest(*, pairs: bool) -> dict[str, Any]:
    inputs = [mf.input_entry(ev.CAPSUM_CENSUS, partition=None)]
    if pairs:
        inputs.append(mf.input_entry(ev.STRENGTHFORM_RESULT, partition="held_out"))
        inputs.append(mf.input_entry(sf.TABLE, partition="held_out"))
    return {
        "sources": [
            {"accession": "data/results/capsumrna_census.json", "version": "f74dca5"},
            {"accession": sf.SOURCE["accession"], "version": "pinned by sha256"},
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38 for the RNA side; n/a for the reporter library",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "form": ev.FORM,
            "lower_end_form": ev.LOWER_FORM,
            "reference_form": ev.REFERENCE_FORM,
            "min_genes": ev.MIN_GENES,
            "identity_draws": ev.IDENTITY_DRAWS,
            "identity_seed": ev.IDENTITY_SEED,
            "identity_tolerance": ev.IDENTITY_TOLERANCE,
            "bootstrap_draws": sf.BOOTSTRAP_DRAWS,
            "bootstrap_seed": sf.BOOTSTRAP_SEED,
            "selection": c2.SELECTION,
        },
        "exclusions": [
            "the RNA side reads no RNA value: the counts are lane-capsumrna's committed census of the"
            " same engine rule at the same conventions",
            "the pair side is descriptive and registers nothing; its outcomes were read and committed"
            " by lane-strengthform (fc328c5)",
        ],
        "partitions": {
            "held_out": "lane-strengthform's combinations, held out from every other lane",
            "training": "none: nothing is fitted here",
        },
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-pairs", action="store_true", help="skip the descriptive pair side")
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    out = census(pairs=not args.no_pairs)
    out["seconds"] = round(time.time() - t0, 1)
    print(json.dumps(out["identity"], indent=1))
    print(out["sentence"])
    print("decision:", out["decision"], "size check:", json.dumps(out["size_check"]))
    if "pair_side_descriptive" in out:
        pair = out["pair_side_descriptive"]
        print("committed intervals reproduced:", pair["committed_intervals_reproduced"]["problems"])
        print(pair["sentence"])
    if not args.no_save:
        print("saved", save_result(RESULT, out, manifest=manifest(pairs=not args.no_pairs)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
