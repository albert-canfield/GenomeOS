"""The AstroREG request list, enumerated and priced, with NOTHING SENT.

This is what the executor reads before any money is spent, so its one job is to be exactly the list
the real run would send and not a second enumeration that could drift from it. The list comes from
`astrorun.plan_requests`, which is the only enumeration there is: a sender, if one is ever
authorised, iterates THAT function, so the thing reviewed here and the thing sent are the same
object. `astrorun.ONE_ENUMERATION_ONLY` says why that matters.

IT TOUCHES NO NETWORK. There is no AlphaGenome import in this file, no HTTP client, and no code path
that could send a request: what it does is read the frozen registration, read the screen's pairs, read
the swept element registry, and print. `--send` is not a flag it has. That is structural and not a
promise -- `tests/test_astrorun.py` reads this file's own source and fails if an AlphaGenome or HTTP
import ever appears in it.

WHAT IT CHECKS, and refuses on:

  (1) the registration is byte-identical to the authorised copy (sha256), so the number the list is
      checked against is the number Albert approved;
  (2) the enumerated count equals the registration's `requests.total_requests_needed` AND its
      `coverage.registry_elements_needed_total`. A list that is one request longer than the authorised
      number is not a rounding difference, it is an unauthorised request, so a mismatch is a refusal
      and not a warning;
  (3) the activity gate's committed verdict. The gate (scripts/astroreg_activity_gate.py) found that
      the frozen model's activity term cannot be supplied for an astrocyte pair in the units its
      weights were fitted on, and the registration's own feasibility rule answers that case with a
      no-go. So this dry run reports `requests_permitted: 0` and says the list may not be bought,
      quoting the gate rather than deciding again. If the gate result is absent it refuses, because an
      unreviewed list is not a reviewed one.

The list is written out in full so it can be diffed between runs: the enumeration is sorted and
deterministic, so two dry runs over the same inputs produce the same bytes.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from astroreg_register import screen_pairs  # noqa: E402

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import astrorun, crispri  # noqa: E402

REGISTRATION = Path("data/results/astroreg_registration.json")
REGISTRATION_SHA256 = "43f2edfab0636e0f07a035aedf5ddf7c57aa1dfcee8d39ea32e9338edd522c7a"

#: The gate whose verdict governs whether the list below may be bought at all.
GATE_RESULT = Path("data/results/astroreg_activity_nogo.json")

#: Where the full list is written for review.
DEFAULT_OUT = Path("data/results/astroreg_request_plan.json")

#: The labels the registered test actually scores. Elements that serve only other labels are counted
#: and reported, never removed: the authorised total counts them.
SCORED_LABELS = ("positive", "negative")


def refuse(message: str) -> int:
    print(f"REFUSED: {message}")
    return 2


def check_the_registration() -> dict[str, Any] | None:
    digest, _, _ = mf.sha256_of(REGISTRATION)
    if digest != REGISTRATION_SHA256:
        print(f"REFUSED: {REGISTRATION} is sha256 {digest}, not the frozen {REGISTRATION_SHA256}")
        return None
    return json.loads(REGISTRATION.read_text())


def the_gates_verdict() -> dict[str, Any] | None:
    """The committed activity gate's verdict, quoted and not decided again."""
    if not GATE_RESULT.exists():
        return None
    g = json.loads(GATE_RESULT.read_text())
    return {
        "result": str(GATE_RESULT),
        "go": g["verdict"]["go"],
        "columns_not_in_frozen_units": g["verdict"]["columns_not_in_frozen_units"],
        "requests_permitted": g["requests_permitted"],
        "registered_rule": g["verdict"]["registered_rule"],
        "reading": g["verdict"]["reading"],
        "decision": g["decision"],
    }


def describe(plan: list[dict[str, Any]], pairs: list[dict[str, Any]]) -> dict[str, Any]:
    """What the list is made of, in the terms a reviewer decides in."""
    per_chrom = Counter(r["chrom"] for r in plan)
    label_sets = Counter()
    serves_a_scored_pair = 0
    for r in plan:
        label_sets[tuple(r["serves_labels"])] += 1
        if any(lab in SCORED_LABELS for lab in r["serves_labels"]):
            serves_a_scored_pair += 1
    pair_labels = Counter(p["label"] for p in pairs)
    return {
        "requests": len(plan),
        "chromosomes": len(per_chrom),
        "requests_per_chromosome": {c: per_chrom[c] for c in sorted(per_chrom, key=lambda c: (len(c), c))},
        "screen_pairs_read": len(pairs),
        "screen_pairs_by_label": dict(sorted(pair_labels.items())),
        "requests_serving_at_least_one_scored_pair": serves_a_scored_pair,
        "requests_serving_no_scored_pair": len(plan) - serves_a_scored_pair,
        "label_combinations_served": {
            "|".join(k): v for k, v in sorted(label_sets.items(), key=lambda kv: -kv[1])
        },
        "elements_are_unique": len({(r["chrom"], r["element"]) for r in plan}) == len(plan),
        "costing_observation": (
            f"{len(plan) - serves_a_scored_pair} of the {len(plan)} requests serve no pair the "
            f"registered test scores (its scored labels are {', '.join(SCORED_LABELS)}): every pair "
            "they cover is one of the 25 increases held apart or one of the 3,154 non-hits that are "
            "not WellPowered at fc 0.25 and are excluded. They are counted here because the "
            "AUTHORISED total counts them -- the registration's own rule is one request per registry "
            "element overlapping a covered pair of any label, and its left_undone already says 'a "
            "cheaper rule may exist and is not explored here'. Removing them would be a cheaper run "
            "than the one registered, which is a change to a registered term and not this script's "
            "to make"
        ),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="enumerate and print the request list. The only mode this script has; required, so "
        "that running it is always a deliberate statement that nothing is being sent",
    )
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT, help="where to write the full list")
    ap.add_argument("--show", type=int, default=10, help="how many rows to print in full")
    args = ap.parse_args()
    if not args.dry_run:
        return refuse(
            "this script has one mode and it must be asked for by name: pass --dry-run. There is no "
            "--send and no code path that sends; the flag is required so that the absence of a spend "
            "is something the operator states rather than something they assume"
        )

    t0 = time.time()
    frozen = check_the_registration()
    if frozen is None:
        return 2
    authorised = frozen["requests"]["total_requests_needed"]
    needed_total = frozen["coverage"]["registry_elements_needed_total"]
    one_request_per = frozen["requests"]["one_request_per"]
    if one_request_per != astrorun.ONE_REQUEST_PER:
        return refuse("the counting rule quoted in code is not the rule in the frozen registration")
    print(f"the frozen registration authorises {authorised} requests, one per {one_request_per}")

    print("reading the screen's pairs (no effect size, p-value or FDR is read)")
    pairs = screen_pairs()
    print(f"  {len(pairs)} pairs")

    print(f"enumerating the request list from {crispri.ELEMENTS} (this reads the element registry)")
    plan = astrorun.plan_requests(pairs, crispri.DeletionTable())
    summary = describe(plan, pairs)

    if len(plan) != authorised or len(plan) != needed_total:
        return refuse(
            f"the enumerated list holds {len(plan)} requests, but the registration authorises "
            f"{authorised} and its coverage table needed {needed_total}. A list that does not match "
            "the authorised number may not be sent: the difference is an unauthorised request or a "
            "missing one, and either way the enumeration and the registration disagree about what "
            "was approved"
        )
    print(f"  {len(plan)} requests, which matches the authorised {authorised} exactly")

    gate = the_gates_verdict()
    if gate is None:
        return refuse(
            f"the activity gate's result {GATE_RESULT} is not present. The list above is the list the "
            "run would send, but whether the registered claim can be computed from it is the gate's "
            "question, and an unreviewed list is not a reviewed one. Run "
            "scripts/astroreg_activity_gate.py first"
        )

    payload = {
        "status": "a dry run. NOTHING WAS SENT: this script has no AlphaGenome import and no code "
        "path that could send a request",
        "lane": "lane-astro",
        "authorisation": astrorun.AUTHORISATION,
        "authorised_requests": authorised,
        "requests_enumerated": len(plan),
        "requests_sent": 0,
        "money_spent": 0,
        "registration": {"path": str(REGISTRATION), "sha256": REGISTRATION_SHA256},
        "one_request_per": one_request_per,
        "one_enumeration_only": astrorun.ONE_ENUMERATION_ONLY,
        "enumerated_by": "genomeos.attribution.astrorun.plan_requests, the same function any sender "
        "must iterate",
        "summary": summary,
        "activity_gate": gate,
        "may_this_list_be_bought": gate["requests_permitted"] > 0,
        "verdict": (
            "the list is complete and matches the authorised number exactly, and it MAY NOT BE "
            f"BOUGHT: the activity gate permits {gate['requests_permitted']} requests, because the "
            "frozen model's activity term cannot be supplied for an astrocyte pair in the units its "
            "weights were fitted on, and the registration's own feasibility rule answers a missing "
            "input with a no-go. Buying it would produce a deletion cache from which the registered "
            "primary claim cannot be computed"
        )
        if gate["requests_permitted"] == 0
        else "the activity gate permits the run; the list above is what it would send",
        "element_registry": {
            "root": str(crispri.ELEMENTS),
            "note": "read by overlap only, one chromosome held at a time. Not hashed here: this is a "
            "review artefact and not a registry result, so it carries no manifest and no number of "
            "its own that a later result could cite",
        },
        "requests": plan,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=1) + "\n")

    print()
    print(f"first {args.show} of {len(plan)} requests:")
    for r in plan[: args.show]:
        print(
            f"  {r['chrom']}:{r['start']}-{r['end']}  {r['element']}  "
            f"genes={','.join(r['serves_genes'][:3])}{'...' if len(r['serves_genes']) > 3 else ''}  "
            f"labels={','.join(r['serves_labels'])}"
        )
    print()
    print(
        json.dumps(
            {
                "authorised_requests": authorised,
                "requests_enumerated": len(plan),
                "matches_the_authorisation": len(plan) == authorised,
                "requests_sent": 0,
                "money_spent": 0,
                "requests_serving_no_scored_pair": summary["requests_serving_no_scored_pair"],
                "activity_gate_requests_permitted": gate["requests_permitted"],
                "may_this_list_be_bought": payload["may_this_list_be_bought"],
            },
            indent=1,
        )
    )
    print(f"({time.time() - t0:.0f} s) -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
