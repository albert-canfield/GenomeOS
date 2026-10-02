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

It is written THROUGH `save_result`, with a manifest declaring every input it read, including the
per-chromosome registry tables reached through `crispri.DeletionTable`. The first version of this
script wrote the file into data/results itself and so put a file there with no manifest at all. That
was a breach of the result contract, and the contract is the thing that makes a number in this
project checkable, so it is fixed here rather than excused. Why the static guard on that contract did
not refuse it is recorded in the commit message, because a guard that exists and did not fire is
worth more attention than the one file it missed.
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

from astroreg_register import TABLE3, screen_pairs  # noqa: E402

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import astrorun, crispri  # noqa: E402
from genomeos.results import save_result  # noqa: E402

REGISTRATION = Path("data/results/astroreg_registration.json")
REGISTRATION_SHA256 = "43f2edfab0636e0f07a035aedf5ddf7c57aa1dfcee8d39ea32e9338edd522c7a"

#: The gate whose verdict governs whether the list below may be bought at all.
GATE_RESULT = Path("data/results/astroreg_activity_nogo.json")

#: The registry name this list is written under, through save_result and its manifest. The first
#: version of this script wrote the file itself, which put a file into data/results with no manifest
#: at all -- the exact breach the result contract exists to prevent. It goes through save_result now.
RESULT = "astroreg_request_plan"

#: This script, as the entry whose transitive import closure is the counting path of its result.
ENTRY = "scripts/astroreg_dry_run.py"

OWN_CODE = (
    "genomeos/attribution/astrorun.py",
    "scripts/astroreg_activity_gate.py",
    "scripts/astroreg_dry_run.py",
    "tests/test_astrorun.py",
)

#: The labels the registered test actually scores, from the one place that defines them.
SCORED_LABELS = astrorun.SCORED_LABELS

#: AstroREG-2's registration, whose scope is the 1,232 requests that serve a pair the test scores.
ASTROREG2 = Path("data/results/astroreg2_registration.json")

#: The number that scope enumerates. A list that does not match it is refused, not trimmed.
ASTROREG2_REQUESTS = 1_232


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


def manifest(chroms: list[str], plan: list[dict[str, Any]], pairs: list[dict[str, Any]]) -> dict[str, Any]:
    """What this list was built from. Declared because every traced read must be declared."""
    element_tables = [crispri.ELEMENTS / f"{c}.json" for c in chroms]
    return {
        "sources": [
            {
                "accession": str(REGISTRATION),
                "version": f"the frozen registration, sha256 {REGISTRATION_SHA256}; its authorised "
                "total and its counting rule are read and checked, never moved",
            },
            {
                "accession": "Green NFO, et al. CRISPRi screening in cultured human astrocytes. "
                "Nature Neuroscience 2025;29(3):703-716 -- Supplementary Table 3",
                "version": "publisher open-access supplementary store; pinned here by sha256",
                "url": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-025-02154-3"
                "/MediaObjects/41593_2025_2154_MOESM5_ESM.xlsx",
            },
            {
                "accession": "the swept all-element deletion registry, "
                f"{crispri.ELEMENTS} ({len(element_tables)} per-chromosome tables)",
                "version": "this project's own sweep; read by overlap only, and no deletion value is "
                "read from it here",
            },
        ],
        "inputs": [
            mf.input_entry(
                REGISTRATION,
                partition=None,
                role="the authorised request total, the coverage table's "
                "registry_elements_needed_total and the counting rule; its sha256 is checked and no "
                "term is moved",
            ),
            mf.input_entry(
                GATE_RESULT,
                partition=None,
                role="the activity gate's committed verdict, quoted rather than decided again: it is "
                "what says whether this list may be bought",
            ),
            mf.input_entry(
                TABLE3,
                partition=None,
                role="the screen's element coordinates and registered labels; no effect size, "
                "p-value, FDR or score was read",
            ),
            mf.files_entry(
                f"{crispri.ELEMENTS}/chr*.json",
                element_tables,
                partition=None,
                role="the swept registry element spans, for the overlap test that decides which "
                "elements a covered pair needs. Element ids and coordinates only; no deletion value",
            ),
        ],
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "coordinates_note": (
            "the registry tables and the overlap test are zero-based half-open. The screen states no "
            "convention for its chrom:start-end strings; the only use made of them here is the same "
            "overlap test crispri.annotate uses for `covered`, so this list is built by the condition "
            "the registration's own coverage table counted with"
        ),
        "parameters": {
            "authorised_requests": astrorun.AUTHORISED_REQUESTS,
            "requests_enumerated": len(plan),
            "one_request_per": astrorun.ONE_REQUEST_PER,
            "chromosomes": chroms,
            "screen_pairs": len(pairs),
            "alphagenome_requests": 0,
            "model_requests": 0,
            "requests_sent": 0,
            "money_spent": 0,
        },
        "exclusions": [
            "no AlphaGenome request is sent: this script has no such import and no code path that "
            "could send one, which tests/test_astrorun.py checks against its syntax tree",
            "no deletion value is read from the registry: only each element's id and span",
            "no effect size, fold change, p-value, FDR or expression level is read from the screen",
            "no astrocyte activity value is read or constructed: whether one can be is the gate's "
            "question and this list does not reopen it",
        ],
        "partitions": {
            "screen_pairs": f"all {len(pairs)} element-gene pairs of the published screen, every "
            "label, because the authorised total counts a covered pair of any label",
        },
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def astroreg2_scope(plan, pairs, summary, authorised, show) -> int:
    """The 1,232 AstroREG-2 would send: the single enumeration, filtered in its one filter."""
    if not ASTROREG2.exists():
        return refuse(f"{ASTROREG2} is absent; AstroREG-2's scope is read from its registration")
    reg_sha, _, _ = mf.sha256_of(ASTROREG2)
    sendable = astrorun.requests_serving_scored_pairs(plan)
    dropped = len(plan) - len(sendable)
    if len(sendable) != ASTROREG2_REQUESTS:
        return refuse(
            f"the sendable list holds {len(sendable)} requests, not the {ASTROREG2_REQUESTS} this "
            "scope enumerates. A list that does not match the number under review may not be sent: "
            "the difference is an unauthorised request or a missing one, and either way the "
            "enumeration and the scope disagree"
        )
    print(f"  {len(sendable)} requests serve a pair the test scores; {dropped} do not and are out")

    payload = {
        "status": "a dry run for AstroREG-2's scope. NOTHING WAS SENT: no AlphaGenome import, no "
        "network, and no code path that could send a request",
        "lane": "lane-astro",
        "scope": "astroreg2-1232",
        "why_1232_and_not_1322": (
            f"the registered counting rule enumerates {authorised} requests, one per registry element "
            f"overlapping a covered pair of ANY label. {dropped} of those serve NO pair the registered "
            f"test scores: every pair they cover is one of the 25 increases held apart or one of the "
            f"3,154 non-hits excluded for not being WellPowered at fc 0.25, and neither class enters "
            f"either arm of the endpoint. Buying them buys nothing the claim can use, so AstroREG-2's "
            f"scope is the remaining {len(sendable)}. This is stated here rather than left for a "
            f"reader to subtract"
        ),
        "authorised_total_is_unchanged": "the ORIGINAL registration's authorised total is "
        f"{authorised} and narrowing it is a cheaper run than the one approved. Only Albert can "
        "authorise that, and his approval named the original registration by hash, so AstroREG-2 "
        "needs its own approval whichever number it names",
        "requests_enumerated_any_label": len(plan),
        "requests_sendable": len(sendable),
        "requests_dropped": dropped,
        "requests_sent": 0,
        "money_spent": 0,
        "registration": {"path": str(ASTROREG2), "sha256": reg_sha, "committed_at": "0f4c372"},
        "one_request_per": astrorun.ONE_REQUEST_PER,
        "one_enumeration_only": astrorun.ONE_ENUMERATION_ONLY,
        "one_list_then_one_filter": astrorun.ONE_LIST_THEN_ONE_FILTER,
        "enumerated_by": "genomeos.attribution.astrorun.plan_requests, then filtered by "
        "astrorun.requests_serving_scored_pairs. A sender iterates the SAME filtered list, so what is "
        "reviewed here and what would be sent cannot drift",
        "summary": summary,
        "label_combinations_dropped": {
            "|".join(sorted(set(r["serves_labels"]))): 1
            for r in plan
            if not any(lab in SCORED_LABELS for lab in r["serves_labels"])
        },
        "requests": sendable,
    }
    man = manifest(sorted({r["chrom"] for r in sendable}, key=lambda c: (len(c), c)), sendable, pairs)
    man["inputs"].append(
        mf.input_entry(
            ASTROREG2,
            partition=None,
            role="AstroREG-2's registration: the scope is read from it rather than from a flag, and "
            "its sha256 is recorded so the list can be tied to the terms it was enumerated under",
        )
    )
    payload[mf.KEY] = man
    out = save_result("astroreg2_request_plan", payload)
    print()
    print(f"first {show} of {len(sendable)}:")
    for r in sendable[:show]:
        print(
            f"  {r['chrom']}:{r['start']}-{r['end']}  {r['element']}  labels={','.join(r['serves_labels'])}"
        )
    print()
    print(
        json.dumps(
            {
                "scope": "astroreg2-1232",
                "enumerated_any_label": len(plan),
                "sendable": len(sendable),
                "dropped": dropped,
                "matches_the_scope": len(sendable) == ASTROREG2_REQUESTS,
                "requests_sent": 0,
                "money_spent": 0,
            },
            indent=1,
        )
    )
    print(f"-> {out}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="enumerate and print the request list. The only mode this script has; required, so "
        "that running it is always a deliberate statement that nothing is being sent",
    )
    ap.add_argument("--show", type=int, default=10, help="how many rows to print in full")
    ap.add_argument(
        "--scope",
        choices=("authorised-1322", "astroreg2-1232"),
        default="authorised-1322",
        help="`authorised-1322` is the original registration's own total, one request per registry "
        "element overlapping a covered pair of ANY label. `astroreg2-1232` is the subset that serves "
        "at least one pair the registered test scores, which is what AstroREG-2 would send",
    )
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

    if args.scope == "astroreg2-1232":
        return astroreg2_scope(plan, pairs, summary, authorised, args.show)

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
            "note": "read by overlap only, one chromosome held at a time, for element ids and spans "
            "and no deletion value. Declared and hashed in this result's manifest like any other "
            "input: an earlier version of this script wrote the list into data/results itself, with "
            "no manifest at all, and that was a breach of the result contract rather than a choice",
        },
        "requests": plan,
    }
    payload[mf.KEY] = manifest(sorted({r["chrom"] for r in plan}, key=lambda c: (len(c), c)), plan, pairs)
    out = save_result(RESULT, payload)

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
    print(f"({time.time() - t0:.0f} s) -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
