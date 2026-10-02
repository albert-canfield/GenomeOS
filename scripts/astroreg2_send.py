"""The ONE authorised AstroREG-2 run: the only script in this project that may buy an astrocyte answer.

    uv run python scripts/astroreg2_send.py --send

IT REFUSES BY DEFAULT AND ON EVERY CLAUSE SEPARATELY. `astrorun.may_send` is called first and checks
each clause of Albert's approval on its own, quoting the clause it enforces, plus the supervisor's item
(f). Today it refuses: no authorisation is recorded, no sign-off is recorded. Recording one does not lift
the others.

ONE PATH, AND THE POINT IS WHAT IT EXCLUDES. Every request goes through
`enhancer_target.score_element`, which is the RECORDING path: it writes the gene-axis record adapter v2
added, including `gene_id`, and each named cell's value multiset before any threshold. A sender that
called the AlphaGenome client directly would buy 1,232 answers and record a summary of them, and adapter
v2 would have bought nothing. That is not a style preference, it is the difference between keeping and
discarding data that has been paid for, so `tests/test_astroreg2_send.py` monkeypatches the client so
that ONLY the recording path can satisfy it and plants the direct call to show the failure.

THE LEDGER LINE IS WRITTEN BEFORE THE REQUEST IS PERMITTED, by `RequestBudget.take`, which is the single
place a charge happens. At the end `mark_run_complete` makes the run a terminal state, so a second run
refuses by name rather than by an exhausted counter.

THE PILOT CHECKPOINT, registered before any send. After the first `PILOT_REQUESTS` answers, every one of
them must carry a gene name and a non-empty effects set. The fault this exists for is undetectable
otherwise: a response whose gene proto omits `name` yields zero effects SILENTLY, which looks exactly
like "no gene moved". Without the checkpoint that is up to 1,232 requests spent for nothing, with nothing
to show it; with it, the cost of finding out is at most ten.

AND THE PRICE OF THE CHECKPOINT, which must not be a surprise. Albert's clause is "One run", and a
completed run is terminal here. So A PILOT STOP CONSUMES HIS ONE RUN: the remaining requests cannot be
resumed on this approval. That is not a budget that ran out, it is a run that finished, and resuming is a
NEW authorisation in his own words rather than a continuation. `PILOT_STOP_CONSUMES_THE_ONE_RUN` says so
and the refusal a later attempt meets quotes it, so the code and what he was told cannot diverge.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import astrorun  # noqa: E402
from genomeos.predict import enhancer_target  # noqa: E402
from genomeos.predict.alphagenome_adapter import ALPHAGENOME_MODEL_VERSION  # noqa: E402

RESULT = "astroreg2_send"
ENTRY = "scripts/astroreg2_send.py"
OWN_CODE = (
    "genomeos/attribution/astrorun.py",
    "scripts/astroreg2_send.py",
    "tests/test_astroreg2_send.py",
)

PLAN = Path("data/results/astroreg2_request_plan.json")
REGISTRATION = Path("data/results/astroreg2_registration.json")
ACTIVITY = Path("data/results/astroreg2_astrocyte_activity.json")
#: The money audit trail. NOT under data/results, and the reason is not a technicality: a ledger is not
#: a RESULT. It is an append-only record of charges, each line written BEFORE its charge precisely so it
#: cannot be reconstructed from outcomes, whereas a result is a computed value written once through
#: save_result with a manifest. The project's own guard -- nothing under genomeos/ or scripts/ writes into
#: data/results except through save_result -- was failing on this file, and the right answer was to move
#: the file rather than to carve an exemption into a claim worth keeping absolute.
LEDGER = Path("data/ledgers/astroreg2.jsonl")

#: The digest of the list the supervisor reviewed, as committed at cc5b097. The sender rebuilds the list
#: and refuses on any difference, so the reviewed list and the sent list are the same object.
REVIEWED_DIGEST = "b96148cab7174cf10407aa07a2d26edf1c1972df5db2c00de3e682ccacf2b519"
REVIEWED_AT = "cc5b097"

#: How many answers the pilot checkpoint inspects before the run may continue.
PILOT_REQUESTS = 10

PILOT_RULE = (
    "after the first 10 requests, EVERY one of them must carry a gene name and a non-empty effects set. "
    "The fault this guards against is silent: a response whose gene proto omits `name` yields zero "
    "effects, which is indistinguishable from 'no gene moved', so up to 1,232 requests could be spent "
    "for nothing with nothing to show it. The checkpoint bounds the cost of finding out at ten"
)

PILOT_STOP_CONSUMES_THE_ONE_RUN = (
    "a pilot stop CONSUMES Albert's one run. His clause is 'One run' and a completed run is terminal "
    "here, so the requests not sent cannot be resumed on this approval. This is NOT a budget that ran "
    "out -- it is a run that finished early, and resuming it is a NEW authorisation in his own words, "
    "not a continuation. The trade is still the right one: ten requests to learn the responses are "
    "usable beats 1,232 to learn they were not. But it is a trade and it reads as one"
)


def reviewed_plan() -> list[dict[str, Any]]:
    """The reviewed list, rebuilt and digest-checked. Refuses on any difference."""
    body = json.loads(PLAN.read_text())
    plan = body["requests"]
    astrorun.check_is_the_reviewed_plan(plan, REVIEWED_DIGEST)
    return plan


def answer_is_usable_note(note: dict[str, Any]) -> dict[str, Any]:
    """Read a ledger answer note back as a usability verdict, so the pilot can span a resume."""
    return {
        "id": note.get("id"),
        "gene_name_present": bool(note.get("gene_name_present")),
        "effects": note.get("effects"),
        "effects_non_empty": bool(note.get("effects_non_empty")),
        "usable": bool(note.get("usable")),
    }


def answer_is_usable(hit: dict[str, Any]) -> dict[str, Any]:
    """Whether one answer carries a gene name and a non-empty effects set."""
    genes = hit.get("genes") or []
    named = [g for g in genes if (g.get("gene") or "").strip()]
    return {
        "id": hit.get("id"),
        "gene_name_present": bool(named),
        "effects": len(genes),
        "effects_non_empty": len(genes) > 0,
        "usable": bool(named) and len(genes) > 0,
    }


def send(
    plan: list[dict[str, Any]],
    budget: astrorun.RequestBudget,
    score: Any,
    pilot_requests: int = PILOT_REQUESTS,
) -> dict[str, Any]:
    """Iterate the reviewed list through the recording path, ledgering before each send.

    `score` is `enhancer_target.score_element` bound to its scorer and fetcher by the caller. It is a
    parameter so a test can substitute a recording-path double; it is NOT a place to substitute a direct
    client call, and the test that plants one shows the failure.
    """
    # (b) iterate the UNCHARGED requests, never the plan. On a fresh ledger these are the same list; on
    # a partial one this is what stops an element being paid for twice. A resume that restarted the plan
    # would re-pay for everything already charged AND, because the cap resumes from the ledger, run out
    # that many requests early -- costing more than was approved while scoring fewer elements.
    todo = astrorun.requests_not_yet_charged(plan, budget.ledger)
    already = len(plan) - len(todo)
    checks: list[dict[str, Any]] = []
    stopped: dict[str, Any] | None = None
    for row in todo:
        budget.take(
            chrom=row["chrom"],
            element=row["element"],
            start=row["start"],
            end=row["end"],
            genes=row.get("serves_genes"),
        )
        try:
            hit = score(chrom=row["chrom"], element_id=row["element"], start=row["start"], end=row["end"])
        except BaseException as e:  # noqa: BLE001 - recorded, then re-raised unchanged
            # SMALL FIX 1: a charged-but-unanswered element is VISIBLE in the ledger rather than merely
            # implied by a missing answer line. Written before the raise propagates.
            budget.note(
                event="error",
                element=row["element"],
                chrom=row["chrom"],
                error=f"{type(e).__name__}: {e}",
            )
            raise
        check = answer_is_usable(hit or {})
        budget.note(event="answer", **check)
        checks.append(check)
        # SMALL FIX 2: the pilot counts the answers OF THE RUN, read from the ledger, not of this
        # process. In memory it would re-pilot after a resume or skip the pilot, and a pilot stop is
        # terminal and consumes Albert's one run, so being wrong either way is expensive.
        answers_so_far = len(astrorun.answer_notes(budget.ledger))
        if answers_so_far == pilot_requests:
            pilot = [answer_is_usable_note(n) for n in astrorun.answer_notes(budget.ledger)]
            bad = [c for c in pilot if not c["usable"]]
            if bad:
                stopped = {
                    "at_request": answers_so_far,
                    "failing": bad,
                    "rule": PILOT_RULE,
                    "consequence": PILOT_STOP_CONSUMES_THE_ONE_RUN,
                }
                budget.note(event="pilot_stop", at_request=answers_so_far, failing=len(bad))
                break
    astrorun.mark_run_complete(
        budget.ledger,
        requests=budget.sent,
        stopped_at_pilot=bool(stopped),
        why=PILOT_STOP_CONSUMES_THE_ONE_RUN if stopped else "the run completed its list",
    )
    return {
        "requests_sent": budget.sent,
        "requests_charged_before_this_pass": already,
        "requests_this_pass": len(checks),
        "answers_checked": len(astrorun.answer_notes(budget.ledger)),
        "pilot": {
            "requests": pilot_requests,
            "rule": PILOT_RULE,
            "passed": stopped is None,
            "stopped": stopped,
        },
        "unusable_answers": [c for c in checks if not c["usable"]],
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--send",
        action="store_true",
        help="required by name. Without it nothing is sent and the clause checks are reported only",
    )
    args = ap.parse_args()
    t0 = time.time()

    plan = reviewed_plan()
    print(f"the reviewed list rebuilt and digest-checked against {REVIEWED_AT}: {len(plan)} requests")

    try:
        permitted = astrorun.may_send(
            activity_result=ACTIVITY,
            registration=REGISTRATION,
            ledger=LEDGER,
            signoff=astrorun.SUPERVISOR_SIGNOFF,
            plan=plan,
            reviewed_digest=REVIEWED_DIGEST,
            committed=committed_in_git,
        )
    except (astrorun.NoAuthorisationError, astrorun.SendRefusedError) as e:
        print(f"REFUSED: {e}")
        return 2
    print(f"every clause satisfied: {permitted}")

    if not args.send:
        print("--send was not given, so nothing is sent. The clause checks above all passed")
        return 0

    budget = astrorun.astroreg2_budget(LEDGER)

    # Built INLINE and exactly as the sweep that produced the frozen K562 values built it
    # (scripts/enhancer_targets_all.py worker_scorer, verified at :98-107): the pinned client with the
    # same timeout, and the scorer with threshold=0.0. The adapter's default is 0.05, and taking it would
    # have censored every |effect| < 0.05 -- values the frozen cache CARRIES -- making the astrocyte
    # deletion term a DIFFERENT FEATURE from the one the weights were frozen on. No helper is written for
    # this: the line that failed before named a function that did not exist, and inventing it now would
    # make the name real while leaving the habit intact.
    from genomeos.genome import IndexedGenome, reference_fasta
    from genomeos.predict import AlphaGenomeAdapter
    from genomeos.predict.alphagenome_adapter import create_client

    adapter = AlphaGenomeAdapter()
    if not adapter.available:
        print(
            "REFUSED: no AlphaGenome key is available, so no request can be made. Reported rather than "
            "defaulted: a sender that proceeded without one would look like a run that bought nothing"
        )
        return 2
    adapter._client = create_client(  # noqa: SLF001 - the pinned version, as the sweep does
        adapter.api_key, timeout=astrorun.FROZEN_CLIENT_TIMEOUT
    )
    scorer = adapter._live_scorer(  # noqa: SLF001 - the only route the adapter exposes
        threshold=astrorun.FROZEN_SCORER_THRESHOLD
    )
    model_version = dict(getattr(adapter, "last_scan", {}) or {})

    genomes: dict[str, Any] = {}

    def fetch(locus: Any) -> str:
        genome = genomes.get(locus.chrom)
        if genome is None:
            fasta = reference_fasta(locus.chrom)
            if not fasta.exists():
                raise FileNotFoundError(
                    f"{locus.chrom} has no reference fasta at {fasta}; run `genomeos data fetch --chrom "
                    f"{locus.chrom}` first. Not defaulted: a sequence this sender invented would be bought"
                )
            genome = genomes[locus.chrom] = IndexedGenome(str(fasta))
        return genome.fetch(locus)

    def score(**kw: Any) -> dict[str, Any]:
        return enhancer_target.score_element(scorer, fetch, **kw)

    out = send(plan, budget, score)
    out["status"] = "the ONE authorised AstroREG-2 run"
    out["lane"] = "lane-astro"
    out["path"] = (
        "every request went through enhancer_target.score_element, the recording path that writes the "
        "gene-axis record and each named cell's value multiset. No client call was made from here"
    )
    out["ledger"] = str(LEDGER)
    out["feature_definition"] = {
        "scorer_threshold": astrorun.FROZEN_SCORER_THRESHOLD,
        "client_timeout": astrorun.FROZEN_CLIENT_TIMEOUT,
        "model_version": ALPHAGENOME_MODEL_VERSION,
        "last_scan": model_version,
        "why_recorded": astrorun.FROZEN_SCORER_PARAMETERS,
        "note": "recorded beside the answers so a later reader can tell WHICH feature definition "
        "produced them without reading the code. A threshold of 0.05 where the frozen feature used 0.0 "
        "is a changed feature definition, which the registration forbids",
    }
    out["charged_in_the_ledger"] = budget.charged()
    out[mf.KEY] = manifest(out, plan)
    from genomeos.results import save_result

    print(f"-> {save_result(RESULT, out)}")
    print(f"({time.time() - t0:.0f} s) sent {budget.sent} of {astrorun.ASTROREG2_CAP}")
    return 0


def committed_in_git(path: Path) -> bool:
    """Whether a path is committed and unmodified. Injected into may_send so it can be tested."""
    import subprocess

    r = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", str(path)], cwd=str(ROOT), check=False)
    if r.returncode != 0:
        return False
    tracked = subprocess.run(
        ["git", "ls-files", "--error-unmatch", str(path)],
        cwd=str(ROOT),
        capture_output=True,
        check=False,
    )
    return tracked.returncode == 0


def manifest(out: dict[str, Any], plan: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "sources": [
            {
                "accession": "AlphaGenome, through genomeos.predict.enhancer_target.score_element",
                "version": "the recording path; the gene-axis record and per-cell multiset of adapter "
                "v2 are written for every answer",
            }
        ],
        "inputs": [
            mf.input_entry(PLAN, partition=None, role="the reviewed request list, digest-checked"),
            mf.input_entry(REGISTRATION, partition=None, role="AstroREG-2's registered terms"),
            mf.input_entry(
                ACTIVITY,
                partition=None,
                role="the astrocyte activity columns, whose presence and rule are a clause of the approval",
            ),
        ],
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "astroreg2_cap": astrorun.ASTROREG2_CAP,
            "requests_in_the_reviewed_list": len(plan),
            "requests_sent": out["requests_sent"],
            "pilot_requests": PILOT_REQUESTS,
            "pilot_passed": out["pilot"]["passed"],
            "reviewed_digest": REVIEWED_DIGEST,
        },
        "exclusions": [
            "no request was made outside enhancer_target.score_element",
            "no request was sent without a ledger line written first",
            "no second run is possible: the ledger records a completed run",
        ],
        "partitions": {"astrocyte": "the reviewed list only"},
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


if __name__ == "__main__":
    raise SystemExit(main())
