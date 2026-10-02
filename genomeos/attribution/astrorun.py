"""The authorised AstroREG run: the money guard, and the gate on the one input the frozen model
needs and the astrocyte side cannot supply.

Three separate things live here on purpose: the money guard, the request list, and the gate on the
one input the frozen model needs. They are separate because each can be wrong on its own, and a
single object that did all three would hide which one had failed.

THE MONEY GUARD. Albert authorised a NUMBER, 1,322 AlphaGenome requests, not a job: "I approve the
1,322 AlphaGenome requests for the AstroREG test as registered (astroreg_registration.json at
e58bcee, exploratory, prevalence-matched power 0.79). Run it once, under the registration, with
every request logged". So the cap is a constant checked before every request is sent, and reaching
it is a REFUSAL that stops the run (`CapRefusedError`), never a warning and never a skipped element: a
run that silently stopped at the cap and reported a result would have chosen its population by how
much money was left. `RequestBudget.take` is the single place a request is charged, and it writes
the ledger line BEFORE it hands back the permission to send, so a request cannot be paid for
without being logged -- the log is the charge, not a side effect of it.

AND THE BUDGET IS RECONSTRUCTED FROM THE LEDGER, NEVER FROM ZERO. The session that registered this
test was killed before it sent a single request, which is the ordinary way a run ends: not with an
error it can catch but with the process gone. A budget that counted from zero on startup would let
the next session spend the whole authorised number again, so the cap would bound a process and not
the money. The ledger is the only durable record of a charge, so `RequestBudget.__init__` counts the
charges already in it and starts there; a restart's remaining cap is `cap` minus what the ledger
holds, and a ledger already holding `cap` charges refuses the very first request of the new process.
That is what `A_RESTART_MAY_NOT_DOUBLE_SPEND` says and what the resume tests prove.

THE REQUEST LIST. `plan_requests` is the only enumeration of what would be sent: one row per registry
element overlapping a covered screen pair, which is the rule the registration costed the run with.
The dry run prints it and any sender must iterate it, so the list reviewed before the money is spent
and the list actually sent cannot drift apart. `ONE_ENUMERATION_ONLY` says why that is structural
rather than a matter of care.

THE ACTIVITY-INPUT GATE. The frozen model is a difference between two arms, 'activity + distance +
deletion' and 'activity + distance', and `activity` is in both. `crispri.Pair.dhs` and
`crispri.Pair.h3k27ac` are read in exactly one place in this project, `crispri.parse`, from the
EPCrisprBenchmark columns `DHS.RPM` and `H3K27ac.RPM`: reads per million mapped reads in the
element. The benchmark publishes them for every pair of all five of its held-out cell types, which
is why the HCT116 arm could be scored on the frozen path at all. The astrocyte screen is not in the
benchmark and its Supplementary Table 3 has no activity sheet, so nothing supplies those two
numbers for an astrocyte pair.

What IS on disk for astrocyte is the registration's own two activity inputs -- ENCFF874OPW DNase
peaks and ENCFF970DKF H3K27ac replicated peaks -- and both are narrowPeak PEAK CALLS whose third
column is `signalValue`, not RPM. A peak call also has no value at all outside a called peak,
whereas the benchmark's RPM is defined for every element and reaches zero continuously. Putting a
signalValue where a weight fitted on an RPM is applied is a substitute for a missing input, and the
registration's feasibility rule decides that case before any money is spent:

    "if an input the frozen model needs is missing, the model is not applied and no substitute is
    put in its place: a no-go is recorded instead"

So this module's gate is a reading of a registered rule, not a new rule. It is deliberately a pure
function over measured availability, so the verdict can be tested without the files and so the
evidence it rests on has to be produced before it will answer.
"""

from __future__ import annotations

import json
import shutil
import sys
import threading
import time
from pathlib import Path
from typing import Any

# --------------------------------------------------------------------------- the money guard

#: Albert's authorisation, 2026-10-02, word for word. Quoted and never summarised, because the
#: number in it is the whole permission.
AUTHORISATION = (
    "I approve the 1,322 AlphaGenome requests for the AstroREG test as registered "
    "(astroreg_registration.json at e58bcee, exploratory, prevalence-matched power 0.79). Run it "
    "once, under the registration, with every request logged."
)

#: The authorised number of AlphaGenome requests. 1,322, not 1,323: the registration's
#: `requests.total_requests_needed`, which is one request per registry element overlapping a covered
#: pair. A number was approved and a number is what may be spent.
AUTHORISED_REQUESTS = 1_322

#: Why reaching the cap stops the run instead of trimming the work to fit.
CAP_IS_A_REFUSAL = (
    "the cap is a refusal, not a budget to be filled: a run that reached 1,322 and carried on with "
    "the elements it had already paid for would have chosen its population by how much money was "
    "left, which is choosing it after seeing which part of the genome was affordable. So the "
    "request that would be the 1,323rd raises CapRefusedError, the run stops, and the count reached and "
    "the reason are what gets reported"
)

#: Why the ledger line is written before the request is permitted rather than after it is answered.
LOG_IS_THE_CHARGE = (
    "the ledger line is written inside the same lock that increments the count, before the caller "
    "is given permission to send, so a request cannot be charged without being logged. A log "
    "written after the answer came back would lose exactly the requests that cost money and "
    "returned nothing"
)


#: Why a new process starts from the ledger's count and not from zero. The brief's own words for the
#: failure this prevents: a lane that can be killed and restarted into a second full spend is worse
#: than one that never runs.
A_RESTART_MAY_NOT_DOUBLE_SPEND = (
    "the cap must bound the money, not the process. A killed run leaves no exception to catch and no "
    "in-memory count to carry, only the ledger it already wrote, so the budget is reconstructed from "
    "the ledger on startup: a new RequestBudget over a ledger holding n charges has n spent and "
    "cap - n remaining, and a ledger holding cap charges refuses the first request of the new "
    "process. A truncated final line -- a process killed mid-write -- is counted AS A CHARGE, "
    "because the line is written before the request is permitted, so a half-written line is "
    "evidence that a charge was being made and the money-safe reading of it is that it was"
)

#: What the ledger does NOT protect, stated here so the guard is not read as more than it is. This is
#: a limitation of the design and not a defect to be fixed by a bigger lock.
WHAT_THE_LEDGER_CANNOT_DO = (
    "the lock makes the count and the ledger line atomic WITHIN one process. It does not make them "
    "atomic ACROSS two. Two processes started against the same ledger would each read the same count, "
    "each believe the same request unspent, and both send it, and no amount of reading the ledger "
    "closes that: the read and the append are two operations and another process can act between "
    "them. So the ledger bounds a RESTART, which is one process after another, and not CONCURRENCY, "
    "which is two at once. The only thing that closes the concurrent case is one executor holding the "
    "budget, which is an operating rule and not a property of this code. A guard that claimed "
    "otherwise would be the more dangerous for being trusted"
)


def ledger_charges(ledger: Path | str) -> dict[str, int]:
    """What the ledger on disk says has already been charged. The count a restart resumes from.

    `charges` counts every `event: request` line plus every line too damaged to read, because the
    ledger line is written before the request is permitted: a line truncated by a kill is evidence
    that a charge was being made, and counting it is the only reading that cannot spend twice. The
    other counts are reported beside it so a damaged ledger is visible rather than silently absorbed.
    """
    path = Path(ledger)
    if not path.exists():
        return {"charges": 0, "request_lines": 0, "unreadable_lines": 0, "other_event_lines": 0}
    request_lines = unreadable = other = 0
    for line in path.read_text(errors="replace").splitlines():
        if not line.strip():
            continue
        try:
            row = json.loads(line)
        except ValueError:
            unreadable += 1
            continue
        if isinstance(row, dict) and row.get("event") == "request":
            request_lines += 1
        else:
            other += 1
    return {
        "charges": request_lines + unreadable,
        "request_lines": request_lines,
        "unreadable_lines": unreadable,
        "other_event_lines": other,
    }


#: A helper that makes a resume cheap is not a permission to resume. Kept beside the function so the
#: two cannot be separated.
A_RESUME_STILL_NEEDS_ALBERTS_WORD = (
    "`requests_not_yet_charged` exists so that IF a resume is ever authorised it cannot re-pay for an "
    "element the ledger already holds. It is not itself an authorisation: a run that failed part way "
    "is reported at the count it reached, and a second run is a second spend that needs Albert's word "
    "again. The function is the money-safe way to carry out a decision someone else has made, not the "
    "decision"
)


def charged_elements(ledger: Path | str) -> set[tuple[str, str]]:
    """(chrom, element) of every charge the ledger records, for a resume that must not re-pay.

    Refuses rather than guesses when the ledger holds a line it cannot read: a damaged line is a
    charge whose element is unknown, so the set of elements already paid for is unknown, and a
    remainder computed from an unknown set could re-send a request that was already bought.
    """
    path = Path(ledger)
    counts = ledger_charges(path)
    if counts["unreadable_lines"]:
        raise ValueError(
            f"{path} holds {counts['unreadable_lines']} line(s) that cannot be read, so which "
            "elements were already charged is unknown and no remainder may be computed from it: "
            f"{A_RESUME_STILL_NEEDS_ALBERTS_WORD}"
        )
    out: set[tuple[str, str]] = set()
    if not path.exists():
        return out
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        if isinstance(row, dict) and row.get("event") == "request":
            out.add((row.get("chrom"), row.get("element")))
    return out


def requests_not_yet_charged(plan: list[dict[str, Any]], ledger: Path | str) -> list[dict[str, Any]]:
    """The rows of `plan` the ledger holds no charge for, in the plan's own order.

    Read `A_RESUME_STILL_NEEDS_ALBERTS_WORD` before calling this.
    """
    already = charged_elements(ledger)
    return [r for r in plan if (r["chrom"], r["element"]) not in already]


class CapRefusedError(RuntimeError):
    """The authorised number of requests is reached. Raised instead of returning, so the run stops.

    The message carries the count reached and the cap, because those two numbers are the report.
    """


class RequestBudget:
    """The one place an AlphaGenome request is charged, and the only thing that may permit one.

    `take` refuses at `cap` by raising `CapRefusedError`. Every charge appends one JSON line to `ledger`
    before it returns, so the committed log and the money spent cannot disagree. Thread safe: the
    count and the ledger line are written under one lock, so two workers cannot both be handed the
    last permitted request.

    Construction RESUMES: the charges already in `ledger` are counted and spent before this object
    permits anything, so a restart cannot spend the authorised number a second time. `sent` is a
    floor a caller may add, not a replacement for the ledger's count -- the larger of the two is
    taken, never their sum, because a caller passing what it read from the ledger must not charge it
    twice over. A resume appends one `event: resume` line naming the count it found, so a restart is
    on the record the committed log makes.
    """

    def __init__(self, ledger: Path | str, cap: int = AUTHORISED_REQUESTS, sent: int = 0) -> None:
        if cap < 0:
            raise ValueError(f"a request cap cannot be negative: {cap}")
        self.cap = cap
        self.ledger = Path(ledger)
        self.on_disk = ledger_charges(self.ledger)
        self.resumed_from_ledger = self.on_disk["charges"]
        self.sent = max(sent, self.resumed_from_ledger)
        self._lock = threading.Lock()
        if self.resumed_from_ledger:
            self._write(
                {
                    "event": "resume",
                    "charges_found_in_the_ledger": self.resumed_from_ledger,
                    **{k: v for k, v in self.on_disk.items() if k != "charges"},
                    "sent": self.sent,
                    "cap": self.cap,
                    "remaining": self.cap - self.sent,
                    "why": A_RESTART_MAY_NOT_DOUBLE_SPEND,
                    "t": time.strftime("%Y-%m-%dT%H:%M:%S"),
                }
            )

    def remaining(self) -> int:
        return self.cap - self.sent

    def take(self, **row: Any) -> int:
        """Charge one request and log it, or refuse. Returns the running total after the charge.

        `row` is what was asked for: the registry element, its coordinates, and whatever else the
        caller wants on the record. The time, the running total and the cap are added here so no
        caller can omit them.
        """
        with self._lock:
            if self.sent >= self.cap:
                raise CapRefusedError(
                    f"the authorised number of AlphaGenome requests is reached: {self.sent} of "
                    f"{self.cap} sent, and this request would be number {self.sent + 1}. "
                    f"{CAP_IS_A_REFUSAL}"
                )
            self.sent += 1
            self._write(
                {
                    **row,
                    "event": "request",
                    "n": self.sent,
                    "cap": self.cap,
                    "t": time.strftime("%Y-%m-%dT%H:%M:%S"),
                }
            )
            return self.sent

    def note(self, **row: Any) -> None:
        """Record something that is not a charge: an answer, a cache hit, a quota refusal, an error.

        Kept apart from `take` so the ledger's `event: request` lines are exactly the requests paid
        for and can be counted without a rule about which other events to subtract.
        """
        self._write({**row, "sent_so_far": self.sent, "t": time.strftime("%Y-%m-%dT%H:%M:%S")})

    def _write(self, row: dict[str, Any]) -> None:
        self.ledger.parent.mkdir(parents=True, exist_ok=True)
        with self.ledger.open("a") as fh:
            fh.write(json.dumps(row) + "\n")

    def charged(self) -> int:
        """The requests the ledger on disk says were charged. The count that is reported.

        Read from the file and not from `self.sent`, so the number reported is the number the
        committed artefact can be checked against.
        """
        return ledger_charges(self.ledger)["charges"]


# ------------------------------------------------------------------- the request list, once

#: The registration's own counting rule, quoted from requests.one_request_per. The list below is
#: this rule executed, and the registration's total is what it is checked against.
ONE_REQUEST_PER = "registry element overlapping a covered pair, as the HCT116 arm counted it"

#: Why the enumeration lives here and not in the dry run and again in a sender.
ONE_ENUMERATION_ONLY = (
    "the dry run exists so the list can be reviewed before it is paid for, and a review is worth "
    "nothing if the thing reviewed is not the thing sent. So the list is built in ONE function, "
    "`plan_requests`, which the dry run prints and which any sender must iterate: not a second "
    "enumeration written to match, because two enumerations drift and the drift is only visible "
    "after the money is gone. The list is deterministic and sorted, so two runs of the dry run on "
    "the same inputs produce the same list in the same order and a reviewer can diff them"
)


def plan_requests(pairs: list[dict[str, Any]], table: Any) -> list[dict[str, Any]]:
    """THE request list: one row per registry element overlapping a covered screen pair.

    `pairs` are the screen's pairs as `astroreg_register.screen_pairs` returns them, each with
    `chrom`, `start`, `end`, `element`, `gene` and `label`. `table` is a `crispri.DeletionTable`, and
    the overlap condition is its own `overlapping`, which is the same condition `crispri.annotate`
    uses for `covered` and the same one the registration's coverage table counted with. A pair that
    overlaps nothing contributes nothing, which is why no row here can be uncovered.

    One row per registry element, NOT per pair: an element overlapping several pairs is one request,
    which is the rule the registration costed the run with (`ONE_REQUEST_PER`) and the rule the
    HCT116 arm counted by. Every pair the element serves is carried on the row, with its label, so a
    reviewer can see what each request is for before it is bought.

    Deliberately NOT filtered by label: the registered total counts elements overlapping a covered
    pair of ANY label, so filtering here would produce a different number from the one authorised.
    The share of elements that serve only labels the test never scores is reported by the dry run
    instead, where it is an observation about the costing and not a change to it.
    """
    rows: dict[tuple[str, str], dict[str, Any]] = {}
    for p in pairs:
        for e in table.overlapping(p["chrom"], p["start"], p["end"]):
            row = rows.setdefault(
                (p["chrom"], e["id"]),
                {
                    "chrom": p["chrom"],
                    "element": e["id"],
                    "start": e["start"],
                    "end": e["end"],
                    "serves_screen_elements": set(),
                    "serves_genes": set(),
                    "serves_labels": set(),
                },
            )
            row["serves_screen_elements"].add(p["element"])
            row["serves_genes"].add(p["gene"])
            row["serves_labels"].add(p["label"])
    out = []
    for key in sorted(rows, key=lambda k: (len(k[0]), k[0], rows[k]["start"], k[1])):
        row = rows[key]
        out.append(
            {
                **{k: v for k, v in row.items() if not isinstance(v, set)},
                "serves_screen_elements": sorted(row["serves_screen_elements"]),
                "serves_genes": sorted(row["serves_genes"]),
                "serves_labels": sorted(row["serves_labels"]),
            }
        )
    return out


#: The labels the registered test actually scores. The other two registered classes -- the 25 increases
#: held apart and the 3,154 non-hits excluded for not being WellPowered at fc 0.25 -- are never
#: negatives and never positives, so an element whose only covered pairs carry those labels cannot
#: enter either arm of the endpoint.
SCORED_LABELS = ("positive", "negative")

#: Why the sendable list is a FILTER of `plan_requests` and not a second enumeration.
ONE_LIST_THEN_ONE_FILTER = (
    "the sendable list is `plan_requests` filtered, in this one function, and a sender iterates THIS. "
    "Re-deriving it from the pairs would be a second enumeration, and two enumerations drift. The "
    "filter is also the only place the scored-label rule lives, so a report and a sender cannot "
    "disagree about which requests are in scope"
)


def requests_serving_scored_pairs(
    plan: list[dict[str, Any]], scored: tuple[str, ...] = SCORED_LABELS
) -> list[dict[str, Any]]:
    """The rows of `plan` that serve at least one pair the registered test scores, in the plan's order.

    An element whose every covered pair is an increase held apart or an underpowered exclusion serves
    neither arm of the endpoint, so buying it buys nothing the claim can use. This is a FILTER of the
    single enumeration and never a re-derivation of it; see `ONE_LIST_THEN_ONE_FILTER`.

    It does NOT change the authorised total. The original registration counted one request per registry
    element overlapping a covered pair of ANY label, and that total is what Albert approved. Narrowing
    it is a cheaper run than the one registered, which only he can authorise.
    """
    return [r for r in plan if any(lab in scored for lab in r["serves_labels"])]


# ------------------------------------------------- AstroREG-2: its own cap, its own authorisation

#: AstroREG-2's scope: the requests that serve a pair the registered test scores. NOT 1,322. The
#: original total counted a covered pair of any label; this is the subset either arm of the endpoint can
#: use, and it is a DIFFERENT number under a DIFFERENT registration.
ASTROREG2_CAP = 1_232

#: Where blob checks resolve paths from. The repository root, two levels above this module.
ROOT_FOR_BLOBS = Path(__file__).resolve().parents[2]

#: Albert's approval for AstroREG-2, to be recorded here VERBATIM in his own words when he gives them.
#: It is None because he has not given it. The approval already on record named the ORIGINAL
#: registration by hash and does not carry: a runner that reused it would spend on the strength of an
#: approval for a different registration.
ASTROREG2_AUTHORISATION = (
    "Albert: I approve 1,232 AlphaGenome requests for AstroREG-2 as registered "
    "(astroreg2_registration at 0f4c372): after the astrocyte activity is computed under Amendment 2's "
    'rule and the supervisor writes "dry run reviewed". One run, every request logged.'
)

#: WHEN HE WROTE IT, bounded by evidence rather than stated to a precision this record does not have.
#: His words reached the coordinator session on 2026-10-02 and were written BEFORE 15:15:13 +0100, which
#: is the commit time of 53a0332 -- the commit encoding his clauses as refusals, which could only be
#: written after they had been relayed. No finer timestamp is available to the session that received
#: them, and none is invented here: a precise time asserted without a source would be the same defect
#: as a relayed text transcribed as an approval.
ASTROREG2_AUTHORISATION_WRITTEN_BEFORE = "2026-10-02 15:15:13 +0100 (evidenced by 53a0332)"

#: Why this slot is filled by the coordinator and not by a lane, and why only now.
WHO_RECORDED_THE_AUTHORISATION = (
    "recorded by the coordinator session that received Albert's words FIRST-HAND, not by a lane. A lane "
    "holds a relay, and a relay is not a first-hand record -- see WHY_THE_RELAYED_TEXT_IS_NOT_THE_APPROVAL, "
    "and note that the lane correctly REFUSED to fill this slot from the relayed text. It is filled only "
    "now, after both of his conditions were met: the astrocyte activity committed at cf1b95f under "
    "Amendment 2's rule, and the supervisor's sign-off written at 17:23:05 +0100. Filling it earlier would "
    "have left the coordinator's own intention as the last barrier to spending, which is not a mechanism"
)

#: Why absence refuses instead of warning.
NO_AUTHORISATION_MEANS_NO_SEND = (
    "the AstroREG-2 scope REFUSES to send while no approval of its own is recorded. Not a warning, a "
    "refusal: a warning can be read past, and the failure it would permit is spending real money under "
    "an approval that covers a different registration. Requiring Albert's actual words to be physically "
    "present in the code before anything can send is the only version of this that cannot drift, "
    "because there is nothing to forget and nothing to infer"
)


class NoAuthorisationError(RuntimeError):
    """No approval is recorded for this scope, so nothing may be sent under it."""


def astroreg2_authorisation() -> str:
    """Albert's AstroREG-2 approval, verbatim, or a refusal. Never a default and never inferred."""
    if not ASTROREG2_AUTHORISATION:
        raise NoAuthorisationError(
            "no AstroREG-2 approval is recorded. The approval on record names the ORIGINAL "
            f"registration by hash and does not carry to it. {NO_AUTHORISATION_MEANS_NO_SEND}"
        )
    return ASTROREG2_AUTHORISATION


# ----------------------------------- one approval per run, and one budget across every run's ledger

A_GUARD_KEYED_TO_WHAT_THE_CALLER_SUPPLIES_IS_NOT_A_GUARD = (
    "the third appearance of one shape on this path: a mechanism that EXISTED and was not WIRED to the "
    "thing it protected. The resume path could re-pay for charged elements while "
    "requests_not_yet_charged sat unused; the blob signature could never be satisfied because a file "
    "held its own hash; and the one-run clause was enforced against a ledger THE RUN ITSELF CHOOSES. "
    "ONE_RUN asked 'is this ledger finished?' when the question is 'has Albert authorised THIS run?'. "
    "So the run id is now the only input: may_send DERIVES the ledger from it instead of accepting "
    "one, which removes the caller's ability to choose what it is judged against"
)

#: Albert's approvals KEYED BY RUN ID, because one text cannot answer a per-run question.
#:
#: The per-run ledger opened this: `--run N` resolves a FRESH, EMPTY file, so run_already_completed is
#: False and charges are 0, so the one-run clause PASSED for run 2, then 3, then 4, under run 1's one
#: consumed approval -- and the cap of 1,232 being per LEDGER left no total bound at all. Demonstrated
#: on 20fe5ee's body before it was changed: 4,928 requests permitted across four further runs with the
#: approval text unchanged, and no stopping point in sight. A per-ledger cap is not a budget.
ASTROREG2_AUTHORISATIONS: dict[int, dict[str, Any]] = {
    1: {
        "words": ASTROREG2_AUTHORISATION,
        "requests": ASTROREG2_CAP,
        "consumed": True,
        "why_consumed": (
            "run 1 happened under it: its ledger records a completed run, and 'One run' is terminal. "
            "So this approval authorises nothing further -- not a second run, and not a resume"
        ),
    },
}

AN_APPROVAL_IS_SPENT_WHEN_ITS_RUN_HAPPENS = (
    "an approval is for ONE run, so a completed run spends it. Marking it consumed in the record rather "
    "than inferring it from the ledger is deliberate: inference is what failed, because the ledger a "
    "later run reads is not the ledger the earlier run wrote. Run N's approval must be Albert's own "
    "words for run N, recorded verbatim by the session that received them FIRST-HAND, and a run with no "
    "approval of its own does not start"
)

A_TOTAL_IS_THE_ONLY_BOUND_ON_MONEY = (
    "the only figure that bounds spending is the sum of the authorised request counts across every run, "
    "checked against the charges in EVERY run ledger on disk. It is enforced twice on purpose: once as "
    "a refusal before a run starts, and once as the cap of the budget object itself, so a run cannot "
    "cross the total part way through. A ledger on disk for a run with no recorded approval refuses "
    "rather than being ignored, because charges outside the scheme are exactly what must not pass"
)


def ledger_path_for_run(run_id: int, root: Path | None = None) -> Path:
    """Run N's ledger as an absolute path. The run id picks the file; nothing else may."""
    base = Path(root) if root is not None else ROOT_FOR_BLOBS
    return base / ledger_for_run(run_id)


def run_ledgers(root: Path | None = None) -> dict[int, Path]:
    """Every run ledger present on disk, by run id, authorised or not."""
    import re as _re

    base = Path(root) if root is not None else ROOT_FOR_BLOBS
    found: dict[int, Path] = {}
    first = base / LEDGER_RUN1
    if first.exists():
        found[1] = first
    folder = base / LEDGER_RUN1.parent
    if folder.is_dir():
        for path in sorted(folder.glob("astroreg2_run*.jsonl")):
            m = _re.fullmatch(r"astroreg2_run(\d+)\.jsonl", path.name)
            if m:
                found[int(m.group(1))] = path
    return dict(sorted(found.items()))


def total_charged_across_runs(root: Path | None = None) -> dict[str, Any]:
    """Charges summed over EVERY run ledger, which is the only number that tracks money."""
    per = {run: ledger_charges(path)["charges"] for run, path in run_ledgers(root).items()}
    return {"per_run": per, "total": sum(per.values())}


def total_authorised_requests() -> int:
    """The sum of the authorised request counts, consumed or not: what has ever been approved in total."""
    return sum(int(rec["requests"]) for rec in ASTROREG2_AUTHORISATIONS.values())


def authorisation_for_run(run_id: int) -> dict[str, Any]:
    """Albert's approval for THIS run, verbatim, or a refusal. Never run 1's, never inferred."""
    rec = ASTROREG2_AUTHORISATIONS.get(int(run_id))
    if rec is None or not rec.get("words"):
        raise NoAuthorisationError(
            f"no approval is recorded for run {run_id}. {AN_APPROVAL_IS_SPENT_WHEN_ITS_RUN_HAPPENS}. "
            f"{NO_AUTHORISATION_MEANS_NO_SEND}"
        )
    if rec.get("consumed"):
        raise NoAuthorisationError(
            f"run {run_id}'s approval is CONSUMED: {rec.get('why_consumed')}. A further run needs a NEW "
            f"approval in Albert's own words, recorded for that run. "
            f"{AN_APPROVAL_IS_SPENT_WHEN_ITS_RUN_HAPPENS}"
        )
    return rec


def check_earlier_runs_are_complete(run_id: int, root: Path | None = None) -> dict[str, Any]:
    """Run N does not start until runs 1..N-1 each END in a run_complete line. No skipping, no parallel."""
    ended = {}
    for earlier in range(1, int(run_id)):
        path = ledger_path_for_run(earlier, root)
        if not run_already_completed(path):
            raise SendRefusedError(
                f"run {run_id} cannot start: run {earlier}'s ledger {ledger_for_run(earlier)} does not "
                f"end in a {RUN_COMPLETE_EVENT} line, so run {earlier} either never happened or is "
                "still going. Starting run "
                f"{run_id} now would skip a run or put two runs on the same approvals at once, and "
                "either way the charges of the run that is still open are not counted against anything"
            )
        ended[earlier] = str(ledger_for_run(earlier))
    return {"earlier_runs_complete": ended}


def check_total_cap(run_id: int, plan_size: int, root: Path | None = None) -> dict[str, Any]:
    """Refuse unless this run's list still fits inside the TOTAL authorised across every run."""
    spent = total_charged_across_runs(root)
    unknown = sorted(set(spent["per_run"]) - set(ASTROREG2_AUTHORISATIONS))
    if unknown:
        raise SendRefusedError(
            f"ledgers exist for run(s) {unknown} that no recorded approval covers, holding "
            f"{sum(spent['per_run'][r] for r in unknown)} charge(s). {A_TOTAL_IS_THE_ONLY_BOUND_ON_MONEY}"
        )
    authorised = total_authorised_requests()
    already = spent["total"] - spent["per_run"].get(int(run_id), 0)
    if already + plan_size > authorised:
        raise CapRefusedError(
            f"run {run_id}'s {plan_size} requests on top of {already} already charged in other runs "
            f"would come to {already + plan_size}, and the total ever authorised is {authorised}. "
            f"{A_TOTAL_IS_THE_ONLY_BOUND_ON_MONEY}"
        )
    return {
        "total_authorised": authorised,
        "charged_in_other_runs": already,
        "charged_per_run": spent["per_run"],
        "remaining_in_total": authorised - already,
    }


def astroreg2_budget(*, run_id: int, root: Path | None = None) -> RequestBudget:
    """A budget for run N, obtainable only once run N's OWN approval is recorded.

    No ledger is accepted: the run id picks it. The cap is the LOWER of this run's approved count and
    what the total leaves, so crossing the total is a refusal part way through a run and not only
    before it starts.
    """
    astroreg2_authorisation()
    rec = authorisation_for_run(run_id)
    totals = check_total_cap(run_id, 0, root)
    cap = min(int(rec["requests"]), int(totals["remaining_in_total"]))
    if cap <= 0:
        raise CapRefusedError(
            f"run {run_id} has no room: {totals['charged_in_other_runs']} of "
            f"{totals['total_authorised']} authorised requests are already charged in other runs. "
            f"{A_TOTAL_IS_THE_ONLY_BOUND_ON_MONEY}"
        )
    return RequestBudget(ledger_path_for_run(run_id, root), cap=cap)


#: The text of Albert's AstroREG-2 approval AS RELAYED to this lane, recorded as relayed and NOT as the
#: approval itself. It is kept apart from `ASTROREG2_AUTHORISATION` on purpose: see
#: `WHY_THE_RELAYED_TEXT_IS_NOT_THE_APPROVAL`.
ASTROREG2_AUTHORISATION_AS_RELAYED = (
    "Albert: I approve 1,232 AlphaGenome requests for AstroREG-2 as registered "
    "(astroreg2_registration at 0f4c372): after the astrocyte activity is computed under Amendment 2's "
    'rule and the supervisor writes "dry run reviewed". One run, every request logged.'
)

#: Why the relayed text does not populate the approval slot.
WHY_THE_RELAYED_TEXT_IS_NOT_THE_APPROVAL = (
    "this lane received the text above through another agent, not from Albert. An agent's account of what "
    "a person approved is not that person's approval, and transcribing it into the slot the send path "
    "reads would manufacture a record of consent this lane cannot verify. So the slot stays empty and the "
    "relayed text is kept beside it, labelled. Nothing is lost operationally: every other clause of the "
    "approval is encoded and refuses on its own, and the supervisor sign-off is absent too, so the send "
    "path refuses on several grounds at once. What is preserved is that the record does not assert "
    "something about Albert that only a relay supports"
)


#: Why the coordinator records this and not the lane, departing from the instruction it was given.
WHY_THE_COORDINATOR_RECORDED_THE_SIGNOFF = (
    "the supervisor asked the LANE to record its words. The coordinator recorded them instead, for the "
    "same reason the lane refused to record Albert's: the supervisor wrote them to the coordinator, so "
    "the coordinator holds them FIRST-HAND and a lane would hold only a relay. Applying the rule to the "
    "party that made it is not an exception to it. The divergence was reported to the supervisor rather "
    "than made silently, and either record would have been checkable -- but only one of them is a "
    "first-hand record"
)

#: Each clause of the approval, quoted, so a refusal can name the clause it enforces rather than saying
#: "not authorised", which teaches nothing and can be cleared by accident.
ASTROREG2_CLAUSES = {
    "scope": '"I approve 1,232 AlphaGenome requests for AstroREG-2 as registered '
    '(astroreg2_registration at 0f4c372)"',
    "activity": '"after the astrocyte activity is computed under Amendment 2\'s rule"',
    "signoff": '"and the supervisor writes \\"dry run reviewed\\""',
    "one_run": '"One run"',
    "logged": '"every request logged"',
}

#: The rule the activity result must have been produced under, as astroreg2_registration at 0f4c372
#: records it. Identified by CONTENT, not by a file existing.
AMENDMENT_2_RULE_FINGERPRINT = {
    "repository": "mayasheth/chrom-annotate",
    "commit": "91cda73ebe3a19153a582cab18cbf7ff70d85cfc",
    "file": "workflow/scripts/neighborhoods.py",
}

#: The registration the approval names, bound by content so its object cannot drift after the fact.
ASTROREG2_REGISTRATION_COMMIT = "0f4c372"

#: WHERE amendment 2's rule is actually registered, established by reading both files rather than
#: assumed. It is spread across two registrations and a reader should not be told it is in one.
AMENDMENT_2_RULE_REGISTERED_IN = {
    "data/results/astroreg2_registration.json at 0f4c372": "names the repository "
    "(mayasheth/chrom-annotate), the commit (91cda73ebe3a19153a582cab18cbf7ff70d85cfc) and the "
    "'NO RESCALING of any kind' term. It does NOT name the file, and it carries no term-by-term "
    "citation block",
    "data/results/astroreg_calibration_registration.json at 230efc8": "carries amendment_2_citations, "
    "the full term-by-term set with a line cited for each of the RPM formula, the numerator, the "
    "sex-chromosome doubling, the denominator and the several-BAM combination, plus the file name",
    "note": "the approval's clause says 'computed under Amendment 2's rule'. Two of the three "
    "fingerprint values the send path binds to are anchored in 0f4c372 and all three are in 230efc8, so "
    "the rule is registered across the pair and the send path's fingerprint is checked against both",
}

#: sha256 of `rpm.PRODUCER` as the registration at 230efc8 records it, verified equal term for term with
#: no key added or missing. Pinned so an edit to the cited rule fails a test instead of passing quietly.
AMENDMENT_2_CITATIONS_SHA256 = "77ff2a0a790d6b5baeacc10fa7ab2965a4ced14e6e64b5598401117c6c97c7db"

#: Run 1's ledger: the completed one, 23 lines including its correction row. Never written again.
LEDGER_RUN1 = Path("data/ledgers/astroreg2.jsonl")

A_NEW_RUN_GETS_ITS_OWN_LEDGER = (
    "a run's ledger is the record of THAT run, so a new run gets a new file rather than appending to a "
    "finished one. Sharing one file would mix two runs' charges in a single sequence and make 'how much "
    "did run 2 cost' unanswerable without a rule about where to cut. It also keeps run 1 exactly as it "
    "is -- 23 lines ending in its correction row -- which is what append-only is for"
)


def ledger_for_run(run_id: int) -> Path:
    """Run 1 keeps its own filename; every later run gets `astroreg2_run<N>.jsonl`."""
    if run_id < 1:
        raise ValueError(f"a run id is 1 or more: {run_id}")
    return LEDGER_RUN1 if run_id == 1 else Path(f"data/ledgers/astroreg2_run{run_id}.jsonl")


#: The ledger event that makes a completed run a TERMINAL state.
RUN_COMPLETE_EVENT = "run_complete"

ONE_RUN_IS_TERMINAL = (
    "'One run' is enforced as a terminal state and not as an exhausted counter. A counter at its cap and "
    "a run that has already happened are different facts, and a reset would look like the former while "
    "being the latter. So a finished run appends an `event: run_complete` line and the send path refuses "
    "BY NAME on finding one, whatever the count says"
)


#: Sign-off item (f). NOT one of Albert's clauses: it is the supervisor's requirement, and it is kept
#: apart so a refusal does not attribute to Albert a condition he did not state.
ADAPTER_V2_IS_A_SUPERVISOR_REQUIREMENT = (
    "item (f) of the supervisor's sign-off checklist, not a clause of Albert's approval. His approval "
    "names the activity precondition, the sign-off, one run, every request logged and the scope; this is "
    "an additional condition the supervisor imposed, and conflating the two would misreport what he "
    "agreed to"
)

WHY_ADAPTER_V2_BEFORE_SENDING = (
    "the prediction chain keeps a summary and discards what the model returned: the adapter reads only "
    "`gene_name` off the response's gene axis and the per-element writer keys by that bare string, so two "
    "response rows under one gene name MERGE. Buying 1,232 NEW answers through that adapter would destroy "
    "the gene identity and the per-track value multiset of data that has been PAID FOR -- permanently, "
    "unrecoverably, in exchange for a summary. Paid data is the one kind this project cannot re-fetch for "
    "free, so this is a refusal BEFORE the purchase rather than a repair after it"
)

#: The wording a merge must be reported in, carried verbatim from lane-generow.
A_MERGE_IS_A_LOST_DISTINCTION = "a merge is a LOST DISTINCTION, NOT A WRONG NUMBER"

#: What adapter v2 must expose for the send to proceed. Checked by CAPABILITY and not by a version
#: string, and default-refusing: if the building lane names these differently the refusal fires and says
#: exactly what it looked for, which is the safe direction to be wrong in.
ADAPTER_MODULE = "genomeos.predict.alphagenome_adapter"

#: Reconciled against adapter v2 as COMMITTED, after this check refused the names this lane had guessed.
#: That refusal was the design working: the guessed names were absent, the send refused, and the fix was
#: to read the adapter rather than to relax the check. GENE_AXIS_COLUMNS_090 is the gene-axis field set
#: taken from the client's own source; CELL_TRACKS and CELL_VALUES_NOTE are the per-cell value multiset
#: before any threshold; MERGE_NOTE is the merge audit that makes a pooled row countable.
ADAPTER_V2_MUST_EXPOSE = (
    "GENE_AXIS_SCHEMA",
    "GENE_AXIS_COLUMNS_090",
    "CELL_TRACKS",
    "CELL_VALUES_NOTE",
    "MERGE_NOTE",
)

#: The gene-axis field without which a purchased answer cannot be tied to a gene identity.
ADAPTER_V2_MUST_RECORD_GENE_ID = "gene_id"


def check_adapter_v2(module_name: str = ADAPTER_MODULE) -> dict[str, Any]:
    """Refuse unless the adapter the runner would write through is v2. See item (f).

    Looks for the capability rather than a declared version, because a version constant can be set
    without the behaviour existing. Refuses by naming each thing it looked for and did not find, so a
    naming disagreement with the lane building it surfaces as a refusal rather than as a silent pass.
    """
    import importlib

    try:
        mod = importlib.import_module(module_name)
    except ModuleNotFoundError as e:
        raise SendRefusedError(
            f"item (f): the adapter module {module_name} is not importable ({e}), so the runner cannot "
            f"be shown to write through adapter v2. {WHY_ADAPTER_V2_BEFORE_SENDING}"
        ) from e
    missing = [name for name in ADAPTER_V2_MUST_EXPOSE if not hasattr(mod, name)]
    if not missing and ADAPTER_V2_MUST_RECORD_GENE_ID not in tuple(getattr(mod, "GENE_AXIS_COLUMNS_090", ())):
        missing = [f"{ADAPTER_V2_MUST_RECORD_GENE_ID} in GENE_AXIS_COLUMNS_090"]
    if missing:
        raise SendRefusedError(
            f"item (f): {module_name} does not expose {missing}, so it is not adapter v2 and the gene "
            f"identity and per-track multiset of 1,232 paid answers would be discarded. "
            f"{WHY_ADAPTER_V2_BEFORE_SENDING} {ADAPTER_V2_IS_A_SUPERVISOR_REQUIREMENT}"
        )
    return {
        "module": module_name,
        "exposes": list(ADAPTER_V2_MUST_EXPOSE),
        "merge_reporting_rule": A_MERGE_IS_A_LOST_DISTINCTION,
    }


# --------------------------------------- the frozen feature's scorer parameters, shared by both sides

#: The threshold the FROZEN K562 deletion values were produced with. scripts/enhancer_targets_all.py's
#: worker_scorer passes `threshold=0.0` EXPLICITLY, so every gene in the window is kept, uncensored.
#: `AlphaGenomeAdapter._live_scorer` DEFAULTS to 0.05, and a sender taking that default would censor every
#: |effect| < 0.05 -- values the frozen cache CARRIES. That would make the astrocyte deletion term a
#: DIFFERENT FEATURE from the one the weights were frozen on, and the registration forbids exactly that:
#: "Nothing is refitted on astrocyte data -- not the weights, not a threshold, not a feature definition".
#: The gain would have been biased by construction and would have read like a result.
FROZEN_SCORER_THRESHOLD = 0.0

#: The client timeout the sweep gave its client, from scripts/enhancer_targets_all.py CALL_TIMEOUT.
FROZEN_CLIENT_TIMEOUT = 300

FROZEN_SCORER_PARAMETERS = (
    "the sweep that produced the frozen K562 values built its client as create_client(api_key, "
    "timeout=300) and its scorer as _live_scorer(threshold=0.0). Any sender buying values for the same "
    "frozen feature must pass the SAME arguments, so they live here as constants and a test reads BOTH "
    "call sites and asserts they agree. A default taken silently is how this diverged: no stub could "
    "catch it, because the stub substitutes the behaviour and the divergence was in a PARAMETER"
)

# --------------------------------------------------- the sign-off is bound to the code it was given for


A_SIGNOFF_IS_FOR_THE_CODE_IT_READ = (
    "the sign-off names a sender and a guard module by blob. If either file's bytes differ from the ones "
    "signed, the send refuses and says WHICH file differs and that a re-sign is required -- not that an "
    "authorisation is missing, which would be a different and misleading fact. This is what stops a "
    "sign-off carrying over to code the reviewer never read, and it bites in the ordinary course of work: "
    "the lane's own fixes invalidate it, which is the intended behaviour and not a defect"
)


def git_blob(path: Path | str, root: Path | None = None) -> str | None:
    """The git blob sha of a file AS IT IS ON DISK, so an uncommitted edit changes it too."""
    import subprocess

    p = Path(path)
    if not p.exists():
        return None
    r = subprocess.run(
        ["git", "hash-object", str(p)],
        cwd=str(root or Path.cwd()),
        capture_output=True,
        text=True,
        check=False,
    )
    return r.stdout.strip() or None


#: The sign-off record. An APPEND-ONLY FILE READ AT RUNTIME, never imported, and never in code.
#:
#: Two defects shaped this. The first version stored a file's own blob INSIDE that file, which is a
#: fixed point with no solution: writing the value changes the file, so no re-sign could ever pass.
#: Keeping the record outside the code dissolved that. The second was the FORMAT: a single JSON object
#: that a new sign-off REPLACED. Replacement destroys the history of what was signed and when, and it makes a
#: withdrawal indistinguishable from never having happened -- so it goes one record per LINE, and ONLY
#: THE LAST LINE GOVERNS.
SIGNOFF_RECORD = Path("data/ledgers/astroreg2_signoff.jsonl")

#: The superseded single-object record. Reverted to its committed bytes and NEVER WRITTEN AGAIN; it is
#: kept only so the history of what was signed is not erased by the move.
SIGNOFF_RECORD_SUPERSEDED = Path("data/ledgers/astroreg2_signoff.json")

A_SIGNOFF_IS_FOR_ONE_RUN = (
    "each sign-off line names the run_id it covers, and the send path requires the governing line to "
    "name THIS run. A sign-off that covered whatever run came next would be reusable, which is the same "
    "defect as a consumed approval authorising a second run: the closure it signed and the words it "
    "carries were written about one run's code and one run's spend"
)

ONLY_THE_LAST_LINE_GOVERNS = (
    "each line is one sign-off, carrying the closure it signed, the digest of that closure, the record "
    "it supersedes and why. ONLY THE LAST LINE AUTHORISES. An earlier line whose closure happens to "
    "match the current code must NOT authorise, because that is exactly how a WITHDRAWN sign-off would "
    "revive if the code were reverted -- and that is not hypothetical here: the 19:14 record signed the "
    "closure of a sender that would have charged 1,232 requests and bought nothing. A format that let "
    "that line speak again on a revert would reauthorise the defect"
)

APPEND_ONLY_IS_ENFORCED = (
    "append-only is a property of the HISTORY, and at send time what is checked is stronger: the record "
    "must be byte-IDENTICAL to its committed copy. A prefix rule is not enough, because HEAD being a "
    "prefix of the working file is exactly what an append produces -- so 'append-only' without "
    "'committed' authorises anybody who can write the file. Requiring identity makes the governing line "
    "a committed line, which is an act in the repository's history rather than an edit to a file"
)


def closure_sha256(closure: dict[str, str]) -> str:
    """The digest a sign-off line records for the closure it signed."""
    import hashlib

    return hashlib.sha256(json.dumps(closure, sort_keys=True).encode()).hexdigest()


def signoff_lines(path: Path | str | None = None) -> list[dict[str, Any]]:
    """Every sign-off record in order. Unreadable lines are kept as a marker, never skipped silently."""
    p = Path(path) if path is not None else ROOT_FOR_BLOBS / SIGNOFF_RECORD
    if not p.exists():
        return []
    out: list[dict[str, Any]] = []
    for line in p.read_text(errors="replace").splitlines():
        if not line.strip():
            continue
        try:
            row = json.loads(line)
        except ValueError:
            out.append({"_unreadable": line[:80]})
            continue
        out.append(row if isinstance(row, dict) else {"_unreadable": str(row)[:80]})
    return out


def git_blobs(paths: list[str], root: Path) -> dict[str, str]:
    """Blob shas for many files in ONE `git hash-object` call, keyed by the path given.

    One call rather than one per file: 52 subprocess spawns per check was both slow and a source of
    load-dependent failure, and a closure that changes under load is not a closure.
    """
    import subprocess

    if not paths:
        return {}
    r = subprocess.run(
        ["git", "hash-object", "--stdin-paths"],
        cwd=str(root),
        input="\n".join(paths) + "\n",
        capture_output=True,
        text=True,
        check=False,
    )
    lines = [x.strip() for x in r.stdout.splitlines() if x.strip()]
    if r.returncode != 0 or len(lines) != len(paths):
        return {}
    return dict(zip(paths, lines, strict=True))


def check_signoff_is_committed(path: Path | str | None = None, root: Path | None = None) -> dict[str, Any]:
    """Refuse unless the sign-off record is BYTE-IDENTICAL to its committed copy.

    Identical, not a prefix. HEAD being a prefix of the working file is exactly what an APPEND produces,
    so a prefix rule lets anyone who can write the file append a governing line and send. The governing
    line must be a COMMITTED line, so a sign-off is an act in the repository's history rather than a
    file someone has edited.

    Three defects in the first version each let an uncommitted sign-off authorise, and each is now a
    refusal: a missing HEAD copy returned "not checked" and was read as a pass; an appended uncommitted
    line satisfied the prefix; and an absolute path broke the `git show` lookup, which landed in the
    first branch, so a tampered file passed. A check that reports "not checked" and is treated as a pass
    is the vacuous-pass shape, and it is worse than no check because it reads as reassurance.
    """
    import subprocess

    base = (Path(root) if root is not None else ROOT_FOR_BLOBS).resolve()
    given = Path(path) if path is not None else SIGNOFF_RECORD
    target = (given if given.is_absolute() else base / given).resolve()
    try:
        rel = target.relative_to(base)
    except ValueError as e:
        raise SendRefusedError(
            f"the sign-off record {target} is outside the repository at {base}, so it is not something "
            "a commit can vouch for and it authorises nothing"
        ) from e

    committed = subprocess.run(
        ["git", "show", f"HEAD:{rel.as_posix()}"], cwd=str(base), capture_output=True, check=False
    )
    if committed.returncode != 0:
        raise SendRefusedError(
            f"the sign-off record {rel.as_posix()} has no committed copy at HEAD. "
            "An uncommitted sign-off authorises nothing: the governing line must be a committed line, "
            "which makes a sign-off an act in the repository's history rather than an edit to a file"
        )
    head_bytes = committed.stdout
    now = target.read_bytes() if target.exists() else b""
    if now != head_bytes:
        extra = len(now) - len(head_bytes)
        how = (
            f"{extra} uncommitted byte(s) beyond HEAD"
            if now.startswith(head_bytes)
            else "its committed bytes differ from the file on disk"
        )
        raise SendRefusedError(
            f"the sign-off record {rel.as_posix()} is not byte-identical to HEAD ({how}), so the line "
            "that would govern is uncommitted and authorises nothing. Append the re-sign AND COMMIT it "
            f"before sending. {APPEND_ONLY_IS_ENFORCED}"
        )
    return {"identical_to_head": True, "bytes": len(head_bytes), "record": rel.as_posix()}


def recorded_signoff(path: Path | str | None = None) -> dict[str, Any] | None:
    """The sign-off that GOVERNS: the last line, or None. Earlier lines never authorise."""
    lines = signoff_lines(path)
    if not lines:
        return None
    last = lines[-1]
    return None if "_unreadable" in last else last


#: The sender whose whole import closure is signed.
SENDER_ENTRY = "scripts/astroreg2_send.py"

A_SIGNOFF_IS_FOR_THE_CODE_IT_READ = (
    "the sign-off covers the sender's WHOLE IMPORT CLOSURE, not two files. Two files were never the right "
    "boundary: an edit to enhancer_target.score_element or to alphagenome_adapter after a sign-off would "
    "have passed unnoticed, and those two decide what is BOUGHT and what is RECORDED. The reviewed unit "
    "is everything the sender can reach. The closure is recomputed from the sender by the shared "
    "mf.counting_path, which parses imports and never imports, so it is the same in any checkout of the "
    "same revision; any file ADDED, REMOVED or CHANGED refuses the send and is NAMED. The refusal says a "
    "re-sign is required and never that an authorisation is missing: those are different facts and "
    "conflating them would misreport what Albert agreed to"
)


#: Where the astrocyte answers are written. NOT the sweep's root, and the reason is a trap rather than
#: tidiness: `score_element` writes a per-element file at `cache/<chrom>/<id>.json`, and `load_cached`
#: prefers that file over the chromosome archive. So a new file under the sweep's root would SHADOW THE
#: ARCHIVE and silently change the FROZEN K562 inputs for every future rebuild -- corrupting the very
#: baseline this test is measured against, invisibly and permanently.
ASTROREG2_CACHE_ROOT = Path("data/knowledge/alphagenome/elements_astroreg2")

WHY_A_SEPARATE_CACHE_ROOT = (
    "the sweep's cache root holds the frozen K562 answers the weights were fitted on. Writing astrocyte "
    "answers there would shadow the archive for those element ids and change the frozen inputs of every "
    "later rebuild, which is worse than a wrong number because it would make the baseline wrong too. The "
    "astrocyte answers go to their own root and `crispri`'s astrocyte annotation reads THAT root"
)

#: Why an already-populated request root is a refusal and not a convenience.
AN_EXISTING_ENTRY_MEANS_NOTHING_WOULD_BE_BOUGHT = (
    "`score_element` returns a CACHED answer without a request. The AstroREG-2 list was built from the "
    "elements the sweep had already scored, so against the sweep's root every one of the 1,232 is a cache "
    "hit: the run would have charged 1,232 to the ledger, sent nothing, and produced a result made "
    "entirely of K562 answers labelled as astrocyte work. It would have looked like a success. So the "
    "send REFUSES unless the request cache root holds NONE of the plan's elements, which is the only "
    "state in which a request is actually made"
)


def plan_elements_already_cached(
    plan: list[dict[str, Any]], cache_root: Path | str | None = None, root: Path | None = None
) -> list[str]:
    """The plan's elements that already have an answer in the request cache root, so would not be bought."""
    from genomeos.predict.enhancer_target import load_cached

    base = Path(root) if root is not None else ROOT_FOR_BLOBS
    cache = Path(cache_root) if cache_root is not None else base / ASTROREG2_CACHE_ROOT
    if not Path(cache).is_absolute():
        cache = base / cache
    present = []
    for row in plan:
        if load_cached(row["chrom"], row["element"], Path(cache)) is not None:
            present.append(row["element"])
    return present


def check_nothing_in_the_plan_is_already_cached(
    plan: list[dict[str, Any]], cache_root: Path | str | None = None, root: Path | None = None
) -> dict[str, Any]:
    """Refuse unless EVERY request would really be a request. See the note above."""
    present = plan_elements_already_cached(plan, cache_root, root)
    if present:
        raise SendRefusedError(
            f"{len(present)} of the {len(plan)} planned elements already have a cached answer in the "
            f"request cache root, so they would be charged and NOT bought (first: {present[:3]}). "
            f"{AN_EXISTING_ENTRY_MEANS_NOTHING_WOULD_BE_BOUGHT}"
        )
    return {"plan": len(plan), "already_cached": 0, "cache_root": str(cache_root or ASTROREG2_CACHE_ROOT)}


def sender_closure(entry: str = SENDER_ENTRY, root: Path | None = None) -> dict[str, str]:
    """`path -> blob` for every repository file the sender can reach by import, plus the sender itself."""
    from genomeos import manifest as mf

    base = Path(root) if root is not None else ROOT_FOR_BLOBS
    rels = sorted(set(mf.counting_path(entry, base)) | {entry})
    present = [r for r in rels if (base / r).exists()]
    blobs = git_blobs(present, base)
    missing = [r for r in present if r not in blobs]
    if missing:
        # Silently dropping a file would change the closure whenever `git hash-object` failed under
        # load, so the signed set would depend on machine conditions rather than on the code. A file
        # that exists and cannot be hashed is an error.
        raise SendRefusedError(
            f"the closure could not be computed: {len(missing)} file(s) exist but could not be hashed "
            f"({missing[:3]}). A closure that silently omitted them would make the sign-off depend on "
            "machine load rather than on the code"
        )
    return blobs


def recorded_signoff_words(path: Path | str | None = None) -> str | None:
    rec = recorded_signoff(path)
    return (rec or {}).get("words")


def check_signoff_closure(
    record: dict[str, Any] | None = None,
    entry: str = SENDER_ENTRY,
    root: Path | None = None,
    path: Path | str | None = None,
    run_id: int | None = None,
) -> dict[str, Any]:
    """Refuse unless the sender's whole closure is byte-for-byte what was signed. See the note above."""
    if record is None:
        check_signoff_is_committed(path, root)
    rec = record if record is not None else recorded_signoff(path)
    if rec is None:
        raise SendRefusedError(
            f"no sign-off record is present at {SIGNOFF_RECORD}, so nothing is signed. "
            f"{A_SIGNOFF_IS_FOR_THE_CODE_IT_READ}"
        )
    # A sign-off is for ONE run. A line that names no run, or names another, covers nothing here:
    # reusing run 2's sign-off for run 3 is the same reuse as spending run 1's approval twice.
    if run_id is not None:
        named = rec.get("run_id")
        if named is None:
            raise SendRefusedError(
                f"the governing sign-off line names no run_id, so it does not cover run {run_id}. "
                f"{A_SIGNOFF_IS_FOR_ONE_RUN}"
            )
        if int(named) != int(run_id):
            raise SendRefusedError(
                f"the governing sign-off line covers run {named}, not run {run_id}. "
                f"{A_SIGNOFF_IS_FOR_ONE_RUN}"
            )

    signed = rec.get("signed_closure") or {}
    if not signed:
        raise SendRefusedError(
            f"the sign-off record at {SIGNOFF_RECORD} names no signed_closure, so it signs nothing"
        )
    now = sender_closure(entry, root)
    added = sorted(set(now) - set(signed))
    removed = sorted(set(signed) - set(now))
    changed = sorted(p for p in set(now) & set(signed) if now[p] != signed[p])
    if added or removed or changed:
        parts = []
        if changed:
            parts.append(
                "CHANGED: " + ", ".join(f"{p} (on disk {now[p]}, signed {signed[p]})" for p in changed)
            )
        if added:
            parts.append("ADDED to the closure: " + ", ".join(added))
        if removed:
            parts.append("REMOVED from the closure: " + ", ".join(removed))
        raise SendRefusedError(
            "the supervisor's sign-off was given for different code and a RE-SIGN is required (this is "
            f"not a missing authorisation): {'; '.join(parts)}. {A_SIGNOFF_IS_FOR_THE_CODE_IT_READ}"
        )
    return {"files_signed": len(signed), "signed_at": rec.get("signed_at"), "words": rec.get("words")}


class SendRefusedError(RuntimeError):
    """A clause of Albert's approval is not satisfied, so nothing may be sent."""


def _refuse(clause: str, detail: str) -> None:
    raise SendRefusedError(
        f"Albert's approval is CONDITIONAL on {ASTROREG2_CLAUSES[clause]} and that is not satisfied: {detail}"
    )


def mark_run_complete(ledger: Path | str, **row: Any) -> None:
    """Record that the one authorised run has happened. Makes a second one refusable by name."""
    path = Path(ledger)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(
            json.dumps({**row, "event": RUN_COMPLETE_EVENT, "t": time.strftime("%Y-%m-%dT%H:%M:%S")}) + "\n"
        )


def run_already_completed(ledger: Path | str) -> bool:
    """Whether a run has already finished, read from the ledger by event NAME."""
    path = Path(ledger)
    if not path.exists():
        return False
    for line in path.read_text(errors="replace").splitlines():
        if not line.strip():
            continue
        try:
            row = json.loads(line)
        except ValueError:
            continue
        if isinstance(row, dict) and row.get("event") == RUN_COMPLETE_EVENT:
            return True
    return False


def answer_notes(ledger: Path | str) -> list[dict[str, Any]]:
    """Every `event: answer` note the ledger holds, in order. What the pilot counts.

    The pilot has to count the answers OF THE RUN and not of this process. Counting in memory would
    re-pilot after a resume, or skip the pilot entirely, and since a pilot stop is terminal and consumes
    Albert's one run, being wrong in either direction is expensive.
    """
    path = Path(ledger)
    if not path.exists():
        return []
    out: list[dict[str, Any]] = []
    for line in path.read_text(errors="replace").splitlines():
        if not line.strip():
            continue
        try:
            row = json.loads(line)
        except ValueError:
            continue
        if isinstance(row, dict) and row.get("event") == "answer":
            out.append(row)
    return out


#: Why a crash stop is not the exception to "One run".
A_CRASH_STOP_IS_NOT_A_CONTINUATION = (
    "a partial ledger means a run started and did not finish. Resuming it is a NEW authorisation in "
    "Albert's own words, not a continuation -- the same terms as a pilot stop, which he has been told "
    "consumes his one run. A crash must not quietly become the exception to that. And a resume that "
    "simply restarted the list would RE-PAY for every element already charged while the cap, which "
    "resumes from the ledger, ran out that many requests early -- so it would cost more than he approved "
    "AND score fewer elements than the registration names. Both halves are breaches of what he agreed to"
)


def may_send(
    *,
    activity_result: Path | str,
    registration: Path | str,
    run_id: int,
    signoff: str | None,
    ledger_root: Path | None = None,
    plan: list[dict[str, Any]],
    reviewed_digest: str,
    committed: Any = None,
    resume_authorisation: str | None = None,
) -> dict[str, Any]:
    """Every clause of Albert's approval as its own refusal. Returns only if ALL of them hold.

    Recording the authorisation text does NOT by itself make this return: the approval is conditional,
    and a guard that lifted on a non-empty string would convert a conditional approval into an
    unconditional one, silently and in the direction of spending money.

    `committed` is a predicate taking a path and saying whether it is committed, injected so the check
    can be tested without a repository.
    """
    if not ASTROREG2_AUTHORISATION:
        raise NoAuthorisationError(
            f"no AstroREG-2 approval is recorded. {NO_AUTHORISATION_MEANS_NO_SEND} "
            f"{WHY_THE_RELAYED_TEXT_IS_NOT_THE_APPROVAL}"
        )

    # THIS run's own approval, before any ledger is looked at. The ledger cannot answer this question:
    # a fresh run resolves a fresh file, so "is this ledger finished?" is False for every new run.
    authorisation_for_run(run_id)
    check_earlier_runs_are_complete(run_id, ledger_root)
    totals = check_total_cap(run_id, len(plan), ledger_root)

    # DERIVED, not accepted: see A_GUARD_KEYED_TO_WHAT_THE_CALLER_SUPPLIES_IS_NOT_A_GUARD.
    ledger = ledger_path_for_run(run_id, ledger_root)

    reg = Path(registration)
    if not reg.exists():
        _refuse("scope", f"the registration {reg} the approval names is not present")
    if len(plan) != ASTROREG2_CAP:
        _refuse(
            "scope",
            f"the list holds {len(plan)} requests and the approval names {ASTROREG2_CAP}",
        )
    check_is_the_reviewed_plan(plan, reviewed_digest)

    act = Path(activity_result)
    if not act.exists():
        _refuse("activity", f"the astrocyte activity result {act} is not present")
    if committed is not None and not committed(act):
        _refuse(
            "activity",
            f"{act} is present but NOT COMMITTED. The clause requires the activity to have been "
            "computed, and a result no commit holds is not a computation anyone can check",
        )
    body = json.loads(act.read_text())
    found = (body.get("rule") or {}).get("producer") or {}
    for key, want in AMENDMENT_2_RULE_FINGERPRINT.items():
        got = found.get(key)
        if got != want:
            _refuse(
                "activity",
                f"the committed activity result records {key}={got!r}, not amendment 2's {want!r}. It "
                "was produced under a different rule",
            )

    check_signoff_closure(root=ROOT_FOR_BLOBS, run_id=run_id)
    check_nothing_in_the_plan_is_already_cached(plan)
    if not signoff or "dry run reviewed" not in signoff:
        _refuse(
            "signoff",
            "no supervisor sign-off containing the words 'dry run reviewed' is recorded. It may not be "
            "anticipated, paraphrased or represented by a flag: the actual text, or the refusal stands",
        )

    if run_already_completed(ledger):
        _refuse("one_run", f"the ledger already records a completed run. {ONE_RUN_IS_TERMINAL}")

    charged = ledger_charges(ledger)["charges"]
    if charged and not resume_authorisation:
        _refuse(
            "one_run",
            f"a partial run exists ({charged} charged) and no completion is recorded; a resume needs "
            f"Albert's new word. {A_CRASH_STOP_IS_NOT_A_CONTINUATION}",
        )

    check_adapter_v2(ADAPTER_MODULE)  # read at CALL time so the module under check is substitutable
    # Item (g), the supervisor's: the recording path must keep every row's whole track vector, and
    # the disk it needs must be there before a paid run starts rather than part way through one.
    check_adapter_writes_full_vectors(ADAPTER_MODULE)
    check_disk_for_full_vectors(len(plan), OBSERVED_ROWS_PER_ANSWER, OBSERVED_TRACKS_PER_ANSWER)

    budget = RequestBudget(ledger, cap=ASTROREG2_CAP)
    if budget.charged() != budget.sent:
        _refuse(
            "logged",
            f"the ledger records {budget.charged()} charges and the budget counts {budget.sent}: the "
            "log and the money disagree, so not every request is logged",
        )
    remaining = requests_not_yet_charged(plan, ledger)
    return {
        "may_send": True,
        "resuming": bool(charged),
        "already_charged": charged,
        "remaining": len(remaining),
        "run_id": int(run_id),
        "ledger": str(ledger),
        "totals": totals,
        "cap": ASTROREG2_CAP,
        "requests": len(plan),
        "digest": reviewed_digest,
        "clauses_checked": sorted(ASTROREG2_CLAUSES),
    }


# ------------------------------------------- the sender must prove it is sending the reviewed list

#: Why the sender recomputes the list and compares a digest rather than trusting the reviewed file.
THE_REVIEWED_LIST_IS_THE_SENT_LIST = (
    "a review of a list the sender could quietly regenerate differently is a review of nothing. So the "
    "sender rebuilds the list through plan_requests and the one filter, takes its digest, and REFUSES "
    "if it differs from the digest that was reviewed. The digest is over the list's own canonical "
    "content rather than the file's bytes, so a reformatting does not trip it and a changed, added or "
    "removed request does"
)


def plan_digest(plan: list[dict[str, Any]]) -> str:
    """A canonical digest of the request list: what was reviewed, and what must be sent.

    Taken over each row's identity and the pairs it serves, in the list's order, so the digest moves if
    any request is changed, added, removed or reordered, and does not move for a change in formatting.
    """
    import hashlib

    canon = [
        [
            r["chrom"],
            r["element"],
            r["start"],
            r["end"],
            sorted(r.get("serves_genes", [])),
            sorted(r.get("serves_labels", [])),
        ]
        for r in plan
    ]
    blob = json.dumps(canon, sort_keys=False, separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()


def check_is_the_reviewed_plan(plan: list[dict[str, Any]], reviewed_digest: str) -> str:
    """Refuse unless the recomputed list is byte-for-byte the list that was reviewed."""
    got = plan_digest(plan)
    if got != reviewed_digest:
        raise ValueError(
            f"the recomputed request list has digest {got}, not the reviewed {reviewed_digest}: "
            f"{len(plan)} requests. It may not be sent. {THE_REVIEWED_LIST_IS_THE_SENT_LIST}"
        )
    return got


# ----------------------------------------------------------------- the activity-input gate

#: The two numbers the frozen activity term is made of, and the one place they are read from.
FROZEN_ACTIVITY_COLUMNS = ("DHS.RPM", "H3K27ac.RPM")

#: The units those two numbers are in. The frozen weight on `log_activity` was fitted against these,
#: so a value in any other unit is a different feature however reasonable it looks.
FROZEN_ACTIVITY_UNITS = (
    "reads per million mapped reads in the element, as EPCrisprBenchmark publishes them per pair "
    "(crispri.parse reads DHS.RPM and H3K27ac.RPM; crispri._annotate_one forms "
    "sqrt(dhs * h3k27ac), then log_activity = log1p(activity) and "
    "activity_over_distance = log_activity - log(distance))"
)

#: The registration's own rule for a missing input, quoted from
#: data/results/astroreg_registration.json, terms.feasibility_gate.rule. The gate below applies this
#: and nothing else.
REGISTERED_MISSING_INPUT_RULE = (
    "if an input the frozen model needs is missing, the model is not applied and no substitute is "
    "put in its place: a no-go is recorded instead"
)

#: Why a peak call cannot stand in for the benchmark's RPM, beyond the units differing. Two separate
#: defects, and the second is the one no rescaling can repair.
WHY_A_PEAK_CALL_IS_NOT_THE_INPUT = (
    "a narrowPeak file differs from the benchmark's RPM in two independent ways. (1) UNITS: column "
    "seven is signalValue -- fold change over control for the ChIP-seq peaks, the peak caller's "
    "signal for the DNase peaks -- and a weight fitted on reads per million applied to it is a "
    "different feature. (2) CENSORING: a peak call carries no value outside a called peak, so an "
    "element in no peak can only be given zero, while the benchmark's RPM is defined for every "
    "element and descends to zero continuously. Rescaling can address the first defect and cannot "
    "address the second, because the values that are missing were never measured into the file"
)


def activity_verdict(available: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """Whether the frozen activity term can be supplied for the astrocyte pairs, and the reading.

    `available` maps each of `FROZEN_ACTIVITY_COLUMNS` to what was found for it, each entry naming
    `source` (where a value would come from), `units` and `in_frozen_units` (whether that source is
    in `FROZEN_ACTIVITY_UNITS`). Every column must be present in the mapping: a column nobody
    looked for is not an absent input, it is an unfinished check, and it is refused as one.

    The verdict is the registration's rule applied: go only when every column is available in its
    own units, and otherwise a no-go that names each column that is not and quotes the rule.
    """
    missing_from_check = [c for c in FROZEN_ACTIVITY_COLUMNS if c not in available]
    if missing_from_check:
        raise ValueError(
            "the activity gate was asked for a verdict without a finding for "
            f"{missing_from_check}: a column nobody looked for is an unfinished check, not an "
            "absent input"
        )
    not_in_units = [c for c in FROZEN_ACTIVITY_COLUMNS if not available[c]["in_frozen_units"]]
    go = not not_in_units
    return {
        "frozen_activity_columns": list(FROZEN_ACTIVITY_COLUMNS),
        "frozen_activity_units": FROZEN_ACTIVITY_UNITS,
        "found": {c: available[c] for c in FROZEN_ACTIVITY_COLUMNS},
        "columns_not_in_frozen_units": not_in_units,
        "go": go,
        "registered_rule": REGISTERED_MISSING_INPUT_RULE,
        "reading": (
            "every input the frozen activity term needs is available in its own units, so the "
            "model may be applied"
            if go
            else "an input the frozen model needs is missing in its own units for "
            + ", ".join(not_in_units)
            + ", so by the registration's own rule the model is not applied, no substitute is put "
            "in its place, and this no-go is recorded instead"
        ),
        "why_a_peak_call_is_not_the_input": WHY_A_PEAK_CALL_IS_NOT_THE_INPUT,
    }


def requests_permitted(verdict: dict[str, Any]) -> int:
    """How many AlphaGenome requests the gate permits: the authorised number, or none.

    The authorisation is for the registered test. A no-go means the registered test cannot be
    computed from what the money would buy, so the number it permits is zero -- not 1,322 spent on a
    deletion cache no registered claim can be made from.
    """
    return AUTHORISED_REQUESTS if verdict["go"] else 0


# ---------------------------------- item (g): the full track vector of every paid answer is kept

#: Sign-off item (g). NOT one of Albert's clauses: the SUPERVISOR's requirement, labelled like item
#: (f) so a refusal never attributes to him a condition he did not state.
ITEM_G_IS_A_SUPERVISOR_REQUIREMENT = (
    "item (g) of the supervisor's sign-off checklist, not a clause of Albert's approval. His approval "
    "names the activity precondition, the sign-off, one run, every request logged and the scope; this "
    "is an additional condition the supervisor imposed, and conflating the two would misreport what he "
    "agreed to"
)

WHY_THE_FULL_TRACK_VECTOR_IS_KEPT = (
    "the recording path keeps each row's whole value multiset for CELL_TRACKS -- K562, HepG2, GM12878 "
    "and IMR-90 -- plus extremes and a mean, and discards every other column. For 1,232 ASTROCYTE "
    "screen elements that means every BRAIN-TISSUE track's value would be BOUGHT AND THROWN AWAY to "
    "keep four non-neural cell lines. Paid data is the one kind this project cannot re-fetch for free, "
    "so the full vector is persisted for every gene row of every paid answer. This is the cache-hit "
    "defect one layer out: there the run would have charged 1,232 and bought nothing, here it would "
    "charge 1,232 and discard the part that answers the question it was authorised for. Both were "
    "invisible to every stub, because a stub returns what the code asks for -- and the code was asking "
    "for four cell lines"
)

WIDENING_IS_FREE_ONLY_BEFORE_THE_RESPONSE_ARRIVES = (
    "keeping more of a response costs nothing while the response has not been received, and a fresh "
    "request once it has: a cached answer no longer holds what it dropped. The 1,232 responses do not "
    "exist yet, which is why item (g) lands BEFORE the send rather than after it"
)

#: float64, not float32. A lossy store of data that cannot be re-fetched is the same mistake one layer
#: down: the values would come back close to what was paid for rather than equal to it.
FULL_VECTOR_DTYPE = "<f8"

#: The key the full vector is written under, on each gene-axis row, beside the existing columns.
FULL_VECTOR_KEY = "track_vector"

#: The free disk this project will not go below, matching scripts/capacity_gate.py's own floor.
DISK_FLOOR_BYTES = 10 * 1024**3


class VectorRefusedError(SendRefusedError):
    """A vector is absent, short, long or unreadable, so the answer does not record what was paid for."""


def pack_track_vector(values: list[float]) -> str:
    """One gene row's whole track vector as base64 float64, in the response's own column order."""
    import array
    import base64

    buf = array.array("d", [float(v) for v in values])
    if sys.byteorder != "little":  # the dtype is declared little-endian, so say so on either machine
        buf.byteswap()
    return base64.b64encode(buf.tobytes()).decode("ascii")


def unpack_track_vector(blob: str, expected: int | None = None) -> list[float]:
    """The vector back out, REFUSING a truncated or over-long one rather than returning what fits."""
    import array
    import base64

    raw = base64.b64decode(blob.encode("ascii"), validate=True)
    if len(raw) % 8:
        raise VectorRefusedError(
            f"a track vector of {len(raw)} bytes is not a whole number of float64 values, so it is "
            "truncated and the answer does not record what was paid for"
        )
    buf = array.array("d")
    buf.frombytes(raw)
    if sys.byteorder != "little":
        buf.byteswap()
    out = list(buf)
    if expected is not None and len(out) != int(expected):
        raise VectorRefusedError(
            f"a track vector carries {len(out)} values and the answer's own track count is {expected}: "
            "a short vector is a partly discarded answer and a long one is not this response's axis. "
            f"{WHY_THE_FULL_TRACK_VECTOR_IS_KEPT}"
        )
    return out


def full_vector_entry(values: list[float]) -> dict[str, Any]:
    """The compact record written on a gene-axis row: the count, the dtype and the packed bytes."""
    return {"n": len(values), "dtype": FULL_VECTOR_DTYPE, "b64": pack_track_vector(values)}


def with_full_track_vectors(
    axis: dict[str, Any], matrix: Any, track_names: list[str], tracks_sha256: str
) -> dict[str, Any]:
    """`axis` with every row's whole track vector added, and the axis's track names recorded once.

    ADDITIVE and nothing else: no existing key is written, removed or reordered, so the registered
    frozen features -- which read the effects and the per-cell values -- cannot move by one bit. The
    names are kept once per answer rather than per row, bound to the `tracks_sha256` the adapter
    computes from the response's own var axis, so 371 NAMED values come back out of one row.
    """
    out = dict(axis)
    rows = []
    for row in axis.get("rows", []):
        gi = int(row.get("row", len(rows)))
        values = [float(matrix[gi, ti]) for ti in range(len(track_names))]
        rows.append({**row, FULL_VECTOR_KEY: full_vector_entry(values)})
    out["rows"] = rows
    out["track_names"] = list(track_names)
    out["tracks_sha256"] = tracks_sha256
    out["full_vector_note"] = WHY_THE_FULL_TRACK_VECTOR_IS_KEPT
    return out


def named_track_values(axis: dict[str, Any], row_index: int) -> dict[str, float]:
    """One gene row's vector as {track name: value}, refusing unless the lengths agree."""
    names = axis.get("track_names") or []
    if not names:
        raise VectorRefusedError(
            "the answer records no track names, so its values cannot be named and a reader cannot tell "
            f"a brain track from a cell line. {WHY_THE_FULL_TRACK_VECTOR_IS_KEPT}"
        )
    for row in axis.get("rows", []):
        if int(row.get("row", -1)) != int(row_index):
            continue
        entry = row.get(FULL_VECTOR_KEY)
        if not entry:
            raise VectorRefusedError(
                f"gene row {row_index} carries no {FULL_VECTOR_KEY}, so this paid row's non-cell-line "
                f"tracks are gone. {WHY_THE_FULL_TRACK_VECTOR_IS_KEPT}"
            )
        values = unpack_track_vector(entry["b64"], expected=len(names))
        return dict(zip(names, values, strict=True))
    raise VectorRefusedError(f"the answer has no gene row {row_index}")


def check_full_vectors_in_answer(answer: dict[str, Any]) -> dict[str, Any]:
    """Refuse unless EVERY gene row of EVERY output of this answer carries its whole track vector."""
    model = answer.get("model") or {}
    axes = model.get("gene_axis_outputs")
    if not axes:
        raise VectorRefusedError(
            "the answer records no gene axis at all, so there is nothing to check and nothing was "
            f"kept. {ITEM_G_IS_A_SUPERVISOR_REQUIREMENT}. {WHY_THE_FULL_TRACK_VECTOR_IS_KEPT}"
        )
    rows = 0
    for axis in axes:
        names = axis.get("track_names") or []
        total = axis.get("tracks_total")
        if not names:
            raise VectorRefusedError(
                f"output {axis.get('output')} records no track names. {ITEM_G_IS_A_SUPERVISOR_REQUIREMENT}"
            )
        if total is not None and len(names) != int(total):
            raise VectorRefusedError(
                f"output {axis.get('output')} names {len(names)} tracks and counts {total}: the names "
                "and the axis disagree, so a named value cannot be trusted to be that track's value"
            )
        if not axis.get("tracks_sha256"):
            raise VectorRefusedError(
                f"output {axis.get('output')} records no tracks_sha256, so the name list it was read "
                "on cannot be identified later"
            )
        for row in axis.get("rows", []):
            entry = row.get(FULL_VECTOR_KEY)
            if not entry:
                raise VectorRefusedError(
                    f"output {axis.get('output')} row {row.get('row')} carries no {FULL_VECTOR_KEY}. "
                    f"{WHY_THE_FULL_TRACK_VECTOR_IS_KEPT}"
                )
            unpack_track_vector(entry["b64"], expected=len(names))
            rows += 1
    return {"outputs": len(axes), "gene_rows_with_full_vectors": rows}


def full_vector_disk_estimate(
    elements: int, rows_per_element: float, tracks: int, free_bytes: int | None = None
) -> dict[str, Any]:
    """What keeping the full vectors costs on disk, measured from the format rather than assumed.

    base64 of float64 is 4 characters per 3 bytes, so a row of `tracks` values is 8*tracks bytes
    packed and ceil(8*tracks/3)*4 characters written, plus the row's existing columns.
    """
    import math

    packed = 8 * int(tracks)
    encoded = math.ceil(packed / 3) * 4
    per_row = encoded + 40  # the {"n": .., "dtype": "<f8", "b64": ".."} wrapper, measured below in tests
    per_answer = per_row * rows_per_element
    total = int(per_answer * elements)
    free = shutil.disk_usage(str(ROOT_FOR_BLOBS)).free if free_bytes is None else int(free_bytes)
    return {
        "elements": int(elements),
        "rows_per_element": rows_per_element,
        "tracks": int(tracks),
        "bytes_per_row": per_row,
        "bytes_per_answer": int(per_answer),
        "bytes_total": total,
        "mb_total": round(total / 1024**2, 1),
        "free_bytes": free,
        "free_gb": round(free / 1024**3, 2),
        "floor_bytes": DISK_FLOOR_BYTES,
        "fits_above_the_floor": free - total >= DISK_FLOOR_BYTES,
    }


def check_disk_for_full_vectors(
    elements: int, rows_per_element: float, tracks: int, free_bytes: int | None = None
) -> dict[str, Any]:
    """Refuse to start a paid run that would leave the machine below the floor it already paged under."""
    est = full_vector_disk_estimate(elements, rows_per_element, tracks, free_bytes)
    if not est["fits_above_the_floor"]:
        raise VectorRefusedError(
            f"keeping the full track vectors for {elements} answers needs about {est['mb_total']} MB "
            f"and {est['free_gb']} GB is free, which would leave less than the "
            f"{DISK_FLOOR_BYTES // 1024**3} GB floor. The machine paged twice tonight, so this is a "
            f"refusal and not a warning. {ITEM_G_IS_A_SUPERVISOR_REQUIREMENT}"
        )
    return est


#: What one answer of the sweep carried, read from adapter v2's own note on the response it recorded:
#: 29 gene-axis rows on a 371-track axis. Cited rather than assumed, and a test reads a cached sweep
#: answer when one is on the machine so a drift in either number fails instead of passing quietly.
OBSERVED_ROWS_PER_ANSWER = 29
OBSERVED_TRACKS_PER_ANSWER = 371
OBSERVED_SHAPE_SOURCE = (
    "adapter v2's own note on the response it recorded: 29 gene-axis rows x 371 tracks, ~12.6 KiB at "
    "the SUMMARY level. The full vector is materially larger, which is why the estimate is computed "
    "from the format and checked against the floor rather than assumed to be small"
)


def check_adapter_writes_full_vectors(module_name: str = ADAPTER_MODULE) -> dict[str, Any]:
    """Refuse unless the adapter the runner writes through keeps each row's WHOLE track vector.

    Item (g). Checked on the module's SOURCE rather than by running it, because the thing that must be
    true is that the recording path persists the vector -- and a stub response cannot show that: a stub
    returns what the code asks for, and the code was asking for four cell lines.
    """
    import ast
    import importlib.util

    spec = importlib.util.find_spec(module_name)
    origin = getattr(spec, "origin", None) if spec is not None else None
    if not origin or not Path(origin).exists():
        raise VectorRefusedError(
            f"item (g): the adapter module {module_name} has no source file to read, so the runner "
            f"cannot be shown to keep the full track vector. {WHY_THE_FULL_TRACK_VECTOR_IS_KEPT}"
        )
    tree = ast.parse(Path(origin).read_text())
    recorder = next(
        (n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == "recorded_axis"),
        None,
    )
    if recorder is None:
        raise VectorRefusedError(
            f"item (g): {module_name} has no recorded_axis, so the gene-axis recording path is not "
            f"where this check believes it is. {WHY_THE_FULL_TRACK_VECTOR_IS_KEPT}"
        )
    body = ast.dump(recorder)
    keeps = FULL_VECTOR_KEY in body or "with_full_track_vectors" in body
    if not keeps:
        raise VectorRefusedError(
            f"item (g): {module_name}.recorded_axis keeps CELL_TRACKS only and writes no "
            f"{FULL_VECTOR_KEY}, so for 1,232 astrocyte elements every brain-tissue track would be "
            f"bought and discarded. {WHY_THE_FULL_TRACK_VECTOR_IS_KEPT}. "
            f"{WIDENING_IS_FREE_ONLY_BEFORE_THE_RESPONSE_ARRIVES}"
        )
    return {"module": module_name, "recorded_axis_keeps_the_full_vector": True}
