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


def astroreg2_budget(ledger: Path | str) -> RequestBudget:
    """A budget for AstroREG-2's 1,232, obtainable only once its own approval is recorded.

    The authorisation is checked BEFORE the budget exists, so there is no object to send with until the
    approval is on record. The cap is `ASTROREG2_CAP` and never `AUTHORISED_REQUESTS`.
    """
    astroreg2_authorisation()
    return RequestBudget(ledger, cap=ASTROREG2_CAP)


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


#: The sign-off record. A FILE READ AT RUNTIME, never imported, and deliberately not in code.
#:
#: The first version of this stored the signed blob of `genomeos/attribution/astrorun.py` INSIDE
#: `genomeos/attribution/astrorun.py`. That is a fixed point with no solution: writing a blob value into
#: the file changes the file, which changes its blob, so the recorded value can never equal the computed
#: one. Demonstrated rather than reasoned about -- writing the current blob in produced a third,
#: different blob. The intended property was "a change invalidates the sign-off"; the property built was
#: "the sign-off can never be valid", and no re-sign by anyone could have passed.
#:
#: Keeping the record OUTSIDE the code dissolves it rather than working around it: a record that
#: describes the code from outside can name the code's hashes without being part of what it hashes. It
#: sits beside the ledger because both are records ABOUT a run rather than parts of one.
SIGNOFF_RECORD = Path("data/ledgers/astroreg2_signoff.json")

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


def sender_closure(entry: str = SENDER_ENTRY, root: Path | None = None) -> dict[str, str]:
    """`path -> blob` for every repository file the sender can reach by import, plus the sender itself."""
    from genomeos import manifest as mf

    base = Path(root) if root is not None else ROOT_FOR_BLOBS
    out: dict[str, str] = {}
    for rel in sorted(set(mf.counting_path(entry, base)) | {entry}):
        blob = git_blob(base / rel, base)
        if blob:
            out[rel] = blob
    return out


def recorded_signoff(path: Path | str | None = None) -> dict[str, Any] | None:
    """The sign-off record as written by whoever holds the supervisor's words first-hand, or None."""
    p = Path(path) if path is not None else ROOT_FOR_BLOBS / SIGNOFF_RECORD
    if not p.exists():
        return None
    try:
        body = json.loads(p.read_text())
    except ValueError:
        return None
    return body if isinstance(body, dict) else None


def recorded_signoff_words(path: Path | str | None = None) -> str | None:
    rec = recorded_signoff(path)
    return (rec or {}).get("words")


def check_signoff_closure(
    record: dict[str, Any] | None = None,
    entry: str = SENDER_ENTRY,
    root: Path | None = None,
    path: Path | str | None = None,
) -> dict[str, Any]:
    """Refuse unless the sender's whole closure is byte-for-byte what was signed. See the note above."""
    rec = record if record is not None else recorded_signoff(path)
    if rec is None:
        raise SendRefusedError(
            f"no sign-off record is present at {SIGNOFF_RECORD}, so nothing is signed. "
            f"{A_SIGNOFF_IS_FOR_THE_CODE_IT_READ}"
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
    ledger: Path | str,
    signoff: str | None,
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

    check_signoff_closure(root=ROOT_FOR_BLOBS)
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
