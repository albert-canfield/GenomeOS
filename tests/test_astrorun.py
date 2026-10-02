"""The money guard and the activity-input gate of the authorised AstroREG run.

The cap tests are written before any request can be sent. They do not check that a warning is
printed at 1,322; they check that the 1,323rd request cannot be obtained at all, that the count
stops at the cap, and that the ledger on disk holds exactly one line per charge.
"""

from __future__ import annotations

import csv
import json
import math
import sys
import threading
from pathlib import Path

import pytest

from genomeos.attribution import astrorun


def grant_run(monkeypatch, run=1, requests=None, consumed=False, words="<the run's words>"):
    """Record an approval for ONE run, so a test of another clause is not stopped by the approval one.

    Run 1's REAL approval is marked consumed, because run 1 happened under it. A test that needs to
    reach a later clause grants its own run an approval here rather than reusing a spent one.
    """
    monkeypatch.setattr(
        astrorun,
        "ASTROREG2_AUTHORISATIONS",
        {
            run: {
                "words": words,
                "requests": astrorun.ASTROREG2_CAP if requests is None else requests,
                "consumed": consumed,
            }
        },
    )


#: `needs_local_data`, the supervisor's AMENDED acceptance rule of 2026-10-02, which I had not read
#: when I converted these two skips into failures. The amended rule is explicit: the verdict is a
#: worktree of the committed tree WITH the git-ignored stores linked read-only, so a test needing local
#: data RUNS there; where the stores are genuinely absent -- CI, or a bare worktree -- it SKIPS BY NAME
#: and never fails. My failures would have RED CI, which is the same break that 21 tests caused in a
#: worktree, 16 of them in tests/test_context_evidence.py on this very file.
#:
#: The marker does not weaken either assertion: both still claim every one of the response's tracks
#: comes back named and equal, and both still read the figure from the file rather than from a
#: constant. It changes WHERE they run, not what they claim -- and the skip names the missing store, so
#: a reader can tell "not run here" from "passed".
needs_track_metadata = pytest.mark.needs_local_data(
    "data/cache/entex/alphagenome_track_metadata_copy.csv",
    how="scripts/entex_feasibility.py writes it (AG_METADATA); data/cache is machine-local by the "
    "data boundary, so a fresh checkout cannot have it",
)


class TestTheRunIsAboutTheTreeItNames:
    """A test, not a belief: the code under test must come from the tree these tests came from.

    THE INCIDENT. An acceptance run in a worktree, using the main checkout's venv, reported item (g)'s
    no-skip test PASSED in a worktree that has no data/cache at all -- because the code it imported was
    the MAIN checkout's, while the status file named the worktree's tree. Nothing in the status file
    could show that.

    THE MECHANISM, measured rather than reasoned about. `tests/` has no __init__.py, so pytest's
    prepend import mode makes `<tree>/tests` sys.path[0] and the TREE ROOT is never on sys.path at all.
    `import genomeos` therefore falls through to site-packages, where the main venv's editable install
    `_editable_impl_genomeos.pth` names the main checkout, and the main checkout wins. Running the same
    interpreter as `python -c` resolves to the worktree instead, because that form puts the working
    directory on sys.path -- which is why the hazard is invisible to a `python -c` probe. Measured, one
    variable apart:

        pytest, no PYTHONPATH : sys.path[0]=<worktree>/tests, root on sys.path False -> MAIN's genomeos
        pytest, PYTHONPATH=wt : sys.path[0]=<worktree>/tests, root on sys.path True  -> WORKTREE's
        python -c             : sys.path[0]='',                                       -> WORKTREE's

    PYTHONPATH=<tree> fixes it, but this assertion is what makes a wrong-tree run impossible to
    publish, because it does not depend on anyone remembering the env var or on which theory of the
    cause is right.
    """

    def test_the_imported_genomeos_comes_from_THIS_tree(self):
        import genomeos

        code = Path(genomeos.__file__).resolve().parents[1]
        tests = Path(__file__).resolve().parents[1]
        print(f"\ngenomeos.__file__ : {genomeos.__file__}")
        print(f"tree under test   : {tests}")
        assert code == tests, (
            f"these tests are {tests} and the genomeos they import is {code}. The run would judge one "
            "tree and report another, which is how an acceptance passes for a tree it never read. If "
            "this is a worktree, export PYTHONPATH=<worktree>; see this class's docstring for the "
            "measured mechanism"
        )

    def test_the_tree_under_test_is_the_one_the_verdict_would_name(self):
        """The same question asked of the verdict's own notion of the repository."""
        import genomeos
        from genomeos.attribution import astrorun

        assert Path(astrorun.ROOT_FOR_BLOBS).resolve() == Path(__file__).resolve().parents[1], (
            "ROOT_FOR_BLOBS is where this module resolves paths from, so if it is not this tree then "
            "every committed-file check in the send path is about another checkout"
        )
        assert Path(genomeos.__file__).resolve().parents[1] == Path(astrorun.ROOT_FOR_BLOBS).resolve()


class TestTheAuthorisedNumber:
    def test_the_cap_is_the_number_albert_approved(self):
        assert astrorun.AUTHORISED_REQUESTS == 1322

    def test_the_authorisation_is_quoted_with_its_number_in_it(self):
        assert "1,322 AlphaGenome requests" in astrorun.AUTHORISATION


class TestTheCapRefuses:
    def test_the_1323rd_request_is_refused_and_the_count_stops_at_1322(self, tmp_path):
        b = astrorun.RequestBudget(tmp_path / "ledger.jsonl")
        for _ in range(astrorun.AUTHORISED_REQUESTS):
            b.take(id="e", chrom="chr1", start=1, end=2)
        assert b.sent == 1322
        assert b.remaining() == 0
        with pytest.raises(astrorun.CapRefusedError) as exc:
            b.take(id="one too many", chrom="chr1", start=1, end=2)
        assert b.sent == 1322, "a refused request must not be charged"
        assert "1322 of 1322 sent" in str(exc.value)
        assert "number 1323" in str(exc.value)

    def test_the_refusal_names_the_count_reached_and_the_cap(self, tmp_path):
        b = astrorun.RequestBudget(tmp_path / "l.jsonl", cap=3)
        for _ in range(3):
            b.take(id="e")
        with pytest.raises(astrorun.CapRefusedError) as exc:
            b.take(id="e")
        assert "3 of 3 sent" in str(exc.value)
        assert "would be number 4" in str(exc.value)

    def test_the_cap_is_read_from_the_constant_not_hard_coded(self, tmp_path):
        """A smaller cap refuses at its own number, so the check is on `cap` and not on 1322."""
        b = astrorun.RequestBudget(tmp_path / "l.jsonl", cap=1)
        b.take(id="e")
        with pytest.raises(astrorun.CapRefusedError):
            b.take(id="e")

    def test_a_cap_of_zero_refuses_the_very_first_request(self, tmp_path):
        b = astrorun.RequestBudget(tmp_path / "l.jsonl", cap=0)
        with pytest.raises(astrorun.CapRefusedError):
            b.take(id="e")
        assert b.sent == 0
        assert not (tmp_path / "l.jsonl").exists(), "a refusal writes no charge"

    def test_a_negative_cap_is_refused_at_construction(self, tmp_path):
        with pytest.raises(ValueError, match="cannot be negative"):
            astrorun.RequestBudget(tmp_path / "l.jsonl", cap=-1)

    def test_requests_already_sent_count_against_the_cap(self, tmp_path):
        """A budget told 1,320 are already spent permits two more, not another 1,322."""
        b = astrorun.RequestBudget(tmp_path / "l.jsonl", sent=1320)
        b.take(id="a")
        b.take(id="b")
        with pytest.raises(astrorun.CapRefusedError):
            b.take(id="c")

    def test_concurrent_workers_cannot_both_take_the_last_request(self, tmp_path):
        b = astrorun.RequestBudget(tmp_path / "l.jsonl", cap=50)
        taken, refused = [], []

        def worker():
            for _ in range(20):
                try:
                    taken.append(b.take(id="e"))
                except astrorun.CapRefusedError:
                    refused.append(1)

        threads = [threading.Thread(target=worker) for _ in range(5)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        assert len(taken) == 50
        assert sorted(taken) == list(range(1, 51)), "every charge has its own running total"
        assert len(refused) == 50
        assert b.charged() == 50


class TestEveryRequestIsLogged:
    def test_each_charge_writes_exactly_one_ledger_line_before_it_returns(self, tmp_path):
        p = tmp_path / "l.jsonl"
        b = astrorun.RequestBudget(p, cap=5)
        for i in range(5):
            n = b.take(id=f"el{i}", chrom="chr21", start=i, end=i + 10)
            assert p.exists()
            assert len(p.read_text().splitlines()) == n, "the line is written before permission"
        rows = [json.loads(line) for line in p.read_text().splitlines()]
        assert [r["n"] for r in rows] == [1, 2, 3, 4, 5]
        assert all(r["event"] == "request" for r in rows)
        assert all(r["cap"] == 5 for r in rows)
        assert all(r["t"] and r["id"] and r["chrom"] for r in rows)

    def test_the_charge_carries_what_was_asked_for(self, tmp_path):
        b = astrorun.RequestBudget(tmp_path / "l.jsonl", cap=1)
        b.take(id="EH38E1234", chrom="chr11", start=5_000, end=5_300, gene="HSPB1")
        row = json.loads((tmp_path / "l.jsonl").read_text().splitlines()[0])
        assert row["id"] == "EH38E1234"
        assert (row["chrom"], row["start"], row["end"]) == ("chr11", 5_000, 5_300)
        assert row["gene"] == "HSPB1"

    def test_notes_are_not_charges_and_do_not_move_the_count(self, tmp_path):
        p = tmp_path / "l.jsonl"
        b = astrorun.RequestBudget(p, cap=2)
        b.take(id="a")
        b.note(event="answer", id="a", status="ok")
        b.note(event="cached", id="b")
        b.note(event="quota", id="c", message="RESOURCE_EXHAUSTED")
        assert b.sent == 1
        assert b.charged() == 1, "only `event: request` lines are requests paid for"
        assert len(p.read_text().splitlines()) == 4
        assert json.loads(p.read_text().splitlines()[1])["sent_so_far"] == 1

    def test_charged_reads_the_committed_ledger_not_the_in_memory_count(self, tmp_path):
        p = tmp_path / "l.jsonl"
        astrorun.RequestBudget(p, cap=3).take(id="a")
        assert astrorun.RequestBudget(p, cap=3).charged() == 1


class TestARestartCannotDoubleSpend:
    """The requirement the death of the registering session added.

    A budget that counted from zero on startup would bound a process and not the money: kill a run
    that had spent 1,322 and the next one could spend 1,322 again. So these tests write a ledger
    holding charges, start a FRESH budget against it, and check the remaining cap and the refusal.
    """

    def ledger_with(self, path, charges, notes=0):
        with path.open("a") as fh:
            for i in range(charges):
                fh.write(json.dumps({"event": "request", "n": i + 1, "id": f"e{i}"}) + "\n")
            for i in range(notes):
                fh.write(json.dumps({"event": "answer", "id": f"e{i}", "status": "ok"}) + "\n")
        return path

    def test_a_fresh_budget_over_a_ledger_of_n_charges_has_1322_minus_n_left(self, tmp_path):
        """The brief's own check: remaining is 1322 - N, and request 1322 - N + 1 is refused."""
        n = 900
        p = self.ledger_with(tmp_path / "l.jsonl", n)
        b = astrorun.RequestBudget(p)
        assert b.resumed_from_ledger == n
        assert b.sent == n
        assert b.remaining() == astrorun.AUTHORISED_REQUESTS - n == 422
        for _ in range(astrorun.AUTHORISED_REQUESTS - n):
            b.take(id="e", chrom="chr1", start=1, end=2)
        assert b.sent == astrorun.AUTHORISED_REQUESTS
        with pytest.raises(astrorun.CapRefusedError) as exc:
            b.take(id="one too many")
        assert "1322 of 1322 sent" in str(exc.value)
        assert b.charged() == astrorun.AUTHORISED_REQUESTS, "the ledger holds 1322, not 1322 + 900"

    def test_a_ledger_already_holding_the_cap_refuses_the_very_first_request(self, tmp_path):
        p = self.ledger_with(tmp_path / "l.jsonl", astrorun.AUTHORISED_REQUESTS)
        b = astrorun.RequestBudget(p)
        assert b.remaining() == 0
        with pytest.raises(astrorun.CapRefusedError):
            b.take(id="the second full spend")
        assert b.charged() == astrorun.AUTHORISED_REQUESTS

    def test_a_ledger_holding_more_than_the_cap_still_refuses(self, tmp_path):
        p = self.ledger_with(tmp_path / "l.jsonl", 7)
        b = astrorun.RequestBudget(p, cap=5)
        assert b.sent == 7
        with pytest.raises(astrorun.CapRefusedError):
            b.take(id="e")

    def test_two_restarts_in_a_row_cannot_between_them_exceed_the_cap(self, tmp_path):
        p = tmp_path / "l.jsonl"
        first = astrorun.RequestBudget(p, cap=10)
        for _ in range(6):
            first.take(id="first")
        second = astrorun.RequestBudget(p, cap=10)
        assert second.remaining() == 4
        for _ in range(4):
            second.take(id="second")
        with pytest.raises(astrorun.CapRefusedError):
            second.take(id="e")
        third = astrorun.RequestBudget(p, cap=10)
        with pytest.raises(astrorun.CapRefusedError):
            third.take(id="e")
        assert astrorun.ledger_charges(p)["request_lines"] == 10

    def test_notes_in_the_ledger_are_not_resumed_as_charges(self, tmp_path):
        p = self.ledger_with(tmp_path / "l.jsonl", 4, notes=9)
        b = astrorun.RequestBudget(p, cap=10)
        assert b.resumed_from_ledger == 4
        assert b.remaining() == 6

    def test_a_line_truncated_by_a_kill_is_counted_as_a_charge(self, tmp_path):
        """The line is written before the request is permitted, so a half-written line was a charge."""
        p = self.ledger_with(tmp_path / "l.jsonl", 3)
        with p.open("a") as fh:
            fh.write('{"event": "request", "n": 4, "id": "e3"')  # killed mid-write, no newline
        counts = astrorun.ledger_charges(p)
        assert counts == {
            "charges": 4,
            "request_lines": 3,
            "unreadable_lines": 1,
            "other_event_lines": 0,
        }
        assert astrorun.RequestBudget(p, cap=10).remaining() == 6

    def test_the_ledger_count_and_an_explicit_sent_are_not_added_together(self, tmp_path):
        """A caller that passes what it read from the ledger must not charge it twice."""
        p = self.ledger_with(tmp_path / "l.jsonl", 5)
        assert astrorun.RequestBudget(p, cap=10, sent=5).remaining() == 5
        assert astrorun.RequestBudget(p, cap=10, sent=8).remaining() == 2, "the larger is taken"

    def test_a_restart_appends_to_the_ledger_rather_than_truncating_it(self, tmp_path):
        p = self.ledger_with(tmp_path / "l.jsonl", 2)
        astrorun.RequestBudget(p, cap=10).take(id="after the restart")
        rows = [json.loads(line) for line in p.read_text().splitlines()]
        assert [r["event"] for r in rows] == ["request", "request", "resume", "request"]
        assert rows[2]["charges_found_in_the_ledger"] == 2
        assert rows[2]["remaining"] == 8
        assert rows[3]["n"] == 3, "the running total continues from the ledger, not from 1"

    def test_a_first_run_with_no_ledger_writes_no_resume_line(self, tmp_path):
        p = tmp_path / "l.jsonl"
        b = astrorun.RequestBudget(p, cap=3)
        assert b.resumed_from_ledger == 0
        assert not p.exists(), "nothing is written until something is charged"
        b.take(id="e")
        assert [json.loads(line)["event"] for line in p.read_text().splitlines()] == ["request"]

    def test_the_reason_a_restart_resumes_is_recorded(self):
        assert "cap must bound the money, not the process" in astrorun.A_RESTART_MAY_NOT_DOUBLE_SPEND

    def test_the_guard_says_what_it_does_not_protect(self):
        """The ledger bounds a restart and not two concurrent processes, and says so."""
        limit = astrorun.WHAT_THE_LEDGER_CANNOT_DO
        assert "does not make them atomic ACROSS two" in limit
        assert "one executor holding the" in limit, "the operating rule is named, not implied"

    def test_two_processes_started_at_once_are_not_protected_and_the_code_says_so(self, tmp_path):
        """By demonstration: interleave the read and the append, and the same element is sent twice.

        This is the limitation above, shown rather than asserted. Both budgets are constructed BEFORE
        either charges, which is what two processes starting together do, and the result is two
        charges for one element against a cap of 1.
        """
        p = tmp_path / "l.jsonl"
        plan = [{"chrom": "chr1", "element": "E0", "start": 1, "end": 2}]
        a = astrorun.RequestBudget(p, cap=1)
        b = astrorun.RequestBudget(p, cap=1)
        assert astrorun.requests_not_yet_charged(plan, p) == plan
        a.take(chrom="chr1", element="E0")
        b.take(chrom="chr1", element="E0")
        charged = [
            json.loads(line)["element"]
            for line in p.read_text().splitlines()
            if json.loads(line).get("event") == "request"
        ]
        assert charged == ["E0", "E0"], (
            "one element, two charges, cap of 1: the ledger cannot close the concurrent case, which "
            "is why one executor holds the budget"
        )


class TestAResumeDoesNotRePayForALedgeredRequest:
    """The count is not enough: a resume must also not buy the same ELEMENT twice."""

    def plan(self, n):
        return [
            {"chrom": "chr1", "element": f"E{i}", "start": i * 100, "end": i * 100 + 50} for i in range(n)
        ]

    def ledger_charging(self, path, plan):
        b = astrorun.RequestBudget(path, cap=len(plan) + 50)
        for r in plan:
            b.take(chrom=r["chrom"], element=r["element"], start=r["start"], end=r["end"])
        return b

    def test_the_remainder_excludes_every_element_the_ledger_already_holds(self, tmp_path):
        p = tmp_path / "l.jsonl"
        plan = self.plan(10)
        self.ledger_charging(p, plan[:4])
        left = astrorun.requests_not_yet_charged(plan, p)
        assert [r["element"] for r in left] == ["E4", "E5", "E6", "E7", "E8", "E9"]
        assert astrorun.charged_elements(p) == {("chr1", f"E{i}") for i in range(4)}

    def test_a_resumed_run_sends_each_element_exactly_once_across_two_processes(self, tmp_path):
        p = tmp_path / "l.jsonl"
        plan = self.plan(10)
        first = astrorun.RequestBudget(p, cap=10)
        for r in astrorun.requests_not_yet_charged(plan, p)[:6]:
            first.take(chrom=r["chrom"], element=r["element"])
        second = astrorun.RequestBudget(p, cap=10)
        for r in astrorun.requests_not_yet_charged(plan, p):
            second.take(chrom=r["chrom"], element=r["element"])
        charged = [
            (json.loads(line)["chrom"], json.loads(line)["element"])
            for line in p.read_text().splitlines()
            if json.loads(line).get("event") == "request"
        ]
        assert len(charged) == 10, "ten elements, ten charges, no element paid for twice"
        assert sorted(charged) == [("chr1", f"E{i}") for i in range(10)]
        assert len(set(charged)) == len(charged), "no element appears in the ledger twice"
        assert second.remaining() == 0

    def test_an_unreadable_ledger_line_refuses_a_remainder_rather_than_guessing(self, tmp_path):
        p = tmp_path / "l.jsonl"
        plan = self.plan(3)
        self.ledger_charging(p, plan[:1])
        with p.open("a") as fh:
            fh.write('{"event": "request", "chrom": "chr1", "element": "E1"')
        with pytest.raises(ValueError, match="cannot be read"):
            astrorun.requests_not_yet_charged(plan, p)

    def test_the_helper_says_it_is_not_an_authorisation_to_resume(self):
        assert "not itself an authorisation" in astrorun.A_RESUME_STILL_NEEDS_ALBERTS_WORD


class TestTheGuardIsAProtectionAndNotADiagnosis:
    """By removal: strip the resume guard from a copy of the real module and show the harm happens.

    The standing rule here is that a guard must be shown to prevent something. A guard whose harm
    cannot be demonstrated is a diagnosis and must be described as one. So this takes the actual
    source of `genomeos/attribution/astrorun.py`, removes the one line that makes the budget resume
    from the ledger, and runs the same scenario against both: the real class refuses, the stripped
    class spends the authorised number a second time.
    """

    #: The guard, as one line of the real source. If this line moves, the test fails loudly rather
    #: than silently passing against a module it did not actually modify.
    GUARD_LINE = "self.sent = max(sent, self.resumed_from_ledger)"

    def stripped_module(self):
        src = Path(astrorun.__file__).read_text()
        assert src.count(self.GUARD_LINE) == 1, "the guard line has moved; this test must be updated"
        stripped = src.replace(self.GUARD_LINE, "self.sent = sent")
        # the module resolves a repository root from __file__, which exec does not provide by default
        ns: dict = {"__file__": astrorun.__file__}
        exec(compile(stripped, "astrorun_without_the_resume_guard", "exec"), ns)
        return ns

    def ledger_holding_a_full_spend(self, path, cap):
        b = astrorun.RequestBudget(path, cap=cap)
        for i in range(cap):
            b.take(chrom="chr1", element=f"E{i}")
        assert b.remaining() == 0
        return path

    def test_without_the_guard_a_restart_spends_the_cap_a_second_time(self, tmp_path):
        cap = 25
        p = self.ledger_holding_a_full_spend(tmp_path / "l.jsonl", cap)
        ns = self.stripped_module()

        unguarded = ns["RequestBudget"](p, cap=cap)
        assert unguarded.sent == 0, "the stripped copy has forgotten the ledger, which is the harm"
        for i in range(cap):
            unguarded.take(chrom="chr1", element=f"second-spend-{i}")

        charges = ns["ledger_charges"](p)["request_lines"]
        assert charges == 2 * cap == 50, (
            "the harm is real and measured: a cap of 25 bought 50 requests across two processes, "
            "which is what a lane killed and restarted would have done to Albert's 1,322"
        )

    def test_with_the_guard_the_same_restart_buys_nothing(self, tmp_path):
        cap = 25
        p = self.ledger_holding_a_full_spend(tmp_path / "l.jsonl", cap)
        guarded = astrorun.RequestBudget(p, cap=cap)
        assert guarded.sent == cap
        with pytest.raises(astrorun.CapRefusedError):
            guarded.take(chrom="chr1", element="second-spend-0")
        assert astrorun.ledger_charges(p)["request_lines"] == cap == 25, "no second spend"

    def test_the_two_halves_differ_only_by_the_guard(self, tmp_path):
        """The comparison is only worth something if the stripped copy is otherwise the real code."""
        ns = self.stripped_module()
        assert ns["AUTHORISED_REQUESTS"] == astrorun.AUTHORISED_REQUESTS == 1322
        assert ns["CAP_IS_A_REFUSAL"] == astrorun.CAP_IS_A_REFUSAL
        b = ns["RequestBudget"](tmp_path / "fresh.jsonl", cap=2)
        b.take(chrom="chr1", element="E0")
        b.take(chrom="chr1", element="E1")
        with pytest.raises(ns["CapRefusedError"]):
            b.take(chrom="chr1", element="E2")


class TestTheRequestList:
    class FakeTable:
        """Stands in for crispri.DeletionTable: the overlap condition, over a handful of elements."""

        def __init__(self, elements):
            self.elements = elements

        def overlapping(self, chrom, start, end):
            return [e for e in self.elements if e["chrom"] == chrom and e["end"] > start and e["start"] < end]

    def pair(self, element, chrom, start, end, gene, label="negative"):
        return {
            "element": element,
            "chrom": chrom,
            "start": start,
            "end": end,
            "gene": gene,
            "label": label,
        }

    def test_one_row_per_registry_element_not_per_pair(self):
        table = self.FakeTable([{"id": "R1", "chrom": "chr1", "start": 100, "end": 400}])
        pairs = [
            self.pair("S1", "chr1", 150, 250, "AAA", "positive"),
            self.pair("S1", "chr1", 150, 250, "BBB", "negative"),
        ]
        plan = astrorun.plan_requests(pairs, table)
        assert len(plan) == 1, "two pairs over one registry element is one request"
        assert plan[0]["element"] == "R1"
        assert (plan[0]["start"], plan[0]["end"]) == (100, 400), "the registry element's own span"
        assert plan[0]["serves_genes"] == ["AAA", "BBB"]
        assert plan[0]["serves_labels"] == ["negative", "positive"]

    def test_a_pair_overlapping_nothing_contributes_no_request(self):
        table = self.FakeTable([{"id": "R1", "chrom": "chr1", "start": 100, "end": 200}])
        plan = astrorun.plan_requests([self.pair("S1", "chr1", 5_000, 5_100, "AAA")], table)
        assert plan == []

    def test_one_pair_overlapping_two_elements_is_two_requests(self):
        table = self.FakeTable(
            [
                {"id": "R1", "chrom": "chr1", "start": 100, "end": 200},
                {"id": "R2", "chrom": "chr1", "start": 180, "end": 300},
            ]
        )
        plan = astrorun.plan_requests([self.pair("S1", "chr1", 150, 250, "AAA")], table)
        assert [r["element"] for r in plan] == ["R1", "R2"]

    def test_the_list_is_deterministic_and_sorted_so_two_dry_runs_can_be_diffed(self):
        table = self.FakeTable(
            [
                {"id": "R2", "chrom": "chr1", "start": 900, "end": 1_000},
                {"id": "R1", "chrom": "chr1", "start": 100, "end": 200},
                {"id": "R3", "chrom": "chr10", "start": 100, "end": 200},
            ]
        )
        pairs = [
            self.pair("S2", "chr1", 950, 960, "BBB"),
            self.pair("S1", "chr1", 150, 160, "AAA"),
            self.pair("S3", "chr10", 150, 160, "CCC"),
        ]
        first = astrorun.plan_requests(pairs, table)
        assert [r["element"] for r in first] == ["R1", "R2", "R3"], "chr1 before chr10, then start"
        assert astrorun.plan_requests(list(reversed(pairs)), table) == first

    def test_elements_are_not_filtered_by_label_because_the_registered_total_is_not(self):
        table = self.FakeTable([{"id": "R1", "chrom": "chr1", "start": 100, "end": 400}])
        plan = astrorun.plan_requests(
            [self.pair("S1", "chr1", 150, 250, "AAA", "excluded_underpowered")], table
        )
        assert len(plan) == 1
        assert plan[0]["serves_labels"] == ["excluded_underpowered"]


class TestTheDryRunCannotSend:
    """The dry run's no-network property is read off its own source, not promised in its docstring."""

    SCRIPTS = (
        Path(__file__).resolve().parents[1] / "scripts" / "astroreg_dry_run.py",
        Path(__file__).resolve().parents[1] / "scripts" / "astroreg_activity_gate.py",
    )

    #: Every top-level module through which this project could reach AlphaGenome or the network. The
    #: check is on what a file IMPORTS, not on what its prose mentions: these files discuss
    #: AlphaGenome at length and must do so, and a word in a docstring cannot send anything. An
    #: import can.
    FORBIDDEN_IMPORTS = frozenset(
        {
            "alphagenome",
            "requests",
            "httpx",
            "urllib",
            "http",
            "socket",
            "aiohttp",
            "grpc",
            "google",
        }
    )

    def imported_roots(self, path):
        """Every top-level module name the file imports, by its syntax tree and not by grep."""
        import ast

        roots = set()
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Import):
                roots.update(a.name.split(".")[0] for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
                roots.add(node.module.split(".")[0])
        return roots

    def test_no_file_in_this_lane_imports_anything_that_could_send_a_request(self):
        for path in (*self.SCRIPTS, Path(astrorun.__file__)):
            offending = self.imported_roots(path) & self.FORBIDDEN_IMPORTS
            assert not offending, f"{path.name} imports {sorted(offending)} and so could spend money"

    def test_the_dry_run_declares_no_send_flag(self):
        src = self.SCRIPTS[0].read_text()
        assert '"--send"' not in src and "'--send'" not in src
        assert '"--dry-run"' in src

    def test_the_check_would_notice_an_import_that_could_send(self, tmp_path):
        """By removal: the check is only worth something if it fails on a file that can send."""
        bad = tmp_path / "sender.py"
        bad.write_text("import alphagenome\n")
        assert self.imported_roots(bad) & self.FORBIDDEN_IMPORTS == {"alphagenome"}


class TestTheActivityInputGate:
    def rpm(self):
        return {
            "source": "EPCrisprBenchmark column",
            "units": astrorun.FROZEN_ACTIVITY_UNITS,
            "in_frozen_units": True,
        }

    def peaks(self):
        return {
            "source": "ENCFF874OPW narrowPeak",
            "units": "narrowPeak signalValue",
            "in_frozen_units": False,
        }

    def test_both_columns_in_their_own_units_is_a_go(self):
        v = astrorun.activity_verdict(dict.fromkeys(astrorun.FROZEN_ACTIVITY_COLUMNS, self.rpm()))
        assert v["go"] is True
        assert v["columns_not_in_frozen_units"] == []
        assert astrorun.requests_permitted(v) == 1322

    def test_one_column_out_of_units_is_a_no_go_that_names_it(self):
        v = astrorun.activity_verdict({"DHS.RPM": self.peaks(), "H3K27ac.RPM": self.rpm()})
        assert v["go"] is False
        assert v["columns_not_in_frozen_units"] == ["DHS.RPM"]
        assert "DHS.RPM" in v["reading"]
        assert astrorun.requests_permitted(v) == 0

    def test_a_no_go_quotes_the_registration_rule_rather_than_a_new_one(self):
        v = astrorun.activity_verdict(dict.fromkeys(astrorun.FROZEN_ACTIVITY_COLUMNS, self.peaks()))
        assert v["go"] is False
        assert v["registered_rule"] == astrorun.REGISTERED_MISSING_INPUT_RULE
        assert "no substitute is put in its place" in v["reading"]
        assert sorted(v["columns_not_in_frozen_units"]) == ["DHS.RPM", "H3K27ac.RPM"]

    def test_a_column_nobody_looked_for_is_refused_not_read_as_absent(self):
        with pytest.raises(ValueError, match="unfinished check"):
            astrorun.activity_verdict({"DHS.RPM": self.rpm()})

    def test_a_no_go_permits_no_spend_at_all(self):
        v = astrorun.activity_verdict(dict.fromkeys(astrorun.FROZEN_ACTIVITY_COLUMNS, self.peaks()))
        assert astrorun.requests_permitted(v) == 0, (
            "a no-go may not buy a deletion cache no registered claim can be made from"
        )


class TestTheSendableSubset:
    """1,232 is a FILTER of the 1,322, in one function, so a sender and a report cannot disagree."""

    def row(self, element, labels):
        return {"chrom": "chr1", "element": element, "start": 1, "end": 2, "serves_labels": labels}

    def test_a_row_serving_a_positive_or_a_negative_is_kept(self):
        plan = [self.row("A", ["positive"]), self.row("B", ["negative"])]
        assert [r["element"] for r in astrorun.requests_serving_scored_pairs(plan)] == ["A", "B"]

    def test_a_row_serving_only_held_apart_or_excluded_labels_is_dropped(self):
        plan = [
            self.row("A", ["excluded_underpowered"]),
            self.row("B", ["increase_held_apart"]),
            self.row("C", ["excluded_underpowered", "increase_held_apart"]),
        ]
        assert astrorun.requests_serving_scored_pairs(plan) == []

    def test_a_mixed_row_is_kept_because_it_serves_a_scored_pair_too(self):
        plan = [self.row("A", ["excluded_underpowered", "negative"])]
        assert len(astrorun.requests_serving_scored_pairs(plan)) == 1

    def test_the_plan_order_is_preserved(self):
        plan = [self.row(c, ["negative"]) for c in "DCBA"]
        assert [r["element"] for r in astrorun.requests_serving_scored_pairs(plan)] == list("DCBA")

    def test_the_filter_never_adds_a_row(self):
        plan = [self.row("A", ["positive"]), self.row("B", ["excluded_underpowered"])]
        assert len(astrorun.requests_serving_scored_pairs(plan)) <= len(plan)

    def test_the_scored_labels_come_from_the_caller(self):
        plan = [self.row("A", ["increase_held_apart"])]
        assert len(astrorun.requests_serving_scored_pairs(plan, ("increase_held_apart",))) == 1

    def test_it_says_it_is_a_filter_and_not_a_second_enumeration(self):
        assert "second enumeration" in astrorun.ONE_LIST_THEN_ONE_FILTER

    def test_narrowing_is_not_presented_as_changing_the_authorised_total(self):
        assert "only he can authorise" in astrorun.requests_serving_scored_pairs.__doc__


class TestAstroreg2HasItsOwnCapAndAuthorisation:
    """The old approval named the original registration by hash and must not be able to fund this one."""

    def test_the_astroreg2_cap_is_1232_and_not_the_original_1322(self):
        assert astrorun.ASTROREG2_CAP == 1232
        assert astrorun.AUTHORISED_REQUESTS == 1322
        assert astrorun.ASTROREG2_CAP != astrorun.AUTHORISED_REQUESTS

    def test_the_authorisation_is_read_only_through_the_function_that_can_refuse(self):
        """Converted from an assertion that the slot was EMPTY, which was a fact about today.

        The slot has since been filled, and that must not change what the code does. What is
        load-bearing is that nothing reads the constant directly: every reader goes through
        astroreg2_authorisation(), which is the only place an absent approval can refuse.
        """
        import inspect

        src = inspect.getsource(astrorun.astroreg2_budget)
        assert "astroreg2_authorisation()" in src, (
            "the budget must obtain the approval through the refusing function, not by reading the "
            "constant, or an empty slot would be permissive"
        )

    def test_an_ABSENT_authorisation_refuses_whatever_the_slot_holds_today(self, monkeypatch):
        """Mechanism, not state: patched to None, the refusal fires. Holds before or after approval."""
        monkeypatch.setattr(astrorun, "ASTROREG2_AUTHORISATION", None)
        with pytest.raises(astrorun.NoAuthorisationError, match="does not carry"):
            astrorun.astroreg2_authorisation()
        monkeypatch.setattr(astrorun, "ASTROREG2_AUTHORISATION", "")
        with pytest.raises(astrorun.NoAuthorisationError):
            astrorun.astroreg2_authorisation()  # an empty string is not an approval either

    def test_no_budget_can_be_obtained_without_an_approval_whatever_the_slot_holds_today(
        self, tmp_path, monkeypatch
    ):
        """The refusal comes BEFORE the budget exists, so there is nothing to send with."""
        monkeypatch.setattr(astrorun, "ASTROREG2_AUTHORISATION", None)
        with pytest.raises(astrorun.NoAuthorisationError):
            astrorun.astroreg2_budget(run_id=1, root=tmp_path)
        assert not astrorun.ledger_path_for_run(1, tmp_path).exists(), "a refusal creates no ledger"

    def test_the_astroreg2_scope_cannot_reach_the_original_1322_cap(self, tmp_path, monkeypatch):
        """With an approval recorded, the cap is 1,232: request 1,233 is refused, 1,322 unreachable."""
        monkeypatch.setattr(astrorun, "ASTROREG2_AUTHORISATION", "<verbatim words would go here>")
        monkeypatch.setattr(
            astrorun,
            "ASTROREG2_AUTHORISATIONS",
            {1: {"words": "<verbatim>", "requests": astrorun.ASTROREG2_CAP, "consumed": False}},
        )
        (tmp_path / astrorun.LEDGER_RUN1.parent).mkdir(parents=True, exist_ok=True)
        b = astrorun.astroreg2_budget(run_id=1, root=tmp_path)
        assert b.cap == 1232
        for _ in range(1232):
            b.take(chrom="chr1", element="e")
        assert b.remaining() == 0
        with pytest.raises(astrorun.CapRefusedError) as exc:
            b.take(chrom="chr1", element="one too many")
        assert "1232 of 1232 sent" in str(exc.value)
        assert b.sent == 1232 < astrorun.AUTHORISED_REQUESTS

    def test_the_reason_absence_refuses_rather_than_warns_is_recorded(self):
        assert "Not a warning, a refusal" in astrorun.NO_AUTHORISATION_MEANS_NO_SEND


class TestTheSenderProvesItSendsTheReviewedList:
    def plan(self, n=5):
        return [
            {
                "chrom": "chr1",
                "element": f"E{i}",
                "start": i * 100,
                "end": i * 100 + 50,
                "serves_genes": ["AAA"],
                "serves_labels": ["negative"],
            }
            for i in range(n)
        ]

    def test_the_same_list_has_the_same_digest(self):
        a, b = self.plan(), self.plan()
        assert astrorun.plan_digest(a) == astrorun.plan_digest(b)
        assert astrorun.check_is_the_reviewed_plan(b, astrorun.plan_digest(a))

    def test_formatting_does_not_move_the_digest(self):
        """The digest is over the list's content, not the file's bytes."""
        a = self.plan()
        b = [dict(reversed(list(r.items()))) for r in a]
        assert astrorun.plan_digest(b) == astrorun.plan_digest(a)

    def test_a_planted_one_element_change_is_refused(self):
        """By planting: change exactly one request and show the refusal fire."""
        reviewed = astrorun.plan_digest(self.plan())
        tampered = self.plan()
        tampered[2]["element"] = "E2-but-different"
        with pytest.raises(ValueError, match="not the reviewed"):
            astrorun.check_is_the_reviewed_plan(tampered, reviewed)

    def test_a_planted_one_coordinate_change_is_refused(self):
        reviewed = astrorun.plan_digest(self.plan())
        tampered = self.plan()
        tampered[0]["end"] += 1
        with pytest.raises(ValueError, match="not the reviewed"):
            astrorun.check_is_the_reviewed_plan(tampered, reviewed)

    def test_an_added_or_removed_request_is_refused(self):
        reviewed = astrorun.plan_digest(self.plan())
        with pytest.raises(ValueError, match="not the reviewed"):
            astrorun.check_is_the_reviewed_plan(self.plan(6), reviewed)
        with pytest.raises(ValueError, match="not the reviewed"):
            astrorun.check_is_the_reviewed_plan(self.plan(4), reviewed)

    def test_reordering_is_refused_because_the_order_is_part_of_the_list(self):
        reviewed = astrorun.plan_digest(self.plan())
        with pytest.raises(ValueError, match="not the reviewed"):
            astrorun.check_is_the_reviewed_plan(list(reversed(self.plan())), reviewed)

    def test_a_changed_label_set_is_refused_because_it_changes_what_the_request_is_for(self):
        reviewed = astrorun.plan_digest(self.plan())
        tampered = self.plan()
        tampered[1]["serves_labels"] = ["positive"]
        with pytest.raises(ValueError, match="not the reviewed"):
            astrorun.check_is_the_reviewed_plan(tampered, reviewed)

    def test_the_refusal_names_the_count_so_a_reviewer_can_see_what_changed(self):
        reviewed = astrorun.plan_digest(self.plan())
        with pytest.raises(ValueError) as exc:
            astrorun.check_is_the_reviewed_plan(self.plan(7), reviewed)
        assert "7 requests" in str(exc.value)

    def test_the_reason_the_sender_recomputes_is_recorded(self):
        assert "a review of nothing" in astrorun.THE_REVIEWED_LIST_IS_THE_SENT_LIST


class TestAlbertsConditionsAreEachTheirOwnRefusal:
    """By PLANTING, with the authorisation string PRESENT in every test.

    The failure this guards against is the one found twice already today: a guard that protects the
    prose and lets the substance through. If any of these passed with the string recorded, Albert's
    conditions would have become decoration.
    """

    DIGEST_SOURCE = [
        {
            "chrom": "chr1",
            "element": f"E{i}",
            "start": i,
            "end": i + 1,
            "serves_genes": ["G"],
            "serves_labels": ["negative"],
        }
        for i in range(astrorun.ASTROREG2_CAP)
    ]

    @pytest.fixture
    def authorised(self, monkeypatch):
        """Albert's words recorded, and the blob check stubbed.

        The blob check is the SUPERVISOR's requirement, not a clause of Albert's, and it refuses while
        this lane holds changes to the signed files. Stubbing it keeps these tests about HIS clauses; it
        has its own planted tests, including the one-character edit.
        """
        monkeypatch.setattr(astrorun, "ASTROREG2_AUTHORISATION", astrorun.ASTROREG2_AUTHORISATION_AS_RELAYED)
        monkeypatch.setattr(astrorun, "check_signoff_closure", lambda *a, **k: {})
        # Run 1's REAL approval is marked consumed, because run 1 happened under it. These tests are
        # about Albert's other clauses, so run 1 is granted an unconsumed approval here; the per-run
        # approval rule has its own planted tests, including run 2 under run 1's words.
        monkeypatch.setattr(
            astrorun,
            "ASTROREG2_AUTHORISATIONS",
            {1: {"words": "<run 1's words>", "requests": astrorun.ASTROREG2_CAP, "consumed": False}},
        )

    @pytest.fixture
    def good(self, tmp_path):
        act = tmp_path / "activity.json"
        act.write_text(json.dumps({"rule": {"producer": dict(astrorun.AMENDMENT_2_RULE_FINGERPRINT)}}))
        reg = tmp_path / "astroreg2_registration.json"
        reg.write_text("{}")
        (tmp_path / astrorun.LEDGER_RUN1.parent).mkdir(parents=True, exist_ok=True)
        return {
            "activity_result": act,
            "registration": reg,
            "run_id": 1,
            "ledger_root": tmp_path,
            # INJECTED, never read from the volume: item (g)'s disk clause still RUNS here, it is
            # simply handed the figure instead of asking the laptop. See
            # astrorun.THE_FIGURE_IS_INJECTED_IN_TESTS_AND_READ_ONLY_IN_A_REAL_RUN.
            "free_bytes": 11 * 1024**3,
            "signoff": 'the supervisor wrote "dry run reviewed" at 2026-10-02T13:00:00',
            "plan": self.DIGEST_SOURCE,
            "reviewed_digest": astrorun.plan_digest(self.DIGEST_SOURCE),
            "committed": lambda p: True,
        }

    def test_with_everything_satisfied_it_returns(self, authorised, good, monkeypatch):
        """The positive control: without it, a refusal proves nothing.

        Item (f), the supervisor's adapter-v2 requirement, is stubbed here because it is a different
        condition from Albert's five and has its own planted tests. Stubbing it is what makes this a
        test of HIS clauses rather than of the adapter.
        """
        monkeypatch.setattr(astrorun, "check_adapter_v2", lambda *a, **k: {"stubbed": True})
        # Item (g) likewise: the supervisor's, with its own planted tests, and it refuses while the
        # recording path still keeps four cell lines. Stubbing it keeps this about Albert's clauses.
        monkeypatch.setattr(astrorun, "check_adapter_writes_full_vectors", lambda *a, **k: {})
        out = astrorun.may_send(**good)
        assert out["may_send"] is True
        assert out["cap"] == 1232
        assert sorted(out["clauses_checked"]) == ["activity", "logged", "one_run", "scope", "signoff"]

    def test_the_relayed_text_is_kept_apart_and_never_consulted_by_the_refusing_function(self):
        """Converted: the old premise (slot empty) is obsolete; these claims are still load-bearing.

        The relayed text remains a SEPARATE constant, the function that can refuse never reads it, and
        the slot's content is byte-equal to it -- which is a cross-check between two transcriptions of
        the same relay, not independent corroboration that Albert said it. It confirms the relay was
        carried faithfully and nothing more.
        """
        import inspect

        assert "I approve 1,232" in astrorun.ASTROREG2_AUTHORISATION_AS_RELAYED
        assert "not that person's approval" in astrorun.WHY_THE_RELAYED_TEXT_IS_NOT_THE_APPROVAL
        src = inspect.getsource(astrorun.astroreg2_authorisation)
        assert "ASTROREG2_AUTHORISATION_AS_RELAYED" not in src, (
            "the refusing function must never fall back to the relayed text"
        )
        assert astrorun.ASTROREG2_AUTHORISATION == astrorun.ASTROREG2_AUTHORISATION_AS_RELAYED, (
            "the recorded approval and this lane's transcription of the relay agree character for character"
        )

    def test_without_an_authorisation_recorded_nothing_sends_whatever_the_slot_holds_today(
        self, good, monkeypatch
    ):
        monkeypatch.setattr(astrorun, "ASTROREG2_AUTHORISATION", None)
        with pytest.raises(astrorun.NoAuthorisationError):
            astrorun.may_send(**good)

    def test_a_FILLED_slot_alone_does_not_permit_a_send_while_another_clause_fails(
        self, good, tmp_path, monkeypatch
    ):
        """What nothing covered until the slot was filled.

        Until an hour ago an empty slot refused first and masked every other clause, so a regression in
        one of them could not have been seen. With the real recorded approval in place and NO patching
        of it, a failing clause must still refuse -- here the sign-off, removed.
        """
        monkeypatch.setattr(astrorun, "check_adapter_v2", lambda *a, **k: {"stubbed": True})
        monkeypatch.setattr(astrorun, "check_signoff_closure", lambda *a, **k: {})
        assert astrorun.ASTROREG2_AUTHORISATION, "this test is about a FILLED slot"
        # The slot stays REAL and unpatched, which is the point; what is granted is run 1's per-run
        # approval, since the real one is consumed and would otherwise mask the clause under test --
        # which is the same masking this test exists to catch.
        grant_run(monkeypatch)
        good["signoff"] = None
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.may_send(**good)
        assert "dry run reviewed" in str(exc.value)

    def test_a_filled_slot_does_not_mask_a_partial_run_either(self, good, monkeypatch):
        """The same masking risk on the clause that guards a double spend."""
        monkeypatch.setattr(astrorun, "check_adapter_v2", lambda *a, **k: {"stubbed": True})
        monkeypatch.setattr(astrorun, "check_signoff_closure", lambda *a, **k: {})
        grant_run(monkeypatch)
        led = astrorun.ledger_path_for_run(good["run_id"], good["ledger_root"])
        astrorun.RequestBudget(led, cap=astrorun.ASTROREG2_CAP).take(chrom="chr1", element="E0")
        with pytest.raises(astrorun.SendRefusedError, match="a partial run exists"):
            astrorun.may_send(**good)

    def test_PLANTED_activity_result_absent(self, authorised, good, tmp_path):
        good["activity_result"] = tmp_path / "missing.json"
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.may_send(**good)
        assert "under Amendment 2" in str(exc.value), "the refusal must quote the clause"
        assert "not present" in str(exc.value)

    def test_PLANTED_activity_result_present_but_not_committed(self, authorised, good):
        good["committed"] = lambda p: False
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.may_send(**good)
        assert "NOT COMMITTED" in str(exc.value)
        assert "under Amendment 2" in str(exc.value)

    def test_PLANTED_activity_produced_under_a_different_rule(self, authorised, good):
        wrong = dict(astrorun.AMENDMENT_2_RULE_FINGERPRINT) | {"commit": "0" * 40}
        Path(good["activity_result"]).write_text(json.dumps({"rule": {"producer": wrong}}))
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.may_send(**good)
        assert "different rule" in str(exc.value)
        assert "0000000" in str(exc.value), "the refusal names what it found"

    def test_PLANTED_activity_naming_no_rule_at_all(self, authorised, good):
        Path(good["activity_result"]).write_text(json.dumps({"rule": {}}))
        with pytest.raises(astrorun.SendRefusedError):
            astrorun.may_send(**good)

    def test_PLANTED_supervisor_signoff_absent(self, authorised, good):
        good["signoff"] = None
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.may_send(**good)
        assert "dry run reviewed" in str(exc.value)
        assert "may not be anticipated" in str(exc.value)

    def test_PLANTED_signoff_paraphrased_rather_than_quoted(self, authorised, good):
        good["signoff"] = "the supervisor approved the dry run"
        with pytest.raises(astrorun.SendRefusedError):
            astrorun.may_send(**good)

    def test_PLANTED_a_second_run_after_a_completed_one(self, authorised, good):
        astrorun.mark_run_complete(
            astrorun.ledger_path_for_run(good["run_id"], good["ledger_root"]), requests=1232
        )
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.may_send(**good)
        assert "One run" in str(exc.value)
        assert "already records a completed run" in str(exc.value)

    def test_a_completed_run_is_terminal_and_not_an_exhausted_counter(self, tmp_path):
        """The two facts are different and a reset would look like the former."""
        led = tmp_path / "l.jsonl"
        b = astrorun.RequestBudget(led, cap=3)
        b.take(chrom="chr1", element="e")
        assert astrorun.run_already_completed(led) is False, "a partial run is not a completed one"
        astrorun.mark_run_complete(led, requests=1)
        assert astrorun.run_already_completed(led) is True
        assert "terminal state and not as an exhausted counter" in astrorun.ONE_RUN_IS_TERMINAL

    def test_PLANTED_the_list_is_the_wrong_length_for_the_approved_scope(self, authorised, good):
        good["plan"] = self.DIGEST_SOURCE[:-1]
        good["reviewed_digest"] = astrorun.plan_digest(good["plan"])
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.may_send(**good)
        assert "1231 requests and the approval names 1232" in str(exc.value)

    def test_PLANTED_one_element_changed_so_it_is_not_the_reviewed_list(self, authorised, good):
        tampered = [dict(r) for r in self.DIGEST_SOURCE]
        tampered[5]["element"] = "E5-but-different"
        good["plan"] = tampered
        with pytest.raises(ValueError, match="not the reviewed"):
            astrorun.may_send(**good)

    def test_PLANTED_the_registration_the_approval_names_is_absent(self, authorised, good, tmp_path):
        good["registration"] = tmp_path / "gone.json"
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.may_send(**good)
        assert "astroreg2_registration at 0f4c372" in str(exc.value)

    def test_every_refusal_quotes_the_clause_it_enforces(self, authorised, good, tmp_path):
        """A refusal saying only 'not authorised' teaches nothing and can be cleared by accident."""
        good["signoff"] = None
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.may_send(**good)
        assert "Albert's approval is CONDITIONAL on" in str(exc.value)


class TestAmendment2sRuleIsUnchangedSinceItWasRegistered:
    """Item (b): shown by hash, not asserted.

    The rule is registered across TWO files and a reader should not be told it is in one: 0f4c372 names
    the repository, the commit and the no-rescaling term; 230efc8 carries the term-by-term citations and
    the file name. Both are recorded.
    """

    def test_the_cited_rule_still_hashes_to_what_the_registration_recorded(self):
        import hashlib

        from genomeos.attribution import rpm

        got = hashlib.sha256(json.dumps(rpm.PRODUCER, sort_keys=True).encode()).hexdigest()
        assert got == astrorun.AMENDMENT_2_CITATIONS_SHA256, (
            "rpm.PRODUCER has changed since amendment 2 was registered. The send path binds to this "
            "rule, so a change here is a change to what Albert's condition means"
        )

    def test_the_fingerprint_the_send_path_binds_to_is_a_subset_of_the_cited_rule(self):
        from genomeos.attribution import rpm

        for key, want in astrorun.AMENDMENT_2_RULE_FINGERPRINT.items():
            assert rpm.PRODUCER[key] == want

    def test_where_the_rule_is_registered_is_recorded_for_both_files(self):
        where = astrorun.AMENDMENT_2_RULE_REGISTERED_IN
        assert any("0f4c372" in k for k in where)
        assert any("230efc8" in k for k in where)
        assert "does NOT name the file" in where["data/results/astroreg2_registration.json at 0f4c372"]


class TestOneApprovalPerRunAndOneTotalAcrossThem:
    """The per-run ledger opened an UNBOUNDED path, and these are the four refusals that close it.

    Shown before it was fixed, on 20fe5ee's body: `--run N` resolves a FRESH, EMPTY ledger, so
    run_already_completed is False and charges are 0, so the one-run clause passed for run 2, then 3,
    then 4, under run 1's single consumed approval -- 4,928 requests permitted across four further runs
    with the approval text unchanged, and no stopping point. The cap of 1,232 was PER LEDGER, so there
    was no total bound at all.

    The shape is the day's third of one kind: a mechanism that existed and was not wired to the thing it
    protected. ONE_RUN asked "is this ledger finished?" when the question is "has Albert authorised THIS
    run?" -- and a guard keyed to a path the caller chooses cannot ask it. So may_send now takes a run
    id and DERIVES the ledger.
    """

    CAP = 3

    def plan(self, n=None):
        return [
            {
                "chrom": "chr21",
                "element": f"e{i}",
                "start": i * 100,
                "end": i * 100 + 50,
                "serves_genes": ["G"],
            }
            for i in range(self.CAP if n is None else n)
        ]

    def kwargs(self, tmp_path, run_id, plan=None, signoff_seen=None):
        plan = self.plan() if plan is None else plan
        act = tmp_path / "activity.json"
        act.write_text(json.dumps({"rule": {"producer": dict(astrorun.AMENDMENT_2_RULE_FINGERPRINT)}}))
        reg = tmp_path / "reg.json"
        reg.write_text("{}")
        (tmp_path / astrorun.LEDGER_RUN1.parent).mkdir(parents=True, exist_ok=True)
        return {
            "activity_result": act,
            "registration": reg,
            "run_id": run_id,
            "ledger_root": tmp_path,
            "free_bytes": 11 * 1024**3,  # injected, so the run rule is not decided by free disk
            "signoff": 'the supervisor wrote "dry run reviewed" at 2026-10-02T20:05:00',
            "plan": plan,
            "reviewed_digest": astrorun.plan_digest(plan),
            "committed": lambda p: True,
        }

    def complete_run(self, tmp_path, run, charges=1, cap=None):
        """A run that ended: charges plus a run_complete line, which is what makes it terminal."""
        led = astrorun.ledger_path_for_run(run, tmp_path)
        led.parent.mkdir(parents=True, exist_ok=True)
        budget = astrorun.RequestBudget(led, cap=self.CAP if cap is None else cap)
        for i in range(charges):
            budget.take(chrom="chr21", element=f"old{i}")
        astrorun.mark_run_complete(led, requests=charges)
        return led

    def approvals(self, monkeypatch, live=(), spent=(1,)):
        rec = {}
        for run in spent:
            rec[run] = {
                "words": f"<run {run}'s words>",
                "requests": self.CAP,
                "consumed": True,
                "why_consumed": f"run {run} happened under it",
            }
        for run in live:
            rec[run] = {"words": f"<run {run}'s words>", "requests": self.CAP, "consumed": False}
        monkeypatch.setattr(astrorun, "ASTROREG2_AUTHORISATIONS", rec)
        monkeypatch.setattr(astrorun, "ASTROREG2_CAP", self.CAP)

    @pytest.fixture(autouse=True)
    def other_clauses(self, monkeypatch):
        """The clauses this class is NOT about, stubbed so a refusal here names the run rule."""
        monkeypatch.setattr(astrorun, "ASTROREG2_AUTHORISATION", astrorun.ASTROREG2_AUTHORISATION_AS_RELAYED)
        monkeypatch.setattr(astrorun, "check_adapter_v2", lambda *a, **k: {"stubbed": True})
        monkeypatch.setattr(astrorun, "check_nothing_in_the_plan_is_already_cached", lambda *a, **k: None)
        # Items (f) and (g) are the supervisor's and have their own planted tests; (g) refuses today,
        # so stubbing it is what lets this class be about the run rule it is named for.
        monkeypatch.setattr(astrorun, "check_adapter_writes_full_vectors", lambda *a, **k: {})

    def test_PLANTED_run_2_with_only_run_1s_approval_REFUSES(self, tmp_path, monkeypatch):
        """The hole itself: a fresh ledger used to make run 2 look like a run nobody had done yet."""
        monkeypatch.setattr(astrorun, "check_signoff_closure", lambda *a, **k: {})
        self.approvals(monkeypatch, live=(), spent=(1,))
        self.complete_run(tmp_path, 1)
        with pytest.raises(astrorun.NoAuthorisationError) as exc:
            astrorun.may_send(**self.kwargs(tmp_path, 2))
        said = str(exc.value)
        assert "no approval is recorded for run 2" in said
        assert "a run with no approval of its own does not start" in said

    def test_PLANTED_run_1s_own_approval_is_CONSUMED_so_even_run_1_cannot_go_again(
        self, tmp_path, monkeypatch
    ):
        """An approval is spent by the run that happened under it, recorded rather than inferred."""
        monkeypatch.setattr(astrorun, "check_signoff_closure", lambda *a, **k: {})
        self.approvals(monkeypatch, live=(), spent=(1,))
        with pytest.raises(astrorun.NoAuthorisationError) as exc:
            astrorun.may_send(**self.kwargs(tmp_path, 1))
        assert "CONSUMED" in str(exc.value)
        assert "A further run needs a NEW approval in Albert's own words" in str(exc.value)

    def test_the_REAL_record_marks_run_1s_approval_spent_and_says_why(self):
        """The fact about today, read from the record rather than from a ledger."""
        rec = astrorun.ASTROREG2_AUTHORISATIONS[1]
        assert rec["consumed"] is True
        assert rec["words"] == astrorun.ASTROREG2_AUTHORISATION, "run 1 holds the existing text"
        assert rec["requests"] == 1232
        assert "not a second run, and not a resume" in rec["why_consumed"]
        assert set(astrorun.ASTROREG2_AUTHORISATIONS) == {1}, (
            "no approval for any later run is recorded, so no later run can start"
        )
        assert astrorun.total_authorised_requests() == 1232

    def test_PLANTED_run_3_after_a_completed_run_2_without_its_own_approval_REFUSES(
        self, tmp_path, monkeypatch
    ):
        """Completing a run must not be what authorises the next one."""
        monkeypatch.setattr(astrorun, "check_signoff_closure", lambda *a, **k: {})
        self.approvals(monkeypatch, live=(), spent=(1, 2))
        self.complete_run(tmp_path, 1)
        self.complete_run(tmp_path, 2)
        with pytest.raises(astrorun.NoAuthorisationError) as exc:
            astrorun.may_send(**self.kwargs(tmp_path, 3))
        assert "no approval is recorded for run 3" in str(exc.value)

    def test_PLANTED_run_3_SKIPPING_run_2_REFUSES(self, tmp_path, monkeypatch):
        """No skipping and no two runs at once: an open run's charges count against nothing."""
        monkeypatch.setattr(astrorun, "check_signoff_closure", lambda *a, **k: {})
        self.approvals(monkeypatch, live=(3,), spent=(1, 2))
        self.complete_run(tmp_path, 1)
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.may_send(**self.kwargs(tmp_path, 3))
        said = str(exc.value)
        assert "run 3 cannot start" in said
        assert "run_complete" in said
        assert "skip a run or put two runs on the same approvals at once" in said

    def test_PLANTED_a_signoff_for_run_2_used_for_run_3_REFUSES(self, tmp_path, monkeypatch):
        """A sign-off covers one run: its closure and its words were written about one run's spend."""
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.check_signoff_closure(record={"run_id": 2, "signed_closure": {"a": "b"}}, run_id=3)
        assert "covers run 2, not run 3" in str(exc.value)
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.check_signoff_closure(record={"signed_closure": {"a": "b"}}, run_id=3)
        assert "names no run_id" in str(exc.value)

    def test_PLANTED_a_total_across_runs_is_what_bounds_money(self, tmp_path, monkeypatch):
        """A per-ledger cap is not a budget. The total is the sum over every run ledger on disk."""
        monkeypatch.setattr(astrorun, "check_signoff_closure", lambda *a, **k: {})
        self.approvals(monkeypatch, live=(2,), spent=(1,))
        self.complete_run(tmp_path, 1, charges=self.CAP)  # run 1 used its whole allowance
        assert astrorun.total_authorised_requests() == 2 * self.CAP
        assert astrorun.total_charged_across_runs(tmp_path)["total"] == self.CAP
        # run 2's own CAP still fits: CAP already + CAP asked == the total
        out = astrorun.may_send(**self.kwargs(tmp_path, 2))
        assert out["may_send"] is True
        # but one more request than the total leaves does not
        monkeypatch.setattr(
            astrorun,
            "ASTROREG2_AUTHORISATIONS",
            {
                1: {"words": "<1>", "requests": self.CAP, "consumed": True, "why_consumed": "ran"},
                2: {"words": "<2>", "requests": self.CAP + 1, "consumed": False},
            },
        )
        monkeypatch.setattr(astrorun, "ASTROREG2_CAP", self.CAP + 1)
        big = self.plan(self.CAP + 2)
        with pytest.raises(astrorun.CapRefusedError) as exc:
            astrorun.may_send(**self.kwargs(tmp_path, 2, plan=big))
        assert "the total ever authorised is" in str(exc.value)

    def test_PLANTED_a_ledger_for_an_unapproved_run_refuses_rather_than_being_ignored(
        self, tmp_path, monkeypatch
    ):
        """Charges outside the scheme are exactly what must not be passed over."""
        monkeypatch.setattr(astrorun, "check_signoff_closure", lambda *a, **k: {})
        self.approvals(monkeypatch, live=(2,), spent=(1,))
        self.complete_run(tmp_path, 1)
        self.complete_run(tmp_path, 9)  # a ledger no approval covers
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.may_send(**self.kwargs(tmp_path, 2))
        assert "[9]" in str(exc.value)
        assert "no recorded approval covers" in str(exc.value)

    def test_the_budget_itself_is_capped_by_what_the_total_leaves(self, tmp_path, monkeypatch):
        """Enforced twice: a run must not cross the total PART WAY THROUGH, not only before it starts."""
        self.approvals(monkeypatch, live=(2,), spent=(1,))
        # run 1 charged MORE than its own allowance would suggest, so the total is what binds run 2
        self.complete_run(tmp_path, 1, charges=self.CAP + 1, cap=self.CAP + 1)
        budget = astrorun.astroreg2_budget(run_id=2, root=tmp_path)
        assert budget.cap == self.CAP - 1, (
            "the lower of this run's allowance and what the total leaves: "
            f"{2 * self.CAP} authorised in all, {self.CAP + 1} charged elsewhere"
        )
        assert str(budget.ledger).endswith(astrorun.ledger_for_run(2).name)

    def test_PLANTED_injected_9_46_GB_REFUSES_at_the_figure_it_is_HANDED(self, tmp_path, monkeypatch):
        """The figure push17's suite actually had, proved as a refusal without reading any volume.

        9.46 GB was the live reading when three tests of Albert's OTHER clauses went red. That was not
        a verdict about the code, and it is the fourth check tonight to assert the state of the live
        machine -- this one inside a MONEY guard, which makes it the most expensive of the four. A guard
        whose test cannot be run on a full disk is one somebody will eventually satisfy by deleting
        files.
        """
        kwargs = self.kwargs(tmp_path, 1)
        kwargs["free_bytes"] = int(9.46 * 1024**3)
        self.approvals(monkeypatch, live=(1,), spent=())
        monkeypatch.setattr(astrorun, "check_signoff_closure", lambda *a, **k: {})
        with pytest.raises(astrorun.VectorRefusedError) as exc:
            astrorun.may_send(**kwargs)
        said = str(exc.value)
        assert "9.46 GB is free" in said, said
        assert "would leave less than the 10 GB floor" in said
        assert "refusal and not a warning" in said
        assert "paged twice tonight" in said
        assert "item (g)" in said, "labelled the supervisor's, never attributed to Albert"

    def test_NEAR_MISS_injected_11_GB_PASSES_and_the_run_proceeds(self, tmp_path, monkeypatch):
        """The other way round, at a figure it is handed: a guard that can never pass is the defect."""
        kwargs = self.kwargs(tmp_path, 1)
        kwargs["free_bytes"] = 11 * 1024**3
        self.approvals(monkeypatch, live=(1,), spent=())
        monkeypatch.setattr(astrorun, "check_signoff_closure", lambda *a, **k: {})
        out = astrorun.may_send(**kwargs)
        assert out["may_send"] is True

    def test_NEAR_MISS_run_2_with_its_OWN_approval_and_its_OWN_signoff_PROCEEDS(self, tmp_path, monkeypatch):
        """A guard that can never pass is the fixed-point defect again, so the positive control.

        The sign-off closure is stubbed because it is the supervisor's requirement and has its own
        planted tests -- but the stub RECORDS what it was handed, so this also shows the run id
        reaching the closure check rather than the check being asked about no run at all.
        """
        seen: dict[str, object] = {}

        def closure(*a, **k):
            seen.update(k)
            return {}

        monkeypatch.setattr(astrorun, "check_signoff_closure", closure)
        self.approvals(monkeypatch, live=(2,), spent=(1,))
        self.complete_run(tmp_path, 1)
        out = astrorun.may_send(**self.kwargs(tmp_path, 2))
        assert out["may_send"] is True
        assert out["run_id"] == 2
        assert out["ledger"].endswith(astrorun.ledger_for_run(2).name)
        assert out["totals"]["total_authorised"] == 2 * self.CAP
        assert seen.get("run_id") == 2, "the sign-off check must be asked about THIS run"


class TestItemGTheFullTrackVectorIsKept:
    """Item (g): every gene row of every paid answer persists its WHOLE track vector.

    The defect: the recording path keeps each row's value multiset for CELL_TRACKS -- K562, HepG2,
    GM12878, IMR-90 -- and discards the rest of the axis. For 1,232 ASTROCYTE elements that buys every
    brain-tissue value and throws it away to keep four non-neural lines, and paid data is the one kind
    this project cannot re-fetch for free.

    It is the cache-hit defect one layer out: there the run would have charged 1,232 and bought nothing,
    here it would charge 1,232 and discard the part that answers the question it was authorised for.
    Neither was visible to a stub, because a stub returns what the code asks for.
    """

    class Matrix:
        """A stand-in for the response's X: rows x tracks, read by [gi, ti] like the real one."""

        def __init__(self, rows, tracks):
            self.values = [
                [round(0.5 - 0.01 * (gi * tracks + ti), 4) for ti in range(tracks)] for gi in range(rows)
            ]

        def __getitem__(self, key):
            gi, ti = key
            return self.values[gi][ti]

    def axis(self, rows=3, tracks=6):
        """An axis as recorded_axis builds one: the four cell lines kept, the rest of the axis gone."""
        return {
            "output": 0,
            "fields": ["gene_name", "gene_id"],
            "rows": [
                {
                    "row": gi,
                    "obs_index": f"ENSG{gi}",
                    "gene_name": f"G{gi}",
                    "gene_id": f"ENSG{gi}",
                    "tracks_total": tracks,
                    "values_emitted": tracks,
                    "cell_values": {"K562": [0.1], "HepG2": [0.2]},
                }
                for gi in range(rows)
            ],
            "tracks_total": tracks,
            "threshold": 0.0,
        }

    def names(self, tracks=6):
        return ["astrocyte", "brain cortex", "brain cerebellum", "K562", "HepG2", "GM12878"][:tracks]

    def test_the_vector_round_trips_and_comes_back_NAMED(self):
        """(c): the values back out, one per track, in the response's own column order."""
        axis = astrorun.with_full_track_vectors(self.axis(), self.Matrix(3, 6), self.names(), "a" * 64)
        pairs = astrorun.named_track_values(axis, 1)
        assert [n for n, _ in pairs] == self.names(), "names in column order"
        assert len(pairs) == 6
        assert pairs[0] == ("astrocyte", self.Matrix(3, 6).values[1][0])
        assert pairs[5] == ("GM12878", self.Matrix(3, 6).values[1][5])

    def test_PLANTED_duplicate_track_names_do_not_COLLAPSE(self):
        """The defect the REAL metadata found: 667 rna_seq rows carry only 316 distinct names.

        Returned as {name: value} a row came back with 316 entries -- 351 paid values dropped,
        silently, under a length that looked plausible. Pairs in column order cannot do that, and no
        stub could show it: a stub's names were tidy and unique.
        """
        dup = ["astrocyte", "astrocyte", "astrocyte", "K562", "K562", "K562"]
        axis = astrorun.with_full_track_vectors(self.axis(), self.Matrix(3, 6), dup, "a" * 64)
        pairs = astrorun.named_track_values(axis, 0)
        assert len(pairs) == 6, "one pair per track, duplicates included"
        assert len({n for n, _ in pairs}) == 2
        assert len(dict(pairs)) == 2, "the collapse, demonstrated rather than described"
        assert "COLLAPSES" in astrorun.NAMES_ARE_NOT_UNIQUE

    def test_PLANTED_a_TRUNCATED_vector_REFUSES_rather_than_returning_what_fits(self):
        """A short vector is a partly discarded answer, and returning its prefix would hide that."""
        full = astrorun.pack_track_vector([1.0, 2.0, 3.0, 4.0])
        with pytest.raises(astrorun.VectorRefusedError, match="carries 4 values"):
            astrorun.unpack_track_vector(full, expected=6)
        # cut mid-float, which is what a truncated write looks like
        import base64

        raw = base64.b64decode(full)[:-3]
        with pytest.raises(astrorun.VectorRefusedError, match="not a whole number of float64"):
            astrorun.unpack_track_vector(base64.b64encode(raw).decode(), expected=4)

    def test_PLANTED_a_row_without_a_vector_REFUSES_and_so_does_a_nameless_axis(self):
        axis = astrorun.with_full_track_vectors(self.axis(), self.Matrix(3, 6), self.names(), "a" * 64)
        stripped = dict(axis)
        stripped["rows"] = [
            {k: v for k, v in r.items() if k != astrorun.FULL_VECTOR_KEY} for r in axis["rows"]
        ]
        with pytest.raises(astrorun.VectorRefusedError, match="carries no track_vector"):
            astrorun.named_track_values(stripped, 0)
        nameless = {**axis, "track_names": []}
        with pytest.raises(astrorun.VectorRefusedError, match="records no track names"):
            astrorun.named_track_values(nameless, 0)

    def test_PLANTED_the_REGISTERED_FROZEN_FEATURES_DO_NOT_MOVE_BY_ONE_BIT(self):
        """(b) the hard constraint: this is additive persistence, not a change to what is computed.

        The frozen deletion term reads an element answer's `genes` rows and their per-cell values
        through crispri.deletion_values/deletion_drop. Those are computed before and after the vectors
        are added and compared EXACTLY, not approximately: a feature that moved would be a different
        feature from the one the weights were frozen on.
        """
        from genomeos.attribution import crispri

        before = self.axis()
        answer = {
            "id": "E1",
            "chrom": "chr1",
            "genes": [{"gene": "G1", "mean_log2fc": -0.4, "n_tracks": 6, "by_cell": {"astrocyte": -0.7}}],
            "predicted": {"gene": "G1"},
            "predicted_by_cell": {"astrocyte": -0.7},
            "model": {"gene_axis_outputs": [before]},
        }
        top_a, vals_a = crispri.deletion_values([answer], "G1", "astrocyte")
        drop_a = crispri.deletion_drop(vals_a)

        after = astrorun.with_full_track_vectors(before, self.Matrix(3, 6), self.names(), "a" * 64)
        widened = {**answer, "model": {"gene_axis_outputs": [after]}}
        top_b, vals_b = crispri.deletion_values([widened], "G1", "astrocyte")
        drop_b = crispri.deletion_drop(vals_b)

        assert (top_a, vals_a, drop_a) == (top_b, vals_b, drop_b), "the frozen feature must not move"
        assert drop_b == 0.7
        # and nothing the old axis carried was rewritten, removed or reordered
        for row_before, row_after in zip(before["rows"], after["rows"], strict=True):
            assert {k: v for k, v in row_after.items() if k != astrorun.FULL_VECTOR_KEY} == row_before
            assert list(row_after)[: len(row_before)] == list(row_before), "appended, not reordered"
        added = set(after) - set(before)
        assert added == {"track_names", "tracks_sha256", "full_vector_note"}

    def test_the_DISK_ESTIMATE_is_stated_and_checked_against_the_floor(self):
        """(d): measured from the format, not assumed, and a refusal rather than a warning."""
        est = astrorun.full_vector_disk_estimate(
            1232,
            astrorun.OBSERVED_ROWS_PER_ANSWER,
            astrorun.OBSERVED_TRACKS_PER_ANSWER,
            free_bytes=50 * 1024**3,
        )
        tracks = astrorun.OBSERVED_TRACKS_PER_ANSWER
        assert est["tracks"] == tracks == 667, "the larger of the two figures"
        assert est["rows_per_element"] == 29
        # base64 of float64: 8 bytes a value, 4 characters a 3 bytes
        assert est["bytes_per_row"] == math.ceil(8 * tracks / 3) * 4 + 40
        assert 200 < est["mb_total"] < 300, est["mb_total"]
        # and on the smaller figure it is smaller, which is why the larger is the one used
        note = astrorun.full_vector_disk_estimate(
            1232, astrorun.ADAPTER_NOTE_ROWS, astrorun.ADAPTER_NOTE_TRACKS, free_bytes=50 * 1024**3
        )
        assert note["mb_total"] < est["mb_total"]
        assert est["fits_above_the_floor"] is True
        with pytest.raises(astrorun.VectorRefusedError) as exc:
            astrorun.check_disk_for_full_vectors(1232, 29, 371, free_bytes=10 * 1024**3 + 1)
        assert "would leave less than the 10 GB floor" in str(exc.value)
        assert "refusal and not a warning" in str(exc.value)

    def test_the_measured_row_size_matches_the_estimate(self):
        """The estimate's own arithmetic, checked against a row actually written."""
        axis = astrorun.with_full_track_vectors(
            self.axis(rows=1, tracks=371), self.Matrix(1, 371), [f"t{i}" for i in range(371)], "a" * 64
        )
        written = len(json.dumps(axis["rows"][0][astrorun.FULL_VECTOR_KEY]))
        est = astrorun.full_vector_disk_estimate(1, 1, 371, free_bytes=50 * 1024**3)
        assert abs(written - est["bytes_per_row"]) <= 40, (written, est["bytes_per_row"])

    @needs_track_metadata
    def test_the_REAL_recording_path_keeps_the_WHOLE_vector(self):
        """(a) The real recorded_axis, run and read back: every value, named, equal.

        A source check cannot say this and a stub cannot either. The response is built with the
        library's own AnnData over the REAL RNA-seq track metadata, so the axis comes from the file
        and the library rather than from what our code asks for.
        """
        import anndata
        import numpy
        import pandas

        from genomeos.predict import alphagenome_adapter as adapter

        meta_csv = Path(astrorun.ROOT_FOR_BLOBS) / "data/cache/entex/alphagenome_track_metadata_copy.csv"
        # The marker guarantees the store is here when this runs, and skips BY NAME when it is not.
        # A bare `assert exists()` was wrong: data/cache is machine-local by the data boundary, so the
        # failure it produced was about the checkout and not about the code -- a verdict determined by
        # the environment, which is the defect class this project spent 2026-10-02 removing.
        with open(meta_csv) as fh:
            meta = [r for r in csv.DictReader(fh) if r["output"] == "rna_seq"]
        assert len(meta) > 300, "the real metadata, not a handful of rows"
        var = pandas.DataFrame(meta, index=[f"{r['name']}|{i}" for i, r in enumerate(meta)])
        obs = pandas.DataFrame(
            {
                "gene_id": ["ENSG1", "ENSG2", "ENSG3"],
                "gene_name": ["AAA", "BBB", "CCC"],
                "strand": ["+", "-", "+"],
            },
            index=["ENSG1", "ENSG2", "ENSG3"],
        )
        matrix = numpy.random.default_rng(7).normal(size=(3, len(meta)))
        matrix[1, 0] = 0.0  # an exact zero: a thresholded writer would drop it and look fine
        adata = anndata.AnnData(X=matrix, obs=obs, var=var)
        tissues = adapter.tissue_names(adata)
        assert len(tissues) == len(meta)

        axis = adapter.recorded_axis(0, adata, tissues, 0.0)
        assert axis["tracks_total"] == len(meta)
        assert axis["track_names"] == tissues
        assert len(axis["tracks_sha256"]) == 64
        for gi in range(3):
            pairs = astrorun.named_track_values(axis, gi)
            assert len(pairs) == len(meta), f"one value per track for row {gi}"
            assert [n for n, _ in pairs] == tissues
            assert [v for _, v in pairs] == [float(matrix[gi, ti]) for ti in range(len(meta))], (
                "equal by value, not close: this is the paid answer"
            )
        assert astrorun.named_track_values(axis, 1)[0][1] == 0.0, "the exact zero survived"
        got = astrorun.check_full_vectors_in_answer({"model": {"gene_axis_outputs": [axis]}})
        assert got["gene_rows_with_full_vectors"] == 3

    def test_the_REAL_adapter_PASSES_the_send_time_check(self):
        """The positive control for item (g)'s gate, run against the module that will do the buying."""
        out = astrorun.check_adapter_writes_full_vectors()
        assert out["recorded_axis_keeps_the_full_vector"] is True
        assert out["how"].startswith("the recording path was RUN")

    def test_PLANTED_a_recorder_that_only_MENTIONS_the_key_in_its_DOCSTRING(self, tmp_path, monkeypatch):
        """(b) Both verdicts side by side: the source check PASSES it, the behavioural check FAILS it.

        This was the defect in the first version of the gate. ast.dump carries docstrings and every
        string constant, so a recorder that merely NAMES the key satisfied it, and it refused
        correctly only while the adapter was untouched. The day it passed would have been the day it
        proved nothing -- and that day was the day the adapter changed.
        """
        import ast

        mod = tmp_path / "recorder_that_only_talks_about_it.py"
        mod.write_text(
            "def recorded_axis(output, adata, tissues, threshold=0.0):\n"
            '    "Keeps each row track_vector, honestly, see astrorun.FULL_VECTOR_KEY."\n'
            "    return {'output': output, 'rows': [{'row': 0, 'gene_name': 'AAA'}],\n"
            "            'tracks_total': len(tissues)}\n"
        )
        tree = ast.parse(mod.read_text())
        fn = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == "recorded_axis")
        source_verdict = astrorun.FULL_VECTOR_KEY in ast.dump(fn)
        assert source_verdict is True, "THE SOURCE CHECK PASSES IT: existence by name"

        monkeypatch.syspath_prepend(str(tmp_path))
        with pytest.raises(astrorun.VectorRefusedError) as exc:
            astrorun.check_adapter_writes_full_vectors("recorder_that_only_talks_about_it")
        assert "does not write the whole track vector" in str(exc.value)
        assert "existence-by-name" in astrorun.EXISTENCE_BY_NAME_IS_NOT_A_CHECK

    def test_PLANTED_a_recorder_that_THRESHOLDS_fails_on_the_VALUES(self, tmp_path, monkeypatch):
        """A vector that is present but not what arrived is the same loss, one step quieter."""
        mod = tmp_path / "recorder_that_filters.py"
        mod.write_text(
            "from genomeos.attribution import astrorun\n\n\n"
            "def recorded_axis(output, adata, tissues, threshold=0.0):\n"
            "    rows = []\n"
            "    for gi in range(adata.X.shape[0]):\n"
            "        vals = [float(adata.X[gi, ti]) for ti in range(len(tissues))]\n"
            "        vals = [v if abs(v) > 0.1 else 0.0 for v in vals]\n"
            "        rows.append({'row': gi, astrorun.FULL_VECTOR_KEY: astrorun.full_vector_entry(vals)})\n"
            "    return {'output': output, 'rows': rows, 'tracks_total': len(tissues),\n"
            "            'track_names': list(tissues), 'tracks_sha256': 'x' * 64}\n"
        )
        monkeypatch.syspath_prepend(str(tmp_path))
        with pytest.raises(astrorun.VectorRefusedError) as exc:
            astrorun.check_adapter_writes_full_vectors("recorder_that_filters")
        said = str(exc.value)
        assert "came back as" in said
        assert "thresholded, rounded or partial vector is not what was paid for" in said

    def test_the_probe_carries_a_zero_and_a_negative_on_purpose(self):
        """A probe whose values all survive a filter would pass a filtering writer."""
        flat = [v for row in astrorun.PROBE_VALUES for v in row]
        assert 0.0 in flat, "an exact zero, which a threshold drops"
        assert any(v < 0 for v in flat), "and a negative, which a sign error flips"

    def test_PLANTED_item_g_is_WIRED_into_may_send_and_a_REGRESSED_recorder_stops_it(
        self, tmp_path, monkeypatch
    ):
        """A mechanism not wired to what it protects is indistinguishable from one that is absent.

        The subject is the wiring. The real adapter passes now, so the plant is the module may_send
        looks at: pointed at a recorder that only mentions the key, the REAL check runs and the send
        stops. Nothing about the check is stubbed here.
        """
        regressed = tmp_path / "regressed_recorder.py"
        regressed.write_text(
            "def recorded_axis(output, adata, tissues, threshold=0.0):\n"
            '    "I keep the track_vector, see astrorun.FULL_VECTOR_KEY."\n'
            "    return {'output': output, 'rows': [{'row': 0}], 'tracks_total': len(tissues)}\n"
        )
        monkeypatch.syspath_prepend(str(tmp_path))
        monkeypatch.setattr(astrorun, "ADAPTER_MODULE", "regressed_recorder")

        act = tmp_path / "activity.json"
        act.write_text(json.dumps({"rule": {"producer": dict(astrorun.AMENDMENT_2_RULE_FINGERPRINT)}}))
        reg = tmp_path / "reg.json"
        reg.write_text("{}")
        (tmp_path / astrorun.LEDGER_RUN1.parent).mkdir(parents=True, exist_ok=True)
        plan = [
            {"chrom": "chr21", "element": f"e{i}", "start": i, "end": i + 1, "serves_genes": ["G"]}
            for i in range(2)
        ]
        monkeypatch.setattr(astrorun, "ASTROREG2_CAP", 2)
        monkeypatch.setattr(astrorun, "ASTROREG2_AUTHORISATION", astrorun.ASTROREG2_AUTHORISATION_AS_RELAYED)
        monkeypatch.setattr(
            astrorun,
            "ASTROREG2_AUTHORISATIONS",
            {1: {"words": "<words>", "requests": 2, "consumed": False}},
        )
        monkeypatch.setattr(astrorun, "check_signoff_closure", lambda *a, **k: {})
        monkeypatch.setattr(astrorun, "check_nothing_in_the_plan_is_already_cached", lambda *a, **k: None)
        monkeypatch.setattr(astrorun, "check_adapter_v2", lambda *a, **k: {"stubbed": True})
        with pytest.raises(astrorun.VectorRefusedError, match=r"item \(g\)"):
            astrorun.may_send(
                activity_result=act,
                registration=reg,
                run_id=1,
                ledger_root=tmp_path,
                signoff='the supervisor wrote "dry run reviewed" at 20:05',
                plan=plan,
                reviewed_digest=astrorun.plan_digest(plan),
                committed=lambda p: True,
                free_bytes=11 * 1024**3,
            )

    def test_may_send_also_checks_the_RECORDER_and_the_DISK_before_a_paid_run(self):
        """Both of item (g)'s gates are called from may_send, before anything is bought."""
        import ast

        tree = ast.parse(Path(astrorun.__file__).read_text())
        fn = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == "may_send")
        called = {c.func.id for c in ast.walk(fn) if isinstance(c, ast.Call) and isinstance(c.func, ast.Name)}
        assert "check_adapter_writes_full_vectors" in called, "item (g)'s recording check"
        assert "check_disk_for_full_vectors" in called, "item (g)'s disk check"

    @needs_track_metadata
    def test_the_two_track_counts_are_both_recorded_and_the_LARGER_is_estimated_on(self):
        """Neither figure is measured from a paid response, so the estimate takes the larger.

        371 x 29 is adapter v2's own note. 667 is the rna_seq row count of the real metadata, checked
        here against the file. No cached answer on this machine records a gene axis, so there is
        nothing to measure either against -- and underestimating the disk a paid run needs is the
        direction that stops a run part way through.
        """
        meta_csv = Path(astrorun.ROOT_FOR_BLOBS) / "data/cache/entex/alphagenome_track_metadata_copy.csv"
        assert astrorun.ADAPTER_NOTE_TRACKS == 371
        assert astrorun.ADAPTER_NOTE_ROWS == 29
        assert (
            max(astrorun.ADAPTER_NOTE_TRACKS, astrorun.RNA_SEQ_TRACKS_IN_THE_METADATA)
            == astrorun.OBSERVED_TRACKS_PER_ANSWER
        )
        with open(meta_csv) as fh:
            rna = [r for r in csv.DictReader(fh) if r["output"] == "rna_seq"]
        assert len(rna) == astrorun.RNA_SEQ_TRACKS_IN_THE_METADATA, "read from the file, not remembered"
        assert astrorun.ADAPTER_NOTE_TRACKS == 371
        assert astrorun.ADAPTER_NOTE_ROWS == 29
        assert max(371, len(rna)) == astrorun.OBSERVED_TRACKS_PER_ANSWER
        assert "Neither" in astrorun.OBSERVED_SHAPE_SOURCE or "neither" in astrorun.OBSERVED_SHAPE_SOURCE

    def test_the_two_figures_straddle_the_floor_so_the_pair_is_a_real_contrast(self):
        """Without this, both plants could sit on one side of the floor and prove nothing."""
        need = astrorun.full_vector_disk_estimate(
            1232,
            astrorun.OBSERVED_ROWS_PER_ANSWER,
            astrorun.OBSERVED_TRACKS_PER_ANSWER,
            free_bytes=11 * 1024**3,
        )
        assert need["fits_above_the_floor"] is True
        tight = astrorun.full_vector_disk_estimate(
            1232,
            astrorun.OBSERVED_ROWS_PER_ANSWER,
            astrorun.OBSERVED_TRACKS_PER_ANSWER,
            free_bytes=int(9.46 * 1024**3),
        )
        assert tight["fits_above_the_floor"] is False
        assert tight["bytes_total"] == need["bytes_total"], "same need, different machine"

    def test_the_LIVE_reading_is_exercised_but_cannot_RED_the_suite(self):
        """One test keeps the live path honest: it returns a number. It does NOT judge the number.

        Asserting the live figure is above or below the floor is what made the suite's verdict depend on
        free disk. So this checks the reading happens and is well formed, and says nothing about how
        much space the machine has.
        """
        live = astrorun.full_vector_disk_estimate(
            1232, astrorun.OBSERVED_ROWS_PER_ANSWER, astrorun.OBSERVED_TRACKS_PER_ANSWER
        )
        assert isinstance(live["free_bytes"], int)
        assert live["free_bytes"] > 0, "the volume was read"
        assert isinstance(live["fits_above_the_floor"], bool), "a verdict was formed, not asserted here"
        assert live["bytes_total"] == 1232 * live["bytes_per_answer"], "the need is machine-independent"
        assert "a PARAMETER" in astrorun.THE_FIGURE_IS_INJECTED_IN_TESTS_AND_READ_ONLY_IN_A_REAL_RUN

    def test_may_send_passes_the_INJECTED_figure_through_rather_than_ignoring_it(self):
        """Planted against a may_send that called the clause with no figure and read the volume anyway."""
        import ast

        tree = ast.parse(Path(astrorun.__file__).read_text())
        fn = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == "may_send")
        assert "free_bytes" in {a.arg for a in fn.args.kwonlyargs}, "may_send takes the figure"
        call = next(
            c
            for c in ast.walk(fn)
            if isinstance(c, ast.Call)
            and isinstance(c.func, ast.Name)
            and c.func.id == "check_disk_for_full_vectors"
        )
        names = [a.id for a in call.args if isinstance(a, ast.Name)]
        assert "free_bytes" in names, (
            "may_send must hand the clause the injected figure; without it the clause reads the volume "
            "and every test of every other clause becomes a test of the laptop"
        )

    def test_item_g_is_labelled_the_supervisors_and_not_alberts(self):
        """Like item (f): a refusal must not attribute to him a condition he did not state."""
        assert "not a clause of Albert's approval" in astrorun.ITEM_G_IS_A_SUPERVISOR_REQUIREMENT
        assert "SUPERVISOR" in astrorun.ITEM_G_IS_A_SUPERVISOR_REQUIREMENT.upper()

    def test_the_vector_is_stored_LOSSLESSLY(self):
        """float64, because a lossy store of data that cannot be re-fetched is the same mistake."""
        values = [0.1234567890123456, -9.87654321098765e-8, 0.0, -0.0, 1e300]
        back = astrorun.unpack_track_vector(astrorun.pack_track_vector(values), expected=5)
        assert back == values, "equal, not close"
        assert astrorun.FULL_VECTOR_DTYPE == "<f8"


class TestItemFAdapterV2:
    """Sign-off item (f), planted. It is the SUPERVISOR's requirement, not a clause of Albert's."""

    def test_adapter_v2_is_now_committed_so_item_f_is_SATISFIED(self):
        """It refused until adapter v2 landed. It passes now, against the real module.

        This check was written before adapter v2 existed and guessed the names it would expose. The
        guesses were wrong, the check refused, and the fix was to READ the committed adapter and
        reconcile to it -- not to relax the check. ADAPTER_V2_MUST_EXPOSE now names what is actually
        there, and `gene_id` is required to be in the gene-axis column set rather than merely implied.
        """
        out = astrorun.check_adapter_v2()
        assert out["module"] == "genomeos.predict.alphagenome_adapter"
        assert "GENE_AXIS_COLUMNS_090" in out["exposes"]

    def test_the_gene_id_field_is_required_and_not_merely_implied(self):
        """Without gene_id a purchased answer cannot be tied to a gene identity at all."""
        import genomeos.predict.alphagenome_adapter as real

        assert astrorun.ADAPTER_V2_MUST_RECORD_GENE_ID in real.GENE_AXIS_COLUMNS_090

    def test_it_is_not_attributed_to_albert(self):
        assert "not a clause of Albert's approval" in astrorun.ADAPTER_V2_IS_A_SUPERVISOR_REQUIREMENT
        assert "adapter" not in " ".join(astrorun.ASTROREG2_CLAUSES.values()).lower()

    def test_a_missing_module_refuses_rather_than_raising_an_import_error(self):
        with pytest.raises(astrorun.SendRefusedError, match="not importable"):
            astrorun.check_adapter_v2("genomeos.predict.no_such_adapter_module")

    def test_a_module_exposing_everything_passes(self, monkeypatch):
        """The positive control: a capability check that can never pass would prove nothing."""
        import types

        mod = types.ModuleType("fake_adapter_v2")
        for name in astrorun.ADAPTER_V2_MUST_EXPOSE:
            setattr(mod, name, True)
        mod.GENE_AXIS_COLUMNS_090 = (astrorun.ADAPTER_V2_MUST_RECORD_GENE_ID, "gene_name")
        monkeypatch.setitem(sys.modules, "fake_adapter_v2", mod)
        out = astrorun.check_adapter_v2("fake_adapter_v2")
        assert out["exposes"] == list(astrorun.ADAPTER_V2_MUST_EXPOSE)

    def test_one_missing_capability_is_enough_to_refuse(self, monkeypatch):
        import types

        mod = types.ModuleType("half_adapter")
        for name in astrorun.ADAPTER_V2_MUST_EXPOSE[:-1]:
            setattr(mod, name, True)
        mod.GENE_AXIS_COLUMNS_090 = (astrorun.ADAPTER_V2_MUST_RECORD_GENE_ID,)
        monkeypatch.setitem(sys.modules, "half_adapter", mod)
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.check_adapter_v2("half_adapter")
        assert astrorun.ADAPTER_V2_MUST_EXPOSE[-1] in str(exc.value)

    def test_a_version_constant_alone_does_not_satisfy_it(self, monkeypatch):
        """Checked by capability, because a version string can be set without the behaviour."""
        import types

        mod = types.ModuleType("version_only")
        mod.GENE_AXIS_SCHEMA = 1
        monkeypatch.setitem(sys.modules, "version_only", mod)
        with pytest.raises(astrorun.SendRefusedError):
            astrorun.check_adapter_v2("version_only")

    def test_the_merge_wording_is_carried_verbatim(self):
        assert astrorun.A_MERGE_IS_A_LOST_DISTINCTION == ("a merge is a LOST DISTINCTION, NOT A WRONG NUMBER")

    def test_PLANTED_a_regressed_adapter_refuses_the_send_with_everything_else_satisfied(
        self, tmp_path, monkeypatch
    ):
        """PLANTED with Albert's authorisation recorded AND every clause of his satisfied.

        Adapter v2 is committed, so item (f) passes against the real module. What must still hold is
        that a FUTURE adapter losing a capability stops the send, so the check is pointed at a module
        missing one and the refusal is shown.
        """
        import types

        half = types.ModuleType("regressed_adapter")
        for name in astrorun.ADAPTER_V2_MUST_EXPOSE[:-1]:
            setattr(half, name, True)
        half.GENE_AXIS_COLUMNS_090 = (astrorun.ADAPTER_V2_MUST_RECORD_GENE_ID,)
        monkeypatch.setitem(sys.modules, "regressed_adapter", half)
        monkeypatch.setattr(astrorun, "ADAPTER_MODULE", "regressed_adapter")
        monkeypatch.setattr(astrorun, "check_signoff_closure", lambda *a, **k: {})
        monkeypatch.setattr(astrorun, "ASTROREG2_AUTHORISATION", astrorun.ASTROREG2_AUTHORISATION_AS_RELAYED)
        plan = [
            {
                "chrom": "chr1",
                "element": f"E{i}",
                "start": i,
                "end": i + 1,
                "serves_genes": ["G"],
                "serves_labels": ["negative"],
            }
            for i in range(astrorun.ASTROREG2_CAP)
        ]
        act = tmp_path / "activity.json"
        act.write_text(json.dumps({"rule": {"producer": dict(astrorun.AMENDMENT_2_RULE_FINGERPRINT)}}))
        reg = tmp_path / "reg.json"
        reg.write_text("{}")
        grant_run(monkeypatch)
        with pytest.raises(astrorun.SendRefusedError, match=r"item \(f\)"):
            astrorun.may_send(
                activity_result=act,
                registration=reg,
                run_id=1,
                ledger_root=tmp_path,
                signoff='the supervisor wrote "dry run reviewed" at 2026-10-02T13:00:00',
                plan=plan,
                reviewed_digest=astrorun.plan_digest(plan),
                committed=lambda p: True,
            )
