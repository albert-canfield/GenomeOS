"""The money guard and the activity-input gate of the authorised AstroREG run.

The cap tests are written before any request can be sent. They do not check that a warning is
printed at 1,322; they check that the 1,323rd request cannot be obtained at all, that the count
stops at the cap, and that the ledger on disk holds exactly one line per charge.
"""

from __future__ import annotations

import json
import sys
import threading
from pathlib import Path

import pytest

from genomeos.attribution import astrorun


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
        ns: dict = {}
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

    def test_no_astroreg2_approval_is_recorded_yet(self):
        assert astrorun.ASTROREG2_AUTHORISATION is None

    def test_asking_for_the_astroreg2_authorisation_refuses_while_it_is_absent(self):
        with pytest.raises(astrorun.NoAuthorisationError, match="does not carry"):
            astrorun.astroreg2_authorisation()

    def test_no_budget_can_be_obtained_for_astroreg2_without_its_own_approval(self, tmp_path):
        """The refusal comes BEFORE the budget exists, so there is nothing to send with."""
        with pytest.raises(astrorun.NoAuthorisationError):
            astrorun.astroreg2_budget(tmp_path / "l.jsonl")
        assert not (tmp_path / "l.jsonl").exists(), "a refusal creates no ledger"

    def test_the_astroreg2_scope_cannot_reach_the_original_1322_cap(self, tmp_path, monkeypatch):
        """With an approval recorded, the cap is 1,232: request 1,233 is refused, 1,322 unreachable."""
        monkeypatch.setattr(astrorun, "ASTROREG2_AUTHORISATION", "<verbatim words would go here>")
        b = astrorun.astroreg2_budget(tmp_path / "l.jsonl")
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
        """Albert's words recorded. Every other clause must still be checked."""
        monkeypatch.setattr(astrorun, "ASTROREG2_AUTHORISATION", astrorun.ASTROREG2_AUTHORISATION_AS_RELAYED)

    @pytest.fixture
    def good(self, tmp_path):
        act = tmp_path / "activity.json"
        act.write_text(json.dumps({"rule": {"producer": dict(astrorun.AMENDMENT_2_RULE_FINGERPRINT)}}))
        reg = tmp_path / "astroreg2_registration.json"
        reg.write_text("{}")
        return {
            "activity_result": act,
            "registration": reg,
            "ledger": tmp_path / "l.jsonl",
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
        out = astrorun.may_send(**good)
        assert out["may_send"] is True
        assert out["cap"] == 1232
        assert sorted(out["clauses_checked"]) == ["activity", "logged", "one_run", "scope", "signoff"]

    def test_the_slot_is_empty_and_the_relayed_text_is_kept_apart(self):
        assert astrorun.ASTROREG2_AUTHORISATION is None
        assert "I approve 1,232" in astrorun.ASTROREG2_AUTHORISATION_AS_RELAYED
        assert "not that person's approval" in astrorun.WHY_THE_RELAYED_TEXT_IS_NOT_THE_APPROVAL

    def test_without_the_authorisation_recorded_nothing_sends(self, good):
        with pytest.raises(astrorun.NoAuthorisationError):
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
        astrorun.mark_run_complete(good["ledger"], requests=1232)
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


class TestItemFAdapterV2:
    """Sign-off item (f), planted. It is the SUPERVISOR's requirement, not a clause of Albert's."""

    def test_it_refuses_today_because_the_adapter_is_still_v1(self):
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.check_adapter_v2()
        msg = str(exc.value)
        assert "item (f)" in msg
        assert "ADAPTER_V2" in msg, "the refusal names what it looked for"
        assert "PAID FOR" in msg

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
        monkeypatch.setitem(sys.modules, "fake_adapter_v2", mod)
        out = astrorun.check_adapter_v2("fake_adapter_v2")
        assert out["exposes"] == list(astrorun.ADAPTER_V2_MUST_EXPOSE)

    def test_one_missing_capability_is_enough_to_refuse(self, monkeypatch):
        import types

        mod = types.ModuleType("half_adapter")
        for name in astrorun.ADAPTER_V2_MUST_EXPOSE[:-1]:
            setattr(mod, name, True)
        monkeypatch.setitem(sys.modules, "half_adapter", mod)
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.check_adapter_v2("half_adapter")
        assert astrorun.ADAPTER_V2_MUST_EXPOSE[-1] in str(exc.value)

    def test_a_version_constant_alone_does_not_satisfy_it(self, monkeypatch):
        """Checked by capability, because a version string can be set without the behaviour."""
        import types

        mod = types.ModuleType("version_only")
        mod.ADAPTER_V2 = True
        monkeypatch.setitem(sys.modules, "version_only", mod)
        with pytest.raises(astrorun.SendRefusedError):
            astrorun.check_adapter_v2("version_only")

    def test_the_merge_wording_is_carried_verbatim(self):
        assert astrorun.A_MERGE_IS_A_LOST_DISTINCTION == ("a merge is a LOST DISTINCTION, NOT A WRONG NUMBER")

    def test_the_send_path_refuses_on_item_f_even_with_everything_else_satisfied(self, tmp_path, monkeypatch):
        """PLANTED with Albert's authorisation recorded AND every clause of his satisfied."""
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
        with pytest.raises(astrorun.SendRefusedError, match=r"item \(f\)"):
            astrorun.may_send(
                activity_result=act,
                registration=reg,
                ledger=tmp_path / "l.jsonl",
                signoff='the supervisor wrote "dry run reviewed" at 2026-10-02T13:00:00',
                plan=plan,
                reviewed_digest=astrorun.plan_digest(plan),
                committed=lambda p: True,
            )
