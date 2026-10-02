"""The sender: one path, ledgered before each send, and a pilot checkpoint that bounds a silent waste.

The two tests that carry the weight are PLANTED, because a test that merely passes against the current
code proves nothing here. The whole risk is a FUTURE sender taking a shortcut, so:

  - the client is monkeypatched so that ONLY the recording path can satisfy it, and a sender that calls
    the client directly is planted and shown failing;
  - a stub response without gene names is planted in the first ten answers and the pilot stop is shown.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

from genomeos.attribution import astrorun

ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location("astroreg2_send", ROOT / "scripts/astroreg2_send.py")
sender = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = sender
_SPEC.loader.exec_module(sender)


def row(i: int) -> dict:
    return {
        "chrom": "chr1",
        "element": f"EH38E{i:07d}",
        "start": i * 1000,
        "end": i * 1000 + 300,
        "serves_genes": ["AAA"],
        "serves_labels": ["negative"],
    }


def good_answer(element_id: str) -> dict:
    return {"id": element_id, "genes": [{"gene": "AAA", "mean_log2fc": -0.4, "n_tracks": 3}]}


def nameless_answer(element_id: str) -> dict:
    """The latent fault: a gene proto without `name` yields ZERO effects, silently."""
    return {"id": element_id, "genes": []}


class RecordingPathOnly:
    """A client double that only the recording path can satisfy.

    `score_element` is the recording path and is the only caller allowed to reach `predict`. Anything
    calling the client directly raises, so a sender that takes the shortcut cannot pass.
    """

    def __init__(self):
        self.calls_through_recording_path = 0
        self.direct_calls = 0

    def predict(self, *a, **k):
        self.direct_calls += 1
        raise AssertionError(
            "the AlphaGenome client was called DIRECTLY. Every request must go through "
            "enhancer_target.score_element, the recording path that writes the gene-axis record and "
            "each named cell's value multiset. A direct call buys the answer and discards what was paid "
            "for"
        )

    def through_recording_path(self, element_id: str, answer=good_answer):
        self.calls_through_recording_path += 1
        return answer(element_id)


class TestOnlyTheRecordingPathMaySend:
    def test_the_sender_sends_through_the_recording_path(self, tmp_path):
        """The positive control: the shape the sender actually uses must work."""
        client = RecordingPathOnly()
        budget = astrorun.RequestBudget(tmp_path / "l.jsonl", cap=20)
        out = sender.send(
            [row(i) for i in range(12)],
            budget,
            lambda **kw: client.through_recording_path(kw["element_id"]),
        )
        assert out["requests_sent"] == 12
        assert client.calls_through_recording_path == 12
        assert client.direct_calls == 0
        assert out["pilot"]["passed"] is True

    def test_PLANTED_a_sender_that_calls_the_client_directly_FAILS(self, tmp_path):
        """By planting the shortcut: the client double makes the bypass unsatisfiable."""
        client = RecordingPathOnly()
        budget = astrorun.RequestBudget(tmp_path / "l.jsonl", cap=20)

        def shortcut(**kw):
            return client.predict(kw["chrom"], kw["start"])  # the direct call being planted

        with pytest.raises(AssertionError, match="called DIRECTLY"):
            sender.send([row(i) for i in range(12)], budget, shortcut)
        assert client.direct_calls == 1
        assert client.calls_through_recording_path == 0

    def test_the_sender_imports_no_client_of_its_own(self):
        """The structural half: the sender has no route to the client except the recording path."""
        import ast

        tree = ast.parse((ROOT / "scripts/astroreg2_send.py").read_text())
        imported = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.update(a.name for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.add(node.module)
                imported.update(f"{node.module}.{a.name}" for a in node.names)
        # The sender DOES import AlphaGenomeAdapter, because a scorer has to be constructed from
        # something, and that scorer is then passed to score_element which does the recording. What must
        # not happen is the sender calling the client itself -- which the planted direct-call test above
        # is the real guarantee of. These assertions cover the structural half: no import of the vendor
        # package, and no import of the client factory.
        assert not any(r == "alphagenome" or r.startswith("alphagenome.") for r in imported), (
            "the sender must not import the vendor package; it reaches the model through the adapter "
            "and the recording path"
        )
        assert "genomeos.predict.alphagenome_adapter.create_client" not in imported, (
            "the client factory is not the sender's to call"
        )
        assert "genomeos.predict.enhancer_target" in imported

    def test_a_ledger_line_is_written_before_each_answer(self, tmp_path):
        led = tmp_path / "l.jsonl"
        client = RecordingPathOnly()
        budget = astrorun.RequestBudget(led, cap=20)
        sender.send(
            [row(i) for i in range(3)],
            budget,
            lambda **kw: client.through_recording_path(kw["element_id"]),
            pilot_requests=99,
        )
        events = [json.loads(x)["event"] for x in led.read_text().splitlines()]
        assert events[:2] == ["request", "answer"], "the charge is logged before its answer"
        assert events.count("request") == 3
        assert budget.charged() == 3

    def test_the_run_is_marked_complete_at_the_end(self, tmp_path):
        led = tmp_path / "l.jsonl"
        client = RecordingPathOnly()
        sender.send(
            [row(i) for i in range(3)],
            astrorun.RequestBudget(led, cap=20),
            lambda **kw: client.through_recording_path(kw["element_id"]),
            pilot_requests=99,
        )
        assert astrorun.run_already_completed(led) is True


class TestThePilotCheckpoint:
    def answers(self, bad_indices):
        def score(**kw):
            i = int(kw["element_id"].replace("EH38E", ""))
            return nameless_answer(kw["element_id"]) if i in bad_indices else good_answer(kw["element_id"])

        return score

    def test_PLANTED_a_nameless_response_in_the_first_ten_stops_the_run(self, tmp_path):
        """At most ten requests spent, the failures recorded, and the run terminal."""
        led = tmp_path / "l.jsonl"
        budget = astrorun.RequestBudget(led, cap=astrorun.ASTROREG2_CAP)
        out = sender.send([row(i) for i in range(50)], budget, self.answers({4}))
        assert out["requests_sent"] == 10, "it stops AT the checkpoint, not after the whole list"
        assert out["pilot"]["passed"] is False
        assert out["pilot"]["stopped"]["at_request"] == 10
        assert [f["id"] for f in out["pilot"]["stopped"]["failing"]] == ["EH38E0000004"]
        assert astrorun.run_already_completed(led) is True, "a pilot stop is terminal"

    def test_the_stop_says_it_consumes_albert_s_one_run_rather_than_exhausting_a_budget(self, tmp_path):
        out = sender.send(
            [row(i) for i in range(50)],
            astrorun.RequestBudget(tmp_path / "l.jsonl", cap=astrorun.ASTROREG2_CAP),
            self.answers({0}),
        )
        said = out["pilot"]["stopped"]["consequence"]
        assert "CONSUMES Albert's one run" in said
        assert "NOT a budget that ran out" in said
        assert "NEW authorisation in his own words" in said
        assert said == sender.PILOT_STOP_CONSUMES_THE_ONE_RUN

    def test_a_clean_first_ten_lets_the_run_continue(self, tmp_path):
        out = sender.send(
            [row(i) for i in range(14)],
            astrorun.RequestBudget(tmp_path / "l.jsonl", cap=50),
            self.answers(set()),
        )
        assert out["requests_sent"] == 14
        assert out["pilot"]["passed"] is True

    def test_a_failure_AFTER_the_checkpoint_does_not_stop_the_run_but_is_recorded(self, tmp_path):
        """The checkpoint is a gate at ten, not a filter over the whole list."""
        out = sender.send(
            [row(i) for i in range(14)],
            astrorun.RequestBudget(tmp_path / "l.jsonl", cap=50),
            self.answers({12}),
        )
        assert out["requests_sent"] == 14
        assert out["pilot"]["passed"] is True
        assert [c["id"] for c in out["unusable_answers"]] == ["EH38E0000012"]

    def test_an_answer_with_a_blank_gene_name_is_unusable(self):
        assert sender.answer_is_usable({"id": "x", "genes": [{"gene": "  "}]})["usable"] is False
        assert sender.answer_is_usable({"id": "x", "genes": []})["usable"] is False
        assert sender.answer_is_usable({"id": "x", "genes": [{"gene": "AAA"}]})["usable"] is True

    def test_the_pilot_rule_names_the_silent_fault_it_guards(self):
        assert "indistinguishable from 'no gene moved'" in sender.PILOT_RULE
        assert sender.PILOT_REQUESTS == 10


class TestTheSenderRefusesToday:
    def test_the_reviewed_digest_is_the_committed_one(self):
        assert sender.REVIEWED_DIGEST == ("b96148cab7174cf10407aa07a2d26edf1c1972df5db2c00de3e682ccacf2b519")
        assert sender.REVIEWED_AT == "cc5b097"

    def test_a_tampered_plan_is_refused_before_anything_is_sent(self, tmp_path, monkeypatch):
        plan = [row(i) for i in range(3)]
        p = tmp_path / "plan.json"
        p.write_text(json.dumps({"requests": plan}))
        monkeypatch.setattr(sender, "PLAN", p)
        monkeypatch.setattr(sender, "REVIEWED_DIGEST", astrorun.plan_digest(plan))
        assert len(sender.reviewed_plan()) == 3
        plan[1]["element"] = "tampered"
        p.write_text(json.dumps({"requests": plan}))
        with pytest.raises(ValueError, match="not the reviewed"):
            sender.reviewed_plan()

    def test_main_refuses_end_to_end_when_no_authorisation_is_recorded(self, tmp_path, monkeypatch, capsys):
        """End to end: the entry point refuses, prints why, and sends nothing.

        Converted from asserting the slots were empty -- a fact about the day -- to patching them empty
        and asserting the behaviour, so it holds now that both are filled and would catch a regression
        that made an empty slot permissive.
        """
        monkeypatch.setattr(astrorun, "ASTROREG2_AUTHORISATION", None)
        monkeypatch.setattr(astrorun, "SUPERVISOR_SIGNOFF", None)
        plan = [row(i) for i in range(3)]
        p = tmp_path / "plan.json"
        p.write_text(json.dumps({"requests": plan}))
        monkeypatch.setattr(sender, "PLAN", p)
        monkeypatch.setattr(sender, "REVIEWED_DIGEST", astrorun.plan_digest(plan))
        monkeypatch.setattr(sys, "argv", ["astroreg2_send.py", "--send"])
        assert sender.main() == 2
        said = capsys.readouterr().out
        assert "REFUSED" in said
        assert "no AstroREG-2 approval is recorded" in said

    def test_the_budget_comes_only_from_astroreg2_budget(self):
        src = (ROOT / "scripts/astroreg2_send.py").read_text()
        assert "astrorun.astroreg2_budget(" in src
        assert "RequestBudget(" not in src, (
            "the sender must not construct a budget of its own: astroreg2_budget is what refuses "
            "without Albert's authorisation"
        )


class TestACrashRestartCannotRePay:
    """The blocker: a process that dies at request k must not pay for elements 1..k again.

    The machinery existed (charged_elements, requests_not_yet_charged) and the intent was written down,
    and NEITHER was wired. A mechanism that exists but is not wired is indistinguishable from one that
    does not exist, which is why these are planted rather than argued.
    """

    def partial_ledger(self, path, plan, k):
        """A ledger holding k charges and NO completion: exactly what a crash leaves behind."""
        budget = astrorun.RequestBudget(path, cap=astrorun.ASTROREG2_CAP)
        for r in plan[:k]:
            budget.take(chrom=r["chrom"], element=r["element"], start=r["start"], end=r["end"])
        assert astrorun.run_already_completed(path) is False
        return path

    def test_PLANTED_a_partial_ledger_makes_main_refuse(self, tmp_path, monkeypatch, capsys):
        plan = [row(i) for i in range(20)]
        p = tmp_path / "plan.json"
        p.write_text(json.dumps({"requests": plan}))
        led = self.partial_ledger(tmp_path / "l.jsonl", plan, 7)
        act = tmp_path / "activity.json"
        act.write_text(json.dumps({"rule": {"producer": dict(astrorun.AMENDMENT_2_RULE_FINGERPRINT)}}))
        reg = tmp_path / "reg.json"
        reg.write_text("{}")
        monkeypatch.setattr(sender, "PLAN", p)
        monkeypatch.setattr(sender, "REVIEWED_DIGEST", astrorun.plan_digest(plan))
        monkeypatch.setattr(sender, "LEDGER", led)
        monkeypatch.setattr(sender, "ACTIVITY", act)
        monkeypatch.setattr(sender, "REGISTRATION", reg)
        monkeypatch.setattr(astrorun, "ASTROREG2_AUTHORISATION", astrorun.ASTROREG2_AUTHORISATION_AS_RELAYED)
        monkeypatch.setattr(astrorun, "SUPERVISOR_SIGNOFF", 'wrote "dry run reviewed" at 13:00')
        monkeypatch.setattr(astrorun, "ASTROREG2_CAP", 20)
        # the activity result is a temp file, so the real git predicate cannot judge it; stubbing it is
        # what lets this test reach the PARTIAL-RUN clause rather than stopping at the activity clause
        monkeypatch.setattr(sender, "committed_in_git", lambda p: True)
        monkeypatch.setattr(sys, "argv", ["astroreg2_send.py", "--send"])
        assert sender.main() == 2
        said = capsys.readouterr().out
        assert "a partial run exists (7 charged)" in said
        assert "a resume needs Albert's new word" in said

    def test_the_partial_refusal_reads_as_a_new_authorisation_not_a_continuation(self):
        said = astrorun.A_CRASH_STOP_IS_NOT_A_CONTINUATION
        assert "NEW authorisation in" in said
        assert "not a continuation" in said
        assert "must not quietly become the exception" in said

    def test_PLANTED_a_resume_sends_exactly_the_uncharged_remainder(self, tmp_path):
        """k charged, len(plan) - k sent, and no element charged twice."""
        plan = [row(i) for i in range(20)]
        led = self.partial_ledger(tmp_path / "l.jsonl", plan, 7)
        client = RecordingPathOnly()
        budget = astrorun.RequestBudget(led, cap=astrorun.ASTROREG2_CAP)
        out = sender.send(
            plan,
            budget,
            lambda **kw: client.through_recording_path(kw["element_id"]),
            pilot_requests=99,
        )
        assert out["requests_charged_before_this_pass"] == 7
        assert out["requests_this_pass"] == 13 == len(plan) - 7
        assert client.calls_through_recording_path == 13, "the first seven are not bought again"

    def test_PLANTED_the_ledger_then_holds_each_element_exactly_once(self, tmp_path):
        plan = [row(i) for i in range(20)]
        led = self.partial_ledger(tmp_path / "l.jsonl", plan, 7)
        client = RecordingPathOnly()
        sender.send(
            plan,
            astrorun.RequestBudget(led, cap=astrorun.ASTROREG2_CAP),
            lambda **kw: client.through_recording_path(kw["element_id"]),
            pilot_requests=99,
        )
        charged = [
            json.loads(x)["element"]
            for x in led.read_text().splitlines()
            if json.loads(x).get("event") == "request"
        ]
        assert len(charged) == 20
        assert len(set(charged)) == 20, "no element appears twice in the ledger"
        assert sorted(charged) == sorted(r["element"] for r in plan)

    def test_a_restart_that_iterated_the_plan_would_have_double_paid(self, tmp_path):
        """The defect, demonstrated on the OLD shape, so the fix is measured and not asserted."""
        plan = [row(i) for i in range(20)]
        led = self.partial_ledger(tmp_path / "l.jsonl", plan, 7)
        budget = astrorun.RequestBudget(led, cap=astrorun.ASTROREG2_CAP)
        for r in plan:  # the old behaviour: iterate the PLAN
            budget.take(chrom=r["chrom"], element=r["element"])
        charged = [
            json.loads(x)["element"]
            for x in led.read_text().splitlines()
            if json.loads(x).get("event") == "request"
        ]
        assert len(charged) == 27, "7 re-paid on top of 20"
        assert len(charged) - len(set(charged)) == 7, "exactly the first seven paid for twice"

    def test_a_fresh_ledger_sends_the_whole_list(self, tmp_path):
        """The positive control: the remainder equals the plan when nothing is charged."""
        plan = [row(i) for i in range(6)]
        client = RecordingPathOnly()
        out = sender.send(
            plan,
            astrorun.RequestBudget(tmp_path / "l.jsonl", cap=50),
            lambda **kw: client.through_recording_path(kw["element_id"]),
            pilot_requests=99,
        )
        assert out["requests_charged_before_this_pass"] == 0
        assert out["requests_this_pass"] == 6


class TestTheErrorNoteAndTheLedgerPilot:
    def test_an_exception_is_noted_in_the_ledger_before_it_propagates(self, tmp_path):
        led = tmp_path / "l.jsonl"

        def boom(**kw):
            raise RuntimeError("quota exhausted")

        with pytest.raises(RuntimeError, match="quota exhausted"):
            sender.send([row(0)], astrorun.RequestBudget(led, cap=5), boom, pilot_requests=99)
        rows = [json.loads(x) for x in led.read_text().splitlines()]
        events = [r["event"] for r in rows]
        assert events == ["request", "error"], "the charge, then the error, both on the record"
        assert "RuntimeError: quota exhausted" in rows[1]["error"]
        assert rows[1]["element"] == "EH38E0000000"

    def test_a_charged_but_unanswered_element_is_visible_rather_than_implied(self, tmp_path):
        led = tmp_path / "l.jsonl"

        def boom(**kw):
            raise OSError("connection reset")

        with pytest.raises(OSError):
            sender.send([row(0)], astrorun.RequestBudget(led, cap=5), boom, pilot_requests=99)
        assert astrorun.ledger_charges(led)["charges"] == 1
        assert len(astrorun.answer_notes(led)) == 0, "no answer, and the error says why"

    def test_the_pilot_counts_answers_from_the_ledger_so_it_spans_a_resume(self, tmp_path):
        """Six answers in a first pass, four in a second: the pilot fires at the tenth overall."""
        plan = [row(i) for i in range(30)]
        led = tmp_path / "l.jsonl"
        client = RecordingPathOnly()
        first = sender.send(
            plan[:6],
            astrorun.RequestBudget(led, cap=astrorun.ASTROREG2_CAP),
            lambda **kw: client.through_recording_path(kw["element_id"]),
        )
        assert first["pilot"]["passed"] is True
        assert len(astrorun.answer_notes(led)) == 6
        out = sender.send(
            plan,
            astrorun.RequestBudget(led, cap=astrorun.ASTROREG2_CAP),
            lambda **kw: client.through_recording_path(kw["element_id"]),
        )
        assert out["answers_checked"] >= 10
        assert out["requests_charged_before_this_pass"] == 6

    def test_the_pilot_does_not_re_fire_on_a_resume_past_ten(self, tmp_path):
        plan = [row(i) for i in range(30)]
        led = tmp_path / "l.jsonl"
        client = RecordingPathOnly()
        sender.send(
            plan[:12],
            astrorun.RequestBudget(led, cap=astrorun.ASTROREG2_CAP),
            lambda **kw: client.through_recording_path(kw["element_id"]),
        )
        assert len(astrorun.answer_notes(led)) == 12
        out = sender.send(
            plan,
            astrorun.RequestBudget(led, cap=astrorun.ASTROREG2_CAP),
            lambda **kw: client.through_recording_path(kw["element_id"]),
        )
        assert out["pilot"]["passed"] is True, "past ten, the checkpoint is behind the run"
        assert out["requests_this_pass"] == 18


class TestEveryLiveEntryPointResolves:
    """The test that was missing, and the only kind that could have caught it.

    `enhancer_target.live_scorer_and_fetch()` was INVENTED: the name never existed. The sender passed 122
    of its own tests, a code review and a dry run, and still failed on that line, because every test
    SUBSTITUTED the entry point and the scorer is built only after the --send check. So the one line that
    only a real send could exercise was the one line that was wrong.

    These tests cost nothing, need no network and no key, and assert that every attribute the sender
    reaches for on a module it does not define actually resolves. If one name was invented, another may
    be, so this is derived from the syntax tree rather than from a list someone maintains.
    """

    SRC = ROOT / "scripts/astroreg2_send.py"

    def imported_modules(self):
        """The module objects the sender imports, by the name it binds them to."""
        import ast
        import contextlib
        import importlib

        tree = ast.parse(self.SRC.read_text())
        bound: dict[str, object] = {}
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
                for alias in node.names:
                    full = f"{node.module}.{alias.name}"
                    # a name imported FROM a module is not a submodule; the parent covers it
                    with contextlib.suppress(ModuleNotFoundError):
                        bound[alias.asname or alias.name] = importlib.import_module(full)
            elif isinstance(node, ast.Import):
                for alias in node.names:
                    with contextlib.suppress(ModuleNotFoundError):
                        bound[alias.asname or alias.name.split(".")[0]] = importlib.import_module(alias.name)
        return bound

    def attribute_uses(self):
        """Every `module.attr` the sender reaches for, from its syntax tree."""
        import ast

        tree = ast.parse(self.SRC.read_text())
        uses: set[tuple[str, str]] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name):
                uses.add((node.value.id, node.attr))
        return uses

    def test_every_attribute_the_sender_reaches_for_resolves(self):
        mods = self.imported_modules()
        missing = [
            f"{name}.{attr}"
            for name, attr in sorted(self.attribute_uses())
            if name in mods and not hasattr(mods[name], attr)
        ]
        assert missing == [], (
            f"the sender names {missing} on an imported module and they do not exist. This is the defect "
            "that stopped the first launch: an invented entry point, invisible because every test "
            "substituted it"
        )

    def test_the_check_would_have_caught_the_invented_name(self):
        """By planting the original defect: the check must fail on it, or it proves nothing."""
        from genomeos.predict import enhancer_target

        assert not hasattr(enhancer_target, "live_scorer_and_fetch"), (
            "the invented name is still absent from enhancer_target, which is why the sender defines its own"
        )
        mods = {"enhancer_target": enhancer_target}
        planted = {("enhancer_target", "live_scorer_and_fetch")}
        missing = [f"{n}.{a}" for n, a in planted if n in mods and not hasattr(mods[n], a)]
        assert missing == ["enhancer_target.live_scorer_and_fetch"]

    def test_the_senders_own_live_entry_point_exists_and_is_callable(self):
        assert callable(sender.live_scorer_and_fetch)

    def test_score_element_is_the_real_recording_path_and_takes_what_the_sender_passes(self):
        """Signature, not just existence: the sender's keywords must be ones it accepts."""
        import inspect

        from genomeos.predict import enhancer_target

        params = inspect.signature(enhancer_target.score_element).parameters
        for kw in ("chrom", "element_id", "start", "end"):
            assert kw in params, f"score_element does not take {kw}"
        positional = [n for n, p in params.items() if p.kind is p.POSITIONAL_OR_KEYWORD][:2]
        assert positional == ["scorer", "fetch"], (
            "the sender passes the scorer and fetcher positionally; if that order changed the send would "
            "break on the one line only a real send exercises"
        )

    def test_the_adapter_names_the_sender_uses_exist(self):
        """The scorer construction reaches into the adapter; those names must resolve too."""
        from genomeos.predict import AlphaGenomeAdapter

        a = AlphaGenomeAdapter.__new__(AlphaGenomeAdapter)
        assert hasattr(type(a), "available")
        assert callable(getattr(type(a), "_live_scorer", None)), (
            "the adapter exposes no public scorer accessor, so the sender uses the private one and says "
            "so; if a public one appears, the sender should move to it"
        )

    def test_the_genome_names_the_fetcher_uses_exist(self):
        from genomeos.genome import IndexedGenome, reference_fasta

        assert callable(reference_fasta)
        assert callable(getattr(IndexedGenome, "fetch", None))


class TestTheLedgerIsNotAResult:
    def test_the_ledger_is_outside_data_results(self):
        """A ledger is an append-only record of charges, not a computed value with a manifest."""
        assert "data/results" not in str(sender.LEDGER)
        assert str(sender.LEDGER) == "data/ledger/astroreg2_requests.jsonl"

    def test_the_reason_it_is_not_a_result_is_recorded(self):
        src = (ROOT / "scripts/astroreg2_send.py").read_text()
        assert "is not" in src and "RESULT" in src
        assert "append-only record of charges" in src
        assert "cannot be reconstructed from outcomes" in src
