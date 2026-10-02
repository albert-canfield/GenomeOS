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
        assert not any("alphagenome_adapter" in r for r in imported), (
            "the sender must reach AlphaGenome only through enhancer_target.score_element"
        )
        assert not any("alphagenome" in r.lower() for r in imported), (
            "no AlphaGenome client import of any kind"
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

    def test_main_refuses_today_and_returns_2_because_nothing_is_authorised(
        self, tmp_path, monkeypatch, capsys
    ):
        """End to end: the entry point refuses, prints why, and sends nothing."""
        assert astrorun.ASTROREG2_AUTHORISATION is None
        assert astrorun.SUPERVISOR_SIGNOFF is None
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
