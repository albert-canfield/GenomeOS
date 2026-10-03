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
import tracked_paths as tp

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


#: A small track axis for the stubs: two brain tracks and two cell lines, which is the shape of the
#: problem item (g) is about -- a writer that kept only the cell lines would look fine on a stub.
STUB_TRACKS = ("astrocyte", "brain cortex", "K562", "HepG2")


def full_axis(rows: int = 1, tracks=STUB_TRACKS, values=None) -> dict:
    """A gene axis that KEEPS the whole track vector, which is what a usable paid answer looks like."""
    return {
        "output": 0,
        "rows": [
            {
                "row": gi,
                "gene_name": "AAA",
                "tracks_total": len(tracks),
                astrorun.FULL_VECTOR_KEY: astrorun.full_vector_entry(
                    values if values is not None else [0.1 * (gi + 1) * (ti + 1) for ti in range(len(tracks))]
                ),
            }
            for gi in range(rows)
        ],
        "tracks_total": len(tracks),
        "track_names": list(tracks),
        "tracks_sha256": "0" * 64,
    }


def good_answer(element_id: str) -> dict:
    """A usable answer: named effects AND every gene row's whole track vector.

    The vectors were added when item (g) landed. Before that, every stub here returned what the code
    asked for -- four cell lines -- which is precisely why no stub could show the defect.
    """
    return {
        "id": element_id,
        "genes": [{"gene": "AAA", "mean_log2fc": -0.4, "n_tracks": 3}],
        "model": {"gene_axis_outputs": [full_axis()]},
    }


def answer_without_full_vectors(element_id: str) -> dict:
    """What the writer produced at 8c6622c: effects, and the brain tracks thrown away."""
    axis = full_axis()
    axis["rows"] = [{k: v for k, v in r.items() if k != astrorun.FULL_VECTOR_KEY} for r in axis["rows"]]
    return {
        "id": element_id,
        "genes": [{"gene": "AAA", "mean_log2fc": -0.4, "n_tracks": 3}],
        "model": {"gene_axis_outputs": [axis]},
    }


def nameless_answer(element_id: str) -> dict:
    """The latent fault: a gene proto without `name` yields ZERO effects, silently."""
    return {"id": element_id, "genes": []}


def writing_score(cache_root: Path, answer=None):
    """A score double that WRITES the cache entry, because that is what score_element does.

    The earlier doubles returned an answer and wrote nothing, so the pilot -- which now reads the written
    entry -- would find nothing. Returning without writing is the substitution that hid the cache-hit
    defect, so the doubles write.
    """

    def score(**kw):
        eid, chrom = kw["element_id"], kw["chrom"]
        hit = (answer or good_answer)(eid)
        d = Path(cache_root) / chrom
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{eid}.json").write_text(json.dumps(hit))
        return {k: v for k, v in hit.items() if k != "genes"}  # score_element strips genes

    return score


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
        root = tmp_path / "c"

        def score(**kw):
            client.through_recording_path(kw["element_id"])
            return writing_score(root)(**kw)

        out = sender.send([row(i) for i in range(12)], budget, score, cache_root=root)
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
        # create_client IS the sanctioned construction: the sweep that produced the frozen K562 values
        # builds its client with it, and the sender must build the SAME pinned client or it is buying a
        # different feature. What must not happen is a request bypassing score_element, which the
        # planted direct-call test above is the guarantee of.
        assert "genomeos.predict.alphagenome_adapter.create_client" in imported, (
            "the sender must build the pinned client the same way the frozen sweep did"
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
    def answers(self, bad_indices, cache_root):
        """A writing double: the pilot reads the WRITTEN entry, so a double must write one."""

        def pick(eid):
            i = int(eid.replace("EH38E", ""))
            return nameless_answer(eid) if i in bad_indices else good_answer(eid)

        return writing_score(cache_root, pick)

    def test_PLANTED_a_nameless_response_in_the_first_ten_stops_the_run(self, tmp_path):
        """At most ten requests spent, the failures recorded, and the run terminal."""
        led = tmp_path / "l.jsonl"
        budget = astrorun.RequestBudget(led, cap=astrorun.ASTROREG2_CAP)
        root = tmp_path / "c"
        out = sender.send([row(i) for i in range(50)], budget, self.answers({4}, root), cache_root=root)
        assert out["requests_sent"] == 10, "it stops AT the checkpoint, not after the whole list"
        assert out["pilot"]["passed"] is False
        assert out["pilot"]["stopped"]["at_request"] == 10
        assert [f["id"] for f in out["pilot"]["stopped"]["failing"]] == ["EH38E0000004"]
        assert astrorun.run_already_completed(led) is True, "a pilot stop is terminal"
        # The stop must be THIS guard's stop and not some other end to the loop: the record carries
        # the pilot's OWN rule and its OWN consequence, word for word, and the ledger holds the
        # `pilot_stop` note. Without these the assertions above are satisfied by any early exit.
        assert out["pilot"]["stopped"]["rule"] == sender.PILOT_RULE
        assert (
            "EVERY one of them must carry a gene name and a non-empty effects set"
            in (out["pilot"]["stopped"]["rule"])
        )
        assert out["pilot"]["stopped"]["consequence"] == sender.PILOT_STOP_CONSUMES_THE_ONE_RUN
        notes = [json.loads(line) for line in led.read_text().splitlines()]
        stops = [n for n in notes if n.get("event") == "pilot_stop"]
        assert len(stops) == 1, "the stop is on the ledger, which is the durable record of it"
        assert stops[0]["at_request"] == 10
        assert stops[0]["failing"] == 1

    def test_PLANTED_a_VECTOR_LESS_response_in_the_first_ten_also_stops_the_run(self, tmp_path):
        """The other half of the rule, which had no send-level plant: effects kept, vectors thrown away.

        `answer_is_usable` was tested on a vector-less answer in isolation, but nothing showed that
        `send` STOPS on one. Item (g)'s whole point is that a paid answer missing its track vector
        cannot be widened without paying again, so the pilot is the last cheap moment to catch it --
        and an answer that is NAMED, with non-empty effects, passes the half of the rule the other
        plant exercises. So this is a separate fact and gets its own plant.
        """
        led = tmp_path / "l.jsonl"
        budget = astrorun.RequestBudget(led, cap=astrorun.ASTROREG2_CAP)
        root = tmp_path / "c"
        out = sender.send(
            [row(i) for i in range(50)],
            budget,
            writing_score(
                root,
                lambda eid: (
                    answer_without_full_vectors(eid)
                    if int(eid.replace("EH38E", "")) == 7
                    else good_answer(eid)
                ),
            ),
            cache_root=root,
        )
        assert out["requests_sent"] == 10, "it stops AT the checkpoint, not after the whole list"
        assert out["pilot"]["passed"] is False
        assert out["pilot"]["stopped"]["at_request"] == 10
        failing = out["pilot"]["stopped"]["failing"]
        assert [f["id"] for f in failing] == ["EH38E0000007"]
        # the cause named: this one is NAMED with effects and fails only on the discarded vector
        assert failing[0]["gene_name_present"] is True
        assert failing[0]["effects_non_empty"] is True
        assert failing[0]["full_vectors_kept"] is False
        assert out["pilot"]["stopped"]["rule"] == sender.PILOT_RULE
        assert out["pilot"]["stopped"]["consequence"] == sender.PILOT_STOP_CONSUMES_THE_ONE_RUN
        assert astrorun.run_already_completed(led) is True, "a pilot stop is terminal"

    def test_the_stop_says_it_consumes_albert_s_one_run_rather_than_exhausting_a_budget(self, tmp_path):
        out = sender.send(
            [row(i) for i in range(50)],
            astrorun.RequestBudget(tmp_path / "l.jsonl", cap=astrorun.ASTROREG2_CAP),
            self.answers({0}, tmp_path / "c"),
            cache_root=tmp_path / "c",
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
            self.answers(set(), tmp_path / "c"),
            cache_root=tmp_path / "c",
        )
        assert out["requests_sent"] == 14
        assert out["pilot"]["passed"] is True

    def test_a_failure_AFTER_the_checkpoint_does_not_stop_the_run_but_is_recorded(self, tmp_path):
        """The checkpoint is a gate at ten, not a filter over the whole list."""
        out = sender.send(
            [row(i) for i in range(14)],
            astrorun.RequestBudget(tmp_path / "l.jsonl", cap=50),
            self.answers({12}, tmp_path / "c"),
            cache_root=tmp_path / "c",
        )
        assert out["requests_sent"] == 14
        assert out["pilot"]["passed"] is True
        assert [c["id"] for c in out["unusable_answers"]] == ["EH38E0000012"]

    def test_an_answer_with_a_blank_gene_name_is_unusable(self):
        blank = {**good_answer("x"), "genes": [{"gene": "  "}]}
        empty = {**good_answer("x"), "genes": []}
        assert sender.answer_is_usable(blank)["usable"] is False
        assert sender.answer_is_usable(empty)["usable"] is False
        assert sender.answer_is_usable(good_answer("x"))["usable"] is True
        # and item (g): named effects are not enough if the paid vector was discarded
        thin = sender.answer_is_usable(answer_without_full_vectors("x"))
        assert thin["usable"] is False
        assert thin["full_vectors_kept"] is False
        assert "carries no track_vector" in thin["full_vectors_refusal"]
        assert "BOUGHT AND THROWN AWAY" in astrorun.WHY_THE_FULL_TRACK_VECTOR_IS_KEPT

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
        monkeypatch.setattr(astrorun, "recorded_signoff_words", lambda *a, **k: None)
        plan = [row(i) for i in range(3)]
        p = tmp_path / "plan.json"
        p.write_text(json.dumps({"requests": plan}))
        monkeypatch.setattr(sender, "PLAN", p)
        monkeypatch.setattr(sender, "REVIEWED_DIGEST", astrorun.plan_digest(plan))
        # main() reassigns the module's LEDGER from --run, and a global set inside a call outlives the
        # call; recording it here is what makes teardown put it back instead of leaking into later tests
        monkeypatch.setattr(sender, "LEDGER", sender.LEDGER)
        # --run is required with --send, and run 1's ledger is complete; naming a run is what lets this
        # test reach the authorisation clause rather than stopping at the missing-run clause
        monkeypatch.setattr(sys, "argv", ["astroreg2_send.py", "--send", "--run", "2"])
        assert sender.main() == 2
        said = capsys.readouterr().out
        assert "REFUSED" in said
        assert "no AstroREG-2 approval is recorded" in said

    def test_main_refuses_when_send_names_no_run(self, monkeypatch, capsys):
        """A run that inherited run 1's complete ledger would refuse for the wrong reason, or resume it."""
        monkeypatch.setattr(sys, "argv", ["astroreg2_send.py", "--send"])
        assert sender.main() == 2
        said = capsys.readouterr().out
        assert "REFUSED: --send requires --run N" in said
        assert "a new run gets its own" in said

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
        # --run resolves the ledger now, so the partial ledger has to arrive through that resolution
        monkeypatch.setattr(astrorun, "ledger_for_run", lambda n: led)
        monkeypatch.setattr(sender, "ACTIVITY", act)
        monkeypatch.setattr(sender, "REGISTRATION", reg)
        monkeypatch.setattr(astrorun, "ASTROREG2_AUTHORISATION", astrorun.ASTROREG2_AUTHORISATION_AS_RELAYED)
        monkeypatch.setattr(
            astrorun, "recorded_signoff_words", lambda *a, **k: 'wrote "dry run reviewed" at 13:00'
        )
        monkeypatch.setattr(astrorun, "ASTROREG2_CAP", 20)
        # Run 1's real approval is CONSUMED, so this test grants the run it uses one; its subject is the
        # partial-run clause, not the approval clause, which has its own planted tests below.
        monkeypatch.setattr(
            astrorun,
            "ASTROREG2_AUTHORISATIONS",
            {1: {"words": "<run 1's words>", "requests": 20, "consumed": False}},
        )
        # the activity result is a temp file, so the real git predicate cannot judge it; stubbing it is
        # what lets this test reach the PARTIAL-RUN clause rather than stopping at the activity clause
        monkeypatch.setattr(sender, "committed_in_git", lambda p: True)
        # the blob check is about a DIFFERENT clause and now refuses first, since this lane changed both
        # signed files; stubbing it is what lets this test reach the partial-run clause it is about
        monkeypatch.setattr(astrorun, "check_signoff_closure", lambda *a, **k: {})
        monkeypatch.setattr(sys, "argv", ["astroreg2_send.py", "--send", "--run", "1"])
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
        root = tmp_path / "c"
        first = sender.send(
            plan[:6],
            astrorun.RequestBudget(led, cap=astrorun.ASTROREG2_CAP),
            writing_score(root),
            cache_root=root,
        )
        assert first["pilot"]["passed"] is True
        assert len(astrorun.answer_notes(led)) == 6
        out = sender.send(
            plan,
            astrorun.RequestBudget(led, cap=astrorun.ASTROREG2_CAP),
            writing_score(root),
            cache_root=root,
        )
        assert out["answers_checked"] >= 10
        assert out["requests_charged_before_this_pass"] == 6

    def test_the_pilot_does_not_re_fire_on_a_resume_past_ten(self, tmp_path):
        plan = [row(i) for i in range(30)]
        led = tmp_path / "l.jsonl"
        root = tmp_path / "c"
        sender.send(
            plan[:12],
            astrorun.RequestBudget(led, cap=astrorun.ASTROREG2_CAP),
            writing_score(root),
            cache_root=root,
        )
        assert len(astrorun.answer_notes(led)) == 12
        out = sender.send(
            plan,
            astrorun.RequestBudget(led, cap=astrorun.ASTROREG2_CAP),
            writing_score(root),
            cache_root=root,
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

    def test_the_sender_defines_no_convenience_helper_of_its_own(self):
        """Ruling: build from what exists. Inventing the name would make it real and keep the habit."""
        assert not hasattr(sender, "live_scorer_and_fetch")
        src = (ROOT / "scripts/astroreg2_send.py").read_text()
        assert "create_client(" in src and "_live_scorer(" in src

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
        # asserted on the DECLARATION, not on sender.LEDGER: main() reassigns that global from --run,
        # so a test reading it would be judging whichever run a previous test last named
        assert str(astrorun.LEDGER_RUN1) == "data/ledgers/astroreg2.jsonl"
        assert str(sender.LEDGER) == "data/ledgers/astroreg2.jsonl"
        for run in (1, 2, 7):
            assert "data/results" not in str(astrorun.ledger_for_run(run))
            assert str(astrorun.ledger_for_run(run)).startswith("data/ledgers/")

    def test_the_reason_it_is_not_a_result_is_recorded(self):
        src = (ROOT / "scripts/astroreg2_send.py").read_text()
        assert "is not" in src and "RESULT" in src
        assert "append-only record of charges" in src
        assert "cannot be reconstructed from outcomes" in src


class TestTheFrozenFeatureParametersCannotDiverge:
    """The third shape of today's lesson: the stub substituted behaviour, and the divergence was a
    PARAMETER.

    The frozen K562 values came from a sweep passing threshold=0.0 explicitly. The adapter DEFAULTS to
    0.05. A sender taking that default would censor every |effect| < 0.05 -- values the frozen cache
    carries -- making the astrocyte deletion term a different feature from the one the weights were frozen
    on. No stub could catch it, so these tests read both call sites instead.
    """

    SWEEP = ROOT / "scripts/enhancer_targets_all.py"
    SENDER = ROOT / "scripts/astroreg2_send.py"

    def threshold_at(self, path, func_names):
        """Every threshold= argument passed to _live_scorer in a file, by literal or by constant name."""
        import ast

        tree = ast.parse(path.read_text())
        out = []
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr in func_names
            ):
                for kw in node.keywords:
                    if kw.arg == "threshold":
                        if isinstance(kw.value, ast.Constant):
                            out.append(kw.value.value)
                        elif isinstance(kw.value, ast.Attribute):
                            out.append(f"{ast.unparse(kw.value)}")
                        else:
                            out.append(ast.unparse(kw.value))
                else:
                    if not any(kw.arg == "threshold" for kw in node.keywords):
                        out.append("<DEFAULT>")
        return out

    def test_the_sweep_passes_the_frozen_threshold_explicitly(self):
        got = self.threshold_at(self.SWEEP, {"_live_scorer"})
        assert got, "the sweep must call _live_scorer somewhere"
        assert all(v == astrorun.FROZEN_SCORER_THRESHOLD for v in got), (
            f"the sweep that produced the FROZEN values passes {got}; the constant is "
            f"{astrorun.FROZEN_SCORER_THRESHOLD}"
        )

    def test_the_sender_passes_the_same_constant_and_never_the_default(self):
        got = self.threshold_at(self.SENDER, {"_live_scorer"})
        assert got, "the sender must call _live_scorer"
        assert "<DEFAULT>" not in got, (
            "the sender takes the adapter's DEFAULT threshold of 0.05, which censors every |effect| "
            "below it and changes the feature definition the weights were frozen on"
        )
        assert got == ["astrorun.FROZEN_SCORER_THRESHOLD"], got

    def test_PLANTED_a_default_call_fails_this_check(self, tmp_path):
        """A constant only one side uses is not a fix, so the check is shown failing on the defect."""
        bad = tmp_path / "defaulted.py"
        bad.write_text("adapter._live_scorer()\n")
        assert self.threshold_at(bad, {"_live_scorer"}) == ["<DEFAULT>"]
        good = tmp_path / "explicit.py"
        good.write_text("adapter._live_scorer(threshold=0.0)\n")
        assert self.threshold_at(good, {"_live_scorer"}) == [0.0]

    def test_PLANTED_a_different_threshold_fails_this_check(self, tmp_path):
        other = tmp_path / "other.py"
        other.write_text("adapter._live_scorer(threshold=0.05)\n")
        got = self.threshold_at(other, {"_live_scorer"})
        assert got == [0.05]
        assert not all(v == astrorun.FROZEN_SCORER_THRESHOLD for v in got)

    def test_the_client_timeout_matches_the_sweeps(self):
        """Ruling 2: prove the construction is the sweep's rather than asserting it."""
        import re

        sweep = self.SWEEP.read_text()
        m = re.search(r"CALL_TIMEOUT\s*=\s*(\d+)", sweep)
        assert m, "the sweep must define CALL_TIMEOUT"
        assert int(m.group(1)) == astrorun.FROZEN_CLIENT_TIMEOUT
        assert "timeout=CALL_TIMEOUT" in sweep
        assert "timeout=astrorun.FROZEN_CLIENT_TIMEOUT" in self.SENDER.read_text()

    def test_both_sides_build_the_client_through_create_client(self):
        for path in (self.SWEEP, self.SENDER):
            assert "create_client(" in path.read_text(), f"{path.name} must use the pinned factory"

    def test_the_threshold_effect_is_demonstrated_not_assumed(self):
        """The counterfactual the cache cannot supply: what the two thresholds DO to a small effect.

        The cached answers store no raw response -- eight gene-row keys and no ninth -- so the stored
        by_cell values cannot be re-derived through any scorer. What CAN be shown at zero cost is the
        mechanism: a |effect| of 0.02 survives a 0.0 threshold and is dropped by a 0.05 one. That is the
        whole of the defect, and it is why the sender must pass 0.0.
        """
        small = 0.02
        assert abs(small) >= astrorun.FROZEN_SCORER_THRESHOLD, "kept by the frozen threshold"
        assert abs(small) < 0.05, "dropped by the adapter's default, which is the defect"

    def test_the_run_records_which_feature_definition_produced_the_answers(self):
        src = self.SENDER.read_text()
        for field in ("scorer_threshold", "client_timeout", "model_version"):
            assert f'"{field}"' in src, f"the result must record {field}"


class TestTheSignoffIsBoundToTheWholeClosure:
    """Rebuilt after a FIXED-POINT defect in the first version, and widened from two files to the closure.

    The first version stored the signed blob of genomeos/attribution/astrorun.py INSIDE that same file.
    Writing a blob value there changes the file, which changes its blob, so the recorded value can never
    equal the computed one: no re-sign by anyone could ever have passed. The record now lives OUTSIDE the
    code, in data/ledgers/astroreg2_signoff.json, read at runtime and never imported, which dissolves the
    fixed point rather than working around it.

    It also signs the sender's WHOLE IMPORT CLOSURE. Two files were never the right boundary: an edit to
    enhancer_target.score_element or to alphagenome_adapter would have passed unnoticed, and those two
    decide what is BOUGHT and what is RECORDED.
    """

    def test_the_fixed_point_is_gone_because_the_record_is_not_in_the_code(self):
        import inspect

        src = inspect.getsource(astrorun)
        assert "SUPERVISOR_SIGNOFF_BLOBS" not in src, "the blob pair must not live in the code it hashes"
        assert str(astrorun.SIGNOFF_RECORD) == "data/ledgers/astroreg2_signoff.jsonl"
        assert "fixed point with no solution" in src, (
            "the reason the record left the code must stay written down beside it"
        )

    def test_a_file_cannot_contain_its_own_hash_which_is_why_this_moved(self, tmp_path):
        """The defect, demonstrated rather than reasoned about."""
        f = tmp_path / "selfref.py"
        f.write_text('BLOB = "x" * 40\n')
        first = astrorun.git_blob(f, root=astrorun.ROOT_FOR_BLOBS)
        f.write_text(f'BLOB = "{first}"\n')
        second = astrorun.git_blob(f, root=astrorun.ROOT_FOR_BLOBS)
        assert second != first, "writing a blob into a file changes that file's blob"

    def test_the_closure_covers_what_decides_what_is_bought_and_recorded(self):
        closure = astrorun.sender_closure()
        assert astrorun.SENDER_ENTRY in closure
        for must in (
            "genomeos/predict/enhancer_target.py",
            "genomeos/predict/alphagenome_adapter.py",
            "genomeos/attribution/astrorun.py",
        ):
            assert must in closure, f"{must} must be signed: it decides what is bought or recorded"
        assert len(closure) > 10, "a closure of two files is the boundary this replaced"

    def test_it_refuses_while_no_record_exists(self):
        """A record outside the repository cannot be vouched for by a commit, so it refuses first."""
        with pytest.raises(astrorun.SendRefusedError, match="outside the repository"):
            astrorun.check_signoff_closure(record=None, path=Path("/nonexistent/signoff.jsonl"))

    def test_an_EMPTY_committed_record_signs_nothing(self):
        """The 'no record' branch, reached with the record supplied rather than read from a file."""
        with pytest.raises(astrorun.SendRefusedError, match="signs nothing"):
            astrorun.check_signoff_closure(record={"words": "dry run reviewed"})

    def test_a_record_with_no_signed_closure_signs_nothing(self):
        with pytest.raises(astrorun.SendRefusedError, match="signs nothing"):
            astrorun.check_signoff_closure(record={"words": "dry run reviewed"})

    def test_a_MATCHING_closure_passes(self):
        """The positive control: an always-refusing check would look safe and be useless."""
        rec = {
            "signed_at": "2026-10-02T17:23:05+01:00",
            "words": "dry run reviewed",
            "signed_closure": astrorun.sender_closure(),
        }
        out = astrorun.check_signoff_closure(record=rec)
        assert out["files_signed"] == len(rec["signed_closure"])
        assert out["words"] == "dry run reviewed"

    def test_PLANTED_an_edit_DEEP_in_the_closure_refuses(self):
        """alphagenome_adapter is not the sender, and changing it must still stop the send."""
        closure = astrorun.sender_closure()
        deep = "genomeos/predict/alphagenome_adapter.py"
        assert deep in closure
        tampered = dict(closure)
        tampered[deep] = "0" * 40
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.check_signoff_closure(record={"signed_closure": tampered})
        said = str(exc.value)
        assert "CHANGED" in said and deep in said
        assert "RE-SIGN is required" in said
        assert "not a missing authorisation" in said

    def test_PLANTED_an_ADDED_module_in_the_closure_refuses(self):
        """The case a hand-maintained path list would miss entirely."""
        closure = astrorun.sender_closure()
        short = {k: v for k, v in closure.items() if k != "genomeos/predict/enhancer_target.py"}
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.check_signoff_closure(record={"signed_closure": short})
        said = str(exc.value)
        assert "ADDED to the closure" in said
        assert "enhancer_target.py" in said

    def test_PLANTED_a_REMOVED_module_refuses_and_is_named(self):
        closure = dict(astrorun.sender_closure())
        closure["genomeos/predict/gone_module.py"] = "0" * 40
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.check_signoff_closure(record={"signed_closure": closure})
        assert "REMOVED from the closure" in str(exc.value)
        assert "gone_module.py" in str(exc.value)

    def test_the_closure_is_computed_by_the_shared_function_not_a_local_copy(self):
        import inspect

        src = inspect.getsource(astrorun.sender_closure)
        assert "mf.counting_path" in src, "eight copies of this closure had already drifted once"

    def test_the_record_is_read_at_runtime_and_never_imported(self):
        import inspect

        src = inspect.getsource(astrorun.signoff_lines)
        assert "json.loads" in src and "read_text" in src
        assert "import_module" not in src
        assert "import_module" not in inspect.getsource(astrorun.recorded_signoff)

    def test_a_malformed_record_is_treated_as_absent_rather_than_trusted(self, tmp_path):
        bad = tmp_path / "signoff.jsonl"
        bad.write_text("{not json")
        assert astrorun.recorded_signoff(bad) is None, "an unreadable governing line authorises nothing"
        assert len(astrorun.signoff_lines(bad)) == 1, "and the bad line is kept as a marker"


class TestTheSweepAndTheSenderReachTheSamePath:
    """Code-path equivalence for the frozen feature, PINNED rather than remembered.

    This replaces a check that could not be run: re-deriving cached by_cell values needs a stored raw
    response and the cache holds none. What can be pinned statically is that the sweep which produced the
    frozen K562 values and the sender which will buy the astrocyte ones both reach the SAME
    `score_element`, with the same cache and the same MIN_EFFECT -- and that the extra arguments the sweep
    passes cannot reach the cached gene rows that `crispri` reads.

    The equivalence is currently a fact about two call sites that a future edit to either could break
    silently, and the whole frozen-feature claim rests on it.
    """

    SWEEP = ROOT / "scripts/enhancer_targets_all.py"
    SENDER = ROOT / "scripts/astroreg2_send.py"
    TARGET = ROOT / "genomeos/predict/enhancer_target.py"

    def test_the_sweep_reaches_score_element_through_context_score(self):
        import ast
        import inspect

        from genomeos.predict.enhancer_target import Context

        sweep = self.SWEEP.read_text()
        assert ".score(" in sweep, "the sweep must call Context.score"
        import textwrap

        src = textwrap.dedent(inspect.getsource(Context.score))
        calls = {
            node.func.id
            for node in ast.walk(ast.parse(src))
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        }
        assert "score_element" in calls, "Context.score must reach score_element"

    def test_the_sender_reaches_score_element_directly(self):
        import ast

        tree = ast.parse(self.SENDER.read_text())
        reached = {
            node.func.attr
            for node in ast.walk(tree)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
        }
        assert "score_element" in reached, "the sender must reach score_element"

    def test_both_use_score_elements_own_cache_and_min_effect(self):
        """Neither side overrides them, so the cached rows are written identically."""
        import inspect

        from genomeos.predict import enhancer_target

        params = inspect.signature(enhancer_target.score_element).parameters
        assert params["cache"].default is enhancer_target.CACHE
        assert params["min_effect"].default == enhancer_target.MIN_EFFECT
        sender = self.SENDER.read_text()
        assert "min_effect=" not in sender, (
            "the sender must not override the effect floor the frozen values used"
        )
        assert "cache=astrorun.ASTROREG2_CACHE_ROOT" in sender, (
            "the sender MUST override the cache root: against the sweep's root every planned element is "
            "a cache hit, so the run would charge and buy nothing, and a file written there would shadow "
            "the frozen archive"
        )
        src = inspect.getsource(enhancer_target.Context.score)
        assert "cache=cache" in src and "min_effect=min_effect" in src
        sig = inspect.signature(enhancer_target.Context.score).parameters
        assert sig["cache"].default is enhancer_target.CACHE
        assert sig["min_effect"].default == enhancer_target.MIN_EFFECT

    def test_the_extra_arguments_the_sweep_passes_cannot_reach_the_cached_rows(self):
        """The substance of the equivalence: they feed predict_target, AFTER the cache is written.

        `crispri` reads the cached gene rows. If inferred_targets, domain_genes or coding could influence
        what is cached, the sweep and the sender would write different rows and the frozen feature would
        not be the same feature. They appear only after the cache write, so they cannot.
        """
        import inspect

        from genomeos.predict import enhancer_target

        src = inspect.getsource(enhancer_target.score_element)
        write = src.index("p.write_text(json.dumps(hit))")
        before = src[:write]
        for name in ("inferred_targets", "domain_genes", "coding"):
            assert name not in before.split("def score_element")[-1].split('"""')[-1], (
                f"{name} appears before the cache is written, so it could change the cached rows that "
                "crispri reads"
            )
        after = src[write:]
        assert "predict_target(" in after, "those arguments feed predict_target, after the cache write"

    def test_the_cached_rows_hold_no_raw_response_so_a_replay_check_is_impossible(self):
        """Why the output-level comparability check was not run, recorded as a fact not an excuse."""
        import gzip
        import json

        arc_path = Path("data/knowledge/alphagenome/elements/chr21.json.gz")
        if not arc_path.exists():
            pytest.skip("the chr21 element archive is not on this machine")
        with gzip.open(arc_path, "rt") as fh:
            arc = json.load(fh)
        keys, gene_keys = set(), set()
        for hit in list(arc.values())[:50]:
            keys |= set(hit)
            for g in (hit.get("genes") or [])[:3]:
                gene_keys |= set(g)
        assert not [
            k for k in keys | gene_keys if any(w in k.lower() for w in ("raw", "response", "adata", "proto"))
        ]
        assert gene_keys <= {
            "by_cell",
            "gene",
            "max_drop_log2fc",
            "max_drop_tissue",
            "max_rise_log2fc",
            "max_rise_tissue",
            "mean_log2fc",
            "n_tracks",
        }


class TestACacheHitIsRefusedNotCharged:
    """The test whose absence let the whole thing through.

    Every earlier test stubbed the scorer, and `score_element`'s cache lookup happens BEFORE the scorer
    is ever called -- so no stub could reach it. Against the sweep's cache root all 1,232 planned
    elements are cache hits, measured: the run would have charged 1,232 to the ledger, sent nothing to
    the service, and produced a result made entirely of K562 answers labelled as astrocyte work. It
    would have reported success.

    These tests use the REAL archive layout, because a stub standing in for the archive is precisely the
    substitution that hid the defect.
    """

    def archive(self, root: Path, chrom: str, entries: dict) -> Path:
        """A chromosome archive in the layout load_cached actually reads."""
        import gzip

        root.mkdir(parents=True, exist_ok=True)
        p = root / f"{chrom}.json.gz"
        with gzip.open(p, "wt") as fh:
            json.dump(entries, fh)
        return p

    def answer(self, element_id: str, genes=1) -> dict:
        return {
            "id": element_id,
            "chrom": "chr1",
            "genes": [
                {"gene": f"G{i}", "mean_log2fc": -0.3, "n_tracks": 4, "by_cell": {"K562": -0.3}}
                for i in range(genes)
            ],
        }

    def test_PLANTED_an_element_already_in_the_archive_is_REFUSED(self, tmp_path):
        """(b) and (d): the real layout, and the refusal names the count."""
        plan = [row(i) for i in range(4)]
        root = tmp_path / "elements_astroreg2"
        self.archive(root, "chr1", {plan[2]["element"]: self.answer(plan[2]["element"])})
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.check_nothing_in_the_plan_is_already_cached(plan, root)
        said = str(exc.value)
        assert "1 of the 4 planned elements already have a cached answer" in said
        assert "charged and NOT bought" in said
        assert plan[2]["element"] in said

    def test_a_loose_per_element_file_is_found_too_because_it_SHADOWS_the_archive(self, tmp_path):
        """The shadowing hazard, in the direction that matters for the pre-flight."""
        plan = [row(0)]
        root = tmp_path / "elements_astroreg2"
        d = root / "chr1"
        d.mkdir(parents=True)
        (d / f"{plan[0]['element']}.json").write_text(json.dumps(self.answer(plan[0]["element"])))
        with pytest.raises(astrorun.SendRefusedError):
            astrorun.check_nothing_in_the_plan_is_already_cached(plan, root)

    def test_an_EMPTY_request_root_passes(self, tmp_path):
        """The positive control: the only state in which a request is really made."""
        plan = [row(i) for i in range(4)]
        out = astrorun.check_nothing_in_the_plan_is_already_cached(plan, tmp_path / "fresh")
        assert out["already_cached"] == 0
        assert out["plan"] == 4

    def test_the_defect_is_measured_on_the_real_sweep_root(self):
        """What would have happened: every planned element already answered, by the K562 sweep."""
        from genomeos.predict.enhancer_target import load_cached

        sweep = Path("data/knowledge/alphagenome/elements")
        if not (sweep / "chr1.json.gz").exists():
            pytest.skip("the sweep archive is not on this machine")
        plan = json.loads(Path("data/results/astroreg2_request_plan.json").read_text())["requests"][:25]
        hits = [r for r in plan if load_cached(r["chrom"], r["element"], sweep) is not None]
        assert len(hits) == len(plan), (
            "every sampled element is already answered in the sweep root, which is why the run would "
            "have bought nothing"
        )
        one = load_cached(plan[0]["chrom"], plan[0]["element"], sweep)
        by_cell = (one["genes"][0].get("by_cell") or {}) if one.get("genes") else {}
        assert "astrocyte" not in by_cell, (
            "and the cached answers carry no astrocyte value, so they could not have served the test"
        )

    def test_the_request_root_is_not_the_sweeps(self):
        assert str(astrorun.ASTROREG2_CACHE_ROOT) == "data/knowledge/alphagenome/elements_astroreg2"
        assert "SHADOW THE" in Path("genomeos/attribution/astrorun.py").read_text()

    def test_the_sender_passes_the_separate_root_to_score_element(self):
        src = (ROOT / "scripts/astroreg2_send.py").read_text()
        assert "cache=astrorun.ASTROREG2_CACHE_ROOT" in src


class TestThePilotReadsTheWrittenEntry:
    """(c) The pilot was reading score_element's return, which has `genes` STRIPPED.

    So it was always an empty list and always failed, whatever the service returned. The failure
    accidentally stopped a run that would have bought nothing -- but it would equally have stopped a run
    that was working.
    """

    def test_score_element_really_does_strip_the_genes_key(self):
        """The premise, checked against the function rather than remembered."""
        import inspect

        from genomeos.predict import enhancer_target

        src = inspect.getsource(enhancer_target.score_element)
        assert 'k != "genes"' in src, "if this stops being true, the pilot's source should be revisited"

    def test_PLANTED_a_written_entry_WITH_genes_passes(self, tmp_path):
        import gzip

        root = tmp_path / "c"
        root.mkdir()
        entry = good_answer("E1")
        with gzip.open(root / "chr1.json.gz", "wt") as fh:
            json.dump({"E1": entry}, fh)
        got = sender.written_entry("chr1", "E1", root)
        assert got is not None
        assert sender.answer_is_usable(got)["usable"] is True

    def test_PLANTED_a_written_entry_WITHOUT_genes_fails(self, tmp_path):
        import gzip

        root = tmp_path / "c"
        root.mkdir()
        with gzip.open(root / "chr1.json.gz", "wt") as fh:
            json.dump({"E1": {"id": "E1", "genes": []}}, fh)
        got = sender.written_entry("chr1", "E1", root)
        assert sender.answer_is_usable(got)["usable"] is False

    def test_a_missing_written_entry_is_unusable_rather_than_an_error(self, tmp_path):
        assert sender.written_entry("chr1", "nope", tmp_path) is None
        assert sender.answer_is_usable({})["usable"] is False

    def test_the_pilot_records_whether_it_read_a_written_entry(self, tmp_path):
        """So a run cannot pass the pilot on answers it never found."""
        src = (ROOT / "scripts/astroreg2_send.py").read_text()
        assert "from_written_entry" in src


class TestCrispriReadsTheSameSeparateRoot:
    """(e) The astrocyte annotation must read the root the answers were written to."""

    def test_an_element_cache_can_be_pointed_at_the_astroreg2_root(self):
        from genomeos.attribution import crispri

        cache = crispri.ElementCache(astrorun.ASTROREG2_CACHE_ROOT)
        assert cache.root == astrorun.ASTROREG2_CACHE_ROOT
        assert cache.root != crispri.ELEMENT_CACHE, "it must not be the sweep's root"

    def test_it_reads_a_value_from_that_root_by_the_frozen_collapse_rule(self, tmp_path):
        import gzip

        from genomeos.attribution import crispri

        root = tmp_path / "elements_astroreg2"
        root.mkdir()
        entry = {
            "id": "E1",
            "genes": [{"gene": "AAA", "n_tracks": 2, "by_cell": {"astrocyte": -0.7, "K562": -0.1}}],
        }
        with gzip.open(root / "chr1.json.gz", "wt") as fh:
            json.dump({"E1": entry}, fh)
        cache = crispri.ElementCache(root)
        assert cache.value("chr1", "E1", "AAA", "astrocyte") == -0.7
        assert cache.value("chr1", "E1", "AAA", "K562") == -0.1
        assert cache.value("chr1", "E1", "AAA", "HepG2") is None

    def test_the_sweeps_root_holds_no_astrocyte_value_so_it_cannot_serve_the_test(self):
        from genomeos.attribution import crispri

        if not (crispri.ELEMENT_CACHE / "chr1.json.gz").exists():
            pytest.skip("the sweep archive is not on this machine")
        plan = json.loads(Path("data/results/astroreg2_request_plan.json").read_text())["requests"][0]
        cache = crispri.ElementCache()
        genes = plan["serves_genes"][:3]
        assert all(cache.value(plan["chrom"], plan["element"], g, "astrocyte") is None for g in genes), (
            "the sweep's answers carry no astrocyte value, which is the whole reason a run is needed"
        )


class TestTheSignoffRecordIsAppendOnly:
    """One record per line, and ONLY THE LAST LINE GOVERNS.

    The single-object format let a new sign-off REPLACE the old one, which destroys the history of what
    was signed and makes a withdrawal indistinguishable from never having happened. The concrete danger
    is not hypothetical: the 19:14 record signed the closure of a sender that would have charged 1,232
    requests and bought nothing, so a format in which that line could speak again on a code revert would
    reauthorise the defect.
    """

    def line(self, closure, why="", supersedes=None, words="dry run reviewed"):
        return {
            "signed_at": "2026-10-02T20:05:00+01:00",
            "words": words,
            "signed_closure": closure,
            "signed_closure_sha256": astrorun.closure_sha256(closure),
            "supersedes": supersedes,
            "why": why,
        }

    def write(self, path, records):
        path.write_text("".join(json.dumps(r) + "\n" for r in records))
        return path

    def test_the_record_is_a_jsonl_and_the_old_object_is_marked_superseded(self):
        assert str(astrorun.SIGNOFF_RECORD) == "data/ledgers/astroreg2_signoff.jsonl"
        assert str(astrorun.SIGNOFF_RECORD_SUPERSEDED) == "data/ledgers/astroreg2_signoff.json"
        assert "NEVER WRITTEN AGAIN" in Path("genomeos/attribution/astrorun.py").read_text()

    def test_only_the_last_line_is_read(self, tmp_path):
        f = self.write(
            tmp_path / "s.jsonl",
            [self.line({"a": "1"}, words="first"), self.line({"b": "2"}, words="second")],
        )
        assert len(astrorun.signoff_lines(f)) == 2
        assert astrorun.recorded_signoff(f)["words"] == "second"
        assert astrorun.recorded_signoff_words(f) == "second"

    def test_PLANTED_an_EARLIER_line_matching_the_current_code_does_NOT_authorise(self, tmp_path):
        """The one that matters: a withdrawn sign-off must not revive when code reverts."""
        current = astrorun.sender_closure()
        f = self.write(
            tmp_path / "s.jsonl",
            [
                self.line(current, why="the cache-defect closure, withdrawn"),
                self.line({"genomeos/x.py": "0" * 40}, why="supersedes the withdrawn one"),
            ],
        )
        governing = astrorun.recorded_signoff(f)
        assert governing["why"] == "supersedes the withdrawn one", "the LAST line governs"
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.check_signoff_closure(record=governing)
        said = str(exc.value)
        assert "RE-SIGN is required" in said
        assert "ONLY THE LAST LINE" in astrorun.ONLY_THE_LAST_LINE_GOVERNS

    def test_the_LAST_line_matching_the_current_code_DOES_authorise(self, tmp_path):
        """Near-miss positive control: the same two lines in the other order."""
        current = astrorun.sender_closure()
        f = self.write(
            tmp_path / "s.jsonl",
            [self.line({"genomeos/x.py": "0" * 40}, why="old"), self.line(current, why="re-signed")],
        )
        governing = astrorun.recorded_signoff(f)
        assert governing["why"] == "re-signed"
        out = astrorun.check_signoff_closure(record=governing)
        assert out["files_signed"] == len(current)

    def commit_repo(self, tmp_path, text):
        import subprocess

        base = tmp_path / "repo"
        base.mkdir()
        subprocess.run(["git", "init", "-q"], cwd=base, check=True)
        rel = Path("rec.jsonl")
        (base / rel).write_text(text)
        subprocess.run(["git", "add", "-A"], cwd=base, check=True)
        subprocess.run(
            ["git", "-c", "user.email=t@t", "-c", "user.name=t", "commit", "-qm", "x"],
            cwd=base,
            check=True,
        )
        return base, rel

    def test_a_COMMITTED_appended_line_passes(self, tmp_path):
        """Appending is the point -- but the appended line must be committed to govern.

        This replaces a test that asserted an UNCOMMITTED appended line passes. That was the defect:
        HEAD being a prefix of the working file is exactly what an append produces, so a prefix rule
        authorised anyone who could write the file.
        """
        import subprocess

        base, rel = self.commit_repo(tmp_path, '{"a": 1}\n')
        with (base / rel).open("a") as fh:
            fh.write('{"a": 2}\n')
        with pytest.raises(astrorun.SendRefusedError, match="not byte-identical to HEAD"):
            astrorun.check_signoff_is_committed(rel, root=base)
        subprocess.run(["git", "add", "-A"], cwd=base, check=True)
        subprocess.run(
            ["git", "-c", "user.email=t@t", "-c", "user.name=t", "commit", "-qm", "append"],
            cwd=base,
            check=True,
        )
        out = astrorun.check_signoff_is_committed(rel, root=base)
        assert out["identical_to_head"] is True, "once committed, the appended line governs"

    def test_PLANTED_an_EDITED_line_refuses(self, tmp_path):
        base, rel = self.commit_repo(tmp_path, '{"a": 1}\n{"a": 2}\n')
        (base / rel).write_text('{"a": 1}\n{"a": 99}\n')
        with pytest.raises(astrorun.SendRefusedError, match="not byte-identical to HEAD"):
            astrorun.check_signoff_is_committed(rel, root=base)

    def test_PLANTED_a_REMOVED_line_refuses_too(self, tmp_path):
        base, rel = self.commit_repo(tmp_path, '{"a": 1}\n{"a": 2}\n')
        (base / rel).write_text('{"a": 1}\n')
        with pytest.raises(astrorun.SendRefusedError, match="not byte-identical to HEAD"):
            astrorun.check_signoff_is_committed(rel, root=base)

    def test_a_record_not_in_HEAD_REFUSES_rather_than_being_reported(self, tmp_path):
        """Replaces a test asserting the opposite. 'Not checked' read as a pass is the vacuous pass."""
        import subprocess

        base = tmp_path / "repo"
        base.mkdir()
        subprocess.run(["git", "init", "-q"], cwd=base, check=True)
        (base / "x.txt").write_text("x")
        subprocess.run(["git", "add", "-A"], cwd=base, check=True)
        subprocess.run(
            ["git", "-c", "user.email=t@t", "-c", "user.name=t", "commit", "-qm", "x"],
            cwd=base,
            check=True,
        )
        (base / "rec.jsonl").write_text('{"a": 1}\n')
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.check_signoff_is_committed(Path("rec.jsonl"), root=base)
        assert "An uncommitted sign-off authorises nothing" in str(exc.value)

    def test_an_unreadable_last_line_does_not_authorise(self, tmp_path):
        f = tmp_path / "s.jsonl"
        f.write_text(json.dumps(self.line({"a": "1"})) + "\n{not json\n")
        assert astrorun.recorded_signoff(f) is None
        assert len(astrorun.signoff_lines(f)) == 2, "the bad line is kept as a marker, not skipped"

    def test_each_line_carries_the_digest_of_what_it_signed(self):
        c = astrorun.sender_closure()
        assert len(astrorun.closure_sha256(c)) == 64
        assert astrorun.closure_sha256(c) == astrorun.closure_sha256(dict(reversed(list(c.items()))))


class TestEachRunHasItsOwnLedger:
    def test_run_one_keeps_its_filename_and_later_runs_do_not(self):
        assert astrorun.ledger_for_run(1) == Path("data/ledgers/astroreg2.jsonl")
        assert astrorun.ledger_for_run(2) == Path("data/ledgers/astroreg2_run2.jsonl")
        assert astrorun.ledger_for_run(3) == Path("data/ledgers/astroreg2_run3.jsonl")

    def test_a_run_id_below_one_is_refused(self):
        with pytest.raises(ValueError, match="1 or more"):
            astrorun.ledger_for_run(0)

    def test_PLANTED_a_second_run_on_the_FIRST_ledger_refuses(self):
        """Run 1 completed, so its ledger carries run_complete and must never be appended to again."""
        first = Path("data/ledgers/astroreg2.jsonl")
        # run 1's ledger is TRACKED: it is the committed record of what was spent, so it is on every
        # machine that has this commit and the skip could never fire. The plant that proves a second
        # run on a completed ledger is refused is now unconditional, which is what a plant must be.
        tp.must_be_present(first, was="run 1's ledger is not on this machine")
        assert astrorun.run_already_completed(first) is True
        assert astrorun.ledger_charges(first)["charges"] == 10

    def test_a_NEW_run_on_a_NEW_ledger_has_nothing_charged_and_is_not_complete(self, tmp_path):
        fresh = tmp_path / "astroreg2_run2.jsonl"
        assert astrorun.run_already_completed(fresh) is False
        assert astrorun.ledger_charges(fresh)["charges"] == 0

    def test_send_requires_run_and_has_no_default(self):
        src = (ROOT / "scripts/astroreg2_send.py").read_text()
        assert '"--run"' in src
        assert "--send requires --run N" in src
        assert "deliberately without a default" in src

    def test_the_reason_a_run_gets_its_own_ledger_is_recorded(self):
        assert "record of THAT run" in astrorun.A_NEW_RUN_GETS_ITS_OWN_LEDGER

    def test_the_result_states_mixed_model_version_provenance(self):
        src = (ROOT / "scripts/astroreg2_send.py").read_text()
        assert "MIXED MODEL-VERSION PROVENANCE" in src
        assert "unrequested" in src and "ALL_FOLDS" in src


class TestAnUncommittedSignoffAuthorisesNothing:
    """Three defects that each let an UNCOMMITTED sign-off authorise a send.

    (a) a missing HEAD copy returned {"checked": False} and the send went on -- a check that reports
        "not checked" and is read as a pass, which is the same vacuous-pass shape as a tracer recording
        zero reads and satisfying "every read is declared";
    (b) HEAD being a byte PREFIX of the working file is exactly what an APPEND produces, so anyone able
        to write the file could append a governing line and send;
    (c) an absolute path made `git show HEAD:<path>` fail, which landed in (a) -- so a tampered file
        passed.

    The send-time rule is therefore IDENTICAL-TO-HEAD, not prefix: the governing line must be a
    COMMITTED line, so a sign-off is an act in the repository's history rather than a file someone has
    edited.
    """

    def repo(self, tmp_path, committed: str | None, working: str):
        import subprocess

        base = tmp_path / "repo"
        base.mkdir()
        subprocess.run(["git", "init", "-q"], cwd=base, check=True)
        rel = Path("rec.jsonl")
        (base / "seed.txt").write_text("seed")
        if committed is not None:
            (base / rel).write_text(committed)
        subprocess.run(["git", "add", "-A"], cwd=base, check=True)
        subprocess.run(
            ["git", "-c", "user.email=t@t", "-c", "user.name=t", "commit", "-qm", "x"],
            cwd=base,
            check=True,
        )
        (base / rel).write_text(working)
        return base, rel

    def test_PLANTED_no_HEAD_copy_REFUSES(self, tmp_path):
        """(a) An uncommitted sign-off authorises nothing."""
        base, rel = self.repo(tmp_path, committed=None, working='{"a": 1}\n')
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.check_signoff_is_committed(rel, root=base)
        assert "An uncommitted sign-off authorises nothing" in str(exc.value)

    def test_PLANTED_an_uncommitted_APPENDED_line_REFUSES(self, tmp_path):
        """(b) Append-only without committed authorises anybody who can write the file."""
        base, rel = self.repo(tmp_path, committed='{"a": 1}\n', working='{"a": 1}\n{"a": 2}\n')
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.check_signoff_is_committed(rel, root=base)
        said = str(exc.value)
        assert "not byte-identical to HEAD" in said
        assert "uncommitted" in said

    def test_PLANTED_an_absolute_path_to_a_tampered_file_REFUSES(self, tmp_path):
        """(c) The path form must not be able to turn a tampered file into a pass."""
        base, rel = self.repo(tmp_path, committed='{"a": 1}\n', working='{"a": 99}\n')
        with pytest.raises(astrorun.SendRefusedError):
            astrorun.check_signoff_is_committed(base / rel, root=base)

    def test_a_path_OUTSIDE_the_repository_refuses(self, tmp_path):
        base, rel = self.repo(tmp_path, committed='{"a": 1}\n', working='{"a": 1}\n')
        outside = tmp_path / "elsewhere.jsonl"
        outside.write_text('{"a": 1}\n')
        with pytest.raises(astrorun.SendRefusedError, match="outside the repository"):
            astrorun.check_signoff_is_committed(outside, root=base)

    def test_NEAR_MISS_committed_and_identical_PASSES(self, tmp_path):
        """A check that can never pass is the defect the fixed-point case already taught us."""
        base, rel = self.repo(tmp_path, committed='{"a": 1}\n', working='{"a": 1}\n')
        out = astrorun.check_signoff_is_committed(rel, root=base)
        assert out["identical_to_head"] is True
        assert out["bytes"] == len('{"a": 1}\n')

    def test_PLANTED_may_send_is_WIRED_to_the_committed_check_and_refuses_today(self, tmp_path, monkeypatch):
        """A mechanism that exists but is not wired is indistinguishable from one that does not exist.

        So the subject here is the WIRING, not the check: may_send is driven to the sign-off clause with
        every earlier clause satisfied, and it must refuse because the sign-off record has no committed
        copy at HEAD. If the call were removed from may_send, this returns and sends.

        THE RECORD IS A PLANT AND NOT THE REAL ONE, and that is the whole correction here. This test
        used to drive the clause against the real `data/ledgers/astroreg2_signoff.jsonl` and relied on
        it being ABSENT FROM THE REPOSITORY -- so the day run 2's sign-off was committed (`f029d47`),
        the committed check started passing, `may_send` walked on to the next clause, and the test went
        red reporting "the governing sign-off line covers run 2, not run 1". It was asserting the
        wiring by way of a fact about the repository's contents, which the first real sign-off was
        always going to falsify. `SIGNOFF_RECORD` is read from the module at call time by both
        `check_signoff_is_committed` and `recorded_signoff`, so pointing it at a repo-relative path
        that no commit holds tests the wiring against a condition this test OWNS.
        """
        reg = tmp_path / "reg.json"
        reg.write_text("{}")
        act = tmp_path / "activity.json"
        act.write_text(json.dumps({"rule": {"producer": dict(astrorun.AMENDMENT_2_RULE_FINGERPRINT)}}))
        plan = [row(i) for i in range(astrorun.ASTROREG2_CAP)]
        monkeypatch.setattr(astrorun, "ASTROREG2_AUTHORISATION", astrorun.ASTROREG2_AUTHORISATION_AS_RELAYED)
        monkeypatch.setattr(
            astrorun,
            "ASTROREG2_AUTHORISATIONS",
            {1: {"words": "<run 1's words>", "requests": astrorun.ASTROREG2_CAP, "consumed": False}},
        )
        absent = Path("data/ledgers/astroreg2_signoff_PLANT_absent.jsonl")
        monkeypatch.setattr(astrorun, "SIGNOFF_RECORD", absent)
        with pytest.raises(astrorun.SendRefusedError) as exc:
            astrorun.may_send(
                plan=plan,
                registration=reg,
                activity_result=act,
                reviewed_digest=astrorun.plan_digest(plan),
                run_id=1,
                ledger_root=tmp_path,
                signoff='wrote "dry run reviewed" at 13:00',
                committed=lambda p: True,
            )
        said = str(exc.value)
        # The PLANTED name, so the refusal is demonstrably about the record this test controls and not
        # about whatever the repository happens to hold today.
        assert absent.as_posix() in said, said
        assert "astroreg2_signoff" in said
        assert "no committed copy at HEAD" in said, said
        assert "An uncommitted sign-off authorises nothing" in said

    def test_the_SENDER_hands_may_send_a_RUN_ID_and_no_ledger(self):
        """The ledger the clauses are checked against must not be a value this script chose.

        Planted against the call site as it stood at 20fe5ee, which passed `ledger=LEDGER` -- and
        LEDGER came from `--run`, so the run picked the file its one-run clause was judged against.
        """
        import ast

        tree = ast.parse((ROOT / "scripts/astroreg2_send.py").read_text())
        calls = [
            c
            for c in ast.walk(tree)
            if isinstance(c, ast.Call) and isinstance(c.func, ast.Attribute) and c.func.attr == "may_send"
        ]
        assert len(calls) == 1, "one gate, called once"
        kw = {k.arg for k in calls[0].keywords}
        assert "run_id" in kw, "the run id is the input"
        assert "ledger" not in kw, (
            "a ledger handed in is the defect: may_send must DERIVE it from the run id, or the caller "
            "chooses what it is judged against"
        )
        budgets = [
            c
            for c in ast.walk(tree)
            if isinstance(c, ast.Call)
            and isinstance(c.func, ast.Attribute)
            and c.func.attr == "astroreg2_budget"
        ]
        assert len(budgets) == 1
        assert {k.arg for k in budgets[0].keywords} == {"run_id"}, "and the budget likewise"
        assert not budgets[0].args, "no positional ledger either"

    def test_may_send_does_NOT_hand_the_check_a_record_which_would_skip_it(self):
        """check_signoff_closure(record=...) skips the committed check by design, for tests.

        That makes `record` a way past the rule, so the call site is pinned: may_send must pass no
        record. Planted against a call site that supplied one.
        """
        import ast

        tree = ast.parse((ROOT / "genomeos/attribution/astrorun.py").read_text())
        fn = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == "may_send")
        calls = [
            c
            for c in ast.walk(fn)
            if isinstance(c, ast.Call)
            and isinstance(c.func, ast.Name)
            and c.func.id == "check_signoff_closure"
        ]
        assert len(calls) == 1, "may_send must check the sign-off closure exactly once"
        kw = {k.arg for k in calls[0].keywords}
        assert "record" not in kw, (
            "a record handed in skips check_signoff_is_committed, so may_send must let the check read "
            "the file it is the rule about"
        )
        assert not calls[0].args, "no positional record either"

    def test_the_governing_line_must_be_a_committed_line(self):
        src = Path("genomeos/attribution/astrorun.py").read_text()
        assert "act in the repository's history" in src
        assert "An uncommitted sign-off authorises nothing" in src
