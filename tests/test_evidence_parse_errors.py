"""A program that fails to parse is reported beside the rows, and the result says it is incomplete.

ROADMAP item 9 recorded that `genomeos evidence` once lost 21 facts (26,845 to 26,824) to one demo
program with a syntax error and said nothing. The report is an addition only: every row, summary and
per-program entry of the programs that do parse must come back byte for byte as it did without the
broken program, and the only new things are the failures (path and message), their count and
`complete: false`. A clean read adds nothing at all, so its payload is what it was before the report
existed; its empty `failed` list already says that nothing failed.
"""

import json
from pathlib import Path

import pytest

from genomeos import evidence
from genomeos.web.server import Api

BASE = """module test.base

param shared_rate = 0.5 h {
  evidence: curated "BioModels BIOMD0000000012"
  confidence: 0.8
}
"""

TOP = """module test.top

import base.bio

gene SOX2 {
  symbol: SOX2
  evidence: experimental "GENCODE v50"
  confidence: 0.95
}

param guessed = 0.35 {
  evidence: inferred "order of magnitude"
  confidence: 0.3
}

rule SOX2 activates SOX2 {
  strength: 1.0
  threshold: 1.0
  evidence: predicted "AlphaGenome"
}
"""

BROKEN = "module test.broken\n\ngene {\n"

# the keys collect returned before the report existed; the report may add keys, never change these
KEYS_BEFORE = {"rows", "files", "failed", "summary", "whole", "weak_line", "evidence_kinds"}
REPORT_KEYS = {"parse_errors", "parse_error_count", "complete"}


@pytest.fixture
def programs(tmp_path, monkeypatch):
    """Two programs that parse, one importing the other, in place of the project's own."""
    d = tmp_path / "data" / "demo"
    d.mkdir(parents=True)
    (d / "base.bio").write_text(BASE)
    (d / "top.bio").write_text(TOP)
    monkeypatch.setattr(evidence, "PROGRAM_DIRS", (Path("data/demo"),))
    evidence._cached.cache_clear()
    yield tmp_path
    evidence._cached.cache_clear()


def _bytes(obj) -> str:
    return json.dumps(obj, sort_keys=True, default=str)


def _collect(root: Path, **kw):
    evidence._cached.cache_clear()
    return evidence.collect(root, **kw)


def test_a_clean_run_adds_nothing_and_names_no_failure(programs):
    out = _collect(programs)

    assert set(out) == KEYS_BEFORE
    assert out["failed"] == []


def test_a_failed_read_adds_exactly_the_report(programs):
    (programs / "data" / "demo" / "broken.bio").write_text(BROKEN)

    out = _collect(programs)

    assert set(out) == KEYS_BEFORE | REPORT_KEYS


def test_a_syntax_error_is_named_with_path_and_message_and_the_result_is_incomplete(programs):
    clean = _collect(programs)
    (programs / "data" / "demo" / "broken.bio").write_text(BROKEN)

    out = _collect(programs)

    assert [e["path"] for e in out["parse_errors"]] == ["data/demo/broken.bio"]
    assert set(out["parse_errors"][0]) == {"path", "message"}
    assert out["parse_errors"][0]["message"].strip()
    assert out["parse_error_count"] == 1
    assert out["complete"] is False
    # the rows, both summaries and every parsed program's entry are what they were without it
    for key in ("rows", "summary", "whole", "weak_line", "evidence_kinds"):
        assert _bytes(out[key]) == _bytes(clean[key]), key
    parsed = [f for f in out["files"] if f["path"] != "data/demo/broken.bio"]
    assert _bytes(parsed) == _bytes(clean["files"])


def test_a_program_importing_a_broken_one_is_reported_too_and_the_rest_are_untouched(programs):
    demo = programs / "data" / "demo"
    (demo / "leaf.bio").write_text(BROKEN)
    (demo / "user.bio").write_text("module test.user\n\nimport leaf.bio\n")
    base_alone = _collect(programs, module="data/demo/base.bio")

    out = _collect(programs)

    assert sorted(e["path"] for e in out["parse_errors"]) == ["data/demo/leaf.bio", "data/demo/user.bio"]
    assert out["parse_error_count"] == 2 and out["complete"] is False
    assert _bytes(_collect(programs, module="data/demo/base.bio")["rows"]) == _bytes(base_alone["rows"])


def test_a_filter_does_not_hide_the_failures(programs):
    """The report is about the whole read, so a query that selects nothing still says it is incomplete."""
    (programs / "data" / "demo" / "broken.bio").write_text(BROKEN)

    out = _collect(programs, query="no fact says this")

    assert out["rows"] == [] and out["parse_error_count"] == 1 and out["complete"] is False


def test_the_api_passes_the_report_through_also_for_the_csv(programs):
    (programs / "data" / "demo" / "broken.bio").write_text(BROKEN)
    evidence._cached.cache_clear()
    api = Api(programs)

    page = api.evidence()
    sheet = api.evidence(csv=True)

    for out in (page, sheet):
        assert [e["path"] for e in out["parse_errors"]] == ["data/demo/broken.bio"]
        assert out["parse_error_count"] == 1 and out["complete"] is False
    assert set(sheet) == {"csv", "matched"} | REPORT_KEYS


def test_the_api_of_a_clean_run_adds_nothing(programs):
    api = Api(programs)

    assert set(api.evidence(csv=True)) == {"csv", "matched"}
    assert not REPORT_KEYS & set(api.evidence())


def test_the_cli_says_incomplete_and_keeps_the_rows_it_prints(programs, monkeypatch, capsys):
    from genomeos.cli import main

    monkeypatch.chdir(programs)
    evidence._cached.cache_clear()
    assert main(["evidence"]) == 0
    clean = capsys.readouterr()
    assert "incomplete" not in clean.out and "incomplete" not in clean.err
    (programs / "data" / "demo" / "broken.bio").write_text(BROKEN)

    evidence._cached.cache_clear()
    assert main(["evidence"]) == 1
    broken = capsys.readouterr()

    err = broken.err.splitlines()
    assert any(
        line.startswith("genomeos evidence: ERROR data/demo/broken.bio did not parse: ") for line in err
    )
    assert "genomeos evidence: incomplete: 1 program(s) failed to parse" in broken.err
    out = broken.out.splitlines()
    assert "incomplete: 1 program(s) failed to parse" in out[1]
    # everything else on stdout is what the clean run printed, the headline less its failure mark
    rest = [out[0].replace(", 1 FAILED TO PARSE (facts missing)", ""), *out[2:]]
    assert rest == clean.out.splitlines()


def test_the_cli_csv_says_incomplete_and_writes_the_same_file(programs, monkeypatch, capsys):
    from genomeos.cli import main

    monkeypatch.chdir(programs)
    evidence._cached.cache_clear()
    assert main(["evidence", "--csv", "clean.csv"]) == 0
    clean = capsys.readouterr()
    (programs / "data" / "demo" / "broken.bio").write_text(BROKEN)

    evidence._cached.cache_clear()
    assert main(["evidence", "--csv", "broken.csv"]) == 1
    broken = capsys.readouterr()

    assert (programs / "broken.csv").read_bytes() == (programs / "clean.csv").read_bytes()
    assert "incomplete: 1 program(s) failed to parse" in broken.out
    assert broken.out.splitlines()[0] == clean.out.splitlines()[0].replace("clean.csv", "broken.csv")
