from pathlib import Path

import pytest

from genomeos import evidence
from genomeos.web.server import Api, ApiError

ROOT = Path(__file__).resolve().parent.parent

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

protein NANOG { evidence: experimental "GENCODE v50"; confidence: 0.9 }

param guessed = 0.35 {
  evidence: inferred "order of magnitude"
  confidence: 0.3
}

rule SOX2 activates NANOG {
  strength: 1.0
  threshold: 1.0
  evidence: predicted "AlphaGenome"
  confidence: 0.4
}
"""


@pytest.fixture
def programs(tmp_path, monkeypatch):
    """Two programs, one importing the other, in place of the project's own."""
    d = tmp_path / "data" / "demo"
    d.mkdir(parents=True)
    (d / "base.bio").write_text(BASE)
    (d / "top.bio").write_text(TOP)
    monkeypatch.setattr(evidence, "PROGRAM_DIRS", (Path("data/demo"),))
    evidence._cached.cache_clear()
    return tmp_path


def test_a_fact_is_credited_to_the_program_that_declares_it(programs):
    out = evidence.collect(programs)
    by_path: dict[str, set[str]] = {}
    for r in out["rows"]:
        by_path.setdefault(r["path"], set()).add(r["label"])
    assert "shared_rate = 0.5 h" in by_path["data/demo/base.bio"]
    # top.bio imports base.bio, so the imported parameter is not counted again under top
    assert not any("shared_rate" in label for label in by_path["data/demo/top.bio"])
    assert {"SOX2", "NANOG", "guessed = 0.35", "SOX2 activates NANOG"} <= by_path["data/demo/top.bio"]
    assert out["whole"]["facts"] == 5


def test_filters_and_summary(programs):
    out = evidence.collect(programs)
    assert out["whole"]["by_evidence"] == {"experimental": 2, "curated": 1, "predicted": 1, "inferred": 1}
    assert out["whole"]["weak"] == 2  # the inferred parameter at 0.3 and the predicted rule at 0.4
    weak = evidence.collect(programs, max_confidence=0.5)["rows"]
    assert [r["confidence"] for r in weak] == [0.3, 0.4]  # weakest first
    kinds = evidence.collect(programs, kinds={"experimental", "curated"})["rows"]
    assert {r["evidence"] for r in kinds} == {"experimental", "curated"}
    assert evidence.collect(programs, query="alphagenome")["rows"][0]["label"] == "SOX2 activates NANOG"
    one = evidence.collect(programs, module="data/demo/base.bio")
    assert one["summary"]["facts"] == 1 and one["whole"]["facts"] == 5
    # the filtered summary describes the selection, the whole one the project
    assert one["summary"]["mean_confidence"] == 0.8


def test_csv_is_the_review_list(programs):
    rows = evidence.collect(programs, max_confidence=0.5)["rows"]
    lines = evidence.to_csv(rows).strip().splitlines()
    assert lines[0].startswith("confidence,evidence,block,label")
    assert lines[1].startswith("0.3,inferred,parameter,guessed = 0.35")
    assert len(lines) == 3


def test_an_unparsable_program_is_reported_not_raised(programs):
    (programs / "data" / "demo" / "broken.bio").write_text("module broken\n\ngene {\n")
    evidence._cached.cache_clear()
    out = evidence.collect(programs)
    broken = next(f for f in out["files"] if f["path"].endswith("broken.bio"))
    assert broken["error"] and broken["facts"] == 0
    assert out["whole"]["facts"] == 5  # the other programs still count


def test_the_real_programs_are_read():
    out = evidence.collect(ROOT)
    assert out["whole"]["facts"] > 5000 and len(out["files"]) >= 20
    assert not [f for f in out["files"] if f.get("error")]
    assert set(out["whole"]["by_evidence"]) <= set(evidence.KINDS)
    order = lambda r: (not r["stated"], r["confidence"], r["path"], r["block"], r["label"])  # noqa: E731
    assert out["rows"] == sorted(out["rows"], key=order)


def test_api_evidence():
    api = Api(ROOT)
    out = api.evidence(kinds="inferred", max_confidence="0.5", limit=5)
    assert out["matched"] > 0 and len(out["rows"]) <= 5
    assert all(r["evidence"] == "inferred" and r["confidence"] <= 0.5 for r in out["rows"])
    assert "csv" not in out and out["whole"]["facts"] > out["matched"]
    assert api.evidence(kinds="inferred", limit=1, csv=True)["csv"].startswith("confidence,evidence")
    with pytest.raises(ApiError):
        api.evidence(kinds="bogus")
    with pytest.raises(ApiError):
        api.evidence(max_confidence="soon")


def test_the_compiled_chromosomes_are_off_by_default_and_on_by_request(tmp_path, monkeypatch) -> None:
    """They hold 940,803 of the project's 967,426 facts and cost twenty seconds to parse.

    A page that loaded them by default would be a page nobody opens, so `programs` takes the flag and
    the parse itself is skipped - filtering them out after parsing would have cost the same time.
    """
    from genomeos import evidence as ev

    (tmp_path / "data" / "demo").mkdir(parents=True)
    (tmp_path / "data" / "demo" / "hand.bio").write_text("module demo.hand\n")
    (tmp_path / ev.COMPILED_DIR).mkdir(parents=True)
    (tmp_path / ev.COMPILED_DIR / "noncoding_chr21.bio").write_text("module human.noncoding.chr21\n")

    assert ev.programs(tmp_path) == ["data/demo/hand.bio"]
    assert str(ev.COMPILED_DIR) in " ".join(ev.programs(tmp_path, compiled=True))


def test_a_missing_compiled_directory_is_not_an_error(tmp_path) -> None:
    """The programs are generated and git-ignored, so a fresh checkout has none."""
    from genomeos import evidence as ev

    (tmp_path / "data" / "demo").mkdir(parents=True)

    assert ev.programs(tmp_path, compiled=True) == []


def test_the_pooled_mean_confidence_is_reported_as_the_mixture_it_is() -> None:
    """A level is a property of a population, so one mean over four kinds is on no scale at all.

    2026-09-17's calibration lane tested the predicted band against four measured populations and
    found the stated confidence outside its own 95% interval in 10 of 12 bands, with the four offsets
    in OPPOSITE directions: +0.472, -0.291, -0.137, +0.518. No constant fixes four disagreeing signs.

    The same applies one level up to `mean_confidence`, which averages curated, inferred, predicted
    and experimental facts. On this repository the pooled figure sits between two populations more
    than half a point apart and describes neither, so the per-kind means are printed beside it.
    """
    rows = [
        {"evidence": "predicted", "confidence": 0.2, "block": "b", "label": "x", "source": ""},
        {"evidence": "predicted", "confidence": 0.3, "block": "b", "label": "y", "source": ""},
        {"evidence": "curated", "confidence": 0.9, "block": "b", "label": "z", "source": "s"},
    ]

    out = evidence.summarise(rows)

    assert out["mean_confidence"] == round((0.2 + 0.3 + 0.9) / 3, 3)
    by = out["mean_confidence_by_evidence"]
    assert by["predicted"] == {"facts": 2, "stated": 2, "unstated": 0, "mean": 0.25}
    assert by["curated"] == {"facts": 1, "stated": 1, "unstated": 0, "mean": 0.9}
    # the pooled figure lies between the two and equals neither: that is the whole point
    assert by["predicted"]["mean"] < out["mean_confidence"] < by["curated"]["mean"]
    assert "not on one scale" in out["mean_confidence_is_a_mixture"]


def test_a_single_kind_still_reports_its_own_mean() -> None:
    """The mixture warning must not imply a mixture where there is none."""
    rows = [{"evidence": "curated", "confidence": 0.9, "block": "b", "label": "z", "source": "s"}]

    out = evidence.summarise(rows)

    assert out["mean_confidence_by_evidence"] == {
        "curated": {"facts": 1, "stated": 1, "unstated": 0, "mean": 0.9}
    }
    assert out["mean_confidence"] == 0.9


def test_a_program_that_fails_to_parse_is_counted_and_the_command_fails(programs, monkeypatch, capsys):
    """2026-09-27: a demo program with a syntax error took 21 facts out of `genomeos evidence`
    (26,845 to 26,824) and the command said nothing. A failure is now named, counted and an error."""
    from genomeos.cli import main

    (programs / "data" / "demo" / "broken.bio").write_text("module broken\n\ngene {\n")
    evidence._cached.cache_clear()
    out = evidence.collect(programs)
    assert [f["path"] for f in out["failed"]] == ["data/demo/broken.bio"]
    assert out["failed"][0]["error"]
    monkeypatch.chdir(programs)
    assert main(["evidence"]) == 1
    cap = capsys.readouterr()
    assert "broken.bio did not parse" in cap.err and "1 program(s) failed to parse" in cap.err
    assert "5 facts in 2 programs, 1 FAILED TO PARSE" in cap.out


def test_a_clean_run_reports_no_failures_and_succeeds(programs, monkeypatch, capsys):
    from genomeos.cli import main

    assert evidence.collect(programs)["failed"] == []
    monkeypatch.chdir(programs)
    assert main(["evidence"]) == 0
    cap = capsys.readouterr()
    assert "FAILED" not in cap.out and "ERROR" not in cap.err


UNSTATED_PROGRAM = """module test.unstated

param judged_zero = 0.1 {
  evidence: inferred "a guess somebody judged worthless"
  confidence: 0.0
}

param judged_weak = 0.2 { evidence: inferred "order of magnitude"; confidence: 0.4 }

param judged_strong = 0.3 { evidence: curated "BioModels"; confidence: 0.9 }

gene NANOG { symbol: NANOG; evidence: curated "GENCODE v50" }

rule NANOG activates NANOG {
  strength: 0.2
  evidence: predicted "AlphaGenome deletion" effect -0.2 log2 fold change, probability unavailable
}
"""


@pytest.fixture
def unstated_program(tmp_path, monkeypatch):
    d = tmp_path / "data" / "demo"
    d.mkdir(parents=True)
    (d / "u.bio").write_text(UNSTATED_PROGRAM)
    monkeypatch.setattr(evidence, "PROGRAM_DIRS", (Path("data/demo"),))
    evidence._cached.cache_clear()
    return tmp_path


def test_an_unstated_confidence_is_not_counted_as_a_weak_one(unstated_program):
    """Review R4 (2026-09-28) took `confidence:` off every predicted compiled fact and the parser read
    the gap as 0.0, so the explorer's weak count rose 812,921 -> 928,094 overnight on facts nobody had
    judged. A stated 0.0 is a judgement and stays weak; a missing one is counted apart."""
    out = evidence.collect(unstated_program)
    w = out["whole"]
    assert w["facts"] == 5
    assert (w["weak"], w["strong"], w["unstated"]) == (2, 1, 2)
    assert w["weak"] + w["strong"] + w["unstated"] == w["facts"]
    by = {r["label"].split(" =")[0]: r for r in out["rows"]}
    assert by["judged_zero"]["stated"] is True and by["judged_zero"]["confidence"] == 0.0
    assert by["NANOG"]["stated"] is False and by["NANOG activates NANOG"]["stated"] is False
    # the mean is over the three stated values, not dragged down by the two absences
    assert w["mean_confidence"] == round((0.0 + 0.4 + 0.9) / 3, 3)
    assert w["bands"]["unstated"] == 2
    assert w["mean_confidence_by_evidence"]["predicted"] == {
        "facts": 1,
        "stated": 0,
        "unstated": 1,
        "mean": None,
    }
    assert out["files"][0]["unstated"] == 2 and out["files"][0]["weak"] == 2


def test_the_weak_filter_and_the_review_list_keep_unstated_apart(unstated_program):
    weak = evidence.collect(unstated_program, max_confidence=0.5)["rows"]
    assert [r["label"].split(" =")[0] for r in weak] == ["judged_zero", "judged_weak"]
    rows = evidence.collect(unstated_program)["rows"]
    assert [r["stated"] for r in rows] == [True, True, True, False, False]  # judged first
    lines = evidence.to_csv(rows).strip().splitlines()
    assert lines[1].startswith("0.0,inferred,parameter,judged_zero")
    assert sum(1 for x in lines[1:] if x.startswith("unstated,")) == 2


def test_the_cli_headline_splits_weak_strong_and_unstated(unstated_program, monkeypatch, capsys):
    from genomeos.cli import main

    monkeypatch.chdir(unstated_program)
    assert main(["evidence", "--by-program"]) == 0
    out = capsys.readouterr().out
    assert "confidence stated on 3: 2 at or below 0.5, 1 above" in out and "2 state none" in out
    assert "unstated" in out.splitlines()[2]


def test_the_parser_marks_a_missing_confidence_and_keeps_its_value_zero():
    from genomeos.ir import UNSTATED, confidence_stated
    from genomeos.lang import parse

    m = parse(UNSTATED_PROGRAM)
    assert not confidence_stated(m.entities["NANOG"].confidence)
    assert m.entities["NANOG"].confidence == 0.0 and f"{m.rules[0].confidence:.2f}" == "0.00"
    assert confidence_stated(m.parameters["judged_zero"].confidence)
    assert m.parameters["judged_zero"].confidence == 0.0
    assert not confidence_stated(UNSTATED) and confidence_stated(0.0)


def test_the_distinction_survives_bioir_json_and_the_explorer_reads_it_back():
    """Item 12 S1 (2026-09-28) replaces the test that accepted the loss: BioIR JSON wrote an unstated
    confidence as 0.0 and from_dict read it back as stated. It now writes null, reads null back as
    UNSTATED, and keeps a stated 0.0 a number; the explorer's rows from the reloaded module say
    exactly what the rows from the parsed program say."""
    import json

    from genomeos.ir import Module, confidence_stated
    from genomeos.lang import parse

    m = parse(UNSTATED_PROGRAM)
    text = json.dumps(m.to_dict())
    d = json.loads(text)
    assert d["records_unstated_confidence"] is True
    by_id = {e["id"]: e for e in d["entities"]}
    assert by_id["NANOG"]["confidence"] is None and d["rules"][0]["confidence"] is None
    assert {p["name"]: p["confidence"] for p in d["parameters"]}["judged_zero"] == 0.0
    back = Module.from_dict(json.loads(text))
    assert confidence_stated(back.parameters["judged_zero"].confidence)
    assert back.parameters["judged_zero"].confidence == 0.0
    assert (
        not confidence_stated(back.entities["NANOG"].confidence) and back.entities["NANOG"].confidence == 0.0
    )
    assert not confidence_stated(back.rules[0].confidence)
    rows = [(r["label"], r["stated"], r["confidence"]) for r in evidence._rows_of(m, "u.bio")]
    assert [(r["label"], r["stated"], r["confidence"]) for r in evidence._rows_of(back, "u.bio")] == rows
    assert evidence.summarise(evidence._rows_of(back, "u.bio")) == evidence.summarise(
        evidence._rows_of(m, "u.bio")
    )
