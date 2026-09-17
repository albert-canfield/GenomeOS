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
    order = lambda r: (r["confidence"], r["path"], r["block"], r["label"])  # noqa: E731
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
