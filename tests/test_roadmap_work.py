from pathlib import Path

import pytest

from genomeos import roadmap, work
from genomeos.web.server import Api

ROOT = Path(__file__).resolve().parent.parent

SAMPLE = """# Plan

Reviewed 2026-09-11.

## 3. Areas

### B. Genome decoding

- **Goal.** Every base accounted for.
- **Code.** `genome/segments.py`, `genome/{hic,mouse}.py`.
- **Missing.** Silencers have no source. Done 2026-09-11: the reader and
  the second mammal.
- **Enhancer to gene (2026-09-11).** Deleting an element moves a gene. The
  node holds 90.2% of them.

  | signals | precision |
  |---|---|
  | matrices | 5.2% |

  Kept side by side.
- **Next.** 1. The `reader` construct in BioLang, version 0.3 and
  22 autosomes. 2. Measured Hi-C, blocked on a key. 3. Done: the parser; next a
  second chromosome. 4. Done 2026-09-12: phasing.
- **Owner.** genomeos-8e (genomeos-fe before).

## 4. Data jobs

| Layer | Done | Job | Notes |
|---|---|---|---|
| Anatomy | 25 of 25 | done | fine |
| Parser | chr21 | `genomeos segments` | open |

| Job | Needs | Owner |
|---|---|---|
| Telomere | **done 2026-09-12**: read | genomeos-8e |
| AlphaGenome live | a key | Albert |

## 6. Milestones

| Version | Milestone | Proof |
|---|---|---|
| **0.9 whole genome** ✅ | every chromosome | reached |
| **1.3 the 98%** ◑ | tiers | budgeted |
| **2.0 BioLang standalone** | package | a program runs |
"""


def test_areas_split_planned_from_done():
    (b,) = roadmap.parse_areas(SAMPLE)
    assert b["letter"] == "B" and b["title"] == "Genome decoding"
    assert [s["state"] for s in b["next"]] == ["planned", "blocked", "partial", "done"]
    assert "22 autosomes" in b["next"][0]["text"]  # a count inside a step is not a step
    assert b["owners"] == ["genomeos-8e", "genomeos-fe"]
    titles = [d["title"] for d in b["done"]]
    assert titles == ["Done since the review", "Enhancer to gene"]
    assert b["missing"] == "Silencers have no source."
    enh = b["done"][1]
    assert enh["date"] == "2026-09-11" and enh["summary"] == "Deleting an element moves a gene."
    assert "| matrices | 5.2% |" in enh["body"].split("\n") and enh["body"].endswith("Kept side by side.")


def test_milestones_and_data_jobs():
    ms = {m["milestone"]: m["state"] for m in roadmap.parse_milestones(SAMPLE)}
    assert ms == {"0.9 whole genome": "done", "1.3 the 98%": "partial", "2.0 BioLang standalone": "planned"}
    jobs = {j["job"]: j["state"] for j in roadmap.parse_data_jobs(SAMPLE)}
    assert jobs == {"Anatomy": "done", "Parser": "planned", "Telomere": "done", "AlphaGenome live": "planned"}


def test_area_of_a_changed_file():
    areas = roadmap.parse_areas(SAMPLE)
    assert roadmap.area_of_path("genomeos/genome/hic.py", areas) == "B"
    assert roadmap.area_of_path("tests/test_segments.py", areas) == "B"
    assert roadmap.area_of_path("genomeos/attribution/closure.py", areas) == "I"  # by where it lives
    assert roadmap.area_of_path("README.md", areas) is None


def test_the_real_roadmap_parses():
    r = roadmap.load(ROOT)
    letters = [a["letter"] for a in r["areas"]]
    assert letters[:3] == ["A", "B", "C"] and len(letters) >= 10
    assert r["totals"]["done"] > 30 and r["milestones"] and r["data_jobs"]
    assert all("_text" not in a for a in r["areas"])


def test_work_board_lifecycle(tmp_path):
    e = work.start(
        tmp_path, "s-1", "the closure rerun", "i", ["a.py, b.py", "c/"], "the 69 candidates", now=1000.0
    )
    assert e["area"] == "I" and e["files"] == ["a.py", "b.py", "c/"] and e["state"] == "working"
    work.update(tmp_path, "s-1", note="half way", state="waiting", now=1100.0)
    work.start(tmp_path, "s-2", "older task", now=0.0)
    board = work.board(tmp_path, now=work.STALE_AFTER + 500.0)
    assert [x["who"] for x in board] == ["s-1", "s-2"]
    assert board[0]["state"] == "waiting" and board[0]["note"] == "half way"
    assert board[1]["state"] == "stale"
    work.done(tmp_path, "s-1", now=2000.0)
    assert work.board(tmp_path, now=2000.0 + 60)[-1]["state"] == "done"
    assert all(x["who"] != "s-1" for x in work.board(tmp_path, now=2000.0 + work.KEEP_DONE + 1))
    with pytest.raises(ValueError):
        work.update(tmp_path, "nobody", note="x")
    with pytest.raises(ValueError):
        work.start(tmp_path, "../escape", "x")


def test_api_work_and_roadmap():
    api = Api(ROOT)
    w = api.work()
    assert {"board", "jobs", "uncommitted", "commits", "areas"} <= set(w)
    assert w["areas"]["B"].startswith("Genome decoding")
    r = api.roadmap()
    assert r["areas"] and "totals" in r
