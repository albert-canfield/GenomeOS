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


def test_a_dead_lane_is_abandoned_and_retired(tmp_path):
    """Two days untouched reads differently from six hours, and retiring moves rather than deletes.

    The board said "working" for eleven lanes on 2026-09-21 whose sessions had ended days before,
    because a restart renames a session and the old entry is never touched again. Stale is a pause;
    abandoned is nobody coming back, and only the second should leave the board.
    """
    work.start(tmp_path, "paused", "a lane between commits", now=0.0)
    work.start(tmp_path, "dead", "a lane whose session ended", now=0.0)
    work.update(tmp_path, "paused", note="still here", now=work.ABANDONED_AFTER)
    board = work.board(tmp_path, now=work.ABANDONED_AFTER + work.STALE_AFTER + 1)
    states = {x["who"]: x["state"] for x in board}
    assert states == {"paused": "stale", "dead": "abandoned"}

    named = work.retire(tmp_path, hours=48, now=work.ABANDONED_AFTER + 1, dry_run=True)
    assert named == ["dead"] and (tmp_path / "data" / "work" / "dead.json").exists()
    assert work.retire(tmp_path, hours=48, now=work.ABANDONED_AFTER + 1) == ["dead"]
    assert [x["who"] for x in work.board(tmp_path, now=work.ABANDONED_AFTER + 1)] == ["paused"]
    assert (tmp_path / "data" / "work" / "retired" / "dead.json").exists(), "moved, not deleted"


def test_api_work_and_roadmap():
    api = Api(ROOT)
    w = api.work()
    assert {"board", "jobs", "uncommitted", "commits", "areas"} <= set(w)
    assert w["areas"]["B"].startswith("Genome decoding")
    r = api.roadmap()
    assert r["areas"] and "totals" in r


def test_the_board_root_can_be_named_by_the_environment(tmp_path, monkeypatch):
    """A session in its own git worktree must be able to claim the SHARED board.

    `genomeos work` resolved the board as `Path.cwd()/data/work/`, and data/work/ is git-ignored, so a
    lane working in an isolated worktree -- which is now the rule for any module on every result's
    counting path -- wrote its claim to a board only it could read. On 2026-10-02 the first lane to work
    that way was invisible on Albert's Progress tab while running. GENOMEOS_WORK_ROOT names the checkout
    whose board to use, so the work stays isolated and the claim does not.

    It also makes probing safe: the same day, a check of which argument order parses was run against the
    live board and left an entry there that had to be removed by hand.
    """
    import subprocess
    import sys

    def run(*args, env_root=None):
        env = {"PATH": "/usr/bin:/bin:/usr/sbin:/sbin", "HOME": str(tmp_path)}
        if env_root is not None:
            env["GENOMEOS_WORK_ROOT"] = str(env_root)
        return subprocess.run(
            [sys.executable, "-c", "from genomeos.cli import main; raise SystemExit(main())", *args],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
            env={**env, "PYTHONPATH": str(ROOT)},
        )

    elsewhere = tmp_path / "a-worktree"
    elsewhere.mkdir()
    r = run("work", "start", "--who", "lane-in-a-worktree", "a task", env_root=elsewhere)
    assert r.returncode == 0, (r.stdout, r.stderr)
    # it went to the named root, and NOT to the live board
    assert (elsewhere / "data" / "work" / "lane-in-a-worktree.json").is_file(), r.stdout
    assert not (ROOT / "data" / "work" / "lane-in-a-worktree.json").exists(), "wrote to the live board"
    # and the listing from that root sees it while the live board does not
    assert "lane-in-a-worktree" in run("work", "list", env_root=elsewhere).stdout
    assert "lane-in-a-worktree" not in run("work", "list").stdout
