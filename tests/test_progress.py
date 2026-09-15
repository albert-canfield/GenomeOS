from pathlib import Path

from genomeos import jobs
from genomeos.web.server import Api


def test_progress_rows_come_from_the_docs():
    p = Api(Path(".")).progress()
    assert p["count"] > 20 and p["done"] > 10
    assert all({"task", "status", "evidence"} <= set(r) for r in p["rows"])


def test_job_catalog_and_status():
    names = {j.name for j in jobs.all_status(Path("."))}
    assert {"anatomy_genome_wide", "signals_learn_chr21", "distil"} <= names
    st = jobs.status("anatomy_genome_wide", Path("."))
    assert (
        st.total == 25 and 0 <= st.done <= 25 and st.state in ("idle", "running", "done", "failed", "unknown")
    )
    d = st.to_dict()
    assert "fraction" in d and "last_lines" in d


def test_unknown_job_is_rejected():
    import pytest

    with pytest.raises(KeyError):
        jobs.start("rm_rf", Path("."))


def test_all_elements_job_reads_its_result_not_its_process(tmp_path, monkeypatch):
    """A chromosome the chain has finished reads done, with the element count its result records.

    The scorer is stopped by the chain rather than exiting on its own, so the job's metadata says
    "started, never finished" for a run that did complete. Before this, finished chromosomes showed
    as stalled and their totals came from whatever the result held when the module was imported.
    """
    import json
    from pathlib import Path

    from genomeos import jobs

    results = tmp_path / "data" / "results"
    results.mkdir(parents=True)
    (results / "enhancer_targets_all_chr21.json").write_text(
        json.dumps({"chrom": "chr21", "scored": 12139, "elements_total": 12139, "complete": True})
    )
    (results / "enhancer_targets_all_chr15.json").write_text(
        json.dumps({"chrom": "chr15", "scored": 3400, "elements_total": 31058, "complete": False})
    )
    meta = tmp_path / "data" / "jobs"
    meta.mkdir(parents=True)
    monkeypatch.setattr(jobs, "JOBS_DIR", meta)
    # a process the chain stopped mid-run: a pid far above any in use, so it reads as gone
    for name in ("enhancer_targets_all_chr21", "enhancer_targets_all_chr15"):
        (meta / f"{name}.json").write_text(json.dumps({"pid": 2**22, "started": 1.0, "finished": None}))

    done = jobs.status("enhancer_targets_all_chr21", Path(tmp_path))
    assert done.state == "done" and done.done == 12139 and done.total == 12139
    part = jobs.status("enhancer_targets_all_chr15", Path(tmp_path))
    assert part.state != "done" and part.done == 3400 and part.total == 31058  # total read per call
    assert jobs.CATALOG["enhancer_targets_all_chr15"]["auto_heal"] is False  # the chain sequences these


def test_two_all_elements_jobs_cannot_hold_the_key_at_once(tmp_path, monkeypatch):
    """One AlphaGenome key, one scorer: the registry refuses the second rather than racing.

    A stray chr1 run was started from the Progress tab (or a supervisor tick) four times in one day
    while the chain was on another chromosome; the chain killed each one, but not before it had spent
    requests. Refusing at the registry costs nothing and names the chain as the way to run them.
    """
    import json
    import os

    import pytest

    from genomeos import jobs

    meta = tmp_path / "data" / "jobs"
    meta.mkdir(parents=True)
    monkeypatch.setattr(jobs, "JOBS_DIR", meta)
    # chr13 is scoring: its recorded pid is this test process, which is certainly alive
    (meta / "enhancer_targets_all_chr13.json").write_text(
        json.dumps({"pid": os.getpid(), "started": 1.0, "finished": None, "code": None})
    )
    assert jobs._key_held_by("enhancer_targets_all_chr1") == "enhancer_targets_all_chr13"
    with pytest.raises(RuntimeError, match="shares one AlphaGenome key"):
        jobs.start("enhancer_targets_all_chr1", tmp_path)
    assert not (meta / "enhancer_targets_all_chr1.json").exists()  # nothing was started
    # a job that does not touch the key is unaffected
    assert jobs._key_held_by("anatomy_genome_wide") is None


def test_a_job_outside_the_registry_can_hold_the_key(tmp_path, monkeypatch):
    """The lock covers any spender, not only registry jobs.

    A Start from the Progress tab during an executor run shared the quota silently for six minutes,
    because the registry could only see its own jobs. Now the executor takes the lock and the
    registry refuses; a holder whose process has died releases it by itself.
    """
    import json
    import os

    import pytest

    from genomeos import jobs

    meta = tmp_path / "data" / "jobs"
    meta.mkdir(parents=True)
    monkeypatch.setattr(jobs, "JOBS_DIR", meta)

    jobs.take_key("executor E1W", "1,000 pairs")
    held = jobs.key_holder()
    assert held["holder"] == "executor E1W" and held["pid"] == os.getpid()
    with pytest.raises(RuntimeError, match=r"shares one AlphaGenome key with executor E1W.*wait for it"):
        jobs.start("enhancer_targets_all_chr1", tmp_path)

    jobs.drop_key("executor E1W")
    assert jobs.key_holder() is None
    assert jobs._key_held_by("enhancer_targets_all_chr1") is None

    # a lock left behind by a process that died is stale, and nobody waits on a corpse
    (meta / f"{jobs.KEY_LOCK}.lock").write_text(
        json.dumps({"holder": "gone", "what": "killed session", "pid": 2**22, "started": 1.0})
    )
    assert jobs.key_holder() is None
