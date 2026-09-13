import json
import os
import subprocess
import sys

from genomeos import jobs


def test_start_does_not_duplicate_a_live_job(tmp_path, monkeypatch):
    monkeypatch.setattr(jobs, "JOBS_DIR", tmp_path)
    monkeypatch.setattr(jobs, "_running", {})
    sleeper = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])
    try:
        (tmp_path / "distil.json").write_text(
            json.dumps({"pid": sleeper.pid, "started": 1, "finished": None, "code": None})
        )
        st = jobs.start("distil", tmp_path)
        assert st.state == "running" and not (tmp_path / "distil.log").exists()  # no twin started
    finally:
        sleeper.kill()
        sleeper.wait()
    assert jobs.status("distil", tmp_path).state == "unknown"  # gone, no exit code recorded


def test_alive_handles_missing_pids():
    assert not jobs._alive(None)
    assert jobs._alive(os.getpid())


def test_stalled_and_heal(tmp_path, monkeypatch):
    import time

    monkeypatch.setattr(jobs, "JOBS_DIR", tmp_path)
    monkeypatch.setattr(jobs, "_running", {})
    monkeypatch.setattr(jobs, "STALL_AFTER", 1)
    root = tmp_path
    (root / "data" / "jobs").mkdir(parents=True)
    (root / "data" / "results").mkdir(parents=True)
    # a live process recorded for the proteome job, but no sign of life for longer than STALL_AFTER
    sleeper = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])
    try:
        (tmp_path / "proteome_genome_wide.json").write_text(
            json.dumps({"pid": sleeper.pid, "started": 1, "finished": None, "code": None})
        )
        log = root / "data" / "jobs" / "proteome_genome_wide.log"
        log.write_text("chr19: 10/100 compiled (5 s)\n")
        old = time.time() - 10
        os.utime(log, (old, old))
        st = jobs.status("proteome_genome_wide", root)
        assert st.state == "stalled" and st.activity is not None and st.activity >= 9
        assert st.detail == "chr19: 10/100" and st.done == 0.1
        # a heartbeat brings it back to running
        jobs.heartbeat("proteome_genome_wide", root)
        assert jobs.status("proteome_genome_wide", root).state == "running"
        assert not jobs.is_complete("proteome_genome_wide", root)
        # a completed job is never healed
        for i in range(25):
            (root / "data" / "results" / f"proteome_chr{i}.json").write_text("{}")
        assert jobs.is_complete("proteome_genome_wide", root)
        assert jobs.supervise_once(root) == []
        sleeper.kill()
        sleeper.wait()
        assert jobs.status("proteome_genome_wide", root).state == "done"  # complete, process gone
    finally:
        sleeper.kill()
        sleeper.wait()
