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
