# SPDX-License-Identifier: AGPL-3.0-or-later
"""The job registry's status against processes started in other places."""


def test_status_ignores_a_stale_handle_for_a_job_restarted_elsewhere(tmp_path, monkeypatch):
    """A server that started a job keeps a handle; when the job is started again from the CLI the handle
    points at the old, finished process and must not mark the new run done."""
    import json
    import subprocess
    import sys

    from genomeos import jobs

    monkeypatch.setattr(jobs, "JOBS_DIR", tmp_path / "jobs")
    (tmp_path / "jobs").mkdir()
    jobs.CATALOG["_stale_test"] = {
        "argv": [sys.executable, "-c", "pass"],
        "describe": "t",
        "total": 1,
        "result": None,
        "count": None,
    }
    try:
        old = subprocess.Popen([sys.executable, "-c", "pass"])
        old.wait()
        jobs._running["_stale_test"] = old  # the handle of a run that has ended
        live = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(5)"])
        try:
            (tmp_path / "jobs" / "_stale_test.json").write_text(
                json.dumps({"pid": live.pid, "started": 0, "finished": None, "code": None})
            )
            (tmp_path / "jobs" / "_stale_test.log").write_text("x")
            st = jobs.status("_stale_test", tmp_path)
            assert st.state == "running"
            assert json.loads((tmp_path / "jobs" / "_stale_test.json").read_text())["finished"] is None
            assert "_stale_test" not in jobs._running
        finally:
            live.kill()
    finally:
        jobs.CATALOG.pop("_stale_test", None)
        jobs._running.pop("_stale_test", None)
