# SPDX-License-Identifier: AGPL-3.0-or-later
"""Adopting a job that was started outside the registry, so the Progress tab can see it."""

import json
import subprocess
import sys

import pytest

from genomeos import jobs


@pytest.fixture
def registry(tmp_path, monkeypatch):
    """A registry writing into a temporary directory, with one job that runs nothing."""
    monkeypatch.setattr(jobs, "JOBS_DIR", tmp_path / "jobs")
    (tmp_path / "jobs").mkdir()
    jobs.CATALOG["_adopt_test"] = {
        "argv": [sys.executable, "-c", "pass"],
        "describe": "t",
        "total": 1,
        "result": None,
        "count": None,
    }
    yield tmp_path
    jobs.CATALOG.pop("_adopt_test", None)
    jobs._running.pop("_adopt_test", None)


@pytest.fixture
def sleeper():
    """A live process to stand in for a job someone started by hand."""
    proc = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])
    yield proc
    proc.kill()
    proc.wait()


def test_an_adopted_job_is_running_and_recorded(registry, sleeper) -> None:
    st = jobs.adopt("_adopt_test", sleeper.pid, registry)

    assert st.state == "running"
    meta = json.loads((registry / "jobs" / "_adopt_test.json").read_text())
    assert meta["pid"] == sleeper.pid
    assert meta["adopted"] is True
    assert meta["finished"] is None


def test_a_dead_pid_is_refused(registry) -> None:
    """Adopting something that is not running would show work where there is none."""
    done = subprocess.Popen([sys.executable, "-c", "pass"])
    done.wait()

    with pytest.raises(RuntimeError, match="not running"):
        jobs.adopt("_adopt_test", done.pid, registry)


def test_a_second_live_process_is_refused(registry, sleeper) -> None:
    """Two processes under one job name would race, which is what the registry exists to prevent."""
    jobs.adopt("_adopt_test", sleeper.pid, registry)
    other = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])
    try:
        with pytest.raises(RuntimeError, match="already has a live process"):
            jobs.adopt("_adopt_test", other.pid, registry)
    finally:
        other.kill()
        other.wait()


def test_re_adopting_the_same_process_is_allowed(registry, sleeper) -> None:
    """Adopting twice is how a restarted session re-attaches; it must not be an error."""
    jobs.adopt("_adopt_test", sleeper.pid, registry)

    assert jobs.adopt("_adopt_test", sleeper.pid, registry).state == "running"


def test_an_unknown_job_name_is_refused(registry, sleeper) -> None:
    with pytest.raises(KeyError):
        jobs.adopt("_not_a_job", sleeper.pid, registry)


def test_the_panel_sweep_is_in_the_catalogue() -> None:
    """The sweep exists so the panel is startable and visible rather than run by hand."""
    spec = jobs.CATALOG["human_panel_sweep"]

    assert spec["total"] == 24
    assert spec["argv"][-1].endswith("human_panel_sweep.py")


def test_a_long_run_reports_progress_from_its_log_not_only_at_the_end(tmp_path) -> None:
    """A job that says nothing for an hour is a job the Progress tab cannot speak about.

    The VISTA run writes its result when it finishes, so until then the only honest source of progress
    is its own log. The result wins once it exists, because the log is the estimate and the result is
    the count.
    """
    from genomeos import jobs

    logs = tmp_path / "data" / "jobs"
    logs.mkdir(parents=True)
    (logs / "satmut_vista_run.log").write_text(
        "[06:11] 100/3200 windows scored\n[06:19] 500/3200 windows scored\n[06:25] 800/3200 windows scored\n"
    )

    assert jobs._windows_scored(tmp_path) == 800.0
    assert not jobs.CATALOG["satmut_vista"]["complete"](tmp_path)


def test_progress_is_zero_before_the_run_writes_anything(tmp_path) -> None:
    from genomeos import jobs

    assert jobs._windows_scored(tmp_path) == 0.0
