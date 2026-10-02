# SPDX-License-Identifier: AGPL-3.0-or-later
"""Each test is its own run, as far as the input tracer is concerned (2026-10-02, lane-tracer).

`genomeos/manifest.py` records every file opened under data/ from the moment it is imported, and
`save_result` reconciles a new registry result against what the open window holds. A pytest session is
one process, so without this fixture a window would span every test since the last `save_result`, and
a result written by one test would be answerable for a file a different test read. A window per test is
the scope that matches what a window means: the reads of one run.

The fixture reopens the window rather than switching the tracer off, so every test still runs with the
hook live and a test can open its own window inside this one.
"""

import extras_lock
import local_data
import pytest

from genomeos import manifest as mf


@pytest.fixture(autouse=True)
def trace_window_per_test():
    mf.trace_begin()
    yield
    mf.trace_begin()


# `requires_extra` (2026-10-02, lane-bare), registered in pyproject's `markers` and given its meaning
# here. A test that genuinely needs an optional extra says which one, and skips by that declaration
# where the extra is not installed. It must never pass because some other extra happened to drag the
# package in -- seven tests did exactly that on pandas, which eight locked packages pull in -- and never
# fail for want of it: tests/extras_lock.py and the `test-bare` CI job are the other two halves. An
# unknown extra name raises rather than skipping, because a marker naming an extra pyproject does not
# provide would otherwise skip, or run, for no stated reason.
# `needs_local_data` (2026-10-02, lane-raceandlocal), registered in pyproject's `markers` and given its
# meaning here and in tests/local_data.py. The amended acceptance rule: the verdict is a worktree of the
# committed tree with the git-ignored stores linked read-only, so a test that needs local data RUNS
# there; where the stores are genuinely absent -- CI, or a bare worktree -- it SKIPS BY NAME and never
# fails. By name means the reason identifies the missing store, so a reader can tell "not run here" from
# "passed". It never weakens an assertion: the marker moves where a test runs, not what it claims. A path
# git does not ignore raises, for the reason `requires_extra` raises on an unknown extra -- a marker on a
# committed path would skip, or run, for no stated reason, and a typo would read as a clean skip forever.
def pytest_runtest_setup(item):
    for mark in item.iter_markers("needs_local_data"):
        absent = local_data.missing(tuple(mark.args))
        if absent:
            pytest.skip(local_data.skip_reason(absent, mark.kwargs.get("how")))
    for mark in item.iter_markers("requires_extra"):
        for extra in mark.args:
            missing = extras_lock.missing_modules_for_extra(extra)
            if missing:
                pytest.skip(
                    f"needs the {extra!r} extra (uv sync --extra {extra}): cannot import {', '.join(missing)}"
                )
