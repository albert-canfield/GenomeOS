# SPDX-License-Identifier: AGPL-3.0-or-later
"""The whole-node scoring job's quota parsing and compact rows."""

import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location("eta", Path("scripts/enhancer_targets_all.py"))
eta = importlib.util.module_from_spec(spec)
spec.loader.exec_module(eta)


def test_retry_seconds_reads_the_quota_answer():
    assert (
        eta.retry_seconds('status = StatusCode.RESOURCE_EXHAUSTED details = "Quota exceeded; retry in 40s"')
        == 40
    )
    assert eta.retry_seconds("RESOURCE_EXHAUSTED") == 60
    assert eta.retry_seconds("ConnectionResetError: peer") is None


def test_compact_keeps_the_fields_the_closure_reads():
    r = {
        "id": "E",
        "start": 1,
        "end": 2,
        "predicted_by_cell": {"K562": -0.1},
        "genes": [1, 2],
        "top_genes": [],
    }
    c = eta.compact(r)
    assert c["predicted_by_cell"] == {"K562": -0.1} and "genes" not in c and c["inferred"] is None


def test_the_watchdog_cuts_off_a_request_that_never_returns():
    import signal
    import time

    signal.signal(signal.SIGALRM, eta._hang)
    signal.setitimer(signal.ITIMER_REAL, 0.05)
    try:
        try:
            time.sleep(5)
        except eta.HungError as ex:
            assert "no answer in" in str(ex)
        else:
            raise AssertionError("the alarm did not fire")
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    # a hang is not a quota answer, so the loop's generic branch retries it
    assert eta.retry_seconds("Hung: no answer in 300s") is None
