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


def test_the_pacer_holds_every_worker_and_escalates():
    import time

    p = eta.Pacer(patience=0.05, long_sleep=7200)
    held = p.hold(1)
    assert 0 < held <= 1 and p.waits == 1
    time.sleep(0.06)
    held = p.hold(1)  # still refusing after the patience: the long sleep takes over
    assert held > 3600 and p.waits == 2
    p.clear()
    assert p.since is None


def test_a_hang_is_not_a_quota_answer():
    assert eta.retry_seconds("DeadlineExceeded: no answer") is None
    assert eta.retry_seconds('RESOURCE_EXHAUSTED details = "Quota exceeded; retry in 30s"') == 30
