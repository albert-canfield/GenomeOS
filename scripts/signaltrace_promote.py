# SPDX-License-Identifier: AGPL-3.0-or-later
"""M1c: promote the registered run's quarantined result, its manifest corrected and nothing else changed.

    uv run python scripts/signaltrace_promote.py data/quarantine/results/signaltrace_result.json

THE DEFECT, reported beside the registered readings. The one registered run (scripts/signaltrace_run.py,
frozen at d21ad05, run from a clean worktree of the registration commit 9637584) computed every stage
and every reading, and `save_result` then refused to write it because the run had read its own
registration (data/results/signaltrace_registration.json) without declaring it among its inputs. The
payload went to the quarantine whole. This script does not recompute anything: it checks that the
quarantined payload was written by the registered code on a clean tree and that the undeclared
registration is its only manifest problem, adds the registration (and the quarantined file it reads)
to the inputs, keeps the run's own code stamp and traced reads under `promotion`, and writes the payload
with every computed field byte-identical, which it checks before writing. No stage, column, population,
metric or reading is changed; the run is not repeated.

No request is sent, no key is read. 0 model requests, no network, no money.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

RESULT = "signaltrace_result"
REGISTRATION = RESULTS_DIR / "signaltrace_registration.json"
REGISTRATION_COMMIT = "9637584698535a9e64ea96c2c0d389c4e2c90c70"
OWN_CODE = (
    "genomeos/attribution/signaltrace.py",
    "scripts/signaltrace_register.py",
    "scripts/signaltrace_run.py",
    "scripts/signaltrace_promote.py",
    "tests/test_signaltrace.py",
)
#: the manifest fields the run's writer gave; everything else in a quarantined manifest is a stamp
WRITER_FIELDS = (
    "sources",
    "inputs",
    "input_count",
    "assembly",
    "coordinates",
    "parameters",
    "exclusions",
    "partitions",
)
ONLY_PROBLEM = "data/results/signaltrace_registration.json"
DEFECT = (
    "the registered run read its own registration (data/results/signaltrace_registration.json) and did"
    " not declare it among its inputs, so save_result quarantined the complete payload; the payload is"
    " written here unchanged with the registration (and the quarantined file read to promote it) declared."
    " No computed field was changed and the run was not repeated (scripts/signaltrace_promote.py)"
)


def checked(q: dict[str, Any]) -> None:
    """Refuse unless the payload is the registered code's, on a clean tree, with only this one problem."""
    m = q.get(mf.KEY) or {}
    code = m.get("code") or {}
    problems = m.get("problems") or []
    if code.get("git_sha") != REGISTRATION_COMMIT or code.get("dirty") is not False:
        raise SystemExit(
            f"refused: written by {code.get('git_sha')} dirty={code.get('dirty')},"
            " not the clean registration commit"
        )
    if code.get("argv") != ["scripts/signaltrace_run.py"]:
        raise SystemExit(f"refused: written by {code.get('argv')}, not the registered run")
    if len(problems) != 1 or ONLY_PROBLEM not in problems[0] or "not declared" not in problems[0]:
        raise SystemExit(f"refused: the quarantined manifest's problems are not the one defect: {problems}")
    clean = m.get("code_cleanliness") or {}
    if not clean.get("own_code_is_committed") or clean.get("foreign_uncommitted_code_on_the_counting_path"):
        raise SystemExit(f"refused: the run's code was not clean: {clean}")


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        raise SystemExit(__doc__)
    src = Path(argv[1])
    q = json.loads(src.read_text())
    checked(q)
    reg = json.loads(REGISTRATION.read_text())
    if reg["result_manifest"]["code"]["git_sha"] != reg["frozen_code"]["commit"]:
        raise SystemExit("refused: the registration was not written at its frozen code")
    run_manifest = q[mf.KEY]
    body = {k: v for k, v in q.items() if k not in (mf.KEY, "quarantine", "result", "date")}
    before = json.dumps(body, sort_keys=True)
    manifest = {k: run_manifest[k] for k in WRITER_FIELDS if k in run_manifest}
    manifest["inputs"] = [
        *run_manifest["inputs"],
        mf.input_entry(REGISTRATION, partition=None, role="the registration the run checked itself against"),
        mf.input_entry(
            src, partition=None, role="the registered run's quarantined payload, promoted unchanged"
        ),
    ]
    manifest["input_count"] = len(manifest["inputs"])
    manifest["code_cleanliness"] = mf.code_cleanliness(__file__, OWN_CODE)
    body["promotion"] = {
        "defect_reported_beside_the_readings": DEFECT,
        "run_code": run_manifest.get("code"),
        "run_code_cleanliness": run_manifest.get("code_cleanliness"),
        "run_traced_inputs": run_manifest.get("traced_inputs"),
        "run_quarantine": q.get("quarantine"),
        "computed_fields_unchanged": True,
    }
    after = json.dumps({k: v for k, v in body.items() if k != "promotion"}, sort_keys=True)
    if after != before:
        raise SystemExit("refused: a computed field would change")
    print("saved", save_result(RESULT, body, manifest=manifest))
    print(json.dumps({"outcome": body["outcome"], "reading": body["reading"]}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
