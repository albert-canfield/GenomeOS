# SPDX-License-Identifier: AGPL-3.0-or-later
"""Amendment 1 to the repression lane's registration: the label two inherited results must carry
(data/results/repress2_registration_amendment.json).

    uv run --frozen python scripts/repress_register_amendment.py

`data/results/context_evidence.json` and `data/results/context_evidence_baserate.json` are suspended
as not yet reproduced: both record grouped inputs through `manifest.files_entry`, which the old
`scripts/manifest_rebuild.py` resolved by path, so a group was reported unavailable with no file
opened and the result was excused rather than compared - 12 of 13 inputs on the first, standing for
873 files. The tool is fixed at `ebbded6` but those two manifests predate its `members` field, so they
cannot be checked from their manifests until their writer re-records them. The suspension is about
reproducibility and not about any known error.

This amendment records that label, says why the suspension does not reach a figure this lane reports,
and names the one result the lane does quote a figure from with its file and its scope. It is written
as its own result rather than into `repress2_registration.json`, so nothing already on the record is
rewritten: the registration stands as committed and this sits beside it. Both are read together.

No count is taken here, no measured pair is opened, no direction is read and no model request is made.
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import repress2 as rp  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "repress2_registration_amendment"
AMENDS = "repress2_registration"
#: The commit the registration this amends was made in. Named so the pair can be read in order.
AMENDS_COMMIT = "7a5dd8b"


def inputs() -> list[dict[str, Any]]:
    """The registration this amends and the two suspended results, each by its own path.

    The two suspended files are digested as they stand on this machine and under the label they must
    be quoted with. Nothing in this lane reads either of them; they are frozen here so that a reader
    can see which bytes carried the label and does not have to take it on trust.
    """
    out: list[dict[str, Any]] = []
    p = Path(f"data/results/{AMENDS}.json")
    if p.exists():
        out.append(mf.input_entry(p, partition=f"the registration this amends, committed in {AMENDS_COMMIT}"))
    for name, label in rp.SUSPENDED_INPUTS.items():
        q = Path(name)
        if q.exists():
            out.append(mf.input_entry(q, partition=label))
    q = Path(rp.QUOTED_FROM["file"])
    if q.exists():
        out.append(mf.input_entry(q, partition=rp.QUOTED_FROM["records_its_inputs"]))
    return out


def main() -> None:
    entries = inputs()
    payload: dict[str, Any] = {
        "result": RESULT,
        "date": date.today().isoformat(),
        "lane": "lane-repress2",
        "amends": f"data/results/{AMENDS}.json",
        "amends_commit": AMENDS_COMMIT,
        "status": (
            "a labelling amendment to inherited figures, registered before any count was taken. It "
            "changes no gate, no floor, no band, no mechanism and no refutation condition: the design "
            f"stands exactly as committed in {AMENDS_COMMIT}, and nothing in this payload is a result"
        ),
        "why_it_is_its_own_file": (
            "written beside the registration rather than into it, so that nothing already on the "
            "record is rewritten. The two are read together, and scripts/repress_population.py "
            "carries both into its own result through repress2.registration()"
        ),
        "suspended_inputs": dict(rp.SUSPENDED_INPUTS),
        "suspension_call": rp.SUSPENSION_CALL,
        "not_inherited_from_a_suspended_file": rp.NOT_INHERITED_FROM_A_SUSPENDED_FILE,
        "quoted_from": dict(rp.QUOTED_FROM),
        "what_does_not_change": [
            "the blinded eligibility count of gate 1 and the fields it may read",
            f"both floors: {rp.LOCUS_FLOOR} independent loci and {rp.POSITIVE_FLOOR} eligible "
            "measured links, imported and not set here",
            "the most-extreme-track mechanism, the prediction it makes and the four comparisons",
            "the refutation condition, word for word",
            "the three sets of reading words, and that nothing is labelled validation",
        ],
        "alphagenome_requests": 0,
        "money": "none: no request is made",
        "code": {
            "module": "genomeos/attribution/repress2.py",
            "amendment": "scripts/repress_register_amendment.py",
            "tests": "tests/test_repress2.py",
        },
    }
    payload["frozen_inputs"] = [{k: e[k] for k in ("path", "sha256", "bytes") if k in e} for e in entries]
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "this lane's own registration and the three results it stands beside",
                "version": "the bytes on this machine, digested in frozen_inputs",
            }
        ],
        "inputs": entries,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "locus_floor": rp.LOCUS_FLOOR,
            "link_floor": rp.POSITIVE_FLOOR,
            "aggregate_function": "none: nothing is aggregated",
            "fill_value": 0,
        },
        "exclusions": ["nothing is excluded: no population is counted here"],
        "partitions": {"amendment": "a label on inherited figures, not a reading of any population"},
    }
    path = save_result(RESULT, payload)
    print(f"{RESULT}: {path}")
    for name, label in rp.SUSPENDED_INPUTS.items():
        print(f"  {name}: {label}")
    print(f"  quoted with its scope: {rp.QUOTED_FROM['file']} ({rp.QUOTED_FROM['scope']})")
    print(f"  inputs digested: {len(entries)}; nothing counted, 0 model requests")


if __name__ == "__main__":
    main()
