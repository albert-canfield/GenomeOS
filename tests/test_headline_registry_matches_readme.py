"""The headline registry and the README must name the same results, in both directions.

Why this file exists, and it is an incident rather than a worry. A promotion gate was aimed at the six
results in `manifest_headlines.json` because that file is the record of what the README's claims rest on.
It is dated 2026-09-28 and the README has moved since: the README's independent-prediction claim now cites
`crispri_published_v2.json` and `crispri_benchmark_v2.json`, and the registry names NEITHER. So the gate
spent a rebuild on `crispri_published.json`, which no README claim depends on any more, and did not look at
two results that one does. A stale registry is how a verification effort ends up verifying something nobody
relies on, which is worse than not verifying at all, because it reports a pass.

The registry's `quoted_in` field is the other half and it drifted too. `node_containment_audit` records
`quoted_in: "README; docs/NODES-READER-WRITER.md; docs/LESSONS.md"` and `quoted_as: "+2.90"`, and the
string `2.90` does not appear in README.md at all. That entry is still relied on -- both docs do quote it --
so the defect is in the registry's claim about the README and not in the result.

The figure test matches the NUMERIC TOKENS of `quoted_as` rather than the whole string, deliberately. The
registry writes `unknown_coverage` as `"160,447 of 30,602,182 bp (0.52%)"` and the README writes "of the
30,602,182 bases of the real unknown, 160,447 (0.52%)". Those are the same claim in different prose, and a
whole-string test would call it drift and teach its reader to skim the failures. Matching the numbers keeps
the test firing only when a quantity is gone.

All three tests failed when this file was written and were strict xfails with the measured counts in the
reason, because the registry and the README are other sessions' files in a shared checkout: this file
records the defect so it cannot be forgotten and flips to a pass the moment either is corrected. It never
adjusts what it expects.

The FIRST of the three was corrected on 2026-10-02 and its marker is gone, which is the mechanism working
rather than an adjustment: `crispri_published_v2` and `crispri_benchmark_v2` were rebuilt and registered in
`manifest_headlines.json`, so the README-to-registry direction holds and the strict xfail would otherwise
XPASS and redden the suite. The other two still fail and keep their markers: `constrained_unknown_targets`
still records "ROADMAP 1.3", and `node_containment_audit` still records `quoted_in: "README"` for a figure
of +2.90 that the README does not contain.
"""

import json
import re
from pathlib import Path

import pytest

RESULTS = Path("data/results")
README = Path("README.md")
RECORD = json.loads((RESULTS / "manifest_headlines.json").read_text())
ENTRIES = RECORD["rebuilt"] + RECORD.get("pending", [])
REGISTERED = {e["result"] for e in ENTRIES}


def readme_result_filenames() -> set[str]:
    """Every `*.json` the README names that is a result file on disk. A name the README mentions which is
    not in data/results is not a result citation and is not this test's business."""
    return {n for n in re.findall(r"[a-z0-9_]+\.json", README.read_text()) if (RESULTS / n).is_file()}


def documents(quoted_in: str) -> list[str]:
    """The file paths a `quoted_in` field names. An entry may name a place that is not a path at all --
    `constrained_unknown_targets` says "docs/ATTRIBUTION.md; ROADMAP 1.3" -- and such a reference is
    returned as given so the test can say it cannot be checked rather than skip it silently."""
    out = []
    for part in (quoted_in or "").split(";"):
        part = part.strip()
        if not part:
            continue
        out.append("README.md" if part.startswith("README") else part)
    return out


def test_every_result_the_readme_names_by_filename_is_registered():
    unregistered = sorted(n for n in readme_result_filenames() if n.removesuffix(".json") not in REGISTERED)
    assert unregistered == [], (
        f"the README rests on these and the registry does not name them: {unregistered}"
    )


@pytest.mark.xfail(
    strict=True,
    reason="constrained_unknown_targets records 'ROADMAP 1.3', which is not a path any test can open",
)
def test_every_registered_entry_names_documents_that_exist():
    missing = [
        f"{e['result']}: {d}" for e in ENTRIES for d in documents(e.get("quoted_in")) if not Path(d).is_file()
    ]
    assert missing == [], f"a registry entry quotes a document that cannot be opened: {missing}"


@pytest.mark.xfail(
    strict=True,
    reason="node_containment_audit records quoted_in 'README' for +2.90 and the README does not contain 2.90",
)
def test_every_registered_figure_is_still_quoted_where_the_registry_says_it_is():
    gone = []
    for e in ENTRIES:
        numbers = re.findall(r"\d[\d,.]*", e.get("quoted_as") or "")
        for d in documents(e.get("quoted_in")):
            path = Path(d)
            if not path.is_file():
                continue  # the test above owns that failure
            text = path.read_text()
            absent = [n for n in numbers if n not in text]
            if absent:
                gone.append(f"{e['result']} -> {d}: {absent}")
    assert gone == [], f"the registry says these figures are quoted there and they are not: {gone}"
