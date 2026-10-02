# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register the number-provenance census before anything is counted
(data/results/number_provenance_registration.json).

    uv run --frozen python scripts/number_provenance_register.py

It reads nothing but its own code and the digests of the programs the census will read. The
population, the scoped fields, the classes, the cascade, the convention tags, the citation-furniture
strip, what will be reported and the denominator all come from `genomeos.lang.number_provenance`, so
the registration and the code that applies it cannot drift apart.

No count is taken here, no class is assigned, no program is written to and no model request is made.
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.lang import number_provenance as np  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "number_provenance_registration"
WILL_WRITE = "data/results/number_provenance_census.json"

#: The measurements this census follows from, by the commits that made them.
DESCRIBES = (
    "d340e5c - the gastrulation program has no mesoderm at all because nothing in it represses "
    "SOX17; SOX17's repression factor is exactly 1.00000",
    "53d9936 / b1e8839 - three different sourced repressors of SOX17 restore a band and all three "
    "agree to the digit, because the module's own `threshold 2.0; hill 3` makes any added inhibitor "
    "a switch and not a brake",
    "5cbce26 - the Tbxt-SOX17 sign correction, which stands",
    "8bb9123 / 19423bd - tuning two free numbers to the unsourced expected_* proportions, refused",
)

#: This lane's own files.
OWN_CODE = (
    "genomeos/lang/number_provenance.py",
    "scripts/number_provenance_register.py",
    "scripts/number_provenance_census.py",
    "tests/test_number_provenance.py",
)


def inputs() -> list[dict[str, Any]]:
    """The programs the census will read, each digested before any count."""
    out: list[dict[str, Any]] = []
    for p in np.programs():
        out.append(mf.input_entry(p, partition=f"hand-authored BioLang program in {p.parent.as_posix()}"))
    return out


def main() -> None:
    payload = np.registration()
    payload["result"] = RESULT
    payload["date"] = date.today().isoformat()
    payload["lane"] = "lane-qualbound"
    payload["describes"] = list(DESCRIBES)
    payload["what_the_run_will_write"] = WILL_WRITE
    entries = inputs()
    payload["frozen_inputs"] = [
        {k: e[k] for k in ("path", "label", "sha256", "bytes") if k in e} for e in entries
    ]
    # The population figure is pinned HERE, as a registered datum with every program digested by
    # name, and not asserted in the test suite. tests/test_number_provenance.py asserts the
    # mechanism instead, because `scripts/package_engine.py` ships `genomeos/lang` - where
    # `number_provenance` lives - together with only part of the corpus, so a file count is a
    # property of this tree and not of the code.
    payload["population_size"] = len(entries)
    payload["population_programs"] = [e["path"] for e in entries]
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "the hand-authored BioLang corpus of this repository: "
                + ", ".join(np.POPULATION_DIRS),
                "version": "the working tree at the code revision stamped below; every program "
                "digested by its own path in inputs",
            },
            {
                "accession": "genomeos/ir/model.py dataclass defaults, read for the value the "
                "language supplies when a field is omitted",
                "version": "the same code revision; the values are transcribed into "
                "number_provenance.SCOPE and a test asserts they still agree",
            },
        ],
        "inputs": entries,
        "assembly": "n/a: no genomic interval is read. Genome coordinates are out of scope and the "
        "reason is in out_of_scope_fields.",
        "coordinates": "n/a: the census addresses declarations by file and line, not by genomic position.",
        "parameters": {
            "population_directories": list(np.POPULATION_DIRS),
            "excluded_corpora": dict(np.EXCLUDED_CORPORA),
            "scoped_fields": payload["scoped_fields"],
            "tiers": payload["tiers"],
            "classes": dict(np.CLASSES),
            "cascade": list(np.CASCADE),
            "convention_tags": dict(np.CONVENTION_TAGS),
            "repeated_in_file_floor": np.REPEATED_IN_FILE_FLOOR,
            "citation_furniture_stripped": payload["citation_furniture_stripped"],
            "comments_are_context_only": payload["comments_are_context_only"],
            "s_derived_procedure": payload["s_derived_procedure"],
            "denominator": payload["denominator"],
            "reports": list(np.REPORTS),
            "model_requests": 0,
            "money_spent": "none; no paid API is called. Web fetches of cited sources are free and "
            "are reported per source with what each one does and does not state about a number.",
        },
        "exclusions": list(np.EXCLUSIONS),
        "partitions": {cls: why for cls, why in np.CLASSES.items()},
        "code_cleanliness": mf.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
    }
    path = save_result(RESULT, payload)
    print(f"{RESULT}: {path}")
    print(f"  population: {len(entries)} hand-authored programs in {', '.join(np.POPULATION_DIRS)}")
    print(
        f"  scoped fields: {len(np.SCOPE)} (kind, field) pairs over "
        f"{len({s.kind for s in np.SCOPE})} block kinds"
    )
    print(f"  classes, fixed now: {', '.join(np.CLASSES)}")
    print(
        f"  convention tags: {', '.join(np.CONVENTION_TAGS)}; repeated-in-file floor "
        f"{np.REPEATED_IN_FILE_FLOOR}"
    )
    print(f"  S_derived adjudications in the table at registration: {len(np.DERIVATIONS)}")
    print(f"  denominator: {payload['denominator']}")
    print(f"  the run will write: {WILL_WRITE}")
    print("  no slot classified, no number changed, no program written to, 0 model requests")


if __name__ == "__main__":
    main()
