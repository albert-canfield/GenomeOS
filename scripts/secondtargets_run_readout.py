# SPDX-License-Identifier: AGPL-3.0-or-later
"""The registered second-targets run on Gasperini 2019's fetched readout, under amendment 1.

    uv run python scripts/secondtargets_run_readout.py --amend  # the amendment: header and ids only
    uv run python scripts/secondtargets_run_readout.py          # the one run, after the amendment

The registration (data/results/secondtargets_registration.json, 18da8cb) was written while the screen's
differential-expression table was not on this machine, and its run (d4dacb3) stopped as
NOT_RUN_NO_READOUT. The coordinator then cleared one fetch of that table. Its header, read before any
value, does not carry the fold-change column the registration recalled ('fold_change.transcript_remap'):
it is 'fold_change.transcript_remaining'. And its pair ids write a candidate enhancer's guide group in
upper case ('CHR1.10703_TOP_TWO:ENSG...') where the pair table writes it in lower case. `--amend` writes
both facts, and the join that follows from them, as an additive amendment before any value is read; the
registration itself is not edited.

The run then does exactly what scripts/secondtargets_run.py does after its NOT_RUN gate, with the
amended column and join: the registration and its frozen code are checked (`run.check_registration`),
the amendment and this file's code are checked, the label-free census is recomputed and must equal the
registered one, the TSS positive control is read first, then the matched tests. The result is written to
a NEW file, data/results/secondtargets_result_run.json, which supersedes the NOT_RUN record; that record
stays as written.

No model request, no key, no money; the one network fetch was made by hand before this file was run and
is recorded in the amendment (URL, bytes, sha256, terms).
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import subprocess
import sys
import time
from collections import Counter
from collections.abc import Iterable
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import secondtargets as s2  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import secondtargets_run as run  # noqa: E402

RESULT = "secondtargets_result_run"
AMENDMENT_NAME = "secondtargets_registration_amendment_1"
AMENDMENT = RESULTS_DIR / f"{AMENDMENT_NAME}.json"
SUPERSEDES = "data/results/secondtargets_result.json"
REGISTRATION_COMMIT = "18da8cb"
OWN_CODE = (*run.OWN_CODE, "scripts/secondtargets_run_readout.py", "tests/test_secondtargets_readout.py")
#: this file, pinned by sha256 in the amendment
AMENDED_CODE = ("scripts/secondtargets_run_readout.py",)

#: the fetch the coordinator cleared on 2026-10-10, made by hand and recorded here
FETCH = {
    "url": s2.READOUT_URL,
    "path": str(s2.READOUT),
    "bytes": 37256050,
    "sha256": "b5da4650330b3afbe61209b12cdf952387e591f7a4a99a5e0ec9bb30e6d98560",
    "fetched_utc": "2026-10-10T13:31Z",
    "how": "curl -f -L to a .part file, gzip -t, then renamed; size matched the server's HEAD (37,256,050)",
    "terms": (
        "NCBI GEO public record GSE120861; the submitters state no licence. NCBI's GEO disclaimer places no"
        " restriction of its own on use or redistribution but notes submitters may claim rights. Cite"
        " Gasperini et al. 2019, Cell 176(1-2):377-390. Kept git-ignored under"
        " data/knowledge/measured_unknown"
    ),
    "authorised_by": "the coordinator (genomeos-99), 2026-10-10, one file only",
}
#: the amended hit columns: the registered rule with the table's real fold-change column name
HIT_AMENDED = {
    **s2.HIT,
    "pair_column": "pairs4merge",
    "fold_column": "fold_change.transcript_remaining",
}
JOIN_AMENDED = (
    "a pair id is compared after upper-casing its guide-group part and removing any Ensembl version from"
    " its gene part, on both sides: the readout table writes candidate-enhancer guide groups in upper case"
    " ('CHR1.10703_TOP_TWO:ENSG...') and the pair table in lower case ('chr1.10703_top_two:ENSG...'); TSS"
    " guide groups are written alike in both"
)
WHAT_WAS_READ = (
    "before this amendment was written, of the fetched table: its header line (20 columns), and its"
    " pairs4merge column alone (identifiers: their pattern with digits masked, whether any carries an"
    " Ensembl version, and how many of the registered ids it holds before and after the case change). No"
    " value column (beta, intercept, fold change, any p-value, outlier_gene) was read"
)
DEPARTURE = (
    "the registration's STOP_RULE allows an amendment that changes column names only. The fold-change"
    " column is such a change. The join is not: without it no candidate-enhancer test joins at all (0 of"
    " the 872 registered ids are in the table as written, 345 after the case change), and it was decided"
    " from identifiers alone, before any value was read. It is stated here so the reader can weigh it;"
    " nothing else registered changes"
)


def norm(key: str) -> str:
    """JOIN_AMENDED: 'chr1.10703_top_two:ENSG00000143862.4' -> 'CHR1.10703_TOP_TWO:ENSG00000143862'."""
    group, _, gene = key.strip('"').rpartition(":")
    return f"{group.upper()}:{s2.ensg_of(gene)}"


def header(path: Path = s2.READOUT) -> list[str]:
    with gzip.open(path, "rt") as fh:
        return fh.readline().rstrip("\n").split("\t")


def missing_columns(head: list[str], hit: dict[str, Any]) -> list[str]:
    return [c for c in (hit["pair_column"], hit["p_column"], hit["fold_column"]) if c not in head]


def identifiers(lines: Iterable[str], pair_column: str) -> set[str]:
    """The table's pair ids as written; only the pair column is split off each row."""
    it = iter(lines)
    head = next(it, "").rstrip("\n").split("\t")
    i = head.index(pair_column)
    return {ln.split("\t", i + 1)[i].strip('"') for ln in it if ln.strip()}


def read_amended(lines: Iterable[str], wanted: dict[str, list[str]]) -> dict[str, bool | None]:
    """HIT_RULE with the amended column and join. `wanted` maps a normalised id to the registered ids it
    stands for; the result is keyed by the registered ids. KeyError before any value if a column lacks."""
    it = iter(lines)
    head = next(it, "").rstrip("\n").split("\t")
    lacking = missing_columns(head, HIT_AMENDED)
    if lacking:
        raise KeyError(f"readout header lacks {lacking}")
    ik, ip, ifc = (head.index(HIT_AMENDED[c]) for c in ("pair_column", "p_column", "fold_column"))
    out: dict[str, bool | None] = {}
    for line in it:
        f = line.rstrip("\n").split("\t")
        if len(f) <= max(ik, ip, ifc):
            continue
        keys = wanted.get(norm(f[ik]))
        if not keys:
            continue
        hit: bool | None
        try:
            p, fold = float(f[ip]), float(f[ifc])
        except ValueError:
            hit = None
        else:
            measured = not (math.isnan(p) or math.isnan(fold))
            hit = (p < HIT_AMENDED["p_below"] and fold < HIT_AMENDED["fold_below"]) if measured else None
        for k in keys:
            out[k] = hit
    return out


def wanted_of(keys: Iterable[str]) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for k in sorted(set(keys)):
        out.setdefault(norm(k), []).append(k)
    return out


def code_sha256(paths: tuple[str, ...]) -> dict[str, str]:
    return {p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths}


def git_head() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True
    ).stdout.strip()


def manifest(entries: list[dict[str, Any]]) -> dict[str, Any]:
    base = run.manifest(entries, own_entry=__file__)
    readout_source = {
        "accession": "GEO GSE120861 (Gasperini et al. 2019, Cell), at-scale differential-expression table",
        "version": f"fetched {FETCH['fetched_utc']}; sha256 {FETCH['sha256']}; {FETCH['terms']}",
        "url": FETCH["url"],
    }
    return {
        **base,
        "sources": [*base["sources"], readout_source],
        "parameters": {**base["parameters"], "hit": HIT_AMENDED, "join": JOIN_AMENDED},
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def amend() -> int:
    """Write the amendment from the header and the identifiers alone."""
    if AMENDMENT.exists():
        raise SystemExit("refused: the amendment has already been written")
    reg = run.check_registration()
    got = mf.input_entry(s2.READOUT, partition="heldout", role="the screen's readout: header and ids only")
    if got["sha256"] != FETCH["sha256"] or got["bytes"] != FETCH["bytes"]:
        raise SystemExit(f"refused: the readout on disk is not the recorded fetch: {got}")
    head = header()
    registered_missing = missing_columns(head, {**s2.HIT, "pair_column": "pairs4merge"})
    if missing_columns(head, HIT_AMENDED):
        raise SystemExit(f"refused: even the amended columns are not in the header: {head}")
    with gzip.open(s2.READOUT, "rt") as fh:
        raw = identifiers(fh, HIT_AMENDED["pair_column"])
    normalised = {norm(k) for k in raw}
    pop = reg["population_ids"]
    reg_path = run.REGISTRATION
    body = {
        "lane": s2.LANE,
        "date": date.today().isoformat(),
        "amends": {
            "registration": str(reg_path),
            "sha256": hashlib.sha256(reg_path.read_bytes()).hexdigest(),
            "registered_in_commit": REGISTRATION_COMMIT,
            "frozen_code_commit": reg["frozen_code"]["commit"],
        },
        "fetch": FETCH,
        "header": head,
        "registered_columns_missing_from_header": registered_missing,
        "hit_registered": s2.HIT,
        "hit_amended": HIT_AMENDED,
        "join_amended": JOIN_AMENDED,
        "identifier_coverage": {
            "readout_distinct_ids": len(raw),
            "registered_matched_ids": len(pop),
            "registered_ids_as_written": sum(1 for k in pop if k in raw),
            "registered_ids_after_the_case_change": sum(1 for k in pop if norm(k) in normalised),
        },
        "what_was_read": WHAT_WAS_READ,
        "departure": DEPARTURE,
        "unchanged": (
            "everything else registered: the population, the shortlist, the matching, the thresholds (p below"
            " 0.1, fold below 1), the treatment of unmeasured tests, the positive control, MIN_MATCHED, the"
            " bootstrap, the decision, the falsifier and every outcome sentence"
        ),
        "result_will_be_written_to": f"data/results/{RESULT}.json",
        "supersedes_when_run": SUPERSEDES,
        "amended_code": {"commit": git_head(), "sha256": code_sha256(AMENDED_CODE)},
        "readout_values_read": 0,
        "alphagenome_requests": 0,
        "money": "none",
        "result_manifest": manifest([got, run.registration_entry()]),
    }
    print("saved", save_result(AMENDMENT_NAME, body))
    print(json.dumps({k: body[k] for k in ("registered_columns_missing_from_header", "identifier_coverage")}))
    return 0


def check_amendment() -> dict[str, Any]:
    am = json.loads(AMENDMENT.read_text())
    if am["amends"]["sha256"] != hashlib.sha256(run.REGISTRATION.read_bytes()).hexdigest():
        raise SystemExit("refused: the amendment does not amend the registration on disk")
    if am["hit_amended"] != json.loads(json.dumps(HIT_AMENDED)) or am["join_amended"] != JOIN_AMENDED:
        raise SystemExit("refused: the amendment and this file disagree on the hit columns or the join")
    if am["amended_code"]["sha256"] != code_sha256(AMENDED_CODE):
        raise SystemExit("refused: this file differs from what the amendment pinned")
    if am["fetch"] != json.loads(json.dumps(FETCH)):
        raise SystemExit("refused: the amendment and this file disagree on the fetch")
    return am


def measured_by_arm(units: list[s2.Unit], hits: dict[str, bool | None]) -> dict[str, int]:
    out: Counter = Counter()
    for u in units:
        out["shortlist"] += hits.get(u.test.test.key) is not None
        out["controls"] += sum(hits.get(c.test.key) is not None for c in u.controls)
    return dict(out)


def main() -> int:
    t0 = time.time()
    if (RESULTS_DIR / f"{RESULT}.json").exists():
        raise SystemExit("refused: the registered run has already been written")
    reg = run.check_registration()
    am = check_amendment()
    built = run.build()
    census = built["census"]
    readout = mf.input_entry(
        s2.READOUT, partition="heldout", role="the screen's readout: hits by its own call"
    )
    census_inputs = run.inputs(built["chroms"])
    entries = [
        *census_inputs,
        run.registration_entry(),
        mf.input_entry(AMENDMENT, partition=None, role="amendment 1: the fold column and the join"),
        readout,
    ]
    base: dict[str, Any] = {
        "lane": s2.LANE,
        "registration": "data/results/secondtargets_registration.json",
        "amendment": f"data/results/{AMENDMENT_NAME}.json",
        "supersedes": (
            f"{SUPERSEDES} (NOT_RUN_NO_READOUT, d4dacb3), written while the readout was not on this machine;"
            " that record stays as written and this file is the registered run on the fetched readout"
        ),
        "question": s2.QUESTION,
        "census": census,
        "shortlist_genome_wide": built["shortlist"],
        "benchmark_rows_by_dataset": built["benchmark_rows_by_dataset"],
        "readout": FETCH,
        "hit": HIT_AMENDED,
        "join": JOIN_AMENDED,
        "limitations": list(s2.LIMITATIONS),
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
        "network": "none during the run; the one fetch is recorded under readout",
    }

    def write(outcome: str, text: str, extra: dict[str, Any] | None = None) -> int:
        out = {
            **base,
            "outcome": outcome,
            "reading": text,
            "falsified": outcome in s2.FALSIFIED_BY,
            "prediction_held": outcome == "MORE_OFTEN",
            **(extra or {}),
            "seconds": round(time.time() - t0, 1),
        }
        print("saved", save_result(RESULT, out, manifest=manifest(entries)))
        print(json.dumps({"outcome": outcome, "reading": text}, indent=1))
        return 0 if outcome != "HARNESS_DEFECT" else 1

    def defect(what: str) -> int:
        return write("HARNESS_DEFECT", s2.reading("HARNESS_DEFECT", what=what))

    pinned = {e["path"]: e["sha256"] for e in reg["result_manifest"]["inputs"]}
    moved = [e["path"] for e in census_inputs if pinned.get(e["path"]) != e["sha256"]]
    if moved:
        return defect(f"pinned inputs moved since the registration: {moved}")
    if readout["sha256"] != FETCH["sha256"]:
        return defect("the readout on disk is not the recorded fetch")
    if census != reg["census"]:
        diff = [k for k in set(census) | set(reg["census"]) if census.get(k) != reg["census"].get(k)]
        return defect(f"the label-free census differs from the registration's on {sorted(diff)}")
    if built["ids"] != reg["population_ids"] or s2.ids_sha256(built["ids"]) != reg["population_sha256"]:
        return defect("the matched shortlist tests' ids or sha256 differ from the registration's")
    if s2.decide(census["matched_shortlist_tests"], True) == "TOO_FEW_PAIRS":
        return write("TOO_FEW_PAIRS", s2.reading("TOO_FEW_PAIRS", **s2.census_numbers(census)))

    # ---- the positive control first, then the matched tests
    try:
        with gzip.open(s2.READOUT, "rt") as fh:
            control_hits = read_amended(fh, wanted_of(t.key for t in built["tss"]))
    except KeyError as e:
        return defect(str(e))
    control = s2.positive_control(built["tss"], control_hits)
    if control["share"] < s2.POSITIVE_CONTROL_MIN:
        what = f"the TSS positive control reads {control}, below the registered {s2.POSITIVE_CONTROL_MIN}"
        return write("HARNESS_DEFECT", s2.reading("HARNESS_DEFECT", what=what), {"positive_control": control})
    units = built["units"]
    keys = {
        k for u in units for k in (u.test.test.key, u.nearest.test.key, *(c.test.key for c in u.controls))
    }
    keys |= {c.test.key for u in units for c in u.compiled}
    with gzip.open(s2.READOUT, "rt") as fh:
        hits = read_amended(fh, wanted_of(keys))
    scored, lost = s2.scored_units(units, hits)
    r = s2.rates(scored)
    extra: dict[str, Any] = {
        "positive_control": control,
        "not_measured": dict(lost),
        "measured_by_arm": measured_by_arm(units, hits),
        "rates": r,
        "amendment_identifier_coverage": am["identifier_coverage"],
        "readout_values_read": len(hits) + len(control_hits),
    }
    if len(scored) < s2.MIN_MATCHED:
        nums = {**s2.census_numbers(census), "matched": f"{len(scored):,} measured"}
        return write("TOO_FEW_PAIRS", s2.reading("TOO_FEW_PAIRS", **nums), extra)
    primary = s2.chromosome_interval(scored, lambda u: u.control_rate)
    nearest = s2.chromosome_interval(scored, lambda u: float(u.nearest_hit))
    outcome = s2.decide(len(scored), True, primary["ci95"])
    text = s2.reading(
        outcome,
        n=f"{r['tests']:,}",
        h=r["hits"],
        r=f"{r['rate_pct']:.2f}",
        rc=f"{r['control_rate_pct']:.2f}",
        d=f"{primary['difference']:+.2f}",
        lo=f"{primary['ci95'][0]:+.2f}",
        hi=f"{primary['ci95'][1]:+.2f}",
        rn=f"{r['nearest_rate_pct']:.2f}",
        dn=f"{nearest['difference']:+.2f}",
        nlo=f"{nearest['ci95'][0]:+.2f}",
        nhi=f"{nearest['ci95'][1]:+.2f}",
    )
    untested = [u for u in scored if u.element_untested]
    compiled = [hits.get(k) for k in sorted({c.test.key for u in units for c in u.compiled})]
    compiled = [c for c in compiled if c is not None]
    by_rank = {b: s2.rates([u for u in scored if s2.rank_band(u.rank) == b]) for b in ("2", "3", "4-5", "6+")}
    extra.update(
        {
            "primary": primary,
            "against_nearest_gene": nearest,
            "descriptives_not_verdicts": {
                "element_independent_subset": {
                    **s2.rates(untested),
                    "difference": s2.chromosome_interval(untested, lambda u: u.control_rate),
                },
                "compiled_targets_at_the_same_elements": {"measured": len(compiled), "hits": sum(compiled)},
                "shortlist_by_rank": by_rank,
            },
        }
    )
    return write(outcome, text, extra)


if __name__ == "__main__":
    if sys.argv[1:] == ["--amend"]:
        raise SystemExit(amend())
    raise SystemExit(main())
