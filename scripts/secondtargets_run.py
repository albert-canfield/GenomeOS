# SPDX-License-Identifier: AGPL-3.0-or-later
"""The registered run: the model's second targets against Gasperini 2019's tests independent of the
CRISPRi benchmark (genomeos/attribution/secondtargets.py).

    uv run python scripts/secondtargets_run.py

Runs only from the committed registration (data/results/secondtargets_registration.json) and refuses if
any registered field differs from the module, or if its own code or the code it calls differs from what
the registration pinned. Order, as registered (secondtargets.DECISION and STOP_RULE):
  1. every pinned input is hashed and compared with the registration's;
  2. the label-free census is recomputed (pair table, independence, shortlist, matching) and must equal
     the registered census and the registered ids of the matched shortlist tests;
  3. TOO_FEW_PAIRS or NOT_RUN_NO_READOUT is decided here, before the readout is opened;
  4. only if the readout table is on this machine: its header is checked, then the TSS positive control
     is read, then the matched tests' readout, and the registered difference is decided.

No model request, no network, no money: every input is on this machine or the run stops at step 3.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import sys
import time
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri  # noqa: E402
from genomeos.attribution import measured_unknown as mu  # noqa: E402
from genomeos.attribution import secondtargets as s2  # noqa: E402
from genomeos.attribution import targetchoice as tc  # noqa: E402
from genomeos.attribution.measured import CRISPRI_SPLIT_OF  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

RESULT = "secondtargets_result"
REGISTRATION = RESULTS_DIR / "secondtargets_registration.json"
OWN_CODE = (
    "genomeos/attribution/secondtargets.py",
    "scripts/secondtargets_register.py",
    "scripts/secondtargets_run.py",
    "tests/test_secondtargets.py",
)
#: the code the run must equal, pinned by sha256 in the registration; the tests are not on the run's path
FROZEN = OWN_CODE[:3]
#: the committed code the census runs, not this lane's; pinned by sha256 in the registration
STAGE_CODE = (
    "genomeos/attribution/targetchoice.py",
    "genomeos/attribution/signaltrace.py",
    "genomeos/attribution/compile.py",
    "genomeos/attribution/crispri.py",
    "genomeos/attribution/measured_unknown.py",
    "genomeos/attribution/targets.py",
    "genomeos/predict/enhancer_target.py",
    "genomeos/genome/annotation.py",
    "genomeos/ir/model.py",
    "genomeos/lang/parser.py",
)
CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]


def code_sha256(paths: tuple[str, ...]) -> dict[str, str]:
    return {p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths}


def program_path(chrom: str) -> Path:
    return Path(COMPILED_DIR) / f"noncoding_{chrom}.bio"


def gencode_path(chrom: str) -> Path:
    """The GENCODE v50 file the project reads for one chromosome (annotation.default_gencode)."""
    from genomeos.genome.annotation import default_gencode

    path = default_gencode({chrom})
    if path is None:
        raise SystemExit(f"refused: no GENCODE annotation for {chrom} on this machine")
    return Path(path)


def coding_ids(chroms: list[str]) -> dict[str, dict[str, str]]:
    """Chromosome -> {Ensembl id without version: GENCODE v50 symbol} of its protein-coding genes."""
    from genomeos.genome import Annotation

    out: dict[str, dict[str, str]] = {}
    for c in chroms:
        ann = Annotation.from_gff3(gencode_path(c), {c})
        out[c] = {s2.ensg_of(g.id): g.symbol for g in ann.protein_coding() if g.locus.chrom == c}
    return out


def lookups_for(
    chroms: list[str], coding: dict[str, dict[str, str]], log: Any = print
) -> dict[str, tc.ChoiceLookup]:
    """Every chromosome's compiled elements with their candidate targets, as targetchoice builds them."""
    from genomeos.attribution.targets import ElementResponses

    t0 = time.time()
    responses = ElementResponses(crispri.ELEMENT_CACHE)
    out: dict[str, tc.ChoiceLookup] = {}
    for chrom in chroms:
        symbols = set(coding.get(chrom, {}).values())
        out[chrom] = tc.chromosome_choice(chrom, program_path(chrom).read_text(), symbols, responses)
        log(chrom, dict(out[chrom].counts), f"{time.time() - t0:.0f} s", flush=True)
    return out


def build(log: Any = print) -> dict[str, Any]:
    """The label-free census: the screen's tests, independence, the shortlist, the matching."""
    chain = mu.Chain.from_file(s2.CHAIN)
    with gzip.open(s2.PAIR_TABLE, "rt") as fh:
        tests, tss, pair_counts = s2.read_pair_table(fh, chain.lift)
    log("pair table", dict(pair_counts), flush=True)
    index = s2.BenchmarkIndex.from_files()
    chroms = [c for c in CHROMS if c in {t.chrom for t in tests}]
    coding = coding_ids(chroms)
    lookups = lookups_for(chroms, coding, log)
    rows, split = s2.rows_of(tests, index, coding, lookups)
    units, unmatched = s2.matched(rows)
    ids = s2.population_ids(units)
    return {
        "chroms": chroms,
        "tss": tss,
        "units": units,
        "ids": ids,
        "census": s2.census(pair_counts, rows, split, units, unmatched),
        "benchmark_rows_by_dataset": dict(index.datasets.most_common()),
        "shortlist": s2.shortlist_census(lookups),
    }


def inputs(chroms: list[str]) -> list[dict[str, Any]]:
    """Every file the census reads, pinned by sha256 (the readout is declared apart, if held)."""
    return [
        mf.input_entry(s2.PAIR_TABLE, partition=None, role="the screen's tests: coordinates and genes only"),
        mf.input_entry(s2.CHAIN, partition=None, role="UCSC hg19 to hg38 chain, for the screen's elements"),
        mf.input_entry(
            crispri.KNOWLEDGE / crispri.TRAINING,
            partition=CRISPRI_SPLIT_OF[crispri.TRAINING],
            role="the benchmark, read by element, gene and dataset only: the independence exclusion",
        ),
        mf.input_entry(
            crispri.KNOWLEDGE / crispri.HELDOUT,
            partition=CRISPRI_SPLIT_OF[crispri.HELDOUT],
            role="the benchmark, read by element, gene and dataset only: the independence exclusion",
        ),
        mf.input_entry(
            crispri.ELEMENT_CACHE,
            partition=None,
            role="ElementResponses: every gene's predicted deletion response at an element (the shortlist)",
        ),
        mf.files_entry(
            "the compiled non-coding programs (targetchoice.chromosome_choice)",
            [str(program_path(c)) for c in chroms],
            partition=None,
        ),
        mf.files_entry(
            "GENCODE v50 per-chromosome GFF3 (protein-coding genes by Ensembl id)",
            [str(gencode_path(c)) for c in chroms],
            partition=None,
        ),
    ]


def registration_entry() -> dict[str, Any]:
    return mf.input_entry(
        REGISTRATION, partition=None, role="this run's registration, read before any reading"
    )


def readout_entry() -> dict[str, Any]:
    return mf.input_entry(s2.READOUT, partition="heldout", role="the screen's readout: hits by its own call")


def manifest(entries: list[dict[str, Any]], own_entry: str = __file__) -> dict[str, Any]:
    return {
        "sources": [
            {
                "accession": "GEO GSE120861 (Gasperini et al. 2019, Cell), at-scale screen pair table",
                "version": "fetched 2026-10-10 by lane-measuredunknown; pinned here by sha256",
                "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120861/suppl/",
            },
            {"accession": "UCSC hg19ToHg38.over.chain", "version": "pinned by sha256"},
            {
                "accession": "EngreitzLab/CRISPR_comparison EPCrisprBenchmark training_K562 and"
                " heldout_5_cell_types (Gschwind et al.)",
                "version": "fetched 2026-09-16; pinned here by sha256",
                "url": crispri.BASE_URL,
            },
            {
                "accession": "ENCODE SCREEN cCREs scored by AlphaGenome deletion"
                " (per-element response cache)",
                "version": "the 2026-09 all-element sweep; pinned here by sha256",
            },
            {"accession": "the compiled non-coding programs (data/knowledge/compiled)", "version": "pinned"},
            {"accession": "GENCODE v50 annotation (protein-coding genes)", "version": "pinned by sha256"},
        ],
        "inputs": entries,
        "input_count": len(entries),
        "assembly": "GRCh38 (the screen's hg19 elements lifted); distances in hg19 from the screen's table",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "min_effect": s2.MIN_EFFECT,
            "guide_group": s2.GUIDE_GROUP,
            "distance_factor": s2.DISTANCE_FACTOR,
            "min_matched": s2.MIN_MATCHED,
            "hit": s2.HIT,
            "positive_control_min": s2.POSITIVE_CONTROL_MIN,
            "bootstrap": s2.BOOTSTRAP,
        },
        "exclusions": [
            s2.POPULATION,
            s2.INDEPENDENCE,
            "experimental rule lines of the compiled programs (signaltrace.PREDICTED_ONLY)",
        ],
        "partitions": {
            "training": "none of this test: the benchmark's rows are read only to exclude the tests"
            " they hold",
            "heldout": "the screen's tests independent of the benchmark; their readout is read once, if held",
        },
        "code_cleanliness": mf.code_cleanliness(own_entry, OWN_CODE),
    }


def check_registration() -> dict[str, Any]:
    reg = json.loads(REGISTRATION.read_text())
    bad = [f for f, v in s2.registered_fields().items() if reg.get(f) != v]
    if bad:
        raise SystemExit(f"refused: the registration and secondtargets.py disagree on {', '.join(bad)}")
    moved = [p for p, h in code_sha256(FROZEN).items() if reg["frozen_code"]["sha256"].get(p) != h]
    moved += [p for p, h in code_sha256(STAGE_CODE).items() if reg["stage_code"]["sha256"].get(p) != h]
    if moved:
        raise SystemExit(f"refused: code differs from what the registration pinned: {moved}")
    return reg


def _pct(x: float) -> str:
    return f"{x:.2f}"


def _signed(x: float) -> str:
    return f"{x:+.2f}"


def main() -> int:
    t0 = time.time()
    if (RESULTS_DIR / f"{RESULT}.json").exists():
        raise SystemExit("refused: the registered run has already been written")
    reg = check_registration()
    built = build()
    entries = [*inputs(built["chroms"]), registration_entry()]
    census = built["census"]
    base: dict[str, Any] = {
        "lane": s2.LANE,
        "registration": "data/results/secondtargets_registration.json",
        "registration_date": reg.get("date"),
        "question": s2.QUESTION,
        "census": census,
        "shortlist_genome_wide": built["shortlist"],
        "benchmark_rows_by_dataset": built["benchmark_rows_by_dataset"],
        "readout": {"path": str(s2.READOUT), "url": s2.READOUT_URL, "held": s2.READOUT.exists()},
        "limitations": list(s2.LIMITATIONS),
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
        "network": "none",
    }

    def write(outcome: str, text: str, extra: dict[str, Any] | None = None, more: list | None = None) -> int:
        out = {
            **base,
            "outcome": outcome,
            "reading": text,
            "falsified": outcome in s2.FALSIFIED_BY,
            "prediction_held": outcome == "MORE_OFTEN",
            **(extra or {}),
            "seconds": round(time.time() - t0, 1),
        }
        print("saved", save_result(RESULT, out, manifest=manifest([*entries, *(more or [])])))
        print(json.dumps({"outcome": outcome, "reading": text}, indent=1))
        return 0 if outcome != "HARNESS_DEFECT" else 1

    def defect(what: str) -> int:
        return write("HARNESS_DEFECT", s2.reading("HARNESS_DEFECT", what=what))

    pinned = {e["path"]: e["sha256"] for e in reg["result_manifest"]["inputs"]}
    moved = [e["path"] for e in entries[:-1] if pinned.get(e["path"]) != e["sha256"]]
    if moved:
        return defect(f"pinned inputs moved since the registration: {moved}")
    if census != reg["census"]:
        diff = [k for k in set(census) | set(reg["census"]) if census.get(k) != reg["census"].get(k)]
        return defect(f"the label-free census differs from the registration's on {sorted(diff)}")
    if built["ids"] != reg["population_ids"] or s2.ids_sha256(built["ids"]) != reg["population_sha256"]:
        return defect("the matched shortlist tests' ids or sha256 differ from the registration's")

    outcome = s2.decide(census["matched_shortlist_tests"], s2.READOUT.exists())
    if outcome in ("TOO_FEW_PAIRS", "NOT_RUN_NO_READOUT"):
        return write(outcome, s2.reading(outcome, **s2.census_numbers(census)), {"readout_values_read": 0})

    # ---- the readout is on this machine: header, positive control, then the matched tests -----------------
    more = [readout_entry()]
    try:
        with gzip.open(s2.READOUT, "rt") as fh:
            control_hits = s2.read_readout(fh, {t.key for t in built["tss"]})
    except KeyError as e:
        return defect(f"{e}; a column-name amendment may be registered before any value is read")
    control = s2.positive_control(built["tss"], control_hits)
    if control["share"] < s2.POSITIVE_CONTROL_MIN:
        return write(
            "HARNESS_DEFECT",
            s2.reading(
                "HARNESS_DEFECT",
                what=f"the TSS positive control reads {control}, below the registered"
                f" {s2.POSITIVE_CONTROL_MIN}",
            ),
            {"positive_control": control},
            more,
        )
    units = built["units"]
    wanted = {
        k for u in units for k in (u.test.test.key, u.nearest.test.key, *(c.test.key for c in u.controls))
    }
    wanted |= {c.test.key for u in units for c in u.compiled}
    with gzip.open(s2.READOUT, "rt") as fh:
        hits = s2.read_readout(fh, wanted)
    scored, lost = s2.scored_units(units, hits)
    primary = s2.chromosome_interval(scored, lambda u: u.control_rate)
    nearest = s2.chromosome_interval(scored, lambda u: float(u.nearest_hit))
    outcome = (
        s2.decide(len(scored), True, primary["ci95"]) if len(scored) >= s2.MIN_MATCHED else "TOO_FEW_PAIRS"
    )
    r = s2.rates(scored)
    if outcome == "TOO_FEW_PAIRS":
        text = s2.reading(outcome, **{**s2.census_numbers(census), "matched": f"{len(scored):,} measured"})
    else:
        text = s2.reading(
            outcome,
            n=f"{r['tests']:,}",
            h=r["hits"],
            r=_pct(r["rate_pct"]),
            rc=_pct(r["control_rate_pct"]),
            d=_signed(primary["difference"]),
            lo=_signed(primary["ci95"][0]),
            hi=_signed(primary["ci95"][1]),
            rn=_pct(r["nearest_rate_pct"]),
            dn=_signed(nearest["difference"]),
            nlo=_signed(nearest["ci95"][0]),
            nhi=_signed(nearest["ci95"][1]),
        )
    untested = [u for u in scored if u.element_untested]
    compiled = [hits.get(k) for k in sorted({c.test.key for u in units for c in u.compiled})]
    compiled = [c for c in compiled if c is not None]
    by_rank: dict[str, Any] = {}
    for band in ("2", "3", "4-5", "6+"):
        sub = [u for u in scored if s2.rank_band(u.rank) == band]
        by_rank[band] = s2.rates(sub)
    extra = {
        "positive_control": control,
        "not_measured": dict(lost),
        "rates": r,
        "primary": primary,
        "against_nearest_gene": nearest,
        "descriptives_not_verdicts": {
            "element_independent_subset": {
                **s2.rates(untested),
                "difference": s2.chromosome_interval(untested, lambda u: u.control_rate),
            },
            "compiled_targets_at_the_same_elements": {
                "measured": len(compiled),
                "hits": sum(compiled),
            },
            "shortlist_by_rank": by_rank,
        },
        "readout_values_read": len(hits) + len(control_hits),
    }
    return write(outcome, text, extra, more)


if __name__ == "__main__":
    raise SystemExit(main())
