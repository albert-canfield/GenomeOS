"""Reconstruct DHS.RPM and H3K27ac.RPM from the ENCODE alignments, and put the result to the gate
registered in data/results/astroreg_calibration_registration.json before any of this was run.

    uv run --with pysam==0.23.0 --python 3.12 python scripts/astroreg_rpm_calibrate.py \
        --column h3k27ac --cell k562

pysam is a transient tool, MIT, fetched by `--with` and deliberately NOT in pyproject, the way
shellcheck is: this project ships no runtime dependency and that is not changed for one script. It
needs Python 3.12 because pysam 0.23.0 does not import on 3.14 (`libcalignedsegment has no attribute
CMATCH`), and the project itself has no dependencies, so the two coexist in one interpreter.

ONE SEQUENTIAL STREAMING PASS PER BAM. The BAM is read over HTTPS, BGZF decoded as it arrives, and
NOTHING is written to disk. No index is needed because there is no random access, and none is
available: neither astrocyte experiment publishes a `.bai`, the conventional URL is 404, and a portal
search for file_format=bai over both returns nothing. The objects are range-readable, which is why
sequential streaming works at all. Size is therefore a cost in time and never a reason to prefer a
file.

THE REPLICATE CHECK IS THE POINT OF THIS SCRIPT'S FILE SELECTION, not a formality. The registered rule
was originally "pool all released filtered GRCh38 alignments", and that was wrong: ENCSR000AKP
publishes seven such BAMs which are three biological replicates reprocessed by three pipeline
versions, so pooling them would have counted two replicates three times each and reweighted the
replicate mixture without making any single number look absurd. The corrected rule takes the
experiment's RELEASED analysis only -- ENCODE's own designation, one BAM per biological replicate -- and
`select_bams` REFUSES a set in which any biological replicate appears twice. A rule is worth less than
the check that makes breaking it impossible.

THE NUMERATOR AND DENOMINATOR COME FROM THE SAME PASS under the same filter, which `rpm.Counter`
enforces by incrementing both. ENCODE's published mapped-read total is fetched and printed BESIDE the
computed denominator with their ratio, and is never used as the denominator.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.request
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from astroreg_calib_register import (  # noqa: E402
    MEDIAN_RATIO_RANGE,
    SPEARMAN_MIN,
    k562_heldout_elements,
)

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri, rpm  # noqa: E402
from genomeos.attribution.measured import CRISPRI_SPLIT_OF  # noqa: E402
from genomeos.results import save_result  # noqa: E402

PORTAL = "https://www.encodeproject.org"
UA = {"Accept": "application/json", "User-Agent": "genomeos"}

REGISTRATION = Path("data/results/astroreg_calibration_registration.json")

#: The experiments, as the registration names them. Nothing here is chosen by size.
EXPERIMENTS = {
    ("h3k27ac", "k562"): "ENCSR000AKP",
    ("h3k27ac", "astrocyte"): "ENCSR000AOQ",
    ("dnase", "k562"): None,  # selected by measurement; see the registration's selection rule
    ("dnase", "astrocyte"): "ENCSR000EPM",
}

COLUMN_OF = {"h3k27ac": "H3K27ac.RPM", "dnase": "DHS.RPM"}

#: Where a long pass checkpoints itself, so a dropped connection costs one BAM and not the run.
CHECKPOINTS = Path(".git/genomeos-rpm")

OWN_CODE = (
    "genomeos/attribution/rpm.py",
    "scripts/astroreg_calib_register.py",
    "scripts/astroreg_rpm_calibrate.py",
    "tests/test_rpm.py",
)


def portal(path: str) -> dict[str, Any]:
    url = path if path.startswith("http") else f"{PORTAL}{path}"
    sep = "&" if "?" in url else "?"
    with urllib.request.urlopen(
        urllib.request.Request(f"{url}{sep}format=json", headers=UA), timeout=120
    ) as r:
        return json.load(r)


def select_bams(experiment: str) -> dict[str, Any]:
    """The filtered GRCh38 alignment BAMs of the experiment's RELEASED analysis, one per replicate.

    Refuses rather than pools when a biological replicate appears more than once, which is the defect
    the registered rule was amended for.
    """
    exp = portal(f"/experiments/{experiment}/")
    released = [a for a in exp.get("analyses", []) if a.get("status") == "released"]
    if len(released) != 1:
        raise SystemExit(
            f"REFUSED: {experiment} has {len(released)} released analyses, not 1 "
            f"({[a.get('title') for a in released]}). The registered rule names the released "
            "analysis, so a different number of them is a case this rule does not cover and must be "
            "ruled on rather than guessed"
        )

    analysis = released[0]
    wanted = set(analysis.get("files", []))
    bams = []
    for f in exp.get("files", []):
        if f.get("@id") not in wanted:
            continue
        if f.get("status") != "released" or f.get("file_format") != "bam":
            continue
        if f.get("output_type") != "alignments" or f.get("assembly") != "GRCh38":
            continue
        bams.append(
            {
                "accession": f["accession"],
                "biological_replicates": f.get("biological_replicates"),
                "technical_replicates": f.get("technical_replicates"),
                "file_size": f.get("file_size"),
                "mapped_run_type": f.get("mapped_run_type"),
                "mapped_read_length": f.get("mapped_read_length"),
                "href": f.get("href"),
            }
        )
    if not bams:
        raise SystemExit(f"REFUSED: {experiment}'s released analysis holds no filtered GRCh38 BAM")
    try:
        seen = rpm.check_one_bam_per_replicate(bams)
    except ValueError as e:
        raise SystemExit(f"REFUSED for {experiment}: {e}") from e
    return {
        "experiment": experiment,
        "analysis": analysis.get("@id"),
        "analysis_title": analysis.get("title"),
        "bams": sorted(bams, key=lambda b: b["accession"]),
        "biological_replicates": sorted(seen),
        "total_bytes": sum(b["file_size"] or 0 for b in bams),
        "selection_rule": "the released analysis only, one BAM per biological replicate, checked",
    }


def metadata_mapped_reads(accession: str) -> dict[str, Any]:
    """ENCODE's own mapped-read total: a CROSS-CHECK printed beside the computed denominator."""
    f = portal(f"/files/{accession}/")
    out: dict[str, Any] = {"accession": accession, "mapped_reads": None, "from": None}
    for qm in f.get("quality_metrics", []) or []:
        ref = qm if isinstance(qm, dict) else portal(qm)
        for key in ("mapped_reads", "mapped", "reads mapped"):
            if ref.get(key) is not None:
                out["mapped_reads"] = ref[key]
                out["from"] = f"{[t for t in ref.get('@type', []) if 'Quality' in t][:1]}:{key}"
                return out
    return out


def stream_one(url: str, elements: list[tuple[str, int, int]], label: str) -> rpm.Counter:
    """One sequential pass. Nothing is written to disk and no index is used."""
    import pysam

    counter = rpm.Counter(elements)
    t0 = time.time()
    save = pysam.set_verbosity(0)  # the "no index" notice is expected and is not news
    try:
        af = pysam.AlignmentFile(url, "rb")
    finally:
        pysam.set_verbosity(save)
    if af.has_index():
        print(f"  note: {label} reports an index; it is not used, the pass is sequential")
    try:
        for read in af:
            counter.add(read)
            if counter.reads_seen % 5_000_000 == 0:
                el = time.time() - t0
                print(
                    f"    {label}: {counter.reads_seen:,} reads, {counter.denominator:,} kept, "
                    f"{counter.reads_seen / max(el, 1e-9):,.0f}/s, {el / 60:.1f} min",
                    flush=True,
                )
    finally:
        af.close()
    print(
        f"  {label}: {counter.reads_seen:,} reads seen, {counter.denominator:,} kept, "
        f"{sum(counter.counts.values()):,} element overlaps, {(time.time() - t0) / 60:.1f} min",
        flush=True,
    )
    return counter


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--column", choices=sorted(COLUMN_OF), required=True)
    ap.add_argument("--cell", choices=("k562", "astrocyte"), required=True)
    ap.add_argument("--experiment", help="override only for the DNase arm, whose source is selected")
    ap.add_argument("--list-only", action="store_true", help="resolve and check the BAMs, stream none")
    args = ap.parse_args()

    if args.cell == "astrocyte":
        return refuse_astrocyte_until_the_gate_passes()

    experiment = args.experiment or EXPERIMENTS[(args.column, args.cell)]
    if experiment is None:
        raise SystemExit(
            "REFUSED: the K562 DNase source is not fixed by the registration, it is SELECTED by "
            "measurement among the 14 preferred_default candidates before any BAM is streamed. Run "
            "the selection first; do not pass --experiment to skip it"
        )

    chosen = select_bams(experiment)
    print(
        f"{experiment} ({chosen['analysis_title']}): {len(chosen['bams'])} BAM(s), "
        f"{chosen['total_bytes'] / 1e9:.2f} GB, biological replicates {chosen['biological_replicates']}"
    )
    for b in chosen["bams"]:
        print(
            f"   {b['accession']} {(b['file_size'] or 0) / 1e9:5.2f}GB rep={b['biological_replicates']} "
            f"{b['mapped_run_type']} {b['mapped_read_length']}bp"
        )
    if args.list_only:
        print("--list-only: nothing streamed")
        return 0

    comparison = k562_heldout_elements()
    published = _published(comparison, COLUMN_OF[args.column])
    elements = sorted(published)
    print(f"{len(elements)} K562 held-out elements; comparing column {COLUMN_OF[args.column]}")

    CHECKPOINTS.mkdir(parents=True, exist_ok=True)
    counters, cross = [], []
    for b in chosen["bams"]:
        url = (
            f"{PORTAL}{b['href']}"
            if b.get("href")
            else f"{PORTAL}/files/{b['accession']}/@@download/{b['accession']}.bam"
        )
        counters.append(stream_one(url, elements, b["accession"]))
        cross.append(metadata_mapped_reads(b["accession"]))
        ck = CHECKPOINTS / f"{experiment}-{b['accession']}.json"
        ck.write_text(
            json.dumps(
                {
                    "summary": counters[-1].summary(),
                    "counts": {"|".join(map(str, k)): v for k, v in counters[-1].counts.items()},
                }
            )
        )

    pooled = rpm.merge(counters)
    computed = [pooled["rpm"][e] for e in elements]
    pub = [published[e] for e in elements]
    verdict = rpm.gate(computed, pub, SPEARMAN_MIN, MEDIAN_RATIO_RANGE)

    meta_total = sum(c["mapped_reads"] or 0 for c in cross)
    out = {
        "column": COLUMN_OF[args.column],
        "cell": args.cell,
        "files": chosen,
        "pooled_denominator": pooled["denominator"],
        "per_pass_denominator": pooled["per_pass_denominator"],
        "metadata_cross_check": {
            "per_file": cross,
            "metadata_total_mapped_reads": meta_total or None,
            "computed_over_metadata": (pooled["denominator"] / meta_total) if meta_total else None,
            "rule": rpm.METADATA_IS_A_CROSS_CHECK,
        },
        "gate": verdict,
        "elements_with_no_read": sum(1 for v in computed if v == 0),
        "read_rule": rpm.READ_RULE,
    }
    print(json.dumps({k: v for k, v in out.items() if k != "files"}, indent=1, default=str))
    out["status"] = (
        "the registered calibration gate, run once on the K562 held-out elements. No AlphaGenome "
        "request was sent, no money was spent, and no BAM was written to disk"
    )
    out["lane"] = "lane-astro"
    out["registration"] = {
        "path": str(REGISTRATION),
        "spearman_min": SPEARMAN_MIN,
        "median_ratio_range": list(MEDIAN_RATIO_RANGE),
        "both_written_before_this_ran": True,
    }
    out["reading"] = (
        "the gate PASSES for this column: the frozen feature is reconstructible from the alignments "
        "to the registered tolerance"
        if verdict["passes"]
        else "the gate FAILS for this column. Under the read rule registered before this ran, the "
        "published column is NOT reproduced to the registered tolerance. By the registration's own "
        "terms the AstroREG no-go STANDS and no astrocyte alignment is fetched. What is established "
        "is that the column is not recoverable UNDER THIS RULE; since the benchmark documents no "
        "read rule, that is as much a finding about the documentation as about the column, and it is "
        "NOT a licence to try other rules until one passes -- a rule chosen after seeing this number "
        "would be fitted to it"
    )
    out["not_claimed"] = [
        "not a demonstration that the published column is wrong: it is a demonstration that this "
        "registered rule does not reproduce it",
        "not an identification of the rule that would reproduce it, and no such rule was searched for",
        "not a result about astrocytes: no astrocyte alignment was read",
    ]
    out[mf.KEY] = _manifest(chosen, comparison, pooled, verdict)
    path = save_result(f"astroreg_calibration_{args.column}_{args.cell}", out)
    dest = CHECKPOINTS / f"calibration-{args.column}-{args.cell}.json"
    dest.write_text(json.dumps(out, indent=1, default=str))
    print(f"-> {path}")
    print(
        "GATE PASSES"
        if verdict["passes"]
        else "GATE FAILS: the no-go stands and nothing is fetched for astrocytes"
    )
    return 0


def _manifest(chosen, comparison, pooled, verdict) -> dict[str, Any]:
    """What this number was computed from. The alignments are SOURCES: they are streamed, never stored."""
    return {
        "sources": [
            {
                "accession": f"ENCODE {chosen['experiment']} ({chosen['analysis_title']}), filtered "
                f"GRCh38 alignments: {', '.join(b['accession'] for b in chosen['bams'])}",
                "version": f"{chosen['total_bytes'] / 1e9:.2f} GB streamed over HTTPS and not stored; "
                "one BAM per biological replicate from the released analysis",
                "url": PORTAL,
            },
            {
                "accession": "EngreitzLab/CRISPR_comparison, EPCrisprBenchmark heldout_5_cell_types "
                "(Gschwind et al. 2025)",
                "version": "fetched 2026-09-16; pinned here by sha256",
                "url": crispri.BASE_URL,
            },
        ],
        "inputs": [
            mf.input_entry(
                crispri.KNOWLEDGE / crispri.HELDOUT,
                partition=CRISPRI_SPLIT_OF[crispri.HELDOUT],
                role="the K562 held-out elements and the published column this is compared against; "
                "already-read development evidence",
            ),
        ],
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "spearman_min": SPEARMAN_MIN,
            "median_ratio_min": MEDIAN_RATIO_RANGE[0],
            "median_ratio_max": MEDIAN_RATIO_RANGE[1],
            "elements": verdict["elements_compared"],
            "pooled_denominator": pooled["denominator"],
            "bams_streamed": len(chosen["bams"]),
            "bytes_written_to_disk": 0,
            "alphagenome_requests": 0,
            "model_requests": 0,
            "requests_sent": 0,
            "money_spent": 0,
        },
        "exclusions": [
            "no AlphaGenome request is sent and this script has no code path that could send one",
            "no BAM is written to disk: each is streamed once, sequentially, and discarded",
            "no MAPQ threshold and no duplicate rule are applied, because the benchmark states none",
            "no astrocyte alignment is read, and none may be until this gate passes",
            "no alternative read rule was tried: the rule was registered before the run and a second "
            "rule chosen after seeing this number would be fitted to it",
        ],
        "partitions": {
            "k562_heldout": f"{comparison['elements']} elements, already-read development evidence",
        },
        "reproducibility": (
            "the published column and the element set are pinned by sha256; the alignments are a "
            "network source streamed by accession and are not stored, so a rebuild re-streams them "
            "from ENCODE rather than reading bytes from this machine"
        ),
        "code_cleanliness": mf.code_cleanliness("scripts/astroreg_rpm_calibrate.py", OWN_CODE, ROOT),
    }


def _published(comparison: dict[str, Any], column: str) -> dict[tuple[str, int, int], float]:
    """The published column per element. Reads the benchmark again rather than trusting a summary."""
    import gzip

    from genomeos.attribution import crispri

    out: dict[tuple[str, int, int], float] = {}
    with gzip.open(crispri.KNOWLEDGE / crispri.HELDOUT, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            row = line.rstrip("\n").split("\t")
            if row[idx["CellType"]] != "K562":
                continue
            raw = (row[idx[column]] or "").strip()
            if raw in ("", "NA", "nan"):
                continue
            out[(row[idx["chrom"]], int(row[idx["chromStart"]]), int(row[idx["chromEnd"]]))] = float(raw)
    return out


def refuse_astrocyte_until_the_gate_passes() -> int:
    print(
        "REFUSED: no astrocyte alignment is streamed until the K562 calibration gate passes. The "
        "registration says so, and the point of a gate is that it is checked before the thing it "
        "guards, not after"
    )
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
