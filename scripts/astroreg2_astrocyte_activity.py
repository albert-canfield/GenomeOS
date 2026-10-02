"""The astrocyte activity columns, reconstructed under amendment 2's rule, and the ZERO RATE BY LABEL.

    uv run --with pysam==0.23.0 --python 3.12 python scripts/astroreg2_astrocyte_activity.py

This computes an INPUT, not a result. No deletion value is read, no AUPRC is formed, no effect size,
fold change, p-value or FDR is read from the screen, and the endpoint needs both arms so nothing about
it is visible from here. What IS read is the registered LABEL of each pair -- positive, negative,
increase held apart, excluded -- because the question this answers is whether a measured zero is
differential by label, and that cannot be answered without the labels. Those labels were already read
and committed by this lane in data/results/astroreg_activity_nogo.json, which reported the peak-call
censoring by label; they are already-read development evidence and no outcome magnitude accompanies
them.

WHY THE ZERO RATE IS THE QUESTION. The original no-go refused a narrowPeak signalValue because it is
CENSORED: 477 of the screen's 957 elements overlap no astrocyte H3K27ac peak, so a peak call could give
them no value but an invented zero -- and that censoring was DIFFERENTIAL BY LABEL, 15.6% of positives
against 53.4% of negatives, which is why a peak-signal substitute would have produced a gain that was an
artefact of the substitution rather than a measurement. A count from the alignments has no such
asymmetry BY CONSTRUCTION, because every element's span is read whether or not anything was called
there, and an element with no reads gets a MEASURED zero. By construction is not the same as in fact, so
this SHOWS it: the zero rate is reported overall and per label, positives and negatives apart. If it is
differential, that is a finding and it bears on the whole route.

AN ASYMMETRY WITH THE K562 VALIDATION, recorded because any astrocyte number must carry it. For K562 the
producer named its own files, so both the counting rule and the file list had the producer's provenance.
mayasheth/chrom-annotate's resources/metadata/epigenetic_datasets.tsv names only GM12878, HCT116,
Jurkat, K562 and WTC11 -- it names NO astrocyte file. So here the COUNTING RULE is the producer's, cited
to repository, commit, file and line, while the FILE SELECTION is this lane's own stated rule: the
experiment's released analysis, one BAM per biological replicate, enforced by
`rpm.check_one_bam_per_replicate`. The reconstruction was validated against published values in K562
only, and in astrocytes there is no published column to check it against.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from astroreg_register import TABLE3, screen_pairs  # noqa: E402
from astroreg_rpm_calibrate import (  # noqa: E402
    PORTAL,
    metadata_mapped_reads,
    select_bams,
    stream_producer,
)

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import rpm  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "astroreg2_astrocyte_activity"
ENTRY = "scripts/astroreg2_astrocyte_activity.py"
OWN_CODE = (
    "genomeos/attribution/rpm.py",
    "scripts/astroreg2_astrocyte_activity.py",
    "scripts/astroreg_rpm_calibrate.py",
    "tests/test_rpm.py",
)

#: The experiments the ORIGINAL registration names for astrocyte. Not chosen here.
EXPERIMENTS = {"DHS.RPM": "ENCSR000EPM", "H3K27ac.RPM": "ENCSR000AOQ"}

#: The labels the registered test actually scores, kept apart in the report.
SCORED = ("positive", "negative")

CACHE = Path(".git/genomeos-rpm/recon-astrocyte.json")

#: The peak-call censoring the original no-go measured, quoted so the comparison is exact.
PEAK_CENSORING = {
    "elements_with_no_h3k27ac_peak": 477,
    "elements": 957,
    "share": 477 / 957,
    "by_label_in_the_no_go": "15.6% of positives against 53.4% of negatives overlapped no H3K27ac "
    "peak, so a peak-signal substitute would have been differentially missing BY LABEL and the gain it "
    "produced would have been an artefact of the substitution",
    "source": "data/results/astroreg_activity_nogo.json",
}


def zero_rate(
    values: dict[tuple[str, int, int], float], label_of: dict[tuple[str, int, int], set[str]]
) -> dict[str, Any]:
    """The share of elements with a MEASURED zero, overall and per label. The thing to be shown."""
    per: dict[str, dict[str, int]] = defaultdict(lambda: {"elements": 0, "zero": 0})
    for el, v in values.items():
        buckets = {"all labels", *label_of.get(el, set())}
        for b in buckets:
            per[b]["elements"] += 1
            if v == 0.0:
                per[b]["zero"] += 1
    out = {
        k: dict(v) | {"zero_share": (v["zero"] / v["elements"]) if v["elements"] else None}
        for k, v in sorted(per.items())
    }
    pos = out.get("positive", {}).get("zero_share")
    neg = out.get("negative", {}).get("zero_share")
    diff = (pos - neg) if (pos is not None and neg is not None) else None
    return {
        "by_label": out,
        "positive_minus_negative": diff,
        "differential_by_label": (abs(diff) > 0.05) if diff is not None else None,
        "reading": (
            "a measured zero means the element's span was read in the alignments and no read was "
            "there. It is not the censored zero the no-go refused, which recorded the absence of a "
            "measurement as though it were one. The shares per label are reported so that 'no "
            "asymmetry by construction' is shown rather than asserted"
        ),
        "threshold_note": "0.05 is this report's own threshold for calling the rate differential and "
        "is not a registered term; the shares themselves are the evidence",
    }


def one_column(column: str, elements: list[tuple[str, int, int]], list_only: bool) -> dict[str, Any]:
    """Stream one experiment's BAMs and average the per-BAM RPMs, as the producer's rule does."""
    experiment = EXPERIMENTS[column]
    chosen = select_bams(experiment)
    print(
        f"\n--- {column}: {experiment} ({chosen['analysis_title']}), {len(chosen['bams'])} BAM(s), "
        f"{chosen['total_bytes'] / 1e9:.2f} GB, biological replicates {chosen['biological_replicates']}"
    )
    for b in chosen["bams"]:
        print(
            f"   {b['accession']} {(b['file_size'] or 0) / 1e9:5.2f}GB rep={b['biological_replicates']} "
            f"{b['mapped_run_type']} {b['mapped_read_length']}bp"
        )
    if list_only:
        return {"chosen": chosen, "listed_only": True}
    per_bam, summaries, cross = [], [], []
    for b in chosen["bams"]:
        url = (
            f"{PORTAL}{b['href']}"
            if b.get("href")
            else f"{PORTAL}/files/{b['accession']}/@@download/{b['accession']}.bam"
        )
        counter = stream_producer(url, elements, b["accession"])
        per_bam.append(counter.rpm())
        summaries.append(counter.summary() | {"accession": b["accession"]})
        cross.append(metadata_mapped_reads(b["accession"]))
    averaged = rpm.mean_of_per_bam_rpm(per_bam)
    computed = sum(s["denominator_reads_passing_the_filter"] for s in summaries)
    meta = sum(c["mapped_reads"] or 0 for c in cross)
    return {
        "chosen": chosen,
        "per_bam": summaries,
        "rpm": averaged,
        "combination": "mean of the per-BAM RPMs, as average_features does",
        "metadata_cross_check": {
            "per_file": cross,
            "computed_denominator_total": computed,
            "metadata_total_mapped_reads": meta or None,
            "computed_over_metadata": (computed / meta) if meta else None,
            "rule": rpm.METADATA_IS_A_CROSS_CHECK,
        },
    }


def manifest(columns: dict[str, Any], pairs: list[dict], elements: list) -> dict[str, Any]:
    names = [b["accession"] for c in columns.values() for b in c["chosen"]["bams"]]
    gb = sum(c["chosen"]["total_bytes"] for c in columns.values()) / 1e9
    return {
        "sources": [
            {
                "accession": "ENCODE astrocyte alignments, released analysis, one per biological "
                f"replicate: {', '.join(names)}",
                "version": f"{gb:.2f} GB streamed over HTTPS and not stored",
                "url": PORTAL,
            },
            {
                "accession": "mayasheth/chrom-annotate, workflow/scripts/neighborhoods.py",
                "version": "commit 91cda73ebe3a19153a582cab18cbf7ff70d85cfc; the COUNTING RULE is read "
                "from this file line by line. It names NO astrocyte file, so the file selection is "
                "this lane's stated rule and not the producer's",
                "url": "https://github.com/mayasheth/chrom-annotate",
            },
            {
                "accession": "Green NFO, et al. CRISPRi screening in cultured human astrocytes. "
                "Nature Neuroscience 2025;29(3):703-716 -- Supplementary Table 3",
                "version": "publisher open-access supplementary store; pinned here by sha256",
            },
        ],
        "inputs": [
            mf.input_entry(
                TABLE3,
                partition=None,
                role="the screen's element coordinates and registered labels. No effect size, fold "
                "change, p-value, FDR or expression level is read",
            ),
        ],
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "elements": len(elements),
            "screen_pairs": len(pairs),
            "bams_streamed": len(names),
            "gb_streamed": round(gb, 2),
            "bytes_written_to_disk": 0,
            "alphagenome_requests": 0,
            "model_requests": 0,
            "requests_sent": 0,
            "money_spent": 0,
            "deletion_values_read": 0,
            "auprc_computed": 0,
        },
        "exclusions": [
            "no AlphaGenome request is sent and this script has no code path that could send one",
            "no deletion value is read or requested, so neither arm of the endpoint exists and nothing "
            "about the result is visible from this work",
            "no effect size, fold change, p-value, FDR or expression level is read from the screen; the "
            "registered LABEL is read, and only to report the zero rate by label",
            "no BAM is written to disk",
            "no rescaling factor is fitted: the columns are as counted",
        ],
        "partitions": {
            "astrocyte_screen": f"all {len(pairs)} pairs over {len(elements)} elements, every label",
        },
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--list-only", action="store_true", help="resolve and check the BAMs, stream none")
    args = ap.parse_args()
    t0 = time.time()

    pairs = screen_pairs()
    label_of: dict[tuple[str, int, int], set[str]] = defaultdict(set)
    for p in pairs:
        label_of[(p["chrom"], p["start"], p["end"])].add(p["label"])
    elements = sorted(label_of)
    print(f"{len(pairs)} screen pairs over {len(elements)} elements")

    columns = {c: one_column(c, elements, args.list_only) for c in ("DHS.RPM", "H3K27ac.RPM")}
    if args.list_only:
        print("--list-only: nothing streamed")
        return 0

    zeros = {c: zero_rate(columns[c]["rpm"], label_of) for c in columns}
    both_zero = sum(
        1 for e in elements if columns["DHS.RPM"]["rpm"][e] == 0.0 and columns["H3K27ac.RPM"]["rpm"][e] == 0.0
    )

    payload = {
        "status": "an INPUT computed, not a result. No deletion value read, no AUPRC formed, no "
        "outcome magnitude read, no AlphaGenome request, no money, nothing written to disk",
        "lane": "lane-astro",
        "what_this_is": "the astrocyte activity columns under amendment 2's rule, and the zero rate by "
        "label, which is the thing the no-go's censoring argument turns on",
        "rule": {
            "counting": "the producer's, cited term by term",
            "producer": rpm.PRODUCER,
            "file_selection": "NOT the producer's: it names no astrocyte file. The released analysis "
            "only, one BAM per biological replicate, enforced by rpm.check_one_bam_per_replicate",
            "asymmetry_with_the_k562_validation": "for K562 both the rule and the files were the "
            "producer's. Here the rule is the producer's and the files are this lane's stated "
            "selection, and there is no published astrocyte column to validate the result against. Any "
            "astrocyte number must carry that: activity reconstructed by the K562-validated rule",
            "no_rescaling": "the columns are as counted; no factor was fitted",
        },
        "columns": {
            c: {
                "experiment": columns[c]["chosen"]["experiment"],
                "analysis": columns[c]["chosen"]["analysis_title"],
                "bams": [b["accession"] for b in columns[c]["chosen"]["bams"]],
                "biological_replicates": columns[c]["chosen"]["biological_replicates"],
                "per_bam": columns[c]["per_bam"],
                "combination": columns[c]["combination"],
                "metadata_cross_check": columns[c]["metadata_cross_check"],
                "zero_rate": zeros[c],
            }
            for c in columns
        },
        "elements_with_zero_in_both_columns": both_zero,
        "the_peak_call_this_replaces": PEAK_CENSORING,
        "a_measured_zero": rpm.A_MEASURED_ZERO,
        "not_claimed": [
            "not validated against a published astrocyte column: none exists. The rule was validated "
            "on K562 and is applied here unchanged",
            "not a result about the deletion gain, and not an input to one until a deletion value is "
            "bought under an approval that names AstroREG-2",
            "not a claim that the activity values are correct in magnitude: the K562 scale gate failed "
            "and stays failed, and what was validated is that the ENDPOINT is insensitive to that",
            "no astrocyte outcome magnitude was read",
        ],
        "requests_sent": 0,
        "money_spent": 0,
    }
    payload[mf.KEY] = manifest(columns, pairs, elements)
    path = save_result(RESULT, payload)

    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(
        json.dumps({c: {"|".join(map(str, k)): v for k, v in columns[c]["rpm"].items()} for c in columns})
    )

    print()
    for c in columns:
        z = zeros[c]["by_label"]
        print(
            f"{c}: zero in {z['all labels']['zero']}/{z['all labels']['elements']} elements "
            f"({z['all labels']['zero_share']:.2%})"
        )
        for lab in SCORED:
            if lab in z:
                print(f"    {lab:24s} {z[lab]['zero']}/{z[lab]['elements']} = {z[lab]['zero_share']:.2%}")
        print(f"    positives minus negatives: {zeros[c]['positive_minus_negative']}")
        print(f"    differential by label: {zeros[c]['differential_by_label']}")
        x = columns[c]["metadata_cross_check"]
        print(
            f"    denominator computed {x['computed_denominator_total']:,} vs ENCODE metadata "
            f"{x['metadata_total_mapped_reads']:,} -> ratio {x['computed_over_metadata']}"
        )
    print(f"zero in BOTH columns: {both_zero}/{len(elements)}")
    print(f"-> {path}")
    print(f"({time.time() - t0:.0f} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
