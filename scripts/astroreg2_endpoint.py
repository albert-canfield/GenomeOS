"""AstroREG-2's gate: does the ENDPOINT move when the activity columns are reconstructed?

    uv run --with pysam==0.23.0 --python 3.12 python scripts/astroreg2_endpoint.py

Registered in data/results/astroreg2_registration.json before either gain was computed: PASS only if
BOTH |gain_reconstructed - gain_published| <= 0.01 AND the Spearman of the per-pair scores between the
two runs >= 0.99.

WHAT IS HELD FIXED, because that is what makes the comparison about the input rather than the fit. The
weights are fitted ONCE, on the K562 TRAINING pairs, with the PUBLISHED activity columns, and reused
unchanged for both runs. The deletion features are unchanged too: re-annotating a pair after swapping its
activity recomputes the activity terms from the new values and reaches the same deletion value through
the same cache, so the only thing that differs between the two runs is the activity term on the HELD-OUT
pairs.

NO RESCALING. The reconstructed columns enter exactly as computed. Fitting any factor, even one derived
from K562, would make the astrocyte values depend on a K562 quantity, and that is a refitting of the
feature definition which the claim forbids.

NOTHING IS SENT AND NO ASTROCYTE VALUE IS COMPUTED. The alignments are streamed once and discarded, as
amendment 2's rule requires; the per-element reconstruction is cached under .git so a re-run of the
endpoint arithmetic does not re-stream 15.85 GB.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from astroreg2_register import (  # noqa: E402
    MAX_ABS_GAIN_SHIFT,
    MIN_SCORE_SPEARMAN,
    ORIGINAL_SHA256,
)
from astroreg_rpm_calibrate import (  # noqa: E402
    PORTAL,
    PRODUCER_BAMS,
    producer_bam_meta,
    stream_producer,
)

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri, rpm  # noqa: E402
from genomeos.attribution.measured import CRISPRI_SPLIT_OF  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "astroreg2_endpoint_k562"
ENTRY = "scripts/astroreg2_endpoint.py"
OWN_CODE = (
    "genomeos/attribution/rpm.py",
    "scripts/astroreg2_endpoint.py",
    "scripts/astroreg2_register.py",
    "scripts/astroreg_rpm_calibrate.py",
    "tests/test_rpm.py",
)

CACHE = Path(".git/genomeos-rpm/recon-k562-heldout.json")

REGISTRATION = Path("data/results/astroreg2_registration.json")


def heldout_k562() -> list[Any]:
    with gzip.open(crispri.KNOWLEDGE / crispri.HELDOUT, "rt") as fh:
        return [p for p in crispri.parse(fh) if p.cell == "K562"]


def training_k562() -> list[Any]:
    with gzip.open(crispri.KNOWLEDGE / crispri.TRAINING, "rt") as fh:
        return list(crispri.parse(fh))


def crispri_column(col: str) -> str:
    return {"dnase": "DHS.RPM", "h3k27ac": "H3K27ac.RPM"}[col]


def reconstruct(elements: list[tuple[str, int, int]], force: bool) -> dict[str, Any]:
    """Both reconstructed columns per element, under amendment 2's rule. Cached, never on the registry."""
    if CACHE.exists() and not force:
        got = json.loads(CACHE.read_text())
        print(
            f"reusing the cached reconstruction from {CACHE} "
            f"({got['streamed_bytes'] / 1e9:.2f} GB was streamed when it was built; nothing re-streamed)"
        )
        return got
    out: dict[str, Any] = {"columns": {}, "streamed_bytes": 0, "bams": []}
    for col in ("dnase", "h3k27ac"):
        experiment, accessions = PRODUCER_BAMS[(col, "k562")]
        metas = [producer_bam_meta(a) for a in accessions]
        print(f"\n--- {crispri_column(col)}: {experiment}, {len(metas)} BAM(s) as the producer names them")
        per_bam = []
        for m in metas:
            url = (
                f"{PORTAL}{m['href']}"
                if m.get("href")
                else f"{PORTAL}/files/{m['accession']}/@@download/{m['accession']}.bam"
            )
            counter = stream_producer(url, elements, m["accession"])
            per_bam.append(counter.rpm())
            out["streamed_bytes"] += m["file_size"] or 0
            out["bams"].append(
                {
                    "column": crispri_column(col),
                    "accession": m["accession"],
                    "denominator": counter.denominator,
                }
            )
        averaged = rpm.mean_of_per_bam_rpm(per_bam)
        out["columns"][crispri_column(col)] = {"|".join(map(str, k)): v for k, v in averaged.items()}
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(out))
    return out


def gain_and_scores(rows: list[Any], weights: dict[str, list[float]]) -> dict[str, Any]:
    """The project's own headline estimator: average precision of each arm, and their difference."""
    with_d = crispri.FEATURES["activity + distance + deletion"]
    without = crispri.FEATURES["activity + distance"]
    lab = [p.regulated for p in rows]
    s_with = crispri.logistic_score(weights["activity + distance + deletion"], crispri.matrix(rows, with_d))
    s_without = crispri.logistic_score(weights["activity + distance"], crispri.matrix(rows, without))
    ap_with = crispri.average_precision(s_with, lab) or 0.0
    ap_without = crispri.average_precision(s_without, lab) or 0.0
    return {
        "auprc_with_deletion": ap_with,
        "auprc_without_deletion": ap_without,
        "gain": ap_with - ap_without,
        "s_with": s_with,
        "s_without": s_without,
    }


def manifest(recon, heldout, training, verdict) -> dict[str, Any]:
    names = [b["accession"] for b in recon["bams"]]
    return {
        "sources": [
            {
                "accession": "ENCODE alignments named by mayasheth/chrom-annotate "
                f"resources/metadata/epigenetic_datasets.tsv: {', '.join(names)}",
                "version": f"{recon['streamed_bytes'] / 1e9:.2f} GB streamed over HTTPS and not stored",
                "url": PORTAL,
            },
            {
                "accession": "mayasheth/chrom-annotate, workflow/scripts/neighborhoods.py",
                "version": "commit 91cda73ebe3a19153a582cab18cbf7ff70d85cfc; the reconstruction rule "
                "is read from this file, line by line",
                "url": "https://github.com/mayasheth/chrom-annotate",
            },
            {
                "accession": "EngreitzLab/CRISPR_comparison, EPCrisprBenchmark training_K562 and "
                "heldout_5_cell_types",
                "version": "fetched 2026-09-16; pinned here by sha256",
                "url": crispri.BASE_URL,
            },
        ],
        "inputs": [
            mf.input_entry(
                Path("data/results/astroreg_registration.json"),
                partition=None,
                role="the frozen original, whose sha256 is checked; no term of it is moved and its "
                "failed calibration gate is not reopened",
            ),
            mf.input_entry(
                REGISTRATION,
                partition=None,
                role="AstroREG-2's registered gate, read rather than restated",
            ),
            mf.input_entry(
                crispri.KNOWLEDGE / crispri.TRAINING,
                partition=CRISPRI_SPLIT_OF[crispri.TRAINING],
                role="the K562 training pairs the frozen weights are fitted on, with the PUBLISHED "
                "columns, once, for both runs; already-read development evidence",
            ),
            mf.input_entry(
                crispri.KNOWLEDGE / crispri.HELDOUT,
                partition=CRISPRI_SPLIT_OF[crispri.HELDOUT],
                role="the K562 held-out pairs, their labels and their published activity columns: the "
                "endpoint comparison set; already-read",
            ),
        ],
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "max_abs_gain_shift": MAX_ABS_GAIN_SHIFT,
            "min_score_spearman": MIN_SCORE_SPEARMAN,
            "heldout_pairs": len(heldout),
            "training_pairs": len(training),
            "gain_shift": verdict["gain_shift"],
            "bams_streamed": len(names),
            "bytes_written_to_disk": 0,
            "alphagenome_requests": 0,
            "model_requests": 0,
            "requests_sent": 0,
            "money_spent": 0,
            "astrocyte_values_computed": 0,
        },
        "exclusions": [
            "no AlphaGenome request is sent and this script has no code path that could send one",
            "no astrocyte value of any kind is read or computed",
            "no rescaling factor was fitted: the reconstructed columns entered as computed",
            "the weights were not refitted for the reconstructed run",
            "the failed scale gate was not re-run and is not reopened",
            "no BAM is written to disk",
        ],
        "partitions": {
            "k562_training": "the frozen weights' source, published columns only",
            "k562_heldout": "the endpoint comparison set, scored twice",
            "astrocyte": "not touched",
        },
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--restream", action="store_true", help="ignore the cached reconstruction")
    args = ap.parse_args()
    t0 = time.time()

    digest, _, _ = mf.sha256_of(Path("data/results/astroreg_registration.json"))
    if digest != ORIGINAL_SHA256:
        print(f"refused: the original registration is sha256 {digest}, not {ORIGINAL_SHA256}")
        return 2
    if not REGISTRATION.exists():
        print(f"refused: {REGISTRATION} is absent. The gate is registered before it is run")
        return 2

    heldout, training = heldout_k562(), training_k562()
    elements = sorted({(p.chrom, p.start, p.end) for p in heldout})
    print(
        f"{len(heldout)} K562 held-out pairs over {len(elements)} elements; "
        f"{len(training)} K562 training pairs"
    )

    table, cache = crispri.DeletionTable(), crispri.ElementCache()
    print("annotating with the PUBLISHED columns and fitting the frozen weights on training")
    crispri.annotate(training, table, cache)
    crispri.annotate(heldout, table, cache)
    labels = [p.regulated for p in training]
    weights = {
        name: crispri.logistic_fit(crispri.matrix(training, cols), labels)
        for name, cols in crispri.FEATURES.items()
    }
    published = gain_and_scores(heldout, weights)
    print(
        f"  published gain {published['gain']:.6f} (with {published['auprc_with_deletion']:.6f}, "
        f"without {published['auprc_without_deletion']:.6f})"
    )

    recon = reconstruct(elements, args.restream)
    dhs = recon["columns"]["DHS.RPM"]
    k27 = recon["columns"]["H3K27ac.RPM"]
    missing = [e for e in elements if "|".join(map(str, e)) not in dhs]
    if missing:
        print(f"refused: {len(missing)} elements have no reconstructed value")
        return 2

    print("swapping in the RECONSTRUCTED columns, no rescaling, and re-annotating the held-out pairs")
    for p in heldout:
        key = "|".join(map(str, (p.chrom, p.start, p.end)))
        p.dhs, p.h3k27ac = dhs[key], k27[key]
    crispri.annotate(heldout, table, cache)
    reconstructed = gain_and_scores(heldout, weights)
    print(
        f"  reconstructed gain {reconstructed['gain']:.6f} "
        f"(with {reconstructed['auprc_with_deletion']:.6f}, "
        f"without {reconstructed['auprc_without_deletion']:.6f})"
    )

    shift = reconstructed["gain"] - published["gain"]
    rho_with = rpm.spearman(published["s_with"], reconstructed["s_with"])
    rho_without = rpm.spearman(published["s_without"], reconstructed["s_without"])
    gain_ok = abs(shift) <= MAX_ABS_GAIN_SHIFT
    rank_ok = rho_with is not None and rho_with >= MIN_SCORE_SPEARMAN
    passes = bool(gain_ok and rank_ok)

    verdict = {
        "gain_published": published["gain"],
        "gain_reconstructed": reconstructed["gain"],
        "gain_shift": shift,
        "max_abs_gain_shift": MAX_ABS_GAIN_SHIFT,
        "gain_shift_passes": gain_ok,
        "score_spearman_deletion_arm": rho_with,
        "score_spearman_no_deletion_arm": rho_without,
        "min_score_spearman": MIN_SCORE_SPEARMAN,
        "score_spearman_passes": rank_ok,
        "passes": passes,
        "rule": "both conditions must hold; either one failing fails the gate",
        "pairs": len(heldout),
        "positives": sum(p.regulated for p in heldout),
        "auprc_published": {
            "with_deletion": published["auprc_with_deletion"],
            "without_deletion": published["auprc_without_deletion"],
        },
        "auprc_reconstructed": {
            "with_deletion": reconstructed["auprc_with_deletion"],
            "without_deletion": reconstructed["auprc_without_deletion"],
        },
    }

    payload = {
        "status": "AstroREG-2's endpoint gate, run once on the K562 held-out pairs. No AlphaGenome "
        "request, no money, no astrocyte value computed, no BAM written to disk",
        "lane": "lane-astro",
        "registration": {"path": str(REGISTRATION), "gate": "registered before either gain existed"},
        "what_was_held_fixed": "the weights, fitted once on the K562 training pairs with the PUBLISHED "
        "columns and reused for both runs; and the deletion features, which come from the same cache "
        "either way. Only the held-out pairs' activity differs",
        "no_rescaling": "the reconstructed columns entered exactly as computed; no factor was fitted",
        "reconstruction": {
            "rule": "amendment 2's, read from the producer's source",
            "producer": rpm.PRODUCER,
            "bams": recon["bams"],
            "streamed_gb": round(recon["streamed_bytes"] / 1e9, 2),
            "bytes_written_to_disk": 0,
        },
        "verdict": verdict,
        "reading": (
            "the endpoint is unchanged within the registered tolerance: the reconstructed activity is "
            "interchangeable with the published activity FOR THIS ENDPOINT, which is what the claim "
            "depends on. This does NOT mean the columns agree -- they do not, and the gate that "
            "measured that failed and stays failed"
            if passes
            else "the endpoint MOVES by more than the registered tolerance allows, or the ranking "
            "disagrees. The reconstructed activity is NOT interchangeable with the published activity "
            "for this endpoint, the route CLOSES, and no astrocyte value may be computed from it"
        ),
        "the_old_gate_is_not_reopened": "the scale gate failed and stays failed. Nothing here rescales "
        "anything, and this result must never be described as that test passing",
        "not_claimed": [
            "not a demonstration that the reconstructed columns equal the published ones: they do not",
            "not a licence to use the reconstruction for anything but this endpoint: a quantity that "
            "leaves one AUPRC difference unchanged is not validated for a threshold, a percentile or "
            "any claim about activity itself",
            "not an astrocyte result: no astrocyte value was read or computed",
            "not a spending approval: Albert's named the original registration by hash",
        ],
        "requests_sent": 0,
        "money_spent": 0,
    }
    payload[mf.KEY] = manifest(recon, heldout, training, verdict)
    path = save_result(RESULT, payload)

    print()
    print(json.dumps({k: v for k, v in verdict.items() if not k.startswith("auprc_")}, indent=1, default=str))
    print(f"-> {path}")
    print("GATE PASSES" if passes else "GATE FAILS: the route closes")
    print(f"({time.time() - t0:.0f} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
