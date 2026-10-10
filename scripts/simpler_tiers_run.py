# SPDX-License-Identifier: AGPL-3.0-or-later
"""Simpler genomes, stage 2: the worm's constraint tiers by range read, registered, then run once.

    uv run python scripts/simpler_tiers_run.py tiers      # cut per track + per-element counts, no RNA read
    uv run python scripts/simpler_tiers_run.py register   # census + registration; needs the tiers committed
    uv run python scripts/simpler_tiers_run.py score      # refuses unless the registration is committed

The phyloP bigWigs are read by HTTP range request and never stored. No model request is sent and no
key is read. docs/SIMPLER-GENOMES.md (stage 2) holds the reasoning.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import simpler  # noqa: E402
from genomeos.attribution import simpler_tiers as st  # noqa: E402
from genomeos.attribution.bigwig import BigWig  # noqa: E402
from genomeos.results import save_result  # noqa: E402

TIERS = "simpler_worm_tiers"
REGISTRATION = "simpler_worm_tiers_registration"
RESULT = "simpler_worm_tiers_result"
TIERS_FILE = f"data/results/{TIERS}.json"
REGISTRATION_FILE = f"data/results/{REGISTRATION}.json"
STAGE1_REGISTRATION = "data/results/simpler_worm_registration.json"
ENTRY = "scripts/simpler_tiers_run.py"
OWN_CODE = (
    "genomeos/attribution/simpler_tiers.py",
    "scripts/simpler_tiers_run.py",
    "tests/test_simpler_tiers.py",
)
TABLE = str(simpler.TABLE_S2)
GFF = str(simpler.GFF3)
#: Range-read budget per track (bytes moved, nothing stored), and the free disk the run needs.
BYTES_BUDGET = 200_000_000
MIN_FREE_GIB = 15

PARAMETERS = {
    "fdr_q": st.FDR_Q,
    "hist_bin": st.HIST_BIN,
    "hist_limit": st.HIST_LIMIT,
    "sample_group_blocks": st.SAMPLE_GROUP,
    "sample_fraction": st.SAMPLE_FRACTION,
    "cut_bootstrap": st.CUT_BOOTSTRAP,
    "constrained_min_fraction": st.CONSTRAINED_MIN,
    "primary_track": st.PRIMARY,
    "secondary_track": st.SECONDARY,
    "element_max_gap_bytes": st.ELEMENT_MAX_GAP,
    "workers": st.WORKERS,
    "min_per_tier": st.MIN_PER_TIER,
    "permutations": st.PERMUTATIONS,
    "bootstrap": st.BOOTSTRAP,
    "cluster_bp": simpler.CLUSTER_BP,
    "distance_bins_bp": list(simpler.DISTANCE_BINS),
    "detection_tpm": simpler.DETECTION_TPM,
    "seed": st.SEED,
    "bytes_budget_per_track": BYTES_BUDGET,
}
EXCLUSIONS = [
    "elements of any class other than putative_enhancer, and elements whose tissue call is not exactly "
    "one of the five tissues (stage 1's population)",
    "exon bases of any transcript in the WBcel235 GFF3, removed from every element before its tier",
    "elements with no scored base left after exon removal (not_measured): no tier, excluded and counted",
    "chrM and every block outside chrI-V and chrX from the genome sample",
    "the source's own element-to-gene columns, never read (simpler.NEVER_READ)",
]


def _sources() -> list[dict[str, Any]]:
    out = [
        {
            "accession": simpler.SOURCE,
            "version": "Genome Res 30(12), December 2020; supplement as served by Europe PMC on 2026-10-10",
        },
        {
            "accession": "Ensembl gene models for WBcel235 (WBcel235.63.gff3.gz)",
            "version": "release 63, the copy held in data/reference/celegans since 2026-09-10",
        },
        {
            "accession": (
                "Christmas et al. 2023, Science 380:eabn3943 (Zoonomia), Methods: phyloP p values to q "
                "values by an FDR correction, q <= 0.05 constrained or accelerated by the score's sign"
            ),
            "version": "the procedure only; no Zoonomia number is used",
        },
    ]
    for name, t in st.TRACKS.items():
        out.append(
            {
                "accession": f"UCSC ce11 {name}: {t['url']}",
                "version": t.get("served", "range-read 2026-10-10"),
            }
        )
    return out


def _manifest(inputs: list[tuple[str, str]], partitions: Any) -> dict:
    return {
        "sources": _sources(),
        "inputs": [mf.input_entry(path, partition=part) for path, part in inputs],
        "assembly": "WBcel235 (ce11): Table S2, the gene models and both phyloP tracks are on the same build",
        "coordinates": {
            "base": 0,
            "interval": "half-open",
            "note": "bigWig is 0-based half-open; the GFF3 (1-based closed) is converted on reading; Table "
            "S2 does not state its base and its start/end are taken as 0-based half-open, at most 1 bp off",
        },
        "parameters": PARAMETERS,
        "exclusions": EXCLUSIONS,
        "partitions": partitions,
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def _committed(rel: str) -> bool:
    logged = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", rel], cwd=ROOT, capture_output=True, text=True
    ).stdout.strip()
    clean = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", rel], cwd=ROOT).returncode == 0
    return bool(logged) and clean


def _commit_of(rel: str) -> str:
    return subprocess.run(
        ["git", "log", "-1", "--format=%h", "--", rel], cwd=ROOT, capture_output=True, text=True
    ).stdout.strip()


def _population() -> list[dict[str, Any]]:
    chroms = simpler.read_genes(simpler.GFF3)
    elements = simpler.read_elements(simpler.TABLE_S2, kinds=(simpler.ELEMENT_TYPE,))
    rows = st.population(chroms, elements)
    stage1 = json.loads(Path(STAGE1_REGISTRATION).read_text())["census"]["elements"]
    if len(rows) != stage1:
        sys.exit(f"refused: population {len(rows)} differs from stage 1's registered {stage1}")
    return rows


def _served(url: str) -> str:
    import urllib.request

    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "GenomeOS/0.9"})
    with urllib.request.urlopen(req, timeout=60) as r:  # noqa: S310
        return f"{r.headers.get('Content-Length')} bytes, Last-Modified {r.headers.get('Last-Modified')}"


def tiers() -> None:
    free = shutil.disk_usage(ROOT).free / 2**30
    if free < MIN_FREE_GIB:
        sys.exit(f"refused: {free:.1f} GiB free, under {MIN_FREE_GIB} GiB")
    mf.trace_begin()
    rows = _population()
    ann = st.read_annotation(simpler.GFF3)
    elems = [(r["chrom"], r["start"], r["end"]) for r in rows]
    per_track: dict[str, Any] = {}
    counts: dict[str, list[tuple[int, int]]] = {}
    for name in (st.PRIMARY, st.SECONDARY):
        t = st.TRACKS[name]
        t["served"] = _served(t["url"])
        t0 = time.time()
        bw = BigWig(t["url"])

        def guard(done: int, total: int, nbytes: int, name: str = name, bw: BigWig = bw) -> None:
            moved = bw.src.bytes_fetched
            if moved > BYTES_BUDGET:
                raise SystemExit(f"stopped: {name} moved {moved} bytes, over the {BYTES_BUDGET} budget")
            if done == total or done % 100 == 0:
                print(
                    f"  {name} {done}/{total} data {nbytes / 1e6:.1f} MB, moved {moved / 1e6:.1f} MB",
                    flush=True,
                )

        sample = st.genome_sample(bw, ann, progress=guard)
        cut = st.derive_cut(sample)
        print(name, json.dumps({k: v for k, v in cut.items() if k != "by_class"}), flush=True)
        if cut["cut"] is None:
            sys.exit(f"refused: no 5% FDR cut on {name}")
        c, ebytes = st.element_counts(bw, elems, ann, cut["cut"], progress=guard)
        counts[name] = c
        sample_keys = ("blocks_total", "groups_total", "groups_sampled", "blocks_sampled", "bases_sampled")
        per_track[name] = {
            "url": t["url"],
            "served": t["served"],
            "aligned": t["aligned"],
            "note": t["note"],
            **cut,
            "sample": {k: sample[k] for k in (*sample_keys, "data_bytes")},
            "abs_hist_by_class": {
                cls: {"pos": sample["pos"][k].tolist(), "neg": sample["neg"][k].tolist()}
                for k, cls in enumerate(st.CLASSES)
            },
            "element_data_bytes": ebytes,
            "bytes_moved_total": bw.src.bytes_fetched,
            "requests": bw.src.requests,
            "seconds": round(time.time() - t0, 1),
        }
        bw.close()
    stored = [
        [r["chrom"], r["start"], r["end"], *counts[st.PRIMARY][i], *counts[st.SECONDARY][i]]
        for i, r in enumerate(rows)
    ]
    out = {
        "result": TIERS,
        "date": time.strftime("%Y-%m-%d"),
        "rna_sheet_read": False,
        "what": (
            "Per track: the 5% FDR phyloP cut re-derived by Zoonomia's procedure on a uniform 10% sample "
            "of the track's data blocks, with |score| histograms by sign and by base class; per element of "
            "stage 1's population: exon-free scored bases and bases at or above the cut, per track."
        ),
        "tracks": per_track,
        "element_columns": list(st.ELEMENT_COLUMNS),
        "elements": stored,
        "bytes_moved_all_tracks": sum(v["bytes_moved_total"] for v in per_track.values()),
        "bytes_stored": 0,
        "alphagenome_requests": 0,
        "money": "none",
    }
    out["result_manifest"] = _manifest(
        [
            (TABLE, "elements and their tissue call only (sheet ATAC_metrics)"),
            (GFF, "gene positions and strands; exon, CDS and gene spans"),
            (STAGE1_REGISTRATION, "stage 1's census, read only to check the population count"),
        ],
        {"rna": "the RNA sheet (RNA_metrics) is not opened by this step"},
    )
    print(save_result(TIERS, out, compact=True))


def _tiered(rows: list[dict[str, Any]], stored: list[list[Any]]):
    prim = st.attach_tiers(rows, stored, st.PRIMARY)
    sec = st.attach_tiers(rows, stored, st.SECONDARY)
    edges = st.length_edges([r["measured"] for r in prim])
    return prim, sec, edges


CUT_KEYS = (
    "cut",
    "cut_interval_95",
    "constrained_share_all_scored",
    "constrained_share_non_exon",
    "instrument_check",
)


def register() -> None:
    if not _committed(TIERS_FILE):
        sys.exit(f"refused: {TIERS_FILE} is not committed and clean")
    mf.trace_begin()
    tiers_doc = json.loads(Path(TIERS_FILE).read_text())
    rows = _population()
    prim, sec, edges = _tiered(rows, tiers_doc["elements"])
    cen = st.census(prim, edges)
    cen_sec = st.census(sec, edges)
    tracks = tiers_doc["tracks"]
    out = {
        "result": REGISTRATION,
        "date": time.strftime("%Y-%m-%d"),
        "written_before_tiers_were_read_against_rna": True,
        "registered_inputs": {
            "table_s2": {"path": TABLE, "sha256": simpler.sha256(simpler.TABLE_S2)},
            "gff3": {"path": GFF, "sha256": simpler.sha256(simpler.GFF3)},
            "tiers": {
                "path": TIERS_FILE,
                "sha256": simpler.sha256(Path(TIERS_FILE)),
                "commit": _commit_of(TIERS_FILE),
            },
        },
        "cuts": {name: {k: tracks[name].get(k) for k in CUT_KEYS} for name in (st.PRIMARY, st.SECONDARY)},
        "registration": st.REGISTRATION,
        "census": cen,
        "census_secondary_beside": cen_sec,
        "alphagenome_requests": 0,
        "money": "none",
    }
    out["result_manifest"] = _manifest(
        [
            (TIERS_FILE, "the committed tiers (no RNA read)"),
            (TABLE, "elements and their tissue call only (sheet ATAC_metrics)"),
            (GFF, "gene positions and strands"),
            (STAGE1_REGISTRATION, "stage 1's census, read only to check the population count"),
        ],
        {"population": "fixed here, before the RNA sheet (RNA_metrics) was read against any tier"},
    )
    print(save_result(REGISTRATION, out))
    print(json.dumps({"cuts": out["cuts"], "census": cen}, indent=1))


def score() -> None:
    for rel in (REGISTRATION_FILE, TIERS_FILE, "genomeos/attribution/simpler_tiers.py"):
        if not _committed(rel):
            sys.exit(f"refused: {rel} is not committed and clean")
    mf.trace_begin()
    reg = json.loads(Path(REGISTRATION_FILE).read_text())
    if reg["registration"] != json.loads(json.dumps(st.REGISTRATION)):
        sys.exit("refused: simpler_tiers.REGISTRATION differs from the committed registration")
    for key in ("table_s2", "gff3", "tiers"):
        if (
            simpler.sha256(Path(reg["registered_inputs"][key]["path"]))
            != reg["registered_inputs"][key]["sha256"]
        ):
            sys.exit(f"refused: {key} differs from the registered input")
    tiers_doc = json.loads(Path(TIERS_FILE).read_text())
    rows = _population()
    prim, sec, edges = _tiered(rows, tiers_doc["elements"])
    cen = st.census(prim, edges)
    if cen != reg["census"]:
        sys.exit("refused: the census differs from the registered census")
    instrument = bool(reg["cuts"][st.PRIMARY]["instrument_check"]["passed"])
    expr = simpler.read_expression(simpler.TABLE_S2)
    res = st.score(prim, sec, expr, edges, cen, instrument)
    out = {
        "result": RESULT,
        "date": time.strftime("%Y-%m-%d"),
        **res,
        "registration_commit": _commit_of(REGISTRATION_FILE),
        "cuts": reg["cuts"],
        "census": cen,
        "alphagenome_requests": 0,
        "money": "none",
    }
    out["result_manifest"] = _manifest(
        [
            (REGISTRATION_FILE, "the committed registration"),
            (TIERS_FILE, "the committed tiers"),
            (TABLE, "elements (ATAC_metrics) and the endpoint (RNA_metrics)"),
            (GFF, "gene positions and strands"),
            (STAGE1_REGISTRATION, "stage 1's census, read only to check the population count"),
        ],
        {"endpoint": "RNA_metrics read once, after the registration was committed"},
    )
    print(save_result(RESULT, out))
    print(out["outcome"], out.get("falsifier", ""))
    print(out["reading"])


if __name__ == "__main__":
    os.chdir(ROOT)
    {"tiers": tiers, "register": register, "score": score}[sys.argv[1]]()
