# SPDX-License-Identifier: AGPL-3.0-or-later
"""Simpler genomes, stage 3: the worm's constraint-dose view registered on a second source, then run once.

    uv run python scripts/simpler_dose_run.py fetch       # Jänes 2018 source data + GEO long-cap tracks
    uv run python scripts/simpler_dose_run.py inventory   # data/results/simpler_dose_inventory.json
    uv run python scripts/simpler_dose_run.py tiers       # phyloP135way counts per element, by range read
    uv run python scripts/simpler_dose_run.py register    # census + registration; no RNA value read
    uv run python scripts/simpler_dose_run.py score       # refuses unless the registration is committed
    uv run python scripts/simpler_dose_run.py register2   # after run 1 stopped at its strand check
    uv run python scripts/simpler_dose_run.py score2      # registration 2, run once
    uv run python scripts/simpler_dose_run.py edges_posthoc   # post hoc, deciding nothing: edge rounding

The phyloP track is read by HTTP range request and never stored. No model request is sent and no key
is read. docs/SIMPLER-GENOMES.md section 10 holds the reasoning.
"""

from __future__ import annotations

import datetime as dt
import json
import os
import shutil
import subprocess
import sys
import time
import urllib.request
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import simpler  # noqa: E402
from genomeos.attribution import simpler_dose as sd  # noqa: E402
from genomeos.attribution import simpler_tiers as st  # noqa: E402
from genomeos.attribution.bigwig import BigWig  # noqa: E402
from genomeos.results import save_result  # noqa: E402

INVENTORY = "simpler_dose_inventory"
TIERS = "simpler_dose_tiers"
REGISTRATION = "simpler_dose_registration"
RESULT = "simpler_dose_result"
TIERS_FILE = f"data/results/{TIERS}.json"
REGISTRATION_FILE = f"data/results/{REGISTRATION}.json"
RESULT_FILE = f"data/results/{RESULT}.json"
#: Registration 2: registration 1 ran once and stopped at its strand check (see sd.REGISTRATION_V2).
REGISTRATION2 = "simpler_dose_v2_registration"
RESULT2 = "simpler_dose_v2_result"
REGISTRATION2_FILE = f"data/results/{REGISTRATION2}.json"
RESULT2_FILE = f"data/results/{RESULT2}.json"
#: Post hoc, deciding nothing: run 2's primary with the census's unrounded tertile edges.
EDGES_POSTHOC = "simpler_dose_v2_edges_posthoc"
STAGE2_TIERS = "data/results/simpler_worm_tiers.json"
ENTRY = "scripts/simpler_dose_run.py"
OWN_CODE = (
    "genomeos/attribution/simpler.py",
    "genomeos/attribution/simpler_tiers.py",
    "genomeos/attribution/simpler_dose.py",
    "scripts/simpler_dose_run.py",
    "tests/test_simpler_dose.py",
)
GFF = str(simpler.GFF3)
DOWNLOAD_BUDGET = 1_500_000_000
RANGE_BUDGET = 200_000_000
MIN_FREE_GIB = 15

#: Bytes moved while choosing the source, outside the fetch: the eLife API record of the article.
#: HEAD requests and directory listings are not counted.
SCOUTING_BYTES = 272_291

CONSIDERED = [
    {
        "source": "Jänes et al. 2018, eLife 7:e37344 (source data) and GEO GSE114483 (long cap RNA)",
        "what": (
            "ATAC of whole wild-type animals at six stages with typed elements, a dynamics call and "
            "per-stage heights; strand-specific long nuclear RNA at the same six stages, same study"
        ),
        "assembly": "ce11 (source data give ce10 and ce11; the GEO tracks are ce11)",
        "licence": sd.LICENCE,
        "used": True,
    },
    {
        "source": "Durham et al. 2021, Genome Res 31:1952",
        "what": "single-cell ATAC at L2 by cell type",
        "assembly": "ce11",
        "licence": "not checked",
        "used": False,
        "why_not": (
            "no RNA measured in the same samples: cell-type RNA would come from another study "
            "(Cao et al. 2017), so the endpoint would not be measured with the accessibility"
        ),
    },
    {
        "source": "Daugherty et al. 2017, Genome Res 27:2096",
        "what": "ATAC at embryo, L3 and young adult, reporter tests of selected elements",
        "assembly": "ce11",
        "licence": "not checked",
        "used": False,
        "why_not": "RNA in the same samples not established when checked; three stages against six",
    },
    {
        "source": "GEO GSE141213 (Serizay 2020's tracks)",
        "what": "per-tissue ATAC and RNA bigWigs",
        "assembly": "ce11",
        "licence": "public GEO record",
        "used": False,
        "why_not": (
            "2,386,268,160-byte archive, over the 1.5 GB budget, and the same source as stages 1 and 2"
        ),
    },
]

PARAMETERS = {
    "stages": list(sd.STAGES),
    "detection_coverage_tpm": sd.DETECTION_TPM,
    "constraint_track": sd.PRIMARY,
    "constraint_cut": sd.STAGE2_CUT,
    "constrained_min_fraction_for_the_5pct_tier": st.CONSTRAINED_MIN,
    "stage2_edges": list(sd.STAGE2_EDGES),
    "min_per_arm": sd.MIN_PER_ARM,
    "min_discordant_secondary": simpler.MIN_DISCORDANT,
    "permutations": st.PERMUTATIONS,
    "bootstrap": st.BOOTSTRAP,
    "control_draws": simpler.CONTROL_DRAWS,
    "cluster_bp": simpler.CLUSTER_BP,
    "distance_bins_bp": list(simpler.DISTANCE_BINS),
    "seed": simpler.SEED,
    "alpha": simpler.ALPHA,
    "range_budget_bytes": RANGE_BUDGET,
    "download_budget_bytes": DOWNLOAD_BUDGET,
}
EXCLUSIONS = [
    "elements of classes other than putative_enhancer (population) and coding_promoter (positive control)",
    "elements the source does not call developmentally dynamic (devel_is_dynamic False)",
    "the ageing (glp-1) samples and their columns",
    "exon bases of any WBcel235 transcript, removed from every element before its constrained fraction",
    "elements with no scored base left after exon removal: no fraction, out of the dose test, counted",
    "the source's own element-to-gene columns, never read (simpler_dose.NEVER_READ)",
    "the RNA-derived flag columns of Figure 2 source data (lcap_*_passed_*), never read",
]


def _sources() -> list[dict[str, Any]]:
    return [
        {
            "accession": sd.SOURCE,
            "version": "eLife source data v2 and GEO GSE114483, as served on 2026-10-10",
        },
        {
            "accession": "Ensembl gene models for WBcel235 (WBcel235.63.gff3.gz)",
            "version": "release 63, the copy held in data/reference/celegans since 2026-09-10",
        },
        {"accession": f"UCSC ce11 {sd.PRIMARY}: {sd.PRIMARY_URL}", "version": "range-read 2026-10-10"},
    ]


def _manifest(inputs: list[tuple[str, str]], partitions: Any) -> dict:
    return {
        "sources": _sources(),
        "inputs": [mf.input_entry(path, partition=part) for path, part in inputs],
        "assembly": "WBcel235 (ce11): the source's ce11 columns, the GEO tracks, the gene models and phyloP",
        "coordinates": {
            "base": 0,
            "interval": "half-open",
            "note": (
                "the source data are bed-style (0-based, half-open); bigWig is 0-based half-open; the GFF3 "
                "(1-based closed) is converted on reading"
            ),
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


def _free_gib() -> float:
    return shutil.disk_usage(ROOT).free / 2**30


# --------------------------------------------------------------------------------------------- fetch


def _urls() -> dict[Path, str]:
    out = {p: f"{sd.ELIFE}/{p.name}" for p in sd.ELIFE_FILES}
    for p in sd.LCAP.values():
        out[p] = f"{sd.GEO}/{p.name}"
    return out


def fetch() -> None:
    free = _free_gib()
    if free < MIN_FREE_GIB:
        sys.exit(f"refused: {free:.1f} GiB free, under {MIN_FREE_GIB} GiB")
    sd.STORE.mkdir(parents=True, exist_ok=True)
    rec: dict[str, Any] = {
        "fetched_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
        "disk_free_bytes_before": shutil.disk_usage(ROOT).free,
        "files": {},
    }
    moved = 0
    for path, url in _urls().items():
        if not path.exists():
            req = urllib.request.Request(url, headers={"User-Agent": "GenomeOS/1.0"})
            with urllib.request.urlopen(req, timeout=600) as r:  # noqa: S310
                body = r.read()
            moved += len(body)
            if moved > DOWNLOAD_BUDGET:
                sys.exit("refused: over the download budget")
            path.write_bytes(body)
        rec["files"][path.name] = {
            "url": url,
            "bytes": path.stat().st_size,
            "sha256": simpler.sha256(path),
            "licence": sd.LICENCE,
        }
        print(path.name, rec["files"][path.name]["bytes"], flush=True)
    rec["download_bytes"] = sum(f["bytes"] for f in rec["files"].values())
    rec["disk_free_bytes_after"] = shutil.disk_usage(ROOT).free
    sd.FETCH_RECORD.write_text(json.dumps(rec, indent=1) + "\n")
    print(json.dumps({k: v for k, v in rec.items() if k != "files"}, indent=1))


def inventory() -> None:
    mf.trace_begin()
    rec = json.loads(sd.FETCH_RECORD.read_text())
    _, funnel = sd.read_elements()
    out = {
        "question": (
            "a second worm source with element accessibility and stage or tissue RNA measured by the same "
            "study: what was fetched, its size, hash and licence, and why this source"
        ),
        "date": time.strftime("%Y-%m-%d"),
        "fetched": [{"name": k, **v} for k, v in rec["files"].items()],
        "fetched_bytes": rec["download_bytes"],
        "scouting_bytes": SCOUTING_BYTES,
        "lane_download_bytes": rec["download_bytes"] + SCOUTING_BYTES,
        "budget_bytes": DOWNLOAD_BUDGET,
        "disk_free_bytes_before": rec["disk_free_bytes_before"],
        "disk_free_bytes_after": rec["disk_free_bytes_after"],
        "element_funnel_no_rna": funnel,
        "considered": CONSIDERED,
        "rna_read": False,
        "money": "none; no model request",
    }
    out["result_manifest"] = _manifest(
        [
            (str(sd.FETCH_RECORD), "the fetch record"),
            (str(sd.FIG1), "coordinates and the six wild-type ATAC heights"),
            (str(sd.FIG2), "coordinates (ce10, ce11) and the class only"),
            (str(sd.FIG4), "coordinates and the developmental-dynamics call only"),
        ],
        "n/a: an inventory has no evaluation partition",
    )
    print(save_result(INVENTORY, out))


# --------------------------------------------------------------------------------------------- tiers


def _population() -> tuple[dict, list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    chroms = simpler.read_genes(simpler.GFF3)
    els, funnel = sd.read_elements()
    rows = sd.attribute(els, chroms)
    enh = [r for r in rows if r["kind"] == sd.ELEMENT_TYPE]
    prom = [r for r in rows if r["kind"] == sd.PROMOTER_TYPE]
    return chroms, enh, prom, funnel


def tiers() -> None:
    free = _free_gib()
    if free < MIN_FREE_GIB:
        sys.exit(f"refused: {free:.1f} GiB free, under {MIN_FREE_GIB} GiB")
    mf.trace_begin()
    stage2 = json.loads(Path(STAGE2_TIERS).read_text())
    s2 = stage2["tracks"][sd.PRIMARY]
    if s2["cut"] != sd.STAGE2_CUT or not s2["instrument_check"]["passed"]:
        sys.exit(
            "refused: stage 2's committed cut or instrument check differs from the one this module names"
        )
    _, enh, _, funnel = _population()
    ann = st.read_annotation(simpler.GFF3)
    elems = [(r["chrom"], r["start"], r["end"]) for r in enh]
    t0 = time.time()
    bw = BigWig(sd.PRIMARY_URL)

    def guard(done: int, total: int, nbytes: int) -> None:
        moved = bw.src.bytes_fetched
        if moved > RANGE_BUDGET:
            raise SystemExit(f"stopped: moved {moved} bytes, over the {RANGE_BUDGET} budget")
        print(f"  {done}/{total} data {nbytes / 1e6:.1f} MB, moved {moved / 1e6:.1f} MB", flush=True)

    counts, ebytes = st.element_counts(bw, elems, ann, sd.STAGE2_CUT, progress=guard)
    out = {
        "result": TIERS,
        "date": time.strftime("%Y-%m-%d"),
        "rna_read": False,
        "what": (
            "Per element of the population (Jänes 2018 putative enhancers called developmentally dynamic, "
            "with a domain-rule gene): exon-free scored bases of UCSC ce11 phyloP135way and bases at or "
            "above stage 2's cut 2.08, read by range request."
        ),
        "track": {
            "name": sd.PRIMARY,
            "url": sd.PRIMARY_URL,
            "cut": sd.STAGE2_CUT,
            "cut_source": f"{STAGE2_TIERS} (commit {_commit_of(STAGE2_TIERS)})",
            "stage2_instrument_check": s2["instrument_check"],
            "element_data_bytes": ebytes,
            "bytes_moved_total": bw.src.bytes_fetched,
            "requests": bw.src.requests,
            "seconds": round(time.time() - t0, 1),
        },
        "element_funnel": funnel,
        "element_columns": ["chrom", "start", "end", "measured", "above"],
        "elements": [[r["chrom"], r["start"], r["end"], *counts[i]] for i, r in enumerate(enh)],
        "bytes_stored": 0,
        "alphagenome_requests": 0,
        "money": "none",
    }
    bw.close()
    out["result_manifest"] = _manifest(
        [
            (str(sd.FIG1), "coordinates and the six wild-type ATAC heights"),
            (str(sd.FIG2), "coordinates (ce10, ce11) and the class only"),
            (str(sd.FIG4), "coordinates and the developmental-dynamics call only"),
            (GFF, "gene positions and strands; exon spans"),
            (STAGE2_TIERS, "stage 2's committed cut and instrument check only"),
        ],
        {"rna": "no RNA file is opened by this step"},
    )
    print(save_result(TIERS, out, compact=True))


# --------------------------------------------------------------------------------------------- register


def _build() -> tuple[dict, list[dict[str, Any]], list[dict[str, Any]], list[bool], dict[str, Any], dict]:
    tiers_doc = json.loads(Path(TIERS_FILE).read_text())
    chroms, enh, prom, funnel = _population()
    rows = sd.attach_counts(enh, tiers_doc["elements"])
    stage2 = json.loads(Path(STAGE2_TIERS).read_text())
    overlap = sd.overlaps_any(rows, [(e[0], e[1], e[2]) for e in stage2["elements"]])
    cen = sd.census(rows, prom, overlap)
    return chroms, rows, prom, overlap, cen, funnel


def _inputs_hashed() -> dict[str, dict[str, str]]:
    out = {
        "fig1": {"path": str(sd.FIG1), "sha256": simpler.sha256(sd.FIG1)},
        "fig2": {"path": str(sd.FIG2), "sha256": simpler.sha256(sd.FIG2)},
        "fig4": {"path": str(sd.FIG4), "sha256": simpler.sha256(sd.FIG4)},
        "gff3": {"path": GFF, "sha256": simpler.sha256(simpler.GFF3)},
        "tiers": {"path": TIERS_FILE, "sha256": simpler.sha256(Path(TIERS_FILE))},
        "stage2_tiers": {"path": STAGE2_TIERS, "sha256": simpler.sha256(Path(STAGE2_TIERS))},
    }
    for (s, d), p in sd.LCAP.items():
        out[f"lcap_{s}_{d}"] = {"path": str(p), "sha256": simpler.sha256(p)}
    return out


def _lcap_inputs(note: str) -> list[tuple[str, str]]:
    return [(str(p), note) for p in sd.LCAP.values()]


def register() -> None:
    for rel in (TIERS_FILE, STAGE2_TIERS):
        if not _committed(rel):
            sys.exit(f"refused: {rel} is not committed and clean")
    mf.trace_begin()
    _, _, _, _, cen, funnel = _build()
    out = {
        "result": REGISTRATION,
        "date": time.strftime("%Y-%m-%d"),
        "written_before_any_rna_value_was_read": True,
        "registered_inputs": _inputs_hashed(),
        "tiers_commit": _commit_of(TIERS_FILE),
        "element_funnel": funnel,
        "registration": sd.REGISTRATION,
        "census": cen,
        "alphagenome_requests": 0,
        "money": "none",
    }
    out["result_manifest"] = _manifest(
        [
            (TIERS_FILE, "the committed per-element counts (no RNA read)"),
            (STAGE2_TIERS, "stage 2's population coordinates, for the overlap count"),
            (str(sd.FIG1), "coordinates and the six wild-type ATAC heights"),
            (str(sd.FIG2), "coordinates (ce10, ce11) and the class only"),
            (str(sd.FIG4), "coordinates and the developmental-dynamics call only"),
            (GFF, "gene positions and strands"),
            *_lcap_inputs("hashed only, so the run can check them; no value is read here"),
        ],
        {"population": "fixed here, before any RNA value of this source was read"},
    )
    print(save_result(REGISTRATION, out))
    print(json.dumps(cen, indent=1))


def register2() -> None:
    """Registration 2, after registration 1's single run stopped at its strand check."""
    for rel in (REGISTRATION_FILE, RESULT_FILE, TIERS_FILE, STAGE2_TIERS):
        if not _committed(rel):
            sys.exit(f"refused: {rel} is not committed and clean")
    mf.trace_begin()
    reg1 = json.loads(Path(REGISTRATION_FILE).read_text())
    run1 = json.loads(Path(RESULT_FILE).read_text())
    if run1["outcome"] != "INSTRUMENT_FAILED":
        sys.exit("refused: registration 2 exists only because run 1 stopped at its strand check")
    _, _, _, _, cen, funnel = _build()
    if json.loads(json.dumps(cen)) != reg1["census"]:
        sys.exit("refused: the census does not recompute to registration 1's")
    out = {
        "result": REGISTRATION2,
        "date": time.strftime("%Y-%m-%d"),
        "written_before_any_gene_value_or_outcome_was_read": True,
        "registration_1_commit": _commit_of(REGISTRATION_FILE),
        "run_1_commit": _commit_of(RESULT_FILE),
        "run_1_outcome": run1["outcome"],
        "run_1_strand_check": run1["strand_check"],
        "rna_seen_before_this_registration": {
            "run_1_strand_check": "the twelve per-stage sums above",
            "value_signs_per_track": sd.sign_census(sd.LCAP),
        },
        "registered_inputs": _inputs_hashed(),
        "tiers_commit": _commit_of(TIERS_FILE),
        "element_funnel": funnel,
        "registration": sd.REGISTRATION_V2,
        "census": cen,
        "alphagenome_requests": 0,
        "money": "none",
    }
    out["result_manifest"] = _manifest(
        [
            (REGISTRATION_FILE, "registration 1, for its census and commit"),
            (RESULT_FILE, "run 1, for its outcome and strand sums"),
            (TIERS_FILE, "the committed per-element counts (no RNA read)"),
            (STAGE2_TIERS, "stage 2's population coordinates, for the overlap count"),
            (str(sd.FIG1), "coordinates and the six wild-type ATAC heights"),
            (str(sd.FIG2), "coordinates (ce10, ce11) and the class only"),
            (str(sd.FIG4), "coordinates and the developmental-dynamics call only"),
            (GFF, "gene positions and strands"),
            *_lcap_inputs("hashed, and the sign of each value counted per track; no gene value formed"),
        ],
        {"population": "registration 1's, recomputed and checked"},
    )
    print(save_result(REGISTRATION2, out))
    print(json.dumps(out["rna_seen_before_this_registration"], indent=1))


# --------------------------------------------------------------------------------------------- score


def score(version: int = 1) -> None:
    reg_name, reg_file, result_name, registration, magnitude = {
        1: (REGISTRATION, REGISTRATION_FILE, RESULT, sd.REGISTRATION, False),
        2: (REGISTRATION2, REGISTRATION2_FILE, RESULT2, sd.REGISTRATION_V2, True),
    }[version]
    before = (RESULT_FILE,) if version == 2 else ()
    for rel in (reg_file, *before, TIERS_FILE, STAGE2_TIERS, *OWN_CODE):
        if not _committed(rel):
            sys.exit(f"refused: {rel} must be committed and unchanged before any RNA value is read")
    if version == 2 and Path(f"data/results/{RESULT2}.json").exists():
        sys.exit("refused: registration 2 has been run; it runs once")
    mf.trace_begin()
    reg = json.loads(Path(reg_file).read_text())
    if reg["registration"] != json.loads(json.dumps(registration)):
        sys.exit(f"refused: the module's registration differs from the committed {reg_name}")
    for key, entry in reg["registered_inputs"].items():
        if simpler.sha256(Path(entry["path"])) != entry["sha256"]:
            sys.exit(f"refused: {key} differs from the registered input")
    chroms, rows, prom, overlap, cen, _ = _build()
    if json.loads(json.dumps(cen)) != reg["census"]:
        sys.exit("refused: the census does not recompute to the registered one")
    t0 = time.time()
    exons = sd.read_gene_exons(simpler.GFF3)
    cov = sd.coverage(sd.LCAP, exons, magnitude=magnitude)
    strand = sd.orientation(cov, exons)
    genes = {g.gid for ch in chroms.values() for g in ch.genes}
    expr = sd.expression(cov, exons, strand["verdict"], genes)
    res = sd.score(rows, prom, overlap, chroms, expr, reg["census"], strand)
    out = {
        "result": result_name,
        "date": time.strftime("%Y-%m-%d"),
        "registration_version": version,
        "registration_commit": _commit_of(reg_file),
        "coverage_read_as_magnitude": magnitude,
        **res,
        "rna_summary": {
            "genes_with_value": len(expr),
            "protein_coding_genes": len(genes),
            "seconds": round(time.time() - t0, 1),
        },
        "census": reg["census"],
        "reading_is_the_registered_wording": (
            "the reading is REGISTRATION['outcome_wording'] with its numbers filled in "
            "and no other word changed"
        ),
        "alphagenome_requests": 0,
        "money": "none",
    }
    out["result_manifest"] = _manifest(
        [
            (reg_file, "the committed registration"),
            *([(RESULT_FILE, "run 1, which must be committed before run 2")] if version == 2 else []),
            (TIERS_FILE, "the committed per-element counts"),
            (STAGE2_TIERS, "stage 2's population coordinates, for the overlap"),
            (str(sd.FIG1), "coordinates and the six wild-type ATAC heights"),
            (str(sd.FIG2), "coordinates (ce10, ce11) and the class only"),
            (str(sd.FIG4), "coordinates and the developmental-dynamics call only"),
            (GFF, "gene positions, strands and exons"),
            *_lcap_inputs("the endpoint: long-cap linear coverage, read once after the registration"),
        ],
        {"endpoint": "the long-cap tracks, read once, after the registration was committed"},
    )
    print(save_result(result_name, out))
    print(out["outcome"], out.get("falsifier", ""))
    print(out["reading"])


def edges_posthoc() -> None:
    """Post hoc, deciding nothing. Run 2 read the tertile edges as the census stored them (rounded to
    six decimals) while the census counted its arms on the exact edges; this computes run 2's primary
    both ways, so the disclosure in the doc rests on a committed result rather than a scratch run."""
    for rel in (REGISTRATION2_FILE, RESULT2_FILE, TIERS_FILE, STAGE2_TIERS, *OWN_CODE):
        if not _committed(rel):
            sys.exit(f"refused: {rel} must be committed and unchanged")
    mf.trace_begin()
    reg = json.loads(Path(REGISTRATION2_FILE).read_text())
    run2 = json.loads(Path(RESULT2_FILE).read_text())
    chroms, rows, _, _, cen, _ = _build()
    if json.loads(json.dumps(cen)) != reg["census"]:
        sys.exit("refused: the census does not recompute to the registered one")
    exact = st.fraction_edges(rows)
    exons = sd.read_gene_exons(simpler.GFF3)
    cov = sd.coverage(sd.LCAP, exons, magnitude=True)
    strand = sd.orientation(cov, exons)
    genes = {g.gid for ch in chroms.values() for g in ch.genes}
    expr = sd.expression(cov, exons, strand["verdict"], genes)
    for r in rows:
        r["c_domain"] = int(sd.concordant(expr.get(r["domain_gene"]), r["tissue"]))
    ledges = cen["length_edges_measured_bases"]
    stored = st.matched_difference(sd.by_edges(rows, cen["fraction_tertile_edges"]), ledges)
    unrounded = st.matched_difference(sd.by_edges(rows, exact), ledges)
    if stored["delta_points"] != run2["primary"]["delta_points"]:
        sys.exit("refused: the stored-edge primary does not reproduce run 2")
    out = {
        "result": EDGES_POSTHOC,
        "date": time.strftime("%Y-%m-%d"),
        "post_hoc_deciding_nothing": True,
        "why": (
            "Run 2 (registration 5b930cb) read the tertile edges as the census stored them, rounded to six "
            "decimals; the census counted its arms on the exact edges. The stored lower edge 0.072848 is "
            "above the exact 11/151 and the stored upper edge 0.19403 above the exact 13/67, so elements "
            "whose fraction equals either exact edge sit in a different arm. The registered outcome is "
            "run 2's; this computes the primary both ways and decides nothing."
        ),
        "run_2_commit": _commit_of(RESULT2_FILE),
        "run_2_outcome": run2["outcome"],
        "edges_stored": cen["fraction_tertile_edges"],
        "edges_exact": exact,
        "stored_edges_as_run": stored,
        "exact_edges_post_hoc": unrounded,
        "exact_edges_label_if_it_had_been_registered": sd.dose_label(*unrounded["delta_ci95_points"]),
        "alphagenome_requests": 0,
        "money": "none",
    }
    out["result_manifest"] = _manifest(
        [
            (REGISTRATION2_FILE, "registration 2, its census"),
            (RESULT2_FILE, "run 2, to check the stored-edge primary reproduces"),
            (TIERS_FILE, "the committed per-element counts"),
            (STAGE2_TIERS, "stage 2's population coordinates, for the overlap"),
            (str(sd.FIG1), "coordinates and the six wild-type ATAC heights"),
            (str(sd.FIG2), "coordinates (ce10, ce11) and the class only"),
            (str(sd.FIG4), "coordinates and the developmental-dynamics call only"),
            (GFF, "gene positions, strands and exons"),
            *_lcap_inputs("the endpoint, read as magnitude, after run 2 was committed"),
        ],
        {"post_hoc": "after run 2 was committed; deciding nothing"},
    )
    print(save_result(EDGES_POSTHOC, out))
    for k in ("stored_edges_as_run", "exact_edges_post_hoc"):
        v = out[k]
        print(
            k,
            v["constrained_in_support"],
            v["not_constrained_in_support"],
            v["delta_points"],
            v["delta_ci95_points"],
        )


if __name__ == "__main__":
    os.chdir(ROOT)
    {
        "fetch": fetch,
        "inventory": inventory,
        "tiers": tiers,
        "register": register,
        "score": score,
        "register2": register2,
        "score2": lambda: score(2),
        "edges_posthoc": edges_posthoc,
    }[sys.argv[1]]()
