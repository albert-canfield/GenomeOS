# SPDX-License-Identifier: AGPL-3.0-or-later
"""Audit C, checkpoint 1: is an EN-TEx allele-coordination test feasible? Depth and counts only.

    uv run python scripts/entex_feasibility.py               # build the depth-only cache once, then count
    uv run python scripts/entex_feasibility.py --cache-only  # stop after the cache
    uv run python scripts/entex_feasibility.py --alphagenome-metadata PATH   # a saved output_metadata copy

The question (lane-entex, 2026-09-29; docs/DATA.md "EN-TEx allele coordination: the feasibility audit").
EN-TEx (Rozowsky et al. 2023, Cell) mapped each donor's ATAC-seq and RNA-seq reads to that donor's own
phased diploid genome and released haplotype read counts per regulatory region and per gene. A later,
separately registered endpoint could ask whether a frozen nominated target's expression leans to the same
haplotype as its element's accessibility more often than distance-, expression- and depth-matched
alternative genes do. That tests coordination, not causality. This script decides, before any signed
allelic outcome is opened, whether such an endpoint would have enough rows, loci and donors. It is a new
question with its own rationale, not a continuation of the stopped coherence pilot.

What it counts, from the C5 probe's two EN-TEx tables and four small public files:

1. donor-tissue samples with both ATAC and RNA allelic counts (not "any chromatin assay");
2. phase: whether an element and its target gene lie in one of the donor's published phased blocks, so the
   element's heterozygous site and the gene's can be put on the same haplotype;
3. total depth per element-gene row: ATAC reads at the element, RNA reads at the gene;
4. how many rows are balanced, per side, by EN-TEx's own unsigned call;
5. exclusions: mapping (ENCODE exclusion list v2), copy number (EN-TEx's SV calls, two donors), imprinting
   (geneimprint), immunoglobulin and T-cell receptor genes and the MHC;
6. matched alternative genes for each row, matched on distance, expression and depth;
7. exposure: where EN-TEx experiments, or anything derived from them, enter the project;
8. independence, stated, not resolved;

and the power sensitivity by depth floor, from total depth only.

The outcome discipline. The two tables carry, per row, hap1_count, hap2_count, hap1_allele_ratio,
p_betabinom and imbalance_significance beside the identifiers. `depth_only` is the only function that
sees a raw line: it keeps the identifiers, sums hap1_count and hap2_count on the line that parses them and
keeps only the sum, reduces imbalance_significance (EN-TEx's 0/1 call) to an unsigned flag, and never
indexes hap1_allele_ratio or p_betabinom. The raw bytes are streamed and hashed, never written to disk;
the cache holds only what `depth_only` returns. The output of `depth_only` is the same when hap1 and hap2
are swapped (tests/test_entex_feasibility.py), so no haplotype direction can leave it. The unsigned flag is
reported as marginals per side only: the joint count (both sides imbalanced in the same row) is not
formed, because against its expectation under independence it is an unsigned coordination statistic that
would pre-empt the endpoint. No flag is computed for the alternative genes.

The go/no-go rule was fixed in this file before the first run (GO_* below). 0 AlphaGenome requests.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import re
import resource
import subprocess
import sys
import tarfile
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.results import save_result  # noqa: E402

_spec = importlib.util.spec_from_file_location("c5_probe", ROOT / "scripts" / "c5_paired_variation_probe.py")
c5 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(c5)  # imported, never edited: its loaders, streaming reader and read formula

RESULT = "entex_feasibility"
CACHE = Path("data/cache/entex")
DEPTH_CACHE = CACHE / "depth_only.tsv.gz"
DEPTH_SIDECAR = CACHE / "depth_only.json"
LEDGER = CACHE / "ledger.jsonl"  # every download and every run of this script, appended
ENTEX = c5.ENTEX
AS_TABLES = {"ccre": "cCREs_default_AS.tsv", "gene": "genes_default_AS.tsv"}
PHASED_URL = ENTEX + "phased_block.tar.gz"
#: EN-TEx names its donors ENC-001..ENC-004 in the AS tables and individual 1..4 in the phased-block and
#: SV files; the correspondence is taken by number
DONOR_INDIVIDUAL = {f"ENC-00{k}": f"ind{k}" for k in (1, 2, 3, 4)}
SV_URLS = {"ENC-001": ENTEX + "individual1_SV.vcf", "ENC-004": ENTEX + "individual4_SV.vcf"}
EXCLUSION_URL = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz"
#: AlphaGenome's track metadata (output_metadata) as a peer session saved it; re-fetching it is a
#: model-service call this lane may not make, so a saved copy is read as found, when one is given
AG_METADATA = CACHE / "alphagenome_track_metadata_copy.csv"
GTEX_TPM = c5.CACHE / Path(c5.GTEX_TPM_URL).name
#: files in which naming an EN-TEx allelic table is expected: this lane, the C5 probe, and prose
ALLOWED_TO_NAME = (
    "scripts/entex_feasibility.py",
    "tests/test_entex_feasibility.py",
    f"data/results/{RESULT}.json",
    "scripts/c5_paired_variation_probe.py",
    "data/results/c5_paired_variation_probe.json",
    "docs/",
)
#: file names of EN-TEx's allelic releases (the data directory listing, 2026-09-29)
ALLELIC_NAMES = (
    "_default_AS",
    "hetSNVs_",
    "_pooled_AS",
    "AS_allhets",
    "AS_ratios_and_eQTL",
    "AlleleSeq2",
    "imprinted_genes_in_ENTEx_ASE",
    "ASB-predictions-on-GTEx",
)

ELEMENT_ASSAY = "ATAC-seq"
SECONDARY_ELEMENT_ASSAY = "HM-ChIP-seq_H3K27ac"
GENE_ASSAY = "RNA-seq"
KEEP_ASSAYS = {"ccre": (ELEMENT_ASSAY, SECONDARY_ELEMENT_ASSAY), "gene": (GENE_ASSAY,)}

#: the outcome discipline, column by column (module docstring)
COLUMNS_READ = ["chr", "start", "end", "region_id", "experiment_accession", "donor", "tissue", "assay"]
COLUMNS_SUMMED_NEVER_KEPT = ["hap1_count", "hap2_count"]
COLUMNS_REDUCED_UNSIGNED = ["imbalance_significance"]
COLUMNS_NEVER_INDEXED = ["hap1_allele_ratio", "p_betabinom"]
UNSIGNED = {"0": 0, "1": 1}  # EN-TEx's call: 1 imbalanced, 0 not; neither value names a haplotype
CACHE_COLUMNS = ["table", *COLUMNS_READ, "total_depth", "imbalanced_unsigned"]

AUTOSOMES = c5.AUTOSOMES
#: GRCh38 MHC (GRC region); IG_* and TR_* gene types are excluded by GENCODE type
MHC = ("chr6", 28_510_120, 33_480_577)
IMPRINT_EXCLUDE = ("Imprinted", "Predicted", "Tissue Dependent", "Provisional Data", "Conflicting Data")
SV_TYPES = ("DEL", "DUP", "CNV")  # copy-number changing; insertions and inversions keep dosage

#: depth floors, both sides: 0 = every row EN-TEx lists; 47, 85, 194 and 783 are the reads one binomial
#: test needs to tell 0.70, 0.65, 0.60 and 0.55 from 0.50 at alpha 0.05 and 80% power (c5.reads_needed)
FLOORS = (0, 10, 20, 47, 85, 194, 783)
AGREEMENT_RATES = (0.55, 0.6, 0.65, 0.7)
LOCUS_GAP = 1_000_000  # target TSSs closer than this share linked variants: the span C5 counted them over
MIN_ALTERNATIVES = 3
ALTERNATIVE_COUNTS = (1, 3, 5)
DIST_FLOOR = 1_000  # bp: distances below this count as this, so a promoter-proximal element can be matched
DIST_LOG2_TOL = 1.0  # an alternative's TSS lies within half to twice the target's distance from the element
EXPR_LOG2_TOL = 1.0  # log2(GTEx v8 median TPM + 1) in the matching tissue, within 1
DEPTH_LOG2_TOL = 1.0  # log2 RNA total depth in the same donor and tissue, within 1

#: the go/no-go rule, fixed before the first run. GO when, at GO_FLOOR reads on both sides, with the element
#: and gene on one phased block, every exclusion applied and at least MIN_ALTERNATIVES matched alternatives,
#: the rows touch at least GO_LOCI independent loci (`freeze_clusters`) and at least GO_DONORS donors carry
#: GO_DONOR_LOCI loci each (so the donor-specific sensitivity can run), and no EN-TEx allelic table or
#: derivative is read anywhere in the project's code or results outside the C5 probe and this lane.
#: GO_LOCI is the number of independent units a one-sample binomial test needs to tell an agreement rate of
#: 0.60 from 0.50; it counts loci, not the rows that will carry an imbalance on both sides, so it is a
#: necessary condition, not a power guarantee. Loci are single-linkage clusters of target TSSs at LOCUS_GAP,
#: frozen once on the parent universe (every ATAC element-target row formed, before any depth, phase,
#: exclusion or matching filter); a filtered row set counts the parent clusters it touches. The first run
#: recomputed the clusters on each filtered set, which split connected components (566 loci in the parent
#: universe became 781 at a depth of 47) without adding independent evidence (review, relayed 2026-09-29).
GO_FLOOR = c5.reads_needed(0.7)
GO_LOCI = c5.reads_needed(0.6)
GO_DONORS = 3
GO_DONOR_LOCI = 30

#: EN-TEx tissue (lower case, underscores) -> GTEx v8 median-TPM column, for the expression match only
TISSUE_TO_GTEX = {
    "adrenal_gland": "Adrenal Gland",
    "body_of_pancreas": "Pancreas",
    "breast_epithelium": "Breast - Mammary Tissue",
    "coronary_artery": "Artery - Coronary",
    "esophagus_muscularis_mucosa": "Esophagus - Muscularis",
    "esophagus_squamous_epithelium": "Esophagus - Mucosa",
    "gastrocnemius_medialis": "Muscle - Skeletal",
    "gastroesophageal_sphincter": "Esophagus - Gastroesophageal Junction",
    "heart_left_ventricle": "Heart - Left Ventricle",
    "lower_leg_skin": "Skin - Sun Exposed (Lower leg)",
    "omental_fat_pad": "Adipose - Visceral (Omentum)",
    "ovary": "Ovary",
    "peyer's_patch": "Small Intestine - Terminal Ileum",
    "peyers_patch": "Small Intestine - Terminal Ileum",  # the AS tables drop the apostrophe
    "ascending_aorta": "Artery - Aorta",
    "prostate_gland": "Prostate",
    "right_atrium_auricular_region": "Heart - Atrial Appendage",
    "right_lobe_of_liver": "Liver",
    "sigmoid_colon": "Colon - Sigmoid",
    "spleen": "Spleen",
    "stomach": "Stomach",
    "subcutaneous_adipose_tissue": "Adipose - Subcutaneous",
    "suprapubic_skin": "Skin - Not Sun Exposed (Suprapubic)",
    "testis": "Testis",
    "thoracic_aorta": "Artery - Aorta",
    "thyroid_gland": "Thyroid",
    "tibial_artery": "Artery - Tibial",
    "tibial_nerve": "Nerve - Tibial",
    "transverse_colon": "Colon - Transverse",
    "upper_lobe_of_left_lung": "Lung",
    "uterus": "Uterus",
    "vagina": "Vagina",
}

#: requests made by hand while this audit was designed (2026-09-29), before the script existed
DESIGN_TIME_REQUESTS = [
    {"what": "EN-TEx data directory listing", "url": ENTEX, "bytes": 16_353},
    {
        "what": "the first 1 MiB of each AS table by HTTP range, to read the imbalance_significance "
        "vocabulary ('0', '1'); no row was printed and both files were deleted",
        "url": ENTEX + "{cCREs,genes}_default_AS.tsv",
        "bytes": 2 * 1_048_576,
    },
    {"what": "phased_block.tar.gz, to read its layout", "url": PHASED_URL, "bytes": 5_388},
    {
        "what": "the first 20 kB of individual1_SV.vcf, header only",
        "url": SV_URLS["ENC-001"],
        "bytes": 20_001,
    },
    {"what": "ENCODE exclusion list v2, reachability", "url": EXCLUSION_URL, "bytes": 5_867},
    {
        "what": "ENCODE portal: the EN-TEx donors and their external ids (HumanDonor search), two "
        "attempts; both answered 504 Gateway Time-out after 60 s",
        "url": "https://www.encodeproject.org/search/?type=HumanDonor&internal_tags=ENTEx",
        "bytes": 0,
    },
    {
        "what": "ENCODE portal home page, reachability (no answer)",
        "url": "https://www.encodeproject.org/",
        "bytes": 0,
    },
]


# ---------------------------------------------------------------- cost


def _peak_bytes() -> int:
    ru = resource.getrusage(resource.RUSAGE_SELF)
    return int(ru.ru_maxrss if sys.platform == "darwin" else ru.ru_maxrss * 1024)  # bytes on macOS, KiB else


def _cpu() -> float:
    ru = resource.getrusage(resource.RUSAGE_SELF)
    return ru.ru_utime + ru.ru_stime


class Cost:
    """Wall time, CPU, peak memory, bytes downloaded and requests made by this process; every download and
    the run itself are appended to the cache's ledger, so the all-in cost covers every run."""

    def __init__(self, stage: str) -> None:
        self.t0 = time.time()
        self.stage = stage
        self.requests = 0
        self.bytes = 0

    def note(self, url: str, nbytes: int) -> None:
        self.requests += 1
        self.bytes += nbytes
        _ledger({"kind": "download", "url": url, "bytes": nbytes, "at": time.strftime("%Y-%m-%dT%H:%M:%S")})

    def report(self) -> dict[str, Any]:
        return {
            "wall_seconds": round(time.time() - self.t0, 1),
            "cpu_seconds": round(_cpu(), 1),
            "peak_memory_bytes": _peak_bytes(),
            "download_bytes": self.bytes,
            "requests": self.requests,
        }

    def close(self) -> dict[str, Any]:
        r = self.report()
        _ledger({"kind": "run", "stage": self.stage, "at": time.strftime("%Y-%m-%dT%H:%M:%S"), **r})
        return r


def _ledger(entry: dict[str, Any]) -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    with LEDGER.open("a") as fh:
        fh.write(json.dumps(entry) + "\n")


def all_in_cost(this_run: dict[str, Any]) -> dict[str, Any]:
    """Every run and download the ledger holds, plus this run, plus the design-time requests."""
    runs, downloads = [], []
    if LEDGER.exists():
        for line in LEDGER.read_text().splitlines():
            e = json.loads(line)
            (runs if e["kind"] == "run" else downloads).append(e)
    runs = [*runs, {"stage": "this run", **this_run}]
    return {
        "runs": len(runs),
        "wall_seconds": round(sum(r["wall_seconds"] for r in runs), 1),
        "cpu_seconds": round(sum(r["cpu_seconds"] for r in runs), 1),
        "peak_memory_bytes": max(r["peak_memory_bytes"] for r in runs),
        "download_bytes_by_the_script": sum(d["bytes"] for d in downloads),
        "requests_by_the_script": len(downloads),
        "download_bytes_by_hand_while_designing": sum(d["bytes"] for d in DESIGN_TIME_REQUESTS),
        "requests_by_hand_while_designing": len(DESIGN_TIME_REQUESTS) + 1,  # the ENCODE search ran twice
        "downloads": downloads,
        "run_log": runs,
    }


def fetch(url: str, dest: Path, cost: Cost) -> Path:
    """A small public file, downloaded once into the git-ignored cache (the probe's fetch), counted."""
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    c5.fetch(url, dest)
    cost.note(url, dest.stat().st_size)
    return dest


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------- the depth-only cache


def _count(s: str) -> int | float:
    try:
        return int(s)
    except ValueError:
        return float(s)


def depth_only(line: str) -> tuple:
    """One raw AS-table line -> (identifiers..., total depth, unsigned imbalance flag).

    The only function that sees a raw line. hap1_count and hap2_count (fields 4 and 5) are summed here and
    neither is kept; hap1_allele_ratio and p_betabinom (fields 10 and 11) are never indexed;
    imbalance_significance (field 12) is EN-TEx's 0/1 call and names no haplotype."""
    f = line.rstrip("\n").split("\t")
    if len(f) != len(c5.ENTEX_HEADER):
        raise ValueError(f"expected {len(c5.ENTEX_HEADER)} fields, found {len(f)}")
    try:
        total = _count(f[4]) + _count(f[5])
    except ValueError:
        raise ValueError("a haplotype count is not a number (value not shown)") from None
    if f[12] not in UNSIGNED:
        raise ValueError("imbalance_significance outside its 0/1 vocabulary (value not shown)")
    return (f[0], int(f[1]), int(f[2]), f[3], f[6], f[7], f[8], f[9], total, UNSIGNED[f[12]])


def build_depth_cache(cost: Cost) -> dict[str, Any]:
    """Stream both AS tables once through `depth_only` into data/cache/entex/depth_only.tsv.gz; the raw
    bytes are hashed as they pass and never written. The experiment accessions of every assay are kept as
    identifiers for the exposure check."""
    CACHE.mkdir(parents=True, exist_ok=True)
    t0, cpu0 = time.time(), _cpu()
    tmp = DEPTH_CACHE.with_name(DEPTH_CACHE.name + ".part")
    tables: dict[str, Any] = {}
    experiments: dict[str, dict] = {}
    with tmp.open("wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as gz:
        gz.write(("\t".join(CACHE_COLUMNS) + "\n").encode())
        for table, name in AS_TABLES.items():
            src, text = c5._entex_rows(name)  # checks the header against the probe's
            rows = kept = 0
            by_assay: Counter = Counter()
            min_depth: dict[str, float] = {}
            exp = experiments.setdefault(table, {})
            buf: list[str] = []
            for line in text:
                r = depth_only(line)
                rows += 1
                by_assay[r[7]] += 1
                exp.setdefault(r[7], {}).setdefault(r[5], {}).setdefault(r[6], set()).add(r[4])
                if r[7] not in KEEP_ASSAYS[table]:
                    continue
                kept += 1
                min_depth[r[7]] = min(min_depth.get(r[7], r[8]), r[8])
                buf.append(table + "\t" + "\t".join(map(str, r)) + "\n")
                if len(buf) >= 100_000:
                    gz.write("".join(buf).encode())
                    buf.clear()
            gz.write("".join(buf).encode())
            text.close()
            cost.note(src.url, src.bytes)
            tables[table] = src.record(
                rows=rows,
                rows_kept=kept,
                rows_by_assay=dict(by_assay.most_common()),
                min_total_depth_by_assay=min_depth,
                streamed_not_stored=True,
            )
    tmp.rename(DEPTH_CACHE)
    side = {
        "built": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "tables": tables,
        "cache": {
            "path": str(DEPTH_CACHE),
            "sha256": sha256_file(DEPTH_CACHE),
            "bytes": DEPTH_CACHE.stat().st_size,
            "columns": CACHE_COLUMNS,
            "assays_kept": KEEP_ASSAYS,
        },
        "experiments": {
            t: {
                a: {d: {ti: sorted(v) for ti, v in dd.items()} for d, dd in ad.items()}
                for a, ad in td.items()
            }
            for t, td in experiments.items()
        },
        "discipline": {
            "read": COLUMNS_READ,
            "summed_never_kept": COLUMNS_SUMMED_NEVER_KEPT,
            "reduced_unsigned": COLUMNS_REDUCED_UNSIGNED,
            "never_indexed": COLUMNS_NEVER_INDEXED,
        },
        "build_cost": {
            "wall_seconds": round(time.time() - t0, 1),
            "cpu_seconds": round(_cpu() - cpu0, 1),
            "peak_memory_bytes": _peak_bytes(),
            "download_bytes": sum(t["bytes"] for t in tables.values()),
            "requests": len(tables),
        },
    }
    DEPTH_SIDECAR.write_text(json.dumps(side, indent=1) + "\n")
    return side


def read_depth_cache():
    """The cache as three maps: ATAC and H3K27ac per (donor, tissue) -> region -> [depth, flag, experiments],
    and RNA per (donor, tissue) -> Ensembl id (no version) -> [depth, flag, experiments]."""
    chrom: dict[str, dict] = {ELEMENT_ASSAY: defaultdict(dict), SECONDARY_ELEMENT_ASSAY: defaultdict(dict)}
    rna: dict[tuple[str, str], dict[str, list]] = defaultdict(dict)
    with gzip.open(DEPTH_CACHE, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header != CACHE_COLUMNS:
            raise SystemExit(f"depth cache columns changed: {header}")
        for line in fh:
            t, c, s, e, rid, _exp, d, ti, assay, depth, flag = line.rstrip("\n").split("\t")
            dt = (d, ti)
            if t == "gene":
                slot = rna[dt].setdefault(rid.split(".")[0], [0.0, 0, 0])
            else:
                slot = chrom[assay][dt].setdefault((c, int(s), int(e)), [0.0, 0, 0])
            slot[0] += float(depth)
            slot[1] |= int(flag)
            slot[2] += 1
    return chrom[ELEMENT_ASSAY], chrom[SECONDARY_ELEMENT_ASSAY], rna


# ---------------------------------------------------------------- phase, copy number, exclusion list


def read_phase_blocks(path: Path) -> dict[str, dict[str, list[tuple[int, int, str, str]]]]:
    """EN-TEx's phased blocks per individual: chrom -> [(start, end, hap1 origin, hap2 origin)], sorted."""
    out: dict[str, dict[str, list]] = defaultdict(lambda: defaultdict(list))
    with tarfile.open(path, "r:gz") as tf:
        for m in tf.getmembers():
            if not m.isfile():
                continue
            ind = Path(m.name).stem.replace("phased_block_", "")
            for line in tf.extractfile(m).read().decode().splitlines():
                f = line.split("\t")
                if len(f) >= 5:
                    out[ind][f[0]].append((int(f[1]), int(f[2]), f[3], f[4]))
    for ind in out:
        for c in out[ind]:
            out[ind][c].sort()
    return {k: dict(v) for k, v in out.items()}


def block_of(blocks: list[tuple[int, int, str, str]], lo: int, hi: int) -> tuple[int | None, str]:
    """The one phased block holding [lo, hi): (index, 'one'); (None, 'none') when no block holds it whole;
    (None, 'ambiguous') when two overlapping blocks both hold it and the haplotype labels could differ."""
    holding = [k for k, (s, e, _, _) in enumerate(blocks) if s <= lo and hi <= e]
    if len(holding) == 1:
        return holding[0], "one"
    return None, "none" if not holding else "ambiguous"


KNOWN_ORIGIN = ("Paternal", "Maternal")


def orientation(blocks, el_block: int | None, gene_block: int | None) -> str:
    """same_block: one block holds both, so hap1 is one physical haplotype across the interval;
    parental_origin: two blocks, each with a parental origin, so the labels can be aligned by parent;
    not_orientable: otherwise."""
    if el_block is None or gene_block is None:
        return "not_orientable"
    if el_block == gene_block:
        return "same_block"
    if blocks[el_block][2] in KNOWN_ORIGIN and blocks[gene_block][2] in KNOWN_ORIGIN:
        return "parental_origin"
    return "not_orientable"


class Intervals:
    """Sorted intervals per chromosome with a prefix maximum of ends: any-overlap in O(log n)."""

    def __init__(self, by_chrom: dict[str, list[tuple[int, int]]]) -> None:
        self.starts: dict[str, list[int]] = {}
        self.maxend: dict[str, list[int]] = {}
        self.n = 0
        for c, iv in by_chrom.items():
            iv = sorted(iv)
            self.starts[c] = [s for s, _ in iv]
            m, acc = -1, []
            for _, e in iv:
                m = max(m, e)
                acc.append(m)
            self.maxend[c] = acc
            self.n += len(iv)

    def hits(self, chrom: str, lo: int, hi: int) -> bool:
        st = self.starts.get(chrom)
        if not st:
            return False
        k = bisect.bisect_left(st, hi)
        return k > 0 and self.maxend[chrom][k - 1] > lo

    def any_of(self, chrom: str, spans) -> bool:
        return any(self.hits(chrom, s, e) for s, e in spans)


def read_sv(path: Path) -> tuple[Intervals, Counter]:
    """Copy-number-changing SVs (DEL, DUP, CNV) with a non-reference genotype, as intervals."""
    by: dict[str, list[tuple[int, int]]] = defaultdict(list)
    types: Counter = Counter()
    with path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t", 10)
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            t = info.get("SVTYPE", "")
            types[t] += 1
            if not t.startswith(SV_TYPES):
                continue
            gt = f[9].split(":", 1)[0] if len(f) > 9 else ""
            if not any(a not in ("0", ".") for a in gt.replace("|", "/").split("/")):
                continue
            s = int(f[1]) - 1
            e = int(info["END"]) if info.get("END", "").isdigit() else s + len(f[3])
            by[f[0]].append((s, max(e, s + 1)))
    return Intervals(by), types


def read_bed_gz(path: Path) -> Intervals:
    by: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 3 and not line.startswith(("#", "track")):
                by[f[0]].append((int(f[1]), int(f[2])))
    return Intervals(by)


def gtex_log_tpm(path: Path, columns: set[str]) -> dict[str, dict[str, float]]:
    """GTEx v8 median TPM as log2(TPM + 1) per stable gene id, for the named tissue columns only."""
    out: dict[str, dict[str, float]] = {}
    with gzip.open(path, "rt") as fh:
        fh.readline()
        fh.readline()
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header) if h in columns}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[0].split(".")[0]] = {h: math.log2(float(f[i]) + 1) for h, i in idx.items()}
    return out


def norm_tissue(t: str) -> str:
    return t.strip().lower().replace(" ", "_")


def norm_gtex(t: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", t.lower()).strip("_")


# ---------------------------------------------------------------- rows


def element_depth(regions_by_dt, elements) -> tuple[dict, dict[str, Any]]:
    """Per (donor, tissue): archive element -> [summed depth over overlapping EN-TEx regions, any-flag].
    EN-TEx names registry-V2 regions; the archive holds V3 elements; they are joined by overlap."""
    overlap: dict[tuple[str, int, int], list[int]] = {}
    out: dict[tuple[str, str], dict[tuple[str, int], list]] = {}
    for dt, regions in regions_by_dt.items():
        acc: dict[tuple[str, int], list] = {}
        for key, (depth, flag, _n) in regions.items():
            c = key[0]
            if c not in elements:
                continue
            if key not in overlap:
                st, en = elements[c]["start"], elements[c]["end"]
                hi = int(np.searchsorted(st, key[2], side="left"))
                overlap[key] = [i for i in range(max(0, hi - 50), hi) if en[i] > key[1]]
            for i in overlap[key]:
                slot = acc.setdefault((c, i), [0.0, 0])
                slot[0] += depth
                slot[1] |= flag
        out[dt] = acc
    info = {
        "distinct_regions": len(overlap),
        "regions_overlapping_an_element": sum(1 for v in overlap.values() if v),
        "regions_overlapping_two_or_more_elements": sum(1 for v in overlap.values() if len(v) > 1),
    }
    return out, info


def build_rows(el_by_dt, rna, elements, genes) -> tuple[list[dict[str, Any]], Counter]:
    """Element-gene rows: an element with element-side depth and its nominated target with RNA depth, in the
    same donor and tissue. The target is the archive's frozen `predicted_coding` gene."""
    rows: list[dict[str, Any]] = []
    skipped: Counter = Counter()
    for dt in sorted(set(el_by_dt) & set(rna)):
        rd = rna[dt]
        for (c, i), (a, af) in el_by_dt[dt].items():
            g = elements[c]["target"][i]
            if not g:
                continue
            if c not in AUTOSOMES:
                skipped["target_on_chrX_or_chrY"] += 1
                continue
            gene = genes[c].get(g)
            if gene is None:
                skipped["target_not_in_gencode_v50"] += 1
                continue
            r = rd.get(gene["id"])
            if r is None:
                skipped["target_has_no_rna_row_in_this_sample"] += 1
                continue
            mid = (int(elements[c]["start"][i]) + int(elements[c]["end"][i])) // 2
            rows.append(
                {
                    "donor": dt[0],
                    "tissue": dt[1],
                    "chrom": c,
                    "element": i,
                    "gene": g,
                    "atac": a,
                    "atac_flag": af,
                    "rna": r[0],
                    "rna_flag": r[1],
                    "mid": mid,
                    "distance": abs(gene["tss"] - mid),
                    "strength": elements[c]["strength"][i],
                }
            )
    return rows, skipped


def gene_span(gene: dict[str, Any]) -> tuple[int, int]:
    ex = gene["exons"]
    return (ex[0][0], max(e for _, e in ex)) if ex else (gene["tss"], gene["tss"] + 1)


def freeze_clusters(tss_by_chrom: dict[str, set[int]], gap: int = LOCUS_GAP) -> dict[tuple[str, int], int]:
    """Single-linkage clusters of target TSSs, computed once on the parent universe: a new cluster starts
    where the next TSS is more than `gap` away. Returns (chrom, TSS) -> cluster id. Filtered row sets are
    counted by the clusters they touch, so a filter can remove a locus but never split one into two."""
    out: dict[tuple[str, int], int] = {}
    k = -1
    for c in sorted(tss_by_chrom):
        last = None
        for p in sorted(tss_by_chrom[c]):
            if last is None or p - last > gap:
                k += 1
            out[(c, p)] = k
            last = p
    return out


def target_tss(rows: list[dict[str, Any]], genes) -> dict[str, set[int]]:
    tss: dict[str, set[int]] = defaultdict(set)
    for r in rows:
        tss[r["chrom"]].add(genes[r["chrom"]][r["gene"]]["tss"])
    return tss


def summarize(rows: list[dict[str, Any]], genes, clusters: dict[tuple[str, int], int]) -> dict[str, Any]:
    """Rows, element-gene pairs, target genes, parent clusters touched (descriptive), donors and
    donor-tissues."""
    touched: set[int] = set()
    per_donor: dict[str, set[int]] = defaultdict(set)
    for r in rows:
        k = clusters[(r["chrom"], genes[r["chrom"]][r["gene"]]["tss"])]
        touched.add(k)
        per_donor[r["donor"]].add(k)
    return {
        "rows": len(rows),
        "element_gene_pairs": len({(r["chrom"], r["element"], r["gene"]) for r in rows}),
        "target_gene_loci": len({(r["chrom"], r["gene"]) for r in rows}),
        "parent_clusters_touched": len(touched),
        "donors": len({r["donor"] for r in rows}),
        "donor_tissues": len({(r["donor"], r["tissue"]) for r in rows}),
        "parent_clusters_touched_by_donor": {d: len(v) for d, v in sorted(per_donor.items())},
    }


def count_alternatives(row, cand, floor: float, el_block: int | None) -> int:
    """Alternative genes for one row: protein-coding, not the target, not excluded, RNA depth at the floor
    in the same donor and tissue, on the element's phased block, and matched on distance, expression and
    depth within the tolerances fixed before the first run."""
    if cand is None or el_block is None or row.get("expr") is None:
        return 0
    d0 = math.log2(max(row["distance"], DIST_FLOOR))
    dist = np.log2(np.maximum(np.abs(cand["tss"] - row["mid"]), DIST_FLOOR))
    m = (
        (cand["depth"] >= floor)
        & (cand["block"] == el_block)
        & ~cand["excluded"]
        & (cand["gene"] != row["gene"])
        & (np.abs(dist - d0) <= DIST_LOG2_TOL)
        & (np.abs(cand["expr"] - row["expr"]) <= EXPR_LOG2_TOL)
        & (np.abs(np.log2(np.maximum(cand["depth"], 1)) - math.log2(max(row["rna"], 1))) <= DEPTH_LOG2_TOL)
    )
    return int(m.sum())


def balanced_marginals(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """EN-TEx's unsigned call, per side only. The joint is not formed (module docstring)."""
    return {
        "rows": len(rows),
        "element_side_balanced": sum(1 for r in rows if not r["atac_flag"]),
        "gene_side_balanced": sum(1 for r in rows if not r["rna_flag"]),
    }


def go_decision(stage: dict[str, Any], exposure_blocks: bool) -> dict[str, Any]:
    by_donor = stage["parent_clusters_touched_by_donor"]
    donors_ok = sorted(d for d, n in by_donor.items() if n >= GO_DONOR_LOCI)
    checks = {
        f"independent loci (parent clusters touched) >= {GO_LOCI}": stage["parent_clusters_touched"]
        >= GO_LOCI,
        f"donors with >= {GO_DONOR_LOCI} loci >= {GO_DONORS}": len(donors_ok) >= GO_DONORS,
        "no EN-TEx allelic table or derivative read outside the C5 probe and this lane": not exposure_blocks,
    }
    return {
        "go": all(checks.values()),
        "checks": checks,
        "parent_clusters_touched": stage["parent_clusters_touched"],
        "donors_meeting_the_per_donor_floor": donors_ok,
    }


def stages(
    rows, genes, clusters, floor, blocks, el_block_of, gene_block_of, excluded_row, cand_of, with_flags
):
    """The row funnel at one depth floor: depth, then phase, then exclusions, then matched alternatives."""
    s0 = [r for r in rows if r["atac"] >= floor and r["rna"] >= floor]
    orient: Counter = Counter()
    s1, s1p = [], []
    for r in s0:
        b = blocks.get(DONOR_INDIVIDUAL[r["donor"]], {}).get(r["chrom"], [])
        o = orientation(b, el_block_of(r), gene_block_of(r))
        orient[o] += 1
        if o == "same_block":
            s1.append(r)
        if o in ("same_block", "parental_origin"):
            s1p.append(r)
    reasons: Counter = Counter()
    s2 = []
    for r in s1:
        why = excluded_row(r)
        for w in why:
            reasons[w] += 1
        if not why:
            s2.append(r)
    alt_counts = [count_alternatives(r, cand_of(r), floor, el_block_of(r)) for r in s2]
    s3 = [r for r, k in zip(s2, alt_counts, strict=True) if k >= MIN_ALTERNATIVES]
    out: dict[str, Any] = {
        "floor": floor,
        "depth": summarize(s0, genes, clusters),
        "phase_orientation_of_depth_rows": dict(orient),
        "phase_same_block": summarize(s1, genes, clusters),
        "phase_same_block_or_parental_origin": summarize(s1p, genes, clusters),
        "exclusion_reasons_among_same_block_rows": dict(reasons.most_common()),
        "after_exclusions": summarize(s2, genes, clusters),
        "rows_without_a_gtex_expression_value": sum(1 for r in s2 if r.get("expr") is None),
        "rows_by_matched_alternatives": {
            f">={k}": sum(1 for n in alt_counts if n >= k) for k in ALTERNATIVE_COUNTS
        },
        "matched_alternatives_per_row": c5.quartiles(alt_counts),
        f"with_{MIN_ALTERNATIVES}_matched_alternatives": summarize(s3, genes, clusters),
    }
    if with_flags:
        out["balanced_marginals"] = {
            "after_exclusions": balanced_marginals(s2),
            f"with_{MIN_ALTERNATIVES}_matched_alternatives": balanced_marginals(s3),
            "joint": "not formed: rows imbalanced on both sides, set against the product of these "
            "marginals, would be an unsigned coordination statistic and would pre-empt the endpoint",
        }
    return out


# ---------------------------------------------------------------- exposure


def _git(*args: str) -> list[str]:
    return subprocess.run(["git", *args], capture_output=True, text=True, cwd=ROOT).stdout.splitlines()


def exposure(side: dict[str, Any]) -> dict[str, Any]:
    """Where EN-TEx, or anything derived from it, enters the project: tracked files that name it or one of
    its allelic releases, every ENCODE experiment accession in a tracked file checked against EN-TEx's own
    experiment list, and a saved copy of AlphaGenome's track metadata checked for EN-TEx's tissue names."""

    def outside(paths: list[str]) -> list[str]:
        return [p for p in paths if not p.startswith(ALLOWED_TO_NAME)]

    named = _git("grep", "-l", "-i", "-e", "entex", "-e", "en-tex")
    grep_args = [x for n in ALLELIC_NAMES for x in ("-e", n)]
    allelic = _git("grep", "-l", "-F", *grep_args)
    entex_exps: dict[str, dict[str, str]] = {}
    for t, by_assay in side["experiments"].items():
        for assay, by_donor in by_assay.items():
            for donor, by_tissue in by_donor.items():
                for tissue, accs in by_tissue.items():
                    for a in accs:
                        entex_exps[a] = {"table": t, "assay": assay, "donor": donor, "tissue": tissue}
    hits: dict[str, set[str]] = defaultdict(set)
    for line in _git("grep", "-o", "-E", "ENCSR[0-9]{3}[A-Z]{3}"):
        path, _, acc = line.rpartition(":")
        if acc in entex_exps and not path.startswith(ALLOWED_TO_NAME):
            hits[acc].add(path)
    shared = {a: {**entex_exps[a], "files": sorted(p)} for a, p in sorted(hits.items())}
    ag: dict[str, Any] = {"available": AG_METADATA.exists()}
    if AG_METADATA.exists():
        tissues = {
            norm_tissue(t)
            for d in side["experiments"].values()
            for a in d.values()
            for dd in a.values()
            for t in dd
        }
        with AG_METADATA.open() as fh:
            meta = list(csv.DictReader(fh))
        enc = [m for m in meta if m.get("data_source") == "encode" and m.get("biosample_type") == "tissue"]
        match = [m for m in enc if norm_tissue(m.get("biosample_name", "")) in tissues]
        gtex_cols = {norm_gtex(TISSUE_TO_GTEX[t]) for t in tissues if t in TISSUE_TO_GTEX}
        gtex = [m for m in meta if m.get("data_source") == "gtex"]
        ag.update(
            {
                "path": str(AG_METADATA),
                "sha256": sha256_file(AG_METADATA),
                "tracks": len(meta),
                "encode_tissue_tracks": len(enc),
                "encode_tissue_tracks_named_as_an_entex_tissue": dict(
                    Counter(m.get("output", "") for m in match).most_common()
                ),
                "entex_tissue_names_among_them": sorted({norm_tissue(m["biosample_name"]) for m in match}),
                "gtex_rna_tracks": len(gtex),
                "gtex_rna_tracks_for_an_entex_tissue": sum(
                    1 for m in gtex if norm_gtex(m.get("gtex_tissue") or "") in gtex_cols
                ),
                "what_it_cannot_say": "the metadata names a biosample, an ontology term and a data source "
                "per track, never an ENCODE experiment or donor accession, so whether a track was built "
                "from an EN-TEx donor cannot be read from it",
            }
        )
    return {
        "tracked_files_naming_entex": named,
        "tracked_files_naming_an_entex_allelic_release": allelic,
        "of_which_outside_the_c5_probe_this_lane_and_docs": outside(allelic),
        "entex_experiment_accessions_known": len(entex_exps),
        "entex_experiments_named_in_tracked_files": shared,
        "alphagenome_track_metadata": ag,
    }


# ---------------------------------------------------------------- the outcome-exposure record


def outcome_exposure_record(side: dict[str, Any]) -> dict[str, Any]:
    """Every file this lane read, how, and which of its columns; frozen with the checkpoint."""
    t = side["tables"]
    as_table = {
        "columns_read": COLUMNS_READ,
        "columns_summed_on_the_parsing_line_never_kept": COLUMNS_SUMMED_NEVER_KEPT,
        "columns_reduced_to_an_unsigned_flag": COLUMNS_REDUCED_UNSIGNED,
        "columns_never_indexed": COLUMNS_NEVER_INDEXED,
        "how": "streamed once through depth_only on 2026-09-29, hashed as it passed, never written to disk; "
        "the C5 probe streamed the same bytes (same sha256) on 2026-09-28 for identifier columns only",
    }
    return {
        "signed_allelic_outcome_opened": False,
        "files": [
            {
                "file": t["ccre"]["url"],
                "sha256": t["ccre"]["sha256"],
                "bytes": t["ccre"]["bytes"],
                **as_table,
            },
            {
                "file": t["gene"]["url"],
                "sha256": t["gene"]["sha256"],
                "bytes": t["gene"]["bytes"],
                **as_table,
            },
            {
                "file": "the first 1 MiB of each AS table, by HTTP range, while designing",
                "how": "a script printed the header and the distinct values of assay, tissue and "
                "imbalance_significance ('0', '1'); no row, count, ratio or p-value was printed; both files "
                "were deleted",
            },
            {
                "file": PHASED_URL,
                "columns_read": ["chrom", "start", "end", "hap1 parental origin", "hap2 parental origin"],
                "outcome": "none: phase blocks",
            },
            *[
                {
                    "file": u,
                    "columns_read": ["CHROM", "POS", "REF (length)", "INFO SVTYPE and END", "FORMAT GT"],
                    "outcome": "none: the donor's structural-variant genotypes, used to exclude copy-number "
                    "changes; while designing, the first 20 kB of individual1_SV.vcf was read for its header "
                    "and three records (GT:DR:DV) were displayed",
                }
                for u in SV_URLS.values()
            ],
            {"file": EXCLUSION_URL, "columns_read": ["chrom", "start", "end"], "outcome": "none"},
            {"file": c5.GENEIMPRINT_URL, "columns_read": ["gene symbol", "status"], "outcome": "none"},
            {
                "file": c5.GTEX_TPM_URL,
                "columns_read": ["Name", "the EN-TEx tissues' median TPM"],
                "outcome": "none",
            },
            {
                "file": str(AG_METADATA),
                "columns_read": ["biosample_name", "biosample_type", "data_source", "output", "gtex_tissue"],
                "outcome": "none: track descriptions",
            },
            {
                "file": str(c5.ARCHIVE),
                "columns_read": ["start", "end", "id", "predicted_coding.gene", "predicted_coding.strength"],
                "outcome": "none: the frozen nominations",
            },
            {
                "file": "data/reference/gencode_v50_chr*.gff3.gz",
                "columns_read": ["genes and exons"],
                "outcome": "none",
            },
        ],
        "not_opened": [
            "every other file in the EN-TEx data directory, among them the per-SNV allelic tables "
            "(hetSNVs_*_AS.tsv), imprinted_genes_in_ENTEx_ASE.tsv, "
            "AS_allhets_alt_allele_ratio_eqtl_intersect.tsv and AS_ratios_and_eQTL_effect.tsv",
        ],
        "unsigned_quantities_reported": "EN-TEx's 0/1 imbalance call, as per-side marginals of the rows "
        "left after exclusions and of the final rows at each depth floor (the after-exclusion marginals were "
        "added after the first run); the joint of the two sides was never formed and no flag was computed "
        "for an alternative gene",
    }


# ---------------------------------------------------------------- the census


def analyse(cost: Cost, side: dict[str, Any]) -> dict[str, Any]:
    atac, h3k27ac, rna = read_depth_cache()
    atac_dt, rna_dt, h_dt = set(atac), set(rna), set(h3k27ac)
    both = sorted(atac_dt & rna_dt)
    donors = sorted({d for d, _ in atac_dt | rna_dt | h_dt})
    exps = side["experiments"]

    def n_exp(table, assay, d, t):
        return len(exps.get(table, {}).get(assay, {}).get(d, {}).get(t, []))

    samples = {
        "donor_tissues_with_atac": len(atac_dt),
        "donor_tissues_with_rna": len(rna_dt),
        "donor_tissues_with_atac_and_rna": len(both),
        "donor_tissues_with_h3k27ac_and_rna": len(h_dt & rna_dt),
        "distinct_tissues_with_atac_and_rna": len({t for _, t in both}),
        "tissues_with_atac_and_rna_in_how_many_donors": {
            str(k): v for k, v in sorted(Counter(Counter(t for _, t in both).values()).items())
        },
        "by_donor": {
            d: {
                "tissues_with_atac": sum(1 for dd, _ in atac_dt if dd == d),
                "tissues_with_rna": sum(1 for dd, _ in rna_dt if dd == d),
                "tissues_with_atac_and_rna": sorted(t for dd, t in both if dd == d),
            }
            for d in donors
        },
        "donor_tissues_with_two_or_more_atac_experiments": sum(
            1 for d, t in both if n_exp("ccre", ELEMENT_ASSAY, d, t) > 1
        ),
        "donor_tissues_with_two_or_more_rna_experiments": sum(
            1 for d, t in both if n_exp("gene", GENE_ASSAY, d, t) > 1
        ),
    }
    print("samples", {k: v for k, v in samples.items() if k != "by_donor"}, flush=True)

    elements = c5.load_elements()
    genes = {c: c5.load_genes(c) for c in c5.CHROMS}
    id_of = {g["id"]: (c, name) for c, gs in genes.items() for name, g in gs.items()}
    unmapped_tissues = sorted({t for _, t in rna_dt if norm_tissue(t) not in TISSUE_TO_GTEX})
    gtex = gtex_log_tpm(GTEX_TPM, set(TISSUE_TO_GTEX.values()))

    blocks = read_phase_blocks(fetch(PHASED_URL, CACHE / "phased_block.tar.gz", cost))
    sv: dict[str, Intervals] = {}
    sv_types: dict[str, dict[str, int]] = {}
    for d, url in SV_URLS.items():
        sv[d], types = read_sv(fetch(url, CACHE / Path(url).name, cost))
        sv_types[d] = dict(types.most_common())
    excl = read_bed_gz(fetch(EXCLUSION_URL, CACHE / "hg38-blacklist.v2.bed.gz", cost))
    imprint_status = {g: s for s, gs in c5.imprinted_genes().items() for g in gs}

    el_block_cache: dict[tuple, int | None] = {}
    gene_block_cache: dict[tuple, int | None] = {}
    phase_status: Counter = Counter()

    def el_block_of(r) -> int | None:
        key = (r["donor"], r["chrom"], r["element"])
        if key not in el_block_cache:
            b = blocks.get(DONOR_INDIVIDUAL[r["donor"]], {}).get(r["chrom"], [])
            el = elements[r["chrom"]]
            el_block_cache[key], st = block_of(
                b, int(el["start"][r["element"]]), int(el["end"][r["element"]])
            )
            phase_status["element_" + st] += 1
        return el_block_cache[key]

    def gene_block(donor, chrom, name) -> int | None:
        key = (donor, chrom, name)
        if key not in gene_block_cache:
            b = blocks.get(DONOR_INDIVIDUAL[donor], {}).get(chrom, [])
            gene_block_cache[key], st = block_of(b, *gene_span(genes[chrom][name]))
            phase_status["gene_" + st] += 1
        return gene_block_cache[key]

    def gene_block_of(r) -> int | None:
        return gene_block(r["donor"], r["chrom"], r["gene"])

    static_cache: dict[tuple, list[str]] = {}

    def static(chrom, name) -> list[str]:
        if (chrom, name) not in static_cache:
            g = genes[chrom][name]
            why = []
            s = imprint_status.get(name)
            if s in IMPRINT_EXCLUDE:
                why.append(f"imprinting: {s}")
            if (g.get("type") or "").startswith(("IG_", "TR_")):
                why.append("immunoglobulin or T-cell receptor gene")
            if chrom == MHC[0] and MHC[1] <= g["tss"] < MHC[2]:
                why.append("MHC (HLA region)")
            if excl.any_of(chrom, g["exons"]):
                why.append("mapping: an exon of the target on the ENCODE exclusion list")
            static_cache[(chrom, name)] = why
        return static_cache[(chrom, name)]

    def cnv(donor, chrom, spans) -> bool:
        return donor in sv and sv[donor].any_of(chrom, spans)

    def excluded_row(r) -> list[str]:
        why = list(static(r["chrom"], r["gene"]))
        el = elements[r["chrom"]]
        espan = (int(el["start"][r["element"]]), int(el["end"][r["element"]]))
        if r["chrom"] == MHC[0] and espan[0] < MHC[2] and espan[1] > MHC[1]:
            why.append("MHC (HLA region)")
        if excl.hits(r["chrom"], *espan):
            why.append("mapping: element on the ENCODE exclusion list")
        if cnv(r["donor"], r["chrom"], [espan]):
            why.append("copy number: a DEL/DUP/CNV call in this donor overlaps the element")
        if cnv(r["donor"], r["chrom"], genes[r["chrom"]][r["gene"]]["exons"]):
            why.append("copy number: a DEL/DUP/CNV call in this donor overlaps an exon of the target")
        return sorted(set(why))

    cand_cache: dict[tuple, dict[str, np.ndarray] | None] = {}

    def cand_of(r):
        key = (r["donor"], r["tissue"], r["chrom"])
        if key not in cand_cache:
            col = TISSUE_TO_GTEX.get(norm_tissue(r["tissue"]))
            names, tss, depth, expr, blk, exc = [], [], [], [], [], []
            for gid, (depth_v, _flag, _n) in rna[(r["donor"], r["tissue"])].items():
                loc = id_of.get(gid)
                if loc is None or loc[0] != r["chrom"]:
                    continue
                g = genes[loc[0]][loc[1]]
                e = gtex.get(gid, {}).get(col) if col else None
                if g.get("type") != "protein_coding" or e is None:
                    continue
                b = gene_block(r["donor"], loc[0], loc[1])
                names.append(loc[1])
                tss.append(g["tss"])
                depth.append(depth_v)
                expr.append(e)
                blk.append(-1 if b is None else b)
                exc.append(bool(static(loc[0], loc[1])) or cnv(r["donor"], loc[0], g["exons"]))
            cand_cache[key] = (
                {
                    "gene": np.array(names, dtype=object),
                    "tss": np.array(tss, dtype=np.int64),
                    "depth": np.array(depth, dtype=float),
                    "expr": np.array(expr, dtype=float),
                    "block": np.array(blk, dtype=np.int64),
                    "excluded": np.array(exc, dtype=bool),
                }
                if names
                else None
            )
        return cand_cache[key]

    result: dict[str, Any] = {"samples_with_atac_and_rna": samples}
    for label, regions, with_flags in (("atac", atac, True), ("h3k27ac_secondary", h3k27ac, False)):
        el_by_dt, ov = element_depth(regions, elements)
        rows, skipped = build_rows(el_by_dt, rna, elements, genes)
        for r in rows:
            col = TISSUE_TO_GTEX.get(norm_tissue(r["tissue"]))
            r["expr"] = gtex.get(genes[r["chrom"]][r["gene"]]["id"], {}).get(col) if col else None
        print(label, "rows", len(rows), dict(skipped), flush=True)
        clusters = freeze_clusters(target_tss(rows, genes))  # the parent universe: every row formed
        by_floor = [
            stages(
                rows,
                genes,
                clusters,
                f,
                blocks,
                el_block_of,
                gene_block_of,
                excluded_row,
                cand_of,
                with_flags,
            )
            for f in FLOORS
        ]
        result[label] = {
            "independent_locus_count": {
                "parent_universe_clusters": len(set(clusters.values())),
                "target_genes_in_the_parent_universe": len({(r["chrom"], r["gene"]) for r in rows}),
                "definition": "single-linkage clusters of target TSSs at 1 Mb over every element-target row "
                "formed (same donor and tissue, any listed depth), frozen before any filter; the counts "
                "under by_floor are the parent clusters a filtered set touches, descriptive only",
            },
            "element_region_overlap": ov,
            "rows_not_formed": dict(skipped),
            "target_types": dict(
                Counter(genes[r["chrom"]][r["gene"]].get("type") for r in rows).most_common()
            ),
            "depth_quartiles_all_rows": {
                "element_side": c5.quartiles([r["atac"] for r in rows]),
                "gene_side_rna": c5.quartiles([r["rna"] for r in rows]),
            },
            "rows_at_or_above_a_depth": {
                str(f): {
                    "element_side": sum(1 for r in rows if r["atac"] >= f),
                    "gene_side_rna": sum(1 for r in rows if r["rna"] >= f),
                    "both": sum(1 for r in rows if r["atac"] >= f and r["rna"] >= f),
                }
                for f in FLOORS
            },
            "by_floor": by_floor,
        }
        key = f"with_{MIN_ALTERNATIVES}_matched_alternatives"
        print(label, [(s["floor"], s[key]["parent_clusters_touched"]) for s in by_floor], flush=True)
    result["phase_lookups_distinct_donor_element_and_donor_gene"] = dict(phase_status)
    result["phase_blocks"] = {
        ind: {
            "blocks": sum(len(v) for v in b.values()),
            "autosomal_blocks": sum(len(v) for c, v in b.items() if c in AUTOSOMES),
            "autosomal_blocks_over_1mb": sum(
                1 for c, v in b.items() if c in AUTOSOMES for s, e, *_ in v if e - s > 1_000_000
            ),
            "autosomal_bases_in_blocks_over_1mb": sum(
                e - s for c, v in b.items() if c in AUTOSOMES for s, e, *_ in v if e - s > 1_000_000
            ),
            "blocks_with_known_parental_origin": sum(
                1 for v in b.values() for x in v if x[2] in KNOWN_ORIGIN
            ),
        }
        for ind, b in sorted(blocks.items())
    }
    result["copy_number_calls"] = {
        "donors_with_published_sv_calls": sorted(SV_URLS),
        "donors_without": sorted(set(DONOR_INDIVIDUAL) - set(SV_URLS)),
        "sv_types_by_donor": sv_types,
        "copy_number_intervals_by_donor": {d: s.n for d, s in sv.items()},
    }
    result["exclusion_list_intervals"] = excl.n
    result["imprinting_statuses_excluded"] = list(IMPRINT_EXCLUDE)
    result["tissues_without_a_gtex_match"] = unmapped_tissues
    return result


def manifest(side: dict[str, Any]) -> dict[str, Any]:
    inputs = [mf.input_entry(c5.ARCHIVE / f"{c}.json", partition=None) for c in c5.CHROMS]
    inputs += [
        mf.input_entry(Path("data/reference") / f"gencode_v50_{c}.gff3.gz", partition=None) for c in c5.CHROMS
    ]
    for t in side["tables"].values():
        inputs.append(
            {
                "path": t["url"],
                "sha256": t["sha256"],
                "bytes": t["bytes"],
                "rows": t["rows"],
                "partition": None,
                "streamed": True,
            }
        )
    inputs += [
        mf.input_entry(DEPTH_CACHE, partition=None, derived_from="the two AS tables, through depth_only"),
        mf.input_entry(CACHE / "phased_block.tar.gz", partition=None, url=PHASED_URL),
        *[mf.input_entry(CACHE / Path(u).name, partition=None, url=u) for u in SV_URLS.values()],
        mf.input_entry(CACHE / "hg38-blacklist.v2.bed.gz", partition=None, url=EXCLUSION_URL),
        mf.input_entry(c5.CACHE / "geneimprint_human.html", partition=None, url=c5.GENEIMPRINT_URL),
        mf.input_entry(GTEX_TPM, partition=None, url=c5.GTEX_TPM_URL),
    ]
    if AG_METADATA.exists():
        inputs.append(mf.input_entry(AG_METADATA, partition=None))
    m = {
        "sources": [
            {
                "accession": "EN-TEx AS catalogue, phased blocks and SV calls (Rozowsky et al. 2023, Cell)",
                "version": "cCREs_default_AS.tsv and genes_default_AS.tsv 2021-12-10; phased_block.tar.gz "
                "2021-12-10; individual1_SV.vcf and individual4_SV.vcf 2022-09-27",
                "url": ENTEX,
            },
            {
                "accession": "ENCODE SCREEN cCREs scored by AlphaGenome deletion (all-element archive)",
                "version": "AlphaGenome as served during the 2026-09 all-element sweep (unpinned); pinned "
                "here by sha256",
                "path": str(c5.ARCHIVE),
            },
            {"accession": "GENCODE human gene annotation", "version": "release 50 (GRCh38.p14)"},
            {
                "accession": "ENCODE exclusion list v2 (Amemiya et al. 2019), hg38",
                "version": "v2",
                "url": EXCLUSION_URL,
            },
            {
                "accession": "geneimprint.com human imprinted genes (symbols only, as a filter)",
                "version": "as cached 2026-09-28",
                "url": c5.GENEIMPRINT_URL,
            },
            {
                "accession": "GTEx v8 gene median TPM (the EN-TEx tissues' columns, for the expression "
                "match)",
                "version": "2017-06-05 v8 RNASeQC 1.1.9",
            },
            {
                "accession": "AlphaGenome track metadata (output_metadata), a peer session's saved copy, "
                "for the exposure check only",
                "version": "saved 2026-09-27; not re-fetched (a model-service call)",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "floors": list(FLOORS),
            "element_assay": ELEMENT_ASSAY,
            "secondary_element_assay": SECONDARY_ELEMENT_ASSAY,
            "gene_assay": GENE_ASSAY,
            "locus_gap_bp": LOCUS_GAP,
            "locus_definition": "single-linkage clusters of target TSSs at locus_gap_bp, frozen on the "
            "parent universe (every ATAC element-target row formed); filtered sets count the clusters they "
            "touch",
            "min_alternatives": MIN_ALTERNATIVES,
            "distance_floor_bp": DIST_FLOOR,
            "distance_log2_tolerance": DIST_LOG2_TOL,
            "expression_log2_tolerance": EXPR_LOG2_TOL,
            "depth_log2_tolerance": DEPTH_LOG2_TOL,
            "go_floor": GO_FLOOR,
            "go_loci": GO_LOCI,
            "go_donors": GO_DONORS,
            "go_donor_loci": GO_DONOR_LOCI,
            "donor_individual": DONOR_INDIVIDUAL,
            "tissue_to_gtex": TISSUE_TO_GTEX,
            "element_overlap": "at least one base (EN-TEx V2 regions to archive V3 elements); depth summed "
            "over the overlapping regions and over experiments of the same donor and tissue",
            "phase": "element and target gene span (first to last exon) inside one published phased block",
        },
        "exclusions": [
            "chrX and chrY (X inactivation; one X in the male donors)",
            f"target genes geneimprint lists as {', '.join(IMPRINT_EXCLUDE)}",
            "immunoglobulin and T-cell receptor genes (GENCODE IG_* and TR_*) and the MHC, "
            "chr6:28,510,120-33,480,577",
            "an element or target exon on the ENCODE exclusion list v2",
            "an element or target exon under a DEL, DUP or CNV call with a non-reference genotype in that "
            "donor (ENC-001 and ENC-004 only: EN-TEx publishes no SV calls for ENC-002 and ENC-003)",
            "rows whose element and target do not lie in one phased block of that donor",
            "every signed allelic quantity: hap1_count and hap2_count only as their sum; hap1_allele_ratio "
            "and p_betabinom never indexed",
        ],
        "partitions": "n/a: a feasibility census; no evaluation partition is scored and no outcome is opened",
    }
    return mf.with_model_dependencies(m, "alphagenome")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cache-only", action="store_true")
    ap.add_argument("--alphagenome-metadata", type=Path, default=None, help="a saved output_metadata CSV")
    args = ap.parse_args()
    cost = Cost("cache" if args.cache_only else "census")
    if DEPTH_CACHE.exists() and DEPTH_SIDECAR.exists():
        side = json.loads(DEPTH_SIDECAR.read_text())
        if sha256_file(DEPTH_CACHE) != side["cache"]["sha256"]:
            raise SystemExit("depth cache does not match its sidecar; delete both and rebuild")
    else:
        side = build_depth_cache(cost)
    print(
        "depth cache",
        {t: {k: v for k, v in r.items() if k != "rows_by_assay"} for t, r in side["tables"].items()},
        flush=True,
    )
    if args.alphagenome_metadata and not AG_METADATA.exists():
        AG_METADATA.write_bytes(args.alphagenome_metadata.read_bytes())
    if args.cache_only:
        cost.close()
        return
    counts = analyse(cost, side)
    exp = exposure(side)
    blocking = bool(exp["of_which_outside_the_c5_probe_this_lane_and_docs"])
    go_stage = next(s for s in counts["atac"]["by_floor"] if s["floor"] == GO_FLOOR)
    decision = go_decision(go_stage[f"with_{MIN_ALTERNATIVES}_matched_alternatives"], blocking)
    payload = {
        "question": "Checkpoint 1 of audit C: can EN-TEx support an allele-coordination test of frozen "
        "nominated targets against matched alternative genes? Depth and counts only; no signed allelic "
        "outcome opened.",
        "donor_names": "the AS tables name donors ENC-001..ENC-004 and the phased-block and SV files "
        "individual 1..4; they are matched by number, which was not checked against the ENCODE portal (it "
        "answered 504 on 2026-09-29)",
        "not_a_continuation": "a new question with its own rationale; the coherence pilot stopped at "
        "197c560 and its phase C does not follow",
        "go_rule": {
            "floor_both_sides": GO_FLOOR,
            "independent_loci": GO_LOCI,
            "donors": GO_DONORS,
            "loci_per_donor": GO_DONOR_LOCI,
            "min_alternatives": MIN_ALTERNATIVES,
            "fixed": "in scripts/entex_feasibility.py before the first run",
            "locus_clusters": "frozen on the parent universe before any filter (review, relayed 2026-09-29): "
            "the first run recomputed single-linkage clusters on each filtered set, which split connected "
            "components (566 in the parent universe became 781 at a depth of 47) without adding independent "
            "evidence",
            "development_runs": "run 1 recomputed clusters per filtered set; run 2 tried a greedy "
            "pairwise-1-Mb count and two looser matchings, withdrawn unread at the reviewer's direction (no "
            "relaxing of matching or thresholds after seeing counts); the committed run is the rule as "
            "fixed before run 1, with the clusters frozen and two EN-TEx tissue spellings mapped to GTEx",
            "what_it_is": "a necessary condition: loci are counted from depth, and the rows that will "
            "carry an imbalance on both sides are a subset not counted here",
        },
        "decision": decision,
        "reading": "a go means the matched-control design fixed before the first run would have the rows, "
        "loci and donors a registration needs; a no-go means that design is infeasible on EN-TEx, not that "
        "EN-TEx lacks usable allelic evidence. A different control or measurement-error design would need "
        "its own rationale and registration; none is started here.",
        "outcome_discipline": {
            "columns_read": COLUMNS_READ,
            "columns_summed_never_kept": COLUMNS_SUMMED_NEVER_KEPT,
            "columns_reduced_to_an_unsigned_flag": COLUMNS_REDUCED_UNSIGNED,
            "columns_never_indexed": COLUMNS_NEVER_INDEXED,
            "how": "depth_only is the only function that sees a raw line; its output is invariant to "
            "swapping hap1 and hap2 and to the ratio and p-value fields (tests); the raw tables were "
            "streamed and hashed, never written; the cache holds only depth_only's output; the unsigned "
            "flag is reported as per-side marginals, never jointly and never for alternative genes",
        },
        "outcome_exposure_record": outcome_exposure_record(side),
        "depth_cache": {k: side[k] for k in ("built", "tables", "cache", "build_cost")},
        **counts,
        "exposure": exp,
        "independence": "unresolved. EN-TEx's four donors are GTEx donors; GTEx v8 eQTLs enter the "
        "attribution layer (eqtl_targets) and GTEx RNA-seq tracks are among AlphaGenome's outputs. Whether "
        "the four donors are in GTEx v8's eQTL cohort, and whether their samples trained AlphaGenome, is "
        "not established here, so a result on EN-TEx is not an independent validation of anything fitted "
        "on GTEx or by AlphaGenome.",
        "power_units_needed": {
            "one_sample_binomial_units_to_tell_rate_from_0_50": {
                str(p): c5.reads_needed(p) for p in AGREEMENT_RATES
            },
            "note": "units are independent loci; matched alternatives add variance (a factor near 1 + 1/k "
            "for k alternatives) and the fraction of rows imbalanced on both sides is not used, so these "
            "are optimistic",
        },
        "design_time_requests": DESIGN_TIME_REQUESTS,
        "alphagenome_requests": 0,
    }
    payload["result_manifest"] = manifest(side)
    this_run = cost.report()
    payload["cost"] = {
        "this_run": this_run,
        "depth_cache_build": side["build_cost"],
        "all_in": all_in_cost(this_run),
    }
    p = save_result(RESULT, payload)
    cost.close()
    print("wrote", p, "go" if decision["go"] else "no-go", decision["checks"], flush=True)


if __name__ == "__main__":
    main()
