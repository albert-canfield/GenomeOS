# SPDX-License-Identifier: AGPL-3.0-or-later
"""Probe the IGVF portal for a second independent CRISPRi benchmark; label-blind, 0 model requests.

    uv run python scripts/indep_igvf_probe.py

lane-indep's first benchmark (the IGVF K562 MHC Perturb-seq, attribution/indep.py) was NOT READABLE,
and it named the IGVF K562 CRISPRi Perturb-seq over many loci (Gersbach, 16 measurement sets,
2026-07-27) as the next 0-request candidate, then without an element-level table. This script records
what the IGVF API (https://api.data.igvf.org) held on 2026-09-28 and makes the two label-blind counts
that decide what can be done now:

1. the only released K562 CRISPRi element-level table outside the MHC, the Engreitz "K562 Random
   DC-TAP-seq" (IGVFDS7288SJVF, file IGVFFI0957PYTA), against the ENCODE EPCrisprBenchmark files the
   project has read. It reads coordinates, gene, type, power and promoter columns only; never the
   effect, p-value or significance columns;
2. the guide library of the 16 new K562 sets (IGVFFI9377CEMM, 350 DHSs), against the stored
   AlphaGenome deletion sweep: elements covered, and element-gene pairs that carry a K562 value.

Nothing is scored and no label is joined to a prediction. Inputs are fetched to the git-ignored
data/cache/indep/ and checked against the portal's md5.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri, indep  # noqa: E402
from genomeos.results import save_result  # noqa: E402

NAME = "indep_igvf_probe"
API = "https://api.data.igvf.org"
LICENCE = "CC BY 4.0 (IGVF data)"

DCTAP = {
    "accession": "IGVFFI0957PYTA",
    "analysis_set": "IGVFDS7288SJVF",
    "lab": "Jesse Engreitz, Stanford (IGVF)",
    "what": "K562 Random DC-TAP-seq, 664 elements drawn at random from accessible peaks in 25 loci",
    "content_type": "differential gene expression quantifications (formerly 'differential element "
    "quantifications' IGVFFI2436YUJO, revoked, and IGVFFI9586LSUT, archived)",
    "bytes": 565_943,
    "md5": "758e116216267969281d1b4c5f0c4b1a",
    "released": "2026-08-08",
    "licence": LICENCE,
    "url": f"{API}/tabular-files/IGVFFI0957PYTA/@@download/IGVFFI0957PYTA.tsv.gz",
}
LIBRARY = {
    "accession": "IGVFFI9377CEMM",
    "construct_library_set": "IGVFDS0835NMEN",
    "lab": "Charles Gersbach, Duke (IGVF)",
    "what": "guide RNA sequences: CRISPRi screen targeting 350 DHS genome-wide with transcriptome-wide "
    "Perturb-seq and growth readout (the library of the 16 K562 sets released 2026-07-27)",
    "bytes": 39_145,
    "md5": "2c1975d5024d8b840343213bd8f42d79",
    "released": "2026-07-27",
    "licence": LICENCE,
    "url": f"{API}/tabular-files/IGVFFI9377CEMM/@@download/IGVFFI9377CEMM.tsv.gz",
}

# The 16 K562 CRISPRi Perturb-seq measurement sets (Gersbach, released 2026-07-27), with the analysis set
# each is `input_for`, and the size of its raw files (fastq + seqspec), read from the API on 2026-09-28.
MANY_LOCI_SETS = {
    "IGVFDS5637QJRP": ("IGVFDS3624DGBH", 35.2),
    "IGVFDS8708IDCC": ("IGVFDS3624DGBH", 28.4),
    "IGVFDS5440YEBR": ("IGVFDS3624DGBH", 32.1),
    "IGVFDS6162UVFJ": ("IGVFDS3624DGBH", 43.1),
    "IGVFDS0567GXFW": ("IGVFDS3624DGBH", 43.6),
    "IGVFDS5791LZOT": ("IGVFDS3624DGBH", 30.7),
    "IGVFDS5632GXBM": ("IGVFDS3624DGBH", 41.4),
    "IGVFDS2463TEJP": ("IGVFDS3624DGBH", 36.5),
    "IGVFDS1341HYXT": ("IGVFDS3624DGBH", 18.1),
    "IGVFDS8094NBFV": ("IGVFDS3624DGBH", 38.4),
    "IGVFDS0146NONB": ("IGVFDS3624DGBH", 48.5),
    "IGVFDS6790CUPS": ("IGVFDS3624DGBH", 45.1),
    "IGVFDS6013FBZQ": ("IGVFDS3481OUHP", 42.6),
    "IGVFDS6260UZSB": ("IGVFDS3624DGBH", 31.3),
    "IGVFDS1474DQOO": ("IGVFDS3624DGBH", 38.3),
    "IGVFDS8491ILFO": ("IGVFDS3624DGBH", 43.8),
}
GRNA_AUX_GB = 93.9  # the 16 gRNA-capture auxiliary sets (fastq), one per measurement set

# Every other released IGVF CRISPRi element-level table (content types 'differential element
# quantifications', 'local/global differential expression per element'), by cell, and why the stored
# sweep cannot score it with 0 requests. The sweep holds K562, HepG2, GM12878 and IMR-90 tracks only.
OTHER_TABLES = [
    {
        "cell": "K562",
        "analysis_set": "IGVFDS7132YVKO",
        "file": "IGVFFI4093WUVB",
        "lab": "Gersbach",
        "why": "the MHC screen: read by lane-indep, NOT READABLE",
    },
    {
        "cell": "K562",
        "analysis_set": "IGVFDS7288SJVF",
        "file": "IGVFFI0957PYTA",
        "lab": "Engreitz",
        "why": "Ray 2025 K562 DC-TAP, already in the ENCODE held-out file (see dctap)",
    },
    {
        "cell": "WTC11",
        "analysis_set": "IGVFDS0523KGZP",
        "file": "IGVFFI9490RJAU (998 kB)",
        "lab": "Gersbach",
        "why": "no WTC11 track in the sweep",
    },
    {
        "cell": "WTC11",
        "analysis_set": "IGVFDS4617SJVM",
        "file": "IGVFFI9246AJEK (15 kB)",
        "lab": "Engreitz",
        "why": "no WTC11 track in the sweep",
    },
    {
        "cell": "WTC11",
        "analysis_set": "IGVFDS3911MOCN",
        "file": "differential element quantifications",
        "lab": "Engreitz",
        "why": "no WTC11 track; WTC11 DC-TAP is also in the ENCODE held-out file",
    },
    {
        "cell": "WTC11",
        "analysis_set": "IGVFDS7340YDHF, IGVFDS6332VCTO, IGVFDS4389OUWU, IGVFDS4003HZAB",
        "file": "local/global differential expression per element",
        "lab": "Hon",
        "why": "no WTC11 track in the sweep (cell line and cardiomyocyte differentiations)",
    },
    {
        "cell": "WTC11-derived",
        "analysis_set": "IGVFDS3959LESA, IGVFDS2406AZOM, IGVFDS1172VSTV",
        "file": "differential element quantifications",
        "lab": "Gersbach, Engreitz",
        "why": "neurons and endothelial cells from WTC11: no track",
    },
    {
        "cell": "HCT116",
        "analysis_set": "12 FlowFISH / TAP-seq sets",
        "file": "differential element quantifications",
        "lab": "Engreitz",
        "why": "no HCT116 track: the 705-request item",
    },
    {
        "cell": "HepG2",
        "analysis_set": "IGVFDS6504OLWV",
        "file": "guide quantifications only",
        "lab": "Sherwood",
        "why": "CRISPRi FACS on a phenotype, no element-gene table",
    },
    {
        "cell": "GM12878, IMR-90",
        "analysis_set": "none",
        "file": "none",
        "lab": "",
        "why": "no released CRISPRi screen on the portal",
    },
]
PENDING_ANALYSIS_SETS = {
    "IGVFDS3624DGBH": "input of 15 of the 16 many-loci K562 sets; exists, HTTP 403 (not released)",
    "IGVFDS3481OUHP": "input of IGVFDS6013FBZQ; exists, HTTP 403 (not released)",
    "IGVFDS1588RHDM": "input of three further K562 sets (IGVFDS1915EGTT, IGVFDS2136PGKM, IGVFDS5321WVPQ, "
    "2026-02-06) on the MHC library IGVFDS2205IYHG; HTTP 403 (not released)",
}


def fetch(src: dict[str, Any]) -> Path:
    path = indep.CACHE / Path(src["url"]).name
    if not path.exists():
        indep.CACHE.mkdir(parents=True, exist_ok=True)
        subprocess.run(["curl", "-sfL", "-o", str(path), src["url"]], check=True)
    digest = hashlib.md5(path.read_bytes()).hexdigest()
    if digest != src["md5"]:
        raise SystemExit(f"{path}: md5 {digest}, the portal says {src['md5']}")
    return path


def encode_spans() -> tuple[dict[str, list[tuple[str, int, int, str]]], dict[str, list[tuple[int, int]]]]:
    """Every ENCODE EPCrisprBenchmark element (training and held-out): per gene symbol, and per chromosome."""
    by_gene: dict[str, list[tuple[str, int, int, str]]] = defaultdict(list)
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for name in (crispri.TRAINING, crispri.HELDOUT):
        with gzip.open(crispri.KNOWLEDGE / name, "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                s, e = int(r["chromStart"]), int(r["chromEnd"])
                by_gene[r["measuredGeneSymbol"]].append((r["chrom"], s, e, r["Dataset"]))
                by_chrom[r["chrom"]].append((s, e))
    return by_gene, by_chrom


def dctap_overlap(path: Path, by_gene: dict[str, list[tuple[str, int, int, str]]]) -> dict[str, Any]:
    """The DC-TAP pairs against the ENCODE files. Reads no effect, p-value or significance column."""
    keep = (
        "type",
        "targeting_chr",
        "targeting_start",
        "targeting_end",
        "gene_symbol",
        "power_at_effect_size_15",
        "promoter_gene",
    )
    by_dataset: Counter[str] = Counter()
    strata: Counter[str] = Counter()
    elements, types = set(), Counter()
    with gzip.open(path, "rt") as fh:
        for full in csv.DictReader(fh, delimiter="\t"):
            r = {k: full[k] for k in keep}
            types[r["type"]] += 1
            if r["type"] != "targeting":
                continue
            c, s, e, g = (
                r["targeting_chr"],
                int(r["targeting_start"]),
                int(r["targeting_end"]),
                r["gene_symbol"],
            )
            elements.add((c, s, e))
            hits = sorted({d for cc, a, b, d in by_gene.get(g, ()) if cc == c and a < e and b > s})
            by_dataset["+".join(hits) or "not in either ENCODE file"] += 1
            try:
                powered = float(r["power_at_effect_size_15"]) >= 0.8
            except ValueError:
                powered = False
            where = "promoter" if r["promoter_gene"] not in ("", "NA", "None") else "distal"
            strata[
                f"{'in ENCODE' if hits else 'not in ENCODE'}, {where}, power15 {'>=' if powered else '<'} 0.8"
            ] += 1
    return {
        "rows_by_type": dict(types),
        "elements": len(elements),
        "pairs_by_encode_dataset": dict(by_dataset),
        "pairs_by_stratum": dict(sorted(strata.items())),
    }


def library_coverage(path: Path, by_chrom: dict[str, list[tuple[int, int]]]) -> dict[str, Any]:
    """The 350-DHS library's targets against the stored sweep: elements covered and K562 pairs in reach."""
    with gzip.open(path, "rt") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    types = Counter(r["type"] for r in rows)
    targets = sorted(
        {
            (r["intended_target_chr"], int(r["intended_target_start"]), int(r["intended_target_end"]))
            for r in rows
            if r["type"] == "targeting" and r["intended_target_chr"]
        }
    )
    table, cache = crispri.DeletionTable(), crispri.ElementCache()
    counts: Counter[str] = Counter()
    pairs, chroms = 0, sorted({c for c, _, _ in targets})
    for chrom in chroms:
        genes = indep.gencode(chrom)
        for c, s, e in (t for t in targets if t[0] == chrom):
            if indep.near_tss(genes["all_tss"], s, e):
                counts["promoter (within 1 kb of a GENCODE v50 TSS)"] += 1
                continue
            counts["distal"] += 1
            if any(a < e and b > s for a, b in by_chrom.get(c, ())):
                counts["distal, overlapping an ENCODE benchmark element"] += 1
            hit = [x for x in table.overlapping(c, s, e) if x["end"] > s and x["start"] < e]
            if not hit:
                continue
            counts["distal, covered by the sweep"] += 1
            k562 = {
                g["gene"]
                for x in hit
                for g in _record(cache, c, x["id"])
                if (g.get("by_cell") or {}).get("K562") is not None
            }
            pairs += len(k562)
    distal = counts["distal"]
    return {
        "guides_by_type": dict(types),
        "targets": len(targets),
        "chromosomes": len(chroms),
        "chromosome_list": chroms,
        "counts": dict(counts),
        "covered_share_of_distal": round(counts["distal, covered by the sweep"] / distal, 3)
        if distal
        else None,
        "element_gene_pairs_with_a_K562_value": pairs,
    }


def _record(cache: crispri.ElementCache, chrom: str, element_id: str) -> list[dict[str, Any]]:
    cache._load(chrom)
    hit = cache._archive.get(element_id)
    if hit is None:
        p = cache.root / chrom / f"{element_id}.json"
        if p.exists():
            hit = json.loads(p.read_text())
    return (hit or {}).get("genes") or []


@mf.depends_on_models("alphagenome")  # the sweep's model, as far as the disk says (R9)
def build_manifest(paths: list[Path], chroms: list[str]) -> dict[str, Any]:
    return {
        "sources": [
            {**DCTAP, "version": f"released {DCTAP['released']}"},
            {**LIBRARY, "version": f"released {LIBRARY['released']}"},
            {
                "accession": "EngreitzLab/CRISPR_comparison EPCrisprBenchmark",
                "version": "main (sha256 pinned)",
            },
            {"accession": "GENCODE", "version": "v50"},
            {"accession": "AlphaGenome all-element deletion sweep", "version": "stored, 2026-09-13..16"},
        ],
        "inputs": [
            mf.input_entry(paths[0], partition="independent (label-blind overlap only)"),
            mf.input_entry(paths[1], partition="independent (label-free guide library)"),
            mf.input_entry(crispri.KNOWLEDGE / crispri.TRAINING, partition="training (overlap check only)"),
            mf.input_entry(crispri.KNOWLEDGE / crispri.HELDOUT, partition="heldout (overlap check only)"),
        ]
        + [
            mf.input_entry(crispri.ELEMENT_CACHE / f"{c}.json.gz", partition=None)
            for c in chroms
            if (crispri.ELEMENT_CACHE / f"{c}.json.gz").exists()
        ],
        "assembly": "GRCh38",
        "exclusions": [
            {"reason": "DC-TAP: effect, p-value and significance columns not read (label-blind)", "pairs": 0},
            {"reason": "library: non-targeting and positive-control guides, not elements", "pairs": 0},
            {
                "reason": "library: targets within 1 kb of a GENCODE v50 TSS, counted apart as promoters",
                "pairs": 0,
            },
        ],
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "promoter_window": indep.PROMOTER_WINDOW,
            "power_floor": 0.8,
            "min_positives": indep.MIN_POSITIVES,
            "min_coverage": indep.MIN_COVERAGE,
        },
        "partitions": {
            "independent": "IGVF files read by no earlier script; no label column read, nothing scored",
            "training": "ENCODE K562 training pairs, read only for overlap",
            "heldout": "ENCODE held-out pairs, read only for overlap",
        },
    }


def main() -> None:
    dctap_path, lib_path = fetch(DCTAP), fetch(LIBRARY)
    by_gene, by_chrom = encode_spans()
    dctap = dctap_overlap(dctap_path, by_gene)
    lib = library_coverage(lib_path, by_chrom)
    raw_gb = round(sum(gb for _, gb in MANY_LOCI_SETS.values()), 1)
    payload = {
        "question": "Has IGVF released a processed element-gene table for a K562/HepG2/WTC11 CRISPRi screen "
        "the project has not read, that the stored sweep can score with 0 requests?",
        "answer": "no. The 16 many-loci K562 sets are raw reads only; their analysis sets exist but are "
        "unreleased. The one other released K562 table is Ray 2025 DC-TAP, already read.",
        "verdict": "NO BENCHMARK NOW: raw single-cell data only; processing is out of this machine's reach",
        "requests": 0,
        "many_loci_k562": {
            "measurement_sets": {
                k: {"analysis_set": a, "raw_gb": gb} for k, (a, gb) in MANY_LOCI_SETS.items()
            },
            "raw_fastq_gb": raw_gb,
            "grna_aux_fastq_gb": GRNA_AUX_GB,
            "released": "2026-07-27",
            "licence": LICENCE,
            "pending_analysis_sets": PENDING_ANALYSIS_SETS,
            "guide_library": LIBRARY,
            "guide_library_coverage": lib,
            "processing_estimate": {
                "download_gb": round(raw_gb + GRNA_AUX_GB, 1),
                "disk_free_gib": 19,
                "tools": "IGVF CRISPR Perturb-seq pipeline (Nextflow): kallisto|bustools or Cell Ranger for "
                "gene and guide counts, guide assignment, then SCEPTRE or PerTurbo "
                "differential expression per element and gene",
                "compute": "16 sets of about 37 GB each; per set a counting step of several "
                "CPU-hours and 32-64 GB RAM; the differential test over 5,676 cis pairs "
                "(or all genes) adds hours",
                "reading": "691 GB of reads against 19 GiB free and a 2 GB cap: not feasible here, "
                "even one set at a time. Wait for IGVFDS3624DGBH to be released, "
                "then pre-register on its table.",
            },
        },
        "dctap": {**DCTAP, "overlap": dctap},
        "other_released_tables": OTHER_TABLES,
        "evidence": "portal metadata (IGVF API, 2026-09-28); label-blind counts only; "
        "no prediction joined to a label",
    }
    manifest = build_manifest([dctap_path, lib_path], lib["chromosome_list"])
    print(save_result(NAME, payload, manifest=manifest))
    print(
        {
            "dctap": dctap["pairs_by_stratum"],
            "library": lib["counts"],
            "pairs": lib["element_gene_pairs_with_a_K562_value"],
        }
    )


if __name__ == "__main__":
    main()
