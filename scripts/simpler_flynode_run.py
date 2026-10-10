# SPDX-License-Identifier: AGPL-3.0-or-later
"""Simpler genomes, the fly's boundary node: fetch, inventory, census and registration, then the run.

    uv run python scripts/simpler_flynode_run.py fetch       # the boundary file into the git-ignored store
    uv run python scripts/simpler_flynode_run.py inventory   # data/results/simpler_flynode_inventory.json
    uv run python scripts/simpler_flynode_run.py register    # census + registration, no RNA value read
    uv run python scripts/simpler_flynode_run.py score       # refuses unless the registration is committed

The population, the endpoint and the RNA are simpler_fly's (registration fd81e37), read from its store
and checked against its committed registration. No model request is sent and no key is read.
docs/SIMPLER-GENOMES.md section 9 holds the reasoning.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import simpler  # noqa: E402
from genomeos.attribution import simpler_fly as sf  # noqa: E402
from genomeos.attribution import simpler_flynode as fn  # noqa: E402
from genomeos.results import save_result  # noqa: E402

INVENTORY = "simpler_flynode_inventory"
REGISTRATION = "simpler_flynode_registration"
RESULT = "simpler_flynode_result"
REGISTRATION_FILE = "data/results/simpler_flynode_registration.json"
FLY_REGISTRATION_FILE = "data/results/simpler_fly_registration.json"
ENTRY = "scripts/simpler_flynode_run.py"
OWN_CODE = (
    "genomeos/attribution/simpler.py",
    "genomeos/attribution/simpler_fly.py",
    "genomeos/attribution/simpler_flynode.py",
    "scripts/simpler_flynode_run.py",
    "tests/test_simpler_flynode.py",
)
BUDGET_BYTES = 200_000_000

#: Bytes moved while choosing the set, outside the fetch: the GEO copy of Stadler 2017's tables
#: (GSE100370, dm3 coordinates only), a scouting copy of the eLife file, and Cavalheiro 2023's
#: supplementary archive from Europe PMC (10-kb boundaries). Page reads and HEAD requests not counted.
SCOUTING_BYTES = 261_022 + 386_317 + 10_991_845

CONSIDERED = [
    {
        "source": "Stadler, Haines, Eisen 2017, eLife 6:e29550, supplementary file 1",
        "what": "Hi-C boundaries of nc14 embryos at 500-bp resolution: manual, computational, representative",
        "assembly": "dm6 as given (dm3 alongside)",
        "licence": "CC BY 4.0",
        "used": True,
    },
    {
        "source": "the same tables as deposited at GEO GSE100370 (Supplemental_tables_fullSub.xlsx)",
        "what": "an earlier version of the same lists",
        "assembly": "dm3 coordinates only",
        "licence": "GEO record; no licence stated by the submitters",
        "used": False,
        "why_not": "dm3 would need lifting; the published eLife file gives dm6",
    },
    {
        "source": "Cavalheiro et al. 2023, Sci Adv 9:eade1085, table S1",
        "what": "boundaries of nc14 embryos scored for insulator depletion, HiCExplorer at 10 kb",
        "assembly": "dm6",
        "licence": "CC BY 4.0",
        "used": False,
        "why_not": (
            "10-kb bins are coarser than many fly intergenic gaps, and the table lists boundaries as "
            "rows of a differential analysis (878) rather than the 772 the text reports"
        ),
    },
    {
        "source": "Ogiyama et al. 2018, Mol Cell 71:73, GEO GSE103625",
        "what": "Hi-C contact matrices over embryo stages",
        "assembly": "dm6",
        "licence": "GEO record; no licence stated by the submitters",
        "used": False,
        "why_not": "no published boundary list: calling one would choose a caller and a threshold here",
    },
    {
        "source": "Hug et al. 2017, Cell 169:216",
        "what": "Hi-C over zygotic genome activation",
        "assembly": "dm6",
        "licence": "not open access",
        "used": False,
        "why_not": "no openly licensed boundary list found",
    },
    {
        "source": "Ramirez et al. 2018, Nat Commun 9:189",
        "what": "high-resolution boundaries",
        "assembly": "not checked",
        "licence": "CC BY 4.0",
        "used": False,
        "why_not": "Kc167 cells, not embryos",
    },
]


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch() -> None:
    fn.STORE.mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(".").free
    req = urllib.request.Request(fn.BOUNDARY_URL, headers={"User-Agent": "GenomeOS/1.0"})
    with urllib.request.urlopen(req, timeout=120) as r:
        body = r.read()
    if len(body) > BUDGET_BYTES:
        sys.exit("refused: over the budget")
    fn.BOUNDARY_FILE.write_bytes(body)
    rec = {
        "fetched_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
        "disk_free_bytes_before": free,
        "files": {
            fn.BOUNDARY_FILE.name: {
                "url": fn.BOUNDARY_URL,
                "what": fn.BOUNDARY_SOURCE,
                "bytes": len(body),
                "sha256": _sha(fn.BOUNDARY_FILE),
                "licence": fn.BOUNDARY_LICENCE,
                "kept": True,
            }
        },
        "download_bytes": len(body),
    }
    fn.FETCH_RECORD.write_text(json.dumps(rec, indent=1) + "\n")
    print(json.dumps(rec, indent=1))


def _manifest(
    inputs: list[tuple[str, str]], parameters: dict, exclusions: list[str], partitions: Any
) -> dict:
    return {
        "sources": [
            {
                "accession": fn.BOUNDARY_SOURCE,
                "version": "eLife supplementary file 1, v2, as served on 2026-10-10",
            },
            {
                "accession": sf.SOURCE,
                "version": "as fetched by lane-fly on 2026-10-10 (registration fd81e37)",
            },
            {"accession": sf.GENE_MODELS, "version": "Ensembl release 104"},
        ],
        "inputs": [mf.input_entry(path, partition=part) for path, part in inputs],
        "assembly": (
            "dm6 (BDGP6.32): the boundaries as the file gives them, the peaks, the RNA and the gene models"
        ),
        "coordinates": {
            "base": 0,
            "interval": "half-open",
            "note": (
                "the peaks are BED (0-based, half-open), the GFF3 1-based closed, and a boundary is the "
                "500-bp bin start + 250; positions move by at most 1 bp between the conventions"
            ),
        },
        "parameters": parameters,
        "exclusions": exclusions,
        "partitions": partitions,
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


PARAMETERS = {
    "boundary_lists": fn.SHEETS,
    "primary_list": fn.PRIMARY_SET,
    "bin_bp": fn.BIN_BP,
    "rotation_draws": fn.ROTATION_DRAWS,
    "tissues": list(sf.TISSUES),
    "detection_cpm": sf.DETECTION_CPM,
    "distance_bins_bp": list(sf.DISTANCE_BINS),
    "control_draws": sf.CONTROL_DRAWS,
    "bootstrap": sf.BOOTSTRAP,
    "cluster_bp": sf.CLUSTER_BP,
    "seed": sf.SEED,
    "min_discordant": fn.MIN_DISCORDANT,
    "alpha": fn.ALPHA,
}
EXCLUSIONS = [
    "simpler_fly's population exclusions, unchanged (registration fd81e37)",
    "primary population: distal elements whose midpoint lies in a protein-coding gene body",
    "boundary rows on chromosomes other than 2L, 2R, 3L, 3R, 4, X, Y",
    "sheets of the boundary file other than Tables S4, S5 and S7",
]


def _record() -> dict:
    return json.loads(fn.FETCH_RECORD.read_text())


def _boundaries() -> tuple[dict[str, dict[str, list[int]]], dict[str, Any]]:
    sets, summaries = {}, {}
    for name in fn.SHEETS:
        sets[name], summaries[name] = fn.read_boundaries(fn.BOUNDARY_FILE, name)
    return sets, summaries


def inventory() -> None:
    mf.trace_begin()
    rec = _record()
    f = rec["files"][fn.BOUNDARY_FILE.name]
    _, summaries = _boundaries()
    out = {
        "question": (
            "one small public boundary set for dm6 embryos: what was fetched, its size, hash, licence and why"
        ),
        "fetched": [{**f, "name": fn.BOUNDARY_FILE.name}],
        "fetched_bytes": rec["download_bytes"],
        "scouting_bytes": SCOUTING_BYTES,
        "lane_download_bytes": rec["download_bytes"] + SCOUTING_BYTES,
        "budget_bytes": BUDGET_BYTES,
        "disk_free_bytes_before": rec["disk_free_bytes_before"],
        "lists": summaries,
        "why_this_set": fn.REGISTRATION["boundary_set"]["why_this_set"],
        "stage": fn.REGISTRATION["boundary_set"]["stage"],
        "considered": CONSIDERED,
        "money": "none; no model request",
    }
    out["result_manifest"] = _manifest(
        [(str(fn.FETCH_RECORD), "the fetch record"), (str(fn.BOUNDARY_FILE), "Tables S4, S5 and S7 only")],
        {"download_budget_bytes": BUDGET_BYTES},
        ["sheets other than Tables S4, S5 and S7 are not read"],
        "n/a: an inventory has no evaluation partition",
    )
    print(save_result(INVENTORY, out))


def _rows() -> tuple[dict, list[dict[str, Any]], dict[str, Any], dict[str, str]]:
    """simpler_fly's population, exactly as scripts/simpler_fly_run.py builds it."""
    chroms, tss, name_to_id = sf.read_annotation(sf.GFF3)
    calls, funnel = sf.tissue_calls(sf.read_markers(sf.SUPP), sf.read_cluster_labels(sf.SUPP))
    els, counts = sf.elements(calls, sf.read_peaks(sf.PEAKS), chroms, tss)
    return chroms, simpler.attribute(els, chroms), {"marker_calls": funnel, "classes": counts}, name_to_id


def _build() -> tuple[dict, list[dict], dict[str, list[dict]], dict, dict, dict, dict[str, str]]:
    fly = json.loads(Path(FLY_REGISTRATION_FILE).read_text())
    for key in ("supplementary_tables", "peaks", "gff3", "rna_by_annotation"):
        entry = fly["registered_inputs"][key]
        if simpler.sha256(Path(entry["path"])) != entry["sha256"]:
            sys.exit(f"refused: simpler_fly's {key} differs from its registered input")
    chroms, rows, funnel, name_to_id = _rows()
    if sf.census(rows) != fly["census"] or funnel != fly["population_funnel"]:
        sys.exit("refused: simpler_fly's census does not recompute to its registered one")
    sets, summaries = _boundaries()
    arm = fn.read_arm_lengths(sf.GFF3)
    main = fn.attribute(rows, chroms, sets[fn.PRIMARY_SET])
    alts = {k: fn.attribute(rows, chroms, v) for k, v in sets.items() if k != fn.PRIMARY_SET}
    return chroms, main, alts, sets[fn.PRIMARY_SET], arm, summaries, name_to_id


def register() -> None:
    mf.trace_begin()
    _, main, alts, bounds, arm, summaries, _ = _build()
    cen = fn.census(main, alts, bounds, arm)
    out = {
        "written_before_any_rna_value_was_read": True,
        "registered_inputs": {
            "boundaries": {"path": str(fn.BOUNDARY_FILE), "sha256": simpler.sha256(fn.BOUNDARY_FILE)},
            "simpler_fly_registration": {
                "path": FLY_REGISTRATION_FILE,
                "sha256": simpler.sha256(Path(FLY_REGISTRATION_FILE)),
                "note": (
                    "its registered inputs (tables, peaks, GFF3, RNA by annotation) are checked by sha256"
                ),
            },
        },
        "boundary_lists": summaries,
        "arm_lengths": arm,
        "registration": fn.REGISTRATION,
        "census": cen,
        "alphagenome_requests": 0,
        "money": "none",
    }
    out["result_manifest"] = _manifest(
        [
            (str(fn.BOUNDARY_FILE), "Tables S4, S5 and S7 (boundary positions)"),
            (FLY_REGISTRATION_FILE, "simpler_fly's committed registration: its inputs' sha256 and census"),
            (str(sf.SUPP), "Table S2 and Table S4 only (ATAC cluster labels and marker calls)"),
            (str(sf.PEAKS), "the peak set"),
            (str(sf.GFF3), "gene and transcript positions and strands, and arm lengths"),
            (str(sf.RNA), "hashed only, against simpler_fly's registered sha256; no value is read here"),
        ],
        PARAMETERS,
        EXCLUSIONS,
        {
            "population": (
                "simpler_fly's, fixed by fd81e37; the primary subset fixed here before any RNA value"
            )
        },
    )
    print(save_result(REGISTRATION, out))
    print(json.dumps(cen, indent=1))


def _committed(rel: str) -> bool:
    logged = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", rel], cwd=ROOT, capture_output=True, text=True
    ).stdout.strip()
    clean = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", rel], cwd=ROOT).returncode == 0
    return bool(logged) and clean


def score() -> None:
    for rel in (REGISTRATION_FILE, FLY_REGISTRATION_FILE, *OWN_CODE):
        if not _committed(rel):
            sys.exit(f"refused: {rel} must be committed and unchanged before any RNA value is read")
    mf.trace_begin()
    reg = json.loads(Path(REGISTRATION_FILE).read_text())
    if reg["registration"] != json.loads(json.dumps(fn.REGISTRATION)):
        sys.exit("refused: simpler_flynode.REGISTRATION differs from the committed registration")
    for key, entry in reg["registered_inputs"].items():
        if simpler.sha256(Path(entry["path"])) != entry["sha256"]:
            sys.exit(f"refused: {key} differs from the registered input")
    chroms, main, alts, bounds, arm, summaries, name_to_id = _build()
    if fn.census(main, alts, bounds, arm) != reg["census"] or summaries != reg["boundary_lists"]:
        sys.exit("refused: the census does not recompute to the registered one")
    expr, rna_summary = sf.read_expression(sf.RNA, name_to_id)
    res = fn.score(chroms, main, alts, bounds, arm, expr)
    reg_sha = subprocess.run(
        ["git", "log", "-1", "--format=%h", "--", REGISTRATION_FILE], cwd=ROOT, capture_output=True, text=True
    ).stdout.strip()
    out = {
        "registration_commit": reg_sha,
        "registered_inputs": reg["registered_inputs"],
        "rna_summary": rna_summary,
        **res,
        "reading_is_the_registered_wording": (
            "the reading is REGISTRATION['outcome_wording'] with its numbers filled in "
            "and no other word changed"
        ),
        "alphagenome_requests": 0,
        "money": "none",
    }
    out["result_manifest"] = _manifest(
        [
            (REGISTRATION_FILE, "the committed registration"),
            (FLY_REGISTRATION_FILE, "simpler_fly's committed registration"),
            (str(fn.BOUNDARY_FILE), "Tables S4, S5 and S7 (boundary positions)"),
            (str(sf.SUPP), "Table S2 and Table S4 only"),
            (str(sf.PEAKS), "the peak set"),
            (str(sf.GFF3), "gene and transcript positions and strands, and arm lengths"),
            (str(sf.RNA), "the endpoint: UMI sums by annotation"),
        ],
        PARAMETERS,
        EXCLUSIONS,
        {"population": "pinned by the committed registration's census", "positive_control": "promoter"},
    )
    print(save_result(RESULT, out))
    print(out["outcome"])
    print(out.get("reading", ""))


if __name__ == "__main__":
    os.chdir(ROOT)
    {"fetch": fetch, "inventory": inventory, "register": register, "score": score}[sys.argv[1]]()
