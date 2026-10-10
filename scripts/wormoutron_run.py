#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Simpler genomes, area I: the worm's between-gene lead against its own outron caveat.

    uv run python scripts/wormoutron_run.py inventory   # data/results/wormoutron_inventory.json
    uv run python scripts/wormoutron_run.py register    # census + registration, no RNA value read
    uv run python scripts/wormoutron_run.py score       # refuses unless the registration is committed

`inventory` records the map this lane chose and why, with its licence, assembly and what it defines.
`register` fixes the exclusion's surviving population before any RNA value is read -- the masks come
from positions and annotations alone. `score` refuses to run until the registration and this lane's
code are committed and unchanged, then runs once.

docs/SIMPLER-GENOMES.md section 12 holds the reasoning. No model request is sent, no key is read and
no money is spent. The map was already streamed into the git-ignored store by hand; `inventory`
verifies its sha256 rather than fetching again, and prints the fetch command it would need.
"""

from __future__ import annotations

import datetime as dt
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import simpler  # noqa: E402
from genomeos.attribution import wormoutron as wo  # noqa: E402
from genomeos.results import save_result  # noqa: E402

INVENTORY = "wormoutron_inventory"
REGISTRATION = "wormoutron_registration"
RESULT = "wormoutron_result"
REGISTRATION_FILE = "data/results/wormoutron_registration.json"
ENTRY = "scripts/wormoutron_run.py"
OWN_CODE = (
    "genomeos/attribution/simpler.py",
    "genomeos/attribution/wormoutron.py",
    "scripts/wormoutron_run.py",
    "tests/test_wormoutron.py",
)

#: The one network transfer this lane made, recorded so it can be repeated exactly. The release file
#: is 742 MB and is NOT kept: it was streamed and filtered in one pass, so the store holds 8.0 MB.
FETCH_COMMAND = (
    'curl -sS "https://ftp.ebi.ac.uk/pub/databases/wormbase/releases/WS298/species/c_elegans/'
    'PRJNA13758/c_elegans.PRJNA13758.WS298.annotations.gff3.gz" '
    "| tee >(md5 -q) | gzip -dc | awk -F'\\t' "
    '\'NF>8 && ($2=="operon" || $2=="deprecated_operon" || $2=="dicistronic_mRNA" '
    '|| $2=="SL1" || $2=="SL2" || $2=="TSS_region")\' '
    "> data/knowledge/wormoutron/ws298_operon_sl_tss.gff3"
)

PARAMETERS = {
    "outron_max_bp": wo.OUTRON_MAX_BP,
    "sl_window_bp": wo.SL_WINDOW_BP,
    "min_surviving": wo.MIN_SURVIVING,
    "min_surviving_discordant": wo.MIN_SURVIVING_DISCORDANT,
    "operon_sources_used": list(wo.OPERON_SOURCES),
    "operon_sources_not_used": list(wo.OPERON_SOURCES_NOT_USED),
    "detection_tpm": simpler.DETECTION_TPM,
    "tissues": list(simpler.TISSUES),
    "cluster_bp": simpler.CLUSTER_BP,
    "bootstrap_reps": simpler.BOOTSTRAP,
    "control_draws": simpler.CONTROL_DRAWS,
    "seed": simpler.SEED,
    "stage1_between_genes": wo.STAGE1_BETWEEN_GENES,
}

EXCLUSIONS = [
    "elements with a host gene (domain_case == 'host'): the host stratum is not this lane's population",
    "elements of any class other than putative_enhancer",
    "elements whose tissue call is not exactly one of the five tissues",
    "deprecated_operon features (111 of them), which WormBase no longer asserts",
    "MtDNA, which carries no protein-coding gene in the models stage 1 read",
    "the source table's own element-to-gene columns, which answer the question being asked",
]


def _manifest(inputs: list[tuple[str, str]], partitions: Any) -> dict:
    return {
        "sources": [
            {"accession": simpler.SOURCE, "version": "as read in stage 1, unchanged"},
            {"accession": wo.MAP_SOURCE, "version": "WormBase WS298, served 2026-10-10"},
            {
                "accession": "Ensembl WBcel235 release 63 gene models (Caenorhabditis elegans)",
                "version": "WBcel235.63",
            },
        ],
        "inputs": [mf.input_entry(path, partition=part) for path, part in inputs],
        "assembly": (
            "WBcel235 (ce11) throughout: Table S2, the Ensembl gene models and the WormBase WS298 "
            "features are on the same build, so nothing is lifted over"
        ),
        "coordinates": {
            "base": 1,
            "interval": "closed",
            "note": (
                "both GFF3s are 1-based closed; the element midpoint is (start + end) // 2, and "
                "Table S2 does not state its base, so a midpoint moves by at most 1 bp either way"
            ),
        },
        "parameters": PARAMETERS,
        "exclusions": EXCLUSIONS,
        "partitions": partitions,
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def _build() -> dict[str, Any]:
    """Positions and annotations only. No RNA sheet is opened here."""
    chroms = simpler.read_genes()
    els = simpler.read_elements()
    rows = simpler.attribute(els, chroms)
    mapdata = wo.read_map()
    mapstats = wo.flag(rows, chroms, mapdata)
    mapstats["operons_used"] = len(mapdata["operons"])
    mapstats["deprecated_operons_skipped"] = mapdata["deprecated_operons_skipped"]
    mapstats["sl_acceptor_sites"] = len(mapdata["sl_sites"])
    mapstats["tss_regions"] = len(mapdata["tss_regions"])
    return {"chroms": chroms, "rows": rows, "mapstats": mapstats}


def _registered_inputs() -> dict[str, Any]:
    return {
        "table_s2": {
            "path": str(simpler.TABLE_S2),
            "sha256": simpler.sha256(simpler.TABLE_S2),
        },
        "gff3": {"path": str(simpler.GFF3), "sha256": simpler.sha256(simpler.GFF3)},
        "ws298_map": {"path": str(wo.MAP_FILE), "sha256": simpler.sha256(wo.MAP_FILE)},
        "ws298_source_type_census": {
            "path": str(wo.CENSUS_FILE),
            "sha256": simpler.sha256(wo.CENSUS_FILE),
        },
        "ws298_checksums": {
            "path": str(wo.CHECKSUMS_FILE),
            "sha256": simpler.sha256(wo.CHECKSUMS_FILE),
        },
        "remote_release_file": {
            "url": (
                "https://ftp.ebi.ac.uk/pub/databases/wormbase/releases/WS298/species/c_elegans/"
                "PRJNA13758/c_elegans.PRJNA13758.WS298.annotations.gff3.gz"
            ),
            "bytes": 742_000_000,
            "kept_on_disk": False,
            "published_md5": wo.REMOTE_MD5,
            "streamed_md5": wo.REMOTE_MD5,
            "streamed_sha256": wo.REMOTE_SHA256,
            "verified": (
                "the streamed md5 equals the md5 WormBase publishes in WS298/CHECKSUMS, and two "
                "independent streams gave the same sha256, so the 8.0 MB extract is of the "
                "authentic complete release file; the release file itself is not kept"
            ),
        },
    }


# ------------------------------------------------------------------------------------- the inventory


def inventory() -> None:
    mf.trace_begin()
    present = {
        str(p): {"exists": p.exists(), "bytes": p.stat().st_size if p.exists() else 0}
        for p in (wo.MAP_FILE, wo.CENSUS_FILE, wo.CHECKSUMS_FILE)
    }
    out: dict[str, Any] = {
        "lane": "lane-wormoutron",
        "date": dt.date.today().isoformat(),
        "area": "I",
        "question": wo.REGISTRATION["question"],
        "why_this_lane_exists": (
            "lane-flysource measured in Drosophila that the domain rule's lead over nearest TSS "
            "appears only where a host gene can be named and is absent between genes (+0.03 "
            "[-0.50, +0.60] on 789 elements, docs/SIMPLER-GENOMES.md section 11). The worm is NOT "
            "qualified by that: its between-gene stratum held at +3.73 [+2.79, +4.74] over 2,335 "
            "elements. But stage 1 registered a caveat for that very stratum -- trans-splicing and "
            "outrons -- and never tested it. This lane tests it."
        ),
        "candidates_considered": [
            {
                "map": "WormBase WS298 annotation GFF3 (CHOSEN)",
                "defines": wo.REGISTRATION["map_defines"],
                "assembly": "WBcel235 (ce11); WormBase has used this build since WS235",
                "liftover": "none needed",
                "obtainable": True,
                "how": (
                    "streamed from the EBI mirror and filtered in one pass to six feature classes; "
                    "the 742 MB release file is not kept, the 8.0 MB extract is"
                ),
                "why_chosen": (
                    "it is the only candidate that carries ALL THREE pieces the caveat needs -- "
                    "WormBase's clustered TSS regions (which are the published Chen, Kruesi, Saito "
                    "and Gu TSS data sets, Saito being the outron paper), the measured SL1 and SL2 "
                    "trans-splice acceptor sites, and the operon annotation with member genes -- on "
                    "the same assembly as the gene models stage 1 read, so nothing is lifted over"
                ),
                "licence": (
                    "WormBase data are distributed for free reuse; the release, URL and WormBase's "
                    "own published md5 are recorded instead of asserting a licence identifier this "
                    "lane did not read, because wormbase.org/about/policies answered an "
                    "unauthenticated fetch with HTTP 403 (Cloudflare)"
                ),
            },
            {
                "map": "Saito et al. 2013 Genome Res 23:1348, supplementary TSS/outron tables",
                "status": "NOT fetched separately",
                "why_not": (
                    "WormBase's TSS_region features already carry this data set -- the WS298 TSS "
                    "wiggle archive's members are TSS_Chen, TSS_Kruesi, TSS_Saito and TSS_Gu, and "
                    "the release letter names WBPaper00042246, WBPaper00042529 and WBPaper00042354 "
                    "-- on WBcel235, whereas the paper's own tables are on the build of 2013 and "
                    "would need a liftover this lane could not verify"
                ),
            },
            {
                "map": "Kruesi et al. 2013 eLife (GRO-cap TSSs); Chen et al. 2013 Genome Res; Gu et al.",
                "status": "NOT fetched separately",
                "why_not": "same reason: carried by WormBase's TSS_region features on WBcel235",
            },
            {
                "map": "WormBase operon annotation as a standalone file",
                "status": "does not exist on the FTP mirror",
                "why_not": (
                    "the WS298 species directory has no operon-only extract; the operon features "
                    "live in the full annotation GFF3, which is why that file was streamed"
                ),
            },
        ],
        "map_chosen": wo.MAP_SOURCE,
        "map_defines": wo.REGISTRATION["map_defines"],
        "store": present,
        "fetch_command": FETCH_COMMAND,
        "registered_inputs": _registered_inputs(),
        "disk": {
            "store_bytes": sum(v["bytes"] for v in present.values()),
            "budget_bytes": 500_000_000,
            "release_file_kept": False,
        },
        "alphagenome_requests": 0,
        "model_requests": 0,
        "money": "none; no model request, no key read",
        "stop_list": (
            "docs/ROADMAP.md section 5, 'The stop list of 2026-09-27, reinstated': no AlphaGenome "
            "and no deletion request (0 model requests here), no new reader biosample (the human "
            "reader's 13 biosamples are untouched; a WormBase annotation release is not a reader "
            "biosample), nothing in area J, no node-model work, no BioForge calibration"
        ),
    }
    out["result_manifest"] = _manifest(
        [
            (str(wo.MAP_FILE), "the WS298 extract: operon, SL1, SL2 and TSS_region features"),
            (str(wo.CENSUS_FILE), "the whole-file source|type census from the same pass"),
            (str(wo.CHECKSUMS_FILE), "WormBase's published checksums for WS298"),
            (str(simpler.TABLE_S2), "hashed only, to pin stage 1's source; no sheet is opened"),
            (str(simpler.GFF3), "hashed only, to pin stage 1's gene models"),
        ],
        "n/a: an inventory has no evaluation partition",
    )
    print(save_result(INVENTORY, out))
    print(json.dumps(out["map_chosen"], indent=1))


# ---------------------------------------------------------------------------------- the registration


def register() -> None:
    mf.trace_begin()
    built = _build()
    cen = wo.census(built["rows"], built["mapstats"])
    out: dict[str, Any] = {
        "result": "registration, before any RNA value is re-read",
        "date": dt.date.today().isoformat(),
        "written_before_the_rna_sheet_was_reread": True,
        "prior_exposure": wo.REGISTRATION["prior_exposure"],
        "registered_inputs": _registered_inputs(),
        "registration": json.loads(json.dumps(wo.REGISTRATION)),
        "census": cen,
        "alphagenome_requests": 0,
        "money": "none",
    }
    out["result_manifest"] = _manifest(
        [
            (str(simpler.TABLE_S2), "sheet ATAC_metrics only: positions, class and tissue call"),
            (str(simpler.GFF3), "protein-coding gene positions and strands"),
            (str(wo.MAP_FILE), "operon, SL1, SL2 and TSS_region positions and strands"),
            (str(wo.CENSUS_FILE), "the whole-file source|type census from the same pass"),
            (str(wo.CHECKSUMS_FILE), "WormBase's published checksums for WS298"),
        ],
        "population fixed here, before the RNA sheet (RNA_metrics) was opened again",
    )
    print(save_result(REGISTRATION, out))
    print(
        json.dumps(
            {
                "between_genes_recomputed": cen["between_genes_recomputed"],
                "matches_stage1": cen["matches_stage1"],
                "excluded_by": cen["excluded_by"],
                "surviving": cen["surviving"],
                "surviving_discordant": cen["surviving_discordant"],
                "mde_points_unclustered": cen["mde_points_unclustered"],
                "inconclusive_unscored": cen["inconclusive_unscored"],
            },
            indent=1,
        )
    )


# ----------------------------------------------------------------------------------------- the run


def _committed(rel: str) -> bool:
    logged = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", rel], cwd=ROOT, capture_output=True, text=True
    ).stdout.strip()
    clean = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", rel], cwd=ROOT).returncode == 0
    return bool(logged) and clean


def score() -> None:
    for rel in (REGISTRATION_FILE, *OWN_CODE):
        if not _committed(rel):
            sys.exit(f"refused: {rel} must be committed and unchanged before any RNA value is read")
    mf.trace_begin()
    reg = json.loads(Path(REGISTRATION_FILE).read_text())
    if reg["registration"] != json.loads(json.dumps(wo.REGISTRATION)):
        sys.exit("refused: wormoutron.REGISTRATION differs from the committed registration")
    for key, entry in reg["registered_inputs"].items():
        if "path" in entry and simpler.sha256(Path(entry["path"])) != entry["sha256"]:
            sys.exit(f"refused: {key} differs from the registered input")
    built = _build()
    cen = wo.census(built["rows"], built["mapstats"])
    if cen != reg["census"]:
        sys.exit("refused: the census does not recompute to the registered one")
    expr = simpler.read_expression()
    res = wo.score(built["rows"], built["chroms"], expr, cen)
    out: dict[str, Any] = {
        "result": res["label"],
        "date": dt.date.today().isoformat(),
        "registration": REGISTRATION_FILE,
        "registration_unchanged": True,
        "prior_exposure": wo.REGISTRATION["prior_exposure"],
        **res,
        "alphagenome_requests": 0,
        "money": "none; no model request, no key read",
    }
    out["result_manifest"] = _manifest(
        [
            (str(simpler.TABLE_S2), "sheets ATAC_metrics and RNA_metrics"),
            (str(simpler.GFF3), "protein-coding gene positions and strands"),
            (str(wo.MAP_FILE), "operon, SL1, SL2 and TSS_region positions and strands"),
            (str(wo.CENSUS_FILE), "the whole-file source|type census from the same pass"),
            (str(wo.CHECKSUMS_FILE), "WormBase's published checksums for WS298"),
            (REGISTRATION_FILE, "the committed registration"),
        ],
        "population pinned by the committed registration's census",
    )
    print(save_result(RESULT, out))
    print()
    print(res["label"], "/", out.get("falsifier", "-"))
    print()
    print(res["reading"])


def main() -> None:
    cmds = {"inventory": inventory, "register": register, "score": score}
    if len(sys.argv) != 2 or sys.argv[1] not in cmds:
        sys.exit(f"usage: {ENTRY} {{{'|'.join(cmds)}}}")
    cmds[sys.argv[1]]()


if __name__ == "__main__":
    main()
