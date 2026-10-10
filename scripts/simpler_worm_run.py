# SPDX-License-Identifier: AGPL-3.0-or-later
"""Simpler genomes, C. elegans first: inventory, census and registration, then the registered run.

    uv run python scripts/simpler_worm_run.py inventory   # data/results/simpler_worm_inventory.json
    uv run python scripts/simpler_worm_run.py fetch       # Table S2 into data/knowledge/simpler_worm
    uv run python scripts/simpler_worm_run.py register    # census + registration, no RNA read
    uv run python scripts/simpler_worm_run.py score       # refuses unless the registration is committed

No model request is sent and no key is read. docs/SIMPLER-GENOMES.md holds the reasoning.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import simpler  # noqa: E402
from genomeos.results import save_result  # noqa: E402

INVENTORY = "simpler_worm_inventory"
REGISTRATION = "simpler_worm_registration"
RESULT = "simpler_worm_result"
REGISTRATION_FILE = "data/results/simpler_worm_registration.json"
ENTRY = "scripts/simpler_worm_run.py"
OWN_CODE = (
    "genomeos/attribution/simpler.py",
    "scripts/simpler_worm_run.py",
    "tests/test_simpler_worm.py",
)

#: What public, measured C. elegans regulatory data exists, and what the project already held.
#: Sizes are bytes as served on 2026-10-10 where they were checked; "not checked" where they were not.
SOURCES = [
    {
        "source": simpler.SOURCE,
        "measures": (
            "ATAC-seq accessibility (RPM) of 47,514 elements in FACS-sorted nuclei of five young-"
            "adult tissues "
            "(germline, neurons, muscle, hypodermis, intestine) with a tissue call per element; "
            "element classes "
            "(20,271 putative_enhancer, 15,277 coding_promoter, 8,565 other_element, 2,244 "
            "unassigned_promoter, "
            "831 non-coding_RNA, 326 pseudogene_promoter); nuclear RNA-seq TPM of 20,222 protein-"
            "coding genes in "
            "the same five tissues with a tissue call per gene"
        ),
        "assembly": "ce11 (WBcel235)",
        "url": simpler.TABLE_S2_URL,
        "bytes": {"bundle_streamed_then_discarded": 28664943, "table_s2_kept": 14777075},
        "licence": "CC BY 4.0 (article and supplements)",
        "fetched": True,
        "used_by": "the registered test (elements and their tissue call; gene TPM as the endpoint)",
    },
    {
        "source": "GEO GSE141213 (Serizay et al. 2020 raw and processed tracks)",
        "measures": (
            "per-tissue ATAC and stranded RNA bigWigs, young adult, plus L1/L3 muscle and L3 germline ATAC"
        ),
        "assembly": "ce11",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141213/suppl/",
        "bytes": {"raw_tar": 2386268160, "single_bigwig_range": [7135176, 236553287]},
        "licence": "public GEO record",
        "fetched": False,
        "used_by": "nothing: Table S2 already holds the per-element and per-gene values distilled from these",
    },
    {
        "source": (
            "Janes J et al. (2018) Chromatin accessibility dynamics across C. elegans development and"
            " ageing. eLife 7:e37344"
        ),
        "measures": (
            "ATAC-seq and long-cap / nuclear RNA across six stages and ageing; 42,245 elements "
            "annotated by type"
        ),
        "assembly": "ce11",
        "url": "https://elifesciences.org/articles/37344",
        "bytes": "not checked",
        "licence": "CC BY 4.0 (eLife)",
        "fetched": False,
        "used_by": (
            "nothing: Table S2 extends the same annotation pipeline (its header says so) to 47,514 elements"
        ),
    },
    {
        "source": "UCSC ce11 phyloP135way (135 nematodes) and phyloP26way / phastCons26way (26 nematodes)",
        "measures": "per-base conservation scores across nematodes",
        "assembly": "ce11",
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/ce11/phyloP135way/ce11.phyloP135way.bw",
        "bytes": {"phyloP135way_bw": 356259120, "phyloP26way_bw": 180170767, "phastCons26way_bw": 159951981},
        "licence": "UCSC Genome Browser download; terms for this track not checked",
        "fetched": False,
        "used_by": (
            "nothing yet: readable by range request with genomeos/attribution/bigwig.py without download, "
            "for the worm constraint tiers (docs/SIMPLER-GENOMES.md)"
        ),
    },
    {
        "source": (
            "Durham TJ et al. (2021) Comprehensive characterization of tissue-specific chromatin "
            "accessibility in L2 C. elegans nematodes. Genome Res 31:1952"
        ),
        "measures": "single-cell (sci-ATAC) accessibility at L2, cell-type resolved",
        "assembly": "ce11",
        "url": "https://genome.cshlp.org/content/31/10/1952",
        "bytes": "not checked",
        "licence": "not checked",
        "fetched": False,
        "used_by": "nothing: a candidate replication of the registered test at a second stage",
    },
    {
        "source": (
            "Daugherty AC et al. (2017) Chromatin accessibility dynamics reveal novel functional "
            "enhancers in C. elegans. Genome Res 27:2096"
        ),
        "measures": (
            "ATAC-seq across embryo, L3 and young adult, with transgenic reporter tests of selected elements"
        ),
        "assembly": "ce11",
        "url": "https://genome.cshlp.org/content/27/12/2096",
        "bytes": "not checked",
        "licence": "not checked",
        "fetched": False,
        "used_by": "nothing",
    },
    {
        "source": (
            "WormBase cis-regulatory features (enhancer, promoter, regulatory_region, TF_binding_site"
            " in the WormBase GFF3)"
        ),
        "measures": "curated and imported regulatory annotations, mostly from reporter and ChIP studies",
        "assembly": "WBcel235",
        "url": "https://downloads.wormbase.org/releases/",
        "bytes": "not checked",
        "licence": "WormBase terms, not checked",
        "fetched": False,
        "used_by": "nothing",
    },
    {
        "source": "Genome-wide reporter assays (STARR-seq, MPRA) in C. elegans",
        "measures": "none found in one web search on 2026-10-10",
        "assembly": None,
        "url": None,
        "bytes": None,
        "licence": None,
        "fetched": False,
        "used_by": "nothing",
    },
    {
        "source": (
            "Further tissue or cell RNA: Kaletsky et al. 2018 (adult tissues), CeNGEN (Taylor et al. "
            "2021, neuron classes), modERN TF ChIP-seq (Kudron et al. 2018)"
        ),
        "measures": "tissue and neuron-class expression; transcription-factor binding",
        "assembly": "ce11",
        "url": None,
        "bytes": "not checked",
        "licence": "not checked",
        "fetched": False,
        "used_by": "nothing",
    },
]

HELD = [
    {
        "path": "data/reference/celegans/WBcel235.dna.toplevel.fa.gz",
        "bytes": 30316631,
        "what": "WBcel235 genome (Ensembl), git-ignored",
    },
    {
        "path": "data/reference/celegans/WBcel235.63.gff3.gz",
        "bytes": 8231538,
        "what": "Ensembl gene models, 19,985 protein-coding genes; the attribution reads this",
    },
    {
        "path": "data/knowledge/celegans/ma2021_atlas.zip",
        "bytes": 88347568,
        "what": "Ma et al. 2021 TF reporter expression per embryonic cell (Zenodo 4737593, CC BY 4.0)",
    },
    {
        "path": "data/results/celegans_*.json",
        "bytes": None,
        "what": "16 results: lineage, Packer 2019 cell types, contacts, fates, TF atlas, time axis",
    },
    {
        "path": "data/organisms/celegans/*.bio",
        "bytes": None,
        "what": "lineage, fates, contacts as BioLang (held by lane-wormcontact)",
    },
]


def _manifest(
    inputs: list[tuple[str, str]], parameters: dict, exclusions: list[str], partitions: Any
) -> dict:
    return {
        "sources": [
            {
                "accession": simpler.SOURCE,
                "version": (
                    "Genome Res 30(12), December 2020; supplement as served by Europe PMC on 2026-10-10"
                ),
            },
            {
                "accession": "Ensembl gene models for WBcel235 (WBcel235.63.gff3.gz)",
                "version": "release 63, the copy held in data/reference/celegans since 2026-09-10",
            },
        ],
        "inputs": [mf.input_entry(path, partition=part) for path, part in inputs],
        "assembly": "WBcel235 (ce11): Table S2 and the gene models are on the same build; nothing is lifted",
        "coordinates": {
            "base": 1,
            "interval": "closed",
            "note": "the GFF3 is 1-based closed; Table S2 does not state its base, and the element "
            "midpoint (start + end) // 2 moves by at most 1 bp either way",
        },
        "parameters": parameters,
        "exclusions": exclusions,
        "partitions": partitions,
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


PARAMETERS = {
    "element_type": simpler.ELEMENT_TYPE,
    "promoter_type": simpler.PROMOTER_TYPE,
    "tissues": list(simpler.TISSUES),
    "detection_tpm": simpler.DETECTION_TPM,
    "distance_bins_bp": list(simpler.DISTANCE_BINS),
    "control_draws": simpler.CONTROL_DRAWS,
    "bootstrap": simpler.BOOTSTRAP,
    "cluster_bp": simpler.CLUSTER_BP,
    "seed": simpler.SEED,
    "min_discordant": simpler.MIN_DISCORDANT,
    "alpha": simpler.ALPHA,
}
EXCLUSIONS = [
    "elements of any class other than putative_enhancer (primary) and coding_promoter (positive control)",
    "elements whose tissue call is not exactly one of the five tissues "
    "(Ubiq., Low, Unclassified, combinations)",
    "genes that are not protein_coding in the GFF3, and the MtDNA genes",
    "the source's own element-to-gene columns, never read: " + ", ".join(simpler.NEVER_READ),
]
TABLE = str(simpler.TABLE_S2)
GFF = str(simpler.GFF3)


def inventory() -> None:
    mf.trace_begin()
    out = {
        "question": (
            "what public, measured C. elegans regulatory data exists, and what the project already holds"
        ),
        "sources": SOURCES,
        "already_held": HELD,
        "held_before_this_lane": (
            "no worm regulatory-element, accessibility, conservation or tissue-RNA data; the worm data "
            "the project held were the lineage, cell types, contacts and the TF reporter atlas"
        ),
        "table_s2_sha256": simpler.sha256(simpler.TABLE_S2),
        "lane_download_bytes": 28664943,
        "money": "none; no model request",
    }
    out["result_manifest"] = _manifest(
        [(TABLE, "the one source fetched, digested here")],
        {"download_budget_bytes": 1_500_000_000},
        ["sources whose size or licence was not checked are listed with 'not checked'"],
        "n/a: an inventory has no evaluation partition",
    )
    print(save_result(INVENTORY, out))


def fetch() -> None:
    if simpler.TABLE_S2.exists():
        print("present", simpler.TABLE_S2, simpler.sha256(simpler.TABLE_S2))
        return
    print(json.dumps(simpler.fetch_table()))


def _rows():
    chroms = simpler.read_genes(simpler.GFF3)
    elements = simpler.read_elements(simpler.TABLE_S2)
    return chroms, simpler.attribute(elements, chroms)


def register() -> None:
    mf.trace_begin()
    chroms, rows = _rows()
    cen = simpler.census(rows)
    out = {
        "written_before_the_rna_sheet_was_read": True,
        "registered_inputs": {
            "table_s2": {"path": TABLE, "sha256": simpler.sha256(simpler.TABLE_S2)},
            "gff3": {"path": GFF, "sha256": simpler.sha256(simpler.GFF3)},
            "protein_coding_genes": sum(len(c.genes) for c in chroms.values()),
        },
        "registration": simpler.REGISTRATION,
        "census": cen,
        "alphagenome_requests": 0,
        "money": "none",
    }
    out["result_manifest"] = _manifest(
        [
            (TABLE, "elements and their tissue call only (sheet ATAC_metrics)"),
            (GFF, "gene positions and strands"),
        ],
        PARAMETERS,
        EXCLUSIONS,
        {"population": "fixed here, before the RNA sheet (RNA_metrics) was opened"},
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
    if not _committed(REGISTRATION_FILE) or not _committed("genomeos/attribution/simpler.py"):
        sys.exit(
            "refused: the registration and simpler.py must be committed and unchanged before the RNA "
            "sheet is read"
        )
    mf.trace_begin()
    reg = json.loads(Path(REGISTRATION_FILE).read_text())
    if reg["registration"] != json.loads(json.dumps(simpler.REGISTRATION)):
        sys.exit("refused: simpler.REGISTRATION differs from the committed registration")
    for key in ("table_s2", "gff3"):
        if (
            simpler.sha256(Path(reg["registered_inputs"][key]["path"]))
            != reg["registered_inputs"][key]["sha256"]
        ):
            sys.exit(f"refused: {key} differs from the registered input")
    chroms, rows = _rows()
    expr = simpler.read_expression(simpler.TABLE_S2)
    res = simpler.score(chroms, rows, expr)
    reg_sha = subprocess.run(
        ["git", "log", "-1", "--format=%h", "--", REGISTRATION_FILE],
        cwd=ROOT,
        capture_output=True,
        text=True,
    ).stdout.strip()
    out = {
        "registration_commit": reg_sha,
        "registered_inputs": reg["registered_inputs"],
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
            (TABLE, "elements (ATAC_metrics) and the endpoint (RNA_metrics)"),
            (GFF, "gene positions and strands"),
        ],
        PARAMETERS,
        EXCLUSIONS,
        {
            "population": "pinned by the committed registration's census",
            "positive_control": "coding_promoter",
        },
    )
    print(save_result(RESULT, out))
    print(out["outcome"])
    print(out.get("reading", ""))


if __name__ == "__main__":
    os.chdir(ROOT)
    {"inventory": inventory, "fetch": fetch, "register": register, "score": score}[sys.argv[1]]()
