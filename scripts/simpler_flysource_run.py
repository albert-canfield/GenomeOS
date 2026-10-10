# SPDX-License-Identifier: AGPL-3.0-or-later
"""Simpler genomes, area I: the fly's between-gene question on a source that reaches the stop rule.

    uv run python scripts/simpler_flysource_run.py fetch       # inputs into the git-ignored store
    uv run python scripts/simpler_flysource_run.py inventory   # data/results/simpler_flysource_inventory.json
    uv run python scripts/simpler_flysource_run.py register    # census + registration, no RNA value read
    uv run python scripts/simpler_flysource_run.py score       # refuses unless the registration is committed

`inventory` answers the lane's first question -- which candidate boundary sets and sources exist,
under what licence, at what stage, and how many between-gene elements each would move -- and settles
route (a) by counting, with the boundary density sweep. `register` fixes the test before any RNA
value is read. docs/SIMPLER-GENOMES.md section 11 holds the reasoning. No model request is sent and
no key is read.
"""

from __future__ import annotations

import datetime as dt
import gzip
import hashlib
import io
import json
import os
import shutil
import subprocess
import sys
import tarfile
import urllib.request
import zipfile
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import simpler  # noqa: E402
from genomeos.attribution import simpler_fly as sf  # noqa: E402
from genomeos.attribution import simpler_flynode as fn  # noqa: E402
from genomeos.attribution import simpler_flysource as fs  # noqa: E402
from genomeos.results import save_result  # noqa: E402

INVENTORY = "simpler_flysource_inventory"
REGISTRATION = "simpler_flysource_registration"
RESULT = "simpler_flysource_result"
REGISTRATION_FILE = "data/results/simpler_flysource_registration.json"
FLYNODE_REGISTRATION_FILE = "data/results/simpler_flynode_registration.json"
ENTRY = "scripts/simpler_flysource_run.py"
OWN_CODE = (
    "genomeos/attribution/simpler.py",
    "genomeos/attribution/simpler_fly.py",
    "genomeos/attribution/simpler_flynode.py",
    "genomeos/attribution/simpler_flysource.py",
    "scripts/simpler_flysource_run.py",
    "tests/test_simpler_flysource.py",
)
BUDGET_BYTES = 1_000_000_000

#: Bytes moved while choosing the source and the sets, outside the fetch: Europe PMC's supplementary
#: archive for Dolsten et al. 2025 (59,258,286) and for Kaushal et al. 2022 (22,355,160), the
#: Janssens peak file (4,774,826) and portal page (87,174). HEAD requests and GEO file lists are not
#: counted.
SCOUTING_BYTES = 59_258_286 + 22_355_160 + 4_774_826 + 87_174

#: Each input: where it comes from, what is read from it, and under what licence.
INPUTS: dict[str, dict[str, Any]] = {
    "cellTypeInfo.tsv": {
        "url": "https://flybrain.aertslab.org/downloads/cellTypeInfo.tsv",
        "what": "the source's own accessibility cell type -> RNA cluster map (`RNAcluster` column)",
        "licence": fs.ELEMENT_LICENCE,
        "member": None,
    },
    "tbl_DARs.txt": {
        "url": "https://flybrain.aertslab.org/downloads/other/tbl_DARs.txt.zip",
        "what": (
            "the source's differentially accessible regions; only the lvl1 / vsRest / WholeBrain / Up "
            "rows and only their region, cell type and selection columns are read"
        ),
        "licence": fs.ELEMENT_LICENCE,
        "member": "tar:tbl_DARs.txt",
    },
    "Davie_Janssens_Koldere_et_al_2018_AdultBrain.loom": {
        "url": (
            "https://flybrain.aertslab.org/downloads/loom/Davie_Janssens_Koldere_et_al_2018_AdultBrain.loom"
        ),
        "what": "the endpoint: gene counts per cell, with each cell's `cell_type`",
        "licence": fs.RNA_LICENCE,
        "member": None,
    },
    "dolsten2025_TableS2.xlsx": {
        "url": (
            "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12648115/supplementaryFiles"
            "?includeInlineImage=false"
        ),
        "what": "Table S2, sheet `boundaries`: the Micro-C atlas's 3,521 boundary spans (beside)",
        "licence": "CC BY-NC-ND 4.0 (Cell Genomics open-access article)",
        "member": "zip:mmc3.zip!TableS2.xlsx",
    },
    "kaushal2022_data_file_s2.xlsx": {
        "url": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9106302/supplementaryFiles",
        "what": "data file S2: the 2-kb bins marked `exist.WT_0` in wild-type 2-6 h embryos (beside)",
        "licence": "CC BY 4.0 (Science Advances)",
        "member": "zip:sciadv.abl8834_data_files_s1_to_s15.zip!sciadv.abl8834_data_file_s2.xlsx",
    },
    "kaushal2021_larval_cns_data_s3.xlsx": {
        "url": (
            "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-21366-2/"
            "MediaObjects/41467_2021_21366_MOESM6_ESM.xlsx"
        ),
        "what": (
            "supplementary data 3: the 2-kb bins marked `exist.wildtype` in wild-type third-instar "
            "larval central nervous system -- the one set whose tissue matches the elements' (beside)"
        ),
        "licence": "CC BY 4.0 (Nature Communications)",
        "member": None,
    },
}

#: The candidate boundary sets. `counted` says whether the between-gene numbers were COMPUTED here
#: from positions; nothing qualifies on a count that could not be made.
BOUNDARY_CANDIDATES = [
    {
        "name": "Stadler, Haines, Eisen 2017, eLife 6:e29550, Table S4 (manual)",
        "stage": "2-3 h embryos (nuclear cycle 14)",
        "assembly": "dm6 as given",
        "licence": "CC BY 4.0",
        "boundaries": 3121,
        "bin_bp": 500,
        "bytes": 386_317,
        "counted": True,
        "key": "stadler_manual_nc14",
        "note": "lane-flynode's set, re-counted here; the primary set of this registration",
    },
    {
        "name": "Dolsten et al. 2025, Cell Genomics 5:101002, Table S2",
        "stage": (
            "a Micro-C atlas of nc1-8, nc14 and stages 10-12; the latest stage is about 6-8 h, so "
            "this is the only counted set that reaches the marker window's lower end"
        ),
        "assembly": "BDGP Release 6 (dm6)",
        "licence": "CC BY-NC-ND 4.0",
        "boundaries": 3521,
        "bin_bp": 400,
        "bytes": 1_087_877,
        "counted": True,
        "key": "dolsten_2025_atlas",
        "note": (
            "the authors' own union list; per-stage calls would need a threshold chosen here, so the "
            "union is used as shipped. The text states 3,673 boundaries and the sheet holds 3,521"
        ),
    },
    {
        "name": "Kaushal et al. 2022, Sci Adv 8:eabl8834, data file S2",
        "stage": "2-6 h embryos",
        "assembly": "dm6",
        "licence": "CC BY 4.0",
        "boundaries": 2334,
        "bin_bp": 2000,
        "bytes": 183_298,
        "counted": True,
        "key": "kaushal_2022_wt_embryo",
        "note": "TopDom; the authors' own `exist.WT_0` call, no threshold chosen here",
    },
    {
        "name": "Kaushal et al. 2021, Nat Commun 12:1011, supplementary data 3",
        "stage": "third-instar larval central nervous system (not embryos)",
        "assembly": "dm6",
        "licence": "CC BY 4.0",
        "boundaries": 3458,
        "bin_bp": 2000,
        "bytes": 140_228,
        "counted": True,
        "key": "kaushal_2021_wt_larval_cns",
        "note": (
            "the only candidate whose TISSUE matches the element source's (brain), so it is the "
            "stage-and-tissue-matched set this lane found; beside here only because the primary set "
            "is carried from f5985f7 unchanged"
        ),
    },
    {
        "name": "Mohana et al. 2023, Cell 186:3826, GEO GSE214707",
        "stage": "6-8 h, 10-12 h AND 14-16 h embryos, FACS-sorted neuroblast, glia and neuron",
        "assembly": "dm6 expected; not confirmed",
        "licence": "the publisher's user licence, not Creative Commons; the GEO data state no terms",
        "boundaries": None,
        "bin_bp": None,
        "bytes": 4_009_000_960,
        "counted": False,
        "why_not": (
            "THE KEY FINDING OF ROUTE (a): this is the one study whose Hi-C covers 6-16 h embryos, in "
            "the very cell types the question is about, but it publishes LOOP anchors and "
            "meta-domains and no boundary or insulation list. Boundaries would have to be called from "
            "its .mcool matrices (about 1.07 GB for the three embryonic stages), which chooses a "
            "caller and a threshold here -- the free parameter lane-flynode refused"
        ),
    },
    {
        "name": "Schuettengruber et al. 2014, Cell Reports 9:219, GEO GSE61471",
        "stage": "16-18 h embryos",
        "assembly": "dm3",
        "licence": "CC BY-NC-ND 3.0 at the article; the GEO file states no terms",
        "boundaries": 1305,
        "bin_bp": None,
        "bytes": 11_396,
        "counted": False,
        "why_not": (
            "a genuine late-embryo author-called domain list, but 1,305 domains is a LOWER density "
            "than every set counted here, and the density sweep's lowest rung (2,000 boundaries) "
            "moves 65 on average and 69 at best, so it cannot reach 200; and its dm3 coordinates "
            "would need a liftover, which is a tool and an error source this lane does not add"
        ),
    },
    {
        "name": "Sexton et al. 2012, Cell 148:458, Table S1",
        "stage": "16-18 h embryos (deep map); 4-12 h pilot",
        "assembly": "dm3",
        "licence": "the publisher's open-archive user licence, not Creative Commons",
        "boundaries": 1169,
        "bin_bp": 10000,
        "bytes": 122_880,
        "counted": False,
        "why_not": "1,169 domains: the same density argument as the row above, and dm3",
    },
    {
        "name": "Pollex et al. 2024, Mol Cell 84:822",
        "stage": "2-4 h, 6-8 h and 16-18 h embryos",
        "assembly": "dm6",
        "licence": "CC BY 4.0",
        "boundaries": None,
        "bin_bp": None,
        "bytes": 107_268,
        "counted": False,
        "why_not": "gene-gene loops only: 1,213 loops across the three stages, no boundary list",
    },
    {
        "name": "Ogiyama et al. 2018, Mol Cell 71:73, GEO GSE103625",
        "stage": "cycles 1-8 to embryonic stage 16 (about 13-16 h)",
        "assembly": "dm3",
        "licence": "the publisher's open-archive user licence; GEO states no terms",
        "boundaries": None,
        "bin_bp": None,
        "bytes": None,
        "counted": False,
        "why_not": (
            "lane-flynode's reason, confirmed in more detail: loops and ChIP peaks only, and nothing "
            "at all for stage 16; the processed files are contact matrices"
        ),
    },
    {
        "name": "Szabo et al. 2018, Sci Adv 4:eaar8082",
        "stage": "16-18 h embryos",
        "assembly": "dm3",
        "licence": "CC BY-NC",
        "boundaries": None,
        "bin_bp": None,
        "bytes": None,
        "counted": False,
        "why_not": "no boundary list: contact matrices and a PDF supplement",
    },
    {
        "name": "Ing-Simmons et al. 2021, Nat Genet 53:487",
        "stage": "late nc14 (about 2.5-3 h) and stage 10 (about 4-5 h)",
        "assembly": "dm6 (FlyBase r6.30)",
        "licence": "CC BY 4.0",
        "boundaries": 1491,
        "bin_bp": 2000,
        "bytes": 510_265,
        "counted": False,
        "why_not": (
            "its 1,491 domains and insulation scores are not deposited as a table; the supplementary "
            "workbook holds enhancer-gene assignments and read statistics only"
        ),
    },
    {
        "name": "Cavalheiro et al. 2023, Sci Adv 9:eade1085, table S1",
        "stage": "2-3 h embryos (nuclear cycle 14)",
        "assembly": "dm6",
        "licence": "CC BY 4.0",
        "boundaries": 878,
        "bin_bp": 10000,
        "bytes": 159_931,
        "counted": False,
        "why_not": (
            "lane-flynode's reason, unchanged: 10-kb bins are coarser than many fly intergenic gaps, "
            "and the table lists 878 rows of a differential analysis against the 772 its text reports"
        ),
    },
    {
        "name": "Hug et al. 2017, Cell 169:216",
        "stage": "nuclear cycles 12-14",
        "assembly": "dm6",
        "licence": "not open access",
        "boundaries": None,
        "bin_bp": None,
        "bytes": None,
        "counted": False,
        "why_not": "lane-flynode's reason, unchanged: no openly licensed boundary list found",
    },
    {
        "name": (
            "Ramirez et al. 2018 (Kc167), Eagen et al. 2017 (Kc167), Chathoth and Zabet 2019 (BG3 and Kc167)"
        ),
        "stage": "n/a: cell lines",
        "assembly": "dm3 and dm6",
        "licence": "CC BY 4.0, the publisher's user licence, and CC BY",
        "boundaries": None,
        "bin_bp": None,
        "bytes": None,
        "counted": False,
        "why_not": (
            "lane-flynode's reason, extended: cell lines, not embryos or tissue. Chathoth and Zabet's "
            "BG3 is larval-CNS-derived and dm6, so it is the nearest cell-line analogue of the "
            "larval-CNS set counted above, which is published from tissue and preferred for that"
        ),
    },
    {
        "name": "4D Nucleome; modENCODE",
        "stage": "n/a",
        "assembly": "n/a",
        "licence": "n/a",
        "boundaries": None,
        "bin_bp": None,
        "bytes": None,
        "counted": False,
        "why_not": (
            "the 4DN portal returns no Drosophila melanogaster experiment sets, and modENCODE "
            "publishes no Hi-C or boundary release"
        ),
    },
]

SOURCE_CANDIDATES = [
    {
        "name": "Calderon et al. 2022, Science 377:eabn5800 (lane-fly's source)",
        "what": "sci-ATAC and sci-RNA over 0-20 h embryos, 14 tissues in both modalities",
        "assembly": "dm6",
        "licence": "not stated",
        "counted": True,
        "note": "the population section 9 could not decide on",
    },
    {
        "name": "Janssens et al. 2022, Nature 601:630 + Davie et al. 2018, Cell 174:982",
        "what": (
            "single-cell accessibility of larval-to-adult brain with the source's own declared map "
            "from each cell type to a cluster of the 2018 brain RNA atlas; 42 shared labels"
        ),
        "assembly": "dm6",
        "licence": "not stated on the article page, the data portal or the GEO record",
        "counted": True,
        "note": "the source this registration uses",
    },
    {
        "name": "Treese et al. 2026, PNAS 123:e2614314123, GEO GSE305940",
        "what": "10x Multiome of the optic lobe, 259 joint clusters, 76,777 consensus peaks",
        "assembly": "BDGP6 (dm6)",
        "licence": "not verified: the publisher page refused the request; not inferred",
        "counted": False,
        "why_not": (
            "its element table (36.6 MB) is affordable, but the registered endpoint needs counts per "
            "label and those are inside a 16.1 GB .rds, over the 1 GB budget; the published marker "
            "tables give marker calls, not counts, which would change the endpoint. On count alone "
            "this is the larger source and it is the first thing to try if the budget is raised"
        ),
    },
    {
        "name": "Wang et al. 2025, Cell 188:4734 (Flysta3D-v2)",
        "what": "whole-organism scATAC (40 clusters) and scRNA (45 clusters), embryo to pupa",
        "assembly": "not stated in the text read; the article page refused the request",
        "licence": "Mendeley Data states CC BY 4.0 for the processed matrices",
        "counted": False,
        "why_not": "no small peak list is published; the element coordinates need 18.7 GB of .Rds",
    },
    {
        "name": "van Nierop y Sanchez et al. 2025, Nat Commun 16:6766",
        "what": "10x Multiome of the adult testis apical tip, 25 clusters, 46,619 regions",
        "assembly": "BDGP6.32 (dm6)",
        "licence": "CC BY 4.0",
        "counted": False,
        "why_not": (
            "the best-licensed candidate, but 25 labels in one tissue; no per-label count table under "
            "1 GB was found, and the per-label element counts would be below the brain's"
        ),
    },
    {
        "name": "Bravo Gonzalez-Blas et al. 2020, Mol Syst Biol 16:e9438",
        "what": "eye-antennal disc scATAC (22 clusters) and scRNA (17 clusters), 38,179 regions",
        "assembly": "fly genome; the exact build not confirmed",
        "licence": "CC BY 4.0",
        "counted": False,
        "why_not": "17-22 labels in one disc: fewer labels and fewer regions than the brain source",
    },
    {
        "name": "Cusanovich et al. 2018, Nature 555:538",
        "what": "sci-ATAC of 2-4, 6-8 and 10-12 h embryos; 35,963 distal elements",
        "assembly": "dm3",
        "licence": "not Creative Commons; the publisher's conditions of use",
        "counted": False,
        "why_not": "no matched RNA on the same labels, and dm3 would need lifting",
    },
    {
        "name": "Reddington et al. 2020, Dev Cell 55:648",
        "what": "DNase and RNA on the same sorted lineages, 2-12 h",
        "assembly": "dm3",
        "licence": "all rights reserved",
        "counted": False,
        "why_not": "about six lineage labels: too few for a one-label-only element call; dm3",
    },
    {
        "name": "Fly Cell Atlas (Li et al. 2022, Science 375:eabk2432)",
        "what": "adult single-nucleus RNA, 15 tissues, over 250 cell types",
        "assembly": "n/a",
        "licence": "the portal states all rights reserved",
        "counted": False,
        "why_not": "RNA only: no matched accessibility is published to pair with it",
    },
]


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _member(body: bytes, member: str) -> bytes:
    kind, rest = member.split(":", 1)
    if kind == "tar":
        with tarfile.open(fileobj=io.BytesIO(gzip.decompress(body))) as tar:
            handle = tar.extractfile(rest)
            if handle is None:
                raise ValueError(f"{rest} is not in the archive")
            return handle.read()
    outer, inner = rest.split("!", 1)
    with zipfile.ZipFile(io.BytesIO(body)) as top:
        nested_bytes = top.read(outer)
    with zipfile.ZipFile(io.BytesIO(nested_bytes)) as nested:
        return nested.read(inner)


def fetch() -> None:
    fs.STORE.mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(".").free
    files: dict[str, Any] = {}
    moved = 0
    for name, spec in INPUTS.items():
        path = fs.STORE / name
        if path.exists():
            files[name] = {
                **spec,
                "bytes": path.stat().st_size,
                "sha256": _sha(path),
                "downloaded_now": False,
            }
            continue
        req = urllib.request.Request(spec["url"], headers={"User-Agent": "GenomeOS/1.0"})
        with urllib.request.urlopen(req, timeout=900) as r:
            body = r.read()
        moved += len(body)
        if moved > BUDGET_BYTES:
            sys.exit("refused: over the budget")
        path.write_bytes(_member(body, spec["member"]) if spec["member"] else body)
        files[name] = {
            **spec,
            "bytes": path.stat().st_size,
            "downloaded_bytes": len(body),
            "sha256": _sha(path),
            "downloaded_now": True,
        }
    record = {
        "fetched_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
        "disk_free_bytes_before": free,
        "files": files,
        "download_bytes": moved,
        "kept_bytes": sum(v["bytes"] for v in files.values()),
    }
    fs.FETCH_RECORD.write_text(json.dumps(record, indent=1) + "\n")
    print(json.dumps({k: v for k, v in record.items() if k != "files"}, indent=1))


# ------------------------------------------------------------------------------------ building the rows


def _elements() -> tuple[dict, list[dict[str, Any]], dict[str, int], dict[str, str], dict[str, int]]:
    chroms, tss, name_to_id = sf.read_annotation(sf.GFF3)
    dars = fs.read_dars(fs.DAR_FILE, fs.read_label_map(fs.CELL_TYPES_FILE))
    els, funnel = fs.elements(dars, chroms, tss)
    return chroms, simpler.attribute(els, chroms), funnel, name_to_id, fn.read_arm_lengths(sf.GFF3)


def _boundary_sets() -> tuple[dict[str, dict[str, list[int]]], dict[str, Any]]:
    stadler, stadler_summary = fn.read_boundaries(fn.BOUNDARY_FILE, "manual")
    dolsten, dolsten_summary = fs.read_boundary_table(fs.read_dolsten())
    embryo, embryo_summary = fs.read_boundary_table(fs.read_kaushal(fs.KAUSHAL_FILE, "exist.WT_0"))
    cns, cns_summary = fs.read_boundary_table(fs.read_kaushal(fs.KAUSHAL_CNS_FILE, "exist.wildtype"))
    return (
        {
            "stadler_manual_nc14": stadler,
            "dolsten_2025_atlas": dolsten,
            "kaushal_2022_wt_embryo": embryo,
            "kaushal_2021_wt_larval_cns": cns,
        },
        {
            "stadler_manual_nc14": stadler_summary,
            "dolsten_2025_atlas": dolsten_summary,
            "kaushal_2022_wt_embryo": embryo_summary,
            "kaushal_2021_wt_larval_cns": cns_summary,
        },
    )


def _build() -> dict[str, Any]:
    chroms, rows, funnel, name_to_id, arm = _elements()
    sets, summaries = _boundary_sets()
    attributed = {
        name: fs.attribute(rows, chroms, bounds, fs.BOUNDARY_SETS[name]["half_bin"])
        for name, bounds in sets.items()
    }
    # A summary's `span_widths_bp` is keyed by the span width as an INTEGER and JSON writes every
    # key as a string, so the summary that goes into a result file is not equal to the one computed
    # here and `score`'s guard on it could never pass, whatever the boundaries were. Normalised once,
    # at the one place the summaries are built, so what is compared is what was written.
    summaries = json.loads(json.dumps(summaries))
    return {
        "chroms": chroms,
        "rows": rows,
        "funnel": funnel,
        "name_to_id": name_to_id,
        "arm": arm,
        "sets": sets,
        "summaries": summaries,
        "main": attributed[fs.PRIMARY_SET],
        "alternatives": {k: v for k, v in attributed.items() if k != fs.PRIMARY_SET},
    }


PARAMETERS = {
    "dar_selection": fs.DAR_SELECTION,
    "never_read": list(fs.NEVER_READ),
    "labels": list(fs.LABELS),
    "boundary_sets": fs.BOUNDARY_SETS,
    "primary_boundary_set": fs.PRIMARY_SET,
    "detection_cpm": fs.DETECTION_CPM,
    "distance_bins_bp": list(fs.DISTANCE_BINS),
    "control_draws": fs.CONTROL_DRAWS,
    "rotation_draws": fs.ROTATION_DRAWS,
    "bootstrap": fs.BOOTSTRAP,
    "cluster_bp": fs.CLUSTER_BP,
    "seed": fs.SEED,
    "stop_rule": fs.STOP_RULE,
    "alpha": fs.ALPHA,
    "sweep_counts": list(fs.SWEEP_COUNTS),
    "sweep_draws": fs.SWEEP_DRAWS,
}
EXCLUSIONS = [
    "accessibility cell types for which the source declares no RNA counterpart, on both sides",
    "DAR rows outside lvl1 / vsRest / WholeBrain / Up",
    "regions the source calls up in more than one label",
    "regions within 1 kb of a transcript 5' end (kept only as the positive control)",
    "regions whose midpoint lies in a protein-coding gene body (reported beside, not the population)",
    "regions on chromosomes other than 2L, 2R, 3L, 3R, 4, X, Y",
    "the DAR table's `closestGene` column, which is the source's own answer to the question asked",
]


def _manifest(inputs: list[tuple[str, str]], partitions: Any) -> dict:
    return {
        "sources": [
            {"accession": fs.ELEMENT_SOURCE, "version": "as served on 2026-10-10"},
            {"accession": fs.RNA_SOURCE, "version": "as served on 2026-10-10"},
            {"accession": fn.BOUNDARY_SOURCE, "version": "eLife supplementary file 1, v2"},
            {"accession": sf.GENE_MODELS, "version": "Ensembl release 104"},
        ],
        "inputs": [mf.input_entry(path, partition=part) for path, part in inputs],
        "assembly": "dm6 (BDGP6.32): the regions, the boundaries, the gene models; RNA by gene name",
        "coordinates": {
            "base": 0,
            "interval": "half-open",
            "note": (
                "the DAR regions are chr:start-end as the source writes them, the GFF3 is 1-based "
                "closed, and a boundary is its span's midpoint; positions move by at most 1 bp"
            ),
        },
        "parameters": PARAMETERS,
        "exclusions": EXCLUSIONS,
        "partitions": partitions,
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


# ---------------------------------------------------------------------------------------- the inventory


def _calderon() -> dict[str, Any]:
    """lane-fly's 917 between-gene elements, and what each measured set moves on them."""
    chroms, tss, _ = sf.read_annotation(sf.GFF3)
    calls, _ = sf.tissue_calls(sf.read_markers(sf.SUPP), sf.read_cluster_labels(sf.SUPP))
    els, _ = sf.elements(calls, sf.read_peaks(sf.PEAKS), chroms, tss)
    rows = simpler.attribute(els, chroms)
    sets, _ = _boundary_sets()
    counts = {
        key: fs.between_gene_moves(rows, chroms, bounds, fs.BOUNDARY_SETS[key]["half_bin"])
        for key, bounds in sets.items()
    }
    return {"chroms": chroms, "rows": rows, "counts": counts}


def inventory() -> None:
    mf.trace_begin()
    built = _build()
    record = json.loads(fs.FETCH_RECORD.read_text())
    here = {
        key: fs.between_gene_moves(built["rows"], built["chroms"], bounds, fs.BOUNDARY_SETS[key]["half_bin"])
        for key, bounds in built["sets"].items()
    }
    calderon = _calderon()
    sweep = fs.density_sweep(calderon["rows"], calderon["chroms"], built["arm"])
    ceiling = fs.tailored_ceiling(calderon["rows"], calderon["chroms"])
    boundaries = []
    for cand in BOUNDARY_CANDIDATES:
        row = dict(cand)
        key = cand.get("key")
        if key:
            row["on_calderon_917_between_gene_elements"] = calderon["counts"][key]
            row["on_this_lanes_brain_elements"] = here[key]
        boundaries.append(row)
    sources = []
    for cand in SOURCE_CANDIDATES:
        row = dict(cand)
        got = None
        if cand["name"].startswith("Calderon"):
            got = calderon["counts"][fs.PRIMARY_SET]
        elif cand["name"].startswith("Janssens"):
            got = here[fs.PRIMARY_SET]
        if got is not None:
            row["between_gene_elements"] = got["between_gene_elements"]
            row["between_gene_elements_that_change_gene"] = got["change_gene"]
            row["counted_under"] = fs.PRIMARY_SET
        sources.append(row)
    out = {
        "question": (
            "a 6-16 h fly boundary set, or a source with more between-gene elements: what exists, "
            "under what licence, at what stage, and how many between-gene elements each would put in "
            "the population and move -- counted from positions only, before any outcome was read"
        ),
        "route_a_boundary_sets": boundaries,
        "route_a_is_closed_by_arithmetic": {
            "why": (
                "under the registered rule a between-gene element changes gene only when its nearest "
                "5' end lies across a boundary AND the domain it leaves holds another 5' end. More "
                "boundaries meet the first condition more often and the second less often, so what "
                "the rule can move is a curve in the boundary count and not a growing function"
            ),
            "density_sweep_on_calderon_917": sweep,
            "tailored_ceiling_on_calderon_917": ceiling,
            "reading": (
                "On Calderon et al. 2022's 917 between-gene elements the number of elements the "
                f"boundary rule names a different gene for peaks at "
                f"{sweep['best_single_draw_change_gene']} across "
                f"{len(fs.SWEEP_COUNTS)} densities from {fs.SWEEP_COUNTS[0]:,} to "
                f"{fs.SWEEP_COUNTS[-1]:,} boundaries, {fs.SWEEP_DRAWS} uniform placements each, and "
                f"never reaches the stop rule's {fs.STOP_RULE}. The four measured sets give "
                + ", ".join(f"{v['change_gene']} ({k})" for k, v in sorted(calderon["counts"].items()))
                + f". A set with one boundary placed inside each element's own gap reaches "
                f"{ceiling['change_gene']}, so the limit is where the boundaries lie and not the rule. "
                "Measured fly boundaries sit mostly at active promoters, which is the case where the "
                "rule does NOT move, so a uniform placement is if anything generous. No measured, "
                "density-plausible fly boundary set decides the between-gene question on this "
                "population, whatever its stage; the source has to change, and it does here"
            ),
        },
        "route_b_sources": sources,
        "counts_on_this_lanes_brain_elements": here,
        "fetched": {
            name: {
                "bytes": entry["bytes"],
                "sha256": entry["sha256"],
                "licence": entry["licence"],
                "what": entry["what"],
            }
            for name, entry in record["files"].items()
        },
        "fetched_bytes": record["download_bytes"],
        "kept_bytes": record["kept_bytes"],
        "scouting_bytes": SCOUTING_BYTES,
        "lane_bytes_moved": record["kept_bytes"] + SCOUTING_BYTES,
        "budget_bytes": BUDGET_BYTES,
        "disk_free_bytes_before": record["disk_free_bytes_before"],
        "boundary_lists": built["summaries"],
        "population_funnel": built["funnel"],
        "verdict": fs.verdict(sources),
        "stop_list": (
            "docs/ROADMAP.md section 5's 'No more biosamples' clause is about the human reader's 13 "
            "biosamples (READER_TERMS). None is added or named here. A Drosophila Hi-C boundary list "
            "and a Drosophila accessibility or RNA source are not reader biosamples"
        ),
        "money": "none; no model request, no key read",
    }
    out["result_manifest"] = _manifest(
        [
            (str(fs.FETCH_RECORD), "the fetch record"),
            (str(fs.CELL_TYPES_FILE), "the `RNAcluster` column only"),
            (str(fs.DAR_FILE), "the selected rows' region and cell type only"),
            (str(fn.BOUNDARY_FILE), "Table S4 (manual boundary positions)"),
            (str(fs.DOLSTEN_FILE), "sheet `boundaries`: chrom, start, end"),
            (str(fs.KAUSHAL_FILE), "chr, position and `exist.WT_0`"),
            (str(fs.KAUSHAL_CNS_FILE), "chr, position and `exist.wildtype`"),
            (str(sf.PEAKS), "Calderon's peak set, for the density sweep's population"),
            (str(sf.SUPP), "Calderon's Table S2 and Table S4, for the density sweep's population"),
            (str(sf.GFF3), "gene and transcript positions and strands, and arm lengths"),
        ],
        "n/a: an inventory has no evaluation partition",
    )
    print(save_result(INVENTORY, out))
    print(json.dumps(out["verdict"], indent=1))
    print(out["route_a_is_closed_by_arithmetic"]["reading"])


# ------------------------------------------------------------------------------------- the registration


def register() -> None:
    mf.trace_begin()
    built = _build()
    cen = fs.census(built["main"], built["alternatives"], built["funnel"])
    record = json.loads(fs.FETCH_RECORD.read_text())
    out = {
        "written_before_any_rna_value_was_read": True,
        "registered_inputs": {
            name: {"path": str(fs.STORE / name), "sha256": record["files"][name]["sha256"]} for name in INPUTS
        }
        | {
            "boundaries_stadler": {
                "path": str(fn.BOUNDARY_FILE),
                "sha256": simpler.sha256(fn.BOUNDARY_FILE),
            },
            "flynode_registration": {
                "path": FLYNODE_REGISTRATION_FILE,
                "sha256": simpler.sha256(Path(FLYNODE_REGISTRATION_FILE)),
            },
            "gff3": {"path": str(sf.GFF3), "sha256": simpler.sha256(sf.GFF3)},
        },
        "boundary_lists": built["summaries"],
        "arm_lengths": built["arm"],
        "registration": fs.REGISTRATION,
        "census": cen,
        "alphagenome_requests": 0,
        "money": "none",
    }
    out["result_manifest"] = _manifest(
        [
            (str(fs.FETCH_RECORD), "each input's registered sha256"),
            (str(fs.CELL_TYPES_FILE), "the `RNAcluster` column only"),
            (str(fs.DAR_FILE), "the selected rows' region and cell type only"),
            (str(fn.BOUNDARY_FILE), "Table S4 (manual boundary positions)"),
            (str(fs.DOLSTEN_FILE), "sheet `boundaries` (beside)"),
            (str(fs.KAUSHAL_FILE), "`exist.WT_0` rows (beside)"),
            (str(fs.KAUSHAL_CNS_FILE), "`exist.wildtype` rows (beside)"),
            (str(sf.GFF3), "gene and transcript positions and strands, and arm lengths"),
            (str(fs.RNA_FILE), "hashed only; no value is read here"),
            (FLYNODE_REGISTRATION_FILE, "the registration this one carries its wording from"),
        ],
        {"population": "fixed here, with the census in view and before any RNA value was read"},
    )
    print(save_result(REGISTRATION, out))
    print(json.dumps(cen["between_genes"], indent=1))


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
    if reg["registration"] != json.loads(json.dumps(fs.REGISTRATION)):
        sys.exit("refused: simpler_flysource.REGISTRATION differs from the committed registration")
    for key, entry in reg["registered_inputs"].items():
        if simpler.sha256(Path(entry["path"])) != entry["sha256"]:
            sys.exit(f"refused: {key} differs from the registered input")
    built = _build()
    if fs.census(built["main"], built["alternatives"], built["funnel"]) != reg["census"]:
        sys.exit("refused: the census does not recompute to the registered one")
    if built["summaries"] != reg["boundary_lists"]:
        sys.exit("refused: a boundary list differs from the registered one")
    expr, rna_summary = fs.read_expression(fs.RNA_FILE, built["name_to_id"])
    res = fs.score(
        built["chroms"],
        built["main"],
        built["alternatives"],
        built["sets"][fs.PRIMARY_SET],
        built["arm"],
        expr,
        built["funnel"],
    )
    reg_sha = subprocess.run(
        ["git", "log", "-1", "--format=%h", "--", REGISTRATION_FILE],
        cwd=ROOT,
        capture_output=True,
        text=True,
    ).stdout.strip()
    out = {
        "registration_commit": reg_sha,
        "registered_inputs": reg["registered_inputs"],
        "rna_summary": rna_summary,
        **res,
        "reading_is_the_registered_wording": (
            "the reading is REGISTRATION['outcome_wording'] with its numbers filled in and no other "
            "word changed"
        ),
        "alphagenome_requests": 0,
        "money": "none",
    }
    out["result_manifest"] = _manifest(
        [
            (REGISTRATION_FILE, "the committed registration"),
            (str(fs.CELL_TYPES_FILE), "the `RNAcluster` column only"),
            (str(fs.DAR_FILE), "the selected rows' region and cell type only"),
            (str(fn.BOUNDARY_FILE), "Table S4 (manual boundary positions)"),
            (str(fs.DOLSTEN_FILE), "sheet `boundaries` (beside)"),
            (str(fs.KAUSHAL_FILE), "`exist.WT_0` rows (beside)"),
            (str(fs.KAUSHAL_CNS_FILE), "`exist.wildtype` rows (beside)"),
            (str(sf.GFF3), "gene and transcript positions and strands, and arm lengths"),
            (str(fs.RNA_FILE), "the endpoint: counts summed by label"),
        ],
        {
            "population": "pinned by the committed registration's census",
            "positive_control": "promoter elements",
        },
    )
    print(save_result(RESULT, out))
    print(out["outcome"])
    print(out.get("reading", ""))


if __name__ == "__main__":
    os.chdir(ROOT)
    {"fetch": fetch, "inventory": inventory, "register": register, "score": score}[sys.argv[1]]()
