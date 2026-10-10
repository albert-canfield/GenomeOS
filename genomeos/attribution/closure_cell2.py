# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""What follows M3: the public element-level fitness screens the cell-level chain could run on next.

M3 (`attribution/closure_cell.py`, result `80bf3db`) ran the chain element -> gene -> what the cell
cannot live without once, on ENCODE's genome-wide CRISPRi proliferation screen in K562, and detected
no difference over the 35 elements whose attributed gene is one K562 depends on and whose nearest
coding gene is not. Area I's Next list asks for what follows it: a larger population or another
screen, registered before it is run.

This module is the SEARCH, and it is deliberately only the search. It holds every public screen the
chain could be run on, and for each one it counts -- from POSITIONS ALONE, before any guide count or
phenotype value is read -- how many of the project's compiled attributed elements the screen's
library covers, and how the four arms of M3's rule would divide them in that screen's own cell. The
gate the count is read against is fixed here (`GATE`) rather than after the counting, because the
question "is there a bigger or better-powered population than M3's 35" has to be answerable with a
no.

Three things the chain needs in ONE cell, and the search is narrow because of the third:

* the project's own claim -- element e drives gene G on that cell's own RNA track. The sweep scored
  four cells (`predicted_coding_by_cell`): K562, HepG2, GM12878 and IMR-90. A screen in any other
  cell cannot carry the chain's first statement, however good the screen is;
* what the cell cannot live without -- a DepMap 24Q4 Chronos column for it. Of those four, K562
  (ACH-000551) and HepG2 (ACH-000739) have one; GM12878 and IMR-90 are not in the release;
* an element-level fitness or growth readout in that cell, public and not M3's own screen.

So the reachable cells are K562 and HepG2 and nothing else, and every other screen in the table
below is recorded with the reason it is out rather than left unsaid.

WHAT THIS CANNOT ESTABLISH. A coverage count is not a measurement: an element the library covers may
still be unmeasured when the counts are read (M3's `MIN_CONTROL_COUNT` gate). The arm sizes here are
therefore an upper bound on the measured arm sizes, which makes them the right thing to gate on --
a population that is too small before any depth is lost is too small afterwards.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import urllib.request
from pathlib import Path
from typing import Any

from genomeos.attribution import closure_cell as cc
from genomeos.therapeutics import depmap

#: Where the design files and the distilled dependency columns are kept; never committed.
STORE = Path("data/knowledge/closure_cell2")

ENCODE_LICENCE = cc.ENCODE_LICENCE
ENCODE_FILE = "https://www.encodeproject.org/files/{a}/@@download/{a}.{ext}"

#: The cells the chain can reach, with the DepMap model whose Chronos column is the dependency.
#: A cell is here only if the sweep scored it AND DepMap 24Q4 carries it.
CELLS = {
    "K562": "ACH-000551",
    "HepG2": "ACH-000739",
}
#: The sweep's other two cells, kept in the record so the absence is stated and not implied.
CELLS_WITHOUT_A_DEPENDENCY_COLUMN = {
    "GM12878": "no model in DepMap 24Q4 Public (a lymphoblastoid line, not a cancer model)",
    "IMR-90": "no model in DepMap 24Q4 Public (a primary fibroblast line, not a cancer model)",
}

#: Every library's DESIGN file: positions, and no measurement of any kind. `rows` says what one row
#: is, because the libraries are not all alike -- a 3 bp guide target in one, a paired-guide deletion
#: interval in another -- and an overlap count means something different in each.
LIBRARIES: dict[str, dict[str, Any]] = {
    "ENCFF232FVA": {
        "ext": "bed.gz",
        "format": "bed6",
        "md5": "6d59b65691282a0abbc54dfad01fa616",
        "rows": "one 3 bp guide target interval per guide; non-targeting guides have no position "
        "and are absent from the file",
        "reference": "ENCSR515QTA (guide locations for the validation CERES screen), GRCh38",
        "guide_type_column": None,
    },
    "ENCFF745ZBH": {
        "ext": "bed.gz",
        "format": "bed3",
        "md5": "aa075503d0a6832b213081d6dcc90f90",
        "rows": "one paired-guide (pgRNA) interval per designed deletion, about 19 bp as published",
        "reference": "ENCSR488FGH (guide locations, elements selected by histone modifications), GRCh38",
        "guide_type_column": None,
    },
    "ENCFF526YZH": {
        "ext": "bed.gz",
        "format": "bed3",
        "md5": "e9cb46ea964276248505be2a47dda896",
        "rows": "one paired-guide (pgRNA) deletion interval per designed deletion, about 1 kb",
        "reference": "ENCSR327RRR (guide locations, MYC and MYB loci), GRCh38",
        "guide_type_column": None,
    },
    "ENCFF135FSS": {
        "ext": "csv",
        "format": "csv",
        "md5": "9da453cd01cbda201b3ac476385fc0df",
        "rows": "one 3 bp guide target interval per guide, with a guideType column; only "
        "`targeting` rows are counted and the negative controls are counted apart",
        "reference": "ENCSR251MMR (guide locations, cCREs around essential genes), GRCh38",
        "guide_type_column": "guideType",
    },
    "ENCFF865MGL": {
        "ext": "csv",
        "format": "csv",
        "md5": "3f9a4e999934e095d4bed70ae44154d3",
        "rows": "one 3 bp guide target interval per guide, with a guideType column",
        "reference": "ENCSR193VJY (guide locations, ~7,500 cCREs around ~450 essential genes), GRCh38",
        "guide_type_column": "guideType",
    },
    "ENCFF271CBI": {
        "ext": "csv",
        "format": "csv",
        "md5": "79b839eae7d79ac25b0891972ccb0dab",
        "rows": "one 3 bp guide target interval per guide, with a guideType column",
        "reference": "ENCSR027FIB (guide locations, CTCF motifs in loop anchors), GRCh38",
        "guide_type_column": "guideType",
    },
    "ENCFF964OOO": {
        "ext": "csv",
        "format": "csv",
        "md5": "22eb43a1353859e6a78dd7b4bfd814fe",
        "rows": "one 3 bp guide target interval per guide, with a guideType column",
        "reference": "ENCSR145WNY (guide locations, GATA1, MYB and ZMYND8 loci), GRCh38",
        "guide_type_column": "guideType",
    },
}

#: Every public screen considered, whether or not the chain can run on it. `cell` is the cell the
#: screen was done in; `reachable` says whether the project can make its own claim in that cell.
#: `readout_files` are NOT fetched here -- they are the outcome, and this search does not open them.
CANDIDATES: tuple[dict[str, Any], ...] = (
    {
        "key": "wgceres_validation_k562",
        "cell": "K562",
        "reachable": True,
        "assay": "proliferation CRISPR screen (CRISPRi, dCas9-KRAB), validation sublibrary",
        "screen": "ENCSR995JLE",
        "control": "ENCSR501FEU",
        "lab": "Tim Reddy, Duke / Charles Gersbach, Duke (ENCODE4 wgCERES)",
        "publication": "Klann et al., Genome-wide annotation of gene regulatory elements linked to "
        "cell fitness",
        "library": "ENCFF232FVA",
        "readout_files": {
            "screen": ["ENCFF237LMB", "ENCFF391APL", "ENCFF656EAR", "ENCFF663WXT"],
            "control": ["ENCFF365SEP", "ENCFF413BFN", "ENCFF112NVI", "ENCFF798BFU"],
            "bytes": 311_641_992,
            "output_type": "guide quantifications, GRCh38, four biological replicates a side",
        },
        "independence": "the validation sublibrary's elements were SELECTED by the discovery "
        "screen's own outcome -- the screen M3 read -- so a population drawn from it is not "
        "independent of M3's measurement, and any reading from it would have to say so",
    },
    {
        "key": "wgceres_validation_ociaml2",
        "cell": "OCI-AML2",
        "reachable": False,
        "not_reachable_because": "the sweep never scored OCI-AML2 (its cells are K562, HepG2, "
        "GM12878 and IMR-90), so the chain's first statement -- the project's own claim on that "
        "cell's track -- cannot be made there; DepMap does carry the line (ACH-000113)",
        "assay": "proliferation CRISPR screen (CRISPRi, dCas9-KRAB), validation sublibrary",
        "screen": "ENCSR689KVW",
        "control": "ENCSR770ENU",
        "lab": "Tim Reddy, Duke / Charles Gersbach, Duke (ENCODE4 wgCERES)",
        "library": "ENCFF232FVA",
        "readout_files": {
            "screen": ["ENCFF549LWY", "ENCFF924VSR", "ENCFF336BMY", "ENCFF963RQE"],
            "control": ["ENCFF693EGA", "ENCFF518YXD", "ENCFF933LBZ", "ENCFF280QIY"],
            "bytes": 312_124_444,
            "output_type": "guide quantifications, GRCh38, four biological replicates a side",
        },
    },
    {
        "key": "ren_histone_selected_k562",
        "cell": "K562",
        "reachable": True,
        "assay": "proliferation CRISPR screen (CRISPRi, paired guides deleting the element)",
        "screen": "ENCSR918PMD",
        "control": "ENCSR469EDZ (the plasmid pool, sequenced as the library's own starting point)",
        "lab": "Bing Ren, UCSD",
        "library": "ENCFF745ZBH",
        "readout_files": {
            "screen": ["ENCFF912HFE"],
            "bytes": 18_147_161,
            "output_type": "guide quantifications, GRCh38, one series-level file",
        },
        "independence": "a different library and a different lab from M3's screen, in the same cell",
    },
    {
        "key": "ren_histone_selected_hct116",
        "cell": "HCT116",
        "reachable": False,
        "not_reachable_because": "the sweep never scored HCT116, so the project's own claim cannot "
        "be made in that cell; DepMap does carry the line (ACH-000971)",
        "assay": "proliferation CRISPR screen (CRISPRi, paired guides deleting the element)",
        "screen": "ENCSR924MQG",
        "lab": "Bing Ren, UCSD",
        "library": "ENCFF745ZBH",
        "readout_files": {
            "screen": ["ENCFF540TVK"],
            "bytes": 20_052_028,
            "output_type": "guide quantifications, GRCh38, one series-level file",
        },
    },
    {
        "key": "ren_myc_myb_k562",
        "cell": "K562",
        "reachable": True,
        "assay": "growth-based CRISPR screen (CRISPRi, paired guides deleting the element), MYC and MYB loci",
        "screen": "ENCSR411TAL",
        "lab": "Bing Ren, UCSD",
        "library": "ENCFF526YZH",
        "readout_files": {
            "screen": ["ENCFF925SXM"],
            "bytes": 4_011_361,
            "output_type": "guide quantifications, GRCh38, one series-level file",
        },
        "independence": "a different library and a different lab from M3's screen, in the same cell",
    },
    {
        "key": "ren_myc_myb_hepg2",
        "cell": "HepG2",
        "reachable": True,
        "assay": "growth-based CRISPR screen (CRISPRi, paired guides deleting the element), MYC and MYB loci",
        "screen": "ENCSR319WEU",
        "lab": "Bing Ren, UCSD",
        "library": "ENCFF526YZH",
        "readout_files": {
            "screen": ["ENCFF375BGA"],
            "bytes": 4_023_030,
            "output_type": "guide quantifications, GRCh38, one series-level file",
        },
        "independence": "the only candidate in a SECOND cell the chain can reach: the sweep scored "
        "HepG2's own track and DepMap 24Q4 carries HepG2 (ACH-000739)",
    },
    {
        "key": "greenleaf_sparse_ccre_k562",
        "cell": "K562",
        "reachable": True,
        "assay": "proliferation CRISPR screen (CRISPRi, several dCas9 effectors), sparse tiling of "
        "~7,500 cCREs around ~450 essential genes",
        "screen": "ENCSR398UPE, ENCSR237VMX, ENCSR546ILA, ENCSR377LWO, ENCSR509FZL",
        "lab": "Will Greenleaf, Stanford",
        "library": "ENCFF865MGL",
        "readout_files": {
            "screen": ["ENCFF432FPQ", "ENCFF875ELG", "ENCFF227RPU", "ENCFF882CKL"],
            "bytes": 125_383_605,
            "output_type": "guide quantifications, GRCh38, two biological replicates per effector",
        },
        "independence": "a different library and a different lab from M3's screen, in the same cell",
    },
    {
        "key": "greenleaf_ccre_tiling_k562",
        "cell": "K562",
        "reachable": True,
        "assay": "proliferation CRISPR screen (CRISPRi and dCas9), peak tiling of a subset of cCREs "
        "around essential genes",
        "screen": "ENCSR839FIS, ENCSR013IWY, ENCSR007UZQ, ENCSR039ZOS, ENCSR294JZL",
        "lab": "Will Greenleaf, Stanford",
        "library": "ENCFF135FSS",
        "readout_files": {
            "screen": ["ENCFF508SFM", "ENCFF101ZBD"],
            "bytes": 4_291_272,
            "output_type": "guide quantifications, GRCh38, two biological replicates",
        },
        "independence": "a different library and a different lab from M3's screen, in the same cell",
    },
    {
        "key": "greenleaf_ctcf_anchors_k562",
        "cell": "K562",
        "reachable": True,
        "assay": "proliferation CRISPR screen (CRISPRi, dCas9 and Cas9), CTCF motifs in loop anchors",
        "screen": "ENCSR371TIR, ENCSR223KWN, ENCSR816GWF, ENCSR842CIE",
        "lab": "Will Greenleaf, Stanford",
        "library": "ENCFF271CBI",
        "readout_files": {
            "screen": ["ENCFF409HZI", "ENCFF411MWK"],
            "bytes": 5_273_870,
            "output_type": "guide quantifications, GRCh38, two biological replicates",
        },
        "independence": "a different library and a different lab from M3's screen, in the same cell",
    },
    {
        "key": "greenleaf_three_loci_k562",
        "cell": "K562",
        "reachable": True,
        "assay": "proliferation CRISPR screen (CRISPRi, CRISPRa, Cas9 and dCas9), full tiling of "
        "the GATA1, MYB and ZMYND8 loci",
        "screen": "ENCSR356RUW, ENCSR403LMA, ENCSR138FGY, ENCSR397CXJ",
        "lab": "Will Greenleaf, Stanford",
        "library": "ENCFF964OOO",
        "readout_files": {
            "screen": ["ENCFF217TZA", "ENCFF124SYO"],
            "bytes": 10_554_135,
            "output_type": "guide quantifications, GRCh38, two biological replicates",
        },
        "independence": "a different library and a different lab from M3's screen, in the same cell",
    },
)

#: The other cells the MYC/MYB growth screen was run in. Each has a DepMap column and none has a
#: track in the sweep, so the chain cannot be run in any of them; listed so the search is complete.
OTHER_CELLS_OF_THE_MYC_MYB_SCREEN = {
    "A549": "ENCSR917YQW",
    "DU 145": "ENCSR112MYY",
    "HCT116": "ENCSR675XBE",
    "MCF-7": "ENCSR147TRG",
    "MDA-MB-231": "ENCSR371XHF",
    "NCI-H460": "ENCSR401HFH",
    "PC-3": "ENCSR703YMV",
    "SW620": "ENCSR172TKR",
}

# ---- the gate, fixed before anything was counted --------------------------------------------------

#: M3's own numbers, which a successor has to beat to be worth running: 35 elements in arm A, and a
#: separation of 0.6217 as the smallest the rank test could have detected.
M3_ARM_A = 35
M3_DETECTABLE = 0.6217
M3_RESULT = "data/results/closure_cell_result.json"

ALPHA = 0.05
POWER = 0.8
Z_ALPHA = 1.6449
Z_BETA = 0.8416

GATE = {
    "arm_a_at_least": M3_ARM_A,
    "smallest_detectable_separation_at_most": M3_DETECTABLE,
    "both": "a candidate qualifies for a registered run only if arm A (the element's attributed "
    "gene is one the cell depends on and its nearest coding gene is not) covers at least as many "
    "elements as M3's 35 AND the rank test against that candidate's own floor arm could detect a "
    "separation at least as small as M3's 0.6217; a candidate that is bigger but worse powered, or "
    "better powered on fewer elements, does not qualify",
    "why_both": "area I's Next list asks for a larger population or a better-decided one; a "
    "successor that is neither repeats M3 with more words",
    "if_none_qualifies": "the lane stops here and this inventory is the result: no second "
    "cell-level run is registered, and what would change the answer is named",
}


def smallest_detectable(n1: int, n2: int) -> float | None:
    """The smallest separation a rank test on these two arm sizes can detect, as a probability.

    M3's formula, unchanged (`scripts/closure_cell_register.py`): the quantity is P(a random element
    of the first arm sits below a random element of the second), 0.5 under the null, with standard
    error sqrt((n1 + n2 + 1) / (12 n1 n2)). Written here rather than imported from the script so
    this module does not depend on a script, and pinned by a test that compares it against M3's
    registered numbers.
    """
    if n1 <= 0 or n2 <= 0:
        return None
    se = math.sqrt((n1 + n2 + 1) / (12 * n1 * n2))
    return round(0.5 + (Z_ALPHA + Z_BETA) * se, 4)


def qualifies(arm_a: int, detectable: float | None) -> bool:
    """Whether a candidate's counted arms meet `GATE`. Read once, over numbers from positions only."""
    if detectable is None:
        return False
    return arm_a >= M3_ARM_A and detectable <= M3_DETECTABLE


# ---- the design files ---------------------------------------------------------------------------


def design_path(accession: str, store: Path = STORE) -> Path:
    return store / f"{accession}.{LIBRARIES[accession]['ext']}"


def fetch_design(accession: str, store: Path = STORE, net: bool = True) -> dict[str, Any]:
    """One library design file on disk, with its md5 checked against what ENCODE published.

    The file carries guide or deletion POSITIONS and no count, which is why it may be read before
    the registration exists. The md5 is compared against the portal's, so a truncated download is a
    failure here and not a quiet change in a coverage count.
    """
    spec = LIBRARIES[accession]
    path = design_path(accession, store)
    if not path.exists():
        if not net:
            raise FileNotFoundError(f"no design file at {path} and network disabled")
        path.parent.mkdir(parents=True, exist_ok=True)
        url = ENCODE_FILE.format(a=accession, ext=spec["ext"])
        request = urllib.request.Request(url, headers={"User-Agent": "GenomeOS/0.1 (closure-cell2)"})
        with urllib.request.urlopen(request, timeout=900) as response:  # noqa: S310
            path.write_bytes(response.read())
    raw = path.read_bytes()
    md5 = hashlib.md5(raw).hexdigest()  # noqa: S324 - compared against ENCODE's published md5
    if md5 != spec["md5"]:
        raise ValueError(f"{accession}: md5 {md5} is not ENCODE's {spec['md5']}")
    return {
        "accession": accession,
        "path": str(path),
        "bytes": len(raw),
        "md5": md5,
        "md5_published_by_encode": spec["md5"],
        "sha256": hashlib.sha256(raw).hexdigest(),
        "format": spec["format"],
        "rows": spec["rows"],
        "reference": spec["reference"],
        "licence": ENCODE_LICENCE,
    }


def read_design(
    accession: str, store: Path = STORE
) -> tuple[dict[str, list[tuple[int, int, None]]], dict[str, int]]:
    """A library's intervals by chromosome, and the row census of the file they came from.

    `counted` is every row that enters the coverage count: all of them for a BED, whose
    non-targeting guides have no position, and the `targeting` rows alone for a CSV that names a
    guide type.
    """
    spec = LIBRARIES[accession]
    path = design_path(accession, store)
    out: dict[str, list[tuple[int, int, None]]] = {}
    census = {"rows": 0, "counted": 0, "not_targeting": 0, "malformed": 0}
    raw = path.read_bytes()
    text = gzip.decompress(raw).decode("utf-8", "replace") if spec["ext"].endswith("gz") else raw.decode()
    lines = text.splitlines()
    if spec["guide_type_column"]:
        for row in csv.DictReader(lines):
            census["rows"] += 1
            chrom = (row.get("chrom") or "").strip()
            if not chrom.startswith("chr"):
                census["malformed"] += 1
                continue
            if (row.get(spec["guide_type_column"]) or "") != cc.TARGETING:
                census["not_targeting"] += 1
                continue
            try:
                start, end = int(row["chromStart"]), int(row["chromEnd"])
            except (KeyError, TypeError, ValueError):
                census["malformed"] += 1
                continue
            census["counted"] += 1
            out.setdefault(chrom, []).append((start, end, None))
    else:
        for line in lines:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            census["rows"] += 1
            f = line.split("\t")
            if len(f) < 3 or not f[0].startswith("chr"):
                census["malformed"] += 1
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                census["malformed"] += 1
                continue
            census["counted"] += 1
            out.setdefault(f[0], []).append((start, end, None))
    return out, census


def median_width(by_chrom: dict[str, list[tuple[int, int, None]]]) -> int | None:
    """The median interval width of a design, so a reader can see what one of its rows covers."""
    widths = sorted(e - s for rows in by_chrom.values() for s, e, _ in rows)
    if not widths:
        return None
    n = len(widths)
    return widths[n // 2] if n % 2 else (widths[n // 2 - 1] + widths[n // 2]) // 2


# ---- the dependency column, for any one DepMap model --------------------------------------------


def gene_effect_row(model: str, genes: set[str], net: bool = True, store: Path = STORE) -> dict[str, Any]:
    """The Chronos gene effect of `genes` in ONE DepMap model, from the cache or streamed once.

    `closure_cell.k562_gene_effect` does this for K562 alone and this is the same procedure with the
    model as an argument, because the second cell the chain can reach is HepG2. One row is kept and
    nothing else; the digest recorded is over the bytes that were streamed, so the record names the
    file that was read.
    """
    cache = store / f"depmap_{model}_gene_effect.json.gz"
    cached = _read_gz(cache)
    wanted = {g.upper() for g in genes}
    if cached is not None and cached.get("release") == depmap.RELEASE and cached.get("model") == model:
        if wanted <= set(cached.get("requested") or ()):
            return {**cached, "from_cache": True}
        wanted |= set(cached.get("requested") or ())
    if not net:
        raise depmap.DepMapUnavailableError(
            f"no cached gene effect for {model} covering {len(wanted)} genes and network disabled"
        )
    file_id, name, what = depmap.FILES["gene_effect"]
    url = depmap.FIGSHARE_FILE.format(file_id)
    digest = hashlib.sha256()
    total = 0
    header: list[str] | None = None
    values: dict[str, float] = {}
    found = False
    request = urllib.request.Request(url, headers=depmap.UA)
    with urllib.request.urlopen(request, timeout=900) as response:  # noqa: S310
        pending = b""
        for chunk in iter(lambda: response.read(1 << 20), b""):
            digest.update(chunk)
            total += len(chunk)
            if found:
                continue  # the row is in hand; the rest is streamed only to finish the digest
            pending += chunk
            lines = pending.split(b"\n")
            pending = lines.pop()
            for raw in lines:
                row = next(csv.reader(io.StringIO(raw.decode("utf-8", "replace"))), [])
                if header is None:
                    header = [depmap.symbol(c) for c in row]
                    continue
                if not row or row[0] != model:
                    continue
                for i, cell in enumerate(row[1:], start=1):
                    if i < len(header) and header[i] in wanted and cell not in ("", "NA"):
                        values[header[i]] = round(float(cell), 4)
                found = True
                break
    payload = {
        "dataset": "gene_effect",
        "file": name,
        "figshare_file": file_id,
        "what": what,
        "model": model,
        "release": depmap.RELEASE,
        "doi": depmap.DOI,
        "licence": depmap.LICENCE,
        "citation": depmap.CITATION,
        "streamed_bytes": total,
        "streamed_sha256": digest.hexdigest(),
        "requested": sorted(wanted),
        "genes": sorted(values),
        "missing": sorted(wanted - set(values)),
        "values": values,
        "row_found": found,
    }
    _write_gz(cache, payload)
    return {**payload, "from_cache": False}


# ---- the counting -------------------------------------------------------------------------------


def effect_in(element: dict[str, Any], cell: str) -> float | None:
    """The element's predicted effect on its target gene on ONE cell's own track, or None.

    `closure_cell.element_effect` is this for K562, with the same sign convention: a fall in the
    gene on deleting the element means the element holds the gene up in that cell, so the effect is
    POSITIVE.
    """
    v = (element.get("predicted_coding_by_cell") or {}).get(cell)
    return None if v is None else -float(v)


def arms_over(
    rows: list[dict[str, Any]],
    intervals: dict[str, list[tuple[int, int, None]]],
    effects: dict[str, float],
    cell: str,
) -> dict[str, Any]:
    """How M3's four arms and its positive control divide the elements one library covers.

    `rows` are the compiled attributed elements, each already carrying its nearest coding gene and
    its per-cell effects; `intervals` is one library's design. Nothing but positions and the
    dependency column is read, and the dependency is the chain's premise, not its outcome.
    """
    sizes: dict[str, int] = {}
    unassigned: dict[str, int] = {}
    counts = {
        "elements_with_a_guide": 0,
        "elements_with_min_guides": 0,
        "no_value_on_the_cells_own_track": 0,
        "promoter_elements": 0,
        "distal_elements": 0,
    }
    index = {c: (rs, cc.bin_index(rs)) for c, rs in intervals.items()}
    for row in rows:
        held = index.get(row["chrom"])
        if held is None:
            continue
        guides = len(cc.overlapping(held[0], held[1], row["start"], row["end"]))
        if guides == 0:
            continue
        counts["elements_with_a_guide"] += 1
        if guides < cc.MIN_GUIDES:
            continue
        counts["elements_with_min_guides"] += 1
        effect = row["effects"].get(cell)
        if effect is None:
            counts["no_value_on_the_cells_own_track"] += 1
            continue
        target = cc.dependency(effects, row["gene"])
        neighbour = cc.dependency(effects, row["nearest_coding_gene"])
        if row["nearest_coding_distance"] <= cc.PROMOTER_WINDOW:
            counts["promoter_elements"] += 1
            if target == "essential":
                sizes[cc.POSITIVE_CONTROL] = sizes.get(cc.POSITIVE_CONTROL, 0) + 1
            else:
                key = "promoter_element_whose_gene_the_cell_does_not_depend_on"
                unassigned[key] = unassigned.get(key, 0) + 1
            continue
        counts["distal_elements"] += 1
        arm = cc.arm_of(target, neighbour, effect)
        if arm is None:
            if abs(effect) < cc.MIN_EFFECT:
                key = "prediction_under_the_magnitude_bar"
            elif target not in ("essential", "not") or neighbour not in ("essential", "not"):
                key = "a_dependency_call_is_borderline_or_unmeasured"
            else:
                key = "negative_direction_with_a_gene_the_cell_does_not_depend_on"
            unassigned[key] = unassigned.get(key, 0) + 1
        else:
            sizes[arm] = sizes.get(arm, 0) + 1
    arm_a = sizes.get("attributed_essential_neighbourhood_not", 0)
    floor = sizes.get("attributed_not_neighbourhood_not", 0)
    detectable = smallest_detectable(arm_a, floor)
    return {
        "counts": counts,
        "sizes": sizes,
        "in_no_arm": unassigned,
        "arm_a": arm_a,
        "floor_arm": floor,
        "smallest_detectable_separation": detectable,
        "qualifies": qualifies(arm_a, detectable),
    }


def _read_gz(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    try:
        with gzip.open(path, "rt") as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError):
        return None


def _write_gz(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        json.dump(payload, fh, separators=(",", ":"))
