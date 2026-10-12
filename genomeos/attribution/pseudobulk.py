# SPDX-License-Identifier: AGPL-3.0-or-later
"""Does pseudobulked single-nucleus RNA reproduce bulk RNA-seq on M2h's endpoint?

WHY THIS LANE EXISTS. M2h (`ac5580b`) found the compiled rules fit measured RNA better when scaled
by one fitted constant -- held-out Delta MSE +0.04339 [98.33% CI +0.0379, +0.04894], scale +0.2146,
about 2% of variance -- and M2i (`f6e71c6`) replicated it with nothing fitted on RNA fresh to this
project, +0.03727 [+0.03202, +0.04299]. lane-trainexp then measured that independence from
AlphaGenome's training data is UNKNOWN FOR 13 OF 13 (`18de375`), and the only way past that is a
cell type OUTSIDE AlphaGenome's training domain. lane-oodinventory (`c441b90` rule, `171f5bf`
result) found exactly ONE such candidate with open paired data -- bile duct, UBERON:0002394 -- and
ITS RNA IS SINGLE-NUCLEUS ONLY. A test there would change the cell type AND the assay modality at
once, so the supervisor ruled it not interpretable as it stands.

This module holds the PREREQUISITE that would make it interpretable, and it needs no biosample
lift: does pseudobulked single-nucleus RNA reproduce bulk RNA-seq ON THIS ENDPOINT, in cell types
the reader ALREADY covers that have both?

THE TWO OUTCOMES ARE SYMMETRIC AND BOTH ARE USEFUL. If pseudobulk reproduces bulk on this endpoint,
the bile-duct test becomes interpretable and the supervisor's recommendation to Albert can change.
If it does not, the out-of-domain test cannot be done with that candidate and the biosample question
closes without a lift.

WORDING. Agreement on this endpoint licenses the modality swap FOR THIS ENDPOINT ONLY. It is never
"single-nucleus RNA equals bulk RNA".

THE MINIMUM IS STATED HERE BEFORE ANYTHING IS COUNTED. `MINIMUM` below fixes how many cell types and
how many genes are required before the comparison is meaningful, and `NO_SOURCE_QUALIFIES` is an
acceptable and complete answer: if the covered cell types do not offer both modalities, the
out-of-domain test cannot be made interpretable this way, and saying so is this lane's whole value.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = ROOT / "data" / "results"

#: M2i's committed record. Its ``cells`` list IS the reader's covered set for this lane, and its
#: sources record names one open ENCODE bulk RNA-seq experiment per cell type, already on disk.
FRESHRNA_RESULT = RESULTS_DIR / "freshrna_result.json"
FRESHRNA_SOURCES = RESULTS_DIR / "freshrna_sources.json"
#: M2h's record, for the endpoint carried verbatim and for the fixed scale.
MAGCAL_RESULT = RESULTS_DIR / "magcal_result.json"
MAGCAL_REGISTRATION = RESULTS_DIR / "magcal_registration.json"
#: lane-oodinventory's record, for the one out-of-domain candidate this lane is the prerequisite of.
OODINVENTORY_RESULT = RESULTS_DIR / "oodinventory_candidates.json"

RULE_PATH = RESULTS_DIR / "pseudobulk_rule.json"
INVENTORY_PATH = RESULTS_DIR / "pseudobulk_inventory.json"
REGISTRATION_PATH = RESULTS_DIR / "pseudobulk_registration.json"
RESULT_PATH = RESULTS_DIR / "pseudobulk_result.json"

#: M2h's fitted scale, used ONLY as this fixed constant and never refitted in this lane.
SCALE = 0.2145911381525326

#: The field the covered set and every candidate are matched on. The same field lane-oodinventory
#: matched on, for the same reason: it is the only field present on every ENCODE experiment and on
#: every row of the covered catalogue.
MATCH_FIELD = "biosample_term_name"

#: ENCODE assay titles that are single-nucleus or single-cell RNA. Learned from the portal's own
#: ``assay_title`` facet over released transcription experiments BEFORE this rule was written; no
#: count of qualifying cell types was read first.
SINGLE_CELL_ASSAYS = ("snRNA-seq", "scRNA-seq", "long read scRNA-seq")

#: The bulk assays M2h and M2i read. Carried from `data/results/freshrna_sources.json`.
BULK_ASSAYS = ("polyA plus RNA-seq", "total RNA-seq")

ENCODE = "https://www.encodeproject.org"
ENCODE_LICENCE = (
    "ENCODE Consortium data: free to download and use without restriction under the ENCODE data "
    "use policy; cite the ENCODE Consortium and the experiment accession "
    "(https://www.encodeproject.org/help/citing-encode/)"
)

WORDING_LIMIT = (
    "Agreement on this endpoint licenses the modality swap FOR THIS ENDPOINT ONLY. It is never "
    "'single-nucleus RNA equals bulk RNA', and it is never a statement about any other endpoint, "
    "any other gene set or any cell type not counted here."
)

#: Stated BEFORE counting. Both numbers must be met or the answer is NO_SOURCE_QUALIFIES.
MINIMUM: dict[str, Any] = {
    "cell_types": 4,
    "cell_types_reason": (
        "The registered control is M2h's WITHIN-GENE CELL PERMUTATION: a gene's values are permuted "
        "among its own cells. With n cells a gene has n! distinct permutations, so n = 4 is the "
        "smallest n for which 4! = 24 >= 20 permutations exist and a control p at or below 0.05 is "
        "attainable AT ALL. With three or fewer cell types the registered control cannot be run, so "
        "the comparison cannot be made as registered and no reading from it would be interpretable."
    ),
    "genes": 1000,
    "genes_reason": (
        "M2h's and M2i's registered floor is min_genes = 100, carried verbatim as the endpoint's own "
        "floor. 1,000 is the STRICTER floor this lane requires, because what must be resolvable here "
        "is not whether the endpoint is positive but whether TWO READINGS OF IT AGREE, and the "
        "difference at issue is the size M2h measured (Delta MSE about +0.04, about 2% of variance). "
        "Below 1,000 genes the gene bootstrap's interval on the bulk reading alone is wide enough to "
        "contain both agreement and disagreement, so the comparison would decide nothing."
    ),
    "genes_clause": (
        "A gene counts toward the 1,000 only if it is in M2h's HELD-OUT population AND resolves in "
        "EVERY qualifying cell type in BOTH modalities: exactly one row for its unversioned gene_id "
        "in that cell type's bulk table (rows suffixed _PAR_Y not counted), and a row for the same "
        "gene in that cell type's pseudobulked single-nucleus matrix. No filter on expression and no "
        "filter on any measured value."
    ),
    "below_either": (
        "NO_SOURCE_QUALIFIES, and the lane stops without reading a value. That is a COMPLETE answer: "
        "it says the out-of-domain test cannot be made interpretable this way, so the biosample "
        "question closes without a lift."
    ),
}

SELECTION_RULE: dict[str, Any] = {
    "version": "pseudobulk-1",
    "question": (
        "Which of the reader's 13 covered cell types have BOTH open bulk RNA-seq AND open "
        "single-nucleus (or single-cell) RNA from a public source?"
    ),
    "covered_set": (
        "The thirteen cell keys of M2h's and M2i's committed `cells` list, read from "
        "data/results/freshrna_result.json and pinned by its sha256. Each key's ENCODE "
        "biosample_term_name is taken from data/results/freshrna_sources.json's `biosamples` map, "
        "also pinned. The covered set is NOT re-derived and NOT widened: a cell type absent from "
        "that list is out of scope for this lane whatever data exists for it."
    ),
    "bulk": (
        "Satisfied by construction and named, not re-searched: M2i already holds one released, open "
        "ENCODE bulk RNA-seq experiment per covered cell type (polyA plus RNA-seq or total RNA-seq) "
        "with its experiment accession, its file accessions and their sha256 in "
        "data/results/freshrna_sources.json, and the files are on disk under "
        "data/cache/freshrna/rna/. The bulk side of the comparison is those exact files and no "
        "others."
    ),
    "single_cell": (
        "A RELEASED HUMAN experiment whose assay_title is one of "
        + ", ".join(SINGLE_CELL_ASSAYS)
        + f", and whose {MATCH_FIELD} is EXACTLY the covered cell type's term name. The biosample "
        "must not be genetically modified and must not be treated, for the same reason "
        "lane-oodinventory excluded shRNA, siRNA, CRISPR and CRISPRi RNA-seq: a knockdown is not "
        "the cell type's own RNA. The experiment must carry at least one OPEN processed file from "
        "which a per-gene value per nucleus or per cell can be read -- a gene-by-barcode sparse "
        "count matrix, or a per-gene count table. A browser track, a bam alone or a clustering "
        "figure does not satisfy this clause."
    ),
    "licence": (
        "Open, no application, no controlled access. This excludes GTEx individual-level data and "
        "anything else behind dbGaP, and it excludes the managed-access tiers of the Human Cell "
        "Atlas. An open redistribution under CC BY 4.0 or the ENCODE data use policy satisfies it."
    ),
    "sources_searched": (
        "ENCODE (www.encodeproject.org) exhaustively through its public JSON API for the three "
        "single-cell assay titles above. Beyond ENCODE, each source lane-oodinventory already "
        "characterised is re-examined against THIS lane's clauses, which differ from that lane's: "
        "this lane needs no matched chromatin accessibility, so a source that failed oodinventory's "
        "accessibility or same-source clause may still qualify here. Every source is recorded with "
        "the clause it passes or fails."
    ),
    "same_source_not_required": (
        "This lane deliberately does NOT require the bulk and the single-nucleus data to come from "
        "one source. The question is whether the MODALITY changes the reading, so a cross-source "
        "pair is admissible and the source difference is recorded as a limitation that any "
        "disagreement must be read against."
    ),
    "minimum": MINIMUM,
    "acceptable_answer": (
        "NO_SOURCE_QUALIFIES is an acceptable and COMPLETE answer, and it is the answer whenever "
        "either number in `minimum` is not met."
    ),
    "procedure_registered_separately": (
        "If and only if sources qualify, the pseudobulk aggregation -- how nuclei are aggregated to "
        "a per-gene value, stated exactly -- and the endpoint carried verbatim from 7d, 7h and M2h "
        "are REGISTERED in data/results/pseudobulk_registration.json and committed as their own "
        "commit BEFORE the run. No measured value is read before that commit."
    ),
    "wording_limit": WORDING_LIMIT,
    "no_biosample": (
        "This lane adds NO biosample. It works only in cell types the reader already covers, makes "
        "no AlphaGenome request and no deletion request, and touches nothing in area J, no node "
        "model and no BioForge calibration."
    ),
}


def sha256_of(path: Path) -> str:
    """The sha256 of a file on disk, for the manifest and for pinning an input."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def covered_cells() -> list[dict[str, Any]]:
    """The thirteen covered cell types with the bulk experiment and files M2i already read.

    Read from the two committed records, never re-derived. Each entry carries the cell key as M2h
    and M2i write it, the ENCODE ``biosample_term_name`` to match a single-cell experiment against,
    the bulk assay and experiment accession, and the bulk files with their sha256.
    """
    cells = json.loads(FRESHRNA_RESULT.read_text())["cells"]
    sources = json.loads(FRESHRNA_SOURCES.read_text())["biosamples"]
    out: list[dict[str, Any]] = []
    for key in cells:
        entry = sources[key]
        out.append(
            {
                "cell": key,
                MATCH_FIELD: entry["biosample"],
                "bulk_assay": entry["assay"],
                "bulk_experiment": entry["experiment"],
                "bulk_files": [
                    {"accession": f["accession"], "path": f["path"], "sha256": f["sha256"]}
                    for f in entry["files"]
                ],
            }
        )
    return out


def verdict(qualifying_cells: int, qualifying_genes: int | None) -> str:
    """``QUALIFIES`` only when BOTH numbers in :data:`MINIMUM` are met, else ``NO_SOURCE_QUALIFIES``.

    ``qualifying_genes`` is ``None`` when the cell-type clause already failed and no gene count was
    taken, which is the ordered form of the rule: the gene count is never read to rescue a cell
    count that did not clear.
    """
    if qualifying_cells < MINIMUM["cell_types"]:
        return "NO_SOURCE_QUALIFIES"
    if qualifying_genes is None or qualifying_genes < MINIMUM["genes"]:
        return "NO_SOURCE_QUALIFIES"
    return "QUALIFIES"
