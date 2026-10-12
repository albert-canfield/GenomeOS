# SPDX-License-Identifier: AGPL-3.0-or-later
"""Open single-cell RNA for the reader's covered cell types OUTSIDE ENCODE, under an adopted rule.

WHY THIS LANE EXISTS. lane-pseudobulk committed its rule and its minimum first (`2b9cd0a`), then
asked which of the reader's thirteen covered cell types have BOTH open bulk RNA-seq and open
single-cell or single-nucleus RNA, so that pseudobulk could be checked against bulk on M2h's
endpoint. It found 2 of 13 against a minimum of 4 -- `NO_SOURCE_QUALIFIES` (`21b58b6` inventory,
`f26c5ae`, `eab497c`) -- and it searched ENCODE exhaustively plus Tabula Sapiens, the Human Cell
Atlas and the GTEx single-nucleus pilot. IT DID NOT SEARCH GEO, CELLxGENE OR THE 10x PUBLIC
DATASETS, AND SAID SO IN ITS OWN RECORD. The coordinator nevertheless told Albert the question was
closed, on an assumption that a wider search would not change the count. This lane searches the
three archives that were named as unsearched, and settles it.

THE RULE IS ADOPTED, NOT RE-DERIVED. Every clause and both numbers come from
`data/results/pseudobulk_rule.json` as committed in `2b9cd0a`: 4 cell types, 1,000 genes, open
access with no application, no genetically modified and no treated biosample, at least one open
processed file from which a per-gene value per cell or per nucleus can be read, the covered set
neither re-derived nor widened, and bulk-and-single-cell from one source NOT required. Nothing here
relaxes a clause and nothing here tightens one. :data:`MINIMUM` IS IMPORTED FROM THAT LANE'S MODULE
rather than copied, so a drift between the two is impossible by construction.

ONE CLAUSE HAS TO BE READ ACROSS ARCHIVES, AND THE READING IS STATED BEFORE ANY COUNT. The adopted
single-cell clause is written in ENCODE's vocabulary: "a RELEASED HUMAN experiment whose
assay_title is one of snRNA-seq, scRNA-seq, long read scRNA-seq, and whose biosample_term_name is
EXACTLY the covered cell type's term name". GEO, CELLxGENE and 10x have no `assay_title` and no
`biosample_term_name` field at all, so the clause cannot be applied as a field comparison without
making every non-ENCODE source fail BY CONSTRUCTION -- which would make this lane's answer vacuous
and would be a tightening. :data:`BIOSAMPLE_IDENTITY_POLICY` and :data:`TERM_MATCH_POLICY` state
how it is read instead, and BOTH counts are reported: the one where the archive's own label is the
IDENTICAL STRING, and the one where it is a different string for the same biosample.

THIS LANE READS NO RNA VALUE, FITS NOTHING AND RUNS NO COMPARISON. Its output is a list and a
count. It adds NO biosample: it looks only at cell types the reader already covers, makes no
AlphaGenome request and no deletion request, and touches nothing in area J, no node model and no
BioForge calibration (the stop list of 2026-09-27, reinstated, ROADMAP section 5).

WHAT A COUNT HERE IS AND IS NOT. An archive listing is what is PUBLISHED, not what is USABLE. No
experimental data file is fetched by this lane: a candidate's openness is established by an HTTP
HEAD -- a status and a byte count -- and its per-gene-per-cell content is taken from the archive's
own declaration of the file (an `.h5ad` whose raw counts the archive names, a filtered
gene-barcode matrix) and NOT from reading the file. A candidate that clears every clause here is a
candidate; whether its matrix resolves 1,000 of M2h's held-out genes is the SECOND number of the
adopted minimum and this lane does not count it.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from genomeos.attribution import pseudobulk as pb

ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = ROOT / "data" / "results"

#: lane-pseudobulk's committed rule and its inventory. Both are read, neither is written.
ADOPTED_RULE = pb.RULE_PATH
PREDECESSOR_INVENTORY = pb.INVENTORY_PATH

INVENTORY_PATH = RESULTS_DIR / "scsweep_inventory.json"

#: The adopted minimum, IMPORTED and never copied, so this lane cannot drift from `2b9cd0a`.
MINIMUM: dict[str, Any] = pb.MINIMUM
ADOPTED_RULE_VERSION = pb.SELECTION_RULE["version"]
WORDING_LIMIT = pb.WORDING_LIMIT

BIOSAMPLE_IDENTITY_POLICY = (
    "The adopted single-cell clause requires the single-cell experiment's BIOSAMPLE to be the "
    "covered cell type. Applied to an archive with no `biosample_term_name` field, that is read as "
    "biosample IDENTITY and nothing wider: the cells the candidate file holds must BE the covered "
    "cell type, not a tissue that contains it and not a related type. A whole dissociated tissue "
    "counts for a covered type that IS a tissue -- `ovary`, `testis` -- which is exactly the "
    "relation of the ENCODE `ovary` snRNA-seq experiments lane-pseudobulk itself ADMITTED, and it "
    "does NOT count for a covered type that is a cell, because selecting that cell type's "
    "annotated subset out of a tissue is the cell-composition decision the adopted rule does not "
    "license. A purified or cultured preparation of the covered cell type DOES count, because its "
    "biosample is that cell type."
)

TERM_MATCH_POLICY = (
    "The adopted clause says EXACTLY, and what that word protects is stated in the rule's own "
    "`covered_set` clause: the covered set is NOT re-derived and NOT widened. It guards against "
    "admitting a cell type that is not in the thirteen; it cannot be read as a requirement that "
    "another archive spell a cell type the way ENCODE spells it, because no other archive has the "
    "field. So BOTH counts are reported and neither is hidden: `exact_label_count` admits a "
    "candidate only when the archive's own label is the IDENTICAL STRING to the covered "
    "biosample_term_name, and `same_biosample_count` also admits one whose label is a different "
    "string for the same biosample (`CD14+ Monocytes` for `CD14-positive monocyte`). The headline "
    "count of this lane is the STRICTER of the two, `exact_label_count`."
)

#: Every archive this lane searched or did not search, with what it is.
ARCHIVES: dict[str, dict[str, str]] = {
    "GEO": {
        "name": "NCBI Gene Expression Omnibus",
        "endpoint": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils (db=gds) and "
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi (targ=self&form=text&view=brief)",
        "searched": "yes",
        "licence": "NCBI GEO holds public submissions available without restriction and without "
        "an application; no controlled-access tier is involved (dbGaP and EGA, which the adopted "
        "licence clause excludes, are separate archives and nothing here comes from them)",
        "limit": "GEO is an archive of independent submissions with no common pipeline and no "
        "common quantification, so it CANNOT be swept exhaustively the way ENCODE's typed JSON "
        "API can. This lane recorded one search count per covered cell type and then read "
        "candidate submissions' own metadata. Its negatives are therefore NOT a claim that no "
        "such submission exists: they are a claim about these named queries.",
    },
    "CELLxGENE": {
        "name": "CZ CELLxGENE Discover",
        "endpoint": "https://api.cellxgene.cziscience.com/curation/v1/datasets",
        "searched": "yes, exhaustively: the whole dataset listing was read and every human "
        "dataset's tissue labels were matched against the covered term names",
        "licence": "CZ CELLxGENE Discover distributes curated datasets openly over HTTPS with no "
        "registration and no application; its submission policy requires CC BY 4.0, and each "
        "dataset carries its own `citation` string naming the publication and the dataset version",
        "limit": "CELLxGENE's `cell_type` field is a CL annotation of cells INSIDE a dataset; a "
        "dataset's biosample is its `tissue`. Under BIOSAMPLE_IDENTITY_POLICY a covered cell type "
        "appearing only as a `cell_type` label inside a tissue dataset does NOT qualify.",
    },
    "10x public datasets": {
        "name": "10x Genomics public example datasets",
        "endpoint": "https://cf.10xgenomics.com/samples/cell-exp/ (the CDN the dataset pages "
        "link to) and https://www.10xgenomics.com/datasets",
        "searched": "yes, for the covered cell types; the catalogue page itself answered HTTP 429 "
        "to this lane, so the per-dataset CDN files were read directly instead",
        "licence": "10x states on its dataset pages that the public example datasets are offered "
        "under CC BY 4.0. THIS LANE COULD NOT FETCH THAT PAGE (HTTP 429) and does not quote it: "
        "what it verified itself is that the matrix downloads over HTTPS with no registration and "
        "no application (HTTP 200 on a HEAD, with a byte count).",
        "limit": "The catalogue is small and is mostly PBMC, tissue and tumour preparations. Only "
        "the FACS-purified immune populations of Zheng et al. 2017 meet a covered cell type.",
    },
    "Tabula Sapiens, Human Cell Atlas, GTEx snRNA pilot": {
        "name": "the three sources lane-pseudobulk characterised beyond ENCODE",
        "endpoint": "n/a",
        "searched": "re-examined, not re-searched: Tabula Sapiens is reached here THROUGH "
        "CELLxGENE, which distributes it per tissue, and that is where this lane's testis "
        "candidate comes from",
        "licence": "as lane-pseudobulk recorded: GTEx individual-level data and the HCA managed "
        "tiers are excluded by the adopted licence clause",
        "limit": "lane-pseudobulk excluded Tabula Sapiens for the covered cell LINES (no match "
        "exists) and for the covered TISSUES (it read a tissue match as needing a "
        "cell-composition decision). The first half holds. THE SECOND HALF IS WHERE THIS LANE "
        "DIFFERS, and the difference is in the reading and not in the rule: pseudobulking a whole "
        "dissociated testis needs no composition decision, and it is the same relation as the "
        "ENCODE whole-ovary snRNA-seq that lane-pseudobulk ADMITTED.",
    },
    "dbGaP, EGA": {
        "name": "controlled-access archives",
        "endpoint": "n/a",
        "searched": "no, and deliberately",
        "licence": "controlled access behind an application",
        "limit": "excluded by the adopted licence clause, which this lane does not relax",
    },
}

#: The candidates this lane names, each with the query that found it. The clause evaluation is NOT
#: stored here: the sweep recomputes it live from each archive's own record, so a rebuild can
#: disagree with this lane and say so.
CANDIDATES: tuple[dict[str, Any], ...] = (
    {
        "cell": "K562",
        "archive": "GEO",
        "accession": "GSM9224207",
        "series": "GSE306708",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9224207",
        "assay": "inCAGE single-cell RNA-seq (droplet; library strategy RNA-Seq)",
        "files": (
            {
                "name": "GSM9224207_Fresh_K562.h5ad",
                "url": "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9224nnn/GSM9224207/suppl/"
                "GSM9224207_Fresh_K562.h5ad",
                "what": "AnnData matrix of the FRESH (unfixed) K562 cells, cells by genes",
            },
        ),
        "found_by": 'GEO esearch db=gds, gsm[Entry Type], "Homo sapiens"[Organism], "K562"[All '
        'Fields] AND ("single cell" OR "single-cell" OR "single nucleus" OR "scRNA" OR "snRNA"), '
        "then the sample's own brief read for its characteristics",
        "label_in_archive": "K562",
        "notes": (
            "The submission's characteristics read `cell line: K562, L1210`: the sample is a 1:1 "
            "human/mouse preparation and the supplementary files are split per line, so the file "
            "admitted here holds the K562 cells alone. The sample's OTHER file, "
            "`GSM9224207_Fixed_K562.h5ad`, is PFA-fixed, which is a treatment of the cells; only "
            "the FRESH file is admitted and the fixed one is recorded as rejected on that clause."
        ),
    },
    {
        "cell": "testis",
        "archive": "CELLxGENE",
        "accession": "ec6ef004-5f30-4a40-b282-980cca20a561",
        "series": "collection e5f58829-1a66-40b5-a624-9046778e74f5 (Tabula Sapiens)",
        "url": "https://api.cellxgene.cziscience.com/curation/v1/collections/"
        "e5f58829-1a66-40b5-a624-9046778e74f5/datasets/ec6ef004-5f30-4a40-b282-980cca20a561",
        "assay": "10x 3' v3 single-cell RNA (suspension type: cell)",
        "files": (),
        "found_by": "the whole CELLxGENE dataset listing read once, then every human dataset whose "
        "`tissue` labels contain a string identical to a covered biosample_term_name",
        "label_in_archive": "testis",
        "notes": (
            "Tabula Sapiens - Testis: a whole dissociated adult testis, `disease: normal`, raw "
            "counts declared at `raw.X`. Its biosample is the testis, which is the covered type "
            "itself, so no cell-composition decision arises. lane-pseudobulk excluded Tabula "
            "Sapiens; this lane reaches it through CELLxGENE's per-tissue distribution and "
            "disagrees with the SECOND half of that exclusion for the reason in ARCHIVES."
        ),
    },
    {
        "cell": "CD14_positive_monocyte",
        "archive": "10x public datasets",
        "accession": "cd14_monocytes (standard 1.1.0)",
        "series": "Zheng et al. 2017, the FACS-purified immune populations",
        "url": "https://www.10xgenomics.com/datasets/cd-14-plus-monocytes-1-standard-1-1-0",
        "assay": "10x Chromium single-cell 3' RNA (Cell Ranger 1.1.0; 2,612 cells in the web "
        "summary this lane read)",
        "files": (
            {
                "name": "cd14_monocytes_filtered_gene_bc_matrices.tar.gz",
                "url": "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd14_monocytes/"
                "cd14_monocytes_filtered_gene_bc_matrices.tar.gz",
                "what": "filtered gene-barcode sparse count matrix (Matrix Market)",
            },
        ),
        "found_by": "the 10x public dataset catalogue, for the covered cell types; the dataset's "
        "own Cell Ranger web summary read from the CDN for its sample id and label",
        "label_in_archive": "CD14+ Monocytes",
        "notes": (
            "FACS-purified CD14-positive monocytes from a healthy donor: the biosample IS the "
            "covered cell type, unmodified and untreated. ITS LABEL IS NOT THE IDENTICAL STRING "
            "to ENCODE's `CD14-positive monocyte`, so under TERM_MATCH_POLICY it counts in "
            "`same_biosample_count` and NOT in the headline `exact_label_count`."
        ),
    },
)

#: Named candidates that FAIL a clause. Recorded because two of them are where the lift was
#: expected to come from, and it does not come from there.
REJECTED: tuple[dict[str, Any], ...] = (
    {
        "cell": "K562",
        "archive": "GEO / ENCODE Perturb-seq",
        "accession": "e.g. GSM9550664 (ENCODE library ENCLB761BDY), GSM9964752-5 (K562 Zim3-dCas9 "
        "pooled dual-guide screen), GSE308682, GSE283574, GSE288996",
        "clause_failed": "not genetically modified",
        "why": (
            "THE SUPERVISOR EXPECTED THE LIFT TO COME FROM PERTURB-SEQ NON-TARGETING CONTROL "
            "CELLS, AND UNDER THE ADOPTED RULE IT CANNOT. A non-targeting control cell in a "
            "Perturb-seq or CRISPRi screen carries a lentivirally integrated guide vector and a "
            "Cas9, dCas9 or Zim3-dCas9 effector; that is a genetically modified biosample, which "
            "the adopted single-cell clause excludes in the same words it excludes shRNA, siRNA, "
            "CRISPR and CRISPRi RNA-seq -- 'a knockdown is not the cell type's own RNA'. This lane "
            "does not relax that clause to admit them. The K562 candidate admitted above is an "
            "unmodified, unfixed K562 culture and owes nothing to Perturb-seq."
        ),
    },
    {
        "cell": "K562",
        "archive": "CELLxGENE",
        "accession": "38a78e91-e361-45c8-bae3-1d6d5108d18e (sci-Plex, tissue label `K-562`)",
        "clause_failed": "not treated",
        "why": (
            "The only CELLxGENE dataset whose tissue label is a covered cell line is 'Massively "
            "multiplex chemical transcriptomics at single-cell resolution': a chemical "
            "perturbation screen in which every arm receives a compound or a vehicle, so there is "
            "no untreated arm to admit. Finding one would mean reading the dataset's own cell "
            "annotations, which are inside a 321 MB matrix this lane does not fetch."
        ),
    },
    {
        "cell": "HepG2",
        "archive": "GEO",
        "accession": "none admitted",
        "clause_failed": "not genetically modified, not treated, or no open per-cell matrix",
        "why": (
            "THE SUPERVISOR EXPECTED HepG2 AND THIS SWEEP DID NOT FIND IT. 689 GEO samples match "
            "HepG2 under the single-cell terms, and the ones these queries surfaced are "
            "overexpression, CRISPRi, vector-control or drug-treatment arms, or preparations of "
            "another cell type in a study that merely mentions HepG2. Three query forms were "
            "tried: the plain single-cell sample sweep, the same narrowed by `untreated OR "
            "parental OR wild-type OR 'wild type'` (24 samples, none an unmodified HepG2 "
            "single-cell library), and the same narrowed by `benchmark OR mixture OR 'cell lines' "
            "OR titration` (364 samples, none an unmodified HepG2 single-cell library). THAT IS A "
            "STATEMENT ABOUT THESE QUERIES AND NOT ABOUT GEO: an archive listing is what is "
            "published, not what is usable, and this lane does not claim no such submission "
            "exists."
        ),
    },
    {
        "cell": "H1",
        "archive": "GEO",
        "accession": "e.g. GSM7976019 (`cell line: H1`, `genotype: WT`)",
        "clause_failed": "not treated",
        "why": (
            "The day-0 sample of a differentiation time course is an exact label match on an "
            "unmodified H1, and its own characteristics read `treatment: BMP4+CHIR` with a "
            "treatment protocol to match, so it is a treated biosample. No untreated H1 "
            "single-cell library with an open per-cell matrix was surfaced by these queries."
        ),
    },
    {
        "cell": "astrocyte, hepatocyte, keratinocyte, cardiac muscle cell, CD14-positive monocyte",
        "archive": "CELLxGENE",
        "accession": "271, 89, 73, 39 and 66 datasets mention these labels respectively",
        "clause_failed": "biosample identity",
        "why": (
            "Each appears as a CL `cell_type` annotation of cells INSIDE a tissue dataset, and NO "
            "human CELLxGENE dataset carries any of them as its `tissue`, which is where the "
            "biosample is stated. Admitting one would mean selecting that cell type's annotated "
            "subset out of a tissue, which is the cell-composition decision the adopted rule does "
            "not license -- the same ground on which lane-pseudobulk excluded Tabula Sapiens and "
            "the HCA, and this lane keeps it."
        ),
    },
    {
        "cell": "IMR_90, SK_N_SH, GM12878, hepatocyte, keratinocyte, cardiac_muscle_cell, astrocyte",
        "archive": "GEO",
        "accession": "none examined beyond the per-cell-type search counts",
        "clause_failed": "not evaluated",
        "why": (
            "The sweep recorded a GEO search count for all thirteen covered cell types, and read "
            "candidate submissions' metadata only for the cell types where one was needed: the "
            "two the supervisor named and the ones that could reach the minimum. GM12878 is "
            "already supplied by ENCODE, so no GEO candidate was sought for it. For the rest this "
            "lane states plainly that it did NOT evaluate individual submissions, so their "
            "absence from the count is UNSEARCHED and not a negative."
        ),
    },
)


def covered_terms() -> dict[str, str]:
    """The thirteen covered cell keys mapped to their ENCODE biosample_term_name.

    Read through lane-pseudobulk's own loader, so the covered set is the one `2b9cd0a` pinned and
    is neither re-derived nor widened here.
    """
    return {c["cell"]: c[pb.MATCH_FIELD] for c in pb.covered_cells()}


def predecessor_qualifying() -> list[str]:
    """The cell types lane-pseudobulk's committed inventory already admits, read from that file."""
    return list(json.loads(PREDECESSOR_INVENTORY.read_text())["qualifying_cells"])


def verdict(exact_label_count: int) -> str:
    """The cell-type half of the adopted minimum, and NOTHING MORE IS CLAIMED.

    `NO_SOURCE_QUALIFIES` is the adopted rule's own complete answer below 4. At or above 4 this
    returns `CELL_TYPE_MINIMUM_REACHED` and NOT `QUALIFIES`: the adopted rule's `QUALIFIES`
    requires BOTH numbers, and the second -- 1,000 of M2h's held-out genes resolving in every
    qualifying cell type in both modalities -- cannot be counted without reading the matrices,
    which this lane does not do. `pb.verdict` is called beside this to show that it agrees.
    """
    if exact_label_count < MINIMUM["cell_types"]:
        return "NO_SOURCE_QUALIFIES"
    return "CELL_TYPE_MINIMUM_REACHED"
