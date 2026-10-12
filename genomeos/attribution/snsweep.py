# SPDX-License-Identifier: AGPL-3.0-or-later
"""Open single-NUCLEUS RNA for the reader's covered cell types, under lane-pseudobulk's rule.

WHY THIS LANE EXISTS. M2h (`ac5580b`) found the compiled rules fit measured RNA better scaled by
one fitted constant, +0.2146; M2i (`f6e71c6`) replicated it with nothing fitted; lane-scalesplit
(`ec9cfd4`) showed it survives a chromosome-held-out split. Its one unresolved limit is that
independence from AlphaGenome's training domain is UNKNOWN FOR 13 OF 13 (`18de375`), and the only
cell type outside that domain with open paired data is bile duct (`171f5bf`) -- whose RNA is
SINGLE-NUCLEUS ONLY. lane-pbrun has just licensed the modality swap on M2h's endpoint
(`ccd5379` registration, `22cc987` result): AGREEMENT, the paired difference +0.00757935
[98.33% CI -0.00243939, +0.0174182] inside the pre-registered margin +-0.0216956, on 8,051 genes
across GM12878, K562, ovary and testis. BUT ONLY ONE OF THOSE FOUR IS SINGLE-NUCLEUS (ovary),
which is bile duct's modality, and the supervisor ruled that n = 1 cannot even run the registered
within-gene permutation, so it licenses nothing on its own. THIS LANE ESTABLISHES WHETHER THE
SINGLE-NUCLEUS LICENCE CAN BE MADE, OR CANNOT.

THE RULE IS ADOPTED, NOT RE-DERIVED. Every clause and both numbers come from
`data/results/pseudobulk_rule.json` as committed in `2b9cd0a`: 4 cell types, 1,000 genes, open
access with no application, no genetically modified and no treated biosample, at least one open
processed file from which a per-gene value per nucleus can be read, the covered set neither
re-derived nor widened, and bulk-and-single-nucleus from one source NOT required. :data:`MINIMUM`
IS IMPORTED from that lane's module rather than copied, as `lane-scsweep` did it, so a drift
between the lanes is impossible by construction. The two readings that lane had to make across
archives -- :data:`BIOSAMPLE_IDENTITY_POLICY` and :data:`TERM_MATCH_POLICY` -- are likewise
IMPORTED from `genomeos/attribution/scsweep.py` and not restated, for the same reason.

THE ONE CLAUSE THIS LANE NARROWS IS THE MODALITY, AND IT NARROWS IT, NEVER RELAXES IT.
lane-pseudobulk's single-cell clause admits three ENCODE assay titles, `snRNA-seq`, `scRNA-seq`
and `long read scRNA-seq`. Here ONLY the single-NUCLEUS one counts: a single-CELL experiment does
not count towards this minimum however well it clears every other clause. That is a STRICTER
reading of the adopted rule and it cannot inflate the count. :data:`MODALITY_POLICY` states, per
archive, how the distinction is made from the archive's own words.

THIS LANE READS NO RNA VALUE, FITS NOTHING AND RUNS NO COMPARISON. Its output is a list and a
count. It adds NO biosample: it looks only at cell types the reader already covers, makes no
AlphaGenome request and no deletion request, and touches nothing in area J, no node model and no
BioForge calibration (the stop list of 2026-09-27, reinstated, ROADMAP section 5).

WHAT A COUNT HERE IS AND IS NOT. An archive listing is what is PUBLISHED, not what is USABLE. No
experimental data file is fetched: a candidate's openness is an HTTP HEAD -- a status and a byte
count -- and its per-gene-per-nucleus content is taken from the archive's own declaration of the
file and NOT from reading it. A candidate that clears every clause here is a candidate; whether
its matrix resolves 1,000 of M2h's held-out genes is the SECOND number of the adopted minimum and
this lane does not count it.
"""

from __future__ import annotations

import json
from typing import Any

from genomeos.attribution import pseudobulk as pb
from genomeos.attribution import scsweep as sc

RESULTS_DIR = pb.RESULTS_DIR

#: lane-pseudobulk's committed rule and inventory, and lane-scsweep's. All read, none written.
ADOPTED_RULE = pb.RULE_PATH
PSEUDOBULK_INVENTORY = pb.INVENTORY_PATH
SCSWEEP_INVENTORY = sc.INVENTORY_PATH

INVENTORY_PATH = RESULTS_DIR / "snsweep_inventory.json"

#: The adopted minimum, IMPORTED and never copied, so this lane cannot drift from `2b9cd0a`.
MINIMUM: dict[str, Any] = pb.MINIMUM
ADOPTED_RULE_VERSION = pb.SELECTION_RULE["version"]
WORDING_LIMIT = pb.WORDING_LIMIT

#: The two cross-archive readings lane-scsweep had to make, IMPORTED and not restated.
BIOSAMPLE_IDENTITY_POLICY = sc.BIOSAMPLE_IDENTITY_POLICY
TERM_MATCH_POLICY = sc.TERM_MATCH_POLICY

#: The single-NUCLEUS member of the adopted rule's three assay titles. The other two are
#: single-CELL and are named here so the exclusion is explicit rather than implied.
SINGLE_NUCLEUS_ASSAY = "snRNA-seq"
SINGLE_CELL_ASSAYS_NOT_COUNTED = ("scRNA-seq", "long read scRNA-seq")

MODALITY_POLICY = (
    "SINGLE-NUCLEUS ONLY, AND A SINGLE-CELL EXPERIMENT DOES NOT COUNT TOWARDS THIS MINIMUM "
    "HOWEVER WELL IT CLEARS EVERY OTHER CLAUSE. This is a NARROWING of the adopted single-cell "
    "clause, never a relaxation of it: every candidate admitted here would also be admitted by "
    "`2b9cd0a`, and candidates `2b9cd0a` admits are dropped here. Each archive states the "
    "modality in its own vocabulary and each is read from that archive's own words: (1) ENCODE -- "
    "`assay_title` EXACTLY `snRNA-seq`; `scRNA-seq` and `long read scRNA-seq` are single-CELL and "
    "are excluded, which is why ENCODE's GM12878 scRNA-seq, admitted by `2b9cd0a` and used by "
    "lane-pbrun, is NOT counted here. (2) CELLxGENE -- the schema's own `suspension_type` field "
    "EXACTLY `nucleus`; `cell` is excluded, which is why Tabula Sapiens - Testis, admitted by "
    "`12b0c4e` and used by lane-pbrun, is NOT counted here. (3) GEO -- no assay vocabulary "
    "exists, so the SAMPLE'S OWN WORDS are read and all three must agree: its "
    "`library_strategy` is RNA-Seq (which excludes an ATAC arm of the same multiome), and its "
    "`extract_protocol`, title or series design says nuclei (`Single nuclei`, `snRNA-seq`, "
    "`single-nucleus`, `DroNc-seq`, `snPATHO-seq`) and does not say the library was made from "
    "whole cells. The free-text word `nuclei` ALONE in a series summary is not taken as a "
    "modality: it matches nuclear extraction and nuclear RNA in bulk assays, and every GEO "
    "candidate's modality here comes from its own sample record."
)

#: Every archive this lane searched or did not search, with what it is.
ARCHIVES: dict[str, dict[str, str]] = {
    "ENCODE": {
        "name": "ENCODE portal",
        "endpoint": "https://www.encodeproject.org/search/ (type=Experiment, status=released, "
        "assay_title, the biosample_ontology.term_name facet)",
        "searched": "yes, exhaustively: the term_name facet of every RELEASED HUMAN experiment of "
        "each of the adopted rule's three assay titles was read, and the single-nucleus one was "
        "intersected with the thirteen covered term names",
        "licence": "ENCODE Consortium data: free to download and use without restriction under "
        "the ENCODE data use policy, no application",
        "limit": "The facet is the portal's own count over `biosample_ontology.term_name`, which "
        "is the adopted rule's match field. A covered term absent from the facet has NO released "
        "human experiment of that assay title, which IS a negative and not an unsearched gap.",
    },
    "CELLxGENE": {
        "name": "CZ CELLxGENE Discover",
        "endpoint": "https://api.cellxgene.cziscience.com/curation/v1/datasets",
        "searched": "yes, exhaustively: the whole dataset listing was read, and every human "
        "dataset whose `suspension_type` includes `nucleus` had its tissue labels matched against "
        "the covered term names",
        "licence": "CZ CELLxGENE Discover distributes curated datasets openly over HTTPS with no "
        "registration and no application; its submission policy requires CC BY 4.0",
        "limit": "CELLxGENE's `cell_type` field is a CL annotation of cells INSIDE a dataset; a "
        "dataset's biosample is its `tissue`. Under BIOSAMPLE_IDENTITY_POLICY a covered cell type "
        "appearing only as a `cell_type` label inside a tissue dataset does NOT qualify.",
    },
    "GEO": {
        "name": "NCBI Gene Expression Omnibus",
        "endpoint": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils (db=gds) and "
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi (targ=self&form=text&view=brief)",
        "searched": "yes, for all thirteen covered cell types at BOTH the sample and the series "
        "level under single-NUCLEUS-only terms, and then candidate submissions' own records read",
        "licence": "NCBI GEO holds public submissions available without restriction and without "
        "an application; dbGaP and EGA, which the adopted licence clause excludes, are separate "
        "archives and nothing here comes from them",
        "limit": "GEO is an archive of independent submissions with no common pipeline and no "
        "common quantification, so it CANNOT be swept exhaustively the way ENCODE's typed JSON "
        "API can. Its negatives here are NOT a claim that no such submission exists: they are a "
        "claim about these named queries and the submissions they surfaced.",
    },
    "Tabula Sapiens, Human Cell Atlas, GTEx snRNA pilot, 10x public datasets": {
        "name": "the sources lane-pseudobulk and lane-scsweep characterised",
        "endpoint": "n/a",
        "searched": "re-examined against the MODALITY clause, not re-searched: Tabula Sapiens is "
        "reached through CELLxGENE, where every one of its datasets is `suspension_type: cell`; "
        "the 10x public CD14-positive monocyte dataset is a whole-CELL 3' library; the HCA open "
        "releases and the GTEx snRNA pilot fail the covered set, the GTEx individual-level data "
        "also the licence clause",
        "licence": "as those lanes recorded",
        "limit": "The three candidates lane-scsweep and lane-pseudobulk admitted for K562, "
        "GM12878 and testis are ALL single-CELL, so none of them counts here. That is the "
        "narrowing doing its work, and it is why this lane's count is not that lane's count.",
    },
    "dbGaP, EGA": {
        "name": "controlled-access archives",
        "endpoint": "n/a",
        "searched": "no, and deliberately",
        "licence": "controlled access behind an application",
        "limit": "excluded by the adopted licence clause, which this lane does not relax",
    },
}

#: The single-NUCLEUS candidates this lane names, each with the query that found it. The clause
#: evaluation is NOT stored here: the sweep recomputes it live from each archive's own record, so a
#: rebuild can disagree with this lane and say so.
CANDIDATES: tuple[dict[str, Any], ...] = (
    {
        "cell": "ovary",
        "archive": "ENCODE",
        "accession": "ENCSR938CZJ, ENCSR297ZTQ, ENCSR876EUH, ENCSR362HUP, ENCSR298MPV",
        "series": "the five released human snRNA-seq experiments on biosample ovary that "
        "lane-pseudobulk's committed inventory (`21b58b6`) admitted",
        "url": "https://www.encodeproject.org/experiments/ENCSR938CZJ/",
        "assay": "snRNA-seq (ENCODE assay_title), 10x multiome",
        "modality_evidence": "ENCODE assay_title is EXACTLY `snRNA-seq`",
        "label_in_archive": "ovary",
        "found_by": "the ENCODE assay_title=snRNA-seq facet over released human experiments, "
        "intersected with the thirteen covered biosample_term_names; the admitted experiments and "
        "their files are carried from lane-pseudobulk's committed inventory and not re-searched",
        "notes": (
            "THE ONE CELL TYPE THAT QUALIFIES, and it is the same one lane-pbrun's single-nucleus "
            "arm rested on. It is not new: it is the n = 1 the supervisor ruled licenses nothing "
            "on its own. Counting it again here changes nothing about that ruling; what this lane "
            "settles is whether three more can be found beside it."
        ),
    },
    {
        "cell": "ovary",
        "archive": "CELLxGENE",
        "accession": "3a702020-cf66-4ef5-9d11-4a5e835b7efb",
        "series": "collection of 'The regulatory landscapes of the human ovarian aging'",
        "url": "https://api.cellxgene.cziscience.com/curation/v1/datasets/"
        "3a702020-cf66-4ef5-9d11-4a5e835b7efb",
        "assay": "10x 3' v3, suspension type nucleus (snRNA-seq)",
        "modality_evidence": "CELLxGENE `suspension_type` is EXACTLY `nucleus`",
        "label_in_archive": "ovary",
        "found_by": "the whole CELLxGENE dataset listing read once, then every human dataset "
        "whose `suspension_type` includes `nucleus` and whose `tissue` labels contain a string "
        "identical to a covered biosample_term_name",
        "notes": (
            "A SECOND SOURCE FOR A CELL TYPE ALREADY COUNTED, AND IT ADDS NOTHING TO THE COUNT. "
            "The adopted minimum counts CELL TYPES, not sources, and ovary is one cell type "
            "whether one archive supplies it or two. It is recorded because it is the honest "
            "answer to 'did the wider search find anything in ovary', and because it would be a "
            "second source for a future ovary-only reading -- not for this count."
        ),
    },
    {
        "cell": "SK_N_SH",
        "archive": "GEO",
        "accession": "GSM8159281",
        "series": "GSE262189, joint snRNA-seq and snATAC-seq of neuroblastoma cell lines",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8159281",
        "assay": "10x Genomics Multiome GEX on isolated nuclei (library strategy RNA-Seq)",
        "modality_evidence": "the sample's own `extract_protocol` reads `Single nuclei` and its "
        "`library_strategy` is RNA-Seq; the ATAC arm of the same cell line is a separate sample "
        "(GSM8159280) and is not this one",
        "label_in_archive": "SK.N.SH",
        "found_by": 'GEO esearch db=gds at the series level, \'"SK-N-SH"[All Fields] AND '
        'gse[Entry Type] AND "Homo sapiens"[Organism]\' narrowed by the single-nucleus terms, '
        "then every sample of GSE262189 read for its cell line, library strategy and files",
        "files": (
            {
                "name": "GSM8159281_SK.N.SH_filtered_feature_bc_matrix.h5",
                "url": "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8159nnn/GSM8159281/suppl/"
                "GSM8159281_SK.N.SH_filtered_feature_bc_matrix.h5",
                "what": "Cell Ranger filtered feature-barcode count matrix, genes by nuclei",
            },
        ),
        "notes": (
            "AN UNMODIFIED, UNTREATED SK-N-SH CULTURE WHOSE NUCLEI WERE SEQUENCED, with an open "
            "per-gene-per-nucleus matrix -- and ITS LABEL IS NOT THE IDENTICAL STRING to ENCODE's "
            "`SK-N-SH`: the submission spells it `SK.N.SH`. Under the IMPORTED TERM_MATCH_POLICY "
            "it therefore counts in `same_biosample_count` and NOT in the headline "
            "`exact_label_count`. That policy is lane-scsweep's and is imported unchanged; this "
            "lane does not reread it to let a near-spelling into the headline. The sample's "
            "characteristics read `tissue: Neuroblastoma tumor, cell line: SK.N.SH, cell type: "
            "Tumor cells`: SK-N-SH IS a neuroblastoma-derived line, which is what the reader's "
            "own covered `SK-N-SH` is, so the `tissue` field describes the line's origin and not "
            "a tumour biopsy."
        ),
    },
)

#: Named candidates that FAIL a clause. Recorded because this is where a lift would have had to
#: come from, and it does not come from any of them.
REJECTED: tuple[dict[str, Any], ...] = (
    {
        "cell": "K562, GM12878, H1",
        "archive": "GEO",
        "accession": "GSE126074 (SNARE-seq), file GSE126074_CellLineMixture_SNAREseq_cDNA_counts.tsv.gz",
        "clause_failed": "biosample identity",
        "why": (
            "THIS IS THE CLOSEST THING IN GEO TO SINGLE-NUCLEUS RNA OF THE COVERED CELL LINES AND "
            "IT CANNOT BE ADMITTED. SNARE-seq profiled 1,047 nuclei of a FOUR-LINE MIXTURE -- BJ, "
            "H1, K562 and GM12878 -- and the archive publishes ONE cDNA count table for the "
            "mixture, with no per-line split declared. Taking the K562 (or H1, or GM12878) "
            "barcodes out of it needs a per-barcode cell-line assignment the archive does not "
            "make, which is exactly the cell-identity decision the adopted rule does not license. "
            "The contrast is lane-scsweep's own K562 candidate, GSM9224207, which WAS admitted "
            "from a human/mouse mixture precisely BECAUSE the submitter split the supplementary "
            "files per line. Nuclei isolation is not a treatment and nothing here is modified: "
            "the biosample is a mixture of four cell lines, and that is the whole reason."
        ),
    },
    {
        "cell": "GM12878",
        "archive": "GEO",
        "accession": "GSE157660 (SNARE-seq2), GSM4772950-61",
        "clause_failed": "not treated",
        "why": (
            "A single-nucleus RNA library of GM12878 alone, which is what this lane needed, and "
            "the series' own design reads `GM12878 cells were freshly fixed with 1% "
            "formaldehyde` and the nuclei `fixed with 1% formaldehyde` before the protocol. "
            "Fixation is a treatment of the material, which is the clause on which lane-scsweep "
            "rejected the PFA-fixed half of its own K562 sample (`12b0c4e`). The same clause is "
            "applied here to the cell type where it costs this lane the most."
        ),
    },
    {
        "cell": "cardiac_muscle_cell",
        "archive": "GEO",
        "accession": "GSE129096 (Drop-seq and DroNc-seq during cardiomyocyte differentiation)",
        "clause_failed": "biosample identity, and the exact label",
        "why": (
            "DroNc-seq IS single-nucleus RNA, and the biosample is an hiPSC culture at one of five "
            "time points of a differentiation, in which the submitters themselves report SIX "
            "distinct cell types. The cardiomyocytes are an annotated subset inside that culture, "
            "so admitting it means the cell-composition decision the adopted rule does not "
            "license. Its label, `cardiomyocyte`, is also not the identical string to ENCODE's "
            "`cardiac muscle cell`. Recorded in full because ENCODE's own covered `cardiac muscle "
            "cell` is itself `in vitro differentiated cells` originated from RUES2, so the "
            "rejection rests on the composition decision and on the label, NOT on a claim that a "
            "differentiated cardiomyocyte can never be that biosample."
        ),
    },
    {
        "cell": "K562",
        "archive": "CELLxGENE",
        "accession": "38a78e91-e361-45c8-bae3-1d6d5108d18e (sci-Plex, tissue label `K-562`)",
        "clause_failed": "not treated, and the exact label",
        "why": (
            "The only CELLxGENE dataset with `suspension_type: nucleus` whose tissue label is a "
            "covered cell line, and it is 'Massively multiplex chemical transcriptomics at "
            "single-cell resolution': a chemical perturbation screen in which every arm receives "
            "a compound or a vehicle, so there is no untreated arm to admit. lane-scsweep "
            "rejected it on the same clause (`12b0c4e`). Its label `K-562` is also not the "
            "identical string to the covered `K562`."
        ),
    },
    {
        "cell": "hepatocyte, astrocyte, cardiac_muscle_cell, CD14_positive_monocyte",
        "archive": "CELLxGENE and GEO",
        "accession": "e.g. GSE223581, GSE212046, GSE189175 (liver snRNA-seq); GSE242197, "
        "GSE268599, GSE180759 (brain snRNA-seq); GSE255992, GSE271229, GSE165838 (heart "
        "snRNA-seq); GSE277496 (blood monocytes)",
        "clause_failed": "biosample identity",
        "why": (
            "SINGLE-NUCLEUS RNA OF THESE CELL TYPES IS PUBLISHED IN QUANTITY AND NONE OF IT IS "
            "ADMISSIBLE HERE. In every case the biosample is a TISSUE -- liver, brain, heart, "
            "blood -- and the covered cell type is a cell-type annotation of nuclei inside it. "
            "Admitting one means selecting that annotated subset, which is the cell-composition "
            "decision the adopted rule does not license, and it is the same ground on which "
            "lane-pseudobulk excluded Tabula Sapiens and the HCA and lane-scsweep excluded these "
            "very five labels. No human CELLxGENE dataset carries any of them as its `tissue`. "
            "Several of the named series are also disease cohorts (MASH, HCC, Huntington's, "
            "multiple sclerosis, atrial fibrillation, Crohn's), which is a second failure."
        ),
    },
    {
        "cell": "testis",
        "archive": "ENCODE, CELLxGENE, GEO",
        "accession": "none admitted",
        "clause_failed": "no single-nucleus experiment on this biosample",
        "why": (
            "THE CELL TYPE WHOSE SINGLE-CELL ARM lane-pbrun USED, AND IT HAS NO SINGLE-NUCLEUS "
            "ARM. ENCODE's released human snRNA-seq facet has no `testis` term. Every CELLxGENE "
            "human dataset whose tissue is `testis` -- Tabula Sapiens - Testis included, which is "
            "the one lane-pbrun read -- is `suspension_type: cell`. The three GEO series these "
            "queries surfaced are a pediatric solid-tumour atlas, a U7 snRNA study and a "
            "replication-timing series, none of them single-nucleus RNA of a normal testis."
        ),
    },
    {
        "cell": "HepG2, IMR_90, keratinocyte",
        "archive": "GEO, ENCODE, CELLxGENE",
        "accession": "none admitted; the GEO series surfaced were read",
        "clause_failed": "no single-nucleus experiment on this biosample, or not modified/treated",
        "why": (
            "No released human snRNA-seq exists on these terms in ENCODE and no human CELLxGENE "
            "nucleus dataset carries them as its tissue. In GEO the series these queries surfaced "
            "are tissue or disease studies that merely mention the line (MASH models, "
            "senescence and fibrosis cohorts, the ENCODE ChIP-seq and replication-timing "
            "deposits), or sample-multiplexing protocols on mixtures (GSE165047, CASB, a mixture "
            "of 16 cell types labelled with concanavalin A). THAT IS A STATEMENT ABOUT THESE "
            "QUERIES AND NOT ABOUT GEO: an archive listing is what is published, not what is "
            "usable, and this lane does not claim no such submission exists."
        ),
    },
)


def covered_terms() -> dict[str, str]:
    """The thirteen covered cell keys mapped to their ENCODE biosample_term_name.

    Read through lane-pseudobulk's own loader, so the covered set is the one `2b9cd0a` pinned and
    is neither re-derived nor widened here.
    """
    return {c["cell"]: c[pb.MATCH_FIELD] for c in pb.covered_cells()}


def predecessor_single_nucleus_cells() -> list[str]:
    """Covered cell types lane-pseudobulk's committed inventory admits on an `snRNA-seq` experiment.

    Read from that committed file, and filtered by the MODALITY clause: a cell type it admitted on
    `scRNA-seq` or `long read scRNA-seq` alone is single-CELL and is NOT returned. This is where
    GM12878 drops out of the count that lane-pbrun used.
    """
    record = json.loads(PSEUDOBULK_INVENTORY.read_text())
    out: list[str] = []
    for entry in record["per_cell"]:
        titles = {a.get("assay_title") for a in entry["single_cell"]["admitted"]}
        if SINGLE_NUCLEUS_ASSAY in titles:
            out.append(entry["cell"])
    return sorted(out)


def predecessor_single_cell_only_cells() -> list[str]:
    """The cell types the two predecessor inventories count that are single-CELL ONLY.

    Named in the record so the difference between this lane's count and lane-pbrun's four is
    visible as a modality exclusion and not as a disagreement about a clause.
    """
    counted = set(json.loads(SCSWEEP_INVENTORY.read_text())["exact_label_cells"])
    return sorted(counted - set(predecessor_single_nucleus_cells()))


def is_single_nucleus_encode(assay_title: str) -> bool:
    """ENCODE's modality, from its own `assay_title` and nothing else."""
    return assay_title == SINGLE_NUCLEUS_ASSAY


def is_single_nucleus_cellxgene(suspension_types: list[str] | tuple[str, ...]) -> bool:
    """CELLxGENE's modality, from the schema's own `suspension_type` field and nothing else."""
    return "nucleus" in set(suspension_types)


def is_single_nucleus_geo(library_strategy: str, words: str) -> bool:
    """GEO's modality, from the sample's own words, with the ATAC arm of a multiome excluded.

    `words` is the sample's own free text -- its extract protocol, title and series design joined.
    The bare word `nuclei` in a series summary is NOT enough on its own, by MODALITY_POLICY: it
    matches nuclear extraction and nuclear RNA in bulk assays.
    """
    if library_strategy.strip().upper() not in ("RNA-SEQ", "RNA-SEQ (SINGLE NUCLEUS)"):
        return False
    low = words.lower()
    if "whole cell" in low or "whole cells" in low:
        return False
    return any(
        term in low
        for term in (
            "single nuclei",
            "single nucleus",
            "single-nucleus",
            "single-nuclei",
            "snrna-seq",
            "snrna seq",
            "dronc-seq",
            "snpatho-seq",
            "isolated nuclei",
        )
    )


def verdict(exact_label_count: int) -> str:
    """The cell-type half of the adopted minimum, read on SINGLE-NUCLEUS cell types alone.

    `NO_SOURCE_QUALIFIES` is the adopted rule's own complete answer below 4. At or above 4 this
    returns `CELL_TYPE_MINIMUM_REACHED` and NOT `QUALIFIES`: the adopted rule's `QUALIFIES`
    requires BOTH numbers, and the second -- 1,000 of M2h's held-out genes resolving in every
    qualifying cell type in both modalities -- cannot be counted without reading the matrices,
    which this lane does not do. `pb.verdict` is called beside this to show that it agrees.
    """
    if exact_label_count < MINIMUM["cell_types"]:
        return "NO_SOURCE_QUALIFIES"
    return "CELL_TYPE_MINIMUM_REACHED"


#: What `NO_SOURCE_QUALIFIES` means for Albert's question, in the words the lane was given. This
#: is a CLOSING answer and a complete one, not a failure to find data.
CLOSES_WITHOUT_A_LIFT = (
    "The modality swap lane-pbrun licensed (`22cc987`) is licensed FOR SINGLE-CELL PSEUDOBULK "
    "ONLY on M2h's endpoint: three of its four cell types are single-CELL and only ovary is "
    "single-NUCLEUS, and n = 1 cannot run the registered within-gene permutation. This lane "
    "searched ENCODE, CELLxGENE and GEO for single-NUCLEUS RNA of the covered cell types under "
    "`2b9cd0a`'s rule imported unchanged, and found fewer than the minimum of 4. So the "
    "single-nucleus-only comparison CANNOT BE REGISTERED OR RUN, bile duct's single-nucleus test "
    "stays UNINTERPRETABLE -- it would change the cell type and the assay modality at once -- and "
    "THE BIOSAMPLE QUESTION CLOSES WITHOUT A LIFT. The minimum was not stretched, no single-CELL "
    "experiment was counted towards it, and no weaker test was proposed to rescue it."
)
