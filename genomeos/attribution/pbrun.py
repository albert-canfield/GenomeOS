# SPDX-License-Identifier: AGPL-3.0-or-later
"""Pseudobulked single-cell/single-nucleus RNA against bulk RNA-seq, on M2h's endpoint.

WHY THIS LANE EXISTS. M2h (`ac5580b`) found the compiled rules fit measured RNA better when scaled
by one fitted constant -- held-out Delta MSE +0.04339 [98.33% CI +0.0379, +0.04894], scale +0.2146,
about 2% of variance -- M2i (`f6e71c6`) replicated it on RNA fresh to this project with nothing
fitted (+0.03727 [+0.03202, +0.04299]), and lane-scalesplit (`ec9cfd4`) showed it survives a
chromosome-held-out split with the scale refitted on the training chromosomes (both folds better,
+0.04478 and +0.04201). Its ONE unresolved limit is that independence from AlphaGenome's training
domain is UNKNOWN FOR 13 OF 13 (`18de375`). The only route past that is a cell type OUTSIDE that
domain, and lane-oodinventory (`171f5bf`) found exactly one -- bile duct, UBERON:0002394 -- whose
RNA is SINGLE-NUCLEUS ONLY. The supervisor ruled a test there uninterpretable while the modality
swap is unlicensed, because it would change the cell type AND the assay at once.

THIS LANE LICENSES, OR REFUSES TO LICENSE, THAT SWAP. lane-pseudobulk (`2b9cd0a`) wrote the rule and
its minimum -- 4 cell types and 1,000 genes -- before counting anything; lane-scsweep (`12b0c4e`)
reached 4 of 4 under that rule verbatim (GM12878 and ovary from ENCODE, K562 from GEO, testis from
CELLxGENE) and returned `CELL_TYPE_MINIMUM_REACHED`, NOT `QUALIFIES`, because the 1,000-gene clause
was uncounted. Counting it is this module's first job, and if it fails the lane ends at
`NO_SOURCE_QUALIFIES`: the modality swap stays unlicensed, the bile-duct test stays uninterpretable,
and Albert's biosample question closes without a lift. That is a complete answer.

NOTHING IS FITTED HERE. M2h's scale enters only as the fixed constant +0.2146 of `ac5580b`, by
import from lane-pseudobulk's module, and `MINIMUM` is lane-pseudobulk's dict itself rather than a
copy of its numbers.

WORDING. Agreement on this endpoint licenses the modality swap FOR THIS ENDPOINT ONLY. It is never
"single-nucleus RNA equals bulk RNA".
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import tarfile
from pathlib import Path
from typing import Any

import h5py
import numpy as np

from genomeos.attribution import pseudobulk as pb

ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = ROOT / "data" / "results"
CACHE_DIR = Path("data/cache/pbrun")

#: lane-scsweep's committed record: the four cell types whose archive label is the IDENTICAL STRING
#: to ENCODE's biosample_term_name. Read, never re-derived and never widened.
SCSWEEP_RESULT = RESULTS_DIR / "scsweep_inventory.json"
#: lane-pseudobulk's committed inventory: the ENCODE single-cell experiments it admitted, and the
#: bulk experiment and files per covered cell type.
PSEUDOBULK_INVENTORY = RESULTS_DIR / "pseudobulk_inventory.json"
#: 7h's committed genes table: 16,509 genes, each with `symbol`, `gene_id`, thirteen P and thirteen D.
MULTRNA_RESULT = RESULTS_DIR / "multrna_result.json"
#: M2h's committed record, for the fixed scale and for the reading this lane's wording is read against.
MAGCAL_RESULT = RESULTS_DIR / "magcal_result.json"
#: M2i's committed records, for the bulk files on disk and the thirteen covered cell keys.
FRESHRNA_RESULT = pb.FRESHRNA_RESULT
FRESHRNA_SOURCES = pb.FRESHRNA_SOURCES

GENE_COUNT_PATH = RESULTS_DIR / "pbrun_gene_count.json"
REGISTRATION_PATH = RESULTS_DIR / "pbrun_registration.json"
RESULT_PATH = RESULTS_DIR / "pbrun_result.json"

#: M2h's fitted scale, as the FIXED CONSTANT of `ac5580b` and never refitted here. Imported from
#: lane-pseudobulk so the two lanes cannot drift apart.
SCALE = pb.SCALE

#: lane-pseudobulk's minimum, IMPORTED and not copied: 4 cell types and 1,000 genes.
MINIMUM = pb.MINIMUM

#: lane-pseudobulk's gene clause, carried word for word as the thing this lane counts.
GENE_CLAUSE = pb.MINIMUM["genes_clause"]

WORDING_LIMIT = pb.WORDING_LIMIT

#: The ONE form of M2h whose reading is "better" and whose scale is +0.2146 (the mean form is
#: inconclusive and the promoter-competition form is a negative), so it is the form whose held-out
#: population this lane counts and compares. Named, not chosen after a value was seen.
FORM = "multiplicative"

#: Which identifier each single-cell matrix carries, read from the file itself before any count.
#: A gene is matched to a matrix on the identifier THAT MATRIX CARRIES: the unversioned Ensembl gene
#: id where it carries Ensembl ids, the gene symbol where it carries symbols. Both come from the
#: same row of 7h's committed genes table, so no new mapping is introduced.
ENSEMBL = "unversioned_ensembl_gene_id"
SYMBOL = "gene_symbol"

SC_SELECTION: dict[str, Any] = {
    "version": "pbrun-1",
    "principle": (
        "One single-cell or single-nucleus source per qualifying cell type, named from the committed "
        "inventories and chosen on METADATA ALONE before any value was read. No source is re-searched "
        "and no cell type is added: the four are lane-scsweep's `exact_label_cells` of `12b0c4e`."
    ),
    "encode": (
        "For a cell type ENCODE supplies (GM12878, ovary), EVERY experiment lane-pseudobulk's "
        "committed inventory admitted for that biosample_term_name contributes, and from each "
        "experiment one file is taken by a fixed order: assembly GRCh38 over any other assembly, "
        "then the highest genome_annotation version present, then output_type by the priority "
        "'sparse gene count matrix of unique reads' (the called-cell STARsolo matrix of unique reads) "
        "> 'gene quantifications' (the per-cell RSEM table), ties broken by the lowest file "
        "accession. For GM12878 ENCODE's scRNA-seq deposits ONE CELL PER FILE, so every qualifying "
        "file of each admitted experiment contributes as one cell."
    ),
    "geo_and_cellxgene": (
        "For K562 and testis exactly one file is admitted by lane-scsweep and no choice arises: "
        "GSM9224207_Fresh_K562.h5ad (the FRESH unfixed preparation; the PFA-fixed file of the same "
        "sample is excluded as a treatment, which is lane-scsweep's judgement carried) and CELLxGENE "
        "dataset ec6ef004-5f30-4a40-b282-980cca20a561 (Tabula Sapiens - Testis), whose raw counts "
        "the archive declares at `raw.X`."
    ),
    "counts_layer": (
        "Integer counts and nothing else: the STARsolo matrix.mtx for ENCODE snRNA-seq, the "
        "`expected_count` column of the per-cell RSEM table for ENCODE scRNA-seq, `X` for the GEO "
        "h5ad (which carries no raw layer and whose X is integral), and `raw.X` for the CELLxGENE "
        "h5ad, which is where that archive's schema declares raw counts. No normalised, imputed, "
        "decontaminated or scaled layer is read: not `layers/Q-Q_normalized` of the GEO file, and "
        "not `X`, `layers/decontXcounts` or `layers/scale_data` of the CELLxGENE file."
    ),
    "modality_note": (
        "OF THE FOUR, ONE IS SINGLE-NUCLEUS: the ovary snRNA-seq. GM12878, K562 and testis are "
        "single-CELL. This is stated here, before any value is read, because the bile-duct source is "
        "single-NUCLEUS: an agreement read here licenses pseudobulked single-cell-or-nucleus RNA "
        "against bulk ON THIS ENDPOINT, and the nucleus-specific part of it rests on one of the four "
        "cell types. The per-cell-type descriptives make the ovary arm readable on its own and "
        "decide nothing."
    ),
}

AGGREGATION: dict[str, Any] = {
    "statement": (
        "SUM THE COUNTS, THEN CPM, THEN LOG2. (1) For a cell type, every contributing file's integer "
        "count matrix is read (SC_SELECTION['counts_layer']) and the counts of a gene are SUMMED over "
        "ALL of that cell type's cells or nuclei and over all of its contributing files, giving one "
        "integer per gene per cell type. (2) That summed vector is scaled to counts per million: "
        "cpm(G, X) = 1e6 * c(G, X) / sum over all genes of c(., X), the sum taken over EVERY row of "
        "the matrix and not only over the counted genes. (3) The measured value is "
        "y_sc(G, X) = log2(cpm(G, X) + 1), with PSEUDOCOUNT 1.0 carried verbatim from M2h."
    ),
    "no_filter": (
        "NO FILTER AT ANY STEP. No expression threshold, no minimum count, no minimum number of "
        "cells, no doublet removal, no mitochondrial-fraction cut, no cell-type subsetting inside a "
        "tissue, no highly-variable-gene selection and no removal of a donor or a batch. Every "
        "barcode present in the chosen matrix contributes and every gene row is kept. This is "
        "lane-pseudobulk's 'no filter on expression and no filter on any measured value' applied to "
        "the single-cell side."
    ),
    "why_sum": (
        "Summing counts over cells is the aggregation that makes a pseudobulk stand for the library a "
        "bulk experiment would have sequenced from the same material, so it is the comparison the "
        "modality question asks. A mean of per-cell normalised values would weight a 1,000-UMI cell "
        "equally with a 20,000-UMI cell and would not be that library; it is not read here."
    ),
    "cpm_not_tpm": (
        "The single-cell side is CPM and the bulk side is RSEM TPM, which is length-normalised. The "
        "difference is a PER-GENE factor, 1/length, and the endpoint centres each gene within itself "
        "across the cell types, so a per-gene factor cancels from both D's up to the pseudocount. "
        "What does NOT cancel is a per-CELL-TYPE term (the ratio of the two normalising sums), which "
        "shifts all of a cell type's genes together; it is declared as a limitation and no "
        "length-normalised variant of the single-cell side was computed or read."
    ),
}

#: The four cell types, in the order lane-scsweep's record lists them. Read from the record.
CELLS_SOURCE = (
    "data/results/scsweep_inventory.json, field `exact_label_cells`, pinned by sha256 in the record"
)


def sha256_of(path: Path) -> str:
    """The sha256 of a file on disk, lane-pseudobulk's helper carried."""
    return pb.sha256_of(path)


def qualifying_cells() -> list[str]:
    """lane-scsweep's four cell types, read from its committed record and never re-derived."""
    record = json.loads(SCSWEEP_RESULT.read_text())
    cells = record["exact_label_cells"]
    return [c["cell"] if isinstance(c, dict) else str(c) for c in cells]


def is_held_out(symbol: str) -> bool:
    """M2h's committed split, carried unchanged: odd first 32 bits of sha256(symbol) are held out."""
    return int(hashlib.sha256(symbol.encode("utf-8")).hexdigest()[:8], 16) % 2 == 1


def unversioned(gene_id: str) -> str:
    """An Ensembl gene id without its version suffix. ``_PAR_Y`` rows keep their suffix and so are
    never equal to a plain id, which is how lane-pseudobulk's clause excludes them."""
    if "_PAR_Y" in gene_id:
        return gene_id
    return gene_id.split(".", 1)[0]


def held_out_population() -> list[dict[str, str]]:
    """M2h's held-out genes of the ``multiplicative`` form: symbol, unversioned Ensembl id.

    7h's committed genes table filtered by M2h's committed split on the gene SYMBOL,
    ``int(sha256(symbol as UTF-8).hexdigest()[:8], 16)`` odd, carried unchanged from `ac5580b` by way
    of lane-pseudobulk. No gene is admitted that 7h did not admit and none it did is dropped.
    """
    table = json.loads(MULTRNA_RESULT.read_text())["genes_table"]
    columns = table["columns"]
    i_symbol, i_gene_id = columns.index("symbol"), columns.index("gene_id")
    out: list[dict[str, str]] = []
    for row in table["rows"]:
        symbol = str(row[i_symbol])
        if not is_held_out(symbol):
            continue
        out.append({"symbol": symbol, "gene_id": unversioned(str(row[i_gene_id]))})
    return out


def row_counts_rsem_tsv(path: Path) -> dict[str, int]:
    """How many rows an ENCODE RSEM table has per unversioned gene id. Spike-ins are not Ensembl
    ids and simply do not match any gene of the population."""
    counts: dict[str, int] = {}
    with path.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        i = header.index("gene_id")
        for line in handle:
            gene_id = line.split("\t", i + 1)[i]
            key = unversioned(gene_id)
            counts[key] = counts.get(key, 0) + 1
    return counts


def _one_member(names: list[str], suffix: str) -> str:
    """The one archive member whose name ends in `suffix`, or a raise. Never a guess among several."""
    hits = [n for n in names if n.endswith(suffix)]
    if len(hits) != 1:
        raise ValueError(f"expected exactly one {suffix} in the archive, found {len(hits)}: {hits}")
    return hits[0]


def row_counts_starsolo_tar(path: Path) -> dict[str, int]:
    """How many rows a STARsolo ``features.tsv`` inside an ENCODE tar has per unversioned gene id."""
    counts: dict[str, int] = {}
    with tarfile.open(path, "r:gz") as tar:
        name = _one_member(tar.getnames(), "features.tsv")
        member = tar.extractfile(name)
        if member is None:
            raise ValueError(f"{name} in {path} is not a readable file")
        for raw in io.TextIOWrapper(member, encoding="utf-8"):
            gene_id = raw.split("\t", 1)[0].strip()
            if not gene_id:
                continue
            key = unversioned(gene_id)
            counts[key] = counts.get(key, 0) + 1
    return counts


def row_counts_h5ad(path: Path, index: str) -> dict[str, int]:
    """How many rows an h5ad's var index has per identifier. `index` is the HDF5 path of that index.

    Read with h5py and not with anndata: only the var index is opened, so no matrix is loaded and no
    value is read.
    """
    counts: dict[str, int] = {}
    with h5py.File(path, "r") as handle:
        for raw in handle[index][:]:
            key = raw.decode() if isinstance(raw, bytes) else str(raw)
            key = unversioned(key) if key.startswith("ENSG") else key
            counts[key] = counts.get(key, 0) + 1
    return counts


def resolves(counts_per_file: list[dict[str, int]], key: str) -> bool:
    """lane-pseudobulk's clause on one side: EXACTLY ONE row for the key in EVERY contributing file."""
    return bool(counts_per_file) and all(c.get(key, 0) == 1 for c in counts_per_file)


def verdict(qualifying_cell_count: int, qualifying_gene_count: int | None) -> str:
    """lane-pseudobulk's own verdict function, imported and not reimplemented."""
    return pb.verdict(qualifying_cell_count, qualifying_gene_count)


def pseudobulk_cpm(counts: np.ndarray) -> np.ndarray:
    """AGGREGATION steps (2) and (3) on one cell type's summed integer counts.

    `counts` is the summed count per gene over every row of the matrix. The library size is the sum
    over EVERY row, not over the counted genes, exactly as AGGREGATION states.
    """
    total = float(counts.sum())
    if total <= 0:
        raise ValueError("a cell type's summed counts are empty, which is a defect and not a value")
    return 1e6 * counts.astype(np.float64) / total


def log2p1(values: np.ndarray) -> np.ndarray:
    """M2h's measured transform, pseudocount 1.0 carried verbatim."""
    return np.log2(values + 1.0)


def centred(matrix: np.ndarray) -> np.ndarray:
    """D(G, X) = y(G, X) minus the mean of y over the cell types in `matrix`'s columns.

    The one deviation from M2h, applied IDENTICALLY to both arms and to P: M2h centres over its
    thirteen cell types and this comparison has four, so every quantity is centred over the same
    four columns. Centring one arm over thirteen and the other over four would make the two numbers
    incomparable, which is the error this function exists to prevent.
    """
    return matrix - matrix.mean(axis=1, keepdims=True)


def delta_mse(d: np.ndarray, p: np.ndarray, scale: float) -> float:
    """M2h's metric at a fixed scale: mean over pairs of 2*a*D*P - a^2*P^2. Nothing is fitted."""
    return float(np.mean(2.0 * scale * d * p - (scale**2) * (p**2)))


def open_gzip_or_plain(path: Path) -> io.TextIOBase:
    """A text handle on a file that may or may not be gzipped, by its magic bytes and not its name."""
    with path.open("rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return path.open(encoding="utf-8")
