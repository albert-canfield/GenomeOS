# SPDX-License-Identifier: AGPL-3.0-or-later
"""Simpler genomes, the fly: fetch Calderon et al. 2022 and distil the RNA to per-annotation sums.

    uv run --no-project --with pyreadr --with pandas --with numpy python scripts/simpler_fly_fetch.py

Runs outside the project environment because the source's cell-to-annotation table is an R data
file and pyreadr is not a project dependency; it imports nothing from genomeos for that reason.

What it does, in order, all into the git-ignored store data/knowledge/simpler_fly:

1. Fetches, if absent: the source's supplementary tables (ATAC cluster annotations, Table S2; ATAC
   marker accessible regions, Table S4), the full peak set, the RNA cell-to-annotation table, the
   GEO RNA count matrix with its row and column names, and the Ensembl 104 (BDGP6.32) gene models.
2. Sums the RNA count matrix over cells by the source's own fine annotation (`manual_annot` in
   rna_meta.rds), the way the source's own script does (make_germ_pseudobulk.R), and writes one
   gene-by-annotation table of integer UMI sums. No tissue grouping is applied here: the grouping is
   registered in genomeos/attribution/simpler_fly.py.
3. Deletes the 698 MB matrix once summed and writes a record of every URL, byte count and sha256.

It prints byte counts, hashes and cells per annotation only. It prints no gene's count: the RNA is
read by `scripts/simpler_fly_run.py score`, after the registration is committed.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STORE = ROOT / "data/knowledge/simpler_fly"
DEAP = "https://shendure-web.gs.washington.edu/content/members/DEAP_website/public"
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190147/suppl"
ENSEMBL = (
    "https://ftp.ensembl.org/pub/release-104/gff3/drosophila_melanogaster/"
    "Drosophila_melanogaster.BDGP6.32.104.gff3.gz"
)

#: name in the store -> (URL, what it is, kept after distillation)
FILES = {
    "SupplementaryTables.xlsx": (
        f"{DEAP}/SupplementaryTables.xlsx",
        "Calderon 2022 supplementary tables S1-S13; read: TableS2 (ATAC cluster annotations) and "
        "TableS4 (ATAC marker accessible regions)",
        True,
    ),
    "peaks.bed": (f"{DEAP}/ATAC/revision/peaks.bed", "the full ATAC peak set, dm6", True),
    "rna_meta.rds": (
        f"{DEAP}/cell_to_annotations/rna_meta.rds",
        "RNA cell identifiers to the annotations used in the source's Fig. 3 (manual_annot)",
        True,
    ),
    "GSE190147_scirnaseq_gene_matrix.rows.txt.gz": (
        f"{GEO}/GSE190147_scirnaseq_gene_matrix.rows.txt.gz",
        "RNA matrix row names (genes)",
        True,
    ),
    "GSE190147_scirnaseq_gene_matrix.columns.txt.gz": (
        f"{GEO}/GSE190147_scirnaseq_gene_matrix.columns.txt.gz",
        "RNA matrix column names (cells) with per-cell QC fields",
        True,
    ),
    "GSE190147_scirnaseq_gene_matrix.mtx.gz": (
        f"{GEO}/GSE190147_scirnaseq_gene_matrix.mtx.gz",
        "RNA UMI count matrix, genes x cells (MatrixMarket); summed by annotation, then deleted",
        False,
    ),
    "Drosophila_melanogaster.BDGP6.32.104.gff3.gz": (
        ENSEMBL,
        "Ensembl 104 gene models for dm6 (BDGP6.32); its 13,968 protein-coding gene names all "
        "appear among the RNA matrix's row names",
        True,
    ),
}
PSEUDOBULK = "rna_umi_by_manual_annot.tsv.gz"
RECORD = "fetch_record.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(name: str, url: str) -> None:
    dest = STORE / name
    if dest.exists():
        return
    part = dest.with_suffix(dest.suffix + ".part")
    subprocess.run(
        ["curl", "-sS", "-L", "--retry", "8", "--retry-delay", "15", "-C", "-", "-o", str(part), url],
        check=True,
    )
    part.rename(dest)


def distil(record: dict) -> None:
    import numpy as np
    import pandas as pd
    import pyreadr

    meta = pyreadr.read_r(str(STORE / "rna_meta.rds"))[None]
    label_of = dict(zip(meta["cell"].astype(str), meta["manual_annot"].astype(str), strict=True))
    labels = sorted(set(label_of.values()))
    li = {lab: i for i, lab in enumerate(labels)}
    with gzip.open(STORE / "GSE190147_scirnaseq_gene_matrix.columns.txt.gz", "rt") as f:
        head = f.readline().rstrip("\n").split("\t")
        ci = head.index("cell")
        cells = [line.rstrip("\n").split("\t")[ci] for line in f]
    with gzip.open(STORE / "GSE190147_scirnaseq_gene_matrix.rows.txt.gz", "rt") as f:
        genes = [line.rstrip("\n") for line in f]
    col_label = np.array([li.get(label_of.get(c, ""), -1) for c in cells], dtype=np.int64)
    n_lab = len(labels)
    mtx = STORE / "GSE190147_scirnaseq_gene_matrix.mtx.gz"
    with gzip.open(mtx, "rt") as f:
        skip = 0
        line = f.readline()
        while line.startswith("%"):
            skip += 1
            line = f.readline()
        n_rows, n_cols, nnz = (int(x) for x in line.split())
        skip += 1
    if (n_rows, n_cols) != (len(genes), len(cells)):
        sys.exit(f"refused: matrix is {n_rows} x {n_cols}, names are {len(genes)} x {len(cells)}")
    acc = np.zeros(n_rows * n_lab, dtype=np.float64)
    seen = 0
    for ch in pd.read_csv(
        mtx,
        sep=" ",
        header=None,
        names=["r", "c", "v"],
        skiprows=skip,
        compression="gzip",
        dtype={"r": np.int64, "c": np.int64, "v": np.float64},
        chunksize=20_000_000,
    ):
        seen += len(ch)
        lab = col_label[ch["c"].to_numpy() - 1]
        keep = lab >= 0
        idx = (ch["r"].to_numpy()[keep] - 1) * n_lab + lab[keep]
        acc += np.bincount(idx, weights=ch["v"].to_numpy()[keep], minlength=n_rows * n_lab)
    if seen != nnz:
        sys.exit(f"refused: read {seen} entries, the header says {nnz}")
    acc = acc.reshape(n_rows, n_lab)
    if not np.all(acc == np.round(acc)):
        sys.exit("refused: the sums are not integers, so the matrix is not raw UMI counts")
    out = STORE / PSEUDOBULK
    with gzip.open(out, "wt") as f:
        f.write("gene\t" + "\t".join(labels) + "\n")
        for g, row in zip(genes, acc.astype(np.int64), strict=True):
            f.write(g + "\t" + "\t".join(str(int(x)) for x in row) + "\n")
    matched = int((col_label >= 0).sum())
    cells_per = {lab: 0 for lab in labels}
    for x in col_label[col_label >= 0]:
        cells_per[labels[x]] += 1
    record["distillation"] = {
        "what": (
            "UMI counts summed over cells by the source's manual_annot (rna_meta.rds), all windows "
            "pooled; one row per matrix row (gene name), one column per annotation"
        ),
        "file": f"data/knowledge/simpler_fly/{PSEUDOBULK}",
        "bytes": out.stat().st_size,
        "sha256": sha256(out),
        "matrix_shape": [n_rows, n_cols],
        "matrix_nonzero_entries": nnz,
        "cells_in_matrix": len(cells),
        "cells_with_an_annotation": matched,
        "annotations": len(labels),
        "cells_per_annotation": cells_per,
        "pyreadr": pyreadr.__version__ if hasattr(pyreadr, "__version__") else "not reported",
        "method_follows": f"{DEAP}/scripts/RNA/misc/make_germ_pseudobulk.R (sums counts by manual_annot)",
    }


def main() -> None:
    STORE.mkdir(parents=True, exist_ok=True)
    rec_path = STORE / RECORD
    record = json.loads(rec_path.read_text()) if rec_path.exists() else {"files": {}}
    for name, (url, what, kept) in FILES.items():
        if name in record["files"] and not (STORE / name).exists() and not kept:
            continue  # summed and deleted on an earlier run
        fetch(name, url)
        p = STORE / name
        record["files"][name] = {
            "url": url,
            "what": what,
            "bytes": p.stat().st_size,
            "sha256": sha256(p),
            "kept": kept,
        }
        print(name, record["files"][name]["bytes"], record["files"][name]["sha256"])
    if not (STORE / PSEUDOBULK).exists():
        distil(record)
    mtx = STORE / "GSE190147_scirnaseq_gene_matrix.mtx.gz"
    if mtx.exists():
        mtx.unlink()
    record["download_bytes"] = sum(f["bytes"] for f in record["files"].values())
    rec_path.write_text(json.dumps(record, indent=1) + "\n")
    d = record["distillation"]
    print("distilled", d["file"], d["bytes"], d["sha256"])
    print(
        "cells", d["cells_in_matrix"], "annotated", d["cells_with_an_annotation"], "labels", d["annotations"]
    )
    print("download bytes", record["download_bytes"])


if __name__ == "__main__":
    main()
