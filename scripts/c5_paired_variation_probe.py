# SPDX-License-Identifier: AGPL-3.0-or-later
"""Item 13 C5, the probe: natural variation as experiments, feasibility only.

    uv run python scripts/c5_paired_variation_probe.py

Different human haplotypes are alternative versions of the same regulatory sequence. Where a genotype
and a molecular measurement are paired in the same samples (allele-specific expression, splicing,
accessibility or binding), a proposed explanation from the attribution layer ("this element regulates
that gene") predicts an allelic difference that could be checked. VCFs alone carry no outcome. This
script establishes, label-blind, how much of the attribution layer's element set the open paired
resources could reach. It measures nothing and registers nothing that measures.

What it counts, per resource (docs/DATA.md, "Natural variation as experiments: the C5 probe"):

- the element universe: every element of the all-element archive (data/knowledge/alphagenome/
  all_elements/chr*.json, 961,227 elements), and the element-gene pairs where the archive names a
  predicted coding target (`predicted_coding.gene`), which are the explanations a paired readout
  could check;
- per phased genome (NA12878 = GM12878 from GIAB HG001 v4.2.1; HG002 from the Q100 dipcall calls
  already cached): elements carrying at least one heterozygous SNV, target genes whose exons carry
  one (the readout side of allele-specific expression), pairs with both, and the heterozygous SNVs
  in the target's cis window (the linked-variant count);
- the Geuvadis samples on chr21 of the 1000 Genomes 30x phased panel (streamed, not stored): per
  element and per pair, how many samples are heterozygous;
- ADASTRA (allele-specific TF binding) and UDACHA (allele-specific accessibility): elements holding
  at least one of the SNVs each release lists as eligible (coverage-passing) in GM12878, K562, HepG2
  and IMR-90, the archive's four context cells;
- EN-TEx: elements and target genes that its catalogue lists as accessible (enough haplotype-resolved
  depth to be tested) per donor and tissue, and pairs accessible together in the same donor and
  tissue;
- the LCL context, since NA12878, HG002 and Geuvadis were all read from lymphoblastoid lines: pairs
  whose element overlaps a GM12878 H3K27ac peak (ENCFF361XMX, already cached) and whose target has a
  GTEx v8 median TPM of at least 1 in EBV-transformed lymphocytes. Neither is an allelic quantity.

Label-blind by construction: from the allelic tables only identifier, position, donor, tissue and
assay columns are read (`COLUMNS_READ`); count, ratio, effect-size and significance columns
(`COLUMNS_NEVER_READ`) are never indexed. Downloads stay under 1 GB each: four small tables and the
HG001 VCF and bed (151 MB) are cached in the git-ignored data/cache/c5/; the panel and the EN-TEx
tables are streamed and discarded, and ADASTRA and UDACHA are read by HTTP range, each with a sha256
of the bytes that passed. 0 AlphaGenome requests: the targets are read from the cached archive.
"""

from __future__ import annotations

import array
import csv
import gzip
import hashlib
import html
import io
import json
import math
import re
import sys
import time
import urllib.parse
import urllib.request
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from statistics import NormalDist
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "c5_paired_variation_probe"
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
ARCHIVE = Path("data/knowledge/alphagenome/all_elements")
CACHE = Path("data/cache/c5")
Q100 = Path("data/cache/q100")
HG002_VCF = Q100 / "GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.dip.vcf.gz"
HG002_BED = Q100 / "HG002_GRCh38_v5.0q_smvar.benchmark.bed"
GIAB = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/"
HG001_VCF_URL = GIAB + "HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
HG001_BED_URL = GIAB + "HG001_GRCh38_1_22_v4.2.1_benchmark.bed"
KGP = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/"
    "20220422_3202_phased_SNV_INDEL_SV/"
    "1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
)
KGP_CHROM = "chr21"  # the one chromosome streamed for the population count (407 MB)
GEUV_SDRF_URL = "https://ftp.ebi.ac.uk/biostudies/fire/E-GEUV-/001/E-GEUV-1/Files/E-GEUV-1.sdrf.txt"
ADASTRA_URL = "https://zenodo.org/records/14174114/files/ADASTRA.v.6.1.Mabel.zip?download=1"
UDACHA_KEY = "https://disk.yandex.ru/d/8Shsn1-OgDuwzQ"
ENTEX = "http://entex.encodeproject.org/data/"
GENEIMPRINT_URL = "https://www.geneimprint.com/site/genes-by-species.Homo+sapiens"
GTEX_TPM_URL = (
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/"
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
)
LCL_TISSUE = "Cells - EBV-transformed lymphocytes"
LCL_TPM_MIN = 1.0  # a target counts as expressed in LCLs at this GTEx median TPM or above
LCL_PEAKS = Path("data/knowledge/epigenome/peaks")  # ENCFF361XMX, GM12878 H3K27ac replicated peaks

#: the archive's four context cells, as each release names its cell-type file
_CELLS = (
    ("GM12878", "GM12878__female_B-cells_lymphoblastoid_cell_line_.tsv"),
    ("K562", "K562__myelogenous_leukemia_.tsv"),
    ("HepG2", "HepG2__hepatoblastoma_.tsv"),
    ("IMR-90", "IMR90__lung_fibroblasts_.tsv"),
)
CELL_FILES = {
    "ADASTRA": {cell: f"CL/{name}" for cell, name in _CELLS},
    "UDACHA": {f"{cell} {assay}": f"{assay}/{name}" for assay in ("dnase", "atac") for cell, name in _CELLS},
}
#: EN-TEx AS tables (cCREs_default_AS.tsv, genes_default_AS.tsv) share this header
ENTEX_HEADER = [
    "chr",
    "start",
    "end",
    "region_id",
    "hap1_count",
    "hap2_count",
    "experiment_accession",
    "donor",
    "tissue",
    "assay",
    "hap1_allele_ratio",
    "p_betabinom",
    "imbalance_significance",
]
COLUMNS_READ = {
    "EN-TEx": ["chr", "start", "end", "region_id", "experiment_accession", "donor", "tissue", "assay"],
    "ADASTRA": ["chr", "start"],
    "UDACHA": ["chr", "start"],
}
COLUMNS_NEVER_READ = {
    "EN-TEx": ["hap1_count", "hap2_count", "hap1_allele_ratio", "p_betabinom", "imbalance_significance"],
    "ADASTRA": ["every column after the position (BAD, coverage, effect sizes, FDR, motif)"],
    "UDACHA": ["every column after the position"],
}
CIS_WINDOW = 1_000_000  # bp either side of the target's TSS: the linked-variant count
GEUV_MIN_SAMPLES = (1, 10, 30)
GROUP_MIN = 10  # a between-sample contrast needs this many readable samples on each side
READ_DEPTH_EFFECTS = (0.55, 0.6, 0.65, 0.7)
UA = {"User-Agent": "genomeos-c5-probe"}


# ---------------------------------------------------------------- small I/O helpers


def fetch(url: str, dest: Path) -> Path:
    """Download a small file once into the git-ignored cache."""
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with (
        urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=300) as r,
        tmp.open("wb") as f,
    ):
        while chunk := r.read(1 << 20):
            f.write(chunk)
    tmp.rename(dest)
    return dest


class HashingReader(io.RawIOBase):
    """A streamed HTTP body, hashed and counted as it passes; nothing is written to disk."""

    def __init__(self, url: str) -> None:
        self.resp = urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=300)
        self.url = self.resp.geturl()
        self.sha = hashlib.sha256()
        self.bytes = 0

    def readable(self) -> bool:
        return True

    def readinto(self, b) -> int:
        data = self.resp.read(len(b))
        n = len(data)
        b[:n] = data
        self.sha.update(data)
        self.bytes += n
        return n

    def close(self) -> None:
        self.resp.close()
        super().close()

    def record(self, **extra: Any) -> dict[str, Any]:
        return {"url": self.url, "sha256": self.sha.hexdigest(), "bytes": self.bytes, **extra}


class RangeFile(io.RawIOBase):
    """A remote file read by HTTP ranges, so zipfile can pull single members out of a release archive."""

    def __init__(self, url: str) -> None:
        req = urllib.request.Request(url, headers={"Range": "bytes=0-0", **UA})
        with urllib.request.urlopen(req, timeout=120) as r:
            self.url = r.geturl()
            self.size = int(r.headers["Content-Range"].split("/")[-1])
        self.pos = 0
        self.fetched = 0

    def seekable(self) -> bool:
        return True

    def readable(self) -> bool:
        return True

    def tell(self) -> int:
        return self.pos

    def seek(self, off: int, whence: int = 0) -> int:
        self.pos = off if whence == 0 else self.pos + off if whence == 1 else self.size + off
        return self.pos

    def readinto(self, b) -> int:
        if self.pos >= self.size or not len(b):
            return 0
        end = min(self.size, self.pos + len(b)) - 1
        req = urllib.request.Request(self.url, headers={"Range": f"bytes={self.pos}-{end}", **UA})
        with urllib.request.urlopen(req, timeout=300) as r:
            data = r.read()
        b[: len(data)] = data
        self.pos += len(data)
        self.fetched += len(data)
        return len(data)


# ---------------------------------------------------------------- the element universe and genes


def load_elements() -> dict[str, dict[str, Any]]:
    """Per chromosome: element coordinates (0-based half-open, sorted by start), ids, and the predicted
    coding target where the archive names one."""
    out = {}
    for c in CHROMS:
        rows = json.loads((ARCHIVE / f"{c}.json").read_text())
        rows.sort(key=lambda e: e["start"])
        out[c] = {
            "start": np.array([e["start"] for e in rows], dtype=np.int64),
            "end": np.array([e["end"] for e in rows], dtype=np.int64),
            "id": [e["id"] for e in rows],
            "target": [(e.get("predicted_coding") or {}).get("gene") for e in rows],
            "strength": [(e.get("predicted_coding") or {}).get("strength") for e in rows],
        }
    return out


def load_genes(chrom: str) -> dict[str, dict[str, Any]]:
    """GENCODE v50 genes on one chromosome by symbol: stable id, TSS (0-based), and the union of every
    transcript's exons merged into disjoint 0-based half-open intervals."""
    path = Path("data/reference") / f"gencode_v50_{chrom}.gff3.gz"
    genes: dict[str, dict[str, Any]] = {}
    exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if len(f) < 9 or f[2] not in ("gene", "exon"):
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].rstrip("\n").split(";") if "=" in kv)
            name = attrs.get("gene_name")
            if not name:
                continue
            gid = attrs.get("gene_id", "").split(".")[0]
            s, e = int(f[3]) - 1, int(f[4])
            if f[2] == "gene":
                if name not in genes:  # a reused symbol keeps its first gene, as Annotation.by_symbol does
                    genes[name] = {
                        "id": gid,
                        "type": attrs.get("gene_type"),
                        "tss": s if f[6] == "+" else e - 1,
                    }
            elif gid == genes.get(name, {}).get("id"):
                exons[name].append((s, e))
    for name, g in genes.items():
        merged: list[list[int]] = []
        for s, e in sorted(exons.get(name, [])):
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        g["exons"] = merged
    return genes


def imprinted_genes() -> dict[str, set[str]]:
    """geneimprint.com's human table, by status; only symbols are kept, and only as a filter."""
    text = fetch(GENEIMPRINT_URL, CACHE / "geneimprint_human.html").read_text(errors="ignore")
    out: dict[str, set[str]] = defaultdict(set)
    for row in re.findall(r"<tr[^>]*>(.*?)</tr>", text, flags=re.S):
        cells = [
            html.unescape(re.sub(r"<[^>]+>", "", c)).strip()
            for c in re.findall(r"<td[^>]*>(.*?)</td>", row, re.S)
        ]
        if len(cells) == 5:
            out[cells[3].replace("\xa0", " ")].add(cells[0])
    return dict(out)


def lcl_expressed(genes_by_chrom) -> set[str]:
    """Target symbols with a GTEx v8 median TPM of LCL_TPM_MIN or more in EBV-transformed lymphocytes,
    matched by stable Ensembl id. Not an allelic quantity: it says whether an LCL readout can exist."""
    id_to_symbol = {g["id"]: name for gs in genes_by_chrom.values() for name, g in gs.items()}
    out: set[str] = set()
    with gzip.open(fetch(GTEX_TPM_URL, CACHE / Path(GTEX_TPM_URL).name), "rt") as fh:
        next(fh), next(fh)
        header = next(fh).rstrip("\n").split("\t")
        col = header.index(LCL_TISSUE)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sym = id_to_symbol.get(f[0].split(".")[0])
            if sym and float(f[col]) >= LCL_TPM_MIN:
                out.add(sym)
    return out


def lcl_active(elements) -> dict[str, np.ndarray]:
    """Per chromosome, which elements overlap a GM12878 H3K27ac replicated peak (ENCFF361XMX)."""
    out = {}
    for c, el in elements.items():
        path = LCL_PEAKS / f"GM12878_H3K27ac_{c}.bed.gz"
        iv: list[tuple[int, int]] = []
        if path.exists():
            with gzip.open(path, "rt") as fh:
                iv = [(int(f[0]), int(f[1])) for f in (ln.split("\t") for ln in fh if not ln.startswith("#"))]
        iv.sort()
        ps = np.array([a for a, _ in iv], dtype=np.int64)
        pe = np.maximum.accumulate(np.array([b for _, b in iv], dtype=np.int64)) if iv else ps
        # an element [s, e) overlaps a peak when the latest-ending peak starting before e ends after s
        i = np.searchsorted(ps, el["end"], side="left") - 1
        ok = i >= 0
        act = np.zeros(len(el["start"]), dtype=bool)
        act[ok] = pe[i[ok]] > el["start"][ok]
        out[c] = act
    return out


# ---------------------------------------------------------------- genotypes


def het_snvs(vcf: Path, keep=None) -> tuple[dict[str, np.ndarray], dict[str, int]]:
    """Heterozygous SNVs of a one-sample VCF as sorted 0-based positions per chromosome, with how many
    of them are written phased."""
    pos: dict[str, list[int]] = defaultdict(list)
    tally: Counter = Counter()
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 10)
            if keep is not None and not keep(f):
                continue
            if len(f[3]) != 1 or any(len(a) != 1 or a == "*" for a in f[4].split(",")):
                continue
            gt = f[9].split(":", 1)[0].rstrip("\n")
            sep = "|" if "|" in gt else "/"
            al = gt.split(sep)
            if len(al) != 2 or "." in al or al[0] == al[1]:
                continue
            pos[f[0]].append(int(f[1]) - 1)
            tally["het_snv"] += 1
            tally["phased"] += sep == "|"
    return {c: np.array(sorted(v), dtype=np.int64) for c, v in pos.items()}, dict(tally)


def read_bed(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    iv: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open() as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            iv[f[0]].append((int(f[1]), int(f[2])))
    out = {}
    for c, v in iv.items():
        v.sort()
        out[c] = (np.array([a for a, _ in v], dtype=np.int64), np.array([b for _, b in v], dtype=np.int64))
    return out


def inside(bed, chrom: str, s: np.ndarray, e: np.ndarray) -> np.ndarray:
    """Which intervals [s, e) lie wholly inside one confident interval."""
    out = np.zeros(len(s), dtype=bool)
    if chrom not in bed:
        return out
    bs, be = bed[chrom]
    i = np.searchsorted(bs, s, side="right") - 1
    ok = i >= 0
    out[ok] = be[i[ok]] >= e[ok]
    return out


def hits(h: np.ndarray | None, s: np.ndarray, e: np.ndarray) -> np.ndarray:
    """Heterozygous SNVs inside each 0-based half-open interval."""
    if h is None or not len(h):
        return np.zeros(len(s), dtype=np.int64)
    return np.searchsorted(h, e, side="left") - np.searchsorted(h, s, side="left")


def quartiles(v) -> dict[str, float] | None:
    if not len(v):
        return None
    q = np.percentile(np.asarray(v, dtype=float), [0, 25, 50, 75, 100])
    return {
        "min": float(q[0]),
        "q1": float(q[1]),
        "median": float(q[2]),
        "q3": float(q[3]),
        "max": float(q[4]),
    }


# ---------------------------------------------------------------- the per-genome count


def target_tables(elements, genes_by_chrom):
    """Per chromosome: exon arrays for every target gene, and each element's gene index (-1: none)."""
    out = {}
    for c, el in elements.items():
        genes = genes_by_chrom[c]
        names = sorted({t for t in el["target"] if t and t in genes})
        gidx = {n: i for i, n in enumerate(names)}
        xs, xe, xg = [], [], []
        for n in names:
            for s, e in genes[n]["exons"]:
                xs.append(s)
                xe.append(e)
                xg.append(gidx[n])
        out[c] = {
            "names": names,
            "tss": np.array([genes[n]["tss"] for n in names], dtype=np.int64),
            "xs": np.array(xs, dtype=np.int64),
            "xe": np.array(xe, dtype=np.int64),
            "xg": np.array(xg, dtype=np.int64),
            "elem_gene": np.array([gidx.get(t, -1) if t else -1 for t in el["target"]], dtype=np.int64),
            "not_found": sum(1 for t in el["target"] if t and t not in genes),
        }
    return out


def gene_hits(h, tb) -> np.ndarray:
    n = len(tb["names"])
    if not n:
        return np.zeros(0, dtype=bool)
    return np.bincount(tb["xg"], weights=hits(h, tb["xs"], tb["xe"]), minlength=n) > 0


def individual(label: str, het, tally, bed, elements, tables, flags, context) -> dict[str, Any]:
    """One phased genome against the element set and its element-gene pairs. `context` is (elements
    active in GM12878, targets expressed in LCLs): both genomes here were sequenced from LCLs."""
    row: Counter = Counter()
    linked: list[int] = []
    genes_both: set[str] = set()
    genes_ctx: set[str] = set()
    per_chrom = {}
    for c in CHROMS:
        el, tb = elements[c], tables[c]
        h = het.get(c)
        n = len(el["start"])
        eh = hits(h, el["start"], el["end"]) > 0
        gh = gene_hits(h, tb)
        has_t = tb["elem_gene"] >= 0
        eg = tb["elem_gene"][has_t]
        both = np.zeros(n, dtype=bool)
        both[has_t] = eh[has_t] & gh[eg]
        row["elements"] += n
        row["elements_in_confident_regions"] += int(inside(bed, c, el["start"], el["end"]).sum())
        row["elements_with_het_snv"] += int(eh.sum())
        row["pairs"] += int(has_t.sum())
        row["pairs_target_not_in_gencode"] += tb["not_found"]
        row["pairs_element_het"] += int((has_t & eh).sum())
        row["pairs_target_exon_het"] += int(gh[eg].sum())
        row["pairs_both_het"] += int(both.sum())
        if c in AUTOSOMES:
            row["pairs_both_het_autosomal"] += int(both.sum())
            active, expressed = context
            for i in np.nonzero(both)[0]:
                gname = tb["names"][tb["elem_gene"][i]]
                genes_both.add(gname)
                for flag, members in flags.items():
                    row[f"pairs_both_het_autosomal_target_{flag}"] += gname in members
                if active[c][i] and gname in expressed:
                    row["pairs_both_het_autosomal_lcl_context"] += 1
                    genes_ctx.add(gname)
            if h is not None and both.any():
                t = tb["tss"][tb["elem_gene"][both]]
                linked.extend(
                    (np.searchsorted(h, t + CIS_WINDOW) - np.searchsorted(h, t - CIS_WINDOW)).tolist()
                )
        per_chrom[c] = {"elements": n, "with_het": int(eh.sum()), "pairs_both_het": int(both.sum())}
    out: dict[str, Any] = {"label": label, **row}
    out["het_snvs"] = tally.get("het_snv", 0)
    out["het_snvs_written_phased"] = tally.get("phased", 0)
    out["fraction_elements_with_het_snv"] = round(out["elements_with_het_snv"] / out["elements"], 4)
    out["distinct_targets_both_het_autosomal"] = len(genes_both)
    out["distinct_targets_both_het_autosomal_lcl_context"] = len(genes_ctx)
    out["het_snvs_in_target_tss_window"] = quartiles(linked)
    out["per_chromosome"] = per_chrom
    return out


# ---------------------------------------------------------------- Geuvadis on the 1000 Genomes panel


def geuvadis_samples() -> set[str]:
    ids: set[str] = set()
    with fetch(GEUV_SDRF_URL, CACHE / "E-GEUV-1.sdrf.txt").open() as fh:
        for line in fh:
            ids.update(re.findall(r"\b((?:HG|NA)\d{5})\b", line))
    return ids


def stream_panel(chrom: str, wanted: set[str]) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Heterozygous SNVs per wanted sample on one chromosome of the phased panel, streamed."""
    src = HashingReader(KGP.format(chrom=chrom))
    gz = gzip.GzipFile(fileobj=io.BufferedReader(src, buffer_size=1 << 22))
    idx = np.zeros(0, dtype=np.int64)
    names: list[str] = []
    per: list[array.array] = []
    nsamp = 0
    tally: Counter = Counter()
    for raw in gz:
        if raw.startswith(b"##"):
            continue
        if raw.startswith(b"#CHROM"):
            cols = raw.rstrip(b"\n").decode().split("\t")[9:]
            names = [s for s in cols if s in wanted]
            idx = np.array([cols.index(s) for s in names], dtype=np.int64)
            per = [array.array("q") for _ in names]
            nsamp = len(cols)
            continue
        f = raw.split(b"\t", 9)
        if len(f[3]) != 1 or len(f[4]) != 1:
            continue
        tally["snv_records"] += 1
        rest = f[9]
        if f[8] == b"GT" and len(rest) == 4 * nsamp:
            a = np.frombuffer(rest, dtype=np.uint8).reshape(nsamp, 4)[idx]
            het = a[:, 0] != a[:, 2]
        else:
            gts = rest.rstrip(b"\n").split(b"\t")
            g = [gts[i].split(b":", 1)[0] for i in idx]
            het = np.array([len(x) == 3 and x[0] != x[2] for x in g])
            tally["slow_path"] += 1
        p = int(f[1]) - 1
        for i in np.nonzero(het)[0]:
            per[i].append(p)
    gz.close()
    info = src.record(chrom=chrom, samples_in_panel=nsamp, **tally)
    return {s: np.frombuffer(v, dtype=np.int64).copy() for s, v in zip(names, per, strict=True)}, info


def population(chrom: str, het: dict[str, np.ndarray], elements, tables, flags, context) -> dict[str, Any]:
    """Per element and per pair on one chromosome: how many samples could read it."""
    el, tb = elements[chrom], tables[chrom]
    n = len(el["start"])
    has_t = tb["elem_gene"] >= 0
    eg = tb["elem_gene"][has_t]
    e_count = np.zeros(n, dtype=np.int64)
    both = np.zeros(len(eg), dtype=np.int64)
    exon_only = np.zeros(len(eg), dtype=np.int64)  # target readable, element homozygous
    linked: list[int] = []
    for h in het.values():
        eh = hits(h, el["start"], el["end"]) > 0
        gh = gene_hits(h, tb)
        e_count += eh
        b = eh[has_t] & gh[eg]
        both += b
        exon_only += (~eh[has_t]) & gh[eg]
        if b.any():
            t = tb["tss"][eg[b]]
            linked.extend((np.searchsorted(h, t + CIS_WINDOW) - np.searchsorted(h, t - CIS_WINDOW)).tolist())
    names = np.array(tb["names"], dtype=object)
    imprinted = np.array([g in flags["imprinted"] for g in names[eg]], dtype=bool)
    contrast = (both >= GROUP_MIN) & (exon_only >= GROUP_MIN)
    active, expressed = context
    ctx = active[chrom][has_t] & np.array([g in expressed for g in names[eg]], dtype=bool)
    return {
        "chrom": chrom,
        "samples": len(het),
        "elements": n,
        "pairs": len(eg),
        "elements_with_het_in_at_least": {str(k): int((e_count >= k).sum()) for k in GEUV_MIN_SAMPLES},
        "median_het_samples_per_element": float(np.median(e_count)) if n else None,
        "pairs_both_het_in_at_least": {str(k): int((both >= k).sum()) for k in GEUV_MIN_SAMPLES},
        "pairs_contrastable": int(contrast.sum()),
        "pairs_contrastable_target_imprinted": int((contrast & imprinted).sum()),
        "distinct_targets_contrastable": len(set(names[eg[contrast]])),
        "pairs_lcl_context": int(ctx.sum()),
        "pairs_contrastable_lcl_context": int((contrast & ctx).sum()),
        "pairs_contrastable_lcl_context_not_imprinted": int((contrast & ctx & ~imprinted).sum()),
        "distinct_targets_contrastable_lcl_context": len(set(names[eg[contrast & ctx]])),
        "het_snvs_in_target_tss_window": quartiles(linked),
        "het_snvs_per_sample": quartiles([len(v) for v in het.values()]),
    }


# ---------------------------------------------------------------- ADASTRA and UDACHA


def udacha_url() -> str:
    q = urllib.parse.urlencode({"public_key": UDACHA_KEY})
    url = f"https://cloud-api.yandex.net/v1/disk/public/resources/download?{q}"
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
        return json.loads(r.read())["href"]


def eligible_snvs(release: str, url: str) -> tuple[dict[str, dict[str, np.ndarray]], dict[str, Any]]:
    """The eligible SNV positions of each context cell's file, pulled out of the release zip by ranges.
    Only the chromosome and the 0-based start column are read."""
    rf = RangeFile(url)
    z = zipfile.ZipFile(io.BufferedReader(rf, buffer_size=1 << 20))
    members = {i.filename: i for i in z.infolist()}
    out: dict[str, dict[str, np.ndarray]] = {}
    info: dict[str, Any] = {"archive_bytes": rf.size, "members_total": len(members), "files": {}}
    for label, name in CELL_FILES[release].items():
        if name not in members:
            info["files"][label] = {"member": name, "missing": True}
            continue
        data = z.read(name)
        reader = csv.reader(io.TextIOWrapper(io.BytesIO(data), encoding="utf-8"), delimiter="\t")
        header = next(reader)
        ci = header.index("chr" if "chr" in header else "#chr")
        si = header.index("start")
        pos: dict[str, list[int]] = defaultdict(list)
        rows = 0
        for r in reader:
            pos[r[ci]].append(int(r[si]))
            rows += 1
        out[label] = {c: np.array(sorted(v), dtype=np.int64) for c, v in pos.items()}
        info["files"][label] = {
            "member": name,
            "bytes": members[name].file_size,
            "crc32": f"{members[name].CRC:08x}",
            "sha256": hashlib.sha256(data).hexdigest(),
            "rows": rows,
            "header": header,
        }
    info["bytes_fetched"] = rf.fetched
    return out, info


def element_reach(snvs: dict[str, dict[str, np.ndarray]], elements) -> dict[str, Any]:
    out = {}
    n_el = sum(len(el["start"]) for el in elements.values())
    for label, per in snvs.items():
        n_hit = n_in = n_snv = 0
        for c in CHROMS:
            el = elements[c]
            h = per.get(c, np.zeros(0, dtype=np.int64))
            k = hits(h, el["start"], el["end"])
            n_hit += int((k > 0).sum())
            n_in += int(k.sum())
            n_snv += len(h)
        out[label] = {
            "eligible_snvs": n_snv,
            "eligible_snvs_inside_elements": n_in,
            "elements_with_eligible_snv": n_hit,
            "fraction_of_elements": round(n_hit / n_el, 5),
        }
    return out


# ---------------------------------------------------------------- EN-TEx


CHROMATIN = ("ATAC", "DNase", "ChIP")  # every EN-TEx chromatin assay name carries one of these


def _entex_rows(name: str) -> tuple[HashingReader, io.TextIOWrapper]:
    src = HashingReader(ENTEX + name)
    text = io.TextIOWrapper(io.BufferedReader(src, buffer_size=1 << 22), encoding="utf-8")
    header = text.readline().rstrip("\n").split("\t")
    if header != ENTEX_HEADER:
        raise SystemExit(f"EN-TEx {name} header changed: {header}")
    return src, text


def entex(elements, genes_by_chrom) -> dict[str, Any]:
    """Accessible regions of EN-TEx's AS catalogue, by donor and tissue, matched to the archive's
    elements by coordinate overlap (EN-TEx names registry-V2 accessions, the archive V3) and to target
    genes by stable Ensembl id. Only the identifier columns are read."""
    id_to_symbol = {g["id"]: name for gs in genes_by_chrom.values() for name, g in gs.items()}
    elem_acc: dict[tuple[str, str], set[tuple[str, int]]] = defaultdict(set)
    assays: Counter = Counter()
    donor_rows: Counter = Counter()
    overlap: dict[tuple[str, int, int], list[int]] = {}
    src, text = _entex_rows("cCREs_default_AS.tsv")
    rows = 0
    for line in text:
        f = line.split("\t", 10)
        rows += 1
        assays[f[9]] += 1
        donor_rows[f[7]] += 1
        c = f[0]
        if c not in elements or not any(k in f[9] for k in CHROMATIN):
            continue
        key = (c, int(f[1]), int(f[2]))
        if key not in overlap:  # archive elements overlapping [s, e): start before e, end after s
            st, en = elements[c]["start"], elements[c]["end"]
            hi = int(np.searchsorted(st, key[2], side="left"))
            overlap[key] = [i for i in range(max(0, hi - 50), hi) if en[i] > key[1]]
        for i in overlap[key]:
            elem_acc[(f[7], f[8])].add((c, i))
    text.close()
    info: dict[str, Any] = {
        "ccre_table": src.record(rows=rows),
        "ccre_rows_by_assay": dict(assays.most_common()),
        "ccre_rows_by_donor": dict(donor_rows.most_common()),
        "distinct_chromatin_regions": len(overlap),
        "distinct_chromatin_regions_overlapping_an_element": sum(1 for v in overlap.values() if v),
    }
    gene_acc: dict[tuple[str, str], set[str]] = defaultdict(set)
    gassays: Counter = Counter()
    src, text = _entex_rows("genes_default_AS.tsv")
    grows = unmapped = 0
    for line in text:
        f = line.split("\t", 10)
        grows += 1
        gassays[f[9]] += 1
        if "RNA" not in f[9]:
            continue
        sym = id_to_symbol.get(f[3].split(".")[0])
        if sym is None:
            unmapped += 1
        else:
            gene_acc[(f[7], f[8])].add(sym)
    text.close()
    info["gene_table"] = src.record(rows=grows)
    info["gene_rows_by_assay"] = dict(gassays.most_common())
    info["gene_rna_rows_not_in_gencode_v50"] = unmapped
    any_elem: set[tuple[str, int]] = set()
    all_pairs: set[tuple[str, int]] = set()
    per_donor = {}
    for d in sorted({d for d, _ in elem_acc} | {d for d, _ in gene_acc}):
        e_d: set[tuple[str, int]] = set()
        pairs_d: set[tuple[str, int]] = set()
        tissues = sorted({t for dd, t in elem_acc if dd == d} | {t for dd, t in gene_acc if dd == d})
        for t in tissues:
            ea, ga = elem_acc.get((d, t), set()), gene_acc.get((d, t), set())
            e_d |= ea
            pairs_d |= {(c, i) for c, i in ea if elements[c]["target"][i] in ga}
        any_elem |= e_d
        all_pairs |= pairs_d
        per_donor[d] = {
            "tissues": len(tissues),
            "tissues_with_chromatin_and_rna": sum(
                1 for t in tissues if (d, t) in elem_acc and (d, t) in gene_acc
            ),
            "elements_accessible_chromatin": len(e_d),
            "pairs_element_and_target_accessible_same_tissue": len(pairs_d),
        }
    n_el = sum(len(el["start"]) for el in elements.values())
    n_pairs = sum(sum(1 for t in el["target"] if t) for el in elements.values())
    info["elements_accessible_any_donor"] = len(any_elem)
    info["fraction_of_elements"] = round(len(any_elem) / n_el, 5)
    info["pairs_accessible_same_donor_and_tissue"] = len(all_pairs)
    info["fraction_of_pairs"] = round(len(all_pairs) / n_pairs, 5)
    info["pairs_strong_accessible"] = sum(1 for c, i in all_pairs if elements[c]["strength"][i] == "strong")
    info["per_donor"] = per_donor
    return info


# ---------------------------------------------------------------- design numbers, no data


def reads_needed(ratio: float, alpha: float = 0.05, power: float = 0.8) -> int:
    """Allelic reads one binomial test needs to tell `ratio` from 0.5 (normal approximation, two-sided);
    overdispersion only raises it."""
    z = NormalDist().inv_cdf
    n = ((z(1 - alpha / 2) * 0.5 + z(power) * math.sqrt(ratio * (1 - ratio))) / (ratio - 0.5)) ** 2
    return math.ceil(n)


# ---------------------------------------------------------------- the candidates, as read 2026-09-28


CANDIDATES = [
    {
        "name": "GTEx v8 allelic expression (Castel et al. 2020, Genome Biology)",
        "paired": "RNA-seq haplotype counts per gene per sample (phASER, with and without WASP) with the "
        "donors' WGS genotypes",
        "samples": "15,253 samples, 54 tissues, 838 donors",
        "access": "open: the phASER haplotype expression matrices on the GTEx portal (adult-gtex bucket, "
        "haplotype-expression/v8/); controlled (dbGaP phs000424.v8): SNP-level ASE, the WGS genotypes and "
        "the read-back phased genotypes. No haplotype expression in v10 or v11 (bulk-qtl and bulk-gex only)",
        "licence": "GTEx portal open-access data; article CC BY 4.0; genotypes under dbGaP terms",
        "size": "phASER_WASP_GTEx_v8_matrix.gw_phased.txt.gz 461.5 MB (the three other matrices 513.6 to "
        "555.0 MB)",
        "format": "one string per gene per sample, HAP_A_COUNT|HAP_B_COUNT; no variant ids, no genotypes",
        "allele_counts": "provided (gene-level haplotype counts)",
        "assembly": "GRCh38",
        "feasible_here": "the open matrix pairs with nothing open: which haplotype carries an element's "
        "allele is in the controlled genotypes, so no element-level prediction can be checked without dbGaP "
        "access (not requested)",
    },
    {
        "name": "Geuvadis (Lappalainen et al. 2013), ArrayExpress E-GEUV-1",
        "paired": "LCL RNA-seq ASE per individual per site, with 1000 Genomes phased genotypes",
        "samples": "462 unique RNA samples (EUR 373, YRI 89); one cell type (LCL)",
        "access": "open (BioStudies E-GEUV-1, no registration); the 1000 Genomes 30x phased panel is open",
        "licence": "EMBL-EBI terms of use; IGSR: 1000 Genomes data available without embargo",
        "size": "GD462.ASE.COV8.ANNOT_PTV.txt.gz 918 MB; phased genotypes 621 MB (chr22) to 3.5 GB (chr1); "
        "one BAM per sample",
        "format": "ASE: INDIVIDUAL, RSID, CHR, POS, ALLELES, REF_COUNT, NONREF_COUNT, TOTAL_COUNT, REF_RATIO "
        "(per-sample null), binomial PVALUE, GENOTYPE; the name says COV8, a coverage floor the README "
        "does not define",
        "allele_counts": "provided",
        "assembly": "GRCh37 (1000 Genomes phase 1); a lift or a GRCh38 re-count is needed",
        "feasible_here": "yes for the table (under 1 GB); the genotype side can come from the GRCh38 30x "
        "phased panel instead, which this probe streams for chr21",
    },
    {
        "name": "EN-TEx (Rozowsky et al. 2023, Cell)",
        "paired": "haplotype-resolved read counts for ATAC, DNase, histone and CTCF/POLR2A ChIP and RNA-seq "
        "on each donor's own phased diploid genome",
        "samples": "4 donors, about 30 tissues, 1,635 datasets; the AS pipeline ran on about 1,000 samples",
        "access": "open: 'fully open-consented and accessible without registration'; the four donors' "
        "personalized genome assemblies on the ENCODE portal (ENCSR792RVV, ENCSR866HYZ, ENCSR780BKN, "
        "ENCSR083ETH; 1.7 to 3.6 GB each)",
        "licence": "ENCODE data use policy ('freely download, analyze and publish ... without "
        "restrictions'); article CC BY 4.0",
        "size": "cCREs_default_AS.tsv 711 MB; genes_default_AS.tsv 94 MB; hetSNVs_default_AS.tsv 2.5 GB",
        "format": "chr, start, end, region_id (cCRE V2 accession or Ensembl gene id), hap1_count, "
        "hap2_count, experiment, donor, tissue, assay, hap1_allele_ratio, p_betabinom, "
        "imbalance_significance",
        "allele_counts": "provided (AlleleSeq2: reads mapped to both personal haplotypes, "
        "beta-binomial test)",
        "assembly": "GRCh38 coordinates in the tables",
        "feasible_here": "yes: both tables are under 1 GB and keyed to cCREs and genes, the element and "
        "readout sides of the archive's pairs",
    },
    {
        "name": "ADASTRA release Mabel v6.1 (Nov 2024; allele-specific TF binding)",
        "paired": "ChIP-seq allelic read counts at heterozygous SNVs called from the reads themselves",
        "samples": "1,073 TFs and 649 cell types from 15,970 GTRD alignments (Bill Cipher v5 base)",
        "access": "open (Zenodo 14174114); the ADASTRA site's downloads page returned HTTP 502 on 2026-09-28",
        "licence": "CC BY 4.0 (Zenodo record)",
        "size": "ADASTRA.v.6.1.Mabel.zip 942.7 MB; the four context-cell files 26 MB compressed",
        "format": "per cell type and per TF: every eligible SNV passing the coverage thresholds, with "
        "aggregated effect sizes and FDR per allele (hg38, 0-based start)",
        "allele_counts": "aggregated statistics per SNV per cell type, not per-sample counts",
        "assembly": "GRCh38",
        "feasible_here": "yes, by HTTP ranges on the zip; binding only, so it reads the element side and "
        "never the gene",
    },
    {
        "name": "UDACHA release IceKing v1.0.3 (June 2023; allele-specific accessibility)",
        "paired": "DNase-, ATAC- and FAIRE-seq allelic read counts at SNVs called from the reads",
        "samples": "5,858 chromatin accessibility datasets from GTRD",
        "access": "open (Yandex Disk link on udacha.autosome.org)",
        "licence": "no data licence found on the site; the article (Nat Commun 2025) is CC BY-NC-ND 4.0",
        "size": "UDACHA_IceKing_May2023_release.zip 826.5 MB; the context-cell files 12 MB compressed",
        "format": "per cell type per assay: #chr, start, end, mean_bad, id, max_cover, ref, alt, n_reps, "
        "then effect sizes, p-values and FDRs; no readme in the zip, so that a file lists every "
        "coverage-passing SNV and not only the significant ones is not stated by this release",
        "allele_counts": "aggregated statistics per SNV per cell type",
        "assembly": "GRCh38",
        "feasible_here": "yes, by HTTP ranges; accessibility only, element side",
    },
    {
        "name": "AlleleDB (Chen et al. 2016)",
        "paired": "ASB and ASE on 1000 Genomes individuals' personal genomes",
        "samples": "382 individuals",
        "access": "open (archive.gersteinlab.org/proj/alleledb/download/), last updated 2016-08",
        "licence": "not stated",
        "size": "accE v2.1 96 MB, accB 6.1 MB, ASE 7.6 MB, ASB 213 kB",
        "format": "accessible SNVs and AS calls",
        "allele_counts": "provided for accessible SNVs",
        "assembly": "GRCh37",
        "feasible_here": "reachable, but for this purpose EN-TEx (same lab, personal genomes, GRCh38, "
        "tissues) and ADASTRA/UDACHA (binding, accessibility) cover it",
    },
    {
        "name": "ENCODE4 Hi-C genophasing (43 biosamples)",
        "paired": "phased variant calls from Hi-C, including HepG2 and IMR-90; its 'allele-specific "
        "variants' output is the Hi-C diploid layout, not a molecular readout",
        "samples": "43 annotations (HepG2 ENCSR777ARQ, IMR-90 ENCSR425DIX; none for K562 or GM12878)",
        "access": "open (ENCODE portal)",
        "licence": "ENCODE data use policy (no restrictions)",
        "size": "phased VCF 122 to 198 MB each",
        "format": "VCF",
        "allele_counts": "none; would be computed from ENCODE RNA, DNase or ATAC BAMs",
        "assembly": "GRCh38",
        "feasible_here": "genotype side only; allele counting from multi-GB BAMs is possible by ranges "
        "(genomeos/genome/bam_range.py) but is a build, not a probe",
    },
    {
        "name": "GIAB HG002 RNA-seq (NIST RNA-seq pilot)",
        "paired": "Illumina mRNA and lncRNA, PacBio and ONT reads on three HG002 cell stocks (GM24385, "
        "GM26105, GM27730), with the Q100 assembly phased by parent",
        "samples": "one individual",
        "access": "open (GIAB FTP)",
        "licence": "NIST data use policy (17 USC 105)",
        "size": "mRNA BAMs 7.3 to 9.1 GB each (.bai 4 to 5 MB); FASTQ 6 to 7 GB per mate",
        "format": "reads only",
        "allele_counts": "none; must be computed",
        "assembly": "GRCh38",
        "feasible_here": "by range reads of the indexed BAM at the heterozygous sites this probe counts; one "
        "individual, so it can check a direction but cannot separate an element from its linked variants",
    },
]


# ---------------------------------------------------------------- main


def main() -> None:
    t0 = time.time()
    CACHE.mkdir(parents=True, exist_ok=True)
    elements = load_elements()
    print(
        f"elements {sum(len(e['start']) for e in elements.values()):,} ({time.time() - t0:.0f}s)", flush=True
    )
    genes_by_chrom = {c: load_genes(c) for c in CHROMS}
    tables = target_tables(elements, genes_by_chrom)
    imp = imprinted_genes()
    flags = {
        "imprinted": imp.get("Imprinted", set()),
        "imprinted_predicted": imp.get("Predicted", set()),
        "immunoglobulin_or_HLA": {
            g for gs in genes_by_chrom.values() for g in gs if g.startswith(("IGH", "IGK", "IGL", "HLA-"))
        },
    }
    context = (lcl_active(elements), lcl_expressed(genes_by_chrom))
    print(f"genes, flags and LCL context ({time.time() - t0:.0f}s)", flush=True)

    hg001_vcf = fetch(HG001_VCF_URL, CACHE / Path(HG001_VCF_URL).name)
    hg001_bed = fetch(HG001_BED_URL, CACHE / Path(HG001_BED_URL).name)
    het, tally = het_snvs(hg001_vcf)
    label = "NA12878 (GM12878), GIAB HG001 v4.2.1"
    na12878 = individual(label, het, tally, read_bed(hg001_bed), elements, tables, flags, context)
    print("NA12878", {k: v for k, v in na12878.items() if k != "per_chromosome"}, flush=True)
    het, tally = het_snvs(HG002_VCF, keep=lambda f: f[6] in (".", "PASS"))
    hg002 = individual(
        "HG002, T2T Q100 v1.1 dipcall", het, tally, read_bed(HG002_BED), elements, tables, flags, context
    )
    print("HG002", {k: v for k, v in hg002.items() if k != "per_chromosome"}, flush=True)
    del het

    wanted = geuvadis_samples()
    phet, pinfo = stream_panel(KGP_CHROM, wanted)
    geuv = population(KGP_CHROM, phet, elements, tables, flags, context)
    geuv["sdrf_sample_ids"] = len(wanted)
    geuv["panel"] = pinfo
    print("Geuvadis", geuv, flush=True)
    del phet

    adastra_snvs, adastra_info = eligible_snvs("ADASTRA", ADASTRA_URL)
    adastra = {"release": adastra_info, "reach": element_reach(adastra_snvs, elements)}
    print("ADASTRA", adastra["reach"], flush=True)
    udacha_snvs, udacha_info = eligible_snvs("UDACHA", udacha_url())
    udacha = {"release": udacha_info, "reach": element_reach(udacha_snvs, elements)}
    print("UDACHA", udacha["reach"], flush=True)
    del adastra_snvs, udacha_snvs

    en = entex(elements, genes_by_chrom)
    print("EN-TEx", {k: v for k, v in en.items() if not k.endswith("by_assay")}, flush=True)

    payload = {
        "question": "Which open resource pairs genotype with an allelic molecular readout, and how much of "
        "the attribution layer's element set could it reach? Feasibility only: counts of overlap, no "
        "allelic outcome read, nothing registered that measures.",
        "universe": {
            "elements": sum(len(e["start"]) for e in elements.values()),
            "pairs_with_predicted_coding_target": sum(
                sum(1 for t in e["target"] if t) for e in elements.values()
            ),
            "pairs_strong": sum(sum(1 for s in e["strength"] if s == "strong") for e in elements.values()),
            "pairs_target_not_in_gencode_v50": sum(t["not_found"] for t in tables.values()),
            "elements_active_in_gm12878_h3k27ac": int(sum(a.sum() for a in context[0].values())),
            "pairs_lcl_context": sum(
                1
                for c, el in elements.items()
                for i, t in enumerate(el["target"])
                if t and context[0][c][i] and t in context[1]
            ),
            "targets_expressed_in_lcl_gtex_v8": len(context[1]),
            "where": str(ARCHIVE),
        },
        "candidates": CANDIDATES,
        "overlap": {
            "NA12878": na12878,
            "HG002": hg002,
            "Geuvadis_1kGP_chr21": geuv,
            "ADASTRA": adastra,
            "UDACHA": udacha,
            "EN-TEx": en,
        },
        "flags": {k: len(v) for k, v in flags.items()},
        "bias_controls_design": {
            "reads_for_80pct_power_single_binomial_test_alpha_0.05": {
                str(r): reads_needed(r) for r in READ_DEPTH_EFFECTS
            },
            "cis_window_bp": CIS_WINDOW,
            "contrast_min_samples_each_side": GROUP_MIN,
        },
        "label_blind": {"columns_read": COLUMNS_READ, "columns_never_read": COLUMNS_NEVER_READ},
        "alphagenome_requests": 0,
        "seconds": round(time.time() - t0),
    }
    save_result(RESULT, payload, manifest=manifest(pinfo, adastra_info, udacha_info, en))
    print(f"saved in {time.time() - t0:.0f}s", flush=True)


def manifest(pinfo, adastra_info, udacha_info, en) -> dict[str, Any]:
    """The provenance contract (review item R9): every file read, streamed or pulled by range."""
    inputs = [mf.input_entry(ARCHIVE / f"{c}.json", partition=None) for c in CHROMS]
    inputs += [
        mf.input_entry(Path("data/reference") / f"gencode_v50_{c}.gff3.gz", partition=None) for c in CHROMS
    ]
    inputs += [
        mf.input_entry(CACHE / Path(HG001_VCF_URL).name, partition=None, url=HG001_VCF_URL),
        mf.input_entry(CACHE / Path(HG001_BED_URL).name, partition=None, url=HG001_BED_URL),
        mf.input_entry(HG002_VCF, partition=None),
        mf.input_entry(HG002_BED, partition=None),
        mf.input_entry(CACHE / "E-GEUV-1.sdrf.txt", partition=None, url=GEUV_SDRF_URL),
        mf.input_entry(CACHE / "geneimprint_human.html", partition=None, url=GENEIMPRINT_URL),
        mf.input_entry(CACHE / Path(GTEX_TPM_URL).name, partition=None, url=GTEX_TPM_URL),
        *[mf.input_entry(p, partition=None) for p in sorted(LCL_PEAKS.glob("GM12878_H3K27ac_chr*.bed.gz"))],
        {**pinfo, "path": pinfo["url"], "partition": None, "streamed": True},
        {**en["ccre_table"], "path": en["ccre_table"]["url"], "partition": None, "streamed": True},
        {**en["gene_table"], "path": en["gene_table"]["url"], "partition": None, "streamed": True},
    ]
    for rel, url, info in (("ADASTRA", ADASTRA_URL, adastra_info), ("UDACHA", UDACHA_KEY, udacha_info)):
        for cell, f in info["files"].items():
            if "sha256" in f:
                inputs.append(
                    {
                        "path": f"{url}#{f['member']}",
                        "sha256": f["sha256"],
                        "bytes": f["bytes"],
                        "partition": None,
                        "release": rel,
                        "cell": cell,
                        "read_by_http_range": True,
                    }
                )
    sources = [
        {
            "accession": "ENCODE SCREEN cCREs scored by AlphaGenome deletion (all-element archive)",
            "version": "AlphaGenome as served during the 2026-09 all-element sweep (unpinned); pinned here "
            "by sha256",
            "path": str(ARCHIVE),
        },
        {"accession": "GENCODE human gene annotation", "version": "release 50 (GRCh38.p14)"},
        {
            "accession": "GIAB HG001 (NA12878) benchmark small variants",
            "version": "NISTv4.2.1 GRCh38",
            "url": GIAB,
        },
        {
            "accession": "GIAB HG002 T2T Q100 dipcall calls and v5.0q small-variant benchmark bed",
            "version": "v5.0q",
        },
        {
            "accession": "1000 Genomes 30x phased panel (NYGC, Byrska-Bishop et al. 2022)",
            "version": "20220422_3202_phased_SNV_INDEL_SV",
            "url": pinfo["url"],
        },
        {
            "accession": "Geuvadis E-GEUV-1 (sample list only)",
            "version": "BioStudies release, modified 2022-11-02",
            "url": GEUV_SDRF_URL,
        },
        {"accession": "ADASTRA (Zenodo 14174114)", "version": "Mabel v6.1, Nov 2024", "url": ADASTRA_URL},
        {"accession": "UDACHA", "version": "IceKing v1.0.3, 2023-06-13", "url": UDACHA_KEY},
        {
            "accession": "EN-TEx AS catalogue (Rozowsky et al. 2023)",
            "version": "cCREs_default_AS.tsv and genes_default_AS.tsv dated 2021-12-10",
            "url": ENTEX,
        },
        {
            "accession": "geneimprint.com human imprinted genes (symbols only, as a filter)",
            "version": "as served 2026-09-28",
            "url": GENEIMPRINT_URL,
        },
    ]
    sources += [
        {"accession": "GTEx v8 gene median TPM (LCL column only)", "version": "2017-06-05 v8 RNASeQC 1.1.9"},
        {
            "accession": "ENCODE ENCFF361XMX (ENCSR000AKC) GM12878 H3K27ac replicated peaks",
            "version": "as cached",
        },
    ]
    m = {
        "sources": sources,
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "cis_window_bp": CIS_WINDOW,
            "population_chromosome": KGP_CHROM,
            "geuvadis_thresholds": list(GEUV_MIN_SAMPLES),
            "contrast_min_samples_each_side": GROUP_MIN,
            "het": "SNV whose GT holds two different alleles; indels and symbolic alleles excluded",
            "element_overlap": "at least one base (EN-TEx V2 regions to archive V3 elements)",
            "target_exons": "union of every GENCODE v50 exon of the target symbol's gene",
            "chromatin_assays": list(CHROMATIN),
            "lcl_context": "element overlaps a GM12878 H3K27ac peak and target GTEx v8 median TPM >= "
            f"{LCL_TPM_MIN} in '{LCL_TISSUE}'",
        },
        "exclusions": [
            "indels and symbolic alleles: allele counting reads SNVs",
            "chrX and chrY from the autosomal pair counts of the two genomes (X inactivation in NA12878, "
            "hemizygous in HG002)",
            "HG002 dipcall calls with a FILTER other than PASS or '.'",
            "the population count covers chr21 only; the panel's other chromosomes were not streamed",
            "every allelic count, ratio, effect-size and significance column (label-blind)",
        ],
        "partitions": "n/a: a feasibility census; no evaluation partition is scored",
    }
    return mf.with_model_dependencies(m, "alphagenome")


if __name__ == "__main__":
    main()
