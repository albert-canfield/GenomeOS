# SPDX-License-Identifier: AGPL-3.0-or-later
"""The second mammal: a mouse chromosome through the same code, and its nodes held against ours.

The node model says regulation is addressed by CTCF-bounded domains. If that is a property of
mammalian genomes and not of one annotation, a mouse chromosome fetched and analysed exactly like a
human one (sequence from UCSC, gene models from GENCODE, candidate elements from ENCODE SCREEN)
should give nodes whose gene content maps onto single human nodes. This module fetches mouse mm10
(the assembly ENCODE's mouse registry is on), infers its domains with the human code path, and
compares node membership across species by gene symbol.

Symbol identity is the orthology here (App ↔ APP): cheap, mostly right for one-to-one orthologues,
wrong for renamed and duplicated genes, so the conservation it reports is a lower bound. Files carry
the `mm10_` prefix and never mix with the human ones; the mouse cCRE subset and the summary are
committed, the sequence and gene models stay local.
"""

from __future__ import annotations

import gzip
import io
import json
import urllib.request
from pathlib import Path
from typing import Any

MM10_FASTA = "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/{chrom}.fa.gz"
GENCODE_MOUSE = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gff3.gz"
)
MOUSE_CCRES = "https://downloads.wenglab.org/V3/mm10-cCREs.bed"
REFERENCE = Path("data/reference")
RESULTS = Path("data/results")
EVIDENCE = {
    "sequence": "curated: UCSC mm10 (GRCm38)",
    "genes": "curated: GENCODE vM25 (GRCm38)",
    "elements": "curated: ENCODE SCREEN registry of mouse cCREs v3 (mm10)",
    "domains": "inferred: CTCF-only elements as boundaries, the same code as for human",
    "orthology": "curated: MGI mouse–human homology classes (symbol identity kept as the fallback)",
}


def _download(url: str, dest: Path, progress=None) -> Path:
    if dest.exists():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "GenomeOS/0.1 (stream)"})
    with urllib.request.urlopen(req, timeout=900) as r, open(dest, "wb") as fh:  # noqa: S310
        while chunk := r.read(1 << 20):
            fh.write(chunk)
    if progress:
        progress(f"{dest.name}: {dest.stat().st_size / 1e6:.1f} MB")
    return dest


def fetch_sequence(chrom: str, progress=None) -> Path:
    return _download(MM10_FASTA.format(chrom=chrom), REFERENCE / f"mm10_{chrom}.fa.gz", progress)


def fetch_gencode(chrom: str, progress=None) -> Path:
    """Stream the mouse GENCODE GFF3 and keep this chromosome's rows."""
    dest = REFERENCE / f"gencode_vM25_{chrom}.gff3.gz"
    if dest.exists():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(GENCODE_MOUSE, headers={"User-Agent": "GenomeOS/0.1 (stream)"})
    kept = 0
    with (
        urllib.request.urlopen(req, timeout=900) as resp,  # noqa: S310
        gzip.open(io.BufferedReader(resp, 1 << 20), "rt") as fh,
        gzip.open(dest, "wt") as out,
    ):
        for line in fh:
            if line.startswith("#"):
                if line.startswith("##gff-version"):
                    out.write(line)
                continue
            if line.split("\t", 1)[0] == chrom:
                out.write(line)
                kept += 1
    if progress:
        progress(f"GENCODE vM25 rows kept for mouse {chrom}: {kept:,}")
    return dest


def ccres_path(chrom: str) -> Path:
    return RESULTS / f"ccres_mm10_{chrom}.bed.gz"


def fetch_ccres(chrom: str, progress=None) -> Path:
    """Stream the mouse registry and keep this chromosome's elements (same columns as the human file)."""
    from genomeos.genome.regulatory import stream_ccres

    dest = ccres_path(chrom)
    if dest.exists():
        return dest
    elements = stream_ccres({chrom}, url=MOUSE_CCRES, progress=None)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(dest, "wt") as fh:
        fh.write(f"# ENCODE mouse cCREs v3 (mm10), {chrom} subset distilled by GenomeOS from {MOUSE_CCRES}\n")
        for c in elements:
            fh.write(f"{c.chrom}\t{c.start}\t{c.end}\t{c.id}\t{c.cls}\t{int(c.ctcf_bound)}\n")
    if progress:
        progress(f"mouse cCREs kept for {chrom}: {len(elements):,}")
    return dest


def load_ccres(chrom: str):
    from genomeos.genome.regulatory import CCRE

    p = ccres_path(chrom)
    out = []
    if p.exists():
        with gzip.open(p, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                out.append(CCRE(f[0], int(f[1]), int(f[2]), f[3], f[4], f[5] == "1"))
    return out


MGI_HOMOLOGY = "https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"


def orthology_path(results_dir: Path = RESULTS) -> Path:
    return results_dir / "mgi_mouse_human_orthology.tsv.gz"


def fetch_orthology(results_dir: Path = RESULTS, progress=None) -> Path:
    """Stream MGI's curated mouse–human homology report once (15 MB) and keep only the symbol pairs
    (one line per pair, a few hundred KB): the orthology the node comparison rests on."""
    dest = orthology_path(results_dir)
    if dest.exists():
        return dest
    req = urllib.request.Request(MGI_HOMOLOGY, headers={"User-Agent": "GenomeOS/0.1 (stream)"})
    groups: dict[str, dict[str, list[str]]] = {}
    with urllib.request.urlopen(req, timeout=600) as resp:  # noqa: S310
        for line in io.TextIOWrapper(io.BufferedReader(resp, 1 << 20), encoding="utf-8", errors="replace"):
            f = line.rstrip("\n").split("\t")
            if len(f) < 4 or f[0] == "DB Class Key":
                continue
            org = "mouse" if f[1].startswith("mouse") else "human" if f[1] == "human" else None
            if org:
                groups.setdefault(f[0], {"mouse": [], "human": []})[org].append(f[3])
    pairs = sorted({(mo, hu) for g in groups.values() for mo in g["mouse"] for hu in g["human"]})
    dest.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(dest, "wt") as fh:
        fh.write(f"# MGI mouse–human homology classes distilled to symbol pairs from {MGI_HOMOLOGY}\n")
        for mo, hu in pairs:
            fh.write(f"{mo}\t{hu}\n")
    if progress:
        progress(f"MGI orthology: {len(pairs):,} mouse–human symbol pairs kept from {len(groups):,} classes")
    return dest


def load_orthology(results_dir: Path = RESULTS) -> dict[str, list[str]]:
    """Mouse symbol → human symbols (curated, MGI); empty if the report has not been fetched."""
    p = orthology_path(results_dir)
    out: dict[str, list[str]] = {}
    if p.exists():
        with gzip.open(p, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                mo, _, hu = line.rstrip("\n").partition("\t")
                out.setdefault(mo, []).append(hu)
    return out


def human_nodes_by_symbol(results_dir: Path = RESULTS) -> tuple[dict[str, str], int]:
    """Gene symbol → human node id from every committed domains_chr*.json; also the node count."""
    index: dict[str, str] = {}
    n = 0
    for p in sorted(results_dir.glob("domains_chr*.json")):
        try:
            d = json.loads(p.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        for dom in d.get("domains", []):
            n += 1
            for g in dom.get("genes", []):
                index.setdefault(g.upper(), dom["id"])
    return index, n


def compare_nodes(
    mouse_domains: list[dict[str, Any]],
    human_index: dict[str, str],
    orthology: dict[str, list[str]] | None = None,
) -> dict[str, Any]:
    """Does a mouse node's gene content land in one human node? Per mouse node with two or more
    coding genes that have a human match: conserved (one human node), split (several), and the
    nodes with fewer than two matches (unmapped). With `orthology` (MGI) a mouse gene maps through
    its curated human orthologues; without it, through symbol identity."""

    def human_node(symbol: str) -> str | None:
        candidates = (orthology or {}).get(symbol) or [symbol.upper()]
        for h in candidates:
            node = human_index.get(h.upper())
            if node:
                return node
        return None

    rows = []
    counts = {"conserved": 0, "split_adjacent": 0, "split_scattered": 0, "unmapped": 0}
    for dom in mouse_domains:
        coding = [g for g in dom.get("coding_symbols", []) if g]
        hits = {g: human_node(g) for g in coding}
        mapped = {g: h for g, h in hits.items() if h}
        targets = sorted(set(mapped.values()), key=_node_key)
        if len(mapped) < 2:
            verdict = "unmapped"
        elif len(targets) == 1:
            verdict = "conserved"
        elif _adjacent(targets):
            verdict = "split_adjacent"  # one human neighbourhood, boundaries drawn differently
        else:
            verdict = "split_scattered"  # genes that are neighbours in mouse live apart in human
        counts[verdict] += 1
        rows.append(
            {
                "id": dom["id"],
                "start": dom["start"],
                "end": dom["end"],
                "coding_genes": len(coding),
                "mapped": len(mapped),
                "human_nodes": targets,
                "verdict": verdict,
                "genes": coding[:12],
            }
        )
    tested = counts["conserved"] + counts["split_adjacent"] + counts["split_scattered"]
    same_place = counts["conserved"] + counts["split_adjacent"]
    return {
        "orthology": "curated: MGI mouse–human homology" if orthology else "inferred: symbol identity",
        "nodes": len(rows),
        "tested": tested,
        **counts,
        "fraction_conserved": round(counts["conserved"] / tested, 3) if tested else None,
        "fraction_same_neighbourhood": round(same_place / tested, 3) if tested else None,
        "rows": rows,
    }


def _node_key(node_id: str) -> tuple[str, int]:
    chrom, _, d = node_id.partition(":")
    return chrom, int(d[1:]) if d[1:].isdigit() else -1


def _adjacent(node_ids: list[str]) -> bool:
    """All on one human chromosome with consecutive node numbers (allowing one gap for an empty node)."""
    keys = sorted(_node_key(n) for n in node_ids)
    if len({c for c, _ in keys}) != 1:
        return False
    nums = [i for _, i in keys]
    return all(b - a <= 2 for a, b in zip(nums, nums[1:], strict=False))


def analyse(chrom: str, progress=None) -> dict[str, Any]:
    """Fetch what is missing, infer the mouse chromosome's nodes with the human code, compare."""
    from genomeos.genome import Annotation, Genome
    from genomeos.genome.domains import infer_domains

    fa = fetch_sequence(chrom, progress)
    gff = fetch_gencode(chrom, progress)
    fetch_ccres(chrom, progress)
    ccres = load_ccres(chrom)
    ann = Annotation.from_gff3(gff, {chrom})
    genome = Genome.from_fasta(fa)
    length = genome.chromosomes[chrom].length
    domains = infer_domains(chrom, length, ccres, ann)
    by_symbol = {g.symbol: g for g in ann.genes.values() if g.locus.chrom == chrom}
    rows = []
    for d in domains:
        dd = d.to_dict()
        dd["coding_symbols"] = [
            s for s in d.genes if s in by_symbol and by_symbol[s].type == "protein_coding"
        ]
        rows.append(dd)
    coding = [g for g in ann.genes.values() if g.locus.chrom == chrom and g.type == "protein_coding"]
    human_index, human_nodes = human_nodes_by_symbol()
    fetch_orthology(progress=progress)
    orthology = load_orthology()
    cmp = compare_nodes(rows, human_index, orthology)
    cmp_symbol = compare_nodes(rows, human_index)
    cls: dict[str, int] = {}
    for c in ccres:
        cls[c.cls] = cls.get(c.cls, 0) + 1
    out = {
        "organism": "mouse",
        "assembly": "mm10",
        "chrom": chrom,
        "length": length,
        "genes": sum(1 for g in ann.genes.values() if g.locus.chrom == chrom),
        "coding_genes": len(coding),
        "elements": len(ccres),
        "element_classes": dict(sorted(cls.items(), key=lambda kv: -kv[1])),
        "domains": len(domains),
        "domains_with_coding_genes": sum(1 for r in rows if r["coding_symbols"]),
        "domain_length_median": sorted(d.length for d in domains)[len(domains) // 2] if domains else None,
        "human_nodes_indexed": human_nodes,
        "human_symbols_indexed": len(human_index),
        "orthology_pairs": sum(len(v) for v in orthology.values()),
        "node_comparison": {k: v for k, v in cmp.items() if k != "rows"},
        "node_comparison_symbol_identity": {k: v for k, v in cmp_symbol.items() if k != "rows"},
        "node_rows": cmp["rows"],
        "mouse_domains": rows,
        "evidence": EVIDENCE,
    }
    return out
