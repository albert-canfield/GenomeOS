# SPDX-License-Identifier: AGPL-3.0-or-later
"""Any person's genome, as a local individual: import a VCF once, then every layer sees it.

HG002 is the built-in test human. A geneticist's own data arrives the same way
HG002 does: one VCF against GRCh38 (plain or gzipped, any caller), streamed
once and split into one small PASS file per chromosome under
data/individuals/<name>/, with a manifest that records where it came from. From
then on `genomeos lookup` says whether that person carries a variant,
`genomeos report` lists their variants inside a gene, `genomeos twin build`
makes their haplotypes, and `genomeos individual genes` walks a chromosome gene
by gene with the consequence of every coding change, all offline.

The directory is git-ignored and nothing here sends a byte off the machine: a
personal genome stays where it was imported. Genotypes are measured evidence
(the caller's), consequences derived by the local trace.
"""

from __future__ import annotations

import gzip
import json
import re
import shutil
import time
from pathlib import Path
from typing import Any

ROOT = Path("data/individuals")
BUILTIN = {"HG002": "GIAB v4.2.1 benchmark calls (a real person, open consent)"}
CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,39}$")
HEADER = (
    "##fileformat=VCFv4.2\n##source={source}, {chrom} PASS subset via GenomeOS\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
)


def _open(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def normalise_chrom(raw: str) -> str | None:
    """`21`, `chr21`, `MT`, `M`, `X` → the hg38 name; anything else (contigs, decoys) → None."""
    c = raw if raw.startswith("chr") else f"chr{raw}"
    if c == "chrMT":
        c = "chrM"
    return c if c in CHROMOSOMES else None


def import_vcf(
    src: str | Path, name: str, note: str = "", root: Path | None = None, replace: bool = False, progress=None
) -> dict[str, Any]:
    """Stream one VCF and keep a PASS file per chromosome plus a manifest. Returns the manifest."""
    root = root or ROOT
    src = Path(src)
    if not _NAME.match(name):
        raise ValueError("name: letters, digits, _ . - only, up to 40 characters")
    if name in BUILTIN:
        raise ValueError(f"{name} is the built-in test human; pick another name")
    if not src.exists():
        raise FileNotFoundError(f"no VCF at {src}")
    d = root / name
    if d.exists():
        if not replace:
            raise FileExistsError(f"{name} already exists under {d}; pass replace to import again")
        shutil.rmtree(d)
    d.mkdir(parents=True)
    counts: dict[str, int] = {}
    handles: dict[str, Any] = {}
    sample = name
    skipped_filter = skipped_contig = 0
    phased = 0
    try:
        with _open(src) as fh:
            for line in fh:
                if line.startswith("#"):
                    if line.startswith("#CHROM"):
                        cols = line.rstrip("\n").split("\t")
                        if len(cols) >= 10:
                            sample = cols[9]
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 8:
                    continue
                if f[6] not in ("PASS", "."):
                    skipped_filter += 1
                    continue
                chrom = normalise_chrom(f[0])
                if chrom is None:
                    skipped_contig += 1
                    continue
                if chrom not in handles:
                    handles[chrom] = open(d / f"{name}_{chrom}.vcf", "w")  # noqa: SIM115
                    handles[chrom].write(HEADER.format(source=src.name, chrom=chrom, sample=name))
                    if progress:
                        progress(f"{chrom}: splitting {name}")
                gt = f[9].split(":")[0] if len(f) >= 10 else "1/1"
                phased += "|" in gt
                handles[chrom].write(f"{chrom}\t{f[1]}\t{f[2]}\t{f[3]}\t{f[4]}\t.\tPASS\t.\tGT\t{gt}\n")
                counts[chrom] = counts.get(chrom, 0) + 1
    finally:
        for h in handles.values():
            h.close()
    manifest = {
        "name": name,
        "source_file": src.name,
        "sample_column": sample,
        "note": note,
        "date": time.strftime("%Y-%m-%d"),
        "assembly": "GRCh38 assumed; positions are used as they are",
        "variants": sum(counts.values()),
        "phased": phased,
        "chromosomes": {c: counts[c] for c in CHROMOSOMES if c in counts},
        "skipped_filtered": skipped_filter,
        "skipped_other_contigs": skipped_contig,
        "evidence": f"measured: {src.name} (the caller's genotypes); nothing leaves this machine",
    }
    (d / "manifest.json").write_text(json.dumps(manifest, indent=1))
    return manifest


def remove(name: str, root: Path | None = None) -> bool:
    root = root or ROOT
    d = root / name
    if name in BUILTIN or not (d / "manifest.json").exists():
        return False
    shutil.rmtree(d)
    return True


def list_individuals(root: Path | None = None) -> list[dict[str, Any]]:
    """The built-in test human first, then every imported person, each with the chromosomes on disk."""
    from genomeos.genome.fetch import individual_vcf_path

    root = root or ROOT

    out = []
    have = [c for c in CHROMOSOMES if individual_vcf_path(c).exists()]
    out.append(
        {
            "name": "HG002",
            "builtin": True,
            "note": BUILTIN["HG002"],
            "chromosomes": have,
            "variants": None,
            "evidence": "measured: " + BUILTIN["HG002"],
        }
    )
    if root.exists():
        for m in sorted(root.glob("*/manifest.json")):
            try:
                d = json.loads(m.read_text())
            except (OSError, json.JSONDecodeError):
                continue
            d["builtin"] = False
            d["chromosomes"] = list(d.get("chromosomes", {}))
            out.append(d)
    return out


def vcf_path(name: str, chrom: str, root: Path | None = None) -> Path | None:
    """Where one person's rows for one chromosome live, or None if that person has none there."""
    root = root or ROOT
    if name == "HG002":
        from genomeos.genome.fetch import individual_vcf_path

        p = individual_vcf_path(chrom)
        return p if p.exists() else None
    p = root / name / f"{name}_{chrom}.vcf"
    return p if p.exists() else None


def sources(chrom: str, root: Path | None = None) -> list[tuple[str, Path, str]]:
    """(name, file, evidence) for every individual with rows on this chromosome, HG002 first."""
    root = root or ROOT
    out = []
    for ind in list_individuals(root):
        p = vcf_path(ind["name"], chrom, root)
        if p is not None:
            out.append((ind["name"], p, ind["evidence"]))
    return out


def rows_in(path: Path, start: int, end: int):
    """VCF rows with 1-based position inside start..end (1-based inclusive); files are sorted."""
    with path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            pos = int(f[1])
            if pos < start:
                continue
            if pos > end:
                break
            yield f


def gene_by_gene(
    name: str, chrom: str, annotation, genome, root: Path | None = None, limit: int = 0
) -> dict[str, Any]:
    """Walk a chromosome's coding genes with one person's variants: counts, and the consequence of
    every coding SNV on the canonical transcript (derived by the local trace)."""
    from genomeos.flow import trace_gene

    path = vcf_path(name, chrom, root)
    if path is None:
        raise FileNotFoundError(f"{name} has no rows on {chrom}")
    positions: list[tuple[int, list[str]]] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            positions.append((int(f[1]), f))
    positions.sort(key=lambda x: x[0])
    import bisect

    keys = [p for p, _ in positions]
    module = annotation.to_module("individual")
    genes = sorted(
        (g for g in annotation.genes.values() if g.locus.chrom == chrom and g.type == "protein_coding"),
        key=lambda g: g.locus.start,
    )
    rows = []
    totals = {"genes": 0, "genes_with_variants": 0, "variants_in_genes": 0, "coding_snvs": 0}
    by_consequence: dict[str, int] = {}
    for g in genes:
        lo = bisect.bisect_left(keys, g.locus.start + 1)
        hi = bisect.bisect_right(keys, g.locus.end)
        inside = positions[lo:hi]
        totals["genes"] += 1
        if not inside:
            continue
        totals["genes_with_variants"] += 1
        totals["variants_in_genes"] += len(inside)
        coding = []
        tr = trace_gene(genome, g, module.entities[g.id].transcripts)
        if tr is not None:
            for pos, f in inside:
                if len(f[3]) == 1 and len(f[4]) == 1:
                    sub = tr.substitute(pos - 1, f[3], f[4])
                    if sub.get("region") == "CDS":
                        c = sub.get("consequence") or "coding"
                        by_consequence[c] = by_consequence.get(c, 0) + 1
                        coding.append(
                            {
                                "pos": pos,
                                "change": f"{f[3]}>{f[4]}",
                                "genotype": f[9].split(":")[0] if len(f) > 9 else "",
                                "consequence": c,
                                "hgvs_p": sub.get("hgvs_p"),
                            }
                        )
        totals["coding_snvs"] += len(coding)
        rows.append(
            {
                "gene": g.symbol,
                "start": g.locus.start,
                "end": g.locus.end,
                "variants": len(inside),
                "coding_snvs": len(coding),
                "protein_changing": sum(
                    1 for c in coding if c["consequence"] not in ("synonymous", "coding")
                ),
                "changes": coding[:12],
            }
        )
    rows.sort(key=lambda r: (-r["protein_changing"], -r["coding_snvs"], -r["variants"]))
    return {
        "individual": name,
        "chrom": chrom,
        "totals": totals,
        "by_consequence": dict(sorted(by_consequence.items(), key=lambda kv: -kv[1])),
        "genes": rows[:limit] if limit else rows,
        "evidence": "measured genotypes; consequence derived by the local trace on the canonical transcript",
    }
