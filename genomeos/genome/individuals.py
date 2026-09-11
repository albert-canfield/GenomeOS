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


MISMATCH_TOLERANCE = (
    0.005  # more than 0.5% of calls disagreeing with the reference base means the wrong assembly
)


def verdict(stats: dict[str, Any]) -> str:
    """What the reference mismatches say about the file: the assembly check every import needs."""
    total = stats.get("variants", 0)
    mism = stats.get("reference_mismatches", 0)
    if total == 0:
        return "no variants to check"
    rate = mism / total
    if rate <= MISMATCH_TOLERANCE:
        return "matches GRCh38: the reference base agrees at the called positions"
    if rate >= 0.2:
        return (
            "does not match GRCh38: most calls disagree with the reference base "
            "(wrong assembly, e.g. GRCh37/hg19?)"
        )
    return "partly disagrees with GRCh38: check the assembly and the chromosome naming of the source file"


def check(
    name: str, chrom: str, root: Path | None = None, reference: Path = Path("data/reference")
) -> dict[str, Any]:
    """Apply one person's variants of one chromosome to the local reference, keep the statistics only
    (no haplotype FASTA), and say whether the file fits GRCh38. Stored under the person's directory."""
    from genomeos.coords import Locus
    from genomeos.genome import IndexedGenome, apply_variants, iter_vcf

    root = root or ROOT
    vcf = vcf_path(name, chrom, root)
    if vcf is None:
        raise FileNotFoundError(f"{name} has no rows on {chrom}")
    fa = reference / f"{chrom}.fa"
    if not fa.exists():
        raise FileNotFoundError(f"no reference sequence for {chrom} (genomeos data fetch --chrom {chrom})")
    t0 = time.time()
    variants = list(iter_vcf(vcf, {chrom}))
    genome = IndexedGenome(str(fa))
    try:
        ref = genome.fetch(Locus(chrom, 0, genome.lengths[chrom]))
    finally:
        genome.close()
    row: dict[str, Any] = {
        "individual": name,
        "chrom": chrom,
        "variants": len(variants),
        "snv": sum(1 for v in variants if v.is_snv),
        "phased": sum(1 for v in variants if v.phased),
    }
    for h in (0, 1):
        _, st = apply_variants(ref, variants, h)
        row[f"hap{h + 1}"] = st
    row["reference_mismatches"] = max(
        row["hap1"]["skipped_ref_mismatch"], row["hap2"]["skipped_ref_mismatch"]
    )
    row["verdict"] = verdict(row)
    row["seconds"] = round(time.time() - t0, 1)
    row["evidence"] = "measured genotypes applied to the local GRCh38 sequence; haplotype FASTA not kept"
    d = root / name
    d.mkdir(parents=True, exist_ok=True)
    (d / f"check_{chrom}.json").write_text(json.dumps(row, indent=1))
    return row


def checks(name: str, root: Path | None = None) -> dict[str, dict[str, Any]]:
    """Every stored reference check for a person, by chromosome."""
    root = root or ROOT
    out = {}
    d = root / name
    if d.exists():
        for p in sorted(d.glob("check_chr*.json")):
            try:
                out[p.stem[len("check_") :]] = json.loads(p.read_text())
            except (OSError, json.JSONDecodeError):
                continue
    return out


TRUNCATING = ("nonsense", "start_lost", "stop_lost")


def knockouts(
    name: str,
    chroms: list[str] | None = None,
    root: Path | None = None,
    reference: Path = Path("data/reference"),
):
    """Genome-wide: the person's SNVs that end, start or extend a canonical protein early or late
    (nonsense, start lost, stop lost), homozygous first: the natural knockouts. Derived by the local
    trace; SNVs only, so frameshift indels are not counted. Stored under the person's directory."""
    from genomeos.flow import trace_gene
    from genomeos.genome import Annotation, IndexedGenome, default_gencode

    root = root or ROOT
    people = {p["name"]: p for p in list_individuals(root)}
    if name not in people:
        raise FileNotFoundError(f"{name} is not a local individual")
    chroms = chroms or people[name]["chromosomes"]
    hits = []
    scanned_genes = 0
    after_reference_stop = 0
    done = []
    for chrom in chroms:
        vcf = vcf_path(name, chrom, root)
        gff = default_gencode({chrom})
        fa = reference / f"{chrom}.fa"
        if vcf is None or not gff or not fa.exists():
            continue
        ann = Annotation.from_gff3(gff, {chrom})
        genome = IndexedGenome(str(fa))
        try:
            module = ann.to_module("knockouts")
            positions: list[tuple[int, list[str]]] = []
            with vcf.open() as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    f = line.rstrip("\n").split("\t")
                    if len(f[3]) == 1 and len(f[4]) == 1:
                        positions.append((int(f[1]), f))
            positions.sort(key=lambda x: x[0])
            keys = [q for q, _ in positions]
            import bisect

            for g in ann.genes.values():
                if g.locus.chrom != chrom or g.type != "protein_coding":
                    continue
                lo = bisect.bisect_left(keys, g.locus.start + 1)
                hi = bisect.bisect_right(keys, g.locus.end)
                if lo == hi:
                    continue
                scanned_genes += 1
                tr = trace_gene(genome, g, module.entities[g.id].transcripts)
                if tr is None:
                    continue
                for pos, f in positions[lo:hi]:
                    sub = tr.substitute(pos - 1, f[3], f[4])
                    if sub.get("region") != "CDS" or sub.get("consequence") not in TRUNCATING:
                        continue
                    residue = sub.get("residue") or 0
                    length = len(tr.protein) or 1
                    if residue > length + 1:
                        # past a stop the reference itself carries: the transcript's CDS runs on after a
                        # premature stop (a reference nonsense allele or a
                        # pseudogene), not this person's doing
                        after_reference_stop += 1
                        continue
                    gt = f[9].split(":")[0] if len(f) > 9 else ""
                    alleles = gt.replace("|", "/").split("/")
                    hits.append(
                        {
                            "chrom": chrom,
                            "pos": pos,
                            "ref": f[3],
                            "alt": f[4],
                            "gene": g.symbol,
                            "transcript": tr.transcript,
                            "consequence": sub["consequence"],
                            "hgvs_p": sub.get("hgvs_p"),
                            "residue": residue,
                            "protein_length": length,
                            "fraction_lost": round(1 - residue / length, 3)
                            if sub["consequence"] == "nonsense"
                            else None,
                            "genotype": gt,
                            "zygosity": "homozygous" if alleles.count("1") >= 2 else "heterozygous",
                        }
                    )
        finally:
            genome.close()
        done.append(chrom)
    hits.sort(key=lambda h: (h["zygosity"] != "homozygous", -(h["fraction_lost"] or 0), h["chrom"], h["pos"]))
    out = {
        "individual": name,
        "chromosomes": done,
        "genes_with_variants": scanned_genes,
        "hits": hits,
        "nonsense": sum(1 for h in hits if h["consequence"] == "nonsense"),
        "start_lost": sum(1 for h in hits if h["consequence"] == "start_lost"),
        "stop_lost": sum(1 for h in hits if h["consequence"] == "stop_lost"),
        "homozygous": sum(1 for h in hits if h["zygosity"] == "homozygous"),
        "skipped_after_reference_stop": after_reference_stop,
        "genes": sorted({h["gene"] for h in hits}),
        "date": time.strftime("%Y-%m-%d"),
        "evidence": "measured genotypes; consequence derived by the local trace on the canonical transcript "
        "(SNVs only)",
        "note": "a truncating SNV is not a proven loss of function: late stops, alternative isoforms and "
        "nonsense-mediated decay escape all soften it; this is the list to look at, not a verdict",
    }
    d = root / name
    d.mkdir(parents=True, exist_ok=True)  # the test human gets a directory too (git-ignored like the rest)
    (d / "knockouts.json").write_text(json.dumps(out, indent=1))
    return out


def dossier(name: str, root: Path | None = None) -> str:
    """One Markdown page for a person from what has been computed and stored under their directory:
    the import, the reference checks, the ClinVar screen and the truncating variants. Nothing is
    recomputed here; a section that has not been run says so. Never written under data/results."""
    root = root or ROOT
    people = {p["name"]: p for p in list_individuals(root)}
    if name not in people:
        raise FileNotFoundError(f"{name} is not a local individual")
    p = people[name]
    d = root / name

    def load(fname: str) -> dict[str, Any] | None:
        f = d / fname
        try:
            return json.loads(f.read_text()) if f.exists() else None
        except (OSError, json.JSONDecodeError):
            return None

    lines = [f"# {name}", ""]
    src = (
        "built-in test human (GIAB HG002)" if p["builtin"] else f"imported from `{p.get('source_file', '?')}`"
    )
    note = f"; {p['note']}" if p.get("note") else ""
    variants = f"{p['variants']:,} PASS variants, " if p.get("variants") else ""
    lines.append(
        f"**Genome** — {src}{note}: {variants}{len(p['chromosomes'])} chromosomes on file. _{p['evidence']}_"
    )
    lines.append("")
    ck = checks(name, root)
    if ck:
        mism = sum(r["reference_mismatches"] for r in ck.values())
        bad = [c for c, r in ck.items() if not r["verdict"].startswith("matches")]
        state = "matches GRCh38" if not bad else f"assembly doubt on {', '.join(bad)}"
        applied = sum(r["hap1"]["applied"] + r["hap2"]["applied"] for r in ck.values())
        lines.append(
            f"**Reference check** — {len(ck)} chromosomes applied to GRCh38 ({applied:,} alleles placed), "
            f"{mism:,} reference mismatches: {state}. _measured genotypes against the local sequence_"
        )
    else:
        lines.append("**Reference check** — not run (`genomeos individual check`).")
    lines.append("")
    sc = load("clinvar_screen.json")
    if sc:
        lines.append(
            f"**ClinVar carrier screen** — {len(sc['hits'])} pathogenic or likely pathogenic alleles carried "
            f"over {sc['variants_scanned']:,} variants on {len(sc['chromosomes'])} chromosomes "
            f"({sc['pathogenic']} pathogenic, {sc['likely_pathogenic']} likely, {sc['homozygous']} "
            f"homozygous, {sc['two_stars_or_more']} with two or more review stars). _{sc['evidence']}_"
        )
        if sc["hits"]:
            lines.append("")
            lines.append("| variant | gene | genotype | significance | stars | condition |")
            lines.append("|---|---|---|---|---|---|")
            for h in sc["hits"][:40]:
                lines.append(
                    f"| {h['chrom']}:{h['pos']:,} {h['ref'][:8]}>{h['alt'][:8]} | "
                    f"{h['gene']} | {h['genotype']} "
                    f"({h['zygosity']}) | {h['significance']} | {h['stars']} | {h['conditions'][:60]} |"
                )
        lines.append("")
        lines.append(f"_{sc['note']}_")
    else:
        lines.append("**ClinVar carrier screen** — not run (`genomeos individual screen`).")
    lines.append("")
    cv = load("coding.json")
    if cv:
        bc = cv["by_consequence"]
        top = [r for r in cv["top"] if r["protein_changing"]][:15]
        lines.append(
            f"**Coding variants** — {sum(bc.values()):,} coding SNVs on canonical transcripts over "
            f"{len(cv['chromosomes'])} chromosomes: "
            + ", ".join(f"{k} {v:,}" for k, v in bc.items())
            + f"; {cv['genes_with_protein_changing']:,} genes carry a protein-changing variant, "
            f"{cv['genes_with_homozygous_changing']:,} a homozygous one. _{cv['evidence']}_"
        )
        if top:
            lines.append("")
            lines.append("Genes with the most protein-changing variants:")
            lines.append("")
            lines.append("| gene | protein-changing | homozygous | coding SNVs | examples |")
            lines.append("|---|---|---|---|---|")
            for r in top:
                lines.append(
                    f"| {r['gene']} | {r['protein_changing']} | "
                    f"{r['homozygous_changing']} | {r['coding_snvs']} | "
                    f"{'; '.join(r['examples'][:4])} |"
                )
        lines.append("")
        lines.append(f"_{cv['note']}_")
    else:
        lines.append("**Coding variants** — not run (`genomeos individual coding`).")
    lines.append("")
    ko = load("knockouts.json")
    if ko:
        lines.append(
            f"**Truncating variants** — {len(ko['hits'])} SNVs that end a canonical "
            f"protein early or remove its "
            f"start or stop, in {len(ko['genes'])} genes over {len(ko['chromosomes'])} chromosomes "
            f"({ko['nonsense']} nonsense, {ko['start_lost']} start lost, {ko['stop_lost']} stop lost; "
            f"{ko['homozygous']} homozygous). _{ko['evidence']}_"
        )
        hom = [h for h in ko["hits"] if h["zygosity"] == "homozygous"][:25]
        if hom:
            lines.append("")
            lines.append("Homozygous, most protein lost first:")
            lines.append("")
            lines.append("| gene | variant | consequence | protein | lost |")
            lines.append("|---|---|---|---|---|")
            for h in hom:
                lost = f"{h['fraction_lost']:.0%}" if h.get("fraction_lost") is not None else ""
                lines.append(
                    f"| {h['gene']} | {h['chrom']}:{h['pos']:,} {h['ref']}>{h['alt']} | "
                    f"{h['consequence'].replace('_', ' ')} | {h['hgvs_p'] or ''} "
                    f"of {h['protein_length']} aa | {lost} |"
                )
        lines.append("")
        lines.append(f"_{ko['note']}_")
    else:
        lines.append("**Truncating variants** — not run (`genomeos individual knockouts`).")
    lines.append("")
    lines.append(
        "_Every section is research annotation of a variant list against public references; none of it is "
        "clinical advice. This page lives under the person's own directory and is never committed._"
    )
    return "\n".join(lines)


def coding_inventory(
    name: str,
    chroms: list[str] | None = None,
    root: Path | None = None,
    reference: Path = Path("data/reference"),
):
    """Genome-wide: every coding SNV of the person on canonical transcripts, counted by consequence and by
    gene; the genes with the most protein-changing variants first. Derived by the local trace, SNVs only.
    Stored under the person's directory."""
    from genomeos.flow import trace_gene
    from genomeos.genome import Annotation, IndexedGenome, default_gencode

    root = root or ROOT
    people = {p["name"]: p for p in list_individuals(root)}
    if name not in people:
        raise FileNotFoundError(f"{name} is not a local individual")
    chroms = chroms or people[name]["chromosomes"]
    by_consequence: dict[str, int] = {}
    per_gene: list[dict[str, Any]] = []
    done = []
    genes_seen = 0
    for chrom in chroms:
        vcf = vcf_path(name, chrom, root)
        gff = default_gencode({chrom})
        fa = reference / f"{chrom}.fa"
        if vcf is None or not gff or not fa.exists():
            continue
        ann = Annotation.from_gff3(gff, {chrom})
        genome = IndexedGenome(str(fa))
        try:
            module = ann.to_module("coding")
            positions: list[tuple[int, list[str]]] = []
            with vcf.open() as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    f = line.rstrip("\n").split("\t")
                    if len(f[3]) == 1 and len(f[4]) == 1:
                        positions.append((int(f[1]), f))
            positions.sort(key=lambda x: x[0])
            keys = [q for q, _ in positions]
            import bisect

            for g in ann.genes.values():
                if g.locus.chrom != chrom or g.type != "protein_coding":
                    continue
                lo = bisect.bisect_left(keys, g.locus.start + 1)
                hi = bisect.bisect_right(keys, g.locus.end)
                if lo == hi:
                    continue
                tr = trace_gene(genome, g, module.entities[g.id].transcripts)
                if tr is None:
                    continue
                genes_seen += 1
                counts: dict[str, int] = {}
                hom_changing = 0
                examples = []
                for pos, f in positions[lo:hi]:
                    sub = tr.substitute(pos - 1, f[3], f[4])
                    if sub.get("region") != "CDS":
                        continue
                    if (sub.get("residue") or 0) > len(tr.protein) + 1:
                        continue  # past a stop the reference itself carries (see knockouts)
                    c = sub.get("consequence") or "coding"
                    counts[c] = counts.get(c, 0) + 1
                    by_consequence[c] = by_consequence.get(c, 0) + 1
                    if c not in ("synonymous", "coding"):
                        gt = f[9].split(":")[0] if len(f) > 9 else ""
                        if gt.replace("|", "/").split("/").count("1") >= 2:
                            hom_changing += 1
                        if len(examples) < 6:
                            examples.append(f"{sub.get('hgvs_p') or c} ({gt})")
                if counts:
                    changing = sum(v for k, v in counts.items() if k not in ("synonymous", "coding"))
                    per_gene.append(
                        {
                            "gene": g.symbol,
                            "chrom": chrom,
                            "coding_snvs": sum(counts.values()),
                            "protein_changing": changing,
                            "homozygous_changing": hom_changing,
                            "missense": counts.get("missense", 0),
                            "synonymous": counts.get("synonymous", 0),
                            "protein_length": len(tr.protein),
                            "examples": examples,
                        }
                    )
        finally:
            genome.close()
        done.append(chrom)
    per_gene.sort(key=lambda r: (-r["protein_changing"], -r["coding_snvs"]))
    out = {
        "individual": name,
        "chromosomes": done,
        "coding_genes_traced": genes_seen,
        "genes_with_coding_snvs": len(per_gene),
        "by_consequence": dict(sorted(by_consequence.items(), key=lambda kv: -kv[1])),
        "genes_with_protein_changing": sum(1 for r in per_gene if r["protein_changing"]),
        "genes_with_homozygous_changing": sum(1 for r in per_gene if r["homozygous_changing"]),
        "top": per_gene[:60],
        "date": time.strftime("%Y-%m-%d"),
        "evidence": "measured genotypes; consequence derived by the local trace on the canonical transcript "
        "(SNVs only)",
        "note": "a count of protein-changing variants per gene is not a measure of harm: long and "
        "polymorphic genes carry many; the ClinVar screen and the truncating list are where to look first",
    }
    d = root / name
    d.mkdir(parents=True, exist_ok=True)
    (d / "coding.json").write_text(json.dumps(out, indent=1))
    return out
