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


def _open(path: Path | str):
    """A local file (plain or gzip) or an http(s) URL streamed straight from the server."""
    s = str(path)
    if s.startswith(("http://", "https://")):
        import io
        import urllib.request

        req = urllib.request.Request(s, headers={"User-Agent": "GenomeOS/0.1 (stream)"})
        resp = urllib.request.urlopen(req, timeout=1800)  # noqa: S310
        raw = io.BufferedReader(resp, 1 << 20)
        return gzip.open(raw, "rt") if s.endswith(".gz") else io.TextIOWrapper(raw)
    return gzip.open(path, "rt") if s.endswith(".gz") else open(path)


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
    remote = str(src).startswith(("http://", "https://"))
    src = str(src) if remote else Path(src)
    if not _NAME.match(name):
        raise ValueError("name: letters, digits, _ . - only, up to 40 characters")
    if name in BUILTIN:
        raise ValueError(f"{name} is the built-in test human; pick another name")
    if not remote and not src.exists():
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
    source_name = src.rsplit("/", 1)[-1] if remote else src.name
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
                    handles[chrom].write(HEADER.format(source=source_name, chrom=chrom, sample=name))
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
        "source_file": source_name,
        "source_url": src if remote else None,
        "sample_column": sample,
        "note": note,
        "date": time.strftime("%Y-%m-%d"),
        "assembly": "GRCh38 assumed; positions are used as they are",
        "variants": sum(counts.values()),
        "phased": phased,
        "chromosomes": {c: counts[c] for c in CHROMOSOMES if c in counts},
        "skipped_filtered": skipped_filter,
        "skipped_other_contigs": skipped_contig,
        "evidence": f"measured: {source_name} (the caller's genotypes); nothing leaves this machine",
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
    pgs = [json.loads(f.read_text()) for f in sorted(d.glob("pgs_*.json")) if f.exists()]
    if pgs:
        lines.append("")
        lines.append("## Polygenic scores")
        lines.append("")
        lines.append(
            "Raw weighted sums over the person's genotypes from PGS Catalog weight tables, with the share of "
            "each score's variants read; not percentiles (the catalog publishes no population distribution), "
            "and only within the ancestries each score was built in. _weights curated; the sum derived_"
        )
        lines.append("")
        lines.append("| score | trait | publication | variants used | raw score |")
        lines.append("|---|---|---|---|---|")
        for r in pgs:
            m = r.get("score") or {}
            pub = m.get("publication") or {}
            who = f"{pub.get('author') or ''} {pub.get('year') or ''}".strip()
            used = f"{r.get('variants_used', 0):,} of {r.get('variants_in_score', 0):,}"
            raw = f"{r.get('raw_score', 0):+.4f}"
            lines.append(f"| {m.get('id', '')} | {m.get('trait') or ''} | {who} | {used} | {raw} |")
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
        pe = cv.get("missense_by_effect") or []
        if pe:
            c = cv["missense_predicted"]
            lines.append("")
            lines.append(
                f"Missense variants by predicted effect (AlphaMissense): {c['likely_pathogenic']:,} likely "
                f"pathogenic ({cv.get('missense_likely_pathogenic_homozygous', 0):,} homozygous), "
                f"{c['ambiguous']:,} ambiguous, {c['likely_benign']:,} likely benign, {c['unscored']:,} "
                f"unscored; strongest first. _{cv.get('missense_evidence', '')}_"
            )
            lines.append("")
            lines.append("| gene | change | genotype | score | class | UniProt says |")
            lines.append("|---|---|---|---|---|---|")
            for m in pe[:15]:
                pr = m.get("predicted") or {}
                lines.append(
                    f"| {m['gene']} | {m['hgvs_p'] or pr.get('protein_variant') or ''} | {m['genotype']} | "
                    f"{pr.get('score', '')} | {(pr.get('class') or 'unscored').replace('_', ' ')} | "
                    f"{'; '.join(m['features'][:2])} |"
                )
        mr = [m for m in cv.get("missense_ranked", []) if m["site"] == "site"][:15]
        if mr:
            lines.append("")
            lines.append(
                f"Missense variants on an annotated UniProt site "
                f"({cv.get('missense_at_annotated_site', 0)} of {cv.get('missense', 0)}; "
                f"{cv.get('missense_in_domain', 0)} more inside a domain), homozygous first:"
            )
            lines.append("")
            lines.append("| gene | change | genotype | site |")
            lines.append("|---|---|---|---|")
            for m in mr:
                lines.append(
                    f"| {m['gene']} | {m['hgvs_p'] or ''} | {m['genotype']} | "
                    f"{'; '.join(m['features'][:2])} |"
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
        rt = cv.get("reference_truncating")
        if rt:
            hidden = [r for r in cv.get("reference_truncating_rows", []) if r["matches_reference"]][:12]
            lines.append("")
            lines.append(
                f"**Reference truncating alleles** — hg38 itself carries a frameshift or nonsense allele in "
                f"{rt['genes']} genes against the curated protein; this person matches the reference in "
                f"{rt['matches_reference']} of the {rt['checked']} on file (so carries that truncation) "
                f"and has "
                f"variants inside the gene in {rt['has_variants_in_gene']}"
                + (
                    ": " + ", ".join(f"{r['gene']} ({r['kind'].split(' allele')[0]})" for r in hidden)
                    if hidden
                    else ""
                )
                + f". _{rt['evidence']}_"
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


SITE_TYPES = {
    "Active site",
    "Binding site",
    "Site",
    "Modified residue",
    "Disulfide bond",
    "Glycosylation",
    "Lipidation",
    "Cross-link",
    "DNA binding",
    "Metal binding",
    "Zinc finger",
    "Motif",
}
DOMAIN_TYPES = {
    "Domain",
    "Region",
    "Repeat",
    "Transmembrane",
    "Coiled coil",
    "Compositional bias",
    "Topological domain",
}


def _predicted_missense(
    name: str, missense: list[dict[str, Any]], root: Path, progress=None
) -> dict[str, Any]:
    """AlphaMissense scores for the person's missense variants (streamed once, cached under the person)."""
    from genomeos.genome import missense as am

    scores = am.scores_for(name, missense, root, progress)
    return am.annotate(missense, scores)


def site_class(features: list[dict[str, Any]]) -> str:
    """How much UniProt says about the residue a missense variant hits: a site, a domain, or nothing."""
    types = {f.get("type") for f in features}
    if types & SITE_TYPES:
        return "site"
    if types & DOMAIN_TYPES:
        return "domain"
    return "none"


def rank_missense(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Missense variants worth reading first: those on an annotated site, then in a domain, homozygous
    before heterozygous within each. A rank from annotation, not a prediction of effect."""
    order = {"site": 0, "domain": 1, "none": 2}
    return sorted(rows, key=lambda m: (order[m["site"]], m["zygosity"] != "homozygous", m["gene"], m["pos"]))


def coding_inventory(
    name: str,
    chroms: list[str] | None = None,
    root: Path | None = None,
    reference: Path = Path("data/reference"),
    predict: bool = False,
    progress=None,
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
    from genomeos.genome.lookup import features_at

    by_consequence: dict[str, int] = {}
    per_gene: list[dict[str, Any]] = []
    missense: list[dict[str, Any]] = []
    changing_all: list[dict[str, Any]] = []
    done = []
    genes_seen = 0
    protein_cache = Path("data/knowledge/proteins")
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
                defn = None
                dp = protein_cache / f"{g.symbol}.json"
                if dp.exists():
                    try:
                        defn = json.loads(dp.read_text())
                    except (OSError, json.JSONDecodeError):
                        defn = None
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
                        hom = gt.replace("|", "/").split("/").count("1") >= 2
                        if hom:
                            hom_changing += 1
                        changing_all.append(
                            {
                                "gene": g.symbol,
                                "chrom": chrom,
                                "pos": pos,
                                "ref": f[3],
                                "alt": f[4],
                                "consequence": c,
                                "hgvs_p": sub.get("hgvs_p"),
                                "genotype": gt,
                                "zygosity": "homozygous" if hom else "heterozygous",
                                "transcript": getattr(tr, "transcript", None),
                            }
                        )
                        if len(examples) < 6:
                            examples.append(f"{sub.get('hgvs_p') or c} ({gt})")
                        if c == "missense":
                            residue = sub.get("residue") or 0
                            feats = [
                                x
                                for x in (features_at(defn, residue) if defn else [])
                                # a bridge or a cross-link is two residues, not everything between them
                                if x["type"] not in ("Disulfide bond", "Cross-link")
                                or residue in (x["start"], x["end"])
                            ]
                            for x in feats:
                                if x["type"] in ("Disulfide bond", "Cross-link") and not x["description"]:
                                    x["description"] = f"Cys{x['start']}–Cys{x['end']}"
                            missense.append(
                                {
                                    "gene": g.symbol,
                                    "chrom": chrom,
                                    "pos": pos,
                                    "ref": f[3],
                                    "alt": f[4],
                                    "transcript": getattr(tr, "transcript", None),
                                    "hgvs_p": sub.get("hgvs_p"),
                                    "residue": sub.get("residue"),
                                    "protein_length": len(tr.protein),
                                    "genotype": gt,
                                    "zygosity": "homozygous" if hom else "heterozygous",
                                    "features": [
                                        f"{x['type']}: {x['description']}"[:60]
                                        for x in feats
                                        if x["type"] != "Chain"
                                    ][:4],
                                    "site": site_class(feats),
                                }
                            )
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
    ranked = rank_missense(missense)
    ref_trunc = reference_alleles_carried(name, root)
    out = {
        "reference_truncating": {k: v for k, v in ref_trunc.items() if k != "rows"},
        "reference_truncating_rows": [r for r in ref_trunc["rows"] if r["chromosome_on_file"]],
        "individual": name,
        "chromosomes": done,
        "coding_genes_traced": genes_seen,
        "genes_with_coding_snvs": len(per_gene),
        "by_consequence": dict(sorted(by_consequence.items(), key=lambda kv: -kv[1])),
        "genes_with_protein_changing": sum(1 for r in per_gene if r["protein_changing"]),
        "genes_with_homozygous_changing": sum(1 for r in per_gene if r["homozygous_changing"]),
        "missense": len(missense),
        "missense_at_annotated_site": sum(1 for m in missense if m["site"] == "site"),
        "missense_in_domain": sum(1 for m in missense if m["site"] == "domain"),
        "missense_ranked": ranked[:60],
        **(_predicted_missense(name, missense, root or ROOT, progress) if predict and missense else {}),
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
    (d / "coding_variants.json").write_text(json.dumps(changing_all))
    return out


def normalise_variant(pos: int, ref: str, alt: str, base_at=None) -> tuple[int, str, str]:
    """One representation for one allele: the common suffix trimmed, then the common prefix (one anchor
    base kept, as in a VCF), then the indel left-aligned along the reference while the base before it
    equals its last base. `base_at(pos1)` returns the reference base at a 1-based position; without it
    the trimming alone is applied. Three files that write TGG>TGGG, T>TG and TGG>T at one position are
    describing an insertion of G and a deletion of GG; after this they compare."""
    ref, alt = ref.upper(), alt.upper()
    if ref == alt:
        return pos, ref, alt
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        pos += 1
    if base_at is not None and len(ref) != len(alt) and ref[0] == alt[0]:
        # an indel with its anchor: shift left while the base before equals the indel's last base
        while pos > 1 and ref[-1] == alt[-1]:
            b = base_at(pos - 1)
            if not b:
                break
            ref, alt = b + ref[:-1], b + alt[:-1]
            pos -= 1
    return pos, ref, alt


def _genotypes(path: Path, base_at=None) -> dict[tuple[int, str, str], str]:
    out: dict[tuple[int, str, str], str] = {}
    with path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            gt = f[9].split(":")[0] if len(f) > 9 else ""
            alleles = gt.replace("|", "/").split("/")
            for i, alt in enumerate(f[4].split(","), 1):
                if str(i) in alleles:
                    key = (int(f[1]), f[3], alt)
                    if base_at is not None and (len(f[3]) > 1 or len(alt) > 1):
                        key = normalise_variant(*key, base_at=base_at)
                    zyg = "hom" if alleles.count(str(i)) >= 2 else "het"
                    out[key] = "hom" if out.get(key) == "hom" or zyg == "hom" else zyg
    return out


def _reference_base_reader(chrom: str, reference: Path = Path("data/reference")):
    """A base_at(pos1) over the local reference, or None when the chromosome is not fetched."""
    fa = reference / f"{chrom}.fa"
    if not fa.exists():
        return None, None
    from genomeos.coords import Locus
    from genomeos.genome import IndexedGenome

    g = IndexedGenome(str(fa))

    def base_at(pos1: int) -> str:
        return str(g.fetch(Locus(chrom, pos1 - 1, pos1))).upper()

    return base_at, g


def trio(
    child: str,
    father: str,
    mother: str,
    chroms: list[str] | None = None,
    root: Path | None = None,
    normalise: bool = True,
    reference: Path = Path("data/reference"),
):
    """Mendelian consistency of a child's variants against both parents, chromosome by chromosome:
    inherited from one or both, present in neither (a de novo candidate), and the impossible ones
    (homozygous in the child, absent from one parent). Where both parents have trusted regions
    (`import_regions`), a child call outside them counts as untrusted, not as an event. Autosomes."""
    root = root or ROOT
    people = {p["name"]: p for p in list_individuals(root)}
    for n in (child, father, mother):
        if n not in people:
            raise FileNotFoundError(f"{n} is not a local individual")
    chroms = chroms or [c for c in people[child]["chromosomes"] if c not in ("chrX", "chrY", "chrM")]
    keys = (
        "child_variants",
        "inherited",
        "in_both_parents",
        "de_novo_candidates",
        "mendelian_errors",
        "outside_a_parent_region",
        "de_novo_snv",
        "de_novo_indel",
    )
    tot = dict.fromkeys(keys, 0)
    normalised_any = False
    per_chrom = {}
    de_novo: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []
    done = []
    regions_used = False
    for chrom in chroms:
        pc, pf, pm = (vcf_path(n, chrom, root) for n in (child, father, mother))
        if not (pc and pf and pm):
            continue
        base_at, genome = _reference_base_reader(chrom, reference) if normalise else (None, None)
        try:
            gc, gf, gm = _genotypes(pc, base_at), _genotypes(pf, base_at), _genotypes(pm, base_at)
        finally:
            if genome is not None:
                genome.close()
        normalised_here = base_at is not None
        rf, rm = load_regions(father, chrom, root), load_regions(mother, chrom, root)
        rc = load_regions(child, chrom, root)  # the child's own trusted regions, when it has them
        sf, sm, sc = [x for x, _ in rf], [x for x, _ in rm], [x for x, _ in rc]
        have_regions = bool(rf and rm)
        regions_used = regions_used or have_regions
        row = dict.fromkeys(keys, 0)
        for key, zyg in gc.items():
            row["child_variants"] += 1
            f_has, m_has = key in gf, key in gm
            pos0 = key[0] - 1
            trusted = (not have_regions) or (
                _inside(rf, sf, pos0) and _inside(rm, sm, pos0) and (not rc or _inside(rc, sc, pos0))
            )
            if f_has and m_has:
                row["in_both_parents"] += 1
                row["inherited"] += 1
            elif f_has or m_has:
                if zyg != "hom":
                    row["inherited"] += 1
                elif trusted:
                    row["mendelian_errors"] += 1
                    if len(errors) < 200:
                        errors.append(
                            {
                                "chrom": chrom,
                                "pos": key[0],
                                "ref": key[1],
                                "alt": key[2],
                                "child": zyg,
                                "father": gf.get(key),
                                "mother": gm.get(key),
                            }
                        )
                else:
                    row["outside_a_parent_region"] += 1
            elif trusted:
                row["de_novo_candidates"] += 1
                row["de_novo_snv" if len(key[1]) == 1 and len(key[2]) == 1 else "de_novo_indel"] += 1
                if len(de_novo) < 500:
                    de_novo.append(
                        {"chrom": chrom, "pos": key[0], "ref": key[1], "alt": key[2], "child": zyg}
                    )
            else:
                row["outside_a_parent_region"] += 1
        for k in keys:
            tot[k] += row[k]
        per_chrom[chrom] = row
        done.append(chrom)
        normalised_any = normalised_any or normalised_here
    cv = tot["child_variants"] or 1
    out = {
        "child": child,
        "father": father,
        "mother": mother,
        "chromosomes": done,
        "totals": tot,
        "fraction_inherited": round(tot["inherited"] / cv, 4),
        "fraction_de_novo_candidates": round(tot["de_novo_candidates"] / cv, 5),
        "fraction_mendelian_errors": round(tot["mendelian_errors"] / cv, 5),
        "regions": (
            "trusted regions applied (both parents, and the child's where present)"
            if regions_used
            else "no trusted regions: every absence counts, most are no-calls"
        ),
        "representation": (
            "alleles normalised: common suffix and prefix trimmed, indels left-aligned on the reference"
            if normalised_any
            else "alleles compared as written"
        ),
        "per_chromosome": per_chrom,
        "de_novo_candidates": de_novo,
        "mendelian_errors": errors,
        "date": time.strftime("%Y-%m-%d"),
        "evidence": "measured genotypes of three people; inheritance derived by comparing calls at the same "
        "position and alleles",
        "note": "a de novo candidate is a child call absent from both parents' files inside both parents' "
        "trusted regions: real de novo variation is about 60 to 100 per genome, so a larger count is "
        "representation differences and residual no-calls; a Mendelian error is a homozygous child call "
        "with a parent lacking the allele, inside the trusted regions",
    }
    d = root / child
    d.mkdir(parents=True, exist_ok=True)
    (d / f"trio_{father}_{mother}.json").write_text(json.dumps(out, indent=1))
    return out


def import_regions(name: str, src: str | Path, root: Path | None = None, progress=None) -> dict[str, int]:
    """The regions a person's calls are trusted in (a BED, file or URL: GIAB's benchmark regions, a
    caller's callable-regions track), kept per chromosome under the person's directory. Without them a
    variant absent from a parent cannot be told from a parent that was never called there."""
    root = root or ROOT
    d = root / name
    if not (d / "manifest.json").exists() and name not in BUILTIN:
        raise FileNotFoundError(f"{name} is not a local individual")
    d.mkdir(parents=True, exist_ok=True)
    counts: dict[str, int] = {}
    handles: dict[str, Any] = {}
    try:
        with _open(src) as fh:
            for line in fh:
                if line.startswith(("#", "track", "browser")):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 3:
                    continue
                chrom = normalise_chrom(f[0])
                if chrom is None:
                    continue
                if chrom not in handles:
                    handles[chrom] = open(d / f"regions_{chrom}.bed", "w")  # noqa: SIM115
                handles[chrom].write(f"{chrom}\t{f[1]}\t{f[2]}\n")
                counts[chrom] = counts.get(chrom, 0) + 1
    finally:
        for h in handles.values():
            h.close()
    m = d / "manifest.json"
    if m.exists():
        try:
            man = json.loads(m.read_text())
            man["regions_source"] = str(src).rsplit("/", 1)[-1]
            man["regions_intervals"] = sum(counts.values())
            m.write_text(json.dumps(man, indent=1))
        except (OSError, json.JSONDecodeError):
            pass
    if progress:
        progress(f"{name}: {sum(counts.values()):,} trusted intervals on {len(counts)} chromosomes")
    return counts


def load_regions(name: str, chrom: str, root: Path | None = None) -> list[tuple[int, int]]:
    """Sorted, merged trusted intervals (0-based, half-open) of a person on a chromosome; empty if none."""
    root = root or ROOT
    p = root / name / f"regions_{chrom}.bed"
    if not p.exists():
        return []
    iv = []
    with p.open() as fh:
        for line in fh:
            f = line.split("\t")
            iv.append((int(f[1]), int(f[2])))
    iv.sort()
    out: list[tuple[int, int]] = []
    for a, b in iv:
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def _inside(intervals: list[tuple[int, int]], starts: list[int], pos0: int) -> bool:
    import bisect

    i = bisect.bisect_right(starts, pos0) - 1
    return i >= 0 and intervals[i][0] <= pos0 < intervals[i][1]


def variants_on_transcript(tr, rows) -> list[dict[str, Any]]:
    """A person's SNVs (VCF rows) read on one traced transcript: the coding ones with their consequence."""
    out = []
    for f in rows:
        if len(f[3]) != 1 or len(f[4]) != 1:
            continue
        sub = tr.substitute(int(f[1]) - 1, f[3], f[4])
        if sub.get("region") != "CDS":
            continue
        gt = f[9].split(":")[0] if len(f) > 9 else ""
        out.append(
            {
                "pos": int(f[1]),
                "ref": f[3],
                "alt": f[4],
                "genotype": gt,
                "consequence": sub.get("consequence"),
                "hgvs_p": sub.get("hgvs_p"),
                "residue": sub.get("residue"),
            }
        )
    return out


def tissue_proteins(
    name: str, gene: str, chrom: str, root: Path | None = None, reference: Path = Path("data/reference")
):
    """Which protein each tissue makes in this person: GTEx's dominant transcript per tissue, traced,
    with the person's coding variants read on that transcript rather than on the canonical one."""
    from genomeos.flow.trace import trace
    from genomeos.genome import Annotation, IndexedGenome, default_gencode
    from genomeos.molecules.rna import gtex_isoforms

    root = root or ROOT
    vcf = vcf_path(name, chrom, root)
    gff = default_gencode({chrom})
    fa = reference / f"{chrom}.fa"
    if vcf is None:
        raise FileNotFoundError(f"{name} has no rows on {chrom}")
    if not gff or not fa.exists():
        raise FileNotFoundError(
            f"{chrom} needs local models and sequence (genomeos data fetch --chrom {chrom})"
        )
    ann = Annotation.from_gff3(gff, {chrom})
    g = ann.gene(gene.upper())
    module = ann.to_module("tissue")
    txs = {t.id.split(".")[0]: t for t in module.entities[g.id].transcripts if t.cds_segments}
    defn_p = Path("data/knowledge/proteins") / f"{g.symbol}.json"
    labels = []
    if defn_p.exists():
        try:
            d = json.loads(defn_p.read_text())
            labels = ((d["sections"].get("genomic_origin") or {}).get("items") or {}).get("transcripts") or []
        except (OSError, json.JSONDecodeError):
            labels = []
    iso = gtex_isoforms(g.symbol, transcripts=labels)
    rows = list(rows_in(vcf, g.locus.start + 1, g.locus.end))
    genome = IndexedGenome(str(fa))
    proteins = []
    try:
        dominant = iso.get("dominant_by_tissue") or {}
        by_tx: dict[str, list[str]] = {}
        for tissue, d in dominant.items():
            by_tx.setdefault(d["transcript"], []).append(tissue)
        canonical = iso.get("canonical")
        if canonical and canonical not in by_tx:
            by_tx[canonical] = []
        for tid, tissues in sorted(by_tx.items(), key=lambda kv: -len(kv[1])):
            tx = txs.get(tid)
            if tx is None:
                proteins.append(
                    {"transcript": tid, "tissues": tissues, "note": "non-coding in the local models"}
                )
                continue
            tr = trace(genome, tx, g.symbol)
            vs = variants_on_transcript(tr, rows)
            changing = [v for v in vs if v["consequence"] not in ("synonymous", "coding")]
            proteins.append(
                {
                    "transcript": tid,
                    "name": (iso.get("isoforms", {}).get(tid) or {}).get("name"),
                    "canonical": tid == canonical,
                    "tissues_dominant": len(tissues),
                    "tissues": tissues[:12],
                    "protein_length": len(tr.protein),
                    "coding_snvs": len(vs),
                    "protein_changing": [
                        {k: v[k] for k in ("hgvs_p", "genotype", "consequence", "residue")} for v in changing
                    ],
                }
            )
    finally:
        genome.close()
    return {
        "individual": name,
        "gene": g.symbol,
        "chrom": chrom,
        "tissues_measured": iso.get("tissues_measured", 0),
        "proteins": proteins,
        "evidence": "measured: GTEx v8 transcript TPM per tissue; measured genotypes; consequence derived "
        "by the local trace on the tissue's dominant transcript",
        "note": "the isoform a tissue makes is GTEx's population median, not this person's own expression",
    }


def reference_truncating_genes(results_dir: Path = Path("data/results")) -> list[dict[str, Any]]:
    """The genes where hg38 itself carries a frameshift or nonsense allele against the curated protein
    (the verified translation disagreements triaged by mechanism): a person who matches the reference
    there carries a truncation the reference hides."""
    out = []
    for f in sorted(results_dir.glob("translation_vs_uniprot_chr*.json")):
        try:
            d = json.loads(f.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        chrom = d.get("chrom")
        if chrom == "chrM":
            continue
        for r in d.get("disagreements", []):
            frac = r.get("translated_fraction_of_cds")
            if r.get("cds_out_of_frame"):
                kind = "frameshift allele in hg38"
            elif frac is not None and frac < 0.7 and r["ours_aa"] < 0.7 * r["uniprot_aa"]:
                kind = "nonsense allele in hg38"
            else:
                continue
            out.append(
                {
                    "gene": r["gene"],
                    "chrom": chrom,
                    "kind": kind,
                    "reference_protein_aa": r["ours_aa"],
                    "curated_aa": r["uniprot_aa"],
                }
            )
    return out


def reference_alleles_carried(name: str, root: Path | None = None, annotation_for=None) -> dict[str, Any]:
    """For each reference truncating allele: does the person differ from the reference inside that gene's
    CDS at all? No variant inside the CDS means the person matches hg38 there, hence carries the
    truncation the curated protein does not have."""
    root = root or ROOT
    genes = reference_truncating_genes()
    rows = []
    by_chrom: dict[str, list[dict]] = {}
    for g in genes:
        by_chrom.setdefault(g["chrom"], []).append(g)
    for chrom, gs in sorted(by_chrom.items()):
        vcf = vcf_path(name, chrom, root)
        ann = annotation_for(chrom) if annotation_for else None
        if ann is None:
            from genomeos.genome import Annotation, default_gencode

            gff = default_gencode({chrom})
            if not gff:
                continue
            ann = Annotation.from_gff3(gff, {chrom})
        for g in gs:
            try:
                gene = ann.gene(g["gene"])
            except KeyError:
                continue
            inside = list(rows_in(vcf, gene.locus.start + 1, gene.locus.end)) if vcf else []
            snvs = sum(1 for f in inside if len(f[3]) == 1 and len(f[4]) == 1)
            rows.append(
                {
                    **g,
                    "variants_in_gene": len(inside),
                    "snvs_in_gene": snvs,
                    "indels_in_gene": len(inside) - snvs,
                    "matches_reference": not inside,
                    "chromosome_on_file": vcf is not None,
                }
            )
    on_file = [r for r in rows if r["chromosome_on_file"]]
    return {
        "individual": name,
        "genes": len(genes),
        "checked": len(on_file),
        "matches_reference": sum(1 for r in on_file if r["matches_reference"]),
        "has_variants_in_gene": sum(1 for r in on_file if not r["matches_reference"]),
        "rows": rows,
        "evidence": "derived: hg38 translation against UniProt (verified per chromosome); measured genotypes",
        "note": "a variant inside the gene does not restore the curated protein on its own (none of the "
        "GIAB trio's do); no variant at all means the person carries hg38's allele there, a truncation the "
        "reference hides",
    }
