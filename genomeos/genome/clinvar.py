# SPDX-License-Identifier: AGPL-3.0-or-later
"""ClinVar carrier screen for a local individual: stream ClinVar once, keep the pathogenic rows, intersect.

ClinVar's GRCh38 VCF (about 190 MB, weekly) is streamed once and distilled to the rows whose
clinical significance is pathogenic or likely pathogenic (no conflicting calls): chromosome,
position, alleles, gene, significance, conditions, review status. That file stays local under
data/knowledge/clinvar; only its counts are committed. A screen then walks a person's per-chromosome
files and reports every ClinVar pathogenic allele the person carries, with the genotype, so a
geneticist sees the recessive carrier state and the rare dominant hit in one table, offline.

This is research annotation of a variant list against a public database, with ClinVar's own
review status carried through; it is not a clinical report, and the output says so. Per-person
results stay under the person's git-ignored directory.
"""

from __future__ import annotations

import gzip
import io
import json
import time
import urllib.request
from pathlib import Path
from typing import Any

CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
KNOWLEDGE = Path("data/knowledge/clinvar")
RESULTS = Path("data/results")
KEEP = ("Pathogenic", "Likely_pathogenic")
EVIDENCE = "curated: ClinVar (NCBI) pathogenic and likely pathogenic assertions, with their review status"
NOTE = "research annotation against a public database, not a clinical report; review status is ClinVar's own"
STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 0,
    "no_assertion_criteria_provided": 0,
    "no_classification_provided": 0,
}


def pathogenic_path(knowledge: Path = KNOWLEDGE) -> Path:
    return knowledge / "pathogenic.tsv.gz"


def _info(field: str) -> dict[str, str]:
    out = {}
    for kv in field.split(";"):
        k, _, v = kv.partition("=")
        out[k] = v
    return out


def distil(knowledge: Path = KNOWLEDGE, url: str = CLINVAR_URL, progress=None) -> dict[str, Any]:
    """Stream the VCF once; keep pathogenic and likely pathogenic rows; return and save the counts."""
    dest = pathogenic_path(knowledge)
    knowledge.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "GenomeOS/0.1 (stream)"})
    kept = 0
    seen = 0
    by_sig: dict[str, int] = {}
    by_chrom: dict[str, int] = {}
    t0 = time.time()
    with (
        urllib.request.urlopen(req, timeout=1800) as resp,  # noqa: S310
        gzip.open(io.BufferedReader(resp, 1 << 20), "rt") as fh,
        gzip.open(dest, "wt") as out,
    ):
        out.write("# ClinVar GRCh38, pathogenic and likely pathogenic rows distilled by GenomeOS\n")
        out.write("#chrom\tpos\tref\talt\tgene\tsignificance\tconditions\treview\tstars\tclinvar_id\trs\n")
        for line in fh:
            if line.startswith("#"):
                continue
            seen += 1
            if progress and seen % 500_000 == 0:
                progress(f"ClinVar: {seen:,} rows read, {kept:,} kept")
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            info = _info(f[7])
            sig = info.get("CLNSIG", "")
            if "Conflicting" in sig or not any(k in sig for k in KEEP):
                continue
            chrom = f[0] if f[0].startswith("chr") else ("chrM" if f[0] == "MT" else f"chr{f[0]}")
            gene = info.get("GENEINFO", "").split("|")[0].split(":")[0]
            conds = [
                c.replace("_", " ") for c in info.get("CLNDN", "").split("|") if c and c != "not_provided"
            ]
            review = info.get("CLNREVSTAT", "")
            out.write(
                f"{chrom}\t{f[1]}\t{f[3]}\t{f[4]}\t{gene}\t{sig}\t{'; '.join(conds[:3])}\t{review}\t"
                f"{STARS.get(review, 0)}\t{f[2]}\t{info.get('RS', '')}\n"
            )
            kept += 1
            key = (
                "pathogenic"
                if sig.startswith("Pathogenic")
                else "likely_pathogenic"
                if sig.startswith("Likely")
                else "other"
            )
            by_sig[key] = by_sig.get(key, 0) + 1
            by_chrom[chrom] = by_chrom.get(chrom, 0) + 1
    summary = {
        "source": url,
        "rows_read": seen,
        "rows_kept": kept,
        "by_significance": by_sig,
        "by_chromosome": by_chrom,
        "seconds": round(time.time() - t0),
        "kept_where": str(dest),
        "evidence": EVIDENCE,
    }
    from genomeos.results import save_result

    save_result("clinvar_pathogenic_summary", summary)
    return summary


def load_chromosome(chrom: str, knowledge: Path = KNOWLEDGE) -> dict[tuple[int, str, str], dict[str, Any]]:
    """Pathogenic rows of one chromosome keyed by (pos, ref, alt)."""
    p = pathogenic_path(knowledge)
    out: dict[tuple[int, str, str], dict[str, Any]] = {}
    if not p.exists():
        return out
    with gzip.open(p, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] != chrom:
                continue
            out[(int(f[1]), f[2], f[3])] = {
                "gene": f[4],
                "significance": f[5].replace("_", " "),
                "conditions": f[6],
                "review": f[7].replace("_", " "),
                "stars": int(f[8]),
                "clinvar_id": f[9],
                "rs": f"rs{f[10]}" if f[10] else "",
            }
    return out


def screen(name: str, chroms: list[str] | None = None, root: Path | None = None, knowledge: Path = KNOWLEDGE):
    """Every ClinVar pathogenic allele the person carries, chromosome by chromosome."""
    from genomeos.genome.individuals import ROOT, list_individuals, vcf_path

    root = root or ROOT
    people = {p["name"]: p for p in list_individuals(root)}
    if name not in people:
        raise FileNotFoundError(f"{name} is not a local individual")
    chroms = chroms or people[name]["chromosomes"]
    hits = []
    scanned = 0
    for chrom in chroms:
        vcf = vcf_path(name, chrom, root)
        if vcf is None:
            continue
        table = load_chromosome(chrom, knowledge)
        if not table:
            continue
        with vcf.open() as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                scanned += 1
                pos = int(f[1])
                gt = f[9].split(":")[0] if len(f) > 9 else ""
                alleles = gt.replace("|", "/").split("/")
                for i, alt in enumerate(f[4].split(","), 1):
                    row = table.get((pos, f[3], alt))
                    if row and (str(i) in alleles or not alleles):
                        hits.append(
                            {
                                "chrom": chrom,
                                "pos": pos,
                                "ref": f[3],
                                "alt": alt,
                                "genotype": gt,
                                "zygosity": "homozygous" if alleles.count(str(i)) >= 2 else "heterozygous",
                                **row,
                            }
                        )
    hits.sort(key=lambda h: (-h["stars"], h["significance"] != "Pathogenic", h["chrom"], h["pos"]))
    out = {
        "individual": name,
        "chromosomes": chroms,
        "variants_scanned": scanned,
        "hits": hits,
        "pathogenic": sum(1 for h in hits if h["significance"].startswith("Pathogenic")),
        "likely_pathogenic": sum(1 for h in hits if h["significance"].startswith("Likely")),
        "homozygous": sum(1 for h in hits if h["zygosity"] == "homozygous"),
        "two_stars_or_more": sum(1 for h in hits if h["stars"] >= 2),
        "date": time.strftime("%Y-%m-%d"),
        "evidence": EVIDENCE,
        "note": NOTE,
    }
    d = root / name
    if d.exists():
        (d / "clinvar_screen.json").write_text(json.dumps(out, indent=1))
    return out
