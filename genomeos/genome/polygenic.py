# SPDX-License-Identifier: AGPL-3.0-or-later
"""Polygenic scores for a person from the PGS Catalog's published weight tables.

A polygenic score is a weighted sum over a person's genotypes at a published set of variants. The
PGS Catalog (EBI) publishes each score's weights, harmonised to GRCh38, with the paper behind it and
the ancestry of the people it was built and tested on. This module streams one score's harmonised
table once (cached under data/knowledge/pgs, git-ignored), reads the person's genotypes chromosome by
chromosome and sums effect-allele dosages: a called variant gives its dosage; a position with no call
inside the person's trusted regions is homozygous reference, so the dosage is 2 when the effect
allele is the reference base and 0 otherwise; a position outside the trusted regions is missing and
counted, never guessed. The weights are `curated` evidence, the sum `derived`; the number is a raw
score with no reference distribution behind it, so it says where the person stands only against
other people scored the same way, and every score's ancestry caveat is carried with it. Results go
under the person's own directory, never under data/results.
"""

from __future__ import annotations

import gzip
import json
import time
import urllib.request
from pathlib import Path
from typing import Any

from genomeos.genome.individuals import (
    ROOT,
    _inside,
    _reference_base_reader,
    list_individuals,
    load_regions,
    vcf_path,
)

REST = "https://www.pgscatalog.org/rest/score/{pgs_id}"
FTP = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{pgs_id}/ScoringFiles/Harmonized/{pgs_id}_hmPOS_GRCh38.txt.gz"
KNOWLEDGE = Path("data/knowledge/pgs")
DEFAULT_SCORES = ("PGS000004", "PGS000115", "PGS000297", "PGS000018")  # breast cancer, LDL, height, CAD
EVIDENCE = (
    "curated: PGS Catalog weights (the cited publication); derived: the sum over this person's genotypes"
)
CONFIDENCE = 0.3


def metadata(pgs_id: str, knowledge: Path = KNOWLEDGE) -> dict[str, Any]:
    """Name, trait, publication, variant count, ancestry distributions and licence, cached."""
    p = knowledge / f"{pgs_id}.json"
    if p.exists():
        return json.loads(p.read_text())
    req = urllib.request.Request(REST.format(pgs_id=pgs_id), headers={"User-Agent": "GenomeOS/0.1"})
    with urllib.request.urlopen(req, timeout=120) as r:  # noqa: S310
        d = json.loads(r.read().decode())
    pub = d.get("publication") or {}
    anc = d.get("ancestry_distribution") or {}
    meta = {
        "id": d.get("id", pgs_id),
        "name": d.get("name"),
        "trait": d.get("trait_reported"),
        "variants": d.get("variants_number"),
        "publication": {
            "author": pub.get("firstauthor"),
            "journal": pub.get("journal"),
            "year": (pub.get("date_publication") or "")[:4],
            "doi": pub.get("doi"),
        },
        "ancestry": {
            stage: {k: v for k, v in (dist.get("dist") or {}).items()} for stage, dist in anc.items() if dist
        },
        "licence": d.get("license"),
    }
    knowledge.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(meta, indent=1))
    return meta


def scoring_file(pgs_id: str, knowledge: Path = KNOWLEDGE) -> Path:
    """The harmonised GRCh38 weight table, streamed once."""
    p = knowledge / f"{pgs_id}_hmPOS_GRCh38.txt.gz"
    if not p.exists():
        knowledge.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(FTP.format(pgs_id=pgs_id), headers={"User-Agent": "GenomeOS/0.1"})
        with urllib.request.urlopen(req, timeout=1800) as r:  # noqa: S310
            p.write_bytes(r.read())
    return p


def weights(fh) -> tuple[dict[str, str], list[tuple[str, int, str, str, float]]]:
    """(header fields, rows) of a scoring file: rows are (chrom, pos, effect, other, weight) on GRCh38."""
    header: dict[str, str] = {}
    rows: list[tuple[str, int, str, str, float]] = []
    cols: dict[str, int] | None = None
    for line in fh:
        if line.startswith("#"):
            if "=" in line and not line.startswith("##"):
                k, _, v = line[1:].rstrip("\n").partition("=")
                header[k] = v
            continue
        f = line.rstrip("\n").split("\t")
        if cols is None:
            cols = {name: i for i, name in enumerate(f)}
            continue
        chrom = f[cols["hm_chr"]] if "hm_chr" in cols else f[cols["chr_name"]]
        pos = f[cols["hm_pos"]] if "hm_pos" in cols else f[cols["chr_position"]]
        if not chrom or not pos:
            continue
        other = f[cols["other_allele"]] if "other_allele" in cols and cols["other_allele"] < len(f) else ""
        if not other and "hm_inferOtherAllele" in cols and cols["hm_inferOtherAllele"] < len(f):
            other = f[cols["hm_inferOtherAllele"]]
        try:
            w = float(f[cols["effect_weight"]])
        except ValueError:
            continue
        rows.append(
            (
                f"chr{chrom}" if not chrom.startswith("chr") else chrom,
                int(pos),
                f[cols["effect_allele"]],
                other,
                w,
            )
        )
    return header, rows


def _dosages(path: Path) -> dict[int, dict[str, int]]:
    """Position → {alt allele: dosage} and the ref allele under key "ref", from one chromosome's calls."""
    out: dict[int, dict[str, Any]] = {}
    with path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            gt = f[9].split(":")[0] if len(f) > 9 else ""
            alleles = gt.replace("|", "/").split("/")
            entry = out.setdefault(int(f[1]), {"ref": f[3], "alts": {}})
            for i, alt in enumerate(f[4].split(","), 1):
                entry["alts"][alt] = alleles.count(str(i))
    return out


def score_person(
    name: str, pgs_id: str, root: Path | None = None, knowledge: Path = KNOWLEDGE, progress=None
) -> dict[str, Any]:
    root = root or ROOT
    people = {p["name"]: p for p in list_individuals(root)}
    if name not in people:
        raise FileNotFoundError(f"{name} is not a local individual")
    meta = metadata(pgs_id, knowledge)
    t0 = time.time()
    with gzip.open(scoring_file(pgs_id, knowledge), "rt") as fh:
        header, rows = weights(fh)
    by_chrom: dict[str, list] = {}
    for r in rows:
        by_chrom.setdefault(r[0], []).append(r)
    total = 0.0
    n = called = assumed_ref = missing = mismatched = 0
    per_chrom: dict[str, dict[str, int]] = {}
    for chrom, items in by_chrom.items():
        vcf = vcf_path(name, chrom, root)
        if vcf is None:
            missing += len(items)
            per_chrom[chrom] = {"variants": len(items), "missing": len(items)}
            continue
        dos = _dosages(vcf)
        regions = load_regions(name, chrom, root)
        starts = [a for a, _ in regions]
        base_at, genome = _reference_base_reader(chrom)
        try:
            c = {"variants": 0, "called": 0, "assumed_reference": 0, "missing": 0, "allele_mismatch": 0}
            for _, pos, eff, _other, w in items:
                n += 1
                c["variants"] += 1
                entry = dos.get(pos)
                if entry is not None:
                    if eff in entry["alts"]:
                        d = entry["alts"][eff]
                    elif eff == entry["ref"]:
                        d = 2 - sum(entry["alts"].values())
                    else:
                        mismatched += 1
                        c["allele_mismatch"] += 1
                        continue
                    called += 1
                    c["called"] += 1
                elif not regions or _inside(regions, starts, pos - 1):
                    ref_base = base_at(pos) if base_at else None
                    if ref_base is None or len(eff) != 1:
                        missing += 1
                        c["missing"] += 1
                        continue
                    d = 2 if eff == ref_base else 0
                    assumed_ref += 1
                    c["assumed_reference"] += 1
                else:
                    missing += 1
                    c["missing"] += 1
                    continue
                total += w * d
            per_chrom[chrom] = c
        finally:
            if genome is not None:
                genome.close()
        if progress:
            progress(f"{pgs_id} {chrom}: {len(items):,} variants, {total:+.3f} so far")
    used = called + assumed_ref
    out = {
        "individual": name,
        "score": meta,
        "genome_build": header.get("HmPOS_build", "GRCh38"),
        "variants_in_score": len(rows),
        "variants_used": used,
        "called": called,
        "assumed_reference": assumed_ref,
        "missing": missing,
        "allele_mismatch": mismatched,
        "coverage": round(used / len(rows), 4) if rows else None,
        "raw_score": round(total, 6),
        "raw_score_per_used_variant": round(total / used, 8) if used else None,
        "per_chromosome": per_chrom,
        "seconds": round(time.time() - t0, 1),
        "evidence": EVIDENCE,
        "confidence": CONFIDENCE,
        "note": "a raw weighted sum, not a percentile: the catalog publishes weights, not the distribution "
        "of the score in any population, so this number places the person only against others scored the "
        "same way; the score was built and tested in the ancestries listed under score.ancestry, and its "
        "effect sizes transfer poorly outside them",
        "date": time.strftime("%Y-%m-%d"),
    }
    d = root / name
    d.mkdir(parents=True, exist_ok=True)
    (d / f"pgs_{pgs_id}.json").write_text(json.dumps(out, indent=1))
    return out


def compare(names: list[str], pgs_id: str, root: Path | None = None) -> list[dict[str, Any]]:
    """The raw scores of several local people on one score, highest first, with a z against the group."""
    root = root or ROOT
    rows = []
    for n in names:
        p = root / n / f"pgs_{pgs_id}.json"
        if p.exists():
            d = json.loads(p.read_text())
            rows.append({"individual": n, "raw_score": d["raw_score"], "coverage": d["coverage"]})
    if len(rows) >= 2:
        xs = [r["raw_score"] for r in rows]
        mu = sum(xs) / len(xs)
        sd = (sum((x - mu) ** 2 for x in xs) / (len(xs) - 1)) ** 0.5
        for r in rows:
            r["z_within_group"] = round((r["raw_score"] - mu) / sd, 3) if sd else 0.0
    rows.sort(key=lambda r: -r["raw_score"])
    return rows
