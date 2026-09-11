"""Storage economics: status, distil, clean.

    genomeos data status        disk use per data directory, largest files
    genomeos data distil        turn raw downloads into small result summaries
    genomeos data clean         delete raw inputs that already have a summary

Nothing under data/results/ is ever deleted. `clean` only removes files whose
distilled summary exists, so no test loses the information it needs.
"""

from __future__ import annotations

import json
import shutil
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from .results import load_result, save_result

DATA = Path("data")


@dataclass(slots=True)
class Distiller:
    name: str  # result name
    inputs: list[Path]  # raw files consumed
    run: Callable[[], dict]  # produces the summary payload
    describe: str = ""


def _size(path: Path) -> int:
    if path.is_file():
        return path.stat().st_size
    return sum(p.stat().st_size for p in path.rglob("*") if p.is_file())


def human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"


def status(data: Path = DATA) -> dict:
    dirs = {p.name: _size(p) for p in sorted(data.iterdir()) if p.is_dir()} if data.is_dir() else {}
    files = sorted((p for p in data.rglob("*") if p.is_file()), key=lambda p: -p.stat().st_size)[:12]
    total, used, free = shutil.disk_usage(data if data.exists() else Path("."))
    return {
        "dirs": dirs,
        "largest": [(str(p), p.stat().st_size) for p in files],
        "disk": {"total": total, "used": used, "free": free},
        "results": len(list((data / "results").glob("*.json"))) if (data / "results").is_dir() else 0,
    }


# ---- distillers ------------------------------------------------------------


def _distil_clock_geo() -> dict:
    from genomeos.twin.clocks import Clock, MethylationMatrix, pearson

    geo = DATA / "knowledge" / "GSE41169_series_matrix.txt.gz"
    h, hn = Clock.horvath(), Clock.hannum()
    m = MethylationMatrix.from_geo_series_matrix(geo, cpgs=set(h.coefficients) | set(hn.coefficients))
    rows = []
    for s in m.samples:
        age = m.metadata.get(s, {}).get("age")
        if not age:
            continue
        betas = m.sample(s)
        rows.append(
            {
                "sample": s,
                "age": float(age),
                "horvath": round(h.predict(betas).age, 2),
                "hannum": round(hn.predict(betas).age, 2),
            }
        )
    ages = [r["age"] for r in rows]
    out = {
        "dataset": "GEO GSE41169 whole blood (Horvath 2013 test set)",
        "samples": len(rows),
        "per_sample": rows,
    }
    for clock in ("horvath", "hannum"):
        preds = [r[clock] for r in rows]
        out[clock] = {
            "r": round(pearson(ages, preds), 4),
            "mae_years": round(sum(abs(a - p) for a, p in zip(ages, preds, strict=True)) / len(rows), 2),
        }
    return out


def _distil_clinvar() -> dict:
    from collections import Counter

    from genomeos.genome import Annotation, IndexedGenome, Variant, default_gencode, iter_vcf
    from genomeos.runtime import classify_all

    ann = Annotation.from_gff3(default_gencode({"chr21"}), {"chr21"})
    m = ann.to_module("chr21")
    genome = IndexedGenome(DATA / "reference" / "chr21.fa.gz")
    coding = {}
    for g in ann.protein_coding():
        txs = [
            t
            for t in m.entities[g.id].transcripts
            if t.cds_segments and t.attrs["transcript_type"] == "protein_coding"
        ]
        if txs:
            coding[g.symbol] = txs
    wanted = {"nonsense", "missense_variant", "synonymous_variant", "frameshift_variant"}
    compatible = {
        "nonsense": {"nonsense", "frameshift_variant"},
        "frameshift_variant": {"frameshift_variant"},
        "missense_variant": {"missense_variant"},
        "synonymous_variant": {"synonymous_variant"},
    }
    total: Counter[str] = Counter()
    agree: Counter[str] = Counter()
    confusion: Counter[str] = Counter()
    clinvar = DATA / "reference" / "clinvar_chr21_MT.vcf"
    if not clinvar.exists():
        clinvar = DATA / "results" / "clinvar_chr21_coding.vcf.gz"
    for v in iter_vcf(clinvar, {"21"}, pass_only=False):
        info = dict(kv.split("=", 1) for kv in v.info.split(";") if "=" in kv)
        labels = {x.split("|")[1] for x in info.get("MC", "").split(",") if "|" in x}
        if len(labels) != 1 or not labels & wanted or len(v.alts) != 1:
            continue
        label = labels.pop()
        genes = [x.split(":")[0] for x in info.get("GENEINFO", "").split("|")]
        txs = [t for g in genes for t in coding.get(g, [])]
        if not txs:
            continue
        got = classify_all(genome, txs, Variant("chr21", v.pos, v.ref, v.alts, gt=(1, 1)))[0].consequence
        total[label] += 1
        if got in compatible[label]:
            agree[label] += 1
        else:
            confusion[f"{label}->{got}"] += 1
    genome.close()
    return {
        "dataset": "ClinVar chr21 (all records with a single molecular consequence)",
        "agreement": {
            k: {"n": total[k], "agree": agree[k], "rate": round(agree[k] / total[k], 4)} for k in total
        },
        "confusion_top": confusion.most_common(10),
    }


def _distil_libraries() -> dict:
    from genomeos.lib import LIBRARIES, KnowledgeBase

    kb = KnowledgeBase()
    members = {lib.id: sorted(kb.members(lib)) for lib in LIBRARIES.values() if lib.go_terms or lib.reactome}
    ver = kb.verify_all()
    return {
        "source": "GO (goa_human) + Reactome (Ensembl2Reactome, all levels) + GENCODE 50 symbols",
        "libraries": len(members),
        "members": members,
        "verification": [
            {
                "library": v.library,
                "members": v.members,
                "checked": len(v.checked),
                "missing": list(v.missing),
            }
            for v in ver
        ],
        "agreement": round(1 - sum(len(v.missing) for v in ver) / sum(len(v.checked) for v in ver), 4),
    }


def _distil_hg002_chr21() -> dict:
    from genomeos.genome import IndexedGenome, Locus, apply_variants, iter_vcf

    variants = list(iter_vcf(DATA / "reference" / "HG002_GRCh38_benchmark.vcf.gz", {"chr21"}))
    ref = IndexedGenome(DATA / "reference" / "chr21.fa.gz").fetch(Locus("chr21", 0, 46_709_983))
    out = {
        "sample": "HG002 GIAB v4.2.1 benchmark, chr21",
        "variants": len(variants),
        "snv": sum(1 for v in variants if v.is_snv),
        "phased": sum(1 for v in variants if v.phased),
    }
    for h in (0, 1):
        _, st = apply_variants(ref, variants, h)
        out[f"hap{h + 1}"] = st
    # keep the chr21 records themselves: a few MB, enough to rebuild the haplotypes without the 156 MB file
    subset = DATA / "results" / "HG002_chr21.vcf"
    with open(subset, "w") as fh:
        fh.write(
            "##fileformat=VCFv4.2\n##source=GIAB HG002 v4.2.1 benchmark, chr21 PASS subset via GenomeOS\n"
        )
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n")
        for v in variants:
            gt = ("|" if v.phased else "/").join(str(x) for x in v.gt) if v.gt else "./."
            fh.write(f"{v.chrom}\t{v.pos + 1}\t{v.id}\t{v.ref}\t{','.join(v.alts)}\t.\tPASS\t.\tGT\t{gt}\n")
    out["subset_file"] = str(subset)
    return out


def _distil_gencode_subset() -> dict:
    import gzip

    src = DATA / "reference" / "gencode.v50.annotation.gff3.gz"
    dst = DATA / "results" / "gencode_v50_chr21_chrM.gff3.gz"
    n = 0
    with gzip.open(src, "rt") as fi, gzip.open(dst, "wt") as fo:
        for line in fi:
            if line.startswith("#") or line.split("\t", 1)[0] in ("chr21", "chrM"):
                fo.write(line)
                n += 1
    return {"subset_file": str(dst), "lines": n, "chromosomes": ["chr21", "chrM"]}


def _distil_ccres_chr21() -> dict:
    from genomeos.genome.regulatory import save_ccres, stream_ccres, summarise

    elements = stream_ccres({"chr21", "chrM"})
    by_chrom: dict[str, list] = {}
    for c in elements:
        by_chrom.setdefault(c.chrom, []).append(c)
    out: dict = {
        "rows_streamed_from": "https://downloads.wenglab.org/V3/GRCh38-cCREs.bed",
        "disk_used_bytes": 0,
    }
    for chrom, els in by_chrom.items():
        p = save_ccres(chrom, els)
        out[chrom] = {**summarise(els), "file": str(p)}
    return out


def _distil_celegans_lineage() -> dict:
    from genomeos.organism.reference import KNOWLEDGE, ReferenceLineage, distil, fetch_wormweb, parse_wormweb

    data = distil(parse_wormweb(fetch_wormweb()))
    KNOWLEDGE.parent.mkdir(parents=True, exist_ok=True)
    KNOWLEDGE.write_text(json.dumps(data, separators=(",", ":")))
    ref = ReferenceLineage.from_dict(data)
    org = DATA / "organisms" / "celegans"
    org.mkdir(parents=True, exist_ok=True)
    (org / "timers.bio").write_text(ref.to_bio_timers())
    (org / "lineage_embryo.bio").write_text(ref.to_bio_program(0.0, 800.0))
    (org / "lineage_larva.bio").write_text(ref.to_bio_program(800.0))
    return {
        **{k: v for k, v in data.items() if k != "cells"},
        **ref.summary(),
        "knowledge_file": str(KNOWLEDGE),
        "generated": [str(org / f) for f in ("timers.bio", "lineage_embryo.bio", "lineage_larva.bio")],
        "cycle_timers": len(ref.cycle_stats()),
    }


def _distil_human_turnover() -> dict:
    from genomeos.organism.human import distil, fetch_milo, read_xlsx, to_bio_tissues

    table = distil(read_xlsx(fetch_milo()))
    org = DATA / "organisms" / "human"
    org.mkdir(parents=True, exist_ok=True)
    from genomeos.organism.haematopoiesis import mutants_bio
    from genomeos.organism.haematopoiesis import to_bio as haematopoiesis_bio

    (org / "tissues.bio").write_text(to_bio_tissues(table))
    (org / "haematopoiesis.bio").write_text(haematopoiesis_bio())
    (org / "haematopoiesis_mutants.bio").write_text(mutants_bio())
    table["generated"] = [
        str(org / f) for f in ("tissues.bio", "haematopoiesis.bio", "haematopoiesis_mutants.bio")
    ]
    return table


def _distil_celegans_tf_atlas() -> dict:
    from genomeos.organism.reference import ReferenceLineage
    from genomeos.organism.tf_atlas import distil, fetch, save_cells, summary, textbook_check, to_bio_reader

    table = distil(fetch())
    save_cells(table)
    ref = ReferenceLineage.load()
    checks = textbook_check(table, ref, load_result("celegans_packer2019"))
    org = DATA / "organisms" / "celegans"
    (org / "reader.bio").write_text(to_bio_reader(table))
    out = summary(table, checks)
    out["generated"] = [str(org / "reader.bio")]
    return out


def _distil_celegans_digital_development() -> dict:
    from genomeos.organism.digital_development import distil

    return distil()


def _distil_celegans_packer() -> dict:
    from genomeos.organism.packer import distil, stream_annotation
    from genomeos.organism.reference import ReferenceLineage

    return distil(stream_annotation(), ReferenceLineage.load())


DISTILLERS: list[Distiller] = [
    Distiller(
        "celegans_lineage",
        [],
        _distil_celegans_lineage,
        "the complete timed C. elegans lineage (WormWeb, CC BY) as a 2,183-cell table plus generated BioLang",
    ),
    Distiller(
        "human_cell_turnover",
        [],
        _distil_human_turnover,
        "Sender & Milo 2021 cell counts, lifespans and turnover per cell type (Summary.xlsx streamed)",
    ),
    Distiller(
        "celegans_tf_atlas",
        [],
        _distil_celegans_tf_atlas,
        "Ma 2021 TF protein atlas (88 MB streamed): factor presence per cell, reader.bio, textbook check",
    ),
    Distiller(
        "celegans_digital_development",
        [],
        _distil_celegans_digital_development,
        "Du 2014 founder fate changes per knockout (17 KB streamed) scored against the program's knockouts",
    ),
    Distiller(
        "celegans_packer2019",
        [],
        _distil_celegans_packer,
        "Packer 2019 lineage -> cell type table (GEO GSE126954, 2.8 MB streamed) checked against the lineage",
    ),
    Distiller(
        "encode_ccres_chr21",
        [],
        _distil_ccres_chr21,
        "ENCODE cCREs streamed (64 MB) and kept only for chr21/chrM: promoters, enhancers, CTCF sites",
    ),
    Distiller(
        "gencode_chr21_chrM",
        [DATA / "reference" / "gencode.v50.annotation.gff3.gz"],
        _distil_gencode_subset,
        "chr21 + chrM rows of GENCODE 50; the 153 MB genome-wide file is then disposable",
    ),
    Distiller(
        "clock_GSE41169",
        [DATA / "knowledge" / "GSE41169_series_matrix.txt.gz"],
        _distil_clock_geo,
        "Horvath/Hannum ages for every GEO sample; the 193 MB matrix is then disposable",
    ),
    Distiller(
        "clinvar_chr21_agreement",
        [DATA / "reference" / "clinvar.vcf.gz"],
        _distil_clinvar,
        "consequence agreement table; the 193 MB genome-wide ClinVar VCF is then disposable",
    ),
    Distiller(
        "library_members",
        [DATA / "knowledge" / "Ensembl2Reactome.txt"],
        _distil_libraries,
        "data-driven membership of all 45 libraries; the 183 MB Reactome mapping is then disposable",
    ),
    Distiller(
        "hg002_chr21",
        [DATA / "reference" / "HG002_GRCh38_benchmark.vcf.gz", DATA / "twins" / "HG002_chr21.fa"],
        _distil_hg002_chr21,
        "haplotype statistics plus a chr21-only VCF; the 156 MB VCF and 94 MB FASTA are then disposable",
    ),
]


def distil(only: list[str] | None = None, force: bool = False) -> list[tuple[str, str]]:
    done = []
    for d in DISTILLERS:
        if only and d.name not in only:
            continue
        if load_result(d.name) and not force:
            done.append((d.name, "kept existing summary"))
            continue
        if not all(p.exists() for p in d.inputs):
            done.append((d.name, "inputs missing: " + ", ".join(str(p) for p in d.inputs if not p.exists())))
            continue
        payload = d.run()
        payload["inputs"] = [str(p) for p in d.inputs]
        payload["describe"] = d.describe
        save_result(d.name, payload)
        done.append((d.name, "distilled"))
    return done


def clean(dry_run: bool = True) -> list[tuple[str, int, str]]:
    """Delete raw inputs whose summary exists. Returns (path, bytes, action)."""
    out = []
    for d in DISTILLERS:
        if not load_result(d.name):
            continue
        for p in d.inputs:
            if p.exists():
                size = p.stat().st_size
                if not dry_run:
                    p.unlink()
                    fai = p.with_name(p.name + ".fai")
                    if fai.exists():
                        fai.unlink()
                out.append((str(p), size, "would delete" if dry_run else "deleted"))
    return out


def manifest() -> dict:
    return {
        d.name: {
            "inputs": [str(p) for p in d.inputs],
            "summary_exists": load_result(d.name) is not None,
            "describe": d.describe,
        }
        for d in DISTILLERS
    }


def to_json(obj) -> str:
    return json.dumps(obj, indent=2, default=str)
