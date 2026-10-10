#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Canonical against tissue isoform: how often a person's coding SNVs read differently, registered first.

    uv run python scripts/flow_tissue_census.py --register   # the registration, before any count
    uv run python scripts/flow_tissue_census.py              # the census

`genomeos flow` traced a variant on the canonical transcript, which a tissue may not make. The trace
now follows the tissue's dominant isoform (genomeos/flow/tissue.py). This census asks, over HG002's own
coding SNVs in the genes the project holds GTEx isoform tables for, how many change consequence class
or residue number between the canonical transcript and each tissue's dominant isoform, per tissue, with
denominators. It is a census of disagreement between two predictions, not an accuracy claim.

The registration names the gene set, the variant set, the tissue rule, the comparison, the outputs and
the wording, pins the code by sha256 and the isoform tables by digest, and holds no count. The census
refuses to start unless the registration is committed, the frozen code on disk is the code it pins and
the tables still hash to what it records. The committed result carries counts only: no position,
allele, genotype or protein change of the person (docs/STORAGE.md).

`--results-dir` and `--individuals-root` exist for a dry run on a synthetic person outside the result
registry; a registry run reads HG002 and nothing else.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from datetime import date
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path and not (Path.cwd() / "genomeos").is_dir():
    sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.flow import tissue as ft  # noqa: E402
from genomeos.genome import Annotation, IndexedGenome, default_gencode, reference_fasta  # noqa: E402
from genomeos.genome.individuals import ROOT as INDIVIDUALS  # noqa: E402
from genomeos.genome.individuals import rows_in, vcf_path  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

REG_NAME = "flow_tissue_census_registration"
RESULT_NAME = "flow_tissue_census_result"
PERSON = "HG002"
CHROMS = tuple(f"chr{i}" for i in range(1, 23))
TABLES = sorted(ft.ISOFORM_CACHE.glob("gtex_isoforms_*.json"))
OWN_CODE = (
    "genomeos/flow/__init__.py",
    "genomeos/flow/tissue.py",
    "genomeos/flow/trace.py",
    "scripts/flow_tissue_census.py",
    "tests/test_flow_tissue.py",
)
#: The code the census must be: pinned by sha256, so an edit after registration stops the run.
FROZEN = ("genomeos/flow/tissue.py", "genomeos/flow/trace.py", "scripts/flow_tissue_census.py")

SOURCES = [
    {
        "accession": "GTEx Portal API v2 medianTranscriptExpression, datasetId gtex_v8",
        "version": "GTEx v8 (RSEM transcript TPM on GENCODE v26 models), cached per gene by "
        "genomeos.molecules.rna.gtex_isoforms",
    },
    {"accession": "GENCODE comprehensive annotation (GFF3)", "version": "v50"},
    {"accession": "UCSC hg38 chromosome FASTA (GRCh38)", "version": "hg38"},
    {
        "accession": "GIAB v4.2.1 benchmark calls (open-consent HG002), PASS rows per chromosome",
        "version": "NISTv4.2.1 GRCh38, chr1-22",
    },
]

REGISTRATION: dict[str, Any] = {
    "question": "Over one person's own coding SNVs, how many change consequence class or residue number "
    "between the canonical transcript and a tissue's dominant isoform, per tissue, with denominators? "
    "A census of disagreement between two predictions of the local trace; not a test and not an "
    "accuracy claim.",
    "person": f"{PERSON}: GIAB v4.2.1 benchmark calls (PASS rows inside the benchmark regions, chr1-22) "
    "imported locally (genomeos.genome.individuals.vcf_path). Its genotypes never enter the committed "
    "result: the result carries counts only, with no position, allele, genotype or protein change.",
    "genes": "every cached GTEx isoform table (data/knowledge/expression/gtex_isoforms_*.json, pinned "
    "by digest in this file's manifest) that names a GTEx gencode_id, matched to the local GENCODE v50 "
    "models by the unversioned Ensembl gene id, on chr1-22 where the person has calls and the reference "
    "sequence is local, and whose canonical transcript codes a protein. Tables without a gencode_id "
    "('not in GTEx') and genes the models or the person's chromosomes lack are counted as exclusions. "
    "A gene with no universe variant counts in the gene tally and adds nothing to the per-tissue counts.",
    "gene_set_caveat": "the cached tables are the project's, not a sample of the genome: 96 of the 99 were "
    "fetched for the translation-disagreement genes (genes whose reference canonical transcript does "
    "not translate to the UniProt protein; data/results/translation_disagreements_isoforms.json) and "
    "the other three are APP, FUT2 and TP53. The census describes these genes only and is not "
    "extrapolated. A genome-wide census needs the GTEx v8 transcript TPM table, a download this lane "
    "was not cleared to make.",
    "variants": "the person's VCF rows inside the gene's locus with a one-base REF, one one-base ALT and "
    "a genotype carrying a non-reference allele, read on the canonical transcript and on every tissue's "
    "established dominant isoform; the gene's universe is those in the CDS of at least one of them. "
    "Indels, multi-allelic rows and rows whose genotype carries no ALT are counted as exclusions.",
    "canonical": "the transcript trace_gene follows (genomeos/flow/trace.py canonical_transcript): the "
    "Ensembl_canonical-tagged coding transcript in the local models, else the longest CDS.",
    "tissues": "every tissue named in the included tables (GTEx tissueSiteDetailId), each read for every "
    "census gene; a tissue a gene's table does not measure is UNKNOWN tissue_not_measured for it.",
    "rule": {
        "min_tpm": ft.MIN_TPM,
        "tie_ratio": ft.TIE_RATIO,
        "established": "the gene has an isoform table; GTEx measured the tissue; the top transcript's "
        "median TPM in the tissue is at least min_tpm; the next transcript's is below tie_ratio times "
        "the top's; the top transcript is in the local models and codes a protein there",
        "reasons": list(ft.REASONS),
        "never": "an UNKNOWN tissue is never read on the canonical transcript in its place",
    },
    "comparison": {
        "same": "same consequence class, residue number, HGVS p. and c. on both, or the tissue's dominant "
        "isoform IS the canonical transcript",
        "class_changed": "a different consequence class (e.g. missense on the canonical, intron on the "
        "tissue's isoform; synonymous on one, missense on the other)",
        "residue_changed": "the same class at a different residue number",
        "change_differs": "the same class and residue number, another amino acid or codon (counted "
        "beside the disagreement, not in it)",
        "disagree": "class_changed + residue_changed",
    },
    "outputs": "per tissue: variants (the denominator, the same for every tissue), unknown by reason, "
    "established, same_transcript, other_transcript, same, class_changed, residue_changed, "
    "change_differs, disagree, genes_established_other_transcript, genes_disagreeing; over all "
    "variant-tissue pairs the same sums (pairs = variants x tissues); the tissues ranked by disagree; "
    "per gene counts only (variants, coding_on_canonical, tissues where another isoform is established, "
    "tissues with a disagreement); the exclusions with their counts",
    "prediction": "none: a census, registered so that the gene set, the variant set, the rule and the "
    "wording are fixed before any count is seen",
    "wording": {
        "headline": "Over {variants} coding SNVs of {person} in {genes} genes with cached GTEx isoform "
        "tables, the consequence predicted on a tissue's dominant isoform differs in class or residue "
        "number from the one predicted on the canonical transcript in {disagree} of {established} "
        "variant-tissue pairs where the tissue's isoform is established ({unknown} of {pairs} pairs "
        "UNKNOWN; per tissue in per_tissue).",
        "reading": "Both consequences are predictions of the local trace on GTEx's population-median "
        "isoform, not measurements and not this person's own expression; the genes are the project's "
        "cached set, not a sample of the genome.",
        "if_no_variant": "No coding SNV of {person} falls in the census genes: the census is empty, and "
        "says nothing about disagreement.",
    },
    "what_this_cannot_say": [
        "which consequence is right: both are predictions; neither is checked against protein data",
        "whether this person makes the isoform: the call is GTEx's median across donors",
        "anything about genes outside the cached tables, or about indels",
        "whether GENCODE v26's model of a transcript (GTEx's quantification) equals the v50 model traced",
    ],
    "what_is_already_known": [
        "the area C reading that the cerebellum's dominant APP isoform is APP695 while the canonical is "
        "APP770 (genomeos individual protein, ROADMAP D)",
        "data/results/translation_disagreements_isoforms.json, published earlier, names each of its 96 "
        "genes' most often dominant transcript and whether it codes in the local models (per gene, not "
        "per tissue, and with no person in it)",
        "before this registration the census code ran only on synthetic transcripts (tests) and on a "
        "synthetic person outside the registry (one SNV at every eleventh base of each census gene, the "
        "reference base read, the next base by A>G, C>T, G>A, T>C as ALT, all heterozygous) to check that "
        "it runs and which files it reads; it has not been run on HG002 or on any real person",
        "what that dry run printed, which reflects the GTEx side and the gene set, not any person: 99 "
        "tables, 86 with a gencode_id, 83 census genes traced; 11,666 synthetic coding SNVs (11,653 "
        "coding on the canonical); over 629,964 variant-tissue pairs (54 tissues), 172,339 established "
        "(122,106 on the canonical itself, 50,233 on another isoform) and 457,625 UNKNOWN (low_expression "
        "420,547, dominant_non_coding 19,650, tie 14,617, dominant_not_in_local_models 2,811); "
        "class_changed 10,104, residue_changed 20,972, change_differs 0. A synthetic SNV every eleventh "
        "base weights genes by length, so these figures are not an expectation for a person's variants",
    ],
}


def sha256(path: str | Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def committed(path: Path) -> bool:
    tracked = subprocess.run(["git", "ls-files", "--error-unmatch", str(path)], capture_output=True)
    clean = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", str(path)])
    return tracked.returncode == 0 and clean.returncode == 0


def is_registry(results_dir: Path) -> bool:
    return Path(results_dir).resolve() == Path(RESULTS_DIR).resolve()


def tables_entry() -> dict[str, Any]:
    return mf.files_entry("data/knowledge/expression/gtex_isoforms_*.json", TABLES)


def manifest(inputs: list[dict[str, Any]], parameters: dict[str, Any], exclusions: list[Any]) -> dict:
    return {
        "sources": SOURCES,
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": parameters,
        "exclusions": exclusions,
        "partitions": "n/a: a census of one person's variants; nothing is fitted, so nothing is held out",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def register(results_dir: Path) -> None:
    payload = {
        "result": REG_NAME,
        "date": date.today().isoformat(),
        "registration": REGISTRATION,
        "tables": len(TABLES),
        "frozen": {p: sha256(p) for p in FROZEN},
    }
    params = {"mode": "register", "min_tpm": ft.MIN_TPM, "tie_ratio": ft.TIE_RATIO, "person": PERSON}
    print("saved", save_result(REG_NAME, payload, results_dir, manifest([tables_entry()], params, [])))


def included_tables() -> tuple[dict[str, dict], dict[str, int]]:
    """Unversioned gene id → its isoform table, and the tables left out with why."""
    keep: dict[str, dict] = {}
    out = {"table_without_gencode_id": 0, "table_unreadable": 0}
    for p in TABLES:
        try:
            d = json.loads(p.read_text())
        except (OSError, json.JSONDecodeError):
            out["table_unreadable"] += 1
            continue
        if not d.get("gencode_id") or not d.get("isoforms"):
            out["table_without_gencode_id"] += 1
            continue
        keep[d["gencode_id"].split(".")[0]] = d
    return keep, out


def person_snvs(vcf: Path, start1: int, end1: int, ex: dict[str, int]) -> list[tuple[int, str, str]]:
    """The person's SNVs inside a 1-based closed span, as (0-based position, ref, alt)."""
    out = []
    for f in rows_in(vcf, start1, end1):
        ref, alt = f[3].upper(), f[4].upper()
        gt = f[9].split(":")[0] if len(f) > 9 else ""
        if "," in alt:
            ex["multi_allelic_rows"] += 1
        elif len(ref) != 1 or len(alt) != 1:
            ex["indel_or_mnv_rows"] += 1
        elif not any(a not in ("0", ".") for a in gt.replace("|", "/").split("/")):
            ex["rows_without_alt_genotype"] += 1
        else:
            out.append((int(f[1]) - 1, ref, alt))
    return out


def census(name: str, root: Path) -> tuple[dict[str, Any], dict[str, list[Path]]]:
    """The census over the included genes; and every data file it read, by kind."""
    keep, excluded = included_tables()
    tissues = sorted({t for d in keep.values() for t in ft.measured_tissues(d)})
    ex = {
        **excluded,
        "gene_not_in_local_models_on_chr1_22": 0,
        "gene_on_a_chromosome_without_calls_or_sequence": 0,
        "gene_without_coding_canonical": 0,
        "multi_allelic_rows": 0,
        "indel_or_mnv_rows": 0,
        "rows_without_alt_genotype": 0,
    }
    found: set[str] = set()
    genes: list[dict[str, Any]] = []
    snvs_in_loci = 0
    read: dict[str, list[Path]] = {"gff": [], "fasta": [], "vcf": []}
    for chrom in CHROMS:
        gff = default_gencode({chrom})
        if not gff:
            continue
        read["gff"].append(Path(gff))
        ann = Annotation.from_gff3(gff, {chrom})
        here = [g for g in ann.genes.values() if g.id.split(".")[0] in keep and g.locus.chrom == chrom]
        if not here:
            continue
        found.update(g.id.split(".")[0] for g in here)
        vcf = vcf_path(name, chrom, root)
        fa = reference_fasta(chrom)
        if vcf is None or not fa.exists():
            ex["gene_on_a_chromosome_without_calls_or_sequence"] += len(here)
            continue
        read["vcf"].append(Path(vcf))
        read["fasta"].append(Path(fa))
        module = ann.to_module("tissue_census")
        genome = IndexedGenome(str(fa))
        try:
            for g in sorted(here, key=lambda x: x.locus.start):
                snvs = person_snvs(vcf, g.locus.start + 1, g.locus.end, ex)
                snvs_in_loci += len(snvs)
                table = keep[g.id.split(".")[0]]
                r = ft.gene_census(genome, g, module.entities[g.id].transcripts, table, snvs, tissues)
                if r is None:
                    ex["gene_without_coding_canonical"] += 1
                    continue
                r["chrom"] = chrom
                genes.append(r)
        finally:
            genome.close()
    ex["gene_not_in_local_models_on_chr1_22"] = len(set(keep) - found)
    with_variants = [g for g in genes if g["variants"]]
    per_tissue = ft.combine(with_variants, tissues)
    pair_keys = ("variants", "established", "same_transcript", "other_transcript", *ft.AGREEMENTS, "disagree")
    totals: dict[str, Any] = {k: sum(c[k] for c in per_tissue.values()) for k in pair_keys}
    totals["pairs"] = totals.pop("variants")
    totals["unknown"] = {r: sum(c["unknown"][r] for c in per_tissue.values()) for r in ft.REASONS}
    variants = sum(g["variants"] for g in with_variants)
    w = REGISTRATION["wording"]
    headline = (
        w["if_no_variant"].format(person=name)
        if not variants
        else w["headline"].format(
            person=name,
            variants=variants,
            genes=len(with_variants),
            disagree=totals["disagree"],
            established=totals["established"],
            unknown=sum(totals["unknown"].values()),
            pairs=totals["pairs"],
        )
    )
    out = {
        "result": RESULT_NAME,
        "date": date.today().isoformat(),
        "person": name,
        "headline": headline,
        "reading": w["reading"],
        "genes": {
            "tables": len(TABLES),
            "tables_included": len(keep),
            "traced": len(genes),
            "with_universe_variants": len(with_variants),
        },
        "variants": {
            "snvs_in_census_gene_loci": snvs_in_loci,
            "universe": variants,
            "coding_on_canonical": sum(g["coding_on_canonical"] for g in with_variants),
            "coding_only_on_a_tissue_isoform": sum(
                g["variants"] - g["coding_on_canonical"] for g in with_variants
            ),
        },
        "tissues": len(tissues),
        "totals": totals,
        "tissues_by_disagreement": [
            [t, c["disagree"], c["established"], c["variants"]]
            for t, c in sorted(per_tissue.items(), key=lambda kv: (-kv[1]["disagree"], kv[0]))
        ],
        "per_tissue": per_tissue,
        "per_gene": [
            {
                "gene": g["gene"],
                "chrom": g["chrom"],
                "variants": g["variants"],
                "coding_on_canonical": g["coding_on_canonical"],
                "tissues_other_isoform_established": sum(
                    1 for c in g["per_tissue"].values() if c["call"]["dominant_is_canonical"] is False
                ),
                "tissues_disagreeing": sum(1 for c in g["per_tissue"].values() if c["disagree"]),
            }
            for g in with_variants
        ],
        "exclusions": ex,
    }
    return out, read


def _sequence_files(fa: Path) -> tuple[Path, ...]:
    """The sequence file read, its indexes, and the other form beside it, which `resolve_fasta`
    opens to ask whether it is blocked."""
    other = fa.with_suffix("") if fa.suffix == ".gz" else fa.with_name(fa.name + ".gz")
    return (fa, fa.with_name(fa.name + ".fai"), fa.with_name(fa.name + ".gzi"), other)


def run(results_dir: Path, name: str, root: Path) -> int:
    reg_path = Path(results_dir) / f"{REG_NAME}.json"
    if not reg_path.exists():
        raise SystemExit(f"{reg_path} does not exist: register first")
    registry = is_registry(results_dir)
    if registry and not committed(reg_path):
        raise SystemExit(f"{reg_path} is not committed: register, commit, then run")
    if registry and name != PERSON:
        raise SystemExit(f"a registry run reads {PERSON} only")
    reg = json.loads(reg_path.read_text())
    drift = [p for p, h in reg["frozen"].items() if sha256(p) != h]
    if drift:
        raise SystemExit(f"frozen code changed since registration: {drift}")
    pinned = (reg.get("result_manifest") or {}).get("inputs") or []
    if not pinned or pinned[0].get("sha256") != tables_entry()["sha256"]:
        raise SystemExit("the GTEx isoform tables are not the ones the registration pins")
    t0 = time.time()
    out, read = census(name, root)
    out["registration"] = {"file": str(reg_path), "date": reg.get("date"), "frozen": reg["frozen"]}
    out["seconds"] = round(time.time() - t0)
    inputs = [mf.input_entry(reg_path), tables_entry()]
    if read["gff"]:
        inputs.append(mf.files_entry("GENCODE v50 annotation, per chromosome", read["gff"]))
    if read["vcf"]:
        inputs.append(mf.files_entry(f"{name} calls, per chromosome", read["vcf"]))
    fasta = [p for f in read["fasta"] for p in _sequence_files(f) if p.exists()]
    if fasta:
        inputs.append(mf.files_entry("hg38 sequence and its indexes, per chromosome", fasta))
    undeclared = mf.undeclared_reads({"inputs": inputs}, mf.traced_reads())
    if undeclared:
        print("undeclared reads:", undeclared, flush=True)
    params = {"min_tpm": ft.MIN_TPM, "tie_ratio": ft.TIE_RATIO, "person": name, "chromosomes": list(CHROMS)}
    exclusions = [f"{k}: {v}" for k, v in out["exclusions"].items()]
    print("saved", save_result(RESULT_NAME, out, results_dir, manifest(inputs, params, exclusions)))
    print(out["headline"])
    print(json.dumps({k: out[k] for k in ("genes", "variants", "totals")}, indent=1))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--register", action="store_true", help="write the registration")
    ap.add_argument("--results-dir", type=Path, default=RESULTS_DIR, help="a scratch directory for a dry run")
    ap.add_argument("--name", default=PERSON, help="the person (a dry run only; a registry run is HG002)")
    ap.add_argument("--individuals-root", type=Path, default=INDIVIDUALS)
    a = ap.parse_args()
    if a.register:
        register(a.results_dir)
        return 0
    return run(a.results_dir, a.name, a.individuals_root)


if __name__ == "__main__":
    sys.exit(main())
