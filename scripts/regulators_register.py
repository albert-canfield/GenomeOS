# SPDX-License-Identifier: AGPL-3.0-or-later
"""M1b registration: one network-level prediction of the compiled regulator rules, fixed before its
endpoint is read.

    uv run python scripts/regulators_register.py [--no-save]

Run after the census (data/results/regulators_census.json, its assignment table) and the inventory
(the ENCODE knockdown metadata), and after the analysis code is committed: it refuses when
genomeos/attribution/regulators.py or scripts/regulators_run.py differ from HEAD, and pins both by
sha256. It reads no expression value: every knockdown file is still undownloaded when it runs. It
writes the git-ignored test inputs (the arms, the universe candidates, the peaks) to
data/cache/regulators/test_inputs.json.gz and pins them by sha256 in the registration.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import holdout as ho  # noqa: E402
from genomeos.attribution import regulators as rg  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.lang.parser import parse  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

RESULT = "regulators_registration"
FROZEN = ("genomeos/attribution/regulators.py", "scripts/regulators_run.py")
TEST_INPUTS = rg.CACHE / "test_inputs.json.gz"
META = rg.KD_DIR / "experiments.json"
ASSIGNMENT = rg.CACHE / "assignment.tsv.gz"
OWN_CODE = (
    "genomeos/attribution/regulators.py",
    "scripts/regulators_inventory.py",
    "scripts/regulators_census.py",
    "scripts/regulators_register.py",
    "scripts/regulators_fetch.py",
    "scripts/regulators_run.py",
    "tests/test_regulators.py",
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def frozen_code() -> dict:
    dirty = subprocess.run(
        ["git", "status", "--porcelain", "--", *FROZEN], capture_output=True, text=True
    ).stdout
    if dirty.strip():
        raise SystemExit(f"refused: the analysis code differs from HEAD; commit it first:\n{dirty}")
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    last = {
        p: subprocess.run(
            ["git", "log", "-1", "--format=%H", "--", p], capture_output=True, text=True
        ).stdout.strip()
        for p in FROZEN
    }
    return {
        "head_at_registration": head,
        "last_commit_per_file": last,
        "sha256": {p: sha256(Path(p)) for p in FROZEN},
    }


def gencode_inputs() -> list[Path]:
    """The files `holdout.genes()` reads: its cache when present, the GENCODE files otherwise."""
    paths = ho.gencode_paths()
    cached = ho._cache_path("genes", ho._files_key(paths))
    return [cached] if cached.exists() else paths


def gene_index() -> dict[str, tuple[str, str, int]]:
    """Symbol -> (Ensembl id, chromosome, TSS) from GENCODE v50, a protein-coding gene preferred."""
    out: dict[str, tuple[int, str, str, int]] = {}
    for g in ho.genes():
        if not g.symbol:
            continue
        rank = 0 if g.gene_type == "protein_coding" else 1
        if g.symbol not in out or rank < out[g.symbol][0]:
            out[g.symbol] = (rank, g.gene_id, g.chrom, g.tss)
    return {s: (v[1], v[2], v[3]) for s, v in out.items()}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args(argv)
    t0 = time.time()
    code = frozen_code()
    census = json.loads((RESULTS_DIR / "regulators_census.json").read_text())
    meta = json.loads(META.read_text())
    exps = meta["experiments"]
    cl = rg.load_clusters()
    chip = cl.factors_in()
    genes = gene_index()

    # the units: TF knockdown experiments whose target ENCODE also profiled by ChIP in the same cell
    units: dict[tuple[str, str], list[dict]] = defaultdict(list)
    not_admitted = Counter()
    for row in meta["listed"]:
        e = exps[row["accession"]]
        cell = rg.context_label(row["cell"])
        if e["target"] not in chip.get(cell, set()):
            not_admitted["no_chip_of_the_target_in_this_cell"] += 1
            continue
        if e["target"] not in genes:
            not_admitted["target_not_in_gencode_v50"] += 1
            continue
        files = rg.choose_files(e, [exps[c] for c in e["possible_controls"] if c in exps])
        if files is None:
            not_admitted["no_control_files_on_a_shared_annotation"] += 1
            continue
        units[(e["target"], cell)].append(
            {
                "accession": e["accession"],
                "assay": row["assay"],
                "lab": e["lab"],
                "controls": e["possible_controls"],
                "annotation": files["annotation"],
                "knockdown_files": [
                    {k: f[k] for k in ("accession", "md5sum", "file_size", "href", "biological_replicates")}
                    for f in files["knockdown"]
                ],
                "control_files": [
                    {k: f[k] for k in ("accession", "md5sum", "file_size", "href", "biological_replicates")}
                    for f in files["control"]
                ],
            }
        )

    # the compiled programme of each test cell, and each unit factor's rule arm in it
    rule_genes: dict[str, set[str]] = defaultdict(set)
    unmapped: Counter = Counter()
    for chrom in rg.CHROMS:
        path = Path(COMPILED_DIR) / f"noncoding_{chrom}.bio"
        for r in rg.compiled_rules(parse(path.read_text(), path.stem), path.stem):
            if r.cell in rg.TEST_CELLS:
                if r.gene in genes:
                    rule_genes[r.cell].add(genes[r.gene][0])
                else:
                    unmapped[r.cell] += 1
    wanted = {(f, c) for f, c in units}
    arm: dict[tuple[str, str], set[str]] = defaultdict(set)
    arm_motif: dict[tuple[str, str], set[str]] = defaultdict(set)
    directed: dict[tuple[str, str], dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    with gzip.open(ASSIGNMENT, "rt") as fh:
        header = next(fh).rstrip("\n").split("\t")
        col = {k: i for i, k in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            key = (f[col["factor"]], f[col["cell"]])
            if key not in wanted or f[col["gene"]] not in genes:
                continue
            gid = genes[f[col["gene"]]][0]
            arm[key].add(gid)
            if f[col["kind"]] == "occupancy_same_cell_with_motif":
                arm_motif[key].add(gid)
            if f[col["direction"]] != "unknown":
                directed[key][f[col["direction"]]].add(gid)
    tss = {gid: (ch, pos) for _s, (gid, ch, pos) in genes.items()}
    needed = set().union(*rule_genes.values()) | {genes[f][0] for f, _c in units}
    peaks = {f"{f} {c}": rg.factor_peaks(cl, f, c) for f, c in sorted(units)}
    inputs = {
        "rule_genes": {c: sorted(v) for c, v in sorted(rule_genes.items())},
        "arms": {f"{f} {c}": sorted(arm[(f, c)]) for f, c in sorted(units)},
        "arms_with_motif": {f"{f} {c}": sorted(arm_motif[(f, c)]) for f, c in sorted(units)},
        "directed": {
            f"{f} {c}": {d: sorted(g) for d, g in directed[(f, c)].items()} for f, c in sorted(units)
        },
        "factor_gene": {f: genes[f][0] for f, _c in sorted(units)},
        "tss": {g: tss[g] for g in sorted(needed)},
        "peaks": peaks,
    }
    TEST_INPUTS.parent.mkdir(parents=True, exist_ok=True)
    TEST_INPUTS.write_bytes(gzip.compress(json.dumps(inputs, sort_keys=True).encode(), mtime=0))
    unit_list = [
        {
            "factor": f,
            "cell": c,
            "experiments": sorted(units[(f, c)], key=lambda x: x["accession"]),
            "rule_arm_genes_before_expression": len(arm[(f, c)]),
            "rule_arm_with_motif_before_expression": len(arm_motif[(f, c)]),
            "peaks_in_cell": sum(len(v) for v in peaks[f"{f} {c}"].values()),
        }
        for f, c in sorted(units)
    ]
    files_total = {
        x["accession"]: x["file_size"]
        for u in unit_list
        for e in u["experiments"]
        for x in e["knockdown_files"] + e["control_files"]
    }
    reg = {
        "lane": rg.LANE,
        "roadmap_item": "section 5 item 7b, M1b",
        "question": (
            "when ENCODE knocked down a transcription factor in K562 or HepG2, did the genes whose compiled "
            "rules in that cell pass through an element the factor occupies there change more than the "
            "other genes of the same compiled programme, more than genes near the factor's peaks, and more "
            "than another factor's assigned genes?"
        ),
        "why_network_level": (
            "the prediction is about genes the compiled programme connects to a regulator only through "
            "its elements: a regulator -> element -> gene path, which no compiled rule held before M1b "
            "(0 of 440,589 rules named a gene as its source). A knockdown changes the regulator, not the "
            "element, so the endpoint reads the path, not any one rule's input observation"
        ),
        "census_facts": {k: census["totals"][k] for k in ("rules", "rules_gaining_a_source", "by_kind")},
        "endpoint_exposure": (
            "none before this registration: the ENCODE knockdown and control quantification files are "
            "undownloaded; the lane read only portal metadata (experiment identities, targets, controls, "
            "file accessions, md5, sizes). No lane of this project has read these experiments (searched "
            "2026-10-10). The compiled rules come from AlphaGenome, trained on ENCODE tracks that may "
            "include K562 and HepG2 RNA-seq, but not on these knockdowns"
        ),
        "population": (
            "factor-cell units: an ENCODE experiment of an assay in kd_assays, in a test cell, whose target "
            "ENCODE classifies as a transcription factor and also profiled by ChIP-seq in the same cell "
            "(the UCSC clusters), with a control experiment sharing a GRCh38 annotation. An experiment is "
            "used when (1) the factor's own gene has >= 1 mean TPM in control and <= 0.70 of it in "
            "knockdown (>= 30% lower; the threshold Replogle et al. 2022 use for a strong perturbation, "
            "adopted here, not ENCODE's), and (2) its universe holds >= 10 rule-arm genes and >= 10 others. "
            "A unit enters with >= 1 used experiment and its AUROCs are the means over them. A cell type "
            "with a single unit is set aside (it cannot be deranged within its cell)"
        ),
        "universe": (
            "per experiment: Ensembl genes with >= 1 compiled rule in the unit's cell (symbols mapped by "
            "GENCODE v50, protein-coding preferred) and >= 1 mean TPM in the control replicates, minus the "
            "factor's own gene and genes whose TSS lies within 10 kb of the factor's TSS"
        ),
        "measured": (
            "response = |log2((mean TPM knockdown + 1) / (mean TPM control + 1))| per gene, the replicates "
            "of the knockdown experiment against all replicates of its listed control experiments"
        ),
        "rule_arm": (
            "universe genes with >= 1 compiled rule in the unit's cell whose element a ChIP-seq peak "
            "cluster of the factor overlaps in that cell (census kinds occupancy_same_cell and "
            "occupancy_same_cell_with_motif)"
        ),
        "baseline": (
            "genes near the factor's peaks, ignoring every compiled rule: the same number of universe genes "
            "as the rule arm, ranked by the factor's peak midpoints in that cell within 50 kb of the TSS, "
            "then by distance to the nearest midpoint, then by gene id"
        ),
        "metric": (
            "per experiment, the stratified AUROC of each arm against the rest of the universe on the "
            "response, pairs taken only within deciles of log2(control TPM + 1); per unit, the mean over "
            "its experiments; the primary is the mean over units of (rule-arm AUROC - baseline AUROC), with "
            "a 95% percentile bootstrap interval over units (10,000 resamples, seed 20261010)"
        ),
        "matched_control": (
            "factor identity: each unit relabelled with another unit-factor's rule arm from the same cell "
            "(a uniform derangement within each cell type, 1,000 draws, seed 20261011), scored inside the "
            "unit's own universe and response; p = (1 + draws whose mean AUROC >= the observed mean "
            "rule-arm AUROC) / 1,001"
        ),
        "falsifier": (
            "control p >= 0.05: the compiled assignment does not pick out a factor's responding genes better "
            "than another factor's assignment does, so the regulator -> element -> gene path carries no "
            "factor-specific prediction on this endpoint, whatever the comparison with the baseline"
        ),
        "decision": (
            "in order: fewer than 20 units -> stopped_below_floor; control p >= 0.05 -> not_factor_specific; "
            "else the interval of the mean paired difference: lower bound > 0 -> rules_beat_peaks, upper "
            "bound < 0 -> peaks_beat_rules, otherwise no_difference_detected"
        ),
        "stop_rule": (
            "the run happens once (the runner refuses a second); fewer than 20 units stops it before any "
            "comparison with no rule relaxed; a knockdown file whose md5 is not ENCODE's stops it before "
            "any value is read; no threshold, universe or arm changes after a value is read"
        ),
        "outcomes": rg.OUTCOMES,
        "descriptives_not_verdicts": [
            "per experiment: knockdown ratio, universe, arm and baseline sizes, assay",
            "the rule arm restricted to elements where the factor's own JASPAR profile also matches",
            "K562 units alone, the paired difference and its interval without the control",
            "for directed candidate rows (curated GO direction agreeing with the element's rule): the "
            "share of genes whose signed change on knockdown has the sign the direction predicts "
            "(activates -> lower, inhibits -> higher), beside the same share among the universe's other "
            "genes",
        ],
        "limitations": [
            "occupancy is not function: a factor bound at an element need not regulate the gene, and a "
            "knockdown also acts through indirect paths, so a response is not a validated connection",
            "the knockdown lowers the factor's RNA; protein loss and timing are undocumented per experiment",
            "two to four replicates per arm: the response is a fold change, not a test statistic",
            "ENCODE profiled factors by ChIP in few cell types, so the test covers K562 and HepG2 only",
            "the clusters are ENCODE 3's uniform peaks (2020), not the ENCODE4 set",
        ],
        "constants": rg.TEST_CONSTANTS,
        "not_admitted_experiments": dict(not_admitted),
        "compiled_genes_per_cell": {c: len(v) for c, v in sorted(rule_genes.items())},
        "compiled_rules_with_unmapped_gene_per_cell": dict(unmapped),
        "units_before_expression": len(unit_list),
        "units": unit_list,
        "files_to_fetch": len(files_total),
        "bytes_to_fetch": sum(v or 0 for v in files_total.values()),
        "test_inputs": {"path": str(TEST_INPUTS), "sha256": sha256(TEST_INPUTS), "committed": False},
        "frozen_code": code,
        "what_the_run_will_write": "data/results/regulators_result.json",
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
        "seconds": round(time.time() - t0, 1),
    }
    man = {
        "sources": [
            {"accession": rg.CLUSTER_SOURCE["url"], "version": rg.CLUSTER_SOURCE["last_modified"]},
            {
                "accession": "https://www.encodeproject.org JSON API",
                "version": "metadata fetched " + meta["fetched"],
            },
            {"accession": "data/results/regulators_census.json", "version": "the committed census"},
        ],
        "inputs": [
            mf.input_entry(rg.CLUSTERS, partition=None),
            mf.input_entry(META, partition=None),
            mf.input_entry(ASSIGNMENT, partition=None),
            mf.input_entry(TEST_INPUTS, partition=None),
            mf.input_entry(RESULTS_DIR / "regulators_census.json", partition=None),
            mf.files_entry(
                "the compiled programs", [Path(COMPILED_DIR) / f"noncoding_{c}.bio" for c in rg.CHROMS]
            ),
            mf.files_entry("GENCODE v50 genes, as attribution.holdout reads them", gencode_inputs()),
        ],
        "input_count": 7,
        "assembly": "GRCh38; GENCODE v50 for symbols and TSS; the knockdown files' own annotation per unit",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": rg.TEST_CONSTANTS,
        "exclusions": [{"reason": k, "count": v} for k, v in sorted(not_admitted.items())],
        "partitions": "n/a: the endpoint is outside every training set this project fitted",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }
    if not args.no_save:
        print("saved", save_result(RESULT, reg, manifest=man))
    print(
        json.dumps(
            {
                k: reg[k]
                for k in (
                    "not_admitted_experiments",
                    "compiled_genes_per_cell",
                    "units_before_expression",
                    "files_to_fetch",
                    "bytes_to_fetch",
                )
            },
            indent=1,
        )
    )
    for u in unit_list:
        print(
            u["factor"],
            u["cell"],
            len(u["experiments"]),
            u["rule_arm_genes_before_expression"],
            u["rule_arm_with_motif_before_expression"],
            u["peaks_in_cell"],
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
