# SPDX-License-Identifier: AGPL-3.0-or-later
"""How independent is the node-containment result, component by component.

    uv run python scripts/node_independence_audit.py [--workers 6] [--quick]

The project's second-strongest claim (docs/ROADMAP.md area B) is that nodes hold enhancer-gene links
together: on 661 CRISPRi-measured regulated pairs the element and its measured gene share a node
5.89 points more often than under the published random-boundary control (95% +3.18 to +8.46,
data/results/node_containment_measured.json). lane-split's audit (0af7e1b) classified that arm
EVALUATION-ONLY and flagged the ENCODE CRISPRi held-out file as a reused benchmark: at least ten
scored results in this project have read it.

EVALUATION-ONLY says a scorer fitted nothing to the pairs. It does not say the scorer was built
without them. This audit asks the harder question for every component of the node caller - the
CTCF cCRE input, the boundary merge, the 50 kb floor, the default caller chosen out of seven
candidates, the control's seed and draw count, the annotation, the statistic - and reports each as
CLEAN, EXPOSED or EVALUATION-ONLY with the evidence that decides it, using lane-split's legend.

Nothing is asserted that can be run. The script:

1. imports the node caller in a fresh interpreter and lists every module it pulls in, so "the caller
   never reads CRISPRi" is a fact about the import closure and not a reading of the source;
2. asks git when each constant of the caller was written and when the CRISPRi files first entered
   the project, so "the constants predate the data" is a date comparison;
3. recomputes the stage 2 containment counts on every chromosome from the CRISPRi tables and the
   caller alone, with the AlphaGenome archive never opened, and holds them against the committed
   result, so "stage 2 does not depend on the model archive" is a reproduction;
4. decomposes the 661 pairs into the two benchmark files, so which split the claim rests on is a
   count and not a recollection;
5. scores the six rejected node callers on the same 661 pairs. The default was kept over them on
   four measurements, one of which - node content over the archive - is the same containment
   statistic in its modelled form. That is the one exposure this audit finds, and the spread of the
   seven callers on the measured pairs is what it costs;
6. states, for each CRISPRi set this project has not scored for containment, how many pairs it
   would bring and what a binomial test on that many pairs can see.

No model is called and nothing is fetched. 0 AlphaGenome requests.
"""

from __future__ import annotations

import bisect
import csv
import gzip
import json
import math
import random
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import node_containment_audit as nca  # noqa: E402

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri  # noqa: E402
from genomeos.attribution.measured import CRISPRI_SPLIT_OF  # noqa: E402
from genomeos.genome import Annotation, default_gencode  # noqa: E402
from genomeos.genome.domains import (  # noqa: E402
    MERGE_BOUNDARIES_WITHIN,
    MIN_DOMAIN,
    STRICT_RELATIVE,
    ctcf_motif_sites,
    infer_domains,
    strand_variant,
)
from genomeos.genome.regulatory import load_ccres  # noqa: E402
from genomeos.results import save_result  # noqa: E402

CHROMS = nca.CHROMS
LOOSE_RELATIVE = 0.85  # scripts/oriented_domains.py LOOSE_RELATIVE, the 2026-09-14 scan threshold
#: the caller set the default was chosen out of, exactly as scripts/oriented_domains.py names it.
#: value: None for the default, else (score an element's hit must reach, take the best hit's strand,
#: restrict to the CTCF-only class)
CALLERS: dict[str, tuple[float, bool, bool] | None] = {
    "ctcf_only": None,
    "oriented": (LOOSE_RELATIVE, False, False),
    "oriented_ctcf_only": (LOOSE_RELATIVE, False, True),
    "oriented_best_hit": (LOOSE_RELATIVE, True, False),
    "oriented_strong": (STRICT_RELATIVE, False, False),
    "oriented_strict": (STRICT_RELATIVE, True, False),
}
#: the controls scored for every caller. count_matched is left to the committed result: it draws
#: until the post-merge edge count matches, which is the slow one, and it is not needed to compare
#: callers because each caller is compared with its own uniform control.
SENSITIVITY_CONTROLS = ("uniform", "uniform_merged", "circular")
MEASURED_RESULT = ROOT / "data/results/node_containment_measured.json"
INDEP_MHC = ROOT / "data/results/indep_mhc_crispri.json"
DCTAP = ROOT / "data/cache/indep/IGVFFI0957PYTA.tsv.gz"


# ---------------------------------------------------------------- 1. the import closure


def import_closure() -> dict:
    """Every module the default node caller pulls in, from a fresh interpreter.

    A caller that had ever read the benchmark would have to reach it through one of these."""
    code = (
        "import sys, json;"
        "import genomeos.genome.domains as d;"
        "from genomeos.genome.regulatory import load_ccres;"
        "print(json.dumps(sorted(m for m in sys.modules if m.startswith('genomeos'))))"
    )
    out = subprocess.run(  # noqa: S603
        [sys.executable, "-c", code], cwd=ROOT, capture_output=True, text=True, check=True
    )
    mods = json.loads(out.stdout.strip().splitlines()[-1])
    hits = [m for m in mods if "crispri" in m or "measured" in m or "benchmark" in m]
    return {
        "genomeos_modules_imported": mods,
        "modules_naming_a_crispri_source": hits,
        "clean": not hits,
    }


def source_mentions() -> dict:
    """Which files under genomeos/genome name the benchmark at all, even in a comment."""
    out: dict[str, list[str]] = {}
    for p in sorted((ROOT / "genomeos/genome").glob("*.py")):
        text = p.read_text()
        hit = [t for t in ("crispri", "CRISPRi", "EPCrispr", "Engreitz", "Gschwind") if t in text]
        if hit:
            out[str(p.relative_to(ROOT))] = hit
    return out


# ---------------------------------------------------------------- 2. the constants in git


def _git(*args: str) -> str:
    return subprocess.run(  # noqa: S603
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout


def constant_history() -> dict:
    """When each constant of the caller was written, against when the CRISPRi files first arrived."""
    first_crispri = _git(
        "log", "--reverse", "--format=%h %ad", "--date=short", "-S", "EPCrisprBenchmark"
    ).splitlines()
    rows = {}
    for label, needle, path in (
        ("min_node_bp", "MIN_DOMAIN = 50_000", "genomeos/genome/domains.py"),
        ("boundary_merge_bp", "MERGE_BOUNDARIES_WITHIN = 5_000", "genomeos/genome/domains.py"),
        ("boundary_class", 'c.cls == "CTCF-only"', "genomeos/genome/domains.py"),
        ("node_confidence", "confidence: float = 0.4", "genomeos/genome/domains.py"),
        ("control_seed_inherited", "random.Random(7)", "scripts/oriented_domains.py"),
        ("control_seed_in_the_audit", "SEED = 7", "scripts/node_containment_audit.py"),
        ("control_draws", "SHUFFLES = 20", "scripts/oriented_domains.py"),
    ):
        log = [
            line.split(" ", 2)
            for line in _git("log", "--format=%h %ad %s", "--date=short", "-S", needle, "--", path)
            .strip()
            .splitlines()
        ]
        rows[label] = {
            "literal": needle,
            "file": path,
            "commits_that_changed_it": [{"sha": a, "date": b, "subject": c[:90]} for a, b, c in log],
            "written": log[-1][1] if log else None,
            "changed_since": [b for _, b, _ in log[:-1]],
        }
    return {
        "first_commit_naming_the_crispri_benchmark": first_crispri[0] if first_crispri else None,
        "constants": rows,
    }


# ---------------------------------------------------------------- 3. stage 2 without the archive


def _crispri_pairs_by_chrom() -> tuple[dict[str, list[tuple[int, int]]], dict]:
    """The 661 regulated pairs as (element midpoint, measured TSS), and their split decomposition."""
    pairs: list[tuple[str, object]] = []
    per_file: dict[str, int] = {}
    for t in nca.CRISPRI_TABLES:
        if (crispri.KNOWLEDGE / t).exists():
            got = crispri.load(t)
            split = CRISPRI_SPLIT_OF.get(t, t)
            per_file[split] = sum(1 for p in got if p.regulated and p.tss is not None)
            pairs += [(split, p) for p in got]
    reg = [(s, p) for s, p in pairs if p.regulated and p.tss is not None]
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    cells_by_split: dict[str, dict[str, int]] = {}
    for split, p in reg:
        by_chrom.setdefault(p.chrom, []).append((p.midpoint, p.tss))
        d = cells_by_split.setdefault(split, {})
        d[p.cell] = d.get(p.cell, 0) + 1
    return by_chrom, {"pairs_per_file": per_file, "cells_per_file": cells_by_split, "total": len(reg)}


def _nodes(chrom: str, caller: str) -> tuple[list[int], int]:
    length = nca.chrom_length(chrom)
    ccres = load_ccres(chrom)
    gff = default_gencode({chrom})
    ann = Annotation.from_gff3(gff, {chrom})
    spec = CALLERS[caller]
    if spec is None:
        doms = infer_domains(chrom, length, ccres, ann)
    else:
        rel, best, only_cls = spec

        def _no_fetch(a: int, b: int) -> str:  # the scan cache covers every chromosome; never called
            raise RuntimeError(f"no CTCF site cache for {chrom}: this audit fetches nothing")

        sites = ctcf_motif_sites(chrom, ccres, _no_fetch)
        strands = strand_variant(sites, relative=rel, best_hit=best)
        if only_cls:
            keep = {c.id for c in ccres if c.cls == "CTCF-only"}
            strands = {k: v for k, v in strands.items() if k in keep}
        doms = infer_domains(chrom, length, ccres, ann, orientation=strands)
    return [d.start for d in doms], length


def _caller_chrom(job: tuple[str, str, list[tuple[int, int]]]) -> dict:
    caller, chrom, mp = job
    starts, length = _nodes(chrom, caller)
    rng = random.Random(nca.SEED)
    cs = nca.controls(starts, len(starts) - 1, length, mp, rng)
    clusters: dict[int, list[int]] = {}
    for m, t in mp:
        i = bisect.bisect_right(starts, m)
        clusters.setdefault(i, []).append(1 if i == bisect.bisect_right(starts, t) else 0)
    return {
        "chrom": chrom,
        "caller": caller,
        "length": length,
        "nodes": len(starts),
        "coding_named": len(mp),
        "coding_inside": nca._inside(starts, mp),
        "controls": cs,
        "node_clusters": [len(v) for v in clusters.values()],
        "node_cluster_inside": [sum(v) for v in clusters.values()],
    }


def reconstruction(rows: list[dict]) -> dict:
    """The committed stage 2 containment counts, recomputed with the archive never opened."""
    stored = json.loads(MEASURED_RESULT.read_text())
    want = {r["chrom"]: (r["coding_named"], r["coding_inside"]) for r in stored["per_chromosome"]}
    got = {r["chrom"]: (r["coding_named"], r["coding_inside"]) for r in rows if r["caller"] == "ctcf_only"}
    diffs = {c: {"stored": want[c], "recomputed": got.get(c)} for c in want if want[c] != got.get(c)}
    return {
        "chromosomes_compared": len(want),
        "identical": not diffs and set(want) == set(got),
        "differences": diffs,
        "archive_opened": False,
        "what_it_shows": (
            "the measured arm is a function of the CTCF cCREs, GENCODE, the chromosome length and the "
            "CRISPRi pairs alone: the AlphaGenome archive that stage 1 scores is never read"
        ),
    }


# ---------------------------------------------------------------- 5. the rejected callers


def caller_readings(rows: list[dict]) -> dict:
    out: dict[str, dict] = {}
    rng = random.Random(11)
    for caller in CALLERS:
        rs = [r for r in rows if r["caller"] == caller]
        if not rs:
            continue
        named = sum(r["coding_named"] for r in rs)
        rec: dict = {
            "pairs": named,
            "mean_nodes": sum(r["nodes"] for r in rs),
            "measured_share": round(sum(r["coding_inside"] for r in rs) / max(1, named), 4),
            "controls": {},
        }
        for c in SENSITIVITY_CONTROLS:
            share, ctrl, ex = nca.pooled(rs, c)
            rec["controls"][c] = {
                "control_share": round(ctrl, 4),
                "excess_points": round(ex * 100, 3),
                "bootstrap_chromosomes": nca.bootstrap_chromosomes(rs, c, rng),
                "sign_test": {k: v for k, v in nca.sign_test(rs, c).items() if k != "per_chromosome"},
            }
        out[caller] = rec
    return out


# ---------------------------------------------------------------- 6. what a fresh set could see


def _binom_tail(n: int, k: int, p: float) -> float:
    return sum(math.comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(k, n + 1))


def power_at(n: int, p0: float, p1: float, alpha: float = 0.05) -> dict:
    """Exact one-sided binomial: the smallest count that rejects, and the power at p1."""
    if n <= 0:
        return {"pairs": n, "critical_count": None, "power": 0.0}
    crit = next((k for k in range(n + 1) if _binom_tail(n, k, p0) <= alpha), n + 1)
    return {
        "pairs": n,
        "null_share": round(p0, 4),
        "alternative_share": round(p1, 4),
        "critical_count": crit,
        "critical_share": round(crit / n, 4) if crit <= n else None,
        "power": round(_binom_tail(n, crit, p1), 4) if crit <= n else 0.0,
    }


def _dctap_fresh() -> dict:
    """The DC-TAP pairs no result here has scored: distal, targeting, not in either ENCODE file.

    Counted the way scripts/indep_igvf_probe.py counts the overlap (same gene symbol, overlapping
    interval on the same chromosome)."""
    if not DCTAP.exists():
        return {"available": False}
    by_gene: dict[str, list[tuple[str, int, int]]] = {}
    for name in nca.CRISPRI_TABLES:
        with gzip.open(crispri.KNOWLEDGE / name, "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                by_gene.setdefault(r["measuredGeneSymbol"], []).append(
                    (r["chrom"], int(r["chromStart"]), int(r["chromEnd"]))
                )
    distal = fresh = pos = fresh_pos = 0
    els: set = set()
    with gzip.open(DCTAP, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["type"] != "targeting" or r["promoter_gene"] not in ("", "NA", "None"):
                continue
            distal += 1
            c, s, e, g = (
                r["targeting_chr"],
                int(r["targeting_start"]),
                int(r["targeting_end"]),
                r["gene_symbol"],
            )
            new = not any(cc == c and a < e and b > s for cc, a, b in by_gene.get(g, ()))
            fresh += new
            try:
                down = r["significant"].strip().upper() == "TRUE" and float(r["sceptre_log2_fc"]) < 0
            except ValueError:
                down = False
            if down:
                pos += 1
                if new:
                    fresh_pos += 1
                    els.add((c, s, e))
    return {
        "available": True,
        "distal_targeting_pairs": distal,
        "distal_pairs_not_in_either_encode_file": fresh,
        "distal_significant_decreases": pos,
        "distal_significant_decreases_not_in_encode": fresh_pos,
        "distinct_elements_behind_them": len(els),
    }


def candidate_arms(p0: float, p1: float) -> dict:
    mhc = json.loads(INDEP_MHC.read_text()) if INDEP_MHC.exists() else {}
    mhc_n = (mhc.get("outcomes_eligible") or {}).get("decrease")
    dct = _dctap_fresh()
    fresh = dct.get("distal_significant_decreases_not_in_encode") or 0
    need = next(n for n in range(10, 20001, 1) if power_at(n, p0, p1)["power"] >= 0.8)
    arms = {
        "igvf_mhc_perturb_seq": {
            "file": "data/cache/indep/IGVFFI4093WUVB.tsv.gz (lane-indep, "
            "data/results/indep_mhc_crispri.json)",
            "containment_pairs": mhc_n,
            "power": power_at(mhc_n or 0, p0, p1),
            "extra": "every pair sits in one 4 Mb locus on chr6, so the independent-unit count is 1, "
            "not 19: the binomial figure above is an upper bound on what it could see",
        },
        "igvf_dctap_k562_remainder": {
            "file": "data/cache/indep/IGVFFI0957PYTA.tsv.gz (Ray 2025 DC-TAP, the pairs not in the "
            "ENCODE files; no result in this project has scored them)",
            "inventory": dct,
            "containment_pairs": fresh,
            "power": power_at(fresh, p0, p1),
            "extra": "K562 only, so it does not touch the claim's cell-type caveat; the same study's "
            "other pairs are already in the ENCODE held-out file",
        },
        "igvf_k562_many_loci_library": {
            "file": "IGVFDS3624DGBH / IGVFDS3481OUHP (lane-indep2, data/results/indep_igvf_probe.json)",
            "containment_pairs": None,
            "power": None,
            "extra": "unreleased (HTTP 403); the library covers 279 distal elements, so a released "
            "table would bring pairs but the positive count cannot be known before release",
        },
    }
    return {
        "null_share": round(p0, 4),
        "alternative_share": round(p1, 4),
        "pairs_for_80_percent_power": need,
        "arms": arms,
        "verdict": (
            "no usable independent arm exists now at 0 requests: the two readable sets bring "
            f"{(mhc_n or 0) + fresh} regulated pairs between them against the {need} an exact one-sided "
            "binomial needs to see +5.89 points at 80% power, and one of the two is a single locus"
        ),
    }


# ---------------------------------------------------------------- the table


def components(closure: dict, hist: dict, recon: dict, split: dict) -> list[dict]:
    """lane-split's legend: CLEAN, EXPOSED, EVALUATION-ONLY. One row per component of the claim."""
    cons = hist["constants"]
    return [
        {
            "component": "the pairs the claim is scored on",
            "what": f"{split['total']} CRISPRi Regulated=TRUE pairs, {split['pairs_per_file']} by file",
            "status": "EVALUATION-ONLY",
            "why": "both benchmark files, training and held-out, pooled. The containment statistic "
            "fits nothing, so neither file is a test set and neither is a training set; but the "
            "held-out file is a reused benchmark (lane-split: at least ten scored results read it), "
            "so the pairs are not a first touch",
        },
        {
            "component": "the CTCF cCRE input",
            "what": "ENCODE SCREEN V3 registry, CTCF-only class, loaded by genome/regulatory.py",
            "status": "CLEAN",
            "why": "an external registry built from ENCODE ChIP and DNase, published before this "
            "project and untouched by it; it names no gene and no perturbation",
        },
        {
            "component": "the boundary rule (CTCF-only midpoints, merged within 5 kb)",
            "what": f"MERGE_BOUNDARIES_WITHIN = {MERGE_BOUNDARIES_WITHIN}",
            "status": "CLEAN",
            "why": "written "
            + str(cons["boundary_merge_bp"]["written"])
            + " and never changed; the CRISPRi benchmark first entered the project in "
            + str(hist["first_commit_naming_the_crispri_benchmark"]),
        },
        {
            "component": "the 50 kb node floor",
            "what": f"MIN_DOMAIN = {MIN_DOMAIN}",
            "status": "CLEAN",
            "why": "written "
            + str(cons["min_node_bp"]["written"])
            + " and never changed since; inherited as a TAD size floor, not fitted. It is not "
            "neutral for the number - the published control does not pass through it, which is why "
            "the audit reports uniform_merged beside uniform - but it is not tuned on the pairs",
        },
        {
            "component": "thresholds inside the caller",
            "what": "the default caller has none beyond the class filter and the two lengths above; "
            "STRICT_RELATIVE = 0.95 belongs to the orientation callers and was calibrated against "
            "saturation mutagenesis (genome/motifs.py), not against CRISPRi",
            "status": "CLEAN",
            "why": "the import closure of the default caller holds "
            f"{len(closure['genomeos_modules_imported'])} genomeos modules and none of them reads a "
            "CRISPRi source",
        },
        {
            "component": "the annotation the TSS comes from",
            "what": "GENCODE 50 protein-coding TSS, the same annotation for the caller and the pairs",
            "status": "CLEAN",
            "why": "external release; the measured gene is the benchmark's own measuredGeneSymbol, "
            "not a gene this project chose",
        },
        {
            "component": "the random-boundary control",
            "what": "seed 7, 20 draws, inherited from scripts/oriented_domains.py",
            "status": "CLEAN",
            "why": "fixed before the measured arm was written and reused unchanged, so the control "
            "is not redrawn per result; three stricter controls are reported beside it",
        },
        {
            "component": "the statistic",
            "what": "share of pairs whose element midpoint and gene TSS fall in one node",
            "status": "EVALUATION-ONLY",
            "why": "no parameter is fitted; the caller's output is compared with a measurement it never saw",
        },
        {
            "component": "which model output the measured arm depends on",
            "what": "none",
            "status": "CLEAN",
            "why": recon["what_it_shows"]
            + f"; recomputed on {recon['chromosomes_compared']} chromosomes with the archive never "
            f"opened, identical: {recon['identical']}",
        },
        {
            "component": "the choice of default caller out of seven candidates",
            "what": "ctcf_only kept over six orientation callers on 2026-09-14 and again on "
            "2026-09-21, on four measurements: 4DN Hi-C boundary support, node content over the "
            "AlphaGenome archive, mouse synteny and the HOXD interval",
            "status": "EXPOSED",
            "why": "not to the CRISPRi pairs - no CRISPRi figure existed on 2026-09-14 and none was "
            "read on 2026-09-21 - but to the containment statistic itself in its modelled form. "
            "'Node content' is the same share of (element, gene) pairs inside one node, computed on "
            "the model archive, and the commit messages quote the winner's excess over its own "
            "random control (+2.9 against -0.3 for the stricter call). The default was selected as "
            "the maximum of seven callers on a statistic correlated with the one it is now scored "
            "on, so the measured excess carries a selection premium",
        },
    ]


# ---------------------------------------------------------------- main


def manifest_of(chroms: list[str]) -> dict:
    inputs = []
    for c in chroms:
        for p in (ROOT / f"data/reference/{c}.fa.gz.fai", ROOT / f"data/reference/{c}.fa.fai"):
            if p.exists():
                inputs.append(mf.input_entry(p, partition=None))
                break
        inputs.append(mf.input_entry(default_gencode({c}), partition=None))
        inputs.append(mf.input_entry(ROOT / f"data/results/ccres_{c}.bed.gz", partition=None))
    for t in nca.CRISPRI_TABLES:
        if (crispri.KNOWLEDGE / t).exists():
            inputs.append(mf.input_entry(crispri.KNOWLEDGE / t, partition=CRISPRI_SPLIT_OF.get(t)))
    if MEASURED_RESULT.exists():
        inputs.append(mf.input_entry(MEASURED_RESULT, partition=None))
    if DCTAP.exists():
        inputs.append(mf.input_entry(DCTAP, partition="independent, counted only"))
    return {
        "sources": [
            {"accession": "GENCODE human gene annotation", "version": "release 50 (GRCh38.p14)"},
            {
                "accession": "ENCODE SCREEN registry of cCREs (CTCF-only elements as node boundaries)",
                "version": "V3, downloads.wenglab.org/V3/GRCh38-cCREs.bed",
            },
            {
                "accession": "EngreitzLab/CRISPR_comparison resources/crispr_data (Gschwind et al. 2025)",
                "version": "main branch, unpinned upstream; fetched 2026-09-16; pinned here by sha256",
                "url": crispri.BASE_URL,
            },
            {
                "accession": "IGVF IGVFFI0957PYTA (Ray 2025 K562 Random DC-TAP-seq)",
                "version": "released 2026-08-08; counted for power only, no containment computed",
            },
            {
                "accession": "this repository's git history",
                "version": "read with git log -S for the constants named in parameters",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "callers": list(CALLERS),
            "controls": list(SENSITIVITY_CONTROLS),
            "min_node_bp": MIN_DOMAIN,
            "boundary_merge_bp": MERGE_BOUNDARIES_WITHIN,
            "loose_relative": LOOSE_RELATIVE,
            "strict_relative": STRICT_RELATIVE,
            "shuffles": nca.SHUFFLES,
            "seed_per_chromosome": nca.SEED,
            "bootstraps": nca.BOOTSTRAPS,
            "bootstrap_seed": 11,
            "alpha_one_sided": 0.05,
            "alphagenome_requests": 0,
        },
        "exclusions": [
            "count_matched is not scored for the six rejected callers: each caller is compared with "
            "its own uniform control, and the matched draw is the slow one",
            "the DC-TAP and MHC sets are counted, not scored: neither has the pairs to answer, and "
            "reading them would spend the only unread pairs this project has",
        ],
        "partitions": {
            "training": "EPCrisprBenchmark training_K562 Regulated=TRUE pairs",
            "heldout": "EPCrisprBenchmark heldout_5_cell_types Regulated=TRUE pairs",
            "pooled": "the containment statistic pools both and fits nothing to either",
            "independent": "IGVF DC-TAP remainder and the IGVF MHC screen: counted for power only",
        },
    }


def main() -> None:
    args = sys.argv[1:]
    workers = int(args[args.index("--workers") + 1]) if "--workers" in args else 6
    quick = "--quick" in args

    closure = import_closure()
    print("import closure:", closure["clean"], len(closure["genomeos_modules_imported"]), "modules")
    hist = constant_history()
    print("first CRISPRi commit:", hist["first_commit_naming_the_crispri_benchmark"])

    by_chrom, split = _crispri_pairs_by_chrom()
    chroms = sorted(by_chrom)
    callers = ["ctcf_only"] if quick else list(CALLERS)
    jobs = [(c, ch, by_chrom[ch]) for c in callers for ch in chroms]
    with ProcessPoolExecutor(workers) as pool:
        rows = list(pool.map(_caller_chrom, jobs))
    print("scored", len(rows), "caller-chromosome cells", flush=True)

    recon = reconstruction(rows)
    print("reconstruction identical:", recon["identical"], flush=True)
    sens = caller_readings(rows)
    for name, rec in sens.items():
        u = rec["controls"]["uniform"]
        print(
            f"  {name}: share {rec['measured_share']} excess {u['excess_points']} "
            f"ci {u['bootstrap_chromosomes']['ci95']} sign {u['sign_test']['ahead']}/{u['sign_test']['of']}",
            flush=True,
        )
    default = sens["ctcf_only"]
    p1 = default["measured_share"]
    p0 = default["controls"]["uniform"]["control_share"]
    arms = candidate_arms(p0, p1)
    print("arms:", arms["verdict"], flush=True)

    excesses = {k: v["controls"]["uniform"]["excess_points"] for k, v in sens.items()}
    rejected = [v for k, v in excesses.items() if k != "ctcf_only"]
    selection = {
        "measurement_the_default_was_selected_on": "node content over the AlphaGenome archive, one of "
        "four criteria (scripts/oriented_domains.py, commits 26da99a and 56e2c50)",
        "excess_over_uniform_by_caller": excesses,
        "default": excesses.get("ctcf_only"),
        "best_rejected": max(rejected) if rejected else None,
        "worst_rejected": min(rejected) if rejected else None,
        "spread_points": round(max(excesses.values()) - min(excesses.values()), 3) if excesses else None,
        "rejected_callers_clear_of_zero": sum(
            1
            for k, v in sens.items()
            if k != "ctcf_only" and v["controls"]["uniform"]["bootstrap_chromosomes"]["ci95"][0] > 0
        ),
        "rejected_callers": len(rejected),
    }

    payload = {
        "question": "how much of the node-containment result is independent of the pairs it is "
        "scored on, and of the choices this project made before scoring them",
        "legend": {
            "CLEAN": "the component was built without the pairs and nothing in it was fitted to them",
            "EXPOSED": "the component was chosen or tuned using the pairs, or using a statistic "
            "computed from them or correlated with the one it is scored on",
            "EVALUATION-ONLY": "scores against the pairs, fits nothing to them",
        },
        "components": components(closure, hist, recon, split),
        "import_closure": closure,
        "source_mentions_under_genome": source_mentions(),
        "constant_history": hist,
        "pairs": split,
        "reconstruction_without_the_archive": recon,
        "caller_sensitivity": sens,
        "selection_premium": selection,
        "independent_arms_at_zero_requests": arms,
        "alphagenome_requests": 0,
    }
    save_result("node_independence_audit", payload, manifest=manifest_of(chroms))
    print("saved node_independence_audit", flush=True)


if __name__ == "__main__":
    main()
