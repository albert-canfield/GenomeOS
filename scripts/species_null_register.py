#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The species-tree null's registration: which area J readings are re-tested, how, and in what words.

    uv run python scripts/species_null_register.py [--fetch]

Written after the code is frozen in a commit and before any J reading is computed under either null.
It names the readings verbatim from GRAMMAR-BY-COMPARISON (and refuses to run if a quote is no longer
in the document word for word), the result files they live in with their digests, the species tree
(source files, digests, the merge rule and the pruned topology's digest), both nulls and every
parameter, the statistic, the falsifier, the decision rule and the exact wording of every outcome.
It holds no figure the run would fill: FORBIDDEN is checked at any depth before the file is written.

`--fetch` downloads the two Ensembl species trees (about 13 KB) when they are absent and checks them
against the digests pinned in `genomeos/knowledge/species_null.py`.

`scripts/species_null_run.py` refuses to start unless this file is committed, the code it pins is
the code on disk, and the tree and inputs still hash to what is recorded here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.knowledge import species_null as sn  # noqa: E402
from genomeos.results import load_result, save_result  # noqa: E402

NAME = "species_null_registration"
ENTRY = "scripts/species_null_register.py"
OWN_CODE = (
    "genomeos/knowledge/species_null.py",
    "scripts/species_null_register.py",
    "scripts/species_null_run.py",
    "tests/test_species_null.py",
)
#: The code the run must be: pinned by blob, so an edit after registration stops the run.
FROZEN = ("genomeos/knowledge/species_null.py", "scripts/species_null_run.py")
DOC = "docs/GRAMMAR-BY-COMPARISON.md"
INPUTS = {
    "origin_presence_genome_wide": "the presence matrix both nulls read: one bit per gene per species",
    "origin_genome_wide": "each gene's origin stratum, which §14's age-matched draw and null A's age "
    "expectation read",
    "library_members": "the computed membership of the libraries, the same for §14 and §19 (checked "
    "identical for all 45 libraries on 2026-10-10)",
    "profiling_genome_wide": "§14's stored pairs, read for the reproduction gate and for which 35 §14 "
    "called; no figure of it is an input to either null",
    "profiling_genes": "§19's stored libraries, read for the reproduction gate and for which 21 §19 "
    "called clear; no figure of it is an input to either null",
}

READINGS = {
    "A_library_pairs": {
        "section": "GRAMMAR-BY-COMPARISON §14 (2026-09-12), with §11's headline pairs through their §14 rows",
        "result_file": "data/results/profiling_genome_wide.json (null.pairs, 860 tested; "
        "null.significant_q05, 35 rows)",
        "quotes": [
            "**What is more alike than age predicts.** 35 of 860 pairs at a false discovery rate of 5%, "
            "led by the immune and signalling libraries: immune against the ligand-receptor protocol "
            "(observed −0.76, age-matched −0.93, z 9.3), degradation and cytoskeleton against the same "
            "protocol (z 8.5), the cell cycle against immunity (z 6.7), blood with the growth-longevity "
            "axis (observed +0.72 against −0.27).",
            "**Both headline findings of §11 are age.** Random age-matched sets already correlate at "
            "0.55 to 0.73 where the developmental libraries sit, so axes with segmentation reads z 1.9, "
            "heart with musculoskeletal 2.2, germ layers with the nervous system 1.6, none significant. "
            "The cores against the ligand-receptor protocol match their null (replication −0.87 against "
            "−0.79, z −0.9).",
        ],
        "why_the_tree_changes_it": "§14's null holds age and not the tree: a random age-matched set's "
        "profile is close to its age class's average, so the null's spread is small, while a library's "
        "own departure from its age expectation is one trait evolving along the tree and agrees across a "
        "clade because the clade shares an ancestor. A z of 9.3 against that spread counts each of 355 "
        "species as an observation.",
        "named_rows_reported_individually": [
            ["systems.immune", "systems.ligand_receptor_protocol"],
            ["core.degradation", "systems.ligand_receptor_protocol"],
            ["core.cytoskeleton", "systems.ligand_receptor_protocol"],
            ["core.cell_cycle", "systems.immune"],
            ["blueprint.organ_blood_immune", "timer.growth_longevity_axis"],
            ["blueprint.axes", "blueprint.segmentation"],
            ["blueprint.organ_heart", "blueprint.organ_musculoskeletal"],
            ["blueprint.germ_layers", "blueprint.organ_nervous"],
            ["core.replication", "systems.ligand_receptor_protocol"],
        ],
    },
    "B_library_coherence": {
        "section": "GRAMMAR-BY-COMPARISON §19 (2026-09-17)",
        "result_file": "data/results/profiling_genes.json (libraries, 42; matched_draws_at_or_above_"
        "observed 0 for 21)",
        "quotes": [
            "| clear of the null (0 draws at or above) | **21** | | marginal (1–9 draws above) | 7 | "
            "| **at or below the null (≥10 draws above)** | **14** |",
            "**The ones that survive are patchy and lineage-restricted**: `timer.developmental_timing` "
            "(+0.0413), `systems.extracellular_matrix_adhesion` (+0.0393), `blueprint.segmentation` and "
            "`timer.segmentation_clock` (+0.0381 each).",
            "Even the matched-gene null draws genes independently of the tree, so a clade-wide loss "
            "counts once in the library and once in each draw but the *correlation structure* between "
            "the drawn genes is whatever the genome supplies rather than what a model of the tree would "
            "say. Felsenstein contrasts or a birth-death model along a species tree is what would "
            "settle a pair, and this project holds no tree.",
        ],
        "why_the_tree_changes_it": "§19's matched genes are matched on prevalence, not on where on the "
        "tree a gene was born: two genes born in the same clade share every absence outside it by "
        "descent, and a library whose members share a birth clade co-occurs more than prevalence-"
        "matched genes of mixed births without any loss in common.",
    },
}
NOT_RETESTED = {
    "§10 origin per gene and age per library": "a placement (the deepest clade with an orthologue), "
    "not agreement between profiles; shared ancestry is what it reads, so a null of shared ancestry "
    "does not apply. The tree's birth nodes are reported beside the strata for the record only.",
    "§11 residual correlations on full membership": "§14 already re-read them on exclusive members "
    "against age; they are re-tested through their §14 rows (named_rows_reported_individually).",
    "§14 pairs far less alike than age (z down to −35)": "recorded and not called by §14; the test "
    "stays one-sided, more alike.",
    "§19 null one, the prevalence shuffle": "§19 itself rejects it as measuring non-independence.",
    "gene pairs": "§19 called no gene pair a dependency, and 188 million gene pairs are a separate "
    "test with its own multiple-testing design; not attempted here.",
    "§12, §13 and §15 to §21": "one locus across species, families and composition, and the enhancer, "
    "assay and sequence readings: none counts agreement of presence across species.",
}

NULLS = {
    "tree": {
        "sources": [
            {"file": sn.PAN_TREE.name, "url": sn.PAN_TREE_URL, "sha256": sn.TREE_SHA256[sn.PAN_TREE.name]},
            {"file": sn.VERT_TREE.name, "url": sn.VERT_TREE_URL, "sha256": sn.TREE_SHA256[sn.VERT_TREE.name]},
        ],
        "why_these": "Ensembl Compara's own species trees, the trees its protein trees were built on, "
        "so the tree is the one the orthologue calls of the presence matrix already assume; release 63 "
        "of Ensembl Genomes is the pan-taxonomic counterpart of Ensembl 116, the release of both "
        "homology dumps the presence matrix was streamed from",
        "merge": "the vertebrate tree's Vertebrata (the common ancestor of human and sea lamprey) "
        "replaces the pan tree's vertebrate clade; Ciona intestinalis, an outgroup only the vertebrate "
        "tree carries, is placed as sister to Ciona savignyi, as the vertebrate tree places it; any "
        "other conflict raises",
        "species_matching": "production name, then without the assembly suffix, then the RENAMES "
        "table (12 renamed taxa), then a unique tip extending the name or the genus and species; two "
        "candidates or two species on one tip raise",
        "renames": sn.RENAMES,
        "branch_lengths": {
            "registered": "Grafen (1989), rho 1: a node's height is its number of descendant tips "
            "minus one, the root's scaled to 1. Neither file has lengths both halves share (the pan "
            "tree has none), so a length is a modelling choice; Grafen's is the conventional one for a "
            "topology",
            "sensitivity": "every branch length 1; reported beside the registered reading and "
            "deciding nothing",
        },
    },
    "A": {
        "statistic": "§14's own: Pearson r over the 355 species between two libraries' residual "
        "profiles (share of the members exclusive to each side with an orthologue in the species, "
        "minus the share of all 19,391 genes), r_exclusive, recomputed and required to match the stored "
        "value",
        "model": "each side's profile P = E + D, with E the age expectation (the exact mean profile of "
        "§14's own matched draw: per origin stratum, the members' count times the mean profile of every "
        "gene of that stratum, over the members with an origin) and D the deviation. Under the null E is "
        "held and D is replaced by an independent Brownian motion along the species tree, D ~ N(mu, "
        "sigma^2 C), with C the tree's covariance (shared root-to-ancestor length) and sigma^2 the REML "
        "estimate from D itself, (D - mu_GLS)' C^-1 (D - mu_GLS) / (n - 1); mu is irrelevant to r. The "
        "two sides' motions are independent; the statistic is recomputed on E + D_sim for each draw",
        "parameters": {
            "draws": sn.A_SIMS,
            "seed": sn.A_SEED,
            "p": "(number of draws with r at or above the observed + 1) / (draws + 1), one-sided, more alike",
            "multiple_testing": "Benjamini-Hochberg over every testable pair, q at or below 0.05, as "
            "§14 did; the draws are shared across pairs (two banks), which keeps each p exact and makes "
            "different pairs' p values dependent",
            "min_exclusive_members": sn.MIN_MEMBERS,
        },
        "calibration": {
            "what": "one random pair of gene sets per testable pair, each set drawn to its side's "
            "origin-stratum make-up (§14's own draw), tested exactly as a library pair",
            "seed": sn.A_CAL_SEED,
            "ceiling": f"more than {sn.CAL_MAX_SHARE:.0%} of them at p at or below 0.05 and the null is "
            "not calibrated: no reading",
        },
        "power": {
            "what": "each library split at random into two halves, tested as a pair: the halves share "
            "the library's history, so a null that never sees them as alike cannot see a pair",
            "seed": sn.A_HALVES_SEED,
            "floor": f"fewer than {sn.HALVES_MIN_SHARE:.0%} of libraries at p at or below 0.05 and a "
            "result with no survivor is read as low power",
        },
        "falsifier": "§14's reading that a pair is more alike than age predicts is re-read as age plus "
        "shared ancestry for every one of its 35 pairs whose q under this null is above 0.05, provided "
        "the null is calibrated; if none of the 860 survives, the 35 are what age and shared ancestry "
        "give on this model",
        "descriptive_only": "the GLS (phylogenetically whitened) correlation of the two deviations, "
        "the null's mean and sd, and the equal-length run; none decides",
        "decision_rule": sn.A_RULE,
        "outcome_wordings": sn.A_OUTCOMES,
        "sensitivity_wording": sn.A_SENSITIVITY,
    },
    "B": {
        "statistic": "§19's per-library mean Jaccard over all member pairs of their presence sets "
        "(exact over every pair, where §19 sampled 60,000 pairs for the largest libraries), less the "
        "mean of the same over the model's draws of the members: the library's excess over its own "
        "tree model",
        "model": "stochastic Dollo along the species tree with human added: a gene is born at the "
        "common ancestor of human and the species carrying it and lost independently on the branches "
        "below, never on the path to human. A branch's loss probability p_e is (losses + "
        f"{sn.LOSS_PSEUDO[0]}) / (genes at risk + {sn.LOSS_PSEUDO[1]}) from every gene's Dollo "
        "reconstruction (at risk: the gene is held at the branch's parent); a gene's survival on a "
        "branch is (1 - p_e)^theta_g, theta_g fitted by bisection so that its expected number of "
        "species equals its observed number",
        "null": "the same excess for random gene sets matched member for member on birth node, drawn "
        "without replacement per node; the model removes birth clade and loss rate, the real genes "
        "remove what the model misses for every gene alike",
        "parameters": {
            "model_draws_per_gene": sn.B_SIMS,
            "model_seed": sn.B_SIM_SEED,
            "random_sets_per_library": sn.B_SETS,
            "set_seed": sn.B_SEED,
            "min_members": sn.MIN_MEMBERS,
            "reproduction_tolerance": sn.B_TOLERANCE,
            "chance_ceiling_for_42": sn.chance_ceiling(42),
        },
        "classes": "clear: (k + 1) / (B_SETS + 1) at or below 0.05, k the random sets whose excess is "
        "at or above the library's; at_or_below: 2k at or above B_SETS; marginal otherwise. With 20 "
        "sets these are exactly §19's classes",
        "falsifier": "§19's reading that a library's members co-occur beyond prevalence-matched genes "
        "is re-read as birth clade and loss rate on the tree for every one of its 21 clear libraries "
        "that is not clear here",
        "descriptive_only": "the number of random sets whose observed mean Jaccard exceeds every model "
        "draw (how far the model misses real genes), and the libraries' excesses; none decides",
        "decision_rule": sn.B_RULE,
        "outcome_wordings": sn.B_OUTCOMES,
    },
}

#: Keys no registration may carry at any depth: each is a figure only the run produces.
FORBIDDEN = (
    "p_tree",
    "q_tree",
    "survivors",
    "n_surv",
    "k_14",
    "k_19",
    "n_clear",
    "class",
    "excess_over_tree_model",
    "outcome",
    "verdict",
    "reading",
    "calibration_result",
    "halves_result",
)


def forbidden_at_any_depth(payload: Any, path: str = "") -> list[str]:
    found: list[str] = []
    if isinstance(payload, dict):
        for key, value in payload.items():
            here = f"{path}.{key}" if path else key
            if key in FORBIDDEN:
                found.append(here)
            found.extend(forbidden_at_any_depth(value, here))
    elif isinstance(payload, list):
        for i, value in enumerate(payload):
            found.extend(forbidden_at_any_depth(value, f"{path}[{i}]"))
    return found


def git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    ).stdout.strip()


def frozen_code() -> dict[str, Any]:
    """The commit and blobs the run must find; refuses when any of them is not committed as it stands."""
    out: dict[str, Any] = {"head": git("rev-parse", "HEAD"), "blobs": {}}
    for path in FROZEN:
        if git("status", "--porcelain", "--", path):
            raise SystemExit(f"{path} is modified or untracked: freeze the code in a commit first")
        out["blobs"][path] = git("rev-parse", f"HEAD:{path}")
        out.setdefault("last_commit", {})[path] = git("log", "-1", "--format=%H", "--", path)
    return out


def check_quotes(doc_text: str) -> None:
    flat = re.sub(r"\s+", " ", doc_text)
    flat_table = re.sub(r"\|\s*\n\s*\|", "| |", doc_text)
    flat_table = re.sub(r"\s+", " ", flat_table)
    for key, reading in READINGS.items():
        for q in reading["quotes"]:
            if q not in flat and q not in flat_table:
                raise SystemExit(f"{key}: a quote is not in {DOC} word for word: {q[:80]}...")


def registration(tree_record: dict[str, Any], doc_blob: str) -> dict[str, Any]:
    return {
        "lane": "lane-speciesnull",
        "item": "ROADMAP §5 item 13 (M6), area J: a species-tree null for the grammar by comparison",
        "registered_in_its_own_commit": "Written after the code was frozen in a commit and before any "
        "J reading was computed under either null. The run reads each reading once, against this file.",
        "question": "Two of area J's readings count agreement of orthologue presence species by species. "
        "Which of them survive a null in which the agreement comes from shared ancestry alone?",
        "readings_retested": READINGS,
        "readings_not_retested": NOT_RETESTED,
        "document": {"path": DOC, "blob_at_registration": doc_blob},
        "nulls": NULLS,
        "tree_as_built": tree_record,
        "gates": [
            "this file committed and unmodified",
            "genomeos/knowledge/species_null.py and scripts/species_null_run.py at the blobs pinned in "
            "frozen_code",
            "the two tree files at their pinned sha256, and the pruned topologies at the digests below",
            "every input result at the sha256 recorded in result_manifest.inputs",
            "the reproduction rule of each null (A_RULE 1_reproduce, B_RULE 1_reproduce)",
        ],
        "frozen_code": frozen_code(),
    }


def tree_record(species: list[str]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for label, with_human in (("without_human_for_A", False), ("with_human_for_B", True)):
        root, mapping = sn.species_tree(species, with_human=with_human)
        nodes = sn.preorder(root)
        out[label] = {
            "tips": len(sn.leaves(root)),
            "nodes": len(nodes),
            "internal_nodes_with_more_than_two_children": sum(1 for v in nodes if len(v.children) > 2),
            "topology_sha256": sn.topology_digest(root),
            "newick": sn.to_newick(root),
        }
        out["species_to_tip"] = mapping
    return out


def result_manifest(input_paths: list[Path]) -> dict[str, Any]:
    return {
        "sources": [
            {
                "accession": "Ensembl Genomes pan-taxonomic Compara species tree",
                "version": "release 63, " + sn.PAN_TREE_URL,
            },
            {
                "accession": "Ensembl Compara vertebrate species tree",
                "version": "release 116, " + sn.VERT_TREE_URL,
            },
            {"accession": "the project's area J results", "version": "as committed, digested in inputs"},
        ],
        "inputs": [
            mf.input_entry(sn.TREE_DIR / sn.PAN_TREE.name, partition="the species tree, pan-taxonomic half"),
            mf.input_entry(sn.TREE_DIR / sn.VERT_TREE.name, partition="the species tree, vertebrate half"),
            *[mf.input_entry(p, partition=INPUTS[p.stem]) for p in input_paths],
        ],
        "assembly": "n/a: species presence of orthologues and species trees; no genomic coordinate",
        "coordinates": "n/a: no interval is computed or read",
        "parameters": {
            "A_SIMS": sn.A_SIMS,
            "A_SEED": sn.A_SEED,
            "A_CAL_SEED": sn.A_CAL_SEED,
            "A_HALVES_SEED": sn.A_HALVES_SEED,
            "CAL_MAX_SHARE": sn.CAL_MAX_SHARE,
            "HALVES_MIN_SHARE": sn.HALVES_MIN_SHARE,
            "B_SIMS": sn.B_SIMS,
            "B_SETS": sn.B_SETS,
            "B_SEED": sn.B_SEED,
            "B_SIM_SEED": sn.B_SIM_SEED,
            "LOSS_PSEUDO": list(sn.LOSS_PSEUDO),
            "B_TOLERANCE": sn.B_TOLERANCE,
            "MIN_MEMBERS": sn.MIN_MEMBERS,
            "Q_ALPHA": sn.Q_ALPHA,
        },
        "exclusions": [
            "every figure either null would produce: none is computed by this writer",
            "the species list is read from the presence matrix to build the tree; no profile is read, "
            "and the other four results are hashed so the run can refuse inputs that moved",
        ],
        "partitions": "n/a: a registration written before the measurement; nothing is evaluated",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fetch", action="store_true", help="download the two species trees if absent")
    args = ap.parse_args(argv)
    if args.fetch:
        print(json.dumps(sn.fetch_trees(), indent=1))
    mf.trace_begin()
    for name in (sn.PAN_TREE.name, sn.VERT_TREE.name):
        path = sn.TREE_DIR / name
        if not path.exists():
            raise SystemExit(f"no {path}: run with --fetch")
        if sn.sha256(path) != sn.TREE_SHA256[name]:
            raise SystemExit(f"{path} does not hash to its pin")
    check_quotes((ROOT / DOC).read_text())
    presence = load_result("origin_presence_genome_wide")
    if not presence:
        raise SystemExit("no origin_presence_genome_wide result")
    record = registration(tree_record(presence["species"]), git("rev-parse", f"HEAD:{DOC}"))
    bad = [
        k for k in forbidden_at_any_depth(record) if "outcome_wordings" not in k and "decision_rule" not in k
    ]
    if bad:
        raise SystemExit(f"a registration must hold no figure the run would fill: {bad}")
    record["result_manifest"] = result_manifest([Path("data/results") / f"{name}.json" for name in INPUTS])
    path = save_result(NAME, record)
    print(f"wrote {path}, sha256 {hashlib.sha256((ROOT / path).read_bytes()).hexdigest()}")
    for label in ("without_human_for_A", "with_human_for_B"):
        t = record["tree_as_built"][label]
        print(f"tree {label}: {t['tips']} tips, {t['nodes']} nodes, topology {t['topology_sha256']}")
    print(f"frozen code: {json.dumps(record['frozen_code']['blobs'])}")
    print(f"readings registered: {sorted(READINGS)}")
    print(f"outcomes: A {sorted(sn.A_OUTCOMES)}, B {sorted(sn.B_OUTCOMES)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
