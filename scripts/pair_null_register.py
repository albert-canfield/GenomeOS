#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The tip-noise pair null's registration: the calibration, and the conditional re-test of §14's pairs.

    uv run python scripts/pair_null_register.py

Written after the code is frozen in a commit and before the calibration or any library pair is
computed under either tip-noise null. It names the readings verbatim from GRAMMAR-BY-COMPARISON (and
refuses to run if a quote is no longer in the document word for word), the null family and its order,
the calibration rule and the confirmation set, the statistic, the multiple-testing rule, the power
check, the decision rule and the exact wording of every outcome, the frozen code by blob, the tree by
topology digest, and every input by sha256. It holds no figure the run would fill: FORBIDDEN is checked
at any depth before the file is written.

`scripts/pair_null_run.py` refuses to start unless this file is committed, the code it pins is the code
on disk, and the tree and inputs still hash to what is recorded here.
"""

from __future__ import annotations

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
from genomeos.knowledge import pair_null as pn  # noqa: E402
from genomeos.knowledge import species_null as sn  # noqa: E402
from genomeos.results import load_result, save_result  # noqa: E402

NAME = "pair_null_registration"
ENTRY = "scripts/pair_null_register.py"
OWN_CODE = (
    "genomeos/knowledge/pair_null.py",
    "scripts/pair_null_register.py",
    "scripts/pair_null_run.py",
    "tests/test_pair_null.py",
)
#: The code the run must be: pinned by blob, so an edit after registration stops the run. species_null
#: is pinned too, because the statistic, the age expectation, the draws and the tree are its functions.
FROZEN = (
    "genomeos/knowledge/pair_null.py",
    "genomeos/knowledge/species_null.py",
    "scripts/pair_null_run.py",
)
DOC = "docs/GRAMMAR-BY-COMPARISON.md"
INPUTS = {
    "origin_presence_genome_wide": "the presence matrix: one bit per gene per species",
    "origin_genome_wide": "each gene's origin stratum, which the age expectation, the matched draws "
    "and the binomial scale read",
    "library_members": "the computed membership of the libraries, as §14 and §22 read it",
    "profiling_genome_wide": "§14's stored pairs, read for the reproduction gate and for which 35 §14 "
    "called; no figure of it is an input to the null",
    "species_null_result": "§22's result, read for one figure only: how many of its random pairs its "
    "Brownian null put at p at or below 0.05 (A_library_pairs.grafen.calibration), set beside the "
    "same figure recomputed here; no figure of it is an input to the null",
}

READINGS = {
    "section_14": {
        "section": "GRAMMAR-BY-COMPARISON §14 (2026-09-12)",
        "result_file": "data/results/profiling_genome_wide.json (null.pairs, 860 tested; "
        "null.significant_q05, 35 rows)",
        "quotes": [
            "**What is more alike than age predicts.** 35 of 860 pairs at a false discovery rate of 5%, "
            "led by the immune and signalling libraries: immune against the ligand-receptor protocol "
            "(observed −0.76, age-matched −0.93, z 9.3), degradation and cytoskeleton against the same "
            "protocol (z 8.5), the cell cycle against immunity (z 6.7), blood with the growth-longevity "
            "axis (observed +0.72 against −0.27).",
        ],
    },
    "section_22": {
        "section": "GRAMMAR-BY-COMPARISON §22 (2026-10-10), the Brownian species-tree null",
        "result_file": "data/results/species_null_result.json (A_library_pairs: outcome null_not_calibrated)",
        "quotes": [
            "No reading: the species-tree null is not calibrated on this data. Of 860 pairs of random "
            "gene sets drawn to the tested pairs' origin-stratum make-up, 205 (23.8%) have a tree-null "
            "p at or below 0.05, above the registered ceiling of 10%, so the null finds random genes "
            "alike more often than it should and a library pair's p under it is not evidence.",
            "In the pair null the simulated deviations swamp the age expectation (null means −0.14 to "
            "+0.04 where the age-expected r is −0.96 to +0.99), which is the likely cause of its "
            "miscalibration; the next pair null adds a tip-noise term and is calibrated on random pairs "
            "before any pair is read.",
        ],
    },
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
    "what_is_already_known": "§22 computed every library pair's p under its Brownian null and "
    "reported nine of them; those p values were read under a null that failed calibration and decide "
    "nothing here. No library pair has been computed under either tip-noise null, and none is until "
    "rules 2 and 3 pass.",
}

NULL = {
    "statistic": "§14's own: Pearson r over the 355 species between two libraries' residual profiles "
    "(share of the members exclusive to each side with an orthologue in the species, minus the share "
    "of all 19,391 genes), r_exclusive, recomputed and required to match the stored value; one-sided, "
    "more alike",
    "decomposition": "each side's profile P = E + D, with E the age expectation (the exact mean "
    "profile of §14's matched draw: per origin stratum, the members' count times the mean profile of "
    "every gene of that stratum, over the members with an origin; species_null.age_expectation) and D "
    "the deviation. The null holds E and replaces D by an independent draw from a model fitted to D",
    "family": {
        "model": "D = h * y, y ~ N(mu 1, sigma^2 V(lambda)), V(lambda) = lambda R + (1 - lambda) I, R "
        "the Brownian covariance of the species tree on Grafen branch lengths (a correlation matrix, "
        "since the Grafen tree is ultrametric with the root at height 1). V is Pagel's lambda and, on "
        "this tree, exactly Brownian motion plus independent white noise at the tips (tip noise)",
        "fit": "per side, on the side's own deviation, by REML on the species where h is positive: "
        "sigma^2 and mu profiled, lambda the maximiser over a grid of 1,001 points on [0, 1] (the "
        "first when two tie). Random sets and libraries are fitted by the same code",
        "simulation": "a side's draw is h * (mu + sigma (sqrt(lambda) b + sqrt(1 - lambda) w)), b a "
        "unit-diagonal Brownian draw along the tree and w unit white noise; side one uses banks one, "
        "side two banks two, each bank P_SIMS draws, shared by every pair (each p exact, different "
        "pairs' p values dependent through the shared draws); p = (draws with r at or above the "
        "observed + 1) / (draws + 1)",
        "members_in_registered_order": {
            "1_lambda_binomial": {
                "name": pn.MODEL_NAMES["lambda_binomial"],
                "h": "the binomial sd of the side's origin-stratum make-up at each species, "
                "sqrt(sum_k n_k p_k (1 - p_k) / sum_k n_k), p_k the share of stratum k's genes with an "
                "orthologue there; zero where every stratum of the make-up is all present or all "
                "absent, and there D is exactly zero for random sets and libraries alike",
                "why_first": "a share cannot move at a species where its strata are fixed (104 species "
                "for a Eukaryota-born set, 341 for a Catarrhini-born one, read from the strata means "
                "before registration) and can move most where they are half present; the scaled model "
                "has the same three parameters as the plain one and draws variance only where a "
                "deviation can occur",
            },
            "2_lambda": {
                "name": pn.MODEL_NAMES["lambda"],
                "h": "1 at every species",
                "why_second": "the null §22 named, without the scale; tried only if the first fails "
                "the first calibration set",
            },
        },
        "brownian_plus_white_noise": "not a separate member: on the ultrametric Grafen tree it is the "
        "same model as Pagel's lambda, with sigma_BM^2 = lambda sigma^2 and tau^2 = (1 - lambda) sigma^2",
    },
    "tree": "the tree §22 built and registered (species_null_registration.json tree_as_built."
    "without_human_for_A): Ensembl Compara's pan tree of Ensembl Genomes r63 with the vertebrate tree "
    "of r116 grafted in, 355 tips; its topology digest is recorded below and checked by the run",
    "branch_lengths": "Grafen, rho 1, the length §22 registered. The equal-length sensitivity of §22 "
    "is not repeated: with a tip-noise term fitted per side, the deep-against-shallow share of "
    "variance is fitted rather than fixed by the lengths",
    "calibration": {
        "what": "for every testable pair, one random pair of gene sets, each set drawn to its side's "
        "origin-stratum make-up without replacement per stratum (§14's draw, species_null.matched_rows), "
        "tested exactly as a library pair. Random pairs only: no library pair is read",
        "first_set": f"seed {pn.CAL1_SEED} in species_null.null_a's draw order: the very random pairs "
        "§22's Brownian null was calibrated on. §22's Brownian null is recomputed on them (its own "
        "banks, seed and draws) and its count at p at or below 0.05 is reported beside §22's; that "
        "figure decides nothing",
        "second_set": f"seed {pn.CAL2_SEED}, the same draw: the confirmation set, computed only for "
        "the selected null",
        "rule": f"a null passes a set when at most {pn.CAL_MAX_SHARE:.0%} of its random pairs have p "
        f"at or below {pn.P_ALPHA} AND at most {pn.CAL_MAX_Q_SHARE:.0%} have a Benjamini-Hochberg q, "
        f"computed over that set's random pairs, at or below {pn.Q_ALPHA}",
        "why_10_percent": "§22's ceiling, kept so the two nulls are judged alike: twice the nominal "
        "rate leaves room for the shared banks and the estimated parameters, and refuses a null like "
        "§22's (23.8%)",
        "why_the_tail_ceiling": "the reading is made at q at or below 0.05 over 860 pairs, which lives "
        "in the null's far tail (p near 0.05 k / 860); a null calibrated at 0.05 but heavy in that tail "
        "would make survivors of random genes. Under a calibrated null with every pair null, any "
        "Benjamini-Hochberg discovery happens at most about 5 times in 100; 1% of the pairs (8 of 860) "
        "leaves room for the dependence of shared draws and refuses a null like §22's, which put 49 "
        "of 860 random pairs at q at or below 0.05",
        "selection": "the first null in the registered order that passes the first set; it must pass "
        "the second set too before any library pair is read. Selection uses random pairs only",
        "failure": "neither passes the first set: no reading (null_not_calibrated). The selected null "
        "fails the second: no reading (not_confirmed). Either is recorded as the result, and no "
        "library pair is computed under the null",
    },
    "multiple_testing": "Benjamini-Hochberg over every testable library pair, q at or below 0.05, as "
    "§14 and §22 did",
    "power": {
        "what": "each library split at random into two halves (§22's split, its seed "
        f"{pn.HALVES_SEED}), tested as a pair under the selected null: the halves share the library's "
        "history, so a null that never sees them as alike cannot see a pair",
        "floor": f"fewer than {pn.HALVES_MIN_SHARE:.0%} of libraries at p at or below 0.05 and a "
        "result with no survivor is read as low power",
        "when": "computed only after rule 3 passes",
    },
    "falsifier": "§14's reading that a pair is more alike than age predicts is not distinguished from "
    "age plus this model of shared ancestry, for every one of its 35 pairs whose q under the selected "
    "null is above 0.05, provided the null passed both calibration sets",
    "descriptive_only": "the fitted lambda and rate of every side, the null's mean and sd per pair, "
    "the lambda quartiles of random and library sides, §22's Brownian null recomputed on the first "
    "calibration set, and the calibration figures of the null not selected; none decides",
    "parameters": {
        "P_SIMS": pn.P_SIMS,
        "P_SEED": pn.P_SEED,
        "CAL1_SEED": pn.CAL1_SEED,
        "CAL2_SEED": pn.CAL2_SEED,
        "HALVES_SEED": pn.HALVES_SEED,
        "LAMBDA_GRID": pn.LAMBDA_GRID,
        "P_ALPHA": pn.P_ALPHA,
        "Q_ALPHA": pn.Q_ALPHA,
        "CAL_MAX_SHARE": pn.CAL_MAX_SHARE,
        "CAL_MAX_Q_SHARE": pn.CAL_MAX_Q_SHARE,
        "HALVES_MIN_SHARE": pn.HALVES_MIN_SHARE,
        "MIN_MEMBERS": sn.MIN_MEMBERS,
        "SCALE_EPS": pn.SCALE_EPS,
        "MODELS": list(pn.MODELS),
        "brownian_recomputed_with": {"A_SIMS": sn.A_SIMS, "A_SEED": sn.A_SEED},
    },
    "decision_rule": pn.P_RULE,
    "outcome_wordings": pn.P_OUTCOMES,
    "calibration_sentence": pn.CAL_SENTENCE,
}

#: Keys no registration may carry at any depth: each is a figure only the run produces.
FORBIDDEN = (
    "p_pair",
    "q_pair",
    "p_tree",
    "q_tree",
    "passes",
    "at_or_below_005",
    "q_at_or_below_005",
    "survivors",
    "n_surv",
    "k_14",
    "selected",
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
    for key in ("section_14", "section_22"):
        for q in READINGS[key]["quotes"]:
            if q not in flat:
                raise SystemExit(f"{key}: a quote is not in {DOC} word for word: {q[:80]}...")


def tree_record(species: list[str]) -> dict[str, Any]:
    root, _ = sn.species_tree(species, with_human=False)
    nodes = sn.preorder(root)
    t = sn.arrays(root, species)
    cov = sn.covariance(t, sn.grafen_lengths(t))
    return {
        "tips": len(sn.leaves(root)),
        "nodes": len(nodes),
        "topology_sha256": sn.topology_digest(root),
        "grafen_diagonal_is_one": bool(abs(cov.diagonal() - 1).max() < 1e-12),
        "same_as": "species_null_registration.json tree_as_built.without_human_for_A",
    }


def registration(tree: dict[str, Any], doc_blob: str) -> dict[str, Any]:
    return {
        "lane": "lane-pairnull",
        "item": "area J, the calibrated pair null: the next step GRAMMAR-BY-COMPARISON §22 named",
        "registered_in_its_own_commit": "Written after the code was frozen in a commit and before the "
        "calibration or any library pair was computed under either tip-noise null. The run reads once, "
        "against this file, in the order of the decision rule, and stops where a rule fails.",
        "question": "Is there a pair null with a tip-noise term that finds random age-matched gene sets "
        "alike no more often than it should, and if there is, which of §14's 35 of 860 library pairs "
        "stay more alike than age and this model of shared ancestry predict?",
        "readings_retested": READINGS,
        "document": {"path": DOC, "blob_at_registration": doc_blob},
        "null": NULL,
        "tree_as_built": tree,
        "gates": [
            "this file committed and unmodified",
            "the files in frozen_code at the blobs pinned there",
            "the two tree files at the sha256 species_null pins, and the pruned topology at the digest "
            "in tree_as_built",
            "every input result at the sha256 recorded in result_manifest.inputs",
            "rule 1 (reproduction) before rule 2, rule 2 (calibration) before rule 3, rule 3 "
            "(confirmation) before any library pair",
        ],
        "frozen_code": frozen_code(),
    }


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
        "assembly": "n/a: species presence of orthologues and a species tree; no genomic coordinate",
        "coordinates": "n/a: no interval is computed or read",
        "parameters": NULL["parameters"],
        "exclusions": [
            "every figure the null would produce: none is computed by this writer",
            "the species list is read from the presence matrix to build the tree; no profile is read, "
            "and the other results are hashed so the run can refuse inputs that moved",
        ],
        "partitions": "n/a: a registration written before the measurement; nothing is evaluated",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    mf.trace_begin()
    for name in (sn.PAN_TREE.name, sn.VERT_TREE.name):
        path = sn.TREE_DIR / name
        if not path.exists():
            raise SystemExit(f"no {path}: run scripts/species_null_register.py --fetch")
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
    t = record["tree_as_built"]
    print(f"tree: {t['tips']} tips, {t['nodes']} nodes, topology {t['topology_sha256']}")
    print(f"frozen code: {json.dumps(record['frozen_code']['blobs'])}")
    print(f"outcomes: {sorted(pn.P_OUTCOMES)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
