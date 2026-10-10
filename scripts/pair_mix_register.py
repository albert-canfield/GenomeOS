#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The mixed calibration's registration: one side a real library, the other random, before any is run.

    uv run python scripts/pair_mix_register.py

Written after the code is frozen in a commit and before any mixed pair is computed. It names §23's
reading and the caveat it is checked against, verbatim from GRAMMAR-BY-COMPARISON (and refuses to run
if a quote is no longer in the document word for word), the mixed pairs and their draw, the null they
are tested under (§23's, by its registered parameters), the two ceilings, the decision rule, the exact
wording of every outcome, what each outcome does to §23's 154 pairs, the frozen code by blob, the tree
by topology digest and every input by sha256. It holds no figure the run would fill: FORBIDDEN is
checked at any depth before the file is written.

`scripts/pair_mix_run.py` refuses to start unless this file is committed, the code it pins is the code
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
from genomeos.knowledge import pair_mix as pm  # noqa: E402
from genomeos.knowledge import pair_null as pn  # noqa: E402
from genomeos.knowledge import species_null as sn  # noqa: E402
from genomeos.results import load_result, save_result  # noqa: E402

NAME = "pair_mix_registration"
ENTRY = "scripts/pair_mix_register.py"
OWN_CODE = (
    "genomeos/knowledge/pair_mix.py",
    "scripts/pair_mix_register.py",
    "scripts/pair_mix_run.py",
    "tests/test_pair_mix.py",
)
#: The code the run must be: pinned by blob, so an edit after registration stops the run. pair_null and
#: species_null are pinned too: the null, its banks and fits, the statistic and the tree are theirs.
FROZEN = (
    "genomeos/knowledge/pair_mix.py",
    "genomeos/knowledge/pair_null.py",
    "genomeos/knowledge/species_null.py",
    "scripts/pair_mix_run.py",
)
DOC = "docs/GRAMMAR-BY-COMPARISON.md"
INPUTS = {
    "origin_presence_genome_wide": "the presence matrix: one bit per gene per species",
    "origin_genome_wide": "each gene's origin stratum, which the age expectation, the matched draws and "
    "the binomial scale read",
    "library_members": "the computed membership of the libraries, as §14, §22 and §23 read it",
    "profiling_genome_wide": "§14's stored pairs, read for the reproduction gate and for which 35 §14 "
    "called; no figure of it is an input to the null",
    "pair_null_registration": "§23's registration, read for the null's registered parameters, which the "
    "registration and the run check equal the code's",
    "pair_null_result": "§23's read: every library pair's stored count of null draws at or above its r "
    "and q, for the reproduction gate, and §23's figures (154 survivors, 6 of §14's 35, 29 lost) that "
    "the outcome wording names; no figure of it is an input to the null",
}
#: §23's registered parameters that define the null being calibrated; each must equal pair_null's.
NULL_PARAMETERS = (
    "P_SIMS",
    "P_SEED",
    "LAMBDA_GRID",
    "SCALE_EPS",
    "P_ALPHA",
    "Q_ALPHA",
    "CAL_MAX_SHARE",
    "CAL_MAX_Q_SHARE",
)

READINGS = {
    "section_23": {
        "section": "GRAMMAR-BY-COMPARISON §23 (2026-10-10), the tip-noise pair null",
        "result_file": "data/results/pair_null_result.json (outcome some_survive; library_pairs, 860 "
        "tested; figures.n_surv 154)",
        "registered_reading_quoted": [
            "154 of 860 testable library pairs are more alike than age and shared ancestry predict under "
            "Pagel's lambda on binomially scaled deviations along the Ensembl Compara species tree "
            "(Grafen branch lengths; lambda, rate and root fitted by REML for each side; "
            "Benjamini-Hochberg q at or below 0.05): 6 of the 35 pairs §14 called more alike than age, "
            "and 148 that §14 did not call. The other 29 of §14's 35 do not survive this null.",
        ],
        "caveat_quoted": [
            "The calibration used random gene sets, whose deviations are small (library sides fit lambda "
            "0.84–0.90, random sides 0.65–0.77), so the null was not checked on two large library "
            "deviations; a mixed calibration, each library against random sets of its partner's "
            "make-up, comes before the 154 are read further.",
        ],
    },
    "what_is_already_known": [
        "§23's read: every library pair's p and q under the null, its 154 pairs at q at or below 0.05, "
        "and the quartiles of the library sides' fitted lambda (0.844, 0.876, 0.898). A library side's "
        "fit depends on its own deviation only, so the library sides of the mixed pairs fit exactly "
        "these values; the random sides' fits are not known",
        "before this registration, three library pairs (testable pairs 1, 401 and 860, all published in "
        "§23's result) were recomputed under the null to check that the computation is deterministic: "
        "each count of null draws at or above r and each fitted lambda equalled §23's stored value. "
        "The stratum pools were read for their sizes (the Homo, Hominoidea and Haplorrhini strata hold "
        "one gene each, Primates three), which is why a draw can fall short. No mixed pair has been "
        "computed, and no random gene set of the mixed draw has been drawn",
    ],
    "what_this_does_not_check": "a pair of two real libraries whose deviations are independent cannot "
    "be built from the data, because whether two real libraries' deviations are independent is the "
    "question §23 asked; the mixed pairs check the null with one library-sized deviation and one random "
    "one, which is the case §23's random pairs did not reach",
}

MIX = {
    "null_calibrated": {
        "model": pm.MODEL,
        "name": pn.MODEL_NAMES[pm.MODEL],
        "what": "the null §23 selected and read the 154 under: Pagel's lambda along §22's tree on each "
        "side's deviation from its age expectation, scaled by the binomial sd of the side's "
        "origin-stratum make-up, lambda, rate and root fitted by REML per side; the same banks (P_SEED, "
        "P_SIMS draws, side one on bank one, side two on bank two) and the same code "
        "(pair_null.pair_test)",
    },
    "mixed_pairs": {
        "what": "for every testable library pair (a, b), with the exclusive members §14 tested: (1) "
        "first_library_real, a's exclusive members on side one against a random gene set on side two "
        "drawn to the origin-stratum make-up of b's exclusive members; (2) second_library_real, a random "
        "gene set on side one drawn to the make-up of a's exclusive members against b's exclusive "
        "members on side two. Each library keeps the side and bank it had in §23's read, and a random "
        "set's age expectation is its partner library's",
        "draw": f"seed {pm.MIX_SEED}, one numpy generator, testable pairs in §14's order; per pair the "
        "random set beside the first library is drawn first, then the one beside the second. Each is "
        "§14's draw (per stratum, strata in sorted order, without replacement) from each stratum's "
        "genes LESS the library side's members, so the two sides share no gene, as a library pair's "
        "exclusive sides share none; where the reduced pool is smaller than the make-up asks, the draw "
        "takes what is left and the shortfall is recorded",
        "why_one_draw_per_mixed_pair": "as §23's calibration drew one random pair per testable pair: "
        "the calibration asks how often the null calls a mixed pair alike across 860 pairs, not "
        "whether any one pair is",
        "count": "two directions, each as large as the number of testable pairs (860 when §14's pairs "
        "reproduce)",
    },
    "statistic": "§14's own, as §23 tested it: Pearson r over the 355 species between the two sides' "
    "residual profiles, one-sided, more alike; p = (null draws with r at or above the observed + 1) / "
    "(draws + 1)",
    "ceilings": {
        "rule": f"in EACH direction, at most {pn.CAL_MAX_SHARE:.0%} of the mixed pairs have p at or "
        f"below {pn.P_ALPHA} AND at most {pn.CAL_MAX_Q_SHARE:.0%} have a Benjamini-Hochberg q, over that "
        f"direction's mixed pairs, at or below {pn.Q_ALPHA}",
        "why_these": "§23's two ceilings, unchanged, so the null is held to the bar it passed on random "
        "pairs; each direction is a set of 860 as each of §23's two calibration sets was, and both must "
        "pass as both of §23's had to",
        "pooled": "the 1,720 mixed pairs pooled, with Benjamini-Hochberg over all of them, are reported "
        "and decide nothing",
    },
    "reproduction_gate": "before any mixed pair is drawn: §14's r_exclusive for every testable pair "
    "(pair_null.reproduce, tolerance 0.0001), and every library pair's count of null draws at or above "
    "its r under the null, which must equal §23's stored count exactly (the computation is "
    "deterministic), and the pairs at q at or below 0.05 recomputed, which must be §23's",
    "descriptive_only": "the pooled figure, the fitted lambda (n, min, quartiles, max) of library and "
    "random sides of the mixed pairs, the shortfall of the random draws, each library's count of its "
    "mixed pairs at p at or below 0.05, and how many of §23's 154 have one or both mixed pairs at p at "
    "or below 0.05 (one random draw each: too noisy to judge any one pair); none decides",
    "parameters": {
        "MODEL": pm.MODEL,
        "MIX_SEED": pm.MIX_SEED,
        "DIRECTIONS": list(pm.DIRECTIONS),
        "P_SIMS": pn.P_SIMS,
        "P_SEED": pn.P_SEED,
        "LAMBDA_GRID": pn.LAMBDA_GRID,
        "SCALE_EPS": pn.SCALE_EPS,
        "P_ALPHA": pn.P_ALPHA,
        "Q_ALPHA": pn.Q_ALPHA,
        "CAL_MAX_SHARE": pn.CAL_MAX_SHARE,
        "CAL_MAX_Q_SHARE": pn.CAL_MAX_Q_SHARE,
        "MIN_MEMBERS": sn.MIN_MEMBERS,
    },
    "decision_rule": pm.MIX_RULE,
    "outcome_wordings": pm.MIX_OUTCOMES,
    "pooled_sentence": pm.POOLED_SENTENCE,
    "what_happens_to_the_154": {
        "calibrated": "§23's 154 of 860 stand as read; no re-reading is needed, and none is made",
        "not_calibrated": "§23's 154 of 860 are withdrawn to no reading, with §23's finding that 29 of "
        "§14's 35 do not survive the null; §14's 35 stand as written, with their caveat",
        "not_reproduced": "nothing is decided about §23; its reading stands as registered, with its "
        "caveat, and the mixed calibration remains to be done",
    },
}

#: Keys no registration may carry at any depth: each is a figure only the run produces.
FORBIDDEN = (
    "p_pair",
    "q_pair",
    "passes",
    "at_or_below_005",
    "q_at_or_below_005",
    "share_at_or_below_005",
    "lambda_random_sides",
    "lambda_library_sides",
    "first_library_real_result",
    "second_library_real_result",
    "pooled_result",
    "outcome",
    "verdict",
    "reading",
    "rows",
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
    s23 = READINGS["section_23"]
    for q in s23["registered_reading_quoted"] + s23["caveat_quoted"]:
        if q not in flat:
            raise SystemExit(f"a quote is not in {DOC} word for word: {q[:80]}...")


def null_parameters_match(registered_23: dict[str, Any]) -> dict[str, Any]:
    """§23's registered parameters against pair_null's constants; refuses on any difference."""
    reg = registered_23["null"]["parameters"]
    out = {}
    for k in NULL_PARAMETERS:
        code = getattr(pn, k)
        if reg.get(k) != code:
            raise SystemExit(f"§23 registered {k} = {reg.get(k)!r}, the code holds {code!r}")
        out[k] = code
    if reg["MODELS"][0] != pm.MODEL:
        raise SystemExit("the null calibrated here is not the first of §23's registered order")
    return {"equal_to_section_23_registration": True, "checked": out}


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
        "same_as": "pair_null_registration.json tree_as_built",
    }


def registration(tree: dict[str, Any], doc_blob: str, params: dict[str, Any]) -> dict[str, Any]:
    return {
        "lane": "lane-pairmix",
        "item": "area J, the mixed calibration of the pair null: the step GRAMMAR-BY-COMPARISON §23 named "
        "before its 154 pairs are read further",
        "registered_in_its_own_commit": "Written after the code was frozen in a commit and before any "
        "mixed pair was drawn or computed. The run reads once, against this file, in the order of the "
        "decision rule, and stops where rule 1 fails.",
        "question": "Does the null §23 read its 154 pairs under find a real library alike to random gene "
        "sets of its partner's origin-stratum make-up no more often than its two ceilings allow, when "
        "one side carries a library-sized deviation and the other a random one?",
        "readings_checked": READINGS,
        "document": {"path": DOC, "blob_at_registration": doc_blob},
        "mixed_calibration": MIX,
        "null_parameters": params,
        "tree_as_built": tree,
        "gates": [
            "this file committed and unmodified",
            "the files in frozen_code at the blobs pinned there",
            "the two tree files at the sha256 species_null pins, and the pruned topology at the digest in "
            "tree_as_built",
            "every input result at the sha256 recorded in result_manifest.inputs",
            "rule 1 (reproduction of §14's statistic and §23's read) before any mixed pair is drawn",
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
        "parameters": MIX["parameters"],
        "exclusions": [
            "every figure the mixed calibration would produce: none is computed by this writer",
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
    reg23 = load_result("pair_null_registration")
    if not reg23:
        raise SystemExit("no pair_null_registration result")
    params = null_parameters_match(reg23)
    record = registration(tree_record(presence["species"]), git("rev-parse", f"HEAD:{DOC}"), params)
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
    print(f"outcomes: {sorted(pm.MIX_OUTCOMES)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
