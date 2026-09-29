# SPDX-License-Identifier: AGPL-3.0-or-later
"""Review item R4b census: every stated `confidence` produced in genomeos/attribution/ outside the
compiler, the measured layer and the CRISPRi module (compile.py, measured.py, crispri.py are the
compile half of R4b and wait for R1).

Each entry names the site, what the number is, how it is computed, whether it rises with the size of
an effect, and every consumer found by search (stored results, compiled BioLang, CLI, web, tests).
`kind` is one of:

- `magnitude`: the number is a function of how far a measured quantity lies past a bar, so it rises
  with an effect's size; this is the pattern review R4 forbids;
- `evidence_quality`: a hand-set constant per rule or per evidence layer present, not a function of
  any effect's size; a label for how good the support is, never a probability;
- `no_call`: a constant 0.0 standing for "nothing was called".

Census taken 2026-09-28 at 19423bd. None of these numbers was ever calibrated against an outcome.

Replaced 2026-09-28: human_panel.block_class and variation.case_of return a `certainty` record
(genomeos.certainty.Certainty: effect in its unit, Poisson spread where there is a count, probability
None with its reason) and no `confidence`. budget.guess keeps `confidence` because compile.py emits it,
now documented as an evidence-quality score and repeated as the record's named model score beside the
constrained fraction as effect. candidates.reading is unchanged: no effect size enters it. Stored
results are not rewritten; the change reaches them on their next run.
"""

from __future__ import annotations

CENSUS: tuple[dict[str, object], ...] = (
    {
        "site": "genomeos/attribution/human_panel.py block_class, lineage_restricted",
        "code": "0.3 + 0.3 * min(1.0, (LINEAGE_MAX - presence) / LINEAGE_MAX)",
        "kind": "magnitude",
        "what": (
            "rises as the block's panel presence falls below 0.5: the size of the absence read as certainty"
        ),
        "consumers": [
            "human_panel.py compact_block -> data/results/human_panel_chr*.json blocks[] and genes[]"
            " .confidence (no reader of the field found)",
            "human_panel.py read_candidates and local_candidates rows carry the whole class dict"
            " (class.confidence, no reader of the field found)",
        ],
    },
    {
        "site": "genomeos/attribution/human_panel.py block_class, polymorphic",
        "code": "0.3 + min(0.3, gap * 3)",
        "kind": "magnitude",
        "what": (
            "rises with the gap below presence 0.95 or above a 0.10 recurring-touch share: the"
            " review's named case"
        ),
        "consumers": ["as for lineage_restricted"],
    },
    {
        "site": "genomeos/attribution/human_panel.py block_class, core",
        "code": "min(0.6, 0.3 + 0.3 * min(1, (0.5 - ratio) / 0.5 + (1 if p < 1e-4 else 0)))",
        "kind": "magnitude",
        "what": (
            "rises as observed/expected recurring events fall below 0.5, plus a step when the Poisson"
            " p < 1e-4"
        ),
        "consumers": ["as for lineage_restricted"],
    },
    {
        "site": "genomeos/attribution/human_panel.py block_class, variable",
        "code": "min(0.6, 0.3 + 0.3 * min(1, |ratio - 0.5| / 0.5))",
        "kind": "magnitude",
        "what": "rises with the distance of observed/expected from 0.5",
        "consumers": ["as for lineage_restricted"],
    },
    {
        "site": "genomeos/attribution/human_panel.py block_class, unplaced (not aligned, too short)",
        "code": "0.0",
        "kind": "no_call",
        "what": "a constant standing for 'no class was called'",
        "consumers": ["as for lineage_restricted"],
    },
    {
        "site": "genomeos/attribution/variation.py case_of",
        "code": "0.3 + 0.3 * min(dm, dh), dm/dh = distance of each constrained fraction from its bar",
        "kind": "magnitude",
        "what": "rises as both the mammalian and the human constrained fractions move away from their bars",
        "consumers": [
            "data/results/variation_chr*.json blocks[].case.confidence and the VISTA genome-wide run",
            "genomeos/cli.py cmd_variation 'conf' column",
            "genomeos/benchmark/loci.py stores the whole case dict in the loci benchmark",
            "tests/test_variation.py pins 0.6 at (1.0, 1.0) and 0.3 at the bars: constant-preserving",
        ],
    },
    {
        "site": "genomeos/attribution/budget.py guess, all tiers",
        "code": "per-rule constants 1.0/0.8/0.7/0.6/0.5/0.4; structural max(class, 0.5)",
        "kind": "evidence_quality",
        "what": (
            "one hand-set constant per rule that fired; the constrained fraction selects the rule and the"
            " label (unconstrained vs weak constraint), and within one label the number does not move with it"
        ),
        "consumers": [
            "data/results/budget_chr*.json blocks[].guess.confidence",
            "budget.tallies guessed_fraction (bp whose score >= 0.5) -> budget_genome_wide, cli budget text,"
            " scripts/budget_genome_wide.py",
            "genomeos/attribution/compile.py:107 emits it as `confidence:` of every `region U_*` block in the"
            " compiled chromosome programs -> genomeos/evidence.py (evidence CLI, web Evidence explorer)",
            "genomeos/attribution/organise.py -> data/results/organised_chr*.json rows[].confidence",
            "tests/test_attribution.py (unmeasured <= 0.3; guessed_fraction)",
        ],
    },
    {
        "site": "genomeos/attribution/budget.py:60 guess, constraint not measured",
        "code": "round(min(confidence, 0.3), 2)",
        "kind": "evidence_quality",
        "what": (
            "a cap on the sequence classifier's own per-rule constant (genomeos/genome/unknown.py"
            " classify_block: curated 0.9, predicted 0.4-0.7, inferred 0.3-0.5, none 0.0) because no"
            " constraint was read; not a function of any effect"
        ),
        "consumers": ["as for budget.py guess"],
    },
    {
        "site": "genomeos/attribution/candidates.py reading",
        "code": (
            "base per class + REGISTRY 0.05 + READER 0.15 + MEASURED 0.15 (booleans), -0.1 if borrowed,"
            " cap 0.7"
        ),
        "kind": "evidence_quality",
        "what": "a sum of hand-set weights for which evidence layers are present; no effect size enters it",
        "consumers": [
            "data/results/syntax_candidates_genome_wide.json candidates[].reading.confidence and"
            " mean_confidence",
        ],
    },
)

OUT_OF_SCOPE: dict[str, str] = {
    "genomeos/attribution/compile.py": "compile half of R4b, after R1 (lane-split holds it)",
    "genomeos/attribution/measured.py": "held by lane-split",
    "genomeos/attribution/crispri.py": "held by lane-split",
    "genomeos/attribution/confidence_calibration.py": (
        "measures compile.py's stated confidence, produces none"
    ),
    "genomeos/attribution/target_calibration.py": "measures the sweep's effect-as-confidence; produces none",
}
