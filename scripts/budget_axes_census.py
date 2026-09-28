# SPDX-License-Identifier: AGPL-3.0-or-later
"""R7 follow-up: who reads the block classifier's and the budget's output, and what they say now.

    uv run python scripts/budget_axes_census.py [--no-save]

Two findings of review item R7 are left open by lane-ontology (329c2f9). `genome/unknown.py`
`classify_block` returns `regulatory` before it reads RepeatMasker, so a regulatory block never had
its repeat coverage recorded; and `attribution/budget.py` writes role strings that read an absence of
constraint as an absence of function ("best guess neutral", "fossil", "dead pseudogene"). This
script changes nothing. It lists, by hand and checked with `git grep`, every reader of the two
outputs and the field each reads, and counts over the stored `unknown_<chrom>` and `budget_<chrom>`
results: blocks per sequence class and whether their repeat coverage was read, what RepeatMasker
says over the blocks it was not read for, and every budget (tier, label) pair with its selection
reading. No requests, no model.
"""

from __future__ import annotations

import argparse
from collections import Counter
from typing import Any

from genomeos import manifest as mf
from genomeos.results import RESULTS_DIR, load_result, save_result

CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]

# Read by hand at d193afe and checked against `git grep` for every consumer.
CONSUMERS: list[dict[str, str]] = [
    {
        "output": "classify_block (class, evidence, confidence, features, hits)",
        "reader": "genome/unknown.py investigate -> unknown_<chrom>.json; cli.py `genomeos unknown`; "
        "scripts/unknown_genome_wide.py",
        "reads": "all five; by_class tallies the class",
    },
    {
        "output": "unknown_<chrom> blocks[].class",
        "reader": "attribution/budget.py build (the tier rule), attribution/compress.py, "
        "attribution/lexicon.py, "
        "attribution/unknown_scoring.py, genome/blocks.py, web/server.py (unknown panel)",
        "reads": "class only; none reads a class value for origin except compile.py",
    },
    {
        "output": "unknown_<chrom> blocks[].features (repeat_coverage, interspersed_coverage)",
        "reader": "attribution/compile.py region_axes (origin axis; falls back to rmsk_<chrom> when the "
        "features lack interspersed_coverage, which is every regulatory block); genome/blocks.py "
        "(f_* columns of classified blocks)",
        "reads": "interspersed_coverage and repeat_coverage",
    },
    {
        "output": "budget guess.tier",
        "reader": "attribution/organise.py (imports TIERS), compress.py, human_panel.py, lexicon.py, "
        "lexicon_axes.py contexts unknown_fossil/unknown_neutral, variation.py, syntax_tiling.py, "
        "element_types.py, unknown_scoring.py, measurability.py, panel_background.py, compile.py "
        "(region role head), cli.py `genomeos budget`, web/server.py budget_wide and index.html colours, "
        "scripts/epigenome.py, scripts/unknown_coverage.py, scripts/syntax_blocks_matched.py",
        "reads": "the tier string as a join key; every one reads budget_<chrom>.json",
    },
    {
        "output": "budget guess.label",
        "reader": "attribution/compile.py _region only (region `role: <tier>, <label>`); "
        "tests/test_attribution.py, tests/test_attribution_certainty.py, tests/test_ontology.py fixture",
        "reads": "the label text, emitted verbatim",
    },
    {
        "output": "budget guess.confidence and certainty",
        "reader": "budget.tallies guessed_fraction, organise.py, compile.py region `confidence:`, "
        "confidence_census.py",
        "reads": "the hand-set score",
    },
    {
        "output": "the loci benchmark",
        "reader": "genomeos/benchmark/loci.py, loci_miss.py and scripts/loci_*.py",
        "reads": "neither unknown_<chrom> nor budget_<chrom>; no figure of it can move",
    },
]

# Budget labels as stored, with what they read into the constraint.
READINGS = {
    "assembly gap: no sequence to attribute": "no constraint reading",
    "transposable-element fossil, unconstrained": "absence of constraint read as a dead fossil",
    "transposable-element fossil, weak constraint": "weak constraint read as a fossil",
    "repeat-derived sequence under constraint: possibly exapted as a regulatory element": (
        "constraint read as a regulatory mechanism (hedged)"
    ),
    "regulatory elements under constraint; target gene unassigned": "constraint read as selection",
    "regulatory elements by the registry, little constraint: lineage-specific or weak": (
        "little constraint read as weak function"
    ),
    "long open reading frame under constraint: unannotated coding or a young pseudogene": (
        "constraint read as coding (two alternatives)"
    ),
    "long open reading frame without constraint: a chance frame or a dead pseudogene": (
        "absence of constraint read as dead"
    ),
    "constrained non-coding sequence, function unknown": "constraint read as selection",
    "unconstrained unique sequence: no evidence of function, best guess neutral": (
        "absence of constraint read as neutral"
    ),
    "weakly constrained unique sequence: mostly neutral": "weak constraint read as neutral",
}


# Registered 2026-09-28 after the census above (b6e6d17) and before the build.
REGISTERED: dict[str, Any] = {
    "classifier": {
        "change": "classify_block reads RepeatMasker right after the gap test, before the telomere and "
        "regulatory returns, and records repeat_coverage, interspersed_coverage and origin "
        "(repeat_derived/<top class> at >= 0.5 interspersed, partly_repeat_derived/<top> above 0, unique "
        "at 0) whenever RepeatMasker was consulted, an empty coverage included; the rules that name a class "
        "are unchanged and keep their order, so a regulatory block is also repeat_derived when it is",
        "must_not_move": "the class, evidence and confidence of every block; stored unknown_<chrom> files",
        "expected": {"regulatory blocks mostly repeat": 8709, "partly": 6275, "unique": 268},
    },
    "budget": {
        "change": "new writes go to budget_axes_<chrom> and budget_axes_genome_wide with a manifest; each "
        "guess states evidence_status (the constraint reading as selection only), origin, the legacy tier "
        "beside the new one, and a label that reads no function into missing constraint",
        "tier_renamed": {"fossil": "repeat_unconstrained", "neutral": "unconstrained_unknown"},
        "labels_renamed": {
            "transposable-element fossil, unconstrained": (
                "repeat-derived sequence; no sign of selection, which is not a sign of no function"
            ),
            "transposable-element fossil, weak constraint": "repeat-derived sequence; weak sign of selection",
            "repeat-derived sequence under constraint: possibly exapted as a regulatory element": (
                "repeat-derived sequence under selection; role unknown"
            ),
            "regulatory elements by the registry, little constraint: lineage-specific or weak": (
                "regulatory elements by the registry; selection weak or not detected; target gene unassigned"
            ),
            "long open reading frame under constraint: unannotated coding or a young pseudogene": (
                "long open reading frame under selection; coding candidate, role unknown"
            ),
            "long open reading frame without constraint: a chance frame or a dead pseudogene": (
                "long open reading frame; selection weak or not detected; role unknown"
            ),
            "unconstrained unique sequence: no evidence of function, best guess neutral": (
                "unique non-coding sequence; no sign of selection, which is not a sign of no function"
            ),
            "weakly constrained unique sequence: mostly neutral": (
                "unique non-coding sequence; weak sign of selection; role unknown"
            ),
        },
        "must_not_move": [
            "26,806 blocks; per tier under the rename: regulatory 15,536, repeat_unconstrained 7,302, "
            "unconstrained_unknown 2,632, constrained_unknown 1,098, structural 238, and their bp",
            "every block's confidence, so guessed_fraction; measured_bp, constrained_bp",
            "labels map one to one, so the count on each label is the count on the old one",
            "budget_<chrom>, budget_genome_wide and unknown_<chrom> byte for byte; every consumer, since all "
            "read budget_<chrom>: no organise, compile, panel, lexicon, web or loci figure moves",
        ],
        "expected": {
            "evidence_status": {
                "selection_not_detected": 21283,
                "selection_weak": 2962,
                "under_selection": 2410,
                "selection_not_measured": 33,
                "none (gap)": 118,
            },
            "new labels saying neutral, fossil, dead, best guess, lineage-specific or exapted": 0,
            "regulatory-tier blocks with origin repeat_derived": 8709,
            "origin equal to compile.py region_axes origin (first group), over all blocks": "all",
        },
        "negatives_registered": [
            "16 blocks without a constraint reading keep the constrained_unknown tier (a tier the rule "
            "gives when phyloP is missing); not moved here",
            "the tier keys fossil and neutral stay in budget_<chrom> and in every consumer until each is "
            "moved to budget_axes_<chrom>",
        ],
    },
}


def selection(fa: float | None) -> str:
    from genomeos.attribution.budget import CONSTRAINED_MIN, NEUTRAL_MAX

    if fa is None:
        return "selection_not_measured"
    if fa >= CONSTRAINED_MIN:
        return "under_selection"
    return "selection_weak" if fa >= NEUTRAL_MAX else "selection_not_detected"


def reading_of(label: str) -> str:
    if label in READINGS:
        return READINGS[label]
    if label.endswith("constraint not measured"):
        return "constraint not measured"
    if "mechanical role" in label:
        return "no constraint reading (sequence class)"
    return "unlisted"


def chrom_census(chrom: str) -> dict[str, Counter]:
    from genomeos.genome.repeats import INTERSPERSED, repeat_index

    unk = {b["start"]: b for b in (load_result(f"unknown_{chrom}") or {}).get("blocks", [])}
    bud = (load_result(f"budget_{chrom}") or {}).get("blocks", [])
    rindex = repeat_index(chrom)
    out: dict[str, Counter] = {
        k: Counter() for k in ("class_read", "unread_origin", "tier_label", "label_selection", "bp")
    }
    for b in bud:
        f = (unk.get(b["start"]) or {}).get("features") or {}
        read = "read" if f.get("interspersed_coverage") is not None else "not read"
        out["class_read"][(b["class"], read)] += 1
        if read == "not read" and b["class"] != "gap":
            cov = rindex.coverage(b["start"], b["end"]) if rindex else {}
            inter = sum(v for k, v in cov.items() if k in INTERSPERSED) / max(1, b["length"])
            state = "repeat_derived" if inter >= 0.5 else "partly_repeat_derived" if inter > 0 else "unique"
            out["unread_origin"][(b["class"], state)] += 1
        g = b["guess"]
        fa = (b.get("phylop") or {}).get("fraction_above")
        out["tier_label"][(g["tier"], g["label"])] += 1
        out["label_selection"][(g["label"], selection(fa) if b["class"] != "gap" else "n/a")] += 1
        out["bp"][g["tier"]] += b["length"]
    return out


def _flat(c: Counter) -> dict[str, int]:
    return {(" | ".join(k) if isinstance(k, tuple) else k): v for k, v in c.most_common()}


def manifest(names: tuple[str, ...] = ("budget", "unknown")) -> dict[str, Any]:
    inputs = [
        mf.input_entry(p, partition=None)
        for ch in CHROMS
        for p in [RESULTS_DIR / f"{n}_{ch}.json" for n in names] + [RESULTS_DIR / f"rmsk_{ch}.bed.gz"]
        if p.exists()
    ]
    return {
        "sources": [
            {"accession": "this repository, stored budget_<chrom> and unknown_<chrom>", "version": "sha256"},
            {"accession": "UCSC rmsk track, hg38", "version": "rmsk_<chrom> as distilled"},
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {"repeat_derived_min_interspersed_fraction": 0.5, "chromosomes": CHROMS},
        "exclusions": [],
        "partitions": "n/a: a census of labels, no evaluation",
    }


def before() -> dict[str, Any]:
    tot: dict[str, Counter] = {}
    for ch in CHROMS:
        for k, v in chrom_census(ch).items():
            tot.setdefault(k, Counter()).update(v)
        print(ch, "done", flush=True)
    by_reading: Counter = Counter()
    for (label, _sel), n in tot["label_selection"].items():
        by_reading[reading_of(label)] += n
    unread = tot["unread_origin"]
    return {
        "review_item": "R7 follow-up (classifier order, budget role strings)",
        "consumers": CONSUMERS,
        "blocks": sum(tot["tier_label"].values()),
        "blocks_by_class_and_repeat_coverage_read": _flat(tot["class_read"]),
        "unread_blocks_by_class_and_rmsk_origin": _flat(unread),
        "regulatory_or_promoter_like_mostly_repeat_unread": sum(
            n
            for (c, o), n in unread.items()
            if c in ("regulatory", "promoter_like") and o == "repeat_derived"
        ),
        "budget_tier_label": _flat(tot["tier_label"]),
        "budget_label_by_selection": _flat(tot["label_selection"]),
        "budget_blocks_by_constraint_reading": _flat(by_reading),
        "budget_bp_by_tier": dict(tot["bp"]),
        "registered": REGISTERED,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args()
    payload = before()
    name = "budget_axes_census"
    for k, v in payload.items():
        if k != "consumers":
            print(k, v if not isinstance(v, dict) else dict(list(v.items())[:12]))
    if not args.no_save:
        print(f"saved {save_result(name, payload, manifest=manifest())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
