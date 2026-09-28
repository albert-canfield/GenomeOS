# SPDX-License-Identifier: AGPL-3.0-or-later
"""The composition budget: every UNKNOWN block with constraint attached and a best guess.

The house argument. A plan that calls half the house pipes is wrong before a brick is
laid, because the quantities do not fit. The genome has quantities too: a tenth of it
is under purifying selection across mammals and the rest has mutated freely for tens
of millions of years. So a block's sequence class (from `genomeos unknown`) and its
constraint (from Zoonomia and the 100-vertebrate elements) together support a guess
that is far better than "no clue":

    structural             gaps, centromeres, satellite arrays: sequence with a mechanical role
    repeat_unconstrained   interspersed-repeat blocks with under 5% of bases constrained
    regulatory             ENCODE-backed element clusters; the target gene is the open question
    constrained_unknown    under selection, not coding, not regulatory by the registry: the real unknown
    unconstrained_unknown  unique or coding-candidate sequence with under 5% of bases constrained

Review R7 (2026-09-28): constraint is evidence of selection only, and its absence is not evidence of
no function. The tiers above are the names new writes use (`budget_axes_<chrom>`); the stored
`budget_<chrom>` results and every consumer of them keep the legacy keys `fossil` and `neutral`
(`TIERS`, `LEGACY_TIER`), and each new guess states its legacy tier beside the new one. Each guess
also states `evidence_status` (the constraint reading: under_selection, selection_weak,
selection_not_detected, selection_not_measured) and `origin` (the classifier's RepeatMasker reading,
or what the sequence class itself states), so a regulatory block that is mostly LINE says both.

Every guess carries an evidence-quality score (`confidence`, hand-set per rule, not a
probability) and a certainty record with the numbers it rests on. The tiers are a
vocabulary for the next steps (attribution of a target gene and a tissue), not a
verdict; `constrained_unknown` is where the work is.
"""

from __future__ import annotations

import time
from pathlib import Path

from genomeos.attribution.constraint import (
    ELEMENTS_TRACK,
    PHYLOP_241_URL,
    PHYLOP_THRESHOLD,
    elements_over_blocks,
    fetch_elements,
    phylop_over_blocks,
)
from genomeos.certainty import Certainty
from genomeos.results import load_result, save_result

#: the tier keys of the stored budget_<chrom> results, which every consumer joins on
TIERS = ("structural", "fossil", "regulatory", "constrained_unknown", "neutral")
#: the tier names new writes use (review R7): no tier reads missing constraint as no function
AXES_TIERS = (
    "structural",
    "repeat_unconstrained",
    "regulatory",
    "constrained_unknown",
    "unconstrained_unknown",
)
LEGACY_TIER = {"repeat_unconstrained": "fossil", "unconstrained_unknown": "neutral"}
AXES_TIER = {v: k for k, v in LEGACY_TIER.items()}
NEUTRAL_MAX = 0.03  # below this constrained fraction a block reads as unconstrained
CONSTRAINED_MIN = 0.05  # from here on the block carries constraint worth a name
STRUCTURAL = {"centromere", "satellite_array", "tandem_repeat"}
REGULATORY = {"regulatory", "promoter_like"}


SCORE_NAME = (
    "evidence-quality score: the hand-set constant of the budget rule that fired, capped at 0.3 when"
    " constraint is unmeasured; within one label it does not move with the constrained fraction; not"
    " a probability"
)
NO_PROBABILITY = (
    "no calibration record: no set of blocks with a known function outcome has been scored against"
    " these tiers, so no probability of the tier is quoted"
)
CONSTRAINT_UNIT = "fraction of measured bases at Zoonomia phyloP >= 2.27 (241 mammals)"
#: sequence classes whose origin the class itself states (as compile.CLASS_ORIGIN)
CLASS_ORIGIN = {
    "gap": "assembly_gap",
    "centromere": "satellite",
    "satellite_array": "satellite",
    "tandem_repeat": "tandem_repeat",
    "telomere": "tandem_repeat",
}


def selection(fraction: float | None) -> str:
    """The R7 evidence_status of a constrained fraction: selection, never a mechanism or its absence."""
    if fraction is None:
        return "selection_not_measured"
    if fraction >= CONSTRAINED_MIN:
        return "under_selection"
    return "selection_weak" if fraction >= NEUTRAL_MAX else "selection_not_detected"


def block_origin(cls: str, features: dict | None) -> str:
    """Origin as compile.region_axes states it: from the class when the class states it, else from the
    classifier's RepeatMasker reading, else unknown."""
    if cls in CLASS_ORIGIN:
        return CLASS_ORIGIN[cls]
    if cls.startswith("interspersed_repeat"):
        top = cls.removeprefix("interspersed_repeat").lstrip("_")
        return f"repeat_derived/{top}" if top else "repeat_derived"
    return (features or {}).get("origin") or "unknown"


def guess(
    cls: str,
    confidence: float,
    phylop: dict | None,
    elements: dict | None,
    features: dict | None = None,
) -> dict:
    """Class plus constraint to a tier, a label, an evidence-quality score and a certainty record.

    With the R7 axes (2026-09-28): `evidence_status` (constraint as selection only), `origin`
    (from `features`, the classifier's output, when given) and the legacy tier the stored results use.

    Review R4b (2026-09-28): `confidence` stays because compile.py emits it as the `confidence:` of
    each region block, but it is a hand-set evidence-quality score per rule (census:
    genomeos/attribution/confidence_census.py), never a probability. The record beside it states the
    constrained fraction as the effect, in its unit, and that no probability exists.
    """
    g = _rule(cls, confidence, phylop, elements)
    fa = phylop.get("fraction_above") if phylop else None
    source = "curated: assembly gap" if cls == "gap" else f"inferred: sequence class {cls} plus constraint"
    g["certainty"] = Certainty(
        evidence_category=source if fa is not None or cls == "gap" else f"{source} not measured",
        effect_estimate=fa,
        effect_unit=CONSTRAINT_UNIT if fa is not None else "",
        uncertainty_note="a block summary; no spread computed" if fa is not None else "constraint not read",
        model_score=g["confidence"],
        model_score_name=SCORE_NAME,
        probability_unavailable=NO_PROBABILITY,
    ).to_dict()
    g["legacy_tier"] = LEGACY_TIER.get(g["tier"], g["tier"])
    g["evidence_status"] = None if cls == "gap" else selection(fa)
    g["origin"] = block_origin(cls, features)
    return g


def _rule(cls: str, confidence: float, phylop: dict | None, elements: dict | None) -> dict:
    """Class plus constraint to a tier, a label and the rule's evidence-quality score."""
    fa = phylop.get("fraction_above") if phylop else None
    n_el = elements.get("n", 0) if elements else 0
    if cls == "gap":
        return {"tier": "structural", "label": "assembly gap: no sequence to attribute", "confidence": 1.0}
    if cls in STRUCTURAL:
        return {
            "tier": "structural",
            "label": f"{cls.replace('_', ' ')}: repeat array with a mechanical role, not read as genes",
            "confidence": round(max(confidence, 0.5), 2),
        }
    if fa is None:
        return {
            "tier": "constrained_unknown" if cls not in REGULATORY else "regulatory",
            "label": f"{cls.replace('_', ' ')}; constraint not measured",
            "confidence": round(min(confidence, 0.3), 2),
        }
    if cls.startswith("interspersed_repeat"):
        if fa >= CONSTRAINED_MIN:
            return {
                "tier": "constrained_unknown",
                "label": "repeat-derived sequence under selection; role unknown",
                "confidence": 0.5,
            }
        if fa < NEUTRAL_MAX:
            return {
                "tier": "repeat_unconstrained",
                "label": "repeat-derived sequence; no sign of selection, which is not a sign of no function",
                "confidence": 0.8,
            }
        return {
            "tier": "repeat_unconstrained",
            "label": "repeat-derived sequence; weak sign of selection",
            "confidence": 0.6,
        }
    if cls in REGULATORY:
        if fa >= CONSTRAINED_MIN:
            return {
                "tier": "regulatory",
                "label": "regulatory elements under constraint; target gene unassigned",
                "confidence": 0.7,
            }
        return {
            "tier": "regulatory",
            "label": (
                "regulatory elements by the registry; selection weak or not detected; target gene unassigned"
            ),
            "confidence": 0.5,
        }
    if cls == "long_orf":
        if fa >= CONSTRAINED_MIN:
            return {
                "tier": "constrained_unknown",
                "label": "long open reading frame under selection; coding candidate, role unknown",
                "confidence": 0.5,
            }
        return {
            "tier": "unconstrained_unknown",
            "label": "long open reading frame; selection weak or not detected; role unknown",
            "confidence": 0.5,
        }
    if fa >= CONSTRAINED_MIN or (n_el >= 3 and fa >= NEUTRAL_MAX):
        return {
            "tier": "constrained_unknown",
            "label": "constrained non-coding sequence, function unknown",
            "confidence": 0.6,
        }
    if fa < NEUTRAL_MAX:
        return {
            "tier": "unconstrained_unknown",
            "label": "unique non-coding sequence; no sign of selection, which is not a sign of no function",
            "confidence": 0.6,
        }
    return {
        "tier": "unconstrained_unknown",
        "label": "unique non-coding sequence; weak sign of selection; role unknown",
        "confidence": 0.4,
    }


def chromosome_length(chrom: str) -> int | None:
    a = load_result("anatomy_hg38_by_chromosome") or {}
    c = (a.get("chromosomes") or {}).get(chrom)
    return c.get("length") if c else None


def composition(chrom: str) -> dict | None:
    a = load_result("anatomy_hg38_by_chromosome") or {}
    c = (a.get("chromosomes") or {}).get(chrom)
    return c.get("composition_bp") if c else None


def build(
    chrom: str,
    threshold: float = PHYLOP_THRESHOLD,
    phylop: bool = True,
    elements: bool = True,
    progress=None,
    unknown: dict | None = None,
    length: int | None = None,
) -> dict:
    """The budget of one chromosome from its UNKNOWN result plus constraint."""
    unknown = unknown or load_result(f"unknown_{chrom}")
    if not unknown:
        raise FileNotFoundError(f"no unknown_{chrom} result; run genomeos unknown --chrom {chrom}")
    blocks = sorted(unknown["blocks"], key=lambda b: b["start"])
    length = length or chromosome_length(chrom) or max(b["end"] for b in blocks)
    t0 = time.time()
    cost: dict = {}
    # phyloP: every block that holds sequence (gaps are N and have no alignment)
    measured = [i for i, b in enumerate(blocks) if b["class"] != "gap"]
    ph: list[dict | None] = [None] * len(blocks)
    if phylop and measured:
        ivs = [(blocks[i]["start"], blocks[i]["end"]) for i in measured]
        stats, cost["phylop"] = phylop_over_blocks(chrom, ivs, threshold, progress=progress)
        for i, st in zip(measured, stats, strict=True):
            ph[i] = st.as_dict()
    el: list[dict | None] = [None] * len(blocks)
    if elements:
        els = fetch_elements(chrom, length)
        cost["elements"] = {"n": len(els), "track": ELEMENTS_TRACK}
        for i, r in zip(
            measured,
            elements_over_blocks(els, [(blocks[i]["start"], blocks[i]["end"]) for i in measured]),
            strict=True,
        ):
            el[i] = r
    rows = []
    for b, p, e in zip(blocks, ph, el, strict=True):
        g = guess(b["class"], b.get("confidence", 0.0), p, e, b.get("features"))
        rows.append(
            {
                "start": b["start"],
                "end": b["end"],
                "length": b["length"],
                "class": b["class"],
                "class_evidence": b.get("evidence"),
                "class_confidence": b.get("confidence"),
                "phylop": p,
                "elements": e,
                "guess": g,
            }
        )
    out = {
        "chrom": chrom,
        "chromosome_length": length,
        "composition_bp": composition(chrom),
        "threshold": threshold,
        "sources": {
            "phylop": f"Zoonomia 241 placental mammals, {PHYLOP_241_URL}" if phylop else None,
            "elements": f"UCSC hg38 {ELEMENTS_TRACK}" if elements else None,
        },
        "unknown_bp": unknown["unknown_bp"],
        "blocks": rows,
        "cost": {**cost, "seconds": round(time.time() - t0, 1)},
    }
    out.update(tallies(rows, length))
    return out


def tallies(rows: list[dict], length: int | None, tiers: tuple[str, ...] = AXES_TIERS) -> dict:
    by_class: dict[str, dict] = {}
    by_tier: dict[str, dict] = {t: {"blocks": 0, "bp": 0, "constrained_bp": 0} for t in tiers}
    constrained_total = 0
    measured_bp = 0
    for r in rows:
        c = by_class.setdefault(r["class"], {"blocks": 0, "bp": 0, "constrained_bp": 0, "measured_bp": 0})
        c["blocks"] += 1
        c["bp"] += r["length"]
        above = (r["phylop"] or {}).get("above") or 0
        bases = (r["phylop"] or {}).get("bases") or 0
        c["constrained_bp"] += above
        c["measured_bp"] += bases
        t = by_tier.setdefault(r["guess"]["tier"], {"blocks": 0, "bp": 0, "constrained_bp": 0})
        t["blocks"] += 1
        t["bp"] += r["length"]
        t["constrained_bp"] += above
        constrained_total += above
        measured_bp += bases
    total = sum(c["bp"] for c in by_class.values()) or 1
    for c in by_class.values():
        c["constrained_fraction"] = (
            round(c["constrained_bp"] / c["measured_bp"], 4) if c["measured_bp"] else None
        )
    for t in by_tier.values():
        t["fraction_of_unknown"] = round(t["bp"] / total, 4)
        t["fraction_of_chromosome"] = round(t["bp"] / length, 4) if length else None
    guessed = sum(r["length"] for r in rows if r["guess"]["confidence"] >= 0.5)
    return {
        "by_class": dict(sorted(by_class.items(), key=lambda kv: -kv[1]["bp"])),
        "by_tier": by_tier,
        "measured_bp": measured_bp,
        "constrained_bp": constrained_total,
        "constrained_fraction": round(constrained_total / measured_bp, 4) if measured_bp else None,
        "guessed_fraction": round(guessed / total, 4),
    }


def distil(results_dir: Path | None = None, name: str = "budget", tiers: tuple[str, ...] = TIERS) -> dict:
    """Every <name>_chr*.json summed: the genome's composition budget. The default reads the stored
    budget_<chrom> results with their legacy tiers; name="budget_axes", tiers=AXES_TIERS the R7 ones."""
    kw = {"results_dir": results_dir} if results_dir else {}
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    per: dict[str, dict] = {}
    by_class: dict[str, dict] = {}
    by_tier: dict[str, dict] = {t: {"blocks": 0, "bp": 0, "constrained_bp": 0} for t in tiers}
    genome = 0
    unknown_bp = 0
    measured = 0
    constrained = 0
    guessed = 0
    cost_mb = 0.0
    for c in chroms:
        r = load_result(f"{name}_{c}", **kw)
        if not r:
            continue
        per[c] = {
            "unknown_bp": r["unknown_bp"],
            "constrained_fraction": r.get("constrained_fraction"),
            "guessed_fraction": r.get("guessed_fraction"),
            "by_tier_bp": {t: v["bp"] for t, v in r["by_tier"].items()},
        }
        genome += r.get("chromosome_length") or 0
        unknown_bp += r["unknown_bp"]
        measured += r.get("measured_bp") or 0
        constrained += r.get("constrained_bp") or 0
        guessed += round((r.get("guessed_fraction") or 0) * r["unknown_bp"])
        cost_mb += ((r.get("cost") or {}).get("phylop") or {}).get("mb_fetched", 0)
        for k, v in r["by_class"].items():
            d = by_class.setdefault(k, {"blocks": 0, "bp": 0, "constrained_bp": 0, "measured_bp": 0})
            for f in ("blocks", "bp", "constrained_bp", "measured_bp"):
                d[f] += v.get(f) or 0
        for t, v in r["by_tier"].items():
            d = by_tier.setdefault(t, {"blocks": 0, "bp": 0, "constrained_bp": 0})
            for f in ("blocks", "bp", "constrained_bp"):
                d[f] += v.get(f) or 0
    for d in by_class.values():
        d["constrained_fraction"] = (
            round(d["constrained_bp"] / d["measured_bp"], 4) if d["measured_bp"] else None
        )
    for t in by_tier.values():
        t["fraction_of_unknown"] = round(t["bp"] / unknown_bp, 4) if unknown_bp else None
        t["fraction_of_genome"] = round(t["bp"] / genome, 4) if genome else None
    return {
        "chromosomes": len(per),
        "per_chromosome": per,
        "genome_bp": genome,
        "unknown_bp": unknown_bp,
        "measured_bp": measured,
        "constrained_bp": constrained,
        "constrained_fraction": round(constrained / measured, 4) if measured else None,
        "guessed_fraction": round(guessed / unknown_bp, 4) if unknown_bp else None,
        "by_class": dict(sorted(by_class.items(), key=lambda kv: -kv[1]["bp"])),
        "by_tier": by_tier,
        "threshold": PHYLOP_THRESHOLD,
        "phylop_mb_fetched": round(cost_mb, 1),
        "note": (
            "Constraint is Zoonomia phyloP >= 2.27 over 241 mammals (5% FDR), read per base from UCSC and "
            "never stored; elements are the 100-vertebrate phastCons set. A tier is a best guess with its "
            "confidence, not a verdict; constrained_unknown is where attribution work remains."
        ),
    }


AXES_NAME = "budget_axes"


def _manifest(chrom: str, names: list[str], phylop: bool = True, results_dir: Path | None = None) -> dict:
    from genomeos import manifest as mf
    from genomeos.results import RESULTS_DIR

    rd = results_dir or RESULTS_DIR
    paths = [rd / f"{n}_{chrom}.json" for n in names] + [rd / f"rmsk_{chrom}.bed.gz"]
    return {
        "sources": [
            {"accession": f"Zoonomia 241-mammal phyloP, {PHYLOP_241_URL}", "version": "as read by the budget"}
            if phylop
            else {"accession": "none: constraint not read", "version": "n/a"},
            {"accession": f"UCSC hg38 {ELEMENTS_TRACK}", "version": "as fetched"},
            {"accession": "UCSC rmsk track, hg38", "version": "rmsk_<chrom> as distilled"},
        ],
        "inputs": [mf.input_entry(p, partition=None) for p in paths if p.exists()],
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "chrom": chrom,
            "phylop_threshold": PHYLOP_THRESHOLD,
            "constrained_min": CONSTRAINED_MIN,
            "neutral_max": NEUTRAL_MAX,
            "tiers": list(AXES_TIERS),
            "legacy_tier": LEGACY_TIER,
        },
        "exclusions": ["gap blocks carry no constraint reading"],
        "partitions": "n/a: a rule per block, no evaluation",
    }


def run_and_save(chrom: str, **kw) -> dict:
    """Build and save under the R7 name `budget_axes_<chrom>`; stored `budget_<chrom>` keeps its values."""
    out = build(chrom, **kw)
    save_result(f"{AXES_NAME}_{chrom}", out, manifest=_manifest(chrom, ["unknown"], kw.get("phylop", True)))
    return out


def restate(chrom: str, results_dir: Path | None = None, rindex=None) -> dict:
    """The stored budget re-read under R7 without a request: each block's stored constraint and
    elements, the classifier's RepeatMasker origin (unknown.repeat_origin, read now for every block,
    the regulatory ones included), and the relabelled guess. Rows keep what the guess rests on and
    drop the per-row certainty record, whose only varying field is the constrained fraction kept here."""
    from genomeos.genome.unknown import repeat_origin

    kw = {"results_dir": results_dir} if results_dir else {}
    budget = load_result(f"budget_{chrom}", **kw)
    if not budget:
        raise FileNotFoundError(f"no budget_{chrom} result")
    feats = {
        b["start"]: b.get("features") or {}
        for b in (load_result(f"unknown_{chrom}", **kw) or {}).get("blocks", [])
    }
    if rindex is None:
        from genomeos.genome.repeats import repeat_index

        rindex = repeat_index(chrom)
    rows = []
    for b in sorted(budget["blocks"], key=lambda b: b["start"]):
        f = dict(feats.get(b["start"]) or {})
        if b["class"] != "gap" and rindex:
            f.update(repeat_origin(rindex.coverage(b["start"], b["end"]), b["length"]))
        g = guess(b["class"], b.get("class_confidence") or 0.0, b.get("phylop"), b.get("elements"), f)
        ph = b.get("phylop") or {}
        rows.append(
            {
                "start": b["start"],
                "end": b["end"],
                "length": b["length"],
                "class": b["class"],
                "phylop": {k: ph[k] for k in ("bases", "above", "fraction_above") if k in ph} or None,
                "interspersed_coverage": f.get("interspersed_coverage"),
                "guess": {
                    k: g[k]
                    for k in ("tier", "legacy_tier", "label", "evidence_status", "origin", "confidence")
                },
            }
        )
    length = budget.get("chromosome_length")
    out = {
        "chrom": chrom,
        "review_item": "R7 follow-up: origin beside role, constraint as selection only",
        "restated_from": f"budget_{chrom}",
        "chromosome_length": length,
        "threshold": budget.get("threshold"),
        "unknown_bp": budget["unknown_bp"],
        "certainty": {
            "model_score_name": SCORE_NAME,
            "probability_unavailable": NO_PROBABILITY,
            "effect": f"phylop.fraction_above of each row, in {CONSTRAINT_UNIT}",
        },
        "blocks": rows,
    }
    out.update(tallies(rows, length))
    return out


def restate_and_save(chrom: str, results_dir: Path | None = None) -> dict:
    out = restate(chrom, results_dir)
    kw = {"results_dir": results_dir} if results_dir else {}
    save_result(
        f"{AXES_NAME}_{chrom}", out, manifest=_manifest(chrom, ["budget", "unknown"], True, results_dir), **kw
    )
    return out
