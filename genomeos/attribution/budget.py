# SPDX-License-Identifier: AGPL-3.0-or-later
"""The composition budget: every UNKNOWN block with constraint attached and a best guess.

The house argument. A plan that calls half the house pipes is wrong before a brick is
laid, because the quantities do not fit. The genome has quantities too: a tenth of it
is under purifying selection across mammals and the rest has mutated freely for tens
of millions of years. So a block's sequence class (from `genomeos unknown`) and its
constraint (from Zoonomia and the 100-vertebrate elements) together support a guess
that is far better than "no clue":

    structural           gaps, centromeres, satellite arrays: sequence with a mechanical role
    fossil               transposable-element remains without constraint
    regulatory           ENCODE-backed element clusters; the target gene is the open question
    constrained_unknown  under selection, not coding, not regulatory by the registry: the real unknown
    neutral              unique sequence with no constraint and no element: best guess, nothing

Every guess carries a confidence and the numbers it rests on. The tiers are a
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
from genomeos.results import load_result, save_result

TIERS = ("structural", "fossil", "regulatory", "constrained_unknown", "neutral")
NEUTRAL_MAX = 0.03  # below this constrained fraction a block reads as unconstrained
CONSTRAINED_MIN = 0.05  # from here on the block carries constraint worth a name
STRUCTURAL = {"centromere", "satellite_array", "tandem_repeat"}
REGULATORY = {"regulatory", "promoter_like"}


def guess(cls: str, confidence: float, phylop: dict | None, elements: dict | None) -> dict:
    """Class plus constraint to a tier, a label and a confidence."""
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
                "label": "repeat-derived sequence under constraint: possibly exapted as a regulatory element",
                "confidence": 0.5,
            }
        if fa < NEUTRAL_MAX:
            return {
                "tier": "fossil",
                "label": "transposable-element fossil, unconstrained",
                "confidence": 0.8,
            }
        return {"tier": "fossil", "label": "transposable-element fossil, weak constraint", "confidence": 0.6}
    if cls in REGULATORY:
        if fa >= CONSTRAINED_MIN:
            return {
                "tier": "regulatory",
                "label": "regulatory elements under constraint; target gene unassigned",
                "confidence": 0.7,
            }
        return {
            "tier": "regulatory",
            "label": "regulatory elements by the registry, little constraint: lineage-specific or weak",
            "confidence": 0.5,
        }
    if cls == "long_orf":
        if fa >= CONSTRAINED_MIN:
            return {
                "tier": "constrained_unknown",
                "label": "long open reading frame under constraint: unannotated coding or a young pseudogene",
                "confidence": 0.5,
            }
        return {
            "tier": "neutral",
            "label": "long open reading frame without constraint: a chance frame or a dead pseudogene",
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
            "tier": "neutral",
            "label": "unconstrained unique sequence: no evidence of function, best guess neutral",
            "confidence": 0.6,
        }
    return {
        "tier": "neutral",
        "label": "weakly constrained unique sequence: mostly neutral",
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
        g = guess(b["class"], b.get("confidence", 0.0), p, e)
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


def tallies(rows: list[dict], length: int | None) -> dict:
    by_class: dict[str, dict] = {}
    by_tier: dict[str, dict] = {t: {"blocks": 0, "bp": 0, "constrained_bp": 0} for t in TIERS}
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
        t = by_tier[r["guess"]["tier"]]
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


def distil(results_dir: Path | None = None) -> dict:
    """Every budget_chr*.json summed: the genome's composition budget."""
    kw = {"results_dir": results_dir} if results_dir else {}
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    per: dict[str, dict] = {}
    by_class: dict[str, dict] = {}
    by_tier: dict[str, dict] = {t: {"blocks": 0, "bp": 0, "constrained_bp": 0} for t in TIERS}
    genome = 0
    unknown_bp = 0
    measured = 0
    constrained = 0
    guessed = 0
    cost_mb = 0.0
    for c in chroms:
        r = load_result(f"budget_{c}", **kw)
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
            for f in ("blocks", "bp", "constrained_bp"):
                by_tier[t][f] += v.get(f) or 0
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


def run_and_save(chrom: str, **kw) -> dict:
    out = build(chrom, **kw)
    save_result(f"budget_{chrom}", out)
    return out
