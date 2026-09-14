# SPDX-License-Identifier: AGPL-3.0-or-later
"""Terminal fates from measured factors (area E, Next step 1).

Until now every terminal fate of the worm program came from the observed lineage lookup
(`fate_<cell>` decisions generated from Sulston's tree), so the fate score was 100% by construction. Here
the fate is decided by the transcription factors the cell carries (Ma et al. 2021 atlas) and the lookup is
only the fallback where no factor rule applies.

What a terminal cell reads:

    instantaneous   the presence set of the cell, or of its nearest ancestor the atlas tracked (what the
                    reader gives it in the Body today): max adjusted expression >= 20% of the factor's max
    integrated      exposure summed along its lineage path from the zygote (minutes at the factor's maximum
                    level, `atlas_levels`), called present at >= 20% of the factor's largest path exposure
                    over the embryonic terminal cells: the same one-parameter threshold, read over time

Two rule sets:

    textbook   the factor -> tissue relations already cited in `tf_atlas.TEXTBOOK` (ELT-2 and ELT-7 for
               intestine, HLH-1 and UNC-120 for body-wall muscle, ELT-1, LIN-26 and NHR-25 for hypodermis);
               no fitting, so their score is out of sample. PHA-4 and CEH-22 name an organ (pharynx), not a
               WormWeb tissue, and END-1 is a specification factor that precedes the gut programme, so they
               are not terminal rules.
    learned    a decision list learned from the atlas (a conjunction of at most two factors per rule), scored
               by leaving one founder sublineage out at a time (ABal, ABar, ABpl, ABpr, MS, E, C, D), so a
               rule is only credited on a lineage it did not see.

Rules are applied in order; the first that applies decides. In the Body, the rules are compiled into
`differentiate` decisions gated on a terminal cell type (the lookup still says which cells stop dividing and
when) and on the embryonic stages (the atlas ends at the bean stage and says nothing about larval cells).
"""

from __future__ import annotations

import random
from collections import Counter, defaultdict
from collections.abc import Iterable
from dataclasses import dataclass

from .reference import HATCH_MIN, TISSUE_CELL_TYPE, ReferenceLineage
from .tf_atlas import PRESENCE_FRACTION, SOURCE

EMBRYONIC_STAGES = "Cleavage|Gastrulation|Morphogenesis"
TERMINAL_TYPES = "|".join(sorted(set(TISSUE_CELL_TYPE.values())))
SUBLINEAGES = ("ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C", "D")
INTEGRATED_SUFFIX = "_integrated"

# (tissue, factors any of which calls it, citation); order is precedence
TEXTBOOK_RULES = [
    (
        "intestine",
        ("ELT-2", "ELT-7"),
        "Fukushige et al. 1998, Dev Biol 198:286; Sommermann et al. 2010, Dev Biol 347:154",
    ),
    (
        "muscle",
        ("HLH-1", "UNC-120"),
        "Krause et al. 1990, Cell 63:907; Fukushige et al. 2006, Genes Dev 20:3395",
    ),
    (
        "hypoderm",
        ("ELT-1", "LIN-26", "NHR-25"),
        "Page et al. 1997, Genes Dev 11:1651; Labouesse et al. 1994, Development 120:2359; "
        "Gissendanner & Sluder 2000, Dev Biol 221:259",
    ),
]


@dataclass(frozen=True, slots=True)
class Rule:
    tissue: str
    factors: tuple[str, ...]  # all present (a conjunction); a textbook rule is one Rule per factor
    source: str = ""
    precision: float = 0.0  # on the data it was learned from
    support: int = 0

    def applies(self, present: set[str]) -> bool:
        return all(f in present for f in self.factors)


def textbook_rules() -> list[Rule]:
    return [Rule(t, (f,), src) for t, fs, src in TEXTBOOK_RULES for f in fs]


# ---- what a terminal cell reads --------------------------------------------


def embryonic_terminal(ref: ReferenceLineage) -> list:
    return [c for c in ref.terminal() if c.born < HATCH_MIN]


def _path(ref: ReferenceLineage, cid: str) -> list[str]:
    out, x = [], ref.cells.get(cid)
    while x is not None:
        out.append(x.id)
        x = ref.cells.get(x.parent) if x.parent else None
    return out


def instantaneous(ref: ReferenceLineage, atlas: dict[str, list[str]], cells: Iterable) -> dict[str, set[str]]:
    """cell -> factors of the cell or of its nearest tracked ancestor (the reader's inheritance)."""
    out = {}
    for c in cells:
        tracked = next((x for x in _path(ref, c.id) if x in atlas), None)
        out[c.id] = set(atlas[tracked]) if tracked else set()
    return out


def presence(levels: dict[str, dict[str, list]], fraction: float = PRESENCE_FRACTION, field: int = 0) -> dict:
    """cell -> factors whose `field` (0 peak, 1 mean over the cell's frames) reaches `fraction` of that
    field's maximum over all cells."""
    out: dict[str, list[str]] = defaultdict(list)
    for tf, per_cell in levels.items():
        top = max((v[field] for v in per_cell.values()), default=0.0)
        for cell, v in per_cell.items():
            if top > 0 and v[field] >= fraction * top:
                out[cell].append(tf)
    return dict(out)


def exposures(ref: ReferenceLineage, levels: dict[str, dict[str, list]], cells: Iterable) -> dict[str, dict]:
    """cell -> factor -> exposure summed over the cell and its ancestors (minutes at the factor's max)."""
    out: dict[str, dict[str, float]] = {}
    for c in cells:
        path = _path(ref, c.id)
        acc: dict[str, float] = {}
        for tf, per_cell in levels.items():
            s = sum(per_cell[x][2] for x in path if x in per_cell)
            if s > 0:
                acc[tf] = s
        out[c.id] = acc
    return out


def integrated(exp: dict[str, dict[str, float]], fraction: float = PRESENCE_FRACTION) -> dict[str, set[str]]:
    top: dict[str, float] = defaultdict(float)
    for acc in exp.values():
        for tf, v in acc.items():
            top[tf] = max(top[tf], v)
    return {cid: {tf for tf, v in acc.items() if v >= fraction * top[tf]} for cid, acc in exp.items()}


# ---- applying and scoring ----------------------------------------------------


def apply_rules(rules: list[Rule], feats: dict[str, set[str]]) -> dict[str, str]:
    pred = {}
    for cid, present in feats.items():
        for r in rules:
            if r.applies(present):
                pred[cid] = r.tissue
                break
    return pred


def score(pred: dict[str, str], labels: dict[str, str]) -> dict:
    """With the lookup as fallback: a cell a rule decides is right or wrong; the rest keep the observed fate.
    Per tissue: cells, decided by a factor rule and right (recall by factors), wrongly taken from it by
    another tissue's rule (lost), and cells the tissue's rules claimed (precision)."""
    per: dict[str, dict[str, int]] = defaultdict(
        lambda: {"cells": 0, "factor_correct": 0, "lost": 0, "claimed": 0}
    )
    for cid, t in labels.items():
        per[t]["cells"] += 1
        p = pred.get(cid)
        if p is None:
            continue
        per[p]["claimed"] += 1
        if p == t:
            per[t]["factor_correct"] += 1
        else:
            per[t]["lost"] += 1
    wrong = sum(1 for cid, p in pred.items() if labels.get(cid) != p)
    n = len(labels)
    tissues = {}
    for t, r in sorted(per.items(), key=lambda kv: -kv[1]["cells"]):
        tissues[t] = {
            **r,
            "after": round((r["cells"] - r["lost"]) / r["cells"], 3) if r["cells"] else None,
            "factor_recall": round(r["factor_correct"] / r["cells"], 3) if r["cells"] else None,
            "precision": round(r["factor_correct"] / r["claimed"], 3) if r["claimed"] else None,
        }
    confusion = Counter(f"{labels[c]}->{p}" for c, p in pred.items() if labels.get(c) != p)
    return {
        "cells": n,
        "decided_by_factors": len(pred),
        "factor_correct": len(pred) - wrong,
        "wrong": wrong,
        "before": 1.0,
        "after": round((n - wrong) / n, 4) if n else None,
        "factor_only": round((len(pred) - wrong) / n, 4) if n else None,
        "precision": round((len(pred) - wrong) / len(pred), 4) if pred else None,
        "tissues": tissues,
        "confusion": dict(confusion.most_common(10)),
    }


# ---- a decision list learned from the atlas -----------------------------------


def learn(
    feats: dict[str, set[str]],
    labels: dict[str, str],
    min_support: int = 3,
    min_precision: float = 0.8,
    max_rules: int = 40,
) -> list[Rule]:
    """Greedy sequential covering: the conjunction of one or two factors with the best Laplace precision for
    any tissue among the cells not yet covered, until no rule reaches `min_precision` on `min_support`."""
    cells = sorted(feats)
    idx = {c: i for i, c in enumerate(cells)}
    bits: dict[str, int] = defaultdict(int)
    for c in cells:
        for f in feats[c]:
            bits[f] |= 1 << idx[c]
    tissue_bits: dict[str, int] = defaultdict(int)
    for c in cells:
        tissue_bits[labels[c]] |= 1 << idx[c]
    factors = sorted(bits)
    remaining = (1 << len(cells)) - 1
    rules: list[Rule] = []
    while len(rules) < max_rules:
        best: tuple[float, int, str, tuple[str, ...], float] | None = None
        singles = [(f, bits[f] & remaining) for f in factors]
        singles = [(f, b) for f, b in singles if b.bit_count() >= min_support]
        candidates: list[tuple[tuple[str, ...], int]] = [((f,), b) for f, b in singles]
        for i, (f, bf) in enumerate(singles):
            for g, bg in singles[i + 1 :]:
                b = bf & bg
                if b.bit_count() >= min_support and b != bf and b != bg:
                    candidates.append(((f, g), b))
        for combo, b in candidates:
            n = b.bit_count()
            for t, tb in tissue_bits.items():
                k = (b & tb).bit_count()
                if k < min_support:
                    continue
                lap = (k + 1) / (n + 2)
                key = (lap, k, t, combo, k / n)
                if best is None or key[:2] > best[:2]:
                    best = key
        if best is None or best[4] < min_precision:
            break
        lap, k, t, combo, prec = best
        covered = remaining
        for f in combo:
            covered &= bits[f]
        rules.append(Rule(t, combo, "learned from Ma 2021 x WormWeb", round(prec, 3), k))
        remaining &= ~covered
    return rules


def sublineage(cid: str) -> str:
    return next((s for s in SUBLINEAGES if cid.startswith(s)), "other")


def lineage_majority(labels: dict[str, str]) -> dict[str, str]:
    """No factors: every cell gets the most common tissue of its founder sublineage (in sample)."""
    by = defaultdict(Counter)
    for c, t in labels.items():
        by[sublineage(c)][t] += 1
    return {c: by[sublineage(c)].most_common(1)[0][0] for c in labels}


def cross_validate(
    feats: dict[str, set[str]],
    labels: dict[str, str],
    before: list[Rule] | None = None,
    tissues: set[str] | None = None,
    **kw,
) -> dict[str, str]:
    """Leave one founder sublineage out: rules learned on the other seven decide the held-out cells. `before`
    are fixed rules with precedence (the textbook), `tissues` keeps only learned rules for those tissues."""
    pred: dict[str, str] = {}
    groups = defaultdict(list)
    for c in feats:
        groups[sublineage(c)].append(c)
    for g, held in groups.items():
        train = {c: feats[c] for c in feats if sublineage(c) != g}
        rules = learn(train, {c: labels[c] for c in train}, **kw)
        if tissues is not None:
            rules = [r for r in rules if r.tissue in tissues]
        pred.update(apply_rules((before or []) + rules, {c: feats[c] for c in held}))
    return pred


def label_null(feats: dict[str, set[str]], labels: dict[str, str], n: int = 20, seed: int = 0, **kw) -> dict:
    """Cross-validated accuracy when tissues are shuffled among cells of the same founder sublineage: keeps
    what the lineage alone says about composition, removes what the factors say about each cell."""
    rng = random.Random(seed)
    by_group = defaultdict(list)
    for c in labels:
        by_group[sublineage(c)].append(c)
    scores = []
    for _ in range(n):
        shuffled = {}
        for cs in by_group.values():
            ts = [labels[c] for c in cs]
            rng.shuffle(ts)
            shuffled.update(zip(cs, ts, strict=True))
        pred = cross_validate(feats, shuffled, **kw)
        scores.append(sum(1 for c, p in pred.items() if shuffled[c] == p) / len(labels))
    return {"permutations": n, "factor_only": sorted(round(s, 4) for s in scores)}


# ---- BioLang -------------------------------------------------------------------


def to_bio_fates(rules: list[Rule], integrated_read: bool, module: str = "organism.celegans.fates") -> str:
    """Rules as `differentiate` decisions with an explicit `priority`: the first rule of the list is the
    highest, and every factor rule outranks the lineage lookup (priority 0). The program that imports them
    declares `regime { fates: first }` (BIOLANG-v0.4-ECONOMY.md §7.3), so one fate is taken per decision
    point and line order carries no meaning."""
    read = "integrated along the lineage path" if integrated_read else "instantaneous presence"
    lines = [
        "# Generated by scripts/celegans_fates.py (genomeos.organism.fate_rules): terminal fates from the",
        f"# measured factors (Ma et al. 2021), read as {read}. Imported after the lineage program, so a",
        "# factor rule overrides the observed fate of a terminal cell and the lookup stays the fallback.",
        "# Precedence is the `priority` of each decision (first rule highest; the lookup has 0) under",
        "# `regime { fates: first }`, declared by the importing program.",
        f"module {module}",
        "import cell_types.bio",
        "",
    ]
    suffix = INTEGRATED_SUFFIX if integrated_read else ""
    for i, r in enumerate(rules):
        cond = ", ".join(f"{f}{suffix} = present" for f in r.factors)
        kind = "experimental" if r.source and "learned" not in r.source else "inferred"
        src = r.source or "learned from Ma 2021 x WormWeb"
        if kind == "inferred":
            src = f"{src} (precision {r.precision}, n {r.support})"
        conf = 0.8 if kind == "experimental" else round(min(0.7, r.precision * 0.75), 2)
        lines.append(
            f"decision factor_fate_{i:02d} {{ action: differentiate; when: cell_type = {TERMINAL_TYPES}, "
            f"stage = {EMBRYONIC_STAGES}, {cond}; to: {TISSUE_CELL_TYPE[r.tissue]}; "
            f"priority: {len(rules) - i}; "
            f'evidence: {kind} "{src}; factors {SOURCE}"; confidence: {conf} }}'
        )
    return "\n".join(lines) + "\n"


def to_bio_exposure(
    feats_integrated: dict[str, set[str]],
    feats_instant: dict[str, set[str]],
    atlas: dict[str, list[str]],
    factors: Iterable[str],
    module: str = "organism.celegans.exposure",
) -> str:
    """One `express` decision per embryonic terminal cell with the integrated calls the rules read. A cell
    the atlas did not track also gets its inherited instantaneous set, because an express decision replaces
    the measured factors a cell inherited."""
    wanted = set(factors)
    lines = [
        "# Generated by scripts/celegans_fates.py: exposure to each rule factor summed along the",
        f"# lineage path (Ma et al. 2021 frames x 1.25 min), present at >= {PRESENCE_FRACTION:.0%} of the",
        "# factor's largest path exposure over the embryonic terminal cells.",
        f"module {module}",
        "",
    ]
    for cid in sorted(feats_integrated):
        sets = sorted(f"{f}{INTEGRATED_SUFFIX}" for f in feats_integrated[cid] & wanted)
        if cid not in atlas:
            sets = sorted(feats_instant[cid]) + sets
        if not sets:
            continue
        lines.append(
            f"decision exposure_{cid} {{ action: express; when: cell = {cid}; sets: {', '.join(sets)}; "
            f'evidence: experimental "{SOURCE}"; confidence: 0.7 }}'
        )
    return "\n".join(lines) + "\n"
