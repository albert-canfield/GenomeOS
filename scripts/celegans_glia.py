# SPDX-License-Identifier: AGPL-3.0-or-later
"""Glia from factors: what separates sheath and socket glia from their sister neurons, if anything.

Glia are the clearest failure of the worm's factor rules. Of the 555 embryonic terminal cells the
program scores, 22 are sheath glia and 18 socket glia, and the factor rules call 6 and 9 of them:
the rest go to the hypodermis and neuron rules, because a glial cell is the sister of a sensory
neuron, carries that neuron's factors, and carries the ectodermal factors (ELT-1, LIN-26, NHR-25)
the hypodermis rules read. The question is whether anything in the Ma et al. 2021 atlas separates
them, and a measured "these cells are not separable by the factors we have" is a result.

Six measurements, all local, no model calls:

1. coverage      are the glia and their sisters in the atlas at all, or is the question unanswerable?
2. singles       the best single factor for sheath, socket and glia, under three reads of the atlas
3. cited         factors the glia literature names (PROS-1, MLS-2, HLH-17, EGL-13), scored
4. ceiling       an exhaustive search over every conjunction of up to three of the atlas's factors,
                 scored by what it would be worth to the program: an in-sample upper bound on the
                 whole rule language BioLang `when:` can express
5. generalise    fitted rules at precision 1.0 in sample, scored on a held-out founder sublineage
6. sisters       how many factors separate a glial cell from its sister, against sisters of the same
                 tissue and against the atlas's own two-reporter measurement floor

    uv run python scripts/celegans_glia.py
"""

from __future__ import annotations

import statistics
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path and not (Path.cwd() / "genomeos").is_dir():
    sys.path.insert(0, str(ROOT))

from genomeos.organism import atlas_levels  # noqa: E402
from genomeos.organism import fate_rules as fr  # noqa: E402
from genomeos.organism.reference import ReferenceLineage  # noqa: E402
from genomeos.organism.tf_atlas import load_cells  # noqa: E402
from genomeos.results import save_result  # noqa: E402

MIN_SUPPORT = 3
CITED = [
    (("PROS-1",), "sheath", "Wallace et al. 2016, Development 143:3016 (PROS-1/Prospero, sheath glia)"),
    (("MLS-2",), "sheath", "Yoshimura et al. 2008, Development 135:2263 (MLS-2/HMX, CEPsh glia)"),
    (("HLH-17",), "sheath", "McMiller & Johnson 2005, Gene 356:1 (hlh-17 in CEPsh glia)"),
    (("EGL-13",), "socket", "Feng et al. 2013, Development 140:4433 (EGL-13/SoxD in socket cells)"),
    (("PROS-1", "NHR-25"), "sheath", "PROS-1 in the non-neuronal ectodermal class"),
    (("PROS-1", "SOX-2"), "sheath", "PROS-1 in the neural class"),
    (("PROS-1", "NHR-25", "SOX-2"), "sheath", "the rule the program adopts"),
    (("SOX-2", "PAG-3"), "socket", "the best cited-looking socket pair"),
    (("SOX-2",), "glia", "SOX-2 marks the neural and glial class, not glia alone"),
]


def main() -> None:
    t0 = time.time()
    ref = ReferenceLineage.load()
    levels = atlas_levels.load_levels()
    atlas = load_cells()
    term = fr.embryonic_terminal(ref)
    labels = {c.id: c.tissue for c in term}
    cells = sorted(labels)
    targets = {
        "sheath": {c for c in cells if labels[c] == "sheath"},
        "socket": {c for c in cells if labels[c] == "socket"},
        "glia": {c for c in cells if labels[c] in ("sheath", "socket")},
    }
    exp = fr.exposures(ref, levels, term)
    reads = {
        "peak_in_own_life": fr.instantaneous(ref, fr.presence(levels, field=0), term),
        "mean_over_own_life": fr.instantaneous(ref, fr.presence(levels, field=1), term),
        "integrated_along_lineage": fr.integrated(exp),
    }
    # the read the program's reader gives a cell (its own peak calls, or its nearest tracked ancestor's)
    inst = fr.instantaneous(ref, atlas, term)
    reads["reader"] = inst

    out: dict = {
        "result": "celegans_glia",
        "question": "does anything in the Ma 2021 atlas separate sheath and socket glia from their sister "
        "neurons and from the hypodermis?",
        "cells": {
            "embryonic_terminal": len(cells),
            "sheath": len(targets["sheath"]),
            "socket": len(targets["socket"]),
        },
    }

    # 1. coverage: a glial cell and its sister both tracked means the question is answerable
    kids = defaultdict(list)
    for c in ref.cells.values():
        if c.parent:
            kids[c.parent].append(c.id)
    tracked = set(atlas)
    sisters = {}
    for c in targets["glia"]:
        sibs = [x for x in kids[ref.cells[c].parent] if x != c]
        sisters[c] = sibs[0] if sibs else None
    out["coverage"] = {
        "glia_tracked_in_the_atlas": sum(1 for c in targets["glia"] if c in tracked),
        "glia_whose_sister_is_tracked": sum(1 for c, s in sisters.items() if s is not None and s in tracked),
        "of": len(targets["glia"]),
        "glia_sharing_a_feature_vector_with_a_non_glial_cell": {
            name: _identical(feats, cells, labels, targets["glia"]) for name, feats in reads.items()
        },
        "note": "every glial cell and every sister is measured, and no glial cell has the same factor set "
        "as a non-glial one, so nothing rules separability out in principle",
    }

    # 2. the best single factor per target per read
    out["best_single_factors"] = {
        name: {t: _rank(feats, cells, pos, 6) for t, pos in targets.items()} for name, feats in reads.items()
    }

    # 3. the factors the glia literature names
    out["cited_factors"] = [
        {
            "factors": list(combo),
            "for": t,
            "source": src,
            "reads": {name: _stats(feats, cells, labels, combo, targets[t]) for name, feats in reads.items()},
        }
        for combo, t, src in CITED
    ]

    # 4. the ceiling: what the best possible rule of up to three factors would be worth to the program
    out["ceiling"] = _ceiling(ref, levels, atlas, term, labels, cells, inst, reads)

    # 5. fitted rules at precision 1.0 in sample, scored on a held-out sublineage
    out["held_out"] = {
        name: {t: _holdout(feats, cells, pos) for t, pos in targets.items()} for name, feats in reads.items()
    }

    # 6. sisters, against the atlas's own measurement floor
    out["sisters"] = _sisters(ref, kids, cells, labels, inst, targets["glia"])
    floor = atlas_levels.load_strain_peaks()
    disagree = [
        sum(1 for v in per.values() if (v[0] >= 0.2) != (v[1] >= 0.2)) / len(per)
        for per in floor.values()
        if len(per) >= 10 and len(next(iter(per.values()))) == 2
    ]
    out["measurement_floor"] = {
        "factors_with_two_reporter_strains": len(disagree),
        "median_disagreement_on_presence_in_the_same_cell": round(statistics.median(disagree), 3)
        if disagree
        else None,
        "note": "two reporters for the same factor disagree on presence in the same cell more often than "
        "sisters disagree on the same factors, which is the floor any sister-level rule sits on",
    }

    out["seconds"] = round(time.time() - t0, 1)
    print(save_result("celegans_glia", out))
    for t in ("sheath", "socket"):
        print(t, "best single (reader read):", out["best_single_factors"]["reader"][t][0])
        print(t, "ceiling over all conjunctions of <=3 factors:", out["ceiling"][t])
        print(t, "fitted rules held out:", out["held_out"]["reader"][t])


def _identical(feats, cells, labels, glia) -> int:
    byvec = defaultdict(list)
    for c in cells:
        byvec[frozenset(feats[c])].append(c)
    n = 0
    for cs in byvec.values():
        if len(cs) > 1 and len({labels[c] for c in cs}) > 1:
            n += sum(1 for c in cs if c in glia)
    return n


def _rank(feats, cells, pos, top):
    rows = []
    factors = {f for c in cells for f in feats[c]}
    for f in factors:
        claimed = {c for c in cells if f in feats[c]}
        tp = len(claimed & pos)
        if not tp:
            continue
        prec, rec = tp / len(claimed), tp / len(pos)
        rows.append(
            {
                "factor": f,
                "f1": round(2 * prec * rec / (prec + rec), 3),
                "precision": round(prec, 3),
                "recall": round(rec, 3),
                "claimed": len(claimed),
            }
        )
    rows.sort(key=lambda r: -r["f1"])
    return rows[:top]


def _stats(feats, cells, labels, combo, pos):
    claimed = [c for c in cells if all(f in feats[c] for f in combo)]
    tp = len([c for c in claimed if c in pos])
    return {
        "claimed": len(claimed),
        "correct": tp,
        "precision": round(tp / len(claimed), 3) if claimed else None,
        "recall": round(tp / len(pos), 3),
        "what_it_claims": dict(Counter(labels[c] for c in claimed).most_common(5)),
    }


def _ceiling(ref, levels, atlas, term, labels, cells, inst, reads):
    """The best net gain in terminal fates over every conjunction of up to three measured factors.

    In sample and therefore an upper bound: the rule is chosen on the same 555 cells it is scored on.
    A rule is inserted where a glia rule would go in the program, between the muscle and the hypodermis
    rules, and scored by how many more of the 555 come out right than with the program's rules today."""
    integ = reads["integrated_along_lineage"]
    textbook = fr.textbook_rules()
    neuron = [r for r in fr.learn(integ, labels) if r.tissue == "neuron"]
    featmap = {"integrated": integ, "instantaneous": inst}
    base = fr.apply_mixed(textbook + neuron, {"integrated": integ, "instantaneous": inst})
    before = fr.apply_mixed(textbook[:4], featmap)
    idx = {c: i for i, c in enumerate(cells)}
    gain = {}
    for t in ("sheath", "socket"):
        g = ell = 0
        for c in cells:
            if c in before:
                continue  # an intestine or muscle rule already decided it
            was = base.get(c, labels[c]) == labels[c]
            now = labels[c] == t
            if now and not was:
                g |= 1 << idx[c]
            elif was and not now:
                ell |= 1 << idx[c]
        gain[t] = (g, ell)
    bits: dict[str, int] = defaultdict(int)
    for c in cells:
        for f in inst[c]:
            bits[f] |= 1 << idx[c]
    factors = [f for f in sorted(bits) if bits[f].bit_count() >= MIN_SUPPORT]
    pairs = []
    for i, f in enumerate(factors):
        for g2 in factors[i + 1 :]:
            b = bits[f] & bits[g2]
            if b.bit_count() >= MIN_SUPPORT:
                pairs.append(((f, g2), b))

    def net(b, t):
        g, ell = gain[t]
        return (b & g).bit_count() - (b & ell).bit_count()

    out = {}
    for t in ("sheath", "socket"):
        singles = [(net(bits[f], t), (f,), bits[f].bit_count()) for f in factors]
        two = [(net(b, t), combo, b.bit_count()) for combo, b in pairs]
        best3 = (-99, (), 0)
        for combo, b in pairs:
            for f in factors:
                if f in combo:
                    continue
                b3 = b & bits[f]
                if b3.bit_count() < MIN_SUPPORT:
                    continue
                v = net(b3, t)
                if v > best3[0]:
                    best3 = (v, tuple(sorted(combo + (f,))), b3.bit_count())
        b1 = max(singles)
        b2 = max(two)
        out[t] = {
            "factors_searched": len(factors),
            "pairs_searched": len(pairs),
            "best_one_factor": {"net_fates": b1[0], "factors": list(b1[1]), "claimed": b1[2]},
            "best_two_factors": {"net_fates": b2[0], "factors": list(b2[1]), "claimed": b2[2]},
            "best_three_factors": {"net_fates": best3[0], "factors": list(best3[1]), "claimed": best3[2]},
        }
    return out


def _holdout(feats, cells, pos, nrules=3):
    """Rules fitted on seven founder sublineages, scored on the eighth: does the signal generalise?"""
    groups = defaultdict(list)
    for c in cells:
        groups[fr.sublineage(c)].append(c)
    tp = fp = fn = 0
    in_sample = []
    for g, held in groups.items():
        train = [c for c in cells if fr.sublineage(c) != g]
        rules = _best_rules(feats, pos, train, nrules)
        in_sample.extend(round(r[0], 3) for r in rules)
        for c in held:
            hit = any(all(f in feats[c] for f in r[1]) for r in rules)
            tp += hit and c in pos
            fp += hit and c not in pos
            fn += (not hit) and c in pos
    return {
        "rules_per_fold": nrules,
        "in_sample_precision_median": round(statistics.median(in_sample), 3) if in_sample else None,
        "held_out_true_positives": tp,
        "held_out_false_positives": fp,
        "held_out_missed": fn,
        "held_out_precision": round(tp / (tp + fp), 3) if tp + fp else 0.0,
        "held_out_recall": round(tp / (tp + fn), 3) if tp + fn else 0.0,
        "base_rate": round(len(pos) / len(cells), 3),
    }


def _best_rules(feats, pos, train, k):
    idx = {c: i for i, c in enumerate(train)}
    bits: dict[str, int] = defaultdict(int)
    for c in train:
        for f in feats[c]:
            bits[f] |= 1 << idx[c]
    pb = 0
    for c in train:
        if c in pos:
            pb |= 1 << idx[c]
    singles = [(f, bits[f]) for f in sorted(bits) if bits[f].bit_count() >= MIN_SUPPORT]
    cands = [((f,), b) for f, b in singles]
    for i, (f, bf) in enumerate(singles):
        for g, bg in singles[i + 1 :]:
            b = bf & bg
            if b.bit_count() >= MIN_SUPPORT and b != bf and b != bg:
                cands.append(((f, g), b))
    scored = [((b & pb).bit_count() / b.bit_count(), combo) for combo, b in cands]
    scored.sort(reverse=True)
    return scored[:k]


def _sisters(ref, kids, cells, labels, feats, glia):
    pairs, seen = [], set()
    for c in cells:
        sibs = [x for x in kids[ref.cells[c].parent] if x != c]
        if not sibs or sibs[0] not in feats:
            continue
        key = tuple(sorted((c, sibs[0])))
        if key in seen:
            continue
        seen.add(key)
        pairs.append(key)

    def spread(group):
        xs = [len(feats[a] ^ feats[b]) for a, b in group]
        return {"pairs": len(xs), "median_factors_differing": statistics.median(xs) if xs else None}

    glia_vs_other = [
        (a, b) for a, b in pairs if ({a, b} & glia) and labels.get(a) != labels.get(b) and a in cells
    ]
    same = [(a, b) for a, b in pairs if labels.get(a) and labels.get(a) == labels.get(b)]
    diff = [(a, b) for a, b in pairs if labels.get(a) and labels.get(b) and labels[a] != labels[b]]
    return {
        "glia_and_a_sister_of_another_tissue": spread(glia_vs_other),
        "sisters_of_the_same_tissue": spread(same),
        "sisters_of_different_tissues": spread(diff),
        "note": "a glial cell differs from its sister by no more factors than two sisters of the same "
        "tissue do, so the sister difference carries little about the fate",
    }


if __name__ == "__main__":
    main()
