# SPDX-License-Identifier: AGPL-3.0-or-later
"""Does CTCF motif orientation say which inferred node boundaries are real?

    uv run python scripts/ctcf_orientation.py [chr21 ...] [--relative 0.85]

The node boundaries are CTCF-only ENCODE elements with no strand. Loop extrusion is blocked
by CTCF in an orientation-dependent way, so loops and domains are anchored by convergent
sites. Each CTCF-only element is scanned with JASPAR MA0139 (the profile set genome/motifs.py
reads), each domain edge takes the strands of its elements (genome/domains.py
orient_boundaries), and the edges are held against 4DN's measured boundary calls within
20 kb, split by class, with two controls:

- the comparison's own random control: as many positions placed uniformly at random;
- a strand shuffle: motif strands permuted among the elements that have a hit, positions and
  motif presence kept, so a convergent subset that beats its shuffles is orientation and not
  the presence of a strong site.

Results: domains_ctcf_orientation (chr21 as asked, and pooled over every chromosome with
measured boundaries). No model is called.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos.coords import Locus  # noqa: E402
from genomeos.genome import IndexedGenome  # noqa: E402
from genomeos.genome import hic as hic_mod  # noqa: E402
from genomeos.genome import motifs as mo  # noqa: E402
from genomeos.genome.domains import CTCF_PROFILE, ORIENTATION_EVIDENCE, orient_boundaries  # noqa: E402
from genomeos.genome.regulatory import load_ccres  # noqa: E402
from genomeos.predict.contact_maps import TOLERANCE, compare, random_control  # noqa: E402
from genomeos.results import save_result  # noqa: E402

BIOSOURCES = ("GM12878", "H1-hESC", "K562", "HepG2", "IMR-90")
CLASSES = ("convergent", "divergent", "motif", "none")
SHUFFLES = 500
MIN_FLIP_GAP = 5_000  # consecutive sites closer than a merge distance are one boundary cluster
FLANK = 20_000  # the positive control looks this far either side of a measured boundary


def element_strands(chrom: str, ccres, motifs, relative: float) -> dict[str, set[str]]:
    g = IndexedGenome(f"data/reference/{chrom}.fa.gz")
    out: dict[str, set[str]] = {}
    try:
        for c in ccres:
            if c.cls != "CTCF-only":
                continue
            seq = str(g.fetch(Locus(chrom, c.start, min(c.end, g.lengths[chrom]))))
            hits = mo.all_hits(motifs, seq, min_relative=relative)
            out[c.id] = {"+" if h[3] == 0 else "-" for h in hits}
    finally:
        g.close()
    return out


def supported(positions: list[int], measured: list[int]) -> list[bool]:
    import bisect

    ms = sorted(measured)
    flags = []
    for p in positions:
        i = bisect.bisect_left(ms, p)
        d = min((abs(ms[j] - p) for j in (i - 1, i) if 0 <= j < len(ms)), default=None)
        flags.append(d is not None and d <= TOLERANCE)
    return flags


def main() -> None:
    args = sys.argv[1:]
    relative = mo.RELATIVE_SCORE
    if "--relative" in args:
        i = args.index("--relative")
        relative = float(args[i + 1])
        del args[i : i + 2]
    chroms = args or [f"chr{i}" for i in range(1, 23)] + ["chrX"]
    motifs = [
        m for m in mo.load_motifs(relative=min(relative, mo.RELATIVE_SCORE)) if m.id.startswith(CTCF_PROFILE)
    ]
    if not motifs:
        raise SystemExit("JASPAR MA0139 not found in the profile set")
    per_chrom: dict = {}
    rng = random.Random(139)
    shuffle_conv: dict[str, list[float]] = {b: [0.0] * SHUFFLES for b in BIOSOURCES}
    shuffle_div: dict[str, list[float]] = {b: [0.0] * SHUFFLES for b in BIOSOURCES}
    shuffle_n: dict[str, list[list[int]]] = {b: [[0, 0, 0, 0] for _ in range(SHUFFLES)] for b in BIOSOURCES}
    elements_total = elements_with_site = 0
    for chrom in chroms:
        ccres = load_ccres(chrom)
        ref = Path(f"data/reference/{chrom}.fa.gz")
        if not ccres or not ref.exists():
            continue
        g = IndexedGenome(ref)
        length = g.lengths[chrom]
        g.close()
        strands = element_strands(chrom, ccres, motifs, relative)
        elements_total += len(strands)
        elements_with_site += sum(1 for s in strands.values() if s)
        rows = orient_boundaries(length, ccres, strands)
        positions = [r["position"] for r in rows]
        measured = {b: hic_mod.load_boundaries(b, chrom) for b in BIOSOURCES}
        flags = {b: supported(positions, m) for b, m in measured.items() if m}
        chrom_out: dict = {
            "boundaries": len(rows),
            "by_class": {k: sum(1 for r in rows if r["class"] == k) for k in CLASSES},
            "by_sites": {k: sum(1 for r in rows if r["sites"] == k) for k in ("+", "-", "both", "none")},
            "biosources": {},
        }
        # strand shuffles among elements with a site, keeping which elements have one and how many strands
        with_site = [i for i, s in strands.items() if s]
        shuffled_rows = []
        for _ in range(SHUFFLES):
            labels = [strands[i] for i in with_site]
            rng.shuffle(labels)
            sh = dict(strands)
            sh.update(zip(with_site, labels, strict=True))
            shuffled_rows.append([r["class"] for r in orient_boundaries(length, ccres, sh)])
        for b, fl in flags.items():
            res = {"all": compare(positions, measured[b])}
            res["all"]["random_control"] = random_control(positions, measured[b], length)
            for k in CLASSES:
                sub = [p for p, r in zip(positions, rows, strict=True) if r["class"] == k]
                res[k] = {
                    "n": len(sub),
                    "supported": sum(f for f, r in zip(fl, rows, strict=True) if r["class"] == k),
                    "random_control": random_control(sub, measured[b], length) if sub else None,
                }
            chrom_out["biosources"][b] = res
            for s_i, classes in enumerate(shuffled_rows):
                for ci, k in enumerate(CLASSES):
                    shuffle_n[b][s_i][ci] += sum(1 for c in classes if c == k)
                shuffle_conv[b][s_i] += sum(f for f, c in zip(fl, classes, strict=True) if c == "convergent")
                shuffle_div[b][s_i] += sum(f for f, c in zip(fl, classes, strict=True) if c == "divergent")
        # positive control for the convention: around each measured boundary, single-strand sites
        # downstream should lean + (anchoring the next domain) and upstream sites - (closing the last)
        single = sorted(
            ((c.start + c.end) // 2, next(iter(strands[c.id])))
            for c in ccres
            if c.cls == "CTCF-only" and len(strands.get(c.id, ())) == 1
        )
        import bisect

        sp = [x for x, _ in single]
        for b, m in measured.items():
            flank = chrom_out.setdefault("measured_flanks", {}).setdefault(
                b, {"upstream": {"+": 0, "-": 0}, "downstream": {"+": 0, "-": 0}}
            )
            for mb in m:
                lo, hi = bisect.bisect_left(sp, mb - FLANK), bisect.bisect_right(sp, mb + FLANK)
                for x, st in single[lo:hi]:
                    flank["upstream" if x < mb else "downstream"][st] += 1
        # an orientation-aware caller, exploratory: a boundary between consecutive single-strand sites
        # where the strand flips from - to + (the pair points away from each other, as sites flanking a
        # measured boundary do); + to - flips (sites pointing at each other, inside a loop) are the control
        flips: dict[str, list[int]] = {"-+": [], "+-": []}
        for (x1, s1), (x2, s2) in zip(single, single[1:], strict=False):
            key = s1 + s2
            if key in flips and x2 - x1 >= MIN_FLIP_GAP:
                flips[key].append((x1 + x2) // 2)
        for b, m in measured.items():
            if not m:
                continue
            fo = chrom_out.setdefault("flips", {}).setdefault(b, {})
            for key, pos in flips.items():
                fo[key] = {
                    "n": len(pos),
                    "supported": sum(supported(pos, m)),
                    "random_control": random_control(pos, m, length) if pos else None,
                }
        per_chrom[chrom] = chrom_out
        print(
            chrom,
            chrom_out["by_class"],
            {
                b: (
                    v["convergent"]["supported"],
                    v["convergent"]["n"],
                    v["all"]["fraction_inferred_supported"],
                )
                for b, v in chrom_out["biosources"].items()
            },
            flush=True,
        )

    def summary(b: str, chrom_list: list[str]) -> dict:
        n = {k: 0 for k in CLASSES}
        s = {k: 0 for k in CLASSES}
        rc_w = 0.0
        tot = sup = 0
        for c in chrom_list:
            v = per_chrom.get(c, {}).get("biosources", {}).get(b)
            if not v:
                continue
            tot += v["all"]["inferred"]
            sup += v["all"]["inferred_on_a_predicted_boundary"]
            rc_w += v["all"]["random_control"] * v["all"]["inferred"]
            for k in CLASSES:
                n[k] += v[k]["n"]
                s[k] += v[k]["supported"]
        out = {
            "boundaries": tot,
            "supported_share": round(sup / tot, 4) if tot else None,
            "random_control": round(rc_w / tot, 4) if tot else None,
            "by_class": {
                k: {"n": n[k], "supported_share": round(s[k] / n[k], 4) if n[k] else None} for k in CLASSES
            },
        }
        return out

    result: dict = {"profile": CTCF_PROFILE, "relative_score": relative, "tolerance": TOLERANCE}
    result["elements_scanned"] = elements_total
    result["elements_with_site"] = elements_with_site
    result["chr21"] = {b: summary(b, ["chr21"]) for b in BIOSOURCES if "chr21" in per_chrom}
    all_chroms = list(per_chrom)
    gw = {}
    for b in BIOSOURCES:
        sm = summary(b, all_chroms)
        if not sm["boundaries"]:
            continue
        nconv = sm["by_class"]["convergent"]["n"]
        ndiv = sm["by_class"]["divergent"]["n"]
        conv_share = sm["by_class"]["convergent"]["supported_share"]
        div_share = sm["by_class"]["divergent"]["supported_share"]
        sh_conv = sorted(
            (shuffle_conv[b][i] / shuffle_n[b][i][0]) if shuffle_n[b][i][0] else 0.0 for i in range(SHUFFLES)
        )
        sh_div = sorted(
            (shuffle_div[b][i] / shuffle_n[b][i][1]) if shuffle_n[b][i][1] else 0.0 for i in range(SHUFFLES)
        )
        sm["strand_shuffle"] = {
            "shuffles": SHUFFLES,
            "convergent_supported_share_median": round(sh_conv[SHUFFLES // 2], 4),
            "convergent_95": [
                round(sh_conv[int(0.025 * SHUFFLES)], 4),
                round(sh_conv[int(0.975 * SHUFFLES) - 1], 4),
            ],
            "p_convergent_at_or_above": round(
                sum(1 for x in sh_conv if x >= (conv_share or 0)) / SHUFFLES, 4
            ),
            "convergent_n_median": sorted(x[0] for x in shuffle_n[b])[SHUFFLES // 2],
            "divergent_supported_share_median": round(sh_div[SHUFFLES // 2], 4),
            "p_divergent_at_or_below": round(sum(1 for x in sh_div if x <= (div_share or 0)) / SHUFFLES, 4),
            "divergent_n_median": sorted(x[1] for x in shuffle_n[b])[SHUFFLES // 2],
        }
        sm["convergent_n"], sm["divergent_n"] = nconv, ndiv
        fl = {"upstream": {"+": 0, "-": 0}, "downstream": {"+": 0, "-": 0}}
        for c in all_chroms:
            f = per_chrom[c].get("measured_flanks", {}).get(b)
            if f:
                for side in fl:
                    for st in "+-":
                        fl[side][st] += f[side][st]
        for side in fl:
            tot_side = fl[side]["+"] + fl[side]["-"]
            fl[side]["plus_share"] = round(fl[side]["+"] / tot_side, 4) if tot_side else None
        sm["positive_control_sites_within_20kb_of_measured_boundaries"] = fl
        flip = {}
        for key in ("-+", "+-"):
            n = sup = 0
            rc = 0.0
            for c in all_chroms:
                f = per_chrom[c].get("flips", {}).get(b, {}).get(key)
                if f and f["n"]:
                    n += f["n"]
                    sup += f["supported"]
                    rc += f["random_control"] * f["n"]
            flip[key] = {
                "n": n,
                "supported_share": round(sup / n, 4) if n else None,
                "random_control": round(rc / n, 4) if n else None,
            }
        sm["strand_flip_boundaries"] = flip
        gw[b] = sm
    result["genome"] = {"chromosomes": all_chroms, "biosources": gw}
    result["per_chromosome"] = per_chrom
    result["evidence"] = {
        "orientation": ORIENTATION_EVIDENCE,
        "measured": hic_mod.EVIDENCE,
        "controls": "uniform random positions (the comparison's own control); strand labels shuffled among "
        "elements with a site, positions and presence kept",
    }
    name = (
        "domains_ctcf_orientation"
        if relative == mo.RELATIVE_SCORE
        else f"domains_ctcf_orientation_{relative}"
    )
    save_result(name, result)
    for b, sm in gw.items():
        print(
            b,
            sm["supported_share"],
            sm["random_control"],
            {k: v for k, v in sm["by_class"].items()},
            sm["strand_shuffle"],
            flush=True,
        )


if __name__ == "__main__":
    main()
