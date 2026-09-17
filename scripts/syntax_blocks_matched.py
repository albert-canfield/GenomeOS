# SPDX-License-Identifier: AGPL-3.0-or-later
"""Is the syntax case's high rate the blocks, or is it how close they sit to a promoter? No requests.

    uv run python scripts/syntax_blocks_matched.py [--chroms chr21,chr22]

The 69 syntax blocks read 0.575 of 40 already-scored elements against a tier average of 0.248, and
that comparison holds nothing fixed. Assembling the pre-registered tiling run (2026-09-17) turned up
the reason to doubt it before a single request was spent: on chr22 the syntax windows sit a median
**2.6 kb** from a coding transcription start, while the relaxed and neutral arms sit at 51 kb and
62 kb. An element near a promoter moves a gene more often whatever else is true of it.

This asks the cheap version of the question, over elements the sweep has already scored: do the
elements inside the 69 syntax blocks move a gene more often than elements matched to them on length,
GC and distance to the nearest coding TSS? Matching is by strata, the convention this project already
uses (`attribution/unknown_scoring`), and the same comparison is run for the relaxed case so the human
axis has its own answer.

If the difference survives the match, the tiling run is worth its 2,400 requests. If it does not, the
0.575 was proximity and the run would have measured proximity more precisely.
"""

from __future__ import annotations

import argparse
import bisect
import json
import time
from pathlib import Path
from typing import Any

from genomeos.attribution import organise
from genomeos.attribution import syntax_tiling as st
from genomeos.compare import Strata, standardised
from genomeos.results import save_result

ELEMENTS = Path("data/knowledge/alphagenome/all_elements")
GC_BINS = (0.35, 0.45, 0.55)
LEN_BINS = (200, 400)
TSS_BINS = (1_000, 5_000, 20_000, 100_000)


def label_of(block: dict[str, Any]) -> str:
    """Every UNKNOWN block labelled the way a comparison needs: tier, and for the real unknown its case.

    The 882-block reading compared "constrained unknown, copies removed" with the neutral tier as one
    number each. Matching showed the ordering between them is not what the raw rates said, so the label
    here is finer: each case of the real unknown separately, the copies apart, and the other tiers by
    name. `arm_blocks` in the tiling module keeps the three arms that registration named; this is the
    whole space.
    """
    tier = block.get("tier")
    if tier == "constrained_unknown":
        if block.get("copy"):
            return "constrained_unknown_copy"
        return f"real_unknown_{block.get('case') or 'uncased'}"
    return str(tier)


def elements_with_features(chrom: str) -> list[dict[str, Any]]:
    """Every scored element of one chromosome with its block case, GC and distance to a coding TSS."""
    p = ELEMENTS / f"{chrom}.json"
    if not p.exists():
        return []
    from genomeos.coords import Locus
    from genomeos.predict.enhancer_target import Context

    try:
        ctx = Context(chrom)
    except FileNotFoundError:
        return []
    coding_tss = sorted(
        (g.locus.end - 1 if g.locus.strand.value == "-" else g.locus.start)
        for g in ctx.annotation.protein_coding()
        if g.locus.chrom == chrom
    )
    cases: list[tuple[int, int, str]] = []
    for b in organise.blocks(chrom):
        cases.append((b["start"], b["end"], label_of(b)))
    cases.sort()
    starts = [c[0] for c in cases]

    rows = []
    for e in json.loads(p.read_text()):
        mid = (e["start"] + e["end"]) // 2
        i = max(0, bisect.bisect_right(starts, mid) - 1)
        arm = "other"
        while i < len(cases) and cases[i][0] <= mid:
            if cases[i][0] <= mid < cases[i][1]:
                arm = cases[i][2]
                break
            i += 1
        j = bisect.bisect_left(coding_tss, e["start"])
        near = min(
            (abs(coding_tss[k] - e["start"]) for k in (j - 1, j) if 0 <= k < len(coding_tss)), default=None
        )
        if near is None:
            continue
        seq = str(ctx.genome.fetch(Locus(chrom, e["start"], e["end"]))).upper()
        acgt = sum(seq.count(x) for x in "ACGT")
        if not acgt:
            continue
        rows.append(
            {
                "arm": arm,
                "length": e["end"] - e["start"],
                "gc": (seq.count("G") + seq.count("C")) / acgt,
                "tss": near,
                "moved": st.moved(e),
            }
        )
    ctx.close()
    return rows


def matched_reading(rows: list[dict], arm: str) -> dict[str, Any]:
    """The arm against every other scored element, standardised by genomeos.compare."""
    target = [r for r in rows if r["arm"] == arm]
    pool = [r for r in rows if r["arm"] == "other"]
    out = standardised(target, pool, strata_spec())
    return {
        "elements": len(target),
        "elements_in_a_shared_stratum": out["targets_matched"],
        "median_tss_distance": out["imbalance"]["tss"]["target_median"],
        "median_tss_distance_of_the_pool": out["imbalance"]["tss"]["control_median"],
        "raw": out["raw"],
        "matched_on_length_gc_and_tss_distance": out["matched"],
        "imbalance": out["imbalance"],
    }


def strata_spec() -> Strata:
    """This script's strata, as the shared comparator understands them."""
    return Strata(length=LEN_BINS, gc=GC_BINS, tss=TSS_BINS)


def matched_pair(rows: list[dict], a: str, b: str) -> dict[str, Any]:
    """One arm against another inside shared strata, which is what a tier-against-tier claim needs.

    The reading of 2026-09-16 put the real unknown at 0.248 against the neutral tier's 0.293 with
    nothing held fixed, and the tiers differ in distance to a promoter by an order of magnitude. This
    is the same comparison with length, GC and that distance held fixed.
    """
    out = standardised([r for r in rows if r["arm"] == a], [r for r in rows if r["arm"] == b], strata_spec())
    return {
        "raw": out["raw"],
        "matched": out["matched"],
        "left_in_a_shared_stratum": out["targets_matched"],
        "dropped_for_want_of_a_control": out["dropped_for_want_of_a_control"],
        "imbalance": out["imbalance"],
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chroms", default="")
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args(argv)
    chroms = args.chroms.split(",") if args.chroms else [f"chr{c}" for c in [*range(1, 23), "X"]]

    t0 = time.time()
    rows: list[dict] = []
    for chrom in chroms:
        got = elements_with_features(chrom)
        rows.extend(got)
        print(f"{chrom}: {len(got)} scored elements placed", flush=True)

    out = {
        "result": "syntax_blocks_matched",
        "chromosomes": chroms,
        "elements": len(rows),
        "syntax": matched_reading(rows, "real_unknown_syntax"),
        "relaxed": matched_reading(rows, "real_unknown_relaxed"),
        "neutral": matched_reading(rows, "neutral"),
        "by_label": {
            lab: matched_reading(rows, lab)
            for lab in sorted({r["arm"] for r in rows} - {"other"})
            if sum(1 for r in rows if r["arm"] == lab) >= 20
        },
        "against_neutral": {
            lab: matched_pair(rows, lab, "neutral")
            for lab in sorted({r["arm"] for r in rows} - {"other", "neutral"})
            if sum(1 for r in rows if r["arm"] == lab) >= 20
        },
        "strata": {"length": LEN_BINS, "gc": GC_BINS, "tss": TSS_BINS},
        "reading": (
            "the raw difference is the comparison that produced the observation; the matched one holds "
            "length, GC and distance to the nearest coding TSS fixed. Where they disagree, the raw one "
            "was measuring what the strata hold fixed"
        ),
        "seconds": round(time.time() - t0, 1),
    }
    if args.no_save:
        print("not saved (--no-save)")
    else:
        print(f"saved {save_result(out['result'], out)}")
    for arm in ("syntax", "relaxed", "neutral"):
        a = out[arm]
        print(
            f"{arm:8} n={a['elements']:5} median TSS {a['median_tss_distance']:>10,} "
            f"(pool {a['median_tss_distance_of_the_pool']:,}) raw {a['raw']['difference']} "
            f"matched {a['matched_on_length_gc_and_tss_distance']['difference']} "
            f"p {a['matched_on_length_gc_and_tss_distance']['p_one_sided']}"
        )
    print("\nevery label with 20 or more scored elements, against the genome's other elements:")
    for lab, a in out["by_label"].items():
        m = a["matched_on_length_gc_and_tss_distance"]
        print(
            f"  {lab:28} n={a['elements']:5} TSS {a['median_tss_distance']:>10,}  "
            f"raw {a['raw']['difference']:+.4f}  matched {m['difference']:+.4f} (p {m['p_one_sided']})"
        )
    print("\nthe same labels against the neutral tier, raw and matched:")
    for lab, pair in out["against_neutral"].items():
        print(
            f"  {lab:28} raw {pair['raw']['difference']:+.4f} -> matched "
            f"{pair['matched']['difference']:+.4f} (p {pair['matched']['p_one_sided']})"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
