# SPDX-License-Identifier: AGPL-3.0-or-later
"""The oligo library that would measure the real unknown, written out rather than proposed.

    uv run python scripts/unknown_library.py [--step 300] [--max-oligos 120000] [--no-save]

`scripts/unknown_coverage.py` established that 0.45% of the unknown space has ever been measured, and
that the part this project most wants measured — the 882 blocks of the real unknown, 30.6 Mb — is
small enough for one MPRA library. This emits that library: every oligo with its coordinates, the
block it tiles, and the arm it belongs to, plus the controls without which the result would be
unreadable.

Four arms, because a reporter library with no controls measures the library:

- **test**: the real unknown, untouched blocks first, tiled end to end;
- **genomic negative**: neutral-tier blocks matched to the test oligos by GC and distance to the
  nearest coding TSS, the two covariates that decided every comparison this project got wrong before
  it standardised on them;
- **positive**: elements lentiMPRA already measured as active, carried so the library calibrates
  against a known answer in the same batch;
- **scrambled**: each test oligo's own dinucleotide-preserving shuffle, which holds composition fixed
  and destroys arrangement — the null the motif lane used, here as a per-oligo partner.

It writes `data/results/unknown_library.json` (the design and its counts) and a manifest TSV under
`data/knowledge/library/`. Sequence is fetched from the local reference; nothing is sent anywhere.
The manifest is a design, not an order: no synthesis, no vendor, no claim that any of it is functional.
"""

from __future__ import annotations

import argparse
import bisect
import random
import time
from collections import defaultdict
from pathlib import Path
from statistics import median
from typing import Any

from genomeos.attribution import mpra, organise
from genomeos.coords import Locus
from genomeos.predict.enhancer_target import Context
from genomeos.results import save_result
from scripts.unknown_coverage import label_of, measured_spans, overlap_bp  # noqa: E402

OUT = Path("data/knowledge/library")
OLIGO = 300
SEED = 17


def dinucleotide_shuffle(seq: str, rng: random.Random) -> str:
    """Altschul-Erikson: keep every dinucleotide count, destroy the arrangement.

    A mononucleotide shuffle would leave CpG content free to change, and CpG is what several of this
    project's readings turned on, so the shuffle has to hold pairs rather than bases.
    """
    if len(seq) < 3:
        return seq
    edges: dict[str, list[str]] = defaultdict(list)
    for a, b in zip(seq, seq[1:], strict=False):
        edges[a].append(b)
    for v in edges.values():
        rng.shuffle(v)
    out = [seq[0]]
    cursor = dict.fromkeys(edges, 0)
    for _ in range(len(seq) - 1):
        here = out[-1]
        i = cursor.get(here, 0)
        if here not in edges or i >= len(edges[here]):
            break
        cursor[here] = i + 1
        out.append(edges[here][i])
    return "".join(out)


def tss_of(ctx: Context, chrom: str) -> list[int]:
    return sorted(
        (g.locus.end - 1 if g.locus.strand.value == "-" else g.locus.start)
        for g in ctx.annotation.protein_coding()
        if g.locus.chrom == chrom
    )


def oligos_of_block(block: dict, step: int) -> list[tuple[int, int]]:
    return [
        (s, s + OLIGO) for s in range(block["start"], max(block["start"], block["end"] - OLIGO + 1), step)
    ]


def build(chroms: list[str], step: int, max_oligos: int) -> dict[str, Any]:
    """Every arm of the library, in the order the budget should be spent."""
    t0 = time.time()
    rng = random.Random(SEED)
    test: list[dict] = []
    pool: list[dict] = []  # neutral-tier oligos, the genomic negatives to match from
    positives: list[dict] = []

    for chrom in chroms:
        try:
            ctx = Context(chrom)
        except FileNotFoundError:
            continue
        coding_tss = tss_of(ctx, chrom)
        spans = measured_spans(chrom)
        active = {(e.start, e.end) for e in mpra.load(chrom) if max(e.activity.values(), default=0) >= 1.0}

        def feature(
            start: int,
            end: int,
            chrom: str = chrom,
            coding_tss: list[int] = coding_tss,
            ctx: Context = ctx,  # bound now: ctx is rebound each chromosome and the closure outlives it
        ):
            seq = str(ctx.genome.fetch(Locus(chrom, start, end))).upper()
            acgt = sum(seq.count(x) for x in "ACGT")
            if acgt < (end - start) * 0.9:
                return None
            i = bisect.bisect_left(coding_tss, start)
            near = min(
                (abs(coding_tss[j] - start) for j in (i - 1, i) if 0 <= j < len(coding_tss)), default=None
            )
            if near is None:
                return None
            return {
                "chrom": chrom,
                "start": start,
                "end": end,
                "gc": round((seq.count("G") + seq.count("C")) / acgt, 4),
                "tss": near,
                "seq": seq,
            }

        for b in organise.blocks(chrom):
            label = label_of(b)
            if label.startswith("real_unknown_"):
                untouched = not any(overlap_bp(v, b["start"], b["end"]) for v in spans.values())
                for start, end in oligos_of_block(b, step):
                    f = feature(start, end)
                    if f:
                        test.append(
                            {
                                **f,
                                "arm": "test",
                                "block": f"{chrom}:{b['start']}-{b['end']}",
                                "case": b.get("case"),
                                "untouched_block": untouched,
                            }
                        )
            elif label == "neutral":
                # the pool is sampled at the same step as the test arm: at a coarser step it ran out and
                # only 54% of test oligos found a matched negative, which is a library with a hole in it
                for start, end in oligos_of_block(b, step):
                    f = feature(start, end)
                    if f:
                        pool.append(
                            {**f, "arm": "genomic_negative", "block": f"{chrom}:{b['start']}-{b['end']}"}
                        )
        for start, end in sorted(active):
            f = feature(start, end if end - start >= OLIGO else start + OLIGO)
            if f:
                positives.append({**f, "arm": "positive", "block": "lentimpra_active"})
        ctx.close()

    # untouched blocks first: the library's purpose is the sequence nothing has measured
    test.sort(key=lambda r: (not r["untouched_block"], r["chrom"], r["start"]))
    test = test[:max_oligos]

    by_key: dict[tuple, list[dict]] = defaultdict(list)
    for r in pool:
        by_key[(round(r["gc"], 1), len(str(r["tss"])))].append(r)
    negatives, used = [], set()
    for r in test:
        bucket = by_key.get((round(r["gc"], 1), len(str(r["tss"]))), [])
        pick = next((c for c in bucket if id(c) not in used), None)
        if pick:
            used.add(id(pick))
            negatives.append(pick)

    scrambled = [
        {
            **r,
            "arm": "scrambled",
            "seq": dinucleotide_shuffle(r["seq"], rng),
            "of": f"{r['chrom']}:{r['start']}",
        }
        for r in test
    ]
    return {
        "test": test,
        "genomic_negative": negatives,
        "positive": positives[: len(test) // 10],
        "scrambled": scrambled,
        "seconds": round(time.time() - t0, 1),
    }


def write_manifest(arms: dict[str, list[dict]]) -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "unknown_library.tsv"
    with open(path, "w") as fh:
        fh.write("id\tarm\tchrom\tstart\tend\tblock\tgc\ttss_distance\tsequence\n")
        n = 0
        for arm in ("test", "genomic_negative", "positive", "scrambled"):
            for r in arms.get(arm, []):
                n += 1
                fh.write(
                    f"olig{n:06d}\t{arm}\t{r['chrom']}\t{r['start']}\t{r['end']}\t{r.get('block', '')}\t"
                    f"{r['gc']}\t{r['tss']}\t{r['seq']}\n"
                )
    return path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chroms", default="")
    ap.add_argument("--step", type=int, default=OLIGO, help="tiling step; OLIGO means end to end")
    ap.add_argument("--max-oligos", type=int, default=120_000, help="cap on the test arm")
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args(argv)
    chroms = args.chroms.split(",") if args.chroms else [f"chr{c}" for c in [*range(1, 23), "X"]]

    arms = build(chroms, args.step, args.max_oligos)
    seconds = arms.pop("seconds")
    path = None if args.no_save else write_manifest(arms)
    test = arms["test"]
    out = {
        "result": "unknown_library",
        "chromosomes": chroms,
        "oligo_bp": OLIGO,
        "step": args.step,
        "arms": {k: len(v) for k, v in arms.items()},
        "total_oligos": sum(len(v) for v in arms.values()),
        "test_from_untouched_blocks": sum(1 for r in test if r["untouched_block"]),
        "blocks_tiled": len({r["block"] for r in test}),
        "by_case": {
            case: sum(1 for r in test if r["case"] == case)
            for case in sorted({r["case"] for r in test if r["case"]})
        },
        "median_gc": {arm: round(median([r["gc"] for r in v]), 4) for arm, v in arms.items() if v},
        "median_tss_distance": {arm: round(median([r["tss"] for r in v]), 1) for arm, v in arms.items() if v},
        "manifest": str(path) if path else None,
        "design": (
            "four arms: the real unknown tiled end to end with its untouched blocks first; neutral-tier "
            "oligos matched on GC and order of magnitude of TSS distance as genomic negatives; lentiMPRA "
            "actives as in-batch positives; and each test oligo's own dinucleotide-preserving shuffle. A "
            "design, not an order: no synthesis, no vendor, and no claim that any of this sequence is "
            "functional"
        ),
        "seconds": seconds,
    }
    if not args.no_save:
        print(f"saved {save_result(out['result'], out)}")
    for arm, n in out["arms"].items():
        print(
            f"{arm:18} {n:7,}  median GC {out['median_gc'].get(arm)}  "
            f"median TSS {out['median_tss_distance'].get(arm)}"
        )
    print(
        f"total {out['total_oligos']:,} oligos over {out['blocks_tiled']:,} blocks; "
        f"manifest {out['manifest']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
