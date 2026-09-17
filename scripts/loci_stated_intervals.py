# SPDX-License-Identifier: AGPL-3.0-or-later
"""Delete the published element itself at the loci that stated one, which no annotation drew.

    uv run python scripts/loci_stated_intervals.py                  # what it would ask, and why
    uv run python scripts/loci_stated_intervals.py --run --quota-handed

docs/LOCI-BENCHMARK.md 11a: `Expect.element` is already a stated interval, and
`Context.score_region` already scores a region as an ad-hoc element when no ENCODE element overlaps
it. The benchmark's deletion layer reads `loci_stated_intervals` first and labels every row
`stated_interval`, so these can never be pooled with registry-derived rates. All that was missing
was the driver, and the key being free.

**Only one of the two stated intervals is worth a request, and that is this run's first result.**
The scorer resizes its input to 1 Mb around the element. SOX9 is 1,450 kb from the Pierre Robin
element, so it is not in the model's input and no number of requests can make the deletion name it;
`read_reach` now says so, and the locus is marked unaskable rather than pending. The same arithmetic
retires a question the benchmark has been grading for months: the ZRS's own deletion named LMBR1,
the gene it sits inside, and that was read as a miss against SHH, which is 979 kb away and was never
a candidate. H19 is the locus where the question survives: the ICR is 3 kb from H19 and 138 kb from
IGF2, both well inside the window.

**The prediction, written here before the request.** Deleting the ICR will name a gene - naming some
gene is not rare, it happens at 92.9% of the panel's matched negative windows - and the question is
which. Three outcomes were fixed before scoring:

  hit      IGF2 or H19 first. The ICR is the CTCF cluster that decides which of the two is read, so
           this is the published answer and would be the first element-level long-range hit the
           panel has.
  trap     MRPL23 first, the nearest coding TSS. That is the heuristic's answer, and the panel keeps
           the trap named in the expectation so the miss can be recognised rather than explained.
  neither  some third gene. The ZRS precedent says this is the likely one: an element-level deletion
           tends to name what it sits in or beside.

The honest prior is `trap` or `neither`. The ICR's effect is parent-of-origin and methylation-borne,
and the benchmark already records that no layer in the project carries parent of origin, so a
sequence model reading one allele has no way to express what the ICR does. A hit would therefore be
the surprise, and a miss is not evidence against the model so much as against the question.
"""

from __future__ import annotations

import argparse
import time
from contextlib import contextmanager
from typing import Any

from genomeos import jobs
from genomeos.benchmark.loci import PANEL, Chromosome, read_reach
from genomeos.predict import AlphaGenomeAdapter, status
from genomeos.predict.enhancer_target import Context
from genomeos.results import save_result

HOLDER = "genomeos-9c"
NAME = "loci_stated_intervals"

#: fixed before the first request; the reading below applies it without interpretation
PREDICTION = {
    "written": "2026-09-17, before any stated interval was scored",
    "expected": "trap or neither, not a hit",
    "why": (
        "the ICR acts through parent-of-origin methylation, which no layer in this project carries"
        " and which a sequence model reading one allele cannot express; and the ZRS precedent is"
        " that an element-level deletion names what it sits in or beside"
    ),
    "outcomes": {
        "hit": "IGF2 or H19 named first",
        "trap": "MRPL23 named first, the nearest coding TSS",
        "neither": "some third gene named first",
    },
}


def say(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


@contextmanager
def the_key(what: str):
    jobs.take_key(HOLDER, what)
    try:
        yield
    finally:
        jobs.drop_key(HOLDER)


def stated() -> list:
    """The panel's stated intervals, in panel order."""
    return [e for e in PANEL if e.element_source == "stated"]


def plan() -> list[dict[str, Any]]:
    """Which stated intervals can be asked, and which are unaskable arithmetic. No request."""
    rows = []
    for expect in stated():
        ch = Chromosome(expect.chrom)
        try:
            reach = read_reach(ch, expect)
            rows.append(
                {
                    "locus": expect.locus,
                    "chrom": expect.chrom,
                    "start": expect.element[0],
                    "end": expect.element[1],
                    "askable": reach["askable"],
                    "targets_in_reach": reach["targets_in_reach"],
                    "targets_out_of_reach": reach["targets_out_of_reach"],
                    "trap": expect.nearest_gene_trap,
                    "targets": list(expect.targets),
                }
            )
        finally:
            ch.close()
    return rows


def verdict(named: list[str], targets: list[str], trap: str | None) -> str:
    """The three registered outcomes, applied without interpretation."""
    if not named:
        return "no gene named: the deletion moved nothing above the floor"
    first = named[0]
    if first in targets:
        return f"hit: {first} named first, and it is published"
    if trap and first == trap:
        return f"trap: {first} named first, which is the nearest coding TSS"
    return f"neither: {first} named first, neither a published target nor the nearest-TSS trap"


def run(rows: list[dict[str, Any]]) -> dict[str, Any]:
    t0 = time.time()
    adapter = AlphaGenomeAdapter()
    scorer = adapter._live_scorer(threshold=0.0)  # noqa: SLF001  (as scripts/loci_score.py does)
    elements: list[dict[str, Any]] = []
    read: list[dict[str, Any]] = []
    asked = [r for r in rows if r["askable"]]
    with the_key(f"known-locus stated intervals, {len(asked)} elements"):
        for row in asked:
            ctx = Context(row["chrom"])
            try:
                say(f"{row['locus']}: deleting {row['chrom']}:{row['start']}-{row['end']}")
                hit = ctx.score_region(scorer, row["start"], row["end"])
            finally:
                ctx.close()
            # score_region substitutes an overlapping ENCODE element when there is one, which would
            # silently answer about a different interval than the published one. Recorded, not hidden.
            substituted = hit.get("id") != f"{row['chrom']}_{row['start']}_{row['end']}"
            elements.append({**hit, "locus": row["locus"], "substituted_annotated_element": substituted})
            named = [g["gene"] for g in (hit.get("top_genes") or []) if g.get("gene")]
            pred = hit.get("predicted_coding") or hit.get("predicted") or {}
            first = [pred["gene"]] if pred.get("gene") else []
            read.append(
                {
                    "locus": row["locus"],
                    "element": f"{row['chrom']}:{row['start']}-{row['end']}",
                    "substituted_annotated_element": substituted,
                    "named_first": pred.get("gene"),
                    "log2_fold_change": pred.get("log2_fold_change"),
                    "tissue": pred.get("tissue"),
                    "action": pred.get("action"),
                    "genes_in_window": hit.get("genes_in_window"),
                    "top_genes": named[:6],
                    "verdict": verdict(first, row["targets"], row["trap"]),
                }
            )
            say(read[-1]["verdict"])
    return {
        "result": NAME,
        "prediction": PREDICTION,
        "elements": elements,
        "read": read,
        "unaskable": [
            {
                "locus": r["locus"],
                "why": (
                    f"{r['targets_out_of_reach'][0]['target']} is"
                    f" {r['targets_out_of_reach'][0]['distance_bp'] // 1000} kb from the element and the"
                    " model's input reaches 524 kb each way: no request can name it"
                ),
            }
            for r in rows
            if not r["askable"] and r["targets_out_of_reach"]
        ],
        "note": (
            "AlphaGenome deletions of intervals the panel STATED from a publication, which no"
            " registry drew an element over (docs/LOCI-BENCHMARK.md 11a). Read first and labelled"
            " `stated_interval` by the benchmark so these rows are never pooled with registry ones"
        ),
        "requests_this_run": len(read),
        "seconds": round(time.time() - t0, 1),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--quota-handed", action="store_true")
    args = ap.parse_args(argv)

    rows = plan()
    for r in rows:
        if r["askable"]:
            say(f"{r['locus']}: askable, targets in reach {r['targets_in_reach']}, trap {r['trap']}")
        else:
            far = r["targets_out_of_reach"][0]
            say(
                f"{r['locus']}: UNASKABLE, {far['target']} is {far['distance_bp'] // 1000} kb away and"
                " the model reaches 524 kb. No request is worth spending here"
            )
    if not args.run:
        say(f"{sum(1 for r in rows if r['askable'])} of {len(rows)} would be asked; pass --run")
        return 0
    if not args.quota_handed:
        raise SystemExit("this spends model requests: pass --quota-handed once the key is yours")
    st = status()
    if not st["enabled"]:
        raise SystemExit(f"AlphaGenome is disabled: {st['reason']}. {st['how']}")
    out = run(rows)
    say(f"saved {save_result(NAME, out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
