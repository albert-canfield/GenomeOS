# SPDX-License-Identifier: AGPL-3.0-or-later
"""The 69 syntax blocks, tiled and deleted against two matched control arms.

On 2026-09-16 the finished sweep was read over the 882 blocks of the real unknown, and the tier came
out as the least active in the genome per element: 0.248 against 0.293 at the neutral tier, P 7e-9
below it. One slice read the other way — the `syntax` case, the 69 blocks held across mammals *and*
constrained among people — at 0.575 of 40 elements, above neutral at P 0.0002. Forty elements, sliced
after the fact, against a tier average rather than matched sequence.

This module is the instrument that can fail that observation. It does not use the sweep's elements: it
tiles the blocks themselves in fixed windows and deletes each one, so the sample is set by the blocks'
length rather than by where ENCODE happened to call a cCRE, and it scores two matched control arms at
the same time.

The pre-registration below was written and committed before the instrument existed and before any
request was made. `PRE_REGISTRATION` is the machine-readable copy; the prose is docs/ATTRIBUTION.md.
"""

from __future__ import annotations

from typing import Any

# re-exported deliberately: `scripts/syntax_tiling.py` calls `st.difference`, and this module used to
# define its own copy of it. `as difference` says the name is part of this module's surface rather
# than an unused import.
from genomeos.compare import difference as difference

# the sweep's own threshold, imported rather than restated: this module's registration text
# interpolates it, so a divergence would silently change what the registration said
from genomeos.predict.enhancer_target import MIN_EFFECT

# The design, fixed before any window was scored.
WINDOW = 300  # bp, close to the median length of an ENCODE cCRE, so the arms are comparable to the sweep
MAX_PER_BLOCK = 12  # evenly spaced, so a 40 kb block cannot outvote a 2 kb one
GC_TOL = 0.04  # the project's matching convention (human panel, unknown_scoring)
TSS_TOL = 0.35  # relative distance to the nearest coding TSS
ARMS = ("syntax", "relaxed", "neutral")

# NOT RUN. Assembling the arms showed they do not overlap in the covariate that decides the outcome:
# the syntax windows sit a median 81 kb from a coding TSS, the relaxed arm 240 kb and the neutral arm
# 418 kb, against 42 kb for the genome's scored elements. The cheaper question - are the 40 elements
# already scored inside these blocks different from matched elements? - answered it for no requests:
# -0.033 at p 0.67 against the genome's own elements. The registration stands as written and unrun;
# see docs/ATTRIBUTION.md, "The syntax observation dissolves under matching" (2026-09-17).
UNRUN = "cancelled 2026-09-17 after scripts/syntax_blocks_matched.py answered the question for free"

PRE_REGISTRATION: dict[str, Any] = {
    "written": (
        "2026-09-17, before the instrument was written and before any window was scored; committed "
        "in this file before the code that runs it"
    ),
    "question": (
        "do deletions inside the 69 syntax blocks move a gene more often than deletions in matched "
        "sequence? The observation that prompts it (0.575 of 40 elements) rests on elements ENCODE "
        "happened to call inside those blocks, compared with a tier average rather than with matched "
        "windows, and was sliced after the fact"
    ),
    "arms": {
        "syntax": (
            "the 69 constrained_unknown, non-copy blocks whose case is syntax: held across mammals "
            "and constrained among people"
        ),
        "relaxed": (
            "the 437 blocks whose case is relaxed: held across mammals, variable among people. This arm "
            "isolates the human axis, since it differs from the syntax arm in that alone"
        ),
        "neutral": "blocks of the neutral tier: the background the whole tier reading was taken against",
    },
    "sampling": (
        f"every block of an arm is tiled in non-overlapping {WINDOW} bp windows from its start; where a "
        f"block yields more than {MAX_PER_BLOCK}, {MAX_PER_BLOCK} are taken evenly spaced across it, so "
        "block length cannot decide the sample. Windows are scored one request each, as the sweep scores "
        "an element"
    ),
    "matching": (
        f"each syntax window is matched, on its own chromosome, to one relaxed and one neutral window "
        f"with GC within {GC_TOL} and distance to the nearest coding TSS within {TSS_TOL} relative. "
        "Length is equal by construction. The arms are NOT matched on constraint: mammalian constraint "
        "is what defines the tier and human constraint is what the relaxed arm isolates, so matching on "
        "either would remove the variable under test. A syntax window with no match in an arm is dropped "
        "from that comparison and counted"
    ),
    "outcome": (
        f"a window moves a gene if deleting it changes a gene's predicted expression by at least "
        f"{MIN_EFFECT} in |log2| and names that gene, which is the sweep's own threshold"
    ),
    "primary": (
        "the difference in the share of windows that move a gene, syntax minus neutral and syntax minus "
        "relaxed, each with a one-sided p and a 95% upper bound"
    ),
    "outcomes": {
        "success": "at least +0.10 against BOTH control arms at one-sided p 0.01 or better",
        "weak": (
            "0.05 to 0.10 against both, or at least +0.10 against one arm only: the observation survives "
            "in a form worth one more instrument, and nothing is claimed for the blocks"
        ),
        "failure": (
            "below 0.05 against the neutral arm: the 0.575 of 40 elements was the sample it came from, "
            "and the syntax case stops being carried as the sharpest candidate in the roadmap"
        ),
    },
    "budget": "2,400 requests, one per window, about 800 per arm",
    "stopping": (
        "one look at half the budget, for futility only: if the 95% upper bound of syntax minus neutral "
        "is below 0.05 the run stops and is reported as stopped for futility. No interim success look"
    ),
    "secondary": (
        "declared now: inside the syntax arm, windows that overlap an element the sweep already scored "
        "should act about as often as those 40 did (0.575), and windows that overlap none are the new "
        "information. If only the overlapping windows act, the tiling has added nothing and the arm's "
        "rate is the old observation in a new coat"
    ),
    "cannot_do": (
        "this is the same model that produced the observation, so it cannot be independent evidence "
        "about it. What it can fix is the sampling: 40 elements chosen by where ENCODE called a cCRE, "
        "compared with a tier average, become a few hundred windows fixed by the blocks themselves and "
        "compared with matched sequence. A positive says the model reacts more to deleting these "
        "sequences than to deleting matched sequence; it does not say the sequence is functional, and "
        "measured evidence on these blocks (lentiMPRA where it covers them) remains the thing that would"
    ),
}


# --- the instrument, written after the registration above was committed ---------------------------

import bisect  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
from pathlib import Path  # noqa: E402

from genomeos.attribution import organise  # noqa: E402

REFERENCE = Path("data/reference")
ELEMENTS = Path("data/knowledge/alphagenome/all_elements")


def tile(start: int, end: int, window: int = WINDOW, cap: int = MAX_PER_BLOCK) -> list[tuple[int, int]]:
    """Non-overlapping windows across a block, at most `cap` of them, evenly spaced.

    Evenly spaced rather than the first `cap`: taking the first twelve of a 40 kb block would sample
    one end of it, and the block is not claimed to be uniform.
    """
    whole = [(s, s + window) for s in range(start, max(start, end - window + 1), window)]
    if len(whole) <= cap:
        return whole
    step = len(whole) / cap
    return [whole[int(i * step)] for i in range(cap)]


def arm_blocks(chrom: str) -> dict[str, list[dict[str, Any]]]:
    """The three arms' blocks on one chromosome, by the organiser's own reading."""
    out: dict[str, list[dict[str, Any]]] = {a: [] for a in ARMS}
    for b in organise.blocks(chrom):
        if b.get("tier") == "neutral":
            out["neutral"].append(b)
        elif (
            b.get("tier") == "constrained_unknown"
            and not b.get("copy")
            and b.get("case") in ("syntax", "relaxed")
        ):
            out[b["case"]].append(b)
    return out


def scored_starts(chrom: str) -> list[tuple[int, int]]:
    """Where the sweep already scored an element, for the declared secondary."""
    p = ELEMENTS / f"{chrom}.json"
    if not p.exists():
        return []
    return sorted((e["start"], e["end"]) for e in json.loads(p.read_text()))


def _overlaps(spans: list[tuple[int, int]], start: int, end: int) -> bool:
    i = bisect.bisect_right(spans, (start, math.inf)) - 1
    return any(0 <= j < len(spans) and spans[j][0] < end and start < spans[j][1] for j in (i, i + 1))


def windows_of(chrom: str, ctx: Any, coding_tss: list[int], spans: list[tuple[int, int]]) -> list[dict]:
    """Every arm's windows on one chromosome, with the two features they are matched on."""
    rows = []
    for arm, blocks in arm_blocks(chrom).items():
        for b in blocks:
            for start, end in tile(b["start"], b["end"]):
                from genomeos.coords import Locus

                seq = str(ctx.genome.fetch(Locus(chrom, start, end))).upper()
                acgt = sum(seq.count(x) for x in "ACGT")
                if acgt < (end - start) * 0.9:  # an assembly gap is not sequence to delete
                    continue
                i = bisect.bisect_left(coding_tss, start)
                near = min(
                    (abs(coding_tss[j] - start) for j in (i - 1, i) if 0 <= j < len(coding_tss)),
                    default=None,
                )
                if near is None:
                    continue
                rows.append(
                    {
                        "arm": arm,
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "block": f"{chrom}:{b['start']}-{b['end']}",
                        "gc": round((seq.count("G") + seq.count("C")) / acgt, 4),
                        "tss": near,
                        "already_scored": _overlaps(spans, start, end),
                    }
                )
    return rows


def match_windows(rows: list[dict]) -> list[dict]:
    """Each syntax window paired with one relaxed and one neutral window, unmatched ones counted.

    Greedy and without reuse: a control window serves one syntax window, so a single well-placed
    control cannot stand in for twenty.
    """
    pools = {a: [r for r in rows if r["arm"] == a] for a in ARMS}
    used: set[int] = set()
    kept: list[dict] = []
    for target in pools["syntax"]:
        picked = {}
        for arm in ("relaxed", "neutral"):
            best, best_d = None, None
            for cand in pools[arm]:
                if id(cand) in used:
                    continue
                dgc = abs(cand["gc"] - target["gc"])
                dtss = abs(cand["tss"] - target["tss"]) / max(target["tss"], 1)
                if dgc > GC_TOL or dtss > TSS_TOL:
                    continue
                d = (dgc / GC_TOL) + (dtss / TSS_TOL)
                if best_d is None or d < best_d:
                    best, best_d = cand, d
            if best is not None:
                used.add(id(best))
                picked[arm] = best
        target["matched"] = sorted(picked)
        kept.append(target)
        kept.extend(picked.values())
    return kept


# `difference` was a line-for-line copy of `compare.difference` and is now that function. It kept the
# old z of 1.96 under the name `upper_95`, which is the two-sided 95% constant and therefore a
# one-sided 97.5% bound beside a one-sided p; the shared helper was corrected to 1.645 on 2026-09-18
# and this copy would have gone on disagreeing with it. Nothing published from this lane - the
# registered run was cancelled unrun, so it has a plan file and no result - and the futility look in
# `scripts/syntax_tiling.py` reads `upper_95`, which is why the copy mattered even unrun.


def moved(row: dict) -> bool:
    """Did this window's deletion move a gene, by the registered threshold?"""
    pred = row.get("predicted") or {}
    log2 = pred.get("log2_fold_change")
    return bool(pred.get("gene")) and log2 is not None and abs(log2) >= MIN_EFFECT
