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

# The design, fixed before any window was scored.
WINDOW = 300  # bp, close to the median length of an ENCODE cCRE, so the arms are comparable to the sweep
MAX_PER_BLOCK = 12  # evenly spaced, so a 40 kb block cannot outvote a 2 kb one
GC_TOL = 0.04  # the project's matching convention (human panel, unknown_scoring)
TSS_TOL = 0.35  # relative distance to the nearest coding TSS
MIN_EFFECT = 0.1  # the sweep's own threshold for "this deletion moved a gene"
ARMS = ("syntax", "relaxed", "neutral")

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
