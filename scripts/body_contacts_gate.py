# SPDX-License-Identifier: AGPL-3.0-or-later
"""The BioLang v0.4 §7.4 gate: lateral inhibition emerging in the Body, scored over seeds.

Grows data/demo/lateral_inhibition.bio (Collier et al. 1996 on two equivalent cells, neighbours from the
grid, the neighbours' Delta read as an amount by each cell's own network) once per seed under three
conditions: noise on, noise off, and the two daughters created in the opposite order. Pass means no
divergence without noise, exactly one anchor cell in nearly every run with it, a winner split near
50:50, and the same winner per seed whatever the creation order (no order smuggled in). The reference
is area E's direct implementation: 198 of 200 pairs diverge, first cell 51%.

    uv run python scripts/body_contacts_gate.py [RUNS]
"""

import sys
import time
from collections import Counter
from pathlib import Path

from genomeos.lang import parse, parse_file
from genomeos.results import save_result
from genomeos.runtime.body import Body

RUNS = int(sys.argv[1]) if len(sys.argv) > 1 else 200
PATH = "data/demo/lateral_inhibition.bio"
UNTIL = 40 * 60
PAIR = ("Z1.ppp", "Z4.aaa")


def outcome(module, seed):
    b = Body(module, seed=seed).run(until=UNTIL)
    types = {n: b.cells[n].cell_type for n in PAIR}
    anchors = [n for n in PAIR if types[n] == "AnchorCell"]
    d = [b.cells[n].levels.get("Dp", 0.0) for n in PAIR]
    return {
        "anchor": anchors[0] if len(anchors) == 1 else None,
        "n_anchor": len(anchors),
        "gap": abs(d[0] - d[1]),
        "types": types,
    }


t0 = time.time()
module = parse_file(PATH)
on = [outcome(module, s) for s in range(RUNS)]
swapped_src = Path(PATH).read_text().replace("daughters: Z1.ppp, Z4.aaa", "daughters: Z4.aaa, Z1.ppp")
swapped = parse(swapped_src, base_dir=None)
sw = [outcome(swapped, s) for s in range(RUNS)]
quiet = parse_file(PATH)
quiet.parameters["noise"].value = 0.0
off = [outcome(quiet, s) for s in range(RUNS)]

diverged = sum(1 for o in on if o["n_anchor"] == 1)
winners = Counter(o["anchor"] for o in on if o["anchor"])
same_as_swapped = sum(1 for a, b in zip(on, sw, strict=True) if a["anchor"] == b["anchor"])
result = {
    "runs": RUNS,
    "noise_on": {
        "exactly_one_anchor": diverged,
        "winners": dict(winners),
        "first_cell_share": round(winners.get("Z1.ppp", 0) / diverged, 3) if diverged else None,
        "both_or_neither": RUNS - diverged,
    },
    "noise_off": {
        "exactly_one_anchor": sum(1 for o in off if o["n_anchor"] == 1),
        "max_delta_gap": max(o["gap"] for o in off),
        "fates": dict(Counter("+".join(sorted(o["types"].values())) for o in off)),
    },
    "daughters_swapped": {
        "same_winner_per_seed": same_as_swapped,
        "winners": dict(Counter(o["anchor"] for o in sw if o["anchor"])),
    },
    "seconds": round(time.time() - t0, 1),
}
result["pass"] = bool(
    result["noise_off"]["exactly_one_anchor"] == 0
    and result["noise_on"]["exactly_one_anchor"] >= 0.95 * RUNS
    and 0.35 <= (result["noise_on"]["first_cell_share"] or 0) <= 0.65
    and result["daughters_swapped"]["same_winner_per_seed"] == RUNS
)
result["reference"] = (
    "area E, genomeos/organism/lateral.py: 198 of 200 diverge, first cell 51% (noise sd 0.02)"
)
result["note"] = (
    "the Body's noise is multiplicative on every network species (Euler-Maruyama, per-cell streams seeded by "
    "run seed and cell name); Collier's reference adds Gaussian noise to Delta only, so the comparison is "
    "qualitative: divergence, balance and order-independence, not the exact rates"
)
print(save_result("body_contacts_gate", result))
print(result["noise_on"], result["noise_off"], result["daughters_swapped"], result["pass"])
