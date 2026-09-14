# SPDX-License-Identifier: AGPL-3.0-or-later
"""AC/VU as an equivalence group: exactly one anchor cell, and never the same one (area E, §7.4).

In the worm, Z1.ppp and Z4.aaa are interchangeable: one becomes the anchor cell, the other a ventral
uterine precursor, and across animals it is about half and half (Kimble & Hirsh 1979; Kimble 1981;
Seydoux & Greenwald 1989). The reference lineage records one outcome, so the falsifier is explicit: a
runtime that picks the same winner every run has smuggled in an order. This script runs
data/organisms/celegans/acvu.bio once per seed under six conditions and reports the split, never one
run:

  noise_on        the program as written: exactly one anchor cell, near 50:50 between the two
  noise_off       the same with the noise parameter zeroed: two equal cells must not diverge
  daughters_swapped   each division creating its daughters in the opposite order: the same winner
                  per seed, or the creation order is deciding it
  mother_inherits the contact between Z1.ppp and Z4.aa (the group member's mother) restored in the
                  contact table: a cell that heard Delta for 56 minutes hands the answer to its
                  daughter, and the contest is over before it starts
  whole_group     LAG-2 presented by all four cells rather than the two of the equivalence group:
                  a four-cell chain alternates, so both central cells lose
  reference       what the run is worth against a reference lineage that records one outcome

    uv run python scripts/celegans_acvu.py [SEEDS]
"""

import sys
import tempfile
import time
from collections import Counter
from pathlib import Path

from genomeos.lang import parse_file
from genomeos.results import save_result
from genomeos.runtime.body import Body

SEEDS = int(sys.argv[1]) if len(sys.argv) > 1 else 200
PROGRAM = "data/organisms/celegans/acvu.bio"
TABLE = "data/organisms/celegans/contacts_acvu.tsv"
UNTIL = 4200.0
GROUP = ["Z1ppa", "Z1ppp", "Z4aaa", "Z4aap"]
PAIR = ("Z1ppp", "Z4aaa")
REFERENCE_ANCHOR = "Z1ppp"  # WormWeb names Z1.ppp `gon herm anch`; the animal does not always agree


def outcomes(module, seeds=SEEDS):
    """One run per seed; who ended up the anchor cell, and how far apart the two cells' Delta got."""
    out = []
    for seed in range(seeds):
        b = Body(module, seed=seed).run(until=UNTIL)
        anchors = [n for n in GROUP if b.cells[n].cell_type == "AnchorCell"]
        d = [b.cells[n].levels.get("Dp", 0.0) for n in PAIR]
        out.append(
            {
                "seed": seed,
                "anchor": anchors[0] if len(anchors) == 1 else None,
                "n_anchor": len(anchors),
                "gap": abs(d[0] - d[1]),
                "fates": {n: b.cells[n].cell_type for n in GROUP},
            }
        )
    return out


def summarise(runs, label):
    exactly_one = sum(1 for r in runs if r["n_anchor"] == 1)
    winners = Counter(r["anchor"] for r in runs if r["anchor"])
    return {
        "condition": label,
        "runs": len(runs),
        "exactly_one_anchor": exactly_one,
        "no_anchor": sum(1 for r in runs if r["n_anchor"] == 0),
        "more_than_one": sum(1 for r in runs if r["n_anchor"] > 1),
        "winners": dict(winners),
        "Z1ppp_share": round(winners.get("Z1ppp", 0) / exactly_one, 3) if exactly_one else None,
        "max_delta_gap": round(max(r["gap"] for r in runs), 6),
        "median_delta_gap": round(sorted(r["gap"] for r in runs)[len(runs) // 2], 6),
    }


t0 = time.time()
module = parse_file(PROGRAM)
on = outcomes(module)
result = {"program": PROGRAM, "table": TABLE, "until_min": UNTIL, "seeds": SEEDS}
result["noise_on"] = summarise(on, "the program as written")

# noise off: two equal cells with a deterministic rule must not diverge at all
quiet = parse_file(PROGRAM)
quiet.parameters["noise"].value = 0.0
off = outcomes(quiet)
result["noise_off"] = summarise(off, "noise = 0")

# creation order: the daughters of each division made in the opposite order
alt = parse_file(PROGRAM)
for d in alt.decisions:
    if d.id in ("div_Z1pp", "div_Z4aa"):
        d.daughters = list(reversed(d.daughters))
sw = outcomes(alt)
result["daughters_swapped"] = summarise(sw, "each division's daughters created in the opposite order")
result["daughters_swapped"]["same_winner_per_seed"] = sum(
    1 for a, b in zip(on, sw, strict=True) if a["anchor"] == b["anchor"]
)

# the mother's contact restored: a Delta heard before the division decides the contest
with tempfile.TemporaryDirectory() as tmp:
    biased = Path(tmp) / "contacts.tsv"
    biased.write_text(
        Path(TABLE).read_text()
        + "\n2742\tZ1ppp\tZ4aa\t1.0\n2750\tZ1ppa\tZ1ppp\t1.0\n2750\tZ1ppp\tZ4aa\t1.0\n"
    )
    mother = parse_file(PROGRAM)
    mother.organism.contacts = str(biased)
    mi = outcomes(mother, seeds=min(SEEDS, 100))
result["mother_inherits"] = summarise(mi, "Z1.ppp touching Z4.aa, the group member's mother, from 2742 min")

# the equivalence group is two cells, not four: with all four presenting LAG-2 the chain alternates
whole = parse_file(PROGRAM)
for d in whole.decisions:
    if d.id == "presents_lag2":
        d.when = {"cell": "|".join(GROUP)}
wg = outcomes(whole, seeds=min(SEEDS, 100))
result["whole_group"] = summarise(wg, "all four cells presenting LAG-2")

# what this is worth against a lineage that records one outcome
hits = sum(1 for r in on if r["anchor"] == REFERENCE_ANCHOR)
result["reference"] = {
    "reference_anchor": REFERENCE_ANCHOR,
    "source": "Sulston & Horvitz 1977 via WormWeb: Z1.pp's posterior daughter is `gon herm anch`",
    "runs_matching_the_reference": hits,
    "share": round(hits / len(on), 3),
    "note": (
        "the reference records one animal's outcome, so this is the ceiling for an honest program "
        "and 100% is the score of a program that has named the winner. The worm program "
        "(lineage_larva.bio) does name it, and takes that 100%: the AC/VU cell is the one place in "
        "the 1,092 terminal cells where the lineage score is not the right instrument"
    ),
}

result["seconds"] = round(time.time() - t0, 1)
result["pass"] = bool(
    result["noise_on"]["exactly_one_anchor"] >= 0.95 * SEEDS
    and 0.35 <= (result["noise_on"]["Z1ppp_share"] or 0) <= 0.65
    and result["noise_off"]["exactly_one_anchor"] == 0
    and result["noise_off"]["max_delta_gap"] < 1e-9
    and result["daughters_swapped"]["same_winner_per_seed"] == SEEDS
)
print(save_result("celegans_acvu", result))
for key in ("noise_on", "noise_off", "daughters_swapped", "mother_inherits", "whole_group"):
    print(key, result[key])
print("reference", result["reference"]["runs_matching_the_reference"], "/", SEEDS)
print("pass", result["pass"])
