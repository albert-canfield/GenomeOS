# SPDX-License-Identifier: AGPL-3.0-or-later
"""What choosing concentrations would cost (BIOLANG-v0.4-ECONOMY.md §10, open decision 2).

Decision 2 holds stage 2, and §9's table files it as *partly measurable*: "run stage 2's burden gate
both ways on cells of different volume; if only one reproduces the published burden scaling, it
decides." That test cannot be run, and finding out why is cheaper than running it — which is the
same shape as decision 1, where a preference turned out to have a countable price.

**The discriminating test is not constructible with what the language holds.** `volume` on a
`compartment` is documented in two places as *a fraction of the cell's volume*
(BIOLANG-v0.4-ECONOMY.md's compartment table and BIOLANG-GRAMMAR.md), and a fraction cannot express
"cells of different volume". A concentration is molecules divided by an absolute volume, and no
construct in v0.4 carries one. So the concentration arm of the gate has nothing to divide by, and the
proposed experiment is blocked on a language change rather than on data.

**What that change would cost is small, and this script counts it rather than asserting it.** The
cost falls only on programs that declare a compartment, because only they can hold an amount that a
volume would convert. The count is the result; the argument is in the result's own fields.

**Why the number will not stay this small.** It is small *because stage 2 has not been built*: pools,
costs and allocations do not exist in any committed program yet, so almost nothing is denominated in
either unit. Every program stage 2 adds raises the price of changing the answer. This is the cheapest
moment the decision will ever have, and that is an argument about timing rather than about which unit
is right.

    uv run python scripts/amounts_or_concentrations.py
"""

from __future__ import annotations

import time
from pathlib import Path

from genomeos.lang import parse_file
from genomeos.results import save_result

PROGRAM_DIRS = ("data/organisms", "data/demo", "genomeos/std")

#: absolute volumes are published for both cell types the committed compartments stand for, so the
#: concentration arm would not require an author to invent one. Recorded here as the citation a
#: program would carry, not as a value this script uses for anything.
PUBLISHED_VOLUMES = {
    "human.erythrocyte": (
        "the mature human erythrocyte, about 90 fL (Alberts MBoC 6e; BioNumbers BNID 101711)"
    ),
    "human.oxphos": (
        "the hepatocyte, whose volume FRACTIONS this program already cites from Alberts MBoC 6e"
        " Table 12-1; the absolute volume is about 3,400 um^3 and comes from the same literature"
    ),
}


def programs(root: Path) -> list[Path]:
    out: list[Path] = []
    for d in PROGRAM_DIRS:
        if (root / d).is_dir():
            out.extend(sorted((root / d).rglob("*.bio")))
    return out


def _absent(v) -> bool:
    """The parser fills a missing field with an UNKNOWN sentinel, not None.

    Written after this script first reported a confident zero: iterating `module.entities` yields
    KEYS, because it is a dict, and `volume` reads `UNKNOWN` rather than `None` when unstated. Both
    produced a plausible count of nothing. A zero from reading nothing and a zero from measuring
    nothing must not be the same output, so the guard below asserts the read found something.
    """
    return v is None or str(v) == "UNKNOWN"


def read(path: Path) -> dict:
    """One program's compartments and whether anything is denominated in a unit at all."""
    module = parse_file(str(path))
    entities = list(module.entities.values())
    compartments = [e for e in entities if type(e).__name__ == "Compartment"]
    with_volume = [c for c in compartments if not _absent(getattr(c, "volume", None))]
    return {
        "path": str(path),
        "module": module.name,
        "compartments": len(compartments),
        "compartments_with_a_volume_fraction": len(with_volume),
        "volume_fractions": {c.id: c.volume for c in with_volume},
        # stage 2's constructs, which are what a unit choice would be baked into
        "entities": len(entities),
        "pools": sum(1 for e in entities if type(e).__name__ == "Pool"),
        "costs": sum(1 for e in entities if not _absent(getattr(e, "costs", None))),
    }


def main() -> int:
    t0 = time.time()
    root = Path(".")
    rows = [read(p) for p in programs(root)]
    # the guard the first run needed: a zero here must mean nothing declares a compartment, not that
    # the read failed. If no program parsed to any entity at all, say so instead of counting zero.
    if rows and not any(r["entities"] for r in rows):
        raise SystemExit("read 0 entities from every program: the reader is broken, not the answer")
    located = [r for r in rows if r["compartments"]]
    volumed = [r for r in located if r["compartments_with_a_volume_fraction"]]
    cited = [r for r in volumed if r["module"] in PUBLISHED_VOLUMES]

    out = {
        "result": "amounts_or_concentrations",
        "decision": "BIOLANG-v0.4-ECONOMY.md section 10, open decision 2",
        "why_the_registered_test_cannot_be_run": (
            "section 9 proposes running stage 2's burden gate both ways on cells of different volume."
            " `volume` on a compartment is documented in two places as a FRACTION of the cell's"
            " volume, and a fraction cannot express cells of different volume. A concentration is"
            " molecules over an absolute volume and v0.4 has no construct carrying one, so the"
            " concentration arm has nothing to divide by. The test is blocked on a language change,"
            " not on data"
        ),
        "programs_read": len(rows),
        "programs_declaring_a_compartment": len(located),
        "programs_with_a_volume_fraction": len(volumed),
        "cell_types_needing_an_absolute_volume": len(volumed),
        "of_those_with_a_published_source": len(cited),
        "published_volumes": PUBLISHED_VOLUMES,
        "stage_two_constructs_committed": {
            "pools": sum(r["pools"] for r in rows),
            "entities_with_a_cost": sum(r["costs"] for r in rows),
        },
        "the_price": (
            "one language field (an absolute volume, which no construct carries) and one citation per"
            " cell type that declares a compartment. Both cell types here have a published absolute"
            " volume, so neither would be invented, and section 2's penalty for an uncited fact does"
            " not apply. That is the whole cost of making concentrations expressible today"
        ),
        "why_it_does_not_stay_this_small": (
            "the price is low because stage 2 is not built: no committed program declares a pool or a"
            " cost, so almost nothing is denominated in either unit yet. Every program stage 2 adds"
            " raises the cost of changing the answer, so this is the cheapest moment the decision"
            " will have. That is an argument about timing, not about which unit is right"
        ),
        "what_this_does_not_decide": (
            "nothing about which unit is correct. Amounts and concentrations are not distinguished by"
            " this count, and the burden gate still cannot separate them until the language can"
            " express an absolute volume. What the count says is that the blocker is one field rather"
            " than a data problem, and that the objection recorded against concentrations - that"
            " volumes vary per cell type - costs two citations at today's scale"
        ),
        "per_program": located,
        "seconds": round(time.time() - t0, 2),
    }
    print(f"programs read: {out['programs_read']}")
    print(f"declaring a compartment: {out['programs_declaring_a_compartment']}")
    print(f"with a volume fraction: {out['programs_with_a_volume_fraction']}")
    print(f"of those with a published absolute volume: {out['of_those_with_a_published_source']}")
    print(f"pools committed: {out['stage_two_constructs_committed']['pools']}")
    print(f"\nsaved {save_result(out['result'], out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
