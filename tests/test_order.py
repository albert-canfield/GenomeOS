# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""BioLang v0.4 §7.6: `order`, a sequence the runtime can be wrong about.

The language could say when a thing happens and where it is, but not that A happens before B, which
is why the known-locus benchmark could not test HOXD colinearity at all. An `order` drives nothing —
it sets no time and fires no rule — so the only thing it can do is be contradicted: by the
coordinates the program already carries (`axis: position`, checked when it compiles), or by the run
(`axis: time`, reported as the order taken and the inversions against it).

The falsifiers, and the honest split between them: the position claim is checked on HOXD, whose
coordinates we hold; the time claim is gated on the worm, whose birth times the reference gives
independently of the program. The motivating claim — that HOXD *activates* in that order — is not
testable here, because the project holds no expression time course for a human HOX cluster.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from genomeos.ir import Order
from genomeos.lang import parse, parse_file
from genomeos.lang.parser import BioLangError
from genomeos.runtime.body import Body

ROOT = Path(__file__).resolve().parent.parent
HOXD = ROOT / "data" / "demo" / "hoxd_order.bio"
WORM = ROOT / "data" / "organisms" / "celegans" / "embryo.bio"

# Sulston's founder cells as the reference orders them: AB and P1 together, then EMS and P2, then the
# four the reference puts at one instant (62 min: it does not order MS against C), then D and P4.
FOUNDERS = ["AB|P1", "EMS|P2", "MS|E|C|P3", "D|P4"]
FLAT = [m for step in FOUNDERS for m in step.split("|")]


# ---- along position, checked when it compiles ---------------------------------------


def test_the_hoxd_cluster_compiles_and_the_order_is_a_block():
    module = parse_file(HOXD)
    (order,) = module.orders()
    assert order.axis == "position" and order.direction == "decreasing"
    assert order.members[0] == "HOXD1" and order.members[-1] == "HOXD13"


def test_a_sequence_that_contradicts_the_genome_is_refused():
    text = HOXD.read_text()
    with pytest.raises(BioLangError, match="comes before"):
        parse(text.replace("HOXD1, HOXD4", "HOXD4, HOXD1"))


def test_the_overlap_that_writing_it_found_is_named_in_the_refusal():
    """HOXD3's GENCODE record spans HOXD4, so the canonical textbook sequence is not monotonic in
    transcription start. The refusal says so, because that is a question about the annotation rather
    than about the sequence, and an author needs to be told which."""
    text = HOXD.read_text()
    with pytest.raises(BioLangError, match="spans overlap"):
        parse(text.replace("HOXD1, HOXD4, HOXD3", "HOXD1, HOXD3, HOXD4"))


def test_an_order_along_position_needs_coordinates_and_one_chromosome():
    base = "module t\ngene A { locus: chr1:100-200(+) }\ngene B { locus: chr1:300-400(+) }\n"
    parse(base + "order o { members: A, B; axis: position; direction: increasing }\n")
    with pytest.raises(BioLangError, match="has no locus"):
        parse(base + "gene C { }\norder o { members: A, C; axis: position; direction: increasing }\n")
    with pytest.raises(BioLangError, match="no common order"):
        parse(
            base
            + "gene D { locus: chr2:100-200(+) }\n"
            + "order o { members: A, D; axis: position; direction: increasing }\n"
        )


def test_the_clauses_it_cannot_check_are_refused_rather_than_ignored():
    base = "module t\ngene A { locus: chr1:100-200(+) }\ngene B { locus: chr1:300-400(+) }\n"
    with pytest.raises(BioLangError, match="chromatin opens"):
        parse(base + "order o { members: A, B; axis: position; direction: opening }\n")
    with pytest.raises(BioLangError, match="direction: increasing or decreasing"):
        parse(base + "order o { members: A, B; axis: position }\n")
    with pytest.raises(BioLangError, match="at least two steps"):
        parse(base + "order o { members: A; axis: position; direction: increasing }\n")
    with pytest.raises(BioLangError, match="only 'observe: birth' is implemented"):
        parse(base + "order o { members: A, B; axis: time; observe: activation }\n")


# ---- along time, checked against a run -----------------------------------------------


def _worm_with_order(steps: list[str]):
    module = parse_file(WORM)
    groups = [step.split("|") for step in steps]
    module.add(Order(id="founders", kind="order", groups=groups, axis="time", observe="birth"))
    return module


def test_the_worm_founders_are_born_in_the_published_order():
    body = Body(_worm_with_order(FOUNDERS), seed=None).run(until=200)
    (report,) = body.orders()
    assert report["inversions"] == 0 and report["never happened"] == []
    assert sorted(report["taken"]) == sorted(FLAT) and report["ok"]


def test_an_order_the_run_takes_backwards_is_reported_and_not_hidden():
    """A run that contradicts the declared sequence is a result, not a crash: it is reported, with the
    order it actually took, the way `ambiguous_fates` is reported."""
    backwards = ["D|P4", "MS|E|C|P3", "EMS|P2", "AB|P1"]
    body = Body(_worm_with_order(backwards), seed=None).run(until=200)
    (report,) = body.orders()
    assert report["inversions"] > 0 and not report["ok"]
    assert report["taken"][0] in ("AB", "P1")  # what the run did, against what the program claimed


def test_a_member_that_never_happens_is_counted_apart_from_one_out_of_order():
    body = Body(_worm_with_order([*FOUNDERS, "ZZZ"]), seed=None).run(until=200)
    (report,) = body.orders()
    assert report["never happened"] == ["ZZZ"] and report["inversions"] == 0 and not report["ok"]


def test_an_order_drives_nothing():
    """The rule that keeps this from being a second `stage`: adding an order must not change a run."""
    plain = Body(parse_file(WORM), seed=None).run(until=200).summary()
    ordered = Body(_worm_with_order(FOUNDERS), seed=None).run(until=200).summary()
    assert {k: v for k, v in ordered.items() if k != "orders"} == {
        k: v for k, v in plain.items() if k != "orders"
    }
