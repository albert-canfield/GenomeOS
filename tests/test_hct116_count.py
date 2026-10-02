# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-hct: the blinded HCT116 count, the pooled gate and the proofs that reading a direction is not
merely avoided here but structurally impossible.

The three proofs are `test_a_poison_value_...`, `test_the_output_is_identical_...` and
`test_the_reader_contains_no_value_read_...`. Together they say: any operation on a cached value other
than an identity test against None raises; a negative, a positive and an exactly-zero value produce
identical output; and the function's own syntax tree carries no comparison, arithmetic or conversion
that could read one.
"""

from __future__ import annotations

import ast
import inspect
import json
from pathlib import Path
from typing import Any

import pytest

from genomeos.attribution import cell2, fresh
from genomeos.attribution import direction_link as dl
from genomeos.attribution import hct116_count as hc

# ---- the poison value ----------------------------------------------------------------------------


class Poison:
    """A stand-in for a cached predicted value that detonates on everything but an identity test.

    `is` and `is not` are the only operations Python performs on an object without consulting it, so a
    reader that only asks `value is None` passes, and a reader that compares it, negates it, converts
    it, formats it, hashes it or looks at an attribute of it raises `hc.BlindingError`. This is the
    structural proof: the blinding is enforced by the value, not by the reader's good manners.
    """

    __slots__ = ()

    def _boom(self, *_: Any, **__: Any) -> Any:
        raise hc.BlindingError("a cached predicted value was read: this lane may only test it against None")

    __eq__ = __ne__ = _boom
    __lt__ = __le__ = __gt__ = __ge__ = _boom
    __add__ = __radd__ = __sub__ = __rsub__ = __mul__ = __rmul__ = _boom
    __truediv__ = __rtruediv__ = __floordiv__ = __mod__ = __pow__ = _boom
    __neg__ = __pos__ = __abs__ = __invert__ = _boom
    __bool__ = __int__ = __float__ = __index__ = __complex__ = _boom
    __round__ = __trunc__ = __floor__ = __ceil__ = _boom
    __str__ = __repr__ = __format__ = __hash__ = _boom
    __len__ = __iter__ = __contains__ = __getitem__ = _boom
    __getattr__ = _boom


def test_the_poison_detonates_on_everything_but_an_identity_test() -> None:
    """The instrument itself, checked first: an identity test passes, and nothing else does."""
    p = Poison()
    assert (p is None) is False
    assert (p is not None) is True
    for call in (
        lambda: p == 0.0,
        lambda: p < 0,
        lambda: p > 0,
        lambda: -p,
        lambda: abs(p),
        lambda: float(p),
        lambda: bool(p),
        lambda: str(p),
        lambda: f"{p}",
        lambda: hash(p),
        lambda: p.real,
        lambda: round(p, 4),
    ):
        with pytest.raises(hc.BlindingError):
            call()


# ---- proof 1: the value cannot be touched --------------------------------------------------------


def _record(value: Any, gene: str = "MYC", cell: str = hc.CELL) -> dict[str, Any]:
    """One cache record shaped like the real one, with `value` where a predicted value lives."""
    return {
        "id": "EH38E0000001",
        "chrom": "chr8",
        "genes": [
            {"gene": gene, "by_cell": {cell: value, "K562": value}, "n_tracks": 3},
            {"gene": "OTHER", "by_cell": {"K562": value}},
        ],
    }


def test_a_poison_value_passes_through_the_reader_untouched() -> None:
    """PROOF 1. With a poison value in the cache, the reader still returns the gene names.

    If `genes_with_an_answer` compared the value with zero, took its sign, its absolute value, its
    float, its bool or its string, this test would raise `BlindingError` instead of passing. It is
    therefore not a statement that this lane was careful: it is a statement that the value was not read.
    """
    rec = _record(Poison())
    assert hc.genes_with_an_answer(rec, hc.CELL) == frozenset({"MYC"})
    assert hc.genes_with_an_answer(rec, hc.PRIMARY_CELL) == frozenset({"MYC", "OTHER"})
    assert hc.genes_with_an_answer(rec, "IMR-90") == frozenset()


def test_the_reader_returns_gene_names_and_nothing_else() -> None:
    """A value has no path out of the reader: every member of the result is a gene name."""
    got = hc.genes_with_an_answer(_record(Poison()), hc.CELL)
    assert isinstance(got, frozenset)
    assert all(isinstance(x, str) for x in got)
    assert got <= {"MYC", "OTHER"}


def test_a_null_entry_is_not_an_answer() -> None:
    """Presence is `the entry exists and is not null`, which is the one identity test allowed."""
    assert hc.genes_with_an_answer(_record(None), hc.CELL) == frozenset()


# ---- proof 2: the output cannot distinguish a sign -----------------------------------------------


def test_the_output_is_identical_for_a_negative_a_positive_and_a_zero_value() -> None:
    """PROOF 2. No sign is recoverable from anything this lane computes, even in principle.

    A cached value of -3.5, +3.5 and exactly 0.0 give identical output from the reader and from
    `answer_presence`, so this lane's result cannot carry a direction however it is read.
    """
    outs = [hc.genes_with_an_answer(_record(v), hc.CELL) for v in (-3.5, 3.5, 0.0, 1e-300, -1e-300)]
    assert len({frozenset(o) for o in outs}) == 1
    assert outs[0] == frozenset({"MYC"})


def test_answer_presence_is_identical_for_every_sign_on_disk(tmp_path: Path) -> None:
    """PROOF 2, through the file reader: three caches differing only in sign are indistinguishable."""
    seen = []
    for i, v in enumerate((-2.0, 2.0, 0.0)):
        root = tmp_path / f"cache{i}"
        (root / "chr8").mkdir(parents=True)
        (root / "chr8" / "EH38E0000001.json").write_text(json.dumps(_record(v)))
        got = hc.answer_presence(root, "chr8", "EH38E0000001", "MYC", hc.CELL)
        assert isinstance(got["has_answer"], bool)
        seen.append({k: v2 for k, v2 in got.items() if k != "opened"})
    assert seen[0] == seen[1] == seen[2]
    assert seen[0]["has_answer"] is True
    assert seen[0]["why"] == hc.PRESENT


def test_answer_presence_reports_its_two_kinds_of_absence(tmp_path: Path) -> None:
    """A missing file and a cached window without the pair's own gene are different facts."""
    root = tmp_path / "cache"
    (root / "chr8").mkdir(parents=True)
    (root / "chr8" / "EH38E0000001.json").write_text(json.dumps(_record(-1.0)))
    gone = hc.answer_presence(root, "chr8", "EH38E0000002", "MYC", hc.CELL)
    assert gone == {"element": "EH38E0000002", "opened": None, "has_answer": False, "why": hc.NO_FILE}
    other = hc.answer_presence(root, "chr8", "EH38E0000001", "NOTHERE", hc.CELL)
    assert other["has_answer"] is False
    assert other["why"] == hc.NO_ENTRY
    assert other["opened"] == (root / "chr8" / "EH38E0000001.json").as_posix()


def test_every_per_element_open_is_reported(tmp_path: Path) -> None:
    """The path opened is on the record, so the result can name every per-element open."""
    root = tmp_path / "cache"
    (root / "chr8").mkdir(parents=True)
    p = root / "chr8" / "EH38E0000001.json"
    p.write_text(json.dumps(_record(1.0)))
    got = hc.answer_presence(root, "chr8", "EH38E0000001", "MYC", hc.CELL)
    assert got["opened"] == p.as_posix()
    assert hc.element_file(root, "chr8", "EH38E0000001") == p


# ---- proof 3: the syntax tree carries no value read ----------------------------------------------

BANNED_CALLS = frozenset(
    {"float", "int", "abs", "round", "str", "repr", "sorted", "min", "max", "sum", "bool", "format"}
)
BANNED_COMPARE = (ast.Eq, ast.NotEq, ast.Lt, ast.LtE, ast.Gt, ast.GtE)


def test_the_reader_contains_no_value_read_in_its_own_syntax_tree() -> None:
    """PROOF 3. An AST check, not a grep: a later edit cannot reintroduce a value read silently.

    `genes_with_an_answer` is asserted to contain no ordering or equality comparison, no arithmetic, and
    no call to a numeric or string conversion. Only `is`/`is not`, `in`/`not in` and `isinstance` remain,
    and none of them can recover a sign.
    """
    tree = ast.parse(inspect.getsource(hc.genes_with_an_answer))
    for node in ast.walk(tree):
        if isinstance(node, ast.Compare):
            bad = [op for op in node.ops if isinstance(op, BANNED_COMPARE)]
            assert not bad, f"an equality or ordering comparison appears in the reader: {ast.dump(node)}"
        assert not isinstance(node, ast.BinOp), "arithmetic appears in the reader"
        assert not isinstance(node, ast.UnaryOp) or isinstance(node.op, ast.Not), (
            "a unary arithmetic operator appears in the reader"
        )
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            assert node.func.id not in BANNED_CALLS, f"the reader calls {node.func.id}"


def test_the_module_never_names_a_value_or_a_sign_field() -> None:
    """The module's source names no direction-carrying field of the cache record it reads.

    The record's gene entries also carry `max_drop_log2fc`, `max_rise_log2fc` and `mean_log2fc`. None of
    them is named anywhere in this lane's module, so no edit can read one by a path already written.
    """
    src = Path(hc.__file__).read_text()
    for field in (
        "max_drop_log2fc",
        "max_rise_log2fc",
        "mean_log2fc",
        "max_drop_tissue",
        "max_rise_tissue",
    ):
        assert field not in src, f"{field} is named in this lane's module"


# ---- the floors and the locus convention are imported, not chosen -------------------------------


def test_the_floors_are_imported_and_not_chosen_here() -> None:
    assert hc.POSITIVE_FLOOR is fresh.POSITIVE_FLOOR
    assert hc.LOCUS_FLOOR is fresh.LOCUS_FLOOR
    assert fresh.LOCUS_FLOOR is cell2.POOLED_LOCUS_FLOOR
    assert (hc.POSITIVE_FLOOR, hc.LOCUS_FLOOR) == (30, 20)
    assert hc.FLOORS["neither_chosen_here"] is True
    assert "fresh.POSITIVE_FLOOR" in hc.FLOORS["links_imported_from"]
    assert "cell2.POOLED_LOCUS_FLOOR" in hc.FLOORS["independent_loci_imported_from"]


def test_the_locus_convention_is_cell2s_own_wording() -> None:
    assert hc.LOCUS_RULE is cell2.INDEPENDENT_LOCUS_RULE
    assert hc.LOCUS_SPAN is cell2.INDEPENDENT_LOCUS_SPAN
    assert "not established biological independence" in hc.LOCUS_RULE
    assert "not established biological independence" in hc.NOT_BIOLOGICAL_INDEPENDENCE
    assert "not established biological independence" in hc.POOLING


def test_the_inherited_wording_is_carried_verbatim() -> None:
    assert hc.NAMED_NOT_READ is dl.CELL_AVAILABILITY
    assert "elements_hct116" in hc.NAMED_NOT_READ
    assert hc.INHERITED_GATE_NO_GO is dl.GATE_NO_GO
    assert hc.FORBIDDEN_OF_A_SHORT_COUNT is dl.FORBIDDEN_OF_A_SHORT_OR_UNDETECTED_OUTCOME


# ---- the pooled gate -----------------------------------------------------------------------------


def _link(arm: str, cell: str, chrom: str, start: int, gene: str) -> dict[str, Any]:
    return {
        "arm": arm,
        "cell": cell,
        "chrom": chrom,
        "start": start,
        "end": start + 500,
        "gene": gene,
        "element": f"E{chrom}{start}",
    }


def test_the_pooled_gate_is_lane_directions_own_gate() -> None:
    """The gate is imported, so its floors and its wording are not re-implemented here."""
    rows = [_link("increases", "K562", "chr1", i * 5_000_000, f"G{i}") for i in range(40)]
    g = hc.pooled_gate(rows, [], "increases")
    assert g["imported_gate"] == dl.gate(rows, "increases")
    assert g["links_with_an_answer"] == 40
    assert g["independent_loci"] == 40
    assert g["meets_both_floors"] is True
    assert g["short_by"] == []


def test_a_short_pooled_arm_carries_its_margin() -> None:
    rows = [_link("increases", "K562", "chr1", i * 5_000_000, f"G{i}") for i in range(21)]
    mine = [_link("increases", hc.CELL, "chr8", 1_000_000 + i * 5_000_000, f"H{i}") for i in range(3)]
    g = hc.pooled_gate(rows, mine, "increases")
    assert g["links_with_an_answer"] == 24
    assert g["from_the_primary_cell"] == 21
    assert g["from_this_lanes_cell"] == 3
    assert g["independent_loci"] == 24
    assert g["meets_both_floors"] is False
    assert g["short_by"] == ["answered links 24, short of 30 by 6"]


def test_pooling_a_second_cell_is_not_additive_in_loci() -> None:
    """cell2.group does not look at the cell, so one link of a second cell can merge two loci.

    This is the reason the pooled loci are grouped over the pooled links together and never added up
    per cell, and the reason a pooled locus count can be lower than the sum of the two.
    """
    k = [
        _link("increases", "K562", "chr1", 1_000, "A"),
        _link("increases", "K562", "chr1", 10_000_000, "B"),
    ]
    assert len(set(dl.loci_of(k))) == 2
    bridge = [_link("increases", hc.CELL, "chr1", 1_000, "B")]
    g = hc.pooled_gate(k, bridge, "increases")
    assert g["loci_in_the_primary_cell_alone"] == 2
    assert g["loci_in_this_lanes_cell_alone"] == 1
    assert g["independent_loci"] == 1, "a bridging link of a second cell merges two loci into one"
    assert g["independent_loci"] < g["loci_in_the_primary_cell_alone"] + g["loci_in_this_lanes_cell_alone"]


def test_the_verdict_has_two_readings_and_no_third() -> None:
    short = hc.pooled_gate([_link("increases", "K562", "chr1", 0, "A")], [], "increases")
    wide = hc.pooled_gate(
        [_link("decreases", "K562", "chr1", i * 5_000_000, f"G{i}") for i in range(40)], [], "decreases"
    )
    no_go = hc.verdict({"decreases": wide, "increases": short})
    assert no_go["both_pooled_arms_meet_both_floors"] is False
    assert no_go["reading"] is hc.GATE_NO_GO
    assert no_go["short_by"]["increases"]
    assert no_go["short_by"]["decreases"] == []
    both = hc.verdict({"decreases": wide, "increases": wide})
    assert both["both_pooled_arms_meet_both_floors"] is True
    assert both["reading"] is hc.GATE_PASS
    assert hc.READINGS == (hc.GATE_PASS, hc.GATE_NO_GO)
    assert hc.THERE_IS_NO_THIRD is True


def test_a_short_count_is_never_softened() -> None:
    """The forbidden words appear only inside the prohibitions that name them."""
    for word in hc.FORBIDDEN_OF_A_SHORT_COUNT:
        assert word in hc.GATE_NO_GO, f"{word} must be forbidden by name in the no-go wording"
    assert "NEVER as close" in hc.GATE_NO_GO


def test_the_increase_wording_is_never_a_mechanism() -> None:
    """Review item R2: an increase on knockdown is an increase on knockdown, and nothing more."""
    src = Path(hc.__file__).read_text()
    assert "increase on knockdown" in src
    for claim in ("is a silencer", "is a repressor", "are silencers", "are repressors"):
        assert claim not in src
    assert "never a silencer" in src and "never a repressor" in src


def test_the_prior_exposure_of_this_cell_is_registered() -> None:
    """HCT116's decreases were scored in a39073d, so it is not a fresh cell for this purpose."""
    assert "a39073d" in hc.PRIOR_EXPOSURE
    assert "NOT a fresh cell" in hc.PRIOR_EXPOSURE
    assert (
        hc.verdict({"increases": hc.pooled_gate([], [], "increases")})["prior_exposure"] is hc.PRIOR_EXPOSURE
    )
    assert "PRIOR_EXPOSURE" in hc.GATE_PASS


def test_presence_is_declared_an_upper_bound_before_any_count() -> None:
    assert "UPPER BOUND" in hc.PRESENCE_IS_AN_UPPER_BOUND
    assert "exactly zero" in hc.PRESENCE_IS_AN_UPPER_BOUND
    g = hc.pooled_gate([], [], "increases")
    assert g["counts_are_presence_not_answered"] is hc.PRESENCE_IS_AN_UPPER_BOUND


def test_the_registration_carries_every_binding_part() -> None:
    r = hc.registration()
    assert r["lane"] == "lane-hct"
    for key in (
        "order",
        "blinding",
        "presence_call",
        "presence_is_an_upper_bound",
        "answer_presence_gap",
        "named_not_read_by_lane_direction",
        "inherited_gate_no_go",
        "floors",
        "pooling",
        "independent_locus_rule",
        "not_biological_independence",
        "prior_exposure",
        "cannot_establish",
        "no_requests",
    ):
        assert r[key], key
    assert r["readings"] == [hc.GATE_PASS, hc.GATE_NO_GO]
    assert "0 model requests" in r["no_requests"]


def test_the_lane_cannot_establish_a_mechanism_or_a_usable_cell() -> None:
    assert "not a measurement" in hc.CANNOT_ESTABLISH
    assert "repression mechanism" in hc.CANNOT_ESTABLISH
    assert "upper bound" in hc.CANNOT_ESTABLISH
