# SPDX-License-Identifier: AGPL-3.0-or-later
"""regdiff.py reads the sweep's whole window beside the one gene cached_prediction keeps.

Wave 2 of section 5 item 10, under the registration at `genomeos/attribution/onetarget2.py`
(24adf33) and its file (aed8ae9). regdiff.py's registered falsifier is one of the three the
registration expects NOT to fire: `with_predicted_target_only_a` and `_only_b` are conditioned on
`predicted` being truthy, which is the yes/no and not the gene, so the window must reproduce them
exactly. A move there refutes the invariant the wave stands on and stops the wave; it is never a
finding. The GENE changing is the expected gain and is counted separately.

Three things are asserted here and each one is a different claim:

  the control -- with no reader every committed field is what it was, key for key, and no `window`
  key exists at all, so a caller of the old shape cannot see this change;
  the gain -- a window that names genes the compact head dropped is reported as `extra`, and does
  not touch either count;
  the refusals -- a with_predicted_target count that moves, and an any-gene head disagreement,
  are both planted and both must raise.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from genomeos.attribution import onetarget2 as ot
from genomeos.attribution.eqtl import Intervals
from genomeos.attribution.targets import prediction_window_reading
from genomeos.genome import regdiff


def _element(eid: str, cls: str, start: int, end: int, gene: str | None) -> SimpleNamespace:
    return SimpleNamespace(
        id=eid,
        cls=cls,
        locus=SimpleNamespace(start=start, end=end),
        targets=[{"gene": gene, "basis": "nearest TSS in domain", "distance": 500}] if gene else [],
    )


def _g(gene: str, drop: float, rise: float) -> dict:
    """One gene of a cached record, in the sweep's own field names.

    BOTH sizes are required and numeric because that is what the sweep writes -- every gene of every
    record under `data/knowledge/alphagenome/elements` carries `max_drop_log2fc` and
    `max_rise_log2fc`, and `predict_target` indexes both without a default. A fixture that left the
    rise at None raised a TypeError here and taught that lesson in this file's first run.
    """
    return {
        "gene": gene,
        "max_drop_log2fc": drop,
        "max_rise_log2fc": rise,
        "max_drop_tissue": "K562",
        "max_rise_tissue": "occipital lobe",
    }


#: The fixture VARIES ON THE AXIS EACH GUARD TURNS ON, which is today's lesson from two guards that
#: passed because every fixture happened to carry the same curated field. The axes varied here are
#: the element's CLASS (asked or not), whether the cache holds it, whether any gene reaches the bar,
#: and whether a gene reaches it by its FALL or by its RISE -- the two halves of `predict_target`'s
#: size rule, which a window made only of falls would never exercise.
#:
#: E1 enhancer, three genes at the bar, and CYYR1 reaches it by its rise (a repressed gene) while
#:    the head APP reaches it by its fall;
#: E2 promoter, a window naming exactly its head and nothing else;
#: E3 insulator, never asked, so no window is read for it;
#: E4 enhancer the cache does not hold at all;
#: E5 enhancer the cache holds with every gene BELOW the bar, so the head is None and the window is
#:    empty -- which is not the same as E4 and must not be counted with it.
RECORDS: dict[str, dict | None] = {
    "E1": {"genes": [_g("APP", -0.8, 0.05), _g("CYYR1", -0.1, 0.4), _g("ADAMTS1", -0.15, 0.02)]},
    "E2": {"genes": [_g("SOD1", -0.5, 0.01)]},
    "E4": None,
    "E5": {"genes": [_g("RUNX1", -0.05, 0.02)]},
}


def _patch(monkeypatch, a: dict, b: dict, records=None) -> None:
    by_id = {
        "E1": _element("E1", "enhancer", 1000, 1300, "APP"),
        "E2": _element("E2", "promoter", 5000, 5400, "SOD1"),
        "E3": _element("E3", "insulator", 9000, 9200, None),
        "E4": _element("E4", "enhancer", 11000, 11400, "RUNX1"),
        "E5": _element("E5", "enhancer", 13000, 13400, "RUNX1"),
    }
    iv = Intervals()
    for e in by_id.values():
        iv.add("chr21", e.locus.start, e.locus.end, e.id)
    monkeypatch.setattr(regdiff, "chromosome_elements", lambda chrom, reference=None: (by_id, iv.freeze()))
    monkeypatch.setattr(regdiff, "_reference_base_reader", lambda chrom, reference=None: (None, None))
    monkeypatch.setattr(regdiff, "vcf_path", lambda name, chrom, root: f"{name}.vcf")
    calls = {"A": a, "B": b}
    monkeypatch.setattr(regdiff, "_genotypes", lambda path, base_at=None: calls[path.split(".")[0]])

    import genomeos.predict.enhancer_target as et

    recs = RECORDS if records is None else records
    monkeypatch.setattr(et, "load_cached", lambda chrom, eid, cache=None: recs.get(eid))

    def pred(chrom, eid, cache=None):
        from genomeos.predict.enhancer_target import predict_target

        r = recs.get(eid)
        return predict_target(r["genes"]) if r else None

    monkeypatch.setattr(et, "cached_prediction", pred)


A_SIDE = {
    (1100, "C", "T"): "het",  # E1, three genes at the bar
    (5100, "G", "A"): "hom",  # E2, only its head
    (11100, "A", "G"): "het",  # E4, not cached
    (13100, "T", "A"): "het",  # E5, cached and nothing at the bar
}
B_SIDE = {(9100, "T", "C"): "het"}  # E3, an insulator: never asked


# --------------------------------------------------------------------------------------------
# the control: without a reader, nothing moved
# --------------------------------------------------------------------------------------------


def test_without_a_reader_every_field_is_what_it_was(tmp_path, monkeypatch):
    _patch(monkeypatch, A_SIDE, B_SIDE)
    old = regdiff.regulatory_diff("A", "B", "chr21", tmp_path)
    new = regdiff.regulatory_diff("A", "B", "chr21", tmp_path, window=False)
    assert old == new
    assert "window" not in old
    for g in old["genes"]:
        for v in g["only_a"] + g["only_b"]:
            assert "window" not in v


def test_with_a_reader_the_two_counts_are_identical_to_without_one(tmp_path, monkeypatch):
    _patch(monkeypatch, A_SIDE, B_SIDE)
    off = regdiff.regulatory_diff("A", "B", "chr21", tmp_path)
    on = regdiff.regulatory_diff("A", "B", "chr21", tmp_path, window=True)
    for k in (
        "elements",
        "variants_in_elements_a",
        "variants_in_elements_b",
        "shared",
        "only_a",
        "only_b",
        "by_class",
        "with_predicted_target_only_a",
        "with_predicted_target_only_b",
    ):
        assert off[k] == on[k], k
    assert set(on) - set(off) == {"window"}


# --------------------------------------------------------------------------------------------
# the gain, and the three populations kept apart
# --------------------------------------------------------------------------------------------


def test_the_window_names_the_genes_the_head_dropped(tmp_path, monkeypatch):
    _patch(monkeypatch, A_SIDE, B_SIDE)
    d = regdiff.regulatory_diff("A", "B", "chr21", tmp_path, window=True)
    rows = {v["element"]: v for g in d["genes"] for v in g["only_a"]}
    assert rows["E1"]["predicted"]["gene"] == "APP"  # the head, carried through unchanged
    assert rows["E1"]["window"]["extra"] == ["CYYR1", "ADAMTS1"]
    assert rows["E1"]["window"]["head_agrees"] is True
    assert rows["E2"]["window"]["extra"] == []  # a window naming exactly its head gains nothing
    assert d["window"]["only_a"]["variants_gaining_a_gene"] == 1
    assert d["window"]["only_a"]["genes_gained"] == 2


def test_an_element_the_cache_does_not_hold_is_not_cached_and_is_not_a_zero(tmp_path, monkeypatch):
    _patch(monkeypatch, A_SIDE, B_SIDE)
    d = regdiff.regulatory_diff("A", "B", "chr21", tmp_path, window=True)
    rows = {v["element"]: v for g in d["genes"] for v in g["only_a"]}
    e4 = rows["E4"]["window"]
    assert e4["not_cached"] is True
    assert e4["at_bar"] == [] and e4["anything_at_bar"] is None and e4["head_agrees"] is None
    assert d["window"]["only_a"]["not_cached"] == 1
    #: E5 IS cached and has nothing at the bar, which is a different answer and is kept apart: the
    #: window is empty rather than absent, so `anything_at_bar` is False and not None.
    e5 = rows["E5"]["window"]
    assert e5["not_cached"] is False
    assert e5["at_bar"] == [] and e5["anything_at_bar"] is False and e5["head_agrees"] is True
    assert rows["E5"]["predicted"] is None  # nothing reaches MIN_EFFECT, so the head is None too


def test_a_class_the_module_never_asks_about_gets_no_window_at_all(tmp_path, monkeypatch):
    """E3 is an insulator: `cached_prediction` is not called for it today, so no window is read for
    it either. Reading one would make the two counts compare different populations of variants."""
    _patch(monkeypatch, A_SIDE, B_SIDE)
    d = regdiff.regulatory_diff("A", "B", "chr21", tmp_path, window=True)
    rows = {v["element"]: v for g in d["genes"] for v in g["only_b"]}
    assert "window" not in rows["E3"]
    assert d["window"]["only_b"]["variants_whose_element_was_asked"] == 0


def test_the_arm_contributes_no_key_when_it_was_not_asked_for():
    assert regdiff._window_arm(False, [], []) == {}


def test_the_registered_limit_is_stated_and_says_the_coding_head_is_not_read():
    assert "coding" in regdiff.WINDOW_DOES_NOT_READ_THE_CODING_HEAD
    d = regdiff._window_arm(True, [], [])
    assert d["limit"] == regdiff.WINDOW_DOES_NOT_READ_THE_CODING_HEAD
    assert d["falsifier"] == regdiff.WINDOW_FALSIFIER


# --------------------------------------------------------------------------------------------
# the refusals, planted
# --------------------------------------------------------------------------------------------


def test_a_count_cannot_move_without_the_head_guard_refusing_first(tmp_path, monkeypatch):
    """PLANT, and a FINDING about the registered falsifier: it is SUBSUMED by the invariant guard.

    Planted here is exactly what "the count moved" means -- a window with a gene at the bar while
    `cached_prediction` says None. Reached through `regulatory_diff` it does NOT produce
    `the_count_moved`: `check_head_invariant` refuses first, because head None against a non-empty
    bar is an any-gene head disagreement. Enumerating the five shapes an asked, cached row can take
    shows this is not an accident of ordering -- every row that could move a count is a row whose
    head disagrees, so no input the sweep or a bug can produce reaches the count branch through the
    public entry point. The count is still computed and still reported, because the registration
    names it as this module's falsifier and the figure belongs in the record; what is recorded here
    is that it is strictly weaker than the refusal above it. The match is on 'THE WAVE STOPS',
    a phrase only `check_head_invariant`'s message uses.
    """
    _patch(monkeypatch, {(1100, "C", "T"): "het"}, {})
    import genomeos.predict.enhancer_target as et

    monkeypatch.setattr(et, "cached_prediction", lambda chrom, eid, cache=None: None)
    with pytest.raises(ot.InvariantRefutedError, match="THE WAVE STOPS"):
        regdiff.regulatory_diff("A", "B", "chr21", tmp_path, window=True)


def test_the_count_moved_branch_itself_reports_a_refutation_when_it_is_reached(tmp_path):
    """The branch above is unreachable through `regulatory_diff`, so it is driven directly.

    This row is self-inconsistent on purpose -- `head_agrees` True beside `anything_at_bar` False
    and a `predicted` gene -- which is a shape only a defect inside `_window_arm` itself could make.
    It is here so the branch is exercised rather than left standing untested, and so a future change
    that made the two counts genuinely divergeable finds a guard already in place.
    """
    row = {
        "element": "E9",
        "predicted": {"gene": "APP"},
        "window": {
            "at_bar": [],
            "extra": [],
            "not_cached": False,
            "head_agrees": True,
            "anything_at_bar": False,
        },
    }
    d = regdiff._window_arm(True, [row], [])
    assert d["only_a"]["with_predicted_target_head"] == 1
    assert d["only_a"]["with_predicted_target_window"] == 0
    assert d["only_a"]["the_count_moved"] is True
    assert "REFUTES THE INVARIANT" in d["verdict_against_its_own_falsifier"]


def test_an_any_gene_head_disagreement_stops_the_wave(tmp_path, monkeypatch):
    """PLANT: the head is not the window's strongest gene, so the invariant itself is refuted.

    The refusal is `onetarget2.check_head_invariant`'s, not this module's, and that is deliberate:
    one refusal for one invariant, never a substitute for a module's own falsifier. The match is on
    'THE WAVE STOPS', a phrase only that guard's message uses.
    """
    _patch(monkeypatch, {(1100, "C", "T"): "het"}, {})
    import genomeos.predict.enhancer_target as et

    monkeypatch.setattr(
        et,
        "cached_prediction",
        lambda chrom, eid, cache=None: {"gene": "ADAMTS1", "log2_fold_change": -0.15, "tissue": "K562"},
    )
    with pytest.raises(ot.InvariantRefutedError, match="THE WAVE STOPS"):
        regdiff.regulatory_diff("A", "B", "chr21", tmp_path, window=True)


def test_the_adaptor_does_not_recompute_the_head_it_is_given():
    """`prediction_window_reading` carries the caller's head through; it never derives one.

    This is the property that keeps a consumer's committed figures its own: if the adaptor
    recomputed the head, a move in `predict_target` would silently rewrite every head in the tree.
    """
    w = prediction_window_reading(RECORDS["E1"], "CYYR1")
    assert w.head == "CYYR1"
    assert w.genes == ("APP", "CYYR1", "ADAMTS1")
    assert w.head_agrees is False  # it disagrees, and the adaptor says so rather than fixing it


def test_the_adaptor_with_no_record_is_the_control():
    w = prediction_window_reading(None, "APP")
    assert w.head == "APP" and w.not_cached is True and w.extra == () and w.genes == ()


def test_the_adaptor_forms_a_coding_window_only_when_asked():
    plain = prediction_window_reading(RECORDS["E1"], "APP")
    assert plain.coding_at_bar is None and plain.coding_genes == ()
    coded = prediction_window_reading(RECORDS["E1"], "APP", coding_head="APP", coding={"APP", "CYYR1"})
    assert coded.coding_genes == ("APP", "CYYR1")
    assert coded.coding_head_agrees is True
