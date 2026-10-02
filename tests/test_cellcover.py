# SPDX-License-Identifier: AGPL-3.0-or-later
"""Targeted tests for the cell-coverage costing (`genomeos.attribution.cellcover`).

The binding ones are: the label rule is the adapter's own rule and cannot drift from it; the
request count is over elements and never over rows or cells; one track means unpurchasable at any
price; and the readings carried from direction_v2 are the same objects, so none can be
strengthened by being retyped here.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import pytest

from genomeos.attribution import cellcover as cc
from genomeos.attribution import direction_v2 as dv
from genomeos.predict import alphagenome_adapter as ad
from genomeos.predict.enhancer_target import CELLS

# ---- the label rule is the adapter's own ------------------------------------------------------


class _Axis:
    """The minimum of a response's var axis: the two columns `tissue_names` reads."""

    def __init__(self, names: list[str], gtex: list[object]) -> None:
        self._d = {"biosample_name": names, "gtex_tissue": gtex}
        self.index = names

    def get(self, key, default=None):  # pragma: no cover - exercised through tissue_names
        return self._d.get(key, default)

    def __getitem__(self, key):
        return self._d[key]

    def __contains__(self, key):
        return key in self._d

    def keys(self):
        return self._d.keys()


class _Response:
    def __init__(self, names: list[str], gtex: list[object]) -> None:
        self.var = _Axis(names, gtex)


@pytest.mark.parametrize(
    "name,gtex",
    [
        ("K562", ""),  # an ENCODE cell line: no GTEx tissue, so the biosample name is the label
        ("Whole_Blood", "Whole_Blood"),  # a GTEx tissue: the GTEx name is the label
        ("liver", "nan"),  # the string "nan", which is not a tissue
        ("CD14-positive monocyte", None),
        ("some biosample", "Artery_Tibial"),  # both present: GTEx wins
    ],
)
def test_track_label_is_the_adapters_own_rule(name, gtex):
    """`track_label` and `alphagenome_adapter.tissue_names` must agree on every branch.

    A difference between them would make the per-cell counts describe a different labelling than
    the one a cell's value actually arrives under, which is the whole basis of the match.
    """
    mine = cc.track_label({"biosample_name": name, "gtex_tissue": gtex})
    theirs = ad.tissue_names(_Response([name], [gtex]))
    assert [mine] == theirs


def test_track_label_never_invents_a_name():
    assert cc.track_label({"biosample_name": "", "gtex_tissue": ""}) == ""


# ---- the request unit -------------------------------------------------------------------------


def _row(element: str, cell: str, gene: str = "G") -> dict[str, object]:
    return {"element": element, "cell": cell, "gene": gene}


def test_requests_are_counted_over_elements_and_not_over_rows():
    """REQUEST_UNIT: one request buys one element, whatever cells or genes the rows name."""
    rows = [_row("E1", "liver"), _row("E1", "placenta"), _row("E1", "liver", gene="H")]
    assert cc.price(rows).requests == 1
    assert cc.price(rows).rows == 3


def test_requests_are_not_multiplied_by_cells():
    rows = [_row("E1", "liver"), _row("E2", "liver"), _row("E2", "placenta")]
    p = cc.price(rows)
    assert p.requests == 2
    assert len(p.cells) == 2
    assert p.requests != p.rows


# ---- purchasable and unpurchasable ------------------------------------------------------------


def _axis(labels: list[str]) -> cc.Axis:
    """A stand-in axis: one distinct track name per entry, labelled by biosample name."""
    tracks = {f"t{i}": {"biosample_name": lab, "gtex_tissue": ""} for i, lab in enumerate(labels)}
    return cc.Axis(tracks=tracks, rows_read=len(labels), source="stub")


def test_one_track_is_unpurchasable_at_any_price():
    """A cell with one column can never yield two values, so no request count resolves it."""
    cover = cc.coverage(Counter({"lonely": 3}), _axis(["lonely"]))
    assert cover[0].rna_tracks == 1
    assert cover[0].track_exists if hasattr(cover[0], "track_exists") else True
    assert cover[0].purchasable is False
    assert cover[0].established_how == "client_track_metadata"


def test_two_tracks_are_purchasable():
    cover = cc.coverage(Counter({"paired": 1}), _axis(["paired", "paired"]))
    assert cover[0].rna_tracks == 2
    assert cover[0].purchasable is True


def test_a_cell_no_track_carries_is_recorded_as_absent_not_as_one():
    """`absent` and `one track` are different findings and must not collapse into each other."""
    cover = cc.coverage(Counter({"nowhere": 2}), _axis(["elsewhere"]))
    assert cover[0].rna_tracks == 0
    assert cover[0].established_how == "absent"
    assert cover[0].purchasable is False
    assert cover[0].to_dict()["track_exists"] is False


def test_split_sends_one_track_rows_to_unpurchasable():
    axis = _axis(["paired", "paired", "lonely"])
    rows = [_row("E1", "paired"), _row("E2", "lonely"), _row("E3", "lonely")]
    cover = cc.coverage(Counter({"paired": 1, "lonely": 2}), axis)
    parts = cc.split(rows, cover)
    assert cc.price(parts["purchasable"]).requests == 1
    assert cc.price(parts["unpurchasable_at_any_price"]).rows == 2


def test_established_how_distinguishes_a_confirmation_from_a_bare_name():
    """HOW_ESTABLISHED: the two tiers must be separate tokens, so a row cannot read as both."""
    assert "client_track_metadata" in cc.HOW_ESTABLISHED
    assert "name_only" in cc.HOW_ESTABLISHED
    assert "absent" in cc.HOW_ESTABLISHED
    tokens = {c.established_how for c in cc.coverage(Counter({"a": 1, "b": 1}), _axis(["a"]))}
    assert tokens == {"client_track_metadata", "absent"}


# ---- the readings are carried, not retyped ----------------------------------------------------


def test_readings_are_the_same_objects_as_direction_v2s():
    """A carried reading must be the same object, so it cannot be strengthened in transit."""
    assert cc.READINGS_CARRIED["unresolved_means"] is dv.UNRESOLVED_MEANS
    assert cc.READINGS_CARRIED["flips_are_structurally_impossible"] is dv.FLIPS_ARE_STRUCTURALLY_IMPOSSIBLE
    assert cc.READINGS_CARRIED["assigned_cell"] is dv.ASSIGNED_CELL
    assert cc.READINGS_CARRIED["cache_is_legacy"] is dv.CACHE_IS_LEGACY


def test_the_reasons_priced_are_members_of_the_registered_tuple():
    assert cc.ONE_VALUE_ONLY in dv.UNRESOLVED_REASONS
    assert cc.ONE_TRACK_SEEN_TWICE in dv.UNRESOLVED_REASONS


def test_values_needed_is_v2s_clause_one_and_not_a_new_threshold():
    assert cc.VALUES_NEEDED == 2
    assert "at least two values are retained" in dv.SIGN_RULE_V2


def test_no_recommendation_is_stated_and_nothing_argues_for_the_purchase():
    for word in ("worth", "should buy", "recommend spending", "cheap", "bargain"):
        assert word not in cc.WHAT_IT_DOES_NOT_BUY.lower()
    assert "no recommendation" in cc.NO_RECOMMENDATION.lower()


def test_what_it_does_not_buy_names_all_three_surviving_clauses():
    t = cc.WHAT_IT_DOES_NOT_BUY.lower()
    assert "distinct" in t and "unanimous" in t and "magnitude_floor" in t


def test_the_premise_correction_names_what_the_code_says():
    t = cc.PREMISE_CORRECTED
    assert "ontology_terms" in t and "371" in t and "retention" in t


# ---- the module sends nothing ------------------------------------------------------------------


#: Names that must never be CALLED by this lane's committed code. A prose mention is allowed and is
#: in fact required -- AD_HOC_READ_INCIDENT names `load_cached` to say why it is not reached -- so the
#: test parses the code and looks at call sites and imports, never at the text.
MUST_NOT_CALL = ("load_cached", "score_variant", "create_client", "_live_scorer", "getenv")
LANE_CODE = (
    "genomeos/attribution/cellcover.py",
    "scripts/cellcover.py",
    "scripts/cellcover_register.py",
)


def _called_names(src: str) -> set[str]:
    """Every name this source calls, by attribute or bare, from its AST and not from its text."""
    import ast

    out: set[str] = set()
    for node in ast.walk(ast.parse(src)):
        if isinstance(node, ast.Call):
            f = node.func
            if isinstance(f, ast.Attribute):
                out.add(f.attr)
            elif isinstance(f, ast.Name):
                out.add(f.id)
    return out


def _imported_names(src: str) -> set[str]:
    import ast

    out: set[str] = set()
    for node in ast.walk(ast.parse(src)):
        if isinstance(node, ast.Import):
            out.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            out.add(node.module or "")
            out.update(a.name for a in node.names)
    return out


@pytest.mark.parametrize("path", LANE_CODE)
def test_this_lane_calls_nothing_that_sends_a_request(path):
    """No committed file of this lane may call the client, build a scorer or read a key."""
    called = _called_names(Path(path).read_text())
    assert not (called & set(MUST_NOT_CALL)), f"{path} calls {sorted(called & set(MUST_NOT_CALL))}"


@pytest.mark.parametrize("path", LANE_CODE)
def test_this_lane_imports_no_client(path):
    imported = _imported_names(Path(path).read_text())
    assert "alphagenome" not in {n.split(".")[0] for n in imported}, path
    assert "os" not in imported, path


@pytest.mark.parametrize("path", LANE_CODE)
def test_the_committed_program_opens_no_chromosome_archive(path):
    """AD_HOC_READ_INCIDENT: `load_cached` inflates whole archives and must not be reached.

    Checked as a call site, because the incident is also DESCRIBED in this lane's prose by name,
    and a text search cannot tell the description from the deed.
    """
    src = Path(path).read_text()
    assert "load_cached" not in _called_names(src), path
    assert ".json.gz" not in src, path


def test_the_axis_witness_is_a_loose_file_and_not_an_archive():
    assert cc.AXIS_WITNESS.suffix == ".json"
    assert ".gz" not in str(cc.AXIS_WITNESS)


def test_the_figures_from_the_forbidden_read_are_marked_unreproducible():
    """A figure the committed program cannot recompute must say so on its own record."""
    assert cc.MEASURED_BY_A_READ_NOW_FORBIDDEN["reproducible"] is False
    assert "not re-run" in cc.AD_HOC_READ_INCIDENT or "not re-run" in str(cc.MEASURED_BY_A_READ_NOW_FORBIDDEN)


def test_retained_cells_today_are_the_four_the_writer_keeps():
    assert CELLS == ("K562", "HepG2", "GM12878", "IMR-90")
    assert str(list(CELLS)) in cc.PREMISE_CORRECTED
