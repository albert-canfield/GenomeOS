# SPDX-License-Identifier: AGPL-3.0-or-later
"""The counting lane's own logic: attachment, the locus-status order, and the imported conventions."""

from __future__ import annotations

import importlib.util
from pathlib import Path

from genomeos.attribution import cell2
from genomeos.attribution import measured as ms

ROOT = Path(__file__).resolve().parents[1]


def _module():
    spec = importlib.util.spec_from_file_location("rmc", ROOT / "scripts" / "response_map_coverage.py")
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


rmc = _module()


def _pair(start: int, end: int, chrom: str = "chr1") -> ms.CrispriPair:
    return ms.CrispriPair(
        chrom=chrom,
        start=start,
        end=end,
        gene="G",
        cell="K562",
        dataset="D",
        reference="R",
        regulated=True,
        significant=True,
        effect_size=-0.3,
        p_adjusted=0.01,
    )


def _elements(spans: list[tuple[int, int]]) -> list[dict]:
    return [{"id": f"E{i}", "start": s, "end": e} for i, (s, e) in enumerate(spans)]


def test_attachment_uses_the_measured_layer_rule_and_not_proximity() -> None:
    elements = _elements([(1_000, 1_500), (900_000, 900_500)])
    starts = [e["start"] for e in elements]
    # the same piece of DNA: reciprocal overlap 1.0
    assert rmc.attaches_to(_pair(1_000, 1_500), starts, elements)["id"] == "E0"
    # a single shared base is never an attachment, however close
    assert rmc.attaches_to(_pair(1_499, 2_600), starts, elements) is None
    # an element 898 kb away is within the locus span but is not a measurement of this interval
    assert rmc.attaches_to(_pair(2_000, 2_500), starts, elements) is None


def test_attachment_window_is_wider_than_the_longest_measured_interval() -> None:
    """A pair must be found however long it is, so the scan window is `measured.REACH` on each side."""
    elements = _elements([(0, ms.REACH)])
    starts = [e["start"] for e in elements]
    assert rmc.attaches_to(_pair(10, ms.REACH), starts, elements)["id"] == "E0"


def test_the_locus_status_order_prefers_the_best_status_a_locus_reaches() -> None:
    order = rmc.STATUS_ORDER
    assert order[rmc.STATUS_ADDABLE] == 0
    assert order[rmc.STATUS_NO_READER] < order[rmc.STATUS_NO_RULE] < order[rmc.STATUS_NO_LINK]
    assert set(rmc.STATUSES) == set(order)


def test_the_locus_rule_is_imported_word_for_word_and_not_restated() -> None:
    """The convention is `cell2`'s; this lane carries it, so the two can never read differently."""
    text = (ROOT / "scripts" / "response_map_coverage.py").read_text()
    assert "cell2.INDEPENDENT_LOCUS_RULE" in text
    assert "1 Mb of one another" not in text  # the rule's own wording is never retyped here
    assert "not established biological independence" in cell2.INDEPENDENT_LOCUS_RULE
    assert "not established biological independence" in text


def test_only_crispri_counts_as_a_perturbation_of_the_native_locus() -> None:
    assert rmc.PERTURBATION_ASSAY == "crispri"
    assert set(rmc.REPORTER_ASSAYS) == set(ms.ASSAYS) - {"crispri"}
    assert rmc.REPORTER_ASSAYS, "the measured layer holds assays other than the perturbation one"
