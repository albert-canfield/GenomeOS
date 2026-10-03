# SPDX-License-Identifier: AGPL-3.0-or-later
"""gwas.py asks whether ANY gene at the bar is the catalogue's mapped gene, not only the head.

Wave 2 of section 5 item 10, under the registration at `genomeos/attribution/onetarget2.py`
(24adf33). gwas.py is the one module in the census whose headline figure can only move ONE WAY, and
its registered falsifier says so: "the agreement share not rising. This is the one module in the
census whose headline figure can only move one way under the window -- an element agrees if ANY gene
at the bar is the mapped gene, which is a superset of agreeing on the head -- so a share that falls
or an element that loses its agreement is a defect in the move, not a finding, and stops the wave."

So `window_agreement` RAISES on an element that agreed on its head and not on its window, and that
refusal is planted below. A share that rose is a result; a share that fell is never published.
"""

from __future__ import annotations

import pytest

from genomeos.attribution import gwas
from genomeos.attribution import onetarget2 as ot


def _g(gene: str, drop: float, rise: float) -> dict:
    return {
        "gene": gene,
        "max_drop_log2fc": drop,
        "max_rise_log2fc": rise,
        "max_drop_tissue": "K562",
        "max_rise_tissue": "occipital lobe",
    }


class _Reader:
    def __init__(self, records):
        self.records = records

    def element(self, chrom, element_id):
        return self.records.get(element_id)

    def at_bar(self, chrom, element_id, min_effect=0.1, genes=None):
        from genomeos.attribution.targets import genes_at_bar

        got = genes_at_bar(self.element(chrom, element_id), min_effect)
        if got is None or genes is None:
            return got
        return [(g, v) for g, v in got if g in genes]


def _e(eid, head, mapped, chrom="chr21"):
    """One element and the catalogue hit for it. `summarise` OVERWRITES e['gwas'] from `hits`, which
    the first form of this fixture did not do, so the hits are built from the same rows below and
    the mapped gene is carried on the element only to build them."""
    return {
        "id": eid,
        "key": f"{chrom}:{eid}",
        "chrom": chrom,
        "start": 1000,
        "end": 1200,
        "predicted": {"gene": head, "strength": "strong"} if head else None,
        "mapped_for_the_fixture": mapped,
    }


def _hits(rows):
    return {r["key"]: [{"mapped_gene": r["mapped_for_the_fixture"], "trait": "height"}] for r in rows}


#: Varied on the axes the guards turn on: A1 agrees on its head; A2 agrees only in the window; A3
#: agrees nowhere; A4 is not cached; A5's mapped field uses the catalogue's " - " intergenic spelling.
RECORDS = {
    "A1": {"genes": [_g("APP", -0.8, 0.05)]},
    "A2": {"genes": [_g("APP", -0.8, 0.05), _g("SOD1", -0.1, 0.4)]},
    "A3": {"genes": [_g("APP", -0.8, 0.05)]},
    "A4": None,
    "A5": {"genes": [_g("APP", -0.8, 0.05), _g("CYYR1", -0.3, 0.02)]},
}
READER = _Reader(RECORDS)


def test_summarise_without_a_reader_has_neither_new_key():
    rows = [_e("A1", "APP", "APP")]
    out = gwas.summarise(rows, _hits(rows))
    assert "catalog_gene_is_at_the_bar" not in out and "window" not in out
    assert gwas.window_agreement([], None) == {}


def test_summarise_with_a_reader_moves_no_committed_field():
    rows = [_e("A1", "APP", "APP"), _e("A2", "APP", "SOD1")]
    off = gwas.summarise([dict(r) for r in rows], _hits(rows))
    on = gwas.summarise([dict(r) for r in rows], _hits(rows), READER)
    for k in off:
        assert off[k] == on[k], k
    assert set(on) - set(off) == {"catalog_gene_is_at_the_bar", "window"}


def test_the_share_rises_and_the_gained_element_is_named():
    rows = [_e("A1", "APP", "APP"), _e("A2", "APP", "SOD1"), _e("A3", "APP", "TP53")]
    out = gwas.summarise(rows, _hits(rows), READER)
    assert out["catalog_gene_is_the_predicted_target"] == round(1 / 3, 3)
    assert out["catalog_gene_is_at_the_bar"] == round(2 / 3, 3)
    w = out["window"]
    assert w["agree_on_the_head"] == 1 and w["agree_anywhere_in_the_window"] == 2
    assert w["elements_gaining_an_agreement"] == 1
    assert w["gained"][0]["head"] == "APP" and w["gained"][0]["mapped"] == ["SOD1"]
    assert w["gained"][0]["rank_in_the_window"] == 2
    assert "FIRES" in w["verdict_against_its_own_falsifier"]


def test_an_element_not_in_the_cache_is_counted_and_keeps_its_head_answer():
    rows = [_e("A4", "APP", "APP")]
    out = gwas.summarise(rows, _hits(rows), READER)
    assert out["window"]["not_cached"] == 1
    #: The head still agrees and the window cannot take that away: an uncached element keeps its
    #: compact answer, which is WindowReading's contract and what keeps the share from falling on
    #: elements that were never scored.
    assert out["catalog_gene_is_the_predicted_target"] == 1.0
    assert out["catalog_gene_is_at_the_bar"] == 1.0
    assert out["window"]["elements_gaining_an_agreement"] == 0


def test_the_catalogues_intergenic_spelling_is_read_by_the_same_rule():
    """A mapped_gene of 'APP - CYYR1' names two genes, which is the rule `summarise` has always
    used and which `_mapped_genes` reproduces for the window."""
    assert gwas._mapped_genes({"mapped_gene": "APP - CYYR1"}) == ["APP", "CYYR1"]
    #: A5's window is [APP, CYYR1] and its head is APP, so the head agrees with NEITHER mapped gene
    #: while the window agrees with CYYR1 -- the gain comes from the second name in the pair.
    rows = [_e("A5", "APP", "CYYR1 - TP53")]
    out = gwas.summarise(rows, _hits(rows), READER)
    assert out["catalog_gene_is_the_predicted_target"] == 0.0
    assert out["catalog_gene_is_at_the_bar"] == 1.0
    assert out["window"]["gained"][0]["mapped"] == ["CYYR1", "TP53"]


def test_an_element_that_loses_its_agreement_refuses(monkeypatch):
    """PLANT: the mapped gene IS the head and is NOT at the bar, which cannot happen because the
    head is always at its own bar. The match is on 'The share is not published', a phrase only this
    guard's message uses."""
    rows = [_e("A1", "TP53", "TP53")]  # TP53 is the head and is absent from A1's window
    with pytest.raises(gwas.WindowLostAnAgreementError, match="The share is not published"):
        gwas.summarise(rows, _hits(rows), READER)


def test_an_any_gene_head_disagreement_stops_the_wave():
    """PLANT: the head is not the window's strongest gene. The refusal is check_head_invariant's."""
    rows = [_e("A2", "SOD1", "SOD1")]  # SOD1 is at the bar but APP is the window's head
    with pytest.raises(ot.InvariantRefutedError, match="THE WAVE STOPS"):
        gwas.summarise(rows, _hits(rows), READER)


def test_the_one_way_rule_and_the_limit_are_carried_in_the_arm():
    out = gwas.summarise([_e("A1", "APP", "APP")], _hits([_e("A1", "APP", "APP")]), READER)
    assert out["window"]["can_only_rise"] == gwas.WINDOW_CAN_ONLY_RISE
    assert out["window"]["limit"] == gwas.WINDOW_READS_THE_ANY_GENE_HEAD
    assert "does NOT fire" in out["window"]["verdict_against_its_own_falsifier"]
