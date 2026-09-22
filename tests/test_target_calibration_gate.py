# SPDX-License-Identifier: AGPL-3.0-or-later
"""The calibration's top-target gate, counted: who is admitted, who is excluded, and why.

Synthetic and local: four pairs on one element, a per-element response cache written into tmp_path,
no network and no cached data. The point of each test is that an excluded pair is either a censoring
the sweep can undo or a named silence it cannot, and never a zero.
"""

import gzip
import io
import json

from genomeos.attribution import crispri
from genomeos.attribution import target_calibration as tc
from genomeos.attribution.targets import NOT_IN_WINDOW, SCORED, ElementResponses

HEADER = (
    "chrom\tchromStart\tchromEnd\tname\tEffectSize\tmeasuredGeneSymbol\tValidConnection\tCellType\t"
    "Regulated\tDataset\tdistanceToTSS\tDHS.RPM\tH3K27ac.RPM\tstartTSS\n"
)


def _row(gene, regulated):
    return (
        f"chr1\t1000\t1500\t{gene}|x\t-0.1\t{gene}\tTRUE\tK562\t"
        f"{'TRUE' if regulated else 'FALSE'}\tD\t30000\t1.0\t1.0\t31000\n"
    )


def _world(tmp_path):
    """One element whose compact table names TOP, and a cache that also holds OTHER and RISES."""
    elements = tmp_path / "elements"
    elements.mkdir()
    pred = {"gene": "TOP", "log2_fold_change": -0.4, "tissue": "K562"}
    (elements / "chr1.json").write_text(
        json.dumps(
            [
                {
                    "id": "EH1",
                    "start": 1_000,
                    "end": 1_500,
                    "inferred": {"gene": "TOP", "distance": 1000, "basis": "nearest TSS in domain"},
                    "predicted": pred,
                    "predicted_by_cell": {"K562": -0.4},
                    "predicted_coding": pred,
                    "predicted_coding_by_cell": {"K562": -0.4},
                }
            ]
        )
    )
    cache = tmp_path / "cache"
    cache.mkdir()
    with gzip.open(cache / "chr1.json.gz", "wt") as fh:
        json.dump(
            {
                "EH1": {
                    "id": "EH1",
                    "genes": [
                        {"gene": "TOP", "by_cell": {"K562": -0.4}, "max_drop_log2fc": -0.4},
                        {"gene": "OTHER", "by_cell": {"K562": -0.25}, "max_drop_log2fc": -0.25},
                        {"gene": "RISES", "by_cell": {"K562": 0.3}, "max_drop_log2fc": -0.01},
                    ],
                }
            },
            fh,
        )
    return elements, cache


def _pairs(tmp_path, cache_on):
    elements, cache = _world(tmp_path)
    text = HEADER + _row("TOP", True) + _row("OTHER", True) + _row("RISES", False) + _row("ABSENT", False)
    pairs = crispri.parse(io.StringIO(text))
    table = crispri.DeletionTable(elements)
    crispri.annotate(pairs, table, crispri.ElementCache(cache) if cache_on else None)
    return pairs, table, ElementResponses(cache)


def test_the_gate_admits_one_pair_of_four_and_the_two_sides_have_different_rates(tmp_path):
    pairs, _, _ = _pairs(tmp_path, cache_on=False)
    got = tc.gate_population(pairs)
    assert got["all"] == {"pairs": 4, "regulated": 2, "rate": 0.5}
    assert got["admitted, the gene is the top target"] == {"pairs": 1, "regulated": 1, "rate": 1.0}
    assert got["excluded by the gate"] == {"pairs": 3, "regulated": 1, "rate": 0.3333}


def test_without_the_cache_every_excluded_pair_carries_a_structural_zero(tmp_path):
    pairs, _, _ = _pairs(tmp_path, cache_on=False)
    got = tc.gate_population(pairs)
    assert got["deletion_drop is zero"] == 3  # the three the compact table has nothing to say about
    assert got["the sweep answered this gene"] == 1  # only the top target


def test_with_the_cache_the_censored_pairs_carry_the_sweep_s_own_number(tmp_path):
    pairs, _, _ = _pairs(tmp_path, cache_on=True)
    got = tc.gate_population(pairs)
    assert got["excluded by the gate"]["pairs"] == 3  # the gate itself has not moved
    assert got["the sweep answered this gene"] == 3  # TOP, OTHER and RISES; ABSENT is not in the window
    # RISES is answered and its answer is a rise, so its drop is still zero and now means it
    assert got["deletion_drop is zero"] == 2
    by_gene = {p.gene: p.features["deletion_drop"] for p in pairs}
    assert by_gene == {"TOP": 0.4, "OTHER": 0.25, "RISES": 0.0, "ABSENT": 0.0}


def test_each_excluded_pair_is_a_censoring_or_a_named_silence_and_never_a_zero(tmp_path):
    pairs, table, responses = _pairs(tmp_path, cache_on=False)
    got = tc.gate_silences(pairs, table, responses)
    assert got["by_side"]["admitted"] == {SCORED: 1}
    assert got["by_side"]["excluded"] == {SCORED: 2, NOT_IN_WINDOW: 1}


def test_the_rank_of_the_measured_gene_says_how_far_down_the_window_it_sits(tmp_path):
    pairs, table, responses = _pairs(tmp_path, cache_on=False)
    rank = tc.gate_silences(pairs, table, responses)["measured_gene_rank_in_the_window"]
    assert rank["median_genes_in_the_window"] == 3
    assert rank["pairs"] == 3 and rank["rank_1"] == 1 and rank["top_3"] == 3  # ABSENT has no rank


def test_the_registration_names_the_three_places_the_gate_acts(tmp_path):
    text = tc.PREREGISTERED_GATE
    assert "matched_element" in text and "deletion_drop" in text and "612,323" in text
    assert "IF RELIABILITY IMPROVES" in text  # the outcome that would tempt a lane to stop checking
