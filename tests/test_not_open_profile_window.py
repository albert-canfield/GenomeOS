# SPDX-License-Identifier: AGPL-3.0-or-later
"""not_open_profile.py moved onto the window, against its own registered falsifier.

Registered in `genomeos/attribution/onetarget2.py` (24adf33, aed8ae9). This module's falsifier: "the
number of predicted rules changing. One element must yield one rule per gene it names, so the rule
COUNT is expected to rise."

The pin this move must not break: `tests/test_not_open_profile.py` holds `rules()`'s loci to
`context_evidence.rule_loci`'s element for element. `rules()` with its default `responses=None`
emits exactly one compiled rule per element as it always has, so that pin stays true, and the window
rules carry their own source and are in no committed program.
"""

from __future__ import annotations

import pytest

from genomeos.attribution import not_open_profile as nop
from genomeos.predict.enhancer_target import STRONG_EFFECT

needs_the_response_cache = pytest.mark.needs_local_data(
    "data/knowledge/alphagenome/elements",
    how="scripts/enhancer_targets_all.py writes it (worker_scorer, threshold=0.0); "
    "data/knowledge is git-ignored machine-local data, so a fresh checkout cannot have it",
)


def test_the_window_source_is_not_one_of_the_committed_sources():
    """A caller counting over SOURCES must be unchanged by this move."""
    assert nop.SOURCE_PREDICTED_WINDOW not in nop.SOURCES
    assert nop.SOURCES == (nop.SOURCE_PREDICTED, nop.SOURCE_MEASURED)
    assert nop.SOURCE_PREDICTED_WINDOW in nop.ALL_SOURCES


def test_the_band_rule_is_predict_targets_own_rule_at_its_own_threshold():
    """Not a second rule: the same threshold, reproduced, because a window gene has no field."""
    assert nop.window_band(-STRONG_EFFECT) == nop.STRONG
    assert nop.window_band(STRONG_EFFECT) == nop.STRONG
    assert nop.window_band(-(STRONG_EFFECT - 1e-9)) == nop.WEAK
    assert nop.window_band(0.1) == nop.WEAK
    assert nop.window_band(-0.9) == nop.STRONG


@needs_the_response_cache
def test_with_no_reader_the_compiled_rules_are_IDENTICAL_and_no_window_rule_exists():
    """The control, on the real chromosome: the element-for-element pin cannot have moved."""
    plain = nop.rules("chr21")
    assert plain, "chr21 has compiled rules on this machine"
    assert not [r for r in plain if r.source == nop.SOURCE_PREDICTED_WINDOW]
    # and window_rules with no reader returns that same list, opening nothing
    assert nop.window_rules("chr21") == plain


@needs_the_response_cache
def test_the_compiled_rules_survive_the_move_unchanged_rule_for_rule():
    from genomeos.attribution.targets import ElementResponses

    plain = [r for r in nop.rules("chr21") if r.source == nop.SOURCE_PREDICTED]
    windowed = nop.window_rules("chr21", responses=ElementResponses(), coding=_coding("chr21"))
    kept = [r for r in windowed if r.source == nop.SOURCE_PREDICTED]
    assert kept == plain


@needs_the_response_cache
def test_the_registered_falsifier_FIRES_the_predicted_rule_count_rises():
    from genomeos.attribution.targets import ElementResponses

    windowed = nop.window_rules("chr21", responses=ElementResponses(), coding=_coding("chr21"))
    extra = [r for r in windowed if r.source == nop.SOURCE_PREDICTED_WINDOW]
    assert extra, "the falsifier fires: elements name more than one coding gene at the bar"
    # a window rule never repeats its element's own compiled gene
    compiled = {(r.element, r.gene) for r in windowed if r.source == nop.SOURCE_PREDICTED}
    assert not {(r.element, r.gene) for r in extra} & compiled


@needs_the_response_cache
def test_every_window_rule_belongs_to_an_element_that_already_has_a_compiled_rule():
    """A window rule locates a further gene at a KNOWN element; it never invents an element."""
    from genomeos.attribution.targets import ElementResponses

    windowed = nop.window_rules("chr21", responses=ElementResponses(), coding=_coding("chr21"))
    elements = {r.element for r in windowed if r.source == nop.SOURCE_PREDICTED}
    for r in windowed:
        if r.source == nop.SOURCE_PREDICTED_WINDOW:
            assert r.element in elements, r.element


@needs_the_response_cache
def test_a_window_rules_band_agrees_with_the_compact_head_on_the_head_gene():
    """The two band rules must agree where both speak, or the reproduction claim is false."""
    from genomeos.attribution import compile as cp
    from genomeos.attribution.targets import RESULTS_DIR

    checked = 0
    for e in cp._attributed("chr21", RESULTS_DIR):
        pc = e["predicted_coding"]
        stated = pc.get("strength") or nop.WEAK
        assert nop.window_band(pc["log2_fold_change"]) == stated, e["id"]
        checked += 1
    assert checked > 100


def _coding(chrom: str) -> set[str]:
    from pathlib import Path

    from genomeos.genome import Annotation

    ann = Annotation.from_gff3(Path("data/reference") / f"gencode_v50_{chrom}.gff3.gz", {chrom})
    return {g.symbol for g in ann.genes.values() if g.type == "protein_coding"}
