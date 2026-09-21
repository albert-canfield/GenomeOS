# SPDX-License-Identifier: AGPL-3.0-or-later
"""Two thresholds that had four and five homes, and what now depends on each.

A duplication sweep found `MIN_EFFECT = 0.1` declared in four modules and the Zoonomia constraint
threshold 2.27 in five. Every copy held the same value, so nothing was wrong — the defect was latent,
and the shape it would have taken is the one this project has already been bitten by: a registration
that interpolates a constant would have gone on quoting its own copy after somebody changed another.

These tests pin the home and name the dependants, so moving either value fails here with a list of
what it moves rather than passing quietly.
"""

from __future__ import annotations

from genomeos.attribution import syntax_tiling, unknown_scoring
from genomeos.attribution.constraint import PHYLOP_THRESHOLD
from genomeos.attribution.lexicon_axes import PHYLOP_THRESHOLD as LEXICON_PHYLOP
from genomeos.knowledge.satmut import PHYLOP_CONSTRAINED
from genomeos.predict import individual_effects
from genomeos.predict.enhancer_target import MIN_EFFECT


def test_min_effect_has_one_home_and_three_dependants() -> None:
    """The deletion scorer owns it; the three that restated it now import it.

    `syntax_tiling` interpolates it into its registration prose and `unknown_scoring` writes it into
    its result as `min_effect_log2`, so a second home would have published two bars under one name.
    """
    assert MIN_EFFECT == 0.1
    assert individual_effects.MIN_EFFECT is MIN_EFFECT
    assert syntax_tiling.MIN_EFFECT is MIN_EFFECT
    assert unknown_scoring.MIN_EFFECT is MIN_EFFECT


def test_the_constraint_threshold_has_one_home() -> None:
    """Zoonomia's 5% FDR bar. `lexicon_axes` still declares its own and is asserted equal instead.

    That module quantises a whole chromosome to one byte per base around this value, so importing it
    would change a stored track's meaning rather than a comparison; equality is the honest pin there,
    and tests/test_lexicon_axes_bins.py holds the rest of that derivation.
    """
    assert PHYLOP_THRESHOLD == 2.27
    assert PHYLOP_CONSTRAINED is PHYLOP_THRESHOLD
    assert LEXICON_PHYLOP == PHYLOP_THRESHOLD


def test_the_argparse_defaults_still_match_the_constant() -> None:
    """Two command lines default to 2.27 as a literal; a drift there is invisible at the call site."""
    import importlib.util
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location("bgw", root / "scripts" / "budget_genome_wide.py")
    assert spec and spec.loader
    text = (root / "scripts" / "budget_genome_wide.py").read_text()

    assert f"default={PHYLOP_THRESHOLD}" in text, "scripts/budget_genome_wide.py's default drifted"
    assert f"default={PHYLOP_THRESHOLD}" in (root / "genomeos" / "cli.py").read_text()
