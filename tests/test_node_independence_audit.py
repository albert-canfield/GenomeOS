# SPDX-License-Identifier: AGPL-3.0-or-later
"""The independence audit's arithmetic and its caller table, pinned.

The audit's claims about history are `git log` output and its claim about the import closure is a
fresh interpreter, so what is left to test is the arithmetic that decides "no usable arm exists"
and the fact that the caller set it scores is the caller set the default was chosen out of.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

_spec = importlib.util.spec_from_file_location("nia", ROOT / "scripts/node_independence_audit.py")
nia = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(nia)


def test_binomial_tail_is_the_exact_tail() -> None:
    assert nia._binom_tail(1, 1, 0.5) == 0.5
    assert abs(nia._binom_tail(2, 0, 0.3) - 1.0) < 1e-12
    assert abs(nia._binom_tail(19, 17, 0.7) - 0.04628) < 1e-4


def test_power_at_nineteen_pairs_cannot_see_six_points() -> None:
    """The MHC screen's 19 regulated pairs, against the claim's own shares."""
    p = nia.power_at(19, 0.699, 0.758)
    assert p["critical_count"] == 17  # 17 of 19, that is 89.5%, to reject at one-sided 0.05
    assert p["power"] < 0.2


def test_power_rises_with_pairs_and_needs_hundreds() -> None:
    small = nia.power_at(36, 0.699, 0.758)["power"]
    big = nia.power_at(400, 0.699, 0.758)["power"]
    assert small < big
    assert nia.power_at(661, 0.699, 0.758)["power"] > 0.8
    need = next(n for n in range(10, 2000) if nia.power_at(n, 0.699, 0.758)["power"] >= 0.8)
    assert 200 < need < 600


def test_the_caller_set_is_the_set_the_default_was_chosen_out_of() -> None:
    """Every candidate scripts/oriented_domains.py compared, minus its declared diagnostic."""
    import oriented_domains as od

    considered = [c for c in od.CALLERS if not c.endswith("0.90")]
    assert sorted(nia.CALLERS) == sorted(considered)
    for name, (rel, best) in od.SITE_CALLS.items():
        if name in nia.CALLERS:
            got = nia.CALLERS[name]
            assert got[:2] == (rel, best)
            assert got[2] == (name == "oriented_ctcf_only")


def test_the_default_caller_carries_no_orientation() -> None:
    assert nia.CALLERS["ctcf_only"] is None


def test_every_component_row_carries_a_legend_status() -> None:
    legend = {"CLEAN", "EXPOSED", "EVALUATION-ONLY"}
    rows = nia.components(
        {"genomeos_modules_imported": [], "modules_naming_a_crispri_source": [], "clean": True},
        {
            "first_commit_naming_the_crispri_benchmark": "x 2026-09-16",
            "constants": dict.fromkeys(
                (
                    "min_node_bp",
                    "boundary_merge_bp",
                    "boundary_class",
                    "node_confidence",
                    "control_seed",
                    "control_draws",
                ),
                {"written": "2026-09-10"},
            ),
        },
        {"what_it_shows": "w", "chromosomes_compared": 23, "identical": True},
        {"total": 661, "pairs_per_file": {}},
    )
    assert {r["status"] for r in rows} <= legend
    assert all(r["why"] for r in rows)
    assert sum(1 for r in rows if r["status"] == "EXPOSED") == 1
