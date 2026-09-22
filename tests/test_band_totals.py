"""The genome-wide band totals, guarded as an arithmetic identity rather than as a remembered figure.

Written 2026-09-22, after "the curve bands all 612,323 targets" travelled into four documents -
ATTRIBUTION.md twice, LESSONS.md, ROADMAP.md and the module's own prose - and was wrong in each.
612,323 is the count of targets the sweep NAMES. 593,765 carry a band; the other 18,558 have no
band because their predicted gene has no GENCODE v50 TSS. Nothing failed when the figure drifted,
because nothing was watching the relation between the two numbers: each document quoted a number
that was true of something, and none of them said which thing.

So these tests pin the IDENTITIES and not the values. A re-banding is allowed to move every count
in the file - the union-axis lane may do exactly that - and is not allowed to leave the parts
disagreeing with the whole, which is the failure that actually occurred.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "data/results/target_calibration.json"

pytestmark = pytest.mark.skipif(not RESULT.exists(), reason="the calibration sweep has not been run here")


def genome_wide() -> dict:
    return json.loads(RESULT.read_text())["genome_wide"]


@pytest.mark.parametrize("arm", ["any gene", "coding gene"])
def test_the_bands_account_for_every_banded_target(arm: str) -> None:
    """`targets` is the sum of the band counts, so a band cannot be dropped without the total moving."""
    g = genome_wide()[arm]
    assert sum(g["predicted_target_bands"].values()) == g["targets"]


@pytest.mark.parametrize("arm", ["any gene", "coding gene"])
def test_the_drop_strata_account_for_every_banded_target(arm: str) -> None:
    """The axis the re-banding of 2026-09-22 registered partitions the same targets, not a subset."""
    g = genome_wide()[arm]
    assert sum(g["drop_bands"].values()) == g["targets"]


def test_banded_plus_unbanded_is_the_number_the_sweep_names() -> None:
    """The identity the four documents lost: banded + no-TSS = named.

    A document is free to quote either number; it is not free to quote one and mean the other, and
    this is the only place the two are written down together.
    """
    g = genome_wide()["any gene"]
    named = g["targets"] + g["missing"]["no_gencode_tss_for_the_predicted_gene"]
    assert named == 612323, "the named-target count moved; every document quoting it must move too"
    assert g["targets"] < named, "a banded count can never exceed the count of targets that exist"


def test_the_coding_arm_is_inside_the_any_gene_arm() -> None:
    """Coding targets are a subset, so every one of its totals is bounded by the wider arm's."""
    g = genome_wide()
    assert g["coding gene"]["targets"] <= g["any gene"]["targets"]


def test_the_registry_classes_do_not_exceed_the_arm_they_cut() -> None:
    """The class cut is a partition of a subset: dELS and pELS together cannot outnumber the arm."""
    g = genome_wide()["any gene"]
    by_class = sum(sum(bands.values()) for bands in g["by_registry_class"].values())
    assert by_class <= g["targets"]
