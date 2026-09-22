# SPDX-License-Identifier: AGPL-3.0-or-later
"""E1's extension: a declared look schedule, and a fresh sample rather than a continuation."""

from genomeos.attribution import executor as ex


def _pairs(n: int, endpoint: str = "E1W_mpra_wide") -> list[dict]:
    """n pairs whose scorer always agrees, so any success look that can fire, will."""
    return [
        {
            "endpoint": endpoint,
            "unit": f"chr1:{i}-{i + 200}",
            "test": {"pos": 1001 + 2 * i, "ref": "A", "alt": "G"},
            "control": {"unit": f"chr1:c{i}", "pos": 2000 + 2 * i, "ref": "C", "alt": "T"},
            "borrowed": {"gene": "GENE1", "tissue": "Liver", "cell": "K562", "measured_sign": 1},
            "order": i,
        }
        for i in range(n)
    ]


def _scorer(chrom, pos, ref, alt):
    """Units move the way the reporter measured; their controls move the other way.

    The test positions are odd and the control positions even (see `_pairs`), so this is the shape a
    run that is working looks like: agreement at the units, none at the matched nulls.
    """
    return [("GENE1", "Liver", 0.5 if pos % 2 else -0.5)]


def test_the_extension_declares_one_futility_look_and_no_success_look() -> None:
    """Alpha 0 at the interim is how 'no success look' is written: no p can meet it."""
    looks = ex.E1_EXTENSION["alpha_spending"]

    assert looks["look_1_pairs"] == 1130
    assert looks["alpha_1"] == 0.0
    assert looks["look_2_pairs"] == 0  # never reached: the first pair is n = 1
    assert looks["final_alpha"] == 0.01


def test_a_declared_schedule_is_the_one_consulted(tmp_path) -> None:
    """The run must read the schedule it was given, not the project-wide one.

    Checked by watching the reads rather than by engineering a stop: whether a look fires is
    `decide`'s business and is tested with `decide`, but which schedule the run asks is this wiring.
    """

    class Watched(dict):
        def __init__(self, *a, **k):
            super().__init__(*a, **k)
            self.read: list[str] = []

        def __getitem__(self, key):
            self.read.append(key)
            return super().__getitem__(key)

    looks = Watched(ex.E1_EXTENSION["alpha_spending"])

    ex.run_pairs(
        _pairs(3),
        _scorer,
        cache=tmp_path,
        looks_at=("E1W_mpra_wide",),
        endpoints=("E1W_mpra_wide",),
        split="cell",
        looks=looks,
    )

    assert "look_1_pairs" in looks.read
    assert "final_alpha" in looks.read


def test_an_alpha_of_zero_cannot_stop_a_run_for_success(tmp_path) -> None:
    """The extension's interim look must be futility-only however good the data look."""
    pairs = _pairs(6)

    out = ex.run_pairs(
        pairs,
        _scorer,
        cache=tmp_path,
        looks_at=("E1W_mpra_wide",),
        endpoints=("E1W_mpra_wide",),
        split="cell",
        looks={"look_1_pairs": 2, "alpha_1": 0.0, "look_2_pairs": 0, "alpha_2": 0.0, "final_alpha": 0.01},
    )

    assert out["pairs_done"] == 6  # ran to the end of the list rather than stopping at the look
    assert "success" not in str(out["stopped"])


def test_the_default_schedule_is_still_the_projects(tmp_path) -> None:
    """Passing no schedule must not change what every earlier run did."""
    out = ex.run_pairs(_pairs(3), _scorer, cache=tmp_path, looks_at=(), endpoints=("E1W_mpra_wide",))

    assert out["pairs_done"] == 3


def test_the_skip_is_exactly_what_the_first_run_scored() -> None:
    """The extension's sample is the complement of cbff5f6's 1,000, not a new selection."""
    assert ex.E1_EXTENSION_SKIP == 1000


def test_the_registration_states_that_pooling_is_not_a_test() -> None:
    """The whole point of the document: the two samples may be described together, never pooled."""
    assert "never be used as a test" in ex.E1_EXTENSION["why"] or "never used as a test" in str(
        ex.E1_EXTENSION
    )
    assert "weaker by construction" in ex.E1_EXTENSION["sample"]
