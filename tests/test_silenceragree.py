# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Guards for the silencer_like independent check.

Each one is planted against a specific way this study could have been wrong, and each was verified
to FAIL against the body it guards before being kept.
"""

from __future__ import annotations

import json

import pytest

from genomeos.attribution import silenceragree as sa
from scripts import silenceragree_register as reg


def _record(description: str, gi_chrom: str, gi_start: int, gi_stop: int, **extra: object) -> dict:
    return {
        "description": description,
        "genomicinfo": [{"chrloc": gi_chrom, "chrstart": gi_start, "chrstop": gi_stop}],
        "summary": "identified as a functional silencer in K562 erythroleukemia cells.",
        **extra,
    }


# ---------------------------------------------------------------- the GRCh37 trap


def test_the_join_uses_grch38_genomicinfo_and_never_the_grch37_description() -> None:
    """The ReSE interval is stated in GRCh37 in the description and GRCh38 in genomicinfo.

    Joining on the description would put every fragment at the wrong GRCh38 address. The record here
    makes the two disagree by a wide margin so a reader of the wrong field cannot land by luck.
    """
    kept, _ = sa.rese_intervals(
        {"1": _record("ReSE screen-validated silencer GRCh37_chr7:2272476-2272854", "7", 2232840, 2233218)}
    )
    assert len(kept) == 1
    assert (kept[0]["start"], kept[0]["end"]) == (2232840, 2233218), "joined on the GRCh37 description"
    assert kept[0]["grch37"] == "chr7:2272476-2272854"


def test_a_record_naming_rese_only_outside_the_description_is_rejected() -> None:
    """The 480 records whose genomicinfo is a PARENT region, median 3.55x the assayed fragment.

    Accepting them would attribute a silencer assay to sequence the assay never carried.
    """
    rec = _record("solute carrier family 6 member 4 gene promoter", "17", 30235480, 30237520)
    rec["otherdesignations"] = "ReSE screen-validated silencer GRCh37_chr17:28562949-28563131"
    kept, accounting = sa.rese_intervals({"1": rec})
    assert kept == []
    assert accounting["rejected"]["not_the_assayed_fragment"] == 1


def test_a_record_unplaced_on_grch38_is_rejected_not_placed_at_chrom_blank() -> None:
    """Two records carry chrloc '', and a reader that trusted it would invent chromosome 'chr'."""
    kept, accounting = sa.rese_intervals(
        {"1": _record("ReSE screen-validated silencer GRCh37_chr1:143460194-143460356", "", 1, 2)}
    )
    assert kept == []
    assert accounting["rejected"]["unplaced_on_grch38"] == 1


# ---------------------------------------------------------------- scope sees no direction


def test_element_index_returns_no_direction_at_all() -> None:
    """scope()'s only door onto the element results must hand it no action field.

    If a direction ever leaks into this tuple, a registration built from scope() could carry an
    outcome, which is the one thing the registration exists to make impossible.
    """
    index = sa.element_index("chr21")
    assert index, "chr21 is not scored; this guard needs the local results"
    for start, end, eid, named in index[:500]:
        assert isinstance(start, int) and isinstance(end, int)
        assert isinstance(eid, str) and isinstance(named, bool)
    flat = {repr(x) for row in index[:2000] for x in row}
    assert not {"'represses'", "'activates'"} & flat


def test_scope_output_holds_no_field_a_run_would_fill() -> None:
    assert reg.forbidden_at_any_depth(sa.scope()) == []


def test_the_forbidden_check_is_recursive_not_only_top_level() -> None:
    """The earlier registrations checked the top level only; a figure one dict down would pass."""
    assert reg.forbidden_at_any_depth({"scope": {"deep": {"agreement_rate": 61.2}}}) == [
        "scope.deep.agreement_rate"
    ]
    assert reg.forbidden_at_any_depth({"a": [{"band": "CONSISTENT"}]}) == ["a[0].band"]


# ---------------------------------------------------------------- the control


def test_the_control_never_draws_an_element_from_the_test_set() -> None:
    """A control that could re-draw the test set would pull its own mean toward the observation."""
    kept, _ = sa.rese_intervals(sa.rese_records())
    rows, _ = sa.join(kept)
    used = {r["element"] for r in rows}
    for drawn in sa.control_draws(rows, 3, sa.CONTROL["seed"]):
        assert not used & set(drawn)


def test_the_control_matches_the_test_set_stratum_by_stratum() -> None:
    """Matched on (chromosome, cCRE class). A control matched on chromosome alone would import the
    genome's class composition, and the join is dELS and pELS only."""
    from collections import Counter

    kept, _ = sa.rese_intervals(sa.rese_records())
    rows, _ = sa.join(kept)
    want = Counter((r["chrom"], r["el_class"]) for r in rows)
    classes = {c: sa.ccre_classes(c) for c in {r["chrom"] for r in rows}}
    where = {}
    for chrom, table in classes.items():
        for eid, cls in table.items():
            where[eid] = (chrom, cls)
    for drawn in sa.control_draws(rows, 2, sa.CONTROL["seed"]):
        assert Counter(where[e] for e in drawn) == want
        assert len(set(drawn)) == len(drawn), "a control set drew the same element twice"


def test_the_control_is_deterministic_under_its_registered_seed() -> None:
    kept, _ = sa.rese_intervals(sa.rese_records())
    rows, _ = sa.join(kept)
    a = sa.control_draws(rows, 2, sa.CONTROL["seed"])
    b = sa.control_draws(rows, 2, sa.CONTROL["seed"])
    assert a == b
    assert sa.control_draws(rows, 2, sa.CONTROL["seed"] + 1) != a


# ---------------------------------------------------------------- the bands


@pytest.mark.parametrize(
    ("diff", "low", "high", "want"),
    [
        (25.0, 18.0, 31.0, "CONSISTENT"),
        (10.0, 1.0, 19.0, "CONSISTENT"),
        (9.9, 1.0, 19.0, "DISPLACED_BUT_SMALL"),
        (3.0, -1.0, 7.0, "NO_SIGNAL"),
        (-25.0, -31.0, -18.0, "INVERTED"),
        (-10.0, -19.0, -1.0, "INVERTED"),
        (-9.9, -19.0, -1.0, "DISPLACED_BUT_SMALL"),
        (0.0, -1.0, 1.0, "NO_SIGNAL"),
    ],
)
def test_the_registered_bands_decide_exactly_as_registered(
    diff: float, low: float, high: float, want: str
) -> None:
    assert sa.classify(diff, low, high, clusters=23) == want


def test_an_interval_touching_zero_is_not_consistent() -> None:
    """The boundary the band rule turns on: an interval whose end IS zero has not excluded zero."""
    assert sa.classify(40.0, 0.0, 80.0, clusters=23) == "INCONCLUSIVE_UNDERPOWERED"
    assert sa.classify(40.0, 0.001, 80.0, clusters=23) == "CONSISTENT"


def test_an_interval_spanning_zero_and_the_threshold_is_underpowered_not_no_signal() -> None:
    """The reading that makes an underpowered check a complete answer instead of a false negative.

    An interval from -1 to +51 contains both 'no signal' and 'consistent', so the population cannot
    tell them apart. Calling that NO_SIGNAL would report an absence the data cannot support.
    """
    assert sa.classify(25.0, -1.0, 51.0, clusters=23) == "INCONCLUSIVE_UNDERPOWERED"
    assert sa.classify(3.0, -1.0, 7.0, clusters=23) == "NO_SIGNAL"
    assert sa.classify(0.0, -12.0, 4.0, clusters=23) == "INCONCLUSIVE_UNDERPOWERED"


def test_fewer_than_ten_clusters_reports_no_band_at_all() -> None:
    """Elements near one another are not independent draws, so below 10 clusters there is no
    interval and therefore no band -- a naive one is not offered in its place."""
    assert sa.classify(40.0, 18.0, 62.0, clusters=9) == "INCONCLUSIVE_TOO_FEW_CLUSTERS"
    assert sa.classify(40.0, 18.0, 62.0, clusters=10) == "CONSISTENT"


def test_the_interval_is_clustered_and_refuses_to_report_below_ten_clusters() -> None:
    few = {f"chr{i}": (40, 100) for i in range(1, 10)}
    out = sa.clustered_interval(few, 37.7, 200, 1)
    assert out["interval_percent_points"] is None
    assert out["clusters"] == 9
    identical = {f"chr{i}": (40, 100) for i in range(1, 24)}
    out = sa.clustered_interval(identical, 37.7, 200, 1)
    assert out["clusters"] == 23
    assert out["interval_percent_points"] == [2.3, 2.3], (
        "clusters that all hold the same rate leave the bootstrap no variation to find"
    )
    varied = {f"chr{i}": (30 + i, 100) for i in range(1, 24)}
    low, high = sa.clustered_interval(varied, 37.7, 500, 1)["interval_percent_points"]
    assert low < high, "clusters that disagree must give a non-degenerate interval"


def test_the_clustered_interval_is_wider_than_a_naive_one_on_clumped_data() -> None:
    """The whole reason the interval is clustered: between-cluster variation has to show up.

    Twenty-three chromosomes that disagree sharply with one another must give a wider interval than
    twenty-three that agree, at the same pooled rate and the same denominator.
    """
    agree = {f"chr{i}": (50, 100) for i in range(1, 24)}
    clumped = {f"chr{i}": ((100 if i % 2 else 0), 100) for i in range(1, 24)}
    wa = sa.clustered_interval(agree, 50.0, 500, 3)["interval_percent_points"]
    wc = sa.clustered_interval(clumped, 50.0, 500, 3)["interval_percent_points"]
    assert (wc[1] - wc[0]) > (wa[1] - wa[0])


# ---------------------------------------------------------------- prior exposure and independence


def test_the_registration_records_prior_exposure_as_attestation_not_audit() -> None:
    """The probe lane committed no artefact, so its not having computed an agreement rate cannot be
    confirmed from the tree. The registration must say that in those terms, not claim blindness."""
    text = sa.PRIOR_EXPOSURE["so_is_this_test_blind"]
    assert "BLIND BY ATTESTATION, NOT BLIND BY AUDIT" in text
    assert "CONFIRM OR REFUTE" in text


def test_independence_from_the_training_data_is_registered_as_not_established() -> None:
    """Assuming independence would be worse than finding nothing, because the number gets quoted."""
    assert sa.INDEPENDENCE["status"].startswith("NOT ESTABLISHED")
    assert "what_cannot_be_established" in sa.INDEPENDENCE
    assert sa.INDEPENDENCE["what_cannot_be_established"].startswith("whether the ReSE screen's data")


# ---------------------------------------------------------------- the negative that stands


def test_no_rese_record_validates_in_gm12878_or_imr90() -> None:
    """The four-cell-line negative, measured on the records rather than asserted.

    The project scores GM12878, HepG2, IMR-90 and K562. The source covers two of them, so 'silencers,
    genome-wide, all four lines' is not fundable by any source that exists.
    """
    kept, _ = sa.rese_intervals(sa.rese_records())
    seen = {c for item in kept for c in item["cells"]}
    assert seen <= {"HepG2", "K562"}
    assert "GM12878" not in seen
    assert "IMR-90" not in seen
    assert set(sa.PROJECT_CELL_LINES) - seen == {"GM12878", "IMR-90"}


def test_the_registration_as_written_holds_no_field_a_run_would_fill() -> None:
    """The committed registration file itself, not just the payload the writer built."""
    from pathlib import Path

    path = Path("data/results/silenceragree_registration.json")
    if not path.exists():
        pytest.skip("the registration is not written yet")
    written = json.loads(path.read_text())
    assert written["result"] == "silenceragree_registration"
    assert [k for k in reg.forbidden_at_any_depth(written) if k != "result"] == []


def test_the_base_rate_excludes_the_joined_elements_it_is_the_base_rate_for() -> None:
    """A base rate containing the observation is not a base rate.

    The joined elements are what is being measured; leaving them inside the stratum population pulls
    the base rate toward the observed rate and shrinks the difference the bands turn on. The
    invariant: the stratum population plus the distinct joined elements must account for every
    named-target element in those strata, with no element in both.
    """
    kept, _ = sa.rese_intervals(sa.rese_records())
    rows, _ = sa.join(kept)
    base = sa.stratum_matched_base_rate(rows)
    joined = {r["element"] for r in rows}
    strata = {(r["chrom"], r["el_class"]) for r in rows}
    total = 0
    for chrom in sorted({c for c, _ in strata}):
        classes = sa.ccre_classes(chrom)
        for eid in sa.directions(chrom):
            if (chrom, classes.get(eid, "unknown")) in strata:
                total += 1
    assert base["elements"] == total - len(joined), (
        "the stratum base-rate population still contains the joined elements it is the base rate for"
    )
    assert base["elements"] < total
