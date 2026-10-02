# SPDX-License-Identifier: AGPL-3.0-or-later
"""The DAP-G range read's own rules: imported constants, the matched baseline, the read's bounds."""

from __future__ import annotations

import importlib.util
import pathlib
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load():
    sys.path.insert(0, str(ROOT))
    sys.path.insert(0, str(ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location("_finemap_dapg", ROOT / "scripts/finemap_dapg.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


fd = _load()


# --- nothing quantitative is chosen here ------------------------------------------------------------


def test_every_parameter_is_the_imported_object_itself() -> None:
    """Identity, not equality: a copied number can drift from its source and an equal one hides it."""
    from genomeos.attribution import cell2, executor, fresh, wiring

    assert fd.PIP is executor.DAPG_PIP
    assert fd.CONTROLS_PER_ELEMENT is executor.CONTROLS_PER_UNIT
    assert fd.CONTROL_WINDOW is executor.CONTROL_WINDOW
    assert fd.SEED is executor.SEED
    assert fd.DENSITY_BINS is wiring.ACTIVITY_BINS
    assert fd.DRAW_TRIES is wiring.DERANGE_TRIES
    assert fd.LINK_FLOOR is fresh.POSITIVE_FLOOR
    assert fd.LOCUS_FLOOR is fresh.LOCUS_FLOOR
    assert fresh.LOCUS_FLOOR is cell2.POOLED_LOCUS_FLOOR


def test_the_floors_are_the_pair_lane_reptest_applied() -> None:
    assert (fd.LINK_FLOOR, fd.LOCUS_FLOOR) == (30, 20)


def test_the_margin_comes_through_increment_three() -> None:
    import response_map_increment3_count as c3

    assert fd.MARGIN is c3.EQTL_MARGIN


def test_the_track_is_the_one_the_panel_reads() -> None:
    from genomeos.attribution import human_panel as hp

    assert hp.BIGBEDS[fd.TRACK_KEY] == "https://hgdownload.soe.ucsc.edu/gbdb/hg38/gtex/eQtl/gtexDapg.bb"


def test_the_retained_rows_do_not_land_on_a_pinned_store() -> None:
    """data/knowledge/gtex and gtex_crispri are sha256-pinned inputs of committed results."""
    from genomeos.attribution import eqtl

    assert fd.KNOWLEDGE.resolve() != eqtl.KNOWLEDGE.resolve()
    assert fd.KNOWLEDGE.as_posix() != "data/knowledge/gtex_crispri"


# --- the gene-density covariate ---------------------------------------------------------------------


def test_max_gene_bp_holds_against_the_annotation_on_disk() -> None:
    """The backward scan's bound, checked against real gene bodies rather than trusted."""
    longest = 0
    for chrom in ("chr16", "chr1", "chr2", "chrX"):
        for start, end in fd.gene_bodies(chrom):
            longest = max(longest, end - start)
    assert longest > 0, "no protein-coding bodies were read at all"
    assert longest < fd.MAX_GENE_BP, f"a body of {longest} would be missed by the scan"


def test_density_counts_every_overlap_including_one_starting_far_before() -> None:
    """The case the first draft of this function got wrong: a long gene opened before the window."""
    bodies = [(0, 2_000_000), (1_000, 2_000), (500_000, 500_100), (9_000_000, 9_001_000)]
    bodies.sort()
    starts = [b[0] for b in bodies]
    assert fd.density(bodies, starts, 1_500_000, 1_500_100) == 1
    assert fd.density(bodies, starts, 1_000, 2_000) == 2
    assert fd.density(bodies, starts, 3_000_000, 3_000_100) == 0
    assert fd.density(bodies, starts, 0, 10_000_000) == 4


def test_density_is_half_open_at_both_ends() -> None:
    bodies = [(100, 200)]
    starts = [100]
    assert fd.density(bodies, starts, 200, 300) == 0
    assert fd.density(bodies, starts, 0, 100) == 0
    assert fd.density(bodies, starts, 199, 300) == 1


# --- the baseline -----------------------------------------------------------------------------------


def _toy() -> dict[str, dict[str, object]]:
    return {
        f"E{i}": {
            "chrom": "chr21",
            "start": 10_000_000 + i * 20_000,
            "end": 10_000_250 + i * 20_000,
            "genes": {"G"},
            "statuses": {"all_three_present"},
        }
        for i in range(12)
    }


def test_controls_match_length_exactly() -> None:
    els = _toy()
    ctrl, _info = fd.controls(els)
    for eid, windows in ctrl.items():
        want = els[eid]["end"] - els[eid]["start"]
        for lo, hi in windows:
            assert hi - lo == want, f"{eid}: control length {hi - lo} against element {want}"


def test_controls_stay_inside_the_imported_window_and_on_the_chromosome() -> None:
    els = _toy()
    ctrl, _info = fd.controls(els)
    for eid, windows in ctrl.items():
        s = els[eid]["start"]
        for lo, _hi in windows:
            assert abs(lo - s) <= fd.CONTROL_WINDOW


def test_controls_overlap_no_element_of_the_frame() -> None:
    """A control inside another element would put a real element on the baseline's own side."""
    els = _toy()
    ctrl, _info = fd.controls(els)
    taken = [(v["start"] - fd.MARGIN, v["end"] + fd.MARGIN) for v in els.values()]
    for windows in ctrl.values():
        for lo, hi in windows:
            for a, b in taken:
                assert not (lo - fd.MARGIN < b and hi + fd.MARGIN > a), "a control overlaps an element"


def test_the_draw_is_deterministic_under_the_imported_seed() -> None:
    els = _toy()
    a, ia = fd.controls(els)
    b, ib = fd.controls(els)
    assert a == b
    assert ia["digest"] == ib["digest"]


def test_an_unmatched_element_is_counted_and_not_dropped_silently() -> None:
    els = _toy()
    _ctrl, info = fd.controls(els)
    total = (
        info["elements_with_the_full_number"]
        + info["elements_with_fewer"]
        + info["elements_with_none_unmatched"]
    )
    assert total == len(els), "every element is in exactly one of the three buckets"


def test_the_read_covers_controls_as_well_as_elements() -> None:
    """A baseline whose windows were not read cannot have an excess computed over it."""
    els = _toy()
    ctrl, _info = fd.controls(els)
    iv = fd.intervals_to_read(els, ctrl)
    flat = {(lo, hi) for v in iv.values() for lo, hi in v}
    for eid, windows in ctrl.items():
        for lo, hi in windows:
            assert (max(0, lo - fd.MARGIN), hi + fd.MARGIN) in flat, f"{eid}'s control is not read"
    for v in els.values():
        assert (max(0, v["start"] - fd.MARGIN), v["end"] + fd.MARGIN) in flat


# --- what the registration must say -----------------------------------------------------------------


def test_assessed_by_construction_is_the_registrations_own_wording() -> None:
    import json

    reg = ROOT / "data/results/eqtl_crispri_frame_registration.json"
    if not reg.exists():
        pytest.skip("the re-distillation registration is not on this machine")
    assert json.loads(reg.read_text())["assessed_by_construction"] == fd.ASSESSED_BY_CONSTRUCTION


def test_the_writer_refuses_to_read_without_a_registration_flag() -> None:
    src = (ROOT / "scripts/finemap_dapg.py").read_text()
    assert "--register is required: nothing is read before the registration exists" in src


def test_the_fired_request_bound_is_recorded_and_never_checked_again() -> None:
    """It fired at 225 against 222. It is not relaxed, not amended and not reused as a live check."""
    assert fd.FIRED_REQUEST_BOUND == 222
    assert fd.FIRED_RUN["measured_requests"] == 225
    assert fd.FIRED_RUN["fired_on"] == "requests"
    assert fd.FIRED_RUN["measured_requests"] > fd.FIRED_REQUEST_BOUND, "the record must show it fired"
    src = (ROOT / "scripts/finemap_dapg.py").read_text()
    body = src[src.index("def read(") : src.index("def load_hits(")]
    assert 'cost["requests"] >=' not in body, "requests must not bound the read any more"


def test_the_live_bound_is_bytes_and_there_is_no_request_ceiling() -> None:
    assert fd.BYTE_BOUND_MB == 170.4
    src = (ROOT / "scripts/finemap_dapg.py").read_text()
    body = src[src.index("def read(") : src.index("def load_hits(")]
    assert "total_mb >= bound_mb" in body, "the only cost check is the cumulative byte figure"


def test_the_geometry_reason_is_recorded_not_just_asserted() -> None:
    """Why bytes and not requests, stated as a fact about the two interval sets."""
    assert "DISJOINT intervals" in fd.WHY_BYTES_AND_NOT_REQUESTS
    assert "total span" in fd.WHY_BYTES_AND_NOT_REQUESTS
    assert "not a rationalisation" in fd.WHY_BYTES_AND_NOT_REQUESTS


def test_the_wrong_quantity_is_attributed_to_the_mis_citation() -> None:
    """Not to the panel results, which claim nothing about which quantity bounds another read."""
    assert "mis-citation" in fd.WHOSE_MISTAKE_THE_WRONG_QUANTITY_WAS
    assert "not the precedent's" in fd.WHOSE_MISTAKE_THE_WRONG_QUANTITY_WAS


# --- the count refuses a partial read --------------------------------------------------------------


def test_the_count_refuses_a_read_that_did_not_finish(tmp_path, monkeypatch) -> None:
    """Counting over the chromosomes that happen to sort first is the error this refusal exists for."""
    monkeypatch.setattr(fd, "KNOWLEDGE", tmp_path)
    (tmp_path / "hits_chr1.tsv").write_text("chrom\tpos\tvariant\tgene\ttissue\tpip\n")
    (tmp_path / "hits_chr2.tsv").write_text("chrom\tpos\tvariant\tgene\ttissue\tpip\n")
    with pytest.raises(fd.RefusedError) as e:
        fd.read_is_complete()
    assert "did not finish" in str(e.value)
    assert "chr1, chr2" in str(e.value), "the refusal names which chromosomes it has"


def test_the_completeness_check_runs_before_any_frame_work(tmp_path, monkeypatch) -> None:
    """A refusal that costs three minutes of element tables to reach is one someone works around."""
    monkeypatch.setattr(fd, "KNOWLEDGE", tmp_path)
    (tmp_path / "hits_chr1.tsv").write_text("chrom\tpos\tvariant\tgene\ttissue\tpip\n")

    def explode() -> None:  # pragma: no cover - must never be reached
        raise AssertionError("the frame was read before the read was known to be complete")

    monkeypatch.setattr(fd, "frame_with_keys", explode)
    monkeypatch.setattr(fd, "controls", explode)
    with pytest.raises(fd.RefusedError):
        fd.count_payload()


def test_read_is_complete_returns_the_summary_when_it_is_there(tmp_path, monkeypatch) -> None:
    import json

    monkeypatch.setattr(fd, "KNOWLEDGE", tmp_path)
    (tmp_path / "hits_chr1.tsv").write_text("chrom\tpos\tvariant\tgene\ttissue\tpip\n")
    (tmp_path / "read_summary.json").write_text(json.dumps({"intervals": 7}))
    assert fd.read_is_complete()["intervals"] == 7


def test_the_cost_bound_is_checked_after_every_chromosome_not_at_the_end() -> None:
    """A bound only checked once the work is done is not a bound."""
    src = (ROOT / "scripts/finemap_dapg.py").read_text()
    body = src[src.index("def read(") : src.index("def load_hits(")]
    assert body.index("total_mb >= bound_mb") < body.index("summary = {")


def test_the_read_refuses_without_the_new_cost_registration() -> None:
    """The fired bound is not reusable, so reading further needs the new registration to exist."""
    src = (ROOT / "scripts/finemap_dapg.py").read_text()
    body = src[src.index("def read(") : src.index("def load_hits(")]
    assert "stands fired; it is not reused and not" in body
    assert body.index("COST_REGISTRATION") < body.index("for chrom in sorted(iv)")


def test_the_resume_skips_chromosomes_already_read() -> None:
    src = (ROOT / "scripts/finemap_dapg.py").read_text()
    body = src[src.index("def read(") : src.index("def load_hits(")]
    assert "already read, not fetched again" in body
    assert "carried_mb" in body, "the byte bound must span both runs, not just this one"


# --- the exposure statement -------------------------------------------------------------------------


def test_the_exposure_answers_both_questions_asked() -> None:
    """Was a carrying figure seen? Was an element-versus-baseline figure seen? Both must be answered."""
    e = fd.THE_EXPOSURE
    assert e["the_answer"] == "no carrying or baseline figure was computed on the partial read"
    assert e["was_a_carrying_figure_computed_or_seen"].startswith("No")
    assert e["was_an_element_versus_baseline_figure_computed_or_seen"].startswith("No")


def test_the_exposure_claim_is_true_of_the_code_it_cites() -> None:
    """The evidence is a mechanism, so the mechanism is checked rather than the sentence describing it."""
    src = (ROOT / "scripts/finemap_dapg.py").read_text()
    body = src[src.index("def read(") : src.index("def load_hits(")]
    for forbidden in (
        "symbols_of",
        "rows_in",
        "cell2.group",
        "symbol_map",
        "LINK_FLOOR",
        "LOCUS_FLOOR",
        "fisher_greater",
    ):
        assert forbidden not in body, (
            f"the exposure statement says {forbidden} appears nowhere in read(), and it does"
        )


def test_the_exposure_discloses_the_prior_source_and_its_figures() -> None:
    """14 and 0 came from the on-disk record before this read; a reader must not have to infer that."""
    e = fd.THE_EXPOSURE["prior_exposure_from_a_DIFFERENT_source_disclosed"]
    assert "finemap_coverage.json" in e
    assert "14 of 1,505" in e and "0 carrying one" in e
    assert "before this read was registered" in e


def test_the_exposure_says_what_the_prior_figures_could_not_have_tuned() -> None:
    e = fd.THE_EXPOSURE["what_that_prior_exposure_could_not_have_tuned"]
    assert "fresh.POSITIVE_FLOOR" in e and "fresh.LOCUS_FLOOR" in e
    assert "before this lane existed" in e


def test_the_bound_change_is_stated_to_be_independent_of_the_rows() -> None:
    e = fd.THE_EXPOSURE["did_the_bound_change_depend_on_anything_the_read_returned"]
    assert e.startswith("No")
    assert "knowable before a single row came back" in e


def test_registered_terms_are_carried_forward_not_recomputed() -> None:
    """Regenerating after the read turned `the_resume` from 11 pending chromosomes into none."""
    fresh = {
        "the_new_bound": {"bound_mb": 170.4},
        "the_resume": {"chromosomes_still_to_read": []},
        "the_exposure": {"the_answer": "no carrying or baseline figure was computed on the partial read"},
    }
    # the function must keep the registered value and name the drift
    out = _carry(fd, fresh)
    assert out["the_resume"]["chromosomes_still_to_read"] == ["chr20", "chr3"], "registered value kept"
    assert out["the_exposure"] == fresh["the_exposure"], "a genuinely new key is added"
    drift = out["terms_that_would_have_moved_on_regeneration"]
    assert drift["keys"] == ["the_resume"]
    assert drift["detail"]["the_resume"]["kept"] == "as_registered"


def _carry(mod, fresh):
    """Run carry_forward_registered_terms against a temporary committed file."""
    import json
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        results = pathlib.Path(d) / "data" / "results"
        results.mkdir(parents=True)
        (results / f"{mod.COST_REGISTRATION}.json").write_text(
            json.dumps(
                {
                    "the_new_bound": {"bound_mb": 170.4},
                    "the_resume": {"chromosomes_still_to_read": ["chr20", "chr3"]},
                }
            )
        )
        real = mod.ROOT
        mod.ROOT = pathlib.Path(d)
        try:
            return mod.carry_forward_registered_terms(fresh)
        finally:
            mod.ROOT = real


def test_carry_forward_is_a_no_op_before_the_first_commit() -> None:
    """With no committed file there is nothing to preserve and the fresh payload stands."""
    import tempfile

    fresh = {"a": 1}
    real = fd.ROOT
    with tempfile.TemporaryDirectory() as d:
        fd.ROOT = pathlib.Path(d)
        try:
            assert fd.carry_forward_registered_terms(fresh) == fresh
        finally:
            fd.ROOT = real
