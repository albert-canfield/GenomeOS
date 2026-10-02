# SPDX-License-Identifier: AGPL-3.0-or-later
"""Increment 3 of the response map: the gates it refuses on, and what it may not call agreement."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from genomeos import response_map as rm
from genomeos import response_map2 as rm2
from genomeos import response_map3 as rm3
from genomeos.attribution import eqtl, mpra
from genomeos.attribution import measured as ms

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "data/results/response_map_increment3.json"


class _Files:
    """A `_Files` that always, or never, has the file an assertion names."""

    def __init__(self, available: bool = True) -> None:
        self.available = available
        self.table: dict[str, dict] = {}

    def ref(self, path, key, record_path, **extra):
        if not self.available:
            return None
        return {
            "path": path,
            "record_key": key,
            "record_path": record_path,
            "sha256": "0" * 64,
            "kind": rm.COMMITTED,
            **extra,
        }


def _map() -> rm._Map:
    return rm._Map({"assembly": "GRCh38", "coordinates": {"base": 0, "interval": "half-open"}})


def _candidate(**over):
    pair = SimpleNamespace(chrom="chr1", start=1000, end=1200, gene="HBE1", cell="K562")
    row = {
        "id": "EH1",
        "start": 1000,
        "end": 1200,
        "assays": ["crispri", "lentimpra"],
        "measured": {
            "lentimpra": {
                "label_by_cell": {"K562": "silent", "HepG2": "active", "WTC11": "silent"},
                "tiles": [
                    {
                        "name": "K562_seq1_F",
                        "start": 1000,
                        "end": 1200,
                        "strand": {"K562": "+", "HepG2": "+", "WTC11": "+"},
                        "overlap": 0.9,
                        "activity": {"K562": -0.4, "HepG2": 1.7, "WTC11": -0.9},
                    }
                ],
            }
        },
    }
    base = {
        "pair": pair,
        "status": rm2.STATUS_ADDABLE,
        "element": {"id": "EH1", "start": 1000, "end": 1200},
        "row": row,
        "overlap": 0.9,
        "state": "open_in_reader",
        "biosample": "K562",
        "key": None,
    }
    base.update(over)
    return SimpleNamespace(**base)


# --- the two gates -----------------------------------------------------------------------------------
def test_the_build_refuses_a_population_of_a_different_size_than_the_count():
    rm3.refuse_unless_population_matches(40, 40)
    with pytest.raises(rm.RefusedError, match="would not be over the set that was counted"):
        rm3.refuse_unless_population_matches(39, 40)
    with pytest.raises(rm.RefusedError):
        rm3.refuse_unless_population_matches(41, 40)


def test_the_build_refuses_a_payload_whose_loci_do_not_reconcile():
    ok = {
        "loci_built": 37,
        "loci_excluded_by_cause": {"a": 2, "b": 1},
        "loci_accounted_for": 40,
        "loci_counted_before_the_build": 40,
        "reconciles": True,
    }
    rm3.refuse_unless_reconciled(ok)
    bad = {**ok, "loci_built": 36, "loci_accounted_for": 39, "reconciles": False}
    with pytest.raises(rm.RefusedError, match="not the 40 the count named"):
        rm3.refuse_unless_reconciled(bad)


def test_reconcile_adds_the_built_to_every_cause():
    loci = [{"id": "l1"}, {"id": "l2"}]
    excluded = {c: [] for c in rm3.CAUSES}
    excluded[rm3.CAUSE_NO_FURTHER_KIND] = [{"locus": "l3"}]
    population = {
        "path": "p",
        "sha256": "x",
        "loci_carrying_at_least_one_further_kind": 3,
        "addable_loci": 110,
    }
    got = rm3._reconcile(loci, population, excluded, dict.fromkeys(rm3.CAUSES, 0), {"by_status": {}})
    assert got["loci_accounted_for"] == 3
    assert got["reconciles"] is True
    assert got["share_of_the_addable_loci"] == round(2 / 110, 4)
    # and it breaks the moment a locus leaves without a cause
    got2 = rm3._reconcile(
        loci,
        {**population, "loci_carrying_at_least_one_further_kind": 4},
        excluded,
        dict.fromkeys(rm3.CAUSES, 0),
        {"by_status": {}},
    )
    assert got2["reconciles"] is False


def test_every_cause_increment_2_could_refuse_under_is_still_a_cause():
    assert set(rm2.CAUSES) <= set(rm3.CAUSES)
    assert rm3.CAUSE_REPORTER_SOURCE in rm3.CAUSES
    assert rm3.CAUSE_EQTL_SOURCE in rm3.CAUSES


# --- the reporter observations ------------------------------------------------------------------------
def test_a_reporter_reading_is_of_the_sequence_and_never_of_the_locus():
    m = _map()
    c = _candidate()
    got, cause = rm3.reporter_observations(c, m, _Files(), "interval:chr1:1000-1200", "basis")
    assert cause is None
    assert len(got) == 3  # one tile, three cells: the measured layer's own (tile, cell) unit
    for a in got:
        assert a["status"] == "observed"
        assert a["of_the_native_locus"] is False
        assert a["basis"] == ms.OUTCOME_KIND["lentimpra"]
        assert a["measurement"]["observation_unit"] == list(ms.REPORTER_UNIT)
        assert a["relation"] == "associated_with_state"
        # the object is a reporter-activity state, never the native gene
        assert m.entities[a["object"]]["kind"] == "cell_state"
        assert m.entities[a["object"]]["of_a_reporter_construct"] is True


def test_the_reporter_label_is_the_aggregation_and_not_the_observation():
    m = _map()
    got, _ = rm3.reporter_observations(_candidate(), m, _Files(), "interval:chr1:1000-1200", "basis")
    k562 = next(a for a in got if a["id"].endswith("|K562"))
    assert k562["measurement"]["value"] == -0.4  # the tile's own log2(RNA/DNA): what was measured
    assert k562["aggregated_label"]["label"] == "silent"
    assert k562["aggregated_label"]["is_the_observation"] is False
    assert k562["aggregated_label"]["rule"] == ms.REPORTER_LABEL_RULE
    assert k562["aggregated_label"]["threshold_log2"] == ms.MPRA_ACTIVE


def test_a_reporter_reading_with_no_source_refuses_under_its_own_cause():
    got, cause = rm3.reporter_observations(
        _candidate(), _map(), _Files(available=False), "interval:chr1:1000-1200", "basis"
    )
    assert got == []
    assert cause == rm3.CAUSE_REPORTER_SOURCE


def test_an_element_with_no_reporter_tiles_is_not_a_refusal():
    c = _candidate(row={"id": "EH1", "start": 1000, "end": 1200, "assays": ["crispri"], "measured": {}})
    got, cause = rm3.reporter_observations(c, _map(), _Files(), "interval:chr1:1000-1200", "basis")
    assert got == [] and cause is None


def test_the_reporter_source_is_the_file_of_the_cell_it_was_read_in():
    m = _map()
    got, _ = rm3.reporter_observations(_candidate(), m, _Files(), "interval:chr1:1000-1200", "basis")
    for a in got:
        cell = a["source"]["cell_of_the_file"]
        assert a["source"]["path"] == (mpra.KNOWLEDGE / f"{mpra.FILES[cell]}.bed.gz").as_posix()


# --- the eQTL observations ----------------------------------------------------------------------------
def _hit(**over):
    h = {
        "tissue": "Whole_Blood",
        "chrom": "chr1",
        "pos": 1100,
        "ref": "G",
        "alt": "A",
        "gene_id": "ENSG00000001",
        "slope": -0.39,
        "pval_nominal": 1.6e-07,
        "file": "data/knowledge/gtex/hits_Whole_Blood.tsv",
    }
    h.update(over)
    return h


def test_an_eqtl_is_observed_but_its_relation_stays_unresolved():
    m = _map()
    got, cause = rm3.eqtl_observations(
        _candidate(), m, _Files(), "interval:chr1:1000-1200", "basis", [_hit()], 500
    )
    assert cause is None and len(got) == 1
    a = got[0]
    assert a["status"] == "observed"
    assert a["of_the_native_locus"] is True
    assert a["an_association_not_a_perturbation"] is True
    assert a["relation"] == rm.UNRESOLVED
    assert set(a["relation_candidates"]) <= set(rm.RELATIONS)
    assert len(a["relation_candidates"]) >= 2
    assert a["measurement"]["is_a_perturbation"] is False
    assert a["context"]["perturbation"] == "none: nothing was done to the locus"


def test_the_eqtl_window_is_the_distillation_s_own_and_excludes_its_low_edge():
    m = _map()
    # the element is 1000-1200; with a margin of 500 the window is `1000 - 500 < pos <= 1200 + 500`
    inside = [_hit(pos=501), _hit(pos=1700)]
    outside = [_hit(pos=500), _hit(pos=1701)]
    got, _ = rm3.eqtl_observations(
        _candidate(), m, _Files(), "interval:chr1:1000-1200", "basis", inside + outside, 500
    )
    assert sorted(a["measurement"]["value"] is not None for a in got) == [True, True]
    assert {int(a["id"].split("|")[1].split(":")[1]) for a in got} == {501, 1700}


def test_one_eqtl_assertion_names_one_row_of_one_file():
    m = _map()
    hits = [_hit(), _hit(tissue="Lung", file="data/knowledge/gtex/hits_Lung.tsv")]
    got, _ = rm3.eqtl_observations(_candidate(), m, _Files(), "interval:chr1:1000-1200", "basis", hits, 500)
    assert len(got) == 2  # no aggregation over tissues: no statistic is introduced here
    assert {a["source"]["tissue"] for a in got} == {"Whole_Blood", "Lung"}
    assert len({a["observation_key"] for a in got}) == 2


def test_two_alleles_at_one_position_are_two_assertions():
    m = _map()
    hits = [_hit(ref="G", alt="A"), _hit(ref="G", alt="T")]
    got, _ = rm3.eqtl_observations(_candidate(), m, _Files(), "interval:chr1:1000-1200", "basis", hits, 500)
    assert len({a["id"] for a in got}) == 2


def test_an_eqtl_with_no_source_refuses_under_its_own_cause():
    got, cause = rm3.eqtl_observations(
        _candidate(), _map(), _Files(available=False), "interval:chr1:1000-1200", "basis", [_hit()], 500
    )
    assert got == [] and cause == rm3.CAUSE_EQTL_SOURCE


def test_the_hits_reader_takes_its_columns_from_where_they_were_registered(tmp_path):
    d = tmp_path / "gtex"
    d.mkdir()
    (d / "hits_Lung.tsv").write_text(
        "\t".join(eqtl.HIT_COLUMNS) + "\n" + "Lung\tchr1\t50\tG\tA\tENSG1\t0.1\t1e-9\tchr1:1-100\n"
    )
    got = rm3.hits_for("chr1", d)
    assert got[0]["gene_id"] == "ENSG1" and got[0]["pos"] == 50
    assert rm3.hits_for("chr2", d) == []
    (d / "hits_Bad.tsv").write_text("wrong\tcolumns\n")
    with pytest.raises(rm.RefusedError, match="are not"):
        rm3.hits_for("chr1", d)


# --- observed is split by what measured it -----------------------------------------------------------
def test_observed_is_split_by_what_measured_it_and_nothing_is_pooled():
    got = rm3.observed_by_assay(
        [
            {"id": "p2|x", "status": "observed"},
            {"id": "s2|x", "status": "observed"},
            {"id": "m3|x", "status": "observed"},
            {"id": "q3|x", "status": "observed"},
            {"id": "r2|x", "status": "predicted"},
        ]
    )
    assert len(got) == 4
    assert got["CRISPRi perturbation (ENCODE benchmark)"] == 1
    assert got[mpra.EVIDENCE] == 1
    assert got[eqtl.EVIDENCE] == 1


def test_an_observed_assertion_nothing_can_name_is_refused():
    with pytest.raises(rm.RefusedError, match="nothing says what measured it"):
        rm3.observed_by_assay([{"id": "zz|x", "status": "observed"}])


# --- the statements the map must carry ----------------------------------------------------------------
def test_increment_2s_limit_is_carried_word_for_word():
    assert "23.26% of one of the measured layer's four assays" in rm3.INCREMENT_2_LIMIT
    assert (
        "they read a sequence, or the bases inside it, not the native locus after a perturbation of it"
        in rm3.INCREMENT_2_LIMIT
    )


def test_more_kinds_is_never_more_established():
    assert "not a locus three or four times established" in rm3.MORE_KINDS_IS_NOT_MORE_ESTABLISHED
    joined = " ".join(rm3.CANNOT_ESTABLISH)
    assert "not agreement, confirmation or validation" in joined
    assert "operational grouping" in joined


def test_the_absent_assays_are_named_with_what_they_measure():
    blk = rm3._kinds_block([], __import__("collections").Counter())
    assert ms.OUTCOME_KIND["vista"] in blk["absent_kinds"]["vista"]
    assert ms.SATMUT_CANNOT_DISAGREE in blk["absent_kinds"]["satmut"]
    assert blk["which_kinds_read_the_native_locus"]["lentimpra_sequence"] is False
    assert blk["which_kinds_read_the_native_locus"]["eqtl_native_locus"] is True


# --- the page ----------------------------------------------------------------------------------------
def _payload(n_loci: int = 20) -> dict:
    loci = [
        {"id": f"l{i}", "assertions": [f"a{i}"], "chains": [f"c{i}"], "kinds_of_evidence": []}
        for i in range(n_loci)
    ]
    return {
        "loci": loci,
        "assertions": [{"id": f"a{i}", "entities": [f"e{i}"], "status": "observed"} for i in range(n_loci)],
        "entities": {f"e{i}": {"kind": "dna_interval"} for i in range(n_loci)},
        "chains": [{"id": f"c{i}", "steps": [f"a{i}"]} for i in range(n_loci)],
        "cycles": [],
        "counts": {"assertions": n_loci},
        "reconciliation": {"reconciles": True, "loci_built": n_loci},
        "why": {"question": "q"},
        "conventions": {"assembly": "GRCh38"},
        "vocabulary": {},
        "kinds_of_evidence": {"loci": n_loci},
        "cannot_establish": ["x"],
        "counted": {},
        "population": {},
        "not_open_in_reader": {},
        "sources": [],
        "title": "t",
    }


def test_a_page_carries_the_whole_maps_counts():
    p = rm3.page(_payload(20), 2, per=8)
    assert p["page"] == 2 and p["pages"] == 3
    assert p["loci_total"] == 20 and p["assertions_total"] == 20
    assert len(p["loci"]) == 8
    assert p["reconciliation"]["loci_built"] == 20
    assert p["counts"]["assertions"] == 20


def test_a_page_after_the_first_does_not_repeat_the_preamble():
    one, two = rm3.page(_payload(20), 1, per=8), rm3.page(_payload(20), 2, per=8)
    for block in rm3.DESCRIBES_THE_MAP:
        if block in _payload(20):
            assert block in one
            assert block not in two
    assert set(two["described_on_page_1"]) <= set(rm3.DESCRIBES_THE_MAP)
    assert "not_repeated" in two
    assert "not_repeated" not in one


def test_the_last_page_is_not_a_truncation():
    p = rm3.page(_payload(20), 3, per=8)
    assert len(p["loci"]) == 4
    assert p["loci_total"] == 20
    assert "of 20" in p["paging"]


def test_a_page_past_the_end_clamps_rather_than_emptying():
    p = rm3.page(_payload(20), 99, per=8)
    assert p["page"] == 3 and len(p["loci"]) == 4


# --- the committed map -------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def built():
    if not RESULT.exists():
        pytest.skip("the genome-wide map has not been built here")
    with RESULT.open() as fh:
        return json.load(fh)


@pytest.mark.skipif(not RESULT.exists(), reason="the genome-wide map has not been built here")
class TestTheCommittedMap:
    def test_it_reconciles_against_the_count_committed_before_it(self, built):
        r = built["reconciliation"]
        assert r["reconciles"] is True
        assert (
            r["loci_built"] + sum(r["loci_excluded_by_cause"].values())
            == (r["loci_counted_before_the_build"])
        )
        assert r["counted_sha256"]

    def test_every_assertion_names_its_source_path_record_key_and_sha256(self, built):
        for a in built["assertions"]:
            s = a["source"]
            assert s["path"] and s["record_key"] and s["sha256"]

    def test_observed_is_split_by_four_assays(self, built):
        got = built["counts"]["observed_by_assay"]
        assert len(got) == 4
        assert sum(got.values()) == built["counts"]["by_status"]["observed"]

    def test_no_predicted_quantity_is_observed(self, built):
        for a in built["assertions"]:
            if a["id"].startswith("r2|"):
                assert a["status"] == "predicted"
                assert a["measurement"]["is_a_measurement"] is False

    def test_no_reporter_reading_claims_the_native_locus(self, built):
        n = 0
        for a in built["assertions"]:
            if a["id"].startswith("m3|"):
                n += 1
                assert a["of_the_native_locus"] is False
        assert n > 0

    def test_no_model_request_and_no_per_element_cache_open(self, built):
        assert built["alphagenome_requests"] == 0
        assert built["per_element_response_cache"]["opens"] == 0

    def test_it_says_what_it_cannot_establish(self, built):
        assert len(built["cannot_establish"]) >= 6

    def test_no_verdict_or_score_is_on_the_map(self, built):
        for a in built["assertions"]:
            assert "verdict" not in a
            assert a["kind"] == "evidence"
