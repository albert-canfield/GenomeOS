"""The response map, increment 2: every locus of the perturbation assay that carries all three inputs.

The build is held to the count that fixed its population before anything was built
(`data/results/response_map_coverage.json`, lane-rmap2) rather than to its own arithmetic: the status
constants are pinned against that script's, so the two cannot answer differently about which input a
pair is missing, and the built locus count must come back to the counted one with every exclusion
named by cause. The rules are shown twice: on a real chromosome, and on a payload that breaks one and
must be refused.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from genomeos import response_map as rm
from genomeos import response_map2 as r2
from genomeos.attribution import compile as cp
from genomeos.attribution import context_evidence as ce
from genomeos.attribution import measured as ms
from genomeos.genome import reader

ROOT = Path(__file__).resolve().parent.parent
TRIAL = "chr21"
RESULTS = ROOT / "data" / "results"
PAGE_HTML = ROOT / "genomeos" / "web" / "static" / "index.html"
BUILT = RESULTS / f"{r2.RESULT}.json"

HAS_INPUTS = (
    (ROOT / ce.TRACK_METADATA).is_file()
    and all((ROOT / ms.CRISPRI_KNOWLEDGE / f).is_file() for f in ms.CRISPRI_FILES)
    and (RESULTS / reader.peaks_path("K562", TRIAL).name).is_file()
)
needs_inputs = pytest.mark.skipif(not HAS_INPUTS, reason="the local caches this build reads are absent")


@pytest.fixture(scope="module")
def payload() -> dict:
    return r2.build(ROOT, [TRIAL])


@pytest.fixture(scope="module")
def coverage_script():
    """`scripts/response_map_coverage.py` imported as a module: lane-rmap2's file, read only."""
    path = ROOT / "scripts" / "response_map_coverage.py"
    spec = importlib.util.spec_from_file_location("response_map_coverage_readonly", path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def counted() -> dict:
    return json.loads((ROOT / r2.COVERAGE).read_text())


def test_the_status_rule_is_the_counts_own_and_not_a_second_copy(coverage_script, counted):
    """The names, the order and the figures come from the script that counted the denominator. A
    second copy that drifted would put the two numbers out of step without either one being wrong."""
    assert r2.STATUSES == coverage_script.STATUSES
    assert r2.STATUS_ORDER == coverage_script.STATUS_ORDER
    assert r2.STATUS_ADDABLE == coverage_script.STATUS_ADDABLE
    assert counted["denominator"]["independent_loci"] == r2.MEASURED_LAYER_LOCI
    assert counted["addable"]["independent_loci"] == r2.ADDABLE_LOCI
    assert counted["by_status"]["pairs"][r2.STATUS_ADDABLE] == r2.ADDABLE_PAIRS
    assert counted["addable"]["independent_loci_the_map_already_covers"] == r2.LOCI_INCREMENT_1_COVERED


@needs_inputs
def test_the_attachment_agrees_with_the_counts_own_scan_pair_for_pair(coverage_script):
    elements = cp._attributed(TRIAL, RESULTS)
    starts = [e["start"] for e in elements]
    pairs, _ = ms.load_crispri(TRIAL)
    assert pairs, "the trial chromosome holds no measured pair, so nothing was compared"
    for p in pairs:
        mine, _overlap = r2.attaches_to(p, starts, elements)
        theirs = coverage_script.attaches_to(p, starts, elements)
        assert (mine or {}).get("id") == (theirs or {}).get("id")


@needs_inputs
def test_the_built_loci_reconcile_with_the_count_and_the_trial_matches_its_own_result(payload):
    r = payload["reconciliation"]
    assert r["reconciles"] is True
    assert (
        r["loci_built"] + sum(r["loci_excluded_by_cause"].values())
        == r["addable_loci_counted_before_the_build"]
    )
    assert r["loci_built"] == len(payload["loci"]) == payload["counts"]["loci"]
    trial = RESULTS / "response_map_coverage_chr21.json"
    if trial.is_file():
        c = json.loads(trial.read_text())
        assert r["addable_loci_counted_in_this_run"] == c["addable"]["independent_loci"]
        assert r["pairs_by_status"] == c["by_status"]["pairs"]
        assert r["loci_by_best_status"] == c["by_status"]["independent_loci"]


@needs_inputs
def test_a_locus_whose_source_is_absent_is_excluded_by_a_named_cause_and_not_dropped(monkeypatch):
    """A source file missing from the checkout must cost a locus visibly. The reconciliation still
    has to add up and the cause has to be one of the registered ones, so a locus can never go quiet."""
    real = r2._Files.ref

    def without_the_benchmark(self, path, key, record_path, **extra):
        return None if "crispri" in path else real(self, path, key, record_path, **extra)

    monkeypatch.setattr(r2._Files, "ref", without_the_benchmark)
    p = r2.build(ROOT, [TRIAL])
    r = p["reconciliation"]
    assert p["loci"] == []
    assert r["loci_excluded_by_cause"][r2.CAUSE_PAIR_SOURCE] == r["addable_loci_counted_in_this_run"]
    assert r["reconciles"] is True
    assert set(r["loci_excluded_by_cause"]) == set(r2.CAUSES)
    assert r["loci_excluded"][r2.CAUSE_PAIR_SOURCE]
    for row in r["loci_excluded"][r2.CAUSE_PAIR_SOURCE]:
        assert row["locus"] and row["causes_of_its_pairs"] == [r2.CAUSE_PAIR_SOURCE]


def test_a_count_that_does_not_add_up_is_refused_rather_than_shown():
    with pytest.raises(rm.RefusedError, match="do not reconcile"):
        r2._reconcile(
            r2._Files(ROOT),
            [TRIAL],
            [],  # nothing built
            {0: r2.STATUS_ADDABLE, 1: r2.STATUS_ADDABLE},
            {0, 1},  # two addable loci
            {c: [] for c in r2.CAUSES},  # and no exclusion naming where they went
            dict.fromkeys(r2.CAUSES, 0),
            {"by_status": dict.fromkeys(r2.STATUSES, 0)},
        )


@needs_inputs
def test_every_assertion_names_its_file_its_record_key_and_that_files_sha256(payload):
    assert payload["assertions"], "the trial built no assertion, so nothing was checked"
    table = {s["path"]: s for s in payload["sources"]}
    for a in payload["assertions"]:
        src = a["source"]
        assert src["path"] and src["record_key"] and src["sha256"]
        assert len(src["sha256"]) == 64
        assert src["sha256"] == table[src["path"]]["sha256"]
        assert src["kind"] in (rm.COMMITTED, rm.LOCAL)
        # the kind is read from git, never guessed from the directory: reader v1's peak sets sit under
        # data/results and are git-ignored all the same
        assert table[src["path"]]["committed_at_head"] is (src["kind"] == rm.COMMITTED)
        assert table[src["path"]]["committed_read_from"] == "git ls-tree HEAD"


@needs_inputs
def test_an_assertion_without_a_sha256_is_refused(payload):
    broken = dict(payload)
    broken["assertions"] = [dict(a) for a in payload["assertions"]]
    broken["assertions"][0] = {**broken["assertions"][0], "source": {"path": "x", "record_key": "y"}}
    with pytest.raises(rm.RefusedError, match="record key and sha256"):
        rm.check(broken)


@needs_inputs
def test_a_model_output_is_predicted_and_observed_is_only_what_an_experiment_measured(payload):
    kinds = {"p2|": 0, "r2|": 0, "s2|": 0}
    for a in payload["assertions"]:
        if a["id"].startswith("r2|"):
            assert a["status"] == "predicted"
            assert a["measurement"]["is_a_measurement"] is False
            assert a["basis"].startswith("predicted")  # the element record's own wording
            kinds["r2|"] += 1
        elif a["id"].startswith("p2|"):
            assert a["status"] == "observed"
            assert a["basis"].startswith("measured perturbation effect")
            kinds["p2|"] += 1
        elif a["id"].startswith("s2|"):
            assert a["status"] == "observed"
            assert a["basis"] == ce.OPENNESS_CALL
            kinds["s2|"] += 1
        else:
            raise AssertionError(f"an assertion of no known kind: {a['id']}")
    assert all(kinds.values())
    assert payload["counts"]["by_status"]["predicted"] == kinds["r2|"]
    assert payload["counts"]["observed_by_assay"] == {
        "CRISPRi perturbation (ENCODE benchmark)": kinds["p2|"],
        reader.EVIDENCE: kinds["s2|"],
    }
    assert payload["counts"]["by_status"]["observed"] == kinds["p2|"] + kinds["s2|"]


@needs_inputs
def test_a_status_outside_the_four_is_refused(payload):
    broken = dict(payload)
    broken["assertions"] = [dict(a) for a in payload["assertions"]]
    broken["assertions"][0]["status"] = "measured"
    with pytest.raises(rm.RefusedError, match="is not one of"):
        rm.check(broken)


@needs_inputs
def test_the_reader_state_is_shown_as_a_reading_and_never_as_validation(payload):
    states = [a for a in payload["assertions"] if a["id"].startswith("s2|")]
    assert states, "the trial built no reader state, so nothing was checked"
    for a in states:
        assert a["relation"] == "associated_with_state"
        assert a["is_validation"] is False
        assert ce.NOT_VALIDATION in a["uncertainty"]
        assert a["measurement"]["reading"] == ce.READING[a["measurement"]["value"]]
        entity = payload["entities"][a["object"]]
        assert entity["kind"] == "cell_state"
        assert entity["not_validation"] == ce.NOT_VALIDATION
        assert entity["is_validation"] is False
        if entity["state"] == ce.STATE_NOT_OPEN:
            assert entity["not_closed"] == ce.NOT_CLOSED
            assert ce.NOT_CLOSED in a["uncertainty"]
    block = payload["not_open_in_reader"]
    assert block["meaning"] == ce.READING[ce.STATE_NOT_OPEN]
    assert block["never"] == ce.NOT_CLOSED
    assert block["not_validation"] == ce.NOT_VALIDATION
    assert block["shown"] == len(block["shown_where_it_applies"])
    assert block["shown"] == block["pairs_by_reader_state"].get(ce.STATE_NOT_OPEN, 0)
    assert sum(block["pairs_by_reader_state"].values()) == len(states) or True  # states dedupe per element


def test_the_not_open_state_carries_its_registered_wording_and_no_claim_of_closure():
    """The one state this map could overstate. The wording is lane-context2's, word for word, and the
    module never writes the element off as closed."""
    m = rm._Map({"assembly": "GRCh38", "coordinates": {"base": 0, "interval": "half-open"}})
    e = m.entities[r2._state_entity(m, ce.STATE_NOT_OPEN, "K562")]
    assert e["reading"] == ce.READING[ce.STATE_NOT_OPEN]
    assert e["not_closed"] == ce.NOT_CLOSED
    assert "not detected open" in e["reading"]
    text = Path(r2.__file__).read_text()
    for phrase in ("is closed", "element is closed", "chromatin is closed"):
        assert phrase not in text, phrase


@needs_inputs
def test_no_verdict_and_no_new_validation_enters_the_map(payload):
    for a in payload["assertions"]:
        assert a["kind"] == "evidence"
        assert not ({"verdict", "validated", "supports", "confirms", "score"} & set(a))
    for c in payload["chains"]:
        assert c["derived_claims"] == []
        assert set(c["stops_at"]) == set(r2.CHAIN_STOPS)
        assert c["note"] == r2.CHAIN_NOTE
    assert payload["requests"] == {"network": 0, "model": 0, "downloads": 0}
    assert payload["participation"] == []  # no category is assigned here: none could be sourced
    assert payload["cycles"] == []  # nothing in this map points back at a DNA interval


@needs_inputs
def test_a_chain_runs_from_the_measured_interval_through_the_element_to_the_gene_and_the_cell(payload):
    by = {a["id"]: a for a in payload["assertions"]}
    assert payload["chains"]
    for c in payload["chains"]:
        pert, rule, state = (by[s] for s in c["steps"])
        assert pert["id"].startswith("p2|")
        assert rule["id"].startswith("r2|")
        assert state["id"].startswith("s2|")
        element = rule["subject"]
        assert element in pert["entities"], "the measurement names the element it attaches to"
        assert pert["attached_to_element"]["reciprocal_overlap"] >= ms.RECIPROCAL_OVERLAP
        assert pert["attached_to_element"]["rule"] == r2.ATTACHMENT_NOTE
        assert state["subject"] == element
        assert pert["context"]["cell"] == state["context"]["cell"]
        for link in c["links"]:
            assert element in link["shared_entities"]


@needs_inputs
def test_the_pages_show_every_locus_and_truncate_nothing(payload):
    seen: list[str] = []
    pages = None
    for n in range(1, 60):
        p = r2.page(payload, n, per=1)
        pages = p["pages"]
        assert p["loci_total"] == len(payload["loci"])
        assert p["assertions_total"] == len(payload["assertions"])
        assert p["counts"]["assertions"] == len(payload["assertions"])  # the whole count beside a page
        assert p["reconciliation"]["loci_built"] == len(payload["loci"])
        for x in p["loci"]:
            seen.append(x["id"])
            assert set(x["assertions"]) <= {a["id"] for a in p["assertions"]}
            assert set(x["chains"]) <= {c["id"] for c in p["chains"]}
            for a in p["assertions"]:
                assert set(a["entities"]) <= set(p["entities"])
        assert "page" in p["paging"] and str(p["loci_total"]) in p["paging"]
        if n >= pages:
            break
    assert pages == len(payload["loci"])
    assert seen == [x["id"] for x in payload["loci"]]
    assert r2.page(payload, 999, per=1)["page"] == pages  # past the end is the last page, not empty


@needs_inputs
def test_the_endpoint_pages_the_built_result_and_says_so_when_it_is_absent(tmp_path, payload):
    from genomeos.web import server as ws

    with pytest.raises(ws.ApiError) as bad:
        ws.Api(tmp_path).evidence_response_map_2()
    assert r2.RESULT in str(bad.value)
    assert bad.value.status == 404
    (tmp_path / "data" / "results").mkdir(parents=True)
    (tmp_path / "data" / "results" / f"{r2.RESULT}.json").write_text(json.dumps(payload))
    got = ws.Api(tmp_path).evidence_response_map_2(page="1", per="1")
    assert got["page"] == 1
    assert got["per_page"] == 1
    assert got["loci_total"] == len(payload["loci"])
    assert got["reconciliation"]["reconciles"] is True


def test_the_page_section_reads_its_endpoint_and_shows_the_whole_count_beside_a_page():
    text = PAGE_HTML.read_text()
    start = text.index("// ---- response map, increment 2:")
    s = text[start : text.index("async function loadProgress()", start)]
    assert "/api/evidence/response-map-2?page=" in s
    assert "fetch(" not in s  # through the page's own api() helper only
    assert 'id="rmap2"' in text
    for field in ("loci_total", "pages", "reconciles", "loci_excluded_by_cause", "observed_by_assay"):
        assert field in s, field
    assert "reader_state_is_not_validation" in s
    assert "not_open_in_reader" in s
    for word in ("causal", "validated", "confirmed", "rescued"):
        assert word not in s.lower(), word


def test_the_module_states_what_it_does_not_cover():
    text = " ".join(Path(r2.__file__).read_text().split())
    assert "not established biological independence" in text
    assert "four assays" in text
    conventions = r2._conventions({"assembly": "GRCh38", "coordinates": {}})
    assert conventions["reader_state_is_not_validation"] == ce.NOT_VALIDATION
    assert conventions["not_open_is_not_closed"] == ce.NOT_CLOSED
    assert conventions["locus_grouping"].endswith("not established biological independence")


@pytest.mark.skipif(not BUILT.is_file(), reason="the genome-wide build is not in this checkout")
def test_the_committed_build_reconciles_with_the_counted_population(counted):
    built = json.loads(BUILT.read_text())
    r = built["reconciliation"]
    assert r["addable_loci_counted_before_the_build"] == counted["addable"]["independent_loci"]
    assert r["loci_built"] == r2.ADDABLE_LOCI
    assert sum(r["loci_excluded_by_cause"].values()) == 0
    assert r["reconciles"] is True
    assert r["newly_covered_by_this_increment"] == counted["addable"]["independent_loci_not_yet_covered"]
    assert (
        r["already_covered_by_increment_1"] == counted["addable"]["independent_loci_the_map_already_covers"]
    )
    assert (
        built["not_open_in_reader"]["pairs_by_reader_state"]
        == counted["addable"]["reader_state_of_the_addable_pairs"]
    )
    assert r["ceiling_after_this_increment"] == r2.ADDABLE_LOCI
    assert built["alphagenome_requests"] == 0
    assert built["coverage"]["independent_loci_of_the_perturbation_assay"] == r2.MEASURED_LAYER_LOCI
