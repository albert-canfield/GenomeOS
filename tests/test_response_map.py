"""The cellular control and response map, increment 1: the beta-globin locus in K562, read-only.

Each test reads the committed sources itself (json.load, gzip, line numbers) and holds the map to them,
so an assertion that drifted from its record fails here. The rules the map enforces are shown twice:
on the real example, and on a small fixture that breaks one rule and must be refused.
"""

import copy
import gzip
import hashlib
import importlib.util
import json
import socket
from pathlib import Path

import pytest

from genomeos import evidence_discovery, response_map
from genomeos.response_map import LOCAL, UNAVAILABLE, RefusedError

ROOT = Path(__file__).resolve().parent.parent
COMMITTED = (
    "data/results/crispri_direction.json",
    "data/results/discovery_review.json",
    "data/results/attribution_correctness_v3.json",
    "data/results/loci_benchmark.json",
    "genomeos/lib/data/proteome.json.gz",
    "data/organisms/human/erythrocyte.bio",
)
REILLY = ("cd|chr11:5253147-5253547|HBG1|Reilly", "cd|chr11:5275847-5276247|HBG1|Reilly")
CLAIM = "v3|EH38E2941908|target|HBE1|HT1080"
# the canonical sha256 of /api/evidence/discovery, recorded at 304a1c8 before the response map existed
DISCOVERY_BEFORE = "5bb466a170fea880e6feaeaca42bab26f540293b2bdf827d5c4ac2b03e061fe3"
HAS_CACHE = (ROOT / response_map.HIC).is_file() and all(
    (ROOT / p).is_file() for p in response_map.CRISPRI_FILES.values()
)


@pytest.fixture(scope="module")
def payload() -> dict:
    return response_map.build(ROOT)


@pytest.fixture(scope="module")
def by_id(payload) -> dict:
    return {a["id"]: a for a in payload["assertions"]}


@pytest.fixture
def bare_root(tmp_path) -> Path:
    """A checkout holding the committed sources only: no local cache."""
    for rel in COMMITTED:
        (tmp_path / rel).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / rel).symlink_to(ROOT / rel)
    return tmp_path


def _load(rel: str):
    raw = (ROOT / rel).read_bytes()
    if rel.endswith(".bio"):
        return raw.decode().splitlines()
    return json.loads(gzip.decompress(raw) if rel.endswith(".gz") else raw)


def _resolve(doc, path: list):
    if path[0] == "line":
        return {"line": doc[path[1] - 1]}
    for step in path:
        doc = doc[step]
    return doc


def _subset(quoted, record) -> bool:
    if isinstance(quoted, dict):
        return isinstance(record, dict) and all(
            k in record and _subset(v, record[k]) for k, v in quoted.items()
        )
    return quoted == record


# --- 1. real-example reconciliation ----------------------------------------------------------------
def test_every_committed_assertion_survives_retrieval_unchanged(payload):
    docs = {rel: _load(rel) for rel in COMMITTED}
    digests = {rel: hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() for rel in COMMITTED}
    committed = [a for a in payload["assertions"] if a["source"]["kind"] == "committed"]
    assert len(committed) >= 90
    for a in committed:
        src = a["source"]
        assert src["sha256"] == digests[src["path"]], a["id"]
        record = _resolve(docs[src["path"]], src["record_path"])
        assert _subset(a["quoted"], record), a["id"]
    listed = {s["path"]: s["sha256"] for s in payload["sources"] if s["kind"] == "committed"}
    assert listed == digests


def test_the_example_reconciles_with_the_brief_where_the_sources_say_so(payload, by_id):
    rows = _load("data/results/crispri_direction.json")["what_it_named_instead"]["rows"]
    for aid, value in zip(REILLY, (-0.7356, -0.6908), strict=True):
        a = by_id[aid]
        rec = rows[a["source"]["record_path"][2]]
        assert a["measurement"]["value"] == rec["measured"] == value
        assert a["context"]["cell"] == "K562" and a["context"]["condition"] == "not recorded"
        assert a["quoted"] == rec and a["source"]["record_key"].endswith("gene=HBG1,dataset=Reilly]")
    # the committed record names HBG2 for the proximal interval and HBE1 for the distal one
    assert by_id[REILLY[0]]["conflicts"][0]["model_named"] == ["HBG2"]
    assert by_id[REILLY[1]]["conflicts"][0]["model_named"] == ["HBE1"]
    assert {i for i in by_id if i.startswith("cd|") and i.endswith("K562_DC_TAP")} == {
        "cd|chr11:5284543-5285042|HBD|K562_DC_TAP",
        "cd|chr11:5288040-5288539|HBD|K562_DC_TAP",
    }
    sch = by_id["dr|training|chr11:5277927-5278054|HBE1|K562"]
    assert sch["context"]["perturbation"].endswith("(Schraivogel2020)")
    assert sch["measurement"]["outcome"] == "significant_decrease"
    both = {by_id[f"dr|training|chr11:5280370-5281170|{g}|K562"]["subject"] for g in ("HBE1", "HBG2")}
    assert both == {"interval:chr11:5280370-5281170"}
    claim = by_id[CLAIM]
    assert claim["claim"]["stated_cell"] == "HT1080"
    assert claim["verdict"]["reason"] == "not_assessed_in_this_context"
    assert claim["entities"] == ["gene:HBE1", "interval:chr11:5272822-5273045"]
    assert [(x["cell"], x["finding"]) for x in claim["other_cell_evidence"]] == [("K562", "agrees")]
    for a in payload["assertions"]:
        e = payload["entities"].get(a["subject"])
        if e and e["kind"] == "dna_interval":
            assert e["assembly"] == "GRCh38" and e["coordinates"] == {"base": 0, "interval": "half-open"}
    assert payload["why"]["examined_before"] and "304a1c8" in payload["why"]["plan"]


def test_the_committed_complex_is_adult_haemoglobin_and_is_not_extended(payload, by_id):
    lines = _load("data/organisms/human/erythrocyte.bio")
    binds = [a for a in payload["assertions"] if a["id"].startswith("bio|")]
    assert sorted((a["subject"], a["object"]) for a in binds) == [
        ("protein:HBA", "complex:HbA"),
        ("protein:HBB", "complex:HbA"),
    ]
    assert all(lines[a["source"]["record_path"][1] - 1] == a["quoted"]["line"] for a in binds)
    for g in ("HBG1", "HBG2", "HBE1"):
        assert not any(
            a["subject"] == f"protein:{g}" and a["object"] == "complex:HbA" for a in payload["assertions"]
        )
        assert by_id[f"prot|{g}|partners"]["status"] == "predicted"
        assert by_id[f"prot|{g}|function"]["status"] == "inferred"
        assert "not measured in K562" in by_id[f"prot|{g}|function"]["basis"]


# --- 2. contact never becomes regulation -----------------------------------------------------------
def _fixture() -> dict:
    """A small example: two intervals, a gene and its protein, each assertion from a named record."""
    src = {
        "path": "fixture.json",
        "record_key": "k",
        "record_path": ["k"],
        "sha256": "0" * 64,
        "kind": "committed",
    }

    def a(aid, relation, subject, obj, **kw):
        base = {
            "id": aid,
            "kind": "evidence",
            "relation": relation,
            "relation_candidates": [],
            "subject": subject,
            "object": obj,
            "entities": sorted({subject, obj}),
            "status": "observed",
            "measurement": None,
            "conflicts": [],
            "source": {**src, "record_key": aid},
        }
        return {**base, **kw}

    entities = {
        "interval:chr1:0-10": {"id": "interval:chr1:0-10", "kind": "dna_interval"},
        "interval:chr1:90-91": {"id": "interval:chr1:90-91", "kind": "dna_interval"},
        "gene:G": {"id": "gene:G", "kind": "gene"},
        "protein:G": {"id": "protein:G", "kind": "protein"},
    }
    assertions = [
        a("contact", "physically_contacts", "interval:chr1:0-10", "interval:chr1:90-91"),
        a(
            "effect",
            "unresolved",
            "interval:chr1:0-10",
            "gene:G",
            relation_candidates=["changes_measured_rna", "changes_measured_protein"],
            measurement={"direction": "decrease", "significant": True},
        ),
        a("encodes", "encodes", "gene:G", "protein:G", status="inferred"),
    ]
    return {"entities": entities, "assertions": assertions, "participation": [], "chains": []}


def _assemble(f: dict) -> dict:
    return response_map.assemble(f["entities"], f["assertions"], f["participation"], f["chains"])


def test_contact_is_never_a_regulatory_assertion(payload):
    contacts = [a for a in payload["assertions"] if a["relation"] == "physically_contacts"]
    assert contacts
    cited = {aid for p in payload["participation"] for aid in p["from_assertions"]}
    for a in contacts:
        assert not a["relation_candidates"] and "mechanism" not in a and a["id"] not in cited
        assert not (a["measurement"] or {}).get("direction")
        assert any("not regulation" in u for u in a["uncertainty"])
    _assemble(_fixture())  # the clean fixture passes
    tagged = _fixture()
    tagged["participation"] = [
        {
            "id": "P",
            "entity": "interval:chr1:0-10",
            "process": "p",
            "categories": ["regulate_gene_responses"],
            "status": "inferred",
            "from_assertions": ["contact"],
        }
    ]
    with pytest.raises(RefusedError, match="contact never becomes a role"):
        _assemble(tagged)
    directed = _fixture()
    directed["assertions"][0]["measurement"] = {"direction": "decrease"}
    with pytest.raises(RefusedError, match="contact carries no direction"):
        _assemble(directed)


# --- 3. other-cell evidence ------------------------------------------------------------------------
def test_other_cell_evidence_is_shown_and_never_support_in_the_claimed_cell(payload, by_id):
    claims = [a for a in payload["assertions"] if a["kind"] == "claim"]
    assert len(claims) == 8 and CLAIM in by_id
    for c in claims:
        assert c["support_in_stated_cell"] == []
        assert c["verdict"]["label"] == "v3 stated-context verdict" and c["verdict"]["authoritative"]
        assert (
            c["verdict"]["verdict"] == c["quoted"]["verdict"]
            and c["verdict"]["reason"] == c["quoted"]["reason"]
        )
        for x in c["other_cell_evidence"]:
            assert x["cell"] != c["claim"]["stated_cell"] and x["counts_as_support_in_stated_cell"] is False
    # historical C4 readings are labelled apart, and a measurement carries no v3 verdict
    for a in payload["assertions"]:
        if "historical_c4" in a:
            assert a["historical_c4"]["label"] == evidence_discovery.VERDICT_NOTE and "verdict" not in a
    one = _fixture()
    claim = {
        "id": "claim",
        "kind": "claim",
        "relation": None,
        "subject": "interval:chr1:0-10",
        "object": "gene:G",
        "entities": ["gene:G", "interval:chr1:0-10"],
        "status": "predicted",
        "conflicts": [],
        "claim": {"stated_cell": "HT1080"},
        "verdict": {
            "reason": "not_assessed_in_this_context",
            "label": "v3 stated-context verdict",
            "authoritative": True,
        },
        "support_in_stated_cell": [],
        "other_cell_evidence": [
            {"cell": "K562", "finding": "agrees", "counts_as_support_in_stated_cell": False}
        ],
        "source": {"path": "fixture.json", "record_key": "claim", "sha256": "0" * 64, "kind": "committed"},
    }
    one["assertions"].append(claim)
    _assemble(copy.deepcopy(one))
    counted = copy.deepcopy(one)
    counted["assertions"][-1]["other_cell_evidence"][0]["counts_as_support_in_stated_cell"] = True
    with pytest.raises(RefusedError, match="never counts as support"):
        _assemble(counted)
    moved = copy.deepcopy(one)
    moved["assertions"][-1]["support_in_stated_cell"] = moved["assertions"][-1]["other_cell_evidence"]
    moved["assertions"][-1]["other_cell_evidence"] = []
    with pytest.raises(RefusedError, match="from the stated cell only"):
        _assemble(moved)


# --- 4. ambiguous candidates -----------------------------------------------------------------------
def test_set_valued_candidates_stay_set_valued_without_duplicating_the_observation(tmp_path):
    obs = {
        "id": "training|chr11:5000000-5001000|HBG1|K562",
        "study": "S",
        "cell": "K562",
        "split": "training",
        "outcome": "significant_decrease",
        "tested_width_bp": 1000,
        "change": "candidates_added_unattached_under_both",
        "multiplicity": {"reciprocal_half": 0, "min_side_half": 3},
        "union_coverage": {"reciprocal_half": 0.0, "min_side_half": 0.9},
        "candidates": [{"id": f"EH{i}", "fraction_covered": 1.0, "current": False} for i in range(3)],
        "attached": [],
        "resolution": {"reciprocal_half": "unattached", "min_side_half": "unattached"},
        "verdict_current_rule": "unattached",
        "status_under_new_rule": "unattached",
        "direction_description_policy": "not_attached",
        "stated_cell_agreement": {"reciprocal_half": "not_attached", "min_side_half": "not_attached"},
    }
    path = tmp_path / response_map.DISCOVERY
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps({"changed_observations": [obs]}))
    src = response_map._Sources(tmp_path)
    m = response_map._Map({"assembly": "GRCh38", "coordinates": {"base": 0, "interval": "half-open"}})
    response_map._discovery(src, m)
    response_map._participation_of_effects(m)
    assert len(m.assertions) == 1 and len(m.participation) == 1
    cands = m.assertions[0]["unresolved_candidates"]
    assert cands["set"] == ["EH0", "EH1", "EH2"] and "not an established cause" in cands["note"]
    assert m.participation[0]["alternatives"] and m.participation[0]["status"] == "inferred"
    assert [e for e in m.entities if e.startswith("interval:")] == [m.assertions[0]["subject"]]
    twin = copy.deepcopy(m.assertions[0])
    twin["id"] = "the same observation again"
    dup = {"entities": m.entities, "assertions": [*m.assertions, twin], "participation": [], "chains": []}
    with pytest.raises(RefusedError, match="is carried by"):
        _assemble(dup)


def test_the_real_candidate_sets_are_kept_whole(payload, by_id):
    data = _load("data/results/discovery_review.json")["changed_observations"]
    for a in payload["assertions"]:
        if a["id"].startswith("dr|"):
            rec = data[a["source"]["record_path"][1]]
            assert a["unresolved_candidates"]["set"] == [c["id"] for c in rec["candidates"]]
    rows = by_id[CLAIM]["other_cell_evidence"][0]["candidate_rows"]
    if HAS_CACHE:
        assert rows["availability"] == LOCAL and len(rows["set"]) == 2
    else:
        assert rows == {"availability": UNAVAILABLE, "set": None}


# --- 5. missing measurements and underpowered observations -----------------------------------------
def test_missing_context_stays_missing_and_underpowered_is_not_negative(payload):
    perturb = [a for a in payload["assertions"] if a.get("observation_key")]
    under = [a for a in perturb if a["measurement"].get("outcome") == "not_significant_underpowered"]
    assert len(under) == 45
    for a in under:
        assert a["status"] == "unknown" and a["negative_evidence"] is False
        assert (
            a["measurement"]["direction"] is None and "not negative evidence" in a["measurement"]["reading"]
        )
    for a in perturb:
        assert a["time_window"] == "not recorded" and a["context"]["condition"] == "not recorded"
        assert a["relation"] == "unresolved"
        assert a["relation_candidates"] == ["changes_measured_rna", "changes_measured_protein"]
    bad = _fixture()
    bad["assertions"][1]["measurement"] = {"outcome": "not_significant_underpowered", "direction": None}
    with pytest.raises(RefusedError, match="not negative evidence"):
        _assemble(bad)


def test_without_the_local_cache_nothing_is_zero_or_negative(bare_root, payload):
    p = response_map.build(bare_root)
    contacts = [a for a in p["assertions"] if a["relation"] == "physically_contacts"]
    assert contacts and all(
        a["status"] == "unknown" and a["measurement"] is None and a["availability"] == UNAVAILABLE
        for a in contacts
    )
    claims = [a for a in p["assertions"] if a["kind"] == "claim"]
    assert all(a["status"] == "unknown" and a["status_from"] == UNAVAILABLE for a in claims)
    stop = next(s for s in p["stop_points"] if s["id"] == "local_rows_not_assembled")
    assert stop["count"] is None and stop["availability"] == UNAVAILABLE
    assert p["dynamics"]["half_lives"]["availability"] == UNAVAILABLE
    assert all(s["kind"] == "committed" or s["available"] is False for s in p["sources"])
    reactome = [a["process"]["name"] for a in p["assertions"] if "|reactome|" in a["id"]]
    assert reactome and set(reactome) == {UNAVAILABLE}
    # the committed assertions are the same with or without the cache
    keep = {a["id"]: a["quoted"] for a in payload["assertions"] if a["source"]["kind"] == "committed"}
    assert {a["id"]: a["quoted"] for a in p["assertions"] if a["source"]["kind"] == "committed"} == keep


# --- 6. contradictions stay visible ---------------------------------------------------------------
def test_the_model_named_instead_conflict_stays_beside_the_measurement(payload):
    stop = next(s for s in payload["stop_points"] if s["id"] == "model_named_another_gene")
    named = [
        a for a in payload["assertions"] if any(c["kind"] == "model_named_instead" for c in a["conflicts"])
    ]
    assert sorted(stop["assertions"]) == sorted(a["id"] for a in named) and len(named) == 4
    for a in named:
        c = a["conflicts"][0]
        assert a["status"] == "observed" and a["object"] == f"gene:{a['quoted']['gene']}"
        assert c["model_named"] == a["quoted"]["model_named_instead"] and c["resolution"] == "unresolved"
    shared = next(s for s in payload["stop_points"] if s["id"] == "one_interval_several_genes")
    assert "chr11:5280370-5281170 (decrease): HBE1, HBG2" in shared["detail"]


# --- 7. multi-step paths ---------------------------------------------------------------------------
def test_a_chain_lists_assertions_and_manufactures_no_claim(payload, by_id):
    assert payload["chains"]
    for c in payload["chains"]:
        assert c["derived_claims"] == [] and all(s in by_id for s in c["steps"] + c["beside"])
        assert all(link["shared_entities"] for link in c["links"])
        first, last = by_id[c["steps"][0]], by_id[c["steps"][-1]]
        # no assertion joins the perturbed interval to the protein or its process
        assert not any(
            a["subject"] == first["subject"] and (a.get("object") or "").startswith("protein:")
            for a in payload["assertions"]
        )
        assert last["relation"] == "participates_in_process" and last["status"] == "inferred"
    assert not any(a.get("derived_from") for a in payload["assertions"])
    f = _fixture()
    f["chains"] = [
        {"id": "c", "title": "t", "steps": ["effect", "encodes"], "beside": [], "derived_claims": []}
    ]
    before = copy.deepcopy(f["assertions"])
    out = _assemble(f)
    assert out["assertions"] == before and out["chains"][0]["links"][0]["shared_entities"] == ["gene:G"]
    concluded = copy.deepcopy(f)
    concluded["chains"][0]["conclusion"] = "the interval regulates protein G"
    with pytest.raises(RefusedError, match="derives no claim"):
        _assemble(concluded)
    derived = _fixture()
    derived["assertions"].append(
        {**derived["assertions"][2], "id": "made", "derived_from": ["effect", "encodes"]}
    )
    with pytest.raises(RefusedError, match="never derived"):
        _assemble(derived)


# --- 8. a cycle without an execution order ----------------------------------------------------------
def test_a_feedback_set_is_represented_without_an_execution_order(payload):
    assert payload["cycles"] == []  # the example's assertions form none
    f = _fixture()
    back = copy.deepcopy(f["assertions"][2])
    back.update(
        id="feedback",
        relation="binds",
        subject="protein:G",
        object="interval:chr1:0-10",
        entities=["interval:chr1:0-10", "protein:G"],
        source={**back["source"], "record_key": "feedback"},
    )
    f["assertions"].append(back)
    out = _assemble(f)
    assert len(out["assertions"]) == 4
    (cy,) = out["cycles"]
    assert cy["execution_order"] is None
    assert cy["assertions"] == ["effect", "encodes", "feedback"]
    assert cy["entities"] == ["gene:G", "interval:chr1:0-10", "protein:G"]
    ordered = copy.deepcopy(out)
    ordered["cycles"][0]["execution_order"] = ["effect", "encodes", "feedback"]
    with pytest.raises(RefusedError, match="no execution order"):
        response_map.check(ordered)


# --- 9. categories, writes and requests ------------------------------------------------------------
def test_categories_live_in_participation_only_and_nothing_is_written(bare_root, monkeypatch):
    def refuse(*a, **k):
        raise AssertionError("the map made a network connection")

    monkeypatch.setattr(socket, "create_connection", refuse)
    monkeypatch.setattr(socket.socket, "connect", refuse)
    before = sorted(str(p) for p in bare_root.rglob("*"))
    p = response_map.build(bare_root)
    assert sorted(str(x) for x in bare_root.rglob("*")) == before
    assert p["requests"] == {"network": 0, "model": 0, "downloads": 0}
    for e in p["entities"].values():
        assert not {"label", "labels", "role", "roles", "categories"} & set(e)
    assert {c for x in p["participation"] for c in x["categories"]} <= set(response_map.CATEGORIES)
    labelled = _fixture()
    labelled["entities"]["gene:G"]["label"] = "regulate_gene_responses"
    with pytest.raises(RefusedError, match="never on an entity"):
        _assemble(labelled)


def test_the_cli_prints_the_payload_and_a_summary(capsys):
    assert response_map.main(["globin_k562"]) == 0
    out = json.loads(capsys.readouterr().out)
    assert out["example"] == "globin_k562" and out["counts"]["assertions"] == len(out["assertions"])
    assert response_map.main(["globin_k562", "--summary"]) == 0
    text = capsys.readouterr().out
    assert "not assessable with current data" in text and "where the evidence stops" in text


# --- 10. existing interfaces stay compatible --------------------------------------------------------
def _discovery_tests():
    spec = importlib.util.spec_from_file_location(
        "_discovery_tests", ROOT / "tests" / "test_evidence_discovery.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_existing_evidence_payloads_are_unchanged(tmp_path, monkeypatch):
    from genomeos import evidence
    from genomeos.web.server import Api, Handler

    dt = _discovery_tests()
    response_map.view(ROOT)
    assert dt._canonical(Api(ROOT).evidence_discovery()) == DISCOVERY_BEFORE
    d = tmp_path / "data" / "demo"
    d.mkdir(parents=True)
    (d / "base.bio").write_text(dt.BASE)
    (d / "top.bio").write_text(dt.TOP)
    monkeypatch.setattr(evidence, "PROGRAM_DIRS", (Path("data/demo"),))
    evidence._cached.cache_clear()
    try:
        srv, saved = dt._serve(tmp_path)
        try:
            for path, digest in dt.EVIDENCE_BEFORE.items():
                status, body = dt._get(srv, path)
                assert status == 200 and dt._canonical(json.loads(body)) == digest, path
        finally:
            srv.shutdown()
            Handler.api = saved
    finally:
        evidence._cached.cache_clear()


def test_the_endpoint_serves_the_map_and_refuses_what_it_cannot_serve(tmp_path):
    from genomeos.web.server import Api, ApiError, Handler

    dt = _discovery_tests()
    srv, saved = dt._serve(ROOT)
    try:
        status, body = dt._get(srv, "/api/evidence/response-map?example=globin_k562")
        unknown, _ = dt._get(srv, "/api/evidence/response-map?example=nothing")
    finally:
        srv.shutdown()
        Handler.api = saved
    assert status == 200 and json.loads(body) == json.loads(json.dumps(response_map.view(ROOT)))
    assert unknown == 400
    with pytest.raises(ApiError) as e:
        Api(tmp_path).evidence_response_map()
    assert e.value.status == 404
