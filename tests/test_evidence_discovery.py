"""The discovery view: set-valued CRISPRi evidence from the committed discovery review, read-only.

Each test reads the committed result independently (json.load, its own id parsing, its own cell
comparison) and holds the served view to it, so a view that drifted from its file fails here rather
than on the page. The existing /api/evidence payload is held to digests recorded before the view was
added (at f4c391e), on a fixed pair of programs, so the addition is shown to change nothing there.
"""

import hashlib
import json
import re
import subprocess
import threading
import urllib.error
import urllib.request
from collections import Counter
from http.server import ThreadingHTTPServer
from pathlib import Path

import pytest

from genomeos import evidence, evidence_discovery
from genomeos.web.server import Api, ApiError, Handler

ROOT = Path(__file__).resolve().parent.parent
SOURCE = "data/results/discovery_review.json"
STRATUM = "Gasperini2019|K562"
CONTEXTS = ("all_attached_claims_match", "some_match", "none_match", "no_assessable_claim_context")
BANNED = ("causal", "identified", "validated", "rescued")

SOME_MATCH = "training|chr12:9764686-9765406|CD69|K562"
LOST_RESOLUTION = "training|chr12:33588022-33588689|SYT10|K562"  # in section (a) and section (b)
NONE_MATCH = "training|chr19:18291477-18292285|JUND|K562"


@pytest.fixture(scope="module")
def data() -> dict:
    return json.loads((ROOT / SOURCE).read_bytes())


@pytest.fixture(scope="module")
def view() -> dict:
    return evidence_discovery.view(ROOT)


def _parse(obs_id: str) -> dict:
    split, interval, gene, cell = obs_id.split("|")
    chrom, span = interval.split(":")
    start, end = span.split("-")
    return {"split": split, "chrom": chrom, "start": int(start), "end": int(end), "gene": gene, "cell": cell}


def _cell(c: str) -> str:
    return re.sub(r"[^a-z0-9]", "", c.lower())


def _claim(claim: dict, assay_cell: str) -> str:
    cells = [c for c in claim["stated_cells"] if _cell(c) not in ("", "unknown")]
    if not cells:
        return "no_stated_cell"
    return "match" if any(_cell(c) == _cell(assay_cell) for c in cells) else "other_cell"


def _row(record: dict) -> str:
    """The row context, recomputed apart from the module."""
    got = [_claim(a, record["cell"]) for a in record["attached"]]
    if not got:
        return "no_attached_claims"
    if set(got) == {"match"}:
        return "all_attached_claims_match"
    if "match" in got:
        return "some_match"
    return "none_match" if "other_cell" in got else "no_assessable_claim_context"


def _dist(records) -> tuple:
    c = Counter(_row(r) for r in records)
    return tuple(c.get(k, 0) for k in CONTEXTS)


def _discovered(data: dict) -> list:
    return [
        r for r in data["changed_observations"] if r["change"] in ("newly_attached", "attached_links_added")
    ]


# --- S1: the source is pinned ------------------------------------------------------------------------


def test_the_source_sha256_equals_the_committed_blob(view):
    blob = subprocess.run(["git", "show", f"HEAD:{SOURCE}"], cwd=ROOT, capture_output=True, check=True).stdout
    assert view["source"]["path"] == SOURCE
    assert view["source"]["sha256"] == hashlib.sha256(blob).hexdigest()


# --- S2 and S3: one record per observation, sections by reference -----------------------------------


def test_every_served_record_equals_the_files_record_on_the_listed_fields(view, data):
    by_id = {r["id"]: r for r in data["changed_observations"]}
    copied = (
        "study", "cell", "split", "outcome", "tested_width_bp", "change", "multiplicity", "union_coverage",
        "candidates", "resolution", "verdict_current_rule", "status_under_new_rule",
        "direction_description_policy", "stated_cell_agreement",
    )  # fmt: skip
    for obs_id, rec in view["records"].items():
        src = by_id[obs_id]
        assert rec["id"] == obs_id
        for field in copied:
            assert rec[field] == src[field], (obs_id, field)
        assert [{k: v for k, v in a.items() if k != "context"} for a in rec["attached"]] == src["attached"]
        iv = _parse(obs_id)
        assert rec["interval"] == {"chrom": iv["chrom"], "start": iv["start"], "end": iv["end"]}
        assert rec["gene"] == iv["gene"] and rec["cell"] == iv["cell"] and rec["split"] == iv["split"]
        assert rec["verdict_rule"] == "C4 any-cell attachment, reciprocal-0.5 rule (b7e4bf0)"


def test_section_a_holds_exactly_the_files_stratum(view, data):
    a = view["sections"]["discovery"]
    stratum = data["screen"]["attached_links_differ"][STRATUM]
    study, cell = STRATUM.split("|")
    expected = {r["id"] for r in _discovered(data) if r["study"] == study and r["cell"] == cell}
    assert set(a["ids"]) == expected and len(a["ids"]) == len(expected)
    assert a["rows"] == {
        "count": stratum["real"],
        "unit": "observations",
        "checked_against_file": stratum["real"],
    }
    assert a["denominator"]["count"] == stratum["observations"] and a["denominator"]["unit"] == "observations"
    assert f"{data['universe']['observations']:,}" in a["denominator"]["definition"]
    assert "attached or not" in a["denominator"]["definition"]
    # links are counted as links, apart from observations
    assert a["attached_links"] == {
        "count": sum(len(view["records"][i]["attached"]) for i in a["ids"]),
        "unit": "links",
    }
    assert a["display"]["rows"].startswith(f"{stratum['real']:,} observations ")
    assert a["display"]["denominator"].startswith(f"denominator: {stratum['observations']:,} observations")
    assert a["display"]["links"].startswith(f"{a['attached_links']['count']:,} attached links ")


def test_section_b_holds_exactly_the_thirteen_newly_ambiguous_observations(view, data):
    b = view["sections"]["resolution_change"]
    each = data["newly_ambiguous"]["each"]
    assert b["ids"] == [e["id"] for e in each] and len(b["ids"]) == 13
    for e in each:
        rec = view["records"][e["id"]]
        assert rec["verdict_current_rule"] == e["historical_verdict_current_rule"]
        assert rec["display"]["verdict"].startswith("historical verdict: ")
        assert rec["display"]["status"] == "unresolved under the new rule"
    # the one observation in both sections is one record, referred to twice
    both = set(b["ids"]) & set(view["sections"]["discovery"]["ids"])
    assert both == {LOST_RESOLUTION} and b["also_in_discovery"] == [LOST_RESOLUTION]
    assert len({r["study"] for r in (view["records"][i] for i in b["ids"])}) > 1  # across studies


def test_record_ids_are_unique_and_every_section_reference_resolves(view):
    records = view["records"]
    a, b = view["sections"]["discovery"]["ids"], view["sections"]["resolution_change"]["ids"]
    assert len(a) == len(set(a)) and len(b) == len(set(b))
    assert all(i in records for i in a + b)
    assert set(records) == set(a) | set(b) and len(records) == len(set(a) | set(b))
    assert all(k == r["id"] for k, r in records.items())
    for group in view["sections"]["other_strata"]["groups"]:
        assert all(
            set(s) >= {"stratum", "n", "observations", "screen"} and "ids" not in s for s in group["strata"]
        )


# --- context categories ---------------------------------------------------------------------------


def test_context_categories_recomputed_match_the_registered_distributions(view, data):
    records = view["records"]
    for rec in records.values():
        assert rec["row_context"] == _row(rec)
        assert [a["context"] for a in rec["attached"]] == [_claim(a, rec["cell"]) for a in rec["attached"]]
    a = [records[i] for i in view["sections"]["discovery"]["ids"]]
    b = [records[i] for i in view["sections"]["resolution_change"]["ids"]]
    assert _dist(a) == (41, 9, 66, 0)
    assert _dist(b)[:3] == (3, 0, 10) and _dist(b)[3] == 0
    discovered = _discovered(data)
    assert len(discovered) == 184 and _dist(discovered) == (54, 20, 110, 0)
    # the module's own function agrees on all 184, not only on the served records
    assert Counter(evidence_discovery.row_context(r) for r in discovered) == Counter(
        _row(r) for r in discovered
    )
    assert view["sections"]["discovery"]["context_counts"]["counts"] == {
        **dict(zip(CONTEXTS, (41, 9, 66, 0), strict=True)),
        "no_attached_claims": 0,
    }
    assert view["sections"]["resolution_change"]["context_counts"]["counts"]["none_match"] == 10


def test_no_attached_claims_is_kept_apart_from_no_assessable_context():
    base = {"cell": "K562"}
    assert evidence_discovery.row_context({**base, "attached": []}) == "no_attached_claims"
    unstated = {**base, "attached": [{"stated_cells": []}, {"stated_cells": ["unknown"]}]}
    assert evidence_discovery.row_context(unstated) == "no_assessable_claim_context"
    mixed = {**base, "attached": [{"stated_cells": ["K562"]}, {"stated_cells": []}]}
    assert evidence_discovery.row_context(mixed) == "some_match"


# --- coverage ---------------------------------------------------------------------------------------


def test_union_coverage_and_candidate_fractions_equal_the_files(view, data):
    by_id = {r["id"]: r for r in data["changed_observations"]}
    sentence = re.compile(
        r"^(\d+) admitted registry candidates; ([\d.]+)% of the tested interval lies inside them; "
        r"the rest is unresolved\.$"
    )
    for obs_id, rec in view["records"].items():
        src = by_id[obs_id]
        assert rec["union_coverage"] == src["union_coverage"]
        assert [(c["id"], c["fraction_covered"]) for c in rec["candidates"]] == [
            (c["id"], c["fraction_covered"]) for c in src["candidates"]
        ]
        m = sentence.match(rec["display"]["coverage"])
        assert m, rec["display"]["coverage"]
        assert int(m[1]) == src["multiplicity"]["min_side_half"] == len(src["candidates"])
        assert abs(float(m[2]) - src["union_coverage"]["min_side_half"] * 100) < 0.005


# --- (c): counts only, registered readings unchanged -------------------------------------------------


def test_other_strata_are_labelled_counts_with_the_registered_readings(view, data):
    c = view["sections"]["other_strata"]
    screen, decision = data["screen"]["attached_links_differ"], data["decision_1_reading"]
    labels = [g["label"] for g in c["groups"]]
    assert labels == [
        "screen passed, sparse (n ≤ 16)",
        "screen failed",
        "screen undefined",
        "nothing discovered",
    ]
    lists = ["strata_passing", "strata_failing", "strata_undefined_shifted_zero", "strata_nothing_discovered"]
    for g, key in zip(c["groups"], lists, strict=True):
        assert [s["stratum"] for s in g["strata"]] == [n for n in decision[key] if n != STRATUM]
        for s in g["strata"]:
            f = screen[s["stratum"]]
            assert (s["n"], s["observations"], s["enrichment"], s["screen"]) == (
                f["real"],
                f["observations"],
                f["enrichment"],
                f["screen"],
            )
            assert s["n_unit"] == "observations whose attached links differ"
    assert all(s["n"] <= 16 for s in c["groups"][0]["strata"])
    failed = c["groups"][1]["strata"]
    assert [s["stratum"] for s in failed] == ["WTC11_DC_TAP|WTC11"] and round(
        failed[0]["enrichment"], 2
    ) == 0.86
    assert [s["stratum"] for s in c["groups"][2]["strata"]] == ["Nasser2021|GM12878"]
    assert (c["reading"], c["meaning"]) == (decision["reading"], decision["meaning"])
    assert c["note"] == "passing the screen is not sufficient evidence for implementation"
    assert c["pooled"]["screen"] == screen["all_observations"]["screen"]


# --- S4: fixed wording and the rendered rows ---------------------------------------------------------


def test_fixed_wording_is_verbatim(view):
    assert view["header"] == "The discovery view adds access to evidence, not new validation."
    assert view["notes"]["context"] == (
        "The stated-cell subset is selected using model predictions and is not representative of "
        "all screened pairs."
    )
    note = (
        "Historical C4 reading (any-cell attachment, reciprocal-0.5 rule, b7e4bf0); "
        "not a v3 stated-context verdict."
    )
    assert view["notes"]["verdict"] == note
    assert all(
        r["verdict_note"] == note and r["display"]["verdict_note"] == note for r in view["records"].values()
    )


def test_the_rule_version_names_a_commit_in_this_repository():
    assert subprocess.run(["git", "cat-file", "-e", "b7e4bf0^{commit}"], cwd=ROOT).returncode == 0


def test_a_mixed_context_row_renders_exactly(view):
    assert view["records"][SOME_MATCH]["display"] == {
        "observation": "tested interval chr12:9764686-9765406 (720 bp) · measured gene CD69 · "
        "assay cell K562 · Gasperini2019 · training split · outcome: significant decrease",
        "coverage": "2 admitted registry candidates; 76.81% of the tested interval lies inside them; "
        "the rest is unresolved.",
        "candidates": "admitted registry candidates: EH38E1591980 (covered fraction 1.0); "
        "EH38E3002772 (covered fraction 1.0)",
        "links": "2 attached links under the new rule",
        "claims": [
            "EH38E1591980 activates CD69 · stated cell: K562 · context: match",
            "EH38E3002772 activates CD69 · stated cell: GM12891 · context: other_cell",
        ],
        "context": "row context: some attached claims match the assay cell (some_match)",
        "verdict": "historical verdict: unattached [C4 any-cell attachment, reciprocal-0.5 rule (b7e4bf0)]",
        "verdict_note": "Historical C4 reading (any-cell attachment, reciprocal-0.5 rule, b7e4bf0); "
        "not a v3 stated-context verdict.",
        "status": "unresolved under the new rule",
        "direction": "direction relative to the attached claims (a description only): agrees",
    }


def test_a_lost_resolution_row_in_both_sections_renders_exactly(view):
    assert view["records"][LOST_RESOLUTION]["display"] == {
        "observation": "tested interval chr12:33588022-33588689 (667 bp) · measured gene SYT10 · "
        "assay cell K562 · Gasperini2019 · training split · outcome: significant decrease",
        "coverage": "2 admitted registry candidates; 75.41% of the tested interval lies inside them; "
        "the rest is unresolved.",
        "candidates": "admitted registry candidates: EH38E1604341 (covered fraction 1.0); "
        "EH38E3011039 (covered fraction 1.0, admitted under the current rule too)",
        "links": "2 attached links under the new rule",
        "claims": [
            "EH38E3011039 activates SYT10 · stated cell: K562 · context: match · "
            "attached under the current rule too",
            "EH38E1604341 activates SYT10 · stated cell: K562 · context: match",
        ],
        "context": "row context: all attached claims match the assay cell (all_attached_claims_match)",
        "verdict": "historical verdict: supported [C4 any-cell attachment, reciprocal-0.5 rule (b7e4bf0)]",
        "verdict_note": "Historical C4 reading (any-cell attachment, reciprocal-0.5 rule, b7e4bf0); "
        "not a v3 stated-context verdict.",
        "status": "unresolved under the new rule",
        "direction": "direction relative to the attached claims (a description only): agrees",
    }


def test_a_none_match_row_renders_exactly(view):
    assert view["records"][NONE_MATCH]["display"] == {
        "observation": "tested interval chr19:18291477-18292285 (808 bp) · measured gene JUND · "
        "assay cell K562 · Gasperini2019 · training split · outcome: significant decrease",
        "coverage": "2 admitted registry candidates; 76.98% of the tested interval lies inside them; "
        "the rest is unresolved.",
        "candidates": "admitted registry candidates: EH38E1944166 (covered fraction 0.8841); "
        "EH38E3295927 (covered fraction 1.0)",
        "links": "2 attached links under the new rule",
        "claims": [
            "EH38E3295927 activates JUND · stated cell: GM12878 · context: other_cell",
            "EH38E1944166 activates JUND · stated cell: HeLa_S3 · context: other_cell",
        ],
        "context": "row context: no attached claim matches the assay cell (none_match)",
        "verdict": "historical verdict: unattached [C4 any-cell attachment, reciprocal-0.5 rule (b7e4bf0)]",
        "verdict_note": "Historical C4 reading (any-cell attachment, reciprocal-0.5 rule, b7e4bf0); "
        "not a v3 stated-context verdict.",
        "status": "unresolved under the new rule",
        "direction": "direction relative to the attached claims (a description only): agrees",
    }


def _page_section() -> str:
    page = (ROOT / "genomeos" / "web" / "static" / "index.html").read_text()
    html = page[page.index('<div class="card" id="evd"') : page.index("</section>", page.index('id="evd"'))]
    js = page[page.index("// ---- evidence explorer: set-valued") : page.index("async function loadProgress")]
    return html + js


def test_the_page_shows_the_served_strings_and_no_banned_word(view):
    section = _page_section()
    assert "/api/evidence/discovery" in section
    used = "d.coverage d.candidates d.claims d.context d.verdict d.verdict_note d.status d.direction"
    for name in used.split() + ["v.header", "v.notes.context", "v.source.sha256", "c.note"]:
        assert name in section, name
    served = json.dumps(view, ensure_ascii=False).lower()
    for word in BANNED:
        assert word not in served, word
        assert word not in section.lower(), word
        assert word not in (ROOT / "genomeos" / "evidence_discovery.py").read_text().lower(), word


def test_neither_count_is_written_into_the_code():
    source = (ROOT / "genomeos" / "evidence_discovery.py").read_text()
    for figure in ("116", "5299", "5,299", "14734", "14,734"):
        assert figure not in source, figure
    assert "116" not in _page_section() and "5,299" not in _page_section()


# --- refusal, the endpoint, and the existing outputs ------------------------------------------------


def test_the_view_refuses_when_section_a_disagrees_with_the_file(tmp_path, data):
    bad = json.loads(json.dumps(data))
    bad["screen"]["attached_links_differ"][STRATUM]["real"] += 1
    (tmp_path / "data" / "results").mkdir(parents=True)
    (tmp_path / SOURCE).write_text(json.dumps(bad))
    with pytest.raises(evidence_discovery.RefusedError, match="section \\(a\\)"):
        evidence_discovery.view(tmp_path)
    with pytest.raises(ApiError) as e:
        Api(tmp_path).evidence_discovery()
    assert e.value.status == 500


def test_a_checkout_without_the_result_gets_a_404(tmp_path):
    with pytest.raises(ApiError) as e:
        Api(tmp_path).evidence_discovery()
    assert e.value.status == 404


def _serve(root: Path):
    saved = getattr(Handler, "api", None)
    Handler.api = Api(root)
    srv = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    srv.verbose = False  # type: ignore[attr-defined]
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    return srv, saved


def _get(srv, path: str) -> tuple[int, bytes]:
    try:
        with urllib.request.urlopen(f"http://127.0.0.1:{srv.server_address[1]}{path}", timeout=60) as r:
            return r.status, r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.read()


def test_the_endpoint_serves_the_view(view):
    srv, saved = _serve(ROOT)
    try:
        status, body = _get(srv, "/api/evidence/discovery")
    finally:
        srv.shutdown()
        Handler.api = saved
    assert status == 200 and json.loads(body) == view


BASE = """module test.base

param shared_rate = 0.5 h {
  evidence: curated "BioModels BIOMD0000000012"
  confidence: 0.8
}
"""

TOP = """module test.top

import base.bio

gene SOX2 {
  symbol: SOX2
  evidence: experimental "GENCODE v50"
  confidence: 0.95
}

protein NANOG { evidence: experimental "GENCODE v50"; confidence: 0.9 }

param guessed = 0.35 {
  evidence: inferred "order of magnitude"
  confidence: 0.3
}

rule SOX2 activates NANOG {
  strength: 1.0
  threshold: 1.0
  evidence: predicted "AlphaGenome"
  confidence: 0.4
}
"""

# sha256 of each payload, recorded at f4c391e before the discovery view existed. The payload's key order
# is not fixed between processes, so the digest is of its canonical form (keys sorted), which fixes every
# key, value and list order.
EVIDENCE_BEFORE = {
    "/api/evidence": "12d2bb29e85d3498b69e1498ba5c2b30b2e806724b5ca4207dc04abd6ec38abe",
    "/api/evidence?kinds=predicted,inferred&max_confidence=0.5&limit=10": (
        "8c3b5e374379bf8d69fdac76c6d9dc664c06131552fc48d6580a74cf46de221c"
    ),
    "/api/evidence?csv=1": "ef89dc25e92d597110dda74c619bcaae92420e027aeff736ec28aff533caf528",
}
# 2026-09-29: `complete: true`, `parse_errors: []` and `parse_error_count: 0` added intentionally to
# every clean /api/evidence payload (ROADMAP item 9: a read says whether every program parsed, and a
# missing field is not a synonym for complete). The digests above are kept as recorded at f4c391e;
# these are of the same three payloads, which differ from them only by those three added keys.
EVIDENCE_BEFORE = {
    "/api/evidence": "8b3c00761010e429a9f418eaffe94cb3f7b0a86c49b2d9baa1d8d9b6a7db645c",
    "/api/evidence?kinds=predicted,inferred&max_confidence=0.5&limit=10": (
        "4108fb4ba1b48292eda5da96ce934cc68e45d2f629b5bf798164156700b43f69"
    ),
    "/api/evidence?csv=1": "abef5feeaa2627b6bd386376a6fc830d6d05ec55e24c26e3ba815d0967e756be",
}


def _canonical(payload: dict) -> str:
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


@pytest.fixture
def two_programs(tmp_path, monkeypatch):
    d = tmp_path / "data" / "demo"
    d.mkdir(parents=True)
    (d / "base.bio").write_text(BASE)
    (d / "top.bio").write_text(TOP)
    monkeypatch.setattr(evidence, "PROGRAM_DIRS", (Path("data/demo"),))
    evidence._cached.cache_clear()
    yield tmp_path
    evidence._cached.cache_clear()


def test_the_existing_evidence_payload_is_unchanged(two_programs):
    api = Api(two_programs)
    assert _canonical(api.evidence()) == EVIDENCE_BEFORE["/api/evidence"]
    filtered = api.evidence(kinds="predicted,inferred", max_confidence="0.5", limit=10)
    assert (
        _canonical(filtered)
        == EVIDENCE_BEFORE["/api/evidence?kinds=predicted,inferred&max_confidence=0.5&limit=10"]
    )
    assert _canonical(api.evidence(csv=True)) == EVIDENCE_BEFORE["/api/evidence?csv=1"]
    srv, saved = _serve(two_programs)
    try:
        for path, digest in EVIDENCE_BEFORE.items():
            status, body = _get(srv, path)
            assert status == 200 and _canonical(json.loads(body)) == digest, path
    finally:
        srv.shutdown()
        Handler.api = saved
