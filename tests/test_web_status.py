"""The second review's three status measures and the discontinued list (roadmap section 5 item 12).

What these tests hold the Progress tab to: each measure is computed from files at request time;
a biological claim is never shown without the qualification its own files attach to it; a
category with nothing validated reads "none" instead of vanishing; and an investigation closed on
a negative is never counted as a capability. None of them reads a git-ignored cache: the CI status
is written into a temporary checkout by the test itself.
"""

import importlib.util
import json
import re
import shutil
from datetime import UTC, datetime
from pathlib import Path

from genomeos.web import server
from genomeos.web.server import Api

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "data" / "results"


def _result(name: str) -> dict:
    return json.loads((RESULTS / f"{name}.json").read_text())


def _checkout(tmp_path: Path, results: tuple[str, ...] = (), docs: tuple[str, ...] = ()) -> Path:
    """A checkout with the plan and the README, and only the result files and docs named."""
    (tmp_path / "docs").mkdir()
    (tmp_path / "data" / "results").mkdir(parents=True)
    shutil.copy(ROOT / "docs" / "ROADMAP.md", tmp_path / "docs" / "ROADMAP.md")
    shutil.copy(ROOT / "README.md", tmp_path / "README.md")
    for d in docs:
        shutil.copy(ROOT / "docs" / d, tmp_path / "docs" / d)
    for r in results:
        shutil.copy(RESULTS / f"{r}.json", tmp_path / "data" / "results" / f"{r}.json")
    return tmp_path


def _cache_script():
    spec = importlib.util.spec_from_file_location("ci_status_cache", ROOT / "scripts" / "ci_status_cache.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _entries(bio: dict) -> list[dict]:
    return [e for c in bio["categories"] for e in c["validated"] + c["not_counted"]]


# -- each measure is computed from files --------------------------------------------------------------


def test_software_is_read_from_the_plan_the_engine_package_and_the_ci_cache(tmp_path):
    from genomeos import roadmap

    root = _checkout(tmp_path, results=("engine_package",))
    written = datetime(2026, 9, 28, 12, 0, tzinfo=UTC)
    cache = {
        "latest_completed": {
            "sha": "a" * 40,
            "conclusion": "failure",
            "event": "schedule",
            "created": "2026-09-28T06:00:00Z",
            "test_summary": "1 failed, 12 passed in 3.00s",
            "failed_tests": ["test_x"],
        },
        "refs": {"main": "b" * 40, "dev": "a" * 40},
        "main_check": {"sha": "b" * 40, "conclusion": None},
        "consecutive_failures": 3,
        "promote_dry_run": {"args": ["--dry-run", "--latest-green"], "exit": 1, "output": "refused"},
    }
    _cache_script().write(cache, root / "data" / "cache" / "ci_status.json", now=written)

    s = Api(root).state()
    sw = s["measures"]["software"]
    parsed = roadmap.parse_milestones((ROOT / "docs" / "ROADMAP.md").read_text())
    assert sw["milestones"]["reached"] == sum(1 for p in parsed if p["state"] == "done")
    assert sw["milestones"]["total"] == len(parsed)
    assert "not knowledge" in sw["definition"]

    first, second = sw["reviews"]
    assert first["source"].endswith("item 11") and second["source"].endswith("item 12")
    assert first["total"] == 9 and second["total"] == 8
    for rv in sw["reviews"]:
        assert rv["done"] + len(rv["reopened"]) + len(rv["closed"]) + len(rv["open"]) == rv["total"]

    engine = _result("engine_package")
    ep = sw["engine_package"]
    assert ep["checks_total"] == len(engine["checks"])
    assert ep["checks_passed"] == sum(1 for c in engine["checks"] if c["ok"])
    assert ep["dirty"] == engine["result_manifest"]["code"]["dirty"]

    # the test suite's status is the cached CI run's, with the cache's age beside it
    assert sw["tests"]["summary"] == cache["latest_completed"]["test_summary"]
    assert sw["tests"]["failed_tests"] == ["test_x"] and sw["tests"]["conclusion"] == "failure"
    assert sw["tests"]["age_seconds"] >= 0

    rel = s["measures"]["release"]
    checks = {c["check"]: c for c in rel["checks"]}
    assert checks["CI's latest completed run is green"]["ok"] is False
    assert checks["main sits on a revision with a green test run"]["ok"] is False
    assert checks["the promotion gate finds a green revision to promote"]["ok"] is False
    assert checks["the promotion gate finds a green revision to promote"]["detail"] == "refused"
    assert rel["ci"]["written"] == "2026-09-28T12:00:00Z" and rel["ci"]["consecutive_failures"] == 3
    assert set(rel["missing"]) | set(rel["met"]) | set(rel["unknown"]) == set(checks)


def test_release_readiness_quotes_the_plan_and_reads_packaging_and_licences_from_files():
    rel = Api(ROOT).state()["measures"]["release"]
    plan = " ".join((ROOT / "docs" / "ROADMAP.md").read_text().split())
    assert rel["criteria"], "the plan's own release criteria are quoted"
    for q in rel["criteria"]:
        assert re.sub(r"[*`]", "", q["text"])[:40] in re.sub(r"[*`]", "", plan)
    checks = {c["check"]: c for c in rel["checks"]}
    version = re.search(r'^version\s*=\s*"([^"]+)"', (ROOT / "pyproject.toml").read_text(), re.M).group(1)
    assert (
        version in checks["the version is tagged"]["detail"] or checks["the version is tagged"]["ok"] is None
    )
    assert checks["the licence files are present and shipped"]["ok"] is all(
        (ROOT / f).exists() for f in ("LICENSE", "LICENSE-APACHE", "NOTICE", "LICENSING.md")
    )
    engine = _result("engine_package")
    assert checks["the engine package builds and passes its own checks"]["ok"] is bool(
        engine["passed"] and all(c["ok"] for c in engine["checks"])
    )
    assert checks["that engine package was built from a committed revision"]["ok"] is (
        engine["result_manifest"]["code"]["dirty"] is False
    )
    assert checks["a PyPI package exists"]["ok"] is None  # not knowable from files; said so


def test_the_validated_figures_are_the_result_files_own():
    bio = Api(ROOT).state()["measures"]["biology"]
    by = {e["id"]: e for e in _entries(bio)}

    crispri = _result("crispri_published")
    f = by["crispri_heldout"]["figures"]
    assert f["gain"] == crispri["heldout_published_pairs"]["deletion_gain"]["gain"]
    assert f["gain_ci95"] == crispri["heldout_published_pairs"]["deletion_gain"]["ci95"]
    assert f["second_replicated"] is crispri["second_cell_type_hct116"]["replicated"]
    # the band field is not the claim's wording, so the validated entry does not carry it
    assert crispri["heldout_published_pairs"]["against_encode_re2g"] not in json.dumps(by["crispri_heldout"])

    audit = _result("node_independence_audit")
    assert by["node_containment"]["figures"]["excess_default"] == audit["selection_premium"]["default"]
    assert by["node_containment"]["figures"]["pairs"] == audit["pairs"]["total"]

    clause2 = _result("clause2_matched_control")
    assert (
        by["clause2"]["figures"]["matched_difference_points"]
        == clause2["primary"]["matched_difference_points"]
    )

    bench = _result("therapeutic_benchmark")
    assert by["therapeutic_benchmark"]["figures"]["targets_recovered"] == bench["targets_recovered"]


def test_no_figure_of_the_measures_is_written_into_the_code():
    s = Api(ROOT).state()
    source = (ROOT / "genomeos" / "web" / "server.py").read_text()
    page = (ROOT / "genomeos" / "web" / "static" / "index.html").read_text()
    values = []

    def walk(v):
        if isinstance(v, bool):
            return
        if isinstance(v, float) and len(repr(v).split(".")[-1]) >= 2:
            values.append(repr(v))
        elif isinstance(v, int) and abs(v) >= 100:
            values.append(str(v))
        elif isinstance(v, dict):
            for x in v.values():
                walk(x)
        elif isinstance(v, list):
            for x in v:
                walk(x)

    for e in _entries(s["measures"]["biology"]):
        walk(e["figures"])
    ep = s["measures"]["software"]["engine_package"]
    walk([ep["checks_passed"], ep["checks_total"]])
    assert values, "the biology entries carry figures read from their files"
    for v in values:
        # a figure standing as a number of its own; "190px" or "BIOMD0000000190" is not one
        figure = re.compile(rf"(?<![\w.]){re.escape(v)}(?![\w])")
        assert not figure.search(source), f"{v} is a result's figure in server.py"
        assert not figure.search(page), f"{v} is a result's figure in the page"


# -- a claim never appears without its qualification --------------------------------------------------


def test_every_shown_claim_carries_its_qualification_and_its_registered_wording():
    bio = Api(ROOT).state()["measures"]["biology"]
    readme = " ".join((ROOT / "README.md").read_text().split())
    shown = [e for e in _entries(bio) if e["class"] != "withheld"]
    assert shown
    for e in shown:
        assert e["wording"] and e["wording_source"]
        assert e["qualifications"] and all(q["text"] and q["source"] for q in e["qualifications"])
        assert e["withheld_because"] is None
    counted = [e for c in bio["categories"] for e in c["validated"]]
    assert all(e["counted"] for e in counted)

    crispri = next(e for e in shown if e["id"] == "crispri_heldout")
    # README's words, which say "in the range of", never the stronger band the file records
    assert crispri["wording"][:60] in readme and "in the range of" in crispri["wording"]
    quals = " ".join(q["text"] for q in crispri["qualifications"])
    verdict = _result("crispri_published")["second_cell_type_hct116"]["verdict"]
    assert verdict in quals and "one cell type" in quals and "not a paired test" in quals
    # lane-split's flag: the held-out file is a reused benchmark, and the panel says so
    audit = _result("crispri_split_audit")
    assert crispri["reused_benchmark"]["flag"] is True
    assert crispri["reused_benchmark"]["text"][:60] == audit["reused_benchmark"][:60]

    node = next(e for e in _entries(bio) if e["id"] == "node_containment")
    assert node["counted"] is False and "may not be called independent" in node["wording"]
    clause2 = next(e for e in _entries(bio) if e["id"] == "clause2")
    assert clause2["counted"] is False and clause2["class"] == "not met"
    assert any("not in prediction accuracy" in q["text"] for q in clause2["qualifications"])
    bench = next(e for e in _entries(bio) if e["id"] == "therapeutic_benchmark")
    assert bench["counted"] is False and bench["class"] == "built knowing the answers"


def test_a_claim_whose_qualification_cannot_be_read_is_withheld_not_shown(tmp_path):
    root = _checkout(tmp_path, results=("crispri_split_audit",), docs=("CRISPRI-RESULT.md",))
    crispri = _result("crispri_published")
    del crispri["second_cell_type_hct116"]  # the replication verdict the claim must travel with
    (root / "data" / "results" / "crispri_published.json").write_text(json.dumps(crispri))

    bio = Api(root).state()["measures"]["biology"]
    e = next(x for x in _entries(bio) if x["id"] == "crispri_heldout")
    assert e["class"] == "withheld" and e["counted"] is False
    assert e["wording"] is None and e["figures"] == {}
    assert "second_cell_type_hct116" in e["withheld_because"]
    target = next(c for c in bio["categories"] if c["category"] == "target")
    assert target["validated"] == [] and target["none"] is True
    assert "crispri_heldout" in bio["withheld"]


# -- a category with nothing validated reads none -----------------------------------------------------


def test_every_category_is_listed_and_an_empty_one_reads_none(tmp_path):
    names = [c for c, _ in server.BIOLOGY_CATEGORIES]
    bio = Api(ROOT).state()["measures"]["biology"]
    assert [c["category"] for c in bio["categories"]] == names
    for c in bio["categories"]:
        assert c["none"] is (not c["validated"])
    assert bio["counted"] == sum(len(c["validated"]) for c in bio["categories"])

    empty = Api(_checkout(tmp_path)).state()["measures"]["biology"]
    assert [c["category"] for c in empty["categories"]] == names  # nothing disappears
    assert all(c["none"] and not c["validated"] for c in empty["categories"])
    assert empty["counted"] == 0 and all(e["class"] == "withheld" for e in _entries(empty))

    page = (ROOT / "genomeos" / "web" / "static" / "index.html").read_text()
    assert "c.none ? '<span class=\"pp-none\">none</span>'" in page


# -- discontinued investigations are never capabilities -----------------------------------------------


def test_every_discontinued_investigation_is_quoted_from_its_own_plan_row():
    d = Api(ROOT).state()["discontinued"]
    assert [r["id"] for r in d["rows"]] == [x["id"] for x in server.DISCONTINUED]
    for r, spec in zip(d["rows"], server.DISCONTINUED, strict=True):
        assert r["state"] == "closed on a negative", f"{r['id']}: its row was not found in the plan"
        assert len(r["closed"]) == len(spec["closed"])
        for q, (_, phrase, *_) in zip(r["closed"], spec["closed"], strict=True):
            assert phrase in q["text"] and q["source"].startswith("docs/ROADMAP.md line ")
    by = {r["id"]: r for r in d["rows"]}
    # R8's transformations: retired AND reopened as item 13's pilot, both facts shown
    assert any("retired" in q["text"] for q in by["r8_transformations"]["closed"])
    assert any("reopened as S3" in q["text"] for q in by["r8_transformations"]["reopened"])
    assert any("pilot" in q["text"] for q in by["r8_transformations"]["reopened"])
    assert by["igvf_mhc"]["result"]["verdict"] == _result("indep_mhc_crispri")["verdict"]


def test_a_discontinued_investigation_is_never_counted_as_done():
    s = Api(ROOT).state()
    ids = [x["id"] for x in server.DISCONTINUED]
    assert s["measures"]["software"]["not_counted"] == ids
    items = {r["item"]: r for r in s["review"]["rows"]}
    for spec in server.DISCONTINUED:
        if spec["review_item"]:
            assert items[spec["review_item"]]["state"] != "done"
    assert items["R8"]["state"] == "reopened"
    assert s["review"]["done"] == sum(1 for r in s["review"]["rows"] if r["state"] == "done")
    # a follow-up row that records a negative never makes its item done by itself
    for r in s["review"]["rows"]:
        if r["state"] == "done":
            assert any("done" in f["state"].lower() for f in r["follow_ups"] if not f["discontinued"])


def test_a_negative_row_alone_leaves_its_item_open():
    anchor = next(x for x in server.DISCONTINUED if x["id"] == "human_rate")["closed"][0][0]
    plan = "\n".join(
        [
            "11. **The external review, for this test.**",
            "   | # | Item | Code | Acceptance | Order |",
            "   | --- | --- | --- | --- | --- |",
            "   | R3 | A bridge | `x.py` | a fixture | now |",
            f"   | ↳ R3 | **{anchor}, lane-x)** | `x.py` | nothing to borrow | done |",
            "   | R8 | An engine | `y.py` | a pretest | later |",
            "   | ↳ R8 | **Pretest: NEGATIVE** | `y.py` | closes | closed |",
            "   | ↳ R8 | *Correction* | — | reopened as S3 |",
            "12. Next item.",
        ]
    )
    rv = server._review(plan)
    by = {r["item"]: r for r in rv["rows"]}
    assert by["R3"]["follow_ups"][0]["discontinued"] is True and by["R3"]["state"] == "open"
    assert by["R8"]["state"] == "reopened" and rv["done"] == 0 and rv["reopened"] == ["R8"]


# -- the CI status cache -------------------------------------------------------------------------------


def test_the_ci_cache_script_reads_github_through_one_runner_and_only_dry_runs_the_gate(tmp_path):
    mod = _cache_script()
    calls: list[list[str]] = []
    log = "\n".join(
        [
            "test\tUNKNOWN STEP\t2026-09-28T13:28:46.2387846Z ________ test_a_partner ________",
            "test\tUNKNOWN STEP\t2026-09-28T13:28:46.5965591Z 1 failed, 20 passed in 5.10s (0:00:05)",
        ]
    )

    def fake(cmd):
        calls.append(cmd)
        if cmd[:3] == ["gh", "repo", "view"]:
            return 0, "owner/repo\n", ""
        if cmd[0] == "scripts/promote_main.sh":
            return 1, "", "promote_main: refused: no green run"
        if cmd[:3] == ["gh", "run", "list"]:
            runs = [
                {"databaseId": 2, "headSha": "d" * 40, "headBranch": "dev", "event": "schedule",
                 "status": "completed", "conclusion": "failure", "createdAt": "2026-09-28", "url": "u2"},
                {"databaseId": 1, "headSha": "c" * 40, "headBranch": "dev", "event": "schedule",
                 "status": "completed", "conclusion": "success", "createdAt": "2026-09-27", "url": "u1"},
            ]  # fmt: skip
            return 0, json.dumps(runs), ""
        if cmd[:3] == ["gh", "run", "view"] and "--json" in cmd:
            return 0, json.dumps({"jobs": [{"name": "test", "conclusion": "failure", "databaseId": 9}]}), ""
        if cmd[:3] == ["gh", "run", "view"] and "--log" in cmd:
            return 0, log, ""
        if cmd[:2] == ["git", "rev-parse"]:
            return 0, ("m" if cmd[-1].endswith("main") else "d") * 40 + "\n", ""
        if cmd[:2] == ["gh", "api"]:
            return 0, json.dumps({"check_runs": []}), ""
        raise AssertionError(cmd)

    status = mod.collect(sh=fake)
    assert all("--push" not in c for c in calls)
    assert ["scripts/promote_main.sh", "--dry-run", "--latest-green"] in calls
    latest = status["latest_completed"]
    assert latest["conclusion"] == "failure" and latest["test_summary"].startswith("1 failed, 20 passed")
    assert latest["failed_tests"] == ["test_a_partner"]
    assert status["latest_success"]["sha"] == "c" * 40 and status["consecutive_failures"] == 1
    assert status["promote_dry_run"]["exit"] == 1 and "refused" in status["promote_dry_run"]["output"]
    assert status["main_check"]["conclusion"] is None and status["errors"] == []

    path = mod.write(
        status, tmp_path / "data" / "cache" / "ci_status.json", now=datetime(2026, 9, 28, tzinfo=UTC)
    )
    ci = server._ci_status(tmp_path)
    assert json.loads(path.read_text())["written"] == "2026-09-28T00:00:00Z"
    assert (
        ci["present"] and ci["green"] is False and ci["tested_dev_tip"] is True and ci["main_green"] is False
    )
    assert ci["age_seconds"] > 0
    assert server._ci_status(tmp_path / "nowhere")["present"] is False
