import json
import threading
import urllib.request
from pathlib import Path

import pytest

from genomeos.web.server import Api, ApiError, make_server

ROOT = Path(__file__).resolve().parent.parent
SRC = (ROOT / "data/demo/repressilator.bio").read_text()


def test_api_files_and_module():
    api = Api(ROOT)
    f = api.files()
    assert any(m["path"].endswith("repressilator.bio") for m in f["modules"])
    assert "module synthetic.repressilator" in api.module_source("data/demo/repressilator.bio")["source"]


def test_api_rejects_paths_outside_data():
    api = Api(ROOT)
    with pytest.raises(ApiError):
        api.module_source("pyproject.toml")
    with pytest.raises(ApiError):
        api.module_source("../../etc/passwd")


def test_api_genome_and_orfs():
    api = Api(ROOT)
    info = api.genome_info("data/demo/demo.fa")
    assert info["total_bp"] == 528
    orfs = api.genome_orfs("data/demo/demo.fa", min_aa=50)
    assert orfs["orfs"][0]["length_aa"] == 61
    seq = api.genome_fetch("data/demo/demo.fa", "demo_synthetic:150-153(+)")
    assert seq["sequence"] == "ATG"


def test_api_compile_and_run():
    api = Api(ROOT)
    c = api.compile(SRC)
    assert c["ok"] and c["report"]["rules"] == 6
    bad = api.compile("module x\ngene a { produces: Nope }")
    assert not bad["ok"] and "undeclared" in bad["error"]
    r = api.run(SRC, hours=60, dt=0.01, initial={"TetR": 10})
    assert r["ok"] and r["peaks"]["TetR"] >= 3 and len(r["times"]) <= 601


def test_api_age():
    a = Api(ROOT).age(cell_type="fibroblast", years=50, cells=100, seed=2)
    assert a["reports"][-1]["telomere_bp"] < a["reports"][0]["telomere_bp"]
    assert all(p["source"] for p in a["evidence"])


def test_http_roundtrip():
    srv = make_server(ROOT, port=0)
    port = srv.server_address[1]
    t = threading.Thread(target=srv.serve_forever, daemon=True)
    t.start()
    try:
        with urllib.request.urlopen(f"http://127.0.0.1:{port}/") as r:
            assert b"<title>GenomeOS</title>" in r.read()
        with urllib.request.urlopen(f"http://127.0.0.1:{port}/api/libs") as r:
            assert len(json.load(r)["libraries"]) > 40
        req = urllib.request.Request(
            f"http://127.0.0.1:{port}/api/compile",
            data=json.dumps({"source": SRC}).encode(),
            headers={"Content-Type": "application/json"},
        )
        with urllib.request.urlopen(req) as r:
            assert json.load(r)["ok"]
        with pytest.raises(urllib.error.HTTPError) as ei:
            urllib.request.urlopen(f"http://127.0.0.1:{port}/api/module?path=pyproject.toml")
        assert ei.value.code == 403
    finally:
        srv.shutdown()
        srv.server_close()


def test_api_grow_reports_the_organism_and_its_space():
    api = Api(ROOT)
    worm = api.grow("data/demo/celegans_lineage.bio", "150", 2)
    assert (
        worm["summary"]["organism"] == "CelegansEarly"
        and worm["space"] is None
        and worm["diff"]["matched"] > 0
    )
    flag = api.grow("data/demo/flag_organism.bio", "200 h", 0)
    rows = flag["space"]["rows"]
    assert (
        flag["space"]["width"] == 30
        and len(rows) == 1
        and rows[0].startswith("BBB")
        and rows[0].endswith("RRR")
    )
    assert [b[0] for b in flag["space"]["bands"]] == ["Blue", "White", "Red"] and flag["asserts"][0]["ok"]


def test_api_budget_wide(tmp_path):
    """The 98% card: per-chromosome budgets summed, the job's state attached, empty when nothing ran."""
    empty = Api(tmp_path).budget_wide()
    assert empty["done"] == 0 and empty["table"] == [] and empty["tiers"]["unconstrained_unknown"] == 0
    assert "neutral" not in empty["tiers"] and "fossil" not in empty["tiers"]
    (tmp_path / "data" / "results").mkdir(parents=True)
    (tmp_path / "data" / "results" / "budget_chr21.json").write_text(
        json.dumps(
            {
                "chrom": "chr21",
                "unknown_bp": 1000,
                "chromosome_length": 4000,
                "constrained_bp": 9,
                "measured_bp": 900,
                "constrained_fraction": 0.01,
                "blocks": [{}, {}],
                "by_tier": {
                    "structural": {"blocks": 1, "bp": 100},
                    "fossil": {"blocks": 0, "bp": 0},
                    "regulatory": {"blocks": 0, "bp": 0},
                    "constrained_unknown": {"blocks": 1, "bp": 300},
                    "neutral": {"blocks": 0, "bp": 600},
                },
                "cost": {"phylop": {"mb_fetched": 3.0}, "seconds": 12.5},
            }
        )
    )
    b = Api(tmp_path).budget_wide()
    assert b["done"] == 1 and b["total"] == 24 and b["unknown_bp"] == 1000 and b["genome_bp"] == 4000
    assert b["tiers"]["constrained_unknown"] == 300 and b["constrained_fraction"] == 0.01
    # a stored budget with no axes record reads under the R7 names, numbers as stored
    assert b["tiers"]["unconstrained_unknown"] == 600 and "neutral" not in b["tiers"]
    row = b["table"][0]
    assert row["blocks"] == 2 and row["constrained_unknown_blocks"] == 1 and row["mb_fetched"] == 3.0
    assert b["job"] is None or b["job"]["name"] == "budget_genome_wide"


def test_api_state_milestones_are_the_roadmap_parser_s():
    """The Progress tab's milestone states are the parser's, not a second reading of the same table."""
    from genomeos import roadmap

    s = Api(ROOT).state()
    parsed = roadmap.parse_milestones((ROOT / "docs" / "ROADMAP.md").read_text())
    rows = s["milestones"]["rows"]
    assert [(r["milestone"], r["state"]) for r in rows] == [(p["milestone"], p["state"]) for p in parsed]
    assert s["milestones"]["reached"] == sum(1 for p in parsed if p["state"] == "done")
    assert s["milestones"]["held"] == [p["milestone"] for p in parsed if p["state"] != "done"]
    # a held milestone carries the reason its own row states; a reached one claims no reason
    for r in rows:
        assert (r["reason"] is None) == (r["state"] == "done")
    held = [r for r in rows if r["state"] != "done"]
    assert all(len(r["reason"]) > 20 for r in held)


def test_api_state_review_reports_the_nine_items_with_their_acceptance_tests():
    s = Api(ROOT).state()["review"]
    assert [r["item"] for r in s["rows"]] == [f"R{i}" for i in range(1, 10)]
    assert all(r["acceptance"] and r["title"] for r in s["rows"])
    assert all(r["state"] in ("done", "open", "reopened", "closed") for r in s["rows"])
    assert s["done"] == sum(1 for r in s["rows"] if r["state"] == "done")
    # an item is only done because a follow-up row in the plan says so
    assert all(r["follow_ups"] for r in s["rows"] if r["state"] == "done")
    # and a row the plan reopened is not done: R8 closed on a negative pretest, then reopened as S3
    for r in s["rows"]:
        if any("reopened" in f["state"].lower() for f in r["follow_ups"]):
            assert r["state"] == "reopened" and r["item"] in s["reopened"]


def test_api_state_claims_read_their_result_files_and_never_restate_a_figure():
    """The honest panel: README's words, the result files' numbers, and no figure written by hand."""
    s = Api(ROOT).state()["claims"]
    readme = (ROOT / "README.md").read_text()
    assert [r["claim"] for r in s["rows"]] == [
        "Sequence",
        "Annotation",
        "Assay coverage",
        "Software correctness",
        "Independent prediction",
    ]
    for r in s["rows"]:
        assert r["text"][:40] in " ".join(readme.split())  # the words are README's, unchanged

    coverage = json.loads((ROOT / "data/results/unknown_coverage.json").read_text())
    fig = next(r["figure"] for r in s["rows"] if r["claim"] == "Assay coverage")
    real = coverage["real_unknown"]
    assert fig["measured_bp"] == real["measured_bp"] and fig["total_bp"] == real["bp"]
    assert fig["percent"] == round(real["measured_bp"] / real["bp"] * 100, 2)
    assert fig["date"] == coverage["date"] and fig["source"].endswith("unknown_coverage.json")
    # README quotes this result, so a README that stops matching its file is reported, not believed
    assert next(r["agrees_with_readme"] for r in s["rows"] if r["claim"] == "Assay coverage") is True

    crispri = json.loads((ROOT / "data/results/crispri_published.json").read_text())
    pred = next(r["figure"] for r in s["rows"] if r["claim"] == "Independent prediction")
    held = crispri["heldout_published_pairs"]
    assert pred["held_out_gain"] == held["deletion_gain"]["gain"]
    assert pred["ci95"] == held["deletion_gain"]["ci95"]
    assert pred["against_encode_re2g"] == held["against_encode_re2g"]
    # the qualification travels with the number: one cell type, and the second one did not replicate
    assert pred["qualification"] == crispri["second_cell_type_hct116"]["verdict"]
    assert pred["replication"]["replicated"] is crispri["second_cell_type_hct116"]["replicated"]

    source = (ROOT / "genomeos" / "web" / "server.py").read_text()
    page = (ROOT / "genomeos" / "web" / "static" / "index.html").read_text()
    for stale in ("160447", "160,447", "30602182", "30,602,182", "0.52%", "0.1095"):
        assert stale not in source, f"{stale} is a result's figure; read it from the file"
        assert stale not in page, f"{stale} is a result's figure; the page may not carry it"


def test_api_state_renders_with_an_empty_work_board(tmp_path):
    """A checkout with the plan and the README but nobody working and nothing computed."""
    (tmp_path / "docs").mkdir()
    (tmp_path / "docs" / "ROADMAP.md").write_text((ROOT / "docs" / "ROADMAP.md").read_text())
    (tmp_path / "README.md").write_text((ROOT / "README.md").read_text())
    api = Api(tmp_path)
    assert api.work()["board"] == [] and api.work()["uncommitted"] == []
    s = api.state()
    assert s["milestones"]["total"] and s["review"]["total"] == 9
    # every claim is still stated; the figures are absent rather than invented
    assert len(s["claims"]["rows"]) == 5
    assert all(r["figure"] is None and r["agrees_with_readme"] is None for r in s["claims"]["rows"])
    assert s["claims"]["caveats"] == []
    assert s["owed"]["board"] == [] and s["owed"]["ledger"] is None
    assert s["owed"]["count"] == len(s["owed"]["data_jobs"]) + len(s["owed"]["steps"])
