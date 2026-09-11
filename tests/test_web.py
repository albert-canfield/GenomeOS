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
