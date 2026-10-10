# SPDX-License-Identifier: AGPL-3.0-or-later
"""The read-by-marks call shown in the Blocks lane (genome/readcall.py, /api/blocks/readcall,
static/readcall.js).

The reader calls a promoter read when a DNase peak covers it, or, where the DNase file has none, when
H3K4me3 and H3K27ac peaks do. The Blocks lane's `<cell>_read` folded the two together, so the call made by
the marks was shown nowhere. These tests plant one gene for each call, and an UNKNOWN for each reason, and
check that the API keeps them apart, that `<cell>_read` is unchanged beside the new `<cell>_read_basis`,
that the card shows UNKNOWN with its reason and never as closed, and that it says "read", not more."""

import json
import shutil
import subprocess
import threading
import urllib.request
from pathlib import Path
from types import SimpleNamespace

import pytest

from genomeos.genome.readcall import CALLS, GLOSS, REASONS, read_calls
from genomeos.web.server import Api, ApiError, make_server

ROOT = Path(__file__).resolve().parents[1]
STATIC = ROOT / "genomeos" / "web" / "static"

#: one coding chromosome of five genes, the reader's own count of them, and a lncRNA beside them
CODING = ["APP", "VPS26C", "BACH1", "TPTE", "DYRK1A"]
TYPES = {**dict.fromkeys(CODING, "protein_coding"), "LINC1": "lncRNA"}


def _reader(cell: str, **lists) -> dict:
    return {
        "cell_type": cell,
        "chrom": "chr21",
        "coding_genes": lists.pop("coding_genes", len(CODING)),
        "top_read": lists.pop("top_read", []),
        "silent_genes": lists.pop("silent_genes", []),
        "poised_genes": lists.pop("poised_genes", []),
        "read_by_marks": lists.pop("read_by_marks", []),
        "marks_used": lists.pop("marks_used", True),
        "node_table": [],
        "silent_node_ids": [],
    }


def _plant(results: Path) -> Path:
    """K562 reads APP on the DNase file (by absence) and VPS26C with a signal, holds BACH1 poised and TPTE
    closed; GM12878 reads APP by its marks; HepG2 counted six coding genes, so absence decides nothing
    there; spinal cord is in the genome-wide panel with no file on chr21; a comparison file is skipped."""
    results.mkdir(parents=True, exist_ok=True)
    (results / "reader_genome_wide.json").write_text(
        json.dumps({"cell_types": ["K562", "GM12878", "HepG2", "spinal cord"]})
    )
    files = {
        "K562": _reader(
            "K562",
            top_read=[{"gene": "VPS26C", "signal": 12.47}],
            silent_genes=["TPTE", "BACH1"],  # poised sits inside silent_genes by design
            poised_genes=["BACH1"],
        ),
        "GM12878": _reader("GM12878", read_by_marks=["APP"], silent_genes=["TPTE"]),
        "HepG2": _reader("HepG2", coding_genes=6, silent_genes=["TPTE"]),
    }
    for cell, r in files.items():
        (results / f"reader_{cell}_chr21.json").write_text(json.dumps(r))
    (results / "reader_K562_vs_HepG2_chr21.json").write_text(json.dumps({"only_in_a": ["APP"]}))
    return results


def _calls(results: Path, symbol: str) -> dict[str, dict]:
    d = read_calls(
        symbol, "chr21", TYPES.get(symbol), len(CODING), results, reader_terms={"K562", "GM12878", "HepG2"}
    )
    return {c["biosample"]: c for c in d["calls"]}


def test_PLANTED_marks_read_dnase_read_poised_closed_and_unknown_are_kept_apart(tmp_path):
    rd = _plant(tmp_path)

    app = _calls(rd, "APP")
    assert list(app) == ["K562", "GM12878", "HepG2", "spinal cord"]  # the panel's order; no comparison file
    # the same gene read in two biosamples on two different grounds
    assert app["GM12878"]["call"] == "read_by_marks" and app["GM12878"]["basis"] == "marks"
    assert app["GM12878"]["read"] is True
    assert app["K562"]["call"] == "open" and app["K562"]["basis"] == "DNase" and app["K562"]["by_absence"]
    # HepG2 counted six coding genes against the annotation's five: absence is not read as open
    assert app["HepG2"]["call"] == "unknown" and app["HepG2"]["reason"] == "coding_set_differs"
    assert app["HepG2"]["read"] is None and "6 coding genes" in app["HepG2"]["reason_text"]
    # no file for spinal cord on chr21: UNKNOWN with its reason, never closed
    assert app["spinal cord"]["call"] == "unknown" and app["spinal cord"]["reason"] == "no_reader_result"
    assert app["spinal cord"]["in_reader_terms"] is False and app["K562"]["in_reader_terms"] is True

    vps = _calls(rd, "VPS26C")["K562"]
    assert vps["call"] == "open" and vps["signal"] == 12.47 and not vps["by_absence"]
    bach1 = _calls(rd, "BACH1")["K562"]
    assert bach1["call"] == "poised" and bach1["read"] is False
    tpte = _calls(rd, "TPTE")
    assert {c["call"] for c in tpte.values() if c["reason"] != "no_reader_result"} == {"closed"}

    linc = _calls(rd, "LINC1")
    assert {c["call"] for c in linc.values()} == {"unknown"}
    assert {c["reason"] for c in linc.values()} == {"not_protein_coding"}
    assert "lncRNA" in linc["K562"]["reason_text"]
    nowhere = _calls(rd, "NOT_A_GENE")
    assert {c["reason"] for c in nowhere.values()} == {"not_in_annotation"}

    d = read_calls("APP", "chr21", "protein_coding", len(CODING), rd, reader_terms={"K562"})
    assert d["counts"] == {"open": 1, "read_by_marks": 1, "poised": 0, "closed": 0, "unknown": 2}
    assert d["read"] == 2 and d["biosamples"] == 4 and set(d["counts"]) == set(CALLS)


def test_PLANTED_marks_not_read_leave_only_the_dnase_call_and_say_so(tmp_path):
    rd = tmp_path / "r"
    rd.mkdir()
    (rd / "reader_K562_chr21.json").write_text(
        json.dumps(_reader("K562", silent_genes=["TPTE"], marks_used=False))
    )
    tpte = _calls(rd, "TPTE")["K562"]
    assert tpte["call"] == "closed" and tpte["basis"] == "DNase" and "not tested" in tpte["note"]
    assert _calls(rd, "APP")["K562"]["basis"] == "DNase"


def test_the_blocks_attribute_says_which_read_rests_on_the_marks_and_leaves_read_unchanged(tmp_path):
    from genomeos.genome.blocks import Block, _apply_reader

    rd = _plant(tmp_path)
    gene = {"gene_type": "protein_coding"}
    blocks = [
        Block("g1", "gene", 10, 500, "+", "APP", None, "curated", 0.9, dict(gene)),
        Block("g2", "gene", 600, 900, "-", "TPTE", None, "curated", 0.9, dict(gene)),
        Block("g3", "gene", 950, 990, "+", "BACH1", None, "curated", 0.9, dict(gene)),
        Block("g4", "gene", 995, 999, "+", "LINC1", None, "curated", 0.9, {"gene_type": "lncRNA"}),
    ]
    assert _apply_reader(blocks, "chr21", rd) == ["GM12878", "HepG2", "K562"]
    app, tpte, bach1, linc = (b.attrs for b in blocks)
    assert app["GM12878_read"] == "read" and app["GM12878_read_basis"] == "marks"
    assert app["K562_read"] == "read" and app["K562_read_basis"] == "DNase"
    assert tpte["K562_read"] == "silent" and "K562_read_basis" not in tpte
    assert bach1["K562_read"] == "poised" and "K562_read_basis" not in bach1
    assert not any(k.endswith("_read_basis") or k.endswith("_read") for k in linc)


def _api(tmp_path, monkeypatch) -> Api:
    _plant(tmp_path / "data" / "results")
    api = Api(tmp_path)
    genes = {
        s: SimpleNamespace(symbol=s, type=t, locus=SimpleNamespace(chrom="chr21")) for s, t in TYPES.items()
    }
    monkeypatch.setattr(api, "_annotation_for", lambda rel, chrom: SimpleNamespace(genes=genes))
    return api


def test_api_blocks_readcall_answers_the_planted_cases_and_refuses_a_gene_it_does_not_have(
    tmp_path, monkeypatch
):
    api = _api(tmp_path, monkeypatch)
    d = api.blocks_readcall("", "chr21", "APP")
    by = {c["biosample"]: c for c in d["calls"]}
    assert by["GM12878"]["call"] == "read_by_marks" and by["K562"]["call"] == "open"
    assert by["spinal cord"]["call"] == "unknown" and by["spinal cord"]["in_reader_terms"] is False
    assert d["gene_type"] == "protein_coding" and d["gloss"] == GLOSS and d["reasons"] == REASONS
    assert api.blocks_readcall("", "chr21", "app")["gene"] == "APP"
    assert {c["reason"] for c in api.blocks_readcall("", "chr21", "LINC1")["calls"]} == {"not_protein_coding"}
    with pytest.raises(ApiError) as e:
        api.blocks_readcall("", "chr21", "NOPE")
    assert e.value.status == 404
    with pytest.raises(ApiError):
        api.blocks_readcall("", "chr21", " ")


def test_http_readcall_on_the_committed_reader_results_for_app():
    """APP on chr21: read by its marks in GM12878 and keratinocyte, on the DNase file in K562; the spinal
    cord reading is shown and flagged as outside the thirteen reader cells. The route and the script are
    served."""
    srv = make_server(ROOT, port=0)
    port = srv.server_address[1]
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    try:
        with urllib.request.urlopen(f"http://127.0.0.1:{port}/api/blocks/readcall?gene=APP&chrom=chr21") as r:
            d = json.load(r)
        with urllib.request.urlopen(f"http://127.0.0.1:{port}/static/readcall.js") as r:
            assert b"window.readCallHtml" in r.read()
    finally:
        srv.shutdown()
        srv.server_close()
    by = {c["biosample"]: c for c in d["calls"]}
    assert d["biosamples"] == 14 and len(by) == 14
    assert by["GM12878"]["call"] == "read_by_marks" and by["keratinocyte"]["call"] == "read_by_marks"
    assert by["K562"]["call"] == "open" and by["K562"]["basis"] == "DNase"
    assert [c for c, x in by.items() if not x["in_reader_terms"]] == ["spinal cord"]
    assert d["counts"]["unknown"] == 0 and d["read"] == d["counts"]["open"] + d["counts"]["read_by_marks"]


def test_the_lane_draws_marks_apart_and_loads_the_card_without_touching_the_page():
    js = (STATIC / "blocks.js").read_text()
    assert "b.attrs[st.reader + '_read_basis'] === 'marks'" in js and "◐ read by marks" in js
    assert "s.src = '/static/readcall.js'" in js and 'id="b-readcall-out"' in js
    assert "/api/blocks/readcall?gene=" in js and "window.readCallHtml(r, st.reader)" in js
    assert "readcall.js" not in (STATIC / "index.html").read_text()


def _render(payloads: list[dict], selected: str = "") -> list[str]:
    """window.readCallHtml from static/readcall.js, run in node."""
    harness = (
        "const vm = require('vm'), fs = require('fs'); const window = {};"
        "vm.runInNewContext(fs.readFileSync(process.argv[1], 'utf8'), {window});"
        "const ps = JSON.parse(fs.readFileSync(0, 'utf8'));"
        "process.stdout.write(JSON.stringify(ps.map(p => window.readCallHtml(p, process.argv[2]))));"
    )
    run = subprocess.run(
        ["node", "-e", harness, str(STATIC / "readcall.js"), selected],
        input=json.dumps(payloads),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert run.returncode == 0, run.stderr
    return json.loads(run.stdout)


@pytest.mark.skipif(shutil.which("node") is None, reason="node is not installed")
def test_the_card_shows_marks_apart_unknown_with_its_reason_and_says_read_not_more(tmp_path):
    rd = _plant(tmp_path)
    terms = {"K562", "GM12878", "HepG2"}
    app = read_calls("APP", "chr21", "protein_coding", len(CODING), rd, reader_terms=terms)
    linc = read_calls("LINC1", "chr21", "lncRNA", len(CODING), rd, reader_terms=terms)
    tpte = read_calls("TPTE", "chr21", "protein_coding", len(CODING), rd, reader_terms=terms)
    html_app, html_linc, html_tpte = _render([app, linc, tpte], "GM12878")
    rows = {r.split("<td>")[1].split(" <span")[0].split("</td>")[0]: r for r in html_app.split("<tr ")[1:]}

    assert "◐</span> read · by marks" in rows["GM12878"] and "<td>marks</td>" in rows["GM12878"]
    assert "background:" in rows["GM12878"]  # the biosample the map shows
    assert "●</span> read · DNase" in rows["K562"] and "<td>DNase</td>" in rows["K562"]
    assert "(<b>1</b> on DNase, <b>1</b> by marks)" in html_app
    for who, reason in (("HepG2", "coding_set_differs"), ("spinal cord", "no_reader_result")):
        assert '<span class="pill none">UNKNOWN</span>' in rows[who]
        assert f'<span class="mono rc-reason">{reason}</span>' in rows[who] and "closed" not in rows[who]
    assert "outside the 13 reader cells" in rows["spinal cord"] and "outside" not in rows["K562"]
    assert html_linc.count("not_protein_coding") == 4 and "closed" not in html_linc.split("</div>", 2)[2]
    assert "closed · not read" in html_tpte
    # accessibility is not function: the card and its source never say active or expressed
    src = (STATIC / "readcall.js").read_text().lower()
    for out in (html_app, html_linc, html_tpte, src):
        assert "active" not in out.lower() and "expressed" not in out.lower()
