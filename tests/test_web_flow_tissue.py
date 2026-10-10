# SPDX-License-Identifier: AGPL-3.0-or-later
"""The Flow tab draws the tissue's own isoform beside the canonical one (static/flow.js, /api/flow?tissue=).

The page carries a tissue selector, a second drawing and an output block; the API gives the selector its
tissues and the tissue's isoform its geometry; the renderer shows an UNKNOWN call's reason verbatim and
never presents the canonical as the tissue's."""

import json
import shutil
import subprocess
from pathlib import Path

import pytest

from genomeos.coords import Locus
from genomeos.flow import trace_gene, trace_gene_in_tissue
from genomeos.flow.tissue import REASONS
from genomeos.genome import Chromosome, Genome, Sequence
from genomeos.ir.model import Gene, Transcript

ROOT = Path(__file__).resolve().parents[1]
STATIC = ROOT / "genomeos" / "web" / "static"

#   0..9 exon1 (GGG UTR + ATG AAA) | 9..19 intron | 19..25 exon2 (GAT CCC) | 25..35 intron |
#   35..44 exon3 (TGG TAA + CCC UTR); the cerebellum's isoform skips exon 2
SEQ = "GGGATGAAA" + "GTTTTTTTAG" + "GATCCC" + "GTTTTTTTAG" + "TGGTAACCC"
E1, E2, E3 = Locus("c", 0, 9), Locus("c", 19, 25), Locus("c", 35, 44)
ISO = {
    "gene": "G1",
    "isoforms": {
        "ENST1": {"tissues": {"Liver": 10.0, "Brain_Cerebellum": 2.0, "Lung": 5.0, "Skin": 0.5}},
        "ENST2": {"tissues": {"Liver": 1.0, "Brain_Cerebellum": 8.0, "Lung": 4.8, "Skin": 0.2}},
        "ENST3": {"tissues": {"Heart": 9.0}},  # not in the local models
        "ENST4": {"tissues": {"Muscle": 7.0}},  # non-coding in the local models
    },
}
#: one synthetic call per reason the API can answer (None: no isoform table for the gene)
UNKNOWN_CASES = {
    "no_isoform_data": ("Liver", None),
    "tissue_not_measured": ("Kidney", ISO),
    "low_expression": ("Skin", ISO),
    "tie": ("Lung", ISO),
    "dominant_not_in_local_models": ("Heart", ISO),
    "dominant_non_coding": ("Muscle", ISO),
}


def _payload(tissue: str, isoforms, variant=(36, "G", "A")) -> dict:
    """What /api/flow returns for the synthetic gene, assembled the way Api.flow assembles it."""
    g = Genome("t")
    g.chromosomes["c"] = Chromosome("c", Sequence(SEQ))
    gene = Gene(id="g1", kind="gene", symbol="G1", locus=Locus("c", 0, 44))
    txs = [
        Transcript(
            id="ENST1.2", kind="transcript", exons=[E1, E2, E3],
            cds_segments=[Locus("c", 3, 9), E2, Locus("c", 35, 41)], tags=["Ensembl_canonical"],
        ),
        Transcript(
            id="ENST2.5", kind="transcript", exons=[E1, E3],
            cds_segments=[Locus("c", 3, 9), Locus("c", 35, 41)],
        ),
        Transcript(id="ENST4.1", kind="transcript", exons=[E1, E2]),
    ]  # fmt: skip
    tr, tt = trace_gene(g, gene, txs), trace_gene_in_tissue(g, gene, txs, tissue, isoforms)
    out = tr.to_dict()
    out["gene_start"], out["gene_end"] = 0, 44
    out["tissue"] = tt.to_dict()
    if tt.isoform is not None:
        out["tissue"]["trace"] = tt.isoform.to_dict()
    out["variant"] = tr.substitute(*variant)
    out["variant_in_tissue"] = tt.substitute(*variant)
    return json.loads(json.dumps(out))


def _render(payloads: list[dict]) -> list[str]:
    """window.flowTissueHtml from static/flow.js, run in node with the DOM stubbed out."""
    harness = (
        "const vm = require('vm'), fs = require('fs');"
        "const window = {};"
        "const document = {addEventListener() {}, querySelector() { return null; }, documentElement: {}};"
        "vm.runInNewContext(fs.readFileSync(process.argv[1], 'utf8'), {window, document});"
        "const ps = JSON.parse(fs.readFileSync(0, 'utf8'));"
        "process.stdout.write(JSON.stringify(ps.map(p => window.flowTissueHtml(p))));"
    )
    run = subprocess.run(
        ["node", "-e", harness, str(STATIC / "flow.js")],
        input=json.dumps(payloads),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert run.returncode == 0, run.stderr
    return json.loads(run.stdout)


def test_the_page_carries_the_tissue_selector_and_the_unknown_reason_hook():
    page = (STATIC / "index.html").read_text()
    js = (STATIC / "flow.js").read_text()
    flow = page[page.index('<section id="flow">') : page.index("<!-- ================= CANCER")]
    assert '<select id="f-tissue"' in flow and '<option value="">canonical only</option>' in flow
    assert 'id="f-tcanvas"' in flow and 'id="f-tissue-out"' in flow
    # the selector feeds the API and is filled from it; the output block is the UNKNOWN renderer
    assert "'&tissue=' + encodeURIComponent(tissue)" in js and "d.tissues_measured" in js
    assert "$('#f-tissue-out').innerHTML = tissueHtml(d)" in js and "window.flowTissueHtml = tissueHtml" in js
    assert 'class="mono f-tissue-reason">${esc(t.reason)}</span>' in js
    # every reason the rule can give has its plain-words gloss beside the verbatim code
    assert all(f"{r}: " in js for r in REASONS)
    # the tissue's drawing comes only from its own trace, and only when the call is ESTABLISHED
    assert "t.status === 'ESTABLISHED' ? t.trace : null" in js


@pytest.mark.skipif(shutil.which("node") is None, reason="node is not installed")
def test_the_renderer_shows_the_isoform_or_the_reason_verbatim_never_the_canonical_in_its_place():
    established = _payload("Brain_Cerebellum", ISO)
    assert established["tissue"]["status"] == "ESTABLISHED" and "trace" in established["tissue"]
    unknown = {r: _payload(*args) for r, args in UNKNOWN_CASES.items()}
    assert {r: p["tissue"]["reason"] for r, p in unknown.items()} == {r: r for r in REASONS}
    assert not any("trace" in p["tissue"] for p in unknown.values())

    html = _render([established, *unknown.values()])
    est = html[0]
    iso, canon = established["tissue"]["transcript"], established["tissue"]["canonical_transcript"]
    assert f"dominant isoform <b>{iso}</b>" in est and "lacks canonical exon 2" in est
    # the stop is residue 3 of the tissue's protein and residue 5 of the canonical's, each named as such
    assert f"on {iso}: c." in est and "p.Trp3Ter" in est
    assert f"on the canonical {canon}: " in est and est.index("p.Trp5Ter") > est.index("on the canonical")
    assert (
        "residue_changed" in est
        and "measured: GTEx v8 median" in est
        and "predicted by the local trace" in est
    )

    for reason, out in zip(unknown, html[1:], strict=True):
        assert (
            f'<span class="pill none">UNKNOWN</span> <span class="mono f-tissue-reason">{reason}</span>'
            in out
        )
        assert "dominant isoform" not in out and "never the canonical's" in out
        # the variant there is UNKNOWN with the same reason; the canonical consequence is labelled as such
        assert f'<b>UNKNOWN</b> <span class="mono f-tissue-reason">{reason}</span>' in out
        assert out.index("p.Trp5Ter") > out.index("on the canonical") and "not the tissue's" in out


@pytest.mark.skipif(
    not (ROOT / "data/reference/chr21.fa.gz").exists()
    or not (ROOT / "data/knowledge/expression/gtex_isoforms_APP.json").exists(),
    reason="needs the local chr21 sequence and the cached GTEx isoform table of APP",
)
def test_api_flow_gives_the_selector_its_tissues_and_the_tissue_isoform_its_geometry():
    from genomeos.web.server import Api

    api = Api(ROOT)
    d = api.flow("APP", "chr21", "25891834 C>A", "Brain_Cerebellum")
    assert "Brain_Cerebellum" in d["tissues_measured"] and len(d["tissues_measured"]) == len(
        set(d["tissues_measured"])
    )
    t = d["tissue"]
    assert t["status"] == "ESTABLISHED" and t["transcript"] == "APP-202" and t["protein_length"] == 695
    assert t["trace"]["transcript"] == "APP-202" and len(t["trace"]["protein"]) == 695
    assert {"exons", "mrna_to_genomic", "cds_start", "cds_end"} <= set(t["trace"])
    v = d["variant_in_tissue"]
    assert v["hgvs_p"] == "p.Gly625Val" and v["canonical"]["hgvs_p"] == "p.Gly700Val"
    assert d["variant"]["hgvs_p"] == "p.Gly700Val"

    x = api.flow("APP", "chr21", "25891834 C>A", "Not_A_Tissue")
    assert x["tissue"]["status"] == "UNKNOWN" and x["tissue"]["reason"] == "tissue_not_measured"
    assert "trace" not in x["tissue"] and x["variant_in_tissue"]["consequence"] == "UNKNOWN"
    assert x["tissues_measured"] == d["tissues_measured"]
