# SPDX-License-Identifier: AGPL-3.0-or-later
"""CRISPRi benchmark parsed, joined to the deletion table, and the predictors compared honestly."""

import io
import json

from genomeos.attribution import crispri

HEADER = (
    "chrom\tchromStart\tchromEnd\tname\tEffectSize\tmeasuredGeneSymbol\tValidConnection\tCellType\t"
    "Regulated\tDataset\tdistanceToTSS\tDHS.RPM\tH3K27ac.RPM\n"
)


def row(chrom, start, end, gene, regulated, distance, dhs=1.0, k27=1.0, valid="TRUE", cell="K562"):
    return (
        f"{chrom}\t{start}\t{end}\t{gene}|x\t-0.1\t{gene}\t{valid}\t{cell}\t"
        f"{'TRUE' if regulated else 'FALSE'}\tD\t{distance}\t{dhs}\t{k27}\n"
    )


def test_parse_drops_invalid_connections():
    text = (
        HEADER
        + row("chr1", 100, 600, "A", True, 5000)
        + row("chr1", 100, 600, "B", False, 900, valid="overlaps potential promoter")
    )
    pairs = crispri.parse(io.StringIO(text))
    assert [(p.gene, p.regulated, p.start, p.distance) for p in pairs] == [("A", True, 100, 5000.0)]


def element(eid, start, end, gene, k562, coding=True):
    pred = {"gene": gene, "log2_fold_change": k562, "tissue": "K562"}
    return {
        "id": eid,
        "start": start,
        "end": end,
        "inferred": None,
        "predicted": pred,
        "predicted_by_cell": {"K562": k562},
        "predicted_coding": pred if coding else None,
        "predicted_coding_by_cell": {"K562": k562} if coding else None,
    }


def test_overlap_and_deletion_features(tmp_path):
    (tmp_path / "chr1.json").write_text(
        json.dumps([element("e2", 900, 1200, "B", 0.3), element("e1", 100, 400, "A", -0.5)])
    )
    table = crispri.DeletionTable(tmp_path)
    assert [e["id"] for e in table.overlapping("chr1", 350, 950)] == ["e2", "e1"]
    assert table.overlapping("chr1", 400, 900) == []  # half-open intervals touch, do not overlap
    assert table.overlapping("chr2", 0, 10) == []
    els = table.overlapping("chr1", 150, 200)
    assert crispri.deletion_for(els, "A", "K562") == (1.0, 0.5)
    assert crispri.deletion_for(els, "Z", "K562") == (0.0, 0.0)
    # a predicted rise is not the activation the screens call
    assert crispri.deletion_for(table.overlapping("chr1", 1000, 1100), "B", "K562") == (1.0, 0.0)


def test_metrics_on_known_rankings():
    assert crispri.average_precision([3, 2, 1], [True, False, True]) == (1 + 2 / 3) / 2
    assert crispri.auroc([3, 2, 1], [True, False, False]) == 1.0
    assert crispri.auroc([1, 1], [True, False]) == 0.5
    assert crispri.average_precision([1, 2], [False, False]) is None


def test_logistic_fit_orders_a_separable_feature():
    x = [[float(i)] for i in range(20)]
    y = [i >= 10 for i in range(20)]
    w = crispri.logistic_fit(x, y)
    assert w[1] > 0
    s = crispri.logistic_score(w, x)
    assert crispri.auroc(s, y) == 1.0


def pairs_for_score(tmp_path):
    """Two genes per element; the near gene is regulated only where the deletion names it."""
    els, text = [], HEADER
    for chrom in ("chr1", "chr2", "chr3"):
        rows = []
        for k in range(12):
            start = 10_000 * (k + 1)
            real = k % 3 == 0
            els_gene = "N" if real else "F"
            rows.append(
                element(f"{chrom}e{k}", start, start + 500, f"{chrom}{els_gene}{k}", -0.4 if real else -0.01)
            )
            text += row(chrom, start, start + 500, f"{chrom}N{k}", real, 5_000, dhs=2.0, k27=2.0)
            text += row(
                chrom, start, start + 500, f"{chrom}F{k}", False, 50_000 + 1_000 * k, dhs=2.0, k27=2.0
            )
        (tmp_path / f"{chrom}.json").write_text(json.dumps(rows))
        els.extend(rows)
    return crispri.parse(io.StringIO(text))


def test_score_runs_the_preregistered_comparison(tmp_path):
    training = pairs_for_score(tmp_path)
    heldout = pairs_for_score(tmp_path)
    for p in heldout[:4]:
        p.cell = "WTC11"
    result = crispri.score(training, heldout, crispri.DeletionTable(tmp_path))
    assert result["coverage"]["training_pairs_on_a_deleted_element"] == len(training)
    assert result["heldout"]["WTC11"]["refused"].startswith("no AlphaGenome line")
    k562 = result["heldout"]["K562"]
    assert k562["models"]["activity + distance + deletion"]["auprc"] == 1.0
    assert k562["passes"] is True and result["verdict"] == "passed"
    # the held-out elements are the training elements here, so nothing is left disjoint
    assert result["heldout_elements_not_in_training"]["K562"]["models"]["distance"]["pairs"] == 0
    calls = result["one_call_per_element_k562"]["0.1"]
    assert calls["all_elements"]["deletion"] == {"calls": 12, "right": 12, "precision": 1.0}
    assert calls["deletion_elements"]["closest gene"]["precision"] == 1.0
