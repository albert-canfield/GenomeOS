# SPDX-License-Identifier: AGPL-3.0-or-later
"""CRISPRi benchmark parsed, joined to the deletion table, and the predictors compared honestly."""

import io
import json
import math

import pytest

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


# --- measured contact in place of 1/distance ------------------------------------------------------

HEADER_TSS = HEADER.replace("EffectSize\t", "EffectSize\tstartTSS\t")


def row_tss(chrom, start, end, gene, regulated, tss, dhs=2.0, k27=2.0, cell="K562"):
    """One benchmark row that carries the gene's TSS, as the real tables do."""
    distance = abs(tss - (start + end) // 2)
    return (
        f"{chrom}\t{start}\t{end}\t{gene}|x\t-0.1\t{tss}\t{gene}\tTRUE\t{cell}\t"
        f"{'TRUE' if regulated else 'FALSE'}\tD\t{distance}\t{dhs}\t{k27}\n"
    )


NEAR_TSS, FAR_TSS = 8_000, 60_000  # the two tested genes of every synthetic element
LOOP = 20_000  # the stub calls a TSS further than this the far one


class StubContacts:
    """A contact source without a matrix: on every second element the far gene is the contacted one.

    The stub answers from the geometry alone, never from the labels, and it is built so that which
    gene an element contacts alternates along the chromosome: 1/distance then carries no information
    about the regulated gene, while the measured contact carries all of it. A model that reads only
    distance cannot do better than the base rate, however it fits its coefficient.
    """

    def __init__(self, cells=("K562",), contacted=60.0, other=1.0, missing=()):
        self.matrices = dict.fromkeys(cells, "STUB")
        self.contacted, self.other = contacted, other
        self.missing = set(missing)
        self.asked = 0

    def contact(self, cell, chrom, pos1, pos2):
        self.asked += 1
        if (chrom, pos1) in self.missing:
            return None  # a bin the balancing left undefined
        distance = abs(pos1 - pos2)
        far_is_contacted = (min(pos1, pos2) // 100_000) % 2 == 0
        observed = self.contacted if (distance >= LOOP) == far_is_contacted else self.other
        return {
            "norm": "KR",
            "binsize": 5_000,
            "same_bin": distance < 5_000,
            "raw": observed,
            "observed": observed,
            "expected": 1.0,
            "oe": observed,
        }

    def provenance(self):
        return {"evidence": "stub"}


def pairs_for_contact(tmp_path, cell="K562"):
    """Two genes per element, and the regulated one is the gene the stub says the element contacts."""
    text = HEADER_TSS
    for chrom in ("chr1", "chr2", "chr3"):
        rows = []
        for k in range(12):
            start = 100_000 * (k + 1)
            far_is_contacted = (start // 100_000) % 2 == 0
            near, far = f"{chrom}N{k}", f"{chrom}F{k}"
            rows.append(element(f"{chrom}e{k}", start, start + 500, far if far_is_contacted else near, -0.4))
            text += row_tss(
                chrom, start, start + 500, near, not far_is_contacted, start + NEAR_TSS, cell=cell
            )
            text += row_tss(chrom, start, start + 500, far, far_is_contacted, start + FAR_TSS, cell=cell)
        (tmp_path / f"{chrom}.json").write_text(json.dumps(rows))
    return crispri.parse(io.StringIO(text))


def test_parse_reads_the_tss_column_when_it_is_there():
    pairs = crispri.parse(io.StringIO(HEADER_TSS + row_tss("chr1", 1_000, 1_500, "A", True, 20_000)))
    assert (pairs[0].tss, pairs[0].midpoint, pairs[0].distance) == (20_000, 1_250, 18_750.0)
    # the older synthetic rows have no TSS column, and a pair without one gets no contact
    assert crispri.parse(io.StringIO(HEADER + row("chr1", 100, 600, "A", True, 5_000)))[0].tss is None


def test_annotate_contact_fills_features_and_counts_what_is_missing(tmp_path):
    pairs = pairs_for_contact(tmp_path)
    crispri.annotate(pairs, crispri.DeletionTable(tmp_path))
    source = StubContacts(missing=[("chr1", 100_250)])
    counts = crispri.annotate_contact(pairs, source)
    assert counts["pairs_with_contact"] == len(pairs) - 2
    assert counts["pairs_without_contact"] == 2 and counts["balanced_by_KR"] == len(pairs) - 2
    contacted = next(p for p in pairs if p.gene == "chr2N0")  # the element at 100 kb contacts the near gene
    assert contacted.contact["observed"] == 60.0 and contacted.regulated is True
    assert contacted.features["log_contact"] == pytest.approx(math.log1p(60.0))
    assert contacted.features["activity_x_contact"] == pytest.approx(math.log1p(2.0) + math.log1p(60.0))
    assert next(p for p in pairs if p.gene == "chr2F0").contact["observed"] == 1.0
    dropped = next(p for p in pairs if p.chrom == "chr1" and p.start == 100_000)
    assert dropped.contact is None and dropped.features["log_contact"] == 0.0
    # a cell line with no matrix is never asked
    for p in pairs:
        p.cell = "HepG2"
    assert crispri.annotate_contact(pairs, source)["pairs_without_contact"] == len(pairs)


def test_score_contact_runs_the_preregistered_contact_comparison(tmp_path):
    training = pairs_for_contact(tmp_path)
    heldout = pairs_for_contact(tmp_path)
    for p in heldout[:4]:
        p.cell = "WTC11"  # no matrix for this cell type
    for p in heldout[4:8]:
        p.cell = "HepG2"  # a deletion line, but no matrix either
    result = crispri.score_contact(
        training, heldout, crispri.DeletionTable(tmp_path), StubContacts(cells=("K562", "GM12878"))
    )
    assert result["preregistered"] == crispri.PREREGISTERED_CONTACT
    assert result["coverage"]["training_pairs_scored"] == len(training)
    assert result["heldout"]["WTC11"]["refused"].startswith("no 4DN")
    assert result["heldout"]["HepG2"]["refused"].startswith("no 4DN")
    k562 = result["heldout"]["K562"]
    assert k562["models"]["activity x contact"]["auprc"] == 1.0
    assert k562["models"]["activity + distance"]["auprc"] < 0.8  # 1/distance cannot tell them apart
    assert k562["contact_gain"]["gain"] > 0 and k562["passes"] is True
    assert result["verdict"] == "passed"
    assert set(result["weights"]) == set(crispri.CONTACT_FEATURES)
    assert "contact_gain" in result["training_leave_chromosome_out"]
    assert result["contact_source"] == {"evidence": "stub"}


def test_score_contact_judges_only_the_pairs_with_a_measured_contact(tmp_path):
    training = pairs_for_contact(tmp_path)
    heldout = pairs_for_contact(tmp_path)
    missing = [("chr1", p.midpoint) for p in heldout if p.chrom == "chr1"]
    result = crispri.score_contact(
        training, heldout, crispri.DeletionTable(tmp_path), StubContacts(missing=missing)
    )
    assert result["coverage"]["heldout_pairs_scored"] == len(heldout) - len(missing)
    assert result["heldout"]["K562"]["models"]["contact"]["pairs"] == len(heldout) - len(missing)


def test_contact_coverage_is_reported_by_arm(tmp_path):
    """The denominator discipline: which arm the contact reaches, and whether that subset is easier."""
    pairs = pairs_for_contact(tmp_path)
    crispri.annotate(pairs, crispri.DeletionTable(tmp_path))
    # drop the contact of four regulated pairs, so the arms are covered differently
    missing = [("chr1", p.midpoint) for p in pairs if p.chrom == "chr1" and p.regulated][:4]
    crispri.annotate_contact(pairs, StubContacts(missing=missing))
    report = crispri.contact_coverage_by_arm(pairs)
    regulated = report["arms"]["regulated"]
    assert regulated["with_a_measured_contact"] == regulated["pairs_on_a_deleted_element"] - 4
    # a dropped bin is the element's bin, so the element's other gene loses its contact as well
    assert report["arms"]["not regulated"]["with_a_measured_contact"] == (
        report["arms"]["not regulated"]["pairs_on_a_deleted_element"] - 4
    )
    assert regulated["fraction"] < 1.0
    # both baselines are read twice, on every covered pair and on the subset that has a contact
    assert set(report["baselines"]) == {"on a deleted element", "and with a measured contact"}
    for reading in report["baselines"].values():
        assert reading["distance"]["auprc"] is not None
    assert report["baselines"]["and with a measured contact"]["distance"]["pairs"] == len(pairs) - 8


def test_a_pair_without_a_contact_is_counted_with_its_reason(tmp_path):
    pairs = pairs_for_contact(tmp_path)
    crispri.annotate(pairs, crispri.DeletionTable(tmp_path))
    for p in pairs[:2]:
        p.tss = None
    for p in pairs[2:4]:
        p.cell = "HepG2"
    counts = crispri.annotate_contact(pairs, StubContacts(missing=[("chr1", pairs[4].midpoint)]))
    assert counts["without_contact_no_tss_column_for_this_pair"] == 2
    assert counts["without_contact_no_matrix_for_this_cell_line"] == 2
    assert counts["without_contact_no_balanced_bin_in_the_matrix"] >= 1
    assert counts["pairs_without_contact"] == 4 + counts["without_contact_no_balanced_bin_in_the_matrix"]
