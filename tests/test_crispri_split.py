# SPDX-License-Identifier: AGPL-3.0-or-later
"""The CRISPRi benchmark's split and power, carried through the measured layer and the compiler.

The split audit of 2026-09-28 (docs/ATTRIBUTION.md, data/results/crispri_split_audit.json) found
that no held-out figure had read held-out pairs through `load_crispri` or the compiled program, but
that 43 held-out-only links sat in the compiled programs unmarked. These tests pin the fix: every
pair knows its file, a negative knows its power, and a held-out link never reaches `targets:` and
never reaches a rule without the mark.
"""

from __future__ import annotations

import gzip
import io

import pytest

from genomeos.attribution import compile as cp
from genomeos.attribution import measured
from genomeos.lang import parse

HEADER = (
    "chrom\tchromStart\tchromEnd\tname\tEffectSize\tmeasuredGeneSymbol\tSignificant\tpValueAdjusted\t"
    "PowerAtEffectSize10\tPowerAtEffectSize15\tPowerAtEffectSize20\tPowerAtEffectSize25\t"
    "PowerAtEffectSize50\tValidConnection\tCellType\tReference\tRegulated\tDataset"
)


def row(start, end, gene, regulated, effect, cell, power="0.9", valid="TRUE"):
    reg = "TRUE" if regulated else "FALSE"
    return (
        f"chrT\t{start}\t{end}\tx\t{effect}\t{gene}\t{reg}\t0.01\t"
        f"{power}\t{power}\t{power}\t{power}\t{power}\t{valid}\t{cell}\tsyn\t{reg}\tSyn2026"
    )


TRAIN_ROWS = [
    row(1000, 1300, "AAA", True, -0.2, "K562"),
    row(1000, 1300, "BBB", False, 0.01, "K562", power="0.05"),
    row(5000, 5300, "CCC", False, 0.0, "K562", valid="FALSE"),
]
HELD_ROWS = [
    row(1000, 1300, "AAA", True, -0.8, "WTC11"),  # stronger than the training pair: must not set the rule
    row(1000, 1300, "DDD", True, -0.5, "GM12878"),  # held-out only
    row(1000, 1300, "EEE", False, 0.0, "GM12878", power="NA"),
]


def table(rows):
    return io.StringIO("\n".join([HEADER, *rows]) + "\n")


@pytest.fixture
def knowledge(tmp_path):
    for name, rows in zip(measured.CRISPRI_FILES, (TRAIN_ROWS, HELD_ROWS), strict=True):
        with gzip.open(tmp_path / name, "wt") as fh:
            fh.write(table(rows).getvalue())
    return tmp_path


def test_the_power_columns_are_the_ones_both_headers_carry():
    """Read from the real headers when the tables are cached, so the names are not a guess."""
    for name in measured.CRISPRI_FILES:
        path = measured.CRISPRI_KNOWLEDGE / name
        if not path.exists():
            pytest.skip("CRISPRi benchmark tables not cached on this machine")
        with gzip.open(path, "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
        assert [c for c in header if c.startswith("PowerAtEffectSize")] == list(measured.POWER_COLUMNS)


def test_a_pair_carries_its_split_and_power_and_an_empty_power_is_none_not_zero():
    got, invalid = measured.parse_crispri(table(TRAIN_ROWS + HELD_ROWS[2:]), split=measured.HELDOUT)
    assert invalid == 1
    by_gene = {p.gene: p for p in got}
    assert all(p.split == measured.HELDOUT for p in got)
    assert by_gene["AAA"].power_at_effect_size_20 == 0.9
    assert by_gene["BBB"].power_at_effect_size_10 == 0.05  # a blind negative, now told apart
    assert by_gene["EEE"].power_at_effect_size_50 is None


def test_parse_without_a_split_keeps_the_old_default():
    got, _ = measured.parse_crispri(table(TRAIN_ROWS))
    assert {p.split for p in got} == {measured.TRAINING}


def test_load_crispri_selects_by_split_and_its_default_is_still_both(knowledge):
    everything, bad = measured.load_crispri("chrT", knowledge)
    assert bad == 1 and len(everything) == 5  # the audit found no exposed caller, so "all" stays
    train, _ = measured.load_crispri("chrT", knowledge, split=measured.TRAINING)
    held, _ = measured.load_crispri("chrT", knowledge, split=measured.HELDOUT)
    assert {p.split for p in train} == {measured.TRAINING} and len(train) == 2
    assert {p.split for p in held} == {measured.HELDOUT} and len(held) == 3
    assert sorted(train + held, key=repr) == sorted(everything, key=repr)
    with pytest.raises(ValueError):
        measured.load_crispri("chrT", knowledge, split="test")


def compiled_row(knowledge):
    pairs, _ = measured.load_crispri("chrT", knowledge)
    layer = measured.Layer(chrom="chrT", crispri=pairs)
    m = layer.for_element(1000, 1300)
    return {
        "id": "E1",
        "start": 1000,
        "end": 1300,
        "class": "enhancer",
        "measured": m,
        "agreement": measured.agreement("AAA", m),
        "predicted_gene": "AAA",
        "confidence": measured.confidence_of(m),
    }


def test_the_layer_keeps_every_link_and_names_the_training_ones(knowledge):
    c = compiled_row(knowledge)["measured"]["crispri"]
    assert c["genes_regulated"] == ["AAA", "DDD"]  # nothing measured is dropped from the layer
    assert c["genes_regulated_training"] == ["AAA"]
    assert {p["split"] for p in c["pairs"] if p["gene"] == "DDD"} == {measured.HELDOUT}


# Since R1 (2026-09-28, docs/ATTRIBUTION.md) a link is one (element, gene, cell): AAA is two rules,
# K562 (training, 0.2) and WTC11 (held-out, 0.8, marked), and the invariant - a held-out pair never
# sets a training rule's number - holds within each cell (also pinned in tests/test_rule_context.py).
def test_a_held_out_measurement_never_sets_a_compiled_number(knowledge):
    assert measured.rule_links(compiled_row(knowledge)) == [
        ("AAA", "activates", 0.2, "K562", measured.TRAINING),  # not WTC11's -0.8
        ("AAA", "activates", 0.8, "WTC11", measured.HELDOUT),  # its own cell's rule, never K562's
        ("DDD", "activates", 0.5, "GM12878", measured.HELDOUT),
    ]
    assert measured.rule_lines(compiled_row(knowledge)) == [
        ("AAA", "activates", 0.2, "K562"),
        ("AAA", "activates", 0.8, "WTC11"),
        ("DDD", "activates", 0.5, "GM12878"),
    ]


# R1: one rule per (element, gene, cell). A held-out link still never reaches `targets:`, and every
# rule resting only on held-out pairs carries the mark, in whichever cell it is gated on.
def test_the_compiled_block_keeps_held_out_links_out_of_targets_and_marks_their_rules(knowledge):
    r = compiled_row(knowledge)
    ident_of = {g: cp.ident(g) for g in ("AAA", "DDD")}
    block = cp._measured_blocks("chrT", r, {}, ident_of)
    text = "\n".join(
        [
            "module t",
            *(
                f'gene {i} {{ symbol: {g}; evidence: curated "x"; confidence: 0.9 }}'
                for g, i in ident_of.items()
            ),
            *block,
        ]
    )
    module = parse(text)
    element = module.entities["E1_measured"]
    targets = [line for line in block if line.strip().startswith("targets:")]
    assert targets == [f"  targets: {ident_of['AAA']}"]  # DDD is held-out only
    own = [r_ for r_ in module.rules if r_.source == "E1_measured"]
    rules = {(r_.target, r_.when.get("cell_type")): r_ for r_ in own}
    assert len(own) == len(rules) == 3  # one rule per gene and cell, none dropped for being held out
    aaa_k562, aaa_wtc11, ddd = (
        (ident_of["AAA"], "K562"),
        (ident_of["AAA"], "WTC11"),
        (ident_of["DDD"], "GM12878"),
    )
    assert set(rules) == {aaa_k562, aaa_wtc11, ddd}
    assert measured.is_heldout(rules[ddd].evidence.source)
    assert measured.is_heldout(rules[aaa_wtc11].evidence.source)
    assert not measured.is_heldout(rules[aaa_k562].evidence.source)
    assert rules[aaa_k562].strength == 0.2  # the training rule's number is its own, not WTC11's 0.8
    assert element.evidence.kind == "experimental"
    assert "held-out split" in [line for line in block if "basis:" in line][0]


# --- R5: provenance on every record, a guard for development readers, the overlap check -------------
def test_every_record_names_its_study_assay_file_and_partition_and_cannot_be_changed(knowledge):
    pairs, _ = measured.load_crispri("chrT", knowledge)
    for p in pairs:
        assert p.study == "Syn2026" and p.assay == "CRISPRi"
        assert p.source_file in measured.CRISPRI_FILES
        assert p.partition == measured.CRISPRI_SPLIT_OF[p.source_file]
    with pytest.raises(AttributeError):
        pairs[0].split = measured.TRAINING  # frozen: a record's partition is fixed at parse time


def test_a_development_reader_refuses_a_held_out_pair_instead_of_dropping_it(knowledge):
    train, _ = measured.load_crispri("chrT", knowledge, split=measured.TRAINING)
    assert measured.development_only(train) == train
    everything, _ = measured.load_crispri("chrT", knowledge)
    with pytest.raises(ValueError, match="held-out"):
        measured.development_only(everything)


def p_(start, end, gene, split, cell="K562"):
    return measured.CrispriPair(
        chrom="chrT",
        start=start,
        end=end,
        gene=gene,
        cell=cell,
        dataset="s",
        reference="r",
        regulated=False,
        significant=False,
        effect_size=0.0,
        p_adjusted=1.0,
        split=split,
    )


def test_the_overlap_check_sees_related_intervals_and_not_only_identical_pairs():
    t, h = measured.TRAINING, measured.HELDOUT
    pairs = [
        p_(1000, 2000, "A", t),
        p_(1000, 2000, "A", h),  # identical pair
        p_(1010, 2000, "A", h, cell="WTC11"),  # near-identical, same gene
        p_(1800, 2800, "A", h),  # overlapping, same gene
        p_(1000, 2000, "B", h),  # the same element on another gene
        p_(1900, 3000, "C", h),  # overlapping, another gene
        p_(9000, 9500, "A", h),  # independent
    ]
    got = measured.split_overlap(pairs)
    assert got["counts"] == {
        "identical_pair": 1,
        "near_identical_same_gene": 1,
        "overlapping_same_gene": 1,
        "near_identical_other_gene": 1,
        "overlapping_other_gene": 1,
        "independent": 1,
    }
    assert got["heldout_related_to_training"] == 5
    assert got["by_heldout_cell"]["WTC11"]["near_identical_same_gene"] == 1


def test_the_committed_overlap_counts_are_what_the_benchmark_files_give():
    """The audit's numbers are recomputed from the files whenever they are on this machine."""
    from pathlib import Path

    from genomeos.results import load_result

    if not all((measured.CRISPRI_KNOWLEDGE / n).exists() for n in measured.CRISPRI_FILES):
        pytest.skip("CRISPRi benchmark tables not cached on this machine")
    audit = load_result("crispri_split_audit")
    assert audit, Path("data/results/crispri_split_audit.json")
    pairs, _ = measured.load_crispri()
    assert measured.split_overlap(pairs) == audit["overlap_between_partitions"]


# --- R2: a CRISPRi outcome is one of five things, and only one of them rejects -----------------------
def q_(gene, significant, effect, power=0.95):
    return measured.CrispriPair(
        chrom="chrT",
        start=1000,
        end=1300,
        gene=gene,
        cell="K562",
        dataset="s",
        reference="r",
        regulated=significant and effect < 0,
        significant=significant,
        effect_size=effect,
        p_adjusted=0.01 if significant else 0.5,
        power_at_effect_size_20=power,
    )


def test_each_pair_is_one_of_five_outcomes():
    assert q_("A", True, -0.3).outcome == measured.DECREASE
    assert q_("A", True, 0.3).outcome == measured.INCREASE
    assert q_("A", False, 0.01, power=0.95).outcome == measured.NULL_INFORMATIVE
    assert q_("A", False, 0.01, power=0.2).outcome == measured.NULL_INCONCLUSIVE
    assert q_("A", False, 0.01, power=None).outcome == measured.NULL_INCONCLUSIVE
    assert q_("A", False, float("nan")).outcome == measured.MISSING
    got, _ = measured.parse_crispri(table([row(1, 2, "Z", False, "NA", "K562")]))
    assert got[0].outcome == measured.MISSING  # an empty effect is not read as zero


def element(*pairs, predicted):
    m = measured.Layer(chrom="chrT", crispri=list(pairs)).for_element(1000, 1300)
    return {"measured": m, "agreement": measured.agreement(predicted, m), "predicted_gene": predicted}


def test_a_significant_increase_never_becomes_no_effect_and_raises_no_rule():
    r = element(q_("UP", True, 0.4), q_("NUL", False, 0.0), predicted="UP")
    c = r["measured"]["crispri"]
    assert c["genes_increased"] == ["UP"]
    assert "UP" not in c["genes_not_regulated"] and c["genes_not_regulated"] == ["NUL"]
    assert r["agreement"]["crispri"] == measured.INCREASED
    assert r["agreement"]["assays_disagreeing"] == 0 and r["agreement"]["assays_agreeing"] == 0
    assert measured.rule_links(r) == []  # an increase is not a silencer, and nothing says it is
    assert "measured no effect on NUL" in measured.basis_text(r)
    assert "no effect on UP" not in measured.basis_text(r)


def test_only_a_well_powered_null_rejects_the_compiled_claim():
    strong = element(q_("G", False, 0.0, power=0.9), predicted="G")
    weak = element(q_("G", False, 0.0, power=0.3), predicted="G")
    assert strong["agreement"]["crispri"] == measured.DISAGREES
    assert weak["agreement"]["crispri"] == measured.UNDERPOWERED
    assert weak["agreement"]["assays_disagreeing"] == 0
    assert weak["measured"]["crispri"]["genes_no_effect_underpowered"] == ["G"]
    assert "underpowered, on G" in measured.basis_text(weak)


def test_a_gene_takes_its_strongest_outcome_across_pairs():
    r = element(q_("G", False, 0.0, power=0.9), q_("G", True, 0.2), predicted="G")
    assert r["measured"]["crispri"]["genes_increased"] == ["G"]
    assert r["measured"]["crispri"]["genes_no_effect_well_powered"] == []
    r = element(q_("G", True, 0.2), q_("G", True, -0.2), predicted="G")
    assert r["agreement"]["crispri"] == measured.AGREES
