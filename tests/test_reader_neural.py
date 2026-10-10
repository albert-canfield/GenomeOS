# SPDX-License-Identifier: AGPL-3.0-or-later
"""The spinal cord sanity check's registered statistic and decision (scripts/reader_neural.py).

The check is registered before any spinal cord reading exists, so what can be tested here is that the
decision does what the registration says: the exact sign test, the three outcomes and their wording,
the falsifier on either ground, and the frozen gene list.
"""

import importlib.util
import os

import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_SPEC = importlib.util.spec_from_file_location(
    "reader_neural", os.path.join(_ROOT, "scripts", "reader_neural.py")
)
rn = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(rn)


def _counts(**kw):
    base = {"n": 13, "a": 10, "k": 2, "b": 8, "c": 0, "sc_read": 11_000, "sc_coding": 20_000}
    base.update({"k_read": 12_775, "k_coding": 20_094})
    base.update(kw)
    return base


def test_the_sign_test_is_the_exact_one_sided_binomial():
    assert rn.sign_test_p(5, 0) == pytest.approx(1 / 32)
    assert rn.sign_test_p(3, 1) == pytest.approx(5 / 16)
    assert rn.sign_test_p(0, 4) == pytest.approx(1.0)
    assert rn.sign_test_p(0, 0) == 1.0  # no discordant gene: nothing to tell from chance


def test_a_clear_neural_excess_passes_and_the_reading_carries_its_denominators():
    v = rn.verdict(_counts())
    assert v["outcome"] == "PASS"
    assert "10 of 13" in v["reading"] and "2 of 13" in v["reading"]
    assert "11,000 of 20,000" in v["reading"] and "12,775 of 20,094" in v["reading"]
    assert "accessibility is not function" in v["reading"]


def test_the_falsifier_fires_when_spinal_cord_reads_no_more_of_the_list_than_K562():
    v = rn.verdict(_counts(a=6, k=6, b=2, c=2))
    assert v["outcome"] == "FAIL"
    assert "does not read more of the listed genes" in v["reading"]
    assert "not to be used in any comparison" in v["reading"]


def test_the_falsifier_fires_when_the_excess_is_only_depth():
    # more of the list read in spinal cord, but no more than its genome-wide share predicts against K562's
    v = rn.verdict(_counts(a=9, k=8, b=1, c=0, sc_read=19_000, sc_coding=20_000))
    assert v["outcome"] == "FAIL"
    assert "excess over the genome-wide read share" in v["reading"]


def test_the_right_direction_on_too_few_discordant_genes_is_undecided():
    v = rn.verdict(_counts(a=5, k=2, b=3, c=0))
    assert v["outcome"] == "UNDECIDED" and v["p_one_sided"] == pytest.approx(0.125)
    assert "cannot be told from chance" in v["reading"]


def test_the_gene_list_is_the_registered_thirteen_on_their_chromosomes():
    assert [(s, c) for s, c, _ in rn.GENES] == [
        ("PAX7", "chr1"),
        ("IRX3", "chr16"),
        ("DBX1", "chr11"),
        ("DBX2", "chr12"),
        ("PAX6", "chr11"),
        ("NKX6-1", "chr4"),
        ("NKX2-2", "chr20"),
        ("EVX1", "chr7"),
        ("EN1", "chr2"),
        ("VSX2", "chr14"),
        ("MNX1", "chr7"),
        ("ISL1", "chr5"),
        ("SIM1", "chr6"),
    ]
    assert set(rn.WORDING) == {"PASS", "FAIL", "UNDECIDED"}


def test_a_gene_is_read_unless_the_reading_lists_it_silent_or_poised():
    reading = {"silent_genes": ["PAX7", "NKX2-2"], "poised_genes": ["NKX2-2"]}
    assert rn.gene_reads(reading, "PAX6") is True
    assert rn.gene_reads(reading, "PAX7") is False
    assert rn.gene_reads(reading, "NKX2-2") is False  # poised genes are listed in silent_genes too
