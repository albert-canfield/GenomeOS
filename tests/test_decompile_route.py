# SPDX-License-Identifier: AGPL-3.0-or-later
"""The decompile route: a gene read back as a program, with its empty layers named."""

from __future__ import annotations

from pathlib import Path

import pytest

from genomeos.web.server import Api, ApiError

ROOT = Path(".")


def test_a_real_gene_comes_back_with_its_program_and_its_gaps() -> None:
    out = Api(ROOT).decompile("BACH1", "chr21")

    assert out["symbol"] == "BACH1"
    assert out["layers"]["gene"]
    assert out["program"].startswith("# BACH1 on chr21")
    # the unknown list is part of the answer: a decompiler that hid its empty layers would read as
    # completeness, which is the opposite of what this view is for
    assert isinstance(out["unknown"], list)


def test_an_unknown_gene_is_a_404_rather_than_an_empty_program() -> None:
    with pytest.raises(ApiError, match="no gene"):
        Api(ROOT).decompile("NOSUCHGENE", "chr21")


def test_naming_no_gene_asks_for_one() -> None:
    with pytest.raises(ApiError, match="name a gene"):
        Api(ROOT).decompile("", "chr21")


def test_the_symbol_is_normalised_so_a_lowercase_click_works() -> None:
    assert Api(ROOT).decompile("bach1", "chr21")["symbol"] == "BACH1"
