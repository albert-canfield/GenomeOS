# SPDX-License-Identifier: AGPL-3.0-or-later
"""The Cell view: one cell type's active rule set, its graph, and what the program does not say.

The view area G's UI track promised and never built. These tests hold it to the one thing that
matters here: it must never let silence read as a finding. A cell type with no `expresses` list
silences nothing, a program with no rule has no graph rather than an empty one, and the network
runtime's own gate is not the same gate the rule table shows -- each of those has to be stated in
words by the endpoint, not left for the reader to infer from a blank panel.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from genomeos.web.server import Api, ApiError

ROOT = Path(__file__).resolve().parent.parent
# the one program in the project whose rules are gated on cell type and whose cell types say what
# they express; every other cell-type program would give the same answer in all its cell types
CONTEXT = "data/demo/cell_context.bio"


def test_the_same_program_reads_differently_in_two_cell_types() -> None:
    """The claim the view exists to make checkable: one rule set, two cells, two active sets."""
    api = Api(ROOT)

    heart = api.cell(CONTEXT, "Cardiomyocyte2")
    brain = api.cell(CONTEXT, "Neuron2")

    assert heart["counts"] == {"total": 4, "active": 3, "gated_out": 0, "silenced_out": 1}
    assert brain["counts"] == {"total": 4, "active": 1, "gated_out": 1, "silenced_out": 2}
    assert heart["genes_silenced"] == ["SYN1"]
    assert brain["genes_silenced"] == ["MYH6", "NKX2-5"]
    assert heart["cell_type"]["expresses"] == ["NKX2-5", "MYH6"]
    assert heart["cell_type"]["lineage"] == ["Cardiomyocyte", "Cell"]


def test_every_dropped_rule_says_why_it_was_dropped() -> None:
    """A rule missing from the graph with no reason attached is the failure mode of this view."""
    rules = {r["id"]: r for r in Api(ROOT).cell(CONTEXT, "Neuron2")["rules"]}

    assert rules["SYN1->Syn1"]["status"] == "active"
    assert rules["MYH6->Myh6"]["status"] == "silenced"
    assert "`expresses`" in rules["MYH6->Myh6"]["reason"]
    assert rules["Nkx2-5 activates MYH6"]["status"] == "gated out"
    assert "`when`" in rules["Nkx2-5 activates MYH6"]["reason"]
    assert all(r["reason"] == "" for r in rules.values() if r["status"] == "active")


def test_a_silenced_gene_is_a_node_the_graph_keeps_and_marks() -> None:
    """Dropping the node would draw a cell that has no such gene; the cell has it and holds it shut."""
    g = Api(ROOT).cell(CONTEXT, "Cardiomyocyte2")["graph"]

    silenced = [n for n in g["nodes"] if n["silenced"]]
    assert [n["id"] for n in silenced] == ["SYN1"]
    assert [n["kind"] for n in g["nodes"] if n["id"] == "Nkx2-5"] == ["protein"]
    # every rule is an edge, active or not, so the drawing can fade the ones this cell does not run
    assert len(g["edges"]) == 4
    assert sorted({e["status"] for e in g["edges"]}) == ["active", "silenced"]


def test_a_cell_type_with_no_expresses_list_is_said_to_be_silent_not_complete() -> None:
    """`expresses` absent means the program says nothing, which is not the same as "expresses all"."""
    out = Api(ROOT).cell(CONTEXT, "Fibroblast")

    assert out["cell_type"]["expresses"] == []
    assert out["genes_silenced"] == []
    assert any("not a claim that it reads every gene" in n for n in out["notes"])


def test_a_program_that_gates_nothing_says_the_rule_set_is_not_cell_specific() -> None:
    """Otherwise a tab headed by a cell type implies that cell type chose these eleven rules."""
    out = Api(ROOT).cell("data/demo/gastrulation.bio", "Neuron")

    assert out["gated_on_cell_type"] == []
    assert out["counts"]["active"] == out["counts"]["total"] == 11
    assert any("would be identical in every one of its" in n for n in out["notes"])


def test_a_program_with_no_rule_has_no_graph_rather_than_an_empty_one() -> None:
    """The lineage programs declare dozens of cell types and not one rule: that must be said."""
    out = Api(ROOT).cell("data/organisms/human/tissues.bio", "Hepatocyte")

    assert out["counts"]["total"] == 0
    assert out["graph"]["nodes"] == [] and out["graph"]["edges"] == []
    assert out["runnable"] is False
    assert any("states no rule at all" in n for n in out["notes"])


def test_the_run_refuses_a_program_with_no_gene_and_names_where_to_go_instead() -> None:
    with pytest.raises(ApiError, match="no network to integrate"):
        Api(ROOT).cell_run("data/organisms/human/tissues.bio", "Hepatocyte")


def test_an_undeclared_cell_type_is_refused_with_the_declared_ones_listed() -> None:
    with pytest.raises(ApiError, match="Cardiomyocyte2"):
        Api(ROOT).cell(CONTEXT, "Hepatocyte")


def test_a_path_outside_data_or_of_the_wrong_kind_is_refused() -> None:
    api = Api(ROOT)
    with pytest.raises(ApiError):
        api.cell("pyproject.toml")
    with pytest.raises(ApiError, match="not a BioLang program"):
        api.cell("data/demo/demo.fa")


def test_the_run_reports_the_runtime_gate_and_the_cell_gate_apart() -> None:
    """The runtime never applies `expresses`, so a silenced gene transcribes unless it is held down.

    The view would otherwise show one rule active in a neuron and then draw MYH6 rising in it.
    """
    api = Api(ROOT)

    loose = api.cell_run(CONTEXT, "Neuron2", hours=20)
    held = api.cell_run(CONTEXT, "Neuron2", hours=20, silence=True)

    assert loose["active_rules"] == 3 and loose["active_rules_in_cell"] == 1
    assert loose["held_at_zero"] == []
    assert loose["levels"]["MYH6.mRNA"][-1] > 1  # the runtime transcribes a gene the cell silences
    assert held["held_at_zero"] == ["MYH6.mRNA", "NKX2-5.mRNA"]
    assert held["levels"]["MYH6.mRNA"][-1] == 0.0
    assert held["levels"]["SYN1.mRNA"][-1] > 1  # what the neuron does express is untouched


def test_the_index_lists_what_it_skipped_rather_than_only_what_it_found() -> None:
    """A program absent from the list for no stated reason reads as a program with no cell types."""
    out = Api(ROOT).cell_programs()

    paths = {p["path"] for p in out["programs"]}
    assert CONTEXT in paths
    assert not any("compiled" in p for p in paths)  # the generated chromosome programs are not read
    assert out["without_cell_types"] > 0
    assert all(e["path"] and e["error"] for e in out["errors"])
    ctx = next(p for p in out["programs"] if p["path"] == CONTEXT)
    assert out["programs"][0] is ctx  # the one program whose answer changes with the cell comes first
    assert ctx["gated_on_cell_type"] == ["Cardiomyocyte2"]
    assert ctx["default_cell_type"] == "Cardiomyocyte2"
    assert {c["id"]: c["rules_active"] for c in ctx["cell_types"]}["Neuron2"] == 1


def test_the_default_cell_type_is_one_the_program_distinguishes() -> None:
    """Most of these programs import `bio.std.cell_types`, so "the first declared" is the generic Cell.

    Opening on it would show the whole rule set under a heading that names a cell type, which is the
    overstatement this view exists to avoid.
    """
    api = Api(ROOT)

    assert api.cell(CONTEXT)["cell_type"]["id"] == "Cardiomyocyte2"
    assert api.cell(CONTEXT)["cell_types"][0] == "Cell"
    # nothing distinguishes a cell type in a program that gates on none: the first declared stands
    assert api.cell("data/demo/gastrulation.bio")["cell_type"]["id"] == "Cell"
