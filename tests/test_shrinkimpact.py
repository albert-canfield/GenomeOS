# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The leaf comparison behind the silent-input-shrink recounts, and the facts it rests on.

Every test here is about the COMPARISON, not about the data: whether a moved leaf is found, whether
a shape difference can be read as agreement, and whether the one-sentence reading can overstate. The
two facts the recount's explanation rests on -- that the CRISPRi library's design carries no chrY
guide, and that `default_gencode` resolves chr21 to the distilled subset under data/results -- are
asserted against TRACKED files only, so nothing here can go inert on a machine without the stores.
"""

from __future__ import annotations

import json
from pathlib import Path

from genomeos.attribution import shrinkimpact as si
from genomeos.genome.annotation import default_gencode
from tests import tracked_paths

ROOT = Path(__file__).resolve().parent.parent
CENSUS = "data/results/closure_cell_inventory.json"
CHR21_SUBSET = "data/results/gencode_v50_chr21_chrM.gff3.gz"


def test_leaves_flattens_nested_dicts_and_lists_to_dotted_paths() -> None:
    flat = si.leaves({"a": {"b": 1}, "c": [{"d": 2}, 3]})
    assert flat == {"a.b": 1, "c[0].d": 2, "c[1]": 3}


def test_an_empty_container_is_itself_a_leaf_so_filling_it_counts_as_movement() -> None:
    # An empty dict flattened to nothing would make `{} -> {"x": 1}` invisible, which is exactly
    # the kind of silent change this lane was opened on.
    assert si.leaves({"a": {}}) == {"a": {}}
    # Filling it changes the SHAPE -- the path `a` disappears and `a.x` appears -- so it is caught
    # as a one-sided leaf rather than as a moved value, and either way it is never agreement.
    report = si.compare({"a": {}}, {"a": {"x": 1}}, label="filled")
    assert report["only_in_committed"] == ["a"]
    assert report["only_in_recomputed"] == ["a.x"]
    assert not si.every_reading_is_no_movement([report])
    # And an empty container that STAYS empty is compared and matches, rather than vanishing.
    assert si.every_reading_is_no_movement([si.compare({"a": {}}, {"a": {}}, label="still empty")])


def test_a_moved_leaf_carries_both_the_old_and_the_new_value() -> None:
    report = si.compare({"n": 35}, {"n": 34}, label="arms.sizes")
    assert report["moved"] == [{"path": "n", "committed": 35, "recomputed": 34}]
    assert "1 of 1 leaves MOVED" in report["reading"]


def test_no_movement_reads_as_no_number_moves_and_names_the_denominator() -> None:
    report = si.compare({"a": 1, "b": 2}, {"a": 1, "b": 2}, label="coverage")
    assert report["moved_count"] == 0
    assert "no number moves" in report["reading"]
    assert "all 2 recomputed leaves" in report["reading"]
    assert si.every_reading_is_no_movement([report])


def test_a_leaf_on_one_side_only_is_never_reported_as_no_number_moves() -> None:
    # The failure this forbids: a recomputation that produced a DIFFERENT SHAPE -- a key the
    # committed result never had, or one it had and the recount lost -- compared clean on the
    # overlap and would have read as agreement.
    report = si.compare({"a": 1}, {"a": 1, "b": 9}, label="coverage")
    assert report["moved_count"] == 0
    assert report["only_in_recomputed"] == ["b"]
    assert "NOT 'no number moves'" in report["reading"]
    assert not si.every_reading_is_no_movement([report])


def test_only_in_committed_is_also_refused_as_agreement() -> None:
    report = si.compare({"a": 1, "b": 2}, {"a": 1}, label="coverage")
    assert report["only_in_committed"] == ["b"]
    assert not si.every_reading_is_no_movement([report])


def test_every_reading_is_no_movement_is_false_on_an_empty_list() -> None:
    # "nothing was compared" must never answer "nothing moved".
    assert not si.every_reading_is_no_movement([])


def test_the_census_declares_22_gencode_members_and_chr21_is_not_one_of_them() -> None:
    """The committed state this lane was briefed on, asserted so the finding cannot drift."""
    tracked_paths.must_be_present(
        ROOT / CENSUS,
        was="skipped when data/results/closure_cell_inventory.json is absent",
    )
    manifest = json.loads((ROOT / CENSUS).read_text())["result_manifest"]
    groups = {entry["path"]: entry for entry in manifest["inputs"] if entry.get("group")}
    gencode = next(e for label, e in groups.items() if label.startswith("GENCODE 50"))
    assert gencode["files"] == 22
    members = [m["path"] for m in gencode["members"]]
    assert not any("chr21" in p for p in members)
    assert not any("chrY" in p for p in members)


def test_chr21s_annotation_is_declared_under_the_sweeps_group_because_it_lives_in_data_results() -> None:
    """The reason the GENCODE group holds 22: a path-prefix grouping, not a missing input."""
    tracked_paths.must_be_present(ROOT / CENSUS, was="skipped without the census result")
    manifest = json.loads((ROOT / CENSUS).read_text())["result_manifest"]
    sweep = next(
        e for e in manifest["inputs"] if e.get("group") and e["path"].startswith("the all-enhancer sweep")
    )
    members = [m["path"] for m in sweep["members"]]
    assert CHR21_SUBSET in members
    assert sweep["files"] == 70  # 3 x 23 sweep files, plus chr21's annotation


def test_default_gencode_resolves_chr21_to_the_distilled_subset_under_data_results() -> None:
    """Why chr21's annotation is read from data/results at all, asserted against a tracked file."""
    tracked_paths.must_be_present(
        ROOT / CHR21_SUBSET,
        was="skipped when the distilled chr21+chrM annotation is absent",
    )
    chosen = default_gencode({"chr21"})
    assert chosen is not None
    assert Path(chosen).name == "gencode_v50_chr21_chrM.gff3.gz"


def test_the_census_itself_states_that_chr21_has_elements_and_chrY_has_none() -> None:
    """chrY's absence from every group is the library's doing, and the census already says so."""
    tracked_paths.must_be_present(ROOT / CENSUS, was="skipped without the census result")
    census = json.loads((ROOT / CENSUS).read_text())
    assert census["coverage"]["per_chromosome"]["chrY"] == 0
    assert census["coverage"]["per_chromosome"]["chr21"] > 0
    assert census["coverage"]["counts"]["chromosomes_with_no_library_guide"] == 1
    assert any("chrY" in line for line in census["result_manifest"]["exclusions"])


THE_FOUR = ("direction_link", "hct116_count", "reptest_answerable", "rescore_cost")


def test_all_four_results_of_2026_10_02_declare_the_same_22_cache_chromosomes() -> None:
    """The committed state the cache-group recount was briefed on, pinned so it cannot drift."""
    declared = {}
    for name in THE_FOUR:
        path = ROOT / f"data/results/{name}.json"
        tracked_paths.must_be_present(path, was=f"skipped when data/results/{name}.json is absent")
        manifest = json.loads(path.read_text())["result_manifest"]
        group = next(
            e
            for e in manifest["inputs"]
            if isinstance(e, dict) and e.get("path") == "per_element_response_cache"
        )
        assert group["files"] == 22
        declared[name] = sorted(Path(m["path"]).name.split(".")[0] for m in group["members"])
    # All four name the SAME 22, and the two they all leave out are chr13 and chrY. If a later
    # change makes them disagree, this lane's single explanation for all four stops being one.
    assert len({tuple(v) for v in declared.values()}) == 1
    missing = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    missing -= set(next(iter(declared.values())))
    assert missing == {"chr13", "chrY"}


def test_the_cache_group_is_built_with_an_exists_filter_in_all_four_entry_scripts() -> None:
    """The mechanism, asserted against the tracked source rather than described in prose."""
    for name in THE_FOUR:
        path = ROOT / f"scripts/{name}.py"
        tracked_paths.must_be_present(path, was=f"skipped when scripts/{name}.py is absent")
        source = path.read_text()
        assert "path.exists()" in source
        assert "per_element_response_cache" in source
