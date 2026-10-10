# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Guards for the ReSE join: assay-validated silencers beside the derived direction call.

Each guard names the way the join could have been wrong. The synthetic ones run anywhere; the ones on
the committed decision and result read tracked files, which RAISE if missing rather than skip.
"""

from __future__ import annotations

import io
import json
import tokenize
from pathlib import Path

import pytest

from genomeos.attribution import measured
from genomeos.attribution import silenceragree as sa
from genomeos.attribution import silencerdoctrine as sd
from genomeos.attribution import silencerjoin as sj
from genomeos.lang import grammar
from genomeos.lang.parser import parse
from scripts import silencerjoin_decide as decide

ROOT = Path(__file__).resolve().parents[1]
OWN_FILES = (
    "genomeos/attribution/silencerjoin.py",
    "scripts/silencerjoin_decide.py",
    "scripts/silencerjoin_run.py",
)


def _record(chrom: str, start: int, stop: int, cells: str) -> dict:
    return {
        "description": f"ReSE screen-validated silencer GRCh37_chr{chrom}:{start + 1}-{stop + 1}",
        "genomicinfo": [{"chrloc": chrom, "chrstart": start, "chrstop": stop}],
        "summary": f"That subregion was identified as a functional silencer in {cells}.",
    }


def _frag(uid: str, start: int, end: int, lines: list[str], chrom: str = "chr1") -> dict:
    return {"fragment": uid, "chrom": chrom, "start": start, "end": end, "lines": lines}


def _readings(frags: list[dict], spans: list[tuple[int, int, str]]) -> dict[str, dict]:
    pairs = sj.pairs_on(frags, spans)
    sj.label_pairs(pairs)
    where = {eid: ("chr1", s, e) for s, e, eid in spans}
    rows = sj.element_rows(pairs, {p["element"]: where[p["element"]] for p in pairs}, {})
    return {r["id"]: r["by_cell_line"] for r in rows}


# ---------------------------------------------------------------- the rule is the project's, not a new one


def test_the_role_rule_is_the_measured_layers_production_rule_at_its_own_value() -> None:
    """A change to the measured layer's constant must break here, not silently re-decide this join."""
    assert sj.ROLE_RULE == measured.RECIPROCAL_HALF
    assert sj.ROLE_OVERLAP == measured.RECIPROCAL_OVERLAP == 0.5


def test_the_cell_line_map_is_identity_on_the_two_lines_rese_names_and_nothing_else() -> None:
    assert set(sj.CELL_LINE_MAP) == set(sa.RESE_CELL_LINES)
    assert all(k == v for k, v in sj.CELL_LINE_MAP.items())
    assert set(sj.COVERED) < set(sj.SCORED_CELL_LINES)
    assert set(sj.SCORED_CELL_LINES) - set(sj.COVERED) == {"GM12878", "IMR-90"}


def test_a_fragment_is_half_open_from_genomicinfo_and_never_from_the_grch37_description() -> None:
    """genomicinfo gives the 0-based first and last base; the description is GRCh37 and 1-based."""
    rec = _record("7", 2232840, 2233218, "K562 erythroleukemia cells")
    rec["description"] = "ReSE screen-validated silencer GRCh37_chr7:2272476-2272854"
    frags, _ = sj.fragments({"1": rec})
    assert [(f["chrom"], f["start"], f["end"]) for f in frags] == [("chr7", 2232840, 2233219)]
    assert frags[0]["end"] - frags[0]["start"] == 2272854 - 2272476 + 1
    assert frags[0]["lines"] == ["K562"]


def test_a_fragment_naming_both_lines_maps_to_both() -> None:
    frags, _ = sj.fragments({"1": _record("1", 100, 299, "K562 and HepG2 cells")})
    assert frags[0]["lines"] == ["HepG2", "K562"]


# ---------------------------------------------------------------- the join sees no direction


def test_the_join_is_handed_no_direction(tmp_path: Path) -> None:
    """element_spans is the join's only door onto the elements, and it must carry no action field."""
    rows = [
        {"id": "E1", "start": 100, "end": 300, "predicted": {"action": "represses", "gene": "G"}},
        {"id": "E2", "start": 400, "end": 600, "predicted": None},
    ]
    (tmp_path / "chr1.json").write_text(json.dumps(rows))
    spans = sj.element_spans("chr1", tmp_path)
    assert spans == [(100, 300, "E1"), (400, 600, "E2")]
    assert "represses" not in repr(spans)


def test_the_module_never_writes_a_file() -> None:
    """The derived call is read and never written: no write or open call anywhere in the module's code."""
    src = (ROOT / "genomeos/attribution/silencerjoin.py").read_text()
    names = {t.string for t in tokenize.generate_tokens(io.StringIO(src).readline) if t.type == tokenize.NAME}
    assert not names & {"write_text", "write_bytes", "open", "dump", "save_result", "unlink", "rename"}


# ---------------------------------------------------------------- the readings


def test_a_fragment_gives_the_role_in_its_own_line_only() -> None:
    got = _readings([_frag("F1", 100, 300, ["K562"])], [(100, 300, "E1")])["E1"]
    assert got["K562"]["molecular_role"] == "silencer"
    assert got["K562"]["evidence_status"] == "measured/ReSE"
    assert got["K562"]["fragments"] == ["F1"] and "Pang and Snyder" in got["K562"]["evidence"]
    assert got["HepG2"] == sj.unknown(sj.ANOTHER_LINE_ONLY)
    assert got["GM12878"] == got["IMR-90"] == sj.unknown(sj.LINE_NOT_COVERED)


def test_a_fragment_inside_a_much_wider_element_is_below_the_rule() -> None:
    """150 bp of a 350 bp element is 0.43 of the element: reported, never upgraded."""
    got = _readings([_frag("F1", 100, 250, ["K562"])], [(50, 400, "E1")])["E1"]
    assert got["K562"] == sj.unknown(sj.BELOW_RULE)


def test_one_shared_base_is_reported_and_sets_nothing() -> None:
    frags = [_frag("F1", 299, 500, ["HepG2"])]
    spans = [(100, 300, "E1")]
    pairs = sj.pairs_on(frags, spans)
    assert len(pairs) == 1 and pairs[0]["shared_bp"] == 1
    assert _readings(frags, spans)["E1"]["HepG2"] == sj.unknown(sj.BELOW_RULE)


def test_touching_ends_do_not_overlap() -> None:
    assert sj.pairs_on([_frag("F1", 300, 500, ["K562"])], [(100, 300, "E1")]) == []


def test_a_fragment_attaching_to_two_elements_sets_the_role_on_neither() -> None:
    """600 bp spanning two whole 300 bp elements: half of the fragment and all of each element."""
    got = _readings([_frag("F1", 100, 700, ["K562"])], [(100, 400, "E1"), (400, 700, "E2")])
    assert got["E1"]["K562"] == got["E2"]["K562"] == sj.unknown(sj.SEVERAL_ELEMENTS)


def test_a_second_fragment_attaching_uniquely_still_gives_the_role() -> None:
    frags = [_frag("F1", 100, 700, ["K562"]), _frag("F2", 120, 380, ["K562"])]
    got = _readings(frags, [(100, 400, "E1"), (400, 700, "E2")])
    assert got["E1"]["K562"]["molecular_role"] == "silencer" and got["E1"]["K562"]["fragments"] == ["F2"]
    assert got["E2"]["K562"] == sj.unknown(sj.SEVERAL_ELEMENTS)


def test_the_join_window_reaches_an_element_starting_before_the_fragment() -> None:
    spans = [(1000, 1350, "E0"), (1300, 1400, "E1")]
    pairs = sj.pairs_on([_frag("F1", 1340, 1360, ["K562"])], spans)
    assert {p["element"] for p in pairs} == {"E0", "E1"}


def test_no_reading_anywhere_is_a_negative() -> None:
    """UNKNOWN with a reason is the only alternative to the role; absence is never `not a silencer`."""
    frags = [
        _frag("F1", 100, 300, ["K562"]),
        _frag("F2", 1000, 1100, ["HepG2"]),
        _frag("F3", 2000, 2600, ["HepG2"]),
    ]
    spans = [(100, 300, "E1"), (900, 1300, "E2"), (2000, 2300, "E3"), (2300, 2600, "E4")]
    for by_line in _readings(frags, spans).values():
        for r in by_line.values():
            assert r["molecular_role"] in {"silencer", "unknown"}
            if r["molecular_role"] == "unknown":
                assert r["evidence_status"] == "unknown" and r["reason"] in sj.REASONS
    assert "not a negative" in sj.REASONS[sj.LINE_NOT_COVERED]
    assert "not a negative" in sj.REASONS[sj.NOT_LISTED]


# ---------------------------------------------------------------- the readings are legal BioLang


def test_the_reading_values_are_in_the_grammar_so_no_grammar_change_is_needed() -> None:
    assert grammar.AXES["molecular_role"]["silencer"].startswith("a repressive element")
    assert sj.ROLE in grammar.AXES["molecular_role"] and sj.UNKNOWN in grammar.AXES["molecular_role"]
    assert sj.EVIDENCE_STATUS.split("/", 1)[0] in grammar.AXES["evidence_status"]
    assert sj.UNKNOWN in grammar.AXES["evidence_status"]


def test_the_bio_block_parses_and_carries_the_role_and_its_evidence() -> None:
    frags = [_frag("106799845", 100, 300, ["K562"])]
    pairs = sj.pairs_on(frags, [(100, 300, "EH38E1")])
    sj.label_pairs(pairs)
    row = sj.element_rows(pairs, {"EH38E1": ("chr1", 100, 300)}, {})[0]
    assert sj.bio_block(row, "HepG2") is None and sj.bio_block(row, "GM12878") is None
    module = parse(sj.bio_block(row, "K562"))
    el = module.entities["EH38E1_rese_K562"]
    assert el.attrs["ontology"]["molecular_role"] == [["silencer"]]
    assert el.attrs["ontology"]["evidence_status"] == [["measured/ReSE"]]
    assert el.evidence.kind == "experimental"


# ---------------------------------------------------------------- the end to end, on synthetic stores


def _synthetic(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict:
    records = {
        "1": _record("1", 100, 299, "K562 erythroleukemia cells"),
        "2": _record("1", 1000, 1199, "HepG2 hepatocellular carcinoma cells"),
        "3": _record("1", 5000, 5199, "K562 and HepG2 cells"),
        "4": _record("2", 100, 299, "K562 erythroleukemia cells"),
        "5": _record("1", 9000, 9199, "K562 erythroleukemia cells"),
    }
    monkeypatch.setattr(sj.sa, "rese_records", lambda: records)
    (tmp_path / "chr1.json").write_text(
        json.dumps(
            [
                {"id": "E1", "start": 100, "end": 300, "predicted": {"action": "represses"}},
                {"id": "E2", "start": 1000, "end": 1200, "predicted": {"action": "activates"}},
                {"id": "E3", "start": 5000, "end": 5200, "predicted": None},
                {"id": "E4", "start": 7000, "end": 7300, "predicted": {"action": "represses"}},
            ]
        )
    )
    (tmp_path / "chr2.json").write_text(
        json.dumps([{"id": "E5", "start": 50, "end": 500, "predicted": None}])
    )
    return sj.build(tmp_path)


def test_the_counts_per_line_on_a_synthetic_store(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    out = _synthetic(tmp_path, monkeypatch)
    lines = out["by_cell_line"]
    assert out["population"]["elements_scored"] == 5
    assert lines["K562"]["elements_given_the_role"] == 2  # E1, E3
    assert lines["HepG2"]["elements_given_the_role"] == 2  # E2, E3
    for line in ("GM12878", "IMR-90"):
        assert lines[line]["source_covers_the_line"] is False
        assert lines[line]["elements_given_the_role"] == 0
        assert lines[line]["unknown_by_reason"] == {sj.LINE_NOT_COVERED: 5}
    k = lines["K562"]
    assert sum(k["unknown_by_reason"].values()) == k["elements_read_unknown"] == 3
    assert k["unknown_by_reason"][sj.ANOTHER_LINE_ONLY] == 1  # E2
    assert k["unknown_by_reason"][sj.BELOW_RULE] == 1  # E5: 200 of 450 bp
    assert k["fragments"] == {
        "placed_on_one_element": 2,
        "attached_to_several_elements": 0,
        "sharing_bases_below_the_rule_only": 1,
        "on_no_scored_element": 1,
    }
    assert k["fragments_validated_in_the_line"] == sum(k["fragments"].values())


def test_the_coverage_gap_names_every_unread_element_and_why(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    gap = _synthetic(tmp_path, monkeypatch)["coverage_gap"]
    assert set(gap) == set(sj.SCORED_CELL_LINES)
    assert gap["GM12878"]["elements_whose_role_is_not_read"] == 5 and gap["GM12878"]["share_not_read"] == 1.0
    assert "validated_fragments_no_scored_element_carries" not in gap["IMR-90"]
    assert gap["K562"]["elements_whose_role_is_not_read"] == sum(gap["K562"]["by_reason"].values()) == 3
    assert gap["K562"]["validated_fragments_no_scored_element_carries"] == 1


def test_the_derived_call_is_described_not_tested(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    out = _synthetic(tmp_path, monkeypatch)["derived_call_on_the_role_elements"]
    assert out["K562"] == {
        "elements_given_the_role": 2,
        "with_a_named_predicted_target": 1,
        "carrying_the_derived_repressive_direction": 1,
        "carrying_the_derived_activating_direction": 0,
        "with_no_named_predicted_target": 1,
    }
    assert out["given_the_role_in_both_covered_lines"] == 1
    assert "not a test" in out["reading"] and "a5a43b3" in out["reading"]
    flat = json.dumps(out).lower()
    for word in ("control", "interval", "band", "p_value", "permutation"):
        assert f'"{word}' not in flat


def test_the_sensitivity_at_the_role_threshold_equals_the_role_count(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    out = _synthetic(tmp_path, monkeypatch)
    by = out["sensitivity"]["by_threshold"]
    assert out["sensitivity"]["applied"] == "reciprocal_0.5"
    for line in sj.COVERED:
        assert by["reciprocal_0.5"][line] == out["by_cell_line"][line]["elements_given_the_role"]
        assert by["any_shared_base"][line] >= by["reciprocal_0.25"][line] >= by["reciprocal_0.5"][line]
        assert by["reciprocal_0.5"][line] >= by["reciprocal_0.75"][line]


def test_the_layer_never_carries_the_derived_keys_token(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The doctrine pins the results carrying the derived key; this layer must not join that population."""
    text = json.dumps(_synthetic(tmp_path, monkeypatch)) + json.dumps(sj.REPORTED) + json.dumps(sj.rule())
    assert f'"{sd.KEY}"' not in text


# ---------------------------------------------------------------- prose


@pytest.mark.parametrize("rel", OWN_FILES)
def test_no_sentence_gives_a_derived_element_the_repressive_function(rel: str) -> None:
    assert sd.scan_text(rel, ROOT / rel) == []


def test_the_decision_refuses_a_field_only_a_join_could_fill() -> None:
    assert decide.forbidden({"rule": {}, "readings": {}}) == []
    assert decide.forbidden({"a": [{"by_cell_line": {}}]}) == ["/a[0]/by_cell_line"]


def test_what_is_reported_names_the_join_fields_only_as_values() -> None:
    """The decision states what the result will report without carrying a field a join fills."""
    assert decide.forbidden({"reported": sj.REPORTED}) == []
    assert {r["field"] for r in sj.REPORTED} >= {"by_cell_line", "coverage_gap", "elements"}
