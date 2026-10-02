# SPDX-License-Identifier: AGPL-3.0-or-later
"""Chromatin evidence for the cell a compiled rule names, or the label saying why there is none.

R1 of the external review (2026-09-28) made every compiled rule carry its cell as an executable
`when:`. That cell is asserted. These tests pin the reading that puts reader v1's own chromatin
beside it, under the conditions it was authorised on:

(a) the states are named by the measurement, so no name may imply support or contradiction;
(b) not-open is not closed, and a never-looked state is a different output from a measured one;
(c) the cell is mapped to a reader biosample by ontology term through a registered table, never by
    comparing names, and an ambiguous mapping is `not_assessable`;
(d) reader v1's registered call is used unchanged and no new cut-off is introduced;
(e) the state survives compile -> parse -> BioIR JSON -> parse, in the style of S1 and R1, and it is
    added beside the rule: no verdict moves and nothing is deleted.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from genomeos.attribution import compile as cp
from genomeos.attribution import context_evidence as ce
from genomeos.attribution import measured
from genomeos.genome import reader
from genomeos.ir import Module
from genomeos.lang import parse

CHROM = "chrT"

# ---- (a) the names carry a measurement, never a verdict ---------------------------------------


def test_no_state_or_reason_name_implies_support_or_contradiction():
    for name in (*ce.STATES, *ce.NOT_ASSESSABLE_REASONS):
        for stem in ce.FORBIDDEN_IN_A_STATE_NAME:
            assert stem not in name.lower(), (
                f"{name!r} contains {stem!r}: a state is a measurement, and this reading issues no "
                "verdict on any rule"
            )


def test_the_states_are_the_three_registered_ones_and_nothing_else():
    assert ce.STATES == ("open_in_reader", "not_open_in_reader", "not_assessable")
    assert ce.NOT_ASSESSABLE_REASONS == ("no_reader_for_this_cell", "region_outside_measured_span")
    with pytest.raises(ValueError):
        ce.value("supported")


def test_every_state_and_reason_has_a_reading_and_the_weakness_is_stated():
    for name in (*ce.STATES, *ce.NOT_ASSESSABLE_REASONS):
        assert ce.READING[name]
    assert "not validation" in ce.NOT_VALIDATION
    assert "trained on ENCODE" in ce.NOT_VALIDATION
    assert "not closed" in ce.NOT_CLOSED
    reg = ce.registration()
    assert reg["not_validation"] == ce.NOT_VALIDATION
    assert reg["alphagenome_requests"] == 0


# ---- (b) not-open is not closed, and the two not-assessable reasons are different outputs ------


def peak_file(tmp_path: Path, cell: str, chrom: str, peaks: list[tuple[int, int, float]]) -> None:
    path = tmp_path / reader.peaks_path(cell, chrom).name
    with gzip.open(path, "wt") as fh:
        fh.write("# synthetic\n")
        for s, e, v in peaks:
            fh.write(f"{s}\t{e}\t{v:.2f}\n")


def test_open_not_open_and_outside_the_span_are_three_different_outputs(tmp_path):
    peak_file(tmp_path, "K562", CHROM, [(1000, 1200, 5.0), (5000, 5200, 4.0)])
    readers = ce.Readers(results_dir=tmp_path)
    assert ce.state_for("K562", CHROM, 1100, 1150, readers) == "open_in_reader, K562"
    # inside the file's span, no peak over the locus: a negative call at the reader's own threshold
    assert ce.state_for("K562", CHROM, 3000, 3100, readers) == "not_open_in_reader, K562"
    # beyond the outermost peak no read was placed, so absence is not a negative call
    assert ce.state_for("K562", CHROM, 90_000, 90_100, readers) == (
        "not_assessable, region_outside_measured_span, K562"
    )


def test_a_cell_with_no_reader_and_a_cell_whose_span_misses_the_locus_differ(tmp_path):
    peak_file(tmp_path, "K562", CHROM, [(1000, 1200, 5.0)])
    readers = ce.Readers(results_dir=tmp_path)
    never_looked = ce.state_for("placenta", CHROM, 1100, 1150, readers)
    looked_elsewhere = ce.state_for("K562", CHROM, 80_000, 80_100, readers)
    assert never_looked == "not_assessable, no_reader_for_this_cell"
    assert looked_elsewhere.startswith("not_assessable, region_outside_measured_span")
    assert never_looked != looked_elsewhere


def test_a_chromosome_the_biosample_has_no_peak_file_for_is_outside_its_measured_span(tmp_path):
    readers = ce.Readers(results_dir=tmp_path)  # no file at all
    assert ce.state_for("K562", CHROM, 1000, 1100, readers) == (
        "not_assessable, region_outside_measured_span, K562"
    )


def test_an_unrecorded_context_is_reported_as_no_reader_and_never_as_not_open(tmp_path):
    """R1 writes `cell_type = unknown` where no cell was recorded. That is a never-looked state."""
    peak_file(tmp_path, "K562", CHROM, [(1000, 1200, 5.0)])
    readers = ce.Readers(results_dir=tmp_path)
    assert ce.state_for(measured.CONTEXT_UNKNOWN, CHROM, 1100, 1150, readers) == (
        "not_assessable, no_reader_for_this_cell"
    )


# ---- (c) the mapping is by ontology term, through a registered table ---------------------------


def test_every_reader_biosample_carries_the_term_the_track_metadata_on_disk_gives_it():
    """The hand-written table cannot drift from the authority it claims to copy."""
    table = ce.mapping()
    for cell, term in ce.READER_TERMS.items():
        assert table.term_of(ce.ident(cell)) == term, (
            f"{cell!r} is registered as {term} but the track metadata gives {table.term_of(ce.ident(cell))}"
        )


def test_the_registered_terms_are_well_formed_distinct_and_in_the_three_ontologies():
    assert len(ce.READER_TERMS) == 13
    assert len(set(ce.READER_TERMS.values())) == 13
    for term in ce.READER_TERMS.values():
        assert ce.TERM_PATTERN.match(term)
        assert term.split(":", 1)[0] in ce.READER_PREFIXES
    assert set(ce.term_to_reader().values()) == set(ce.READER_TERMS)


def test_two_different_names_for_one_term_both_reach_the_reader():
    """`Ovary` is a GTEx tissue name and `ovary` an ENCODE biosample name; the term is the same one."""
    assert ce.reader_for("Ovary")[0] == "ovary"
    assert ce.reader_for("ovary")[0] == "ovary"
    assert ce.reader_for("female_gonad")[0] == "ovary"
    assert ce.reader_for("H1_hESC")[0] == "H1"


def test_a_name_that_merely_contains_a_reader_s_name_does_not_map_to_it():
    """What a string comparison would get wrong: a different cell type with an overlapping name."""
    for label in ("foreskin_keratinocyte", "hair_follicular_keratinocyte", "regular_cardiac_myocyte"):
        biosample, term, why = ce.reader_for(label)
        assert biosample is None, f"{label} must not be read as a reader biosample"
        assert term and term not in ce.READER_TERMS.values()
        assert why


def test_a_cell_whose_label_resolves_to_more_than_one_term_is_not_assessable(tmp_path):
    csv_path = tmp_path / "metadata.csv"
    csv_path.write_text(
        "name,ontology_curie,biosample_name,gtex_tissue\n"
        "a,EFO:0002067,twoways,\n"
        "b,CL:0000127,twoways,\n"
        "c,EFO:0002067,K562,\n"
    )
    table = ce.mapping(csv_path)
    assert table.ambiguous["twoways"] == ("CL:0000127", "EFO:0002067")
    biosample, term, why = ce.reader_for("twoways", table)
    assert (biosample, term) == (None, None)
    assert "more than one ontology term" in why
    peak_file(tmp_path, "K562", CHROM, [(1000, 1200, 5.0)])
    readers = ce.Readers(results_dir=tmp_path)
    assert ce.state_for("twoways", CHROM, 1100, 1150, readers, table) == (
        "not_assessable, no_reader_for_this_cell"
    )
    ce.mapping(ce.TRACK_METADATA)  # the module-level cache goes back to the real table


def test_the_mapping_refuses_to_guess_when_the_table_is_not_on_this_machine(tmp_path):
    with pytest.raises(FileNotFoundError, match="never guessed from names"):
        ce.mapping(tmp_path / "absent.csv")
    ce.mapping(ce.TRACK_METADATA)


# ---- (d) reader v1's own call, no new cut-off --------------------------------------------------


def test_the_openness_call_is_reader_v1_s_registered_evidence_unchanged():
    assert reader.EVIDENCE in ce.OPENNESS_CALL
    assert "No new cut-off is introduced" in ce.OPENNESS_CALL
    assert ce.registration()["openness_call"] == ce.OPENNESS_CALL


def test_openness_is_peak_overlap_and_nothing_here_thresholds_a_signal_value(tmp_path):
    """A peak with the lowest signal in the file still opens the locus: no signal is thresholded."""
    peak_file(tmp_path, "K562", CHROM, [(1000, 1200, 0.01), (5000, 5200, 900.0)])
    readers = ce.Readers(results_dir=tmp_path)
    assert ce.state_for("K562", CHROM, 1100, 1150, readers) == "open_in_reader, K562"
    source = Path("genomeos/attribution/context_evidence.py").read_text()
    for forbidden in ("signal >", "signal_value", "> 0.0", ">= 0.0"):
        assert forbidden not in source


def test_a_locus_touching_a_peak_edge_counts_as_overlapping_as_the_reader_s_index_does(tmp_path):
    # a far peak so the whole region is inside the file's measured span and the edge is what is tested
    peak_file(tmp_path, "K562", CHROM, [(1000, 1200, 5.0), (90_000, 90_200, 5.0)])
    readers = ce.Readers(results_dir=tmp_path)
    assert ce.state_for("K562", CHROM, 1199, 1400, readers) == "open_in_reader, K562"
    assert ce.state_for("K562", CHROM, 1200, 1400, readers) == "not_open_in_reader, K562"


# ---- (e) the state survives the round trip, beside the rule ------------------------------------


def pair(gene, effect, cell):
    return measured.CrispriPair(
        chrom=CHROM,
        start=1000,
        end=1300,
        gene=gene,
        cell=cell,
        dataset="Syn2026",
        reference="synthetic",
        regulated=effect < 0,
        significant=True,
        effect_size=effect,
        p_adjusted=0.01,
        split=measured.TRAINING,
        power_at_effect_size_20=0.95,
    )


def row_of(*pairs):
    m = measured.Layer(chrom=CHROM, crispri=list(pairs)).for_element(1000, 1300)
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


def program(row, context_state=None, genes=("AAA",)):
    ident_of = {g: cp.ident(g) for g in genes}
    block = cp._measured_blocks(CHROM, row, {}, ident_of, None, context_state)
    stubs = [
        f'gene {cp.ident(g)} {{ symbol: {g}; evidence: curated "GENCODE v50" symbol only, a stub }}'
        for g in genes
    ]
    return "\n".join([f"module human.noncoding.{CHROM}", "", *stubs, *block]) + "\n"


def round_trip(text: str) -> Module:
    """compile -> parse -> BioIR JSON -> parse, as S1 and R1 check it."""
    return Module.from_dict(json.loads(json.dumps(parse(text).to_dict())))


def test_the_state_survives_compile_parse_bioir_json_and_parse(tmp_path):
    peak_file(tmp_path, "K562", CHROM, [(1000, 1200, 5.0)])
    # HepG2 was read over this element's region - its span brackets it - and called no peak on it
    peak_file(tmp_path, "HepG2", CHROM, [(100, 200, 5.0), (8000, 8200, 5.0)])
    readers = ce.Readers(results_dir=tmp_path)

    def state(label, start, end):
        return ce.state_for(label, CHROM, start, end, readers)

    row = row_of(pair("AAA", -30.0, "K562"), pair("AAA", -25.0, "HepG2"))
    text = program(row, state)
    for module in (parse(text), round_trip(text)):
        got = {r.when["cell_type"]: r.context_evidence for r in module.rules}
        assert got == {
            "K562": "open_in_reader, K562",
            "HepG2": "not_open_in_reader, HepG2",
        }, got


def test_the_written_state_parses_back_to_its_state_and_its_qualifiers():
    for state in ce.STATES:
        assert ce.parse_value(ce.value(state, "K562")) == (state, ("K562",))
    assert ce.parse_value("not_assessable, region_outside_measured_span, K562") == (
        "not_assessable",
        ("region_outside_measured_span", "K562"),
    )
    with pytest.raises(ValueError):
        ce.parse_value("open in K562")


def test_the_state_is_added_beside_the_rule_and_deletes_nothing(tmp_path):
    """The only difference the reading makes to a program is the field it adds."""
    peak_file(tmp_path, "K562", CHROM, [(1000, 1200, 5.0)])
    readers = ce.Readers(results_dir=tmp_path)
    row = row_of(pair("AAA", -30.0, "K562"))
    without = program(row)
    with_state = program(row, lambda label, s, e: ce.state_for(label, CHROM, s, e, readers))
    assert "context_evidence" not in without
    assert with_state.replace("context_evidence: open_in_reader, K562; ", "") == without


def test_the_state_gates_nothing_a_rule_not_detected_open_still_applies_in_its_cell(tmp_path):
    peak_file(tmp_path, "HepG2", CHROM, [(100, 200, 5.0), (8000, 8200, 5.0)])
    readers = ce.Readers(results_dir=tmp_path)
    row = row_of(pair("AAA", -30.0, "HepG2"))
    text = program(row, lambda label, s, e: ce.state_for(label, CHROM, s, e, readers))
    (rule,) = parse(text).rules
    assert rule.context_evidence == "not_open_in_reader, HepG2"
    assert rule.applies({"cell_type": "HepG2"})
    assert not rule.applies({"cell_type": "K562"})


def test_the_grammar_names_the_key_so_the_written_doc_does_too():
    from genomeos.lang import grammar

    assert "context_evidence" in grammar.BLOCKS["rule"]["props"]
    assert "context_evidence" in Path("docs/BIOLANG-GRAMMAR.md").read_text()


# ---- the denominator: every compiled rule ------------------------------------------------------


@pytest.mark.skipif(
    not (Path("data/results/budget_chr21.json").exists() and ce.TRACK_METADATA.exists()),
    reason="chr21's budget or the track metadata is not on this machine",
)
def test_every_rule_of_a_real_compiled_chromosome_carries_a_state():
    text = cp.compile_chromosome("chr21")
    rules = [ln for ln in text.splitlines() if ln.startswith("rule ")]
    assert rules
    for line in rules:
        assert "context_evidence: " in line, line[:120]
        value = line.split("context_evidence: ", 1)[1].split(";", 1)[0]
        state, _rest = ce.parse_value(value)
        assert state in ce.STATES
