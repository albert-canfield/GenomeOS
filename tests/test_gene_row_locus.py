# SPDX-License-Identifier: AGPL-3.0-or-later
"""Can a cached gene row be tied to a locus? The reading, on constructed caches.

Every test here builds its own archive, so none depends on the project's data being present. The
ones that matter are the two that distinguish a property of a SYMBOL from a property of a LOCUS:
`test_pinning_distinguishes_two_loci_of_the_same_name` and
`test_coding_name_is_a_property_of_the_symbol_and_pins_nothing`.
"""

from __future__ import annotations

import gzip
import json

import pytest

from genomeos.attribution import gene_row_locus as gr
from genomeos.attribution.gene_identity import GeneRow


def gene_row(name: str, n_tracks: int = 371, **extra):
    """A cached gene row in the chain's own field order, so the streaming reader sees what it sees."""
    return {
        "gene": name,
        "n_tracks": n_tracks,
        "mean_log2fc": 0.0,
        "max_drop_log2fc": -0.01,
        "max_drop_tissue": "heart",
        "max_rise_log2fc": 0.01,
        "max_rise_tissue": "K562",
        "by_cell": {"K562": 0.0},
        **extra,
    }


def element(eid: str, start: int, end: int, rows: list[dict], chrom: str = "chr1"):
    return {
        "id": eid,
        "chrom": chrom,
        "start": start,
        "end": end,
        "length": end - start,
        "genes_in_window": len(rows),
        "tracks": rows[0]["n_tracks"] if rows else 0,
        "seconds": 1.0,
        "genes": rows,
    }


def locus(gene_id: str, start: int, end: int, gene_type: str = "lncRNA", symbol: str = "N"):
    return GeneRow(gene_id=gene_id, symbol=symbol, gene_type=gene_type, start=start, end=end, tss=start)


def archive(tmp_path, name: str, records: list[dict], gz: bool = True):
    p = tmp_path / (f"{name}.json.gz" if gz else f"{name}.json")
    body = json.dumps({r["id"]: r for r in records})
    if gz:
        with gzip.open(p, "wt") as fh:
            fh.write(body)
    else:
        p.write_text(body)
    return p


# --- the gene row's contents -------------------------------------------------------------------


def test_gene_row_keys_are_the_eight_and_carry_no_locus():
    assert len(gr.GENE_ROW_KEYS) == 8
    for forbidden in ("gene_id", "ensembl", "strand", "transcript", "gene_type", "start", "end", "tss"):
        assert forbidden not in gr.GENE_ROW_KEYS


def test_row_key_vocabulary_finds_a_ninth_field_if_one_exists(tmp_path):
    plain = archive(tmp_path, "a", [element("E1", 100, 200, [gene_row("N")])])
    assert set(gr.row_key_vocabulary([plain])) >= set(gr.GENE_ROW_KEYS)
    assert "gene_id" not in gr.row_key_vocabulary([plain])

    withid = archive(tmp_path, "b", [element("E1", 100, 200, [gene_row("N", gene_id="ENSG1")])])
    assert gr.row_key_vocabulary([withid])["gene_id"] == 1


def test_row_key_vocabulary_spans_chunk_boundaries(tmp_path):
    """The reader streams, so a key split across two chunks must still be seen -- and seen once."""
    rows = [gene_row(f"G{i}") for i in range(4000)]
    p = archive(tmp_path, "big", [element("E1", 100, 200, rows)])
    assert p.stat().st_size > 0
    v = gr.row_key_vocabulary([p])
    assert v["gene"] == len(rows), v["gene"]
    assert v["n_tracks"] == len(rows)


def test_row_key_vocabulary_reads_a_plain_json_table_too(tmp_path):
    p = archive(tmp_path, "c", [element("E1", 100, 200, [gene_row("N")])], gz=False)
    assert "max_drop_tissue" in gr.row_key_vocabulary([p])


# --- the merge ---------------------------------------------------------------------------------


def test_a_row_above_the_column_count_is_a_merge(tmp_path):
    p = archive(
        tmp_path,
        "m",
        [element("E1", 100, 200, [gene_row("A", 740), gene_row("B", 371), gene_row("C", 370)])],
    )
    s = gr.row_scan([p])
    assert s["gene_rows"] == 3
    assert s["merged_rows"] == 1
    assert s["merged_names"] == {"A": 1}


def test_merged_at_least_is_a_lower_bound_not_a_count(tmp_path):
    """ceil(n_tracks / TRACKS) bands the row; it cannot be the number merged, because threshold=0.0
    drops a value of exactly zero and a constituent row may contribute fewer than TRACKS values."""
    p = archive(tmp_path, "m2", [element("E1", 1, 2, [gene_row("A", 485), gene_row("B", 1109)])])
    s = gr.row_scan([p])
    assert s["merged_at_least"] == {2: 1, 3: 1}
    assert s["merged_rows"] == 2
    assert "lower bound" in " ".join(gr.row_scan.__doc__.lower().split())


def test_exactly_the_column_count_is_not_a_merge(tmp_path):
    p = archive(tmp_path, "m3", [element("E1", 1, 2, [gene_row("A", gr.TRACKS)])])
    assert gr.row_scan([p])["merged_rows"] == 0


def test_ensembl_named_rows_are_counted(tmp_path):
    p = archive(tmp_path, "e", [element("E1", 1, 2, [gene_row("ENSG00000123"), gene_row("MN1")])])
    assert gr.row_scan([p])["ensembl_named_rows"] == 1


def test_the_two_drop_points_are_named_with_their_lines():
    files = {d["file"]: d["line"] for d in gr.DROPPED_AT}
    assert files["genomeos/predict/alphagenome_adapter.py"] == 227
    assert files["genomeos/predict/enhancer_target.py"] == 198


# --- the window --------------------------------------------------------------------------------


def test_loci_in_window_takes_the_gene_body_not_the_tss():
    """A gene whose TSS is outside the window but whose body reaches into it is a candidate: that is
    lane-identity's rule, and taking the TSS instead would drop it."""
    far = locus("ENSG1", 2_000_000, 2_000_100)
    reaching = locus("ENSG2", 2_000_000, 3_000_000)
    got = gr.loci_in_window([far, reaching], 2_500_000, 2_500_100, half=100_000)
    assert [g.gene_id for g in got] == ["ENSG2"]


# --- the classification, and the symbol-versus-locus test ---------------------------------------


def two_loci():
    return {"N": [locus("ENSG_A", 1_000, 2_000), locus("ENSG_B", 1_500, 2_500)]}


def test_one_locus_in_window_is_not_an_ambiguity():
    rec = element("E1", 1_200, 1_300, [gene_row("N")])
    assert gr.classify_element(rec, {"N": [locus("ENSG_A", 1_000, 2_000)]}, set()) == []


def test_a_merged_row_is_not_settled():
    rec = element("E1", 1_700, 1_800, [gene_row("N", 740)])
    got = gr.classify_element(rec, two_loci(), set())
    assert [p["verdict"] for p in got] == ["merged"]
    assert got[0]["the_symbol_row_is_then"] is None


def test_no_signal_at_all_is_nothing():
    rec = element("E1", 1_700, 1_800, [gene_row("N")])
    got = gr.classify_element(rec, two_loci(), set())
    assert [p["verdict"] for p in got] == ["nothing"]


def test_pinning_distinguishes_two_loci_of_the_same_name():
    """The test that matters. Both loci carry the symbol N; the cache names one of them by its bare
    Ensembl id, so the row named N can only be the other. A property of the symbol could not do
    this, because the symbol is the same for both."""
    rec = element("E1", 1_700, 1_800, [gene_row("N"), gene_row("ENSG_A")])
    got = gr.classify_element(rec, two_loci(), set())
    assert [p["verdict"] for p in got] == ["pinned"]
    assert got[0]["the_symbol_row_is_then"] == "ENSG_B"
    assert got[0]["named_by_bare_ensembl_id_in_this_element"] == ["ENSG_A"]


def test_pinning_needs_every_other_locus_named_not_just_one():
    three = {
        "N": [
            locus("ENSG_A", 1_000, 2_000),
            locus("ENSG_B", 1_500, 2_500),
            locus("ENSG_C", 1_600, 2_600),
        ]
    }
    rec = element("E1", 1_700, 1_800, [gene_row("N"), gene_row("ENSG_A")])
    got = gr.classify_element(rec, three, set())
    assert [p["verdict"] for p in got] == ["partly_pinned"]
    assert got[0]["the_symbol_row_is_then"] is None

    rec2 = element("E2", 1_700, 1_800, [gene_row("N"), gene_row("ENSG_A"), gene_row("ENSG_C")])
    got2 = gr.classify_element(rec2, three, set())
    assert [p["verdict"] for p in got2] == ["pinned"]
    assert got2[0]["the_symbol_row_is_then"] == "ENSG_B"


def test_a_merged_row_is_merged_even_with_an_id_sibling_present():
    """The merge wins: whatever else the element names, pooled numbers stay pooled."""
    rec = element("E1", 1_700, 1_800, [gene_row("N", 740), gene_row("ENSG_A")])
    assert [p["verdict"] for p in gr.classify_element(rec, two_loci(), set())] == ["merged"]


def test_coding_name_is_a_property_of_the_symbol_and_pins_nothing():
    """The named trap. One of the two loci is protein_coding, which is exactly lane-identity's 295 of
    547: it is recorded, it changes no verdict, and it names no locus."""
    mixed = {"N": [locus("ENSG_A", 1_000, 2_000, "protein_coding"), locus("ENSG_B", 1_500, 2_500)]}
    rec = element("E1", 1_700, 1_800, [gene_row("N")])
    got = gr.classify_element(rec, mixed, {"N"})
    assert got[0]["verdict"] == "nothing"
    assert got[0]["the_symbol_row_is_then"] is None
    assert got[0]["name_is_protein_coding_somewhere_on_the_chromosome"] is True
    # and the same row with no coding locus gets the same verdict: the flag moves nothing
    plain = gr.classify_element(rec, two_loci(), set())
    assert plain[0]["verdict"] == got[0]["verdict"]


# --- the tally ---------------------------------------------------------------------------------


def test_the_verdicts_sum_to_the_ambiguous_pairs_with_nothing_in_a_residue():
    pairs = []
    for v in ("merged", "pinned", "nothing", "partly_pinned", "pinned"):
        pairs.append(
            {
                "verdict": v,
                "name": "N",
                "loci_in_window": 2,
                "name_is_protein_coding_somewhere_on_the_chromosome": v == "nothing",
            }
        )
    t = gr.tally(pairs)
    assert t["ambiguous_element_name_pairs"] == 5
    assert sum(t["by_verdict"].values()) == 5
    assert t["settled_by_what_is_on_disk"] == 2
    assert t["distinct_pinned_names"] == 1
    coding = t["restricted_to_names_with_a_protein_coding_locus_on_the_chromosome"]
    assert coding["ambiguous"] == 1
    assert coding["pinned"] == 0


def test_the_coding_restriction_is_a_subset_of_the_whole():
    pairs = [
        {
            "verdict": "pinned",
            "name": "N",
            "loci_in_window": 2,
            "name_is_protein_coding_somewhere_on_the_chromosome": True,
        }
    ]
    t = gr.tally(pairs)
    c = t["restricted_to_names_with_a_protein_coding_locus_on_the_chromosome"]
    assert c["pinned"] == t["by_verdict"]["pinned"] == 1
    assert c["ambiguous"] == t["ambiguous_element_name_pairs"] == 1


def test_an_empty_tally_reports_zero_and_not_an_error():
    t = gr.tally([])
    assert t["ambiguous_element_name_pairs"] == 0
    assert t["settled_by_what_is_on_disk"] == 0
    assert set(t["by_verdict"]) == set(gr.CLASSES)


def test_every_verdict_has_a_stated_meaning():
    assert set(gr.CLASSES) == {"merged", "pinned", "partly_pinned", "nothing"}
    for v in gr.CLASSES.values():
        assert len(v) > 20


# --- what the figures are not -------------------------------------------------------------------


def test_what_this_is_not_names_the_denominator_and_the_547():
    w = gr.what_this_is_not()
    assert "547" in w["the_547"]
    assert "440,589" in w["the_denominator"] and "NOT" in w["the_denominator"]
    assert "lost distinction" in w["merged_rows"]
    assert "MATR3" in w["pinned"]
    for v in w.values():
        assert len(v) > 40


def test_a_record_without_coordinates_is_skipped_not_guessed():
    rec = {"id": "E1", "chrom": "chr1", "genes": [gene_row("N")]}
    assert gr.classify_element(rec, two_loci(), set()) == []


def test_rows_without_a_name_are_skipped():
    rec = element("E1", 1_700, 1_800, [{"gene": "", "n_tracks": 371}])
    assert gr.classify_element(rec, two_loci(), set()) == []


@pytest.mark.parametrize(
    "name,hit", [("ENSG00000123456", True), ("ENSG1", True), ("MN1", False), ("ENSGX", False)]
)
def test_the_bare_ensembl_shape(name, hit):
    assert bool(gr.ENSG.match(name)) is hit
