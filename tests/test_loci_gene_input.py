"""The versioned gene-input reader: v1's identity, v2's repair, and the plant that separates them.

Three things are proved here and they are deliberately not the same thing.

1. **v1 is the reader as it stood**, and that is checked two ways that fail independently. The
   pre-change loop is transcribed verbatim from `genomeos/benchmark/loci.py` as it stood before
   `reading` existed and is asserted equal to `credit_rows(rows, GENE_INPUT_V1)` over an exhaustive
   enumeration of head shapes and a seeded random sweep, so the claim is about ALL inputs and not
   about the ones that happened to be convenient. Separately, every committed benchmark frame's
   saved `gene_input` reading is read off disk and compared: the key set, the one numeric invariant
   that only v1 can satisfy, and - for the loci whose element rows are reachable without opening a
   bulk per-chromosome archive - the whole dict field for field.
2. **v2 is opt-in and changes nothing by default.** `score_locus` asks for no reading, so it gets
   v1; a v1 dict carries no version key at all.
3. **the plant**, on a synthetic fixture and never on a benchmark frame: rows built to the shape of
   the fourth frame's deletion-layer overlaps, read under both versions, and classified by
   `loci_miss.classify_gene` - bf33233's own classifier, imported and not restated. The fixture
   also records why the four eQTL overlaps section 25 names CANNOT be moved by any version of this
   reader, which is a finding about the repair's reach.

Every input is a tracked file under `data/results` or is built in the test. No git-ignored store is
required, no bulk chromosome archive is opened, no cached-element loader is called, 0 AlphaGenome
requests are spent and no rate is computed anywhere in this file.
"""

from __future__ import annotations

import json
import os
import random
import resource
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest

from genomeos.benchmark import loci as lb
from genomeos.benchmark import loci_gene_input as gi
from genomeos.benchmark import loci_miss

RESULTS = Path("data/results")

#: the six committed frames section 24 and section 25 report over, named by their result stem.
FRAMES = (
    "loci_benchmark",
    "loci_candidates",
    "loci_third",
    "loci_fourth",
    "loci_noncoding",
    "loci_fourth_rejoined",
)

#: the result files `loci._deletion_rows` reads that are TRACKED in git. The per-chromosome element
#: tables of the `enhancer_targets_all` sweep are not here on purpose: they are 30 MB+ archives in a
#: git-ignored store, and this file is forbidden to open one. A results directory holding only these
#: makes `_deletion_rows` return the rows it can see and no others, with no archive touched.
TRACKED_ROW_SOURCES = (
    "loci_stated_intervals",
    "loci_candidate_intervals",
    "loci_third_intervals",
    "loci_fourth_intervals",
    "loci_deletions",
)

#: every key a v1 reading may carry. A committed frame holding any other key would mean either that
#: a v2 reading was saved into it or that the reader's output shape moved; both are failures here.
V1_KEYS = frozenset(
    {
        "layer",
        "provenance",
        "window",
        "elements_scored",
        "naming_a_coding_gene",
        "genes",
        "published_targets",
        "target",
        "rank_of_first_published_target",
        "pending",
        "evidence",
    }
)

#: keys only a v2 reading emits. None may appear in a committed frame.
V2_ONLY_KEYS = (
    "reading",
    "naming_a_gene",
    "coding_head_genes",
    "coding_head_target",
    "coding_head_rank_of_first_published_target",
    "coding_head_reading",
)

#: Registered before the comparison was run: this file reads the six frames (about 4.3 MB of JSON)
#: and `loci_deletions.json` (1.8 MB) one at a time and holds no chromosome table at all, so its
#: own peak resident set has no business exceeding this. The check RAISES; it does not warn.
RSS_CEILING_MB = 900


def _rss_mb() -> float:
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss / 1e6 if rss > 1e7 else rss / 1024  # Darwin reports bytes, Linux kibibytes


def _under_the_ceiling(where: str) -> None:
    got = _rss_mb()
    if got > RSS_CEILING_MB:
        raise AssertionError(f"{where}: peak RSS {got:.0f} MB over the registered {RSS_CEILING_MB} MB")


def _frame(stem: str) -> dict[str, Any]:
    p = RESULTS / f"{stem}.json"
    if not p.exists():  # tracked in git, so this is a broken checkout rather than a missing store
        pytest.skip(f"{p} is absent from this checkout")
    return json.loads(p.read_text())


def _saved_gene_input_readings() -> list[tuple[str, str, dict[str, Any], dict[str, Any]]]:
    """(frame, locus, the saved gene_input reading, the locus's `expected` block) for all six."""
    out = []
    for stem in FRAMES:
        for row in _frame(stem).get("loci") or []:
            reading = (row.get("readings") or {}).get("gene_input")
            if reading is not None:
                out.append((stem, row["locus"], reading, row["expected"]))
    _under_the_ceiling("reading the six committed frames")
    return out


# ----------------------------------------------------------------------------------------------
# 1. v1 is the reader as it stood, proved against the implementation it replaced
# ----------------------------------------------------------------------------------------------


def _reader_before_the_change(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """`loci.read_gene_input`'s accumulator loop, transcribed verbatim from the tree as it stood
    before `reading` existed. Not imported, not refactored, not shared with the code under test:
    this is a second copy on purpose, so that a change to `credit_rows` cannot change it too."""
    by_gene: dict[str, dict[str, Any]] = {}
    for e in rows:
        p = e.get("predicted_coding") or {}
        if not p.get("gene"):
            continue
        g = by_gene.setdefault(p["gene"], {"gene": p["gene"], "elements": 0, "activating": 0, "summed": 0.0})
        g["elements"] += 1
        g["activating"] += int(p["action"] == "activates")
        g["summed"] = round(g["summed"] + abs(p["log2_fold_change"]), 3)
    return by_gene


def _head(gene: str | None, action: str, lfc: float) -> dict[str, Any] | None:
    return None if gene is None else {"gene": gene, "action": action, "log2_fold_change": lfc}


def _exhaustive_row_space() -> list[list[dict[str, Any]]]:
    """Every shape one element's two heads can take, and then every pair of them.

    The axes are the ones the two readings can differ on: whether each head names a gene at all,
    whether the two heads name the SAME gene, the action word, and the sign of the effect.
    """
    genes = (None, "CODA", "NONC")
    singles: list[dict[str, Any]] = []
    for ck in genes:
        for pk in genes:
            for action in ("activates", "represses"):
                for lfc in (-1.691, 0.1433):
                    e: dict[str, Any] = {"id": f"E{len(singles)}", "start": 0, "end": 100}
                    if ck is not None:
                        e["predicted_coding"] = _head(ck, action, lfc)
                    if pk is not None:
                        e["predicted"] = _head(pk, "activates", -1.691)
                    singles.append(e)
    space = [[e] for e in singles]
    space += [[a, b] for a in singles[:12] for b in singles[:12]]
    space.append([])
    space.append([{"id": "bare", "start": 0, "end": 1}])  # an element with neither head
    return space


def test_v1_is_the_loop_it_replaced_on_every_shape_two_heads_can_take() -> None:
    space = _exhaustive_row_space()
    assert len(space) > 150, "the enumeration shrank; it is meant to cover both heads' whole space"
    for rows in space:
        assert gi.credit_rows(rows, gi.GENE_INPUT_V1) == _reader_before_the_change(rows)
        assert gi.credit_rows(rows) == _reader_before_the_change(rows), "v1 must be the DEFAULT"


def test_v1_is_the_loop_it_replaced_over_a_seeded_random_sweep() -> None:
    rng = random.Random(20261003)
    names = ["AAA", "BBB", "CCC", "DDD", "ENSG00000240739"]
    for _ in range(4000):
        rows = []
        for i in range(rng.randrange(0, 6)):
            e: dict[str, Any] = {"id": f"E{i}", "start": i * 10, "end": i * 10 + 5}
            for key in gi.HEADS:
                if rng.random() < 0.8:
                    e[key] = _head(
                        rng.choice(names),
                        rng.choice(("activates", "represses")),
                        round(rng.uniform(-3, 3), 4),
                    )
            rows.append(e)
        assert gi.credit_rows(rows, gi.GENE_INPUT_V1) == _reader_before_the_change(rows)


def test_no_committed_frame_carries_a_v2_key_and_every_one_satisfies_v1s_own_invariant() -> None:
    """The frames' actual numbers, not the code's appearance.

    Two claims per locus. First the key set: a v1 reading carries no version stamp, so a committed
    frame holding one would mean a v2 reading had been saved into it. Second the invariant that
    separates the readings arithmetically - under v1 each element credits AT MOST ONE gene, so the
    element credits can never exceed the elements scored, while under v2 an element whose two heads
    name different genes credits twice and the sum can exceed it. Every committed frame satisfying
    the v1 bound is evidence from the published numbers that they are v1 readings.
    """
    saved = _saved_gene_input_readings()
    assert len(saved) == 147, f"the six frames held 147 gene_input readings, found {len(saved)}"
    for stem, locus, reading, _expected in saved:
        extra = set(reading) - V1_KEYS
        assert not extra, f"{stem}/{locus} carries non-v1 keys {sorted(extra)}"
        for key in V2_ONLY_KEYS:
            assert key not in reading, f"{stem}/{locus} carries the v2 key {key}"
        assert reading["naming_a_coding_gene"] <= reading["elements_scored"], (
            f"{stem}/{locus} credits {reading['naming_a_coding_gene']} over"
            f" {reading['elements_scored']} elements, which v1 cannot do"
        )


def test_v1s_assembly_reproduces_every_committed_readings_target_and_rank() -> None:
    """`target` and `rank_of_first_published_target` recomputed from each saved `genes` list.

    `genes` is the top eight, so a locus whose first published target ranks ninth or lower cannot be
    checked from the saved dict at all: those are counted and named, not asserted around. 114 of the
    147 are checkable and all 114 must agree.
    """
    saved = _saved_gene_input_readings()
    checked = beyond_the_eight = 0
    for stem, locus, reading, expected in saved:
        ranked = reading["genes"]
        rank = reading["rank_of_first_published_target"]
        if rank is not None and rank > len(ranked):
            beyond_the_eight += 1
            continue
        targets = expected["targets"]
        mine = [g for g in ranked if g["gene"] in targets]
        assert gi.rank_of_first_target(ranked, targets) == rank, f"{stem}/{locus} rank"
        want = mine[0]["gene"] if mine else (ranked[0]["gene"] if ranked else None)
        assert reading["target"] == want, f"{stem}/{locus} target"
        checked += 1
    assert checked == 114, f"expected 114 checkable loci, got {checked}"
    assert checked + beyond_the_eight == len(saved)


def _tracked_only_results_dir(tmp: Path) -> Path:
    for name in TRACKED_ROW_SOURCES:
        src = RESULTS / f"{name}.json"
        if src.exists():
            os.symlink(src.resolve(), tmp / f"{name}.json")
    return tmp


def test_every_committed_reading_whose_rows_are_reachable_replays_byte_for_byte_under_v1() -> None:
    """The strongest form available without opening a bulk archive: the whole dict, re-derived.

    `loci._deletion_rows` is called against a results directory holding only the TRACKED row
    sources, so the 30 MB+ per-chromosome element tables of the sweep are not on the path and are
    never opened. Where that directory yields the same number of rows the frame recorded, the row
    set is the one the frame read and the whole reading is rebuilt and compared field for field.
    Where it yields fewer, the locus's rows live in the git-ignored archive: it is counted and
    SKIPPED BY NAME rather than asserted about, because a fresh checkout has no archive at all.
    """
    saved = _saved_gene_input_readings()
    with tempfile.TemporaryDirectory() as raw:
        results_dir = _tracked_only_results_dir(Path(raw))
        replayed: list[str] = []
        unreachable: list[str] = []
        for stem, locus, reading, expected in saved:
            lo, hi = reading["window"]
            rows = lb._deletion_rows(expected["chrom"], lo, hi, results_dir)
            if len(rows) != reading["elements_scored"]:
                unreachable.append(f"{stem}/{locus}")
                continue
            ranked = gi.rank(gi.credit_rows(rows, gi.GENE_INPUT_V1))
            mine = [g for g in ranked if g["gene"] in expected["targets"]]
            rebuilt = {
                "layer": "gene_input",
                "provenance": "derived",
                "window": list(reading["window"]),
                "elements_scored": len(rows),
                "naming_a_coding_gene": sum(g["elements"] for g in ranked),
                "genes": ranked[:8],
                "published_targets": mine,
                "target": mine[0]["gene"] if mine else (ranked[0]["gene"] if ranked else None),
                "rank_of_first_published_target": gi.rank_of_first_target(ranked, expected["targets"]),
                "pending": None if rows else "no already-scored element anywhere in the locus window",
                "evidence": (
                    "predicted: AlphaGenome deletions already computed, summed |log2| per coding gene"
                ),
            }
            assert rebuilt == reading, f"{stem}/{locus} does not replay identically under v1"
            replayed.append(f"{stem}/{locus}")
    _under_the_ceiling("replaying the reachable loci")
    assert len(replayed) == 5, (
        "5 of the 147 loci have their rows wholly in tracked results and all 5 must replay;"
        f" replayed {sorted(replayed)}, unreachable {len(unreachable)}"
    )
    assert len(unreachable) == 142


# ----------------------------------------------------------------------------------------------
# 2. v2 is opt-in, and asking for nothing gets v1
# ----------------------------------------------------------------------------------------------


@dataclass
class _Ch:
    """The one attribute `read_gene_input` uses of a chromosome, and nothing else."""

    chrom: str


def _expect(chrom: str, window: tuple[int, int], targets: tuple[str, ...]) -> Any:
    return lb.Expect(
        locus="SYNTHETIC",
        chrom=chrom,
        element=(window[0] + 1000, window[0] + 1300),
        window=window,
        classes=("synthetic",),
        targets=targets,
        direction="activates",
        tissues=("none",),
    )


def test_no_call_anywhere_in_the_tree_asks_this_reader_for_v2() -> None:
    """Parsed, not grepped, so a docstring naming v2 does not read as a caller asking for it.

    v2 is opt-in from outside the benchmark. Every `read_gene_input` call in the repository is
    inspected and none may pass `reading`; `score_locus`'s own call is checked by name.
    """
    import ast

    asking: list[str] = []
    for path in sorted(Path("genomeos").rglob("*.py")) + sorted(Path("scripts").rglob("*.py")):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            fn = node.func
            name = fn.attr if isinstance(fn, ast.Attribute) else getattr(fn, "id", None)
            if name != "read_gene_input":
                continue
            if any(k.arg == "reading" for k in node.keywords) or len(node.args) > 5:
                asking.append(f"{path}:{node.lineno}")
    assert asking == [], f"these callers ask for a non-default reading: {asking}"
    src = Path("genomeos/benchmark/loci.py").read_text()
    assert 'readings["gene_input"] = read_gene_input(ch, e, results_dir)' in src, (
        "score_locus's call to read_gene_input has moved; if it now names a reading, this test is"
        " what should have stopped that"
    )


def test_the_default_signature_is_v1_and_a_v1_reading_carries_no_version_key() -> None:
    import inspect

    sig = inspect.signature(lb.read_gene_input)
    assert sig.parameters["reading"].default == gi.GENE_INPUT_V1


def test_an_unknown_reading_raises_rather_than_falling_back() -> None:
    with pytest.raises(ValueError, match="unknown gene-input reading"):
        gi.credit_rows([], "v3")


# ----------------------------------------------------------------------------------------------
# 3. the plant, on a synthetic fixture
# ----------------------------------------------------------------------------------------------


@dataclass
class _Locus:
    start: int
    end: int


@dataclass
class _Gene:
    symbol: str
    locus: _Locus


class _Annotation:
    """The attributes `loci_miss._record` and `classify_gene` use of a GENCODE annotation.

    Synthetic, so the plant needs no `data/reference` store and classifies no real coordinate. The
    classifier itself is `loci_miss`'s, imported: the overlap test and the distance bands are
    bf33233's and are not restated here.
    """

    def __init__(self, bodies: dict[str, tuple[int, int]]) -> None:
        self.genes = {s: _Gene(s, _Locus(a, b)) for s, (a, b) in bodies.items()}


#: the shape of the fourth frame's SLC2A3 overlap, on synthetic coordinates: a 315 bp non-coding
#: gene lying wholly inside the published target's body, and a coding neighbour 200 kb away. Those
#: are the two numbers section 25's own table gives for that locus (315 bp, one body inside the
#: other); the distance is chosen to fall outside every tolerated band so that the v1 arm's class is
#: unambiguous. Nothing here is read off a benchmark frame.
PLANT_BODIES = {
    "SLC2A3": (7_980_000, 8_000_000),
    "ENSG00000240739": (7_990_000, 7_990_315),
    "NEIGHBOUR": (8_200_000, 8_250_000),
}


def _plant_rows() -> list[dict[str, Any]]:
    """One already-scored element whose two heads disagree in exactly the way the repair is about.

    `predicted` names the non-coding gene inside the target's body at the larger effect;
    `predicted_coding` names the coding neighbour at a twelfth of it. The two magnitudes are the
    ones `read_deletion`'s docstring records for the chr8 CCDC26 element that caught this defect
    with a positive - -1.691 against -0.1433 - so the fixture is the documented shape and not a
    convenient one.
    """
    return [
        {
            "id": "PLANT1",
            "chrom": "chrS",
            "start": 7_990_100,
            "end": 7_990_300,
            "predicted_coding": {"gene": "NEIGHBOUR", "action": "activates", "log2_fold_change": -0.1433},
            "predicted": {"gene": "ENSG00000240739", "action": "activates", "log2_fold_change": -1.691},
        }
    ]


def test_the_plant_v1_names_the_distant_coding_gene_and_v2_names_the_gene_inside_the_target() -> None:
    rows = _plant_rows()
    v1 = gi.rank(gi.credit_rows(rows, gi.GENE_INPUT_V1))
    v2 = gi.rank(gi.credit_rows(rows, gi.GENE_INPUT_V2))
    assert [g["gene"] for g in v1] == ["NEIGHBOUR"]
    assert [g["gene"] for g in v2] == ["ENSG00000240739", "NEIGHBOUR"]

    ann = _Annotation(PLANT_BODIES)
    targets = ["SLC2A3"]
    first_v1 = loci_miss.classify_gene(ann, v1[0]["gene"], targets, [])
    first_v2 = loci_miss.classify_gene(ann, v2[0]["gene"], targets, [])
    # bf33233's classifier, on bf33233's classes. The v1 arm's class is a statement about where the
    # gene is; nothing in this test turns either class into a rate.
    assert first_v1["class"] == "elsewhere"
    assert first_v1["gap_bp"] == 200_000
    assert first_v2["class"] == "overlaps_the_target_body"
    assert first_v2["gap_bp"] == 0
    assert first_v2["nearest_target"] == "SLC2A3"
    # and both are misses under the benchmark's own rule, which is symbol equality and is untouched.
    assert first_v1["class"] != "exact" and first_v2["class"] != "exact"


def test_the_plant_shows_a_coding_target_can_be_displaced_and_never_promoted() -> None:
    """`REPAIR`'s registered direction, demonstrated rather than asserted in prose.

    The fourth frame's drawing rule requires a protein-coding published target, so the only thing
    v2 can do to that frame's gene-input layer is take a hit away. Both arms of that are planted.
    """
    # the published target IS the coding head's gene: v1 names it, v2 is displaced by the stronger
    # non-coding effect. A loss, and it is the expected direction.
    rows = _plant_rows()
    v1 = gi.rank(gi.credit_rows(rows, gi.GENE_INPUT_V1))
    v2 = gi.rank(gi.credit_rows(rows, gi.GENE_INPUT_V2))
    assert v1[0]["gene"] == "NEIGHBOUR"
    assert v2[0]["gene"] != "NEIGHBOUR"
    # the reverse never happens: a gene the coding head names can only lose rank under v2, because
    # v2's ranking is v1's with further entries added at non-negative credit.
    rng = random.Random(7)
    for _ in range(500):
        rs = []
        for i in range(rng.randrange(1, 5)):
            e: dict[str, Any] = {"id": f"E{i}"}
            for key in gi.HEADS:
                if rng.random() < 0.85:
                    e[key] = _head(
                        rng.choice(("P", "Q", "R", "S")),
                        rng.choice(("activates", "represses")),
                        round(rng.uniform(-2, 2), 3),
                    )
            rs.append(e)
        one = gi.credit_rows(rs, gi.GENE_INPUT_V1)
        two = gi.credit_rows(rs, gi.GENE_INPUT_V2)
        for gene, got in one.items():
            assert gene in two, "v2 dropped a gene v1 credited, which it cannot do"
            assert two[gene]["summed"] >= got["summed"] - 1e-9
        assert set(one) <= set(two)


def test_v2_is_carried_beside_v1_in_one_call_so_the_two_movements_stay_separable() -> None:
    """The `window_reading` precedent: added beside the layer, not in place of it."""
    rows = _plant_rows()
    with tempfile.TemporaryDirectory() as raw:
        d = Path(raw)
        (d / "loci_deletions.json").write_text(json.dumps({"elements": [{**rows[0]}]}))
        expect = _expect("chrS", (7_900_000, 8_300_000), ("SLC2A3",))
        v1 = lb.read_gene_input(_Ch("chrS"), expect, d)
        v2 = lb.read_gene_input(_Ch("chrS"), expect, d, reading=gi.GENE_INPUT_V2)
    assert set(v1) == V1_KEYS
    assert v1["target"] == "NEIGHBOUR"
    assert v2["reading"] == gi.GENE_INPUT_V2
    assert v2["target"] == "ENSG00000240739"
    # the pre-repair answer travels in the same dict, so no rate published under v1 needs the git
    # history to be recomputed.
    assert v2["coding_head_target"] == "NEIGHBOUR"
    assert v2["coding_head_genes"] == v1["genes"]
    assert v2["coding_head_rank_of_first_published_target"] == v1["rank_of_first_published_target"]
    assert v2["naming_a_coding_gene"] == v1["naming_a_coding_gene"] == 1
    assert v2["naming_a_gene"] == 2
    assert v2["elements_scored"] == v1["elements_scored"] == 1


def test_the_four_eqtl_overlaps_cannot_be_moved_by_any_version_of_this_reader() -> None:
    """Section 25's sentence, checked against the committed result rather than taken.

    The item says four of the fourth frame's seven overlaps come from the eQTL layer. They do - and
    that is precisely why they are out of this repair's reach: `read_eqtl` ranks GTEx p-values and
    reads neither stored head, so no reading of `read_gene_input` can change an eQTL-layer answer.
    The gene-input layer contributes NONE of the seven. The plant above is synthetic for this
    reason; there is no arrangement of these frames in which the repair moves those four.
    """
    miss = _frame("loci_miss")
    overlaps = miss["frames"]["loci_fourth"]["overlapping_misses"]
    by_layer: dict[str, int] = {}
    for o in overlaps:
        by_layer[o["layer"]] = by_layer.get(o["layer"], 0) + 1
    assert len(overlaps) == 7
    assert by_layer == {"eqtl": 4, "deletion": 3}
    assert "gene_input" not in by_layer
    src = Path("genomeos/benchmark/loci.py").read_text()
    body = src.split("def read_eqtl(")[1].split("\ndef ")[0]
    assert "predicted" not in body, "read_eqtl now reads a predicted head; this finding would change"


def test_this_file_imports_the_governing_prohibition_and_computes_no_rate() -> None:
    """bf33233's prohibition is imported, not paraphrased, and it binds this work.

    The distance classes are description only and feed no rate, ever. Nothing in this file or in
    `loci_gene_input` divides one count by another, and `REPAIR` names the registration it is bound
    by rather than restating its rule.
    """
    assert "feed no rate, ever" in loci_miss.PREREGISTRATION["only_overlap_can_be_tolerated"]
    assert (
        "score the nearest-gene baseline as a hit"
        in loci_miss.PREREGISTRATION["only_overlap_can_be_tolerated"]
    )
    assert loci_miss.TOLERATED == ("exact", "overlaps_the_target_body")
    assert "bf33233" in gi.REPAIR["the_governing_rules_are_imported_not_restated"]


def test_the_registration_quotes_the_item_it_answers_and_claims_no_measurement() -> None:
    assert "still reads `predicted_coding` alone" in gi.REPAIR["the_item"]
    assert "no rate is recomputed" in gi.REPAIR["written"]
    assert (
        "restriction of the same computation rather than a better"
        in (gi.REPAIR["why_this_reading_and_not_the_favourable_one"])
    )
    assert "two movements behind one number" in gi.REPAIR["why_this_reading_and_not_the_favourable_one"]
    assert gi.REPAIR["cost"].startswith("0 AlphaGenome requests")
