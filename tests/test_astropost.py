"""AstroREG-2's post-run analyses: the gate, item (h)'s registered statistics, the primary claim's tail.

EVERY ANSWER HERE IS SYNTHETIC. No test reads the real ledger or the real request root: the ledgers and
answers are written under a repository root of the test's own, and a tripwire on the cache loader shows
that a refused run is refused before anything of the outcome is opened.
"""

from __future__ import annotations

import importlib.util
import json
import random
import sys
from pathlib import Path

import pytest

from genomeos.attribution import astropost, astrorun

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
_SPEC = importlib.util.spec_from_file_location("astroreg2_analyse", ROOT / "scripts/astroreg2_analyse.py")
analyse = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = analyse
_SPEC.loader.exec_module(analyse)

BANDS = astropost.item_h_registration()["bands"]
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def plan_row(i: int, chrom: str) -> dict:
    return {
        "chrom": chrom,
        "element": f"EH38E{i:07d}",
        "start": i * 1000,
        "end": i * 1000 + 300,
        "serves_genes": ["AAA"],
        "serves_labels": ["negative"],
    }


def synthetic_answer(eid: str, chrom: str, tissue: str | None) -> dict:
    """A written entry's gene rows in the recording path's schema. `tissue=None` moves no gene."""
    drop, rise = (-0.5, 0.02) if tissue else (-0.01, 0.01)
    return {
        "id": eid,
        "chrom": chrom,
        "genes": [
            {
                "gene": "AAA",
                "n_tracks": 3,
                "mean_log2fc": -0.1,
                "max_drop_log2fc": drop,
                "max_drop_tissue": tissue or "K562",
                "max_rise_log2fc": rise,
                "max_rise_tissue": "HepG2",
                "by_cell": {"K562": -0.1},
            }
        ],
        "model": {"gene_axis_outputs": [{"tracks_total": 371, "rows": []}]},
    }


class World:
    """A repository root of the test's own: a run-2 ledger, a request root and a reviewed list."""

    def __init__(self, tmp_path, monkeypatch, n=46, tissue_of=None):
        self.repo = tmp_path / "repo"
        monkeypatch.setattr(astrorun, "ROOT_FOR_BLOBS", self.repo)
        self.cache = self.repo / astrorun.ASTROREG2_CACHE_ROOT
        self.ledger = astrorun.ledger_path_for_run(2)
        self.plan = [plan_row(i, CHROMS[i % len(CHROMS)]) for i in range(n)]
        p = tmp_path / "plan.json"
        p.write_text(json.dumps({"requests": self.plan}))
        monkeypatch.setattr(analyse.astroreg2_send, "PLAN", p)
        monkeypatch.setattr(analyse.astroreg2_send, "REVIEWED_DIGEST", astrorun.plan_digest(self.plan))
        self.tissue_of = tissue_of or (lambda i: "astrocyte" if i % 3 == 0 else "K562")
        self.saved: dict = {}

        def save_result(name, out):
            self.saved[name] = out
            return Path("not-written")

        import genomeos.results as gr

        monkeypatch.setattr(gr, "save_result", save_result)
        monkeypatch.setattr(analyse, "manifest", lambda ledger, plan: {"stubbed": True})

    def run_to(self, charged, answered, complete, **complete_row):
        budget = astrorun.RequestBudget(self.ledger, cap=len(self.plan))
        for i, r in enumerate(self.plan[:charged]):
            budget.take(chrom=r["chrom"], element=r["element"])
            if i < answered:
                d = self.cache / r["chrom"]
                d.mkdir(parents=True, exist_ok=True)
                hit = synthetic_answer(r["element"], r["chrom"], self.tissue_of(i))
                (d / f"{r['element']}.json").write_text(json.dumps(hit))
                budget.note(event="answer", id=r["element"], usable=True)
        if complete:
            astrorun.mark_run_complete(self.ledger, requests=charged, **complete_row)
        return self


def run_main(monkeypatch, capsys, *args):
    monkeypatch.setattr(sys, "argv", ["astroreg2_analyse.py", *args])
    code = analyse.main()
    return code, capsys.readouterr().out


@pytest.fixture
def no_answer_may_be_opened(monkeypatch):
    """The tripwire: a refused run must be refused before any answer is loaded."""
    import genomeos.predict.enhancer_target as et

    def tripped(*a, **k):
        raise AssertionError("an answer was opened before the gate refused")

    monkeypatch.setattr(et, "load_cached", tripped)


# ------------------------------------------------------------------------------------------ the gate


class TestTheGateRefusesAnIncompleteRun:
    def test_PLANTED_run_2s_own_shape_REFUSES_and_opens_no_answer(
        self, tmp_path, monkeypatch, capsys, no_answer_may_be_opened
    ):
        """Request/answer pairs, the last request unanswered, no stop event: run 2 as it stands."""
        World(tmp_path, monkeypatch).run_to(charged=14, answered=13, complete=False)
        code, said = run_main(monkeypatch, capsys, "--run", "2")
        assert code == 2
        assert "REFUSED: the ledger" in said
        assert "does not show the run complete" in said
        assert "it has no `run_complete` line" in said
        assert "14 charge(s) and 13 answer(s)" in said
        assert "its last line is `request` n=14" in said
        assert "decide nothing on their own" in said

    def test_no_ledger_at_all_REFUSES(self, tmp_path, monkeypatch, capsys, no_answer_may_be_opened):
        World(tmp_path, monkeypatch)
        code, said = run_main(monkeypatch, capsys, "--run", "2")
        assert code == 2
        assert "there is no ledger" in said

    def test_PLANTED_a_pilot_stop_REFUSES(self, tmp_path, monkeypatch, capsys, no_answer_may_be_opened):
        World(tmp_path, monkeypatch).run_to(10, 10, True, stopped_at_pilot=True)
        code, said = run_main(monkeypatch, capsys, "--run", "2")
        assert code == 2
        assert "PILOT STOP" in said

    def test_PLANTED_lines_after_run_complete_REFUSE(
        self, tmp_path, monkeypatch, capsys, no_answer_may_be_opened
    ):
        w = World(tmp_path, monkeypatch).run_to(12, 12, True)
        astrorun.RequestBudget(w.ledger, cap=99).take(chrom="chr1", element="late")
        code, said = run_main(monkeypatch, capsys, "--run", "2")
        assert code == 2
        assert "is not the last line" in said

    def test_PLANTED_a_cached_answer_the_ledger_never_answered_REFUSES(self, tmp_path, monkeypatch, capsys):
        w = World(tmp_path, monkeypatch).run_to(40, 40, True)
        extra = w.plan[44]
        d = w.cache / extra["chrom"]
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{extra['element']}.json").write_text(
            json.dumps(synthetic_answer(extra["element"], "c", "K562"))
        )
        code, said = run_main(monkeypatch, capsys, "--run", "2")
        assert code == 2
        assert "the request root and the ledger disagree" in said
        assert extra["element"] in said

    def test_ALLOWED_a_complete_synthetic_run_is_analysed(self, tmp_path, monkeypatch, capsys):
        w = World(tmp_path, monkeypatch).run_to(charged=46, answered=45, complete=True)
        code, said = run_main(monkeypatch, capsys, "--run", "2")
        assert code == 0, said
        out = w.saved[analyse.RESULT]
        assert out["gate"]["complete"] is True
        assert out["gate"]["charged_without_an_answer"] == [w.plan[45]["element"]]
        h = out["item_h"]
        assert h["population"]["elements_read"] == 45
        assert h["deciding_set"] == "CNS_PRIMARY"
        assert h["the_reading"] == {
            "reading": h["by_set"]["CNS_PRIMARY"]["reading"],
            "words": h["by_set"]["CNS_PRIMARY"]["words"],
        }
        assert set(h["by_set"]) == {
            "CNS_PRIMARY",
            "BRAIN_ONLY",
            "CNS_AND_CNS_LINES",
            "CNS_AND_ALL_NEURAL_LINES",
        }
        assert h["axis_columns_observed"] == [371]
        assert out["primary_claim"]["refused"] is True

    def test_the_gate_reads_the_ledger_only(self):
        """No cache path is a parameter of the gate, so it cannot open an answer."""
        import inspect

        assert list(inspect.signature(astropost.check_the_run_is_complete).parameters) == ["ledger"]


# ------------------------------------------------------------------------------------------ item (h)


#: 23 chromosome sizes spread like the registered design's (5 to 137 requests). Equal sizes give a
#: resampled ratio few distinct values -- a sum of whole counts over a fixed total -- which the
#: registered rule rightly calls degenerate; the real clusters are unequal, and so are these.
VARIED = [5, 137, 40, 66, 23, 91, 12, 58, 77, 31, 104, 47, 19, 85, 62, 28, 120, 36, 70, 15, 53, 99, 44]


def per_chrom(rates: list[float], size, seed: int = 1):
    """`in_set` and `chrom_of` with each chromosome's share of hits set by `rates`."""
    rng = random.Random(seed)
    in_set, chrom_of = {}, {}
    sizes = size if isinstance(size, list) else [size] * len(rates)
    for ci, (rate, size) in enumerate(zip(rates, sizes, strict=True)):
        hits = round(rate * size)
        flags = [True] * hits + [False] * (size - hits)
        rng.shuffle(flags)
        for j, f in enumerate(flags):
            e = f"c{ci}_{j}"
            in_set[e], chrom_of[e] = f, f"chr{ci + 1}"
    return in_set, chrom_of


class TestItemHIsTheRegisteredStatistic:
    def test_the_registration_is_read_and_hash_checked(self):
        terms = astropost.item_h_registration()
        assert terms["sha256"] == "941fe4d18ad5aa8516aa6f29375d3b15e883d589930faf5327bef50eda416dcc"
        assert terms["controls"]["CNS_PRIMARY"] == 0.07623431
        assert len(terms["sets"]["CNS_PRIMARY"]) == 27
        assert "CL:0000127" in terms["sets"]["CNS_PRIMARY"]
        assert terms["bands"]["bootstrap_seed"] == 20261003
        assert terms["bands"]["bootstrap_draws"] == 2000

    def test_PLANTED_a_changed_registration_REFUSES(self, tmp_path):
        for rel in (astropost.ASTROARGMAX_REGISTRATION, astropost.ASTROARGMAX_AMENDMENT_1):
            (tmp_path / rel).parent.mkdir(parents=True, exist_ok=True)
            (tmp_path / rel).write_bytes((ROOT / rel).read_bytes())
        astropost.item_h_registration(tmp_path)  # the copy is the registration and passes
        reg = tmp_path / astropost.ASTROARGMAX_REGISTRATION
        reg.write_bytes(
            reg.read_bytes().replace(b'"ratio_enriched_floor": 2.0', b'"ratio_enriched_floor": 1.5')
        )
        with pytest.raises(astropost.AnalysisRefusedError, match="registers nothing"):
            astropost.item_h_registration(tmp_path)

    def test_enriched_needs_the_ratio_AND_the_level(self):
        rates = [0.35, 0.4, 0.45, 0.5, 0.38, 0.42, 0.47, 0.36, 0.44, 0.41, 0.39, 0.46] * 2
        in_set, chrom_of = per_chrom(rates[:23], VARIED)
        out = astropost.rate_against_control(in_set, chrom_of, 0.07623431, BANDS)
        assert out["bootstrap"]["degenerate"] is False
        assert out["reading"] == "enriched", out
        assert out["words"] == BANDS["readings"]["enriched"]
        assert out["decided_by"] == "the cluster bootstrap"

    def test_PLANTED_a_ratio_without_the_level_is_NOT_enriched(self):
        """A 3x excess at a 6% level: the band's own example of a ratio that misleads."""
        rates = [0.05, 0.075] * 12
        in_set, chrom_of = per_chrom(rates[:23], VARIED)
        out = astropost.rate_against_control(in_set, chrom_of, 0.02, BANDS)
        assert out["bootstrap"]["degenerate"] is False
        assert out["observed_rate"] < BANDS["absolute_level_floor"]
        assert out["bootstrap"]["ci95_ratio"][0] >= BANDS["ratio_enriched_floor"]
        assert out["reading"] != "enriched"

    def test_depleted(self):
        in_set, chrom_of = per_chrom(([0.0, 0.03] * 12)[:23], VARIED)
        out = astropost.rate_against_control(in_set, chrom_of, 0.07623431, BANDS)
        assert out["bootstrap"]["degenerate"] is False
        assert out["reading"] == "depleted", out
        assert out["bootstrap"]["identical_resample_share"] < 1

    def test_PLANTED_a_degenerate_bootstrap_goes_to_the_named_exact_bound(self):
        """Every chromosome at the same share: every resample identical, so the bootstrap decides nothing."""
        in_set, chrom_of = per_chrom([0.075] * 23, 200)
        out = astropost.rate_against_control(in_set, chrom_of, 0.075, BANDS)
        assert out["bootstrap"]["degenerate"] is True
        assert out["bootstrap"]["distinct_resamples"] == 1
        assert out["bootstrap"]["identical_resample_share"] == 1.0
        assert "Wilson" in out["exact_bound"]["method"]
        assert "NOT the cluster bootstrap's" in out["exact_bound"]["not_the_cluster_bootstrap"]
        assert out["decided_by"].startswith("the exact bound alone")
        lo, hi = out["exact_bound"]["ratio_interval"]
        assert 0.8 <= lo < 1 < hi <= 1.25
        assert out["reading"] == "equivalent"

    def test_all_in_the_set_is_degenerate_too(self):
        in_set, chrom_of = per_chrom([1.0] * 12, 10)
        out = astropost.rate_against_control(in_set, chrom_of, 0.07623431, BANDS)
        assert out["bootstrap"]["degenerate"] is True
        assert out["decided_by"].startswith("the exact bound alone")

    def test_PLANTED_fewer_than_ten_clusters_is_inconclusive_and_says_why(self):
        in_set, chrom_of = per_chrom([0.5] * 5, 40)
        out = astropost.rate_against_control(in_set, chrom_of, 0.07623431, BANDS)
        assert out["clusters"] == 5
        assert out["reading"] == "inconclusive"
        assert "THE DATA CANNOT TELL" in out["words"]
        assert "fixes no rule for grouping them" in out["why_no_interval"]
        assert "bootstrap" not in out

    def test_the_readings_are_the_registered_words(self):
        assert astropost.reading_h(2.0, 3.0, 0.1, BANDS)["reading"] == "enriched"
        assert astropost.reading_h(2.0, 3.0, 0.0999, BANDS)["reading"] == "inconclusive"
        assert astropost.reading_h(0.3, 0.8, 0.02, BANDS)["reading"] == "depleted"
        assert astropost.reading_h(0.8, 1.25, 0.08, BANDS)["reading"] == "equivalent"
        assert astropost.reading_h(0.79, 1.25, 0.08, BANDS)["reading"] == "inconclusive"
        for key in ("enriched", "depleted", "equivalent", "inconclusive"):
            assert BANDS["readings"][key] in {
                astropost.reading_h(lo, hi, r, BANDS)["words"]
                for lo, hi, r in ((2, 3, 0.1), (0.3, 0.8, 0.02), (0.8, 1.25, 0.08), (0.5, 3, 0.1))
            }

    def test_one_label_per_element_joined_through_the_census(self, monkeypatch):
        tissues = ["astrocyte", "K562", "brain", "Not_A_Census_Label", None, "neural cell"] * 4
        answers, chrom_of = {}, {}
        for i, t in enumerate(tissues):
            e = f"E{i}"
            answers[e] = synthetic_answer(e, CHROMS[i % 12], t)
            chrom_of[e] = CHROMS[i % 12]
        out = astropost.item_h(answers, chrom_of)
        pop = out["population"]
        assert pop["elements_read"] == 24
        assert pop["elements_with_no_target"] == 4, "counted beside the denominator, never dropped"
        assert pop["elements_with_a_target"] == 20
        assert out["labels_the_census_cannot_join"] == ["Not_A_Census_Label"]
        primary = out["by_set"]["CNS_PRIMARY"]
        # astrocyte, brain and neural cell are in CNS_PRIMARY: 12 of the 20 with a target
        assert (primary["in_the_set"], primary["elements"]) == (12, 20)
        assert set(out["not_computed_here"]) == {
            "comparison_b_the_K562_arm",
            "the_667_column_and_name_collapsed_split",
            "secondary_all_gene_rows",
        }


# ---------------------------------------------------------------------------------- the primary claim


class TestThePrimaryClaimsTailAndItsRefusedAssembly:
    def scored(self, chroms: int, per: int, seed: int = 3):
        rng = random.Random(seed)
        pairs, labels, with_d, without = [], [], [], []
        for c in range(chroms):
            for j in range(per):
                pos = j % 7 == 0
                pairs.append(
                    {
                        "chrom": f"chr{c + 1}",
                        "start": j * 2_000_000,
                        "end": j * 2_000_000 + 500,
                        "gene": f"G{c}_{j}",
                    }
                )
                labels.append(pos)
                base = rng.random()
                without.append(base + (0.2 if pos else 0.0))
                with_d.append(base + (0.8 if pos else 0.0))
        return with_d, without, labels, pairs

    def test_the_tail_quotes_the_wider_interval_and_the_registered_reading(self):
        with_d, without, labels, pairs = self.scored(12, 21)
        out = astropost.primary_tail(with_d, without, labels, pairs)
        assert out["gain"] > 0
        assert out["quoted_interval"] in (out["locus_interval"], out["chromosome_interval"])
        spans = [i["ci95"][1] - i["ci95"][0] for i in (out["locus_interval"], out["chromosome_interval"])]
        quoted = out["quoted_interval"]["ci95"]
        assert quoted[1] - quoted[0] == max(spans)
        assert out["reading"] == "the deletion gain is detected in cultured fetal astrocytes"
        assert out["required_quote"].startswith("exploratory (prevalence-matched power 0.79")

    def test_too_few_clusters_has_no_interval(self):
        with_d, without, labels, pairs = self.scored(4, 1)
        out = astropost.primary_tail(with_d, without, labels, pairs)
        assert "ci95" not in out["chromosome_interval"]
        assert out["reading"].startswith("interval unreliable")

    def test_the_assembly_is_refused_with_each_gap_named(self):
        out = astropost.primary_assembly()
        assert out["refused"] is True
        assert set(out["gaps"]) == {"astrocyte_deletion_value", "distance", "coverage_controls"}

    def test_PLANTED_the_registered_path_finds_no_astrocyte_value_on_a_recorded_answer(self, tmp_path):
        """The first gap shown rather than asserted: the recording path's own aggregate, then the
        registered reader, on an answer whose astrocyte track moved by -0.7."""
        from genomeos.attribution import crispri
        from genomeos.predict import alphagenome_adapter, enhancer_target

        rows = enhancer_target.aggregate([("AAA", "astrocyte", -0.7), ("AAA", "K562", -0.1)])
        d = tmp_path / "chr1"
        d.mkdir()
        (d / "E1.json").write_text(json.dumps({"id": "E1", "genes": rows}))
        cache = crispri.ElementCache(tmp_path)
        assert cache.value("chr1", "E1", "AAA", "K562") == -0.1
        assert cache.value("chr1", "E1", "AAA", "astrocyte") is None
        assert rows[0]["max_drop_tissue"] == "astrocyte", "the value was there, and it was not kept"
        assert "astrocyte" not in enhancer_target.CELLS
        assert "astrocyte" not in alphagenome_adapter.CELL_TRACKS
