# SPDX-License-Identifier: AGPL-3.0-or-later
"""Many human genomes at once: the MAF distiller, events, value domains and the classes."""

import gzip

import pytest

from genomeos.attribution import human_panel as hp

MAF = b"""##maf version=1

a score=0.0
s hg38.chrT          100 10 + 1000 ACGTACGTAC
s hs1.chrT           50 10 + 900 ACGTACGTAC
i hs1.chrT           C 0 C 0
s GCA_000000001.1.c1 10 10 + 800 ACTTACGTAC
i GCA_000000001.1.c1 C 0 C 0
s GCA_000000002.1.c7 20 8 - 700 ACG--CGTAC
i GCA_000000002.1.c7 I 5 C 0
s GCA_000000003.1.c2 30 10 + 600 ACTTACGTAC

a score=0.0
s hg38.chrT          110 6 + 1000 GG--ATCC
s hs1.chrT           60 8 + 900 GGTTATCC
s GCA_000000001.1.c1 20 8 + 800 GGTTATCC
s GCA_000000002.1.c7 27 6 - 700 GG--ATCC
e GCA_000000003.1.c2 40 0 + 600 C

"""


@pytest.fixture()
def panel(tmp_path):
    d = hp.Distiller("chrT")
    for para in hp.split_blocks(MAF):
        d.add(hp.parse_block(para))
    d.save(tmp_path / "chrT")
    return hp.Panel(tmp_path / "chrT")


def test_assembly_names():
    assert hp.assembly_of("GCA_018472595.1.JAHBCB010000033.1") == "GCA_018472595.1"
    assert hp.assembly_of("hs1.chr21") == "hs1"


def test_distiller_keeps_blocks_sites_and_insertions(panel):
    assert panel.n == 4
    assert panel.block_start == [100, 110] and panel.block_end == [110, 116]
    # the second block: three assemblies present, one deleted with contiguous flanks
    assert panel.block_status[1].count(hp.PRESENT) == 3 and panel.block_status[1].count(hp.DELETED) == 1
    by_pos = dict(zip(panel.site_pos, panel.site_alleles, strict=True))
    assert by_pos[102] == b"=T=T"  # G>T in two assemblies
    assert by_pos[103] == b"==-=" and by_pos[104] == b"==-="  # a two-base deletion in one
    # two bases inserted after hg38 position 111 by the assemblies that carry TT
    assert panel.ins_aligned[111] == {0: 2, 1: 2}
    assert panel.ins_unaligned[99] == {2: 5}
    assert panel.gc[1] == (5 + 4, 10 + 6)


def test_events_collapse_deletion_runs(panel):
    evs = hp.events(panel, [(100, 116)])
    kinds = [(e.kind, e.start, e.end) for e in evs]
    assert ("snv", 102, 103) in kinds
    assert ("deletion", 103, 105) in kinds
    assert ("structural", 110, 116) in kinds
    ins = [e for e in evs if e.kind == "insertion" and e.start == 111][0]
    assert (
        ins.tokens == [2, 2, 0, "*"] and ins.minor == 1
    )  # the absent assembly is its own state, not a minor allele
    snv = [e for e in evs if e.kind == "snv"][0]
    assert snv.minor == 2 and snv.informative == 4


def test_measure_presence_identity_and_structure(panel):
    m = hp.measure(panel, [(100, 110)])
    assert m["bases"] == 10 and m["aligned_share"] == 1.0
    assert m["presence"] == pytest.approx(38 / 40)
    assert m["identity"] == pytest.approx(1 - 2 / 38)
    assert m["events"]["snv"] == {"events": 1, "recurring": 1}
    assert m["recurring_events"] == 1
    m2 = hp.measure(panel, [(110, 116)])
    assert m2["presence"] == pytest.approx(18 / 24) and m2["deleted_share"] == pytest.approx(6 / 24)


def test_domain_counts_values(panel):
    d = hp.domain(panel, [(100, 110)])
    assert d["recurring_events"] == 1 and d["values"] == 2 and d["recurring_values"] == 2
    assert d["top_share"] == 0.5 and d["effective_values"] == 2.0
    assert d["class"] == "storage"
    described = hp.describe_values(d)
    assert {tuple(v["differs_from_hg38"]) for v in described} == {("hg38",), ("103>T",)}


def _stats(**kw):
    base = {
        "aligned_share": 1.0,
        "informative": 89,
        "top_share": 1.0,
        "recurring_cover": 1.0,
        "recurring_values": 1,
    }
    base.update(kw)
    return base


def test_unit_class_rules():
    assert hp.unit_class(_stats(), 89) == "fixed"
    assert hp.unit_class(_stats(top_share=0.6, recurring_values=3, recurring_cover=0.95), 89) == "storage"
    assert (
        hp.unit_class(_stats(top_share=0.2, recurring_values=20, recurring_cover=0.95), 89) == "hypervariable"
    )
    assert (
        hp.unit_class(_stats(top_share=0.2, recurring_values=3, recurring_cover=0.4), 89) == "hypervariable"
    )
    assert hp.unit_class(_stats(informative=60), 89) == "unplaced"
    assert hp.unit_class(_stats(aligned_share=0.3), 89) == "unplaced"


def test_tile_joins_short_tail():
    assert hp.tile([(0, 450)], 200) == [[(0, 200)], [(200, 400), (400, 450)]]
    assert hp.tile([(0, 150), (300, 400)], 200) == [[(0, 150), (300, 350), (350, 400)]]
    assert hp.tile([(0, 520)], 200) == [[(0, 200)], [(200, 400)], [(400, 520)]]


def test_byte_runs_join_consecutive_blocks_and_skip_unwanted():
    index = [(0, 10, 0), (10, 20, 100), (20, 30, 250), (500, 510, 400)]
    assert hp.byte_runs(index, None, 600) == [(0, 600)]
    assert hp.byte_runs(index, [(12, 25)], 600) == [(100, 300)]
    assert hp.byte_runs(index, None, 600, limit=260) == [(0, 250), (250, 150), (400, 200)]


def test_poisson_tails_and_block_classes():
    lo, hi = hp.poisson_tails(0, 10.0)
    assert lo == pytest.approx(4.54e-5, rel=1e-2) and hi == pytest.approx(1.0)
    base = {"aligned_share": 1.0, "missing_share": 0.0, "presence": 0.999, "touched_recurring_share": 0.01}
    assert hp.block_class(base, 2, 20.0)["class"] == "core"
    assert hp.block_class(base, 18, 20.0)["class"] == "variable"
    assert hp.block_class(base, 0, 2.0)["reason"] == "too short to call"
    assert hp.block_class({**base, "presence": 0.9}, 18, 20.0)["class"] == "polymorphic"
    assert hp.block_class({**base, "presence": 0.3}, 18, 20.0)["class"] == "lineage_restricted"
    assert hp.block_class({**base, "aligned_share": 0.2}, 18, 20.0)["class"] == "unplaced"


def test_lift_point_through_a_chain(tmp_path):
    p = tmp_path / "c.chain.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("chain 100 chrT 1000 + 100 300 chrQ 2000 + 500 700 1\n50 10 10\n140\n\n")
        fh.write("chain 50 chrT 1000 + 400 450 chrR 1000 - 100 150 2\n50\n\n")
    chain = hp.load_chain("chrT", p)
    starts = [c[0] for c in chain]
    assert hp.lift_point(chain, starts, 120) == ("chrQ", 520)
    assert hp.lift_point(chain, starts, 155) is None  # inside the gap
    assert hp.lift_point(chain, starts, 170) == ("chrQ", 570)
    assert hp.lift_point(chain, starts, 410) == ("chrR", 1000 - 110 - 1)
    assert hp.rt_stratum(50.0, [30.0, 55.0]) == 1 and hp.rt_stratum(20.0, [30.0, 55.0]) == 0


def test_background_rate_and_expectation(panel):
    evs = hp.events(panel, [(0, 1000)])
    bg = hp.build_background(panel, [], evs)
    assert bg.overall == pytest.approx(0.0)  # one bin of 16 aligned bases is under half a kilobase
    bg.density = {None: 0.01}
    bg.overall = 0.01
    assert bg.expected(panel, [(100, 116)]) == pytest.approx(0.16)


def test_small_helpers():
    assert hp.entropy_bits([5, 5]) == pytest.approx(1.0)
    assert hp.spearman([1, 2, 3, 4], [10, 20, 30, 40]) == 1.0
    assert hp.merge_intervals([(5, 10), (1, 6), (20, 30)]) == [(1, 10), (20, 30)]
    assert hp.gc_stratum(0.42) == 2 and hp.gc_stratum(None) is None
    assert hp.sv_row_carriers("hprc_v21_sv", {"AC": "17", "alleleNumber": "464", "alleleFreq": "0.0366"}) == (
        17,
        464,
        0.0366,
    )
    assert hp.sv_row_carriers("hprc_inv", {"label": "2"}) == (2, 90, 0.0222)
    assert hp.sv_row_carriers("hprc_del", {"label": "5:5512bp"}) == (5, 90, 0.0556)
    assert hp.Overlaps([{"chromStart": 7, "chromEnd": 7}]).over(7, 8)[0]["chromEnd"] == 8


def test_gnomad_join_keeps_the_two_counts_apart():
    ev = hp.Event("snv", 102, 103, ["=", "T", "=", "T"], "=")
    rows = hp.Overlaps(
        [
            {
                "chromStart": 102,
                "chromEnd": 103,
                "ref": "G",
                "alt": "T",
                "AF": "0.25",
                "AN": "1000",
                "rsId": "rs1",
            }
        ]
    )
    out = hp.gnomad_for_events(rows, [ev], 4)
    assert out[0]["panel"] == {"haplotypes": 4, "minor": 2, "minor_share": 0.5}
    assert out[0]["gnomad"]["af"] == 0.25 and out[0]["gnomad"]["an"] == 1000


def test_event_index_slices_spans(panel):
    idx = hp.EventIndex(hp.events(panel, [(0, 1000)]))
    kinds = {e.kind for e in idx.over([(104, 105)])}
    assert kinds == {"deletion"}
    assert {e.kind for e in idx.over([(112, 113)])} == {"structural"}


def test_length_class_reads_copy_numbers():
    # DRD4-like: 4R in most, 7R in a quarter, 2R in a few
    assert hp.length_class({0: 58, 144: 22, -96: 4, 49: 2, 97: 1}, 89)["class"] == "storage"
    assert hp.length_class({0: 89}, 89)["class"] == "fixed"
    assert hp.length_class({i: 1 for i in range(89)}, 89)["class"] == "hypervariable"


def test_depletion_null_counts_chance_and_dispersion():
    bg = hp.Background(
        bins={0: [1000, 0, 0.4], 1: [1000, 20, 0.4], 2: [1000, 10, 0.4]}, density={2: 0.01}, overall=0.01
    )
    out = hp.depletion_null({"_bg": (bg, bg)})
    assert out["kilobases"] == 3 and out["depleted_observed"] == 1
    assert out["depleted_expected_under_poisson"] == 0  # three bins at 3% each
    assert out["pearson_dispersion"] > 1


class _FakeResponse:
    def __init__(self, status, data, total):
        self.status, self._data = status, data
        self.headers = {"Content-Range": f"bytes 0-{len(data) - 1}/{total}"}

    def read(self):
        return self._data


class _FakeConn:
    def __init__(self, host, plan):
        self.host, self.plan, self.calls = host, plan, 0

    def request(self, method, path, headers=None):
        self.calls += 1
        if self.plan.get(self.host) == "fail":
            raise TimeoutError("timed out")

    def getresponse(self):
        return _FakeResponse(206, b"abcd", self.plan.get(self.host, 100))

    def close(self):
        pass


def test_keepalive_source_retries_on_the_mirror_and_checks_sizes(monkeypatch):
    plan = {"hgdownload.soe.ucsc.edu": "fail", "hgdownload2.soe.ucsc.edu": 100}
    monkeypatch.setattr(hp.time, "sleep", lambda s: None)
    monkeypatch.setattr(hp.http.client, "HTTPSConnection", lambda host, timeout=None: _FakeConn(host, plan))
    hp._POOL.conns = {}
    src = hp.KeepAliveSource("https://hgdownload.soe.ucsc.edu/gbdb/hg38/x.bb")
    assert src.read(0, 4) == b"abcd" and src.mirror_requests == 1 and src.total == 100
    plan["hgdownload2.soe.ucsc.edu"] = 101  # a mirror with another file
    with pytest.raises(OSError):
        src.read(0, 4)
    hp._POOL.conns = {}


def _fake_result(cds_ratio, storage_late, fixed_late, core_constrained, var_constrained):
    """The fields claim_numbers reads, with everything else left out."""
    return {
        "pooled": {
            "cds": {"ratio_gc": 0.35, "ratio_gc_rt": cds_ratio},
            "neutral": {"ratio_gc": 0.96, "ratio_gc_rt": 0.91},
        },
        "units": {
            "cds": {"placed_shares": {"fixed": 0.70}, "expected_gc_rt_matched": {"fixed": 0.44}},
            "neutral": {"placed_shares": {"fixed": 0.48}, "expected_gc_rt_matched": {"fixed": 0.47}},
            "fossil": {"placed_shares": {"fixed": 0.44}, "expected_gc_rt_matched": {"fixed": 0.47}},
            "background": {
                "placed_shares": {"storage": 0.5},
                "classes": {"fixed": {"n": 10}, "storage": {"n": 10}, "hypervariable": {"n": 0}},
            },
        },
        "genes": {"core_share_of_callable": 0.58},
        "matched_windows": {"cds_genes": {"core_share_of_callable": 0.14}},
        "timing_by_unit_class": {"fixed": {"late": fixed_late}, "storage": {"late": storage_late}},
        "against_gnocchi": {
            "per_kilobase": {
                "spearman_z_vs_ratio": -0.08,
                "poisson_null": {"pearson_dispersion": 15.6},
                "gnocchi_unscored_panel_depleted": {"duplicated_share": 0.58},
            },
            "blocks_core_or_variable_by_gnocchi": {
                "core_gnocchi_constrained": core_constrained,
                "core_gnocchi_free": 10 - core_constrained,
                "variable_gnocchi_constrained": var_constrained,
                "variable_gnocchi_free": 100 - var_constrained,
            },
        },
        "catalogues": {"unknown_space": {"fixed": {"bp": 100}, "storage": {"bp": 200}}},
        "storage": {"storage_units": 3, "with_either": 2, "recurring_values_histogram": {"2": 3}},
    }


def test_claim_numbers_reads_the_gaps_and_the_timing_survival():
    n = hp.claim_numbers(_fake_result(0.41, 0.34, 0.33, 3, 30))
    assert n["cds_pooled_ratio_gc_rt"] == 0.41
    assert n["cds_units_fixed_excess"] == pytest.approx(0.26)
    assert n["neutral_units_fixed_excess"] == pytest.approx(0.01)
    assert n["genes_core_excess_over_matched"] == pytest.approx(0.44)
    # the coding-neutral gap is 0.61 before the timing match and 0.50 after
    assert n["timing_match_survival"] == pytest.approx(0.82, abs=0.01)
    assert n["unit_classes_even_over_timing"] == pytest.approx(0.01)
    # core blocks constrained 30% of the time, variable blocks 30%: no agreement
    assert n["core_blocks_gnocchi_agreement"] == pytest.approx(0.0)
    # the worst of the unknown tiers present: fossil is 3 points under its matched share, neutral 1 over
    assert n["no_unknown_tier_above_matched"] == pytest.approx(0.01)


def test_claims_across_names_the_chromosomes_a_claim_fails_on():
    per = {
        "chr21": hp.claim_numbers(_fake_result(0.41, 0.34, 0.33, 3, 30)),
        "chr22": hp.claim_numbers(_fake_result(0.42, 0.37, 0.31, 18, 28)),
        "chr19": hp.claim_numbers(_fake_result(0.44, 0.34, 0.34, 5, 30)),
    }
    out = hp.claims_across(per)
    assert out["cds_pooled_ratio_gc_rt"]["verdict"] == "genome-wide"
    assert out["cds_pooled_ratio_gc_rt"]["spread"] == {"min": 0.41, "median": 0.42, "max": 0.44}
    assert out["unit_classes_even_over_timing"]["fails_on"] == ["chr22"]
    assert out["unit_classes_even_over_timing"]["holds_on"] == 2


def test_genome_wide_pools_the_autosomes_and_keeps_chry_apart(tmp_path, monkeypatch):
    import json as _json

    for chrom, ratio in (("chr21", 0.41), ("chr22", 0.42), ("chrY", 1.66)):
        (tmp_path / f"human_panel_{chrom}.json").write_text(
            _json.dumps(_fake_result(ratio, 0.34, 0.33, 3, 30))
        )
    monkeypatch.setattr("genomeos.results.RESULTS_DIR", tmp_path)
    out = hp.genome_wide(("chr21", "chr22", "chrY"), results_dir=tmp_path)
    assert out["chromosomes_read"] == ["chr21", "chr22", "chrY"]
    assert out["pooled_over"] == ["chr21", "chr22"] and out["not_pooled"] == ["chrY"]
    assert out["claims"]["cds_pooled_ratio_gc_rt"]["chromosomes"] == 2  # chrY's 1.66 is not in it
    assert out["unknown_space_bp"] == {"fixed": 200, "storage": 400}
    assert out["storage"]["storage_units"] == 6 and out["storage"]["with_either"] == 4
    assert out["chrY"]["storage"]["storage_units"] == 3
