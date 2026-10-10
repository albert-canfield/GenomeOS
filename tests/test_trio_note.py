# SPDX-License-Identifier: AGPL-3.0-or-later
"""The note `trio()` writes into every trio result says what area D's phasing and read allele
fractions found, and each number in it is read back from the committed result file it cites, so the
note cannot drift from the record again. Until 2026-10-10 it said a count above 60 to 100 was
"representation differences and residual no-calls", which the phasing (4d63910) and the read
fractions (0125926) no longer support after normalisation (docs/DATA.md, 2026-09-27 and 2026-09-28)."""

import json
import re
from pathlib import Path

from genomeos.genome import individuals as ind

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "data" / "results"
NOTE = ind.TRIO_NOTE


def _load(name):
    return json.loads((RESULTS / name).read_text())


def test_the_old_reading_is_gone_and_the_cited_files_exist():
    assert "so a larger count is representation differences" not in NOTE
    assert "residual no-calls" not in NOTE
    cited = re.findall(r"data/results/[\w.]+\.json", NOTE)
    assert sorted(set(cited)) == [
        "data/results/trio_q100_phase.json",
        "data/results/trio_q100_population.json",
        "data/results/trio_vaf.json",
    ]
    for path in cited:
        assert (REPO / path).is_file(), path
    for sha in ("4d63910", "0125926"):
        assert sha in NOTE


def test_the_phasing_numbers_are_the_phasing_result():
    d = _load("trio_q100_phase.json")
    c = d["candidates"]
    assessable = c["total"] - c["outside_assembly"]
    one = c["father"] + c["mother"]
    assert f"the {c['total']:,} candidates" in NOTE
    assert f"of {assessable:,} assessable, {one:,} sit on one parental chromosome" in NOTE
    assert f"{c['father']:,} paternal to {c['mother']:,} maternal" in NOTE
    assert one == round(c["father"] / d["candidates_paternal_share"])


def test_the_population_numbers_are_the_population_result():
    groups = _load("trio_q100_population.json")["groups"]
    one_parent = [v for k, v in groups.items() if k.startswith(("cand_father_", "cand_mother_"))]
    snvs = sum(g["snvs"] for g in one_parent)
    seen = sum(g["seen_in_panel"] for g in one_parent)
    ctrl = groups["control_one_parent"]
    by_type = _load("trio_q100_phase.json")["candidates_by_type"]["snv"]
    assert snvs == by_type["father"] + by_type["mother"]
    assert f"{seen} of those {snvs:,} SNVs ({seen / snvs:.1%})" in NOTE
    assert f"{ctrl['seen_in_panel']:,} of {ctrl['snvs']:,} inherited calls" in NOTE
    assert f"({ctrl['seen_in_panel'] / ctrl['snvs']:.1%}; exploratory" in NOTE


def test_the_read_fraction_numbers_are_the_read_fraction_result():
    d = _load("trio_vaf.json")
    assert d["verdicts"]["void"] is False
    assert d["candidates_total"] == _load("trio_q100_phase.json")["candidates"]["total"]
    cand = d["groups"]["candidates_one_parent"]
    ctrl = d["groups"]["control_snv"]
    assert (
        f"the {cand['assessed']:,} assessed sit at a median read allele fraction of {cand['median_vaf']}"
        in NOTE
    )
    assert f"against {ctrl['median_vaf']} for {ctrl['assessed']:,} inherited heterozygous SNVs" in NOTE
    # a heterozygous mutation in a fraction c of the cells sits at a read fraction c / 2
    assert f"about {round(200 * cand['median_vaf'])}% of cells" in NOTE


def test_the_germline_bound_is_the_upper_tail_over_the_controls_share_there():
    d = _load("trio_vaf.json")
    edges = d["histogram_edges"]
    i = edges.index(0.5)
    cand = d["groups"]["candidates_one_parent"]["vaf_histogram"]
    ctrl = d["groups"]["control"]["vaf_histogram"]
    assert sum(cand[edges.index(0.6) :]) == 0  # the candidates' upper tail ends below 0.60
    upper = sum(cand[i:])
    share = (ctrl[i] + ctrl[i + 1]) / sum(ctrl)  # the control's share at 0.50 to 0.60 (42.7%)
    assert f"{upper} sit at 0.50 or above, so at most about {int(upper / share)} can be germline" in NOTE


def test_trio_writes_the_note(tmp_path):
    root = tmp_path / "individuals"
    for name, rows in {
        "kid": [(100, "A", "G", "0/1"), (300, "G", "A", "0/1")],
        "dad": [(100, "A", "G", "0/1")],
        "mum": [(200, "C", "T", "1/1")],
    }.items():
        src = tmp_path / f"{name}.vcf"
        src.write_text(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tX\n"
            + "".join(f"chr21\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\tGT\t{g}\n" for p, r, a, g in rows)
        )
        ind.import_vcf(src, name, root=root)
    r = ind.trio("kid", "dad", "mum", root=root)
    assert r["totals"]["de_novo_candidates"] == 1
    assert r["note"] == NOTE
    assert json.loads((root / "kid" / "trio_dad_mum.json").read_text())["note"] == NOTE
