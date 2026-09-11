"""The Body runtime: one cell to an organism by decisions, timers and signals, scored against truth."""

from pathlib import Path

import pytest

from genomeos.lang import parse
from genomeos.organism import REFERENCES, ReferenceLineage
from genomeos.organism.diff import compare
from genomeos.runtime.body import Body, evaluate_assert

TOY = """
module toy.body
import bio.std.development

organism Toy {
  root: Z; cell_type: Zygote; factors: M; tempo: 1.0
  assert: count at 10 min = 2
  assert: count at 100 min in 2..2
  assert: deaths at 100 min = 1
  assert: fate Neuron at 100 min >= 1
  assert: lineage A at 100 min = 1
  evidence: curated "toy"; confidence: 0.5
}
stage Early { from: 0 min; to: 30 min }
stage Late { from: 30 min }
timer fast { duration: 10 min; when: lineage = A; lengthening: 2.0; }
timer slow { duration: 40 min; when: cell = Z; }
signal S { mode: contact; from: cell = B; to: cell = Aa; sets: Sig = received }
cell_type Neuron2 { parent: PostMitotic; }
decision zyg {
  action: divide; when: cell = Z; daughters: A, B; lineages: A = A, B = B; asymmetric: M -> B; after: 5 min
}
decision a_div { action: divide; when: lineage = A, generation = 0; }
decision a_late { action: divide; when: cell = Aa, stage = Late; }
decision aa_neuron { action: differentiate; when: cell = Aa, Sig = received; to: Neuron2; name: N1; }
decision ap_dies { action: die; when: cell = Ap; after: 20 min; }
decision b_germ { action: differentiate; when: cell = B, M = present; to: GermCell; }
decision b_quiet { action: quiesce; when: cell = B; }
"""


def test_toy_body_divides_signals_differentiates_and_dies():
    b = Body(parse(TOY)).run(until=100)
    assert set(b.cells) == {"Z", "A", "B", "Aa", "Ap"}
    z, a, bb, aa, ap = (b.cells[n] for n in ("Z", "A", "B", "Aa", "Ap"))
    assert z.divides_at == 5 and a.lineage == "A" and a.generation == 0 and bb.lineage == "B"
    assert "M" in bb.factors and "M" not in a.factors  # asymmetric inheritance
    assert a.divides_at == 15 and aa.generation == 1 and aa.born == 15  # timer `fast` for lineage A
    assert bb.cell_type == "GermCell" and bb.quiescent and bb.divides_at is None
    assert aa.factors.get("Sig") == "received" and aa.cell_type == "Neuron2" and aa.terminal_name == "N1"
    assert aa.divides_at is None  # stage Late never applied at birth (15 min); no division without a decision
    assert ap.dies_at == 35 and b.deaths_by(100) == 1 and b.count_at(100) == 2  # B and Aa remain
    assert b.fates_at(100) == {"Neuron2": 1, "GermCell": 1}
    checks = b.check_asserts()
    assert [c["ok"] for c in checks] == [True, True, True, False, True], checks
    assert evaluate_assert(b, "nonsense text")["ok"] is False
    rep = b.uncertainty().to_dict()
    assert rep["organism"]["items"] >= 5 and rep["cellular"]["items"] == 1
    assert b.unknown == {} and not any("UNKNOWN" in line for line in b.tree(depth=3))
    assert b.summary()["cells_born"] == 5


def test_tempo_and_seeded_noise():
    src = TOY.replace("tempo: 1.0", "tempo: 2.0")
    b = Body(parse(src)).run(until=100)
    assert b.cells["Z"].divides_at == 10 and b.cells["A"].divides_at == 30
    noisy = TOY.replace("timer fast { duration: 10 min;", "timer fast { duration: 10 min; sd: 3;")
    times = {Body(parse(noisy), seed=s).run(until=100).cells["A"].divides_at for s in range(5)}
    assert len(times) > 1 and all(5 < t < 40 for t in times)


def test_default_naming_and_lengthening():
    src = """
module toy.naming
import bio.std.development
organism T { root: R; cell_type: Zygote }
timer c { duration: 10 min; lengthening: 2.0; when: lineage = R }
decision start { action: divide; when: cell = R; lineages: Ra = R, Rp = R; after: 0 min }
decision go { action: divide; when: lineage = R, generation = <=2 }
"""
    b = Body(parse(src)).run(until=1000)
    assert {"Ra", "Rp", "Raa", "Rap", "Raal", "Raar", "Raala", "Raalp"} <= set(b.cells)  # a/p, l/r, a/p ...
    assert b.cells["Ra"].divides_at == 10 and b.cells["Raa"].divides_at == 10 + 10
    assert b.cells["Raal"].divides_at == 20 + 20  # lengthening applies from generation 2
    assert b.unknown == {}


def test_no_organism_block_is_an_error():
    with pytest.raises(ValueError, match="no organism"):
        Body(parse("module x\nimport bio.std.development\ntimer t { duration: 1 min }"))


REF = Path(REFERENCES["celegans"])


@pytest.mark.skipif(not REF.exists(), reason="run `genomeos data distil --only celegans_lineage`")
def test_celegans_grows_to_the_adult_and_matches_the_reference():
    from genomeos.lang import parse_file

    m = parse_file("data/organisms/celegans/embryo.bio")
    body = Body(m, means=True).run(until=6000)  # timers at their means: the topology and fates exactly
    ref = ReferenceLineage.load()
    d = compare(body, ref)
    assert d.reference_cells == 2183 and d.matched == 2183 and not d.missing and not d.extra
    assert not d.parent_mismatch and d.fate_accuracy == 1.0 and d.deaths_matched == d.deaths_expected == 131
    assert d.timing_median is not None and d.timing_median < 20
    assert body.count_at(6000) == 961 and body.deaths_by(6000) == 131
    # with the measured spread (the program's own seed) the early wave lands nearer the reference
    seeded = Body(m).run(until=8000)
    assert m.organism.seed == 0 and 230 <= seeded.count_at(350) <= 300 and seeded.count_at(8000) == 961
    assert all(c["ok"] for c in seeded.check_asserts()), seeded.check_asserts()
    # mechanism decided the founders: Wnt reached EMS, Notch reached ABp, PIE-1 stayed in the germline
    assert body.cells["EMS"].factors.get("Wnt") == "received" and body.cells["ABp"].factors.get("Notch")
    assert body.cells["P4"].cell_type == "GermCell" and "PIE-1" not in body.cells["AB"].factors
    assert body.cells["MS"].factors.get("POP-1") and "POP-1" not in body.cells["E"].factors
    embryo = compare(Body(m, means=True).run(until=800), ref, until=800)
    assert embryo.deaths_matched == 110 and embryo.fates_checked == 555
    rep = body.uncertainty().to_dict()
    assert rep["organism"]["label"] == "high" and rep["cellular"]["items"] > 500
