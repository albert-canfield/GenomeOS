"""Experiments: knockouts against the wild type reproduce the classic C. elegans founder phenotypes."""

from genomeos.lang import parse, parse_file
from genomeos.organism.experiment import run_experiment, run_experiments
from genomeos.runtime.body import Body


def test_absent_and_knockouts_in_a_toy():
    src = """
module toy.ko
import bio.std.development
organism T { root: Z; cell_type: Zygote; factors: F }
cell_type WithF { parent: Blastomere }
cell_type WithoutF { parent: Blastomere }
signal S { mode: contact; ligand: L; receptor: R; from: cell = B; to: cell = A; sets: Sig = received }
cell_type Signalled { parent: Blastomere }
decision z { action: divide; when: cell = Z; daughters: A, B; after: 1 min; asymmetric: F -> A }
decision with { action: differentiate; when: cell = A, F = present; to: WithF }
decision without { action: differentiate; when: cell = A, F = absent; to: WithoutF }
decision sig { action: differentiate; when: cell = B, Sig = absent; to: Signalled }
experiment noF {
  knockout: F; until: 10 min; expect: "A takes the F-less fate"; assert: type WithoutF at 5 min = 1
}
experiment noL { knockout: L; until: 10 min }
"""
    m = parse(src)
    wt = Body(m).run(until=10)
    assert wt.cells["A"].cell_type == "WithF" and wt.cells["B"].cell_type == "Blastomere"
    r_f, r_l = run_experiments(m)
    assert r_f.mutant.cells["A"].cell_type == "WithoutF" and [c.cell for c in r_f.founders] == ["A"]
    assert r_f.founders[0].lost == ["with"] and r_f.founders[0].gained == ["without"]
    assert r_f.asserts()[0]["ok"] and "A takes the F-less fate" in r_f.format()
    # knocking out the ligand silences the signal: B never receives it and takes the Sig-absent fate
    assert (
        r_l.mutant.cells["B"].cell_type == "Signalled" and r_l.wild_type.cells["B"].cell_type == "Blastomere"
    )
    d = r_f.to_dict()
    assert d["knockouts"] == ["F"] and d["changed_cells"] == 1 and d["evidence"]["kind"] == "none"


def test_celegans_mutants_reproduce_the_published_founder_phenotypes():
    m = parse_file("data/organisms/celegans/mutants.bio")
    assert [e.name for e in m.experiments] == ["pop1", "skn1", "pie1", "apx1", "glp1", "pal1"]
    by = {e.name: e for e in m.experiments}

    def types(r, _t=None):
        return {
            c: r.mutant.cells[c].cell_type for c in ("ABa", "ABp", "EMS", "MS", "E", "C", "D", "P2", "P4")
        }

    wt = Body(m, means=True).run(until=800)
    assert wt.cells["MS"].cell_type == "MSPrecursor" and wt.cells["E"].cell_type == "EPrecursor"
    assert wt.cells["ABp"].cell_type == "ABpPrecursor" and wt.cells["C"].cell_type == "CPrecursor"
    pop1 = run_experiment(m, by["pop1"])
    assert types(pop1, 0)["MS"] == "EPrecursor" and [f.cell for f in pop1.founders] == ["MS"]
    assert all(a["ok"] for a in pop1.asserts()), pop1.asserts()
    skn1 = run_experiment(m, by["skn1"])
    t = types(skn1, 0)
    assert t["MS"] == "CPrecursor" and t["E"] == "CPrecursor" and t["C"] == "CPrecursor"
    assert all(a["ok"] for a in skn1.asserts()), skn1.asserts()
    pie1 = run_experiment(m, by["pie1"])
    assert types(pie1, 0)["P2"] == "EMSPrecursor" and types(pie1, 0)["P4"] != "GermCell"
    assert all(a["ok"] for a in pie1.asserts()), pie1.asserts()
    for name in ("apx1", "glp1"):
        r = run_experiment(m, by[name])
        assert types(r, 0)["ABp"] == "ABaPrecursor" and r.founders[0].cell == "ABp"
        assert r.descendants_changed("ABp") > 500 and all(a["ok"] for a in r.asserts())
    pal1 = run_experiment(m, by["pal1"])
    assert {f.cell for f in pal1.founders} == {"C", "D"} and all(a["ok"] for a in pal1.asserts())
    # the wild type is untouched by the experiments and terminal fates stay the observed ones
    assert pop1.wild_type.count_at(800) == pal1.mutant.count_at(800) == 610
    assert "MS " in pop1.format() and "192 descendants" in pop1.format()
