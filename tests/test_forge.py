"""Phase 5: BioForge finds a design meeting a target under constraints, labelled predicted."""

from genomeos.forge import Knob, forge, period_of
from genomeos.ir import EvidenceKind
from genomeos.lang import parse_file
from genomeos.runtime.grn import NetworkRuntime


def test_forge_retunes_the_repressilator_period():
    m = parse_file("data/demo/repressilator.bio")
    base = NetworkRuntime(m).run(hours=60, dt=0.02, initial={"TetR": 10}, record_every=5)
    p0 = period_of(base, "TetR")
    assert p0 is not None
    target = p0 * 1.8
    knobs = [
        Knob("param:mrna_half_life", 0.2, 3.0),
        Knob("param:protein_half_life", 0.05, 1.0),
        Knob("gene:tetR.max", 20, 1000),
        Knob("gene:lacI.max", 20, 1000),
        Knob("gene:cI.max", 20, 1000),
    ]

    def loss(traj):
        p = period_of(traj, "TetR")
        return abs(p - target) / target if p else 10.0

    def still_oscillates(traj):
        return traj.peaks("TetR") >= 3

    res = forge(
        m,
        knobs,
        loss,
        [still_oscillates],
        hours=120,
        dt=0.02,
        initial={"TetR": 10},
        iterations=40,
        restarts=2,
        seed=1,
    )
    assert res.best.feasible and res.best.loss < 0.15, (res.best.loss, res.changes())
    assert res.changes()  # something was tuned
    changed_params = [
        p for p in res.best.module.parameters.values() if p.evidence.kind is EvidenceKind.PREDICTED
    ]
    changed_genes = [g for g in res.best.module.genes() if g.evidence.kind is EvidenceKind.PREDICTED]
    assert changed_params or changed_genes
    assert all(p.confidence <= 0.3 for p in changed_params) and all(
        g.confidence <= 0.3 for g in changed_genes
    )
    assert res.evaluations > 40 and res.history[-1] <= res.history[0]
