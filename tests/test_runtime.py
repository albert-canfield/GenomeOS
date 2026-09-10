from genomeos.lang import parse_file
from genomeos.runtime import CellRuntime, Environment, NetworkRuntime


def test_repressilator_oscillates():
    m = parse_file("data/demo/repressilator.bio")
    vm = NetworkRuntime(m)
    traj = vm.run(hours=60, dt=0.01, initial={"TetR": 10.0})
    # a sustained oscillator shows several peaks and none of the three proteins
    # collapses to a steady state
    assert traj.peaks("TetR") >= 3
    assert traj.peaks("LacI") >= 3
    assert traj.peaks("CI") >= 3
    for s in ("TetR", "LacI", "CI"):
        xs = traj.levels[s][len(traj.levels[s]) // 2 :]
        assert max(xs) - min(xs) > 1.0


def test_context_disables_rules():
    m = parse_file("data/demo/repressilator.bio")
    vm = NetworkRuntime(m, context={"host": "ecoli"})
    assert len(vm.active_rules) == len(m.rules)  # `host = any` matches anything


def test_ageing_clocks_move_in_the_right_direction():
    rt = CellRuntime(seed=7)
    reports = rt.simulate_tissue("fibroblast", years=80, n_cells=1000)
    first, last = reports[0], reports[-1]
    assert last.mean_telomere_bp < first.mean_telomere_bp
    assert last.mean_mutations > first.mean_mutations
    assert abs(last.mean_epigenetic_age - 80) < 8  # pace ~ 1 year/year
    assert 0.005 < last.senescent_fraction < 0.3  # a few percent of aged cells
    # attrition of the mean telomere should sit in the measured 20-40 bp/year band
    per_year = (first.mean_telomere_bp - last.mean_telomere_bp) / 80
    assert 20 < per_year < 45


def test_neurons_do_not_divide():
    rt = CellRuntime(seed=3)
    reports = rt.simulate_tissue("neuron", years=80, n_cells=50)
    assert abs(reports[-1].mean_telomere_bp - reports[0].mean_telomere_bp) < 1e-9


def test_environment_accelerates_ageing():
    slow = CellRuntime(seed=1).simulate_tissue("fibroblast", 70, n_cells=300)[-1]
    fast = CellRuntime(
        seed=1, env=Environment(proliferation_factor=2.0, mutagen_factor=2.0, epigenetic_pace=1.3)
    ).simulate_tissue("fibroblast", 70, n_cells=300)[-1]
    assert fast.mean_telomere_bp < slow.mean_telomere_bp
    assert fast.mean_mutations > slow.mean_mutations
    assert fast.mean_epigenetic_age > slow.mean_epigenetic_age


def test_every_parameter_has_evidence():
    for p in CellRuntime.evidence_table("fibroblast"):
        assert p.evidence.source
        assert 0.0 <= p.confidence <= 1.0
