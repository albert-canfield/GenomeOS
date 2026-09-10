"""Task 3.3: BioTwin fork / run / diff, and timer-gene variants shift the trajectory."""

from genomeos.twin import Twin, diff_runs, timer_effects_for_variants
from genomeos.twin.twin import MeasuredState


def test_variant_rules():
    from genomeos.twin.twin import Variant

    mods = timer_effects_for_variants([Variant("TERT", "nonsense")])
    assert mods["telomerase_compensation"] == 0.0 and mods["telomere_loss_per_division_bp"] > 1
    assert timer_effects_for_variants([Variant("TERT", "missense_variant")]) == {}
    assert "mutations_per_year" in timer_effects_for_variants([Variant("BRCA1", "frameshift_variant")])


def test_fork_run_diff_from_measured_state(tmp_path):
    base = Twin("HG002", measured=MeasuredState(chronological_age=45, telomere_bp=7200, epigenetic_age=47))
    p = base.save(tmp_path)
    loaded = Twin.load(p)
    assert loaded.measured.telomere_bp == 7200
    tert = loaded.fork("HG002_TERT_LoF", note="hypothetical TERT loss of function")
    applied = tert.add_variant("TERT", "nonsense", "p.R100*")
    assert "telomerase_compensation" in applied
    a = loaded.run("hematopoietic_stem", years=30, cells=300, seed=3)
    b = tert.run("hematopoietic_stem", years=30, cells=300, seed=3)
    assert a.reports[0].age_years == 45 and abs(a.reports[0].mean_telomere_bp - 7200) < 150
    assert abs(a.reports[0].mean_epigenetic_age - 47) < 1e-9
    d = diff_runs(a, b)
    assert d["telomere_bp"]["delta"] < -300  # TERT loss: faster attrition in stem cells
    assert d["senescent_fraction"]["delta"] > 0
    assert b.uncertainty.levels["cellular"].confidence < a.uncertainty.levels["cellular"].confidence
    assert b.uncertainty.levels["cellular"].weakest


def test_p53_loss_reduces_senescence():
    base = Twin("X", measured=MeasuredState(chronological_age=50, telomere_bp=6500))
    p53 = base.fork("X_TP53")
    p53.add_variant("TP53", "frameshift_variant")
    a = base.run("fibroblast", years=30, cells=400, seed=5)
    b = p53.run("fibroblast", years=30, cells=400, seed=5)
    assert b.final.senescent_fraction < a.final.senescent_fraction
