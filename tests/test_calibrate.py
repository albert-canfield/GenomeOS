from genomeos.results import load_result
from genomeos.twin.calibrate import calibrate_attrition, expected_telomere


def test_calibration_recovers_measurement(tmp_path, monkeypatch):
    import genomeos.twin.calibrate as c

    monkeypatch.setattr(c, "save_result", lambda name, payload: tmp_path / f"{name}.json")
    cal = calibrate_attrition("hematopoietic_stem", age_years=45, measured_bp=5900, source="test")
    assert abs(cal.predicted_bp - 5900) < 1e-6
    assert cal.parameter.evidence.kind.value == "inferred" and cal.parameter.confidence == 0.4
    assert 0 < cal.fitted_divisions_per_year < 20
    # a default model run with the fitted value lands near the measurement
    assert abs(expected_telomere("hematopoietic_stem", 45, cal.fitted_divisions_per_year) - 5900) < 1e-6


def test_saved_calibration_if_present():
    r = load_result("calibration_hematopoietic_stem_attrition")
    if r is None:
        return
    assert 20 < r["net_bp_per_year"] < 200
