"""Task 3.1: epigenetic clocks from published coefficients; real-data check on GEO if downloaded."""

import math
from pathlib import Path

import pytest

from genomeos.twin.clocks import Clock, MethylationMatrix, _inverse_horvath, horvath_transform, pearson

HORVATH = Path("data/knowledge/Horvath1.csv")
HANNUM = Path("data/knowledge/Hannum.csv")
GEO = Path("data/knowledge/GSE41169_series_matrix.txt.gz")


def test_horvath_transform_roundtrip():
    for age in (0.5, 5, 20, 21, 45, 90):
        assert abs(_inverse_horvath(horvath_transform(age)) - age) < 1e-9
    assert _inverse_horvath(0.0) == 20.0


@pytest.mark.skipif(not (HORVATH.exists() and HANNUM.exists()), reason="fetch clock coefficient files")
def test_clocks_load_and_respond_to_methylation():
    h = Clock.horvath()
    assert len(h.coefficients) == 353 and h.intercept != 0.0
    hn = Clock.hannum()
    assert len(hn.coefficients) == 71
    base = dict.fromkeys(h.coefficients, 0.5)
    r0 = h.predict(base)
    assert r0.coverage == 1.0 and 0 < r0.age < 120
    # raising methylation at positively weighted sites must raise predicted age
    older = {c: (0.9 if k > 0 else 0.1) for c, k in h.coefficients.items()}
    assert h.predict(older).age > r0.age
    partial = h.predict({c: 0.5 for c in list(h.coefficients)[:100]})
    assert partial.missing == 253 and partial.parameter.confidence < r0.parameter.confidence


@pytest.mark.skipif(not (GEO.exists() and HORVATH.exists()), reason="GEO series matrix not downloaded")
def test_horvath_on_real_blood_samples_tracks_age():
    h = Clock.horvath()
    m = MethylationMatrix.from_geo_series_matrix(GEO, cpgs=h.coefficients)
    assert len(m.samples) >= 50
    ages, preds = [], []
    for s in m.samples:
        age = m.metadata.get(s, {}).get("age")
        if not age:
            continue
        r = h.predict(m.sample(s))
        if not math.isnan(r.age):
            ages.append(float(age))
            preds.append(r.age)
    assert len(ages) >= 50
    r = pearson(ages, preds)
    mae = sum(abs(a - p) for a, p in zip(ages, preds, strict=True)) / len(ages)
    assert r > 0.8, r  # Horvath reports r≈0.96 on training data; blood test sets are >0.8
    assert mae < 10, mae
