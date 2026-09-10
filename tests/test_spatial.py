"""Task 4.1: a gradient-reading rule produces a French-flag pattern."""

from genomeos.runtime.spatial import Field2D, french_flag, gradient_decay_length


def test_field_diffuses_and_conserves_without_decay():
    f = Field2D("m", 21, 1, diffusion=0.2, decay=0.0)
    f.add(10, 0, 100.0)
    for _ in range(50):
        f.step(1.0)
    assert abs(sum(f.profile_x(0)) - 100.0) < 1e-6
    assert f.at(10, 0) < 100 and f.at(5, 0) > 0


def test_french_flag_forms_three_bands_in_order():
    rt = french_flag(width=60, height=6, hours=150)
    bands = rt.bands_along_x()
    types = [b[0] for b in bands]
    assert types == ["blue", "white", "red"], bands
    census = rt.census()
    assert all(census[t] > 0 for t in ("blue", "white", "red"))
    # the gradient is monotonic along x
    prof = rt.fields["morphogen"].profile_x()
    assert all(prof[i] >= prof[i + 1] for i in range(len(prof) - 1))
    assert 3 < gradient_decay_length(0.4, 0.02) < 6
