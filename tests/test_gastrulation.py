"""Task 4.3: three germ layers in the right order along the NODAL gradient."""

from genomeos.runtime.gastrulation import run_gastrulation


def test_three_layers_in_order_and_plausible_proportions():
    r = run_gastrulation(cells=60, hours=30, dt=0.05)
    props = r.proportions()
    assert all(props[f] > 0.05 for f in props), props
    # endoderm nearest the NODAL source, ectoderm farthest
    first = r.fates[0]
    last = r.fates[-1]
    assert first == "endoderm" and last == "ectoderm", (first, last)
    order = [f for i, f in enumerate(r.fates) if i == 0 or f != r.fates[i - 1]]
    assert order == ["endoderm", "mesoderm", "ectoderm"], order
    exp = {k: r.module.parameters[f"expected_{k}"].value for k in ("ectoderm", "mesoderm", "endoderm")}
    for k in exp:
        assert abs(props[k] - exp[k]) < 0.25, (
            k,
            props[k],
            exp[k],
        )  # stated as order-of-magnitude, confidence 0.3
    # half of this module's rules are inferred, and the report must say so
    unc = r.uncertainty().to_dict()["molecular"]
    assert unc["label"] == "low" and unc["items"] > 5
    assert any(rule.evidence.kind.value == "inferred" for rule in r.module.rules)
