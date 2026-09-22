# SPDX-License-Identifier: AGPL-3.0-or-later
"""The BioLang v0.4 §7.2a gate: the published plasticity series, and what it costs to remove a construct.

Grows data/organisms/celegans/plasticity.bio under each experiment, then re-runs the same arms with
`competence` stripped from the program and with `commitment` stripped, so the table says what each
construct is worth rather than that the program passes. Pass means: the four published outcomes hold
with both constructs; removing `competence` loses the closed window (a factor forced after it still
converts the embryo); removing `commitment` loses the protection of cells that have already
differentiated. If a construct can be removed with no change to any arm, it is not carrying the result
and this gate says so.

    uv run python scripts/plasticity_gate.py
"""

import copy

from genomeos.lang import parse_file
from genomeos.results import save_result
from genomeos.runtime.body import Body, evaluate_assert

PATH = "data/organisms/celegans/plasticity.bio"


def arm(module, ex):
    body = Body(
        module,
        seed=None,
        knockouts=set(ex.knockouts),
        adds=set(ex.adds),
        add_at=dict(ex.add_at),
        environment=ex.environment,
    ).run(until=ex.until)
    types: dict[str, float] = {}
    for c in body.alive_at(ex.until):
        types[c.cell_type] = types.get(c.cell_type, 0.0) + c.count
    checks = [evaluate_assert(body, a) for a in ex.asserts]
    summary = body.summary()
    return {
        "types": dict(sorted(types.items(), key=lambda kv: -kv[1])),
        "muscle": types.get("Muscle", 0.0),
        "asserts": [{"assert": c["assert"], "value": c.get("value"), "ok": bool(c["ok"])} for c in checks],
        "passed": all(c["ok"] for c in checks),
        "outside_competence": summary["outside_competence"],
        "refused_committed": summary["refused_committed"],
        "revised_fates": summary["revised_fates"],
    }


def run(module, label):
    return {ex.name: arm(module, ex) for ex in module.experiments} | {"program": label}


def main() -> None:
    module = parse_file(PATH)
    no_competence = copy.deepcopy(module)
    no_competence.competences.clear()
    for d in no_competence.decisions:
        d.competence = ""
    no_commitment = copy.deepcopy(module)
    no_commitment.commitments.clear()

    conditions = {
        "as written": run(module, "competence and commitment"),
        "competence removed": run(no_competence, "commitment only"),
        "commitment removed": run(no_commitment, "competence only"),
    }
    names = [ex.name for ex in module.experiments]
    width = max(len(n) for n in names)
    print(f"{'arm':<{width}}  " + "  ".join(f"{c:>20}" for c in conditions))
    for name in names:
        cells = []
        for cond in conditions.values():
            a = cond[name]
            cells.append(f"{a['muscle']:>7.0f} muscle {'ok ' if a['passed'] else 'FAIL'}")
        print(f"{name:<{width}}  " + "  ".join(f"{c:>20}" for c in cells))
    verdict = {
        "all arms pass as written": all(conditions["as written"][n]["passed"] for n in names),
        "competence is load-bearing": any(
            conditions["as written"][n]["passed"] != conditions["competence removed"][n]["passed"]
            for n in names
        ),
        "commitment is load-bearing": any(
            conditions["as written"][n]["passed"] != conditions["commitment removed"][n]["passed"]
            for n in names
        ),
    }
    for k, v in verdict.items():
        print(f"{k}: {v}")
    path = save_result(
        "plasticity_gate",
        {
            "program": PATH,
            "series": "Fukushige & Krause 2005, Development 132:1795; Yuzyuk et al. 2009, Dev Cell 16:699",
            "conditions": conditions,
            "verdict": verdict,
        },
    )
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
