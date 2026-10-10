# SPDX-License-Identifier: AGPL-3.0-or-later
"""Published MaBoSS figures, and the score of a Boolean engine against them.

ROADMAP §5 item 13 (area C): the MaBoSS adapter (genomeos/runtime/maboss.py) is scored on a published
MaBoSS model whose paper prints the probabilities, with the quantities, the tolerance and what each
engine is expected to do fixed here before any scored run, as the libRoadRunner lane did
(genomeos/molecules/sbml_reference.py).

THE MODEL. Cohen DPA, Martignetti L, Robine S, Barillot E, Zinovyev A, Calzone L (2015) Mathematical
Modelling of Molecular Pathways Enabling Tumour Cell Invasion and Migration. PLoS Comput Biol 11(11):
e1004571, doi:10.1371/journal.pcbi.1004571 (PMC4631357, CC BY 4.0). Its S1 File carries the master
model in MaBoSS format; data/models/cohen2015_metastasis.bnd and .cfg are those two files byte for
byte (SuppMat_Model_MaBoSS_MasterModel.bnd and .cfg; sha256 below). 32 nodes, every rate 1 ("All rates
are set to be 1 in this model since it is difficult to estimate them", Methods); four inputs
(ECMicroenv, DNAdamage, GF, TGFbeta) start random because their istate lines are commented out; six
output nodes; 50,000 trajectories, seed 100, max_time 50.

THE FIGURE. Fig 3 prints, in its four pie charts, the probability of each phenotype at the end of the
MaBoSS simulation, in whole percent (read from the published image, sha256 below):

    A. Wild type            HS 25, Metastasis 25, EMT 1, Apoptosis 49
    B. p53 LoF              HS 23, Metastasis 27, EMT 3, Apoptosis 47
    C. NICD GoF             Invasion 18, Metastasis 36, EMT 46      ("apoptosis is no longer reachable")
    D. p53 LoF NICD GoF     Metastasis 100                          ("only metastasis is observed")

THE PHENOTYPE OF A FINAL STATE (`label`). Table 2 lists the wild type's nine stable states with their
output nodes: HS has none on; the four apoptotic states have Apoptosis and CellCycleArrest; the two
EMT states EMT and CellCycleArrest; the two metastatic states Metastasis, Migration, Invasion, EMT and
CellCycleArrest. Fig 3C names an Invasion phenotype the table does not list. So a state is named by
the first rule that applies: no output on -> HS; Metastasis -> Metastasis; Apoptosis -> Apoptosis;
Migration -> Migration; Invasion -> Invasion; EMT -> EMT; CellCycleArrest alone -> CCA. Every Table 2
state gets its own column's name under this rule; the order among the rest is a choice, fixed here.

THE MUTANTS (`CONDITIONS`). The paper does not say how a LoF or GoF was imposed. Each is taken as the
node held at 0 (LoF) or 1 (GoF) from time 0: `MaBoSSModel.held` sets both of its rates to 0 and its
initial state to the value. This is what MaBoSS's own mutant convention (rate forced to the maximum
toward the value, the High/Low form pyMaBoSS's `mutate` writes) does once the first instant passes,
and what the stand-in does with a constant rule.

THE CHECKS. Per panel, each printed phenotype's probability, and one more: the total probability of
every phenotype the panel does not print. 4 + 1, 4 + 1, 3 + 1, 1 + 1 = 16 checks.

THE TOLERANCE (`TOLERANCE`), fixed now: a printed phenotype passes when the engine's probability is
within 0.02 of the printed percent / 100; the unprinted total passes when it is at most 0.02. Why 0.02:
the print rounds to whole percent (<= 0.005), and two independent 50,000-trajectory estimates of a
probability differ by sqrt(2) x sqrt(p(1-p)/50000) <= 0.0032 at one SE, so 3 SE is <= 0.0095; the sum,
0.0145, is rounded up. A model passes on an engine only when all 16 checks pass.

THE ENGINES, and what each is expected to do (`EXPECTED`), registered before any run:

  maboss          runtime/maboss.py on the published .bnd/.cfg as published (seed 100, 50,000
                  trajectories, max_time 50, the last time window [49.9, 50)). Expected: 16 of 16.
  standin_async   runtime/boolean.py's `step_async` (one node drawn uniformly, re-evaluated) under an
                  exponential clock of rate 32 (nodes) x 1 (the common rate), from the cfg's initial
                  distribution, to t = 50; 50,000 trajectories, seed 100; the network from
                  maboss.to_network. With every rate equal this is the same Markov chain as MaBoSS's,
                  so it is expected to pass 16 of 16 as well -- and to fail where the rates differ,
                  which to_network refuses rather than run.
  standin         runtime/boolean.py's committed analysis: the synchronous `attractor` from each of
                  the 16 initial states the four random inputs allow, weighted equally (a cycle's
                  states share its weight). No value is predicted. Synchronous update is not the
                  semantics the paper used, and it is expected not to pass all 16.

The two engines that run the same chain are also compared with each other (`AGREEMENT`): the total
variation distance between MaBoSS's and standin_async's final distributions over the output states,
per panel, is expected to be at most 0.02.
"""

from __future__ import annotations

ROOT_FILES = {
    "bnd": "data/models/cohen2015_metastasis.bnd",
    "cfg": "data/models/cohen2015_metastasis.cfg",
}

SOURCE = {
    "paper": (
        "Cohen DPA, Martignetti L, Robine S, Barillot E, Zinovyev A, Calzone L (2015) Mathematical Modelling "
        "of Molecular Pathways Enabling Tumour Cell Invasion and Migration. PLoS Comput Biol 11(11): e1004571"
    ),
    "doi": "10.1371/journal.pcbi.1004571",
    "pmcid": "PMC4631357",
    "licence": "CC BY 4.0 (PLoS Computational Biology)",
    "model_files": "https://journals.plos.org/ploscompbiol/article/file?type=supplementary&id=10.1371/journal.pcbi.1004571.s017",
    "model_zip_sha256": "b2cc7908f2e290c53b094defc519c20cf9d720157a48b63d9e835a8e15d54816",
    "model_members": {
        "SuppMat_models/SuppMat_Model_MaBoSS_MasterModel.bnd": "data/models/cohen2015_metastasis.bnd",
        "SuppMat_models/SuppMat_Model_MaBoSS_MasterModel.cfg": "data/models/cohen2015_metastasis.cfg",
    },
    "bnd_sha256": "02fc71367b480cef46e7bcc858c756e7fb7e510ccbbf1d7b5b3cc074b7780cc7",
    "cfg_sha256": "3acdbc7efb5d60942b16f40db238e0451945fa552e500b96adad796881b0bcb5",
    "figure": "Fig 3, https://journals.plos.org/ploscompbiol/article/figure/image?size=large&id=10.1371/journal.pcbi.1004571.g003",
    "figure_sha256": "91c8e1018222700f5590a4e90c4ccef4a1ecf240c12d13ccbf2280e9e8911687",
}

OUTPUTS = ("Metastasis", "Migration", "Invasion", "EMT", "Apoptosis", "CellCycleArrest")

#: the first rule that applies names the state; HS when no output is on
PRECEDENCE = ("Metastasis", "Apoptosis", "Migration", "Invasion", "EMT")

#: panel -> the nodes held, and their values
CONDITIONS: dict[str, dict[str, bool]] = {
    "A. Wild type": {},
    "B. p53 LoF": {"p53": False},
    "C. NICD GoF": {"NICD": True},
    "D. p53 LoF NICD GoF": {"p53": False, "NICD": True},
}

#: panel -> phenotype -> probability, as Fig 3 prints it (whole percent / 100)
FIGURE: dict[str, dict[str, float]] = {
    "A. Wild type": {"HS": 0.25, "Metastasis": 0.25, "EMT": 0.01, "Apoptosis": 0.49},
    "B. p53 LoF": {"HS": 0.23, "Metastasis": 0.27, "EMT": 0.03, "Apoptosis": 0.47},
    "C. NICD GoF": {"Invasion": 0.18, "Metastasis": 0.36, "EMT": 0.46},
    "D. p53 LoF NICD GoF": {"Metastasis": 1.00},
}

TOLERANCE = {"printed_phenotype_abs": 0.02, "unprinted_total_max": 0.02}

#: total variation distance between maboss and standin_async final output-state distributions, per panel
AGREEMENT = 0.02

ENGINES = ("maboss", "standin_async", "standin")

SETTINGS = {
    "maboss": "the published .cfg as published: seed_pseudorandom 100, sample_count 50000, max_time 50, "
    "time_tick 0.1, thread_count 4; the last window [49.9, 50)",
    "standin_async": {"trajectories": 50_000, "seed": 100, "t_end": 50.0},
    "standin": "synchronous attractor from each of the 2^4 initial states, equal weights; a cycle's "
    "states share its weight",
}

EXPECTED = {
    "maboss": "16 of 16",
    "standin_async": "16 of 16",
    "standin": "fewer than 16; no value predicted",
    "agreement": "maboss and standin_async within 0.02 total variation in all four panels",
}


def label(state: tuple[str, ...] | frozenset[str] | set[str]) -> str:
    """The phenotype of a final state, from the output nodes on in it (other nodes are ignored)."""
    on = set(state) & set(OUTPUTS)
    if not on:
        return "HS"
    for name in PRECEDENCE:
        if name in on:
            return name
    return "CCA"


def phenotypes(states: dict[tuple[str, ...], float]) -> dict[str, float]:
    """Sum state probabilities by phenotype."""
    out: dict[str, float] = {}
    for state, p in states.items():
        out[label(state)] = out.get(label(state), 0.0) + float(p)
    return out


def score(panel: str, probs: dict[str, float]) -> list[dict]:
    """The checks of one panel against the phenotype probabilities an engine produced."""
    printed = FIGURE[panel]
    tol = TOLERANCE["printed_phenotype_abs"]
    checks = []
    for name, ref in printed.items():
        got = float(probs.get(name, 0.0))
        checks.append(
            {
                "check": f"{name} probability",
                "figure": ref,
                "engine": round(got, 5),
                "difference": round(got - ref, 5),
                "bound": tol,
                "pass": abs(got - ref) <= tol,
            }
        )
    rest = {k: v for k, v in probs.items() if k not in printed and v > 0}
    total = float(sum(rest.values()))
    checks.append(
        {
            "check": "phenotypes the panel does not print, total",
            "figure": 0.0,
            "engine": round(total, 5),
            "which": {k: round(v, 5) for k, v in sorted(rest.items())},
            "bound": TOLERANCE["unprinted_total_max"],
            "pass": total <= TOLERANCE["unprinted_total_max"],
        }
    )
    return checks


def total_variation(p: dict, q: dict) -> float:
    keys = set(p) | set(q)
    return 0.5 * sum(abs(float(p.get(k, 0.0)) - float(q.get(k, 0.0))) for k in keys)
