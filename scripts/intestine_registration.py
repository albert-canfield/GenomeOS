# SPDX-License-Identifier: AGPL-3.0-or-later
"""The registration for the small-intestinal mechanism module, committed BEFORE its numbers are fitted.

Area E's second human mechanism module after haematopoiesis (docs/AREA-E-NEXT.md holds the written
choice between this and the worm's next number). The list of gating factors below is CLOSED: every
entry names a factor, the step it gates and a published loss-of-function phenotype, each citation
fetched and read in this lane rather than recalled. Anything without a source is in `EXCLUDED` with
the reason, and `check_factors` refuses a module that gates on a factor this file does not name, so
the list cannot grow between the registration and the run without a commit saying so.

`MUTANTS` is what the knockouts must reproduce, written before the module existed. Two entries are
registered as PARTIAL in advance, because their published phenotype is a re-allocation (more
enteroendocrine cells without GFI1; proliferative cells converted to goblet cells without RBPJ) and
a population module made of flows can delete a lineage but cannot re-specify one. `FITTED` names the
one quantity fitted to a measurement, so that no part of the result can be read as a prediction it
is not: the knockout pattern is the prediction, and nothing about any mutant was fitted.

Run `python3 scripts/intestine_registration.py` to print the registration.
"""

from __future__ import annotations

import json

ORGAN = "human small-intestinal epithelium (crypt-villus)"
MODULE = "data/organisms/human/intestine.bio"
MUTANT_MODULE = "data/organisms/human/intestine_mutants.bio"

# The closed list. `gates` are the module steps each factor is required at; the module may require a
# factor only at the steps named here.
GATES = [
    {
        "id": "S1",
        "factor": "TCF7L2",
        "gates": ["ISC"],
        "phenotype": "no proliferative compartment is maintained in the prospective crypt regions; the "
        "neonatal epithelium is composed entirely of differentiated, non-dividing villus cells",
        "source": "Korinek et al. 1998, Nat Genet 19:379",
    },
    {
        "id": "S2",
        "factor": "ASCL2",
        "gates": ["ISC"],
        "phenotype": "induced deletion in adult small intestine makes the Lgr5 stem cells disappear "
        "within days",
        "source": "van der Flier et al. 2009, Cell 136:903",
    },
    {
        "id": "S3",
        "factor": "RBPJ",
        "gates": ["TA"],
        "phenotype": "conditional removal converts proliferative crypt cells into post-mitotic goblet cells",
        "source": "van Es et al. 2005, Nature 435:959",
    },
    {
        "id": "S4",
        "factor": "ATOH1",
        "gates": ["SecretoryProgenitor"],
        "phenotype": "goblet, enteroendocrine and Paneth cells are depleted; enterocytes are unaffected",
        "source": "Yang et al. 2001, Science 294:2155",
    },
    {
        "id": "S5",
        "factor": "NEUROG3",
        "gates": ["EnteroendocrineCell"],
        "phenotype": "null mice generate no intestinal endocrine cells and endocrine progenitors are "
        "lacking, while the other main intestinal epithelial cell types differentiate properly",
        "source": "Jenny et al. 2002, EMBO J 21:6338",
    },
    {
        "id": "S6",
        "factor": "GFI1",
        "gates": ["GobletCell", "PanethCell"],
        "phenotype": "Gfi1-/- mice lack Paneth cells, have fewer goblet cells and supernumerary "
        "enteroendocrine cells",
        "source": "Shroyer et al. 2005, Genes Dev 19:2412",
    },
    {
        "id": "S7",
        "factor": "SPDEF",
        "gates": ["GobletCell", "PanethCell"],
        "phenotype": "loss severely impairs the maturation of goblet and Paneth cells and accumulates "
        "immature secretory progenitors",
        "source": "Gregorieff et al. 2009, Gastroenterology 137:1333",
    },
    {
        "id": "S8",
        "factor": "SOX9",
        "gates": ["PanethCell"],
        "phenotype": "in the absence of SOX9 Paneth cells are not formed and the other epithelial cell "
        "types are unaffected",
        "source": "Mori-Akiyama et al. 2007, Gastroenterology 133:539; "
        "Bastide et al. 2007, J Cell Biol 178:635",
    },
]

EXCLUDED = [
    {
        "factor": "POU2F3",
        "source": "Gerbe et al. 2016, Nature 529:226",
        "phenotype": "Pou2f3-/- mice lack intestinal tuft cells",
        "reason": "the factor is sourced but the quantities are not: no measurement of the tuft share of "
        "daily output or of a tuft residence time was found, so the lineage would be a population with "
        "invented numbers",
    },
    {
        "factor": "KLF4",
        "source": "Katz et al. 2002, Development 129:2619",
        "phenotype": "Klf4-/- mice show a 90% decrease in the number of goblet cells in the colon",
        "reason": "colon-specific; this module is the small intestine",
    },
    {
        "factor": "HES1",
        "source": "Jensen et al. 2000, Nat Genet 24:36",
        "phenotype": "precocious and excessive differentiation of endocrine cell types in the developing "
        "stomach and gut",
        "reason": "the phenotype is an increase, and a flow module cannot re-specify a fate, so gating the "
        "absorptive branch on HES1 would predict a loss no paper reports",
    },
    {
        "factor": "CDX2",
        "source": None,
        "phenotype": None,
        "reason": "considered; no citation fetched in this lane, and its published phenotype is a change of "
        "organ identity rather than the loss of a lineage inside the epithelium",
    },
    {
        "factor": "LGR5",
        "source": None,
        "phenotype": None,
        "reason": "considered; no citation fetched in this lane: it marks the stem cell without a "
        "loss-of-function phenotype that removes the compartment",
    },
]

MATURE = ["Enterocyte", "GobletCell", "PanethCell", "EnteroendocrineCell"]

# What each knockout must reproduce. `lost` must reach zero, `kept` must stay above `KEPT_FLOOR`.
MUTANTS = [
    {
        "name": "tcf7l2",
        "knockout": "TCF7L2",
        "lost": MATURE,
        "kept": [],
        "match": "full",
        "note": "the stem compartment does not exist, so no lineage is produced",
    },
    {
        "name": "ascl2",
        "knockout": "ASCL2",
        "lost": MATURE,
        "kept": [],
        "match": "full",
        "note": "the stem compartment does not exist, so no lineage is produced",
    },
    {
        "name": "rbpj",
        "knockout": "RBPJ",
        "lost": MATURE,
        "kept": [],
        "match": "partial",
        "note": "registered partial in advance: the loss of the proliferative compartment is the match, the "
        "published conversion of those cells into goblet cells is not expressible in a flow module",
    },
    {
        "name": "atoh1",
        "knockout": "ATOH1",
        "lost": ["GobletCell", "PanethCell", "EnteroendocrineCell"],
        "kept": ["Enterocyte"],
        "match": "full",
        "note": "the headline prediction: the secretory branch disappears and the absorptive branch does not",
    },
    {
        "name": "neurog3",
        "knockout": "NEUROG3",
        "lost": ["EnteroendocrineCell"],
        "kept": ["Enterocyte", "GobletCell", "PanethCell"],
        "match": "full",
        "note": "one secretory type only",
    },
    {
        "name": "gfi1",
        "knockout": "GFI1",
        "lost": ["GobletCell", "PanethCell"],
        "kept": ["Enterocyte", "EnteroendocrineCell"],
        "match": "partial",
        "note": "registered partial in advance: 'lack Paneth cells' is a match, 'fewer goblet cells' is "
        "overstated by a module that can only delete, and the supernumerary enteroendocrine cells cannot "
        "be expressed",
    },
    {
        "name": "spdef",
        "knockout": "SPDEF",
        "lost": ["GobletCell", "PanethCell"],
        "kept": ["Enterocyte", "EnteroendocrineCell"],
        "match": "full",
        "extra": "SecretoryProgenitor above its wild-type pool (the published accumulation of immature "
        "secretory progenitors)",
        "note": "the only mutant with a registered expectation about a progenitor pool",
    },
    {
        "name": "sox9",
        "knockout": "SOX9",
        "lost": ["PanethCell"],
        "kept": ["Enterocyte", "GobletCell", "EnteroendocrineCell"],
        "match": "full",
        "note": "one secretory type only",
    },
]
KEPT_FLOOR = 1.0e6  # a population the knockout does not touch must stay above this

# Sourced quantities, used as given.
MEASURED = [
    {
        "what": "small-intestine epithelial cells in the reference adult",
        "value": 1.2e11,
        "source": "Sender & Milo 2021, Nat Med 27:45 (distilled: data/results/human_cell_turnover.json)",
    },
    {
        "what": "small-intestine epithelial cells replaced per day",
        "value": 2.2e10,
        "source": "Sender & Milo 2021, Nat Med 27:45",
    },
    {
        "what": "goblet cell turnover time, mouse duodenum (days)",
        "value": 3.0,
        "source": "Cheng 1974, Am J Anat 141:481 (II. Mucous cells)",
    },
    {
        "what": "Paneth cell turnover time, mouse (days)",
        "value": 15.0,
        "source": "Cheng 1974, Am J Anat 141:521 (IV. Paneth cells)",
    },
    {
        "what": "stem-cell division interval, mouse Lgr5 crypt base columnar cell (days)",
        "value": 1.0,
        "source": "Snippert et al. 2010, Cell 143:134",
    },
    {
        "what": "enterocyte share of villus epithelial cells (bound, mouse)",
        "value": 0.95,
        "source": "Cheng & Leblond 1974, Am J Anat 141:461 (I. Columnar cells): >95% of villus cells",
    },
    {
        "what": "enteroendocrine share of epithelial cells (bound)",
        "value": 0.01,
        "source": "Cheng & Leblond 1974, Am J Anat 141:461: enteroendocrine cells are <1% of the "
        "epithelium of stomach, small intestine and colon",
    },
]

# Fitted, and to what. Nothing else in the module is fitted, and no mutant number is fitted at all.
FITTED = [
    {
        "what": "villus_residence_days (shared by enterocytes and enteroendocrine cells)",
        "fitted_to": "the published human epithelial count of 1.2e11 cells: the residence that makes the "
        "sum of every compartment's pool equal it",
        "against": "the mouse-measured absorptive turnover of about 3 days",
    },
]
INFERRED = [
    "the lineage shares of daily output (enterocyte 0.95, goblet 0.04, enteroendocrine 0.008, Paneth the "
    "remaining 0.002), taken from the sourced bounds above and from the ~4% duodenal goblet share of the "
    "Karam 1999 review (Front Biosci 4:D286), not fitted to anything",
    "the transit-amplifying transit time (2.5 d) and division count (5), with no primary citation fetched "
    "in this lane",
    "the secretory progenitor residence (1 d), with no primary citation fetched in this lane",
]
NOT_CLAIMED = [
    "the human stem-cell pool and the crypt count it implies: both are consequences of the flows, "
    "reported and checked against nothing",
]

# Bounds the wild-type run must clear. The first two are fit closure, not predictions.
BOUNDS = {
    "total_epithelium_at_1_yr": {"low": 1.02e11, "high": 1.38e11, "fit_closure": True},
    "turnover_per_day": {"low": 1.87e10, "high": 2.53e10, "fit_closure": True},
    "every_mature_population_present": {"fit_closure": False},
}

REGISTRATION = {
    "registered": "2026-10-10",
    "lane": "lane-organism",
    "organ": ORGAN,
    "module": MODULE,
    "mutant_module": MUTANT_MODULE,
    "gates": GATES,
    "excluded": EXCLUDED,
    "mature": MATURE,
    "mutants": MUTANTS,
    "kept_floor": KEPT_FLOOR,
    "measured": MEASURED,
    "fitted": FITTED,
    "inferred": INFERRED,
    "not_claimed": NOT_CLAIMED,
    "bounds": BOUNDS,
}


def factors() -> set[str]:
    """The closed set of factors the module may gate on."""
    return {g["factor"] for g in GATES}


def check_factors(used: set[str] | list[str]) -> None:
    """Refuse a module that gates on a factor this registration does not name."""
    extra = sorted(set(used) - factors())
    if extra:
        raise ValueError(f"not in the registered closed list: {', '.join(extra)}")


def check_gate(factor: str, step: str) -> None:
    """Refuse a gate at a step the registration did not name for that factor."""
    for g in GATES:
        if g["factor"] == factor:
            if step not in g["gates"]:
                raise ValueError(f"{factor} is registered at {g['gates']}, not at {step}")
            return
    raise ValueError(f"not in the registered closed list: {factor}")


def main() -> None:
    print(json.dumps(REGISTRATION, indent=2))


if __name__ == "__main__":
    main()
