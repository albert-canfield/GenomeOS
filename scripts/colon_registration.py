# SPDX-License-Identifier: AGPL-3.0-or-later
"""The registration for the colonic-crypt mechanism module, committed BEFORE its numbers are fitted.

Area E's third human mechanism module, beside haematopoiesis and the small-intestinal crypt
(scripts/intestine_registration.py, whose discipline this file follows; docs/AREA-E-NEXT.md holds the
written comparison of the colon with the small intestine). The list of gating factors below is
CLOSED: every entry names a factor, the step it gates and a published loss-of-function phenotype IN
THE COLON, with the fetched sentence it rests on quoted. Every citation was fetched in this lane
(PubMed E-utilities, PMC full text) and not recalled. A factor without a fetched colon phenotype is
in EXCLUDED with that as the reason, however well sourced it is in the small intestine, and
`check_factors` refuses a module that gates on a factor this file does not name.

Three kinds of agreement are registered in advance:

- "full": the published colon phenotype is a loss of lineages a flow module can state;
- "partial": the published colon phenotype is a reduction or a conversion, and a module whose gates
  are presence-or-absence can only delete, so the registered check is the deletion and the rest of
  the phenotype is named here as not expressible;
- "contested": two fetched papers disagree about the colon and the module's gate can agree with at
  most one. Such a knockout holding its registered pattern is NOT counted as reproducing the colon;
  it is counted as agreeing with the paper in `agrees_with` and disagreeing with the other.

READINGS are statements about the formalism registered before any run, one of them a number computed
from the registered shares alone, so that what the run shows cannot be explained afterwards. FITTED
names the one quantity fitted to a measurement; nothing about any mutant is fitted.

Run `python3 scripts/colon_registration.py` to print the registration.
"""

from __future__ import annotations

import json

ORGAN = "human colonic epithelium (crypt plus surface epithelium: no villus, no Paneth cells)"
MODULE = "data/organisms/human/colon.bio"
MUTANT_MODULE = "data/organisms/human/colon_mutants.bio"

# The module's compartments. The names are distinct from the small-intestinal module's so that no
# inventory across programs can read a colonic population as an intestinal one.
STEM = "ColonStem"
TA = "ColonTA"
SECRETORY = "ColonSecretoryProgenitor"
COLONOCYTE = "Colonocyte"
GOBLET = "ColonGobletCell"
EEC = "ColonEnteroendocrineCell"
MATURE = [COLONOCYTE, GOBLET, EEC]

SRC_MILO = "Sender & Milo 2021, Nat Med 27:45"
SRC_VANES2012 = "van Es et al. 2012, Mol Cell Biol 32:1918"
SRC_ANGUSHILL = "Angus-Hill et al. 2011, Proc Natl Acad Sci USA 108:4914"
SRC_PELLEGRINET = "Pellegrinet et al. 2011, Gastroenterology 140:1230"
SRC_LO = "Lo et al. 2016, Cell Mol Gastroenterol Hepatol 3:51"
SRC_MELLITZER = "Mellitzer et al. 2010, J Clin Invest 120:1708"
SRC_JENNY = "Jenny et al. 2002, EMBO J 21:6338"
SRC_SHROYER = "Shroyer et al. 2005, Genes Dev 19:2412"
SRC_BASTIDE = "Bastide et al. 2007, J Cell Biol 178:635"
SRC_KATZ = "Katz et al. 2002, Development 129:2619"
SRC_SASAKI = "Sasaki et al. 2016, Proc Natl Acad Sci USA 113:E5399"

# The closed list. `gates` are the module steps each factor is required at; the module may require a
# factor only at the steps named here. `quote` is the fetched sentence the gate rests on.
GATES = [
    {
        "id": "C1",
        "factor": "TCF7L2",
        "gates": [STEM],
        "match": "contested",
        "phenotype": "an adult inducible deletion abolishes the proliferating cells of the colon within "
        "7 days; a mosaic adult deletion instead enlarges the colonic crypts and leaves proliferative and "
        "differentiated cells present",
        "source": SRC_VANES2012,
        "quote": "the proliferating/Ki67+ cells in the colon derived from the Villin Creert2_Tcf4 LoxP/LoxP "
        "mice were abolished 7 days after cre induction ... Tcf4 ... played an essential role in the "
        "homeostasis of the adult colon",
        "contradicted_by": SRC_ANGUSHILL,
        "contradicting_quote": "loss of Tcf4 early in development and in adult colon results in increased "
        "cell proliferation ... adult colonic crypts lacking Tcf4 were enlarged compared with adjacent "
        "heterozygous crypts and ... there were proliferative and differentiated cell types present",
        "agrees_with": SRC_VANES2012,
    },
    {
        "id": "C2",
        "factor": "RBPJ",
        "gates": [TA],
        "match": "partial",
        "phenotype": "an inducible gut-specific deletion loses the proliferative crypt compartment "
        "completely and converts it into goblet cells, in the colon as in the proximal intestine",
        "source": SRC_PELLEGRINET,
        "quote": "Alcian blue (top panels) and Ki67 (bottom panels) staining on sections derived from the "
        "proximal intestine (A) and colon (B) of the indicated mice. Note the complete loss of the "
        "proliferative crypt compartment and the complete conversion into goblet cells in RBP-J and RBP-J "
        "Klf4 mutant intestine",
    },
    {
        "id": "C3",
        "factor": "ATOH1",
        "gates": [SECRETORY],
        "match": "full",
        "phenotype": "an adult inducible deletion leaves the secretory lineages nearly absent in the colon, "
        "the absorptive colonocytes remaining",
        "source": SRC_LO,
        "quote": "ATOH1 is required for all secretory lineages in ileum and colon ... After tamoxifen "
        "injection for 3 consecutive days, Alcian blue staining showed that the secretory lineages were "
        "nearly absent in both the ileum and colon compared with wild-type control Atoh1 lox/WT;VilCreERT2 "
        "mice",
    },
    {
        "id": "C4",
        "factor": "NEUROG3",
        "gates": [EEC],
        "match": "full",
        "phenotype": "an intestine-specific deletion leaves no enteroendocrine cells anywhere along the "
        "intestinal tract, the colon included, while goblet cells remain",
        "source": SRC_MELLITZER,
        "quote": "Sections of adult duodenum, jejunum, ileum, and colon were examined for the presence of "
        "endocrine cells ... This analysis revealed a complete lack of chromogranin A- and CCK/gastrin-"
        "positive enteroendocrine cells all along the intestinal tract",
        "also": SRC_JENNY,
        "also_quote": "Goblet cells were clearly identified in the small and large intestine by the "
        "production of mucus in mice lacking ngn3 ... The colon did not show a similar increase in goblet "
        "cell density",
    },
    {
        "id": "C5",
        "factor": "GFI1",
        "gates": [GOBLET],
        "match": "partial",
        "phenotype": "Gfi1-/- mice have about a third the normal number of goblet cells, a reduction also "
        "seen in the colon, and more enteroendocrine cells (PYY sixfold up in the colon)",
        "source": SRC_SHROYER,
        "quote": "Gfi1-/- mice had approximately one-third the normal number of goblet cells throughout the "
        "small intestine (p < 0.0001) ... This reduction in goblet cell number was also observed in the "
        "colon ... a significant increase (p < 0.0001) in all enteroendocrine sub-types examined across the "
        "small intestine and colon (... PYY, sixfold increase in colon)",
    },
    {
        "id": "C6",
        "factor": "SOX9",
        "gates": [GOBLET],
        "match": "partial",
        "phenotype": "an epithelium-wide deletion strongly decreases the colonic goblet population and "
        "leaves the enteroendocrine cells and the absorptive cells comparable to the control",
        "source": SRC_BASTIDE,
        "quote": "Alcian blue staining (a and b) and Muc2 immunoreactivity (c and d) reveal a strong "
        "reduction of the goblet cell population in the colon epithelium of Sox9 flox/flox-vil-Cre mice. "
        "Sox9 flox/flox and Sox9 flox/flox-vil-Cre mice have comparable populations of chromogranin "
        "A-positive enteroendocrine cells (e and f) and Cdx-2-positive enterocytes (g and h)",
    },
    {
        "id": "C7",
        "factor": "KLF4",
        "gates": [GOBLET],
        "match": "contested",
        "phenotype": "Klf4-/- newborns have 90% fewer colonic goblet cells, the colonocytes and "
        "enteroendocrine cells being normal; but an adult inducible deletion changes no goblet cell number "
        "in the colon at all",
        "source": SRC_KATZ,
        "quote": "Klf4-/- mice showed a 90% decrease in the ratio of mature goblet cells to total epithelial "
        "cells compared with wild-type littermates ... All other epithelial cell types are present in the "
        "colon of Klf4-/- mice",
        "contradicted_by": SRC_PELLEGRINET,
        "contradicting_quote": "histological analysis with Alcian blue staining did not reveal any "
        "differences in goblet cell numbers between control and Klf4 mutant animals in either the small "
        "intestine or the colon",
        "agrees_with": SRC_KATZ,
    },
]

# Excluded, each with the fetched reason. The list is closed and none of these enters the module.
EXCLUDED = [
    {
        "factor": "ASCL2",
        "source": "van der Flier et al. 2009, Cell 136:903",
        "phenotype": "an induced deletion in the adult small intestine makes the Lgr5 stem cells disappear "
        "within days",
        "reason": "the fetched record states the small intestine and says nothing about the colon, so "
        "carrying the small-intestinal gate across would be an assumption and not a source. Excluded for "
        "want of a colon phenotype, not for want of genetics",
    },
    {
        "factor": "SPDEF",
        "source": "Gregorieff et al. 2009, Gastroenterology 137:1333",
        "phenotype": "loss severely impairs the maturation of goblet and Paneth cells and accumulates "
        "immature secretory progenitors",
        "reason": "the gate of the small-intestinal module, and its abstract says 'intestinal epithelium' "
        "without naming the colon; the publisher's full text could not be fetched in this lane (HTTP 403), "
        "so no colon sentence exists to quote. Lo et al. 2016 use Spdef-null colonic crypts for gene "
        "expression, which is not a lineage loss. Excluded for want of a fetched colon phenotype",
    },
    {
        "factor": "REG4",
        "source": SRC_SASAKI,
        "phenotype": "diphtheria-toxin ablation of the Reg4+ deep crypt secretory cells makes the colonic "
        "Lgr5+ stem cells disappear and disrupts homeostasis",
        "reason": "the colon's Paneth-cell equivalent, and the one population this module would need to add "
        "to the small-intestinal pattern; but the experiment ablates the cells rather than removing a "
        "factor, and no share of daily output or residence time for deep crypt secretory cells was found, "
        "so the lineage would be a population with invented numbers. Its measured 12-14 Reg4+ cells per "
        "crypt in the caecum and distal colon is recorded in MEASURED and used for nothing",
    },
    {
        "factor": "POU2F3",
        "source": "Gerbe et al. 2016, Nature 529:226",
        "phenotype": "Pou2f3-/- mice lack intestinal tuft cells",
        "reason": "the same exclusion the small-intestinal registration made and for the same reason: the "
        "factor is sourced, the quantities are not",
    },
    {
        "factor": "HES1",
        "source": "Jensen et al. 2000, Nat Genet 24:36",
        "phenotype": "precocious and excessive differentiation of endocrine cell types in the developing "
        "stomach and gut",
        "reason": "the phenotype is an increase, and a flow module cannot re-specify a fate, so gating the "
        "absorptive branch on HES1 would predict a loss no paper reports (carried from the small-intestinal "
        "registration, which excluded it for this reason)",
    },
    {
        "factor": "NOX1",
        "source": "Coant et al. 2010, Mol Cell Biol 30:2636",
        "phenotype": "NOX1-deficient mice show a massive conversion of progenitor cells into post-mitotic "
        "goblet cells at the cost of colonocytes",
        "reason": "colon-specific and well sourced, but the phenotype is a conversion towards goblet cells: "
        "a module made of flows would have to gate the colonocyte branch on it and would then predict a "
        "loss of colonocytes with no compensating goblet gain, which is not what the paper reports",
    },
    {
        "factor": "DLL1 with DLL4",
        "source": SRC_PELLEGRINET,
        "phenotype": "single Dll1 loss raises goblet numbers 1.83-fold; Dll1 with Dll4 converts the "
        "proliferative compartment into goblet cells, in the colon too",
        "reason": "the single-ligand phenotype is an increase, which the formalism cannot express, and the "
        "double-ligand phenotype is the RBPJ phenotype already gated at C2, so adding the ligands would "
        "duplicate a gate rather than add a mechanism",
    },
]

# What each knockout must reproduce, written before the module existed. `lost` must reach zero,
# `kept` must stay above KEPT_FLOOR, and `match` is the agreement registered in advance.
MUTANTS = [
    {
        "name": "tcf7l2",
        "knockout": "TCF7L2",
        "lost": MATURE,
        "kept": [],
        "match": "contested",
        "agrees_with": SRC_VANES2012,
        "disagrees_with": SRC_ANGUSHILL,
        "note": "the stem compartment does not exist, so no lineage is produced. Registered contested in "
        "advance: this agrees with the abolished colonic proliferation of van Es et al. 2012 and disagrees "
        "with the enlarged, still-differentiating crypts of Angus-Hill et al. 2011",
    },
    {
        "name": "rbpj",
        "knockout": "RBPJ",
        "lost": MATURE,
        "kept": [],
        "match": "partial",
        "note": "registered partial in advance: the complete loss of the proliferative compartment is the "
        "match, and the published conversion of those cells into colonic goblet cells is not expressible in "
        "a module made of flows, which deletes the goblet lineage with the rest",
    },
    {
        "name": "atoh1",
        "knockout": "ATOH1",
        "lost": [GOBLET, EEC],
        "kept": [COLONOCYTE],
        "match": "full",
        "note": "the headline prediction: the secretory branch disappears from the colon and the absorptive "
        "branch does not",
    },
    {
        "name": "neurog3",
        "knockout": "NEUROG3",
        "lost": [EEC],
        "kept": [COLONOCYTE, GOBLET],
        "match": "full",
        "note": "one secretory type only, and the colonic goblet cells remain, as Jenny et al. 2002 report",
    },
    {
        "name": "gfi1",
        "knockout": "GFI1",
        "lost": [GOBLET],
        "kept": [COLONOCYTE, EEC],
        "match": "partial",
        "note": "registered partial in advance: the published colonic phenotype is a reduction in goblet "
        "cells, and a presence-or-absence gate can only delete them; the published rise in colonic "
        "enteroendocrine cells (PYY sixfold) has the direction READINGS R2 predicts and no fold this lane "
        "can compare with",
    },
    {
        "name": "sox9",
        "knockout": "SOX9",
        "lost": [GOBLET],
        "kept": [COLONOCYTE, EEC],
        "match": "partial",
        "note": "registered partial in advance: 'a strong reduction of the goblet cell population' is stated "
        "as a deletion, while the enteroendocrine and absorptive populations Bastide et al. found comparable "
        "are the registered kept lineages",
    },
    {
        "name": "klf4",
        "knockout": "KLF4",
        "lost": [GOBLET],
        "kept": [COLONOCYTE, EEC],
        "match": "contested",
        "agrees_with": SRC_KATZ,
        "disagrees_with": SRC_PELLEGRINET,
        "note": "the one colon-specific gate, excluded from the small-intestinal module for being "
        "colon-specific. Registered contested in advance: the 90% decrease of Katz et al. 2002 in the "
        "newborn colon is the agreement, and the unchanged adult goblet numbers of Pellegrinet et al. 2011 "
        "are the disagreement; the module has no age axis and therefore cannot hold both",
    },
]
KEPT_FLOOR = 1.0e6  # a population the knockout does not touch must stay above this

# Sourced quantities, used as given.
MEASURED = [
    {
        "what": "colon epithelial cells in the reference adult",
        "value": 6.0e10,
        "source": f"{SRC_MILO} (distilled: data/results/human_cell_turnover.json, row 'Colon epithelia')",
    },
    {
        "what": "colon epithelial cells replaced per day",
        "value": 1.7e10,
        "source": SRC_MILO,
    },
    {
        "what": "colon epithelial lifespan, published mean residence (days)",
        "value": 3.4,
        "source": SRC_MILO,
    },
    {
        "what": "mature goblet cells as a share of colonic epithelial cells, wild-type mouse colon at P1",
        "value": 0.20,
        "source": f"{SRC_KATZ}: 'In wild-type colons, 20% of epithelial cells were mature goblet cells; "
        "this fraction was reduced to 2% in Klf4-/- mice'",
    },
    {
        "what": "Lgr5+ stem cells per colonic crypt (mouse)",
        "value": 20.0,
        "source": f"{SRC_SASAKI}: 'approximately 20 Lgr5+ stem cells per crypt were detected in each region'",
    },
    {
        "what": "Reg4+ deep crypt secretory cells per crypt, caecum and distal colon (mouse)",
        "value": 13.0,
        "source": f"{SRC_SASAKI}: '12-14 Reg4+ cells in cecum, distal colon'",
        "used_for": "nothing: the lineage is excluded, and this number is recorded so the exclusion can be "
        "read as a choice and not an oversight",
    },
    {
        "what": "enteroendocrine share of the epithelium (upper bound)",
        "value": 0.01,
        "source": "the bound registered by the small-intestinal module (scripts/intestine_registration.py, "
        "460c46f), itself citing Cheng & Leblond 1974, Am J Anat 141:461 for <1% of the epithelium of "
        "stomach, small intestine and colon; not re-fetched in this lane",
    },
]

# Fitted, and to what. Nothing else in the module is fitted, and no mutant number is fitted at all.
FITTED = [
    {
        "what": "surface_residence_days, shared by colonocytes, goblet cells and enteroendocrine cells",
        "fitted_to": "the published colonic epithelial count of 6.0e10 cells: the residence that makes the "
        "sum of every compartment's pool equal it",
        "against": "the published mean residence of 3.4 days and the 3.53 days the two measured numbers "
        "imply on their own",
    },
]
INFERRED = [
    "the goblet share of daily output is taken to equal its measured standing share of 0.20 (Katz et al. "
    "2002). That holds only if goblet cells and colonocytes have the same residence, which the colon makes "
    "more defensible than the small intestine did: there are no Paneth cells and no long-lived resident "
    "lineage, every mature colonic type rides the same crypt-to-surface conveyor and is shed at the "
    "surface. It is an inference and not a fit, and FITTED's one knob is what closes the count",
    "the enteroendocrine share of daily output, 0.004, inside the <1% bound above and fitted to nothing",
    "the colonocyte share of output, 0.796, as the remainder",
    "the stem-cell division interval, 1 day, with no primary colon citation fetched in this lane",
    "the transit-amplifying transit time (2.5 d) and division count (5), with no primary citation fetched; "
    "the same inference the small-intestinal module made",
    "the secretory progenitor residence (1 d), with no primary citation fetched",
]
NOT_CLAIMED = [
    "the colonic stem-cell pool and the crypt count it implies at the measured 20 Lgr5+ cells per crypt: "
    "both are consequences of the flows, reported and checked against nothing",
    "the 20% goblet standing share: it comes back by construction, because the share was inferred as the "
    "output share and all three mature residences are the one fitted number. It is not a check and must "
    "never be reported as one",
]

# Statements about the formalism, registered before the module ran. R2 carries a number computed from
# the registered shares alone, so it can be wrong in writing.
READINGS = [
    {
        "id": "R1",
        "says": "two of the seven gates are contested in the literature, so no run can settle them: a "
        "knockout that holds its registered pattern agrees with the paper in `agrees_with` and disagrees "
        "with the other, and the module's count of reproduced phenotypes must say so rather than count "
        "seven of seven",
    },
    {
        "id": "R2",
        "says": "blocking the goblet route (gfi1, sox9 or klf4) will leave the enteroendocrine pool at "
        "about 51 times its wild-type value, by mass conservation and not by fate re-specification: the "
        "secretory progenitor backs up until its one open outflow carries the whole inflow, so the factor "
        "is 1 / (0.004 / 0.204) = 51.0. The direction is the published one for Gfi1 (PYY sixfold up in the "
        "colon) and the magnitude is a property of the registered shares, claimed of nothing",
        "number": 51.0,
    },
    {
        "id": "R3",
        "says": "without ATOH1 the colonocytes will stand ABOVE their wild-type pool, because each "
        "compartment's growth term was solved for the wild-type outflows and the secretory outflow is gone. "
        "Only the direction is registered: the transit-amplifying compartment carries a growth term and a "
        "cap, so the magnitude is not predictable from the shares. Lo et al. 2016 report the absorptive "
        "lineage remaining, not expanding",
    },
    {
        "id": "R4",
        "says": "the published RBPJ phenotype is a conversion into goblet cells and the module will delete "
        "every lineage instead, the goblet cells included. That is the formalism and not the colon",
    },
    {
        "id": "R5",
        "says": "three of the seven published colonic phenotypes are reductions (a third of the goblet "
        "cells for Gfi1, 'a strong reduction' for Sox9, 90% for Klf4) and a presence-or-absence gate can "
        "only delete, so those three overstate their phenotype by construction. It is worse here than in "
        "the small intestine, where the deleted lineage was 4% of the epithelium: in the colon the goblet "
        "lineage is 20% of it",
    },
    {
        "id": "R6",
        "says": "the colon needs no Paneth compartment, and its Paneth equivalent (Reg4+ deep crypt "
        "secretory cells) is excluded for want of a quantity, so this module has three mature lineages "
        "where the small-intestinal one has four. With the villus gone, that is the whole structural "
        "difference between the two programs",
    },
]

# Bounds the wild-type run must clear. The first two are fit closure, not predictions.
BOUNDS = {
    "total_epithelium_at_1_yr": {"low": 5.1e10, "high": 6.9e10, "fit_closure": True},
    "turnover_per_day": {"low": 1.445e10, "high": 1.955e10, "fit_closure": True},
    "every_mature_population_present": {"fit_closure": False},
}

REGISTRATION = {
    "registered": "2026-10-10",
    "lane": "lane-colon",
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
    "readings": READINGS,
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
