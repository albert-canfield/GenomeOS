# SPDX-License-Identifier: AGPL-3.0-or-later
"""The registration for the gastric-gland mechanism module, committed BEFORE its numbers are fitted.

Area E's fourth human mechanism module, after haematopoiesis, the small-intestinal crypt
(scripts/intestine_registration.py) and the colonic crypt (scripts/colon_registration.py), whose
discipline this file follows; docs/AREA-E-NEXT.md holds the written comparison of the stomach with
both gut modules. The list of gating factors below is CLOSED: every entry names a factor, the step it
gates and a published loss-of-function phenotype IN THE STOMACH, with the fetched sentence it rests on
quoted. Every citation was fetched in this lane (PubMed E-utilities, PMC full text, EBI OLS4 for the
Cell Ontology ids) and not recalled. A factor without a fetched stomach phenotype is in EXCLUDED with
that as the reason, however well sourced it is in the intestine or the colon - ATOH1, the secretory
gate of both gut modules, is excluded here for exactly that reason - and `check_factors` refuses a
module that gates on a factor this file does not name.

Three kinds of agreement are registered in advance, as in the colonic registration:

- "full": the published stomach phenotype is a loss of a lineage a flow module can state;
- "partial": the published stomach phenotype is a reduction, a conversion, a maturation defect or an
  increase, and a module whose gates are presence-or-absence can only delete, so the registered check
  is the deletion and the rest of the phenotype is named here as not expressible;
- "contested": two fetched stomach EXPERIMENTS disagree and the module's single gate can satisfy at
  most one of them. The colonic registration required two papers; this one WIDENS the category to two
  experiments, because the stomach's Notch result is a paper whose own genetic and pharmacological
  arms differ in severity, and a widening that makes a gate count for less is the safe direction. A
  contested knockout holding its registered pattern is NOT counted as reproducing the stomach; it is
  counted as agreeing with the experiment in `agrees_with` and disagreeing with the other.

The stomach's own finding, registered here before any run: unlike the intestine and the colon, most
well-sourced gastric transcription-factor phenotypes are MATURATION defects with the cell still
present (XBP1, HNF4A, SPDEF, FOXQ1, KLF4 below), not lineage losses. That is why this closed list has
five gates and eight exclusions, and why the lineage that carries most of the stomach's turnover - the
pit (foveolar) cell - has NO gate in this module at all.

READINGS are statements about the formalism registered before any run, three of them numbers computed
from the registered inputs alone, so that what the run shows cannot be explained afterwards. FITTED
names the one quantity fitted to a measurement; nothing about any mutant is fitted.

Run `python3 scripts/stomach_registration.py` to print the registration.
"""

from __future__ import annotations

import json
import math

ORGAN = (
    "human gastric epithelium (the oxyntic/corpus gland unit: pit, isthmus, neck and base, with the "
    "antral gland's endocrine output folded into the one endocrine lineage)"
)
MODULE = "data/organisms/human/stomach.bio"
MUTANT_MODULE = "data/organisms/human/stomach_mutants.bio"

# The module's compartments. The names are distinct from the intestinal and colonic modules' so that
# no inventory across programs can read a gastric population as a gut one.
STEM = "GastricIsthmusStem"
NECK = "GastricNeckCell"
PIT = "GastricPitCell"
PARIETAL = "GastricParietalCell"
CHIEF = "GastricChiefCell"
EEC = "GastricEnteroendocrineCell"
MATURE = [PIT, PARIETAL, CHIEF, EEC]

SRC_MILO = "Sender & Milo 2021, Nat Med 27:45"
SRC_KARAM1 = "Karam & Leblond 1993, Anat Rec 236:259 (part I: proliferative cell types and the stem cell)"
SRC_KARAM2 = "Karam & Leblond 1993, Anat Rec 236:280 (part II: outward migration of pit cells)"
SRC_KARAM3 = "Karam & Leblond 1993, Anat Rec 236:297 (part III: neck cells into zymogenic cells)"
SRC_KARAM4 = "Karam 1993, Anat Rec 236:314 (part IV: bidirectional migration of parietal cells)"
SRC_KARAM5 = "Karam & Leblond 1993, Anat Rec 236:333 (part V: entero-endocrine and caveolated cells)"
SRC_KIM = "Kim & Shivdasani 2011, J Exp Med 208:677"
SRC_ADKINS = "Adkins-Threats et al. 2024, Cell Stem Cell 31:886"
SRC_RAMSEY = "Ramsey et al. 2007, Development 134:211"
SRC_LEE = "Lee et al. 2002, Genes Dev 16:1488"
SRC_HIBDON = "Hibdon et al. 2023, Am J Physiol Gastrointest Liver Physiol 325:G458"
SRC_HUH = "Huh et al. 2010, Gastroenterology 139:2038"
SRC_MOORE = "Moore et al. 2016, Am J Physiol Gastrointest Liver Physiol 311:G267"
SRC_HORST = "Horst et al. 2010, J Biol Chem 285:35047"
SRC_KATZ = "Katz et al. 2005, Gastroenterology 128:935"
SRC_VERZI = "Verzi et al. 2008, Gastroenterology 135:591"
SRC_QUE = "Que et al. 2007, Development 134:2521"
SRC_DELAFOREST = "DeLaForest et al. 2021, Cell Mol Gastroenterol Hepatol 12:1391"
SRC_BARKER = "Barker et al. 2010, Cell Stem Cell 6:25"

# The closed list. `gates` are the module steps each factor is required at; the module may require a
# factor only at the steps named here. `quote` is the fetched sentence the gate rests on.
GATES = [
    {
        "id": "S1",
        "factor": "RBPJ",
        "gates": [STEM],
        "match": "contested",
        "phenotype": "pan-Notch inhibition virtually abolishes proliferation in corpus and antral glands, "
        "leaving no Ki67+ cells in most glands; but the genetic deletion of the Notch effector RBP-Jk in "
        "the stomach epithelium leaves its gland units alive and merely at a competitive disadvantage, 47% "
        "fewer than in controls",
        "source": SRC_KIM,
        "quote": "Although stomach epithelial proliferation was only mildly reduced at the end of the 5-d "
        "treatment, proliferation in corpus and antral glands was virtually abolished 2 d later, leaving "
        "few Ki67+ cells in some isthmi and no Ki67+ cells in most glands",
        "contradicted_by": f"{SRC_KIM}, its own genetic arm (Osr1-Cre; RBP-Jk Fl/Fl)",
        "contradicting_quote": "Similarly, YFP+ gastric gland units were reduced by 47% compared with "
        "control animals ... Thus, RBP-Jk-null intestinal and stomach epithelia both face a significant "
        "competitive disadvantage, allowing rapid dominance by clonal units that are proficient in Notch "
        "signaling",
        "agrees_with": f"{SRC_KIM}, its pharmacological arm (the pan-Notch inhibitor DBZ)",
    },
    {
        "id": "S2",
        "factor": "ESRRG",
        "gates": [PARIETAL],
        "match": "full",
        "phenotype": "deleting Esrrg in cells already committed to the parietal lineage leaves whole "
        "swaths of gastric units with no parietal cells at all, from birth onward",
        "source": SRC_ADKINS,
        "quote": "Atp4bCre; Esrrg flox/flox mice had whole swaths of gastric units completely lacking PCs "
        "... Unlike other characterized PC-specific genes, deletion of Esrrg leads to the absence of PCs "
        "from birth onward",
        "also": SRC_ADKINS,
        "also_quote": "When we crossed Fgf20Cre with Esrrg flox/flox mice, we saw dramatic loss of PCs in "
        "glands lacking ERRgamma in adult mice",
        "not_expressible": "the same mice show 'a significant, steady-state increase in proliferation', "
        "which a module made of flows cannot state; only the parietal loss is checked",
    },
    {
        "id": "S3",
        "factor": "MIST1",
        "gates": [CHIEF],
        "match": "partial",
        "phenotype": "Mist1-/- neck cells differentiate normally, but the mature chief (zymogenic) cells "
        "of the gland base are uniformly structurally abnormal and the immature transitional cells "
        "upstream of them accumulate markedly; the gene is Bhlha15",
        "source": SRC_RAMSEY,
        "quote": "By contrast, Mist1(-/-) neck cells differentiate normally, but ZCs in the mature, basal "
        "portion of the gastric unit uniformly exhibit multiple apical cytoplasmic structural "
        "abnormalities. This defect in terminal ZC differentiation is also associated with markedly "
        "increased abundance of TCs, especially in late-stage TCs that predominantly have features of "
        "immature ZCs",
        "not_expressible": "the mature chief cell is PRESENT and abnormal, not absent, and this module "
        "deletes it. The accumulation of immature transitional cells is the part the formalism does state, "
        "because the neck compartment has no second outflow once the chief route is gated off (READINGS "
        "R4)",
    },
    {
        "id": "S4",
        "factor": "NEUROG3",
        "gates": [EEC],
        "match": "partial",
        "phenotype": "in ngn3-/- mice the glucagon-, somatostatin- and gastrin-secreting endocrine cells "
        "are absent from the epithelium of the glandular stomach and the serotonin-expressing "
        "enterochromaffin cells are dramatically decreased; the same mice show intestinal metaplasia of "
        "the gastric epithelium",
        "source": SRC_LEE,
        "quote": "In ngn3(-/-) mice, glucagon secreting A-cells, somatostatin secreting D-cells, and "
        "gastrin secreting G-cells are absent from the epithelium of the glandular stomach, whereas the "
        "number of serotonin-expressing enterochromaffin (EC) cells is decreased dramatically. In "
        "addition, ngn3(-/-) mice display intestinal metaplasia of the gastric epithelium",
        "not_expressible": "two parts. The EC subtype is a dramatic decrease and not an absence, and this "
        "module has one endocrine lineage and no subtypes; and the intestinal metaplasia means the paper "
        "does NOT report the other gastric lineages as normal, so no kept lineage is registered for this "
        "knockout and the pit, parietal and chief cells the module leaves standing are not a claim",
    },
    {
        "id": "S5",
        "factor": "ASCL1",
        "gates": [EEC],
        "match": "partial",
        "phenotype": "newborn Ascl1-null mice lose the markers of the neurogenin-3-dependent "
        "enteroendocrine cells of the stomach while the enterochromaffin-like cells, mucous cells, chief "
        "cells and parietal cells are normal",
        "source": SRC_HIBDON,
        "quote": "Newborn Ascl1 null mice showed a loss of expression of markers of neurogenin-3-dependent "
        "enteroendocrine cells, with normal expression of enterochromaffin-like cells, mucous cells, chief "
        "cells, and parietal cells",
        "not_expressible": "the enterochromaffin-like cells are explicitly normal and this module's one "
        "endocrine lineage has no subtypes, so deleting it overstates the phenotype; and the readout is a "
        "newborn one, while the module runs to a year",
    },
]

# Excluded, each with the fetched reason. The list is closed and none of these enters the module. The
# first five are the stomach's own pattern: a well-sourced factor whose loss spoils a cell that is
# still there.
EXCLUDED = [
    {
        "factor": "XBP1",
        "source": SRC_HUH,
        "phenotype": "an inducible deletion shrinks the chief (zymogenic) cells fourfold and disrupts "
        "their rough endoplasmic reticulum and granulogenesis",
        "quote": "Tamoxifen-induced deletion of Xbp1 (Xbp1delta) did not affect survival of ZCs but "
        "prevented formation of their structure",
        "reason": "the paper explicitly denies the loss a presence-or-absence gate would predict: the "
        "chief cells survive. Gating the chief route on XBP1 would state an absence the fetched sentence "
        "rules out",
    },
    {
        "factor": "HNF4A",
        "source": SRC_MOORE,
        "phenotype": "an inducible deletion in the adult stomach reduces chief-cell size and ER content, "
        "phenocopying the Xbp1 deletion, and also increases proliferation in the isthmal stem cell zone "
        "and alters mucous neck cell migration",
        "quote": "Loss of HNF4alpha in adult mouse stomach led to reduced ZC size and ER content, "
        "phenocopying previously characterized effects of Xbp1 deletion. However, HNF4alpha(delta/delta) "
        "stomachs also exhibited additional phenotypes including increased proliferation in the isthmal "
        "stem cell zone and altered mucous neck cell migration",
        "reason": "a maturation defect with the cell present, plus an INCREASE in proliferation; the "
        "formalism can state neither. Excluded for the same reason as XBP1 and for HES1's reason in both "
        "gut registrations",
    },
    {
        "factor": "SPDEF",
        "source": SRC_HORST,
        "phenotype": "Spdef-null mice fail to mature the antral mucous gland cells - less Muc6 and Tff2, "
        "fewer secretory granules - and nearly half develop antral mucosal hyperplasia within 1.5 years",
        "quote": "The absence of Spdef impaired terminal maturation of antral mucous gland cells, as "
        "reflected in reduced expression of Muc6 and Tff2 and reduced numbers of secretory granules",
        "reason": "the candidate the task named for the stomach, and it is a maturation defect with the "
        "cell present plus a hyperplasia (an increase). The colonic registration excluded SPDEF for want "
        "of a colon sentence; here the sentence exists and says the wrong kind of thing",
    },
    {
        "factor": "FOXQ1",
        "source": SRC_VERZI,
        "phenotype": "Foxq1 coding mutants lose the MUC5AC mucin backbone and most of the foveolar cell's "
        "secretory vesicles, and the stomach's mucus with them",
        "quote": "Mice carrying Foxq1 coding mutations show virtual absence of mRNA and protein for the "
        "backbone of the major stomach mucin MUC5AC. These observations correspond to a paucity of "
        "foveolar cell secretory vesicles and notable loss of stomach but not intestinal mucus",
        "reason": "stomach-restricted and well sourced, and the best candidate there was for a pit "
        "(foveolar) gate; but the pit cell is present and empty, not absent. This exclusion is why the "
        "lineage that carries most of the stomach's daily turnover has no gate in this module",
    },
    {
        "factor": "KLF4",
        "source": SRC_KATZ,
        "phenotype": "a stomach-specific deletion increases proliferation and alters differentiation, with "
        "aberrant acidic mucins and TFF2/SP-positive cells, and no lineage reported missing",
        "quote": "Klf4 mutant mice survive to adulthood and show increased proliferation and altered "
        "differentiation of their gastric epithelia",
        "reason": "the colon's one colon-specific gate, and in the stomach its phenotype is an increase "
        "with no lineage loss. The formalism cannot express an increase (the reason HES1 is excluded in "
        "both gut registrations)",
    },
    {
        "factor": "ATOH1",
        "source": "none fetched for the stomach",
        "phenotype": "the secretory gate of both gut modules: its loss leaves the secretory lineages "
        "nearly absent in ileum and colon (Lo et al. 2016, carried from the colonic registration)",
        "reason": "no loss-of-function sentence about the STOMACH was found in this lane for ATOH1. "
        "Carrying the gut gate across would be an assumption and not a source - exactly the ground on "
        "which the colonic registration excluded ASCL2. That the secretory gate of both gut modules has "
        "no stomach sentence to quote is a finding about the record, not about the stomach",
    },
    {
        "factor": "SOX2",
        "source": SRC_QUE,
        "phenotype": "Sox2 hypomorphs have an abnormal forestomach with reduced keratinization, ectopic "
        "mucus cells and columnar epithelium, and an esophagus that resembles the trachea",
        "quote": "In all hypomorphic embryos the forestomach has an abnormal phenotype, with reduced "
        "keratinization, ectopic mucus cells and columnar epithelium ... Sox2 plays a second role in "
        "establishing the boundary between the keratinized, squamous esophagus/forestomach and glandular "
        "hindstomach",
        "reason": "the candidate the task named first, and the fetched phenotype is a boundary and "
        "patterning phenotype of the squamous forestomach in hypomorphic embryos, not the loss of a "
        "lineage of the adult glandular stomach this module describes. No adult glandular loss-of-function "
        "sentence was fetched",
    },
    {
        "factor": "GATA4",
        "source": SRC_DELAFOREST,
        "phenotype": "GATA4 controls epithelial morphogenesis in the developing stomach and the "
        "establishment of the glandular columnar epithelium of the hindstomach",
        "quote": "Differential GATA4 expression within these domains coincides with the onset of distinct "
        "tissue morphogenetic events, suggesting a role for GATA4 in diversifying foregut endoderm into "
        "discrete esophageal/FS and HS epithelial tissues",
        "reason": "the second candidate the task named, and the fetched phenotype is developmental "
        "morphogenesis before the gland's lineages exist. This module is an adult steady state and has no "
        "morphogenesis to lose",
    },
    {
        "factor": "LGR5 with AXIN2",
        "source": SRC_BARKER,
        "phenotype": "Lgr5+ cells drive self-renewal in the stomach and build long-lived gastric units in "
        "vitro",
        "reason": "the third candidate the task named. What was fetched is a lineage-tracing and organoid "
        "result: a marker of the antral stem cell, not a loss of function. Deleting Lgr5 is not reported "
        "in what was fetched to remove a gastric lineage, so there is nothing for a gate to state",
    },
    {
        "factor": "POU2F3 (the caveolated/tuft lineage)",
        "source": SRC_KARAM5,
        "phenotype": "caveolated (brush) cells arise in the isthmus and migrate away in both directions",
        "quote": "Little information is available on the dynamics of caveolated cells",
        "reason": "the stomach's tuft lineage, excluded for the reason the tuft lineage is excluded from "
        "both gut modules, and here the fetched sentence says so outright: no share of output and no "
        "residence time exist for it. A lineage with invented numbers is worse than none",
    },
]

# What each knockout must reproduce, written before the module existed. `lost` must reach zero,
# `kept` must stay above KEPT_FLOOR, and `match` is the agreement registered in advance.
MUTANTS = [
    {
        "name": "rbpj",
        "knockout": "RBPJ",
        "lost": MATURE,
        "kept": [],
        "match": "contested",
        "agrees_with": f"{SRC_KIM}, its pharmacological arm (DBZ)",
        "disagrees_with": f"{SRC_KIM}, its genetic arm (Osr1-Cre; RBP-Jk Fl/Fl)",
        "note": "the isthmus compartment does not exist, so no lineage is produced. Registered contested "
        "in advance: this agrees with the virtually abolished proliferation of corpus and antral glands "
        "under the pan-Notch inhibitor DBZ and disagrees with the same paper's genetic arm, where "
        "RBP-Jk-null gastric units are alive and merely 47% fewer. One module, one gate, two experiments "
        "it cannot both satisfy",
    },
    {
        "name": "esrrg",
        "knockout": "ESRRG",
        "lost": [PARIETAL],
        "kept": [PIT, CHIEF, EEC],
        "match": "full",
        "note": "the headline prediction and the one registered full match in this list: the acid-secreting "
        "lineage disappears and the other three do not. The kept lineages rest on the tool, not on an "
        "assumption: Atp4bCre deletes Esrrg only in cells already committed to the parietal lineage, so "
        "the other lineages are untouched in that experiment. The paper's steady-state increase in "
        "proliferation is not expressible and is not checked",
    },
    {
        "name": "mist1",
        "knockout": "MIST1",
        "lost": [CHIEF],
        "kept": [PIT, PARIETAL, EEC],
        "match": "partial",
        "note": "registered partial in advance: Ramsey et al. 2007 report the mature chief cell PRESENT "
        "and structurally abnormal, not absent, and a presence-or-absence gate can only delete it. What "
        "the module does state is the other half of that sentence - the immature transitional cells "
        "accumulate - because the neck compartment loses its only outflow (READINGS R4)",
    },
    {
        "name": "neurog3",
        "knockout": "NEUROG3",
        "lost": [EEC],
        "kept": [],
        "match": "partial",
        "note": "registered partial in advance, and with NO kept lineage: the glandular stomach's A, D and "
        "G cells are absent, which is the match, but the enterochromaffin cells are dramatically decreased "
        "rather than absent and the same mice show intestinal metaplasia of the gastric epithelium. Lee et "
        "al. 2002 therefore do not report the other lineages as normal, so this lane registers none of "
        "them as kept and the pit, parietal and chief cells the module leaves standing are not a claim",
    },
    {
        "name": "ascl1",
        "knockout": "ASCL1",
        "lost": [EEC],
        "kept": [PIT, PARIETAL, CHIEF],
        "match": "partial",
        "note": "the second gate on the one endocrine step, and the one whose paper states the kept "
        "lineages outright: mucous, chief and parietal cells normal. Registered partial in advance because "
        "the enterochromaffin-like cells are normal too and this module has no endocrine subtypes to "
        "spare, and because the readout is a newborn one",
    },
]
KEPT_FLOOR = 1.0e6  # a population the knockout does not touch must stay above this

# Sourced quantities, used as given. The per-unit standing counts are the module's composition and the
# residences are its loss rates; both come from the five-part Karam & Leblond series on the mouse
# corpus, and the human totals come from the body program's own source.
MEASURED = [
    {
        "what": "stomach epithelial cells in the reference adult",
        "value": 1.7e10,
        "source": f"{SRC_MILO} (distilled: data/results/human_cell_turnover.json, row 'Stomach epithelia')",
    },
    {
        "what": "stomach epithelial cells replaced per day",
        "value": 1.3e9,
        "source": SRC_MILO,
    },
    {
        "what": "stomach epithelial lifespan, published mean residence (days)",
        "value": 13.0,
        "source": f"{SRC_MILO}; it is 1.7e10 / 1.3e9 = 13.08 and therefore NOT independent of the two "
        "numbers above",
    },
    {
        "what": "pit (surface mucous) cells per gastric unit, mouse corpus",
        "value": 37.0,
        "source": f"{SRC_KARAM2}: 'The pit region includes an average of 37 pit cells'",
    },
    {
        "what": "pit cell residence (days): the time the migration to the surface takes",
        "value": 3.1,
        "source": f"{SRC_KARAM2}: 'The continuous 3H-thymidine infusion results indicate that this "
        "migration takes 3.1 days on the average'",
    },
    {
        "what": "mucous neck cells per gastric unit, mouse corpus",
        "value": 12.6,
        "source": f"{SRC_KARAM3}: 'The neck region includes an average of 12.6 neck cells'",
    },
    {
        "what": "zymogenic (chief) cells per gastric unit, mouse corpus",
        "value": 67.5,
        "source": f"{SRC_KARAM3}: 'Zymogenic cells, which average 67.5 per unit, further migrate inward'",
    },
    {
        "what": "parietal cells per gastric unit, mouse corpus (3.6 pit + 6.2 isthmus + 5.6 neck + "
        "10.6 base)",
        "value": 26.0,
        "source": f"{SRC_KARAM4}: 'On the average 3.6 cells are in the pit, 6.2 in the isthmus, 5.6 in the "
        "neck, and 10.6 in the base'",
    },
    {
        "what": "parietal cell residence (days)",
        "value": 54.0,
        "source": f"{SRC_KARAM4}: 'The overall turnover time of parietal cells averages 54 days'",
    },
    {
        "what": "isthmal granule-free cell turnover time (days): an UPPER BOUND, and what the one fitted "
        "number is compared with",
        "value": 2.6,
        "source": f"{SRC_KARAM1}: 'These three cell types, as well as pre-parietal cells, are rapidly "
        "renewed, with the turnover times estimated at 3.0 days for pre-neck and pre-parietal cells and "
        "less than 2.6 days for granule-free and pre-pit cells'",
        "used_for": "nothing in the solve: it is the bound the one fitted number is held against, fetched "
        "before the fit",
    },
    {
        "what": "gastric endocrine cell half-life (days), by lineage tracing in the adult antrum",
        "value": 60.0,
        "source": f"{SRC_HIBDON}: 'ASCL1-expressing endocrine cells persisted for several weeks "
        "posttamoxifen labeling with a half-life of approximately 2 months ... confirming that gastric "
        "endocrine cells are long-lived'",
    },
    {
        "what": "mucous neck cell inward migration time (days): a fetched RANGE, 1 to 2 weeks",
        "value": [7.0, 14.0],
        "source": f"{SRC_KARAM3}: 'Neck cells migrate inward in 1-2 weeks from the high through the mid and "
        "low neck segments'",
    },
    {
        "what": "pre-zymogenic cells per gastric unit, mouse corpus",
        "value": 5.0,
        "source": f"{SRC_KARAM3}: 'the pre-zymogenic cells, of which there are 5.0 per unit on the average'",
        "used_for": "nothing: the module folds the pre-zymogenic stage into the neck-to-chief step rather "
        "than giving it a compartment, and this number is recorded so that choice can be read as a choice",
    },
    {
        "what": "pre-parietal cells per gastric unit, mouse corpus",
        "value": 0.6,
        "source": f"{SRC_KARAM4}: 'The pre-parietal cells, of which there are 0.6 per unit on the average'",
        "used_for": "nothing: the module has no pre-parietal compartment",
    },
    {
        "what": "enteroendocrine share of the epithelium (upper bound)",
        "value": 0.01,
        "source": "the bound registered by the small-intestinal module (scripts/intestine_registration.py, "
        "460c46f), itself citing Cheng & Leblond 1974, Am J Anat 141:461 for <1% of the epithelium of "
        "stomach, small intestine and colon; carried, not re-fetched in this lane, exactly as the colonic "
        "registration carried it",
    },
]

# Fitted, and to what. Nothing else in the module is fitted, and no mutant number is fitted at all.
FITTED = [
    {
        "what": "isthmus_residence_days, the isthmal stem (granule-free) cell's turnover time",
        "fitted_to": "the published stomach epithelial count of 1.7e10 cells: the residence that makes the "
        "sum of every compartment's pool equal it",
        "against": f"the fetched UPPER BOUND of 2.6 days for the granule-free cell's turnover time "
        f"({SRC_KARAM1}), fetched before the fit was run. READINGS R2 states in advance the value the "
        "registered inputs imply",
    },
]
INFERRED = [
    "the mouse corpus unit's measured standing counts are used as the human stomach's standing "
    "COMPOSITION, scaled by the published human total. The species carry-over is an inference, and so is "
    "the use of a corpus composition for a stomach that also has an antrum with no parietal or chief "
    "cells; both are registered here and neither is fitted",
    "the mucous neck cell residence, 10.5 days, the midpoint of the fetched 1-to-2-week range",
    "the endocrine cell residence, 86.6 days, obtained from the fetched 2-month half-life as "
    "60 / ln 2, which is the mean residence of the exponential loss law this module actually applies. "
    "The conversion is the inference; the half-life is fetched",
    "the chief (zymogenic) cell residence, 190 days. NO number for it was fetched in this lane: it is the "
    "round figure this module uses for the longest-lived lineage, and READINGS R5 registers in advance "
    "that the designed total and the one fitted number barely move for any value between 50 and 1000 "
    "days, so nothing in this module rests on it",
    "the endocrine standing count, 1.0 cell per gastric unit, inside the <1% bound above and fitted to "
    "nothing. It is the one lineage of the five with no measured count",
    "that the module needs no transit-amplifying compartment: Karam & Leblond part I puts the division in "
    "the isthmus itself (the granule-free cell's labelling index is 32.4%), and no division count per "
    "transit was fetched, so the root divides and there is no amplifying step",
]
NOT_CLAIMED = [
    "the standing composition of the module: pit, parietal, chief and neck pools come back in the "
    "measured per-unit proportions BY CONSTRUCTION, because those proportions are what the module was "
    "built from. They are not a check and must never be reported as one",
    "the published total of 1.7e10 cells: it is the target the one fitted number closes, so the count "
    "assert is the fit closing and not a prediction",
    "the implied number of gastric units: a consequence of the flows, reported and checked against "
    "nothing, because no per-unit cell total was fetched in this lane",
    "the published mean lifespan of 13.0 days: it is the ratio of the two measured numbers above and is "
    "not independent of them",
]

# Statements about the formalism, registered before the module ran. R2, R3 and R5 carry numbers computed
# from the registered inputs alone, so they can be wrong in writing.
READINGS = [
    {
        "id": "R1",
        "says": "one of the five gates is registered full, three partial and one contested, so no run can "
        "report 'five of five'. The result must count the one full match, name the three partial ones with "
        "the part of the phenotype the formalism cannot state, and count the contested one as agreement "
        "with one experiment of one paper and disagreement with another experiment of the same paper",
    },
    {
        "id": "R2",
        "says": "the one fitted number, the isthmal stem residence, will come out at about 2.50 days. That "
        "is arithmetic on the registered inputs and not insight: it is "
        "(1.7e10 - the pools the measured composition and residences already account for) / 1.3e9. The "
        "content is the comparison and not the number: the fetched upper bound for the granule-free cell's "
        "turnover time is 2.6 days, so a knob with one degree of freedom will land inside a bound fetched "
        "before the fit, or it will not",
        "number": 2.50,
        "tolerance": 0.02,
        "bound": 2.6,
    },
    {
        "id": "R3",
        "says": "the stomach's standing composition and its daily output composition diverge by an order "
        "of magnitude, which is the whole arithmetic difference from the colon. The pit cells will be "
        "about 22% of the standing epithelium and about 93% of the daily output; the chief cells about 40% "
        "of the standing epithelium and about 2.8% of the output. In the colonic module the two were "
        "identical by construction, because all three mature residences were one fitted number; here four "
        "residences differ by a factor of 60, so the measured standing share of a lineage says almost "
        "nothing about its share of the turnover",
        "pit_standing": 0.221,
        "pit_output": 0.934,
        "chief_standing": 0.404,
        "chief_output": 0.028,
        "tolerance": 0.01,
    },
    {
        "id": "R4",
        "says": "blocking the chief route (mist1) is NOT the colon's backing-up case. The colonic "
        "secretory progenitor had a second outflow and settled at 51 times its wild-type pool; the gastric "
        "neck compartment has exactly one outflow, so with it gated off nothing leaves and the compartment "
        "accumulates at the wild-type chief inflow rate for as long as the run lasts. Only a DIRECTION is "
        "registered - the neck compartment stands above its wild-type pool - and the magnitude is a "
        "property of the run length and not of the stomach. The direction has the sign of Ramsey et al.'s "
        "markedly increased transitional cells, and the sign is all that is claimed",
    },
    {
        "id": "R5",
        "says": "no number in this module rests on the one residence that was not fetched. Holding "
        "everything else, moving the chief residence from 50 days to 1000 days will move the fitted "
        "isthmal residence by less than 1% (2.504 to 2.497 days) and the designed total by less than 0.1%, "
        "because a longer chief residence raises the chief pool and lowers the chief output and the neck "
        "pool in the same proportion",
        "fitted_at_50_days": 2.504,
        "fitted_at_1000_days": 2.497,
        "tolerance": 0.01,
    },
    {
        "id": "R6",
        "says": "a knockout's effect on the standing count is not its effect on the turnover, which did "
        "not arise in either gut module. Deleting the parietal lineage (esrrg) removes about 16% of the "
        "standing epithelium and about 4% of the daily output; deleting the endocrine lineage (neurog3, "
        "ascl1) removes about 0.6% of the standing epithelium and about 0.09% of the output. The two "
        "phenotypes rank the same way round on either measure, but by factors four times apart",
    },
    {
        "id": "R7",
        "says": "the lineage that carries most of the stomach's turnover has no gate. Every fetched "
        "loss-of-function sentence for a pit (foveolar) factor - FOXQ1, SPDEF, KLF4 - leaves the cell "
        "present and empty or hyperplastic, so no knockout in this module removes the pit cells, and the "
        "93% of daily output they carry is untouched by every mutant. That is a statement about the "
        "stomach's record, not about the stomach",
    },
    {
        "id": "R8",
        "says": "this module has no corpus/antrum axis and no gland geometry. The parietal and chief "
        "lineages are corpus-only in life and unrestricted here; the G cell, which is antral and is what "
        "NEUROG3's fetched sentence names, is inside the one endocrine lineage together with the corpus "
        "endocrine cells; and the isthmus, neck and base are compartments of a flow and not places. The "
        "one antral quantity in MEASURED (the endocrine half-life) is used for the whole endocrine lineage "
        "for that reason",
    },
]

# Bounds the wild-type run must clear. The first is fit closure, not a prediction.
BOUNDS = {
    "total_epithelium_at_1_yr": {"low": 1.445e10, "high": 1.955e10, "fit_closure": True},
    "turnover_per_day": {"low": 1.105e9, "high": 1.495e9, "fit_closure": False},
    "every_mature_population_present": {"fit_closure": False},
}

REGISTRATION = {
    "registered": "2026-10-10",
    "lane": "lane-stomach",
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


def endocrine_residence_days() -> float:
    """The mean residence the fetched 2-month half-life implies under exponential loss."""
    return 60.0 / math.log(2.0)


def main() -> None:
    print(json.dumps(REGISTRATION, indent=2))


if __name__ == "__main__":
    main()
