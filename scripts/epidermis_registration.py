# SPDX-License-Identifier: AGPL-3.0-or-later
"""The registration for the interfollicular-epidermis mechanism module, committed BEFORE its numbers are set.

Area E's fifth human mechanism module and the first outside blood and the gut: after haematopoiesis,
the small-intestinal crypt (scripts/intestine_registration.py), the colonic crypt
(scripts/colon_registration.py) and the gastric gland (scripts/stomach_registration.py), whose discipline
this file follows. docs/AREA-E-NEXT.md holds the written criteria and the comparison that chose the
epidermis, also committed before the module existed.

The list of gating factors below is CLOSED: every entry names a factor, the compartment it gates and a
published loss-of-function phenotype IN THE EPIDERMIS, with the fetched sentence it rests on quoted.
Every citation was fetched in this lane (PubMed E-utilities abstracts, PMC full text, EBI OLS4 for the
Cell Ontology ids, the milo-lab/cellular_turnover repository for the tabulated kinetics, BioNumbers for
Hoath & Leahy's Table I) and not recalled. A factor without a fetched epidermal sentence that a
presence-or-absence gate can state is in EXCLUDED with that as the reason, and `check_factors` refuses a
module that gates on a factor this file does not name.

Three kinds of agreement are registered in advance, as in the colonic and gastric registrations:

- "full": the published epidermal phenotype is the loss of a compartment (a layer) a flow module can
  state; marker details inside a present layer that a population cannot carry are named, not checked;
- "partial": part of the published phenotype is a reduction, a conversion or a sequence in an existing
  tissue, which a module whose gates are presence-or-absence and whose knockouts start at the first cell
  cannot state, so the registered check is the deletion and the rest is named here as not expressible;
- "contested": two fetched EXPERIMENTS disagree and the module's single gate can satisfy at most one of
  them. A contested knockout holding its registered pattern is NOT counted as reproducing the epidermis;
  it is counted as agreeing with the experiment in `agrees_with` and disagreeing with the other.

READINGS are statements about the formalism and the sources registered before any run, three of them
numbers computed from the registered inputs alone, so that what the run shows cannot be explained
afterwards. FITTED names the one quantity fitted to a measurement; nothing about any mutant is fitted.

Licence-compatible use: every number below is a fact quoted with its citation; no text beyond short
quotations and no data file from any source is redistributed, and no run of the module reads the network.

Run `python3 scripts/epidermis_registration.py` to print the registration.
"""

from __future__ import annotations

import json

ORGAN = "human interfollicular epidermis (basal, spinous, granular, cornified)"
MODULE = "data/organisms/human/epidermis.bio"
MUTANT_MODULE = "data/organisms/human/epidermis_mutants.bio"

# The four compartments, bottom to top. The names are the module's cell types; the Cell Ontology ids were
# fetched from EBI OLS4 (exact label search in `cl`).
BASAL = "BasalKeratinocyte"  # CL:0002187 basal cell of epidermis
SPINOUS = "SpinousKeratinocyte"  # CL:0000649 spinous cell of epidermis
GRANULAR = "GranularKeratinocyte"  # CL:0002189 granular cell of epidermis
CORNEOCYTE = "Corneocyte"  # CL:0002153 corneocyte
LAYERS = [BASAL, SPINOUS, GRANULAR, CORNEOCYTE]
NUCLEATED = [BASAL, SPINOUS, GRANULAR]  # what Sender & Milo count; corneocytes are dead and not counted
CELL_ONTOLOGY = {
    BASAL: "CL:0002187",
    SPINOUS: "CL:0000649",
    GRANULAR: "CL:0002189",
    CORNEOCYTE: "CL:0002153",
}
# The step a gate may sit at: a compartment (its existence and everything it does) or, for the one gate on
# self-renewal, the basal compartment's divisions alone.
BASAL_DIVIDE = f"{BASAL}:divide"

SRC_MILO = "Sender & Milo 2021, Nat Med 27:45"
SRC_MILO_SKIN = (
    "Sender & Milo 2021, Nat Med 27:45: github.com/milo-lab/cellular_turnover, Skin_cells.ipynb and "
    "Skin_cells_data.xlsx (sheet epidermal_turnoevr)"
)
SRC_WEINSTEIN = "Weinstein, McCullough & Ross 1984, J Invest Dermatol 82:623"
SRC_HALPRIN = "Halprin 1972, Br J Dermatol 86:14"
SRC_HOATH = "Hoath & Leahy 2003, J Invest Dermatol 121:1440 (Table I, via BioNumbers 109888)"
SRC_BERGSTRESSER = "Bergstresser, Pariser & Taylor 1978, J Invest Dermatol 70:280"
SRC_MILLS = "Mills et al. 1999, Nature 398:708"
SRC_YANG = "Yang et al. 1999, Nature 398:714"
SRC_SEN2010 = "Sen et al. 2010, Nature 463:563"
SRC_BLANPAIN = "Blanpain et al. 2006, Genes Dev 20:3022"
SRC_INGRAHAM = "Ingraham et al. 2006, Nat Genet 38:1335"
SRC_RICHARDSON = "Richardson et al. 2006, Nat Genet 38:1329"
SRC_LEE = "Lee, Kong & Weatherbee 2013, Dev Biol 383:201"
SRC_LI = "Li et al. 1999, Genes Dev 13:1322"
SRC_TAKEDA = "Takeda et al. 1999, Science 284:313"
SRC_GAREUS = "Gareus et al. 2007, Nat Cell Biol 9:461"
SRC_SEN2012 = "Sen et al. 2012, Dev Cell 22:669"
SRC_BUTERA = "Butera et al. 2023, Sci Adv 9:eadg5423"

# The closed list. `gates` are the module steps each factor is required at; the module may require a
# factor only at the steps named here. `quote` is the fetched sentence the gate rests on.
GATES = [
    {
        "id": "E1",
        "factor": "TP63",
        "gates": [BASAL],
        "match": "full",
        "phenotype": "no stratified epidermis at birth: no layer is produced",
        "quote": [
            "The skin of p63-deficient mice does not progress past an early developmental stage: it lacks "
            "stratification and does not express differentiation markers.",
            "The embryonic epidermis of p63-/- mice undergoes an unusual process of non-regenerative "
            "differentiation, culminating in a striking absence of all squamous epithelia and their "
            "derivatives",
        ],
        "source": f"{SRC_MILLS} (PMID 10227293, abstract); {SRC_YANG} (PMID 10227294, abstract)",
    },
    {
        "id": "E2",
        "factor": "DNMT1",
        "gates": [BASAL_DIVIDE],
        "match": "partial",
        "phenotype": "human epidermis regenerated from DNMT1-depleted keratinocytes loses self-renewal and "
        "is "
        "eventually lost, after premature differentiation in the progenitor compartment",
        "quote": [
            "In tissue, DNMT1 depletion led to exit from the progenitor cell compartment, premature "
            "differentiation and eventual tissue loss.",
            "By two weeks, DNMT1i epidermis was markedly hypoplastic ( Fig. 2a ) without increased cell "
            "death "
            "( Supplementary Fig. 2 ), and only 34% of DNMT1-deficient epidermal tissues survived to week 3 "
            "compared to 100% for controls ( Fig. 2a,b ).",
        ],
        "source": f"{SRC_SEN2010} (PMID 20081831; PMC3050546, full text): shRNA in human keratinocytes "
        "regenerating epidermis on immune-deficient mice",
        "not_expressible": "the premature differentiation (K1 in the basal layer) and the sequence itself: a "
        "knockout here acts from the first cell, so the tissue never forms rather than forming and being "
        "lost; "
        "the knockdown is partial and 34% of tissues survived to week 3, so total loss at one year is the "
        "module's endpoint, not the paper's",
    },
    {
        "id": "E3",
        "factor": "RBPJ",
        "gates": [SPINOUS],
        "match": "partial",
        "phenotype": "spinous layers absent; fewer granular cells; basal markers relatively normal",
        "quote": [
            "Spinous layers are absent in RBP-J conditional null epidermis and expanded when Notch1 "
            "signaling is "
            "elevated transgenically in epidermis.",
            "Additionally, KO epidermis displayed fewer granular layer cells, which were also distinguished "
            "by a "
            "reduction in keratohyalin granules.",
            "In contrast to the repression of spinous and granular markers, basal markers appeared "
            "relatively "
            "normal in RBP-J cKO epidermis ( Fig. 3 ).",
        ],
        "source": f"{SRC_BLANPAIN} (PMID 17079689; PMC1620020, full text): K14-Cre conditional null",
        "not_expressible": "'fewer granular layer cells' is overstated by a linear chain that can only "
        "delete: "
        "with the spinous compartment gone, the granular and cornified compartments above it are gone too",
    },
    {
        "id": "E4",
        "factor": "IRF6",
        "gates": [GRANULAR],
        "match": "full",
        "phenotype": "basal layer present, spinous layer greatly expanded, granular and cornified layers "
        "absent",
        "quote": [
            "In the epidermis from null embryos, the basal layer was present along with a greatly expanded "
            "spinous layer, but the outer two layers appeared absent.",
        ],
        "source": f"{SRC_INGRAHAM} (PMID 17041601; PMC2082114, full text): Irf6 null; agreeing, "
        f"{SRC_RICHARDSON} (PMID 17041603, abstract): a homozygous missense allele with 'a "
        f"hyperproliferative "
        "epidermis that fails to undergo terminal differentiation'",
        "not_expressible": "the proliferation inside the expanded spinous layer and the ectopic K14 and p63 "
        "there; the spinous compartment here does not divide, and only the direction of its expansion is "
        "registered",
    },
    {
        "id": "E5",
        "factor": "KDF1",
        "gates": [GRANULAR],
        "match": "full",
        "phenotype": "granular and cornified markers not induced; basal and spinous markers expanded",
        "quote": [
            "Markers of granular and cornified layers (D, Loricrin in red; E, Filaggrin in red) are not "
            "induced "
            "in Kdf1 GT homozygote embryos and in shd-Kdf1 GT transheterozygotes.",
            "Basal (B, Keratin 5 in red) and spinous (C, Keratin 10 in red) cell markers are expanded in "
            "keratinocytes in Kdf1 GT homozygote embryos as well as in Kdf1 GT/shd transheterozygotes.",
            "However, in shd mutants, the most distal keratinocytes still retained their nuclei and the "
            "epidermis showed no signs of a cornified layer ( Figure 1G ).",
        ],
        "source": f"{SRC_LEE} (PMID 24075906; PMC3841015, full text): the gene-trap allele "
        "1810019J16Rik tm1a(EUCOMM)Wtsi, which the paper uses to show shd is a loss of function",
        "not_expressible": "the suprabasal proliferation and the mixed basal/spinous identity",
    },
    {
        "id": "E6",
        "factor": "CHUK",
        "gates": [GRANULAR],
        "match": "contested",
        "agrees_with": f"{SRC_LI}: germline IKK1 (IKKalpha, CHUK) null",
        "disagrees_with": f"{SRC_GAREUS}: keratinocyte-restricted ablation",
        "phenotype": "germline null: granular and cornified layers absent, single basal layer unchanged, "
        "keratin-10-positive suprabasal cells increased; keratinocyte-restricted ablation: a normally "
        "differentiated stratified epidermis",
        "quote": [
            "Dead cells devoid of nuclei, characteristic of the cornified layer (Fig. 3 C), were not "
            "observed in "
            "mutant mouse skin. The granular layer was also not obvious as the granular cells with "
            "distinctive "
            "keratohyalin granules were not present.",
            "Although anti-keratin 14 immunostaining showed no difference of the single basal layer in both "
            "wild-type and mutant skin (Fig. 3 C,D), keratin 10-positive cells increased dramatically in the "
            "mutant suprabasal layer (SB) of skin compared with littermates (Fig. 3 E,F).",
            "mice with epidermal keratinocyte-specific IKK1 ablation (hereafter referred to as IKK1(EKO)) "
            "develop a normally differentiated stratified epidermis, demonstrating that the function of "
            "IKK1 in "
            "inducing epidermal differentiation is not keratinocyte-autonomous.",
        ],
        "source": f"{SRC_LI} (PMID 10346820; PMC316728, full text); agreeing, {SRC_TAKEDA} (PMID 10195895, "
        "abstract: 'highly proliferative with dysregulated epidermal differentiation'); disagreeing, "
        f"{SRC_GAREUS} (PMID 17351639, abstract)",
        "not_expressible": "where the factor acts: a knockout here removes it from the whole module from the "
        "first cell, which is the germline null's setting; a keratinocyte-restricted ablation cannot be "
        "stated",
    },
    {
        "id": "E7",
        "factor": "ZNF750",
        "gates": [GRANULAR],
        "match": "contested",
        "agrees_with": f"{SRC_SEN2012}: knockdown in regenerated organotypic human epidermis",
        "disagrees_with": f"{SRC_BUTERA}: mouse Zfp750 null",
        "phenotype": "human knockdown: keratins 1 and 10 kept, the late (granular) differentiation genes not "
        "induced; mouse null: granular layer and stratum corneum present, filaggrin granules altered, "
        "loricrin "
        "unchanged, barrier lost",
        "quote": [
            "Although leaving early differentiation gene markers keratin 1 and 10 unaffected ( Figure "
            "S1B-S1C ), "
            "ZNF750 loss abolished induction of late differentiation genes, including FLG , LOR , LCE3A , "
            "LCE3D , "
            "and SPRR1A ( Figures 1D-1E ).",
            "ZFP750 −/− image reveals a more densely packed SC and perturbation of filaggrin (F) granule "
            "distribution.",
            "Loricrin levels were unchanged both at protein and mRNA level",
        ],
        "source": f"{SRC_SEN2012} (PMID 22364861; PMC3306457, full text); disagreeing, {SRC_BUTERA} "
        "(PMID 37115925; PMC10146900, full text)",
        "not_expressible": "a gene-expression readout is read here as the absence of the layer those genes "
        "define; the human tissue is regenerated, not a body",
    },
]

EXCLUDED = [
    {
        "factor": "NOTCH1",
        "source": "Rangarajan et al. 2001, EMBO J 20:3427 (PMID 11432830, abstract); Blanpain et al. 2006 "
        "(PMC1620020) describing it",
        "quote": "This contrasts strongly with Notch1 cKO epidermis, which is typified by an increase rather "
        "than decrease of granular marker expression and no noticeable defect in spinous layers",
        "reason": "the phenotype is hyperplasia and an increase; no layer is lost, so a gate would state a "
        "loss "
        "no paper reports. The spinous gate is RBPJ, whose loss does remove the layer",
    },
    {
        "factor": "ADAM10",
        "source": "Weber et al. 2011, Development 138:495 (PMID 21205794, abstract)",
        "quote": "A reduction of spinous layers, not associated with differences in either proliferation or "
        "apoptosis, indicates that loss of Adam10 triggers a premature differentiation of spinous "
        "keratinocytes.",
        "reason": "a reduction with premature differentiation, not a loss: gating the spinous compartment on "
        "ADAM10 would delete the granular and cornified layers that premature differentiation produces",
    },
    {
        "factor": "RIPK4",
        "source": "Holland et al. 2002, Curr Biol 12:1424 (PMID 12194825, abstract); Urwyler-Rosselet et al. "
        "2018, J Invest Dermatol 138:1268 (PMID 29317263, abstract)",
        "quote": "In contrast to RIPK4 full KO epidermis, RIPK4EKO epidermis was normally stratified",
        "reason": "the full knockout's layer phenotype is stated in no sentence fetched in this lane (the "
        "2002 "
        "abstract says only 'abnormal epidermal differentiation' and the full text was not open to it), and "
        "the epidermis-specific knockout stratifies normally; a contested gate needs two fetched sentences",
    },
    {
        "factor": "KLF4",
        "source": "Segre, Bauer & Fuchs 1999, Nat Genet 22:356 (PMID 10431239, abstract)",
        "quote": "Loss of the barrier occurs without morphological and biochemical alterations to the "
        "well-known structural features of epidermis that are essential for mechanical integrity.",
        "reason": "the layers are present and the cornified envelope is spoiled; a presence gate has nothing "
        "to say about a spoiled cell",
    },
    {
        "factor": "GRHL3",
        "source": "Ting et al. 2005, Science 308:411 (PMID 15831758, abstract); Yu et al. 2006, Dev Biol "
        "299:122 (PMID 16949565, abstract)",
        "quote": "The epidermis of Get-1(-/-) mice showed a severe barrier function defect associated with "
        "impaired differentiation of the epidermis, including defects of the stratum corneum, extracellular "
        "lipid composition and cell adhesion in the granular layer.",
        "reason": "a defective stratum corneum and granular layer, both present: a maturation defect",
    },
    {
        "factor": "TGM1",
        "source": "Matsuki et al. 1998, Proc Natl Acad Sci USA 95:1044 (PMID 9448282, abstract)",
        "quote": "In their stratum corneum, degradation of nuclei and keratohyalin F-granules was "
        "incomplete and "
        "cell envelope assembly was defective.",
        "reason": "the stratum corneum is present and defective",
    },
    {
        "factor": "PRDM1",
        "source": "Magnusdottir et al. 2007, Proc Natl Acad Sci USA 104:14988 (PMID 17846422, abstract)",
        "quote": "conditional knockout mice have delayed barrier formation as embryos, enlarged granular "
        "layer "
        "cells and corneocytes, and a morphologically abnormal cornified layer.",
        "reason": "the granular-to-cornified transition is delayed and abnormal, not absent",
    },
    {
        "factor": "SFN",
        "source": "Herron et al. 2005, Nat Genet 37:1210 (PMID 16200063, abstract)",
        "quote": "Er/+ mice expressed a truncated Sfn protein, which probably contributes to the defects in "
        "Er/Er "
        "and Er/+ epidermis",
        "reason": "the repeated-epilation allele makes a truncated protein with a phenotype in "
        "heterozygotes, so "
        "it is not a loss of function, and no null was fetched",
    },
    {
        "factor": "OVOL1",
        "source": "Nair et al. 2006, J Cell Biol 173:253 (PMID 16636146, abstract)",
        "quote": "In its absence, developing epidermis fails to properly restrict the proliferative "
        "potential of "
        "progenitor cells",
        "reason": "an increase in proliferation, not a loss of a layer",
    },
    {
        "factor": "CEBPA with CEBPB; MAF with MAFB; HDAC1 with HDAC2",
        "source": None,
        "quote": None,
        "reason": "considered and not fetched in this lane: each published loss is of a PAIR, and a gate of "
        "one "
        "factor cannot state redundancy, so no citation is claimed",
    },
]

# What each knockout must reproduce. `lost` must be exactly 0 at one year, `below_one` must be below one
# cell, `kept` must stay at or above KEPT_FLOOR; `extra` is a registered direction checked against the
# wild-type design.
MUTANTS = [
    {
        "name": "tp63",
        "knockout": "TP63",
        "lost": [SPINOUS, GRANULAR, CORNEOCYTE],
        "seed_only": BASAL,
        "kept": [],
        "match": "full",
        "note": "the basal compartment never forms: it holds the single seed cell and produces nothing",
    },
    {
        "name": "dnmt1",
        "knockout": "DNMT1",
        "lost": [],
        "below_one": LAYERS,
        "kept": [],
        "match": "partial",
        "note": "registered partial in advance: without self-renewal the seed cell's descendants drain away "
        "and every compartment is below one cell at one year; the published premature differentiation and "
        "the "
        "loss of an already-formed tissue are not expressible",
    },
    {
        "name": "rbpj",
        "knockout": "RBPJ",
        "lost": [SPINOUS, GRANULAR, CORNEOCYTE],
        "kept": [BASAL],
        "match": "partial",
        "note": "registered partial in advance: 'spinous layers are absent' is the match; 'fewer granular "
        "layer "
        "cells' is overstated as none by a chain that can only delete",
    },
    {
        "name": "irf6",
        "knockout": "IRF6",
        "lost": [GRANULAR, CORNEOCYTE],
        "kept": [BASAL, SPINOUS],
        "extra": f"{SPINOUS} above twice its wild-type pool (the published greatly expanded spinous layer)",
        "match": "full",
        "note": "the headline pattern: basal present, spinous expanded, the outer two layers absent",
    },
    {
        "name": "kdf1",
        "knockout": "KDF1",
        "lost": [GRANULAR, CORNEOCYTE],
        "kept": [BASAL, SPINOUS],
        "extra": f"{SPINOUS} above twice its wild-type pool (the published expanded spinous markers)",
        "match": "full",
        "note": "the same pattern as IRF6, from an independent gene",
    },
    {
        "name": "chuk",
        "knockout": "CHUK",
        "lost": [GRANULAR, CORNEOCYTE],
        "kept": [BASAL, SPINOUS],
        "extra": f"{SPINOUS} above twice its wild-type pool (the published increase of keratin-10-positive "
        "suprabasal cells)",
        "match": "contested",
        "note": "registered contested in advance: agrees with the germline null, disagrees with the "
        "keratinocyte-restricted ablation; counted as agreement, not reproduction",
    },
    {
        "name": "znf750",
        "knockout": "ZNF750",
        "lost": [GRANULAR, CORNEOCYTE],
        "kept": [BASAL, SPINOUS],
        "extra": f"{SPINOUS} above twice its wild-type pool (a consequence of the formalism, not a published "
        "phenotype of this knockdown: the human paper reports keratins 1 and 10 'unaffected')",
        "match": "contested",
        "note": "registered contested in advance: agrees with the human knockdown, disagrees with the mouse "
        "null; counted as agreement, not reproduction",
    },
]
KEPT_FLOOR = 1.0e6  # a population the knockout does not touch must stay above this
EXPANDED_FACTOR = 2.0  # the registered `extra`: the spinous pool above this multiple of its wild-type design

# Sourced quantities, used as given.
MEASURED = [
    {
        "what": "nucleated epidermal cells in the reference adult (corneocytes not counted)",
        "value": 1.2e11,
        "source": f"{SRC_MILO} (distilled: data/results/human_cell_turnover.json); Skin_cells.ipynb: 'As "
        "corneocytes are essentially dead cells, we don't count them.'",
    },
    {
        "what": "epidermal cells replaced per day",
        "value": 3.7e9,
        "source": f"{SRC_MILO}: the count over the 31-day lifespan of the viable compartments",
    },
    {
        "what": "lifespan in the viable (nucleated) compartments, days",
        "value": 31.0,
        "source": f"{SRC_MILO_SKIN}: the mean of Halprin 1972 (minimal), Bergstresser & Taylor 1977 and "
        "Weinstein et al. 1984, Halprin's maximal estimate removed",
    },
    {
        "what": "basal (germinative) compartment time, days: the proliferative compartment's cell-cycle "
        "duration of 311 h",
        "value": 311.0 / 24.0,
        "source": f"{SRC_WEINSTEIN} (PMID 6725985, abstract: 'The proliferative compartment of normal "
        "epidermis has a cell cycle duration (Tc) of 311 h'), tabulated as 12.958 d for the germinative cell "
        f"layer in {SRC_MILO_SKIN}",
    },
    {
        "what": "stratum corneum turnover time, days",
        "value": 14.0,
        "source": f"{SRC_HOATH}: 'Average turnover time for human stratum corneum' 14 d (Rothberg et al. "
        "1961; Weinstein & Van Scott 1965; Jansen et al. 1974); 14 d in every row of the Sender & Milo table",
    },
    {
        "what": "one flux through every compartment: the turnover of the proliferative, "
        "viable-differentiated "
        "and cornified compartments, cells/day/mm2",
        "value": [1246, 1417, 1490],
        "source": f"{SRC_WEINSTEIN}, abstract: 'Expected kinetic homeostasis in the epidermis is confirmed "
        f"by "
        'the very similar "turnover" rates in each of the compartments that are, respectively, 1246, 1417, '
        "and 1490 cells/day/mm2 surface area. The mean epidermal turnover time of the entire tissue is 39 "
        "days.'",
    },
    {
        "what": "transit through the viable differentiated compartment (spinous plus granular), days: the "
        "published span the fitted number is compared with, not an input",
        "value": {"weinstein_1984": 12.0, "halprin_1972_minimal": 26.0, "halprin_1972_maximal": 42.0},
        "source": f"{SRC_MILO_SKIN} tabulating {SRC_WEINSTEIN} and {SRC_HALPRIN}",
    },
    {
        "what": "corneocytes per mm2 of human skin",
        "value": 17778,
        "source": SRC_HOATH,
    },
    {
        "what": "cellular rate of formation of stratum corneum (human), cells per mm2 per second",
        "value": 0.0147,
        "source": SRC_HOATH,
    },
    {
        "what": "total body skin surface area, m2 (head and trunk 0.8, upper arms and thighs 0.55, forearms "
        "and lower legs 0.55)",
        "value": 1.9,
        "source": f"ICRP Publication 89 (2002), as tabulated in {SRC_MILO_SKIN} (sheet Skin_properties)",
    },
    {
        "what": "nucleated epidermal cells per mm2: not independent of the count above, which Sender & Milo "
        "built from these densities and the area",
        "value": {"bergstresser_1978": [45200, 48700], "bauer_2001": 75346},
        "source": f"{SRC_BERGSTRESSER} (abstract: 'The mean densities of nucleated cells per 10(4) mu2 "
        f"surface "
        f"area were 452, 483, and 487'); Bauer et al. 2001, J Invest Dermatol 116:313, via {SRC_HOATH}",
    },
]

# Fitted, and to what. Nothing else in the module is fitted, and no mutant number is fitted at all.
FITTED = [
    {
        "what": "suprabasal_residence_days: the time a cell spends in the viable suprabasal layers (spinous "
        "plus granular)",
        "fitted_to": "the published 1.2e11 nucleated epidermal cells: the residence that makes the basal, "
        "spinous and granular pools sum to it at the published output of 3.7e9 per day",
        "against": "the published transits of the viable differentiated compartment, 12 d (Weinstein et al. "
        "1984) to 26-42 d (Halprin 1972), fetched before the fit (READINGS R2)",
    },
]
INFERRED = [
    "the granular residence, 2.0 days, with no primary citation fetched in this lane: about two cell layers "
    "at the 1.15 corneum layers formed per day that Hoath & Leahy tabulate, an analogy and not a "
    "measurement. "
    "The spinous residence is the fitted suprabasal residence less this. READINGS R5 registers that nothing "
    "but the split of the suprabasal pool between the two layers depends on it",
    "no death in the viable layers: every cell born in the basal layer leaves through the stratum corneum, "
    "on the strength of Weinstein et al.'s three near-equal compartment turnovers",
    "the basal layer's renewal is capped near its designed pool, because it is a single layer on a fixed "
    "area of basement membrane; READINGS R6 registers what that cap does under RBPJ",
]
NOT_CLAIMED = [
    "the basal pool (the output times the 311 h compartment time, about 4.8e10, 40% of the nucleated "
    "cells): a consequence of the flows, checked against no basal density, because none was fetched",
    "the published 31-day viable lifespan: it is the ratio of the two counted numbers above and is not "
    "independent of them; the module's implied mean is their ratio, 32.4 days",
]

# Statements about the formalism and the sources, registered before the module ran. R2, R3 and R4 carry
# numbers computed from the registered inputs alone, so they can be wrong in writing.
READINGS = [
    {
        "id": "R1",
        "says": "three of the seven knockouts are registered full (TP63, IRF6, KDF1), two partial (DNMT1, "
        "RBPJ) and two contested (CHUK, ZNF750), so no run can report 'seven of seven reproduced'. The "
        "result "
        "must count the full ones as reproductions, name the part of each partial one the formalism cannot "
        "state, and count each contested one as agreement with one experiment and disagreement with another",
    },
    {
        "id": "R2",
        "says": "the one fitted number, the suprabasal residence, will come out at about 19.47 days. That is "
        "arithmetic on the registered inputs and not insight: it is 1.2e11 / 3.7e9 - 311 / 24. The content "
        "is "
        "the comparison: the fetched transits of the viable differentiated compartment are 12 d (Weinstein) "
        "and 26 to 42 d (Halprin), so the fit lands between the two studies and inside neither; a residence "
        "outside 12 to 42 d would have meant the published count and output contradict the published "
        "transits",
        "number": 19.47,
        "tolerance": 0.02,
        "bound": [12.0, 42.0],
    },
    {
        "id": "R3",
        "says": "the module's corneocyte pool is not fitted and will disagree with the corneocyte density. "
        "At "
        "the published output and the 14-day stratum corneum it is 3.7e9 x 14 = 5.18e10 corneocytes, while "
        "17,778 per mm2 over 1.9 m2 is 3.38e10: the module stands 1.53 times above. That ratio is exactly "
        "the "
        "ratio of Sender & Milo's 3.7e9 per day to the 0.0147 per mm2 per second (2.41e9 per day over 1.9 "
        "m2) "
        "of Hoath & Leahy, because both routes use 14 days. It is a disagreement between two sources about "
        "the "
        "flux, inherited, and not a property of the epidermis the module found",
        "corneocyte_pool": 5.18e10,
        "corneocyte_from_density": 3.38e10,
        "ratio": 1.53,
        "tolerance": 0.02,
    },
    {
        "id": "R4",
        "says": "blocking the granular route (IRF6, KDF1, CHUK, ZNF750) leaves the spinous compartment with "
        "no "
        "outflow, so it accumulates at the wild-type inflow, about 3.7e9 per day, for as long as the run "
        "lasts. "
        "Only a DIRECTION is registered - the spinous pool stands above twice its wild-type design at one "
        "year "
        "- and that direction is the sign of the published expanded spinous layer (IRF6, KDF1, CHUK). The "
        "size, "
        "of the order of 18 times, is a property of the run's length and not of the epidermis. For ZNF750 "
        "the "
        "same expansion contradicts the human paper's 'unaffected' keratins 1 and 10 and is recorded as such",
        "spinous_wild_type": 6.47e10,
        "order_of_fold_at_one_year": 18.0,
    },
    {
        "id": "R5",
        "says": "no knockout outcome and no total rests on the one inferred residence. Moving the granular "
        "residence from 1 to 4 days moves only the split of the suprabasal pool between the spinous and "
        "granular layers; the fitted suprabasal residence, the nucleated total, the corneocyte pool and "
        "every "
        "presence-or-absence check are identical, because the fit is of the sum",
    },
    {
        "id": "R6",
        "says": "under RBPJ the basal compartment loses its only outflow. Because its renewal is capped "
        "near its "
        "designed pool (INFERRED), it will stand within 10% of its wild-type pool instead of growing without "
        "bound. That has the sign of the published 'basal markers appeared relatively normal' and of the "
        "published 'reduction in proliferation', but the cap is this module's assumption, stated here, so "
        "the "
        "agreement is not evidence",
    },
    {
        "id": "R7",
        "says": "TP63 and DNMT1 cannot be told apart by what is left: in both the tissue never forms, "
        "because a "
        "knockout here acts from the first cell. They differ in what they gate (the compartment, or only its "
        "divisions), not in any count this module can read",
    },
    {
        "id": "R8",
        "says": "this module has no appendages (hair follicle, sebaceous gland), no regional thickness, no "
        "space and no barrier. The layers are compartments of a flow, not places; the limb and craniofacial "
        "phenotypes of TP63, IRF6 and CHUK, and every barrier phenotype, are outside it",
    },
]

# Bounds the wild-type run must clear. The first is fit closure, not a prediction.
BOUNDS = {
    "nucleated_at_1_yr": {"low": 1.02e11, "high": 1.38e11, "fit_closure": True},
    "turnover_per_day": {"low": 3.145e9, "high": 4.255e9, "fit_closure": False},
    "every_layer_present": {"fit_closure": False},
}

FALSIFIER = (
    "a registered knockout that fails its pass rule is reported as failed, verbatim, and not tuned; the wild "
    "type fails if the nucleated total misses 1.2e11 by more than 15% or the daily turnover misses 3.7e9 by "
    "more than 15%; R2 fails if the fitted residence leaves 12 to 42 days; R3 fails if the corneocyte ratio "
    "is not 1.53 within 0.02"
)

REGISTRATION = {
    "registered": "2026-10-10",
    "lane": "lane-module3",
    "organ": ORGAN,
    "module": MODULE,
    "mutant_module": MUTANT_MODULE,
    "layers": LAYERS,
    "nucleated": NUCLEATED,
    "cell_ontology": CELL_ONTOLOGY,
    "gates": GATES,
    "excluded": EXCLUDED,
    "mutants": MUTANTS,
    "kept_floor": KEPT_FLOOR,
    "expanded_factor": EXPANDED_FACTOR,
    "measured": MEASURED,
    "fitted": FITTED,
    "inferred": INFERRED,
    "not_claimed": NOT_CLAIMED,
    "readings": READINGS,
    "bounds": BOUNDS,
    "falsifier": FALSIFIER,
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


def suprabasal_residence_days() -> float:
    """R2's number, from the registered inputs alone: the residence the fit must return."""
    return 1.2e11 / 3.7e9 - 311.0 / 24.0


def main() -> None:
    print(json.dumps(REGISTRATION, indent=2))


if __name__ == "__main__":
    main()
