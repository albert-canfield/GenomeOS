# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""What the literature holds for each of the 21 rule numbers of `data/demo/gastrulation.bio`.

Written AFTER every source was fetched and read, which is why it is a module of its own. The
question, the population, the classes, the cascade, the two other axes, the keys every finding must
carry and the search plan are in `rule_number_sources`, and they were committed before one source
was fetched (`data/results/rule_number_sources_registration.json` at commit 4392d24, code revision
7c80d18). That module's `FINDINGS`, `SOURCES_FETCHED`, `ATTRIBUTION_FINDINGS`, `CONTRADICTIONS` and
`VERDICT` are the tables AS REGISTERED - empty - and they stay empty there, so that
`registration()` keeps reproducing the committed artefact and the emptiness it recorded stays
checkable. The filled tables live here, and the census writer loads them.

Nothing here is taken from memory. Every entry names the URL that was fetched, what was retrieved
from it, what that source claims and - separately, because this is where a predecessor was weakest -
what it does NOT claim. A WebSearch result is never evidence for a class: a search only pointed at a
source, which was then fetched and read.

Nothing here changes a number. A value change would be a MODEL change and would need its own
registration and its own re-run.
"""

from __future__ import annotations

from typing import Any

from genomeos.provenance.rule_number_sources import CLASSES, slots

__all__ = [
    "CONTRADICTIONS",
    "DECLINED_MAPPINGS",
    "FINDINGS",
    "SEARCHES_RUN",
    "SOURCES_FETCHED",
    "VERDICT",
    "attribution_findings",
    "findings",
]

# --- the fetches ---------------------------------------------------------------------------------

_EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

#: Every source reached, one entry per URL, with what was retrieved from it. A WebSearch result is
#: NOT a source here: only these fetches are used as evidence, and the queries are kept separately.
SOURCES_FETCHED: tuple[dict[str, Any], ...] = (
    {
        "key": "conlon1994_epmc",
        "url": f"{_EPMC}?query=AUTH:%22Conlon+FL%22+AND+PUB_YEAR:1994+AND+JOURNAL:%22Development%22"
        "&resultType=core&format=json",
        "retrieved": (
            "the Europe PMC core record and the FULL abstract of Conlon et al., `A primary "
            "requirement for nodal in the formation and maintenance of the primitive streak in the "
            "mouse`, Development 120:1919-1928, PMID 7924997, DOI 10.1242/dev.120.7.1919. The "
            "citation in the program is correct to the volume and first page."
        ),
        "claims": (
            "the 413.d proviral insertion is a loss-of-function mutation of nodal; `413.d mutant "
            "embryos show no morphological evidence for the formation of a primitive streak`; `about "
            "25% of mutant embryos do form randomly positioned patches of cells of a posterior "
            "mesodermal character`; nodal RNA expression is described by stage and position."
        ),
        "does_not_claim": (
            "the abstract does not name Brachyury, Bra or the T gene at all, and gives no dose, "
            "concentration, EC50, half-maximal value, Hill coefficient or quantified strength. It is "
            "a loss-of-function and expression study, not a dose-response."
        ),
    },
    {
        "key": "conlon1994_publisher",
        "url": "https://journals.biologists.com/dev/article/120/7/1919/38439/"
        "A-primary-requirement-for-nodal-in-the-formation",
        "retrieved": (
            "the publisher's article page, reached through DOI 10.1242/dev.120.7.1919. ONLY THE "
            "ABSTRACT is available: the full text is behind a paywall, so the full text of this "
            "citation was NOT read and nothing is claimed from it."
        ),
        "claims": "the same abstract as the Europe PMC record, retrieved independently.",
        "does_not_claim": (
            "neither Brachyury nor the T gene appears in what was retrieved, and no number of any "
            "kind relating nodal level to mesoderm or to Brachyury appears in it."
        ),
    },
    {
        "key": "kanaiazuma2002_epmc",
        "url": f"{_EPMC}?query=AUTH:%22Kanai-Azuma+M%22+AND+PUB_YEAR:2002&resultType=core&format=json",
        "retrieved": (
            "the Europe PMC core record and the FULL abstract of Kanai-Azuma et al., `Depletion of "
            "definitive gut endoderm in Sox17-null mutant mice`, Development 129:2367-2379, PMID "
            "11973269, DOI 10.1242/dev.129.10.2367. The citation in the program is correct to the "
            "volume and first page. The fetch was asked explicitly whether `Nodal` or `nodal` "
            "appears anywhere in the abstract: it does not."
        ),
        "claims": (
            "`Sox17(-/-) mutant embryos are deficient of gut endoderm`; the earliest defect is "
            "reduced occupancy of the definitive endoderm in the prospective mid- and hindgut; "
            "Sox17-null ES cells contribute to ectoderm and mesoderm but are excluded from mid- and "
            "hindgut endoderm; Sox17's upstream factors in Xenopus and zebrafish are Mixer and "
            "Casanova."
        ),
        "does_not_claim": (
            "the abstract does not contain the word Nodal at all, so it says nothing about Nodal "
            "activating SOX17 - the relation the rule it is cited for states. It is evidence that "
            "Sox17 is REQUIRED for gut endoderm, which is a statement about what Sox17 does "
            "downstream, not about what induces it. It gives no dose, threshold, EC50 or Hill "
            "coefficient."
        ),
    },
    {
        "key": "kanaiazuma2002_publisher",
        "url": "https://journals.biologists.com/dev/article/129/10/2367/41701/"
        "Depletion-of-definitive-gut-endoderm-in-Sox17-null",
        "retrieved": (
            "the publisher's article page, reached through DOI 10.1242/dev.129.10.2367. ONLY THE "
            "ABSTRACT is available: the full text is behind a paywall, so the full text of this "
            "citation was NOT read and nothing is claimed from it. `Nodal` does not appear in what "
            "was retrieved, which is the second independent fetch to find that."
        ),
        "claims": "the same abstract as the Europe PMC record, retrieved independently.",
        "does_not_claim": "any numerical dose, threshold, EC50 or Hill coefficient for Nodal on Sox17.",
    },
    {
        "key": "thomson2011_epmc",
        "url": f"{_EPMC}?query=AUTH:%22Thomson+M%22+AND+PUB_YEAR:2011+AND+JOURNAL:%22Cell%22"
        "&resultType=core&format=json",
        "retrieved": (
            "the Europe PMC core record and the FULL abstract of Thomson et al., `Pluripotency "
            "factors in embryonic stem cells regulate differentiation into germ layers`, Cell "
            "145:875-889, PMID 21663792, PMCID PMC5603300. The citation in the program is correct "
            "to the volume and first page."
        ),
        "claims": (
            "`Oct4 suppresses neural ectodermal differentiation and promotes mesendodermal "
            "differentiation; Sox2 inhibits mesendodermal differentiation and promotes neural "
            "ectodermal differentiation`; differentiation signals modulate Oct4 and Sox2 protein "
            "levels asymmetrically."
        ),
        "does_not_claim": (
            "the abstract names neither Brachyury nor the T gene, and contains no Hill coefficient, "
            "EC50 or half-maximal value."
        ),
    },
    {
        "key": "thomson2011_fulltext",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5603300/",
        "retrieved": (
            "the full text in PMC, read for Brachyury, for the direction of each repression, and for "
            "any number."
        ),
        "claims": (
            "Brachyury appears throughout as a mesendoderm marker and readout; `Sox2 specifically "
            "represses Brachyury and the mesendodermal lineage`; `Sox2 specifically represses only "
            "the mesendodermal fate`; and the repressor of the neural program is OCT4: `Oct4 "
            "specifically represses Sox1 and the NE lineage`."
        ),
        "does_not_claim": (
            "it does NOT state that Brachyury represses Sox2 - the rule it is cited for on line 21. "
            "It reports no Hill coefficient, EC50, half-maximal concentration, dissociation constant "
            "or quantified fold-repression for Sox2 on mesendoderm genes, and its protein levels are "
            "in `units normalized to the population mean`, which is arbitrary units, not absolute "
            "concentration. No dose-response curve in absolute concentration is reported."
        ),
    },
    {
        "key": "lolas2014_epmc",
        "url": f"{_EPMC}?query=AUTH:%22Lolas+M%22+AND+PUB_YEAR:2014&resultType=core&format=json",
        "retrieved": (
            "the Europe PMC core record and the FULL abstract of Lolas et al., `Charting "
            "Brachyury-mediated developmental pathways during early mouse embryogenesis`, PNAS "
            "111:4478-4483, PMID 24616493, PMCID PMC3970479. The citation in the program is correct "
            "to the volume and first page."
        ),
        "claims": (
            "`Brachyury functions primarily as a transcriptional activator genome-wide`; `an "
            "unexpected gene-regulatory feedback loop consisting of Brachyury, Foxa2, and Sox17 "
            "directs proper stem-cell lineage commitment during streak formation`."
        ),
        "does_not_claim": (
            "the abstract gives no dose, Hill coefficient, EC50 or half-maximal value: the words "
            "Hill, EC50, dose and half-maximal do not appear in it."
        ),
    },
    {
        "key": "lolas2014_fulltext",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3970479/",
        "retrieved": (
            "the full text in PMC, read for the two figure panels the program cites, for the "
            "direction of each link of the feedback loop, and for any number."
        ),
        "claims": (
            "`Brachyury first targets and promotes the expression of Foxa2 and Sox17`, and `high "
            "levels of Sox17 would, in turn, directly or indirectly repress the expression of "
            "Brachyury`. Fig. 3A is a western blot of `Foxa2 and Sox17 protein levels after "
            "Brachyury knockdown`. Fig. 3C shows that `Sox17 represses the expression of Brachyury "
            "during in vitro streak formation` after Sox17 overexpression. Both signs the program "
            "states for the Tbxt/Sox17 pair are supported."
        ),
        "does_not_claim": (
            "no numerical dose, strength, fold-change threshold, Hill coefficient, EC50 or "
            "half-maximal value is given for the Brachyury-Sox17 relationship anywhere. The Sox17 "
            "overexpression effect is reported as `significant down-regulation of Brachyury at EB "
            "day 4` with NO quantitative measure. It also does not claim a direct edge: the "
            "repression is `directly or indirectly`. And Fig. 1A, which the program cites for the "
            "ChIP-seq, shows `representative tracks for Brachyury-bound seq-regions` for Fgf8, "
            "Foxa2, Foxj1 and Dusp6 - Sox17 is not among the tracks displayed in that panel."
        ),
    },
    {
        "key": "saka_smith_2007",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC1885807/",
        "retrieved": (
            "the full text in PMC of Saka & Smith, `A mechanism for the sharp transition of "
            "morphogen gradient interpretation in Xenopus`, BMC Dev Biol 2007;7:47, read for every "
            "numerical parameter it states and for whether those values were measured or chosen. "
            "This is the nearest thing in the literature to a stated cooperativity for a TGF-beta "
            "input onto a Brachyury orthologue."
        ),
        "claims": (
            "the model's parameter values as stated: `ka = 5.5, kb = 5.4` (synthesis), `alpha = 6, "
            "beta = 3` (cooperativities of repression), `mu = 3` (cooperativity of induction), `kda "
            "= kdb = 1`; and an alternative set `(ka, kb) = (6, 14)` with `(kda, kdb) = (1, 5)` and "
            "`mu = 2`. The mutual repression modelled is between Xbra (the Xenopus Brachyury "
            "orthologue) and Gsc, with Xom mediating Xbra's repression of Gsc."
        ),
        "does_not_claim": (
            "none of these values is a measurement: they are chosen for the simulations, and the "
            "paper states simplifying assumptions such as `for the sake of simplicity we assume "
            "that M stays constant throughout the simulation`. The same cooperativity of induction "
            "is given as 3 in one parameter set and 2 in another, which is itself evidence that it "
            "is a choice and not an observation. Sox2, Sox17 and Nodal are NOT part of the model, "
            "the ligand is activin, and the repressed partner is Gsc, not SOX2."
        ),
    },
    {
        "key": "green_smith_dose_thresholds",
        "url": f"{_EPMC}?query=%28AUTH:%22Green+JB%22+AND+AUTH:%22Smith+JC%22%29+AND+"
        "%28PUB_YEAR:1990+OR+PUB_YEAR:1992%29&resultType=core&format=json",
        "retrieved": (
            "the Europe PMC core records and FULL abstracts of Green & Smith, Nature 1990;347:391-394 "
            "(PMID 1699129); Green, New & Smith, Cell 1992;71:731-739 (PMID 1423628); and Green et "
            "al., Development 1990;108:173-183 (PMID 2351061). Fetched because a measured ligand "
            "dose threshold, if one exists anywhere for this kind of interaction, is here."
        ),
        "claims": (
            "measured dose behaviour in Xenopus animal-pole tissue: `responding cells distinguish "
            "sharply between doses of pure XTC-MIF differing by less than 1.5-fold`; `Two different "
            "response thresholds have been found, defining three cell states`; `Each state is "
            "induced in a narrow dose range bounded by sharp thresholds`; and `a similar minimum "
            "active concentration (0.1-0.2 ng ml-1)` for XTC-MIF and XbFGF, with muscle induced by "
            "XbFGF only `above 24 ng ml-1`."
        ),
        "does_not_claim": (
            "none of the three abstracts names Xbra or Brachyury, and none reports a Hill "
            "coefficient or a half-maximal concentration for a named target gene. A minimum active "
            "concentration and a 1.5-fold discrimination are not the quantity the program's "
            "`threshold` field is - the regulator level at which the Hill term is one half - and "
            "turning the 1.5-fold into a Hill exponent would be a modelling step, not a reading. "
            "The ligand is activin/XTC-MIF in Xenopus, not Nodal, and the concentration is of "
            "ligand in the medium, not of the regulator protein the program integrates."
        ),
    },
    {
        "key": "damour2005",
        "url": f"{_EPMC}?query=TITLE:%22Efficient+differentiation+of+human+embryonic+stem+cells+to+"
        "definitive+endoderm%22&resultType=core&format=json",
        "retrieved": (
            "the Europe PMC core record and FULL abstract of D'Amour et al., Nat Biotechnol "
            "2005;23:1534-1541, PMID 16258519, fetched as the standard source for activin A dosing "
            "onto definitive endoderm."
        ),
        "claims": (
            "`Differentiation of hES cells in the presence of activin A and low serum produced "
            "cultures consisting of up to 80% definitive endoderm cells`."
        ),
        "does_not_claim": (
            "the abstract gives no activin A concentration, no dose threshold, no half-maximal dose, "
            "no fold-induction for SOX17 and no Hill coefficient; the only number in it is the 80%. "
            "SOX17 is not named in the abstract."
        ),
    },
    {
        "key": "tosic2019",
        "url": f"{_EPMC}?query=AUTH:%22Tosic+J%22+AND+PUB_YEAR:2019&resultType=core&format=json",
        "retrieved": (
            "the Europe PMC core record and FULL abstract of Tosic et al., `Eomes and Brachyury "
            "control pluripotency exit and germ-layer segregation by changing the chromatin state`, "
            "Nat Cell Biol 2019;21:1518-1531, PMID 31792383. Fetched because a search indicated it "
            "as the source for Brachyury repressing Sox2 - the rule the program attributes to "
            "Thomson 2011."
        ),
        "claims": (
            "`the induction of ME gene programs critically relies on the T-box transcription factors "
            "Eomesodermin ... and Brachyury, which concomitantly repress pluripotency and NE gene "
            "programs`; cells deficient in these factors `retain pluripotency and differentiate to "
            "NE lineages despite the presence of ME-inducing signals transforming growth factor beta "
            "(TGF-beta)/Nodal and Wnt`."
        ),
        "does_not_claim": (
            "the ABSTRACT does not name Sox2 - it says `pluripotency and NE gene programs` - so the "
            "specific Brachyury-Sox2 edge is taken from the curated Reactome entry below and not "
            "from this abstract. It reports no Hill coefficient, EC50, half-maximal value, "
            "fold-repression or dose."
        ),
    },
    {
        "key": "reactome_eomes_tbxt_sox2",
        "url": "https://reactome.org/content/detail/R-MMU-9756665",
        "retrieved": (
            "the Reactome curated event `Eomes and Tbxt bind the Sox2 gene` with its summation and "
            "its literature reference."
        ),
        "claims": (
            "verbatim: `Eomes and Tbxt positively regulate the exit from pluripotency by binding and "
            "repressing the expression of Sox2, a gene that maintains pluripotency (Tosic et al. "
            "2019)`, citing Tosic et al. 2019, Nat Cell Biol, PMID 31792383. This names both factors "
            "of the program's `rule Tbxt inhibits SOX2` and supports its sign."
        ),
        "does_not_claim": (
            "no binding affinity, fold-repression, Hill coefficient or concentration is given. A "
            "curated existence statement is not a value for any of the three fields."
        ),
    },
    {
        "key": "vincent2003",
        "url": f"{_EPMC}?query=TITLE:%22Cell+fate+decisions+within+the+mouse+organizer+are+governed+"
        "by+graded+Nodal+signals%22&resultType=core&format=json",
        "retrieved": (
            "the Europe PMC core record and FULL abstract of Vincent et al., Genes Dev "
            "2003;17:1646-1662, PMID 12842913, PMCID PMC196136, fetched as the mouse paper on GRADED "
            "Nodal, which is what a threshold on a Nodal rule would need."
        ),
        "claims": (
            "`graded Nodal/Smad2 signals govern allocation of the axial mesendoderm precursors that "
            "selectively give rise to the ADE and PCP mesoderm`, by conditional Smad2 inactivation "
            "in the epiblast and by a novel allele lacking the proximal epiblast enhancer."
        ),
        "does_not_claim": (
            "the abstract names neither Sox17 nor Brachyury, and gives no numerical dose, threshold, "
            "fold-change, EC50 or Hill coefficient. `Graded` here is established by allele series "
            "and tissue outcome, not by a measured concentration scale."
        ),
    },
    {
        "key": "sox17_nodal_open_access",
        "url": f"{_EPMC}?query=%28ABSTRACT:%22SOX17%22+AND+ABSTRACT:%22Nodal%22%29+AND+OPEN_ACCESS:y"
        "&resultType=core&format=json",
        "retrieved": (
            "the Europe PMC core records and FULL abstracts of four open-access papers naming both "
            "SOX17 and Nodal, fetched to find any source at all that names both factors of `rule "
            "Nodal activates SOX17`, which the program's own citation does not."
        ),
        "claims": (
            "iScience 2023, PMID 37502260, PMCID PMC10368912: `Activin/Nodal signaling is the "
            "primary driver of definitive endoderm formation`, and Activin/Nodal-SMAD2 `participates "
            "in the initiation of mesendoderm and DE specification` for `SOX17+ DE`. Cell & "
            "Bioscience 2022, PMID 36333732, PMCID PMC9636699: `TET proteins stimulated SOX17 "
            "through the NODAL signaling pathway`. Both name both factors and support the sign the "
            "rule states, in human pluripotent cells in vitro."
        ),
        "does_not_claim": (
            "neither gives any numerical dose, concentration, threshold, EC50, fold-change or Hill "
            "coefficient. And the relation is context-dependent rather than universal: Oncotarget "
            "2016, PMID 27283990, PMCID PMC5216926, reports the opposite association in a germ-cell "
            "tumour, where `induction of NODAL signaling` accompanies `downregulation of seminoma "
            "markers, like SOX17` and strong induction of SOX2."
        ),
    },
)

#: The free web searches that were run, kept apart from the fetches. A search snippet is never used
#: as evidence for a class: it only pointed at a source, which was then fetched and read.
SEARCHES_RUN = (
    "Nodal Activin dose-response Hill coefficient Brachyury T induction measured EC50",
    "SOX17 endoderm Activin dose threshold half-maximal induction human embryonic stem cells quantitative",
    "published mathematical model germ layer tristable switch SOX17 Brachyury SOX2 Hill coefficient "
    "parameter values gastrulation",
    '"SOX17" Activin Nodal dose-dependent induction definitive endoderm open access eLife '
    "quantitative Smad2 high dose",
    '"Brachyury" represses "Sox2" OR "SOX17 represses SOX2" direct repression evidence',
)

# --- the findings --------------------------------------------------------------------------------

_NO_VALUE_ANYWHERE = (
    "no source fetched in this lane states a value for this field, by measurement or by fitting"
)

#: Per rule: the class and the record that is the same for all three of its fields.
_PER_RULE: dict[str, dict[str, Any]] = {
    "Nodal activates TBXT": {
        "class": "E_existence_only",
        "fetched": ("conlon1994_epmc", "conlon1994_publisher", "tosic2019", "saka_smith_2007"),
        "source_claims": (
            "Conlon et al. 1994, the program's own citation, establishes that nodal is required for "
            "the primitive streak and that a quarter of mutant embryos still form patches of "
            "posterior mesodermal character. Tosic et al. 2019 names both TGF-beta/Nodal and "
            "Brachyury and places Nodal among the ME-inducing signals upstream of T-box factor "
            "induction, which is the existence of the rule's relation."
        ),
        "attribution": "names_one_factor_only",
        "attribution_note": (
            "Conlon's abstract names nodal and does not name Brachyury, Bra or the T gene. The full "
            "text is paywalled and was not read, so this is a statement about the abstract, twice "
            "retrieved. What the abstract supports is nodal's requirement for the primitive streak "
            "and for posterior mesoderm - one step removed from the TBXT gene the rule names."
        ),
    },
    "Nodal activates SOX17": {
        "class": "E_existence_only",
        "fetched": (
            "kanaiazuma2002_epmc",
            "kanaiazuma2002_publisher",
            "sox17_nodal_open_access",
            "vincent2003",
        ),
        "source_claims": (
            "The rule's relation exists and is supported, but by sources the program does not cite: "
            "iScience 2023 (`Activin/Nodal signaling is the primary driver of definitive endoderm "
            "formation`, for SOX17+ DE) and Cell & Bioscience 2022 (`TET proteins stimulated SOX17 "
            "through the NODAL signaling pathway`). Vincent et al. 2003 establishes that Nodal/Smad2 "
            "signalling is GRADED in the mouse organizer, which is the premise a threshold needs."
        ),
        "attribution": "names_one_factor_only",
        "attribution_note": (
            "Kanai-Azuma et al. 2002 does not contain the word Nodal in its abstract - checked by "
            "two independent fetches, Europe PMC and the publisher page, the second asked the "
            "question explicitly. It is a Sox17-null mouse study: evidence that Sox17 is REQUIRED "
            "for gut endoderm, which is about what Sox17 does, not about what induces it. This is "
            "the same shape of weakness the predecessor lane was corrected for."
        ),
    },
    "Tbxt inhibits SOX2": {
        "class": "E_existence_only",
        "fetched": ("thomson2011_epmc", "thomson2011_fulltext", "tosic2019", "reactome_eomes_tbxt_sox2"),
        "source_claims": (
            "The rule's relation exists, and the Reactome curated event R-MMU-9756665 states it "
            "naming both factors: `Eomes and Tbxt positively regulate the exit from pluripotency by "
            "binding and repressing the expression of Sox2`, citing Tosic et al. 2019. The program "
            "does not cite either."
        ),
        "attribution": "names_both_but_not_this_direction",
        "attribution_note": (
            "Thomson et al. 2011, the program's own citation, names both factors in its full text "
            "and states the OPPOSITE direction: `Sox2 specifically represses Brachyury and the "
            "mesendodermal lineage`. The factor it names as the repressor of the neural ectoderm "
            "program is OCT4 - `Oct4 specifically represses Sox1 and the NE lineage` - and Oct4 is "
            "not in this program at all. So this rule's citation supports the rule BELOW it (line "
            "22) and not this one."
        ),
    },
    "Sox2 inhibits TBXT": {
        "class": "E_existence_only",
        "fetched": ("thomson2011_epmc", "thomson2011_fulltext"),
        "source_claims": (
            "Thomson et al. 2011's full text states this rule's relation in the program's own "
            "direction and names both factors: `Sox2 specifically represses Brachyury and the "
            "mesendodermal lineage`, and `Sox2 specifically represses only the mesendodermal fate`. "
            "Of the program's seven rules this is the one whose own citation squarely supports it."
        ),
        "attribution": "supports_declaration",
        "attribution_note": (
            "the citation on this rule is the bare `Thomson et al. 2011`, without the journal and "
            "pages the rule above carries; the fuller citation on line 21 is the one that does NOT "
            "support its own rule. Both point at the same work."
        ),
    },
    "Sox17 inhibits TBXT": {
        "class": "E_existence_only",
        "fetched": ("lolas2014_epmc", "lolas2014_fulltext"),
        "source_claims": (
            "Lolas et al. 2014 states this rule's relation and the cited panel shows it: Fig. 3C "
            "demonstrates that `Sox17 represses the expression of Brachyury during in vitro streak "
            "formation` after Sox17 overexpression, and the paper's model says `high levels of Sox17 "
            "would, in turn, directly or indirectly repress the expression of Brachyury`."
        ),
        "attribution": "supports_declaration",
        "attribution_note": (
            "with one qualification the program does not carry: the source says `directly or "
            "indirectly`, while a rule in this language is a direct regulatory edge. The effect is "
            "reported as `significant down-regulation of Brachyury at EB day 4` with no number."
        ),
    },
    "Sox17 inhibits SOX2": {
        "class": "U_unsourced",
        "fetched": ("thomson2011_fulltext", "sox17_nodal_open_access"),
        "source_claims": (
            'nothing: the declaration cites no work - its evidence is `inferred "endoderm excludes '
            "ectoderm program\"`, which BioLang's own grammar defines as not measured - and no "
            "source fetched in this lane supports SOX17 repressing SOX2 in either a measured or a "
            "modelled form."
        ),
        "attribution": "no_source_cited",
        "attribution_note": (
            "the one fetched source that names both factors, Oncotarget 2016, runs the other way: "
            "SOX2-deficient cells MAINTAIN SOX17 expression, which bears on SOX2 repressing SOX17, "
            "not on this rule's direction, and it is a germ-cell tumour model rather than "
            "gastrulation. A search for `SOX17 represses SOX2` returned nothing that states it. "
            "This rule is unsourced on every axis, and all three of its numbers with it."
        ),
    },
    "Tbxt activates SOX17": {
        "class": "E_existence_only",
        "fetched": ("lolas2014_epmc", "lolas2014_fulltext"),
        "source_claims": (
            "Lolas et al. 2014 states this rule's relation: `Brachyury first targets and promotes "
            "the expression of Foxa2 and Sox17`, Fig. 3A is a western blot of `Foxa2 and Sox17 "
            "protein levels after Brachyury knockdown`, and the abstract reports that `Brachyury "
            "functions primarily as a transcriptional activator genome-wide`. The sign correction at "
            "5cbce26 is what this source supports, and it stands."
        ),
        "attribution": "supports_declaration",
        "attribution_note": (
            "one pointer in the evidence string is inexact: the rule cites `Fig. 1A ChIP-seq`, and "
            "Fig. 1A shows representative Brachyury-bound tracks for Fgf8, Foxa2, Foxj1 and Dusp6 - "
            "Sox17 is not among the tracks displayed in that panel. The Sox17 evidence in the paper "
            "is Fig. 3A and the model text, both of which do support the rule."
        ),
    },
}

#: Per field: the commensurability verdict and what would have had to exist for an M.
_PER_FIELD: dict[str, dict[str, Any]] = {
    "strength": {
        "commensurable": "dimensionless_and_expressible",
        "what_an_M_would_need": (
            "a measured saturating fold-change of the TARGET gene's transcription rate caused by "
            "this one regulator alone, in the cell being modelled, expressed as the fraction of the "
            "target's maximum rate the regulator removes (or drives). Every one of this program's "
            "seven rules writes a strength, six of them 1.0, and 1.0 on an inhibitor is complete "
            "shutoff: R = prod(1 - s H) reaches exactly 0 at saturation, which is a switch and not a "
            "brake. Nothing fetched reports such a fold-change for any of the seven."
        ),
    },
    "threshold": {
        "commensurable": "arbitrary_units_no_conversion",
        "what_an_M_would_need": (
            "a measured regulator concentration at half-maximal target response, in units this "
            "program could convert to its own. It could not: the run is in arbitrary units by the "
            "runtime's own statement, this program's `translation_rate` and `mrna_half_life` are "
            '`evidence: inferred "dimensionless"`, and the NODAL level the two Nodal rules are '
            "read against is set by `run_gastrulation`'s unsourced `nodal_max = 6.0`. Writing the "
            "measured unit onto the field instead is refused: the unlocated GRN runtime raises on "
            "any rule with a `threshold_unit`. The nearest measured quantities found - a minimum "
            "active concentration of 0.1-0.2 ng/ml and a 1.5-fold dose discrimination in Xenopus "
            "animal caps - are ligand concentrations in the medium for fate outcomes, not the "
            "regulator level at half-maximal induction of a named target."
        ),
    },
    "hill": {
        "commensurable": "dimensionless_and_expressible",
        "what_an_M_would_need": (
            "a directly observed exponent, which is not a thing the literature reports: a Hill "
            "coefficient is obtained by fitting a Hill function to a dose-response curve, or chosen "
            "inside a model. This was declared at registration, before the search. Nothing fetched "
            "reports a Hill coefficient for any of this program's seven rules by either route; the "
            "one cooperativity found for a corresponding interaction is a modeller's choice."
        ),
    },
}

#: The one slot where a fetched source does state a value for the quantity - by choosing it inside a
#: model, in another species, for another ligand. It is named as fitted and is NOT counted as
#: measured. The strict and the loose measured counts are both unaffected by it.
_MODELLED_SLOTS: dict[tuple[str, str], dict[str, Any]] = {
    ("Nodal activates TBXT", "hill"): {
        "class": "F_fitted_or_modelled",
        "value_in_the_source": "mu = 3, and mu = 2 in a second parameter set of the same paper",
        "fitted_how": (
            "chosen for the simulations, not fitted to data and not measured: Saka & Smith 2007 "
            "report `mu = 3` as the cooperativity of induction of Xbra and Gsc by the morphogen in "
            "one parameter set and `mu = 2` in another, alongside the assumption that the morphogen "
            "`stays constant throughout the simulation`. One paper giving the same constant two "
            "values is itself the evidence that it is a choice."
        ),
        "why_not_measured": (
            "it is a modeller's choice, so by the rule registered before the search it may not be "
            "labelled measured. Three further mismatches are recorded rather than smoothed over: "
            "the ligand is activin, not Nodal; the organism is Xenopus; and the target is Xbra, the "
            "Brachyury orthologue, in a model whose other half is Gsc, which this program does not "
            "contain. It is reported because it is the closest the literature came to this field, "
            "not because it fills it."
        ),
    }
}

#: Declined on purpose: Saka & Smith also state cooperativities of REPRESSION, `alpha = 6` and `beta
#: = 3`. They are for Gsc repressing Xbra through Xom. None of this program's four inhibiting rules
#: is that pair, so mapping those exponents onto `Tbxt inhibits SOX2`, `Sox2 inhibits TBXT`, `Sox17
#: inhibits TBXT` or `Sox17 inhibits SOX2` would be stretching a source to fill a cell. They are
#: recorded here and used for nothing.
DECLINED_MAPPINGS = (
    {
        "source": "saka_smith_2007",
        "values": "alpha = 6 and beta = 3, cooperativities of repression",
        "declined_for": (
            "all four inhibiting rules of this program, and for the hill field of each of them",
        ),
        "why": (
            "the modelled repression is Gsc on Xbra by way of Xom. Neither Gsc nor Xom is in this "
            "program, and none of its four repressing pairs is that pair. A value for a different "
            "pair is not a value for these, even in a model."
        ),
    },
    {
        "source": "green_smith_dose_thresholds",
        "values": "a minimum active concentration of 0.1-0.2 ng/ml, and discrimination of doses "
        "differing by less than 1.5-fold",
        "declined_for": ("the threshold field of all seven rules",),
        "why": (
            "these are measured, and they are measurements of something else: a ligand "
            "concentration in the medium at which a fate outcome changes, in Xenopus, for "
            "activin/XTC-MIF. The field is the level of the regulator the program integrates, at "
            "half-maximal induction of the named target, in the program's own arbitrary units. "
            "Converting the one into the other needs a scale this program does not have, and "
            "turning the 1.5-fold into a Hill exponent would be a fit, not a reading."
        ),
    },
)


def _finding(rule: str, field: str) -> dict[str, Any]:
    per_rule = _PER_RULE[rule]
    per_field = _PER_FIELD[field]
    modelled = _MODELLED_SLOTS.get((rule, field))
    urls = [s["url"] for s in SOURCES_FETCHED if s["key"] in per_rule["fetched"]]
    fetched_what = " | ".join(
        f"{s['key']}: {s['retrieved']}" for s in SOURCES_FETCHED if s["key"] in per_rule["fetched"]
    )
    does_not_claim = " | ".join(
        f"{s['key']}: {s['does_not_claim']}" for s in SOURCES_FETCHED if s["key"] in per_rule["fetched"]
    )
    record: dict[str, Any] = {
        "class": modelled["class"] if modelled else per_rule["class"],
        "fetched_url": urls,
        "fetched_what": fetched_what,
        "source_claims": per_rule["source_claims"],
        "source_does_not_claim": (
            f"for this field specifically: {_NO_VALUE_ANYWHERE}. Per source: {does_not_claim}"
            if not modelled
            else f"per source: {does_not_claim}"
        ),
        "cited_source_attribution": per_rule["attribution"],
        "attribution_note": per_rule["attribution_note"],
        "commensurable": per_field["commensurable"],
        "what_an_M_would_need": per_field["what_an_M_would_need"],
        "searches": list(SEARCHES_RUN),
        "measured_elsewhere": False,
        "measured_elsewhere_why_not": (
            "no fetched source reports a measurement of this quantity for this interaction in any "
            "organism or cell type, so the loose count does not differ from the strict one here."
        ),
    }
    if modelled:
        record.update(
            {
                "value_in_the_source": modelled["value_in_the_source"],
                "fitted_how": modelled["fitted_how"],
                "why_not_measured": modelled["why_not_measured"],
                "fetched_url": urls + ["https://pmc.ncbi.nlm.nih.gov/articles/PMC1885807/"],
            }
        )
    return record


#: What the program's own citations do and do not establish, one entry per citation.
ATTRIBUTION_FINDINGS: tuple[dict[str, Any], ...] = tuple(
    {
        "rule": rule,
        "cited_work": _PER_RULE[rule]["fetched"][0],
        "attribution": _PER_RULE[rule]["attribution"],
        "note": _PER_RULE[rule]["attribution_note"],
    }
    for rule in _PER_RULE
)

#: What the code or the literature contradicts in what was assumed going in. Recorded because being
#: corrected is worth more than being agreed with.
CONTRADICTIONS: tuple[dict[str, Any], ...] = (
    {
        "about": "data/demo/gastrulation.bio line 3, `Mutual repression makes the choice sharp`",
        "sharper_than_reported": (
            "the line is false of the three-fate choice but NOT of every pair: SOX2 and TBXT DO "
            "repress each other (lines 21 and 22). What is false is the rest - SOX17 is repressed by "
            "nothing, SOX17 to SOX2 is one-way, and the SOX17/TBXT pair is repression one way and "
            "ACTIVATION the other, which is a negative feedback loop and not mutual repression. So "
            "one of the three pairs is mutually repressing and the sharpness of the three-way choice "
            "is not what the line says it is."
        ),
        "edited": False,
    },
    {
        "about": "the program's citation for `rule Tbxt inhibits SOX2` (line 21)",
        "finding": (
            "Thomson et al. 2011 states the opposite direction in its full text - `Sox2 specifically "
            "represses Brachyury and the mesendodermal lineage` - and names OCT4, not Brachyury, as "
            "the repressor of the neural ectoderm program. A source that does support the rule "
            "exists and the program does not cite it: Tosic et al. 2019, as curated at Reactome "
            "R-MMU-9756665. Reported, not edited: the numbers and the evidence strings are "
            "untouched."
        ),
    },
    {
        "about": "the corpus-wide reading that `lateral_inhibition.bio` is the one place with a "
        "traceable strength, threshold and hill",
        "finding": (
            "true as stated, and it does not make that number measured. Its source is Collier et al. "
            "1996's own dimensionless Hill functions, `f(d) = d^2 / (0.01 + d^2)` and `g(n) = 1 / (1 "
            "+ 100 n^2)`, which is a parameterised theoretical model. So the count of MEASURED rule "
            "numbers in the hand-authored corpus is 0, not 1 program's worth."
        ),
    },
    {
        "about": "whether the search could have succeeded for `threshold` at all",
        "finding": (
            "it could not have, as the program stands, and that is a property of the code rather "
            "than of the literature: the unlocated GRN runtime refuses any rule carrying a "
            "`threshold_unit`, so a measured EC50 could not be written into these seven fields even "
            "if one existed for the right regulator in the right cell. A test in "
            "tests/test_rule_number_sources.py establishes the refusal."
        ),
    },
    {
        "about": "the premise that the module's own numbers make any added inhibitor a switch",
        "finding": (
            "confirmed from the runtime rather than assumed: R = prod(1 - s H) with s = 1.0 reaches "
            "exactly 0 at saturating regulator level, and all four of this program's inhibiting "
            "rules carry `strength: 1.0`. A brake would need a strength below 1, and no source "
            "fetched here supplies one for any of them."
        ),
    },
)

#: The reading, in the words the result carries.
VERDICT = (
    "Of the 21 rule numbers of data/demo/gastrulation.bio, 0 have a measured value with a quoted "
    "source - 0 strictly, and 0 counting a measurement of the same quantity in any organism or cell "
    "type. 20 of the 21 are E_existence_only or U_unsourced: for 17 the interaction the rule states "
    "is supported by a source that was fetched and read, while the number is not; for 3, all of "
    "them on `rule Sox17 inhibits SOX2`, there is no source for the interaction either. 1 slot, the "
    "`hill` of `rule Nodal activates TBXT`, is F_fitted_or_modelled: Saka & Smith 2007 state a "
    "cooperativity of induction of mu = 3 for activin onto Xbra, and mu = 2 in a second parameter "
    "set of the same paper, chosen for a simulation in another species for another ligand. It is "
    "named as modelled and counted as measured nowhere. No number in the program was changed."
)


#: The filled findings table, one entry per rule number, keyed as the registration keys the slots.
FINDINGS: dict[tuple[str, str], dict[str, Any]] = {
    (s["rule"], s["field"]): _finding(s["rule"], s["field"]) for s in slots()
}


def findings() -> dict[tuple[str, str], dict[str, Any]]:
    """The filled table, which the census writer loads into the registered-empty one."""
    return dict(FINDINGS)


def attribution_findings() -> tuple[dict[str, Any], ...]:
    """What the program's own citation does for each rule it sits on, one entry per rule."""
    return ATTRIBUTION_FINDINGS


def counts() -> dict[str, int]:
    """The classes recounted from the findings, never asserted."""
    out = {c: 0 for c in CLASSES}
    for record in FINDINGS.values():
        out[record["class"]] += 1
    return out
