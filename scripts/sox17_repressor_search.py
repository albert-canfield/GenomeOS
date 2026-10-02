#!/usr/bin/env python3
"""Set-one-in over a closed list of SOURCED repressors of SOX17, each added alone.

lane-mesoderm's diagnosis (d340e5c) is the input: SOX17's repression factor in the runtime is
exactly 1.00000 because nothing in data/demo/gastrulation.bio represses SOX17, while TBXT carries
two repressors, one of them Sox17. This lane asks the one question that follows: does the
literature name a factor that represses SOX17, or antagonises the endoderm program, IN THE
MESODERM LINEAGE -- and if a sourced one is set into the program alone, at the program's own
default numbers, does a mesoderm band come back?

The reading is one of exactly two: "a sourced repressor of SOX17 restores a mesoderm band", or
"no sourced repressor on the list does". If the literature yielded no citable candidate at all,
that would itself be the result. It yielded two.

The difference from lane-mesoderm: that lane had no network, so its citations were bibliography
recalled from its own knowledge. EVERY citation below was FETCHED THIS SESSION and read, and each
entry records the URL, what was retrieved (abstract or full text) and what the source does and
does not claim. A factor shown to repress SOX17 in a different lineage or at a different stage is
EXCLUDED and listed as excluded, not stretched into a candidate.

Nothing is tuned. The bound, the tail rule, the run and the number convention are lane-mesoderm's,
imported from its committed module rather than restated, so that they are the same objects and not
merely the same words. No strength, threshold, Hill coefficient, basal, max, nodal_max or
decay_length is varied in any run, and the sign correction at 5cbce26 stands in every run.
"""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

from genomeos import manifest as mf
from genomeos.results import save_result
from genomeos.runtime.gastrulation import CENSUS_MODEL_RUN

# The bound, the tail rule, the recording interval, the number convention and the variant builder
# are lane-mesoderm's committed ones, imported so that "the same weak bound and the same tail rule"
# is an identity and not a transcription.
from scripts.mesoderm_diagnosis import (
    BOUND_LAST_FRACTION,
    BOUND_SHARE,
    MODULE,
    RECORD_EVERY,
    evaluate,
    variant_text,
)
from scripts.mesoderm_diagnosis import REGISTRATION as PRIOR_REGISTRATION

ENTRY = "scripts/sox17_repressor_search.py"
OWN_CODE = (ENTRY, "tests/test_sox17_repressor_search.py")

#: Europe PMC's REST search and full-text endpoints, the one route used for every fetch below.
FETCH_ROUTE = "https://www.ebi.ac.uk/europepmc/webservices/rest/"

#: THE REGISTRATION, in committed code, written and committed BEFORE any run of this file.
REGISTRATION = {
    "registered": "2026-10-02",
    "lane": "lane-sox17rep",
    "builds_on": {
        "registration": "0790d51",
        "result": "d340e5c",
        "row": "65d283d",
        "figure_correction": "99b86de",
        "what_it_established": "the gastrulation model has no mesoderm at all (0.000 at all 41 "
        "recorded samples) and SOX17's repression factor is exactly 1.00000 "
        "because nothing in the program represses SOX17, while TBXT carries "
        "two repressors, one of them Sox17 itself",
        "readings_carried_unchanged": [
            "M1, a MIXL1 node, CLEARS THE WEAK BOUND: sustained 0.0167 against a bound of "
            "0.01. It does not restore mesoderm.",
            "R1 recovers a band by LEAVING OUT a SOURCED rule (Sox17 inhibits TBXT, Lolas "
            "2014 Fig. 3C), so it LOCATES WHERE THE BAND RESTS. It is not a fix, and the "
            "rule stays.",
            "The finding is STRUCTURAL: the program lacks any repressor of SOX17. It names "
            "no biological factor.",
        ],
    },
    "question": "Does the literature name a factor that represses SOX17, or antagonises the "
    "endoderm program, in the mesoderm lineage; and does setting a sourced one into "
    "the program alone, at the program's own default numbers, clear the weak mesoderm "
    "bound? The answer is one of exactly two readings: 'a sourced repressor of SOX17 "
    "restores a mesoderm band', or 'no sourced repressor on the list does'.",
    "what_this_is_not": list(PRIOR_REGISTRATION["what_this_is_not"])
    + [
        "Not a naming of the missing factor. Whichever entry clears the bound, the reading "
        "stays 'the program lacked a repressor of SOX17' and the structural finding names no "
        "biological factor. A candidate clearing a deliberately weak bound in a four-gene "
        "module is not evidence that the factor specifies mesoderm in an embryo.",
        "Not a citation of anything unfetched. Every entry's source was retrieved this "
        "session through "
        + FETCH_ROUTE
        + " and is quoted below at the level actually read; what each source does NOT claim is "
        "recorded beside what it does.",
    ],
    "candidate_rules": {
        "1": "a factor the literature shows represses SOX17, or antagonises the endoderm "
        "program, IN THE MESODERM LINEAGE, at or near gastrulation",
        "2": "the citation must have been FETCHED AND READ this session; a candidate without a "
        "fetched citation is EXCLUDED and listed as excluded",
        "3": "a source showing repression in a different lineage or at a different stage is "
        "EXCLUDED with that reason stated, not stretched",
        "4": "each candidate is added ALONE, with the source-faithful sign and the program's "
        "default strengths; nothing is fitted or varied",
        "5": "the list is closed at the count below and may not grow after any run",
    },
    # the bound, the run and the number convention are the prior lane's objects, not copies
    "bound": PRIOR_REGISTRATION["bound"],
    "run": PRIOR_REGISTRATION["run"],
    "number_convention": PRIOR_REGISTRATION["number_convention"],
    "citation_verification": {
        "network": "USED. This lane is free -- 0 model requests, no money, no paid API -- and the "
        "web is free. Every citation below was fetched through " + FETCH_ROUTE + ".",
        "difference_from_the_prior_lane": "lane-mesoderm had no network, so every one of its "
        "citations was bibliography given from its own knowledge "
        "and was never fetched. Its seven are re-checked in "
        "prior_lane_verification below, against fetched records.",
        "level_read": "for each entry: 'abstract' means the Europe PMC abstract was retrieved and "
        "read; 'full text' means the Europe PMC fullTextXML was retrieved and the "
        "quoted sentence read in it.",
    },
    "closed_list": [
        {
            "id": "C1",
            "name": "a MESP1 node repressing SOX17, with the module's NODAL drive on it",
            "kind": "missing species and rules",
            "factor": "MESP1",
            "citation": "Costello I, Pimeisl IM, Draeger S, Bikoff EK, Robertson EJ, Arnold SJ "
            "(2011) The T-box transcription factor Eomesodermin acts upstream of Mesp1 "
            "to specify cardiac mesoderm during mouse gastrulation. Nat Cell Biol "
            "13:1084-1091. PMID 21822279, PMCID PMC4531310, doi "
            "10.1038/ncb2304.",
            "fetched": {
                "level": "full text",
                "url": FETCH_ROUTE + "PMC4531310/fullTextXML",
                "also": FETCH_ROUTE
                + "search?query=TITLE:%22Mesp1%20acts%20as%20a%20master%20regulator...%22 "
                "(Bondue et al. 2008, abstract only)",
            },
            "source_claims": "quoted verbatim from the fetched full text: 'Within this early "
            "sub-set Mesp1 directly represses genes required for formation of DE "
            "including Foxa2, Gsc and Sox17.' The same full text reports, for the "
            "drive: 'we manipulated ActivinA concentrations in differentiating ES "
            "cells and confirmed that low levels (5ng/ml) are sufficient for robust "
            "Mesp1 expression, while conversely maintaining cultures in high "
            "ActivinA concentrations (50ng/ml) leads to induction of Sox17.' Mesp1 "
            "marks the earliest cardiac mesoderm within the primitive streak at "
            "mouse gastrulation, which is the lineage and the stage rule 1 asks for.",
            "source_does_not_claim": "no dose, strength, threshold or Hill coefficient, so every "
            "number used here is the module's own. The sentence names "
            "Sox17 among the directly repressed DE genes; the "
            "panel-level ChIP evidence for the Sox17 locus is in Bondue "
            "et al. 2008 (Cell Stem Cell 3:69-84, PMID 18593560), whose "
            "ABSTRACT WAS FETCHED and says only 'Mesp1 also directly "
            "represses the expression of key genes regulating other "
            "early mesoderm and endoderm cell fates' WITHOUT NAMING "
            "Sox17; its full text is paywalled and was NOT read, so the "
            "E-box panel is recorded here as unverified. Activin is a "
            "NODAL-family ligand read by the same SMAD2/3 transducer, "
            "not NODAL itself; the Nodal->MESP1 rule is faithful to the "
            "source's sign, not to its ligand identity.",
            "species": [
                {
                    "gene": "MESP1",
                    "protein": "Mesp1",
                    "note": "sourced SOX17 repressor C1/C2; level unsourced, module convention",
                }
            ],
            "rules": [
                {
                    "src": "Nodal",
                    "dst": "MESP1",
                    "action": "activates",
                    "cite": "Costello et al. 2011, Nat Cell Biol 13:1084 (ActivinA induces Mesp1)",
                },
                {
                    "src": "Mesp1",
                    "dst": "SOX17",
                    "action": "inhibits",
                    "cite": "Costello et al. 2011, Nat Cell Biol 13:1084 "
                    "(Mesp1 directly represses Foxa2, Gsc and Sox17)",
                },
            ],
            "predicted": "the first brake SOX17 has ever had in this program; graded, because its "
            "carrier is driven by the clamped NODAL through the module's own "
            "activating numbers",
            "runnable": True,
        },
        {
            "id": "C2",
            "name": "a MESP1 node repressing SOX17, unregulated, in the module's W1/F1 shape",
            "kind": "missing species and rule",
            "factor": "MESP1",
            "citation": "Costello I, Pimeisl IM, Draeger S, Bikoff EK, Robertson EJ, Arnold SJ "
            "(2011) Nat Cell Biol 13:1084-1091. PMID 21822279, PMCID PMC4531310.",
            "fetched": {
                "level": "full text",
                "url": FETCH_ROUTE + "PMC4531310/fullTextXML",
            },
            "source_claims": "the same fetched sentence as C1 for the repression of Sox17.",
            "source_does_not_claim": "nothing about the carrier's level. C2 adds NO activator, so "
            "the runtime drives MESP1 at full max_rate, exactly the "
            "shape lane-mesoderm registered for W1 and F1 ('alone the "
            "added gene has no activator, so the runtime drives it at "
            "full max_rate and the drive is not gradient-shaped'). C1 "
            "and C2 are registered as two entries because the choice of "
            "carrier drive is a choice the source does not make, and "
            "picking one and reporting it would be picking a number.",
            "species": [
                {
                    "gene": "MESP1",
                    "protein": "Mesp1",
                    "note": "sourced SOX17 repressor C1/C2; level unsourced, module convention",
                }
            ],
            "rules": [
                {
                    "src": "Mesp1",
                    "dst": "SOX17",
                    "action": "inhibits",
                    "cite": "Costello et al. 2011, Nat Cell Biol 13:1084",
                },
            ],
            "predicted": "a brake on SOX17 everywhere along the axis, not only where NODAL is high",
            "runnable": True,
        },
        {
            "id": "C3",
            "name": "a HAND1 node repressing SOX17, unregulated, in the module's W1/F1 shape",
            "kind": "missing species and rule",
            "factor": "HAND1",
            "citation": "Lynch AT, Phillips N, Douglas M, Dorgnach M, Lin IH, Adamson AD, Darieva "
            "Z, Whittle J, Hanley NA, Bobola N, Birket MJ (2025) HAND1 level controls "
            "the specification of multipotent cardiac and extraembryonic progenitors "
            "from human pluripotent stem cells. EMBO J 44:2541-2565. PMID 40164946, "
            "PMCID PMC12048643.",
            "fetched": {
                "level": "full text",
                "url": FETCH_ROUTE + "PMC12048643/fullTextXML",
            },
            "source_claims": "quoted verbatim from the fetched full text: 'Conversely, genes "
            "potentially repressed by HAND1 were associated with processes including "
            "neural development, WNT signalling and endoderm development (including "
            "SOX17 and NODAL).' And: 'Regions of decreased accessibility were also "
            "enriched in motifs for GSC and SOX17, consistent with the "
            "downregulation of these genes', and 'SOX17 overlapped with the HAND1[-] "
            "regulon.' The system is human pluripotent stem cells differentiating to "
            "cardiac and extraembryonic mesoderm, which is the mesoderm lineage rule "
            "1 asks for.",
            "source_does_not_claim": "the word in the source is 'potentially repressed'. The "
            "evidence is expression and chromatin accessibility with "
            "motif enrichment and a regulon overlap; NO ChIP-seq "
            "binding of HAND1 at the SOX17 locus is shown. So C3's sign "
            "is sourced at the level of a HAND1-dependent "
            "downregulation, NOT of a demonstrated direct repression, "
            "and that is weaker than C1's. The same sentence names "
            "NODAL among the repressed genes, which this program cannot "
            "represent at all -- see excluded entry Y1.",
            "species": [
                {
                    "gene": "HAND1",
                    "protein": "Hand1",
                    "note": "sourced SOX17 repressor C3; level unsourced, module convention",
                }
            ],
            "rules": [
                {
                    "src": "Hand1",
                    "dst": "SOX17",
                    "action": "inhibits",
                    "cite": "Lynch et al. 2025, EMBO J 44:2541 (genes potentially repressed by "
                    "HAND1 include SOX17)",
                },
            ],
            "predicted": "the same shape as C2 with a different carrier: if C2 and C3 agree to the "
            "digit, the run is reporting the module's inhibiting numbers and not the "
            "identity of the factor, and that is worth knowing",
            "runnable": True,
        },
    ],
    "closed_list_n": 3,
    "no_pairwise": "This lane is set-one-in only, as briefed: each candidate added ALONE. No pair "
    "is run and no pair is registered.",
    "excluded": [
        {
            "id": "Y1",
            "item": "HAND1 (or CER1/LEFTY1, or SOX17 itself) acting through NODAL",
            "citation": "Lynch et al. 2025, EMBO J 44:2541-2565 (fetched, full text) names NODAL "
            "among the genes potentially repressed by HAND1; Perea-Gomez et al. 2002, "
            "Dev Cell 3:745-756 (fetched, abstract) names Cerberus-like and Lefty1 as "
            "the NODAL antagonists.",
            "reason": "EXCLUDED for the structural reason lane-mesoderm registered for A1 and "
            "carried unchanged: NODAL is not a modelled species in this program. Its "
            "gene declares max 0 and the runtime holds Nodal at a clamped external "
            "level in every Runge-Kutta stage, so no rule the module can state will "
            "move it. Any route to SOX17 that runs through NODAL is unrunnable here "
            "without changing the program's input semantics.",
            "excluded_from": "the closed candidate list; it cannot clear the bound",
        },
        {
            "id": "Y2",
            "item": "FOXP1 repressing Sox17",
            "citation": "Zhang Y, Li S, Yuan L, Tian Y, Weidenfeld J, Yang J, Liu F, Chokas AL, "
            "Morrisey EE (2010) Foxp1 coordinates cardiomyocyte proliferation through "
            "both cell-autonomous and nonautonomous mechanisms. Genes Dev 24:1746-1757. "
            "PMID 20713518, PMCID PMC2922503. Fetched: Europe PMC record, found by a "
            "full-text phrase search for 'represses Sox17' and for 'repression of "
            "Sox17', both of which return it.",
            "reason": "EXCLUDED by rule 3, wrong lineage and wrong stage. Foxp1 targets Sox17 in "
            "the ENDOCARDIUM, to control myocardial proliferation in an already formed "
            "heart. That is a mesoderm DERIVATIVE at organogenesis, not the mesoderm "
            "lineage at gastrulation, and the model has no heart. Named here rather "
            "than stretched into a candidate.",
            "excluded_from": "the closed candidate list",
        },
        {
            "id": "Y3",
            "item": "KLF5 repressing Sox17",
            "citation": "Lin SC, Wani MA, Whitsett JA, Wells JM (2010) Klf5 regulates lineage "
            "formation in the pre-implantation mouse embryo. Development 137:3953-3963. "
            "PMID 20980403, PMCID PMC2976279. Fetched: Europe PMC record via the "
            "'represses Sox17' full-text phrase search.",
            "reason": "EXCLUDED by rule 3, wrong stage. The pre-implantation embryo is before "
            "gastrulation and the SOX17 at issue there is extraembryonic/primitive "
            "endoderm, not the definitive endoderm this module's SOX17 stands for. Klf5 "
            "is not a mesoderm-lineage factor.",
            "excluded_from": "the closed candidate list",
        },
        {
            "id": "Y4",
            "item": "SOX2 repressing SOX17",
            "citation": "none found at the level rule 1 requires. Teo AKK et al. (2011) Genes Dev "
            "25:238-250 (fetched, abstract) places SOX2 with NANOG and OCT4 as "
            "pluripotency factors that ACTIVELY DIRECT endoderm differentiation by "
            "controlling EOMES -- the opposite sign to a repression of SOX17.",
            "reason": "EXCLUDED twice over. By rule 2, no fetched source states a SOX2 repression "
            "of SOX17; by rule 1, SOX2 is the ectoderm fate in this module and not a "
            "mesoderm-lineage factor. Named because it is the one candidate that would "
            "have needed no new gene, and because a lane under pressure to find a brake "
            "would reach for it.",
            "excluded_from": "the closed candidate list",
        },
        {
            "id": "Y5",
            "item": "SNAI1, TBX6, CDX2 and MSGN1 repressing SOX17",
            "citation": "NONE FOUND. Each was searched for a SOX17 or endoderm-program repression "
            "and each search returned a different target instead: Snail1 downregulates "
            "ectodermal genes within the mesoderm and reciprocally represses Sox3; Tbx6 "
            "represses Sox2 by inactivating enhancer N1; Cdx2 interacts repressively "
            "with Sox2 at the stomach-intestine border; Msgn1 switches off progenitor "
            "maintenance genes. None of the four returned a source for repressing "
            "SOX17.",
            "reason": "EXCLUDED by rule 2: candidates without a citation are excluded and listed "
            "as excluded. They are mesoderm-lineage repressors, which is why they were "
            "searched; their published target is SOX2 or SOX3, not SOX17. A repression "
            "of SOX2 is a rule this module already carries from Tbxt.",
            "excluded_from": "the closed candidate list",
        },
        {
            "id": "Y6",
            "item": "Vox / Vent repressing endoderm development",
            "citation": "the zebrafish ventral repressor Vox, reported to repress endoderm "
            "development by interacting with Casanova and Pou2 (Development 140:1090). "
            "Found by search; the record was seen but the paper's full text was not "
            "fetched, because the entry is excluded before that matters.",
            "reason": "EXCLUDED: Vox/Vent have no human orthologue, and every gene in this module "
            "is human (evidence 'curated GENCODE 50', proteins UniProt P48431, O15178, "
            "Q9H6I2, Q96S42). Adding a gene with no human counterpart to a human module "
            "would be an unsourced act whatever the paper says. Noted as the one "
            "excluded entry whose exclusion is about species rather than lineage or "
            "stage.",
            "excluded_from": "the closed candidate list",
        },
        {
            "id": "Y7",
            "item": "ZFP281 and NANOG restricting SOX17",
            "citation": "ZFP281 recruiting polycomb repressive complex 2 to restrict "
            "extraembryonic endoderm potential in embryonic stem cells (Europe PMC "
            "record seen via the 'represses Sox17' phrase search); Li Z et al. (2022) "
            "The balance between NANOG and SOX17 mediated by TET proteins regulates "
            "specification of human primordial germ cell fate. Cell Biosci 12:181, PMID "
            "36333732 (Europe PMC record seen).",
            "reason": "EXCLUDED by rule 3, wrong lineage. Both act in pluripotency or the "
            "primordial germ cell program, not in the mesoderm lineage. A repression of "
            "SOX17 by a pluripotency factor is a brake on differentiation, not a brake "
            "the mesoderm lineage applies to the endoderm one.",
            "excluded_from": "the closed candidate list",
        },
        {
            "id": "Y8",
            "item": "every exclusion lane-mesoderm registered: X1 to X5",
            "citation": "0790d51, carried unchanged.",
            "reason": "EXCLUDED, inherited without re-litigation: the unsourced thresholds, the "
            "expected_* proportions, nodal_max and decay_length, the runtime's a = 1.0 "
            "convention, and leaving out Nodal activates SOX17. No run in this file "
            "touches any of them.",
            "excluded_from": "the closed candidate list",
        },
    ],
    "prior_lane_verification": {
        "what_was_done": "each of lane-mesoderm's seven citations was looked up this session "
        "through " + FETCH_ROUTE + " and compared, field by field, with the "
        "string its registration carries.",
        "bibliography_verdict": "ALL SEVEN CONFIRMED. Every author list, year, journal, volume and "
        "page range matches the fetched record exactly, including the page "
        "ranges 13:3185-3190, 15:121-133, 13:4469-4481, 1:37-49, 5:62-67, "
        "129:3597-3608, 135:501-511, 25:238-250, 13:1084-1091 and "
        "3:745-756. Nothing was mis-cited and nothing was not found.",
        "entries": [
            {
                "id": "W1",
                "verdict": "CONFIRMED, and the abstract states the claim",
                "fetched": "abstract",
                "record": "Yamaguchi TP, Takada S, Yoshikawa Y, Wu N, McMahon AP (1999) Genes Dev "
                "13:3185-3190, PMID 10617567",
                "supports_the_rule": "yes: 'Transgenic analysis of the T promoter identifies T "
                "(Brachyury) as a direct transcriptional target of the Wnt "
                "signaling pathway', in cells fated to form paraxial "
                "mesoderm.",
            },
            {
                "id": "W2",
                "verdict": "CONFIRMED bibliographically; the claim is weaker than 'direct'",
                "fetched": "abstract",
                "record": "Martin BL, Kimelman D (2008) Dev Cell 15:121-133, PMID 18606146",
                "supports_the_rule": "the sign, yes: Ntl and Bra 'are required for and can induce "
                "expression of the canonical Wnts wnt8 and wnt3a'. The "
                "abstract does NOT say direct, and the system is zebrafish.",
            },
            {
                "id": "F1",
                "verdict": "CONFIRMED, both papers",
                "fetched": "abstracts",
                "record": "Isaacs HV, Pownall ME, Slack JM (1994) EMBO J 13:4469-4481, PMID "
                "7925289; Ciruna B, Rossant J (2001) Dev Cell 1:37-49, PMID 11703922",
                "supports_the_rule": "yes: 'eFGF and Xbra can activate each other's expression' "
                "(Xenopus), and FGFR1 acts 'in mesoderm cell fate "
                "specification through positive regulation of Brachyury and "
                "Tbx6' (mouse streak).",
            },
            {
                "id": "F2",
                "verdict": "CONFIRMED",
                "fetched": "abstract",
                "record": "Schulte-Merker S, Smith JC (1995) Curr Biol 5:62-67, PMID 7535172",
                "supports_the_rule": "yes: 'Brachyury activates FGF gene expression, and FGF "
                "reciprocally maintains Brachyury expression'.",
            },
            {
                "id": "M1",
                "verdict": "CONFIRMED as a paper; ITS TWO RULES ARE NOT CONFIRMED BY IT",
                "fetched": "abstract (Europe PMC) and the publisher's article page; the full text "
                "was NOT reachable without a subscription",
                "record": "Hart AH, Hartley L, Sourris K, Stadler ES, Li R, Stanley EG, Tam PP, "
                "Elefanty AG, Robb L (2002) Development 129:3597-3608, PMID 12117810",
                "supports_the_rule": "NO, not at the level the registration states. M1 attributes "
                "TWO rules to this paper, 'Nodal activates MIXL1' and "
                "'Mixl1 inhibits SOX2'. The fetched abstract mentions "
                "NEITHER Nodal NOR Sox2: it is a null-mutant morphogenesis "
                "paper reporting a thickened primitive streak, absent heart "
                "and gut, deficient paraxial mesoderm, and that 'Mixl1 "
                "activity is most crucial for ENDODERMAL differentiation'. "
                "It sources Mixl1 as a required mesendoderm factor; it does "
                "not source either of M1's two rules, and its own emphasis "
                "is endodermal. THIS MATTERS BECAUSE M1 IS THE ENTRY THAT "
                "CLEARS THE WEAK BOUND. The prior row's M1 reading should "
                "carry that its two rules are inferred from the factor's "
                "role and are not stated by the cited paper.",
            },
            {
                "id": "E1",
                "verdict": "CONFIRMED, all three papers; one of the two rules is supported, the "
                "other only indirectly",
                "fetched": "abstracts of Arnold 2008 and Teo 2011; FULL TEXT of Costello 2011",
                "record": "Arnold SJ, Hofmann UK, Bikoff EK, Robertson EJ (2008) Development "
                "135:501-511, PMID 18171685; Teo AK, Arnold SJ, Trotter MW, Brown S, "
                "Ang LT, Chng Z, Robertson EJ, Dunn NR, Vallier L (2011) Genes Dev "
                "25:238-250, PMID 21245162; Costello I, Pimeisl IM, Draeger S, Bikoff "
                "EK, Robertson EJ, Arnold SJ (2011) Nat Cell Biol 13:1084-1091, PMID "
                "21822279",
                "supports_the_rule": "'Eomes activates SOX17' is well sourced: Teo reports EOMES "
                "'marks the onset of endoderm specification' and 'interacts "
                "with SMAD2/3 to initiate the transcriptional network "
                "governing endoderm formation', and Arnold reports 'Eomes is "
                "required for specification of the definitive endoderm "
                "lineage'. 'Nodal activates EOMES' is NOT stated by either "
                "abstract; the nearest fetched support is Teo's EOMES-SMAD2/3 "
                "interaction, SMAD2/3 being the NODAL transducer, and "
                "Arnold's Eomes/Nodal double heterozygotes. Costello's "
                "fetched full text adds, in the lane's favour, 'Eomes "
                "expression within the primitive streak marks the earliest "
                "cardiac mesoderm and promotes formation of cardiovascular "
                "progenitors by directly activating the bHLH transcription "
                "factor Mesp1'.",
            },
            {
                "id": "A1",
                "verdict": "CONFIRMED, and the antagonists are named exactly as the lane named them",
                "fetched": "abstract",
                "record": "Perea-Gomez A, Vella FD, Shawlot W, Oulad-Abdelghani M, Chazaud C, Meno "
                "C, Pfister V, Chen L, Robertson E, Hamada H, Behringer RR, Ang SL "
                "(2002) Dev Cell 3:745-756, PMID 12431380",
                "supports_the_rule": "yes, and the structural non-run stands: 'Cerberus-like and "
                "Lefty1 in the anterior visceral endoderm restrict primitive "
                "streak formation to the posterior end of mouse embryos by "
                "antagonizing Nodal signaling'. Both antagonise the LIGAND, "
                "which this program clamps.",
            },
        ],
        "in_repo_check_confirmed": "the prior lane's statement that the repository holds no "
        "regulatory citation for WNT3A->TBXT, FGF->TBXT, MIXL1 or "
        "EOMES was not re-litigated; this lane checked the published "
        "record, not the repository.",
    },
    "header_consequence": {
        "line": 'data/demo/gastrulation.bio\'s header says "Mutual repression makes the choice sharp"',
        "status_before_this_lane": "FALSE as the program stands: SOX17 represses two genes and is "
        "repressed by none, so the repression is one-way and not mutual.",
        "this_lane_does_not_change_it": "no run in this file edits data/demo/gastrulation.bio. "
        "Every variant is written to a temporary file built from "
        "the committed module text, so the header stays as it is "
        "and stays false. Changing it is a change to a shared "
        "demo module and is reported to the coordinator rather "
        "than taken here.",
    },
}

_LIST = {e["id"]: e for e in REGISTRATION["closed_list"]}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.parse_args()

    base = MODULE.read_text()
    run = dict(CENSUS_MODEL_RUN)
    results: dict[str, dict] = {}

    with tempfile.TemporaryDirectory() as tmp:

        def measure(name: str, chosen: list[dict]) -> dict:
            p = Path(tmp) / f"{name}.bio"
            p.write_text(variant_text(base, chosen))
            return evaluate(p, run)

        results["baseline"] = measure("baseline", [])
        for cid, entry in _LIST.items():
            if not entry.get("runnable", True):
                results[cid] = {"no_run": entry["no_run_reason"]}
                continue
            results[cid] = measure(cid, [entry])

    cleared = sorted(k for k, v in results.items() if v.get("clears_bound"))
    reading = (
        "a sourced repressor of SOX17 restores a mesoderm band"
        if [c for c in cleared if c != "baseline"]
        else "no sourced repressor on the list does"
    )
    payload = {
        "pass": "single",
        "run": run,
        "record_every": RECORD_EVERY,
        "bound": {
            "statement": REGISTRATION["bound"]["statement"],
            "share": BOUND_SHARE,
            "last_fraction": BOUND_LAST_FRACTION,
            "not_the_census_range": REGISTRATION["bound"]["deliberately_weak"],
        },
        "registration": REGISTRATION,
        "runs": results,
        "cleared_the_bound": cleared,
        "cleared_count": len(cleared),
        "reported_count": len(results),
        "reading": reading,
        "reading_was_one_of_exactly_two": [
            "a sourced repressor of SOX17 restores a mesoderm band",
            "no sourced repressor on the list does",
        ],
        "carried_verbatim": {
            "census_bound_d55cb19": (
                "a sampled census of one embryo is not a 1-D axis, and the model names no stage"
            ),
            "cs7_verdict": "falsified on all three layers, unchanged by anything in this file",
            "sign_correction": "5cbce26 stands: Tbxt activates SOX17 in every run here, and "
            "Sox17 inhibits TBXT is removed in none of them",
            "m1_reading": REGISTRATION["builds_on"]["readings_carried_unchanged"][0],
            "r1_reading": REGISTRATION["builds_on"]["readings_carried_unchanged"][1],
            "structural_finding": REGISTRATION["builds_on"]["readings_carried_unchanged"][2],
        },
        "wording": (
            "any entry that restores the band means the program lacked a repressor of SOX17; "
            "it does not mean that factor specifies mesoderm"
        ),
        "alphagenome_requests": 0,
        "money": "none: 0 model requests, no paid API. The web was used, which is free, for "
        "citation fetching only; no run below reads the network.",
    }
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "data/demo/gastrulation.bio",
                "version": "in-repo demo module at this revision",
            },
            {
                "accession": "Costello et al. 2011, Nat Cell Biol 13:1084-1091 (PMC4531310)",
                "version": "Europe PMC fullTextXML, fetched 2026-10-02",
            },
            {
                "accession": "Lynch et al. 2025, EMBO J 44:2541-2565 (PMC12048643)",
                "version": "Europe PMC fullTextXML, fetched 2026-10-02",
            },
        ],
        "inputs": [mf.input_entry(MODULE, partition=None)],
        "assembly": "n/a: a simulated regulatory network, no genome coordinates",
        "coordinates": "n/a: a normalised axis x in [0, 1], not genomic intervals",
        "parameters": {
            "pass": "single",
            "model_run": run,
            "record_every": RECORD_EVERY,
            "bound_share": BOUND_SHARE,
            "bound_last_fraction": BOUND_LAST_FRACTION,
            "added_gene": REGISTRATION["number_convention"]["added_gene"],
            "added_activating_rule": REGISTRATION["number_convention"]["added_activating_rule"],
            "added_inhibiting_rule": REGISTRATION["number_convention"]["added_inhibiting_rule"],
            "inherited_from": "scripts/mesoderm_diagnosis.py at 0790d51, imported rather than "
            "restated so the bound and the convention are the same objects",
        },
        "exclusions": [x["item"] + " -- " + x["reason"] for x in REGISTRATION["excluded"]],
        "partitions": "n/a: no evaluation split; every entry of the closed list is reported",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE),
    }
    path = save_result("sox17_repressor_single", payload)
    print(f"single: {len(results)} reported, cleared {len(cleared)}: {', '.join(cleared) or 'none'}")
    for k, v in results.items():
        if "no_run" in v:
            print(f"  {k:10s} no run: {v['no_run'][:64]}")
            continue
        f = v["final"]
        print(
            f"  {k:10s} ecto {f['ectoderm']:.3f} meso {f['mesoderm']:.3f} endo {f['endoderm']:.3f}"
            f"  tail_min {v['mesoderm_tail_min']:.4f}  {'CLEARS' if v['clears_bound'] else '-'}"
        )
    print(f"  reading: {reading}")
    print(f"  -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
