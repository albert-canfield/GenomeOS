# SPDX-License-Identifier: AGPL-3.0-or-later
"""Area E, the pie-1 discrepancy: a closed, sourced list of candidate corrections, registered in this
file before any rule of the worm program changes, then run against the wild type and every knockout.

THE DISCREPANCY (lane-wormcontact, 6bdecd4 registration, 09be446 result). In embryo_nolookup.bio the
pie-1 knockout loses the germ line, as published, but instead of the published extra intestine it makes
no intestine at all: `presents_mom2` in founders_contacts.bio needs PIE-1 in P2, so without PIE-1 no Wnt
reaches EMS at the contact and E takes the MS fate. The same gate is the program's only route to the
E-to-MS change Du et al. 2014 observe in par-2, where P2 keeps no PIE-1 (div_P1_no_par2).

THE DISCIPLINE (the mesoderm lane's, 0790d51 registration, b1e8839 result). `REGISTRATION` below is the
registration: a closed list, each candidate with its sources as fetched on 2026-10-10 (quoted, never
recalled), the exact text edits it makes, and what it is predicted to do to the wild type and to every
knockout at the founder and the terminal level. It is committed before the first run; the run refuses
unless this file is committed and clean, and `tests/test_pie1.py` pins the registration's digest, so the
list cannot grow, shrink or change its predictions without a commit that says so. Candidates not on the
list are refused by the runner. Exclusions are named with the reason, not dropped.

The falsifier is fixed in advance: a candidate that changes the wild type's 56 / 56 / 0, or changes ANY
other knockout at the founder or the terminal level, or loses any of Du et al. 2014's observed founder
changes the baseline reproduces, is rejected whatever it does for pie-1. Matching pie-1 by tuning is not
a result: every edit is the mechanism its source states, and no number is chosen.

What was read before writing this: the program files, the runtime's division, express and contact
semantics (genomeos/runtime/body.py, to derive the predictions), the measured contact table's rows for
C and P3, and the baseline pie-1 knockout's founder types, which 09be446 already published (P2
EMSPrecursor, ABp ABaPrecursor, E MSPrecursor). No candidate was grown before this file was committed.

    uv run python scripts/pie1_registration.py --digest   # the registration's digest
    uv run python scripts/pie1_registration.py            # the run (refuses unless committed and clean)

No model calls, no network, about a minute on one core.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
from datetime import date
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.lang import parse_file  # noqa: E402
from genomeos.organism import digital_development as dd  # noqa: E402
from genomeos.organism import provenance as pv  # noqa: E402
from genomeos.organism.diff import compare  # noqa: E402
from genomeos.organism.reference import ReferenceLineage  # noqa: E402
from genomeos.results import save_result  # noqa: E402
from genomeos.runtime.body import Body  # noqa: E402
from scripts.celegans_nolookup import (  # noqa: E402
    EXPECTED_FOUNDERS,
    KNOCKOUT_FOUNDERS,
    KNOCKOUT_TERMINAL,
    PUBLISHED_TEST,
    TISSUE_RULES,
    WT_COUNTS,
    counts,
    decided,
    founder_types,
    terminal_cells,
)

ORG = "data/organisms/celegans"
FOUNDERS = "founders_contacts.bio"
NOLOOKUP = "embryo_nolookup.bio"
CONTACTS = "embryo_contacts.bio"
#: the last commit that touched the worm's programs before this registration: every arm, the baseline
#: included, is grown from these files, so the run stays the registered one after a candidate is applied
BASE_COMMIT = "1c2ae677a7e0eecbe36ae47d1e8a3c947ae171e0"
TABLE = Path("data/results/celegans_digital_development.json")
RESULT_NAME = "pie1_corrections"
HATCH = 800.0
SEEDS = list(range(20))
OWN_CODE = ("scripts/pie1_registration.py", "tests/test_pie1.py")
#: the knockouts every arm runs: the eight registered at 6bdecd4, the two PAR knockouts the contact
#: program's Du et al. 2014 comparison covers, and POS-1, reported and not scored (see the registration)
KNOCKOUTS = (
    "POP-1",
    "SKN-1",
    "PIE-1",
    "APX-1",
    "GLP-1",
    "PAL-1",
    "MOM-2",
    "MOM-5",
    "PAR-2",
    "PAR-3",
    "POS-1",
)
TERMINAL_KNOCKOUTS = tuple(KNOCKOUT_TERMINAL)  # the eight
NOT_SCORED = ("POS-1",)
FOUNDER_NAMES = ("ABa", "ABp", "EMS", "P2", "MS", "E", "C", "D", "P4")

# ---- the sources, as fetched on 2026-10-10 ------------------------------------------------------------

ABSTRACTS = "abstract fetched from NCBI E-utilities (efetch db=pubmed, rettype=abstract), 2026-10-10"
SOURCES: dict[str, dict[str, Any]] = {
    "Mello1992": {
        "cite": "Mello CC, Draper BW, Krause M, Weintraub H, Priess JR (1992) Cell 70:163-176, PMID 1623520, "
        "doi:10.1016/0092-8674(92)90542-k",
        "fetched": ABSTRACTS + "; full text NOT retrievable (cell.com and sciencedirect.com 403; Unpaywall "
        "and Semantic Scholar list no open copy)",
        "quote": "In pie-1 mutants one additional posterior blastomere adopts an MS-like fate, and in mex-1 "
        "mutants four additional anterior blastomeres adopt an MS-like fate.",
    },
    "Seydoux1996": {
        "cite": "Seydoux G, Mello CC, Pettitt J, Wood WB, Priess JR, Fire A (1996) Nature 382:713-716, "
        "PMID 8751441",
        "fetched": ABSTRACTS,
        "quote": "after each asymmetrical division, embryonically transcribed RNAs are detected in somatic, "
        "but not germline, blastomeres. This asymmetry depends on the activity of the germline specific "
        "factor, PIE-1. In the absence of PIE-1, embryonically transcribed RNAs are detected in both "
        "somatic and germline blastomeres.",
    },
    "Mello1996": {
        "cite": "Mello CC, Schubert C, Draper B, Zhang W, Lobel R, Priess JR (1996) Nature 382:710-712, "
        "PMID 8751440",
        "fetched": ABSTRACTS,
        "quote": "Totipotent germline blastomeres in Caenorhabditis elegans contain, but do not respond to, "
        "factors that promote somatic differentiation in other embryonic cells. Mutations in the maternal "
        "gene pie-1 result in the germline blastomeres adopting somatic cell fates.",
    },
    "Tenenhaus1998": {
        "cite": "Tenenhaus C, Schubert C, Seydoux G (1998) Dev Biol 200:212-224, PMID 9705228",
        "fetched": ABSTRACTS,
        "quote": "mutants that disrupt the asymmetric divisions of germline blastomeres mislocalize PIE-1, "
        "and activate RNAPII-H5 expression in the germ lineage. In contrast, mutants that alter somatic "
        "cell identities do not affect PIE-1 localization or RNAPII-H5 expression. Our observations "
        "suggest that PIE-1 represses mRNA transcription in each germline blastomere in a "
        "concentration-dependent manner.",
    },
    "Tenenhaus2001": {
        "cite": "Tenenhaus C, Subramaniam K, Dunn MA, Seydoux G (2001) Genes Dev 15:1031-1040, PMID 11316796",
        "fetched": ABSTRACTS,
        "quote": "PIE-1 has a second function in germ cells; it is required for efficient expression of the "
        "maternally encoded Nanos homolog NOS-2.",
    },
    "Thorpe1997": {
        "cite": "Thorpe CJ, Schlesinger A, Carter JC, Bowerman B (1997) Cell 90:695-705, PMID 9288749",
        "fetched": ABSTRACTS,
        "quote": "mom-1, mom-2, and mom-3 are required in the signaling cell, P2, while mom-4 is required in "
        "EMS. P2 signaling downregulates an HMG domain protein, POP-1, in one EMS daughter.",
    },
    "Rocheleau1997": {
        "cite": "Rocheleau CE, Downs WD, Lin R, Wittmann C, et al. (1997) Cell 90:707-716, PMID 9288750",
        "fetched": ABSTRACTS,
        "quote": "We present evidence that there may be partially redundant inputs into endoderm "
        "specification",
    },
    "Oldenbroek2013": {
        "cite": "Oldenbroek M, Robertson SM, Guven-Ozkan T, Spike C, Greenstein D, Lin R (2013) Development "
        "140:4614-4623, PMID 24131629, PMC3817945, doi:10.1242/dev.096313",
        "fetched": "abstract from NCBI efetch db=pmc; results text from journals.biologists.com through a "
        "fetch tool that returns quoted sentences (each quote below was asked for verbatim and re-asked "
        "once; one paraphrase the tool produced on the first pass was rejected and re-checked), 2026-10-10",
        "quote": [
            "Here, we characterize the translational regulation of maternally supplied mom-2 mRNA, which "
            "encodes a Wnt ligand essential for two separate cell-cell interactions in early embryos.",
            "Depletion of any of the three germline blastomere-enriched proteins, MEX-1, PIE-1 and POS-1, "
            "resulted in a reduction, but not abolishment, of GFP::H2Bmom-2 expression: mex-1(zu121), "
            "50±3% reduction in P3 (n=13); pie-1(zu154), 54±4% reduction in P3 (n=13); pos-1(zu148), "
            "33.1±4% reduction in P3 (n=12).",
            "pos-1(RNAi) in either the pie-1(zu154) or mex-1(zu121) strain resulted in no detectable "
            "GFP::H2Bmom-2 in embryos (n=73 and n=59, respectively).",
            "Combined depletion of pie-1 and mex-1 did not result in a more severe defect in GFP::H2Bmom-2 "
            "expression than either single mutant.",
            "We conclude that PIE-1, MEX-1 and POS-1 function together to promote the translation of "
            "gfp::h2bmom-2 in germline blastomeres.",
            "The P2 blastomere also has high levels of MEX-1, PIE-1 and POS-1.",
            "P2 and P3 have a high level of the negative regulator SPN-4, as well as the three positive "
            "regulators PIE-1, MEX-1 and POS-1",
            "the reporter: nuclear GFP::histone H2B under the control of the germline-specific pie-1 "
            "promoter and the mom-2 3'UTR; not detected until the 4-cell stage, initially at a very low "
            "level in the P2 blastomere",
            "All pos-1(zu148) embryos fail to produce endoderm, whereas mex-1(zu121) mutant embryos have a "
            "cold-sensitive, incompletely penetrant defect in endoderm production.",
        ],
    },
    "Mango1994": {
        "cite": "Mango SE, Thorpe CJ, Martin PR, Chamberlain SH, Bowerman B (1994) Development "
        "120:2305-2315, "
        "PMID 7925031",
        "fetched": ABSTRACTS,
        "quote": "in pie-1 embryos, APp [sic, ABp] fails to produce its characteristic cell types and "
        "instead adopts a fate similar to that of ABa. We conclude that the changed identity of P2 in "
        "pie-1 embryos prevents the P2/ABp interaction.",
    },
    "Tabara1999": {
        "cite": "Tabara H, Hill RJ, Mello CC, Priess JR, Kohara Y (1999) Development 126:1-11, PMID 9834181",
        "fetched": ABSTRACTS,
        "quote": "During early embryogenesis, pos-1 mRNA and POS-1 protein are present predominantly in the "
        "germ precursors. POS-1 is a novel protein with two copies of a CCCH finger motif previously "
        "described in the germline proteins PIE-1 and MEX-1 ... The earliest defect detected in pos-1 "
        "mutants is the failure to express APX-1 protein from maternally provided apx-1 mRNA",
    },
    "Guedes1997": {
        "cite": "Guedes S, Priess JR (1997) Development 124:731-739, PMID 9043088",
        "fetched": ABSTRACTS,
        "quote": "We show here that MEX-1 also is present in germline blastomeres and is a P granule "
        "component ... We further show that MEX-1 is required to restrict PIE-1 expression and activity "
        "to the germline blastomeres during the early embryonic cleavages.",
    },
    "Bei2002": {
        "cite": "Bei Y, Hogan J, Berkowitz LA, Soto M, Rocheleau CE, et al. (2002) Dev Cell 3:113-125, "
        "PMID 12110172",
        "fetched": ABSTRACTS,
        "quote": "Although Wnt signaling contributes to this polarizing interaction, no mutants identified "
        "to date abolish P2/EMS signaling. ... src-1 and mes-1 mutants strongly enhance endoderm and EMS "
        "spindle rotation defects associated with Wnt pathway mutants. ... Wnt and Src signaling function "
        "in parallel to control developmental outcomes within a single responding cell.",
    },
    "Goldstein1992": {
        "cite": "Goldstein B (1992) Nature 357:255-257, PMID 1589023",
        "fetched": ABSTRACTS,
        "quote": "Gut differentiation is rescued by recombining EMS with its posterior neighbour (P2), but "
        "not by recombining EMS with one or both of the other two cells of the four-cell stage.",
    },
    "Du2014": {
        "cite": "Du Z, Santella A, He F, Tiongson M, Bao Z (2014) Cell 156:359-372, PMID 24439388; the "
        "Fate_regulation table distilled in data/results/celegans_digital_development.json",
        "fetched": "full text from NCBI efetch db=pmc (PMC3998820); table rows read from the committed "
        "distillation, 2026-10-10",
        "quote": [
            "pie-1: {ABp: ABa, P2: EMS} (no change of E is listed)",
            "par-2: {ABp: ABa, ABar: ABal, P2: EMS, E: MS}",
            "mom-2: {ABar: ABal, ABpra: ABpla, E: MS}",
            "the Wnt pathway components functioning in the early embryo are provided maternally",
        ],
    },
    "Boyd1996": {
        "cite": "Boyd L, Guo S, Levitan D, Stinchcomb DT, Kemphues KJ (1996) Development 122:3075-3084, "
        "PMID 8898221",
        "fetched": ABSTRACTS,
        "quote": "although par-2 activity is not required for posterior localization of P granules at the "
        "one-cell stage, it is required for proper cortical association of P granules in P1.",
    },
    "Cheeks2004": {
        "cite": "Cheeks RJ, Canman JC, Gabriel WN, Meyer N, Strome S, Goldstein B (2004) Curr Biol "
        "14:851-862, "
        "PMID 15186741",
        "fetched": ABSTRACTS,
        "quote": "PAR-1 functions to refine the asymmetric localization of germline ribonucleoprotein "
        "complexes by selectively stabilizing only those complexes that reach the PAR-1-enriched "
        "posterior cell cortex",
    },
    "Lin1998": {
        "cite": "Lin R, Hill RJ, Priess JR (1998) Cell 92:229-239, PMID 9458047",
        "fetched": ABSTRACTS,
        "quote": "Each anterior cell resulting from an anterior-posterior division appears to have a higher "
        "level of nuclear POP-1 protein than does its posterior sister. Genes in the C. elegans Wnt "
        "pathway are required for this inequality in POP-1 levels.",
    },
    "Bowerman1993": {
        "cite": "Bowerman B, Draper BW, Mello CC, Priess JR (1993) Cell 74:443-452, PMID 8348611",
        "fetched": ABSTRACTS,
        "quote": "pie-1(+) activity may function to regulate the activity of Skn-1 protein in the "
        "descendants of the posterior blastomere P1.",
    },
}

# ---- the candidates: exact text edits, each `old` present exactly once when it is applied -------------

PRESENTS_LAG2 = "decision presents_lag2 {"
POS1_COMMENT = (
    "# MOM-2 is maternal mRNA translated in the germline blastomeres, and PIE-1 is one of three proteins\n"
    "# there that raise its translation, with MEX-1 and POS-1 (Oldenbroek et al. 2013): losing PIE-1 alone\n"
    "# reduces it by about half, losing POS-1 as well leaves none. So P2 presents MOM-2 while it holds\n"
    "# either PIE-1 or POS-1. MEX-1 is not modelled: it is in the same cells as PIE-1 (Guedes & Priess\n"
    "# 1997), and losing PIE-1 and MEX-1 together is no worse than losing one (Oldenbroek et al. 2013).\n"
)
CANDIDATES: dict[str, dict[str, Any]] = {
    "C1": {
        "what": "presents_mom2 loses its PIE-1 condition: P2 presents MOM-2 whatever it holds",
        "mechanism": "mom-2 is maternal mRNA (Oldenbroek et al. 2013; Du et al. 2014, 'the Wnt pathway "
        "components ... are provided maternally'); PIE-1 represses NEW transcription in the germ line "
        "(Seydoux et al. 1996), so it cannot be what makes a maternal ligand present; and without PIE-1 "
        "the mom-2 3'UTR reporter is reduced, 'not abolishment' (Oldenbroek et al. 2013: pie-1(zu154) 54% "
        "in P3). Read as presence, reduced is present. The same paper's double depletion contradicts the "
        "candidate for pie-1 with pos-1, which no registered knockout is",
        "sources": ["Oldenbroek2013", "Seydoux1996", "Du2014", "Thorpe1997", "Mello1992"],
        "edits": [
            (
                FOUNDERS,
                "decision presents_mom2 { action: express; when: cell = P2, PIE-1 = present; sets: MOM-2;",
                "decision presents_mom2 { action: express; when: cell = P2; sets: MOM-2;",
            ),
            (
                FOUNDERS,
                '(MOM-2 from P2)"; confidence: 0.9 }',
                "(MOM-2 from P2); Oldenbroek et al. 2013, Development 140:4614 (mom-2 is maternal mRNA; "
                "without PIE-1 its 3'UTR reporter is reduced by 54% in P3, not abolished)\"; "
                "confidence: 0.9 }",
            ),
        ],
    },
    "C2": {
        "what": "a second, Wnt-parallel route from P2 to EMS: P2 presents MES-1, EMS reads it at the "
        "contact through SRC-1, and either read polarises EMS",
        "mechanism": "Bei et al. 2002: SRC-1 and MES-1 signal between P2 and EMS in parallel with Wnt, and "
        "no single mutant abolishes the P2/EMS signal; Rocheleau et al. 1997: 'partially redundant inputs "
        "into endoderm specification'. Neither abstract gates MES-1 on PIE-1, so none is written",
        "sources": ["Bei2002", "Rocheleau1997", "Goldstein1992"],
        "edits": [
            (
                FOUNDERS,
                PRESENTS_LAG2,
                "decision presents_mes1 { action: express; when: cell = P2; sets: MES-1; evidence: "
                "experimental "
                '"Bei et al. 2002, Dev Cell 3:113 (MES-1 and SRC-1 signal between P2 and EMS, in '
                "parallel with "
                'Wnt)"; confidence: 0.7 }\n' + PRESENTS_LAG2,
            ),
            (
                FOUNDERS,
                "signal Notch2_at_contact {",
                "signal Src_at_contact { mode: contact; ligand: MES-1; receptor: SRC-1; reads: amount; from: "
                'MES-1 = present; sets: Src_amount; evidence: experimental "Bei et al. 2002, Dev '
                'Cell 3:113"; '
                "confidence: 0.7 }\nsignal Notch2_at_contact {",
            ),
            (
                FOUNDERS,
                "decision fate_EMS { ",
                "decision fate_EMS_polarised_src { action: differentiate; when: cell_type = EMSPrecursor, "
                'Src_amount = >=0.5; to: EMSPolarised; evidence: experimental "Bei et al. 2002, Dev '
                "Cell 3:113 "
                "(Wnt and Src act in parallel in EMS); Rocheleau et al. 1997, Cell 90:707 "
                "(partially redundant "
                'inputs into endoderm specification)"; confidence: 0.7 }\ndecision fate_EMS { ',
            ),
        ],
    },
    "C3": {
        "what": "P2 presents MOM-2 while it holds PIE-1 OR POS-1; POS-1 is a maternal factor segregated to "
        "the germ precursors as PIE-1 is, and lost from both daughters of P1 without PAR-2 as PIE-1 is",
        "mechanism": "Oldenbroek et al. 2013: MEX-1, PIE-1 and POS-1 'function together to promote the "
        "translation' of mom-2 in germline blastomeres; each alone is a 'reduction, but not abolishment'; "
        "pos-1(RNAi) in pie-1(zu154) leaves 'no detectable' reporter (n=73); P2 and P3 hold high POS-1. "
        "Tabara et al. 1999: POS-1 is 'present predominantly in the germ precursors'. The par-2 branch is a "
        "PREMISE, stated as one: the program's div_P1_no_par2 already says the germline determinants are "
        "not retained on the posterior side without PAR-2 (Boyd et al. 1996; Cheeks et al. 2004), and C3 "
        "applies that to POS-1, the second such protein; no fetched source measures POS-1 in a par-2 P2",
        "sources": [
            "Oldenbroek2013",
            "Tabara1999",
            "Guedes1997",
            "Boyd1996",
            "Cheeks2004",
            "Tenenhaus1998",
            "Du2014",
        ],
        "edits": [
            (
                FOUNDERS,
                PRESENTS_LAG2,
                "decision presents_mom2_pos1 { action: express; when: cell = P2, POS-1 = present; "
                "sets: MOM-2; "
                'evidence: experimental "Oldenbroek et al. 2013, Development 140:4614 (mom-2 is '
                "maternal mRNA "
                "translated in the germline blastomeres; MEX-1, PIE-1 and POS-1 each raise its translation; "
                "losing PIE-1 reduces it by 54% in P3 and does not abolish it, and pos-1(RNAi) in "
                "pie-1(zu154) "
                'leaves none detectable)"; confidence: 0.7 }\n' + PRESENTS_LAG2,
            ),
            (
                FOUNDERS,
                "# --- which cells present a ligand: measured, and independent of who hears it "
                "---------------------\n",
                "# --- which cells present a ligand: measured, and independent of who hears it "
                "---------------------\n" + POS1_COMMENT,
            ),
            (
                FOUNDERS,
                "asymmetric: PIE-1 -> P1, SKN-1 -> P1,",
                "asymmetric: PIE-1 -> P1, POS-1 -> P1, SKN-1 -> P1,",
            ),
            (
                FOUNDERS,
                '(GLP-1 protein in AB only)"',
                "(GLP-1 protein in AB only); Tabara et al. 1999, Development 126:1 (pos-1 mRNA and POS-1 "
                'protein are predominantly in the germ precursors)"',
            ),
            (
                FOUNDERS,
                "asymmetric: PIE-1 -> P1, PAR-2 -> P1,",
                "asymmetric: PIE-1 -> P1, POS-1 -> P1, PAR-2 -> P1,",
            ),
            (FOUNDERS, "asymmetric: PIE-1 -> P2;", "asymmetric: PIE-1 -> P2, POS-1 -> P2;"),
            (FOUNDERS, "asymmetric: PIE-1 -> none;", "asymmetric: PIE-1 -> none, POS-1 -> none;"),
            (
                FOUNDERS,
                'all three of which follow from P2 losing its identity)"',
                "all three of which follow from P2 losing its identity; POS-1 is taken to follow PIE-1 here, "
                "as the second germline-blastomere protein that raises mom-2 translation (Oldenbroek et al. "
                '2013), a premise no source in hand measures in par-2)"',
            ),
            (
                FOUNDERS,
                "asymmetric: PIE-1 -> P3, SKN-1 -> P3;",
                "asymmetric: PIE-1 -> P3, POS-1 -> P3, SKN-1 -> P3;",
            ),
            (
                FOUNDERS,
                "asymmetric: PIE-1 -> P4, SKN-1 -> P4;",
                "asymmetric: PIE-1 -> P4, POS-1 -> P4, SKN-1 -> P4;",
            ),
            (NOLOOKUP, "factors: SKN-1, PIE-1, PAL-1,", "factors: SKN-1, PIE-1, POS-1, PAL-1,"),
            (CONTACTS, "factors: SKN-1, PIE-1, PAL-1,", "factors: SKN-1, PIE-1, POS-1, PAL-1,"),
        ],
    },
}
#: run because two of the three act on the same event (EMS's polarisation): every pair, every one reported
PAIRS = (("C1", "C2"), ("C1", "C3"), ("C2", "C3"))

EXCLUDED = {
    "X1": {
        "what": "MOM-2 presented at a reduced amount without PIE-1 (about 0.46 of the wild type, from the "
        "54% reduction)",
        "why": "a reporter's fluorescence in P3 at the 8-cell stage has no mapping to the contact-weighted "
        "amount and the 0.5 threshold the EMS read uses; choosing the number is tuning a free parameter",
    },
    "X2": {
        "what": "presents_apx1 loses its PIE-1 condition",
        "why": "contradicted: Mango et al. 1994 ('the changed identity of P2 in pie-1 embryos prevents the "
        "P2/ABp interaction') and Du et al. 2014's pie-1 row (ABp adopts ABa)",
    },
    "X3": {
        "what": "a pie-1 P2 divides as EMS does: its daughters take MS-like and E-like fates, which is where "
        "the published EXTRA intestine would come from",
        "why": "no fetched source states it at the level a rule needs: Mello et al. 1992's abstract says one "
        "additional posterior blastomere is MS-like and does not name it or say which cell makes the extra "
        "intestine (the full text was not retrievable, 403); a POP-1 split at P2's division would rest on "
        "Lin et al. 1998's general statement about anterior-posterior divisions, and no fetched text says "
        "P2's is one; C keeping SKN-1 in pie-1 has no fetched source (Bowerman et al. 1993's abstract says "
        "only that pie-1 regulates SKN-1 activity in P1's descendants). Excluded for want of a source",
    },
    "X4": {
        "what": "a MOM-2 gate on PAR-2 itself, or on MEX-1, chosen because it separates par-2 from pie-1",
        "why": "no fetched source measures MOM-2 or mom-2 translation in a par-2 embryo; a gate chosen "
        "because it separates the two rows would be fitted to them. MEX-1 is not modelled (C3 says why)",
    },
    "X5": {
        "what": "C3 with POS-1 kept by P2 (or shared by both daughters) in par-2",
        "why": "no fetched source says that either; the program's par-2 rule already states the class "
        "premise and C3 follows it. Named so the dependence is visible: under X5, par-2's E-to-MS would be "
        "lost exactly as under C1",
    },
}

#: what each candidate is predicted to do, derived from the rules and the runtime before any run; every
#: knockout not named is predicted to be identical to the baseline at the founder and the terminal level
PREDICTIONS: dict[str, dict[str, Any]] = {
    "baseline": {
        "wild_type": "56 / 56 / 0 at the means and in seeds 0-19 (by_founder D Muscle 20, E Intestine 34, P4 "
        "Reproductive 2); the eight founder types as expected",
        "PIE-1": "founders P2 EMSPrecursor, ABp ABaPrecursor, E MSPrecursor; terminal 0/20/0/20",
        "PAR-2": "E MSPrecursor, ABp ABaPrecursor, P2 EMSPrecursor",
        "du2014": "11 reproduced, 0 missed on embryo_contacts.bio (as committed) and on embryo_nolookup.bio",
    },
    "C1": {
        "wild_type": "unchanged: P2 holds PIE-1 in the wild type, so the dropped condition was true there",
        "PIE-1": "E EPrecursor; terminal 34/20/0/54; P2 EMSPrecursor and ABp ABaPrecursor unchanged; P2's "
        "daughters do not carry MOM-2 on, because reader.bio's read_C and read_P3 replace the measured set",
        "PAR-2": "E EPrecursor: Du et al. 2014's par-2 E-to-MS is MISSED on both programs (10 of 11)",
        "verdict": "rejected: it changes PAR-2",
    },
    "C2": {
        "wild_type": "unchanged: EMS is polarised by both routes",
        "PIE-1": "E EPrecursor; terminal 34/20/0/54",
        "MOM-2": "E EPrecursor (Src polarises EMS without Wnt): founder E-to-MS lost, terminal 34/20/2/56",
        "MOM-5": "as MOM-2",
        "PAR-2": "E EPrecursor: par-2's E-to-MS missed",
        "verdict": "rejected: it changes MOM-2, MOM-5 and PAR-2",
    },
    "C3": {
        "wild_type": "unchanged: P2 holds both PIE-1 and POS-1",
        "PIE-1": "E EPrecursor; terminal 34/20/0/54; P2 EMSPrecursor and ABp ABaPrecursor unchanged",
        "PAR-2": "unchanged: P2 holds neither PIE-1 nor POS-1, so E MSPrecursor as before",
        "PAR-3": "unchanged: POS-1 goes to P1 in the par-3 branch as PIE-1 does",
        "POS-1": "reported, not scored: P2 keeps PIE-1, so MOM-2 is presented and E is endoderm; Oldenbroek "
        "et al. 2013 report that every pos-1(zu148) embryo lacks endoderm, which this program does not "
        "model either before or after (the knockout of an absent factor is the wild type at baseline)",
        "du2014": "11 reproduced, 0 missed; pie-1's predicted-but-unobserved E-to-MS disappears",
        "verdict": "survives and moves pie-1: E is no longer lost, Intestine 0 -> 34; the published EXTRA "
        "intestine is still not reproduced (34 is the wild type's)",
    },
    "C1+C2": {"verdict": "rejected (C2's and C1's changes to MOM-2, MOM-5 and PAR-2)"},
    "C1+C3": {"verdict": "rejected: MOM-2 is ungated as in C1, so PAR-2 changes"},
    "C2+C3": {"verdict": "rejected (C2's changes to MOM-2, MOM-5 and PAR-2)"},
}

#: the published pie-1 phenotype as readings the program can be held to; a candidate must hold every
#: reading the baseline holds and at least one more
PIE1_READINGS = {
    "R1 P2 takes the EMS fate": "Mello et al. 1992; Du et al. 2014's pie-1 row",
    "R2 ABp takes the ABa fate": "Mango et al. 1994; Du et al. 2014's pie-1 row",
    "R3 E is not MS-like": "Du et al. 2014's pie-1 row lists no change of E; Mello et al. 1992's abstract "
    "counts ONE additional MS-like blastomere",
    "R4 no germ line (Reproductive 0)": "Mello et al. 1996",
    "R5 extra intestine (Intestine above the wild type's 34)": "Mello et al. 1992 as registered at 6bdecd4 "
    "(the abstract does not say it; the registered reading is kept, never loosened)",
}

REGISTRATION: dict[str, Any] = {
    "question": "Which sourced correction, if any, lets the worm program's pie-1 knockout keep E's intestine "
    "(or make extra) without changing the wild type or any other knockout?",
    "registered": "2026-10-10, lane-pie1, before any candidate was grown",
    "base": {"commit": BASE_COMMIT, "programs": [f"{ORG}/{NOLOOKUP}", f"{ORG}/{CONTACTS}"]},
    "sources": SOURCES,
    "candidates": CANDIDATES,
    "pairs": [list(p) for p in PAIRS],
    "excluded": EXCLUDED,
    "predictions": PREDICTIONS,
    "pie1_readings": PIE1_READINGS,
    "measured_per_arm": [
        "wild type: decided / right / wrong among the 555 at the means and in seeds 0-19, by founder; the "
        "eight founder types; the tree matched to 800 min; the Body's own compare",
        f"each knockout of {list(KNOCKOUTS)}: the types of {list(FOUNDER_NAMES)} (the last founder decision "
        "fired, as scripts/celegans_nolookup.py reads them) under the measured table",
        f"each of the eight registered knockouts {list(TERMINAL_KNOCKOUTS)}: Intestine / Muscle / "
        "Reproductive / decided among the 555, and the registered published reading",
        "Du et al. 2014: genomeos.organism.digital_development.compare on embryo_contacts.bio and on "
        "embryo_nolookup.bio, every observed change and every predicted-but-unobserved one",
    ],
    "falsifier": "a candidate is REJECTED if, against the baseline arm grown from the same base commit: "
    "(a) the wild type is not 56 / 56 / 0 at the means and in every seed, or the eight founder types are "
    "not the expected ones; (b) any knockout other than PIE-1 and POS-1 differs in any founder type, or "
    "any of the seven other registered knockouts differs in any terminal count; (c) any observed change of "
    "Du et al. 2014 that the baseline reproduces is missed, on either program; or (d) the arm is void. "
    "Whatever it does for pie-1",
    "apply_rule": "a candidate that is not rejected SURVIVES; it is applied only if it holds every pie-1 "
    "reading the baseline holds and at least one more. If several single candidates survive, C3 is "
    "preferred to C1 because C1 is contradicted by its own source's double depletion; any other tie is "
    "recorded and nothing is applied. A pair is applied only if no single candidate survives and exactly "
    "one pair does. Applying = making the registered edits verbatim to the files at the base commit; "
    "then mutants_nolookup.bio's pie1_terminal and the tests that pin PIE-1's terminal row are re-pinned "
    "to the measured values, with a dated comment",
    "void_if": [
        "V1 an edit's old text is not present exactly once when it is applied",
        "V2 the arm raises (a division naming collision, a parse error)",
        "V3 the wild type does not grow the reference tree to 800 min (missing cells or parent mismatches)",
        "V4 the Body's compare counts a different number of correct fates than the instrument",
        "V5 a decided wild-type cell is not attributed to tissue_E, tissue_D or tissue_germline",
    ],
    "wording": {
        "applied": "{c} is applied: {what}. With it the pie-1 knockout keeps E's intestine (Intestine "
        "{i}), loses the germ line, and P2 still takes the EMS fate and ABp the ABa fate; the published "
        "EXTRA intestine is still not reproduced. The wild type stays 56 / 56 / 0 and no other knockout "
        "changes. pie-1's E outcome follows from a rule written from Oldenbroek et al. 2013's measurement "
        "of mom-2 translation, not from the pie-1 phenotype; it is consistency with that measurement, not "
        "a prediction of the phenotype",
        "none": "No registered candidate survives; nothing is changed. The program's pie-1 still has no "
        "intestine",
    },
}


def registration_digest() -> str:
    return hashlib.sha256(json.dumps(REGISTRATION, sort_keys=True, ensure_ascii=True).encode()).hexdigest()


# ---- applying edits -----------------------------------------------------------------------------------


def base_files() -> dict[str, str]:
    """The worm's program directory at the base commit, file name -> text."""
    names = subprocess.run(
        ["git", "ls-tree", "--name-only", f"{BASE_COMMIT}:{ORG}"], capture_output=True, text=True, check=True
    ).stdout.split()
    return {
        n: subprocess.run(
            ["git", "show", f"{BASE_COMMIT}:{ORG}/{n}"], capture_output=True, text=True, check=True
        ).stdout
        for n in names
    }


def apply_edits(files: dict[str, str], arm: str) -> tuple[dict[str, str], list[str]]:
    """The files with the registered edits of an arm ('baseline', 'C1', 'C1+C3'), and any V1 failure."""
    out, void = dict(files), []
    if arm == "baseline":
        return out, void
    for cid in arm.split("+"):
        if cid not in CANDIDATES:
            raise SystemExit(f"{cid} is not a registered candidate: refused")
        for fname, old, new in CANDIDATES[cid]["edits"]:
            n = out[fname].count(old)
            if n != 1:
                void.append(f"V1 {cid} {fname}: {old[:60]!r} found {n} times")
                continue
            out[fname] = out[fname].replace(old, new)
    return out, void


def materialise(files: dict[str, str], where: Path) -> Path:
    org = where / "celegans"
    org.mkdir(parents=True)
    for n, text in files.items():
        (org / n).write_text(text)
    return org


# ---- one arm ------------------------------------------------------------------------------------------


def run_arm(arm: str, files: dict[str, str], ref: ReferenceLineage, cells, table) -> dict[str, Any]:  # noqa: ANN001
    edited, void = apply_edits(files, arm)
    out: dict[str, Any] = {"arm": arm, "void": void}
    if void:
        return out
    with tempfile.TemporaryDirectory() as tmp:
        org = materialise(edited, Path(tmp))
        try:
            module = parse_file(org / NOLOOKUP)
            targets = {d.id: d.to for d in module.decisions if d.action == "differentiate"}
            body = Body(module, means=True).run(until=HATCH)
            diff = compare(body, ref, until=HATCH)
            wt = decided(body, cells, ref, targets)
            ref_by_800 = sum(1 for c in ref.cells.values() if c.born <= HATCH)
            if diff.missing or diff.parent_mismatch or diff.matched != ref_by_800:
                void.append(f"V3 missing {len(diff.missing)}, parent mismatches {len(diff.parent_mismatch)}")
            if diff.fates_correct != wt["right"]:
                void.append(f"V4 compare {diff.fates_correct} against {wt['right']}")
            if wt["unattributed"] or set(wt["by_rule"]) - set(TISSUE_RULES):
                void.append(f"V5 unattributed {wt['unattributed'][:5]}, rules {wt['by_rule']}")
            founders = founder_types(body, targets, EXPECTED_FOUNDERS)
            seeds = {}
            for s in SEEDS:
                d = decided(Body(module, seed=s).run(until=HATCH), cells, ref, targets)
                seeds[s] = [d["decided"], d["right"], d["wrong"]]
            out["wild_type"] = {
                "decided": wt["decided"],
                "right": wt["right"],
                "wrong": wt["wrong"],
                "by_founder": wt["by_founder"],
                "by_rule": wt["by_rule"],
                "founders": founders,
                "founders_correct": sum(founders[n] == t for n, t in EXPECTED_FOUNDERS.items()),
                "seeds_all_56_56_0": all(v == [56, 56, 0] for v in seeds.values()),
                "seeds": seeds,
            }
            kos = {}
            for factor in KNOCKOUTS:
                b = Body(module, means=True, knockouts={factor}).run(until=HATCH)
                row: dict[str, Any] = {"founders": founder_types(b, targets, FOUNDER_NAMES)}
                if factor in TERMINAL_KNOCKOUTS:
                    d = decided(b, cells, ref, targets)
                    c = counts(d)
                    test = PUBLISHED_TEST.get(factor)
                    wanted = KNOCKOUT_FOUNDERS[factor]
                    row |= {
                        "terminal": c,
                        "terminal_by_founder": d["by_founder"],
                        "founder_phenotype_reproduced": all(
                            row["founders"].get(n) == t for n, t in wanted.items()
                        ),
                        "published_terminal_phenotype_reproduced": test(c, d["by_founder"]) if test else None,
                    }
                kos[factor] = row
            out["knockouts"] = kos
            out["du2014"] = {}
            for label, path in (("embryo_contacts", org / CONTACTS), ("embryo_nolookup", org / NOLOOKUP)):
                r = dd.compare(parse_file(path), table)
                out["du2014"][label] = {
                    "reproduced": r["reproduced"],
                    "missed": r["missed"],
                    "not_modelled": r["not_modelled"],
                    "observed": {
                        g["gene"]: {o["cell"]: o["status"] for o in g["observed"]} for g in r["genes"]
                    },
                    "predicted_not_observed": {
                        g["gene"]: sorted(
                            f"{p['cell']}:{p['mutant']}"
                            for p in g["predicted"]
                            if p["cell"] not in {o["cell"] for o in g["observed"]}
                        )
                        for g in r["genes"]
                    },
                }
        except Exception as e:  # noqa: BLE001 -- a void arm is reported, never dropped
            void.append(f"V2 {type(e).__name__}: {e}")
    return out


# ---- the verdicts -------------------------------------------------------------------------------------


def pie1_readings(arm: dict[str, Any]) -> dict[str, bool]:
    k = arm["knockouts"]["PIE-1"]
    f, t = k["founders"], k["terminal"]
    return {
        "R1 P2 takes the EMS fate": f.get("P2") == "EMSPrecursor",
        "R2 ABp takes the ABa fate": f.get("ABp") == "ABaPrecursor",
        "R3 E is not MS-like": f.get("E") != "MSPrecursor",
        "R4 no germ line (Reproductive 0)": t["Reproductive"] == 0,
        "R5 extra intestine (Intestine above the wild type's 34)": t["Intestine"] > WT_COUNTS["Intestine"],
    }


def breaks(arm: dict[str, Any], base: dict[str, Any]) -> list[str]:
    """The falsifier, clause by clause (REGISTRATION['falsifier'])."""
    if arm["void"]:
        return [f"(d) void: {v}" for v in arm["void"]]
    out = []
    w = arm["wild_type"]
    if (w["decided"], w["right"], w["wrong"]) != (56, 56, 0) or not w["seeds_all_56_56_0"]:
        out.append(f"(a) wild type {w['decided']}/{w['right']}/{w['wrong']}, seeds {w['seeds_all_56_56_0']}")
    if w["founders_correct"] != len(EXPECTED_FOUNDERS):
        out.append(f"(a) wild-type founders {w['founders_correct']} of {len(EXPECTED_FOUNDERS)}")
    for factor in KNOCKOUTS:
        if factor in ("PIE-1", *NOT_SCORED):
            continue
        a, b = arm["knockouts"][factor], base["knockouts"][factor]
        moved = {
            n: [b["founders"].get(n), a["founders"].get(n)]
            for n in FOUNDER_NAMES
            if a["founders"].get(n) != b["founders"].get(n)
        }
        if moved:
            out.append(f"(b) {factor} founders changed {moved}")
        if factor in TERMINAL_KNOCKOUTS and a["terminal"] != b["terminal"]:
            out.append(f"(b) {factor} terminal {b['terminal']} -> {a['terminal']}")
    for prog, d in arm["du2014"].items():
        for gene, cells in base["du2014"][prog]["observed"].items():
            for cell, status in cells.items():
                if status == "reproduced" and d["observed"][gene][cell] != "reproduced":
                    out.append(f"(c) Du et al. 2014 {gene} {cell} missed on {prog}")
    return out


def verdicts(arms: dict[str, dict[str, Any]]) -> dict[str, Any]:
    base = arms["baseline"]
    if base["void"]:
        raise SystemExit(f"the baseline arm is void: {base['void']}")
    base_held = pie1_readings(base)
    out = {}
    for name, arm in arms.items():
        if name == "baseline":
            continue
        b = breaks(arm, base)
        readings = pie1_readings(arm) if not arm["void"] else {}
        moves = (
            bool(readings)
            and all(readings[r] for r, v in base_held.items() if v)
            and sum(readings.values()) > sum(base_held.values())
        )
        out[name] = {"breaks": b, "survives": not b, "pie1_readings": readings, "moves_pie1": moves}
    singles = [c for c in CANDIDATES if out[c]["survives"] and out[c]["moves_pie1"]]
    pairs = ["+".join(p) for p in PAIRS if out["+".join(p)]["survives"] and out["+".join(p)]["moves_pie1"]]
    if len(singles) == 1:
        applied = singles[0]
    elif set(singles) == {"C1", "C3"}:
        applied = "C3"
    elif not singles and len(pairs) == 1:
        applied = pairs[0]
    else:
        applied = None
    return {
        "arms": out,
        "survivors": singles,
        "surviving_pairs": pairs,
        "apply": applied,
        "base_pie1": base_held,
    }


# ---- the run ------------------------------------------------------------------------------------------


def committed(path: Path) -> bool:
    tracked = subprocess.run(["git", "ls-files", "--error-unmatch", str(path)], capture_output=True)
    clean = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", str(path)])
    return tracked.returncode == 0 and clean.returncode == 0


def run() -> int:
    me = Path(__file__).resolve().relative_to(ROOT)
    if not committed(me):
        raise SystemExit(f"{me} is not committed and clean: the registration is the committed file; refused")
    files = base_files()
    entries = pv.inputs(atlas=False, extra=[TABLE])
    entries += [
        pv.bytes_entry(f"git {BASE_COMMIT[:12]}:{ORG}/{n}", files[n].encode())
        for n in sorted(files)
        if n.endswith(".bio") or n.endswith(".tsv")
    ]
    ref = ReferenceLineage.load()
    cells = terminal_cells(ref)
    table = json.loads(TABLE.read_text())["table"]
    arms = {}
    for arm in ["baseline", *CANDIDATES, *("+".join(p) for p in PAIRS)]:
        arms[arm] = run_arm(arm, files, ref, cells, table)
        print(arm, "void" if arms[arm]["void"] else "ran", arms[arm]["void"])
    v = verdicts(arms)
    applied = v["apply"]
    if applied:
        k = arms[applied]["knockouts"]["PIE-1"]["terminal"]
        what = CANDIDATES[applied]["what"] if applied in CANDIDATES else applied
        headline = REGISTRATION["wording"]["applied"].format(c=applied, what=what, i=k["Intestine"])
    else:
        headline = REGISTRATION["wording"]["none"]
    out = {
        "result": RESULT_NAME,
        "date": date.today().isoformat(),
        "registration": str(me),
        "registration_digest": registration_digest(),
        "base_commit": BASE_COMMIT,
        "base_sha256": {
            n: hashlib.sha256(files[n].encode()).hexdigest() for n in (FOUNDERS, NOLOOKUP, CONTACTS)
        },
        "headline": headline,
        "verdicts": v,
        "arms": arms,
    }
    m = pv.manifest(
        entries,
        {"until_min": HATCH, "seeds": SEEDS, "timers": "means", "base_commit": BASE_COMMIT},
        exclusions=[{k: e["what"] for k, e in EXCLUDED.items()}],
        partitions=pv.NO_SPLIT,
        atlas=False,
    )
    m["sources"] += [{"accession": s["cite"], "version": s["fetched"]} for s in SOURCES.values()]
    m["code_cleanliness"] = mf.code_cleanliness(__file__, OWN_CODE)
    print("saved", save_result(RESULT_NAME, out, manifest=m))
    summary = {
        "headline": headline,
        "survivors": v["survivors"],
        "apply": applied,
        "breaks": {k: a["breaks"] for k, a in v["arms"].items()},
    }
    print(json.dumps(summary, indent=1))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--digest", action="store_true", help="print the registration's digest")
    a = ap.parse_args()
    if a.digest:
        print(registration_digest())
        return 0
    return run()


if __name__ == "__main__":
    raise SystemExit(main())
