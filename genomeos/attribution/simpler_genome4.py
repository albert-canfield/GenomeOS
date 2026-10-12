# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""simpler_genome4: is the worm's between-gene lead THE WORM, or is it operons and trans-splicing?

THE THREE-SPECIES PICTURE, carried verbatim from the records that registered it and not restated.

* **Worm** (`595d7c2` registered, docs/SIMPLER-GENOMES.md): ahead of nearest TSS by **+3.73 points
  [+2.79, +4.74]** on the 2,335 elements that lie BETWEEN genes, where no host gene can be named --
  a lead that **survives** the worm's own outron and trans-splicing caveat at **+4.47 [+3.16,
  +5.78]** on the 1,745 survivors of that exclusion.
* **Fly** (`ff063e6` registered): between genes **+0.12 [-0.25, +0.61]** on Janssens et al. 2022's
  3,265 brain elements.
* **Mouse** (docs/SIMPLER-GENOMES.md section 12): **+0.26 points [-0.01, +0.53]**,
  `NO_DIFFERENCE_DETECTED`, on a population powered in advance to about **0.74 points, eleven times
  the fly's differing**, with the primary interval's lower edge at -0.01, touching zero.

So the worm is the outlier of three, and its between-gene lead reads as A PROPERTY OF THE WORM ON
THESE SOURCES, not of the domain rule. That wording is the supervisor's constraint and is carried
word for word: **on these sources, and nothing wider.**

THE QUESTION THIS LANE ASKS, and it is the only one it asks: the worm is not merely a worm. It has
**operons and SL trans-splicing**, which the fly and the mouse do not. So "the worm" and "operons
and trans-splicing" are two different explanations of the same outlier, and the three genomes
measured so far cannot separate them, because the only genome with the mechanism is also the only
genome with the lead.

WHAT WOULD SEPARATE THEM, and the two candidate families answer DIFFERENT questions. This is the
design point and it is registered before any source is read:

* **another nematode** (C. briggsae, C. remanei, Pristionchus pacificus) keeps the mechanism AND
  the lineage, so it separates **C. elegans specifically** from **nematodes generally**. If the
  lead reappears, "the worm" narrows to "nematodes"; it does not reach the mechanism.
* **a non-nematode genome with operons and trans-splicing** -- a tunicate such as Ciona, or a
  platyhelminth -- keeps the mechanism and DROPS the lineage, so it separates **operons and
  trans-splicing as a mechanism** from **the nematode lineage**. This is the sharper test of the
  question as asked, and it is the one a reader of the three-species picture would want.

Neither is run here. THIS LANE IS AN INVENTORY and nothing else: it asks whether any candidate has
a source that qualifies, by criteria fixed below before a single source is opened.

**NO_SOURCE_QUALIFIES IS A COMPLETE ANSWER**, and for a fourth genome it is the likely one. The
worm's source (Serizay et al. 2020) and the fly's (Janssens et al. 2022) are unusual: tissue- or
cell-type-resolved accessibility AND matched RNA in ONE source, with a usable gene annotation. For
a non-model operon-bearing genome that combination is rare, and reporting its absence is the
lane's result rather than its failure.

**THE STOP LIST OF 2026-09-27 WAS READ, NOT IGNORED**, and this lane touches none of its six items.
Its clauses are: nothing more on the node model until wave 0's measured test lands (met since
2026-09-27); **no more biosamples**; nothing is built in area J; no AlphaGenome requests on the 882
blocks; area E's fate ceiling is not chased; BioForge's confidence is not calibrated. The
biosample clause's subject is the **human reader's 13 biosamples** -- `READER_TERMS` in
`genomeos/genome/context_evidence.py` -- and this lane adds none, names none and widens none, which
is the reading `lane-flysource` and `lane-mousedomain` both recorded for their own sources. A
fourth genome's tissue source is not a human reader biosample. This lane makes no model request,
spends nothing, and reads no person's genome.
"""

from __future__ import annotations

from genomeos.attribution import simpler

#: Carried from the shared module rather than restated, so it cannot drift from the worm's and the
#: fly's and the mouse's. Same stop rule, same number, one definition.
MIN_DISCORDANT = simpler.MIN_DISCORDANT

#: The question, in one sentence, so a later reader cannot widen it.
QUESTION = (
    "is the worm's between-gene lead a property of THE WORM, or of OPERONS AND TRANS-SPLICING?"
    " The three genomes measured so far cannot separate those, because the only genome with the"
    " mechanism is also the only genome with the lead"
)

#: What each candidate family can and cannot settle. Registered before any source is opened.
CANDIDATE_FAMILIES = (
    {
        "family": "another nematode",
        "examples": ("C. briggsae", "C. remanei", "Pristionchus pacificus"),
        "keeps": "the mechanism AND the lineage",
        "separates": "C. elegans specifically from nematodes generally",
        "cannot_settle": "the mechanism. If the lead reappears, 'the worm' narrows to 'nematodes'"
        " and the operon explanation is untouched",
    },
    {
        "family": "a non-nematode genome with operons and trans-splicing",
        "examples": ("a tunicate such as Ciona", "a platyhelminth"),
        "keeps": "the mechanism, and DROPS the lineage",
        "separates": "operons and trans-splicing as a mechanism from the nematode lineage",
        "cannot_settle": "whether some third property of that genome carries the lead instead;"
        " a single fourth genome cannot, and the reading must say so",
    },
)

#: The qualification criteria. FIXED HERE BEFORE ANY SOURCE IS OPENED, so a source cannot be
#: admitted by a criterion written once it was seen.
QUALIFIES = (
    "a candidate qualifies only if ALL FOUR hold: (1) ACCESSIBILITY resolved by tissue or cell"
    " type, from a named source; (2) RNA for THE SAME tissues or cell types, FROM THE SAME SOURCE"
    " -- the worm's and the fly's both were, and two sources would introduce a between-study"
    " difference the endpoint cannot separate from the organism; (3) a USABLE GENE ANNOTATION with"
    " positions and strands for protein-coding genes, which is what the domain rule reads; and (4)"
    " enough elements that at least MIN_DISCORDANT of them get DIFFERENT genes under the two rules"
    " -- the same 200 the worm, the fly and the mouse were held to"
)

#: And what does not qualify, named so the inventory cannot quietly relax.
DISQUALIFIES = (
    "accessibility and RNA from DIFFERENT studies, however well matched; whole-animal or"
    " whole-embryo RNA with no tissue resolution, because the endpoint compares a gene's"
    " expression ACROSS tissues; an annotation without strands, since the domain rule needs the"
    " face of the 5' end; a genome with operons but NO demonstrated SL trans-splicing, or the"
    " reverse, because the hypothesis names both; and any candidate whose element count cannot"
    " reach 200 discordant, which is UNDERPOWERED and is reported as such rather than scored"
)

#: The complete answer, stated in advance so it is not read later as a failure to find one.
NO_SOURCE_QUALIFIES_IS_COMPLETE = (
    "if no candidate meets all four criteria the lane reports NO_SOURCE_QUALIFIES and stops. That"
    " is a COMPLETE ANSWER and the likely one: the worm's and the fly's sources are unusual in"
    " carrying tissue-resolved accessibility and matched RNA together, and for a non-model"
    " operon-bearing genome that combination is rare. Reporting its absence tells the next reader"
    " that the mechanism-versus-lineage question is NOT CURRENTLY DECIDABLE on public data, which"
    " is worth more than a strained population that cannot carry the endpoint"
)

#: The wording constraint, carried word for word from the supervisor.
WORDING_CONSTRAINT = (
    "say 'on these sources', nothing wider. The three-species picture licenses 'the worm is the"
    " outlier of three on these sources', and it does not license any statement about the domain"
    " rule in general, about nematodes as a group, or about operons, until a fourth genome is"
    " measured"
)

#: What an inventory cannot claim, however it comes out.
DOES_NOT_CLAIM = (
    "not that a qualifying source would settle the question -- one fourth genome separates two"
    " explanations and leaves every other difference between it and C. elegans uncontrolled; not"
    " that NO_SOURCE_QUALIFIES means no such data exists, only that none was found by these"
    " criteria on the sources checked, each named; and nothing at all about the worm's own"
    " reading, which stands as registered at +4.47 [+3.16, +5.78] on its 1,745 survivors"
)
