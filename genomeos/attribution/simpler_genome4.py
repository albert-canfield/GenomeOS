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


# ----------------------- the limit, as a LIMIT and not as prose -- registered before any source
# The coordinator asked for this in the shape of globcensus's
# `test_PLANTED_a_group_constant_at_a_shrunken_value_is_INVISIBLE_here`: a boundary a later summary
# cannot drop, because a test fails if the result does not carry it.

#: The key every result of this lane must carry, whatever it finds.
LIMIT_KEY = "what_this_cannot_reach"

#: And its text, fixed here, so the result cannot restate it more weakly.
SECOND_NEMATODE_LIMIT = (
    "IF THE QUALIFYING SOURCE IS A SECOND NEMATODE, THE READING SAYS THE ANSWER NARROWS 'THE"
    " WORM' AND DOES NOT REACH 'OPERONS AND TRANS-SPLICING'. A nematode keeps both the mechanism"
    " and the lineage, so a lead that reappears there is consistent with either explanation and"
    " separates neither. Only a non-nematode genome with operons and trans-splicing drops the"
    " lineage while keeping the mechanism, and only that comparison reaches the question as asked"
)

#: The archives to be searched, the number of queries, and the species. DECLARED BEFORE SEARCHING,
#: because a negative's boundary is the negative: "no source qualifies" means nothing without them.
SEARCH_PROTOCOL = {
    "species": (
        "Caenorhabditis briggsae",
        "Caenorhabditis remanei",
        "Pristionchus pacificus",
        "Ciona robusta / Ciona intestinalis",
        "Schistosoma mansoni",
    ),
    "archives": (
        "GEO / SRA",
        "ArrayExpress / ENA",
        "ENCODE",
        "WormBase and WormBase ParaSite",
        "ANISEED (the Ciona-specific database)",
        "Ensembl Metazoa (annotation only)",
    ),
    "queries_per_species": "at most 4, and the count actually run is recorded per species",
    "what_a_search_can_and_cannot_show": "a search can show that a source EXISTS and meets the"
    " four criteria; it cannot show that none exists. So a negative is reported as 'none found by"
    " these queries in these archives for these species', never as 'no such data exists'",
}

#: The rule for carrying this lane's negative upward. The coordinator's standing correction,
#: recorded here because this lane's negative is the kind that gets widened.
NEGATIVE_CARRIES_ITS_BOUNDARY = (
    "when this lane's negative is carried upward it travels WITH ITS OWN BOUNDARY IN THE SAME"
    " SENTENCE: which archives were searched, how many queries were run, and which species. A"
    " negative without those three is a different and much larger claim than the one measured,"
    " and the coordinator records having inflated a lane's negative three times in two days"
)


def result_carries_its_limit(result: dict) -> bool:
    """Whether a result of this lane carries :data:`LIMIT_KEY` with the limit's own wording.

    The point of the function is that a test can call it: a result that drops the limit, or
    softens it, fails rather than passing quietly with a weaker sentence.
    """
    text = result.get(LIMIT_KEY)
    if not isinstance(text, str):
        return False
    return "DOES NOT REACH 'OPERONS AND TRANS-SPLICING'" in text and "separates neither" in text


# ============================ THE INVENTORY'S ANSWER: NO_SOURCE_QUALIFIES -- 2026-10-12
# Five species, six archives, nine queries. The answer is the one declared likely in advance, and
# it is a complete answer. NOTE THERE IS NO RESULT JSON: nothing was measured, no RNA value was
# read, and a result record would imply a measurement. The inventory is a SOURCE ASSESSMENT and
# lives here, in committed code, which also keeps it clear of the ordering trap this project
# recorded today -- a test committed before the record it reads reddens every sha in between.

#: The answer.
INVENTORY_ANSWER = "NO_SOURCE_QUALIFIES"

#: Per species: the queries actually run, what was found, and WHICH criterion it fails. The
#: criterion is named because "no source" without it is unfalsifiable.
INVENTORY_BY_SPECIES = (
    {
        "species": "Caenorhabditis briggsae",
        "queries": 1,
        "found": "nothing tissue-resolved. The searches return C. elegans work instead",
        "fails": "criterion (1), accessibility resolved by tissue or cell type",
    },
    {
        "species": "Caenorhabditis remanei",
        "queries": 1,
        "found": "nothing. Same searches, same C. elegans results",
        "fails": "criterion (1)",
    },
    {
        "species": "Pristionchus pacificus",
        "queries": 1,
        "found": "ATAC-seq EXISTS, as one input to a ChromHMM chromatin-state annotation beside"
        " histone modifications (Genome Research gr.242362.118, 2019). But it is whole-animal"
        " material for state annotation, not a tissue or cell-type panel, and no matched"
        " tissue RNA from the same source was found",
        "fails": "criteria (1) and (2)",
    },
    {
        "species": "Ciona robusta / intestinalis / savignyi",
        "queries": 2,
        "found": "BOTH assays exist and are integrated, which is why this was the most promising"
        " candidate -- but not by tissue. C. robusta: an accessome of 56,090 non-overlapping"
        " accessible regions covering 9.25% of the genome, from FACS-purified cells, resolved by"
        " CARDIOPHARYNGEAL FATE CHOICE -- one lineage, not a tissue panel. C. savignyi: ATAC-seq"
        " and RNA-seq jointly across LARVAL DEVELOPMENTAL STAGES. A stage is not a tissue, and"
        " the endpoint compares a gene's expression ACROSS TISSUES",
        "fails": "criterion (1), and so criterion (2) with it",
    },
    {
        "species": "Schistosoma mansoni",
        "queries": 1,
        "found": "rich single-cell RNA atlases (a miracidium whole-body atlas of about 365 cells"
        " in 19 transcriptionally distinct clusters; a separate intra-mammalian stage atlas), and"
        " an ATAC-seq PROJECT at Aberystwyth running October 2021 to July 2022 -- but no"
        " published tissue-resolved accessibility matched to RNA in one source was found",
        "fails": "criterion (2), RNA and accessibility from the same source",
    },
)

#: The negative WITH ITS BOUNDARY IN THE SAME SENTENCE, as the coordinator's standing correction
#: requires. This is the sentence to carry upward; it is not to be shortened.
INVENTORY_NEGATIVE = (
    "NO SOURCE QUALIFIES among five species -- C. briggsae, C. remanei, Pristionchus pacificus,"
    " Ciona robusta/intestinalis/savignyi and Schistosoma mansoni -- on SIX queries across"
    " GEO/SRA, ArrayExpress/ENA, ENCODE, WormBase and WormBase ParaSite, ANISEED and Ensembl"
    " Metazoa. THAT IS NOT 'no such data exists': it is none found by these queries in these"
    " archives for these species, and a deeper search of any one of them could overturn it"
)

#: What the inventory therefore says about the question.
INVENTORY_READING = (
    "the mechanism-versus-lineage question is NOT CURRENTLY DECIDABLE on the public data found"
    " here. The worm's lead stands as registered, and whether it belongs to the worm or to"
    " operons and trans-splicing cannot be separated on these sources. The failure is almost"
    " always the SAME criterion: accessibility resolved by TISSUE with RNA for the same tissues"
    " from the SAME source. Ciona has both assays and integrates them, and still fails, because"
    " its resolution is a lineage or a developmental stage rather than a tissue panel"
)

#: And the by-product, which is worth more than the negative and is NOT this lane's to run.
WITHIN_WORM_REPLICATION_INSTEAD = (
    "THE SEARCHES FOUND SOMETHING BETTER THAN WHAT THEY WERE LOOKING FOR. Durham et al. 2021"
    " (Genome Research 31:1952, gr.271791.120) is sci-ATAC-seq in L2 C. elegans built"
    " deliberately TO MATCH available sci-RNA-seq, tissue by tissue -- an INDEPENDENT C. elegans"
    " source from Serizay et al. 2020, which is the worm's registered source. That does not"
    " answer the mechanism question, because it is the same species and so changes neither the"
    " mechanism nor the lineage. It answers a DIFFERENT and currently unasked one: is the worm's"
    " +4.47 a property of THE WORM or of SERIZAY 2020? A source replication inside the worm is"
    " cheaper than a fourth genome, needs no new organism, and is the only one of these routes"
    " for which qualifying data was actually found. Proposed, not run: it needs its own"
    " registration, and this lane's registration is for a fourth genome"
)


def stated_query_count() -> int:
    """The query count in :data:`INVENTORY_NEGATIVE`, parsed back out of the sentence.

    It exists so a test can compare the sentence against the per-species records. The first
    version of that sentence said NINE while the records summed to SIX -- a number written rather
    than counted, which is the error this whole day was about, and the records caught it. A
    boundary that can drift from its own evidence is not a boundary.
    """
    import re

    words = {
        "one": 1,
        "two": 2,
        "three": 3,
        "four": 4,
        "five": 5,
        "six": 6,
        "seven": 7,
        "eight": 8,
        "nine": 9,
        "ten": 10,
    }
    match = re.search(r"on ([A-Za-z]+) queries", INVENTORY_NEGATIVE)
    if match is None:
        raise ValueError("INVENTORY_NEGATIVE no longer states a query count")
    return words[match.group(1).lower()]
