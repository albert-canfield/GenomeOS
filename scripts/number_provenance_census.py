# SPDX-License-Identifier: AGPL-3.0-or-later
"""Count what every number in the hand-authored BioLang corpus rests on
(data/results/number_provenance_census.json).

    uv run --frozen python scripts/number_provenance_census.py

Pre-registered in data/results/number_provenance_registration.json. The population, the scoped
fields, the classes, the cascade, the convention tags, the citation-furniture strip and the
denominator all come from `genomeos.lang.number_provenance` and were committed before this count.
This script adds only what a count produces: the tallies, the per-rule table for the module the
lane was asked about, the hand adjudication of every member of the sourced class, and the record of
which cited sources were fetched and what each one does NOT claim.

It changes no number. It writes no .bio file. It makes no model request and spends no money.
"""

from __future__ import annotations

import collections
import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.lang import number_provenance as np  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "number_provenance_census"
REGISTRATION = "number_provenance_registration"

OWN_CODE = (
    "genomeos/lang/number_provenance.py",
    "scripts/number_provenance_register.py",
    "scripts/number_provenance_census.py",
    "tests/test_number_provenance.py",
)

#: The lane was asked about this module in particular, so its rules are tabulated one by one.
FOCUS = "data/demo/gastrulation.bio"

# --- the hand adjudication of the sourced class --------------------------------------------------
# The registration fixed the rule: a member of S_quoted is demoted when the matched digits fall in
# citation furniture the strip missed. Applying that rule to all 30 members demotes NONE of them.
# Two members would be demoted under a WIDER rule, and widening a test after seeing the data is the
# error this lane exists to avoid, so they keep their registered class and are recorded here as
# contested instead. The registered figure and the contested pair are both reported.
CONTESTED_S_QUOTED = (
    {
        "file": "data/demo/lateral_inhibition.bio",
        "line": 39,
        "declaration": "Np inhibits D",
        "field": "strength",
        "value": 1.0,
        "matched": "the `1` in `g(n) = 1 / (1 + 100 n^2)`",
        "contest": "the matched 1 is the numerator of a Hill function and not, as written, a "
        "statement about the strength of the interaction, so the match looked accidental.",
        "the_contest_is_withdrawn": "checking the engine's arithmetic refutes it. The inhibiting "
        "branch is `rep *= 1 - strength * x^h/(threshold^h + x^h)`, which equals the quoted "
        "`1 / (1 + 100 n^2)` only when strength is exactly 1.0: any smaller strength leaves a "
        "floor of 1 - strength. So the quoted formula does constrain this slot's strength to 1.0, "
        "and the match is not accidental. S_quoted stands on its merits and not only on the "
        "registered rule.",
        "registered_class": "S_quoted",
        "class_if_the_rule_were_widened": "S_quoted; the widening does not reach it",
    },
    {
        "file": "genomeos/std/human_stages.bio",
        "line": 33,
        "declaration": "human_cleavage",
        "field": "duration",
        "value": 20.0,
        "matched": "the `20` in `about one division per 20 h`",
        "contest": "the 20 h is the file's own arithmetic from two cited milestones (2-cell at "
        "~24 h, 8-cell at day 3). The source states the milestones; the 20 is the author's "
        "division of them, written inside a curated evidence string.",
        "registered_class": "S_quoted",
        "class_if_the_rule_were_widened": "E_existence_only",
    },
)

#: The strongest case in the corpus that the registered S_derived test does NOT admit, reported in
#: full with its arithmetic so a reader can judge it. `DERIVATIONS` stays empty and S_derived is 0:
#: the transform below runs through the engine's own Hill form in genomeos/runtime/grn.py and not
#: through anything the evidence string says, and the registered test required the latter. On top of
#: that, Collier et al. 1996 is not open access, so the parameter values the evidence string
#: attributes to it could not be fetched and confirmed at all.
STRONGEST_NON_S = (
    {
        "file": "data/demo/lateral_inhibition.bio",
        "line": 37,
        "declaration": "Dext activates N",
        "field": "threshold",
        "value": 0.1,
        "evidence": 'curated "Collier et al. 1996: f(d) = d^2 / (0.01 + d^2)"',
        "arithmetic": "the engine's activating Hill form is x^h / (threshold^h + x^h), so "
        "threshold 0.1 with hill 2 gives d^2 / (0.1^2 + d^2) = d^2 / (0.01 + d^2), the quoted "
        "formula exactly. 0.1 = sqrt(0.01).",
        "registered_class": "E_existence_only",
        "why_not_s_derived": "the transform is in genomeos/runtime/grn.py, not in the evidence "
        "string, and the registered S_derived test required the transform to be written in the "
        "string. Promoting it would be widening a test after seeing the data.",
    },
    {
        "file": "data/demo/lateral_inhibition.bio",
        "line": 39,
        "declaration": "Np inhibits D",
        "field": "threshold",
        "value": 0.1,
        "evidence": 'curated "Collier et al. 1996: g(n) = 1 / (1 + 100 n^2)"',
        "arithmetic": "the engine's inhibiting Hill form is threshold^h / (threshold^h + x^h) = "
        "1 / (1 + (x/threshold)^h), so threshold 0.1 with hill 2 gives 1 / (1 + 100 n^2), the "
        "quoted formula exactly. 0.1 = 1/sqrt(100).",
        "registered_class": "E_existence_only",
        "why_not_s_derived": "the same reason.",
    },
)

# --- what was fetched, and what each source does NOT claim ---------------------------------------
# A citation that supports a rule's EXISTENCE is not a citation for its STRENGTH. These are the
# checks that back that sentence for this corpus. Abstracts were fetched through the Europe PMC REST
# service (resultType=core) and one full text through PubMed Central; where a paper is not open
# access the full text could not be fetched and the entry says so rather than implying a check.
SOURCE_CHECKS = (
    {
        "cited_as": "Conlon et al. 1994, Development 120:1919",
        "cited_for": "data/demo/gastrulation.bio line 18, `rule Nodal activates TBXT`, and its "
        "strength 1.0, threshold 1.0 and hill 2",
        "fetched": "abstract, Europe PMC REST resultType=core, pmid 7924997",
        "full_text_fetched": False,
        "full_text_why_not": "not open access (Europe PMC isOpenAccess=N, no PMCID); the "
        "fullTextXML endpoint returns nothing.",
        "what_it_states": "A primary requirement for nodal in the formation and maintenance of the "
        "primitive streak in the mouse. The 413.d insertional proviral insertion is a loss of "
        "function mutation; nodal RNA is detected from day 5.5; 413.d mutant embryos show no "
        "morphological evidence for a primitive streak.",
        "what_it_does_not_state": "no dose-response, no half-maximal concentration, no Hill "
        "coefficient and no numeric strength for the effect of nodal on Brachyury/T. The only "
        "number in the abstract is a penetrance, `about 25% of mutant embryos do form randomly "
        "positioned patches of cells of a posterior mesodermal character`, which is a fraction of "
        "embryos and not a parameter of a regulatory function.",
    },
    {
        "cited_as": "Kanai-Azuma et al. 2002, Development 129:2367",
        "cited_for": "data/demo/gastrulation.bio line 19, `rule Nodal activates SOX17`, and its "
        "strength 1.0, threshold 3.0 and hill 4",
        "fetched": "abstract, Europe PMC REST resultType=core, pmid 11973269",
        "full_text_fetched": False,
        "full_text_why_not": "not open access (isOpenAccess=N, no PMCID).",
        "what_it_states": "Depletion of definitive gut endoderm in Sox17-null mutant mice: "
        "Sox17(-/-) embryos are deficient of gut endoderm; reduced occupancy of the definitive "
        "endoderm; elevated apoptosis in the foregut; Sox17-null ES cells are excluded from mid- "
        "and hindgut endoderm in chimeras.",
        "what_it_does_not_state": "it is a loss-of-function study of Sox17, not a measurement of "
        "nodal dose. It states no threshold at which NODAL activates SOX17, no Hill coefficient and "
        "no strength. It does not even address the direction `Nodal activates SOX17` "
        "quantitatively; it establishes that Sox17 is required for endoderm.",
    },
    {
        "cited_as": "Thomson et al. 2011, Cell 145:875",
        "cited_for": "data/demo/gastrulation.bio lines 21 and 22, `rule Tbxt inhibits SOX2` "
        "(strength 1.0, threshold 2.0, hill 3) and `rule Sox2 inhibits TBXT` (strength 1.0, "
        "threshold 4.0, hill 3)",
        "fetched": "abstract, Europe PMC REST resultType=core, pmid 21663792",
        "full_text_fetched": False,
        "full_text_why_not": "PMCID PMC5603300 exists but isOpenAccess=N; the fullTextXML endpoint "
        "returns nothing.",
        "what_it_states": "Pluripotency factors in embryonic stem cells regulate differentiation "
        "into germ layers: Oct4 suppresses neural ectodermal and promotes mesendodermal "
        "differentiation; Sox2 inhibits mesendodermal and promotes neural ectodermal "
        "differentiation; differentiation signals `continuously and asymmetrically modulate Oct4 "
        "and Sox2 protein levels`. This supports the EXISTENCE and the SIGN of a Sox2 brake on the "
        "mesendodermal programme.",
        "what_it_does_not_state": "no threshold, no Hill coefficient, no strength. The abstract "
        "states that the regulation is level-dependent and gives no level. The two rules carry "
        "DIFFERENT thresholds, 2.0 and 4.0, from the same citation, and nothing in the citation "
        "distinguishes them.",
    },
    {
        "cited_as": "Lolas et al. 2014, PNAS 111:4478",
        "cited_for": "data/demo/gastrulation.bio line 23, `rule Sox17 inhibits TBXT` (Fig. 3C) and "
        "line 27, `rule Tbxt activates SOX17` (Fig. 1A ChIP-seq, Fig. 3A knockdown)",
        "fetched": "abstract, Europe PMC REST resultType=core, pmid 24616493",
        "full_text_fetched": False,
        "full_text_why_not": "PMCID PMC3970479 exists but isOpenAccess=N; the fullTextXML endpoint "
        "returns nothing, so Fig. 3C itself was not read here.",
        "what_it_states": "Charting Brachyury-mediated developmental pathways: ChIP-seq and "
        "ChIP-exo localisation of Brachyury; `Brachyury functions primarily as a transcriptional "
        "activator genome-wide`; `an unexpected gene-regulatory feedback loop consisting of "
        "Brachyury, Foxa2, and Sox17`. This supports the EXISTENCE of both rules and the sign of "
        "`Tbxt activates SOX17`, which is the `5cbce26` correction, and that correction stands.",
        "what_it_does_not_state": "no dose, no strength, no threshold and no Hill coefficient. "
        "data/demo/gastrulation.bio says so itself in the comment at line 28: `Lolas 2014 supports "
        "the sign of the two Tbxt/Sox17 rules above and reports no dose, strength, threshold or "
        "Hill coefficient`. The census agrees with that comment and classifies all six numbers of "
        "those two rules E_existence_only.",
        "a_false_match_this_citation_would_have_caused": "the evidence string contains `Fig. 3C`. "
        "Without the citation-furniture strip the matcher would have read that 3 as a quoted number "
        "and called `hill: 3` sourced. tests/test_number_provenance.py pins the demotion.",
    },
    {
        "cited_as": "Zhou et al. 2018, Nat Genet 50:591",
        "cited_for": "genomeos/std/methylation.bio, `param methylation.neighbour_window_bp = 35 "
        "bp`, and (in data/demo/methylation_erosion.bio) the mechanism behind an unsourced "
        "per-division fidelity",
        "fetched": "abstract (Europe PMC, pmid 29610480) and full text (PubMed Central "
        "PMC5893360, open access)",
        "full_text_fetched": True,
        "what_it_states": "the solo-WCGW context is defined by `CpGs with the combination of zero "
        'neighboring CpGs ("solo") and the WCGW motif`, solo meaning no further CpG in a 35 bp '
        "window on either side. The 35 is the paper's own definition, so "
        "`methylation.neighbour_window_bp = 35 bp` is a number this census confirms as sourced.",
        "what_it_does_not_state": "no per-division maintenance fidelity as a number. Asked "
        "directly, the full text gives none: it describes methylation loss accumulating over "
        "divisions and a re-methylation-window model without quantifying fidelity per division. "
        "data/demo/methylation_erosion.bio already says this - `Zhou et al. 2018 show the loss, not "
        "the rate` - and the census classifies its 0.90 A_asserted_names_a_work, agreeing with the "
        "file's own note.",
    },
    {
        "cited_as": "Matsuda et al. 2020, Science 369:1450",
        "cited_for": "data/demo/segmentation_clock.bio `param period_h = 5.0 h` and "
        "genomeos/std/human_stages.bio `timer human_segmentation_clock duration: 5 h`",
        "fetched": "abstract, Europe PMC REST resultType=core, pmid 32943519",
        "full_text_fetched": False,
        "full_text_why_not": "not open access (isOpenAccess=N, no PMCID).",
        "what_it_states": "`we recapitulate murine and human segmentation clocks that display 2- to "
        "3-hour and 5- to 6-hour oscillation periods, respectively`. The source states a number for "
        "the quantity the field names, so both slots are sourced.",
        "what_it_does_not_state": "the source states a RANGE, 5 to 6 h, and both files write its "
        "lower bound, 5. The census counts them S_quoted and records here that the written value is "
        "the bottom of the stated interval and not a point estimate the paper gives.",
    },
    {
        "cited_as": "Sender & Milo 2021, Nat Med 27:45",
        "cited_for": "genomeos/std/human_turnover.bio, all 16 `timer *_lifespan` durations, which "
        "are 16 of the 30 members of the sourced class",
        "fetched": "abstract, Europe PMC REST resultType=core, pmid 33432173",
        "full_text_fetched": False,
        "full_text_why_not": "not open access (isOpenAccess=N, no PMCID). The per-cell-type "
        "lifespans live in the paper's table, which was not fetched.",
        "what_it_states": "`a total cellular mass turnover of 80 +/- 20 grams per day` and `close "
        "to 90% of the (0.33 +/- 0.02) x 10^12 cells per day turnover was blood cells`. The paper "
        "is explicitly a quantitative integration of `ubiquity, mass and lifespan of all major cell "
        "types`, so it is the right KIND of source for a per-type lifespan.",
        "what_it_does_not_state": "the abstract states none of the 16 values the file attributes to "
        "it. Those 16 are the largest block of the sourced class and this census did NOT verify "
        "any of them against the paper. `S_quoted` for them means only that the FILE quotes a "
        "number equal to the written value; a lane with access to the table should check them.",
    },
    {
        "cited_as": "Collier et al. 1996, J Theor Biol 183:429",
        "cited_for": "data/demo/lateral_inhibition.bio, the two rules' strength, threshold and hill "
        "and the two proteins' initial values",
        "fetched": "abstract, Europe PMC REST resultType=core, pmid 9015458",
        "full_text_fetched": False,
        "full_text_why_not": "not open access (isOpenAccess=N, no PMCID).",
        "what_it_states": "`we construct and analyse a simple and general mathematical model of "
        "such contact-mediated lateral inhibition`, with the postulate that receipt of inhibition "
        "diminishes the ability to deliver it. It is the one cited work in this corpus whose "
        "subject is a parameterised model, which is why the file can quote formulae at all.",
        "what_it_does_not_state": "the abstract states no parameter value: not the 0.01, not the "
        "100, not the exponents and not the initial conditions. So the figures the file attributes "
        "to Collier 1996 could not be confirmed here, and the two threshold slots that reproduce "
        "its quoted formulae exactly stay E_existence_only for that reason as well as for the "
        "registered one.",
    },
)

#: The gastrulation header claims a circuit the file does not contain. Reported, never edited.
HEADER_CLAIM = dict(np.KNOWN_WRONG_LINES[0])

#: Precedent for the existence/strength distinction that sits OUTSIDE this census's population,
#: recorded so a reader does not take it as a counted finding.
OUTSIDE_THE_POPULATION = (
    {
        "what": "the `M1` entry of the closed candidate list in scripts/mesoderm_diagnosis.py "
        "(registered at 0790d51) attributes two rules, `Nodal activates MIXL1` and `Mixl1 inhibits "
        "SOX2`, to Hart et al. 2002, Development 129:3597.",
        "why_it_is_not_counted": "it is a Python registration table, not a .bio program, so it is "
        "outside this census's population.",
        "what_it_shows_anyway": "the same distinction, and the prior lane was already honest about "
        "the level: the entry's own species note reads `level unsourced, module convention`.",
        "fetched_here": "Hart AH et al. 2002, Development 129:3597-3608, pmid 12117810, abstract "
        "via Europe PMC REST resultType=core. It is a real paper and it is the paper cited.",
        "what_the_abstract_does_not_claim": "it names neither Nodal nor Sox2. It is a Mixl1-null "
        "mouse study: thickened primitive streak, arrest at the early somite stage, deficient "
        "paraxial mesoderm, and in chimeras Mixl1(-/-) cells excluded from the hindgut. So it "
        "supports that Mixl1 acts in the mesendoderm split, which is what the entry's "
        "`acts_in_mesoderm_specification` field claims, and it supports neither of the two rules "
        "the entry attributes to it - not their existence and not their strength. Confirmed as a "
        "paper; its rules not confirmed by it. The full text is not open access (no PMCID), so a "
        "figure inside it could still name one of the two; the abstract does not.",
    },
)


def rows_for(path: Path) -> list[dict[str, Any]]:
    return np.census_file(path)


def tally(rows: list[dict[str, Any]], key: str) -> dict[str, int]:
    c: collections.Counter[str] = collections.Counter(r[key] for r in rows)
    return dict(sorted(c.items()))


def summarise(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Counts per class, the denominator, and what sits outside it."""
    outside = ("D_no_number", "N_non_numeric")
    inside = [r for r in rows if r["cls"] not in outside]
    conv: collections.Counter[str] = collections.Counter()
    for r in inside:
        for t in r.get("convention_tags", []):
            conv[f"{r['cls']}:{t}"] += 1
    return {
        "slots_scoped": len(rows),
        "denominator_numbers": len(inside),
        "classes": tally(rows, "cls"),
        "written": sum(1 for r in rows if r["written"]),
        "outside_the_denominator": {k: sum(1 for r in rows if r["cls"] == k) for k in outside},
        "convention_tags_by_class": dict(sorted(conv.items())),
    }


def main() -> None:
    programs = np.programs()
    all_rows: list[dict[str, Any]] = []
    per_file: dict[str, Any] = {}
    for p in programs:
        rows = rows_for(p)
        all_rows.extend(rows)
        per_file[p.as_posix()] = {
            "all": summarise(rows),
            "tier_1_dynamics_core": summarise([r for r in rows if r["tier"] == 1]),
            "tier_2_other_declared_numbers": summarise([r for r in rows if r["tier"] == 2]),
        }

    focus_rows = [r for r in all_rows if r["file"] == FOCUS]
    focus_rules: dict[str, Any] = {}
    for r in focus_rows:
        if r["kind"] != "rule":
            continue
        d = focus_rules.setdefault(
            r["declaration"], {"line": r["line"], "evidence_kind": r.get("evidence_kind"), "fields": {}}
        )
        d["fields"][r["field"]] = {
            "value": r["value"],
            "class": r["cls"],
            "convention_tags": r.get("convention_tags", []),
        }

    sourced = [r for r in all_rows if r["cls"].startswith("S_")]
    candidates = [r for r in all_rows if r.get("derivation_candidate")]
    unclassified = [r for r in all_rows if r["cls"] not in np.CLASSES]

    payload: dict[str, Any] = {
        "result": RESULT,
        "date": date.today().isoformat(),
        "lane": "lane-qualbound",
        "registration": f"data/results/{REGISTRATION}.json",
        "question": np.registration()["question"],
        "headline": {
            "corpus": summarise(all_rows),
            "tier_1_dynamics_core": summarise([r for r in all_rows if r["tier"] == 1]),
            "tier_2_other_declared_numbers": summarise([r for r in all_rows if r["tier"] == 2]),
            "programs": len(programs),
            "programs_with_no_scoped_slot": [p.as_posix() for p in programs if not rows_for(p)],
        },
        "per_file": per_file,
        "focus_module": {
            "file": FOCUS,
            "summary": summarise(focus_rows),
            "rules": focus_rules,
            "reading": "not one number in this module is traceable to a source that states it. "
            "Every rule number is either a citation for the rule's existence that states no number "
            "(E_existence_only) or an assertion with no source (A_asserted_bare), and the same "
            "holds for every gene max and basal, every protein half_life and every param. The "
            "module is qualitative throughout, by count and not by assertion.",
        },
        "sourced_class_in_full": [
            {
                k: r[k]
                for k in (
                    "file",
                    "line",
                    "kind",
                    "declaration",
                    "field",
                    "value",
                    "cls",
                    "evidence_kind",
                    "evidence_text",
                    "numbers_quoted_after_strip",
                    "convention_tags",
                )
                if k in r
            }
            for r in sourced
        ],
        "hand_adjudication": {
            "registered_rule": "a member of S_quoted is demoted to E_existence_only when the "
            "matched digits fall in citation furniture the strip missed.",
            "demoted_under_the_registered_rule": 0,
            "contested_but_kept": list(CONTESTED_S_QUOTED),
            "why_they_were_kept": "demoting them needs a wider rule than the one registered, and "
            "widening a test after seeing the data is the error this lane exists to avoid. One of "
            "the two contests was then refuted by the engine's own arithmetic and is withdrawn "
            "inside its entry, so the registered count of 30 stands and a reader who accepts the "
            "wider rule should read 29, not 28.",
            "s_derived_promotions": len(np.DERIVATIONS),
            "derivation_candidates_the_matcher_exposed": len(candidates),
            "strongest_case_the_registered_test_does_not_admit": list(STRONGEST_NON_S),
        },
        "source_checks": list(SOURCE_CHECKS),
        "reported_not_edited": [HEADER_CLAIM],
        "outside_the_population": list(OUTSIDE_THE_POPULATION),
        "unclassified": unclassified,
        "limits_of_this_census": [
            "it is CONSERVATIVE in one direction and the direction is known. A file that cites a "
            "quantitative paper without quoting its number lands in E_existence_only even if the "
            "paper does state it. The clearest instance is data/demo/repressilator.bio: its three "
            'rules carry `strength: 1.0; threshold: 1.0; hill: 2` against `experimental "Elowitz '
            '& Leibler 2000"`, and the genes carry `max: 200; basal: 0.2`, whose ratio of 1e-3 '
            "looks like the paper's own dimensionless leakiness. The paper's Box 1 may well state "
            "the Hill coefficient and that ratio, but the evidence strings quote nothing, the "
            "paper is not open access (pmid 10659856, isOpenAccess=N) and its abstract states no "
            "parameter value, so the census leaves all fifteen of those slots in E. A reader should "
            "treat E_existence_only as `the FILE does not carry the number's provenance`, which is "
            "what a census of a repository can establish, and not as `no source states it`.",
            "it cannot distinguish an honest assertion from a careless one. "
            "data/demo/concentration_threshold.bio and data/demo/stage4_division.bio say in every "
            "evidence string that their numbers are modelling devices and carry confidence 0.1 to "
            "0.2; they count A_asserted alongside anything else unsourced. The `confidence` field "
            "would separate them and is deliberately out of scope, because classifying the "
            "provenance of a self-reported confidence is classifying an opinion about provenance.",
            "`S_quoted` is a claim about the file, not about the paper, except where source_checks "
            "says a source was fetched. 16 of the 30 members are the Sender & Milo lifespans and "
            "none of those 16 was verified against the paper.",
        ],
        "observations_the_counts_make_visible": [
            "bio.std.methylation writes six of its params as `= unknown`, which is why they count "
            "N_non_numeric. That is not a missing number: the engine refuses to run until a program "
            "binds them, and data/demo/methylation_erosion.bio binds them as its own assumptions "
            "with `inferred` evidence, which is why that file carries nine A_asserted slots and one "
            "E. The corpus's most honest file is also its least sourced, and the census cannot tell "
            "the two apart - that is a limit of this census, not a finding about the file.",
            "data/demo/hoxd_order.bio contributes nine numbers and all nine are omitted gene "
            "`basal` slots at the language default 0.0, because the program declares loci and no "
            "kinetics at all. A per-file denominator keeps that from reading as nine conventional "
            "kinetic numbers.",
            "genomeos/std/signalling.bio contributes 148 tier-2 slots and not one tier-1 slot: its "
            "149 signals declare ligand, receptor and mode and almost never a threshold. The one "
            "large file in the corpus carries no kinetics.",
            "data/demo/concentration_threshold.bio is 12 A_asserted_bare slots out of 14 and says "
            "so in every evidence string (`a modelling device`), with confidence 0.1. An asserted "
            "number that declares itself asserted is the class working as intended.",
        ],
        "corrections_to_the_premises_this_lane_was_given": [
            "the module's inhibiting convention is not `threshold 2.0; hill 3` across the "
            "inhibiting rules. `Sox2 inhibits TBXT` writes threshold 4.0, so threshold 2.0 is on "
            "three of the four inhibiting rules, while hill 3 is on five of the seven rules - the "
            "activating `Tbxt activates SOX17` included. The convention is wider than the "
            "inhibiting rules and narrower in the threshold than the brief said.",
            "`max: 10` is not a language default. The IR supplies NO max at all when a gene omits "
            "it (genomeos/runtime/grn.py warns `gene_parameter_missing` and reads zero), so the 10 "
            "on three of gastrulation's four genes is the file's own repeated value. The language's "
            "rule defaults are strength 1.0, threshold 1.0 and hill 2.0, which is why the module's "
            "written `threshold: 2.0; hill: 3` counts as a repeated-in-file convention and not as "
            "an inherited one.",
            "`strength: 1.0` on six of the seven rules IS exactly the language default, so for "
            "strength the module's convention and the language's coincide; that is the one field "
            "where leaving it out would have changed nothing.",
        ],
        "claims_not_made": list(np.EXCLUSIONS),
        "result_manifest": {
            "sources": [
                {
                    "accession": "the hand-authored BioLang corpus of this repository: "
                    + ", ".join(np.POPULATION_DIRS),
                    "version": "the working tree at the code revision stamped below; every program "
                    "digested by its own path in inputs",
                },
                {
                    "accession": "genomeos/ir/model.py dataclass defaults, read for the value the "
                    "language supplies when a field is omitted",
                    "version": "the same code revision; transcribed into number_provenance.SCOPE "
                    "and asserted equal by tests/test_number_provenance.py",
                },
                {
                    "accession": "Europe PMC REST search (resultType=core) and PubMed Central, for "
                    "the abstracts and the one open full text in source_checks",
                    "version": "fetched 2026-10-02; each entry names what was fetched and what could not be",
                },
            ],
            "inputs": [
                mf.input_entry(p, partition=f"hand-authored BioLang program in {p.parent.as_posix()}")
                for p in programs
            ]
            + (
                [
                    mf.input_entry(
                        Path(f"data/results/{REGISTRATION}.json"),
                        partition="the registration this count was pre-registered in",
                    )
                ]
                if Path(f"data/results/{REGISTRATION}.json").exists()
                else []
            ),
            "assembly": "n/a: no genomic interval is read; genome coordinates are out of scope and "
            "the reason is in the registration's out_of_scope_fields.",
            "coordinates": "n/a: slots are addressed by file and line, not by genomic position.",
            "parameters": {
                "population_directories": list(np.POPULATION_DIRS),
                "scoped_field_count": len(np.SCOPE),
                "classes": dict(np.CLASSES),
                "cascade": list(np.CASCADE),
                "repeated_in_file_floor": np.REPEATED_IN_FILE_FLOOR,
                "denominator": np.registration()["denominator"],
                "model_requests": 0,
                "money_spent": "none; no paid API was called. The source checks are free web "
                "fetches of Europe PMC and PubMed Central.",
                "numbers_changed": 0,
                "bio_files_written": 0,
            },
            "exclusions": list(np.EXCLUSIONS),
            "partitions": dict(np.CLASSES),
            "code_cleanliness": mf.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
        },
    }

    path = save_result(RESULT, payload)
    h = payload["headline"]
    print(f"{RESULT}: {path}")
    for label in ("corpus", "tier_1_dynamics_core", "tier_2_other_declared_numbers"):
        s = h[label]
        print(f"  {label}: {s['denominator_numbers']} numbers of {s['slots_scoped']} scoped slots")
        for cls, n in s["classes"].items():
            print(f"      {cls:26s} {n}")
    f = payload["focus_module"]["summary"]
    print(
        f"  {FOCUS}: {f['denominator_numbers']} numbers; "
        + ", ".join(f"{k} {v}" for k, v in f["classes"].items())
    )
    print(
        f"  sourced class: {len(sourced)} members, {len(CONTESTED_S_QUOTED)} contested and kept, "
        f"{len(np.DERIVATIONS)} S_derived"
    )
    print(
        f"  source checks: {sum(1 for s in SOURCE_CHECKS if s['full_text_fetched'])} full texts "
        f"fetched, {len(SOURCE_CHECKS)} abstracts"
    )
    print(f"  unclassified slots: {len(unclassified)}")
    print("  0 numbers changed, 0 .bio files written, 0 model requests, no money")


if __name__ == "__main__":
    main()
