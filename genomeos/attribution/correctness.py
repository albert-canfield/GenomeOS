# SPDX-License-Identifier: AGPL-3.0-or-later
"""What a correct attribution means, and a scorer that keeps its three quantities apart (item 12 S4).

Registered 2026-09-29 before the scorer was applied to any real labelling (docs/ATTRIBUTION.md, "S4
registered"). The second external review (2026-09-28) asked for five things, each fixed here:

1. **Five claims, kept apart** (`AXES`): where the sequence came from (origin), what it is biochemically
   (molecular role), what it was seen to do (activity), which gene it acts on (target) and in which cell
   (context). The first three are R7's axes as the grammar states them (`genomeos.lang.grammar.AXES`);
   target is the gene a rule names, context the cell its `when` names (R1). `evidence_status` states
   what a claim rests on, not a claim about the biology, and is never judged.
2. **What can establish, refute or say nothing about each claim** (`TABLE`, one `Cell` per observation
   kind and axis). An observation decides a verdict only where the table says it establishes, refutes
   or cannot be read under its own observation model; a cell marked `CANNOT` is never consulted to
   decide a verdict, and one marked `SUGGESTS` is counted beside and never decides one.
3. **Not naming a coding gene is not absence of function.** No axis has a value meaning "no function";
   `unknown` and `unassigned` are not claims and are never judged (they are counted beside, `not_claims`);
   a well-powered CRISPRi null refutes the named gene in the cell screened and nothing else.
4. **Unresolved and competing explanations are kept** (`EXPLANATIONS`), including "the observation model
   is inadequate", which is also a verdict of its own (`MODEL_INADEQUATE`) where an assay's reading rule
   cannot produce a label at all.
5. **Target accuracy, role accuracy and coverage are reported apart** (`Report.quantities`), each a
   `Share` with its own numerator and denominator. A `Share` refuses addition, and nothing in this module
   combines two of them into one number.

The scorer uses item 13 C4's split discipline and nothing else: every held-out source whose units judge
a labelling is first passed to `holdout.check_provenance`, whose `LeakError` is never caught, and an
observation of a kind holdout does not load (an annotation, the registry) is refused with the same error
when the labelling read its source. The pilot (item 12 S3, lane-pilot) calls
`judge(claims, labels, sources=[S])` with a labelling built from `holdout.evidence(without=S)`.

An internal development benchmark only (`holdout.STATUS`): every source here has been read by this
project before, and no count this module returns is a fresh or external validation.
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from genomeos.attribution import holdout as ho
from genomeos.attribution import measured as ms
from genomeos.lang.grammar import AXES as R7_AXES

# ==================================================================================================
# The registration (2026-09-29, item 12 S4, lane-s4). Written before the scorer was run on any real
# labelling; docs/ATTRIBUTION.md carries the same text under "S4 registered".
# ==================================================================================================
REGISTERED = "2026-09-29"
STATUS = ho.STATUS
REUSE = ho.REUSE

ORIGIN, ROLE, ACTIVITY, TARGET, CONTEXT = "origin", "molecular_role", "activity", "target", "context"
AXES = (ORIGIN, ROLE, ACTIVITY, TARGET, CONTEXT)
AXIS_TEXT = {
    ORIGIN: "where the sequence came from (R7 `origin`)",
    ROLE: "what the element is biochemically (R7 `molecular_role`)",
    ACTIVITY: "what it was seen or predicted to do: the direction of its effect on its target in place, "
    "or its activity in a reporter (R7 `activity`)",
    TARGET: "which gene it acts on (the gene a rule names; R7 `target_relation` says how it was named)",
    CONTEXT: "the cell in which the relation holds (the rule's `when: cell_type`, R1)",
}
NEVER_JUDGED = {
    "evidence_status": "states what a claim rests on (a registry, a model, a measurement, constraint), not "
    "a claim about the biology; constraint lives here and is never read as a role",
}
#: the three quantities, and the axes each is read on. Origin and context accuracy are reported under
#: their own names beside them, never folded into target or role accuracy
QUANTITIES = {
    "target_accuracy": (TARGET,),
    "role_accuracy": (ROLE, ACTIVITY),
    "coverage": AXES,
}
ALSO_REPORTED = {"origin_accuracy": (ORIGIN,), "context_accuracy": (CONTEXT,)}
NO_SUM = (
    "target accuracy, role accuracy and coverage are three quantities with three denominators; they are "
    "never added, averaged or weighted into one score, and role accuracy is reported per axis "
    "(molecular role and activity apart), never pooled across them"
)

ANY = "*"  # target and context cells judge any gene or any cell
DIRECTION = ("activates_target", "represses_target")
REPORTER = ("active_in_reporter", "inactive_in_reporter")
SIGNATURE = ("promoter_like", "enhancer_like", "insulator_like", "open_chromatin")
ORIGIN_VALUES = tuple(v for v in R7_AXES[ORIGIN] if v != "unknown")
NOT_A_CLAIM = frozenset({"unknown", "unassigned", ""})

# --- the table: observation kind x axis ------------------------------------------------------------
ESTABLISHES, REFUTES, READ_INADEQUATE = "establishes", "refutes", "cannot_be_read"
JUDGES, SUGGESTS, CANNOT = "judges", "suggests_only", "cannot_establish"


@dataclass(frozen=True)
class Cell:
    """What one kind of observation can say about one axis, value by value.

    `establishes` and `refutes` list the axis values the observation can make correct or incorrect
    (`ANY` for any gene or cell); `inadequate` the values for which its reading rule can fail to produce
    a label; `suggests` the values it is consistent with and never decides. Everything else: cannot."""

    why: str
    establishes: tuple[str, ...] = ()
    refutes: tuple[str, ...] = ()
    inadequate: tuple[str, ...] = ()
    suggests: tuple[str, ...] = ()

    @property
    def kind(self) -> str:
        if self.establishes or self.refutes or self.inadequate:
            return JUDGES
        return SUGGESTS if self.suggests else CANNOT

    def reading(self, value: str) -> str:
        if value in self.establishes:
            return ESTABLISHES
        if value in self.refutes:
            return REFUTES
        if value in self.inadequate:
            return READ_INADEQUATE
        if value in self.suggests:
            return SUGGESTS
        return CANNOT

    def code(self) -> str:
        """The short form the documented table prints: E, R, M, S or -."""
        parts = [c for c, v in (("E", self.establishes), ("R", self.refutes), ("M", self.inadequate)) if v]
        if parts:
            return "/".join(parts)
        return "S" if self.suggests else "-"


def _no(why: str) -> Cell:
    return Cell(why)


_NULL_FUNCTION = (
    "a null on the genes tested is not absence of function: other genes, other cells, a redundant "
    "element, and roles no expression readout sees"
)
_NO_ORIGIN = "says nothing about where the sequence came from"
_MODEL = "a model output is a claim to be judged, never an observation that judges one"

#: the observation kinds, in the table's row order
OBSERVATIONS = (
    "crispri_decrease",
    "crispri_increase",
    "crispri_null_well_powered",
    "crispri_null_underpowered",
    "crispri_missing",
    "combinatorial_joint_change",
    "combinatorial_joint_null",
    "lentimpra_active",
    "lentimpra_inactive",
    "lentimpra_tiles_conflict",
    "vista_positive",
    "vista_negative",
    "satmut_functional_bases",
    "satmut_bases_inert",
    "gtex_associated",
    "gtex_not_associated",
    "chromatin_contact_present",
    "chromatin_contact_absent",
    "conservation_constrained",
    "conservation_not_detected",
    "allelic_imbalance",
    "allelic_balance",
    "alphagenome_deletion_prediction",
    "sequence_annotation",
    "registry_biochemical",
)
OBSERVATION_TEXT = {
    "crispri_decrease": "CRISPRi, significant decrease of the measured gene on silencing the element (R2)",
    "crispri_increase": "CRISPRi, significant increase of the measured gene (R2)",
    "crispri_null_well_powered": "CRISPRi, not significant with PowerAtEffectSize20 >= 0.8 (R2)",
    "crispri_null_underpowered": "CRISPRi, not significant and underpowered (R2)",
    "crispri_missing": "CRISPRi, no effect recorded (R2)",
    "combinatorial_joint_change": "two elements silenced together, the gene changes significantly",
    "combinatorial_joint_null": "two elements silenced together, a well-powered null",
    "lentimpra_active": "lentiMPRA, one cell, the matched tiles read active by R6's share rule",
    "lentimpra_inactive": "lentiMPRA, one cell, the matched tiles read silent by R6's share rule",
    "lentimpra_tiles_conflict": "lentiMPRA, one cell, the matched tiles split evenly about the threshold",
    "vista_positive": "VISTA transgenic mouse embryo e11.5, positive in some tissue",
    "vista_negative": "VISTA transgenic mouse embryo e11.5, negative",
    "satmut_functional_bases": "saturation mutagenesis, some measured base of the element is functional",
    "satmut_bases_inert": "saturation mutagenesis, every measured base of the element is inert",
    "gtex_associated": "GTEx v8 cis-eQTL, a significant variant inside the element for this gene",
    "gtex_not_associated": "GTEx v8, a protein-coding gene in the cis window with no such variant",
    "chromatin_contact_present": "chromatin contact (Hi-C, Micro-C, ChIA-PET) between element and gene",
    "chromatin_contact_absent": "no contact detected between element and gene",
    "conservation_constrained": "constraint across mammals (phyloP)",
    "conservation_not_detected": "no constraint detected",
    "allelic_imbalance": "allele-specific readout (item 13 C5): the two haplotypes differ",
    "allelic_balance": "allele-specific readout (item 13 C5): no imbalance",
    "alphagenome_deletion_prediction": "AlphaGenome's predicted expression change on deleting the element",
    "sequence_annotation": "RepeatMasker, the assembly and the segmental duplication track",
    "registry_biochemical": "the ENCODE cCRE registry's biochemical classes",
}

TABLE: dict[str, dict[str, Cell]] = {
    "crispri_decrease": {
        ORIGIN: _no(_NO_ORIGIN),
        ROLE: Cell(
            "the element activates a gene in place, consistent with an enhancer or a promoter; it cannot "
            "establish the biochemical signature R7's values name, or tell an enhancer from a promoter",
            suggests=("enhancer_like", "promoter_like"),
        ),
        ACTIVITY: Cell(
            "for the gene that fell, in the cell screened: removing the element lowers it. Says nothing "
            "about a reporter",
            establishes=("activates_target",),
            refutes=("represses_target", "no_effect_measured"),
        ),
        TARGET: Cell(
            "the gene responds when the element is silenced in place, in the cell screened. KRAB silences "
            "about a kilobase, so an element inside the tested interval is not told from its neighbour",
            establishes=(ANY,),
        ),
        CONTEXT: Cell("the relation holds in the cell screened", establishes=(ANY,)),
    },
    "crispri_increase": {
        ORIGIN: _no(_NO_ORIGIN),
        ROLE: Cell(
            "an increase is an effect, not evidence of a silencer (R2): silencer, insulator and competing "
            "promoter stay alternatives, and nothing here chooses among them",
            suggests=("silencer", "insulator_like", "competing_promoter"),
        ),
        ACTIVITY: Cell(
            "for the gene that rose, in the cell screened: removing the element raises it",
            establishes=("represses_target",),
            refutes=("activates_target", "no_effect_measured"),
        ),
        TARGET: Cell("the gene responds when the element is silenced in place", establishes=(ANY,)),
        CONTEXT: Cell("the relation holds in the cell screened", establishes=(ANY,)),
    },
    "crispri_null_well_powered": {
        ORIGIN: _no(_NO_ORIGIN),
        ROLE: _no(_NULL_FUNCTION),
        ACTIVITY: Cell(
            "establishes only `no_effect_measured` for that gene in that cell. A null has no direction, so "
            "it never judges activates_target or represses_target: the direction is judged only where the "
            "gene responded, and the null is read on the target axis instead",
            establishes=("no_effect_measured",),
        ),
        TARGET: Cell(
            "refutes the named gene only in the cell screened and only for an effect of 20% or more; a "
            "redundant element can hide a real target, which only a combinatorial perturbation separates; "
            + _NULL_FUNCTION,
            refutes=(ANY,),
        ),
        CONTEXT: Cell(
            "refutes the stated cell when the gene responds to the element in another cell",
            refutes=(ANY,),
        ),
    },
    "crispri_null_underpowered": {
        a: _no("a screen too weak to have seen a 20% effect says nothing on any axis") for a in AXES
    },
    "crispri_missing": {a: _no("no effect was recorded for the pair") for a in AXES},
    "combinatorial_joint_change": {
        ORIGIN: _no(_NO_ORIGIN),
        ROLE: _no("a joint effect says two elements matter together, not what either is biochemically"),
        ACTIVITY: Cell(
            "the direction of the joint effect where single silencing left the gene unchanged: the route to "
            "additive, redundant and cooperative explanations single perturbations cannot separate (C6)",
            establishes=DIRECTION,
            refutes=("no_effect_measured",),
        ),
        TARGET: Cell("the gene responds to the pair, in the cell screened", establishes=(ANY,)),
        CONTEXT: Cell("the joint relation holds in the cell screened", establishes=(ANY,)),
    },
    "combinatorial_joint_null": {
        ORIGIN: _no(_NO_ORIGIN),
        ROLE: _no(_NULL_FUNCTION),
        ACTIVITY: Cell("no effect of the pair on that gene", establishes=("no_effect_measured",)),
        TARGET: Cell(
            "refutes the gene for the pair in the cell screened, removing redundancy between the two as "
            "the explanation of a single null",
            refutes=(ANY,),
        ),
        CONTEXT: _no("a joint null in one cell says nothing about another"),
    },
    "lentimpra_active": {
        ORIGIN: _no(_NO_ORIGIN),
        ROLE: Cell(
            "the sequence acts autonomously outside its locus, consistent with an enhancer or promoter; it "
            "is not the biochemical signature and not the element in place",
            suggests=("enhancer_like", "promoter_like"),
        ),
        ACTIVITY: Cell(
            "active in a reporter in that cell. Never activates_target or represses_target: a 200 bp copy "
            "driving its own reporter measures activity, not regulation",
            establishes=("active_in_reporter",),
            refutes=("inactive_in_reporter",),
        ),
        TARGET: _no("a reporter tile can establish activity, never a target: it drives its own reporter"),
        CONTEXT: Cell(
            "the sequence can act in that cell's trans environment; the relation in place is not measured",
            suggests=(ANY,),
        ),
    },
    "lentimpra_inactive": {
        ORIGIN: _no(_NO_ORIGIN),
        ROLE: _no("a sequence silent in a reporter can carry its signature and act in place"),
        ACTIVITY: Cell(
            "inactive in a reporter in that cell only; says nothing about the element in place",
            establishes=("inactive_in_reporter",),
            refutes=("active_in_reporter",),
        ),
        TARGET: _no("a reporter never measures a target"),
        CONTEXT: _no("silence in a reporter does not refute a relation in place in that cell"),
    },
    "lentimpra_tiles_conflict": {
        ORIGIN: _no(_NO_ORIGIN),
        ROLE: _no("an unreadable reporter says nothing about the signature"),
        ACTIVITY: Cell(
            "R6's share rule cannot read a label when the tiles split evenly: the observation model is "
            "inadequate for this element in this cell",
            inadequate=REPORTER,
        ),
        TARGET: _no("a reporter never measures a target"),
        CONTEXT: _no("an unreadable reporter says nothing about a context"),
    },
    "vista_positive": {
        ORIGIN: _no(_NO_ORIGIN),
        ROLE: Cell(
            "in vivo reporter activity, consistent with an enhancer; not the biochemical signature",
            suggests=("enhancer_like",),
        ),
        ACTIVITY: Cell(
            "active in a transgenic reporter, mouse e11.5, in the tissues named; never regulation in place",
            establishes=("active_in_reporter",),
            refutes=("inactive_in_reporter",),
        ),
        TARGET: _no("a transgenic reporter drives its own minimal promoter; it never names a target"),
        CONTEXT: _no("a mouse embryo at e11.5 is not a context any compiled rule states"),
    },
    "vista_negative": {
        ORIGIN: _no(_NO_ORIGIN),
        ROLE: _no("inactive at one stage in a mouse embryo says nothing about the signature"),
        ACTIVITY: Cell(
            "inactive in the transgenic reporter at e11.5 only",
            establishes=("inactive_in_reporter",),
            refutes=("active_in_reporter",),
        ),
        TARGET: _no("a transgenic reporter never names a target"),
        CONTEXT: _no("one developmental stage in another species"),
    },
    "satmut_functional_bases": {
        ORIGIN: _no(_NO_ORIGIN),
        ROLE: Cell(
            "bases whose identity changes reporter activity: motifs, consistent with a regulatory element, "
            "not a role",
            suggests=("enhancer_like", "promoter_like"),
        ),
        ACTIVITY: Cell(
            "the bases matter in a reporter; the element's own activity is not what was scored",
            suggests=("active_in_reporter",),
        ),
        TARGET: _no("single substitutions in a reporter never measure a gene in place"),
        CONTEXT: _no("the reporter's cell is not the relation's context"),
    },
    "satmut_bases_inert": {
        a: _no("a measurement of those bases only, never of the element (SATMUT_CANNOT_DISAGREE)")
        for a in AXES
    },
    "gtex_associated": {
        ORIGIN: _no(_NO_ORIGIN),
        ROLE: _no("association says a variant tracks expression, nothing about the biochemical signature"),
        ACTIVITY: Cell(
            "a variant inside the element tracks the gene's expression: association, not perturbation, and "
            "the causal variant may be a linked one",
            suggests=DIRECTION,
        ),
        TARGET: Cell("association with the gene, not perturbation; linkage is not excluded", suggests=(ANY,)),
        CONTEXT: _no("the GTEx units read here pool 49 tissues (holdout.gtex_units): they carry no context"),
    },
    "gtex_not_associated": {
        a: _no("not significant and not tested are not told apart (holdout's GTEx endpoint)") for a in AXES
    },
    "chromatin_contact_present": {
        ORIGIN: _no(_NO_ORIGIN),
        ROLE: _no("contact cannot tell an enhancer from an insulator or a structural anchor"),
        ACTIVITY: _no("proximity is not activity"),
        TARGET: Cell(
            "the element and the gene are near in three dimensions; neither sufficient nor necessary",
            suggests=(ANY,),
        ),
        CONTEXT: Cell("the conformation exists in the cell measured", suggests=(ANY,)),
    },
    "chromatin_contact_absent": {
        a: _no("an undetected contact does not refute regulation, nor any other claim") for a in AXES
    },
    "conservation_constrained": {
        a: _no(
            "constraint suggests selection, recorded on evidence_status; it never establishes a role, an "
            "activity, a target, a context or an origin (an exapted repeat can be constrained)"
        )
        for a in AXES
    },
    "conservation_not_detected": {a: _no("no sign of selection is not a sign of no function") for a in AXES},
    "allelic_imbalance": {
        ORIGIN: _no(_NO_ORIGIN),
        ROLE: _no("a haplotype difference says nothing about the signature"),
        ACTIVITY: Cell(
            "the haplotypes differ in the donor tissue; with a median 1,580 to 2,012 heterozygous SNVs "
            "within 1 Mb of a TSS (the C5 probe), a within-person imbalance cannot be credited to one "
            "element",
            suggests=DIRECTION,
        ),
        TARGET: Cell("a haplotype, not an element, is tied to the gene", suggests=(ANY,)),
        CONTEXT: Cell("the donor tissue measured", suggests=(ANY,)),
    },
    "allelic_balance": {
        a: _no("the element may carry no variant, or a neutral one; balance says nothing") for a in AXES
    },
    "alphagenome_deletion_prediction": {a: _no(_MODEL) for a in AXES},
    "sequence_annotation": {
        ORIGIN: Cell(
            "states the origin; independent only of a labelling that did not read it (the compiled labels "
            "read it, so it is their input and never their check)",
            establishes=ORIGIN_VALUES,
            refutes=ORIGIN_VALUES,
        ),
        ROLE: _no("a repeat class is not a biochemical signature"),
        ACTIVITY: _no("sequence class is not activity"),
        TARGET: _no("sequence class names no gene"),
        CONTEXT: _no("sequence class is the same in every cell"),
    },
    "registry_biochemical": {
        ORIGIN: _no(_NO_ORIGIN),
        ROLE: Cell(
            "states the biochemical signature; independent only of a labelling that did not read it (the "
            "compiled labels read it). Nothing in this table establishes a silencer, a competing promoter, "
            "a structural role or a coding candidate",
            establishes=SIGNATURE,
            refutes=SIGNATURE,
        ),
        ACTIVITY: _no("a signature is not a measured effect"),
        TARGET: _no("a signature names no gene"),
        CONTEXT: _no("the registry's classes are pooled over biosamples"),
    },
}

_LOADED_ASSAYS = ("crispri", "lentimpra", "vista", "satmut", "gtex")
#: the kinds holdout's sources supply, and why every other kind is not loaded by this scorer
LOADED = tuple(k for k in OBSERVATIONS if k.split("_", 1)[0] in _LOADED_ASSAYS)
NOT_LOADED = {
    "combinatorial_joint_change": "no paired perturbation is among the cached screens (C6 is phase C)",
    "combinatorial_joint_null": "no paired perturbation is among the cached screens (C6 is phase C)",
    "chromatin_contact_present": "no holdout source carries contact",
    "chromatin_contact_absent": "no holdout source carries contact",
    "conservation_constrained": "an input of the compiled labels (evidence_status), and it judges nothing",
    "conservation_not_detected": "an input of the compiled labels (evidence_status), and it judges nothing",
    "allelic_imbalance": "the C5 probe counted reach only; no allelic outcome has been read",
    "allelic_balance": "the C5 probe counted reach only; no allelic outcome has been read",
    "alphagenome_deletion_prediction": "the model output is the unchanged labels' own input; 0 requests",
    "sequence_annotation": "an input of the compiled labels' origin, so never their check",
    "registry_biochemical": "an input of the compiled labels' molecular role, so never their check",
}
_PAIR_PREFIXES = ("crispri", "combinatorial", "gtex", "chromatin", "allelic")
#: kinds whose observations are about one gene; the rest are about the element
PAIR_KINDS = frozenset(k for k in OBSERVATIONS if k.startswith(_PAIR_PREFIXES))
#: perturbation kinds: an element they measured was screened
SCREEN_KINDS = frozenset(k for k in OBSERVATIONS if k.startswith(("crispri", "combinatorial")))
#: kinds that carry a stated value to compare with the claim's (establish when equal, refute otherwise)
VALUE_KINDS = frozenset({"sequence_annotation", "registry_biochemical"})
CRISPRI_KIND = {
    ms.DECREASE: "crispri_decrease",
    ms.INCREASE: "crispri_increase",
    ms.NULL_INFORMATIVE: "crispri_null_well_powered",
    ms.NULL_INCONCLUSIVE: "crispri_null_underpowered",
    ms.MISSING: "crispri_missing",
}
VISTA_CONTEXT = "mouse_e11.5"
MATCH_RULE = (
    "an observation is of the element when its interval meets the measured layer's rule "
    "(measured.measures: reciprocal overlap >= 0.5) for CRISPRi pairs, lentiMPRA tiles and VISTA "
    "elements; a saturation-mutagenesis base when it lies inside the element; a GTEx unit when its "
    "element shares a base (holdout's and the ablation's rules). lentiMPRA is read per cell from the "
    "matched tiles by R6's share rule (measured.reporter_label). Genes by symbol; cells equal after "
    "lowercasing and dropping every character that is not a letter or a digit (the ablation's rule)"
)

# --- verdicts, reasons, explanations ---------------------------------------------------------------
CORRECT, INCORRECT, UNRESOLVED = "correct", "incorrect", "unresolved"
MODEL_INADEQUATE = "observation_model_inadequate"
NOT_JUDGED = "not_judged"
VERDICTS = (CORRECT, INCORRECT, UNRESOLVED, MODEL_INADEQUATE, NOT_JUDGED)
JUDGED = (CORRECT, INCORRECT, UNRESOLVED, MODEL_INADEQUATE)
#: an explanation that stays live beside every verdict an observation decides: the assay's observation
#: model (one effect per element, gene and cell; one label per reporter cell) may not describe what it saw
OBSERVATION_MODEL_INADEQUATE = "the observation model is inadequate"
EXPLANATIONS = {
    CORRECT: ("the claim holds", OBSERVATION_MODEL_INADEQUATE),
    INCORRECT: ("the claim is wrong", OBSERVATION_MODEL_INADEQUATE),
    UNRESOLVED: (
        "the claim holds and the refuting observation is wrong",
        "the claim is wrong and the establishing observation is wrong",
        "the answer differs between studies, conditions or contexts",
        OBSERVATION_MODEL_INADEQUATE,
    ),
    MODEL_INADEQUATE: (OBSERVATION_MODEL_INADEQUATE,),
}
#: explanations added on one axis
AXIS_EXPLANATIONS = {
    (TARGET, INCORRECT): ("a redundant element compensates; only a combinatorial perturbation separates it",),
    (CONTEXT, INCORRECT): ("a redundant element compensates in the stated cell",),
    (ACTIVITY, INCORRECT): ("the effect is indirect, through another gene",),
}
#: why a claim was not judged: one reason per claim, the first that applies in this order
OUTSIDE_SCOPE = "judging_observations_outside_the_claims_scope"
CONDITION_UNMET = "not_judged_until_the_target_responds"
SUGGESTS_ONLY = "only_observations_that_suggest"
CANNOT_ONLY = "only_observations_that_cannot_judge_this_axis"
NO_OBSERVATION = "no_observation"
REASONS = (OUTSIDE_SCOPE, CONDITION_UNMET, SUGGESTS_ONLY, CANNOT_ONLY, NO_OBSERVATION)
NULL_ELSEWHERE = "well_powered_null_only_in_another_context"
GENE_NOT_TESTED = "element_screened_this_gene_not_tested"

AXIS_RULE = {
    TARGET: "established by a significant change of the named gene on silencing the element in any cell; "
    "refuted by a well-powered null of that gene in the cell the claim states, and only when the gene "
    "responds in no cell. A response in another cell with a null in the stated cell leaves the target "
    "established and refutes the context; a response and a null both in the stated cell is unresolved. "
    f"A null in another cell is not a refutation (`{NULL_ELSEWHERE}`, counted); a claim with no stated "
    "cell can be established and never refuted. `refutable` counts the claims whose stated cell was "
    "screened on the named gene with an outcome that can establish or refute",
    CONTEXT: "judged only where the gene responds to the element in some cell: established by a response in "
    "the stated cell, refuted by a well-powered null there. A claim whose target never responded is not "
    "judged on context, so one failure is never counted on both axes",
    ACTIVITY: "the direction (activates_target, represses_target) is judged only where the named gene "
    "responded, in the stated cell when it responded there and otherwise in every cell it responded in; "
    "both directions in one cell is `observation_model_inadequate`, in different cells `unresolved`. "
    "Reporter values are judged by reporter observations in the stated cell when one was measured there, "
    "otherwise as R7 defines them (active in at least one reporter context; inactive in every one)",
    ROLE: "judged only by an observation the table lets establish a biochemical signature (the registry); "
    "activity in a reporter or in place suggests and never decides",
    ORIGIN: "judged only by a sequence annotation the labelling did not read",
}
MAY_BE_CALLED = (
    "an internal development benchmark reading: of N claims a labelling states on axis A, J can be judged "
    "by an observation the S4 table lets establish or refute A (coverage J / N); of the J, C are "
    "established and I refuted (accuracy C / (C + I)), with U unresolved and M observation-model-"
    "inadequate beside"
)
MAY_NOT_BE_CALLED = (
    "a validation, fresh or external: every source here has been read by this project before",
    "the accuracy of the labelling on claims no observation judges: the judged claims are the ones screens "
    "chose to test, near the genes they chose, mostly in K562",
    "a single score, or a comparison between axes: each quantity has its own denominator and evidence",
    "evidence that an unjudged element has no function",
)
EXPECTED = (
    "written before the run, after reading the ablation's published counts (93 of 128 measured links "
    "survive, 28 contradicted of which 23 by a well-powered null and 5 by the opposite sign, 7 "
    "inconclusive; 39 rules survive in their own cell and none is contradicted there): most claims on "
    "every axis are not judged; origin and molecular role are judged by no loaded source, 0 claims by "
    "construction; target is judged on the order of 100 of 440,377 rules and refuted on about none, "
    "because a refutation needs a null in the stated cell, with the claims whose only null is in another "
    "cell (about 20) counted beside; activity is judged on about the same claims as target, a few of them "
    "refuted by the opposite sign; context is established on about 40 and refuted on about none"
)

# ==================================================================================================
# The direction rule, versioned (registered 2026-09-29, item 12 S4 follow-up, lane-judge2), before any
# claim was judged under v2. docs/ATTRIBUTION.md carries the same text under "The direction rule,
# versioned". Everything above is v1 and stays as registered; v1 is the rule S4's committed result was
# judged by, and v2 is the rule for every call from this registration on.
# ==================================================================================================
RULE_REGISTERED = "2026-09-29"
RULE_V1, RULE_V2 = "v1", "v2"
RULES = (RULE_V1, RULE_V2)
#: the rule a call uses when it names none; a script that reproduces a result judged under v1 names RULE_V1
DEFAULT_RULE = RULE_V2
#: v2's reason for a direction claim the stated cell does not decide (a new reason; no v1 reason is reused)
NOT_ASSESSED = "not_assessed_in_this_context"
#: what another cell's response says about the claimed direction, kept beside the verdict under v2
AGREES, DISAGREES = "agrees", "disagrees"
RULE_DECISION = (
    'Albert, 2026-09-29: "Judge rule: preserve the historical result, correct future judging. The '
    "current rule falls back to another cell when the claimed cell lacks measurements. That cannot "
    "establish that a cell-specific direction is wrong. Registration makes the analysis traceable; it "
    "does not make that interpretation valid. Introduce a versioned rule: judge direction in the stated "
    "cell; otherwise return 'not assessed in this context.' Keep cross-cell disagreement as a separate "
    "finding. Apply this symmetrically: the trace says 54 correct and five incorrect verdicts came from "
    "other cells. All 59 should become unassessed under the revised rule—not just the five errors. "
    "Preserve the original results beside the correction. ID1 remains a separate, useful finding: its "
    "cached K562 prediction disagrees with the K562 experiment. That does not directly refute the "
    'compiled whole-blood claim."'
)
DIRECTION_RULE = {
    RULE_V1: "v1 (S4 as registered, AXIS_RULE[activity], unchanged): " + AXIS_RULE[ACTIVITY],
    RULE_V2: "v2 (from 2026-09-29): the direction (activates_target, represses_target) is judged only in "
    "the cell the claim states, and only by the named gene's significant changes there: a change the "
    "claimed way establishes it, a change the other way refutes it, and both in the stated cell is "
    "`observation_model_inadequate`. When the claim states no cell, or the named gene responded in some "
    "cell but has no significant change in the stated cell, the claim is not judged, with the reason "
    "`not_assessed_in_this_context`, whatever the gene did in any other cell: agreement in another cell "
    "no longer establishes a direction, just as disagreement in another cell no longer refutes one. Each "
    "significant change of the named gene in another cell is kept on the verdict as a separate cross-cell "
    "finding (`cross_cell`), marked `agrees` or `disagrees` with the claimed direction; it never decides "
    "the verdict. A claim whose gene responded in no cell keeps v1's reason, since no cell decides it "
    "under either rule. Reporter values and no_effect_measured, and the target, context, origin and "
    "molecular-role axes, are judged exactly as under v1",
}
REFUTABLE_RULE = {
    RULE_V1: "v1: computed on target and context only (the stated cell was screened on the named gene with "
    "an outcome that can establish or refute); every other verdict carries the default false, which on "
    "activity means 'not computed' (the repression trace, 2026-09-29)",
    RULE_V2: "v2: on a direction claim, true when the stated cell was screened on the named gene with an "
    "outcome that can establish or refute the direction (a significant change there; the table lets no "
    "null decide a direction), false otherwise; target and context keep v1's definition; where the field "
    "is not computed (origin, molecular role, and the activity values that are not a direction) v2 writes "
    "null, never false",
}
RULE_RECORD = (
    "a result judged under v2 names its rule: `judge_rule` in its manifest parameters and `rule` in its "
    "body. A result that names no rule was judged under v1, the only rule before 2026-09-29; the scripts "
    "that reproduce such results pin RULE_V1 and add nothing to their output, so their committed files "
    "still reproduce (scripts/s4_correctness_run.py without `--rule v2`, scripts/prior_only_test.py, "
    "scripts/pilot_biological_gate.py, scripts/repression_trace.py)"
)
SEEN_BEFORE_V2 = (
    "seen before this registration, and used to write EXPECTED_V2: the repression trace "
    "(data/results/repression_trace.json, commits 54c47f1..641909e) had already counted the committed S4 "
    "activity verdicts by where they were decided: 39 in the stated cell, all established; 54 established "
    "and 5 refuted only in another cell. It had also computed that judging direction only in the stated "
    "cell would give 39 judged, 39 established, with 59 moved to not judged. The 39 / 59 split and the "
    "S4 counts expected below were therefore known before this rule was registered. Also read before "
    "registration: the committed S4-beside blocks of prior_only_test.json and pilot_biological_gate.json "
    "(their activity, target-refutable and context counts), from which the expected moves there are "
    "derived. The registration makes the re-judge traceable; it does not make the interpretation valid"
)
EXPECTED_V2 = (
    "on the committed S4 claims and sources: activity is judged on 39 of 440,377 claims, all 39 "
    "established (0 refuted, 0 unresolved, 0 observation-model-inadequate); the 59 verdicts v1 decided only "
    "in another cell (54 established, 5 refuted) are not judged, reason not_assessed_in_this_context, each "
    "with its other-cell observations kept as a cross-cell finding, 54 agreeing and 5 disagreeing with the "
    "claimed direction; every other activity reason keeps its v1 count (1,377 outside the claim's scope, "
    "23 until the target responds, 3,194 only suggest, 18,747 only cannot judge, 416,938 no observation); "
    "activity refutable is 39; no count moves on target (98 of 98), context (39 of 39), origin or "
    "molecular role (0 judged). ID1's represses_target claim in Whole_Blood (EH38E3426791) is not assessed "
    "in this context, with Gasperini2019's K562 decrease of ID1 kept as a disagreeing cross-cell finding. "
    "How many of the 39 judged verdicts also carry a cross-cell finding was not counted before "
    "registration. Beside S4, derived and not rerun: in the prior-only test's unchanged_on_the_same_blocks "
    "(Gasperini2019) activity would go from 41 of 43 to 24 of 24, 19 moving (17 established, 2 refuted); "
    "in the pilot gate's unchanged_on_same_blocks from 27 of 28 to 12 of 12 on Gasperini2019 (16 moving: "
    "15 established, 1 refuted), from 3 of 3 to 0 of 0 on Schraivogel2020 (3 moving, all established) and "
    "1 of 1 unchanged on Xie; the pilot's and the prior's own labellings state no activity claim, so none "
    "of their counts move"
)


# --- claims and observations ------------------------------------------------------------------------
@dataclass(frozen=True, slots=True)
class Claim:
    """One statement a labelling makes about one element on one axis.

    `value` is an R7 value on origin, molecular role and activity, the gene's symbol on target, the cell
    on context. `gene` is the gene an activity or context claim is relative to (a target claim repeats
    its value there); `cell` the cell a target or activity claim is stated in ("" or "unknown": none)."""

    element: str
    chrom: str
    start: int  # 0-based, half-open
    end: int
    axis: str
    value: str
    gene: str = ""
    cell: str = ""
    block: str = "element"


@dataclass(frozen=True, slots=True)
class Observation:
    """One reading of one source about one element: its kind (a row of `TABLE`), for pair kinds the
    gene, the cell it was made in, the CRISPRi split, and for annotation and registry kinds the value(s)
    stated, comma-joined. The interval is used only to index `extra` observations."""

    kind: str
    source: str
    gene: str = ""
    cell: str = ""
    split: str = ""
    value: str = ""
    chrom: str = ""
    start: int = 0
    end: int = 0

    def to_dict(self) -> dict[str, Any]:
        d = {"kind": self.kind, "source": self.source}
        d.update({k: getattr(self, k) for k in ("gene", "cell", "split", "value") if getattr(self, k)})
        return d


def norm_cell(c: str) -> str:
    return re.sub(r"[^a-z0-9]", "", (c or "").lower())


def stated(cell: str) -> bool:
    return norm_cell(cell) not in ("", norm_cell(ms.CONTEXT_UNKNOWN))


def reading(kind: str, axis: str, value: str) -> str:
    """What the table lets one observation of `kind` say about `value` on `axis`."""
    return TABLE[kind][axis].reading(ANY if axis in (TARGET, CONTEXT) else value.split("/", 1)[0])


# --- one claim's verdict ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Verdict:
    claim: Claim
    verdict: str
    reason: str = ""  # NOT_JUDGED only: one of REASONS
    detail: str = ""
    explanations: tuple[str, ...] = ()
    deciding: tuple[Observation, ...] = ()
    refutable: bool = False  # target and context: the stated cell was screened on the gene

    def to_dict(self) -> dict[str, Any]:
        c = self.claim
        d: dict[str, Any] = {
            "element": c.element,
            "locus": f"{c.chrom}:{c.start}-{c.end}",
            "block": c.block,
            "axis": c.axis,
            "value": c.value,
            "gene": c.gene,
            "cell": c.cell,
            "verdict": self.verdict,
            "reason": self.reason,
            "detail": self.detail,
            "explanations": list(self.explanations),
            "deciding": [o.to_dict() for o in self.deciding],
            "refutable": self.refutable,
        }
        return d

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> Verdict:
        m = ho.LOCUS_RE.match(d["locus"])
        if not m:
            raise ValueError(f"bad locus {d['locus']!r}")
        if d["verdict"] not in VERDICTS:
            raise ValueError(f"unknown verdict {d['verdict']!r}")
        claim = Claim(
            d["element"],
            m.group(1),
            int(m.group(2)),
            int(m.group(3)),
            d["axis"],
            d["value"],
            d.get("gene", ""),
            d.get("cell", ""),
            d.get("block", "element"),
        )
        return cls(
            claim,
            d["verdict"],
            d.get("reason", ""),
            d.get("detail", ""),
            tuple(d.get("explanations", ())),
            tuple(Observation(**o) for o in d.get("deciding", ())),
            bool(d.get("refutable", False)),
        )


def _decided(
    claim: Claim, verdict: str, deciding: list[Observation], detail: str = "", refutable: bool = False
) -> Verdict:
    for o in deciding:  # the table decides: an observation it says cannot judge never decides
        r = reading(o.kind, claim.axis, claim.value)
        if r not in (ESTABLISHES, REFUTES, READ_INADEQUATE):
            raise AssertionError(f"{o.kind} is {r} on {claim.axis}={claim.value}; it cannot decide a verdict")
    expl = EXPLANATIONS[verdict] + AXIS_EXPLANATIONS.get((claim.axis, verdict), ())
    return Verdict(claim, verdict, "", detail, expl, tuple(deciding), refutable)


def _not_judged(claim: Claim, reason: str, detail: str = "", refutable: bool = False) -> Verdict:
    return Verdict(claim, NOT_JUDGED, reason, detail, (), (), refutable)


def touching(claim: Claim, obs: Iterable[Observation]) -> list[Observation]:
    """The observations about this claim's element, and for a pair kind about this claim's gene."""
    return [o for o in obs if o.kind not in PAIR_KINDS or not claim.gene or o.gene == claim.gene]


def _fallback(claim: Claim, obs: list[Observation]) -> Verdict:
    """Not judged, and nothing in scope: say what touched the claim, suggestive first."""
    near = touching(claim, obs)
    if not near:
        return _not_judged(claim, NO_OBSERVATION)
    kinds = ",".join(sorted({o.kind for o in near}))
    if any(reading(o.kind, claim.axis, claim.value) == SUGGESTS for o in near):
        return _not_judged(claim, SUGGESTS_ONLY, kinds)
    return _not_judged(claim, CANNOT_ONLY, kinds)


def _on_gene(obs: list[Observation], gene: str, want: str, axis: str = TARGET) -> list[Observation]:
    return [o for o in obs if o.kind in PAIR_KINDS and o.gene == gene and reading(o.kind, axis, gene) == want]


def _in(obs: list[Observation], cell: str) -> list[Observation]:
    return [o for o in obs if norm_cell(o.cell) == norm_cell(cell)] if stated(cell) else []


def _screened(obs: list[Observation], gene: str | None = None) -> bool:
    """A perturbation measured the element (on `gene`, when given), with any outcome."""
    return any(o.kind in SCREEN_KINDS and (gene is None or o.gene == gene) for o in obs)


def _unscreened(claim: Claim, obs: list[Observation], gene: str) -> Verdict:
    if _screened(obs) and not _screened(obs, gene):
        return _not_judged(claim, OUTSIDE_SCOPE, GENE_NOT_TESTED)
    return _fallback(claim, obs)


def _target(claim: Claim, obs: list[Observation]) -> Verdict:
    gene, cell = claim.value, claim.cell
    est = _on_gene(obs, gene, ESTABLISHES)
    ref = _on_gene(obs, gene, REFUTES)
    est_here, ref_here = _in(est, cell), _in(ref, cell)
    refutable = bool(est_here or ref_here)
    if est and ref_here:
        if est_here:
            return _decided(
                claim, UNRESOLVED, est_here + ref_here, "response and null in the stated cell", True
            )
        return _decided(
            claim, CORRECT, est, "responds in another cell; the stated cell's null is read on context", True
        )
    if est:
        return _decided(claim, CORRECT, est, "", refutable)
    if ref_here:
        return _decided(claim, INCORRECT, ref_here, "well-powered null in the stated cell", True)
    if ref:
        return _not_judged(claim, OUTSIDE_SCOPE, NULL_ELSEWHERE)
    return _unscreened(claim, obs, gene)


def _context(claim: Claim, obs: list[Observation]) -> Verdict:
    gene, cell = claim.gene, claim.value
    responds = _on_gene(obs, gene, ESTABLISHES)
    est = [o for o in _in(responds, cell) if reading(o.kind, CONTEXT, cell) == ESTABLISHES]
    ref = _in(_on_gene(obs, gene, REFUTES, CONTEXT), cell)
    if not responds:
        if _on_gene(obs, gene, REFUTES):
            return _not_judged(claim, CONDITION_UNMET, "the gene responded in no cell tested", bool(ref))
        return _unscreened(claim, obs, gene)
    if est and ref:
        return _decided(claim, UNRESOLVED, est + ref, "response and null in the stated cell", True)
    if est:
        return _decided(claim, CORRECT, est, "", True)
    if ref:
        return _decided(
            claim, INCORRECT, ref, "responds in another cell, well-powered null in the stated one", True
        )
    return _not_judged(claim, OUTSIDE_SCOPE, "responds only in other cells; the stated cell is untested")


def _direction(claim: Claim, obs: list[Observation]) -> Verdict:
    v = claim.value
    responds = _on_gene(obs, claim.gene, ESTABLISHES)
    if not responds:
        if _on_gene(obs, claim.gene, REFUTES):
            return _not_judged(claim, CONDITION_UNMET, "the gene responded in no cell tested")
        return _unscreened(claim, obs, claim.gene)
    scope = _in(responds, claim.cell) or responds
    est = [o for o in scope if reading(o.kind, ACTIVITY, v) == ESTABLISHES]
    ref = [o for o in scope if reading(o.kind, ACTIVITY, v) == REFUTES]
    if est and ref:
        if {norm_cell(o.cell) for o in est} & {norm_cell(o.cell) for o in ref}:
            return _decided(claim, MODEL_INADEQUATE, est + ref, "significant in both directions in one cell")
        return _decided(claim, UNRESOLVED, est + ref, "the direction differs between cells")
    if est:
        return _decided(claim, CORRECT, est)
    if ref:
        return _decided(claim, INCORRECT, ref, "the gene moved the other way")
    return _fallback(claim, obs)


def _no_effect(claim: Claim, obs: list[Observation]) -> Verdict:
    pool = [o for o in obs if o.kind in PAIR_KINDS and (not claim.gene or o.gene == claim.gene)]
    pool = _in(pool, claim.cell) if stated(claim.cell) else pool
    est = [o for o in pool if reading(o.kind, ACTIVITY, claim.value) == ESTABLISHES]
    ref = [o for o in pool if reading(o.kind, ACTIVITY, claim.value) == REFUTES]
    if est and ref:
        return _decided(claim, UNRESOLVED, est + ref, "a response and a well-powered null")
    if est or ref:
        return _decided(claim, CORRECT if est else INCORRECT, est or ref)
    return _fallback(claim, obs)


def _reporter(claim: Claim, obs: list[Observation]) -> Verdict:
    v = claim.value
    rep = [o for o in obs if reading(o.kind, ACTIVITY, v) in (ESTABLISHES, REFUTES, READ_INADEQUATE)]
    scoped = _in(rep, claim.cell)
    pool = scoped or rep
    est = [o for o in pool if reading(o.kind, ACTIVITY, v) == ESTABLISHES]
    ref = [o for o in pool if reading(o.kind, ACTIVITY, v) == REFUTES]
    inad = [o for o in pool if reading(o.kind, ACTIVITY, v) == READ_INADEQUATE]
    if scoped:
        if est and ref:
            return _decided(claim, UNRESOLVED, est + ref, "two reporters disagree in the stated cell")
        order = ((est, CORRECT), (ref, INCORRECT), (inad, MODEL_INADEQUATE))
    elif v == "active_in_reporter":  # R7: a reporter read it active in at least one cell or tissue
        order = ((est, CORRECT), (inad, MODEL_INADEQUATE), (ref, INCORRECT))
    else:  # R7: a reporter read it inactive everywhere it was tested
        order = ((ref, INCORRECT), (inad, MODEL_INADEQUATE), (est, CORRECT))
    for found, verdict in order:
        if found:
            return _decided(claim, verdict, found)
    return _fallback(claim, obs)


def _states(stated_values: str, value: str) -> bool:
    """An annotation states a value when it names it, or names it with a qualifier (`repeat_derived`
    is stated by `repeat_derived/LINE`; `repeat_derived/SINE` is not)."""
    named = [v.strip() for v in stated_values.split(",") if v.strip()]
    return value in named or value in {v.split("/", 1)[0] for v in named}


def _valued(claim: Claim, obs: list[Observation]) -> Verdict:
    """Origin and molecular role: observations that state a value (annotation, registry)."""
    pool = [
        o
        for o in obs
        if o.kind in VALUE_KINDS
        and o.value
        and reading(o.kind, claim.axis, claim.value) in (ESTABLISHES, REFUTES)
    ]
    est = [o for o in pool if _states(o.value, claim.value)]
    ref = [o for o in pool if not _states(o.value, claim.value)]
    if est and ref:
        return _decided(claim, UNRESOLVED, est + ref, "two annotations disagree")
    if est:
        return _decided(claim, CORRECT, est)
    if ref:
        return _decided(claim, INCORRECT, ref, "the annotation states another value")
    return _fallback(claim, obs)


def verdict_of(claim: Claim, observations: Iterable[Observation]) -> Verdict:
    """One claim against the observations of its element, under `TABLE` and `AXIS_RULE`."""
    if claim.axis not in AXES:
        raise ValueError(f"{claim.axis!r} is not a judged axis ({', '.join(AXES)})")
    if claim.value in NOT_A_CLAIM or "|" in claim.value:
        raise ValueError(f"{claim.value!r} is not a claim: unknown and unchosen alternatives are counted")
    obs = list(observations)
    if claim.axis == TARGET:
        return _target(claim, obs)
    if claim.axis == CONTEXT:
        return _context(claim, obs)
    if claim.axis == ACTIVITY:
        if claim.value in DIRECTION:
            return _direction(claim, obs)
        if claim.value == "no_effect_measured":
            return _no_effect(claim, obs)
        if claim.value in REPORTER:
            return _reporter(claim, obs)
        return _fallback(claim, obs)
    return _valued(claim, obs)


# --- the evidence: holdout's units as observations of an element ------------------------------------
class Evidence:
    """The observations each element has, from the held-out sources the caller may judge with, plus any
    observations of kinds holdout does not load (`extra`, for tests and later sources)."""

    def __init__(
        self, units: dict[str, tuple[ho.Unit, ...]], sources: Iterable[str], extra: Iterable[Observation] = ()
    ) -> None:
        self.sources = list(sources)

        def of(assay: str) -> list[ho.Unit]:
            return [u for s in self.sources if ho.assay_of(s) == assay for u in units[s]]

        self.crispri = ho._Index(of("crispri"))
        self.mpra = {s: ho._Index(units[s]) for s in self.sources if ho.assay_of(s) == "lentimpra"}
        self.vista = ho._Index(of("vista"))
        self.satmut = ho._Index(of("satmut"))
        self.gtex = ho._Index(of("gtex"))
        self.extra = ho._Index(list(extra))
        self._memo: dict[tuple[str, int, int], list[Observation]] = {}

    def _reporter(self, chrom: str, start: int, end: int) -> Iterator[Observation]:
        for s, idx in self.mpra.items():
            vals = [
                t.value
                for t in idx.near(chrom, start, end)
                if t.value is not None and ms.measures(start, end, t.start, t.end)
            ]
            if vals:
                lab = ms.reporter_label(vals)
                kind = {ms.LABEL_ACTIVE: "lentimpra_active", ms.LABEL_SILENT: "lentimpra_inactive"}.get(
                    lab, "lentimpra_tiles_conflict"
                )
                yield Observation(kind, s, cell=s.split(":", 1)[1])

    def at(self, chrom: str, start: int, end: int) -> list[Observation]:
        key = (chrom, start, end)
        hit = self._memo.get(key)
        if hit is not None:
            return hit
        out = [
            Observation(CRISPRI_KIND[u.outcome], u.source, u.gene, u.cell, u.split)
            for u in self.crispri.near(chrom, start, end)
            if ms.measures(start, end, u.start, u.end)
        ]
        out.extend(self._reporter(chrom, start, end))
        for v in self.vista.near(chrom, start, end):
            if ms.measures(start, end, v.start, v.end):
                kind = "vista_positive" if v.outcome == "positive" else "vista_negative"
                out.append(Observation(kind, "vista", cell=VISTA_CONTEXT))
        bases = [b for b in self.satmut.near(chrom, start, end) if b.start >= start and b.end <= end]
        if bases:
            functional = any(b.outcome == "functional" for b in bases)
            out.append(
                Observation("satmut_functional_bases" if functional else "satmut_bases_inert", "satmut")
            )
        for g in self.gtex.near(chrom, start, end):
            kind = "gtex_associated" if g.outcome == "associated" else "gtex_not_associated"
            out.append(Observation(kind, "gtex", g.gene))
        out.extend(self.extra.near(chrom, start, end))
        if len(self._memo) > 4096:
            self._memo.clear()
        self._memo[key] = out
        return out


# --- the three quantities -------------------------------------------------------------------------
class NotSummableError(TypeError):
    """Raised when anything tries to add two quantities: they have different denominators (NO_SUM)."""


@dataclass(frozen=True)
class Share:
    """A count over its own denominator. It refuses addition: target accuracy, role accuracy and
    coverage are never combined, and neither are two axes."""

    numerator: int
    denominator: int
    counts: str  # what the numerator counts
    over: str  # what the denominator counts
    beside: dict[str, int] = field(default_factory=dict)

    @property
    def value(self) -> float | None:
        return round(self.numerator / self.denominator, 4) if self.denominator else None

    def __add__(self, other: Any) -> Any:
        raise NotSummableError(NO_SUM)

    __radd__ = __add__

    def to_dict(self) -> dict[str, Any]:
        return {
            "numerator": self.numerator,
            "denominator": self.denominator,
            "value": self.value,
            "counts": self.counts,
            "over": self.over,
            "beside": dict(self.beside),
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> Share:
        return cls(d["numerator"], d["denominator"], d["counts"], d["over"], dict(d.get("beside", {})))


@dataclass
class AxisTally:
    claims: int = 0
    verdicts: Counter = field(default_factory=Counter)
    not_judged_by_reason: Counter = field(default_factory=Counter)
    not_judged_detail: Counter = field(default_factory=Counter)
    decided_by_kind: Counter = field(default_factory=Counter)
    touched_by_kind: Counter = field(default_factory=Counter)
    by_block: Counter = field(default_factory=Counter)
    refutable: int = 0
    decided_with_heldout_file_pair: int = 0

    def add(self, v: Verdict, touched: set[str]) -> None:
        self.claims += 1
        self.verdicts[v.verdict] += 1
        self.by_block[v.claim.block] += 1
        self.refutable += v.refutable
        for k in touched:
            self.touched_by_kind[k] += 1
        if v.verdict == NOT_JUDGED:
            self.not_judged_by_reason[v.reason] += 1
            if v.reason in (OUTSIDE_SCOPE, CONDITION_UNMET) and v.detail:
                self.not_judged_detail[v.detail] += 1
        else:
            for k in {o.kind for o in v.deciding}:
                self.decided_by_kind[k] += 1
            self.decided_with_heldout_file_pair += any(o.split == ms.HELDOUT for o in v.deciding)

    def judged(self) -> int:
        return sum(self.verdicts[x] for x in JUDGED)

    def accuracy(self, axis: str) -> Share:
        beside = {UNRESOLVED: self.verdicts[UNRESOLVED], MODEL_INADEQUATE: self.verdicts[MODEL_INADEQUATE]}
        if axis in (TARGET, CONTEXT):
            beside["refutable"] = self.refutable
        return Share(
            self.verdicts[CORRECT],
            self.verdicts[CORRECT] + self.verdicts[INCORRECT],
            f"{axis} claims established",
            f"{axis} claims established or refuted",
            beside,
        )

    def coverage(self, axis: str) -> Share:
        return Share(
            self.judged(),
            self.claims,
            f"{axis} claims an observation the table allows establishes, refutes or cannot read",
            f"{axis} claims stated with one definite value",
            {NOT_JUDGED: self.verdicts[NOT_JUDGED], **{r: self.not_judged_by_reason[r] for r in REASONS}},
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "claims": self.claims,
            "verdicts": {v: self.verdicts[v] for v in VERDICTS},
            "not_judged_by_reason": {r: self.not_judged_by_reason[r] for r in REASONS},
            "not_judged_detail": dict(sorted(self.not_judged_detail.items())),
            "decided_by_kind": dict(sorted(self.decided_by_kind.items())),
            "touched_by_kind": dict(sorted(self.touched_by_kind.items())),
            "by_block": dict(sorted(self.by_block.items())),
            "refutable": self.refutable,
            "decided_with_heldout_file_pair": self.decided_with_heldout_file_pair,
        }

    @classmethod
    def from_dict(cls, t: dict[str, Any]) -> AxisTally:
        return cls(
            claims=t["claims"],
            verdicts=Counter(t["verdicts"]),
            not_judged_by_reason=Counter(t["not_judged_by_reason"]),
            not_judged_detail=Counter(t["not_judged_detail"]),
            decided_by_kind=Counter(t["decided_by_kind"]),
            touched_by_kind=Counter(t["touched_by_kind"]),
            by_block=Counter(t["by_block"]),
            refutable=t["refutable"],
            decided_with_heldout_file_pair=t["decided_with_heldout_file_pair"],
        )


@dataclass
class Report:
    """What `judge` returns: per axis tallies, the judged claims, and the quantities kept apart."""

    labelling: str
    reads: list[str]
    built_without: str | None
    sources: list[str]
    axes: dict[str, AxisTally]
    judged: list[Verdict]
    not_claims: dict[str, dict[str, int]] = field(default_factory=dict)

    def quantities(self) -> dict[str, dict[str, Share]]:
        """The three quantities, each per axis, never combined (NO_SUM)."""
        return {q: {a: self._share(q, a) for a in axes} for q, axes in QUANTITIES.items()}

    def also_reported(self) -> dict[str, dict[str, Share]]:
        return {q: {a: self.axes[a].accuracy(a) for a in axes} for q, axes in ALSO_REPORTED.items()}

    def _share(self, quantity: str, axis: str) -> Share:
        t = self.axes[axis]
        return t.coverage(axis) if quantity == "coverage" else t.accuracy(axis)

    def to_dict(self) -> dict[str, Any]:
        def shares(block: dict[str, dict[str, Share]]) -> dict[str, Any]:
            return {q: {a: s.to_dict() for a, s in v.items()} for q, v in block.items()}

        return {
            "labelling": self.labelling,
            "reads": self.reads,
            "built_without": self.built_without,
            "sources": self.sources,
            "quantities": shares(self.quantities()),
            "also_reported": shares(self.also_reported()),
            "no_sum": NO_SUM,
            "axes": {a: t.to_dict() for a, t in self.axes.items()},
            "not_claims": self.not_claims,
            "judged": [v.to_dict() for v in self.judged],
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> Report:
        return cls(
            d["labelling"],
            list(d["reads"]),
            d["built_without"],
            list(d["sources"]),
            {a: AxisTally.from_dict(t) for a, t in d["axes"].items()},
            [Verdict.from_dict(v) for v in d["judged"]],
            d.get("not_claims", {}),
        )


# --- the call ------------------------------------------------------------------------------------
def judge(
    claims: Iterable[Claim],
    labels: ho.Labels,
    *,
    units: dict[str, tuple[ho.Unit, ...]] | None = None,
    sources: Iterable[str] | None = None,
    extra: Iterable[Observation] = (),
    references: dict[str, frozenset[str]] | None = None,
    not_claims: dict[str, Counter] | None = None,
) -> Report:
    """Judge a labelling's claims against held-out sources, under holdout's split discipline.

    `labels` carries the labelling's provenance (`reads`, `built_without`); its `predict` is not called.
    Every source in `sources` (default: all holdout sources present) is passed to
    `holdout.check_provenance` first, and its LeakError propagates: a labelling built from
    `holdout.evidence(without=S)` is judged with `sources=[S]`. An `extra` observation whose source the
    labelling read is refused with the same error. `not_claims` (counts of what the labelling states
    that is not a claim) is read after the claims are consumed, so a generator may fill it."""
    units = units if units is not None else ho.all_units()
    srcs = list(sources) if sources is not None else ho.sources(units)
    for s in srcs:
        if s not in units:
            raise KeyError(f"unknown source {s!r}")
        ho.check_provenance(labels, s, references)
    extra = list(extra)
    for o in extra:
        if o.kind in LOADED:
            raise ValueError(f"{o.kind} comes from holdout's sources, never from `extra`")
        if o.source in labels.reads:
            raise ho.LeakError(f"{labels.name} read {o.source!r}, so it cannot be judged by it")
    ev = Evidence(units, srcs, extra)
    tallies = {a: AxisTally() for a in AXES}
    judged: list[Verdict] = []
    for c in claims:
        obs = ev.at(c.chrom, c.start, c.end)
        v = verdict_of(c, obs)
        tallies[c.axis].add(v, {o.kind for o in touching(c, obs)})
        if v.verdict != NOT_JUDGED:
            judged.append(v)
    return Report(
        labels.name,
        sorted(labels.reads),
        labels.built_without,
        srcs,
        tallies,
        judged,
        {a: dict(sorted(n.items())) for a, n in (not_claims or {}).items()},
    )


# --- the unchanged compiled labels as claims -------------------------------------------------------
#: what the compiled programs' stated labels read: none of it is a holdout source
UNCHANGED_READS = frozenset(
    {
        ho.MODEL,
        "encode_ccre_registry",
        "repeatmasker",
        "segmental_duplications",
        "phylop_constraint",
        "gencode_v50",
        "unknown_block_classes",
    }
)
UNCHANGED_SCOPE = (
    "the compiled programs as they stand (data/knowledge/compiled/noncoding_*.bio): every predicted "
    "`element` block (origin and molecular role), every `region` block (origin and molecular role), and "
    "every `rule` with `evidence: predicted` (target, activity and context). The measured twins "
    "(`*_measured`) are the evidence, never the labelling, and are not read. A `,`-joined axis value gives "
    "one claim per value; a `|` group (unchosen alternatives) and `unknown` or `unassigned` are counted "
    "as not claims"
)


def _no_prediction(u: ho.Unit, endpoint: str) -> float | None:
    raise NotImplementedError("the S4 scorer judges claims; it never asks a labelling to predict a unit")


def unchanged_labels() -> ho.Labels:
    return ho.Labels(ho.UNCHANGED, _no_prediction, UNCHANGED_READS, note=UNCHANGED_SCOPE)


def _groups(value: str) -> list[list[str]]:
    return [g.split("|") for g in (x.strip() for x in value.split(",")) if g]


def compiled_claims(
    compiled: Path = ho.COMPILED_DIR,
    chroms: Iterable[str] | None = None,
    not_claims: dict[str, Counter] | None = None,
) -> Iterator[Claim]:
    """The unchanged labels' claims, streamed program by program (`UNCHANGED_SCOPE`). `not_claims`, when
    given, counts per axis what is stated but is not a claim."""
    nc = not_claims if not_claims is not None else defaultdict(Counter)
    only = set(chroms) if chroms is not None else None
    for p in ho.compiled_programs(compiled):
        if only is not None and p.stem.removeprefix("noncoding_") not in only:
            continue
        loci: dict[str, tuple[str, int, int]] = {}
        for kind, name, f in ho.program_blocks(p):
            if kind == "rule":
                if f["evidence"] != "predicted":
                    continue
                if name not in loci:
                    nc["rules"]["predicted_rule_without_its_element"] += 1
                    continue
                c, s, e = loci[name]
                gene, cell = f["gene"], f["cell"]
                yield Claim(name, c, s, e, TARGET, gene, gene, cell)
                act = "activates_target" if f["action"] == "activates" else "represses_target"
                yield Claim(name, c, s, e, ACTIVITY, act, gene, cell)
                if stated(cell):
                    yield Claim(name, c, s, e, CONTEXT, cell, gene, cell)
                else:
                    nc[CONTEXT]["unknown"] += 1
                continue
            if name.endswith("_measured"):
                continue
            m = ho.LOCUS_RE.match(f.get("locus", ""))
            if not m:
                nc["blocks"]["no_locus"] += 1
                continue
            c, s, e = m.group(1), int(m.group(2)), int(m.group(3))
            if kind == "element":
                loci[name] = (c, s, e)
            for axis in (ORIGIN, ROLE):
                for g in _groups(f.get(axis, "")) or [["unknown"]]:
                    if len(g) == 1 and g[0] not in NOT_A_CLAIM:
                        yield Claim(name, c, s, e, axis, g[0], block=kind)
                    else:
                        nc[axis]["unknown" if len(g) == 1 else "unchosen_alternatives"] += 1
            if kind == "region":
                for axis, key in ((ACTIVITY, "activity"), (TARGET, "target_relation")):
                    val = f.get(key, "unknown")
                    nc[axis][f"region_{val}" if val in NOT_A_CLAIM else "region_value_not_read"] += 1


def table_rows() -> list[dict[str, Any]]:
    """The table as rows, for a result file and the documented table."""
    return [
        {
            "observation": k,
            "text": OBSERVATION_TEXT[k],
            "loaded": k in LOADED,
            "not_loaded_because": NOT_LOADED.get(k),
            "cells": {
                a: {
                    "code": TABLE[k][a].code(),
                    "kind": TABLE[k][a].kind,
                    "establishes": list(TABLE[k][a].establishes),
                    "refutes": list(TABLE[k][a].refutes),
                    "cannot_be_read": list(TABLE[k][a].inadequate),
                    "suggests": list(TABLE[k][a].suggests),
                    "why": TABLE[k][a].why,
                }
                for a in AXES
            },
        }
        for k in OBSERVATIONS
    ]


def registration() -> dict[str, Any]:
    """The registered constants, as a result file states them."""
    return {
        "registered": REGISTERED,
        "status": STATUS,
        "reuse": REUSE,
        "axes": AXIS_TEXT,
        "never_judged": NEVER_JUDGED,
        "quantities": {k: list(v) for k, v in QUANTITIES.items()},
        "also_reported": {k: list(v) for k, v in ALSO_REPORTED.items()},
        "no_sum": NO_SUM,
        "axis_rules": AXIS_RULE,
        "match_rule": MATCH_RULE,
        "verdicts": list(VERDICTS),
        "explanations": {k: list(v) for k, v in EXPLANATIONS.items()},
        "axis_explanations": {f"{a}:{v}": list(x) for (a, v), x in AXIS_EXPLANATIONS.items()},
        "not_judged_reasons_in_order": list(REASONS),
        "observation_model_inadequate": OBSERVATION_MODEL_INADEQUATE,
        "may_be_called": MAY_BE_CALLED,
        "may_not_be_called": list(MAY_NOT_BE_CALLED),
        "expected": EXPECTED,
        "unchanged_scope": UNCHANGED_SCOPE,
        "unchanged_reads": sorted(UNCHANGED_READS),
        "table": table_rows(),
        "holdout_registration": ho.REGISTERED,
    }
