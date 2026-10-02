# SPDX-License-Identifier: AGPL-3.0-or-later
"""Source census of the rule numbers of `data/demo/gastrulation.bio`.

For every number the program's `rule` declarations carry - `strength`, `threshold`, `hill`, one slot
each - this asks the literature a single question: is there a MEASURED value with a QUOTED source,
or is the number unsourced? The answer is allowed to be "unsourced" for all of them; that is a
finding and not a failure, and nothing here changes a number in any program.

Three axes are kept apart, because the predecessor census (`number_provenance`, registration
cd43fce) showed that conflating them is how a weakness hides:

1. `class` - what the LITERATURE holds for that quantity: `M_measured_quoted`,
   `F_fitted_or_modelled`, `E_existence_only`, `U_unsourced`.
2. `cited_source_attribution` - what the PROGRAM's own citation does for the declaration it is
   attached to. A citation can be correct to the page and still not mention one of the two factors.
3. `commensurable` - whether a value in the literature's units could be written into this field at
   all. The GRN runtime is in arbitrary concentration units, and it REFUSES a rule whose threshold
   carries a concentration unit unless the program is located in a compartment with a volume
   (`genomeos/runtime/grn.py`, the `molar` check). A measured EC50 in nM is therefore not merely
   absent from this program: as the program is written it cannot be expressed in it.

The rule that decides the headline count: a parameter obtained by FITTING a model to data is not a
measurement of that parameter. It is recorded, named as fitted, and never counted as measured. The
one traceable strength/threshold/hill in the corpus, in `data/demo/lateral_inhibition.bio`, is of
exactly this kind - its source is Collier et al. 1996's own dimensionless Hill functions, a
parameterised theoretical model - so it is not a measurement either.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from genomeos.lang import number_provenance as np

#: The one program this census is about.
PROGRAM = "data/demo/gastrulation.bio"

#: The three numbers a `rule` declaration carries, one slot each. Every rule of the program writes
#: all three explicitly, so the population is 3 x (the number of rules) and no slot is a default.
RULE_FIELDS = ("strength", "threshold", "hill")

#: What each field IS in the runtime, so that it is fixed BEFORE the search what a measurement of it
#: would have to be a measurement of. Transcribed from `genomeos/runtime/grn.py`:
#: transcription = basal + max_rate * A * R, A = mean_i s_i H(x_i), R = prod_j (1 - s_j H(x_j)),
#: H(x) = x^n / (K^n + x^n).
QUANTITY_DEFINITIONS = {
    "strength": (
        "s, dimensionless. For an inhibitor, R = prod (1 - s H), so s is the FRACTION of the "
        "target's max_rate that this regulator alone removes at saturating regulator level: s = 1 "
        "is complete shutoff, a switch and not a brake. For an activator under the mean rule it is "
        "the fraction of max_rate this regulator alone can drive at saturation. A measurement of it "
        "would be the saturating fold-change of the TARGET'S TRANSCRIPTION RATE caused by that one "
        "regulator in the cell being modelled, with the regulator at saturation."
    ),
    "threshold": (
        "K, in the same units as the regulator's protein level, which this runtime states are "
        "ARBITRARY concentration units. It is the regulator level at which the Hill term is one "
        "half. A measurement of it would be the regulator concentration at half-maximal target "
        "response (an EC50 or IC50 read on the regulator's own concentration), in units the program "
        "can convert to its own - and this program provides no such conversion."
    ),
    "hill": (
        "n, the exponent of H(x) = x^n / (K^n + x^n): the sharpness of the response. Values for n "
        "in the literature are obtained by FITTING a Hill function to a dose-response curve, or by "
        "choosing n inside a model. A fit is not a measurement of n, so a literature Hill "
        "coefficient is classed F_fitted_or_modelled by construction. This consequence is declared "
        "here, before the search, so that it cannot later look like a convenience."
    ),
}

#: Axis 1: what the literature holds for that quantity. Exactly one per slot, first that applies.
CLASSES: dict[str, str] = {
    "M_measured_quoted": (
        "a source that was FETCHED AND READ in this lane reports a DIRECT MEASUREMENT of the "
        "quantity this field is (as QUANTITY_DEFINITIONS fixes it), the measured value is quoted "
        "from that source, and the value is in units commensurable with the field. Not a value "
        "obtained by fitting a model, not a value chosen inside a model, not a value for a "
        "different quantity that shares the field's name."
    ),
    "F_fitted_or_modelled": (
        "a fetched source gives a value for this quantity, but by FITTING a model to data, or by "
        "assuming or choosing it inside a model. It is recorded with its value and its source and "
        "is NOT a measurement of the parameter. Hill coefficients reach this class by construction."
    ),
    "E_existence_only": (
        "a fetched source supports the EXISTENCE of the interaction the rule states, naming both "
        "factors, but states no value for this quantity by either route. A citation that supports "
        "an interaction is not a citation for its strength."
    ),
    "U_unsourced": (
        "no source, fetched or cited, gives a value for this quantity by either route, and the "
        "slot's own citation does not even establish the interaction. Also the class when the "
        "declaration cites nothing at all."
    ),
}

#: The order the classes are tried in. The first that applies wins.
CASCADE = ("M_measured_quoted", "F_fitted_or_modelled", "E_existence_only", "U_unsourced")

#: Axis 2: what the program's OWN citation does for the declaration it sits on, decided only from
#: what was fetched. This is the axis on which the predecessor lane was weakest: two rules were
#: attributed to a paper whose abstract mentions neither factor.
ATTRIBUTION = {
    "supports_declaration": (
        "what was fetched names both factors of the rule and supports a regulatory relation of the "
        "stated sign between them"
    ),
    "supports_sign_in_another_system": (
        "what was fetched names both factors and supports the sign, but in an organism, cell type "
        "or stage other than the one the program models"
    ),
    "names_one_factor_only": "what was fetched names one of the two factors and not the other",
    "names_neither_factor": "what was fetched names neither factor of the rule",
    "no_source_cited": (
        "the declaration's evidence kind is `inferred` or `none`, which BioLang's own grammar "
        "defines as not measured, so there is no citation to check"
    ),
}

#: Axis 3: whether a literature value could be written into this field as the program stands.
COMMENSURABILITY = {
    "dimensionless_and_expressible": (
        "the field is dimensionless (strength, hill), so a literature value needs no conversion to "
        "be written into it"
    ),
    "arbitrary_units_no_conversion": (
        "the field is in the runtime's arbitrary concentration units and this program supplies no "
        "conversion to physical units: its `translation_rate` and `mrna_half_life` are "
        '`evidence: inferred "dimensionless"`, and the NODAL level the two Nodal rules are read '
        "against is set by `run_gastrulation`'s unsourced `nodal_max = 6.0`. A measured EC50 "
        "therefore has no scale to be written on."
    ),
    "molar_threshold_refused_by_the_runtime": (
        "writing the measured unit onto the field is refused at run time: `NetworkRuntime.__init__` "
        "raises when any active rule has a `threshold_unit`, because a concentration threshold "
        "needs a compartment with an `absolute_volume`, and this program is not located. The "
        "language parses `threshold: 3.0 nM`; this program cannot run with it."
    ),
}

#: Each finding must carry these keys, so that no entry can be a bare verdict.
REQUIRED_FINDING_KEYS = (
    "rule",
    "field",
    "written_value",
    "class",
    "cited_source",
    "fetched_url",
    "fetched_what",
    "source_claims",
    "source_does_not_claim",
    "cited_source_attribution",
    "commensurable",
    "searches",
)

#: What will be searched for, fixed before the search so the search cannot be narrowed to fit the
#: answer. Every entry is a query, not a conclusion.
SEARCH_PLAN = (
    "each of the four works the program cites, fetched at its publisher or PubMed record: Conlon "
    "et al. 1994 Development 120:1919; Kanai-Azuma et al. 2002 Development 129:2367; Thomson et al. "
    "2011 Cell 145:875; Lolas et al. 2014 PNAS 111:4478",
    "for each rule, a query naming both factors together with each of `EC50`, `dose-response`, "
    "`Hill coefficient`, `half-maximal` and `fold repression`",
    "a query for a measured or fitted dose-response of NODAL or Activin signalling on T/Brachyury "
    "and on SOX17, since those are the two rules whose regulator is the gradient",
    "a query for any published parameterised model of the germ-layer switch that states a strength, "
    "threshold or Hill coefficient for these interactions, which would be F and not M",
)

#: Reported, not acted on.
REPORTS = (
    "`data/demo/gastrulation.bio` line 3 says `Mutual repression makes the choice sharp`, and that "
    "is false of the program as written: SOX17 represses TBXT and SOX2 and is repressed by nothing. "
    "`data/demo` is shared, so this is a report and the line is not edited.",
    "no number in any `.bio` file is changed by this lane. The area has already refused to tune two "
    "free numbers to unsourced targets (8bb9123, 19423bd) and that refusal is inherited. The "
    "Tbxt-SOX17 sign correction (5cbce26) stands.",
    "the headline count is reported twice over the same 21 slots: strictly, counting only a "
    "measurement of the quantity the field is in commensurable units; and loosely, counting a "
    "measurement of that quantity in any organism or cell type. Both denominators are 21.",
)

#: What is NOT claimed, whichever way the count falls.
WHICHEVER_WAY_IT_FALLS = (
    "if none of the 21 slots reaches M_measured_quoted, that is the result and it is stated plainly: "
    "the program's rule numbers are unsourced. No source is stretched to fill a cell, no slot is "
    "promoted because its paper is about the right pathway, and no number is adjusted toward any "
    "value found. If some slot does reach M, the finding says so with the quoted value and does NOT "
    "write it into the program: a census does not edit its population."
)

EXCLUSIONS = (
    "the program's `gene`, `protein` and `param` numbers are out of scope: this census is the 21 "
    "rule numbers only. The predecessor census covers all 42 scoped slots of the file.",
    "no judgement is passed on whether a rule is true, only on whether its numbers have a source.",
    "no accuracy, proportion or fate figure is computed, and no run of the program is made for this census.",
)


def slots(path: str | Path = PROGRAM) -> list[dict[str, Any]]:
    """Every rule number of the program, read from the program itself.

    The population is enumerated from the file and never asserted, so the count cannot drift from
    what the file holds.
    """
    out: list[dict[str, Any]] = []
    for d in np.read_program(path):
        if d.kind != "rule":
            continue
        kind, text = np.evidence_kind_and_text(d.props)
        for f in RULE_FIELDS:
            written = d.props.get(f)
            out.append(
                {
                    "rule": d.name,
                    "field": f,
                    "line": d.line,
                    "written_value": written,
                    "evidence_kind": kind,
                    "cited_source": text,
                    "confidence": d.props.get("confidence"),
                }
            )
    return out


#: The findings table, filled in by hand after every citation has been fetched and read. It is EMPTY
#: at registration, and a test asserts that the registration payload assigns no class.
FINDINGS: dict[tuple[str, str], dict[str, Any]] = {}

#: Every source reached, one entry per URL, with what was retrieved from it. Empty at registration.
SOURCES_FETCHED: tuple[dict[str, Any], ...] = ()

#: What the program's own citations do and do not establish, per citation. Empty at registration.
ATTRIBUTION_FINDINGS: tuple[dict[str, Any], ...] = ()

#: Anything the code or the literature contradicts in what this lane was told or assumed. Empty at
#: registration: a contradiction may only be recorded after it is found.
CONTRADICTIONS: tuple[dict[str, Any], ...] = ()

#: The reading, in the words the result will carry. Empty at registration.
VERDICT = ""


def registration() -> dict[str, Any]:
    """The whole question, population, classes and search plan, before anything is searched."""
    pop = slots()
    return {
        "question": (
            "For each of the rule numbers of data/demo/gastrulation.bio - strength, threshold and "
            "hill on every rule - is there a MEASURED value with a QUOTED source, or is the number "
            "unsourced? A value obtained by fitting a model is never counted as measured."
        ),
        "why_now": (
            "The predecessor census (registration cd43fce, result 00e7028, row 2e1f676) established "
            "by count, 0 slots unclassified, that this file has 42 scoped slots and ZERO sourced "
            "numbers, and that corpus-wide not one strength, threshold or hill outside "
            "lateral_inhibition.bio is sourced. That is why three different sourced repressors of "
            "SOX17 gave answers identical to the digit: the module's own inhibiting numbers make "
            "any added inhibitor a switch and not a brake. This lane asks the complementary "
            "question, of the literature rather than of the file."
        ),
        "program": PROGRAM,
        "rule_fields": list(RULE_FIELDS),
        "quantity_definitions": dict(QUANTITY_DEFINITIONS),
        "population_size": len(pop),
        "population": [{k: s[k] for k in ("rule", "field", "line", "written_value")} for s in pop],
        "classes": dict(CLASSES),
        "cascade": list(CASCADE),
        "cited_source_attribution_values": dict(ATTRIBUTION),
        "commensurability_values": dict(COMMENSURABILITY),
        "required_finding_keys": list(REQUIRED_FINDING_KEYS),
        "search_plan": list(SEARCH_PLAN),
        "every_citation_is_fetched": (
            "Nothing is taken from memory. Each finding records the URL fetched and what was "
            "retrieved from it - at minimum an abstract actually retrieved - and states separately "
            "what the source claims and what it does not. A source that could not be fetched is "
            "recorded as not fetched and supports nothing."
        ),
        "fitted_is_never_measured": (
            "A parameter estimated by fitting a model to data is not a measurement of that "
            "parameter. Such a value goes to F_fitted_or_modelled, is named as fitted, and is never "
            "counted in the measured total. lateral_inhibition.bio's threshold and hill, the only "
            "traceable ones in the corpus, are of this kind: Collier et al. 1996 is itself a "
            "parameterised theoretical model."
        ),
        "denominator": (
            "21 rule numbers: 3 fields x 7 rule declarations, enumerated from the program. Both the "
            "strict and the loose measured counts are reported against this same 21."
        ),
        "findings_at_registration": len(FINDINGS),
        "sources_fetched_at_registration": len(SOURCES_FETCHED),
        "verdict_at_registration": VERDICT,
        "reports": list(REPORTS),
        "whichever_way_it_falls": WHICHEVER_WAY_IT_FALLS,
        "exclusions": list(EXCLUSIONS),
    }
