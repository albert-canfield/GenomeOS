# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The defaults a compiled BioLang program leans on, written into the program that leans on them.

ADOPTED 2026-10-03 on Albert's decision, item (7): "adopt the header fix, not_assessed, visible
defaults and per-slot evidence tier; defer cell provenance". The registered proposal this adopts is
`genomeos/provenance/rule_evidence_tier.py` (design `0d4247b`, artefact `709aebc`); its vocabulary is
IMPORTED here and not re-spelled, so there is one axis in the language and not two.

Three things are made visible, and all three were invisible before:

1. **The numeric defaults.** A compiled rule line states `strength:` and nothing else of the Hill
   function. `threshold` and `hill` are then supplied by `genomeos/ir/model.py`'s `Rule` dataclass -
   1.0 and 2.0 - and a reader of the program cannot see that. A default that is applied and not
   printed is a number in the dynamics that no line of the program carries. `declaration_lines()`
   prints them, and it reads them off `dataclasses.fields(Rule)` rather than transcribing them, so
   the declaration cannot drift away from the engine while both still look right.

2. **The evidence tier, PER SLOT.** One rule carries three numbers and the census found tiers
   DIFFERING WITHIN one rule across them (`4d4003d`), so the mark is keyed by field name and never
   by rule. A compiled program names each of `rule_evidence_tier.TIERED_FIELDS` with its own value.

3. **`not_assessed` as the absent value, never `unsourced`.** A slot nobody has assessed and a slot
   assessed and found to rest on nothing are different facts. `U_unsourced` is a finding - a search
   was made and came back empty - and the compiler has made no search, so it may not write it. The
   absent value is `rule_evidence_tier.NOT_ASSESSED`, which is the absence of a tier and not the
   lowest one; `rule_evidence_tier.rank` refuses to order it against a tier, and that refusal is
   what keeps the two apart downstream.

The same distinction governs `context_evidence`: the empty string stored on a `Rule` whose program
stated none READS `not_assessed` here, and `not_assessed` is not `not_assessable` - the latter is a
reading that was attempted and could not be taken, and it carries a reason. A guard refuses a
program that gives `not_assessed` a reason, because a reason is the mark of an attempt.

WHAT THIS DOES NOT DO. It gates nothing, it changes no number, and it does not add a per-rule
`evidence_tier:` property to the grammar: writing a mark onto every compiled rule LINE would move
all 5,176 rule lines of chr21 on the next cut, and a re-cut of a published program is Albert's
decision. The declaration is a program-level statement, so the rule lines of a cut are untouched by
it. Cell provenance (`genomeos/provenance/rule_cell_provenance.py`) is deferred and is not read here.

AND ONE PIECE WAS BUILT AND WITHDRAWN, named rather than dropped. Writing `context_evidence:
not_assessed` onto the rule line where no reading was taken is blocked by two artefacts that belong
to another lane: `tests/test_context_evidence.py`'s
`test_the_state_is_added_beside_the_rule_and_deletes_nothing` pins that a program compiled with no
reader carries no `context_evidence` at all, and `attribution.context_evidence.parse_value` refuses
any value not beginning with one of its three states, which `not_assessed` does not. Neither was
weakened. What is delivered instead is the READING - `context_evidence_value` and the header's own
sentence - which says the same thing from the program alone without touching a peer's pin.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Iterable, Mapping
from typing import Any

from genomeos.ir.model import Rule
from genomeos.provenance import rule_evidence_tier as ret

#: The absent value on the evidence-tier axis, the proposal's own constant and not a second copy.
ABSENT: str = ret.NOT_ASSESSED

#: The three numbers one rule carries, one mark each. The same tuple object as the census's.
TIERED_FIELDS: tuple[str, ...] = ret.TIERED_FIELDS

#: The four tiers. `ABSENT` is not among them and that is the whole point of it.
TIERS: tuple[str, ...] = ret.TIERS

#: Every value a mark may hold: the four tiers plus the absence of one.
AXIS_VALUES: tuple[str, ...] = ret.AXIS_VALUES

#: The `Rule` fields whose default the engine applies where a program line is silent, and which a
#: compiled program must therefore print. `strength` is included although every compiled rule states
#: it: a HAND-authored rule may leave it out, and the declaration speaks for the language and not
#: only for what one generator happens to emit.
DECLARED_NUMBERS: tuple[str, ...] = ("strength", "threshold", "hill")

#: The non-numeric field whose default is the absent reading rather than a value.
DECLARED_ABSENT_FIELD = "context_evidence"

#: Every field named in the declaration, numeric and not.
DECLARED: tuple[str, ...] = (*DECLARED_NUMBERS, DECLARED_ABSENT_FIELD)

#: The word the absent value may never be spelled with. `U_unsourced` is a finding about a search
#: that was made; the compiler has made none, so a compiled program that writes it is claiming one.
FORBIDDEN_ABSENT_SPELLING = "unsourced"

#: The marker that opens the declaration block, so a guard can find it and a reader can see where it
#: begins. Kept here because the guard and the writer must not hold two spellings of it.
BLOCK_MARKER = "# ---- the defaults this program leans on, printed because an applied default that is"


def rule_defaults(rule_type: type = Rule) -> dict[str, Any]:
    """The default each declared field gets, read off the dataclass and never transcribed.

    `rule_type` is injected so a test can show the reading - and the guard built on it - refuse a
    dataclass whose fields have moved, without mutating the live `Rule` that every parse uses.
    """
    if not dataclasses.is_dataclass(rule_type):
        raise ValueError(f"{rule_type!r} is not a dataclass, so it declares no defaults to read")
    by_name = {f.name: f for f in dataclasses.fields(rule_type)}
    missing = [n for n in DECLARED if n not in by_name]
    if missing:
        raise ValueError(
            f"{rule_type.__name__} has no field {missing}: this module declares a default for a "
            "field the engine no longer has, so the declaration would print a number nothing uses"
        )
    out: dict[str, Any] = {}
    for name in DECLARED:
        default = by_name[name].default
        if default is dataclasses.MISSING:
            raise ValueError(
                f"{rule_type.__name__}.{name} has no default, so no default of it can be declared"
            )
        out[name] = default
    return out


def context_evidence_value(stored: str) -> str:
    """What a rule's `context_evidence` READS. The empty string reads `not_assessed`, never "".

    The empty string is what the field holds when no chromatin reading was taken at all. Read as
    "" it is indistinguishable from a field this version of the language does not have; read as
    `not_assessed` it says the one true thing - nobody looked - and stays distinct from
    `not_assessable`, a reading that was attempted and could not be taken.
    """
    text = (stored or "").strip()
    return text or ABSENT


def slot_tiers(rule: Any, tiered_fields: Iterable[str] = TIERED_FIELDS) -> dict[str, str]:
    """The evidence tier of each of a rule's numbers, one mark per number, absent reading `ABSENT`.

    Duck-typed on an optional `evidence_tier` mapping, so this works on today's `Rule`, which has
    none and therefore reads `not_assessed` on every slot, and on a later one that carries marks
    without this module changing. A value off the axis is refused rather than counted: an unknown
    word in a mark is not a quieter kind of unassessed.
    """
    marks = getattr(rule, "evidence_tier", None) or {}
    if not isinstance(marks, Mapping):
        raise ValueError(f"evidence_tier must be a mapping of field to tier, got {type(marks).__name__}")
    fields = tuple(tiered_fields)
    out: dict[str, str] = {}
    for name in fields:
        value = str(marks.get(name, ABSENT))
        if value not in AXIS_VALUES:
            raise ValueError(f"{name}: {value!r} is not on the axis; the axis is {list(AXIS_VALUES)}")
        out[name] = value
    stray = sorted(set(marks) - set(fields))
    if stray:
        raise ValueError(f"evidence_tier marks a field that carries no number: {stray}")
    return out


def _number(value: Any) -> str:
    """A default rendered as the declaration prints it, so writer and guard cannot disagree."""
    if isinstance(value, float) and value == int(value):
        return f"{value:.1f}"
    return str(value)


def declaration_lines(rule_type: type = Rule, tiered_fields: Iterable[str] = TIERED_FIELDS) -> list[str]:
    """The `#` comment block a compiled program carries so its defaults are readable from it alone."""
    defaults = rule_defaults(rule_type)
    return [
        "#",
        BLOCK_MARKER,
        "# invisible to a reader is a number in the dynamics that no line of this program carries.",
        "# Where a rule line is silent the engine supplies, from genomeos/ir/model.py's own Rule:",
        *(f"#   `{n}` is {_number(defaults[n])}" for n in DECLARED_NUMBERS),
        "# Read off that dataclass when this program was written, not transcribed, so the two cannot",
        "# drift apart while both still look right. A rule that states no `confidence:` states none",
        "# (review R4) and that is not a confidence of 0.",
        "#",
        "# ---- evidence tier, one mark per NUMBER and not per rule "
        "(genomeos/provenance/rule_evidence_tier.py)",
        "# One rule carries three numbers and the source census found tiers differing within one rule",
        "# across them (4d4003d), so each number is marked on its own:",
        *(f"#   {n} = {ABSENT}" for n in tiered_fields),
        f"# `{ABSENT}` is the ABSENCE of a tier and not the lowest one. No adjudication has classed any",
        f"# number in this program, which is a different fact from `U_{FORBIDDEN_ABSENT_SPELLING}`: a search",
        "# that was made and came back empty. Nothing here may be read as measured, as fitted or as",
        f"# {FORBIDDEN_ABSENT_SPELLING}; it has not been assessed. The mark gates nothing - no rule is",
        "# refused, skipped or down-weighted by it and no integration moves by a digit.",
        "#",
        f"# A rule that states no `{DECLARED_ABSENT_FIELD}:` reads `{ABSENT}` - no chromatin reading was",
        "# taken at all - which is not `not_assessable`, a reading that was attempted and could not be",
        "# taken and which carries its reason. `not_assessed` carries none, because nothing was tried.",
    ]


def check_declaration(
    text: str, rule_type: type = Rule, tiered_fields: Iterable[str] = TIERED_FIELDS
) -> None:
    """Refuse a compiled program whose declaration is missing, stale, or claims a tier.

    Raises `ValueError` naming what is wrong. Three ways to fail, each one a way a reader could be
    misled: the block is absent, so the defaults are invisible again; a declared default disagrees
    with the dataclass the engine actually applies, so the program prints a number the run does not
    use; or a slot is marked with one of the four tiers, which asserts an adjudication the compiler
    never made - `U_unsourced` above all, which claims a search.
    """
    if BLOCK_MARKER not in text:
        raise ValueError(
            "this program declares none of the defaults it leans on: the block opening "
            f"{BLOCK_MARKER!r} is absent, so `threshold` and `hill` are applied and not printed"
        )
    defaults = rule_defaults(rule_type)
    block = text.split(BLOCK_MARKER, 1)[1]
    header = block.split("\nelement ", 1)[0].split("\nrule ", 1)[0]
    for name in DECLARED_NUMBERS:
        want = f"`{name}` is {_number(defaults[name])}"
        if want not in header:
            raise ValueError(
                f"the declaration does not say {want!r}: the engine applies "
                f"{name}={defaults[name]!r} where a line is silent and this program says otherwise"
            )
    for name in tuple(tiered_fields):
        # the tier claim is tested FIRST: a slot that asserts a tier is also a slot whose `not_assessed`
        # is gone, and naming the missing mark there would report the smaller of the two faults
        for tier in TIERS:
            if f"{name} = {tier}" in header:
                raise ValueError(
                    f"{name} is marked {tier!r}, a tier: a compiled program asserts no adjudication, "
                    f"so the only value it may write is {ABSENT!r}"
                )
        if f"{name} = {ABSENT}" not in header:
            raise ValueError(
                f"the declaration marks no tier for {name!r}: every number the census tiers must "
                f"carry its own mark, and an unmarked number reads as nothing at all"
            )
    if f"{DECLARED_ABSENT_FIELD}:` reads `{ABSENT}`" not in header:
        raise ValueError(
            f"the declaration does not say what an absent `{DECLARED_ABSENT_FIELD}:` reads; it reads "
            f"{ABSENT!r} and must say so, or the empty field is readable as a missing feature"
        )


def check_not_assessed_carries_no_reason(text: str) -> None:
    """Refuse a program that gives `not_assessed` a reason, which would make it `not_assessable`.

    A reason is the mark of an attempt: `not_assessable, no_reader_for_this_cell` says a reading was
    tried and could not be taken. `not_assessed` says nothing was tried, so it has no reason to
    give, and a value that carries one has collapsed the very distinction this module exists for.
    """
    needle = f"{DECLARED_ABSENT_FIELD}: {ABSENT},"
    for n, line in enumerate(text.splitlines(), start=1):
        if needle in line:
            raise ValueError(
                f"line {n}: `{ABSENT}` carries a reason - {line.strip()[:120]!r}. A reason means a "
                "reading was attempted, which is `not_assessable`; nothing was attempted here"
            )


def adoption() -> dict[str, Any]:
    """What was adopted, for a result or a registration to carry without re-spelling it."""
    return {
        "adopted": "rule_evidence_tier, per slot, with not_assessed as the absent value",
        "decided": "Albert, 2026-10-03, item (7)",
        "proposal": "genomeos/provenance/rule_evidence_tier.py (design 0d4247b, artefact 709aebc)",
        "absent_value": ABSENT,
        "absent_value_is_not_a_tier": ABSENT not in TIERS,
        "tiered_fields": list(TIERED_FIELDS),
        "marks_are_per_slot_not_per_rule": True,
        "declared_defaults": {
            k: _number(v) if k in DECLARED_NUMBERS else v for k, v in rule_defaults().items()
        },
        "declared_defaults_are_read_off_the_dataclass": "dataclasses.fields(genomeos.ir.model.Rule)",
        "gates_nothing": True,
        "grammar_unchanged": (
            "no per-rule `evidence_tier:` property: a mark on every rule LINE would move all of a "
            "published program's rule lines on the next cut, which is a re-cut and not this lane's "
            "decision. The declaration is program-level and leaves every rule line byte-identical."
        ),
        "cell_provenance": "deferred, untouched",
    }
