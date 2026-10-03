# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""What every number in a hand-authored BioLang program rests on.

This module holds the population, the scoped fields and the classification of a census of number
provenance, so that the registration (`scripts/number_provenance_register.py`) and the count
(`scripts/number_provenance_census.py`) cannot drift apart. It reads programs and classifies the
slots in them. It changes no strength, threshold, Hill coefficient, basal or max, and it says nothing
about whether any number is right: the census reports what a number RESTS ON, never whether it holds.

Why it exists. Two runs on 2026-10-02 established, by measurement, that the gastrulation module
cannot discriminate mechanisms, and that the cause is its number convention and not its biology:

* `d340e5c` - the program has no mesoderm at all (0.000 at all 41 samples) because nothing in it
  represses SOX17: SOX17's repression factor is exactly 1.00000 while TBXT carries two repressors.
* `53d9936` / `b1e8839` - three different sourced repressors of SOX17 (MESP1 twice, HAND1) each
  restore a band and all three agree to the digit, because an added inhibiting rule takes the
  module's own `threshold 2.0; hill 3`, the carrier settles ten to twenty times above that
  threshold, the repression factor collapses to 1.09e-4, and the repressed gene sits at its basal
  floor. The module's own inhibiting numbers make any added inhibitor a switch, not a brake.

So telling one mechanism from another in this module needs measured kinetics that are not on disk.
The census draws that boundary by counting rather than asserting it.

Four refusals this module inherits and does not revisit: `8bb9123` and `19423bd` refused tuning two
free numbers to meet unsourced targets; the `5cbce26` sign correction stands; and `data/demo` is a
shared demo directory, so a wrong line in a program is reported, never edited.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from genomeos.lang.parser import _parse_blocks

# --- the population -----------------------------------------------------------------------------

#: The hand-authored BioLang corpus, as two directories. Every program a person wrote is in scope;
#: the counts stay separable per file and are never pooled across the two directories in a headline.
POPULATION_DIRS = ("data/demo", "genomeos/std")

#: Named and excluded, with the reason, so the scope cannot be read as "everything".
EXCLUDED_CORPORA = {
    "data/knowledge/compiled/noncoding_chr*.bio": "machine-written from the attribution pipeline, "
    "hundreds of thousands of rules; every number in it comes from one generator, so a per-rule "
    "census would count the generator once per rule. Its provenance is the subject of the separate "
    "context-evidence census (data/results/context_evidence_census*.json), not of this one.",
    "data/organisms/**/*.bio": "machine-written per-organism programs, the same reason.",
    "build/biolang/*.bio": "build artefacts of the engine's own smoke tests, not authored biology.",
}

# --- the scoped fields --------------------------------------------------------------------------

#: A sentinel for "the language supplies no number here": the field is in the schema, the program
#: leaves it out, and the engine holds UNKNOWN or warns rather than inventing a value.
NO_NUMBER = object()


@dataclass(frozen=True)
class Slot:
    """One scoped field of one block kind: its name, its tier and what the language supplies."""

    kind: str
    field: str
    tier: int
    language_default: Any  # a float, or NO_NUMBER
    note: str = ""


#: Tier 1, the dynamics core: every numeric field of a `rule`, a `gene` or a `protein`, and every
#: `param` value. These are exactly the numbers the GRN integrator reads
#: (`genomeos/runtime/grn.py`: basal + max_rate * A * R, with A and R Hill functions of strength,
#: threshold and hill), so they are the numbers a mechanism claim would rest on. `param` is taken
#: whole rather than by a chosen list of names, because choosing which params matter is a judgement
#: and this census makes none.
#:
#: Tier 2, the rest of the declared numbers: timers, stages, compartments, pools, fields, signals,
#: events, regimes and decisions. Counted with the same classifier and reported with its own
#: denominator, never pooled into tier 1's.
SCOPE: tuple[Slot, ...] = (
    Slot("rule", "strength", 1, 1.0),
    Slot("rule", "threshold", 1, 1.0),
    Slot("rule", "hill", 1, 2.0),
    Slot("gene", "max", 1, NO_NUMBER, "omitted: the engine reads 0 and warns `gene_parameter_missing`"),
    Slot("gene", "basal", 1, 0.0),
    Slot("protein", "half_life", 1, NO_NUMBER, "omitted: the IR holds UNKNOWN; no turnover is applied"),
    Slot("protein", "initial", 1, 0.0),
    Slot("param", "value", 1, NO_NUMBER, "a param exists only when written"),
    Slot("timer", "duration", 2, NO_NUMBER, "required by the grammar"),
    Slot("timer", "lengthening", 2, 1.0),
    Slot("stage", "from", 2, 0.0),
    Slot("stage", "to", 2, NO_NUMBER, "omitted: the stage has no end"),
    Slot("compartment", "volume", 2, NO_NUMBER),
    Slot("compartment", "absolute_volume", 2, NO_NUMBER),
    Slot("pool", "size", 2, NO_NUMBER),
    Slot("pool", "regenerates", 2, NO_NUMBER),
    Slot("field", "diffusion", 2, 0.2),
    Slot("field", "decay", 2, 0.01),
    Slot("field", "source", 2, NO_NUMBER),
    Slot("signal", "threshold", 2, 0.0),
    Slot("event", "rate", 2, 0.0),
    Slot("regime", "threshold", 2, 50.0),
    Slot("decision", "after", 2, NO_NUMBER),
    Slot("decision", "priority", 2, 0.0),
)

#: Out of scope, by field, with the reason. Stated so that the denominator is a decision and not an
#: oversight.
OUT_OF_SCOPE_FIELDS = {
    "confidence": "a declaration's self-report about itself, not a quantity the engine integrates; "
    "classifying its provenance would be classifying an opinion about provenance.",
    "locus / exons / cds / members": "genome coordinates. They are numbers, and GENCODE does state "
    "them, so counting them would fill the sourced class with annotation coordinates and hide the "
    "kinetic boundary this census is for.",
    "seed / replicates / space / origin / placement / sense / cell_network / tempo (organism, "
    "regime)": "run control: how long to run, how many times, on what grid. Not biology, so no "
    "source could state them.",
    "assert / test lines": "the bars a program sets itself. `19423bd` already established that "
    "gastrulation's three expected_* proportions cite no source; a bar is not a model number.",
}

# --- the classification, decided before counting ------------------------------------------------

CLASSES: dict[str, str] = {
    "S_quoted": "the field is written; the declaration's own evidence kind is experimental, curated "
    "or predicted; and the written value appears as a literal number in the evidence string once "
    "citation furniture is stripped. Every member is listed with the matched substring and is "
    "adjudicated by hand in the result: a match that falls in furniture is demoted to E.",
    "S_derived": "the same, except the written value is an exact arithmetic transform of a number "
    "the evidence string quotes, AND the transform is written in that same string (for example an "
    "evidence string quoting `f(d) = d^2 / (0.01 + d^2)` for `threshold: 0.1`, where 0.01 = 0.1^2). "
    "Every member is listed with its arithmetic so a reader can check it.",
    "E_existence_only": "the field is written; the declaration cites a source (evidence kind "
    "experimental, curated or predicted with a non-empty source) and that source supports the "
    "declaration's EXISTENCE; and no number equal to the written value is quoted anywhere in the "
    "evidence string. A citation that supports a rule's existence is not a citation for its "
    "strength. This class is expected to be the largest and may not collapse into either neighbour.",
    "A_asserted_names_a_work": "the field is written and the declaration's evidence kind is "
    "`inferred` or `none` - which BioLang's own grammar defines as not measured - but its note names "
    "a publication, a database or a reported range. The number is asserted; the prose is context.",
    "A_asserted_bare": "the field is written, the evidence kind is `inferred`, `none` or absent, and "
    "the note names no work at all.",
    "D_language_default": "the field is NOT written and the language supplies a number for it "
    "(`genomeos/ir/model.py` dataclass defaults: rule strength 1.0, threshold 1.0, hill 2.0, gene "
    "basal 0.0, protein initial 0.0, field diffusion 0.2, decay 0.01, signal threshold 0.0, event "
    "rate 0.0, regime threshold 50.0, stage from 0.0, timer lengthening 1.0, decision priority 0).",
    "D_no_number": "the field is NOT written and the language supplies no number: the IR holds "
    "UNKNOWN, or the engine warns and reads zero. Counted, and reported OUTSIDE the denominator of "
    "numbers, because no number exists to have a provenance.",
    "N_non_numeric": "the field is written but its value is not a number (`regenerates: from "
    "Glycolysis, OXPHOS`). Counted and reported outside the denominator.",
}

#: The order the cascade is applied in. It is a cascade and not a set of tests, so the classes
#: partition the slots: every slot falls in exactly one.
CASCADE = (
    "1. the field is not written -> D_language_default if the language supplies a number for it, "
    "else D_no_number.",
    "2. the field is written but holds no number -> N_non_numeric.",
    "3. the evidence kind is `inferred`, `none` or absent -> A_asserted_names_a_work if the source "
    "or note names a publication, a database or a reported range, else A_asserted_bare. An "
    "`inferred` declaration cannot be a source for its own number whatever it mentions, because "
    "BioLang's grammar defines the kind as not measured; this is decided here and not while "
    "counting.",
    "4. the evidence kind claims a source -> S_quoted if the written value appears literally in the "
    "stripped evidence string, else S_derived if a stated transform of a quoted number, else "
    "E_existence_only.",
)

#: Orthogonal to the classes, reported per class rather than as a class of its own, because a
#: convention and a provenance are different questions and a slot has both.
CONVENTION_TAGS = {
    "equals_language_default": "the written value is exactly what the language would have supplied "
    "had the field been left out.",
    "repeated_in_file": "the same field carries this identical value on at least three declarations "
    "of the same kind in the same file. This is the tag that makes `the module's default "
    "convention` visible: gastrulation's four inhibiting rules all write `threshold: 2.0; hill: 3`, "
    "which is why `53d9936` found three different repressors agreeing to the digit.",
}

REPEATED_IN_FILE_FLOOR = 3

#: What the census reports. Fixed before counting.
REPORTS = (
    "the count per class, per file, per tier, with the denominator stated as written + "
    "D_language_default slots;",
    "the D_no_number and N_non_numeric counts, outside that denominator;",
    "the convention tags per class, per file;",
    "every S_quoted and S_derived member in full, with its matched substring or arithmetic and a "
    "hand adjudication, because the sourced class is the one a reader must be able to audit;",
    "for data/demo/gastrulation.bio alone, the per-rule table: the seven rules by name with each "
    "one's strength, threshold and Hill class;",
    "which cited sources were fetched and which could not be, and what a fetched source does NOT claim;",
    "every slot that could not be classified, with the reason.",
)

#: Claims this census does not make, fixed before counting.
EXCLUSIONS = (
    "no number is changed: no strength, threshold, Hill coefficient, basal or max is edited, and no "
    "program file is written to. The census reads.",
    "no claim about correctness: the census says what a number rests on, never whether it is right. "
    "No validation, no accuracy figure, no verdict on any rule's truth.",
    "no claim that a cited paper is wrong, and none that it is right. Where a source was fetched, "
    "what it states and what it does not state about a NUMBER is reported, and nothing more.",
    "no tuning. `8bb9123` and `19423bd` refused fitting two free numbers to unsourced targets and "
    "this lane inherits the refusal.",
    "S_quoted means the FILE quotes a number equal to the written value. It does not by itself mean "
    "the cited work states that number; only a fetch can say that, and the fetch is reported "
    "separately per source.",
    "no number of this corpus is read as evidence about the machine-written corpora excluded above.",
)

#: Reported, never edited: data/demo is shared and a lane does not rewrite a peer's demo.
KNOWN_WRONG_LINES = (
    {
        "file": "data/demo/gastrulation.bio",
        "line": 3,
        "text": "Mutual repression makes the choice sharp.",
        "why_it_is_false": "the repression in the program is not mutual. SOX17 represses TBXT and "
        "SOX2 and is repressed by nothing: `d340e5c` measured SOX17's repression factor at exactly "
        "1.00000. The header describes a circuit the file does not contain.",
        "action": "reported, not edited: data/demo is a shared demo directory.",
    },
)

# --- reading a program --------------------------------------------------------------------------

_NUM = re.compile(r"(?<![\w.])[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")

#: Stripped from an evidence string before any number is extracted, so that a citation's year,
#: volume, page, figure label or accession digits cannot be read as a quoted value. A `hill: 3`
#: against `Fig. 3C` is exactly the false match this guards: without the strip it would be counted
#: sourced, and it is not.
_FURNITURE = (
    re.compile(r"\b(19|20)\d{2}[a-z]?\b"),  # a year
    re.compile(r"\b\d+\s*:\s*\d+(?:\s*[-–]\s*\d+)?"),  # volume:page, volume:page-page
    re.compile(r"\b(?:Fig(?:ure)?s?|Table|Panel|Suppl\w*|Extended Data)\.?\s*[A-Za-z]?\d+[A-Za-z]?", re.I),
    re.compile(r"\b(?:v|version)\s*\.?\s*\d+(?:\.\d+)*", re.I),  # GENCODE v50
    re.compile(r"\b[A-Z]{1,6}\d{3,}(?:\.\d+)?\b"),  # accessions: P48431, ENST00000329198.5
    re.compile(r"\bGENCODE\s*\d+\b", re.I),
    re.compile(r"\bBNID\s*\d+\b", re.I),
    re.compile(r"\bE-[A-Z]+-\d+\b"),
    re.compile(r"\bCL:\d+\b"),
    re.compile(r"\bdoi:[^\s;]+", re.I),
    re.compile(r"\bPMC\d+\b"),
    re.compile(r"\b\d+\s*e\b", re.I),  # "MBoC 6e"
)

#: A publication, a database or a reported range named in an `inferred` note. Matching any of these
#: separates A_asserted_names_a_work from A_asserted_bare.
_NAMES_A_WORK = (
    re.compile(r"\b[A-Z][A-Za-zÀ-ɏ'’-]+\s+(?:et al\.|&|and)\s"),  # Zhou et al. / X & Y
    re.compile(r"\b(?:GENCODE|UniProt|BioNumbers|BNID|Cell Ontology|Reactome|Ensembl|MBoC)\b", re.I),
    re.compile(r"\b(19|20)\d{2}\b"),
)

_EVIDENCE_KINDS_WITH_A_SOURCE = ("experimental", "curated", "predicted")

#: Hand adjudications promoting an E_existence_only slot to S_derived, by (file, line, field), each
#: with the arithmetic that justifies it. The registration fixes the TEST and the candidate set -
#: only a slot whose evidence string quotes some number can be promoted, and `classify_slot` marks
#: every such candidate - so this table can grow only inside a set the registration already exposed,
#: and every entry is a sentence a reader can check. It is empty in the registration commit and
#: filled in the result commit, where each entry is also printed in full.
DERIVATIONS: dict[tuple[str, int, str], str] = {}


def strip_furniture(s: str) -> str:
    """An evidence string with its citation digits removed, ready for number extraction."""
    for pat in _FURNITURE:
        s = pat.sub(" ", s)
    return s


def numbers_in(s: str) -> list[float]:
    out: list[float] = []
    for m in _NUM.finditer(s):
        try:
            out.append(float(m.group()))
        except ValueError:  # pragma: no cover - the regex only matches parseable literals
            continue
    return out


def leading_number(value: str) -> float | None:
    """The number a field's written value starts with, or None when it is not a number.

    A unit after the number (`0.5 h`, `94 fL`, `100 nM`, `20 min`) does not make the slot
    non-numeric: the quantity is still written down. `from Glycolysis, OXPHOS` is non-numeric.
    """
    m = _NUM.match(value.strip())
    if m is None:
        return None
    try:
        return float(m.group())
    except ValueError:  # pragma: no cover
        return None


def evidence_kind_and_text(props: dict[str, str]) -> tuple[str, str]:
    """A declaration's evidence kind and the rest of its evidence line, verbatim."""
    raw = props.get("evidence", "").strip()
    if not raw:
        return "", ""
    head, _, rest = raw.partition(" ")
    kind = head.strip().lower()
    if kind not in ("experimental", "curated", "predicted", "inferred", "none"):
        return "", raw
    return kind, rest.strip()


@dataclass
class Decl:
    """One declaration of a program, with the lines around it kept for context only."""

    file: str
    kind: str
    name: str
    line: int
    props: dict[str, str] = field(default_factory=dict)
    header: str = ""
    comments: list[str] = field(default_factory=list)


def read_program(path: str | Path) -> list[Decl]:
    """Every declaration of one program, nested blocks included, in source order."""
    p = Path(path)
    lines = p.read_text().splitlines()
    _, top = _parse_blocks(lines)
    flat: list[Decl] = []

    def walk(blocks: list[Any]) -> None:
        for b in blocks:
            flat.append(
                Decl(
                    file=p.as_posix(),
                    kind=b.kind,
                    name=b.header.strip(),
                    line=b.line,
                    props=dict(b.props),
                    header=b.header.strip(),
                )
            )
            walk(b.children)

    walk(top)
    flat.sort(key=lambda d: d.line)
    # the contiguous comment run immediately above a declaration, and the trailing comments inside
    # its span. Context only: a comment NEVER promotes a class, because a run of comment lines above
    # four rules belongs to all four and attaching it to the first would be an invented convention.
    starts = [d.line for d in flat] + [len(lines) + 1]
    for i, d in enumerate(flat):
        above: list[str] = []
        j = d.line - 2  # 0-based index of the line above the header
        while j >= 0 and lines[j].lstrip().startswith("#"):
            above.append(lines[j].strip())
            j -= 1
        inside = [
            ln[ln.index("#") :].strip()
            for ln in lines[d.line - 1 : starts[i + 1] - 1]
            if "#" in ln and not ln.lstrip().startswith("#")
        ]
        d.comments = list(reversed(above)) + inside
    return flat


def param_value(header: str) -> str:
    """The written value of a `param NAME = VALUE [unit]` header."""
    _, _, rest = header.partition("=")
    return rest.strip()


def written_value(d: Decl, slot: Slot) -> str | None:
    """What the program writes for this slot, or None when it writes nothing."""
    if slot.kind == "param" and slot.field == "value":
        v = param_value(d.header)
        return v or None
    return d.props.get(slot.field)


# --- classifying ---------------------------------------------------------------------------------


def repeated_values(decls: list[Decl], slot: Slot) -> dict[float, int]:
    """How many declarations of this kind in this file write each value of this field."""
    counts: dict[float, int] = {}
    for d in decls:
        if d.kind != slot.kind:
            continue
        raw = written_value(d, slot)
        if raw is None:
            continue
        n = leading_number(raw)
        if n is None:
            continue
        counts[n] = counts.get(n, 0) + 1
    return counts


def classify_slot(d: Decl, slot: Slot, repeats: dict[float, int]) -> dict[str, Any]:
    """The class of one (declaration, field) slot under CASCADE, with its evidence."""
    out: dict[str, Any] = {
        "file": d.file,
        "line": d.line,
        "kind": d.kind,
        "declaration": d.name,
        "field": slot.field,
        "tier": slot.tier,
    }
    raw = written_value(d, slot)
    if raw is None:
        out["written"] = False
        out["value"] = None
        if slot.language_default is NO_NUMBER:
            out["cls"] = "D_no_number"
            out["why"] = slot.note or "the language supplies no number for this field"
        else:
            out["cls"] = "D_language_default"
            out["value"] = slot.language_default
            out["why"] = f"omitted; the language supplies {slot.language_default}"
        return out

    out["written"] = True
    out["value_text"] = raw.strip()
    n = leading_number(raw)
    if n is None:
        out["cls"] = "N_non_numeric"
        out["why"] = "written, but the value is not a number"
        return out
    out["value"] = n
    kind, text = evidence_kind_and_text(d.props)
    out["evidence_kind"] = kind or "absent"
    out["evidence_text"] = text
    tags = []
    if slot.language_default is not NO_NUMBER and n == slot.language_default:
        tags.append("equals_language_default")
    if repeats.get(n, 0) >= REPEATED_IN_FILE_FLOOR:
        tags.append("repeated_in_file")
    out["convention_tags"] = tags
    out["repeated_in_file_count"] = repeats.get(n, 0)

    if kind not in _EVIDENCE_KINDS_WITH_A_SOURCE:
        names_a_work = any(p.search(text) for p in _NAMES_A_WORK)
        out["cls"] = "A_asserted_names_a_work" if names_a_work else "A_asserted_bare"
        out["why"] = f"evidence kind `{kind or 'absent'}`, which the grammar defines as not measured; " + (
            "the note names a work, which is context and not a source for the number"
            if names_a_work
            else "the note names no work"
        )
        return out

    stripped = strip_furniture(text)
    out["evidence_stripped"] = " ".join(stripped.split())
    quoted = numbers_in(stripped)
    out["numbers_quoted_after_strip"] = quoted
    if any(abs(q - n) <= 1e-12 for q in quoted):
        out["cls"] = "S_quoted"
        out["why"] = "the written value appears literally in the stripped evidence string"
        out["adjudication"] = "pending hand adjudication"
        return out
    out["cls"] = "E_existence_only"
    out["why"] = "a source is cited and quotes no number equal to the written value"
    # The only slots that can be promoted to S_derived by hand are those whose evidence string
    # quotes some number: the candidate set is fixed mechanically here, by the registration's own
    # matcher, so the promotion cannot reach a slot the registration did not expose.
    out["derivation_candidate"] = bool(quoted)
    key = (d.file, d.line, slot.field)
    if key in DERIVATIONS:
        out["cls"] = "S_derived"
        out["arithmetic"] = DERIVATIONS[key]
        out["why"] = "an exact transform of a number the same evidence string quotes, stated in it"
    return out


def census_file(path: str | Path) -> list[dict[str, Any]]:
    """Every scoped slot of one program, classified."""
    decls = read_program(path)
    repeats = {(s.kind, s.field): repeated_values(decls, s) for s in SCOPE}
    rows: list[dict[str, Any]] = []
    for d in decls:
        for slot in SCOPE:
            if slot.kind != d.kind:
                continue
            rows.append(classify_slot(d, slot, repeats[(slot.kind, slot.field)]))
    return rows


def programs(root: Path | None = None) -> list[Path]:
    """The population, as files, in a fixed order."""
    base = root or Path()
    out: list[Path] = []
    for d in POPULATION_DIRS:
        out.extend(sorted((base / d).glob("*.bio")))
    return out


def registration() -> dict[str, Any]:
    """Everything decided before the count, as the registration writes it."""
    return {
        "question": "which numbers in the hand-authored BioLang corpus are traceable to a source "
        "that states a number, which inherit a default convention, and which are asserted - so that "
        "the model's qualitative/quantitative boundary is counted rather than asserted",
        "why_now": "`d340e5c` and `53d9936` established by measurement that the gastrulation module "
        "cannot discriminate mechanisms and that its number convention, not its biology, is the "
        "cause: three different sourced repressors of SOX17 restore a band and agree to the digit, "
        "because the module's own `threshold 2.0; hill 3` makes any added inhibitor a switch rather "
        "than a brake.",
        "population": {
            "directories": list(POPULATION_DIRS),
            "rule": "every .bio program a person wrote in these two directories. Counts stay "
            "separable per file and the two directories are never pooled in a headline.",
            "excluded": dict(EXCLUDED_CORPORA),
        },
        "scoped_fields": [
            {
                "kind": s.kind,
                "field": s.field,
                "tier": s.tier,
                "language_default": "no number" if s.language_default is NO_NUMBER else s.language_default,
                "note": s.note,
            }
            for s in SCOPE
        ],
        "tiers": {
            "1": "the dynamics core: every numeric field of a rule, a gene or a protein, and every "
            "param value. These are the numbers genomeos/runtime/grn.py integrates.",
            "2": "the rest of the declared numbers: timers, stages, compartments, pools, fields, "
            "signals, events, regimes, decisions. Own denominator, never pooled into tier 1's.",
        },
        "out_of_scope_fields": dict(OUT_OF_SCOPE_FIELDS),
        "classes": dict(CLASSES),
        "cascade": list(CASCADE),
        "convention_tags": dict(CONVENTION_TAGS),
        "repeated_in_file_floor": REPEATED_IN_FILE_FLOOR,
        "citation_furniture_stripped": [p.pattern for p in _FURNITURE],
        "why_furniture_is_stripped": "without it a `hill: 3` whose evidence string reads `Fig. 3C` "
        "would be counted as a quoted number and called sourced. It is not.",
        "comments_are_context_only": "a comment above or inside a declaration is recorded and never "
        "promotes a class. A comment run above four rules belongs to all four, and attaching it to "
        "the first would be a convention invented while counting.",
        "s_derived_procedure": "S_derived is the one class no matcher can decide, so the "
        "registration fixes its test and its candidate set instead of its members: a slot may be "
        "promoted to S_derived only if it classified E_existence_only AND its evidence string "
        "quotes at least one number after the strip, which `classify_slot` marks as "
        "`derivation_candidate`. Each promotion is one entry in `DERIVATIONS` carrying the "
        "arithmetic, committed with the result and printed in full. The table is empty here.",
        "s_derived_candidates_at_registration": "unknown at registration by construction; the count "
        "of candidates is reported with the result.",
        "reports": list(REPORTS),
        "denominator": "written slots + D_language_default slots, stated per file and per tier. "
        "D_no_number and N_non_numeric are counted and reported outside it.",
        "exclusions": list(EXCLUSIONS),
        "reported_not_edited": list(KNOWN_WRONG_LINES),
        "whichever_way_it_falls": "a census finding almost nothing sourced is the finished result, "
        "not a failure: it draws the qualitative boundary in counts and names the rules a "
        "measurement would have to supply before any mechanism claim is possible.",
    }
