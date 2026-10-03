# SPDX-License-Identifier: AGPL-3.0-or-later
"""Chromatin evidence for the cell a compiled rule names, or a label saying why there is none.

R1 of the external review (2026-09-28) made every compiled rule carry its cell as an executable
`when: cell_type = ...`. That cell is **asserted**: it is the tissue of the winning AlphaGenome
track in `predict_target`, or the cell a CRISPRi screen silenced in, and until now nothing asked
whether the element is open in that cell at all. Reader v1 already holds the chromatin: ENCODE
DNase-seq narrowPeak sets for thirteen biosamples, genome-wide, cached under `data/results`.

This module reads one state per rule and the compiler writes it beside the rule, so an asserted
context becomes either an evidenced one or a labelled one.

**The states are named by the measurement, never by a verdict.**

    open_in_reader       the element's locus overlaps a DNase narrowPeak of the reader biosample
                         this rule's cell maps to
    not_open_in_reader   no narrowPeak of that biosample overlaps the locus
    not_assessable       the reading could not be taken, with the reason

**A fourth value, which is not a state, was adopted on 2026-10-03** (Albert, item (7)): `not_assessed`,
the absence of a reading. The three above are readings and this one says none was ATTEMPTED - the
mapping table was not on the machine that compiled the program, so nothing was tried on any rule. Since
that date the compiler WRITES it onto the rule line in that case instead of omitting the field, so
"nobody looked" is readable from the rule itself rather than inferred from a silence, and `parse_value`
accepts it. It is kept out of `STATES`, which stays the three readings, and out of `value()`, which
writes a reading: it is never `not_assessable`, a reading that was attempted and could not be taken and
which carries its reason, so `not_assessed` is accepted bare and refused with one.

`not_open_in_reader` is **not** "closed", and it is not "contradicted". Peak absence is *not
detected open at the reader's registered call*: a peak set is a call set, so a locus with no peak
may be shut, may be open below the caller's sensitivity in that experiment, or may be open in a
state of the cell the experiment did not sample. Nothing here converts absence into closure.

**`not_assessable` is split in two, because a never-looked state and a measured one are different
outputs.**

    no_reader_for_this_cell        reader v1 never read this cell: the rule's cell carries no
                                   registered ontology term, or its term is not one of the thirteen
                                   biosamples, or its label resolves to more than one term
    region_outside_measured_span   the cell does have a reader, but this chromosome's peak set for
                                   it covers nothing at this locus: no peak file, no peak in it, or
                                   a locus beyond the outermost peak of the file

**The cell is mapped by ontology term, never by comparing names.** A rule's cell label is the
`ident`-ed name of an AlphaGenome track biosample (or of the GTEx tissue a track names), and
AlphaGenome's own track metadata, already on disk at `TRACK_METADATA`, gives each of those names its
`ontology_curie` (EFO, CL or UBERON). Each reader biosample carries the same kind of term in
`READER_TERMS`, and a rule is assessable exactly when its label resolves to one term and that term
is a reader biosample's. A label the metadata gives two or more terms is ambiguous and is
`not_assessable`; so is a label the metadata does not carry at all.

**No new cut-off is introduced.** The openness call is reader v1's own: `reader.EVIDENCE`, the
released ENCODE DNase-seq narrowPeak set for the biosample, as the reader loads it. The reader
applies that call to a promoter window; this module applies the same call to the element's own
interval. The only other quantity here is the peak file's covered span, which is a statement about
where reads were placed and not a threshold on signal.

**The weakness, stated here and in every result this module writes.** Reader v1's openness comes
from ENCODE DNase and histone peaks, and AlphaGenome was trained on ENCODE. `open_in_reader` is
therefore a **consistency check between two readings of the same chromatin** and never independent
evidence that the rule is right. It must never be quoted as validation. The informative direction is
the other one: a `not_open_in_reader` rule is a *located candidate defect*, a place where the
compiled program asserts activity in a cell whose own chromatin, read directly, shows no open
element there.
"""

from __future__ import annotations

import csv
import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from genomeos import manifest as _mf
from genomeos.genome import reader
from genomeos.provenance import rule_evidence_tier as ret

# ---- the states and the reasons, named by the measurement -------------------------------------

STATE_OPEN = "open_in_reader"
STATE_NOT_OPEN = "not_open_in_reader"
STATE_NOT_ASSESSABLE = "not_assessable"
STATES = (STATE_OPEN, STATE_NOT_OPEN, STATE_NOT_ASSESSABLE)

#: The FOURTH value a rule's `context_evidence` may carry, adopted 2026-10-03 on Albert's item (7)
#: ("adopt the header fix, not_assessed, visible defaults and per-slot evidence tier; defer cell
#: provenance"). It is the only one of the four that is NOT a reading: `not_assessed` says nothing was
#: tried - the mapping table was not on the machine that compiled the program, so no reading was
#: attempted on any rule. It is NOT `not_assessable`, which is a reading that WAS attempted and could
#: not be taken, and which carries its reason; so `not_assessed` may never carry one, and
#: `genomeos.provenance.compiled_defaults.check_not_assessed_carries_no_reason` and `parse_value`
#: below both refuse a value that gives it one.
#:
#: The spelling is IMPORTED and not re-typed. `rule_evidence_tier.NOT_ASSESSED` is the one place this
#: project spells the absence of an assessment, and a second spelling is how two levels of the
#: language start disagreeing about the same fact.
#:
#: It is deliberately NOT added to `STATES`, and that is the whole shape of this adoption. `STATES` is
#: the set of READINGS: `value()` writes a reading and still refuses this one, `registration()`
#: registers the readings and is byte-unmoved, and `tests/test_context_evidence.py` pins `STATES` to
#: its three members. Only `parse_value`, which reads what a program may legitimately CARRY, knows the
#: fourth - because a program may now carry it, and a reader of a program must not choke on it.
STATE_NOT_ASSESSED = ret.NOT_ASSESSED

#: Every value `parse_value` accepts: the three readings and the absence of one. `STATES` is unchanged.
#:
#: It is a DECLARATION and not the membership test. `parse_value` takes the fourth value in a clause
#: that returns before the committed three-state test it leaves intact, so this tuple could in
#: principle drift away from what that function really accepts. It cannot drift silently:
#: `tests/test_context_evidence.py` feeds every member of this tuple through `parse_value` and then
#: shows that a value outside it is refused, so the declaration is held to the behaviour.
VALUES = (*STATES, STATE_NOT_ASSESSED)

REASON_NO_READER = "no_reader_for_this_cell"
REASON_OUTSIDE_SPAN = "region_outside_measured_span"
NOT_ASSESSABLE_REASONS = (REASON_NO_READER, REASON_OUTSIDE_SPAN)

#: Words a state name may never contain. A state is a measurement; "supported" and "contradicted"
#: are verdicts, and this module's whole point is that it issues none. tests/test_context_evidence.py
#: fails if any name above picks one of these up.
FORBIDDEN_IN_A_STATE_NAME = (
    "support",
    "contradict",
    "validat",
    "confirm",
    "refut",
    "verif",
    "consistent",
    "agree",
    "disagree",
    "correct",
    "wrong",
    "true",
    "false",
    "good",
    "bad",
    "pass",
    "fail",
)

READING = {
    STATE_OPEN: (
        "the element's locus overlaps a DNase narrowPeak of the reader biosample this rule's cell "
        "maps to by ontology term"
    ),
    STATE_NOT_OPEN: (
        "no DNase narrowPeak of that biosample overlaps the locus: not detected open at the "
        "reader's registered call, which is not the same as closed"
    ),
    STATE_NOT_ASSESSABLE: "no reading was taken; the reason is carried with the state",
    REASON_NO_READER: (
        "reader v1 never read this cell: its label carries no single registered ontology term, or "
        "that term is not one of the thirteen biosamples"
    ),
    REASON_OUTSIDE_SPAN: (
        "the cell has a reader, but that biosample's peak set for this chromosome covers nothing at "
        "this locus: no file, no peak in it, or a locus beyond its outermost peak"
    ),
}

NOT_VALIDATION = (
    "a consistency check between two readings of the same ENCODE chromatin, not validation: reader "
    "v1's openness comes from ENCODE peaks and AlphaGenome was trained on ENCODE, so open_in_reader "
    "is never independent evidence that the rule is right and must never be quoted as validation"
)

NOT_CLOSED = (
    "not_open_in_reader means not detected open at the reader's registered threshold, not closed: a "
    "narrowPeak set is a call set, and absence is absence of a call"
)

OPENNESS_CALL = (
    "reader v1's own call, unchanged: " + reader.EVIDENCE + ", as genomeos.genome.reader loads it. "
    "The reader tests a promoter window for an overlapping peak; this tests the element's own "
    "interval for one, which is the element-level call genomeos.attribution.closure has used since "
    "it was written ('an element is active in the cell when it overlaps a DNase peak there'). No new "
    "cut-off is introduced and no signal value is thresholded."
)

SPAN_CALL = (
    "a biosample's measured span on a chromosome is the interval from the first peak's start to the "
    "last peak's end in its cached narrowPeak file for that chromosome. Beyond it no read was "
    "placed, so an absent peak there is not a negative call. This is a statement about coverage, "
    "not a threshold on signal."
)

# ---- the mapping table ------------------------------------------------------------------------

#: AlphaGenome's own track metadata, already on disk: `biosample_name`, `gtex_tissue` and the
#: `ontology_curie` each carries. This is the authority for the rule side of the mapping, because a
#: rule's cell label is the `ident`-ed name of one of these tracks.
TRACK_METADATA = Path("data/cache/entex/alphagenome_track_metadata_copy.csv")

#: The thirteen biosamples reader v1 read genome-wide, by the ENCODE `biosample_ontology.term_name`
#: `genomeos.genome.reader` queried, with the ontology term each one is. Every value is checked
#: against TRACK_METADATA by a test, so this hand-written table cannot drift from the authority.
READER_TERMS = {
    "K562": "EFO:0002067",
    "HepG2": "EFO:0001187",
    "GM12878": "EFO:0002784",
    "IMR-90": "EFO:0001196",
    "H1": "EFO:0003042",
    "SK-N-SH": "EFO:0003072",
    "astrocyte": "CL:0000127",
    "hepatocyte": "CL:0000182",
    "keratinocyte": "CL:0000312",
    "cardiac muscle cell": "CL:0000746",
    "CD14-positive monocyte": "CL:0001054",
    "ovary": "UBERON:0000992",
    "testis": "UBERON:0000473",
}

#: Any well-formed CURIE is read from the metadata, so a label keeps the term it actually carries and
#: lands in `not_assessable` for the right reason. Only the three ontologies a reader biosample can be
#: named in may ever match one; the metadata also uses NTR (ENCODE new-term requests) and CLO (Cell
#: Line Ontology), and a label carrying one of those can never be a reader biosample.
TERM_PATTERN = re.compile(r"^[A-Z]{2,10}:\d{4,9}$")
READER_PREFIXES = ("EFO", "CL", "UBERON")


def term_to_reader() -> dict[str, str]:
    """Ontology term -> reader biosample. One term may name only one biosample, or the table is wrong."""
    out: dict[str, str] = {}
    for cell, term in READER_TERMS.items():
        if term in out:
            raise ValueError(f"two reader biosamples claim {term}: {out[term]!r} and {cell!r}")
        out[term] = cell
    return out


def ident(s: str) -> str:
    """The identifier form the compiler writes a cell name in (`genomeos.attribution.compile.ident`)."""
    out = re.sub(r"[^A-Za-z0-9_]", "_", s)
    return out if not out[:1].isdigit() else f"g_{out}"


@dataclass(frozen=True)
class Mapping:
    """Cell label -> ontology term, from the track metadata, with what could not be mapped.

    `terms` holds the labels that resolve to exactly one term. `ambiguous` holds the labels the
    metadata gives two or more terms, which are `not_assessable` by the registered rule and are
    reported rather than silently resolved.
    """

    terms: dict[str, str]
    ambiguous: dict[str, tuple[str, ...]]
    rows: int
    source: str
    prefixes: dict[str, int] = field(default_factory=dict)

    def term_of(self, label: str) -> str | None:
        return self.terms.get(label)


_MAPPING: Mapping | None = None


def mapping(path: Path = TRACK_METADATA) -> Mapping:
    """The registered label-to-term table, read once from the track metadata on disk."""
    global _MAPPING
    if _MAPPING is not None and _MAPPING.source == str(path):
        return _MAPPING
    if not path.exists():
        raise FileNotFoundError(
            f"the AlphaGenome track metadata is not on this machine: {path}. The cell-to-biosample "
            "mapping is read from it and is never guessed from names."
        )
    seen: dict[str, set[str]] = {}
    rows = 0
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh):
            rows += 1
            term = (row.get("ontology_curie") or "").strip()
            if not TERM_PATTERN.match(term):
                continue
            for name in (row.get("biosample_name"), row.get("gtex_tissue")):
                name = (name or "").strip()
                if name:
                    seen.setdefault(ident(name), set()).add(term)
    terms = {k: next(iter(v)) for k, v in seen.items() if len(v) == 1}
    ambiguous = {k: tuple(sorted(v)) for k, v in seen.items() if len(v) > 1}
    prefixes: dict[str, int] = {}
    for term in terms.values():
        prefix = term.split(":", 1)[0]
        prefixes[prefix] = prefixes.get(prefix, 0) + 1
    _MAPPING = Mapping(terms=terms, ambiguous=ambiguous, rows=rows, source=str(path), prefixes=prefixes)
    return _MAPPING


def reader_for(label: str, table: Mapping | None = None) -> tuple[str | None, str | None, str]:
    """(reader biosample, ontology term, why not) for a rule's cell label.

    The biosample is None exactly when no reading can be taken for the cell, and the third value
    then says which of the registered reasons applies.
    """
    table = table or mapping()
    if label in table.ambiguous:
        return (
            None,
            None,
            "the label resolves to more than one ontology term: " + ", ".join(table.ambiguous[label]),
        )
    term = table.term_of(label)
    if term is None:
        return None, None, "no registered ontology term for this cell label"
    cell = term_to_reader().get(term)
    if cell is None:
        return None, term, f"{term} is not one of the thirteen biosamples reader v1 read"
    return cell, term, ""


# ---- the chromatin ----------------------------------------------------------------------------


@dataclass
class Readers:
    """Reader v1's DNase peak sets, loaded per biosample and chromosome on first use."""

    results_dir: Path = reader.RESULTS
    _index: dict[tuple[str, str], reader.PeakIndex | None] = field(default_factory=dict)
    _span: dict[tuple[str, str], tuple[int, int] | None] = field(default_factory=dict)

    def _load(self, biosample: str, chrom: str) -> None:
        key = (biosample, chrom)
        path = self.results_dir / reader.peaks_path(biosample, chrom).name
        peaks: list[tuple[int, int, float]] = []
        if path.exists():
            with gzip.open(path, "rt") as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    s, e, v = line.rstrip("\n").split("\t")
                    peaks.append((int(s), int(e), float(v)))
        if not peaks:
            self._index[key] = None
            self._span[key] = None
            return
        self._index[key] = reader.PeakIndex(peaks)
        self._span[key] = (min(s for s, _, _ in peaks), max(e for _, e, _ in peaks))

    def index(self, biosample: str, chrom: str) -> reader.PeakIndex | None:
        if (biosample, chrom) not in self._index:
            self._load(biosample, chrom)
        return self._index[(biosample, chrom)]

    def span(self, biosample: str, chrom: str) -> tuple[int, int] | None:
        if (biosample, chrom) not in self._span:
            self._load(biosample, chrom)
        return self._span[(biosample, chrom)]

    def cells_loaded(self) -> list[str]:
        return sorted({b for b, _ in self._index})


# ---- the state --------------------------------------------------------------------------------


def value(state: str, *parts: str) -> str:
    """The property value the compiler writes: the state first, then what qualifies it."""
    if state not in STATES:
        raise ValueError(f"{state!r} is not one of {STATES}")
    return ", ".join([state, *parts])


# TAUGHT THE FOURTH VALUE on 2026-10-03, Albert's item (7), and taught it by ADDITION: the docstring
# below, the split, the strip, the tuple, the three-state membership test and the `ValueError` it
# raises are the lines `bd0e56a` committed, unchanged and in the same order. `not_assessed` is handled
# by a clause ABOVE them and returns before they are reached, so nothing about the three readings can
# have moved - there is no line left for it to have moved in. (It is also why the ValueError below
# still names three states and not four: that line is the committed one. The fourth value never
# reaches it.)
#
# Two reasons for the shape rather than one. The first is the project's: this lane may not remove a
# committed line, `--force` is denied to a lane, and `dbb5d4a` established that the route left open is
# to add. The second is the better one: "the three readings behave exactly as before" is a claim, and
# a clause that returns before their code is reached turns it into a fact a reader can check by eye.
#
# `not_assessed` is accepted BARE only. With anything after the comma it is refused, for the same
# reason `compiled_defaults.check_not_assessed_carries_no_reason` refuses it in a whole program: a
# reason is the mark of an attempt, and nothing was attempted. The value that was attempted and could
# not be taken is `not_assessable`, and `not_assessable, <reason>, <cell>` is accepted as it always was.
def parse_value(text: str) -> tuple[str, tuple[str, ...]]:
    """(state, the rest) from a written value. The inverse of `value`."""
    bits = [b.strip() for b in text.split(",")]
    if bits and bits[0] == STATE_NOT_ASSESSED:
        if bits[1:]:
            raise ValueError(
                f"{text!r} gives {STATE_NOT_ASSESSED!r} a reason. A reason is the mark of an attempt "
                f"and nothing was attempted; the value that was attempted and not takeable is "
                f"{STATE_NOT_ASSESSABLE!r}, and that is the one which carries its reason."
            )
        return bits[0], ()
    if not bits or bits[0] not in STATES:
        raise ValueError(f"{text!r} does not begin with one of {STATES}")
    return bits[0], tuple(bits[1:])


def state_for(
    label: str, chrom: str, start: int, end: int, readers: Readers, table: Mapping | None = None
) -> str:
    """The context-evidence value for one rule: its cell label and its element's locus.

    `label` is the rule's `when: cell_type` value exactly as the compiler writes it, so `unknown`
    (R1's never-recorded context) arrives here as a label with no ontology term and is reported as
    `no_reader_for_this_cell`, which is what it is.
    """
    biosample, _term, _why = reader_for(label, table)
    if biosample is None:
        return value(STATE_NOT_ASSESSABLE, REASON_NO_READER)
    span = readers.span(biosample, chrom)
    if span is None or not (end > span[0] and start < span[1]):
        return value(STATE_NOT_ASSESSABLE, REASON_OUTSIDE_SPAN, biosample)
    idx = readers.index(biosample, chrom)
    assert idx is not None  # a span exists only where peaks were loaded
    if idx.overlapping(start, end):
        return value(STATE_OPEN, biosample)
    return value(STATE_NOT_OPEN, biosample)


# ---- the rules, enumerated the way the compiler emits them ------------------------------------


def rule_loci(chrom: str, results_dir: Path | None = None, layer: Any = None) -> list[tuple[str, int, int]]:
    """(cell label, start, end) for every rule `compile_chromosome` emits for this chromosome.

    The compiled text carries each rule's cell but not its locus, so a reading that needs both - the
    base-rate comparison, for one - enumerates the rules here instead. It mirrors the compiler block
    for block, and `tests/test_context_evidence.py` holds it to the compiler's own output on chr21,
    so the two cannot drift apart.
    """
    from genomeos.attribution import compile as cp
    from genomeos.attribution.measured import rule_links

    results_dir = results_dir if results_dir is not None else cp.RESULTS_DIR
    out: list[tuple[str, int, int]] = []
    elements = cp._attributed(chrom, results_dir)
    for e in elements:
        pc = e["predicted_coding"]
        out.append((cp.context(pc.get("tissue")), e["start"], e["end"]))
    _layer, measured_rows = cp._measured_rows(chrom, elements, results_dir, layer)
    for row in measured_rows:
        for _gene, _action, _strength, cell, _split in rule_links(row):
            out.append((cp.context(cell), row["start"], row["end"]))
    return out


# ---- the stamp: which uncommitted code could have reached a count -----------------------------

REPO_ROOT = Path(__file__).resolve().parents[2]

#: Re-exported, not redefined. This module held one of the ten copies of `code_cleanliness` and one of
#: the eight copies of the closure behind it; the copies had drifted into three algorithms that gave
#: three different answers about the same tree, so there is now one implementation in
#: `genomeos/manifest.py` (tests/test_code_cleanliness_shared.py fails if a second one appears). The
#: names stay bound here because five writers reach them through this module
#: (scripts/context_evidence_baserate.py, scripts/context_evidence_register.py,
#: scripts/not_open_profile.py, scripts/repress_population.py, scripts/response_map_increment2.py) and
#: the shared signatures are the ones this module already used: `(entry, own_code)` and `(entry)`.
code_cleanliness = _mf.code_cleanliness
counting_path = _mf.counting_path


def registration() -> dict[str, Any]:
    """Everything a result or a pre-registration states about this reading, from the code itself."""
    table = mapping()
    return {
        "states": {s: READING[s] for s in STATES},
        "not_assessable_reasons": {r: READING[r] for r in NOT_ASSESSABLE_REASONS},
        "state_names_carry_no_verdict": {
            "rule": "no state or reason name may contain any of these stems",
            "forbidden_stems": list(FORBIDDEN_IN_A_STATE_NAME),
            "tested_by": "tests/test_context_evidence.py",
        },
        "not_closed": NOT_CLOSED,
        "not_validation": NOT_VALIDATION,
        "openness_call": OPENNESS_CALL,
        "measured_span": SPAN_CALL,
        "mapping": {
            "by": "ontology term (EFO, CL or UBERON), never by comparing cell names",
            "rule_side_authority": str(TRACK_METADATA),
            "rule_side_columns": ["biosample_name", "gtex_tissue", "ontology_curie"],
            "metadata_rows": table.rows,
            "labels_with_one_term": len(table.terms),
            "labels_per_ontology": dict(sorted(table.prefixes.items())),
            "ontologies_a_reader_biosample_can_be_named_in": list(READER_PREFIXES),
            "labels_with_more_than_one_term": len(table.ambiguous),
            "ambiguous_labels": {k: list(v) for k, v in sorted(table.ambiguous.items())},
            "reader_biosamples": dict(READER_TERMS),
            "ambiguous_is": STATE_NOT_ASSESSABLE + ", " + REASON_NO_READER,
        },
        "denominator": "every compiled rule, counted as the compiler emits it",
        "alphagenome_requests": 0,
    }
