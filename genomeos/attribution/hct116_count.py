# SPDX-License-Identifier: AGPL-3.0-or-later
"""A blinded count of the HCT116 links the partial `elements_hct116` cache can answer, and the pooled
arms against the two imported floors. No direction, no value and no sign is read anywhere.

**What came before.** lane-direction registered a link-level direction test
(`data/results/direction_link_registration.json`) and its gate refused the comparison
(`data/results/direction_link.json`): of the K562 links on an attributed element, the decreases arm met
both imported floors with 174 answered links over 112 independent loci, and the increases arm did not,
with 21 answered links against a floor of 30 and 13 independent loci against a floor of 20. That lane
read the K562 track only, substituted nothing, and named a second, partial cache in its registration as
present and deliberately not read. Its own words, carried verbatim here as `NAMED_NOT_READ`, are the
reason this lane exists and the reason this lane counts instead of measuring.

**What this lane is.** The only reason to look at a second cache is that one arm missed a floor. A lane
that read predicted directions while counting whether it has enough of them would be choosing its
population after seeing the outcome, which is the thing the floors exist to prevent. So this module
reads **answer presence** and nothing else: whether the cached window of a link's own element holds an
entry for the link's own gene on the link's own cell track. Whether that entry is negative, positive or
exactly zero is not read, is not returned, is not recorded and is not recoverable from anything this
lane computes.

**How that is made structurally impossible rather than merely avoided.** Exactly one function in this
module is ever handed a parsed cache record, `genes_with_an_answer`, and its return type is a
`frozenset` of gene **names**: the value has no path out of it. Inside it the only operation applied to
the value is an identity test against `None`. `tests/test_hct116_count.py` pins that three ways. First a
poison object, which raises `BlindingError` from every comparison, every arithmetic operation,
every conversion to a number, string or bool, every attribute access and every hash, is put where a
cached value would be: if the reader did anything to the value other than test it against `None`, the
test fails with that violation instead of passing. Second, the same records are read with a negative
value, a positive value and an exactly-zero value, and the reader's output is asserted **identical** in
all three, so no sign is distinguishable in this lane's output even in principle. Third, an AST check
over this module's own source asserts that `genes_with_an_answer` contains no ordering or equality
comparison, no arithmetic, and no call to `float`, `int`, `abs`, `round`, `str`, `repr` or `sorted`, so
a later edit cannot reintroduce a value read without failing the test.

**Presence is an upper bound, and that is stated before the count.** lane-direction's *answered* links
are those whose cached value is present **and not exactly zero** (`direction_link.SIGN_CALL`): an
exactly-zero value carries no predicted direction and is excluded from every denominator. Testing a
value against zero would be reading the value, so this lane does not do it. Every count here is
therefore **answer presence**, an upper bound on answered links, and the pooled counts this lane tests
against the floors are upper bounds. A count short of a floor at its upper bound is short of it, so an
upper bound that fails a floor settles the gate; an upper bound that cleared one would not, and this
module says so in `PRESENCE_IS_AN_UPPER_BOUND` rather than leaving it to a reader. lane-direction
reported `predicted_zero_excluded` 0 on all 230 of its in-cell links, so on the K562 side presence and
answered coincide there, and this lane reconciles its own K562 presence count against that lane's
answered count as a check.

**The locus convention.** cell2's grouping, imported with its own wording that it is an operational
grouping and **not** established biological independence. The gate itself is lane-direction's own
`gate`, imported and not re-implemented, so both floors come from the files that fixed them before
either lane existed: `fresh.POSITIVE_FLOOR` = 30 links and `fresh.LOCUS_FLOOR`, which is
`cell2.POOLED_LOCUS_FLOOR` = 20 independent loci. cell2's `group` does not look at the cell, so pooling
a second cell's links into an arm is **not** additive in loci: a link of either cell that shares a gene
with another, or lies within 1 Mb of it on the same chromosome, is the same locus. The pooled loci are
therefore grouped over the pooled links together and never added up per cell.

**The wording that is binding (review item R2).** A measured increase on knockdown is *an increase on
knockdown*. It is never a silencer, never a repressor and never evidence of a repression mechanism. A
count short of a floor is a no-go and is never close, promising, nearly enough, a good start or enough
for a pilot.

**HCT116 is not a fresh cell for this purpose.** Its decreases were scored in a39073d. Any later
registration that uses these links must say so, and `PRIOR_EXPOSURE` is the sentence it must carry.

**What this lane does not do.** It does not register a direction test, does not run one, does not read
or re-read a direction, does not reopen lane-direction's no-go, does not move a floor, does not widen
lane-direction's answered set and does not pool the two arms to reach a floor. It makes no model
request, downloads nothing and spends nothing.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from genomeos.attribution import cell2, fresh
from genomeos.attribution import direction_link as dl

# ---- what this lane inherits, carried word for word with its file --------------------------------

#: lane-direction's own sentence naming this cache as present and not read, verbatim. It is the whole
#: authority for opening it here, and it is quoted rather than paraphrased.
NAMED_NOT_READ = dl.CELL_AVAILABILITY

#: lane-direction's registered gate no-go, verbatim. This lane does not reopen it and does not weaken
#: it: it stands as that lane's reading of that lane's population.
INHERITED_GATE_NO_GO = dl.GATE_NO_GO

INHERITED_FROM = {
    "registration": "data/results/direction_link_registration.json",
    "result": "data/results/direction_link.json",
    "module": "genomeos/attribution/direction_link.py",
    "commits": {
        "code": "b3705c5",
        "registration": "0779c97",
        "no_go_wording": "d65e803",
        "result": "de0c879",
        "attribution_section": "232aa19",
    },
}

#: Every figure of lane-direction's this lane pools onto, with the file it is read from. None of them is
#: recomputed as a new finding here and none may be quoted under a noun its own breakdown does not
#: support. The script reads these same numbers out of the result file and reconciles against them, so a
#: retyped figure cannot drift from the record.
INHERITED_FIGURES = {
    "file": "data/results/direction_link.json",
    "eligible_links": 260,
    "eligible_by_arm": {"decreases": 212, "increases": 48},
    "eligible_by_cell": {"GM12878": 4, "HCT116": 7, "Jurkat": 5, "K562": 230, "WTC11": 14},
    "eligible_by_arm_and_cell": {
        "decreases": {"GM12878": 4, "HCT116": 4, "Jurkat": 5, "K562": 191, "WTC11": 8},
        "increases": {"HCT116": 3, "K562": 39, "WTC11": 6},
    },
    "cell_read": "K562",
    "answered": {"decreases": 174, "increases": 21},
    "predicted_zero_excluded": {"decreases": 0, "increases": 0},
    "absent": {"decreases": 17, "increases": 18},
    "independent_loci_answered": {"decreases": 112, "increases": 13},
    "increases_short_by": [
        "answered links 21, short of 30 by 9",
        "independent loci 13, short of 20 by 7",
    ],
}

# ---- the cells, the cache and what it is permitted to give up -------------------------------------

CELL = "HCT116"
PRIMARY_CELL = dl.PRIMARY_CELL
CACHE_ROOT = "data/knowledge/alphagenome/elements_hct116"
CACHE = f"{CACHE_ROOT}/<chrom>/<element_id>.json"
CACHE_SHAPE = (
    f"the partial cache is one JSON file per element, {CACHE}, and not the per-chromosome gzip archive "
    f"{dl.CACHE} the finished sweep wrote. A link's own element is therefore one small file, opened by "
    "name, and every open this lane makes is recorded in the result by path. Where the main archives "
    "are read for the K562 side they are STREAMED and never loaded: the reader is lane-direction's own "
    "`cached_records`, imported from scripts/direction_link.py rather than copied, which decompresses "
    "in chunks and materialises one record at a time, so peak memory is one record and not one "
    "chromosome"
)

PRESENCE_CALL = (
    "the only thing read out of a cached record is ANSWER PRESENCE: whether the element's cached window "
    "holds an entry for the link's own gene on the link's own cell track, under "
    f"`genes[*].by_cell[{CELL!r}]` for the gene whose name is exactly the link's gene, by the same path "
    "lane-direction's `PREDICTED_CALL` names. Presence is `the entry exists and is not null`. The entry "
    "itself is never returned, never stored, never printed and never compared with anything but None, "
    "and the result of this lane is identical whether that entry is negative, positive or exactly zero"
)

BLINDING = (
    "no predicted direction, no predicted value and no sign is read in this lane, and that is "
    "structural rather than careful. Exactly one function is handed a parsed cache record, "
    "`genes_with_an_answer`, and it returns a frozenset of gene NAMES, so the value has no path out of "
    "it; inside it the only operation applied to the value is an identity test against None. "
    "tests/test_hct116_count.py pins this three ways: a poison object that raises BlindingError "
    "from every comparison, arithmetic operation, numeric or string conversion, attribute access and "
    "hash is read in place of a cached value, so any use of the value other than an identity test "
    "fails the test; the same record is read with a negative, a positive and an exactly-zero value and "
    "the output is asserted identical in all three, so no sign is distinguishable in this lane's "
    "output even in principle; and an AST check over this module's own source asserts that the "
    "function contains no ordering or equality comparison, no arithmetic and no call to float, int, "
    "abs, round, str, repr or sorted, so a later edit cannot reintroduce a value read silently"
)

PRESENCE_IS_AN_UPPER_BOUND = (
    "every count here is ANSWER PRESENCE and is an UPPER BOUND on lane-direction's `answered`, because "
    "that lane's answered links are those whose cached value is present AND not exactly zero "
    "(direction_link.SIGN_CALL), and testing a value against zero would be reading the value. A count "
    "short of a floor at its upper bound is short of that floor, so an upper bound settles a gate it "
    "fails; an upper bound that CLEARED a floor would not settle the gate, and a lane reading it would "
    "have to establish the zero exclusion first. lane-direction reported predicted_zero_excluded 0 on "
    "all 230 of its in-cell links, so presence and answered coincide on its K562 side, and this lane "
    "reconciles its own K562 presence count against that lane's answered count as a check"
)

#: The reasons a link has no answer, kept apart because they are different facts about the cache.
NO_FILE = "no cached element file for this element in this partial cache"
NO_ENTRY = "the element's cached window holds no entry for this gene on this cell's track"
PRESENT = "an entry is present for this gene on this cell's track"

ANSWER_PRESENCE_GAP = (
    "the gap between eligible links and links with an answer is reported as its own number and split "
    "into its two reasons, because they are different facts about the cache: an element with NO FILE in "
    f"this partial cache ({NO_FILE!r}) and an element whose cached window holds NO ENTRY for the "
    f"pair's own gene ({NO_ENTRY!r}). lane-direction found that gap large on its own arm and said so: "
    "of 48 measured increase links, 21 were answerable in their own cell, 9 having no cached track for "
    "their cell at all and 18 having an element whose cached window carried no entry for their own "
    "gene. A count is not a measurement of the thing you want (docs/LESSONS.md), and an eligible link "
    "is not an answerable one"
)

# ---- the floors and the locus convention, imported and not chosen here ---------------------------

POSITIVE_FLOOR = fresh.POSITIVE_FLOOR
LOCUS_FLOOR = fresh.LOCUS_FLOOR
LOCUS_RULE = cell2.INDEPENDENT_LOCUS_RULE
LOCUS_SPAN = cell2.INDEPENDENT_LOCUS_SPAN
NOT_BIOLOGICAL_INDEPENDENCE = dl.NOT_BIOLOGICAL_INDEPENDENCE

FLOORS = {
    "links": POSITIVE_FLOOR,
    "links_imported_from": "genomeos.attribution.fresh.POSITIVE_FLOOR",
    "independent_loci": LOCUS_FLOOR,
    "independent_loci_imported_from": (
        "genomeos.attribution.fresh.LOCUS_FLOOR, which is genomeos.attribution.cell2.POOLED_LOCUS_FLOOR"
    ),
    "neither_chosen_here": True,
    "applied": "to each pooled arm separately, never to the two arms pooled together",
    "gate_implementation": (
        "genomeos.attribution.direction_link.gate, imported and not re-implemented, so the floors and "
        "their wording come from the files that fixed them before either lane existed"
    ),
}

POOLING = (
    "the pooled count of an arm is lane-direction's K562 links of that arm with an answer plus this "
    "lane's HCT116 links of that arm with an answer, and the pooled independent loci are cell2's "
    f"grouping over those pooled links TOGETHER at a span of {LOCUS_SPAN} bp. cell2.group does not look "
    "at the cell, so pooling a second cell into an arm is NOT additive in loci: a link of either cell "
    "that shares its gene with another, or lies within 1 Mb of it on the same chromosome, is the same "
    "locus, and two K562 loci bridged by one HCT116 link become one. Adding links can therefore lower a "
    "locus count as well as raise it, so the pooled loci are counted over the pooled links and never "
    "added up per cell. The grouping is cell2's operational grouping and not established biological "
    "independence"
)

# ---- what this lane may and may not conclude ------------------------------------------------------

PRIOR_EXPOSURE = (
    "HCT116 is NOT a fresh cell for this purpose: its decreases were scored in a39073d, so these links "
    "carry prior exposure. Any later registration that uses them must state that exposure, and may not "
    "describe this cell as held out, fresh or unseen"
)

GATE_PASS = (
    "both pooled arms are at or above both imported floors on their links with an answer: a direction "
    "test may be registered as its own lane, with cell as a stratum. That registration must state "
    "HCT116's prior exposure in the words of PRIOR_EXPOSURE, and must establish the exactly-zero "
    "exclusion this lane did not read before it treats a presence count as an answered count. No "
    "direction is read in this lane and no test is registered or run here"
)

GATE_NO_GO = (
    "a pooled arm is below an imported floor: a NO-GO. It is recorded with the count and the margin it "
    "falls short by, because a reader is owed it and not because a small margin is better than a large "
    "one. No floor is moved, the two arms are not pooled together to reach a floor, no third cache is "
    "opened to widen an arm, no cell, gene or value is substituted, and no direction test is registered "
    "or run. A count short of a floor is reported as a no-go and NEVER as close, promising, nearly "
    "enough, a good start or enough for a pilot"
)

READINGS = (GATE_PASS, GATE_NO_GO)
THERE_IS_NO_THIRD = True

#: The words no summary, result or commit message of this lane may use of a short count. lane-increase
#: forbade them by name, lane-direction carried the prohibition, and this lane carries it again.
FORBIDDEN_OF_A_SHORT_COUNT = dl.FORBIDDEN_OF_A_SHORT_OR_UNDETECTED_OUTCOME

CANNOT_ESTABLISH = (
    "this lane is a count of answer presence and is not a measurement of anything a direction test "
    "would measure. It does not establish that any direction is readable, does not establish that any "
    "increase on knockdown is a repression mechanism, is not a silencer or repressor call, and is not "
    "a verdict on any compiled rule, because no compiled rule was consulted. It does not establish "
    "that HCT116 is a usable second cell: it establishes only how many links of each arm this partial "
    "cache holds an entry for, as an upper bound on how many could be answered. Clearing a count floor "
    "is not a statement that the effects behind the count are large, and the floors are on counts and "
    "not on effect sizes"
)

NO_REQUESTS = (
    "0 model requests, no money, nothing downloaded and no network access. Every input was already on "
    "disk, and the partial cache was written by an earlier run which this lane neither extends nor "
    "re-runs"
)

ORDER = (
    "the order is the whole point and it is fixed here before the run: (1) the eligibility join, by "
    "lane-direction's own committed predicate, blind to every predicted value; (2) the blinded presence "
    "pass over this cache, which can return presence and never a value; (3) the pooled arms against "
    "both imported floors, per arm; (4) the verdict in the words registered above, and nothing beyond "
    "it. No direction is read at any step, and no step is taken out of order to see an outcome first"
)


# ---- the one function that is ever handed a parsed cache record ----------------------------------


class BlindingError(RuntimeError):
    """Raised by the test's poison value when anything but an identity test is applied to it."""


def genes_with_an_answer(record: Any, cell: str) -> frozenset[str]:
    """The gene NAMES this element's cached window holds an entry for on `cell`'s track.

    The only function of this lane that is handed a parsed cache record, and the only operation it
    applies to a cached value is an identity test against None. Its return type carries gene names and
    nothing else, so a predicted value, a predicted direction and a sign have no path out of it. See
    `BLINDING` for the three tests that hold this, and `tests/test_hct116_count.py` for them.
    """
    if not isinstance(record, dict):
        return frozenset()
    genes = record.get("genes")
    if not isinstance(genes, list):
        return frozenset()
    names: set[str] = set()
    for g in genes:
        if not isinstance(g, dict):
            continue
        name = g.get("gene")
        if not isinstance(name, str):
            continue
        by_cell = g.get("by_cell")
        if not isinstance(by_cell, dict):
            continue
        if cell in by_cell and by_cell[cell] is not None:
            names.add(name)
    return frozenset(names)


def element_file(root: Path, chrom: str, element: str) -> Path:
    """Where this partial cache keeps one element: one file per element, under its chromosome."""
    return Path(root) / chrom / f"{element}.json"


def answer_presence(root: Path, chrom: str, element: str, gene: str, cell: str) -> dict[str, Any]:
    """Whether this cache holds an entry for one link's own gene on its own cell track.

    Returns presence and the reason, never a value: `has_answer` is a `bool` and the record's values
    are reachable only through `genes_with_an_answer`, which returns gene names. The path opened is
    reported so every per-element open this lane makes is on the record.
    """
    path = element_file(root, chrom, element)
    if not path.exists():
        return {"element": element, "opened": None, "has_answer": False, "why": NO_FILE}
    record = json.loads(path.read_text())
    present = gene in genes_with_an_answer(record, cell)
    return {
        "element": element,
        "opened": path.as_posix(),
        "has_answer": present,
        "why": PRESENT if present else NO_ENTRY,
    }


# ---- the pooled gate ------------------------------------------------------------------------------


def pooled_gate(
    k562_rows: list[dict[str, Any]], hct116_rows: list[dict[str, Any]], arm: str
) -> dict[str, Any]:
    """One arm's pooled links against both imported floors, by lane-direction's own `gate`.

    The links of the two cells are grouped TOGETHER, which is what cell2's rule requires and what makes
    the pooled locus count something other than the sum of two per-cell counts.
    """
    pooled = [*k562_rows, *hct116_rows]
    imported = dl.gate(pooled, arm)
    return {
        "arm": arm,
        "links_with_an_answer": len(pooled),
        "from_the_primary_cell": len(k562_rows),
        "from_this_lanes_cell": len(hct116_rows),
        "independent_loci": imported["independent_loci"],
        "loci_in_the_primary_cell_alone": len(set(dl.loci_of(k562_rows))) if k562_rows else 0,
        "loci_in_this_lanes_cell_alone": len(set(dl.loci_of(hct116_rows))) if hct116_rows else 0,
        "floors": {"links": POSITIVE_FLOOR, "independent_loci": LOCUS_FLOOR},
        "meets_both_floors": imported["meets_both_floors"],
        "short_by": imported["short_by"],
        "counts_are_presence_not_answered": PRESENCE_IS_AN_UPPER_BOUND,
        "imported_gate": imported,
        "pooling": POOLING,
    }


def verdict(gates: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """The verdict in the words registered above, with no third reading."""
    met = all(g["meets_both_floors"] for g in gates.values())
    return {
        "both_pooled_arms_meet_both_floors": met,
        "reading": GATE_PASS if met else GATE_NO_GO,
        "short_by": {a: gates[a]["short_by"] for a in sorted(gates)},
        "prior_exposure": PRIOR_EXPOSURE,
        "forbidden_of_a_short_count": FORBIDDEN_OF_A_SHORT_COUNT,
        "there_is_no_third_reading": THERE_IS_NO_THIRD,
        "no_direction_read": BLINDING,
        "cannot_establish": CANNOT_ESTABLISH,
    }


def registration() -> dict[str, Any]:
    """Everything this lane fixed before its run, as the result records it."""
    return {
        "lane": "lane-hct",
        "order": ORDER,
        "blinding": BLINDING,
        "presence_call": PRESENCE_CALL,
        "presence_is_an_upper_bound": PRESENCE_IS_AN_UPPER_BOUND,
        "answer_presence_gap": ANSWER_PRESENCE_GAP,
        "cache": CACHE,
        "cache_shape": CACHE_SHAPE,
        "named_not_read_by_lane_direction": NAMED_NOT_READ,
        "inherited_gate_no_go": INHERITED_GATE_NO_GO,
        "inherited_from": INHERITED_FROM,
        "inherited_figures": INHERITED_FIGURES,
        "floors": FLOORS,
        "pooling": POOLING,
        "independent_locus_rule": LOCUS_RULE,
        "not_biological_independence": NOT_BIOLOGICAL_INDEPENDENCE,
        "eligibility": dl.ELIGIBILITY,
        "contest_rule": dl.CONTEST_RULE,
        "binding_wording_r2": dl.registration()["binding_wording_r2"],
        "prior_exposure": PRIOR_EXPOSURE,
        "readings": list(READINGS),
        "there_is_no_third_reading": THERE_IS_NO_THIRD,
        "forbidden_of_a_short_count": FORBIDDEN_OF_A_SHORT_COUNT,
        "cannot_establish": CANNOT_ESTABLISH,
        "no_requests": NO_REQUESTS,
    }
