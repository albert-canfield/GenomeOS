# SPDX-License-Identifier: AGPL-3.0-or-later
"""Does a compiled rule's argmax cell carry cell-type information?

`lane-clause1` established (registration `ba8c41f`, result `b4622ad`) that on its own population of
166 predicted rows the rule's cell is an argmax on every row - 150 by `max_drop`, 16 by `max_rise`,
0 assigned any other way. `direction_v2.ASSIGNED_CELL` is `predicted_coding['tissue']`, and
`predict.enhancer_target.predict_target` sets that field to `max_drop_tissue` or `max_rise_tissue`:
the tissue of the larger of the biggest fall and the biggest rise over every track the scorer
returned, chosen jointly with the gene. So a rule's "cell" is the place an extreme happened to fall
and not a context anyone chose before the values were read.

That leaves the question this module asks and does not answer by assertion: **is the argmax cell
nevertheless informative about where the element acts?** It is answerable because some element-gene
pairs were MEASURED in a named cell. For those pairs the rate at which the compiled cell IS the
measured cell can be set beside the rate at which that same cell label appears over every compiled
rule genome-wide - the base rate. A rate at or near the base rate means the cell label is not
evidence of where the rule acts; a rate well above it means the label carries something.

The asymmetry in `CONFOUND` is registered before any count and is the reason the two readings are
not symmetric: the CRISPRi benchmark chose which elements to test, and it tested elements in regions
active in the cell it was testing, so a positive difference is explained by that selection as
readily as by the model. A difference at or near zero is NOT explained by it, because the selection
pushes the other way. This module can therefore refute the proposition cleanly and can only fail to
refute it, never confirm it.

Nothing here edits a rule, a threshold or a `.bio` file: `direction_v2`, `compile` and
`enhancer_target` are imported and read.

READ DISCIPLINE. `enhancer_target.load_cached` falls back to the per-chromosome archive and caches
every archive it opens in a module-level dict that never evicts; a peer's inline loop over 93
elements reached 23 GB RSS on 2026-10-02 and drove the machine's free disk under the floor. No
function here opens an archive or calls that loader. The compact per-chromosome element table is
read by `stream_table`, which decodes ONE element at a time with `json.JSONDecoder.raw_decode` and
retains four scalars per element, so the largest chromosome's 59 MB table is never held as a parsed
object graph. `check_rss` raises at `RSS_CEILING_BYTES` after every chromosome rather than warning.
"""

from __future__ import annotations

import json
import math
import resource
import sys
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from genomeos.attribution import compile as cp
from genomeos.attribution import crispri as cr
from genomeos.attribution import direction_v2 as dv2
from genomeos.attribution import targets as tg
from genomeos.predict import enhancer_target as et

# ---- the files this lane reads, and only these -------------------------------------------------

#: The committed genome-wide census of every compiled rule's `when: cell_type` label. The base rate
#: comes from here and from nowhere else, so it is a published figure this lane did not compute.
CENSUS = Path("data/results/context_evidence.json")
CENSUS_RELATIVE = "data/results/context_evidence.json"

#: The two cached ENCODE CRISPRi benchmark tables. `crispri.load` parses them; this lane reads a
#: pair's cell, locus and gene and reads no outcome column, no DHS value and no power column.
TRAINING = cr.TRAINING
HELDOUT = cr.HELDOUT

#: The compact per-chromosome element tables the all-enhancer sweep wrote, pointed at by each
#: committed `enhancer_targets_all_<chrom>` summary's `elements_where`. Plain JSON, never an archive.
COMPACT = cr.ELEMENTS

CHROMS = tuple(f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"])

# ---- the question, the population and the denominator ------------------------------------------

QUESTION = (
    "A compiled rule's cell is an argmax over the scorer's tracks: `predicted_coding['tissue']` is "
    "`max_drop_tissue` or `max_rise_tissue`, the place the largest movement happened to fall. Does "
    "that argmax carry cell-type information at all? For element-gene pairs MEASURED in a named "
    "cell, at what rate is the compiled cell label equal to the measured cell, against the rate at "
    "which that label appears over every compiled rule genome-wide?"
)

POPULATION = (
    "One arm per measured cell. An arm's population is every DISTINCT attributed element that "
    "overlaps at least one valid pair of the two cached ENCODE CRISPRi benchmark tables whose "
    "`CellType` is that cell. `attributed` is `genomeos.attribution.targets.attributed`'s own "
    "population - the elements of the three deletion runs that name a coding target gene, "
    "deduplicated by element id with the whole-chromosome run first - which is exactly the "
    "population `compile.compile_chromosome` emits one rule for, and therefore exactly the "
    "population the base rate is taken over. A pair is valid by the benchmark's own "
    "`ValidConnection` column, which `crispri.parse` applies; promoter and exon overlaps are "
    "already dropped there. Overlap is `OVERLAP_RULE`."
)

DENOMINATOR = (
    "DISTINCT ELEMENTS, not pairs. The compiled cell is a property of the element's rule and not of "
    "the pair, so counting an element once per tested gene would weight elements by how many genes "
    "the screen happened to assay near them. The pair-level figure is reported beside it as a "
    "secondary count and is never the headline. Both are stated with the arm's cell name."
)

DENOMINATOR_IS_SEPARATE = (
    "This lane's denominators are its own and are added to nothing. They are NOT to be added to or "
    "compared with the 112, 93, 73, 67, 26, 19, 18, 14 or 10 of tonight's other lanes: several of "
    "those overlap, several are different populations, and two were wrongly added together earlier "
    "today. Every figure below is reported with the arm and the denominator it was taken over."
)

WILL_REPORT = (
    "per arm: the measured cell; the number of valid pairs in that cell; how many of them overlap "
    "an attributed element; the number of DISTINCT attributed elements they overlap (the "
    "denominator); how many of those elements carry the arm's cell as their compiled label (the "
    "numerator); the observed rate; its Wilson 95% interval; the registered base rate; the "
    "difference; and which of the three registered readings the difference selects. Plus, over the "
    "same elements, the ten most frequent compiled labels with each one's own genome-wide base rate "
    "(the CROSS_ARM check), and the cross-arm matrix of every arm's rate on every other arm's cell."
)

# ---- the base rate, and how it is computed -----------------------------------------------------

BASE_RATE_RULE = (
    "base_rate(cell) = census['per_cell'][label]['rules'] / census['rules'], where census is the "
    "COMMITTED data/results/context_evidence.json and label is `compile.context(cell)` - the same "
    "function the compiler used to write the `when: cell_type` value the census counted. The census "
    "denominator is 440,589: every rule `compile.compile_chromosome` emits over the 24 chromosomes, "
    "counted as the compiler emits it, and its own `per_cell` entries sum to exactly that. No part "
    "of this base rate is computed by this lane and none of it is derived from the measured pairs, "
    "so it cannot have been chosen after seeing an observed figure."
)

BASE_RATE_IS_NOT_UNIFORM = (
    "The base rate is NOT 1/371. K562 is the single most frequent compiled label genome-wide, so "
    "its base rate is about twenty-three times what a uniform draw over the scorer's tracks would "
    "give. That is precisely why the base rate and not the track count is the control: a high rate "
    "of `label == K562` on K562-measured elements is what the genome-wide label distribution "
    "produces on its own, before any cell-type information is supposed."
)

SECOND_BASE_RATE = (
    "A second base rate is computed in the SAME streaming pass as the observed rate, from this "
    "lane's own reading of the element tables: the share of all attributed elements of all 24 "
    "chromosomes whose compiled label is the arm's cell. It exists to check the committed census "
    "rather than to replace it. The registered reading uses the COMMITTED base rate; the second is "
    "reported beside it with the difference, and a difference between the two is reported as a "
    "difference and not reconciled away."
)

#: The arms, in the order they are reported. Every one is a cell a benchmark table measured in.
ARMS = ("K562", "GM12878", "HCT116", "Jurkat Clone E6-1", "WTC11")

# ---- the falsifier, stated before any count ----------------------------------------------------

TOLERANCE = 0.02  # absolute, on the difference between the observed rate and the committed base rate
USABLE = 0.25  # the share of an arm's elements that must carry the arm's cell for the label to be usable
MIN_ELEMENTS = 30  # fewer than this and the arm reports counts and NO rate reading

FALSIFIER = (
    "Stated before any count is taken. d = observed_rate - committed_base_rate, both over the arm's "
    "own denominator.\n"
    "(1) d <= -0.02: the argmax lands on the measured cell LESS often than the label's own "
    "genome-wide frequency. The label is not evidence of where the rule acts, and on this arm it "
    "points away from it.\n"
    "(2) |d| < 0.02: NO cell-type information is detected. On this arm the compiled cell is the "
    "label's background frequency and nothing more, so a rule's `cell` is NOT evidence of where it "
    "acts - and that is the finding, stated plainly.\n"
    "(3) d >= +0.02: a difference IS detected. Read with CONFOUND, which says in advance that this "
    "branch does not establish that the argmax carries cell-type information, because the "
    "benchmark's choice of which elements to test produces the same sign.\n"
    "A SECOND, separate threshold, on usability rather than on detection: observed_rate >= 0.25 - "
    "the compiled cell names the measured cell on at least a quarter of the arm's elements, so a "
    "reader could use it; observed_rate < 0.25 - whatever d says, the label names the measured cell "
    "on a minority of elements and may not be read as the place the rule acts.\n"
    "An arm with fewer than 30 elements reports its counts and NO rate reading at all. An interval "
    "is a Wilson 95% interval on the observed rate; whether it excludes the base rate is reported, "
    "and it is a binomial interval over elements, which is NOT a clustered interval and is labelled "
    "as such wherever it appears."
)

CONFOUND = (
    "Registered before any count, and it makes the two readings asymmetric. The CRISPRi benchmark "
    "chose which elements to test, and it tested candidate elements in regions active in the cell "
    "it was testing. So the elements of the K562 arm are not a random draw from the compiled "
    "population: they are enriched for K562-active sequence before any model touched them. That "
    "enrichment pushes the observed rate UP. Therefore:\n"
    "- a difference at or near zero is NOT explained by the selection, because the selection pushes "
    "the other way, and it is clean evidence that the argmax carries no cell-type information;\n"
    "- a positive difference is explained by the selection at least as readily as by the model, and "
    "this lane cannot separate the two with the data it reads. A positive difference is therefore "
    "reported as a difference and MUST NOT be reported as the argmax carrying cell-type "
    "information.\n"
    "This lane can refute the proposition and can fail to refute it. It cannot confirm it."
)

CROSS_ARM = (
    "The one check that bears on the confound without new data, reported for every pair of arms: "
    "the rate of `label == C` over the elements of arm A, for every arm cell C. If the argmax "
    "tracked the cell a measurement was taken in, arm A's own cell would stand out among the "
    "columns of row A. If instead each column is near that column's genome-wide base rate whatever "
    "the row, the compiled label is tracking the genome-wide label distribution and not the "
    "measured cell. This is a comparison of measured rates and sets no threshold of its own."
)

# ---- what a name match is, and is not ----------------------------------------------------------

LABEL_RULE = (
    "An element's compiled label is `compile.context(predicted_coding['tissue'])`, which is "
    "`compile.ident` on a non-empty name and `unknown` on an absent one. The arm's cell is passed "
    "through the same function before comparison, so both sides of the equality are the string the "
    "compiler would write into `when: cell_type`. Two raw tissue names that differ only in "
    "characters `ident` replaces collapse to one label; that is a property of the census this base "
    "rate comes from, so the observed side must share it or the two would not be comparable."
)

NAME_MATCH_IS_NOT_A_CONFIRMATION = (
    "Equality of two labels is equality of two strings. It does not establish that the scorer's "
    "track and the screen's cells are the same material, the same passage or the same laboratory, "
    "and no file this lane reads could establish any of that."
)

METADATA_COPY_LIMITATION = (
    "Carried word for word from the correction a peer established, because it bounds what any cell "
    "statement of this lane may say: the saved copy of the model client's track metadata carries NO "
    "biosample accession and NO experiment accession. `one biosample` is therefore unsupportable "
    "from it; `one biosample NAME` is the most that can be said, and a two-track cell's two tracks "
    "are polyA plus RNA-seq and total RNA-seq of one biosample name, `not known to be biological "
    "replicates` in `_cell_summary`'s own words. This lane reads that metadata copy not at all and "
    "makes no claim that needs it."
)

CELL2_GROUP_CAVEAT = "an operational grouping, NOT established biological independence"

GROUPING = (
    "No grouping is applied. Elements are counted individually and the interval is binomial over "
    f"elements. If a later run groups them, `cell2.group`'s caveat travels with it: {CELL2_GROUP_CAVEAT}."
)

# ---- what this lane may not claim --------------------------------------------------------------

VALIDATES_NOTHING = (
    "This lane validates nothing about the direction rule in either direction. A peer's registered "
    "reading stands and is carried word for word: v2 can withhold a direction and never reverse "
    "one, and an unresolved class is NOT absence of regulation. Nothing measured here is evidence "
    "for or against v1's directions, and no reading below may be cited as such. The question here "
    "is about the CELL a rule names, not about the direction it calls."
)

NO_RECOMMENDATION = (
    "No rule, threshold, clause or `.bio` file is changed by this lane, and none is proposed for "
    "change. This is a measurement. If the finding is that a rule's cell is not evidence of where "
    "it acts, that is a finding for Albert and the supervisor; changing the cell assignment would "
    "invalidate every count already taken under it."
)

CANNOT_ESTABLISH = (
    "whether the scorer's track for a cell name and the screen's cells of that name are the same "
    "material: no file read here carries an accession on either side",
    "whether a compiled label that is NOT the measured cell is wrong: the element may act in both, "
    "and a screen in one cell measures nothing about the other",
    "whether the argmax carries cell-type information when the observed rate is ABOVE the base "
    "rate: the benchmark's selection of tested elements produces the same sign (CONFOUND)",
    "any direction, magnitude or effect: no outcome column of the benchmark tables is read, and no "
    "predicted value is read beyond the tissue field of the compact table",
    "anything about the 98% non-coding space outside the attributed elements: an element with no "
    "predicted coding target emits no rule and is in neither the base rate nor any arm",
)

COUNTS_NAMED = (
    "pairs_in_cell",
    "pairs_on_an_attributed_element",
    "elements",
    "elements_whose_label_is_the_measured_cell",
    "observed_rate",
    "wilson95",
    "committed_base_rate",
    "difference",
    "second_base_rate_this_run",
    "label_distribution_top10",
    "cross_arm_matrix",
    "attributed_elements_streamed",
    "chromosomes_streamed",
    "max_rss_bytes",
)

REFUSALS = (
    "the committed census is absent or untracked: raise. A committed artefact's absence is a "
    "defect, not a difference, and a guard that skips on it would be measuring the machine",
    "the census `per_cell` entries do not sum to its own `rules`: raise, because the base rate's "
    "denominator would then not be the population it claims",
    "a compact element table named by a committed summary is not on this machine: raise. "
    "`targets.run_elements` already raises for this reason and the same rule holds here",
    "resident memory passes RSS_CEILING_BYTES after any chromosome: raise, naming the chromosome "
    "and the figure. Never widen the ceiling to make a run pass",
    "any call that would open a per-element archive or the cached-element loader: there is no such "
    "call, and the test suite asserts the module's source contains neither name",
    "an arm with fewer than 30 elements: counts are reported and NO rate reading is",
)

# ---- the read discipline -----------------------------------------------------------------------

#: 2 GiB. Measured against a gate sample before the full run; the full run's own figure is reported.
RSS_CEILING_BYTES = 2 * 1024**3

READ_DISCIPLINE = (
    "No archive is opened and the cached-element loader is never called. The compact "
    "per-chromosome element table is decoded ONE element at a time by `stream_table` with "
    "`json.JSONDecoder.raw_decode`, and four scalars per element are retained, so the parsed object "
    "graph of a 59 MB table is never held. One chromosome is in hand at a time and its text is "
    "released before the next is opened. `check_rss` raises at RSS_CEILING_BYTES after every "
    f"chromosome. The ceiling is {RSS_CEILING_BYTES} bytes and the run reports its own peak."
)

#: The overlap rule, copied from `crispri.DeletionTable.overlapping` so the arm's elements are the
#: ones that module would have found. `test_argmaxcell` asserts the two agree element for element.
OVERLAP_RULE = (
    "`crispri.DeletionTable.overlapping`'s own rule, reimplemented over four scalars per element "
    "and asserted equal to it by test: an element overlaps [start, end) when element.end > start "
    f"and element.start > start - {cr.REACH}, searched from bisect_left(starts, end) - 1 downwards. "
    f"{cr.REACH} is `crispri.REACH`, the longest element in the table, and is imported rather than "
    "restated."
)

ARGMAX_RULE = getattr(dv2, "ASSIGNED_CELL", "")

ARGMAX_IS_A_SELECTION = (
    "`direction_v2.ASSIGNED_CELL` is `predicted_coding['tissue']`, and "
    "`enhancer_target.predict_target` sets that field to `max_drop_tissue` or `max_rise_tissue` - "
    "the tissue of the larger of the biggest fall and the biggest rise over every track, chosen "
    "jointly with the gene. An argmax is a selection: the cell is where an extreme fell, not a "
    "context chosen before the values were read. `lane-clause1` measured this on 166 of 166 rows of "
    "its own population (150 by max_drop, 16 by max_rise, 0 otherwise), registration ba8c41f, "
    "result b4622ad. That count is carried here, not recomputed, and this lane's denominators are "
    "separate from its 166."
)

EXPLORATION_PRECEDED_THIS = (
    "Stated plainly rather than implied. Before this registration was written this lane read: "
    "`genomeos/predict/enhancer_target.py` for CELLS, the cached-element loader and predict_target; "
    "`genomeos/attribution/direction_v2.py` for the assigned-cell field; "
    "`genomeos/attribution/targets.py` and `compile.py` for the attributed population and the label "
    "function; `genomeos/attribution/crispri.py` for the pair parser, the overlap rule and REACH; "
    "and the committed `data/results/context_evidence.json` for its `per_cell` block. From that "
    "last file it had SEEN the base rates before writing them down: K562 27,445/440,589 = 0.062292, "
    "HepG2 14,779 = 0.033544, GM12878 6,406 = 0.014540, IMR_90 2,942 = 0.006677, HCT116 594 = "
    "0.001348, WTC11 96 = 0.000218, Jurkat 5 = 0.000011. It had seen NO observed rate, NO arm "
    "denominator and NO element of any arm: not one CRISPRi pair had been joined to an element and "
    "not one element table had been opened when this was written. What this registration binds is "
    "the population, the denominator, the base rate and its formula, the falsifier, the confound "
    "and the count list, all of them before the observed side exists."
)


# ---- the base rate ------------------------------------------------------------------------------


def census(path: Path = CENSUS) -> dict[str, Any]:
    """The committed genome-wide label census, with the refusals that make it a base rate applied.

    Absence RAISES. A committed artefact's absence is a defect, not a difference: a guard that
    skipped here would report a pass while measuring nothing, which is the defect class
    `must_be_committed()` was written against at 4f44dbf.
    """
    if not path.exists():
        raise FileNotFoundError(
            f"{path} is absent from this checkout. It is a committed result and the base rate of "
            "this lane is read from it; its absence is a defect and never a reason to skip"
        )
    with path.open() as fh:
        d = json.load(fh)
    per_cell = d.get("per_cell") or {}
    total = d.get("rules")
    if not per_cell or not isinstance(total, int) or total <= 0:
        raise ValueError(f"{path} carries no usable per_cell census (REFUSALS)")
    summed = sum(int(v["rules"]) for v in per_cell.values())
    if summed != total:
        raise ValueError(
            f"{path}: per_cell rules sum to {summed} but the census reports {total} rules, so the "
            "base rate's denominator is not the population it claims (REFUSALS)"
        )
    return d


def base_rate(cell: str, d: dict[str, Any]) -> dict[str, Any]:
    """`BASE_RATE_RULE`, applied. The label is the compiler's own, on both sides of the lookup."""
    label = cp.context(cell)
    row = (d.get("per_cell") or {}).get(label) or {}
    rules = int(row.get("rules") or 0)
    total = int(d["rules"])
    return {
        "cell": cell,
        "label": label,
        "rules_with_this_label": rules,
        "rules_total": total,
        "rate": round(rules / total, 6),
        "how": BASE_RATE_RULE,
        "source": CENSUS_RELATIVE,
        "ontology_term": row.get("ontology_term"),
    }


# ---- the streaming element reader ---------------------------------------------------------------

_DECODER = json.JSONDecoder()
_WS = " \t\r\n,"


def _results_dir() -> Path:
    from genomeos.results import RESULTS_DIR

    return RESULTS_DIR


def table_path(chrom: str, results_dir: Path | None = None) -> Path:
    """The compact element table a committed `enhancer_targets_all_<chrom>` summary points at.

    A summary whose table is not on this machine RAISES, which is `targets.run_elements`'s own rule.
    """
    from genomeos.results import load_result

    rd = results_dir or _results_dir()
    summary = load_result(f"enhancer_targets_all_{chrom}", rd) or {}
    where = summary.get("elements_where")
    if not where:
        raise FileNotFoundError(
            f"enhancer_targets_all_{chrom} names no elements_where, so the compact table cannot be "
            "located; this lane will not fall back to an archive (READ_DISCIPLINE)"
        )
    p = Path(where)
    if not p.is_absolute() and not p.exists():
        p = rd.parent.parent / where
    if not p.exists():
        raise FileNotFoundError(
            f"enhancer_targets_all_{chrom} points at {where}, which is not on this machine; rerun "
            "the job. This lane will not read the chromosome archive instead (READ_DISCIPLINE)"
        )
    return p


def stream_table(path: Path) -> Iterator[tuple[str, int, int, str]]:
    """(id, start, end, compiled label) for every element of one compact table that names a gene.

    One element is decoded at a time and four scalars are retained; the element's own dict is
    released before the next is decoded. No archive is opened and no cached-element loader is used.
    """
    text = path.read_text()
    i = text.find("[")
    if i < 0:
        raise ValueError(f"{path} is not a JSON list of elements")
    i += 1
    n = len(text)
    while True:
        while i < n and text[i] in _WS:
            i += 1
        if i >= n or text[i] == "]":
            return
        e, i = _DECODER.raw_decode(text, i)
        pc = e.get("predicted_coding") or {}
        if not pc.get("gene"):
            continue
        yield (str(e["id"]), int(e["start"]), int(e["end"]), cp.context(pc.get("tissue")))


def small_run_elements(chrom: str, results_dir: Path | None = None) -> Iterator[tuple[str, int, int, str]]:
    """The two sampled runs' attributed elements, which are inline in their committed summaries.

    `targets.RUNS` order is kept: the whole-chromosome run first (streamed by `stream_table`), then
    `constrained_targets`, then `enhancer_targets`. Deduplication by id is the caller's, exactly as
    `targets.attributed` does it, so the population here is that function's population.
    """
    for name in tg.RUNS[1:]:
        for e in tg.run_elements(name, chrom, results_dir or _results_dir()):
            pc = e.get("predicted_coding") or {}
            if not pc.get("gene"):
                continue
            yield (str(e["id"]), int(e["start"]), int(e["end"]), cp.context(pc.get("tissue")))


def attributed_table(chrom: str, results_dir: Path | None = None) -> list[tuple[int, int, str, str]]:
    """`targets.attributed(chrom)`'s population as (start, end, id, label), sorted by start.

    Four scalars per element and nothing else, so a chromosome of 45,000 elements costs a few
    megabytes instead of the table's 59.
    """
    seen: set[str] = set()
    out: list[tuple[int, int, str, str]] = []
    for eid, start, end, label in stream_table(table_path(chrom, results_dir)):
        if eid in seen:
            continue
        seen.add(eid)
        out.append((start, end, eid, label))
    for eid, start, end, label in small_run_elements(chrom, results_dir):
        if eid in seen:
            continue
        seen.add(eid)
        out.append((start, end, eid, label))
    out.sort(key=lambda t: t[0])
    return out


def overlapping(
    table: list[tuple[int, int, str, str]], starts: list[int], start: int, end: int
) -> list[tuple[int, int, str, str]]:
    """`OVERLAP_RULE`: `crispri.DeletionTable.overlapping`'s rule over the lean table."""
    import bisect

    out = []
    j = bisect.bisect_left(starts, end) - 1
    while j >= 0 and table[j][0] > start - cr.REACH:
        if table[j][1] > start:
            out.append(table[j])
        j -= 1
    return out


# ---- resident memory ----------------------------------------------------------------------------


def max_rss_bytes() -> int:
    """Peak resident set size of this process, in BYTES on every platform.

    macOS reports `ru_maxrss` in bytes and Linux in kilobytes; the unit is normalised here rather
    than left to the caller, because a ceiling compared against the wrong unit is a ceiling that
    never fires.
    """
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(raw) if sys.platform == "darwin" else int(raw) * 1024


def check_rss(where: str, ceiling: int = RSS_CEILING_BYTES) -> int:
    """The peak, or a RuntimeError naming it. Raises; never warns and never widens the ceiling."""
    peak = max_rss_bytes()
    if peak > ceiling:
        raise RuntimeError(
            f"resident memory reached {peak} bytes after {where}, over the registered ceiling of "
            f"{ceiling}. The run stops. Widening the ceiling to make it pass is not a fix "
            "(REFUSALS, READ_DISCIPLINE)"
        )
    return peak


# ---- the arms -----------------------------------------------------------------------------------


def wilson(k: int, n: int, z: float = 1.959964) -> list[float] | None:
    """A Wilson 95% interval on k/n, or None when n is 0. Binomial over elements, NOT clustered."""
    if n <= 0:
        return None
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(max(0.0, centre - half), 6), round(min(1.0, centre + half), 6)]


INTERVAL_IS_BINOMIAL = (
    "a Wilson 95% interval over ELEMENTS, assuming independent elements. It is NOT a clustered "
    "interval: elements of one locus are not independent, so it is narrower than a chromosome- or "
    "locus-clustered interval would be. It is reported because the readings turn on a difference "
    "from a base rate and a reader is entitled to the sampling width; it is labelled here and "
    "wherever it is printed, and no reading of this lane rests on it alone."
)


@dataclass
class Arm:
    """One measured cell's arm: its pairs, its elements and the counts registered in COUNTS_NAMED.

    `element_label` is the arm's whole population: one entry per DISTINCT attributed element that
    any of the arm's pairs overlaps, holding that element's compiled label. Every count below is
    derived from it, so the element denominator and the label tally cannot disagree.
    """

    cell: str
    label: str
    pairs: int = 0
    pairs_on_an_element: int = 0
    pair_level_matching: int = 0
    element_label: dict[str, str] = field(default_factory=dict)

    def add(self, elements: list[tuple[int, int, str, str]]) -> bool:
        """Record one pair's overlapping elements; True when any of them carries the arm's label."""
        self.pairs += 1
        if not elements:
            return False
        self.pairs_on_an_element += 1
        hit = False
        for _start, _end, eid, label in elements:
            self.element_label[eid] = label
            hit = hit or label == self.label
        if hit:
            self.pair_level_matching += 1
        return hit

    @property
    def elements(self) -> int:
        return len(self.element_label)

    @property
    def elements_matching(self) -> int:
        return sum(1 for v in self.element_label.values() if v == self.label)

    @property
    def labels(self) -> dict[str, int]:
        out: dict[str, int] = {}
        for v in self.element_label.values():
            out[v] = out.get(v, 0) + 1
        return out

    @property
    def rate(self) -> float | None:
        return round(self.elements_matching / self.elements, 6) if self.elements else None


def reading(arm: Arm, base: float) -> dict[str, Any]:
    """`FALSIFIER`, applied to one arm's measured figures. No branch is chosen by hand."""
    n, k = arm.elements, arm.elements_matching
    if n < MIN_ELEMENTS:
        return {
            "rate_reported": False,
            "why": (
                f"{n} elements is below the registered floor of {MIN_ELEMENTS}, so this arm reports "
                "its counts and no rate reading. The counts stand; the reading does not exist"
            ),
            "counts_only": {"elements": n, "elements_whose_label_is_the_measured_cell": k},
        }
    rate = k / n
    d = rate - base
    ci = wilson(k, n)
    if d <= -TOLERANCE:
        detected = (
            "the argmax lands on the measured cell LESS often than that label's own genome-wide "
            "frequency. On this arm the compiled cell is not evidence of where the rule acts, and "
            "it points away from it (FALSIFIER branch 1)"
        )
        carries: bool | None = False
    elif abs(d) < TOLERANCE:
        detected = (
            "NO cell-type information is detected. On this arm the rate at which the compiled cell "
            "is the measured cell is that label's background frequency and nothing more, so a "
            "rule's `cell` is NOT evidence of where it acts (FALSIFIER branch 2). The confound "
            "pushes the other way, so this reading is not explained by the benchmark's selection "
            "of tested elements"
        )
        carries = False
    else:
        detected = (
            "a difference IS detected above the registered tolerance. It does NOT establish that "
            "the argmax carries cell-type information: the benchmark chose elements active in the "
            "cell it was testing, which produces this sign on its own, and this lane cannot "
            "separate the two (FALSIFIER branch 3, CONFOUND)"
        )
        carries = None
    return {
        "rate_reported": True,
        "observed_rate": round(rate, 6),
        "committed_base_rate": round(base, 6),
        "difference": round(d, 6),
        "tolerance": TOLERANCE,
        "wilson95": ci,
        "wilson95_is": INTERVAL_IS_BINOMIAL,
        "wilson95_excludes_the_base_rate": None if ci is None else not (ci[0] <= base <= ci[1]),
        "detection": detected,
        "argmax_carries_cell_type_information": carries,
        "usable": rate >= USABLE,
        "usability": (
            f"observed rate {rate:.4f} is at or above the registered {USABLE}: the compiled cell "
            "names the measured cell on at least a quarter of this arm's elements"
            if rate >= USABLE
            else f"observed rate {rate:.4f} is below the registered {USABLE}: whatever the "
            "difference says, the compiled cell names the measured cell on a MINORITY of this "
            "arm's elements and may not be read as the place the rule acts"
        ),
    }


def top_labels(arm: Arm, d: dict[str, Any], n: int = 10) -> list[dict[str, Any]]:
    """The arm's most frequent compiled labels, each beside its own genome-wide base rate."""
    total = arm.elements or 1
    out = []
    for label, count in sorted(arm.labels.items(), key=lambda kv: (-kv[1], kv[0]))[:n]:
        row = (d.get("per_cell") or {}).get(label) or {}
        bt = int(row.get("rules") or 0) / int(d["rules"])
        out.append(
            {
                "label": label,
                "elements": count,
                "share_of_arm": round(count / total, 6),
                "genome_wide_base_rate": round(bt, 6),
                "ratio_to_base_rate": round((count / total) / bt, 3) if bt else None,
            }
        )
    return out


#: Names that would mean an archive had been opened: the cached-element loader, the archive helper
#: it calls, and the decompression module either would need. `names_used` asserts their absence from
#: this module's own CODE rather than from its text, so the sentence in the docstring above that has
#: to name the loader in order to say it is not used does not make the check pass or fail.
ARCHIVE_NAMES = ("load" + "_cached", "_archive", "gzip", "GzipFile")

_SOURCE = Path(__file__)


def names_used(source: str | None = None) -> set[str]:
    """Every identifier this module's code refers to: imports, names and attribute names.

    Read with `ast`, so a string literal and a docstring are not identifiers. A text scan would have
    to decide whether the docstring sentence that names the loader in order to say it is never
    called counts as a use, and whichever way it decided it would be measuring prose.
    """
    import ast

    tree = ast.parse(source if source is not None else _SOURCE.read_text())
    out: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute):
            out.add(node.attr)
        elif isinstance(node, ast.Name):
            out.add(node.id)
        elif isinstance(node, ast.Import):
            for a in node.names:
                out.add(a.name.split(".")[0])
                out.add((a.asname or a.name).split(".")[0])
        elif isinstance(node, ast.ImportFrom):
            out.add((node.module or "").split(".")[0])
            for a in node.names:
                out.add(a.name)
                out.add(a.asname or a.name)
    return out


def source_opens_no_archive() -> bool:
    """True when this module's code refers to no archive loader. Asserted by test, not claimed."""
    return not (names_used() & set(ARCHIVE_NAMES))


# Named so a reader can see the module knows what it is not using. et is imported for CELLS only.
RETAINED_CELLS_TODAY = tuple(et.CELLS)
RETENTION_IS_NOT_THIS_QUESTION = (
    "`enhancer_target.CELLS` is the four cell lines whose own track value is kept per gene. It "
    "bounds which cells a `by_cell` value exists for and so which rows clause (1) can reach - a "
    "peer's registered reading. It does NOT bound this question: the argmax tissue is over every "
    "track the scorer returned, so a compiled label may be any of the census's 317 and HCT116's 594 "
    "rules are assessable here although HCT116 carries no `by_cell` column at all."
)
