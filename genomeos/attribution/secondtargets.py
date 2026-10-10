# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The model's second targets as a shortlist, held against an independent perturbation (ROADMAP §3,
area I, Next 3; lane-secondtargets, 2026-10-10).

lane-targetchoice (`targetchoice.py`, b746a23 registered, 715d830 result) found that keeping every
protein-coding gene the compiler's bar names at an element ranks the 1,918 held-out ENCODE CRISPRi pairs
in K562 better than keeping one target per element. Those pairs are part of the CRISPRi benchmark the
project has read for weeks, and the off-gate shortlist the roadmap names (the 52 and 67 pairs above a
drop of 0.1) was found on the benchmark's training pairs. This module asks the shortlist's question on
a perturbation the benchmark does not hold.

**The shortlist.** For each compiled element, its second and further targets: the candidate targets
`targetchoice.CANDIDATES` lists (compiled target first, then every protein-coding gene whose predicted
change on deletion reaches `enhancer_target.MIN_EFFECT`, strongest first) other than the compiled
target itself. Nothing is re-tuned; the list is `targetchoice.chromosome_choice`'s.

**The independent perturbation.** Gasperini et al. 2019 (GEO GSE120861), the at-scale K562 CRISPRi
screen with a single-cell RNA readout, tested each candidate enhancer against every gene within 1 Mb.
The CRISPRi benchmark (EPCrisprBenchmark) pools screens and its training file carries 5,299 of
Gasperini's pairs, so a Gasperini pair is independent here only when no row of either benchmark file
names the same gene at an overlapping element (`independence`).

**What is on this machine.** The screen's pair table (every gRNA group by gene pair tested, coordinates
in hg19) is held; its differential-expression table, which holds the readout, is not
(`READOUT`, never fetched). `read_pair_table` reads coordinates, identifiers and the design's selection
label only. The readout is read by `read_readout`, after the registration, and only if the file exists.

In the order written and committed: this code, then the registration
(data/results/secondtargets_registration.json, with the label-free census of the population), then one
run. No model request, no key, no money, no network.
"""

from __future__ import annotations

import bisect
import csv
import gzip
import hashlib
import json
import math
import random
from collections import Counter, defaultdict
from collections.abc import Callable, Iterable, Sequence
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from genomeos.attribution import crispri
from genomeos.attribution import measured_unknown as mu
from genomeos.attribution import targetchoice as tc
from genomeos.predict.enhancer_target import MIN_EFFECT  # noqa: F401 - bound by BOUND, read via globals()

REGISTERED = "2026-10-10"
LANE = "lane-secondtargets"
CELL = "K562"

# ---- inputs
PAIR_TABLE = mu.STORE / "GSE120861_gene_gRNAgroup_pair_table.at_scale.txt.gz"
CHAIN = mu.STORE / "hg19ToHg38.over.chain.gz"
#: the screen's differential-expression table; not on this machine on 2026-10-10 (never fetched)
READOUT = mu.STORE / "GSE120861_all_deg_results.at_scale.txt.gz"
READOUT_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120861/suppl/GSE120861_all_deg_results.at_scale.txt.gz"
BENCHMARK = (crispri.TRAINING, crispri.HELDOUT)

#: the guide group of a candidate enhancer each test is read from; 'second_two' repeats a subset
GUIDE_GROUP = "top_two"
#: design selection labels whose elements entered the screen because a pilot screen found a hit there
SELECTED_FOR_A_HIT = ("not_model_selected_",)
#: the pair table's columns read, by name; every other field is split and dropped unread
PAIR_COLUMNS = (
    "gRNAgroup.chr",
    "gRNAgroup.start",
    "gRNAgroup.stop",
    "gRNAgroup",
    "general_group",
    "start.targetgene",
    "ENSG.targetgene",
    "targetgene_short_name",
)
#: the benchmark columns read, by name (no effect, significance or label column)
BENCHMARK_COLUMNS = (
    "chrom",
    "chromStart",
    "chromEnd",
    "measuredGeneSymbol",
    "measuredGeneEnsemblId",
    "Dataset",
)

# ---- the registered test
#: a control gene is at a matched distance when its distance is within this factor of the shortlist gene's
DISTANCE_FACTOR = 2.0
CALIPER = math.log10(DISTANCE_FACTOR)
#: below this many matched shortlist tests the test is not run (decided before any readout is read)
MIN_MATCHED = 100
#: a screen hit, by the screen's own call: column names as published in the readout table
HIT = {
    "pair_column": ("pairs4merge", "pairs"),
    "p_column": "pvalue.empirical.adjusted",
    "p_below": 0.1,
    "fold_column": "fold_change.transcript_remap",
    "fold_below": 1.0,
}
#: the positive control: TSS-targeting guide groups against the gene whose TSS they target
POSITIVE_CONTROL_MIN = 0.5
BOOTSTRAP = {
    "unit": "chromosome (the matched tests of one chromosome are drawn together)",
    "resamples": crispri.BOOTSTRAPS,
    "seed": 0,
    "interval": "percentile, 2.5 and 97.5",
}

#: the arms of a test at an element overlapping a compiled element
COMPILED, SHORT, OTHER = "compiled_target", "shortlist", "other"

QUESTION = (
    "Is a compiled element's second or further target at the compiler's bar (a protein-coding gene whose"
    " predicted change on deletion reaches enhancer_target.MIN_EFFECT, 0.1 log2, other than the compiled"
    " target; targetchoice.CANDIDATES, nothing re-tuned) more often a measured hit in Gasperini et al."
    " 2019's at-scale K562 CRISPRi screen than the genes the same screen tested at the same element, at a"
    " matched distance, that the shortlist does not name, on the screen's tests that the CRISPRi benchmark"
    " does not hold?"
)
SHORTLIST = (
    "for each compiled element (the compiled programme with its experimental rule lines removed,"
    " signaltrace.PREDICTED_ONLY), the genes targetchoice.chromosome_choice lists as its candidate targets"
    " other than its compiled target(s), ranked as listed (strongest predicted change first): its second"
    " and further targets at the compiler's bar. Read as built by targetchoice; no parameter of it is set"
    " here"
)
POPULATION = (
    "the at-scale screen's tests in its pair table (GSE120861_gene_gRNAgroup_pair_table.at_scale), one per"
    " candidate-enhancer guide group 'top_two' and gene, lifted from hg19 to GRCh38 by"
    " measured_unknown.Chain at liftOver's default (95% of bases); excluded and counted: non-targeting"
    " groups, TSS and other control groups, the 'second_two' repeat groups, elements the design labels"
    " not_model_selected_* (entered because a pilot screen found a hit there), tests that do not lift, and"
    " repeated tests. Then INDEPENDENCE. Then the gene must be protein-coding in GENCODE v50 by its Ensembl"
    " id, and the element must overlap at least one compiled element (half-open, one base)"
)
INDEPENDENCE = (
    "a test is excluded, and counted, when any row of either EPCrisprBenchmark file (training_K562 and"
    " heldout_5_cell_types, every dataset and cell type, read by coordinates, gene and dataset only) names"
    " the same gene (Ensembl id without version, or symbol) at an element overlapping the test's lifted"
    " element by at least one base. The remaining tests are independent of the benchmark pair by pair; the"
    " subset whose element overlaps no benchmark element at all is independent element by element and is"
    " read beside the primary"
)
ARMS = {
    COMPILED: "the gene is a compiled target of an overlapping compiled element",
    SHORT: (
        "the gene is a second or further candidate target of an overlapping compiled element and the"
        " compiled target of none; its rank is its best position among the candidates (2 = second)"
    ),
    OTHER: "the gene is on no overlapping compiled element's candidate list",
}
MATCHING = (
    "each shortlist test is matched to every 'other' test of the same guide group (the same element) whose"
    " distance (hg19 |TSS - element midpoint|, from the pair table, at least 1 bp) is within a factor of"
    " DISTANCE_FACTOR of the shortlist gene's; its control rate is the share of those controls that are"
    " screen hits. A shortlist test with no such control is unmatched and counted. The nearest gene of the"
    " test's element is the independent protein-coding gene the screen tested there at the smallest"
    " distance (ties to the smaller Ensembl id), whatever its arm"
)
HIT_RULE = (
    "a test is a screen hit when the readout table's empirical p-value adjusted by the screen"
    " (pvalue.empirical.adjusted) is below 0.1 and its fold change (fold_change.transcript_remap) is below"
    " 1: the screen's own call, a decrease. Rows join on the pair id 'gRNAgroup:ENSG' (pairs4merge, else"
    " pairs). A test with no readout row or a missing value is not measured: a shortlist test not measured"
    " leaves the comparison, a control not measured leaves its control set, and a shortlist test left with"
    " no measured control leaves; all counted"
)
METRIC = (
    "over the matched shortlist tests that remain: the shortlist's hit rate minus the mean of their control"
    " rates, in percentage points, with a 95% interval from BOOTSTRAP (whole chromosomes resampled). Beside"
    " it and deciding nothing: the shortlist's rate minus the rate of the nearest gene at each test's"
    " element; the same difference on the element-independent subset; the compiled targets' rate at the"
    " same elements; the shortlist's rate by rank (2, 3, 4-5, 6+)"
)
DECISION = (
    "TOO_FEW_PAIRS if fewer than MIN_MATCHED shortlist tests are matched (decided on the label-free census"
    " before any readout is opened); else NOT_RUN_NO_READOUT if the readout table is not on this machine;"
    " else HARNESS_DEFECT on any STOP_RULE condition; else, with [lo, hi] the 95% interval of the primary"
    " difference: MORE_OFTEN if lo > 0, LESS_OFTEN if hi < 0, otherwise NO_DIFFERENCE; TOO_FEW_PAIRS as well"
    " if fewer than MIN_MATCHED matched tests remain measured under HIT_RULE. The prediction holds"
    " on MORE_OFTEN only. The falsifier is NO_DIFFERENCE or LESS_OFTEN: then second targets were not"
    " detected to be measured hits more often than distance-matched genes the shortlist does not name. No"
    " margin beyond the interval is applied"
)
STOP_RULE = (
    "one run, from this registration's committed code, on the files it pins by sha256. Nothing registered"
    " is changed after a readout value is read; a defect found afterwards is reported beside the reading"
    " and corrected only under a new registration. If the readout table's header lacks a registered"
    " column the run stops before any value is read, and a registration amendment changing column names"
    " only, committed before any value is read, is the one change allowed. The run stops with"
    " HARNESS_DEFECT before any test pair's readout is read if: a pinned input moved; the matched"
    " shortlist tests' ids or sha256 differ from the registration's; or, read first, fewer than"
    " POSITIVE_CONTROL_MIN of the screen's TSS guide groups are hits against the gene whose TSS they"
    " target (the readout is then being read wrongly)"
)
OUTCOMES = {
    "TOO_FEW_PAIRS": (
        "too few independent pairs: of the at-scale screen's {total} candidate-enhancer tests in K562,"
        " {excluded} are already in the CRISPRi benchmark (the same gene at an overlapping element) and"
        " {independent} are independent; {coding} of those test a protein-coding gene at an element that"
        " overlaps a compiled element, {short} of them a second or further target, and {matched} of those"
        " have a non-shortlist gene at a matched distance at the same element, fewer than the {min} the"
        " registered test needs; no readout was read"
    ),
    "NOT_RUN_NO_READOUT": (
        "not run, too few independent pairs with a readout: 0 of the {matched} independent Gasperini 2019"
        " K562 tests of a second or further target that have a distance-matched control at the same element"
        " carry a readout on this machine, because the screen's differential-expression table ({readout})"
        " was never fetched; {excluded} of the screen's {total} candidate-enhancer tests were excluded as"
        " already in the CRISPRi benchmark, {independent} are independent, {coding} of those test a"
        " protein-coding gene at a compiled element and {short} a second or further target. The registered"
        " test waits on that one table; nothing was read"
    ),
    "HARNESS_DEFECT": "stopped before any test pair's readout was read: {what}",
    "MORE_OFTEN": (
        "second targets were measured hits more often: of {n} independent Gasperini 2019 K562 tests of a"
        " second or further target of a compiled element, {h} ({r}%) are screen hits, against {rc}% for the"
        " genes the screen tested at the same element at a matched distance that the shortlist does not"
        " name, {d} percentage points [95% CI {lo}, {hi}]; the nearest tested gene at those elements is a"
        " hit for {rn}%, {dn} points from the shortlist [95% CI {nlo}, {nhi}]. A screen hit is a measured"
        " decrease in K562, not function elsewhere"
    ),
    "NO_DIFFERENCE": (
        "no difference detected: of {n} independent Gasperini 2019 K562 tests of a second or further target"
        " of a compiled element, {h} ({r}%) are screen hits, against {rc}% for the genes the screen tested at"
        " the same element at a matched distance that the shortlist does not name, {d} percentage points"
        " [95% CI {lo}, {hi}]; the nearest tested gene at those elements is a hit for {rn}%, {dn} points"
        " from the shortlist [95% CI {nlo}, {nhi}]. No difference detected is not no difference"
    ),
    "LESS_OFTEN": (
        "second targets were measured hits less often: of {n} independent Gasperini 2019 K562 tests of a"
        " second or further target of a compiled element, {h} ({r}%) are screen hits, against {rc}% for the"
        " genes the screen tested at the same element at a matched distance that the shortlist does not"
        " name, {d} percentage points [95% CI {lo}, {hi}]; the nearest tested gene at those elements is a"
        " hit for {rn}%, {dn} points from the shortlist [95% CI {nlo}, {nhi}]"
    ),
}
RAN = ("MORE_OFTEN", "NO_DIFFERENCE", "LESS_OFTEN")
FALSIFIED_BY = ("NO_DIFFERENCE", "LESS_OFTEN")
LIMITATIONS = (
    "the benchmark kept the screen's pairs that passed its own filters (power, expression, distance), so"
    " the independent remainder is what those filters left out: on average less well powered, so its hit"
    " rates are lower than the benchmark's and a difference is harder to detect. Both arms share that",
    "the shortlist names genes AlphaGenome predicts to respond, and a gene that responds in a predicted"
    " track is likely expressed; a single-cell readout detects a change more easily in a well-expressed"
    " gene. The matched control shares the element and the distance, not the expression level",
    "AlphaGenome was trained on ENCODE data that include K562; the screen itself is not in its training"
    " as far as this lane knows, but its cell line is",
    "a screen hit is a measured decrease of a gene's RNA in K562 after CRISPRi at the element; it is not"
    " function in another cell and not a mechanism",
    "the screen's coordinates are hg19 and are lifted; distances are read in hg19 from the screen's own"
    " table",
    "the design's selection label (general_group) is read only to exclude elements selected for a pilot"
    " hit; its other labels may encode pilot-screen status and are not used",
)
DESCRIPTIVES = (
    "the shortlist genome-wide: compiled elements, those with at least one second target, and the size of"
    " each element's shortlist",
    "every exclusion count, the independence split (same gene in the benchmark, element in the benchmark"
    " for another gene, element not in the benchmark) and the benchmark datasets the excluded tests match",
    "per arm: tests, elements, and (after the readout) hits; unmatched shortlist tests; the shortlist by"
    " rank",
)
#: the registration fields the run refuses to differ from: registration field -> constant of this module
BOUND = {
    "question": "QUESTION",
    "shortlist": "SHORTLIST",
    "population": "POPULATION",
    "independence": "INDEPENDENCE",
    "arms": "ARMS",
    "matching": "MATCHING",
    "distance_factor": "DISTANCE_FACTOR",
    "min_matched": "MIN_MATCHED",
    "min_effect": "MIN_EFFECT",
    "guide_group": "GUIDE_GROUP",
    "selected_for_a_hit": "SELECTED_FOR_A_HIT",
    "hit": "HIT",
    "hit_rule": "HIT_RULE",
    "positive_control_min": "POSITIVE_CONTROL_MIN",
    "metric": "METRIC",
    "bootstrap": "BOOTSTRAP",
    "decision": "DECISION",
    "stop_rule": "STOP_RULE",
    "outcomes": "OUTCOMES",
    "limitations": "LIMITATIONS",
    "descriptives_not_verdicts": "DESCRIPTIVES",
}


def registered_fields() -> dict[str, Any]:
    """Every BOUND constant as JSON stores it, so a registration and this module compare equal or not."""
    return {f: json.loads(json.dumps(globals()[c])) for f, c in BOUND.items()}


def ensg_of(text: str) -> str:
    """An Ensembl gene id without prefix or version ('ENSEMBL:ENSG00000108179.12' -> 'ENSG00000108179')."""
    i = text.find("ENSG")
    return text[i:].split(".", 1)[0] if i >= 0 else ""


# ==== the screen's tests, label-free


@dataclass(frozen=True)
class Tested:
    """One test of the screen: a guide group's element (GRCh38, lifted) against one gene."""

    group: str
    chrom: str
    start: int
    end: int
    ensg: str
    symbol: str  # as the pair table names it
    distance: int  # hg19 |TSS - element midpoint|, at least 1

    @property
    def key(self) -> str:
        """The pair table's own pair id, 'gRNAgroup:ENSG'."""
        return f"{self.group}:{self.ensg}"


Lift = Callable[[str, int, int], "tuple[str, int, int] | None"]


def read_pair_table(lines: Iterable[str], lift: Lift) -> tuple[list[Tested], list[Tested], Counter]:
    """The candidate-enhancer tests (POPULATION's first filters), the TSS groups' own-gene tests (the
    positive control), and the exclusion counts. Reads PAIR_COLUMNS by name and nothing else."""
    it = iter(lines)
    head = next(it, "").rstrip("\n").split("\t")
    col = {n: i for i, n in enumerate(head)}
    missing = [n for n in PAIR_COLUMNS if n not in col]
    if missing:
        raise KeyError(f"pair table lacks {missing}")
    ic, i0, i1, ig, isel, itss, ie, isym = (col[n] for n in PAIR_COLUMNS)
    tests: list[Tested] = []
    tss: list[Tested] = []
    seen: set[str] = set()
    counts: Counter = Counter()
    for line in it:
        f = line.rstrip("\n").split("\t")
        if len(f) < len(head):
            counts["short_rows"] += 1
            continue
        counts["rows"] += 1
        group = f[ig]
        if f[ic] == "NTC" or not f[ic].startswith("chr"):
            counts["non_targeting"] += 1
            continue
        suffix = group.split("_", 1)[1] if "_" in group else ""
        ensg = ensg_of(f[ie])
        if suffix == "TSS":
            counts["tss_control"] += 1
            if group[: -len("_TSS")] == f[isym]:
                tss.append(Tested(group, f[ic], int(f[i0]), int(f[i1]), ensg, f[isym], 1))
            continue
        if suffix != GUIDE_GROUP:
            counts["other_guide_group" if suffix == "second_two" else "other_control"] += 1
            continue
        counts["candidate_enhancer_tests"] += 1
        if f[isel].startswith(SELECTED_FOR_A_HIT):
            counts["selected_for_a_pilot_hit"] += 1
            continue
        start, end = int(f[i0]), int(f[i1])
        got = lift(f[ic], start, end)
        if got is None:
            counts["not_lifted"] += 1
            continue
        key = f"{group}:{ensg}"
        if key in seen:
            counts["repeated"] += 1
            continue
        seen.add(key)
        distance = max(1, abs(int(f[itss]) - (start + end) // 2))
        tests.append(Tested(group, got[0], got[1], got[2], ensg, f[isym], distance))
    counts["tests"] = len(tests)
    return tests, tss, counts


# ==== independence from the benchmark


class BenchmarkIndex:
    """Every row of the benchmark files by element, gene and dataset (no effect, label or significance)."""

    def __init__(self) -> None:
        self.rows: dict[str, list[tuple[int, int, str, str, str]]] = defaultdict(list)
        self._starts: dict[str, list[int]] = {}
        self._longest: dict[str, int] = {}
        self.datasets: Counter = Counter()

    @classmethod
    def from_tables(cls, tables: Iterable[Iterable[str]]) -> BenchmarkIndex:
        self = cls()
        for lines in tables:
            for r in csv.DictReader(lines, delimiter="\t"):
                c, s, e, sym, ens, ds = (r[n] for n in BENCHMARK_COLUMNS)
                self.rows[c].append((int(s), int(e), sym, ensg_of(ens), ds))
                self.datasets[ds] += 1
        for c, rows in self.rows.items():
            rows.sort()
            self._starts[c] = [r[0] for r in rows]
            self._longest[c] = max(r[1] - r[0] for r in rows)
        return self

    @classmethod
    def from_files(cls, names: Sequence[str] = BENCHMARK, root: Path = crispri.KNOWLEDGE) -> BenchmarkIndex:
        with ExitStack() as stack:
            return cls.from_tables([stack.enter_context(gzip.open(root / n, "rt")) for n in names])

    def overlapping(self, chrom: str, start: int, end: int) -> list[tuple[int, int, str, str, str]]:
        rows = self.rows.get(chrom)
        if not rows:
            return []
        lo = bisect.bisect_left(self._starts[chrom], start - self._longest[chrom])
        hi = bisect.bisect_left(self._starts[chrom], end)
        return [r for r in rows[lo:hi] if r[1] > start]


def independence(t: Tested, index: BenchmarkIndex) -> tuple[str, list[str]]:
    """'same_gene' (excluded), 'element_other_gene' or 'element_untested', with the matching datasets."""
    hits = index.overlapping(t.chrom, t.start, t.end)
    same = sorted({r[4] for r in hits if (t.ensg and r[3] == t.ensg) or r[2] == t.symbol})
    if same:
        return "same_gene", same
    return ("element_other_gene" if hits else "element_untested"), []


# ==== the shortlist at the screen's elements


@dataclass(frozen=True)
class Row:
    """An independent protein-coding test at an element overlapping a compiled element."""

    test: Tested
    symbol: str  # GENCODE v50
    arm: str
    rank: int | None  # the best position among the candidates, for the shortlist arm
    element_untested: bool  # the element overlaps no benchmark element at all


def arm_of(t: Tested, symbol: str, cl: tc.ChoiceLookup | None) -> tuple[str | None, int | None]:
    """The test's arm (ARMS), or None when no compiled element overlaps its element."""
    els = cl.overlapping(t.start, t.end) if cl is not None else []
    if not els:
        return None, None
    heads = {lk.symbol for el in els for lk in cl.look.links.get(el, [])}
    if symbol in heads:
        return COMPILED, 1
    ranks = [cl.candidates[el].index(symbol) + 1 for el in els if symbol in cl.candidates.get(el, ())]
    return (SHORT, min(ranks)) if ranks else (OTHER, None)


def shortlist_census(lookups: dict[str, tc.ChoiceLookup]) -> dict[str, Any]:
    """The shortlist genome-wide: how many second and further targets each compiled element carries."""
    sizes = []
    for cl in lookups.values():
        for el, cands in cl.candidates.items():
            heads = {lk.symbol for lk in cl.look.links.get(el, [])}
            sizes.append(sum(1 for g in cands if g not in heads))
    return {
        "compiled_elements": len(sizes),
        "with_a_second_target": sum(1 for s in sizes if s > 0),
        "second_and_further_targets": sum(sizes),
        "shortlist_size_histogram": tc.size_histogram(sizes),
    }


def rows_of(
    tests: Sequence[Tested],
    index: BenchmarkIndex,
    coding: dict[str, dict[str, str]],
    lookups: dict[str, tc.ChoiceLookup],
) -> tuple[list[Row], dict[str, Any]]:
    """POPULATION after its first filters: independence, protein-coding, a compiled element. Label-free.

    `coding` is chromosome -> {Ensembl id: GENCODE v50 symbol} of its protein-coding genes."""
    split: Counter = Counter()
    by_dataset: Counter = Counter()
    counts: Counter = Counter()
    rows: list[Row] = []
    for t in tests:
        kind, datasets = independence(t, index)
        split[kind] += 1
        if kind == "same_gene":
            by_dataset[",".join(datasets)] += 1
            continue
        symbol = coding.get(t.chrom, {}).get(t.ensg)
        if symbol is None:
            counts["not_protein_coding_in_gencode_v50"] += 1
            continue
        arm, rank = arm_of(t, symbol, lookups.get(t.chrom))
        if arm is None:
            counts["no_compiled_element_overlaps"] += 1
            continue
        rows.append(Row(t, symbol, arm, rank, kind == "element_untested"))
    counts["rows"] = len(rows)
    arms = Counter(r.arm for r in rows)
    return rows, {
        "independence": dict(split),
        "excluded_by_benchmark_dataset": dict(by_dataset.most_common()),
        "after_independence": dict(counts),
        "arms": {a: arms.get(a, 0) for a in (COMPILED, SHORT, OTHER)},
        "guide_groups_at_compiled_elements": len({r.test.group for r in rows}),
    }


# ==== matching


@dataclass(frozen=True)
class Unit:
    """One matched shortlist test with its controls and its element's nearest gene."""

    test: Row
    controls: tuple[Row, ...]
    nearest: Row
    compiled: tuple[Row, ...]


def _log(d: int) -> float:
    return math.log10(max(1, d))


def matched(rows: Sequence[Row]) -> tuple[list[Unit], int]:
    """MATCHING: the matched shortlist tests, in a fixed order, and how many found no control."""
    by_group: dict[str, list[Row]] = defaultdict(list)
    for r in rows:
        by_group[r.test.group].append(r)
    units: list[Unit] = []
    unmatched = 0
    for g in sorted(by_group):
        rs = sorted(by_group[g], key=lambda r: (r.test.distance, r.test.ensg))
        nearest = rs[0]
        others = [r for r in rs if r.arm == OTHER]
        compiled = tuple(r for r in rs if r.arm == COMPILED)
        for s in (r for r in rs if r.arm == SHORT):
            ctrl = tuple(
                o for o in others if abs(_log(o.test.distance) - _log(s.test.distance)) <= CALIPER + 1e-12
            )
            if ctrl:
                units.append(Unit(s, ctrl, nearest, compiled))
            else:
                unmatched += 1
    return units, unmatched


def population_ids(units: Sequence[Unit]) -> list[str]:
    return [u.test.test.key for u in units]


def ids_sha256(ids: Sequence[str]) -> str:
    return hashlib.sha256("\n".join(ids).encode()).hexdigest()


def rank_band(rank: int | None) -> str:
    if rank is None:
        return "none"
    return str(rank) if rank <= 3 else "4-5" if rank <= 5 else "6+"


def census(
    pair_counts: Counter, rows: Sequence[Row], split: dict[str, Any], units: Sequence[Unit], unmatched: int
) -> dict[str, Any]:
    """The label-free census the registration pins and the run must reproduce."""
    ind = split["independence"]
    total = pair_counts["tests"]
    return {
        "pair_table": dict(pair_counts),
        "candidate_enhancer_tests_read": total,
        "excluded_same_gene_in_benchmark": ind.get("same_gene", 0),
        "independent": ind.get("element_other_gene", 0) + ind.get("element_untested", 0),
        **split,
        "coding_tests_at_compiled_elements": len(rows),
        "shortlist_tests": sum(1 for r in rows if r.arm == SHORT),
        "shortlist_by_rank": dict(Counter(rank_band(r.rank) for r in rows if r.arm == SHORT)),
        "matched_shortlist_tests": len(units),
        "matched_elements": len({u.test.test.group for u in units}),
        "matched_chromosomes": len({u.test.test.chrom for u in units}),
        "matched_element_untested": sum(1 for u in units if u.test.element_untested),
        "unmatched_shortlist_tests": unmatched,
        "controls_per_matched_test": _quantiles([len(u.controls) for u in units]),
    }


def _quantiles(xs: Sequence[int]) -> dict[str, float] | None:
    if not xs:
        return None
    s = sorted(xs)
    return {"min": s[0], "median": s[len(s) // 2], "max": s[-1], "mean": round(sum(s) / len(s), 2)}


# ==== the readout, after the registration


def read_readout(lines: Iterable[str], wanted: set[str]) -> dict[str, bool | None]:
    """HIT_RULE for the wanted pair ids: True a hit, False not, None not measured. Raises KeyError,
    before any value is read, when the header lacks a registered column."""
    it = iter(lines)
    head = next(it, "").rstrip("\n").split("\t")
    col = {n: i for i, n in enumerate(head)}
    pair = next((c for c in HIT["pair_column"] if c in col), None)
    missing = [c for c in (HIT["p_column"], HIT["fold_column"]) if c not in col]
    if pair is None or missing:
        raise KeyError(f"readout header lacks {missing or list(HIT['pair_column'])}")
    ik, ip, ifc = col[pair], col[HIT["p_column"]], col[HIT["fold_column"]]
    out: dict[str, bool | None] = {}
    for line in it:
        f = line.rstrip("\n").split("\t")
        if len(f) <= max(ik, ip, ifc):
            continue
        key = f[ik].strip('"')
        if key not in wanted:
            continue
        try:
            p, fold = float(f[ip]), float(f[ifc])
        except ValueError:
            out[key] = None
            continue
        out[key] = (
            None if math.isnan(p) or math.isnan(fold) else (p < HIT["p_below"] and fold < HIT["fold_below"])
        )
    return out


def positive_control(tss: Sequence[Tested], hits: dict[str, bool | None]) -> dict[str, Any]:
    measured = [hits[t.key] for t in tss if hits.get(t.key) is not None]
    share = sum(measured) / len(measured) if measured else 0.0
    return {"groups": len(tss), "measured": len(measured), "hits": sum(measured), "share": round(share, 4)}


@dataclass(frozen=True)
class Scored:
    chrom: str
    hit: bool
    control_rate: float
    nearest_hit: bool
    element_untested: bool
    rank: int | None


def scored_units(units: Sequence[Unit], hits: dict[str, bool | None]) -> tuple[list[Scored], Counter]:
    """HIT_RULE's handling of unmeasured tests, then one Scored per remaining matched shortlist test."""
    out: list[Scored] = []
    lost: Counter = Counter()
    for u in units:
        h = hits.get(u.test.test.key)
        if h is None:
            lost["shortlist_test_not_measured"] += 1
            continue
        ctrl = [hits.get(c.test.key) for c in u.controls]
        ctrl = [c for c in ctrl if c is not None]
        lost["controls_not_measured"] += len(u.controls) - len(ctrl)
        if not ctrl:
            lost["no_measured_control"] += 1
            continue
        near = bool(hits.get(u.nearest.test.key))
        out.append(
            Scored(u.test.test.chrom, h, sum(ctrl) / len(ctrl), near, u.test.element_untested, u.test.rank)
        )
    return out, lost


def _diff(units: Sequence[Scored], other: Callable[[Scored], float]) -> float:
    return 100 * (sum(u.hit for u in units) - sum(other(u) for u in units)) / len(units)


def chromosome_interval(
    units: Sequence[Scored], other: Callable[[Scored], float], seed: int = 0, n: int = crispri.BOOTSTRAPS
) -> dict[str, Any]:
    """The difference in percentage points with a 95% interval from resampling whole chromosomes."""
    if not units:
        return {"difference": None, "ci95": None, "chromosomes": 0}
    by: dict[str, list[Scored]] = defaultdict(list)
    for u in units:
        by[u.chrom].append(u)
    chroms = sorted(by)
    rng = random.Random(seed)
    diffs = []
    for _ in range(n):
        draw = [u for c in (rng.choice(chroms) for _ in chroms) for u in by[c]]
        diffs.append(_diff(draw, other))
    diffs.sort()
    lo, hi = diffs[int(0.025 * n)], diffs[min(n - 1, int(0.975 * n))]
    return {
        "difference": round(_diff(units, other), 3),
        "ci95": [round(lo, 3), round(hi, 3)],
        "chromosomes": len(chroms),
        "resamples": n,
    }


def rates(units: Sequence[Scored]) -> dict[str, Any]:
    n = len(units)
    return {
        "tests": n,
        "hits": sum(u.hit for u in units),
        "rate_pct": round(100 * sum(u.hit for u in units) / n, 3) if n else None,
        "control_rate_pct": round(100 * sum(u.control_rate for u in units) / n, 3) if n else None,
        "nearest_rate_pct": round(100 * sum(u.nearest_hit for u in units) / n, 3) if n else None,
    }


def decide(matched_n: int, readout_held: bool, ci95: Sequence[float] | None = None) -> str:
    if matched_n < MIN_MATCHED:
        return "TOO_FEW_PAIRS"
    if not readout_held:
        return "NOT_RUN_NO_READOUT"
    if ci95 is None:
        return "HARNESS_DEFECT"
    lo, hi = ci95
    if lo > 0:
        return "MORE_OFTEN"
    if hi < 0:
        return "LESS_OFTEN"
    return "NO_DIFFERENCE"


def reading(kind: str, **numbers: Any) -> str:
    """The registered sentence of an outcome, its numbers filled in and nothing else."""
    return OUTCOMES[kind].format(**numbers)


def census_numbers(c: dict[str, Any]) -> dict[str, Any]:
    """The numbers the label-free outcomes name, from a census."""
    return {
        "total": f"{c['candidate_enhancer_tests_read']:,}",
        "excluded": f"{c['excluded_same_gene_in_benchmark']:,}",
        "independent": f"{c['independent']:,}",
        "coding": f"{c['coding_tests_at_compiled_elements']:,}",
        "short": f"{c['shortlist_tests']:,}",
        "matched": f"{c['matched_shortlist_tests']:,}",
        "min": MIN_MATCHED,
        "readout": READOUT.name,
    }
