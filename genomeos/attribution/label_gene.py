"""Is a compiled rule's `when: cell_type` a property of the ELEMENT, the GENE, or the LOCUS?

A compiled rule line reads `rule <element> activates <gene> { ... when: cell_type = X ... }`.
Syntactically the cell sits on the element's rule, which reads as a property of the element. But
X is the argmax over the scorer's RNA-seq track names: it is the tissue whose predicted
expression OF THE TARGET GENE moved most when the element was deleted. If the label is a
property of the GENE, then the cell context is scoped wrongly in BioLang - a finding about the
language and not about a number.

This module holds the registration text and the statistics. It reads one committed text file,
takes no model request, and spends no money.

WHAT THIS POPULATION CANNOT ANSWER, registered in advance so that no later reading borrows it.
Every element of `data/organisms/human/noncoding_chr21.bio` appears in exactly one predicted
rule, paired with exactly one gene. The SYMMETRIC test - within-ELEMENT label concordance across
one element's different target genes, which is what would directly establish an ELEMENT property
- therefore CANNOT BE RUN on this population at all. A result at the control level here means
"NOT a gene property and not a locus property". It does NOT mean "an element property". That
second claim is a reading this data cannot deliver, and it may not be made from this lane's
numbers in any document.
"""

from __future__ import annotations

import math
import random
import re
from dataclasses import dataclass
from typing import Any

try:  # pragma: no cover - numpy is a hard dependency of the count script only
    import numpy as _np
except ImportError:  # pragma: no cover
    _np = None

QUESTION = (
    "A compiled rule line claims `when: cell_type = X`, which reads as a property of the ELEMENT. "
    "X is the argmax over the model's RNA-seq track names, so it is the tissue whose predicted "
    "EXPRESSION OF THE TARGET GENE changed most when the element was deleted. Is the label a "
    "property of the element, of the gene, or of the locus?"
)

POPULATION = (
    "the committed data/organisms/human/noncoding_chr21.bio, restricted to the rule lines "
    'carrying `predicted "AlphaGenome deletion`. The two rule lines carrying `experimental` are '
    "EXCLUDED by the registration: they are measured CRISPRi rows whose cell is the cell the "
    "assay was run in and not an argmax, so they are not this lane's object. The declared input "
    "is hashed in the registration and the hash is re-checked by the count script."
)

ELEMENT_APPEARS_ONCE = (
    "Verified by this lane rather than taken on trust: within the predicted population every "
    "element id appears in exactly one rule, paired with exactly one gene. The check is in the "
    "count script and it RAISES if it ever fails, because every reading below depends on it."
)

# ---------------------------------------------------------------------------
# Controls
# ---------------------------------------------------------------------------

CONTROL_ONE = (
    "CONTROL 1, the independent-draw control: the concordance expected if each element's label "
    "were drawn independently from chr21's OWN compiled label distribution over this population, "
    "i.e. sum of p squared over labels. It is NOT one over the track count. A label's genome-wide "
    "frequency already reflects properties of its tracks - scale, dynamic range, how often that "
    "track carries the extreme - before any cell information is supposed, so the track count "
    "would be a control against a uniform draw nobody claims. The argmaxcell registration says "
    'this in writing: "that is precisely why the base rate and not the track count is the '
    "control: a high rate of `label == K562` on K562-measured elements is what the genome-wide "
    'label distribution produces on its own, before any cell-type information is supposed." '
    "Control 1 is quoted beside every ratio this lane reports."
)

CONTROL_ONE_IS_A_COMMITTED_CONSTANT = (
    "Control 1 is computed ONCE over the whole population and then held fixed through every "
    "bootstrap resample. It is the stated null, not a second estimate, so it contributes no width "
    "to any interval. This is an ASSUMPTION and it is printed beside every interval: an interval "
    "on a ratio to a fixed denominator is narrower than one that also carried the denominator's "
    "sampling error, and no reading of this lane may be quoted without this sentence."
)

CONTROL_TWO = (
    "CONTROL 2, the DECIDING control for reading (a): pairs of elements of DIFFERENT genes at "
    "MATCHED element-to-element distance. On one chromosome locality and gene identity are "
    "confounded by construction - the elements of one gene ARE near each other - so a raw "
    "within-gene excess over control 1 is equally well explained by a regional effect. Only the "
    "distance-matched different-gene arm separates them."
)

# ---------------------------------------------------------------------------
# The matched-distance rule, fixed here before any concordance is computed
# ---------------------------------------------------------------------------

DISTANCE_IS = (
    "element-to-element distance is the absolute difference of the two elements' MIDPOINTS in base "
    "pairs, each midpoint being (start + end) // 2 of the element block's `locus:` field. Both "
    "elements of every pair are on chr21, so no cross-chromosome pair exists in this population."
)

DISTANCE_BIN_EDGES_BP: tuple[float, ...] = (
    0.0,
    1_000.0,
    10_000.0,
    100_000.0,
    1_000_000.0,
    10_000_000.0,
    math.inf,
)

BIN_LABELS: tuple[str, ...] = (
    "0-1kb",
    "1kb-10kb",
    "10kb-100kb",
    "100kb-1Mb",
    "1Mb-10Mb",
    ">10Mb",
)

MATCHING_RULE = (
    "STRATIFIED EXACT REWEIGHTING at bin resolution, and not a sampled match, because a sampled "
    "match throws away pairs and adds its own sampling noise. The bin edges are decade bins in "
    f"base pairs, fixed here: {[int(e) if e != math.inf else 'inf' for e in DISTANCE_BIN_EDGES_BP]}, "
    "each bin half-open as (lower, upper] except the first, which includes 0. The TOLERANCE is the "
    "bin width itself: two pairs are 'matched in distance' when they fall in the same decade bin, "
    "and no finer claim of matching is made anywhere in this lane. The matched different-gene "
    "concordance is the within-bin different-gene concordance combined with weights equal to the "
    "WITHIN-GENE pair share of each bin, so the two arms are compared on one distance "
    "distribution. The within-gene arm is reported twice: the PRIMARY over all its own bins, and "
    "a reweighted within_matched over the contributing bins with exactly the weights the "
    "different-gene arm is given, and only that second figure is ever compared with the "
    "different-gene arm."
)

MIN_PAIRS_PER_BIN = 10
MAX_DROPPED_WEIGHT = 0.20

BIN_CONTRIBUTION_RULE = (
    f"A bin CONTRIBUTES to the matched comparison only if it holds at least {MIN_PAIRS_PER_BIN} "
    f"within-gene pairs AND at least {MIN_PAIRS_PER_BIN} different-gene pairs over the FULL "
    "population. The set of contributing bins is fixed once on the full population and held "
    "through every resample, so the estimand does not move between resamples. The within-gene "
    "weight of a non-contributing bin is dropped and the remaining weights renormalised, and the "
    "DROPPED WEIGHT is reported. If the dropped weight exceeds "
    f"{MAX_DROPPED_WEIGHT:.2f} the matched comparison is INCONCLUSIVE BY RULE and prints no "
    "reading, because most of the within-gene distance distribution would then have no "
    "different-gene counterpart to compare against."
)

# ---------------------------------------------------------------------------
# The primary statistic
# ---------------------------------------------------------------------------

PRIMARY = (
    "WITHIN-GENE LABEL CONCORDANCE: over the unordered PAIRS of elements that share one target "
    "gene, the share of pairs whose two rules carry the same `when: cell_type` label. Pooled over "
    "genes as (concordant pairs summed over genes) / (pairs summed over genes), so a gene with "
    "many elements carries more pairs - which is why the INTERVAL, and not the point, is clustered "
    "by gene."
)

UNIT_OF_ANALYSIS = (
    "THE UNIT OF ANALYSIS IS GENES, NOT PAIRS. A gene with 27 elements contributes 351 pairs that "
    "share one gene's whole biology; those pairs are not 351 independent observations. A "
    "pair-level interval would be spuriously narrow - the same mechanism that made another lane's "
    "bootstrap anti-conservative. Every interval in this lane resamples whole GENES with "
    "replacement, every power statement is in units of GENES, and a comparison with fewer than "
    "MIN_CLUSTERS gene clusters prints NO interval and is INCONCLUSIVE BY RULE rather than by "
    "judgement."
)

MIN_CLUSTERS = 10
BOOTSTRAPS = 2000
MIN_RESAMPLES = 1000
SEED = 20261002
Z95 = 1.959964

CLUSTER_RULE = (
    "GENES ARE THE CLUSTERS. One resample draws G genes with replacement from the G genes of the "
    "population, G being the number drawn each time. A WITHIN-GENE pair enters with the "
    "multiplicity of its gene. A DIFFERENT-GENE pair enters with the PRODUCT of its two genes' "
    "multiplicities - the dyadic cluster bootstrap - because such a pair is only observed when "
    "both of its genes are in the sample. Elements within a gene are NOT assumed independent; "
    "that is the whole reason for clustering."
)

DEGENERACY_RULE = (
    "At 0 or n successes a cluster bootstrap is DEGENERATE: every resample returns the same rate "
    "and the interval decides nothing. Such a comparison is routed to a WILSON bound on an "
    "EFFECTIVE sample size and the degeneracy is NAMED in the result file. Because a degenerate "
    "arm cannot estimate its own design effect, the effective sample size is taken as the number "
    "of GENE CLUSTERS - the most conservative unit available - and that substitution is printed "
    "beside the bound. The IDENTICAL-RESAMPLE SHARE is reported for EVERY interval, degenerate or "
    "not, so a near-degenerate interval is visible to a reader without this lane inventing a "
    "threshold for 'near'."
)

# ---------------------------------------------------------------------------
# Bands: relative to control 1, never absolute
# ---------------------------------------------------------------------------

BAND_AT_CONTROL = (0.80, 1.25)
BAND_FAR_ABOVE = 1.50

BANDS_ARE_RELATIVE = (
    "EVERY comparison is on a RELATIVE band: a concordance is read as a MULTIPLE of its control, "
    "never as an absolute share. An absolute equivalence band at a low base rate is close to a "
    "vacuous claim - if control 1 is 0.02, an absolute band of +/-0.02 spans a factor of two - and "
    "that defect cost a result. The bands: a ratio is AT CONTROL when its whole interval lies "
    f"inside [{BAND_AT_CONTROL[0]:.2f}, {BAND_AT_CONTROL[1]:.2f}]; it is FAR ABOVE ('>>') when its "
    f"interval's LOWER bound exceeds {BAND_FAR_ABOVE:.2f}; anything else CROSSES A BAND EDGE. "
    "The same two bands are applied to the ratio of the two arms, within_matched / "
    "different_gene_matched, so 'within-gene >> different-gene' and 'within-gene ~ "
    "different-gene' mean exactly what they mean against control 1."
)

READINGS = {
    "a": (
        "GENE PROPERTY. within_matched >> different_gene_matched AND within-gene >> control 1. The "
        "label travels with the target gene and not merely with the neighbourhood. Two sub-cases "
        "are reported and neither denies (a): (a-i) different_gene_matched ~ control 1, a gene "
        "property alone; (a-ii) different_gene_matched >> control 1 as well, the coordinator's "
        "literal chain within-gene >> matched different-gene >> control 1, which is a gene "
        "property WITH a regional component on top. A pure gene property would put the "
        "different-gene arm AT control 1, so requiring the second link of that chain would have "
        "denied the very reading it was written to describe; it is recorded as a sub-case instead."
    ),
    "a_prime": (
        "LOCUS / REGIONAL PROPERTY. within_matched ~ different_gene_matched, and BOTH >> control "
        "1. Neighbouring genes are often co-expressed, so this is a real third outcome and not a "
        "weak form of (a): the label is a property of the neighbourhood, not of the gene - and it "
        "is not evidence of an element property either."
    ),
    "b": (
        "NOT A GENE PROPERTY AND NOT A LOCUS PROPERTY. Both arms ~ control 1. This is an "
        "equivalence result on the relative band and not a failure to reject: the comparison has "
        "excluded any ratio outside the band. It does NOT mean the label is an element property; "
        "see the registration's `cannot_be_answered` clause."
    ),
    "c": (
        "INCONCLUSIVE. Any deciding interval crosses a band edge, or a comparison is inconclusive "
        "by rule (under MIN_CLUSTERS gene clusters, dropped matching weight above "
        "MAX_DROPPED_WEIGHT, or an achieved relative half-width wider than the band itself). Its "
        "words are THE DATA CANNOT TELL. Never 'no information': those are different claims and "
        "eliding them is how overstatements got into the record."
    ),
}

POWER_RULE = (
    "POWER PER COMPARISON, STATED IN GENES AND IN ADVANCE. (1) The number of gene clusters "
    f"available, against the floor of {MIN_CLUSTERS}; under the floor, no interval. (2) The "
    "ACHIEVED relative half-width of the deciding interval, against the band's own half-width "
    f"{BAND_AT_CONTROL[1] - 1.0:.2f} on the high side. If the achieved relative half-width exceeds "
    "that, the comparison cannot separate 'at control' from 'far above' at all and reads (c) by "
    "construction, whatever the point estimate. This is the rule that stops an underpowered null "
    "being presented as a finding. (3) The genes needed for a half-width of "
    f"{BAND_AT_CONTROL[1] - 1.0:.2f} at the observed design effect, so a reader can see how far "
    "short a comparison falls."
)

# ---------------------------------------------------------------------------
# Secondary: GTEx
# ---------------------------------------------------------------------------

SECONDARY = (
    "SECONDARY, AND LABELLED SECONDARY WHEREVER IT APPEARS: the share of rules whose label equals "
    "the target gene's top GTEx v8 median-TPM tissue, over the rules whose label is a GTEx tissue "
    "and whose gene is in the GTEx matrix. The declared input is "
    "data/cache/c5/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz, hashed in "
    "this registration. The label-to-column rule, fixed here: both strings are normalised by "
    "lowercasing, replacing every run of non-alphanumeric characters with a single underscore and "
    "stripping leading and trailing underscores; a label matches a GTEx column when the normalised "
    "forms are equal. A rule whose normalised label matches no column is EXCLUDED from the "
    "denominator and the excluded count is reported, because a non-GTEx label cannot agree or "
    "disagree with a GTEx ranking. Genes are matched on the gct `Description` symbol; a symbol "
    "appearing on more than one row is dropped and the count reported. This arm is clustered by "
    "gene like every other."
)

TRAINING_EXPOSURE_CAVEAT = (
    "THE TRAINING-EXPOSURE CAVEAT, registered in advance and the strongest limitation on this "
    "whole lane: AlphaGenome's GTEx RNA-seq tracks were trained on this same GTEx RNA-seq. "
    "Agreement between a rule's argmax label and the gene's top GTEx tissue may therefore be "
    "RECALL OF SEEN DATA rather than biology, and the secondary arm cannot distinguish the two. "
    "No number from this arm may be read as evidence about biology, and this sentence must travel "
    "with it word for word. It must not be softened, shortened or moved to a footnote in any "
    "later document."
)

# ---------------------------------------------------------------------------
# Model consequence - a proposal only
# ---------------------------------------------------------------------------

MODEL_CONSEQUENCE = (
    "PROPOSAL ONLY, carried for the owner and not recommended by this lane. If reading (a) holds, "
    "the cell context belongs on the GENE rather than on the element's rule: the rule line's "
    "`when: cell_type` would then be asserting of the element a fact about the target gene's "
    "expression profile. This lane makes NO recommendation about adopting that change - adoption "
    "is the owner's decision. Nothing in this lane touches a .bio file, the compiler, the IR, the "
    "grammar or the parser, and no such change is implemented here."
)

# ---------------------------------------------------------------------------
# Resource ceiling
# ---------------------------------------------------------------------------

RSS_CEILING_BYTES = 4 * 1024**3

RSS_RULE = (
    f"The count script checks its own peak resident set size against a ceiling of "
    f"{RSS_CEILING_BYTES // 1024**3} GiB and RAISES - it does not warn - if the ceiling is "
    "reached. It reads ONE committed element table, line by line; it opens no bulk chromosome "
    "archive and calls no cached-element loader, because an inline loop over 93 elements through "
    "such a loader reached 23 GB resident. The pair aggregation is over gene-by-gene-by-bin "
    "counters, never over a materialised list of pairs."
)

NO_SPEND = (
    "0 model requests, no network, no AlphaGenome request of any kind, no money. Two committed "
    "files are read from disk and nothing else."
)

NEVER_WEAKENED = (
    "No test, threshold, pin or registration is weakened, exempted or allowlisted by this lane. No "
    ".bio file, compiler, IR, grammar or parser file is touched."
)


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

_ELEMENT_RE = re.compile(r"^element (\S+) \{")
_LOCUS_RE = re.compile(r"^\s*locus:\s*(\w+):(\d+)-(\d+)\s*$")
_RULE_RE = re.compile(r"^rule (\S+) (?:activates|inhibits) (\S+) \{[^}]*?when: cell_type = ([^;]+);")
_PREDICTED = 'predicted "AlphaGenome deletion'
_EXPERIMENTAL = "evidence: experimental "


@dataclass(frozen=True)
class RuleRow:
    """One compiled predicted rule line: its element, its gene, its argmax label, its midpoint."""

    element: str
    gene: str
    label: str
    chrom: str
    start: int
    end: int

    @property
    def midpoint(self) -> int:
        return (self.start + self.end) // 2


def parse_population(path: str) -> tuple[list[RuleRow], dict[str, int]]:
    """Read one committed .bio text file line by line and return the predicted rule rows.

    The second return value is the line census this lane VERIFIES rather than assumes:
    `rule_lines`, `predicted_rules`, `experimental_rules`, `element_blocks`.
    """
    loci: dict[str, tuple[str, int, int]] = {}
    rows: list[RuleRow] = []
    census = {
        "rule_lines": 0,
        "predicted_rules": 0,
        "experimental_rules": 0,
        "element_blocks": 0,
    }
    current: str | None = None
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("element "):
                match = _ELEMENT_RE.match(line)
                if match:
                    current = match.group(1)
                    census["element_blocks"] += 1
                continue
            if current is not None and "locus:" in line:
                match = _LOCUS_RE.match(line)
                if match:
                    loci[current] = (match.group(1), int(match.group(2)), int(match.group(3)))
                    current = None
                continue
            if not line.startswith("rule "):
                continue
            census["rule_lines"] += 1
            if _EXPERIMENTAL in line:
                census["experimental_rules"] += 1
                continue
            if _PREDICTED not in line:
                continue
            match = _RULE_RE.match(line)
            if match is None:
                raise ValueError(f"a predicted rule line did not parse: {line[:160]!r}")
            census["predicted_rules"] += 1
            element, gene, label = match.group(1), match.group(2), match.group(3).strip()
            where = loci.get(element)
            if where is None:
                raise ValueError(f"predicted rule {element} has no element block locus")
            rows.append(RuleRow(element, gene, label, where[0], where[1], where[2]))
    return rows, census


def verify_one_gene_per_element(rows: list[RuleRow]) -> dict[str, Any]:
    """RAISE if any element carries more than one predicted rule: every reading depends on it."""
    seen: dict[str, set[str]] = {}
    for row in rows:
        seen.setdefault(row.element, set()).add(row.gene)
    many = {element: sorted(genes) for element, genes in seen.items() if len(genes) > 1}
    if many:
        raise AssertionError(
            f"{len(many)} elements carry more than one target gene; the within-element "
            f"symmetric test would become possible and the registration's "
            f"`cannot_be_answered` clause would be wrong: {sorted(many)[:5]}"
        )
    if len(seen) != len(rows):
        raise AssertionError(
            f"{len(rows)} predicted rules over {len(seen)} distinct elements: an element appears "
            "more than once, so the pair construction below would double-count it"
        )
    return {"distinct_elements": len(seen), "rules": len(rows)}


# ---------------------------------------------------------------------------
# Controls and bins
# ---------------------------------------------------------------------------


def label_distribution(rows: list[RuleRow]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        counts[row.label] = counts.get(row.label, 0) + 1
    return counts


def control_one(rows: list[RuleRow]) -> float:
    """Sum of p squared over chr21's own compiled label distribution. NOT one over the tracks."""
    counts = label_distribution(rows)
    total = sum(counts.values())
    if total == 0:
        raise ValueError("control 1 is undefined over an empty population")
    return sum((count / total) ** 2 for count in counts.values())


def distance_bin(distance: float) -> int:
    """The decade bin of a distance in bp, by the edges fixed in the registration."""
    for index in range(len(DISTANCE_BIN_EDGES_BP) - 1):
        lower, upper = DISTANCE_BIN_EDGES_BP[index], DISTANCE_BIN_EDGES_BP[index + 1]
        if (distance == 0.0 and index == 0) or (lower < distance <= upper):
            return index
    raise ValueError(f"distance {distance} fell in no bin")


# ---------------------------------------------------------------------------
# Intervals
# ---------------------------------------------------------------------------


def wilson(k: int, n: int, z: float = Z95) -> list[float] | None:
    """A Wilson 95% interval on k/n, or None when n is 0. Binomial, unclustered."""
    if n <= 0:
        return None
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(max(0.0, centre - half), 6), round(min(1.0, centre + half), 6)]


def percentile_interval(values: list[float]) -> dict[str, Any]:
    """A 2.5/97.5 percentile interval plus the identical-resample share, for any resampled series.

    The identical-resample share is the share of resamples returning the MODAL value. It is
    reported for every interval this lane prints, so a near-degenerate interval is visible
    without this lane inventing a threshold for 'near'.
    """
    if not values:
        return {"ci95": None, "identical_share": None, "resamples": 0}
    ordered = sorted(values)
    lo = ordered[int(0.025 * (len(ordered) - 1))]
    hi = ordered[int(0.975 * (len(ordered) - 1))]
    modal = max(set(ordered), key=ordered.count)
    share = ordered.count(modal) / len(ordered)
    return {
        "ci95": [round(lo, 6), round(hi, 6)],
        "identical_share": round(share, 4),
        "degenerate": share >= 1.0,
        "resamples": len(values),
        "method": "percentile bootstrap resampling whole GENES with replacement",
    }


def bootstrap_degenerate(k: int, n: int) -> bool:
    """At 0 or n successes every cluster resample returns the same rate."""
    return n > 0 and (k == 0 or k == n)


def band_of(interval: list[float] | None) -> str:
    """Which registered relative band an interval on a RATIO falls in."""
    if interval is None:
        return "no_interval"
    lo, hi = interval
    if BAND_AT_CONTROL[0] <= lo and hi <= BAND_AT_CONTROL[1]:
        return "at_control"
    if lo > BAND_FAR_ABOVE:
        return "far_above"
    return "crosses_a_band_edge"


def relative_half_width(interval: list[float] | None) -> float | None:
    if interval is None:
        return None
    return round((interval[1] - interval[0]) / 2.0, 6)


def read_the_comparison(within_vs_diff: str, within_vs_control: str, diff_vs_control: str) -> tuple[str, str]:
    """Which of (a)/(a')/(b)/(c) the three banded intervals read, by the registered rule alone."""
    if "crosses_a_band_edge" in (within_vs_diff, within_vs_control, diff_vs_control):
        return "c", (
            "INCONCLUSIVE - THE DATA CANNOT TELL. At least one deciding interval crosses a band "
            "edge, so the registered rule selects no substantive reading. This is not 'no "
            "information'."
        )
    if "no_interval" in (within_vs_diff, within_vs_control, diff_vs_control):
        return "c", (
            "INCONCLUSIVE BY RULE - THE DATA CANNOT TELL. A deciding interval was not printed, so "
            "no reading is selected."
        )
    if within_vs_diff == "far_above" and within_vs_control == "far_above":
        if diff_vs_control == "far_above":
            return "a", READINGS["a"] + " Sub-case (a-ii) holds."
        return "a", READINGS["a"] + " Sub-case (a-i) holds."
    if within_vs_diff == "at_control" and within_vs_control == "far_above" and diff_vs_control == "far_above":
        return "a_prime", READINGS["a_prime"]
    if within_vs_control == "at_control" and diff_vs_control == "at_control":
        return "b", READINGS["b"]
    return "c", (
        "INCONCLUSIVE - THE DATA CANNOT TELL. The three banded intervals fall in a combination "
        "the registration does not assign to (a), (a') or (b)."
    )


def genes_needed(
    observed_ratio: float, design_effect: float, half_width: float = BAND_AT_CONTROL[1] - 1.0
) -> int | None:
    """Gene clusters needed for a half-width of `half_width` on the ratio, at a measured design
    effect. Normal approximation on the pair-level rate scaled by the design effect; reported so a
    reader can see how far short a comparison falls, and never used to license a reading."""
    if observed_ratio <= 0 or design_effect <= 0:
        return None
    return max(1, math.ceil(Z95 * Z95 * design_effect / (half_width / observed_ratio) ** 2))


def cluster_multiplicities(count: int, rng: random.Random) -> list[int]:
    """One gene resample: draw `count` genes with replacement, return per-gene multiplicities."""
    multiplicity = [0] * count
    for _ in range(count):
        multiplicity[rng.randrange(count)] += 1
    return multiplicity


def peak_rss_bytes() -> int:
    """This process's peak resident set size in bytes, from resource.getrusage."""
    import resource
    import sys

    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(peak) if sys.platform == "darwin" else int(peak) * 1024


def check_rss(ceiling: int = RSS_CEILING_BYTES) -> int:
    """RAISE - never warn - if peak RSS has reached the registered ceiling."""
    peak = peak_rss_bytes()
    if peak >= ceiling:
        raise MemoryError(
            f"peak RSS {peak} bytes reached the registered ceiling of {ceiling} bytes; "
            "the registration says this RAISES and does not warn"
        )
    return peak


def normalise_tissue(name: str) -> str:
    """The registered label-to-GTEx-column rule: lowercase, non-alphanumeric runs to one `_`."""
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def registration() -> dict[str, Any]:
    """The whole registration, as the dict the register script commits before any count."""
    return {
        "result": "label_gene_registration",
        "date": "2026-10-02",
        "lane": "lane-labelgene",
        "registered_in_its_own_commit": (
            "this file is committed ALONE, before any concordance is computed and before "
            "scripts/label_gene_count.py is run. A registration committed beside its result rests "
            "on a sentence; this one rests on the history, and the result names this file's sha."
        ),
        "question": QUESTION,
        "population": POPULATION,
        "element_appears_once": ELEMENT_APPEARS_ONCE,
        "cannot_be_answered": __doc__.split("WHAT THIS POPULATION CANNOT ANSWER")[1].strip(),
        "primary": PRIMARY,
        "unit_of_analysis": UNIT_OF_ANALYSIS,
        "control_1": CONTROL_ONE,
        "control_1_is_a_committed_constant": CONTROL_ONE_IS_A_COMMITTED_CONSTANT,
        "control_2": CONTROL_TWO,
        "distance_is": DISTANCE_IS,
        "distance_bin_edges_bp": [int(e) if e != math.inf else "inf" for e in DISTANCE_BIN_EDGES_BP],
        "distance_bin_labels": list(BIN_LABELS),
        "matching_rule": MATCHING_RULE,
        "bin_contribution_rule": BIN_CONTRIBUTION_RULE,
        "bands_are_relative": BANDS_ARE_RELATIVE,
        "band_at_control": list(BAND_AT_CONTROL),
        "band_far_above": BAND_FAR_ABOVE,
        "readings": READINGS,
        "cluster_rule": CLUSTER_RULE,
        "degeneracy_rule": DEGENERACY_RULE,
        "power_rule": POWER_RULE,
        "min_clusters_genes": MIN_CLUSTERS,
        "bootstraps": BOOTSTRAPS,
        "min_resamples": MIN_RESAMPLES,
        "seed": SEED,
        "secondary": SECONDARY,
        "training_exposure_caveat": TRAINING_EXPOSURE_CAVEAT,
        "model_consequence_proposal_only": MODEL_CONSEQUENCE,
        "rss_rule": RSS_RULE,
        "rss_ceiling_bytes": RSS_CEILING_BYTES,
        "no_spend": NO_SPEND,
        "never_weakened": NEVER_WEAKENED,
        "declared_inputs": {
            "data/organisms/human/noncoding_chr21.bio": (
                "778615512beaa46ddbb02f3cbee8502ceb1ff8d07dd43ce553e086b29e035192"
            ),
            "data/cache/c5/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz": (
                "ee7201ff2f280b0de5657d4b08e9a240362d9757efed7f7bd5dba35a5f8617b8"
            ),
        },
        "declared_inputs_are_rechecked": (
            "the count script recomputes both sha256 digests and RAISES on a mismatch, so a "
            "result can never be read against a different input than the one registered."
        ),
    }
