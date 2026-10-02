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


# ---------------------------------------------------------------------------
# Aggregation: gene-by-gene-by-bin counters, never a materialised pair list
# ---------------------------------------------------------------------------

#: How many pair records are buffered before they are folded into the counters. The buffer bounds
#: the script's resident set: 13.4 million pairs held at once would be 107 MB of indices alone, and
#: the registration's read discipline is that nothing in this lane scales with the pair count.
FLUSH_AT = 2_000_000


def vector_bins(distances: Any) -> Any:
    """The registered decade bin of an array of distances, agreeing with `distance_bin` elementwise.

    `searchsorted(..., side='left') - 1` puts an exact edge in the lower bin, which is the
    half-open `(lower, upper]` the registration fixed; the clip at 0 is the `distance == 0` case,
    which the registration puts in the first bin.
    """
    edges = _np.asarray(DISTANCE_BIN_EDGES_BP, dtype=_np.float64)
    return _np.maximum(_np.searchsorted(edges, distances, side="left") - 1, 0)


@dataclass
class PairCounts:
    """Concordant and total element pairs, counted per gene pair and per distance bin.

    `within` is indexed [gene, bin]; `diff_*` are indexed [bin, gene * genes + gene] so that one
    resample's dyadic weights can be applied as a single matrix-vector product.
    """

    genes: list[str]
    within_concordant: Any
    within_total: Any
    diff_concordant: Any
    diff_total: Any

    @property
    def gene_count(self) -> int:
        return len(self.genes)


def aggregate_pairs(rows: list[RuleRow]) -> PairCounts:
    """Every unordered element pair of the population, folded into gene-by-gene-by-bin counters.

    Nothing here scales with the number of pairs: the pair records are produced one source element
    at a time and folded into fixed-size counters every FLUSH_AT records.
    """
    if _np is None:  # pragma: no cover
        raise RuntimeError("numpy is required to aggregate pairs")
    genes = sorted({row.gene for row in rows})
    index = {gene: i for i, gene in enumerate(genes)}
    order = sorted(range(len(rows)), key=lambda i: rows[i].midpoint)
    gene_of = _np.array([index[rows[i].gene] for i in order], dtype=_np.int64)
    midpoint = _np.array([rows[i].midpoint for i in order], dtype=_np.float64)
    labels = sorted({row.label for row in rows})
    label_index = {label: i for i, label in enumerate(labels)}
    label_of = _np.array([label_index[rows[i].label] for i in order], dtype=_np.int64)

    n_genes, n_bins = len(genes), len(BIN_LABELS)
    cells = n_genes * n_genes * n_bins
    total = _np.zeros(cells, dtype=_np.int64)
    concordant = _np.zeros(cells, dtype=_np.int64)
    buffer_flat: list[Any] = []
    buffer_same: list[Any] = []
    buffered = 0

    def flush() -> None:
        nonlocal buffered
        if not buffer_flat:
            return
        flat = _np.concatenate(buffer_flat)
        same = _np.concatenate(buffer_same)
        total[:] += _np.bincount(flat, minlength=cells)
        concordant[:] += _np.bincount(flat[same], minlength=cells)
        buffer_flat.clear()
        buffer_same.clear()
        buffered = 0

    for i in range(len(order) - 1):
        others = slice(i + 1, len(order))
        bins = vector_bins(_np.abs(midpoint[others] - midpoint[i]))
        low = _np.minimum(gene_of[others], gene_of[i])
        high = _np.maximum(gene_of[others], gene_of[i])
        buffer_flat.append((low * n_genes + high) * n_bins + bins)
        buffer_same.append(label_of[others] == label_of[i])
        buffered += len(order) - i - 1
        if buffered >= FLUSH_AT:
            flush()
    flush()

    total = total.reshape(n_genes, n_genes, n_bins)
    concordant = concordant.reshape(n_genes, n_genes, n_bins)
    diagonal = _np.arange(n_genes)
    within_total = total[diagonal, diagonal, :].copy()
    within_concordant = concordant[diagonal, diagonal, :].copy()
    total[diagonal, diagonal, :] = 0
    concordant[diagonal, diagonal, :] = 0
    return PairCounts(
        genes=genes,
        within_concordant=within_concordant,
        within_total=within_total,
        diff_concordant=_np.ascontiguousarray(
            concordant.reshape(n_genes * n_genes, n_bins).T.astype(_np.float64)
        ),
        diff_total=_np.ascontiguousarray(total.reshape(n_genes * n_genes, n_bins).T.astype(_np.float64)),
    )


def contributing_bins(counts: PairCounts) -> dict[str, Any]:
    """Which decade bins contribute to the matched comparison, by the registered rule, fixed once
    over the FULL population and held through every resample."""
    within = counts.within_total.sum(axis=0)
    different = counts.diff_total.sum(axis=1)
    keep = [
        b
        for b in range(len(BIN_LABELS))
        if within[b] >= MIN_PAIRS_PER_BIN and different[b] >= MIN_PAIRS_PER_BIN
    ]
    all_within = float(within.sum())
    dropped = 1.0 - (float(within[keep].sum()) / all_within if all_within else 0.0)
    return {
        "contributing": keep,
        "contributing_labels": [BIN_LABELS[b] for b in keep],
        "within_gene_pairs_per_bin": {BIN_LABELS[b]: int(within[b]) for b in range(len(BIN_LABELS))},
        "different_gene_pairs_per_bin": {BIN_LABELS[b]: int(different[b]) for b in range(len(BIN_LABELS))},
        "dropped_within_gene_weight": round(dropped, 6),
        "maximum_dropped_weight": MAX_DROPPED_WEIGHT,
        "inconclusive_by_rule": dropped > MAX_DROPPED_WEIGHT,
        "rule": BIN_CONTRIBUTION_RULE,
    }


def _matched(
    within_c: Any, within_t: Any, diff_c: Any, diff_t: Any, keep: list[int]
) -> tuple[float | None, float | None, float]:
    """The two reweighted arms on one distance distribution, plus the weight actually used.

    A contributing bin whose resampled denominator is empty in either arm takes weight 0 and the
    remaining weights are renormalised; the share of resamples in which that happened is reported
    beside the interval, because it is a mechanical necessity of resampling and not a registered
    choice.
    """
    usable = [b for b in keep if within_t[b] > 0 and diff_t[b] > 0]
    mass = float(sum(within_t[b] for b in usable))
    if not usable or mass <= 0:
        return None, None, 0.0
    weights = {b: float(within_t[b]) / mass for b in usable}
    within = sum(weights[b] * float(within_c[b]) / float(within_t[b]) for b in usable)
    different = sum(weights[b] * float(diff_c[b]) / float(diff_t[b]) for b in usable)
    used = mass / float(sum(within_t[b] for b in keep)) if sum(within_t[b] for b in keep) else 0.0
    return within, different, used


def point_estimates(counts: PairCounts, keep: list[int]) -> dict[str, Any]:
    """The primary and the two matched arms over the full population, before any resampling."""
    within_c = counts.within_concordant.sum(axis=0)
    within_t = counts.within_total.sum(axis=0)
    diff_c = counts.diff_concordant.sum(axis=1)
    diff_t = counts.diff_total.sum(axis=1)
    primary_k, primary_n = int(within_c.sum()), int(within_t.sum())
    within, different, used = _matched(within_c, within_t, diff_c, diff_t, keep)
    return {
        "within_gene_concordant_pairs": primary_k,
        "within_gene_pairs": primary_n,
        "within_gene_concordance_primary": (primary_k / primary_n) if primary_n else None,
        "within_gene_concordance_matched": within,
        "different_gene_concordance_matched": different,
        "different_gene_concordant_pairs": int(diff_c.sum()),
        "different_gene_pairs": int(diff_t.sum()),
        "matched_weight_used": round(used, 6),
    }


def resample(
    counts: PairCounts, keep: list[int], draws: int = BOOTSTRAPS, seed: int = SEED
) -> dict[str, Any]:
    """`draws` gene resamples, returning each resample's arms so every ratio shares one resample.

    GENES ARE THE CLUSTERS. A within-gene pair carries its gene's multiplicity; a different-gene
    pair carries the PRODUCT of its two genes', because such a pair is only observed when both of
    its genes are in the sample.
    """
    rng = random.Random(seed)
    n_genes = counts.gene_count
    series: dict[str, list[float]] = {
        "within_primary": [],
        "within_matched": [],
        "different_matched": [],
    }
    empty_bin_resamples = 0
    for _ in range(draws):
        multiplicity = _np.array(cluster_multiplicities(n_genes, rng), dtype=_np.float64)
        within_c = multiplicity @ counts.within_concordant.astype(_np.float64)
        within_t = multiplicity @ counts.within_total.astype(_np.float64)
        dyadic = _np.outer(multiplicity, multiplicity).ravel()
        diff_c = counts.diff_concordant @ dyadic
        diff_t = counts.diff_total @ dyadic
        primary = float(within_c.sum() / within_t.sum()) if within_t.sum() > 0 else None
        within, different, used = _matched(within_c, within_t, diff_c, diff_t, keep)
        if used < 1.0:
            empty_bin_resamples += 1
        if primary is None or within is None or different is None:
            continue
        series["within_primary"].append(primary)
        series["within_matched"].append(within)
        series["different_matched"].append(different)
    return {
        "series": series,
        "resamples_requested": draws,
        "resamples_used": len(series["within_primary"]),
        "resamples_with_an_empty_contributing_bin": empty_bin_resamples,
        "seed": seed,
        "clusters": n_genes,
        "cluster_minimum": MIN_CLUSTERS,
        "met_minimum": draws >= MIN_RESAMPLES and n_genes >= MIN_CLUSTERS,
        "cluster_rule": CLUSTER_RULE,
    }


def ratio_interval(numerator: list[float], denominator: list[float] | float, name: str) -> dict[str, Any]:
    """An interval on a RATIO, from the resamples both arms shared, with its band and assumptions.

    `denominator` is a list when the denominator is itself resampled (arm against arm) and a float
    when it is control 1, which the registration holds fixed. The band, the identical-resample
    share and the instrument's assumptions travel with the interval and are never printed apart
    from it.
    """
    if isinstance(denominator, list):
        if len(numerator) != len(denominator):
            raise ValueError("a ratio of two arms needs both arms from the same resamples")
        values = [n / d for n, d in zip(numerator, denominator, strict=True) if d > 0]
        assumption = "both arms from the SAME gene resample, so their common cluster draw cancels"
    else:
        if denominator <= 0:
            raise ValueError("control 1 must be positive")
        values = [n / denominator for n in numerator]
        assumption = CONTROL_ONE_IS_A_COMMITTED_CONSTANT
    interval = percentile_interval(values)
    return {
        "ratio": name,
        "ci95": interval["ci95"],
        "identical_resample_share": interval["identical_share"],
        "degenerate": interval.get("degenerate"),
        "resamples": interval["resamples"],
        "band": band_of(interval["ci95"]),
        "relative_half_width": relative_half_width(interval["ci95"]),
        "band_half_width": round(BAND_AT_CONTROL[1] - 1.0, 6),
        "wider_than_the_band": (
            None
            if interval["ci95"] is None
            else bool((relative_half_width(interval["ci95"]) or 0.0) > (BAND_AT_CONTROL[1] - 1.0))
        ),
        "instrument_assumptions": (
            "GENES are the clusters; elements within a gene are NOT assumed independent. "
            + assumption
            + ". A percentile bootstrap over whole genes, "
            f"{interval['resamples']} usable resamples, seed {SEED}. The identical-resample share "
            "above is the share of resamples returning the modal value: at 1.0 the interval is "
            "degenerate and decides nothing, and the reading routes to a Wilson bound on the "
            "cluster count."
        ),
    }


# ---------------------------------------------------------------------------
# The SECONDARY arm: agreement with the gene's top GTEx v8 median-TPM tissue
# ---------------------------------------------------------------------------


def gtex_top_tissue(path: str, wanted: set[str]) -> dict[str, Any]:
    """The top median-TPM tissue of each wanted gene symbol, by the registered normalisation rule.

    Only the wanted symbols are kept, so nothing here scales with the 56,000 rows of the matrix. A
    symbol appearing on more than one row is DROPPED and counted, because two rows cannot both be
    that symbol's profile.
    """
    import gzip

    top: dict[str, str] = {}
    seen: dict[str, int] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        handle.readline()
        handle.readline()
        header = handle.readline().rstrip("\n").split("\t")
        columns = header[2:]
        normalised = [normalise_tissue(c) for c in columns]
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            symbol = fields[1]
            if symbol not in wanted:
                continue
            seen[symbol] = seen.get(symbol, 0) + 1
            values = [float(v) for v in fields[2:]]
            top[symbol] = normalised[max(range(len(values)), key=values.__getitem__)]
    duplicated = sorted(s for s, n in seen.items() if n > 1)
    for symbol in duplicated:
        top.pop(symbol, None)
    return {
        "top_tissue": top,
        "gtex_columns": set(normalised),
        "column_count": len(columns),
        "duplicated_symbols_dropped": duplicated,
        "wanted": len(wanted),
        "matched": len(top),
    }


def gtex_arm(rows: list[RuleRow], gtex: dict[str, Any], control1: float) -> dict[str, Any]:
    """SECONDARY: the share of rules whose argmax label is the gene's top GTEx tissue.

    Clustered by gene like every other interval in this lane, and carrying the training-exposure
    caveat, which is the strongest limitation on this whole lane.
    """
    top: dict[str, str] = gtex["top_tissue"]
    columns: set[str] = gtex["gtex_columns"]
    per_gene: dict[str, list[bool]] = {}
    excluded_label = excluded_gene = 0
    chance = []
    counts = label_distribution(rows)
    total_rules = sum(counts.values())
    by_normalised: dict[str, int] = {}
    for label, count in counts.items():
        key = normalise_tissue(label)
        by_normalised[key] = by_normalised.get(key, 0) + count
    for row in rows:
        key = normalise_tissue(row.label)
        if key not in columns:
            excluded_label += 1
            continue
        if row.gene not in top:
            excluded_gene += 1
            continue
        per_gene.setdefault(row.gene, []).append(key == top[row.gene])
        chance.append(by_normalised.get(top[row.gene], 0) / total_rules)
    flags = [flag for gene in sorted(per_gene) for flag in per_gene[gene]]
    k, n = sum(flags), len(flags)
    rng = random.Random(SEED)
    keys = sorted(per_gene)
    values = []
    if len(keys) >= MIN_CLUSTERS:
        for _ in range(BOOTSTRAPS):
            multiplicity = cluster_multiplicities(len(keys), rng)
            num = den = 0
            for index, mult in enumerate(multiplicity):
                if not mult:
                    continue
                arm = per_gene[keys[index]]
                num += mult * sum(arm)
                den += mult * len(arm)
            if den:
                values.append(num / den)
    interval = percentile_interval(values)
    expected = sum(chance) / len(chance) if chance else None
    return {
        "labelled": "SECONDARY. Never the headline, and never read as evidence about biology.",
        "rules_in_the_arm": n,
        "agree": k,
        "share": (k / n) if n else None,
        "genes_in_the_arm": len(keys),
        "clustered_ci95": interval["ci95"] if len(keys) >= MIN_CLUSTERS else None,
        "identical_resample_share": interval["identical_share"],
        "degenerate": bootstrap_degenerate(k, n),
        "why_no_interval": (
            None
            if len(keys) >= MIN_CLUSTERS
            else f"{len(keys)} gene clusters is below the floor of {MIN_CLUSTERS}"
        ),
        "wilson_unclustered_secondary": wilson(k, n),
        "expected_share_if_labels_were_drawn_from_this_population": (
            None if expected is None else round(expected, 6)
        ),
        "how_the_expected_share_is_computed": (
            "the mean over the arm's rules of the population frequency of that gene's top GTEx "
            "tissue as a compiled label. It is the chance level for THIS statistic, reported "
            "beside the observed share; it is not a registered band and decides nothing"
        ),
        "control_1_quoted_beside_it": round(control1, 6),
        "excluded_label_not_a_gtex_tissue": excluded_label,
        "excluded_gene_not_in_the_gtex_matrix": excluded_gene,
        "duplicated_symbols_dropped": gtex["duplicated_symbols_dropped"],
        "training_exposure_caveat": TRAINING_EXPOSURE_CAVEAT,
        "rule": SECONDARY,
    }


# ---------------------------------------------------------------------------
# AMENDMENT_1: the estimand, on the coordinator's direction, still before any count
# ---------------------------------------------------------------------------

TOP_GENES_FOR_CONCENTRATION = 5

AMENDMENT_1 = (
    "AMENDED on the coordinator's direction. DISCLOSED HONESTLY: the direction arrived AFTER the "
    "registration file was committed at 47835d1, so this text is an APPEND and not part of that "
    "commit. It is nonetheless BLIND and pre-count, and git shows why rather than this sentence "
    "asserting it: at the moment this amendment is committed no concordance has been computed, "
    "scripts/label_gene_count.py does not exist in the tree, and no figure of this lane exists "
    "anywhere. An amendment written after the figures are seen is NOT a pre-registration; this "
    "one is written before any figure exists. It is purely ADDITIVE - no band, no floor, no "
    "reading and no threshold of the registration is changed, and nothing is made easier to "
    "reach.\n"
    "THE DEFECT IT ADDRESSES IS THE ESTIMAND, NOT THE UNCERTAINTY. The primary pools PAIRS across "
    "genes, so a gene with k elements contributes k(k-1)/2 of them. At a mean near 27 elements a "
    "handful of large genes can dominate the number. The gene-clustered interval says how UNCERTAIN "
    "the estimate is; it does not change WHAT is estimated. Two things are therefore registered "
    "beside the primary:\n"
    f"(1) THE CONCENTRATION, PRINTED: the share of all within-gene pairs contributed by the top "
    f"{TOP_GENES_FOR_CONCENTRATION} genes by pair count. One number and NO threshold, so a reader "
    "can see at a glance whether the primary is a statement about 180 genes or about five.\n"
    "(2) A GENE-EQUAL-WEIGHT CONCORDANCE, labelled a SENSITIVITY that DECIDES NOTHING: each "
    "gene's own distance-matched concordance averaged with equal weight over genes, and the SAME "
    "equal weighting applied to the matched different-gene arm - each GENE PAIR weighted equally - "
    "so the two arms stay comparable. A unit's distance reweighting uses the global within-gene "
    "bin weights restricted to the contributing bins in which that unit actually has pairs, "
    "renormalised; a unit with pairs in no contributing bin is EXCLUDED and counted. The reading "
    "is still selected by the PAIR-POOLED primary and by nothing else.\n"
    "(3) THE BINDING SENTENCE, registered here so it binds before either number is seen: IF THE "
    "PAIR-POOLED PRIMARY AND THE GENE-EQUAL-WEIGHT SENSITIVITY DISAGREE, THE RESULT MUST SAY SO "
    "IN THOSE WORDS, and may not report only whichever of the two reads better. They DISAGREE "
    "when the registered band of the gene-equal-weight ratio differs from the band of the "
    "pair-pooled ratio against the same control, or when the gene-equal-weight point estimate "
    "falls outside the pair-pooled ratio's own interval. Either condition is enough."
)


def concentration(counts: PairCounts, top: int = TOP_GENES_FOR_CONCENTRATION) -> dict[str, Any]:
    """`AMENDMENT_1` (1): the share of all within-gene pairs contributed by the top genes.

    Printed with NO threshold attached. It tells a reader whether the pair-pooled primary is a
    statement about the genes of the chromosome or about a handful of large ones.
    """
    per_gene = counts.within_total.sum(axis=1)
    order = sorted(range(len(per_gene)), key=lambda g: -int(per_gene[g]))
    total = int(per_gene.sum())
    chosen = order[:top]
    return {
        "within_gene_pairs": total,
        "top_n": top,
        "top_genes": [counts.genes[g] for g in chosen],
        "top_gene_pairs": [int(per_gene[g]) for g in chosen],
        "share_of_pairs_from_the_top_genes": (
            round(float(sum(int(per_gene[g]) for g in chosen)) / total, 6) if total else None
        ),
        "genes_with_at_least_one_pair": int((per_gene > 0).sum()),
        "no_threshold": "printed as a fact about the estimand; no reading turns on it",
    }


def _unit_matched(conc: Any, tot: Any, keep: list[int], weights: Any) -> Any:
    """Each unit's distance-matched concordance: the global within-gene bin weights restricted to
    the contributing bins where that unit has pairs, renormalised. NaN where a unit has none."""
    sub_c = conc[..., keep].astype(_np.float64)
    sub_t = tot[..., keep].astype(_np.float64)
    weight = _np.where(sub_t > 0, weights[None, :], 0.0)
    mass = weight.sum(axis=1)
    rate = _np.divide(sub_c, sub_t, out=_np.zeros_like(sub_c), where=sub_t > 0)
    with _np.errstate(invalid="ignore", divide="ignore"):
        return _np.where(mass > 0, (weight * rate).sum(axis=1) / _np.where(mass > 0, mass, 1.0), _np.nan)


def gene_equal_weight(counts: PairCounts, keep: list[int]) -> dict[str, Any]:
    """`AMENDMENT_1` (2): the equal-weight sensitivity, point estimates and its clustered interval.

    A SENSITIVITY. It decides nothing: the reading is selected by the pair-pooled primary alone.
    """
    within_t_total = counts.within_total.sum(axis=0).astype(_np.float64)[keep]
    weights = within_t_total / within_t_total.sum() if within_t_total.sum() > 0 else within_t_total
    per_gene = _unit_matched(counts.within_concordant, counts.within_total, keep, weights)
    n_genes = counts.gene_count
    pair_c = counts.diff_concordant.T.reshape(n_genes, n_genes, len(BIN_LABELS))
    pair_t = counts.diff_total.T.reshape(n_genes, n_genes, len(BIN_LABELS))
    upper = _np.triu_indices(n_genes, k=1)
    per_pair = _unit_matched(pair_c[upper], pair_t[upper], keep, weights)
    gene_ok = ~_np.isnan(per_gene)
    pair_ok = ~_np.isnan(per_pair)
    rng = random.Random(SEED + 1)
    within_series: list[float] = []
    diff_series: list[float] = []
    for _ in range(BOOTSTRAPS):
        multiplicity = _np.array(cluster_multiplicities(n_genes, rng), dtype=_np.float64)
        gene_weight = _np.where(gene_ok, multiplicity, 0.0)
        if gene_weight.sum() > 0:
            within_series.append(float((gene_weight * _np.nan_to_num(per_gene)).sum() / gene_weight.sum()))
        dyadic = _np.outer(multiplicity, multiplicity)[upper]
        pair_weight = _np.where(pair_ok, dyadic, 0.0)
        if pair_weight.sum() > 0:
            diff_series.append(float((pair_weight * _np.nan_to_num(per_pair)).sum() / pair_weight.sum()))
    return {
        "labelled": "SENSITIVITY (AMENDMENT_1 item 2). It decides nothing.",
        "within_gene_equal_weight": (float(_np.nanmean(per_gene)) if gene_ok.any() else None),
        "different_gene_equal_weight": (float(_np.nanmean(per_pair)) if pair_ok.any() else None),
        "genes_contributing": int(gene_ok.sum()),
        "genes_excluded_no_pairs_in_a_contributing_bin": int((~gene_ok).sum()),
        "gene_pairs_contributing": int(pair_ok.sum()),
        "gene_pairs_excluded_no_pairs_in_a_contributing_bin": int((~pair_ok).sum()),
        "series": {"within": within_series, "different": diff_series},
        "rule": AMENDMENT_1,
    }


def disagreement(pooled: dict[str, Any], equal: dict[str, Any], point: float | None) -> dict[str, Any]:
    """`AMENDMENT_1` (3): whether the pair-pooled primary and the equal-weight sensitivity disagree.

    The words are fixed by the amendment and are printed whenever either condition holds.
    """
    interval = pooled.get("ci95")
    outside = (
        None if (interval is None or point is None) else bool(point < interval[0] or point > interval[1])
    )
    bands_differ = pooled.get("band") != equal.get("band")
    disagree = bool(bands_differ or outside)
    return {
        "bands_differ": bands_differ,
        "equal_weight_point_outside_the_pooled_interval": outside,
        "disagree": disagree,
        "statement": (
            "THE PAIR-POOLED PRIMARY AND THE GENE-EQUAL-WEIGHT SENSITIVITY DISAGREE."
            if disagree
            else "the pair-pooled primary and the gene-equal-weight sensitivity agree: same "
            "registered band, and the equal-weight point lies inside the pooled interval"
        ),
        "the_reading_is_still_the_pooled_one": (
            "the registration selects the reading on the PAIR-POOLED primary and on nothing else; "
            "the sensitivity decides nothing either way"
        ),
    }


def amendment_1_payload() -> dict[str, Any]:
    """AMENDMENT_1 as its own registration file, the way this repository records an amendment.

    The landed registration at 47835d1 is NOT rewritten: an amendment is a separate file naming
    what it amends, which is how gene_identity recorded its three. This one is additive and
    pre-count.
    """
    return {
        "result": "label_gene_registration_amendment_1",
        "date": "2026-10-02",
        "lane": "lane-labelgene",
        "amends": (
            "data/results/label_gene_registration.json as committed at 47835d1, whose design code "
            "is e13992a. Nothing in that file is rewritten or withdrawn; this file adds to it."
        ),
        "amendment_1": AMENDMENT_1,
        "additive_only": (
            "the bands, the floors, the four readings and (c)'s words are exactly as registered at "
            "47835d1 and are repeated here so a reader can check that none moved."
        ),
        "band_at_control": list(BAND_AT_CONTROL),
        "band_far_above": BAND_FAR_ABOVE,
        "min_clusters_genes": MIN_CLUSTERS,
        "reading_c_words": "THE DATA CANNOT TELL",
        "top_genes_for_concentration": TOP_GENES_FOR_CONCENTRATION,
        "the_reading_is_selected_by_the_pooled_primary": (
            "unchanged by this amendment: the registration selects (a), (a'), (b) or (c) on the "
            "PAIR-POOLED within-gene concordance and on nothing else. The gene-equal-weight figure "
            "is a SENSITIVITY and decides nothing."
        ),
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }
