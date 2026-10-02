"""Reads per million mapped reads in an element, computed from alignments by its own definition.

This is the quantity `DHS.RPM` and `H3K27ac.RPM` are, and that is the whole reason this module is
allowed to exist. The AstroREG no-go refused a narrowPeak `signalValue` as a SUBSTITUTE for reads per
million: different units, and censored outside called peaks, so an element lying in no peak could
only be handed an invented zero. An RPM counted from the alignments is not a substitute. It is the
definition, and an element with no reads over it gets a MEASURED zero, which is a number about the
cell rather than a number about the file format. Everything here is registered in
data/results/astroreg_calibration_registration.json before any read was counted.

Nothing in this module opens a network connection or knows what a BAM is. It takes read spans and
flags, which `scripts/astroreg_rpm_calibrate.py` feeds it from a streaming pysam iterator, so the
counting rule can be tested exhaustively on handfuls of synthetic reads with no file and no network.
That separation is the point: the part that is easy to get wrong is the part that is cheap to test.

THE DENOMINATOR IS COUNTED HERE TOO, in the same pass and under the same filters as the numerator.
That is not a convenience. A numerator from this pass divided by a denominator from somebody else's
filters is wrong by exactly the ratio between the two filter sets, and no test of either number on its
own would show it: both would look reasonable. So `Counter.denominator` is incremented by the same
`keep` that admits a read to the numerator, and the published metadata figure is carried separately as
`METADATA_IS_A_CROSS_CHECK` says, printed beside it and never substituted for it.
"""

from __future__ import annotations

import bisect
import math
from collections.abc import Iterable
from typing import Any, Protocol

#: The read rule, as registered. Quoted here so the code and the registration can be compared.
READ_RULE = (
    "a read is counted when it is not unmapped, not secondary and not supplementary. Duplicates are "
    "whatever the filtered BAM delivers, because excluding them would be a rule the benchmark does "
    "not state. A paired-end fragment is counted ONCE, by counting read1 only. Overlap means any "
    "base: a read overlaps an element when its reference span intersects the element's span by at "
    "least one base, both half-open"
)

#: Why no MAPQ threshold is applied here.
NO_THRESHOLD_OF_OUR_OWN = (
    "the benchmark states no MAPQ threshold, so none is applied. ENCODE's filtered alignments already "
    "carry that experiment's own pipeline thresholds, and taking them as given is the one choice that "
    "adds nothing of this lane's invention. A threshold chosen here would be presented as the "
    "published definition while being ours"
)

#: Why the metadata read count is never the denominator.
METADATA_IS_A_CROSS_CHECK = (
    "ENCODE publishes a mapped-read total in each alignment's quality metrics. It is printed beside "
    "the computed denominator with their ratio, and it is NEVER used as the denominator, because it "
    "was counted under ENCODE's filters and the numerator here is counted under these. Two filter "
    "sets silently divided into one another is an error that looks like a result"
)

#: The registered term, carried word for word from astroreg2_registration at 0f4c372, gate.what_is_computed.
NO_RESCALING_OF_ANY_KIND = (
    "the registration requires the activity be reconstructed under amendment 2's rule 'and NO RESCALING "
    "of any kind'. That is a registered term and it is carried here word for word. It is enforced "
    "structurally rather than promised: `Counter.denominator` is a property with no setter, so the "
    "published mapped-read total cannot be put there; it is fetched by the runner, printed BESIDE the "
    "computed denominator with their ratio, and has no path into any number. A ratio of 1.0 is a CHECK "
    "that the two agree and is not a scaling by one of them, and those are different things"
)

#: A measured zero is not a censored absence. The sentence the whole amendment turns on.
A_MEASURED_ZERO = (
    "an element with no reads over it has RPM exactly 0, and that 0 is MEASURED: the alignments were "
    "read across that span and nothing was there. It is categorically different from the zero a peak "
    "call forces on an element outside any called peak, which is the absence of a measurement being "
    "written down as if it were one. The first is admissible in the frozen feature and the second is "
    "what the no-go refused"
)


#: Why one BAM per biological replicate is checked and not merely intended.
ONE_BAM_PER_REPLICATE = (
    "'released' does not mean 'one per replicate'. An ENCODE experiment keeps every pipeline version "
    "it has ever run, so ENCSR000AKP publishes seven released filtered GRCh38 alignment BAMs that are "
    "three biological replicates reprocessed three times. Pooling them would count two replicates "
    "three times each, and because the duplicated reads inflate the numerator and the denominator "
    "together it would reweight the replicate mixture without making any single number look absurd. "
    "So the set is CHECKED, not just described: a replicate appearing twice is refused"
)


def check_one_bam_per_replicate(bams: list[dict[str, Any]]) -> dict[int, str]:
    """Refuse a pooling set in which any biological replicate appears more than once.

    `bams` are dicts with `accession` and `biological_replicates`. Returns the replicate-to-accession
    map when the set is sound, so the caller records which replicate each pass came from.

    A BAM naming no replicate is refused too: it cannot be shown not to duplicate another, and a file
    whose provenance is unknown is exactly what this check exists to keep out of a pooled total.
    """
    seen: dict[int, str] = {}
    for bam in sorted(bams, key=lambda b: str(b.get("accession"))):
        reps = bam.get("biological_replicates")
        if not reps:
            raise ValueError(
                f"{bam.get('accession')} names no biological replicate, so it cannot be shown not to "
                f"duplicate another file in the pool. {ONE_BAM_PER_REPLICATE}"
            )
        for rep in reps:
            if rep in seen:
                raise ValueError(
                    f"biological replicate {rep} appears in both {seen[rep]} and "
                    f"{bam.get('accession')}: pooling them would count the same reads twice. "
                    f"{ONE_BAM_PER_REPLICATE}"
                )
            seen[rep] = bam["accession"]
    if not seen:
        raise ValueError("no BAM to pool")
    return seen


class Read(Protocol):
    """What the counter needs of a read. pysam's AlignedSegment satisfies this."""

    reference_name: str | None
    reference_start: int
    is_unmapped: bool
    is_secondary: bool
    is_supplementary: bool
    is_paired: bool
    is_read1: bool

    @property
    def reference_end(self) -> int | None: ...


def keep(read: Read) -> bool:
    """The registered read rule, and nothing else. See `READ_RULE`."""
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return False
    # a paired fragment is counted once, by its read1; single-ended data has no read2 to drop
    return not (read.is_paired and not read.is_read1)


class Counter:
    """Counts, in one pass, the reads over each element and the reads passing the filter at all.

    Elements are keyed by `(chrom, start, end)`, half-open, and an element given twice is one element:
    the activity columns are a property of the element, so counting it twice would double a value
    rather than add one.

    The overlap search is a bisect over each chromosome's elements sorted by start, walked back by at
    most the longest element on that chromosome. That bound is computed from the elements themselves,
    so it cannot be wrong for a longer element than the author imagined.
    """

    def __init__(self, elements: Iterable[tuple[str, int, int]]) -> None:
        self.counts: dict[tuple[str, int, int], int] = {}
        self._by_chrom: dict[str, tuple[list[int], list[tuple[int, int, int]]]] = {}
        self._reach: dict[str, int] = {}
        self._denominator = 0
        self.reads_seen = 0
        self.reads_rejected = 0
        self.reads_on_unknown_chrom = 0
        by: dict[str, list[tuple[int, int]]] = {}
        for chrom, start, end in elements:
            if end <= start:
                raise ValueError(f"an element must be non-empty and half-open: {(chrom, start, end)}")
            key = (chrom, start, end)
            if key in self.counts:
                continue
            self.counts[key] = 0
            by.setdefault(chrom, []).append((start, end))
        for chrom, spans in by.items():
            spans.sort()
            self._by_chrom[chrom] = (
                [s for s, _ in spans],
                [(s, e, i) for i, (s, e) in enumerate(spans)],
            )
            self._reach[chrom] = max(e - s for s, e in spans)

    @property
    def denominator(self) -> int:
        """Reads passing the filter, countable only by `add`. There is deliberately NO setter.

        The registration forbids rescaling of any kind, and the obvious way to breach that is to put
        ENCODE's published mapped-read total here instead of the count from this pass. A ratio of 1.0
        between the two cannot be the evidence against it, because a ratio of exactly 1.0 is also what
        substituting the published total would produce -- the observation does not distinguish the
        compliant case from the breach. A property with no setter does: assigning to it raises
        AttributeError, so the substitution is not a thing this class can be made to do.
        """
        return self._denominator

    def __len__(self) -> int:
        return len(self.counts)

    def add(self, read: Read) -> int:
        """Offer one read. Returns how many elements it was counted in (0 if it was not kept).

        A kept read increments the denominator exactly once however many elements it overlaps, because
        the denominator is reads passing the filter and not reads-times-elements.
        """
        self.reads_seen += 1
        if not keep(read):
            self.reads_rejected += 1
            return 0
        self._denominator += 1
        chrom = read.reference_name
        if chrom is None or chrom not in self._by_chrom:
            if chrom is not None:
                self.reads_on_unknown_chrom += 1
            return 0
        start = read.reference_start
        end = read.reference_end
        if end is None or end <= start:
            return 0
        return self._count_span(chrom, start, end)

    def _count_span(self, chrom: str, start: int, end: int) -> int:
        starts, spans = self._by_chrom[chrom]
        reach = self._reach[chrom]
        hit = 0
        j = bisect.bisect_left(starts, end) - 1
        while j >= 0 and spans[j][0] > start - reach - 1:
            s, e, _ = spans[j]
            if e > start and s < end:
                self.counts[(chrom, s, e)] += 1
                hit += 1
            j -= 1
        return hit

    def rpm(self) -> dict[tuple[str, int, int], float]:
        """Reads per million per element. Every element gets a number; no element is dropped."""
        if self.denominator == 0:
            raise ValueError(
                "no read passed the filter, so reads per million has no denominator. A pass that "
                "kept nothing is a failed pass and is not reported as a column of zeros"
            )
        scale = 1e6 / self.denominator
        return {k: v * scale for k, v in self.counts.items()}

    def summary(self) -> dict[str, Any]:
        return {
            "elements": len(self.counts),
            "denominator_reads_passing_the_filter": self.denominator,
            "reads_seen": self.reads_seen,
            "reads_rejected_by_the_read_rule": self.reads_rejected,
            "reads_on_a_chromosome_with_no_element": self.reads_on_unknown_chrom,
            "elements_with_no_read": sum(1 for v in self.counts.values() if v == 0),
            "total_element_overlaps": sum(self.counts.values()),
            "read_rule": READ_RULE,
        }


def merge(counters: Iterable[Counter]) -> dict[str, Any]:
    """Pool several passes, which is what the registered replicate rule means by pooling.

    Numerators add and denominators add, which is arithmetically the same as having streamed the
    concatenation of the files. Pooling ratios instead would weight a shallow replicate equally with a
    deep one, so it is the counts that are pooled and the division happens once at the end.
    """
    counters = list(counters)
    if not counters:
        raise ValueError("nothing to pool")
    keys = {k for c in counters for k in c.counts}
    if any(set(c.counts) != keys for c in counters):
        raise ValueError("pooled passes must be over the same element set")
    num = {k: sum(c.counts[k] for c in counters) for k in keys}
    den = sum(c.denominator for c in counters)
    if den == 0:
        raise ValueError("pooled passes kept no read at all")
    scale = 1e6 / den
    return {
        "rpm": {k: num[k] * scale for k in keys},
        "numerator": num,
        "denominator": den,
        "passes": len(counters),
        "per_pass_denominator": [c.denominator for c in counters],
        "pooling_rule": "numerators and denominators are summed and divided once, which equals "
        "streaming the concatenation; ratios are never averaged",
    }


# ----------------------------------------------------------------------- comparing with published


def rank(xs: list[float]) -> list[float]:
    """Ranks with ties averaged, which is what a Spearman over a column holding many zeros needs."""
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    ranks = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        shared = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[order[k]] = shared
        i = j + 1
    return ranks


def spearman(xs: list[float], ys: list[float]) -> float | None:
    """Spearman rank correlation, ties averaged. None when it is undefined rather than 0.0."""
    if len(xs) != len(ys):
        raise ValueError(f"paired series must be the same length: {len(xs)} vs {len(ys)}")
    if len(xs) < 3:
        return None
    rx, ry = rank(xs), rank(ys)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry, strict=True))
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    if sxx <= 0 or syy <= 0:
        return None
    return sxy / math.sqrt(sxx * syy)


def median_ratio(computed: list[float], published: list[float]) -> dict[str, Any]:
    """The median of computed/published, over elements whose PUBLISHED value is strictly positive.

    A ratio to zero is undefined, so those elements cannot enter the median. The number excluded is
    reported, because an exclusion that is not counted is a filter nobody can see.
    """
    if len(computed) != len(published):
        raise ValueError("paired series must be the same length")
    ratios = [c / p for c, p in zip(computed, published, strict=True) if p > 0]
    excluded = len(published) - len(ratios)
    if not ratios:
        return {
            "median_ratio": None,
            "elements_used": 0,
            "elements_excluded_published_not_positive": excluded,
            "why": "no element has a positive published value, so no ratio is defined",
        }
    ratios.sort()
    n = len(ratios)
    med = ratios[n // 2] if n % 2 else (ratios[n // 2 - 1] + ratios[n // 2]) / 2
    return {
        "median_ratio": med,
        "elements_used": n,
        "elements_excluded_published_not_positive": excluded,
        "p05": ratios[max(0, int(0.05 * n) - 1)],
        "p95": ratios[min(n - 1, int(0.95 * n))],
    }


def gate(
    computed: list[float],
    published: list[float],
    spearman_min: float,
    ratio_range: tuple[float, float],
) -> dict[str, Any]:
    """The registered gate for one column: BOTH conditions, each reported whether or not it passes."""
    rho = spearman(computed, published)
    ratio = median_ratio(computed, published)
    med = ratio["median_ratio"]
    rho_ok = rho is not None and rho >= spearman_min
    ratio_ok = med is not None and ratio_range[0] <= med <= ratio_range[1]
    return {
        "spearman": rho,
        "spearman_min": spearman_min,
        "spearman_passes": rho_ok,
        **ratio,
        "median_ratio_range": list(ratio_range),
        "median_ratio_passes": ratio_ok,
        "passes": bool(rho_ok and ratio_ok),
        "rule": "both conditions must hold; either one failing fails the column",
        "elements_compared": len(computed),
    }


# ------------------------------------------------------- the producer's own rule, read from its code

#: Where every term below comes from. The rule is taken from the code that PRODUCED the columns, not
#: from a search for a rule that fits: the difference is not one of degree, and this records which was
#: done. `mayasheth/chrom-annotate` is the annotator whose config_CRISPR.yml declares
#: `RPM_assays: [CTCF, DHS, H3K27ac, H3K27me3, H3K4me1]` over the benchmark's own
#: `chrom/chromStart/chromEnd` columns for the validation set, which is the held-out five cell types.
PRODUCER = {
    "repository": "mayasheth/chrom-annotate",
    "commit": "91cda73ebe3a19153a582cab18cbf7ff70d85cfc",
    "file": "workflow/scripts/neighborhoods.py",
    "rpm_formula": "L579: df[feature_name + '.RPM'] = 1e6 * df[featurecount] / float(total_counts)",
    "numerator": "L433-445 count_bam: pysam.AlignmentFile(bam).count(chr, start, end). pysam's "
    "default read_callback='all' skips BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL and BAM_FDUP, and does "
    "NOT skip BAM_FSUPPLEMENTARY, so a duplicate is excluded from the numerator and a supplementary "
    "alignment is counted",
    "sex_chromosomes": "L425-431 double_sex_chrom_counts, called from run_count_reads: the count of "
    "any region whose contig name ENDS in X or Y is multiplied by 2, 'to make it seem like they have "
    "2 copies'",
    "denominator": "L662-673 count_bam_mapped: `samtools idxstats`, summing column 3 over every "
    "reference. That is mapped read-segments including duplicates, secondaries and supplementaries, "
    "so the denominator is NOT the numerator's filter applied to the whole file",
    "several_bams": "L596-602 average_features: df[feature + '.RPM'] = df[feature_RPM_cols].mean("
    "axis=1). The per-BAM RPMs are AVERAGED; the counts are not pooled",
    "files": "resources/metadata/epigenetic_datasets.tsv names them per biosample and assay. K562 "
    "DNase-seq ENCSR000EOT: ENCFF205FNC, ENCFF860XAE. K562 H3K27ac ChIP-seq ENCSR000AKP: ENCFF790GFL, "
    "ENCFF817HMW. Both pairs are ENCODE `unfiltered alignments`",
    "region": "config/config_CRISPR.yml: for the `validation` set the chr/start/end columns are the "
    "benchmark's own chrom, chromStart and chromEnd, so the region is the element span as published "
    "and is not resized. `RPM_expanded_assays` is a SEPARATE set of columns (.expandedRegion) and is "
    "not the column the frozen model reads",
}

#: The asymmetry is the producer's, and is implemented rather than corrected.
THE_ASYMMETRY_IS_THEIRS = (
    "the numerator excludes duplicates (pysam's count default) while the denominator counts them "
    "(samtools idxstats). That is what the code does, so that is what is implemented. Making the two "
    "agree would be a third rule of this lane's own invention, and the point of reading the producer's "
    "source was to stop inventing rules"
)


def keep_numerator_producer(read: Any) -> bool:
    """pysam `count()`'s default read_callback='all': drops unmapped, secondary, QC-fail, duplicate."""
    return not (
        read.is_unmapped
        or read.is_secondary
        or getattr(read, "is_qcfail", False)
        or getattr(read, "is_duplicate", False)
    )


def keep_denominator_producer(read: Any) -> bool:
    """`samtools idxstats` column 3: every mapped read-segment, duplicates and all."""
    return not read.is_unmapped


def doubles_on_sex_chromosome(chrom: str) -> bool:
    """The producer's own test: the contig NAME ends in X or Y."""
    return bool(chrom) and chrom[-1] in ("X", "Y")


class ProducerCounter(Counter):
    """`Counter` under the producer's rule: its two filters, and its sex-chromosome doubling.

    Kept as a subclass rather than a flag on `Counter` so the first rule's result stays readable as
    what it was. Two rules were run; both are on the record.
    """

    def add(self, read: Any) -> int:
        self.reads_seen += 1
        if keep_denominator_producer(read):
            self._denominator += 1
        else:
            self.reads_rejected += 1
        if not keep_numerator_producer(read):
            return 0
        chrom = read.reference_name
        if chrom is None or chrom not in self._by_chrom:
            if chrom is not None:
                self.reads_on_unknown_chrom += 1
            return 0
        start, end = read.reference_start, read.reference_end
        if end is None or end <= start:
            return 0
        return self._count_span(chrom, start, end)

    def rpm(self) -> dict[tuple[str, int, int], float]:
        """RPM with the producer's sex-chromosome doubling applied to the COUNT, as its code does."""
        if self.denominator == 0:
            raise ValueError("no mapped read, so reads per million has no denominator: a failed pass")
        scale = 1e6 / self.denominator
        return {k: (v * 2 if doubles_on_sex_chromosome(k[0]) else v) * scale for k, v in self.counts.items()}


def mean_of_per_bam_rpm(
    per_bam: list[dict[tuple[str, int, int], float]],
) -> dict[tuple[str, int, int], float]:
    """The producer's `average_features`: the MEAN of the per-BAM RPMs, element by element."""
    if not per_bam:
        raise ValueError("nothing to average")
    keys = set(per_bam[0])
    if any(set(d) != keys for d in per_bam):
        raise ValueError("averaged passes must be over the same element set")
    n = len(per_bam)
    return {k: sum(d[k] for d in per_bam) / n for k in keys}
