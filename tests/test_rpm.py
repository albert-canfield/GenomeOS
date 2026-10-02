"""The RPM reconstruction: the read rule, the overlap, the pooling and the registered gate.

Every test here runs on synthetic reads with no file and no network, because the counting rule is the
part that is easy to get wrong and expensive to re-run. The streaming is checked separately against a
real BAM; what is checked here is what a read means once it has arrived.
"""

from __future__ import annotations

import pytest

from genomeos.attribution import rpm


class FakeRead:
    """The minimum of pysam's AlignedSegment that the counter uses."""

    def __init__(
        self,
        chrom="chr1",
        start=100,
        end=136,
        unmapped=False,
        secondary=False,
        supplementary=False,
        paired=False,
        read1=True,
    ):
        self.reference_name = chrom
        self.reference_start = start
        self._end = end
        self.is_unmapped = unmapped
        self.is_secondary = secondary
        self.is_supplementary = supplementary
        self.is_paired = paired
        self.is_read1 = read1

    @property
    def reference_end(self):
        return self._end


E = ("chr1", 1_000, 1_500)


class TestTheReadRule:
    def test_a_plain_mapped_read_is_kept(self):
        assert rpm.keep(FakeRead()) is True

    @pytest.mark.parametrize("flag", ["unmapped", "secondary", "supplementary"])
    def test_unmapped_secondary_and_supplementary_are_rejected(self, flag):
        assert rpm.keep(FakeRead(**{flag: True})) is False

    def test_a_paired_fragment_is_counted_once_by_keeping_read1_only(self):
        assert rpm.keep(FakeRead(paired=True, read1=True)) is True
        assert rpm.keep(FakeRead(paired=True, read1=False)) is False

    def test_an_unpaired_read_is_kept_whatever_read1_says(self):
        """Single-ended data must not be silently halved by the fragment rule."""
        assert rpm.keep(FakeRead(paired=False, read1=False)) is True

    def test_a_duplicate_is_not_excluded_because_the_benchmark_states_no_such_rule(self):
        """A duplicate flag is deliberately not consulted; the registration says why."""
        r = FakeRead()
        r.is_duplicate = True
        assert rpm.keep(r) is True
        assert "does not state" in rpm.READ_RULE

    def test_no_mapq_threshold_is_applied_and_the_reason_is_recorded(self):
        r = FakeRead()
        r.mapping_quality = 0
        assert rpm.keep(r) is True
        assert "states no MAPQ threshold" in rpm.NO_THRESHOLD_OF_OUR_OWN


class TestTheOverlap:
    def counter(self):
        return rpm.Counter([E])

    @pytest.mark.parametrize(
        "start,end,hits",
        [
            (1_000, 1_036, 1),  # at the start
            (1_464, 1_500, 1),  # at the end
            (999, 1_001, 1),  # one base inside, from the left
            (1_499, 1_600, 1),  # one base inside, from the right
            (964, 1_000, 0),  # abuts the start, half-open, no overlap
            (1_500, 1_536, 0),  # abuts the end, half-open, no overlap
            (500, 2_000, 1),  # spans the whole element
            (1_100, 1_200, 1),  # inside it
        ],
    )
    def test_overlap_is_any_base_and_both_ends_are_half_open(self, start, end, hits):
        c = self.counter()
        assert c.add(FakeRead(start=start, end=end)) == hits
        assert c.counts[E] == hits

    def test_a_read_on_another_chromosome_is_counted_in_the_denominator_only(self):
        c = self.counter()
        assert c.add(FakeRead(chrom="chr2", start=1_100, end=1_136)) == 0
        assert c.denominator == 1, "it passed the filter, so it is part of the million"
        assert c.counts[E] == 0
        assert c.reads_on_unknown_chrom == 1

    def test_a_rejected_read_moves_neither_the_numerator_nor_the_denominator(self):
        c = self.counter()
        assert c.add(FakeRead(start=1_100, end=1_136, unmapped=True)) == 0
        assert c.denominator == 0
        assert c.counts[E] == 0
        assert c.reads_rejected == 1

    def test_a_read_overlapping_two_elements_counts_in_both_but_once_in_the_denominator(self):
        a, b = ("chr1", 1_000, 1_200), ("chr1", 1_150, 1_400)
        c = rpm.Counter([a, b])
        assert c.add(FakeRead(start=1_160, end=1_196)) == 2
        assert c.counts[a] == 1 and c.counts[b] == 1
        assert c.denominator == 1, "the denominator is reads, not reads times elements"

    def test_elements_much_longer_than_a_read_are_still_found(self):
        """The backward walk is bounded by the longest element, computed from the elements."""
        long = ("chr1", 0, 500_000)
        c = rpm.Counter([long, ("chr1", 900_000, 900_100)])
        assert c.add(FakeRead(start=499_000, end=499_036)) == 1
        assert c.counts[long] == 1

    def test_a_duplicate_element_is_one_element(self):
        c = rpm.Counter([E, E, ("chr1", 1_000, 1_500)])
        assert len(c) == 1
        c.add(FakeRead(start=1_100, end=1_136))
        assert c.counts[E] == 1, "an element given twice must not double its value"

    def test_an_empty_or_reversed_element_is_refused(self):
        with pytest.raises(ValueError, match="half-open"):
            rpm.Counter([("chr1", 500, 500)])
        with pytest.raises(ValueError, match="half-open"):
            rpm.Counter([("chr1", 600, 500)])

    def test_a_read_with_no_reference_end_is_not_counted_in_an_element(self):
        c = self.counter()
        r = FakeRead(start=1_100)
        r._end = None
        assert c.add(r) == 0
        assert c.denominator == 1


class TestReadsPerMillion:
    def test_rpm_is_a_million_times_the_share_of_reads(self):
        c = rpm.Counter([E])
        for _ in range(3):
            c.add(FakeRead(start=1_100, end=1_136))
        for _ in range(999_997):
            c.add(FakeRead(chrom="chrX", start=10, end=46))
        assert c.denominator == 1_000_000
        assert c.rpm()[E] == pytest.approx(3.0)

    def test_an_element_with_no_read_gets_a_measured_zero_and_is_not_dropped(self):
        c = rpm.Counter([E, ("chr1", 50_000, 50_500)])
        c.add(FakeRead(start=1_100, end=1_136))
        out = c.rpm()
        assert len(out) == 2
        assert out[("chr1", 50_000, 50_500)] == 0.0
        assert "MEASURED" in rpm.A_MEASURED_ZERO

    def test_a_pass_that_kept_nothing_is_a_failed_pass_and_not_a_column_of_zeros(self):
        c = rpm.Counter([E])
        c.add(FakeRead(unmapped=True))
        with pytest.raises(ValueError, match="failed pass"):
            c.rpm()

    def test_the_metadata_figure_is_never_the_denominator(self):
        assert "NEVER used as the denominator" in rpm.METADATA_IS_A_CROSS_CHECK


class TestPooling:
    def counters(self):
        a = rpm.Counter([E])
        b = rpm.Counter([E])
        a.add(FakeRead(start=1_100, end=1_136))
        for _ in range(9):
            a.add(FakeRead(chrom="chrX", start=1, end=37))
        for _ in range(3):
            b.add(FakeRead(start=1_100, end=1_136))
        for _ in range(97):
            b.add(FakeRead(chrom="chrX", start=1, end=37))
        return a, b

    def test_pooling_sums_counts_and_divides_once(self):
        a, b = self.counters()
        pooled = rpm.merge([a, b])
        assert pooled["numerator"][E] == 4
        assert pooled["denominator"] == 110
        assert pooled["rpm"][E] == pytest.approx(4 / 110 * 1e6)

    def test_pooling_is_not_an_average_of_the_two_rpms(self):
        """A shallow pass must not weigh as much as a deep one; that is why counts pool, not ratios."""
        a, b = self.counters()
        pooled = rpm.merge([a, b])["rpm"][E]
        averaged = (a.rpm()[E] + b.rpm()[E]) / 2
        assert pooled != pytest.approx(averaged)
        assert pooled < averaged, "the deeper pass should pull the pooled value down here"

    def test_pooling_equals_one_pass_over_the_concatenation(self):
        a, b = self.counters()
        both = rpm.Counter([E])
        for _ in range(4):
            both.add(FakeRead(start=1_100, end=1_136))
        for _ in range(106):
            both.add(FakeRead(chrom="chrX", start=1, end=37))
        assert rpm.merge([a, b])["rpm"][E] == pytest.approx(both.rpm()[E])

    def test_pooling_different_element_sets_is_refused(self):
        with pytest.raises(ValueError, match="same element set"):
            rpm.merge([rpm.Counter([E]), rpm.Counter([("chr2", 1, 2)])])

    def test_pooling_nothing_is_refused(self):
        with pytest.raises(ValueError, match="nothing to pool"):
            rpm.merge([])


class TestSpearmanAndTheRatio:
    def test_a_perfect_monotone_relation_is_one_even_when_not_linear(self):
        xs = [1.0, 2.0, 3.0, 4.0, 5.0]
        assert rpm.spearman(xs, [1.0, 4.0, 9.0, 16.0, 25.0]) == pytest.approx(1.0)

    def test_a_reversed_relation_is_minus_one(self):
        assert rpm.spearman([1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 1.0]) == pytest.approx(-1.0)

    def test_ties_are_averaged_so_a_column_of_many_zeros_is_handled(self):
        assert rpm.rank([0.0, 0.0, 0.0, 5.0]) == [2.0, 2.0, 2.0, 4.0]
        rho = rpm.spearman([0.0, 0.0, 1.0, 2.0], [0.0, 0.0, 3.0, 4.0])
        assert rho == pytest.approx(1.0)

    def test_a_constant_series_has_no_correlation_rather_than_zero(self):
        assert rpm.spearman([1.0, 1.0, 1.0, 1.0], [1.0, 2.0, 3.0, 4.0]) is None

    def test_too_few_points_is_none(self):
        assert rpm.spearman([1.0, 2.0], [1.0, 2.0]) is None

    def test_mismatched_lengths_are_refused(self):
        with pytest.raises(ValueError, match="same length"):
            rpm.spearman([1.0], [1.0, 2.0])

    def test_a_constant_factor_shows_up_in_the_median_ratio_not_the_correlation(self):
        """The two conditions catch different failures, which is why the gate needs both."""
        published = [1.0, 2.0, 3.0, 4.0, 5.0]
        computed = [2 * p for p in published]
        assert rpm.spearman(computed, published) == pytest.approx(1.0)
        assert rpm.median_ratio(computed, published)["median_ratio"] == pytest.approx(2.0)

    def test_elements_with_a_published_zero_are_excluded_from_the_ratio_and_counted(self):
        out = rpm.median_ratio([1.0, 2.0, 3.0], [0.0, 2.0, 3.0])
        assert out["elements_used"] == 2
        assert out["elements_excluded_published_not_positive"] == 1
        assert out["median_ratio"] == pytest.approx(1.0)

    def test_all_published_zero_gives_no_ratio_rather_than_a_wrong_one(self):
        out = rpm.median_ratio([1.0, 2.0], [0.0, 0.0])
        assert out["median_ratio"] is None
        assert out["elements_used"] == 0


class TestTheRegisteredGate:
    MIN, RANGE = 0.98, (0.9, 1.1)

    def test_a_faithful_reconstruction_passes_both_conditions(self):
        published = [float(i) for i in range(1, 51)]
        computed = [p * 1.01 for p in published]
        g = rpm.gate(computed, published, self.MIN, self.RANGE)
        assert g["passes"] is True
        assert g["spearman_passes"] and g["median_ratio_passes"]

    def test_a_constant_factor_out_fails_on_the_ratio_despite_a_perfect_correlation(self):
        published = [float(i) for i in range(1, 51)]
        g = rpm.gate([p * 3 for p in published], published, self.MIN, self.RANGE)
        assert g["spearman"] == pytest.approx(1.0)
        assert g["spearman_passes"] is True
        assert g["median_ratio_passes"] is False
        assert g["passes"] is False, "the ratio alone must be able to fail the column"

    def test_right_on_average_and_wrong_element_by_element_fails_on_the_correlation(self):
        published = [float(i) for i in range(1, 51)]
        computed = list(reversed(published))
        g = rpm.gate(computed, published, self.MIN, self.RANGE)
        assert g["median_ratio"] == pytest.approx(1.0, rel=0.1)
        assert g["median_ratio_passes"] is True
        assert g["spearman_passes"] is False
        assert g["passes"] is False, "the correlation alone must be able to fail the column"

    def test_the_gate_reports_both_numbers_whether_or_not_it_passes(self):
        published = [float(i) for i in range(1, 51)]
        g = rpm.gate([p * 3 for p in published], published, self.MIN, self.RANGE)
        for k in ("spearman", "median_ratio", "elements_compared", "spearman_min"):
            assert k in g
        assert g["elements_compared"] == 50

    def test_the_thresholds_are_taken_from_the_caller_and_not_hard_coded(self):
        published = [float(i) for i in range(1, 51)]
        computed = [p * 3 for p in published]
        assert rpm.gate(computed, published, 0.5, (2.0, 4.0))["passes"] is True


class TestOneBamPerReplicate:
    """The guard that stands between a pooled total and a silent threefold double count.

    Shown to FIRE on the real shape that caused it: ENCSR000AKP's seven released filtered GRCh38
    alignment BAMs, which are three biological replicates reprocessed by three pipeline versions.
    """

    #: The real seven, as the portal reports them.
    SEVEN = [
        {"accession": "ENCFF301TVL", "biological_replicates": [1]},  # ENCODE3
        {"accession": "ENCFF121RHF", "biological_replicates": [1]},  # ENCODE4 v1.5.1
        {"accession": "ENCFF600THN", "biological_replicates": [1]},  # ENCODE4 v1.8.0
        {"accession": "ENCFF879BWC", "biological_replicates": [2]},
        {"accession": "ENCFF907MNY", "biological_replicates": [2]},
        {"accession": "ENCFF704LGA", "biological_replicates": [2]},
        {"accession": "ENCFF232RQF", "biological_replicates": [3]},
    ]

    #: The released analysis alone, which is what the corrected rule selects.
    RELEASED = [
        {"accession": "ENCFF600THN", "biological_replicates": [1]},
        {"accession": "ENCFF704LGA", "biological_replicates": [2]},
        {"accession": "ENCFF232RQF", "biological_replicates": [3]},
    ]

    def test_the_real_seven_bams_are_refused(self):
        with pytest.raises(ValueError, match="appears in both"):
            rpm.check_one_bam_per_replicate(self.SEVEN)

    def test_the_released_analysis_alone_is_accepted_and_maps_each_replicate(self):
        assert rpm.check_one_bam_per_replicate(self.RELEASED) == {
            1: "ENCFF600THN",
            2: "ENCFF704LGA",
            3: "ENCFF232RQF",
        }

    def test_the_refusal_names_both_files_so_the_duplicate_can_be_found(self):
        with pytest.raises(ValueError) as exc:
            rpm.check_one_bam_per_replicate(self.SEVEN)
        msg = str(exc.value)
        assert "ENCFF121RHF" in msg and "ENCFF301TVL" in msg

    def test_a_bam_naming_no_replicate_is_refused_rather_than_assumed_unique(self):
        with pytest.raises(ValueError, match="names no biological replicate"):
            rpm.check_one_bam_per_replicate([{"accession": "ENCFFxxxxxx", "biological_replicates": []}])

    def test_a_bam_pooling_two_replicates_still_conflicts_with_either_of_them(self):
        """A pooled BAM plus one of its own replicates is the same double count, less obviously."""
        with pytest.raises(ValueError, match="appears in both"):
            rpm.check_one_bam_per_replicate(
                [
                    {"accession": "ENCFFaaaaaa", "biological_replicates": [1, 2]},
                    {"accession": "ENCFFbbbbbb", "biological_replicates": [2]},
                ]
            )

    def test_an_empty_set_is_refused(self):
        with pytest.raises(ValueError, match="no BAM to pool"):
            rpm.check_one_bam_per_replicate([])

    def test_the_reason_is_recorded_in_the_module(self):
        assert "does not mean" in rpm.ONE_BAM_PER_REPLICATE


class ProducerRead(FakeRead):
    """A read with the two flags the producer's filters consult and the first rule ignored."""

    def __init__(self, duplicate=False, qcfail=False, **kw):
        super().__init__(**kw)
        self.is_duplicate = duplicate
        self.is_qcfail = qcfail


class TestTheProducersOwnRule:
    """The rule read out of mayasheth/chrom-annotate, not searched for.

    Each test pins one term to the producer's code so a later reader can check the implementation
    against the citation rather than against an intention.
    """

    E = ("chr1", 1_000, 1_500)
    X = ("chrX", 1_000, 1_500)

    def test_the_citation_names_repository_commit_and_file(self):
        assert rpm.PRODUCER["repository"] == "mayasheth/chrom-annotate"
        assert len(rpm.PRODUCER["commit"]) == 40
        assert rpm.PRODUCER["file"].endswith("neighborhoods.py")
        for term in ("rpm_formula", "numerator", "denominator", "sex_chromosomes", "several_bams"):
            assert "L" in rpm.PRODUCER[term], f"{term} must cite a line"

    def test_a_duplicate_leaves_the_numerator_but_stays_in_the_denominator(self):
        """The producer's asymmetry: pysam count drops duplicates, idxstats counts them."""
        c = rpm.ProducerCounter([self.E])
        assert c.add(ProducerRead(start=1_100, end=1_136, duplicate=True)) == 0
        assert c.counts[self.E] == 0, "a duplicate is not in the numerator"
        assert c.denominator == 1, "a duplicate IS in the denominator"
        assert "what the code does, so that is what is implemented" in rpm.THE_ASYMMETRY_IS_THEIRS

    def test_a_qc_fail_read_is_dropped_from_the_numerator_only(self):
        c = rpm.ProducerCounter([self.E])
        assert c.add(ProducerRead(start=1_100, end=1_136, qcfail=True)) == 0
        assert c.counts[self.E] == 0
        assert c.denominator == 1

    def test_a_supplementary_read_IS_counted_because_pysam_does_not_drop_it(self):
        """BAM_FSUPPLEMENTARY is absent from read_callback='all', so it counts. Easy to get wrong."""
        c = rpm.ProducerCounter([self.E])
        assert c.add(ProducerRead(start=1_100, end=1_136, supplementary=True)) == 1
        assert c.counts[self.E] == 1

    def test_a_secondary_read_is_dropped_from_the_numerator(self):
        c = rpm.ProducerCounter([self.E])
        assert c.add(ProducerRead(start=1_100, end=1_136, secondary=True)) == 0
        assert c.denominator == 1

    def test_an_unmapped_read_is_in_neither(self):
        c = rpm.ProducerCounter([self.E])
        assert c.add(ProducerRead(start=1_100, end=1_136, unmapped=True)) == 0
        assert c.denominator == 0

    def test_read2_of_a_pair_is_counted_unlike_the_first_rule(self):
        """The producer does not count fragments once; it counts read-segments."""
        c = rpm.ProducerCounter([self.E])
        assert c.add(ProducerRead(start=1_100, end=1_136, paired=True, read1=False)) == 1

    def test_counts_on_chrX_and_chrY_are_doubled(self):
        for chrom in ("chrX", "chrY"):
            el = (chrom, 1_000, 1_500)
            c = rpm.ProducerCounter([el])
            c.add(ProducerRead(chrom=chrom, start=1_100, end=1_136))
            for _ in range(999_999):
                c.add(ProducerRead(chrom="chr1", start=10, end=46))
            assert c.rpm()[el] == pytest.approx(2.0), f"{chrom} must be doubled"

    def test_an_autosome_is_not_doubled(self):
        c = rpm.ProducerCounter([self.E])
        c.add(ProducerRead(start=1_100, end=1_136))
        for _ in range(999_999):
            c.add(ProducerRead(chrom="chr2", start=10, end=46))
        assert c.rpm()[self.E] == pytest.approx(1.0)

    def test_the_doubling_test_is_on_the_name_ending_as_the_producer_wrote_it(self):
        assert rpm.doubles_on_sex_chromosome("chrX") is True
        assert rpm.doubles_on_sex_chromosome("chrY") is True
        assert rpm.doubles_on_sex_chromosome("chr1") is False
        assert rpm.doubles_on_sex_chromosome("") is False

    def test_several_bams_are_averaged_and_not_pooled(self):
        a = {self.E: 10.0}
        b = {self.E: 20.0}
        assert rpm.mean_of_per_bam_rpm([a, b])[self.E] == pytest.approx(15.0)

    def test_averaging_differs_from_pooling_when_depths_differ(self):
        """The two combination rules give different answers, which is why the citation matters."""
        shallow = rpm.ProducerCounter([self.E])
        deep = rpm.ProducerCounter([self.E])
        shallow.add(ProducerRead(start=1_100, end=1_136))
        for _ in range(9):
            shallow.add(ProducerRead(chrom="chr2", start=1, end=37))
        for _ in range(3):
            deep.add(ProducerRead(start=1_100, end=1_136))
        for _ in range(997):
            deep.add(ProducerRead(chrom="chr2", start=1, end=37))
        averaged = rpm.mean_of_per_bam_rpm([shallow.rpm(), deep.rpm()])[self.E]
        pooled = rpm.merge([shallow, deep])["rpm"][self.E]
        assert averaged != pytest.approx(pooled)
        assert averaged > pooled, "averaging lets the shallow pass carry equal weight"

    def test_averaging_different_element_sets_is_refused(self):
        with pytest.raises(ValueError, match="same element set"):
            rpm.mean_of_per_bam_rpm([{self.E: 1.0}, {self.X: 1.0}])

    def test_averaging_nothing_is_refused(self):
        with pytest.raises(ValueError, match="nothing to average"):
            rpm.mean_of_per_bam_rpm([])


class TestNoRescalingOfAnyKind:
    """The registered term, enforced structurally rather than observed.

    A ratio of 1.0 between the computed denominator and ENCODE's published total cannot be the evidence
    that the published total was not used, because substituting it would also produce exactly 1.0. These
    tests make the substitution impossible instead.
    """

    E = ("chr1", 1_000, 1_500)

    def test_the_registered_term_is_carried_word_for_word(self):
        assert "NO RESCALING of any kind" in rpm.NO_RESCALING_OF_ANY_KIND

    def test_the_published_total_cannot_be_assigned_to_the_denominator(self):
        """By planting the breach: assigning ENCODE's figure raises instead of taking effect."""
        c = rpm.Counter([self.E])
        c.add(FakeRead(start=1_100, end=1_136))
        with pytest.raises(AttributeError):
            c.denominator = 219_538_839
        assert c.denominator == 1, "the count from this pass is untouched"

    def test_the_producer_counter_denominator_is_equally_unassignable(self):
        c = rpm.ProducerCounter([self.E])
        c.add(ProducerRead(start=1_100, end=1_136))
        with pytest.raises(AttributeError):
            c.denominator = 10_925_509
        assert c.denominator == 1

    def test_only_add_can_move_the_denominator(self):
        c = rpm.Counter([self.E])
        assert c.denominator == 0
        c.add(FakeRead(chrom="chr9", start=1, end=37))
        assert c.denominator == 1
        c.add(FakeRead(unmapped=True))
        assert c.denominator == 1, "a rejected read does not move it either"

    def test_rpm_is_computed_from_the_counted_denominator_and_nothing_else(self):
        """The arithmetic is pinned, so any substituted denominator would change the number."""
        c = rpm.Counter([self.E])
        for _ in range(7):
            c.add(FakeRead(start=1_100, end=1_136))
        for _ in range(93):
            c.add(FakeRead(chrom="chr9", start=1, end=37))
        assert c.denominator == 100
        assert c.rpm()[self.E] == pytest.approx(7 / 100 * 1e6)

    def test_the_rpm_module_never_fetches_the_published_total(self):
        """The metadata lookup lives in the runner; the counting module has no route to it."""
        import inspect

        src = inspect.getsource(rpm)
        assert "metadata_mapped_reads" not in src
        assert "urllib" not in src and "requests" not in src
