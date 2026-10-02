"""Amendment 1 to the AstroREG registration: the calibration gate for RECONSTRUCTING the frozen
activity columns from ENCODE alignments, written before a single read is counted.

The original registration (data/results/astroreg_registration.json, sha256
43f2edfab0636e0f07a035aedf5ddf7c57aa1dfcee8d39ea32e9338edd522c7a) stays frozen and is cited by that
hash. This is ADDITIVE: it adds a way to obtain an input the original needs and could not supply, and
it is written before any astrocyte deletion value, any reconstructed RPM and any outcome effect is
read. That ordering is the whole reason it is an amendment and not a term moved after a result.

WHAT THE NO-GO ESTABLISHED, and what this changes. The registered claim's activity term is DHS.RPM
and H3K27ac.RPM, reads per million mapped reads in the element. No astrocyte pair carries either, and
the only astrocyte activity on disk is a narrowPeak signalValue -- wrong units, and censored outside
called peaks. A peak call is a SUBSTITUTE and the registration forbids substitutes. An RPM computed
from the alignments is not a substitute: it is the quantity itself, by its own definition. That is the
distinction this amendment rests on, and it is the only thing that makes the path admissible.

THE GATE IS WRITTEN HERE, BEFORE THE FIRST COUNT. Reconstruct the columns from the K562 alignments and
compare with the PUBLISHED DHS.RPM and H3K27ac.RPM on the K562 held-out elements. Pass only if, for
EACH column, Spearman >= 0.98 AND the median ratio of computed to published lies in [0.9, 1.1]. Both
numbers are in this file before anything is computed and neither moves afterwards for any reason. If
it fails, the no-go STANDS, nothing is fetched for astrocytes, and the finding is about the two
published columns themselves -- that they are not recoverable from the alignments -- which bears on
every result in this project that rests on them. That outcome is the more informative of the two.

WHAT THE BENCHMARK DOES NOT DOCUMENT, established and recorded rather than filled in silently. Its
README defines the columns ("DNase-seq RPM in perturbed element", "H3K27ac ChIP-seq RPM in perturbed
element") and cites Gschwind et al. 2025. It names NO ENCODE accession, NO MAPQ threshold, NO
duplicate rule and NO replicate rule, and neither does the upstream repository. So this file states the
rule it uses and why, and the same rule is applied to K562 and to astrocyte without exception. A rule
chosen here and applied to both is auditable; a rule guessed for one of them is not.
"""

from __future__ import annotations

import gzip
import json
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri  # noqa: E402
from genomeos.attribution.measured import CRISPRI_SPLIT_OF  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "astroreg_calibration_registration"
ENTRY = "scripts/astroreg_calib_register.py"
OWN_CODE = (
    "genomeos/attribution/astrorun.py",
    "scripts/astroreg_calib_register.py",
    "tests/test_astrorun.py",
)

ORIGINAL = Path("data/results/astroreg_registration.json")
ORIGINAL_SHA256 = "43f2edfab0636e0f07a035aedf5ddf7c57aa1dfcee8d39ea32e9338edd522c7a"

BENCHMARK_README = crispri.KNOWLEDGE / "README.md"

# --------------------------------------------------------------------------- the gate, fixed here

#: The pass condition, per column. Written before any count and not moved afterwards.
SPEARMAN_MIN = 0.98
MEDIAN_RATIO_RANGE = (0.9, 1.1)

TOLERANCE_RULE = (
    "per column, BOTH must hold: Spearman rank correlation of computed against published >= 0.98, "
    "AND the median of (computed / published) in [0.9, 1.1]. The correlation alone would pass a "
    "reconstruction that was a constant factor out, and the ratio alone would pass one that was "
    "right on average and wrong element by element. The pair of them is the condition; either one "
    "failing fails the column, and either column failing fails the gate"
)

RATIO_DENOMINATOR_RULE = (
    "the median ratio is taken over elements whose PUBLISHED value is strictly positive, because a "
    "ratio to zero is undefined. The count of elements excluded for that reason is reported beside "
    "the median, and the Spearman is computed over ALL elements including those, so the exclusion "
    "cannot quietly become a filter on the comparison"
)

# ------------------------------------------------------------------- the read rule, fixed here

READ_RULE = {
    "alignments": "the filtered GRCh38 `alignments` BAMs ENCODE publishes, used AS THEY ARE",
    "why_as_they_are": "the benchmark states no MAPQ threshold and no duplicate rule, so any "
    "threshold applied here would be this lane's invention presented as the published definition. "
    "ENCODE's filtered alignments already apply that experiment's own pipeline thresholds, and "
    "taking them as given is the one choice that adds nothing of ours. The unfiltered alignments are "
    "NOT used",
    "reads_counted": "a read is counted when it is not unmapped, not secondary and not supplementary. "
    "Nothing else is excluded here: duplicates are whatever the filtered BAM delivers, because "
    "removing them would be a rule the benchmark does not state",
    "fragments": "a paired-end fragment is counted ONCE. Where reads are paired, only read1 is "
    "counted, so a fragment overlapping an element contributes 1 and not 2; single-ended data is "
    "unaffected, and which case each experiment is is recorded with the result",
    "overlap": "any base. A read overlaps an element when its reference span intersects the element's "
    "span by at least one base, half-open on both sides",
    "denominator": "counted IN THE SAME PASS, under the SAME filters, as the number of reads passing "
    "them. The ENCODE metadata figure is a cross-check printed beside it and is NEVER the "
    "denominator: a numerator from this pass divided by a denominator from someone else's filters "
    "would be wrong in a way no test of either number alone could catch",
    "rpm": "reads per million = 1e6 * reads overlapping the element / reads passing the filters",
    "nothing_on_disk": "one sequential HTTPS streaming pass per BAM, BGZF decoded as it arrives, "
    "memory bounded by the element list. No index is required because there is no random access, and "
    "no BAM is written to disk at any point",
    "size_is_not_a_criterion": "no file is chosen or rejected for its size. Streaming makes size a "
    "cost in time and not in disk",
}

REPLICATE_RULE = (
    "the benchmark documents no replicate rule, so one is stated here and applied identically to "
    "K562 and to astrocyte: ALL of an experiment's released filtered GRCh38 `alignments` BAMs are "
    "streamed and POOLED -- every one of them contributes to both the numerator and the denominator, "
    "which is the same thing as concatenating them. Pooling is chosen over picking one replicate "
    "because picking one would need a reason to prefer it that the benchmark does not give. The "
    "number of BAMs pooled per experiment is recorded with the result"
)

# --------------------------------------------------------------- which files, decided by measurement

H3K27AC_K562 = {
    "experiment": "ENCSR000AKP",
    "why": "the only released K562 H3K27ac ChIP-seq experiment on GRCh38 in the ENCODE portal, so "
    "there is no choice to make and none is made",
    "selection": "forced",
}

DNASE_K562_SELECTION = {
    "candidates": "the K562 DNase-seq experiments on GRCh38 whose files ENCODE itself flags "
    "preferred_default (14 experiments at the time of writing, among 62 released K562 DNase-seq "
    "experiments). ENCODE's own flag narrows the field; this lane does not",
    "problem": "the benchmark names no accession, so which of the 14 produced the published DHS.RPM "
    "is unknown. Picking one would make a gate failure uninterpretable: it could mean the "
    "reconstruction is wrong, or only that the wrong source file was read, and those have opposite "
    "consequences for every result that uses the column",
    "rule": "the source experiment is IDENTIFIED BY MEASUREMENT before any BAM is streamed. For each "
    "candidate, read its `read-depth normalized signal` bigWig over the K562 held-out elements by "
    "RANGED remote reads, and take the Spearman of that signal against the published DHS.RPM. The "
    "candidate with the highest Spearman is the selected experiment, and every candidate's figure is "
    "recorded so the margin between first and second is visible",
    "what_the_bigwig_does_and_does_not_do": "it selects a FILE and supplies NO VALUE. No RPM in this "
    "amendment is derived from a bigWig; the reconstructed column comes from the alignments alone. "
    "This is a third use of the track, distinct from the diagnostic arm below",
    "if_no_candidate_is_clear": "if the best Spearman does not exceed the second best by at least "
    "0.02, the source is recorded as UNIDENTIFIED and the DNase arm is reported as such rather than "
    "run against a guess. That threshold is fixed here, before any correlation is computed",
    "margin_required": 0.02,
}

BIGWIG_DIAGNOSTIC = (
    "the DNase `read-depth normalized signal` bigWig is also run BESIDE the BAM pass on K562 as a "
    "DIAGNOSTIC ONLY, through the same tolerance for information. The BAM pass is PRIMARY because it "
    "is the definition of reads per million. If both pass, the BAM is used. If ONLY the bigWig "
    "passes, this lane STOPS and reports without choosing: that would say something real about what "
    "'read-depth normalized signal' is normalised by, and it is not this lane's to decide"
)

NOT_CLAIMED = [
    "not a demonstration that the published columns were produced this way: a reconstruction that "
    "agrees to this tolerance is consistent with the published method and does not establish it",
    "not an identification of the benchmark's source files as a documented fact: the DNase experiment "
    "is identified by correlation, which is evidence about which file matches and not a record of "
    "which file was used",
    "not a licence to apply any other substitute: this amendment admits ONE thing, an RPM computed "
    "from alignments by its own definition, and the original registration's bar on substitutes stands "
    "for everything else",
    "no claim about astrocytes at all unless the gate passes, and none about the deletion gain until "
    "the registered test is run under the original terms",
]


def k562_heldout_elements() -> dict[str, Any]:
    """The comparison set: the K562 held-out elements and their published activity values.

    Counted here so the registration states the size of the comparison before it is made. The values
    themselves are read because they are the thing being compared against and they are already-read
    development evidence; no astrocyte value is touched.
    """
    path = crispri.KNOWLEDGE / crispri.HELDOUT
    elements: dict[tuple[str, int, int], dict[str, float | None]] = {}
    pairs = 0
    disagreements = 0
    for row in _rows(path):
        if row["CellType"] != "K562":
            continue
        pairs += 1
        key = (row["chrom"], int(row["chromStart"]), int(row["chromEnd"]))
        vals = {c: _num(row.get(c)) for c in ("DHS.RPM", "H3K27ac.RPM")}
        if key in elements and elements[key] != vals:
            disagreements += 1
        elements[key] = vals
    both = [k for k, v in elements.items() if all(x is not None for x in v.values())]
    positive = {c: sum(1 for k in elements if (elements[k][c] or 0) > 0) for c in ("DHS.RPM", "H3K27ac.RPM")}
    return {
        "source": str(path),
        "partition": CRISPRI_SPLIT_OF[crispri.HELDOUT],
        "cell_type": "K562",
        "pairs": pairs,
        "elements": len(elements),
        "elements_with_both_columns": len(both),
        "elements_with_a_positive_published_value": positive,
        "elements_whose_rows_disagreed_on_a_value": disagreements,
        "note": "the comparison is per ELEMENT, not per pair: the activity columns are a property of "
        "the element, and several pairs share one element. A disagreement between two rows of the "
        "same element would mean the columns are not an element property after all, so the count is "
        "reported and is expected to be nil",
    }


def _rows(path: Path):
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(header, line.rstrip("\n").split("\t"), strict=False))


def _num(raw: str | None) -> float | None:
    raw = (raw or "").strip()
    if raw in ("", "NA", "nan"):
        return None
    return float(raw)


def benchmark_documentation() -> dict[str, Any]:
    """What the benchmark's own documentation says about the two columns, read from the file."""
    text = BENCHMARK_README.read_text()
    lines = [ln.strip() for ln in text.splitlines() if ".RPM" in ln and "percentile" not in ln]
    import re

    return {
        "readme": str(BENCHMARK_README),
        "cites": "Gschwind et al., 2025",
        "column_definitions_found": [ln.lstrip("- ") for ln in lines][:6],
        "encode_accessions_named": bool(re.search(r"ENCFF\w+|ENCSR\w+", text)),
        "mapq_mentioned": bool(re.search(r"\bMAPQ\b|mapping quality", text, re.I)),
        "duplicates_mentioned": bool(re.search(r"duplicat", text, re.I)),
        "replicates_mentioned": bool(re.search(r"replicat", text, re.I)),
        "upstream_repository": "EngreitzLab/CRISPR_comparison, main: 29 files, none naming an ENCODE "
        "accession for either column (config/config.yml, resources/crispr_data/README.md and the "
        "genomic_features and genome_annotations directories were the only candidates)",
        "finding": "the documentation DEFINES the columns and does not say how they were made. It "
        "names no ENCODE accession, no MAPQ threshold, no duplicate rule and no replicate rule. That "
        "is recorded here as the reason this amendment states its own rule, rather than being filled "
        "in from what a reader might assume",
    }


def manifest(comparison: dict[str, Any]) -> dict[str, Any]:
    return {
        "sources": [
            {
                "accession": str(ORIGINAL),
                "version": f"the frozen original registration, sha256 {ORIGINAL_SHA256}; amended "
                "additively and not edited",
            },
            {
                "accession": "EngreitzLab/CRISPR_comparison, EPCrisprBenchmark heldout_5_cell_types "
                "(Gschwind et al. 2025)",
                "version": "fetched 2026-09-16; pinned here by sha256",
                "url": crispri.BASE_URL,
            },
            {
                "accession": "ENCODE portal metadata for the candidate K562 and astrocyte DNase-seq "
                "and H3K27ac ChIP-seq experiments",
                "version": "queried 2026-10-02; no file was downloaded and no alignment was read by "
                "this script",
            },
        ],
        "inputs": [
            mf.input_entry(
                ORIGINAL,
                partition=None,
                role="the registration this amends: its hash is recorded and no term of it is moved",
            ),
            mf.input_entry(
                crispri.KNOWLEDGE / crispri.HELDOUT,
                partition=CRISPRI_SPLIT_OF[crispri.HELDOUT],
                role="the K562 held-out elements and their published DHS.RPM and H3K27ac.RPM, which "
                "are the comparison the gate is against; already-read development evidence",
            ),
            mf.input_entry(
                BENCHMARK_README,
                partition=None,
                role="the benchmark's own column documentation, read to establish what it does NOT "
                "state about accessions, MAPQ, duplicates and replicates",
            ),
        ],
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "spearman_min": SPEARMAN_MIN,
            "median_ratio_min": MEDIAN_RATIO_RANGE[0],
            "median_ratio_max": MEDIAN_RATIO_RANGE[1],
            "dnase_selection_margin_required": DNASE_K562_SELECTION["margin_required"],
            "k562_heldout_elements": comparison["elements"],
            "alphagenome_requests": 0,
            "model_requests": 0,
            "requests_sent": 0,
            "money_spent": 0,
            "bams_streamed": 0,
            "bytes_written_to_disk": 0,
        },
        "exclusions": [
            "nothing is computed here: no read is counted, no BAM is opened, no bigWig is read and no "
            "correlation is taken. This file is the design, written first",
            "no AlphaGenome request is sent and this script has no code path that could send one",
            "no astrocyte value of any kind is read",
        ],
        "partitions": {
            "k562_heldout": "the comparison set for the gate, already-read development evidence",
            "astrocyte": "not touched by this script, and not touched at all unless the gate passes",
        },
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    t0 = time.time()
    digest, _, _ = mf.sha256_of(ORIGINAL)
    if digest != ORIGINAL_SHA256:
        print(f"refused: {ORIGINAL} is sha256 {digest}, not the frozen {ORIGINAL_SHA256}")
        return 2

    print("reading the K562 held-out comparison set (published values only)")
    comparison = k562_heldout_elements()
    print(f"  {comparison['elements']} elements over {comparison['pairs']} pairs")
    docs = benchmark_documentation()

    payload = {
        "status": "a registration, written before anything is computed. No read counted, no BAM "
        "opened, no bigWig read, no correlation taken, no request sent, no money spent",
        "lane": "lane-astro",
        "amends": {
            "registration": str(ORIGINAL),
            "sha256": ORIGINAL_SHA256,
            "kind": "additive amendment 1",
            "what_it_adds": "a way to obtain DHS.RPM and H3K27ac.RPM for a cell type the benchmark "
            "does not cover, by computing them from the ENCODE alignments rather than substituting a "
            "different quantity for them",
            "what_it_does_not_change": "no threshold, label, floor, interval rule, reading or power "
            "class of the original moves. The primary claim, the label rule, the two bootstraps, the "
            "three readings and the not_claimed list stand exactly as frozen",
            "why_it_is_an_amendment_and_not_a_term_moved": "it is written before any reconstructed "
            "RPM, any astrocyte deletion value and any outcome effect has been computed or read",
            "approval": "the original authorisation named a specific hash, so the inputs behind the "
            "claim changing is a matter for Albert. This lane does not treat this file as approval to "
            "spend any of the 1,322, and spends none",
        },
        "the_distinction_this_rests_on": (
            "a narrowPeak signalValue is a SUBSTITUTE for reads per million: different units, and "
            "censored outside called peaks, so an element in no peak could only be given an invented "
            "zero. An RPM computed from the alignments is not a substitute, it is the quantity's own "
            "definition, and an element with no reads gets a MEASURED zero. The original "
            "registration forbids the first and says nothing against the second"
        ),
        "benchmark_documentation": docs,
        "gate": {
            "scope": "the K562 held-out elements, per column, against the published values",
            "spearman_min": SPEARMAN_MIN,
            "median_ratio_range": list(MEDIAN_RATIO_RANGE),
            "rule": TOLERANCE_RULE,
            "ratio_denominator_rule": RATIO_DENOMINATOR_RULE,
            "written_before_any_count": True,
            "on_failure": "the no-go STANDS, nothing is fetched for astrocytes, and the finding is "
            "about the two published columns themselves: that they are not recoverable from the "
            "alignments. That bears on every result in this project that rests on them, so it is "
            "reported as a finding and not as a dead end",
        },
        "read_rule": READ_RULE,
        "replicate_rule": REPLICATE_RULE,
        "files": {
            "h3k27ac_k562": H3K27AC_K562,
            "dnase_k562_selection": DNASE_K562_SELECTION,
            "astrocyte": "the experiments the original registration already names, ENCSR000AOQ for "
            "H3K27ac and ENCSR000EPM for DNase, read through the IDENTICAL read and replicate rules "
            "as K562. Neither is touched unless the gate passes",
            "no_index_exists": "neither astrocyte experiment publishes a BAM index: a portal search "
            "for file_format=bai over both returns nothing, the conventional .bam.bai URL is 404, and "
            "the file records name no index. The objects ARE range-readable (HTTP 206 from "
            "encode-public.s3.amazonaws.com, BGZF magic confirmed), which is why the mechanism is one "
            "sequential streaming pass and not an indexed region fetch",
        },
        "bigwig_diagnostic": BIGWIG_DIAGNOSTIC,
        "comparison_set": comparison,
        "metadata_cross_check": {
            "rule": "the ENCODE metadata mapped-read figure is printed beside the computed denominator "
            "with their ratio, and is never used as the denominator",
            "astrocyte_h3k27ac_mapped_reads_ENCFF312ZIP": 10_925_509,
            "astrocyte_dnase_mapped_ENCFF586NXB": 229_650_233,
            "note": "these two were read from the portal's quality metrics and are recorded as the "
            "cross-check values, not as inputs to any RPM",
        },
        "not_claimed": NOT_CLAIMED,
        "requests_sent": 0,
        "money_spent": 0,
        "left_undone": [
            "nothing is computed by this file: the streaming pass, the selection correlation and the "
            "gate are all still to run",
            "which ENCODE files produced the published columns remains undocumented upstream; the "
            "DNase selection identifies a best match and does not establish provenance",
            "the original registration's stale left_undone entry about the astrocyte accessions is "
            "reported in data/results/astroreg_activity_nogo.json and is still not edited, because "
            "the artefact is frozen",
        ],
    }
    payload[mf.KEY] = manifest(comparison)
    path = save_result(RESULT, payload)

    print(
        json.dumps(
            {
                "spearman_min": SPEARMAN_MIN,
                "median_ratio_range": list(MEDIAN_RATIO_RANGE),
                "k562_heldout_elements": comparison["elements"],
                "elements_with_both_columns": comparison["elements_with_both_columns"],
                "rows_disagreeing_on_an_element_value": comparison[
                    "elements_whose_rows_disagreed_on_a_value"
                ],
                "benchmark_names_encode_accessions": docs["encode_accessions_named"],
                "benchmark_states_a_mapq_threshold": docs["mapq_mentioned"],
                "requests_sent": 0,
                "money_spent": 0,
            },
            indent=1,
        )
    )
    print(f"({time.time() - t0:.0f} s) -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
