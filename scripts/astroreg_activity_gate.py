"""Whether the frozen model can be applied to the astrocyte pairs at all, decided before any money
is spent.

The authorised run buys one thing: an astrocyte deletion value per registry element. The registered
primary claim needs two arms, 'activity + distance + deletion' and 'activity + distance', and
`activity` is in both. This script asks where the astrocyte side's `activity` would come from, and
answers it from the files rather than from the registration's intent.

What it establishes, each from disk:

  (a) the frozen activity term is `DHS.RPM` and `H3K27ac.RPM`, read in exactly one place
      (`crispri.parse`), and the EPCrisprBenchmark publishes both for every pair of all five of its
      held-out cell types -- which is why the HCT116 arm could be scored on the frozen path;
  (b) the astrocyte screen is not in that benchmark and its Supplementary Table 3 has no activity
      sheet, so no astrocyte pair carries either column;
  (c) the only astrocyte activity on disk is the registration's own two inputs, and both are
      narrowPeak PEAK CALLS whose third column is signalValue, not reads per million;
  (c2) an uncensored continuous astrocyte H3K27ac track does exist in the signal store, so the
      censoring argument alone would be answerable; it is counted, and the verdict is unchanged
      because that track is fold change over control and the store holds no DNase profile for any
      biosample, leaving DHS.RPM with no continuous source at all;
  (d) the peak calls are also censored: an element in no called peak can only be given zero, and
      this counts how many of the screen's elements that is;
  (e) the two ENCODE accessions ARE recorded in the repository, in each per-chromosome file's own
      first line, so the registration's `left_undone` entry saying they are not is the stale one.

By the registration's own feasibility rule -- "if an input the frozen model needs is missing, the
model is not applied and no substitute is put in its place: a no-go is recorded instead" -- the
verdict is a no-go, and the number of requests it permits is zero. No request is sent by this
script and it has no code path that could send one.
"""

from __future__ import annotations

import gzip
import json
import statistics
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from astroreg_register import TABLE3, screen_pairs  # noqa: E402

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import astroreg, astrorun, crispri  # noqa: E402
from genomeos.attribution.measured import CRISPRI_SPLIT_OF  # noqa: E402
from genomeos.genome import epigenome, reader  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "astroreg_activity_nogo"
REGISTRATION = Path("data/results/astroreg_registration.json")
REGISTRATION_SHA256 = "43f2edfab0636e0f07a035aedf5ddf7c57aa1dfcee8d39ea32e9338edd522c7a"
CELL = "astrocyte"
MARK = "H3K27ac"

#: This script, as the entry whose transitive import closure is the counting path of its result.
ENTRY = "scripts/astroreg_activity_gate.py"

OWN_CODE = (
    "genomeos/attribution/astrorun.py",
    "scripts/astroreg_activity_gate.py",
    "tests/test_astrorun.py",
)


def quantiles(xs: list[float]) -> dict[str, Any]:
    """The shape of a distribution, enough to see a scale from. Never a conversion between two."""
    if not xs:
        return {"n": 0}
    s = sorted(xs)

    def at(p: float) -> float:
        return s[min(len(s) - 1, int(p * len(s)))]

    return {
        "n": len(s),
        "min": round(s[0], 4),
        "p05": round(at(0.05), 4),
        "p25": round(at(0.25), 4),
        "median": round(at(0.50), 4),
        "p75": round(at(0.75), 4),
        "p95": round(at(0.95), 4),
        "max": round(s[-1], 4),
        "mean": round(statistics.fmean(s), 4),
        "exactly_zero": sum(1 for v in s if v == 0.0),
    }


def benchmark_activity() -> dict[str, Any]:
    """(a) The frozen activity columns in the benchmark, per cell type: how many pairs carry them.

    This is the comparison the astrocyte side is measured against: the benchmark hands the frozen
    model both numbers for every pair it covers, in every one of its cell types.
    """
    out: dict[str, Any] = {"per_file": {}}
    dists: dict[str, list[float]] = {}
    for name in (crispri.TRAINING, crispri.HELDOUT):
        path = crispri.KNOWLEDGE / name
        per_cell: dict[str, dict[str, int]] = defaultdict(lambda: {"pairs": 0, "with_both": 0})
        with gzip.open(path, "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idx = {c: i for i, c in enumerate(header)}
            missing = [c for c in astrorun.FROZEN_ACTIVITY_COLUMNS if c not in idx]
            if missing:
                raise KeyError(f"{name} has no column {missing}")
            cell_col = idx.get("CellType")
            for line in fh:
                row = line.rstrip("\n").split("\t")
                cell = row[cell_col] if cell_col is not None else "unknown"
                per_cell[cell]["pairs"] += 1
                vals: list[float | None] = []
                for c in astrorun.FROZEN_ACTIVITY_COLUMNS:
                    raw = (row[idx[c]] or "").strip()
                    vals.append(None if raw in ("", "NA", "nan") else float(raw))
                if all(v is not None for v in vals):
                    per_cell[cell]["with_both"] += 1
                if name == crispri.TRAINING and cell == "K562":
                    for c, v in zip(astrorun.FROZEN_ACTIVITY_COLUMNS, vals, strict=True):
                        if v is not None:
                            dists.setdefault(c, []).append(v)
        out["per_file"][name] = {
            c: dict(v) | {"share_with_both": v["with_both"] / v["pairs"] if v["pairs"] else None}
            for c, v in sorted(per_cell.items())
        }
    out["k562_training_distribution"] = {c: quantiles(v) for c, v in dists.items()}
    out["reading"] = (
        "the benchmark hands the frozen model both activity columns for every pair of every cell "
        "type it covers, which is the input the HCT116 arm had and an astrocyte pair does not"
    )
    return out


def screen_has_no_activity_column() -> dict[str, Any]:
    """(b) The screen's own workbook: its sheets, and whether any of them could carry activity."""
    import re
    import zipfile

    with zipfile.ZipFile(TABLE3) as z:
        names = re.findall(r'name="([^"]+)" sheetId', z.read("xl/workbook.xml").decode())
    return {
        "workbook": str(TABLE3),
        "sheets": names,
        "sheet_that_could_carry_activity": None,
        "reading": (
            "Supplementary Table 3 holds pair differential expression, negative and positive "
            "controls, Nanostring probes and results, and power calculations. There is no chromatin "
            "or activity sheet, so neither frozen activity column can be read from the screen "
            "itself, whatever the permitted-column list allowed"
        ),
    }


def astrocyte_peaks_on_disk(chroms: list[str]) -> dict[str, Any]:
    """(c) What astrocyte activity is on disk: the two registered accessions, their units, and (e)
    the accession each per-chromosome file records in its own first line."""
    out: dict[str, Any] = {}
    groups = (
        ("dnase", [reader.peaks_path(CELL, c) for c in chroms]),
        ("h3k27ac", [epigenome.peaks_path(CELL, MARK, c) for c in chroms]),
    )
    for key, paths in groups:
        headers = {}
        for p in paths:
            with gzip.open(p, "rt") as fh:
                headers[p.name] = fh.readline().strip()
        accessions = sorted({h.split()[1] for h in headers.values() if h.startswith("#")})
        registered = astroreg.ACTIVITY_INPUTS[key]
        out[key] = {
            "registered_file_accession": registered["file_accession"],
            "registered_experiment": registered["experiment"],
            "files": len(paths),
            "files_whose_first_line_records_an_accession": sum(
                1 for h in headers.values() if h.startswith("#")
            ),
            "accessions_found_in_first_lines": accessions,
            "matches_the_registration": accessions == [registered["file_accession"]],
            "one_example_first_line": headers[paths[0].name],
            "units": "narrowPeak column 7, signalValue "
            + (
                "(the peak caller's signal for a DNase peak)"
                if key == "dnase"
                else "(fold change over control for a ChIP-seq replicated peak)"
            ),
            "in_frozen_units": False,
            "written_by": "genomeos.genome.reader.fetch_peaks (float(f[6]))"
            if key == "dnase"
            else "genomeos.genome.epigenome.fetch_mark_peaks",
        }
    out["no_rpm_anywhere"] = (
        "nothing in this repository computes reads per million: `DHS.RPM` and `H3K27ac.RPM` are read "
        "from the benchmark's own columns and never derived. The epigenome layer's continuous track "
        f"is `{epigenome.SIGNAL_OUTPUT}`, not RPM, and it holds no DNase profile for any biosample "
        f"(its marks are {sorted(epigenome.MARKS)}), so there is no depth-normalised astrocyte "
        "coverage to normalise either. What continuous astrocyte signal there is is counted "
        "separately under `c2_the_continuous_track_does_not_rescue_it`, because a censoring argument "
        "that ignored it would be answerable"
    )
    out["accession_finding"] = (
        "the ENCODE accessions ARE recorded in the repository: every one of these per-chromosome "
        "files carries its accession in its own first line. The registration's left_undone entry "
        "saying they 'are not recorded anywhere in the repository' is the stale one, and "
        "terms.activity_inputs is correct; terms.activity_inputs.correction already says so inside "
        "the same frozen artefact, so the two disagree with each other and the correction is right"
    )
    return out


def continuous_signal_on_disk() -> dict[str, Any]:
    """(c2) The one answer to the censoring argument, measured and then answered on units instead.

    The peak files the registration names are censored, but this project also keeps a CONTINUOUS
    200 bp-binned track for some biosample-mark pairs, and astrocyte H3K27ac is one of them. An
    uncensored astrocyte H3K27ac quantity therefore does exist on disk, and a gate that rested only
    on censoring would be answerable by reaching for it. It is counted here so the record says so,
    and the verdict does not move, for two reasons it also measures: that track is
    `fold change over control`, not reads per million, so the units defect survives; and the signal
    store holds no DNase profile for ANY biosample, so the other half of the frozen activity term has
    no continuous source at all, censored or otherwise.
    """
    coverage = epigenome.signal_coverage()
    astro = {m: len(cs) for m, cs in sorted(coverage.get(CELL, {}).items()) if cs}
    biosamples_with_any = sorted(c for c, marks in coverage.items() if any(marks.values()))
    return {
        "store": str(epigenome.signal_path(CELL, MARK, "chr1").parent),
        "units": epigenome.SIGNAL_OUTPUT,
        "bin_bp": epigenome.BIN,
        "marks_the_store_can_hold": sorted(epigenome.MARKS),
        "dnase_profiles_in_the_store": 0,
        "astrocyte_chromosomes_per_mark": astro,
        "astrocyte_h3k27ac_chromosomes": len(coverage.get(CELL, {}).get(MARK, [])),
        "biosamples_with_any_continuous_track": len(biosamples_with_any),
        "in_frozen_units": False,
        "reading": (
            f"an uncensored continuous astrocyte {MARK} track IS on disk, "
            f"{len(coverage.get(CELL, {}).get(MARK, []))} chromosomes of 200 bp-binned "
            f"{epigenome.SIGNAL_OUTPUT}, so the censoring defect could be avoided for one of the two "
            "columns. It changes nothing, on two counts measured here. UNITS: fold change over "
            "control is not reads per million, and the frozen weight was fitted on reads per "
            "million, so the substitution is still a different feature. DNASE: the signal store "
            f"holds only {sorted(epigenome.MARKS)} and no DNase profile for any of its "
            f"{len(biosamples_with_any)} biosamples, so DHS.RPM has no continuous astrocyte source "
            "to be rescaled from at all -- its only candidate on disk remains a censored peak call. "
            "The missing input is missing in every form the repository holds"
        ),
    }


def censoring(pairs: list[dict], chroms: list[str]) -> dict[str, Any]:
    """(d) How many of the screen's elements a peak call could only give zero to.

    The benchmark's RPM is defined for every element and descends to zero continuously. A peak call
    has no value outside a called peak, so an element overlapping no peak can only be given zero.
    This counts those elements, by label, so the size of the fabrication is on the record.
    """
    per_label: dict[str, dict[str, int]] = defaultdict(
        lambda: {"elements": 0, "no_dnase_peak": 0, "no_h3k27ac_peak": 0, "neither": 0}
    )
    signal: dict[str, list[float]] = {"dnase": [], "h3k27ac": []}
    seen: set[tuple[str, str]] = set()
    for chrom in chroms:
        dnase = reader.PeakIndex(reader.load_peaks(CELL, chrom))
        marks = epigenome.load_mark_peaks(CELL, MARK, chrom)
        if marks is None:
            raise FileNotFoundError(f"no astrocyte {MARK} peaks for {chrom}")
        k27 = reader.PeakIndex(marks)
        for p in pairs:
            if p["chrom"] != chrom:
                continue
            key = (chrom, p["element"])
            if key in seen:
                continue
            seen.add(key)
            d = dnase.overlapping(p["start"], p["end"])
            h = k27.overlapping(p["start"], p["end"])
            signal["dnase"].extend(v for _, _, v in d)
            signal["h3k27ac"].extend(v for _, _, v in h)
            for bucket in (per_label[p["label"]], per_label["all labels"]):
                bucket["elements"] += 1
                if not d:
                    bucket["no_dnase_peak"] += 1
                if not h:
                    bucket["no_h3k27ac_peak"] += 1
                if not d and not h:
                    bucket["neither"] += 1
    return {
        "by_label": {
            k: v | {"share_with_no_h3k27ac_peak": round(v["no_h3k27ac_peak"] / v["elements"], 4)}
            for k, v in sorted(per_label.items())
        },
        "astrocyte_peak_signal_at_the_screens_elements": {k: quantiles(v) for k, v in signal.items()},
        "reading": (
            "an element overlapping no called peak can be given no value but zero, and zero in the "
            "frozen feature means log_activity = 0 and activity_over_distance = -log(distance). "
            "That is a number invented for the arm the gain is subtracted by, not a measurement of "
            "it. The benchmark's own K562 elements carry a continuous low tail instead, and the "
            "share that is exactly zero is reported beside this"
        ),
    }


def manifest(chroms: list[str], pairs: list[dict]) -> dict[str, Any]:
    dnase_files = [reader.peaks_path(CELL, c) for c in chroms]
    mark_files = [epigenome.peaks_path(CELL, MARK, c) for c in chroms]
    return {
        "sources": [
            {
                "accession": "data/results/astroreg_registration.json",
                "version": f"the frozen registration, sha256 {REGISTRATION_SHA256}; its terms are "
                "quoted, never moved",
            },
            {
                "accession": "Green NFO, et al. CRISPRi screening in cultured human astrocytes. "
                "Nature Neuroscience 2025;29(3):703-716 -- Supplementary Table 3",
                "version": "publisher open-access supplementary store; pinned here by sha256",
                "url": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-025-02154-3"
                "/MediaObjects/41593_2025_2154_MOESM5_ESM.xlsx",
            },
            {
                "accession": "EngreitzLab/CRISPR_comparison, EPCrisprBenchmark training_K562 and "
                "heldout_5_cell_types (Gschwind et al.)",
                "version": "fetched 2026-09-16; pinned here by sha256",
                "url": crispri.BASE_URL,
            },
            {
                "accession": "ENCODE ENCFF874OPW (ENCSR000EPM) astrocyte DNase-seq peaks and "
                "ENCFF970DKF (ENCSR000AOQ) astrocyte H3K27ac replicated peaks",
                "version": "released GRCh38 narrowPeak; each per-chromosome derivation records its "
                "accession in its own first line, and this result checks that it does",
            },
        ],
        "inputs": [
            mf.input_entry(
                REGISTRATION,
                partition=None,
                role="the frozen registration: its feasibility rule and activity_inputs are read "
                "and its sha256 is checked; no term is moved",
            ),
            mf.input_entry(
                TABLE3,
                partition=None,
                role="the screen's element coordinates and registered labels, and the workbook's "
                "sheet names; no effect size, p-value, FDR or score was read",
            ),
            mf.input_entry(
                crispri.KNOWLEDGE / crispri.TRAINING,
                partition=CRISPRI_SPLIT_OF[crispri.TRAINING],
                role="the frozen activity columns DHS.RPM and H3K27ac.RPM, for the availability "
                "count and the K562 distribution; already-read development evidence",
            ),
            mf.input_entry(
                crispri.KNOWLEDGE / crispri.HELDOUT,
                partition=CRISPRI_SPLIT_OF[crispri.HELDOUT],
                role="the same two columns per held-out cell type, the comparison the astrocyte "
                "side is measured against; already-read",
            ),
            mf.files_entry(
                "data/results/dnase_astrocyte_chr*.bed.gz",
                dnase_files,
                partition=None,
                role="ENCFF874OPW astrocyte DNase narrowPeak per chromosome: signalValue and the "
                "accession in each first line",
            ),
            mf.files_entry(
                "data/knowledge/epigenome/peaks/astrocyte_H3K27ac_chr*.bed.gz",
                mark_files,
                partition=None,
                role="ENCFF970DKF astrocyte H3K27ac replicated-peak narrowPeak per chromosome: "
                "signalValue and the accession in each first line",
            ),
        ],
        "assembly": "GRCh38",
        "coordinates": (
            "half-open, zero-based for the peak files as the narrowPeak derivations store them. The "
            "screen states no convention for its chrom:start-end strings; the only use made of them "
            "here is an overlap test against peak calls, whose outcome a one-base shift cannot "
            "change for any element this result counts"
        ),
        "parameters": {
            "frozen_activity_columns": list(astrorun.FROZEN_ACTIVITY_COLUMNS),
            "authorised_requests": astrorun.AUTHORISED_REQUESTS,
            "requests_permitted_by_this_gate": 0,
            "chromosomes": chroms,
            "screen_pairs": len(pairs),
            "alphagenome_requests": 0,
            "model_requests": 0,
            "astrocyte_pairs_scored": 0,
            "money_spent": 0,
        },
        "exclusions": [
            "no AlphaGenome request is sent and this script has no code path that could send one",
            "no astrocyte deletion value is computed, requested or read",
            "no effect size, fold change, p-value, FDR or expression level is read from the screen",
            "no substitute activity value is constructed: the whole point of the result is that "
            "constructing one is what the registration forbids",
            "the continuous signal store is counted by file name only: no profile is opened and no "
            "value is read from one, because the question asked of it is which biosamples and marks "
            "it covers and in what units, not what it says at any element",
            "the 1.4 GB all-element deletion table is not read: coverage is the registration's "
            "already-recorded figure and the activity question applies to covered and uncovered "
            "pairs alike",
        ],
        "partitions": {
            "screen_pairs": f"all {len(pairs)} element-gene pairs of the published screen",
            "k562": "already-read development evidence, used only for the activity-column "
            "availability count and the K562 distribution",
        },
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    t0 = time.time()
    digest, _, _ = mf.sha256_of(REGISTRATION)
    if digest != REGISTRATION_SHA256:
        print(f"refused: {REGISTRATION} is sha256 {digest}, not the frozen {REGISTRATION_SHA256}")
        return 2
    frozen = json.loads(REGISTRATION.read_text())
    rule = frozen["terms"]["feasibility_gate"]["rule"]
    if rule != astrorun.REGISTERED_MISSING_INPUT_RULE:
        print("refused: the rule quoted in code is not the rule in the frozen registration")
        return 2

    print("reading the screen's elements and registered labels")
    pairs = screen_pairs()
    chroms = sorted({p["chrom"] for p in pairs}, key=lambda c: (len(c), c))
    print(f"  {len(pairs)} pairs over {len(chroms)} chromosomes")

    print("(a) the frozen activity columns in the benchmark")
    bench = benchmark_activity()
    print("(b) the screen's own workbook")
    screen = screen_has_no_activity_column()
    print("(c,e) the astrocyte activity on disk, and its accessions")
    disk = astrocyte_peaks_on_disk(chroms)
    print("(c2) the continuous signal store, the one answer to the censoring argument")
    cont = continuous_signal_on_disk()
    print("(d) what a peak call could only give zero to")
    cens = censoring(pairs, chroms)

    verdict = astrorun.activity_verdict(
        {
            "DHS.RPM": {
                "source": f"{disk['dnase']['registered_file_accession']} astrocyte DNase narrowPeak "
                "(the registration's own activity input)",
                "units": disk["dnase"]["units"],
                "in_frozen_units": False,
                "why_not": "a peak signalValue is not reads per million, and no element outside a "
                "called peak carries any value at all",
            },
            "H3K27ac.RPM": {
                "source": f"{disk['h3k27ac']['registered_file_accession']} astrocyte H3K27ac "
                "replicated-peak narrowPeak (the registration's own activity input)",
                "units": disk["h3k27ac"]["units"],
                "in_frozen_units": False,
                "why_not": "fold change over control is not reads per million, and no element "
                "outside a called peak carries any value at all. The continuous astrocyte H3K27ac "
                "track on disk lifts the censoring half of that and not the units half: it is also "
                f"{epigenome.SIGNAL_OUTPUT}",
            },
        }
    )
    permitted = astrorun.requests_permitted(verdict)

    payload = {
        "status": "a gate on one input, decided before any money was spent. No AlphaGenome request "
        "was sent, no astrocyte pair was scored, nothing was bought",
        "lane": "lane-astro",
        "question": "can the frozen model's activity term be supplied for the astrocyte pairs in "
        "the units its weights were fitted on, so that the registered primary claim can be computed "
        "from the 1,322 requests Albert authorised",
        "authorisation": astrorun.AUTHORISATION,
        "authorised_requests": astrorun.AUTHORISED_REQUESTS,
        "requests_sent": 0,
        "money_spent": 0,
        "registration": {
            "path": str(REGISTRATION),
            "sha256": REGISTRATION_SHA256,
            "byte_identical_to_the_authorised_copy": True,
            "primary_claim": frozen["terms"]["primary_claim"],
            "feasibility_rule": rule,
            "not_claimed": frozen["terms"]["not_claimed"],
        },
        "frozen_activity_term": {
            "columns": list(astrorun.FROZEN_ACTIVITY_COLUMNS),
            "units": astrorun.FROZEN_ACTIVITY_UNITS,
            "read_in_exactly_one_place": "genomeos/attribution/crispri.py, crispri.parse "
            "(dhs=float(r['DHS.RPM'] or 0), h3k27ac=float(r['H3K27ac.RPM'] or 0))",
            "features_that_depend_on_it": [
                f for f in crispri.FEATURES["activity + distance"] if f != "log_distance"
            ],
            "in_both_arms": "'activity + distance + deletion' and 'activity + distance' both carry "
            "log_activity and activity_over_distance, so the term is in the minuend and the "
            "subtrahend of the only registered claim",
        },
        "a_benchmark_supplies_it": bench,
        "b_the_screen_does_not": screen,
        "c_what_is_on_disk_instead": disk,
        "c2_the_continuous_track_does_not_rescue_it": cont,
        "d_and_it_is_censored": cens,
        "verdict": verdict,
        "requests_permitted": permitted,
        "decision": (
            "NO-GO. The run is not made and nothing is spent. This is the registration's own "
            "feasibility rule applied, not a new rule and not a term moved: the frozen model needs "
            "DHS.RPM and H3K27ac.RPM, the astrocyte side has neither, and the only candidates on "
            "disk are a peak signalValue and a fold change over control, which would be substitutes "
            "for a missing input"
        ),
        "what_the_money_would_and_would_not_have_bought": (
            "the 1,322 requests would have bought exactly what they were costed for, an astrocyte "
            "deletion value per registry element. They would not have bought the registered claim, "
            "because the gain is a difference between two arms that both need an activity term the "
            "astrocyte side cannot supply. Spending first and discovering this afterwards would "
            "have left a cache and no reading"
        ),
        "what_would_clear_it": [
            "an activity input in the frozen units for astrocyte: a depth-normalised coverage of "
            "the two ENCODE experiments, from which reads per million in an element is computable. "
            "That is a different file from either registered accession and is not on disk, so it is "
            "an amendment to the registration's activity_inputs, which only the registering lane "
            "and Albert may make",
            "or a registered second claim that does not contain an activity term at all, whose "
            "arms are both scorable from what is on disk. That is a new registration written before "
            "any astrocyte score exists, not a reading of this one",
        ],
        "not_claimed": list(astroreg.NOT_CLAIMED),
        "left_undone": [
            "no AlphaGenome request was sent, so nothing is known about the astrocyte deletion "
            "value beyond the registration's feasibility gate: the track exists and 224 chr21 rows "
            "name it as a winning track",
            "whether a depth-normalised astrocyte coverage exists on the ENCODE portal for these "
            "two experiments is not established here; no portal request was made",
            "the registration's stale left_undone entry is reported and NOT edited: the artefact is "
            "frozen and the wording is the registering lane's to correct",
        ],
    }
    payload[mf.KEY] = manifest(chroms, pairs)
    path = save_result(RESULT, payload)

    print(
        json.dumps(
            {
                "verdict_go": verdict["go"],
                "columns_not_in_frozen_units": verdict["columns_not_in_frozen_units"],
                "requests_permitted": permitted,
                "requests_sent": 0,
                "money_spent": 0,
                "elements_with_no_h3k27ac_peak": cens["by_label"]["all labels"]["no_h3k27ac_peak"],
                "elements_total": cens["by_label"]["all labels"]["elements"],
                "continuous_astrocyte_h3k27ac_chromosomes": cont["astrocyte_h3k27ac_chromosomes"],
                "continuous_dnase_profiles_anywhere": cont["dnase_profiles_in_the_store"],
                "accessions_match_the_registration": {
                    k: disk[k]["matches_the_registration"] for k in ("dnase", "h3k27ac")
                },
            },
            indent=1,
        )
    )
    print(f"({time.time() - t0:.0f} s) -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
