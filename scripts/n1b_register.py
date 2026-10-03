#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Write N1b's registration, which is committed ALONE and before any byte of `X` is read.

The file this writes carries no field a run would fill. There is no `result`, no `verdict`, no
`figures`, no `passed` and no key holding null awaiting a number: those fields do not EXIST in it,
rather than sitting empty, so that nothing in it can read as a preview of an answer. The file
authorises no conclusion. It states what will be measured, the rule that will decide, and what a
fail means, and that is the whole of its content.
"""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import n1b_calibration as n1b  # noqa: E402
from genomeos.results import save_result  # noqa: E402
from scripts import n1b_fetch  # noqa: E402

NAME = "n1b_calibration_registration"
OUT = ROOT / f"data/results/{NAME}.json"
#: This script, as the entry whose transitive import closure is the counting path of its result.
ENTRY = "scripts/n1b_register.py"
OWN_CODE = (
    "genomeos/attribution/n1b_calibration.py",
    "scripts/n1b_fetch.py",
    "scripts/n1b_register.py",
    "scripts/n1b_run.py",
    "tests/test_n1b_calibration.py",
)
FROZEN = ("genomeos/attribution/n1b_calibration.py", "tests/test_n1b_calibration.py")
PLUMBING = ("scripts/n1b_fetch.py", "scripts/n1b_register.py", "scripts/n1b_run.py")

ALBERT = "(8) open N1b calibration-only."


def sha256_text(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def last_commit(path: str) -> str | None:
    out = subprocess.run(  # noqa: S603
        ["git", "log", "-1", "--format=%H", "--", path],  # noqa: S607
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    return out.stdout.strip() or None


def uncommitted(paths: tuple[str, ...]) -> list[str]:
    out = subprocess.run(  # noqa: S603
        ["git", "diff", "HEAD", "--name-only", "--", *paths],  # noqa: S607
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    return [p for p in out.stdout.split() if p]


def power(clusters: int, genes_per_stratum: int) -> dict[str, Any]:
    """Stated IN ADVANCE, in the units of the cluster, which is the non-targeting ROW.

    The honest bound is a range, because the design effect -- how much genes within one row move
    together -- is not known before the run. Both ends are given, and the bootstrap reports which end
    the data sits at. One leg is strong and the other is weak, and saying which before the run is the
    point of this field.
    """
    m = clusters
    per_gene_sd = math.sqrt(2.0 / (m - 1))  # SD of a per-gene variance estimate at m rows, null
    median_sd_independent = per_gene_sd * 1.2533 / math.sqrt(genes_per_stratum)
    tail_se_independent = math.sqrt(n1b.NOMINAL_TAIL * (1 - n1b.NOMINAL_TAIL) / (m * genes_per_stratum))
    tail_se_worst = math.sqrt(n1b.NOMINAL_TAIL * (1 - n1b.NOMINAL_TAIL) / m)
    return {
        "unit": "non-targeting row (the cluster); cells pool within a row and genes correlate within one",
        "clusters": m,
        "genes_per_stratum": genes_per_stratum,
        "values_per_stratum": m * genes_per_stratum,
        "variance_leg": {
            "per_gene_sd_of_the_variance_estimate": per_gene_sd,
            "deciding_summary_sd_if_genes_were_independent": median_sd_independent,
            "deciding_summary_sd_if_genes_moved_together_perfectly": per_gene_sd,
            "band_half_widths": [1.0 - n1b.VAR_LO, n1b.VAR_HI - 1.0],
            "reading": "even at the pessimistic end the band edges sit 3.2 and 4.0 SD from 1, so a "
            "median per-gene variance genuinely at 0.80 or 1.25 is distinguished from 1 with power "
            "above 0.999. What this study CANNOT resolve, at that end, is a deviation smaller than "
            "about 0.06 in the deciding summary.",
        },
        "tail_leg": {
            "se_if_values_were_independent": tail_se_independent,
            "se_if_genes_moved_together_perfectly": tail_se_worst,
            "band_half_width": 0.01,
            "distance_to_band_edge_in_worst_case_se": 0.01 / tail_se_worst,
            "power_at_the_band_edge_in_the_worst_case": 0.36,
            "reading": "THIS IS THE WEAK LEG, said before the run. At the pessimistic end the band "
            "edge is about 1.0 SE away, so a tail share genuinely at 0.04 or 0.06 would be told from "
            "0.05 with power near 0.36. The cluster bootstrap measures where between the two ends the "
            "data sits; if the interval spans the band, the reading is that THE DATA CANNOT TELL, and "
            "that is not the same as the band holding.",
        },
    }


def registration(probe: dict[str, Any]) -> dict[str, Any]:
    nt = probe["non_targeting_rows"]
    genes = probe["genes_total"]
    per_stratum = genes // n1b.STRATA
    return {
        "study": n1b.STUDY,
        "registered": "2026-10-03",
        "lane": "lane-n1b",
        "question": "Is the derived per-gene statistic T = X * sqrt(num_cells_filtered) CALIBRATED? "
        "Nothing about factors. Nothing about whether anything responds. A calibration check on the "
        "null, and the whole of the study.",
        "albert_approved": ALBERT,
        "not_a_resumption_of_n1": (
            "This is a NEW registered study. It is deliberately NOT a resumption of N1. N1 was closed "
            "on 2026-10-01 under its own registered stop rule, on the external reviewer's advice "
            "relayed by Albert, and NOTHING here reopens that closure or revisits its finding. N1's "
            "registration, its two amendments and its committed result are not read by this study's "
            "code, are not altered by it, and none of its thresholds is inherited. The one thing "
            "shared is the SOURCE constant naming the file, so the project holds one copy of the file "
            "identity. N1's own gate read |T| >= 3 against a 1.5x ratio ceiling over ten deciles; "
            "this study's rule is a variance band and a 1.96 tail band over five quintiles, and "
            "neither rule's verdict carries to the other."
        ),
        "authorises_no_conclusion": (
            "This file states a plan and a rule. It holds no measurement, no figure and no verdict, "
            "and the fields a run would fill DO NOT EXIST in it rather than sitting empty, so that "
            "nothing in it can read as a preview. It authorises no conclusion of any kind. The check "
            "that says so runs twice: once on the payload, and once on the file AS WRITTEN, because "
            "genomeos/results.py adds the registry's own stamp after the payload leaves the writer. "
            "That stamp uses one of the forbidden spellings -- `result` -- and it holds this "
            "registration's own NAME, 'n1b_calibration_registration', and not an outcome. It is "
            "allowed at exactly that value and nothing else, so a future stamp putting a figure there "
            "would stop the write rather than pass unnoticed. `date` is the same stamp's write date."
        ),
        "data": {
            "study": "Replogle et al. 2022, genome-scale Perturb-seq in K562",
            "article": "Figshare+ 20029387",
            "doi": "10.25452/figshare.plus.20029387.v1",
            "paper_doi": "10.1016/j.cell.2022.05.013",
            "licence": "CC BY 4.0",
            "cost": "none. A public CC BY 4.0 dataset read by HTTP range request. 0 AlphaGenome "
            "requests, 0 model requests, no money.",
            "files": {
                "normalized": {
                    **n1b.source_file("normalized"),
                    "role": "X, the statistic's numerator, and the obs columns",
                    "read_how": "HTTP range requests only; the file is never downloaded whole",
                    "md5_whole_file": n1b.MD5_NOT_VERIFIABLE,
                },
                "raw": {
                    **n1b.source_file("raw"),
                    "role": "the stratifier only: per-gene control expression. No figure of the pass "
                    "rule is computed from it.",
                    "read_how": "already present in data/cache/n1 from N1's download; read locally in "
                    "blocks of 64 rows",
                    "md5_verified_against_figshare": True,
                    "local_path": "data/cache/n1/K562_gwps_raw_bulk_01.h5ad",
                    "local_sha256": n1b_fetch.sha256_file(n1b_fetch.RAW_LOCAL),
                },
            },
            "cache_rule": "every byte written lands under data/cache/, which is NEVER committed: "
            "Albert has ruled 'data/cache stays local'. The result declares it by path and sha256.",
        },
        "rows": {
            "measured_from": "the normalized file's own index (obs/gene_transcript) and its own obs "
            "columns, read by range request before this registration was written. Nothing counted "
            "here was taken on anyone's word.",
            "rows_in_file": probe["rows_total"],
            "genes_in_file": probe["genes_total"],
            "non_targeting": probe["non_targeting_count"],
            "core_controls_among_them": probe["core_control_count_among_non_targeting"],
            "core_controls_in_whole_file": probe["core_control_count_overall"],
            "usable": probe["usable_non_targeting"],
            "unusable": probe["unusable_non_targeting"],
            "against_the_counts_i_was_given": (
                "I was told about 585 non-targeting rows of which about 514 are the core set. Measured: "
                f"{probe['non_targeting_count']} non-targeting and "
                f"{probe['core_control_count_among_non_targeting']} core. Both match. A third fact was "
                "not in the brief and is measured here: obs/num_cells_filtered is non-finite on exactly "
                f"the {probe['unusable_non_targeting']} non-core rows and finite on all "
                f"{probe['usable_non_targeting']} core ones, so the usable set IS the core set and is "
                "determined by an obs field rather than chosen by this lane. core_control is true "
                "nowhere outside the non-targeting rows."
            ),
            "selection": "every row whose label carries 'non-targeting', from the index alone. No row "
            "is selected or dropped by anything measured on it. NO FACTOR ROW IS READ, AT ALL.",
        },
        "statistic": {
            "definition": "T = X * sqrt(num_cells_filtered), per (non-targeting row, gene)",
            "x": "the normalized pseudobulk value: per the producers, a row's mean over its cells of "
            "each cell's gemgroup z-score",
            "under_the_null": "T ~ N(0, 1) exactly, if the three instrument assumptions hold",
        },
        "pass_rule": {
            "text": n1b.PASS_RULE,
            "variance_band": [n1b.VAR_LO, n1b.VAR_HI],
            "variance_nominal": n1b.NOMINAL_VAR,
            "tail_z": n1b.TAIL_Z,
            "tail_band": [n1b.TAIL_LO, n1b.TAIL_HI],
            "tail_nominal": n1b.NOMINAL_TAIL,
            "both_must_hold": "either failing, in any one scope, fails the whole rule",
            "deciding_variance_summary": n1b.PRIMARY_VAR_SUMMARY,
            "deciding_tail_summary": n1b.PRIMARY_TAIL_SUMMARY,
            "summaries_reported_beside_and_deciding_nothing": [
                *n1b.SECONDARY_VAR_SUMMARIES,
                *n1b.SECONDARY_TAIL_SUMMARIES,
            ],
            "why_a_summary_had_to_be_named": "'the per-gene empirical variance of T' is a quantity per "
            "gene and a band is read against one number. Naming the deciding summary after the figures "
            "were in hand would be choosing the one that passes.",
            "a_known_bias_the_band_carries": "at ddof=1 the per-gene variance of a perfectly calibrated "
            f"T has expectation 1, but a balanced two-point T gives m/(m-1) = "
            f"{probe['usable_non_targeting'] / (probe['usable_non_targeting'] - 1):.5f} at "
            f"{probe['usable_non_targeting']} rows. It is inside the band by three orders of magnitude "
            "and is recorded rather than corrected.",
        },
        "stratification": {
            "registered": f"quintiles ({n1b.STRATA}) of control expression, genes ranked by it, sizes "
            "differing by at most one; stratum 0 holds the least expressed genes",
            "why": "the derivation is EXPECTED to fail on sparse genes, where X is a mean of z-scores "
            "over few non-zero counts and the normal approximation is weakest. An unstratified pass "
            "would hide the thing most worth knowing, so the rule is read in every quintile and an "
            "overall pass does not excuse a failing one.",
            "the_name_i_was_given_and_what_it_turned_out_to_be": (
                "I was told to stratify by quintiles of `control_expr`, that its NAME was known to me "
                "and its VALUES were not outcome rows. Stating it explicitly rather than leaving a "
                "reader to wonder: obs/control_expr IS in the file, it IS a per-ROW field and not a "
                "per-gene one -- it is the control expression of the row's own target gene -- and it is "
                f"NaN on ALL {probe['non_targeting_count']} non-targeting rows, measured, because a "
                "non-targeting row has no target gene. "
                f"{probe['control_expr_finite_among_non_targeting']} of them carry a finite value. "
                "Quintiles of it therefore CANNOT BE FORMED over these rows, and this is reported as a "
                "difference from the brief rather than worked around silently."
            ),
            "what_is_used_instead": (
                "per-gene control expression, measured as the cell-weighted mean of RAW X over the "
                "usable non-targeting rows: e_j = sum_i(raw_x[i][j] * n_i) / sum_i(n_i). It is control "
                "expression in the sense the stratification needs -- the level at which a gene is seen "
                "in control cells -- and it is per gene, which is the axis sparse genes live on. It is "
                "measured on non-targeting rows only and carries no factor response."
            ),
            "it_is_not_an_outcome": "the stratifier is an expression level of control cells. It is not "
            "a response, not a factor row, and no leg of the pass rule is computed from the raw file.",
            "genes_that_leave": "a gene with any non-finite T, or a non-finite control expression, "
            "cannot have a band read on it: it leaves before any figure is formed and is counted.",
        },
        "bands_are_relative": {
            "kind": "relative, to the derivation's own nominal value",
            "reason": "the nominal is a THEORETICAL CONSTANT -- exactly 1 for the variance and "
            f"{n1b.NOMINAL_TAIL:.7f} for the tail -- with no sampling error of its own, so relative and "
            "absolute coincide here. The bands are stated relative because that is what a reader checks.",
            "control": n1b.CONTROL_KIND,
            "no_control_group_is_possible": "the null IS the reference. A second null would be another "
            "draw from the same rows, so there is nothing a control group could add.",
        },
        "every_ratio_beside_its_level": (
            "Every ratio this study reports is emitted in ONE object with its absolute level and its "
            "control, and cannot be separated from them by a reader or by a later summary. The reason "
            "is a result of 2026-10-02 in this project that reported a 2.9-4.1x relative excess whose "
            "absolute level was 6.0%, so 94% of pairs disagreed and the ratio read as a finding the "
            "level did not support."
        ),
        "instrument_assumptions": list(n1b.ASSUMPTIONS),
        "assumptions_are_printed_beside_every_figure": (
            "all four travel in the `assumptions` field of every interval the run emits. The first is "
            "UNCONFIRMED: the producers document X as a per-row mean of per-cell gemgroup z-scores but "
            "do not document that the divisor is the same num_cells_filtered the obs column reports. A "
            "figure that fails only because that is false would look exactly like a derivation that is "
            "wrong, and this study CANNOT tell those apart. It says so rather than choosing."
        ),
        "intervals": {
            "cluster": "the non-targeting ROW",
            "why_the_row": "cells pool within a row and genes correlate within a row. A resample of "
            f"(row, gene) values would treat {probe['usable_non_targeting'] * genes:,} correlated "
            "values as independent and return an interval the data has not earned.",
            "method": "95% percentile interval over "
            f"{n1b.BOOTSTRAP_B} resamples of the rows, with replacement, seed {n1b.BOOTSTRAP_SEED}",
            "identical_resample_share": "reported for EVERY interval, degenerate or not. It is the one "
            "number that says whether an interval came from variation in the data or from its absence.",
            "degenerate_route": "at 0 or all successes a cluster bootstrap returns the identical value "
            "every time and DECIDES NOTHING. The tail share then routes to a Wilson score interval on "
            "the count of ROW-CLUSTERS carrying any exceedance -- never on the (row, gene) count, which "
            "would shrink the bound by a factor the data does not earn. A variance has no exact "
            f"binomial analogue, so a degenerate variance interval means '{n1b.INCONCLUSIVE_WORDING}' "
            "and the figure stands with no interval.",
        },
        "power": power(probe["usable_non_targeting"], per_stratum),
        "inconclusive_reading": {
            "wording": n1b.INCONCLUSIVE_WORDING,
            "never": "no information",
            "when": "where an interval spans the band, the reading is that the data cannot tell whether "
            "the band holds. That is not the band holding, and it is not the absence of information "
            "either: the figure and its interval are both reported.",
        },
        "cannot_establish": {
            "statements": list(n1b.CANNOT_ESTABLISH),
            "pass_licence": n1b.PASS_LICENCE,
            "why_the_licence_is_narrow": (
                "the producers' normalisation computes each per-cell z-score per gemgroup AGAINST THE "
                "NON-TARGETING CELLS -- the mean and SD of the control population. So every "
                "non-targeting row is a per-guide SUBSET of the very reference that defined the scale, "
                "and a pass tests the sqrt(n) scaling, cell independence and normality WITHIN that "
                "reference. It does NOT show that T is calibrated for perturbed rows, whose cell "
                "counts, variances and knockdown-induced shifts all differ. A pass is therefore weaker "
                "than it would read to someone skimming, and this is registered in advance rather than "
                "left to be noticed later."
            ),
            "the_caveat_travels": "decide() cannot emit a verdict without pass_licence and "
            "cannot_establish attached, so the sentence sits beside every figure rather than once in "
            "a file nobody re-reads. If step 2 ever runs, it states it again beside every "
            "perturbed-row figure.",
            "it_moves_no_band": "the variance band, the tail band, every-quintile-plus-overall, a "
            "stratum that cannot be formed counting as a failure, the row as the cluster and the four "
            "printed assumptions all stand exactly as frozen at d04261a. This bounds what a pass "
            "LICENSES, not what counts as one.",
            "perturbed_row_cell_counts": {
                **probe["cell_counts"],
                "per_factor_counts_are_not_read": "NOT in N1b, at all: not sampled, not spot-checked, "
                "not printed for a sanity check. A perturbed row's num_cells_filtered is "
                "OUTCOME-ADJACENT -- knocking down an essential gene lowers its cell count, which is a "
                "fitness phenotype and not metadata -- so reading per-factor counts would read an "
                "outcome and would cost exactly the blindness this registration exists to protect. "
                "Only a distribution over ALL perturbed rows is read, as quantiles, with no identity.",
            },
        },
        "disclosure": {
            "read_before_this_registration_was_written": [
                "the normalized file's superblock and object headers (where X lives, and its shape)",
                "obs/gene_transcript and var/gene_id: the file's own row and gene index",
                "obs/num_cells_filtered, obs/core_control and obs/control_expr, whole columns",
                probe["cell_counts"]["disclosure"],
            ],
            "when": "in scripts/n1b_fetch.py probe, by 33 HTTP range requests over "
            f"{probe['bytes_fetched']:,} bytes, before this file was written and before any byte of X",
            "never_read": "no byte of X at this point; no perturbed row's own cell count and no "
            "perturbed row's index persisted anywhere; no expression value of any perturbed row; no "
            "factor row at all",
            "standard": "disclosure, not abstinence, for the file's shape -- and abstinence for a "
            "perturbed row's identity, because that one is outcome-adjacent.",
        },
        "what_a_fail_means": {
            "written_before_the_run": True,
            "verdict": "the derivation T = X * sqrt(n) is UNUSABLE.",
            "what_this_study_then_does": "says so, and stops. It does not propose a repair, a wider "
            "band or a second statistic.",
            "fallbacks": "the authors' own per-gene results, or the single-cell file. Both are "
            "ALBERT'S DECISION and not this lane's.",
            "either_way": "nothing about factors is read. No factor row is opened on a pass or a fail.",
        },
        "outcome_exposure": {
            "recorded": "null-calibration exposure, no outcome row",
            "why": "the non-targeting rows are the NULL and carry no factor response, by construction "
            "of the experiment: the sgRNA targets nothing. Reading them spends no confirmatory outcome "
            "and leaves this project's independent outcome data untouched, because no row that could "
            "show a response is opened.",
        },
        "access_plan": {
            "x_rows": n1b_fetch.row_ranges(probe["x_layout"], nt),
            "x_layout": probe["x_layout"],
            "rows_to_read": nt,
            "why_all_585_and_not_the_514": "the 71 unusable rows are read too, at a cost of "
            f"{71 * probe['x_layout']['row_stride_bytes']:,} further bytes, so that their exclusion is "
            "auditable from fetched bytes rather than applied before the fetch. No figure is computed "
            "from them.",
            "logging": "every range request is logged with its URL, its exact byte range and the "
            "sha256 of the bytes returned.",
            "headers_already_read": {
                "bytes": probe["bytes_fetched"],
                "requests": len(probe["fetch_log"]),
                "what": "superblock, object headers, obs/gene_transcript, var/gene_id and three obs "
                "columns. No byte of X.",
            },
            "memory": "rows arrive one range request at a time and are written into a "
            f"{probe['usable_non_targeting']} x {genes} float32 matrix of "
            f"{probe['usable_non_targeting'] * genes * 4 / 1e6:.1f} MB. No file is opened whole and the "
            "raw file is read in blocks of 64 rows. A memory ceiling, if wanted, is INJECTED in tests; "
            "a live reading is recorded only in a real run, and it is valid only in a process that has "
            "done nothing substantial first.",
        },
        "frozen_code": {
            path: {"sha256": sha256_text(ROOT / path), "last_commit": last_commit(path)}
            for path in (*FROZEN, *PLUMBING)
        },
        "frozen_code_uncommitted": uncommitted((*FROZEN, *PLUMBING)),
        "cost": {
            "alphagenome_requests": 0,
            "model_requests": 0,
            "money": "none",
            "network": "the public Figshare range reads and nothing else",
        },
    }


FORBIDDEN = ("result", "verdict", "passed", "figures", "reading", "overall", "strata_figures")


def result_manifest() -> dict[str, Any]:
    """What produced these numbers, so the registration is reproducible like any other result.

    It has inputs although it measures nothing: the row counts and the byte ranges it states come
    from the probe's record of the file's own index, and from the local raw file whose sha256 it
    carries. Both live under data/cache, which is never committed, so a rebuild on a fresh checkout
    reports them ABSENT rather than differing -- which is the honest reading and not a gap to paper
    over, since Albert has ruled "data/cache stays local".
    """
    return {
        "sources": [
            {
                "accession": "Replogle et al. 2022 processed Perturb-seq, Figshare+ 20029387 "
                "(K562_gwps_normalized_bulk_01.h5ad: headers, index and three obs columns only, by "
                "HTTP range request; no byte of X)",
                "version": "10.25452/figshare.plus.20029387.v1, CC BY 4.0",
            },
            {
                "accession": "Replogle et al. 2022 processed Perturb-seq, Figshare+ 20029387 "
                "(K562_gwps_raw_bulk_01.h5ad: present locally from N1's download, md5 checked against "
                "Figshare's published one; read here for its sha256 alone)",
                "version": "10.25452/figshare.plus.20029387.v1, CC BY 4.0",
            },
            {
                "accession": "Replogle et al. 2022, Cell (the documentation of what X is)",
                "version": "10.1016/j.cell.2022.05.013",
            },
        ],
        "inputs": [
            mf.input_entry(n1b_fetch.PROBE, partition="the file's own index and obs columns"),
            mf.input_entry(n1b_fetch.RAW_LOCAL, partition="the stratifier's file, hashed not read here"),
        ],
        "assembly": "n/a: a registration over pseudobulk rows and gene columns; no coordinate is read",
        "coordinates": "n/a: no genomic interval is read",
        "parameters": {
            **n1b.CONSTANTS,
            "primary_var_summary": n1b.PRIMARY_VAR_SUMMARY,
            "primary_tail_summary": n1b.PRIMARY_TAIL_SUMMARY,
            "pass_rule": n1b.PASS_RULE,
            "x_row_bytes_planned": 585 * 32992,
        },
        "exclusions": [
            "every row that is not non-targeting. NO FACTOR ROW IS READ, AT ALL, on a pass or a fail",
            "the 71 non-core non-targeting rows are read but carry no finite num_cells_filtered, so no "
            "figure is computed from them; they are fetched so the exclusion is auditable from bytes",
            "obs/control_expr, which the brief named as the stratifier: it is NaN on all 585 "
            "non-targeting rows and quintiles of it cannot be formed",
        ],
        "partitions": "n/a: a registration committed before any measurement; nothing is evaluated in it",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    mf.trace_begin()
    probe = json.loads(n1b_fetch.PROBE.read_text())
    rec = registration(probe)
    present = set(FORBIDDEN) & set(rec)
    if present:
        raise SystemExit(f"the registration must hold no field a run would fill: {sorted(present)}")
    rec["result_manifest"] = result_manifest()
    path = save_result(NAME, rec)
    # The check runs AGAIN on the file as WRITTEN, not only on the payload, because save_result adds
    # the registry's own stamp after the payload leaves here. `result` is the one forbidden spelling
    # the stamp uses, and it holds this registration's own NAME -- not an outcome. It is allowed only
    # at exactly that value, so a future stamp that put anything else there would stop the write.
    written = json.loads((ROOT / path).read_text())
    for key in FORBIDDEN:
        if key == "result" and written.get(key) == NAME:
            continue
        if key in written:
            raise SystemExit(f"the written registration holds {key!r}, a field a run would fill")
    print(f"wrote {path}, sha256 {sha256_text(ROOT / path)}")
    print(f"fields a run would fill: none of {sorted(FORBIDDEN)} is in the written file")
    print(f"the one exception, checked by value: result == {written['result']!r}, the registry's name stamp")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
