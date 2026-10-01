# SPDX-License-Identifier: AGPL-3.0-or-later
"""Step 1 of a second cell type for the CRISPRi deletion result: a blinded eligibility count.

Writes data/results/cell2_eligibility.json.

    uv run --frozen python scripts/cell2_eligibility.py

Counts, per held-out cell type, the pairs, the measured positives and the independent loci among
those positives, by the convention registered in genomeos/attribution/cell2.py and committed before
this script was first run. Blinded: it reads the benchmark's labels, coordinates, measured genes and
cell types through the project's own loader and nothing else. No prediction, no deletion value, no
model score, no AUPRC; no gain is computed, no model is fitted and no model request is made.

The decision the count serves is whether any second cell type can support a stratified test at all.
Pooled independent loci below cell2.POOLED_LOCUS_FLOOR is a no-go, recorded and stopped there. The
request counts a frozen deletion feature would need per cell type are quoted from the committed
data/results/crispri_published.json (second_cell_type.cost), not re-derived here.
"""

from __future__ import annotations

import ast
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import cell2, crispri  # noqa: E402
from genomeos.attribution.measured import CRISPRI_SPLIT_OF  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "cell2_eligibility"
PUBLISHED = Path("data/results/crispri_published.json")
#: The files this lane wrote, so the stamp can say which uncommitted code is somebody else's. Several
#: sessions work in this one checkout, so a stamp may be dirty through no act of this lane; a referee
#: can accept a named foreign file that the counting path does not import, but not a bare dirty flag.
OWN_CODE = (
    "genomeos/attribution/cell2.py",
    "scripts/cell2_eligibility.py",
    "tests/test_cell2_eligibility.py",
)
#: The repository's own package, for computing the counting path rather than listing it by hand.
PACKAGE = "genomeos"
ROOT = Path(__file__).resolve().parents[1]


def counting_path(entry: Path | None = None) -> list[str]:
    """Every module of this repository that the entry script can reach by import, computed as the
    transitive closure of its import statements. A hand-written list cannot be checked and goes stale;
    this is read from the syntax of the files themselves at write time."""
    entry = entry or Path(__file__).resolve()
    seen: dict[str, Path] = {}
    queue = [entry]
    while queue:
        path = queue.pop()
        rel = str(path.relative_to(ROOT))
        if rel in seen:
            continue
        seen[rel] = path
        try:
            tree = ast.parse(path.read_text())
        except (OSError, SyntaxError):
            continue
        names: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names.update(a.name for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module and not node.level:
                names.add(node.module)
                names.update(f"{node.module}.{a.name}" for a in node.names)
        for name in names:
            if not (name == PACKAGE or name.startswith(PACKAGE + ".")):
                continue
            parts = name.split(".")
            for candidate in (
                ROOT.joinpath(*parts).with_suffix(".py"),
                ROOT.joinpath(*parts, "__init__.py"),
            ):
                if candidate.is_file():
                    queue.append(candidate)
    return sorted(seen)


def code_cleanliness() -> dict:
    """Which uncommitted code the stamp names, split into this lane's and other lanes', and whether any
    of it is on the counting path. Read from git and from the imports at write time, not asserted."""
    rev = mf.code_revision()
    path = counting_path()
    dirty = list(rev["dirty_code_paths"])
    own = [p for p in dirty if p in OWN_CODE]
    foreign = [p for p in dirty if p not in OWN_CODE]
    return {
        "git_sha": rev["git_sha"],
        "dirty": rev["dirty"],
        "own_uncommitted_code": own,
        "own_code_is_committed": not own,
        "foreign_uncommitted_code": foreign,
        "foreign_uncommitted_code_on_the_counting_path": [p for p in foreign if p in path],
        "counting_path": path,
        "counting_path_is_computed": (
            "the transitive import closure of this script over the repository's own package, computed "
            "from the files' import statements at write time, not a hand-written list"
        ),
        "note": (
            "several sessions work in this one checkout. A file listed under foreign_uncommitted_code "
            "belongs to another lane; this lane did not write it and did not commit it. The counting "
            "path is the computed closure above, so a foreign file outside it cannot have entered the "
            "count, and foreign_uncommitted_code_on_the_counting_path names any that could"
        ),
    }


def quoted_cost() -> dict:
    """The per-cell-type request counts already in the committed CRISPRi result, quoted as they stand."""
    payload = json.loads(PUBLISHED.read_text())
    cost = payload["second_cell_type"]["cost"]
    return {
        "source": f"{PUBLISHED} second_cell_type.cost, quoted not re-derived",
        "per_cell_type": cost,
        "budget_in_that_result": payload["second_cell_type"]["budget"],
        "why_that_result_did_not_run": payload["second_cell_type"]["why"],
        "total_requests_for_all_three": sum(v["requests_frozen_feature"] for v in cost.values()),
    }


def registration_draft(counted: dict, verdict: dict, cost: dict, exposed: dict, power_gate: dict) -> dict:
    """The registration this lane would propose, drafted and not run. The request count goes to Albert
    for authorisation; this lane makes no request."""
    candidates = counted["candidate_second_cell_types"]
    with_feature = [c for c in candidates if counted["per_cell_type"][c]["deletion_value_available"]]
    return {
        "status": "draft only, not registered and not run; the coordinator reviews before anything else",
        "question": (
            "does the deletion feature's contribution to the CRISPRi element-to-gene benchmark hold "
            "outside K562, in a second held-out cell type measured on its own positives"
        ),
        "stratum": (
            "cell type: each candidate second cell type is its own stratum, scored and reported alone, "
            "never pooled with K562 into one number"
        ),
        "unit_of_independence": counted["convention"],
        "unit_of_independence_is_operational": counted["convention_is_operational"],
        "preset_floor": {
            "pooled_independent_loci": counted["floor"],
            "verdict_of_step_1": verdict,
            "note": "the floor was committed in genomeos/attribution/cell2.py before this count was taken",
        },
        "cell_types_that_could_carry_a_deletion_gain_at_all": with_feature,
        "cell_types_that_could_not": [c for c in candidates if c not in with_feature],
        "loci_where_a_gain_could_be_measured_at_all": counted["pooled_over_candidates_with_a_deletion_value"],
        "the_constraint_this_puts_on_the_draft": (
            "the pooled count the floor is read against is carried mostly by cell types that have no "
            "deletion value, where no deletion gain can be measured. Only the cell types listed under "
            "cell_types_that_could_carry_a_deletion_gain_at_all can be a stratum of this test as the "
            "feature stands, and their own locus count is the one that limits it"
        ),
        "cost_table_does_not_cover_every_candidate": (
            "the quoted second_cell_type.cost covers HCT116, WTC11 and Jurkat only. It carries no figure "
            "for a candidate outside those three, so no request count for such a candidate is stated "
            "here: this lane quotes the committed figures and derives none"
        ),
        "availability_rule": (
            "a cell type outside crispri.MODEL_CELLS has no deletion value, so no deletion gain is "
            "measured there: crispri.gain_where_available refuses a number rather than report the "
            "difference between two model forms as a measurement (crispri.UNAVAILABLE_GAIN)"
        ),
        "falsifier": (
            "the claim is that the deletion feature adds to the frozen baseline in a second cell type's "
            "own positives. It is falsified if, in every candidate stratum that has enough independent "
            "loci, the measured gain's 95% interval over loci includes zero or lies below it. A gain "
            "that holds only when the strata are pooled with K562, or only in a stratum whose largest "
            "locus holds most of its positives, does not count as support"
        ),
        "outcomes": {
            "interval_above_zero_in_a_stratum": (
                "the feature's contribution is reproduced outside K562 in that cell type's own "
                "positives, on the loci named in this count"
            ),
            "interval_includes_zero": (
                "no contribution is demonstrated in that cell type at this number of loci; the stratum "
                "is reported as uninformative, not as evidence of no effect"
            ),
            "interval_below_zero": (
                "the feature costs performance outside K562, which contradicts the single-cell-type "
                "reading and must be reported as such"
            ),
        },
        "model_requests_needed_per_cell_type": cost,
        "prior_exposure": exposed,
        "resampling_unit": {
            "primary": (
                "a cluster bootstrap over the independent loci registered here, resampling whole loci "
                "with replacement, with cell type as a stratum so each stratum is resampled within "
                "itself and never pooled with K562"
            ),
            "also_reported": (
                "the chromosome-cluster interval, resampling whole chromosomes, which is the estimator "
                "the K562 claim used (crispri.gain_interval)"
            ),
            "which_is_quoted": (
                "the wider of the two intervals is the one quoted; the narrower one is reported beside "
                "it and never in its place"
            ),
            "why": (
                "the locus is the unit the eligibility count treats as independent, so it is the unit "
                "the interval must resample; the chromosome interval is kept for comparability with "
                "the published K562 figure, and quoting the wider of the two cannot flatter either"
            ),
        },
        "power_gate": power_gate,
        "readings_fixed_in_advance": cell2.READINGS,
        "cost_measured_not_estimated": {
            "source": (
                "the comparable HCT116 fetch of 705 requests on 2026-09-28, "
                "data/jobs/crispri_hct116_fetch.log, 10:23:23 to 10:26:38, refused_quota: 0"
            ),
            "rate": "3.6 requests per second over the client's five workers",
            "fetching_for_725_requests": "about 3.3 minutes",
            "scoring": "about a minute",
            "cache": "0.57 MB, at 793 bytes per answer measured over 135,829 cached answers",
            "money": (
                "none: AlphaGenome is free under non-commercial terms, with an unpublished daily quota"
            ),
            "local_compute": "negligible",
            "quoted_as": "measured, with the source named; this lane re-derived none of it",
        },
        "cost_that_is_not_a_resource": {
            "what": (
                "WTC11 and Jurkat are the only cell types in this benchmark never scored here, holding "
                "15 positives in 9 loci and 7 in 4"
            ),
            "why_it_is_a_cost": (
                "sending the requests converts them permanently into examined data, so this is the "
                "project's one unexposed look, and the power gate decides whether it is worth it"
            ),
        },
        "authorisation": {
            "albert_approved_wording": (
                "I approve the 725 AlphaGenome requests for the second cell type (WTC11 614, Jurkat "
                "111), sent only after the registration is committed with genomeos-3b's four terms and "
                "the power figure is at least 0.5, and only with genomeos-3b's written sign-off."
            ),
            "requests_approved": 725,
            "per_cell_type": {"WTC11": 614, "Jurkat": 111},
            "cached_values_reused_without_a_request": ["GM12878", "HCT116"],
            "binds_before_a_single_request_is_sent": [
                "the registration is committed, carrying all four terms",
                "the power figure is at least 0.5, against the line registered before it is computed",
                "the supervisor's written sign-off, which the coordinator obtains",
            ],
            "below_the_line": (
                "below 0.5 nothing is spent and the reading goes on record as 'not powered to decide', "
                "which is a finished lane and a useful answer, not a failure. The subsampling, the "
                "structure, the effect size and the line are not adjusted to lift the figure over it"
            ),
            "this_lane_sent_nothing": (
                "this lane made no request, holds none, and sends none; it did not obtain the sign-off "
                "and did not run the comparison"
            ),
        },
        "what_this_lane_did_not_do": [
            "computed no gain and fitted no model",
            "read no prediction, deletion value, model score or AUPRC",
            "made no model request and sent none",
            "did not compute the power share, for the two reasons recorded under power_gate",
        ],
    }


def power_gate(heldout: list[crispri.Pair], counted: dict) -> dict:
    """The gate's registered design, and the one half of it that can be settled blind: whether the
    primary cell type's held-out positives can supply the second-cell structure at all.

    The share itself is not computed here. Two reasons, both recorded: the exact profile the gate
    specifies does not exist in the primary cell type's held-out positives, and computing a share of
    intervals at an observed effect requires reading measured deletion values, which this lane's
    blinding forbids."""
    pooled = counted["pooled_over_candidates"]
    shapes = pooled["locus_shapes_genome_wide"]
    primary = counted["primary_cell_type"]
    primary_keys = [
        cell2.locus_key(p)._replace(cell="") for p in heldout if p.regulated and p.cell == primary
    ]
    match = cell2.profile_match(shapes["positives_per_locus_counts"], primary_keys)
    return {
        "design": cell2.POWER_GATE_DESIGN,
        "share_floor": cell2.POWER_GATE_SHARE_FLOOR,
        "share_floor_registered_before_the_computation": True,
        "structure_to_match": {
            "population": pooled["population"],
            "independent_loci": pooled["independent_loci_genome_wide"],
            "positives": pooled["positives"],
            "positives_per_locus_counts": shapes["positives_per_locus_counts"],
        },
        "source_population": (
            f"the measured positives of the {len([p for p in heldout if p.cell == primary])} held-out "
            f"{primary} pairs, grouped by the same locus convention"
        ),
        "profile_feasibility": match,
        "share": None,
        "share_not_computed_because": [
            (
                "the structure the gate specifies cannot be drawn from the primary cell type's "
                "held-out positives: they hold no locus of 6, 7, 8 or 10 positives, which the "
                "second-cell profile requires, so no subsample reproduces it. The shortfall is given "
                "under profile_feasibility, computed from labels alone"
            ),
            (
                "a share of intervals at the primary cell type's observed effect cannot be computed "
                "without reading measured deletion values and a gain, which this lane is blinded "
                "against: its stop conditions are labels, coordinates, genes and cells only, and to "
                "compute no gain and fit no model"
            ),
        ],
        "what_this_means_for_the_approved_requests": (
            "Albert's approval requires the power figure to be at least 0.5. No figure exists yet, so "
            "the condition is not met and nothing may be spent. This is not a figure below the line "
            "and must not be read as one: the gate as registered is not computable on this source, and "
            "whether to relax the profile, or to assign the share to a lane that is not blinded, is a "
            "decision for the coordinator and the supervisor, taken before the figure is computed"
        ),
        "not_adjusted": (
            "the subsampling, the structure, the effect size and the line were not adjusted; the "
            "shortfall is reported as found"
        ),
    }


def manifest(heldout: list[crispri.Pair]) -> dict:
    return {
        "sources": [
            {
                "accession": "EngreitzLab/CRISPR_comparison resources/crispr_data, EPCrisprBenchmark "
                "heldout_5_cell_types (Gschwind et al.)",
                "version": "main branch, unpinned upstream; fetched 2026-09-16; pinned here by sha256",
                "url": crispri.BASE_URL,
            },
            {
                "accession": "data/results/crispri_published.json, second_cell_type.cost",
                "version": "the committed result of 2026-09-27; quoted, not re-derived",
            },
        ],
        "inputs": [
            mf.input_entry(
                crispri.KNOWLEDGE / crispri.HELDOUT,
                partition=CRISPRI_SPLIT_OF[crispri.HELDOUT],
                pairs=len(heldout),
                role="labels, coordinates, measured genes and cell types only",
            ),
            mf.input_entry(
                PUBLISHED,
                partition=None,
                role="second_cell_type.cost: the request counts quoted in the draft",
            ),
        ],
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "independent_locus_span_bp": cell2.INDEPENDENT_LOCUS_SPAN,
            "pooled_locus_floor": cell2.POOLED_LOCUS_FLOOR,
            "primary_cell": cell2.PRIMARY_CELL,
            "blinded_fields": list(cell2.BLINDED_FIELDS),
            "model_cells": list(crispri.MODEL_CELLS),
            "alphagenome_requests": 0,
            "performance_metrics_computed": 0,
            "gains_computed": 0,
        },
        "exclusions": [
            "pairs the benchmark marks ValidConnection=FALSE are dropped by crispri.parse, as in every "
            "other use of this table",
            "K562 is not a candidate second cell type: the committed result already covers it; its "
            "counts are reported for comparison only",
            "no training pair is read: the question is about the held-out cell types",
        ],
        "partitions": {
            "heldout": "the 4,378 valid held-out pairs across five cell types; the only pairs read here"
        },
        "code_cleanliness": code_cleanliness(),
    }


def main() -> int:
    heldout = crispri.load(crispri.HELDOUT)
    counted = cell2.eligibility(heldout, crispri.MODEL_CELLS)
    verdict = cell2.verdict(counted)
    cost = quoted_cost()
    exposed = cell2.exposure(heldout)
    gate = power_gate(heldout, counted)

    payload = {
        "status": (
            "step 1 only: a blinded eligibility count of independent loci among the measured positives "
            "of each held-out cell type. No gain, no model, no request"
        ),
        "lane": "lane-cell2",
        "date": "2026-10-01",
        "question": (
            "can any second held-out cell type support a stratified test of the CRISPRi deletion "
            "result, judged by independent loci among its measured positives rather than by pair counts"
        ),
        "blinding": (
            "labels, coordinates, measured genes and cell types only, through genomeos.attribution."
            "crispri.load; no prediction, deletion value, model score or AUPRC was read or computed"
        ),
        "code_cleanliness": code_cleanliness(),
        "counts": counted,
        "prior_exposure": exposed,
        "power_gate": gate,
        "verdict": verdict,
        "cost_quoted": cost,
        "alphagenome_requests": 0,
    }
    if verdict["meets_floor"]:
        payload["registration_draft"] = registration_draft(counted, verdict, cost, exposed, gate)
    else:
        payload["no_go"] = {
            "decision": "no-go: the pooled independent loci are below the pre-set floor; the lane stops",
            "floor": verdict["floor"],
            "counted": verdict["pooled_independent_loci_genome_wide"],
            "note": "no registration is drafted, no gain is computed and no model request is made",
        }
    payload[mf.KEY] = manifest(heldout)
    path = save_result(RESULT, payload)

    for cell, row in counted["per_cell_type"].items():
        shapes = row["locus_shapes"]
        print(
            f"{cell:8s} pairs={row['pairs']:5d} positives={row['positives']:4d} "
            f"loci={row['independent_loci']:4d} largest={shapes['largest_locus_positives']} "
            f"({shapes['largest_locus_share_of_positives']}) deletion={row['deletion_value_available']}"
            f"{'' if row['candidate_second_cell_type'] else '  [primary, not a candidate]'}"
        )
    pooled = counted["pooled_over_candidates"]
    print(
        f"pooled over {', '.join(counted['candidate_second_cell_types'])}: "
        f"pairs={pooled['pairs']} positives={pooled['positives']} "
        f"loci_genome_wide={pooled['independent_loci_genome_wide']} "
        f"loci_summed={pooled['independent_loci_summed_per_cell_type']}"
    )
    subset = counted["pooled_over_candidates_with_a_deletion_value"]
    print(
        f"of those, in candidates with a deletion value ({', '.join(subset['cell_types']) or 'none'}): "
        f"positives={subset['positives']} loci_genome_wide={subset['independent_loci_genome_wide']} "
        "(reported beside the floor, not read against it)"
    )
    print(f"verdict: {verdict['decision']} (floor {verdict['floor']})")
    print(
        f"exposure under the pooled convention: {exposed['already_exposed_loci']} of "
        f"{exposed['pooled_loci']} loci already scored "
        f"({exposed['already_exposed_positives']} of {exposed['pooled_positives']} positives), "
        f"{exposed['never_exposed_loci']} never scored; "
        f"{exposed['loci_in_more_than_one_cell_type']} loci appear in more than one cell type"
    )
    print(
        f"power gate: share floor {gate['share_floor']}, share not computed; exact profile available "
        f"in {counted['primary_cell_type']}: {gate['profile_feasibility']['exact_profile_available']}"
    )
    print(f"-> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
