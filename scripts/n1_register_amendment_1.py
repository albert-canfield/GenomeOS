# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register amendment 1 to N1 before any measurement (data/results/n1_registration_amendment_1.json).

    uv run --frozen python scripts/n1_register_amendment_1.py

The original registration (data/results/n1_registration.json, f9a9130) and its code (b359867) stay as
they were; this records five changes the owner's reviewer required, each with what it replaces and why,
and freezes the sha256 of the code scripts/n1_run_v2.py checks before opening any measurement. It reads
no data file. docs/ATTRIBUTION.md, "N1 ... registered before any measurement", subsection "Amendment 1".
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import n1_perturb_response as n1  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

ORIGINAL = RESULTS_DIR / "n1_registration.json"
CODE = {
    "module": Path("genomeos/attribution/n1_perturb_response.py"),
    "runner": Path("scripts/n1_run_v2.py"),
    "original_runner": Path("scripts/n1_run.py"),
    "tests": Path("tests/test_n1_perturb_response.py"),
    "register_amendment": Path("scripts/n1_register_amendment_1.py"),
}
PMC = "https://europepmc.org/article/PMC/PMC9380471"
FIGSHARE = "https://doi.org/10.25452/figshare.plus.20029387.v1"
PGI = "https://github.com/thomasmaxwellnorman/Perturbseq_GI/blob/3b25109aeb9c0c2026bd70abd50304a0ad4e5395"


def last_commit(path: Path) -> str | None:
    out = subprocess.run(["git", "log", "-1", "--format=%h", "--", str(path)], capture_output=True, text=True)
    return out.stdout.strip() or None


def knockdown_evidence() -> list[dict]:
    return [
        {
            "source": "Replogle et al. 2022, STAR Methods, 'Leverage scores for quantifying perturbation "
            "penetrance and variability'",
            "link": PMC,
            "quote": "Knockdown was computed as the ratio of mean (unnormalized) expression of the target "
            "gene "
            "within perturbed cells vs. that in cells with non-targeting sgRNAs.",
            "establishes": "knockdown as a ratio of mean unnormalized expression, perturbed over "
            "non-targeting",
            "does_not_establish": "that fold_expr or control_expr hold these quantities; whether "
            "'unnormalized' "
            "means raw or depth-adjusted counts; all or core controls; pooled or per gemgroup; filtered or "
            "unfiltered cells",
        },
        {
            "source": "Replogle et al. 2022, STAR Methods, 'Filtering and internal normalization of gene "
            "expression measurements'",
            "link": PMC,
            "quote": "We first computed scale factors to adjust for variable sequencing depths across "
            "gemgroups: "
            "we examined all core control cells (which make up ~4% of all cells), computed factors that "
            "equalized the mean UMI counts within these cells across gemgroups, and then applied these "
            "factors "
            "to all cells in the gemgroup to produce adjusted UMI counts.",
            "establishes": "depth-adjusted (non-integer) UMI counts exist in the pipeline, besides raw "
            "counts "
            "from reads downsampled per gemgroup",
            "does_not_establish": "which of them the knockdown fields use",
        },
        {
            "source": "Replogle et al. 2022, Figure S3 legend",
            "link": PMC,
            "quote": "The fractional change in expression is defined as the expression in the targeted cells "
            "minus the expression in non-targeting cells, relative to the expression in the non-targeting "
            "cell "
            "population (-1 implies 100% knockdown).",
            "establishes": "a fractional change, (perturbed - control) / control, used in the paper's "
            "figures",
            "does_not_establish": "that pct_expr is this quantity; pct_expr = fold_expr - 1 would be "
            "algebraic "
            "consistency only, not meaning or units",
        },
        {
            "source": "Figshare+ 20029387 description",
            "link": FIGSHARE,
            "quote": "In the anndata format, the .var annotation details genes while the .obs annotation "
            "details "
            "single-cells/pseudobulk populations.",
            "establishes": "nothing about individual .obs fields: no field is defined",
            "does_not_establish": "any field's meaning or units",
        },
        {
            "source": "Norman et al. 2019 producer code (the predecessor codebase the 2022 paper cites), "
            "GI_generate_populations.ipynb, cell 27",
            "link": f"{PGI}/GI_generate_populations.ipynb",
            "quote": "'control_first_expr': lambda meta, expr: expr.loc[meta['first_id'], 'mean'] ... "
            "mean_pop.cells['fold_first_expr'] = "
            "mean_pop.cells['first_expr']/mean_pop.cells['control_first_expr']",
            "establishes": "in the 2019 CRISPRa data, a control expression taken as the control population's "
            "per-gene 'mean' and a fold as the perturbed mean over it, with fields named in the same pattern",
            "does_not_establish": "the 2022 computation: that code is not public, the 2022 pipeline added "
            "per-gemgroup depth adjustment, and pct_expr has no 2019 counterpart",
        },
        {
            "source": "Norman et al. 2019 producer code, perturbseq/cell_population.py",
            "link": f"{PGI}/perturbseq/cell_population.py",
            "quote": "def metaapply(self, function_dict, axis=1, normalized=False, **kwargs): if not "
            "normalized: "
            "df = self.matrix ... gene_list['mean'] = matrix.mean()",
            "establishes": "in 2019, both the perturbed and the control means come from the unnormalized "
            "matrix, pooled over all control cells",
            "does_not_establish": "the same for the 2022 file",
        },
    ]


def changes() -> list[dict]:
    return [
        {
            "item": 1,
            "title": "AUROC of a constant predictor",
            "was": "auroc() returned (None, 'no_score_variation') when an arm's scores did not vary, so "
            "such a "
            "factor left the paired estimate as undefined",
            "now": "auroc_v2: the Mann-Whitney probability that a responder outscores a non-responder, ties "
            "counting half; a constant score gives 0.5; undefined only when a class is empty",
            "code": "auroc_v2, used by stratified_auroc and factor_result_v2",
            "tests": [
                "test_auroc_v2_counts_ties_half_and_a_constant_score_gives_half",
                "test_the_analysis_after_the_stop_once_a_documented_amendment_supports_it",
            ],
        },
        {
            "item": 2,
            "title": "Knockdown units",
            "finding": "unsupported",
            "was": "eligibility required num_cells_filtered * control_expr >= 10, read as expected target "
            "UMIs "
            "with a Poisson argument, and checked pct_expr = fold_expr - 1 as the fields' meaning",
            "now": "the knockdown rule is marked unsupported: no producer documentation or public 2022 code "
            "defines control_expr, fold_expr or pct_expr, or the units of control_expr. "
            "num_cells_filtered * control_expr is not established as an expected UMI count, so the Poisson "
            "reading is withdrawn, and pct_expr = fold_expr - 1 would show algebraic consistency only. "
            "ELIGIBILITY_SUPPORTED = False and run_v2 stops before reading any candidate row. No threshold "
            "replaces the rule.",
            "evidence": knockdown_evidence(),
            "what_would_resolve_it": [
                "the 2022 producer code, or a statement by the authors defining the three fields and their "
                "units",
                "a separately reviewed amendment deriving knockdown from documented quantities (for example "
                "the "
                "raw pseudobulk's target column against its non-targeting rows), with a threshold that "
                "documentation supports",
            ],
            "tests": ["test_run_v2_stops_before_the_knockdown_rule"],
        },
        {
            "item": 3,
            "title": "The gate's claim, narrowed",
            "framing": "N1 is an exploratory comparison of two prediction arms against an operational "
            "response "
            "label, |T| >= 3. T = X * sqrt(num_cells_filtered) is a proposed statistic, not a calibrated "
            "test. "
            "Passing the gate does not validate per-gene biological responses.",
            "what_the_gate_cannot_establish": [
                "that cells are independent",
                "the uncertainty in the control mean and sd that define each z-score",
                "calibration across cell counts: the gate sees only the non-targeting rows' own counts",
                "calibration of pooled rows",
                "the variance of perturbed cells",
            ],
            "correction": "The registration of f9a9130 gave as a reason for using all 585 non-targeting rows "
            "that the 514 core controls define the normalization's control mean and sd and would flatter the "
            "null. The 585 include those 514, so using them does not avoid that role; it only adds the 71 "
            "guides outside the core set.",
        },
        {
            "item": 4,
            "title": "Expression matching",
            "departure": "The proposal (docs/ROADMAP.md, N1 row of 2026-10-01) listed control-based "
            "expression "
            "matching among what step 2 freezes. The registration of f9a9130 used control-expression deciles "
            "only in the gate, not in the AUROC comparison, so the matching was not implemented there.",
            "now": "the primary per-factor AUROC compares responders only with non-responders of the same "
            "control-expression decile (the gate's deciles over the calibrated universe genes), each stratum "
            "weighted by its responder x non-responder pairs (stratified_auroc). The unstratified auroc_v2 "
            "is "
            "reported beside it with no criterion. The estimate, floor, bootstrap and the two criteria are "
            "those of f9a9130, applied to the stratified differences.",
            "runs_under_this_amendment": False,
            "why_not_run": "run_v2 stops before eligibility (item 2); the comparison is frozen and tested so "
            "that a documented amendment changes one constant, not the analysis",
            "tests": [
                "test_stratified_auroc_compares_within_expression_strata",
                "test_the_analysis_after_the_stop_once_a_documented_amendment_supports_it",
            ],
        },
        {
            "item": 5,
            "title": "Execution matches the freeze",
            "now": [
                "scripts/n1_run_v2.py checks, before opening any measurement, the sha256 of the original "
                "registration, of the analysis module and of itself, and CONSTANTS_V2, against this "
                "amendment, "
                "and refuses on any difference (freeze_problems)",
                "H5adPseudobulkV2 decodes X at the selected columns only, one row at a time (X[i, sorted "
                "columns]); the operating system may still read surrounding bytes from disk",
                "scripts/n1_run.py refuses with a pointer to scripts/n1_run_v2.py (one added statement)",
            ],
            "correction": "the registration of f9a9130 said the gate reads X 'universe genes only'; its "
            "reader "
            "(H5adPseudobulk) decoded whole rows of 8,248 genes and kept the 860, so that description was "
            "inaccurate for it",
            "tests": [
                "test_freeze_problems_refuse_any_difference",
                "test_the_v2_run_script_checks_the_freeze_first_then_runs_once",
                "test_reader_v2_decodes_only_the_selected_columns",
                "test_the_original_run_script_refuses_and_points_to_v2",
            ],
        },
    ]


def main() -> int:
    original = json.loads(ORIGINAL.read_text())
    frozen_code = {
        k: {"path": str(p), "sha256": n1.sha256_file(p), "last_commit": last_commit(p)}
        for k, p in CODE.items()
    }
    n_universe = len(original["plan"]["universe"])
    payload = {
        "status": "amendment registered before any measurement; not run",
        "lane": "lane-n1",
        "amends": {
            "result": "n1_registration",
            "commit": last_commit(ORIGINAL),
            "code_commit": original["result_manifest"]["code"]["git_sha"][:7],
            "sha256": n1.sha256_file(ORIGINAL),
            "unchanged": "the original registration and its code stay as committed; this amendment adds",
        },
        "framing": changes()[2]["framing"],
        "changes": changes(),
        "unchanged": [
            "the shared universe (860 genes), the candidates (56 factors, 59 rows, 39 units) and the ledger",
            "both arms' scores; their positive counts differ by design (for example CGGBP1, 40 attribution "
            "and "
            "598 proximity genes) and are not equalised",
            "the operational label |T| >= 3, the gate's rows and pass rule, the multi-row rule, aliases, "
            "families, the floor of 30, the cluster bootstrap and the two separate criteria",
        ],
        "run_v2": {
            "needs": "the owner's separate authorisation, given to scripts/n1_run_v2.py --authorisation",
            "reads_in_order": [
                "the freeze check: no data file is opened",
                "both files' md5 against Figshare, and their sha256",
                "both files' identities: obs/gene_transcript and var/gene_id",
                f"the gate: normalized and raw X at the 585 non-targeting rows, decoded at the {n_universe} "
                "universe columns only, and obs/num_cells_filtered at those rows",
            ],
            "then": "stops with status eligibility_unsupported; no candidate row, knockdown field or "
            "response "
            "is read",
            "result": "data/results/n1_result_amendment_1.json, written once",
        },
        "frozen_code": frozen_code,
        "constants": json.loads(json.dumps(n1.CONSTANTS_V2)),
    }
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "GenomeOS n1_registration (the registration amended)",
                "version": payload["amends"]["commit"],
            },
            {
                "accession": "Replogle et al. 2022, Cell, STAR Methods and Figure S3 (documentation only)",
                "version": "doi 10.1016/j.cell.2022.05.013, PMC9380471",
            },
            {"accession": "Figshare+ 20029387 description (documentation only)", "version": "v1"},
            {"accession": "Perturbseq_GI (Norman et al. 2019) producer code", "version": "3b25109"},
        ],
        "inputs": [mf.input_entry(ORIGINAL, partition=None)]
        + [mf.input_entry(p, partition=None) for k, p in CODE.items() if k in ("module", "runner")],
        "assembly": "n/a: an amendment; no coordinate is read",
        "coordinates": "n/a: an amendment; no coordinate is read",
        "parameters": dict(n1.CONSTANTS_V2),
        "exclusions": [],
        "partitions": "n/a: an amendment registered before any measurement; nothing is evaluated",
    }
    p = save_result(n1.AMENDMENT, payload)
    print(f"wrote {p}: amends {payload['amends']['commit']} ({payload['amends']['sha256'][:12]})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
