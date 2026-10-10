# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register lane-magcal before it is run (data/results/magcal_registration.json).

    uv run python scripts/magcal_register.py

Every binding word is imported from `genomeos.attribution.magcal` (its BOUND constants), and
scripts/magcal_run.py refuses to run if any registered field differs from that module or any pinned
file's sha256 differs from the file on disk.

This script computes no scale, no split count, no held-out value and no outcome. It reads 7d's and
7h's committed results only for their gene and pair counts, their committed uncalibrated Delta MSE
(the number the run's harness must reproduce) and their registered readings, which are quoted
verbatim as context. No request is sent, no key is read.
"""

from __future__ import annotations

import json
import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import magcal as mc  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from magcal_run import C2_RESULT, COMMITTED, MR_RESULT, OWN_CODE, SOURCE_FILES  # noqa: E402

RESULT = "magcal_registration"
WILL_WRITE = "data/results/magcal_result.json"

EXPLORATION_PRECEDED_THIS = (
    "stated so that nothing here is read as a blind prediction. Before writing it this lane read:"
    " genomeos/attribution/combine2.py, multrna.py, capsumrna.py and saturate.py and their register and"
    " run scripts; scripts/saturable_gate.py for the split convention; ROADMAP §5 items 7, 7b, 7c, 7d and"
    " 7h and the stop list of 2026-09-27; and, from 7d's and 7h's committed results, each form's"
    " registered reading, outcome, gene and pair counts, Delta MSE, interval, control p, baseline MSE and"
    " form MSE, and the COLUMN LAYOUT and FIRST ROW of each genes table. It therefore knows the three"
    " committed uncalibrated Delta MSE values (mean -0.4435, competition -0.8102, multiplicative"
    " -0.5361), which is the harness it must reproduce, and it has seen two genes' P and D rows (ADAR,"
    " AGO3, A3GALT2). It has computed NO scale, NO split, NO training value and NO held-out value of any"
    " kind: this script reads the committed tables' lengths and the committed summary numbers, and the"
    " first held-out number of this design will be computed by scripts/magcal_run.py after this"
    " registration is committed. PRIOR EXPOSURE, which is larger than this paragraph, is stated"
    " separately and binds the reading"
)
DECISIONS = (
    "design: 7d's and 7h's, with one scale per form added and nothing else changed. The endpoint (ENCODE"
    " polyA plus RNA-seq of the thirteen reader biosamples, log2(TPM+1), centred within gene across the"
    " thirteen cells), the no-rule baseline, Delta MSE, the gene bootstrap, MIN_GENES and the reading"
    " vocabulary are theirs, carried word for word; the populations are theirs, gene for gene",
    "one number per form, and no more: the calibration is a single real scale a_F, fitted in closed form"
    " on the training pairs. There is no per-gene and no per-cell free parameter, no intercept and no"
    " transformation of P. a_F = 1 is the uncalibrated form 7d and 7h read",
    "reuse, not re-run: 7d's and 7h's committed genes tables carry the very P and D their registered"
    " readings were computed from, so the calibration is read off them. That is why the harness"
    " reproduces each form's committed Delta MSE before anything is split: it is the check that this lane"
    " reads the same population, endpoint and baseline. No compiled program, narrowPeak file or RNA"
    " quantification is opened by this lane, and none of 7d's, 7h's or capsumrna's files is edited",
    "split: lane-saturable's convention, taken unchanged (sha256 of the gene symbol, odd held out), on"
    " the SYMBOL so that a gene falls the same side in all three forms, with all thirteen of its pairs",
    "intervals: 7d Bonferroni-adjusted two forms from 95% to 97.5%; three forms are judged here, so the"
    " per-form interval is 98.33% on the same convention. The held-out interval resamples held-out genes"
    " with the scale fixed; the scale interval resamples training genes and refits",
    "control: a within-gene cell permutation of this lane's own, because 7d's and 7h's shuffled-reader"
    " control permutes element-level open/held rows that the committed genes tables do not carry, and"
    " re-deriving them would be a re-run of 7d and 7h rather than a calibration. It is named differently"
    " for that reason, their registered control p values are carried as context, and the control changes"
    " no verdict",
    "the decision procedure is the supervisor's science ruling of 2026-10-10 (magcal.WHICH_GOVERNS),"
    " carried verbatim: the held-out interval decides, and the scale interval only splits the straddle"
    " case. magcal.NOT_VACUOUS states in advance why that is not vacuous: a dead scale can manufacture a"
    " straddle but it cannot manufacture a win",
    "nothing in area J, no node model, no BioForge calibration, no new biosample (the thirteen reader"
    " cells only, READER_TERMS unchanged), no deletion request and 0 AlphaGenome requests: the stop list"
    " of 2026-09-27 as ROADMAP §5 reinstated it. No stop item is touched by this item",
)


def committed_facts() -> dict[str, Any]:
    """7d's and 7h's committed numbers the run's harness must reproduce, and their readings verbatim."""
    c2res = json.loads(C2_RESULT.read_text())
    mres = json.loads(MR_RESULT.read_text())
    return {
        "harness_targets": {f: dict(v) for f, v in COMMITTED.items()},
        "seven_d_verbatim": {
            "source": "data/results/combine2_result.json",
            "reading": c2res["reading"],
            "outcomes": c2res["outcomes"],
            "control_p": {f: c2res["per_form"][f]["control_p"] for f in ("mean", "competition")},
        },
        "seven_h_verbatim": {
            "source": "data/results/multrna_result.json",
            "reading": mres["reading"],
            "outcome": mres["outcome"],
            "control_p": mres["multiplicative"]["control_p"],
        },
    }


def inputs() -> list[dict[str, Any]]:
    return [
        mf.input_entry(
            p,
            partition=(
                "the source of one or two forms' P and D: pinned by sha256; its committed summary numbers"
                " and registered reading read, no scale, split or held-out value computed"
                if p in (C2_RESULT, MR_RESULT)
                else "pinned by sha256; the design it carries is quoted, not re-decided"
            ),
        )
        for p in SOURCE_FILES
    ]


def payload(entries: list[dict[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {
        "lane": mc.LANE,
        "date": date.today().isoformat(),
        "roadmap_item": (
            "§5 item M2h: the magnitude calibration of the compiled rules against measured RNA, registered"
        ),
        **mc.registered_fields(),
    }
    out.update(
        {
            "committed_facts": committed_facts(),
            "exploration_preceded_this": EXPLORATION_PRECEDED_THIS,
            "decisions": list(DECISIONS),
            "pinned_sha256": {x["path"]: x["sha256"] for x in entries},
            "what_the_run_will_write": WILL_WRITE,
            "the_code_that_runs_it": (
                "genomeos/attribution/magcal.py and scripts/magcal_run.py, committed before this"
                " registration, with the committed combine2.py and multrna.py they import for the words"
                " and constants carried verbatim"
            ),
            "alphagenome_requests": 0,
            "money": "none: no request is sent and no key is read",
        }
    )
    return out


def manifest(entries: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "sources": [
            {"accession": "data/results/combine2_result.json (7d, lane-combine2)", "version": "pinned"},
            {"accession": "data/results/multrna_result.json (7h, lane-multrna)", "version": "pinned"},
            {
                "accession": "ENCODE polyA plus RNA-seq gene quantifications, V29, 13 reader biosamples",
                "version": "read by 7d and 7h; this lane reads their committed rows, not the files",
            },
        ],
        "inputs": entries,
        "input_count": len(entries),
        "assembly": "GRCh38; genes as 7d and 7h matched them, by gene_id",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "forms": list(mc.FORMS),
            "split": mc.SPLIT,
            "min_genes": mc.MIN_GENES,
            "form_interval": list(mc.FORM_INTERVAL),
            "bootstrap": mc.BOOTSTRAP,
            "permutations": mc.PERMUTATIONS,
            "control_alpha": mc.CONTROL_ALPHA,
            "tolerance": mc.TOLERANCE,
            "free_parameters": "one real scale per form; no per-gene and no per-cell parameter",
        },
        "exclusions": ["none of this lane's own: each form's population is 7d's or 7h's, unchanged"],
        "partitions": "n/a: a registration; no split computed and no outcome read",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    if (RESULTS_DIR / "magcal_result.json").exists():
        raise SystemExit("refused: data/results/magcal_result.json already exists; the run is not re-run")
    entries = inputs()
    print("saved", save_result(RESULT, payload(entries), manifest=manifest(entries)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
