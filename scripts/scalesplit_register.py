# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register lane-scalesplit before it is run (data/results/scalesplit_registration.json).

    uv run python scripts/scalesplit_register.py

Every binding word is imported from `genomeos.attribution.scalesplit` (its BOUND constants), and
scripts/scalesplit_run.py refuses to run if any registered field differs from that module or any
pinned file's sha256 differs from the file on disk.

This script computes NO chromosome-fold scale, NO chromosome-fold Delta MSE, NO held-out value and NO
outcome. It reads 7h's and M2h's committed results only for their counts, their committed numbers (the
numbers the harness must reproduce) and their registered readings, which are quoted verbatim as
context. No request is sent, no key is read.
"""

from __future__ import annotations

import json
import os
import sys
from datetime import date
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import scalesplit as ss  # noqa: E402
from genomeos.results import save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scalesplit_harness import MAGCAL_RESULT, MR_RESULT, OWN_CODE, SOURCE_FILES  # noqa: E402

RESULT = "scalesplit_registration"
WILL_WRITE = "data/results/scalesplit_result.json"
HARNESS_RESULT = "data/results/scalesplit_harness.json"

EXPLORATION_PRECEDED_THIS = (
    "stated so that nothing here is read as a blind prediction. Before writing it this lane read:"
    " genomeos/attribution/magcal.py in full and scripts/magcal_run.py and scripts/magcal_register.py;"
    " the committed M2h result's multiplicative block (its split counts, fitted scale, scale interval,"
    " held-out Delta MSE, held-out interval, control p, outcome and the first lines of its descriptives);"
    " the shape of data/results/multrna_result.json and 7h's committed gene count, pair count and"
    " uncalibrated Delta MSE; and the GENE COUNT PER CHROMOSOME of 7h's committed table (chr1 1,700 down"
    " to chrY 13), which is what told it that both parity sets clear MIN_GENES and that chrY alone falls"
    " below the floor as a single chromosome. It therefore knows, before registering: 7h's -0.5361, M2h's"
    " +0.2146 and +0.04339, and the 24 chromosomes' gene counts. It has computed NO scale on any"
    " chromosome set, NO Delta MSE on any chromosome subset and NO held-out value of this design; this"
    " lane's step 1 (data/results/scalesplit_harness.json, committed before this registration) states"
    " that in its own words and reproduced M2h's numbers only. The first chromosome-fold number of this"
    " design will be computed by scripts/scalesplit_run.py after this registration is committed. PRIOR"
    " EXPOSURE, which is larger than this paragraph, is stated separately and binds the reading"
)
DECISIONS = (
    "WHY THIS LANE EXISTS: M2h's one positive result is split by sha256 of the gene SYMBOL -- a RANDOM"
    " split. Genes that share a compiled element are not independent, which M2i's own post hoc"
    " acknowledged by reporting a 1-Mb locus-cluster interval beside its gene-level one, so a training"
    " gene and a held-out gene may be fitted through the same element and 'held out' is weaker than it"
    " sounds. This lane re-splits the same rows by CHROMOSOME, which cannot leak that way, and reads M2h's"
    " endpoint again",
    "design: M2h's, with ONE thing changed and nothing else. The endpoint, the no-rule baseline, Delta"
    " MSE, the closed-form fit, the gene bootstrap, the within-gene cell permutation control, MIN_GENES,"
    " the four-way decision table and every reading word are imported from genomeos/attribution/magcal.py"
    " and carried verbatim; the population is 7h's, gene for gene. The ONE change is the split",
    "one form, not three: the multiplicative composition, which is M2h's arm that read 'better'. The two"
    " arms that read 'worse' at their best scale are not re-split, because there is no positive result"
    " there to audit for leakage",
    "THE SCALE IS REFITTED on each fold's training chromosomes and is NOT taken as M2h's constant"
    " +0.2146. That is deliberate: refitting tests whether the SCALE TRANSPORTS between chromosome sets,"
    " a property of the compiled rules, rather than whether one number fitted on a random half happens to"
    " hold, which is M2h's own question and is already answered",
    "split: by chromosome, by the PARITY OF THE CHROMOSOME NAME, with TWO COMPLEMENTARY FOLDS both read"
    " so that neither side is a chosen one. Parity needs no gene count, no RNA value and no fitted number"
    " to state, and puts no chromosome in both sides. The sex chromosomes go with the even autosomes by"
    " the one rule that they are not odd-numbered -- 646 of 16,509 genes, declared rather than tuned",
    "interval: 98.33%, magcal.FORM_INTERVAL unchanged, decided now. It keeps the number directly"
    " comparable to M2h's committed +0.04339 [98.33% CI +0.0379, +0.04894]; it is what magcal's bootstrap"
    " functions compute, and this lane calls them read-only rather than editing magcal or reimplementing"
    " its bootstrap; and it is WIDER than the 97.5% a Bonferroni over this lane's two folds would require"
    " (7d's two-arm convention), so it cannot manufacture an exclusion of 0. It is labelled a deliberate"
    " over-correction and never a Bonferroni over two",
    "control: M2h's own within-gene cell permutation, applied per fold with the scale refitted on the"
    " permuted training chromosomes. It is reported in every fold sentence and changes no verdict",
    "the decision procedure is M2h's WHICH_GOVERNS, the supervisor's science ruling of 2026-10-10,"
    " carried verbatim and applied by magcal.decide: the held-out interval decides, and the scale interval"
    " only splits the straddle case -- a scale interval including 0 with a straddling held-out interval"
    " reads 'calibration removes the rules', A NEGATIVE. Above it, this lane's own TRANSPORT_RULE is fixed"
    " here: both folds BETTER -> SURVIVES_BOTH_FOLDS, one fold -> SURVIVES_ONE_FOLD, neither ->"
    " DOES_NOT_SURVIVE, and the transport sentence is read FIRST so that a failure to transport cannot be"
    " buried",
    "the magnitude comparison against M2h is registered as DECIDING NOTHING, with its rule fixed in"
    " advance: overlapping intervals read CONSISTENT WITH THE SAME EFFECT and neither 'stronger' nor"
    " 'weaker' may be said. The two populations are the same 16,509 genes repartitioned, so no test of the"
    " difference is registered and none will be computed",
    "the per-chromosome scales are POST HOC and decide nothing: no interval is placed on one and no"
    " comparison between chromosomes is registered. They are reported so a reader can see whether the"
    " magnitude is uniform or regional, with chrY's 13 genes flagged as below MIN_GENES",
    "nothing in area J, no node model, no BioForge calibration, no new reader biosample (the thirteen"
    " reader cells only, READER_TERMS unchanged), no AlphaGenome request and no deletion request: the stop"
    " list of 2026-09-27 as ROADMAP §5 reinstated it. No stop item is touched by this item, and no file of"
    " magcal's, multrna's, combine2's, freshrna's, trainexp's, emit7g's or scalefit's is edited",
)


def committed_facts() -> dict[str, Any]:
    """7h's and M2h's committed numbers the harness must reproduce, and their readings verbatim."""
    mres = json.loads(Path(MR_RESULT).read_text())
    gres = json.loads(Path(MAGCAL_RESULT).read_text())
    return {
        "harness_targets": {
            "seven_h": dict(ss.SEVEN_H_COMMITTED),
            "m2h": dict(ss.M2H_COMMITTED),
        },
        "seven_h_verbatim": {
            "source": MR_RESULT,
            "reading": mres["reading"],
            "outcome": mres["outcome"],
        },
        "m2h_verbatim": {
            "source": MAGCAL_RESULT,
            "reading": gres["per_form"]["multiplicative"]["reading"],
            "outcome": gres["per_form"]["multiplicative"]["outcome"],
            "split_rule": gres["per_form"]["multiplicative"]["split"]["rule"],
        },
        "note": (
            "quoted and not re-decided. M2h's reading is the one this lane audits for leakage through"
            " shared compiled elements; 7h's is the uncalibrated form it calibrates"
        ),
    }


def inputs() -> list[dict[str, Any]]:
    return [
        mf.input_entry(
            p,
            partition=(
                "pinned by sha256; the source of the committed P and D. Its counts, committed"
                " uncalibrated Delta MSE and registered reading are read; no scale, split or held-out"
                " value of this lane is computed"
                if p == MR_RESULT
                else "pinned by sha256; the committed numbers the harness must reproduce and the design"
                " this lane audits, quoted and not re-decided"
            ),
        )
        for p in (*SOURCE_FILES, HARNESS_RESULT)
    ]


def payload(entries: list[dict[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {
        "lane": ss.LANE,
        "date": date.today().isoformat(),
        "roadmap_item": (
            "§5, area I: does M2h's one scale survive a split that cannot leak through shared compiled"
            " elements -- held out by CHROMOSOME, registered"
        ),
        **ss.registered_fields(),
    }
    out.update(
        {
            "committed_facts": committed_facts(),
            "exploration_preceded_this": EXPLORATION_PRECEDED_THIS,
            "decisions": list(DECISIONS),
            "pinned_sha256": {x["path"]: x["sha256"] for x in entries},
            "what_the_run_will_write": WILL_WRITE,
            "the_code_that_runs_it": (
                "genomeos/attribution/scalesplit.py, scripts/scalesplit_harness.py and"
                " scripts/scalesplit_run.py, committed before this registration, with the committed"
                " magcal.py, multrna.py and combine2.py they import read-only for the words, constants"
                " and functions carried verbatim"
            ),
            "stop_list_checked": (
                "docs/ROADMAP.md §5, 'The stop list of 2026-09-27, reinstated', read before this lane was"
                " briefed: no AlphaGenome request and no deletion request; no new reader biosample (the 13"
                " covered cells only, READER_TERMS untouched); nothing in area J; no node-model work; no"
                " BioForge calibration. None is touched here"
            ),
            "alphagenome_requests": 0,
            "money": "none: no request is sent and no key is read",
        }
    )
    return out


def manifest(entries: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "sources": [
            {"accession": "data/results/multrna_result.json (7h, lane-multrna)", "version": "pinned"},
            {"accession": "data/results/magcal_result.json (M2h, lane-magcal)", "version": "pinned"},
            {
                "accession": "ENCODE polyA plus RNA-seq gene quantifications, V29, 13 reader biosamples",
                "version": "read by 7d, 7h and M2h; this lane reads their committed rows, not the files",
            },
        ],
        "inputs": entries,
        "input_count": len(entries),
        "assembly": "GRCh38; genes as 7h matched them, by gene_id",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "form": ss.FORM,
            "folds": list(ss.FOLDS),
            "split": ss.SPLIT,
            "refitted": ss.REFITTED,
            "min_genes": ss.MIN_GENES,
            "form_interval": list(ss.FORM_INTERVAL),
            "interval_decision": ss.INTERVAL,
            "bootstrap": ss.BOOTSTRAP,
            "permutations": ss.PERMUTATIONS,
            "control_alpha": ss.CONTROL_ALPHA,
            "tolerance": ss.TOLERANCE,
            "m2h_tolerance": ss.M2H_TOLERANCE,
            "free_parameters": "one real scale per fold; no per-gene and no per-cell parameter",
        },
        "exclusions": ["none of this lane's own: the population is 7h's, unchanged"],
        "partitions": "n/a: a registration; no chromosome-fold scale computed and no outcome read",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    if Path(WILL_WRITE).exists():
        raise SystemExit(f"refused: {WILL_WRITE} already exists; the run is not re-run")
    if not Path(HARNESS_RESULT).exists():
        raise SystemExit(f"refused: {HARNESS_RESULT} does not exist; run and commit the harness first")
    entries = inputs()
    print("saved", save_result(RESULT, payload(entries), manifest=manifest(entries)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
