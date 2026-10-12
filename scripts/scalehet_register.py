# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register lane-scalehet before it is run (data/results/scalehet_registration.json).

    uv run python scripts/scalehet_register.py

Every binding word is imported from `genomeos.attribution.scalehet` (its BOUND constants), and
scripts/scalehet_run.py refuses to run if any registered field differs from that module or any pinned
file's sha256 differs from the file on disk.

This script computes NO heterogeneity statistic, NO range, NO permutation, NO p value and NO
per-chromosome scale. It reads 7h's, M2h's and scalesplit's committed results only for their counts,
the committed numbers the harness must reproduce, and their registered readings, which are quoted
verbatim as context. No request is sent, no key is read.
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
from genomeos.attribution import scalehet as sh  # noqa: E402
from genomeos.results import save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scalehet_harness import (  # noqa: E402
    MAGCAL_RESULT,
    MR_RESULT,
    OWN_CODE,
    READERS,
    SCALESPLIT_RESULT,
    SOURCE_FILES,
)

RESULT = "scalehet_registration"
WILL_WRITE = "data/results/scalehet_result.json"
HARNESS_RESULT = "data/results/scalehet_harness.json"

EXPLORATION_PRECEDED_THIS = (
    "stated so that nothing here is read as a blind prediction. Before writing it this lane read:"
    " genomeos/attribution/magcal.py's constants and its fit, bootstrap, control, decision and reading"
    " functions in full; genomeos/attribution/scalesplit.py's registration constants, its POST_HOC"
    " paragraph and its committed M2H_COMMITTED and SEVEN_H_COMMITTED numbers;"
    " scripts/scalesplit_harness.py in full, which this lane's harness follows; the GENE COUNT PER"
    " CHROMOSOME of 7h's committed table (chr1 1,700, chr2 1,089, chr3 930, chr4 642, chr5 732, chr6 886,"
    " chr7 773, chr8 559, chr9 654, chr10 631, chr11 983, chr12 871, chr13 280, chr14 522, chr15 509,"
    " chr16 703, chr17 969, chr18 238, chr19 1,199, chr20 464, chr21 156, chr22 373, chrX 633, chrY 13),"
    " which is what told it that chrY alone falls below MIN_GENES; and the shape of"
    " data/results/scalesplit_result.json. IT ALSO KNEW, FROM ITS OWN BRIEF AND BEFORE REGISTERING,"
    " scalesplit's POST HOC per-chromosome spread: the range +0.1232 (chr13) to +0.2874 (chr14), the"
    " median +0.2176, the quartiles +0.2041 and +0.2376, the absence of a sign change, that 10 of the 24"
    " lie inside M2h's committed interval, and that chrY was flagged at 13 genes. THAT IS WHY THE PRIOR"
    " EXPOSURE PARAGRAPH SAYS THIS IS A TEST OF A PATTERN THAT HAS BEEN LOOKED AT. It has computed NO"
    " heterogeneity statistic, NO range statistic, NO permutation of the chromosome labels and NO p value"
    " of any kind; this lane's step 1 (data/results/scalehet_harness.json, committed before this"
    " registration) states in its own words that it wrote no per-chromosome scale value and formed no"
    " dispersion. The first statistic of this design will be computed by scripts/scalehet_run.py after"
    " this registration is committed"
)
DECISIONS = (
    "WHY THIS LANE EXISTS: M2h's one scale is ONE NUMBER for the whole genome, and that one number is"
    " what makes the compiled rules read better than no rule (about 2% of the variance). scalesplit then"
    " showed, post hoc, within-chromosome scales from +0.1232 to +0.2874. If that spread is real, the"
    " defect is not one mis-specified constant but an error that depends on the locus, which is a"
    " different defect and a different repair. NOBODY HAS TESTED IT: a spread that wide is also what"
    " sampling noise looks like at a few hundred genes per chromosome",
    "design: M2h's rows, M2h's fit, M2h's floor, M2h's 98.33% width and M2h's alpha, every one imported"
    " read-only from genomeos/attribution/magcal.py, with ONE thing added -- a heterogeneity statistic on"
    " the per-chromosome fitted scales and a null for it. No new endpoint, no new measurement, no new"
    " reader biosample, no new composition form",
    "the PRIMARY statistic is a PRECISION-WEIGHTED DISPERSION, Q = sum over admitted chromosomes of pp_c"
    " * (a_c - a_g)^2, Cochran's homogeneity form, where pp_c is that chromosome's sum of P^2 and is"
    " exactly the precision of its fitted scale. It uses every admitted chromosome and names none, so it"
    " is not a statistic about chr13 or chr14",
    "the SECONDARY statistic is the RANGE, max minus min, because that is the statistic the post hoc"
    " sentence used and this lane should test the pattern that was looked at rather than a kinder one."
    " It cannot declare either answer by itself: its only registered force is the MIXED outcome",
    "THE NULL IS A PERMUTATION OF THE CHROMOSOME LABELS AMONG THE GENES, holding EVERY CHROMOSOME'S GENE"
    " COUNT EXACTLY at its observed value, with each gene moving as one block carrying all thirteen of"
    " its pairs. That is the one null that answers the actual objection: it has the same per-chromosome"
    " gene counts and therefore the same sampling noise as the observed spread, it redistributes the same"
    " genes under the same single global fit, it prices the upward bias of fitting 23 scales because it"
    " fits 23 scales the same way, and it assumes no distribution. Its assumption -- gene"
    " exchangeability between chromosomes -- is written as a LIMITATION, so a rejection is a finding that"
    " the scale depends on something about the group and NOT an identification of what",
    "NO TRAIN-AND-HELD-OUT SPLIT, deliberately: the question is about the dispersion of FITTED scales,"
    " which is in-sample by construction, and the null prices that exactly. The held-out question of"
    " whether the calibrated rules beat no rule is M2h's and scalesplit's and is quoted, never re-decided",
    "THE GENE FLOOR IS FIXED NOW at MIN_GENES = 100, magcal's own, which excludes chrY (13 genes) and"
    " admits the other 23 chromosomes. It is not tuned here and it was scalesplit's floor before this"
    " lane existed. Fewer than MIN_CHROMOSOMES = 10 admitted would read INCONCLUSIVE",
    "THE THRESHOLD IS ALPHA = 0.05, one-sided, magcal's own, AND IS THE LAXER OF THE TWO DEFENSIBLE"
    " CHOICES. This lane expects 'consistent with one scale', which is a failure to reject, so it takes"
    " the threshold that gives heterogeneity the easier passage and makes its own expected answer the"
    " harder to obtain. No Bonferroni is applied over the two statistics, for the same reason; the"
    " multiplicity is handled by the MIXED outcome, which claims neither answer when they disagree",
    "the 98.33% width, magcal.FORM_INTERVAL unchanged, is used for DESCRIPTION ONLY -- the null Q"
    " distribution, the null range distribution, and each admitted chromosome's own null band -- and the"
    " verdict rests on the permutation p alone. The 23 per-chromosome bands carry no correction over 23"
    " and are declared NOT a test of any chromosome",
    "1000 permutations at seed 20261012, fixed now. 1000 is magcal's own count; the attainable floor"
    " 0.000999 is far below ALPHA, so the test can reject",
    "THE EXPOSURE IS DOUBLE AND IS STATED FIRST, NOT AS A CAVEAT: this RNA was read by items 7, 7d, 7h,"
    " the capped-sum census, M2h, M2i, scalesplit and lane-recovered, so this is a re-analysis of a"
    " measurement already read; AND the per-chromosome spread tested here was ALREADY SEEN POST HOC by"
    " scalesplit, so this is a test of a pattern that has been looked at, which is weaker than a test of"
    " one that has not",
    "IF THE OUTCOME IS HETEROGENEOUS, THAT VARIATION EXISTS IS THE WHOLE FINDING. No regional model is"
    " proposed, no per-region scale is fitted, no cause is named and no chromosome is read as a verdict"
    " about itself, in the reading, in the result or in the ROADMAP row. Fitting per-region scales is a"
    " separate registration and a separate decision",
    "nothing in area J, no node model, no BioForge calibration, no new reader biosample (the thirteen"
    " reader cells only, READER_TERMS unchanged), no AlphaGenome request and no deletion request: the"
    " stop list of 2026-09-27 as ROADMAP §5 reinstated it. No stop item is touched by this item, and no"
    " file of magcal's, scalesplit's, recovered's, scalefit's, freshrna's, trainexp's or multrna's is"
    " edited",
)


def committed_facts() -> dict[str, Any]:
    """The committed numbers the harness must reproduce, and the readings quoted verbatim."""
    mres = json.loads(Path(MR_RESULT).read_text())
    gres = json.loads(Path(MAGCAL_RESULT).read_text())
    sres = json.loads(Path(SCALESPLIT_RESULT).read_text())
    return {
        "harness_targets": {
            "seven_h": dict(sh.SEVEN_H_COMMITTED),
            "m2h": dict(sh.M2H_COMMITTED),
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
        "scalesplit_verbatim": {
            "source": SCALESPLIT_RESULT,
            "outcome": sres["outcome"],
            "transport_verdict": sres["transport_verdict"],
            "reading": sres["reading"],
            "post_hoc_rule": sres["post_hoc_deciding_nothing"].get("rule"),
        },
        "note": (
            "quoted and not re-decided. M2h's reading is the one whose MAGNITUDE this lane asks about;"
            " scalesplit's post hoc per-chromosome spread is the PATTERN this lane tests, and that it was"
            " already seen is why the exposure paragraph calls this a test of a looked-at pattern"
        ),
    }


def inputs() -> list[dict[str, Any]]:
    return [
        mf.input_entry(
            p,
            partition=(
                "pinned by sha256; the source of the committed P and D. Its counts, committed"
                " uncalibrated Delta MSE and registered reading are read; no statistic of this lane is"
                " computed"
                if p == MR_RESULT
                else "pinned by sha256; committed numbers and readings, quoted and not re-decided"
            ),
        )
        for p in (*SOURCE_FILES, HARNESS_RESULT)
    ]


def payload(entries: list[dict[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {
        "lane": sh.LANE,
        "date": date.today().isoformat(),
        "roadmap_item": (
            "§5, area I: is M2h's magnitude ONE global scale that the compiled rules are wrong by, or"
            " does it vary by region? Registered, with the null and the exposure"
        ),
        **sh.registered_fields(),
    }
    out.update(
        {
            "committed_facts": committed_facts(),
            "exploration_preceded_this": EXPLORATION_PRECEDED_THIS,
            "decisions": list(DECISIONS),
            "readers": READERS,
            "pinned_sha256": {x["path"]: x["sha256"] for x in entries},
            "what_the_run_will_write": WILL_WRITE,
            "the_code_that_runs_it": (
                "genomeos/attribution/scalehet.py, scripts/scalehet_harness.py and"
                " scripts/scalehet_run.py, committed before this registration, with the committed"
                " magcal.py and scalesplit.py they import read-only for the words, constants and"
                " functions carried verbatim"
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
                "accession": "data/results/scalesplit_result.json (lane-scalesplit)",
                "version": "pinned; its post hoc per-chromosome spread is the pattern this lane tests",
            },
            {
                "accession": "ENCODE polyA plus RNA-seq gene quantifications, V29, 13 reader biosamples",
                "version": "read by 7d, 7h and M2h; this lane reads their committed rows, not the files",
            },
        ],
        "inputs": entries,
        "input_count": len(entries),
        "readers": READERS,
        "reader_count": len(READERS),
        "assembly": "GRCh38; genes as 7h matched them, by gene_id",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "form": sh.FORM,
            "grouping": sh.GROUPING,
            "statistic": list(sh.STATISTIC),
            "null": list(sh.NULL),
            "no_split": sh.NO_SPLIT,
            "min_genes": sh.MIN_GENES,
            "min_chromosomes": sh.MIN_CHROMOSOMES,
            "form_interval": list(sh.FORM_INTERVAL),
            "interval_decision": sh.INTERVAL,
            "permutations": sh.PERMUTATIONS,
            "control_alpha": sh.CONTROL_ALPHA,
            "alpha_decision": sh.ALPHA,
            "tolerance": sh.TOLERANCE,
            "m2h_tolerance": sh.M2H_TOLERANCE,
            "free_parameters": (
                "one real scale per admitted chromosome and one global scale, all in closed form; no"
                " per-gene and no per-cell parameter, and nothing is fitted that the null does not refit"
            ),
        },
        "exclusions": [
            "chromosomes carrying fewer than MIN_GENES = 100 genes in 7h's committed table, by the floor"
            " fixed in this registration before any value of this lane was read; on 7h's table that is"
            " chrY alone, at 13 genes, which removes 13 of 16,509 genes"
        ],
        "partitions": "n/a: a registration; no statistic computed, no permutation drawn, no outcome read",
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
