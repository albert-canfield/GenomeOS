#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The M3 registration: what the cell-level closure predicts, written before a count is read.

This file carries no figure a run would fill. The arm SIZES are in it, and they are not outcomes:
an arm's size is fixed by the project's own prediction, by K562's dependency on a gene and by
whether the library covers the element, all of which were read by `closure_cell_inventory.py`
before this was written and none of which is the screen's measurement. The sha256 of the frozen
population is pinned here for the same reason: the set that gets measured must be the set that was
chosen, not the set the measurement suggested.

The run that follows (`scripts/closure_cell_run.py`) refuses to start unless this file is in a
commit and the population on disk still hashes to the digest below.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import closure_cell as cc  # noqa: E402
from genomeos.results import load_result, save_result  # noqa: E402
from genomeos.therapeutics import depmap  # noqa: E402

NAME = "closure_cell_registration"
ENTRY = "scripts/closure_cell_register.py"
OWN_CODE = (
    "genomeos/attribution/closure_cell.py",
    "scripts/closure_cell_inventory.py",
    "scripts/closure_cell_register.py",
    "scripts/closure_cell_run.py",
    "tests/test_closure_cell.py",
)
INVENTORY = "closure_cell_inventory"
#: Fewer than this in either of the two arms the primary comparison needs, or in the positive
#: control, and nothing is compared: the result records the stop and the reason.
MIN_PER_ARM = 30

#: The exact wording of every outcome, fixed here. The run fills the numbers into the sentence its
#: decision rule selects and changes no other word of it. Each one states its denominators, and
#: none of them says "equal": a comparison that does not separate two arms has not shown them to
#: be the same.
OUTCOMES = {
    "chain_supported": "The cell-level closure holds in K562 for the elements the attribution "
    "carries on its own: over {n_a} elements whose attributed gene is one K562 depends on and "
    "whose nearest coding gene is not, the median log2 relative guide abundance in the CRISPRi "
    "screen against the control is {median_a}, against {median_d} over the {n_d} elements where "
    "neither gene is one K562 depends on (one-sided Mann-Whitney p {p_ad}); {rate_a} of {n_a} fall "
    "at or below the floor arm's 5th percentile, which is 5% of {n_d} by construction. The "
    "positive control works in the same population: {n_p} promoter elements of genes K562 depends "
    "on have median {median_p} (p {p_pd} against the same floor).",
    "no_difference_detected": "No difference was detected at the cell level: over {n_a} elements "
    "whose attributed gene is one K562 depends on and whose nearest coding gene is not, the median "
    "log2 relative guide abundance in the CRISPRi screen against the control is {median_a}, "
    "against {median_d} over the {n_d} elements where neither gene is one K562 depends on "
    "(one-sided Mann-Whitney p {p_ad}); {rate_a} of {n_a} fall at or below the floor arm's 5th "
    "percentile, which is 5% of {n_d} by construction. The positive control does separate in the "
    "same population ({n_p} promoter elements of genes K562 depends on, median {median_p}, p "
    "{p_pd}), so the readout is responsive and this is a negative result about the attribution, "
    "not about the assay. No difference detected is not no difference.",
    "proximity_explains_it": "What the screen separates is proximity and not the attribution: the "
    "{n_c} elements whose NEAREST coding gene is one K562 depends on while their attributed gene "
    "is not have median {median_c} against the floor's {median_d} over {n_d} elements (p {p_cd}), "
    "while the {n_a} elements the attribution carries on its own have median {median_a} (p "
    "{p_ad}). The positive control works in the same population ({n_p} promoter elements, median "
    "{median_p}, p {p_pd}). On this readout the attribution adds nothing to the neighbourhood.",
    "readout_not_shown_to_work": "The result is uninterpretable as registered: the positive "
    "control does not separate. Over {n_p} promoter elements of genes K562 depends on the median "
    "log2 relative guide abundance is {median_p} against the floor arm's {median_d} over {n_d} "
    "elements (p {p_pd}), so this pairing of the screen with these element boundaries is not shown "
    "to detect a cell-level cost even where one must exist. The primary arms are reported for the "
    "record and decide nothing: {n_a} elements at median {median_a}, p {p_ad}.",
    "stopped_for_power": "Stopped at the registered stop rule, with nothing compared: the arms "
    "the primary comparison needs hold {n_a} and {n_d} measured elements and the positive control "
    "{n_p}, against the registered minimum of {min_per_arm}. No comparison was run and no reading "
    "is offered.",
}


def decision_rule() -> dict[str, str]:
    """Which outcome the run must select, written as the order the run applies it in."""
    return {
        "1_stop": "if any of arm A (attributed_essential_neighbourhood_not), arm D "
        "(attributed_not_neighbourhood_not) or the positive control holds fewer than "
        f"{MIN_PER_ARM} MEASURED elements, select stopped_for_power and compare nothing",
        "2_instrument": "else if the positive control's median is not below arm D's, or its "
        "one-sided p is at or above 0.05, select readout_not_shown_to_work",
        "3_chain": "else if arm A's median is below arm D's and the one-sided Mann-Whitney p for A "
        "below D is under 0.05, select chain_supported",
        "4_proximity": "else if arm C (attributed_not_neighbourhood_essential) separates from arm D "
        "on the same two conditions while arm A does not, select proximity_explains_it",
        "5_negative": "else select no_difference_detected",
        "ties": "the thresholds are read as written: 'under 0.05' is strict, 'below' is strict, and "
        "a p that is exactly 0.05 does not pass",
    }


def registration(inventory: dict[str, Any], population_sha: str) -> dict[str, Any]:
    arms = inventory["arms"]["sizes"]
    return {
        "lane": "lane-closurecell",
        "item": "ROADMAP section 5 item 14, M3: closure at the cell level (area I, Next 5)",
        "registered_in_its_own_commit": "This registration and the code that reads the outcome are "
        "committed before the first count of the screen is read. The inventory commit that precedes "
        "it holds the census and the frozen population; this commit holds the prediction, the "
        "baselines, the falsifier, the stop rule and the exact wording of every outcome. Nothing "
        "below is a number the run produces.",
        "question": "The gene-level closure asked whether a gene's attributed elements go with the "
        "RNA the cell makes. This asks the next link: when the project says element e drives gene G "
        "on K562's own track, and K562 is a cell that depends on G, does silencing e in K562 cost "
        "the cell its share of a pooled culture?",
        "chain_under_test": [
            "the project's claim: e's deletion changes G on the K562 RNA track "
            "(predicted_coding_by_cell['K562'], signed so that a positive effect means e holds G up)",
            f"the premise about the cell: K562 ({cc.K562_MODEL}) has a Chronos gene effect at or "
            f"below {cc.ESSENTIAL} for G in {depmap.RELEASE}",
            "the consequence measured: the guides targeting e are less abundant in ENCODE's CRISPRi "
            "proliferation screen in K562 than in ENCODE's control screen for it",
        ],
        "prediction": "Elements whose attributed gene is one K562 depends on are depleted in the "
        "CRISPRi screen relative to the control more than elements whose attributed gene is one "
        "K562 does not depend on, WHEN the element's nearest coding gene is in neither case a gene "
        "K562 depends on. The qualifier is the test: it is what makes the comparison about the "
        "attribution rather than about the neighbourhood.",
        "measured_readout": "Per element, log2 of its share of the CRISPRi library over its share "
        "of the control library: counts per million of the summed counts of the targeting guides "
        f"inside the element, pooled over the four replicates of each side, with a pseudocount of "
        f"{cc.PSEUDOCOUNT} on the element. Per-replicate-pair values are reported beside the pooled "
        "one and decide nothing.",
        "readout_is_not": "This is relative guide abundance at harvest in CRISPRi-modified cells "
        "against unmodified cells carrying the same library. It is not a growth rate, not a count "
        "of doublings and not a measurement of gene G.",
        "population": {
            "frozen_at": str(cc.POPULATION),
            "sha256": population_sha,
            "elements": inventory["frozen_population"]["elements"],
            "how_it_was_built": "every element of the project's all-enhancer deletion sweep with a "
            f"named coding target, at least {cc.MIN_GUIDES} of the screen library's guides inside "
            "it and a value on K562's own track; the census that built it read no count",
            "digest_is_a_gate": "the run recomputes this sha256 and refuses to measure if it moved",
        },
        "arms": {
            "sizes_known_in_advance": arms,
            "meaning": cc.ARM_MEANING,
            "positive_control": {
                "name": cc.POSITIVE_CONTROL,
                "size": arms.get(cc.POSITIVE_CONTROL),
                "what_it_is": "elements within "
                f"{cc.PROMOTER_WINDOW} bp of the TSS of a protein-coding gene K562 depends on. "
                "CRISPRi at a promoter silences that gene directly, so these must deplete; if "
                "they do not, the readout is not shown to work in this population. They are "
                "excluded from every other arm.",
            },
            "direction_arm": {
                "name": cc.DIRECTION_ARM,
                "size": arms.get(cc.DIRECTION_ARM),
                "what_it_is": "elements whose predicted effect on K562's own track has the "
                "opposite sign, attributed to a gene K562 depends on. The sign makes no depletion "
                "prediction, so this arm is reported apart and is in no comparison that decides.",
            },
        },
        "baselines": {
            "the_gene_s_dependency_alone": "arm C (attributed_not_neighbourhood_essential): the "
            "element's NEAREST coding gene is one K562 depends on and its attributed gene is not. "
            "This is what distance-to-nearest-gene plus essentiality predicts, with the project's "
            "attribution contributing nothing, and it is the baseline the chain has to beat.",
            "the_floor": "arm D (attributed_not_neighbourhood_not): neither gene is one K562 "
            "depends on. Every arm is read against this one.",
            "both_at_once": "arm B (attributed_essential_neighbourhood_essential): the claim holds "
            "but proximity would make it too, so a separation here is not attributable to either.",
            "the_assay_s_own_null": "the non-targeting guides of the same library, pooled per file. "
            "Their pooled log2 ratio is reported as the assay's zero and is not an arm.",
        },
        "matched_control": "Arms A and D are matched by construction on the three things that "
        "could otherwise carry the contrast: the element is distal (more than "
        f"{cc.PROMOTER_WINDOW} bp from any coding TSS), its nearest coding gene is one K562 does "
        "not depend on, and its prediction clears the same magnitude bar on the same cell's track. "
        "They differ in one thing: whether the gene the project attributes the element to is one "
        "K562 depends on. The guide count per element is reported per arm so a coverage difference "
        "cannot be mistaken for the effect.",
        "falsifier": "The chain is refuted on this readout if arm A's median log2 ratio is not "
        "below arm D's, or if the one-sided Mann-Whitney p for A below D is at or above 0.05, while "
        "the positive control does separate from D on the same two conditions. It is also refuted "
        "as an ADDITION if arm C separates from D and arm A does not: that is the neighbourhood "
        "doing the work the attribution claimed. Nothing here can establish that a separation is "
        "mediated by gene G: a pooled fitness screen cannot say through what an element acts.",
        "stop_rule": f"If arm A, arm D or the positive control holds fewer than {MIN_PER_ARM} "
        "MEASURED elements -- measured meaning at least "
        f"{cc.MIN_GUIDES} of its guides seen in the files and at least {cc.MIN_CONTROL_COUNT} "
        "pooled control reads -- no comparison is run, the stopped_for_power wording is the result, "
        "and the counts are reported as they stand. The run also stops, with nothing written, if "
        "the population's digest has moved, if a streamed file's md5 does not match the md5 ENCODE "
        "publishes, or if more than 0.1% of a file's rows do not fit the layout.",
        "thresholds": {
            "min_guides_per_element": cc.MIN_GUIDES,
            "min_pooled_control_count": cc.MIN_CONTROL_COUNT,
            "min_prediction_magnitude_on_the_cell_s_track": cc.MIN_EFFECT,
            "essential_at_or_below": cc.ESSENTIAL,
            "not_essential_at_or_above": cc.NOT_ESSENTIAL,
            "promoter_window_bp": cc.PROMOTER_WINDOW,
            "depletion_rate_cut": "the floor arm's own "
            f"{cc.DEPLETION_PERCENTILE}th percentile, so the floor arm's rate is "
            f"{cc.DEPLETION_PERCENTILE}% by construction",
            "alpha": 0.05,
            "sided": "one-sided, for depletion; the direction is part of the prediction",
            "min_per_arm": MIN_PER_ARM,
        },
        "decision_rule": decision_rule(),
        "outcome_wordings": OUTCOMES,
        "what_will_be_reported_whatever_happens": [
            "every arm's size, measured count, mean, median and depletion rate, each with its "
            "denominator beside it",
            "the mean guide count and the mean pooled control count per arm, so coverage is visible",
            "the non-targeting guides' pooled ratio as the assay's zero",
            "the per-replicate-pair readings beside the pooled one",
            "the row census of each streamed file, the malformed-row count and the md5 check",
        ],
        "cannot_establish": [
            "that a depleted element acts through the gene the project attributes it to: a pooled "
            "fitness screen measures the culture, not the gene",
            "that an element with no depletion is not regulatory: the screen's power per element is "
            "set by its guide count and the gene's dependency, and absence here is not absence",
            "anything about cells other than K562, or about any gene K562 does not depend on",
            "any direction claim about elements whose predicted sign is negative: that arm is "
            "reported apart and decides nothing",
        ],
        "prior_exposure": "The proliferation readout of this screen has never been read in this "
        "project. Klann 2021 appears in two places and neither touches it: the ENCODE CRISPRi "
        "benchmark's expression-effect pairs under data/knowledge/crispri, and "
        "scripts/prior_only_test.py, which classified that benchmark slice unscorable. The eight "
        "guide-quantification files were opened once each for their first two megabytes, to fix the "
        "column layout, and no count value was kept; the census result records that probe.",
        "money": "none: ENCODE and DepMap are public, and no model request of any kind is made",
        "alphagenome_requests": 0,
        "inventory_result": INVENTORY,
    }


def result_manifest(inventory_path: Path, population: Path) -> dict[str, Any]:
    return {
        "sources": [
            {"accession": cc.SCREEN_SOURCE, "version": "ENCODE portal, released; " + cc.ENCODE_LICENCE},
            {
                "accession": "DepMap CRISPRGeneEffect (Chronos), K562 = " + cc.K562_MODEL,
                "version": f"{depmap.RELEASE}, doi:{depmap.DOI}; {depmap.LICENCE}",
            },
            {
                "accession": "the project's own census and frozen population for this test",
                "version": f"data/results/{INVENTORY}.json as committed, digested in inputs",
            },
        ],
        "inputs": [
            mf.input_entry(
                inventory_path,
                partition="the census this registration's arm sizes are read from; it contains no outcome",
            ),
            mf.input_entry(
                population,
                partition="the frozen population, pinned by digest so the measured set cannot move",
            ),
        ],
        "assembly": "GRCh38 throughout; nothing is lifted",
        "coordinates": {
            "base": 0,
            "interval": "half-open",
            "note": "the frozen population's own element intervals and the files' 3 bp guide "
            "intervals; this registration computes no coordinate of its own",
        },
        "parameters": {
            "min_guides": cc.MIN_GUIDES,
            "min_pooled_control_count": cc.MIN_CONTROL_COUNT,
            "min_effect": cc.MIN_EFFECT,
            "essential": cc.ESSENTIAL,
            "not_essential": cc.NOT_ESSENTIAL,
            "promoter_window": cc.PROMOTER_WINDOW,
            "pseudocount": cc.PSEUDOCOUNT,
            "depletion_percentile": cc.DEPLETION_PERCENTILE,
            "min_per_arm": MIN_PER_ARM,
            "alpha": 0.05,
        },
        "exclusions": [
            "every element the census left out: fewer than MIN_GUIDES guides, or no value on K562's "
            "own track",
            "every element whose dependency call, for the attributed or the nearest coding gene, is "
            "borderline or unmeasured",
            "promoter elements from every arm but the positive control",
            "every count of the screen: none is read by this writer",
        ],
        "partitions": "n/a: a registration written before the measurement; nothing is evaluated",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


#: Keys no registration of this test may carry: each one is a figure only the run can produce.
FORBIDDEN = (
    "median_a",
    "median_d",
    "median_c",
    "median_p",
    "p_ad",
    "p_cd",
    "p_pd",
    "rate_a",
    "readout",
    "log2_ratio",
    "per_element",
    "outcome",
    "verdict",
    "reading",
    "depletion_rate",
    "library_total",
)


def forbidden_at_any_depth(payload: Any, path: str = "") -> list[str]:
    """Every forbidden key anywhere in the payload, named by its path. Recursive on purpose."""
    found: list[str] = []
    if isinstance(payload, dict):
        for key, value in payload.items():
            here = f"{path}.{key}" if path else key
            if key in FORBIDDEN:
                found.append(here)
            found.extend(forbidden_at_any_depth(value, here))
    elif isinstance(payload, list):
        for i, value in enumerate(payload):
            found.extend(forbidden_at_any_depth(value, f"{path}[{i}]"))
    return found


def main() -> int:
    mf.trace_begin()
    inventory = load_result(INVENTORY)
    if inventory is None:
        raise SystemExit(f"no {INVENTORY} result: run scripts/closure_cell_inventory.py first")
    inventory_path = ROOT / "data/results" / f"{INVENTORY}.json"
    if not cc.POPULATION.exists():
        raise SystemExit(f"no frozen population at {cc.POPULATION}")
    population_sha = hashlib.sha256(cc.POPULATION.read_bytes()).hexdigest()
    pinned = inventory["frozen_population"]["sha256"]
    if population_sha != pinned:
        raise SystemExit(
            f"the population on disk ({population_sha}) is not the one the census froze ({pinned})"
        )
    record = registration(inventory, population_sha)
    bad = forbidden_at_any_depth(record)
    # The outcome wordings are templates and hold the placeholders by name, so they are allowed
    # under `outcome_wordings` and nowhere else.
    bad = [k for k in bad if not k.startswith("outcome_wordings")]
    if bad:
        raise SystemExit(f"a registration must hold no figure the run would fill: {bad}")
    record["result_manifest"] = result_manifest(inventory_path, cc.POPULATION)
    path = save_result(NAME, record)
    print(f"wrote {path}, sha256 {hashlib.sha256((ROOT / path).read_bytes()).hexdigest()}")
    print(f"population pinned at {population_sha}")
    print(f"arm sizes registered: {json.dumps(record['arms']['sizes_known_in_advance'], sort_keys=True)}")
    print(f"outcome wordings registered: {sorted(OUTCOMES)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
