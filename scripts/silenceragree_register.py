#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Write the silencer_like independent-check registration, committed ALONE and before any agreement
rate exists.

The file this writes carries no field a run would fill. There is no agreement count, no rate, no
band, no control figure and no key holding null awaiting a number: those fields do not EXIST in it
rather than sitting empty. The check that says so runs RECURSIVELY, at every depth of the payload
and again on the file as written, because a forbidden key nested one level down would otherwise pass.

The scope this registration states is built through `silenceragree.scope()`, which takes its element
data only through `element_index()` and is handed no action field at all, so the registration cannot
carry a direction even by mistake. The one direction count in the file is the genome-wide BASE RATE,
which is the test set's prior and not its agreement rate, and which the earlier probe already
published.
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
from genomeos.attribution import silenceragree as sa  # noqa: E402
from genomeos.results import save_result  # noqa: E402

NAME = "silenceragree_registration"
ENTRY = "scripts/silenceragree_register.py"
OWN_CODE = (
    "genomeos/attribution/silenceragree.py",
    "scripts/silenceragree_fetch.py",
    "scripts/silenceragree_register.py",
    "scripts/silenceragree_run.py",
    "tests/test_silenceragree.py",
)


def forbidden_at_any_depth(payload: Any, path: str = "") -> list[str]:
    """Every forbidden key anywhere in a nested payload, named by its path.

    Recursive on purpose: the earlier registrations checked the top level only, and a run's figure
    nested one dict down would have passed unseen. This supersedes that check additively -- it
    refuses everything the flat check refused and more.
    """
    found: list[str] = []
    if isinstance(payload, dict):
        for key, value in payload.items():
            here = f"{path}.{key}" if path else key
            if key in sa.FORBIDDEN_IN_REGISTRATION:
                found.append(here)
            found.extend(forbidden_at_any_depth(value, here))
    elif isinstance(payload, list):
        for i, value in enumerate(payload):
            found.extend(forbidden_at_any_depth(value, f"{path}[{i}]"))
    return found


def _stratum_base() -> dict[str, Any]:
    """The informative base rate, stated in advance because a bare agreement rate is uninterpretable.

    This reads a direction, so it is NOT built through `scope()`. It is the base rate of the
    population the ReSE elements were drawn from, with the ReSE elements themselves excluded, so it
    carries no information about their own directions.
    """
    kept, _ = sa.rese_intervals(sa.rese_records())
    rows, _ = sa.join(kept)
    return sa.stratum_matched_base_rate(rows)


def registration() -> dict[str, Any]:
    return {
        "study": "silenceragree",
        "registered": "2026-10-03",
        "lane": "lane-silenceragree",
        "question": "Does the model's signed direction -- `silencer_like`, the winner-take-all claim "
        "that the largest-moving gene in a deleted element's 1 Mb window moved UP -- agree with "
        "assay-validated repressive ability on the elements where both instruments have an opinion? "
        "This is the FIRST independent check this layer has ever had: 230,839 of 612,323 named "
        "elements carry silencer_like across 75 committed results and nothing outside the model has "
        "ever been consulted about a single one of them.",
        "the_two_quantities_are_not_the_same_thing": {
            "the_derived_call": "a DIRECTION over a whole 1 Mb window: of every gene in the window, "
            "the one that moves most in one deterministic AlphaGenome run, and whether its move was "
            "up. One gene, one run, no spread, no probability (review R4), ties to activation.",
            "the_assay": "a FRAGMENT'S repressive ability in a lentiviral reporter: ReSE clones the "
            "fragment beside an apoptosis-inducing fusion and scores cell survival under "
            "transcriptional repression. The fragment is out of its locus and the readout is one "
            "reporter, not the native neighbourhood.",
            "so": "a fragment can repress a reporter it is cloned beside while the model's largest "
            "window mover in the native locus is a different gene going the other way, and neither "
            "instrument is wrong when that happens. Agreement is a point of contact, not a shared "
            "definition.",
            "licenses": sa.WHAT_AGREEMENT_LICENSES,
            "does_not_license": sa.WHAT_AGREEMENT_DOES_NOT_LICENSE,
        },
        "may_not_be_concluded": "No element is given molecular_role: silencer by this study and no "
        "lane proposes it here. genomeos/lang/grammar.py defines that value as 'a repressive "
        "element; never inferred from a direction of effect alone', and a direction of effect is "
        "exactly what is read. genomeos/attribution/increase_links.py forbids the words silencer, "
        "repressor, represses, repression and represser of an increase link, and nothing here "
        "weakens that. Whether an assay-validated source can satisfy evidence_status and license "
        "that value is a SEPARATE registered decision and it is not this lane's.",
        "source": sa.SOURCE,
        "silencerdb_not_used": "SilencerDB's validated tier shares all 4,998 ReSE intervals at "
        "offset exactly (0,0) -- the same data -- and SilencerDB states no licence, terms or "
        "copyright anywhere, which makes vendoring it a decision for Albert. Its genuinely "
        "additional 3,599 intervals are not needed for this check and taking them would manufacture "
        "a licence question the work does not require. It is not fetched, not read and not vendored.",
        "the_grch37_trap": {
            "what": "every one of the 4,998 ReSE records states its interval in its DESCRIPTION in "
            "GRCh37, and the GRCh38 coordinates are in `genomicinfo`, lifted by RefSeq.",
            "so": "the description's own numbers must never be joined against this project's GRCh38 "
            "elements, and this study joins on genomicinfo alone. A guard holds it.",
            "and": "480 further records name ReSE only elsewhere in the record. For those the record "
            "is a larger functional element -- a promoter, an enhancer -- that merely CONTAINS the "
            "assayed fragment, so genomicinfo is the parent's span: measured median 3.55x the "
            "fragment and up to 138.9x. They are excluded, and a guard holds that too.",
        },
        "prior_exposure": sa.PRIOR_EXPOSURE,
        "independence_of_the_source_from_the_training_data": sa.INDEPENDENCE,
        "nearest_prior_art_in_the_project": sa.PRIOR_ART,
        "join_rule": sa.JOIN_RULE,
        "scope": sa.scope(),
        "the_denominator": {
            "statement": "whatever this study finds is a statement about the joined elements and "
            "must be worded as one. The joined set is a fraction of one percent of the elements that "
            "carry the call, in two cell lines of the four scored, and in two cCRE classes of the "
            "nine the project holds.",
            "of_the_elements_carrying_silencer_like": "the scope block's `denominator` against the "
            "base_rate block's `silencer_like`; both are stated so the ratio is read off the file "
            "rather than asserted in it",
            "cell_lines": "K562 and HepG2 only",
            "element_classes": "dELS and pELS only",
            "not_covered": "403 ReSE fragments overlap a scored element on which NO gene moved by "
            "MIN_EFFECT, so the model names no target and there is no direction to compare. They are "
            "outside the denominator and are counted in the join block, not discarded quietly.",
        },
        "base_rates_stated_in_advance": sa.BASE_RATES,
        "base_rate_global_for_comparison": sa.base_rate(),
        "base_rate_stratum_matched": _stratum_base(),
        "why_a_control_is_not_optional": "about 37.7% of named elements are already silencer_like, so "
        "a matched random set of scored elements agrees at roughly that rate BY CONSTRUCTION. An "
        "agreement rate quoted without its control is uninterpretable, and every ratio in the result "
        "is to be published beside its absolute level and its control.",
        "control": sa.CONTROL,
        "bands": sa.BANDS,
        "bands_are_relative_not_absolute": sa.BANDS_ARE_RELATIVE,
        "band_threshold_points": sa.BAND_THRESHOLD_POINTS,
        "interval": sa.INTERVAL,
        "control_interval": sa.CONTROL_INTERVAL,
        "power": sa.POWER,
        "an_inconclusive_is_a_complete_reading": "If the registered power run cannot tell the bands "
        "apart at this denominator, the reading is INCONCLUSIVE_UNDERPOWERED and that IS the "
        "finding. 1,181 elements clustered in 23 chromosomes may simply not distinguish a 10-point "
        "displacement from none, and reporting that is worth more than a point estimate quoted "
        "without its interval.",
        "readings": {
            "1_as_shipped": sa.READING_1,
            "2_cell_matched": sa.READING_2,
            "3_tie_bias_descriptive": sa.READING_3,
        },
        "tie_bias": sa.TIE_BIAS,
        "four_cell_line_negative": sa.FOUR_LINE_NEGATIVE,
        "falsifiers": sa.FALSIFIERS,
        "authorises_no_conclusion": "This file states a plan and a rule. It holds no agreement count, "
        "no agreement rate, no band and no control figure, and the fields a run would fill DO NOT "
        "EXIST in it rather than sitting empty, so that nothing in it can read as a preview. The "
        "check that says so is recursive and runs twice: once on the payload and once on the file AS "
        "WRITTEN, because genomeos/results.py adds the registry's own stamp after the payload leaves "
        "the writer. That stamp uses one forbidden spelling -- `result` -- and it holds this "
        "registration's own NAME and not an outcome, so it is allowed at exactly that value and "
        "nothing else.",
        "cost": {
            "alphagenome_requests": 0,
            "paid_requests": 0,
            "money": "none",
            "network": "free NCBI eutils, no API key, and nothing else",
        },
    }


#: The contract's shape for `coordinates`, which `save_result` requires as {base, interval} and
#: which `result_manifest()` first stated as prose. ADDED rather than edited in place: the guard
#: refused removing the prose line -- it is work already in HEAD -- and this project supersedes
#: additively. The prose sentence is kept verbatim below as the note, so nothing it said is lost and
#: the declaration a reader of the written file sees is this one.
COORDINATES = {
    "base": 0,
    "interval": "half-open",
    "note": "0-based half-open for the project's elements as stored; the NCBI genomicinfo span is "
    "used as given, which costs at most one base at each end of a fragment whose median length is "
    "192 bp, and overlap is tested as element_start < rese_end and element_end > rese_start. The "
    "overlap in bases is recorded on every joined row, so a stricter rule can be applied afterwards "
    "without re-joining.",
}

#: The 23 committed cCRE subsets the control matches its STRATA on. Declared as a group because the
#: first write was quarantined for reading all 23 without declaring one of them: the file that
#: decides every stratum boundary was invisible to a rebuild.
def ccre_group():
    return mf.files_entry(
        "ENCODE cCRE v3 per-chromosome subsets, read for the element CLASS the control matches on",
        sorted((ROOT / "data/results").glob("ccres_chr*.bed.gz")),
        partition="the registry's own committed cCRE subsets",
    )


def result_manifest() -> dict[str, Any]:
    return {
        "sources": [
            {
                "accession": "ReSE screen-validated silencers in RefSeq Functional Elements, "
                "NCBI Gene esummary, query 'ReSE[All Fields] AND human[orgn]'",
                "version": "fetched 2026-10-03; 5,478 records; NCBI places no use restrictions",
            },
            {
                "accession": "Pang B, Snyder MP. Systematic identification of silencers in human "
                "cells. Nat Genet 2020;52(3):254-263 (what the assay measures)",
                "version": "10.1038/s41588-020-0578-5",
            },
            {
                "accession": "ENCODE cCRE v3 element classes (dELS, pELS), per-chromosome subsets "
                "already committed as data/results/ccres_chr*.bed.gz",
                "version": "as committed",
            },
        ],
        "inputs": [
            mf.input_entry(sa.RESE_CACHE, partition="the external source; data/cache is never committed"),
            mf.input_entry(
                sa.ELEMENTS_DIR,
                partition="the project's own scored elements; data/knowledge is never committed",
            ),
        ],
        "assembly": "GRCh38. The ReSE intervals are joined on NCBI genomicinfo, which is GRCh38; the "
        "GRCh37 coordinates in each record's description are carried for audit and never joined on.",
        "coordinates": "0-based half-open for the project's elements as stored; the NCBI genomicinfo "
        "span is used as given and overlap is tested as est < rese_end and een > rese_start",
        "parameters": {
            "control_draws": sa.CONTROL["draws"],
            "control_seed": sa.CONTROL["seed"],
            "band_threshold_points": sa.BAND_THRESHOLD_POINTS,
            "near_zero": sa.NEAR_ZERO,
            "rese_cell_lines": list(sa.RESE_CELL_LINES),
            "project_cell_lines": list(sa.PROJECT_CELL_LINES),
        },
        "exclusions": [
            "480 records naming ReSE only outside the description: the record is a larger element "
            "containing the fragment, so its genomicinfo is a parent span",
            "4 records with no genomicinfo and 2 unplaced on GRCh38",
            "3,408 ReSE fragments overlapping no scored element at all",
            "403 fragments whose overlapping element has no named predicted target",
            "SilencerDB entirely: no licence is stated anywhere in it",
        ],
        "partitions": "n/a: a registration committed before any measurement; nothing is evaluated",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    mf.trace_begin()
    rec = registration()
    bad = forbidden_at_any_depth(rec)
    if bad:
        raise SystemExit(f"the registration must hold no field a run would fill: {bad}")
    rec["result_manifest"] = result_manifest()
    # Both lines SUPERSEDE what result_manifest() stated, additively: the guard refuses removing a
    # line that is already in HEAD, so the prose `coordinates` and the short input list stay where
    # they are and the correct values are written over them here. What lands in the file is this.
    rec["result_manifest"]["coordinates"] = COORDINATES
    rec["result_manifest"]["inputs"].append(ccre_group())
    path = save_result(NAME, rec)
    written = json.loads((ROOT / path).read_text())
    stamp = written.get("result")
    if stamp != NAME:
        raise SystemExit(f"the registry's name stamp holds {stamp!r}, not this registration's name")
    bad = [k for k in forbidden_at_any_depth(written) if k != "result"]
    if bad:
        raise SystemExit(f"the written registration holds fields a run would fill: {bad}")
    print(f"wrote {path}, sha256 {hashlib.sha256((ROOT / path).read_bytes()).hexdigest()}")
    print(f"no forbidden key at any depth, out of {sorted(sa.FORBIDDEN_IN_REGISTRATION)}")
    print(f"the one exception, checked by value: result == {stamp!r}, the registry's name stamp")
    print(f"denominator registered: {rec['scope']['denominator']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
