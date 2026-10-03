# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register the remaining one-target consumers before any of them is moved.

    uv run --frozen python scripts/onetarget2_register.py

Writes data/results/onetarget2_registration.json. Every binding word is imported from
`genomeos.attribution.onetarget2`, the module that applies them, so the registration and the code
cannot drift apart; and the census's declared lines are re-derived from the tree by
`check_census_matches_tree()` before the file is written, so a registration naming a line that has
moved refuses to be written at all.

The payload is checked against `check_no_preview()` before it is saved: it holds no measurement of
any consumer's output, no verdict, no float and no field a run would fill. The only figures in it
are the SCOPE (how many modules read a compact head, and at which lines) and PRIOR EVIDENCE quoted
from a committed result with its date and commit.

No consumer is run, no element table and no response cache is opened, no request is sent and no key
is read. 0 AlphaGenome requests, no network, no money.
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import onetarget2 as ot  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "onetarget2_registration"


def inputs() -> list[dict[str, Any]]:
    """The one file this registration reads: the committed result THE_INVARIANT is quoted from."""
    p = Path(__file__).resolve().parents[1] / ot.PRIOR_EVIDENCE
    return [
        mf.input_entry(
            p,
            partition="its `candidates` arm's head_control only: the 849,469 elements compared and "
            "the 0 head disagreements THE_INVARIANT quotes. Read to CHECK the quote, not to take a "
            "new measurement; no other field of it is read",
        )
    ]


def payload() -> dict[str, Any]:
    """The registration. Nothing here is measured and nothing here is a preview."""
    entries = inputs()
    return {
        "result": RESULT,
        "date": date.today().isoformat(),
        "registered": ot.REGISTERED,
        "lane": ot.LANE,
        "wave": ot.WAVE,
        "row": ot.ROW,
        "question": ot.QUESTION,
        "authorises_no_conclusion": ot.AUTHORISES_NO_CONCLUSION,
        "exploration_preceded_this": ot.EXPLORATION_PRECEDED_THIS,
        "the_shared_reader": dict(ot.THE_SHARED_READER),
        "the_invariant_already_measured": ot.THE_INVARIANT,
        "the_quote_is_checked_against_the_file": (
            f"check_the_quoted_invariant() loads {ot.PRIOR_EVIDENCE} and refuses unless its "
            f"`{ot.PRIOR_EVIDENCE_ARM}` arm's head_control holds compared "
            f"{ot.PRIOR_EVIDENCE_COMPARED:,}, head_disagreements {ot.PRIOR_EVIDENCE_DISAGREEMENTS} "
            f"and coding_head_disagreements {ot.PRIOR_EVIDENCE_DISAGREEMENTS}. The figure is the "
            "committed result's, re-read here and not recalled"
        ),
        "what_the_projection_loses": ot.WHAT_THE_PROJECTION_LOSES,
        "the_count_the_row_states_and_the_count_the_tree_holds": (ot.ROW_SAYS_FIVE_THE_TREE_SAYS_OTHERWISE),
        "scope_counts": dict(ot.SCOPE_COUNTS),
        "classes": list(ot.CLASSES),
        "statuses": list(ot.STATUSES),
        "census": [dict(c) for c in ot.CENSUS],
        "record_only": [dict(m) for m in ot.RECORD_ONLY],
        "already_moved_each_under_its_own_registration": [dict(m) for m in ot.ALREADY_MOVED],
        "excluded_by_name": [dict(m) for m in ot.EXCLUDED_BY_NAME],
        "scripts_out_of_scope": list(ot.SCRIPTS_OUT_OF_SCOPE),
        "closure_frozen": list(ot.CLOSURE_FROZEN),
        "closure_rule": (
            "the paid study's authorisation is bound to a 52-file import closure and "
            "sender_closure() hashes the WORKING TREE, so even an UNCOMMITTED edit to one of those "
            "files refuses the frozen send. The two identity-class consumers inside it keep their "
            "falsifier registered and unfired, and whether they move is Albert's call"
        ),
        "peer_held_rule": (
            "a status of peer_held is not a decision and not a blocker of this wave's reasoning: "
            "the commit gate refuses a staged path a live peer lane holds on the work board, and "
            "naming the holder is the coordinator's business. No hold is released by this lane"
        ),
        "a_falsifier_per_module_not_one_for_the_set": (
            "each census entry carries its own falsifier, able to fail on its own, and "
            "check_one_falsifier_per_module() refuses two modules sharing one verbatim. One "
            "falsifier for the set would let most of the set pass on the strength of one member, "
            "which is what the row's wording exists to prevent"
        ),
        "superseding_is_additive": (
            "no test, pin or registration is weakened by this wave. Where a consumer's reading is "
            "superseded, the old assertion is kept as a strict expected failure with its reason "
            "beside it and the new one added, as dbb5d4a did"
        ),
        "nothing_is_moved_by_this_file": (
            "this registration edits no consumer. It is committed ALONE, before any code that "
            "produces a figure, and the first consumer moves in a later commit"
        ),
        "cannot_establish": list(ot.CANNOT_ESTABLISH),
        "refusals": list(ot.REFUSALS),
        "stop_list_not_approached": (
            "nothing is built on the node model, no biosample is added, nothing is built in area J, "
            "no AlphaGenome request is made on the 882 blocks, area E's fate ceiling is not chased "
            "and BioForge's confidence is not calibrated. This wave reads code and, later, files "
            "already on this machine"
        ),
        "what_the_run_will_write": ot.WILL_WRITE,
        "requests": {
            "alphagenome_requests": 0,
            "money": "none: no request is sent and no key is read",
            "network": "none: the census is taken from the source tree",
        },
        "result_manifest": {
            "sources": [
                {
                    "accession": "the GenomeOS source tree itself",
                    "version": "the working tree this file was written from; the declared lines are "
                    "re-derived by check_census_matches_tree() at write time",
                },
                {
                    "accession": ot.PRIOR_EVIDENCE,
                    "version": "the committed 2026-09-27 result (commit d717b28), quoted for "
                    "the_invariant_already_measured, re-read to check the quote and not recomputed",
                },
            ],
            "inputs": entries,
            "input_count": len(entries),
            "assembly": "n/a: no coordinate is read, written or compared. The unit of this "
            "registration is a source file and a line number",
            "coordinates": "n/a: no interval is read, written or compared; the unit is a source "
            "file and a line number",
            "parameters": {
                "head_keys": list(ot.HEAD_KEYS),
                "package_scanned": "genomeos",
                "site_form": 'a .get("<head key>") or a ["<head key>"] on one line',
                "no_threshold_is_set_or_moved": "MIN_EFFECT is neither set nor read here; the bar a "
                "gene must reach stays predict_target's own",
            },
            "exclusions": [
                "scripts/: 24 head-reading entry points, named in scripts_out_of_scope. Each writes "
                "a committed result and re-running one is its own decision",
                "tests/: a test's head read is a fixture, not a consumer",
                "every module in excluded_by_name: its `predicted` spelling is not a compact "
                "element head, and the reason is given per module",
            ],
            "partitions": {
                "identity": "the modules whose output depends on which or how many genes an "
                "element targets: this wave's population",
                "invariant": "the modules that read a head and cannot move, registered as NOT to be moved",
                "record": "the modules that read a head without being a live consumer of one",
            },
            "code_cleanliness": mf.code_cleanliness(__file__, ot.OWN_CODE),
        },
    }


def main() -> int:
    ot.check_all()
    p = payload()
    ot.check_no_preview(p)
    path = save_result(RESULT, p)
    print(f"registered: {path}")
    print(f"census: {len(ot.CENSUS)} modules, {len(ot.by_class('identity'))} of class identity")
    print(f"free to move today: {len(ot.free_to_move())}")
    print(f"the run will write: {ot.WILL_WRITE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
