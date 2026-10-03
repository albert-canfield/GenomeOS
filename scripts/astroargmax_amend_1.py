# SPDX-License-Identifier: AGPL-3.0-or-later
"""Amendment 1 to item (h)'s registration: THE UNIT IS AN ELEMENT.

    uv run --frozen python scripts/astroargmax_amend_1.py

Writes data/results/astroargmax_registration_amendment_1.json, committed ALONE and ADDITIVE.
data/results/astroargmax_registration.json is NOT rewritten: it stands as committed at b5fd1b6 and
this amendment names and checks the sha256 of that committed blob. A registration edited in place
cannot be told from one that was always right.

THE DEFECT the supervisor's review found before any request was sent, verified here against the
code: the registration left the unit undefined, and as written it compared two populations. Its
question spoke of a GENE ROW's argmax, while control (a) counts COMPILED RULES -- one label per
ELEMENT, `compile.py:851 label = context(pc.get("tissue"))` on the `predicted_coding` head chosen by
`enhancer_target.predict_target` -- and comparison (b) reads that same per-element tissue.

THE THREE PARTS, all required and all here:
  1. ONE LABEL PER ELEMENT, by the IMPORTED `predict_target` -> tissue -> `compile.context`.
  2. Elements with NO predicted coding target EXCLUDED, and the excluded count REPORTED beside the
     denominator.
  3. The all-gene-row argmax demoted to a LABELLED SECONDARY that decides nothing.

AND IT IS BLIND: not one of the 1,232 requests has been sent, no astrocyte answer exists, and the
sender still refuses on two of Albert's five clauses, so the outcome is unobtainable and not merely
unseen. 0 AlphaGenome requests, no network to the service, no key read, no money.
"""

from __future__ import annotations

import subprocess
import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import astroargmax as aa  # noqa: E402

RESULT = "astroargmax_registration_amendment_1"
REGISTRATION = Path("data/results/astroargmax_registration.json")

OWN_CODE = (
    "genomeos/attribution/astroargmax.py",
    "scripts/astroargmax_register.py",
    "scripts/astroargmax_amend_1.py",
    "tests/test_astroargmax.py",
)


def amended_blob() -> dict[str, Any]:
    """The committed registration this amends, checked by the sha256 of its blob at that commit.

    Read from git and not from the working copy: the point of an amendment is that the thing it
    corrects has not moved, and only the object store can say so.
    """
    import hashlib

    root = Path(__file__).resolve().parents[1]
    shown = subprocess.run(
        ["git", "show", f"{aa.AMENDS_COMMIT}:{REGISTRATION.as_posix()}"],
        cwd=str(root),
        capture_output=True,
        check=True,
    )
    digest = hashlib.sha256(shown.stdout).hexdigest()
    if digest != aa.AMENDS_BLOB_SHA256:
        raise ValueError(
            f"the registration committed at {aa.AMENDS_COMMIT} hashes {digest} and this amendment "
            f"names {aa.AMENDS_BLOB_SHA256}. An amendment that cannot identify what it amends "
            "amends nothing"
        )
    return {
        "result": aa.AMENDS,
        "path": REGISTRATION.as_posix(),
        "commit": aa.AMENDS_COMMIT,
        "blob_sha256": digest,
        "bytes": len(shown.stdout),
        "not_edited": aa.AMENDMENT_1_THE_REGISTRATION_IS_NOT_EDITED,
    }


def the_evidence_in_the_committed_file() -> dict[str, int]:
    """Counted on the COMMITTED bytes, so the defect is shown rather than asserted."""
    root = Path(__file__).resolve().parents[1]
    text = subprocess.run(
        ["git", "show", f"{aa.AMENDS_COMMIT}:{REGISTRATION.as_posix()}"],
        cwd=str(root),
        capture_output=True,
        check=True,
    ).stdout.decode("utf-8")
    return {
        term: text.count(term) for term in ("gene row", "predict_target", "predicted_coding", "per element")
    }


def inputs() -> list[dict[str, Any]]:
    """The one file this amendment reads: the registration it corrects, as committed."""
    out: list[dict[str, Any]] = []
    if REGISTRATION.exists():
        out.append(
            mf.input_entry(
                REGISTRATION,
                partition="the registration this amends. Read for the four term counts that EVIDENCE "
                "the defect, and identified by the sha256 of its blob at b5fd1b6. Not rewritten",
            )
        )
    return out


def payload() -> dict[str, Any]:
    from genomeos.attribution.compile import context
    from genomeos.predict.enhancer_target import MIN_EFFECT, predict_target

    entries = inputs()
    return {
        "result": RESULT,
        "date": date.today().isoformat(),
        "lane": "lane-astroargmax",
        "amends": aa.AMENDS,
        "amended_registration": amended_blob(),
        "status": (
            "AN AMENDMENT AND NOTHING ELSE, committed ALONE and ADDITIVE. No count, no result, no "
            "answer, no rate, no interval, no verdict. 0 AlphaGenome requests, no network to the "
            "service, no key read, no money"
        ),
        "why": aa.AMENDMENT_1,
        "the_evidence_counted_on_the_committed_bytes": the_evidence_in_the_committed_file(),
        "the_label_path_in_the_code": {
            "compiled_rule_label": "genomeos/attribution/compile.py: `label = context(pc.get("
            '"tissue"))`, where `pc` is the element\'s `predicted_coding` head',
            "head_chosen_by": "genomeos.predict.enhancer_target.predict_target(rows, "
            "min_effect=MIN_EFFECT): ONE dict per element, the gene and the tissue chosen JOINTLY, "
            "or None when no gene moves by MIN_EFFECT",
            "label_function": "genomeos.attribution.compile.context",
            "min_effect": float(MIN_EFFECT),
            "min_effect_is_imported_not_copied": "astroargmax.min_effect() reads it from "
            "enhancer_target. No literal 0.1 exists in this lane's code, because a copied literal "
            "and a re-implemented rule are the same defect",
            "verified_here_not_relayed": "both functions are imported in this script and their "
            "identities are recorded below, so the amendment rests on the code and not on a review's "
            "account of it",
            "imported_identities": {
                "predict_target": f"{predict_target.__module__}.{predict_target.__qualname__}",
                "context": f"{context.__module__}.{context.__qualname__}",
            },
        },
        # ------------------------------------------------------------- the three parts, all required
        "part_1_the_unit": aa.AMENDMENT_1_PART_1,
        "part_2_the_denominator": aa.AMENDMENT_1_PART_2,
        "part_3_the_secondary": aa.AMENDMENT_1_PART_3,
        "the_one_place_the_unit_is_computed": {
            "function": "genomeos.attribution.astroargmax.element_label",
            "it_holds_no_argmax_of_its_own": "it calls predict_target and compile.context and "
            "returns None when the head is None. A test substitutes predict_target and shows the "
            "label follows it, which proves delegation rather than agreement",
            "population_builder": "genomeos.attribution.astroargmax.label_population, which returns "
            "elements_with_a_target as the denominator and elements_with_no_target beside it",
        },
        "secondary_all_gene_rows": aa.SECONDARY_ALL_GENE_ROWS,
        "one_population_definition_for_all_three": (
            "the astrocyte arm, control (a) and comparison (b) are now the SAME population "
            "definition: one label per element, from its predicted coding target, by these two "
            "functions. That is the property the registration lacked and the only thing this "
            "amendment adds"
        ),
        # -------------------------------------------------------------------- what does not move
        "what_does_not_move": aa.AMENDMENT_1_WHAT_DOES_NOT_MOVE,
        "bands_unchanged": {
            "ratio_enriched_floor": aa.RATIO_ENRICHED_FLOOR,
            "absolute_level_floor": aa.ABSOLUTE_LEVEL_FLOOR,
            "equivalence_band_relative": aa.EQUIVALENCE_BAND,
            "ratio_depleted_ceiling": aa.RATIO_DEPLETED_CEILING,
            "minimum_clusters_for_an_interval": aa.MIN_CLUSTERS,
            "bootstrap_draws": aa.BOOTSTRAP_DRAWS,
            "bootstrap_seed": aa.BOOTSTRAP_SEED,
            "minimum_distinct_resamples": aa.MIN_DISTINCT_RESAMPLES,
            "label_set_CNS_PRIMARY": [c for c, _ in aa.REGISTERED_SET],
            "note": "byte-identical to the committed registration. This amendment fixes WHAT IS "
            "COUNTED ONCE, not how much of it is needed to read anything",
        },
        # ------------------------------------------------------------------------- blindness
        "this_amendment_is_blind": aa.AMENDMENT_1_IS_BLIND,
        "what_was_read_of_the_outcome_when_this_was_written": (
            "nothing, because nothing exists. 0 of the 1,232 requests sent, 0 astrocyte answers, 0 "
            "argmax labels of either arm, 0 rates, 0 intervals. Not a disclosure of restraint: a "
            "statement that the quantity has no bytes anywhere on this machine"
        ),
        "and_it_could_not_have_been_obtained": (
            "the sender refuses on two of Albert's five clauses -- today's CI green and the "
            "supervisor's 'dry run reviewed' -- so there was no path by which a figure could have "
            "been looked at before this was written. Unobtainable, not merely unseen"
        ),
        "which_kind_of_amendment_this_is": (
            "THE BLIND KIND. The other kind is on record: a lane's amendment on 2026-10-02 was "
            "written after its first count had been seen, and it cost that result its clean "
            "negative. This one is the first kind and the evidence is in the repository -- the "
            "registration's commit, this amendment's commit, and an empty ledger between them -- "
            "rather than in this sentence"
        ),
        "the_registration_is_not_edited": aa.AMENDMENT_1_THE_REGISTRATION_IS_NOT_EDITED,
        "this_file_does_not_authorise_the_run": aa.NOT_AN_AUTHORISATION,
        "alphagenome_requests": 0,
        "network_requests": 0,
        "money": "none: no request is sent and no key is read",
        "result_manifest": {
            "sources": [
                {
                    "accession": "data/results/astroargmax_registration.json, the registration this amends",
                    "version": f"the blob committed at {aa.AMENDS_COMMIT}, sha256 "
                    f"{aa.AMENDS_BLOB_SHA256}; read from the object store and not from the working "
                    "copy, because an amendment has to show that what it corrects has not moved",
                }
            ],
            "inputs": entries,
            "input_count": len(entries),
            "assembly": "n/a: this amendment reads no coordinate and no sequence",
            "coordinates": "n/a: no interval is read. The amendment is about which label belongs to "
            "an element, not about where the element is",
            "parameters": {
                "unit": "one label per paid-answer ELEMENT",
                "label_path": "enhancer_target.predict_target -> head['tissue'] -> compile.context, "
                "all IMPORTED",
                "min_effect": float(MIN_EFFECT),
                "denominator": "elements_with_a_target",
                "excluded_and_counted": "elements_with_no_target",
                "secondary_that_decides_nothing": "the all-gene-row argmax, denominator gene rows",
                "no_band_or_threshold_moves": "every number of the committed registration is "
                "reproduced above and is byte-identical. No threshold, floor, cap, pin, test or "
                "registration of any lane is moved, weakened, exempted or allowlisted",
            },
            "exclusions": [
                "every element whose predict_target returns None: no gene moves by MIN_EFFECT, so it "
                "emits no compiled rule and is in neither the control's numerator nor its "
                "denominator. COUNTED and reported, never silently dropped",
                "the all-gene-row argmax, from the reading: reported as a labelled secondary only",
                "every astrocyte answer: none exists and none is read",
            ],
            "partitions": {
                "elements_with_a_target": "the registered denominator of every arm",
                "elements_with_no_target": "reported beside it, summing to the answers read",
                "gene_rows": "the secondary's denominator, and never the reading's",
            },
            "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
        },
    }


def main() -> int:
    from genomeos.results import save_result

    body = payload()
    path = save_result(RESULT, body)
    print(f"amended: {path}")
    a = body["amended_registration"]
    print(
        f"  amends {a['result']} at {a['commit']}, blob {a['blob_sha256'][:12]}..., "
        f"{a['bytes']} bytes, NOT edited"
    )
    print(f"  evidence on the committed bytes: {body['the_evidence_counted_on_the_committed_bytes']}")
    print(f"  unit: one label per ELEMENT, min_effect {body['the_label_path_in_the_code']['min_effect']}")
    print("  bands unchanged: " + ", ".join(f"{k}={v}" for k, v in list(body["bands_unchanged"].items())[:4]))
    print("  blind: 0 of 1232 requests sent, 0 astrocyte answers, and 2 of 5 clauses still refuse")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
