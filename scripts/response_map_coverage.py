# SPDX-License-Identifier: AGPL-3.0-or-later
"""Count the loci a second increment of the response map could cover, and the denominator.

    uv run --frozen python scripts/response_map_coverage.py
    uv run --frozen python scripts/response_map_coverage.py --chroms chr21 \
        --result response_map_coverage_chr21

The response map (`genomeos/response_map.py`) covers one example, and the question before a second
increment is how much more there is to cover. This script counts only. It builds no increment, adds
no locus to the map, scores nothing, fits nothing, calls no model and downloads nothing: every input
was already on disk.

**The three inputs a locus needs to be addable**, each taken from the record rather than defined here:

- *a measured perturbation* - a CRISPRi element-gene pair from the ENCODE benchmark
  (`genomeos.attribution.measured`), which is the only perturbation assay in the measured layer: the
  other three (lentiMPRA, VISTA, saturation mutagenesis) read a sequence or the bases inside it, not
  the native locus after a perturbation of it;
- *a compiled link* - a rule the compiler emits for the element. Every element in `_attributed`
  carries exactly one predicted rule (`compile_chromosome`'s own test line, `rules == elements +
  measured rules`), so the compiled link is present for every attached element and absent exactly for
  a measured pair that attaches to no compiled element under the measured layer's own overlap rule,
  `RECIPROCAL_OVERLAP = 0.5`;
- *a reader state* - `genomeos.attribution.context_evidence.state_for`, lane-context2's own call,
  unchanged, on the rule's cell label and the element's locus. `not_assessable` is not a reader state:
  it is the label saying no reading could be taken, with the reason.

**Independent loci.** The grouping is imported from `genomeos.attribution.cell2`, where it was
registered, and is not restated here: `cell2.INDEPENDENT_LOCUS_RULE` is carried into the result word
for word. It is an operational grouping for deciding whether a population holds enough distinct
places in the genome to be worth measuring, **not established biological independence**, and it pools
cell types: `cell2.group` joins two keys on a shared measured gene or on 1 Mb proximity on one
chromosome and never looks at the cell. Every count below is therefore genome-wide and cell-pooled.

**Nothing here is a verdict on a rule.** `open_in_reader` is not validation, `not_open_in_reader` is
not the element being closed, and a count of addable loci says nothing about whether any one of them
would turn out to hold an explanation.
"""

from __future__ import annotations

import argparse
import bisect
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos import response_map as rm  # noqa: E402
from genomeos.attribution import cell2  # noqa: E402
from genomeos.attribution import compile as cp  # noqa: E402
from genomeos.attribution import context_evidence as ce  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.genome import reader  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

RESULT = "response_map_coverage"
CHROMS = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y"])
ROOT = Path(__file__).resolve().parents[1]
#: This lane's own files. Everything else uncommitted in this shared checkout belongs to another lane.
#: This script, as the entry whose transitive import closure is the counting path of its results
#: (genomeos.manifest.counting_path, the one implementation). Named rather than derived from
#: __file__ so the closure is the same however the script is invoked.
ENTRY = "scripts/response_map_coverage.py"

OWN_CODE = ("scripts/response_map_coverage.py", "tests/test_response_map_coverage.py")

#: The result files the compiler reads per chromosome, as globs, so the manifest names its inputs
#: without a hand-written list of several hundred paths going stale. Recorded one file per entry
#: rather than as groups, because `scripts/manifest_rebuild.py` resolves an input by its `path` and a
#: grouped entry carries a label there, so a grouped input would read as absent without being checked.
INPUT_GLOBS = (
    "budget_chr*.json",
    "budget_axes_chr*.json",
    "variation_chr*.json",
    "duplication_chr*.json",
    "domains_chr*.json",
    "unknown_chr*.json",
    "reader_*_chr*.json",
    "ccres_chr*.bed.gz",
    "enhancer_targets_chr*.json",
    "enhancer_targets_all_chr*.json",
    "constrained_targets_chr*.json",
)

STATUS_ADDABLE = "all_three_present"
STATUS_NO_READER = "no_reader_state"
STATUS_NO_RULE = "no_measured_rule_the_screen_measured_no_regulation"
STATUS_NO_LINK = "no_compiled_link_attaches_to_no_compiled_element"
STATUSES = (STATUS_ADDABLE, STATUS_NO_READER, STATUS_NO_RULE, STATUS_NO_LINK)
STATUS_ORDER = {s: i for i, s in enumerate(STATUSES)}

#: Where a locus sits when its pairs do not agree: the best status any one of them reaches. A locus is
#: a group of pairs, so a single status per locus needs a stated rule.
LOCUS_STATUS_RULE = (
    "a locus takes the best status any one of its measured pairs reaches, in the order "
    f"{' > '.join(STATUSES)}. A locus counted as addable may hold pairs that are not, and the four "
    "locus counts therefore sum to the total while the pair counts inside them do not partition a "
    "locus"
)

PERTURBATION_ASSAY = "crispri"
REPORTER_ASSAYS = tuple(a for a in ms.ASSAYS if a != PERTURBATION_ASSAY)


def attaches_to(
    pair: ms.CrispriPair, starts: list[int], elements: list[dict[str, Any]]
) -> dict[str, Any] | None:
    """The first compiled element the measured layer's overlap rule makes this pair a measurement of.

    The scan window is `measured.REACH` wide on each side, which is the longest measured interval any
    assay holds, so no element that could meet the rule lies outside it.
    """
    lo = bisect.bisect_left(starts, pair.start - ms.REACH)
    hi = bisect.bisect_right(starts, pair.end + ms.REACH)
    for e in elements[lo:hi]:
        if ms.measures(e["start"], e["end"], pair.start, pair.end):
            return e
    return None


def map_now(root: Path) -> dict[str, Any]:
    """The loci the response map covers now, by the map's own definition of a covered locus."""
    payload = rm.view(root)
    entities = payload["entities"]
    rows = []
    for a in payload["assertions"]:
        if not str(a.get("basis", "")).startswith("measured perturbation effect"):
            continue
        iv = entities[a["subject"]]
        rows.append(
            {
                "assertion": a["id"],
                "cell": a["context"]["cell"],
                "chrom": iv["chrom"],
                "start": int(iv["start"]),
                "end": int(iv["end"]),
                "gene": entities[a["object"]]["symbol"],
            }
        )
    keys = [cell2.LocusKey(r["cell"], r["chrom"], r["start"], r["end"], r["gene"]) for r in rows]
    return {
        "example": payload["example"],
        "examples_the_map_holds": list(rm.EXAMPLES),
        "increment": payload["increment"],
        "covered_locus_is": (
            "the map's own unit is an example, and increment 1 holds one: the beta-like globin genes "
            f"{', '.join(rm.LOCUS_GENES)} on {rm.LOCUS} in {rm.CELL}, selected by "
            "`response_map.WHY['selection_rule']`. Its measured perturbation assertions are counted "
            "here as loci under the imported convention so that this figure and the addable one are "
            "the same kind of number"
        ),
        "selection_rule": rm.WHY["selection_rule"],
        "measured_perturbation_assertions": len(rows),
        "independent_loci": cell2.count_loci(keys),
        "assertions": rows,
        "keys": [list(k) for k in keys],
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", default=",".join(CHROMS))
    ap.add_argument(
        "--result",
        default=RESULT,
        help=(
            "the result name to write. A run over fewer than all 24 chromosomes writes under a name "
            "of its own, so a partial count can never overwrite the genome-wide one"
        ),
    )
    args = ap.parse_args()
    chroms = [c for c in args.chroms.split(",") if c]
    name = args.result
    if name == RESULT and tuple(chroms) != CHROMS:
        raise SystemExit(
            f"{RESULT} is the genome-wide count; a run over {len(chroms)} chromosomes needs --result"
        )

    started = time.time()
    readers = ce.Readers()
    table = ce.mapping()

    pair_status: Counter[str] = Counter()
    reason: Counter[str] = Counter()
    state_of_addable: Counter[str] = Counter()
    per_chrom: dict[str, dict[str, Any]] = {}
    keys: list[cell2.LocusKey] = []  # one per measured pair, in the order the statuses were taken
    statuses: list[str] = []
    invalid_rows = 0
    compiled_elements = 0
    attached_elements = 0
    perturbation_elements = 0
    reporter_only_elements = 0
    assay_elements: Counter[str] = Counter()
    measured_rule_links = 0

    for chrom in chroms:
        elements = cp._attributed(chrom, cp.RESULTS_DIR)
        starts = [e["start"] for e in elements]
        _layer, measured_rows = cp._measured_rows(chrom, elements, cp.RESULTS_DIR)
        compiled_elements += len(elements)
        attached_elements += len(measured_rows)
        rule_of: dict[tuple[str, str, str], dict[str, Any]] = {}
        for r in measured_rows:
            for a in r["assays"]:
                assay_elements[a] += 1
            if PERTURBATION_ASSAY in r["assays"]:
                perturbation_elements += 1
            else:
                reporter_only_elements += 1
            for gene, _action, _strength, cell, _split in ms.rule_links(r):
                rule_of[(r["id"], gene, cell)] = r
                measured_rule_links += 1

        pairs, invalid = ms.load_crispri(chrom)
        invalid_rows += invalid
        here: Counter[str] = Counter()
        for p in pairs:
            e = attaches_to(p, starts, elements)
            if e is None:
                status, locus = STATUS_NO_LINK, (p.chrom, p.start, p.end)
            else:
                locus = (chrom, e["start"], e["end"])
                r = rule_of.get((e["id"], p.gene, p.cell))
                if r is None:
                    status = STATUS_NO_RULE
                else:
                    raw = ce.state_for(cp.context(p.cell), chrom, r["start"], r["end"], readers, table)
                    state, rest = ce.parse_value(raw)
                    if state == ce.STATE_NOT_ASSESSABLE:
                        status = STATUS_NO_READER
                        reason[rest[0]] += 1
                    else:
                        status = STATUS_ADDABLE
                        state_of_addable[state] += 1
            keys.append(cell2.LocusKey(p.cell, locus[0], locus[1], locus[2], p.gene))
            statuses.append(status)
            pair_status[status] += 1
            here[status] += 1
        per_chrom[chrom] = {
            "compiled_elements": len(elements),
            "elements_with_an_attached_measurement": len(measured_rows),
            "measured_perturbation_pairs": len(pairs),
            **{s: here[s] for s in STATUSES},
        }
        print(
            f"  {chrom}: {len(pairs):,} pairs, {here[STATUS_ADDABLE]:,} with all three"
            f" ({len(elements):,} compiled elements)",
            flush=True,
        )

    # --- loci: one grouping over the measured pairs and the map's own keys together --------------
    covered = map_now(ROOT)
    map_keys = [cell2.LocusKey(*k) for k in covered["keys"]]
    groups = cell2.group([*keys, *map_keys])
    mine, theirs = groups[: len(keys)], set(groups[len(keys) :])
    best: dict[int, str] = {}
    for g, s in zip(mine, statuses, strict=True):
        if g not in best or STATUS_ORDER[s] < STATUS_ORDER[best[g]]:
            best[g] = s
    loci_total = len(best)
    loci_by_status = Counter(best.values())
    addable_groups = {g for g, s in best.items() if s == STATUS_ADDABLE}
    already = len(addable_groups & theirs)

    payload: dict[str, Any] = {
        "question": (
            "how many loci could a second increment of the response map cover, out of how many the "
            "measured layer holds, and which of the three inputs is the limit?"
        ),
        "counting_only": (
            "no increment is built, no locus is added to the map, nothing is scored and nothing is "
            "fitted. 0 model requests, no downloads, no money: every input was already on disk"
        ),
        "chromosomes": chroms,
        "scope": (
            "every CRISPRi pair of the measured layer on all 24 chromosomes"
            if tuple(chroms) == CHROMS
            else f"the CRISPRi pairs of {len(chroms)} of 24 chromosomes: {', '.join(chroms)}"
        ),
        "map_now": covered,
        "denominator": {
            "population": (
                "every element-gene pair of the ENCODE CRISPRi benchmark on these chromosomes that "
                "the measured layer parses as valid, both arms (training and held-out) and every cell "
                "type it holds: this is the measured layer's perturbation assay, the only one of its "
                "four that perturbs the native locus"
            ),
            "measured_perturbation_pairs": sum(pair_status.values()),
            "invalid_rows_the_parser_rejected": invalid_rows,
            "independent_loci": loci_total,
            "independent_locus_rule": cell2.INDEPENDENT_LOCUS_RULE,
            "not_biological_independence": (
                "the grouping is operational and pools cell types: `cell2.group` joins two keys on a "
                "shared measured gene or on 1 Mb proximity on one chromosome and never reads the "
                "cell. It is not established biological independence and nothing may cite it as one"
            ),
            "compiled_elements": compiled_elements,
            "elements_with_an_attached_measurement": attached_elements,
            "elements_with_an_attached_perturbation": perturbation_elements,
            "elements_with_a_reporter_or_base_assay_only": reporter_only_elements,
            "elements_by_assay": dict(sorted(assay_elements.items())),
            "reporter_and_base_assays_are_not_perturbations": (
                f"{', '.join(REPORTER_ASSAYS)} read a sequence, or the bases inside it, and not the "
                "native locus after a perturbation of it, so an element carrying only those has no "
                "measured perturbation for this count"
            ),
            "measured_rule_links": measured_rule_links,
            "measured_rule_link_unit": list(ms.RULE_UNIT),
        },
        "by_status": {
            "pairs": {s: pair_status[s] for s in STATUSES},
            "independent_loci": {s: loci_by_status[s] for s in STATUSES},
            "locus_status_rule": LOCUS_STATUS_RULE,
        },
        "addable": {
            "population": (
                "loci holding at least one measured pair with all three inputs present: a measured "
                "CRISPRi perturbation, a compiled rule for the element it attaches to, and a reader "
                "state that is not `not_assessable`"
            ),
            "independent_loci": loci_by_status[STATUS_ADDABLE],
            "independent_loci_the_map_already_covers": already,
            "independent_loci_not_yet_covered": loci_by_status[STATUS_ADDABLE] - already,
            "share_of_the_measured_layer": (
                round(loci_by_status[STATUS_ADDABLE] / loci_total, 4) if loci_total else None
            ),
            "reader_state_of_the_addable_pairs": dict(sorted(state_of_addable.items())),
            "a_reader_state_is_not_a_verdict": (
                "open_in_reader is never validation, support or confirmation of a rule, and "
                "not_open_in_reader is never the element being closed or the rule contradicted: both "
                "count equally as a reader state here, which is what makes the locus addable"
            ),
            "addable_is_not_informative": (
                "a locus with all three inputs is one the map could arrange; it is not a locus where "
                "an explanation is known to exist, and this count says nothing about what any one of "
                "them would show"
            ),
        },
        "not_addable_by_what_is_missing": {
            STATUS_NO_LINK: {
                "pairs": pair_status[STATUS_NO_LINK],
                "independent_loci": loci_by_status[STATUS_NO_LINK],
                "meaning": (
                    "the pair's tested interval meets the measured layer's overlap rule "
                    f"(RECIPROCAL_OVERLAP = {ms.RECIPROCAL_OVERLAP}) against no compiled element, so "
                    "there is no compiled rule for it to be read beside"
                ),
            },
            STATUS_NO_RULE: {
                "pairs": pair_status[STATUS_NO_RULE],
                "independent_loci": loci_by_status[STATUS_NO_RULE],
                "meaning": (
                    "the pair attaches to a compiled element, which carries a compiled rule, but the "
                    "screen measured no regulation of that gene in that cell, so the compiler writes "
                    "no measured rule for the pair. A measured null is a measurement and is carried "
                    "in the element's `_measured` block; what is missing is a measured perturbation "
                    "*link*, which is what the map's measured perturbation assertions are made of"
                ),
            },
            STATUS_NO_READER: {
                "pairs": pair_status[STATUS_NO_READER],
                "independent_loci": loci_by_status[STATUS_NO_READER],
                "by_reason": {r: reason[r] for r in ce.NOT_ASSESSABLE_REASONS},
                "meaning": (
                    "the pair has a measured rule and a compiled link, but no reading could be taken "
                    "for the rule's cell at that locus; the reason is one of the two registered ones"
                ),
            },
        },
        "reader_state_source": (
            "genomeos.attribution.context_evidence.state_for, lane-context2's own call with its own "
            "registered threshold, applied unchanged to this population. No new cut-off is "
            "introduced, and the genome-wide census in data/results/context_evidence.json is not "
            "re-derived: that result holds counts per state, per cell and per chromosome and not a "
            "state per element, so a state had to be read per candidate from the same function"
        ),
        "openness_call": ce.OPENNESS_CALL,
        "no_interval_is_reported": (
            "these are exhaustive counts over every measured pair the layer holds, not an estimate "
            "from a sample, so no confidence interval is computed and none would mean anything here"
        ),
        "per_chromosome": per_chrom,
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
        "alphagenome_requests": 0,
        "money": "none: every input was already on disk",
        "seconds": round(time.time() - started, 1),
    }

    inputs = [mf.input_entry(ce.TRACK_METADATA, partition=None)]
    for f in ms.CRISPRI_FILES:
        p = ms.CRISPRI_KNOWLEDGE / f
        if p.exists():
            inputs.append(mf.input_entry(p, partition=ms.CRISPRI_SPLIT_OF[f]))
    seen = {str(i["path"]) for i in inputs}
    paths = [
        p
        for cell in ce.READER_TERMS
        for chrom in chroms
        if (p := RESULTS_DIR / reader.peaks_path(cell, chrom).name).exists()
    ]
    for glob in INPUT_GLOBS:
        for chrom in chroms:
            paths += sorted(RESULTS_DIR.glob(glob.replace("chr*", chrom)))
    for p in paths:
        if p.as_posix() not in seen:
            seen.add(p.as_posix())
            inputs.append(mf.input_entry(p, partition=None))

    payload["result_manifest"] = {
        "sources": [
            {
                "accession": rm.CRISPRI_ACCESSION,
                "version": "the two local benchmark tables under data/knowledge/crispri",
            },
            {
                "accession": "ENCODE DNase-seq narrowPeak, GRCh38, released, 13 biosamples",
                "version": "as reader v1 cached them under data/results/dnase_*_chr*.bed.gz",
            },
            {"accession": "AlphaGenome output track metadata", "version": str(ce.TRACK_METADATA)},
            {
                "accession": "the compiled non-coding programs of the 24 chromosomes",
                "version": (
                    "the compiler's own element and measured-layer enumeration, run in this process "
                    "from the results on disk, not read back from a .bio file"
                ),
            },
        ],
        "inputs": inputs,
        "inputs_are_recorded_one_file_per_entry": (
            "scripts/manifest_rebuild.py resolves an input by its `path`, and a grouped entry carries "
            "a label there, so a grouped input would read as absent without any file being checked"
        ),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "independent_locus_span": cell2.INDEPENDENT_LOCUS_SPAN,
            "reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
            "openness_call": ce.OPENNESS_CALL,
            "states": list(ce.STATES),
            "not_assessable_reasons": list(ce.NOT_ASSESSABLE_REASONS),
            "chromosomes": chroms,
        },
        "exclusions": [
            "no pair is excluded: the denominator is every CRISPRi pair the measured layer parses as "
            "valid on these chromosomes, both arms and every cell type",
            "the three reporter and base-level assays of the measured layer are counted but are not "
            "measured perturbations, and are named as such rather than dropped silently",
            "no new cut-off is introduced: the overlap rule, the locus span and the openness call are "
            "all imported from where they were registered",
        ],
        "partitions": {
            "by_status": "which of the three inputs a pair, and a locus, has",
            "per_chromosome": "the counts of each chromosome, summing to the whole",
        },
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }
    path = save_result(name, payload)
    print(f"{name}: {path}")
    print(f"  map now: {covered['independent_loci']} independent loci ({covered['example']})")
    print(f"  measured layer: {loci_total:,} independent loci over {sum(pair_status.values()):,} pairs")
    for s in STATUSES:
        print(f"  {s}: {loci_by_status[s]:,} loci, {pair_status[s]:,} pairs")
    print(f"  addable loci the map already covers: {already}")


if __name__ == "__main__":
    main()
