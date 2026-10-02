# SPDX-License-Identifier: AGPL-3.0-or-later
"""GTEx v8 re-distilled against the CRISPRi element set, to turn "unassessed" into an answer.

    uv run python scripts/eqtl_crispri_frame.py --register   # the frame and the rule, before any byte
    uv run python scripts/eqtl_crispri_frame.py --distil     # stream 1.56 GB once, keep nothing else
    uv run python scripts/eqtl_crispri_frame.py --count      # what the new frame assessed, and found

WHY. `response_map_increment4_count` found that the existing GTEx distillation was indexed against
three sampled element sets, so of the 1,505 compiled elements a CRISPRi pair attaches to only 61
lay inside its frame at all and 1,444 were UNASSESSED -- no eQTL could have been seen at them
whatever is really there. That is not an absence of eQTLs, and the only thing that turns it into an
answer is a distillation indexed against the CRISPRi element set itself. This does that.

WHAT IS REGISTERED BEFORE ANY BYTE IS STREAMED, and all of it is imported rather than chosen:
  * the frame: every compiled element a valid CRISPRi pair of the benchmark attaches to, by
    `response_map2.candidates`'s own attachment, widened by MARGIN bases each side;
  * MARGIN = 500, `scripts/eqtl_targets.MARGIN`, the margin the existing retained set was distilled
    under, imported via `response_map_increment3_count.EQTL_MARGIN`;
  * what counts as a retained hit: `genomeos.attribution.eqtl.distil`'s own rule, unchanged -- a
    pair of GTEx's OWN significant variant-gene set whose variant position falls inside an indexed
    interval -- written in `eqtl.HIT_COLUMNS` order so `eqtl.load_hits` reads it;
  * the floor the deepening arm must clear to be worth reporting as a candidate again, below.

NOTHING OVERWRITES THE EXISTING DISTILLATION. The hits go to `data/knowledge/gtex_crispri`, not to
`data/knowledge/gtex`, because the latter is pinned by sha256 in the manifests of increment 2, of
increment 3, of its count and of `eqtl_targets`, and re-indexing it in place would break every one
of them while every number in them still looked plausible.

NO MODEL REQUEST and no paid call: GTEx v8 is a public bucket. 0 AlphaGenome requests.

CLEARING THE FLOOR DOES NOT AUTHORISE A BUILD. It authorises a report. Counting answerability in
the same breath as building on it is the error the direction test avoided by counting it
separately, which is how it noticed that only 21 of 48 links were answerable in their own cell.
"""

from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import response_map_increment3_count as c3  # noqa: E402

from genomeos import manifest as mf  # noqa: E402
from genomeos import response_map2 as rm2  # noqa: E402
from genomeos.attribution import cell2, eqtl  # noqa: E402
from genomeos.attribution import context_evidence as ce  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.jobs import heartbeat  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

REGISTRATION = "eqtl_crispri_frame_registration"
RESULT = "eqtl_crispri_frame"
ENTRY = "scripts/eqtl_crispri_frame.py"
OWN_CODE = (
    "scripts/response_map_increment4_count.py",
    "scripts/eqtl_crispri_frame.py",
    "tests/test_response_map_increment4.py",
)

#: A directory of its own. `data/knowledge/gtex` is a sha256-pinned input of four committed results
#: and is not touched.
KNOWLEDGE = Path("data/knowledge/gtex_crispri")

#: `scripts/eqtl_targets.MARGIN`, imported rather than restated, via increment 3's count.
MARGIN = c3.EQTL_MARGIN

#: The floor, fixed here before any byte is streamed. Increment 3 built 6 addable loci carrying a
#: native-locus eQTL on an element one of its chains stands on. A deepening arm is worth proposing
#: again when the new frame puts that figure at or above this, which is more than triple the 6 and
#: about a fifth of the 110 addable loci. Below it, the arm stays refused and nothing is proposed.
FLOOR_ADDABLE_LOCI_ON_A_SHOWN_ELEMENT = 20

#: What a retained hit is. Quoted from the code that decides it, not restated in this lane's words.
RETAINED_HIT_RULE = (
    "a pair of GTEx v8's OWN single-tissue significant variant-gene set whose variant position "
    "falls inside one of the indexed intervals, kept by `genomeos.attribution.eqtl.distil` "
    "unchanged and written in `eqtl.HIT_COLUMNS` order so that `eqtl.load_hits` reads it. No "
    "p-value, effect size or tissue is filtered here: GTEx's significance call is the one used, and "
    "this lane sets no cut-off of its own"
)

FRAME_RULE = (
    "every compiled element that a valid CRISPRi pair of the ENCODE benchmark attaches to, by "
    "`response_map2.candidates`'s own attachment at `measured.RECIPROCAL_OVERLAP = 0.5`, over all "
    f"24 chromosomes and every cell type, each widened by MARGIN = {MARGIN} bases on each side"
)

ASSESSED_BY_CONSTRUCTION = (
    "every element of the frame is assessed BY CONSTRUCTION, because the frame is built from them. "
    "That is not a result and may not be reported as one: the figure that is a measurement is how "
    "many of them carry a retained hit, and the two are kept apart below for exactly that reason"
)


class RefusedError(ValueError):
    """Raised instead of streaming or saving when a precondition does not hold."""


def frame(root: Path = ROOT) -> dict[str, Any]:
    """The CRISPRi element set, its intervals, and the locus grouping it belongs to."""
    cands, counted = rm2.candidates(list(rm2.CHROMS), RESULTS_DIR, ce.Readers(), ce.mapping())
    covered = rm2.map1_keys(root)
    groups = cell2.group([c.key for c in cands] + covered)
    mine = groups[: len(cands)]
    best: dict[int, str] = {}
    for g, c in zip(mine, cands, strict=True):
        if g not in best or rm2.STATUS_ORDER[c.status] < rm2.STATUS_ORDER[best[g]]:
            best[g] = c.status
    addable = {g for g, s in best.items() if s == rm2.STATUS_ADDABLE}
    att_iv: dict[str, tuple[str, int, int]] = {}
    att_groups: dict[str, set[int]] = collections.defaultdict(set)
    shown: dict[int, set[str]] = collections.defaultdict(set)
    shown_iv: dict[str, tuple[str, int, int]] = {}
    for g, c in zip(mine, cands, strict=True):
        if c.element is None:
            continue
        eid = c.element["id"]
        att_iv[eid] = (c.pair.chrom, c.element["start"], c.element["end"])
        att_groups[eid].add(g)
        if g in addable and c.status == rm2.STATUS_ADDABLE:
            shown[g].add(eid)
            shown_iv[eid] = (c.pair.chrom, c.element["start"], c.element["end"])
    return {
        "cands": cands,
        "counted": counted,
        "mine": mine,
        "best": best,
        "addable": addable,
        "att_iv": att_iv,
        "att_groups": att_groups,
        "shown": shown,
        "shown_iv": shown_iv,
    }


def intervals(att_iv: dict[str, tuple[str, int, int]]) -> eqtl.Intervals:
    """The frame as `eqtl.Intervals`, keyed the way `eqtl_targets.py` keys them."""
    iv = eqtl.Intervals()
    seen = set()
    for _eid, (chrom, s, e) in sorted(att_iv.items()):
        key = (chrom, s, e)
        if key in seen:
            continue
        seen.add(key)
        iv.add(chrom, max(0, s - MARGIN), e + MARGIN, f"{chrom}:{s}-{e}")
    return iv.freeze()


def registration_payload(f: dict[str, Any], iv: eqtl.Intervals) -> dict[str, Any]:
    by_status = collections.Counter(f["best"].values())
    payload: dict[str, Any] = {
        "result": REGISTRATION,
        "date": "2026-10-02",
        "registered_before": "any byte of the GTEx archive was streamed",
        "question": (
            "how many of the compiled elements a CRISPRi pair attaches to carry a GTEx cis-eQTL, "
            "once the distillation is indexed against those elements instead of three sampled sets"
        ),
        "why": (
            "the existing retained set was distilled against 9,286 intervals of three sampled "
            "element sets, so only 61 of the 1,505 CRISPRi-attached compiled elements lay inside it "
            "and 1,444 were UNASSESSED. Unassessed is not absent, and only a re-distillation indexed "
            "against these elements can say which it is"
        ),
        "frame_rule": FRAME_RULE,
        "retained_hit_rule": RETAINED_HIT_RULE,
        "margin_bp": MARGIN,
        "margin_source": (
            "scripts/eqtl_targets.MARGIN, imported through "
            "scripts/response_map_increment3_count.EQTL_MARGIN, not chosen here"
        ),
        "imported_not_chosen": [
            "genomeos.attribution.eqtl.distil - what a retained hit is, and the file format",
            "genomeos.attribution.eqtl.load_hits - how it is read back",
            "genomeos.response_map2.candidates - the element a pair attaches to, and the status",
            "genomeos.attribution.cell2.group - the independent-locus grouping",
            "scripts/eqtl_targets.MARGIN = 500 - the margin",
        ],
        "frame": {
            "elements": len(f["att_iv"]),
            "intervals_indexed": iv.n,
            "chromosomes": sorted({c for c, _s, _e in f["att_iv"].values()}),
            "crispri_pairs": f["counted"]["pairs"],
            "crispri_loci": len(set(f["mine"])),
            "crispri_loci_by_best_status": dict(by_status),
            "addable_loci": len(f["addable"]),
            "shown_elements_of_the_built_chains": len(f["shown_iv"]),
        },
        "assessed_by_construction": ASSESSED_BY_CONSTRUCTION,
        "the_two_figures_kept_apart": {
            "assessed": "how many of the frame's elements could have been seen at all",
            "carrying_a_hit": "how many of them carry at least one retained hit",
            "why": (
                "conflating them is how an unassessed element becomes an element with no eQTL. The "
                "first is a property of the frame and the second is a measurement"
            ),
        },
        "floor": {
            "quantity": (
                "addable loci carrying a retained eQTL inside an element one of the built chains "
                "stands on, which is increment 3's own attachment rule and not a relaxed one"
            ),
            "floor": FLOOR_ADDABLE_LOCI_ON_A_SHOWN_ELEMENT,
            "increment_3_built": 6,
            "why_this_number": (
                "more than triple what increment 3 built, and about a fifth of the 110 addable "
                "loci, so a deepening arm clearing it would be a material share of the population "
                "rather than a handful"
            ),
            "clearing_it_does_not_authorise_a_build": (
                "it authorises a REPORT. Counting answerability in the same breath as building on "
                "it is the error the direction test avoided by counting it separately, which is how "
                "it noticed that only 21 of 48 links were answerable in their own cell"
            ),
        },
        "where_the_hits_go": KNOWLEDGE.as_posix(),
        "nothing_is_overwritten": (
            "data/knowledge/gtex is a sha256-pinned input of response_map_increment2, "
            "response_map_increment3, its count and eqtl_targets. Re-indexing it in place would "
            "break all four while every number in them still looked plausible, so this writes to a "
            "directory of its own"
        ),
        "streaming": (
            "the 1.56 GB archive is read member by member over HTTP ranges and discarded; only the "
            "retained hits reach disk, a few megabytes"
        ),
        "alphagenome_requests": 0,
        "money": "none: GTEx v8 is a public bucket and no paid call is made",
    }
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "GTEx v8 single-tissue cis-eQTLs (significant variant-gene pairs), hg38",
                "version": eqtl.GTEX_EQTL_TAR,
            },
            {
                "accession": "EngreitzLab/CRISPR_comparison (Gschwind et al. 2025)",
                "version": "the two local benchmark tables under data/knowledge/crispri",
            },
        ],
        "inputs": _inputs(),
        "inputs_opened_beside_the_declared_results": {
            "why_they_are_listed": (
                "an input is what the writer opened, not what it said it read. The audit hook in "
                "`genomeos.manifest` records every open under `data/` and `save_result` refuses a "
                "result whose traced reads exceed its declared inputs"
            ),
            "groups": {k: len(v) for k, v in c3.OPENED_BESIDE_THE_RESULTS.items()},
            "the_pointer_case": (
                "`response_map2.candidates` opens the 24 per-chromosome element tables through a "
                "result file's `elements_where` POINTER field - about 1.4 GB that no reading of the "
                "source can name"
            ),
        },
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "eqtl_margin_bp": MARGIN,
            "reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
            "independent_locus_span": cell2.INDEPENDENT_LOCUS_SPAN,
            "floor_addable_loci_on_a_shown_element": FLOOR_ADDABLE_LOCI_ON_A_SHOWN_ELEMENT,
            "chromosomes": list(rm2.CHROMS),
        },
        "exclusions": [
            "no element of the frame is excluded: it is every element a valid CRISPRi pair attaches "
            "to, over all 24 chromosomes and every cell type",
            "GTEx's own significance set is the population; no further filter is applied here",
        ],
        "partitions": {"frame": "the CRISPRi element set the archive is indexed against"},
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }
    return payload


def _inputs() -> list[dict[str, Any]]:
    """Everything this writer opens under data/, including what it reaches through a pointer field."""
    from genomeos.genome import reader

    inputs = [mf.input_entry(ce.TRACK_METADATA, partition=None)]
    for f in ms.CRISPRI_FILES:
        p = ms.CRISPRI_KNOWLEDGE / f
        if p.exists():
            inputs.append(mf.input_entry(p, partition=ms.CRISPRI_SPLIT_OF[f]))
    seen = {str(i["path"]) for i in inputs}
    paths = [
        p
        for cell in ce.READER_TERMS
        for chrom in rm2.CHROMS
        if (p := RESULTS_DIR / reader.peaks_path(cell, chrom).name).exists()
    ]
    for g in (
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
        "vista_chr*.json",
    ):
        for chrom in rm2.CHROMS:
            paths += sorted(RESULTS_DIR.glob(g.replace("chr*", chrom)))
    for group in c3.OPENED_BESIDE_THE_RESULTS.values():
        paths += [Path(x) for x in group]
    paths += sorted(KNOWLEDGE.glob("hits_*.tsv"))
    p = KNOWLEDGE / "distil_summary.json"
    if p.exists():
        paths.append(p)
    for q in paths:
        if q.as_posix() not in seen:
            seen.add(q.as_posix())
            inputs.append(mf.input_entry(q, partition=None))
    return inputs


def count_payload(f: dict[str, Any], iv: eqtl.Intervals) -> dict[str, Any]:
    """What the new frame assessed, and what it found, with the two kept apart."""
    hits = c3.eqtl_hits_by_chrom(KNOWLEDGE)
    summary_path = KNOWLEDGE / "distil_summary.json"
    distil_summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}

    att_iv: dict[str, tuple[str, int, int]] = f["att_iv"]
    with_hit: dict[str, list[dict[str, Any]]] = {}
    for eid, (chrom, s, e) in att_iv.items():
        got = c3.hits_in(hits.get(chrom, []), s, e)
        if got:
            with_hit[eid] = got

    hit_groups = {g for eid in with_hit for g in f["att_groups"][eid]}
    shown_hit = {eid: v for eid, v in with_hit.items() if eid in f["shown_iv"]}
    shown_hit_loci = {g for g, eids in f["shown"].items() if eids & set(shown_hit)}
    addable_hit_any = {g for g in hit_groups if f["best"][g] == rm2.STATUS_ADDABLE}
    records = sum(len(v) for v in with_hit.values())
    neg = sum(1 for v in with_hit.values() for h in v if h["slope"] < 0)
    pos = records - neg

    # TWO UNITS, and the first version of this block subtracted one from the other and published
    # -1,220. `increment3_count.hits_in` widens the element by EQTL_MARGIN itself, so the window it
    # attaches in is the SAME window the distillation retained under -- which is increment 3's
    # registered wording, "inside the compiled element widened by 500 bases each side". An earlier
    # draft of this comment said the attachment was the element ITSELF; that was wrong, and the
    # margin-only figure of 0 below is that identity rather than a finding. What does differ is the
    # unit: a row inside two overlapping compiled elements attaches TWICE, so `records` above is in
    # ATTACHMENTS and the three figures below are in ROWS. Only rows may be subtracted from rows.
    def row_key(chrom: str, h: dict[str, Any]) -> tuple[str, str, int, str]:
        return (h["tissue"], chrom, h["pos"], h["gene_id"])

    retained_rows = {row_key(c, h) for c, v in hits.items() for h in v}
    attached_rows = {row_key(att_iv[eid][0], h) for eid, v in with_hit.items() for h in v}
    in_the_margin_only = len(retained_rows) - len(attached_rows)

    # the before/after, taken from the committed count rather than restated
    with open(ROOT / "data/results/response_map_increment4_count.json") as fh:
        prior = json.load(fh)
    before = next(c for c in prior["candidates"] if c["candidate"] == "eqtl")[
        "it_splits_into_two_arms_and_both_are_refused"
    ]["deepening"]

    cleared = len(shown_hit_loci) >= FLOOR_ADDABLE_LOCI_ON_A_SHOWN_ELEMENT
    payload: dict[str, Any] = {
        "result": RESULT,
        "date": "2026-10-02",
        "question": (
            "how many of the compiled elements a CRISPRi pair attaches to carry a GTEx cis-eQTL, "
            "with the new frame indexed against those elements"
        ),
        "registered_at": f"data/results/{REGISTRATION}.json",
        "registered_sha256": mf.input_entry(ROOT / f"data/results/{REGISTRATION}.json")["sha256"],
        "frame_rule": FRAME_RULE,
        "retained_hit_rule": RETAINED_HIT_RULE,
        "assessed": {
            "elements_of_the_frame": len(att_iv),
            "assessed": len(att_iv),
            "unassessed": 0,
            "by_construction": ASSESSED_BY_CONSTRUCTION,
            "what_changed": (
                f"{before['inside_the_distillation_frame_at_all']} of "
                f"{before['crispri_attached_compiled_elements']} elements were inside the old frame "
                f"and {before['outside_the_frame_unassessed']} were unassessed. All "
                f"{len(att_iv)} are now assessed, and that part of the change is the frame's "
                "definition rather than a finding"
            ),
        },
        "found": {
            "elements_carrying_at_least_one_retained_hit": len(with_hit),
            "elements_assessed_and_carrying_none": len(att_iv) - len(with_hit),
            "attachments_element_by_record": records,
            "distinct_rows_attached_inside_an_element": len(attached_rows),
            "distinct_rows_retained_by_the_distillation": len(retained_rows),
            "distinct_rows_retained_in_the_margin_only": in_the_margin_only,
            "two_units_and_one_window": (
                "the attachment window is increment 3's registered one, imported and not chosen: "
                "'a retained hit whose variant position falls inside the compiled element widened "
                f"by {MARGIN} bases each side, which is the margin the retained set was distilled "
                "under (scripts/eqtl_targets.MARGIN)'. It is therefore the SAME window the "
                "distillation retained under, so distinct_rows_attached_inside_an_element equals "
                "distinct_rows_retained_by_the_distillation and the margin-only figure is 0 BY "
                "CONSTRUCTION - an identity between two imported rules, and not a finding. What "
                "does differ is the UNIT: a row inside two overlapping compiled elements attaches "
                "TWICE, so attachments_element_by_record is in ATTACHMENTS and the other three are "
                "in ROWS, and the slope split is of attachments at that same denominator. The "
                "first version of this block subtracted a row count from an attachment count and "
                "published a margin-only figure of -1,220; a smaller overlap would have made it "
                "positive and plausible, which is why the units are named here"
            ),
            "slope_negative": neg,
            "slope_positive": pos,
            "tissues": len({h["tissue"] for v in with_hit.values() for h in v}),
            "distinct_egenes": len({h["gene_id"] for v in with_hit.values() for h in v}),
            "a_count_is_not_a_measurement": (
                "the record count counts variant-gene-tissue rows and not places in the genome, and "
                "the slope split is of records rather than of loci. A cis-eQTL is an association "
                "over donors besides: linkage carries the signal some distance, so the variant "
                "inside the element need not be causal and the eGene need not be regulated by it"
            ),
        },
        "loci": {
            "crispri_loci_total": len(set(f["mine"])),
            "loci_with_an_eqtl_inside_an_attached_element": len(hit_groups),
            "loci_with_a_hit_by_best_status": dict(collections.Counter(f["best"][g] for g in hit_groups)),
            "addable_loci_with_a_hit_on_ANY_attached_element": len(addable_hit_any),
            "addable_loci_with_a_hit_on_a_SHOWN_element": len(shown_hit_loci),
            "shown_elements": len(f["shown_iv"]),
            "shown_elements_carrying_a_hit": len(shown_hit),
            "the_rule_that_counts": (
                "increment 3's own attachment, carried in its registered wording: a retained hit "
                f"whose variant position falls inside the compiled element widened by {MARGIN} "
                "bases each side - and, for the figure that counts, inside an element one of the "
                "BUILT CHAINS stands on. The ANY-element figure is reported beside it and is NOT "
                "the same quantity: 86 of the 92 the grouping offered were a different element "
                "within 1 Mb, which is the number this lane withdrew before quoting it"
            ),
        },
        "the_floor": {
            "quantity": "addable loci with a hit on a shown element",
            "floor": FLOOR_ADDABLE_LOCI_ON_A_SHOWN_ELEMENT,
            "reached": len(shown_hit_loci),
            "increment_3_built": 6,
            "cleared": cleared,
            "registered_before_any_byte_was_streamed": True,
            "what_follows": (
                "a report and nothing else. Clearing the floor authorises proposing the deepening "
                "arm again; it does not authorise building increment 4, and nothing is built here"
            )
            if cleared
            else (
                "the deepening arm stays refused, now on a measurement rather than on an "
                "unassessed frame. That is the substantive change: the refusal is no longer "
                "'nothing could have been seen' but 'this is what there is'"
            ),
        },
        "before_and_after": {
            "old_frame": {
                "intervals": before["frame_intervals"],
                "indexed_against": "three sampled element sets",
                "elements_inside_it": before["inside_the_distillation_frame_at_all"],
                "elements_unassessed": before["outside_the_frame_unassessed"],
                "elements_with_a_hit": before["elements_with_at_least_one_retained_hit"],
                "attachments_element_by_record": before["retained_records_on_them"],
                "loci_with_a_hit": before["loci_with_an_eqtl_inside_an_attached_element"],
                "addable_loci_on_a_shown_element": before["addable_loci_with_a_hit_on_a_SHOWN_element"],
            },
            "new_frame": {
                "intervals": iv.n,
                "indexed_against": "the CRISPRi element set",
                "elements_inside_it": len(att_iv),
                "elements_unassessed": 0,
                "elements_with_a_hit": len(with_hit),
                "attachments_element_by_record": records,
                "loci_with_a_hit": len(hit_groups),
                "addable_loci_on_a_shown_element": len(shown_hit_loci),
            },
            "not_a_comparison_of_biology": (
                "the two frames index different element sets, so the difference between them is "
                "coverage and not an effect. Nothing here says the old numbers were wrong: they "
                "were about the elements they could see"
            ),
        },
        "stream": {
            "tissues": distil_summary.get("tissues"),
            "pairs_scanned": distil_summary.get("pairs_scanned_this_run"),
            "hits_kept": distil_summary.get("hits_this_run"),
            "mb_streamed": distil_summary.get("mb_streamed_this_run"),
            "seconds": distil_summary.get("seconds"),
            "bytes_of_the_archive_written_to_disk": 0,
            "kept_where": KNOWLEDGE.as_posix(),
        },
        "cannot_establish": [
            "that an element assessed and carrying no retained hit has no eQTL. It has none in "
            "GTEx v8's significant single-tissue set, which is a call set over 49 tissues and a "
            "finite number of donors, and the CRISPRi cell lines are not among those tissues",
            "that an element carrying an eQTL is regulated by it, or that the eGene is its target. "
            "A cis-eQTL is an association over donors and linkage carries the signal some distance",
            "any agreement between a GTEx eQTL and a CRISPRi effect. Nothing here reads one against "
            "the other, and GTEx tissues are not the CRISPRi cell lines",
            "that the assessed figure is a finding. It is the frame's own definition",
            "the independent-locus grouping is operational and pools cell types, NOT established "
            "biological independence",
        ],
        "alphagenome_requests": 0,
        "money": "none: GTEx v8 is a public bucket and no paid call was made",
    }
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "GTEx v8 single-tissue cis-eQTLs (significant variant-gene pairs), hg38",
                "version": (
                    f"streamed from {eqtl.GTEX_EQTL_TAR} and re-distilled against the CRISPRi "
                    f"element set; the retained hits are under {KNOWLEDGE.as_posix()}"
                ),
            },
            {
                "accession": "EngreitzLab/CRISPR_comparison (Gschwind et al. 2025)",
                "version": "the two local benchmark tables under data/knowledge/crispri",
            },
            {
                "accession": "the response map's fourth-increment count",
                "version": "data/results/response_map_increment4_count.json, by sha256",
            },
        ],
        "inputs": _inputs()
        + [
            mf.input_entry(ROOT / f"data/results/{REGISTRATION}.json"),
            mf.input_entry(ROOT / "data/results/response_map_increment4_count.json"),
        ],
        "inputs_opened_beside_the_declared_results": {
            "why_they_are_listed": (
                "an input is what the writer opened, not what it said it read. The audit hook in "
                "`genomeos.manifest` records every open under `data/` and `save_result` refuses a "
                "result whose traced reads exceed its declared inputs"
            ),
            "groups": {k: len(v) for k, v in c3.OPENED_BESIDE_THE_RESULTS.items()},
            "the_pointer_case": (
                "`response_map2.candidates` opens the 24 per-chromosome element tables through a "
                "result file's `elements_where` POINTER field - about 1.4 GB that no reading of the "
                "source can name"
            ),
        },
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "eqtl_margin_bp": MARGIN,
            "reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
            "independent_locus_span": cell2.INDEPENDENT_LOCUS_SPAN,
            "floor_addable_loci_on_a_shown_element": FLOOR_ADDABLE_LOCI_ON_A_SHOWN_ELEMENT,
            "chromosomes": list(rm2.CHROMS),
        },
        "exclusions": [
            "no element of the frame is excluded, and an element assessed with no hit is reported "
            "as such rather than dropped",
            "GTEx's own significance set is the population; no further filter is applied here",
        ],
        "partitions": {
            "assessed": "what the frame could see at all, by construction",
            "found": "what carries a retained hit, which is the measurement",
        },
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }
    return payload


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--register", action="store_true", help="write the frame and rule, before any byte")
    ap.add_argument("--distil", action="store_true", help="stream the archive once against the frame")
    ap.add_argument("--count", action="store_true", help="what the new frame assessed, and found")
    ap.add_argument("--tissues", default="", help="comma list of GTEx tissues (default: all 49)")
    args = ap.parse_args()
    if not (args.register or args.distil or args.count):
        raise SystemExit("one of --register, --distil or --count is required")

    f = frame()
    iv = intervals(f["att_iv"])

    if args.register:
        payload = registration_payload(f, iv)
        path = save_result(REGISTRATION, payload)
        print(f"{REGISTRATION}: {path}")
        print(f"  frame: {payload['frame']['elements']} elements, {iv.n} intervals indexed")
        print(f"  floor: {FLOOR_ADDABLE_LOCI_ON_A_SHOWN_ELEMENT} (increment 3 built 6)")
        print(f"  hits go to {KNOWLEDGE.as_posix()}, nothing overwrites data/knowledge/gtex")
        return 0

    if args.distil:
        reg = ROOT / f"data/results/{REGISTRATION}.json"
        if not reg.exists():
            raise RefusedError(
                f"{reg} does not exist: the frame and the retained-hit rule are registered before "
                "any byte is streamed, never after"
            )
        registered = json.loads(reg.read_text())
        if registered["frame"]["intervals_indexed"] != iv.n:
            raise RefusedError(
                "the frame has moved since it was registered: "
                f"{registered['frame']['intervals_indexed']} intervals registered, {iv.n} now. "
                "Nothing is streamed against an unregistered frame"
            )
        if eqtl.KNOWLEDGE.resolve() == KNOWLEDGE.resolve():
            raise RefusedError(
                "refusing to distil into data/knowledge/gtex: it is a sha256-pinned input of four "
                "committed results"
            )
        KNOWLEDGE.mkdir(parents=True, exist_ok=True)
        tissues = [t for t in args.tissues.split(",") if t]

        def progress(msg: str) -> None:
            heartbeat("eqtl_crispri_frame")
            print(msg, flush=True)

        s = eqtl.distil(iv, knowledge=KNOWLEDGE, tissues=tissues or None, progress=progress)
        print(json.dumps(s, indent=1))
        return 0

    if not any(KNOWLEDGE.glob("hits_*.tsv")):
        raise RefusedError(f"no retained hits under {KNOWLEDGE}: run --distil first")
    payload = count_payload(f, iv)
    path = save_result(RESULT, payload)
    print(f"{RESULT}: {path}")
    a, g = payload["assessed"], payload["found"]
    print(f"  assessed {a['assessed']} of {a['elements_of_the_frame']} BY CONSTRUCTION, 0 unassessed")
    print(
        f"  found {g['elements_carrying_at_least_one_retained_hit']} elements with a hit, "
        f"{g['elements_assessed_and_carrying_none']} assessed with none, "
        f"{g['distinct_rows_attached_inside_an_element']} rows attached "
        f"({g['attachments_element_by_record']} attachments) over {g['tissues']} tissues"
    )
    lo = payload["loci"]
    print(
        f"  loci: {lo['loci_with_an_eqtl_inside_an_attached_element']} of "
        f"{lo['crispri_loci_total']}; addable on a shown element "
        f"{lo['addable_loci_with_a_hit_on_a_SHOWN_element']} "
        f"(any attached element {lo['addable_loci_with_a_hit_on_ANY_attached_element']})"
    )
    fl = payload["the_floor"]
    print(f"  floor {fl['floor']}: reached {fl['reached']}, cleared {fl['cleared']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
