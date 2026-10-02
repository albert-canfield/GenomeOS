# SPDX-License-Identifier: AGPL-3.0-or-later
"""GTEx v8 DAP-G read BY HTTP RANGE against the 1,505 CRISPRi element intervals, and its controls.

`data/results/finemap_coverage.json` found that the fine-mapped record already on disk answers for
almost none of the element set: DAP-G was read over 554 of the 1,505 elements, its posterior survives
on disk for 48, 14 carried a variant at PIP >= 0.5 and 0 carried one for the element's own linked gene.
That made the 1,262 of the GTEx re-distillation an LD-TAGGED count and left the localised question
unanswered rather than answered in the negative.

This writer asks it properly. The posterior lives in a remote bigBed, UCSC `gtexEqtlDapg`, and is read
HERE BY RANGE - request by request over the intervals, never by downloading the file - exactly as
`genomeos/attribution/human_panel.py` already reads it over the panel's units, through that module's own
`track_rows`, so one implementation fetches this track everywhere.

WHAT THE READ DOES AND DOES NOT CHANGE. It widens the frame a posterior could be seen over from 48
elements to the whole 1,505, because every element is read. That is a property of the frame and NOT a
result: the re-distillation's registration says so in its own words and this writer carries them. It
does not make any element carry anything.

EVERYTHING QUANTITATIVE IS IMPORTED. Not one number below is chosen here:

  * PIP 0.5 is `executor.DAPG_PIP`;
  * the margin is `scripts/eqtl_targets.MARGIN` through `response_map_increment3_count.EQTL_MARGIN`;
  * the frame is `scripts/finemap_coverage.elements`, the same 1,505 the coverage count used, which is
    `response_map2.candidates`'s own attachment;
  * the floors are `fresh.POSITIVE_FLOOR` (30) and `fresh.LOCUS_FLOOR` (20, which is
    `cell2.POOLED_LOCUS_FLOOR`), the pair lane-reptest applied, fixed before this lane existed;
  * the baseline's controls are `executor.CONTROLS_PER_UNIT` (2) per element, each within
    `executor.CONTROL_WINDOW` (250 kb) on the same chromosome, drawn under `executor.SEED`;
  * the density strata are deciles by `wiring.cuts` and `wiring.bin_of` with `wiring.ACTIVITY_BINS`
    bins, and the draw gets `wiring.DERANGE_TRIES` attempts;
  * the locus convention is `cell2.group`'s, with its wording that it is an operational grouping and
    NOT established biological independence.

THE CONTROLS ARE DRAWN AND REGISTERED BEFORE THE READ, and the read covers elements AND controls in one
pass. A baseline drawn after the hits are seen is not a baseline, and a second fetch for the controls
would let this lane choose when to stop.
"""

from __future__ import annotations

import argparse
import bisect
import collections
import hashlib
import random
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import finemap_coverage as fc  # noqa: E402
import response_map_increment3_count as c3  # noqa: E402

from genomeos import manifest as mf  # noqa: E402
from genomeos import response_map2 as rm2  # noqa: E402
from genomeos.attribution import (  # noqa: E402
    cell2,
    fresh,
)
from genomeos.attribution import executor as ex  # noqa: E402
from genomeos.attribution import human_panel as hp  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.attribution import wiring as wr  # noqa: E402
from genomeos.results import save_result  # noqa: E402

REGISTRATION = "finemap_dapg_registration"
RESULT = "finemap_dapg"
ENTRY = "scripts/finemap_dapg.py"
OWN_CODE = (ENTRY, "scripts/finemap_coverage.py")

#: Where the retained rows land. NOT data/knowledge/gtex and NOT data/knowledge/gtex_crispri, both of
#: which are sha256-pinned inputs of committed results.
KNOWLEDGE = Path("data/knowledge/gtex_dapg_crispri")

TRACK_KEY = "gtex_dapg"
PIP = ex.DAPG_PIP
MARGIN = c3.EQTL_MARGIN
CONTROLS_PER_ELEMENT = ex.CONTROLS_PER_UNIT
CONTROL_WINDOW = ex.CONTROL_WINDOW
DENSITY_BINS = wr.ACTIVITY_BINS
DRAW_TRIES = wr.DERANGE_TRIES
SEED = ex.SEED
LINK_FLOOR = fresh.POSITIVE_FLOOR
LOCUS_FLOOR = fresh.LOCUS_FLOOR

#: The precedent this read's cost is held against, from the 24 panel results already on disk. If the
#: read is NOT a fraction of this it stops and reports the figure instead of continuing.
PRECEDENT_REQUESTS = 222
PRECEDENT_MB = 170.4

LD_TAGGED_LABEL = "LD-tagged association, not localised"

#: Carried word for word from the re-distillation's registration, and a test compares it with the file.
ASSESSED_BY_CONSTRUCTION = (
    "every element of the frame is assessed BY CONSTRUCTION, because the frame is built from them. "
    "That is not a result and may not be reported as one: the figure that is a measurement is how "
    "many of them carry a retained hit, and the two are kept apart below for exactly that reason"
)


class RefusedError(ValueError):
    """Raised instead of reading or saving when a precondition does not hold."""


# --------------------------------------------------------------------------- the frame


def frame() -> dict[str, dict[str, Any]]:
    """The 1,505 elements, from `finemap_coverage` so the two results cannot index different sets."""
    return fc.elements()


def gene_bodies(chrom: str) -> list[tuple[int, int]]:
    """Protein-coding gene bodies on one chromosome, from the one GENCODE resolution already used."""
    from genomeos.genome import Annotation
    from genomeos.genome.annotation import default_gencode

    path = default_gencode({chrom})
    if path is None:
        return []
    ann = Annotation.from_gff3(path, {chrom})
    out = [
        (g.locus.start, g.locus.end)
        for g in ann.genes.values()
        if g.type == "protein_coding" and g.locus.chrom == chrom
    ]
    out.sort()
    return out


#: The longest protein-coding gene body on the human genome is under this, so a backward scan that
#: stops here cannot miss a gene that starts before the window and reaches into it. Checked by a test
#: against the annotation actually on disk rather than taken on trust.
MAX_GENE_BP = 3_000_000


def density(bodies: list[tuple[int, int]], starts: list[int], lo: int, hi: int) -> int:
    """Protein-coding gene bodies overlapping [lo, hi). The gene-density covariate, one definition.

    `bodies` is sorted by start, so every gene that can overlap the window starts in
    [lo - MAX_GENE_BP, hi). A gene longer than that would be missed, which is why a test holds the
    constant against the annotation on disk instead of trusting the comment above it.
    """
    n = 0
    i = bisect.bisect_left(starts, max(0, lo - MAX_GENE_BP))
    while i < len(bodies) and bodies[i][0] < hi:
        if bodies[i][1] > lo:
            n += 1
        i += 1
    return n


def controls(els: dict[str, dict[str, Any]]) -> tuple[dict[str, list[tuple[int, int]]], dict[str, Any]]:
    """Matched control windows per element: same length exactly, same gene-density decile, same
    chromosome, within CONTROL_WINDOW, overlapping no element of the frame.

    Drawn here, BEFORE the read, under the imported seed, and their digest goes in the registration.
    An element for which no matched control is found after DRAW_TRIES draws is reported as UNMATCHED
    and left out of the baseline; it is never silently dropped and never matched on fewer covariates.
    """
    by_chrom: dict[str, list[tuple[int, int, str]]] = collections.defaultdict(list)
    for eid, v in els.items():
        by_chrom[v["chrom"]].append((v["start"], v["end"], eid))
    taken: dict[str, list[tuple[int, int]]] = {}
    starts_taken: dict[str, list[int]] = {}
    for c, rows in by_chrom.items():
        iv = sorted((max(0, s - MARGIN), e + MARGIN) for s, e, _ in rows)
        merged: list[list[int]] = []
        for s, e in iv:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        taken[c] = [(a, b) for a, b in merged]
        starts_taken[c] = [a for a, _ in taken[c]]

    def clashes(c: str, lo: int, hi: int) -> bool:
        t, st = taken[c], starts_taken[c]
        i = bisect.bisect_right(st, hi) - 1
        while i >= 0 and t[i][1] > lo:
            if t[i][0] < hi and t[i][1] > lo:
                return True
            i -= 1
        return False

    # density deciles over the frame's OWN elements, so the strata are the population's
    bodies: dict[str, list[tuple[int, int]]] = {}
    starts_b: dict[str, list[int]] = {}
    for c in by_chrom:
        bodies[c] = gene_bodies(c)
        starts_b[c] = [b[0] for b in bodies[c]]
    el_density = {
        eid: density(bodies[v["chrom"]], starts_b[v["chrom"]], v["start"] - MARGIN, v["end"] + MARGIN)
        for eid, v in els.items()
    }
    edges = wr.cuts([float(x) for x in el_density.values()], DENSITY_BINS)

    out: dict[str, list[tuple[int, int]]] = {}
    unmatched: list[str] = []
    partial: list[str] = []
    for eid in sorted(els):
        v = els[eid]
        c, s, e = v["chrom"], v["start"], v["end"]
        length = e - s
        want = wr.bin_of(float(el_density[eid]), edges)
        rng = random.Random(f"{SEED}:{eid}")
        got: list[tuple[int, int]] = []
        for _ in range(DRAW_TRIES):
            if len(got) >= CONTROLS_PER_ELEMENT:
                break
            offset = rng.randint(-CONTROL_WINDOW, CONTROL_WINDOW)
            lo = s + offset
            if lo < 0:
                continue
            hi = lo + length
            if clashes(c, lo - MARGIN, hi + MARGIN):
                continue
            if any(abs(lo - a) < length + 2 * MARGIN for a, _ in got):
                continue
            d = density(bodies[c], starts_b[c], lo - MARGIN, hi + MARGIN)
            if wr.bin_of(float(d), edges) != want:
                continue
            got.append((lo, hi))
        if not got:
            unmatched.append(eid)
        elif len(got) < CONTROLS_PER_ELEMENT:
            partial.append(eid)
        if got:
            out[eid] = got
    info = {
        "controls_per_element": CONTROLS_PER_ELEMENT,
        "elements_with_the_full_number": sum(1 for v in out.values() if len(v) == CONTROLS_PER_ELEMENT),
        "elements_with_fewer": len(partial),
        "elements_with_none_unmatched": len(unmatched),
        "unmatched_elements": unmatched[:50],
        "control_windows_drawn": sum(len(v) for v in out.values()),
        "matched_on": [
            "length, exactly equal to the element's",
            f"gene-density decile, protein-coding bodies overlapping the window widened by {MARGIN}",
            "the same chromosome",
            f"within {CONTROL_WINDOW} bases of the element",
            "overlapping no element of the frame widened by the margin, and not each other",
        ],
        "density_bins": DENSITY_BINS,
        "density_bin_edges": [round(x, 3) for x in edges],
        "seed": SEED,
        "draw_tries_per_control": DRAW_TRIES,
        "digest": _digest(out),
        "an_unmatched_element_is_not_a_zero": (
            "an element with no matched control contributes to neither side of the baseline and is "
            "reported by count. Matching it on fewer covariates would make the excess a comparison of "
            "unlike windows, which is the thing the baseline exists to prevent"
        ),
    }
    return out, info


def _digest(ctrl: dict[str, list[tuple[int, int]]]) -> str:
    h = hashlib.sha256()
    for eid in sorted(ctrl):
        h.update(eid.encode())
        for lo, hi in ctrl[eid]:
            h.update(f":{lo}-{hi}".encode())
    return h.hexdigest()


def intervals_to_read(
    els: dict[str, dict[str, Any]], ctrl: dict[str, list[tuple[int, int]]]
) -> dict[str, list[tuple[int, int]]]:
    """Every interval the read covers, elements and controls together, per chromosome."""
    out: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    for eid, v in els.items():
        out[v["chrom"]].append((max(0, v["start"] - MARGIN), v["end"] + MARGIN))
        for lo, hi in ctrl.get(eid, []):
            out[v["chrom"]].append((max(0, lo - MARGIN), hi + MARGIN))
    return {c: sorted(set(v)) for c, v in out.items()}


def registration_payload(
    els: dict[str, dict[str, Any]], ctrl_info: dict[str, Any], iv: dict[str, list[tuple[int, int]]]
) -> dict[str, Any]:
    return {
        "result": REGISTRATION,
        "date": "2026-10-02",
        "registered_before": "any byte of the DAP-G track was read",
        "question": (
            "how many of the 1,505 compiled CRISPRi elements carry a FINE-MAPPED cis-eQTL (DAP-G PIP "
            "at or above the imported threshold) for the element's OWN linked gene, and how much of "
            "that is in excess of length- and gene-density-matched windows"
        ),
        "why": (
            "data/results/finemap_coverage.json established that the fine-mapped record on disk "
            "covered 48 of the 1,505 for the threshold, so the localised question was UNANSWERED "
            "rather than answered. The 1,262 of the re-distillation counts GTEx's significant "
            f"single-tissue set and is a {LD_TAGGED_LABEL}"
        ),
        "the_read": {
            "track": hp.BIGBEDS[TRACK_KEY],
            "source_note": hp.SOURCE_NOTES[TRACK_KEY],
            "how": (
                "BY HTTP RANGE through `genomeos.attribution.human_panel.track_rows`, the one "
                "implementation that already reads this track over the panel's units. The file is "
                "NEVER downloaded"
            ),
            "chromosomes": sorted(iv),
            "intervals": sum(len(v) for v in iv.values()),
            "cost_precedent": {
                "requests": PRECEDENT_REQUESTS,
                "mb_fetched": PRECEDENT_MB,
                "from": (
                    "the cost.storage_tracks.gtex_dapg fields of the 24 human_panel_chr*.json "
                    "results, which covered the panel's millions of unit intervals"
                ),
                "stop_rule": (
                    "this read covers a few thousand intervals and must therefore be a FRACTION of "
                    "the precedent. If it is not, the run stops and reports the figure rather than "
                    "continuing"
                ),
            },
            "retained_hit_rule": (
                f"a DAP-G row whose variant base (`human_panel._eqtl_pos`) falls inside an indexed "
                f"interval and whose PIP is at or above {PIP}. No tissue, gene class or effect size is "
                "filtered and this lane sets no cut-off of its own"
            ),
            "where_the_rows_go": KNOWLEDGE.as_posix(),
            "nothing_is_overwritten": (
                "not data/knowledge/gtex and not data/knowledge/gtex_crispri, both sha256-pinned "
                "inputs of committed results"
            ),
        },
        "the_two_figures_kept_apart": {
            "assessed": "how many of the frame's elements could have been seen at all",
            "carrying_a_hit": "how many of them carry at least one retained hit",
            "why": (
                "conflating them is how an unassessed element becomes an element with no eQTL. The "
                "first is a property of the frame and the second is a measurement"
            ),
            "quoted_from": "data/results/eqtl_crispri_frame_registration.json",
            "what_the_read_changes": (
                "the read makes all 1,505 assessed, and "
                f"{ASSESSED_BY_CONSTRUCTION}. So 1,505 is the honest denominator for 'is a "
                "fine-mapped variant here', and it is NOT a denominator for 'is this element the "
                "cause': the first is answered by the read and the second is not answered at all"
            ),
            "the_limit_the_read_cannot_remove": (
                "DAP-G fine-maps cis-eQTLs per gene, so an element whose linked gene GTEx never "
                "tested cannot carry a hit for it whatever is there. This read sees only the genes "
                "with a row in its own windows and therefore CANNOT count those elements. The limit "
                "is named, not counted, and no ratio is adjusted for it"
            ),
        },
        "conditions": {
            "a_gene_match_is_mandatory": (
                "the retained row's geneName must be the element's own linked gene, resolved through "
                "`scripts/eqtl_targets.py symbol_map`, the one Ensembl-to-symbol resolution. An eQTL "
                "for some other gene inside the window is not evidence about this chain"
            ),
            "b_the_reportable_quantity_is_the_excess": (
                "the raw gene-matched count is NOT the result. The result is its excess over the "
                "matched control windows, because with a dense track a window of any kind catches "
                "something and a count cannot tell signal from density"
            ),
            "floors": {
                "gene_matched_elements_at_least": LINK_FLOOR,
                "imported_from": "genomeos.attribution.fresh.POSITIVE_FLOOR",
                "independent_loci_at_least": LOCUS_FLOOR,
                "loci_imported_from": (
                    "genomeos.attribution.fresh.LOCUS_FLOOR, which is "
                    "genomeos.attribution.cell2.POOLED_LOCUS_FLOOR"
                ),
                "applied_to": (
                    "the GENE-MATCHED carrying count and its independent loci, never the raw "
                    "carrying count and never the assessed count"
                ),
                "neither_is_chosen_here": (
                    "both were fixed before this lane existed and lane-reptest applied the same pair. "
                    "Neither moves after the count is seen, not down to admit a population and not "
                    "sidestepped by dropping the gene match or by widening the margin"
                ),
                "if_the_floors_are_missed": (
                    "no baseline is computed and that is the result. A matched-window excess over a "
                    "handful of elements is a number nobody can read"
                ),
            },
        },
        "the_baseline": ctrl_info,
        "frame": {
            "elements": len(els),
            "intervals_indexed": sum(len(v) for v in intervals_to_read(els, {}).values()),
            "from": (
                "scripts/finemap_coverage.elements, which is `response_map2.candidates`'s own "
                f"attachment at `measured.RECIPROCAL_OVERLAP = {ms.RECIPROCAL_OVERLAP}`"
            ),
            "margin_bp": MARGIN,
            "margin_source": (
                "scripts/eqtl_targets.MARGIN through "
                "scripts/response_map_increment3_count.EQTL_MARGIN, not chosen here"
            ),
        },
        "locus_convention": (
            "genomeos.attribution.cell2.group's, carrying its wording that it is an operational "
            "grouping and NOT established biological independence"
        ),
        "alphagenome_requests": 0,
        "money": (
            "none. UCSC's gbdb is a public host and a range read of a public bigBed is not a paid "
            "call. If anything in this lane would cost money the run stops: only Albert can "
            "authorise that"
        ),
        "result_manifest": {
            "sources": [
                {"accession": hp.SOURCE_NOTES[TRACK_KEY], "version": hp.BIGBEDS[TRACK_KEY]},
                {
                    "accession": "EngreitzLab/CRISPR_comparison (Gschwind et al. 2025)",
                    "version": "the two local benchmark tables under data/knowledge/crispri",
                },
            ],
            "inputs": fc._inputs(),
            "inputs_opened_beside_the_declared_results": {
                "why_they_are_listed": (
                    "an input is what the writer opened, not what it said it read. The audit hook in "
                    "`genomeos.manifest` records every open under `data/` and `save_result` refuses "
                    "a result whose traced reads exceed its declared inputs"
                ),
                "groups": {k: len(v) for k, v in c3.OPENED_BESIDE_THE_RESULTS.items()},
                "the_pointer_case": (
                    "`response_map2.candidates` opens the 24 per-chromosome element tables through a "
                    "result file's `elements_where` POINTER field - about 1.4 GB that no reading of "
                    "the source can name, and how 17 published results read a table none declared"
                ),
            },
            "assembly": "GRCh38",
            "coordinates": {"base": 0, "interval": "half-open"},
            "parameters": {
                "dapg_pip": PIP,
                "eqtl_margin_bp": MARGIN,
                "reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
                "independent_locus_span": cell2.INDEPENDENT_LOCUS_SPAN,
                "controls_per_element": CONTROLS_PER_ELEMENT,
                "control_window_bp": CONTROL_WINDOW,
                "density_bins": DENSITY_BINS,
                "draw_tries": DRAW_TRIES,
                "seed": SEED,
                "gene_matched_floor": LINK_FLOOR,
                "independent_loci_floor": LOCUS_FLOOR,
                "chromosomes": list(rm2.CHROMS),
            },
            "exclusions": [
                "no element of the frame is excluded: it is every element a valid CRISPRi pair attaches to",
                "no DAP-G row is excluded other than by the imported PIP threshold",
                "an element with no matched control is excluded from the BASELINE only, by count",
            ],
            "partitions": {
                "frame": "the CRISPRi element set the track is read against",
                "controls": "length- and gene-density-matched windows drawn before the read",
            },
            "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--register", action="store_true", help="write the frame, controls and rule first")
    args = ap.parse_args()
    if not args.register:
        raise SystemExit("--register is required: nothing is read before the registration exists")
    els = frame()
    ctrl, info = controls(els)
    iv = intervals_to_read(els, ctrl)
    payload = registration_payload(els, info, iv)
    path = save_result(REGISTRATION, payload)
    print(f"{REGISTRATION}: {path}")
    print(f"  frame: {len(els)} elements; controls drawn: {info['control_windows_drawn']}")
    print(
        f"  full number {info['elements_with_the_full_number']}, "
        f"fewer {info['elements_with_fewer']}, none {info['elements_with_none_unmatched']}"
    )
    print(f"  intervals to read: {sum(len(v) for v in iv.values())} over {len(iv)} chromosomes")
    print(f"  floors: {LINK_FLOOR} gene-matched elements, {LOCUS_FLOOR} independent loci (both imported)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
