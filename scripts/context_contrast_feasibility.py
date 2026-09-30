# SPDX-License-Identifier: AGPL-3.0-or-later
"""Audit B, cached context contrast, checkpoint 1: feasibility first (lane-contrast, 2026-09-29). No score.

The audit asks whether the per-element cached AlphaGenome deletion answer supports one added feature: the
same-cell signed deletion minus the median other-cell signed deletion, per element, with a frozen
aggregation. This script answers only whether it can, and stops:

1. semantics: what the cache holds per element and gene (signed or absolute, which tracks, whether values
   below a recording threshold are stored or missing), from the code that wrote it and from the answers
   the admissible pairs match;
2. completeness on the admissible endpoints (lane-prior's admissibility, a48e9e8): Gasperini2019 on both of
   the pilot's partitions and Schraivogel2020 on the seen partition, per k other cells, with contact
   coverage, abstentions and independent loci;
3. duplication: whether the contrast is nearly a function of the compiled target, or of the same-cell
   deletion that the baseline already carries.

GO_RULE, written in this file before any count below was made, decides the reading. No labelling is
scored. The only outcome read is the harness's floor check on the complete pairs (positives and
negatives, as `prior_only_test.partition_counts` counts them). No model request is made, the 4DN contact
answers are read from their local cache only (a miss is counted, never fetched), and nothing is downloaded.

    uv run python scripts/context_contrast_feasibility.py

writes data/results/context_contrast_feasibility.json.
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import hashlib
import importlib.util
import json
import resource
import statistics
import sys
import time
from collections import Counter, defaultdict
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import holdout as ho  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.genome import hic_contact as hc  # noqa: E402
from genomeos.predict import enhancer_target as et  # noqa: E402
from genomeos.results import save_result  # noqa: E402

_SPEC = importlib.util.spec_from_file_location("prior_only_test", ROOT / "scripts" / "prior_only_test.py")
pot = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(pot)  # lane-prior's admissibility and partitions, read and never changed

NAME = "context_contrast_feasibility"
CHECKPOINT = "2026-09-29"
STATUS = "feasibility checkpoint: no labelling scored"
MAY_BE_CALLED = "a feasibility count on development evidence"
MAY_NOT_BE_CALLED = ("a score", "a test of the feature", "a validation", *ho.MAY_NOT_BE_CALLED)
ENDPOINT = pot.ENDPOINT

# ==================================================================================================
# Written before any count of this script was made (2026-09-29, lane-contrast).
# ==================================================================================================
SAME_CELL = "K562"  # both admissible studies screened K562; asserted on the units
OTHER_CELLS = tuple(c for c in et.CELLS if c != SAME_CELL)  # HepG2, GM12878, IMR-90: all the cache keeps
KS = (1, 2, 3)
#: the sweep's track count, from track_fingerprint over every answer (genomeos/manifest.py ALPHAGENOME_SWEEP)
FULL_TRACKS = mf.ALPHAGENOME_SWEEP["track_metadata"]["observed"]["tracks_per_element_modal"]
LOW = 0.05  # the adapter's default recording threshold, which docs/ATTRIBUTION.md says the cache used
ROUNDING = 1e-4  # by_cell values are rounded to four places, so a difference of two is known to +/- 1e-4
HALF_WINDOW = 524_288  # the scorer's reach each way from the deleted element (1 Mb window)
CONTACT_CELL = SAME_CELL
#: (source, partition, chromosomes, role): lane-prior's admissible endpoints, as prior_only_test fixes them
ENDPOINT_SETS = (
    ("crispri:Gasperini2019", "fresh", pot.FRESH, "primary"),
    ("crispri:Gasperini2019", "seen", pot.SEEN, "beside"),
    ("crispri:Schraivogel2020", "seen", pot.SEEN, "beside"),
)
PRIMARY = ENDPOINT_SETS[0]
USEFUL_COVERAGE = pot.USEFUL_COVERAGE  # 0.80, lane-prior's useful coverage
DUPLICATE_R2 = 0.80
DUPLICATE_SPEARMAN = 0.95
UNCENSORED_SHARE = 0.99
RESOLUTION_SHARE = 0.50

SEMANTICS_FROM_CODE = {
    "writer": "scripts/enhancer_targets_all.py worker_scorer builds the scorer with "
    "alphagenome_adapter._live_scorer(threshold=0.0); the scorer emits (gene, track name, value) for "
    "every gene-track value with abs(value) > threshold, the track name being the track's GTEx tissue or, "
    "when it has none, its biosample name",
    "aggregation": "enhancer_target.aggregate keeps per gene: n_tracks (the values emitted), their mean, "
    "the largest drop and the largest rise with the track names they sat on, and by_cell[name] for each "
    "of K562, HepG2, GM12878 and IMR-90: the value of the LAST emitted track carrying that name, rounded "
    "to four places, signed. Not a replicate aggregate, not the largest, and the track's identity is not "
    "stored. Genes are keyed by the model's gene name, so two genes of one name share one row",
    "what_is_not_on_disk": "the track table (how many tracks carry each cell name and in which order), "
    "every other track of a cell name, the per-track values, and the gene's Ensembl id",
    "docs_statement_checked": "docs/ATTRIBUTION.md ('Found but not fixed: by_cell keeps one value per cell "
    "name') says the kept value is the last track above the scorer's 0.05 recording threshold. That is the "
    "adapter's default; the sweep passed 0.0. Which one the matched answers were written at is measured "
    "below (an answer that stores any by_cell value with |v| <= 0.05 was not written at 0.05)",
}
RELIABLE_RULE = (
    "a pair's same-cell value is reliable when the gene's row (i) comes from an uncensored answer, one that "
    "stores at least one by_cell value with |v| <= 0.05 and so was not written at the 0.05 threshold; "
    f"(ii) has n_tracks equal to the sweep's full track count ({FULL_TRACKS}), so no track of the gene "
    "returned exactly 0.0 and was dropped and no two genes of one name share the row, which makes each "
    "by_cell value the last track of that name in the model's order, one fixed track per cell; and (iii) "
    "holds a K562 value. Rows short of the full count may hold a cell's earlier track in place of its last "
    "one, or no value; they are counted, never repaired"
)
CONTRAST_RULE = (
    "per element and gene: contrast = by_cell['K562'] - median(by_cell[c] for c in HepG2, GM12878, IMR-90 "
    "that the row holds), defined only when K562 is held and at least k other cells are; the median of two "
    "is their mean. A missing cell is missing: it is never read as zero, and no consensus of the others "
    "stands in for it. The other cells are a context, not negative controls: the feature claims no cell's "
    "value is null"
)
PAIR_RULE = (
    "a pair's element is the cached answer sharing at least one base with the pair's element whose row for "
    "the pair's gene (matched by the benchmark's symbol, else its Ensembl id) is reliable; among several, "
    "the one sharing the most bases with the pair's element, then the lower start, then the lower id. "
    "Loose per-element files take precedence over the chromosome archive, as enhancer_target.load_cached "
    "reads them. A pair with no such answer abstains, with the reason of its best-overlapping matched row"
)
CONTACT_RULE = (
    f"K562 in-situ Hi-C ({hc.MATRICES[CONTACT_CELL]}, {hc.RESOLUTION:,} bp, the file's own balancing), the "
    "contact between the pair's midpoint bin and its TSS bin (the benchmark's startTSS), read through "
    "hic_contact.ContactSource from data/knowledge/hic_contact only; no matrix is opened, a miss is counted"
)
ACTIVITY_RULE = (
    "the benchmark's own DHS.RPM and H3K27ac.RPM columns of the pair's row (K562), present when both are "
    "non-empty and not NA; the input of crispri.py's 'activity + distance (+ deletion)' models"
)
DUPLICATION_RULE = (
    "the contrast is nearly a function of an existing feature, and the checkpoint a no-go, if on the "
    "reliable complete pairs (k = 3) of the three admissible endpoint sets pooled: (D1) a least-squares fit "
    "of the contrast on the compiled target's description of the pair (not named; or named, by action x "
    "whether the link's model cell is K562, each with an intercept and a slope on strength) explains at "
    f"least {DUPLICATE_R2:.2f} of its variance; or (D2) its Spearman correlation with the same-cell value, "
    f"which the deletion baseline already carries, is at least {DUPLICATE_SPEARMAN:.2f} in absolute value"
)
GO_RULE = (
    "a clear go needs all three: (G1, semantics) at least "
    f"{UNCENSORED_SHARE:.2f} of the answers the pairs match are uncensored, and fewer than "
    f"{RESOLUTION_SHARE:.2f} of the reliable complete pairs carry a contrast inside the rounding bound "
    f"(|c| <= {ROUNDING}); (G2, completeness) on the primary endpoint (Gasperini2019, fresh partition), the "
    "reliable complete pairs at k = 3 with activity and distance present clear the harness's floors (20 "
    f"positives, 20 negatives, 10 loci) and cover at least {USEFUL_COVERAGE:.2f} of the endpoint's units; "
    "(G3, duplication) neither D1 nor D2 fires. Contact coverage decides only whether comparison (b) would "
    "be on the same pairs or beside, never the reading. Anything short of a clear go ends the lane here"
)
READINGS = {
    "go": "the cache supports the feature on the admissible endpoints: register before any score",
    "no_go_semantics": "the cache cannot give a faithful signed same-cell value and other-cell median",
    "no_go_duplicate": "the contrast is nearly a function of a feature the comparison already holds",
    "no_go_completeness": "too few admissible pairs carry a reliable contrast for a test without coverage "
    "loss",
}
#: what earlier lanes read on these endpoints and this cache, by commit; every one is development evidence
SEEN_BEFORE = (
    {
        "what": "activity + distance (+ deletion, from 2026-09-22 read from this per-element cache) on the "
        "pooled K562 training file, which holds Gasperini2019 and Schraivogel2020; the same with measured "
        "K562 contact",
        "result": "crispri_benchmark, crispri_contact",
        "commit": "138824f, bda8a83, 5a31c39",
    },
    {
        "what": "hold-one-chromosome-out on all 10,356 training pairs with the same-cell deletion drop",
        "result": "crispri_published",
        "commit": "42d7b1b, a39073d",
    },
    {
        "what": "a calibrated probability fitted on the covered K562 training pairs from the K562 drop",
        "result": "target_calibration",
        "commit": "bdc2364",
    },
    {
        "what": "the signed K562 by_cell value against the measured sign (held-out arm)",
        "result": "crispri_direction, crispri_direction_both",
        "commit": "5c842ca, d492834",
    },
    {
        "what": "distance, unchanged labels and rest on every CRISPRi endpoint, and the AlphaGenome ablation",
        "result": "c4_holdout_scores, c4_alphagenome_ablation",
        "commit": "b7e4bf0",
    },
    {
        "what": "the pilot's labellings on the seven validation chromosomes (the prior-only lead)",
        "result": "pilot_biological_gate",
        "commit": "197c560",
    },
    {
        "what": "the prior-only test and its decomposition on both partitions of Gasperini2019 and the seen "
        "partition of Schraivogel2020 (compiled target +0.1825 [+0.1298, +0.2280]; H3K27ac over distance "
        "+0.0098 [-0.0188, +0.0349])",
        "result": "prior_only_test",
        "commit": "a48e9e8",
    },
    {
        "what": "the repression trace, which read by_cell and found one value kept per cell name",
        "result": "repression_trace",
        "commit": "641909e",
    },
    {
        "what": "S4 judged on these screens under rules v1, v2 and v3",
        "result": "attribution_correctness, _v2, _v3",
        "commit": "ab62d99, 4567219, d1e09ae",
    },
    {
        "what": "this lane, before writing the rules above: chr21's archive and the 3,209 loose answer files "
        "read for their recording threshold and track counts, no pair joined, nothing scored",
        "result": "none (a probe)",
        "commit": "n/a",
    },
)
BUDGET = "one process, cached inputs, at most one cache reader at a time; no model request, no download"

# ==================================================================================================
# Added 2026-09-29 after the provisional run (code 7870427, stamped 0183ef5), at the external reviewer's
# reading relayed by the coordinator. Wording and locus accounting only: the gate, the rules above and
# the reading are unchanged, and no other design is launched.
# ==================================================================================================
LOCI_NOTE = (
    "loci are frozen on the parent universe (every unit of the endpoint set, or of all three sets) before "
    "any filtering, and a filtered set's independent-locus count is the number of parent loci it falls "
    "in. Recomputing connected components after filtering splits them (on the primary endpoint 329 parent "
    "loci became 351), which is not new independent evidence; that recount is kept as descriptive only. "
    "clears_floors is left as registered"
)
READING_NOTES = {
    "scope": "a no-go of this registered complete-case design on coverage. It does not say that a "
    "context contrast cannot help",
    "reliable_rule_reading": f"n_tracks = {FULL_TRACKS} is a consistency check: every value of the gene "
    "was emitted, so an exact-zero track cannot have moved a cell's by_cell value to an earlier track of "
    "that name. It is not independent proof of which track by_cell holds, nor that one model version "
    "answered every element (no answer read records a requested version)",
    "possible_future_design": "a baseline with the feature optional (the baseline's own score wherever the "
    "contrast is missing) could keep the baseline's coverage. A possibility only: not registered, not "
    "launched, and nothing here tests it",
}


# --- the cache, streamed ---------------------------------------------------------------------------
_WS = " \t\n\r,"


def _skip(buf: str, i: int) -> int:
    n = len(buf)
    while i < n and buf[i] in _WS:
        i += 1
    return i


def stream_answers(path: Path, chunk: int = 1 << 23) -> Iterator[tuple[str, dict[str, Any]]]:
    """(element id, answer) for every entry of one gzipped archive (a JSON object written by `pack`), read
    a chunk at a time, so an archive is never held whole: the same entries `json.load` gives, in order."""
    dec = json.JSONDecoder()
    with gzip.open(path, "rt") as fh:
        buf = fh.read(chunk)
        eof = not buf
        i = buf.index("{") + 1
        while True:
            j = _skip(buf, i)
            if j < len(buf) and buf[j] == "}":
                return
            try:
                if j >= len(buf):
                    raise IndexError
                key, k = dec.raw_decode(buf, j)
                k = _skip(buf, k)
                if k >= len(buf):
                    raise IndexError
                if buf[k] != ":":
                    raise ValueError(f"{path}: expected ':' at {k}")
                k = _skip(buf, k + 1)
                if k >= len(buf):
                    raise IndexError
                val, end = dec.raw_decode(buf, k)
            except (json.JSONDecodeError, IndexError):
                if eof:
                    raise ValueError(f"{path}: truncated or malformed archive") from None
                if len(buf) - j > 64 * chunk:
                    raise ValueError(f"{path}: an entry larger than {64 * chunk} characters") from None
                more = fh.read(chunk)
                eof = not more
                buf, i = buf[j:] + more, 0
                continue
            yield key, val
            i = end


def cached_answers(
    chrom: str, keep: Callable[[dict[str, Any]], bool], cache: Path = et.CACHE
) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    """The answers of one chromosome that `keep(answer)` accepts: the archive streamed, then the loose
    per-element files, which replace an archive entry of the same id (enhancer_target.load_cached)."""
    kept: dict[str, dict[str, Any]] = {}
    stats: dict[str, Any] = {
        "archive_answers": 0,
        "loose_answers": 0,
        "kept_from_archive": 0,
        "kept_loose": 0,
    }
    loose_paths: list[str] = []
    arch = et.archive_path(chrom, cache)
    if arch.exists():
        for eid, hit in stream_answers(arch):
            stats["archive_answers"] += 1
            if keep(hit):
                kept[eid] = {**hit, "_from": "archive"}
                stats["kept_from_archive"] += 1
    d = cache / chrom
    for f in sorted(d.glob("*.json")) if d.exists() else []:
        stats["loose_answers"] += 1
        loose_paths.append(str(f))
        hit = json.loads(f.read_text())
        if keep(hit):
            kept[f.stem] = {**hit, "_from": "loose"}
            stats["kept_loose"] += 1
    stats["loose_paths"] = loose_paths
    return kept, stats


def overlap_filter(spans: list[tuple[int, int]]) -> Callable[[dict[str, Any]], bool]:
    """keep(answer): whether the answer's element shares at least one base with any of `spans`."""
    spans = sorted(spans)
    starts = [s for s, _ in spans]
    reach = max((e - s for s, e in spans), default=0)

    def keep(hit: dict[str, Any]) -> bool:
        s0, s1 = int(hit["start"]), int(hit["end"])
        lo = bisect.bisect_left(starts, s0 - reach)
        hi = bisect.bisect_left(starts, s1)
        return any(e > s0 and s < s1 for s, e in spans[lo:hi])

    return keep


# --- per answer and row ----------------------------------------------------------------------------
def by_cell_values(hit: dict[str, Any]) -> list[float]:
    return [float(v) for g in hit.get("genes") or [] for v in (g.get("by_cell") or {}).values()]


def uncensored(hit: dict[str, Any]) -> bool:
    """Whether an answer was written below the 0.05 threshold: it stores a by_cell value with |v| <= 0.05."""
    return any(abs(v) <= LOW for v in by_cell_values(hit))


def row_for(hit: dict[str, Any], gene: str, gene_id: str) -> tuple[dict[str, Any] | None, str]:
    """The answer's row for a gene: by symbol, else by Ensembl id (the model names a gene with no symbol
    by its id); (None, '') when the answer does not list the gene."""
    rows = hit.get("genes") or []
    for g in rows:
        if gene and g.get("gene") == gene:
            return g, "symbol"
    for g in rows:
        if gene_id and g.get("gene") == gene_id:
            return g, "ensembl_id"
    return None, ""


def row_status(hit: dict[str, Any], row: dict[str, Any]) -> str:
    """RELIABLE_RULE for one matched row: 'reliable', or the first condition it fails."""
    if not uncensored(hit):
        return "answer_censored"
    n = int(row.get("n_tracks") or 0)
    if n > FULL_TRACKS:
        return "row_merges_two_genes"
    if n < FULL_TRACKS:
        return "row_dropped_exact_zero_tracks"
    if SAME_CELL not in (row.get("by_cell") or {}):
        return "same_cell_missing"
    return "reliable"


def others_of(row: dict[str, Any]) -> list[float]:
    bc = row.get("by_cell") or {}
    return [float(bc[c]) for c in OTHER_CELLS if c in bc]


def contrast(row: dict[str, Any], k: int) -> float | None:
    """CONTRAST_RULE for one row: None when K562 is missing or fewer than k (at least one) other cells
    are held."""
    bc = row.get("by_cell") or {}
    if SAME_CELL not in bc:
        return None
    others = others_of(row)
    if len(others) < max(k, 1):
        return None
    return round(float(bc[SAME_CELL]) - statistics.median(others), 6)


def shared(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


# --- one pair --------------------------------------------------------------------------------------
def pair_record(unit: ho.Unit, hits: list[dict[str, Any]]) -> dict[str, Any]:
    """PAIR_RULE for one pair over the answers overlapping it: the chosen reliable row's values, the most
    other cells any matched row with a K562 value holds (`any_k`, -1 for none), or the abstention's reason."""
    rec: dict[str, Any] = {"answers_overlapping": len(hits), "any_k": -1}
    if not hits:
        rec["reason"] = "no_cached_answer_overlaps"
        return rec
    matched = []
    for h in hits:
        row, how = row_for(h, unit.gene, unit.gene_id)
        if row is not None:
            ov = shared(unit.start, unit.end, int(h["start"]), int(h["end"]))
            matched.append((ov, int(h["start"]), str(h["id"]), h, row, how))
    if not matched:
        if unit.tss is None:
            rec["reason"] = "gene_not_listed_no_tss"
        elif all(abs((int(h["start"]) + int(h["end"])) // 2 - unit.tss) > HALF_WINDOW for h in hits):
            rec["reason"] = "gene_not_listed_out_of_reach"
        else:
            rec["reason"] = "gene_not_listed_in_reach"
        return rec
    matched.sort(key=lambda m: (-m[0], m[1], m[2]))
    rec["matched_by"] = matched[0][5]
    rec["answers_matched"] = len(matched)
    rec["matched_uncensored"] = sum(uncensored(m[3]) for m in matched)
    rec["any_k"] = max(
        (len(others_of(m[4])) for m in matched if SAME_CELL in (m[4].get("by_cell") or {})), default=-1
    )
    reliable = [m for m in matched if row_status(m[3], m[4]) == "reliable"]
    if not reliable:
        rec["reason"] = row_status(matched[0][3], matched[0][4])
        return rec
    _, _, eid, hit, row, how = reliable[0]
    others = others_of(row)
    rec.update(
        {
            "reason": "reliable",
            "element": eid,
            "from": hit.get("_from", ""),
            "matched_by_chosen": how,
            "same": float(row["by_cell"][SAME_CELL]),
            "others": len(others),
            "other_median": statistics.median(others) if others else None,
            "contrast": {k: contrast(row, k) for k in KS},
        }
    )
    return rec


# --- the other inputs ------------------------------------------------------------------------------
def activity_index(rows: dict[str, list[dict[str, str]]]) -> dict[tuple, bool]:
    """ACTIVITY_RULE per (source, chrom, start, end, gene)."""
    out: dict[tuple, bool] = {}
    for src, rs in rows.items():
        for r in rs:
            ok = all(r.get(c) not in (None, "", "NA") for c in ("DHS.RPM", "H3K27ac.RPM"))
            key = (src, r["chrom"], int(r["chromStart"]), int(r["chromEnd"]), r["measuredGeneSymbol"])
            out[key] = out.get(key, False) or ok
    return out


def compiled_description(idx: Any, unit: ho.Unit) -> dict[str, Any]:
    """The compiled target's description of a pair: whether a predicted link sharing a base with the
    pair's element names its gene, and, for the strongest such link, its action, strength and cell."""
    named = [x for x in idx.near(unit.chrom, unit.start, unit.end) if x.gene == unit.gene]
    if not named:
        return {"named": False}
    x = max(named, key=lambda x: (x.strength, x.action, x.cell))
    return {"named": True, "action": x.action, "strength": x.strength, "cell_is_same": x.cell == SAME_CELL}


def floors(units: list[ho.Unit]) -> dict[str, Any]:
    """The harness's floor check, counted as prior_only_test.partition_counts counts it."""
    pos = int(sum(ho.target(u, ENDPOINT) or 0 for u in units))
    n_loci = max(ho.loci(units)) + 1 if units else 0
    return {
        "units": len(units),
        "positives": pos,
        "negatives": len(units) - pos,
        "loci": n_loci,
        "clears_floors": pos >= ho.MIN_POSITIVES
        and len(units) - pos >= ho.MIN_NEGATIVES
        and n_loci >= ho.MIN_LOCI,
    }


def parent_loci(fl: dict[str, Any], units: list[ho.Unit], parent: dict[ho.Unit, int]) -> dict[str, Any]:
    """A floor check with LOCI_NOTE applied: the independent-locus count is the number of parent-universe
    loci the units fall in; the count `floors` recomputes on the filtered units is kept as descriptive.
    `clears_floors` is left exactly as registered."""
    out = {k: v for k, v in fl.items() if k != "loci"}
    out["independent_loci_parent_universe"] = len({parent[u] for u in units})
    out["loci_recomputed_after_filtering_descriptive"] = fl["loci"]
    return out


def with_parent_loci(
    summary: dict[str, Any], us: list[ho.Unit], recs: list[dict[str, Any]]
) -> dict[str, Any]:
    """LOCI_NOTE and the loss steps applied to one endpoint set's summary after `set_summary` built it:
    the floor dicts count independent loci on the set's own units, and `where_coverage_is_lost` lists
    the pairs left after each step. Nothing the gate reads changes."""
    parent = dict(zip(us, ho.loci(us), strict=True))
    full_recs = [
        r
        for r in recs
        if r["reason"] == "reliable" and r["contrast"][3] is not None and r["activity"] and r["distance"]
    ]
    full = [r["unit"] for r in full_recs]
    full_c = [r["unit"] for r in full_recs if r["contact"]]
    key = "reliable_k3_with_activity_and_distance"
    fold = summary[key]["fold_sha256"]
    summary[key] = {**parent_loci(summary[key], full, parent), "fold_sha256": fold}
    key_c = "reliable_k3_with_activity_distance_and_contact"
    summary[key_c] = parent_loci(summary[key_c], full_c, parent)
    reasons = Counter(r["reason"] for r in recs)
    summary["where_coverage_is_lost"] = [
        {"step": "the endpoint's units", "pairs": len(us)},
        {
            "step": "a cached answer overlaps the pair's element",
            "pairs": len(us) - reasons["no_cached_answer_overlaps"],
            "lost": {"no_cached_answer_overlaps": reasons["no_cached_answer_overlaps"]},
        },
        {
            "step": "an overlapping answer lists the pair's gene, with a same-cell value",
            "pairs": summary["same_cell_value_any_row"],
            "lost": {k: v for k, v in reasons.items() if k.startswith("gene_not_listed")},
        },
        {
            "step": "a reliable contrast at k = 3, with activity and distance present",
            "pairs": len(full),
            "lost": {
                k: v
                for k, v in reasons.items()
                if k not in ("reliable", "no_cached_answer_overlaps") and not k.startswith("gene_not_listed")
            },
        },
    ]
    return summary


def unit_key(u: ho.Unit) -> str:
    return f"{u.chrom}\t{u.start}\t{u.end}\t{u.gene_id or u.gene}"


def fold_sha(units: list[ho.Unit]) -> str:
    """sha256 of a fold's pair identities (chromosome, element, gene), sorted; no outcome enters it."""
    return hashlib.sha256("\n".join(sorted(unit_key(u) for u in units)).encode()).hexdigest()


# --- duplication -----------------------------------------------------------------------------------
def r_squared(y: Any, x: Any) -> float | None:
    import numpy as np

    y = np.asarray(y, dtype=float)
    if len(y) < 3 or float(np.var(y)) == 0.0:
        return None
    x = np.asarray(x, dtype=float)
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    res = y - x @ beta
    return float(1.0 - (res @ res) / ((y - y.mean()) @ (y - y.mean())))


def compiled_design(descs: list[dict[str, Any]]) -> list[list[float]]:
    """DUPLICATION_RULE D1's design: an intercept, then per named class (action x model cell is K562) an
    indicator and a slope on strength."""
    classes = [(a, s) for a in ("activates", "inhibits") for s in (True, False)]
    rows = []
    for d in descs:
        r = [1.0]
        for a, s in classes:
            on = bool(d.get("named")) and d.get("action") == a and d.get("cell_is_same") == s
            r += [float(on), float(d["strength"]) if on else 0.0]
        rows.append(r)
    return rows


def spearman(a: list[float], b: list[float]) -> float | None:
    if len(a) < 3:
        return None
    v = float(ho.spearman(list(a), list(b), [1.0] * len(a))[0])
    return None if v != v else v


def _r(x: float | None, n: int = 4) -> float | None:
    return None if x is None else round(x, n)


def duplication(recs: list[dict[str, Any]]) -> dict[str, Any]:
    """DUPLICATION_RULE on the reliable complete pairs (k = 3), and what the contrast looks like."""
    import numpy as np

    rows = [r for r in recs if r.get("reason") == "reliable" and r["contrast"].get(3) is not None]
    c = [r["contrast"][3] for r in rows]
    s = [r["same"] for r in rows]
    m = [r["other_median"] for r in rows]
    descs = [r["compiled"] for r in rows]
    named = [i for i, d in enumerate(descs) if d.get("named")]
    r2 = r_squared(c, compiled_design(descs)) if rows else None
    r2_named = r_squared([c[i] for i in named], compiled_design([descs[i] for i in named])) if named else None
    rho = spearman(c, s)
    d1 = bool(r2 is not None and r2 >= DUPLICATE_R2)
    d2 = bool(rho is not None and abs(rho) >= DUPLICATE_SPEARMAN)
    ca = np.asarray(c, dtype=float)
    q = np.quantile(ca, [0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99]).tolist() if c else []
    k562_act = [i for i in named if descs[i].get("cell_is_same") and descs[i].get("action") == "activates"]
    return {
        "pairs": len(rows),
        "d1_r2_on_compiled_target": _r(r2),
        "d1_fires": d1,
        "r2_within_named_pairs_beside": _r(r2_named),
        "named_pairs": len(named),
        "named_share": round(len(named) / len(rows), 4) if rows else None,
        "named_by_class": dict(
            Counter(
                f"{descs[i]['action']}, model cell {'K562' if descs[i]['cell_is_same'] else 'other'}"
                for i in named
            )
        ),
        "d2_spearman_contrast_same_cell": _r(rho),
        "d2_fires": d2,
        "spearman_other_median_same_cell": _r(spearman(m, s)),
        "contrast_quantiles_1_10_25_50_75_90_99": [round(x, 5) for x in q],
        "median_abs_contrast": round(float(np.median(np.abs(ca))), 5) if c else None,
        "median_abs_same_cell": round(float(np.median(np.abs(np.asarray(s)))), 5) if s else None,
        "share_within_rounding": round(float(np.mean(np.abs(ca) <= ROUNDING)), 4) if c else None,
        "share_contrast_exactly_zero": round(float(np.mean(ca == 0)), 4) if c else None,
        "k562_activating_links_with_negative_contrast": (
            f"{sum(c[i] < 0 for i in k562_act)} of {len(k562_act)}" if k562_act else "none"
        ),
        "duplicate": d1 or d2,
    }


# --- the assessment --------------------------------------------------------------------------------
def assess(semantics: dict[str, Any], sets: dict[str, dict[str, Any]], dup: dict[str, Any]) -> dict[str, Any]:
    """GO_RULE, as written above."""
    p = sets.get(f"{PRIMARY[0]} {PRIMARY[1]}") or {}
    within = dup.get("share_within_rounding")
    g1 = bool(
        (semantics.get("matched_answers_uncensored_share") or 0) >= UNCENSORED_SHARE
        and within is not None
        and within < RESOLUTION_SHARE
    )
    full = p.get("reliable_k3_with_activity_and_distance") or {}
    g2 = bool(full.get("clears_floors") and (p.get("coverage_reliable_k3_full") or 0) >= USEFUL_COVERAGE)
    g3 = not dup.get("duplicate", True)
    order = (("no_go_semantics", g1), ("no_go_duplicate", g3), ("no_go_completeness", g2))
    fails = [n for n, ok in order if not ok]
    reading = fails[0] if fails else "go"
    return {
        "g1_semantics": g1,
        "g2_completeness": g2,
        "g3_not_duplicate": g3,
        "contact_on_the_same_pairs": bool(
            (p.get("coverage_reliable_k3_full_contact") or 0) >= USEFUL_COVERAGE
        ),
        "failed": fails,
        "reading": reading,
        "reading_text": READINGS[reading],
        "clear_go": reading == "go",
    }


# --- the run ---------------------------------------------------------------------------------------
def endpoint_units() -> dict[str, list[ho.Unit]]:
    units = ho.crispri_units()
    out: dict[str, list[ho.Unit]] = {}
    for src, part, chroms, _role in ENDPOINT_SETS:
        assert src in pot.ADMISSIBLE_SOURCES, src
        keep = set(chroms)
        us = [u for u in units.get(src, []) if u.chrom in keep and ho.target(u, ENDPOINT) is not None]
        out[f"{src} {part}"] = sorted(us, key=lambda u: (u.chrom, u.start, u.end, u.gene))
    return out


def run() -> tuple[dict[str, Any], dict[str, Any]]:
    t0, c0 = time.time(), time.process_time()
    sets = endpoint_units()
    for us in sets.values():
        assert all(u.cell == SAME_CELL for u in us), "an admissible unit outside K562"
    links = ho.compiled_links()
    idx = ho._Index(x for v in links.values() for x in v)  # the harness's shared-base index
    wanted = {x[0] for x in ENDPOINT_SETS}
    act = activity_index({s: r for s, r in pot.benchmark_rows().items() if s in wanted})
    contact_src = hc.ContactSource(readers={CONTACT_CELL: None})  # cache only: no matrix is ever opened
    by_chrom: dict[str, list[tuple[str, ho.Unit]]] = defaultdict(list)
    for key, us in sets.items():
        for u in us:
            by_chrom[u.chrom].append((key, u))
    records: dict[str, list[dict[str, Any]]] = {k: [] for k in sets}
    cache_stats: dict[str, Any] = {}
    loose_read: list[str] = []
    archives_read: list[str] = []
    answers_seen: dict[str, dict[str, Any]] = {}
    self_check: dict[str, Any] = {}
    for chrom in sorted(by_chrom, key=ho.CHROMS.index):
        pairs = by_chrom[chrom]
        hits, st = cached_answers(chrom, overlap_filter([(u.start, u.end) for _, u in pairs]))
        loose_read += st.pop("loose_paths")
        if et.archive_path(chrom).exists():
            archives_read.append(str(et.archive_path(chrom)))
        cache_stats[chrom] = st
        if chrom == "chr21" and et.archive_path(chrom).exists():  # the smallest archive: stream == json.load
            with gzip.open(et.archive_path(chrom), "rt") as fh:
                whole = json.load(fh)
            same = dict(stream_answers(et.archive_path(chrom))) == whole
            self_check = {"archive": "chr21", "entries": len(whole), "stream_equals_json_load": same}
            del whole
        ordered = sorted(hits.values(), key=lambda h: int(h["start"]))
        hstarts = [int(h["start"]) for h in ordered]
        hreach = max((int(h["end"]) - int(h["start"]) for h in ordered), default=0)
        for key, u in pairs:
            lo = bisect.bisect_left(hstarts, u.start - hreach)
            hi = bisect.bisect_left(hstarts, u.end)
            near = [h for h in ordered[lo:hi] if int(h["end"]) > u.start and int(h["start"]) < u.end]
            rec = pair_record(u, near)
            rec["unit"] = u
            rec["activity"] = act.get((key.split(" ")[0], u.chrom, u.start, u.end, u.gene), False)
            rec["distance"] = u.tss is not None
            c = contact_src.contact(CONTACT_CELL, u.chrom, u.mid, u.tss) if u.tss is not None else None
            rec["contact"] = c is not None
            rec["contact_zero"] = bool(c is not None and c.get("observed", 0) <= 0)
            rec["compiled"] = compiled_description(idx, u)
            records[key].append(rec)
            for h in near:
                answers_seen[f"{chrom}:{h['id']}"] = {
                    "uncensored": uncensored(h),
                    "from": h.get("_from", ""),
                    "model_version": et.model_version_of(h),
                    "tracks": h.get("tracks"),
                }
        del hits, ordered
    body = summarise(sets, records, cache_stats, answers_seen, self_check)
    body["contact_cache_misses_not_fetched"] = contact_src.misses
    body["contact_matrices_opened"] = sum(m is not None for m in contact_src._open.values())
    return body, {
        "archives_read": archives_read,
        "loose_read": loose_read,
        "wall_seconds": round(time.time() - t0, 1),
        "cpu_seconds": round(time.process_time() - c0, 1),
    }


def set_summary(
    src: str, part: str, chroms: tuple[str, ...], role: str, us: list[ho.Unit], recs: list[dict]
) -> dict:
    n = len(us)
    same_any = sum(r["any_k"] >= 0 for r in recs)
    rel = [r for r in recs if r["reason"] == "reliable"]
    full = [r for r in rel if r["contrast"][3] is not None and r["activity"] and r["distance"]]
    full_c = [r for r in full if r["contact"]]
    return {
        "source": src,
        "partition": part,
        "role": role,
        "chromosomes": list(chroms),
        "fold_sha256": fold_sha(us),
        "units": n,
        "loci": max(ho.loci(us)) + 1 if us else 0,
        "units_by_split": dict(Counter(u.split for u in us)),
        "units_by_cell": dict(Counter(u.cell for u in us)),
        "with_an_overlapping_answer": sum(r["answers_overlapping"] > 0 for r in recs),
        "gene_listed_in_an_overlapping_answer": sum(bool(r.get("answers_matched")) for r in recs),
        "same_cell_value_any_row": same_any,
        "same_cell_and_at_least_k_others_any_row": {k: sum(r["any_k"] >= k for r in recs) for k in KS},
        "reliable_same_cell_value": len(rel),
        "reliable_and_at_least_k_others": {k: sum(r["contrast"][k] is not None for r in rel) for k in KS},
        "abstention_reasons": dict(Counter(r["reason"] for r in recs)),
        "activity_present": sum(r["activity"] for r in recs),
        "distance_present": sum(r["distance"] for r in recs),
        "contact_present": sum(r["contact"] for r in recs),
        "contact_observed_zero": sum(r["contact_zero"] for r in recs),
        "compiled_target_names_the_gene": sum(r["compiled"]["named"] for r in recs),
        "reliable_k3_with_activity_and_distance": {
            **floors([r["unit"] for r in full]),
            "fold_sha256": fold_sha([r["unit"] for r in full]),
        },
        "reliable_k3_with_activity_distance_and_contact": floors([r["unit"] for r in full_c]),
        "coverage_reliable_k3_full": round(len(full) / n, 4) if n else None,
        "coverage_reliable_k3_full_contact": round(len(full_c) / n, 4) if n else None,
        "coverage_same_cell_any_row": round(same_any / n, 4) if n else None,
        "coverage_contact": round(sum(r["contact"] for r in recs) / n, 4) if n else None,
    }


def summarise(
    sets: dict[str, list[ho.Unit]],
    records: dict[str, list[dict[str, Any]]],
    cache_stats: dict[str, Any],
    answers_seen: dict[str, dict[str, Any]],
    self_check: dict[str, Any],
) -> dict[str, Any]:
    all_recs = [r for rs in records.values() for r in rs]
    matched = [r for r in all_recs if r.get("answers_matched")]
    n_matched = sum(r["answers_matched"] for r in matched)
    n_matched_unc = sum(r["matched_uncensored"] for r in matched)
    semantics = {
        **SEMANTICS_FROM_CODE,
        "answers_overlapping_a_pair": len(answers_seen),
        "answers_uncensored": sum(a["uncensored"] for a in answers_seen.values()),
        "answers_from": dict(Counter(a["from"] for a in answers_seen.values())),
        "answers_by_model_version": dict(Counter(a["model_version"] for a in answers_seen.values())),
        "answers_by_track_count_top5": [
            list(x) for x in Counter(a["tracks"] for a in answers_seen.values()).most_common(5)
        ],
        "matched_answers": n_matched,
        "matched_answers_uncensored": n_matched_unc,
        "matched_answers_uncensored_share": round(n_matched_unc / n_matched, 4) if n_matched else None,
        "pairs_by_status_of_the_best_overlapping_matched_row": dict(Counter(r["reason"] for r in matched)),
        "pairs_matched_by": dict(Counter(r["matched_by"] for r in matched)),
        "stream_self_check": self_check,
    }
    out_sets = {
        key: set_summary(src, part, chroms, role, sets[key], records[key])
        for (src, part, chroms, role), key in zip(ENDPOINT_SETS, records, strict=True)
    }
    for key in records:  # LOCI_NOTE, added after the provisional run
        with_parent_loci(out_sets[key], sets[key], records[key])
    complete = [
        r["unit"]
        for r in all_recs
        if r["reason"] == "reliable" and r["contrast"][3] is not None and r["activity"] and r["distance"]
    ]
    codes = ho.loci(complete)
    by_locus: dict[int, set[str]] = defaultdict(set)
    for u, code in zip(complete, codes, strict=True):
        by_locus[code].add(u.source)
    all_units = [u for us in sets.values() for u in us]
    refs = ho.crispri_references()
    provenance = {
        "references": {s: sorted(refs.get(s, ())) for s in sorted({x[0] for x in ENDPOINT_SETS})},
        "siblings": {s: list(ho.siblings(s, refs)) for s in sorted({x[0] for x in ENDPOINT_SETS})},
        "pairs_all_sets": len(all_units),
        "independent_loci_all_sets": max(ho.loci(all_units)) + 1 if all_units else 0,
        "pairs_reliable_k3_all_sets": len(complete),
        "independent_loci_reliable_k3_all_sets": max(codes) + 1 if codes else 0,
        "loci_holding_both_studies": sum(len(v) > 1 for v in by_locus.values()),
        "locus_rule": ho.LOCUS_RULE,
    }
    # LOCI_NOTE, added after the provisional run: the counts above that were recomputed on the filtered
    # pairs move to a descriptive entry, and the independent-locus counts are taken on the parent universe
    parent = dict(zip(all_units, ho.loci(all_units), strict=True))
    by_parent: dict[int, set[str]] = defaultdict(set)
    for u in all_units:
        by_parent[parent[u]].add(u.source)
    complete_parents = {parent[u] for u in complete}
    provenance["descriptive_recomputed_after_filtering"] = {
        "loci_reliable_k3_all_sets": provenance.pop("independent_loci_reliable_k3_all_sets"),
        "loci_holding_both_studies_among_reliable_k3_pairs": provenance.pop("loci_holding_both_studies"),
    }
    provenance["independent_loci_all_sets_parent_universe"] = len(by_parent)
    provenance["independent_loci_reliable_k3_all_sets"] = len(complete_parents)
    provenance["loci_holding_both_studies"] = sum(len(v) > 1 for v in by_parent.values())
    provenance["loci_holding_both_studies_with_a_reliable_k3_pair"] = sum(
        len({u.source for u in complete if parent[u] == c}) > 1 for c in complete_parents
    )
    provenance["locus_note"] = LOCI_NOTE
    dup = duplication(all_recs)
    return {
        "semantics": semantics,
        "cache_reads": cache_stats,
        "completeness": out_sets,
        "provenance_and_loci": provenance,
        "duplication": dup,
        "verdict": assess(semantics, out_sets, dup),
        "reading_notes": reading_notes(out_sets, dup, assess(semantics, out_sets, dup), semantics),
    }


def reading_notes(
    sets: dict[str, dict[str, Any]], dup: dict[str, Any], verdict: dict[str, Any], semantics: dict[str, Any]
) -> dict[str, Any]:
    """READING_NOTES with this run's numbers. Beside the verdict; the gate is not read from here."""
    p = sets[f"{PRIMARY[0]} {PRIMARY[1]}"]
    full = p["reliable_k3_with_activity_and_distance"]["units"]
    steps = {x["step"]: x for x in p["where_coverage_is_lost"]}
    covered = p["same_cell_value_any_row"]
    lost_before = {k: v for x in p["where_coverage_is_lost"][1:3] for k, v in x.get("lost", {}).items() if v}
    out = dict(READING_NOTES)
    if not verdict["clear_go"] and "no_go_completeness" in verdict["failed"]:
        out["coverage_wording"] = (
            f"insufficient coverage for this registered complete-case design: {full:,} of {p['units']:,} "
            f"primary pairs ({p['coverage_reliable_k3_full']:.2%}) carry a reliable complete contrast, "
            f"against the registered {USEFUL_COVERAGE:.0%}. {full:,} of the {covered:,} pairs that hold a "
            "same-cell value also hold a reliable contrast, so most of the loss predates the contrast: "
            + ", ".join(
                f"{v:,} {k.replace('_', ' ')}" for k, v in sorted(lost_before.items(), key=lambda kv: -kv[1])
            )
        )
    out["where_coverage_is_lost_primary"] = [steps[k]["pairs"] for k in steps]
    out["d1_margin"] = (
        f"D1 did not fire: the compiled target's description explains {dup['d1_r2_on_compiled_target']} of "
        f"the contrast's variance against the registered {DUPLICATE_R2:.2f}. The variance sits in the tail "
        f"the compiled target names ({dup['named_pairs']:,} of {dup['pairs']:,} pairs; "
        f"{dup['k562_activating_links_with_negative_contrast']} K562-activating links have a negative "
        f"contrast), and the Spearman correlation with the same-cell value is "
        f"{dup['d2_spearman_contrast_same_cell']} (D2 bar {DUPLICATE_SPEARMAN:.2f})"
    )
    out["answers_by_model_version"] = semantics["answers_by_model_version"]
    return out


def manifest(parameters: dict[str, Any], reads: dict[str, Any], versions: dict[str, int]) -> dict[str, Any]:
    inputs = []
    for name in ms.CRISPRI_FILES:
        p = ms.CRISPRI_KNOWLEDGE / name
        if p.exists():
            inputs.append(mf.input_entry(p, partition=ms.CRISPRI_SPLIT_OF[name]))
    inputs += [mf.input_entry(p, partition=None) for p in reads["archives_read"]]
    if reads["loose_read"]:
        inputs.append(
            mf.files_entry("data/knowledge/alphagenome/elements/<chrom>/*.json", reads["loose_read"])
        )
    cpath = hc.KNOWLEDGE / f"{CONTACT_CELL}_{hc.MATRICES[CONTACT_CELL]}_{hc.RESOLUTION}.json"
    if cpath.exists():
        inputs.append(mf.input_entry(cpath, partition=None))
    inputs += [mf.input_entry(p, partition=None) for p in ho.compiled_programs()]
    return {
        "sources": [
            {
                "accession": "EngreitzLab/CRISPR_comparison EPCrisprBenchmark (Gschwind et al., Nature 2026)",
                "version": "pinned by sha256",
            },
            {
                "accession": "AlphaGenome per-element deletion answers (local cache)",
                "version": "pinned by sha256",
            },
            {
                "accession": f"4DN in-situ Hi-C {hc.MATRICES[CONTACT_CELL]} (K562), cached answers only",
                "version": "pinned by sha256",
            },
            {
                "accession": "this repository: data/knowledge/compiled/noncoding_<chrom>.bio",
                "version": "pinned by sha256",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": parameters,
        "exclusions": [
            "every study of the CRISPRi held-out file and Nasser2021: inadmissible under lane-prior's rule "
            "(a48e9e8); the held-out file is parsed by the harness's loader and no unit of it is used",
            "Schraivogel2020 on the fresh partition: it has no unit there",
            "CRISPRi increases, underpowered nulls and missing effects: outside the decrease endpoint",
            "no labelling scored; the only outcome read is the floor check on the complete pairs",
            "no AlphaGenome request; 4DN contacts read from the local cache only, misses counted and not "
            "fetched",
        ],
        "partitions": {
            "fresh": pot.PARTITIONS["fresh"],
            "seen": pot.PARTITIONS["seen"],
            ms.TRAINING: "the CRISPRi benchmark's training file (K562), which holds both admissible studies",
            ms.HELDOUT: "the benchmark's held-out file: parsed by the loader, no unit used",
        },
        # the answers this result reads, per requested model version (review R9 follow-up)
        "model_dependencies": [mf.answers_model_dependency("alphagenome", versions)],
    }


def registration() -> dict[str, Any]:
    return {
        "written": CHECKPOINT,
        "same_cell": SAME_CELL,
        "other_cells": list(OTHER_CELLS),
        "ks": list(KS),
        "full_tracks": FULL_TRACKS,
        "endpoint_sets": [
            {"source": s, "partition": p, "chromosomes": list(c), "role": r} for s, p, c, r in ENDPOINT_SETS
        ],
        "admissibility": "lane-prior, prior_only_test.STUDIES (a48e9e8), read and not changed",
        "reliable_rule": RELIABLE_RULE,
        "contrast_rule": CONTRAST_RULE,
        "pair_rule": PAIR_RULE,
        "contact_rule": CONTACT_RULE,
        "activity_rule": ACTIVITY_RULE,
        "duplication_rule": DUPLICATION_RULE,
        "go_rule": GO_RULE,
        "readings": READINGS,
        "seen_before": list(SEEN_BEFORE),
        "budget": BUDGET,
    }


def peak_rss_mb() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return round(r / (1024 * 1024) if sys.platform == "darwin" else r / 1024, 1)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.parse_args(argv)
    body, reads = run()
    versions = body["semantics"]["answers_by_model_version"]
    payload = {
        "question": "can the per-element cached AlphaGenome deletion answer support one added feature, the "
        "same-cell signed deletion minus the median other-cell signed deletion, on lane-prior's admissible "
        "CRISPRi decrease endpoints, without duplicating the compiled target or the same-cell deletion",
        "status": STATUS,
        "may_be_called": MAY_BE_CALLED,
        "may_not_be_called": list(MAY_NOT_BE_CALLED),
        "registration": registration(),
        **body,
        "every_dataset_read_is_development_evidence": True,
        "alphagenome_requests": 0,
        "downloads": 0,
        "compute": {
            "wall_seconds": reads["wall_seconds"],
            "cpu_seconds": reads["cpu_seconds"],
            "peak_rss_mb": peak_rss_mb(),
            "cache_archives_streamed": len(reads["archives_read"]),
            "loose_answer_files_read": len(reads["loose_read"]),
            "measured_with": "time.time, time.process_time and getrusage, this process; the archives are "
            "streamed a chunk at a time by one reader, and chr21 is also loaded whole once for the "
            "self-check",
        },
    }
    params = {
        "checkpoint": CHECKPOINT,
        "same_cell": SAME_CELL,
        "other_cells": list(OTHER_CELLS),
        "ks": list(KS),
        "full_tracks": FULL_TRACKS,
        "low": LOW,
        "rounding": ROUNDING,
        "half_window": HALF_WINDOW,
        "useful_coverage": USEFUL_COVERAGE,
        "duplicate_r2": DUPLICATE_R2,
        "duplicate_spearman": DUPLICATE_SPEARMAN,
        "uncensored_share": UNCENSORED_SHARE,
        "resolution_share": RESOLUTION_SHARE,
        "fresh_chromosomes": list(pot.FRESH),
        "seen_chromosomes": list(pot.SEEN),
        "floors": {"positives": ho.MIN_POSITIVES, "negatives": ho.MIN_NEGATIVES, "loci": ho.MIN_LOCI},
    }
    p = save_result(NAME, payload, manifest=manifest(params, reads, versions))
    v = payload["verdict"]
    print("reading", v["reading"], "failed", v["failed"], flush=True)
    print("saved", p, "peak_rss_mb", payload["compute"]["peak_rss_mb"], flush=True)


if __name__ == "__main__":
    main()
