# SPDX-License-Identifier: AGPL-3.0-or-later
"""Re-read every scored frame through the current readers, and the registration that licenses it.

Two defects were confirmed by two different frames and neither was fixed, because both sit in code
the whole benchmark is scored through and the sessions that found them said so rather than reaching
for it.

**One: `loci.read_deletion` breaks at the first key that holds a gene.** The sweep stores two views
of one deletion - `predicted`, the strongest expression change over any gene, and
`predicted_coding`, the strongest one restricted to protein-coding genes. The reader walked
``("predicted_coding", "predicted")`` and ``break``-ed at the first key with a gene, so the
any-gene answer was visible only where no coding answer existed at all. Section 18 of
docs/LOCI-BENCHMARK.md recorded this at the H19 ICR as a caveat that changed no verdict. Section 23
caught it with a positive: at the chr8 CCDC26 element the stored row names the published lncRNA at
-1.691 and the reader reported GSDMC, a protein-coding gene 204 kb away, at -0.1433 - a twelfth of
the effect - and the locus scored as a miss.

**Two: the draw of the fourth frame matches gene symbols against GENCODE symbols**, and the ENCODE
benchmark's `measuredGeneSymbol` column is as old as the screens that filled it. Section 23 joined
the `measuredGeneEnsemblId` the same file already carries and found SSFA2 is ITPRID2, SARS is SARS1
and WDR61 is SKIC8. Eight held-out elements pass every step of `loci_fourth`'s own rule under the
current symbols and are missing from its draw, so "the frame is exhausted" was true of the rule as
implemented and not of the rule as written.

**Why this file exists rather than a patch.** `read_deletion` is the reader every frame is scored
through, so changing it moves sections 9, 19, 21, 22 and 23 at once, and the session that makes the
change is the session that decides which of the new numbers get reported. That is the shape of
error this project has made before. So the reading, the expected direction of every affected rate,
and what would make the change wrong rather than right are written into `PREREGISTRATION` below and
committed before a single number is recomputed.

**What it costs: zero AlphaGenome requests, by construction.** `read_deletion` makes no request and
no network call; it ranks rows the finished all-element sweep had already bought. The eight elements
the corrected draw recovers were checked against `loci._deletion_rows` before this was written and
all eight are already covered by `enhancer_targets_all`. A re-read that spends a request is a bug in
this module, and `reread_frame` asserts it read no new interval.

**The re-read is not a re-score.** It recomputes one layer from rows that are already on disk and
runs `loci.score_locus`, the shipped scorer, over the result. `historical_deletion_reading` keeps the
pre-2026-09-21 loop verbatim, and `reread_frame` refuses to report anything unless that loop
reproduces every stored reading field for field. If it does, the harness is reading the same rows the
run read and the only difference between before and after is the fix.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from genomeos.benchmark import loci
from genomeos.benchmark.loci import RESULTS_DIR, Expect
from genomeos.results import load_result

#: the frames scored through `loci.read_deletion`, with the doc section each one is reported in.
FRAMES: dict[str, str] = {
    "loci_benchmark": "9 - the seventeen panel",
    "loci_candidates": "19 - the nine candidates",
    "loci_third": "21 - the third set",
    "loci_fourth": "22 - the fourth frame, drawn by a rule",
    "loci_noncoding": "23 - the non-coding frame",
}

PREREGISTRATION: dict[str, Any] = {
    "written": (
        "2026-09-21, before any frame was re-read and before the corrected draw was taken. Committed"
        " on its own, ahead of the one-line change it licenses"
    ),
    "the_two_changes": (
        "ONE, the reader: `loci.read_deletion` ranks both stored keys together instead of breaking"
        " at the first one that holds a gene. TWO, the draw: `loci_fourth.select` joins the"
        " benchmark file's own `measuredGeneEnsemblId` to GENCODE before testing whether a"
        " published target is protein coding"
    ),
    "registered_together_measured_apart": (
        "one registration, two stages, and the reason is that they move different things. The reader"
        " changes the READINGS over a draw that is held fixed, so every denominator in sections 9,"
        " 19, 21, 22 and 23 is the same before and after and a difference is a difference in what a"
        " layer named. The Ensembl join changes the DRAW, so the fourth frame's denominator itself"
        " moves and a rate computed after it is over a different set of loci. Applying both at once"
        " would leave no way to say which of them moved the fourth frame's 0.200, which is the one"
        " number this whole exercise is most likely to be accused of having been aimed at. So the"
        " reader is applied first and every frame is re-read with the draw untouched; only then is"
        " the join applied and the frame re-drawn, and the fourth frame is reported three times:"
        " as it stands, after the reader alone, and after the reader and the join"
    ),
    "what_the_reader_will_do": (
        "for every already-scored element inside the window, both `predicted` and `predicted_coding`"
        " are entered into the same ranking, a gene named by both keys of one element counts once,"
        " and the layer's answer is the gene with the largest absolute predicted log2 fold change"
        " whatever its biotype. On an exact tie the coding entry keeps its place, because it is"
        " inserted first and the sort is stable - the tie-break is deliberately the one that"
        " preserves the old answer. The pre-2026-09-21 reading is computed alongside and returned in"
        " the same dict as `coding_first_target` and `coding_first_targets`, so every rate this"
        " benchmark reported before today stays computable from the shipped reader rather than only"
        " from the git history"
    ),
    "why_this_reading_and_not_the_favourable_one": (
        "FIRST, it is the question the layer says it is answering. `read_deletion`'s own first line"
        " is 'which gene moves when an element here is deleted', and the answer to that question"
        " does not depend on the biotype of the gene that moves. A reader that filters by biotype is"
        " answering 'which CODING gene moves most' while being scored as though it had answered the"
        " first, and the benchmark has no line that says so."
        " SECOND, `predicted_coding` is not a better estimate than `predicted`, it is a restriction"
        " of it: the same deletion, the same scorer, the same run, with the candidate set cut down."
        " It exists because the panel's first loci all had coding targets. Preferring a restriction"
        " of a computation over the computation is a choice that needs a reason, and there is none"
        " on the record."
        " THIRD, and this is what makes it the right reading rather than the convenient one, the"
        " change is expected to be NET UNFAVOURABLE to this benchmark's headline. Because"
        " |predicted| >= |predicted_coding| holds for one element by construction, a published"
        " CODING target that used to rank first can now be displaced by a stronger non-coding"
        " effect, and it can never be promoted. Only a published NON-CODING target can gain. Three"
        " of the five frames contain coding targets only, and the fourth frame's drawing rule"
        " REQUIRES the target to be protein coding, so in four of the five frames the reader fix can"
        " lower the derived rate and cannot raise it. The one frame it can raise is n = 2."
        " FOURTH, the readings that would have been favourable were available and are not taken. The"
        " lenient `among` reading - a published target anywhere in the six a layer names - rises"
        " monotonically under this change and is already excluded from every headline by section 9's"
        " own rule; it stays excluded. Choosing per locus whichever of the two keys names the"
        " published target would read the answer before reporting it and is not a reader at all"
    ),
    "scope_what_is_not_being_changed": (
        "`read_gene_input` reads `predicted_coding` and nothing else and is left alone. Section 23"
        " registered that as a separate property of a separate layer - it is a summed-window reading"
        " rather than an element answer - and folding a second layer into this change would put two"
        " movements behind one number. `summed_by_cell` keeps reading `predicted_coding_by_cell`,"
        " which is a coding-specific field by name. The drawing rules of the first three frames are"
        " curated panels and are not re-drawn here"
    ),
    "expected_direction_per_section": {
        "reading": (
            "written before the re-read. Each entry is the direction the rate is expected to move and"
            " by roughly how much. The headline `target_derived` is a UNION over the derived layers,"
            " so a deletion-layer loss only reaches the headline at a locus where the deletion layer"
            " was the only derived layer that hit"
        ),
        "9 - the seventeen panel": (
            "derived target 15/17 (0.882) and 13/15 where the model could answer. Expected to FALL by"
            " 0 to 2 loci, most likely 0 or 1. The panel's targets are coding except at the H19 ICR,"
            " which is where section 18 first saw the defect and where a gain is possible"
        ),
        "19 - the nine candidates": (
            "derived target 8/9 (0.889). Expected to FALL by 0 or 1 locus. Every published target"
            " here is coding, so a rise is not available"
        ),
        "21 - the third set": (
            "derived target 8/9 (0.889). Expected to FALL by 0 or 1 locus, same reasoning"
        ),
        "22 - the fourth frame": (
            "derived target 10/60 (0.167) and 10/50 (0.200) where the model could answer; the"
            " deletion layer alone is 3/50. Expected to FALL by 0 to 3 loci, most likely 1: the"
            " deletion layer is the only derived hit at some of the ten, and the rule requires a"
            " coding target so no locus here can gain. 0.200 is expected to land between 0.140 and"
            " 0.200. The 'named the nearest coding TSS' count of 26/50 is expected to FALL, because"
            " a non-coding gene displacing the nearest coding TSS moves a locus out of that column"
            " and into 'another gene'"
        ),
        "23 - the non-coding frame": (
            "derived target 0/2 rises to 1/2. This is NOT a prediction: section 23 computed it with"
            " `loci_noncoding.counterfactual` and disclosed that the CCDC26 row had been read before"
            " its own registration was written. It is listed here so that it is not counted twice"
        ),
        "the matched-window comparisons in 8, 19, 21 and 22": (
            "expected to be EXACTLY unchanged, targets and controls alike. `loci.claims` asks whether"
            " a layer named a gene at all, not which one; the fix never empties a reading that held a"
            " gene and never fills one that was empty, because the candidate set only ever grows."
            " `elements_scored` and `summed_by_cell` are likewise expected to be identical. Any"
            " movement in a claim rate, in an element count or in a p-value is a defect in the change"
            " and not a finding"
        ),
        "cell and direction": (
            "both read the rank-1 gene's tissue and action, so both can move in EITHER direction and"
            " no direction is registered for them. They are reported before and after whatever they"
            " do"
        ),
    },
    "known_before_the_registration": (
        "what the fix does at the chr8 CCDC26 element, because section 23 measured it and published"
        " it; and that all eight elements the corrected draw recovers are already covered by the"
        " finished sweep, because that was checked with `loci._deletion_rows` while this was being"
        " written, in order to state the budget. Neither is a prediction this re-read confirms"
    ),
    "falsifiers": {
        "the harness": (
            "`historical_deletion_reading` must reproduce the stored deletion reading of every locus"
            " and every negative in all five result files, field for field. If it does not, the"
            " re-read is not reading what the runs read, and NO number out of it is admissible -"
            " including the ones that would look good"
        ),
        "the mechanism": (
            "a frame whose published targets are all coding must not RISE. |predicted| >="
            " |predicted_coding| for one element, so a coding rank-1 can only be displaced, never"
            " promoted; a rise there would mean the mechanism argued above is false and that the"
            " change is doing something not understood. If one appears, the run stops, the defect"
            " stays pinned, and the rise is reported as an unexplained result rather than as a gain"
        ),
        "double counting": (
            "`elements_scored` must be identical before and after at every window. The fix enters two"
            " keys per element, and an element whose two keys name the same gene must still count"
            " once. A changed element count means the ranking is inflated and every rate from it is"
            " void"
        ),
        "the claim tables": (
            "the four claim rates and their control rates must not move by one window. They are"
            " existence claims and the fix cannot touch them"
        ),
        "unfalsifiability": (
            "if after the fix the deletion layer names a non-coding gene at most loci, then 'names"
            " the published target' stops discriminating in the four coding frames and the layer has"
            " been turned into a different measurement rather than repaired. The threshold registered"
            " in advance: if the rank-1 gene changes at more than half the graded loci of the fourth"
            " frame, the fix is reported as a change of measurement and the old reading keeps the"
            " headline until a frame is registered for the new one"
        ),
    },
    "if_the_fourth_frames_0.200_RISES": (
        "it is reported in the same place, in the same table and in the same sentence as the fall it"
        " would partly undo, and with the same emphasis section 22 gave the fall. BOTH readings stay"
        " in the record permanently: section 22's numbers are not rewritten, the new ones are added"
        " beside them, and every table that carries one carries the other. The registered position is"
        " that a rise here is WEAK evidence and a fall is weak evidence too - n = 50 over one screen"
        " - and that the interesting number is not the headline but which gene the deletion layer"
        " names instead of the published one. If the rise came from the Ensembl join rather than the"
        " reader, it is reported as a change of denominator, because the eight recovered elements are"
        " seven perturbations of one gene in one cell line plus one other, and a rate over 58 loci"
        " that contains seven elements of ITPRID2 is less independent than a rate over 50 that does"
        " not"
    ),
    "if_the_fix_is_wrong": (
        "if the reader fix makes the frames worse in a way the mechanism does not explain, or trips"
        " any falsifier above, the change is reverted, the defect stays pinned by its test, and the"
        " outcome is written up as a failed attempt. That is a legitimate result and it is the reason"
        " this is registered first"
    ),
    "budget": (
        "0 AlphaGenome requests for the whole exercise, and the number is not an estimate. The reader"
        " makes no request. All eight elements the corrected draw recovers already have a scored"
        " deletion inside them from `enhancer_targets_all`. If the re-draw turns out to need more"
        " than 10 requests the run stops and the number goes back to the coordinator unspent"
    ),
}


def historical_deletion_reading(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """`read_deletion`'s ranking loop exactly as it stood until 2026-09-21, kept so it can be checked.

    This is not a fallback and nothing scores through it. It exists so that `reread_frame` can prove
    it is reading the same rows the original runs read: if this loop reproduces every stored reading
    field for field, then the difference between the stored numbers and the new ones is the change to
    `loci.read_deletion` and nothing else - not a different results directory, not a sweep that has
    grown since, not a window read at different coordinates.
    """
    genes: dict[str, dict[str, Any]] = {}
    for e in rows:
        for key in ("predicted_coding", "predicted"):
            p = e.get(key) or {}
            if not p.get("gene"):
                continue
            g = genes.setdefault(
                p["gene"],
                {"gene": p["gene"], "elements": 0, "best": 0.0, "action": p["action"], "tissue": p["tissue"]},
            )
            g["elements"] += 1
            if abs(p["log2_fold_change"]) > abs(g["best"]):
                g["best"] = p["log2_fold_change"]
                g["action"] = p["action"]
                g["tissue"] = p["tissue"]
            break
    ranked = sorted(genes.values(), key=lambda g: -abs(g["best"]))
    return {
        "elements_scored": len(rows),
        "targets": ranked[:6],
        "target": ranked[0]["gene"] if ranked else None,
        "action": ranked[0]["action"] if ranked else None,
        "tissue": ranked[0]["tissue"] if ranked else None,
    }


def _expect_from(stored: dict[str, Any]) -> Expect:
    """The `Expect` a saved locus row was scored against. `as_dict` is field for field but `length`."""
    return Expect(**{k: v for k, v in stored.items() if k != "length"})


def _comparable(reading: dict[str, Any]) -> dict[str, Any]:
    """The fields of a deletion reading that the historical loop is responsible for."""
    return {
        "elements_scored": reading.get("elements_scored"),
        "targets": reading.get("targets"),
        "target": reading.get("target"),
        "action": reading.get("action"),
        "tissue": reading.get("tissue"),
    }


def named_instead(rows: list[dict[str, Any]], layer: str = "deletion") -> dict[str, int]:
    """Section 22's own table for one layer: the published target, the trap, another gene, nothing.

    The categories are the frame's, not new ones: `nearest_gene_trap` is the element's nearest coding
    TSS, which the fourth frame's drawing rule guarantees is NOT a published target.
    """
    tally = {"the published target": 0, "the nearest coding TSS": 0, "another gene": 0, "nothing": 0}
    for r in rows:
        named = (r["score"]["scored"]["target"]["by_layer"].get(layer) or {}).get("named") or []
        first = named[0] if named else None
        want = set(r["expected"]["targets"])
        trap = r["expected"].get("nearest_gene_trap")
        if first is None:
            tally["nothing"] += 1
        elif first in want:
            tally["the published target"] += 1
        elif trap and first == trap:
            tally["the nearest coding TSS"] += 1
        else:
            tally["another gene"] += 1
    return tally


def layer_hits(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Per layer, how often it ranks a published target first. The same `hit` the headline uses."""
    out: dict[str, dict[str, Any]] = {}
    for r in rows:
        for layer, v in r["score"]["scored"]["target"]["by_layer"].items():
            slot = out.setdefault(layer, {"k": 0, "n": 0, "provenance": v.get("provenance")})
            slot["n"] += 1
            slot["k"] += int(bool(v.get("hit")))
    for slot in out.values():
        slot["rate"] = round(slot["k"] / slot["n"], 3) if slot["n"] else None
    return out


def neighbourhood(results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """Where each gene the repaired reader promoted sits relative to the target it displaced.

    Read off the saved `loci_reread` movements and GENCODE, free. It answers the question the
    movements raise on sight: the reader stopped naming LCT, SOST, MC1R and HBA2 and started naming
    Y_RNA, LINC02594, ENSG00000259006 and ENSG00000290010, and whether that is the model pointing
    somewhere else or the model pointing at the same place under a non-coding name is a fact about
    GENCODE coordinates rather than a matter of opinion.
    """
    saved = load_result("loci_reread", results_dir) or {}
    chroms: dict[str, loci.Chromosome] = {}
    per: list[dict[str, Any]] = []
    try:
        for frame, v in (saved.get("frames") or {}).items():
            where_each = {
                r["locus"]: r["expected"]["chrom"]
                for r in (load_result(frame, results_dir) or {}).get("loci", [])
            }
            for m in v.get("movements") or []:
                chrom = where_each.get(m["locus"])
                if chrom is None:
                    continue
                ch = chroms.get(chrom) or chroms.setdefault(chrom, loci.Chromosome(chrom, results_dir))
                copies = [x for x in ch.annotation.genes.values() if x.symbol == m["now"]]
                promoted = ch.annotation.genes.get(m["now"]) or (copies[0] if copies else None)
                where = "no GENCODE record on this chromosome for the promoted name"
                gap = None
                for t in m["published"]:
                    pub = ch.annotation.genes.get(t) or next(
                        (x for x in ch.annotation.genes.values() if x.symbol == t), None
                    )
                    if promoted is None or pub is None:
                        continue
                    if promoted.locus.start < pub.locus.end and promoted.locus.end > pub.locus.start:
                        where, gap = "overlaps the published target's gene body", 0
                        break
                    d = min(
                        abs(promoted.locus.start - pub.locus.end), abs(pub.locus.start - promoted.locus.end)
                    )
                    if gap is None or d < gap:
                        where, gap = "outside it", d
                per.append(
                    {
                        "frame": frame,
                        "locus": m["locus"],
                        "was": m["was"],
                        "now": m["now"],
                        "published": m["published"],
                        "where": where,
                        "gap_bp": gap,
                        # Y_RNA, and symbols like it, name many genes on one chromosome. A distance
                        # to "the first one" is arithmetic about the wrong copy, so it is flagged
                        # rather than quoted, and such a row is excluded from the counts below.
                        "copies_of_this_symbol_on_the_chromosome": len(copies),
                        "distance_is_meaningful": bool(promoted is not None and len(copies) <= 1),
                    }
                )
    finally:
        for ch in chroms.values():
            ch.close()
    usable = [p for p in per if p["distance_is_meaningful"]]
    overlapping = [p for p in usable if p["where"].startswith("overlaps")]
    within_100k = [p for p in usable if p["gap_bp"] is not None and 0 < p["gap_bp"] <= 100_000]
    return {
        "movements": len(per),
        "locatable_unambiguously": len(usable),
        "not_locatable": [
            {"locus": p["locus"], "now": p["now"], "copies": p["copies_of_this_symbol_on_the_chromosome"]}
            for p in per
            if not p["distance_is_meaningful"]
        ],
        "the_promoted_gene_overlaps_the_published_target": len(overlapping),
        "outside_it_but_within_100_kb": len(within_100k),
        "per_movement": per,
        "reading": (
            "a promoted gene that overlaps the target it displaced means the model named the right"
            " PLACE and the coding-only reader was reporting a different gene at that place; one far"
            " away means the repaired reader moved the answer somewhere else. The two are not the"
            " same finding and the split is what says which happened. The counts are over the"
            " movements whose promoted gene can be located unambiguously: a symbol GENCODE uses at"
            " many loci on one chromosome, Y_RNA above all, gives a distance to whichever copy comes"
            " first, which is arithmetic about the wrong gene, and those rows are listed apart"
        ),
    }


def reread_frame(name: str, results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """Re-read one saved frame's deletion layer and re-score it. No request, no network.

    Returns the stored numbers, the new numbers, the per-locus movements and the fidelity check.
    Raises if the historical loop fails to reproduce a stored reading, because a harness that is not
    reading the same rows cannot be used to report a change in them.
    """
    out = load_result(name, results_dir)
    if not out:
        raise FileNotFoundError(f"no {name} result to re-read")
    chroms: dict[str, loci.Chromosome] = {}
    moved: list[dict[str, Any]] = []
    unfaithful: list[dict[str, Any]] = []
    after: list[dict[str, Any]] = []
    try:
        for row in out["loci"]:
            expect = _expect_from(row["expected"])
            ch = chroms.get(expect.chrom) or chroms.setdefault(
                expect.chrom, loci.Chromosome(expect.chrom, results_dir)
            )
            s, en = expect.element
            rows = loci._deletion_rows(expect.chrom, s, en, results_dir)
            stored = row["readings"]["deletion"]
            replay = historical_deletion_reading(rows)
            if _comparable(replay) != _comparable(stored):
                unfaithful.append({"locus": row["locus"], "stored": _comparable(stored), "replay": replay})
            fresh = loci.read_deletion(ch, s, en, results_dir)
            # `loci.build` marks a reading unaskable when the published target lies outside the
            # model's input. That is a property of the panel's geometry, not of the reader, and
            # dropping it here would silently change the denominator of every
            # "where the model could answer" rate.
            if stored.get("unaskable"):
                fresh["unaskable"] = stored["unaskable"]
                if fresh.get("pending"):
                    fresh["pending"] = stored["unaskable"]
            readings = {**row["readings"], "deletion": fresh}
            score = loci.score_locus(expect, ch, readings)
            after.append({**row, "readings": readings, "score": score})
            if stored.get("target") != fresh.get("target"):
                moved.append(
                    {
                        "locus": row["locus"],
                        "published": list(expect.targets),
                        "was": stored.get("target"),
                        "was_log2": (stored.get("targets") or [{}])[0].get("best"),
                        "now": fresh.get("target"),
                        "now_log2": (fresh.get("targets") or [{}])[0].get("best"),
                        "hit_before": bool(row["score"]["scored"]["target"]["by_layer"]["deletion"]["hit"]),
                        "hit_after": bool(score["scored"]["target"]["by_layer"]["deletion"]["hit"]),
                    }
                )
        negatives = _reread_negatives(out, chroms, results_dir)
    finally:
        for ch in chroms.values():
            ch.close()
    if unfaithful:
        raise AssertionError(
            f"{name}: the historical loop did not reproduce {len(unfaithful)} stored deletion"
            f" reading(s); the re-read is not reading the rows the run read. First: {unfaithful[0]}"
        )
    return {
        "frame": name,
        "section": FRAMES.get(name, "?"),
        "loci": len(out["loci"]),
        "fidelity": {
            "loci_checked": len(out["loci"]),
            "negatives_checked": negatives["checked"],
            "mismatches": 0,
            "reading": (
                "the pre-2026-09-21 loop reproduced every stored deletion reading, so the rows read"
                " here are the rows the run read and the difference below is the change to the reader"
            ),
        },
        "before": loci.aggregate(out["loci"]),
        "after": loci.aggregate(after),
        # the positive control goes through the same readers, so `loci.build` aggregates it in; every
        # frame that has one reports its headline over the DRAWN loci, and the split is kept here.
        "drawn_only_before": loci.aggregate(
            [r for r in out["loci"] if not r["locus"].startswith("CONTROL_")]
        ),
        "drawn_only_after": loci.aggregate([r for r in after if not r["locus"].startswith("CONTROL_")]),
        "by_layer_before": layer_hits(out["loci"]),
        "by_layer_after": layer_hits(after),
        "named_instead_before": named_instead(out["loci"]),
        "named_instead_after": named_instead(after),
        "rank_one_changed": len(moved),
        "movements": moved,
        "negatives": negatives,
        "loci_after": after,
    }


def _reread_negatives(
    out: dict[str, Any], chroms: dict[str, loci.Chromosome], results_dir: Path
) -> dict[str, Any]:
    """The matched negatives, which carry claims rather than readings. Existence claims only."""
    changed: list[dict[str, Any]] = []
    counts: list[str] = []
    for w in out.get("negatives") or []:
        ch = chroms.get(w["chrom"]) or chroms.setdefault(w["chrom"], loci.Chromosome(w["chrom"], results_dir))
        fresh = loci.read_deletion(ch, w["start"], w["end"], results_dir)
        # the `target` claim is a union with the eQTL layer, so it is tested through `named_a_gene`
        # below rather than against `claims["target"]`, which a negative can hold on the eQTL side
        # alone.
        if bool(fresh.get("action")) != bool(w["claims"]["direction"]):
            changed.append({"window": f"{w['chrom']}:{w['start']}-{w['end']}", "claim": "direction"})
        if bool(fresh.get("elements_scored")) != bool(w["claims"]["deletion_scored"]):
            counts.append(f"{w['chrom']}:{w['start']}-{w['end']}")
        if (w.get("deletion_target") is None) != (fresh.get("target") is None):
            changed.append({"window": f"{w['chrom']}:{w['start']}-{w['end']}", "claim": "named_a_gene"})
    return {
        "checked": len(out.get("negatives") or []),
        "claims_changed": changed,
        "element_counts_changed": counts,
        "reading": (
            "the registration predicted that no claim over a negative window can move, because a"
            " claim asks whether a layer named a gene at all and the fix only ever enlarges the"
            " candidate set. An entry here falsifies that"
        ),
    }


def reread_all(results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """Every frame scored through `loci.read_deletion`, re-read and re-scored. Free."""
    frames = {}
    for name in FRAMES:
        if load_result(name, results_dir) is None:
            continue
        frames[name] = reread_frame(name, results_dir)
    return {
        "result": "loci_reread",
        "preregistration": PREREGISTRATION,
        "frames": {k: {x: y for x, y in v.items() if x != "loci_after"} for k, v in frames.items()},
        "requests_spent": 0,
        "note": (
            "Every frame this benchmark has scored, re-read through the repaired"
            " `loci.read_deletion` with its draw held fixed. No request, no network: the rows were"
            " bought by the finished all-element sweep and are read off disk. The registration that"
            " licensed the change is carried in this file and was committed before any number below"
            " was computed. " + loci.NOTE
        ),
    }
