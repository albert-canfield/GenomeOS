# SPDX-License-Identifier: AGPL-3.0-or-later
"""The last one-target consumers, read from the sweep's per-element cache beside the compact table.

    uv run python scripts/onetarget_consumers.py [--no-save]

Four modules read one target gene per element from the compact `all_elements` tables: the motif
transfer's sweep lead (`motif_transfer.tier_blocks`), the syntax tiling's `moved`, the syntax
candidates' node context (`candidates.node_context`) and the locus benchmark's gene-input layer
(`loci.read_gene_input`). Each is read here twice -- as it reads today and through
`attribution.targets.ElementResponses`, where every gene in the scorer's window is kept -- and the
registration below, committed before any of this was run, says which figures must come back to the
digit, which may move and in which direction, and what would falsify it. 0 AlphaGenome requests: the
cache is on disk, and a chromosome's archive is held one at a time.
"""

from __future__ import annotations

from genomeos.predict.enhancer_target import MIN_EFFECT

REGISTERED = "2026-09-27"

# A gene is "at the bar" at an element when the sweep's larger predicted change for it, the fall or the
# rise over all 371 tracks, reaches MIN_EFFECT. That is `predict_target`'s own rule applied to every
# gene of the window instead of to its head only, so an element has a gene at the bar exactly when the
# compact table names one; what the window adds is the other genes, never a different yes or no.
BAR = MIN_EFFECT

PRE_REGISTRATION = {
    "registered": REGISTERED,
    "one_registration_for": [
        "motif_transfer.tier_blocks (the sweep lead per block)",
        "syntax_tiling.moved (the tiling's own scored windows)",
        "candidates.node_context (node_targets_named, deleted_elements_on_block)",
        "loci.read_gene_input (the gene-input layer of the locus benchmark)",
    ],
    "shared_premise": (
        "The roadmap line says an element that moves a non-top gene reads as 'did not move'. It is "
        "registered here as expected FALSE: the compact table's head is the maximum over the window, so "
        "'some gene in the window is at the bar' and 'the head is at the bar' are the same statement. "
        "Measured on chr21 during the survey: 7,846 elements yes/yes, 4,293 no/no, 0 disagreements. "
        "What a one-target reader loses is the count and identity of the other genes, not the yes/no."
    ),
    "shared_controls": [
        "every element read from both sources: the window's head at the bar names the compact table's "
        "`predicted` gene and its coding head names `predicted_coding` -- 0 disagreements required",
        "an element the cache does not hold keeps its compact reading and is counted by name "
        "(not_cached), never scored as zero",
    ],
    "motif_transfer": {
        "reproduce_exactly": {
            "real_unknown.with_a_tested_element": 531,
            "real_unknown.with_a_sweep_lead": 331,
            "neutral.with_a_tested_element": 1181,
            "neutral.with_a_sweep_lead": 768,
            "by_case.with_a_sweep_lead": {
                "recent": 11,
                "relaxed": 167,
                "syntax": 19,
                "tolerant": 131,
                "unmeasured": 3,
            },
        },
        "reproduce_by_construction": (
            "the agreement AUROC 0.5422 [0.4937, 0.5934]: its label is with_a_sweep_lead per block, so if "
            "every block's label reproduces the AUROC cannot move; the count model is not re-run"
        ),
        "may_change": (
            "new field only, genes at the bar per block over its elements; expected >= the distinct head "
            "genes of the block, about 2.3 genes per moving element (chr21 survey 2.34)"
        ),
        "falsifier": "any block whose lead or tested count differs from the stored result",
    },
    "syntax_tiling": {
        "reproduce_exactly": "nothing: the registered run was cancelled unrun and has no result",
        "not_rerun": "scoring the tiles is model requests, and the key is held elsewhere",
        "change": (
            "`moved` is unchanged; the tiling's rows gain the genes at the bar read from the cache the "
            "scorer writes, so its result will not be one-target when it runs"
        ),
        "falsifier": "a cached element where `moved` and 'some gene at the bar' disagree (unit test)",
    },
    "candidates": {
        "reproduce_exactly": (
            "`node_targets_named`, `deleted_elements_in_node` and `deleted_elements_on_block` of the one "
            "candidate on chr22, the only candidate chromosome whose whole-chromosome run was complete "
            "when syntax_candidates_genome_wide was written on 2026-09-13"
        ),
        "expected_to_differ_from_the_stored_file_by_input": (
            "the other 68, whose chromosomes completed after 2026-09-13: the old reading recomputed today "
            "is the control for the window reading, and its difference from the stored file is attributed "
            "to the input and counted, not to the reader (the confound lane-scoring measured on 7f7c8c1)"
        ),
        "may_change": {
            "node_targets_named_window": (
                "for every gene, window count >= table count (a gene keeps every vote it had as head and "
                "gains one wherever it is at the bar below the head); strictly monotone, any gene whose "
                "window count is lower falsifies the reader"
            ),
            "top_node_gene": "two-sided, no direction registered: counted how many of the 69 change",
        },
        "not_rerun": "the reading, labels and set tests of syntax_candidates_genome_wide stand as stored",
    },
    "loci_gene_input": {
        "reproduce_exactly": (
            "every stored gene_input field of the 17 loci in loci_benchmark (2026-09-17, after the sweep "
            "finished): elements_scored, naming_a_coding_gene, the top-8 genes with their sums, target and "
            "rank_of_first_published_target; stored strict hits 10 of 17 (SHH_ZRS, HERC2_OCA2, MCM6_LCT, "
            "HBB_LCR, FTO_IRX3, BCL11A_enhancer, ABO, APP, HOXD, SOX9_PierreRobin), among 15 of 17"
        ),
        "window_rule": (
            "each element credits |effect| to every coding gene at the bar in its window, not only to its "
            "coding head; elements the cache does not hold (stated intervals, VISTA and lentiMPRA rows) "
            "keep their compact credit"
        ),
        "expected_direction": (
            "strict hits hold or rise (>= 10 of 17): a published target that is second at many elements "
            "gains credit (MYC_8q24 rank 2, PMP22_CMT1A rank 2 are the likeliest to flip); the median "
            "rank of the first published target over the 7 misses holds or improves"
        ),
        "falsifier": (
            "strict hits fall below 10, or any of the four loci gene_input carries alone in the headline "
            "(SHH_ZRS, HERC2_OCA2, FTO_IRX3, SOX9_PierreRobin) loses its hit: crediting every gene then "
            "crowds the published target, and the one-target sum was the better reader"
        ),
        "headline": (
            "target_derived 15/17 is NOT recomputed in place: swapping a layer of the shared benchmark is "
            "the owner's call. The counterfactual rate with the window layer in place of gene_input is "
            "reported beside it"
        ),
    },
}


# ---------------------------------------------------------------- the run (after the registration)

import argparse  # noqa: E402
import time  # noqa: E402
from collections import Counter  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any  # noqa: E402

from genomeos.results import RESULTS_DIR, load_result, save_result  # noqa: E402

REFERENCE = Path("data/reference")
STATED_AS_OF_THE_STORED_BENCHMARK = ("loci_stated_intervals",)


def coding_of(chrom: str) -> tuple[set[str], list[tuple[int, str]]]:
    """The chromosome's protein-coding symbols and their TSSs, from the per-chromosome GENCODE file the
    candidates reading uses (the element control below checks it names the sweep's coding heads)."""
    from genomeos.genome import Annotation

    ann = Annotation.from_gff3(REFERENCE / f"gencode_v50_{chrom}.gff3.gz", {chrom})
    coding = [g for g in ann.genes.values() if g.type == "protein_coding"]
    tss = sorted(
        ((g.locus.end - 1 if g.locus.strand.value == "-" else g.locus.start), g.symbol) for g in coding
    )
    return {g.symbol for g in coding}, tss


def head_control(chrom: str, rows: list[dict[str, Any]], responses: Any, coding: set[str]) -> Counter:
    """The shared control: the window's head at the bar is the table's `predicted` gene and its coding
    head is `predicted_coding`, element by element, wherever both sources speak."""
    c: Counter = Counter()
    for e in rows:
        if "predicted" not in e and "predicted_coding" not in e:
            continue
        got = responses.at_bar(chrom, e["id"], BAR) if e.get("id") else None
        if got is None:
            c["not_cached"] += 1
            continue
        c["compared"] += 1
        head = got[0][0] if got else None
        coding_head = next((g for g, _ in got if g in coding), None)
        c["head_disagreements"] += head != (e.get("predicted") or {}).get("gene")
        if "predicted_coding" not in e:
            # a row that never carried a coding field (the VISTA rows): the table is silent by
            # construction, which is a different thing from disagreeing with the cache
            c["no_coding_field"] += 1
            c["no_coding_field_but_the_cache_names_one"] += coding_head is not None
            continue
        c["coding_head_disagreements"] += coding_head != (e.get("predicted_coding") or {}).get("gene")
    return c


def motif(results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    from genomeos.attribution import motif_transfer as mt
    from genomeos.attribution.targets import ElementResponses

    responses = ElementResponses()
    tiers = mt.tier_blocks(mt.CHROMS, results_dir, responses=responses)
    stored = (load_result("motif_transfer", results_dir) or {}).get("real_unknown", {})
    cov = stored["measured"]["coverage"]
    got = {
        t: {
            "tested": sum(1 for b in blocks if b["tested_elements"]),
            "lead": sum(1 for b in blocks if b["moving_elements"]),
        }
        for t, blocks in tiers.items()
    }
    by_case = Counter(b["case"] for b in tiers["real_unknown"] if b["moving_elements"])
    want_case = {k: v["with_a_sweep_lead"] for k, v in stored["by_case"].items()}
    checks = {
        "real_unknown_tested": [
            got["real_unknown"]["tested"],
            cov["real_unknown"]["blocks_with_an_element_the_sweep_tested"],
        ],
        "real_unknown_lead": [got["real_unknown"]["lead"], cov["real_unknown"]["blocks_with_a_sweep_lead"]],
        "neutral_tested": [
            got["neutral"]["tested"],
            cov["neutral"]["blocks_with_an_element_the_sweep_tested"],
        ],
        "neutral_lead": [got["neutral"]["lead"], cov["neutral"]["blocks_with_a_sweep_lead"]],
        "leads_by_case": [dict(sorted(by_case.items())), dict(sorted(want_case.items()))],
    }
    window = {}
    for t, blocks in tiers.items():
        led = [b for b in blocks if b["moving_elements"]]
        movers = sum(b["moving_elements"] for b in led)
        window[t] = {
            "blocks_where_window_and_table_disagree_on_movers": sum(
                1 for b in blocks if b["moving_elements_window"] != b["moving_elements"]
            ),
            "elements_not_cached": sum(b["elements_not_cached"] for b in blocks),
            "head_genes_over_led_blocks": sum(len(b["head_genes"]) for b in led),
            "genes_at_bar_over_led_blocks": sum(len(b["genes_at_bar_window"]) for b in led),
            "moving_elements": movers,
            "blocks_whose_window_names_more_genes": sum(
                1 for b in led if len(b["genes_at_bar_window"]) > len(b["head_genes"])
            ),
        }
    return {
        "reproduced": {k: {"now": a, "stored": b, "equal": a == b} for k, (a, b) in checks.items()},
        "window": window,
        "caches_decompressed": len(responses._loaded),  # noqa: SLF001
    }


def candidates(results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    from genomeos.attribution import candidates as cm
    from genomeos.attribution.targets import ElementResponses

    stored = (load_result("syntax_candidates_genome_wide", results_dir) or {}).get("candidates", [])
    responses = ElementResponses()
    rows, control = [], Counter()
    fields = ("node", "deleted_elements_in_node", "node_targets_named", "deleted_elements_on_block")
    for chrom in sorted({c["chrom"] for c in stored}, key=lambda c: (len(c), c)):
        coding, tss = coding_of(chrom)

        class View:  # the three things node_context reads from a Chromosome, without its heavy loads
            pass

        ch = View()
        ch.chrom = chrom
        ch.coding_tss = tss
        ch.domains = sorted(
            (load_result(f"domains_{chrom}", results_dir) or {}).get("domains", []), key=lambda d: d["start"]
        )
        ch.domain_starts = [d["start"] for d in ch.domains]
        ch.node_at = lambda pos, ch=ch: cm.Chromosome.node_at(ch, pos)
        ch.deleted = cm.deletion_elements(chrom, results_dir)
        complete_then = chrom in ("chr21", "chr22")
        control.update(head_control(chrom, ch.deleted, responses, coding))
        for c in (x for x in stored if x["chrom"] == chrom):
            old = cm.node_context(ch, c["start"], c["end"])
            new = cm.node_context(ch, c["start"], c["end"], responses, coding)
            same_as_stored = all(old.get(f) == c.get(f) for f in fields)
            top_t = (old.get("node_targets_named") or [{}])[0].get("gene")
            top_w = (new.get("node_targets_named_window") or [{}])[0].get("gene")
            rows.append(
                {
                    "block": f"{chrom}:{c['start']}-{c['end']}",
                    "input_complete_on_2026_09_13": complete_then,
                    "old_reading_equals_stored": same_as_stored,
                    "elements_in_node_stored": c.get("deleted_elements_in_node"),
                    "elements_in_node_now": old.get("deleted_elements_in_node"),
                    "genes_named_table": new.get("node_genes_named_table"),
                    "genes_at_bar_window": new.get("node_genes_at_bar_window"),
                    "top_table": top_t,
                    "top_window": top_w,
                    "genes_below_table_count": new.get("genes_below_table_count"),
                    "elements_not_cached": new.get("elements_not_cached"),
                    "on_block": len(old.get("deleted_elements_on_block") or []),
                    "old_equals_new_on_old_fields": all(
                        old.get(f) == new.get(f) for f in fields if f != "deleted_elements_on_block"
                    ),
                }
            )
    then = [r for r in rows if r["input_complete_on_2026_09_13"]]
    later = [r for r in rows if not r["input_complete_on_2026_09_13"]]
    with_node_genes = [r for r in rows if r["genes_named_table"]]
    return {
        "candidates": len(rows),
        "control_stored_where_input_unchanged": {
            "candidates": len(then),
            "reproduced": sum(r["old_reading_equals_stored"] for r in then),
        },
        "stored_differs_by_input": {
            "candidates": len(later),
            "old_reading_equals_stored": sum(r["old_reading_equals_stored"] for r in later),
            "elements_in_node_stored": sum(r["elements_in_node_stored"] or 0 for r in later),
            "elements_in_node_now": sum(r["elements_in_node_now"] or 0 for r in later),
        },
        "head_control": dict(control),
        "falsifier_genes_below_table_count": sum(len(r["genes_below_table_count"] or []) for r in rows),
        "candidates_with_a_named_node_gene": len(with_node_genes),
        "genes_named_table_total": sum(r["genes_named_table"] or 0 for r in rows),
        "genes_at_bar_window_total": sum(r["genes_at_bar_window"] or 0 for r in rows),
        "top_node_gene_changed": sum(1 for r in with_node_genes if r["top_table"] != r["top_window"]),
        "candidates_with_elements_on_block": sum(1 for r in rows if r["on_block"]),
        "changed": [
            {
                k: r[k]
                for k in ("block", "top_table", "top_window", "genes_named_table", "genes_at_bar_window")
            }
            for r in with_node_genes
            if r["top_table"] != r["top_window"]
        ],
    }


def loci(results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    from genomeos.attribution.targets import ElementResponses
    from genomeos.benchmark import loci as lb

    stored = {x["locus"]: x for x in (load_result("loci_benchmark", results_dir) or {}).get("loci", [])}
    responses = ElementResponses()
    out, control = [], Counter()
    keep = (
        "elements_scored",
        "naming_a_coding_gene",
        "genes",
        "published_targets",
        "target",
        "rank_of_first_published_target",
    )
    panel = sorted(
        (e for e in lb.PANEL if e.locus in stored), key=lambda e: ((len(e.chrom), e.chrom), e.locus)
    )
    coding_by: dict[str, set[str]] = {}
    for e in panel:
        if e.chrom not in coding_by:
            coding_by.clear()
            coding_by[e.chrom] = coding_of(e.chrom)[0]
        coding = coding_by[e.chrom]

        class View:
            pass

        ch = View()
        ch.chrom = e.chrom
        got = lb.read_gene_input(ch, e, results_dir, responses, coding)
        # the stored benchmark was written at 7a1c751 (2026-09-17 09:53); two stated-interval files came
        # later (loci_third_intervals 2026-09-21, loci_fourth_intervals 2026-09-21) and one after it the
        # same morning (loci_candidate_intervals 10:02). The reader is checked on the input it had
        saved = lb.STATED_INTERVAL_RESULTS
        lb.STATED_INTERVAL_RESULTS = STATED_AS_OF_THE_STORED_BENCHMARK
        try:
            as_of = lb.read_gene_input(ch, e, results_dir)
        finally:
            lb.STATED_INTERVAL_RESULTS = saved
        rows = lb._deletion_rows(e.chrom, e.window[0], e.window[1], results_dir)  # noqa: SLF001
        control.update(head_control(e.chrom, rows, responses, coding))
        before = stored[e.locus]["readings"]["gene_input"]
        by_layer = stored[e.locus]["score"]["scored"]["target"]["by_layer"]
        others = [
            k for k, v in by_layer.items() if v["provenance"] == "derived" and v["hit"] and k != "gene_input"
        ]
        w = got["window_reading"]
        out.append(
            {
                "locus": e.locus,
                "reproduced": all(got.get(k) == before.get(k) for k in keep),
                "differing_fields": [k for k in keep if got.get(k) != before.get(k)],
                "reproduced_on_the_stored_input": all(as_of.get(k) == before.get(k) for k in keep),
                "hit_table": bool(got["genes"]) and got["genes"][0]["gene"] in e.targets,
                "hit_window": w["hit"],
                "rank_table": got["rank_of_first_published_target"],
                "rank_window": w["rank_of_first_published_target"],
                "top_table": got["genes"][0]["gene"] if got["genes"] else None,
                "top_window": w["target"],
                "genes_table": len({(r.get("predicted_coding") or {}).get("gene") for r in rows} - {None})
                if rows
                else 0,
                "genes_window": w["genes_credited"],
                "elements_not_cached": w["elements_not_cached"],
                "derived_hit_stored": bool(stored[e.locus]["score"]["scored"]["target"]["hit_derived"]),
                "derived_hit_with_window": bool(others) or w["hit"],
            }
        )
    misses = [r for r in out if not r["hit_table"]]
    med = lambda xs: sorted(xs)[len(xs) // 2] if xs else None  # noqa: E731
    return {
        "loci": len(out),
        "reproduced": sum(r["reproduced"] for r in out),
        "reproduced_on_the_stored_input": sum(r["reproduced_on_the_stored_input"] for r in out),
        "differing_by_input": {r["locus"]: r["differing_fields"] for r in out if not r["reproduced"]},
        "head_control": dict(control),
        "strict_hits_table": sum(r["hit_table"] for r in out),
        "strict_hits_window": sum(r["hit_window"] for r in out),
        "gained": [r["locus"] for r in out if r["hit_window"] and not r["hit_table"]],
        "lost": [r["locus"] for r in out if r["hit_table"] and not r["hit_window"]],
        "misses_median_rank_table": med([r["rank_table"] or 99 for r in misses]),
        "misses_median_rank_window": med([r["rank_window"] or 99 for r in misses]),
        "target_derived_stored": sum(r["derived_hit_stored"] for r in out),
        "target_derived_counterfactual_with_window": sum(r["derived_hit_with_window"] for r in out),
        "rows": out,
    }


def verdicts(res: dict[str, Any]) -> dict[str, Any]:
    """The registration read without interpretation, negatives first."""
    m, c, lo = res["motif_transfer"], res["candidates"], res["loci_gene_input"]
    sole = ("SHH_ZRS", "HERC2_OCA2", "FTO_IRX3", "SOX9_PierreRobin")
    return {
        "motif_transfer": (
            "reproduced: every block count equals the stored result"
            if all(v["equal"] for v in m["reproduced"].values())
            and all(w["blocks_where_window_and_table_disagree_on_movers"] == 0 for w in m["window"].values())
            else "FALSIFIED: a block count differs"
        ),
        "candidates": (
            "FALSIFIED: a gene's window count is below its table count"
            if c["falsifier_genes_below_table_count"]
            else "held: window counts >= table counts for every gene"
        )
        + (
            "; the stored control reproduced"
            if c["control_stored_where_input_unchanged"]["reproduced"]
            == c["control_stored_where_input_unchanged"]["candidates"]
            else "; the stored control did NOT reproduce"
        ),
        "loci_gene_input": (
            f"FALSIFIED: strict hits {lo['strict_hits_table']} -> {lo['strict_hits_window']} of {lo['loci']}"
            f" (lost {', '.join(lo['lost']) or 'none'}; sole-carrier losses "
            f"{', '.join(x for x in sole if x in lo['lost']) or 'none'})"
            if lo["strict_hits_window"] < 10 or any(x in lo["lost"] for x in sole)
            else "held: strict hits held or rose and no sole-carrier locus lost its hit"
        )
        + (
            # the registered control is the reader on the input the stored benchmark had (the stated
            # intervals as of 7a1c751); a locus that differs only on today's input is named, not hidden
            "; the stored fields reproduced on the stored input"
            if lo["reproduced_on_the_stored_input"] == lo["loci"]
            else "; stored fields did NOT reproduce on the stored input"
        )
        + (
            f"; on today's input {', '.join(sorted(lo['differing_by_input']))} differ(s)"
            " by later stated intervals"
            if lo["differing_by_input"]
            else ""
        ),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-save", action="store_true")
    ap.add_argument("--only", default="", help="comma-separated: loci,candidates,motif")
    args = ap.parse_args(argv)
    only = set(args.only.split(",")) if args.only else {"loci", "candidates", "motif"}
    t0 = time.time()
    res: dict[str, Any] = {"result": "onetarget_consumers", "pre_registration": PRE_REGISTRATION}
    if "loci" in only:
        res["loci_gene_input"] = loci()
        print("loci", {k: v for k, v in res["loci_gene_input"].items() if k != "rows"}, flush=True)
    if "candidates" in only:
        res["candidates"] = candidates()
        print("candidates", {k: v for k, v in res["candidates"].items() if k != "changed"}, flush=True)
    if "motif" in only:
        res["motif_transfer"] = motif()
        print("motif", res["motif_transfer"], flush=True)
    res["syntax_tiling"] = "not re-run: unrun by registration, and its tiles cost model requests"
    if only == {"loci", "candidates", "motif"}:
        res["verdicts"] = verdicts(res)
        print(res["verdicts"])
    res["alphagenome_requests"] = 0
    res["seconds"] = round(time.time() - t0, 1)
    if not args.no_save:
        print(f"saved {save_result(res['result'], res)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
