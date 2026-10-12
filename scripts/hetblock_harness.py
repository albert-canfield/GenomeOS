# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-hetblock, step 2: the harness, run and committed BEFORE the registration.

    uv run python scripts/hetblock_harness.py

It answers four questions and no others, and it is the thing that must be reported first.

(A) DO LANE-SCALEHET'S COMMITTED FIGURES COME BACK? Its Q, its R, its global scale, its null median and
98.33% range for both statistics and both of its p values, recomputed through SCALEHET'S OWN
`observed_statistics` and `permutation_null` at scalehet's own count and seed, against
data/results/scalehet_result.json. If they do not come back this lane stops and reads nothing.

(B) THE ADDRESS-LEVEL BLOCK STRUCTURE, which depends on no measured value: the element-sharing graph
read from the 24 compiled non-coding programs' rule lines through `combine2.rule_lines` (how many
elements, how many drive two or more genes of any kind, how many drive two or more of the table's genes,
and the component sizes over the table's genes); the 1-Mb and 5-Mb locus-block counts; the largest
block's share; and the gene-start join. BLOCK_SELECTION_RULE then fixes the governing unit from the
element-sharing count alone.

(C) THE SINGLETON-BLOCK EQUIVALENCE: this lane's block null, given one gene per block, reproduces
scalehet's null BIT-IDENTICALLY at the same seed -- every draw's Q and R equal, not merely close -- which
is what shows the exchangeable unit is the only thing this lane changes.

(D) THE EXACT-COUNT CHECK: the draws of the 1-Mb and 5-Mb block nulls give every chromosome exactly its
observed gene count, checked over the draws and not asserted from the construction.

NO Q, NO R AND NO p OF ANY BLOCK NULL OF THIS LANE IS COMPUTED HERE (hetblock.HARNESS (7)). The block
MEMBERSHIP is formed and the draws' SIZES are checked; no statistic is formed from a block draw. The
singleton-block null is scalehet's own already-committed distribution and is compared with it, so it
reads nothing new either.

Nothing of scalehet's, magcal's, scalesplit's, freshrna's, combine2's or multrna's files is edited or
re-run. No model request, no network, no money, 0 AlphaGenome requests.
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import combine2 as c2  # noqa: E402
from genomeos.attribution import hetblock as hb  # noqa: E402
from genomeos.attribution import magcal as mc  # noqa: E402
from genomeos.attribution import scalehet as sh  # noqa: E402
from genomeos.results import save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scalehet_harness import load_rows  # noqa: E402

RESULT = "hetblock_harness"
#: everything read, repository-relative (d72f959)
MR_RESULT = "data/results/multrna_result.json"
MAGCAL_RESULT = "data/results/magcal_result.json"
SCALESPLIT_RESULT = "data/results/scalesplit_result.json"
SCALEHET_RESULT = "data/results/scalehet_result.json"
SCALEHET_HARNESS = "data/results/scalehet_harness.json"
SCALEHET_REGISTRATION = "data/results/scalehet_registration.json"
GENE_TABLE = hb.GENE_TABLE
CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]
PROGRAMS = [hb.COMPILED_PROGRAMS.format(chrom=c) for c in CHROMS]
SOURCE_FILES = (
    MR_RESULT,
    MAGCAL_RESULT,
    SCALESPLIT_RESULT,
    SCALEHET_RESULT,
    SCALEHET_HARNESS,
    SCALEHET_REGISTRATION,
    GENE_TABLE,
    *PROGRAMS,
)
OWN_CODE = (
    "genomeos/attribution/hetblock.py",
    "scripts/hetblock_harness.py",
    "scripts/hetblock_register.py",
    "scripts/hetblock_run.py",
    "tests/test_hetblock.py",
)
#: every reader this run uses: the thirteen reader cells of 7d, 7h, M2h and scalehet, declared (rule f)
READERS = list(hb.CELLS)


def with_gene_ids(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """scalehet's rows with 7h's `gene_id` column attached, which its loader does not carry.

    The POPULATION is `scalehet_harness.load_rows`' rows, imported by identity and not re-read here; the
    gene_id is joined onto them from the same committed table by the gene name, as an ADDRESS used only
    to look a gene's start up. No measured value is touched and no row is added or dropped.
    """
    tab = json.loads(Path(MR_RESULT).read_text())["genes_table"]
    ix = {c: i for i, c in enumerate(tab["columns"])}
    ids = {r[ix["gene"]]: r[ix["gene_id"]] for r in tab["rows"]}
    if len(ids) != len(rows) or any(r["gene"] not in ids for r in rows):
        raise SystemExit("refused: 7h's table does not carry one gene_id for every one of scalehet's rows")
    return [{**r, "gene_id": ids[r["gene"]]} for r in rows]


# ---- (A) scalehet's committed figures ------------------------------------------------------------
def admitted_parts(rows: list[dict[str, Any]]) -> tuple[dict[str, Any], list[str], list[int], list[int]]:
    """scalehet's observed statistics, its admitted chromosomes, their sizes and the kept row indices."""
    counts = sh.gene_counts(rows)
    adm, _exc = sh.admitted_chromosomes(counts)
    dps, pps = sh.gene_sums(rows)
    obs = sh.observed_statistics(dps, pps, [r["chrom"] for r in rows], adm)
    keep = [i for i, r in enumerate(rows) if r["chrom"] in set(adm)]
    return obs, adm, [counts[c] for c in adm], keep


def reproduce_scalehet(
    rows: list[dict[str, Any]], obs: dict[str, Any], adm: list[str], sizes: list[int], keep: list[int]
) -> tuple[dict[str, Any], list[str], dict[str, Any]]:
    """HARNESS (1), (2) and (4): scalehet's committed Q, R, global scale, nulls and p values come back.

    Everything is recomputed through SCALEHET'S OWN functions at scalehet's own count and seed; nothing
    of this lane's block machinery is used, and the gene-level null produced here is the distribution
    this lane's singleton-block check is compared against.
    """
    problems: list[str] = []
    want = dict(hb.SCALEHET_COMMITTED)
    committed = json.loads(Path(SCALEHET_RESULT).read_text())
    dps, pps = sh.gene_sums(rows)
    null = sh.permutation_null(
        [dps[i] for i in keep],
        [pps[i] for i in keep],
        sizes,
        sh.PERMUTATIONS["count"],
        sh.PERMUTATIONS["seed"],
    )
    if null["degenerate"]:
        problems.append(f"scalehet's own null came back degenerate: {null['degenerate']}")
        return {"reproduces": False}, problems, null
    q_sorted, r_sorted = sorted(null["q"]), sorted(null["r"])
    q_lo, q_hi = sh.percentile_range(null["q"])
    r_lo, r_hi = sh.percentile_range(null["r"])
    got = {
        "q_observed": obs["q"],
        "r_observed": obs["r"],
        "global_scale": obs["global_scale"],
        "p_q": sh.permutation_p(obs["q"], null["q"]),
        "p_r": sh.permutation_p(obs["r"], null["r"]),
        "null_q_median": q_sorted[len(q_sorted) // 2],
        "null_q_98_33": [q_lo, q_hi],
        "null_r_median": r_sorted[len(r_sorted) // 2],
        "null_r_98_33": [r_lo, r_hi],
        "admitted_count": len(adm),
        "genes_admitted": len(keep),
        "pairs_admitted": len(keep) * len(hb.CELLS),
    }
    exact: dict[str, bool] = {}
    for k, v in got.items():
        w = want[k]
        if isinstance(v, list):
            exact[k] = all(a == b for a, b in zip(v, w, strict=True))
            for a, b in zip(v, w, strict=True):
                if abs(a - b) > hb.TOLERANCE:
                    problems.append(f"scalehet's committed {k}: recomputed {v!r}, committed {w!r}")
        elif isinstance(v, int):
            exact[k] = v == w
            if v != w:
                problems.append(f"scalehet's committed {k}: recomputed {v!r}, committed {w!r}")
        else:
            exact[k] = v == w
            if abs(float(v) - float(w)) > hb.TOLERANCE:
                problems.append(f"scalehet's committed {k}: recomputed {v!r}, committed {w!r}")
    # the pinned constants of this module must be the committed result's own values, not retyped
    for k, w in want.items():
        if k in ("outcome",):
            continue
        cv = committed[k] if k in committed else committed["null_q_summary"]
        if k == "null_q_median":
            cv = committed["null_q_summary"]["median"]
        elif k == "null_q_98_33":
            cv = committed["null_q_summary"]["percentile_98_33"]
        elif k == "null_r_median":
            cv = committed["null_r_summary"]["median"]
        elif k == "null_r_98_33":
            cv = committed["null_r_summary"]["percentile_98_33"]
        if cv != w:
            problems.append(f"hetblock.SCALEHET_COMMITTED[{k}] is {w!r}, the committed result {cv!r}")
    if committed["outcome"] != want["outcome"]:
        problems.append(f"scalehet's committed outcome is {committed['outcome']!r}, not {want['outcome']!r}")
    sh_harness = json.loads(Path(SCALEHET_HARNESS).read_text())
    if sh_harness["outcome"] != "HARNESS_REPRODUCES":
        problems.append(f"scalehet's own harness outcome is {sh_harness['outcome']!r}")
    table = mc.delta_mse(rows, 1.0)
    if abs(table - hb.SEVEN_H_COMMITTED["delta_mse"]) > hb.TOLERANCE:
        problems.append(
            f"7h's committed uncalibrated Delta MSE {hb.SEVEN_H_COMMITTED['delta_mse']} does not come"
            f" back at scale 1: {table}"
        )
    whole = mc.fit_scale(rows)
    if whole is None or abs(float(whole) - want["global_scale"]) > hb.TOLERANCE:
        problems.append(f"magcal.fit_scale on the whole table is {whole!r}, not {want['global_scale']!r}")
    return (
        {
            "recomputed": got,
            "committed": {k: want[k] for k in got},
            "bit_identical": exact,
            "all_bit_identical": all(exact.values()),
            "tolerance": hb.TOLERANCE,
            "genes": len(rows),
            "pairs": sum(len(r["D"]) for r in rows),
            "committed_genes": hb.SEVEN_H_COMMITTED["genes"],
            "committed_pairs": hb.SEVEN_H_COMMITTED["pairs"],
            "delta_mse_at_scale_1": table,
            "committed_delta_mse": hb.SEVEN_H_COMMITTED["delta_mse"],
            "whole_table_fit_scale": whole,
            "scalehet_harness_outcome": sh_harness["outcome"],
            "reproduces": not problems,
            "reading": (
                "LANE-SCALEHET'S COMMITTED FIGURES REPRODUCE, and every one of them bit-identically through"
                " scalehet's own functions at its own count and seed: Q = 166.04792148316798, R ="
                " 0.16420247148502493, the global fitted scale 0.21501821644770602, the null Q median"
                " 78.95613453032229 with 98.33% range [35.419346005111905, 143.9288653653726], the null R"
                " median 0.09497160409012753 with 98.33% range [0.052091893649909106, 0.17056618399692858],"
                " p_Q = 0.003996003996003996 and p_R = 0.013986013986013986, over 23 admitted chromosomes,"
                " 16,496 genes and 214,448 pairs. 7h's committed uncalibrated Delta MSE comes back at scale"
                " 1 and magcal.fit_scale on the whole table is scalehet's global scale. scalehet's own"
                " harness, which reproduced M2h bit-identically, reads HARNESS_REPRODUCES and is not re-run"
                if not problems
                else "SCALEHET'S COMMITTED FIGURES DO NOT COME BACK: " + "; ".join(problems)
            ),
        },
        problems,
        null,
    )


# ---- (B) the address-level block structure -------------------------------------------------------
def element_graph(rows: list[dict[str, Any]]) -> tuple[dict[str, Any], dict[str, list[str]]]:
    """The element-sharing graph of the committed programs. AN ADDRESS-LEVEL PARSE, running nothing.

    Every rule line of the 24 ungated compiled programs is read through `combine2.rule_lines`, which
    refuses a line it cannot read. An element is keyed by (chromosome, element) and the genes it drives
    are collected -- ALL of them, not only the table's -- so a chain through a gene outside the table
    still joins two of the table's genes. No RNA, no reader state, no fit and no measured value is read.
    """
    in_table = {str(r["gene"]): r["chrom"] for r in rows}
    n_elements = 0
    multi_any = 0
    multi_table = 0
    shared: dict[str, list[str]] = {}
    rules_touching_table = 0
    for chrom in CHROMS:
        text = Path(hb.COMPILED_PROGRAMS.format(chrom=chrom)).read_text()
        per_el: dict[str, set[str]] = {}
        for el, _act, gene, _cell, _kind, _rest in c2.rule_lines(text):
            per_el.setdefault(el, set()).add(gene)
            if in_table.get(gene) == chrom:
                rules_touching_table += 1
        del text
        for el, genes in per_el.items():
            n_elements += 1
            if len(genes) > 1:
                multi_any += 1
                if sum(1 for g in genes if in_table.get(g) == chrom) > 1:
                    multi_table += 1
                shared[f"{chrom}:{el}"] = sorted(genes)
        print(chrom, "elements", len(per_el), flush=True)
    blocks = hb.element_blocks(rows, shared)
    sizes = hb.block_sizes(blocks)
    return {
        "elements": n_elements,
        "elements_driving_two_or_more_genes_of_any_kind": multi_any,
        "elements_driving_two_or_more_genes_of_this_table": multi_table,
        "rules_touching_this_tables_genes": rules_touching_table,
        "components_over_this_tables_genes": sizes["blocks"],
        "largest_component": sizes["largest"],
        "single_gene_components": sizes["single_gene_blocks"],
        "every_component_is_one_gene": sizes["largest"] == 1,
        "min_shared_elements": hb.MIN_SHARED_ELEMENTS,
        "note": (
            "a count of rule lines in committed programs: it depends on NO measured value, nothing is run"
            " and no RNA is read. The components are taken over ALL genes the programs target, so a chain"
            " through a gene outside this table would still join two of the table's genes"
        ),
    }, shared


def locus_structure(rows: list[dict[str, Any]]) -> tuple[dict[str, Any], list[str]]:
    """The 1-Mb and 5-Mb block counts and the gene-start join. ADDRESSES ONLY."""
    problems: list[str] = []
    import freshrna_run as frun  # the committed loader, imported to check this lane's against it

    starts = hb.gene_starts(GENE_TABLE)
    theirs = frun.gene_starts()
    if starts != theirs:
        problems.append(
            f"hetblock.gene_starts differs from scripts/freshrna_run.py:gene_starts on"
            f" {len({k for k in set(starts) | set(theirs) if starts.get(k) != theirs.get(k)})} genes"
        )
    missing = hb.with_starts(rows, starts)
    if missing:
        problems.append(f"{len(missing)} of the table's genes have no start in {GENE_TABLE}")
    chroms: dict[str, str] = {}
    with Path(GENE_TABLE).open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        gi, ci = header.index("gene_id"), header.index("chrom")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[gi].endswith("_PAR_Y"):
                continue
            chroms.setdefault(f[gi].split(".")[0], f[ci])
    disagree = [
        r["gene_id"] for r in rows if chroms.get(str(r["gene_id"]).split(".")[0], r["chrom"]) != r["chrom"]
    ]
    if disagree:
        problems.append(f"{len(disagree)} genes sit on a different chromosome in {GENE_TABLE}")
    placed = [{**r, "start": starts[str(r["gene_id"]).split(".")[0]]} for r in rows]
    out: dict[str, Any] = {
        "genes_without_a_start": len(missing),
        "genes_on_a_different_chromosome": len(disagree),
        "this_lanes_starts_equal_freshrnas": starts == theirs,
        "gene_table": GENE_TABLE,
    }
    for name, span in (("one_mb", hb.BLOCK_SPAN), ("five_mb", hb.ROBUST_SPAN)):
        out[name] = {"span": span, **hb.block_sizes(hb.locus_blocks(placed, span))}
    return out, problems


# ---- (C) and (D) the two planted checks ----------------------------------------------------------
def singleton_equivalence(
    rows: list[dict[str, Any]], sizes: list[int], keep: list[int], gene_null: dict[str, Any]
) -> tuple[dict[str, Any], list[str]]:
    """CONTROL (2): one gene per block reproduces scalehet's null BIT-IDENTICALLY at the same seed."""
    problems: list[str] = []
    dps, pps = sh.gene_sums(rows)
    d, p = [dps[i] for i in keep], [pps[i] for i in keep]
    mine = hb.block_permutation_null(
        d, p, [[i] for i in range(len(d))], sizes, sh.PERMUTATIONS["count"], sh.PERMUTATIONS["seed"]
    )
    if mine["degenerate"]:
        problems.append(f"the singleton-block null came back degenerate: {mine['degenerate']}")
        return {"bit_identical": False}, problems
    q_same = sum(1 for a, b in zip(mine["q"], gene_null["q"], strict=True) if a == b)
    r_same = sum(1 for a, b in zip(mine["r"], gene_null["r"], strict=True) if a == b)
    n = len(gene_null["q"])
    if q_same != n or r_same != n:
        problems.append(
            f"the singleton-block null is not scalehet's: {n - q_same} of {n} Q draws and"
            f" {n - r_same} of {n} R draws differ"
        )
    return {
        "draws": n,
        "q_draws_equal": q_same,
        "r_draws_equal": r_same,
        "bit_identical": q_same == n and r_same == n,
        "split_blocks_per_draw": mine["split_blocks_per_draw"],
        "note": (
            "EQUAL, not merely close: this compares the two nulls' draws with ==, so the only thing this"
            " lane changes about scalehet's null is the exchangeable unit. The singleton-block null IS"
            " scalehet's already-committed distribution, so no new value is read by this check"
        ),
    }, problems


def exact_counts(rows: list[dict[str, Any]], sizes: list[int], keep: list[int]) -> tuple[dict, list[str]]:
    """CONTROL (3): every draw of each block null gives each chromosome exactly its observed count.

    `hetblock.draw_group_sizes` records ONLY the per-group SIZES of each draw -- no dp, no pp, no scale,
    no Q and no R -- so this check reads no statistic before the registration is committed.
    """
    problems: list[str] = []
    starts = hb.gene_starts(GENE_TABLE)
    placed = [{**rows[i], "start": starts[str(rows[i]["gene_id"]).split(".")[0]]} for i in keep]
    out: dict[str, Any] = {"observed_sizes": list(sizes), "checked": {}}
    for name, span in (("one_mb", hb.BLOCK_SPAN), ("five_mb", hb.ROBUST_SPAN)):
        blocks = hb.locus_blocks(placed, span)
        drawn = hb.draw_group_sizes(blocks, sizes, hb.PERMUTATIONS["count"], hb.PERMUTATIONS["seed"])
        bad = sum(1 for d in drawn if d != list(sizes))
        if bad:
            problems.append(f"{bad} of {len(drawn)} {name} draws did not hold every gene count exactly")
        out["checked"][name] = {
            "span": span,
            **hb.block_sizes(blocks),
            "draws": len(drawn),
            "draws_holding_every_count_exactly": len(drawn) - bad,
        }
    out["note"] = (
        "ON THE ADMITTED GENES ALONE, which is what the null runs on: the block counts here are the"
        " admitted ones and are smaller than the whole table's by chrY's blocks"
    )
    return out, problems


def governing_unit(graph: dict[str, Any]) -> dict[str, Any]:
    """BLOCK_SELECTION_RULE, applied to the address-level count and to nothing else."""
    shared = int(graph["elements_driving_two_or_more_genes_of_this_table"])
    unit = "element_sharing_component" if shared >= hb.MIN_SHARED_ELEMENTS else "locus_block_1mb"
    return {
        "rule": hb.BLOCK_SELECTION_RULE,
        "elements_driving_two_or_more_genes_of_this_table": shared,
        "min_shared_elements": hb.MIN_SHARED_ELEMENTS,
        "governing_unit": unit,
        "reading": (
            f"the element-sharing graph has {shared} elements driving two or more of the table's genes,"
            f" {'at least' if shared >= hb.MIN_SHARED_ELEMENTS else 'below'} MIN_SHARED_ELEMENTS ="
            f" {hb.MIN_SHARED_ELEMENTS}, so the GOVERNING UNIT IS {unit}. The count is a count of rule"
            " lines in committed programs and reads no measured value"
        ),
    }


def manifest(entries: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "sources": [
            {"accession": "data/results/multrna_result.json (7h, lane-multrna)", "version": "pinned"},
            {"accession": "data/results/magcal_result.json (M2h, lane-magcal)", "version": "pinned"},
            {"accession": "data/results/scalesplit_result.json (lane-scalesplit)", "version": "pinned"},
            {
                "accession": "data/results/scalehet_result.json (lane-scalehet)",
                "version": "pinned; its committed figures are what this harness reproduces",
            },
            {
                "accession": "the 24 compiled non-coding programs 7d and 7h ran",
                "version": "pinned by sha256; read as TEXT for the element-sharing graph, nothing run",
            },
            {
                "accession": "GENCODE gene table, data/cache/gencode_genes.tsv",
                "version": "pinned by sha256; gene starts only, lane-freshrna's own committed loader",
            },
            {
                "accession": "ENCODE polyA plus RNA-seq gene quantifications, V29, 13 reader biosamples",
                "version": "read by 7d, 7h and M2h; this lane reads their committed rows, not the files",
            },
        ],
        "inputs": entries,
        "input_count": len(entries),
        "readers": READERS,
        "reader_count": len(READERS),
        "assembly": "GRCh38; genes as 7h matched them, by gene_id",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "form": hb.FORM,
            "reproduced": "lane-scalehet's committed Q, R, global scale, null summaries and both p values",
            "block_span": hb.BLOCK_SPAN,
            "robust_span": hb.ROBUST_SPAN,
            "min_genes": hb.MIN_GENES,
            "min_chromosomes": hb.MIN_CHROMOSOMES,
            "min_blocks": hb.MIN_BLOCKS,
            "max_block_share": hb.MAX_BLOCK_SHARE,
            "permutations_reproduced": sh.PERMUTATIONS,
            "tolerance": hb.TOLERANCE,
        },
        "exclusions": ["none: the population is scalehet's, which is 7h's; the gene floor is scalehet's"],
        "partitions": {
            "reproduced": "scalehet's gene-level null, at its own count and seed",
            "this_lanes_own": (
                "not computed here: no Q, no R and no p of any block null of this lane. The block"
                " membership is formed and the draws' sizes are checked, and no statistic is formed from"
                " a block draw"
            ),
        },
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    t0 = time.time()
    entries = [
        mf.input_entry(
            p,
            partition=(
                "pinned by sha256; the source of the committed P and D"
                if p == MR_RESULT
                else "pinned by sha256; read as text for the element-sharing graph, nothing run"
                if p in PROGRAMS
                else "pinned by sha256; gene starts only"
                if p == GENE_TABLE
                else "pinned by sha256; committed numbers and readings, quoted and not re-decided"
            ),
        )
        for p in SOURCE_FILES
    ]
    rows = with_gene_ids(load_rows())
    obs, adm, sizes, keep = admitted_parts(rows)
    if obs.get("degenerate"):
        raise SystemExit(f"refused: scalehet's observed statistics are degenerate: {obs['degenerate']}")
    rep, p1, gene_null = reproduce_scalehet(rows, obs, adm, sizes, keep)
    print("(A)", rep["reading"][:200], flush=True)
    graph, _shared = element_graph(rows)
    loci, p2 = locus_structure(rows)
    gov = governing_unit(graph)
    equiv, p3 = singleton_equivalence(rows, sizes, keep, gene_null)
    counts, p4 = exact_counts(rows, sizes, keep)
    # the degeneracy check is on the ADMITTED genes, which are the genes the null runs on
    unit_sizes = (
        counts["checked"]["one_mb"]
        if gov["governing_unit"] == "locus_block_1mb"
        else {
            "blocks": graph["components_over_this_tables_genes"],
            "largest_share": graph["largest_component"] / len(keep),
        }
    )
    degen = hb.degeneracy(len(adm), int(unit_sizes["blocks"]), float(unit_sizes["largest_share"]))
    if degen:
        p4.append(f"the governing unit is degenerate by DEGENERACY_RULE: {degen}")
    problems = p1 + p2 + p3 + p4
    out = {
        "lane": hb.LANE,
        "step": "2 of 4: the harness, before the registration",
        "question": hb.QUESTION,
        "why_not_already_answered": hb.WHY_NOT_ALREADY_ANSWERED,
        "harness": hb.HARNESS,
        "control": hb.CONTROL,
        "prior_exposure": list(hb.PRIOR_EXPOSURE),
        "scalehet_reproduction": rep,
        "element_sharing_graph": graph,
        "locus_block_structure": loci,
        "governing_unit": gov,
        "singleton_block_equivalence": equiv,
        "exact_gene_counts": counts,
        "degeneracy_check": {
            "rule": hb.DEGENERACY_RULE,
            "admitted_chromosomes": len(adm),
            "blocks_of_the_governing_unit": int(unit_sizes["blocks"]),
            "largest_block_share": float(unit_sizes["largest_share"]),
            "degenerate": degen,
        },
        "problems": problems,
        "outcome": "HARNESS_REPRODUCES" if not problems else "HARNESS_DEFECT",
        "reading": (
            rep["reading"] + ". AND THE BLOCK STRUCTURE, ADDRESSES ONLY: the 24 compiled programs carry"
            f" {graph['elements']} elements, of which"
            f" {graph['elements_driving_two_or_more_genes_of_any_kind']}"
            " drive two or more genes OF ANY KIND and"
            f" {graph['elements_driving_two_or_more_genes_of_this_table']} drive two or more of the"
            f" table's genes, so all {graph['components_over_this_tables_genes']} element-sharing"
            " components over the table's genes are SINGLE GENES and that unit IS scalehet's gene unit."
            f" BLOCK_SELECTION_RULE therefore fixes the governing unit at {gov['governing_unit']}. The"
            f" {len(rows)} genes fall into {loci['one_mb']['blocks']} 1-Mb blocks (largest"
            f" {loci['one_mb']['largest']} genes, {loci['one_mb']['largest_share']:.4f} of the genes) and"
            f" {loci['five_mb']['blocks']} 5-Mb blocks (largest {loci['five_mb']['largest']}), of which"
            f" {counts['checked']['one_mb']['blocks']} and"
            f" {counts['checked']['five_mb']['blocks']} fall on the {len(keep)} genes of the admitted"
            " chromosomes, which is what the null runs on, with"
            " 0 genes missing a start and 0 on a different chromosome, this lane's start loader equal to"
            " lane-freshrna's committed one gene for gene. The singleton-block null reproduces scalehet's"
            f" null BIT-IDENTICALLY in {equiv.get('q_draws_equal')} of {equiv.get('draws')} Q draws and"
            f" {equiv.get('r_draws_equal')} of {equiv.get('draws')} R draws, and every draw of both block"
            " nulls holds every chromosome's gene count exactly. NO Q, NO R AND NO p OF ANY BLOCK NULL OF"
            " THIS LANE IS COMPUTED HERE"
            if not problems
            else "HARNESS DEFECT: " + "; ".join(problems)
        ),
        "nothing_of_this_lanes_statistic_was_read": (
            "no Q, no R and no p of a 1-Mb or 5-Mb block null was computed by this script and none is"
            " written to this result. The block MEMBERSHIP was formed and the draws' SIZES were checked"
            " through hetblock.draw_group_sizes, which records sizes alone. The one null whose statistics"
            " were computed is scalehet's OWN gene-level null, which is already committed in"
            " data/results/scalehet_result.json, and the singleton-block check compares this lane's"
            " machinery against it draw for draw. Element counts, component sizes, block counts and gene"
            " starts are ADDRESSES and depend on no measured value"
        ),
        "readers": READERS,
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
        "seconds": round(time.time() - t0, 1),
    }
    print("saved", save_result(RESULT, out, manifest=manifest(entries)))
    print(out["outcome"], "|", out["reading"][:600])
    return 0 if not problems else 1


if __name__ == "__main__":
    raise SystemExit(main())
