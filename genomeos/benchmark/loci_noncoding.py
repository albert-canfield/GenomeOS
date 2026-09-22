# SPDX-License-Identifier: AGPL-3.0-or-later
"""The elements the fourth frame set aside because their published target is not protein coding.

`loci_fourth.assess` applies the drawing rule to one perturbed element, step by step, and the FIRST
step it applies is this one::

    non_coding = [t for t in row["targets"] if t not in coding]
    if non_coding:
        verdict["reason"] = "a regulated target is not protein coding in GENCODE"

Fourteen elements of the held-out arm stopped there, and `loci_fourth` counted them apart rather than
folding them into its geometry question. Section 22 says why: the deletion layer ranks
`predicted_coding` first, so those elements test the defect section 18 found at the H19 ICR - a
published NON-coding target hidden behind a coding read-through - and not the geometry the fourth
frame was measuring. Section 21 had said in as many words that curation could not produce a published
non-coding target; the held-out arm produced these by rule, without anybody choosing them.

This module takes that set as its frame. It adds no curation: membership is `loci_fourth`'s own
branch, re-derived from the same file and the same GENCODE release, and `dropped_for_non_coding`
is checked against the committed verdicts in `data/results/loci_fourth.json` by the tests.

**What the frame turns out to be, and it is settled before any scoring.** `coding` in the branch above
is a set of SYMBOLS, and the ENCODE benchmark's `measuredGeneSymbol` column is as old as the screen
that filled it. The benchmark file also carries `measuredGeneEnsemblId`, which does not go stale, and
joining that column to GENCODE's own table says what each of the fourteen targets actually is:

===============  =================  ============================  ==========
published symbol  Ensembl id        GENCODE now                   elements
===============  =================  ============================  ==========
SSFA2            ENSG00000138434    ITPRID2, protein_coding       10
SARS             ENSG00000031698    SARS1, protein_coding          1
WDR61            ENSG00000140395    SKIC8, protein_coding          1
LINC00885        ENSG00000224652    LINC00885, lncRNA              1
CCDC26           ENSG00000229140    CCDC26, lncRNA                 1
===============  =================  ============================  ==========

So **twelve of the fourteen are stale gene symbols, not non-coding targets**, and the frame this
module scores is **n = 2**. That is the first registered number here, and it is arithmetic over two
files rather than anything a run discovers. `corrected_rule` states the other half of the arithmetic:
under the current symbol, eight of the twelve pass every step of `loci_fourth`'s rule and would have
been DRAWN into the fourth frame, three are rejected one step later because the nearest coding TSS is
the target after all, and one falls to the keep-out. That is reported to the fourth frame's owner as a
finding about its draw; nothing here re-cuts a frame that has been scored.

**Whether the layers can express a non-coding answer at all**, which is the question that decides
whether a rate here is a measurement or a division:

* ``gene_input`` (the summed window) reads ``predicted_coding`` and nothing else (`loci.py` 1569).
  It can never name a non-coding gene. **Zero by construction.**
* ``node`` (the heuristic) is the nearest coding TSS inside the CTCF node. **Zero by construction.**
* ``coding`` is a protein-change layer and is not built for an ``activates`` locus at all.
* ``deletion`` walks ``("predicted_coding", "predicted")`` and **breaks at the first key that has a
  gene** (`loci.py` 1103-1116). ``predicted`` names the gene that moves most whatever its type, so the
  layer CAN return a non-coding gene - but only where the element has no coding prediction at all.
  Where it has one, the non-coding answer is discarded unread. That is the H19 defect exactly, and it
  is a conditional zero rather than a flat one.
* ``eqtl`` maps GTEx gene ids through `loci.symbols`, which is the whole GENCODE table, lncRNAs
  included, so it can express a non-coding answer. It is the only derived layer here that can do so
  freely. This frame does not distil GTEx for it: a pass over the 49-tissue archive for two windows is
  not proportionate, and `loci.stream_gtex` clears the directory whose manifest changed, so sharing
  `data/knowledge/loci_fourth` would delete a peer frame's cache. The layer is therefore **unasked,
  reported as pending and never as a miss** - registered here, not discovered later.

So two of the four derived paths are zero before a locus is read, one is unasked, and exactly one can
be right. A derived rate over this frame is very nearly arithmetic, and the registration says so
rather than letting a low number read as a measurement of the model.

**Requests: zero, and that is registered as a hard budget.** Both elements are already covered by the
finished all-element sweep, so every reading below is a lookup into data this project has already
bought. `REQUEST_BUDGET` is 0 and `plan` refuses to proceed if the frame would cost one.

**Disclosure.** Establishing that ``deletion``'s fallback path is reachable at all meant reading the
sweep row over the chr8 element, and that row was read BEFORE this registration was written. What the
deletion layer does at CCDC26 was therefore known in advance and is **not** a prediction this run
confirms; it is stated in `PREREGISTRATION["known_before_the_registration"]` and reported as a
disclosure in the section. The chr3 element's row was read in the same call. Nothing else about either
locus - no other layer, no control, no aggregate - was read before this file was committed.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from genomeos.attribution import crispri
from genomeos.benchmark import loci, loci_fourth
from genomeos.benchmark.loci import Expect
from genomeos.results import RESULTS_DIR, save_result

NAME = "loci_noncoding"
#: the file the fourteen come from, through `loci_fourth`, so there is one source of truth.
SOURCE = loci_fourth.SOURCE
#: GENCODE's own gene table, which carries `gene_type` and the current symbol per Ensembl id.
GENE_TABLE = Path("data/cache/gencode_genes.tsv")
#: this frame spends nothing. Both elements are in the finished sweep; see the module docstring.
REQUEST_BUDGET = 0
#: half the locus window, as in every other frame.
WINDOW_HALF = loci_fourth.WINDOW_HALF
#: this frame's own GTEx directory. Left EMPTY on purpose: see the docstring. It is named rather than
#: reused so that no run here can clear `data/knowledge/loci_fourth` by writing a different manifest.
GTEX_DIR = Path("data/knowledge/loci_noncoding")

PREREGISTRATION: dict[str, Any] = {
    "written": "2026-09-21, before any reading was scored and with a budget of zero requests",
    "the_frame": (
        "every perturbed element of the ENCODE CRISPR benchmark's held-out arm that"
        " `loci_fourth.assess` stopped at its FIRST branch - a regulated target absent from the"
        " chromosome's protein-coding symbols in GENCODE. Membership is that branch's own output,"
        " re-derived from the same file and checked against the committed verdicts in"
        " data/results/loci_fourth.json. No locus is curated, added, dropped or renamed by hand"
    ),
    "n_as_the_branch_reports_it": 14,
    "n_after_the_ensembl_join": 2,
    "why_n_falls": (
        "the branch compares SYMBOLS, and the benchmark's `measuredGeneSymbol` column is as old as"
        " the screen that filled it. Joining the file's own `measuredGeneEnsemblId` to GENCODE's"
        " table shows twelve of the fourteen are protein-coding genes under a current symbol:"
        " SSFA2 is ITPRID2 (10 elements), SARS is SARS1 (1), WDR61 is SKIC8 (1). Two are genuinely"
        " non-coding: LINC00885 (chr3, alongside the coding co-target TFRC) and CCDC26 (chr8, sole"
        " target). This is arithmetic over two files, settled before any scoring, not a result"
    ),
    "request_budget": 0,
    "what_plan_costs": (
        "zero. Both elements are already covered by the finished all-element sweep, so every"
        " reading is a lookup into data already bought. `plan` raises rather than spend one, and a"
        " run that needed a request would be abandoned and reported as unaffordable, not funded"
    ),
    "predictions": {
        "reach_fatalities": (
            "0 of 2. This is not a hope: CCDC26's element at chr8:129,581,954 lies INSIDE the"
            " CCDC26 gene body, and the chr3 element is about 1 kb from LINC00885's 3' end, so both"
            " published targets are inside the scorer's 1,048,576 bp input by inspection of the"
            " annotation. If `loci.read_reach` says otherwise, the reach reader is wrong and that"
            " is the finding"
        ),
        "gene_input": "0 of 2, BY CONSTRUCTION: the layer reads `predicted_coding` only",
        "node": "0 of 2, BY CONSTRUCTION: the layer is the nearest coding TSS inside the CTCF node",
        "eqtl": "unasked at both, reported as pending. Not a miss and never counted as one",
        "deletion": (
            "the only derived layer that could be right here. It can return a non-coding gene only"
            " where `predicted_coding` is empty, because it breaks at the first key with a gene."
            " See `known_before_the_registration`: the answer at chr8 was already read"
        ),
        "the_derived_rate": (
            "expected 0 of 2, and it is to be read as ARITHMETIC, not as a measurement of the"
            " model. Two of the four derived paths are zero before a locus is read, one is unasked,"
            " and the fourth discards a non-coding answer whenever a coding one exists. A rate"
            " built on that says what the reader is allowed to say, not what the model knows"
        ),
    },
    "the_right_baseline": (
        "NOT the nearest coding TSS. That rule cannot name a non-coding gene, so it scores 0 of 2"
        " here for the same reason `gene_input` does, and quoting a beaten baseline that is zero by"
        " construction would be the same error twice. The baseline this frame is read against is"
        " THE MODEL'S OWN ANY-GENE PREDICTION - the `predicted` field the sweep already computed and"
        " `read_deletion` declines to look at. It is the right baseline because it is what the"
        " pipeline would answer if the one line that ranks coding genes first were removed, so the"
        " gap between it and the derived rate measures the READER's defect and nothing about the"
        " model's biology. `counterfactual` computes it, and it changes no shipped code"
    ),
    "outcomes": {
        "the_reader_is_the_limit": (
            "the derived rate is 0 of 2 while the any-gene baseline names a published target at one"
            " or both. Then the H19 defect is confirmed with a positive - the model had the answer"
            " and the reader threw it away - and the fix is one line in `loci.read_deletion`,"
            " named but not made here, because changing a shared reader is its owner's call"
        ),
        "the_model_is_the_limit": (
            "the derived rate is 0 of 2 AND the any-gene baseline is also 0 of 2. Then the defect"
            " is real but costs nothing at these two loci, the H19 finding stands untested for"
            " want of a case, and this frame has failed to test it"
        ),
        "failure": (
            "this frame FAILS if it reports a derived rate as though it measured the model. It also"
            " fails if n is quoted as 14 anywhere, if the eqtl pending is counted as a miss, if a"
            " request is spent, or if a locus is added, dropped or renamed to make a number move."
            " A frame of two is a frame of two and is reported as one"
        ),
    },
    "known_before_the_registration": (
        "the sweep rows over BOTH elements were read while establishing that `read_deletion`'s"
        " fallback to `predicted` is reachable at all, which is the arithmetic this registration"
        " turns on. At chr8 the row reads `predicted` = CCDC26 at -1.691 and `predicted_coding` ="
        " GSDMC at -0.1433; at chr3 both keys read ZDHHC19. The deletion-layer outcome at these two"
        " loci was therefore KNOWN when this was written and is not a prediction the run confirms."
        " It is disclosed here and in the section rather than presented as foresight"
    ),
    "what_this_frame_cannot_answer": (
        "whether the pipeline can find non-coding targets in general. n is 2, both are lncRNAs,"
        " both come from one screen, and one of them has a coding co-target. Nothing here"
        " generalises, and the section says so before it says anything else"
    ),
}


# ------------------------------------------------------------------ the frame, free, no request
def gene_table(path: Path = GENE_TABLE) -> dict[str, dict[str, str]]:
    """GENCODE's gene table keyed by unversioned Ensembl id: current symbol and `gene_type`."""
    out: dict[str, dict[str, str]] = {}
    with path.open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            out[r["gene_id"].split(".")[0]] = {
                "symbol": r["symbol"],
                "gene_type": r["gene_type"],
                "chrom": r["chrom"],
            }
    return out


def coding_symbols(path: Path = GENE_TABLE) -> dict[str, set[str]]:
    """Per chromosome, the protein-coding symbols - the same set `loci_fourth.assess` tests against.

    `loci_fourth` builds it from the per-chromosome GFF3 through `loci.Chromosome.coding`. This reads
    the cached genome-wide table instead, so the fourteen can be re-derived without opening a FASTA
    for every chromosome. The tests assert the two agree by comparing the result with the committed
    verdicts of the fourth frame's own draw.
    """
    out: dict[str, set[str]] = {}
    with path.open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["gene_type"] == "protein_coding":
                out.setdefault(r["chrom"], set()).add(r["symbol"])
    return out


def ensembl_ids(source: str = SOURCE, knowledge: Path = crispri.KNOWLEDGE) -> dict[str, str]:
    """`measuredGeneSymbol` to `measuredGeneEnsemblId`, read from the benchmark file's own columns.

    This is what makes the stale-symbol finding a join rather than a curation: no alias table is
    written here, and no symbol is corrected by anybody's memory.
    """
    import gzip

    out: dict[str, str] = {}
    with gzip.open(knowledge / source, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            gid = (r.get("measuredGeneEnsemblId") or "").strip()
            if gid:
                out.setdefault(r["measuredGeneSymbol"], gid.split(".")[0])
    return out


def dropped_for_non_coding(source: str = SOURCE, path: Path = GENE_TABLE) -> list[dict[str, Any]]:
    """The elements `loci_fourth.assess` stops at its first branch. Free: the file and GENCODE."""
    coding = coding_symbols(path)
    rows = loci_fourth.group_elements(loci_fourth.heldout(source))
    out = []
    for r in rows:
        non_coding = [t for t in r["targets"] if t not in coding.get(r["chrom"], set())]
        if non_coding:
            out.append({**r, "non_coding_targets": non_coding})
    return out


def resolve(
    rows: list[dict[str, Any]], source: str = SOURCE, path: Path = GENE_TABLE
) -> list[dict[str, Any]]:
    """Say what each published target actually is, by Ensembl id rather than by symbol."""
    ids = ensembl_ids(source)
    table = gene_table(path)
    out = []
    for r in rows:
        resolved = []
        for t in r["non_coding_targets"]:
            gid = ids.get(t)
            g = table.get(gid or "", {})
            resolved.append(
                {
                    "published_symbol": t,
                    "ensembl_id": gid,
                    "current_symbol": g.get("symbol"),
                    "gene_type": g.get("gene_type"),
                    "stale_symbol": bool(g.get("gene_type") == "protein_coding"),
                }
            )
        out.append({**r, "resolved": resolved})
    return out


def classify(rows: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    """Split the branch's output into stale symbols and targets that really are non-coding."""
    non_coding, stale = [], []
    for r in rows:
        (non_coding if any(not x["stale_symbol"] for x in r["resolved"]) else stale).append(r)
    return {"non_coding": non_coding, "stale_symbol": stale}


def corrected_rule(stale: list[dict[str, Any]]) -> dict[str, Any]:
    """What `loci_fourth`'s rule does to a stale-symbol element once the symbol is current.

    Reported to the fourth frame's owner. Nothing here re-cuts a frame that has already been scored:
    the fourth frame's n, denominator and headline are its owner's to change or to leave alone.
    """
    taken = [(e.chrom, (e.element[0] + e.element[1]) // 2) for e in loci_fourth.earlier_frames()]
    per: list[dict[str, Any]] = []
    for r in stale:
        current = {x["current_symbol"] or x["published_symbol"] for x in r["resolved"]}
        current |= {t for t in r["targets"] if t not in {x["published_symbol"] for x in r["resolved"]}}
        mid = r["mid"]
        nearest = r.get("nearest_coding")
        if nearest in current:
            verdict = "the nearest coding TSS IS a published target - the shortcut is right here"
        elif [c for c, m in taken if c == r["chrom"] and abs(m - mid) < loci_fourth.KEEP_OUT]:
            verdict = f"within {loci_fourth.KEEP_OUT} bp of a locus in an earlier frame"
        else:
            verdict = "would have been DRAWN into the fourth frame"
        per.append(
            {
                "chrom": r["chrom"],
                "element": [r["start"], r["end"]],
                "cell": r["cell"],
                "published": sorted(r["targets"]),
                "current": sorted(current),
                "nearest_coding": nearest,
                "verdict": verdict,
            }
        )
    tally: dict[str, int] = {}
    for p in per:
        tally[p["verdict"]] = tally.get(p["verdict"], 0) + 1
    return {
        "n": len(per),
        "by_verdict": tally,
        "per_element": per,
        "reading": (
            "the fourth frame's first branch rejected these for being non-coding when they are not."
            " Under the current symbol most of them pass every step of the rule, so the fourth"
            " frame's draw is smaller than its own rule specifies. The correction is the fourth"
            " frame owner's to make; it is stated here as a row, not applied"
        ),
    }


def with_nearest_coding(rows: list[dict[str, Any]], results_dir: Path = RESULTS_DIR) -> list[dict[str, Any]]:
    """Attach the nearest coding TSS to each row, the way `loci_fourth.select` does. Free."""
    out = []
    for chrom in sorted({r["chrom"] for r in rows}, key=lambda c: loci_fourth.CHROM_ORDER.get(c, 99)):
        ch = loci.Chromosome(chrom, results_dir)
        try:
            for r in (x for x in rows if x["chrom"] == chrom):
                near = ch.nearest_coding(r["mid"])
                out.append(
                    {
                        **r,
                        "nearest_coding": near[0] if near else None,
                        "nearest_coding_distance": near[1] if near else None,
                        "length": ch.length,
                    }
                )
        finally:
            ch.close()
    return out


def frame(rows: list[dict[str, Any]]) -> tuple[Expect, ...]:
    """The genuinely non-coding elements as `Expect`s, through `loci_fourth.as_expect`."""
    out = []
    for r in rows:
        verdict = {
            "locus": loci_fourth.locus_name(r),
            "nearest_coding": r["nearest_coding"],
            "nearest_coding_distance": r["nearest_coding_distance"],
        }
        e = loci_fourth.as_expect(r, verdict, r["length"])
        kinds = ", ".join(f"{x['published_symbol']} ({x['gene_type']})" for x in r["resolved"])
        out.append(
            Expect(
                **{
                    **{k: getattr(e, k) for k in e.__dataclass_fields__},
                    "note": (
                        f"{e.note}. This locus is in the NON-CODING frame, not the geometry one: its"
                        f" published target is {kinds}. Two of the four derived layers cannot name a"
                        " non-coding gene at all, so a miss here is a fact about the reader before"
                        " it is a fact about the model"
                    ),
                }
            )
        )
    return tuple(out)


def plan(panel: tuple[Expect, ...], results_dir: Path = RESULTS_DIR) -> list[dict[str, Any]]:
    """Reach and cost per locus. Raises rather than let the frame spend a request it did not budget."""
    rows: list[dict[str, Any]] = []
    for e in panel:
        ch = loci.Chromosome(e.chrom, results_dir)
        try:
            reach = loci.read_reach(ch, e)
            s, en = e.element
            scored = loci._deletion_rows(e.chrom, s, en, results_dir)  # noqa: SLF001
            rows.append(
                {
                    "locus": e.locus,
                    "chrom": e.chrom,
                    "element": [s, en],
                    "targets": list(e.targets),
                    "trap": e.nearest_gene_trap,
                    "distance": e.distance,
                    "askable": reach["askable"],
                    "targets_in_reach": reach["targets_in_reach"],
                    "targets_out_of_reach": reach["targets_out_of_reach"],
                    "already_scored_over_element": len(scored),
                    "graded": reach["askable"],
                    "requests": 0 if (not reach["askable"] or scored) else 1,
                }
            )
        finally:
            ch.close()
    spend = sum(r["requests"] for r in rows)
    if spend > REQUEST_BUDGET:
        raise RuntimeError(
            f"this frame registered a budget of {REQUEST_BUDGET} requests and would spend {spend}."
            " The registration says a run that needed a request is abandoned and reported as"
            " unaffordable, not funded. Nothing has been scored"
        )
    return rows


# ------------------------------------------------------------------------------ the readings
def layer_universe(panel: tuple[Expect, ...], results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """Per layer, whether it CAN name a non-coding gene, checked against the data rather than asserted.

    The registration calls two of these zero by construction. This recomputes that from the rows the
    layers actually read, so the claim is checkable by anyone who opens the result file.
    """
    table = gene_table()
    by_symbol = {v["symbol"]: v["gene_type"] for v in table.values()}
    out: dict[str, Any] = {}
    for e in panel:
        s, en = e.element
        rows = loci._deletion_rows(e.chrom, s, en, results_dir)  # noqa: SLF001
        any_gene = [(x.get("predicted") or {}).get("gene") for x in rows]
        coding_only = [(x.get("predicted_coding") or {}).get("gene") for x in rows]
        out[e.locus] = {
            "elements_over_the_element": len(rows),
            "predicted_any_gene": [{"gene": g, "gene_type": by_symbol.get(g)} for g in any_gene if g],
            "predicted_coding": [{"gene": g, "gene_type": by_symbol.get(g)} for g in coding_only if g],
            "deletion_reaches_the_any_gene_key": not any(coding_only),
        }
    return {
        "per_locus": out,
        "gene_input": {
            "can_name_a_non_coding_gene": False,
            "why": "`loci.read_gene_input` reads `predicted_coding` and nothing else",
        },
        "node": {
            "can_name_a_non_coding_gene": False,
            "why": "the heuristic is the nearest coding TSS inside the CTCF node",
        },
        "deletion": {
            "can_name_a_non_coding_gene": "only where `predicted_coding` is empty",
            "why": "`loci.read_deletion` walks ('predicted_coding', 'predicted') and breaks at the"
            " first key with a gene, so a coding prediction hides the non-coding one",
        },
        "eqtl": {
            "can_name_a_non_coding_gene": True,
            "why": "GTEx gene ids go through `loci.symbols`, the whole GENCODE table, lncRNAs"
            " included. Unasked in this frame by registration: see the module docstring",
        },
    }


def counterfactual(panel: tuple[Expect, ...], results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """The registered baseline: what the deletion layer would name WITHOUT the coding-first break.

    Not a scorer and not a change to one. It replays `loci.read_deletion`'s own loop over the same
    rows with the `break` removed, which is the one-line fix section 18 left to the reader's owner.
    """
    per = []
    hit = 0
    for e in panel:
        s, en = e.element
        rows = loci._deletion_rows(e.chrom, s, en, results_dir)  # noqa: SLF001
        genes: dict[str, float] = {}
        for x in rows:
            for key in ("predicted_coding", "predicted"):
                p = x.get(key) or {}
                if not p.get("gene"):
                    continue
                g, v = p["gene"], p["log2_fold_change"]
                if abs(v) > abs(genes.get(g, 0.0)):
                    genes[g] = v
        ranked = sorted(genes, key=lambda g: -abs(genes[g]))
        named = ranked[0] if ranked else None
        ok = named in e.targets
        hit += int(ok)
        per.append(
            {
                "locus": e.locus,
                "targets": list(e.targets),
                "would_name": named,
                "effect": genes.get(named) if named else None,
                "hit": ok,
                "ranked": [{"gene": g, "effect": genes[g]} for g in ranked[:5]],
            }
        )
    return {
        "k": hit,
        "n": len(panel),
        "rate": round(hit / len(panel), 4) if panel else None,
        "per_locus": per,
        "reading": (
            "this is the model's own any-gene prediction, which the sweep computed and stored and"
            " which `read_deletion` declines to look at whenever a coding prediction exists. The gap"
            " between this and the derived rate is the READER's defect, measured. It is reported as"
            " a baseline, not as a score: no hit here was earned by any code that ships"
        ),
    }


def what_the_layers_named(result: dict[str, Any]) -> dict[str, Any]:
    """Per layer: how often it named the published target, and what it named instead. As section 22."""
    rows = result["loci"]
    layers: dict[str, dict[str, Any]] = {}
    for r in rows:
        for name, v in (
            ((r.get("score") or {}).get("scored") or {}).get("target", {}).get("by_layer", {}).items()
        ):
            layers.setdefault(name, {"provenance": v["provenance"], "hit": 0, "n": 0, "pending": 0})
            layers[name]["n"] += 1
            layers[name]["hit"] += int(v["hit"])
            layers[name]["pending"] += int(bool(v.get("pending")))
    instead: dict[str, dict[str, int]] = {}
    for layer in ("deletion", "eqtl", "gene_input", "node"):
        tally = {"the published target": 0, "the nearest coding TSS": 0, "another gene": 0, "nothing": 0}
        for r in rows:
            # the gene the HIT RULE reads, not the layer's own `target` field. For `deletion` and
            # `node` the two are the same. For `gene_input` they are not: its `target` prefers a
            # published target whenever one is anywhere in the summed window, so reading it here
            # would report a hit the scorer did not give - it named TFRC at the chr3 locus while
            # ranking another gene first. `by_layer[...]["named"]` is `score_target`'s own list.
            scored = ((r.get("score") or {}).get("scored") or {}).get("target", {}).get("by_layer", {})
            named_list = (scored.get(layer) or {}).get("named") or []
            named = named_list[0] if named_list else None
            trap = r["expected"]["nearest_gene_trap"]
            if named is None:
                tally["nothing"] += 1
            elif named in r["expected"]["targets"]:
                tally["the published target"] += 1
            elif named == trap:
                tally["the nearest coding TSS"] += 1
            else:
                tally["another gene"] += 1
        instead[layer] = tally
    return {
        "by_layer": layers,
        "what_it_named_instead": instead,
        "reading": (
            "read this against `layer_universe` and not on its own. `gene_input` and `node` cannot"
            " name a non-coding gene at all, so their zeros are the denominator speaking, and"
            " `eqtl` is unasked here by registration. Only `deletion` had a path to a hit, and only"
            " where the element has no coding prediction"
        ),
    }


def run(results_dir: Path = RESULTS_DIR, network: bool = True, progress=None) -> dict[str, Any]:
    """Draw the frame, check it costs nothing, then read and score it through `loci.build`."""
    say = progress or (lambda _m: None)
    branch = dropped_for_non_coding()
    say(f"the branch stops at {len(branch)} element(s)")
    resolved = resolve(branch)
    split = classify(resolved)
    say(f"{len(split['non_coding'])} genuinely non-coding, {len(split['stale_symbol'])} stale symbols")
    rows = with_nearest_coding(split["non_coding"], results_dir)
    panel = frame(rows)
    cost = plan(panel, results_dir)
    say(f"{sum(r['requests'] for r in cost)} request(s) to spend")
    with loci_fourth.keep_out_all_four_sets(panel):
        out = loci.build(
            panel=panel,
            results_dir=results_dir,
            gtex_dir=GTEX_DIR,
            network=network,
            negatives=False,
            progress=progress,
        )
    out["result"] = NAME
    out["preregistration"] = PREREGISTRATION
    out["branch"] = {
        "n": len(branch),
        "non_coding": len(split["non_coding"]),
        "stale_symbol": len(split["stale_symbol"]),
        "resolved": [
            {
                "chrom": r["chrom"],
                "element": [r["start"], r["end"]],
                "cell": r["cell"],
                "targets": sorted(r["targets"]),
                "resolved": r["resolved"],
            }
            for r in resolved
        ],
    }
    out["corrected_rule"] = corrected_rule(with_nearest_coding(split["stale_symbol"], results_dir))
    out["plan"] = cost
    out["layer_universe"] = layer_universe(panel, results_dir)
    out["baseline"] = counterfactual(panel, results_dir)
    out["layers"] = what_the_layers_named(out)
    out["requests_spent"] = sum(r["requests"] for r in cost)
    out["note"] = (
        "The FIFTH registered frame, and the first whose published targets are non-coding. It is"
        " n = 2, not the 14 the fourth frame's branch counted: twelve of those fourteen are stale"
        " gene symbols, which the benchmark file's own Ensembl ids settle. Two of the four derived"
        " layers cannot express a non-coding answer at all and one is unasked, so the derived rate"
        " here is very nearly arithmetic and is reported as such. Zero requests. " + loci.NOTE
    )
    save_result(NAME, out, results_dir)
    return out
