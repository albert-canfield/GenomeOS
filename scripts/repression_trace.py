# SPDX-License-Identifier: AGPL-3.0-or-later
"""The four repressions S4 judged wrong, traced to where their direction came from (item 12, S4 follow-up).

    uv run python scripts/repression_trace.py

S4 (`data/results/attribution_correctness.json`, `genomeos/attribution/correctness.py`) judged every
`represses_target` claim any held evidence reached. There were four, and all four came out `incorrect`
with the detail "the gene moved the other way". Each is judged in a cell other than the one the claim
names, and each carries `refutable: false`. This script traces each claim back to the compiled rule, the
run row and the cached model answer that produced it, and traces the deciding observation back to its raw
benchmark row. It does the same for the one `activates_target` claim S4 judged wrong (CCND1, HCT116), as a
comparison.

An internal development trace of existing claims. It is not a validation, it builds no attribution
machinery and no scorer, and it changes no label and no verdict. It makes no model request. It reads the
per-element response cache (the model's answers already on disk) and two rows of the CRISPRi benchmark's
held-out file, only to trace verdicts S4 already decided with them. Nothing is fitted, selected or
relabelled from either. Writes data/results/repression_trace.json.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

NAME = "repression_trace"

# ==================================================================================================
# The registration (2026-09-29, lane-repress). The cause classes and the rule that picks the primary
# class are fixed here, before any trace code was written and before the result was computed.
# ==================================================================================================
REGISTERED = "2026-09-29"
STATUS = "an internal development trace of existing claims, not a validation"
WRITTEN_AFTER = (
    "The brief asked for the classes before any tracing. This lane had already read the five compiled "
    "rules, their rows in the whole-chromosome run, the five cached model answers and the matching raw "
    "CRISPRi rows when it wrote them down, so this registration is not blind to the data. Classes (a) to "
    "(e) are the brief's five, with its definitions. Class (f) and the primary-class rule were "
    "written after that first reading, and are marked as such."
)
QUESTION = (
    "for each represses_target claim S4 judged incorrect (four) and the one activates_target claim it judged "
    "incorrect (CCND1), where did the claimed direction come from, and why does it disagree with the "
    "observation that decided the verdict"
)

#: the claims traced, exactly as S4's result states them (element, locus, gene, the claim's cell, and the
#: deciding observation's source, gene, cell and split)
CASES = (
    {
        "element": "EH38E2447945",
        "locus": "chr6:14095446-14095736",
        "value": "represses_target",
        "gene": "CD83",
        "cell": "K562",
        "deciding": {"source": "crispri:Nasser2021", "gene": "CD83", "cell": "GM12878", "split": "heldout"},
    },
    {
        "element": "EH38E3896696",
        "locus": "chr9:97960870-97961216",
        "value": "represses_target",
        "gene": "HEMGN",
        "cell": "Brain_Cerebellar_Hemisphere",
        "deciding": {"source": "crispri:Gasperini2019", "gene": "HEMGN", "cell": "K562", "split": "training"},
    },
    {
        "element": "EH38E3426791",
        "locus": "chr20:31608511-31608860",
        "value": "represses_target",
        "gene": "ID1",
        "cell": "Whole_Blood",
        "deciding": {"source": "crispri:Gasperini2019", "gene": "ID1", "cell": "K562", "split": "training"},
    },
    {
        "element": "EH38E3940111",
        "locus": "chrX:103233290-103233585",
        "value": "represses_target",
        "gene": "BEX4",
        "cell": "hair_follicular_keratinocyte",
        "deciding": {"source": "crispri:Gasperini2019", "gene": "BEX4", "cell": "K562", "split": "training"},
    },
)
#: the comparison case: not one of the four
COMPARISON = {
    "element": "EH38E1548926",
    "locus": "chr11:69588295-69588640",
    "value": "activates_target",
    "gene": "CCND1",
    "cell": "CD8_positive__alpha_beta_memory_T_cell",
    "deciding": {"source": "crispri:HCT116", "gene": "CCND1", "cell": "HCT116", "split": "heldout"},
}

SIGN, CELL, CONTRADICTION, INDIRECT, JUDGE, SPLIT = "a", "b", "c", "d", "e", "f"
PRESENT, ABSENT, UNKNOWN = "present", "absent", "unknown"
CLASSES = {
    SIGN: {
        "name": "sign or convention error",
        "definition": "a sign flipped, or a convention misread, somewhere between the source and the label, "
        "or between the raw measurement and the observation. On the label side the steps are: the model's "
        "log2 fold change (alternate over reference, where the alternate allele is the deletion), the "
        "choice between the gene's largest drop and its largest rise (enhancer_target.predict_target), the "
        "rule verb the compiler writes (activates or inhibits), and the claim value S4 reads from the verb. "
        "On the measurement side the steps are: the benchmark's EffectSize and Significant columns, the "
        "pair's outcome (measured.CrispriPair.outcome), the observation kind (correctness.CRISPRI_KIND) "
        "and the S4 table's reading of that kind on the claim's value",
        "decided_by": "step by step: present when a step's stored sign or direction is not the one the "
        "previous step's value implies under that step's own stated convention. The label side is "
        "recomputed from the cached answer with predict_target and the compiler's verb rule, and compared "
        "with the stored run row, the compiled rule line and the S4 claim. The measurement side is "
        "recomputed from the raw benchmark row and compared with the observation S4 recorded",
        "added_after_first_reading": False,
    },
    CELL: {
        "name": "cell mismatch",
        "definition": "the claim's stated cell and the deciding observation's cell differ, so by the S4 "
        "table the observation may not bear on the claim at all. The table's CRISPRi activity and context "
        "cells speak of 'the cell screened'",
        "decided_by": "present when no deciding observation's cell equals the claim's cell after S4's own "
        "normalisation (correctness.norm_cell)",
        "added_after_first_reading": False,
    },
    CONTRADICTION: {
        "name": "genuine contradiction",
        "definition": "the same reading in the same cell has the opposite direction: the model's own value "
        "for the named gene on the measured cell's track has the opposite sign to the measured effect",
        "decided_by": "the cached answer's value for the gene on the measured cell's track: `by_cell`, and "
        "the largest drop or rise where the track it came from is named for that cell. Present when a value "
        "for that cell is opposite in sign to the measurement. Absent when every such value has the "
        "measurement's sign. Unknown when the cache holds no value for that cell. The cache keeps one "
        "`by_cell` value per cell name, so where a cell has several tracks the others are not seen, and "
        "this is stated beside each case",
        "added_after_first_reading": False,
    },
    INDIRECT: {
        "name": "indirect or ambiguous",
        "definition": "another gene may carry the effect: another gene in the window, a shared promoter or "
        "a multi-target element",
        "decided_by": "present when any of: the model's any-gene prediction names a different gene from its "
        "coding prediction; the run's compact verdict says 'another gene in the same domain'; or the same "
        "screen changed another gene significantly on silencing this element (a multi-target element). "
        "Each of the three is reported separately, and the annotated genes nearest the element are listed "
        "as description only",
        "added_after_first_reading": False,
    },
    JUDGE: {
        "name": "judge inconsistency",
        "definition": "the verdict `incorrect` alongside `refutable: false`: the judge gave a verdict its "
        "own registered rules (correctness.TABLE, correctness.AXIS_RULE) do not give, or its rules "
        "contradict each other on this claim",
        "decided_by": "correctness.verdict_of re-run on the claim with the observations S4 held for the "
        "element must give the committed verdict and the committed deciding observations. What `refutable` "
        "means is read from correctness.py and stated. Separately, the committed S4 result is counted per "
        "axis by verdict and by `refutable`, to show where `incorrect` can occur with `refutable` false",
        "added_after_first_reading": False,
    },
    SPLIT: {
        "name": "direction split across tracks",
        "definition": "the model moves the named gene both ways on different tracks, and the label takes "
        "the direction and the cell of the single most extreme track (enhancer_target.predict_target: the "
        "larger of the largest drop and the largest rise over every track wins)",
        "decided_by": "present when the cached answer records both a drop and a rise for the named gene "
        "(max_drop_log2fc < 0 < max_rise_log2fc). The live scorer records only values beyond 0.05, so both "
        "being non-zero means both passed that recording threshold. The winning margin is reported",
        "added_after_first_reading": True,
    },
}
PRIMARY_RULE = (
    "written after the first reading. The primary class is the first that applies, in this order: (a) "
    "when present, because a flipped sign explains any disagreement; (e) when present; (c) when present, "
    "because the disagreement then survives aligning the cells; (b) when present and (c) is absent or "
    "unknown, because the claim and the observation are then in different cells and the model's own "
    "reading of the measured cell does not contradict the measurement, or cannot be read. (d) and (f) are "
    "recorded beside and are never primary: S4 matches the deciding observation to the named gene, so "
    "another gene cannot be what disagrees, and (f) says where a direction came from, not why a "
    "measurement disagrees with it"
)
NO_CHANGE = (
    "no label and no verdict is changed by this trace. A bug, if one is found, is fixed with a test that "
    "fails before the fix. When the fix would move a committed result's figures, the before-and-after "
    "counts are written and the result is not regenerated: that is the coordinator's call"
)
MAY_BE_CALLED = (
    "an internal development trace of five existing claims: where each claimed direction came from and "
    "which registered class explains its disagreement with the observation that decided it"
)
MAY_NOT_BE_CALLED = (
    "a validation of the labels, the model or the judge",
    "a rate: five claims are traced, chosen because S4 judged them wrong, and no share of anything is "
    "estimated from them",
    "evidence about claims S4 did not judge",
)


def registration() -> dict[str, Any]:
    """The registered constants, as the result file states them."""
    return {
        "registered": REGISTERED,
        "status": STATUS,
        "written_after": WRITTEN_AFTER,
        "question": QUESTION,
        "cases": list(CASES),
        "comparison": COMPARISON,
        "classes": CLASSES,
        "primary_rule": PRIMARY_RULE,
        "no_change": NO_CHANGE,
        "may_be_called": MAY_BE_CALLED,
        "may_not_be_called": list(MAY_NOT_BE_CALLED),
    }


def primary(classes: dict[str, str]) -> str:
    """The primary class under PRIMARY_RULE, from each class's reading (present, absent or unknown)."""
    for k in (SIGN, JUDGE, CONTRADICTION):
        if classes.get(k) == PRESENT:
            return k
    if classes.get(CELL) == PRESENT:
        return CELL
    raise ValueError(f"no registered class is primary for {classes}")


# ==================================================================================================
# The trace. Everything below reads existing files and calls existing functions; nothing is fitted.
# ==================================================================================================
S4_RESULT = Path("data/results/attribution_correctness.json")
PILOT_RESULT = Path("data/results/pilot_biological_gate.json")
RUN = "enhancer_targets_all"
NEAR_BP = 5_000  # genes listed as description only when their annotated body lies this close
RULE_RE = re.compile(
    r"^rule (?P<name>\S+) (?P<verb>\w+) (?P<gene>\S+) \{ strength: (?P<strength>[0-9.]+); "
    r"when: cell_type = (?P<cell>\S+); evidence: (?P<evidence>\w+)"
)
EFFECT_RE = re.compile(r"effect (?P<effect>[+-][0-9.]+) log2 fold change on deletion")
CONVENTIONS = {
    "model": "AlphaGenome RECOMMENDED_VARIANT_SCORERS['RNA_SEQ'] is GeneMaskLFCScorer, whose docstring reads "
    "'the log fold change of gene expression level between ALT and REF alleles'; "
    "genomeos.predict.alphagenome_adapter.EFFECT_UNIT says 'alternate allele over reference'",
    "deletion": "enhancer_target.score_element passes the reference as the anchor base plus the element and "
    "the alternate as the anchor base alone, so the alternate allele is the deletion and a positive log2 "
    "fold change is a rise on deleting the element",
    "direction": "enhancer_target.predict_target: 'Loss of expression on deletion means the element "
    "activates the gene; a rise means it represses it. The larger of the two wins; ties go to activation'",
    "verb": "attribution/compile.py writes `activates` when the prediction's action is activates and "
    "`inhibits` otherwise, with `when: cell_type` = compile.context(the prediction's tissue)",
    "claim": "correctness.compiled_claims reads `activates_target` from the verb activates and "
    "`represses_target` from any other verb",
    "crispri": "measured.EFFECT_UNIT: EffectSize is the fractional change in the measured gene on "
    "silencing, so -0.2 is a 20% decrease; CrispriPair.outcome is a significant decrease when Significant is "
    "TRUE and EffectSize < 0; the benchmark's Regulated means a significant decrease only",
}


def _locus(s: str) -> tuple[str, int, int]:
    chrom, rest = s.split(":")
    a, b = rest.split("-")
    return chrom, int(a), int(b)


def _sign(x: float) -> int:
    return (x > 0) - (x < 0)


def s4_wrong_directions(s4: dict[str, Any]) -> list[dict[str, Any]]:
    """Every activity verdict S4 committed as incorrect."""
    return [v for v in s4["judged"] if v["axis"] == "activity" and v["verdict"] == "incorrect"]


# --- the label side ----------------------------------------------------------------------------------
def compiled_lines(program: Path, element: str) -> dict[str, Any]:
    """The element block, its predicted rule and its measured rules, with 1-based line numbers."""
    out: dict[str, Any] = {"program": str(program), "element_block": [], "rule": None, "measured_rules": []}
    inside = False
    with open(program) as fh:
        for i, line in enumerate(fh, 1):
            s = line.rstrip("\n")
            if s == f"element {element} {{":
                inside = True
            if inside:
                out["element_block"].append([i, s])
                inside = s != "}"
                continue
            if s.startswith(f"rule {element} "):
                out["rule"] = [i, s]
            elif s.startswith(f"rule {element}_measured "):
                out["measured_rules"].append([i, s])
    return out


def coding_genes(chrom: str) -> set[str]:
    """Protein-coding symbols in GENCODE v50 for the chromosome (the compact run's coding restriction)."""
    return {g["gene"] for g in annotated_genes(chrom) if g["type"] == "protein_coding"}


_GENES: dict[str, list[dict[str, Any]]] = {}


def annotated_genes(chrom: str) -> list[dict[str, Any]]:
    """GENCODE v50 gene records of one chromosome, 0-based half-open."""
    if chrom in _GENES:
        return _GENES[chrom]
    out = []
    with gzip.open(ROOT / f"data/reference/gencode_v50_{chrom}.gff3.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            out.append(
                {
                    "gene": attrs.get("gene_name", attrs.get("gene_id", "")),
                    "type": attrs.get("gene_type", ""),
                    "start": int(f[3]) - 1,
                    "end": int(f[4]),
                    "strand": f[6],
                }
            )
    _GENES[chrom] = out
    return out


def genes_near(chrom: str, start: int, end: int, named: str) -> dict[str, Any]:
    """Description only: the named gene's annotated body and TSS, and every gene within NEAR_BP."""

    def gap(g: dict[str, Any]) -> int:
        return max(0, g["start"] - end, start - g["end"])

    genes = annotated_genes(chrom)
    named_rows = [g for g in genes if g["gene"] == named]
    near = sorted((g for g in genes if gap(g) <= NEAR_BP), key=gap)
    out: dict[str, Any] = {
        "near_bp": NEAR_BP,
        "within_near_bp": [{**g, "gap_bp": gap(g)} for g in near],
    }
    if named_rows:
        g = named_rows[0]
        tss = g["start"] if g["strand"] == "+" else g["end"] - 1
        out["named_gene"] = {**g, "tss": tss, "element_to_tss_bp": (start + end) // 2 - tss, "gap_bp": gap(g)}
    return out


def label_trace(case: dict[str, Any]) -> dict[str, Any]:
    """Where the claimed direction came from: cached answer -> run row -> compiled rule -> S4 claim."""
    from genomeos.attribution import compile as cp
    from genomeos.attribution import correctness as co
    from genomeos.attribution import holdout as ho
    from genomeos.attribution import targets as tg
    from genomeos.predict import enhancer_target as et

    chrom, start, end = _locus(case["locus"])
    eid, gene = case["element"], case["gene"]
    rows = [e for e in tg.attributed(chrom) if e["id"] == eid]
    if len(rows) != 1:
        raise ValueError(f"{eid}: {len(rows)} attributed rows")
    row = rows[0]
    pc = row["predicted_coding"]
    hit = et.load_cached(chrom, eid)
    if hit is None:
        raise FileNotFoundError(f"{eid} is not in the response cache")
    cache_file = et.cache_path(chrom, eid)
    genes = hit["genes"]
    g = next(r for r in genes if r["gene"] == gene)
    coding = coding_genes(chrom)
    again_coding = et.predict_target([r for r in genes if r["gene"] in coding])
    again_any = et.predict_target(genes)
    lines = compiled_lines(ho.COMPILED_DIR / f"noncoding_{chrom}.bio", eid)
    m = RULE_RE.match(lines["rule"][1]) if lines["rule"] else None
    if m is None:
        raise ValueError(f"{eid}: no predicted rule line")
    eff = EFFECT_RE.search(lines["rule"][1])
    claims = [
        c
        for c in co.compiled_claims(ho.COMPILED_DIR, [chrom])
        if c.element == eid and c.gene == gene and c.axis == co.ACTIVITY
    ]
    if len(claims) != 1:
        raise ValueError(f"{eid}: {len(claims)} activity claims")
    claim = claims[0]
    drop, rise = g["max_drop_log2fc"], g["max_rise_log2fc"]
    implied = "activates" if -drop >= rise else "represses"
    verb_expected = "activates" if pc["action"] == "activates" else "inhibits"
    value_expected = "activates_target" if m["verb"] == "activates" else "represses_target"
    steps = [
        {
            "step": "L1 cached model answer",
            "where": f"{cache_file if cache_file.exists() else et.archive_path(chrom)} [{eid}] genes[{gene}]",
            "stored": {k: g[k] for k in g if k != "gene"},
            "implies": f"{implied}: the larger of the largest drop ({drop:+.4f}, {g['max_drop_tissue']}) and "
            f"the largest rise ({rise:+.4f}, {g['max_rise_tissue']})",
            "consistent": True,
            "convention": CONVENTIONS["model"] + "; " + CONVENTIONS["deletion"],
        },
        {
            "step": "L2 run row (predicted_coding)",
            "where": f"{RUN}_{chrom} -> {row.get('origin')} run table, element {eid}",
            "stored": {k: pc.get(k) for k in ("gene", "action", "log2_fold_change", "tissue", "strength")},
            "recomputed": {
                k: (again_coding or {}).get(k) for k in ("gene", "action", "log2_fold_change", "tissue")
            },
            "consistent": bool(again_coding)
            and all(again_coding[k] == pc[k] for k in ("gene", "action", "tissue"))
            and abs(again_coding["log2_fold_change"] - pc["log2_fold_change"]) < 1e-9
            and pc["action"] == implied
            and (pc["action"] == "represses") == (pc["log2_fold_change"] > 0),
            "convention": CONVENTIONS["direction"],
        },
        {
            "step": "L3 compiled rule",
            "where": f"{lines['program']}:{lines['rule'][0]}",
            "stored": {
                "verb": m["verb"],
                "gene": m["gene"],
                "cell": m["cell"],
                "evidence": m["evidence"],
                "effect": float(eff["effect"]) if eff else None,
            },
            "consistent": m["verb"] == verb_expected
            and m["gene"] == cp.ident(pc["gene"])
            and m["cell"] == cp.context(pc.get("tissue"))
            and m["evidence"] == "predicted"
            and eff is not None
            and _sign(float(eff["effect"])) == _sign(pc["log2_fold_change"]),
            "convention": CONVENTIONS["verb"],
        },
        {
            "step": "L4 S4 claim",
            "where": "correctness.compiled_claims",
            "stored": {"value": claim.value, "gene": claim.gene, "cell": claim.cell},
            "consistent": claim.value == value_expected
            and claim.value == case["value"]
            and claim.cell == m["cell"] == case["cell"],
            "convention": CONVENTIONS["claim"],
        },
    ]
    return {
        "source_run": row.get("origin"),
        "run_row": {
            k: row.get(k)
            for k in (
                "id",
                "start",
                "end",
                "predicted",
                "predicted_coding",
                "predicted_by_cell",
                "predicted_coding_by_cell",
                "verdict_coding",
                "inferred",
                "model_version",
            )
        },
        "cached_gene_row": g,
        "cached_tracks": hit.get("tracks"),
        "cached_genes_in_window": hit.get("genes_in_window"),
        "any_gene_prediction_recomputed": again_any,
        "compiled": lines,
        "steps": steps,
        "genes_near": genes_near(chrom, start, end, gene),
    }


# --- the measurement side ----------------------------------------------------------------------------
def benchmark_rows(chrom: str, start: int, end: int) -> list[dict[str, Any]]:
    """Every valid-or-not benchmark row the element's interval measures (measured.measures), in either
    file, with its 1-based line number in the decompressed table (the header is line 1)."""
    from genomeos.attribution import measured as ms

    out = []
    for name in ms.CRISPRI_FILES:
        p = ms.CRISPRI_KNOWLEDGE / name
        with gzip.open(p, "rt") as fh:
            header = fh.readline()
            cols = header.rstrip("\n").split("\t")
            for i, line in enumerate(fh, 2):
                if not line.startswith(chrom + "\t"):
                    continue
                r = dict(zip(cols, line.rstrip("\n").split("\t"), strict=False))
                if r["chrom"] != chrom or not ms.measures(
                    start, end, int(r["chromStart"]), int(r["chromEnd"])
                ):
                    continue
                pairs, invalid = ms.parse_crispri([header, line], None, ms.CRISPRI_SPLIT_OF[name], name)
                keep = (
                    "chrom",
                    "chromStart",
                    "chromEnd",
                    "name",
                    "measuredGeneSymbol",
                    "measuredGeneEnsemblId",
                    "EffectSize",
                    "Significant",
                    "pValueAdjusted",
                    "Regulated",
                    "PowerAtEffectSize20",
                    "ValidConnection",
                    "CellType",
                    "Dataset",
                    "Reference",
                )
                out.append(
                    {
                        "file": name,
                        "split": ms.CRISPRI_SPLIT_OF[name],
                        "line": i,
                        "row": {k: r.get(k) for k in keep},
                        "outcome": pairs[0].outcome if pairs else None,
                        "invalid": bool(invalid),
                    }
                )
    return out


def measurement_trace(case: dict[str, Any], rows: list[dict[str, Any]]) -> dict[str, Any]:
    """The deciding observation back to its raw row: sign, size, significance and cell at each step."""
    from genomeos.attribution import correctness as co
    from genomeos.attribution import measured as ms

    d = case["deciding"]
    dataset = d["source"].split(":", 1)[1]
    hits = [
        r
        for r in rows
        if r["row"]["measuredGeneSymbol"] == d["gene"]
        and r["row"]["CellType"] == d["cell"]
        and r["row"]["Dataset"] == dataset
        and r["split"] == d["split"]
        and not r["invalid"]
    ]
    if len(hits) != 1:
        raise ValueError(f"{case['element']}: {len(hits)} raw rows for the deciding observation")
    h = hits[0]
    es = float(h["row"]["EffectSize"])
    sig = h["row"]["Significant"] == "TRUE"
    expected_outcome = (ms.DECREASE if es < 0 else ms.INCREASE) if sig else None
    kind = co.CRISPRI_KIND[h["outcome"]]
    reading = co.reading(kind, co.ACTIVITY, case["value"])
    others = [
        {
            "gene": r["row"]["measuredGeneSymbol"],
            "line": r["line"],
            "effect_size": float(r["row"]["EffectSize"]),
            "p_adjusted": float(r["row"]["pValueAdjusted"]),
            "outcome": r["outcome"],
        }
        for r in rows
        if r["row"]["measuredGeneSymbol"] != d["gene"]
        and r["row"]["CellType"] == d["cell"]
        and r["row"]["Dataset"] == dataset
        and not r["invalid"]
    ]
    steps = [
        {
            "step": "M1 raw benchmark row",
            "where": f"{h['file']} line {h['line']} (name {h['row']['name']})",
            "stored": h["row"],
            "consistent": True,
            "convention": CONVENTIONS["crispri"],
        },
        {
            "step": "M2 pair outcome",
            "where": "measured.parse_crispri -> CrispriPair.outcome",
            "stored": h["outcome"],
            "consistent": h["outcome"] == expected_outcome
            and (h["row"]["Regulated"] == "TRUE") == (h["outcome"] == ms.DECREASE),
            "convention": "significant and EffectSize < 0 is a decrease, > 0 an increase; Regulated TRUE "
            "only "
            "for a significant decrease",
        },
        {
            "step": "M3 observation kind",
            "where": "correctness.CRISPRI_KIND",
            "stored": kind,
            "consistent": kind in {"crispri_decrease", "crispri_increase"},
            "convention": "one kind per outcome",
        },
        {
            "step": "M4 S4 table reading on the claim",
            "where": f"correctness.TABLE[{kind}][activity]",
            "stored": reading,
            "consistent": reading == co.REFUTES,
            "convention": "crispri_decrease refutes represses_target; crispri_increase refutes "
            "activates_target",
        },
    ]
    return {
        "deciding_row": h,
        "effect_size": es,
        "significant": sig,
        "p_adjusted": float(h["row"]["pValueAdjusted"]),
        "cell": h["row"]["CellType"],
        "steps": steps,
        "other_genes_same_screen": others,
        "all_rows_at_element": [
            {
                "file": r["file"],
                "line": r["line"],
                "gene": r["row"]["measuredGeneSymbol"],
                "cell": r["row"]["CellType"],
                "dataset": r["row"]["Dataset"],
                "effect_size": r["row"]["EffectSize"],
                "significant": r["row"]["Significant"],
                "outcome": r["outcome"],
                "valid": not r["invalid"],
            }
            for r in rows
        ],
    }


# --- the judge -----------------------------------------------------------------------------------------
def rerun(case: dict[str, Any], evidence: Any) -> dict[str, Any]:
    """correctness.verdict_of on the element's three rule claims with the observations S4 held."""
    from genomeos.attribution import correctness as co

    chrom, start, end = _locus(case["locus"])
    obs = evidence.at(chrom, start, end)
    g, cell = case["gene"], case["cell"]
    claims = {
        co.TARGET: co.Claim(case["element"], chrom, start, end, co.TARGET, g, g, cell),
        co.ACTIVITY: co.Claim(case["element"], chrom, start, end, co.ACTIVITY, case["value"], g, cell),
        co.CONTEXT: co.Claim(case["element"], chrom, start, end, co.CONTEXT, cell, g, cell),
    }
    out = {a: co.verdict_of(c, obs, co.RULE_V1).to_dict() for a, c in claims.items()}  # S4's rule, pinned
    out["observations_touching"] = sorted(
        {(o.kind, o.source, o.gene, o.cell) for o in co.touching(claims[co.ACTIVITY], obs)}
    )
    return out


def judge_block(s4: dict[str, Any]) -> dict[str, Any]:
    """What `refutable` means, where `incorrect` occurs with it false, and the counts a cell-scoped
    direction rule would give (descriptive only: no verdict is changed)."""
    from collections import Counter

    from genomeos.attribution import correctness as co

    table = Counter((v["axis"], v["verdict"], v["refutable"]) for v in s4["judged"])
    by_axis: dict[str, dict[str, int]] = {}
    for (axis, verdict, refutable), n in sorted(table.items()):
        by_axis.setdefault(axis, {})[f"{verdict}, refutable {str(refutable).lower()}"] = n
    act = [v for v in s4["judged"] if v["axis"] == co.ACTIVITY]

    def here(v: dict[str, Any]) -> bool:
        return any(co.norm_cell(o.get("cell", "")) == co.norm_cell(v["cell"]) for o in v["deciding"])

    decided_here = Counter(v["verdict"] for v in act if here(v))
    decided_elsewhere = Counter(v["verdict"] for v in act if not here(v))
    return {
        "refutable_means": "correctness.Verdict.refutable, commented 'target and context: the stated cell "
        "was screened on the gene'. AXIS_RULE[target] ends: '`refutable` counts the claims whose stated cell "
        "was screened on the named gene with an outcome that can establish or refute'. _target sets it to "
        "whether an establishing or refuting observation of the gene exists in the stated cell; _context "
        "sets it True on every decided verdict and to whether a refuting null exists there when not judged; "
        "_direction never sets it, so every activity verdict carries the default False. On the activity "
        "axis the field is not computed: false there means 'not computed', not 'could not have been "
        "refuted'. The registered quantities report it only beside target and context (AxisTally.accuracy); "
        "the per-verdict records and the per-axis tallies carry the default on every axis",
        "incorrect_reachable_with_refutable_false": {
            "target": "no: _target returns incorrect only with a well-powered null in the stated cell, and "
            "then passes refutable True",
            "context": "no: _context passes refutable True on every decided verdict",
            "activity": "yes, by the registered rule: AXIS_RULE[activity] judges the direction "
            "'in the stated "
            "cell when it responded there and otherwise in every cell it responded in', so a response in "
            "another cell can refute it, and refutable is never computed on this axis",
        },
        "committed_s4_counts_by_axis": by_axis,
        "activity_verdicts_decided_in_the_stated_cell": dict(decided_here),
        "activity_verdicts_decided_only_in_another_cell": dict(decided_elsewhere),
        "tension": "the S4 table's CRISPRi activity cells read "
        "'for the gene that fell, in the cell screened' "
        "and 'for the gene that rose, in the cell screened', and the target rule refutes only in the stated "
        "cell; the activity rule applies a response in another cell to a claim about the stated cell. The "
        "judge follows its registration exactly; whether the direction should be judged only in the stated "
        "cell is a registration question, not a bug fix, and is left to the coordinator",
        "if_direction_were_refuted_only_in_the_stated_cell": {
            "as_the_target_rule_refutes": "the incorrect activity verdicts decided only in another cell "
            "would "
            "become not judged; the correct ones would stay",
            "incorrect_that_would_move": decided_elsewhere.get(co.INCORRECT, 0),
            "correct_that_would_move": 0,
        },
        "if_direction_were_judged_only_in_the_stated_cell": {
            "judged_would_be": sum(decided_here.values()),
            "correct_would_be": decided_here.get(co.CORRECT, 0),
            "incorrect_would_be": decided_here.get(co.INCORRECT, 0),
            "moved_to_not_judged": sum(decided_elsewhere.values()),
        },
    }


# --- the classes ---------------------------------------------------------------------------------------
def classify(
    case: dict[str, Any], label: dict[str, Any], meas: dict[str, Any], again: dict[str, Any]
) -> dict:
    from genomeos.attribution import correctness as co

    reading: dict[str, str] = {}
    why: dict[str, str] = {}
    bad = [s["step"] for s in label["steps"] + meas["steps"] if not s["consistent"]]
    reading[SIGN] = PRESENT if bad else ABSENT
    why[SIGN] = (
        f"inconsistent steps: {bad}"
        if bad
        else "every step on both sides carries the sign the one before implies"
    )

    cell = meas["cell"]
    reading[CELL] = PRESENT if co.norm_cell(cell) != co.norm_cell(case["cell"]) else ABSENT
    why[CELL] = f"the claim names {case['cell']}; the deciding observation is in {cell}"

    g = label["cached_gene_row"]
    vals = [
        (f"by_cell[{k}]", v)
        for k, v in (g.get("by_cell") or {}).items()
        if co.norm_cell(k) == co.norm_cell(cell)
    ]
    if co.norm_cell(g.get("max_drop_tissue", "")) == co.norm_cell(cell):
        vals.append(("max_drop_log2fc", g["max_drop_log2fc"]))
    if co.norm_cell(g.get("max_rise_tissue", "")) == co.norm_cell(cell):
        vals.append(("max_rise_log2fc", g["max_rise_log2fc"]))
    ms_sign = _sign(meas["effect_size"])
    if not vals:
        reading[CONTRADICTION] = UNKNOWN
        why[CONTRADICTION] = (
            f"model on {cell}: the cache keeps no value for {case['gene']} on a track named for this cell"
        )
    else:
        opposite = [f"{k} {v:+.4f}" for k, v in vals if _sign(v) == -ms_sign]
        reading[CONTRADICTION] = PRESENT if opposite else ABSENT
        shown = ", ".join(f"{k} {v:+.4f}" for k, v in vals)
        why[CONTRADICTION] = f"model on {cell}: {shown}; measured EffectSize {meas['effect_size']:+.4f}" + (
            "" if len(vals) > 1 else "; one value kept, other tracks named for this cell are not seen"
        )

    rr = label["run_row"]
    any_gene = (rr.get("predicted") or {}).get("gene")
    multi = [
        o
        for o in meas["other_genes_same_screen"]
        if o["outcome"] in ("significant_decrease", "significant_increase")
    ]
    parts = {
        "any_gene_prediction_differs": any_gene != case["gene"],
        "verdict_another_gene_in_domain": rr.get("verdict_coding") == "another gene in the same domain",
        "multi_target_in_same_screen": bool(multi),
    }
    reading[INDIRECT] = PRESENT if any(parts.values()) else ABSENT
    why[INDIRECT] = (
        f"any-gene prediction {any_gene}; compact verdict '{rr.get('verdict_coding')}'; other genes changed "
        f"significantly in the same screen: {[(o['gene'], o['effect_size']) for o in multi] or 'none'}"
    )

    committed_deciding = [
        {k: o.get(k, "") for k in ("kind", "source", "gene", "cell", "split")}
        for o in case["committed"]["deciding"]
    ]
    rerun_deciding = [
        {k: o.get(k, "") for k in ("kind", "source", "gene", "cell", "split")}
        for o in again["activity"]["deciding"]
    ]
    same = (
        again["activity"]["verdict"] == case["committed"]["verdict"] and rerun_deciding == committed_deciding
    )
    reading[JUDGE] = ABSENT if same else PRESENT
    why[JUDGE] = (
        "verdict_of reproduces the committed verdict and deciding observation; the activity rule as "
        "registered lets a response in another cell decide, and refutable is not computed on activity"
        if same
        else f"re-run gives {again['activity']['verdict']} with {rerun_deciding}"
    )

    drop, rise = g["max_drop_log2fc"], g["max_rise_log2fc"]
    reading[SPLIT] = PRESENT if drop < 0 < rise else ABSENT
    margin = round(abs(rise - (-drop)), 4)
    why[SPLIT] = (
        f"largest drop {drop:+.4f} ({g['max_drop_tissue']}), "
        f"largest rise {rise:+.4f} ({g['max_rise_tissue']}); "
        f"the {'rise' if rise > -drop else 'drop'} wins by {margin}"
    )
    return {
        "reading": reading,
        "why": why,
        "indirect_parts": parts,
        "split_margin": margin,
        "primary": primary(reading),
    }


def one_line(case: dict[str, Any], label: dict[str, Any], meas: dict[str, Any], cls: dict[str, Any]) -> str:
    g = label["cached_gene_row"]
    pc = label["run_row"]["predicted_coding"]
    verb = "represses" if pc["action"] == "represses" else "activates"
    moved = "fell" if meas["effect_size"] < 0 else "rose"
    return (
        f"{case['gene']}: the model's largest rise ({g['max_rise_log2fc']:+.4f}, {g['max_rise_tissue']}) and "
        f"largest drop ({g['max_drop_log2fc']:+.4f}, {g['max_drop_tissue']}) disagree, and the "
        f"{'rise' if verb == 'represses' else 'drop'} won by {cls['split_margin']}, "
        f"so the rule says {verb} in "
        f"{case['cell']}. The screen silenced the element in {meas['cell']} and {case['gene']} {moved} "
        f"(EffectSize {meas['effect_size']:+.4f}, adjusted p {meas['p_adjusted']:.2g}). "
        f"The {cls['why'][CONTRADICTION]}. Primary class ({cls['primary']}) "
        f"{CLASSES[cls['primary']]['name']}"
    )


def split_tracks(traces: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Cases where a largest drop or rise sits on a track named for a cell whose by_cell value differs."""
    from genomeos.attribution import correctness as co

    out = []
    for t in traces:
        g = t["label_side"]["cached_gene_row"]
        for key, tissue in (("max_drop_log2fc", "max_drop_tissue"), ("max_rise_log2fc", "max_rise_tissue")):
            for cell, v in (g.get("by_cell") or {}).items():
                if co.norm_cell(cell) == co.norm_cell(g.get(tissue, "")) and v != g[key]:
                    out.append({"gene": t["gene"], "cell": cell, key: g[key], "by_cell": v})
    return out


# --- context -------------------------------------------------------------------------------------------
def predicted_verbs() -> dict[str, int]:
    """Predicted rules by verb across the compiled programs: how many repressions the labels state."""
    from collections import Counter

    from genomeos.attribution import holdout as ho

    n: Counter = Counter()
    for p in ho.compiled_programs(ho.COMPILED_DIR):
        with open(p) as fh:
            for line in fh:
                if line.startswith("rule ") and "evidence: predicted" in line:
                    m = RULE_RE.match(line)
                    if m and not m["name"].endswith("_measured"):
                        n[m["verb"]] += 1
    return dict(sorted(n.items()))


def where_ccnd1_is() -> dict[str, Any]:
    text = PILOT_RESULT.read_text()
    pilot = json.loads(text)
    return {
        "found_in": str(S4_RESULT),
        "pilot_biological_gate_mentions": {"CCND1": text.count("CCND1"), "HCT116": text.count("HCT116")},
        "pilot_sources": sorted(pilot.get("sources", {})),
        "reading": "the CCND1 activates_target claim judged wrong by an HCT116 pair of the held-out file is "
        "in the S4 result, not in the pilot's result: the pilot's sources do not include crispri:HCT116 and "
        "its result names neither CCND1 nor HCT116",
    }


# --- the manifest and the run ----------------------------------------------------------------------
def manifest(chroms: list[str], parameters: dict[str, Any]) -> dict[str, Any]:
    import s4_correctness_run as s4

    from genomeos import manifest as mf
    from genomeos.attribution import holdout as ho
    from genomeos.attribution import targets as tg
    from genomeos.predict import enhancer_target as et

    programs = [ho.COMPILED_DIR / f"noncoding_{c}.bio" for c in chroms]
    m = s4.manifest(programs, parameters)
    extra = [mf.input_entry(S4_RESULT, partition=None), mf.input_entry(PILOT_RESULT, partition=None)]
    for c in chroms:
        r = json.loads((ROOT / f"data/results/{RUN}_{c}.json").read_text())
        extra.append(mf.input_entry(ROOT / f"data/results/{RUN}_{c}.json", partition=None))
        if r.get("elements_where"):
            extra.append(mf.input_entry(r["elements_where"], partition=None))
        if et.archive_path(c).exists():
            extra.append(mf.input_entry(et.archive_path(c), partition=None))
    m["inputs"] = extra + m["inputs"]
    m["sources"] = [
        {
            "accession": "this repository: data/results/attribution_correctness.json (S4, ab62d99)",
            "version": "pinned by sha256",
        },
        {
            "accession": f"this repository: the {RUN} run tables and the per-element response cache "
            f"({tg.ELEMENT_CACHE}), local and untracked",
            "version": "pinned by sha256",
        },
        *m["sources"],
    ]
    m["exclusions"] = [
        "no model request: the per-element response cache is read, never asked",
        "no label and no verdict is changed, and attribution_correctness.json is not regenerated",
        "the cache keeps one `by_cell` value per cell name and only the largest drop and rise over every "
        "track, so the model's other tracks for a cell are not seen",
        "HCT116 has no `by_cell` value in the cache (only K562, HepG2, GM12878 and IMR-90 are kept)",
        *m["exclusions"],
    ]
    m["partitions"] = {
        "training": "the three Gasperini2019 K562 rows traced (HEMGN, ID1, BEX4) come from the benchmark's "
        "training file",
        "heldout": "two rows of the held-out file (Nasser2021 CD83 in GM12878, HCT116 CCND1) are read "
        "only to "
        "trace verdicts S4 already decided with them; nothing is fitted, selected or relabelled from them",
    }
    return m


def trace_case(case: dict[str, Any], evidence: Any) -> dict[str, Any]:
    chrom, start, end = _locus(case["locus"])
    label = label_trace(case)
    meas = measurement_trace(case, benchmark_rows(chrom, start, end))
    again = rerun(case, evidence)
    cls = classify(case, label, meas, again)
    return {
        "element": case["element"],
        "locus": case["locus"],
        "gene": case["gene"],
        "claim": {"value": case["value"], "cell": case["cell"]},
        "committed_s4_verdict": case["committed"],
        "one_line": one_line(case, label, meas, cls),
        "classes": cls["reading"],
        "primary": cls["primary"],
        "why": cls["why"],
        "indirect_parts": cls["indirect_parts"],
        "label_side": label,
        "measurement_side": meas,
        "s4_rerun": again,
    }


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--dry", action="store_true", help="print the traces, write nothing")
    a = ap.parse_args(argv)
    t0 = time.time()
    from genomeos.attribution import correctness as co
    from genomeos.attribution import holdout as ho
    from genomeos.results import save_result

    s4 = json.loads(S4_RESULT.read_text())
    wrong = s4_wrong_directions(s4)
    by_key = {(v["element"], v["gene"]): v for v in wrong}
    cases = []
    for c in (*CASES, COMPARISON):
        v = by_key.get((c["element"], c["gene"]))
        if v is None or v["value"] != c["value"] or v["cell"] != c["cell"] or v["locus"] != c["locus"]:
            raise ValueError(f"{c['element']} is not committed in S4 as registered")
        d = v["deciding"]
        if len(d) != 1 or any(d[0][k] != c["deciding"][k] for k in ("source", "gene", "cell", "split")):
            raise ValueError(f"{c['element']}: the committed deciding observation is not the registered one")
        cases.append({**c, "committed": v})
    if len(wrong) != len(cases):
        raise ValueError(f"S4 holds {len(wrong)} wrong directions, the registration names {len(cases)}")
    units = ho.all_units()
    evidence = co.Evidence(units, ho.sources(units))
    traces = [trace_case(c, evidence) for c in cases]
    for t in traces:
        print(t["one_line"], json.dumps(t["classes"]), flush=True)
    if a.dry:
        return
    four, comparison = traces[:-1], traces[-1]
    chroms = sorted({_locus(c["locus"])[0] for c in cases})
    payload = {
        "question": QUESTION,
        "status": STATUS,
        "may_be_called": MAY_BE_CALLED,
        "may_not_be_called": list(MAY_NOT_BE_CALLED),
        "registration": registration(),
        "conventions": CONVENTIONS,
        "summary": [
            {
                "element": t["element"],
                "gene": t["gene"],
                "claim": f"{t['claim']['value']} in {t['claim']['cell']}",
                "deciding": f"{t['measurement_side']['deciding_row']['row']['Dataset']}, "
                f"{t['measurement_side']['cell']}, {t['measurement_side']['deciding_row']['split']}",
                "classes": t["classes"],
                "primary": t["primary"],
                "one_line": t["one_line"],
            }
            for t in traces
        ],
        "four": four,
        "comparison": comparison,
        "comparison_shares_a_cause_with_the_four": {
            k: comparison["classes"][k] == PRESENT and any(t["classes"][k] == PRESENT for t in four)
            for k in CLASSES
        },
        "judge": judge_block(s4),
        "where_the_comparison_is": where_ccnd1_is(),
        "predicted_rules_by_verb": predicted_verbs(),
        "bugs_fixed": [],
        "other_findings": {
            "by_cell_keeps_one_track_per_cell_name": {
                "reading": "enhancer_target.aggregate writes by_cell[name] = value for every track above "
                "the scorer's 0.05 recording threshold, so where several tracks carry one cell's name the "
                "last one read is kept and the others are dropped. The cases below show a largest drop or "
                "rise on a track named for a cell that differs from that cell's by_cell value, so at least "
                "two tracks carry the name. It changes no class here (both values have one sign in each "
                "case), and no S4 verdict reads by_cell. Readers that take by_cell as 'the cell's track' "
                "(crispri_direction, closure) read one of several. Not fixed: which track or combination a "
                "cell should keep is a decision, and the cached answers cannot be changed without new model "
                "requests",
                "evidence": split_tracks(traces),
            },
        },
        "unknown": [
            "the model's value for CCND1 on an HCT116 track: the cache keeps only K562, HepG2, GM12878 and "
            "IMR-90",
            "the model's other tracks for a cell where only one by_cell value is kept (CD83 in GM12878, ID1 "
            "in K562): another K562 track for ID1 could carry the other sign",
            "how many of the 156,925 predicted repressions share class (f): the run tables keep only the "
            "winning direction, and reading the losing one needs every chromosome's response cache",
            "why the model's ID1 prediction swings from -0.74 (liver) to +0.81 (whole blood) at an element "
            "1.3 kb from MIR3193 and 1.8 kb past ID1's annotated end",
        ],
        "labels_or_verdicts_changed": 0,
        "alphagenome_requests": 0,
        "per_element_response_cache_opened": True,
        "seconds": round(time.time() - t0, 1),
    }
    p = save_result(NAME, payload, manifest=manifest(chroms, {"registration": registration(), "run": RUN}))
    print("saved", p, round(time.time() - t0, 1), "s", flush=True)


if __name__ == "__main__":
    main()
