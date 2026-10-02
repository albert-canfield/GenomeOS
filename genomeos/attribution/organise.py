# SPDX-License-Identifier: AGPL-3.0-or-later
"""The block organiser: every UNKNOWN block read with all the evidence there is, copies first.

The budget gave each block a tier from its sequence class and mammalian constraint.
Two more readings have since been measured per block: the human axis (area J,
gnomAD Gnocchi: is the block depleted of variation among people, and which of the
four cases does the pair of axes make), and the copy flag (curated segmental
duplications: is the block a copy of sequence elsewhere, where a mammalian
alignment is paralogous and human variants do not map). This module joins the
three per block and re-reads the tier in that order:

1. a copy (half or more of the block duplicated) is read as a copy before anything
   else is said about it, whatever its tier;
2. a constrained_unknown block that is not a copy is read by its case: constrained
   on both axes is the sharpest candidate for something unannotated; held across
   mammals but variable among people is a frame whose value varies, or lost
   function; measured on neither is unmeasured, not neutral;
3. the other tiers keep their label, with the copy flag and the case attached.

Nothing is recomputed: the three inputs are `budget_<chrom>`, `variation_<chrom>` and
`duplication_<chrom>`. The result keeps the per-tier tallies, the copies and the
candidates, not a second copy of every block; `blocks()` returns the full join on
demand. The real unknown of a chromosome is what remains of the constrained_unknown
tier after the copies, split by case.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any

from genomeos.attribution.budget import TIERS, read_axes
from genomeos.results import RESULTS_DIR, load_result, save_result

COPY_MIN = 0.5  # duplicated fraction from which a block is read as a copy first
CASE_READING = {
    "syntax": "constrained across mammals and among people, unannotated: the sharpest candidate",
    "relaxed": "held across mammals, variable among people: a frame whose value varies, or lost function",
    "recent": "free across mammals, constrained among people: recent function or selection",
    "tolerant": "free on both scales at kilobase resolution: check the alignment before attributing",
    "unmeasured": "no human-axis coverage (the block is under a kilobase or on chrY): mammals only",
}


def blocks(
    chrom: str, results_dir: Path = RESULTS_DIR, budget: dict[str, Any] | None = None
) -> list[dict[str, Any]]:
    """Every UNKNOWN block of the chromosome with its tier, both constraint axes, the copy flag
    and the attributed elements inside it.

    `budget` is the `read_axes` record, for a caller that has already read it and needs to know
    WHICH of the two budget files answered (`axes_source`): passing it back avoids a second read of
    up to 2.7 MB and, more to the point, keeps the declaration and the join on one reading instead
    of two that could disagree."""
    budget = budget if budget is not None else read_axes(chrom, results_dir)
    if not budget:
        raise FileNotFoundError(f"no budget_{chrom} result; run genomeos budget --chrom {chrom}")
    var = {
        (b["start"], b["end"]): b
        for b in (load_result(f"variation_{chrom}", results_dir) or {}).get("blocks", [])
    }
    dup = {
        (b["start"], b["end"]): b
        for b in (load_result(f"duplication_{chrom}", results_dir) or {}).get("blocks", [])
    }
    elements = _attributed(chrom, results_dir)
    out = []
    for b in sorted(budget["blocks"], key=lambda x: x["start"]):
        key = (b["start"], b["end"])
        v, d = var.get(key), dup.get(key)
        case = (v or {}).get("case") or {}
        gn = (v or {}).get("gnocchi") or {}
        dup_frac = (d or {}).get("duplicated_fraction")
        inside = [e for e in elements if b["start"] <= e["start"] and e["end"] <= b["end"]]
        row = {
            "start": b["start"],
            "end": b["end"],
            "length": b["length"],
            "class": b["class"],
            "tier": b["guess"]["legacy_tier"],  # the stored key organise_<chrom> joins on
            "confidence": b["guess"]["confidence"],
            "mammal_fraction": (b.get("phylop") or {}).get("fraction_above"),
            "conserved_elements": (b.get("elements") or {}).get("n"),
            "human_fraction": gn.get("fraction_above"),
            "case": case.get("case") or ("unmeasured" if v is not None else None),
            "duplicated_fraction": dup_frac,
            "partners": (d or {}).get("pairs"),
            "copy": bool(dup_frac is not None and dup_frac >= COPY_MIN),
            "attributed_elements": len(inside),
            "targets": sorted({e["target"] for e in inside}),
        }
        row["reading"] = reading(row)
        out.append(row)
    return out


def inputs(chrom: str, results_dir: Path = RESULTS_DIR) -> list[Path]:
    """The files `blocks(chrom)` reads, for a result manifest (review item R9): the three per-block
    readings and every attribution run with the local table a summary points at.

    A POINTER THAT DOES NOT RESOLVE RAISES, naming the result that carries the pointer and the path
    it names. Until 2026-10-02 the pointed-at table was appended only `if w.exists()`, so a summary
    whose table was not on this machine left the declaration one entry shorter and said nothing.
    That is worse than a declaration that never followed the pointer at all, because this one still
    claims to be complete: the manifest goes on naming the summary, `targets.run_elements` goes on
    reading the run's elements out of the table the summary points at, and nothing in the result
    distinguishes a short declaration from a short run. The tables are git-ignored and run to about
    1.4 GB, so the machine that wrote the result is the only one that has them -- which is exactly
    why the pointer has to be declared and hashed rather than dropped.

    The per-chromosome results in `names` are NOT the same case and are deliberately not required
    one by one. `blocks()` reads `variation_<chrom>` and `duplication_<chrom>` through `or {}` and
    reaches the attribution runs through `run_elements`, each of which tolerates absence, so a
    chromosome that has not had one of them run is a smaller join and not a failed one; this list
    enumerates candidates rather than following pointers, and raising per name would make a partial
    genome unwritable. The pair of budget files is the one exception, because `read_axes` returns
    None only when BOTH are absent and `blocks()` then raises: this function raises on that same
    condition and names both paths, rather than handing a writer a declaration for a join that
    cannot be built.
    """
    from genomeos.attribution.targets import RUNS

    names = [f"budget_{chrom}", f"variation_{chrom}", f"duplication_{chrom}", *(f"{r}_{chrom}" for r in RUNS)]
    names.insert(0, f"budget_axes_{chrom}")  # read_axes opens both budget files
    budget = [results_dir / f"budget_axes_{chrom}.json", results_dir / f"budget_{chrom}.json"]
    if not any(b.exists() for b in budget):
        raise FileNotFoundError(
            f"neither budget file for {chrom} is on this machine, so there is no declaration to make "
            f"for a join read_axes cannot build: {', '.join(str(b) for b in budget)}"
        )
    out = []
    for name in names:
        p = results_dir / f"{name}.json"
        if not p.exists():
            continue
        out.append(p)
        where = (load_result(name, results_dir) or {}).get("elements_where")
        if not where:
            continue
        w = Path(where)
        if not w.is_absolute() and not w.exists():
            w = results_dir.parent.parent / where
        if not w.exists():
            raise FileNotFoundError(
                f"{p} declares elements_where={where} and the table it points at is not on this "
                f"machine (resolved to {w}); this run reads that run's elements out of that table, so "
                f"a declaration without it would understate what the result depends on"
            )
        out.append(w)
    return out


def _attributed(chrom: str, results_dir: Path) -> list[dict[str, Any]]:
    from genomeos.attribution.targets import attributed

    return [
        {"id": e["id"], "start": e["start"], "end": e["end"], "target": e["predicted_coding"]["gene"]}
        for e in attributed(chrom, results_dir)
    ]


def reading(row: dict[str, Any]) -> str:
    """One line per block, copies first."""
    if row["copy"]:
        n = row.get("partners") or 0
        return (
            f"copy: {row['duplicated_fraction'] * 100:.0f}% of the block is a curated segmental duplication "
            f"({n} partner{'s' if n != 1 else ''}); read the source first, attribute second"
        )
    if row["tier"] == "constrained_unknown":
        return CASE_READING.get(row["case"] or "unmeasured", CASE_READING["unmeasured"])
    base = row["tier"]
    if row["case"] and row["case"] != "unmeasured":
        base += f", {row['case']} on the two axes"
    if row["attributed_elements"]:
        n = row["attributed_elements"]
        base += f", {n} attributed element{'s' if n != 1 else ''}"
    return base


def _declared(chrom: str, results_dir: Path, budget: dict[str, Any] | None) -> dict[str, Any]:
    """What `organise(chrom)` ACTUALLY read, measured, for the result's `inputs` block.

    Until 2026-10-02 this block was four f-strings built without looking at anything, and a
    declaration decoupled from the read is worse than one that shrinks, because it cannot even be
    wrong in a way the filesystem would show. Three things it got wrong on every committed
    `organised_<chrom>`, all three measured on chr21 before the change:

    * it named `budget_<chrom>` as "the budget", but `read_axes` prefers `budget_axes_<chrom>` when
      it is there and restates the tiers from it -- which it did on all 24 chromosomes, because both
      files are present for all 24. The record named the file that did NOT drive the join, and
      `axes_source`, which says which did, was thrown away. That is the shape where a declaration
      does not merely report a smaller run but names the WRONG FILE;
    * it named two of the three attribution runs in `RUNS` and omitted `enhancer_targets_all`, which
      is the largest by far (12,139 of the 12,439 elements joined on chr21) and the only one whose
      summary is a POINTER at a local table of about 1.4 GB. The biggest input, and the one that
      cannot be recovered from the repository, was the one left out;
    * it named `variation_<chrom>` and `duplication_<chrom>` unconditionally, whether or not they
      were there. `blocks()` reads both through `or {}`, so on a chromosome without them every
      block's human-axis and copy field is null and the result would still have claimed them. The
      mirror of a shrinking declaration: one that will not shrink when it should.

    So: the budget branch that ran is named and `budget_branch` records which it was; every run in
    `RUNS` is asked for; a conditional input that is absent is recorded under
    `inputs_optional_absent` rather than dropped or claimed; and `declared` carries the path list
    from `inputs()`, which follows each pointer and raises on one that does not resolve."""
    from genomeos.attribution.targets import RUNS

    declared = inputs(chrom, results_dir)  # raises if a pointer does not resolve
    present = {p.name for p in declared}
    absent: dict[str, str] = {}

    def named(stem: str, note: str) -> str | None:
        if f"{stem}.json" in present:
            return stem
        absent[stem] = note
        return None

    axes = f"budget_axes_{chrom}"
    branch = (budget or {}).get("axes_source")
    return {
        "budget": axes if branch == axes else f"budget_{chrom}",
        "budget_branch": branch,
        "budget_pair_considered": [axes, f"budget_{chrom}"],
        "variation": named(
            f"variation_{chrom}",
            "not on this machine: blocks() reads it through `or {}`, so every block's human_fraction "
            "and case is null and no block carries a human axis",
        ),
        "duplication": named(
            f"duplication_{chrom}",
            "not on this machine: blocks() reads it through `or {}`, so every block's "
            "duplicated_fraction is null and no block is read as a copy",
        ),
        "elements": [
            n
            for r in RUNS
            if (n := named(f"{r}_{chrom}", "not on this machine: its elements are not in the join"))
        ],
        "declared": [str(p) for p in declared],
        "inputs_optional_absent": absent,
    }


def organise(chrom: str, results_dir: Path = RESULTS_DIR, top: int = 25) -> dict[str, Any]:
    t0 = time.time()
    budget = read_axes(chrom, results_dir)
    rows = blocks(chrom, results_dir, budget)
    by_tier: dict[str, dict[str, Any]] = {}
    for t in TIERS:
        sel = [r for r in rows if r["tier"] == t]
        copies = [r for r in sel if r["copy"]]
        cases: dict[str, dict[str, int]] = {}
        for r in sel:
            if r["copy"]:
                continue
            c = r["case"] or "unmeasured"
            cases.setdefault(c, {"blocks": 0, "bp": 0})
            cases[c]["blocks"] += 1
            cases[c]["bp"] += r["length"]
        by_tier[t] = {
            "blocks": len(sel),
            "bp": sum(r["length"] for r in sel),
            "copies": len(copies),
            "copies_bp": sum(r["length"] for r in copies),
            "after_copies": len(sel) - len(copies),
            "after_copies_bp": sum(r["length"] for r in sel if not r["copy"]),
            "cases_after_copies": cases,
        }
    cu = [r for r in rows if r["tier"] == "constrained_unknown" and not r["copy"]]
    rank = {"syntax": 0, "relaxed": 1, "recent": 2, "unmeasured": 3, "tolerant": 4, None: 3}
    cu.sort(key=lambda r: (rank.get(r["case"], 3), -(r["mammal_fraction"] or 0)))
    copies = sorted((r for r in rows if r["copy"]), key=lambda r: -r["length"])
    return {
        "chrom": chrom,
        "blocks": len(rows),
        "with_human_axis": sum(1 for r in rows if r["case"] and r["case"] != "unmeasured"),
        "with_copy_flag": sum(1 for r in rows if r["duplicated_fraction"] is not None),
        "copies": len(copies),
        "copies_bp": sum(r["length"] for r in copies),
        "copy_min_fraction": COPY_MIN,
        "by_tier": by_tier,
        "real_unknown": {
            "blocks": len(cu),
            "bp": sum(r["length"] for r in cu),
            "by_case": by_tier["constrained_unknown"]["cases_after_copies"],
        },
        "candidates": cu[:top],
        "largest_copies": copies[:top],
        "inputs": _declared(chrom, results_dir, budget),
        "evidence": {
            "tier": "inferred: sequence class plus Zoonomia constraint (the budget)",
            "case": "inferred: Zoonomia phyloP against gnomAD Gnocchi (area J); unmeasured below a kilobase",
            "copy": "curated: UCSC genomicSuperDups, >= 1 kb at >= 90% identity; copy at >= 50% of the block",
        },
        "cost": {"seconds": round(time.time() - t0, 2)},
    }


def distil(results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    per: dict[str, dict[str, Any]] = {}
    by_tier: dict[str, dict[str, Any]] = {
        t: {"blocks": 0, "bp": 0, "copies": 0, "copies_bp": 0, "cases_after_copies": {}} for t in TIERS
    }
    real = {"blocks": 0, "bp": 0, "by_case": {}}
    for c in chroms:
        r = load_result(f"organised_{c}", results_dir)
        if not r:
            continue
        per[c] = {
            "blocks": r["blocks"],
            "copies": r["copies"],
            "real_unknown_blocks": r["real_unknown"]["blocks"],
            "real_unknown_bp": r["real_unknown"]["bp"],
            "syntax_candidates": (r["real_unknown"]["by_case"].get("syntax") or {}).get("blocks", 0),
        }
        for t, v in r["by_tier"].items():
            d = by_tier[t]
            for k in ("blocks", "bp", "copies", "copies_bp"):
                d[k] += v[k]
            for case, cv in v["cases_after_copies"].items():
                e = d["cases_after_copies"].setdefault(case, {"blocks": 0, "bp": 0})
                e["blocks"] += cv["blocks"]
                e["bp"] += cv["bp"]
        real["blocks"] += r["real_unknown"]["blocks"]
        real["bp"] += r["real_unknown"]["bp"]
        for case, cv in r["real_unknown"]["by_case"].items():
            e = real["by_case"].setdefault(case, {"blocks": 0, "bp": 0})
            e["blocks"] += cv["blocks"]
            e["bp"] += cv["bp"]
    return {
        "chromosomes": len(per),
        "per_chromosome": per,
        "by_tier": by_tier,
        "real_unknown": real,
        "copy_min_fraction": COPY_MIN,
        "note": (
            "The budget's constrained_unknown tier re-read with copies set apart and the human axis "
            "attached: the real unknown is what remains, split by the case the two axes make."
        ),
    }


def run_and_save(chrom: str, results_dir: Path = RESULTS_DIR, **kw) -> dict[str, Any]:
    out = organise(chrom, results_dir, **kw)
    save_result(f"organised_{chrom}", out, results_dir)
    return out
