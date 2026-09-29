# SPDX-License-Identifier: AGPL-3.0-or-later
"""Review item R7, step 1: every single-label class or role the compiled programs carry, and how often.

    uv run python scripts/ontology_census.py [--no-save]

The external review of 2026-09-28 found that `attribution/compile.py` writes `class: enhancer` on
every attributed element, even when its rule says the element represses its gene, and that UNKNOWN
regions get one broad role string that joins sequence origin, a guessed function and a reading of
constraint. This script changes nothing. It lists where each such label is assigned and read (by hand,
checked with `git grep`), counts the values over the 24 compiled programs in data/knowledge/compiled,
and counts, from the local ENCODE registry, RepeatMasker and `genomeos unknown` files, what the same
elements and regions would carry if origin, role, activity, target relation and evidence status were
kept apart. No requests, no model.
"""

from __future__ import annotations

import argparse
import bisect
import re
from collections import Counter
from pathlib import Path
from typing import Any

from genomeos import manifest as mf
from genomeos.genome.regulatory import load_ccres
from genomeos.genome.repeats import INTERSPERSED, load_repeats
from genomeos.results import RESULTS_DIR, load_result, save_result

CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]
COMPILED = Path("data/knowledge/compiled")

# Read by hand at 3a1fc0e and checked against `git grep` for every consumer.
ASSIGNMENTS: list[dict[str, str]] = [
    {
        "where": "genomeos/attribution/compile.py:390 (compile_chromosome, attributed elements)",
        "label": "element `class:`",
        "assigns": "the constant `enhancer` for every predicted element, whatever the rule's action",
    },
    {
        "where": "genomeos/attribution/compile.py:229 (_measured_blocks)",
        "label": "element `class:` of `<id>_measured`",
        "assigns": "row['class'] or `enhancer`; no measured row carries a class, so always `enhancer`",
    },
    {
        "where": "genomeos/attribution/compile.py:149-152 (_region)",
        "label": "region `role:`",
        "assigns": "`<tier>, <budget label>`, prefixed `copy, ` for a duplicated block; `unknown` for "
        "the constrained_unknown tier (that is what Module.unknowns counts)",
    },
    {
        "where": "genomeos/attribution/budget.py:81-155 (_rule)",
        "label": "budget tier and label",
        "assigns": "one tier from the block's one sequence class plus constraint; the label reads "
        "constraint as function or its absence as none (see constraint_readings)",
    },
    {
        "where": "genomeos/genome/unknown.py:272-383 (classify_block)",
        "label": "block sequence class",
        "assigns": "first matching rule wins: a block dense in ENCODE cCREs is `regulatory` before its "
        "RepeatMasker coverage is read, so a repeat-derived regulatory block keeps only `regulatory`",
    },
    {
        "where": "genomeos/attribution/organise.py:125-140 (reading)",
        "label": "block reading",
        "assigns": "copy first, else the tier; reads the budget, assigns no new class",
    },
    {
        "where": "genomeos/genome/regulation.py:34-41 (CLS)",
        "label": "RegulatoryElement.cls from ENCODE",
        "assigns": "PLS promoter, pELS/dELS enhancer, CTCF-only insulator, DNase-H3K4me3 open_chromatin; "
        "a different path from compile.py, not used by the compiled programs",
    },
    {
        "where": "genomeos/lang/parser.py:477-500, genomeos/lang/grammar.py:56-69",
        "label": "grammar",
        "assigns": "region `role: text | unknown` (free text); element `class: promoter | enhancer | "
        "insulator | open_chromatin` (one value, not checked); no key for origin, activity or status",
    },
]

CONSUMERS: list[dict[str, str]] = [
    {
        "what": "RegulatoryElement.cls (parsed `class:`)",
        "readers": "genome/regulation.py:136-153, genome/regdiff.py:74-81, predict/enhancer_target.py:352 "
        "(all on elements built from ENCODE by regulation.py, none on a compiled program); "
        "tests/test_regulation.py; test fixtures in tests/test_grn_bridge.py and "
        "tests/test_measured_layer.py write `class: enhancer` by hand",
    },
    {
        "what": "Region.role (parsed `role:`)",
        "readers": "ir/model.py Module.unknowns (role is UNKNOWN), the `# test: unknowns ==` line of every "
        "compiled program, ir/model.py from_dict",
    },
    {
        "what": "the runtime",
        "readers": "runtime/grn.py reads rules, strength and `when`; it reads neither `class:` nor `role:`",
    },
    {
        "what": "evidence explorer and web",
        "readers": "evidence.py and web/server.py read evidence kinds and confidences, not `class:`; "
        "web/static/blocks.js and index.html show the unknown classes and cCRE classes from results, "
        "not from programs",
    },
    {
        "what": "stored results",
        "readers": "budget_<chrom>.json (tier, label), organise_<chrom>.json (tier), "
        "the committed data/organisms/human/noncoding_chr21.bio",
    },
]

# Budget labels that read constraint as a mechanism, or its absence as no function.
CONSTRAINT_READINGS = {
    "transposable-element fossil, unconstrained": "absence of constraint read as a dead fossil",
    "unconstrained unique sequence, no evidence of function, best guess neutral": "absence read as neutral",
    "long open reading frame without constraint, a chance frame or a dead pseudogene": "absence read as dead",
    "weakly constrained unique sequence, mostly neutral": "weak constraint read as neutral",
    "regulatory elements by the registry, little constraint, lineage-specific or weak": (
        "little constraint read as weak function"
    ),
    "repeat-derived sequence under constraint, possibly exapted as a regulatory element": (
        "constraint read as a regulatory mechanism (hedged)"
    ),
    "long open reading frame under constraint, unannotated coding or a young pseudogene": (
        "constraint read as coding (two alternatives)"
    ),
}

ELEMENT = re.compile(r"^element (\S+) \{")
RULE = re.compile(r"^rule (\S+) (activates|inhibits) ")


def program_census(path: Path) -> dict[str, Any]:
    cls: Counter = Counter()
    heads: Counter = Counter()
    roles: Counter = Counter()
    actions: Counter = Counter()
    el_class: dict[str, str] = {}
    loci: dict[str, tuple[int, int]] = {}
    cur = None
    for line in path.open():
        m = ELEMENT.match(line)
        if m:
            cur = m.group(1)
            continue
        if cur and line.startswith("  class: "):
            el_class[cur] = line[9:].strip()
            cls[("measured" if cur.endswith("_measured") else "predicted", el_class[cur])] += 1
        elif cur and line.startswith("  locus: "):
            a, b = line.split(":")[2].split("-")
            loci[cur] = (int(a), int(b))
        elif line.startswith("}"):
            cur = None
        elif line.startswith("  role: "):
            r = line[8:].strip()
            roles[r] += 1
            heads[r.split(",")[0]] += 1
        else:
            m = RULE.match(line)
            if m:
                kind = "measured" if m.group(1).endswith("_measured") else "predicted"
                actions[(kind, el_class.get(m.group(1), "?"), m.group(2))] += 1
    return {"class": cls, "heads": heads, "roles": roles, "actions": actions, "loci": loci}


def interspersed(repeats: list, starts: list[int], max_len: int, s: int, e: int) -> tuple[float, str]:
    cover: Counter = Counter()
    for r in repeats[bisect.bisect_left(starts, s - max_len) :]:
        if r.start >= e:
            break
        ov = min(r.end, e) - max(r.start, s)
        if ov > 0 and r.cls in INTERSPERSED:
            cover[r.cls] += ov
    if not cover:
        return 0.0, ""
    top, _ = cover.most_common(1)[0]
    return min(1.0, sum(cover.values()) / max(1, e - s)), top


def element_census(chrom: str, loci: dict[str, tuple[int, int]]) -> dict[str, Counter]:
    reg = {c.id: c for c in load_ccres(chrom)}
    reps = load_repeats(chrom)
    starts = [r.start for r in reps]
    max_len = max((r.end - r.start for r in reps), default=0)
    ccre: Counter = Counter()
    origin: Counter = Counter()
    for eid, (s, e) in loci.items():
        if eid.endswith("_measured"):
            continue
        c = reg.get(eid)
        ccre[(c.cls if c else "not_in_registry") + (" CTCF-bound" if c and c.ctcf_bound else "")] += 1
        f, top = interspersed(reps, starts, max_len, s, e)
        origin[
            "repeat_derived (>=0.5 interspersed)" if f >= 0.5 else "partly (>0, <0.5)" if f > 0 else "none"
        ] += 1
        if f >= 0.5:
            origin[f"repeat_derived {top}"] += 1
    return {"ccre": ccre, "origin": origin}


def region_census(chrom: str) -> Counter:
    unk = {b["start"]: b for b in (load_result(f"unknown_{chrom}") or {}).get("blocks", [])}
    out: Counter = Counter()
    for b in (load_result(f"budget_{chrom}") or {}).get("blocks", []):
        f = (unk.get(b["start"]) or {}).get("features", {})
        ic = f.get("interspersed_coverage")
        state = "not read" if ic is None else ">=0.5" if ic >= 0.5 else ">0" if ic > 0 else "0"
        out[f"{b['class']} | interspersed {state}"] += 1
    return out


def inputs() -> list[dict[str, Any]]:
    out = []
    for ch in CHROMS:
        for p in (
            COMPILED / f"noncoding_{ch}.bio",
            RESULTS_DIR / f"ccres_{ch}.bed.gz",
            RESULTS_DIR / f"rmsk_{ch}.bed.gz",
            RESULTS_DIR / f"budget_{ch}.json",
            RESULTS_DIR / f"unknown_{ch}.json",
        ):
            if p.exists():
                out.append(mf.input_entry(p, partition=None))
    return out


def manifest() -> dict[str, Any]:
    return {
        "sources": [
            {"accession": "this repository, the 24 compiled programs", "version": "pinned by sha256"},
            {"accession": "ENCODE SCREEN cCREs v3", "version": "ccres_<chrom> as distilled"},
            {"accession": "UCSC rmsk track, hg38", "version": "rmsk_<chrom> as distilled"},
        ],
        "inputs": inputs(),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {"repeat_derived_min_interspersed_fraction": 0.5, "chromosomes": CHROMS},
        "exclusions": [],
        "partitions": "n/a: a census of labels, no evaluation",
    }


def _flat(c: Counter) -> dict[str, int]:
    return {(" | ".join(k) if isinstance(k, tuple) else k): v for k, v in c.most_common()}


AXIS_KEYS = ("origin", "molecular_role", "activity", "target_relation", "evidence_status")


def axis_census() -> dict[str, Any]:
    """After the build: what the programs state on each axis, per block kind, counted per value (a
    `|` group counted as one alternatives value), plus the checks the pre-registration named."""
    by: dict[str, Counter] = {}
    checks: Counter = Counter()
    for ch in CHROMS:
        p = COMPILED / f"noncoding_{ch}.bio"
        if not p.exists():
            continue
        kind = None
        block: dict[str, str] = {}
        for line in p.open():
            if line.startswith(("element ", "region ")) and line.rstrip().endswith("{"):
                name = line.split()[1]
                kind = (
                    "region"
                    if line.startswith("region")
                    else ("measured" if name.endswith("_measured") else "predicted")
                )
                block = {}
            elif kind and line.startswith("}"):
                checks[f"{kind} blocks"] += 1
                if all(k in block for k in AXIS_KEYS):
                    checks[f"{kind} blocks stating all five axes"] += 1
                role = block.get("molecular_role", "")
                if "silencer" in [v.strip() for v in role.split(",")]:
                    checks[f"{kind} with silencer as a label on its own"] += 1
                if "represses_target" in block.get("activity", ""):
                    checks[f"{kind} represses_target"] += 1
                    if "silencer|" in role:
                        checks[f"{kind} represses_target with the role left open"] += 1
                if kind == "region" and "repeat_derived" in block.get("origin", "") and "_like" in role:
                    checks["regions repeat-derived and with a registry role"] += 1
                if kind == "predicted" and "repeat_derived" in block.get("origin", "") and "_like" in role:
                    checks["predicted elements repeat-derived and with a registry role"] += 1
                kind = None
            elif kind:
                key, _, value = line.strip().partition(": ")
                if key in AXIS_KEYS or key == "class":
                    block[key] = value
                    for v in value.split(", "):
                        base = v if "|" in v else v.split("/")[0]
                        by.setdefault(f"{kind} {key}", Counter())[base] += 1
    return {"values": {k: dict(c.most_common()) for k, c in sorted(by.items())}, "checks": dict(checks)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-save", action="store_true")
    ap.add_argument("--after", action="store_true", help="count the five axes in rebuilt programs")
    args = ap.parse_args()
    if args.after:
        payload = {"review_item": "R7", "stage": "after the build", **axis_census()}
        for k, v in payload["checks"].items():
            print(f"{k}: {v}")
        if not args.no_save:
            print(f"saved {save_result('ontology_census_after', payload, manifest=manifest())}")
        return 0
    tot: dict[str, Counter] = {k: Counter() for k in ("class", "heads", "roles", "actions", "ccre", "origin")}
    regions: Counter = Counter()
    programs = 0
    for ch in CHROMS:
        p = COMPILED / f"noncoding_{ch}.bio"
        if not p.exists():
            continue
        programs += 1
        pc = program_census(p)
        for k in ("class", "heads", "roles", "actions"):
            tot[k].update(pc[k])
        ec = element_census(ch, pc["loci"])
        tot["ccre"].update(ec["ccre"])
        tot["origin"].update(ec["origin"])
        regions.update(region_census(ch))
        print(ch, "done", flush=True)
    readings = {
        label: {"reading": why, "regions": sum(n for r, n in tot["roles"].items() if label in r)}
        for label, why in CONSTRAINT_READINGS.items()
    }
    payload = {
        "review_item": "R7",
        "programs": programs,
        "assignments": ASSIGNMENTS,
        "consumers": CONSUMERS,
        "element_class": _flat(tot["class"]),
        "rules_by_element_class_and_action": _flat(tot["actions"]),
        "region_role_head": _flat(tot["heads"]),
        "region_role": _flat(tot["roles"]),
        "constraint_readings": readings,
        "predicted_elements_by_registry_class": _flat(tot["ccre"]),
        "predicted_elements_by_repeat_origin": _flat(tot["origin"]),
        "regions_by_sequence_class_and_interspersed_coverage": _flat(regions),
        "reading": (
            "every compiled element says `class: enhancer`, including every element whose rule inhibits; "
            "the class is a constant, not a reading of the registry, the prediction or the sequence. Region "
            "roles join origin, function and a constraint reading in one string, and several of them read "
            "the absence of constraint as the absence of function"
        ),
    }
    for k in ("element_class", "rules_by_element_class_and_action", "region_role_head"):
        print(k, payload[k])
    print("registry", payload["predicted_elements_by_registry_class"])
    print("origin", payload["predicted_elements_by_repeat_origin"])
    if not args.no_save:
        print(f"saved {save_result('ontology_census', payload, manifest=manifest())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
