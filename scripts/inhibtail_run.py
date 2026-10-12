# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-inhibtail: is the compiled layer's above-bound inhibitor tail a property of the PREDICTION
INSTRUMENT? TWO STAGES, and the order between them is the whole point.

    uv run python scripts/inhibtail_run.py --stage baserate   # -> data/results/inhibtail_registration.json
    uv run python scripts/inhibtail_run.py --stage tail       # -> data/results/inhibtail_result.json

STAGE `baserate` RUNS AND IS COMMITTED FIRST, BEFORE ANY TAIL NUMBER IS READ. It conditions on nothing
but RHO > 1: it emits the evidence-tier composition of the inhibiting elements of the eligible
multi-element genes and of the non-inhibiting elements of the same genes, the rule-line counts by kind,
and the control-arm counts -- and NO strength, NO fold threshold, NO above-bound count and NO per-element
row. `_refuse_a_tail_number` is the check that holds that promise: every key the stage writes is matched
against `BASE_RATE_KEYS` and the stage REFUSES to save if it has produced anything else.

STAGE `tail` runs only after the registration is committed. It reads the registration back for the base
rate and the dead-field verdicts, computes the registered axes A1-A8, selects the outcome by
`inhibtail.select_outcome` (OUTCOME_RULE as code) and writes the result.

0 MODEL REQUESTS, NO KEY READ, NO NETWORK, $0.00. NO RNA VALUE IS OPENED and nothing about RNA is
claimed. NOTHING IS FITTED and NO FIT CHANGE IS PROPOSED. The stop list of 2026-09-27 is read and none
of its six items is touched (the check is at the head of genomeos/attribution/inhibtail.py).
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import bridge  # noqa: E402
from genomeos.attribution import combine2 as c2  # noqa: E402
from genomeos.attribution import compile as cp  # noqa: E402
from genomeos.attribution import inhibtail as it  # noqa: E402
from genomeos.attribution import readergate as rg  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.results import save_result  # noqa: E402

REGISTRATION = "inhibtail_registration"
RESULT = "inhibtail_result"
REGISTRATION_PATH = Path("data/results/inhibtail_registration.json")
CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]
PROGRAM_GROUP = "the 24 compiled non-coding programs"
OWN_CODE = (
    "genomeos/attribution/inhibtail.py",
    "scripts/inhibtail_run.py",
    "tests/test_inhibtail.py",
)

#: The ONLY keys stage `baserate` may write beside the registration's own frozen text. Anything else is
#: a tail number leaking into the registration and the stage refuses rather than saves.
BASE_RATE_KEYS = frozenset(
    {
        "base_rate",
        "dead_fields",
        "A9_control_arm",
        "populations_measured",
        "programs",
        "group_declaration",
        "pinned",
        "seconds",
    }
)

#: Key substrings forbidden inside stage `baserate`'s COMPUTED blocks, because each names a quantity
#: that could only come from the bound.
#: "rho" alone is NOT forbidden: `how_rho_was_read` is a base-rate field name and names no threshold.
#: What is forbidden is a key that could only exist if a bound had been applied.
FORBIDDEN_IN_BASERATE = (
    "tail",
    "the_bound",
    "above_bound",
    "strength",
    "offend",
    "rho_above",
    "log2_fold",
    "gap",
    "regime",
    "quartile",
    "enrich",
)


def program_paths(chroms: list[str]) -> list[Path]:
    """The declared input group, with the member count this lane DECLARED in `inhibtail.INPUT_GROUPS`.

    The count is asserted against the declaration and NOT discovered: no `.exists()` filter stands
    between the declaration and the run, so a missing program makes the run REFUSE instead of silently
    shrinking what it claims to have read.
    """
    want = it.INPUT_GROUPS[PROGRAM_GROUP]
    paths = [Path(COMPILED_DIR) / f"noncoding_{c}.bio" for c in chroms]
    if len(paths) != want:
        raise SystemExit(
            f"REFUSED: {PROGRAM_GROUP} is declared as {want} members in inhibtail.INPUT_GROUPS and this"
            f" run was given {len(paths)}. The declaration is the commit's, not the filesystem's."
        )
    missing = [str(p) for p in paths if not p.exists()]
    if missing:
        raise SystemExit(
            f"REFUSED: {PROGRAM_GROUP} declares {want} members and {len(missing)} are absent:"
            f" {missing}. A declared group is not filtered by existence."
        )
    return paths


def pinned_check() -> dict[str, Any]:
    """Every inherited figure read back FROM THE OBJECT IT WAS READ FROM, against `inhibtail.PINNED`."""
    cache: dict[str, Any] = {}
    rows, bad = [], []
    for label, path, keys, expected in it.PINNED:
        if path not in cache:
            cache[path] = json.loads(Path(path).read_text())
        node: Any = cache[path]
        for k in keys:
            node = node[k]
        ok = node == expected
        rows.append(
            {
                "label": label,
                "source": path,
                "json_path": "/".join(keys),
                "declared": expected,
                "read": node,
                "agrees": ok,
            }
        )
        if not ok:
            bad.append(f"{label}: declared {expected!r}, {path} {'/'.join(keys)} holds {node!r}")
    return {
        "note": (
            "Each inherited figure is pinned AS THE SAME OBJECT it was read from: the file, the json"
            " path and the value this lane declared, read back and compared. A disagreement is a"
            " refusal, not a warning, so this lane's summary cannot drift from the committed record."
        ),
        "rows": rows,
        "all_agree": not bad,
        "disagreements": bad,
    }


def _element_rows(chrom: str, genes: set[str]) -> dict[str, Any]:
    """One program: every element of every eligible multi-element gene, with its provenance.

    NOTHING here conditions on a bound, a strength or a fold: the same function serves both stages and
    the stages differ only in what they then compute. The rows carry RHO because RHO is what the gate
    read; the BASE-RATE STAGE NEVER LOOKS AT IT except through `rho > 1`, which is the definition of
    `inhibiting` and of nothing else.
    """
    from genomeos.lang.parser import parse

    text = (Path(COMPILED_DIR) / f"noncoding_{chrom}.bio").read_text()
    gated, info = cp.gate_program(text, chrom, selection=c2.SELECTION)
    module = parse(gated, f"noncoding_{chrom}_reader")
    prov: dict[tuple[str, str], tuple[str, str, str]] = {}
    for r in module.rules:
        if r.action in c2.REGULATORY and r.target in genes:
            ob = bridge.observation_of(r, module)
            how = ob.how if isinstance(ob, bridge.Observation) else f"unresolved: {ob.reason}"
            prov[(r.target, r.through)] = (bridge._kind(r), how, str(r.evidence.source))
    cells = rg.compiled_cells(text)
    out: dict[str, Any] = {
        "chrom": chrom,
        "program": {k: info[k] for k in ("base_sha256", "gated_sha256", "elements", "rules_gated")},
        "eligible_genes": 0,
        "rows": [],
        "unreadable_elements": 0,
    }
    for gr in c2.gene_rules(module, genes, chrom, cells):
        _fits, why = c2.classify(gr)
        if why:
            continue
        out["eligible_genes"] += 1
        n_inhibiting = sum(1 for r in gr.rho.values() if r > 1.0)
        for el in gr.elements:
            kind, how, source = prov.get((gr.gene, el), ("", "", ""))
            value = gr.rho.get(el)
            if value is None:
                out["unreadable_elements"] += 1
            out["rows"].append(
                {
                    "symbol": gr.symbol,
                    "chrom": chrom,
                    "element": el,
                    "rho": value,
                    "kind": kind,
                    "how": how,
                    "source_prefix": source.split(",")[0],
                    "cell": gr.compiled.get(el, ""),
                    "gene_n_inhibiting": n_inhibiting,
                }
            )
    return out


def _walk(chroms: list[str]) -> tuple[list[dict[str, Any]], dict[str, Any], int]:
    """The one traversal both stages use: POPULATION (a) then every element of every eligible gene."""
    texts = {c: (Path(COMPILED_DIR) / f"noncoding_{c}.bio").read_text() for c in chroms}
    chrom_of, why = c2.multi_rule_genes(texts)
    del texts
    by_chrom: dict[str, set[str]] = {c: set() for c in chroms}
    for g, c in chrom_of.items():
        if c in by_chrom:
            by_chrom[c].add(g)
    rows: list[dict[str, Any]] = []
    programs: dict[str, Any] = {}
    eligible = 0
    for c in chroms:
        part = _element_rows(c, by_chrom[c])
        programs[c] = part["program"]
        programs[c]["eligible_genes"] = part["eligible_genes"]
        programs[c]["unreadable_elements"] = part["unreadable_elements"]
        rows.extend(part["rows"])
        eligible += part["eligible_genes"]
        print(f"  {c}: {part['eligible_genes']} eligible genes, {len(part['rows'])} elements", flush=True)
    programs["multi_rule_exclusions"] = dict(why)
    return rows, programs, eligible


# ---------------------------------------------------------------------------------------------------
# STAGE 1: the base rate. Conditions on nothing but RHO > 1.
# ---------------------------------------------------------------------------------------------------


def base_rate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """The evidence-tier composition of the inhibiting elements and of the non-inhibiting elements.

    `inhibiting` is `rho > 1`. That is the only comparison in this function, and it is the definition
    of the population rather than a threshold on the fit.
    """
    out: dict[str, Any] = {}
    for which, keep in (
        ("inhibiting", lambda r: r["rho"] is not None and r["rho"] > 1.0),
        ("non_inhibiting", lambda r: r["rho"] is not None and r["rho"] <= 1.0),
        ("unreadable", lambda r: r["rho"] is None),
    ):
        sub = [r for r in rows if keep(r)]
        out[which] = {
            "n": len(sub),
            "evidence_kind": dict(Counter(r["kind"] for r in sub)),
            "how_rho_was_read": dict(Counter(r["how"] for r in sub)),
            "evidence_source_prefix": dict(Counter(r["source_prefix"] for r in sub)),
        }
    out["note"] = it.BASE_RATE_DECLARATION
    out["population_clause"] = it.POPULATIONS["inhibiting_elements"]
    return out


def control_arm(rows: list[dict[str, Any]], programs: dict[str, Any]) -> dict[str, Any]:
    """AXES['A9_control_arm'], whose expected value 0 and its reason are declared in the registration."""
    gated: Counter = Counter()
    for _c, p in programs.items():
        if isinstance(p, dict) and "rules_gated" in p:
            gated.update(p["rules_gated"])
    inhibiting = [r for r in rows if r["rho"] is not None and r["rho"] > 1.0]
    not_predicted = [r for r in inhibiting if r["kind"] != "predicted"]
    by_kind = Counter(r["kind"] for r in not_predicted)
    which = "no_control_arm" if not not_predicted else "control_arm"
    if which == "no_control_arm":
        text = it.reading(
            "A9_control_arm",
            which,
            rules=sum(gated.values()),
            experimental=gated.get("experimental", 0),
            inhibiting=len(inhibiting),
        )
    else:
        text = it.reading(
            "A9_control_arm",
            which,
            n=len(not_predicted),
            inhibiting=len(inhibiting),
            by_kind=dict(by_kind),
        )
    return {
        "question": it.AXES["A9_control_arm"]["question"],
        "rule": it.AXES["A9_control_arm"]["rule"],
        "rule_lines_gated_by_kind": dict(gated),
        "elements_of_eligible_genes": len(rows),
        "inhibiting_elements": len(inhibiting),
        "inhibiting_elements_not_predicted": len(not_predicted),
        "not_predicted_by_kind": dict(by_kind),
        "expected_value_was": 0,
        "expected_value_held": not not_predicted,
        "decision": which,
        "reading": text,
        "decides_the_outcome": False,
    }


def _refuse_a_tail_number(payload: dict[str, Any]) -> None:
    """The registration stage's own guard: it must not have produced a tail number.

    Two checks. (a) Every computed key is in BASE_RATE_KEYS. (b) No KEY anywhere inside those computed
    blocks matches a FORBIDDEN_IN_BASERATE substring; the registration's frozen TEXT may name those
    quantities, because text is not a measurement, and the walk looks only at keys.
    """
    computed = set(payload["computed_keys"])
    extra = sorted(computed - BASE_RATE_KEYS)
    if extra:
        raise SystemExit(f"REFUSED: the base-rate stage computed keys it did not declare: {extra}")
    bad: list[str] = []

    def walk(node: Any, path: str) -> None:
        if isinstance(node, dict):
            for k, v in node.items():
                low = str(k).lower()
                if any(f in low for f in FORBIDDEN_IN_BASERATE):
                    bad.append(f"{path}/{k}")
                walk(v, f"{path}/{k}")
        elif isinstance(node, list):
            for i, v in enumerate(node):
                walk(v, f"{path}/{i}")

    for key in sorted(computed):
        walk(payload.get(key), key)
    if bad:
        raise SystemExit(
            "REFUSED: the base-rate stage's computed blocks carry keys naming a bound-conditioned"
            f" quantity, so a tail number would have entered the registration: {bad}"
        )


def registration(chroms: list[str]) -> dict[str, Any]:
    t0 = time.time()
    paths = program_paths(chroms)
    rows, programs, eligible = _walk(chroms)
    br = base_rate(rows)
    payload: dict[str, Any] = {
        "result": "the registration",
        "date": it.DATE,
        "lane": it.LANE,
        "area": it.AREA,
        "question": it.QUESTION,
        "populations": it.POPULATIONS,
        "disclosure": it.DISCLOSURE,
        "dead_by_definition": it.DEAD_BY_DEFINITION,
        "base_rate_declaration": it.BASE_RATE_DECLARATION,
        "base_rate_fields": list(it.BASE_RATE_FIELDS),
        "enrichment_statement": it.ENRICHMENT_STATEMENT,
        "composition_axis_dead_rule": it.COMPOSITION_AXIS_DEAD_RULE,
        "unanswerable_answer": it.UNANSWERABLE_ANSWER,
        "axes": it.AXES,
        "not_computed": it.NOT_COMPUTED,
        "outcomes": it.OUTCOMES,
        "outcome_rule": it.OUTCOME_RULE,
        "conventions": it.CONVENTIONS,
        "inherited_readings": it.INHERITED_READINGS,
        "declared_readers": list(it.DECLARED_READERS),
        "input_groups": it.INPUT_GROUPS,
        "limitations": list(it.LIMITATIONS),
        "rna_read": "none",
        "alphagenome_requests": 0,
        "money": "none: no request is sent, no key is read and no network is used",
        # ---- everything below this line is COMPUTED by this stage ------------------------------------
        "computed_keys": sorted(BASE_RATE_KEYS),
        "base_rate": br,
        "dead_fields": it.dead_fields(br),
        "A9_control_arm": control_arm(rows, programs),
        "populations_measured": {
            "eligible_multi_element_genes": eligible,
            "elements_of_those_genes": len(rows),
            "inhibiting_elements": br["inhibiting"]["n"],
            "non_inhibiting_elements": br["non_inhibiting"]["n"],
            "elements_with_an_unreadable_observation": br["unreadable"]["n"],
        },
        "programs": programs,
        "group_declaration": {
            "label": PROGRAM_GROUP,
            "declared_members": it.INPUT_GROUPS[PROGRAM_GROUP],
            "members_read": len(paths),
            "note": (
                "The member count is this lane's own declaration in inhibtail.INPUT_GROUPS. It is"
                " asserted, never discovered: no .exists() filter stands between the declaration and"
                " the run, so a missing program refuses the run instead of shrinking what it claims."
            ),
        },
        "pinned": pinned_check(),
        "seconds": round(time.time() - t0, 1),
    }
    _refuse_a_tail_number(payload)
    if not payload["pinned"]["all_agree"]:
        raise SystemExit(f"REFUSED: inherited figures disagree: {payload['pinned']['disagreements']}")
    return payload


# ---------------------------------------------------------------------------------------------------
# STAGE 2: the tail. Runs only after the registration is committed.
# ---------------------------------------------------------------------------------------------------


def _concentration(
    groups: dict[str, tuple[int, int]],
    pooled: float,
    divisor: int,
    total: int,
    axis: str,
    label_key: str,
) -> dict[str, Any]:
    """A4 and A5's shared rule: tested groups, exact two-sided binomial p, rate ratio, verdict."""
    tested = {g: (k, n) for g, (k, n) in groups.items() if n >= it.MIN_GROUP}
    per: dict[str, Any] = {}
    hits = []
    for g, (k, n) in sorted(tested.items()):
        rate = k / n
        pval = it.binom_two_sided(k, n, pooled)
        ratio = (rate / pooled) if pooled > 0 else None
        per[g] = {"k": k, "n": n, "rate": rate, "rate_ratio": ratio, "p_two_sided": pval}
        if pval < it.ALPHA / divisor and ratio is not None and ratio >= it.RATE_RATIO:
            hits.append((g, k, n, rate, ratio, pval))
    rates = [v["rate"] for v in per.values()]
    spread = (max(rates) / min(rates)) if rates and min(rates) > 0 else None
    if hits:
        g, k, n, rate, ratio, pval = max(hits, key=lambda h: h[4])
        kw: dict[str, Any] = {
            label_key: g,
            "k": k,
            "n": n,
            "rate": rate,
            "p": pooled,
            "total": total,
            "ratio": ratio,
            "pval": pval,
        }
        if axis == "A5_compiled_cell":
            kw["labels"] = divisor
        text = it.reading(axis, "concentrated", **kw)
        decision = "concentrated"
    else:
        kw = {
            "tested": len(tested),
            "lo": min(rates) if rates else None,
            "hi": max(rates) if rates else None,
            "p": pooled,
            "spread": spread,
        }
        if axis == "A5_compiled_cell":
            kw["labels"] = divisor
        text = it.reading(axis, "not_concentrated", **kw)
        decision = "not_concentrated"
    return {
        "question": it.AXES[axis]["question"],
        "rule": it.AXES[axis]["rule"],
        "pooled_rate": pooled,
        "groups_total": len(groups),
        "groups_tested": len(tested),
        "min_group": it.MIN_GROUP,
        "bonferroni_divisor": divisor,
        "threshold": it.ALPHA / divisor,
        "per_group": per,
        "groups_not_tested": sorted(set(groups) - set(tested)),
        "max_over_min_rate_ratio": spread,
        "hits": [
            {"group": h[0], "k": h[1], "n": h[2], "rate": h[3], "rate_ratio": h[4], "p_two_sided": h[5]}
            for h in hits
        ],
        "decision": decision,
        "reading": text,
        "decides_the_outcome": it.AXES[axis]["decides_the_outcome"],
    }


def tail(chroms: list[str]) -> dict[str, Any]:
    t0 = time.time()
    paths = program_paths(chroms)
    reg = json.loads(REGISTRATION_PATH.read_text())
    rows, programs, eligible = _walk(chroms)
    inh = [r for r in rows if r["rho"] is not None and r["rho"] > 1.0]
    rhos = [r["rho"] for r in inh]
    n = len(inh)
    strengths = [it.implied_strength(r) for r in rhos]
    scaled = [it.implied_strength(1.0 + it.SCALE * (r - 1.0)) for r in rhos]
    in_tail = [it.above_bound(s) for s in strengths]
    in_tail_scaled = [it.above_bound(s) for s in scaled]
    tail_n = sum(in_tail)
    body_n = n - tail_n
    tail_rhos = [r for r, t in zip(rhos, in_tail, strict=True) if t]
    body_rhos = [r for r, t in zip(rhos, in_tail, strict=True) if not t]

    # ---- A8: the tautology, verified as an identity ---------------------------------------------------
    mismatches = sum(
        1 for s, r in zip(strengths, rhos, strict=True) if it.above_bound(s) != (r > it.BOUND_RHO_UNSCALED)
    )
    a8_which = "tautological" if mismatches == 0 else "broken"
    a8 = {
        "question": it.AXES["A8_tautology"]["question"],
        "rule": it.AXES["A8_tautology"]["rule"],
        "mismatches": mismatches,
        "min_rho_in_the_tail": min(tail_rhos) if tail_rhos else None,
        "max_rho_in_the_body": max(body_rhos) if body_rhos else None,
        "bound_rho_unscaled": it.BOUND_RHO_UNSCALED,
        "bound_log2_fold_unscaled": math.log2(it.BOUND_RHO_UNSCALED),
        "bound_rho_at_the_scale": it.BOUND_RHO_AT_THE_SCALE,
        "bound_log2_fold_at_the_scale": math.log2(it.BOUND_RHO_AT_THE_SCALE),
        "decision": a8_which,
        "reading": it.reading("A8_tautology", a8_which, mismatches=mismatches, n=n),
        "decides_the_outcome": False,
    }

    # ---- A1: the composition, against the registration's committed base rate -------------------------
    dead = reg["dead_fields"]
    a1_fields: dict[str, Any] = {}
    any_live_excludes_1 = False
    for field, key in (
        ("evidence_kind", "kind"),
        ("how_rho_was_read", "how"),
        ("evidence_source_prefix", "source_prefix"),
    ):
        d = dead[field]
        value = d["most_common_value"]
        if d["dead"]:
            a1_fields[field] = {
                "dead": True,
                "base_rate_fraction": d["fraction"],
                "most_common_value": value,
                "decision": "dead",
                "reading": it.reading(
                    "A1_composition",
                    "dead",
                    field=field,
                    frac=d["fraction"],
                    n=d["n"],
                    value=value,
                    dead=it.DEAD_FRACTION,
                ),
            }
            continue
        tv = sum(1 for r, t in zip(inh, in_tail, strict=True) if t and r[key] == value)
        bv = sum(1 for r, t in zip(inh, in_tail, strict=True) if not t and r[key] == value)
        ci = it.katz_interval(tv, tail_n, bv, body_n)
        excl = bool(not ci["degenerate"] and ci["lo"] is not None and (ci["lo"] > 1.0 or ci["hi"] < 1.0))
        any_live_excludes_1 = any_live_excludes_1 or excl
        which = "enriched" if excl else "not_enriched"
        a1_fields[field] = {
            "dead": False,
            "base_rate_fraction": d["fraction"],
            "most_common_value": value,
            "tail_with_the_value": tv,
            "body_with_the_value": bv,
            "interval": ci,
            "interval_excludes_1": excl,
            "decision": which,
            "reading": it.reading(
                "A1_composition",
                which,
                field=field,
                value=value,
                rr=ci["ratio"],
                lo=ci["lo"],
                hi=ci["hi"],
                tail_n=tail_n,
                body_n=body_n,
            ),
        }
    a1 = {
        "question": it.AXES["A1_composition"]["question"],
        "rule": it.AXES["A1_composition"]["rule"],
        "base_rate_from_the_registration": {f: dead[f] for f in it.BASE_RATE_FIELDS},
        "composition_axis_is_dead": bool(dead["composition_axis_is_dead"]),
        "fields": a1_fields,
        "any_live_interval_excludes_1": any_live_excludes_1,
        "decision": "dead" if dead["composition_axis_is_dead"] else "measured",
        "decides_the_outcome": True,
    }

    # ---- A2: the gap at the bound --------------------------------------------------------------------
    whole = it.largest_multiplicative_gap(rhos)
    straddling = it.largest_multiplicative_gap(rhos, *it.GAP_STRADDLING_WINDOW)
    control = it.largest_multiplicative_gap(rhos, *it.GAP_CONTROL_WINDOW)
    at = whole["at"]
    near = at is not None and (
        it.BOUND_RHO_UNSCALED / it.GAP_NEAR_BAND <= at <= it.BOUND_RHO_UNSCALED * it.GAP_NEAR_BAND
    )
    w_ratio = (straddling["gap"] / control["gap"]) if straddling["gap"] and control["gap"] else None
    w_which = "windowed_gap" if (w_ratio is not None and w_ratio >= it.GAP_RATIO) else "windowed_no_gap"
    a2 = {
        "question": it.AXES["A2_gap_at_the_bound"]["question"],
        "rule": it.AXES["A2_gap_at_the_bound"]["rule"],
        "largest_gap_overall": whole,
        "straddling_window": {"window": list(it.GAP_STRADDLING_WINDOW), **straddling},
        "control_window": {"window": list(it.GAP_CONTROL_WINDOW), **control},
        "windowed_ratio": w_ratio,
        "decision": "separated" if near else "not_separated",
        "reading": it.reading(
            "A2_gap_at_the_bound",
            "separated" if near else "not_separated",
            n=n,
            gap=whole["gap"],
            at=whole["at"],
        ),
        "windowed_decision": w_which,
        "windowed_reading": it.reading(
            "A2_gap_at_the_bound",
            w_which,
            straddling=straddling["gap"],
            control=control["gap"],
            ratio=w_ratio,
        ),
        "decides_the_outcome": False,
    }

    # ---- A3: one regime or two -----------------------------------------------------------------------
    below = it.survival_slope(rhos, *it.REGIME_WINDOW_BELOW, include_hi=False)
    above = it.survival_slope(rhos, *it.REGIME_WINDOW_ABOVE, include_hi=True)
    ratio = (
        (above["slope"] / below["slope"])
        if below.get("slope") not in (None, 0.0) and above.get("slope") is not None
        else None
    )
    one = ratio is not None and abs(ratio - 1.0) <= it.REGIME_TOLERANCE
    a3_which = "one_regime" if one else "two_regimes"
    a3 = {
        "question": it.AXES["A3_one_regime_or_two"]["question"],
        "rule": it.AXES["A3_one_regime_or_two"]["rule"],
        "below": below,
        "above": above,
        "slope_ratio": ratio,
        "decision": a3_which,
        "reading": it.reading(
            "A3_one_regime_or_two",
            a3_which,
            below=below.get("slope"),
            above=above.get("slope"),
            ratio=ratio,
            n=n,
        ),
        "decides_the_outcome": False,
    }

    # ---- A4 and A5: the concentration axes -----------------------------------------------------------
    pooled = tail_n / n if n else 0.0
    by_chrom: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    by_cell: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for r, t in zip(inh, in_tail, strict=True):
        by_chrom[r["chrom"]][1] += 1
        by_cell[r["cell"]][1] += 1
        if t:
            by_chrom[r["chrom"]][0] += 1
            by_cell[r["cell"]][0] += 1
    chrom_groups = {k: (v[0], v[1]) for k, v in by_chrom.items()}
    cell_groups = {k: (v[0], v[1]) for k, v in by_cell.items()}
    a4 = _concentration(chrom_groups, pooled, it.N_PROGRAMS, n, "A4_chromosome", "chrom")
    tested_labels = sum(1 for _k, (_a, b) in cell_groups.items() if b >= it.MIN_GROUP)
    a5 = _concentration(cell_groups, pooled, max(tested_labels, 1), n, "A5_compiled_cell", "cell")

    # ---- A6: the gene's inhibitor count --------------------------------------------------------------
    labels = [r["gene_n_inhibiting"] for r in inh]
    q1, q3 = it.quantile(labels, 0.25), it.quantile(labels, 0.75)
    bottom = [(r, t) for r, t in zip(inh, in_tail, strict=True) if r["gene_n_inhibiting"] <= q1]
    top = [(r, t) for r, t in zip(inh, in_tail, strict=True) if r["gene_n_inhibiting"] > q3]
    b_rate = (sum(1 for _r, t in bottom if t) / len(bottom)) if bottom else None
    t_rate = (sum(1 for _r, t in top if t) / len(top)) if top else None
    a6_ratio = (t_rate / b_rate) if (b_rate and t_rate is not None) else None
    dep = a6_ratio is not None and (a6_ratio >= it.QUARTILE_RATIO or a6_ratio <= 1.0 / it.QUARTILE_RATIO)
    a6_which = "dependent" if dep else "independent"
    tail_labels = [r["gene_n_inhibiting"] for r, t in zip(inh, in_tail, strict=True) if t]
    body_labels = [r["gene_n_inhibiting"] for r, t in zip(inh, in_tail, strict=True) if not t]
    a6 = {
        "question": it.AXES["A6_gene_inhibitor_count"]["question"],
        "rule": it.AXES["A6_gene_inhibitor_count"]["rule"],
        "quartiles_of_the_label": {"q1": q1, "q3": q3},
        "bottom_quartile": {"elements": len(bottom), "tail_rate": b_rate},
        "top_quartile": {"elements": len(top), "tail_rate": t_rate},
        "rate_ratio_top_over_bottom": a6_ratio,
        "median_label_in_the_tail": it.quantile(tail_labels, 0.5) if tail_labels else None,
        "median_label_in_the_body": it.quantile(body_labels, 0.5) if body_labels else None,
        "decision": a6_which,
        "reading": it.reading(
            "A6_gene_inhibitor_count",
            a6_which,
            top=t_rate,
            top_n=len(top),
            bottom=b_rate,
            bottom_n=len(bottom),
            ratio=a6_ratio,
        ),
        "decides_the_outcome": False,
    }

    # ---- A7: sharing ---------------------------------------------------------------------------------
    genes_of: dict[str, set[str]] = defaultdict(set)
    for r in inh:
        genes_of[r["element"]].add(r["symbol"])
    shared = sum(1 for g in genes_of.values() if len(g) > 1)
    a7_which = "one_per_gene" if shared == 0 else "shared"
    a7 = {
        "question": it.AXES["A7_sharing"]["question"],
        "rule": it.AXES["A7_sharing"]["rule"],
        "distinct_element_ids": len(genes_of),
        "element_ids_under_more_than_one_gene": shared,
        "expected_value_was": 0,
        "expected_value_held": shared == 0,
        "max_genes_per_element": max((len(g) for g in genes_of.values()), default=0),
        "decision": a7_which,
        "reading": it.reading("A7_sharing", a7_which, elements=len(genes_of), shared=shared),
        "decides_the_outcome": False,
    }

    # ---- the outcome, by OUTCOME_RULE as code --------------------------------------------------------
    outcome = it.select_outcome(
        bool(dead["composition_axis_is_dead"]), a4["decision"], a5["decision"], any_live_excludes_1
    )

    return {
        "result": "the tail",
        "date": it.DATE,
        "lane": it.LANE,
        "area": it.AREA,
        "question": it.QUESTION,
        "registration": {
            "path": str(REGISTRATION_PATH),
            "note": (
                "The registration was written and COMMITTED by `--stage baserate` before this stage was"
                " run; the commit shas and the clock between them are in this lane's commit messages."
                " This stage reads the registration's `dead_fields` and does not recompute the base rate"
                " from its own traversal, so the base rate this result is read against is the committed"
                " one."
            ),
            "base_rate_is_dead": bool(dead["composition_axis_is_dead"]),
        },
        "populations_measured": {
            "eligible_multi_element_genes": eligible,
            "elements_of_those_genes": len(rows),
            "inhibiting_elements": n,
            "the_tail_unscaled": tail_n,
            "the_body_unscaled": body_n,
            "the_tail_at_the_scale": sum(in_tail_scaled),
        },
        "populations": it.POPULATIONS,
        "A1_composition": a1,
        "A2_gap_at_the_bound": a2,
        "A3_one_regime_or_two": a3,
        "A4_chromosome": a4,
        "A5_compiled_cell": a5,
        "A6_gene_inhibitor_count": a6,
        "A7_sharing": a7,
        "A8_tautology": a8,
        "A9_control_arm": reg["A9_control_arm"],
        "rho_distribution": {
            "n": n,
            "min": min(rhos) if rhos else None,
            "max": max(rhos) if rhos else None,
            "median": it.quantile(rhos, 0.5) if rhos else None,
            "q90": it.quantile(rhos, 0.9) if rhos else None,
            "q99": it.quantile(rhos, 0.99) if rhos else None,
            "the_tails_min": min(tail_rhos) if tail_rhos else None,
            "the_tails_median": it.quantile(tail_rhos, 0.5) if tail_rhos else None,
            "the_tails_max": max(tail_rhos) if tail_rhos else None,
            "log2_fold_of_the_tails_min": math.log2(min(tail_rhos)) if tail_rhos else None,
            "log2_fold_of_the_tails_max": math.log2(max(tail_rhos)) if tail_rhos else None,
            "population": it.POPULATIONS["inhibiting_elements"],
        },
        "outcome": outcome,
        "outcome_text": it.OUTCOMES[outcome],
        "outcome_rule": it.OUTCOME_RULE,
        "unanswerable_answer": it.UNANSWERABLE_ANSWER,
        "dead_by_definition": it.DEAD_BY_DEFINITION,
        "inherited_readings": it.INHERITED_READINGS,
        "not_computed": it.NOT_COMPUTED,
        "conventions": it.CONVENTIONS,
        "limitations": list(it.LIMITATIONS),
        "programs": programs,
        "group_declaration": {
            "label": PROGRAM_GROUP,
            "declared_members": it.INPUT_GROUPS[PROGRAM_GROUP],
            "members_read": len(paths),
        },
        "pinned": pinned_check(),
        "rna_read": "none",
        "alphagenome_requests": 0,
        "money": "none: no request is sent, no key is read and no network is used",
        "seconds": round(time.time() - t0, 1),
    }


def manifest(chroms: list[str], stage: str) -> dict[str, Any]:
    inputs = [
        mf.files_entry(PROGRAM_GROUP, program_paths(chroms)),
        mf.input_entry(it.CENSUS_RESULT, partition=None),
        mf.input_entry(it.EMIT7G_RESULT, partition=None),
    ]
    if stage == "tail":
        inputs.append(mf.input_entry(str(REGISTRATION_PATH), partition=None))
    return {
        "sources": [
            {
                "accession": "the 24 compiled non-coding programs, pinned by sha256 as one declared"
                " group of 24 members",
                "version": "pinned by the group digest",
            },
            {"accession": it.CENSUS_RESULT + " (lane-inhibitor, 85ca300)", "version": "dev"},
            {"accession": it.EMIT7G_RESULT + " (lane-emit7g, 589b6ca)", "version": "dev"},
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {"stage": stage, "selection": c2.SELECTION, **it.CONVENTIONS},
        "exclusions": [
            "combine2.multi_rule_genes' own: a gene with one rule, a gene with rules in several"
            " programs, and A GENE WITH ANY RULE THAT IS NOT PREDICTED -- the third is why the"
            " composition axis is dead by definition and it is this lane's finding, not a caveat",
            "combine2.classify's own: an unreadable or exactly-1 observation, an action that disagrees"
            " with it, a gene with no activating rule",
            "an element with RHO < 1 is an activator and is in the non-inhibiting population, never in"
            " the tail or the body",
            "A4 and A5 test only groups with at least MIN_GROUP = 1000 inhibiting elements; a"
            " concentration confined to a smaller group is invisible here",
        ],
        "partitions": "n/a: a predictor-side census; no outcome is read and no RNA is opened",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=("baserate", "tail"), required=True)
    ap.add_argument("--chroms", nargs="*", default=CHROMS)
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args()
    held_on = mf.trace_begin()
    if args.stage == "baserate":
        out, name = registration(args.chroms), REGISTRATION
    else:
        if not REGISTRATION_PATH.exists():
            raise SystemExit(f"REFUSED: {REGISTRATION_PATH} must be written AND COMMITTED first.")
        out, name = tail(args.chroms), RESULT
    man = manifest(args.chroms, args.stage)
    held = mf.trace_close()
    man["traced_inputs"] = mf.traced_inputs(man, held)
    man["traced_inputs"]["tracer_was_armed"] = bool(held_on)
    out["result_manifest"] = man
    if args.stage == "baserate":
        print("base rate, inhibiting:", out["base_rate"]["inhibiting"]["evidence_kind"], flush=True)
        print("base rate, non-inhibiting:", out["base_rate"]["non_inhibiting"]["evidence_kind"], flush=True)
        print("dead fields:", {k: out["dead_fields"][k]["dead"] for k in it.BASE_RATE_FIELDS}, flush=True)
        print("control arm:", out["A9_control_arm"]["decision"], flush=True)
    else:
        print("outcome:", out["outcome"], flush=True)
        for a in (
            "A1_composition",
            "A2_gap_at_the_bound",
            "A3_one_regime_or_two",
            "A4_chromosome",
            "A5_compiled_cell",
            "A6_gene_inhibitor_count",
            "A7_sharing",
            "A8_tautology",
        ):
            print(" ", a, "->", out[a]["decision"], flush=True)
    print("reconciliation:", mf.reconciliation_problems(man["traced_inputs"]) or "clean", flush=True)
    if args.no_save:
        print("--no-save: nothing written", flush=True)
        return 0
    path = save_result(name, out, manifest=man)
    print("written:", path, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
