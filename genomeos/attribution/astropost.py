# SPDX-License-Identifier: AGPL-3.0-or-later
"""AstroREG-2 after the run: the registered analyses, and nothing before the run is complete.

Three things, in the order they are applied.

1. THE GATE. Nothing here reads an answer unless the run's ledger shows the run COMPLETE: its last line
   is `run_complete`, written by the sender and by nothing else, and the run did not end at the pilot.
   Run 2 today holds 364 charges and 363 answers and no such line, so the gate refuses it and says so.
   The 363 answers decide nothing on their own (ROADMAP §5 item 14, finding 1), and an analysis that ran
   on them would be the outcome read before the run it belongs to had finished.

2. ITEM (h), the astrocyte-argmax question, exactly as `astroargmax_registration.json` (b5fd1b6) and its
   amendment 1 fix it: one label per element through `astroargmax.element_label`, the elements with no
   target counted beside the denominator, the label joined to the registered CURIE sets through the
   census, the rate against each set's registered control, a cluster bootstrap over whole chromosomes at
   the registered draws and seed, the identical-resample share for every interval, the degenerate route
   to an exact bound on an effective sample size, and the four registered readings in their own words.
   CNS_PRIMARY decides; the three alternatives are reported and may not be promoted.

3. THE PRIMARY CLAIM's registered tail -- the deletion gain, the locus and chromosome bootstraps, the
   wider quoted, the registered reading and the required quote -- over scored pairs. The ASSEMBLY of
   those pairs is REFUSED here, with each unregistered term named (`PRIMARY_GAPS`): the registrations
   do not fix them, and choosing them now, with the paid answers on disk, is a decision for a registered
   amendment and not for a runner.

Every number below is read from a committed registration at run time, never retyped, and the
registration's bytes are checked against the sha256 its amendment names. Nothing here sends a request,
reads a key, or writes a ledger.
"""

from __future__ import annotations

import collections
import hashlib
import json
import math
import random
from pathlib import Path
from typing import Any

from genomeos.attribution import astroargmax, astroreg, astrorun

#: The repository root. Registrations are read from here, not from the working directory.
REPO = Path(__file__).resolve().parents[2]

ASTROARGMAX_REGISTRATION = astrorun.ASTROARGMAX_REGISTRATION
ASTROARGMAX_AMENDMENT_1 = Path("data/results/astroargmax_registration_amendment_1.json")
ASTROREG2_REGISTRATION = Path("data/results/astroreg2_registration.json")
PRIMARY_REGISTRATION = Path("data/results/astroreg_registration.json")

#: The set whose reading is THE reading. The registration's own words: "CNS_PRIMARY ... decides the
#: reading. The three alternatives are reported ... none of them may be promoted to the reading".
DECIDING_SET = "CNS_PRIMARY"


class AnalysisRefusedError(RuntimeError):
    """A registered analysis may not run on what is on disk. The message says which condition failed."""


ONLY_A_COMPLETE_RUN_IS_ANALYSED = (
    "the registered post-run order is: the run completes, the supervisor reads every answer, the "
    "registered analyses run from their committed registrations, then the ledger is committed. A partial "
    "run's answers decide nothing on their own -- run 2's 363 cover 5 of 23 chromosome clusters, below "
    "the 10-cluster floor -- and reading them before the run they belong to is complete would put an "
    "outcome in view while the rest of the run could still be shaped by it. So no answer is opened until "
    "the ledger's last line is the run_complete line the sender writes"
)


# ----------------------------------------------------------------------------------------- the gate


def _ledger_rows(path: Path) -> list[dict[str, Any]]:
    rows = []
    for line in path.read_text(errors="replace").splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        rows.append(row if isinstance(row, dict) else {"event": None, "raw": row})
    return rows


def check_the_run_is_complete(ledger: Path | str) -> dict[str, Any]:
    """Refuse unless the ledger's LAST line is `run_complete` and the run did not stop at the pilot.

    Reads the ledger only. No answer and no cache entry is touched by this function, so a refusal is
    reached without anything of the outcome having been opened.
    """
    path = Path(ledger)
    if not path.exists():
        raise AnalysisRefusedError(
            f"there is no ledger at {path}, so no run is on record and nothing can be analysed. "
            f"{ONLY_A_COMPLETE_RUN_IS_ANALYSED}"
        )
    counts = astrorun.ledger_charges(path)
    if counts["unreadable_lines"]:
        raise AnalysisRefusedError(
            f"the ledger {path} holds {counts['unreadable_lines']} unreadable line(s), so whether the run "
            f"completed cannot be read from it. {ONLY_A_COMPLETE_RUN_IS_ANALYSED}"
        )
    rows = _ledger_rows(path)
    answered = astrorun.answered_elements(path)
    last = rows[-1] if rows else {}
    if not rows or last.get("event") != astrorun.RUN_COMPLETE_EVENT:
        where = (
            f"its last line is `{last.get('event')}`"
            + (f" n={last['n']}" if "n" in last else "")
            + (f" at {last['t']}" if "t" in last else "")
        )
        complete_earlier = astrorun.run_already_completed(path)
        raise AnalysisRefusedError(
            f"the ledger {path} does not show the run complete: "
            + (
                "a run_complete line exists but is not the last line, so the record continued past it"
                if complete_earlier
                else f"it has no `{astrorun.RUN_COMPLETE_EVENT}` line"
            )
            + f". It records {counts['charges']} charge(s) and {len(answered)} answer(s), and "
            f"{where}. {ONLY_A_COMPLETE_RUN_IS_ANALYSED}"
        )
    if last.get("stopped_at_pilot"):
        raise AnalysisRefusedError(
            f"the ledger {path} ends in a run_complete line that records a PILOT STOP. That run ended "
            "because its first answers were unusable; the registered analyses describe a run that bought "
            f"its list, and they are not run on the pilot's answers. {ONLY_A_COMPLETE_RUN_IS_ANALYSED}"
        )
    charged = astrorun.charged_elements(path)
    return {
        "ledger": str(path),
        "complete": True,
        "charges": counts["charges"],
        "answered": len(answered),
        "charged_without_an_answer": sorted(e for _, e in charged - answered),
        "run_complete_line": {k: v for k, v in last.items() if k != "why"},
    }


def load_the_answers(
    plan: list[dict[str, Any]], ledger: Path | str, cache_root: Path | str
) -> dict[str, dict[str, Any]]:
    """The answers of a COMPLETE run, one per answered plan element, after the gate.

    Refuses when the request root and the ledger disagree in either direction: a cached plan element the
    ledger never answered came from somewhere else, and an answered element with no entry means the root
    is not the one the run wrote to.
    """
    from genomeos.predict.enhancer_target import load_cached

    check_the_run_is_complete(ledger)
    answered = astrorun.answered_elements(ledger)
    planned = {(r["chrom"], r["element"]) for r in plan}
    stray = sorted(answered - planned)
    if stray:
        raise AnalysisRefusedError(
            f"the ledger answers {len(stray)} element(s) outside the reviewed list "
            f"(first: {[e for _, e in stray[:3]]}), so it is not this list's ledger"
        )
    out: dict[str, dict[str, Any]] = {}
    unanswered_but_cached, answered_but_absent = [], []
    for r in plan:
        key = (r["chrom"], r["element"])
        hit = load_cached(r["chrom"], r["element"], Path(cache_root))
        if key in answered:
            if hit is None:
                answered_but_absent.append(r["element"])
            else:
                out[r["element"]] = hit
        elif hit is not None:
            unanswered_but_cached.append(r["element"])
    if unanswered_but_cached or answered_but_absent:
        raise AnalysisRefusedError(
            f"the request root and the ledger disagree: {len(unanswered_but_cached)} cached element(s) "
            f"the ledger never answered (first: {unanswered_but_cached[:3]}) and "
            f"{len(answered_but_absent)} answered element(s) with no entry (first: "
            f"{answered_but_absent[:3]}). An analysis over either would not be an analysis of this run"
        )
    return out


# ----------------------------------------------------------------------------------- item (h)


def item_h_registration(root: Path | None = None) -> dict[str, Any]:
    """The committed item-(h) terms, read from the registration whose bytes amendment 1 names."""
    base = Path(root) if root is not None else REPO
    raw = (base / ASTROARGMAX_REGISTRATION).read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    amendment = json.loads((base / ASTROARGMAX_AMENDMENT_1).read_text())
    named = (amendment.get("amended_registration") or {}).get("blob_sha256")
    if digest != astroargmax.AMENDS_BLOB_SHA256 or digest != named:
        raise AnalysisRefusedError(
            f"{ASTROARGMAX_REGISTRATION} has sha256 {digest}; amendment 1 names {named} and the code "
            f"names {astroargmax.AMENDS_BLOB_SHA256}. A registration that is not the one registered "
            "registers nothing"
        )
    reg = json.loads(raw)
    bands = reg["bands"]
    for key, value in (amendment.get("bands_unchanged") or {}).items():
        if key in bands and bands[key] != value:
            raise AnalysisRefusedError(
                f"amendment 1 says {key} is unchanged at {value!r} and the registration holds "
                f"{bands[key]!r}. The two records disagree, so neither can be read as the band"
            )
    sets = {DECIDING_SET: [e["curie"] for e in reg["label_set_CNS_PRIMARY"]["included"]]}
    for name, block in reg["label_sets_registered_as_alternatives"].items():
        sets[name] = [e["curie"] for e in block["included"]]
    controls = {DECIDING_SET: float(reg["control"]["rate"])}
    for name, block in reg["control_for_each_alternative_set"].items():
        controls[name] = float(block["rate"])
    if set(sets) != set(controls):
        raise AnalysisRefusedError(f"label sets {sorted(sets)} and controls {sorted(controls)} differ")
    return {"bands": bands, "sets": sets, "controls": controls, "sha256": digest}


def wilson(rate: float, n: float, z: float = 1.959963984540054) -> list[float]:
    """The Wilson score interval for a proportion on (possibly fractional) `n`. The bound is NAMED."""
    if n <= 0:
        return [0.0, 1.0]
    denom = 1 + z * z / n
    centre = (rate + z * z / (2 * n)) / denom
    half = z * math.sqrt(rate * (1 - rate) / n + z * z / (4 * n * n)) / denom
    return [max(0.0, centre - half), min(1.0, centre + half)]


def reading_h(lo: float, hi: float, rate: float, bands: dict[str, Any]) -> dict[str, str]:
    """The registered reading of an interval on the ratio, in the registration's own words."""
    words = bands["readings"]
    eq_lo, eq_hi = bands["equivalence_interval_on_the_ratio"]
    if lo >= bands["ratio_enriched_floor"] and rate >= bands["absolute_level_floor"]:
        key = "enriched"
    elif hi <= bands["ratio_depleted_ceiling"]:
        key = "depleted"
    elif lo >= eq_lo and hi <= eq_hi:
        key = "equivalent"
    else:
        key = "inconclusive"
    return {"reading": key, "words": words[key]}


def rate_against_control(
    in_set: dict[str, bool], chrom_of: dict[str, str], control: float, bands: dict[str, Any]
) -> dict[str, Any]:
    """One label set's rate, its ratio to the control, the clustered interval and the reading."""
    members: dict[str, list[int]] = collections.defaultdict(list)
    for element, hit in in_set.items():
        members[chrom_of[element]].append(int(hit))
    keys = sorted(members)
    n = sum(len(v) for v in members.values())
    k = sum(sum(v) for v in members.values())
    rate = k / n if n else 0.0
    out: dict[str, Any] = {
        "elements": n,
        "in_the_set": k,
        "observed_rate": rate,
        "control": control,
        "ratio": rate / control if control else None,
        "clusters": len(keys),
        "clusters_are": "whole chromosomes",
    }
    if not n:
        out.update(reading="inconclusive", words=bands["readings"]["inconclusive"], why="no element")
        return out
    floor = int(bands["minimum_clusters_for_an_interval"])
    if len(keys) < floor:
        out["why_no_interval"] = (
            f"{len(keys)} chromosome clusters, below the registered {floor}. The registration names 1 Mb "
            "loci as the fallback cluster but fixes no rule for grouping them, so this runner does not "
            "choose one, and by the registered rule fewer than 10 clusters reads inconclusive"
        )
        out.update(reading="inconclusive", words=bands["readings"]["inconclusive"])
        return out
    draws, seed = int(bands["bootstrap_draws"]), int(bands["bootstrap_seed"])
    rng = random.Random(seed)
    stats = []
    for _ in range(draws):
        chosen = [rng.choice(keys) for _ in keys]
        succ = sum(sum(members[c]) for c in chosen)
        tot = sum(len(members[c]) for c in chosen)
        stats.append((succ / tot) / control)
    stats.sort()
    lo, hi = stats[int(0.025 * (len(stats) - 1))], stats[int(0.975 * (len(stats) - 1))]
    most = collections.Counter(stats).most_common(1)[0][1]
    distinct = len(set(stats))
    degenerate = k in (0, n) or distinct < int(bands["minimum_distinct_resamples"])
    out["bootstrap"] = {
        "ci95_ratio": [lo, hi],
        "draws": draws,
        "seed": seed,
        "identical_resample_share": most / len(stats),
        "distinct_resamples": distinct,
        "degenerate": degenerate,
    }
    if not degenerate:
        out.update(reading_h(lo, hi, rate, bands))
        out["decided_by"] = "the cluster bootstrap"
        return out
    sizes = {c: len(v) for c, v in members.items()}
    design = astroargmax.design(sizes)
    n_eff = n / design["design_effect_from_sizes"] if design["design_effect_from_sizes"] else n
    r_lo, r_hi = wilson(rate, n_eff)
    exact = reading_h(r_lo / control, r_hi / control, rate, bands)
    out["exact_bound"] = {
        "method": "Wilson score, 95%, on n_eff = n / design_effect (the registration permits "
        "Clopper-Pearson or Wilson, named which: this is Wilson)",
        "n_eff": n_eff,
        "design_effect_from_sizes": design["design_effect_from_sizes"],
        "rate_interval": [r_lo, r_hi],
        "ratio_interval": [r_lo / control, r_hi / control],
        "not_the_cluster_bootstrap": "this bound is NOT the cluster bootstrap's, which was degenerate",
    }
    # "the reading is INCONCLUSIVE unless the band is cleared by the exact bound alone"
    out.update(exact)
    out["decided_by"] = "the exact bound alone, the bootstrap being degenerate"
    return out


def item_h(
    answers: dict[str, dict[str, Any]], chrom_of: dict[str, str], root: Path | None = None
) -> dict[str, Any]:
    """Item (h) over a complete run's answers: one label per element, four sets, CNS_PRIMARY decides."""
    terms = item_h_registration(root)
    population = astroargmax.label_population({e: a.get("genes") or [] for e, a in answers.items()})
    per_cell = astroargmax.census()["per_cell"]
    labels = population["labels_by_element"]
    unjoinable = sorted(
        {lab for lab in labels.values() if not (per_cell.get(lab) or {}).get("ontology_term")}
    )
    sets = {}
    for name, curies in terms["sets"].items():
        members = set(curies)
        in_set = {e: (per_cell.get(lab) or {}).get("ontology_term") in members for e, lab in labels.items()}
        sets[name] = rate_against_control(in_set, chrom_of, terms["controls"][name], terms["bands"])
    widths = sorted(
        {
            int(axis.get("tracks_total") or 0)
            for a in answers.values()
            for axis in ((a.get("model") or {}).get("gene_axis_outputs") or [])
        }
    )
    return {
        "registration_sha256": terms["sha256"],
        "population": {k: v for k, v in population.items() if k != "labels_by_element"},
        "labels_the_census_cannot_join": unjoinable,
        "unjoinable_rule": (
            "a label the census carries no CURIE for cannot be placed in any set by the registered join, "
            "so it counts in the denominator and not the numerator, as a non-set label does in the "
            "control. Named here so the count is on the record"
        ),
        "by_set": sets,
        "deciding_set": DECIDING_SET,
        "the_reading": {"reading": sets[DECIDING_SET]["reading"], "words": sets[DECIDING_SET]["words"]},
        "alternatives_may_not_be_promoted": True,
        "axis_columns_observed": widths,
        "not_computed_here": ITEM_H_NOT_COMPUTED,
    }


#: What item (h)'s registration asks for that this runner does NOT compute, each with its reason.
ITEM_H_NOT_COMPUTED = {
    "comparison_b_the_K562_arm": (
        "registered to be computed at reading time from the deletion sweep's per-chromosome element "
        "tables restricted to the ENCODE K562 CRISPRi valid pairs, on the name-collapsed axis. Not built "
        "here: it needs the sweep's own collapse, which no committed function exposes"
    ),
    "the_667_column_and_name_collapsed_split": (
        "the registration expects each paid answer to carry 667 columns against the sweep's 371 names, "
        "while adapter v2's own note records 371 columns (astrorun.ADAPTER_NOTE_TRACKS). "
        "`axis_columns_observed` reports what the answers carry; where it is not 667 the registration's "
        "premise for the split does not hold, and settling that is an amendment's job, not a runner's"
    ),
    "secondary_all_gene_rows": (
        "amendment 1 part 3's labelled secondary decides nothing and fixes no per-row rule for a row "
        "with no movement; not computed rather than computed under a rule chosen here"
    ),
}


# ------------------------------------------------------------------------------ the primary claim

#: The terms the primary analysis needs and no committed registration or function fixes. Each one is a
#: choice that changes the number, which is why it is a REFUSAL here and not a default.
PRIMARY_GAPS = {
    "astrocyte_deletion_value": (
        "the registered path is crispri.annotate with `cells` extended to astrocyte, which reads each "
        "gene's `by_cell[cell]` through crispri.ElementCache. The recording path writes `by_cell` (and "
        "`cell_values`) for enhancer_target.CELLS / alphagenome_adapter.CELL_TRACKS only -- K562, HepG2, "
        "GM12878, IMR-90 -- so on the paid answers `by_cell['astrocyte']` is ABSENT and the registered "
        "path returns no astrocyte value for any pair. The astrocyte values exist only in the full track "
        "vectors kept by item (g). The registration's track_roster_rule says the frozen legacy per-name "
        "collapse applies, but no committed code applies it to a vector, so the rule that turns a vector "
        "into the deletion value must be registered before any is read"
    ),
    "distance": (
        "the frozen log_distance weight was fitted on the benchmark's own distances; the screen's pairs "
        "(astroreg_register.screen_pairs) carry coordinates and a gene symbol and no distance, and no "
        "committed function fixes the TSS source or the element point the distance is measured from"
    ),
    "coverage_controls": (
        "the registration requires the coverage-indicator control and the coverage-matched arm of the "
        "K562 result 'of any run of this design'; both are computed from the assembled pairs and wait on "
        "the two terms above"
    ),
}


def primary_tail(
    with_deletion: list[float],
    without_deletion: list[float],
    labels: list[bool],
    pairs: list[dict[str, Any]],
    root: Path | None = None,
) -> dict[str, Any]:
    """The registered tail over SCORED pairs: gain, both bootstraps, the wider quoted, the reading.

    `pairs` carry chrom, start, end and gene. The locus clusters are built exactly as the registered
    power figure built them (scripts/astroreg_register.py power_figure): fresh.locus_keys then
    cell2.group over every scored pair.
    """
    from genomeos.attribution import cell2, fresh

    base = Path(root) if root is not None else REPO
    reg2 = json.loads((base / ASTROREG2_REGISTRATION).read_text())
    keys = fresh.locus_keys(
        "astrocyte", [{"c": f"{p['chrom']}:{p['start']}-{p['end']}", "g": p["gene"]} for p in pairs], "c", "g"
    )
    loci = cell2.group(keys)
    locus = astroreg.cluster_bootstrap(with_deletion, without_deletion, labels, loci)
    chrom = astroreg.cluster_bootstrap(with_deletion, without_deletion, labels, [p["chrom"] for p in pairs])
    quoted = astroreg.wider(locus, chrom)
    return {
        "gain": locus["point"],
        "pairs": len(labels),
        "positives": sum(labels),
        "locus_interval": locus,
        "chromosome_interval": chrom,
        "quoted_interval": quoted,
        "reading": astroreg.reading(quoted),
        "interval_rule": astroreg.INTERVAL_RULE,
        "required_quote": reg2["claim"]["required_quote"],
        "travels_beside_it": reg2["claim"]["travels_beside_it"],
    }


def primary_assembly() -> dict[str, Any]:
    """The primary claim's assembly, refused with each unregistered term named."""
    return {
        "refused": True,
        "why": "the registrations do not fix every term the scored pairs need; each gap is named, and "
        "a blind amendment registered before any answer is read is what closes it",
        "gaps": dict(PRIMARY_GAPS),
        "tail_available": "primary_tail, tested on synthetic scores, runs once the assembly exists",
    }
