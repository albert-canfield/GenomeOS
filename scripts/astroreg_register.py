# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register the astrocyte CRISPRi test, and compute everything that can be known at 0 requests.

Writes `data/results/astroreg_registration.json`: the registered terms
(`genomeos/attribution/astroreg.py`), plus the four things the supervisor asked be computed before
Albert rules on the design --

  (a) the identity check: 20 K562 pairs recomputed through the extended cell path, shown unchanged;
  (b) the genome-wide coverage table, per chromosome, of the screen's elements against the swept
      element registry, with whether a cached astrocyte value exists for the measured gene;
  (c) the coverage rates of the screen's POSITIVES and of its NEGATIVES, separately, because a
      difference between them is a bias and not a detail;
  and the request figures, quoted against the HCT116 precedent.

No AlphaGenome request is made. No astrocyte pair is scored. The power figure is computed from the
K562 held-out pairs, whose deletion values are already in the local cache, so it too costs nothing.
"""

from __future__ import annotations

import ast
import gzip
import json
import random
import sys
from collections import defaultdict
from pathlib import Path

import openpyxl

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import astroreg, cell2, crispri, fresh  # noqa: E402
from genomeos.attribution.measured import CRISPRI_SPLIT_OF  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "astroreg_registration"
ROOT = Path(__file__).resolve().parents[1]
PACKAGE = "genomeos"
CACHE = Path("data/cache/fresh")
TABLE3 = CACHE / "astro_MOESM5.xlsx"
PAIR_SHEET = "3A_EnhancerGenePairDE"
POWER_SHEET = "3F_PowerCalculations"
IDENTITY_PAIRS = 20

OWN_CODE = (
    "genomeos/attribution/astroreg.py",
    "genomeos/attribution/fresh.py",
    "scripts/astroreg_register.py",
    "scripts/fresh_crispri_eligibility.py",
    "tests/test_astroreg.py",
    "tests/test_fresh_crispri.py",
)


def _is_file_exactly(root: Path, parts: list[str]) -> bool:
    """A file at `parts` below `root`, spelled as the directories spell it (carried from cd263bc)."""
    node = root
    for part in parts:
        try:
            if part not in {p.name for p in node.iterdir()}:
                return False
        except OSError:
            return False
        node = node / part
    return node.is_file()


def counting_path(entry: Path | None = None) -> list[str]:
    """The transitive import closure of this script over the repository's own package, from the
    files' own import statements. A hand-written list cannot be checked and goes stale."""
    entry = entry or Path(__file__).resolve()
    seen: dict[str, Path] = {}
    queue = [entry]
    while queue:
        path = queue.pop()
        rel = str(path.relative_to(ROOT))
        if rel in seen:
            continue
        seen[rel] = path
        try:
            tree = ast.parse(path.read_text())
        except (OSError, SyntaxError):
            continue
        names: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names.update(a.name for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module and not node.level:
                names.add(node.module)
                names.update(f"{node.module}.{a.name}" for a in node.names)
        for name in names:
            if not (name == PACKAGE or name.startswith(PACKAGE + ".")):
                continue
            parts = name.split(".")
            for candidate in ([*parts[:-1], f"{parts[-1]}.py"], [*parts, "__init__.py"]):
                if _is_file_exactly(ROOT, candidate):
                    queue.append(ROOT.joinpath(*candidate))
    return sorted(seen)


def code_cleanliness() -> dict:
    rev = mf.code_revision()
    path = counting_path()
    dirty = list(rev["dirty_code_paths"])
    own = [p for p in dirty if p in OWN_CODE]
    foreign = [p for p in dirty if p not in OWN_CODE]
    return {
        "git_sha": rev["git_sha"],
        "dirty": rev["dirty"],
        "own_uncommitted_code": own,
        "own_code_is_committed": not own,
        "foreign_uncommitted_code": foreign,
        "foreign_uncommitted_code_on_the_counting_path": [p for p in foreign if p in path],
        "counting_path": path,
        "counting_path_is_computed": (
            "the transitive import closure of this script over the repository's own package, "
            "computed from the files' import statements at write time; every path component is "
            "matched against what its directory actually lists, so the closure is the same on a "
            "case-sensitive filesystem (carried from cd263bc)"
        ),
        "note": (
            "several lanes work in this one checkout. A file under foreign_uncommitted_code belongs "
            "to another lane; this lane did not write it and did not commit it. The counting path is "
            "the computed closure above, so a foreign file outside it cannot have entered any number "
            "here, and foreign_uncommitted_code_on_the_counting_path names any that could"
        ),
    }


# ---------------------------------------------------------------- the screen's pairs and labels


def read_sheet(path: Path, sheet: str, columns: list[str]) -> list[dict]:
    """Rows of `sheet`, restricted to `columns`, each of which must pass `fresh.permit` first."""
    fresh.permit(columns)
    workbook = openpyxl.load_workbook(path, read_only=True)
    try:
        worksheet = workbook[sheet]
        rows = worksheet.iter_rows(values_only=True)
        header = list(next(rows))
        missing = [c for c in columns if c not in header]
        if missing:
            raise KeyError(f"{path.name}:{sheet} has no column {missing}")
        index = {c: header.index(c) for c in columns}
        out, blank = [], 0
        for row in rows:
            if not any(cell is not None for cell in row):
                blank += 1
                continue
            if len(row) <= max(index.values()):
                raise ValueError(f"{path.name}:{sheet} has a short non-blank row: {row!r}")
            out.append({c: row[index[c]] for c in columns})
    finally:
        workbook.close()
    print(f"  {path.name}:{sheet}: {len(out)} rows, {blank} wholly blank skipped")
    return out


def _true(v) -> bool:
    return str(v).strip().lower() == "true"


def screen_pairs() -> list[dict]:
    """Every screen pair with its registered label, coordinates and measured gene.

    The two sheets are joined on `Pair`, and the join must be total in both directions: a pair in one
    sheet and not the other would mean one of them is being misread.
    """
    a = read_sheet(TABLE3, PAIR_SHEET, ["Pair", "Enhancer", "GeneSymbol", "Hit", "EnhancerCoord"])
    f = read_sheet(
        TABLE3,
        POWER_SHEET,
        ["Pair", "Hit", "Hit&Downregulated", "WellPowered_at_FC_0.15", "WellPowered_at_FC_0.25"],
    )
    power = {r["Pair"]: r for r in f}
    if len(power) != len(f):
        raise ValueError("the power sheet repeats a pair id")
    missing = [r["Pair"] for r in a if r["Pair"] not in power]
    if missing:
        raise ValueError(f"{len(missing)} pairs of the pair sheet are absent from the power sheet")
    if len(a) != len(f):
        raise ValueError(f"the two sheets disagree on pair count: {len(a)} vs {len(f)}")

    out = []
    for r in a:
        p = power[r["Pair"]]
        if _true(r["Hit"]) != _true(p["Hit"]):
            raise ValueError(f"the two sheets disagree on the hit call for {r['Pair']}")
        chrom, start, end = fresh.parse_coord(r["EnhancerCoord"])
        out.append(
            {
                "pair": r["Pair"],
                "element": r["Enhancer"],
                "gene": str(r["GeneSymbol"]),
                "chrom": chrom,
                "start": start,
                "end": end,
                "label": astroreg.label_of(
                    _true(r["Hit"]),
                    _true(p["Hit&Downregulated"]),
                    _true(p["WellPowered_at_FC_0.25"]),
                ),
                "well_powered_015": _true(p["WellPowered_at_FC_0.15"]),
            }
        )
    return out


# ---------------------------------------------------------------- (a) the identity check


def identity_check() -> dict:
    """20 K562 pairs recomputed with astrocyte added to `cells`, and shown unchanged.

    `crispri._annotate_one` reaches the deletion features through `if p.cell in cells`, so extending
    the tuple is a pure widening of a membership gate. This runs it both ways and compares the
    features leaf by leaf, rather than asserting it from the source.
    """
    with gzip.open(crispri.KNOWLEDGE / crispri.HELDOUT, "rt") as fh:
        lines = fh.read().splitlines()
    pairs = [p for p in crispri.parse(lines) if p.cell == "K562"]
    if not pairs:
        # the held-out file carries no K562; the training file is the K562 set
        with gzip.open(crispri.KNOWLEDGE / crispri.TRAINING, "rt") as fh:
            lines = fh.read().splitlines()
        pairs = [p for p in crispri.parse(lines) if p.cell == "K562"]
    sample = sorted(pairs, key=lambda p: (p.chrom, p.start, p.gene))[:IDENTITY_PAIRS]
    table, cache = crispri.DeletionTable(), crispri.ElementCache()

    def features(cells):
        rows = [crispri.Pair(**{k: getattr(p, k) for k in _pair_fields(p)}) for p in sample]
        crispri.annotate(rows, table, cache, cells=cells)
        return [dict(r.features) for r in rows], [r.covered for r in rows]

    frozen_f, frozen_c = features(crispri.MODEL_CELLS)
    extended_f, extended_c = features(crispri.MODEL_CELLS + ("astrocyte",))
    differing = [
        {"pair": sample[i].gene + "|" + sample[i].chrom, "frozen": frozen_f[i], "extended": extended_f[i]}
        for i in range(len(sample))
        if frozen_f[i] != extended_f[i] or frozen_c[i] != extended_c[i]
    ]
    return {
        "pairs_compared": len(sample),
        "cells_frozen": list(crispri.MODEL_CELLS),
        "cells_extended": list(crispri.MODEL_CELLS + ("astrocyte",)),
        "feature_leaves_compared": sum(len(f) for f in frozen_f) + len(frozen_c),
        "pairs_that_differ": len(differing),
        "differences": differing,
        "unchanged": not differing,
        "why": astroreg.SAME_CODE_PATH,
    }


def _pair_fields(p) -> list[str]:
    """The constructor fields of a `crispri.Pair`, so a fresh copy can be made without sharing
    mutable feature dictionaries between the two runs."""
    import dataclasses

    return [f.name for f in dataclasses.fields(p) if f.init]


# ------------------------------------------------- (b) and (c) genome-wide coverage and bias


def coverage(pairs: list[dict]) -> dict:
    """Per chromosome: overlap with the swept registry, cached astrocyte values, requests needed.

    A pair is `covered` when a swept registry element overlaps its element, which is the same
    condition `crispri.annotate` uses for `covered`. An element with no swept counterpart cannot be
    requested at all, so the overlap half of this table decides what is even purchasable.

    `astrocyte_value_cached` asks only whether a value is present, never what it is.
    """
    table, cache = crispri.DeletionTable(), crispri.ElementCache()
    per_chrom: dict[str, dict] = {}
    by_label_cov: dict[str, dict[str, int]] = defaultdict(lambda: {"covered": 0, "total": 0})
    needed_elements: set[tuple[str, str]] = set()
    cached_astro = 0
    covered_positive_loci: list = []

    for chrom in sorted({p["chrom"] for p in pairs}, key=lambda c: (len(c), c)):
        rows = [p for p in pairs if p["chrom"] == chrom]
        stat = {
            "pairs": len(rows),
            "covered_pairs": 0,
            "elements": len({r["element"] for r in rows}),
            "covered_elements": 0,
            "registry_elements_needed": 0,
            "astrocyte_value_cached": 0,
        }
        covered_els: set[str] = set()
        for r in rows:
            els = table.overlapping(chrom, r["start"], r["end"])
            cov = bool(els)
            r["covered"] = cov
            by_label_cov[r["label"]]["total"] += 1
            if cov:
                by_label_cov[r["label"]]["covered"] += 1
                stat["covered_pairs"] += 1
                covered_els.add(r["element"])
                for e in els:
                    needed_elements.add((chrom, e["id"]))
                    # presence only: the value itself is never read into a number here
                    if cache.value(chrom, e["id"], r["gene"], "astrocyte") is not None:
                        cached_astro += 1
                        stat["astrocyte_value_cached"] += 1
                if r["label"] == "positive":
                    covered_positive_loci.append(r)
        stat["covered_elements"] = len(covered_els)
        stat["registry_elements_needed"] = len({e for c, e in needed_elements if c == chrom})
        per_chrom[chrom] = stat
        print(
            f"  {chrom}: {stat['covered_pairs']}/{stat['pairs']} pairs covered, "
            f"{stat['registry_elements_needed']} registry elements"
        )

    rates = {}
    for label, d in by_label_cov.items():
        rates[label] = {
            "total": d["total"],
            "covered": d["covered"],
            "covered_share": (d["covered"] / d["total"]) if d["total"] else None,
        }
    pos, neg = rates.get("positive", {}), rates.get("negative", {})
    ps, ns = pos.get("covered_share"), neg.get("covered_share")
    bias = {
        "positive_covered_share": ps,
        "negative_covered_share": ns,
        "difference_positive_minus_negative": (ps - ns) if (ps is not None and ns is not None) else None,
        "reading": (
            "the swept elements cover the screen's positives and its negatives at different rates, "
            "which is a bias in what the test can see and is reported as one"
            if (ps is not None and ns is not None and abs(ps - ns) > 0.01)
            else "the covered shares of positives and negatives are within one percentage point"
        ),
        "controls_that_apply": (
            "the coverage-indicator control and the coverage-matched arm already in use for the K562 "
            "result apply here unchanged, and both are required of any run of this design"
        ),
    }

    keys = cell2.group(
        fresh.locus_keys(
            "astrocyte",
            [{"c": f"{r['chrom']}:{r['start']}-{r['end']}", "g": r["gene"]} for r in covered_positive_loci],
            "c",
            "g",
        )
    )
    return {
        "per_chromosome": per_chrom,
        "by_label": rates,
        "coverage_bias": bias,
        "astrocyte_values_already_cached": cached_astro,
        "astrocyte_values_already_cached_note": (
            "counted as presence only; no cached value was read as a number. Expected to be nil, "
            "because the sweep retained four cell lines' values per gene and astrocyte is not one"
        ),
        "registry_elements_needed_total": len(needed_elements),
        "covered_positive_loci": len(set(keys)),
    }


def requests(cov: dict) -> dict:
    """What the information costs, in the terms the decision is made in."""
    total = cov["registry_elements_needed_total"] - cov["astrocyte_values_already_cached"]
    pos = cov["by_label"].get("positive", {})
    neg = cov["by_label"].get("negative", {})
    covered_pos = pos.get("covered") or 0
    per_positive = (total / covered_pos) if covered_pos else None
    hct = astroreg.HCT116_PRECEDENT
    return {
        "total_requests_needed": total,
        "one_request_per": "registry element overlapping a covered pair, as the HCT116 arm counted it",
        "covered_share_of_positives": pos.get("covered_share"),
        "covered_positives": covered_pos,
        "total_positives": pos.get("total"),
        "covered_share_of_negatives": neg.get("covered_share"),
        "covered_negatives": neg.get("covered"),
        "total_negatives": neg.get("total"),
        "expected_requests_per_positive": per_positive,
        "covered_independent_loci_among_positives": cov["covered_positive_loci"],
        "against_hct116_precedent": {
            "hct116_requests": hct["requests"],
            "hct116_covered_pairs": hct["covered_pairs"],
            "hct116_requests_per_covered_pair": hct["requests"] / hct["covered_pairs"],
            "ours_requests_per_covered_positive": per_positive,
            "ratio_of_totals": (total / hct["requests"]) if hct["requests"] else None,
        },
        "framing": astroreg.COST_FRAMING,
        "requests_sent": 0,
    }


# ---------------------------------------------------------------- the power figure


def power_figure() -> dict:
    """K562's own gain, and the share of resamples to 86 loci whose interval excludes zero.

    Costs nothing: K562's deletion values are already in the local cache. The frozen weights are
    refitted on the K562 TRAINING pairs exactly as the committed result does, and no astrocyte pair
    is touched.
    """
    table, cache = crispri.DeletionTable(), crispri.ElementCache()
    with gzip.open(crispri.KNOWLEDGE / crispri.TRAINING, "rt") as fh:
        train = crispri.parse(fh.read().splitlines())
    with gzip.open(crispri.KNOWLEDGE / crispri.HELDOUT, "rt") as fh:
        held = crispri.parse(fh.read().splitlines())
    print(f"  annotating {len(train)} training and {len(held)} held-out pairs from the local cache")
    crispri.annotate(train + held, table, cache)
    tr = [p for p in train if p.covered]
    labels = [p.regulated for p in tr]
    weights = {n: crispri.logistic_fit(crispri.matrix(tr, c), labels) for n, c in crispri.FEATURES.items()}

    k562 = [p for p in held if p.cell == "K562" and p.covered] or tr
    scope = "held-out K562" if any(p.cell == "K562" for p in held) else "K562 training (no held-out K562)"
    lab = [p.regulated for p in k562]
    with_d = crispri.logistic_score(
        weights["activity + distance + deletion"],
        crispri.matrix(k562, crispri.FEATURES["activity + distance + deletion"]),
    )
    without = crispri.logistic_score(
        weights["activity + distance"], crispri.matrix(k562, crispri.FEATURES["activity + distance"])
    )

    keys = fresh.locus_keys(
        "K562",
        [{"c": f"{p.chrom}:{p.start}-{p.end}", "g": p.gene} for p in k562],
        "c",
        "g",
    )
    loci = cell2.group(keys)
    observed = astroreg.cluster_bootstrap(with_d, without, lab, loci)

    members: dict[int, list[int]] = defaultdict(list)
    for i, g in enumerate(loci):
        members[g].append(i)
    pool = sorted(members)
    k562_prevalence = sum(lab) / len(lab)
    target = astroreg.ASTROCYTE_PREVALENCE
    draws = 200  # outer resamples; each carries its own inner interval

    arms: dict[str, dict] = {}
    for k in astroreg.ATTENUATIONS:
        attenuated = astroreg.attenuate(with_d, without, k)
        rng = random.Random(0)  # the same draws at every k, so only the attenuation differs
        excludes, usable, kept_pos, kept_n, gains = 0, 0, [], [], []
        for _ in range(draws):
            chosen = [rng.choice(pool) for _ in range(astroreg.REGISTERED_LOCI)]
            idx, cl = [], []
            for c in chosen:
                for i in members[c]:
                    idx.append(i)
                    cl.append(c)
            of_cluster = dict(zip(idx, cl, strict=True))
            thinned = astroreg.match_prevalence(idx, lab, target, rng)
            sub_lab = [lab[i] for i in thinned]
            if not any(sub_lab):
                continue
            usable += 1
            kept_pos.append(sum(sub_lab))
            kept_n.append(len(thinned))
            sub = astroreg.cluster_bootstrap(
                [attenuated[i] for i in thinned],
                [without[i] for i in thinned],
                sub_lab,
                [of_cluster[i] for i in thinned],
                seed=usable,
                n=200,
            )
            gains.append(sub["point"])
            if "ci95" in sub and (sub["ci95"][0] > 0 or sub["ci95"][1] < 0):
                excludes += 1
        share = (excludes / usable) if usable else 0.0
        gains.sort()
        arms[f"{k:g}x"] = {
            "attenuation": k,
            "share_of_resamples_excluding_zero": share,
            "usable_outer_draws": usable,
            "intervals_excluding_zero": excludes,
            "median_achieved_gain": (gains[len(gains) // 2] if gains else None),
            "mean_positives_kept": (sum(kept_pos) / len(kept_pos)) if kept_pos else None,
            "mean_pairs_kept": (sum(kept_n) / len(kept_n)) if kept_n else None,
            "achieved_prevalence": ((sum(kept_pos) / sum(kept_n)) if kept_n and sum(kept_n) else None),
        }

    full = arms[f"{1.0:g}x"]
    return {
        "method": astroreg.POWER_RULE,
        "amendment_1": astroreg.POWER_AMENDMENT_1,
        "k562_scope": scope,
        "k562_pairs": len(k562),
        "k562_positives": sum(lab),
        "k562_prevalence": k562_prevalence,
        "astrocyte_prevalence_target": target,
        "prevalence_mismatch_ratio": (k562_prevalence / target) if target else None,
        "k562_independent_loci": len(set(loci)),
        "k562_observed_gain_at_k562_prevalence": observed,
        "resampled_to_loci": astroreg.REGISTERED_LOCI,
        "outer_draws": draws,
        "by_attenuation": arms,
        "class_set_by": "the prevalence-matched figure at the full effect (1x)",
        "verdict": astroreg.power_class(full["share_of_resamples_excluding_zero"]),
        "travels_beside_the_class": {
            "0.5x_share": arms[f"{0.5:g}x"]["share_of_resamples_excluding_zero"],
            "rule": "the 0.5x figure travels beside the class in every quote of it",
        },
        "superseded": {
            "earlier_share": 0.96,
            "why": "computed at K562's own prevalence of 6.5% and at the full K562 gain; it does "
            "not describe a test run at 2.9% prevalence",
        },
        "alphagenome_requests": 0,
    }


# ---------------------------------------------------------------- manifest and main


def manifest(pairs: list[dict]) -> dict:
    url = (
        "https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-025-02154-3"
        "/MediaObjects/41593_2025_2154_MOESM5_ESM.xlsx"
    )
    return {
        "sources": [
            {
                "accession": "Green NFO, et al. CRISPRi screening in cultured human astrocytes. "
                "Nature Neuroscience 2025;29(3):703-716 -- Supplementary Table 3",
                "version": "publisher open-access supplementary store; pinned here by sha256",
                "url": url,
            },
            {
                "accession": "EngreitzLab/CRISPR_comparison, EPCrisprBenchmark training_K562 and "
                "heldout_5_cell_types (Gschwind et al.)",
                "version": "fetched 2026-09-16; pinned here by sha256",
                "url": crispri.BASE_URL,
            },
            {
                "accession": "data/results/fresh_crispri_eligibility.json",
                "version": "this lane's committed eligibility result; its counts are carried, not re-derived",
            },
        ],
        "inputs": [
            mf.input_entry(
                TABLE3,
                partition=None,
                role="the screen's hit calls, direction labels, power flags, element coordinates and "
                "measured gene symbols; no effect size, p-value, FDR or score was read",
            ),
            mf.input_entry(
                crispri.KNOWLEDGE / crispri.TRAINING,
                partition=CRISPRI_SPLIT_OF[crispri.TRAINING],
                role="the K562 pairs the frozen weights are fitted on, for the identity check and "
                "the power figure; already-read development evidence",
            ),
            mf.input_entry(
                crispri.KNOWLEDGE / crispri.HELDOUT,
                partition=CRISPRI_SPLIT_OF[crispri.HELDOUT],
                role="the held-out pairs, for the power figure's K562 scope; already-read",
            ),
        ],
        "assembly": "GRCh38",
        "coordinates": (
            "n/a: the screen states no convention for its chrom:start-end strings. The uses made of "
            "them here are the 1 Mb locus grouping and an overlap test against the registry, "
            "neither of whose outcome a one-base shift can change; the build itself is confirmed "
            "hg38 (see terms.build)"
        ),
        "parameters": {
            "registered_loci": astroreg.REGISTERED_LOCI,
            "registered_positives": astroreg.REGISTERED_POSITIVES,
            "independent_locus_span_bp": cell2.INDEPENDENT_LOCUS_SPAN,
            "bootstraps": crispri.BOOTSTRAPS,
            "min_clusters_for_an_interval": crispri.MIN_CLUSTERS_FOR_AN_INTERVAL,
            "identity_check_pairs": IDENTITY_PAIRS,
            "model_cells_frozen": list(crispri.MODEL_CELLS),
            "permitted_columns": sorted(fresh.PERMITTED_COLUMNS),
            "alphagenome_requests": 0,
            "model_requests": 0,
            "astrocyte_pairs_scored": 0,
            "money_spent": 0,
        },
        "exclusions": [
            "no effect size, fold change, p-value, FDR, Z statistic or expression level was read "
            "from the screen: fresh.permit refuses the column before the workbook is opened",
            "the 224 cached rows whose winning track is astrocyte are NOT reused as values; "
            "astroreg.CACHE_REUSE_REFUSED gives the reason",
            "no astrocyte pair is scored and no model request is made: this is a registration",
            "a hit whose change is an increase is neither a positive nor a negative and is never "
            "folded into the negatives",
            "a non-hit that is not WellPowered at fc 0.25 is excluded rather than counted as a negative",
        ],
        "partitions": {
            "screen_pairs": f"all {len(pairs)} element-gene pairs of the published screen",
            "k562": "already-read development evidence, used only for the identity check and the "
            "power figure",
        },
        "code_cleanliness": code_cleanliness(),
    }


def main() -> int:
    print("Reading the screen's labels and coordinates (no effect size is read)")
    pairs = screen_pairs()
    counts = defaultdict(int)
    for p in pairs:
        counts[p["label"]] += 1
    print(f"  labels: {dict(counts)}")

    print("(a) identity check: 20 K562 pairs through the extended cell path")
    identity = identity_check()
    print(f"  pairs that differ: {identity['pairs_that_differ']} of {identity['pairs_compared']}")

    print("(b,c) genome-wide coverage against the swept registry")
    cov = coverage(pairs)
    req = requests(cov)

    print("power figure from the K562 pairs already in the cache")
    power = power_figure()

    payload = {
        "terms": astroreg.terms(),
        "label_counts": dict(counts),
        "label_counts_note": (
            "positives are the authors' significant decreases; increases are held apart; a non-hit "
            "that is not WellPowered at fc 0.25 is excluded rather than counted as a negative"
        ),
        "identity_check": identity,
        "coverage": cov,
        "requests": req,
        "power": power,
        "left_undone": [
            "no astrocyte pair is scored and no request is sent: the design is registered, not run",
            "the ENCODE accessions of the astrocyte DNase and H3K27ac files already on disk are not "
            "recorded anywhere in the repository; they are present but unprovenanced",
            "whether AlphaGenome's astrocyte output is one track or several is not established; the "
            "cache records a winning name, not the roster",
            "the coverage table's requests assume one request per registry element overlapping a "
            "covered pair, the HCT116 counting rule; a cheaper rule may exist and is not explored here",
        ],
    }
    path = save_result(RESULT, payload, manifest=manifest(pairs))
    print(
        json.dumps(
            {
                "identity_unchanged": identity["unchanged"],
                "total_requests": req["total_requests_needed"],
                "covered_share_positives": req["covered_share_of_positives"],
                "covered_share_negatives": req["covered_share_of_negatives"],
                "requests_per_positive": req["expected_requests_per_positive"],
                "power": power["verdict"],
            },
            indent=1,
        )
    )
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
