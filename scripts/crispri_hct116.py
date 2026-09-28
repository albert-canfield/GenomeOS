# SPDX-License-Identifier: AGPL-3.0-or-later
"""The HCT116 arm of the CRISPRi result: the frozen model in a second cell type, as registered.

    uv run python scripts/crispri_hct116.py --dry-run            # the request count, no request
    uv run python scripts/crispri_hct116.py --fetch --cap 760    # resumable; stops at the cap
    uv run python scripts/crispri_hct116.py --score              # no request; writes the result

The registration is crispri.PREREGISTERED_PUBLISHED['second_cell_type'] (fe61c36, 2026-09-27):
"the same code path with HCT116 added to the cells kept per gene, weights unchanged; passes if the
gain on the 363 covered HCT116 pairs is above zero; stated as replicated only if its
chromosome-bootstrap interval excludes zero".

Fetch. One deletion request per registry element overlapping a held-out HCT116 pair, made exactly
as the sweep made it (enhancer_target.score_element's deletion: anchor base plus element, the
RNA_SEQ gene scorer on a 1 Mb window, threshold 0.0), through alphagenome_adapter.create_client, so
every request asks for ALL_FOLDS. The answer is written in the sweep's cache format (one JSON per
element, `genes` rows with `by_cell`) with HCT116 added to the cells kept per gene, under its own
root (ELEMENT_CACHE_HCT116), never over the sweep's cache: the sweep's answers are the ones the
frozen weights were fitted on and stay as they are. Each answer carries its run record (`model`)
and `model_version`. Every attempt is appended to a ledger before it is sent, so a restart never
asks an answered element again and the count of requests spent is read from disk.

Score. The frozen weights are refitted exactly as score_published fits them (covered K562 training
pairs, the sweep's cache), then applied unchanged to the HCT116 pairs, whose deletion values come
from the HCT116 cache on HCT116's own track. The registration names the pairs (the 363 covered),
the sign rule and the chromosome bootstrap but not the AUPRC estimator; this script uses the
headline estimator of the K562 claim it replicates (average precision on the covered pairs,
crispri.gain_interval), and reports score_published's per-cell block (all 396 pairs, weighted,
the benchmark's estimator) beside it. That choice was written here before any HCT116 answer was
read.
"""

from __future__ import annotations

import argparse
import json
import sys
import threading
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from genomeos import manifest as mf
from genomeos.attribution import crispri
from genomeos.predict.enhancer_target import CELLS, aggregate, model_version_of, predict_target

CELL = "HCT116"
KEPT_CELLS = (*CELLS, CELL)  # the cells kept per gene: the sweep's four plus HCT116
ELEMENT_CACHE_HCT116 = Path("data/knowledge/alphagenome/elements_hct116")
LEDGER = ELEMENT_CACHE_HCT116 / "ledger.jsonl"
COSTED = 705
CAP = 760  # owner's approval, 2026-09-28: 705 costed plus a retry margin
WORKERS = 4
CALL_TIMEOUT = 300
RESULT = Path("data/results/crispri_published.json")
KEY = "second_cell_type_hct116"
INSERT_BEFORE = "post_hoc_positive_filter"  # a key inserted before an existing one changes no line


def to_score(heldout: list[crispri.Pair], table: crispri.DeletionTable) -> dict[str, dict[str, Any]]:
    """Every registry element overlapping a held-out HCT116 pair, by id (the frozen feature's inputs)."""
    out: dict[str, dict[str, Any]] = {}
    for p in heldout:
        if p.cell != CELL:
            continue
        for e in table.overlapping(p.chrom, p.start, p.end):
            out[e["id"]] = {"id": e["id"], "chrom": p.chrom, "start": e["start"], "end": e["end"]}
    return out


def answer_path(chrom: str, element_id: str) -> Path:
    return ELEMENT_CACHE_HCT116 / chrom / f"{element_id}.json"


def ledger() -> list[dict[str, Any]]:
    if not LEDGER.exists():
        return []
    return [json.loads(line) for line in LEDGER.read_text().splitlines() if line.strip()]


def spent(rows: list[dict[str, Any]]) -> dict[str, int]:
    """Requests sent, split by outcome; a quota refusal is sent but not answered."""
    c = Counter(r["event"] for r in rows)
    return {
        "sent": c["request"],
        "answered": c["answer"],
        "refused_quota": c["quota"],
        "failed_other": c["error"],
    }


class Fetcher:
    """One client (ALL_FOLDS), a few requests in flight, every attempt on the ledger before it is sent."""

    def __init__(self, cap: int) -> None:
        from genomeos.predict import AlphaGenomeAdapter
        from genomeos.predict.alphagenome_adapter import create_client

        self.cap = cap
        self.lock = threading.Lock()
        self.genomes: dict[str, Any] = {}
        a = AlphaGenomeAdapter()
        a._client = create_client(a.api_key, timeout=CALL_TIMEOUT)  # noqa: SLF001  (asks for ALL_FOLDS)
        self.scorer = a._live_scorer(threshold=0.0)  # noqa: SLF001  (the sweep's threshold)
        self.sent = spent(ledger())["sent"]
        self.quota_until = 0.0

    def log(self, **row: Any) -> None:
        row["t"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        LEDGER.parent.mkdir(parents=True, exist_ok=True)
        with LEDGER.open("a") as fh:
            fh.write(json.dumps(row) + "\n")

    def seq(self, chrom: str, start: int, end: int) -> str:
        from genomeos.genome import IndexedGenome, Locus, reference_fasta

        with self.lock:
            if chrom not in self.genomes:
                self.genomes[chrom] = IndexedGenome(str(reference_fasta(chrom)))
            return str(self.genomes[chrom].fetch(Locus(chrom, start - 1, end))).upper()

    def one(self, e: dict[str, Any]) -> str:
        """One element: answered already, or one request; returns what happened."""
        path = answer_path(e["chrom"], e["id"])
        while True:
            if path.exists():
                return "cached"
            wait = self.quota_until - time.time()
            if wait > 0:
                time.sleep(min(wait, 10))
                continue
            with self.lock:
                if self.sent >= self.cap:
                    return "cap"
                self.sent += 1
                self.log(event="request", id=e["id"], n=self.sent)
            seq = self.seq(e["chrom"], e["start"], e["end"])
            t0 = time.time()
            try:
                effects = self.scorer(e["chrom"], e["start"], seq, seq[0])  # as score_element deletes
            except Exception as ex:  # noqa: BLE001
                msg = str(ex)
                if "RESOURCE_EXHAUSTED" in msg or "Quota" in msg:
                    with self.lock:
                        self.log(event="quota", id=e["id"], message=msg[:200])
                        self.quota_until = time.time() + 120
                    continue
                with self.lock:
                    self.log(event="error", id=e["id"], message=f"{type(ex).__name__}: {msg[:200]}")
                return "error"
            rows = aggregate(effects, KEPT_CELLS)
            per_gene: dict[str, list[float]] = {}
            for g, t, v in effects:
                if t == CELL:
                    per_gene.setdefault(g, []).append(round(float(v), 4))
            hit = {
                "id": e["id"],
                "chrom": e["chrom"],
                "start": e["start"],
                "end": e["end"],
                "length": e["end"] - e["start"],
                "genes_in_window": len(rows),
                "tracks": rows[0]["n_tracks"] if rows else 0,
                "seconds": round(time.time() - t0, 1),
                "genes": rows,
                "model": dict(self.scorer.model),  # type: ignore[attr-defined]
                "cells_kept": list(KEPT_CELLS),
                # every HCT116 track's value per gene in track order; by_cell keeps the last one, as
                # the sweep's aggregate does for K562 (the registered code path)
                "hct116_tracks_by_gene": per_gene,
            }
            hit["model_version"] = model_version_of(hit)
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp = path.with_suffix(".part")
            tmp.write_text(json.dumps(hit))
            tmp.rename(path)
            with self.lock:
                self.log(
                    event="answer", id=e["id"], model_version=hit["model_version"], hct116_genes=len(per_gene)
                )
            return "answered"


def fetch(elements: dict[str, dict[str, Any]], cap: int, workers: int, limit: int | None) -> int:
    from genomeos.predict import status

    st = status()
    if not st["enabled"]:
        print(f"AlphaGenome is disabled: {st['reason']}")
        return 2
    todo = [e for e in elements.values() if not answer_path(e["chrom"], e["id"]).exists()]
    todo.sort(key=lambda e: (e["chrom"], e["start"]))
    if limit:
        todo = todo[:limit]
    f = Fetcher(cap)
    print(f"{len(todo)} to ask; {f.sent} sent before this run; cap {cap}", flush=True)
    outcome: Counter[str] = Counter()
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(f.one, e) for e in todo]
        for i, fut in enumerate(as_completed(futures), 1):
            outcome[fut.result()] += 1
            if i % 25 == 0 or i == len(futures):
                print(f"  {i}/{len(todo)} {dict(outcome)}; sent so far {f.sent}", flush=True)
    print(f"done: {dict(outcome)}; ledger {spent(ledger())}")
    return 0 if not outcome["cap"] and not outcome["error"] else 1


def compare_with_sweep(elements: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """The same element asked twice, once with no version requested (the sweep) and once for ALL_FOLDS.

    Reads the sweep's cache one chromosome at a time (one reader) and the HCT116 answers; no request.
    Agreement on the sweep's four cell lines says whether the two requests reached the same model.
    """
    old_cache = crispri.ElementCache()
    diffs: list[float] = []
    exact = n_elements = top_same = top_n = 0
    for e in sorted(elements.values(), key=lambda e: (e["chrom"], e["start"])):
        p = answer_path(e["chrom"], e["id"])
        if not p.exists():
            continue
        new = json.loads(p.read_text())
        old_cache._load(e["chrom"])  # noqa: SLF001
        old = old_cache._archive.get(e["id"])  # noqa: SLF001
        if old is None:
            q = old_cache.root / e["chrom"] / f"{e['id']}.json"
            old = json.loads(q.read_text()) if q.exists() else None
        if not old:
            continue
        n_elements += 1
        o = {g["gene"]: g for g in old.get("genes") or []}
        for g in new["genes"]:
            h = o.get(g["gene"])
            if not h:
                continue
            for cell in CELLS:
                a, b = (g.get("by_cell") or {}).get(cell), (h.get("by_cell") or {}).get(cell)
                if a is None or b is None:
                    continue
                diffs.append(abs(a - b))
                exact += a == b
        pn, po = predict_target(new["genes"]), predict_target(old.get("genes") or [])
        top_n += 1
        top_same += (pn or {}).get("gene") == (po or {}).get("gene")
    diffs.sort()
    return {
        "elements_in_both": n_elements,
        "gene_cell_values_compared": len(diffs),
        "identical_to_4_places": exact,
        "median_abs_difference": diffs[len(diffs) // 2] if diffs else None,
        "p99_abs_difference": diffs[int(0.99 * (len(diffs) - 1))] if diffs else None,
        "max_abs_difference": diffs[-1] if diffs else None,
        "same_top_target": top_same,
        "top_target_compared": top_n,
        "cells": list(CELLS),
    }


def score(
    training: list[crispri.Pair], heldout: list[crispri.Pair], elements: dict[str, Any]
) -> dict[str, Any]:
    table = crispri.DeletionTable()
    hct = [p for p in heldout if p.cell == CELL]
    others = [p for p in heldout if p.cell != CELL]
    # 1. the frozen weights, exactly as score_published fits them (the sweep's cache, one reader)
    crispri.annotate(training + others, table, crispri.ElementCache())
    for p in training + heldout:
        p.features["covered"] = float(p.covered)
    train = [p for p in training if p.covered]
    sweep_answers = len({e["id"] for p in train for e in table.overlapping(p.chrom, p.start, p.end)})
    weights = {
        name: crispri.logistic_fit(crispri.matrix(train, cols), [p.regulated for p in train])
        for name, cols in crispri.FEATURES.items()
    }
    with_d, without = "activity + distance + deletion", "activity + distance"
    k562 = [p for p in others if p.cell == "K562" and p.covered]
    check = crispri.gain_interval(
        crispri.logistic_score(weights[with_d], crispri.matrix(k562, crispri.FEATURES[with_d])),
        crispri.logistic_score(weights[without], crispri.matrix(k562, crispri.FEATURES[without])),
        k562,
    )
    # 2. HCT116 added to the cells kept per gene, its values from its own cache; weights unchanged
    crispri.annotate(
        hct, table, crispri.ElementCache(ELEMENT_CACHE_HCT116), cells=(*crispri.MODEL_CELLS, CELL)
    )
    for p in hct:
        p.features["covered"] = float(p.covered)
    covered = [p for p in hct if p.covered]
    lab = [p.regulated for p in covered]
    s = {
        n: crispri.logistic_score(weights[n], crispri.matrix(covered, c)) for n, c in crispri.FEATURES.items()
    }
    registered = crispri.gain_interval(s[with_d], s[without], covered)
    ci = registered["ci95"]
    passed = registered["gain"] > 0
    replicated = bool(passed and ci and ci[0] > 0)
    verdict = (
        "replicated: the gain is above zero and its chromosome-bootstrap interval excludes zero"
        if replicated
        else (
            "passes, not replicated: the gain is above zero and its interval includes zero"
            if passed
            else "fails: the gain is at or below zero, so the result is a K562 result and is described as one"
        )
    )
    # 3. beside it, score_published's per-cell block for HCT116 now that it has a deletion value
    s_all = {
        n: crispri.logistic_score(weights[n], crispri.matrix(hct, c)) for n, c in crispri.FEATURES.items()
    }
    per_cell = {
        "pairs": len(hct),
        "models": {n: crispri.bench_metrics(v, hct, weighted=True) for n, v in s_all.items()},
        "deletion_gain": crispri.weighted_gain(s_all[with_d], s_all[without], hct, weighted=True),
    }
    answered = [p for p in covered if p.features.get("deletion_answered")]
    versions: Counter[str] = Counter()
    for e in elements.values():
        p = answer_path(e["chrom"], e["id"])
        if p.exists():
            versions[model_version_of(json.loads(p.read_text()))] += 1
    return {
        "registered": crispri.PREREGISTERED_PUBLISHED["second_cell_type"],
        "estimator": "crispri.gain_interval: average precision on the covered HCT116 pairs, 'activity + "
        "distance + deletion' minus 'activity + distance', 95% interval from 200 resamples of whole "
        "chromosomes (seed 0); the headline estimator of the K562 claim this arm replicates",
        "verdict": verdict,
        "passes": passed,
        "replicated": replicated,
        "covered_pairs": len(covered),
        "regulated": sum(lab),
        "chromosomes": len({p.chrom for p in covered}),
        "models": {n: crispri.metrics(v, lab) for n, v in s.items()},
        "deletion_gain": registered,
        "pairs_with_a_deletion_answer": len(answered),
        "regulated_with_a_deletion_answer": sum(p.regulated for p in answered),
        "pairs_with_a_nonzero_drop": sum(p.features["deletion_drop"] > 0 for p in covered),
        "regulated_with_a_nonzero_drop": sum(p.features["deletion_drop"] > 0 for p in covered if p.regulated),
        "gene_is_top_target": sum(bool(p.features["top_target"]) for p in covered),
        "frozen_weights_check_k562_covered": {
            "gain": check,
            "expected": "+0.14 on 1,744 pairs, 114 regulated (the K562 headline), from the same weights",
            "pairs": len(k562),
        },
        "per_cell_block_weighted_all_pairs": per_cell,
        "hct116_answers_by_model_version": dict(versions),
        "sweep_answers_behind_the_weights": sweep_answers,
        "requests": spent(ledger()),
    }


def manifest(training: list[crispri.Pair], heldout: list[crispri.Pair], versions: dict[str, int]) -> dict:
    """What this arm read; the model dependency counts every answer it reads by requested version."""
    from genomeos.attribution.measured import CRISPRI_SPLIT_OF

    inputs = [
        mf.input_entry(crispri.KNOWLEDGE / name, partition=CRISPRI_SPLIT_OF[name], pairs=len(pairs))
        for name, pairs in ((crispri.TRAINING, training), (crispri.HELDOUT, heldout))
    ]
    inputs.append(mf.input_entry(crispri.ELEMENTS, partition=None, role="DeletionTable: elements by overlap"))
    inputs.append(
        mf.input_entry(crispri.ELEMENT_CACHE, partition=None, role="ElementCache: the frozen weights' values")
    )
    inputs.append(
        mf.input_entry(ELEMENT_CACHE_HCT116, partition=None, role="ElementCache: HCT116 deletion values")
    )
    return {
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "cells_kept_per_gene": list(KEPT_CELLS),
            "estimator": "crispri.gain_interval (average precision, chromosome bootstrap)",
            "bootstraps": crispri.BOOTSTRAPS,
            "bootstrap_seed": 0,
            "fit": "crispri.logistic_fit on the covered K562 training pairs, unchanged",
            "threshold": 0.0,
            "request_cap": CAP,
        },
        "partitions": {
            crispri.TRAINING: CRISPRI_SPLIT_OF[crispri.TRAINING],
            crispri.HELDOUT: CRISPRI_SPLIT_OF[crispri.HELDOUT],
            "heldout": "evaluation only: the HCT116 pairs, frozen weights",
        },
        # the weights (and top_target) read the sweep's unrequested answers; HCT116's deletion values
        # read the ALL_FOLDS answers: a result reading both says "mixed"
        "model_dependencies": [mf.answers_model_dependency("alphagenome", versions)],
    }


def write_additively(block: dict[str, Any]) -> None:
    """Insert KEY before INSERT_BEFORE in the committed file, so no existing line changes."""
    text = RESULT.read_text()
    if f'\n  "{KEY}": ' in text:
        raise SystemExit(f"{KEY} is already in {RESULT}; not overwritten")
    i = text.index(f'\n  "{INSERT_BEFORE}": ')
    body = json.dumps({KEY: block}, indent=2)[2:-2]  # the key and its value, at the top level's indent
    new = text[: i + 1] + body + ",\n" + text[i + 1 :]
    json.loads(new)
    RESULT.write_text(new)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--fetch", action="store_true")
    ap.add_argument("--score", action="store_true")
    ap.add_argument("--cap", type=int, default=CAP)
    ap.add_argument("--workers", type=int, default=WORKERS)
    ap.add_argument("--limit", type=int)
    args = ap.parse_args()
    training, heldout = crispri.load(crispri.TRAINING), crispri.load(crispri.HELDOUT)
    elements = to_score(heldout, crispri.DeletionTable())
    have = sum(answer_path(e["chrom"], e["id"]).exists() for e in elements.values())
    print(f"{len(elements)} elements to score ({COSTED} costed), {have} answered, ledger {spent(ledger())}")
    if len(elements) > CAP:
        print(f"refused: {len(elements)} is above the cap {CAP}")
        return 2
    if args.dry_run:
        return 0
    if args.fetch:
        return fetch(elements, min(args.cap, CAP), args.workers, args.limit)
    if args.score:
        t0 = time.time()
        r = score(training, heldout, elements)
        r["versus_the_sweep"] = compare_with_sweep(elements)
        versions = {**r["hct116_answers_by_model_version"]}
        versions[mf.MODEL_VERSION_UNREQUESTED] = (
            versions.get(mf.MODEL_VERSION_UNREQUESTED, 0) + r["sweep_answers_behind_the_weights"]
        )
        r["result_manifest"] = manifest(training, heldout, versions)
        r["date"] = time.strftime("%Y-%m-%d")
        write_additively(r)
        print(
            json.dumps({k: v for k, v in r.items() if k not in ("registered", "result_manifest")}, indent=1)
        )
        print(f"({time.time() - t0:.0f} s) -> {RESULT} [{KEY}]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
