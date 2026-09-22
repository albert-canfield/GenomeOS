# SPDX-License-Identifier: AGPL-3.0-or-later
"""Does §17's motif-count positive transfer off the episomal reporter? Two readings, both pre-registered.

    uv run python scripts/motif_transfer.py --stage sites            # cache the VISTA site counts
    uv run python scripts/motif_transfer.py --stage preregister      # training numbers only, no held-out
    uv run python scripts/motif_transfer.py --stage score            # the held-out sets, scored once

Section 17 of docs/GRAMMAR-BY-COMPARISON.md failed the arrangement claim and left one large positive:
strict family-collapsed JASPAR site counts lift Spearman with lentiMPRA activity by +0.234 (K562) over
composition alone. This asks whether that survives a different assay (VISTA's transgenic mouse embryos)
and a different question (the 882 constrained-unknown blocks of the real unknown), or whether it is a
property of a 200 bp plasmid.

Stage order matters. `sites` writes the cached scan under data/knowledge/motif_transfer (git-ignored),
`preregister` writes `motif_transfer_preregistration` with the training-fold numbers and the
pre-registered claims and reads nothing from the held-out chromosomes, and `score` writes
`motif_transfer`. The pre-registration is committed before `score` runs.

No network: VISTA's loci and annotation, the lentiMPRA rows and site scan, the JASPAR matrices, the
deletion sweep's element tables and the reference FASTAs are all on disk already.
"""

from __future__ import annotations

import argparse
import time

from genomeos.attribution import motif_transfer as mt


def cache_sites(chroms: list[str]) -> int:
    total = 0
    for chrom in chroms:
        t0 = time.time()
        n = mt.build_vista_counts(chrom)
        total += n
        print(f"{chrom}: {n} VISTA elements scanned in {time.time() - t0:.1f}s", flush=True)
    return total


def report_vista(out: dict) -> None:
    v = out.get("vista")
    if not v:
        return
    print(f"\nVISTA: {v['training_elements']} training, {v['held_out_elements']} held out")
    print(f"  coverage: {v['coverage']}")
    cov = v["covariates_held_out"]["positives_above_negatives_auroc"]
    print(f"  unmatched covariates, positives above negatives (AUROC): {cov}")
    held = v.get("held_out")
    if not held:
        print(
            "  inner folds only: "
            + ", ".join(f"{m} {r['inner_auroc']}" for m, r in v["inner_cross_validation"].items())
        )
        return
    for name in v["models"]:
        r = held[name]
        print(f"  ({name}) AUROC {r['auroc']} {r['ci95']}  AUPRC {r.get('auprc')}")
    for key in ("counts - composition", "counts - counts_shuffled", "both - conservation", "both - counts"):
        if key in held:
            print(f"  {key}: {held[key]['auroc']} {held[key]['ci95']}")
    print(f"  verdict: {v['verdict']['statement']}")
    for g, r in (v.get("tissue") or {}).get("by_group", {}).items():
        print(f"  tissue {g}: {r.get('auroc')} {r.get('ci95')} {r.get('skipped', '')}")


def report_blocks(out: dict) -> None:
    b = out.get("real_unknown")
    if not b:
        return
    if "measured" not in b:  # the pre-registration stage: coverage only, nothing scored
        print("\nreal unknown, coverage before scoring:")
        for tier, rec in b["coverage"].items():
            print(f"  {tier}: {rec}")
        return
    m = b["measured"]
    print(f"\nreal unknown: {b['blocks']} blocks, {b['windows']}")
    print(f"  measured coverage: {m['coverage']}")
    for cell, r in (m.get("level") or {}).items():
        print(f"  measured level {cell}: AUROC {r['auroc']} z {r['z']} n {r['n']} medians {r['median']}")
    for cell, r in (m.get("transfer") or {}).items():
        d = r["b-a"]
        print(
            f"  measured transfer {cell}: a {r['a']['point']} b {r['b']['point']} "
            f"b-a {d['point']} {d['ci95']}"
        )
    for cell in ("K562", "HepG2", "WTC11"):
        r = (b["predicted_level"] or {}).get(cell)
        if r:
            print(f"  predicted level {cell}: vs matched neutral {r['real_unknown_vs_matched_neutral']}")
    a = b["agreement_with_the_sweep"].get("tested_blocks_only")
    if a:
        print(f"  sweep agreement (two model readings): {a['count_model_mean']}")
        print(
            f"    length alone {a['block_length']['auroc']}, tested elements {a['tested_elements']['auroc']}"
        )
    print(f"  verdict: {b['verdict']['measured_statement']}")
    print(f"  verdict: {b['verdict']['transfer_statement']}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--stage", choices=("sites", "preregister", "score"), default="score")
    ap.add_argument("--chroms", default="", help="comma-separated, for --stage sites; default all")
    ap.add_argument("--readings", default="vista,blocks", help="comma-separated: vista, blocks")
    ap.add_argument("--rescan", action="store_true", help="rescan the block windows instead of reusing")
    ap.add_argument("--no-save", action="store_true", help="print without writing the result")
    ap.add_argument("--name", default="", help="write under this result name instead of the default")
    args = ap.parse_args(argv)
    t0 = time.time()
    if args.stage == "sites":
        chroms = args.chroms.split(",") if args.chroms else list(mt.CHROMS)
        n = cache_sites([c for c in chroms if (mt.VISTA_KNOWLEDGE / f"rows_{c}.json").exists()])
        print(f"\n{n} VISTA elements cached under {mt.CACHE} in {time.time() - t0:.0f}s")
        return 0
    readings = tuple(r for r in args.readings.split(",") if r)
    score = args.stage == "score"
    out = (
        mt.run(readings=readings, score_held_out=score, reuse=not args.rescan, progress=print)
        if args.no_save
        else mt.run_and_save(
            readings=readings,
            score_held_out=score,
            reuse=not args.rescan,
            name=args.name or None,
            progress=print,
        )
    )
    report_vista(out)
    report_blocks(out)
    print(f"\n{'not saved' if args.no_save else 'saved'} in {time.time() - t0:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
