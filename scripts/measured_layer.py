# SPDX-License-Identifier: AGPL-3.0-or-later
"""The experimental layer of the compiled genome, measured and counted. No requests, no model.

    uv run python scripts/measured_layer.py [--chroms chr21,chr22] [--no-save] [--write-programs]

The compiled per-chromosome programs state 940,803 facts and not one of them is `experimental`.
This attaches the measurements the project already holds - CRISPRi enhancer-gene screens, ENCODE4
lentiMPRA, VISTA transgenic assays - to the compiled elements they were made over, under the overlap
rule in `attribution/measured.py`, and reports:

- the **census**: how many compiled facts a measurement now sits beside, per chromosome and pooled,
  as a count and a share, with `coverage_elements_measured` as the named coverage number;
- the **agreement and the disagreement**, separately, because the second is the interesting one, and
  with the denominator that is *not* the subset the prediction chose;
- the **overlap rule at three values**, and the containments it refuses to upgrade on;
- two standardised comparisons through `genomeos.compare`, each printing length, GC and median
  distance to a coding TSS for both arms, so a biased sample is a number and not an oversight.

`--write-programs` recompiles the named chromosomes into `data/knowledge/compiled/` (git-ignored)
and refreshes the committed `data/organisms/human/noncoding_chr21.bio`, which is the copy `bio test`
and the Evidence explorer read without being asked for the compiled tree.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import Any

from genomeos.attribution import measured
from genomeos.attribution.compile import write_program
from genomeos.attribution.targets import attributed
from genomeos.compare import Strata, standardised
from genomeos.evidence import COMPILED_DIR, KINDS, WEAK, collect
from genomeos.results import save_result

CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]
# a prediction "moves a gene" at this magnitude; the bar both arms of every comparison are scored on
MOVES = 0.2
STRATA = Strata(length=(200.0, 400.0), gc=(0.35, 0.45, 0.55), nearest_coding_tss=(1e3, 5e3, 2e4, 1e5))


def covariates(chrom: str, elements: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Length, GC and distance to the nearest coding TSS per element, from the scoring lane's reader."""
    from genomeos.attribution.unknown_scoring import annotate

    return {r["id"]: r for r in annotate(chrom, elements)}


def arm(elements: list[dict[str, Any]], cov: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    """One comparison arm: the covariates the strata need and the hit both arms are scored on."""
    out = []
    for e in elements:
        c = cov.get(e["id"])
        if not c or c.get("gc") is None or c.get("nearest_coding_tss") is None:
            continue
        pc = e.get("predicted_coding") or {}
        out.append(
            {
                "id": e["id"],
                "length": float(c["length"]),
                "gc": float(c["gc"]),
                "nearest_coding_tss": float(c["nearest_coding_tss"]),
                "predicted_moves": abs(float(pc.get("log2_fold_change") or 0.0)) >= MOVES,
            }
        )
    return out


def comparisons(chrom: str, elements: list[dict[str, Any]], rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Is the measured slice a biased sample, and are the disagreements where the model shouts?"""
    cov = covariates(chrom, elements)
    seen = {r["id"] for r in rows}
    disagreed = {r["id"] for r in rows if r["agreement"]["verdict"] == measured.DISAGREES}
    agreed = {r["id"] for r in rows if r["agreement"]["verdict"] == measured.AGREES}
    out: dict[str, Any] = {"hit": f"the deletion moves the target by at least {MOVES} log2"}
    out["measured_vs_unmeasured"] = standardised(
        arm([e for e in elements if e["id"] in seen], cov),
        arm([e for e in elements if e["id"] not in seen], cov),
        STRATA,
        hit="predicted_moves",
    )
    out["disagrees_vs_agrees"] = standardised(
        arm([e for e in elements if e["id"] in disagreed], cov),
        arm([e for e in elements if e["id"] in agreed], cov),
        STRATA,
        hit="predicted_moves",
    )
    return out


def read_chromosome(chrom: str, with_comparison: bool = True) -> dict[str, Any]:
    elements = attributed(chrom)
    layer = measured.Layer.load(chrom)
    rows = measured.rows(chrom, elements, layer)
    out: dict[str, Any] = {
        "census": measured.census(chrom, rows, elements, layer),
        "sensitivity": measured.sensitivity(elements, layer),
    }
    if with_comparison and rows:
        try:
            out["comparison"] = comparisons(chrom, elements, rows)
        except (OSError, KeyError, ValueError) as exc:
            out["comparison"] = {"refused": f"{type(exc).__name__}: {exc}"[:160]}
    return out


def read_evidence(root: Path = Path(".")) -> dict[str, Any]:
    """The Evidence explorer over the compiled programs, per chromosome and pooled.

    `compiled=True` is not optional here and the reason is a defect this project has been bitten by:
    `evidence.collect` skips the compiled tree *before* parsing it, so a census that forgot the flag
    would report zero experimental facts and the zero would read as a finding rather than as a bug.
    `compiled_facts_read` is therefore in the result beside every count, so a zero from reading
    nothing and a zero from measuring nothing can never look the same.
    """
    got = collect(root, compiled=True)
    rows = [r for r in got["rows"] if r["path"].startswith(str(COMPILED_DIR))]
    per: dict[str, dict[str, int]] = {}
    for r in rows:
        chrom = Path(r["path"]).stem.replace("noncoding_", "")
        d = per.setdefault(chrom, {"facts": 0, "weak": 0, **dict.fromkeys(KINDS, 0)})
        d["facts"] += 1
        d[r["evidence"]] = d.get(r["evidence"], 0) + 1
        d["weak"] += r["confidence"] <= WEAK
    pooled = {"facts": len(rows), "weak": sum(1 for r in rows if r["confidence"] <= WEAK)}
    for k in KINDS:
        pooled[k] = sum(1 for r in rows if r["evidence"] == k)
    return {
        "compiled_facts_read": len(rows),
        "programs_read": sorted(per),
        "per_chromosome": per,
        "pooled": pooled,
        "experimental_share": round(pooled["experimental"] / len(rows), 5) if rows else None,
    }


def print_chromosome(chrom: str, got: dict[str, Any]) -> None:
    c = got["census"]
    # the three denominators in order: what there is, what could have been found, what was
    print(
        f"{chrom}: of {c['compiled_elements']:,} compiled elements, "
        f"{c['elements_in_an_assay_footprint']:,} were eligible (an assay's footprint touched them) "
        f"and {c['elements_with_any_measurement']:,} were raised "
        f"({c['coverage_elements_measured'] * 100:.2f}% of all, "
        f"{(c['raised_share_of_the_eligible'] or 0) * 100:.1f}% of the eligible); "
        f"{c['experimental_facts_added']:,} experimental facts "
        f"({c['experimental_element_blocks']:,} elements + {c['experimental_rule_blocks']:,} rules); "
        f"agree {c['agrees']}, disagree {c['disagrees']}",
        flush=True,
    )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chroms", default="", help="comma separated; default every chromosome")
    ap.add_argument("--no-save", action="store_true")
    ap.add_argument("--no-comparison", action="store_true", help="skip the standardised comparisons")
    ap.add_argument("--write-programs", action="store_true", help="recompile the programs with the layer")
    ap.add_argument("--no-evidence", action="store_true", help="skip the 30 s read over the compiled tree")
    args = ap.parse_args(argv)
    chroms = [c.strip() for c in args.chroms.split(",") if c.strip()] or CHROMS

    t0 = time.time()
    per: dict[str, Any] = {}
    refused: dict[str, str] = {}
    for chrom in chroms:
        try:
            per[chrom] = read_chromosome(chrom, not args.no_comparison)
        except (FileNotFoundError, KeyError, ValueError) as exc:
            refused[chrom] = f"{type(exc).__name__}: {exc}"[:160]
            print(f"{chrom}: not read ({refused[chrom]})", flush=True)
            continue
        print_chromosome(chrom, per[chrom])

    pooled = measured.pool([v["census"] for v in per.values()])
    written: dict[str, str] = {}
    if args.write_programs:
        (Path(".") / COMPILED_DIR).mkdir(parents=True, exist_ok=True)
        for chrom in per:
            p = write_program(chrom, Path(COMPILED_DIR) / f"noncoding_{chrom}.bio")
            written[chrom] = str(p)
            if chrom == "chr21":  # the copy the Evidence explorer and bio test read by default
                written["chr21_committed"] = str(write_program(chrom))
            print(f"{chrom}: program written to {written[chrom]}", flush=True)

    evidence: dict[str, Any] = {}
    if not args.no_evidence:
        evidence = read_evidence()
        p = evidence["pooled"]
        print(
            f"\nevidence over {evidence['compiled_facts_read']:,} compiled facts read from "
            f"{len(evidence['programs_read'])} programs: experimental {p['experimental']:,}, "
            f"curated {p['curated']:,}, predicted {p['predicted']:,}, inferred {p['inferred']:,}"
        )

    payload = {
        "rule": {
            "name": "reciprocal overlap",
            "value": measured.RECIPROCAL_OVERLAP,
            "reported_at": list(measured.OVERLAP_SENSITIVITY),
            "statement": (
                "a measurement counts as a measurement of a compiled element only when the overlap "
                "covers at least this fraction of the element and of the tested interval; proximity, "
                "similarity and mere containment never raise a fact"
            ),
        },
        "denominators": (
            "read in this order, because a small census is a statement about the search and not about "
            "the genome: compiled_facts_read (what was parsed at all), "
            "eligibility.elements_in_an_assay_footprint (what any assay could have measured, at one "
            "shared base), then elements_with_any_measurement, agrees and disagrees"
        ),
        "per_chromosome": per,
        "pooled": pooled,
        "evidence": evidence,
        "chromosomes_not_read": refused,
        "programs_written": written,
        "reading": (
            "the experimental layer is thin and that is the finding: a compiled element is measured "
            "only where an assay happened to be pointed, and the assays cover a few per cent of what "
            "the genome is compiled from. lentiMPRA measures the sequence in a reporter and not the "
            "locus, so its disagreements are weaker evidence against a prediction than a CRISPRi "
            "negative on the predicted gene itself"
        ),
        "seconds": round(time.time() - t0, 1),
    }
    el = pooled["eligibility"]
    print(
        f"\npooled: of {pooled['compiled_elements']:,} compiled elements, "
        f"{el['elements_in_an_assay_footprint']:,} were eligible for a measurement at all "
        f"({el['share_eligible'] * 100:.2f}%) and {el['elements_no_assay_ever_covered']:,} were never "
        f"covered by any assay"
    )
    print(
        f"  raised: {pooled['elements_with_any_measurement']:,} "
        f"({pooled['coverage_elements_measured'] * 100:.2f}% of all, "
        f"{(pooled['raised_share_of_the_eligible'] or 0) * 100:.1f}% of the eligible), "
        f"{pooled['experimental_facts_added']:,} experimental facts beside "
        f"{pooled['facts_with_a_measured_counterpart']:,} predicted ones"
    )
    print(
        f"agreement: {pooled['agrees']} agree, {pooled['disagrees']} disagree; "
        f"CRISPRi rate over all matched elements "
        f"{pooled['agreement_rate_over_all_matched_elements']} "
        f"(n={pooled['agreement_denominator_all_matched_elements']}), over the pairs where the "
        f"predicted gene was tested {pooled['agreement_rate_where_the_predicted_gene_was_tested']} "
        f"(n={pooled['agreement_denominator_predicted_gene_tested']})"
    )
    if not args.no_save:
        p = save_result("measured_layer_genome", payload)
        print(f"saved {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
