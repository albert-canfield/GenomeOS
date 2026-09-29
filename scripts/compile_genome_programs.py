# SPDX-License-Identifier: AGPL-3.0-or-later
"""Compile the non-coding genome to BioLang, one program per chromosome, and read its evidence.

    uv run python scripts/compile_genome_programs.py [--chroms chr21,chr22] [--no-save]

Area I compiles a chromosome's non-coding space into a program whose every region carries a role, an
evidence kind and a confidence (`attribution/compile.py`). Only chr21 was ever compiled, because the
inputs for the rest were not finished; both genome-wide sweeps closed on 2026-09-16 and 2026-09-17,
so the whole genome can be compiled now.

The programs go to `data/knowledge/compiled/`, git-ignored: 24 chromosomes are about 120 MB of
generated text, and this project commits summaries and rebuilds text. `genomeos.evidence` reads that
directory when it exists, so the Evidence explorer covers them without another switch.

What it is for. The Evidence explorer has been reading the demos, the organism programs and the
prelude — which are curated and hand-written, and therefore mostly strong. The compiled chromosomes
are where the weak half sits: a region whose role came from a tier is `inferred`, an element's target
from a model is `predicted`, and a program that is 99% predicted says so only if someone counts. This
counts, per chromosome and pooled, and writes `evidence_compiled_genome`.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

from genomeos.attribution.compile import compile_chromosome
from genomeos.evidence import COMPILED_DIR, KINDS, WEAK, collect
from genomeos.results import save_result

CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]


def compile_all(chroms: list[str], root: Path = Path(".")) -> dict[str, dict]:
    """Write one program per chromosome, reporting the ones whose inputs are not there."""
    out: dict[str, dict] = {}
    (root / COMPILED_DIR).mkdir(parents=True, exist_ok=True)
    for chrom in chroms:
        t0 = time.time()
        try:
            text = compile_chromosome(chrom)
        except (FileNotFoundError, KeyError, ValueError) as exc:
            out[chrom] = {"written": False, "why": f"{type(exc).__name__}: {exc}"[:160]}
            print(f"{chrom}: not compiled ({out[chrom]['why']})", flush=True)
            continue
        path = root / COMPILED_DIR / f"noncoding_{chrom}.bio"
        path.write_text(text)
        out[chrom] = {
            "written": True,
            "bytes": len(text),
            "lines": text.count("\n"),
            "seconds": round(time.time() - t0, 2),
        }
        print(f"{chrom}: {len(text) / 1e6:.1f} MB, {text.count(chr(10)):,} lines", flush=True)
    return out


STATED_COUNTS = ("stated", "unstated", "weak", "strong", "weak_before_r4_counting_unstated")
WEAK_MEANS = (
    f"a STATED confidence at or below {WEAK}; an unstated one (the parser's UNSTATED, since R4 every "
    "predicted fact) is counted apart. weak_before_r4_counting_unstated is the old count, which "
    "compared the unstated 0.0 with the line, kept beside it for comparison"
)


def count_stated(tally: dict[str, int], row: dict) -> None:
    """Add one evidence row to stated/unstated and stated-weak/stated-strong (evidence.py's `stated`)."""
    stated = row.get("stated", True)
    low = row["confidence"] <= WEAK
    tally["stated"] += stated
    tally["unstated"] += not stated
    tally["weak"] += stated and low
    tally["strong"] += stated and not low
    tally["weak_before_r4_counting_unstated"] += low


def read_evidence(root: Path = Path(".")) -> dict:
    """The evidence reading over the compiled programs alone, per chromosome and pooled."""
    got = collect(root, compiled=True)
    rows = [r for r in got["rows"] if r["path"].startswith(str(COMPILED_DIR))]
    per_chrom: dict[str, dict] = {}
    pooled = {"facts": 0, **dict.fromkeys(STATED_COUNTS, 0), **dict.fromkeys(KINDS, 0)}
    for r in rows:
        chrom = Path(r["path"]).stem.replace("noncoding_", "")
        d = per_chrom.setdefault(
            chrom, {"facts": 0, **dict.fromkeys(STATED_COUNTS, 0), **dict.fromkeys(KINDS, 0)}
        )
        for tally in (d, pooled):
            tally["facts"] += 1
            tally[r["evidence"]] = tally.get(r["evidence"], 0) + 1
            count_stated(tally, r)
    return {"weak_means": WEAK_MEANS, "per_chromosome": per_chrom, "pooled": pooled}


def root_program(chrom: str, root: Path = Path(".")) -> Path:
    return root / COMPILED_DIR / f"noncoding_{chrom}.bio"


def manifest(chroms: list[str], root: Path = Path(".")) -> dict:
    """The provenance contract: the compiled programs read, by sha256."""
    from genomeos import manifest as mf

    return {
        "sources": [
            {
                "accession": "this repository, the compiled programs noncoding_<chrom>.bio",
                "version": "pinned by sha256",
            }
        ],
        "inputs": [
            mf.input_entry(root_program(c, root), partition=None)
            for c in chroms
            if root_program(c, root).exists()
        ],
        "assembly": "GRCh38",
        "coordinates": "n/a: counts of stated facts, no interval is read",
        "parameters": {"weak_at_or_below": WEAK, "compiled": True},
        "exclusions": ["facts in programs outside the compiled tree"],
        "partitions": "n/a: no evaluation partition",
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chroms", default="")
    ap.add_argument("--no-save", action="store_true")
    ap.add_argument("--result-name", default="evidence_compiled_genome", help="the name to save under")
    ap.add_argument(
        "--evidence-only", action="store_true", help="read the compiled tree as it is; compile nothing"
    )
    args = ap.parse_args(argv)
    chroms = args.chroms.split(",") if args.chroms else CHROMS

    t0 = time.time()
    if args.evidence_only:
        present = {c: root_program(c) for c in chroms}
        written = {
            c: {"written": False, "why": "not compiled by this run (--evidence-only)"}
            if p.exists()
            else {"written": False, "why": "no program on disk"}
            for c, p in present.items()
        }
    else:
        written = compile_all(chroms)
    print("reading the evidence over the compiled programs", flush=True)
    evidence = read_evidence()
    pooled = evidence["pooled"]
    out = {
        "result": args.result_name,
        "chromosomes_compiled": sorted(c for c, v in written.items() if v["written"]),
        "programs_read": sorted(evidence["per_chromosome"]),
        "chromosomes_not_compiled": {c: v["why"] for c, v in written.items() if not v["written"]},
        "megabytes": round(sum(v.get("bytes", 0) for v in written.values()) / 1e6, 1),
        "evidence": evidence,
        "weak_share": round(pooled["weak"] / pooled["facts"], 4) if pooled["facts"] else None,
        "unstated_share": round(pooled["unstated"] / pooled["facts"], 4) if pooled["facts"] else None,
        "weak_share_of_the_stated": round(pooled["weak"] / pooled["stated"], 4) if pooled["stated"] else None,
        "reading": (
            "the compiled chromosomes are the project's weak half by construction: a region's role "
            "comes from its tier (inferred) and an element's target from a model (predicted), while the "
            "hand-written organism programs are curated. Reported so the Evidence explorer's totals "
            "cannot be read as if the whole project were as well evidenced as its demos"
        ),
        "seconds": round(time.time() - t0, 1),
    }
    if args.no_save:
        print("not saved (--no-save)")
    else:
        print(
            f"saved {save_result(out['result'], out, manifest=manifest(sorted(evidence['per_chromosome'])))}"
        )
    print(f"\n{pooled['facts']:,} facts over {len(out['programs_read'])} programs read")
    for k in KINDS:
        share = pooled[k] / pooled["facts"] if pooled["facts"] else 0
        print(f"  {k:14} {pooled[k]:9,}  {share:6.1%}")
    print(f"  {'weak (<= ' + str(WEAK) + ')':14} {pooled['weak']:9,}  {out['weak_share']:.1%}  stated only")
    print(f"  {'unstated':14} {pooled['unstated']:9,}  {out['unstated_share']:.1%}")
    print(
        f"  {'weak, old rule':14} {pooled['weak_before_r4_counting_unstated']:9,}  unstated counted as weak"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
