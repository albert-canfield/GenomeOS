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


def read_evidence(root: Path = Path(".")) -> dict:
    """The evidence reading over the compiled programs alone, per chromosome and pooled."""
    got = collect(root, compiled=True)
    rows = [r for r in got["rows"] if r["path"].startswith(str(COMPILED_DIR))]
    per_chrom: dict[str, dict] = {}
    for r in rows:
        chrom = Path(r["path"]).stem.replace("noncoding_", "")
        d = per_chrom.setdefault(chrom, {"facts": 0, "weak": 0, **{k: 0 for k in KINDS}})
        d["facts"] += 1
        d[r["evidence"]] = d.get(r["evidence"], 0) + 1
        if r["confidence"] <= WEAK:
            d["weak"] += 1
    pooled = {"facts": len(rows), "weak": sum(1 for r in rows if r["confidence"] <= WEAK)}
    for k in KINDS:
        pooled[k] = sum(1 for r in rows if r["evidence"] == k)
    return {"per_chromosome": per_chrom, "pooled": pooled}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chroms", default="")
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args(argv)
    chroms = args.chroms.split(",") if args.chroms else CHROMS

    t0 = time.time()
    written = compile_all(chroms)
    print("reading the evidence over the compiled programs", flush=True)
    evidence = read_evidence()
    pooled = evidence["pooled"]
    out = {
        "result": "evidence_compiled_genome",
        "chromosomes_compiled": sorted(c for c, v in written.items() if v["written"]),
        "chromosomes_not_compiled": {c: v["why"] for c, v in written.items() if not v["written"]},
        "megabytes": round(sum(v.get("bytes", 0) for v in written.values()) / 1e6, 1),
        "evidence": evidence,
        "weak_share": round(pooled["weak"] / pooled["facts"], 4) if pooled["facts"] else None,
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
        print(f"saved {save_result(out['result'], out)}")
    print(f"\n{pooled['facts']:,} facts over {len(out['chromosomes_compiled'])} chromosomes")
    for k in KINDS:
        share = pooled[k] / pooled["facts"] if pooled["facts"] else 0
        print(f"  {k:14} {pooled[k]:9,}  {share:6.1%}")
    print(f"  {'weak (<= ' + str(WEAK) + ')':14} {pooled['weak']:9,}  {out['weak_share']:.1%}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
