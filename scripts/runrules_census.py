# SPDX-License-Identifier: AGPL-3.0-or-later
"""M1 census: which compiled rules can run, per cell context, and why the rest cannot.

    uv run python scripts/runrules_census.py [--chroms chr21] [--no-save]

Reads the 24 compiled programs and the rate table only. No measured RNA is opened, no model request
is sent, nothing is downloaded. For every cell context a compiled rule is gated on it reports the
rules, mechanisms and gene-cell pairs, how many pairs `attribution.bridge.parameterize` resolves on
the measured tier (the gene's own mouse rate) and the borrowed tier (relative units), and every
unresolved reason by name; then the same pooled over the 54 contexts named for GTEx tissues, which
are the panel the registered test (genomeos/attribution/runrules.py) is run on. It also counts the
compiled rules whose source is a gene, which decides whether a run is a network at all.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import bridge  # noqa: E402
from genomeos.attribution import runrules as rr  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.lang.parser import parse  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

RESULT = "runrules_census"
CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]
RATE_TABLE = RESULTS_DIR / "gene_rates.json"
OWN_CODE = (
    "genomeos/attribution/runrules.py",
    "scripts/runrules_census.py",
    "scripts/runrules_run.py",
    "scripts/runrules_register.py",
    "tests/test_runrules.py",
)


def load_rates() -> tuple[dict[str, float], float, dict[str, dict], dict[str, str]]:
    """The rate table as `scripts/bridge_audit.py` reads it: measured T, the median, provenance."""
    d = json.loads(RATE_TABLE.read_text())
    rates = {g: float(r["transcription_rate"]) for g, r in d["rates"].items()}
    prov = {
        g: {
            "source": str(r.get("source", "unstated")),
            "species_cell": str(r.get("species_cell", "unstated")),
        }
        for g, r in d["rates"].items()
    }
    sources = sorted({v["source"] for v in prov.values()})
    cells = sorted({v["species_cell"] for v in prov.values()})
    lent = {"source": "genome median of " + ", ".join(sources), "species_cell": "; ".join(cells)}
    return rates, float(d["borrowed_median_transcription_rate"]), prov, lent


def manifest(paths: list[Path]) -> dict:
    inputs = [mf.input_entry(p, partition=None) for p in paths] + [mf.input_entry(RATE_TABLE, partition=None)]
    return {
        "sources": [
            {
                "accession": "the compiled non-coding programs, scripts/compile_genome_programs.py",
                "version": "the files on this machine, pinned by sha256 in inputs",
            },
            {
                "accession": "data/results/gene_rates.json (lane-rates, e921fae): mouse NIH 3T3 rates",
                "version": "the committed file, pinned by sha256 in inputs",
            },
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38, the assembly the compiled programs were written on; no interval is read",
        "coordinates": "n/a: counts of rules, mechanisms and gene-cell pairs per cell context",
        "parameters": {
            "tiers": {
                "measured": bridge.RATE_TIERS["measured"],
                "borrowed": bridge.RATE_TIERS["borrowed_median"],
            },
            "rate_species": bridge.RATE_SPECIES,
            "panel": "the 54 compiled contexts named for GTEx v8 tissues (runrules.PANEL)",
            "no_measured_rna_read": "the census opens no expression value of any kind",
        },
        "exclusions": [],
        "partitions": "n/a: every compiled rule is counted, held-out or not",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chroms", default="")
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args(argv)
    chroms = args.chroms.split(",") if args.chroms else CHROMS
    t0 = time.time()
    rates, borrowed, prov, lent = load_rates()
    per: dict[str, dict] = {}
    paths = []
    for chrom in chroms:
        path = Path(COMPILED_DIR) / f"noncoding_{chrom}.bio"
        paths.append(path)
        per[chrom] = rr.census_module(parse(path.read_text(), path.stem), rates, borrowed, prov, lent)
        print(chrom, per[chrom]["rules"], len(per[chrom]["contexts"]), flush=True)
    pooled = rr.pool_census(per)
    contexts = pooled["contexts"]
    panel = set(rr.PANEL)
    out = {
        "lane": rr.LANE,
        "question": "which compiled rules can run, per cell context, and why the rest cannot",
        "structure": rr.STRUCTURE,
        "totals": pooled["totals"],
        "all_contexts": rr.summarise(contexts),
        "panel_contexts": rr.summarise(contexts, panel),
        "panel_labels_present": sorted(panel & set(contexts)),
        "panel_labels_absent": sorted(panel - set(contexts)),
        "per_context": contexts,
        "per_chromosome_totals": {
            c: {"rules": v["rules"], "contexts": len(v["contexts"])} for c, v in per.items()
        },
        "reading": (
            "a gene-cell pair is resolved when the bridge fits its one mechanism; the borrowed tier"
            " resolves in relative units only and is never a measured count. Pairs are counted per"
            " program, so a gene whose rules sit in two programs (the X/Y pseudoautosomal genes) counts"
            " once in each; the registered run pools mechanisms over programs before it admits a pair"
        ),
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
        "seconds": round(time.time() - t0, 1),
    }
    if not args.no_save:
        print("saved", save_result(RESULT, out, manifest=manifest(paths)))
    print(json.dumps({k: out[k] for k in ("totals", "all_contexts", "panel_contexts")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
