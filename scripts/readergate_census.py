# SPDX-License-Identifier: AGPL-3.0-or-later
"""M2b census: what the compiler's reader layer emits over the 24 programs, and what can be run.

    uv run python scripts/readergate_census.py [--chroms chr21] [--no-save] [--emit DIR]

Gates each of the 24 generated programs M1 ran (data/knowledge/compiled, pinned by sha256) with
`compile.gate_program(selection='reader_cells')` in memory, parses the gated text, and counts per
reader cell type how many elements the reader states open, closed and unknown, how many rules the
reader holds, and the sha256 of every gated program. Then the single-rule genes (POPULATION (a)) are
bridged and run in all thirteen cells by `readergate.run_chromosome`, with the closed-form and gate
checks the registered run will apply. No RNA is opened, no model request is sent, nothing is fetched.

The gated programs are NOT written into the generated store: scripts there (when_census,
gene_identity, holdout, evidence, the rule-evidence-tier cost census) glob data/knowledge/compiled and
data/**/*.bio, and a second copy of every rule would be counted by all of them. `--emit DIR` writes
them to a directory outside data/ for anyone who wants to `bio check` one; it refuses a path inside.
"""

from __future__ import annotations

import argparse
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import compile as cp  # noqa: E402
from genomeos.attribution import readergate as rg  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readergate_fetch import OWN_CODE  # noqa: E402
from runrules_census import RATE_TABLE, load_rates  # noqa: E402

RESULT = "readergate_census"
CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]
ROOT = Path(__file__).resolve().parents[1]


def peak_files(chroms: list[str]) -> list[Path]:
    from genomeos.attribution import context_evidence as ce
    from genomeos.genome import reader

    return [RESULTS_DIR / reader.peaks_path(b, c).name for b in ce.READER_TERMS for c in chroms]


def census(chroms: list[str], emit: Path | None = None) -> dict:
    rates, borrowed, prov, lent = load_rates()
    texts = {c: (Path(COMPILED_DIR) / f"noncoding_{c}.bio").read_text() for c in chroms}
    singles, why = rg.single_rule_genes(texts)
    programs, genes, counts = {}, [], Counter()
    held: dict[str, Counter] = {cid: Counter() for cid in rg.CELLS}
    states: dict[str, Counter] = {cid: Counter() for cid in rg.CELLS}
    for chrom in chroms:
        t = time.time()
        info, rows, c, module = rg.run_chromosome(chrom, texts[chrom], singles, rates, borrowed, prov, lent)
        for cid in rg.CELLS:
            states[cid].update(info["states_by_cell"][cid])
            held[cid].update(s for _r, s, _w in module.reader_gate({"cell_type": cid}))
        if emit is not None:
            gated, again = cp.gate_program(texts[chrom], chrom, selection=rg.SELECTION)
            assert again["gated_sha256"] == info["gated_sha256"]
            (emit / f"noncoding_{chrom}_reader.bio").write_text(gated)
        programs[chrom] = {k: v for k, v in info.items() if k != "states_by_cell"}
        counts.update(c)
        genes += rows
        print(chrom, info["elements"], dict(c), round(time.time() - t, 1), flush=True)
    bad_cf = [g["gene"] for g in genes if abs(g["m"] - g["closed_form"]) > rg.CLOSED_FORM_TOLERANCE]
    bad_gate = [g["gene"] for g in genes if g["gate_mismatch"]]
    opened = Counter(sum(1 for s in g["states"].values() if s == "open") for g in genes)
    informative = [
        g for g in genes if 0 < sum(1 for s in g["states"].values() if s == "open") < len(rg.CELLS)
    ]
    ids = rg.symbol_gene_ids()
    return {
        "lane": rg.LANE,
        "date": time.strftime("%Y-%m-%d"),
        "question": "what the reader layer emits, and which genes a gated run can run (no RNA read)",
        "selection": rg.SELECTION,
        "structure": rg.STRUCTURE,
        "emitted": {
            "programs": len(programs),
            "elements": sum(p["elements"] for p in programs.values()),
            "rules_gated": dict(sum((Counter(p["rules_gated"]) for p in programs.values()), Counter())),
            "element_states_by_cell": {cid: dict(states[cid]) for cid in rg.CELLS},
            "rules_by_gate_state_by_cell": {cid: dict(held[cid]) for cid in rg.CELLS},
            "per_program": programs,
            "written_to_disk": str(emit) if emit is not None else "no: built in memory (module docstring)",
        },
        "population_census": {
            **dict(why),
            **dict(counts),
            "closed_form_beyond_tolerance": len(bad_cf),
            "gate_mismatch_genes": len(bad_gate),
            "genes_by_open_cell_count": {str(k): v for k, v in sorted(opened.items())},
            "informative_genes": len(informative),
            "informative_by_action": dict(Counter(g["action"] for g in informative)),
            "informative_symbol_on_one_gene_id": sum(1 for g in informative if g["symbol"] in ids),
            "informative_compiled_cell_is_a_reader_cell": sum(
                1 for g in informative if g["compiled_cell"] in rg.CELLS
            ),
            "max_abs_m_minus_closed_form": max((abs(g["m"] - g["closed_form"]) for g in genes), default=None),
            "rna_read": "none",
        },
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }


def manifest(chroms: list[str]) -> dict:
    progs = [mf.input_entry(Path(COMPILED_DIR) / f"noncoding_{c}.bio", partition=None) for c in chroms]
    inputs = [
        *progs,
        mf.input_entry(RATE_TABLE, partition=None),
        mf.input_entry(rg.GENE_TABLE, partition=None),
        mf.files_entry("reader v1 DNase narrowPeak files", [p for p in peak_files(chroms) if p.exists()]),
    ]
    return {
        "sources": [
            {"accession": "the 24 compiled non-coding programs M1 ran", "version": "pinned by sha256"},
            {"accession": "reader v1 ENCODE DNase-seq narrowPeak files, 13 biosamples", "version": "pinned"},
            {"accession": "data/results/gene_rates.json (mouse NIH 3T3, e921fae)", "version": "pinned"},
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38",
        "coordinates": "element loci as the compiled programs write them, 0-based half-open",
        "parameters": {"selection": rg.SELECTION, "runtime": rg.RUNTIME, "reader_state": cp.READER_STATE},
        "exclusions": ["genes with several rules or one measured rule", "genes the bridge does not resolve"],
        "partitions": "n/a: a census; no outcome read",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", nargs="*", default=CHROMS)
    ap.add_argument("--no-save", action="store_true")
    ap.add_argument("--emit", type=Path, default=None, help="write the gated programs here (outside data/)")
    args = ap.parse_args()
    if args.emit is not None:
        emit = args.emit.resolve()
        if emit == ROOT / "data" or (ROOT / "data") in emit.parents:
            raise SystemExit(f"refused: {emit} is inside data/, where other censuses would count it")
        emit.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    out = census(args.chroms, args.emit)
    out["seconds"] = round(time.time() - t0, 1)
    print(out["population_census"])
    if not args.no_save:
        print("saved", save_result(RESULT, out, manifest=manifest(args.chroms)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
