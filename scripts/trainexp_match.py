# SPDX-License-Identifier: AGPL-3.0-or-later
"""POST HOC: M2i's independence from AlphaGenome's training data, and the eleven-cell subset.

    uv run python scripts/trainexp_match.py

Writes data/results/trainexp_independence.json. Everything in it is labelled post hoc and decides
nothing: it QUALIFIES M2i (`f6e71c6`), it does not revise it. No model request, no key, no money, no
network: the AlphaGenome track metadata is read from the copy already on disk, and if that copy is
absent the run reports its absence rather than fetching anything.

Two things it computes.

1. The matching. M2i's thirteen fresh experiments and items 7/7d/7h's twenty-seven already-read
   quantifications (thirteen experiments) against AlphaGenome's recorded output-track metadata, by
   ENCODE accession where the metadata carries one and otherwise by biosample ontology CURIE plus
   assay title, with the method named per experiment and its strength stated.

2. The eleven-cell subset. M2i's comparison with cardiac muscle cell and hepatocyte dropped -- the
   two whose biosamples an already-read experiment also used -- at M2i's constant scale, with M2i's
   own interval method. NOTHING IS FITTED: `freshrna.fit_scale` raises by design and is never called.
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import freshrna as fr  # noqa: E402
from genomeos.attribution import readergate as rg  # noqa: E402
from genomeos.attribution import trainexp as tx  # noqa: E402
from genomeos.results import save_result  # noqa: E402

GENE_TABLE = Path("data/cache/gencode_genes.tsv")
M2I_RESULT = Path("data/results/freshrna_result.json")
OWN_CODE = (
    "genomeos/attribution/trainexp.py",
    "scripts/trainexp_match.py",
    "tests/test_trainexp.py",
)


def m2i_run() -> Any:
    """M2i's own run script, imported read-only so the subset is built by M2i's code, not a copy."""
    path = Path("scripts/freshrna_run.py")
    spec = importlib.util.spec_from_file_location("freshrna_run_readonly", path)
    if spec is None or spec.loader is None:  # pragma: no cover - the file is in the commit
        raise SystemExit(f"refused: {path} cannot be imported")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def matching() -> dict[str, Any]:
    """Both sets against the recorded track metadata, each match naming its method and strength."""
    rows = tx.read_tracks()
    rna = tx.rna_tracks(rows)
    accessions = tx.accession_index(rows)
    by_ba = tx.biosample_assay_index(rna)
    by_curie = tx.curie_index(rna)
    terms = tx.cell_terms()

    fresh = json.loads(tx.FRESH_SOURCES.read_text())
    read = json.loads(tx.READ_SOURCES.read_text())

    fresh_entries = []
    for cell, entry in fresh["biosamples"].items():
        name, curie = terms[cell]
        fresh_entries.append(
            tx.classify_one(
                cell,
                entry["experiment"],
                name,
                curie,
                fresh["assays"][cell],
                accessions,
                by_ba,
                by_curie,
                files=len(entry["files"]),
            )
        )

    # 7/7d/7h's selection rule admits ONE assay for all thirteen; it is read from that file's rule,
    # not assumed, and the run refuses rather than guessing if the wording ever changes.
    read_assay = "polyA plus RNA-seq"
    if f"assay '{read_assay}'" not in read["selection_rule"]:
        raise SystemExit(
            "refused: readergate_rna_sources.json's selection rule no longer names"
            f" assay '{read_assay}', so the already-read assay cannot be taken from it"
        )
    read_entries = []
    for cell, entry in read["biosamples"].items():
        name, curie = terms[cell]
        read_entries.append(
            tx.classify_one(
                cell,
                entry["experiment"],
                name,
                curie,
                read_assay,
                accessions,
                by_ba,
                by_curie,
                files=len(entry["files"]),
            )
        )

    return {
        "label": "POST HOC (1): the two sets against AlphaGenome's recorded output-track metadata",
        "metadata_file": tx.METADATA.as_posix(),
        "metadata_fields": sorted(rows[0].keys()) if rows else [],
        "metadata_rows": len(rows),
        "rna_seq_rows": len(rna),
        "rna_seq_distinct_track_names": len({r["name"] for r in rna}),
        "rna_seq_rows_carrying_an_accession": sum(
            1 for r in rna if any(tx.ACCESSION.match(f) for f in (r.get("name") or "").split())
        ),
        "accessions_anywhere_in_the_metadata": {k: v for k, v in sorted(accessions.items())},
        "accession_matching_available_for_rna_seq": False,
        "classes": dict(tx.CLASSES),
        "training_set_record": tx.TRAINING_SET_RECORD,
        "fresh_thirteen": {
            "what": "M2i's thirteen fresh experiments (freshrna_sources non_overlap.this_lane)",
            "entries": fresh_entries,
            "tally": tx.tally(fresh_entries),
        },
        "already_read_twenty_seven": {
            "what": "items 7/7d/7h's twenty-seven already-read quantifications, thirteen experiments",
            "assay": read_assay,
            "entries": read_entries,
            "tally": tx.tally(read_entries),
        },
        "sentence_for_every_m2i_better": tx.independence_sentence(fresh_entries),
    }


def eleven_cells(mod: Any) -> dict[str, Any]:
    """M2i's comparison on eleven cells, built with M2i's own code at M2i's constant scale."""
    sources = json.loads(fr.SOURCES.read_text())
    seven_h = mod._seven_h_rows()  # noqa: SLF001 - M2i's committed P rows, read-only
    per_cell, counts = mod.fresh_tpm(sources)
    all_rows, attrition = mod.assemble(seven_h, per_cell)
    primary = [r for r in all_rows if fr.is_held_out(r["symbol"])]

    keep = tx.subset_indices(fr.CELLS, tx.SHARED_BIOSAMPLE_CELLS)
    dropped = tx.excluded_cells(fr.CELLS, tx.SHARED_BIOSAMPLE_CELLS)
    if len(keep) != 11 or len(dropped) != 2:
        raise SystemExit(f"refused: the subset is {len(keep)} cells, not eleven (dropped {dropped})")

    eleven = [
        {
            "symbol": r["symbol"],
            "chrom": r["chrom"],
            "start": r["start"],
            "P": [r["P"][i] for i in keep],
            "D": rg._centre([r["y"][i] for i in keep]),  # noqa: SLF001 - M2i's within-gene centring
        }
        for r in primary
    ]
    thirteen = mod.side(primary, "M2i's thirteen cells, recomputed here for comparison")
    sub = mod.side(eleven, "POST HOC: the eleven cells, cardiac muscle cell and hepatocyte dropped")

    committed = json.loads(M2I_RESULT.read_text())
    reproduces = abs(thirteen["delta_mse"] - committed["primary"]["delta_mse"]) < 1e-12

    return {
        "label": "POST HOC (2): the eleven-cell subset",
        "note": tx.ELEVEN_CELL_NOTE,
        "nothing_is_fitted": fr.NOTHING_IS_FITTED,
        "scale": fr.SCALE,
        "interval_method": (
            f"M2i's gene-level percentile bootstrap: {fr.BOOTSTRAP['resamples']} resamples, seed"
            f" {fr.BOOTSTRAP['seed']}, percentiles {list(fr.FORM_INTERVAL)}, the scale held fixed"
        ),
        "cells_kept": [fr.CELLS[i] for i in keep],
        "cells_dropped": dropped,
        "attrition": attrition,
        "fresh_tables": len(counts),
        "thirteen_cells": thirteen,
        "eleven_cells": sub,
        "difference_in_delta_mse": sub["delta_mse"] - thirteen["delta_mse"],
        "thirteen_reproduces_the_committed_result": reproduces,
        "committed_delta_mse": committed["primary"]["delta_mse"],
        "reading": (
            "POST HOC and deciding nothing. With the two cells whose biosamples an already-read"
            f" experiment also used dropped, Delta MSE is {sub['delta_mse']:+.5f} ["
            f"{fr.FORM_INTERVAL[1] - fr.FORM_INTERVAL[0]:.2f}% CI {sub['ci'][0]:+.5f},"
            f" {sub['ci'][1]:+.5f}] on {sub['genes']:,} genes and {sub['pairs']:,} pairs, against"
            f" {thirteen['delta_mse']:+.5f} on all thirteen. This isolates the one non-independence"
            " M2i disclosed and nothing else; it does not revise M2i's reading and it says nothing"
            " about independence from AlphaGenome's training data, which the record cannot settle."
        ),
    }


def main() -> int:
    t0 = time.time()
    mod = m2i_run()
    sources = json.loads(fr.SOURCES.read_text())
    inputs = [
        mf.input_entry(p, partition=None)
        for p in (
            tx.METADATA.as_posix(),
            tx.FRESH_SOURCES.as_posix(),
            tx.READ_SOURCES.as_posix(),
            M2I_RESULT.as_posix(),
            fr.P_SOURCE,
            fr.SCALE_SOURCE,
            GENE_TABLE.as_posix(),
        )
    ]
    inputs += [
        mf.input_entry(f["path"], partition=None) for c in sources["biosamples"].values() for f in c["files"]
    ]
    manifest = {
        "sources": [
            {"accession": "f6e71c6", "version": "M2i, the registered result this qualifies"},
            {"accession": "0119387", "version": "M2i, the post-hoc stage"},
            {"accession": "ac5580b", "version": "M2h, the fixed scale"},
            {
                "accession": "alphagenome_track_metadata_copy.csv",
                "version": "AlphaGenome's recorded output-track metadata, the copy on disk",
            },
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38, GENCODE V29 gene models",
        "coordinates": "n/a: gene-level quantifications keyed by gene_id",
        "parameters": {
            "scale": fr.SCALE,
            "nothing_is_fitted": fr.NOTHING_IS_FITTED,
            "interval": list(fr.FORM_INTERVAL),
            "resamples": fr.BOOTSTRAP["resamples"],
            "seed": fr.BOOTSTRAP["seed"],
            "cells_dropped": list(tx.SHARED_BIOSAMPLE_CELLS),
        },
        "exclusions": [
            "non-rna_seq tracks, when matching a measured RNA quantification",
            "cardiac muscle cell and hepatocyte, in the eleven-cell subset alone",
        ],
        "partitions": {
            "fresh_thirteen": "M2i's thirteen fresh experiments",
            "already_read_twenty_seven": "items 7/7d/7h's twenty-seven already-read quantifications",
            "eleven_cells": "M2i's primary population on eleven of the thirteen cells",
        },
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }
    payload = {
        "result": "trainexp_independence",
        "date": time.strftime("%Y-%m-%d"),
        "lane": tx.LANE,
        "post_hoc": tx.POST_HOC,
        "what": (
            "M2i's independence from AlphaGenome's training data, measured from the recorded track"
            " metadata, and the eleven-cell subset. Both post hoc; both decide nothing"
        ),
        "matching": matching(),
        "eleven_cell_subset": eleven_cells(mod),
        "alphagenome_requests": 0,
        "money": "none: metadata and public files already on disk, no key read",
    }
    payload["seconds"] = round(time.time() - t0, 1)
    print(payload["matching"]["sentence_for_every_m2i_better"])
    print(payload["eleven_cell_subset"]["reading"])
    print("saved", save_result("trainexp_independence", payload, manifest=manifest))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
