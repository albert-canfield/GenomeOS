#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Does a cached scorer answer carry enough to tie a gene row to a locus?

    uv run python scripts/gene_row_locus.py [--result gene_row_locus] [--chrom chr21 ...]

Reads only what is already on disk: the per-element response cache the deletion sweep wrote, the
compact tables its committed results point at, and GENCODE v50. No model request, no network, no
money. It reads the compiled rules not at all and changes nothing.

The three stages, each of which can be read on its own:

1. the key vocabulary of every cached answer, by streaming, which is the test for a ninth gene-row
   field -- a gene id, a coordinate, a strand, a transcript, anything locus-identifying;
2. every gene row's ``n_tracks``, which exposes the rows that pooled more than one response row,
   since one response row cannot exceed the response's own recorded column count;
3. every (element, gene name) pair whose name has two or more GENCODE v50 loci with a gene body
   inside the scorer window -- the kind the 547 are -- classified by whether anything on disk tells
   those loci apart.

A chromosome is held one at a time and dropped before the next, so the bill is the largest
chromosome and not the tree.
"""

from __future__ import annotations

import argparse
import collections
import gc
import json
import re
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import gene_row_locus as gr  # noqa: E402
from genomeos.attribution import targets  # noqa: E402
from genomeos.attribution.gene_identity import HALF_WINDOW, annotation_path  # noqa: E402
from genomeos.results import save_result  # noqa: E402

CHROMS = tuple(f"chr{i}" for i in list(range(1, 23))) + ("chrX", "chrY")

#: The caches this lane reads, all of them already on disk.
ARCHIVES = Path("data/knowledge/alphagenome/elements")
COMPACT = Path("data/knowledge/alphagenome/all_elements")
HCT116 = Path("data/knowledge/alphagenome/elements_hct116")

OWN_CODE = ("genomeos/attribution/gene_row_locus.py", "scripts/gene_row_locus.py")

#: A key is a schema key rather than an element id when it occurs far more often than once. The
#: archives are dicts keyed by element id, so every id appears exactly once and no schema key does.
SCHEMA_KEY_FLOOR = 100

#: The shapes the cached element ids take, so a rare key that is NOT an id can be shown separately:
#: that is where a field present on only a few answers would hide.
ID_SHAPES = re.compile(r"^(EH38[A-Z0-9]+|hs\d+|satmut_\S+|chr[0-9XY]+[_:]\S*|[A-Z][A-Z0-9_.\-]*)$")


def schema_keys(vocab: dict[str, int], floor: int = SCHEMA_KEY_FLOOR) -> dict[str, int]:
    return {k: v for k, v in sorted(vocab.items(), key=lambda kv: -kv[1]) if v >= floor}


def not_id_shaped(vocab: dict[str, int], floor: int = SCHEMA_KEY_FLOOR) -> dict[str, int]:
    """The low-frequency keys that are NOT element ids, which is where a rare field would hide."""
    return {k: v for k, v in vocab.items() if v < floor and not ID_SHAPES.match(k)}


def main() -> int:
    started = time.time()
    ap = argparse.ArgumentParser()
    ap.add_argument("--result", default="gene_row_locus")
    ap.add_argument("--chrom", nargs="*", default=list(CHROMS))
    args = ap.parse_args()
    chroms = [c for c in args.chrom if (ARCHIVES / f"{c}.json.gz").exists()]

    arch_paths = [ARCHIVES / f"{c}.json.gz" for c in chroms]
    compact_paths = [COMPACT / f"{c}.json" for c in chroms if (COMPACT / f"{c}.json").exists()]
    loose = sorted(ARCHIVES.glob("*/*.json"))
    hct = sorted(HCT116.glob("*/*.json"))

    # --- 1. the key vocabulary: is there a ninth gene-row field anywhere? ---------------------
    vocab = gr.row_key_vocabulary(arch_paths)
    compact_vocab = gr.row_key_vocabulary(compact_paths)
    loose_vocab = gr.row_key_vocabulary(loose)
    hct_vocab = gr.row_key_vocabulary(hct)

    # --- 2. every gene row's n_tracks ---------------------------------------------------------
    scan = gr.row_scan(arch_paths)

    # --- 3. the same-named-locus ambiguities --------------------------------------------------
    per_chrom: dict[str, dict] = {}
    pairs: list[dict] = []
    examples: list[dict] = []
    elements = 0
    for c in chroms:
        repeated = gr.repeated_names(c)
        coding = gr.coding_names(c)
        r = targets.ElementResponses(root=ARCHIVES)
        r._load(c)  # noqa: SLF001  # one chromosome held, then dropped
        got: list[dict] = []
        for rec in r._archive.values():  # noqa: SLF001
            elements += 1
            got.extend(gr.classify_element(rec, repeated, coding))
        per_chrom[c] = gr.tally(got)
        for p in got:
            if p["verdict"] == "pinned" and len(examples) < 8:
                examples.append(p)
        pairs.extend(got)
        del r
        gc.collect()
        print(f"{c}: {per_chrom[c]['ambiguous_element_name_pairs']} ambiguous pairs", flush=True)

    overall = gr.tally(pairs)
    coding_only = overall["restricted_to_names_with_a_protein_coding_locus_on_the_chromosome"]

    merged_by_type: collections.Counter[str] = collections.Counter()
    for p in pairs:
        if p["verdict"] == "merged":
            merged_by_type["+".join(sorted(set(p["gene_types_in_window"])))] += 1

    payload = {
        "result": args.result,
        "lane": "lane-generow",
        "question": (
            "does a cached scorer answer already record enough to tie a gene row to a locus, which is "
            "what would let lane-identity's gene-identity comparison fire on the 547 at all"
        ),
        "status": (
            "descriptive, and a reading question only: nothing was scored, re-scored, re-resolved or "
            "changed. It reads the caches already on disk and the committed chain's own code"
        ),
        "answer": (
            "the cached answers carry NO gene id, coordinate, strand, transcript or gene type per gene "
            "row -- eight keys and no ninth, over every cached answer. But the response did carry the "
            "distinction and the chain destroyed it: a cached row whose n_tracks exceeds the response's "
            "own 371 columns pooled two or more response rows under one gene_name, and there are "
            f"{scan['merged_rows']} such rows. One field does settle the same-named-locus question "
            "where it applies -- the row's own `gene`, when AlphaGenome returned a bare Ensembl id "
            "because its annotation gives that locus no symbol -- and it settles "
            f"{overall['settled_by_what_is_on_disk']} of the "
            f"{overall['ambiguous_element_name_pairs']} ambiguous (element, name) pairs. It settles "
            f"{coding_only['pinned']} of the pairs whose name has a protein_coding locus on the "
            "chromosome, which is the class every compiled `predicted_coding` target belongs to"
        ),
        "chromosomes": chroms,
        "cached_elements_read": elements,
        "gene_row_contents": {
            "keys": list(gr.GENE_ROW_KEYS),
            "quoted_from_the_data": {
                "note": (
                    "a row of data/knowledge/alphagenome/elements/chr22/"
                    "satmut_chr22_27842324_27842349.json, as it is on disk"
                ),
                "row": {
                    "gene": "MN1",
                    "n_tracks": 371,
                    "mean_log2fc": 0.0049,
                    "max_drop_log2fc": -0.0237,
                    "max_drop_tissue": "heart",
                    "max_rise_log2fc": 0.043,
                    "max_rise_tissue": "HepG2",
                    "by_cell": {"HepG2": 0.0397, "IMR-90": 0.033, "K562": -0.0055, "GM12878": 0.0024},
                },
            },
            "nothing_locus_identifying": (
                "no gene_id, no coordinate, no strand, no transcript, no gene_type. `gene` is "
                "AlphaGenome's gene_name column and nothing else from its gene axis is kept"
            ),
        },
        "key_vocabulary": {
            "per_element_archives": {
                "schema_keys": schema_keys(vocab),
                "low_frequency_keys_that_are_not_element_ids": not_id_shaped(vocab),
                "distinct_keys": len(vocab),
            },
            "compact_tables_the_results_point_at": {
                "all_keys": schema_keys(compact_vocab, 1),
                "distinct_keys": len(compact_vocab),
            },
            "loose_per_element_files": {"schema_keys": schema_keys(loose_vocab, 1), "files": len(loose)},
            "the_newer_hct116_answers": {"schema_keys": schema_keys(hct_vocab), "files": len(hct)},
        },
        "the_merge": {
            "response_columns": gr.TRACKS,
            "how_that_is_known": (
                "the chain records the response's own track table per answer: model.tracks = 371 under "
                "a single tracks_sha256 over all 705 answers that carry a run record "
                "(data/knowledge/alphagenome/elements_hct116, alphagenome client 0.9.0, "
                "variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ'], 2026-09-28)"
            ),
            **{k: v for k, v in scan.items() if k != "merged_names"},
            "merged_names_count": len(scan["merged_names"]),
            "merged_names_top": dict(list(scan["merged_names"].items())[:20]),
            "merged_by_gene_type_of_the_window_loci": dict(merged_by_type.most_common()),
        },
        "dropped_at": [dict(d) for d in gr.DROPPED_AT],
        "same_named_locus_ambiguities": overall,
        "per_chromosome": per_chrom,
        "examples_of_a_pair_the_cache_settles": examples,
        "the_symbol_versus_locus_test": {
            "what_must_be_shown": (
                "that whatever is proposed distinguishes two loci of the SAME NAME, not merely that it "
                "exists. lane-identity's 295 of 547 with exactly one protein_coding locus among the "
                "window candidates failed this: `predicted_coding` restricts by symbols that are "
                "protein_coding SOMEWHERE on the chromosome, a property of the symbol, so it resolved "
                "none of the 547"
            ),
            "how_the_pinning_field_passes_it": (
                "the two loci carry one GENCODE symbol and the cache gives them DIFFERENT names -- one "
                "the symbol, the other its bare Ensembl id. No property of a symbol can differ between "
                "two loci that have the same symbol, so the distinction is locus-level by construction"
            ),
            "and_why_it_still_settles_none_of_the_547": (
                "a compiled target comes from a `predicted_coding` block, so its name has a "
                f"protein_coding locus on the chromosome. Of the {coding_only['ambiguous']} ambiguous "
                f"pairs whose name has one, {coding_only['pinned']} are pinned. AlphaGenome emits a "
                "bare id only where its annotation gives the locus no symbol, which does not happen "
                "for the loci of a protein-coding name"
            ),
            "the_condition_on_pinning": (
                "pinning names the remaining locus by elimination and assumes AlphaGenome's gene set "
                "holds the same loci of that name in that window as GENCODE v50 does. MATR3, NOX5 and "
                "ZNF724 have one GENCODE v50 locus each and a cached row that merges two response "
                "rows, so that assumption fails for them"
            ),
        },
        "what_a_future_scoring_run_would_have_to_record": [
            "the Ensembl gene id of each response row, read from the gene axis beside gene_name at "
            "genomeos/predict/alphagenome_adapter.py:227, and carried through aggregate as the "
            "accumulator key in place of the name (genomeos/predict/enhancer_target.py:198), so that "
            "two loci of one name stay two rows",
            "with it, the annotation release the service's gene axis comes from, because an id the "
            "service returns is only comparable to a GENCODE v50 id when the two sets hold the same "
            "loci, and MATR3, NOX5 and ZNF724 show they do not",
            "the merge cannot be undone after the fact: n_tracks says that a row pooled two loci, and "
            "nothing in the pooled numbers says which contributed what",
        ],
        "what_this_cannot_establish": gr.what_this_is_not(),
        "nothing_was_changed": (
            "no compiled rule, no cached answer, no result and no part of the scoring chain was "
            "altered. data/results/gene_identity.json is untouched"
        ),
        "alphagenome_requests": 0,
        "network_requests": 0,
        "money": "none: every input was already on disk",
        "seconds": round(time.time() - started, 1),
    }

    inputs = [mf.input_entry(p) for p in arch_paths]
    inputs += [mf.input_entry(p) for p in compact_paths]
    inputs += [mf.input_entry(annotation_path(c)) for c in chroms if annotation_path(c).exists()]
    inputs += [mf.input_entry(p) for p in loose]
    inputs += [mf.input_entry(p) for p in hct]
    inputs.sort(key=lambda e: e["path"])

    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "the deletion sweep's per-element response cache, as the project wrote it",
                "version": (
                    "data/knowledge/alphagenome/elements/chr*.json.gz plus the loose per-element files "
                    "under data/knowledge/alphagenome/elements/<chrom>/, every gene the scorer read in "
                    "each element's 1 Mb window"
                ),
            },
            {
                "accession": "the compact tables the enhancer_targets_all_chr* results point at",
                "version": (
                    "data/knowledge/alphagenome/all_elements/chr*.json, declared here as bytes because "
                    "they are reached through an elements_where POINTER field inside a result"
                ),
            },
            {
                "accession": "the HCT116 deletion answers, the only ones carrying a run record",
                "version": (
                    "data/knowledge/alphagenome/elements_hct116/<chrom>/*.json, read for model.tracks "
                    "and tracks_sha256, which fix the response's own column count at 371"
                ),
            },
            {"accession": "GENCODE", "version": "v50, data/reference/gencode_v50_chr*.gff3.gz"},
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "half_window": HALF_WINDOW,
            "half_window_imported_from": (
                "genomeos.attribution.not_open_profile.SCORER_HALF_WINDOW, via gene_identity, so both "
                "lanes read the same window"
            ),
            "candidate_rule": (
                "a GENCODE v50 gene is a candidate when its gene BODY overlaps the element midpoint "
                "plus or minus the half-window: lane-identity's own rule"
            ),
            "response_columns": gr.TRACKS,
            "schema_key_floor": SCHEMA_KEY_FLOOR,
            "gene_row_keys": list(gr.GENE_ROW_KEYS),
            "verdicts": dict(gr.CLASSES),
            "chromosomes": chroms,
        },
        "exclusions": [
            "no cached element is excluded: every element in every archive of every named chromosome "
            "is classified, and the verdicts sum to the ambiguous pairs with nothing in a residue",
            "no new cut-off is introduced: the half-window and the candidate rule are imported from "
            "where the project already fixed them, and 371 is read from the chain's own record of the "
            "response rather than chosen",
            "the compiled rules are not read at all, so no figure here is at the rule denominator and "
            "the 547 are neither recounted nor revised",
        ],
        "partitions": {
            "by_chromosome": "the archive the element was cached in",
            "by_verdict": "merged, pinned, partly_pinned, nothing, as the module defines them",
            "by_name_coding": (
                "whether the name has a protein_coding locus anywhere on the chromosome, which is the "
                "property `predicted_coding` restricts by and is a property of the symbol, not a locus"
            ),
            "by_gene_type_of_the_window_loci": "the GENCODE gene_types of the loci sharing the name",
        },
        "code_cleanliness": mf.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
    }

    path = save_result(args.result, payload)
    print(f"wrote {path}")
    print(json.dumps(overall["by_verdict"], indent=1))
    print(f"ambiguous pairs: {overall['ambiguous_element_name_pairs']}")
    print(f"settled from disk: {overall['settled_by_what_is_on_disk']}")
    print(f"of names with a protein_coding locus: {coding_only}")
    print(f"merged gene rows: {scan['merged_rows']} of {scan['gene_rows']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
