# SPDX-License-Identifier: AGPL-3.0-or-later
"""M1b census: how many compiled rules gain a factor as their source, on what evidence, per cell type.

    uv run python scripts/regulators_census.py --motifs-only   # the JASPAR scan, cached per chromosome
    uv run python scripts/regulators_census.py [--chroms chr21] [--no-save]

Reads the 24 compiled programs, the ENCODE TF ChIP-seq peak clusters
(data/knowledge/regulators/, see genomeos/attribution/regulators.py CLUSTER_SOURCE), the cached JASPAR
2026 CORE vertebrate profiles and the cached GO annotations. No expression value, perturbation or model
output is read and nothing is downloaded. Writes data/results/regulators_census.json and the git-ignored
assignment table data/cache/regulators/assignment.tsv.gz (one row per compiled rule and factor occupying
its element in the rule's own cell), which is what compile.py would read.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import regulators as rg  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.lang.parser import parse  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "regulators_census"
CHUNK = 2000
OWN_CODE = (
    "genomeos/attribution/regulators.py",
    "scripts/regulators_inventory.py",
    "scripts/regulators_census.py",
    "scripts/regulators_register.py",
    "scripts/regulators_fetch.py",
    "scripts/regulators_run.py",
    "tests/test_regulators.py",
)


def program_rules(chrom: str) -> list[rg.CompiledRule]:
    path = Path(COMPILED_DIR) / f"noncoding_{chrom}.bio"
    return rg.compiled_rules(parse(path.read_text(), path.stem), path.stem)


def motif_cache(chrom: str) -> Path:
    return rg.CACHE / f"motifs_{chrom}.json.gz"


def scan_chrom(chrom: str, rules: list[rg.CompiledRule], keep: list[str]) -> dict:
    """Per element: how many factor names any JASPAR profile matches, and which of `keep` match."""
    from genomeos.coords import Locus
    from genomeos.genome import motifs as mo
    from genomeos.genome.fetch import REFERENCE
    from genomeos.genome.genome import Genome

    els = sorted({(r.element, r.chrom, r.start, r.end) for r in rules})
    genome = Genome.from_fasta(REFERENCE / f"{chrom}.fa.gz")
    profiles = mo.load_motifs()
    index = {f: i for i, f in enumerate(keep)}
    out: dict[str, list] = {}
    t0 = time.time()
    for c in range(0, len(els), CHUNK):
        part = els[c : c + CHUNK]
        seqs = [str(genome.fetch(Locus(ch, s, e))) for _, ch, s, e in part]
        names: list[set[str]] = [set() for _ in part]
        for prof, hits in mo.scan(profiles, seqs).items():
            fs = rg.motif_factor_names(prof)
            for i in hits:
                names[i].update(fs)
        for (el, *_), ns in zip(part, names, strict=True):
            out[el] = [len(ns), sorted(index[f] for f in ns if f in index)]
        print(f"  {chrom} {c + len(part)}/{len(els)} {time.time() - t0:.0f}s", flush=True)
    return {
        "chrom": chrom,
        "threshold": mo.RELATIVE_SCORE,
        "profiles": len(profiles),
        "keep": keep,
        "elements": out,
    }


def load_motifs(chrom: str, rules: list[rg.CompiledRule], keep: list[str]) -> dict:
    p = motif_cache(chrom)
    if p.exists():
        d = json.loads(gzip.decompress(p.read_bytes()))
        if d["keep"] == keep:
            return d
    d = scan_chrom(chrom, rules, keep)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(gzip.compress(json.dumps(d).encode()))
    return d


def manifest(chroms: list[str]) -> dict:
    from genomeos.genome import motifs as mo

    paths = [Path(COMPILED_DIR) / f"noncoding_{c}.bio" for c in chroms]
    inputs = [mf.input_entry(p, partition=None) for p in paths] + [
        mf.input_entry(rg.CLUSTERS, partition=None),
        mf.input_entry(mo.JASPAR_PATH, partition=None),
        mf.input_entry(rg.GOA, partition=None),
        mf.files_entry("JASPAR matches per compiled element, cached", [motif_cache(c) for c in chroms]),
    ]
    return {
        "sources": [
            {
                "accession": "the compiled non-coding programs, scripts/compile_genome_programs.py",
                "version": "the files on this machine, pinned by sha256 in inputs",
            },
            {
                "accession": rg.CLUSTER_SOURCE["url"],
                "version": "Last-Modified " + rg.CLUSTER_SOURCE["last_modified"],
            },
            {"accession": mo.JASPAR_URL, "version": "JASPAR 2026 CORE vertebrates, non-redundant"},
            {
                "accession": "Gene Ontology annotations, goa_human.gaf",
                "version": "the cached file, by sha256",
            },
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38 for the compiled elements, the peak clusters and the reference sequence",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "motif_relative_score": mo.RELATIVE_SCORE,
            "cell_aliases": rg.CELL_ALIASES,
            "go_activator": sorted(rg.GO_ACTIVATOR),
            "go_repressor": sorted(rg.GO_REPRESSOR),
            "no_expression_read": "the census opens no expression, perturbation or model value",
        },
        "exclusions": [],
        "partitions": "n/a: every compiled rule is counted",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chroms", default="")
    ap.add_argument("--motifs-only", action="store_true")
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args(argv)
    chroms = args.chroms.split(",") if args.chroms else list(rg.CHROMS)
    t0 = time.time()
    cl = rg.load_clusters()
    keep = sorted(cl.factors)
    print("clusters", len(cl.factors), "factors", f"{time.time() - t0:.0f}s", flush=True)
    if args.motifs_only:
        for chrom in chroms:
            load_motifs(chrom, program_rules(chrom), keep)
        return 0
    go = rg.load_go_directions()
    overall, per_cell = rg.empty_tally(), {}
    per_chrom = {}
    elements_seen: set[str] = set()
    occupied_elements = 0
    rg.CACHE.mkdir(parents=True, exist_ok=True)
    table = rg.CACHE / "assignment.tsv.gz"
    example = None
    with gzip.open(table, "wt") as fh:
        fh.write(
            "program\telement\tgene\tcell\trule_action\tfactor\tkind\tmotif\toccupied_in_cells"
            "\tdirection\tdirection_basis\n"
        )
        for chrom in chroms:
            rules = program_rules(chrom)
            mot = load_motifs(chrom, rules, keep)["elements"]
            ct = rg.empty_tally()
            for r in rules:
                occ = rg.element_occupancy(cl, r.chrom, r.start, r.end)
                count, idx = mot[r.element]
                cands = rg.candidates(r, occ, frozenset(keep[i] for i in idx), go)
                for t in (overall, ct, per_cell.setdefault(r.cell, rg.empty_tally())):
                    rg.tally(t, r, cands, count)
                if r.element not in elements_seen:
                    elements_seen.add(r.element)
                    occupied_elements += bool(occ)
                for c in cands:
                    if c["kind"] not in rg.SOURCE_KINDS:
                        continue
                    fh.write(
                        f"{r.program}\t{r.element}\t{r.gene}\t{r.cell}\t{r.action}\t{c['factor']}\t"
                        f"{c['kind']}\t{int(c['motif'])}\t{c['occupied_in']}\t{c['direction']}\t"
                        f"{c['direction_basis']}\n"
                    )
                    if example is None and c["direction"] != "unknown":
                        example = rg.biolang_rule(c["factor"], r, c["direction"], c["direction_basis"])
            per_chrom[chrom] = {"rules": ct["rules"], "rules_gaining_a_source": ct["rules_gaining_a_source"]}
            print(chrom, per_chrom[chrom], f"{time.time() - t0:.0f}s", flush=True)
    chip_cells = {c: sorted(fs) for c, fs in cl.factors_in().items()}
    cells = {c: rg.finish_tally(t) for c, t in sorted(per_cell.items(), key=lambda kv: -kv[1]["rules"])}
    for c, v in cells.items():
        v["factors_with_chip_in_this_cell"] = len(chip_cells.get(c, []))
    with_chip = {c: v for c, v in cells.items() if v["factors_with_chip_in_this_cell"]}
    fclass = Counter(rg.factor_class(go, f) for f in keep)
    tot = rg.finish_tally(overall, top=25)
    out = {
        "lane": rg.LANE,
        "question": "how many compiled rules gain a factor as their source, on what evidence, per cell type",
        "evidence_kinds": {
            "occupancy_same_cell_with_motif": "measured occupancy in the rule's own cell, with the "
            "factor's own JASPAR profile matching inside the element; occupancy is not function",
            "occupancy_same_cell": "measured occupancy in the rule's own cell, no match of the factor's "
            "profile (or no profile); possibly indirect; occupancy is not function",
            "occupancy_other_cell": "measured occupancy only in other cell types: not a source in the "
            "rule's cell and never emitted there",
            "motif_only": "predicted: a JASPAR match and no occupancy anywhere; a motif is not occupancy",
            "unknown": "no evidence of either kind: the source stays UNKNOWN",
        },
        "totals": {
            "rules": tot["rules"],
            "rules_gaining_a_source": tot["rules_gaining_a_source"],
            "share_gaining_a_source": tot["rules_gaining_a_source"] / tot["rules"] if tot["rules"] else None,
            "by_kind": tot["by_kind"],
            "candidate_pairs_in_own_cell": tot["candidate_pairs_in_own_cell"],
            "candidate_pairs_with_motif": tot["candidate_pairs_with_motif"],
            "direction_of_candidate_pairs": tot["direction_of_candidate_pairs"],
            "rules_with_a_directed_candidate": tot["rules_with_a_directed_candidate"],
            "distinct_factors_used": tot["distinct_factors"],
            "top_factors_by_rules": tot["top_factors_by_rules"],
            "sources_per_rule_median": tot["sources_per_rule_median"],
            "sources_per_rule_max": tot["sources_per_rule_max"],
            "motif_matched_factors_per_rule_median": tot["motif_matched_factors_per_rule_median"],
            "rules_with_any_motif_match": tot["rules_with_any_motif_match"],
            "elements": len(elements_seen),
            "elements_with_any_occupancy": occupied_elements,
        },
        "contexts": {
            "compiled_contexts": len(cells),
            "contexts_with_chip": len(with_chip),
            "rules_in_contexts_with_chip": sum(v["rules"] for v in with_chip.values()),
            "rules_in_contexts_without_chip": sum(v["rules"] for c, v in cells.items() if c not in with_chip),
            "chip_cell_labels_not_compiled": sorted(set(chip_cells) - set(cells)),
        },
        "per_cell_with_chip": with_chip,
        "per_cell_without_chip": {
            c: {k: v[k] for k in ("rules", "by_kind", "motif_matched_factors_per_rule_median")}
            for c, v in cells.items()
            if c not in with_chip
        },
        "per_chromosome": per_chrom,
        "factors": {
            "with_any_cluster": len(keep),
            "curated_direction_class": dict(sorted(fclass.items())),
        },
        "compile_spec": rg.COMPILE_SPEC,
        "example_directed_rule": example,
        "assignment_table": {
            "path": str(table),
            "committed": False,
            "rows": tot["candidate_pairs_in_own_cell"],
            "columns": "program element gene cell rule_action factor kind motif occupied_in_cells "
            "direction direction_basis",
        },
        "reading": (
            "a rule gains a source when at least one factor has a ChIP-seq peak cluster over its element "
            "in the rule's own cell type: measured occupancy, which is not function. Only the cell types "
            "ENCODE profiled can gain one; every other rule keeps an UNKNOWN source whatever its motifs"
        ),
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
        "seconds": round(time.time() - t0, 1),
    }
    if not args.no_save:
        print("saved", save_result(RESULT, out, manifest=manifest(chroms)))
    print(json.dumps(out["totals"], indent=1)[:3000])
    print(json.dumps(out["contexts"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
