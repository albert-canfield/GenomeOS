# SPDX-License-Identifier: AGPL-3.0-or-later
"""The orientation-aware node caller against the CTCF-only one, on four measurements.

    uv run python scripts/oriented_domains.py [--workers 6]

Three node sets, same chromosomes, same code downstream (genome/domains.py infer_domains):

- ctcf_only   the current caller: CTCF-only ENCODE elements as boundaries
- oriented    boundaries where consecutive single-strand CTCF motifs (JASPAR MA0139 inside elements
              with CTCF ChIP support) flip from reverse to forward
- oriented_ctcf_only  the same rule on CTCF-only elements alone, so the site set's effect shows

And the stricter site call registered on 2026-09-21 (docs/NODES-READER-WRITER.md), with the two
halves of it scored separately so that any movement can be attributed to one change:

- oriented_best_hit  the strand of the best-scoring hit, at the same 0.85 threshold
- oriented_strong    the set of hit strands, at the 0.95 threshold motifs.py calibrated
- oriented_strict    both: the best hit's strand at 0.95, the call under test

One caller is not registered and cannot decide anything: oriented_strict_0.90 was added after the
registered run showed that 0.95 keeps a small minority of the sites. It is the diagnostic that tells
"strictness does not help" apart from "0.95 kept nothing", and it is reported as an observation.

Measured on each:

1. Hi-C support: interior edges within 20 kb of a 4DN boundary call (five biosources), against as
   many random positions, and the share of measured boundaries an edge reaches.
2. Node content: of the elements the deletion archive scored (data/knowledge/alphagenome/all_elements),
   the share whose most-moved coding gene has its TSS in the element's node; against as many boundaries
   placed at random, and with both callers cut to the same boundary count (fewer boundaries keep more
   pairs together by construction).
3. Mouse synteny: mouse chr19 and chr11 nodes called the same way, held against human nodes through
   MGI orthology (genome/mouse.py compare_nodes).
4. HOXD: is there an edge in the published HOXD11 to HOXD13 boundary interval, and how many nodes do the
   nine HOXD genes fall into.

No model is called; everything is read from local caches and committed results.
"""

from __future__ import annotations

import bisect
import json
import random
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos.genome import (  # noqa: E402
    Annotation,
    IndexedGenome,
    default_gencode,
    mouse,
)
from genomeos.genome import hic as hic_mod  # noqa: E402
from genomeos.genome.domains import (  # noqa: E402
    STRICT_RELATIVE,
    ctcf_motif_sites,
    infer_domains,
    strand_variant,
)
from genomeos.genome.regulatory import load_ccres  # noqa: E402
from genomeos.predict.contact_maps import TOLERANCE, random_control  # noqa: E402
from genomeos.results import save_result  # noqa: E402

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
LOOSE_RELATIVE = 0.85  # the scan threshold the 2026-09-14 caller reads
CALLERS = (
    "ctcf_only",
    "oriented",
    "oriented_ctcf_only",
    "oriented_best_hit",
    "oriented_strong",
    "oriented_strict",
    "oriented_strict_0.90",  # a diagnostic, not a candidate: see the note in main()
)
# the pair whose boundary count sets the matched-resolution draw (unchanged from 2026-09-14)
MATCHED_PAIR = ("ctcf_only", "oriented")
BIOSOURCES = ("GM12878", "H1-hESC", "K562", "HepG2", "IMR-90")
ARCHIVE = Path("data/knowledge/alphagenome/all_elements")
SHUFFLES = 20
HOXD_BOUNDARY = (176_096_240, 176_109_754)  # Rodriguez-Carballo et al. 2017, as genomeos/benchmark/loci.py
HOXD_GENES = ("HOXD1", "HOXD3", "HOXD4", "HOXD8", "HOXD9", "HOXD10", "HOXD11", "HOXD12", "HOXD13")


def _tss(g) -> int:
    return g.locus.end - 1 if g.locus.strand.value == "-" else g.locus.start


#: caller name -> (score an element's hit must reach, take the best hit's strand)
SITE_CALLS = {
    "oriented": (LOOSE_RELATIVE, False),
    "oriented_ctcf_only": (LOOSE_RELATIVE, False),
    "oriented_best_hit": (LOOSE_RELATIVE, True),
    "oriented_strong": (STRICT_RELATIVE, False),
    "oriented_strict": (STRICT_RELATIVE, True),
    "oriented_strict_0.90": (0.90, True),
}


def site_call_counts(sites: dict) -> dict[str, int]:
    """The registered descriptive numbers: how many elements carry a site, and how many the strand
    rule can orient, under each call. They tell "improved nothing" from "kept nothing"."""
    out = {"elements": len(sites)}
    for name, (rel, best) in SITE_CALLS.items():
        if name == "oriented_ctcf_only":
            continue
        v = strand_variant(sites, relative=rel, best_hit=best)
        out[f"{name}_with_site"] = sum(1 for s in v.values() if s)
        out[f"{name}_orientable"] = sum(1 for s in v.values() if len(s) == 1)
    return out


def node_sets(chrom: str, length: int, ccres, ann, fetch, tag: str = "") -> tuple[dict[str, list], dict]:
    sites = ctcf_motif_sites(chrom, ccres, fetch, tag=tag)
    only = {c.id for c in ccres if c.cls == "CTCF-only"}
    sets = {"ctcf_only": infer_domains(chrom, length, ccres, ann)}
    for name, (rel, best) in SITE_CALLS.items():
        strands = strand_variant(sites, relative=rel, best_hit=best)
        if name == "oriented_ctcf_only":
            strands = {k: v for k, v in strands.items() if k in only}
        sets[name] = infer_domains(chrom, length, ccres, ann, orientation=strands)
    return sets, site_call_counts(sites)


def _edges(doms) -> list[int]:
    return [d.start for d in doms[1:]]


def _hits(positions: list[int], measured: list[int]) -> int:
    ms = sorted(measured)
    n = 0
    for p in positions:
        i = bisect.bisect_left(ms, p)
        if any(abs(ms[j] - p) <= TOLERANCE for j in (i - 1, i) if 0 <= j < len(ms)):
            n += 1
    return n


def _inside_share(starts: list[int], elements: list[tuple[int, int]]) -> tuple[int, int]:
    """(inside, total) for (element midpoint, gene TSS) pairs over nodes starting at `starts`."""
    inside = sum(1 for m, t in elements if bisect.bisect_right(starts, m) == bisect.bisect_right(starts, t))
    return inside, len(elements)


def human_chrom(chrom: str) -> dict:
    ref = Path(f"data/reference/{chrom}.fa.gz")
    gff = default_gencode({chrom})
    ccres = load_ccres(chrom)
    if not ref.exists() or not gff or not ccres:
        return {"chrom": chrom, "skipped": True}
    from genomeos.coords import Locus

    g = IndexedGenome(ref)
    length = g.lengths[chrom]
    ann = Annotation.from_gff3(gff, {chrom})
    try:
        sets, site_counts = node_sets(
            chrom, length, ccres, ann, lambda s, e: str(g.fetch(Locus(chrom, s, min(e, length))))
        )
    finally:
        g.close()
    coding = {
        x.symbol: x for x in ann.genes.values() if x.locus.chrom == chrom and x.type == "protein_coding"
    }
    anyg = {x.symbol: x for x in ann.genes.values() if x.locus.chrom == chrom}
    out: dict = {"chrom": chrom, "length": length, "callers": {}, "site_calls": site_counts}
    measured = {b: hic_mod.load_boundaries(b, chrom) for b in BIOSOURCES}
    # the archive's scored elements: (midpoint, TSS of the gene the deletion moves most)
    pairs_coding: list[tuple[int, int]] = []
    pairs_any: list[tuple[int, int]] = []
    arch = ARCHIVE / f"{chrom}.json"
    archive_rows = 0
    if arch.exists():
        for e in json.loads(arch.read_text()):
            archive_rows += 1
            mid = (e["start"] + e["end"]) // 2
            pc = (e.get("predicted_coding") or {}).get("gene")
            if pc in coding:
                pairs_coding.append((mid, _tss(coding[pc])))
            pa = (e.get("predicted") or {}).get("gene")
            if pa in anyg:
                pairs_any.append((mid, _tss(anyg[pa])))
    rng = random.Random(7)
    for name, doms in sets.items():
        edges = _edges(doms)
        lengths = [d.length for d in doms]
        row: dict = {
            "nodes": len(doms),
            "edges": len(edges),
            "median_length": sorted(lengths)[len(lengths) // 2] if lengths else None,
            "hic": {},
        }
        for b, m in measured.items():
            if not m:
                continue
            row["hic"][b] = {
                "edges_supported": _hits(edges, m),
                "random_control": random_control(edges, m, length) if edges else None,
                "measured": len(m),
                "measured_reached": _hits(m, edges),
            }
        starts = [d.start for d in doms]
        if pairs_coding:
            ins, tot = _inside_share(starts, pairs_coding)
            ins_any, tot_any = _inside_share(starts, pairs_any)
            # the control the Hi-C comparison uses: as many boundaries placed at random. Fewer
            # boundaries keep more pairs together by construction; the control absorbs that.
            null = []
            for _ in range(SHUFFLES):
                st = [0, *sorted(rng.randint(1, length - 1) for _ in edges)]
                null.append(_inside_share(st, pairs_coding)[0])
            row["content"] = {
                "coding_named": tot,
                "coding_inside": ins,
                "coding_inside_random_mean": round(sum(null) / len(null), 1),
                "any_named": tot_any,
                "any_inside": ins_any,
            }
        # genes per node for the mouse comparison's human index
        row["index"] = {s: doms_i for doms_i, d in enumerate(doms) for s in d.genes}
        row["hoxd"] = None
        if chrom == "chr2":
            lo, hi = HOXD_BOUNDARY
            near = min(edges, key=lambda x: min(abs(x - lo), abs(x - hi))) if edges else None
            nodes_of = {s: bisect.bisect_right(starts, _tss(coding[s])) for s in HOXD_GENES if s in coding}
            groups: dict[int, list[str]] = {}
            for s, n in nodes_of.items():
                groups.setdefault(n, []).append(s)
            row["hoxd"] = {
                "edge_in_published_interval": any(lo <= x <= hi for x in edges),
                "nearest_edge": near,
                "nearest_edge_distance": (
                    0 if near is not None and lo <= near <= hi else min(abs(near - lo), abs(near - hi))
                )
                if near is not None
                else None,
                "nodes": len(groups),
                "groups": [sorted(v, key=lambda s: int(s[4:])) for _, v in sorted(groups.items())],
            }
        out["callers"][name] = row
    # placement against resolution: the callers cut to the same number of boundaries
    if pairs_coding:
        k = min(len(_edges(sets[n])) for n in MATCHED_PAIR)
        for name in CALLERS:
            edges = _edges(sets[name])
            if len(edges) < k or not k:
                out["callers"][name]["content"]["matched"] = None
                continue
            shares = []
            for _ in range(SHUFFLES):
                st = [0, *sorted(rng.sample(edges, k))]
                shares.append(_inside_share(st, pairs_coding)[0])
            out["callers"][name]["content"]["matched"] = {
                "boundaries": k,
                "coding_inside_mean": round(sum(shares) / len(shares), 1),
            }
    out["archive_rows"] = archive_rows
    return out


def mouse_chrom(chrom: str) -> dict:
    from genomeos.coords import Locus

    fa = mouse.fetch_sequence(chrom)
    gff = mouse.fetch_gencode(chrom)
    ccres = mouse.load_ccres(chrom)
    ann = Annotation.from_gff3(gff, {chrom})
    g = IndexedGenome(fa)
    length = g.lengths[chrom]
    try:
        sets, _ = node_sets(
            chrom, length, ccres, ann, lambda s, e: str(g.fetch(Locus(chrom, s, min(e, length)))), tag="mm10_"
        )
    finally:
        g.close()
    by_symbol = {x.symbol: x for x in ann.genes.values() if x.locus.chrom == chrom}
    out = {"chrom": chrom, "callers": {}}
    for name, doms in sets.items():
        rows = []
        for d in doms:
            dd = {"id": d.id, "start": d.start, "end": d.end}
            dd["coding_symbols"] = [
                s for s in d.genes if s in by_symbol and by_symbol[s].type == "protein_coding"
            ]
            rows.append(dd)
        out["callers"][name] = {
            "nodes": len(doms),
            "rows": rows,
            "median_length": sorted(d.length for d in doms)[len(doms) // 2],
        }
    return out


def main() -> None:
    args = sys.argv[1:]
    workers = int(args[args.index("--workers") + 1]) if "--workers" in args else 6
    with ProcessPoolExecutor(workers) as pool:
        human = [h for h in pool.map(human_chrom, CHROMS) if not h.get("skipped")]
        mice = list(pool.map(mouse_chrom, ["chr19", "chr11"]))
    print("chromosomes:", [h["chrom"] for h in human], flush=True)
    result: dict = {"callers": {}, "chromosomes": [h["chrom"] for h in human]}
    site_calls: dict[str, int] = {}
    for h in human:
        for k, v in h.get("site_calls", {}).items():
            site_calls[k] = site_calls.get(k, 0) + v
    result["site_calls"] = site_calls
    print("site calls", json.dumps(site_calls), flush=True)
    orthology = mouse.load_orthology()
    for name in CALLERS:
        agg: dict = {"nodes": 0, "edges": 0, "hic": {}, "content": {}, "mouse": {}}
        lengths_all = []
        for h in human:
            r = h["callers"][name]
            agg["nodes"] += r["nodes"]
            agg["edges"] += r["edges"]
            lengths_all.append(r["median_length"])
            for b, v in r["hic"].items():
                a = agg["hic"].setdefault(
                    b, {"edges": 0, "supported": 0, "rc_weighted": 0.0, "measured": 0, "reached": 0}
                )
                n_edges = r["edges"]
                a["edges"] += n_edges
                a["supported"] += v["edges_supported"]
                a["rc_weighted"] += (v["random_control"] or 0.0) * n_edges
                a["measured"] += v["measured"]
                a["reached"] += v["measured_reached"]
            if "content" in r:
                c = agg["content"]
                for k, v in r["content"].items():
                    if k == "matched":
                        if v:
                            c["matched_boundaries"] = c.get("matched_boundaries", 0) + v["boundaries"]
                            c["matched_inside"] = c.get("matched_inside", 0.0) + v["coding_inside_mean"]
                            c["matched_named"] = c.get("matched_named", 0) + r["content"]["coding_named"]
                        continue
                    c[k] = c.get(k, 0) + v
                c.setdefault("chromosomes", []).append(h["chrom"])
        agg["median_of_chromosome_median_lengths"] = sorted(lengths_all)[len(lengths_all) // 2]
        for a in agg["hic"].values():
            a["supported_share"] = round(a["supported"] / a["edges"], 4) if a["edges"] else None
            a["random_control"] = round(a.pop("rc_weighted") / a["edges"], 4) if a["edges"] else None
            a["enrichment"] = (
                round(a["supported_share"] / a["random_control"], 3) if a["random_control"] else None
            )
            a["measured_reached_share"] = round(a["reached"] / a["measured"], 4) if a["measured"] else None
        c = agg["content"]
        if c.get("coding_named"):
            c["coding_inside_share"] = round(c["coding_inside"] / c["coding_named"], 4)
            c["coding_inside_random_share"] = round(c["coding_inside_random_mean"] / c["coding_named"], 4)
            c["coding_excess_over_random"] = round(
                c["coding_inside_share"] - c["coding_inside_random_share"], 4
            )
            if c.get("matched_named"):
                c["coding_inside_share_at_matched_resolution"] = round(
                    c["matched_inside"] / c["matched_named"], 4
                )
            c["any_inside_share"] = round(c["any_inside"] / c["any_named"], 4)
        # the human index for this caller: symbol -> chrom:D(n)
        index: dict[str, str] = {}
        for h in human:
            for s, i in h["callers"][name]["index"].items():
                index.setdefault(s.upper(), f"{h['chrom']}:D{i + 1}")
        for m in mice:
            rows = m["callers"][name]["rows"]
            cmp = mouse.compare_nodes(rows, index, orthology)
            agg["mouse"][m["chrom"]] = {
                k: cmp[k]
                for k in (
                    "tested",
                    "conserved",
                    "split_adjacent",
                    "split_scattered",
                    "fraction_conserved",
                    "fraction_same_neighbourhood",
                )
            } | {
                "mouse_nodes": m["callers"][name]["nodes"],
                "mouse_median_length": m["callers"][name]["median_length"],
            }
        hoxd = next((h["callers"][name]["hoxd"] for h in human if h["chrom"] == "chr2"), None)
        agg["hoxd"] = hoxd
        result["callers"][name] = agg
        print(name, json.dumps({k: v for k, v in agg.items() if k != "hic"})[:1500], flush=True)
        print("  hic", json.dumps(agg["hic"]), flush=True)
    result["evidence"] = {
        "ctcf_only": "inferred: CTCF-only ENCODE elements as boundaries (the committed nodes)",
        "oriented": "inferred: reverse-to-forward CTCF motif flips (genome/domains.py oriented_boundaries)",
        "oriented_strict": (
            "inferred: the same flips, with each element's strand taken from its best MA0139 hit at a "
            "relative score of 0.95 (genome/domains.py strand_variant); registered 2026-09-21"
        ),
        "oriented_strict_0.90": (
            "inferred: the same call at 0.90; not registered, a diagnostic for how much of the strict "
            "call's movement is the sites it drops rather than the strand it picks"
        ),
        "hic": hic_mod.EVIDENCE,
        "content": "predicted: AlphaGenome deletion archive, most-moved coding gene; no new model call",
        "mouse": "curated: MGI mouse-human homology; mouse nodes on mm10 by the same caller",
    }
    result["controls"] = {
        "hic": "as many uniformly random positions (contact_maps.random_control)",
        "content": f"as many boundaries placed uniformly at random, {SHUFFLES} draws; and the callers cut to "
        f"the smaller boundary count of ctcf_only and oriented ({SHUFFLES} draws), placement against "
        "resolution",
    }
    save_result("domains_oriented_comparison", result)
    print("saved", flush=True)


if __name__ == "__main__":
    main()
