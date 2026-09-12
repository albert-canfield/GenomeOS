# SPDX-License-Identifier: AGPL-3.0-or-later
"""Constraint on two axes: across 241 mammals and within 76,156 people, block by block.

The budget (attribution/budget.py) reads one axis, Zoonomia phyloP across placental
mammals: a block that passes is under selection over tens of millions of years. The
second axis is variation among people. gnomAD's Gnocchi score (Chen et al. 2024) is a
Z score per kilobase of observed against expected variants in 76,156 genomes; Z >= 2.18
is the top decile of constrained non-coding sequence, Z >= 4.0 the top percentile. Read
together, each block and each attributed element falls in one of four cases:

    mammals      humans       case       reading
    constrained  constrained  syntax     held still over deep time and among people
    free         free         tolerant   free to change on both scales, or a value slot
    free         constrained  recent     primate or human function, or recent selection
    constrained  free         relaxed    held across mammals, variable among people: a
                                         frame whose *value* varies, or function lost

Both tracks are read per base over HTTP range requests with the in-house bigWig reader
and never stored; a chromosome's blocks cost a few hundred kilobytes because Gnocchi is
one value per kilobase. Nothing here reads a case as a verdict: within-human constraint
at one kilobase says little about one base, and a "tolerant" block may hold a value slot
(eye colour's rs12913832 sits in a kilobase Gnocchi does not score at all). The cases are
evidence with a confidence, reported next to the tier, and the controls (coding exons
must read as constrained, VISTA positives more than negatives) say how far to trust the
axis on each chromosome. Design and numbers: docs/GRAMMAR-BY-COMPARISON.md.
"""

from __future__ import annotations

import glob
import time
from pathlib import Path
from typing import Any

from genomeos.attribution.bigwig import BigWig, IntervalStats
from genomeos.attribution.budget import CONSTRAINED_MIN
from genomeos.attribution.constraint import PHYLOP_241_URL, PHYLOP_THRESHOLD, phylop_over_blocks
from genomeos.results import RESULTS_DIR, load_result, save_result

GNOCCHI_URL = "https://hgdownload.soe.ucsc.edu/gbdb/hg38/gnomAD/mutConstraint/mutConstraint.bw"
GNOCCHI_THRESHOLD = 2.18  # top decile of constrained non-coding kilobases (Chen et al. 2024)
GNOCCHI_STRONG = 4.0  # top percentile
HUMAN_MIN_FRACTION = 0.25  # a quarter of the measured kilobases in the top decile: constrained in humans
ELEMENT_MAMMAL_MIN = 0.20  # the bar the constrained-targets run set for an element
CASES = {
    (True, True): "syntax",
    (False, False): "tolerant",
    (False, True): "recent",
    (True, False): "relaxed",
}
CASE_ORDER = ("syntax", "relaxed", "recent", "tolerant")
CASE_LABELS = {
    "syntax": "constrained across mammals and among people",
    "relaxed": "constrained across mammals, variable among people: a frame whose value varies, or lost",
    "recent": "free across mammals, constrained among people: recent function or selection",
    "tolerant": "free on both scales: tolerant, or a value slot",
}
KNOWN_VALUES = [
    {
        "rsid": "rs12913832",
        "chrom": "chr15",
        "pos": 28120472,
        "note": "HERC2 intron 86 enhancer of OCA2; A brown, G blue (Sturm 2008)",
    }
]
EVIDENCE = {
    "mammals": f"curated: Zoonomia phyloP over 241 placental mammals, constrained at >= {PHYLOP_THRESHOLD}",
    "humans": (
        f"curated: gnomAD Gnocchi, Z per kb from 76,156 genomes, constrained at >= {GNOCCHI_THRESHOLD}"
        f" (top decile), strong at >= {GNOCCHI_STRONG}"
    ),
    "case": "inferred: the two axes read together; a best guess with a confidence, never a verdict",
}


def merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Sorted, non-overlapping union (the reader requires intervals that do not overlap)."""
    out: list[tuple[int, int]] = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def stats_dict(s: IntervalStats | None) -> dict | None:
    """The kilobase summary of one interval, or None where the track has no value."""
    if s is None or not s.bases:
        return None
    return {
        "bases": s.bases,
        "mean": round(s.mean, 3) if s.mean is not None else None,
        "maximum": round(s.maximum, 3),
        "fraction_above": round(s.above / s.bases, 4),
        "strong": bool(s.maximum >= GNOCCHI_STRONG),
    }


def case_of(
    mammal_fraction: float | None,
    human_fraction: float | None,
    mammal_min: float = CONSTRAINED_MIN,
    human_min: float = HUMAN_MIN_FRACTION,
) -> dict | None:
    """One of four cases from the two constrained fractions, with a confidence.

    Confidence rises with the distance of each fraction from its bar and is capped at 0.6:
    two summary statistics over a block are a reading, not a measurement of any base.
    """
    if mammal_fraction is None or human_fraction is None:
        return None
    m = mammal_fraction >= mammal_min
    h = human_fraction >= human_min
    dm = min(1.0, abs(mammal_fraction - mammal_min) / mammal_min)
    dh = min(1.0, abs(human_fraction - human_min) / human_min)
    return {
        "case": CASES[(m, h)],
        "mammals": "constrained" if m else "free",
        "humans": "constrained" if h else "free",
        "confidence": round(0.3 + 0.3 * min(dm, dh), 2),
    }


def gnocchi_over(
    chrom: str, intervals: list[tuple[int, int]], progress=None, url: str = GNOCCHI_URL
) -> tuple[list[IntervalStats], dict]:
    """Gnocchi summaries for non-overlapping intervals of one chromosome, and what it cost."""
    bw = BigWig(url)
    try:
        t0 = time.time()
        stats = bw.summarise(chrom, intervals, GNOCCHI_THRESHOLD, progress=progress)
        cost = {
            "requests": bw.src.requests,
            "mb_fetched": round(bw.src.bytes_fetched / 1e6, 2),
            "seconds": round(time.time() - t0, 1),
        }
    finally:
        bw.close()
    return stats, cost


def classify_blocks(blocks: list[dict], stats: list[IntervalStats | None]) -> list[dict]:
    """One row per budget block: the tier it has, the human axis, and the case."""
    rows = []
    for b, s in zip(blocks, stats, strict=True):
        ph = b.get("phylop") or {}
        mf = ph.get("fraction_above")
        g = stats_dict(s)
        rows.append(
            {
                "start": b["start"],
                "end": b["end"],
                "length": b["length"],
                "class": b["class"],
                "tier": b["guess"]["tier"],
                "mammal_fraction": mf,
                "gnocchi": g,
                "case": case_of(mf, g["fraction_above"] if g else None),
            }
        )
    return rows


def elements_of(chrom: str, results_dir: Path = RESULTS_DIR) -> list[dict]:
    """Every element the target runs scored on this chromosome, constrained run first, by id."""
    seen: set[str] = set()
    out: list[dict] = []
    for name, origin in (("constrained_targets", "constrained"), ("enhancer_targets", "uniform")):
        r = load_result(f"{name}_{chrom}", results_dir) or {}
        for e in r.get("elements", []):
            if e["id"] in seen:
                continue
            seen.add(e["id"])
            pc = e.get("predicted_coding") or {}
            out.append(
                {
                    "id": e["id"],
                    "start": e["start"],
                    "end": e["end"],
                    "origin": origin,
                    "target": pc.get("gene"),
                    "log2_fold_change": pc.get("log2_fold_change"),
                    "mammal_fraction": e.get("constrained_fraction"),
                }
            )
    out.sort(key=lambda e: (e["start"], e["end"]))
    return out


def non_overlapping(items: list[dict]) -> tuple[list[dict], int]:
    """Keep items in order, dropping any that overlaps the one kept before it."""
    kept: list[dict] = []
    last = -1
    dropped = 0
    for it in items:
        if it["start"] < last:
            dropped += 1
            continue
        kept.append(it)
        last = it["end"]
    return kept, dropped


def classify_elements(elements: list[dict], stats: list[IntervalStats | None]) -> list[dict]:
    rows = []
    for e, s in zip(elements, stats, strict=True):
        g = stats_dict(s)
        rows.append(
            {
                **e,
                "gnocchi": g,
                "case": case_of(
                    e.get("mammal_fraction"), g["fraction_above"] if g else None, ELEMENT_MAMMAL_MIN
                ),
            }
        )
    return rows


def tally_blocks(rows: list[dict]) -> dict:
    """Blocks and bases per tier and case, and the bp-weighted human-constrained fraction per tier."""
    out: dict[str, dict] = {}
    for r in rows:
        t = out.setdefault(
            r["tier"],
            {"blocks": 0, "measured_blocks": 0, "bp": 0, "measured_bp": 0, "human_constrained_bp": 0},
        )
        t["blocks"] += 1
        t["bp"] += r["length"]
        g = r["gnocchi"]
        if g:
            t["measured_blocks"] += 1
            t["measured_bp"] += g["bases"]
            t["human_constrained_bp"] += round(g["fraction_above"] * g["bases"])
        if r["case"]:
            c = t.setdefault("cases", {k: {"blocks": 0, "bp": 0} for k in CASE_ORDER})[r["case"]["case"]]
            c["blocks"] += 1
            c["bp"] += r["length"]
    for t in out.values():
        t["human_constrained_fraction"] = (
            round(t["human_constrained_bp"] / t["measured_bp"], 4) if t["measured_bp"] else None
        )
    return out


def tally_elements(rows: list[dict]) -> dict:
    """Per case: how many elements, how many name a gene, how many are strong in humans."""
    out = {k: {"elements": 0, "name_a_gene": 0, "strong": 0} for k in CASE_ORDER}
    unmeasured = 0
    for r in rows:
        if not r["case"]:
            unmeasured += 1
            continue
        c = out[r["case"]["case"]]
        c["elements"] += 1
        c["name_a_gene"] += 1 if r.get("target") else 0
        c["strong"] += 1 if r["gnocchi"] and r["gnocchi"]["strong"] else 0
    for c in out.values():
        c["name_a_gene_share"] = round(c["name_a_gene"] / c["elements"], 3) if c["elements"] else None
    return {"by_case": out, "unmeasured": unmeasured}


def aggregate(stats: list[IntervalStats]) -> dict:
    """Track-weighted summary of one set of intervals (a control set, a VISTA class).

    Gnocchi is one value per kilobase, and the reader counts every kilobase an interval
    touches, so `track_bases` exceeds the summed length for intervals shorter than 1 kb
    (exons, small elements); `fraction_above` is then the share of touched kilobases.
    """
    bases = sum(s.bases for s in stats)
    above = sum(s.above for s in stats)
    maxima = [s.maximum for s in stats if s.bases]
    return {
        "intervals": len(stats),
        "track_bases": bases,
        "fraction_above": round(above / bases, 4) if bases else None,
        "share_intervals_constrained": (
            round(sum(m >= GNOCCHI_THRESHOLD for m in maxima) / len(maxima), 3) if maxima else None
        ),
        "share_intervals_strong": (
            round(sum(m >= GNOCCHI_STRONG for m in maxima) / len(maxima), 3) if maxima else None
        ),
    }


def coding_control(chrom: str) -> tuple[list[tuple[int, int]], list[tuple[int, int]]] | None:
    """Merged canonical CDS segments and introns of the chromosome's coding genes, if GENCODE is local."""
    from genomeos.genome.annotation import Annotation, default_gencode

    gff = default_gencode({chrom})
    if gff is None:
        return None
    ann = Annotation.from_gff3(gff, {chrom})
    cds: list[tuple[int, int]] = []
    introns: list[tuple[int, int]] = []
    for g in ann.protein_coding():
        if g.locus.chrom != chrom:
            continue
        ts = list(g.transcripts.values())
        canon = [t for t in ts if "Ensembl_canonical" in t.tags] or sorted(
            ts, key=lambda t: -(t.locus.end - t.locus.start)
        )
        if not canon:
            continue
        t = canon[0]
        cds += [(c.start, c.end) for c, _ in t.cds]
        ex = sorted((e.start, e.end) for e in t.exons)
        introns += [(a[1], b[0]) for a, b in zip(ex, ex[1:], strict=False) if b[0] > a[1]]
    return merge(cds), merge(introns)


def vista_rows(chrom: str, results_dir: Path = RESULTS_DIR) -> list[dict]:
    r = load_result(f"vista_{chrom}", results_dir) or {}
    rows = [
        {"id": x["id"], "start": x["start"], "end": x["end"], "status": x["status"]}
        for x in r.get("rows", [])
        if x.get("status") in ("positive", "negative")
    ]
    rows.sort(key=lambda x: (x["start"], x["end"]))
    return rows


def vista_by_status(rows: list[dict], stats: list[IntervalStats]) -> dict:
    by: dict[str, list[IntervalStats]] = {"positive": [], "negative": []}
    for r, s in zip(rows, stats, strict=True):
        by[r["status"]].append(s)
    return {k: aggregate(v) for k, v in by.items()}


def known_value(kv: dict, ccres: list | None = None) -> dict:
    """Both axes at one named value: the base, its kilobase, and the registry elements next to it.

    A textbook value slot should read constrained across mammals at the base (the frame) and
    variable among people (the slot); the report says what each track has there.
    """
    chrom, pos = kv["chrom"], kv["pos"]
    near = sorted(
        (c for c in ccres or [] if c.end >= pos - 2000 and c.start <= pos + 2000), key=lambda c: c.start
    )
    # the base sits inside its flank, and the reader wants non-overlapping intervals: two reads
    base_iv = [(pos - 1, pos)]
    rest_iv = [(pos - 501, pos + 500)] + [(c.start, c.end) for c in near]
    labels = ["base", "flank_1kb"] + [f"{c.id} ({c.cls}, {c.start - pos:+,} bp)" for c in near]
    ps_base, _ = phylop_over_blocks(chrom, base_iv, PHYLOP_THRESHOLD)
    ps_rest, pcost = phylop_over_blocks(chrom, rest_iv, PHYLOP_THRESHOLD)
    gs_base, _ = gnocchi_over(chrom, base_iv)
    gs_rest, gcost = gnocchi_over(chrom, rest_iv)
    windows: dict[str, dict] = {}
    for label, p, g in zip(labels, [*ps_base, *ps_rest], [*gs_base, *gs_rest], strict=True):
        windows[label] = {
            "phylop": {
                "bases": p.bases,
                "mean": round(p.mean, 3) if p.mean is not None else None,
                "maximum": round(p.maximum, 3) if p.bases else None,
                "fraction_above": round(p.fraction_above, 3) if p.fraction_above is not None else None,
            },
            "gnocchi": stats_dict(g),
        }
    base = windows["base"]
    mammals = (
        "constrained across mammals at the base"
        if (base["phylop"]["fraction_above"] or 0) >= 1
        else "not constrained across mammals at the base"
    )
    humans = (
        "its kilobase is not scored by Gnocchi"
        if base["gnocchi"] is None
        else f"its kilobase reads Z {base['gnocchi']['maximum']} among people"
    )
    return {
        **kv,
        "windows": windows,
        "reading": f"{mammals}; {humans}",
        "cost": {"phylop": pcost, "gnocchi": gcost},
    }


def build(chrom: str, progress=None, results_dir: Path = RESULTS_DIR, controls: bool = True) -> dict:
    """Both axes over one chromosome's UNKNOWN blocks, attributed elements and VISTA elements."""
    budget = load_result(f"budget_{chrom}", results_dir)
    if not budget:
        raise FileNotFoundError(f"no budget_{chrom} result; run genomeos budget --chrom {chrom}")
    t0 = time.time()
    cost: dict[str, Any] = {}
    blocks = sorted(budget["blocks"], key=lambda b: b["start"])
    measured = [i for i, b in enumerate(blocks) if b["class"] != "gap"]
    stats: list[IntervalStats | None] = [None] * len(blocks)
    got, cost["blocks"] = gnocchi_over(
        chrom, [(blocks[i]["start"], blocks[i]["end"]) for i in measured], progress
    )
    for i, s in zip(measured, got, strict=True):
        stats[i] = s
    block_rows = classify_blocks(blocks, stats)

    elements, dropped = non_overlapping(elements_of(chrom, results_dir))
    element_rows: list[dict] = []
    if elements:
        need = [i for i, e in enumerate(elements) if e["mammal_fraction"] is None]
        if need:
            ps, cost["elements_phylop"] = phylop_over_blocks(
                chrom, [(elements[i]["start"], elements[i]["end"]) for i in need], PHYLOP_THRESHOLD
            )
            for i, p in zip(need, ps, strict=True):
                elements[i]["mammal_fraction"] = (
                    round(p.fraction_above, 4) if p.fraction_above is not None else None
                )
        es, cost["elements"] = gnocchi_over(chrom, [(e["start"], e["end"]) for e in elements])
        element_rows = classify_elements(elements, es)

    vista: dict[str, Any] = {}
    vrows, vdropped = non_overlapping(vista_rows(chrom, results_dir))
    if vrows:
        vs, cost["vista"] = gnocchi_over(chrom, [(r["start"], r["end"]) for r in vrows])
        vista = {**vista_by_status(vrows, vs), "overlapping_dropped": vdropped}

    ctrl: dict[str, Any] = {}
    if controls:
        cc = coding_control(chrom)
        if cc:
            cds, introns = cc
            cs, cost["coding_control"] = gnocchi_over(chrom, cds)
            ins, cost["intron_control"] = gnocchi_over(chrom, introns)
            ctrl = {"canonical_cds": aggregate(cs), "canonical_introns": aggregate(ins)}
    per_tier = tally_blocks(block_rows)
    ctrl["unknown_by_tier"] = {t: v["human_constrained_fraction"] for t, v in per_tier.items()}

    known = [known_value(kv, _ccres(chrom)) for kv in KNOWN_VALUES if kv["chrom"] == chrom]
    out = {
        "chrom": chrom,
        "thresholds": {
            "phylop": PHYLOP_THRESHOLD,
            "mammal_min_fraction_blocks": CONSTRAINED_MIN,
            "mammal_min_fraction_elements": ELEMENT_MAMMAL_MIN,
            "gnocchi": GNOCCHI_THRESHOLD,
            "gnocchi_strong": GNOCCHI_STRONG,
            "human_min_fraction": HUMAN_MIN_FRACTION,
        },
        "sources": {"mammals": PHYLOP_241_URL, "humans": GNOCCHI_URL},
        "evidence": EVIDENCE,
        "cases": CASE_LABELS,
        "blocks": block_rows,
        "by_tier": per_tier,
        "elements": element_rows,
        "elements_overlapping_dropped": dropped,
        "by_element_case": tally_elements(element_rows),
        "vista": vista,
        "controls": ctrl,
        "known_values": known,
        "cost": {**cost, "seconds": round(time.time() - t0, 1)},
    }
    return out


def _ccres(chrom: str):
    try:
        from genomeos.genome.regulatory import load_ccres

        return load_ccres(chrom)
    except (FileNotFoundError, ImportError, TypeError):
        return None


def run_and_save(chrom: str, progress=None, results_dir: Path = RESULTS_DIR, controls: bool = True) -> dict:
    out = build(chrom, progress=progress, results_dir=results_dir, controls=controls)
    save_result(f"variation_{chrom}", out, results_dir)
    return out


def vista_genome(results_dir: Path = RESULTS_DIR, progress=None) -> dict:
    """VISTA positives against negatives on the human axis over every chromosome with a result."""
    t0 = time.time()
    by: dict[str, list[IntervalStats]] = {"positive": [], "negative": []}
    cost = {"requests": 0, "mb_fetched": 0.0}
    chroms = 0
    for f in sorted(glob.glob(str(results_dir / "vista_chr*.json"))):
        chrom = Path(f).stem.split("_", 1)[1]
        rows, _ = non_overlapping(vista_rows(chrom, results_dir))
        if not rows:
            continue
        chroms += 1
        st, c = gnocchi_over(chrom, [(r["start"], r["end"]) for r in rows])
        cost["requests"] += c["requests"]
        cost["mb_fetched"] = round(cost["mb_fetched"] + c["mb_fetched"], 2)
        for r, s in zip(rows, st, strict=True):
            by[r["status"]].append(s)
        if progress:
            progress(f"{chrom}: {len(rows)} VISTA elements")
    out = {
        "chromosomes": chroms,
        **{k: aggregate(v) for k, v in by.items()},
        "thresholds": {"gnocchi": GNOCCHI_THRESHOLD, "gnocchi_strong": GNOCCHI_STRONG},
        "evidence": EVIDENCE["humans"],
        "cost": {**cost, "seconds": round(time.time() - t0, 1)},
    }
    return out
