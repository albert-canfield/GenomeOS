# SPDX-License-Identifier: AGPL-3.0-or-later
"""One locus across species: does the module reconstruct in mouse, chicken and fish?

The last step of ROADMAP area J is the experiment the conversation proposed: take one
developmental locus with every layer available, read it in human and in other vertebrates,
and see whether the same grammar is there. Two loci are the first candidates. The ZRS
(VISTA hs2496, chr7:156,791,086-156,791,875) is the limb enhancer of SHH, a megabase away
inside an intron of LMBR1, conserved to fish and the textbook long-range pointer. The
HERC2 intron-86 enhancer of OCA2 carries rs12913832, Albert's eye-colour value slot.

For each locus and each species, Ensembl Compara's pairwise LASTZ alignment of the human
region is fetched (one call, under a second), giving the other species' locus, the
identity over aligned columns and the share of the human region that aligns at all. The
JASPAR scan (genome/motifs.py) then runs over both sides of each alignment, and a factor's
motif counts as conserved only where it hits the *same aligned site* in both: identity alone
makes nearly every motif "present" in a close species (a first version called 170 factors
conserved across the ZRS), so each human hit is mapped through the alignment columns. The
factors holding their site in every species are the locus's grammar. Both
constraint axes (Zoonomia phyloP, gnomAD Gnocchi) are read over the human region as in
`genomeos variation`. Everything fetched is cached under data/knowledge/across.

Evidence: the alignment is `curated` (Ensembl LASTZ_NET), a motif hit `predicted`, a
"conserved grammar" call `inferred` from the two. What the alignment cannot say is
whether the element is active in the other species; VISTA's mouse assay says it for the
ZRS (positive, limb), and for the rest that is the reader's job.
"""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

from genomeos.results import RESULTS_DIR, save_result

CACHE = Path("data/knowledge/across")
REST = (
    "https://rest.ensembl.org/alignment/region/homo_sapiens/{region}"
    "?method=LASTZ_NET&species_set=homo_sapiens&species_set={species}"
)
SPECIES = ("mus_musculus", "monodelphis_domestica", "gallus_gallus", "xenopus_tropicalis", "danio_rerio")
LOCI = {
    "ZRS": {
        "chrom": "chr7",
        "start": 156_791_086,
        "end": 156_791_875,
        "gene": "SHH",
        "note": (
            "VISTA hs2496, limb positive; the SHH limb enhancer inside LMBR1 intron 5, 1 Mb from its target"
        ),
    },
    "HERC2_OCA2": {
        "chrom": "chr15",
        "start": 28_120_106,
        "end": 28_121_063,
        "gene": "OCA2",
        "note": (
            "the HERC2 intron-86 enhancer of OCA2 around rs12913832 (chr15:28,120,472), eye colour's slot"
        ),
    },
}
EVIDENCE = {
    "alignment": "curated: Ensembl Compara LASTZ_NET pairwise alignment of the human region",
    "motif": "predicted: JASPAR 2026 log-odds hit at >= 85% of the matrix range",
    "grammar": "inferred: a factor whose motif holds the same aligned site in human and every species",
}


def _get(url: str, cache_key: str, timeout: int = 180) -> Any:
    CACHE.mkdir(parents=True, exist_ok=True)
    p = CACHE / f"{cache_key}.json"
    if p.exists():
        return json.loads(p.read_text())
    req = urllib.request.Request(
        url, headers={"Content-Type": "application/json", "User-Agent": "GenomeOS/0.9 (across)"}
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310
        d = json.load(r)
    p.write_text(json.dumps(d))
    return d


def fetch_alignment(chrom: str, start: int, end: int, species: str) -> list[dict]:
    """The LASTZ_NET blocks of the human region against one species (1-based inclusive region)."""
    region = f"{chrom.removeprefix('chr')}:{start + 1}-{end}"
    url = REST.format(region=urllib.parse.quote(region), species=species)
    d = _get(url, f"{chrom}_{start}_{end}_{species}")
    return d if isinstance(d, list) else [d]


def summarise_alignment(blocks: list[dict], species: str, length: int) -> dict[str, Any]:
    """Identity over aligned columns, share of the human region covered, the other locus, the sequences."""
    human_cols = other_cols = matches = 0
    human_seq_parts: list[str] = []
    other_seq_parts: list[str] = []
    human_blocks: list[tuple[int, int]] = []  # (offset in the concatenated human bases, 0-based genome start)
    offset = 0
    loci: list[str] = []
    for b in blocks:
        h = next((a for a in b.get("alignments", []) if a["species"] == "homo_sapiens"), None)
        o = next((a for a in b.get("alignments", []) if a["species"] == species), None)
        if not h or not o:
            continue
        for x, y in zip(h["seq"], o["seq"], strict=False):
            if x != "-":
                human_cols += 1
            if y != "-":
                other_cols += 1
            if x != "-" and y != "-" and x.upper() == y.upper():
                matches += 1
        # gapped, column for column: the same-site test maps positions through these columns
        human_seq_parts.append(h["seq"])
        other_seq_parts.append(o["seq"])
        human_blocks.append((offset, h["start"] - 1))
        offset += len(h["seq"].replace("-", ""))
        loci.append(f"{o['seq_region']}:{o['start']}-{o['end']}({'+' if o['strand'] == 1 else '-'})")
    aligned = min(human_cols, other_cols)
    return {
        "species": species,
        "blocks": len(loci),
        "locus": ", ".join(loci) or None,
        "human_bases_aligned": human_cols,
        "coverage": round(human_cols / length, 3) if length else None,
        "identity": round(matches / human_cols, 3) if human_cols else None,
        "other_sequence": "".join(other_seq_parts),
        "human_aligned_sequence": "".join(human_seq_parts),
        "human_blocks": human_blocks,
        "aligned_columns": aligned,
    }


SAME_SITE_WINDOW = 15  # a motif counts as the same site when its start maps within this many bases


def column_map(gapped: str) -> list[int]:
    """For each ungapped base of a gapped sequence, its column index."""
    return [i for i, c in enumerate(gapped) if c != "-"]


def base_at_column(gapped: str) -> dict[int, int]:
    """Column index to ungapped base index."""
    return {i: n for n, i in enumerate(column_map(gapped))}


def to_genome(pos: int, blocks: list[tuple[int, int]]) -> int | None:
    """A position in the concatenated human bases of an alignment, as a 0-based genome coordinate."""
    for off, start0 in reversed(blocks):
        if pos >= off:
            return start0 + (pos - off)
    return None


MIN_CORE_COVERAGE = 0.5  # a species aligning over at least this share of the element must hold a site


def cluster_sites(held: list[tuple[int, int, str, str, str]], core: list[str]) -> list[dict]:
    """Held hits (genome start, end, factor, family, species) clustered into sites where they overlap.

    A site records every species holding a hit in it; it is held *everywhere* when all core species do.
    Hits of one site may come from different matrices in different species (a HOX site read by HOXA9 in
    mouse and CDX1 in chicken): what is conserved is the site, not the matrix.
    """
    sites: list[dict] = []
    for a, b, tf, fam, sp in sorted(held, key=lambda t: (t[0], t[1], t[2], t[4])):
        if sites and a < sites[-1]["end"]:
            st = sites[-1]
            st["end"] = max(st["end"], b)
        else:
            st = {"start": a, "end": b, "factors": {}, "families": set(), "species": set()}
            sites.append(st)
        st["factors"].setdefault(tf, set()).add(sp)
        st["families"].add(fam)
        st["species"].add(sp)
    out = []
    for st in sites:
        everywhere = bool(core) and set(core) <= st["species"]
        out.append(
            {
                "start": st["start"],
                "end": st["end"],
                "families": sorted(st["families"]),
                "factors": sorted(st["factors"]),
                "species": sorted(st["species"]),
                "everywhere": everywhere,
                # a factor held at this site by every core species on its own
                "factors_everywhere": sorted(f for f, sps in st["factors"].items() if set(core) <= sps),
            }
        )
    return out


def conserved_grammar(
    pairs: dict[str, tuple],
    motifs,
    window: int = SAME_SITE_WINDOW,
    families: dict[str, dict[str, str]] | None = None,
    core: list[str] | None = None,
) -> dict[str, Any]:
    """Factors hitting the human element, and the sites whose motif the other species keep in place.

    Sequence identity alone makes almost every motif "present" in a close species, so a hit counts as
    held only where the other species' hit sits at the aligned position: the human hit's start is mapped
    through the gapped alignment columns and must land within `window` bases of the other hit. Each held
    human hit is placed on the genome through that species' own alignment blocks, and held hits from all
    species are clustered into sites. A site is held everywhere when every `core` species holds it (by
    default every species compared). Only each factor's best hit per sequence is scanned, so a site held
    by a weaker hit can be missed: the call is conservative, never inflated.
    """
    from genomeos.genome.motifs import family_unit, load_families, scan

    families = families if families is not None else load_families()
    widths = {m.name.upper(): m.width for m in motifs}
    per_species: dict[str, dict[str, int]] = {}
    human_hits: dict[str, tuple[float, int]] = {}
    held: list[tuple[int, int, str, str, str]] = []
    for sp, pair in pairs.items():
        human_gapped, other_gapped = pair[0], pair[1]
        blocks = pair[2] if len(pair) > 2 else None
        h_seq, o_seq = human_gapped.replace("-", ""), other_gapped.replace("-", "")
        if not h_seq or not o_seq:
            continue
        hits = scan(motifs, [h_seq, o_seq])
        h_cols, o_base = column_map(human_gapped), base_at_column(other_gapped)
        for tf, per_seq in hits.items():
            if 0 not in per_seq:
                continue
            score, pos, _ = per_seq[0]
            if tf not in human_hits or score > human_hits[tf][0]:
                human_hits[tf] = (score, pos)
            if 1 not in per_seq:
                continue
            col = h_cols[pos] if pos < len(h_cols) else None
            if col is None:
                continue
            near = min(o_base, key=lambda c: abs(c - col)) if o_base else None
            expected = o_base.get(near) if near is not None else None
            if expected is None or abs(per_seq[1][1] - expected) > window:
                continue
            per_species.setdefault(tf, {})[sp] = per_seq[1][1]
            g = to_genome(pos, blocks) if blocks else pos
            if g is not None:
                held.append((g, g + widths.get(tf, 10), tf, family_unit(tf, families), sp))
    species = [sp for sp, pr in pairs.items() if pr[0].replace("-", "") and pr[1].replace("-", "")]
    core = [sp for sp in (core if core is not None else species) if sp in species]
    sites = cluster_sites(held, core)
    held_factor = {f for st in sites if st["everywhere"] for f in st["factors_everywhere"]}
    rows = []
    for tf, (score, pos) in human_hits.items():
        rows.append(
            {
                "factor": tf,
                "family": family_unit(tf, families),
                "human_score": score,
                "human_position": pos,
                "same_site_in": sorted(per_species.get(tf, {})),
                "everywhere": tf in held_factor,
            }
        )
    rows.sort(key=lambda r: (-len(r["same_site_in"]), -r["human_score"]))
    everywhere_sites = [st for st in sites if st["everywhere"]]
    held_families: dict[str, list[str]] = {}
    for st in everywhere_sites:
        for f in st["factors_everywhere"]:
            held_families.setdefault(family_unit(f, families), [])
            if f not in held_families[family_unit(f, families)]:
                held_families[family_unit(f, families)].append(f)
    return {
        "core_species": core,
        "species_compared": species,
        "same_site_window": window,
        "factors_in_human": len(rows),
        "families_in_human": len({r["family"] for r in rows}),
        "sites_held": sites,
        "sites_everywhere": everywhere_sites,
        "sites_everywhere_count": len(everywhere_sites),
        "families_everywhere": dict(sorted(held_families.items(), key=lambda kv: (-len(kv[1]), kv[0]))),
        "everywhere": [r["factor"] for r in rows if r["everywhere"]],
        "rows": rows,  # every factor that hit human: a capped list made absent and unranked look alike
    }


def constraint(chrom: str, start: int, end: int) -> dict[str, Any]:
    from genomeos.attribution.constraint import PHYLOP_THRESHOLD, phylop_over_blocks
    from genomeos.attribution.variation import gnocchi_over, stats_dict

    ps, _ = phylop_over_blocks(chrom, [(start, end)], PHYLOP_THRESHOLD)
    gs, _ = gnocchi_over(chrom, [(start, end)])
    p = ps[0]
    return {
        "phylop": {
            "bases": p.bases,
            "mean": round(p.mean, 3) if p.mean is not None else None,
            "maximum": round(p.maximum, 3) if p.bases else None,
            "fraction_above": round(p.fraction_above, 3) if p.fraction_above is not None else None,
        },
        "gnocchi": stats_dict(gs[0]),
    }


def human_sequence(chrom: str, start: int, end: int) -> str | None:
    from genomeos.coords import Locus
    from genomeos.genome.fetch import REFERENCE
    from genomeos.genome.genome import Genome

    fa = REFERENCE / f"{chrom}.fa.gz"
    if not fa.exists():
        return None
    return str(Genome.from_fasta(fa).fetch(Locus(chrom, start, end)))


def vista_loci(results_dir: Path = RESULTS_DIR) -> dict[str, dict]:
    """Every VISTA element read genome-wide, as a locus: id, place, status, tissues, constraint."""
    import glob

    from genomeos.results import load_result

    out: dict[str, dict] = {}
    for f in sorted(glob.glob(str(results_dir / "vista_chr*.json"))):
        r = load_result(Path(f).stem, results_dir) or {}
        for x in r.get("rows", []):
            if x.get("status") not in ("positive", "negative"):
                continue
            target = (x.get("predicted") or {}).get("gene") or (x.get("inferred") or {}).get("gene")
            out[x["id"]] = {
                "chrom": r["chrom"],
                "start": x["start"],
                "end": x["end"],
                "gene": target,
                "status": x["status"],
                "tissues": x.get("tissues") or [],
                "constrained_fraction": x.get("constrained_fraction"),
                "note": f"VISTA {x['id']}, {x['status']}"
                + (f" ({', '.join(x.get('tissues') or [])})" if x.get("tissues") else ""),
            }
    return out


def build(
    name: str | dict,
    species: tuple[str, ...] = SPECIES,
    motifs=None,
    progress=None,
    with_constraint: bool = True,
) -> dict[str, Any]:
    if isinstance(name, dict):
        loc = name
        name = loc.get("id") or f"{loc['chrom']}_{loc['start']}"
    else:
        if name not in LOCI:
            raise KeyError(f"unknown locus {name}; known: {', '.join(LOCI)}")
        loc = LOCI[name]
    chrom, start, end = loc["chrom"], loc["start"], loc["end"]
    t0 = time.time()
    length = end - start
    per_species: dict[str, dict] = {}
    for sp in species:
        try:
            blocks = fetch_alignment(chrom, start, end, sp)
            per_species[sp] = summarise_alignment(blocks, sp, length)
        except OSError as e:
            per_species[sp] = {"species": sp, "error": str(e)[:120], "coverage": 0.0, "identity": None}
        if progress:
            v = per_species[sp]
            progress(f"{name} against {sp}: coverage {v.get('coverage')} identity {v.get('identity')}")
    human = human_sequence(chrom, start, end)
    grammar: dict[str, Any] = {}
    if human:  # the local reference confirms the locus; the scan reads the aligned sequences
        if motifs is None:
            from genomeos.genome.motifs import load_motifs

            motifs = load_motifs()
        pairs = {
            sp: (v["human_aligned_sequence"], v["other_sequence"], v.get("human_blocks") or [])
            for sp, v in per_species.items()
            if v.get("other_sequence") and v.get("human_aligned_sequence")
        }
        core = [sp for sp, v in per_species.items() if (v.get("coverage") or 0) >= MIN_CORE_COVERAGE]
        grammar = conserved_grammar(pairs, motifs, core=core)
    out = {
        "locus": name,
        **loc,
        "length": length,
        "species": {
            sp: {
                k: v
                for k, v in d.items()
                if k not in ("other_sequence", "human_aligned_sequence", "human_blocks")
            }
            for sp, d in per_species.items()
        },
        "aligned_in": [sp for sp, v in per_species.items() if (v.get("coverage") or 0) > 0],
        "grammar": grammar,
        "constraint": constraint(chrom, start, end) if with_constraint else None,
        "evidence": EVIDENCE,
        "cost": {"seconds": round(time.time() - t0, 1)},
    }
    return out


def run_and_save(name: str, results_dir: Path = RESULTS_DIR, progress=None) -> dict:
    out = build(name, progress=progress)
    save_result(f"across_{name}", out, results_dir)
    return out


# ------------------------------------------------------------------------------------------ panel
NEURAL_TISSUES = {"fb", "mb", "hb", "nt"}
PANEL_CACHE = CACHE / "panel"
MIN_CORE_SPECIES = 2  # "held" must mean more than one other genome


def fisher_greater(a: int, b: int, c: int, d: int) -> float:
    """One-sided Fisher exact p for a 2x2 table [[a, b], [c, d]]: P(X >= a) with margins fixed."""
    from math import comb

    n, row1, col1 = a + b + c + d, a + b, a + c
    denom = comb(n, col1)
    lo = a
    hi = min(row1, col1)
    return min(1.0, sum(comb(row1, k) * comb(n - row1, col1 - k) for k in range(lo, hi + 1)) / denom)


def mann_whitney_greater(x: list[float], y: list[float]) -> float | None:
    """One-sided Mann–Whitney p that x tends to exceed y (normal approximation with tie-averaged ranks)."""
    from statistics import NormalDist

    if not x or not y:
        return None
    allv = sorted((v, 0) for v in x) + sorted((v, 1) for v in y)
    allv.sort(key=lambda t: t[0])
    ranks = [0.0] * len(allv)
    i = 0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1][0] == allv[i][0]:
            j += 1
        for k in range(i, j + 1):
            ranks[k] = (i + j) / 2 + 1
        i = j + 1
    r1 = sum(r for r, (_, g) in zip(ranks, allv, strict=True) if g == 0)
    n1, n2 = len(x), len(y)
    u = r1 - n1 * (n1 + 1) / 2
    mu, sd = n1 * n2 / 2, (n1 * n2 * (n1 + n2 + 1) / 12) ** 0.5
    return round(1 - NormalDist().cdf((u - mu) / sd), 6) if sd else None


def bh(pvalues: list[float]) -> list[float]:
    """Benjamini–Hochberg q-values in the input order."""
    m = len(pvalues)
    order = sorted(range(m), key=lambda k: pvalues[k])
    q = [1.0] * m
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        k = order[rank]
        prev = min(prev, pvalues[k] * m / (rank + 1))
        q[k] = prev
    return q


def locus_group(tissues: list[str], status: str) -> str | None:
    if status == "negative":
        return "negative"
    t = set(tissues)
    if "lb" in t and not t & NEURAL_TISSUES:
        return "limb"
    if t & NEURAL_TISSUES and "lb" not in t:
        return "neural"
    return None


def select_panel(loci: dict[str, dict], per_group: int = 40, seed: int = 5) -> list[dict]:
    """Limb-only and neural-only positives, and negatives matched to the positives' constraint bins."""
    import random

    rng = random.Random(seed)

    def ok(v: dict) -> bool:
        cf = v.get("constrained_fraction") or 0
        return 0.3 <= cf <= 0.9 and 500 <= v["end"] - v["start"] <= 2000

    edges = (0.3, 0.45, 0.6, 0.75, 0.9001)

    def cbin(v: dict) -> int:
        cf = v["constrained_fraction"]
        return next(i for i in range(len(edges) - 1) if edges[i] <= cf < edges[i + 1])

    chosen: list[dict] = []
    for grp in ("limb", "neural"):
        pool = sorted(
            (
                {**v, "id": k, "group": grp}
                for k, v in loci.items()
                if ok(v) and locus_group(v["tissues"], v["status"]) == grp
            ),
            key=lambda v: v["id"],
        )
        rng.shuffle(pool)
        chosen += pool[:per_group]
    need: dict[int, int] = {}
    for v in chosen:
        need[cbin(v)] = need.get(cbin(v), 0) + 1
    negatives = sorted(
        ({**v, "id": k, "group": "negative"} for k, v in loci.items() if ok(v) and v["status"] == "negative"),
        key=lambda v: v["id"],
    )
    rng.shuffle(negatives)
    for bin_, k in sorted(need.items()):
        chosen += [v for v in negatives if cbin(v) == bin_][:k]
    return chosen


def panel_record(locus: dict, motifs) -> dict:
    """One locus through `build`, kept compact and cached: coverage, core species, strict held sites."""
    PANEL_CACHE.mkdir(parents=True, exist_ok=True)
    cache = PANEL_CACHE / f"{locus['id']}.json"
    if cache.exists():
        return json.loads(cache.read_text())
    r = build(locus, motifs=motifs, with_constraint=False)
    g = r.get("grammar") or {}
    strict = [st for st in g.get("sites_held", []) if st["everywhere"] and st["factors_everywhere"]]
    from genomeos.genome.motifs import family_unit, load_families

    fam = load_families()
    rec = {
        "id": locus["id"],
        "group": locus["group"],
        "chrom": locus["chrom"],
        "start": locus["start"],
        "end": locus["end"],
        "length": r["length"],
        "constrained_fraction": locus.get("constrained_fraction"),
        "coverage": {sp: v.get("coverage") for sp, v in r["species"].items()},
        "core_species": g.get("core_species", []),
        "strict_sites": len(strict),
        "families_held": sorted({family_unit(f, fam) for st in strict for f in st["factors_everywhere"]}),
        "families_hit": sorted({row["family"] for row in g.get("rows", [])}),
    }
    cache.write_text(json.dumps(rec))
    return rec


def panel_summary(records: list[dict], min_core: int = MIN_CORE_SPECIES) -> dict[str, Any]:
    """Per positive group against the negatives: held-site density, and which families are held more often."""
    usable = [r for r in records if len(r["core_species"]) >= min_core]
    groups: dict[str, list[dict]] = {}
    for r in usable:
        groups.setdefault(r["group"], []).append(r)
    neg = groups.get("negative", [])
    out: dict[str, Any] = {
        "loci": len(records),
        "usable": {k: len(v) for k, v in groups.items()},
        "min_core_species": min_core,
        "groups": {},
    }
    density = {k: [1000 * r["strict_sites"] / r["length"] for r in v] for k, v in groups.items()}
    for grp in ("limb", "neural"):
        pos = groups.get(grp, [])
        if not pos or not neg:
            continue
        fams = sorted({f for r in pos + neg for f in r["families_held"]})
        rows = []
        for f in fams:
            a = sum(1 for r in pos if f in r["families_held"])
            c = sum(1 for r in neg if f in r["families_held"])
            if a + c < 3:
                continue
            rows.append(
                {
                    "family": f,
                    "positives_held": a,
                    "negatives_held": c,
                    "positive_share": round(a / len(pos), 3),
                    "negative_share": round(c / len(neg), 3),
                    "p": fisher_greater(a, len(pos) - a, c, len(neg) - c),
                }
            )
        for row, q in zip(rows, bh([r["p"] for r in rows]), strict=True):
            row["q"] = round(q, 6)
        rows.sort(key=lambda r: r["p"])
        dp, dn = density.get(grp, []), density.get("negative", [])
        out["groups"][grp] = {
            "positives": len(pos),
            "negatives": len(neg),
            "density_median_positive": sorted(dp)[len(dp) // 2] if dp else None,
            "density_median_negative": sorted(dn)[len(dn) // 2] if dn else None,
            "density_p_greater": mann_whitney_greater(dp, dn),
            "families_tested": len(rows),
            "families_q05": [r for r in rows if r["q"] <= 0.05],
            "families": rows[:40],
        }
    return out
