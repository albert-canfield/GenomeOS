# SPDX-License-Identifier: AGPL-3.0-or-later
"""Phylogenetic profiling between libraries: which modules were gained and lost together.

The conversation behind area J asked for cross-references between libraries inferred from
covariance across the tree: two modules whose members are present in the same species and
absent from the same species depend on each other, or were built at the same time. The
presence matrix from the Compara stream (`origin_presence_genome_wide`: one bit per gene
per species, 375 species ordered by clade) makes that a small computation.

Per library: the share of its placed members with an orthologue in each species, a profile
over the species; and its gain curve, the share of members whose origin is at or below each
stratum. Between libraries: Pearson correlation of the two profiles over the species, a
dependency candidate where it exceeds what the genome's own profile explains. Everything is
`inferred`; a high correlation says two libraries have the same history, not that one calls
the other, and the Reactome hierarchy the libraries were built from is the check.
"""

from __future__ import annotations

import math
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

from genomeos.knowledge.homology import RANK, STRATA
from genomeos.results import RESULTS_DIR, load_result, save_result

MIN_MEMBERS = 8
NULL_DRAWS = 200  # age-matched random pairs per library pair
NULL_SEED = 7
EVIDENCE = {
    "profile": "curated: presence of an orthologue per species from Ensembl Compara (release 116)",
    "correlation": (
        "inferred: two libraries' presence profiles over the species agree beyond the genome's own"
    ),
}


def unpack(mask: str, n: int) -> list[int]:
    """The hex bitmask of one gene as a list of 0/1 over the n species, left to right."""
    if not mask:
        return [0] * n
    bits = bin(int(mask, 16))[2:].zfill(len(mask) * 4)
    return [int(b) for b in bits[-n:]] if n else []


def profile(masks: list[str], n: int) -> list[float]:
    """Share of the given genes present in each species."""
    if not masks:
        return [0.0] * n
    tot = [0] * n
    for m in masks:
        for i, b in enumerate(unpack(m, n)):
            tot[i] += b
    return [round(t / len(masks), 4) for t in tot]


def pearson(a: list[float], b: list[float]) -> float | None:
    n = len(a)
    if n < 3:
        return None
    ma, mb = sum(a) / n, sum(b) / n
    sa = math.sqrt(sum((x - ma) ** 2 for x in a))
    sb = math.sqrt(sum((x - mb) ** 2 for x in b))
    if sa == 0 or sb == 0:
        return None
    return round(sum((x - ma) * (y - mb) for x, y in zip(a, b, strict=True)) / (sa * sb), 4)


def gain_curve(origins: list[str]) -> dict[str, float]:
    """Cumulative share of genes whose origin is at or before each stratum, deep to shallow."""
    if not origins:
        return {}
    ranks = sorted(RANK[o] for o in origins if o in RANK)
    out: dict[str, float] = {}
    for s in STRATA:
        out[s] = round(sum(1 for r in ranks if r <= RANK[s]) / len(origins), 4)
    return out


def build(
    presence: dict, origin: dict, members: dict[str, list[str]], min_members: int = MIN_MEMBERS
) -> dict[str, Any]:
    species: list[str] = presence["species"]
    n = len(species)
    masks: dict[str, str] = presence["genes"]
    genes: dict[str, dict] = origin.get("genes") or {}
    genome_profile = profile(list(masks.values()), n)
    libraries: dict[str, dict] = {}
    for lib, syms in members.items():
        here = [s for s in syms if s in masks]
        if len(here) < min_members:
            continue
        prof = profile([masks[s] for s in here], n)
        libraries[lib] = {
            "members_placed": len(here),
            "profile": prof,
            "gain_curve": gain_curve([genes[s]["origin"] for s in here if s in genes]),
            "against_genome": pearson(prof, genome_profile),
        }
    # residual profiles: what a library does beyond the genome's own presence pattern
    resid = {
        lib: [p - g for p, g in zip(v["profile"], genome_profile, strict=True)]
        for lib, v in libraries.items()
    }
    pairs = []
    names = sorted(libraries)
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            r = pearson(resid[a], resid[b])
            if r is not None:
                pairs.append({"libraries": [a, b], "residual_correlation": r})
    pairs.sort(key=lambda x: -x["residual_correlation"])
    # per stratum: which species carry what share of the genome, for the reader
    by_stratum: dict[str, list[str]] = defaultdict(list)
    for sp in species:
        by_stratum[presence["strata"].get(sp) or "unplaced"].append(sp)
    return {
        "species": n,
        "species_by_stratum": {s: len(v) for s, v in by_stratum.items()},
        "genome_profile": genome_profile,
        "libraries": libraries,
        "pairs": pairs,
        "top_pairs": pairs[:25],
        "bottom_pairs": pairs[-10:],
        "min_members": min_members,
        "evidence": EVIDENCE,
    }


def run_and_save(results_dir: Path = RESULTS_DIR, null: bool = True, draws: int = NULL_DRAWS) -> dict:
    presence = load_result("origin_presence_genome_wide", results_dir)
    origin = load_result("origin_genome_wide", results_dir)
    if not presence or not origin:
        raise FileNotFoundError("run scripts/origin_genome_wide.py first (origin and presence results)")
    from genomeos.knowledge.homology import library_members

    members = library_members() or {}
    out = build(presence, origin, members)
    out["pairs"] = out["pairs"][:400]  # the full list is derivable; keep the file small
    if null:
        t0 = time.time()
        nl = pair_null(presence, origin, members, draws=draws)
        nl["seconds"] = round(time.time() - t0, 1)
        out["null"] = nl
    save_result("profiling_genome_wide", out, results_dir)
    return out


# --------------------------------------------------------------------------------------------- null


def species_masks(presence: dict) -> tuple[list[str], list[int], dict[str, int]]:
    """Each species as one integer bitmask over genes: bit g is set when gene g has an orthologue there.

    A gene set is then one integer too, and a set's presence count in a species is a popcount of the
    AND, so a profile over 355 species costs 355 popcounts instead of a loop over thousands of genes.
    """
    species: list[str] = presence["species"]
    n = len(species)
    genes = sorted(presence["genes"])
    index = {g: i for i, g in enumerate(genes)}
    masks = [0] * n
    for g, i in index.items():
        bits = unpack(presence["genes"][g], n)
        for j, b in enumerate(bits):
            if b:
                masks[j] |= 1 << i
    return species, masks, index


def set_mask(symbols, index: dict[str, int]) -> int:
    m = 0
    for s_ in symbols:
        i = index.get(s_)
        if i is not None:
            m |= 1 << i
    return m


def mask_profile(gene_set: int, masks: list[int], keep: list[int] | None = None) -> list[float]:
    size = gene_set.bit_count()
    cols = keep if keep is not None else range(len(masks))
    return [(masks[j] & gene_set).bit_count() / size if size else 0.0 for j in cols]


def matched_sampler(origins: dict[str, str], index: dict[str, int], rng):
    """Draw a random gene set with the same origin-stratum make-up as a given set."""
    pool: dict[str, list[int]] = defaultdict(list)
    for g, i in index.items():
        if g in origins:
            pool[origins[g]].append(i)

    def draw(symbols) -> int:
        need: dict[str, int] = defaultdict(int)
        for s_ in symbols:
            if s_ in origins and s_ in index:
                need[origins[s_]] += 1
        m = 0
        for stratum, k in need.items():
            for i in rng.sample(pool[stratum], min(k, len(pool[stratum]))):
                m |= 1 << i
        return m

    return draw


def pair_null(
    presence: dict,
    origin: dict,
    members: dict[str, list[str]],
    min_members: int = MIN_MEMBERS,
    draws: int = NULL_DRAWS,
    seed: int = NULL_SEED,
) -> dict[str, Any]:
    """Library pairs held against age-matched random gene sets, on members exclusive to each library.

    Two libraries' presence profiles correlate for three reasons besides a shared history: shared
    members (the same genes under two names), shared age (two ancient libraries look alike because
    both are ancient) and the bacterial species dominating the profile. For each pair the observed
    residual correlation is taken on the members exclusive to each side, over all species and over
    eukaryotic species alone, and compared with `draws` pairs of random sets matching each side's
    origin-stratum make-up; z is (observed − null mean) / null sd, and q the Benjamini–Hochberg
    false-discovery rate of the normal upper-tail p over all pairs.
    """
    import random
    from itertools import combinations
    from statistics import NormalDist

    rng = random.Random(seed)
    species, masks, index = species_masks(presence)
    strata = presence.get("strata") or {}
    euk = [j for j, sp in enumerate(species) if strata.get(sp) not in (None, "Life")]
    origins = {g: v["origin"] for g, v in (origin.get("genes") or {}).items()}
    everyone = set_mask(index, index)
    genome_all = mask_profile(everyone, masks)
    genome_euk = [genome_all[j] for j in euk]
    sets = {lib: {s_ for s_ in syms if s_ in index} for lib, syms in members.items()}
    sets = {lib: s_ for lib, s_ in sets.items() if len(s_) >= min_members}
    draw = matched_sampler(origins, index, rng)

    def residual(profile: list[float], base: list[float]) -> list[float]:
        return [p - b for p, b in zip(profile, base, strict=True)]

    rows = []
    for a, b in combinations(sorted(sets), 2):
        only_a, only_b = sets[a] - sets[b], sets[b] - sets[a]
        shared = len(sets[a] & sets[b])
        row: dict[str, Any] = {
            "libraries": [a, b],
            "shared_members": shared,
            "exclusive_members": [len(only_a), len(only_b)],
        }
        if len(only_a) < min_members or len(only_b) < min_members:
            row["note"] = "too few exclusive members: the pair is mostly the same genes"
            rows.append(row)
            continue
        ma, mb = set_mask(only_a, index), set_mask(only_b, index)
        obs = pearson(
            residual(mask_profile(ma, masks), genome_all), residual(mask_profile(mb, masks), genome_all)
        )
        obs_euk = pearson(
            residual(mask_profile(ma, masks, euk), genome_euk),
            residual(mask_profile(mb, masks, euk), genome_euk),
        )
        null: list[float] = []
        for _ in range(draws):
            ra, rb = draw(only_a), draw(only_b)
            r = pearson(
                residual(mask_profile(ra, masks), genome_all), residual(mask_profile(rb, masks), genome_all)
            )
            if r is not None:
                null.append(r)
        mean = sum(null) / len(null) if null else None
        sd = (sum((x - mean) ** 2 for x in null) / (len(null) - 1)) ** 0.5 if len(null) > 1 else None
        z = (obs - mean) / sd if obs is not None and mean is not None and sd else None
        row.update(
            {
                "r_exclusive": obs,
                "r_exclusive_eukaryotes": obs_euk,
                "null_mean": round(mean, 4) if mean is not None else None,
                "null_sd": round(sd, 4) if sd else None,
                "z": round(z, 2) if z is not None else None,
                "p": round(1 - NormalDist().cdf(z), 6) if z is not None else None,
            }
        )
        rows.append(row)
    tested = sorted((r for r in rows if r.get("p") is not None), key=lambda r: r["p"])
    m = len(tested)
    prev = 1.0
    for k in range(m - 1, -1, -1):  # Benjamini–Hochberg, monotone from the largest p down
        q = min(prev, tested[k]["p"] * m / (k + 1))
        tested[k]["q"] = round(q, 6)
        prev = q
    return {
        "pairs": rows,
        "tested": m,
        "significant_q05": [r for r in tested if r["q"] <= 0.05],
        "draws": draws,
        "seed": seed,
        "species": len(species),
        "eukaryotic_species": len(euk),
        "evidence": (
            "inferred: residual presence-profile correlation on exclusive members, against "
            "age-matched random gene sets; Benjamini–Hochberg over all pairs"
        ),
    }
