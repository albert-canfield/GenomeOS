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
from collections import defaultdict
from pathlib import Path
from typing import Any

from genomeos.knowledge.homology import RANK, STRATA
from genomeos.results import RESULTS_DIR, load_result, save_result

MIN_MEMBERS = 8
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


def run_and_save(results_dir: Path = RESULTS_DIR) -> dict:
    presence = load_result("origin_presence_genome_wide", results_dir)
    origin = load_result("origin_genome_wide", results_dir)
    if not presence or not origin:
        raise FileNotFoundError("run scripts/origin_genome_wide.py first (origin and presence results)")
    from genomeos.knowledge.homology import library_members

    out = build(presence, origin, library_members() or {})
    out["pairs"] = out["pairs"][:400]  # the full list is derivable; keep the file small
    save_result("profiling_genome_wide", out, results_dir)
    return out
