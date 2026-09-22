# SPDX-License-Identifier: AGPL-3.0-or-later
"""Phylogenetic profiling at the gene level, against a null that holds prevalence fixed. No requests.

    uv run python scripts/profiling_genes.py [--libraries core.replication,core.splicing] [--no-save]

Area J profiled the libraries against each other and left the gene level open, with the condition
attached: a shuffled null before any pair is called a dependency. The condition is the whole
difficulty. Two genes present in all 355 species co-occur perfectly and mean nothing by it, and a null
that shuffles presence freely would make every such pair look extraordinary. So the null here holds
each gene's PREVALENCE fixed — the same number of species, drawn at random — and destroys only which
species they are. A pair beats it or it does not.

Read from `origin_presence_genome_wide` (a hex bitmask per gene over 355 species) and the library
catalogue. Within a library, every member pair scores by Jaccard over their presence sets; the library
scores by its mean pair Jaccard; and the same arithmetic runs over `DRAWS` prevalence-preserving
shuffles to say what that mean is worth.

**What this cannot do, stated before the numbers.** Species are not independent: two genes both lost
in the same clade share that loss once, not in every species of it, and a null that draws species
independently does not know that. It therefore *overstates* how surprising a co-occurrence is, and
every figure here is a ceiling on the evidence rather than the evidence. A species-tree-aware null
(Felsenstein contrasts, or a birth-death model along the tree) is what would settle a pair, and this
project does not hold one. Nothing here calls a pair a dependency; it says which pairs would survive
the weakest honest test, so a stronger one has somewhere to start.
"""

from __future__ import annotations

import argparse
import random
import time
from itertools import combinations
from statistics import mean
from typing import Any

from genomeos.lib import LIBRARIES, KnowledgeBase
from genomeos.results import load_result, save_result

DRAWS = 20  # prevalence-preserving shuffles per library
MIN_MEMBERS = 8  # the library-level run's own floor, kept so the two are comparable
MAX_PAIRS = 60_000  # a library of 350 members is 61k pairs; beyond that the mean is sampled
SEED = 11


def masks(presence: dict[str, Any]) -> tuple[dict[str, int], int]:
    """Gene -> presence bitmask as an int, and the number of species."""
    n = presence.get("species")
    n = n if isinstance(n, int) else len(n or [])
    return {g: int(h, 16) for g, h in (presence.get("genes") or {}).items() if h}, n


def jaccard(a: int, b: int) -> float:
    union = (a | b).bit_count()
    return ((a & b).bit_count() / union) if union else 0.0


def mean_pair_score(members: list[int], rng: random.Random) -> float | None:
    """Mean Jaccard over member pairs, sampled when the library is large enough to be quadratic."""
    n = len(members)
    if n < 2:
        return None
    total = n * (n - 1) // 2
    if total <= MAX_PAIRS:
        return mean(jaccard(a, b) for a, b in combinations(members, 2))
    return mean(jaccard(members[rng.randrange(n)], members[rng.randrange(n)]) for _ in range(MAX_PAIRS))


def shuffled(mask: int, species: int, rng: random.Random) -> int:
    """A mask with the same number of species set, drawn at random.

    Prevalence is what has to be held: a gene in every species co-occurs perfectly with any other such
    gene, and a null free to change how many species a gene is in would call that a discovery.
    """
    out = 0
    for i in rng.sample(range(species), mask.bit_count()):
        out |= 1 << i
    return out


def matched_draw(members: list[int], pool: list[int], rng: random.Random) -> list[int]:
    """A random gene set of the same size, each gene matched on prevalence to one of the members.

    This is the null the roadmap asked for and the one that matters: the drawn genes are REAL genes
    with real phylogenetic profiles, so the tree structure that makes species non-independent is
    present in the null as well as in the library. The prevalence shuffle cannot do that, which is why
    every library beats it.
    """
    out = []
    for m in members:
        want = m.bit_count()
        lo = max(0, int(want * 0.9) - 1)
        hi = int(want * 1.1) + 1
        near = [p for p in pool if lo <= p.bit_count() <= hi]
        out.append(rng.choice(near) if near else rng.choice(pool))
    return out


def profile_library(
    name: str,
    genes: list[str],
    by_gene: dict[str, int],
    species: int,
    pool: list[int] | None = None,
) -> dict | None:
    rng = random.Random(SEED)
    members = [by_gene[g] for g in genes if g in by_gene]
    if len(members) < MIN_MEMBERS:
        return None
    observed = mean_pair_score(members, rng)
    if observed is None:
        return None
    null = []
    for _ in range(DRAWS):
        drawn = [shuffled(m, species, rng) for m in members]
        score = mean_pair_score(drawn, rng)
        if score is not None:
            null.append(score)
    null_mean = mean(null) if null else None
    beats = sum(1 for x in null if x >= observed)
    matched: list[float] = []
    if pool:
        for _ in range(DRAWS):
            score = mean_pair_score(matched_draw(members, pool, rng), rng)
            if score is not None:
                matched.append(score)
    matched_mean = mean(matched) if matched else None
    matched_beats = sum(1 for x in matched if x >= observed)
    return {
        "library": name,
        "members": len(genes),
        "members_with_a_profile": len(members),
        "mean_pair_jaccard": round(observed, 4),
        "null_mean": round(null_mean, 4) if null_mean is not None else None,
        "excess_over_null": round(observed - null_mean, 4) if null_mean is not None else None,
        "null_draws_at_or_above_observed": beats,
        "matched_gene_null_mean": round(matched_mean, 4) if matched_mean is not None else None,
        "excess_over_matched_genes": (
            round(observed - matched_mean, 4) if matched_mean is not None else None
        ),
        "matched_draws_at_or_above_observed": matched_beats if matched else None,
        "draws": len(null),
        "median_prevalence": (round(mean(m.bit_count() for m in members) / species, 4) if species else None),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--libraries", default="")
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args(argv)

    t0 = time.time()
    presence = load_result("origin_presence_genome_wide")
    if not presence:
        raise SystemExit("origin_presence_genome_wide has not been computed")
    by_gene, species = masks(presence)
    pool = list(by_gene.values())
    wanted = set(args.libraries.split(",")) if args.libraries else None

    # the computed membership, not the curated gene list: core.replication has 201 members and 10
    # named genes, and profiling ten of them would be profiling the catalogue rather than the library
    kb = KnowledgeBase() if KnowledgeBase.available() else None
    rows = []
    for name, lib in sorted(LIBRARIES.items()):
        if wanted and name not in wanted:
            continue
        members = sorted(kb.members(lib)) if kb else sorted(lib.genes)
        row = profile_library(name, members, by_gene, species, pool)
        if row:
            rows.append(row)
            print(
                f"{name:34} {row['members_with_a_profile']:4} genes  "
                f"J {row['mean_pair_jaccard']:.4f} vs null {row['null_mean']:.4f}  "
                f"excess {row['excess_over_null']:+.4f}  | matched genes "
                f"{row['matched_gene_null_mean']} excess {row['excess_over_matched_genes']:+.4f} "
                f"({row['matched_draws_at_or_above_observed']} draws above)",
                flush=True,
            )

    rows.sort(key=lambda r: -(r["excess_over_null"] or 0))
    out = {
        "result": "profiling_genes",
        "species": species,
        "genes_with_a_profile": len(by_gene),
        "draws_per_library": DRAWS,
        "min_members": MIN_MEMBERS,
        "libraries": rows,
        "beaten_by_the_null": [r["library"] for r in rows if r["null_draws_at_or_above_observed"]],
        "null": (
            "each gene keeps the NUMBER of species it is present in and loses which ones. A null free "
            "to change prevalence would make two ubiquitous genes look like a discovery"
        ),
        "cannot_do": (
            "species are not phylogenetically independent, so a null that draws them independently "
            "overstates how surprising a co-occurrence is: every excess here is a ceiling on the "
            "evidence, not the evidence. A species-tree-aware null is what would settle a pair, and "
            "this project does not hold one. Nothing here calls a pair a dependency"
        ),
        "seconds": round(time.time() - t0, 1),
    }
    if args.no_save:
        print("\nnot saved (--no-save)")
    else:
        print(f"\nsaved {save_result(out['result'], out)}")
    beaten = out["beaten_by_the_null"]
    print(f"{len(rows)} libraries profiled; {len(beaten)} did not beat their own null: {beaten}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
