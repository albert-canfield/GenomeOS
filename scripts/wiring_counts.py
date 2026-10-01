"""List the wiring feasibility counts the result file keeps only as aggregates.

    uv run --frozen python scripts/wiring_counts.py

Read-only and descriptive: it loads the committed pairs and the frozen matching, then counts links,
alternatives, genes and locus clusters. It fits no model, computes no metric and writes no file. It exists
so the coordinator's record can be traced to a command. The counts are for the training pairs under the
primary rewiring (`element_kept`), under scheme S1 (distance bins 0.20 log10 wide) and S3 (0.30), both
without expression classes, as `wiring.SCHEMES` fixes them.

Three quantities it keeps apart, because they are not the same set:
  * movable: links with at least one admissible alternative;
  * moved: links a single valid rewiring actually moves, which is fewer, because each stratum must admit a
    derangement and the distance bins take a random offset per draw. Its mean, its smallest and its largest
    are observed over the registered draws; no maximum over all rewirings is proved anywhere;
  * the distinct genes, elements and chromosomes the movable links fall in, and the locus clusters they
    chain into. Chaining elements within 1 Mb is operational grouping; it establishes no biological
    independence, so nothing here is labelled independent.
"""

from __future__ import annotations

from collections import defaultdict

from genomeos.attribution import crispri, wiring

SCHEMES = (("S1", 0.20), ("S3", 0.30))
LOCUS_GAP = (
    1_000_000  # elements chained within this distance form one locus cluster: operational, not independence
)


def locus_clusters(elements: set[tuple[str, int, int]]) -> int:
    """Elements chained while consecutive starts lie within LOCUS_GAP of the running end.

    An operational grouping of nearby elements. It does not establish that the clusters are independent.
    """
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for chrom, start, end in elements:
        by_chrom[chrom].append((start, end))
    total = 0
    for spans in by_chrom.values():
        spans.sort()
        end: int | None = None
        for start, stop in spans:
            if end is None or start - end > LOCUS_GAP:
                total += 1
                end = stop
            else:
                end = max(end, stop)
    return total


def block(pair: crispri.Pair) -> tuple[float, ...]:
    return tuple(pair.features[column] for column in wiring.BLOCK)


def main() -> int:
    training = crispri.load(crispri.TRAINING)
    crispri.annotate(training, crispri.DeletionTable(), crispri.ElementCache())
    links = wiring.links_of(training, {})
    by_index = {link.index: link for link in links}
    groups, distinct = wiring.keys_for(links, "element_kept", expression_classes=False)
    regulated_genes: dict[tuple[str, int, int], set[str]] = defaultdict(set)
    for link in links:
        if link.regulated:
            regulated_genes[link.element].add(link.gene)
    regulated = [link for link in links if link.regulated]
    print(
        f"population links {len(links)}; regulated {len(regulated)}; "
        f"genes {len({link.gene for link in regulated})}; "
        f"elements {len({link.element for link in regulated})}; "
        f"locus clusters {locus_clusters({link.element for link in regulated})}"
    )

    for name, width in SCHEMES:
        alternatives = wiring.alternatives(links, groups, distinct, width)
        movable = [k for k, link in enumerate(links) if link.regulated and alternatives[k]]
        elements = {links[k].element for k in movable}
        moved_per_draw: list[int] = []
        unchanged_per_draw: list[int] = []
        ever_moved: set[int] = set()
        for draw in range(wiring.FEASIBILITY_DRAWS):
            rewiring = wiring.rewire(
                links, groups, distinct, width, wiring.rng_for("training", "element_kept", name, draw)
            )
            moved = [t for t in rewiring if by_index[t].regulated]
            ever_moved.update(moved)
            moved_per_draw.append(len(moved))
            unchanged_per_draw.append(
                sum(1 for t in moved if block(training[t]) == block(training[rewiring[t]]))
            )
        same_class = [k for k in movable if alternatives[k] <= regulated_genes[links[k].element]]
        assignments = [(k, gene) for k in same_class for gene in alternatives[k]]
        identical = 0
        for k, gene in assignments:
            link = links[k]
            other = next(x for x in links if x.element == link.element and x.gene == gene)
            identical += block(training[link.index]) == block(training[other.index])
        print(
            f"\n{name} (bins {width} log10, no expression classes)"
            f"\n  movable, at least one alternative: {len(movable)} links"
            f"\n  moved by one rewiring, over the {wiring.FEASIBILITY_DRAWS} registered draws: "
            f"mean {sum(moved_per_draw) / len(moved_per_draw):.1f}, smallest observed {min(moved_per_draw)}, "
            f"largest observed {max(moved_per_draw)} (observed, not a proved maximum)"
            f"\n  unique links moved across those draws: {len(ever_moved)}"
            f"\n  of the moved, attached feature block identical to their own (unchanged inputs): "
            f"mean {sum(unchanged_per_draw) / len(unchanged_per_draw):.1f} per draw"
            f"\n  distinct: genes {len({links[k].gene for k in movable})}, elements {len(elements)}, "
            f"chromosomes {len({links[k].chrom for k in movable})}; "
            f"locus clusters {locus_clusters(elements)} (1 Mb chaining, operational, not independence)"
            f"\n  every alternative a regulated link of the same element: {len(same_class)} links, "
            f"{len(assignments)} assignments, of which the attached feature block is identical in "
            f"{identical} (feature identity only; a differing block is not a measured score change)"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
