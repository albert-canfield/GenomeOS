# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Simpler genomes, area I: can ANY fly boundary set decide the between-gene question, and on what source.

Area I, after the fly's boundary node (genomeos/attribution/simpler_flynode.py,
docs/SIMPLER-GENOMES.md section 9). That run registered 917 between-gene elements as the question and
found the boundary rule changed the nearest-TSS gene for only **90** of them, under the stop rule's
200, so the between-gene question could not be decided. Its "what follows" named two routes: a
boundary set from the marker calls' own stages (6-16 h embryos), and a source with more between-gene
elements.

THIS MODULE IS THE INVENTORY AND THE STOP-RULE DECISION FOR THOSE TWO ROUTES. It reads positions
only -- the peak midpoints, the gene 5' ends and candidate boundary positions. No RNA value, no
outcome and no endpoint is read anywhere in it, and no request is sent to any model.

WHAT IT ADDS TO THE RECORD. Route (a) was treated in section 9 as a search problem: find a later
boundary set. It is also an arithmetic one, and the arithmetic can be settled without any new data.
Under the registered boundary rule a between-gene element changes gene only when its nearest 5' end
lies across a boundary AND the domain it leaves holds another 5' end. Add boundaries and the first
condition is met more often; add more and the domains shrink until they hold no other 5' end, the
rule falls back to nearest TSS (its EMPTY case) and the count falls again. So the number of
between-gene elements the rule can move is bounded by a curve in the boundary count, and
`density_sweep` traces it on the real element and gene positions. `tailored_ceiling` gives the other
end: what one boundary placed inside each element's own gap would reach, which no measured global set
can be expected to do.

Boundaries placed uniformly at random are, for this count, if anything GENEROUS: measured fly
boundaries sit mostly at active promoters (Ramirez et al. 2018, Nat Commun 9:189), and a boundary at
the nearest 5' end puts that 5' end inside the element's domain, which is the case where the rule
does NOT move. So a density the sweep cannot carry to 200 is not carried to 200 by a measured set of
that density either.

WHAT IT CANNOT ESTABLISH. Nothing here is a result about the fly's regulation, about any gene or
about any boundary. It is a count of how often two attribution rules name different genes, and a
licence-and-stage inventory of the sources that could supply either side. Accessibility is not
function; a boundary is not a regulatory unit by itself. Where a reading says a gene is attributed,
it says the attributed gene is expressed where the element is open, never that the element acts on
the gene.

THE STOP LIST. docs/ROADMAP.md section 5 reinstates "No more biosamples: 1.1 needs its readings to
mean what their names say." That clause is about the HUMAN READER's biosamples -- the 13 reader cells
of READER_TERMS. This module adds none and names none. A Drosophila Hi-C boundary list, and a
Drosophila accessibility or RNA source, are not reader biosamples and are not covered by that clause.
"""

from __future__ import annotations

import bisect
import csv
import math
import random
from pathlib import Path
from typing import Any

from genomeos.attribution import simpler
from genomeos.attribution import simpler_fly as sf
from genomeos.attribution import simpler_flynode as fn
from genomeos.attribution.simpler import Chromosome, Gene

#: The stop rule the record already set (docs/SIMPLER-GENOMES.md section 9): fewer than this many
#: between-gene elements changing gene under the boundary rule leaves the question undecidable.
STOP_RULE = fn.MIN_DISCORDANT

#: lane-flynode's registered counts on Stadler 2017's manual list, read first-hand from
#: data/results/simpler_flynode_registration.json (registration f5985f7).
FLYNODE = {
    "registration": "f5985f7",
    "boundaries_manual": 3121,
    "between_gene_elements": 917,
    "between_gene_elements_that_change_gene": 90,
    "primary_elements": 1634,
    "primary_elements_that_change_gene": 234,
    "median_domain_bp_around_a_primary_element": 84000,
    "outcome": "NO_DIFFERENCE_DETECTED",
}

#: Boundary counts the density sweep traces. The measured lists of section 9 hold 3,121 (manual),
#: 2,666 (computational) and 954 (representative), so the sweep brackets them by a factor of 30 up.
SWEEP_COUNTS = (2_000, 3_121, 4_500, 5_500, 7_000, 8_000, 12_000, 20_000, 35_000, 60_000, 100_000)
#: Draws per count. Each draw is one whole-genome placement, so the spread is the placement's.
SWEEP_DRAWS = 5
#: Fixed before any count was read.
SEED = 20261010

#: A tailored boundary sits midway in an element's gap; gaps under this are inside a boundary bin
#: either way, so no boundary is placed in them.
MIN_GAP_BP = 2 * fn.HALF_BIN + 100


# ------------------------------------------------------------------- the registered rule, at any bin size


def domain_of(bounds: list[int], pos: int, half_bin: int) -> tuple[float, float]:
    """`simpler_flynode.domain_of` with the source's own bin, not Stadler's 500 bp.

    Every candidate list calls boundaries in bins of its own resolution (500 bp, 800 bp, 2 kb), and
    the registered rule counts a 5' end inside either flanking boundary's bin as in the domain. With
    `half_bin == simpler_flynode.HALF_BIN` this is that function, which `tests/` checks row by row.
    """
    i = bisect.bisect_right(bounds, pos)
    lo = bounds[i - 1] - half_bin if i > 0 else -math.inf
    hi = bounds[i] + half_bin if i < len(bounds) else math.inf
    return lo, hi


def boundary_rule(chrom: Chromosome, bounds: list[int], pos: int, half_bin: int) -> tuple[Gene | None, str]:
    """`simpler_flynode.boundary_rule` with the source's own bin. The node itself is unchanged."""
    nearest = chrom.nearest_tss(pos)
    if nearest is None:
        return None, "none"
    lo, hi = domain_of(bounds, pos, half_bin)
    if lo <= nearest.tss <= hi:
        return nearest, fn.IN_DOMAIN
    a = bisect.bisect_left(chrom.tss, lo)
    b = bisect.bisect_right(chrom.tss, hi)
    cands = chrom.by_tss[a:b]
    if not cands:
        return nearest, fn.EMPTY
    return min(cands, key=lambda g: (abs(g.tss - pos), g.gid)), fn.MOVES


# --------------------------------------------------------------------------------- position-only counts


def between_gene_moves(
    rows: list[dict[str, Any]],
    chroms: dict[str, Chromosome],
    bounds: dict[str, list[int]],
    half_bin: int = fn.HALF_BIN,
) -> dict[str, int]:
    """How many between-gene elements change gene under the boundary rule, and how it splits by case.

    Reads the element midpoints, the gene 5' ends and `bounds` only. `rows` are
    `simpler.attribute` rows, so `domain_case` already says which elements lie between genes.
    """
    distal = [r for r in rows if r["kind"] == sf.DISTAL]
    moved_all = 0
    bg: list[tuple[str, bool]] = []
    for r in distal:
        gene, case = boundary_rule(chroms[r["chrom"]], bounds.get(r["chrom"], []), r["mid"], half_bin)
        if gene is None:
            raise ValueError(f"no gene on {r['chrom']}")
        moves = gene.gid != r["nearest_gene"]
        moved_all += int(moves)
        if r["domain_case"] != fn.HOST:
            bg.append((case, moves))
    return {
        "between_gene_elements": len(bg),
        "change_gene": sum(1 for _, m in bg if m),
        "nearest_5p_end_in_domain": sum(1 for c, _ in bg if c == fn.IN_DOMAIN),
        "empty_domain_falls_back_to_nearest_tss": sum(1 for c, _ in bg if c == fn.EMPTY),
        "all_distal_elements": len(distal),
        "all_distal_change_gene": moved_all,
    }


def read_boundary_table(
    rows: list[tuple[str, int, int]], arms: tuple[str, ...] = sf.CHROMOSOMES
) -> tuple[dict[str, list[int]], dict[str, Any]]:
    """(chrom, start, end) boundary spans -> the peaks' chromosome names and the spans' midpoints.

    The candidate lists name the arm without the `chr` the peak file carries, and give the bin as a
    half-open span; the position used is the span's midpoint and the half-bin half its width, so the
    source's own call is read and no threshold is chosen here.
    """
    out: dict[str, list[int]] = {}
    widths: dict[int, int] = {}
    skipped = 0
    for chrom, start, end in rows:
        name = chrom if chrom.startswith("chr") else "chr" + chrom
        if name[3:] not in arms:
            skipped += 1
            continue
        out.setdefault(name, []).append((start + end) // 2)
        widths[end - start] = widths.get(end - start, 0) + 1
    for name in out:
        out[name] = sorted(out[name])
    return out, {
        "boundaries": sum(len(v) for v in out.values()),
        "by_chromosome": {k: len(v) for k, v in sorted(out.items())},
        "span_widths_bp": dict(sorted(widths.items())),
        "rows_off_the_used_arms": skipped,
    }


def place_uniform(arm: dict[str, int], count: int, rnd: random.Random) -> dict[str, list[int]]:
    """`count` boundaries spread over the arms in proportion to their length, uniformly within each."""
    total = sum(arm.values())
    out: dict[str, list[int]] = {}
    for chrom, length in sorted(arm.items()):
        k = max(1, round(count * length / total))
        lo, hi = fn.HALF_BIN, length - fn.HALF_BIN
        out[chrom] = sorted(rnd.randrange(lo, hi) for _ in range(k))
    return out


def density_sweep(
    rows: list[dict[str, Any]],
    chroms: dict[str, Chromosome],
    arm: dict[str, int],
    counts: tuple[int, ...] = SWEEP_COUNTS,
    draws: int = SWEEP_DRAWS,
    seed: int = SEED,
) -> dict[str, Any]:
    """The curve: how many between-gene elements change gene, against the number of boundaries.

    Uniform placement, `draws` placements per count. Reports the mean, the minimum and the maximum
    over the draws, so the best draw at every density is in the record and not only the average.
    """
    rnd = random.Random(seed)
    total = sum(arm.values())
    rung = []
    for count in counts:
        moved, empty = [], []
        for _ in range(draws):
            got = between_gene_moves(rows, chroms, place_uniform(arm, count, rnd))
            moved.append(got["change_gene"])
            empty.append(got["empty_domain_falls_back_to_nearest_tss"])
        rung.append(
            {
                "boundaries": count,
                "mean_spacing_bp": total // count,
                "change_gene_mean": round(sum(moved) / len(moved), 1),
                "change_gene_min": min(moved),
                "change_gene_max": max(moved),
                "empty_domain_mean": round(sum(empty) / len(empty), 1),
            }
        )
    best = max(rung, key=lambda r: r["change_gene_max"])
    return {
        "placement": "uniform within each arm, arms weighted by length",
        "draws_per_count": draws,
        "seed": seed,
        "rungs": rung,
        "best_rung": best,
        "best_single_draw_change_gene": best["change_gene_max"],
        "reaches_stop_rule_at_any_density": best["change_gene_max"] >= STOP_RULE,
    }


def tailored_ceiling(rows: list[dict[str, Any]], chroms: dict[str, Chromosome]) -> dict[str, Any]:
    """What a boundary set placed inside each element's own gap would reach: the rule's own ceiling.

    One boundary midway between each between-gene element and its nearest 5' end. No measured set is
    placed this way, and none is expected to be; the number says what the rule could do at best, so
    that the sweep's plateau is read as a property of measured densities and not of the rule.
    """
    bounds: dict[str, set[int]] = {}
    for r in rows:
        if r["kind"] != sf.DISTAL or r["domain_case"] == fn.HOST:
            continue
        gene = chroms[r["chrom"]].nearest_tss(r["mid"])
        if gene is None or abs(gene.tss - r["mid"]) < MIN_GAP_BP:
            continue
        bounds.setdefault(r["chrom"], set()).add((r["mid"] + gene.tss) // 2)
    placed = {c: sorted(v) for c, v in bounds.items()}
    got = between_gene_moves(rows, chroms, placed)
    return {
        "boundaries_placed": sum(len(v) for v in placed.values()),
        "min_gap_bp": MIN_GAP_BP,
        **got,
    }


# --------------------------------------------------------------------------------- the stop-rule verdict


def mde_points(discordant: int, population: int) -> float | None:
    """The 80%-power, two-sided 5% difference in points before clustering, as simpler.census bounds it."""
    return round(100 * 2.8 * math.sqrt(discordant) / population, 3) if population else None


def verdict(candidates: list[dict[str, Any]]) -> dict[str, Any]:
    """NO_SOURCE_QUALIFIES, or the candidate whose counted between-gene moves reach the stop rule.

    A candidate qualifies only on a COUNT, never on a count it was not possible to make: a candidate
    whose boundary list or element set could not be obtained carries `counted: False` and cannot
    qualify, and the reason it could not be counted is kept with it.
    """
    counted = [c for c in candidates if c.get("counted")]
    passing = [c for c in counted if (c.get("between_gene_elements_that_change_gene") or 0) >= STOP_RULE]
    best = max(
        counted,
        key=lambda c: c.get("between_gene_elements_that_change_gene") or 0,
        default=None,
    )
    return {
        "stop_rule": STOP_RULE,
        "stop_rule_is": (
            "at least this many between-gene elements must change gene under the boundary rule for "
            "the test to be decidable; it is the record's, not this lane's (section 9)"
        ),
        "candidates_counted": len(counted),
        "candidates_not_counted": len(candidates) - len(counted),
        "qualifying": [c["name"] for c in passing],
        "best_counted": None
        if best is None
        else {
            "name": best["name"],
            "between_gene_elements": best.get("between_gene_elements"),
            "between_gene_elements_that_change_gene": best.get("between_gene_elements_that_change_gene"),
        },
        "outcome": "SOURCE_QUALIFIES" if passing else "NO_SOURCE_QUALIFIES",
    }


# ------------------------------------------------------------------------- the source that qualified
#
# Route (a) is closed by the counts above: no measured boundary set, at any density, carries the
# between-gene count on Calderon 2022's 917 elements to the stop rule's 200. Route (b) is open --
# a source with more between-gene elements. The one that reaches the stop rule on a count is the fly
# brain: Janssens et al. 2022's single-cell accessibility, whose own cell-type table declares, per
# cell type, which cluster of Davie et al. 2018's brain RNA atlas is its counterpart. That declared
# map is what makes one label vocabulary cover both modalities, which is lane-fly's design.

STORE = Path("data/knowledge/simpler_flysource")
DAR_FILE = STORE / "tbl_DARs.txt"
CELL_TYPES_FILE = STORE / "cellTypeInfo.tsv"
RNA_FILE = STORE / "Davie_Janssens_Koldere_et_al_2018_AdultBrain.loom"
DOLSTEN_FILE = STORE / "dolsten2025_TableS2.xlsx"
KAUSHAL_FILE = STORE / "kaushal2022_data_file_s2.xlsx"
KAUSHAL_CNS_FILE = STORE / "kaushal2021_larval_cns_data_s3.xlsx"
FETCH_RECORD = STORE / "fetch_record.json"

ELEMENT_SOURCE = (
    "Janssens J, Aibar S, Taskiran II, et al. (2022) Decoding gene regulation in the fly brain. "
    "Nature 601:630-636, doi:10.1038/s41586-021-04262-z; the authors' differentially accessible "
    "region table (tbl_DARs.txt) and cell-type table (cellTypeInfo.tsv) from flybrain.aertslab.org; "
    "GEO GSE163697 holds the reads"
)
ELEMENT_LICENCE = (
    "not stated: the article is subscription-only (Europe PMC records no licence and no PMCID) and "
    "the data portal flybrain.aertslab.org states no licence, no terms of use and no copyright line; "
    "the GEO record states no terms either. Recorded, not inferred. This is the same standing as "
    "Calderon et al. 2022, which lane-fly used and recorded the same way (section 8)"
)
RNA_SOURCE = (
    "Davie K, Janssens J, Koldere D, et al. (2018) A single-cell transcriptome atlas of the ageing "
    "Drosophila brain. Cell 174:982-998; the loom served by flybrain.aertslab.org as the RNA atlas "
    "whose clusters Janssens et al. 2022's cellTypeInfo.tsv names as the counterpart of each "
    "accessibility cell type"
)
RNA_LICENCE = "not stated by the portal that serves the loom; recorded, not inferred"

#: Which rows of the DAR table are read. The source's own call, at its own resolution and contrast;
#: no threshold, caller or cut-off is chosen here.
DAR_SELECTION = {
    "clResolution": "lvl1",
    "DAR_contrast": "vsRest",
    "DAR_within": "WholeBrain",
    "DAR_type": "Up",
}
#: Columns of the DAR table that carry the source's own region-to-gene link. NEVER read: the lane
#: tests attribution rules, so a source's own answer would be the answer.
NEVER_READ = ("closestGene",)

#: The one vocabulary for both modalities: the RNA clusters that Janssens et al. 2022's
#: cellTypeInfo.tsv names in its `RNAcluster` column, which is the source's own declared map from an
#: accessibility cell type to a cluster of the RNA atlas. Two accessibility cell types (T4 and T5)
#: name one RNA cluster (T4/T5), so the label is the RNA cluster and both map onto it; a region
#: accessible in T4 and in T5 and in nothing else therefore carries one label, not two.
LABELS = (
    "15a",
    "15b",
    "26",
    "90",
    "97",
    "98",
    "A/B*-KC",
    "A/B-KC",
    "AST",
    "C2",
    "C3",
    "CSM",
    "CTX",
    "Dm3",
    "Dm8/Dm11",
    "ENS",
    "G-KC",
    "L2/L3",
    "LLPC1",
    "Mi1",
    "Mi15",
    "Mi9",
    "PB4/5",
    "PLM",
    "PNG",
    "PR",
    "Pm4",
    "Poxn",
    "SUB",
    "T1",
    "T2",
    "T2a",
    "T3",
    "T4/T5",
    "Tm1/TmY8",
    "Tm2",
    "Tm20",
    "Tm3",
    "Tm4",
    "Tm9",
    "TmY3",
    "TmY5a",
)
#: Accessibility cell types with no declared RNA counterpart are excluded on both sides.
ACCESSIBILITY_CELL_TYPES_WITH_A_DECLARED_RNA_CLUSTER = 43

# simpler_fly's thresholds and statistics, carried unchanged.
DETECTION_CPM = sf.DETECTION_CPM
DISTANCE_BINS = sf.DISTANCE_BINS
CONTROL_DRAWS = sf.CONTROL_DRAWS
BOOTSTRAP = sf.BOOTSTRAP
CLUSTER_BP = sf.CLUSTER_BP
ALPHA = sf.ALPHA
ROTATION_DRAWS = fn.ROTATION_DRAWS
#: Labels are shown beside when they hold at least this many primary elements.
MIN_LABEL = fn.MIN_TISSUE

#: The primary boundary set is lane-flynode's, carried unchanged: only the element source changes, so
#: that a disagreement with section 9 is attributable to the source and not to the boundaries too
#: (simpler_fly.py's own reason, section 8 item 1). The other two measured sets are beside.
BOUNDARY_SETS = {
    "stadler_manual_nc14": {"half_bin": fn.HALF_BIN, "primary": True},
    "dolsten_2025_atlas": {"half_bin": 400, "primary": False},
    "kaushal_2022_wt_embryo": {"half_bin": 1000, "primary": False},
    "kaushal_2021_wt_larval_cns": {"half_bin": 1000, "primary": False},
}
PRIMARY_SET = "stadler_manual_nc14"


# ------------------------------------------------------------------------------------------- reading


def read_label_map(path: Path = CELL_TYPES_FILE) -> dict[str, str]:
    """The source's own accessibility cell type -> RNA cluster map, from its `RNAcluster` column."""
    out: dict[str, str] = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rna = (row.get("RNAcluster") or "").strip()
            if rna and rna != "-":
                out[row["CellType_lvl1"].strip()] = rna
        if sorted(set(out.values())) != sorted(LABELS):
            raise ValueError("the source's declared RNA clusters are not the frozen label set")
    return out


def read_dars(path: Path = DAR_FILE, label_of: dict[str, str] | None = None) -> dict[str, set[str]]:
    """Region -> the set of labels the source calls it differentially accessible UP in.

    Only the rows of `DAR_SELECTION` are read, and only the region, the cell type and the selection
    columns. `closestGene` is never read.
    """
    label_of = label_of if label_of is not None else read_label_map()
    out: dict[str, set[str]] = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for key in NEVER_READ:
            if key not in (reader.fieldnames or []):
                raise ValueError(f"the DAR table no longer has the never-read column {key}")
        for row in reader:
            if any(row[k] != v for k, v in DAR_SELECTION.items()):
                continue
            label = label_of.get(row["cellType"])
            if label is None:
                continue
            out.setdefault(row["region"], set()).add(label)
    return out


def elements(
    dars: dict[str, set[str]], chroms: dict[str, Chromosome], tss: dict[str, list[int]]
) -> tuple[list[simpler.Element], dict[str, int]]:
    """The population: regions the source calls up in EXACTLY one label, classified as simpler_fly does."""
    out: list[simpler.Element] = []
    counts = {
        "regions_up_in_at_least_one_label": len(dars),
        "more_than_one_label": 0,
        "chromosome_not_used": 0,
        "neither_class": 0,
        sf.DISTAL: 0,
        sf.PROMOTER: 0,
    }
    for region, labels in sorted(dars.items()):
        if len(labels) != 1:
            counts["more_than_one_label"] += 1
            continue
        label = next(iter(labels))
        chrom, span = region.rsplit(":", 1)
        start, end = (int(x) for x in span.split("-"))
        name = chrom if chrom.startswith("chr") else "chr" + chrom
        if name.removeprefix("chr") not in sf.CHROMOSOMES or name not in chroms:
            counts["chromosome_not_used"] += 1
            continue
        probe = simpler.Element(name, start, end, "", label)
        kind = sf.classify(name, probe.mid, chroms, tss)
        if kind is None:
            counts["neither_class"] += 1
            continue
        counts[kind] += 1
        out.append(simpler.Element(name, start, end, kind, label))
    return out, counts


def read_expression(
    path: Path = RNA_FILE, name_to_id: dict[str, str] | None = None
) -> tuple[dict[str, tuple[float, ...]], dict[str, Any]]:
    """The endpoint: per protein-coding gene (FBgn), CPM in each of `LABELS`. Read only by `score`.

    The loom's counts are summed over the cells of each label, exactly as lane-fly pooled Calderon's
    UMIs over the cells of each tissue, and turned into counts per million of that label's own total.
    Cells whose `cell_type` is not one of `LABELS` are not read into any label's sum; they are counted
    so that the share of the atlas the labels cover is in the record.
    """
    import h5py  # a loom is HDF5; imported here so the module loads without it

    name_to_id = name_to_id or {}
    with h5py.File(path, "r") as f:
        genes = [_text(x) for x in f["row_attrs/Gene"][:]]
        cell_types = [_text(x) for x in f["col_attrs/cell_type"][:]]
        index = {lab: i for i, lab in enumerate(LABELS)}
        column = [index.get(c, -1) for c in cell_types]
        cells_per_label = [0] * len(LABELS)
        for j in column:
            if j >= 0:
                cells_per_label[j] += 1
        matrix = f["matrix"]
        sums: dict[str, list[float]] = {}
        totals = [0.0] * len(LABELS)
        for i, gene in enumerate(genes):
            row = matrix[i, :]
            v = [0.0] * len(LABELS)
            for j, x in zip(column, row, strict=True):
                if j >= 0:
                    v[j] += float(x)
            for j in range(len(LABELS)):
                totals[j] += v[j]
            gid = name_to_id.get(gene)
            if gid is not None:
                sums[gid] = v
    expr = {
        g: tuple(1e6 * v[j] / totals[j] if totals[j] else 0.0 for j in range(len(LABELS)))
        for g, v in sums.items()
    }
    return expr, {
        "rows_in_the_rna_atlas": len(genes),
        "protein_coding_genes_matched_to_an_fbgn": len(expr),
        "cells_in_the_rna_atlas": len(cell_types),
        "cells_under_a_label": sum(cells_per_label),
        "cells_per_label": dict(zip(LABELS, cells_per_label, strict=True)),
        "label_count_totals": {k: round(v, 1) for k, v in zip(LABELS, totals, strict=True)},
    }


def _text(value: Any) -> str:
    return value.decode() if isinstance(value, bytes) else str(value)


def read_dolsten(path: Path = DOLSTEN_FILE) -> list[tuple[str, int, int]]:
    """Dolsten et al. 2025 Table S2, sheet `boundaries`: the atlas's boundary spans, as given."""
    return [
        (_text(r[1]), int(r[2]), int(r[3]))
        for r in _rows_of(path, "boundaries")
        if r[1] is not None and str(r[1]) != "chrom"
    ]


def read_kaushal(
    path: Path = KAUSHAL_FILE, exist_column: str = "exist.WT_0", half_bin: int = 1000
) -> list[tuple[str, int, int]]:
    """A Gambetta-lab TopDom table: the 2-kb bins the authors mark as present in wild type.

    Kaushal et al. 2022 (2-6 h embryos) names that column `exist.WT_0`; Kaushal et al. 2021 (larval
    central nervous system) names it `exist.wildtype`. Only the authors' own call is read.
    """
    rows = list(_rows_of(path, None))
    header = [_text(x) for x in rows[0]]
    exists = header.index(exist_column)
    out = []
    for r in rows[1:]:
        if r[0] is None or _text(r[exists]).strip().lower() != "true":
            continue
        centre = int(r[1])
        out.append((_text(r[0]), centre - half_bin, centre + half_bin))
    return out


def _rows_of(path: Path, sheet: str | None) -> list[tuple]:
    import openpyxl

    book = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = book[sheet] if sheet else book.worksheets[0]
    return list(ws.iter_rows(values_only=True))


# -------------------------------------------------------------------------- attribution and the census


def attribute(
    rows: list[dict[str, Any]],
    chroms: dict[str, Chromosome],
    bounds: dict[str, list[int]],
    half_bin: int,
) -> list[dict[str, Any]]:
    """`simpler_flynode.attribute` with the set's own bin: the boundary gene, its case and distance."""
    out = []
    for r in rows:
        ch = chroms[r["chrom"]]
        gene, case = boundary_rule(ch, bounds.get(r["chrom"], []), r["mid"], half_bin)
        if gene is None:
            raise ValueError(f"no gene on {r['chrom']}")
        hosts = {h.gid for h in ch.hosts(r["mid"])}
        out.append(
            {
                **r,
                "boundary_gene": gene.gid,
                "boundary_case": case,
                "boundary_dist": abs(gene.tss - r["mid"]),
                "boundary_gene_contains_element": gene.gid in hosts,
                "nearest_gene_contains_element": r["nearest_gene"] in hosts,
            }
        )
    return out


def distal(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [r for r in rows if r["kind"] == sf.DISTAL]


def between_genes(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """THE POPULATION: distal elements whose midpoint lies in no protein-coding gene body.

    This is the brief's population and section 9's open question, kept as the primary one here --
    unlike section 9, which had to move to a wider population because this one held only 90 elements
    on which the rules differ. Here it holds 216, over the stop rule, so no move is needed.
    """
    return [r for r in distal(rows) if r["domain_case"] != fn.HOST]


def _population(rows: list[dict[str, Any]]) -> dict[str, Any]:
    n = len(rows)
    moves = [r for r in rows if r["boundary_gene"] != r["nearest_gene"]]
    k_gs = sum(1 for r in rows if r["boundary_gene"] != r["domain_gene"])
    return {
        "elements": n,
        "elements_by_label": {
            lab: sum(1 for r in rows if r["tissue"] == lab)
            for lab in LABELS
            if any(r["tissue"] == lab for r in rows)
        },
        "boundary_cases": {c: sum(1 for r in rows if r["boundary_case"] == c) for c in fn.CASES},
        "differ_from_nearest_tss": len(moves),
        "differ_from_nearest_tss_where_the_nearest_gene_contains_the_element": sum(
            1 for r in moves if r["nearest_gene_contains_element"]
        ),
        "differ_from_gene_structure": k_gs,
        "gene_structure_names_a_host": sum(1 for r in rows if r["domain_case"] == fn.HOST),
        "clusters_1mb": len({(r["chrom"], r["mid"] // CLUSTER_BP) for r in rows}),
        "mde_points_unclustered_vs_nearest_tss": mde_points(len(moves), n),
        "mde_points_unclustered_vs_gene_structure": mde_points(k_gs, n),
    }


def census(
    main: list[dict[str, Any]],
    alternatives: dict[str, list[dict[str, Any]]],
    funnel: dict[str, int],
) -> dict[str, Any]:
    """Counts and power in advance, from positions, the element-side label and the boundaries only."""
    out = {
        "population_funnel": funnel,
        "between_genes": _population(between_genes(main)),
        "all_distal": _population(distal(main)),
        "alternative_boundary_sets_between_genes": {
            name: _population(between_genes(rows)) for name, rows in sorted(alternatives.items())
        },
        "positive_control_promoters": sum(1 for r in main if r["kind"] == sf.PROMOTER),
        "labels": len(LABELS),
    }
    out["stop_rule_met"] = out["between_genes"]["differ_from_nearest_tss"] < STOP_RULE
    return out


# ------------------------------------------------------------------------------------------ statistics


def concordant(cpm: tuple[float, ...] | None, label: str) -> bool:
    """At least DETECTION_CPM in the element's own label and higher there than in every other."""
    if cpm is None:
        return False
    i = LABELS.index(label)
    v = cpm[i]
    return v >= DETECTION_CPM and all(v > cpm[j] for j in range(len(LABELS)) if j != i)


class ControlIndex:
    """`simpler_fly.ControlIndex` over `LABELS`: the same matched-distance draw, exactly.

    A matched-control gene is drawn uniformly from the protein-coding genes of the element's
    chromosome whose 5' end lies in the same distance bin as the attributed gene's, the attributed
    gene excluded; the chance it is concordant is computed exactly from prefix counts.
    """

    def __init__(self, chroms: dict[str, Chromosome], expr: dict[str, tuple[float, ...]]):
        self.chroms = chroms
        self.expr = expr
        self.prefix: dict[tuple[str, str], list[int]] = {}
        for name, ch in chroms.items():
            for lab in LABELS:
                acc, p = 0, [0]
                for g in ch.by_tss:
                    acc += concordant(expr.get(g.gid), lab)
                    p.append(acc)
                self.prefix[(name, lab)] = p

    def _count(self, chrom: str, label: str, lo_pos: float, hi_pos: float) -> tuple[int, int]:
        tss = self.chroms[chrom].tss
        i, j = bisect.bisect_left(tss, lo_pos), bisect.bisect_left(tss, hi_pos)
        p = self.prefix[(chrom, label)]
        return j - i, p[j] - p[i]

    def probability(self, chrom: str, label: str, mid: int, gene: str, dist: int) -> float:
        b = simpler.distance_bin(dist)
        lo = 0 if b == 0 else DISTANCE_BINS[b - 1]
        hi = DISTANCE_BINS[b] if b < len(DISTANCE_BINS) else math.inf
        n1, c1 = self._count(chrom, label, mid - hi + 1, mid - lo + 1)
        n2, c2 = self._count(chrom, label, mid + lo, mid + hi)
        if lo == 0:
            n0, c0 = self._count(chrom, label, mid, mid + 1)
            n1, c1 = n1 - n0, c1 - c0
        total, conc = n1 + n2 - 1, c1 + c2 - int(concordant(self.expr.get(gene), label))
        if total > 0:
            return conc / total
        ch = self.chroms[chrom]
        return sum(concordant(self.expr.get(g.gid), label) for g in ch.genes) / max(len(ch.genes), 1)


def compare(rows: list[dict[str, Any]], expr: dict[str, tuple[float, ...]], a: str, b: str) -> dict[str, Any]:
    """`simpler_flynode.compare`, unchanged but for the label vocabulary."""
    ca = [int(concordant(expr.get(r[a]), r["tissue"])) for r in rows]
    cb = [int(concordant(expr.get(r[b]), r["tissue"])) for r in rows]
    diffs = [{"chrom": r["chrom"], "mid": r["mid"], "d": x - y} for r, x, y in zip(rows, ca, cb, strict=True)]
    bs = simpler.cluster_bootstrap(diffs, "d") if rows else None
    n = len(rows)
    return {
        "n": n,
        "a": a,
        "b": b,
        "differ": sum(1 for r in rows if r[a] != r[b]),
        "a_concordant": sum(ca),
        "b_concordant": sum(cb),
        "a_rate_pct": simpler._rate(sum(ca), n),
        "b_rate_pct": simpler._rate(sum(cb), n),
        "delta_points": round(100 * bs["point"], 3) if bs else None,
        "delta_ci95_points": [round(100 * bs["lo"], 3), round(100 * bs["hi"], 3)] if bs else None,
        "clusters": bs["clusters"] if bs else 0,
    }


def rotation_null(
    rows: list[dict[str, Any]],
    chroms: dict[str, Chromosome],
    bounds: dict[str, list[int]],
    arm: dict[str, int],
    expr: dict[str, tuple[float, ...]],
    observed: int,
    half_bin: int,
    draws: int = ROTATION_DRAWS,
    seed: int = SEED,
) -> dict[str, Any]:
    """`simpler_flynode.rotation_null`: the same boundaries rotated along each arm by one offset."""
    rng = random.Random(seed)
    sims = []
    for _ in range(draws):
        rotated = {}
        for name in sorted(bounds):
            off = rng.randrange(arm[name])
            rotated[name] = sorted((b + off) % arm[name] for b in bounds[name])
        k = 0
        for r in rows:
            gene, _ = boundary_rule(chroms[r["chrom"]], rotated.get(r["chrom"], []), r["mid"], half_bin)
            k += concordant(expr.get(gene.gid) if gene else None, r["tissue"])
        sims.append(k)
    ge = sum(1 for s in sims if s >= observed)
    n = len(rows)
    return {
        "draws": draws,
        "observed": observed,
        "mean_rotated": round(sum(sims) / draws, 3),
        "mean_rotated_rate_pct": round(100 * sum(sims) / draws / n, 3) if n else None,
        "draws_at_or_above_observed": ge,
        "p_one_sided": round((1 + ge) / (1 + draws), 5),
    }


def positive_control(
    rows: list[dict[str, Any]], idx: ControlIndex, expr: dict[str, tuple[float, ...]]
) -> dict[str, Any]:
    """simpler_fly's positive control: nearest TSS on the promoter elements, against its own control."""
    prom = [r for r in rows if r["kind"] == sf.PROMOTER]
    k = sum(int(concordant(expr.get(r["nearest_gene"]), r["tissue"])) for r in prom)
    probs = [
        idx.probability(r["chrom"], r["tissue"], r["mid"], r["nearest_gene"], r["nearest_dist"]) for r in prom
    ]
    out: dict[str, Any] = {
        "n": len(prom),
        "nearest_concordant": k,
        "nearest_rate_pct": simpler._rate(k, len(prom)),
    }
    out["control"] = simpler.control_test(probs, k)
    out["responsive"] = out["control"]["p_one_sided"] < ALPHA
    return out


# ---------------------------------------------------------------------------------- the registration


def _pts(x: float) -> str:
    return f"{x:+.2f}"


def _fill(res: dict[str, Any]) -> dict[str, str]:
    p = res["primary"]
    g = res["primary_vs_gene_structure"]
    s = res["all_distal"]
    sg = res["all_distal_vs_gene_structure"]
    ctl = res["control"]
    rot = res["rotation_beside"]
    cen = res["census"]["between_genes"]
    return {
        "n": f"{p['n']:,}",
        "kb": f"{p['a_concordant']:,}",
        "kn": f"{p['b_concordant']:,}",
        "cb": f"{p['a_rate_pct']:.1f}",
        "cn": f"{p['b_rate_pct']:.1f}",
        "k": f"{p['differ']:,}",
        "d": _pts(p["delta_points"]),
        "lo": _pts(p["delta_ci95_points"][0]),
        "hi": _pts(p["delta_ci95_points"][1]),
        "mde": f"{cen['mde_points_unclustered_vs_nearest_tss']:.1f}",
        "kg": f"{g['b_concordant']:,}",
        "cg": f"{g['b_rate_pct']:.1f}",
        "k2": f"{g['differ']:,}",
        "d2": _pts(g["delta_points"]),
        "lo2": _pts(g["delta_ci95_points"][0]),
        "hi2": _pts(g["delta_ci95_points"][1]),
        "N": f"{s['n']:,}",
        "kA": f"{s['differ']:,}",
        "kbA": f"{s['a_concordant']:,}",
        "cbA": f"{s['a_rate_pct']:.1f}",
        "knA": f"{s['b_concordant']:,}",
        "cnA": f"{s['b_rate_pct']:.1f}",
        "dA": _pts(s["delta_points"]),
        "loA": _pts(s["delta_ci95_points"][0]),
        "hiA": _pts(s["delta_ci95_points"][1]),
        "kgA": f"{sg['b_concordant']:,}",
        "cgA": f"{sg['b_rate_pct']:.1f}",
        "dgA": _pts(sg["delta_points"]),
        "logA": _pts(sg["delta_ci95_points"][0]),
        "higA": _pts(sg["delta_ci95_points"][1]),
        "cc": f"{100 * ctl['expected_rate']:.1f}",
        "p": simpler._fmt_p(ctl["p_one_sided"]),
        "rm": f"{rot['mean_rotated_rate_pct']:.1f}",
        "pr": simpler._fmt_p(rot["p_one_sided"]),
        "t": f"{len(LABELS)}",
        "nb": f"{res['boundaries']:,}",
        "miss": f"{res['attributed_genes_missing_from_rna_rows']:,}",
        "fly": f"{FLYNODE['between_gene_elements_that_change_gene']:,}",
    }


def reading(res: dict[str, Any]) -> str:
    """The registered outcome wording, filled in. The templates are REGISTRATION['outcome_wording']."""
    w = REGISTRATION["outcome_wording"]
    fill = _fill(res)
    parts = [
        w["head"],
        w[res["comparison"]],
        w["meaning_" + res["comparison"]],
        w["gene_structure"],
        w["gene_structure_" + res["comparison_vs_gene_structure"]],
        w["FALSIFIED" if res["falsifier"] == "FALSIFIED" else "NOT_FALSIFIED"],
        w["rotation_" + ("ahead" if res["rotation_beside"]["p_one_sided"] < ALPHA else "not_ahead")],
        w["all_distal"],
        w["missing"],
        w["tail"],
    ]
    return " ".join(x.format(**fill) for x in parts)


#: Written and committed before `score` read any RNA value. Not to be changed afterwards.
REGISTRATION: dict[str, Any] = {
    "question": (
        "Between genes, where no rule can name the gene whose body the element sits in, does the gene "
        "with the nearest 5' end inside the element's boundary-delimited domain match the element's "
        "own cell type more often than the gene with the nearest 5' end on the arm? Section 9 could "
        "not ask it: on Calderon et al. 2022's 917 between-gene elements the two rules named the same "
        "gene for all but 90, under the stop rule's 200. This lane's counts show that no measured "
        "boundary set reaches 200 on that population at any density, so the source changes, not the "
        "boundary set."
    ),
    "route": (
        "route (b) of section 9's 'what follows': a source with more between-gene elements. Route (a), "
        "a 6-16 h boundary set, is closed by the counts in this registration's inventory."
    ),
    "element_source": ELEMENT_SOURCE,
    "element_licence": ELEMENT_LICENCE,
    "rna_source": RNA_SOURCE,
    "rna_licence": RNA_LICENCE,
    "boundary_set": {
        "primary": (
            "Stadler, Haines, Eisen 2017 (eLife 6:e29550, CC BY 4.0), supplementary file 1, Table S4: "
            "the authors' manual call on nuclear-cycle-14 embryos at 500-bp resolution, dm6 as given. "
            "lane-flynode's set, carried unchanged (registration f5985f7)"
        ),
        "why_unchanged": (
            "only the element source changes here, so a disagreement with section 9 is attributable to "
            "the source and not to the boundaries as well; this is simpler_fly.py's own reason "
            "(section 8 item 1) applied the other way round"
        ),
        "beside_deciding_nothing": [
            "Dolsten et al. 2025 (Cell Genomics 5:101002, CC BY-NC-ND 4.0), Table S2: 3,521 boundary "
            "spans of a Micro-C atlas covering nc1-8, nc14 and stages 10-12 (about 6-8 h), 800-bp "
            "spans at 400-bp insulation bins, Release 6. Its text states 3,673 boundaries and the "
            "shipped sheet holds 3,521 rows; the sheet is used as shipped and the difference is noted",
            "Kaushal et al. 2022 (Sci Adv 8:eabl8834, CC BY 4.0), data file S2: the 2,334 2-kb bins "
            "the authors mark `exist.WT_0` in wild-type 2-6 h embryos, dm6",
            "Kaushal et al. 2021 (Nat Commun 12:1011, CC BY 4.0), supplementary data 3: the 3,458 "
            "2-kb bins the authors mark `exist.wildtype` in wild-type third-instar larval central "
            "nervous system, dm6. This is the only set whose TISSUE matches the elements', and it is "
            "beside rather than primary only because the primary set is carried from f5985f7 "
            "unchanged; it is the obvious primary of a next registration",
        ],
        "stage": (
            "the primary boundaries come from 2-3 h embryos while the elements and the RNA come from "
            "larval-to-adult brain. That is a WIDER stage gap than section 9's, not a narrower one, "
            "and it is the price of a population large enough to ask the question. It is stated in the "
            "reading, it is not tested here, and the larval-CNS set reported beside is the one whose "
            "tissue matches"
        ),
    },
    "node": (
        "simpler_flynode's boundary rule, unchanged but for each set's own bin: an element's domain "
        "runs from the boundary at or left of its midpoint to the next boundary right of it, or to the "
        "arm end; a gene is in the domain when its 5' end lies there or inside either boundary's bin; "
        "the rule names the in-domain gene with the nearest 5' end, and falls back to nearest TSS when "
        "the domain holds no 5' end"
    ),
    "population": (
        "the between-gene elements: regions Janssens et al. 2022 calls differentially accessible UP, "
        "at cluster resolution lvl1 against the rest of the whole brain, in EXACTLY one of the 42 RNA "
        "clusters its own cellTypeInfo.tsv declares a counterpart for, whose midpoint lies more than "
        "1 kb from any transcript 5' end and in no protein-coding gene body"
    ),
    "labels": LABELS,
    "endpoint": (
        "simpler_fly's, unchanged: the attributed gene is counted concordant when it reaches "
        "5 counts per million in the element's own label and is higher there than in every other "
        "label, in RNA summed over that label's cells. A gene absent from the RNA atlas's rows is "
        "not concordant, and the number of such genes is reported"
    ),
    "outcome": (
        "the boundary rule against nearest TSS on the between-gene elements; percentile cluster "
        "bootstrap over 1-Mb windows; BOUNDARY_RULE_AHEAD, NEAREST_AHEAD or NO_DIFFERENCE_DETECTED"
    ),
    "second_registered_comparison": (
        "the boundary rule against the worm's gene-structure rule on the same elements, with its own "
        "label, not the outcome. Two comparisons on one population, each at 95%, no correction"
    ),
    "matched_control_and_falsifier": (
        "simpler_fly's matched-distance control applied to the boundary rule's gene on the population: "
        "FALSIFIED when the rule's genes match the element's label no more often than genes drawn at a "
        "matched distance (one-sided p >= 0.05)"
    ),
    "positive_control": (
        "simpler_fly's: nearest TSS on the promoter elements (within 500 bp of a 5' end) against its "
        "own matched control. If it does not separate, the readout is declared unresponsive and the "
        "primary reading is withheld"
    ),
    "beside_deciding_nothing": [
        "the same boundaries rotated along each arm by one uniform offset, number and spacing kept, "
        "1,000 draws",
        "the two other measured boundary sets on the same elements",
        "the elements where the two rules differ, alone",
        "all distal elements, and the gene-structure rule against nearest TSS",
        "labels holding 50 or more population elements",
    ],
    "power_in_advance": (
        "216 of 3,265 population elements differ under the primary boundary set: about 1.3 points at "
        "80% power, two-sided 5%, before clustering, against the 2.9 points section 9 could bound on "
        "its 90"
    ),
    "stop_rule": (
        "one run; fewer than 200 differing between-gene elements stops it unscored. The census counts "
        "216 under the primary set and, beside, 209 under Dolsten 2025, 174 under Kaushal 2021's "
        "larval CNS and 162 under Kaushal 2022. Only the primary set and Dolsten 2025 reach 200, so "
        "the two sets reported beside that do not are read as beside and nothing else"
    ),
    "money": "none; no model request, no key read",
    "outcome_wording": {
        "head": (
            "On {n} D. melanogaster brain distal elements called differentially accessible in exactly "
            "one of {t} cell types (Janssens et al. 2022, larval to adult brain, the source's own "
            "lvl1 call against the rest of the brain) and lying in no protein-coding gene body, the "
            "gene with the nearest 5' end inside the element's domain, delimited by {nb} Hi-C "
            "boundaries of nuclear-cycle-14 embryos (Stadler et al. 2017, the authors' manual call), "
            "is expressed at 5 counts per million or more and highest in the element's own cell type, "
            "in RNA summed over that cell type's cells (Davie et al. 2018, the counterpart cluster "
            "the accessibility source itself declares), for {kb} of {n} ({cb}%), against {kn} of {n} "
            "({cn}%) for the gene with the nearest 5' end on the arm; the rules differ on {k} of {n} "
            "elements; difference {d} points, 95% {lo} to {hi} over 1-Mb clusters."
        ),
        "BOUNDARY_RULE_AHEAD": "On this endpoint the boundary rule is ahead of nearest TSS.",
        "NEAREST_AHEAD": "On this endpoint nearest TSS is ahead of the boundary rule.",
        "NO_DIFFERENCE_DETECTED": ("No difference detected between the boundary rule and nearest TSS."),
        "meaning_BOUNDARY_RULE_AHEAD": (
            "No rule names a host gene on this population, so the lead does not come from the element "
            "sitting inside the named gene's transcribed body, which is what section 8 could not rule "
            "out and section 9 had too few elements to ask."
        ),
        "meaning_NEAREST_AHEAD": (
            "Where the nearest 5' end lies across a measured boundary, the gene across it matches the "
            "element's cell type more often than the gene inside the domain."
        ),
        "meaning_NO_DIFFERENCE_DETECTED": (
            "With {k} elements differing, a difference smaller than about {mde} points could go "
            "undetected; section 9 could only bound it at 2.9 points, on {fly} differing elements."
        ),
        "gene_structure": (
            "Against the gene-structure rule on the same elements ({kg} of {n}, {cg}%; the two differ "
            "on {k2}, and neither names a host here), the difference is {d2} points, 95% {lo2} to "
            "{hi2}:"
        ),
        "gene_structure_BOUNDARY_RULE_AHEAD": ("the boundary rule is ahead of the gene-structure rule."),
        "gene_structure_GENE_STRUCTURE_AHEAD": ("the gene-structure rule is ahead of the boundary rule."),
        "gene_structure_NO_DIFFERENCE_DETECTED": "no difference detected between the two nodes.",
        "NOT_FALSIFIED": (
            "The boundary rule's genes match the element's cell type more often than genes at a "
            "matched distance ({cc}%, one-sided p {p}), so position carries cell-type information here."
        ),
        "FALSIFIED": (
            "The boundary rule's genes match the element's cell type no more often than genes at a "
            "matched distance ({cc}%, one-sided p {p}): on this endpoint the boundary rule carries no "
            "cell-type information beyond distance, and it is falsified as a carrier of it."
        ),
        "rotation_ahead": (
            "With the same boundaries rotated along each arm the rule averages {rm}% on these "
            "elements, below the measured boundaries' (one-sided p {pr}), so where the boundaries lie, "
            "not only how many there are, carries information here."
        ),
        "rotation_not_ahead": (
            "With the same boundaries rotated along each arm the rule averages {rm}% on these elements "
            "(one-sided p {pr}), so the measured positions are not shown to carry more than the "
            "boundaries' number and spacing."
        ),
        "all_distal": (
            "On all {N} distal elements (beside, deciding nothing; the rules differ on {kA}) the "
            "boundary rule is concordant for {kbA} ({cbA}%) against {knA} ({cnA}%) for nearest TSS, "
            "{dA} points [{loA}, {hiA}], and against {kgA} ({cgA}%) for the gene-structure rule, "
            "{dgA} points [{logA}, {higA}]."
        ),
        "missing": (
            "{miss} of the attributed genes on this population are absent from the RNA atlas's rows "
            "and are counted as not concordant."
        ),
        "tail": (
            "Accessibility is not function, and a boundary is not a regulatory unit by itself: this "
            "says the attributed gene is expressed where the element is open, not that the element "
            "acts on it. The boundaries come from 2-3 h embryos while the elements and the RNA come "
            "from larval-to-adult brain, a wider stage gap than section 9's; accessibility and RNA "
            "come from different nuclei and different studies; and the element source states no "
            "licence."
        ),
    },
}


# ------------------------------------------------------------------------------------------- the run


def score(
    chroms: dict[str, Chromosome],
    main: list[dict[str, Any]],
    alternatives: dict[str, list[dict[str, Any]]],
    bounds: dict[str, list[int]],
    arm: dict[str, int],
    expr: dict[str, tuple[float, ...]],
    funnel: dict[str, int],
) -> dict[str, Any]:
    """The registered run. `main` and `alternatives` come from `attribute`; `expr` from this module."""
    cen = census(main, alternatives, funnel)
    if cen["stop_rule_met"]:
        return {"outcome": "UNDERPOWERED", "census": cen}
    pop = between_genes(main)
    b_, n_, g_ = "boundary_gene", "nearest_gene", "domain_gene"
    half = BOUNDARY_SETS[PRIMARY_SET]["half_bin"]
    idx = ControlIndex(chroms, expr)
    first = compare(pop, expr, b_, n_)
    probs = [idx.probability(r["chrom"], r["tissue"], r["mid"], r[b_], r["boundary_dist"]) for r in pop]
    control = simpler.control_test(probs, first["a_concordant"])
    beside: dict[str, Any] = {
        "population_where_the_rules_differ": compare([r for r in pop if r[b_] != r[n_]], expr, b_, n_),
        "population_gene_structure_vs_nearest_tss": compare(pop, expr, g_, n_),
    }
    for name, alt in sorted(alternatives.items()):
        beside[f"{name}_vs_nearest_tss"] = compare(between_genes(alt), expr, b_, n_)
    for lab in LABELS:
        sub = [r for r in pop if r["tissue"] == lab]
        if len(sub) >= MIN_LABEL:
            beside["label_" + lab] = compare(sub, expr, b_, n_)
    res: dict[str, Any] = {
        "census": cen,
        "boundaries": sum(len(v) for v in bounds.values()),
        "primary": first,
        "primary_vs_gene_structure": compare(pop, expr, b_, g_),
        "all_distal": compare(distal(main), expr, b_, n_),
        "all_distal_vs_gene_structure": compare(distal(main), expr, b_, g_),
        "control": control,
        "positive_control": positive_control(main, idx, expr),
        "rotation_beside": rotation_null(pop, chroms, bounds, arm, expr, first["a_concordant"], half),
        "beside": beside,
        "attributed_genes_missing_from_rna_rows": sum(
            1 for r in pop for k in (b_, n_, g_) if r[k] not in expr
        ),
        "flynode": FLYNODE,
    }
    res["comparison"] = fn.label(first["delta_ci95_points"], "NEAREST_AHEAD")
    res["comparison_vs_gene_structure"] = fn.label(
        res["primary_vs_gene_structure"]["delta_ci95_points"], "GENE_STRUCTURE_AHEAD"
    )
    res["falsifier"] = "FALSIFIED" if control["p_one_sided"] >= ALPHA else "NOT_FALSIFIED"
    if not res["positive_control"]["responsive"]:
        res["outcome"] = "READOUT_UNRESPONSIVE"
        res["reading"] = (
            "The positive control did not separate from its matched control, so the readout is "
            "declared unresponsive and the primary reading is withheld."
        )
    else:
        res["outcome"] = res["comparison"]
        res["reading"] = reading(res)
    return res
