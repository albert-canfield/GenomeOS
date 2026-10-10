# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Closure at the cell level: element -> gene -> what the cell cannot live without.

The gene-level closure (`attribution/closure.py`) asked whether the attributed elements of a gene
go with the RNA the cell makes. This module asks the next question down the chain, the one area I's
Next list names: when the project says element e drives gene G in cell C, and G is a gene cell C
cannot live without, does silencing e cost the cell its growth?

Three independent statements are joined, and only the first is this project's:

* the project's own claim -- element e's deletion changes gene G on K562's own RNA track
  (`enhancer_targets_all_<chrom>`, the per-cell value `predicted_coding_by_cell["K562"]`, read
  through `attribution/targets.py` and signed the way `closure.effect_of` signs it);
* what K562 cannot live without -- the Chronos gene effect of G in K562 (ACH-000551) from the
  DepMap 24Q4 CRISPR knockout screens, under the release `therapeutics/depmap.py` pins;
* what happens when e is silenced -- the abundance of the guide RNAs that target e in ENCODE's
  genome-wide CRISPRi proliferation screen in K562 (ENCSR162HMH) against the same library in K562
  cells without the silencing machinery (ENCSR625TAF, ENCODE's own control for that screen).

The third is the cell-level readout and no part of this project has ever read it. The chain is
falsifiable in the direction that matters: a gene's essentiality is not a property of the element,
and proximity is not attribution, so the arms below hold the essentiality of the element's NEAREST
coding gene fixed while the attributed gene's essentiality varies. An element whose neighbourhood
contains a gene K562 depends on would deplete under any assignment; an element whose attribution
names such a gene where its neighbourhood does not is where the attribution does the work.

WHAT THIS CANNOT ESTABLISH, before anything is counted. A guide's abundance falling in a pooled
screen is a measurement over the whole culture, not a measurement of gene G; it cannot say the
silenced element acts through G rather than through anything else within reach. CRISPRi spreads
over roughly a kilobase of chromatin, so an element's measurement is not separable from its
immediate neighbours'. And the readout is relative guide abundance at harvest in CRISPRi-modified
cells against unmodified cells carrying the same library: the screen's own contrast, not a growth
rate and not a count of doublings.
"""

from __future__ import annotations

import bisect
import csv
import gzip
import hashlib
import io
import json
import math
import urllib.request
from pathlib import Path
from typing import Any

from genomeos.genome.annotation import Annotation, default_gencode
from genomeos.therapeutics import depmap

#: The cell line. One line, because the chain needs the prediction, the dependency and the screen
#: to be about the SAME cell, and K562 is the only cell with all three in this project's reach.
CELL = "K562"
#: DepMap's id for K562 (`StrippedCellLineName == "K562"` in the 24Q4 model table).
K562_MODEL = "ACH-000551"

#: The CRISPRi proliferation screen and ENCODE's own control for it. Both are K562 and both carry
#: the same genome-wide guide library against K562 DNase peaks; the screen's cells carry dCas9-KRAB
#: (`applied_modifications`: CRISPRi) and the control's do not.
SCREEN = {
    "experiment": "ENCSR162HMH",
    "assay": "proliferation CRISPR screen (CRISPRi), K562",
    "files": {  # biological replicate -> ENCODE file accession, GRCh38 guide quantifications
        1: "ENCFF945GIJ",
        2: "ENCFF995GEZ",
        3: "ENCFF141MGA",
        4: "ENCFF100ABB",
    },
}
CONTROL = {
    "experiment": "ENCSR625TAF",
    "assay": "Control CRISPR screen, K562: the same library in cells without dCas9-KRAB",
    "files": {
        1: "ENCFF902VHH",
        2: "ENCFF743IPS",
        3: "ENCFF482NVE",
        4: "ENCFF209ITV",
    },
}
#: The portal's md5 for every file above, so a streamed copy is checked against what ENCODE
#: published rather than against itself.
FILE_MD5 = {
    "ENCFF945GIJ": "85e514674741eac7ae48a48f220da67b",
    "ENCFF995GEZ": "5343e43b5969bd1e1cb8608a31ef87a9",
    "ENCFF141MGA": "5bf218f98ff466ea3ff8e9fdf544e6b2",
    "ENCFF100ABB": "04ffc8a74b54a0f889db4a002780ea4b",
    "ENCFF902VHH": "ee1aa18adb6024e2264ba3a526724894",
    "ENCFF743IPS": "2105a046ddd9b4e101349b5213d4e2ed",
    "ENCFF482NVE": "e64121069019df41eae981269e7ce0c6",
    "ENCFF209ITV": "2caa266ad796924082f030adffd094a9",
}
ENCODE_FILE = "https://www.encodeproject.org/files/{a}/@@download/{a}.tsv"
ENCODE_LICENCE = "ENCODE data use policy (no restrictions)"
SCREEN_SOURCE = (
    "ENCODE4 functional characterisation (wgCERES, Gersbach and Reddy labs): ENCSR162HMH, a "
    "CRISPRi proliferation screen in K562, and ENCSR625TAF, ENCODE's control screen for it"
)

#: The library's design as ENCODE publishes it: one 3 bp target interval per guide, GRCh38.
#: Positions only -- it carries no measurement, which is why the population may be built from it
#: before the registration is committed.
GUIDE_BED = Path("data/knowledge/closure_cell/ENCFF630CYG.bed.gz")
GUIDE_BED_MD5 = "ab78d6efd788631ff3e359860405ef33"
GUIDE_BED_SOURCE = (
    "ENCODE ENCSR366ESN (guide locations for the whole-genome CERES screen in K562), file ENCFF630CYG, GRCh38"
)

#: Where the distilled K562 dependency column and the frozen population are kept. Both sit under
#: data/knowledge, which this repository never commits.
STORE = Path("data/knowledge/closure_cell")
DEPMAP_K562 = STORE / "depmap_k562_gene_effect.json.gz"
POPULATION = STORE / "closure_cell_population.json.gz"

# ---- the thresholds, every one of them fixed here and read from here ---------------------------

#: An element enters the population when at least this many of the library's guides fall in it.
MIN_GUIDES = 3
#: The magnitude on K562's own track at which the project's prediction counts as a claim at all.
#: The bar `targets.genes_at_bar` already uses, so this module invents no new one.
MIN_EFFECT = 0.1
#: Chronos gene effect: at or below ESSENTIAL, K562 is called dependent on the gene; at or above
#: NOT_ESSENTIAL, not dependent. The gap is deliberate -- a gene between the two is in no arm,
#: because a borderline dependency would blur the only contrast this test has.
ESSENTIAL = -0.5
NOT_ESSENTIAL = -0.1
#: An element within this distance of any protein-coding TSS is out of the main population: a
#: CRISPRi guide at a promoter silences the gene directly and the chain under test is the distal
#: one. Those elements are kept apart, as the readout's positive control.
PROMOTER_WINDOW = 2000
#: A pseudocount on each element's pooled count, so an element the screen depleted to zero keeps a
#: finite log ratio instead of dropping out of the arm it belongs to.
PSEUDOCOUNT = 1.0
#: An element is measured only if the control side, pooled over its four replicates, saw at least
#: this many reads for its guides. Below that the ratio is a statement about sequencing depth. Set
#: before any count was read, from the shape of the library (about a million guides) rather than
#: from its depth, which is why it is deliberately low.
MIN_CONTROL_COUNT = 20
#: The depletion rate's cut is the floor arm's own 5th percentile, so the floor arm's rate is 5%
#: BY CONSTRUCTION and every other arm's rate is read against that. A rate is never reported
#: without its denominator, and never without the mean and median readout beside it.
DEPLETION_PERCENTILE = 5.0

ARMS = (
    "attributed_essential_neighbourhood_not",
    "attributed_not_neighbourhood_not",
    "attributed_essential_neighbourhood_essential",
    "attributed_not_neighbourhood_essential",
)
#: What each arm holds, in words, so a result's table cannot be read the wrong way round.
ARM_MEANING = {
    "attributed_essential_neighbourhood_not": "the project attributes the element to a gene K562 "
    "depends on, and the element's nearest coding gene is one K562 does not depend on: the "
    "attribution carries the claim and proximity does not",
    "attributed_not_neighbourhood_not": "neither the attributed gene nor the nearest coding gene "
    "is one K562 depends on: the floor",
    "attributed_essential_neighbourhood_essential": "both the attributed gene and the nearest "
    "coding gene are genes K562 depends on: the claim holds, but proximity would make it too",
    "attributed_not_neighbourhood_essential": "the nearest coding gene is one K562 depends on and "
    "the attributed gene is not: proximity carries the claim and the attribution does not",
}
POSITIVE_CONTROL = "promoter_of_a_gene_k562_depends_on"
DIRECTION_ARM = "negative_direction_attributed_essential"


# ---- the genome side: elements, guides, the nearest coding gene ---------------------------------


def element_effect(element: dict[str, Any]) -> float | None:
    """The element's predicted effect on its target gene on K562's OWN track, or None.

    Signed as `closure.effect_of` signs it: a fall in the gene on deleting the element (a negative
    log2 change) means the element holds the gene up in that cell, so the effect is POSITIVE. The
    per-cell value is used and never the tissue-agnostic magnitude -- that distinction is what
    turned the gene-level closure's verdict on 2026-09-13, and this arm needs the cell's own.
    """
    v = (element.get("predicted_coding_by_cell") or {}).get(CELL)
    return None if v is None else -float(v)


def bin_index(rows: list[tuple[int, int, Any]], bin_size: int = 10_000) -> dict[int, list[int]]:
    """Row indices by 10 kb bin: an overlap test over a million guides without an interval tree."""
    index: dict[int, list[int]] = {}
    for i, (start, end, _) in enumerate(rows):
        for b in range(start // bin_size, end // bin_size + 1):
            index.setdefault(b, []).append(i)
    return index


def overlapping(
    rows: list[tuple[int, int, Any]],
    index: dict[int, list[int]],
    start: int,
    end: int,
    bin_size: int = 10_000,
) -> list[int]:
    """Every row of `rows` whose interval overlaps [start, end), half-open, both sides 0-based."""
    hits: list[int] = []
    for b in range(start // bin_size, end // bin_size + 1):
        for i in index.get(b, ()):
            s, e, _ = rows[i]
            if s < end and e > start and i not in hits:
                hits.append(i)
    return hits


def guides_by_chrom(path: Path = GUIDE_BED) -> dict[str, list[tuple[int, int, str]]]:
    """The library's guide target intervals, by chromosome, from the design BED.

    The BED is the library as designed and carries no count, so reading it is not reading an
    outcome. The name is kept only so a guide can be traced back to the file. Read once and
    grouped, because the alternative -- one pass per chromosome -- reads 11.8 MB twenty-four times.
    """
    out: dict[str, list[tuple[int, int, str]]] = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line or line.startswith(("#", "track")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            out.setdefault(f[0], []).append((int(f[1]), int(f[2]), f[3]))
    return out


def index_from_elements(rows: list[dict[str, Any]]) -> GuideIndex:
    """The overlap index of a frozen population, built the same way by the census and by the run.

    One implementation, because an element that the census counted as covered and the run did not
    find would be a population that moved between the registration and the measurement.
    """
    by_chrom: dict[str, list[tuple[int, int, Any]]] = {}
    for row in rows:
        by_chrom.setdefault(row["chrom"], []).append((int(row["start"]), int(row["end"]), row["id"]))
    return {c: (rs, bin_index(rs)) for c, rs in by_chrom.items()}


def coding_tss(chrom: str) -> tuple[list[int], list[str]]:
    """Every protein-coding gene's TSS on one chromosome, as two parallel lists sorted by position.

    Two lists rather than a list of pairs because `nearest_coding` bisects the positions: twenty
    thousand elements against two thousand TSS is a scan this module cannot afford to do linearly
    on all twenty-four chromosomes.
    """
    path = default_gencode({chrom})
    if path is None:
        raise FileNotFoundError(
            f"no GENCODE annotation for {chrom} (data/reference/gencode_v50_{chrom}.gff3.gz)"
        )
    ann = Annotation.from_gff3(path, {chrom})
    pairs: list[tuple[int, str]] = []
    for g in ann.protein_coding():
        locus = g.locus
        plus = getattr(locus.strand, "value", "+") == "+"
        pairs.append((int(locus.start if plus else locus.end), g.symbol))
    pairs.sort()
    return [p for p, _ in pairs], [s for _, s in pairs]


def nearest_coding(positions: list[int], symbols: list[str], start: int, end: int) -> tuple[str | None, int]:
    """The nearest protein-coding TSS to an element, as (symbol, distance in bp).

    The distance is 0 when a TSS sits inside the element, half-open on the element's own interval,
    and otherwise counts from the element's nearest BASE to the TSS, so a TSS one base past the
    element's last base is at distance 1.
    """
    if not positions:
        return None, 1 << 62
    i = bisect.bisect_left(positions, start)
    best: tuple[str | None, int] = (None, 1 << 62)
    for j in (i - 1, i, i + 1):
        if 0 <= j < len(positions):
            pos = positions[j]
            d = 0 if start <= pos < end else min(abs(pos - start), abs(pos - (end - 1)))
            if d < best[1]:
                best = (symbols[j], d)
    return best


# ---- the dependency side: one row of DepMap, streamed and distilled ----------------------------


def k562_gene_effect(genes: set[str], net: bool = True, cache: Path = DEPMAP_K562) -> dict[str, Any]:
    """The Chronos gene effect of `genes` in K562 alone, from the cache or streamed once.

    `therapeutics/depmap.py` holds the release pin, the figshare id, the DOI and the licence, and
    its `matrix` keeps the requested columns for EVERY cell line; fifteen thousand genes across
    1,178 models would be gigabytes in memory for the one row this needs. So this reads the same
    pinned file through the same figshare id and keeps ONE row, K562's. Nothing else is held and
    nothing else is written. The digest recorded is over the bytes that were streamed, so the
    record names the file that was read and not a file with the same name.
    """
    cached = _read_gz(cache)
    wanted = {g.upper() for g in genes}
    if cached is not None and cached.get("release") == depmap.RELEASE:
        if wanted <= set(cached.get("requested") or ()):
            return {**cached, "from_cache": True}
        wanted |= set(cached.get("requested") or ())
    if not net:
        raise depmap.DepMapUnavailableError(
            f"no cached K562 gene effect covering {len(wanted)} genes and network disabled ({cache})"
        )
    file_id, name, what = depmap.FILES["gene_effect"]
    url = depmap.FIGSHARE_FILE.format(file_id)
    digest = hashlib.sha256()
    total = 0
    header: list[str] | None = None
    values: dict[str, float] = {}
    found = False
    request = urllib.request.Request(url, headers=depmap.UA)
    with urllib.request.urlopen(request, timeout=900) as response:  # noqa: S310
        pending = b""
        for chunk in iter(lambda: response.read(1 << 20), b""):
            digest.update(chunk)
            total += len(chunk)
            if found:
                continue  # the row is in hand; the rest is streamed only to finish the digest
            pending += chunk
            lines = pending.split(b"\n")
            pending = lines.pop()
            for raw in lines:
                row = next(csv.reader(io.StringIO(raw.decode("utf-8", "replace"))), [])
                if header is None:
                    header = [depmap.symbol(c) for c in row]
                    continue
                if not row or row[0] != K562_MODEL:
                    continue
                for i, cell in enumerate(row[1:], start=1):
                    if i < len(header) and header[i] in wanted and cell not in ("", "NA"):
                        values[header[i]] = round(float(cell), 4)
                found = True
                break
    payload = {
        "dataset": "gene_effect",
        "file": name,
        "figshare_file": file_id,
        "what": what,
        "model": K562_MODEL,
        "cell_line": CELL,
        "release": depmap.RELEASE,
        "doi": depmap.DOI,
        "licence": depmap.LICENCE,
        "citation": depmap.CITATION,
        "streamed_bytes": total,
        "streamed_sha256": digest.hexdigest(),
        "requested": sorted(wanted),
        "genes": sorted(values),
        "missing": sorted(wanted - set(values)),
        "values": values,
        "row_found": found,
    }
    _write_gz(cache, payload)
    return {**payload, "from_cache": False}


def dependency(effects: dict[str, float], gene: str | None) -> str:
    """`"essential"`, `"not"`, `"borderline"` or `"unmeasured"` for one gene in K562."""
    if not gene:
        return "unmeasured"
    v = effects.get(gene.upper())
    if v is None:
        return "unmeasured"
    if v <= ESSENTIAL:
        return "essential"
    if v >= NOT_ESSENTIAL:
        return "not"
    return "borderline"


def arm_of(attributed: str, neighbourhood: str, effect: float) -> str | None:
    """The arm an element belongs to, or None when it belongs to none of them.

    An element is in no arm when either dependency call is borderline or unmeasured, or when the
    prediction's magnitude on K562's own track is under `MIN_EFFECT`. The direction arm is separate
    and reported apart: a prediction whose sign says the element lowers the gene in K562 makes no
    depletion prediction, and reading it beside the others would mix two claims.
    """
    if attributed not in ("essential", "not") or neighbourhood not in ("essential", "not"):
        return None
    if abs(effect) < MIN_EFFECT:
        return None
    if effect < 0:
        return DIRECTION_ARM if attributed == "essential" else None
    return f"attributed_{attributed}_neighbourhood_{neighbourhood}"


# ---- the readout, once the counts exist ---------------------------------------------------------


def log2_ratio(screen: float, screen_total: float, control: float, control_total: float) -> float | None:
    """An element's readout: log2 of its share of the CRISPRi library over its share of the control.

    Counts per million on each side with `PSEUDOCOUNT` added to the element, so an element the
    screen emptied keeps a finite, strongly negative value rather than leaving its arm. None when
    either library total is zero, which can only mean a file was not read.
    """
    if screen_total <= 0 or control_total <= 0:
        return None
    a = (screen + PSEUDOCOUNT) / screen_total * 1e6
    b = (control + PSEUDOCOUNT) / control_total * 1e6
    if a <= 0 or b <= 0:
        return None
    return math.log2(a / b)


def mean(xs: list[float]) -> float | None:
    return sum(xs) / len(xs) if xs else None


def median(xs: list[float]) -> float | None:
    if not xs:
        return None
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2


def percentile(xs: list[float], q: float) -> float | None:
    """The q-th percentile by linear interpolation between order statistics, or None if empty."""
    if not xs:
        return None
    s = sorted(xs)
    if len(s) == 1:
        return s[0]
    pos = (len(s) - 1) * q / 100.0
    low = int(math.floor(pos))
    high = min(low + 1, len(s) - 1)
    return s[low] + (s[high] - s[low]) * (pos - low)


def mann_whitney_u(a: list[float], b: list[float]) -> dict[str, Any]:
    """A rank test of `a` against `b`, one-sided for a below b, by the normal approximation.

    Written out rather than imported: scipy is an optional extra in this repository and a result
    whose number needs an extra cannot be produced in the bare environment. Ties take mid-ranks and
    the variance carries the tie correction, so a readout with many equal values is not given a
    p-value computed as though it had none.
    """
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0:
        return {"u": None, "p_one_sided_a_below_b": None, "n_a": n1, "n_b": n2}
    pooled = sorted([(v, 0) for v in a] + [(v, 1) for v in b])
    ranks: list[float] = [0.0] * len(pooled)
    ties: list[int] = []
    i = 0
    while i < len(pooled):
        j = i
        while j + 1 < len(pooled) and pooled[j + 1][0] == pooled[i][0]:
            j += 1
        r = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[k] = r
        ties.append(j - i + 1)
        i = j + 1
    r1 = sum(r for r, (_, which) in zip(ranks, pooled, strict=True) if which == 0)
    u1 = r1 - n1 * (n1 + 1) / 2
    n = n1 + n2
    mu = n1 * n2 / 2
    tie_term = sum(t**3 - t for t in ties)
    var = n1 * n2 / 12 * ((n + 1) - tie_term / (n * (n - 1))) if n > 1 else 0.0
    if var <= 0:
        return {"u": u1, "p_one_sided_a_below_b": None, "n_a": n1, "n_b": n2}
    z = (u1 - mu) / math.sqrt(var)
    return {
        "u": u1,
        "z": round(z, 4),
        # P(U1 <= u1) under the null: small when the values of `a` rank below those of `b`
        "p_one_sided_a_below_b": 0.5 * math.erfc(-z / math.sqrt(2)),
        "n_a": n1,
        "n_b": n2,
    }


# ---- the screen's own files, streamed -----------------------------------------------------------

#: The guide-quantification layout, by position, as the files themselves show it: seventeen
#: tab-separated fields, no header line. The chromosome and a 3 bp target interval first, then the
#: BED name, the replicate's count, the strand, the perturbation target's id, six unused fields,
#: the spacer and guide sequences, the guide's type and a note. A row that does not fit is counted,
#: never guessed at.
#:
#: These indices were WRONG when this module was first written, and the census' layout probe is
#: what found it: the BED name is `NA|chr1:605582-605585:+`, one field holding a pipe, and reading
#: the file through a `tr '\t' '|'` pipeline had split it in two, which pushed every index past it
#: up by one. The probe reported 17 fields on 13,978 of 13,978 rows and 0 rows matching the
#: declared layout, before a single count had been read. The malformed-row gate in the run would
#: have refused the file too, since the strand sat where the count was expected and does not parse
#: as a number -- but it would have refused after streaming 173 MB, and the probe cost 2 MB.
FIELDS = {
    "chrom": 0,
    "start": 1,
    "end": 2,
    "name": 3,
    "count": 4,
    "strand": 5,
    "target_id": 6,
    "spacer": 13,
    "guide": 14,
    "guide_type": 15,
    "notes": 16,
}
FIELD_COUNT = 17
TARGETING = "targeting"
GuideIndex = dict[str, tuple[list[tuple[int, int, Any]], dict[int, list[int]]]]


def parse_guide_row(line: str) -> tuple[str, int, int, float, str] | None:
    """One guide-quantification row as (chrom, start, end, count, guide_type), or None if malformed."""
    f = line.rstrip("\n").split("\t")
    if len(f) != FIELD_COUNT:
        return None
    try:
        start = int(f[FIELDS["start"]])
        end = int(f[FIELDS["end"]])
        count = float(f[FIELDS["count"]])
    except ValueError:
        return None
    chrom = f[FIELDS["chrom"]]
    if not chrom.startswith("chr") or end <= start:
        return None
    return chrom, start, end, count, f[FIELDS["guide_type"]]


def new_census() -> dict[str, Any]:
    return {
        "rows": 0,
        "malformed": 0,
        "targeting": 0,
        "non_targeting": 0,
        "outside_the_population": 0,
        "unplaced_chromosome": 0,
        "count_outside_the_population": 0.0,
    }


def fold_row(
    line: str,
    index: GuideIndex,
    per_element: dict[str, float],
    guides_per_element: dict[str, int],
    census: dict[str, Any],
    non_targeting: dict[str, float],
) -> None:
    """One row folded into the aggregates. Every row is counted somewhere, the bad ones included."""
    if not line.strip():
        return
    census["rows"] += 1
    parsed = parse_guide_row(line)
    if parsed is None:
        census["malformed"] += 1
        return
    chrom, start, end, count, guide_type = parsed
    if guide_type != TARGETING:
        census["non_targeting"] += 1
        non_targeting["guides"] += 1
        non_targeting["count"] += count
        return
    census["targeting"] += 1
    held = index.get(chrom)
    if held is None:
        census["unplaced_chromosome"] += 1
        census["count_outside_the_population"] += count
        return
    rows, bins = held
    hits = overlapping(rows, bins, start, end)
    if not hits:
        census["outside_the_population"] += 1
        census["count_outside_the_population"] += count
        return
    for i in hits:
        element_id = rows[i][2]
        per_element[element_id] = per_element.get(element_id, 0.0) + count
        guides_per_element[element_id] = guides_per_element.get(element_id, 0) + 1


def stream_guide_counts(accession: str, index: GuideIndex, log=None) -> dict[str, Any]:
    """Stream one guide-quantification file and keep only what the population needs.

    Per element: the summed count of the targeting guides inside it, and how many guides that was.
    Per file: the library total over every guide, the non-targeting guides pooled as the assay's
    own null, the row census, and the sha256 and md5 of the bytes that were streamed. The rows are
    discarded as they are read; nothing of these files reaches disk.
    """
    url = ENCODE_FILE.format(a=accession)
    per_element: dict[str, float] = {}
    guides_per_element: dict[str, int] = {}
    census = new_census()
    non_targeting = {"guides": 0.0, "count": 0.0}
    sha = hashlib.sha256()
    md5 = hashlib.md5()  # noqa: S324 - compared against the md5 ENCODE publishes; not a defence
    streamed = 0
    request = urllib.request.Request(url, headers={"User-Agent": "GenomeOS/0.1 (closure-cell)"})
    with urllib.request.urlopen(request, timeout=1800) as response:  # noqa: S310
        pending = b""
        for chunk in iter(lambda: response.read(1 << 20), b""):
            sha.update(chunk)
            md5.update(chunk)
            streamed += len(chunk)
            pending += chunk
            lines = pending.split(b"\n")
            pending = lines.pop()
            for raw in lines:
                fold_row(
                    raw.decode("utf-8", "replace"),
                    index,
                    per_element,
                    guides_per_element,
                    census,
                    non_targeting,
                )
        if pending.strip():
            fold_row(
                pending.decode("utf-8", "replace"),
                index,
                per_element,
                guides_per_element,
                census,
                non_targeting,
            )
    library_total = (
        sum(per_element.values()) + non_targeting["count"] + census["count_outside_the_population"]
    )
    if log:
        print(f"  {accession}: {census['rows']} rows, {streamed} bytes streamed", file=log, flush=True)
    return {
        "accession": accession,
        "url": url,
        "streamed_bytes": streamed,
        "sha256": sha.hexdigest(),
        "md5": md5.hexdigest(),
        "md5_matches_encode": md5.hexdigest() == FILE_MD5.get(accession),
        "census": census,
        "library_total": library_total,
        "non_targeting": non_targeting,
        "per_element": per_element,
        "guides_per_element": guides_per_element,
    }


# ---- the frozen population ----------------------------------------------------------------------


def _read_gz(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    try:
        with gzip.open(path, "rt") as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError):
        return None


def _write_gz(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        json.dump(payload, fh, separators=(",", ":"))


def save_population(rows: list[dict[str, Any]], meta: dict[str, Any], path: Path = POPULATION) -> str:
    """Write the frozen population and return its sha256, which the registration pins."""
    _write_gz(path, {**meta, "elements": rows})
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_population(path: Path = POPULATION) -> dict[str, Any]:
    held = _read_gz(path)
    if held is None:
        raise FileNotFoundError(f"no frozen population at {path}; run scripts/closure_cell_inventory.py")
    return held
