# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Simpler genomes, area I: the domain rule in a MAMMAL, and the between-gene stratum above all.

The domain rule attributes a non-coding element to a gene. Area I has asked it in two
invertebrates and they DISAGREE about the one case that matters:

* **Worm** (docs/SIMPLER-GENOMES.md sections 4 and 12, `eede7d9`, `b97a8d5`): the domain rule is
  ahead of nearest TSS by **+2.89 points [+1.83, +3.91]** overall, and ahead by
  **+3.73 points [+2.79, +4.74]** on the 2,335 elements that lie BETWEEN genes, where no host gene
  can be named -- a lead that **survives** the worm's own outron and trans-splicing caveat at
  **+4.47 [+3.16, +5.78]** on the 1,745 survivors of that exclusion.
* **Fly** (sections 8, 9 and 11, `6bfc956`, `b0d05f9`, `ce168cb`): **+2.70** overall, but between
  genes **+0.98 [-0.42, +2.31]** on Calderon et al. 2022's elements and **+0.12 [-0.25, +0.61]** on
  Janssens et al. 2022's 3,265 brain elements -- and the lead appears only where a host gene CAN be
  named (on all 10,052 distal elements the gene-structure rule is ahead of the boundary rule by
  -2.25 [-2.88, -1.66]), which is consistent with the fly's lead coming from the host gene's own
  transcription rather than from regulation.

**Nobody has asked a mammal.** This module asks it on the mouse, on one source for both
modalities, with the worm's and the fly's registration carried verbatim: the same endpoint, the
same matched-distance control, the same rotated-boundary control, the same 1-Mb cluster bootstrap,
the same `MIN_DISCORDANT = 200` stop rule, the same reading wording. Only the organism changes.

**The stop list of 2026-09-27 was read, not ignored.** docs/ROADMAP.md section 5 reinstates "No
more biosamples: 1.1 needs its readings to mean what their names say." That clause's subject is the
**human reader's 13 biosamples** -- `READER_TERMS` in `genomeos/genome/context_evidence.py`. This
lane adds none, names none and widens none. **A mouse tissue source is not a human reader
biosample**, exactly as `lane-flysource` recorded for its fly source (section 11). The mouse
accessibility, RNA and boundary sets named here are a comparative-genomics population, not a
reading of a human element.

Wording, from the worm's and the fly's registrations and not loosened here: conservation and
accessibility are **not** function. The endpoint says **the attributed gene is expressed where the
element is open** -- never that the element regulates the gene.

Imports `simpler.py` (the worm node, the bootstrap, the matched control), `simpler_fly.py` (the
endpoint and its CPM threshold) and `simpler_flysource.py` (the bin-parameterised boundary rule and
the position-only move count), all unchanged.
"""

from __future__ import annotations

import bisect
import gzip
import json
import math
import random
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from genomeos.attribution import simpler
from genomeos.attribution import simpler_fly as sf
from genomeos.attribution import simpler_flynode as fn
from genomeos.attribution import simpler_flysource as fs
from genomeos.attribution.simpler import Chromosome, Gene

# ------------------------------------------------------------------------------------- what is read

STORE = Path("data/knowledge/simpler_mouse")
GFF3 = STORE / "gencode.vM21.annotation.gff3.gz"
CCRE_FILE = STORE / "mm10-cCREs.bed"
DOMAIN_FILE = STORE / "bonev2017_TableS2.xlsx"
FETCH_RECORD = STORE / "fetch_record.json"

ASSEMBLY = "mm10 (GRCm38): the gene models, the cCREs, the accessibility peaks and the boundaries"
LIFTOVER = "none: every input is distributed in mm10/GRCm38 and no coordinate is converted here"

GENE_MODELS = (
    "GENCODE vM21 (GRCm38/mm10), gene records of type protein_coding -- M21 and not a later "
    "release because that is the annotation ENCODE quantified this panel's RNA on, so the gene "
    "positions and the RNA rows come from ONE annotation and no id has to be mapped between two"
)
GENE_MODELS_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.annotation.gff3.gz"
)
GENE_MODELS_LICENCE = (
    "GENCODE/Ensembl: no restriction on use of the annotation (EMBL-EBI terms of use); recorded "
    "from the GENCODE data-release page, not inferred"
)

ELEMENT_SOURCE = (
    "ENCODE SCREEN registry of candidate cis-regulatory elements, mouse mm10 v3 "
    "(ENCODE Project Consortium 2020, Nature 583:699; Moore et al.), 368,121 cCREs"
)
ELEMENT_URL = "https://downloads.wenglab.org/V3/mm10-cCREs.bed"
ELEMENT_LICENCE = (
    "ENCODE data-use policy as encodeproject.org states it: ENCODE data are free of restriction "
    "on use; recorded from the portal, not inferred"
)

#: The 21 chromosomes the population is drawn from. mm10's unplaced and alternate contigs are
#: excluded because a boundary set and a gene order are both undefined there.
CHROMOSOMES = tuple(f"chr{c}" for c in list(range(1, 20)) + ["X", "Y"])

# ------------------------------------------------------------------- the source, chosen before any
# ------------------------------------------------------------------- outcome and recorded as such

ACCESSIBILITY_SOURCE = (
    "ENCODE ATAC-seq, Mus musculus embryo (14.5 days), 12 tissues, Bing Ren lab (UCSD), ENCODE3; "
    "the file read per tissue is the experiment's own `replicated peaks` narrowPeak on mm10 "
    "(the ENCODE mouse developmental matrix; Gorkin et al. 2020, Nature 583:744)"
)
RNA_SOURCE = (
    "ENCODE polyA plus RNA-seq, Mus musculus embryo (14.5 days), the SAME 12 tissues, Barbara Wold "
    "lab (Caltech), ENCODE3; the file read per tissue is the experiment's own preferred-default "
    "`gene quantifications` TSV on mm10 with GENCODE M21, and the column read is TPM"
)
BOUNDARY_SOURCE = (
    "Bonev et al. 2017 (Cell 171:557-572.e24, DOI 10.1016/j.cell.2017.09.043, GEO GSE96107), "
    "Table S2 'Genomic Coordinates and Compartment Identity of the Topological Domains in Different "
    "Cell Types': the authors' own insulation-score domain call on mm10, read as the domain edges. "
    "The primary sheet is `ncx NPC` -- neural progenitors FACS-purified from EMBRYONIC NEOCORTEX, "
    "in vivo -- 2,398 domains. The other four sheets and three ENCODE e14.5 CTCF sets are beside"
)
BOUNDARY_URL = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417311376-mmc2.xlsx"
BOUNDARY_LICENCE = (
    "CC BY 4.0, verbatim from the article's front matter: 'This is an open access article under the "
    "CC BY license (http://creativecommons.org/licenses/by/4.0/)'"
)
BOUNDARY_METHOD = (
    "the authors' own, verbatim: 'The insulation score was computed ... at 1Kb resolution within a "
    "region of +/- 250Kb ... Domain boundaries were then defined as the local 2Kb maxima in "
    "regions, where the insulation score is above the 90% quantile of the genome-wide distribution. "
    "Boundaries within 10Kb of each other were merged and only boundaries, which were detected in "
    "at least 3 biological replicates, were retained. Domains shorter than 50Kb were discarded.' "
    "No caller and no threshold is chosen here"
)

#: Bonev's five sheets with the row count each must reproduce, read before any RNA value. A reader
#: that does not reproduce its sheet's count refuses, so a changed workbook cannot pass silently.
BONEV_SHEETS = {
    "bonev_2017_ncx_npc": ("ncx NPC", 2398, "neural progenitors, embryonic neocortex, in vivo"),
    "bonev_2017_ncx_cn": ("ncx CN", 2307, "cortical neurons, embryonic neocortex, in vivo"),
    "bonev_2017_es": ("ES", 2712, "embryonic stem cells, in vitro"),
    "bonev_2017_npc": ("NPC", 2522, "neural progenitors, in vitro"),
    "bonev_2017_cn": ("CN", 2442, "cortical neurons, in vitro"),
}

ENCODE_LICENCE = (
    "ENCODE data-use policy as encodeproject.org states it: released ENCODE data carry no "
    "restriction on use or publication; recorded from the portal, not inferred"
)
ENCODE_FILE_URL = "https://www.encodeproject.org/files/{acc}/@@download/{acc}.{ext}"

#: THE LABEL VOCABULARY: 12 ENCODE biosample term names at embryonic day 14.5, each of which has
#: BOTH an ATAC-seq and a polyA plus RNA-seq experiment at that same stage. That is the whole
#: vocabulary and it is the source's own: no tissue is grouped, renamed or split here, and no tissue
#: is nested inside another (ENCODE's e14.5 `brain` experiment is EXCLUDED because it contains
#: forebrain, midbrain and hindbrain, which are three of the twelve).
LABELS = (
    "embryonic facial prominence",
    "forebrain",
    "heart",
    "hindbrain",
    "intestine",
    "kidney",
    "limb",
    "liver",
    "lung",
    "midbrain",
    "neural tube",
    "stomach",
)
LABEL_KIND = "ENCODE embryonic-day-14.5 mouse tissues"
LABEL_NOUN = "tissue"

#: Accessibility, one experiment per label; the portal offers exactly one ATAC-seq experiment for
#: each of the twelve at e14.5, so no experiment is chosen here. file -> (experiment, file).
ATAC_FILES = {
    "embryonic facial prominence": ("ENCSR876SYO", "ENCFF426FDN"),
    "forebrain": ("ENCSR810HQR", "ENCFF922LUD"),
    "heart": ("ENCSR068YGC", "ENCFF923BLK"),
    "hindbrain": ("ENCSR798FDL", "ENCFF343ZQT"),
    "intestine": ("ENCSR150EOO", "ENCFF964BLN"),
    "kidney": ("ENCSR758IRM", "ENCFF049CMP"),
    "limb": ("ENCSR460BUL", "ENCFF899IHU"),
    "liver": ("ENCSR032HKE", "ENCFF133VPS"),
    "lung": ("ENCSR335VJW", "ENCFF991RKC"),
    "midbrain": ("ENCSR384JBF", "ENCFF870PSK"),
    "neural tube": ("ENCSR700QBR", "ENCFF522NJW"),
    "stomach": ("ENCSR618HDK", "ENCFF179ZMB"),
}

#: RNA, one experiment per label. Four tissues have more than one e14.5 polyA experiment on the
#: portal (heart 2, liver 5, limb 3); the rule that picks one is stated in advance and is not a
#: tiebreak: **the experiment of the panel's own series** -- Barbara Wold, Caltech, ENCODE3 -- which
#: exists exactly once for every one of the twelve. Within an experiment the file is the portal's
#: own `preferred_default` gene quantification on GENCODE M21.
RNA_FILES = {
    "embryonic facial prominence": ("ENCSR823VEE", "ENCFF932AJG"),
    "forebrain": ("ENCSR185LWM", "ENCFF270DCV"),
    "heart": ("ENCSR727FHP", "ENCFF554PCI"),
    "hindbrain": ("ENCSR559TRB", "ENCFF304ILZ"),
    "intestine": ("ENCSR932TRU", "ENCFF609AAY"),
    "kidney": ("ENCSR504GEG", "ENCFF449ACY"),
    "limb": ("ENCSR216NEG", "ENCFF394CND"),
    "liver": ("ENCSR867YNV", "ENCFF578AKP"),
    "lung": ("ENCSR039ADS", "ENCFF722NZT"),
    "midbrain": ("ENCSR343YLB", "ENCFF954BEO"),
    "neural tube": ("ENCSR928OXI", "ENCFF457PAJ"),
    "stomach": ("ENCSR290RRR", "ENCFF361QMV"),
}

#: The boundary sets. Every CTCF ChIP-seq experiment ENCODE has released at embryonic day 14.5 --
#: three of them -- and nothing else, so no stage gap is opened and no set is dropped on an outcome.
CTCF_FILES = {
    "ctcf_e14_5_liver": ("ENCSR397RHW", "ENCFF123IQA", "liver, Richard Myers lab (HAIB), ENCODE3"),
    "ctcf_e14_5_lung": ("ENCSR677SIH", "ENCFF558ZIE", "lung, Richard Myers lab (HAIB), ENCODE3"),
    "ctcf_e14_5_brain": ("ENCSR000CEH", "ENCFF756OJJ", "brain, Bing Ren lab (UCSD), ENCODE2-Mouse"),
}

#: The primary boundary set, fixed before any RNA value is read: the ONE author-called mouse Hi-C
#: domain set that is mm10 (no liftover), carries a stated licence, and was measured in EMBRYONIC
#: tissue -- Bonev's `ncx NPC`, neural progenitors FACS-purified from embryonic neocortex. Of its
#: five sheets it is the in-vivo progenitor population of the developing cortex, which is the
#: closest of the five to this panel's stage and to three of its twelve tissues (forebrain,
#: midbrain, hindbrain). The other four sheets are reported beside, so nothing is dropped.
#:
#: Section 11 registered its stage gap -- 2-3 h fly embryos against larval-to-adult brain elements
#: -- as the price of a large enough population, and named a tissue-matched registration as the
#: first thing that should follow. Here the boundaries and the elements are both embryonic mouse.
PRIMARY_SET = "bonev_2017_ncx_npc"

#: The bin, carried verbatim from `simpler_flynode` (500 bp, half 250): the registered rule counts a
#: gene's 5' end inside either flanking boundary's bin as inside the domain. NOT retuned for the
#: mouse, so the only thing that differs from the fly's run is the organism and the source.
BIN_BP = fn.BIN_BP
HALF_BIN = fn.HALF_BIN

#: What a boundary IS here, recorded rather than chosen: Bonev publishes DOMAINS, and the boundary
#: positions are their start and end coordinates, deduplicated and sorted. That is the only
#: conversion, it has no parameter, and the CTCF sets beside take each peak span's midpoint instead
#: because a CTCF peak is already a point-like call.
BOUNDARY_CAVEAT = (
    "Bonev's domains are read as their edges -- every domain start and end, deduplicated -- which "
    "is the only conversion and has no parameter. The three ENCODE e14.5 CTCF sets reported beside "
    "take every conservative-IDR peak as a boundary, which is the brief's CTCF option and not a "
    "Hi-C domain call. 4D Nucleome's 110 mm10 boundary BEDs are NOT used: they are DCIC-pipeline "
    "calls rather than author calls, the portal requires a login, and the portal states no licence"
)

#: SCREEN's own cCRE classes. The element side reads this column and nothing else from the registry:
#: the classes are the source's call, not a threshold chosen here.
PROMOTER_CLASSES = ("PLS", "pELS", "DNase-H3K4me3")
DISTAL = sf.DISTAL
PROMOTER = sf.PROMOTER

TRANSCRIPT_TYPES = frozenset({"protein_coding"})

# ---------------------------------------------------------------- carried verbatim, nothing retuned

#: The worm's threshold and the worm's unit, not the fly's: ENCODE quantifies this panel's RNA as
#: TPM, which is what `simpler.DETECTION_TPM` is, so the endpoint is the WORM's registered one
#: (5 TPM or more and highest in the element's own label) with nothing retuned. The fly's identical
#: 5.0 was in CPM because its source counted UMIs; the number is the same and the unit is named.
DETECTION_TPM = simpler.DETECTION_TPM
DISTANCE_BINS = simpler.DISTANCE_BINS
CONTROL_DRAWS = simpler.CONTROL_DRAWS
BOOTSTRAP = simpler.BOOTSTRAP
CLUSTER_BP = simpler.CLUSTER_BP
SEED = simpler.SEED
MIN_DISCORDANT = simpler.MIN_DISCORDANT
ALPHA = simpler.ALPHA
ROTATION_DRAWS = fn.ROTATION_DRAWS
MIN_LABEL = fn.MIN_TISSUE
STOP_RULE = MIN_DISCORDANT

#: What the two invertebrates found, for the reading to name. Registered figures, word for word
#: from the commits that carry them; none is recomputed here and none is strengthened.
WORM = {
    "overall_delta_points": 2.89,
    "overall_ci95_points": [1.83, 3.91],
    "between_gene_delta_points": 3.73,
    "between_gene_ci95_points": [2.79, 4.74],
    "between_gene_elements": 2335,
    "outron_excluded_delta_points": 4.47,
    "outron_excluded_ci95_points": [3.16, 5.78],
    "outron_excluded_elements": 1745,
    "registration": "595d7c2",
    "commits": ["eede7d9", "b97a8d5"],
}
FLY = {
    "overall_delta_points": 2.70,
    "between_gene_delta_points": 0.98,
    "between_gene_ci95_points": [-0.42, 2.31],
    "brain_source_delta_points": 0.12,
    "brain_source_ci95_points": [-0.25, 0.61],
    "brain_source_elements": 3265,
    "all_distal_gene_structure_vs_boundary_points": -2.25,
    "all_distal_gene_structure_vs_boundary_ci95_points": [-2.88, -1.66],
    "registrations": ["fd81e37", "f5985f7", "ff063e6"],
    "commits": ["6bfc956", "b0d05f9", "ce168cb"],
}

# ----------------------------------------------------------------------------------- the gene models


_ATTR = re.compile(r"([^=;]+)=([^;]*)")


def _attrs(field: str) -> dict[str, str]:
    return {k: v for k, v in _ATTR.findall(field)}


def read_annotation(
    path: Path = GFF3, chromosomes: tuple[str, ...] = CHROMOSOMES
) -> tuple[dict[str, Chromosome], dict[str, list[int]], dict[str, str]]:
    """GENCODE vM21 -> the protein-coding genes per chromosome, every transcript 5' end, the names.

    `simpler_fly.read_annotation` for the mouse: the same two products, built from GENCODE's `gene`
    and `transcript` records instead of Ensembl's. Versioned GENCODE ids are reduced to the bare
    ENSMUSG accession, because the RNA side names genes that way.
    """
    genes: dict[str, list[Gene]] = {}
    tss: dict[str, list[int]] = {}
    names: dict[str, str] = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:  # type: ignore[operator]
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[0] not in chromosomes:
                continue
            kind = parts[2]
            if kind not in ("gene", "transcript"):
                continue
            attrs = _attrs(parts[8])
            if kind == "gene":
                if attrs.get("gene_type") != "protein_coding":
                    continue
                gid = attrs["gene_id"].split(".")[0]
                gene = Gene(
                    gid=gid,
                    name=attrs.get("gene_name", gid),
                    chrom=parts[0],
                    start=int(parts[3]),
                    end=int(parts[4]),
                    strand=parts[6],
                )
                genes.setdefault(parts[0], []).append(gene)
                names[gid] = gene.name
            elif attrs.get("transcript_type") in TRANSCRIPT_TYPES:
                start, end = int(parts[3]), int(parts[4])
                tss.setdefault(parts[0], []).append(start if parts[6] == "+" else end)
    chroms = {c: Chromosome(g) for c, g in sorted(genes.items())}
    return chroms, {c: sorted(v) for c, v in tss.items()}, names


def chrom_lengths(path: Path = CCRE_FILE, chromosomes: tuple[str, ...] = CHROMOSOMES) -> dict[str, int]:
    """The rotation control's arm lengths: the last cCRE end on each chromosome.

    `simpler_flynode.read_arm_lengths` takes the fly's arm lengths from the GFF3's region records;
    GENCODE's GFF3 carries no `##sequence-region` line, so the span the rotation turns is taken from
    the element registry itself -- the same quantity, read from a file this lane already declares.
    """
    out: dict[str, int] = {}
    with open(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if parts[0] not in chromosomes:
                continue
            out[parts[0]] = max(out.get(parts[0], 0), int(parts[2]))
    return out


# -------------------------------------------------------------------------------------- the elements


@dataclass(frozen=True)
class CCRE:
    chrom: str
    start: int
    end: int
    eid: str
    klass: str

    @property
    def mid(self) -> int:
        return (self.start + self.end) // 2

    @property
    def kind(self) -> str:
        """SCREEN's own class, reduced to this lane's two: promoter-proximal or distal."""
        return PROMOTER if any(k in self.klass for k in PROMOTER_CLASSES) else DISTAL


def read_ccres(path: Path = CCRE_FILE, chromosomes: tuple[str, ...] = CHROMOSOMES) -> list[CCRE]:
    """The registry, positions and SCREEN's own class column only. No activity value is read here."""
    out: list[CCRE] = []
    with open(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6 or parts[0] not in chromosomes:
                continue
            out.append(CCRE(parts[0], int(parts[1]), int(parts[2]), parts[4], parts[5]))
    return out


def peak_index(spans: dict[str, list[tuple[int, int]]]) -> dict[str, tuple[list[int], list[int]]]:
    """Merged, sorted peak spans per chromosome, as (starts, ends) for bisect."""
    out: dict[str, tuple[list[int], list[int]]] = {}
    for chrom, raw in spans.items():
        merged: list[list[int]] = []
        for s, e in sorted(raw):
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[chrom] = ([m[0] for m in merged], [m[1] for m in merged])
    return out


def overlaps(index: dict[str, tuple[list[int], list[int]]], chrom: str, start: int, end: int) -> bool:
    """Does [start, end) touch any span of `index[chrom]`? One base of overlap is enough."""
    pair = index.get(chrom)
    if not pair:
        return False
    starts, ends = pair
    i = bisect.bisect_right(starts, end - 1) - 1
    while i >= 0 and ends[i] > start:
        if starts[i] < end:
            return True
        i -= 1
    return False


def single_label_elements(
    ccres: list[CCRE],
    peaks: dict[str, dict[str, tuple[list[int], list[int]]]],
    labels: tuple[str, ...],
) -> tuple[dict[str, str], dict[str, int]]:
    """The fly's element rule, on mouse tissues: open in EXACTLY ONE label, by the source's own calls.

    `simpler_flysource.read_dars` keeps a region called up in exactly one of 42 cell types and drops
    the 18,129 called in more than one. Here the same rule runs over the accessibility source's own
    replicated-peak calls: a cCRE is kept when its span touches a peak in exactly one of `labels`,
    and dropped when it touches peaks in two or more, or in none. Nothing is thresholded here -- the
    peak lists are the source's.

    Returns the element id -> its one label, and the funnel that produced it.
    """
    label_of: dict[str, str] = {}
    funnel = {
        "ccres_on_the_chromosomes": len(ccres),
        "open_in_no_label": 0,
        "open_in_more_than_one_label_excluded": 0,
        "open_in_exactly_one_label": 0,
    }
    for c in ccres:
        hits = [lab for lab in labels if overlaps(peaks[lab], c.chrom, c.start, c.end)]
        if not hits:
            funnel["open_in_no_label"] += 1
        elif len(hits) > 1:
            funnel["open_in_more_than_one_label_excluded"] += 1
        else:
            funnel["open_in_exactly_one_label"] += 1
            label_of[c.eid] = hits[0]
    return label_of, funnel


def elements(
    ccres: list[CCRE],
    label_of: dict[str, str],
    chroms: dict[str, Chromosome],
) -> tuple[list[simpler.Element], dict[str, int]]:
    """Single-label cCREs as `simpler.Element`s, split promoter / distal by SCREEN's own class.

    `simpler_fly.elements`' funnel, with one difference registered as such: the promoter / distal
    split is SCREEN's class column, not a basepair distance chosen here, so `simpler_fly.PROMOTER_BP`
    and `DISTAL_BP` are NOT used and no element is dropped for sitting between the two.
    """
    out: list[simpler.Element] = []
    funnel = {
        "single_label_ccres": len(label_of),
        "promoter_class_positive_control": 0,
        "distal": 0,
        "no_protein_coding_gene_on_the_chromosome": 0,
    }
    for c in ccres:
        lab = label_of.get(c.eid)
        if lab is None:
            continue
        if c.chrom not in chroms or not chroms[c.chrom].genes:
            funnel["no_protein_coding_gene_on_the_chromosome"] += 1
            continue
        kind = c.kind
        funnel["promoter_class_positive_control" if kind == PROMOTER else "distal"] += 1
        out.append(simpler.Element(chrom=c.chrom, start=c.start, end=c.end, kind=kind, tissue=lab))
    out.sort(key=lambda e: (e.chrom, e.start, e.end))
    return out, funnel


# ------------------------------------------------------------------------------------- the boundaries


def read_bed_spans(path: Path, chromosomes: tuple[str, ...] = CHROMOSOMES) -> list[tuple[str, int, int]]:
    """(chrom, start, end) of a gzipped or plain BED, the first three columns only."""
    opener = gzip.open if str(path).endswith(".gz") else open
    rows: list[tuple[str, int, int]] = []
    with opener(path, "rt") as handle:  # type: ignore[operator]
        for line in handle:
            if line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] not in chromosomes:
                continue
            rows.append((parts[0], int(parts[1]), int(parts[2])))
    return rows


def read_atac_peaks(
    paths: dict[str, Path], chromosomes: tuple[str, ...] = CHROMOSOMES
) -> tuple[dict[str, dict[str, tuple[list[int], list[int]]]], dict[str, int]]:
    """Each label's own replicated-peak narrowPeak -> a merged, bisect-ready span index.

    Columns 1 to 3 only. No score, no q-value and no summit is read: the peak list IS the source's
    call, and reading its score would be choosing a threshold here.
    """
    index: dict[str, dict[str, tuple[list[int], list[int]]]] = {}
    counts: dict[str, int] = {}
    for label, path in paths.items():
        spans: dict[str, list[tuple[int, int]]] = {}
        n = 0
        for chrom, start, end in read_bed_spans(path, chromosomes):
            spans.setdefault(chrom, []).append((start, end))
            n += 1
        index[label] = peak_index(spans)
        counts[label] = n
    return index, counts


def read_expression(
    paths: dict[str, Path], labels: tuple[str, ...] = LABELS
) -> tuple[dict[str, tuple[float, ...]], dict[str, Any]]:
    """Each label's ENCODE gene quantification -> gene id (version stripped) -> TPM per label.

    The `gene_id` and `TPM` columns are found by HEADER NAME, so a column order change cannot be
    read as a value. Spike-in and ERCC rows, and any id that is not an ENSMUSG accession, are
    dropped. A gene absent from a label's rows is 0.0 there; a gene absent from EVERY label's rows
    has no entry at all and `concordant` returns False for it, which is how the fly's run counted
    its 3,911 missing slots and is not changed here.
    """
    per: dict[str, dict[str, float]] = {}
    rows_read: dict[str, int] = {}
    for label in labels:
        path = paths[label]
        opener = gzip.open if str(path).endswith(".gz") else open
        values: dict[str, float] = {}
        with opener(path, "rt") as handle:  # type: ignore[operator]
            header = handle.readline().rstrip("\n").split("\t")
            gi, ti = header.index("gene_id"), header.index("TPM")
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(gi, ti):
                    continue
                gid = parts[gi].split(".")[0]
                if not gid.startswith("ENSMUSG"):
                    continue
                values[gid] = max(values.get(gid, 0.0), float(parts[ti]))
        per[label] = values
        rows_read[label] = len(values)
    every = sorted(set().union(*(set(v) for v in per.values()))) if per else []
    expr = {g: tuple(per[label].get(g, 0.0) for label in labels) for g in every}
    summary = {
        "genes_with_a_row_in_at_least_one_label": len(expr),
        "rows_per_label": rows_read,
        "unit": "TPM, as ENCODE's RSEM gene quantification reports it",
        "detection": DETECTION_TPM,
    }
    return expr, summary


def boundaries_from_spans(
    rows: list[tuple[str, int, int]], chromosomes: tuple[str, ...] = CHROMOSOMES
) -> tuple[dict[str, list[int]], dict[str, Any]]:
    """`simpler_flysource.read_boundary_table` for the mouse: span midpoints, sorted, deduplicated."""
    bounds: dict[str, set[int]] = {}
    widths: list[int] = []
    for chrom, start, end in rows:
        if chrom not in chromosomes:
            continue
        bounds.setdefault(chrom, set()).add((start + end) // 2)
        widths.append(end - start)
    out = {c: sorted(v) for c, v in sorted(bounds.items())}
    total = sum(len(v) for v in out.values())
    spacing: list[int] = []
    for v in out.values():
        spacing.extend(b - a for a, b in zip(v, v[1:], strict=False))
    summary = {
        "boundaries": total,
        "chromosomes": len(out),
        "median_span_bp": sorted(widths)[len(widths) // 2] if widths else None,
        "median_spacing_bp": sorted(spacing)[len(spacing) // 2] if spacing else None,
        "mean_spacing_bp": round(sum(spacing) / len(spacing), 1) if spacing else None,
    }
    return out, summary


def read_domain_sheet(
    path: Path, sheet: str, expect_rows: int, chromosomes: tuple[str, ...] = CHROMOSOMES
) -> list[tuple[str, int, int]]:
    """One sheet of Bonev's Table S2 -> (chrom, start, end) domains. Refuses a changed workbook.

    The sheet's header is `chrom | start | end | Compartment`; the `Compartment` column (A or B) is
    NEVER read, because which compartment a domain sits in is not what the rule uses. `expect_rows`
    is the count recorded in BONEV_SHEETS before any RNA value was read: a sheet that does not
    reproduce it raises, so a reissued workbook cannot pass as the registered one.
    """
    import openpyxl

    book = openpyxl.load_workbook(path, read_only=True, data_only=True)
    if sheet not in book.sheetnames:
        raise ValueError(f"{path} has no sheet {sheet!r}; it has {book.sheetnames}")
    rows: list[tuple[str, int, int]] = []
    seen = 0
    for chrom, start, end, *_rest in book[sheet].iter_rows(values_only=True):
        if chrom == "chrom":
            continue
        if chrom is None or start is None or end is None:
            continue
        seen += 1
        if str(chrom) in chromosomes:
            rows.append((str(chrom), int(start), int(end)))
    if seen != expect_rows:
        raise ValueError(f"{path}:{sheet} holds {seen} domains, the registration recorded {expect_rows}")
    return rows


def boundaries_from_domains(
    domains: list[tuple[str, int, int]], chromosomes: tuple[str, ...] = CHROMOSOMES
) -> tuple[dict[str, list[int]], dict[str, Any]]:
    """Author-called DOMAINS -> boundary positions: every domain start and end, deduplicated, sorted.

    The only conversion this lane makes to a boundary set, and it has no parameter. A domain's two
    edges ARE its boundaries; where two domains abut, the shared edge is one boundary and not two.
    """
    bounds: dict[str, set[int]] = {}
    sizes: list[int] = []
    for chrom, start, end in domains:
        if chrom not in chromosomes:
            continue
        bounds.setdefault(chrom, set()).update((start, end))
        sizes.append(end - start)
    out = {c: sorted(v) for c, v in sorted(bounds.items())}
    spacing: list[int] = []
    for v in out.values():
        spacing.extend(b - a for a, b in zip(v, v[1:], strict=False))
    summary = {
        "domains": len(sizes),
        "boundaries": sum(len(v) for v in out.values()),
        "chromosomes": len(out),
        "median_domain_bp": sorted(sizes)[len(sizes) // 2] if sizes else None,
        "median_spacing_bp": sorted(spacing)[len(spacing) // 2] if spacing else None,
        "mean_spacing_bp": round(sum(spacing) / len(spacing), 1) if spacing else None,
    }
    return out, summary


def domain_of(bounds: list[int], pos: int, half_bin: int) -> tuple[float, float]:
    """`simpler_flysource.domain_of`, unchanged. Kept as a name so the tests can pin the identity."""
    return fs.domain_of(bounds, pos, half_bin)


def boundary_rule(chrom: Chromosome, bounds: list[int], pos: int, half_bin: int) -> tuple[Gene | None, str]:
    """`simpler_flysource.boundary_rule`, unchanged: the node is the fly's, not a new one."""
    return fs.boundary_rule(chrom, bounds, pos, half_bin)


# ------------------------------------------------------------------ attribution and position counts


def attribute(
    rows: list[dict[str, Any]],
    chroms: dict[str, Chromosome],
    bounds: dict[str, list[int]],
    half_bin: int,
) -> list[dict[str, Any]]:
    """`simpler_flysource.attribute`, unchanged."""
    return fs.attribute(rows, chroms, bounds, half_bin)


def distal(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [r for r in rows if r["kind"] == DISTAL]


def between_genes(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """THE POPULATION: distal elements whose midpoint lies in no protein-coding gene body.

    The stratum the worm and the fly disagree about, and the one where neither rule can win on a
    host gene's own transcription, because neither can name a host.
    """
    return [r for r in distal(rows) if r["domain_case"] != fn.HOST]


def mde_points(discordant: int, population: int) -> float | None:
    return fs.mde_points(discordant, population)


def _population(rows: list[dict[str, Any]], labels: tuple[str, ...]) -> dict[str, Any]:
    """`simpler_flysource._population` over this lane's labels."""
    n = len(rows)
    moves = [r for r in rows if r["boundary_gene"] != r["nearest_gene"]]
    k_gs = sum(1 for r in rows if r["boundary_gene"] != r["domain_gene"])
    return {
        "elements": n,
        "elements_by_label": {
            lab: sum(1 for r in rows if r["tissue"] == lab)
            for lab in labels
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
    funnel: dict[str, Any],
    labels: tuple[str, ...],
) -> dict[str, Any]:
    """Counts and power in advance, from POSITIONS, the element-side calls and the boundaries only.

    No RNA value is read by anything this function reaches. `stop_rule_met` is True when the
    population's between-gene elements that change gene reach `STOP_RULE`; the run stops unscored
    when it is False, and the inventory is then the result (`NO_SOURCE_QUALIFIES`).
    """
    out = {
        "population_funnel": funnel,
        "between_genes": _population(between_genes(main), labels),
        "all_distal": _population(distal(main), labels),
        "alternative_boundary_sets_between_genes": {
            name: _population(between_genes(rows), labels) for name, rows in sorted(alternatives.items())
        },
        "positive_control_promoters": sum(1 for r in main if r["kind"] == PROMOTER),
        "labels": len(labels),
    }
    out["between_gene_elements_that_change_gene"] = out["between_genes"]["differ_from_nearest_tss"]
    out["stop_rule"] = STOP_RULE
    out["stop_rule_met"] = out["between_gene_elements_that_change_gene"] >= STOP_RULE
    out["verdict"] = "DECIDABLE" if out["stop_rule_met"] else "NO_SOURCE_QUALIFIES"
    return out


# ------------------------------------------------------------------------------------- the statistics


def concordant(tpm: tuple[float, ...] | None, label: str, labels: tuple[str, ...]) -> bool:
    """`simpler.concordant`: at least DETECTION_TPM in the element's own label, and highest there."""
    if tpm is None:
        return False
    i = labels.index(label)
    v = tpm[i]
    return v >= DETECTION_TPM and all(v > tpm[j] for j in range(len(labels)) if j != i)


class ControlIndex:
    """`simpler_flysource.ControlIndex` over this lane's labels: the same matched-distance draw.

    A matched-control gene is drawn uniformly from the protein-coding genes of the element's
    chromosome whose 5' end lies in the same distance bin as the attributed gene's, the attributed
    gene excluded; the Bernoulli probability is computed exactly from prefix counts.
    """

    def __init__(
        self,
        chroms: dict[str, Chromosome],
        expr: dict[str, tuple[float, ...]],
        labels: tuple[str, ...],
    ):
        self.chroms = chroms
        self.labels = labels
        self.expr = expr
        self.prefix: dict[tuple[str, str], list[int]] = {}
        for c, ch in chroms.items():
            for lab in labels:
                acc, p = 0, [0]
                for g in ch.by_tss:
                    acc += concordant(expr.get(g.gid), lab, labels)
                    p.append(acc)
                self.prefix[(c, lab)] = p

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
        total, conc = n1 + n2, c1 + c2
        own = concordant(self.expr.get(gene), label, self.labels)
        total, conc = total - 1, conc - int(own)
        if total > 0:
            return conc / total
        ch = self.chroms[chrom]
        others = [g for g in ch.genes if g.gid != gene]
        if not others:
            return 0.0
        sub = min(others, key=lambda g: (abs(abs(g.tss - mid) - dist), g.gid))
        return float(concordant(self.expr.get(sub.gid), label, self.labels))


def compare(
    rows: list[dict[str, Any]],
    expr: dict[str, tuple[float, ...]],
    a: str,
    b: str,
    labels: tuple[str, ...],
) -> dict[str, Any]:
    """`simpler_flysource.compare`: rule `a` against rule `b`, percentile bootstrap over 1-Mb clusters."""
    for r in rows:
        r[f"c_{a}"] = int(concordant(expr.get(r[f"{a}_gene"]), r["tissue"], labels))
        r[f"c_{b}"] = int(concordant(expr.get(r[f"{b}_gene"]), r["tissue"], labels))
        r[f"d_{a}_{b}"] = r[f"c_{a}"] - r[f"c_{b}"]
    n = len(rows)
    ka = sum(r[f"c_{a}"] for r in rows)
    kb = sum(r[f"c_{b}"] for r in rows)
    bs = simpler.cluster_bootstrap(rows, f"d_{a}_{b}") if rows else None
    differ = sum(1 for r in rows if r[f"{a}_gene"] != r[f"{b}_gene"])
    ci = [round(100 * bs["lo"], 3), round(100 * bs["hi"], 3)] if bs else None
    return {
        "n": n,
        "differ": differ,
        f"{a}_concordant": ka,
        f"{b}_concordant": kb,
        f"{a}_rate_pct": round(100 * ka / n, 3) if n else None,
        f"{b}_rate_pct": round(100 * kb / n, 3) if n else None,
        "delta_points": round(100 * bs["point"], 3) if bs else None,
        "delta_ci95_points": ci,
        "clusters": bs["clusters"] if bs else 0,
        "mde_points_unclustered": mde_points(differ, n),
        "label": fn.label(ci, b),
    }


def rotation_null(
    rows: list[dict[str, Any]],
    chroms: dict[str, Chromosome],
    bounds: dict[str, list[int]],
    lengths: dict[str, int],
    expr: dict[str, tuple[float, ...]],
    labels: tuple[str, ...],
    half_bin: int,
    draws: int = ROTATION_DRAWS,
    seed: int = SEED,
) -> dict[str, Any]:
    """The registered rotated-boundary control: the same boundaries turned along each chromosome.

    `simpler_flysource.rotation_null`'s design, carried verbatim -- the NUMBER and the SPACING of
    the boundaries are kept and only their measured positions are destroyed, by one uniform offset
    per chromosome per draw. It asks whether the rule uses where the boundaries were measured or
    only how many there are and how far apart.
    """
    observed = sum(concordant(expr.get(r["boundary_gene"]), r["tissue"], labels) for r in rows)
    rnd = random.Random(seed)
    sims: list[int] = []
    for _ in range(draws):
        shifted = {
            c: sorted((b + rnd.randrange(lengths.get(c, 1))) % max(lengths.get(c, 1), 1) for b in v)
            for c, v in bounds.items()
        }
        k = 0
        for r in rows:
            gene, _case = boundary_rule(chroms[r["chrom"]], shifted.get(r["chrom"], []), r["mid"], half_bin)
            if gene is not None and concordant(expr.get(gene.gid), r["tissue"], labels):
                k += 1
        sims.append(k)
    n = len(rows)
    ge = sum(1 for s in sims if s >= observed)
    return {
        "n": n,
        "observed_concordant": observed,
        "observed_rate_pct": round(100 * observed / n, 3) if n else None,
        "rotated_mean_concordant": round(sum(sims) / len(sims), 2) if sims else None,
        "rotated_mean_rate_pct": round(100 * sum(sims) / len(sims) / n, 3) if sims and n else None,
        "draws": draws,
        "draws_at_or_above_observed": ge,
        "p_one_sided": round((1 + ge) / (1 + draws), 5),
    }


def positive_control(
    rows: list[dict[str, Any]],
    chroms: dict[str, Chromosome],
    expr: dict[str, tuple[float, ...]],
    labels: tuple[str, ...],
) -> dict[str, Any]:
    """`simpler_flysource.positive_control`: nearest TSS on the promoter-class elements, matched distance.

    It answers whether the readout is responsive at all. If the nearest-TSS gene of a
    promoter-proximal element is not expressed in that element's own tissue more often than a gene
    at a matched distance, nothing on the distal population could be believed.
    """
    index = ControlIndex(chroms, expr, labels)
    probs, observed = [], 0
    for r in rows:
        observed += int(concordant(expr.get(r["nearest_gene"]), r["tissue"], labels))
        probs.append(
            index.probability(r["chrom"], r["tissue"], r["mid"], r["nearest_gene"], r["nearest_dist"])
        )
    test = simpler.control_test(probs, observed)
    n = len(rows)
    return {
        "n": n,
        "concordant": observed,
        "rate_pct": round(100 * observed / n, 3) if n else None,
        **test,
    }


def matched_distance_control(
    rows: list[dict[str, Any]],
    chroms: dict[str, Chromosome],
    expr: dict[str, tuple[float, ...]],
    labels: tuple[str, ...],
    rule: str = "boundary",
) -> dict[str, Any]:
    """THE FALSIFIER, carried verbatim: the rule's genes against genes at a matched distance.

    If the attributed gene is not expressed in the element's own tissue more often than a
    protein-coding gene drawn at the same distance from the element, then position carries no
    tissue information on this population and the rule is falsified on this endpoint.
    """
    index = ControlIndex(chroms, expr, labels)
    probs, observed = [], 0
    for r in rows:
        observed += int(concordant(expr.get(r[f"{rule}_gene"]), r["tissue"], labels))
        probs.append(
            index.probability(r["chrom"], r["tissue"], r["mid"], r[f"{rule}_gene"], r[f"{rule}_dist"])
        )
    test = simpler.control_test(probs, observed)
    n = len(rows)
    return {
        "rule": rule,
        "n": n,
        "concordant": observed,
        "rate_pct": round(100 * observed / n, 3) if n else None,
        **test,
        "falsifier": "FALSIFIED" if test["p_one_sided"] > ALPHA else "NOT_FALSIFIED",
    }


# ----------------------------------------------------------------------------------------- the reading


def _pts(x: float | None) -> str:
    return "n/a" if x is None else f"{x:+.2f}"


def _rate(k: int, n: int) -> str:
    return f"{k:,} of {n:,} ({100 * k / n:.1f}%)" if n else "n/a"


def _verdict_sentence(ci: list[float] | None) -> str:
    """The three registered outcome words, as a sentence. `simpler.comparison_label` decides."""
    name = "NO_DIFFERENCE_DETECTED" if ci is None else simpler.comparison_label(*ci)
    return {
        "DOMAIN_RULE_AHEAD": "The boundary rule is ahead of nearest TSS",
        "NEAREST_AHEAD": "Nearest TSS is ahead of the boundary rule",
        "NO_DIFFERENCE_DETECTED": ("No difference detected between the boundary rule and nearest TSS"),
    }[name]


def reading(res: dict[str, Any]) -> str:
    """The registered reading, filled from the run and printed verbatim whatever it says.

    The template is the worm's and the fly's, sentence for sentence (`595d7c2`, `f5985f7`,
    `ff063e6`), with the organism, the source names and the two invertebrates' figures in it. It is
    written here, before any RNA value is read, and `score` prints it filled without editing it.
    """
    bg = res["between_genes"]["boundary_vs_nearest"]
    gs = res["between_genes"]["boundary_vs_gene_structure"]
    ad = res["all_distal"]["boundary_vs_nearest"]
    ags = res["all_distal"]["boundary_vs_gene_structure"]
    fal = res["between_genes"]["falsifier"]
    rot = res["between_genes"]["rotation"]
    pc = res["positive_control"]
    cfg = res["configuration"]
    return (
        f"On {bg['n']:,} Mus musculus distal candidate cis-regulatory elements "
        f"(ENCODE SCREEN mm10 registry v3) open in exactly one of {cfg['labels']} "
        f"{cfg['label_kind']} by the accessibility source's own replicated-peak calls "
        f"({cfg['element_source_short']}) and lying in no protein-coding gene body, the gene with "
        f"the nearest 5' end inside the element's domain, delimited by "
        f"{cfg['boundaries']:,} {cfg['boundary_source_short']}, is expressed at "
        f"{DETECTION_TPM:g} transcripts per million or more and highest in the element's own "
        f"{cfg['label_noun']}, in RNA from that same {cfg['label_noun']} "
        f"({cfg['rna_source_short']}), for {_rate(bg['boundary_concordant'], bg['n'])}, against "
        f"{_rate(bg['nearest_concordant'], bg['n'])} for the gene with the nearest 5' end on the "
        f"chromosome; the rules differ on {bg['differ']:,} of {bg['n']:,} elements; difference "
        f"{_pts(bg['delta_points'])} points, 95% {_pts(bg['delta_ci95_points'][0])} to "
        f"{_pts(bg['delta_ci95_points'][1])} over 1-Mb clusters. "
        f"{_verdict_sentence(bg['delta_ci95_points'])}. "
        f"With {bg['differ']:,} elements differing, a difference smaller than about "
        f"{bg['mde_points_unclustered']:.2f} points could go undetected. "
        f"Against the gene-structure rule on the same elements "
        f"({_rate(gs['domain_concordant'], gs['n'])}; the two differ on {gs['differ']:,}, and "
        f"neither names a host here), the difference is {_pts(gs['delta_points'])} points, 95% "
        f"{_pts(gs['delta_ci95_points'][0])} to {_pts(gs['delta_ci95_points'][1])}. "
        f"The boundary rule's genes match the element's {cfg['label_noun']} "
        f"{'more often than' if fal['falsifier'] == 'NOT_FALSIFIED' else 'NO more often than'} "
        f"genes at a matched distance ({fal['expected_rate'] * 100:.1f}%, one-sided p "
        f"{fal['p_one_sided']:.3f}), so position "
        f"{'carries' if fal['falsifier'] == 'NOT_FALSIFIED' else 'is NOT shown to carry'} "
        f"{cfg['label_noun']} information here. With the same boundaries rotated along each "
        f"chromosome the rule averages {rot['rotated_mean_rate_pct']:.1f}% on these elements "
        f"(one-sided p {rot['p_one_sided']:.3f}), so the measured positions "
        f"{'carry more than' if rot['p_one_sided'] <= ALPHA else 'are not shown to carry more than'} "
        f"the boundaries' number and spacing. "
        f"On all {ad['n']:,} distal elements (beside, deciding nothing; the rules differ on "
        f"{ad['differ']:,}) the boundary rule is concordant for "
        f"{_rate(ad['boundary_concordant'], ad['n'])} against "
        f"{_rate(ad['nearest_concordant'], ad['n'])} for nearest TSS, {_pts(ad['delta_points'])} "
        f"points [{_pts(ad['delta_ci95_points'][0])}, {_pts(ad['delta_ci95_points'][1])}], and "
        f"against {_rate(ags['domain_concordant'], ags['n'])} for the gene-structure rule, "
        f"{_pts(ags['delta_points'])} points [{_pts(ags['delta_ci95_points'][0])}, "
        f"{_pts(ags['delta_ci95_points'][1])}]. "
        f"The positive control separates: on {pc['n']:,} promoter-class elements the nearest-TSS "
        f"gene is concordant for {pc['rate_pct']:.1f}% against {pc['expected_rate'] * 100:.1f}% at "
        f"a matched distance (one-sided p {pc['p_one_sided']:.3f}). "
        f"{res['attributed_genes_absent_from_the_rna_rows']:,} of "
        f"{res['attributed_gene_slots']:,} attributed gene slots on this population name a gene "
        f"absent from the RNA source's rows and are counted as not concordant. "
        f"The mammal is asked here because the two invertebrates disagree about this very stratum: "
        f"the worm's domain rule is ahead of nearest TSS between genes by "
        f"+{WORM['between_gene_delta_points']:.2f} points "
        f"[+{WORM['between_gene_ci95_points'][0]:.2f}, +{WORM['between_gene_ci95_points'][1]:.2f}] "
        f"on {WORM['between_gene_elements']:,} elements, surviving its outron exclusion at "
        f"+{WORM['outron_excluded_delta_points']:.2f} "
        f"[+{WORM['outron_excluded_ci95_points'][0]:.2f}, "
        f"+{WORM['outron_excluded_ci95_points'][1]:.2f}] on "
        f"{WORM['outron_excluded_elements']:,}, while the fly's is "
        f"{FLY['brain_source_delta_points']:+.2f} "
        f"[{FLY['brain_source_ci95_points'][0]:+.2f}, {FLY['brain_source_ci95_points'][1]:+.2f}] "
        f"on {FLY['brain_source_elements']:,} brain elements, no difference detected. "
        f"Accessibility is not function, and a boundary is not a regulatory unit by itself: this "
        f"says the attributed gene is expressed where the element is open, not that the element "
        f"acts on it. "
        f"{cfg['what_is_worse_here']}"
    )


def load(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text())


# --------------------------------------------------------------------------------- the registration

REGISTRATION: dict[str, Any] = {
    "question": (
        "The domain rule attributes a non-coding element to a gene. Between genes, where no rule can "
        "name the gene whose body the element sits in, does the gene with the nearest 5' end inside "
        "the element's boundary-delimited domain match the element's own tissue more often than the "
        "gene with the nearest 5' end on the chromosome? Two invertebrates DISAGREE about exactly "
        "this stratum -- the worm's domain rule is ahead by +3.73 points [+2.79, +4.74] on 2,335 "
        "between-gene elements and still ahead at +4.47 [+3.16, +5.78] after its own outron "
        "exclusion, while the fly's is +0.12 [-0.25, +0.61] on 3,265 brain elements, no difference "
        "detected -- and NOBODY HAS ASKED A MAMMAL. This is that question in the mouse."
    ),
    "why_a_mammal": (
        "the fly's lead appears only where a host gene CAN be named (on all 10,052 distal elements "
        "the gene-structure rule is ahead of the boundary rule by -2.25 [-2.88, -1.66]), which is "
        "consistent with the fly's lead coming from host-gene transcription rather than regulation. "
        "The worm's lead survives where no host can be named. A mammal has 100 times the worm's "
        "intergenic space per gene, so whichever pattern it shows is the informative one."
    ),
    "stop_list_of_2026_09_27": (
        "read, not ignored. 'No more biosamples: 1.1 needs its readings to mean what their names "
        "say' has as its subject the HUMAN reader's 13 biosamples, READER_TERMS in "
        "genomeos/genome/context_evidence.py. This lane adds none, names none and widens none. A "
        "MOUSE TISSUE SOURCE IS NOT A HUMAN READER BIOSAMPLE, as lane-flysource recorded for its "
        "fly source (docs/SIMPLER-GENOMES.md section 11)."
    ),
    "assembly": ASSEMBLY,
    "liftover": LIFTOVER,
    "gene_models": GENE_MODELS,
    "gene_models_licence": GENE_MODELS_LICENCE,
    "element_registry": ELEMENT_SOURCE,
    "element_registry_licence": ELEMENT_LICENCE,
    "accessibility_source": ACCESSIBILITY_SOURCE,
    "rna_source": RNA_SOURCE,
    "encode_licence": ENCODE_LICENCE,
    "one_source_for_both_modalities": (
        "ENCODE, one stage (embryonic day 14.5) and the SAME TWELVE TISSUE NAMES on both sides. The "
        "vocabulary is the source's own: the twelve are exactly the biosample term names that have "
        "both an ATAC-seq and a polyA plus RNA-seq experiment at e14.5, and no tissue is grouped, "
        "renamed or split here. ENCODE's e14.5 `brain` experiment is excluded because it contains "
        "forebrain, midbrain and hindbrain, which are three of the twelve. The accessibility and the "
        "RNA come from different animals of the same strain and stage, not the same nuclei, and that "
        "is stated in the reading."
    ),
    "boundary_set": {
        "primary": BOUNDARY_SOURCE,
        "primary_url": BOUNDARY_URL,
        "primary_licence": BOUNDARY_LICENCE,
        "primary_method": BOUNDARY_METHOD,
        "what_a_boundary_is_here": BOUNDARY_CAVEAT,
        "why_this_one": (
            "it is the only author-called mouse Hi-C domain set found that is distributed in mm10 "
            "(so no liftover), states a licence in prose (CC BY 4.0), and was measured in EMBRYONIC "
            "tissue rather than a cell line. Of its five sheets, `ncx NPC` is the in-vivo "
            "progenitor population of the developing neocortex, the closest of the five to this "
            "panel's stage and to three of its twelve tissues. Section 11 registered its stage gap "
            "-- 2-3 h fly embryos against larval-to-adult brain elements -- as the price of a large "
            "enough population, and named a tissue-matched registration as the first thing that "
            "should follow; here the boundaries and the elements are both embryonic mouse."
        ),
        "what_was_searched_and_rejected": [
            "Gorkin et al. 2020 (Nature 583:744), the ENCODE mouse developmental matrix this "
            "panel's accessibility and RNA come from, publishes NO TAD or boundary table and "
            "generated no Hi-C: its own TADs are borrowed from Dixon et al. 2012 mESC. Its `d-TAC` "
            "is developmental transposase-accessible chromatin, not a TAD",
            "Dixon et al. 2012 (Nature 485:376, GSE35156), 2,200 mESC domains at 40-kb bins, is "
            "mm9 and would need a liftover; its article states no licence, and its supplementary "
            "legends could not be read from the 6.3 MB SI PDF, so which table holds the domains is "
            "inference rather than a read statement. Not used",
            "Vietri Rudan et al. 2015 (Cell Rep 10:1297, GSE65126) is mm10 adult liver but deposits "
            "only 50-kb interaction matrices; its Table S1 is 'Genomic Regions of Interest' and not "
            "a domain list. No author-called domain table exists",
            "Nora et al. 2017 (Cell 169:930, GSE98671) deposits a raw tar and 5C primers; no "
            "boundary or domain table, and no CC licence is stated",
            "Rao et al. 2014 (Cell 159:1665, GSE63525) does publish 2,927 Arrowhead contact domains "
            "for CH12-LX, but that is a B-lymphoblast CELL LINE and the coordinates are mm9",
            "ENCODE serves no TAD or boundary file for any organism (0 files of output type "
            "'topologically associated domains'). Its 105 mm10 `contact domains` bedpe files are "
            "Arrowhead calls on a 5-kb grid with a median length near 130 kb -- sub-TAD scale -- "
            "and ENCODE has NO embryonic or fetal mouse Hi-C at all, so none is tissue-matched to "
            "this panel",
            "4D Nucleome has 110 mm10 boundary BEDs, but they are DCIC-pipeline calls and not "
            "author calls, the portal requires a login to download, and the portal states no "
            "licence anywhere",
        ],
        "beside_deciding_nothing": [
            "Bonev's other four sheets on the same elements: `ncx CN` (2,307 domains, cortical "
            "neurons of embryonic neocortex, in vivo), `ES` (2,712, in vitro), `NPC` (2,522, in "
            "vitro) and `CN` (2,442, in vitro). None is dropped and none is chosen on an outcome",
            "ENCODE CTCF ChIP-seq, embryonic day 14.5 -- the panel's own stage -- conservative IDR "
            "thresholded peaks on mm10: LIVER (ENCSR397RHW, ENCFF123IQA) and LUNG (ENCSR677SIH, "
            "ENCFF558ZIE), both of which are among the panel's twelve tissues, and BRAIN "
            "(ENCSR000CEH, ENCFF756OJJ), a tissue that is the union of three of them. These three "
            "are EVERY CTCF ChIP-seq experiment ENCODE has released at e14.5",
        ],
        "bin": (
            f"{BIN_BP} bp ({HALF_BIN} bp each side), carried verbatim from simpler_flynode and NOT "
            "retuned for the mouse: a gene's 5' end inside either flanking boundary's bin counts as "
            "inside the domain"
        ),
    },
    "node": (
        "simpler_flysource.boundary_rule, unchanged: an element's domain runs from the boundary at or "
        "left of its midpoint to the next boundary right of it, or to the chromosome end; a gene is "
        "in the domain when its 5' end lies there or inside either boundary's bin; the rule names the "
        "in-domain gene with the nearest 5' end, and falls back to nearest TSS when the domain holds "
        "no 5' end. No new node is proposed here."
    ),
    "population": (
        "THE BETWEEN-GENE ELEMENTS: ENCODE SCREEN mm10 cCREs whose SCREEN class is not PLS, pELS or "
        "DNase-H3K4me3, which overlap an ATAC-seq replicated peak in EXACTLY ONE of the twelve "
        "labels, and whose midpoint lies in no protein-coding gene body. On none of them can either "
        "rule name a host, which is the whole reason for this population."
    ),
    "labels": LABELS,
    "endpoint": (
        "the worm's, unchanged: the attributed gene is counted concordant when it reaches 5 "
        "transcripts per million in the element's own tissue and is higher there than in every other "
        "tissue. A gene absent from the RNA source's rows is not concordant, and the number of such "
        "genes is reported. Accessibility is not function; this endpoint says the attributed gene is "
        "expressed where the element is open, never that the element regulates it."
    ),
    "outcome": (
        "the boundary rule against nearest TSS on the between-gene elements; percentile cluster "
        "bootstrap over 1-Mb windows; BOUNDARY_RULE_AHEAD, NEAREST_AHEAD or NO_DIFFERENCE_DETECTED"
    ),
    "second_registered_comparison": (
        "the boundary rule against the worm's gene-structure rule (simpler.domain_rule) on the same "
        "elements, with its own label, not the outcome. Two comparisons on one population, each at "
        "95%, no correction for the pair."
    ),
    "matched_control_and_falsifier": (
        "the matched-distance control applied to the boundary rule's gene on the population: "
        "FALSIFIED when the rule's genes match the element's tissue no more often than genes drawn "
        "at a matched distance from the element (one-sided p >= 0.05)"
    ),
    "positive_control": (
        "nearest TSS on the promoter-class elements (SCREEN's PLS, pELS and DNase-H3K4me3) against "
        "its own matched control. If it does not separate, the readout is declared unresponsive and "
        "the primary reading is withheld."
    ),
    "beside_deciding_nothing": [
        "the same boundaries rotated along each chromosome by one uniform offset, number and spacing "
        "kept, 1,000 draws",
        "the two other e14.5 CTCF sets on the same elements",
        "the elements where the two rules differ, alone",
        "all distal elements, and the gene-structure rule against nearest TSS",
        "labels holding 50 or more population elements",
    ],
    "power_in_advance": (
        "filled by scripts/simpler_mouse_run.py register from the COMMITTED inventory, which counts "
        "it from positions, the ATAC peak calls and the boundaries alone"
    ),
    "stop_rule": (
        f"one run; fewer than {STOP_RULE} between-gene elements changing gene under the boundary rule "
        "stops it unscored and the inventory becomes the result (NO_SOURCE_QUALIFIES). The worm's and "
        "the fly's rule, unchanged."
    ),
    "what_is_worse_here": (
        "The accessibility and the RNA come from different animals of the same strain and stage, "
        "not the same nuclei; a tissue is a mixture of cell types, so an element open in one tissue "
        "and a gene highest in that tissue can belong to different cells of it; the boundaries come "
        "from neural progenitors of embryonic neocortex and are applied to all twelve tissues, so "
        "they are stage-matched and tissue-matched to only three of them; and a domain's edges are "
        "read as its boundaries, which is the only conversion made to the boundary set."
    ),
    "money": "none; no model request, no key read, no AlphaGenome and no deletion request",
    "outcome_wording": {
        "BOUNDARY_RULE_AHEAD": "the boundary-delimited domain picks the expressed gene better",
        "NEAREST_AHEAD": "the nearest 5' end on the chromosome picks the expressed gene better",
        "NO_DIFFERENCE_DETECTED": "no difference detected between the two rules on this endpoint",
    },
}
