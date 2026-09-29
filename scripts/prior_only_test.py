# SPDX-License-Identifier: AGPL-3.0-or-later
"""The coherence pilot's prior-only labelling, tested on its own (item 13 follow-up, lane-prior).

Gate 2 of the coherence pilot (scripts/pilot_biological_gate.py, result 197c560) was discontinued by its
stop rule and left one descriptive lead: its prior-only labelling, which reads no holdable source (the
pilot's distance power law, the compiled target and H3K27ac), scored above distance to TSS on
Gasperini2019 (+0.193 [+0.134, +0.249]) and Xie (+0.086 [+0.022, +0.181]) on the pilot's seven
validation chromosomes. The pilot recorded its own caveat: the CRISPRi benchmark's held-out file, which
holds Morris and Xie, chose its positives by chromatin at the element, and the prior reads H3K27ac.

This script does three things, in the order they were committed:

1. `STUDIES`: for each CRISPRi study of the C4 harness, whether its tested elements or its positives were
   chosen using chromatin at the element, with the evidence (the benchmark's own documentation and each
   study's methods, quoted with their URLs) and the local counts that support it (`local_counts`).
   Committed before any score of any labelling of this script was computed.
2. The registration (`registration()`): the prior exactly as the pilot defined it, decomposed on
   identical pairs into distance alone and distance plus activity; the admissible endpoints; the
   baselines; the metric (C4's `holdout.compare`); the pass rule, the falsifier and each outcome's
   reading; what a pass would and would not mean; and what had already been seen. Written before any
   score.
3. The run, once:

       uv run python scripts/prior_only_test.py                 # writes data/results/prior_only_test.json
       uv run python scripts/prior_only_test.py --admissibility # step 1's local counts only; no score

The pilot's code is reused unchanged (genomeos/attribution/pilot.py, pilot_bio.py): each labelling is
`pilot_bio.solve(mode="prior")` over the ENCODE cCREs the held-out units overlap, read by
`pilot_bio.labels_from`; the three variants differ only in which of the prior's terms are present (no
compiled links, no H3K27ac peaks), never in a weight. An internal development result on withheld
sources, never a validation. No model request; the per-element response cache is not opened.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import resource
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import correctness as cx  # noqa: E402
from genomeos.attribution import holdout as ho  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.attribution import pilot as pl  # noqa: E402
from genomeos.attribution import pilot_bio as pb  # noqa: E402
from genomeos.results import save_result  # noqa: E402

NAME = "prior_only_test"
STATUS = pb.STATUS  # "internal development result on withheld sources"
ENDPOINT = ho.DECREASE_EP

# ==================================================================================================
# Step 1, admissibility (2026-09-29, lane-prior). Recorded before any score of any labelling here.
# ==================================================================================================
ADMISSIBILITY_DATE = "2026-09-29"
#: the marks the prior reads at an element: H3K27ac in its activity term (ENCODE replicated peaks), and
#: DNase only through the ENCODE cCRE registry that defines its blocks and so where it abstains
PRIOR_READS_MARKS = ("H3K27ac", "DNase (through the cCRE registry only)")
INADMISSIBLE, ADMISSIBLE_STATED, ADMISSIBLE = "inadmissible", "admissible_with_selection_stated", "admissible"
ADMISSIBILITY_RULE = (
    "an endpoint is inadmissible when its positives were selected or filtered using H3K27ac at the "
    "element, the mark the prior's activity term reads: the prior would then read the rule that made a "
    "pair positive, and a gain over distance could not be told from the selection. It is admissible, "
    "with the selection stated, when its tested elements (positives and negatives alike) were chosen "
    "using chromatin before the outcome was measured: a selection that does not depend on the outcome "
    "cannot make H3K27ac predict the outcome by construction, but it restricts the population, and the "
    "result then holds for chromatin-selected candidates only. It is admissible without qualification "
    "when neither holds. Only the decrease endpoint is considered: the prior states a link and no "
    "direction, so it is never asked about increases"
)
HELDOUT_FILE = "EPCrisprBenchmark_combined_data.heldout_5_cell_types.GRCh38.tsv.gz"
TRAINING_FILE = "EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz"
URL_NATURE = "https://www.nature.com/articles/s41586-026-10781-4"
URL_SI = (
    "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-026-10781-4/"
    "MediaObjects/41586_2026_10781_MOESM1_ESM.pdf"
)
URL_PAPER_REPO = (
    "https://github.com/EngreitzLab/ENCODE-rE2G-Paper/blob/main/crispr_analyses/"
    "CRISPR_comparison_resources/crispr_data/README.md"
)
URL_BENCHMARK_REPO = "https://github.com/EngreitzLab/CRISPR_comparison"
#: the benchmark's own filters, quoted (Gschwind et al., Nature 2026, Methods and Supplementary Methods)
BENCHMARK_EVIDENCE: dict[str, Any] = {
    "held_out_file": {
        "filters": "positives kept only with H3K27ac at the element and at least a 5% decrease; "
        "negatives kept by power only, never by chromatin",
        "quotes": [
            {
                "text": "To obtain a clean set of positives likely to represent real enhancer interactions, "
                "we applied filters based on effect size and chromatin state at the tested element.",
                "where": "Supplementary Methods 22.1",
                "url": URL_SI,
            },
            {
                "text": "we removed positives with elements in the CTCF element, H3K27me3 element, or no "
                "H3K27ac categories.",
                "where": "Supplementary Methods 22.1.4",
                "url": URL_SI,
            },
            {
                "text": "Positives lacking H3K27ac signal (CTCF, H3K27me3 or no H3K27ac categories) were "
                "removed to focus on enhancer-mediated regulation.",
                "where": "Methods, Creating the combined held-out CRISPR dataset",
                "url": URL_NATURE,
            },
            {
                "text": "Elements were categorized hierarchically based on H3K27ac peak overlap and RPM "
                "percentile as High H3K27ac, H3K27ac, CTCF, H3K27me3 or No H3K27ac",
                "where": "Methods",
                "url": URL_NATURE,
            },
            {
                "text": "EPCrisprBenchmark.combined_heldout.annotated.tsv.gz: Combined held-out dataset "
                "before filtering for non-enhancer elements (see Methods)",
                "where": "ENCODE-rE2G-Paper, crispr_data README",
                "url": URL_PAPER_REPO,
            },
        ],
        "reproduced_by_this_lane": "the unfiltered combined held-out file of the paper's repository "
        "(EPCrisprBenchmark.combined_heldout.annotated.tsv.gz, sha256 "
        "2390273590c5bc35090f78662f645e815249942bc41166af8f60651520b78ddd, 4,437 rows, read in a scratch "
        "copy and not added to the project) becomes the project's held-out file row for row, 4,378 of "
        "4,378 with no exception, under one rule: keep every negative, and keep a positive only if its "
        "elementChromatinCategory is H3K27ac or High H3K27ac and its EffectSize is at most -0.05. The 59 "
        "rows removed are all positives (25 outside the H3K27ac categories, 34 inside them with a smaller "
        "effect); none was relabelled a negative. Morris lost 1 and Xie 2, all on validation chromosomes",
    },
    "training_file": {
        "filters": "distance 1 kb to 1 Mb, promoters and target gene bodies excluded, negatives by power; "
        "no chromatin filter on positives or negatives",
        "quotes": [
            {
                "text": "Datasets were filtered to retain element–gene pairs with distances of 1 kb–1 Mb "
                "from TSS, excluding elements overlapping GENCODE v.29-annotated promoters or target gene "
                "bodies. Negative pairs required ≥80% power to detect a 15% decrease in expression.",
                "where": "Methods, Creating the combined K562 CRISPR training dataset",
                "url": URL_NATURE,
            },
            {
                "text": "Combined K562 training dataset with same filters applied as for the filtered "
                "held-out dataset (see Methods)",
                "where": "ENCODE-rE2G-Paper, crispr_data README: a separate filtered copy "
                "(EPCrisprBenchmark.combined_training.annotated.filtered_positives.tsv.gz), used for one "
                "comparison figure; the project holds the unfiltered training file",
                "url": URL_PAPER_REPO,
            },
        ],
        "local_check": "the project's training file keeps positives outside the H3K27ac categories "
        "(Gasperini2019: 3 CTCF element, 2 No H3K27ac; Nasser2021: 2 CTCF element), which the held-out "
        "file's rule would have removed",
    },
    "element_universe": {
        "text": "the benchmark's re-analysed elements are DNase-seq peaks: 'Candidate elements were defined "
        "by extending the summits of the top 150,000 K562 DNase-seq peaks' (Gasperini2019, Supplementary "
        "Methods 18.2) and 'we defined a unified set of candidate elements using the ABC pipeline 1 from "
        "K562 DNase-seq data' (Xie, Klann, Morris; Supplementary Methods 19.1). This touches positives and "
        "negatives alike",
        "url": URL_SI,
    },
}
#: one entry per CRISPRi source of the C4 harness (the benchmark's `Dataset`, both files); the class is
#: ADMISSIBILITY_RULE applied to the evidence quoted. Quotes are verbatim from the URL given
STUDIES: dict[str, dict[str, Any]] = {
    "crispri:Gasperini2019": {
        "file": TRAINING_FILE,
        "cells": ["K562"],
        "study": "Gasperini et al. 2019, Cell (CRISPRi with single-cell RNA-seq), re-analysed by the "
        "benchmark",
        "tested_elements": "chosen on chromatin, H3K27ac included: DNase peaks intersected with H3K27ac, "
        "EP300, GATA1 and RNA Pol II ChIP-seq (pilot); the top DNase peaks by a classifier over 170 ChIP-seq "
        "tracks (at scale); DNase peaks chosen for epigenomic diversity (948 exploratory). The benchmark's "
        "elements are the top 150,000 K562 DNase-seq peaks",
        "positives": "by the expression readout only; no chromatin filter in the training file",
        "quotes": [
            {
                "text": "The 1,119 candidate enhancers were all intergenic DNase I hypersensitive sites "
                "(DHSs) representing various combinations of H3K27 acetylation, p300, GATA1, and RNA Pol II "
                "binding (Figure 2A).",
                "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6690346/",
            },
            {
                "text": "The remaining regions were largely taken from intersections with K562 GATA1 "
                "ChIP-seq narrowPeaks (ENCSR000EFT, lifted to hg19), H3K27ac ChIP-seq narrowPeaks "
                "(ENCSR000AKP, lifted to hg19), RNA Pol II ChIP-seq narrowPeaks (ENCSR000AKY), and EP300 "
                "ChIP-seq narrowPeaks (ENCSR000EHI) (Figure 2A).",
                "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6690346/",
            },
            {
                "text": "A logistic regression classifier built using the 145 enhancer-gene pairs originally "
                "identified in the pilot experiment ... was used to select the top 5,000 intergenic open "
                "chromatin regions in K562s (as defined by DNase-seq narrowPeaks (ENCSR000EKS)).",
                "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6690346/",
            },
        ],
        "note": "ENCSR000AKP, the H3K27ac experiment Gasperini intersected, is the experiment whose "
        "replicated peaks (ENCFF532MMV) the prior reads for K562",
        "class": ADMISSIBLE_STATED,
    },
    "crispri:Schraivogel2020": {
        "file": TRAINING_FILE,
        "cells": ["K562"],
        "study": "Schraivogel et al. 2020, Nature Methods (TAP-seq), chr8 and chr11 only",
        "tested_elements": "chosen on chromatin: DNase hotspots overlapping GenoSTAN active-enhancer "
        "states, a chromatin-state model fitted on histone marks, p300 and DNase; whether H3K27ac is "
        "among its marks was not confirmed from a quote",
        "positives": "by the expression readout only; no chromatin filter in the training file",
        "quotes": [
            {
                "text": "DNase HS hotspots that overlapped with GenoSTAN active enhancer annotations "
                "(Enh.15, EnhWF.2, EnhF.10) by at least 50 bases were included as candidate enhancers.",
                "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7610614/",
            },
            {
                "text": "we perturbed all 1,778 putatively active enhancers predicted based on ENCODE data "
                "in two regions on chromosome 8 and 11",
                "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7610614/",
            },
        ],
        "note": "every unit lies on chr8 or chr11, both validation chromosomes, where the pilot already "
        "scored the prior-only labelling",
        "class": ADMISSIBLE_STATED,
    },
    "crispri:Nasser2021": {
        "file": f"{TRAINING_FILE} (K562) and {HELDOUT_FILE} (GM12878, Jurkat): one harness source",
        "cells": ["K562", "GM12878", "Jurkat"],
        "study": "Nasser et al. 2021, Nature: a compendium, mostly Fulco et al. 2019 (all K562 DHS within "
        "450 kb of 30 genes), with Fulco 2016 (tiling) and Klann 2017 (DHS); new FlowFISH in GM12878 and "
        "Jurkat",
        "tested_elements": "K562: DNase peaks (Fulco 2019, Klann 2017) or tiling (Fulco 2016); GM12878 and "
        "Jurkat: accessible regions (ATAC or DNase peaks) around PPIF, and ABC-selected elements for 12 "
        "genes",
        "positives": "in the training file by the readout only; in the held-out file filtered on H3K27ac "
        "by the benchmark",
        "quotes": [
            {
                "text": "tested all DNase I hypersensitive (DHS) elements in K562 cells within 450 kb of any "
                "of the genes",
                "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6886585/",
            },
            {
                "text": "we designed gRNAs tiling across all accessible regions (here, defined as the union "
                "of the MACS2 narrow peaks and 250-bp regions on either side of the MACS2 summit)",
                "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9153265/",
            },
        ],
        "note": "the harness pools both files under one source, and its decrease endpoint cannot be scored "
        "anywhere: 111 positives and 0 well-powered negatives",
        "class": INADMISSIBLE,
    },
    "crispri:HCT116": {
        "file": HELDOUT_FILE,
        "cells": ["HCT116"],
        "study": "Guckelberger et al. 2024, bioRxiv (CRUDO-FlowFISH, HCT116)",
        "tested_elements": "DNase peaks around five genes chosen with ABC predictions",
        "positives": "filtered on H3K27ac twice: the study's own enhancer call and the benchmark's filter",
        "quotes": [
            {
                "text": "We identified elements as enhancers if they fulfilled the following criteria: ... "
                "v) an H3K27ac level above the median of all tested elements.",
                "url": "https://www.biorxiv.org/content/10.1101/2024.07.12.603288v1.full",
            },
        ],
        "note": "34 positives and no well-powered negative: unscorable in the harness",
        "class": INADMISSIBLE,
    },
    "crispri:K562_DC_TAP": {
        "file": HELDOUT_FILE,
        "cells": ["K562"],
        "study": "Ray et al. 2025, bioRxiv (DC-TAP-seq)",
        "tested_elements": "DNase peaks chosen at random, of any signal strength",
        "positives": "filtered on H3K27ac by the benchmark",
        "quotes": [
            {
                "text": "Within these loci, we randomly selected ~40 DHS elements of any DHS signal strength "
                "to be included in the target set.",
                "url": "https://www.biorxiv.org/content/10.1101/2025.09.16.676677v1.full",
            },
        ],
        "note": "12 positives genome-wide: below the harness's floor",
        "class": INADMISSIBLE,
    },
    "crispri:WTC11_DC_TAP": {
        "file": HELDOUT_FILE,
        "cells": ["WTC11"],
        "study": "Ray et al. 2025, bioRxiv (DC-TAP-seq)",
        "tested_elements": "DNase peaks chosen at random, of any signal strength",
        "positives": "filtered on H3K27ac by the benchmark",
        "quotes": [
            {
                "text": "Our goal for this design was to select putative cis-regulatory elements with as "
                "little prior bias as possible.",
                "url": "https://www.biorxiv.org/content/10.1101/2025.09.16.676677v1.full",
            },
        ],
        "note": "15 decreases genome-wide: below the harness's floor",
        "class": INADMISSIBLE,
    },
    "crispri:Klann": {
        "file": HELDOUT_FILE,
        "cells": ["K562"],
        "study": "Klann et al. 2021, bioRxiv (wgCERES)",
        "tested_elements": "every K562 DNase peak genome-wide",
        "positives": "filtered on H3K27ac by the benchmark",
        "quotes": [
            {
                "text": "For the initial genome-wide library, 1,092,706 gRNAs were selected, targeting "
                "111,756 DHSs",
                "url": "https://www.biorxiv.org/content/10.1101/2021.03.08.434470v1.full",
            },
        ],
        "note": "21 positives and no well-powered negative: unscorable",
        "class": INADMISSIBLE,
    },
    "crispri:Morris": {
        "file": HELDOUT_FILE,
        "cells": ["K562"],
        "study": "Morris et al. 2023, Science (STING-seq)",
        "tested_elements": "fine-mapped GWAS variants inside K562 cCREs defined from DNase, H3K27ac and "
        "ATAC peaks (41 of 507 deliberately in closed chromatin)",
        "positives": "filtered on H3K27ac by the benchmark",
        "quotes": [
            {
                "text": "we used DNase I hypersensitive sites (DHS) from ENCODE (65), H3K27ac ChIP-seq peak "
                "calls from ENCODE, and ATAC-seq peak calls that we generated previously (81) to identify "
                "candidate cis-regulatory elements (cCREs).",
                "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10518238/",
            },
        ],
        "note": "one of the pilot's four primary endpoints; on the fresh partition 14 positives, below the "
        "floor even with its one removed positive restored",
        "class": INADMISSIBLE,
    },
    "crispri:Reilly": {
        "file": HELDOUT_FILE,
        "cells": ["K562"],
        "study": "Reilly et al. 2021, Nature Genetics (HCR-FlowFISH)",
        "tested_elements": "mostly tiling; every DHS at the GATA1 locus",
        "positives": "filtered on H3K27ac by the benchmark",
        "quotes": [
            {
                "text": "For each library, 52,500 guides were designed to target > 1 Mb surrounding each "
                "gene.",
                "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8925018/",
            },
        ],
        "note": "8 positives: unscorable",
        "class": INADMISSIBLE,
    },
    "crispri:Xie": {
        "file": HELDOUT_FILE,
        "cells": ["K562"],
        "study": "Xie et al. 2019, Cell Reports (Mosaic-seq)",
        "tested_elements": "DNase peaks with H3K4me1, about half chosen for strong p300; H3K27ac not used",
        "positives": "filtered on H3K27ac by the benchmark",
        "quotes": [
            {
                "text": "The putative enhancers targeted in this study are defined by DNase-seq peaks which "
                "are at least 2kb away from any annotated TSS and that also harbor H3K4me1 signal",
                "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6904118/",
            },
        ],
        "note": "one of the pilot's four primary endpoints and half of its lead (+0.086); on the fresh "
        "partition 16 positives, below the floor even with its two removed positives restored",
        "class": INADMISSIBLE,
    },
}
#: the held-out file's filter makes every study in it inadmissible, whatever its own design
HELDOUT_REASON = (
    "its positives are the benchmark's held-out positives, kept only with H3K27ac at the element "
    "(BENCHMARK_EVIDENCE['held_out_file']); the prior's activity term reads H3K27ac"
)
ADMISSIBLE_SOURCES = tuple(s for s, v in STUDIES.items() if v["class"] != INADMISSIBLE)


def benchmark_rows(knowledge: Path = ms.CRISPRI_KNOWLEDGE) -> dict[str, list[dict[str, str]]]:
    """The benchmark's valid rows, per harness source, both files."""
    out: dict[str, list[dict[str, str]]] = defaultdict(list)
    for name in ms.CRISPRI_FILES:
        p = knowledge / name
        if not p.exists():
            continue
        with gzip.open(p, "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("ValidConnection") == "TRUE":
                    out[f"crispri:{r.get('Dataset', '')}"].append({**r, "_file": name})
    return dict(out)


H3K27AC_CATEGORIES = ("H3K27ac", "High H3K27ac")


def local_counts(rows: dict[str, list[dict[str, str]]]) -> dict[str, Any]:
    """Per source and file: the benchmark's chromatin category of each positive (Regulated) and of every
    pair, and how many positives lie outside the H3K27ac categories. Read to establish selection; this
    counts the benchmark's own annotation and scores no labelling."""
    out: dict[str, Any] = {}
    for src, rs in sorted(rows.items()):
        per: dict[str, Any] = {}
        for f in sorted({r["_file"] for r in rs}):
            sub = [r for r in rs if r["_file"] == f]
            pos = [r for r in sub if r.get("Regulated") == "TRUE"]
            per[f] = {
                "pairs": len(sub),
                "positives": len(pos),
                "positives_by_category": dict(Counter(r["elementChromatinCategory"] for r in pos)),
                "pairs_by_category": dict(Counter(r["elementChromatinCategory"] for r in sub)),
                "positives_outside_h3k27ac_categories": sum(
                    r["elementChromatinCategory"] not in H3K27AC_CATEGORIES for r in pos
                ),
                "pairs_outside_h3k27ac_categories": sum(
                    r["elementChromatinCategory"] not in H3K27AC_CATEGORIES for r in sub
                ),
            }
        out[src] = per
    return out


def partition_counts(units: dict[str, tuple[ho.Unit, ...]], chroms: tuple[str, ...]) -> dict[str, Any]:
    """Units, positives, negatives and loci of each CRISPRi decrease endpoint on a partition: counts,
    no score."""
    out: dict[str, Any] = {}
    keep = set(chroms)
    for src in sorted(s for s in units if ho.assay_of(s) == "crispri"):
        kept = [u for u in units[src] if u.chrom in keep and ho.target(u, ENDPOINT) is not None]
        pos = int(sum(ho.target(u, ENDPOINT) or 0 for u in kept))
        n_loci = max(ho.loci(kept)) + 1 if kept else 0
        out[src] = {
            "units": len(kept),
            "positives": pos,
            "negatives": len(kept) - pos,
            "loci": n_loci,
            "clears_floors": pos >= ho.MIN_POSITIVES
            and len(kept) - pos >= ho.MIN_NEGATIVES
            and n_loci >= ho.MIN_LOCI,
        }
    return out


def admissibility() -> dict[str, Any]:
    """Step 1, as recorded: the rule, the benchmark's evidence, each study's class and the local counts."""
    rows = benchmark_rows()
    return {
        "date": ADMISSIBILITY_DATE,
        "rule": ADMISSIBILITY_RULE,
        "prior_reads_marks": list(PRIOR_READS_MARKS),
        "benchmark": BENCHMARK_EVIDENCE,
        "held_out_file_reason": HELDOUT_REASON,
        "studies": STUDIES,
        "admissible_sources": list(ADMISSIBLE_SOURCES),
        "local_counts": local_counts(rows),
    }


# ==================================================================================================
# Step 2, the registration (2026-09-29, lane-prior). Written before any score of any labelling here;
# docs/ATTRIBUTION.md carries the same text under its dated heading.
# ==================================================================================================
REGISTERED = "2026-09-29"
#: the fresh partition: the 17 chromosomes on which no labelling of the pilot, the prior-only one
#: included, was ever scored (the pilot was timed and debugged there without a held-out score)
FRESH = pb.DEVELOPMENT
#: the seen partition: the pilot's seven validation chromosomes, where the lead was read
SEEN = pb.VALIDATION
PARTITIONS = {
    "fresh": "the 17 development chromosomes of the pilot ("
    + ", ".join(FRESH)
    + "): no labelling of the pilot, the prior-only one included, was ever scored here. The pass rule "
    "is read here and nowhere else",
    "seen": "the pilot's 7 validation chromosomes ("
    + ", ".join(SEEN)
    + "): the lead was read here (197c560). Reported beside, for the reproduction check and the "
    "decomposition of the frozen lead, and carrying no weight",
}
#: the endpoints the pass rule reads: admissible, and clearing the harness's floors on the fresh partition
PRIMARY = (("crispri:Gasperini2019", ENDPOINT),)
#: admissible endpoints reported on the seen partition only
BESIDE_SEEN = (("crispri:Gasperini2019", ENDPOINT), ("crispri:Schraivogel2020", ENDPOINT))
#: the three prior labellings: the pilot's prior, and the same prior with terms removed. They share their
#: blocks, their units and their abstentions, so they differ by one term at a time on identical pairs
PRIOR_DISTANCE, PRIOR_ACTIVITY, PRIOR = "prior_distance", "prior_distance_activity", pb.PRIOR
DISTANCE, UNCHANGED = ho.DISTANCE, ho.UNCHANGED
LABELLINGS = (PRIOR, PRIOR_ACTIVITY, PRIOR_DISTANCE, DISTANCE, UNCHANGED)
#: which of the prior's inputs each variant keeps: compiled links (the compiled target, and its model
#: cell's activity bonus) and H3K27ac peaks. No weight of pilot.W is changed
VARIANTS = {
    PRIOR_DISTANCE: {"compiled_links": False, "h3k27ac_peaks": False},
    PRIOR_ACTIVITY: {"compiled_links": False, "h3k27ac_peaks": True},
    PRIOR: {"compiled_links": True, "h3k27ac_peaks": True},
}
VARIANT_READS = {
    PRIOR_DISTANCE: frozenset({"gencode_v50", "encode_ccre_registry"}),
    PRIOR_ACTIVITY: frozenset({"gencode_v50", "encode_ccre_registry", "encode_h3k27ac"}),
    PRIOR: pb.READS_ALWAYS,
}
LABELLING_TEXT = {
    PRIOR: "(c) the prior-only labelling exactly as the pilot defined it: pilot_bio.solve(mode='prior') on "
    "pilot_bio.Fixed.load, pilot.W unchanged. A pair's score is sigmoid(t + 3) x sigmoid(a), the largest "
    "over the ENCODE cCREs the unit overlaps, None outside every cCRE. t = -ln(1 + d / 5000) with d from "
    "the block's midpoint to the pair's TSS (the benchmark's), plus 1 + min(strength, 1) where the block's "
    "compiled link names the pair's gene; a = -1, plus 2 with a K562 ENCODE H3K27ac replicated peak over "
    "the block or -1 without (another cell: 2 x the share of 13 biosamples with a peak), plus 0.5 where "
    "the compiled link's model cell is the pair's cell",
    PRIOR_ACTIVITY: "(b) the activity-and-distance baseline: the same prior with no compiled links, so its "
    "score is sigmoid(t + 3) x sigmoid(a) with t the distance term alone and a the H3K27ac term alone. "
    "H3K27ac only: the pilot's prior reads no DNase except through the cCRE registry that defines its "
    "blocks. It is not ABC (see ABC_FORM)",
    PRIOR_DISTANCE: "(a) distance alone, within the pilot's prior: no compiled links and no H3K27ac peaks, "
    "so a = -1 for every block and the score, sigmoid(-ln(1 + d / 5000) + 3) x sigmoid(-1), is monotone "
    "in the block's distance to the TSS. Same blocks and abstentions as (b) and (c)",
    DISTANCE: "distance to TSS, C4's baseline: " + ho.DISTANCE_RULE,
    UNCHANGED: "the unchanged compiled labels, C4's baseline: " + ho.UNCHANGED_RULE,
}
ABC_FORM = (
    "not computed. ABC's score is activity x contact normalised over every candidate element within 5 Mb "
    "of the gene, with activity the geometric mean of DNase-seq and H3K27ac signal at the element. The "
    "project holds ENCODE H3K27ac peaks and signal and K562 Hi-C, but no experimental DNase-seq signal "
    "(only AlphaGenome's predicted DNase, a model output), so ABC's activity term cannot be formed, and an "
    "H3K27ac-only substitute would not be the ABC form. Labelling (b) is therefore called an "
    "activity-and-distance baseline and never ABC"
)
#: (a, b, what the paired difference asks); `holdout.compare` computes each
COMPARISONS = (
    (PRIOR, DISTANCE, "the frozen lead: the pass rule"),
    (PRIOR_ACTIVITY, PRIOR_DISTANCE, "activity on identical pairs"),
    (PRIOR, PRIOR_ACTIVITY, "the compiled target on identical pairs"),
    (
        PRIOR_DISTANCE,
        DISTANCE,
        "the pilot's machinery: its distance form over cCRE blocks and its abstentions",
    ),
    (PRIOR_ACTIVITY, DISTANCE, "the activity-and-distance baseline against distance to TSS"),
    (PRIOR, UNCHANGED, "the prior against the unchanged compiled labels"),
)
METRIC = (
    "C4's harness unchanged: average precision of the CRISPRi significant decrease (1) against the "
    "well-powered nulls (0), abstentions ranked below every stated score, coverage beside; each paired "
    "difference from holdout.compare (both labellings on the same units and the same 1,000 locus "
    "resamples, seed 20260928), with its 95% percentile interval"
)
USEFUL_COVERAGE = 0.80
PASS_RULE = (
    "the lead replicates if, on every primary endpoint (an admissible CRISPRi decrease endpoint that "
    "clears the harness's floors of 20 positives, 20 negatives and 10 loci on the fresh partition: "
    "Gasperini2019 alone), the prior-only labelling exactly as the pilot defined it has a paired "
    "difference in average precision against distance to TSS whose 95% locus-bootstrap interval lies "
    "above zero, at a coverage of at least 0.80 of the endpoint's units"
)
DECOMPOSITION_RULE = (
    "on identical pairs (the same units, blocks and abstentions): activity carries the gain when (b) minus "
    "(a) lies above zero; the compiled target adds to it when (c) minus (b) lies above zero; the pilot's "
    "machinery carries it when (a) minus distance to TSS lies above zero. Each is reported whatever the pass"
)
LEAD = {
    "source": "crispri:Gasperini2019",
    "endpoint": ENDPOINT,
    "a": PRIOR,
    "b": DISTANCE,
    "difference": 0.1927,
    "interval": [0.1335, 0.2491],
    "a_value": 0.6601,
    "b_value": 0.4674,
    "commit": "197c560",
    "partition": "seen",
}
FALSIFIER = (
    "a replication is not an activity-and-distance reading (i) if (b) minus (a) does not lie above zero: "
    "activity does not carry the gain; (ii) if (c) minus (b) lies above zero while (b) minus distance to "
    "TSS does not: the gain over distance is the compiled target's, an AlphaGenome-derived term, not "
    "activity's; (iii) if (a) minus distance to TSS lies above zero while (b) minus (a) does not: the gain "
    "is the machinery's (cCRE coverage and the distance form). The whole run is void, and nothing else is "
    "read, if the seen partition does not reproduce the registered lead to the fourth decimal "
    "(Gasperini2019 prior-only minus distance +0.1927 [+0.1335, +0.2491], 0.6601 against 0.4674, "
    "197c560): the prior, its inputs or the harness would then have changed since the lead was read. Any "
    "LeakError is fatal and never caught"
)
READINGS = {
    "replicated_activity": "an internal development result on withheld sources: on chromosomes where it "
    "was never scored, the pilot's frozen prior predicts Gasperini2019's CRISPRi decreases better than "
    "distance to TSS, and on identical pairs H3K27ac activity adds to distance. The project's prior "
    "reproduces the known activity-and-distance relationship on CRISPRi decreases whose positives were "
    "not chosen by chromatin, among candidates that were. Not a new finding; not a validation of the "
    "attribution layer",
    "replicated_not_activity": "the lead replicates where it was never scored, but a clause of the "
    "falsifier fires: activity is not what carries it, and no activity-and-distance claim may be made. "
    "What carries it is named from the decomposition",
    "not_replicated": "the lead does not replicate where it was never scored: it stays a descriptive "
    "observation of the validation chromosomes, and may not be used",
    "reversed": "the prior is worse than distance to TSS where it was never scored: the lead was specific "
    "to the chromosomes it was read on",
    "no_admissible_endpoint": "no admissible CRISPRi endpoint clears the harness's floors on the fresh "
    "partition: the lead cannot be tested with what the project holds and stays descriptive",
    "void_reproduction": "the seen partition did not reproduce the registered lead: the run is void and "
    "only the discrepancy is reported",
}
MEANING = {
    "a_pass_would_mean": "the project's prior, its weights registered 2026-09-29 before any score, "
    "reproduces the known activity-and-distance relationship (an element carrying H3K27ac near a gene's "
    "TSS is a CRISPRi hit for that gene more often than its distance alone predicts) on K562 CRISPRi "
    "decreases whose positives were not chosen by chromatin, on chromosomes where this prior was never "
    "scored. The tested candidates were chosen on chromatin (DNase peaks, and for Gasperini2019 also "
    "H3K27ac, p300, GATA1 and Pol II ChIP), so it holds for chromatin-selected candidates only",
    "a_pass_would_not_mean": [
        "a new finding: the relationship is published on these screens (Fulco et al. 2019; Nasser et al. "
        "2021; Gschwind et al. 2026), and this project read it on the pooled training file on 2026-09-16 "
        "(crispri_benchmark, 138824f: activity over distance 0.519 against distance 0.441 in average "
        "precision)",
        "a validation of the attribution layer, of the compiled labels, or of the coherence pilot, which "
        "stays a discontinued investigation",
        "a reproduction of ABC: (b) is an activity-and-distance baseline, not activity x contact "
        "normalised over candidates",
        "an external, independent or fresh validation: Gasperini2019 is development evidence that many "
        "earlier results read",
        "anything about unselected elements, about cells other than K562, or about the held-out file's "
        "studies, which are inadmissible here",
    ],
    "a_failure_would_mean": READINGS["not_replicated"],
}
SEEN_BEFORE = (
    {
        "what": "activity over distance against distance on the pooled training file (Gasperini2019, "
        "Nasser2021 and Schraivogel2020, K562): average precision 0.519 against 0.441, 9,237 pairs, "
        "451 positives; the same in crispri_contact",
        "result": "crispri_benchmark, crispri_contact",
        "commit": "138824f, 5a31c39 (2026-09-16)",
    },
    {
        "what": "the held-out file read by at least ten scored results (crispri_benchmark, crispri_contact, "
        "crispri_published, crispri_direction, crispri_direction_both, target_calibration and others)",
        "result": "crispri_split_audit",
        "commit": "0af7e1b (2026-09-28)",
    },
    {
        "what": "distance to TSS, the unchanged labels and rest on every CRISPRi endpoint, all "
        "chromosomes: Gasperini2019 decrease distance 0.493 [0.439, 0.560], unchanged 0.425 at coverage "
        "0.55; Schraivogel2020 distance 0.390, unchanged 0.162; Morris distance 0.752; Xie distance 0.589",
        "result": "c4_holdout_scores (lane-c4)",
        "commit": "b7e4bf0 (2026-09-28)",
    },
    {
        "what": "on the 7 validation chromosomes, the prior-only labelling, the pilot, the "
        "independent-block variant, the unchanged labels and distance on the four decrease endpoints: "
        "prior-only minus distance +0.1927 [+0.1335, +0.2491] on Gasperini2019 (0.6601 against 0.4674), "
        "+0.0858 [+0.0219, +0.1814] on Xie, -0.0345 on Morris, +0.0165 [-0.4385, +0.1934] on "
        "Schraivogel2020. No labelling of the pilot was scored on the development chromosomes",
        "result": "pilot_biological_gate (lane-pilot)",
        "commit": "197c560 (2026-09-29)",
    },
    {
        "what": "this lane, before registering: the benchmark's chromatin category of every positive and "
        "pair per study (Gasperini2019: 355 of 360 positives in an H3K27ac category against 4,911 of "
        "5,299 pairs), and unit counts per partition. No score of any labelling",
        "result": "admissibility (step 1)",
        "commit": "d19ccac (2026-09-29)",
    },
    {
        "what": "the published relationship on these very screens: ABC was built on Fulco 2019's screens "
        "and ENCODE-rE2G trained on the training file",
        "result": "the literature",
        "commit": "n/a",
    },
)
S4_BESIDE = (
    "descriptive, beside and not in any rule: the prior-only labelling's (c) labels on the blocks it "
    "solved for a primary endpoint on the fresh partition, as target and context claims "
    "(pilot_bio.claims_of), and the unchanged labels' claims on the same blocks, each judged by "
    "correctness.judge with sources=[S]; target accuracy, context accuracy and coverage reported apart. "
    "(a) and (b) name no target and state no claim"
)
BUDGET_RULE = (
    "one process, cached inputs, at most 0.5 CPU-hour; no model request, and the per-element response "
    "cache is not opened"
)
MAY_BE_CALLED = "an internal development result on withheld sources"
MAY_NOT_BE_CALLED = ("a validation", "a reproduction of ABC", *ho.MAY_NOT_BE_CALLED)


def registration() -> dict[str, Any]:
    return {
        "registered": REGISTERED,
        "status": STATUS,
        "partitions": PARTITIONS,
        "primary_endpoints": [list(x) for x in PRIMARY],
        "beside_on_the_seen_partition": [list(x) for x in BESIDE_SEEN],
        "admissible_sources": list(ADMISSIBLE_SOURCES),
        "labellings": LABELLING_TEXT,
        "variants": VARIANTS,
        "abc_form": ABC_FORM,
        "comparisons": [{"a": a, "b": b, "asks": w} for a, b, w in COMPARISONS],
        "metric": METRIC,
        "useful_coverage": USEFUL_COVERAGE,
        "pass_rule": PASS_RULE,
        "decomposition_rule": DECOMPOSITION_RULE,
        "lead": LEAD,
        "falsifier": FALSIFIER,
        "readings": READINGS,
        "meaning": MEANING,
        "seen_before": list(SEEN_BEFORE),
        "s4_beside": S4_BESIDE,
        "budget_rule": BUDGET_RULE,
        "pilot_model": pl.registration(),
    }


# --- the labellings ---------------------------------------------------------------------------------
def variant_fixed(fixed: pb.Fixed, name: str) -> pb.Fixed:
    """The pilot's fixed inputs with the variant's terms removed: no compiled links, no H3K27ac peaks."""
    spec = VARIANTS[name]
    chroms = list(fixed.blocks)
    return pb.Fixed(
        fixed.blocks,
        fixed.peaks if spec["h3k27ac_peaks"] else {c: pl.Peaks({}) for c in chroms},
        fixed.tss,
        fixed.links if spec["compiled_links"] else {c: {} for c in chroms},
    )


def prior_labels(
    name: str, fixed: pb.Fixed, held: list[ho.Unit], src: str, chroms: tuple[str, ...]
) -> tuple[ho.Labels, dict[str, pb.Solved]]:
    """One prior variant as the pilot builds its prior-only labelling: `pilot_bio.solve(mode='prior')`
    over the blocks the held-out units overlap, read by `pilot_bio.labels_from`. Its own solve cache."""
    fx = variant_fixed(fixed, name)
    cache: dict[str, pl.Outcome] = {}
    solved = {}
    for chrom in chroms:
        rel = pb.relevant_blocks([u for u in held if u.chrom == chrom], fx.blocks[chrom])
        solved[chrom] = pb.solve("prior", chrom, rel, [], fx, None, cache)
    lab = pb.labels_from(name, solved, fx, VARIANT_READS[name], src, note=LABELLING_TEXT[name])
    return lab, solved


# --- the verdict ------------------------------------------------------------------------------------
def side(c: dict[str, Any] | None) -> int:
    """+1 when a paired interval lies above zero, -1 when below, 0 otherwise or unscored."""
    if not c or c.get("status") != "scored" or not c.get("interval"):
        return 0
    lo, hi = c["interval"]
    return 1 if lo > 0 else (-1 if hi < 0 else 0)


def reproduces(c: dict[str, Any] | None, lead: dict[str, Any] = LEAD) -> bool:
    if not c or c.get("status") != "scored":
        return False
    return (
        round(c["difference"], 4) == lead["difference"]
        and [round(x, 4) for x in c["interval"]] == lead["interval"]
        and round(c["a_value"], 4) == lead["a_value"]
        and round(c["b_value"], 4) == lead["b_value"]
    )


def assess(
    comparisons: list[dict[str, Any]], scores: list[dict[str, Any]], reproduced: bool
) -> dict[str, Any]:
    """PASS_RULE, DECOMPOSITION_RULE, FALSIFIER and READINGS on the primary endpoints, as registered."""
    by = {(c["source"], c["endpoint"], c["a"], c["b"]): c for c in comparisons}
    cov = {(s["source"], s["endpoint"], s["labels"]): (s.get("coverage") or {}).get("share") for s in scores}
    status = {(s["source"], s["endpoint"], s["labels"]): s.get("status") for s in scores}
    rows = []
    for src, ep in PRIMARY:

        def sd(a: str, b: str, src: str = src, ep: str = ep) -> int:
            return side(by.get((src, ep, a, b)))

        scored = status.get((src, ep, PRIOR)) == "scored"
        coverage = cov.get((src, ep, PRIOR))
        row = {
            "source": src,
            "endpoint": ep,
            "scored": scored,
            "prior_coverage": coverage,
            "useful_coverage": bool(coverage is not None and coverage >= USEFUL_COVERAGE),
            "lead": sd(PRIOR, DISTANCE),
            "activity": sd(PRIOR_ACTIVITY, PRIOR_DISTANCE),
            "compiled_target": sd(PRIOR, PRIOR_ACTIVITY),
            "machinery": sd(PRIOR_DISTANCE, DISTANCE),
            "activity_and_distance_over_distance": sd(PRIOR_ACTIVITY, DISTANCE),
            "over_unchanged": sd(PRIOR, UNCHANGED),
        }
        row["passes"] = bool(scored and row["useful_coverage"] and row["lead"] > 0)
        row["falsifier"] = {
            "i_activity_does_not_carry": row["activity"] <= 0,
            "ii_compiled_target_carries": row["compiled_target"] > 0
            and row["activity_and_distance_over_distance"] <= 0,
            "iii_machinery_carries": row["machinery"] > 0 and row["activity"] <= 0,
        }
        rows.append(row)
    scored_rows = [r for r in rows if r["scored"]]
    if not reproduced:
        reading = "void_reproduction"
    elif not scored_rows:
        reading = "no_admissible_endpoint"
    elif all(r["passes"] for r in rows):
        fired = any(any(r["falsifier"].values()) for r in rows)
        reading = "replicated_not_activity" if fired else "replicated_activity"
    elif any(r["lead"] < 0 for r in scored_rows):
        reading = "reversed"
    else:
        reading = "not_replicated"
    return {
        "endpoints": rows,
        "reproduced_the_registered_lead": reproduced,
        "reading": reading,
        "reading_text": READINGS[reading],
        "replicated": reading in ("replicated_activity", "replicated_not_activity"),
    }


# --- the run ----------------------------------------------------------------------------------------
def score_endpoint(
    src: str,
    chroms: tuple[str, ...],
    labels: dict[str, ho.Labels],
    units: dict[str, tuple[ho.Unit, ...]],
    refs: dict[str, frozenset[str]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    res = ho.evaluate(
        [labels[n] for n in LABELLINGS], src, ENDPOINT, units=units, chroms=chroms, references=refs
    )
    comps = []
    for a, b, asks in COMPARISONS:
        c = ho.compare(labels[a], labels[b], src, ENDPOINT, units=units, chroms=chroms, references=refs)
        comps.append({**c, "asks": asks})
    return res, comps


def s4_beside(
    src: str,
    solved: dict[str, dict[str, pb.Solved]],
    labels: dict[str, ho.Labels],
    unchanged_claims: list[Any],
    units: dict[str, tuple[ho.Unit, ...]],
    refs: dict[str, frozenset[str]],
) -> dict[str, Any]:
    v1 = cx.RULE_V1  # the rule this committed result was judged by, pinned (correctness.RULE_RECORD)

    def brief(rep: Any) -> dict[str, Any]:
        d = rep.to_dict()
        return {
            "claims": {a: d["axes"][a]["claims"] for a in (cx.TARGET, cx.CONTEXT)},
            "target_accuracy": d["quantities"]["target_accuracy"]["target"]
            if "target_accuracy" in d["quantities"]
            else None,
            "quantities": d["quantities"],
            "also_reported": d["also_reported"],
            "judged": d["judged"][:20],
        }

    out: dict[str, Any] = {}
    for name in (PRIOR, PRIOR_ACTIVITY, PRIOR_DISTANCE):
        claims = pb.claims_of(solved[name])
        out[name] = {"claims_stated": len(claims)}
        if claims:
            out[name].update(
                brief(cx.judge(claims, labels[name], units=units, sources=[src], references=refs, rule=v1))
            )
    same = pb.unchanged_claims_on(solved[PRIOR], unchanged_claims)
    out["unchanged_on_the_same_blocks"] = {"claims_stated": len(same)}
    if same:
        rep = cx.judge(same, cx.unchanged_labels(), units=units, sources=[src], references=refs, rule=v1)
        out["unchanged_on_the_same_blocks"].update(brief(rep))
    return out


@mf.depends_on_models("alphagenome")  # the prior's compiled target and the unchanged labels
def manifest(parameters: dict[str, Any]) -> dict[str, Any]:
    inputs = [mf.input_entry(p, partition=None) for p in ho.compiled_programs()]
    for name in ms.CRISPRI_FILES:
        p = ms.CRISPRI_KNOWLEDGE / name
        if p.exists():
            inputs.append(mf.input_entry(p, partition=ms.CRISPRI_SPLIT_OF[name]))
    inputs += [mf.input_entry(p, partition=None) for p in ho.gencode_paths()]
    for part, chroms in (("fresh", FRESH), ("seen", SEEN)):
        for c in chroms:
            inputs.append(mf.input_entry(ROOT / pb.CCRE_TEMPLATE.format(chrom=c), partition=part))
            for s in pl.BREADTH_BIOSAMPLES:
                p = ROOT / pb.PEAK_TEMPLATE.format(sample=s, chrom=c)
                if p.exists():
                    inputs.append(mf.input_entry(p, partition=part))
    return {
        "sources": [
            {
                "accession": "EngreitzLab/CRISPR_comparison EPCrisprBenchmark (Gschwind et al., Nature 2026)",
                "version": "pinned by sha256",
            },
            {"accession": "GENCODE v50 (per chromosome)", "version": "pinned by sha256"},
            {
                "accession": "ENCODE cCREs v3 (downloads.wenglab.org/V3/GRCh38-cCREs.bed), per chromosome",
                "version": "pinned by sha256",
            },
            {
                "accession": "ENCODE H3K27ac replicated peaks, 13 biosamples, per chromosome",
                "version": "pinned by sha256",
            },
            {
                "accession": "this repository: data/knowledge/compiled/noncoding_<chrom>.bio",
                "version": "pinned by sha256",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": parameters,
        "exclusions": [
            "every study of the CRISPRi held-out file, and Nasser2021 as the harness pools it: inadmissible "
            "(their positives were filtered on H3K27ac at the element); not scored",
            "CRISPRi increases: the prior states no direction; not scored",
            "Schraivogel2020 has no unit on the fresh partition; it is reported on the seen partition only",
            "CRISPRi increases, underpowered nulls and missing effects are excluded from the decrease "
            "endpoint by the harness and counted",
            "ABC's form is not computed: no experimental DNase-seq signal is held",
            "the AlphaGenome per-element response cache is not opened; the model is not called",
        ],
        "partitions": {**PARTITIONS, ms.TRAINING: "the CRISPRi benchmark's training file (K562)"},
    }


def run() -> dict[str, Any]:
    t0, c0 = time.time(), time.process_time()
    fixed = pb.Fixed.load((*FRESH, *SEEN))
    units = ho.all_units()
    refs = ho.crispri_references()
    base = {DISTANCE: ho.distance_labels(), UNCHANGED: ho.unchanged_labels()}
    load_cpu = time.process_time() - c0
    for src in {s for s, _ in (*PRIMARY, *BESIDE_SEEN)}:
        assert src in ADMISSIBLE_SOURCES, src
    counts = {
        "fresh": partition_counts(units, FRESH),
        "seen": partition_counts(units, SEEN),
    }
    out: dict[str, Any] = {"fresh": {}, "seen": {}}
    all_scores: list[dict[str, Any]] = []
    all_comps: list[dict[str, Any]] = []
    s4: dict[str, Any] = {}
    reproduced = False
    lead_seen: dict[str, Any] | None = None
    # the seen partition first: the reproduction check decides whether anything else is read
    for part, chroms, endpoints in (("seen", SEEN, BESIDE_SEEN), ("fresh", FRESH, PRIMARY)):
        for src, _ep in endpoints:
            held = [u for u in units[src] if u.chrom in chroms]
            labels = dict(base)
            solved: dict[str, dict[str, pb.Solved]] = {}
            for name in (PRIOR, PRIOR_ACTIVITY, PRIOR_DISTANCE):
                labels[name], solved[name] = prior_labels(name, fixed, held, src, chroms)
            scores, comps = score_endpoint(src, chroms, labels, units, refs)
            for r in scores:
                r["partition"] = part
            for c in comps:
                c["partition"] = part
            out[part][src] = {"scores": scores, "comparisons": comps}
            if part == "seen" and src == LEAD["source"]:
                lead_seen = next(c for c in comps if (c["a"], c["b"]) == (LEAD["a"], LEAD["b"]))
                reproduced = reproduces(lead_seen)
                if not reproduced:
                    break
            if part == "fresh":
                all_scores += scores
                all_comps += comps
                unchanged_claims = [
                    c
                    for c in cx.compiled_claims(chroms=chroms)
                    if c.axis in (cx.TARGET, cx.CONTEXT, cx.ACTIVITY)
                ]
                s4[src] = s4_beside(src, solved, labels, unchanged_claims, units, refs)
            brief = [(r["labels"], r.get("value"), r["coverage"]["share"]) for r in scores]
            print(part, src, brief, flush=True)
        if not reproduced:
            out["fresh"] = {}
            break
    verdict = assess(all_comps, all_scores, reproduced)
    peak_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)
    return {
        "question": "does the coherence pilot's prior-only labelling, exactly as the pilot defined it, "
        "predict CRISPRi decreases better than distance to TSS on admissible endpoints on chromosomes where "
        "it was never scored, and which of its terms carries the difference on identical pairs",
        "status": STATUS,
        "may_be_called": MAY_BE_CALLED,
        "may_not_be_called": list(MAY_NOT_BE_CALLED),
        "registration": registration(),
        "admissibility": admissibility(),
        "partition_counts": counts,
        "reproduction": {"lead": LEAD, "seen": lead_seen, "reproduced": reproduced},
        "verdict": verdict,
        "fresh": out["fresh"],
        "seen": out["seen"],
        "s4_beside": s4,
        "alphagenome_requests": 0,
        "per_element_response_cache_opened": False,
        "compute": {
            "loading_cpu_seconds": round(load_cpu, 1),
            "total_cpu_seconds": round(time.process_time() - c0, 1),
            "wall_seconds": round(time.time() - t0, 1),
            "peak_rss_mb": round(peak_mb, 1),
            "measured_with": "time.process_time and getrusage, this process",
        },
    }


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--admissibility", action="store_true", help="step 1's record and local counts only")
    args = ap.parse_args(argv)
    if args.admissibility:
        units = ho.crispri_units()
        rec = admissibility()
        rec["partition_counts"] = {
            "validation": partition_counts({k: tuple(v) for k, v in units.items()}, SEEN),
            "development": partition_counts({k: tuple(v) for k, v in units.items()}, FRESH),
        }
        print(
            json.dumps(
                {k: rec[k] for k in ("admissible_sources", "local_counts", "partition_counts")}, indent=1
            )
        )
        return
    payload = run()
    params = {
        "registered": REGISTERED,
        "fresh_chromosomes": list(FRESH),
        "seen_chromosomes": list(SEEN),
        "primary_endpoints": [list(x) for x in PRIMARY],
        "labellings": list(LABELLINGS),
        "variants": VARIANTS,
        "useful_coverage": USEFUL_COVERAGE,
        "resamples": ho.RESAMPLES,
        "seed": ho.SEED,
        "pilot_weights": dict(pl.W),
    }
    p = save_result(NAME, payload, manifest=manifest(params))
    v = payload["verdict"]
    print("reading", v["reading"], "reproduced", v["reproduced_the_registered_lead"], flush=True)
    print("saved", p, flush=True)


if __name__ == "__main__":
    main()
