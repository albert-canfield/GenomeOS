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
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomeos.attribution import holdout as ho  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.attribution import pilot_bio as pb  # noqa: E402

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


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--admissibility", action="store_true", help="step 1's record and local counts only")
    args = ap.parse_args(argv)
    if args.admissibility:
        units = ho.crispri_units()
        rec = admissibility()
        rec["partition_counts"] = {
            "validation": partition_counts({k: tuple(v) for k, v in units.items()}, pb.VALIDATION),
            "development": partition_counts({k: tuple(v) for k, v in units.items()}, pb.DEVELOPMENT),
        }
        print(
            json.dumps(
                {k: rec[k] for k in ("admissible_sources", "local_counts", "partition_counts")}, indent=1
            )
        )
        return
    raise SystemExit("the registered run is added with the registration (step 2)")


if __name__ == "__main__":
    main()
