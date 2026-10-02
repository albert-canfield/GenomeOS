# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register N1 before any measurement: the protocol, the frozen inputs' hashes, the shared universe, the
candidate ledger and the plan `scripts/n1_run.py` applies once (data/results/n1_registration.json).

    uv run --frozen python scripts/n1_register.py
    uv run --frozen python scripts/n1_register.py --result n1_registration_v2 --supersedes

`--result` changes the name written and nothing else; `--supersedes` records which committed result the
run is written beside. The committed data/results/n1_registration.json is a byte freeze: amendment 1
declares its sha256 and tests/test_n1_perturb_response.py checks the file on disk against it, so a
re-record goes under a new name and leaves it alone.

Reads only frozen predictions (data/results/motifs_chr*.json), JASPAR 2026 and its TFClass families
(data/knowledge/jaspar/), GENCODE v50 (data/reference/gencode_v50_chr*.gff3.gz), and the identities of the
Perturb-seq pseudobulk file: its row labels and measured-gene IDs and names, fetched once by HTTP range
requests from a fixed list of four datasets and cached in data/cache/n1/. No expression value, knockdown
value, cell count or response statistic is read. docs/ATTRIBUTION.md, "N1, the Perturb-seq response test,
registered before any measurement".
"""

from __future__ import annotations

import argparse
import io
import json
import subprocess
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import n1_perturb_response as n1  # noqa: E402
from genomeos.genome import motifs  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

RESULT = "n1_registration"
#: This script, as the entry whose transitive import closure is the counting path of its result
#: (genomeos.manifest.counting_path). Named rather than taken from __file__ so the closure is the same
#: however the script is invoked.
ENTRY = "scripts/n1_register.py"
ROOT = Path(__file__).resolve().parents[1]
#: This lane's own files. Everything else uncommitted in this shared checkout belongs to another lane.
OWN_CODE = (
    "genomeos/attribution/n1_perturb_response.py",
    "scripts/n1_register.py",
    "scripts/n1_run.py",
    "tests/test_n1_perturb_response.py",
)
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
CACHE = Path("data/cache/n1")
IDENTITIES = CACHE / "K562_gwps_normalized_bulk_01.identities.json"
GENCODE = Path("data/reference")
IDENTITY_DATASETS = ("obs/gene_transcript", "var/gene_id", "var/gene_name", "var/__categories/gene_name")
BLOCK = 1 << 16

PAPER = "https://doi.org/10.1016/j.cell.2022.05.013"
PMC = "https://europepmc.org/article/PMC/PMC9380471"
FIGSHARE = "https://doi.org/10.25452/figshare.plus.20029387.v1"
PRODUCER = (
    "https://github.com/thomasmaxwellnorman/Perturbseq_GI/blob/3b25109aeb9c0c2026bd70abd50304a0ad4e5395/"
    "perturbseq/cell_population.py"
)
ANNDATA = "https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.html"


# --- identities, by range requests over four named datasets -----------------------------------------------


class RangeFile(io.RawIOBase):
    """A remote file read in 64 KiB blocks by HTTP range requests, counting every byte fetched."""

    def __init__(self, url: str, size: int):
        self.url, self.size, self.pos, self.fetched, self.cache = url, size, 0, 0, {}

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def seek(self, off: int, whence: int = 0) -> int:
        self.pos = off if whence == 0 else self.pos + off if whence == 1 else self.size + off
        return self.pos

    def tell(self) -> int:
        return self.pos

    def _block(self, i: int) -> bytes:
        if i not in self.cache:
            start, end = i * BLOCK, min(self.size, (i + 1) * BLOCK) - 1
            req = urllib.request.Request(
                self.url, headers={"Range": f"bytes={start}-{end}", "User-Agent": "GenomeOS/n1-identities"}
            )
            with urllib.request.urlopen(req, timeout=120) as r:  # noqa: S310
                if r.status != 206:
                    raise RuntimeError(f"the server did not honour the range request (HTTP {r.status})")
                data = r.read()
            self.fetched += len(data)
            self.cache[i] = data
        return self.cache[i]

    def readinto(self, b) -> int:
        n = min(len(b), self.size - self.pos)
        if n <= 0:
            return 0
        out, p = bytearray(), self.pos
        while len(out) < n:
            blk = self._block(p // BLOCK)
            off = p % BLOCK
            take = blk[off : off + n - len(out)]
            out += take
            p += len(take)
        b[:n] = out
        self.pos += n
        return n


def identities() -> dict:
    """The file's row labels and measured genes, fetched once and cached; nothing else is decoded."""
    if IDENTITIES.exists():
        return json.loads(IDENTITIES.read_text())
    import h5py

    src = n1.SOURCE["files"]["normalized"]
    rf = RangeFile(src["url"], src["bytes"])
    got = {}
    with h5py.File(io.BufferedReader(rf, buffer_size=BLOCK), "r") as h:
        for name in IDENTITY_DATASETS:
            got[name] = [x.decode() if isinstance(x, bytes) else x for x in h[name][()]]
    cats = got["var/__categories/gene_name"]
    rec = {
        "source": src["url"],
        "file_bytes": src["bytes"],
        "bytes_fetched": rf.fetched,
        "datasets_read": list(IDENTITY_DATASETS),
        "note": "identity datasets only; whole 64 KiB blocks are fetched, and nothing outside the four named "
        "datasets is decoded",
        "obs_index": [str(x) for x in got["obs/gene_transcript"]],
        "var_gene_id": [str(x) for x in got["var/gene_id"]],
        "var_gene_name": [cats[int(c)] if int(c) >= 0 else None for c in got["var/gene_name"]],
    }
    CACHE.mkdir(parents=True, exist_ok=True)
    IDENTITIES.write_text(json.dumps(rec))
    return rec


# --- GENCODE v50, the annotation the predictions used -----------------------------------------------------


def gencode_index() -> tuple[dict, dict, dict]:
    """Per chromosome, symbol -> protein-coding gene IDs (unversioned); genome-wide symbol -> IDs; and each
    gene's canonical TSS, chosen as `motifs.promoter_loci` chooses it (Ensembl_canonical, else longest)."""
    from genomeos.coords import Strand
    from genomeos.genome.annotation import Annotation

    per: dict[str, dict[str, list[str]]] = {}
    every: dict[str, set[str]] = defaultdict(set)
    tss: dict[str, tuple[str, int]] = {}
    for c in CHROMS:
        ann = Annotation.from_gff3(GENCODE / f"gencode_v50_{c}.gff3.gz", {c})
        index: dict[str, list[str]] = defaultdict(list)
        for g in ann.protein_coding():
            if g.locus.chrom != c:
                continue
            gid = g.id.split(".")[0]
            index[g.symbol].append(gid)
            every[g.symbol.upper()].add(gid)
            ts = list(g.transcripts.values())
            canon = [t for t in ts if "Ensembl_canonical" in t.tags] or sorted(
                ts, key=lambda t: -(t.locus.end - t.locus.start)
            )
            if canon and gid not in tss:
                t = canon[0]
                tss[gid] = (c, t.locus.start if t.locus.strand == Strand.PLUS else t.locus.end)
        per[c] = dict(index)
    return per, dict(every), tss


def jaspar_names() -> list[str]:
    return [
        line[1:].split()[1] for line in motifs.JASPAR_PATH.read_text().splitlines() if line.startswith(">")
    ]


def last_commit(path: Path) -> str | None:
    out = subprocess.run(["git", "log", "-1", "--format=%h", "--", str(path)], capture_output=True, text=True)
    return out.stdout.strip() or None


def unmodified(paths: list[Path]) -> bool:
    out = subprocess.run(
        ["git", "status", "--porcelain", "--", *map(str, paths)], capture_output=True, text=True
    )
    return out.returncode == 0 and not out.stdout.strip()


# --- the protocol, word for word ---------------------------------------------------------------------------


def field_semantics() -> list[dict]:
    return [
        {
            "field": "X (normalized pseudobulk)",
            "meaning": "per row, the mean over the row's cells of each cell's gemgroup z-score: z = (x - "
            "mean_control) / sd_control per gene within the cell's gemgroup, after each cell's UMI total is "
            "scaled to the median of core-control cells",
            "status": "documented: the normalization (paper, Figshare); the mean as the aggregation "
            "(producer "
            "code of 2019, and the paper's 'mean normalized expression profile'); the 2022 code is not "
            "public",
            "sources": [
                {
                    "where": f"{PMC}, STAR Methods, 'Filtering and internal normalization of gene expression "
                    "measurements'",
                    "quote": "Within each gemgroup, for each gene, we compute the mean and standard "
                    "deviation of "
                    "expression within control cells and use these to z -normalize expression.",
                },
                {
                    "where": FIGSHARE,
                    "quote": "gemgroup Z-normalized pseudo-bulk expression data for genes expressed at "
                    ">0.01 UMI "
                    "per cell (named $pop_normalized_bulk_01.h5ad)",
                },
                {
                    "where": f"{PMC}, STAR Methods, 'Global analysis and clustering of strong perturbations'",
                    "quote": "We represented perturbations by their mean normalized expression profile",
                },
                {
                    "where": f"{PRODUCER}, CellPopulation.average (lines 663-709)",
                    "quote": "the expression of each gene is equal to the average expression over all cells "
                    "in "
                    "the parent population belonging to that category",
                },
                {"where": ANNDATA, "quote": "X: A #observations x #variables data matrix."},
            ],
        },
        {
            "field": "num_cells_filtered",
            "meaning": "the number of the row's cells passing the quality filters (>= 2,000 adjusted UMIs "
            "and "
            "< 25% mitochondrial RNA in the genome-scale screen), over which X averages",
            "status": "inferred from the name and the methods; unconfirmed: the field is float64 while "
            "num_cells_unfiltered is int64",
            "sources": [
                {
                    "where": f"{PMC}, STAR Methods, 'Filtering and internal normalization'",
                    "quote": "Finally, we computed a normalized gene expression matrix for cells passing the "
                    "quality filters",
                },
                {
                    "where": f"{PRODUCER}, CellPopulation.average",
                    "quote": "'num_cells': lambda x: x.shape[0]",
                },
            ],
        },
        {
            "field": "control_expr, fold_expr, pct_expr",
            "meaning": "control_expr: the mean unnormalized expression of the perturbed gene in "
            "non-targeting "
            "control cells (the knockdown denominator); fold_expr: its mean unnormalized expression in the "
            "row's cells over control_expr; pct_expr: fold_expr - 1, the fractional change (-1 is complete "
            "loss)",
            "status": "inferred: no source defines these fields by name. The run checks pct_expr = "
            "fold_expr - 1 "
            "(or the same in percent) on every assessable candidate row and stops on any mismatch. Not "
            "known: "
            "whether the controls are all non-targeting cells or the core ones, pooled or per gemgroup, and "
            "whether 'unnormalized' means raw or depth-scaled UMIs",
            "sources": [
                {
                    "where": f"{PMC}, STAR Methods, 'Leverage scores for quantifying perturbation "
                    "penetrance'",
                    "quote": "Knockdown was computed as the ratio of mean (unnormalized) expression of the "
                    "target "
                    "gene within perturbed cells vs. that in cells with non-targeting sgRNAs.",
                },
                {
                    "where": f"{PMC}, Figure S3 legend",
                    "quote": "The fractional change in expression is defined as the expression in the "
                    "targeted "
                    "cells minus the expression in non-targeting cells, relative to the expression in the "
                    "non-targeting cell population (-1 implies 100% knockdown).",
                },
            ],
        },
        {
            "field": "energy_test_p_value",
            "meaning": "a permutation test of whether the row's cells differ from 5,000 control cells in the "
            "top 20 principal components: perturbation-level, it cannot label a gene as responding",
            "status": "documented; never read by this test",
            "sources": [
                {
                    "where": f"{PMC}, STAR Methods, 'Energy distance test'",
                    "quote": "each cell is represented by a vector composed of its top 20 principal "
                    "component "
                    "scores, and we compare whether the distribution of these 20-dimensional vectors is "
                    "equal "
                    "or not between unperturbed control cells and cells bearing each perturbation",
                }
            ],
        },
        {
            "field": "anderson_darling_counts, mann_whitney_counts",
            "meaning": "per row, counts of genes the producers' per-gene tests called (Benjamini-Hochberg); "
            "perturbation-level. The per-gene test results are not distributed in the file",
            "status": "inferred from the names (int64) and the methods; never read by this test",
            "sources": [
                {
                    "where": f"{PMC}, STAR Methods, 'Gene-level differential expression testing'",
                    "quote": "for each gene test whether the distribution of normalized expression is "
                    "identical "
                    "between control cells bearing non-targeting sgRNAs and cells bearing each perturbation",
                }
            ],
        },
        {
            "field": "mean_leverage_score, std_leverage_score",
            "meaning": "per-cell outlier scores summarised per row; perturbation-level",
            "status": "documented; never read by this test",
            "sources": [
                {"where": f"{PMC}, STAR Methods, 'Leverage scores'", "quote": "a scalar single-cell score"}
            ],
        },
        {
            "field": ".var mean, std, cv, gini; clean_mean, clean_std, clean_cv",
            "meaning": "in the 2019 producer code, statistics across the pseudobulk rows of the "
            "unnormalized means, "
            "not per-cell dispersions; the clean_* fields are undocumented",
            "status": "never read by this test",
            "sources": [
                {
                    "where": f"{PRODUCER}, MeanPopulation.__init__ (lines 994-1022)",
                    "quote": "gene_list['std'] = matrix.std()",
                }
            ],
        },
        {
            "field": "file structure",
            "meaning": "X dense float32, 11,258 rows x 8,248 genes; no layers, uns, obsm or varm: no "
            "per-gene "
            "standard error, statistic or p-value is stored anywhere in the file",
            "status": "read: names, shapes and dtypes only, by range requests (1,163,794 bytes), 2026-10-01",
            "sources": [
                {
                    "where": ANNDATA,
                    "quote": "layers: Key-indexed multi-dimensional arrays aligned to dimensions of X.",
                }
            ],
        },
        {
            "field": "X (raw pseudobulk, K562_gwps_raw_bulk_01.h5ad)",
            "meaning": "per row, the mean raw UMI count per cell of each gene; used only for the calibration "
            "gate's expression strata, at the non-targeting rows",
            "status": "documented as raw pseudobulk (Figshare); the mean as the aggregation from the 2019 "
            "code",
            "sources": [
                {
                    "where": FIGSHARE,
                    "quote": "Raw, pseudo-bulk expression data for genes expressed at >0.01 UMI per cell",
                }
            ],
        },
    ]


def protocol(n_cands: int, n_rows: int, n_universe: int) -> dict:
    return {
        "question": "Do the committed motif-to-gene predictions anticipate which genes respond to a factor's "
        "CRISPRi knockdown in K562 better than promoter proximity does?",
        "claim_scope": "The test concerns the sampled predictions: the committed scan of 300 elements per "
        "chromosome (7,063 elements), at most 66 candidate factors before quality filtering, not genome-wide "
        "attribution. A knockdown can act indirectly, so a response anticipated by a motif does not validate "
        "an enhancer-gene connection.",
        "decision_on_the_file": "The pseudobulk file alone gives a per-gene effect size (X) but no per-gene "
        "uncertainty. The endpoint adds one derivation, T = X * sqrt(num_cells_filtered), whose assumptions "
        "are stated below and checked by the calibration gate where nothing was perturbed.",
        "endpoint": {
            "statistic": "T = X * sqrt(n), n = num_cells_filtered of the row; for a factor with several "
            "eligible "
            "rows, X is their cell-weighted mean and n their summed cells (the multi-row rule)",
            "assumptions": [
                "(a) X averages exactly the num_cells_filtered cells of the row (inferred, unconfirmed)",
                "(b) where a gene does not respond, its per-cell z-scores in the row follow the control "
                "distribution, so their variance is 1 (true by construction within each gemgroup)",
                "(c) cells are independent",
                "(d) the mean of n per-cell z-scores is close enough to normal at |T| = 3; weakest for "
                "sparse "
                "genes, where calls would be too liberal and fall on low-expression genes, which may "
                "correlate "
                "with the motif scores",
            ],
            "responds": "|T| >= 3, two-sided: a gene that rises or falls responds",
            "threshold_rationale": "the two-sided normal tail at 3 is 0.27%, about two chance responders per "
            f"factor over ~{n_universe} genes; 2.58 (1%) would add about eight more, diluting the labels, "
            "and 4 "
            "(0.006%) would leave many factors without responders and so with undefined AUROCs. No "
            "multiplicity "
            "correction: the labels are the outcome of a ranking metric, not discoveries, and chance "
            "responders "
            "that are independent of the scores pull both arms toward 0.5 alike. Chosen before any value was "
            "seen; the one primary rule, with no alternative tried later.",
        },
        "calibration_gate": {
            "when": "first, once measurement access is authorised, before any factor row is read",
            "rows": "all 585 non-targeting rows, selected from the identity index alone",
            "why_585_not_514": "selecting the 514 core controls would need the core_control value column; "
            "the "
            "producers chose the core set for showing few differentially expressed genes, and it defines the "
            "normalization's control mean and sd, so it would flatter the null. The other 71 carry the "
            "guide-level effects any factor guide can also carry, so all 585 err toward stopping.",
            "values_read": [
                "normalized X at the 585 non-targeting rows, universe genes only",
                "normalized obs/num_cells_filtered at the 585 non-targeting rows",
                "raw X at the 585 non-targeting rows, universe genes only (expression strata)",
            ],
            "strata": "deciles of control-row expression over the universe genes: the cell-weighted mean of "
            "raw X "
            "over the usable non-targeting rows",
            "pass_rule": "the share of |T| >= 3 over all non-targeting rows x universe genes is at most 1.5 "
            "times "
            "the nominal 0.0026998, overall and in each of the 10 deciles, and each decile expects at least "
            "20 "
            "exceedances. Anything else fails.",
            "housekeeping": "a row with a non-finite or sub-1 cell count is left out and counted; a gene "
            "with any "
            "non-finite T or expression leaves the universe before any factor row is read, and is listed",
            "on_failure": "the experiment stops and the failure is recorded. No fallback endpoint inside "
            "this "
            "run: the authors' per-gene results or the single-cell file would be a new proposal.",
            "cannot_test": "whether a perturbed row's cells have a different variance from controls",
        },
        "eligibility": {
            "when": "after the gate passes; reads obs/num_cells_filtered, control_expr, fold_expr and "
            "pct_expr "
            f"at the {n_rows} candidate rows only",
            "rules_in_order": [
                "num_cells_filtered finite and >= 25, else too_few_cells (the producers' minimum: 'at least "
                "25 "
                "cells that passed our quality filters')",
                "control_expr finite and > 0 with num_cells_filtered * control_expr >= 10 expected target "
                "UMIs, "
                "else knockdown_unassessable: a target too weakly expressed to measure its knockdown is "
                "ineligible, never assumed knocked down (unlike the producers' rule, which admits undetected "
                "targets). At 10 expected UMIs a non-functional guide shows <= 4 with Poisson probability "
                "0.029",
                "fold_expr finite and >= 0, else knockdown_unassessable",
                "pct_expr = fold_expr - 1 (or the same in percent), else semantics_mismatch, which stops "
                "the run",
                "fold_expr <= 0.40 (>= 60% knockdown, the producers' threshold for analyses needing a "
                "functional "
                "perturbation: 'an on-target knockdown of at least 60%'; median knockdown 85.5% in K562), "
                "else "
                "insufficient_knockdown",
            ],
            "denominator": "control_expr, read as the mean unnormalized expression of the target in "
            "non-targeting cells (an inference; see field_semantics)",
            "floor_before_responses": "fewer than 30 knockdown-eligible factors stop the run before any "
            "response row is read",
        },
        "multi_row_rule": "every knockdown-eligible row of a factor is pooled, weighted by its cells; an "
        "ineligible row is never used; no row is chosen by its response",
        "aliases": "a JASPAR monomer name that is a GENCODE v50 protein-coding gene name with one gene ID "
        "maps to "
        "that ID, matched to the perturbation index by Ensembl ID; a name GENCODE does not know falls back "
        "to "
        "the perturbation index's own symbol when it names one gene ID; a name GENCODE resolves is never "
        "re-matched by symbol",
        "heterodimers": "a heterodimer profile (A::B) is attributed to neither partner: its hits score for "
        "no "
        "factor in either arm, and it is no candidate",
        "families": "factors of one TFClass family stay separate factors with equal weight; the uncertainty "
        "is "
        "clustered by TFClass unit (motifs.family_unit: the family from JASPAR's TRANSFAC file; a C2H2 "
        "zinc-finger factor or an unannotated one is its own unit)",
        "scores": {
            "proximity": "the factor's hit score in the gene's promoter requires list (TSS +- 1,000 bp of "
            "the "
            "canonical transcript, GENCODE v50); 0 when the promoter was scanned and the factor is not "
            "listed",
            "attribution": "the largest of the factor's hit scores over the sampled elements whose target "
            "is the "
            "gene; 0 when such elements were scanned and none lists the factor",
            "zero_means": "the factor is not in the frozen top-8 requires list of any declared scanned "
            "region of "
            "the gene: no hit at 85% of the matrix range, or a hit ranked below the top 8 enriched factors, "
            "or a "
            "factor not enriched on that chromosome. Never 'not scanned': an unscanned gene is outside the "
            "universe",
        },
        "universe": "measured genes (by Ensembl ID) with a scanned promoter and at least one sampled element "
        "attributed to them, the same for both arms; a gene with a non-finite control T leaves it at the "
        "gate; "
        "per factor, minus the perturbed gene and every gene whose canonical TSS lies within 10 kb of the "
        "factor's (CRISPRi neighbour silencing, Replogle 2022 Figure S3)",
        "candidates": "a perturbed monomer factor with >= 10 nonzero genes in each arm of its universe: "
        f"{n_cands} "
        "factors, at most 66 before quality filtering",
        "auroc": "per factor and arm, the Mann-Whitney AUROC of the score for responders over the factor's "
        "universe, midranks for ties. Undefined when no gene responds, every gene responds, or an arm's "
        "scores "
        "do not vary; an undefined factor is counted in coverage with its reason and never dropped silently",
        "comparison": "per factor, d = AUROC(attribution) - AUROC(proximity); the estimate is the mean of d "
        "over "
        "factors with both AUROCs defined, each factor weighted equally",
        "uncertainty": {
            "method": "cluster bootstrap over TFClass units: 10,000 replicates (seed 20261001), each "
            "drawing as "
            "many units as there are, with replacement, and averaging every factor in the drawn units with "
            "equal weight; the 95% interval is the 251st and 9,750th of the sorted means",
            "assessable_only_with": ">= 10 units",
            "limits": [
                "with few units a percentile interval undercovers",
                "factors in different families share downstream programmes and the same universe genes, "
                "which "
                "the bootstrap does not resample",
                "the noise of the responder labels is not propagated",
                "the element sample is fixed, not resampled: the result is conditional on it",
            ],
        },
        "criteria": {
            "point_gain": "the estimate is at least +0.02 AUROC (the proposed useful gain, with no claim the "
            "study can detect it)",
            "lower_bound": "the 95% cluster-bootstrap lower bound is above 0",
            "kept_separate": "each is reported on its own and never merged into one verdict",
        },
        "floor": "30 defined paired differences: a feasibility floor, not shown power; fewer stop the "
        "analysis "
        "with no estimate, and no rule is relaxed",
        "limits": [
            "the element arm is a sample of 300 elements per chromosome, not the 440,377 compiled elements; "
            "a "
            "full scan would be a new frozen prediction set",
            "requires lists are truncated at 8 factors per region, so a zero includes hits ranked below the "
            "top 8",
            "one cell line (K562), one time point (day 8), CRISPRi knockdown",
            "the field semantics of fold_expr, control_expr, pct_expr and num_cells_filtered are inferred",
            "the gate cannot test whether perturbed cells have a different variance from controls",
            "no prior project exposure to this dataset was found; whether upstream data or models saw "
            "related "
            "data cannot be established, so independence is not complete",
        ],
    }


def access_required(n_rows: int, n_universe: int) -> dict:
    files = n1.SOURCE["files"]
    return {
        "authorisation": "the owner's separate authorisation, recorded by scripts/n1_run.py --authorisation",
        "in_order": [
            f"download {files['normalized']['name']} and {files['raw']['name']} "
            f"({files['normalized']['bytes']:,} "
            f"bytes each), check Figshare's md5 ({files['normalized']['md5']}, {files['raw']['md5']}) and "
            "record "
            "each file's sha256",
            "read both files' identities (obs/gene_transcript, var/gene_id); they must match the frozen "
            "digest",
            f"calibration gate: normalized X and raw X at the 585 non-targeting rows over the {n_universe} "
            "universe genes, and obs/num_cells_filtered at those rows",
            "only if the gate passes: obs/num_cells_filtered, control_expr, fold_expr, pct_expr at the "
            f"{n_rows} "
            "candidate rows",
            "only if every assessable candidate row agrees with the field semantics and >= 30 factors are "
            "eligible: normalized X at the eligible candidate rows over the calibrated universe genes",
        ],
        "never_read": [
            "energy_test_p_value",
            "anderson_darling_counts",
            "mann_whitney_counts",
            "mean_leverage_score",
            "std_leverage_score",
            "core_control",
            "UMI_count_unfiltered",
            "num_cells_unfiltered",
            "z_gemgroup_UMI",
            "mitopercent",
            "TE_ratio",
            "cnv_score_z",
            "any .var statistic",
            "any row that is neither non-targeting nor a candidate's",
        ],
        "once": "scripts/n1_run.py refuses to run when data/results/n1_result.json exists",
    }


# --- assembling the registration ------------------------------------------------------------------------


#: Why this registration is re-recorded under a name of its own rather than over the committed file.
#: The committed manifest records its 24 GENCODE files as one group under a label with no member list,
#: so a rebuild in a second environment has no path to open or hash: `files_entry` did not name the
#: members of a group until 2026-10-02. The committed bytes cannot be replaced to fix it: the sha256 of
#: data/results/n1_registration.json is a declared input of three other committed results
#: (n1_registration_amendment_1, n1_registration_amendment_2, n1_result_amendment_2), and
#: tests/test_n1_perturb_response.py checks the file on disk against the sha256 amendment 1 froze, so
#: new bytes under the same name would make all three unrebuildable and fail that test, correctly.
WHY_A_NEW_NAME = (
    "re-recorded under a new name on 2026-10-02 so that every group of input files names its members, "
    "each with its own sha256 and byte count, which is what a rebuild in a second environment needs to "
    "open and hash them; no figure of the run differs. The committed result is kept unchanged because "
    "its sha256 is a declared input of three other committed results and is the freeze amendment 1 "
    "checks"
)


def supersedes(root: Path = ROOT) -> dict[str, Any]:
    """The committed result this run is written beside, never over, with the bytes it is kept at."""
    p = root / RESULTS_DIR / f"{RESULT}.json"
    entry = mf.input_entry(p, partition=None)
    old = json.loads(p.read_text())
    return {
        "file": (RESULTS_DIR / f"{RESULT}.json").as_posix(),
        "date": old.get("date"),
        "sha256": entry["sha256"],
        "bytes": entry["bytes"],
        "kept": "unchanged; this run is written beside it under a new name, not over it",
        "why": WHY_A_NEW_NAME,
        "declared_as_an_input_by": [
            "data/results/n1_registration_amendment_1.json",
            "data/results/n1_registration_amendment_2.json",
            "data/results/n1_result_amendment_2.json",
        ],
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--result",
        default=RESULT,
        help=(
            "the result name to write, and nothing else: no figure, no input and no parameter of the "
            "run is read from it. A re-record under a new name leaves the committed file alone"
        ),
    )
    ap.add_argument(
        "--supersedes",
        action="store_true",
        help=(
            f"record, beside the manifest, that this run is written beside the committed "
            f"{RESULT} rather than over it, with that file's date and sha256"
        ),
    )
    args = ap.parse_args()
    result_name = args.result

    motif_paths = [RESULTS_DIR / f"motifs_{c}.json" for c in CHROMS]
    if not unmodified(motif_paths):
        print(
            "refused: a motifs_chr*.json file differs from the commit; the predictions must be the "
            "committed ones"
        )
        return 2
    results = {c: json.loads(p.read_text()) for c, p in zip(CHROMS, motif_paths, strict=True)}
    ident = identities()
    per_chrom, every, tss = gencode_index()
    pred = n1.predictions(results, per_chrom)
    uni = n1.shared_universe(ident["var_gene_id"], pred)
    universe = uni["genes"]
    rows = [{**n1.parse_row(label), "index": i} for i, label in enumerate(ident["obs_index"])]
    control_rows = [r["index"] for r in rows if r["non_targeting"]]
    fams = motifs.load_families()
    names = jaspar_names()
    cluster_of = {
        nm: motifs.family_unit(nm, fams)
        for nm in {x for x in (n1.monomer_name(j) for j in names) if x is not None}
    }
    ledger = n1.candidate_ledger(names, every, rows, pred, universe, tss, cluster_of)
    cands = [line for line in ledger if line["status"] == "candidate"]
    ids = [c["gene_id"] for c in cands]
    if len(ids) != len(set(ids)):
        print("refused: two candidate names share one perturbed gene")
        return 2
    n_rows = sum(len(c["rows"]) for c in cands)
    step_counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for line in ledger:
        for step, outcome in line["steps"].items():
            step_counts[step][outcome.split(":")[0].split(" (")[0]] += 1
    element_ge10 = sum(
        1 for line in ledger if line.get("positives", {}).get("attribution", 0) >= n1.MIN_POSITIVES
    )
    plan = {
        "identity_digest": n1.identity_digest(ident["obs_index"], ident["var_gene_id"]),
        "universe": universe,
        "control_rows": control_rows,
        "candidates": [
            {k: c[k] for k in ("factor", "cluster", "gene_id", "rows", "excluded_genes", "scores")}
            for c in cands
        ],
    }
    gencode_files = [GENCODE / f"gencode_v50_{c}.gff3.gz" for c in CHROMS]
    inputs = [mf.input_entry(p, partition=None, last_commit=last_commit(p)) for p in motif_paths]
    inputs += [
        mf.input_entry(motifs.JASPAR_PATH, partition=None, url=motifs.JASPAR_URL),
        mf.input_entry(motifs.JASPAR_TRANSFAC_PATH, partition=None, url=motifs.JASPAR_TRANSFAC_URL),
        mf.files_entry("data/reference/gencode_v50_chr*.gff3.gz (24 files)", gencode_files, partition=None),
        mf.input_entry(IDENTITIES, partition=None, url=n1.SOURCE["files"]["normalized"]["url"]),
    ]
    frozen = [
        {k: e[k] for k in ("path", "sha256", "bytes") if k in e}
        | ({"files": e["files"]} if "files" in e else {})
        for e in inputs
    ]
    payload = {
        "status": "registered before any measurement; not run",
        "lane": "lane-n1",
        "protocol": protocol(len(cands), n_rows, len(universe)),
        "field_semantics": field_semantics(),
        "source": {
            **{k: v for k, v in n1.SOURCE.items() if k != "files"},
            "files": {
                k: {
                    **v,
                    "sha256": None,
                    "sha256_note": "computed at the authorised download, then recorded by the run",
                }
                for k, v in n1.SOURCE["files"].items()
            },
        },
        "identities_read": {k: ident[k] for k in ("source", "file_bytes", "bytes_fetched", "datasets_read")}
        | {
            "rows": len(ident["obs_index"]),
            "non_targeting_rows": len(control_rows),
            "measured_genes": len(ident["var_gene_id"]),
            "identity_digest": plan["identity_digest"],
            "structure_read_bytes": 1_163_794,
        },
        "coverage": {
            **uni["coverage"],
            "exclusions": uni["exclusions"],
            "prediction_mapping_drops": dict(sorted(pred.dropped.items())),
            "promoters_scanned_genes": len(pred.promoter_scanned),
            "element_target_genes": len(pred.element_attributed),
        },
        "ledger_summary": {
            "jaspar_names": len(ledger),
            "by_step": {s: dict(v) for s, v in step_counts.items()},
            "perturbed_with_element_ge_10": element_ge10,
            "candidates": len(cands),
            "candidate_rows": n_rows,
            "candidates_with_several_rows": sorted(c["factor"] for c in cands if len(c["rows"]) > 1),
            "clusters": len({c["cluster"] for c in cands}),
        },
        "candidates": [
            {
                "factor": c["factor"],
                "gene_id": c["gene_id"],
                "via": c["via"],
                "cluster": c["cluster"],
                "rows": [r["label"] for r in c["rows"]],
                "positives": c["positives"],
                "universe_genes": c["universe_genes"],
                "excluded_genes": len(c["excluded_genes"]),
            }
            for c in cands
        ],
        "measurement_access_required": access_required(n_rows, len(universe)),
        "frozen_inputs": frozen,
        "constants": n1.CONSTANTS,
        "code": {
            "module": "genomeos/attribution/n1_perturb_response.py",
            "tests": "tests/test_n1_perturb_response.py (synthetic inputs only)",
            "register": "scripts/n1_register.py",
            "run": "scripts/n1_run.py",
        },
        "plan": plan,
        "ledger": [{k: v for k, v in line.items() if k != "scores"} for line in ledger],
    }
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "Replogle et al. 2022 processed Perturb-seq, Figshare+ 20029387 (identities "
                "only)",
                "version": "v1",
                "doi": n1.SOURCE["doi"],
                "license": "CC BY 4.0",
            },
            {
                "accession": "JASPAR 2026 CORE vertebrates non-redundant",
                "version": "2026",
                "url": motifs.JASPAR_URL,
            },
            {"accession": "GENCODE human gene annotation", "version": "release 50 (GRCh38.p14)"},
            {
                "accession": "GenomeOS motifs_chr*.json (committed predictions)",
                "version": ", ".join(sorted({e["last_commit"] for e in inputs[:24] if e.get("last_commit")})),
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": dict(n1.CONSTANTS),
        "exclusions": [
            {"what": k, "count": v} for k, v in {**uni["exclusions"], **dict(pred.dropped)}.items()
        ]
        + [{"what": "JASPAR names not candidates, by the step they left", "count": len(ledger) - len(cands)}],
        "partitions": {
            "calibration": "the 585 non-targeting rows, read first",
            "evaluation": "the knockdown-eligible candidate rows, read only after the gate passes",
        },
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }
    if args.supersedes:
        payload["result_manifest"]["supersedes"] = supersedes()
    # Both branches are live: the name --result gave, or this writer's own default. A ternary, which
    # ruff would prefer as SIM108, re-indents the second call, and scripts/check_staged.py holds that
    # exact line because b359867 added it on 2026-10-01 13:39. Keeping it at its own indentation is
    # what lets this commit go in with no --force. After 2026-10-03 13:39 the two collapse into one
    # call on result_name and this annotation goes with them.
    if result_name != RESULT:  # noqa: SIM108
        p = save_result(result_name, payload)
    else:
        p = save_result(RESULT, payload)
    print(
        f"wrote {p}: universe {len(universe)} genes, {len(cands)} candidates ({n_rows} rows, "
        f"{payload['ledger_summary']['clusters']} clusters), element >= 10: {element_ge10}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
