# A predicted-deletion feature on the ENCODE CRISPRi enhancer–gene benchmark

*GenomeOS, 2026-09-27. Status: internal. Nothing here has been published or submitted.*

## The claim

On the ENCODE CRISPRi enhancer–gene benchmark, adding one feature to a simple activity-and-distance
model raises held-out K562 AUPRC from 0.55 to 0.69 (+0.14, 95% +0.08 to +0.23). The feature is the
expression drop AlphaGenome predicts when the enhancer is deleted, and the model was frozen before
the held-out pairs were scored. The result has been tested in **one cell type**.

## The setup

- **Data.** The ENCODE enhancer–gene benchmark (EngreitzLab/CRISPR_comparison; Gschwind et al.,
  *Nature* 2026). The training set is 10,356 K562 element–gene pairs (471 regulated). The held-out
  set is 4,378 pairs in five cell types (190 regulated). A pair is regulated when silencing the
  element significantly lowers the gene.
- **Baseline.** A logistic model on log distance, log activity and activity over distance, where
  activity is sqrt(DNase × H3K27ac) at the element in the screen's own cell.
- **Added feature.** For every registry enhancer overlapping the tested element, AlphaGenome's
  predicted log2 fold change of the measured gene when that enhancer is deleted, on the cell's own
  RNA-seq track. It enters as the largest predicted drop, plus a flag for "this gene is the
  element's top predicted target".
- **Protocol.** Fitted on K562 training pairs only. The feature set and the pass condition were
  committed before any held-out pair was scored, and nothing was refitted afterwards.

## The numbers

**Headline, held-out K562:** pairs on a scored element, average precision, 95% intervals from
resampling whole chromosomes.

| model | AUPRC | n |
|---|---|---|
| activity + distance | 0.550 | 1,744 pairs, 114 regulated |
| + predicted deletion | **0.691** | same |
| gain | **+0.141 [+0.082, +0.231]** | |

**Beside the published models, on the same pairs and with the benchmark's own estimator.** The
estimator reproduces the paper's distance-to-TSS figure exactly on both sets: 0.4359 and 0.3631.

| set | this project | published (Gschwind et al. 2026, Suppl. Table 3) |
|---|---|---|
| training, 10,356 pairs, leave-one-chromosome-out | baseline 0.507; **+ deletion 0.724** | ABC 0.565 [0.511, 0.610]; ENCODE-rE2G 0.662 [0.616, 0.706]; ENCODE-rE2G Extended 0.737 [0.693, 0.775] |
| held-out, 4,378 pairs, five cells pooled, weighted | DNase-only baseline 0.476; **+ deletion 0.639**\* | ABC 0.465 [0.378, 0.541]; ENCODE-rE2G 0.556 [0.468, 0.631] |

\* The DNase-only row was a post hoc choice; the reason is under *Limits*. The pre-registered model,
which also reads H3K27ac, scores 0.567 and 0.677 on this set.

**Read plainly.** This project's own baseline is weaker than published ABC. With the deletion
feature it reaches the range of the published ENCODE-rE2G models on the same pairs. The comparison
is against published intervals, not a paired test, so it supports "in the range of", not
"better than".

**The gain does not come from which pairs were scored.** Regulated pairs were more often on a
scored element (96.6% against 90.6% in held-out K562), so three checks were registered:

| check | gain |
|---|---|
| every held-out K562 pair, unscored ones at zero | +0.136 [+0.081, +0.225] |
| coverage matched, 1,000 random draws | median +0.145, all 1,000 above zero |
| a "was it scored" flag in place of the deletion | +0.001 [−0.001, +0.003] |

## Limits

1. **One cell type.** The second held-out cell type with usable n is HCT116 (396 pairs, 34
   regulated). AlphaGenome has HCT116 tracks, but scoring it costs 705 model requests, and that has
   not been done. GM12878 (68 pairs, 16 regulated) gives +0.015 [−0.094, +0.205]: uninformative.
   WTC11 (15 regulated) and Jurkat (7) are too small to carry a result.
2. **The held-out positives were selected on H3K27ac.** All 190 sit in H3K27ac elements, and 1,438
   of the 4,188 negatives do not. A model that reads H3K27ac, as the pre-registered baseline does,
   is flattered on that set. ENCODE-rE2G's published held-out model reads DNase only. That is why the
   held-out comparison above is quoted with DNase in place of sqrt(DNase × H3K27ac). With that
   change the gain grows rather than shrinks: +0.164 pooled, +0.179 on K562.
3. **Feature selection.** The deletion features were chosen after looking at single predictors on
   the training pairs. The held-out figures are free of this; the training comparison is not.
   ENCODE-rE2G's features were also selected on the training set.
4. **Not a new idea.** The AlphaGenome preprint (Avsec et al. 2025, bioRxiv 10.1101/2025.06.25.661532,
   Fig. 4j) added an AlphaGenome input-gradient score to ENCODE-rE2G on this dataset. What this
   project adds is the deletion form of the feature and a frozen held-out test.
5. **The model is not sequence-naive.** AlphaGenome was trained on ENCODE tracks genome-wide,
   including K562 RNA-seq. It was not trained on CRISPR outcomes, so the labels are held out but
   the genomic regions are not.

## Where everything is

- Code: `genomeos/attribution/crispri.py` (`score`, `score_published`, `PREREGISTERED`,
  `PREREGISTERED_PUBLISHED`, `PUBLISHED`).
- Scripts: `scripts/crispri_score.py`, `scripts/crispri_published.py`.
- Results: `data/results/crispri_benchmark.json`, `data/results/crispri_published.json`.
- Full account, with the registrations and their dates: `docs/ATTRIBUTION.md`, the CRISPRi
  sections of 2026-09-16, 2026-09-22 and 2026-09-27.
- Published figures: Gschwind et al. 2026, doi:10.1038/s41586-026-10781-4,
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC13471189/>, Supplementary Table 3.
