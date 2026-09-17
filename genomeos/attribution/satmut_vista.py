# SPDX-License-Identifier: AGPL-3.0-or-later
"""In-silico mutagenesis over the VISTA panel: are the bases the model cares about the conserved ones?

Area J's open step asks for a different readout on the panel of §15 — 1,133 VISTA in-vivo positives
and 1,090 negatives, matched — where the finding was that motif sites held across species do not
separate an enhancer from inactive conserved sequence. Saturation mutagenesis answered the
measurement half on lentiMPRA (§16). This asks the model half, and it asks it as a comparison rather
than as a description, because "the model is sensitive at conserved bases" is not a claim that can
fail on its own: elements selected for conservation are conserved everywhere.

The pre-registration below was written and committed before the instrument existed and before any
window was scored. `PRE_REGISTRATION` is the machine-readable copy; the prose is
docs/GRAMMAR-BY-COMPARISON.md.
"""

from __future__ import annotations

from typing import Any

WINDOW = 25  # bp per perturbed window: finer than a motif, coarser than a base
WINDOWS_PER_ELEMENT = 20  # evenly spaced across the element, so length does not decide the sample
ELEMENTS_PER_ARM = 80  # the panel's own arms: 40 limb + 40 neural positives, 80 negatives
DRAWS = 200  # within-element shuffles of the constraint values

# AMENDED 2026-09-17, before any window was scored and before the instrument was written, with the
# reason: the registration below assumed the panel pairs each positive with a matched negative one to
# one. It does not. `across_panel_vista` matches as GROUPS - 40 limb, 40 neural, 80 negatives, matched
# on length, GC, constrained fraction and distance, with the medians reported per group - so a 1:1
# pairing would be this run's invention rather than the panel's. The arms are therefore the panel's
# own groups, every element of them, and the budget rises from 2,400 requests to 3,200. Nothing else
# changes: the same windows, the same outcome, the same bars, the same null. An amendment before the
# first request is a design decision; the same words after the first number would be a different thing.
AMENDMENT = (
    "arms are the panel's groups (80 positives across limb and neural, 80 negatives), not 1:1 pairs, "
    "because the panel matches at group level; 3,200 requests rather than 2,400"
)

PRE_REGISTRATION: dict[str, Any] = {
    "written": (
        "2026-09-17, before the instrument existed and before any window was scored; committed in "
        "this file before the code that runs it"
    ),
    "question": (
        "does the model's sensitivity to perturbation concentrate on conserved bases MORE in VISTA "
        "in-vivo positives than in their matched negatives? The panel's own finding (§15) is that "
        "conserved sequence and enhancer activity come apart: motif sites held across species do not "
        "separate a positive from an inactive conserved negative. If that holds, the two arms should "
        "look alike here too"
    ),
    "why_it_is_a_comparison": (
        "'the model is sensitive where the sequence is conserved' cannot fail on this panel: VISTA "
        "elements are selected for conservation, so a positive correlation is guaranteed by the "
        "sampling. The matched negatives are the only thing that makes the question answerable, and "
        "the quantity under test is the DIFFERENCE between the arms, not either arm's value"
    ),
    "sampling": (
        "every element of the panel's own arms, as amended above: 80 positives (40 limb, 40 neural) "
        f"and 80 negatives from `across_panel_vista`, each tiled with {WINDOWS_PER_ELEMENT} windows of "
        f"{WINDOW} bp, evenly spaced across it so a long element does not outvote a short one. One "
        "request per window"
    ),
    "outcome": (
        "per window: the model's predicted effect of substituting that window, and the window's mean "
        "Zoonomia phyloP from the constraint reader. Per element: the Spearman correlation between the "
        "two across its windows. Per arm: the mean of those correlations"
    ),
    "primary": (
        "the difference in mean within-element Spearman, positives minus negatives, with a one-sided p "
        "from a within-element shuffle of the constraint values (200 draws, which preserves both "
        "distributions and destroys only the pairing)"
    ),
    "outcomes": {
        "success": (
            "positives exceed negatives by at least 0.10 in mean Spearman at one-sided p 0.01: the "
            "model's sensitivity tracks conservation more tightly where the element is really an "
            "enhancer, which would be the first thing on this panel to separate the arms"
        ),
        "weak": "0.05 to 0.10, or 0.10 without the p: worth one more instrument, claims nothing",
        "failure": (
            "below 0.05: conservation and the model's sensitivity are related the same way in both "
            "arms, which is §15's finding arriving by a second route and is the outcome to expect"
        ),
    },
    "budget": f"{2 * ELEMENTS_PER_ARM * WINDOWS_PER_ELEMENT} requests, one per window (amended)",
    "stopping": (
        "one look at half the budget, futility only: if the 95% upper bound of the difference is below "
        "0.05 the run stops. No interim success look"
    ),
    "cannot_do": (
        "one model, and predicted effects rather than measured ones: this says where AlphaGenome is "
        "sensitive, not which bases matter in an embryo. Window resolution is 25 bp, so a single "
        "decisive base is averaged with 24 neighbours. And the negatives are VISTA's own negatives - "
        "sequence that failed to drive expression in the assay - which is a measured class but not a "
        "sequence-matched control, a caveat this panel has carried since §15"
    ),
}
