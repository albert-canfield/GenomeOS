# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-recovered (ROADMAP §5, area I, 2026-10-12): do the genes scaling makes fittable read better?

`lane-scalefit` (`8469df9` rule and floor committed before counting, `7f65544` census) established, with
NO RNA OPENED, that 7d excluded 12,500 of the 17,983 multi-element genes as `basal_not_positive` -- their
summed single-removal drops reach or exceed the gene's whole intact level -- and that 8,792 of those 12,500
(70.34%) get a positive basal rate back when every drop is multiplied by M2h's fixed constant 0.2146
(`ac5580b`), of which 8,535 pass 7d's full runnable test. Its registered reading claimed NOTHING about
expression and said the RNA test is a separate registration.

THIS MODULE IS THAT REGISTRATION. It asks one question about measured RNA on that recovered population,
and it is built in two committed steps: the population first (POPULATION, `scripts/recovered_population.py`,
`data/results/recovered_population.json`), then the registration (`scripts/recovered_register.py`,
`data/results/recovered_registration.json`), then one run.

The stop list of 2026-09-27 (docs/ROADMAP.md §5) is carried forward: no AlphaGenome request and no
deletion request, no new reader biosample (the 13 covered cells only), nothing in area J, no node-model
work, no BioForge calibration.
"""

from __future__ import annotations

import hashlib
from typing import Any

from genomeos.attribution import combine2 as c2
from genomeos.attribution import magcal as mc
from genomeos.attribution import multrna as mr
from genomeos.attribution import scalefit as sf

LANE = "lane-recovered"
REGISTERED = "2026-10-12"

# ==== step 1: the population, established and committed before the registration ==========================

#: the form whose exclusion defines the population: the MEAN, 7d's runtime default and scalefit's primary
POPULATION_FORM = "mean"

#: the form whose prediction is read against RNA: the MULTIPLICATIVE composition, M2h's arm that read
#: "better". The population form and the prediction form are DIFFERENT and that is deliberate: the
#: scaling recovers genes under the mean's fit, and the question is whether M2h's calibrated
#: multiplicative run -- the one arm of M2h that read better -- is also closer to measured RNA than the
#: no-rule baseline on exactly those genes.
FORM = mr.FORM

POPULATION = (
    "THE GENES SCALING MAKES FITTABLE, re-derived here and not copied from any count. Starting from the"
    " 17,983 multi-element genes 7d/7e found eligible for a fit (`combine2`'s own eligibility, recomputed"
    " from the same 24 compiled non-coding programs), a gene is in this population when all four hold:"
    " (a) under the MEAN form, 7d's fit on the unscaled single-removal observations fails with the reason"
    " `basal_not_positive` (combine2.fit's `if not q < 1.0`), which is the 12,500; (b) under the MEAN"
    " form, 7d's SAME code on observations whose every single-removal effect is multiplied by"
    " scalefit.SCALE (`RHO_i -> 1 - SCALE * (1 - RHO_i)`) returns a positive basal rate, which is"
    " scalefit's PRIMARY `recovered`; (c) that scaled fit also passes 7d's FULL runnable test, i.e. every"
    " fitted strength lies in (0, 1], which is scalefit's SECONDARY `runnable_after`; (d) the gene has a"
    " row in 7h's committed multiplicative genes table (data/results/multrna_result.json, genes_table),"
    " which is what supplies this lane's P and D. (a) to (c) are predictor-side and are recomputed by"
    " scalefit's own functions; (d) is a lookup of 7h's population by SYMBOL ONLY and reads no value."
    " The ENDPOINT is read on the HELD-OUT half of this population under magcal.is_held_out, so that no"
    " gene whose RNA entered M2h's fit of the scale is in the number that decides"
)

#: the genes of (a)+(b)+(c), 7h's table not yet intersected: scalefit's committed `runnable_after`
POPULATION_SOURCE_COUNTS = {
    "source": "data/results/scalefit_census.json, per_form.mean, commit 7f65544 (lane-scalefit)",
    "eligible_for_a_fit": 17_983,
    "excluded_basal_not_positive": 12_500,
    "recovered": 8_792,
    "runnable_after": 8_535,
    "note": (
        "these are scalefit's committed counts, quoted so that the re-derivation has something to"
        " disagree with. THE POPULATION THIS LANE USES IS THE ONE IT RE-DERIVES, and a disagreement with"
        " any of these four numbers stops the lane (POPULATION_HARNESS)"
    ),
}

POPULATION_HARNESS = (
    "before any registration is written and before any RNA value is read: the re-derivation must"
    " reproduce all four of POPULATION_SOURCE_COUNTS' counts exactly, from the same compiled programs and"
    " the same code; the recovered genes must be a subset of the excluded genes and the runnable a subset"
    " of the recovered; every symbol must be non-empty and occur once. A miss in any of them is"
    " POPULATION_DEFECT and nothing else is computed"
)

#: what the overlap with M2h's own population is measured against, and why it must be stated
OVERLAP_QUESTION = (
    "M2h read the multiplicative composition on 16,509 genes, 8,181 training and 8,328 held out, and its"
    " held-out reading of +0.04339 is over those 8,328. The brief that opened this lane stated that the"
    " recovered genes 'were in no previous reading'. THAT IS A CLAIM ABOUT POPULATIONS AND IT IS MEASURED"
    " HERE, NOT ASSUMED: the overlap of the re-derived recovered population with M2h's 16,509, with its"
    " 8,181 training genes and with its 8,328 held-out genes is counted and reported, and whichever way"
    " it comes out is what the registration's exposure sentence says. 7d's mean population (4,725 genes)"
    " is disjoint from the recovered population BY CONSTRUCTION -- a gene 7d excluded is not a gene 7d"
    " fitted -- but 7h's multiplicative population is NOT selected by the mean's basal rate and so may"
    " contain the recovered genes in full"
)


# ==== the population, re-derived ========================================================================


def is_recovered(rho: dict[str, float]) -> tuple[bool, bool, bool]:
    """(excluded, recovered, runnable) for one gene under POPULATION_FORM, by scalefit's own code.

    `excluded` is POPULATION (a), `recovered` is (b) and `runnable` is (c). Nothing here is fitted:
    `scalefit.SCALE` is M2h's constant and `scalefit.recovery` calls `combine2.fit` twice.
    """
    rec = sf.recovery(POPULATION_FORM, rho)
    excluded = rec.before_reason == "basal_not_positive"
    return excluded, rec.recovered, rec.runnable


def held_out_split(symbols: list[str]) -> tuple[list[str], list[str]]:
    """(training, held out) under magcal's SPLIT, carried unchanged: sha256 parity of the symbol."""
    held = sorted(s for s in symbols if mc.is_held_out(s))
    train = sorted(s for s in symbols if not mc.is_held_out(s))
    return train, held


def sha256_of_symbols(symbols: list[str]) -> str:
    """A pin on the population itself, so the registration and the run cannot drift from step 1."""
    h = hashlib.sha256()
    for s in sorted(symbols):
        h.update(s.encode("utf-8"))
        h.update(b"\n")
    return h.hexdigest()


def overlap(mine: set[str], theirs: set[str]) -> dict[str, Any]:
    """The measured relation between two symbol sets, with the word for it."""
    both = mine & theirs
    return {
        "mine": len(mine),
        "theirs": len(theirs),
        "in_both": len(both),
        "mine_only": len(mine - theirs),
        "theirs_only": len(theirs - mine),
        "share_of_mine_in_theirs": (len(both) / len(mine)) if mine else None,
        "disjoint": not both,
        "mine_is_a_subset_of_theirs": mine <= theirs,
    }


def overlap_sentence(ov: dict[str, Any], theirs_name: str) -> str:
    """The plain words for an overlap, used in the registration's exposure sentence."""
    if ov["disjoint"]:
        return (
            f"the recovered population is DISJOINT from {theirs_name} ({ov['mine']} against {ov['theirs']})"
        )
    if ov["mine_is_a_subset_of_theirs"]:
        return (
            f"the recovered population is WHOLLY CONTAINED in {theirs_name}: all {ov['in_both']} of its"
            f" genes are among those {ov['theirs']}"
        )
    return (
        f"the recovered population OVERLAPS {theirs_name}: {ov['in_both']} of its {ov['mine']} genes are"
        f" among those {ov['theirs']} ({100 * ov['share_of_mine_in_theirs']:.2f}%), and"
        f" {ov['mine_only']} are not"
    )


#: the engine constants the population depends on, recorded so a change is visible
CONVENTIONS = {
    "scale": sf.SCALE,
    "scale_rounded": sf.SCALE_ROUNDED,
    "scale_is_refitted": False,
    "population_form": POPULATION_FORM,
    "prediction_form": FORM,
    "theta": c2.THETA,
    "selection": c2.SELECTION,
    "cells": list(c2.CELLS),
    "split_rule": mc.SPLIT,
}
