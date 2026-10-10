# SPDX-License-Identifier: AGPL-3.0-or-later
"""The small-intestinal crypt as a mechanism module: one stem compartment to four mature lineages.

Area E's second human mechanism module after haematopoiesis, and the second place in the human body
program where a lineage choice is decided by a named factor instead of a share. The gates come from
the closed, sourced list registered in scripts/intestine_registration.py BEFORE any number here was
fitted (docs/AREA-E-NEXT.md holds the written choice): eight factors, each with a published
loss-of-function phenotype, and five more excluded with the reason.

What is measured, what is inferred and what is fitted is kept apart on purpose:

- measured and used as given: the epithelial count and the daily output (Sender & Milo 2021), the
  goblet and Paneth turnover times (Cheng 1974, mouse), the stem-cell division interval (Snippert
  et al. 2010, mouse);
- inferred from sourced bounds, fitted to nothing: the lineage shares of daily output;
- inferred with no primary citation: the transit-amplifying transit time and division count, and the
  secretory progenitor residence;
- FITTED, and the only fitted quantity: `villus_residence_days`, solved so that the sum of every
  compartment's pool equals the published epithelial count. Nothing about any mutant is fitted, so
  the knockout pattern is the module's prediction and the count bound is only the fit closing.

The stem-cell pool is a consequence of the flows, not a claim: it is reported together with the
number of crypts it would imply at the measured 14 Lgr5 cells per crypt, and checked against nothing.
"""

from __future__ import annotations

import math

SRC_MILO = "Sender & Milo 2021, Nat Med 27:45"

# Measured, from the distilled Sender & Milo table (data/results/human_cell_turnover.json).
TOTAL_EPITHELIUM = 1.2e11  # small intestine epithelia in the reference adult
DAILY_OUTPUT = 2.2e10  # small intestine epithelia replaced per day
EPITHELIAL_LIFESPAN_DAYS = 5.3  # the published mean residence of all of them

STEP_DAYS = 0.25  # flows and growth are applied every quarter day, as in the haematopoiesis module
ISC_DIVISION_DAYS = 1.0  # Snippert et al. 2010, Cell 143:134: Lgr5 crypt base columnar cells divide daily
LGR5_PER_CRYPT = 14.0  # Snippert et al. 2010: about 14 Lgr5 stem cells per crypt (mouse)
TA_TRANSIT_DAYS = 2.5  # inferred, no primary citation fetched in this lane
TA_DIVISIONS = 5  # inferred, no primary citation fetched in this lane
SECRETORY_PROGENITOR_DAYS = 1.0  # inferred, no primary citation fetched in this lane

# The four mature lineages. `share` is the share of daily output; `residence` is the measured turnover
# time where one exists and None where the villus residence is fitted instead.
MATURE = {
    "Enterocyte": {
        "share": 0.95,
        "residence": None,
        "cl": "CL:0000584",
        "from": "TA",
        "requires": [],
        "share_evidence": "Cheng & Leblond 1974, Am J Anat 141:461: enterocytes are >95% of villus cells",
        "share_tier": "inferred",
        "residence_evidence": "fitted: the villus residence that makes the compartment pools sum to the "
        f"published {TOTAL_EPITHELIUM:.3g} cells of {SRC_MILO}",
    },
    "GobletCell": {
        "share": 0.04,
        "residence": 3.0,
        "cl": "CL:0000160",
        "from": "SecretoryProgenitor",
        "requires": ["GFI1", "SPDEF"],
        "share_evidence": "Karam 1999, Front Biosci 4:D286 (review): about 4% of epithelial cells in the "
        "duodenum",
        "share_tier": "inferred",
        "residence_evidence": "Cheng 1974, Am J Anat 141:481 (II. Mucous cells): turnover about 3 days in "
        "the mouse",
    },
    "EnteroendocrineCell": {
        "share": 0.008,
        "residence": None,
        "cl": "CL:0000164",
        "from": "SecretoryProgenitor",
        "requires": ["NEUROG3"],
        "share_evidence": "Cheng & Leblond 1974, Am J Anat 141:461: enteroendocrine cells are <1% of the "
        "epithelium",
        "share_tier": "inferred",
        "residence_evidence": "fitted with the enterocytes: no measured enteroendocrine residence was "
        "found, so the villus residence is used and the number is fitted, not measured",
    },
    "PanethCell": {
        "share": 0.002,
        "residence": 15.0,
        "cl": "CL:0000510",
        "from": "SecretoryProgenitor",
        "requires": ["GFI1", "SPDEF", "SOX9"],
        "share_evidence": "the remainder of the output once the three shares above are taken; no measured "
        "Paneth share of output was found",
        "share_tier": "inferred",
        "residence_evidence": "Cheng 1974, Am J Anat 141:521 (IV. Paneth cells): turnover about 15 days in "
        "the mouse",
    },
}

# The progenitor compartments, top-down. `requires` is the registered gate of each step.
COMPARTMENTS = [
    {
        "id": "ISC",
        "parent": None,
        "requires": ["TCF7L2", "ASCL2"],
        "residence": ISC_DIVISION_DAYS,
        "divisions": 0,
        "cl": "CL:0002250",
        "evidence": "Korinek et al. 1998, Nat Genet 19:379 (Tcf4: no proliferative compartment is "
        "maintained in the prospective crypts); van der Flier et al. 2009, Cell 136:903 (Ascl2: the Lgr5 "
        "stem cells disappear within days of an induced deletion)",
    },
    {
        "id": "TA",
        "parent": "ISC",
        "requires": ["RBPJ"],
        "residence": TA_TRANSIT_DAYS,
        "divisions": TA_DIVISIONS,
        "cl": None,  # no Cell Ontology id is stated here: none was looked up in this lane
        "evidence": "van Es et al. 2005, Nature 435:959 (conditional removal of the Notch effector "
        "RBP-J/CSL converts the proliferative crypt cells into post-mitotic goblet cells)",
    },
    {
        "id": "SecretoryProgenitor",
        "parent": "TA",
        "requires": ["ATOH1"],
        "residence": SECRETORY_PROGENITOR_DAYS,
        "divisions": 0,
        "cl": None,  # no Cell Ontology id is stated here: none was looked up in this lane
        "evidence": "Yang et al. 2001, Science 294:2155 (Math1/Atoh1: goblet, enteroendocrine and Paneth "
        "cells are depleted and enterocytes are unaffected, so the secretory branch has its own gate)",
    },
]


def plan() -> dict:
    """Solve the flows from the measured output upward, then fit the one villus residence.

    Every compartment carries `output` (cells leaving per day), `inflow` (what it needs from upstream
    per day, the output divided by its amplification), `pool` (the standing count the module designs
    for) and the per-step fractions the runtime applies. The amplifying compartment's pool uses the
    standing count of an exponentially amplifying transit, `output * tau / ln 2`, and not
    `output * transit`, which would put a third of the epithelium in the crypts.
    """
    out = {m: DAILY_OUTPUT * d["share"] for m, d in MATURE.items()}
    p: dict[str, dict] = {}
    for m in MATURE:
        p[m] = {"id": m, "output": out[m], "inflow": out[m], "amplification": 1, "mature": True}

    secretory = [m for m, d in MATURE.items() if d["from"] == "SecretoryProgenitor"]
    absorptive = [m for m, d in MATURE.items() if d["from"] == "TA"]

    sp = {"id": "SecretoryProgenitor", "output": sum(out[m] for m in secretory), "amplification": 1}
    sp["inflow"] = sp["output"] / sp["amplification"]
    sp["pool"] = sp["output"] * SECRETORY_PROGENITOR_DAYS
    p["SecretoryProgenitor"] = sp

    ta = {
        "id": "TA",
        "output": sum(out[m] for m in absorptive) + sp["inflow"],
        "amplification": 2**TA_DIVISIONS,
    }
    ta["inflow"] = ta["output"] / ta["amplification"]
    tau = TA_TRANSIT_DAYS / TA_DIVISIONS  # the doubling time implied by the transit and the divisions
    ta["pool"] = ta["output"] * tau / math.log(2.0)
    ta["doubling_days"] = tau
    p["TA"] = ta

    isc = {"id": "ISC", "output": ta["inflow"], "amplification": 1}
    isc["inflow"] = 0.0  # the root: self-renewal, not an inflow
    isc["pool"] = isc["output"] * ISC_DIVISION_DAYS
    isc["implied_crypts"] = isc["pool"] / LGR5_PER_CRYPT
    p["ISC"] = isc

    # The one fitted quantity: the villus residence of the lineages without a measured one, solved so
    # that every pool together equals the published epithelial count.
    fixed = isc["pool"] + ta["pool"] + sp["pool"]
    fixed += sum(out[m] * MATURE[m]["residence"] for m in MATURE if MATURE[m]["residence"] is not None)
    fitted_out = sum(out[m] for m in MATURE if MATURE[m]["residence"] is None)
    villus = (TOTAL_EPITHELIUM - fixed) / fitted_out
    if villus <= 0:
        raise ValueError(f"the fixed pools already exceed {TOTAL_EPITHELIUM:.3g}: no villus residence fits")
    for m, d in MATURE.items():
        residence = d["residence"] if d["residence"] is not None else villus
        p[m]["residence"] = residence
        p[m]["pool"] = out[m] * residence
        p[m]["loss_fraction"] = STEP_DAYS / residence

    # Per-step flows. The outflows are written as BioLang `share` lines, which take their part of the
    # population as it stood at the start of the step, so three outflows due together retain
    # 1 - s1 - s2 - s3 and not the product of their complements (docs/BIOLANG-v0.3.md, `share`). The
    # growth term is then solved for the discrete balance (1 + g) * retained = 1 - o / A per step,
    # which is zero for a compartment that only passes cells through.
    children = {
        "ISC": ["TA"],
        "TA": [*absorptive, "SecretoryProgenitor"],
        "SecretoryProgenitor": secretory,
    }
    for kids in children.values():
        total = sum(p[k]["inflow"] for k in kids) or 1.0
        for k in kids:
            p[k]["share"] = p[k]["inflow"] / total
    for cid, kids in children.items():
        c = p[cid]
        c["out_fraction"] = c["output"] / c["pool"]  # per day
        c["step_out"] = c["out_fraction"] * STEP_DAYS
        c["retained"] = 1.0 - sum(c["step_out"] * p[k]["share"] for k in kids)
        if cid == "ISC":
            c["growth_fraction"] = 1.0 / c["retained"] - 1.0  # the root: self-renewal balances the outflow
        else:
            c["growth_fraction"] = (1.0 - c["step_out"] / c["amplification"]) / c["retained"] - 1.0
    p["_fit"] = {
        "villus_residence_days": villus,
        "total_designed": sum(p[k]["pool"] for k in [*MATURE, *children]),
        "implied_mean_residence_days": TOTAL_EPITHELIUM / DAILY_OUTPUT,
        "published_mean_residence_days": EPITHELIAL_LIFESPAN_DAYS,
    }
    return p


def _req(factors: list[str]) -> str:
    return "".join(f", {f} = present" for f in sorted(factors))


def _fill_lines(
    cell_type: str, pool: float, retained: float, requires: list[str], growth: float
) -> list[str]:
    """Fill a compartment to its designed pool during the Fill stage, as the haematopoiesis module does:
    doubling to half, then 10% steps to nine tenths, then 1% steps, so the pool lands within 1%.

    A `count` condition is read at the start of a step, before that step's growth, while the standing
    count the run reports is read after it, so a compartment whose growth term is large would settle a
    factor of (1 + growth) above its designed pool. The caps are therefore the pre-growth count that
    leaves the reported pool where the module designed it: `pool / (1 + growth)`.
    """
    ev = 'evidence: inferred "a module of the adult epithelium starts at its steady pool"; confidence: 0.3 }'
    base = pool / (1.0 + growth)
    lines = []
    for tag, step, upto in (("2x", 1.0, 0.5), ("10pc", 0.1, 0.9), ("1pc", 0.01, 1.0)):
        lines.append(
            f"decision {cell_type}_fill_{tag} {{ action: divide; "
            + f"when: cell_type = {cell_type}{_req(requires)}, stage = Fill, count = <{base * upto:.4g}; "
            + f"fraction: {(1.0 + step) ** STEP_DAYS / retained - 1.0:.6g}; "
            + f"after: {STEP_DAYS:g} d; {ev}"
        )
    return lines


def to_bio(module_name: str = "organism.human.intestine") -> str:
    p = plan()
    fit = p["_fit"]
    factors = sorted(
        {f for c in COMPARTMENTS for f in c["requires"]} | {f for d in MATURE.values() for f in d["requires"]}
    )
    lines = [
        "# Generated by genomeos.organism.intestine: the small-intestinal crypt-villus lineages from one",
        "# stem compartment, gated by the factors of the closed sourced list in",
        "# scripts/intestine_registration.py, which was committed before any number below was fitted.",
        "# Measured and used as given: the epithelial count and daily output (Sender & Milo 2021), the",
        "# goblet and Paneth turnover times (Cheng 1974, mouse), the stem-cell division interval",
        "# (Snippert et al. 2010, mouse). Inferred from sourced bounds: the lineage shares of output.",
        "# Inferred with no primary citation: the transit-amplifying transit and division count, and the",
        "# secretory progenitor residence. FITTED, and the only fitted number: the villus residence,",
        f"# {fit['villus_residence_days']:.3g} d, solved so the pools sum to the published",
        f"# {TOTAL_EPITHELIUM:.3g} cells; the mouse-measured absorptive turnover is about 3 d and the",
        f"# published human mean residence is {EPITHELIAL_LIFESPAN_DAYS:g} d. The count asserts below are",
        "# therefore the closure of that fit and not a prediction; the knockouts in",
        "# intestine_mutants.bio are the prediction, because no mutant number is fitted.",
        f"module {module_name}",
        "import bio.std.development",
        "",
        "organism Intestine {",
        "  species: Homo sapiens",
        "  genome: GRCh38",
        "  root: ISC",
        "  resolution: populations",
        "  cell_type: ISC",
        f"  factors: {', '.join(factors)}",
        "  observe: count, fates",
        f"  assert: count at 1 yr in {TOTAL_EPITHELIUM * 0.85:.3g}..{TOTAL_EPITHELIUM * 1.15:.3g}",
    ]
    for m in MATURE:
        lines.append(f"  assert: type {m} at 1 yr in {p[m]['pool'] * 0.85:.3g}..{p[m]['pool'] * 1.15:.3g}")
    lines += [
        f'  evidence: curated "{SRC_MILO}; Cheng 1974, Am J Anat 141:461/481/521; Snippert et al. 2010, '
        + 'Cell 143:134"; confidence: 0.4',
        "}",
        "",
        'stage Fill { from: 0 d; to: 90 d; evidence: inferred "every compartment fills to its adult size"; '
        + "confidence: 0.3 }",
        'stage Steady { from: 90 d; evidence: curated "adult steady state: production balances loss"; '
        + "confidence: 0.7 }",
        "",
    ]
    for c in COMPARTMENTS:
        parent = "Blastomere" if c["parent"] is None else c["parent"]
        onto = f"ontology: {c['cl']}; " if c["cl"] else ""
        lines.append(
            f"cell_type {c['id']} {{ parent: {parent}; {onto}"
            + f'evidence: experimental "{c["evidence"]}"; confidence: 0.8 }}'
        )
    for m, d in MATURE.items():
        lines.append(
            f"cell_type {m} {{ parent: PostMitotic; ontology: {d['cl']}; "
            + 'evidence: curated "Cell Ontology"; confidence: 0.9 }'
        )
    isc = p["ISC"]
    lines += [
        "",
        f"# the stem compartment: {isc['output']:.3g} cells leave per day and each one leaves every",
        f"# {ISC_DIVISION_DAYS:g} day (Snippert et al. 2010, mouse), so the pool is {isc['pool']:.3g} cells.",
        "# That pool is a consequence of the flows and is claimed of nothing: at the measured",
        f"# {LGR5_PER_CRYPT:g} Lgr5 cells per crypt it would imply {isc['implied_crypts']:.3g} crypts,",
        "# a number this module reports and does not check against a source.",
        *_fill_lines(
            "ISC", isc["pool"], isc["retained"], COMPARTMENTS[0]["requires"], isc["growth_fraction"]
        ),
        f"decision isc_renew {{ action: divide; when: cell_type = ISC{_req(COMPARTMENTS[0]['requires'])}, "
        + f"count = >={isc['pool'] * 0.99 / (1.0 + isc['growth_fraction']):.3g}; "
        + f"fraction: {isc['growth_fraction']:.6g}; "
        + f"after: {STEP_DAYS:g} d; "
        + 'evidence: inferred "self-renewal balances the outflow: one division per '
        + f'{ISC_DIVISION_DAYS:g} day per stem cell"; confidence: 0.3 }}',
        "",
    ]
    # the progenitor steps, each gated by the union of the source's and the target's requirements
    for c in COMPARTMENTS:
        if c["parent"] is None:
            continue
        parent_requires = next(x["requires"] for x in COMPARTMENTS if x["id"] == c["parent"])
        req = _req(sorted(set(parent_requires) | set(c["requires"])))
        frac = p[c["parent"]]["step_out"] * p[c["id"]]["share"]
        lines.append(
            f"decision {c['parent']}_to_{c['id']} {{ action: differentiate; "
            + f"when: cell_type = {c['parent']}{req}; to: {c['id']}; share: {frac:.6g}; "
            + f'after: {STEP_DAYS:g} d; evidence: experimental "{c["evidence"]}"; confidence: 0.7 }}'
        )
    lines += [
        "",
        f"# amplification inside the transit-amplifying compartment ({TA_DIVISIONS} inferred divisions over",
        f"# an inferred {TA_TRANSIT_DAYS:g} d transit, a doubling every {p['TA']['doubling_days']:.3g} d),",
        "# capped at three times the designed pool so a blocked outflow backs up instead of growing",
        "# without bound",
    ]
    for c in COMPARTMENTS:
        if c["parent"] is None:
            continue
        pc = p[c["id"]]
        lines.extend(_fill_lines(c["id"], pc["pool"], pc["retained"], c["requires"], pc["growth_fraction"]))
        if c["divisions"]:
            lines.append(
                f"decision {c['id']}_amplify {{ action: divide; "
                + f"when: cell_type = {c['id']}{_req(c['requires'])}, "
                + f"count = <{pc['pool'] * 3 / (1.0 + pc['growth_fraction']):.3g}; "
                + f"fraction: {pc['growth_fraction']:.6g}; after: {STEP_DAYS:g} d; "
                + f'evidence: inferred "{c["divisions"]} divisions per transit of {c["residence"]:g} days; '
                + f'pool {pc["pool"]:.3g}"; confidence: 0.3 }}'
            )
    lines.append("")
    # the mature outputs, each gated by the registered factors of that lineage
    for m, d in MATURE.items():
        source = d["from"]
        req = _req(
            sorted(set(next(x["requires"] for x in COMPARTMENTS if x["id"] == source)) | set(d["requires"]))
        )
        frac = p[source]["step_out"] * p[m]["share"]
        ev = d["share_evidence"]
        lines.append(
            f"decision {source}_to_{m} {{ action: differentiate; when: cell_type = {source}{req}; "
            + f"to: {m}; share: {frac:.6g}; after: {STEP_DAYS:g} d; "
            + f'evidence: {d["share_tier"]} "share {d["share"]:g} of {DAILY_OUTPUT:.3g} cells per day '
            + f'({SRC_MILO}): {ev}"; confidence: 0.3 }}'
        )
    lines.append("")
    for m, d in MATURE.items():
        tier = "experimental" if d["residence"] is not None else "inferred"
        conf = "0.6" if d["residence"] is not None else "0.2"
        lines.append(
            f"decision {m}_loss {{ action: die; when: cell_type = {m}; "
            + f"fraction: {p[m]['loss_fraction']:.6g}; after: {STEP_DAYS:g} d; "
            + f'evidence: {tier} "residence {p[m]["residence"]:.3g} d: {d["residence_evidence"]}"; '
            + f"confidence: {conf} }}"
        )
    return "\n".join(lines) + "\n"


# The knockouts, exactly as scripts/intestine_registration.py registered them before this module's
# numbers were fitted. `expect` carries the published phenotype, including the part a module made of
# flows cannot express, which the registration named in advance.
MUTANTS = [
    (
        "tcf7l2",
        "TCF7L2",
        "Korinek et al. 1998, Nat Genet 19:379",
        "no proliferative compartment in the prospective crypts, so no lineage is produced",
        ["Enterocyte", "GobletCell", "PanethCell", "EnteroendocrineCell"],
        [],
    ),
    (
        "ascl2",
        "ASCL2",
        "van der Flier et al. 2009, Cell 136:903",
        "the stem cells disappear within days of the deletion, so no lineage is produced",
        ["Enterocyte", "GobletCell", "PanethCell", "EnteroendocrineCell"],
        [],
    ),
    (
        "rbpj",
        "RBPJ",
        "van Es et al. 2005, Nature 435:959",
        "the proliferative compartment is lost (published as a conversion of those cells into goblet "
        "cells, which a module made of flows cannot express: registered partial in advance)",
        ["Enterocyte", "GobletCell", "PanethCell", "EnteroendocrineCell"],
        [],
    ),
    (
        "atoh1",
        "ATOH1",
        "Yang et al. 2001, Science 294:2155",
        "goblet, enteroendocrine and Paneth cells are depleted; enterocytes are unaffected",
        ["GobletCell", "PanethCell", "EnteroendocrineCell"],
        ["Enterocyte"],
    ),
    (
        "neurog3",
        "NEUROG3",
        "Jenny et al. 2002, EMBO J 21:6338",
        "no intestinal endocrine cells; the other epithelial types differentiate properly",
        ["EnteroendocrineCell"],
        ["Enterocyte", "GobletCell", "PanethCell"],
    ),
    (
        "gfi1",
        "GFI1",
        "Shroyer et al. 2005, Genes Dev 19:2412",
        "no Paneth cells and fewer goblet cells, with supernumerary enteroendocrine cells (the module "
        "can only delete a lineage, so 'fewer goblet' is overstated: registered partial in advance)",
        ["GobletCell", "PanethCell"],
        ["Enterocyte", "EnteroendocrineCell"],
    ),
    (
        "spdef",
        "SPDEF",
        "Gregorieff et al. 2009, Gastroenterology 137:1333",
        "mature goblet and Paneth cells are lost and immature secretory progenitors accumulate",
        ["GobletCell", "PanethCell"],
        ["Enterocyte", "EnteroendocrineCell"],
    ),
    (
        "sox9",
        "SOX9",
        "Mori-Akiyama et al. 2007, Gastroenterology 133:539; Bastide et al. 2007, J Cell Biol 178:635",
        "Paneth cells are not formed; the other epithelial types are unaffected",
        ["PanethCell"],
        ["Enterocyte", "GobletCell", "EnteroendocrineCell"],
    ),
]
KEPT_FLOOR = 1.0e6  # the floor the registration set for a population a knockout does not touch


def mutants_bio(module_name: str = "organism.human.intestine_mutants") -> str:
    lines = [
        "# The eight knockouts of the closed sourced list in scripts/intestine_registration.py, against",
        "# the published (mouse) phenotypes: each removes one factor and states what must disappear and",
        "# what must remain. Two are registered PARTIAL in advance, because their published phenotype is a",
        "# re-allocation and not a deletion. Run:",
        "# genomeos grow data/organisms/human/intestine_mutants.bio --experiments",
        f"module {module_name}",
        "import intestine.bio",
        "",
    ]
    for name, factor, source, expect, lost, kept in MUTANTS:
        asserts = "".join(f" assert: type {t} at 1 yr = 0;" for t in lost)
        asserts += "".join(f" assert: type {t} at 1 yr >= {KEPT_FLOOR:.3g};" for t in kept)
        lines.append(
            f'experiment {name} {{ knockout: {factor}; until: 1 yr; expect: "{expect}";{asserts} '
            + f'evidence: experimental "{source}"; confidence: 0.8 }}'
        )
    return "\n".join(lines) + "\n"
