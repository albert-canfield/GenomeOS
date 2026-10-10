# SPDX-License-Identifier: AGPL-3.0-or-later
"""The colonic crypt as a mechanism module: one stem compartment to three mature lineages.

Area E's third human mechanism module, after haematopoiesis and the small-intestinal crypt. The gates
come from the closed, sourced list registered in scripts/colon_registration.py BEFORE any number here
was fitted, and the written comparison of the colon with the small intestine is the dated colon
section of docs/AREA-E-NEXT.md. Seven factors, each with a loss-of-function phenotype fetched FOR THE
COLON, and seven more excluded with the reason - ASCL2 and SPDEF among them, because they are gates
of the small-intestinal module with no colon sentence to quote.

What differs from the small-intestinal module, and why the code differs with it:

- three mature lineages, not four. The colon has no Paneth cells, and its Paneth equivalent (Reg4+
  deep crypt secretory cells) is excluded for want of a quantity, so there is no long-lived resident
  lineage at the crypt base;
- no villus. Every mature colonic cell rides the same crypt-to-surface conveyor and is shed at the
  surface, so ONE residence is shared by all three and that residence is the single fitted number;
- the goblet lineage is 20% of the epithelium (Katz et al. 2002, measured) and not about 4%, which
  makes every consequence of deleting it five times larger;
- KLF4 is a gate here. The small-intestinal registration excluded it for being colon-specific.

What is measured, what is inferred and what is fitted is kept apart:

- measured and used as given: the colon's epithelial count, daily output and lifespan (Sender & Milo
  2021), the 20% goblet share of colonic epithelial cells (Katz et al. 2002, mouse), about 20 Lgr5+
  stem cells per colonic crypt (Sasaki et al. 2016, mouse);
- inferred from the measured share: the lineage shares of daily output, which equal the standing
  shares because the three mature residences are one number;
- inferred with no primary citation: the stem division interval, the transit-amplifying transit and
  division count, and the secretory progenitor residence;
- FITTED, and the only fitted quantity: `surface_residence_days`, solved so that the sum of every
  compartment's pool equals the published colonic epithelial count. Nothing about any mutant is
  fitted, so the knockout pattern is what the module predicts and the count asserts are the fit
  closing.

Two of the seven gates are CONTESTED in the literature (TCF7L2 and KLF4) and the registration says so
in advance: those knockouts are counted as agreeing with one named paper, never as reproducing the
colon.

The step machinery is imported from the small-intestinal generator rather than copied, so that the
two programs cannot drift apart in how a `share` outflow and a `divide` growth term are balanced.
"""

from __future__ import annotations

import math

from genomeos.organism.intestine import STEP_DAYS, _fill_lines, _req

SRC_MILO = "Sender & Milo 2021, Nat Med 27:45"
SRC_KATZ = "Katz et al. 2002, Development 129:2619"
SRC_SASAKI = "Sasaki et al. 2016, Proc Natl Acad Sci USA 113:E5399"

# Measured, from the distilled Sender & Milo table (data/results/human_cell_turnover.json), row
# "Colon epithelia".
TOTAL_EPITHELIUM = 6.0e10  # colon epithelia in the reference adult
DAILY_OUTPUT = 1.7e10  # colon epithelia replaced per day
EPITHELIAL_LIFESPAN_DAYS = 3.4  # the published mean residence of all of them

LGR5_PER_CRYPT = 20.0  # Sasaki et al. 2016: about 20 Lgr5+ stem cells per colonic crypt (mouse)
STEM_DIVISION_DAYS = 1.0  # inferred, no primary colon citation fetched in this lane
TA_TRANSIT_DAYS = 2.5  # inferred, no primary citation fetched in this lane
TA_DIVISIONS = 5  # inferred, no primary citation fetched in this lane
SECRETORY_PROGENITOR_DAYS = 1.0  # inferred, no primary citation fetched in this lane

GOBLET_SHARE = 0.20  # measured as a standing share (Katz et al. 2002), used as the output share
EEC_SHARE = 0.004  # inferred inside the <1% bound of the small-intestinal registration

# The three mature lineages. `residence` is None for every one of them: the colon has no long-lived
# resident type, so the one fitted surface residence is shared.
MATURE = {
    "Colonocyte": {
        "share": 1.0 - GOBLET_SHARE - EEC_SHARE,
        "residence": None,
        "cl": "CL:1000347",
        "from": "ColonTA",
        "requires": [],
        "share_evidence": "the remainder once the measured goblet share and the inferred enteroendocrine "
        "share are taken",
        "share_tier": "inferred",
        "residence_evidence": "fitted: the surface residence that makes the compartment pools sum to the "
        f"published {TOTAL_EPITHELIUM:.3g} cells of {SRC_MILO}",
    },
    "ColonGobletCell": {
        "share": GOBLET_SHARE,
        "residence": None,
        "cl": "CL:0009039",
        "from": "ColonSecretoryProgenitor",
        "requires": ["GFI1", "KLF4", "SOX9"],
        "share_evidence": f"{SRC_KATZ}: in wild-type colons 20% of epithelial cells were mature goblet "
        "cells; used as the share of daily output, which follows because all three mature residences are "
        "the one fitted number",
        "share_tier": "experimental",
        "residence_evidence": "fitted with the colonocytes: no measured colonic goblet residence was found, "
        "and the colon has no resident lineage, so the surface residence is shared and the number is fitted",
    },
    "ColonEnteroendocrineCell": {
        "share": EEC_SHARE,
        "residence": None,
        "cl": "CL:0009042",
        "from": "ColonSecretoryProgenitor",
        "requires": ["NEUROG3"],
        "share_evidence": "inside the <1% of the epithelium bound registered by the small-intestinal "
        "module (scripts/intestine_registration.py, citing Cheng & Leblond 1974 for stomach, small "
        "intestine and colon); fitted to nothing",
        "share_tier": "inferred",
        "residence_evidence": "fitted with the colonocytes: no measured colonic enteroendocrine residence "
        "was found",
    },
}

# The progenitor compartments, top-down. `requires` is the registered gate of each step.
COMPARTMENTS = [
    {
        "id": "ColonStem",
        "parent": None,
        "requires": ["TCF7L2"],
        "residence": STEM_DIVISION_DAYS,
        "divisions": 0,
        "cl": "CL:0009043",
        "evidence": "van Es et al. 2012, Mol Cell Biol 32:1918 (an adult inducible deletion of Tcf4 "
        "abolished the proliferating Ki67+ cells of the colon within 7 days); REGISTERED CONTESTED: "
        "Angus-Hill et al. 2011, PNAS 108:4914 report adult colonic crypts lacking Tcf4 enlarged, with "
        "proliferative and differentiated cells present",
    },
    {
        "id": "ColonTA",
        "parent": "ColonStem",
        "requires": ["RBPJ"],
        "residence": TA_TRANSIT_DAYS,
        "divisions": TA_DIVISIONS,
        "cl": "CL:0009011",
        "evidence": "Pellegrinet et al. 2011, Gastroenterology 140:1230 (in the colon as in the proximal "
        "intestine, the RBP-J mutant loses the proliferative crypt compartment completely and converts it "
        "into goblet cells; the conversion is not expressible here, registered partial in advance)",
    },
    {
        "id": "ColonSecretoryProgenitor",
        "parent": "ColonTA",
        "requires": ["ATOH1"],
        "residence": SECRETORY_PROGENITOR_DAYS,
        "divisions": 0,
        "cl": None,  # no Cell Ontology id is stated here: none was looked up in this lane
        "evidence": "Lo et al. 2016, Cell Mol Gastroenterol Hepatol 3:51 (an adult inducible deletion of "
        "Atoh1 leaves the secretory lineages nearly absent in the colon, the absorptive colonocytes "
        "remaining, so the secretory branch has its own gate)",
    },
]


def plan() -> dict:
    """Solve the flows from the measured output upward, then fit the one surface residence.

    The structure is the small-intestinal module's, with the four-lineage villus replaced by a
    three-lineage surface: every compartment carries `output` (cells leaving per day), `inflow` (what
    it needs from upstream, the output divided by its amplification), `pool` (the standing count the
    module designs for) and the per-step fractions the runtime applies. The amplifying compartment's
    pool uses the standing count of an exponentially amplifying transit, `output * tau / ln 2`.
    """
    out = {m: DAILY_OUTPUT * d["share"] for m, d in MATURE.items()}
    p: dict[str, dict] = {}
    for m in MATURE:
        p[m] = {"id": m, "output": out[m], "inflow": out[m], "amplification": 1, "mature": True}

    secretory = [m for m, d in MATURE.items() if d["from"] == "ColonSecretoryProgenitor"]
    absorptive = [m for m, d in MATURE.items() if d["from"] == "ColonTA"]

    sp = {
        "id": "ColonSecretoryProgenitor",
        "output": sum(out[m] for m in secretory),
        "amplification": 1,
    }
    sp["inflow"] = sp["output"] / sp["amplification"]
    sp["pool"] = sp["output"] * SECRETORY_PROGENITOR_DAYS
    p["ColonSecretoryProgenitor"] = sp

    ta = {
        "id": "ColonTA",
        "output": sum(out[m] for m in absorptive) + sp["inflow"],
        "amplification": 2**TA_DIVISIONS,
    }
    ta["inflow"] = ta["output"] / ta["amplification"]
    tau = TA_TRANSIT_DAYS / TA_DIVISIONS  # the doubling time implied by the transit and the divisions
    ta["pool"] = ta["output"] * tau / math.log(2.0)
    ta["doubling_days"] = tau
    p["ColonTA"] = ta

    stem = {"id": "ColonStem", "output": ta["inflow"], "amplification": 1}
    stem["inflow"] = 0.0  # the root: self-renewal, not an inflow
    stem["pool"] = stem["output"] * STEM_DIVISION_DAYS
    stem["implied_crypts"] = stem["pool"] / LGR5_PER_CRYPT
    p["ColonStem"] = stem

    # The one fitted quantity: the surface residence, shared by all three mature lineages because the
    # colon has no resident type, solved so that every pool together equals the published count.
    fixed = stem["pool"] + ta["pool"] + sp["pool"]
    fixed += sum(out[m] * MATURE[m]["residence"] for m in MATURE if MATURE[m]["residence"] is not None)
    fitted_out = sum(out[m] for m in MATURE if MATURE[m]["residence"] is None)
    surface = (TOTAL_EPITHELIUM - fixed) / fitted_out
    if surface <= 0:
        raise ValueError(f"the fixed pools already exceed {TOTAL_EPITHELIUM:.3g}: no surface residence fits")
    for m, d in MATURE.items():
        residence = d["residence"] if d["residence"] is not None else surface
        p[m]["residence"] = residence
        p[m]["pool"] = out[m] * residence
        p[m]["loss_fraction"] = STEP_DAYS / residence

    children = {
        "ColonStem": ["ColonTA"],
        "ColonTA": [*absorptive, "ColonSecretoryProgenitor"],
        "ColonSecretoryProgenitor": secretory,
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
        if cid == "ColonStem":
            c["growth_fraction"] = 1.0 / c["retained"] - 1.0  # the root: self-renewal balances the outflow
        else:
            c["growth_fraction"] = (1.0 - c["step_out"] / c["amplification"]) / c["retained"] - 1.0
    mature_pool = sum(p[m]["pool"] for m in MATURE)
    p["_fit"] = {
        "surface_residence_days": surface,
        "total_designed": sum(p[k]["pool"] for k in [*MATURE, *children]),
        "implied_mean_residence_days": TOTAL_EPITHELIUM / DAILY_OUTPUT,
        "published_mean_residence_days": EPITHELIAL_LIFESPAN_DAYS,
        # The measured 20% goblet share returns as a share of the MATURE cells by construction, and is
        # not a check (scripts/colon_registration.py, NOT_CLAIMED). As a share of every epithelial cell,
        # crypt compartments included, it is smaller, and Katz et al. counted total epithelial cells in
        # sections, so neither figure may be reported as a match.
        "goblet_share_of_mature": p["ColonGobletCell"]["pool"] / mature_pool,
        "goblet_share_of_all": p["ColonGobletCell"]["pool"] / TOTAL_EPITHELIUM,
        # Registered reading R2, computed from the shares alone before the run.
        "eec_fold_if_goblet_blocked": 1.0 / (p["ColonEnteroendocrineCell"]["share"]),
    }
    return p


def to_bio(module_name: str = "organism.human.colon") -> str:
    p = plan()
    fit = p["_fit"]
    factors = sorted(
        {f for c in COMPARTMENTS for f in c["requires"]} | {f for d in MATURE.values() for f in d["requires"]}
    )
    lines = [
        "# Generated by genomeos.organism.colon: the colonic crypt's three mature lineages from one stem",
        "# compartment, gated by the factors of the closed sourced list in scripts/colon_registration.py,",
        "# which was committed before any number below was fitted. Each of its seven gates rests on a",
        "# loss-of-function phenotype fetched FOR THE COLON; ASCL2 and SPDEF, gates of the small-intestinal",
        "# module, are excluded here for want of one, and the colon's Paneth equivalent (Reg4+ deep crypt",
        "# secretory cells) is excluded for want of a quantity, so there are three mature lineages and no",
        "# Paneth compartment. Two gates, TCF7L2 and KLF4, are REGISTERED CONTESTED: the literature",
        "# disagrees about the colon and this module can agree with one paper only.",
        "# Measured and used as given: the colonic count, daily output and lifespan (Sender & Milo 2021),",
        "# the 20% goblet share of colonic epithelial cells (Katz et al. 2002, mouse), about 20 Lgr5+ cells",
        "# per colonic crypt (Sasaki et al. 2016, mouse). Inferred with no primary citation: the stem",
        "# division interval, the transit-amplifying transit and division count, the secretory progenitor",
        "# residence. FITTED, and the only fitted number: the surface residence,",
        f"# {fit['surface_residence_days']:.3g} d, shared by all three mature lineages and solved so the",
        f"# pools sum to the published {TOTAL_EPITHELIUM:.3g} cells; the published mean residence is",
        f"# {EPITHELIAL_LIFESPAN_DAYS:g} d and the two measured numbers imply",
        f"# {fit['implied_mean_residence_days']:.3g} d on their own. The count asserts below are therefore",
        "# the closure of that fit and not a prediction; the knockouts in colon_mutants.bio are the",
        "# prediction,",
        "# because no mutant number is fitted. The measured 20% goblet share returns as a share of the",
        f"# mature cells ({fit['goblet_share_of_mature']:.3g}) BY CONSTRUCTION and is not a check.",
        f"module {module_name}",
        "import bio.std.development",
        "",
        "organism Colon {",
        "  species: Homo sapiens",
        "  genome: GRCh38",
        "  root: ColonStem",
        "  resolution: populations",
        "  cell_type: ColonStem",
        f"  factors: {', '.join(factors)}",
        "  observe: count, fates",
        f"  assert: count at 1 yr in {TOTAL_EPITHELIUM * 0.85:.3g}..{TOTAL_EPITHELIUM * 1.15:.3g}",
    ]
    for m in MATURE:
        lines.append(f"  assert: type {m} at 1 yr in {p[m]['pool'] * 0.85:.3g}..{p[m]['pool'] * 1.15:.3g}")
    lines += [
        f'  evidence: curated "{SRC_MILO}; {SRC_KATZ}; {SRC_SASAKI}"; confidence: 0.4',
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
    stem = p["ColonStem"]
    lines += [
        "",
        f"# the stem compartment: {stem['output']:.3g} cells leave per day and each one leaves every",
        f"# {STEM_DIVISION_DAYS:g} day (inferred, no colon citation fetched), so the pool is",
        f"# {stem['pool']:.3g} cells. That pool is a consequence of the flows and is claimed of nothing:",
        "# at the measured",
        f"# {LGR5_PER_CRYPT:g} Lgr5+ cells per colonic crypt (Sasaki et al. 2016) it would imply",
        f"# {stem['implied_crypts']:.3g} crypts, a number this module reports and checks against nothing.",
        *_fill_lines(
            "ColonStem", stem["pool"], stem["retained"], COMPARTMENTS[0]["requires"], stem["growth_fraction"]
        ),
        "decision colonstem_renew { action: divide; "
        + f"when: cell_type = ColonStem{_req(COMPARTMENTS[0]['requires'])}, "
        + f"count = >={stem['pool'] * 0.99 / (1.0 + stem['growth_fraction']):.3g}; "
        + f"fraction: {stem['growth_fraction']:.6g}; "
        + f"after: {STEP_DAYS:g} d; "
        + 'evidence: inferred "self-renewal balances the outflow: one division per '
        + f'{STEM_DIVISION_DAYS:g} day per stem cell"; confidence: 0.3 }}',
        "",
    ]
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
        f"# an inferred {TA_TRANSIT_DAYS:g} d transit, a doubling every",
        f"# {p['ColonTA']['doubling_days']:.3g} d),",
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
    for m, d in MATURE.items():
        source = d["from"]
        req = _req(
            sorted(set(next(x["requires"] for x in COMPARTMENTS if x["id"] == source)) | set(d["requires"]))
        )
        frac = p[source]["step_out"] * p[m]["share"]
        lines.append(
            f"decision {source}_to_{m} {{ action: differentiate; when: cell_type = {source}{req}; "
            + f"to: {m}; share: {frac:.6g}; after: {STEP_DAYS:g} d; "
            + f'evidence: {d["share_tier"]} "share {d["share"]:g} of {DAILY_OUTPUT:.3g} cells per day '
            + f'({SRC_MILO}): {d["share_evidence"]}"; confidence: 0.3 }}'
        )
    lines.append("")
    for m, d in MATURE.items():
        lines.append(
            f"decision {m}_loss {{ action: die; when: cell_type = {m}; "
            + f"fraction: {p[m]['loss_fraction']:.6g}; after: {STEP_DAYS:g} d; "
            + f'evidence: inferred "residence {p[m]["residence"]:.3g} d: {d["residence_evidence"]}"; '
            + "confidence: 0.2 }"
        )
    return "\n".join(lines) + "\n"


# The knockouts, exactly as scripts/colon_registration.py registered them before this module's
# numbers were fitted. `expect` carries the published colon phenotype, including the part a module
# made of flows cannot express and, for the two contested gates, the paper the module disagrees with.
MUTANTS = [
    (
        "tcf7l2",
        "TCF7L2",
        "van Es et al. 2012, Mol Cell Biol 32:1918",
        "the proliferating cells of the colon are abolished, so no lineage is produced (REGISTERED "
        "CONTESTED in advance: Angus-Hill et al. 2011, PNAS 108:4914 instead report adult colonic crypts "
        "lacking Tcf4 enlarged, with proliferative and differentiated cells present)",
        ["Colonocyte", "ColonGobletCell", "ColonEnteroendocrineCell"],
        [],
    ),
    (
        "rbpj",
        "RBPJ",
        "Pellegrinet et al. 2011, Gastroenterology 140:1230",
        "the proliferative crypt compartment of the colon is lost completely (published as a conversion of "
        "those cells into goblet cells, which a module made of flows cannot express: registered partial in "
        "advance)",
        ["Colonocyte", "ColonGobletCell", "ColonEnteroendocrineCell"],
        [],
    ),
    (
        "atoh1",
        "ATOH1",
        "Lo et al. 2016, Cell Mol Gastroenterol Hepatol 3:51",
        "the secretory lineages are nearly absent in the colon; the absorptive colonocytes remain",
        ["ColonGobletCell", "ColonEnteroendocrineCell"],
        ["Colonocyte"],
    ),
    (
        "neurog3",
        "NEUROG3",
        "Mellitzer et al. 2010, J Clin Invest 120:1708; Jenny et al. 2002, EMBO J 21:6338",
        "no enteroendocrine cells anywhere along the intestinal tract, the colon included; the colonic "
        "goblet cells remain",
        ["ColonEnteroendocrineCell"],
        ["Colonocyte", "ColonGobletCell"],
    ),
    (
        "gfi1",
        "GFI1",
        "Shroyer et al. 2005, Genes Dev 19:2412",
        "about a third the normal number of goblet cells, the reduction also seen in the colon, with more "
        "enteroendocrine cells (colonic PYY sixfold up); a presence-or-absence gate can only delete the "
        "goblet lineage, so the reduction is overstated: registered partial in advance",
        ["ColonGobletCell"],
        ["Colonocyte", "ColonEnteroendocrineCell"],
    ),
    (
        "sox9",
        "SOX9",
        "Bastide et al. 2007, J Cell Biol 178:635",
        "a strong reduction of the colonic goblet population, the enteroendocrine and absorptive "
        "populations comparable to the control (the reduction is stated as a deletion: registered partial "
        "in advance)",
        ["ColonGobletCell"],
        ["Colonocyte", "ColonEnteroendocrineCell"],
    ),
    (
        "klf4",
        "KLF4",
        "Katz et al. 2002, Development 129:2619",
        "90% fewer colonic goblet cells in the newborn, the colonocytes and enteroendocrine cells normal "
        "(REGISTERED CONTESTED in advance: Pellegrinet et al. 2011 found an adult inducible deletion "
        "changed no goblet cell number in either the small intestine or the colon, and this module has no "
        "age axis to hold both)",
        ["ColonGobletCell"],
        ["Colonocyte", "ColonEnteroendocrineCell"],
    ),
]
KEPT_FLOOR = 1.0e6  # the floor the registration set for a population a knockout does not touch


def mutants_bio(module_name: str = "organism.human.colon_mutants") -> str:
    lines = [
        "# The seven knockouts of the closed sourced list in scripts/colon_registration.py, against the",
        "# published COLON phenotypes: each removes one factor and states what must disappear and what must",
        "# remain. Three are registered PARTIAL in advance, because their published phenotype is a",
        "# reduction or a conversion and not a deletion, and TWO are registered CONTESTED, because the",
        "# literature disagrees about the colon: tcf7l2 and klf4 can agree with one paper only and are not",
        "# counted as reproducing the colon. Run:",
        "# genomeos grow data/organisms/human/colon_mutants.bio --experiments",
        f"module {module_name}",
        "import colon.bio",
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
