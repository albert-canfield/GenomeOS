# SPDX-License-Identifier: AGPL-3.0-or-later
"""The gastric gland as a mechanism module: one isthmus compartment to four mature lineages.

Area E's fourth human mechanism module, after haematopoiesis, the small-intestinal crypt and the
colonic crypt. The gates come from the closed, sourced list registered in
scripts/stomach_registration.py BEFORE any number here was fitted, and the written comparison of the
stomach with both gut modules is the dated stomach section of docs/AREA-E-NEXT.md. Five factors, each
with a loss-of-function phenotype fetched FOR THE STOMACH, and ten excluded with the reason - ATOH1,
the secretory gate of both gut modules, among them, for want of any stomach sentence.

What differs from the gut modules, and why the code differs with it:

- the proliferative zone is in the MIDDLE of the unit. Karam & Leblond 1993 put the division in the
  isthmus, from which cells migrate outward to the pit and inward to the neck and base, so the root
  of this module is also its amplifying compartment and carries the cap the colonic
  transit-amplifying compartment carried. There is no separate transit-amplifying step;
- the four residences differ by a factor of 60 (pit 3.1 d and parietal 54 d measured, endocrine
  86.6 d from a fetched half-life, chief 190 d inferred), where the colon had ONE shared residence.
  So the standing composition and the daily-output composition are not the same thing here, and the
  module is solved from the measured standing composition rather than from assumed output shares;
- the mucous neck cell is an intermediate on the chief route with exactly ONE outflow, where the
  colonic secretory progenitor had two. Blocking the chief route therefore has nowhere to send the
  flow, which is registered reading R4;
- every mature lineage gets Fill decisions, which the gut modules needed only for their progenitors:
  a lineage that lives 190 days cannot reach its adult pool inside a 90-day Fill stage by inflow.

What is measured, what is inferred and what is fitted is kept apart:

- measured and used as given: the stomach's epithelial count and daily output (Sender & Milo 2021);
  the standing count per gastric unit of pit, neck, zymogenic and parietal cells and the pit and
  parietal residences (Karam & Leblond 1993, parts I-IV, mouse corpus); the gastric endocrine cell's
  half-life of about 2 months (Hibdon et al. 2023);
- inferred: that a mouse corpus composition can stand for a human stomach that also has an antrum;
  the neck residence as the midpoint of the fetched 1-to-2-week range; the endocrine residence as the
  mean the fetched half-life implies under exponential loss; the chief residence, for which NO number
  was fetched in this lane; the endocrine standing count, inside the carried <1% bound;
- FITTED, and the only fitted quantity: `isthmus_residence_days`, solved so that the sum of every
  compartment's pool equals the published stomach epithelial count. It is held against a bound FETCHED
  BEFORE THE FIT: less than 2.6 days for the isthmal granule-free cells. Nothing about any mutant is
  fitted.

One of the five gates is registered full (ESRRG), three partial and one CONTESTED (RBPJ, whose paper's
own genetic and pharmacological arms differ in severity), so this module can never report five of five.

The step machinery is imported from the small-intestinal generator rather than copied, so that the
programs cannot drift apart in how a `share` outflow and a `divide` growth term are balanced.
"""

from __future__ import annotations

import math

from genomeos.organism.intestine import STEP_DAYS, _fill_lines, _req

SRC_MILO = "Sender & Milo 2021, Nat Med 27:45"
SRC_KARAM = "Karam & Leblond 1993, Anat Rec 236:259-340 (parts I-V, mouse corpus)"
SRC_HIBDON = "Hibdon et al. 2023, Am J Physiol Gastrointest Liver Physiol 325:G458"

# Measured, from the distilled Sender & Milo table (data/results/human_cell_turnover.json), row
# "Stomach epithelia".
TOTAL_EPITHELIUM = 1.7e10  # stomach epithelia in the reference adult
DAILY_OUTPUT = 1.3e9  # stomach epithelia replaced per day
EPITHELIAL_LIFESPAN_DAYS = 13.0  # the published mean residence; it is the ratio of the two above

# Measured per gastric unit in the mouse corpus, used as the standing COMPOSITION of the module.
NECK_UNIT_COUNT = 12.6  # Karam & Leblond part III
NECK_RESIDENCE_DAYS = 10.5  # inferred: the midpoint of the fetched "1-2 weeks"
CHIEF_RESIDENCE_DAYS = 190.0  # inferred, NO number fetched in this lane; see READINGS R5
ENDOCRINE_HALF_LIFE_DAYS = 60.0  # Hibdon et al. 2023: a half-life of about 2 months
# The mean residence that half-life implies under the exponential loss law the module applies.
ENDOCRINE_RESIDENCE_DAYS = ENDOCRINE_HALF_LIFE_DAYS / math.log(2.0)
GRANULE_FREE_TURNOVER_BOUND_DAYS = 2.6  # Karam & Leblond part I, an UPPER bound: what the fit is held to

# The four mature lineages. `unit_count` is the measured standing count per gastric unit (the endocrine
# one inferred), `residence` the mean residence in days. Both are given, so each lineage's share of the
# daily output follows as unit_count / residence and is not assumed.
MATURE = {
    "GastricPitCell": {
        "unit_count": 37.0,
        "residence": 3.1,
        "cl": "CL:0002179",
        "from": "GastricIsthmusStem",
        "requires": [],
        "count_evidence": f"{SRC_KARAM} part II: 'The pit region includes an average of 37 pit cells'",
        "count_tier": "experimental",
        "residence_evidence": f"{SRC_KARAM} part II: the migration to the gastric surface 'takes 3.1 days "
        "on the average'",
        "gate_note": "NO gate: every fetched loss-of-function sentence for a pit factor (FOXQ1, SPDEF, "
        "KLF4) leaves the cell present and empty or hyperplastic, so no knockout here removes the lineage "
        "that carries most of the stomach's turnover",
    },
    "GastricParietalCell": {
        "unit_count": 26.0,
        "residence": 54.0,
        "cl": "CL:0000162",
        "from": "GastricIsthmusStem",
        "requires": ["ESRRG"],
        "count_evidence": f"{SRC_KARAM} part IV: '3.6 cells are in the pit, 6.2 in the isthmus, 5.6 in the "
        "neck, and 10.6 in the base'",
        "count_tier": "experimental",
        "residence_evidence": f"{SRC_KARAM} part IV: 'The overall turnover time of parietal cells averages "
        "54 days'",
        "gate_note": "Adkins-Threats et al. 2024, Cell Stem Cell 31:886: deleting Esrrg in cells already "
        "committed to the parietal lineage leaves whole gastric units with no parietal cells at all",
    },
    "GastricChiefCell": {
        "unit_count": 67.5,
        "residence": CHIEF_RESIDENCE_DAYS,
        "cl": "CL:0000155",
        "from": "GastricNeckCell",
        "requires": ["MIST1"],
        "count_evidence": f"{SRC_KARAM} part III: 'Zymogenic cells, which average 67.5 per unit'",
        "count_tier": "experimental",
        "residence_evidence": "INFERRED: no number for the chief cell's residence was fetched in this "
        "lane. READINGS R5 registers in advance that the fit moves by less than 1% for any value between "
        "50 and 1000 days",
        "gate_note": "Ramsey et al. 2007, Development 134:211: the mature chief cells are uniformly "
        "structurally abnormal and the immature transitional cells accumulate - a maturation defect stated "
        "here as a deletion, registered partial in advance",
    },
    "GastricEnteroendocrineCell": {
        "unit_count": 1.0,
        "residence": ENDOCRINE_RESIDENCE_DAYS,
        "cl": "CL:1000222",
        "from": "GastricIsthmusStem",
        "requires": ["NEUROG3", "ASCL1"],
        "count_evidence": "INFERRED: the one lineage of the four with no measured count. 1.0 cell per unit "
        "is inside the <1% of the epithelium bound carried from the small-intestinal registration (Cheng & "
        "Leblond 1974 for stomach, small intestine and colon), and is fitted to nothing",
        "count_tier": "inferred",
        "residence_evidence": f"{SRC_HIBDON}: gastric endocrine cells have 'a half-life of approximately 2 "
        f"months'; {ENDOCRINE_RESIDENCE_DAYS:.3g} d is the mean residence that implies under the "
        "exponential loss law this module applies, and the conversion is the inference",
        "gate_note": "two gates on the one step: Lee et al. 2002 (the glandular stomach's A, D and G cells "
        "absent) and Hibdon et al. 2023 (the newborn Ascl1 null, with the other lineages normal)",
    },
}

# The progenitor compartments, top-down. `requires` is the registered gate of each step.
COMPARTMENTS = [
    {
        "id": "GastricIsthmusStem",
        "parent": None,
        "requires": ["RBPJ"],
        "residence": None,  # FITTED: the one fitted number in this module
        "cl": "CL:0002183",
        "evidence": "Kim & Shivdasani 2011, J Exp Med 208:677 (pan-Notch inhibition leaves no Ki67+ cells "
        "in most corpus and antral glands); REGISTERED CONTESTED: the same paper's genetic arm (Osr1-Cre; "
        "RBP-Jk Fl/Fl) leaves the units alive and merely 47% fewer by competition. The isthmus is also the "
        "amplifying compartment: Karam & Leblond part I give the granule-free cell the highest labelling "
        "index of the unit, 32.4%",
    },
    {
        "id": "GastricNeckCell",
        "parent": "GastricIsthmusStem",
        "requires": [],
        "residence": NECK_RESIDENCE_DAYS,
        "cl": "CL:0002181",
        "evidence": "Karam & Leblond 1993, Anat Rec 236:297 (part III): the mucous neck cell is the chief "
        "cell's progenitor - 'neck cells migrate inward in 1-2 weeks' and become pre-zymogenic and then "
        "zymogenic cells - and this module gives it exactly one outflow for that reason. No gate: no "
        "loss-of-function sentence removing the neck cell was fetched, and Ramsey et al. 2007 report the "
        "Mist1-null neck cells differentiating normally",
    },
]
NECK = "GastricNeckCell"
STEM = "GastricIsthmusStem"
# The cap the root carries, as the colonic transit-amplifying compartment did: a blocked outflow backs
# the compartment up to this multiple of its designed pool instead of growing without bound.
ROOT_CAP = 3.0


def plan(chief_residence_days: float | None = None) -> dict:
    """Solve the flows from the measured standing composition, then fit the one isthmus residence.

    The colon was solved downward from assumed output shares; the stomach is solved from measured
    standing counts, because its residences are measured too and differ by a factor of 60. For each
    mature lineage the steady state gives output = pool / residence, so the OUTPUT shares follow from
    the measured counts and residences as (n_i / r_i) normalised, and are not assumed. The published
    daily output is then distributed by those shares, every pool follows as output * residence, the
    neck compartment carries the chief route's whole throughput, and the one fitted number - the
    isthmus residence - is solved so that every pool together equals the published count.

    `chief_residence_days` overrides the one residence for which no number was fetched; it exists so
    READINGS R5's sensitivity statement can be checked rather than asserted.
    """
    residence = {m: d["residence"] for m, d in MATURE.items()}
    if chief_residence_days is not None:
        residence["GastricChiefCell"] = chief_residence_days

    rate = {m: MATURE[m]["unit_count"] / residence[m] for m in MATURE}
    total_rate = sum(rate.values())
    out = {m: DAILY_OUTPUT * rate[m] / total_rate for m in MATURE}

    p: dict[str, dict] = {}
    for m in MATURE:
        p[m] = {
            "id": m,
            "output": out[m],
            "inflow": out[m],
            "amplification": 1,
            "mature": True,
            "residence": residence[m],
            "pool": out[m] * residence[m],
            "loss_fraction": STEP_DAYS / residence[m],
        }

    from_neck = [m for m, d in MATURE.items() if d["from"] == NECK]
    from_stem = [m for m, d in MATURE.items() if d["from"] == STEM]

    neck = {
        "id": NECK,
        "output": sum(out[m] for m in from_neck),
        "amplification": 1,
        "residence": NECK_RESIDENCE_DAYS,
    }
    neck["inflow"] = neck["output"] / neck["amplification"]
    neck["pool"] = neck["output"] * NECK_RESIDENCE_DAYS
    p[NECK] = neck

    stem = {
        "id": STEM,
        "output": sum(out[m] for m in from_stem) + neck["inflow"],
        "amplification": 1,
    }
    stem["inflow"] = 0.0  # the root: self-renewal, not an inflow
    # The one fitted quantity: the isthmus residence, solved so that every pool together equals the
    # published stomach epithelial count. It is compared with the fetched upper bound of 2.6 days.
    fixed = neck["pool"] + sum(p[m]["pool"] for m in MATURE)
    fitted = (TOTAL_EPITHELIUM - fixed) / stem["output"]
    if fitted <= 0:
        raise ValueError(
            f"the measured composition already accounts for {fixed:.4g} of {TOTAL_EPITHELIUM:.3g}: "
            "no isthmus residence fits"
        )
    stem["residence"] = fitted
    stem["pool"] = stem["output"] * fitted
    p[STEM] = stem

    children = {STEM: [*from_stem, NECK], NECK: from_neck}
    for kids in children.values():
        total = sum(p[k]["inflow"] for k in kids) or 1.0
        for k in kids:
            p[k]["share"] = p[k]["inflow"] / total
    for cid, kids in children.items():
        c = p[cid]
        c["out_fraction"] = c["output"] / c["pool"]  # per day
        c["step_out"] = c["out_fraction"] * STEP_DAYS
        c["retained"] = 1.0 - sum(c["step_out"] * p[k]["share"] for k in kids)
        if cid == STEM:
            c["growth_fraction"] = 1.0 / c["retained"] - 1.0  # the root: self-renewal balances the outflow
        else:
            c["growth_fraction"] = (1.0 - c["step_out"] / c["amplification"]) / c["retained"] - 1.0

    designed = sum(p[k]["pool"] for k in [*MATURE, NECK, STEM])
    standing = {k: p[k]["pool"] / designed for k in [*MATURE, NECK, STEM]}
    p["_fit"] = {
        "isthmus_residence_days": fitted,
        "granule_free_turnover_bound_days": GRANULE_FREE_TURNOVER_BOUND_DAYS,
        "total_designed": designed,
        "implied_mean_residence_days": TOTAL_EPITHELIUM / DAILY_OUTPUT,
        "published_mean_residence_days": EPITHELIAL_LIFESPAN_DAYS,
        # Registered reading R3, computed from the measured inputs alone: the standing composition and
        # the output composition are not the same thing in the stomach, which is the whole arithmetic
        # difference from the colonic module, where one shared residence made them identical.
        "standing_share": standing,
        "output_share": {m: out[m] / DAILY_OUTPUT for m in MATURE},
        # The per-unit composition comes back BY CONSTRUCTION (it is what the module was built from) and
        # is not a check; the implied unit count is a consequence of the flows, claimed of nothing.
        "implied_units": p["GastricPitCell"]["pool"] / MATURE["GastricPitCell"]["unit_count"],
    }
    return p


def to_bio(module_name: str = "organism.human.stomach") -> str:
    p = plan()
    fit = p["_fit"]
    factors = sorted(
        {f for c in COMPARTMENTS for f in c["requires"]} | {f for d in MATURE.values() for f in d["requires"]}
    )
    lines = [
        "# Generated by genomeos.organism.stomach: the gastric gland's four mature lineages from one",
        "# isthmus compartment, gated by the factors of the closed sourced list in",
        "# scripts/stomach_registration.py, which was committed before any number below was fitted. Each of",
        "# its five gates rests on a loss-of-function phenotype fetched FOR THE STOMACH; ATOH1, the",
        "# secretory gate of both gut modules, is excluded here for want of one, and XBP1, HNF4A, SPDEF,",
        "# FOXQ1 and KLF4 are excluded because the stomach's best-sourced phenotypes leave the cell PRESENT",
        "# and spoiled, which a presence-or-absence gate cannot state. One gate, RBPJ, is REGISTERED",
        "# CONTESTED: the one paper's genetic and pharmacological arms differ in severity and this module",
        "# can satisfy one of them.",
        "# Measured and used as given: the stomach's count and daily output (Sender & Milo 2021); per",
        "# gastric unit in the mouse corpus, 37 pit, 12.6 neck, 67.5 zymogenic and 26.0 parietal cells, and",
        "# the pit (3.1 d) and parietal (54 d) residences (Karam & Leblond 1993, parts I-IV); the gastric",
        "# endocrine half-life of about 2 months (Hibdon et al. 2023). Inferred: the mouse-to-human carry-",
        "# over of the composition, the neck residence (the midpoint of a fetched 1-2 week range), the",
        "# endocrine residence (the half-life converted under exponential loss), the chief residence (no",
        "# number fetched), and the endocrine standing count. FITTED, and the only fitted number: the",
        f"# isthmus residence, {fit['isthmus_residence_days']:.3g} d, solved so the pools sum to the",
        f"# published {TOTAL_EPITHELIUM:.3g} cells - and held against the fetched UPPER BOUND of",
        f"# {GRANULE_FREE_TURNOVER_BOUND_DAYS:g} d for the isthmal granule-free cells, fetched before the",
        "# fit. The count assert below is therefore the closure of that fit and not a prediction; the",
        "# knockouts in stomach_mutants.bio are the prediction, because no mutant number is fitted.",
        "# The standing composition comes back in the measured per-unit proportions BY CONSTRUCTION and is",
        f"# not a check. The daily output is used as given and is NOT fitted: {DAILY_OUTPUT:.3g} cells per",
        f"# day, of which the pit lineage carries {fit['output_share']['GastricPitCell']:.3g} while standing",
        f"# at {fit['standing_share']['GastricPitCell']:.3g} of the epithelium.",
        f"module {module_name}",
        "import bio.std.development",
        "",
        "organism Stomach {",
        "  species: Homo sapiens",
        "  genome: GRCh38",
        f"  root: {STEM}",
        "  resolution: populations",
        f"  cell_type: {STEM}",
        f"  factors: {', '.join(factors)}",
        "  observe: count, fates",
        f"  assert: count at 1 yr in {TOTAL_EPITHELIUM * 0.85:.3g}..{TOTAL_EPITHELIUM * 1.15:.3g}",
    ]
    for m in MATURE:
        lines.append(f"  assert: type {m} at 1 yr in {p[m]['pool'] * 0.85:.3g}..{p[m]['pool'] * 1.15:.3g}")
    lines += [
        f'  evidence: curated "{SRC_MILO}; {SRC_KARAM}; {SRC_HIBDON}"; confidence: 0.4',
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
        lines.append(
            f"cell_type {c['id']} {{ parent: {parent}; ontology: {c['cl']}; "
            + f'evidence: experimental "{c["evidence"]}"; confidence: 0.8 }}'
        )
    for m, d in MATURE.items():
        lines.append(
            f"cell_type {m} {{ parent: PostMitotic; ontology: {d['cl']}; "
            + 'evidence: curated "Cell Ontology, looked up in EBI OLS4"; confidence: 0.9 }'
        )
    stem = p[STEM]
    lines += [
        "",
        f"# the isthmus: {stem['output']:.3g} cells leave it per day, which is the whole daily output, and",
        f"# the ONE FITTED NUMBER is its residence, {stem['residence']:.4g} d, solved so that every pool",
        f"# together equals the published {TOTAL_EPITHELIUM:.3g} cells. Karam & Leblond part I put the",
        f"# granule-free cell's turnover at less than {GRANULE_FREE_TURNOVER_BOUND_DAYS:g} d, fetched before",
        "# the fit, so the knob had a bound to miss. The pool it implies is",
        f"# {stem['pool']:.3g} cells, and at the measured 37 pit cells per gastric unit the module would",
        f"# imply {fit['implied_units']:.3g} units - a consequence of the flows, checked against nothing.",
        "# The root carries the cap the colonic transit-amplifying compartment carried, because in the",
        "# stomach the root IS the amplifying compartment: a blocked outflow backs it up to",
        f"# {ROOT_CAP:g} times its designed pool instead of growing without bound.",
        *_fill_lines(
            STEM, stem["pool"], stem["retained"], COMPARTMENTS[0]["requires"], stem["growth_fraction"]
        ),
        "decision isthmus_renew { action: divide; "
        + f"when: cell_type = {STEM}{_req(COMPARTMENTS[0]['requires'])}, "
        + f"count = <{stem['pool'] * ROOT_CAP / (1.0 + stem['growth_fraction']):.3g}; "
        + f"fraction: {stem['growth_fraction']:.6g}; "
        + f"after: {STEP_DAYS:g} d; "
        + 'evidence: inferred "self-renewal balances the outflow: the isthmus residence is the one fitted '
        + f'number, {stem["residence"]:.4g} d"; confidence: 0.3 }}',
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
        "# the neck compartment, the chief cell's progenitor, with exactly ONE outflow: with the chief",
        "# route gated off nothing leaves it (registered reading R4), which is not the colonic secretory",
        "# progenitor's case, where a second outflow carried the whole inflow",
    ]
    for c in COMPARTMENTS:
        if c["parent"] is None:
            continue
        pc = p[c["id"]]
        lines.extend(_fill_lines(c["id"], pc["pool"], pc["retained"], c["requires"], pc["growth_fraction"]))
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
            + f'evidence: {d["count_tier"]} "standing count {d["unit_count"]:g} per gastric unit and '
            + f"residence {p[m]['residence']:.4g} d give this share {p[m]['share']:.4g} of the inflow: "
            + f'{d["count_evidence"]}"; confidence: 0.3 }}'
        )
    lines += [
        "",
        "# every mature lineage is filled to its adult pool during the Fill stage, which neither gut module",
        "# needed: a lineage that lives 190 days cannot reach its pool inside 90 days by inflow alone",
    ]
    for m, d in MATURE.items():
        lines.extend(_fill_lines(m, p[m]["pool"], 1.0, d["requires"], 0.0))
    lines.append("")
    for m, d in MATURE.items():
        lines.append(
            f"decision {m}_loss {{ action: die; when: cell_type = {m}; "
            + f"fraction: {p[m]['loss_fraction']:.6g}; after: {STEP_DAYS:g} d; "
            + f'evidence: inferred "residence {p[m]["residence"]:.4g} d: {d["residence_evidence"]}"; '
            + "confidence: 0.2 }"
        )
    return "\n".join(lines) + "\n"


# The knockouts, exactly as scripts/stomach_registration.py registered them before this module's
# numbers were fitted. `expect` carries the published stomach phenotype, including the part a module
# made of flows cannot express and, for the contested gate, the experiment the module disagrees with.
MUTANTS = [
    (
        "rbpj",
        "RBPJ",
        "Kim & Shivdasani 2011, J Exp Med 208:677",
        "pan-Notch inhibition leaves no Ki67+ cells in most corpus and antral glands, so no lineage is "
        "produced (REGISTERED CONTESTED in advance: the same paper's genetic arm, Osr1-Cre; RBP-Jk Fl/Fl, "
        "instead leaves the gastric units alive and merely 47% fewer by competition, and this module has "
        "one gate and cannot satisfy both)",
        ["GastricPitCell", "GastricParietalCell", "GastricChiefCell", "GastricEnteroendocrineCell"],
        [],
    ),
    (
        "esrrg",
        "ESRRG",
        "Adkins-Threats et al. 2024, Cell Stem Cell 31:886",
        "whole swaths of gastric units completely lack parietal cells, from birth onward; the other "
        "lineages are untouched because Atp4bCre deletes Esrrg only in cells already committed to the "
        "parietal lineage (the paper's steady-state increase in proliferation is not expressible here)",
        ["GastricParietalCell"],
        ["GastricPitCell", "GastricChiefCell", "GastricEnteroendocrineCell"],
    ),
    (
        "mist1",
        "MIST1",
        "Ramsey et al. 2007, Development 134:211",
        "the mature chief cells of the gland base are uniformly structurally abnormal and the immature "
        "transitional cells upstream accumulate markedly, while the neck cells differentiate normally; a "
        "presence-or-absence gate can only delete the mature cell, so the phenotype is overstated: "
        "registered partial in advance",
        ["GastricChiefCell"],
        ["GastricPitCell", "GastricParietalCell", "GastricEnteroendocrineCell"],
    ),
    (
        "neurog3",
        "NEUROG3",
        "Lee et al. 2002, Genes Dev 16:1488",
        "the glandular stomach's glucagon-, somatostatin- and gastrin-secreting cells are absent and the "
        "enterochromaffin cells dramatically decreased; registered partial in advance and with NO kept "
        "lineage, because the same mice show intestinal metaplasia of the gastric epithelium and the paper "
        "therefore does not report the other lineages as normal",
        ["GastricEnteroendocrineCell"],
        [],
    ),
    (
        "ascl1",
        "ASCL1",
        "Hibdon et al. 2023, Am J Physiol Gastrointest Liver Physiol 325:G458",
        "newborn Ascl1 nulls lose the neurogenin-3-dependent endocrine cells of the stomach with normal "
        "mucous, chief and parietal cells; registered partial in advance because the enterochromaffin-like "
        "cells are normal too and this module has no endocrine subtypes, and because the readout is a "
        "newborn one",
        ["GastricEnteroendocrineCell"],
        ["GastricPitCell", "GastricParietalCell", "GastricChiefCell"],
    ),
]
KEPT_FLOOR = 1.0e6  # the floor the registration set for a population a knockout does not touch


def mutants_bio(module_name: str = "organism.human.stomach_mutants") -> str:
    lines = [
        "# The five knockouts of the closed sourced list in scripts/stomach_registration.py, against the",
        "# published STOMACH phenotypes: each removes one factor and states what must disappear and what",
        "# must remain. Only ESRRG is registered FULL. Three are registered PARTIAL in advance, because",
        "# their published phenotype is a maturation defect or a partial loss and not a deletion, and ONE",
        "# is registered CONTESTED, because the one paper's genetic and pharmacological arms differ in",
        "# severity: rbpj can agree with one of them only and is not counted as reproducing the stomach.",
        "# neurog3 carries NO kept lineage, because its paper reports intestinal metaplasia and so does not",
        "# report the other lineages as normal. Run:",
        "# genomeos grow data/organisms/human/stomach_mutants.bio --experiments",
        f"module {module_name}",
        "import stomach.bio",
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
