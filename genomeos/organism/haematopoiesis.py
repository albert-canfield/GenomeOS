# SPDX-License-Identifier: AGPL-3.0-or-later
"""Haematopoiesis as a mechanism module: from the stem cell to the blood lineages.

The first human module in which lineage choices are decided by named factors rather
than by shares alone. Each compartment requires the transcription factors that the
genetics showed to be necessary (knockouts in mice, cited per compartment), so a
knockout in an `experiment` block removes the lineages that depend on them. The
quantities are fitted to the measured daily outputs and lifespans of Sender & Milo
2021: a compartment's outflow is its output, its residence time sets its pool, its
amplification (divisions per transit) sets what it needs from upstream, and the sum
of everything upstream sets what the stem-cell pool has to deliver per day. The
implied stem-cell division interval is reported and compared with the published
estimates; the amplifications are inferred, and say so.
"""

from __future__ import annotations

import math

# Mature outputs per day and lifespans (Sender & Milo 2021 table, grouped as in human.py).
MATURE = {
    "Erythrocyte": {"output": 2.14e11, "lifespan": 119.0, "cl": "CL:0000232"},
    "Neutrophil": {"output": 6.0e10, "lifespan": 6.6, "cl": "CL:0000775"},
    "Monocyte": {"output": 1.52e9, "lifespan": 3.5, "cl": "CL:0000576"},
    "BCell": {"output": 5.4e9, "lifespan": 63.0, "cl": "CL:0000236"},
    "TCell": {"output": 2.2e9, "lifespan": 320.0, "cl": "CL:0000084"},
}

# Compartments: parent, the factors each requires (genetic necessity, cited), residence in days and
# amplification (2^divisions per transit; inferred), and the mature type a terminal compartment feeds.
COMPARTMENTS = [
    {
        "id": "HSC",
        "parent": None,
        "requires": ["TAL1", "RUNX1", "GATA2"],
        "residence": None,
        "amplification": 1,
        "cl": "CL:0000037",
        "evidence": "Shivdasani et al. 1995 (TAL1); Okuda et al. 1996 (RUNX1); Tsai et al. 1994 (GATA2)",
    },
    {
        "id": "MPP",
        "parent": "HSC",
        "requires": ["GATA2", "MYB"],
        "residence": 7.0,
        "amplification": 2**14,
        "cl": "CL:0000837",
        "evidence": "Busch et al. 2015 (a self-renewing progenitor tier, not HSCs, supplies the steady state)"
        "; Mucenski et al. 1991 (MYB)",
    },
    {
        "id": "CMP",
        "parent": "MPP",
        "requires": ["GATA2"],
        "residence": 3.0,
        "amplification": 2**6,
        "cl": "CL:0000049",
        "evidence": "Akashi et al. 2000, Nature 404:193 (CMP)",
    },
    {
        "id": "CLP",
        "parent": "MPP",
        "requires": ["IKZF1", "TCF3"],
        "residence": 3.0,
        "amplification": 2**6,
        "cl": "CL:0000051",
        "evidence": "Kondo et al. 1997 (CLP); Georgopoulos et al. 1994 (Ikaros: no lymphoid lineages)",
    },
    {
        "id": "MEP",
        "parent": "CMP",
        "requires": ["GATA1"],
        "residence": 3.0,
        "amplification": 2**6,
        "cl": "CL:0000050",
        "evidence": "Pevny et al. 1991; Fujiwara et al. 1996 (GATA1: no erythroid or megakaryocytic lineage)",
    },
    {
        "id": "GMP",
        "parent": "CMP",
        "requires": ["SPI1", "CEBPA"],
        "residence": 3.0,
        "amplification": 2**6,
        "cl": "CL:0000557",
        "evidence": "Scott et al. 1994 (PU.1: no myeloid or B lineages); Zhang et al. 1997 (C/EBPa)",
    },
    {
        "id": "Erythroblast",
        "parent": "MEP",
        "requires": ["GATA1", "KLF1"],
        "residence": 7.0,
        "amplification": 2**5,
        "cl": "CL:0000765",
        "feeds": "Erythrocyte",
        "evidence": "Nuez et al. 1995; Perkins et al. 1995 (KLF1/EKLF: no mature erythrocytes); "
        "~5 divisions from proerythroblast to reticulocyte",
    },
    {
        "id": "Myeloblast",
        "parent": "GMP",
        "requires": ["CEBPA", "GFI1"],
        "residence": 7.0,
        "amplification": 2**4,
        "cl": "CL:0000835",
        "feeds": "Neutrophil",
        "evidence": "Hock et al. 2003 (GFI1: neutropenia); ~4 divisions myeloblast to metamyelocyte",
    },
    {
        "id": "Monoblast",
        "parent": "GMP",
        "requires": ["SPI1", "IRF8"],
        "residence": 3.0,
        "amplification": 2**3,
        "cl": "CL:0000040",
        "feeds": "Monocyte",
        "evidence": "Holtschke et al. 1996 (IRF8/ICSBP: monocyte and macrophage deficiency)",
    },
    {
        "id": "ProB",
        "parent": "CLP",
        "requires": ["SPI1", "TCF3", "EBF1", "PAX5"],
        "residence": 5.0,
        "amplification": 2**5,
        "cl": "CL:0000826",
        "feeds": "BCell",
        "evidence": "Bain et al. 1994 (E2A); Lin & Grosschedl 1995 (EBF1); "
        "Urbanek et al. 1994 (PAX5: block at pro-B)",
    },
    {
        "id": "Thymocyte",
        "parent": "CLP",
        "requires": ["NOTCH1", "TCF7", "GATA3"],
        "residence": 21.0,
        "amplification": 2**7,
        "cl": "CL:0000893",
        "feeds": "TCell",
        "evidence": "Radtke et al. 1999 (Notch1: no T cells); Verbeek et al. 1995 (TCF-1); "
        "Ting et al. 1996 (GATA-3)",
    },
]
HSC_POOL = 1.0e5  # Lee-Six et al. 2018, Nature 561:473: 50,000 to 200,000 human HSCs
STEP_DAYS = 0.25  # flows and growth are applied every quarter day: the steady state of a compartment with a
# large amplification is a small difference between large flows, and a daily step leaves it 10-20% off
HSC_DIVISION_DAYS_PUBLISHED = 280.0  # Catlin et al. 2011, Blood 117:4460: about once per 40 weeks


def plan() -> dict:
    """Solve the flows: required daily inflow per compartment from the mature outputs upward."""
    by_id = {c["id"]: dict(c) for c in COMPARTMENTS}
    children = {c["id"]: [x["id"] for x in COMPARTMENTS if x["parent"] == c["id"]] for c in COMPARTMENTS}
    order = [c["id"] for c in reversed(COMPARTMENTS)]  # leaves before parents (table is top-down)
    for cid in order:
        c = by_id[cid]
        if c.get("feeds"):
            c["output"] = MATURE[c["feeds"]]["output"]
        else:
            c["output"] = sum(by_id[k]["inflow"] for k in children[cid])
        c["inflow"] = c["output"] / c["amplification"]
        if c["residence"]:
            c["pool"] = c["output"] * c["residence"]
            c["out_fraction"] = 1.0 / c["residence"]
    hsc = by_id["HSC"]
    hsc["pool"] = HSC_POOL
    hsc["out_fraction"] = hsc["output"] / HSC_POOL  # per day
    hsc["implied_division_days"] = 1.0 / hsc["out_fraction"]
    for _cid, kids in children.items():
        total = sum(by_id[k]["inflow"] for k in kids) or 1.0
        for k in kids:
            by_id[k]["share"] = by_id[k]["inflow"] / total
    # The runtime applies each day's growth and outflows one after another, so they multiply. The
    # growth step is therefore solved for the discrete balance: (1 + g) * prod(1 - f_i) = 1 - o / A,
    # where the f_i are the day's outflow fractions and o / A the share the inflow must supply.
    for cid, c in by_id.items():
        c["step_out"] = c["out_fraction"] * STEP_DAYS  # outflow fraction per step
        outs = [c["step_out"] * by_id[k]["share"] for k in children[cid]] or [c["step_out"]]
        c["retained"] = math.prod(1.0 - f for f in outs)
        if cid == "HSC":
            c["growth_fraction"] = 1.0 / c["retained"] - 1.0  # self-renewal exactly balances the outflow
        else:
            c["growth_fraction"] = (1.0 - c["step_out"] / c["amplification"]) / c["retained"] - 1.0
    return by_id


def _req(factors: list[str]) -> str:
    return "".join(f", {f} = present" for f in factors)


def _fill_lines(cell_type: str, pool: float, retained: float, requires: list[str]) -> list[str]:
    """Fill a compartment to its steady pool during the Fill stage: doubling to half, then 10% steps to
    nine tenths, then 1% steps, so the pool lands within 1% and the slow steady-state dynamics start
    from the right place (an adult module describes an organism whose compartments already exist)."""
    ev = 'evidence: inferred "an adult compartment starts at its steady pool"; confidence: 0.3 }'
    out = []
    for tag, step, upto in (("2x", 1.0, 0.5), ("10pc", 0.1, 0.9), ("1pc", 0.01, 1.0)):
        out.append(  # the step is on top of the daily outflow, so the pool really grows by it
            f"decision {cell_type}_fill_{tag} {{ action: divide; "
            + f"when: cell_type = {cell_type}{_req(requires)}, stage = Fill, count = <{pool * upto:.4g}; "
            + f"fraction: {(1.0 + step) ** STEP_DAYS / retained - 1.0:.6g}; "
            + f"after: {STEP_DAYS:g} d; {ev}"
        )
    return out


def _pool_range(mid: str, low: float = 0.85, high: float = 1.15) -> str:
    steady = MATURE[mid]["output"] * MATURE[mid]["lifespan"]
    return f"{steady * low:.3g}..{steady * high:.3g}"


def to_bio(module_name: str = "organism.human.haematopoiesis") -> str:
    p = plan()
    src = "Sender & Milo 2021, Nat Med 27:45"
    factors = sorted({f for c in COMPARTMENTS for f in c["requires"]})
    hsc = p["HSC"]
    lines = [
        "# Generated by genomeos.organism.haematopoiesis: the blood lineages from the stem cell,",
        "# decided by the factors the genetics showed to be necessary; flows fitted to the measured",
        "# daily outputs and lifespans. Amplifications are inferred; requirements are cited per compartment.",
        f"module {module_name}",
        "import bio.std.development",
        "",
        "organism Haematopoiesis {",
        "  species: Homo sapiens",
        "  genome: GRCh38",
        "  root: HSC",
        "  resolution: populations",
        "  cell_type: HSC",
        f"  factors: {', '.join(factors)}",
        "  observe: count, fates",
        f"  assert: type Erythrocyte at 3 yr in {_pool_range('Erythrocyte')}",
        f"  assert: type Neutrophil at 3 yr in {_pool_range('Neutrophil')}",
        f"  assert: type TCell at 3 yr in {_pool_range('TCell', 0.8)}",
        f'  evidence: curated "{src}; Lee-Six et al. 2018; Catlin et al. 2011"; confidence: 0.5',
        "}",
        "",
        'stage Fill { from: 0 d; to: 90 d; evidence: inferred "every compartment fills to its adult size"; '
        + "confidence: 0.3 }",
        'stage Steady { from: 90 d; evidence: curated "adult steady state"; confidence: 0.7 }',
        "",
    ]
    for c in COMPARTMENTS:
        parent = "Blastomere" if c["parent"] is None else c["parent"]
        lines.append(
            f"cell_type {c['id']} {{ parent: {parent}; ontology: {c['cl']}; "
            + f'evidence: experimental "{c["evidence"]}"; confidence: 0.8 }}'
        )
    for mid, m in MATURE.items():
        lines.append(
            f"cell_type {mid} {{ parent: PostMitotic; ontology: {m['cl']}; "
            + 'evidence: curated "Cell Ontology"; confidence: 0.9 }'
        )
    lines += [
        "",
        f"# stem-cell pool: {HSC_POOL:.3g} cells (Lee-Six 2018); output {hsc['output']:.3g}/day, so each HSC",
        f"# leaves the pool every {hsc['implied_division_days']:.0f} days and one self-renewal division",
        "# balances the outflow",
        f"# (Catlin 2011: about one division per {HSC_DIVISION_DAYS_PUBLISHED:.0f} days)",
        *_fill_lines("HSC", HSC_POOL, hsc["retained"], COMPARTMENTS[0]["requires"]),
        f"decision hsc_renew {{ action: divide; when: cell_type = HSC{_req(COMPARTMENTS[0]['requires'])}, "
        + f"count = >={HSC_POOL * 0.99:.3g}; "
        + f"fraction: {hsc['growth_fraction']:.6g}; after: {STEP_DAYS:g} d; "
        + 'evidence: inferred "self-renewal balances the outflow: one division per '
        + f'{hsc["implied_division_days"]:.0f} days per HSC"; confidence: 0.3 }}',
        "",
    ]
    for c in COMPARTMENTS:
        if c["parent"] is None:
            continue
        pc = p[c["id"]]
        parent_requires = next(x["requires"] for x in COMPARTMENTS if x["id"] == c["parent"])
        req = _req(
            sorted(set(parent_requires) | set(c["requires"]))
        )  # the source must exist, the target be allowed
        out_frac = p[c["parent"]]["step_out"] * pc["share"]
        lines.append(
            f"decision {c['parent']}_to_{c['id']} {{ action: differentiate; "
            + f"when: cell_type = {c['parent']}{req}; to: {c['id']}; fraction: {out_frac:.6g}; "
            + f"after: {STEP_DAYS:g} d; "
            + f'evidence: experimental "{c["evidence"]}"; confidence: 0.7 }}'
        )
    lines += [
        "",
        "# amplification inside each compartment (inferred divisions per transit), capped at three times the",
        "# steady pool so that a blocked outflow backs up instead of growing without bound",
    ]
    for c in COMPARTMENTS:
        if c["parent"] is None:
            continue
        pc = p[c["id"]]
        divisions = int(math.log2(c["amplification"]))
        lines.extend(_fill_lines(c["id"], pc["pool"], pc["retained"], c["requires"]))
        lines.append(
            f"decision {c['id']}_amplify {{ action: divide; "
            + f"when: cell_type = {c['id']}{_req(c['requires'])}, count = <{pc['pool'] * 3:.3g}; "
            + f"fraction: {pc['growth_fraction']:.6g}; after: {STEP_DAYS:g} d; "
            + f'evidence: inferred "{divisions} divisions per transit of {c["residence"]:g} days; '
            + f'pool {pc["pool"]:.3g}"; confidence: 0.3 }}'
        )
    lines.append("")
    for c in COMPARTMENTS:
        if not c.get("feeds"):
            continue
        pc = p[c["id"]]
        out = MATURE[c["feeds"]]["output"]
        lines.append(
            f"decision {c['id']}_to_{c['feeds']} {{ action: differentiate; when: cell_type = {c['id']}; "
            + f"to: {c['feeds']}; fraction: {pc['step_out']:.6g}; after: {STEP_DAYS:g} d; "
            + f'evidence: experimental "{src}: {out:.3g} {c["feeds"]} per day"; confidence: 0.7 }}'
        )
    lines.append("")
    for mid, m in MATURE.items():
        loss = STEP_DAYS / m["lifespan"]
        lines.append(
            f"decision {mid}_loss {{ action: die; when: cell_type = {mid}; fraction: {loss:.6g}; "
            + f"after: {STEP_DAYS:g} d; "
            + f'evidence: experimental "{src}: lifespan {m["lifespan"]:g} days"; confidence: 0.7 }}'
        )
    return "\n".join(lines) + "\n"


MUTANTS = [
    (
        "gata1",
        "GATA1",
        "Pevny et al. 1991, Nature 349:257; Fujiwara et al. 1996",
        "no erythrocytes (and no megakaryocytes); other lineages intact",
        ["Erythrocyte"],
        ["Neutrophil", "BCell", "TCell"],
    ),
    (
        "klf1",
        "KLF1",
        "Nuez et al. 1995, Nature 375:316; Perkins et al. 1995",
        "no mature erythrocytes (fatal anaemia); myeloid and lymphoid intact",
        ["Erythrocyte"],
        ["Neutrophil", "TCell"],
    ),
    (
        "spi1",
        "SPI1",
        "Scott et al. 1994, Science 265:1573",
        "no neutrophils, monocytes or B cells; erythrocytes and T cells present",
        ["Neutrophil", "Monocyte", "BCell"],
        ["Erythrocyte", "TCell"],
    ),
    (
        "cebpa",
        "CEBPA",
        "Zhang et al. 1997, PNAS 94:569",
        "no neutrophils; monocytes are not lost this way (they come from GMP through IRF8)",
        ["Neutrophil"],
        ["Erythrocyte", "BCell", "TCell"],
    ),
    (
        "irf8",
        "IRF8",
        "Holtschke et al. 1996, Cell 87:307",
        "monocyte and macrophage deficiency",
        ["Monocyte"],
        ["Neutrophil", "Erythrocyte"],
    ),
    (
        "pax5",
        "PAX5",
        "Urbanek et al. 1994, Cell 79:901; Nutt et al. 1999",
        "no B cells: a block at the pro-B stage",
        ["BCell"],
        ["TCell", "Erythrocyte", "Neutrophil"],
    ),
    (
        "notch1",
        "NOTCH1",
        "Radtke et al. 1999, Immunity 10:547",
        "no T cells",
        ["TCell"],
        ["BCell", "Erythrocyte", "Neutrophil"],
    ),
    (
        "ikzf1",
        "IKZF1",
        "Georgopoulos et al. 1994, Cell 79:143",
        "no B or T cells (no lymphoid lineages)",
        ["BCell", "TCell"],
        ["Erythrocyte", "Neutrophil"],
    ),
    (
        "tal1",
        "TAL1",
        "Shivdasani et al. 1995, Nature 373:432",
        "no blood at all: no haematopoietic stem cells",
        ["Erythrocyte", "Neutrophil", "Monocyte", "BCell", "TCell"],
        [],
    ),
]


def mutants_bio() -> str:
    lines = [
        "# Haematopoietic knockouts against the published (mouse) phenotypes: each removes one factor from",
        "# the module and states what should disappear and what should remain. Run:",
        "# genomeos grow data/organisms/human/haematopoiesis_mutants.bio --experiments",
        "module organism.human.haematopoiesis_mutants",
        "import haematopoiesis.bio",
        "",
    ]
    for name, factor, source, expect, lost, kept in MUTANTS:
        asserts = "".join(f" assert: type {t} at 3 yr = 0;" for t in lost)
        asserts += "".join(f" assert: type {t} at 3 yr >= 1e9;" for t in kept)
        lines.append(
            f'experiment {name} {{ knockout: {factor}; until: 3 yr; expect: "{expect}";{asserts} '
            + f'evidence: experimental "{source}"; confidence: 0.8 }}'
        )
    return "\n".join(lines) + "\n"
