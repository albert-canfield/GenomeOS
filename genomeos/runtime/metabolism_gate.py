# SPDX-License-Identifier: Apache-2.0
# Part of the BioLang engine (language, IR, VM, standard library); see LICENSING.md.
"""Pre-registration for stage 3's gate, written before the instrument and before any run.

§6 states the stage and the gate (docs/BIOLANG-v0.4-ECONOMY.md): *"A reduced metabolic module
(glycolysis, TCA, oxidative phosphorylation) in `bio.std.core_metabolism`, producing the ATP the pools
spend, with stoichiometries from the textbook and yields with their spread (about 2.5 ATP per NADH and
1.5 per FADH2, Hinkle 2005, BBA 1706:1; about 30 to 32 ATP per glucose aerobically against 2
anaerobically). Mitochondria are compartments with `genome: chrM` and `copies: N`. Each copy may carry
its own genome variant, which makes heteroplasmy expressible: `heteroplasmy: m.3243A>G = 0.7`, with the
biochemical threshold effect (Rossignol et al. 2003, Biochem J 370:751) as the check."* Its gate: *"(i)
The ATP budget of a modelled cell type lands within the published order of magnitude, and its split
across consumers matches the hierarchy of Buttgereit & Brand 1995 (Biochem J 312:163: protein synthesis
and the Na⁺/K⁺-ATPase first). (ii) Lowering oxygen shifts production to glycolysis and raises glucose
consumption (the Pasteur effect); removing glucose collapses a glycolysis-only cell. (iii) The red blood
cell (no mitochondria) runs on glycolysis alone and keeps its ATP and 2,3-BPG within the ranges of
Mulquiney & Kuchel 1999 (Biochem J 342:581)."*

**Whether stage 3 may be built at all (checked first, 2026-10-10).** The spec was a draft "for Albert's
review" and says nothing past stage 1 is built until he has steered the constructs. Of §10's ten
decisions the one stage 3 touches is 4, "whose cost accounting", and §10's own table records it as
settled BY this stage rather than gating it: *"4. whose cost accounting | measurable | the ATP budget
gate (Buttgereit & Brand's hierarchy, Lynch & Marinov's totals): the accounting that lands in the
published order of magnitude wins"*, and §10's closing note: *"4 waits on stage 3's ATP budget, which
is not implemented"*. Decision 3 (membranes) is preference and both shapes passed stage 1; 6 is stage
4's and 7 the control layer's. ROADMAP-HISTORY item 5, "Albert's language decisions": *"No decision
blocks code any more"*; and ROADMAP area A lists stage 3 as next item 2 with no steer pending (the
control layer, item 5, is the one that "waits for Albert's steer"). So nothing blocks it, and decision
4 stays open unless clause (i)'s magnitude can be asked — which, below, it cannot.

**What the three sources actually say, read before anything was built** (PubMed abstracts, fetched
2026-10-10; the bodies of the 1995 and 1999 papers are scans this project does not hold):

- Buttgereit & Brand 1995 (PMID 7492307): *"There was a clear hierarchy of the responses of different
  energy-consuming reactions to changes in energy supply: pathways of macromolecule biosynthesis
  (protein synthesis and RNA/DNA synthesis) were most sensitive to energy supply, followed by sodium
  cycling and then calcium cycling across the plasma membrane. Mitochondrial proton leak was the least
  sensitive to energy supply."* The hierarchy is an ORDER OF SENSITIVITY TO SUPPLY, not a split of the
  budget. §6's parenthesis, "protein synthesis and the Na⁺/K⁺-ATPase first", reads it as a split; the
  abstract states no shares and puts sodium cycling AFTER macromolecule synthesis.
- Hinkle 2005 (PMID 15620362): *"Values of about 2.5 with NADH-linked substrates and 1.5 with succinate
  are consistent with most reports"*, and *"P/O ratios of 2.3 with NADH and 1.4 with succinate, values
  that are also possible"*. That is the spread the yields carry here.
- Rossignol et al. 2003 (PMID 12467494): *"it was possible to inhibit considerably the activity of a
  respiratory chain complex, up to a critical value, without affecting the rate of mitochondrial
  respiration or ATP synthesis. This phenomenon was called the 'biochemical threshold effect'."*
- Mulquiney, Bubb & Kuchel 1999 (Biochem J 342:567, PMID 10477268), the first of the three papers whose
  second is §6's citation: *"the percentage of glycolytic carbon flux that passes through the 2,3-BPG
  shunt in the normal in vivo steady state was estimated to be 19%."* None of the three abstracts states
  the ATP or 2,3-BPG concentration ranges §6 names.

**The design, so the run cannot choose it.** One energy currency, refilled per hour by the processes a
pool names in `regenerates: from` — Glycolysis, TCA and OXPHOS, with every stoichiometric number read
from the program (`bio.std.core_metabolism`), none held by the engine. Supply follows demand: pyruvate is
oxidised as far as oxygen, mitochondrial capacity and glucose allow, and the rest ends as lactate. That
routing is the reduced module's ONE structural claim, and it is the Pasteur-competent cell by
construction; the Warburg effect (lactate with oxygen present) is therefore not expressible in stage 3,
and clause 3 passing is a property of the routing, not evidence about cells. Every yield-dependent
number is computed at both ends of the spread: the HIGH arm (malate–aspartate shuttle, P/O 2.5 and 1.5)
and the LOW arm (glycerol 3-phosphate shuttle, P/O 2.3 and 1.4), and a clause holds only if it holds in
both. Heteroplasmy removes the OXPHOS capacity of the copies that carry the declared variant; the engine
holds no table of which mtDNA variants are pathogenic, so a variant is declared there because the
program models it as disabling.

**The arithmetic that can be written down now.** Per glucose: glycolysis nets 2 ATP, 2 cytosolic NADH
and 2 pyruvate; each pyruvate oxidised gives 4 NADH, 1 FADH2 and 1 GTP; each electron pair reduces half
an O2. So a glucose oxidised in full yields 2 + 2 + 8·P/O(NADH) + 2·P/O(FADH2) + 2·P/O(shuttle) = **32.0**
in the high arm and **28.0** in the low arm (**30.0** with P/O 2.5 and 1.5 through the glycerol
3-phosphate shuttle, which is the low end of §6's "30 to 32"), and needs 12 electron pairs and exactly
**6 O2** in both arms: the combustion stoichiometry, which is how the module's constants are checked.
A glycolysis-only cell gets 2 − 2f per glucose, where f is the fraction of triose flux through the
2,3-BPG shunt (which bypasses phosphoglycerate kinase): **1.62** at the 19% above. With a demand D met
in full, glucose uptake is D/Y and O2 uptake 6·D/Y; with no oxygen it is D/(2 − 2f), so anoxia multiplies
glucose uptake by Y/(2 − 2f): **16.0** (high) and **14.0** (low) when uptake does not limit. The
biochemical threshold is where the copies left unmutated can no longer carry the electron pairs the
demand needs: h* = 1 − 12·D/(Y·C) for a capacity C, so it is set by the program's declared reserve —
and with no reserve there is no threshold at all, which is the falsifier.

**What clause 7 can and cannot show.** Stage 2's finding returns here: `proportional` and `competitive`
hand every demander of a pool one scalar (§9.2), so under either of them every ATP consumer falls by the
same factor and NO hierarchy can appear, whatever the supply. `priority` is the only policy that can
carry an order, and it carries it because the program writes it. So the verdict this registration
commits to in advance is: **the hierarchy is reproduced only when written in as a priority order; no
policy in the engine predicts it.** And priority carries it as cut-offs — a consumer is served in full
until the next one is cut — where the abstract describes graded sensitivities; the order is the part
that is comparable.
"""

from __future__ import annotations

from typing import Any

#: the two demonstration programs, committed with the language surface (they use `heteroplasmy:`, which
#: does not parse until then). Every number below that they declare is registered HERE, and the gate
#: script refuses to run when a program disagrees with this file.
PROGRAM = "data/demo/stage3_energy.bio"
RBC_PROGRAM = "data/demo/stage3_erythrocyte.bio"

#: the two ends of the yields' spread (Hinkle 2005): P/O for NADH, for FADH2, and for cytosolic NADH,
#: which reaches the chain as NADH through the malate-aspartate shuttle and as FADH2 through the
#: glycerol 3-phosphate shuttle.
ARMS: dict[str, dict[str, float]] = {
    "high": {"po_nadh": 2.5, "po_fadh2": 1.5, "po_cytosolic": 2.5},
    "low": {"po_nadh": 2.3, "po_fadh2": 1.4, "po_cytosolic": 1.4},
}

#: registered yields, ATP per glucose, computed by hand from the textbook stoichiometry above
PREDICTED_YIELD_AEROBIC = {"high": 32.0, "low": 28.0}
PREDICTED_YIELD_SPEC_LOW_END = 30.0  # P/O 2.5 and 1.5 through the glycerol 3-phosphate shuttle
PREDICTED_YIELD_ANAEROBIC = 2.0
PREDICTED_O2_PER_GLUCOSE = 6.0
PREDICTED_PAIRS_PER_GLUCOSE = 12.0

# ---- the nucleated cell (data/demo/stage3_energy.bio): every number DECLARED, none measured ----------

#: four ATP consumers at an equal 2.5e9 ATP per hour each. Equal on purpose: no share is held for any
#: cell type (see not_askable), so the program claims no split, and the hierarchy clause tests an order
#: only. Per-event costs are the textbook ones: 4 ATP per peptide bond and 2 per nucleotide (§5.2,
#: polymerisation only), 1 ATP per Na+/K+-ATPase cycle, 1 per Ca2+ pumped by the plasma-membrane ATPase.
CONSUMERS: dict[str, dict[str, float]] = {
    "protein_synthesis": {"rate_per_h": 6.25e8, "atp_per_event": 4.0},
    "rna_dna_synthesis": {"rate_per_h": 1.25e9, "atp_per_event": 2.0},
    "na_k_pumping": {"rate_per_h": 2.5e9, "atp_per_event": 1.0},
    "ca_pumping": {"rate_per_h": 2.5e9, "atp_per_event": 1.0},
}
DEMAND_PER_H = 1.0e10  # ATP per hour, the sum of the four

#: the priority order the program declares, least sensitive first (served first), from the abstract:
#: calcium cycling, then sodium cycling, then macromolecule synthesis (RNA/DNA before protein, an order
#: the abstract does not make inside "macromolecule biosynthesis"; the clause tests the groups only)
PRIORITY_ORDER = ("ca_pumping", "na_k_pumping", "rna_dna_synthesis", "protein_synthesis")
MACROMOLECULE = ("protein_synthesis", "rna_dna_synthesis")

GLUCOSE_UPTAKE_PER_H = 2.0e9  # the glucose pool's `regenerates`, the uptake ceiling
O2_SUPPLY_PER_H = 3.0e9  # the O2 pool's `regenerates`
ATP_STOCK = 3.0e9  # the ATP pool's size (Milo 2013, as stage 2's demo states it)
MITO_COPIES = 1000
CAPACITY_PER_COPY = 2.0e7  # electron pairs per hour per mitochondrial copy (`oxphos_capacity_per_copy`)
CAPACITY = MITO_COPIES * CAPACITY_PER_COPY  # 2e10 electron pairs per hour
DEMO_HETEROPLASMY = 0.7  # §6's own example, m.3243A>G = 0.7: below either arm's threshold

#: the energy-supply sweep, fixed before the run: oxygen lowered at full glucose, then glucose lowered
#: at no oxygen. Each step is (fraction of O2_SUPPLY_PER_H, fraction of GLUCOSE_UPTAKE_PER_H).
SUPPLY_SWEEP: tuple[tuple[float, float], ...] = (
    (1.0, 1.0),
    (0.75, 1.0),
    (0.5, 1.0),
    (0.25, 1.0),
    (0.1, 1.0),
    (0.0, 1.0),
    (0.0, 0.5),
    (0.0, 0.0),
)

#: registered, by hand: ATP supply per hour at each step of SUPPLY_SWEEP, per arm
PREDICTED_SUPPLY = {
    "high": (1e10, 1e10, 1e10, 7.75e9, 5.5e9, 4e9, 2e9, 0.0),
    "low": (1e10, 1e10, 1e10, 7.25e9, 5.3e9, 4e9, 2e9, 0.0),
}
#: glucose uptake per hour at each step
PREDICTED_GLUCOSE = {
    "high": (3.125e8, 3.125e8, 1.25e9, 2e9, 2e9, 2e9, 1e9, 0.0),
    "low": (1e10 / 28.0, 1e10 / 28.0, 1.75e9, 2e9, 2e9, 2e9, 1e9, 0.0),
}
#: the share of ATP made by OXPHOS at each step (0 where nothing is made)
PREDICTED_OXPHOS_SHARE = {
    "high": (0.875, 0.875, 0.7, 3.5e9 / 7.75e9, 1.4e9 / 5.5e9, 0.0, 0.0, 0.0),
    "low": (24.0 / 28.0, 24.0 / 28.0, 0.6, 3.0e9 / 7.25e9, 1.2e9 / 5.3e9, 0.0, 0.0, 0.0),
}
#: priority factors at the anoxic step with full glucose (supply 4e9 in both arms): exact
PREDICTED_ANOXIC_FACTORS = {
    "ca_pumping": 1.0,
    "na_k_pumping": 0.6,
    "rna_dna_synthesis": 0.0,
    "protein_synthesis": 0.0,
}

#: with uptake not limiting, anoxia multiplies glucose uptake by Y_aerobic / Y_anaerobic
UNCAPPED_GLUCOSE_PER_H = 1.0e11
PREDICTED_PASTEUR_RATIO = {"high": 16.0, "low": 14.0}

#: the heteroplasmy sweep at full oxygen and glucose, fixed before the run; h* lies between two points
#: in each arm and on neither
HETEROPLASMY_SWEEP = (0.0, 0.2, 0.4, 0.6, 0.7, 0.75, 0.78, 0.8, 0.82, 0.85, 0.9, 0.95, 1.0)
PREDICTED_THRESHOLD = {"high": 0.8125, "low": 1.0 - 12.0 * DEMAND_PER_H / (28.0 * CAPACITY)}  # 0.785714...

#: stage 4's resource checkpoint: its demo's own division cost, against this cell over one interval
DIVISION_COST_ATP = 1.0e10  # data/demo/stage4_division.bio, `cost: ATP 1e10 per division`
INTERVAL_H = 24.0
PREDICTED_MAX_SUPPLY = {"high": 1.9e10, "low": 1.7e10}  # ATP per hour at full O2 and glucose
PREDICTED_AVAILABLE = {"high": 2.19e11, "low": 1.71e11}  # stock + (max supply - demand) * 24 h

# ---- the glycolysis-only cell (data/demo/stage3_erythrocyte.bio) ----------------------------------

RBC_CONSUMERS: dict[str, dict[str, float]] = {
    "na_k_pumping": {"rate_per_h": 1.0e8, "atp_per_event": 1.0},
    "ca_pumping": {"rate_per_h": 1.0e8, "atp_per_event": 1.0},
}
RBC_DEMAND_PER_H = 2.0e8
RBC_GLUCOSE_UPTAKE_PER_H = 5.0e8
BPG_SHUNT = 0.19  # Mulquiney, Bubb & Kuchel 1999, Biochem J 342:567, abstract
PREDICTED_RBC_YIELD = 1.62
PREDICTED_RBC_GLUCOSE = RBC_DEMAND_PER_H / PREDICTED_RBC_YIELD  # 1.2345679e8 per hour
PREDICTED_RBC_LACTATE = 2.0 * PREDICTED_RBC_GLUCOSE
PREDICTED_RBC_SHUNT_FLUX = 2.0 * BPG_SHUNT * PREDICTED_RBC_GLUCOSE  # 2,3-BPG made per hour

#: relative tolerance where a run is compared with a closed form: both sides are the same arithmetic
#: in floating point, so anything looser would hide a real difference
CLOSED_FORM_RTOL = 1e-9

PRE_REGISTRATION: dict[str, Any] = {
    "written": (
        "2026-10-10, before genomeos/runtime/metabolism.py, bio.std.core_metabolism, the two demo "
        "programs and scripts/metabolism_gate.py existed and before any number was produced; committed "
        "in this file before any of them"
    ),
    "decision_check": (
        "§10 decision 4 (whose cost accounting) is the one open decision stage 3 touches, and §10 "
        "records it as measurable BY stage 3's ATP budget gate, not as a precondition; ROADMAP-HISTORY "
        "item 5: 'No decision blocks code any more'. Stage 3 is built; decision 4 stays open unless "
        "the magnitude clause is askable, and it is not (not_askable.magnitude)"
    ),
    "claim": (
        "a pool regenerated `from Glycolysis, TCA, OXPHOS` is refilled per hour by a reduced module "
        "whose every stoichiometric number is the program's; supply follows demand up to what glucose, "
        "oxygen and mitochondrial capacity allow; heteroplasmy removes the capacity of the copies that "
        "carry a declared variant; and what the cell cannot supply is divided among its consumers by "
        "the allocation it names"
    ),
    "counter_claims_registered_in_advance": (
        "(1) the hierarchy of Buttgereit & Brand cannot appear under `proportional` or `competitive`, "
        "which hand every consumer one scalar; under `priority` it appears because it is written in. "
        "The verdict will say 'reproduced as an input, not predicted'. (2) The Pasteur clause holds by "
        "construction of the routing; its passing is not evidence about cells. (3) The threshold's "
        "position is the program's declared reserve, and with no reserve there is none. (4) Every cell "
        "here collapses without glucose, because the reduced module has no other fuel"
    ),
    "design": (
        f"{PROGRAM}: a nucleated cell with {MITO_COPIES} mitochondrial copies of capacity "
        f"{CAPACITY_PER_COPY:.0e} electron pairs per hour each, four ATP consumers at an equal "
        f"{DEMAND_PER_H / 4:.1e} ATP per hour (sum {DEMAND_PER_H:.0e}), glucose uptake ceiling "
        f"{GLUCOSE_UPTAKE_PER_H:.0e}/h, oxygen {O2_SUPPLY_PER_H:.0e}/h, ATP stock {ATP_STOCK:.0e}, "
        f"priority order {', '.join(PRIORITY_ORDER)}, heteroplasmy m.3243A>G = {DEMO_HETEROPLASMY}. "
        f"{RBC_PROGRAM}: no mitochondrion, two pumps at {RBC_DEMAND_PER_H / 2:.0e} ATP per hour each, "
        f"glucose ceiling {RBC_GLUCOSE_UPTAKE_PER_H:.0e}/h, 2,3-BPG shunt {BPG_SHUNT}. Both import "
        "bio.std.core_metabolism. Arms: high and low ends of the yields' spread (ARMS). Sweeps: "
        "SUPPLY_SWEEP and HETEROPLASMY_SWEEP, fixed here"
    ),
    "clauses": {
        "1_yields": (
            f"per glucose, from the module's constants: aerobic {PREDICTED_YIELD_AEROBIC['high']} (high) "
            f"and {PREDICTED_YIELD_AEROBIC['low']} (low) exactly, {PREDICTED_YIELD_SPEC_LOW_END} with P/O "
            f"2.5 and 1.5 through the glycerol 3-phosphate shuttle, anaerobic "
            f"{PREDICTED_YIELD_ANAEROBIC}, and {PREDICTED_O2_PER_GLUCOSE} O2 for "
            f"{PREDICTED_PAIRS_PER_GLUCOSE} electron pairs in both arms (the combustion stoichiometry)"
        ),
        "2_budget": (
            f"at full oxygen and glucose: supply equals the demand {DEMAND_PER_H:.0e}/h exactly, glucose "
            "uptake is demand / aerobic yield, O2 uptake six times glucose, lactate 0, and the OXPHOS "
            f"share 28/32 (high) and 24/28 (low), to a relative {CLOSED_FORM_RTOL:g}"
        ),
        "3_pasteur": (
            "along SUPPLY_SWEEP's oxygen steps the OXPHOS share never rises and reaches exactly 0 with no "
            "oxygen, glucose uptake never falls, and supply, glucose uptake and OXPHOS share equal "
            "PREDICTED_SUPPLY, PREDICTED_GLUCOSE and PREDICTED_OXPHOS_SHARE at every step in both arms"
        ),
        "4_pasteur_uncapped": (
            f"with glucose uptake raised to {UNCAPPED_GLUCOSE_PER_H:.0e}/h, anoxia leaves supply equal to "
            "demand and multiplies glucose uptake by exactly "
            f"{PREDICTED_PASTEUR_RATIO['high']} (high) and {PREDICTED_PASTEUR_RATIO['low']} (low)"
        ),
        "5_glycolysis_only": (
            f"the red cell: ATP from TCA and OXPHOS exactly 0; yield {PREDICTED_RBC_YIELD}; glucose uptake "
            f"{PREDICTED_RBC_GLUCOSE:.7e}/h, lactate twice that, 2,3-BPG made at {2 * BPG_SHUNT} per "
            "glucose; supply equal to demand and every reading identical, to the bit, at every oxygen "
            "level of SUPPLY_SWEEP (falsifier: a glycolysis that read oxygen); with the shunt set to 0, "
            f"glucose uptake falls by exactly {PREDICTED_RBC_YIELD / 2}"
        ),
        "6_glucose_removed": (
            "with glucose uptake 0 the red cell's supply is exactly 0 and every consumer's factor is 0 "
            "under its allocation: the collapse. The nucleated cell also supplies exactly 0, reported as "
            "the module's limitation (no other fuel), not as the gate's contrast"
        ),
        "7_hierarchy": (
            "(a) under the program's priority order, along SUPPLY_SWEEP, the first step at which each "
            "consumer's factor falls below 1 comes for both macromolecule consumers before it comes for "
            "na_k_pumping, and for na_k_pumping before ca_pumping, in both arms; at the anoxic step with "
            "full glucose the factors are exactly PREDICTED_ANOXIC_FACTORS. (b) The same program with "
            "the allocation set to proportional, and to competitive: at every step the four factors are "
            "identical to one another. (c) Hence the verdict: reproduced as an input, not predicted"
        ),
        "8_threshold": (
            "along HETEROPLASMY_SWEEP at full oxygen and glucose: OXPHOS ATP equals its h = 0 value "
            "exactly at every h below the arm's PREDICTED_THRESHOLD, is strictly lower at every h above "
            "it, glucose uptake is strictly higher above it, and at h = 1 supply equals the anoxic "
            "step's. Falsifier: with capacity set to exactly what the demand needs (12 * D / Y), every "
            "h > 0 lowers OXPHOS ATP — no plateau, so the plateau is the reserve's"
        ),
        "9_resource_checkpoint": (
            f"stage 4's demo division ({DIVISION_COST_ATP:.0e} ATP) against this cell: refused under the "
            f"standing-stock reading (capacity = pool size {ATP_STOCK:.0e}, as §9.4 found), paid when the "
            f"pool is regenerated over {INTERVAL_H:g} h (available {PREDICTED_AVAILABLE['high']:.3g} high, "
            f"{PREDICTED_AVAILABLE['low']:.3g} low), and refused again with glucose removed (available = "
            "the stock)"
        ),
        "10_refusals": (
            "compile errors: heteroplasmy on a compartment whose genome does not read chrM; a fraction "
            "outside 0..1; fractions summing above 1; a variant not written m.<position><ref>><alt>. "
            "Refused by the runtime, never defaulted: a process outside Glycolysis, TCA and OXPHOS; OXPHOS "
            "with no compartment reading chrM; TCA without OXPHOS; Glycolysis with no Glucose pool; "
            "OXPHOS with no O2 pool or no stated capacity; a stoichiometric constant the program does not "
            "state (the engine holds none)"
        ),
    },
    "not_askable": {
        "magnitude": (
            "whether the budget lands within the published order of magnitude: this project holds no "
            "measured per-cell ATP turnover for any cell type (the three abstracts give none), so the "
            "demo's consumer rates are declared and no clause compares them with a measurement. §10 "
            "decision 4, which this magnitude was to settle, therefore stays open"
        ),
        "split": (
            "whether the budget's split matches Buttgereit & Brand: the abstract states a hierarchy of "
            "sensitivity to supply, not shares; the shares are in the paper's body, a scan not held. "
            "§6's parenthesis 'protein synthesis and the Na⁺/K⁺-ATPase first' reads the hierarchy as a "
            "split, and the abstract does not support that reading: sodium cycling comes after "
            "macromolecule synthesis"
        ),
        "rbc_ranges": (
            "whether the red cell keeps ATP and 2,3-BPG within Mulquiney & Kuchel's ranges: the ranges "
            "are in the papers' bodies, not held, and a steady-state concentration is a kinetic quantity "
            "the reduced module (stoichiometry and capacities, no rate laws) cannot produce. What is held "
            "and used is the shunt's 19% of glycolytic carbon flux"
        ),
        "threshold_values": (
            "whether the threshold sits where Rossignol's measured thresholds sit: the review's abstract "
            "gives no value, and the position here is the program's declared reserve, not a prediction"
        ),
        "pasteur_magnitude": (
            "how far real cells raise glucose uptake under hypoxia: they also lower demand, and no "
            "measured Pasteur ratio is held; 16 and 14 are the module's ceiling, not a prediction"
        ),
        "warburg": (
            "aerobic glycolysis: the routing oxidises pyruvate first whenever oxygen and capacity allow, "
            "so a cell that ferments with oxygen present is not expressible in stage 3"
        ),
        "other_fuels": (
            "fatty acids, glutamine and ketone bodies are not in the reduced module, so 'removing glucose "
            "collapses a glycolysis-only cell' cannot be contrasted with a cell that survives on them"
        ),
        "proton_leak": (
            "the abstract's least sensitive process spends the proton gradient, not ATP, and the module "
            "does not model leak; that rung of the hierarchy is not asked"
        ),
    },
    "must_not_move": {
        "bio_test": (
            "the 47 files and 288 checks `bio test data/demo genomeos/std data/organisms` passed at "
            "d32fe44 pass unchanged; the new files add to the count and none is subtracted"
        ),
        "committed_programs": (
            "programs that write `regenerates: from Glycolysis, OXPHOS` with no mitochondrion "
            "(stage 2's and stage 4's demos) still compile: the refusal is the runtime's, when it is asked "
            "to run them"
        ),
        "stage_2": (
            "Economy with no capacity override behaves exactly as before: burden factors "
            "1.0,1.0,1.0,1.0,0.50,0.25,0.05; abundance MAE 0.952005"
        ),
        "stage_4": "runtime/division.py untouched: 8 of 8 clauses",
        "body_runtime": "runtime/body.py untouched: area E's figures stand",
    },
}
