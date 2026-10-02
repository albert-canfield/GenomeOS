#!/usr/bin/env python3
"""Leave-one-out and set-one-in over the closed candidate list registered for the
mesoderm band the Tbxt->SOX17 sign correction removed (5cbce26).

The list is closed: it is read from the registration JSON and this script refuses a
candidate the registration does not name. Nothing here edits data/demo/gastrulation.bio;
each variant is written to a temporary file built from the committed module text.

No number is chosen per candidate. Every added gene takes the `max` and `basal` that TBXT
and SOX17 already carry in the module, every added activating rule the strength, threshold
and hill of the module's `Nodal activates TBXT`, and every added inhibiting rule those of
its `Tbxt inhibits SOX2`. Nothing is swept and no proportion is fitted: the bound is
"mesoderm share >= 0.01, sustained over the last 10% of simulated time", not the CS7 range.
"""

from __future__ import annotations

import argparse
import itertools
import math
import tempfile
import warnings
from pathlib import Path

from genomeos import manifest as mf
from genomeos.lang import parse_file
from genomeos.results import save_result
from genomeos.runtime.gastrulation import CENSUS_MODEL_RUN
from genomeos.runtime.grn import NetworkRuntime, UnresolvedModelWarning

ENTRY = "scripts/mesoderm_diagnosis.py"
OWN_CODE = (ENTRY, "tests/test_mesoderm_diagnosis.py")

#: THE REGISTRATION. It lives in committed code, not in a data file, so the closed list provably
#: cannot grow between the single-component pass and the pairwise pass without a commit saying so.
REGISTRATION = {
    "registered": "2026-10-02",
    "lane": "lane-mesoderm",
    "question": "The Tbxt->SOX17 sign correction at 5cbce26 made the gastrulation model match its "
    "literature source (Lolas et al. 2014, PNAS 111:4478, Fig. 1A ChIP-seq and Fig. 3A "
    "knockdown) and the mesoderm band went 0.10 to 0.00. The asymmetry is that the program "
    "became more faithful and lost the layer the measured embryo is mostly made of, so "
    "something in the program was being compensated for by the wrong sign: another rule wrong "
    "in a cancelling direction, a required species or input absent, or runtime semantics the "
    "existing tests did not reach.",
    "what_this_is_not": [
        "Not a tuning. 8bb9123 and 19423bd record that the module's expected_* "
        "proportions cite no source and that tuning two free numbers to meet them was "
        "REFUSED; that refusal is inherited. No proportion is fitted and no free number "
        "is adjusted to approach the census.",
        "Not a reversal. The sign correction at 5cbce26 stands. Tbxt activates SOX17 in "
        "every run below. A rule matching its source does not yield to a proportion.",
        "Not a validation. The CS7 comparison stays falsified on all three layers "
        "whatever is found here, with d55cb19's bound carried verbatim: 'a sampled "
        "census of one embryo is not a 1-D axis, and the model names no stage.'",
        "Not a threshold fix. d55cb19 records that no measured NODAL threshold maps to "
        "model units and that Dubrulle et al. 2015 argue induction kinetics rather than "
        "a clamped level set fates. A threshold the model cannot source is a finding, "
        "not a fix.",
        "Not a discovery. Whichever candidate restores the band, the result says 'the "
        "program lacked X', never 'X specifies mesoderm'. This diagnoses a program.",
    ],
    "bound": {
        "statement": "mesoderm share >= 0.01 at the program's final stage, sustained over the last "
        "10% of simulated time, under the default numbers",
        "share": 0.01,
        "last_fraction": 0.1,
        "cells": 120,
        "cells_needed": 2,
        "samples_in_tail": 5,
        "deliberately_weak": "Clearing the bound means 'mesoderm becomes non-zero at all'. It is "
        "NOT the measured CS7 human range 0.694-0.773 for mesoderm; matching "
        "that range would be tuning to the target and is not the goal.",
    },
    "run": {
        "numbers": "genomeos.runtime.gastrulation.CENSUS_MODEL_RUN, unchanged: 120 cells, 40 h, dt "
        "0.05, nodal_max 6.0, decay_length 0.35",
        "record_every": 20,
        "fate": "argmax over the final Sox2 / Tbxt / Sox17 protein levels, as the program already decides it",
    },
    "number_convention": {
        "why": "so that no number is chosen per candidate and nothing is swept",
        "added_gene": "max: 10; basal: 0.05 -- the values TBXT and SOX17 already carry "
        "in data/demo/gastrulation.bio",
        "added_protein": "half_life: 2 -- the value all four proteins already carry",
        "added_activating_rule": "strength: 1.0; threshold: 1.0; hill: 2 -- the "
        "module's own 'Nodal activates TBXT'",
        "added_inhibiting_rule": "strength: 1.0; threshold: 2.0; hill: 3 -- the "
        "module's own 'Tbxt inhibits SOX2'",
        "not_done": "no strength, threshold, hill, basal, max, nodal_max or "
        "decay_length is varied in any run",
    },
    "citation_verification": {
        "network": "none; this lane is free, 0 model requests, no money, no network",
        "consequence": "every citation below is given from the lane's own knowledge "
        "and was NOT re-fetched or re-read this session. Volume and "
        "page numbers are stated as recalled and are not verified. "
        "Figure-level claims are made only where the brief or the "
        "module already carried them (Lolas et al. 2014) and are "
        "otherwise stated at the level of the paper's result, not a "
        "panel.",
        "in_repo_check": "the repository holds no regulatory citation for "
        "WNT3A->TBXT, FGF->TBXT, MIXL1 or EOMES: "
        "genomeos/std/signalling.bio carries WNT3A only as "
        "ligand-receptor pairs (CellPhoneDB v5), and "
        "genomeos/lib/catalog.py carries EOMES, MIXL1 and WNT3A "
        "only as marker-gene names. So the sources could not be "
        "confirmed against anything already committed.",
    },
    "closed_list": [
        {
            "id": "W1",
            "name": "the program lacks a canonical-Wnt drive on TBXT",
            "kind": "missing species and rule",
            "citation": "Yamaguchi TP, Takada S, Yoshikawa Y, Wu N, McMahon AP (1999) T "
            "(Brachyury) is a direct target of Wnt3a during paraxial mesoderm "
            "specification. Genes Dev 13:3185-3190.",
            "acts_in_mesoderm_specification": "Wnt3a is required for paraxial mesoderm and T is "
            "a direct Wnt3a target there.",
            "species": [
                {
                    "gene": "WNT3A",
                    "protein": "Wnt3a",
                    "note": "mesoderm candidate W1/W2; level unsourced, module convention",
                }
            ],
            "rules": [
                {
                    "src": "Wnt3a",
                    "dst": "TBXT",
                    "action": "activates",
                    "cite": "Yamaguchi et al. 1999, Genes Dev 13:3185",
                }
            ],
            "predicted": "raises Tbxt; alone the added gene has no activator, so the runtime "
            "drives it at full max_rate and the drive is not gradient-shaped",
            "runnable": True,
        },
        {
            "id": "W2",
            "name": "the program lacks the Brachyury limb that closes the Wnt loop",
            "kind": "missing species and rule",
            "citation": "Martin BL, Kimelman D (2008) Regulation of canonical Wnt signaling by "
            "Brachyury is essential for posterior mesoderm formation. Dev Cell "
            "15:121-133.",
            "acts_in_mesoderm_specification": "Brachyury maintains canonical Wnt signalling and "
            "that autoregulation is required for posterior "
            "mesoderm.",
            "species": [
                {
                    "gene": "WNT3A",
                    "protein": "Wnt3a",
                    "note": "mesoderm candidate W1/W2; level unsourced, module convention",
                }
            ],
            "rules": [
                {
                    "src": "Tbxt",
                    "dst": "WNT3A",
                    "action": "activates",
                    "cite": "Martin & Kimelman 2008, Dev Cell 15:121",
                }
            ],
            "predicted": "alone, nothing reads Wnt3a, so the fates cannot move: a null control "
            "the pass must reproduce",
            "runnable": True,
        },
        {
            "id": "F1",
            "name": "the program lacks an FGF drive on TBXT",
            "kind": "missing species and rule",
            "citation": "Isaacs HV, Pownall ME, Slack JMW (1994) eFGF regulates Xbra expression "
            "during Xenopus gastrulation. EMBO J 13:4469-4481; Ciruna B, Rossant J "
            "(2001) FGF signaling regulates mesoderm cell fate specification and "
            "morphogenetic movement at the primitive streak. Dev Cell 1:37-49.",
            "acts_in_mesoderm_specification": "FGF signalling is required to maintain Brachyury "
            "expression and Fgfr1 loss misspecifies mesoderm "
            "at the streak.",
            "species": [
                {
                    "gene": "FGF4",
                    "protein": "Fgf4",
                    "note": "mesoderm candidate F1/F2; level unsourced, module convention",
                }
            ],
            "rules": [
                {
                    "src": "Fgf4",
                    "dst": "TBXT",
                    "action": "activates",
                    "cite": "Isaacs et al. 1994, EMBO J 13:4469; Ciruna & Rossant 2001, Dev Cell 1:37",
                }
            ],
            "predicted": "same shape as W1: a second TBXT activator the runtime drives at full max_rate",
            "runnable": True,
        },
        {
            "id": "F2",
            "name": "the program lacks the Brachyury limb that closes the FGF loop",
            "kind": "missing species and rule",
            "citation": "Schulte-Merker S, Smith JC (1995) Mesoderm formation in response to "
            "Brachyury requires FGF signalling. Curr Biol 5:62-67; Isaacs HV, "
            "Pownall ME, Slack JMW (1994) EMBO J 13:4469-4481.",
            "acts_in_mesoderm_specification": "Brachyury's mesoderm-inducing activity requires "
            "FGF signalling, and Brachyury induces eFGF, "
            "closing the loop.",
            "species": [
                {
                    "gene": "FGF4",
                    "protein": "Fgf4",
                    "note": "mesoderm candidate F1/F2; level unsourced, module convention",
                }
            ],
            "rules": [
                {
                    "src": "Tbxt",
                    "dst": "FGF4",
                    "action": "activates",
                    "cite": "Schulte-Merker & Smith 1995, Curr Biol 5:62",
                }
            ],
            "predicted": "alone, nothing reads Fgf4: a second null control",
            "runnable": True,
        },
        {
            "id": "M1",
            "name": "the program lacks a MIXL1 mesendoderm node",
            "kind": "missing species and rules",
            "citation": "Hart AH, Hartley L, Sourris K, Stadler ES, Li R, Stanley EG, Tam PPL, "
            "Elefanty AG, Robb L (2002) Mixl1 is required for axial mesendoderm "
            "morphogenesis and patterning in the murine embryo. Development "
            "129:3597-3608.",
            "acts_in_mesoderm_specification": "Mixl1-null embryos fail axial mesendoderm "
            "morphogenesis and patterning, so Mixl1 acts in "
            "the mesendoderm split.",
            "species": [
                {
                    "gene": "MIXL1",
                    "protein": "Mixl1",
                    "note": "mesoderm candidate M1; level unsourced, module convention",
                }
            ],
            "rules": [
                {
                    "src": "Nodal",
                    "dst": "MIXL1",
                    "action": "activates",
                    "cite": "Hart et al. 2002, Development 129:3597",
                },
                {
                    "src": "Mixl1",
                    "dst": "SOX2",
                    "action": "inhibits",
                    "cite": "Hart et al. 2002, Development 129:3597",
                },
            ],
            "predicted": "weakens the ectoderm competitor where NODAL is present, without "
            "touching TBXT's own drive",
            "runnable": True,
        },
        {
            "id": "E1",
            "name": "the program lacks an EOMES node, so SOX17's endoderm drive has no other carrier",
            "kind": "missing species and rules",
            "citation": "Arnold SJ, Hofmann UK, Bikoff EK, Robertson EJ (2008) Pivotal roles for "
            "eomesodermin during axis formation, epithelium-to-mesenchyme transition "
            "and endoderm specification in the mouse. Development 135:501-511; Teo "
            "AKK, Arnold SJ, Trotter MWB, Brown S, Ang LT, Chng Z, Robertson EJ, "
            "Dunn NR, Vallier L (2011) Pluripotency factors regulate definitive "
            "endoderm specification through eomesodermin. Genes Dev 25:238-250; "
            "Costello I, Pimeisl IM, Drager S, Bikoff EK, Robertson EJ, Arnold SJ "
            "(2011) The T-box transcription factor Eomesodermin acts upstream of "
            "Mesp1 to specify cardiac mesoderm during mouse gastrulation. Nat Cell "
            "Biol 13:1084-1091.",
            "acts_in_mesoderm_specification": "Eomes is required for the EMT that moves cells "
            "out of the streak and acts upstream of Mesp1 to "
            "specify cardiac mesoderm, while being biased to "
            "definitive endoderm; it is a candidate because it "
            "is a second carrier of the endoderm drive the "
            "model routes entirely through SOX17.",
            "species": [
                {
                    "gene": "EOMES",
                    "protein": "Eomes",
                    "note": "mesoderm candidate E1; level unsourced, module convention",
                }
            ],
            "rules": [
                {
                    "src": "Nodal",
                    "dst": "EOMES",
                    "action": "activates",
                    "cite": "Arnold et al. 2008, Development 135:501",
                },
                {
                    "src": "Eomes",
                    "dst": "SOX17",
                    "action": "activates",
                    "cite": "Teo et al. 2011, Genes Dev 25:238",
                },
            ],
            "predicted": "raises endoderm, not mesoderm; registered because the brief names it "
            "and because every entry is reported whether or not it helps",
            "runnable": True,
        },
        {
            "id": "A1",
            "name": "the program lacks a NODAL antagonist",
            "kind": "missing species, registered and NOT RUNNABLE",
            "citation": "Perea-Gomez A, Vella FDJ, Shawlot W, Oulad-Abdelghani M, Chazaud C, "
            "Meno C, Pfister V, Chen L, Robertson E, Hamada H, Behringer RR, Ang SL "
            "(2002) Nodal antagonists in the anterior visceral endoderm prevent the "
            "formation of multiple primitive streaks. Dev Cell 3:745-756.",
            "acts_in_mesoderm_specification": "Cer1 and Lefty1 restrict where the primitive "
            "streak and therefore nascent mesoderm form.",
            "runnable": False,
            "no_run_reason": "Cer1 and Lefty1 antagonise the NODAL ligand, and in this program "
            "NODAL is not a modelled species: its gene declares max 0 and the "
            "runtime holds Nodal at a clamped external level in every "
            "Runge-Kutta stage. No rule the module can state will move a "
            "clamped value, so the candidate cannot be set in without changing "
            "the program's input semantics. Reported as registered and "
            "unrunnable rather than implemented as an unsourced rule on SOX17.",
        },
        {
            "id": "R1",
            "name": "leave out Sox17 inhibits TBXT",
            "kind": "named rule already in the program",
            "citation": "Lolas M, Valenzuela PDT, Tjian R, Liu Z (2014) Charting "
            "Brachyury-mediated developmental pathways during early mouse "
            "embryogenesis. PNAS 111:4478-4483, Fig. 3C (Sox17 overexpression lowers "
            "Brachyury).",
            "acts_in_mesoderm_specification": "Sox17 overexpression lowers Brachyury, so the "
            "rule is a repression of the mesoderm program.",
            "remove": ["rule Sox17 inhibits TBXT"],
            "predicted": "the supervisor's fourth candidate: the corrected sign may have "
            "unmasked this repression without a balancing term",
            "runnable": True,
        },
        {
            "id": "R2",
            "name": "leave out Sox2 inhibits TBXT",
            "kind": "named rule already in the program",
            "citation": "Thomson M, Liu SJ, Zou LN, Smith Z, Meissner A, Ramanathan S (2011) "
            "Pluripotency factors in embryonic stem cells regulate differentiation "
            "into germ layers. Cell 145:875-889.",
            "acts_in_mesoderm_specification": "Sox2 level biases cells towards ectoderm and away "
            "from mesendoderm, so the rule is a repression of "
            "the mesoderm program.",
            "remove": ["rule Sox2 inhibits TBXT"],
            "runnable": True,
        },
        {
            "id": "R3",
            "name": "leave out Tbxt inhibits SOX2",
            "kind": "named rule already in the program",
            "citation": "Thomson M, Liu SJ, Zou LN, Smith Z, Meissner A, Ramanathan S (2011) "
            "Pluripotency factors in embryonic stem cells regulate differentiation "
            "into germ layers. Cell 145:875-889.",
            "acts_in_mesoderm_specification": "the same germ-layer choice: the mesendoderm "
            "program lowers Sox2.",
            "remove": ["rule Tbxt inhibits SOX2"],
            "predicted": "removes the only brake TBXT holds on its competitor, so mesoderm "
            "should fall or stay at zero",
            "runnable": True,
        },
    ],
    "closed_list_n": 10,
    "pairs": {
        "n": 10,
        "n_choose_2": 45,
        "runnable": 36,
        "not_runnable": 9,
        "not_runnable_reason": "the nine pairs containing A1, for A1's stated structural reason",
        "when": "the pairwise pass runs ONLY if no single component clears the bound, and every "
        "pair is reported",
        "why": "a two-component compensation is exactly what single-component runs miss, and it is "
        "this lane's leading hypothesis: the program has no positive feedback on TBXT at "
        "all, so a sourced loop needs both of its limbs",
        "list_may_not_grow": "the pairwise pass is over the SAME closed list; the list may not grow "
        "between passes or after any run",
    },
    "eliminated_in_advance": [
        {
            "item": "the activator combination rule",
            "shas": ["93caf61", "267cc88", "a6fc5c4"],
            "withdrawn_attribution": "3f7b0f5",
            "finding_carried": "mesoderm is absent under all four activator rules "
            "(mean, sum_capped, max, or): 0.633/0.000/0.367, "
            "0.558/0.000/0.442, 0.575/0.000/0.425, "
            "0.558/0.000/0.442; only SOX17 has two activators among "
            "the 41 hand-written programs, and 'the averaging did "
            "not remove the' band.",
            "runs_spent": 0,
        }
    ],
    "excluded": [
        {
            "id": "X1",
            "item": "leave out Nodal activates SOX17",
            "reason": "EXCLUDED by rule 2. Its citation, Kanai-Azuma et al. 2002 (Development "
            "129:2367), reports depletion of definitive gut endoderm in Sox17-null mice "
            "and says nothing about mesoderm specification. Kept as a characterisation "
            "run outside the closed list, labelled as not a candidate.",
            "excluded_from": "the closed candidate list; it cannot clear the bound",
        },
        {
            "id": "X2",
            "item": "the runtime's a = 1.0 for a gene with no activators",
            "reason": "EXCLUDED by rule 2: no paper about mesoderm specification can source a "
            "runtime convention. It is reported as a structural observation about the "
            "program, as arithmetic, not as a candidate that can clear the bound. Note "
            "that rule 1 lists runtime semantics as an allowed entry kind while rule 2 "
            "excludes it; the exclusion is taken as the stricter rule and is reported as "
            "a contradiction in the brief.",
            "excluded_from": "the closed candidate list; it cannot clear the bound",
        },
        {
            "id": "X3",
            "item": "the module's unsourced thresholds (Sox2->TBXT 4.0, Sox17->TBXT 2.0, Tbxt->SOX2 "
            "2.0, Tbxt->SOX17 4.0)",
            "reason": "EXCLUDED as an intervention: adjusting one is tuning a free number, refused "
            "at 8bb9123 and 19423bd. Per d55cb19, a threshold the model cannot source is "
            "a finding, not a fix. Reported as a finding if the diagnosis lands on one.",
            "excluded_from": "the closed candidate list",
        },
        {
            "id": "X4",
            "item": "the module's expected_ectoderm / expected_mesoderm / expected_endoderm",
            "reason": "EXCLUDED: not a mechanism, and they cite no source (8bb9123, 19423bd). No "
            "run touches them.",
            "excluded_from": "the closed candidate list",
        },
        {
            "id": "X5",
            "item": "nodal_max and decay_length",
            "reason": "EXCLUDED: 19423bd records that no gradient shape makes the expected "
            "proportions meetable, and the two are free numbers in model units. Varying "
            "either is tuning.",
            "excluded_from": "the closed candidate list",
        },
    ],
    "characterisation_runs": [
        {
            "id": "X1",
            "note": "outside the closed list, reported for the record only; cannot clear the bound",
            "remove": ["rule Nodal activates SOX17"],
        }
    ],
    "prior_observation": {
        "disclosed": "measured before this registration was written, on the UNMODIFIED "
        "committed module, with no candidate set in or left out",
        "what": "Tbxt's final protein level along the NODAL gradient never exceeds its "
        "competitors anywhere. At the clamped levels 6.0, 4.0, 3.0, 2.5, 2.0, "
        "1.5, 1.0, 0.5 the final Tbxt reads 0.244, 0.274, 0.414, 0.775, 1.206, "
        "0.256, 0.225, 0.215; its highest value, 1.206 at Nodal 2.0, stands "
        "against Sox2 5.667 and Sox17 4.013 at the same point.",
        "why_it_shaped_the_hypothesis": "TBXT's only activator is Nodal and both of its "
        "repressors run high, so the lane's leading "
        "hypothesis is that the program holds no "
        "positive feedback on TBXT at all. That is why "
        "the pairwise pass is registered in advance: a "
        "sourced autoregulatory loop has two limbs and "
        "neither limb alone is the loop.",
        "not_a_candidate_run": "no entry of the closed list was set in or left out to obtain these numbers",
    },
}
MODULE = Path("data/demo/gastrulation.bio")
#: the bound registered before any run: not the CS7 range, which stays falsified
BOUND_SHARE = 0.01
BOUND_LAST_FRACTION = 0.10
#: one sample per simulated hour, so the last 10% of a 40 h run holds five samples
RECORD_EVERY = 20
_LAYERS = {"ectoderm": "Sox2", "mesoderm": "Tbxt", "endoderm": "Sox17"}

#: numbers copied from the module, never chosen for a candidate
_EV = 'evidence: experimental "{cite}"; confidence: 0.4'
GENE_TEMPLATE = (
    "gene {gene} {{ max: 10; basal: 0.05; produces: {protein}; "
    'evidence: inferred "{note}"; confidence: 0.4 }}'
)
PROTEIN_TEMPLATE = 'protein {protein} {{ half_life: 2; evidence: inferred "{note}"; confidence: 0.4 }}'
ACTIVATES_TEMPLATE = "rule {src} activates {dst} {{ strength: 1.0; threshold: 1.0; hill: 2; " + _EV + " }}"
INHIBITS_TEMPLATE = "rule {src} inhibits {dst} {{ strength: 1.0; threshold: 2.0; hill: 3; " + _EV + " }}"


def variant_text(base: str, entries: list[dict]) -> str:
    """The module text with every entry's additions appended and removals deleted."""
    lines = base.splitlines()
    removals = [r for e in entries for r in e.get("remove", [])]
    for target in removals:
        hits = [i for i, line in enumerate(lines) if line.strip().startswith(target)]
        if len(hits) != 1:
            raise ValueError(f"removal {target!r} matched {len(hits)} lines, expected exactly 1")
        lines[hits[0]] = "# removed for this run: " + lines[hits[0]].strip()
    added: list[str] = []
    seen: set[str] = set()
    for e in entries:
        for sp in e.get("species", []):
            if sp["gene"] in seen:
                continue
            seen.add(sp["gene"])
            added.append(GENE_TEMPLATE.format(**sp))
            added.append(PROTEIN_TEMPLATE.format(note=sp["note"], protein=sp["protein"]))
        for rule in e.get("rules", []):
            tpl = ACTIVATES_TEMPLATE if rule["action"] == "activates" else INHIBITS_TEMPLATE
            added.append(tpl.format(**rule))
    return "\n".join(lines + added) + "\n"


def shares_over_time(path: Path, run: dict) -> tuple[list[float], list[dict[str, float]]]:
    """Fate shares at every recorded sample: one NetworkRuntime per cell, as the model does."""
    module = parse_file(path)
    per_cell: list[dict[str, list[float]]] = []
    times: list[float] = []
    for i in range(run["cells"]):
        x = (i + 0.5) / run["cells"]
        nodal = run["nodal_max"] * math.exp(-x / run["decay_length"])
        vm = NetworkRuntime(module)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UnresolvedModelWarning)
            traj = vm.run(
                hours=run["hours"],
                dt=run["dt"],
                initial={"Sox2": 1.0},
                record_every=RECORD_EVERY,
                clamp={"Nodal": nodal},
            )
        times = traj.times
        per_cell.append({sp: traj.levels[sp] for sp in _LAYERS.values()})
    out = []
    for k in range(len(times)):
        counts = dict.fromkeys(_LAYERS, 0)
        for cell in per_cell:
            counts[max(_LAYERS, key=lambda lay: cell[_LAYERS[lay]][k])] += 1
        out.append({lay: counts[lay] / run["cells"] for lay in _LAYERS})
    return times, out


def evaluate(path: Path, run: dict) -> dict:
    times, shares = shares_over_time(path, run)
    cut = run["hours"] * (1.0 - BOUND_LAST_FRACTION)
    tail = [s["mesoderm"] for t, s in zip(times, shares, strict=True) if t >= cut]
    return {
        "final": shares[-1],
        "mesoderm_tail_min": min(tail),
        "mesoderm_tail_samples": len(tail),
        "mesoderm_tail": tail,
        "clears_bound": bool(tail) and min(tail) >= BOUND_SHARE,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pass", dest="which", choices=("single", "pairs"), default="single")
    args = ap.parse_args()

    reg = REGISTRATION
    entries = {e["id"]: e for e in reg["closed_list"]}
    base = MODULE.read_text()
    run = dict(CENSUS_MODEL_RUN)
    results: dict[str, dict] = {}

    with tempfile.TemporaryDirectory() as tmp:

        def measure(name: str, chosen: list[dict]) -> dict:
            p = Path(tmp) / f"{name.replace('+', '_and_')}.bio"
            p.write_text(variant_text(base, chosen))
            return evaluate(p, run)

        results["baseline"] = measure("baseline", [])
        if args.which == "single":
            for cid, e in entries.items():
                if not e.get("runnable", True):
                    results[cid] = {"no_run": e["no_run_reason"]}
                    continue
                results[cid] = measure(cid, [e])
            for e in reg["characterisation_runs"]:
                results[e["id"] + "_characterisation"] = measure(e["id"], [e])
        else:
            for a, b in itertools.combinations(list(entries), 2):
                name = f"{a}+{b}"
                blocked = [c for c in (a, b) if not entries[c].get("runnable", True)]
                results[name] = (
                    {"no_run": entries[blocked[0]]["no_run_reason"]}
                    if blocked
                    else measure(name, [entries[a], entries[b]])
                )

    cleared = sorted(k for k, v in results.items() if v.get("clears_bound"))
    payload = {
        "pass": args.which,
        "run": run,
        "record_every": RECORD_EVERY,
        "bound": {
            "statement": reg["bound"]["statement"],
            "share": BOUND_SHARE,
            "last_fraction": BOUND_LAST_FRACTION,
            "not_the_census_range": reg["bound"]["deliberately_weak"],
        },
        "registration": reg,
        "runs": results,
        "cleared_the_bound": cleared,
        "cleared_count": len(cleared),
        "reported_count": len(results),
        "carried_verbatim": {
            "census_bound_d55cb19": (
                "a sampled census of one embryo is not a 1-D axis, and the model names no stage"
            ),
            "cs7_verdict": "falsified on all three layers, unchanged by anything in this file",
            "sign_correction": "5cbce26 stands: Tbxt activates SOX17 in every run here",
        },
        "wording": (
            "any entry that restores the band means the program lacked it; "
            "it does not mean it specifies mesoderm"
        ),
        "alphagenome_requests": 0,
        "money": "none: no network, no dataset read, one in-repo module text",
    }
    payload["result_manifest"] = {
        "sources": [
            {"accession": "data/demo/gastrulation.bio", "version": "in-repo demo module at this revision"},
            {
                "accession": "Lolas et al. 2014, PNAS 111:4478",
                "version": "the sign the module already states",
            },
        ],
        "inputs": [mf.input_entry(MODULE, partition=None)],
        "assembly": "n/a: a simulated regulatory network, no genome coordinates",
        "coordinates": "n/a: a normalised axis x in [0, 1], not genomic intervals",
        "parameters": {
            "pass": args.which,
            "model_run": run,
            "record_every": RECORD_EVERY,
            "bound_share": BOUND_SHARE,
            "bound_last_fraction": BOUND_LAST_FRACTION,
            "added_gene": REGISTRATION["number_convention"]["added_gene"],
            "added_activating_rule": REGISTRATION["number_convention"]["added_activating_rule"],
            "added_inhibiting_rule": REGISTRATION["number_convention"]["added_inhibiting_rule"],
        },
        "exclusions": [x["item"] + " -- " + x["reason"] for x in reg["excluded"]]
        + [
            "the activator combination rule: eliminated in advance by 93caf61, 267cc88, a6fc5c4; 0 runs spent"
        ],
        "partitions": "n/a: no evaluation split; every entry of the closed list is reported",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE),
    }
    path = save_result(f"mesoderm_diagnosis_{args.which}", payload)
    print(f"{args.which}: {len(results)} reported, cleared {len(cleared)}: {', '.join(cleared) or 'none'}")
    for k, v in results.items():
        if "no_run" in v:
            print(f"  {k:18s} no run: {v['no_run'][:64]}")
        else:
            f = v["final"]
            print(
                f"  {k:18s} ecto {f['ectoderm']:.3f} meso {f['mesoderm']:.3f} endo {f['endoderm']:.3f}"
                f"  tail_min {v['mesoderm_tail_min']:.4f}  {'CLEARS' if v['clears_bound'] else '-'}"
            )
    print(f"  -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
