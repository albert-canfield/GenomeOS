# SPDX-License-Identifier: AGPL-3.0-or-later
"""Network knockouts as an experiment: hold nodes or edges of a network at zero and read the consequence.

The `experiment` block perturbs organisms; a network (a Boolean model or a rule module run by the GRN
runtime) had knockouts only inside hand-written tests. Here an experiment names nodes (a species held
at zero for the whole run) and edges (`A>B`: B reads A as zero, everything else reads A as it is), the
perturbed run is compared against the unperturbed one from the same start, and the size of the change
is ranked against knockouts of the same size drawn at random, so a consequence is read against chance.
"""

from __future__ import annotations

# Registered 2026-09-27, before the code below was written.
EDGE_SEPARATOR = ">"  # "A>B": the edge from A into B's rule
CONTROL_DRAWS = 200  # random matched knockouts; all of them when fewer matched sets exist
CONTROL_SEED = 0
# consequence: mean over the nodes knocked out in NEITHER run of |fraction of the attractor (Boolean) or
# of the second half of the run (continuous) spent ON / the mean level, knockout minus unperturbed|, on
# a 0..1 scale for Boolean and as |log2((ko + 1e-3) / (wt + 1e-3))| for continuous levels; plus whether
# the dynamics changed kind (fixed point <-> cycle; oscillating <-> not).
# chance: p = (1 + #matched random knockouts with a consequence >= the observed) / (1 + #draws). Matched
# = the same number of nodes and the same number of edges, drawn from nodes and edges not chosen.
NETWORK_KNOCKOUT_SPEC = (
    "run(network, knockouts) returns the unperturbed and the perturbed dynamics from the same start, "
    "the per-node change, whether the dynamics changed kind, and the consequence ranked against "
    "matched random knockouts that never include the chosen nodes or edges."
)
NETWORK_KNOCKOUT_FALSIFIER = (
    "A knockout of a node no rule reads must have consequence 0 and chance p = 1.0; if the control "
    "calls it exceptional, or ever draws a chosen node, the control is broken. And the known case "
    "below must hold."
)
# the known case, from an existing fixture: data/models/mammalian_cell_cycle.bnet (Faure, Naldi,
# Chaouiya & Thieffry 2006, Bioinformatics 22:e124) and tests/test_boolean.py, which pins that with
# CycD (growth factor) on the network cycles and with CycD off it rests in a G1 fixed point with Rb,
# p27 and Cdh1 on and CycA, CycB off.
NETWORK_KNOCKOUT_KNOWN_CASE = {
    "network": "data/models/mammalian_cell_cycle.bnet",
    "start": "CycD on, every other node off (the fixture's start)",
    "knockout": "CycD",
    "unperturbed": "a cyclic attractor (length > 1)",
    "perturbed": "a fixed point with Rb, p27, Cdh1 on and CycA, CycB off; kind changed",
    "source": "Faure et al. 2006, as pinned by tests/test_boolean.py::test_faure_mammalian_cell_cycle",
}
