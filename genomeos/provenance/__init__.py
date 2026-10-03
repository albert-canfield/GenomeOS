# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Number and evidence provenance of hand-authored BioLang programs: what each number RESTS ON.

Five censuses live here. Each reads programs and classifies the slots in them -- what the literature
holds for a quantity, what the program's own citation does for the declaration it is attached to,
whether a measured value could be expressed in the field at all -- and none of them changes a
strength, a threshold, a Hill coefficient, a basal or a max. They report; they decide nothing.

WHY THIS PACKAGE EXISTS, 2026-10-03. All five were written under `genomeos/lang/`, which
`scripts/package_engine.py` copies verbatim into a wheel whose pyproject declares **Apache-2.0**
(PACKAGES = lang, ir, runtime, std). Each file carried `SPDX-License-Identifier: AGPL-3.0-or-later`,
so the engine package shipped AGPL files inside an Apache package -- the defect
`tests/test_engine_licence_headers.py` pinned on 2026-10-02 and Albert resolved by moving them here.
A census of what a program's numbers rest on is application work, not language, IR, VM or standard
library: nothing in the toolchain imports it, and the engine never did.

Two consequences the move buys rather than patches. These modules may now import the application
directly -- `genomeos/attribution/measured.py` above all, which `rule_cell_provenance` had to take as
a parameter from `scripts/rule_cell_provenance_cost.py` because decision D40 forbade the import. And
their suites leave the packaged engine's test selection on their own, because
`package_engine.engine_tests()` picks tests by the imports it reads: a test importing
`genomeos.provenance` is no longer one that imports only the engine.
"""
