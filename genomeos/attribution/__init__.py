# SPDX-License-Identifier: AGPL-3.0-or-later
"""The 98%: attributing a function, or the honest absence of one, to every UNKNOWN block.

The coding genome is decoded and the rest is classified by sequence class
(``genomeos.genome.unknown``). This package attaches *evidence about function* to
those blocks, starting with evolutionary constraint (Zoonomia phyloP over 241
mammals and the 100-vertebrate conserved elements), and turns the class plus the
evidence into a best guess with a confidence: structural, fossil, regulatory,
constrained-unknown or neutral. docs/ATTRIBUTION.md is the design document.
"""
