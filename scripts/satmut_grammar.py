# SPDX-License-Identifier: AGPL-3.0-or-later
"""Measured grammar: saturation-mutagenesis MPRA against conservation and motif sites.

    uv run python scripts/satmut_grammar.py

Downloads (once, git-ignored) the 44,658 GRCh38 measurements of Kircher et al. 2019, calls the bases
whose substitution changes activity, and asks what marks them: mammalian phyloP, JASPAR sites across a
sweep of thresholds, both together, and which factor families' sites carry effects. Result:
`satmut_grammar.json`. No AlphaGenome quota is used, so it runs beside the deletion-scoring chain.
"""

from __future__ import annotations

import sys

from genomeos.knowledge.satmut import run_and_save


def main() -> int:
    r = run_and_save(progress=lambda m: print(f"  {m}", flush=True))
    for key in ("functional", "strong"):
        p = r["pooled"][key]
        strict = p["sweep"][-1]
        print(
            f"{key}: {p['bases_functional']:,} of {p['bases_measured']:,} bases over {p['loci']} loci; "
            f"AUC phyloP {p['phylop_auc_median']}, motif score {p['motif_auc_median']}; "
            f"at threshold {strict['threshold']} sites cover {strict['bases_in_sites']:.0%} of bases "
            f"and are {strict['enrichment']}x enriched; conserved bases in a site "
            f"{p['conserved']['in_site']['share']:.1%} functional against "
            f"{p['conserved']['outside']['share']:.1%} outside; "
            f"families at q<=0.05: {len(p['families_q05'])} of {p['families_tested']}",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
