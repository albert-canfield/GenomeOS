# SPDX-License-Identifier: AGPL-3.0-or-later
"""Digital Development (Du et al. 2014, Cell 156:359; Du et al. 2015, Nucleic Acids Res 44:D781): 204
conserved essential genes knocked down, 1,368 embryos lineaged, and for each gene the founder cells whose
fate changed and which other founder's identity they adopted (the "Fate_regulation" table).

This is the published knockout set the C. elegans founder rules are checked against: for every gene of
the table that names a factor or a signal in the organism program, the program's own knockout is run and
its founder-level changes are compared with the observed transformations. The table is 17 KB and is
distilled to data/results; the comparison is recomputed from the current program each time.
"""

from __future__ import annotations

import urllib.request
from pathlib import Path

from genomeos.ir import Experiment, Module

from .experiment import run_experiment
from .human import read_xlsx

FATE_URL = "http://www.digital-development.org/Batch-download/Fate_regulation.xlsx"
SOURCE = "Du et al. 2014, Cell 156:359 (Digital Development, Fate_regulation.xlsx)"
# founder identity -> the precursor type the program gives a cell with that identity
IDENTITY_TYPE = {
    "ABa": "ABaPrecursor",
    "ABp": "ABpPrecursor",
    "EMS": "EMSPrecursor",
    "MS": "MSPrecursor",
    "E": "EPrecursor",
    "C": "CPrecursor",
    "D": "DPrecursor",
}
# table gene -> factor or signal component knocked out in the program
GENE_TO_FACTOR = {
    "pop-1": "POP-1",
    "skn-1": "SKN-1",
    "pie-1": "PIE-1",
    "pal-1": "PAL-1",
    "apx-1": "APX-1",
    "glp-1": "GLP-1",
    "mom-2": "MOM-2",
    "mom-5": "MOM-5",
    "par-2": "PAR-2",
    "par-3": "PAR-3",
    "lag-2": "LAG-2",
}


def fetch(url: str = FATE_URL, timeout: float = 120.0) -> bytes:
    with urllib.request.urlopen(url, timeout=timeout) as r:  # noqa: S310 (fixed public URL)
        return r.read()


def parse_table(archive: bytes) -> dict[str, dict[str, str]]:
    """gene -> {cell: identity it adopted} from the gene-centric sheet."""
    sheets = read_xlsx(archive)
    rows = sheets.get("Gene Centric View") or next(iter(sheets.values()))
    header = rows[0]
    out: dict[str, dict[str, str]] = {}
    for r in rows[1:]:
        if not r or not r[0]:
            continue
        changes = {header[i]: r[i].strip() for i in range(1, min(len(r), len(header))) if r[i].strip()}
        if changes:
            out[r[0].strip()] = changes
    return out


def _descendants(module_cells: dict, cell: str, depth: int = 2) -> list[str]:
    out, frontier = [], [cell]
    for _ in range(depth):
        nxt = []
        for c in frontier:
            kids = module_cells[c].children if c in module_cells else []
            nxt.extend(kids)
        out.extend(nxt)
        frontier = nxt
    return out


def compare(module: Module, table: dict[str, dict[str, str]], until: float = 300.0) -> dict:
    """Run the program's knockout for every gene the program models and score it against the table.

    An observed change `X adopts Y` counts as reproduced when the mutant's X (or, if X is not typed, a
    daughter of X) carries Y's precursor type while the wild type does not; a change whose identity the
    program has no type for (an 8-cell AB granddaughter) is reported as not modelled, not as a miss."""
    rows = []
    for gene, changes in sorted(table.items()):
        factor = GENE_TO_FACTOR.get(gene)
        if factor is None:
            continue
        ex = Experiment(name=gene, knockouts=[factor], until=until)
        r = run_experiment(module, ex)
        wt, mu = r.wild_type.cells, r.mutant.cells
        observed = []
        for cell, identity in changes.items():
            want = IDENTITY_TYPE.get(identity)
            if want is None or cell not in mu:
                observed.append({"cell": cell, "adopts": identity, "status": "not modelled"})
                continue
            candidates = [cell] + _descendants(mu, cell)
            hit = any(mu[c].cell_type == want and wt[c].cell_type != want for c in candidates if c in wt)
            observed.append({"cell": cell, "adopts": identity, "status": "reproduced" if hit else "missed"})
        predicted = [
            {"cell": ch.cell, "wild_type": ch.wild_type, "mutant": ch.mutant}
            for ch in r.early()
            if ch.wild_type != ch.mutant and ch.mutant in IDENTITY_TYPE.values()
        ]
        rows.append({"gene": gene, "factor": factor, "observed": observed, "predicted": predicted})
    reproduced = sum(1 for row in rows for o in row["observed"] if o["status"] == "reproduced")
    missed = sum(1 for row in rows for o in row["observed"] if o["status"] == "missed")
    not_modelled = sum(1 for row in rows for o in row["observed"] if o["status"] == "not modelled")
    return {
        "source": SOURCE,
        "url": FATE_URL,
        "genes_in_table": len(table),
        "genes_modelled": len(rows),
        "observed_changes": reproduced + missed + not_modelled,
        "reproduced": reproduced,
        "missed": missed,
        "not_modelled": not_modelled,
        "rate_over_modelled": round(reproduced / (reproduced + missed), 3) if reproduced + missed else None,
        "genes": rows,
        "unmodelled_genes": sorted(g for g in table if g not in GENE_TO_FACTOR),
    }


def distil(program: str | Path = "data/organisms/celegans/embryo.bio") -> dict:
    from genomeos.lang import parse_file

    table = parse_table(fetch())
    out = compare(parse_file(program), table)
    out["table"] = table
    return out
