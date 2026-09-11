"""Structural context, and the epitope a mutation may create on the outside.

The highest-value target class in the whole pipeline is a tumour-specific
change that is physically exposed:

    normal protein   A-B-C-D-E
    tumour protein   A-B-X-D-E
                         ^ changed, and outside the cell

That is a binder's dream, because the wild-type protein on healthy cells
carries a different surface there. It is also rare, and the usual reason it
fails is that the changed residue is inside the cell, or buried, or in a
domain nobody has ever resolved.

So four separate questions are answered separately and never merged:
extracellular, surface-accessible, structurally resolved, predicted bindable.
Coordinates are never invented; PDB residue ranges come from the deposited
entries and the solvent accessibility is reported as not computed.
"""

from __future__ import annotations

import re
from typing import Any

from .evidence import Evidence, database, derived
from .model import Epitope, Localisation, Measure, Region, StructureContext

CHAIN_RANGE = re.compile(r"=(\d+)-(\d+)")


def _covered(entry: dict[str, Any]) -> tuple[int, int] | None:
    m = CHAIN_RANGE.search(str(entry.get("chains") or ""))
    return (int(m.group(1)), int(m.group(2))) if m else None


def _domains(defn: dict[str, Any]) -> list[Region]:
    items = (defn.get("sections", {}).get("domains") or {}).get("items") or {}
    return [
        Region("domain", f.get("start"), f.get("end"), f.get("description") or "")
        for f in items.get("features") or []
        if f.get("type") == "Domain"
    ]


def context(
    gene: str,
    defn: dict[str, Any],
    localisation: Localisation,
    residue: int | None,
    structures: dict[str, Any] | None,
) -> StructureContext:
    """Structures, domains and coverage around the altered residue."""
    sections = defn.get("sections", {})
    ident = (sections.get("identity") or {}).get("items") or {}
    exp = (structures or {}).get("experimental") or []
    pred = (structures or {}).get("predicted") or []
    ctx = StructureContext(
        uniprot_id=ident.get("accession"),
        predicted_structure_id=(pred[0].get("id") if pred else None),
        mutation_position=residue,
    )
    ctx.pdb_ids = [
        {
            "id": e.get("id"),
            "method": e.get("method"),
            "resolution_A": e.get("resolution_A"),
            "residue_range": list(_covered(e) or ()),
        }
        for e in sorted(exp, key=lambda x: (x.get("resolution_A") is None, x.get("resolution_A") or 99))[:12]
    ]

    ecto = localisation.extracellular_regions
    if ecto and exp:
        biggest = max(ecto, key=lambda r: (r.end or 0) - (r.start or 0))
        covering = [
            e
            for e in exp
            if (c := _covered(e)) and c[0] >= (biggest.start or 0) - 5 and c[1] <= (biggest.end or 10**9) + 5
        ]
        ctx.extracellular_structure_available = bool(covering)
        if covering:
            best = min(covering, key=lambda x: x.get("resolution_A") or 99)
            ctx.evidence.append(
                database(
                    "PDB (via UniProt cross-references)",
                    f"{len(covering)} experimental structures cover the extracellular region of {gene} "
                    f"({biggest.start}-{biggest.end}); best {best.get('id')} {best.get('method')}",
                    0.95,
                    "in_vitro",
                    identifier=best.get("id"),
                    url=f"https://www.ebi.ac.uk/pdbe/entry/pdb/{str(best.get('id', '')).lower()}",
                )
            )
        else:
            ctx.evidence.append(
                derived(
                    "GenomeOS structure",
                    f"{len(exp)} structures exist for {gene} but none is confined to the extracellular "
                    f"region {biggest.start}-{biggest.end}; a binder's epitope is not resolved",
                    0.6,
                )
            )
    elif ecto and not exp:
        ctx.extracellular_structure_available = False
        ctx.evidence.append(
            derived(
                "GenomeOS structure",
                f"{gene} has an extracellular region but no deposited experimental structure",
                0.5,
            )
        )

    domains = _domains(defn)
    if residue is not None:
        home = next((d for d in domains if d.contains(residue)), None)
        if home:
            ctx.target_domain = home.description
            ctx.residue_range = [home.start, home.end]
            ctx.evidence.append(
                database(
                    "UniProtKB/InterPro domains",
                    f"{gene} residue {residue} sits in the {home.description} domain "
                    f"({home.start}-{home.end})",
                    0.9,
                    "human",
                )
            )
        covering = [e for e in exp if (c := _covered(e)) and c[0] <= residue <= c[1]]
        if covering:
            best = min(covering, key=lambda x: x.get("resolution_A") or 99)
            ctx.evidence.append(
                database(
                    "PDB (via UniProt cross-references)",
                    f"residue {residue} of {gene} is resolved in {len(covering)} structures; best "
                    f"{best.get('id')} ({best.get('method')})",
                    0.95,
                    "in_vitro",
                    identifier=best.get("id"),
                )
            )
        elif exp:
            ctx.evidence.append(
                derived(
                    "GenomeOS structure",
                    f"residue {residue} of {gene} is not covered by any deposited structure",
                    0.6,
                )
            )
    if pred:
        ctx.evidence.append(
            Evidence(
                "AlphaFold DB",
                "prediction",
                f"predicted structure {pred[0].get('id')} is available for {gene}; a prediction is not an "
                "experimental structure and pLDDT is not solvent accessibility",
                "computational",
                0.7,
                identifier=pred[0].get("id"),
                url=f"https://alphafold.ebi.ac.uk/entry/{pred[0].get('id')}",
            )
        )
    ctx.surface_accessibility = Measure.unavailable(
        "solvent accessibility is not computed; it needs coordinates and a per-residue calculation"
    )
    return ctx


def mutation_epitope(
    gene: str,
    protein_change: str,
    residue: int | None,
    sequence: str | None,
    localisation: Localisation,
    ctx: StructureContext,
    topology_evidence: Evidence | None,
    flank: int = 7,
) -> Epitope | None:
    """A candidate mutation-specific epitope, with each claim kept separate."""
    if residue is None or not protein_change:
        return None
    outside = localisation.extracellular_residue(residue)
    ep = Epitope(
        epitope_type="mutation_specific_surface" if outside else "mutation_specific_internal",
        protein=gene,
        protein_position=residue,
        extracellular=outside,
        tumour_specific_change=protein_change,
    )
    if topology_evidence:
        ep.evidence.append(topology_evidence)
    if sequence and residue <= len(sequence) and len(protein_change) >= 3:
        wild, mutant = protein_change[0], protein_change[-1]
        if sequence[residue - 1] == wild:
            a, b = max(0, residue - 1 - flank), min(len(sequence), residue + flank)
            ep.reference_context = sequence[a:b]
            ep.mutant_context = sequence[a : residue - 1] + mutant + sequence[residue:b]
            ep.region = Region("epitope_window", a + 1, b, f"{flank} residues either side of {residue}")
    home = next(
        (r for r in localisation.extracellular_regions if r.contains(residue)),
        None,
    )
    if home and ctx.pdb_ids:
        resolved = any(
            r.get("residue_range") and r["residue_range"][0] <= residue <= r["residue_range"][1]
            for r in ctx.pdb_ids
        )
        ep.structurally_resolved = resolved
    elif ctx.pdb_ids:
        ep.structurally_resolved = any(
            r.get("residue_range") and r["residue_range"][0] <= residue <= r["residue_range"][1]
            for r in ctx.pdb_ids
        )
    ep.surface_accessible = None if outside else False
    ep.predicted_bindable = None
    ep.experimentally_validated = None
    ep.normal_human_match_risk = (
        "the wild-type protein carries the reference residue at this position on every healthy cell; a "
        "mutation-specific binder must discriminate a single-residue difference"
    )
    if outside:
        ep.evidence.append(
            derived(
                "GenomeOS epitope reasoning",
                f"{gene} {protein_change} changes a residue on the outside of the cell: a candidate "
                "mutation-specific surface epitope, subject to structural exposure being shown",
                0.7,
            )
        )
    elif outside is False:
        ep.evidence.append(
            derived(
                "GenomeOS epitope reasoning",
                f"{gene} {protein_change} changes a residue that faces the inside of the cell, so it "
                "cannot be recognised by a circulating binder however accessible the protein is",
                0.9,
            )
        )
    else:
        ep.evidence.append(
            derived(
                "GenomeOS epitope reasoning",
                f"topology of {gene} around residue {residue} is not curated; whether the change is "
                "exposed cannot be established",
                0.0,
            )
        )
    return ep


def structural_bindability(ctx: StructureContext, localisation: Localisation) -> tuple[float | None, str]:
    """How much structural information a binder designer would actually have."""
    if not localisation.reachable:
        return None, "protein is not surface-reachable, so an extracellular epitope does not arise"
    if ctx.extracellular_structure_available:
        return 0.9, "an experimental structure of the extracellular region is deposited"
    if ctx.pdb_ids:
        return 0.5, (
            f"{len(ctx.pdb_ids)} experimental structures exist but none is confined to the "
            "extracellular region"
        )
    if ctx.predicted_structure_id:
        return 0.3, "only a predicted structure is available; a prediction is not an experimental epitope"
    return None, "no experimental or predicted structure found"
