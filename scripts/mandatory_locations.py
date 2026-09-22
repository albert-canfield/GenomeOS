# SPDX-License-Identifier: AGPL-3.0-or-later
"""What making locations mandatory would cost (BIOLANG-v0.4-ECONOMY.md §10, open decision 1).

Decision 1 is filed as a preference — "no experiment decides it" — but its price is measurable, and a
price is what the decision actually needs. A located program must give every gene and every protein a
compartment (§4.2). So over every committed program: how many entities would need one they do not
have, and for how many of those does a source we already hold supply it, rather than an author
inventing one? A fact we must invent is exactly what the evidence rules exist to prevent, so the
answer is the cost of mandatory, in facts.

The source is the packaged human proteome (`genomeos.lib.proteome`), whose `location` field is
UniProt's subcellular location. A gene is located by the compartment whose genome holds its
chromosome, so a gene with a locus needs no fact at all: chrM genes are mitochondrial and the rest
nuclear.

    uv run python scripts/mandatory_locations.py
"""

from collections import Counter
from pathlib import Path

from genomeos.lang import parse_file
from genomeos.results import save_result

# UniProt's vocabulary against the compartments v0.4 programs declare; the left-hand side is what the
# packaged proteome says, and a location that maps to nothing is a fact the author would have to invent.
KNOWN = {
    "Cytoplasm": "Cytosol",
    "Cytosol": "Cytosol",
    "Nucleus": "Nucleus",
    "Mitochondrion": "Mitochondrion",
    "Mitochondrion inner membrane": "Mitochondrion",
    "Mitochondrion matrix": "Mitochondrion",
    "Mitochondrion outer membrane": "Mitochondrion",
    "Endoplasmic reticulum": "ER",
    "Endoplasmic reticulum membrane": "ER",
    "Cell membrane": "PlasmaMembrane",
    "Membrane": "PlasmaMembrane",
    "Secreted": "Extracellular",
}


def symbol_of(entity_id: str) -> str:
    """Programs name a protein after its gene with a `p` suffix (SDHB -> SDHBp)."""
    return entity_id[:-1] if entity_id.endswith("p") and len(entity_id) > 1 else entity_id


def main() -> None:
    from genomeos.lib.proteome import load

    proteome = load()
    paths = sorted(Path("data").rglob("*.bio")) + sorted(Path("genomeos/std").rglob("*.bio"))
    rows, totals = [], Counter()
    for path in paths:
        try:
            module = parse_file(path)
        except Exception:  # noqa: BLE001  (a program that does not compile is not this measurement's business)
            continue
        genes, proteins = module.genes(), module.proteins()
        if not genes and not proteins:
            totals["programs unaffected (no genes or proteins)"] += 1
            continue
        need_gene = [g for g in genes if not g.location and g.locus is None]
        need_protein = [p for p in proteins if not p.location]
        # a gene with no locus but a human symbol we hold is a lookup, not an invention
        gene_lookup = [g for g in need_gene if g.id in proteome or symbol_of(g.id) in proteome]
        covered, mapped, unmapped, no_source, not_human = [], [], [], [], []
        for p in need_protein:
            record = proteome.get(symbol_of(p.id)) or proteome.get(p.id)
            if record is None:
                not_human.append(p.id)  # no human symbol, so no human source can ever locate it
                continue
            where = record.get("location") or []
            if not where:
                no_source.append(p.id)  # a human protein UniProt does not place
                continue
            covered.append(p.id)
            (mapped if any(w in KNOWN for w in where) else unmapped).append(p.id)
        rows.append(
            {
                "program": str(path),
                "located": module.located,
                "genes": len(genes),
                "proteins": len(proteins),
                "genes needing a fact": len(need_gene),  # no locus, so no chromosome to place them by
                "of those, a symbol we hold (a lookup)": len(gene_lookup),
                "proteins needing a fact": len(need_protein),
                "of those, UniProt has a location": len(covered),
                "of those, it maps to a declared compartment": len(mapped),
                "UniProt location we cannot map": sorted(unmapped)[:10],
                "human, but UniProt places it nowhere": sorted(no_source)[:10],
                "not a human symbol, so no source exists": sorted(not_human)[:10],
            }
        )
        totals["programs with genes or proteins"] += 1
        totals["located already"] += int(module.located)
        for key in ("genes needing a fact", "proteins needing a fact", "of those, UniProt has a location"):
            totals[key] += rows[-1][key]
        totals["genes that are a lookup, not an invention"] += len(gene_lookup)
        totals["of those, it maps to a declared compartment"] += len(mapped)
        totals["human proteins UniProt places nowhere"] += len(no_source) + len(unmapped)
        totals["proteins with no possible source (not human symbols)"] += len(not_human)

    width = max(len(r["program"]) for r in rows) if rows else 10
    print(f"{'program':<{width}} {'genes':>6} {'prot':>6} {'need':>6} {'UniProt':>8} {'mappable':>9}")
    for r in rows:
        print(
            f"{r['program']:<{width}} {r['genes']:>6} {r['proteins']:>6} "
            f"{r['genes needing a fact'] + r['proteins needing a fact']:>6} "
            f"{r['of those, UniProt has a location']:>8} "
            f"{r['of those, it maps to a declared compartment']:>9}"
        )
    print()
    for key, value in totals.items():
        print(f"  {key}: {value}")
    path = save_result(
        "mandatory_locations",
        {
            "decision": "BIOLANG-v0.4-ECONOMY.md §10 open decision 1: locations opt-in or mandatory",
            "source": "packaged human proteome (UniProt subcellular location) via genomeos.lib.proteome",
            "mapping": KNOWN,
            "totals": dict(totals),
            "programs": rows,
        },
    )
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
