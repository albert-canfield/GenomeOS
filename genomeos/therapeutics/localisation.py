"""Where the protein is, and which part of it a binder could actually touch.

A mutated gene is not a surface target. RUNX1 is mutated in leukaemia and is
a nuclear transcription factor: no circulating binder reaches it. ERBB2 is a
single-pass receptor with 600 residues outside the cell. The difference is
topology, and topology is curated, so it is read rather than guessed.

The strongest statement available wins, and each one carries its own
confidence: UniProt topological domains (curated, per-residue) beat UniProt
subcellular-location text, which beats Human Protein Atlas immunofluorescence,
which beats a Gene Ontology cellular-component term.

The stage also answers a question that decides the whole recognition strategy:
is the *mutated residue itself* on the outside? A surface protein whose
mutation sits in the cytoplasmic tail cannot be targeted mutation-specifically
from outside, however accessible the protein is.
"""

from __future__ import annotations

from typing import Any

from .evidence import Evidence, database, derived
from .model import Localisation, Region

#: A "Cell membrane" annotation on its own does not say which side of the
#: membrane the protein is on. SRC sits on the inner leaflet on a myristoyl
#: anchor and is annotated exactly the same way as a receptor whose ectodomain
#: faces the blood. So a bare membrane annotation is held below the
#: reachability threshold until an external signal raises it.
AMBIGUOUS_MEMBRANE = 0.45

#: Confidence a reachable compartment needs before a circulating binder is said
#: to be able to engage the protein.
REACHABLE_THRESHOLD = 0.6

#: Lipid anchors that hold a protein against the cytoplasmic face.
INNER_LEAFLET_ANCHORS = ("myristoyl", "palmitoyl", "farnesyl", "geranylgeranyl", "prenyl")

#: UniProt subcellular-location phrases mapped onto the pipeline's compartments.
#: Ordered: the first match on a phrase wins, so specific phrases come first.
LOCATION_TERMS: tuple[tuple[str, str, float], ...] = (
    ("cell membrane", "plasma_membrane", AMBIGUOUS_MEMBRANE),
    ("plasma membrane", "plasma_membrane", AMBIGUOUS_MEMBRANE),
    ("cell surface", "cell_surface", 0.9),
    ("apical cell membrane", "plasma_membrane", AMBIGUOUS_MEMBRANE),
    ("basolateral cell membrane", "plasma_membrane", AMBIGUOUS_MEMBRANE),
    ("secreted", "secreted", 0.9),
    ("extracellular", "extracellular", 0.85),
    ("nucleus", "nuclear", 0.9),
    ("nucleoplasm", "nuclear", 0.85),
    ("nucleolus", "nuclear", 0.85),
    ("mitochond", "mitochondrial", 0.9),
    ("cytoplasm", "cytoplasmic", 0.85),
    ("cytosol", "cytoplasmic", 0.85),
    ("endoplasmic reticulum", "other_intracellular", 0.8),
    ("golgi", "other_intracellular", 0.8),
    ("endosome", "other_intracellular", 0.8),
    ("lysosome", "other_intracellular", 0.8),
    ("peroxisome", "other_intracellular", 0.8),
    ("membrane", "transmembrane", 0.5),
)

KEYWORD_TERMS: dict[str, tuple[str, float]] = {
    "cell membrane": ("plasma_membrane", AMBIGUOUS_MEMBRANE),
    "membrane": ("transmembrane", 0.45),
    "transmembrane": ("transmembrane", 0.85),
    "transmembrane helix": ("transmembrane", 0.9),
    "gpi-anchor": ("gpi_anchored", 0.9),
    "secreted": ("secreted", 0.85),
    "nucleus": ("nuclear", 0.85),
    "cytoplasm": ("cytoplasmic", 0.8),
    "mitochondrion": ("mitochondrial", 0.85),
}

HPA_TERMS: dict[str, tuple[str, float]] = {
    "plasma membrane": ("plasma_membrane", AMBIGUOUS_MEMBRANE),
    "nucleoplasm": ("nuclear", 0.7),
    "nucleoli": ("nuclear", 0.7),
    "nuclear speckles": ("nuclear", 0.7),
    "nuclear bodies": ("nuclear", 0.7),
    "cytosol": ("cytoplasmic", 0.7),
    "mitochondria": ("mitochondrial", 0.7),
    "golgi apparatus": ("other_intracellular", 0.65),
    "endoplasmic reticulum": ("other_intracellular", 0.65),
    "vesicles": ("other_intracellular", 0.6),
}

GO_TERMS: dict[str, tuple[str, float, str]] = {
    "GO:0005886": ("plasma_membrane", 0.55, "plasma membrane"),
    "GO:0009986": ("cell_surface", 0.7, "cell surface"),
    "GO:0005576": ("extracellular", 0.7, "extracellular region"),
    "GO:0005634": ("nuclear", 0.7, "nucleus"),
    "GO:0005737": ("cytoplasmic", 0.65, "cytoplasm"),
    "GO:0005739": ("mitochondrial", 0.7, "mitochondrion"),
    "GO:0031225": ("gpi_anchored", 0.7, "anchored component of membrane"),
}

OUTSIDE = ("extracellular", "lumenal", "lumenal, vesicle", "periplasmic", "exoplasmic")
INSIDE = ("cytoplasmic", "nuclear", "mitochondrial intermembrane", "mitochondrial matrix")


def _features(defn: dict[str, Any], section: str) -> list[dict[str, Any]]:
    items = (defn.get("sections", {}).get(section) or {}).get("items")
    if isinstance(items, dict):
        return items.get("features") or []
    return items or []


def localise(
    gene: str, defn: dict[str, Any], go_terms: frozenset[str] | set[str] | None = None
) -> Localisation:
    """Build the localisation of one protein from its compiled definition."""
    loc = Localisation()
    sections = defn.get("sections", {})
    ident = (sections.get("identity") or {}).get("items") or {}
    accession = ident.get("accession")
    fn = (sections.get("function") or {}).get("items") or {}

    external: list[str] = []  # positive evidence that part of the protein faces outwards

    def vote(compartment: str, confidence: float) -> None:
        loc.compartments[compartment] = max(loc.compartments.get(compartment, 0.0), confidence)

    # 1. curated per-residue topology: the strongest statement there is
    topo = [f for f in _features(defn, "domains") if f.get("type") == "Topological domain"]
    tm = [f for f in _features(defn, "domains") if f.get("type") == "Transmembrane"]
    for f in topo:
        desc = (f.get("description") or "").lower()
        region = Region("topological_domain", f.get("start"), f.get("end"), f.get("description") or "")
        if any(desc.startswith(o) for o in OUTSIDE):
            loc.extracellular_regions.append(region)
        elif any(desc.startswith(i) for i in INSIDE):
            loc.intracellular_regions.append(region)
    for f in tm:
        loc.transmembrane_regions.append(
            Region("transmembrane", f.get("start"), f.get("end"), f.get("description") or "")
        )
    if loc.transmembrane_regions:
        vote("transmembrane", 0.95)
        n = len(loc.transmembrane_regions)
        loc.topology = "single-pass" if n == 1 else f"multi-pass ({n} segments)"
        loc.evidence.append(
            database(
                "UniProtKB/Swiss-Prot topology",
                f"{gene} has {n} curated transmembrane segment{'s' if n != 1 else ''}"
                + (
                    f" and {len(loc.extracellular_regions)} extracellular topological domain(s)"
                    if loc.extracellular_regions
                    else ""
                ),
                0.95,
                "human",
                identifier=accession,
            )
        )
    if loc.extracellular_regions:
        external.append("curated extracellular topological domain")
        vote("plasma_membrane", 0.9 if loc.transmembrane_regions else 0.6)
        longest = max(loc.extracellular_regions, key=lambda r: (r.end or 0) - (r.start or 0))
        span = (longest.end or 0) - (longest.start or 0) + 1
        loc.evidence.append(
            database(
                "UniProtKB/Swiss-Prot topology",
                f"{gene} exposes residues {longest.start}-{longest.end} ({span} aa) outside the cell",
                0.95,
                "human",
                identifier=accession,
            )
        )
        if loc.transmembrane_regions and loc.intracellular_regions:
            first_tm = min(loc.transmembrane_regions, key=lambda r: r.start or 0)
            outside_first = any((r.start or 0) < (first_tm.start or 0) for r in loc.extracellular_regions)
            loc.orientation = "N-terminus outside" if outside_first else "N-terminus inside"

    # 2. curated subcellular-location text
    for phrase in fn.get("location") or []:
        low = phrase.lower()
        if "cell surface" in low:
            external.append("curated cell-surface annotation")
        for term, compartment, conf in LOCATION_TERMS:
            if term in low:
                vote(compartment, conf)
                break
    if fn.get("location"):
        loc.evidence.append(
            database(
                "UniProtKB/Swiss-Prot subcellular location",
                f"{gene} is annotated at: " + "; ".join(sorted(set(fn["location"]))[:8]),
                0.85,
                "human",
                identifier=accession,
            )
        )

    # 3. keywords, including the GPI anchor and the signal peptide
    for kw in ident.get("keywords") or []:
        hit = KEYWORD_TERMS.get(kw.lower())
        if hit:
            vote(hit[0], hit[1])
    lipids = [
        f
        for f in _features(defn, "modifications")
        if f.get("type") == "Lipidation" and "gpi" in (f.get("description") or "").lower()
    ]
    if lipids:
        external.append("GPI anchor")
        vote("gpi_anchored", 0.9)
        loc.evidence.append(
            database(
                "UniProtKB/Swiss-Prot lipidation",
                f"{gene} carries a GPI anchor ({lipids[0].get('description')}), so it sits on the outer "
                "leaflet with no intracellular segment",
                0.9,
                "human",
                identifier=accession,
            )
        )
    signal = next((f for f in _features(defn, "processing") if f.get("type") == "Signal"), None)
    if signal:
        loc.signal_peptide = Region("signal_peptide", signal.get("start"), signal.get("end"))
        external.append("signal peptide")
        vote("transmembrane", 0.4)
        loc.evidence.append(
            database(
                "UniProtKB/Swiss-Prot processing",
                f"{gene} has a signal peptide ({signal.get('start')}-{signal.get('end')}): it enters the "
                "secretory pathway",
                0.9,
                "human",
                identifier=accession,
            )
        )

    # 4. Human Protein Atlas immunofluorescence
    hpa = (sections.get("expression") or {}).get("items") or {}
    hpa_locs = list(hpa.get("subcellular_main") or []) + list(hpa.get("subcellular_additional") or [])
    for h in hpa_locs:
        hit = HPA_TERMS.get(str(h).lower())
        if hit:
            vote(hit[0], hit[1])
    if hpa_locs:
        loc.evidence.append(
            database(
                "Human Protein Atlas (immunofluorescence)",
                f"{gene} imaged in: " + ", ".join(str(x) for x in hpa_locs[:6]),
                0.7,
                "human",
            )
        )

    # 5. Gene Ontology cellular component, as a second opinion
    if go_terms:
        hits = [(t, GO_TERMS[t]) for t in GO_TERMS if t in go_terms]
        for t, (compartment, conf, _name) in hits:
            if t in ("GO:0009986", "GO:0031225"):
                external.append(f"Gene Ontology {t}")
            vote(compartment, conf)
        if hits:
            loc.evidence.append(
                database(
                    "Gene Ontology cellular component (goa_human)",
                    f"{gene} is annotated to " + ", ".join(f"{n} ({t})" for t, (_c, _f, n) in hits),
                    0.7,
                    "human",
                )
            )

    if not loc.compartments:
        loc.compartments["unknown"] = 0.0
        loc.evidence.append(
            derived(
                "GenomeOS localisation",
                f"no curated localisation found for {gene}; compartment unknown",
                0.0,
            )
        )
    _resolve_orientation(gene, loc, defn, external)
    loc.plasma_membrane = _plasma_membrane(loc)
    loc.confidence = max(loc.compartments.values()) if loc.compartments else 0.0
    return loc


def _resolve_orientation(gene: str, loc: Localisation, defn: dict[str, Any], external: list[str]) -> None:
    """Decide whether anything about this protein actually faces outwards.

    Three cases have to be told apart, and a database that says "Cell membrane"
    tells them apart for none of them:

      * a receptor with an ectodomain, which a circulating binder can engage;
      * a protein tethered to the inner leaflet by a lipid anchor, which it
        cannot, however membrane-bound the annotation says it is;
      * a peripheral or junctional protein whose side is simply not stated.
    """
    if external:
        if loc.compartments.get("plasma_membrane", 0.0) < 0.85:
            loc.compartments["plasma_membrane"] = max(loc.compartments.get("plasma_membrane", 0.0), 0.85)
        loc.orientation = loc.orientation if loc.orientation != "unknown" else "ectodomain outside"
        loc.evidence.append(
            derived(
                "GenomeOS accessibility",
                f"{gene} has positive evidence of an outward-facing part: "
                + ", ".join(sorted(set(external))),
                0.85,
            )
        )
        return

    anchors = [
        f.get("description", "")
        for f in _features(defn, "modifications")
        if f.get("type") == "Lipidation"
        and any(a in (f.get("description") or "").lower() for a in INNER_LEAFLET_ANCHORS)
    ]
    membrane = loc.compartments.get("plasma_membrane", 0.0) > 0.0
    if membrane and anchors and not loc.transmembrane_regions:
        loc.compartments["plasma_membrane"] = min(loc.compartments["plasma_membrane"], 0.25)
        loc.compartments["cytoplasmic"] = max(loc.compartments.get("cytoplasmic", 0.0), 0.85)
        loc.orientation = "cytoplasmic face (lipid anchor)"
        loc.evidence.append(
            derived(
                "GenomeOS accessibility",
                f"{gene} is held against the membrane by a lipid anchor ({anchors[0]}) with no "
                "transmembrane segment and no extracellular domain: it sits on the cytoplasmic face "
                "and a circulating binder cannot reach it, despite the 'Cell membrane' annotation",
                0.85,
            )
        )
        return
    if membrane and not loc.transmembrane_regions:
        loc.compartments["plasma_membrane"] = min(loc.compartments["plasma_membrane"], AMBIGUOUS_MEMBRANE)
        loc.orientation = "membrane side not established"
        loc.evidence.append(
            derived(
                "GenomeOS accessibility",
                f"{gene} is annotated at the membrane but with no transmembrane segment, no "
                "extracellular topological domain, no signal peptide and no cell-surface annotation: "
                "which side of the membrane it occupies is not established, so external accessibility "
                "cannot be claimed",
                0.0,
            )
        )


def _plasma_membrane(loc: Localisation) -> bool | None:
    """True only with positive evidence of an outward-facing part."""
    surface = max(
        (loc.compartments.get(c, 0.0) for c in ("plasma_membrane", "cell_surface", "gpi_anchored")),
        default=0.0,
    )
    internal = max(
        (
            loc.compartments.get(c, 0.0)
            for c in ("nuclear", "cytoplasmic", "mitochondrial", "other_intracellular")
        ),
        default=0.0,
    )
    if surface >= REACHABLE_THRESHOLD:
        return True
    if internal >= REACHABLE_THRESHOLD and surface < AMBIGUOUS_MEMBRANE:
        return False
    return None


def mutation_topology(
    loc: Localisation, residue: int | None, gene: str, change: str
) -> tuple[bool | None, Evidence | None]:
    """Is the altered residue itself outside the cell?

    This is the difference between a mutation-specific extracellular epitope
    and a change a circulating binder can never see.
    """
    if residue is None:
        return None, None
    outside = loc.extracellular_residue(residue)
    if outside is None:
        return None, None
    where = next(
        (
            r
            for r in (loc.extracellular_regions + loc.intracellular_regions + loc.transmembrane_regions)
            if r.contains(residue)
        ),
        None,
    )
    label = (where.description or where.kind) if where else ("extracellular" if outside else "internal")
    claim = (
        f"{gene} {change} falls in the {label.lower()} segment {where.start}-{where.end}"
        if where
        else f"{gene} {change} is {label}"
    )
    return outside, database("UniProtKB/Swiss-Prot topology", claim, 0.95, "human")
