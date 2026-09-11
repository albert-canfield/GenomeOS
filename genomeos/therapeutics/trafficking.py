"""What happens after the binder attaches.

Binding is the beginning, not the end:

    binding -> internalisation? -> endosome? -> recycling? -> lysosome?

The answer decides the mechanism. Payload delivery needs the receptor to be
swallowed and routed somewhere the cargo can be released; recruiting an NK
cell needs the opposite, a target that stays on the surface long enough to be
seen. A receptor that sheds its ectodomain floods the circulation with decoy
antigen and blunts everything.

None of this is guessed. Curated UniProt keywords, subcellular locations and
Gene Ontology terms are read; where the annotation is silent the field stays
None with the reason attached. Absence of a keyword is not evidence of
absence, and internalisation *rate* is essentially never in a public database,
so it is reported unavailable rather than approximated.
"""

from __future__ import annotations

from typing import Any

from .evidence import database, derived
from .model import Localisation, Trafficking
from .providers import Answer

ENDOCYTOSIS_KEYWORDS = ("endocytosis", "coated pit")
ENDOCYTOSIS_GO = {
    "GO:0006898": "receptor-mediated endocytosis",
    "GO:0031623": "receptor internalization",
    "GO:0072583": "clathrin-dependent endocytosis",
}
ENDOSOME_GO = {"GO:0005768": "endosome"}
LYSOSOME_GO = {"GO:0005764": "lysosome"}
RECYCLING_GO = {"GO:0055037": "recycling endosome", "GO:0001881": "receptor recycling"}

SHEDDING_HINTS = ("sheddase", "ectodomain", "shedding", "soluble form", "released")


def assess(gene: str, answer: Answer, localisation: Localisation) -> Trafficking:
    """Read trafficking properties out of curated annotation."""
    t = Trafficking()
    if not answer.available:
        t.reason = answer.reason
        return t
    data: dict[str, Any] = answer.data
    keywords = set(data.get("keywords") or [])
    locations = set(data.get("locations") or [])
    go = set(data.get("go_terms") or [])
    found: list[str] = []

    kw_hits = sorted(k for k in keywords if k in ENDOCYTOSIS_KEYWORDS)
    go_hits = sorted(g for g in ENDOCYTOSIS_GO if g in go)
    if kw_hits or go_hits:
        t.internalises = True
        parts = [f"UniProt keyword '{k}'" for k in kw_hits]
        parts += [f"GO {g} ({ENDOCYTOSIS_GO[g]})" for g in go_hits]
        found.append("internalisation")
        t.evidence.append(
            database(
                "UniProtKB and Gene Ontology",
                f"{gene} is annotated for receptor internalisation: " + "; ".join(parts),
                0.85,
                "human",
            )
        )

    if any("endosome" in x for x in locations) or (ENDOSOME_GO.keys() & go):
        t.endosomal = True
        found.append("endosomal routing")
        t.evidence.append(
            database(
                "UniProtKB subcellular location / Gene Ontology",
                f"{gene} is found in endosomes after uptake",
                0.8,
                "human",
            )
        )
    if any("lysosome" in x for x in locations) or (LYSOSOME_GO.keys() & go):
        t.lysosomal = True
        found.append("lysosomal routing")
        t.evidence.append(
            database(
                "UniProtKB subcellular location / Gene Ontology",
                f"{gene} reaches the lysosome, the compartment where an acid- or protease-cleaved "
                "payload would be released",
                0.8,
                "human",
            )
        )
    if any("recycling" in x for x in locations) or (RECYCLING_GO.keys() & go):
        t.recycling = True
        found.append("recycling")
        t.evidence.append(
            database(
                "UniProtKB subcellular location / Gene Ontology",
                f"{gene} recycles back to the surface; recycling competes with lysosomal delivery but "
                "raises the number of payload rounds per receptor",
                0.75,
                "human",
            )
        )

    membrane = localisation.reachable
    secreted = localisation.compartments.get("secreted", 0.0) >= 0.6
    chains = [f for f in data.get("processing") or [] if f.get("type") in ("Chain", "Peptide")]
    text = " ".join(str(f.get("description", "")).lower() for f in chains)
    if membrane and secreted:
        t.shedding = True
        t.soluble_antigen = True
        found.append("shedding")
        t.evidence.append(
            database(
                "UniProtKB subcellular location",
                f"{gene} is annotated both at the membrane and as secreted: an ectodomain or cleaved "
                "fragment circulates, which can absorb a binder before it reaches the cell",
                0.75,
                "human",
            )
        )
    elif membrane and (len(chains) > 1 or any(h in text for h in SHEDDING_HINTS)):
        t.shedding = True
        t.soluble_antigen = None
        found.append("proteolytic processing")
        t.evidence.append(
            database(
                "UniProtKB processing features",
                f"{gene} is proteolytically processed into {len(chains)} chains; ectodomain release is "
                "possible but not stated",
                0.55,
                "human",
            )
        )

    if membrane:
        t.membrane_residence = "surface-resident" if t.internalises is None else "internalising"

    if found:
        t.evidence.append(
            derived(
                "GenomeOS trafficking",
                f"{gene} trafficking established for: " + ", ".join(found),
                0.6,
            )
        )
    else:
        t.reason = (
            "no curated endocytosis, endosomal, lysosomal or recycling annotation for this protein; "
            "internalisation is unknown, not absent"
        )
        t.evidence.append(derived("GenomeOS trafficking", f"{gene}: {t.reason}", 0.0))
    if not membrane and t.internalises is None:
        t.reason = t.reason or (
            "protein is not at the cell surface, so post-binding trafficking does not apply to a "
            "circulating binder"
        )
    return t


def internalisation_score(t: Trafficking) -> tuple[float | None, str]:
    """0-1 for payload-delivery suitability, or None when nothing is known."""
    if t.internalises is None:
        return None, t.reason or "no internalisation evidence found"
    if not t.internalises:
        return 0.0, "curated annotation states the receptor is not internalised"
    score = 0.6
    parts = ["curated internalisation annotation (0.60)"]
    if t.endosomal:
        score += 0.1
        parts.append("endosomal routing (+0.10)")
    if t.lysosomal:
        score += 0.2
        parts.append("lysosomal routing (+0.20)")
    if t.recycling:
        score -= 0.1
        parts.append("recycling competes with degradative routing (-0.10)")
    if t.shedding:
        score -= 0.15
        parts.append("shedding produces soluble decoy antigen (-0.15)")
    return max(0.0, min(1.0, score)), " ".join(parts) + (
        "; kinetics unavailable, so this is a qualitative ranking, not a rate"
    )


def lysosomal_score(t: Trafficking) -> tuple[float | None, str]:
    """Separate from internalisation: reaching the lysosome is its own question."""
    if t.lysosomal is True:
        return 1.0, "curated lysosomal localisation"
    if t.endosomal is True:
        return 0.5, "endosomal but no curated lysosomal annotation"
    if t.internalises is None:
        return None, t.reason or "no trafficking annotation"
    return None, "internalisation annotated but the downstream compartment is not"
