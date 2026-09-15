"""The `cancer.*` library layer: which genes break, how, and in which cancers.

The BioLib catalogue names the toolkits a healthy genome uses. This layer
names the ones a tumour breaks, and it is the only layer whose membership is
a *frequency* rather than a function: a gene belongs to `cancer.amplified`
because tumours amplify it, not because of what its product does.

Membership is computed from the two committed distillations
(`cancer_msk_impact_2017` for mutations, `cancer_alterations_msk_impact_2017`
for copy number and structural variants), never hand-written, so a library
cannot drift from the evidence and every member carries the number that put
it there. Nothing here needs the network.

Three limits are part of the answer and travel with it. The distillation
covers a 110-gene driver panel of a 468-gene sequencing study, so a gene
absent from a library may simply never have been looked at: silence here is
not evidence of a healthy gene, and `breaks_in` says which of the two it is.
A frequency over all 10,945 tumours hides the cancer type where the gene
actually breaks — ERBB2 amplification is 4.0% overall and 14.1% of breast
carcinomas — so the type-level numbers are reported beside the overall one,
except for hotspots, where the distillation counts changes study-wide and no
per-type number exists to report. And a recurrent alteration is evidence that
a gene is *selected*, which is not evidence that it drives, still less that
it can be treated.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from genomeos.results import load_result

LAYER = "cancer"
MUTATIONS = "cancer_msk_impact_2017"
ALTERATIONS = "cancer_alterations_msk_impact_2017"

#: A gene joins a frequency-defined library at or above this share of tumours.
#: 1% of 10,945 tumours is about 110 samples, which is where a per-type
#: frequency stops being a handful of cases.
THRESHOLD = 0.01

#: A fusion is rarer than a mutation by an order of magnitude, so it gets its
#: own, lower bar; EML4-ALK is in 0.4% of all tumours and defines a disease.
FUSION_THRESHOLD = 0.001


@dataclass(frozen=True, slots=True)
class CancerLibrary:
    """One way a gene breaks, and the genes that break that way."""

    id: str
    purpose: str
    kind: str  # mutation | amplification | deep_deletion | fusion | hotspot
    threshold: float
    layer: str = LAYER
    note: str = ""

    @property
    def source(self) -> str:
        return f"cBioPortal {MUTATIONS if self.kind in ('mutation', 'hotspot') else ALTERATIONS}"


CANCER_LIBRARIES: dict[str, CancerLibrary] = {}


def _add(lib: CancerLibrary) -> None:
    CANCER_LIBRARIES[lib.id] = lib


_add(
    CancerLibrary(
        "cancer.mutated",
        "genes a tumour breaks by point mutation",
        "mutation",
        THRESHOLD,
        note="a mutation frequency alone under-counts: CCND1 and MYC are here at under 1% and "
        "amplified in over 4%",
    )
)
_add(
    CancerLibrary(
        "cancer.hotspots",
        "genes broken at the same residue again and again",
        "hotspot",
        THRESHOLD,
        note="membership is the share of tumours carrying the gene's single commonest protein change; "
        "recurrence at one residue is the strongest single argument that a change is selected rather "
        "than carried, and it says nothing about whether the residue is reachable. The distillation "
        "counts recurrent changes study-wide and not per cancer type, so no per-type breakdown of a "
        "hotspot exists here",
    )
)
_add(
    CancerLibrary(
        "cancer.amplified",
        "oncogenes a tumour multiplies",
        "amplification",
        THRESHOLD,
        note="amplification raises how much product a cell could display and never shows that it does",
    )
)
_add(
    CancerLibrary(
        "cancer.deleted",
        "tumour suppressors a tumour removes outright",
        "deep_deletion",
        THRESHOLD,
        note="a deleted gene makes no product, so it is never a binder target; what it is worth is the "
        "dependency the loss creates",
    )
)
_add(
    CancerLibrary(
        "cancer.rearranged",
        "genes a tumour joins to another gene",
        "fusion",
        FUSION_THRESHOLD,
        note="a fusion protein exists in no healthy cell, which is the cleanest tumour specificity "
        "there is, and the junction sequence is not reconstructed here",
    )
)


# ---- membership, computed from the distilled tables ---------------------------------


@dataclass(slots=True)
class Knowledge:
    """The two committed tables, read once."""

    mutations: dict[str, Any] = field(default_factory=dict)
    alterations: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def load(cls, mutations: dict | None = None, alterations: dict | None = None) -> Knowledge:
        return cls(
            mutations if mutations is not None else (load_result(MUTATIONS) or {}),
            alterations if alterations is not None else (load_result(ALTERATIONS) or {}),
        )

    @property
    def available(self) -> bool:
        return bool(self.mutations or self.alterations)

    @property
    def samples(self) -> int:
        return int(self.mutations.get("samples") or self.alterations.get("samples") or 0)

    @property
    def study(self) -> str:
        return str(self.mutations.get("study") or self.alterations.get("study") or "")

    def panel(self) -> set[str]:
        return set(self.mutations.get("genes") or {}) | set(self.alterations.get("genes") or {})

    def frequency(self, gene: str, kind: str) -> float | None:
        """How often tumours carry this kind of break in this gene, or None if not called."""
        if kind == "mutation":
            entry = (self.mutations.get("genes") or {}).get(gene)
            return entry["frequency"] if entry else None
        if kind == "hotspot":
            entry = (self.mutations.get("genes") or {}).get(gene)
            if not entry or not entry.get("hotspots"):
                return None
            top = entry["hotspots"][0]
            # the share of all tumours carrying the single commonest change
            return top[1] / self.samples if self.samples else None
        entry = ((self.alterations.get("genes") or {}).get(gene) or {}).get(kind)
        return entry["frequency"] if entry else None

    def by_cancer_type(self, gene: str, kind: str) -> dict[str, float]:
        """Per-type frequency of this break, where the distillation records one.

        Deliberately empty for a hotspot. The mutation table counts recurrent
        protein changes over the whole study and not per cancer type, so the
        only per-type numbers available are the gene's *mutation* frequencies.
        Reporting those beside a hotspot frequency would read as the hotspot's
        distribution, which nothing here measured.
        """
        if kind == "hotspot":
            return {}
        if kind == "mutation":
            entry = (self.mutations.get("genes") or {}).get(gene) or {}
            return dict(sorted((entry.get("by_cancer_type") or {}).items(), key=lambda kv: -kv[1])[:6])
        entry = ((self.alterations.get("genes") or {}).get(gene) or {}).get(kind) or {}
        return dict(entry.get("by_cancer_type") or {})

    def hotspots(self, gene: str) -> list[tuple[str, int]]:
        entry = (self.mutations.get("genes") or {}).get(gene) or {}
        return [tuple(h) for h in entry.get("hotspots", []) if h[0] not in ("?", "")]

    def partners(self, gene: str) -> list[tuple[str, int]]:
        entry = (self.alterations.get("genes") or {}).get(gene) or {}
        return [tuple(p) for p in entry.get("fusion_partners", [])]


def library_members(lib: CancerLibrary, k: Knowledge | None = None) -> dict[str, float]:
    """The genes in this library, each with the frequency that put it there."""
    k = k or Knowledge.load()
    out = {}
    for gene in sorted(k.panel()):
        f = k.frequency(gene, lib.kind)
        if f is not None and f >= lib.threshold:
            out[gene] = round(f, 5)
    return dict(sorted(out.items(), key=lambda kv: -kv[1]))


def cancer_libraries_of(gene: str, k: Knowledge | None = None) -> list[str]:
    """Which cancer libraries a gene belongs to."""
    k = k or Knowledge.load()
    gene = gene.upper()
    return [lib.id for lib in CANCER_LIBRARIES.values() if gene in library_members(lib, k)]


def breaks_in(gene: str, k: Knowledge | None = None, top: int = 5) -> dict[str, Any]:
    """How this gene breaks, in which cancers, with the evidence for each claim.

    The question the layer exists to answer. Every way the study records the
    gene breaking is listed with its overall frequency, the cancer types that
    carry it, and what the alteration means for a therapy; a way the study
    does not record is reported as not called rather than as zero.
    """
    k = k or Knowledge.load()
    gene = gene.upper()
    if not k.available:
        return {"gene": gene, "available": False, "reason": "no distilled cancer knowledge"}
    on_panel = gene in k.panel()
    ways = []
    for lib in CANCER_LIBRARIES.values():
        f = k.frequency(gene, lib.kind)
        if f is None or f <= 0:
            continue
        types = k.by_cancer_type(gene, lib.kind)
        way: dict[str, Any] = {
            "library": lib.id,
            "kind": lib.kind,
            "frequency": round(f, 5),
            "member": f >= lib.threshold,
            "by_cancer_type": dict(list(types.items())[:top]),
            "meaning": lib.note,
            "evidence": {
                "kind": "experimental",
                "source": f"{lib.source} ({k.samples} tumours)",
                "claim": (
                    f"{gene} carries this break in {f:.2%} of tumours"
                    + (
                        f", most often in {next(iter(types))} ({next(iter(types.values())):.1%})"
                        if types
                        else ""
                    )
                ),
            },
            "confidence": 0.8,
        }
        if lib.kind == "hotspot":
            way["recurrent_changes"] = k.hotspots(gene)[:5]
        if lib.kind == "fusion":
            way["partners"] = k.partners(gene)[:5]
        ways.append(way)
    ways.sort(key=lambda w: -w["frequency"])
    return {
        "gene": gene,
        "available": True,
        "on_panel": on_panel,
        "study": k.study,
        "samples": k.samples,
        "libraries": [w["library"] for w in ways if w["member"]],
        "breaks": ways,
        "not_called": [lib.kind for lib in CANCER_LIBRARIES.values() if k.frequency(gene, lib.kind) is None],
        "note": (
            f"{gene} is on the {k.study} panel and every absence below is an absence of a call"
            if on_panel
            else f"{gene} is not on the {k.study} panel: this study cannot say whether it breaks, and "
            "the silence is not evidence that it does not"
        ),
    }


def summary(k: Knowledge | None = None) -> dict[str, Any]:
    """The whole layer: each library, its size, and its largest members."""
    k = k or Knowledge.load()
    return {
        "layer": LAYER,
        "study": k.study,
        "samples": k.samples,
        "panel": len(k.panel()),
        "libraries": {
            lib.id: {
                "purpose": lib.purpose,
                "kind": lib.kind,
                "threshold": lib.threshold,
                "source": lib.source,
                "note": lib.note,
                "members": len(m := library_members(lib, k)),
                "top": dict(list(m.items())[:10]),
            }
            for lib in CANCER_LIBRARIES.values()
        },
        "limits": [
            "a panel study calls only the genes on its panel; a gene absent from a library may be "
            "absent from the panel, and silence is not evidence of a healthy gene",
            "a frequency over every tumour hides the cancer type where the gene actually breaks; "
            "breaks_in() reports the per-type numbers",
            "a recurrent alteration shows a gene is selected, not that it drives, and not that it "
            "can be treated",
        ],
    }


def register(libraries: dict[str, Any], layers: tuple[str, ...]) -> tuple[str, ...]:
    """Fold this layer into the BioLib catalogue, for whoever owns it.

    Kept explicit rather than run on import: `genomeos.lib` is another area's
    file, and a layer that appears in its catalogue as a side effect of
    importing `genomeos.cancer` would be a surprise. One call wires it.
    """
    from genomeos.lib.catalog import Library

    k = Knowledge.load()
    for lib in CANCER_LIBRARIES.values():
        genes = tuple(library_members(lib, k))
        libraries[lib.id] = Library(
            id=lib.id,
            layer=LAYER,
            purpose=lib.purpose,
            genes=genes,
            scale=f"{len(genes)} genes of the {len(k.panel())}-gene panel at >= {lib.threshold:g}",
            source=lib.source,
            note=lib.note,
        )
    return layers if LAYER in layers else (*layers, LAYER)
