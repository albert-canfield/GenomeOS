# SPDX-License-Identifier: AGPL-3.0-or-later
"""The alterations that are not point mutations: copy number and structure.

A tumour breaks a gene in three ways and a VCF of point mutations carries
only one of them. The other two are among the most actionable findings in a
real tumour: a deep deletion removes a tumour suppressor outright, and an
amplification multiplies an oncogene until the cell displays far more of its
product than any healthy cell does. ERBB2 is the case that settles the
argument — mutated in 3% of MSK-IMPACT tumours, amplified in 14% of the
breast carcinomas, and it is the amplification that trastuzumab was approved
against.

This module gives those events the same standing a variant has: a typed
record, a frequency from the distilled cohort, a score on the same scale as
the variant score, and an evidence line per claim. Two things it refuses to
do. It never converts a copy-number call into an amount of protein: copy
number bounds what a cell could make and never shows that it does. And it
never reads a deep deletion as a reason to aim a binder at the gene — a
homozygously deleted gene makes no product, so the actionable direction is
the opposite one, and `removes_product` is how the rest of the pipeline is
told.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from genomeos.results import load_result

#: The alteration kinds, ordered from the strongest gain to the strongest loss.
KINDS: tuple[str, ...] = ("amplification", "gain", "fusion", "shallow_deletion", "deep_deletion")

#: Absolute copy-number thresholds. A focal amplification in a real tumour
#: runs to dozens of copies; 4 is the usual line between a gain and an
#: amplification, and 0.5 between one lost allele and both.
AMPLIFIED_COPIES = 4.0
GAINED_COPIES = 2.5
SHALLOW_LOSS_COPIES = 1.5
DEEP_LOSS_COPIES = 0.5

#: GISTIC discrete calls, as cBioPortal stores them.
GISTIC: dict[int, str] = {2: "amplification", 1: "gain", -1: "shallow_deletion", -2: "deep_deletion"}

#: Base score per kind, on the same scale as the variant severity in
#: genomeos.cancer.tumour: 1.0 is what a truncating mutation scores.
KIND_SCORE: dict[str, float] = {
    "deep_deletion": 1.0,
    "amplification": 1.0,
    "fusion": 0.9,
    "gain": 0.3,
    "shallow_deletion": 0.2,
}

#: A fusion partner seen this often in the cohort is recurrent, which is to a
#: fusion what a hotspot residue is to a missense change.
RECURRENT_PARTNER = 3

KNOWLEDGE = "cancer_alterations_msk_impact_2017"


@dataclass(slots=True)
class GeneAlteration:
    """One copy-number or structural event on one gene in one tumour."""

    gene: str
    kind: str
    copies: float | None = None  # absolute copies, when the source gives them
    gistic: int | None = None  # discrete call, when the source gives one
    partner: str = ""  # the other side of a fusion
    detail: str = ""  # DELETION, INVERSION, TRANSLOCATION, ...
    in_frame: str = ""
    #: Which end of the fusion product this gene contributes: "5'" or "3'",
    #: empty when the source does not say. It decides whether the gene's own
    #: N-terminus — its signal peptide, and for a type-I receptor its whole
    #: ectodomain — is in the product at all, which is the difference between
    #: a surface target and a cytoplasmic one.
    orientation: str = ""
    #: The junction as the source describes it, e.g. "EML4 exons 1-20 with
    #: ALK exons 20-29", kept verbatim so the claim can be checked.
    junction: str = ""
    source: str = ""
    cohort_frequency: float | None = None
    cohort_types: dict[str, float] = field(default_factory=dict)
    cohort_study: str = ""
    recurrent_partner: bool = False
    score: float = 0.0
    evidence: list[str] = field(default_factory=list)

    @property
    def removes_product(self) -> bool:
        """Both copies gone: there is no protein left for a binder to find."""
        return self.kind == "deep_deletion"

    @property
    def raises_product(self) -> bool:
        return self.kind in ("amplification", "gain")

    def label(self) -> str:
        if self.kind == "fusion":
            return f"{self.gene}-{self.partner} fusion" if self.partner else f"{self.gene} rearrangement"
        if self.copies is not None:
            return f"{self.gene} {self.kind.replace('_', ' ')} ({self.copies:g} copies)"
        return f"{self.gene} {self.kind.replace('_', ' ')}"

    def to_dict(self) -> dict[str, Any]:
        return {
            "gene": self.gene,
            "kind": self.kind,
            "label": self.label(),
            "copies": self.copies,
            "gistic": self.gistic,
            "partner": self.partner or None,
            "detail": self.detail or None,
            "in_frame": self.in_frame or None,
            "source": self.source,
            "cohort_frequency": self.cohort_frequency,
            "cohort_by_cancer_type": self.cohort_types,
            "cohort_study": self.cohort_study,
            "recurrent_partner": self.recurrent_partner,
            "removes_product": self.removes_product,
            "score": round(self.score, 3),
            "evidence": self.evidence,
        }


# ---- reading a patient's table ----------------------------------------------------


def detect_cna_format(values: list[float]) -> tuple[str, str]:
    """GISTIC discrete calls or absolute copies? The two disagree about 2.

    A GISTIC 2 is an amplification; 2 absolute copies is a normal diploid gene.
    Reading one table as the other would invent an amplification on every
    untouched gene, so the rule is stated rather than guessed at: a negative
    value can only be a discrete call, and anything above 2 can only be a copy
    count. A table whose values are all 0, 1 or 2 is genuinely ambiguous and is
    read as copies, which is the reading that invents no amplification.
    """
    if not values:
        return "copies", "the table is empty"
    if any(v < 0 for v in values):
        return "gistic", "negative values are present, which only a discrete GISTIC call can be"
    if any(v > 2 for v in values):
        return "copies", "values above 2 are present, which only an absolute copy count can be"
    if any(v != int(v) for v in values):
        return "copies", "fractional values are present, which a discrete call cannot be"
    return "copies", (
        "every value is 0, 1 or 2, which both conventions produce; read as absolute copies because that "
        "is the reading that invents no amplification. Pass an explicit format if this table is GISTIC"
    )


def kind_from_copies(copies: float) -> str | None:
    if copies >= AMPLIFIED_COPIES:
        return "amplification"
    if copies > GAINED_COPIES:
        return "gain"
    if copies <= DEEP_LOSS_COPIES:
        return "deep_deletion"
    if copies < SHALLOW_LOSS_COPIES:
        return "shallow_deletion"
    return None


def rows(path: Path) -> list[list[str]]:
    """Non-comment, non-blank lines of a small table, split on tabs or commas."""
    out = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = [p.strip() for p in line.replace(",", "\t").split("\t") if p.strip()]
        if parts:
            out.append(parts)
    return out


def read_cna_table(path: str, fmt: str = "auto") -> list[GeneAlteration]:
    """Gene and copy-number call, from a two-column table.

    `fmt` is "auto", "gistic" or "copies"; auto applies the rule in
    `detect_cna_format` and writes the reading it chose, with its reason, into
    the evidence of every alteration it produces.
    """
    p = Path(path)
    if not p.exists():
        return []
    pairs: list[tuple[str, float]] = []
    for parts in rows(p):
        if len(parts) < 2:
            continue
        try:
            pairs.append((parts[0].upper(), float(parts[1])))
        except ValueError:
            continue  # a header row
    if fmt == "auto":
        fmt, reason = detect_cna_format([v for _, v in pairs])
    else:
        reason = f"format given as {fmt}"
    out: list[GeneAlteration] = []
    for gene, value in pairs:
        if fmt == "gistic":
            kind = GISTIC.get(int(value))
            a = GeneAlteration(gene, kind or "", gistic=int(value)) if kind else None
        else:
            kind = kind_from_copies(value)
            a = GeneAlteration(gene, kind, copies=value) if kind else None
        if a is None:
            continue
        a.source = f"patient copy-number table ({p.name}), read as {fmt}"
        a.evidence.append(f"{a.label()} in this tumour's copy-number table; {reason}")
        out.append(a)
    return out


#: The column names a structural-variant table may carry, normalised to
#: lower-case alphanumerics, against the field each one fills. Two naming
#: families are here on purpose. `gene`/`partner` is this project's own
#: four-column form, where the gene order is a naming convention; `site1`/`site2`
#: is how a caller and a cBioPortal export name the two ends of the event, and
#: there the order IS the measurement — site1 is the 5' partner and site2 the 3'
#: (genomeos/cancer/cbioportal.py says the same of the API's rows).
SV_COLUMNS: dict[str, str] = {
    "gene": "gene",
    "gene1": "gene",
    "site1": "site1",
    "site1gene": "site1",
    "site1hugosymbol": "site1",
    "partner": "partner",
    "gene2": "partner",
    "site2": "site2",
    "site2gene": "site2",
    "site2hugosymbol": "site2",
    "class": "detail",
    "svclass": "detail",
    "variantclass": "detail",
    "eventinfo": "event",
    "end": "orientation",
    "orientation": "orientation",
    "inframe": "in_frame",
    "site2effectonframe": "in_frame",
    "annotation": "annotation",
    "site1chromosome": "site1chr",
    "site1position": "site1pos",
    "site2chromosome": "site2chr",
    "site2position": "site2pos",
    "ncbibuild": "build",
}

#: What a caller writes where it has no value. Treated as absent, never as text.
SV_ABSENT = frozenset({"", "na", "n/a", "nan", "none", "null", "."})


def _sv_field(name: str) -> str:
    return SV_COLUMNS.get("".join(ch for ch in name.lower() if ch.isalnum()), "")


def _sv_header(cells: list[str]) -> dict[str, int] | None:
    """The row read as a header, or None when it is not one.

    A header is a row that names at least one gene column and at least one
    further column this reader understands. Anything else is data, and the
    positional reading applies — which is what keeps a file with no header
    working, and what keeps a one-word header from being read as a table.
    """
    found: dict[str, int] = {}
    for i, cell in enumerate(cells):
        field = _sv_field(cell)
        if field and field not in found:
            found[field] = i
    if not ({"gene", "site1"} & set(found)) or len(found) < 2:
        return None
    return found


def _sv_value(cells: list[str], header: dict[str, int], field: str) -> str:
    i = header.get(field)
    if i is None or i >= len(cells):
        return ""
    value = cells[i].strip()
    return "" if value.lower() in SV_ABSENT else value


def _sv_junction(cells: list[str], header: dict[str, int]) -> str:
    """The junction as the source states it: its own annotation, else the breakpoints."""
    annotation = _sv_value(cells, header, "annotation")
    if annotation:
        return annotation
    left = _sv_value(cells, header, "site1chr"), _sv_value(cells, header, "site1pos")
    right = _sv_value(cells, header, "site2chr"), _sv_value(cells, header, "site2pos")
    if not (all(left) and all(right)):
        return ""
    build = _sv_value(cells, header, "build")
    where = f"{left[0]}:{left[1]} to {right[0]}:{right[1]}"
    return f"{where} ({build})" if build else where


def _sv_named(lines: list[str], header: dict[str, int], p: Path) -> list[GeneAlteration]:
    """A table whose header names its columns, one record per gene it names as a site.

    A `site1`/`site2` pair carries the orientation in the column names: site1 is
    the 5' partner and site2 the 3', so each row yields a record for each end
    with the other as its partner. A `gene`/`partner` pair does not — the gene
    order in a fusion's name is a naming convention — so it yields the one
    record the positional reading yields, and its orientation comes only from an
    explicit `end` or `orientation` column.
    """
    out: list[GeneAlteration] = []
    for line in lines:
        cells = line.split("\t") if "\t" in line else line.split(",")
        if not any(c.strip() for c in cells):
            continue
        sides: list[tuple[str, str, str]] = []  # gene, partner, orientation
        if "site1" in header or "site2" in header:
            first = _sv_value(cells, header, "site1").upper()
            second = _sv_value(cells, header, "site2").upper()
            if first:
                sides.append((first, second, "5'"))
            if second:
                sides.append((second, first, "3'"))
        else:
            gene = _sv_value(cells, header, "gene").upper()
            if gene:
                raw = _sv_value(cells, header, "orientation").rstrip("'").rstrip("’")
                sides.append(
                    (gene, _sv_value(cells, header, "partner").upper(), {"5": "5'", "3": "3'"}.get(raw, ""))
                )
        detail = _sv_value(cells, header, "detail") or _sv_value(cells, header, "event")
        junction = _sv_junction(cells, header)
        for gene, partner, orientation in sides:
            a = GeneAlteration(
                gene,
                "fusion",
                partner=partner,
                detail=detail,
                in_frame=_sv_value(cells, header, "in_frame"),
                orientation=orientation,
                junction=junction,
                source=f"patient structural-variant table ({p.name})",
            )
            a.evidence.append(f"{a.label()} in this tumour's structural-variant table")
            out.append(a)
    return out


def read_sv_table(path: str) -> list[GeneAlteration]:
    """Gene, partner, optionally the event class, and optionally the orientation.

    The fourth column is "5" or "3": which end of the fusion product this gene
    contributes. It is optional because most callers do not report it, and a
    fusion whose orientation is unknown must stay unknown rather than be
    guessed — the gene order in a fusion's name is a naming convention, not a
    measurement.

    A table whose first row names its columns is read by those names instead
    (`SV_COLUMNS`), which is what lets a caller's own export be read rather than
    rewritten into four columns by hand. Positionally, a cBioPortal
    structural-variant export reads its sample identifier as the gene and its
    chromosome as the orientation, inventing a candidate named after the sample:
    a record the pipeline already reads correctly from the REST API could not be
    read from the file the same portal exports.
    """
    p = Path(path)
    if not p.exists():
        return []
    lines = [ln for ln in p.read_text().splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        return []
    first = lines[0].split("\t") if "\t" in lines[0] else lines[0].split(",")
    header = _sv_header(first)
    if header is not None:
        return _sv_named(lines[1:], header, p)
    out: list[GeneAlteration] = []
    for parts in rows(p):
        gene = parts[0].upper()
        if gene in ("GENE", "GENE1", "SITE1"):
            continue
        partner = parts[1].upper() if len(parts) > 1 else ""
        orientation = ""
        if len(parts) > 3:
            raw = parts[3].strip().rstrip("'").rstrip("’")
            orientation = {"5": "5'", "3": "3'"}.get(raw, "")
        a = GeneAlteration(
            gene,
            "fusion",
            partner=partner,
            detail=parts[2] if len(parts) > 2 else "",
            orientation=orientation,
            source=f"patient structural-variant table ({p.name})",
        )
        a.evidence.append(f"{a.label()} in this tumour's structural-variant table")
        out.append(a)
    return out


def from_cbioportal(sample: dict[str, Any]) -> list[GeneAlteration]:
    """The alterations of one tumour as `CBioPortal.sample_alterations` returns them."""
    out = []
    for row in sample.get("alterations", []):
        a = GeneAlteration(
            gene=(row.get("gene") or "").upper(),
            kind=row.get("kind", ""),
            gistic=row.get("gistic"),
            partner=row.get("partner", ""),
            detail=row.get("detail", ""),
            in_frame=row.get("in_frame", ""),
            orientation=row.get("orientation", ""),
            junction=row.get("junction", ""),
            source=row.get("source", f"cBioPortal {sample.get('study', '')} {sample.get('sample', '')}"),
        )
        if not a.gene or a.kind not in KIND_SCORE:
            continue
        a.evidence.append(f"{a.label()} called in {sample.get('sample', 'this sample')} ({a.source})")
        out.append(a)
    return out


# ---- grading against the distilled cohort -------------------------------------------


def grade(alterations: list[GeneAlteration], knowledge: dict | None = None) -> list[GeneAlteration]:
    """Attach the cohort frequency of each event and score it.

    Deliberately the same shape as `genomeos.cancer.tumour.grade`: severity of
    the event, plus how often tumours carry it, plus a recurrence bonus. A
    fusion with a partner the cohort keeps producing is to a rearrangement what
    a hotspot residue is to a missense change.
    """
    knowledge = knowledge if knowledge is not None else (load_result(KNOWLEDGE) or {})
    kgenes = knowledge.get("genes", {})
    study = knowledge.get("study", "")
    samples = knowledge.get("samples")
    for a in alterations:
        entry = kgenes.get(a.gene) or {}
        k = entry.get(a.kind) or {}
        if k:
            a.cohort_frequency = k.get("frequency")
            a.cohort_types = dict(list((k.get("by_cancer_type") or {}).items())[:6])
            a.cohort_study = study
            a.evidence.append(
                f"{a.gene} carries this {a.kind.replace('_', ' ')} in {a.cohort_frequency:.2%} of "
                f"{samples} tumours (cBioPortal {study})"
            )
            if a.cohort_types:
                top, freq = next(iter(a.cohort_types.items()))
                a.evidence.append(f"most often in {top} ({freq:.1%} of them)")
        elif entry:
            a.evidence.append(
                f"{a.gene} is on the panel of {study} and this event is not recorded there; the "
                "frequency is unknown rather than zero"
            )
        if a.kind == "fusion" and a.partner:
            partners = dict(entry.get("fusion_partners") or [])
            seen = partners.get(a.partner, 0)
            if seen >= RECURRENT_PARTNER:
                a.recurrent_partner = True
                a.evidence.append(
                    f"{a.gene}-{a.partner} is a recurrent pairing: {seen} tumours of {study} carry it"
                )
        if a.removes_product:
            a.evidence.append(
                "a deep deletion removes both copies, so this gene makes no product: it is a reason to "
                "look for what the tumour now depends on, never a target for a binder"
            )
        a.score = score(a)
    alterations.sort(key=lambda a: -a.score)
    return alterations


def score(a: GeneAlteration) -> float:
    s = KIND_SCORE.get(a.kind, 0.1)
    if a.cohort_frequency:
        s += min(1.0, math.log10(1 + 100 * a.cohort_frequency))
    if a.recurrent_partner:
        s += 1.0
    return s


def by_gene(alterations: list[GeneAlteration]) -> dict[str, list[GeneAlteration]]:
    out: dict[str, list[GeneAlteration]] = {}
    for a in alterations:
        out.setdefault(a.gene, []).append(a)
    return out


def copy_number_map(alterations: list[GeneAlteration]) -> dict[str, float]:
    """Gene to absolute copies, for the stages that want a number.

    Only alterations that carry an absolute count are included. A discrete
    GISTIC call is a call, not a copy count, and turning one into a number
    would be inventing a measurement.
    """
    return {a.gene: a.copies for a in alterations if a.copies is not None}


def suggest_cancer_type(
    alterations: list[GeneAlteration], knowledge: dict | None = None, top: int = 5
) -> list[dict[str, Any]]:
    """Rank cancer types by how well this tumour's events match per-type frequencies.

    The same enrichment the mutated-driver profile gives, over the other two
    kinds of alteration. An ERBB2 amplification is far more informative about
    the cancer type than an ERBB2 mutation, because the amplification is
    concentrated in two cancer types and the mutation is not.
    """
    knowledge = knowledge if knowledge is not None else (load_result(KNOWLEDGE) or {})
    kgenes = knowledge.get("genes", {})
    types = knowledge.get("cancer_types", {})
    scores: dict[str, float] = {}
    for a in alterations:
        k = (kgenes.get(a.gene) or {}).get(a.kind) or {}
        base = k.get("frequency")
        if not base:
            continue
        for ct, f in (k.get("by_cancer_type") or {}).items():
            if ct in types:
                scores[ct] = scores.get(ct, 0.0) + math.log((f + 0.005) / (base + 0.005))
    ranked = sorted(scores.items(), key=lambda kv: -kv[1])[:top]
    return [
        {
            "cancer_type": ct,
            "log_enrichment": round(s, 3),
            "samples_in_study": types[ct],
            "evidence": (
                "inferred from per-type copy-number and structural-variant frequencies (cBioPortal); "
                "suggestion, not diagnosis"
            ),
        }
        for ct, s in ranked
    ]
