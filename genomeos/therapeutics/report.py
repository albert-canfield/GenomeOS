"""The therapeutic analysis report, for a person and for a machine.

Two outputs from the same analysis. The text report is written to be argued
with: every claim is followed by where it came from and what is missing. The
machine report is the compact structure a later autonomous agent consumes.

Both carry the same disclaimer, and neither ever says a therapy will work.
"""

from __future__ import annotations

import textwrap
from typing import Any

from .design import design_readiness, recognition_mode, specification
from .logic import describe
from .mechanisms import MECHANISMS
from .model import TherapeuticTargetCandidate
from .pipeline import DISCLAIMER

RULE = "=" * 78
THIN = "-" * 78


def _wrap(text: str, indent: str = "  ", width: int = 78) -> list[str]:
    """Wrap to the report width, keeping continuation lines aligned under the text."""
    return textwrap.wrap(text, width=width, initial_indent=indent, subsequent_indent=" " * len(indent)) or [
        indent.rstrip()
    ]


def _bullets(items: list[str], indent: str = "  - ", empty: str = "  (none)") -> list[str]:
    out: list[str] = []
    for x in items:
        out.extend(_wrap(str(x), indent))
    return out or [empty]


def _pct(v: float | None) -> str:
    return "n/a" if v is None else f"{v:.0%}"


def _num(v: float | None) -> str:
    return "unknown" if v is None else f"{v:.2f}"


def group(candidates: list[TherapeuticTargetCandidate]) -> dict[str, list[TherapeuticTargetCandidate]]:
    out: dict[str, list[TherapeuticTargetCandidate]] = {
        "direct_surface": [],
        "neoantigen_hla": [],
        "pathway_induced_surface": [],
        "secreted": [],
        "intracellular_only": [],
        "unknown": [],
        "unsuitable": [],
    }
    for c in candidates:
        out.setdefault(c.target_class, []).append(c)
    return out


def _plural(n: int, one: str, many: str = "") -> str:
    return f"{n} {one if n == 1 else (many or one + 's')}"


def overview(analysis: dict[str, Any]) -> list[str]:
    g = group(analysis["candidates"])
    candidates = analysis["candidates"]
    variants = sum(len(c.origins) for c in candidates)
    altered = sum(1 for c in candidates if c.origins)
    indirect = len(candidates) - altered
    neo = [c for c in candidates if c.neoantigen and c.neoantigen.applicable and c.neoantigen.peptides]
    return [
        f"{_plural(variants, 'somatic alteration')} in {_plural(altered, 'gene')}"
        + (f", plus {indirect} gene(s) reached indirectly" if indirect else ""),
        _plural(len(g["direct_surface"]), "potential direct surface target"),
        f"{_plural(len(neo), 'gene')} yielding mutation-derived peptides for the HLA route",
        _plural(len(g["pathway_induced_surface"]), "indirect / pathway-induced candidate"),
        f"{len(g['intracellular_only']) + len(g['unknown']) + len(g['secreted'])} "
        "unsuitable, secreted or unknown",
    ]


def candidate_block(c: TherapeuticTargetCandidate) -> list[str]:
    mode, mode_why = recognition_mode(c)
    lines = [
        THIN,
        f"TARGET: {c.gene}" + (f"  ({c.protein})" if c.protein else ""),
        THIN,
        "",
        "Classification:",
        f"  {c.target_class.replace('_', ' ')}",
        f"  {c.class_reason}",
        "",
        "Why interesting:",
        f"  {c.why_interesting}",
        "",
        "Tumour evidence:",
    ]
    for o in c.origins:
        lines.append(
            f"  {o.chromosome}:{o.position} {o.reference}>{o.alternate}  {o.variant_type}"
            + (f"  {o.protein_change}" if o.protein_change else "")
            + (f"  driver in {o.driver_frequency:.1%} of tumours" if o.driver_frequency else "")
            + ("  hotspot" if o.hotspot else "")
        )
    if not c.origins:
        lines.append("  not altered in this tumour; proposed indirectly (see limitations)")
    lines.append(f"  tumour expression: {_measure(c.tumour.expression)}")
    lines.append(
        f"  clonality: {c.tumour.clonality}"
        + (f" (VAF {c.tumour.vaf:.2f})" if c.tumour.vaf is not None else f" - {c.tumour.clonality_reason}")
    )
    lines += [
        "",
        "Normal tissue concern:",
        f"  risk: {c.normal_tissue.on_target_off_tumour_risk}",
        f"  {c.normal_tissue.risk_basis or c.normal_tissue.summary}",
    ]
    if c.normal_tissue.tissues_at_risk:
        lines.append("  tissues above the concern threshold: " + ", ".join(c.normal_tissue.tissues_at_risk))
    lines += [
        "",
        "Surface accessibility:",
        f"  {c.localization.primary.replace('_', ' ')}, confidence {c.localization.confidence:.2f}"
        + (f", {c.localization.topology}" if c.localization.topology != "unknown" else ""),
    ]
    for r in c.localization.extracellular_regions[:3]:
        lines.append(f"  extracellular {r.start}-{r.end}")
    for r in c.localization.transmembrane_regions[:3]:
        lines.append(f"  transmembrane {r.start}-{r.end}")
    lines += [
        "",
        "Internalisation:",
        f"  internalises: {_tri(c.trafficking.internalises)}"
        + (f"  endosomal: {_tri(c.trafficking.endosomal)}" if c.trafficking.endosomal is not None else "")
        + (f"  lysosomal: {_tri(c.trafficking.lysosomal)}" if c.trafficking.lysosomal is not None else "")
        + (f"  shedding: {_tri(c.trafficking.shedding)}" if c.trafficking.shedding is not None else ""),
    ]
    if c.trafficking.reason:
        lines.append(f"  {c.trafficking.reason}")
    lines += [
        "",
        "Structure:",
        f"  experimental entries: {len(c.structure.pdb_ids)}"
        + (f"  predicted: {c.structure.predicted_structure_id}" if c.structure.predicted_structure_id else "")
        + (
            f"  extracellular region resolved: {_tri(c.structure.extracellular_structure_available)}"
            if c.structure.extracellular_structure_available is not None
            else ""
        ),
    ]
    for e in c.structure.candidate_epitopes[:2]:
        lines.append(
            f"  epitope {e.tumour_specific_change} at residue {e.protein_position}: "
            f"extracellular={_tri(e.extracellular)}, structurally resolved={_tri(e.structurally_resolved)}"
        )
        if e.reference_context:
            lines.append(f"    reference {e.reference_context} -> mutant {e.mutant_context}")
    if c.neoantigen:
        lines += [
            "",
            "Neoantigen / HLA route:",
            f"  novel peptide sequence: {_tri(c.neoantigen.novel_peptide_sequence)}",
            f"  {c.neoantigen.reason}",
        ]
        if c.neoantigen.peptides:
            best = [p for p in c.neoantigen.peptides if p.length == 9] or c.neoantigen.peptides
            lines.append(
                f"  {len(c.neoantigen.peptides)} mutation-spanning peptides; example "
                f"{best[0].wild_type_sequence} -> {best[0].sequence}"
            )
        lines.append("  HLA alleles: " + (", ".join(c.neoantigen.hla_alleles) or "not supplied"))
        lines.append(
            f"  predicted binding: {'available' if c.neoantigen.predicted_binding else 'unavailable'}; "
            f"observed presentation: "
            f"{'yes' if c.neoantigen.observed_immunopeptidomics else 'not supplied'}"
        )
    lines += ["", "Scores (research priority, not clinical):"]
    for key, comp in c.scores.components.items():
        shown = "unknown" if comp.value is None else f"{comp.value:.2f}"
        lines.append(f"  {key:<26} {shown}")
        why = comp.unknown_reason if comp.value is None else comp.basis
        if why:
            lines.extend(_wrap(why, "      "))
    lines.append(
        f"  {'OVERALL':<26} {_num(c.scores.overall)}       "
        f"coverage {_pct(c.scores.coverage)}, confidence {_num(c.scores.confidence)}"
    )
    for a in c.scores.adjustments:
        lines.append(f"  adjustment: {a['adjustment']} {a['from']} -> {a['to']} ({a['reason']})")

    viable = [m for m in c.therapeutic_mechanisms if m.viable and m.compatibility > 0]
    established = [m for m in viable if m.status != "experimental"]
    experimental = [m for m in viable if m.status == "experimental"]
    lines += ["", "Compatible mechanisms (computational compatibility score):"]
    for i, m in enumerate(established[:6], 1):
        lines.append(
            f"  {i}. {m.mechanism:<26} {m.compatibility:.0%}  "
            f"cargo={m.cargo.kind}, effector={m.effector.kind}"
        )
    if not established:
        lines.append("  (none passed their hard requirements)")
    lines += ["", "Experimental possibilities:"]
    for m in experimental[:4]:
        lines.append(f"  - {m.mechanism:<26} {m.compatibility:.0%}  EXPERIMENTAL")
    if not experimental:
        lines.append("  (none)")
    blocked = [m for m in c.therapeutic_mechanisms if not m.viable]
    if blocked:
        lines += ["", "Ruled out, and why:"]
        seen: set[str] = set()
        for m in blocked:
            reason = m.gates_failed[0]
            if reason in seen:
                continue
            seen.add(reason)
            others = [x.mechanism for x in blocked if x.gates_failed and x.gates_failed[0] == reason]
            lines.append(f"  - {', '.join(others)}:")
            lines.extend(_wrap(reason, "      "))
    lines += ["", "Major uncertainties:"] + _bullets(c.limitations[:6])
    readiness = design_readiness(c)
    lines += [
        "",
        f"Design readiness: {readiness['overall']:.2f} (coverage {_pct(readiness['component_coverage'])})",
    ]
    lines += _bullets(readiness["blocking_unknowns"][:5], "  blocked by: ")
    lines += ["", "Evidence:"]
    for level, items in c.ledger.by_level().items():
        lines.append(f"  {level}:")
        for e in items[:4]:
            lines.extend(_wrap(f"[{e.source}] {e.claim}", "    "))
    lines += ["", f"Recognition mode for design: {mode}", f"  {mode_why}", ""]
    return lines


def _tri(v: bool | None) -> str:
    return "unknown" if v is None else ("yes" if v else "no")


def _measure(m: Any) -> str:
    if m.value is not None:
        return f"{m.value:g} {m.unit} ({m.source})"
    if m.qualitative:
        return f"{m.qualitative} ({m.source})"
    return f"unavailable - {m.reason}"


def text_report(analysis: dict[str, Any], detail: int = 5) -> str:
    """The full human-readable therapeutic analysis."""
    candidates: list[TherapeuticTargetCandidate] = analysis["candidates"]
    g = group(candidates)
    lines = [
        RULE,
        "THERAPEUTIC TARGET ANALYSIS",
        RULE,
        "",
        DISCLAIMER,
        "",
        f"Sample: {analysis['sample']}   ({analysis['sample_id']})",
        f"Data level reached: {analysis['data_level']['level_reached']} of 9",
        *_wrap(analysis["data_level"]["note"], "  "),
        "",
        "TUMOUR SUMMARY",
        THIN,
    ]
    lines += _bullets(overview(analysis))
    lines += ["", "TOP DIRECT SURFACE TARGETS", THIN]
    lines += _summary_table(g["direct_surface"])
    lines += ["", "TOP NEOANTIGEN / HLA CANDIDATES", THIN]
    neo = [c for c in candidates if c.neoantigen and c.neoantigen.peptides]
    if neo:
        for c in neo:
            lines.append(
                f"  {c.gene:<10} {len(c.neoantigen.peptides)} peptides, "
                f"HLA {'supplied' if c.neoantigen.hla_alleles else 'MISSING'}, "
                f"presentation {'predicted' if c.neoantigen.predicted_binding else 'unestablished'}"
            )
    else:
        lines.append("  (no variant in this sample yields a mutation-derived peptide)")
    for c in candidates:
        if c.neoantigen and c.neoantigen.applicable is False:
            lines.append(f"  {c.gene:<10} excluded: {c.neoantigen.reason}")
    lines += ["", "INDIRECT / PATHWAY-INDUCED SURFACE TARGETS", THIN]
    lines += (
        _summary_table(g["pathway_induced_surface"])
        if g["pathway_induced_surface"]
        else ["  (none proposed)"]
    )
    lines += ["", "REJECTED OR LOW-PRIORITY TARGETS", THIN]
    rejected = g["intracellular_only"] + g["unknown"] + g["secreted"] + g["unsuitable"]
    for c in rejected:
        lines.append(f"  {c.gene:<10} {c.target_class:<22} {c.class_reason}")
    if not rejected:
        lines.append("  (none)")
    lines += ["", "NORMAL-TISSUE SAFETY CONCERNS", THIN]
    any_risk = False
    for c in candidates:
        if c.normal_tissue.tissues_at_risk:
            any_risk = True
            lines.append(
                f"  {c.gene:<10} {c.normal_tissue.on_target_off_tumour_risk:<9} "
                + ", ".join(c.normal_tissue.tissues_at_risk[:5])
            )
        elif c.normal_tissue.known:
            lines.append(
                f"  {c.gene:<10} {c.normal_tissue.on_target_off_tumour_risk:<9} {c.normal_tissue.summary}"
            )
        else:
            lines.append(f"  {c.gene:<10} unknown   healthy-tissue expression unavailable")
    if not any_risk:
        lines.append("  no candidate exceeds the concern threshold in a queried healthy tissue")
    lines += ["", "TRAFFICKING / INTERNALISATION", THIN]
    for c in candidates:
        lines.append(
            f"  {c.gene:<10} internalises {_tri(c.trafficking.internalises):<8}"
            f"endosomal {_tri(c.trafficking.endosomal):<8}"
            f"lysosomal {_tri(c.trafficking.lysosomal):<8}"
            f"shedding {_tri(c.trafficking.shedding)}"
        )
    lines += ["", "THERAPEUTIC MECHANISM RECOMMENDATIONS", THIN]
    ranked = _mechanism_ranking(candidates)
    for gene, mechanism, score, status, cargo, effector in ranked[:10]:
        flag = "  EXPERIMENTAL" if status == "experimental" else ""
        lines.append(f"  {gene:<10} {mechanism:<26} {score:.0%}  cargo={cargo:<22} effector={effector}{flag}")
    if not ranked:
        lines.append("  (no mechanism passed its hard requirements for any candidate)")
    lines += [
        "",
        "  Scores above are computational compatibility scores. They are not response probabilities.",
    ]
    precedent_rows = []
    for c in candidates:
        hits = [m for m in c.therapeutic_mechanisms if m.viable and (m.precedent or {}).get("approved")]
        if not hits:
            continue
        examples = sorted({e for m in hits for e in m.precedent["examples"]})
        mechanisms = sorted({m.mechanism for m in hits})
        precedent_rows.append((c.gene, examples, mechanisms))
    if precedent_rows:
        lines += ["", "  Existing therapeutic precedent:"]
        for gene, examples, mechanisms in precedent_rows[:6]:
            lines.extend(
                _wrap(
                    f"{gene}: approved agents already engage this target ({', '.join(examples)}), which "
                    f"supports {', '.join(mechanisms)} as reachable modalities. An existing drug against "
                    "this gene does not mean it suits this patient's tumour, whose alteration, "
                    "expression and disease context may differ entirely.",
                    "  - ",
                )
            )
    lines += ["", "COMBINATION TARGET LOGIC", THIN]
    combos = analysis.get("combinations") or []
    for x in combos:
        lines.append(f"  {describe(x)}: {x.rationale}")
        if x.normal_tissues_excluded:
            lines.append(f"    healthy tissues excluded: {', '.join(x.normal_tissues_excluded)}")
    if not combos:
        lines.append(
            "  (no combination assessed: needs two or more surface candidates with healthy-tissue data)"
        )
    lines += ["", "EXPERIMENTAL MECHANISMS", THIN]
    for mechanism, spec in MECHANISMS.items():
        if spec.status != "experimental":
            continue
        best = max(
            (
                (c.gene, m.compatibility)
                for c in candidates
                for m in c.therapeutic_mechanisms
                if m.mechanism == mechanism and m.viable
            ),
            key=lambda t: t[1],
            default=None,
        )
        where = f"best fit {best[0]} at {best[1]:.0%}" if best else "no candidate satisfies its requirements"
        lines.append(f"  {mechanism:<30} {where}")
        lines.append(f"    {spec.summary}")
        for n in spec.notes[:2]:
            lines.append(f"    - {n}")
    lines += ["", "REQUIRED MISSING DATA", THIN]
    for m in analysis["missing_data"]:
        lines.append(f"  - {m['input']}: {m['why']}")
    lines += ["", "CANDIDATE DETAIL", RULE, ""]
    for c in candidates[:detail]:
        lines += candidate_block(c)
    lines += [
        RULE,
        "LIMITATIONS",
        RULE,
        "",
    ]
    lines += _bullets(global_limitations(analysis))
    lines += ["", DISCLAIMER, ""]
    return "\n".join(lines)


def _summary_table(candidates: list[TherapeuticTargetCandidate]) -> list[str]:
    if not candidates:
        return ["  (none)"]
    out = [
        f"  {'gene':<10}{'score':<8}{'access':<8}{'select':<8}{'safety':<8}{'intern':<8}{'evid':<7}"
        "best mechanism"
    ]
    for c in candidates:
        best = c.best_mechanism
        out.append(
            f"  {c.gene:<10}"
            f"{_num(c.scores.overall):<8}"
            f"{_num(c.scores.value('surface_accessibility')):<8}"
            f"{_num(c.scores.value('tumour_selectivity')):<8}"
            f"{_num(c.scores.value('normal_tissue_safety')):<8}"
            f"{_num(c.scores.value('internalisation')):<8}"
            f"{_num(c.scores.value('evidence_strength')):<7}"
            + (f"{best.mechanism} {best.compatibility:.0%}" if best else "none")
        )
    return out


def _mechanism_ranking(candidates: list[TherapeuticTargetCandidate]) -> list[tuple]:
    rows = []
    for c in candidates:
        for m in c.therapeutic_mechanisms:
            if m.viable and m.compatibility > 0:
                rows.append((c.gene, m.mechanism, m.compatibility, m.status, m.cargo.kind, m.effector.kind))
    rows.sort(key=lambda r: -r[2])
    return rows


def global_limitations(analysis: dict[str, Any]) -> list[str]:
    level = analysis["data_level"]["level_reached"]
    out = [
        f"This analysis ran at data level {level} of 9. Conclusions that need a higher level are "
        "reported as unavailable, not estimated.",
        "Localisation, topology, trafficking and healthy-tissue expression are curated population-level "
        "knowledge about the gene, not measurements of this patient's tumour.",
        "No matched normal sample means somatic status is estimated from population allele frequency; "
        "rare inherited variants cannot be separated from somatic ones.",
        "Mechanism compatibility scores describe how well the biology fits a mechanism class. They are "
        "not probabilities that a patient responds.",
        "Correcting or replacing a single gene does not restore a normal cell: established malignancies "
        "carry several cooperating alterations.",
        "GenomeOS produces no therapeutic sequence, construct, vector, formulation or protocol.",
    ]
    if not any(c.tumour.expression.known for c in analysis["candidates"]):
        out.insert(
            1,
            "No tumour RNA or protein data was supplied, so for every candidate the chain is: DNA "
            "evidence present, RNA evidence unavailable, protein evidence unavailable, and surface "
            "expression cannot be established.",
        )
    return out


def target_specification_text(c: TherapeuticTargetCandidate, providers: Any = None) -> str:
    """The interface object between cancer genomics and molecular design."""
    s = specification(c, providers)
    lead = c.origins[0] if c.origins else None
    r = s["design_readiness"]
    lines = [
        RULE,
        f"TARGET SPECIFICATION: {c.gene}",
        RULE,
        "",
        "Patient tumour alteration:",
        f"  {lead.chromosome}:{lead.position} {lead.reference}>{lead.alternate} "
        f"{lead.variant_type} {lead.protein_change or ''}"
        if lead
        else "  not altered in this tumour",
        "",
        "Cancer protein state:",
        f"  {c.protein or c.gene} ({c.uniprot or 'no accession'}), "
        + c.localization.primary.replace("_", " "),
        "",
        "Cancer-specific molecular feature:",
        f"  {s['recognition_specification']['recognition_rationale']}",
        "",
        "Physical accessibility:",
        f"  {c.class_reason}",
        "",
        "Required recognition region:",
        f"  {s['positive_targets']['region'] or 'no delimited extracellular region'}",
        "",
        "What must NOT be recognised:",
    ]
    for n in s["negative_targets"][:10]:
        lines.extend(_wrap(f"{n['label']}: {n['reason']}", "  - "))
    lines += [
        "",
        "Normal tissues at risk:",
    ]
    lines += _bullets(
        [f"{e['tissue']} ({e['value']:g} {e['unit']})" for e in s["normal_tissue_exclusions"]["avoid"][:8]]
    )
    lines += [
        "",
        f"Cancer expression: {_measure(c.tumour.expression)}",
        f"Normal expression: {c.normal_tissue.summary}",
        f"Structural evidence: {len(c.structure.pdb_ids)} PDB entries, "
        f"predicted {c.structure.predicted_structure_id or 'none'}",
        f"Internalisation: {_tri(c.trafficking.internalises)}",
        f"Tumour clonality: {c.tumour.clonality}",
        "",
        "Compatible therapeutic mechanisms:",
    ]
    viable = [m for m in c.therapeutic_mechanisms if m.viable and m.compatibility > 0]
    lines += _bullets([f"{m.mechanism} {m.compatibility:.0%} ({m.status})" for m in viable[:6]])
    best = c.best_mechanism
    lines += [
        "",
        f"Preferred mechanism: {best.mechanism if best else 'none'}"
        + (f" - cargo {best.cargo.kind}, effector {best.effector.kind}" if best else ""),
        "",
        f"Design readiness: {r['overall']:.2f}",
    ]
    lines += _bullets(r["blocking_unknowns"], "  blocked by: ")
    lines += ["", "Missing evidence:"]
    lines += _bullets([f"{m.what}: {m.reason}" for m in c.ledger.missing[:8]])
    lines += [
        "",
        f"Confidence: {_num(c.scores.confidence)} over {_pct(c.scores.coverage)} of the scoring dimensions",
        "",
        DISCLAIMER,
        "",
    ]
    return "\n".join(lines)


def machine_report(analysis: dict[str, Any]) -> dict[str, Any]:
    """The compact structure for a later autonomous agent."""
    return {
        "sample": analysis["sample"],
        "sample_id": analysis["sample_id"],
        "disclaimer": DISCLAIMER,
        "data_level": analysis["data_level"]["level_reached"],
        "therapeutic_candidates": [
            {
                "gene": c.gene,
                "protein": c.protein,
                "uniprot": c.uniprot,
                "target_class": c.target_class,
                "why": c.why_interesting,
                "origin_variants": [o.to_dict() for o in c.origins],
                "scores": {
                    "overall": c.scores.overall,
                    **{
                        k: v.value
                        for k, v in c.scores.components.items()
                        if k
                        in (
                            "tumour_selectivity",
                            "surface_accessibility",
                            "normal_tissue_safety",
                            "internalisation",
                            "evidence_strength",
                        )
                    },
                    "coverage": c.scores.coverage,
                    "confidence": c.scores.confidence,
                },
                "recommended_mechanisms": [
                    {
                        "mechanism": m.mechanism,
                        "compatibility": m.compatibility,
                        "payload_required": m.payload_required,
                        "cargo": m.cargo.kind,
                        "effector": m.effector.kind,
                        "status": m.status,
                    }
                    for m in c.therapeutic_mechanisms
                    if m.viable and m.compatibility > 0
                ][:6],
                "neoantigen": (
                    {
                        "applicable": c.neoantigen.applicable,
                        "novel_peptide_sequence": c.neoantigen.novel_peptide_sequence,
                        "peptides": len(c.neoantigen.peptides),
                        "hla_alleles": c.neoantigen.hla_alleles,
                        "presentation": "unestablished"
                        if not c.neoantigen.predicted_binding
                        else "predicted",
                        "reason": c.neoantigen.reason,
                    }
                    if c.neoantigen
                    else None
                ),
                "evidence": [e.to_dict() for e in c.ledger.items],
                "limitations": c.limitations,
            }
            for c in analysis["candidates"]
        ],
        "combinations": [x.to_dict() for x in analysis.get("combinations", [])],
        "missing_data": analysis["missing_data"],
        "providers": analysis["providers"],
    }
