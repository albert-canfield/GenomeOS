# SPDX-License-Identifier: AGPL-3.0-or-later
"""The attributions as a BioLang program: a chromosome's non-coding space made executable.

`compile_chromosome("chr21")` reads what the attribution lane has computed for one
chromosome and writes it as blocks the engine parses and checks:

- every UNKNOWN block of the budget as a `region` whose role is its tier and label,
  with constraint as evidence; the `constrained_unknown` tier keeps `role: unknown`,
  so the program's own count of unknowns is the chromosome's real unknown;
- every regulatory element with a predicted coding target (the constrained-target
  and the uniform enhancer-deletion runs) as an `element` with its `domain`, its
  `targets` and the basis, plus one `rule` (activates or inhibits) with the predicted
  magnitude as strength; the target genes declared once as stubs;
- the CTCF domains those elements sit in, each with a comment saying in which cell
  types the reader finds the node open.

Everything is labelled with the evidence kind it has (curated, inferred, predicted)
and a confidence; predicted evidence is capped at 0.7 as everywhere in GenomeOS. The
program carries `# test:` lines so `bio test` checks that what was compiled is what
was meant. Generated files say so in their header and are not edited by hand.
"""

from __future__ import annotations

import re
import time
from pathlib import Path

from genomeos.results import RESULTS_DIR, load_result

PREDICTED_CAP = 0.7
EVIDENCE_BY_TIER = {
    "structural": ("curated", "RepeatMasker and the assembly, classified by genomeos unknown"),
    "fossil": ("curated", "RepeatMasker family, Zoonomia phyloP over 241 mammals"),
    "regulatory": ("curated", "ENCODE cCRE registry, Zoonomia phyloP over 241 mammals"),
    "constrained_unknown": ("inferred", "Zoonomia phyloP over 241 mammals, 100-vertebrate elements"),
    "neutral": ("inferred", "Zoonomia phyloP over 241 mammals, no registry element"),
}


def ident(s: str) -> str:
    """A BioLang identifier from any label: word characters only, never starting with a digit."""
    out = re.sub(r"[^A-Za-z0-9_]", "_", s)
    return out if not out[:1].isdigit() else f"g_{out}"


def _text(s: str) -> str:
    """Property text the block syntax accepts: no semicolons, braces, quotes or colons."""
    return re.sub(r"[;{}\"]", ",", str(s)).replace(":", ",")


def _human_axis(chrom: str, results_dir: Path = RESULTS_DIR) -> dict[int, dict]:
    """The human constraint axis per block start, from variation_<chrom> when it has been read."""
    r = load_result(f"variation_{chrom}", results_dir) or {}
    return {blk["start"]: blk for blk in r.get("blocks", []) if blk.get("gnocchi")}


def _copies(chrom: str, results_dir: Path = RESULTS_DIR) -> dict[int, dict]:
    """The copy flag per block start, from duplication_<chrom> when it has been read (area J)."""
    r = load_result(f"duplication_{chrom}", results_dir) or {}
    return {blk["start"]: blk for blk in r.get("blocks", []) if blk.get("duplicated_fraction") is not None}


def _region(
    chrom: str, b: dict, human: dict[int, dict] | None = None, copies: dict[int, dict] | None = None
) -> list[str]:
    tier = b["guess"]["tier"]
    kind, source = EVIDENCE_BY_TIER.get(tier, ("inferred", "genomeos budget"))
    ph = b.get("phylop") or {}
    el = b.get("elements") or {}
    facts = [f"class {b['class']}"]
    if ph.get("fraction_above") is not None:
        facts.append(f"constrained {ph['fraction_above'] * 100:.1f}% of {ph['bases']:,} bases")
    if el.get("n"):
        facts.append(f"{el['n']} conserved elements")
    h = (human or {}).get(b["start"])
    if h:
        # the second axis (area J): variation among people, gnomAD Gnocchi per kilobase
        facts.append(f"people {h['gnocchi']['fraction_above'] * 100:.0f}% of kilobases constrained")
        if h.get("case"):
            facts.append(f"case {h['case']['case']}")
    c = (copies or {}).get(b["start"])
    is_copy = bool(c and c["duplicated_fraction"] >= 0.5)
    if c and c["duplicated_fraction"] > 0:
        facts.append(f"duplicated {c['duplicated_fraction'] * 100:.0f}% with {c.get('pairs', 0)} partners")
    role = "unknown" if tier == "constrained_unknown" else _text(f"{tier}, {b['guess']['label']}")
    if is_copy:
        # a copy is read as a copy first, whatever its tier; an unknown that is a copy stays unknown
        role = role if role == "unknown" else _text(f"copy, {role}")
    return [
        f"region U_{chrom}_{b['start']} {{",
        f"  locus: {chrom}:{b['start']}-{b['end']}",
        f"  role: {role}",
        f'  evidence: {kind} "{_text(source)}" {_text(", ".join(facts))}',
        f"  confidence: {min(b['guess']['confidence'], 1.0):.2f}",
        "}",
    ]


def _attributed(chrom: str, results_dir: Path = RESULTS_DIR) -> list[dict]:
    """Elements with a predicted coding target, constrained-target run first, deduplicated by id."""
    seen: set[str] = set()
    out: list[dict] = []
    for name, origin in (("constrained_targets", "constrained"), ("enhancer_targets", "uniform")):
        r = load_result(f"{name}_{chrom}", results_dir) or {}
        for e in r.get("elements", []):
            pc = e.get("predicted_coding") or {}
            if not pc.get("gene") or e["id"] in seen:
                continue
            seen.add(e["id"])
            out.append({**e, "origin": origin})
    out.sort(key=lambda e: e["start"])
    return out


def _open_in(chrom: str, domain_id: str, readers: list[dict]) -> str:
    rows = []
    for r in readers:
        for n in r.get("node_table", []):
            if n["id"] == domain_id:
                rows.append((r["cell_type"], n.get("open_fraction") or 0.0))
    if not rows:
        return ""
    rows.sort(key=lambda x: -x[1])
    opened = [f"{c} ({f:.2f})" for c, f in rows if f >= 0.05]
    silent = [c for c, f in rows if f < 0.05]
    parts = []
    if opened:
        parts.append("open in " + ", ".join(opened[:6]))
    if silent:
        parts.append("silent in " + ", ".join(silent[:6]))
    return "; ".join(parts)


def compile_chromosome(chrom: str, results_dir: Path = RESULTS_DIR) -> str:
    budget = load_result(f"budget_{chrom}", results_dir)
    human = _human_axis(chrom, results_dir)
    copies = _copies(chrom, results_dir)
    if not budget:
        raise FileNotFoundError(f"no budget_{chrom} result; run genomeos budget --chrom {chrom}")
    domains = {d["id"]: d for d in (load_result(f"domains_{chrom}", results_dir) or {}).get("domains", [])}
    readers = [
        load_result(p.stem, results_dir)
        for p in sorted(results_dir.glob(f"reader_*_{chrom}.json"))
        if "_vs_" not in p.stem
    ]
    readers = [r for r in readers if r and r.get("node_table")]
    elements = _attributed(chrom, results_dir)

    lines = [
        f"module human.noncoding.{chrom}",
        "",
        f"# The non-coding space of {chrom} as attributions, generated on {time.strftime('%Y-%m-%d')} by",
        f"# `genomeos budget --chrom {chrom} --bio` from budget_{chrom}, constrained_targets_{chrom},",
        f"# enhancer_targets_{chrom}, domains_{chrom} and reader_*_{chrom}. Do not edit; rerun the",
        "# command. docs/ATTRIBUTION.md explains the tiers and the evidence.",
        "#",
        "# A region's role is its budget tier; the constrained_unknown tier keeps role unknown, so this",
        "# program's own count of unknowns is the chromosome's real unknown. An element carries the gene",
        "# AlphaGenome says it reaches when deleted, the tissue and the magnitude, as predicted evidence.",
    ]
    regions = [b for b in sorted(budget["blocks"], key=lambda b: b["start"])]
    n_unknown = sum(1 for b in regions if b["guess"]["tier"] == "constrained_unknown")
    genes = sorted({e["predicted_coding"]["gene"] for e in elements})
    used_domains = sorted({e["domain"] for e in elements if e.get("domain") in domains})
    n_entities = len(regions) + len(genes) + len(elements) + len(used_domains)
    lines += [
        "#",
        f"# test: entities >= {n_entities}",
        f"# test: unknowns == {n_unknown}",
        f"# test: rules == {len(elements)}",
        "",
        f"# ---- the budget: {len(regions)} UNKNOWN blocks, {budget['unknown_bp'] / 1e6:.1f} Mb, "
        f"constrained {((budget.get('constrained_fraction') or 0) * 100):.2f}% of measured bases",
    ]
    for b in regions:
        lines += _region(chrom, b, human, copies)
    if used_domains:
        lines += [
            "",
            f"# ---- the nodes the attributed elements sit in ({len(used_domains)}), with the reader's view",
        ]
        for did in used_domains:
            d = domains[did]
            note = _open_in(chrom, did, readers)
            if note:
                lines.append(f"# {did}: {note}")
            lines.append(
                f"domain {ident(did)} {{ locus: {chrom}:{d['start']}-{d['end']}; "
                f'evidence: inferred "CTCF-only ENCODE elements as boundary proxies, no Hi-C"; '
                f"confidence: {d.get('confidence', 0.4):.2f} }}"
            )
    if genes:
        lines += [
            "",
            f"# ---- target genes ({len(genes)}), declared once so elements and rules can name them",
        ]
        for g in genes:
            lines.append(
                f'gene {ident(g)} {{ symbol: {g}; evidence: curated "GENCODE v50"; confidence: 0.9 }}'
            )
    if elements:
        lines += ["", f"# ---- attributed elements ({len(elements)}): deletion in AlphaGenome names the gene"]
        for e in elements:
            pc = e["predicted_coding"]
            conf = round(
                min(PREDICTED_CAP, max(0.05, float(pc.get("confidence") or abs(pc["log2_fold_change"])))), 2
            )
            action = "activates" if pc.get("action") == "activates" else "inhibits"
            basis = (
                f"predicted, deleting the element moves {pc['gene']} by {pc['log2_fold_change']:+.2f} log2 "
                f"in {pc.get('tissue') or 'the strongest track'} ({pc.get('strength') or 'weak'})"
            )
            if e.get("constrained_fraction") is not None:
                basis += f", constrained {e['constrained_fraction'] * 100:.0f}% of bases (Zoonomia)"
            if e.get("verdict_coding"):
                basis += f", {e['verdict_coding']}"
            props = ["class: enhancer", f"locus: {chrom}:{e['start']}-{e['end']}"]
            if e.get("domain") in domains:
                props.append(f"domain: {ident(e['domain'])}")
            props += [
                f"targets: {ident(pc['gene'])}",
                f"basis: {_text(basis)}",
                'evidence: predicted "AlphaGenome RNA-seq gene scorer, expression change on deletion"',
                f"confidence: {conf:.2f}",
            ]
            lines.append(f"element {e['id']} {{")
            lines += [f"  {p}" for p in props]
            lines.append("}")
            strength = round(min(1.0, abs(float(pc["log2_fold_change"]))), 3)
            lines.append(
                f"rule {e['id']} {action} {ident(pc['gene'])} {{ strength: {strength}; "
                f'evidence: predicted "AlphaGenome deletion, {_text(pc.get("tissue") or "strongest track")}";'
                " "
                f"confidence: {conf:.2f} }}"
            )
    return "\n".join(lines) + "\n"


def write_program(chrom: str, out: Path | None = None, results_dir: Path = RESULTS_DIR) -> Path:
    out = out or Path("data/organisms/human") / f"noncoding_{chrom}.bio"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(compile_chromosome(chrom, results_dir))
    return out
