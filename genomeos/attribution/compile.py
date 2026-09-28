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

Everything is labelled with the evidence kind it has (curated, inferred, predicted).
Since review R4 (2026-09-28) a predicted element or rule states no `confidence:` (the parser
reads that as 0.0, unstated): it carries its effect in its unit and says its probability is
unavailable, and the full `genomeos.certainty.Certainty` record is written once in the
header. A region keeps the budget's hand-set evidence-quality score, labelled as such in its
evidence note; gene stubs and domains state no constant. The program carries `# test:`
lines so `bio test` checks that what was compiled is what was meant. Generated files say so
in their header and are not edited by hand.

Since 2026-09-17 the program also carries an **experimental layer**
(`attribution/measured.py`): where a CRISPRi screen, a lentiMPRA library, a VISTA
transgenic assay or a saturation-mutagenesis experiment measured the same piece of DNA as
a compiled element, a second block `<id>_measured` is written beside it, carrying
`evidence: experimental` and the source.
Two names, never one: the predicted block is not touched, and a measured negative - "the
screen found no effect on this gene", "none of the bases measured here matters" - is
stated in the measured block rather than dropped or turned back into UNKNOWN. The
base-level assay is stated as a fact about bases and never as a verdict on the element.
The layer is thin on purpose; the census in `measured_layer_genome` says how thin.
"""

from __future__ import annotations

import bisect
import gzip
import re
import time
from pathlib import Path
from typing import Any

from genomeos.certainty import Certainty
from genomeos.genome.repeats import INTERSPERSED as _INTERSPERSED
from genomeos.results import RESULTS_DIR, load_result

#: retired as a compiled number by review R4 (2026-09-28); kept because
#: attribution/confidence_calibration.py recomputes the historical formula to read old programs
PREDICTED_CAP = 0.7
EFFECT_UNIT = (
    "log2 fold change of the target gene's predicted RNA-seq expression on deleting the element"
    " (AlphaGenome gene scorer, one track)"
)
MODEL_SCORE_NAME = (
    "the target run's `confidence` field, equal to |log2 fold change|: the magnitude the run ranked"
    " by, not a probability"
)
NO_PROBABILITY = (
    "no calibration record: no compiled prediction has been scored against a measured outcome"
    " population by a stated method, so no probability that deleting this element moves this gene is"
    " quoted"
)
REGION_SCORE_NOTE = "confidence is the budget's hand-set evidence-quality score per rule, not a probability"
EVIDENCE_BY_TIER = {
    "structural": ("curated", "RepeatMasker and the assembly, classified by genomeos unknown"),
    "fossil": ("curated", "RepeatMasker family, Zoonomia phyloP over 241 mammals"),
    "regulatory": ("curated", "ENCODE cCRE registry, Zoonomia phyloP over 241 mammals"),
    "constrained_unknown": ("inferred", "Zoonomia phyloP over 241 mammals, 100-vertebrate elements"),
    "neutral": ("inferred", "Zoonomia phyloP over 241 mammals, no registry element"),
}


# ---- review R7 (2026-09-28): five axes kept apart (vocabulary: genomeos.lang.grammar.AXES) ----
#: the biochemical signature the ENCODE registry states, per cCRE class
CCRE_ROLE = {
    "PLS": "promoter_like",
    "pELS": "enhancer_like",
    "dELS": "enhancer_like",
    "CTCF-only": "insulator_like",
    "DNase-H3K4me3": "open_chromatin",
}
#: `class:` is kept as a derived, labelled summary of the registry role (the first one), never of
#: the activity; an element the registry does not hold says `class: unknown`
ROLE_CLASS = {
    "promoter_like": "promoter",
    "enhancer_like": "enhancer",
    "open_chromatin": "open_chromatin",
    "insulator_like": "insulator",
}
#: what a predicted or measured increase on removal can mean; none is chosen, so none is forced
REPRESSION_ALTERNATIVES = ("silencer", "insulator_like", "competing_promoter", "unknown")
#: interspersed-repeat share of an interval from which it reads as repeat-derived (below: partly)
REPEAT_DERIVED_MIN = 0.5
#: duplicated share of a block from which it reads as a segmental duplication (as the copy flag)
DUPLICATED_MIN = 0.5
#: sequence classes of `genomeos unknown` whose origin the class itself states
CLASS_ORIGIN = {
    "gap": "assembly_gap",
    "centromere": "satellite",
    "satellite_array": "satellite",
    "tandem_repeat": "tandem_repeat",
    "telomere": "tandem_repeat",
}
CLASS_ROLE = {
    "promoter_like": [["promoter_like", "unknown"]],
    "centromere": [["structural"]],
    "satellite_array": [["structural", "unknown"]],
    "tandem_repeat": [["structural", "unknown"]],
    "long_orf": [["coding_candidate", "unknown"]],
}
AXIS_ORDER = ("origin", "molecular_role", "activity", "target_relation", "evidence_status")


def _ccres(chrom: str, results_dir: Path) -> dict[str, tuple[str, bool]]:
    """The ENCODE registry's class and CTCF flag per cCRE id, from ccres_<chrom> when present."""
    p = results_dir / f"ccres_{chrom}.bed.gz"
    if not p.exists():
        return {}
    out = {}
    with gzip.open(p, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                f = line.rstrip("\n").split("\t")
                out[f[3]] = (f[4], f[5] == "1")
    return out


class _Interspersed:
    """Interspersed-repeat coverage of any interval, and its largest class, from rmsk_<chrom>."""

    def __init__(self, chrom: str, results_dir: Path) -> None:
        rows: list[tuple[int, int, str]] = []
        p = results_dir / f"rmsk_{chrom}.bed.gz"
        if p.exists():
            with gzip.open(p, "rt") as fh:
                for line in fh:
                    f = line.split("\t", 3)
                    if f[2] in _INTERSPERSED:
                        rows.append((int(f[0]), int(f[1]), f[2]))
        rows.sort()
        self.rows = rows
        self.starts = [r[0] for r in rows]
        self.max_len = max((r[1] - r[0] for r in rows), default=0)

    def __bool__(self) -> bool:
        return bool(self.rows)

    def cover(self, start: int, end: int) -> tuple[float, str]:
        by: dict[str, int] = {}
        for s, e, c in self.rows[bisect.bisect_left(self.starts, start - self.max_len) :]:
            if s >= end:
                break
            ov = min(e, end) - max(s, start)
            if ov > 0:
                by[c] = by.get(c, 0) + ov
        if not by:
            return 0.0, ""
        return min(1.0, sum(by.values()) / max(1, end - start)), max(by, key=lambda k: by[k])


def _origin(fraction: float | None, top: str) -> str:
    if fraction is None:
        return "unknown"
    if fraction >= REPEAT_DERIVED_MIN:
        return f"repeat_derived/{top}" if top else "repeat_derived"
    if fraction > 0:
        return f"partly_repeat_derived/{top}" if top else "partly_repeat_derived"
    return "unique"


def selection(fraction: float | None) -> str:
    """Constraint as evidence of selection, never as a mechanism; its absence is not no function."""
    from genomeos.attribution.budget import CONSTRAINED_MIN, NEUTRAL_MAX

    if fraction is None:
        return "selection_not_measured"
    if fraction >= CONSTRAINED_MIN:
        return "under_selection"
    return "selection_weak" if fraction >= NEUTRAL_MAX else "selection_not_detected"


def registry_roles(cls: str | None, ctcf_bound: bool = False) -> list[str]:
    roles = [CCRE_ROLE[cls]] if cls in CCRE_ROLE else []
    if ctcf_bound and "insulator_like" not in roles:
        roles.append("insulator_like")
    return roles


def derived_class(axes: dict[str, list[list[str]]]) -> str:
    """`class:` as a summary of the first registry role the element holds for certain, else unknown."""
    for group in axes.get("molecular_role", []):
        if len(group) == 1 and group[0] in ROLE_CLASS:
            return ROLE_CLASS[group[0]]
    return "unknown"


def axis_lines(axes: dict[str, list[list[str]]]) -> list[str]:
    """`axis: a, b|c` per stated axis, in the grammar's order: `,` all hold, `|` unresolved."""
    return [f"{k}: {', '.join('|'.join(g) for g in axes[k])}" for k in AXIS_ORDER if axes.get(k)]


def element_axes(e: dict, ccre: dict[str, tuple[str, bool]], rep: _Interspersed | None) -> dict:
    """The five axes of a predicted element (R7): sequence facts from the registry and RepeatMasker,
    the direction from the prediction, and a repression left as unresolved alternatives."""
    cls, ctcf = ccre.get(e["id"], (None, False))
    frac, top = rep.cover(e["start"], e["end"]) if rep else (None, "")
    roles = [[r] for r in registry_roles(cls, ctcf)]
    pc = e["predicted_coding"]
    represses = pc.get("action") != "activates"
    if represses:
        roles.append(list(REPRESSION_ALTERNATIVES))
    relation = [["predicted_deletion_target"]]
    if e.get("verdict_coding") == "agrees with nearest TSS in domain":
        relation.append(["nearest_tss_in_domain"])
    status = (["registry_biochemical"] if cls else []) + ["predicted_model"]
    status.append(selection(e.get("constrained_fraction")))
    return {
        "origin": [[_origin(frac, top)]],
        "molecular_role": roles or [["unknown"]],
        "activity": [["represses_target" if represses else "activates_target"]],
        "target_relation": relation,
        "evidence_status": [[s] for s in status],
    }


def measured_axes(row: dict, sequence: dict | None) -> dict:
    """The five axes of a `_measured` block: the sequence facts of its predicted twin, and activity,
    relation and status from the assays alone."""
    from genomeos.attribution import measured as ms

    seq = sequence or {}
    roles = [g for g in seq.get("molecular_role", []) if g != list(REPRESSION_ALTERNATIVES)]
    m = row["measured"]
    activity: list[list[str]] = []
    relation: list[list[str]] = []
    status = ["measured"]
    c = m.get("crispri")
    if c:
        actions = {a for _, a, _, _, _ in ms.rule_links(row)}
        if "activates" in actions:
            activity.append(["activates_target"])
        if "inhibits" in actions or c.get("genes_increased"):
            activity.append(["represses_target"])
            roles.append(list(REPRESSION_ALTERNATIVES))
        if not activity and c.get("genes_no_effect_well_powered"):
            activity.append(["no_effect_measured"])
        if c.get("genes_regulated_training"):
            relation.append(["measured_perturbation_target"])
        elif c.get("genes_tested") and not c.get("genes_regulated") and not c.get("genes_increased"):
            relation.append(["tested_no_effect"])
        if c.get("genes_not_regulated"):
            status.append("measured_negative")
    lt = m.get("lentimpra")
    if lt:
        if lt.get("cells_active"):
            activity.append(["active_in_reporter"])
        elif lt.get("cells_conflicting"):
            activity.append(["active_in_reporter", "inactive_in_reporter"])
        elif lt.get("cells_silent"):
            activity.append(["inactive_in_reporter"])
    v = m.get("vista")
    if v:
        if v.get("positive"):
            activity.append(["active_in_reporter"])
        elif v.get("negative"):
            activity.append(["inactive_in_reporter"])
    if row.get("agreement", {}).get("verdict") == ms.DISAGREES:
        status.append("conflicting")
    seen: list[list[str]] = []
    for g in activity:
        if g not in seen:
            seen.append(g)
    return {
        "origin": seq.get("origin") or [["unknown"]],
        "molecular_role": roles or [["unknown"]],
        "activity": seen or [["unknown"]],
        "target_relation": relation or [["unassigned"]],
        "evidence_status": [[s] for s in status],
    }


def region_axes(b: dict, features: dict | None, rep: _Interspersed | None, duplicated: float | None) -> dict:
    """The five axes of an UNKNOWN block: origin from its sequence class or RepeatMasker, roles from
    the registry classes it contains, and constraint only as a reading of selection."""
    cls = b["class"]
    f = features or {}
    status: list[str] = []
    if cls in CLASS_ORIGIN or cls.startswith("interspersed_repeat"):
        top = cls.removeprefix("interspersed_repeat").lstrip("_")
        origin = CLASS_ORIGIN.get(cls) or (f"repeat_derived/{top}" if top else "repeat_derived")
        if b.get("class_evidence") == "curated":
            status.append("curated_annotation")
    elif f.get("interspersed_coverage") is not None:
        cover = {k: v for k, v in (f.get("repeat_coverage") or {}).items() if k in _INTERSPERSED}
        origin = _origin(f["interspersed_coverage"], max(cover, key=lambda k: cover[k]) if cover else "")
        status.append("curated_annotation")
    elif rep:
        origin = _origin(*rep.cover(b["start"], b["end"]))
        status.append("curated_annotation")
    else:
        origin = "unknown"
    origins = [[origin]]
    if duplicated is not None and duplicated >= DUPLICATED_MIN:
        origins.append(["segmental_duplication"])
    if cls in ("regulatory",) and f.get("ccre"):
        counts = f["ccre"]
        roles: list[str] = []
        for k in sorted(counts, key=lambda k: -counts[k]):
            r = CCRE_ROLE.get(k)
            if counts[k] and r and r not in roles:
                roles.append(r)
        role_groups = [[r] for r in roles] or [["unknown"]]
        status.insert(0, "registry_biochemical")
    else:
        role_groups = CLASS_ROLE.get(cls, [["unknown"]])
    if cls != "gap":
        status.append(selection((b.get("phylop") or {}).get("fraction_above")))
    return {
        "origin": origins,
        "molecular_role": role_groups,
        "activity": [["unknown"]],
        "target_relation": [["unassigned"]],
        "evidence_status": [[s] for s in (status or ["unknown"])],
    }


def ident(s: str) -> str:
    """A BioLang identifier from any label: word characters only, never starting with a digit."""
    out = re.sub(r"[^A-Za-z0-9_]", "_", s)
    return out if not out[:1].isdigit() else f"g_{out}"


def _text(s: str) -> str:
    """Property text the block syntax accepts: no semicolons, braces, quotes or colons."""
    return re.sub(r"[;{}\"]", ",", str(s)).replace(":", ",")


def context(name: str | None) -> str:
    """The `when: cell_type = ...` value for a cell or biosample name, or `unknown` when none was
    recorded (R1, 2026-09-28): a missing context stays explicit and never makes a rule universal."""
    from genomeos.attribution.measured import CONTEXT_UNKNOWN

    return ident(name) if name and name.strip() else CONTEXT_UNKNOWN


def _human_axis(chrom: str, results_dir: Path = RESULTS_DIR) -> dict[int, dict]:
    """The human constraint axis per block start, from variation_<chrom> when it has been read."""
    r = load_result(f"variation_{chrom}", results_dir) or {}
    return {blk["start"]: blk for blk in r.get("blocks", []) if blk.get("gnocchi")}


def _copies(chrom: str, results_dir: Path = RESULTS_DIR) -> dict[int, dict]:
    """The copy flag per block start, from duplication_<chrom> when it has been read (area J)."""
    r = load_result(f"duplication_{chrom}", results_dir) or {}
    return {blk["start"]: blk for blk in r.get("blocks", []) if blk.get("duplicated_fraction") is not None}


def element_certainty(pc: dict[str, Any]) -> Certainty:
    """What a predicted element-to-gene link rests on, R4 (2026-09-28): effect in its unit, the run's
    ranking score named, no probability. Nothing here converts the effect into a certainty."""
    score = pc.get("confidence")
    return Certainty(
        evidence_category="predicted: AlphaGenome deletion, one model run",
        effect_estimate=float(pc["log2_fold_change"]),
        effect_unit=EFFECT_UNIT,
        uncertainty_note="one deterministic model run; no spread computed",
        model_score=None if score is None else float(score),
        model_score_name=MODEL_SCORE_NAME,
        probability_unavailable=NO_PROBABILITY,
    )


def _effect_note(pc: dict[str, Any]) -> str:
    """The per-block form of the record: the effect with its unit's short name, no probability."""
    return (
        f"effect {float(pc['log2_fold_change']):+.3g} log2 fold change on deletion, probability unavailable"
    )


def _region(
    chrom: str,
    b: dict,
    human: dict[int, dict] | None = None,
    copies: dict[int, dict] | None = None,
    features: dict | None = None,
    rep: _Interspersed | None = None,
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
        f'  evidence: {kind} "{_text(source)}" {_text(", ".join([*facts, REGION_SCORE_NOTE]))}',
        f"  confidence: {min(b['guess']['confidence'], 1.0):.2f}",
        *(f"  {ln}" for ln in axis_lines(region_axes(b, features, rep, c and c["duplicated_fraction"]))),
        "}",
    ]


def _attributed(chrom: str, results_dir: Path = RESULTS_DIR) -> list[dict]:
    """Elements with a predicted coding target, the whole-chromosome run first, deduplicated by id."""
    from genomeos.attribution.targets import attributed

    out = attributed(chrom, results_dir)
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


def _measured_rows(
    chrom: str, elements: list[dict], results_dir: Path, layer: Any = None
) -> tuple[Any, list[dict]]:
    """The layer and its rows, or an empty layer if no assay cache is on this machine.

    A missing cache must read as "no measurement was available", never as "nothing was measured":
    the caller gets an empty list and the program then says, in its own header, that it carries no
    experimental facts.
    """
    from genomeos.attribution.measured import Layer, rows

    layer = layer if layer is not None else Layer.load(chrom)
    return layer, ([] if layer.empty else rows(chrom, elements, layer, results_dir=results_dir))


def link_pairs(row: dict) -> dict[tuple[str, str], int]:
    """How many regulated pairs each experimental rule's strength was taken from, per (gene, cell).

    `measured.rule_links` takes a link's strength as the largest |EffectSize| among the regulated
    training pairs of one (element, gene, cell), or among the held-out ones when there is no training
    pair (lane-assay's census A8). With one pair that is an observation; with more it is a maximum. The
    count is written on the measured element so a reader, and `attribution.bridge`, can tell (R3).
    """
    from genomeos.attribution import measured as ms

    pairs = row["measured"].get("crispri", {}).get("pairs", [])
    out: dict[tuple[str, str], int] = {}
    for gene, cell in sorted({(p["gene"], p["cell"]) for p in pairs if p["regulated"]}):
        hit = [p for p in pairs if p["gene"] == gene and p["cell"] == cell and p["regulated"]]
        train = [p for p in hit if p.get("split", ms.TRAINING) == ms.TRAINING]
        out[(gene, cell)] = len(train or hit)
    return out


def _measured_blocks(
    chrom: str, row: dict, domains: dict, ident_of: dict, sequence: dict | None = None
) -> list[str]:
    """One `<id>_measured` element and one rule per measured regulated link, all experimental."""
    from genomeos.attribution import measured as ms

    axes = measured_axes(row, sequence)
    props = [f"class: {derived_class(axes)}", f"locus: {chrom}:{row['start']}-{row['end']}"]
    if row.get("domain") in domains:
        props.append(f"domain: {ident(row['domain'])}")
    c = row["measured"].get("crispri", {})
    # a held-out-only link never reaches `targets:`, the field a feature reads; it stays as a marked
    # rule below (the split audit of 2026-09-28 in docs/ATTRIBUTION.md)
    regulated = c.get("genes_regulated_training", c.get("genes_regulated", []))
    if regulated:
        props.append("targets: " + ", ".join(ident_of[g] for g in regulated))
    props += [
        f"basis: {_text(ms.basis_text(row))}",
        f'evidence: experimental "{_text(ms.sources_of(row["measured"]))}"'
        + "".join(
            (" links " if i == 0 else ", ") + f"{ident_of[g]} in {context(cell)} from {n} pairs"
            for i, ((g, cell), n) in enumerate(link_pairs(row).items())
        ),
        f"confidence: {row['confidence']:.2f}",
    ]
    props += axis_lines(axes)
    lines = [f"element {row['id']}_measured {{}}".replace("{}", "{")]
    lines += [f"  {p}" for p in props]
    lines.append("}")
    for gene, action, strength, cell, split in ms.rule_links(row):
        mark = f", {ms.HELDOUT_MARK}" if split == ms.HELDOUT else ""
        lines.append(
            f"rule {row['id']}_measured {action} {ident_of[gene]} {{ strength: {strength}; "
            f"when: {ms.CONTEXT_KEY} = {context(cell)}; "
            f'evidence: experimental "{_text(ms.SOURCES["crispri"])}, silenced in {_text(cell)}{mark}"; '
            f"confidence: {row['confidence']:.2f} }}"
        )
    return lines


def compile_chromosome(chrom: str, results_dir: Path = RESULTS_DIR, layer: Any = None) -> str:
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
    ccre = _ccres(chrom, results_dir)
    rep = _Interspersed(chrom, results_dir)
    features = {
        blk["start"]: blk.get("features") or {}
        for blk in (load_result(f"unknown_{chrom}", results_dir) or {}).get("blocks", [])
    }
    sequence: dict[str, dict] = {}
    layer, measured_rows = _measured_rows(chrom, elements, results_dir, layer)
    from genomeos.attribution.measured import rule_links

    # one experimental rule per (element, gene, cell), R1 of 2026-09-28
    n_measured_rules = sum(len(rule_links(r)) for r in measured_rows)

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
        "#",
        "# Two blocks per measured element, never one: `<id>` is what the model predicts and",
        "# `<id>_measured` is what an assay measured over the same DNA, with evidence: experimental.",
        "# Neither overwrites the other, and a measured negative - no effect on this gene in this",
        "# screen - is stated in the measured block, because that is a measurement and not an unknown.",
        "#",
        "# One of the four assays is base-level: saturation mutagenesis says which BASES inside an",
        "# element matter. It can support the compiled claim and cannot contradict it, so an element",
        "# whose measured bases are all inert is stated as a fact about those bases and never as a",
        "# verdict on the element, and it raises no rule, because it names no gene.",
        "#",
        "# Every rule is gated on the cell it was measured or predicted in (`when: cell_type = K562`), one",
        "# rule per element, gene and cell, so a run in HepG2 integrates none of K562's. A rule whose",
        "# cell was not recorded says `cell_type = unknown`, which matches no cell: it is never universal.",
        "#",
        "# Certainty (review R4, 2026-09-28). A predicted element or rule states no `confidence:`: the",
        "# size of a predicted effect is not how sure anyone is. Each carries its effect in its unit and",
        "# says its probability is unavailable; the record behind that note, for every predicted link:",
        *(
            "# " + ln[2:]
            for ln in Certainty(
                evidence_category="predicted: AlphaGenome deletion, one model run",
                effect_unit=EFFECT_UNIT,
                uncertainty_note="one deterministic model run; no spread computed",
                model_score_name=MODEL_SCORE_NAME,
                probability_unavailable=NO_PROBABILITY,
            ).comment_lines()
        ),
        "# with the effect and the score per element. A region's `confidence:` is the budget's hand-set",
        "# evidence-quality score for the rule that fired, not a probability; a `_measured` block's is the",
        "# hand-set rank of its strongest assay kind (perturbation above reporter), not a probability.",
        "# Gene stubs and domains state no confidence.",
        "#",
        "# Executable annotation, not a simulation (review R3, 2026-09-28). A rule's `strength` is its",
        "# observation's magnitude in that observation's unit: |log2 fold change| clipped at 1 for a",
        "# predicted deletion, the largest |EffectSize| of the link's CRISPRi pairs for a measured one (the",
        "# measured element says how many pairs). It is not a rate constant. Target genes are stubs with",
        "# no basal or max rate and an element is not a species, so a runtime run of this program reports",
        "# both as unresolved. genomeos.attribution.bridge.parameterize turns one observation per gene and",
        "# cell into one fitted strength, given declared rates, and names every input it lacks.",
        "#",
        "# Five axes, not one label (review R7, 2026-09-28). Elements and regions state `origin`,",
        "# `molecular_role`, `activity`, `target_relation` and `evidence_status` (vocabulary:",
        "# genomeos.lang.grammar.AXES). `,` joins values that all hold, `|` alternatives none of which is",
        "# chosen. A predicted or measured increase on removal is `represses_target` with the role left",
        "# as silencer|insulator_like|competing_promoter|unknown, never a silencer label. Constraint is",
        "# read only as evidence of selection; its absence is `selection_not_detected`, not no function.",
        "# `class:` is a summary of the registry role and a region's `role:` the budget's tier summary.",
    ]
    regions = [b for b in sorted(budget["blocks"], key=lambda b: b["start"])]
    n_unknown = sum(1 for b in regions if b["guess"]["tier"] == "constrained_unknown")
    measured_genes = {
        g for r in measured_rows for g in r["measured"].get("crispri", {}).get("genes_regulated", [])
    }
    genes = sorted({e["predicted_coding"]["gene"] for e in elements} | measured_genes)
    ident_of = {g: ident(g) for g in genes}
    used_domains = sorted({e["domain"] for e in elements if e.get("domain") in domains})
    n_entities = len(regions) + len(genes) + len(elements) + len(used_domains) + len(measured_rows)
    lines += [
        "#",
        f"# test: entities >= {n_entities}",
        f"# test: unknowns == {n_unknown}",
        f"# test: rules == {len(elements) + n_measured_rules}",
        "",
        f"# ---- the budget: {len(regions)} UNKNOWN blocks, {budget['unknown_bp'] / 1e6:.1f} Mb, "
        f"constrained {((budget.get('constrained_fraction') or 0) * 100):.2f}% of measured bases",
    ]
    for b in regions:
        lines += _region(chrom, b, human, copies, features.get(b["start"]), rep)
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
                'evidence: inferred "CTCF-only ENCODE elements as boundary proxies, no Hi-C" '
                "no confidence stated, the domain call has no measured outcome }"
            )
    if genes:
        lines += [
            "",
            f"# ---- target genes ({len(genes)}), declared once so elements and rules can name them",
        ]
        for g in genes:
            lines.append(
                f'gene {ident(g)} {{ symbol: {g}; evidence: curated "GENCODE v50" symbol only, a stub }}'
            )
    if elements:
        lines += ["", f"# ---- attributed elements ({len(elements)}): deletion in AlphaGenome names the gene"]
        for e in elements:
            pc = e["predicted_coding"]
            action = "activates" if pc.get("action") == "activates" else "inhibits"
            basis = (
                f"predicted, deleting the element moves {pc['gene']} by {pc['log2_fold_change']:+.2f} log2 "
                f"in {pc.get('tissue') or 'the strongest track'} ({pc.get('strength') or 'weak'})"
            )
            if e.get("constrained_fraction") is not None:
                basis += f", constrained {e['constrained_fraction'] * 100:.0f}% of bases (Zoonomia)"
            if e.get("verdict_coding"):
                basis += f", {e['verdict_coding']}"
            axes = element_axes(e, ccre, rep)
            sequence[e["id"]] = axes
            props = [f"class: {derived_class(axes)}", f"locus: {chrom}:{e['start']}-{e['end']}"]
            if e.get("domain") in domains:
                props.append(f"domain: {ident(e['domain'])}")
            props += [
                f"targets: {ident(pc['gene'])}",
                f"basis: {_text(basis)}",
                'evidence: predicted "AlphaGenome RNA-seq gene scorer, expression change on deletion" '
                + _text(_effect_note(pc)),
                *axis_lines(axes),
            ]
            lines.append(f"element {e['id']} {{")
            lines += [f"  {p}" for p in props]
            lines.append("}")
            strength = round(min(1.0, abs(float(pc["log2_fold_change"]))), 3)
            lines.append(
                f"rule {e['id']} {action} {ident(pc['gene'])} {{ strength: {strength}; "
                f"when: cell_type = {context(pc.get('tissue'))}; "
                f'evidence: predicted "AlphaGenome deletion, {_text(pc.get("tissue") or "strongest track")}" '
                f"{_text(_effect_note(pc))} }}"
            )
    from genomeos.attribution.measured import AGREES, DISAGREES, RECIPROCAL_OVERLAP, eligibility

    agree = sum(1 for r in measured_rows if r["agreement"]["verdict"] == AGREES)
    disagree = sum(1 for r in measured_rows if r["agreement"]["verdict"] == DISAGREES)
    negatives = sum(
        len(r["measured"].get("crispri", {}).get("genes_not_regulated", [])) for r in measured_rows
    )
    share = len(measured_rows) / len(elements) if elements else 0.0
    n_eligible = eligibility(elements, layer)["elements_in_an_assay_footprint"]
    # the section is written even when it is empty: a program that says nothing about the experimental
    # layer cannot be told apart from one whose assays were never read, and this project has been
    # bitten twice by a zero that came from not looking
    lines += [
        "",
        f"# ---- the experimental layer ({len(measured_rows)} of {len(elements)} elements, "
        f"{share * 100:.2f}%): what an assay measured over the same DNA",
        f"# The rule is reciprocal overlap >= {RECIPROCAL_OVERLAP}; nothing is raised on proximity or",
        f"# on an interval that merely contains the element. {agree} elements agree with the",
        f"# prediction, {disagree} disagree, and {negatives} element-gene pairs were measured as not",
        "# regulated - those are experimental facts about the absence of an effect, kept as such.",
    ]
    lines.append(
        f"# Eligible first, raised second: {n_eligible} of {len(elements)} elements lie in some "
        f"assay's footprint at all, {len(elements) - n_eligible} were never covered by one."
    )
    sat = [r["measured"]["satmut"] for r in measured_rows if "satmut" in r["measured"]]
    if sat:
        inert = sum(1 for s in sat if s["bases_measured"] and not s["bases_functional"])
        lines.append(
            f"# Base level: saturation mutagenesis measured {sum(s['bases_measured'] for s in sat)} "
            f"bases of {len(sat)} of these elements, {sum(s['bases_functional'] for s in sat)} of them "
            f"functional; {inert} elements had every measured base inert, which is a fact about those "
            "bases and is counted in neither the agreements nor the disagreements."
        )
    if not measured_rows:
        lines.append(
            "# No assay in data/knowledge measured any of these elements under the rule, so this "
            "program states no experimental fact."
        )
    for r in measured_rows:
        lines += _measured_blocks(chrom, r, domains, ident_of, sequence.get(r["id"]))
    return "\n".join(lines) + "\n"


def write_program(
    chrom: str, out: Path | None = None, results_dir: Path = RESULTS_DIR, layer: Any = None
) -> Path:
    out = out or Path("data/organisms/human") / f"noncoding_{chrom}.bio"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(compile_chromosome(chrom, results_dir, layer))
    return out
