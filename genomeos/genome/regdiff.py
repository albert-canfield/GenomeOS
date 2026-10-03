# SPDX-License-Identifier: AGPL-3.0-or-later
"""Genome diff at the regulatory level: two people's variants inside the elements that reach genes.

The coding diff says which proteins differ. Most of a person's variants are not in proteins; the ones
that matter next sit in the promoters and enhancers the registry knows, whose target the node model
infers and, where the deletion job has been, AlphaGenome predicts. This module reads a chromosome's
elements once (ENCODE cCREs, the CTCF nodes, the targets each element reaches, the cached deletion
prediction where one exists), places two people's normalised variants inside them, and writes the
difference the way the coding diff does: per target gene, the element and the variant only the first
person carries (+), only the second (-), with the element's class, its inferred target and basis, the
predicted target where cached, and, when asked, the variant's own phyloP constraint read by range.
Everything stays under the people's directories; the Markdown is printed or written where asked.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from genomeos.attribution.constraint import PHYLOP_THRESHOLD
from genomeos.genome.individuals import ROOT, _genotypes, _reference_base_reader, vcf_path

REFERENCE = "reference"
CLASS_RANK = {"promoter": 0, "enhancer": 1, "open_chromatin": 2, "insulator": 3, "unknown": 4}
CONSTRAINED = PHYLOP_THRESHOLD  # the attribution's own threshold for a constrained base


def chromosome_elements(chrom: str, reference: Path = Path("data/reference")):
    """The chromosome's regulatory elements with their targets, and an interval index over them."""
    from genomeos.attribution.eqtl import Intervals
    from genomeos.genome import Annotation, IndexedGenome, default_gencode, reference_fasta
    from genomeos.genome.domains import infer_domains
    from genomeos.genome.regulation import assign_targets
    from genomeos.genome.regulatory import load_ccres

    ccres = load_ccres(chrom)
    fasta = reference_fasta(chrom, reference)
    gff = default_gencode({chrom})
    if not ccres or not fasta.exists() or gff is None:
        raise FileNotFoundError(f"{chrom} is not fetched; run `genomeos data fetch --chrom {chrom}` first")
    genome = IndexedGenome(str(fasta))
    try:
        length = genome.lengths[chrom]
    finally:
        genome.close()
    ann = Annotation.from_gff3(gff, {chrom})
    domains = infer_domains(chrom, length, ccres, ann)
    elements = assign_targets(chrom, ccres, ann, domains)
    by_id = {e.id: e for e in elements}
    iv = Intervals()
    for e in elements:
        iv.add(chrom, e.locus.start, e.locus.end, e.id)
    return by_id, iv.freeze()


#: The element classes whose cached deletion prediction this module reads. The window below is read
#: for EXACTLY this population and no other, because a window formed where no head was asked for
#: would make the two counts compare different sets of variants.
PREDICTED_CLASSES = ("enhancer", "promoter")


def variants_in_elements(
    name: str, chrom: str, by_id, iv, root: Path, base_at=None, window: bool = False, coding=None
) -> dict[tuple, dict[str, Any]]:
    """The person's variants (normalised) that fall inside an element, keyed (pos, ref, alt).

    With `window=False`, which is the default and what every committed figure was produced under,
    nothing extra is opened and each row is exactly what it was: one `predicted` gene per element
    from `cached_prediction`. With `window=True` the SAME cached record `cached_prediction` already
    opened one call down is read again from the same per-process archive -- no new reader, no
    request -- and each row gains one new key, `window`, beside its unchanged `predicted`.
    """
    from genomeos.attribution.targets import prediction_window_reading
    from genomeos.predict.enhancer_target import cached_prediction, load_cached

    if name == REFERENCE:
        return {}
    vcf = vcf_path(name, chrom, root)
    if vcf is None:
        return {}
    out: dict[tuple, dict[str, Any]] = {}
    for (pos, ref, alt), zyg in _genotypes(vcf, base_at).items():
        ids = iv.at(chrom, pos)
        if not ids:
            continue
        e = by_id[ids[0]]
        target = e.targets[0] if e.targets else None
        asked = e.cls in PREDICTED_CLASSES
        pred = cached_prediction(chrom, e.id) if asked else None
        win = None
        if window and asked:
            w = prediction_window_reading(load_cached(chrom, e.id), (pred or {}).get("gene"), coding=coding)
            win = {
                "at_bar": [list(t) for t in (w.at_bar or ())],
                "extra": list(w.extra),
                "not_cached": w.not_cached,
                "head_agrees": w.head_agrees,
                "anything_at_bar": None if w.at_bar is None else bool(w.at_bar),
            }
        out[(pos, ref, alt)] = {
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "zygosity": "homozygous" if zyg == "hom" else "heterozygous",
            "element": e.id,
            "class": e.cls,
            "target": target["gene"] if target else None,
            "basis": target["basis"] if target else None,
            "distance": target["distance"] if target else None,
            "predicted": (
                {"gene": pred["gene"], "log2_fold_change": pred["log2_fold_change"], "tissue": pred["tissue"]}
                if pred
                else None
            ),
        }
        if win is not None:  # only with a reader, so the committed row shape is unchanged without one
            out[(pos, ref, alt)]["window"] = win
    return out


def _constrained(v: dict[str, Any]) -> bool:
    return (v.get("phylop") or 0) >= CONSTRAINED


def _rank(v: dict[str, Any]) -> tuple:
    return (
        0 if (v.get("phylop") or 0) >= CONSTRAINED else 1,
        0 if v.get("predicted") else 1,
        CLASS_RANK.get(v["class"], 5),
        v["zygosity"] != "homozygous",
        v["pos"],
    )


def add_constraint(chrom: str, rows: list[dict[str, Any]], progress=None) -> dict[str, Any]:
    """phyloP at each variant's base, one range pass over the sorted, distinct positions."""
    from genomeos.attribution.constraint import phylop_over_blocks

    positions = sorted({r["pos"] for r in rows})
    if not positions:
        return {}
    intervals = [(p - 1, p) for p in positions]
    stats, cost = phylop_over_blocks(chrom, intervals, threshold=CONSTRAINED, progress=progress)
    value = {p: (round(s.mean, 2) if s.bases else None) for p, s in zip(positions, stats, strict=False)}
    for r in rows:
        r["phylop"] = value.get(r["pos"])
    return cost


#: What `_window_arm` reports when it is asked for, and what each figure is for. Registered in
#: `attribution/onetarget2.py` for this module as: "either `with_predicted_target_only_a` or
#: `with_predicted_target_only_b` changing ... by THE_INVARIANT they must NOT move; a move means
#: `cached_prediction` and the window disagree about whether anything was predicted at all. The GENE
#: may change and that is the expected gain."
WINDOW_FALSIFIER = (
    "the two with_predicted_target counts are conditioned on `predicted` being truthy, which is the "
    "yes/no and not the gene, so the window must reproduce them exactly. The gene is where the gain "
    "is: `extra` names every gene the sweep put at the bar at that element and the one-target "
    "projection dropped."
)
#: The limit the invariant travels with, and whether it bites here. It does NOT: every figure below
#: is read off the ANY-GENE head, which is the head `cached_prediction` itself names, so the
#: 2026-09-27 coding-head exception (10 elements of 4,794 with no `predicted_coding` where the cache
#: named a coding gene) cannot reach this module's counts. A caller that passes `coding` gets the
#: coding window BESIDE the any-gene one and must read `coding_head_agrees` before using it.
WINDOW_DOES_NOT_READ_THE_CODING_HEAD = (
    "regdiff reads cached_prediction's own any-gene head, so the registered coding-head limit is not "
    "load-bearing for these counts; the coding window is formed only when a caller asks for it"
)
#: MEASURED, not assumed, and reported rather than used to drop the check: this module's registered
#: falsifier is SUBSUMED by `onetarget2.check_head_invariant`. An asked, cached row can take five
#: shapes -- head None with an empty bar, head G with G at the head of the bar, head None with a
#: non-empty bar, head G with an empty bar, head G with another gene at the head of the bar -- and
#: the first two leave the two counts equal while the last three are any-gene head disagreements,
#: which the invariant guard refuses before the counts are returned. So no input reachable through
#: `regulatory_diff` can set `the_count_moved`, and `tests/test_regdiff_window.py` records that by
#: driving `_window_arm` directly for the branch. The count is still computed and still reported,
#: because the registration names it as this module's own falsifier and the figure belongs in the
#: record; what is NOT claimed is that it adds a refusal the invariant guard does not already make.
WINDOW_FALSIFIER_IS_SUBSUMED_BY_THE_INVARIANT_GUARD = (
    "the two counts cannot diverge without an any-gene head disagreement, which check_head_invariant "
    "refuses first; this module's falsifier is therefore strictly weaker than that refusal and is "
    "reported as a figure rather than relied on as a second guard"
)


def _window_arm(window: bool, only_a: list[dict[str, Any]], only_b: list[dict[str, Any]]) -> dict[str, Any]:
    """The module's own falsifier, computed only when a window was read; `{}` otherwise.

    `{}` is why every committed figure is untouched without a reader: this contributes no key at all
    to `regulatory_diff`'s result unless it was asked for.
    """
    if not window:
        return {}
    from genomeos.attribution import onetarget2 as ot

    out: dict[str, Any] = {
        "falsifier": WINDOW_FALSIFIER,
        "limit": WINDOW_DOES_NOT_READ_THE_CODING_HEAD,
        "the_falsifier_is_subsumed": WINDOW_FALSIFIER_IS_SUBSUMED_BY_THE_INVARIANT_GUARD,
    }
    disagreements = 0
    for side, rows in (("only_a", only_a), ("only_b", only_b)):
        asked = [v for v in rows if "window" in v]
        head_yes = sum(1 for v in asked if v["predicted"])
        bar_yes = sum(1 for v in asked if v["window"]["anything_at_bar"])
        extra = [v for v in asked if v["window"]["extra"]]
        disagreements += sum(1 for v in asked if v["window"]["head_agrees"] is False)
        out[side] = {
            "variants_whose_element_was_asked": len(asked),
            "with_predicted_target_head": head_yes,
            "with_predicted_target_window": bar_yes,
            "the_count_moved": head_yes != bar_yes,
            "not_cached": sum(1 for v in asked if v["window"]["not_cached"]),
            "variants_gaining_a_gene": len(extra),
            "genes_gained": sum(len(v["window"]["extra"]) for v in asked),
            "genes_gained_at": [
                {
                    "element": v["element"],
                    "head": (v["predicted"] or {}).get("gene"),
                    "extra": v["window"]["extra"],
                }
                for v in extra[:20]
            ],
        }
    out["head_disagreements"] = disagreements
    ot.check_head_invariant(out, "genome/regdiff.py window arm")
    moved = out["only_a"]["the_count_moved"] or out["only_b"]["the_count_moved"]
    gained = out["only_a"]["variants_gaining_a_gene"] + out["only_b"]["variants_gaining_a_gene"]
    out["verdict_against_its_own_falsifier"] = (
        "REFUTES THE INVARIANT: a with_predicted_target count moved, so cached_prediction and the "
        "window disagree about whether anything was predicted. The wave stops."
        if moved
        else f"does NOT fire, as registered: both counts reproduced exactly; {gained} variant(s) "
        f"gained a dropped gene, which is the registered expected gain and not a count moving"
    )
    return out


def regulatory_diff(
    a: str,
    b: str,
    chrom: str,
    root: Path | None = None,
    constraint: bool = False,
    reference: Path = Path("data/reference"),
    progress=None,
    window: bool = False,
    coding=None,
) -> dict[str, Any]:
    root = root or ROOT
    by_id, iv = chromosome_elements(chrom, reference)
    base_at, genome = _reference_base_reader(chrom, reference)
    try:
        va = variants_in_elements(a, chrom, by_id, iv, root, base_at, window=window, coding=coding)
        vb = variants_in_elements(b, chrom, by_id, iv, root, base_at, window=window, coding=coding)
    finally:
        if genome is not None:
            genome.close()
    only_a = [v for k, v in va.items() if k not in vb]
    only_b = [v for k, v in vb.items() if k not in va]
    shared = sum(1 for k in va if k in vb)
    cost = add_constraint(chrom, only_a + only_b, progress) if constraint else None
    genes: dict[str, dict[str, Any]] = {}
    for side, rows in (("only_a", only_a), ("only_b", only_b)):
        for v in rows:
            g = genes.setdefault(
                v["target"] or "(no target in the node)", {"gene": v["target"], "only_a": [], "only_b": []}
            )
            g[side].append(v)
    for g in genes.values():
        g["only_a"].sort(key=_rank)
        g["only_b"].sort(key=_rank)
    ordered = sorted(genes.values(), key=lambda g: min(_rank(v) for v in g["only_a"] + g["only_b"]))

    def count(rows, pred) -> int:
        return sum(1 for v in rows if pred(v))

    by_class = {}
    for side, rows in (("only_a", only_a), ("only_b", only_b)):
        by_class[side] = {}
        for v in rows:
            by_class[side][v["class"]] = by_class[side].get(v["class"], 0) + 1
    return {
        "a": a,
        "b": b,
        "chrom": chrom,
        "elements": len(by_id),
        "variants_in_elements_a": len(va),
        "variants_in_elements_b": len(vb),
        "shared": shared,
        "only_a": len(only_a),
        "only_b": len(only_b),
        "by_class": by_class,
        "with_predicted_target_only_a": count(only_a, lambda v: v["predicted"]),
        "with_predicted_target_only_b": count(only_b, lambda v: v["predicted"]),
        **({"window": _window_arm(window, only_a, only_b)} if window else {}),
        "constrained_only_a": count(only_a, _constrained) if constraint else None,
        "constrained_only_b": count(only_b, _constrained) if constraint else None,
        "phylop_cost": cost,
        "genes": ordered,
        "evidence": "measured genotypes; curated: ENCODE cCRE classes; inferred: the CTCF node's nearest TSS "
        "as target; predicted: AlphaGenome deletion target where the element was scored; measured: Zoonomia "
        "phyloP at the base when asked",
        "note": "a variant inside an element is a candidate, not an effect: the element's activity in the "
        "cell, the base's role in it and the direction of the change are not read here",
    }


def _line(sign: str, v: dict[str, Any]) -> str:
    tags = []
    if v.get("phylop") is not None:
        tags.append(f"phyloP {v['phylop']:+.1f}" + (" constrained" if v["phylop"] >= CONSTRAINED else ""))
    if v.get("predicted"):
        p = v["predicted"]
        tags.append(f"deletion moves {p['gene']} {p['log2_fold_change']:+.2f} in {p['tissue'][:20]}")
    where = f"{v['class']} {v['element']}"
    if v.get("target"):
        where += f" → {v['target']} ({v['basis']}, {v['distance']:,} bp)"
    return f"{sign} {v['pos']:,} {v['ref']}>{v['alt']} ({v['zygosity']}) in {where}" + (
        "  [" + "; ".join(tags) + "]" if tags else ""
    )


def render(d: dict[str, Any], top: int = 40) -> str:
    a, b = d["a"], d["b"]
    against = "the reference" if b == REFERENCE else b
    cons = ""
    if d["constrained_only_a"] is not None:
        cons = f" Constrained bases (phyloP ≥ {CONSTRAINED}): {d['constrained_only_a']} only {a}" + (
            f", {d['constrained_only_b']} only {b}." if b != REFERENCE else "."
        )
    lines = [
        f"# {a} against {against}, {d['chrom']}: variants inside regulatory elements",
        "",
        f"{d['variants_in_elements_a']:,} of {a}'s variants sit inside one of the chromosome's "
        f"{d['elements']:,} elements"
        + (f", {d['variants_in_elements_b']:,} of {b}'s, {d['shared']:,} shared" if b != REFERENCE else "")
        + f". Only {a}: {d['only_a']:,} ({d['with_predicted_target_only_a']} in an element the deletion "
        f"job has scored)"
        + (
            f"; only {b}: {d['only_b']:,} ({d['with_predicted_target_only_b']} scored)."
            if b != REFERENCE
            else "."
        )
        + cons,
        "",
        f"_{d['evidence']}._ {d['note']}",
        "",
    ]
    for g in d["genes"][:top]:
        lines.append(f"## {g['gene'] or 'no coding target in the node'}")
        lines.append("")
        lines.append("```diff")
        for v in g["only_a"]:
            lines.append(_line("+", v))
        for v in g["only_b"]:
            lines.append(_line("-", v))
        lines.append("```")
        lines.append("")
    if len(d["genes"]) > top:
        lines.append(f"… and {len(d['genes']) - top:,} more target genes.")
    return "\n".join(lines)
