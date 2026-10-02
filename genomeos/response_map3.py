# SPDX-License-Identifier: AGPL-3.0-or-later
"""Increment 3 of the response map: one chain read through several kinds of evidence.

Increment 2 built a chain over 110 loci and stated its own limit, carried here word for word: the
evidenced chain runs over 110 loci, and that is "23.26% of one of the measured layer's four assays";
lentiMPRA, VISTA and saturation mutagenesis contributed NOTHING to it, because "they read a sequence,
or the bases inside it, not the native locus after a perturbation of it". So increment 2 is one assay
deep. This increment reads a chain through more than one kind of evidence, over the loci the count at
`data/results/response_map_increment3_count.json` found carry one.

**What is added, and what is not.** Nothing about increment 2's three assertions changes: they are
built by `response_map2`'s own makers, imported, not reimplemented. Two kinds of observation are added
beside them:

- *a lentiMPRA tile reading* - one tile of the ENCODE4 joint library in one cell, with its own
  log2(RNA/DNA). The observation is of a **200 bp copy driving an integrated lentiviral reporter,
  outside its native locus**, so it is an observation of the element's SEQUENCE and never of the
  native locus after a perturbation of it. The per-cell active/silent label is the measured layer's
  declared aggregation, not the observation, and is carried as such with its rule word for word.
- *a GTEx cis-eQTL* - one retained significant variant-gene pair in one tissue. The variant sits on
  the chromosome in place and the expression is the donor's own, so this one does read the native
  locus; but it is natural variation and an association, never an experimental perturbation. Its
  relation therefore stays unresolved between the two registered alternatives.

**No verdict, no score, no new cut-off.** `response_map.check` is the gate, as in increment 2. Several
kinds of evidence at one locus are not agreement, confirmation or validation of one another, and the
map says so in its own section. Every rule, span, threshold and call is imported from where it was
registered: the overlap rule and the assay labels from `attribution.measured`, the openness call from
`attribution.context_evidence`, the locus grouping from `attribution.cell2`, and the eQTL window from
the margin the retained set was distilled under.

**Reconciliation.** `loci_built + sum(loci_excluded_by_cause)` must equal the population the committed
count named, and `build` refuses to return a payload where it does not.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

from genomeos import response_map as rm
from genomeos import response_map2 as rm2
from genomeos.attribution import cell2, eqtl, mpra
from genomeos.attribution import context_evidence as ce
from genomeos.attribution import measured as ms
from genomeos.genome import reader
from genomeos.results import RESULTS_DIR

INCREMENT = 3
VIEW = rm.VIEW
RESULT = "response_map_increment3"
COUNT_RESULT = "data/results/response_map_increment3_count.json"
CHROMS = rm2.CHROMS
#: Fewer loci a page than increment 2's ten, because a locus here carries more assertions. The whole
#: map's counts sit beside every page all the same.
PAGE = 8

HEADER = (
    "A read-only map of existing evidence and of where it stops, for the loci of the measured "
    "perturbation layer that carry a further kind of evidence beside the perturbation and the "
    "reader's own reading. It adds no biological validation, no prediction and no verdict, and a "
    "second kind of evidence at a locus is not agreement with the first."
)
TITLE = "The loci read through several kinds of evidence"

#: Increment 2's stated limit, word for word, and the one this increment starts from.
INCREMENT_2_LIMIT = (
    "the evidenced chain runs over 110 loci, and that is \"23.26% of one of the measured layer's four "
    'assays". lentiMPRA, VISTA and saturation mutagenesis contributed NOTHING to it, because "they '
    'read a sequence, or the bases inside it, not the native locus after a perturbation of it"'
)

#: Increment 2's four causes, kept as they are, and the three this increment can refuse under.
CAUSE_REPORTER_SOURCE = "reporter_source_file_not_in_this_checkout"
CAUSE_EQTL_SOURCE = "eqtl_hits_file_not_in_this_checkout"
CAUSE_NO_FURTHER_KIND = "no_further_kind_of_evidence_survived_with_a_full_source"
CAUSES = (*rm2.CAUSES, CAUSE_REPORTER_SOURCE, CAUSE_EQTL_SOURCE, CAUSE_NO_FURTHER_KIND)

#: Why a reporter reading is not an observation of the locus, stated on every one of them.
REPORTER_NOT_THE_LOCUS = (
    "this reading is of a 200 bp copy of the sequence driving an integrated lentiviral reporter, "
    "outside its native locus. It is an observation of the SEQUENCE and not of the native locus after "
    "a perturbation of it; it neither supports nor contradicts the measured perturbation beside it, "
    "and the two are not two measurements of one thing"
)
#: Why a silent reporter beside a measured decrease is not a conflict, carried on the locus.
SILENT_IS_NOT_A_CONTRADICTION = (
    "a reporter calling the sequence silent and a CRISPRi screen measuring a decrease at the native "
    "locus are not in conflict: the reporter read a copy outside the locus and the screen perturbed "
    "the locus in place. Neither checks the other, and no agreement or disagreement is recorded"
)
#: The reporter's own relation. The object is a reporter-activity state, never the native gene.
REPORTER_RELATION_NOTE = (
    "the tile is associated with a measured reporter-activity state of the construct it was read in. "
    "The state is of the reporter, not of the native locus, and no relation to any gene is asserted"
)
#: Why an eQTL's relation is left unresolved between the two registered alternatives.
EQTL_RELATION_NOTE = (
    "the retained record says this variant's alleles go with the gene's measured RNA across donors in "
    "this tissue. Whether the interval the variant sits in changes that RNA, or merely lies beside "
    "something that does, is not established by the record, so the relation stays unresolved between "
    "the two registered alternatives"
)
EQTL_UNCERTAINTY = (
    "an association over donors, not a perturbation: nothing was done to the locus",
    "linkage carries the signal some distance, so the variant inside the interval need not be the causal one",
    "the eGene need not be regulated by this interval",
    "GTEx tissues are not the CRISPRi cell lines, so this reading and the perturbation beside it are "
    "observations in different contexts",
)
#: What the retained eQTL set is, and what its absence at an interval does and does not mean.
EQTL_FRAME = (
    "the retained hits were distilled against the intervals of three sampled element sets, so a "
    "variant outside that frame was discarded before this map existed. Where no eQTL is shown at an "
    "interval inside the frame the retained set holds none; where the interval is outside the frame "
    "nothing could have been seen at all, and the absence is UNASSESSED and not an absence of eQTLs"
)
#: The map's own statement that more kinds is not more established.
MORE_KINDS_IS_NOT_MORE_ESTABLISHED = (
    "a locus carrying three or four kinds of evidence is not a locus three or four times established. "
    "The kinds measure different things in different contexts - a perturbation of the locus, the "
    "locus's accessibility, a copy of its sequence in a reporter, an association over donor tissues - "
    "and none of them is a replication or a check of another"
)

CHAIN_STOPS = rm2.CHAIN_STOPS
CHAIN_NOTE = rm.CHAIN_NOTE


# --- the further kinds of evidence -------------------------------------------------------------------
def _reporter_state(m: rm._Map, cell: str, label: str) -> str:
    """The measured reporter-activity state a tile reading is an association with, as its own entity."""
    return m.entity(
        f"state:reporter_activity:{label}:{cell}",
        "cell_state",
        state=f"reporter_activity_{label}",
        biosample=cell,
        reading=(
            "the share of this element's matched tiles at or above the reporter's own activity "
            f"threshold in {cell}, by the measured layer's declared aggregation"
        ),
        assay=mpra.EVIDENCE,
        of_a_reporter_construct=True,
        not_the_native_locus=REPORTER_NOT_THE_LOCUS,
        is_validation=False,
    )


def reporter_observations(
    c: rm2.Candidate,
    m: rm._Map,
    files: rm2._Files,
    element_iv: str,
    basis: str,
    knowledge: Path = mpra.KNOWLEDGE,
) -> tuple[list[dict[str, Any]], str | None]:
    """One assertion per (tile, cell), the measured layer's own reporter unit, or a refusal cause.

    The observation is the tile's log2(RNA/DNA) in one cell, which is what was measured. The per-cell
    active/silent label is the declared aggregation over tiles and is carried beside it with its rule,
    never as the observation.
    """
    row = c.row or {}
    blk = (row.get("measured") or {}).get("lentimpra") or {}
    tiles = blk.get("tiles") or []
    labels = blk.get("label_by_cell") or {}
    if not tiles:
        return [], None
    out: list[dict[str, Any]] = []
    for tile in tiles:
        for cell, value in sorted((tile.get("activity") or {}).items()):
            acc = mpra.FILES.get(cell)
            if acc is None:
                continue
            path = (knowledge / f"{acc}.bed.gz").as_posix()
            key = f"{tile['name']} {c.pair.chrom}:{tile['start']}-{tile['end']} {cell}"
            ref = files.ref(
                path,
                key,
                ["bed", c.pair.chrom, [tile["start"], tile["end"]], tile["name"]],
                accession=acc,
                cell_of_the_file=cell,
            )
            if ref is None:
                return [], CAUSE_REPORTER_SOURCE
            label = labels.get(cell, rm.NOT_RECORDED)
            tile_iv = m.interval(
                c.pair.chrom,
                int(tile["start"]),
                int(tile["end"]),
                basis,
                reporter_tile=tile["name"],
                read_as_a_copy_outside_its_locus=True,
            )
            out.append(
                {
                    "id": (f"m3|{c.pair.chrom}:{tile['start']}-{tile['end']}|{tile['name']}|{cell}"),
                    "relation": "associated_with_state",
                    "relation_note": REPORTER_RELATION_NOTE,
                    "subject": tile_iv,
                    "object": _reporter_state(m, cell, str(label)),
                    "entities": [element_iv],
                    "status": "observed",
                    "basis": ms.OUTCOME_KIND["lentimpra"],
                    "context": {
                        "cell": cell,
                        "condition": "an integrated lentiviral reporter construct",
                        "perturbation": (
                            "none of the native locus: a 200 bp copy of the sequence was placed in a reporter"
                        ),
                    },
                    "measurement": {
                        "quantity": "log2(RNA/DNA) of this tile in this cell",
                        "value": value,
                        "unit": "log2 ratio",
                        "observation_unit": list(ms.REPORTER_UNIT),
                        "strand": (tile.get("strand") or {}).get(cell, rm.NOT_RECORDED),
                        "overlap_with_the_compiled_element": tile.get("overlap"),
                        "uncertainty": ms.REPORTER_UNCERTAINTY,
                    },
                    "aggregated_label": {
                        "cell": cell,
                        "label": label,
                        "rule": ms.REPORTER_LABEL_RULE,
                        "threshold_log2": ms.MPRA_ACTIVE,
                        "is_the_observation": False,
                        "note": (
                            "the label is the measured layer's declared aggregation over this "
                            "element's matched tiles, not a reading of this tile"
                        ),
                    },
                    "of_the_native_locus": False,
                    "not_the_native_locus": REPORTER_NOT_THE_LOCUS,
                    "is_validation": False,
                    "mechanism": "not established: reporter activity names no mechanism",
                    "uncertainty": [REPORTER_NOT_THE_LOCUS, ms.REPORTER_UNCERTAINTY],
                    "source": ref,
                    "assembly_basis": basis,
                    "observation_key": f"lentimpra|{key}",
                }
            )
    return out, None


def hits_for(chrom: str, knowledge: Path = eqtl.KNOWLEDGE) -> list[dict[str, Any]]:
    """Every retained GTEx hit of one chromosome, with the file it came from and its own row.

    The columns are taken from `eqtl.HIT_COLUMNS`, where the distillation registered them, so the
    reading of the file is not this module's own guess at its shape.
    """
    cols = list(eqtl.HIT_COLUMNS)
    out: list[dict[str, Any]] = []
    for p in sorted(knowledge.glob("hits_*.tsv")):
        with p.open() as fh:
            header = (fh.readline() or "").rstrip("\n").split("\t")
            if header != cols:
                raise rm.RefusedError(f"{p}: columns {header} are not {cols}")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[1] != chrom:
                    continue
                row = dict(zip(cols, f, strict=True))
                row["pos"] = int(row["pos"])
                row["slope"] = float(row["slope"])
                row["pval_nominal"] = float(row["pval_nominal"])
                row["file"] = p.as_posix()
                out.append(row)
    out.sort(key=lambda r: (r["pos"], r["gene_id"], r["tissue"]))
    return out


def eqtl_observations(
    c: rm2.Candidate,
    m: rm._Map,
    files: rm2._Files,
    element_iv: str,
    basis: str,
    hits: list[dict[str, Any]],
    margin: int,
) -> tuple[list[dict[str, Any]], str | None]:
    """One assertion per retained variant-gene pair inside the element, or a refusal cause.

    One record, one assertion: no aggregation over tissues or variants, so no statistic is introduced
    here and every assertion names exactly one row of one file.
    """
    e = c.element
    assert e is not None
    lo, hi = e["start"] - margin, e["end"] + margin
    out: list[dict[str, Any]] = []
    for h in hits:
        if not lo < h["pos"] <= hi:
            continue
        key = f"{h['tissue']} {h['chrom']}:{h['pos']} {h['ref']}>{h['alt']} {h['gene_id']}"
        ref = files.ref(h["file"], key, ["hits", h["tissue"], h["pos"], h["gene_id"]], tissue=h["tissue"])
        if ref is None:
            return [], CAUSE_EQTL_SOURCE
        variant = m.interval(
            h["chrom"],
            h["pos"] - 1,
            h["pos"],
            basis,
            variant=f"{h['chrom']}:{h['pos']}:{h['ref']}>{h['alt']}",
            inside_element=e["id"],
        )
        out.append(
            {
                "id": (f"q3|{h['chrom']}:{h['pos']}|{h['ref']}>{h['alt']}|{h['gene_id']}|{h['tissue']}"),
                "relation": rm.UNRESOLVED,
                "relation_candidates": ["changes_measured_rna", "associated_with_state"],
                "relation_note": EQTL_RELATION_NOTE,
                "subject": variant,
                "object": m.gene(h["gene_id"]),
                "entities": [element_iv],
                "status": "observed",
                "basis": eqtl.EVIDENCE,
                "context": {
                    "cell": f"GTEx tissue: {h['tissue']}",
                    "condition": "natural variation across donors",
                    "perturbation": "none: nothing was done to the locus",
                },
                "measurement": {
                    "quantity": "slope of the gene's expression on the alternate allele (GTEx)",
                    "value": round(h["slope"], 6),
                    "unit": "normalised expression per allele",
                    "direction": "increase" if h["slope"] > 0 else "decrease",
                    "p_nominal": h["pval_nominal"],
                    "significance_call": (
                        "GTEx v8's own significant variant-gene pair set for this tissue; no threshold "
                        "is applied here"
                    ),
                    "is_a_perturbation": False,
                },
                "of_the_native_locus": True,
                "an_association_not_a_perturbation": True,
                "is_validation": False,
                "mechanism": "not established: an association names no mechanism",
                "uncertainty": list(EQTL_UNCERTAINTY),
                "source": ref,
                "assembly_basis": basis,
                "observation_key": f"eqtl|{key}",
            }
        )
    return out, None


# --- the build ---------------------------------------------------------------------------------------
def counted_population(root: Path = Path()) -> dict[str, Any]:
    """The population this increment is built over, read from the committed count, with its digest.

    The number is not restated here. It is read out of the count's own
    `kinds_per_locus.loci_carrying_at_least_one_further_kind`, and the reconciliation below is against
    that figure, so the build cannot quietly be over a different set than the one that was counted.
    """
    path = Path(root) / COUNT_RESULT
    files = rm2._Files(Path(root))
    digest = files.hash(COUNT_RESULT)
    if not path.is_file():
        raise rm.RefusedError(
            f"{COUNT_RESULT} is not in this checkout: the population was counted before the build and "
            "the build is reconciled against that count, so it cannot run without it"
        )
    with path.open() as fh:
        d = json.load(fh)
    block = d["kinds_per_locus"]
    return {
        "result": d["result"],
        "path": COUNT_RESULT,
        "sha256": digest["sha256"],
        "addable_loci": d["reconciliation"]["addable_loci"],
        "loci_carrying_at_least_one_further_kind": block["loci_carrying_at_least_one_further_kind"],
        "beyond_those_two": block["beyond_those_two"],
        "question_it_answered": d["question"],
        "increment_2_stated_limit": d["starting_point"]["increment_2_stated_limit"],
    }


def build(  # noqa: PLR0915 - one pass over the loci, kept in the order the payload is assembled in
    root: Path = Path(),
    chroms: tuple[str, ...] | list[str] = CHROMS,
    results_dir: Path = RESULTS_DIR,
    readers: ce.Readers | None = None,
    table: ce.Mapping | None = None,
    margin: int = 500,
) -> dict[str, Any]:
    """The map of every locus carrying a further kind of evidence, reconciled against the count."""
    chroms = list(chroms)
    root = Path(root)
    readers = readers if readers is not None else ce.Readers()
    table = table if table is not None else ce.mapping()
    population = counted_population(root)
    cands, counted = rm2.candidates(chroms, results_dir, readers, table)

    covered = rm2.map1_keys(root)
    groups = cell2.group([c.key for c in cands] + covered)
    mine, theirs = groups[: len(cands)], set(groups[len(cands) :])
    best: dict[int, str] = {}
    for g, c in zip(mine, cands, strict=True):
        if g not in best or rm2.STATUS_ORDER[c.status] < rm2.STATUS_ORDER[best[g]]:
            best[g] = c.status
    addable = {g for g, s in best.items() if s == rm2.STATUS_ADDABLE}
    members: dict[int, list[int]] = {}
    for i, g in enumerate(mine):
        members.setdefault(g, []).append(i)

    hits_of: dict[str, list[dict[str, Any]]] = {}

    def hits(chrom: str) -> list[dict[str, Any]]:
        if chrom not in hits_of:
            hits_of[chrom] = hits_for(chrom)
        return hits_of[chrom]

    # which of the addable loci carry a further kind at all: the population, re-derived
    def carries_further(g: int) -> bool:
        for i in members[g]:
            c = cands[i]
            if c.status != rm2.STATUS_ADDABLE or c.element is None:
                continue
            if "lentimpra" in ((c.row or {}).get("assays") or ()):
                return True
            lo, hi = c.element["start"] - margin, c.element["end"] + margin
            if any(lo < h["pos"] <= hi for h in hits(c.pair.chrom)):
                return True
        return False

    population_groups = {g for g in sorted(addable) if carries_further(g)}
    refuse_unless_population_matches(
        len(population_groups), population["loci_carrying_at_least_one_further_kind"]
    )

    files = rm2._Files(root)
    assembly = rm._assembly(rm._Sources(root))
    basis = ", ".join(assembly["from"])
    m = rm._Map(assembly)
    loci: list[dict[str, Any]] = []
    chains: list[dict[str, Any]] = []
    excluded: dict[str, list[dict[str, Any]]] = {c: [] for c in CAUSES}
    pairs_excluded: dict[str, int] = dict.fromkeys(CAUSES, 0)
    seen_keys: set[str] = set()
    added: set[str] = set()
    state_counts: Counter[str] = Counter()
    not_open: list[dict[str, Any]] = []
    kinds_per_locus: Counter[str] = Counter()

    for g in sorted(population_groups):
        idx = members[g]
        keys = [cands[i].key for i in idx]
        lid = rm2._locus_id([cands[i].key for i in idx if cands[i].status == rm2.STATUS_ADDABLE])
        built: list[dict[str, Any]] = []
        causes: list[str] = []
        for i in idx:
            c = cands[i]
            if c.status != rm2.STATUS_ADDABLE:
                continue
            assert c.element is not None and c.row is not None
            element_iv = rm._iv(c.pair.chrom, c.element["start"], c.element["end"])
            pert = rm2._perturbation(c, m, files, element_iv, basis)
            rule = rm2._compiled_rule(c, m, files, c.pair.chrom, results_dir, basis)
            state = rm2._reader_state(c, m, files, c.pair.chrom, readers, element_iv, basis)
            reporter, rcause = reporter_observations(c, m, files, element_iv, basis)
            eq, qcause = eqtl_observations(c, m, files, element_iv, basis, hits(c.pair.chrom), margin)
            if pert is None:
                cause: str | None = rm2.CAUSE_PAIR_SOURCE
            elif rule is None:
                cause = rm2.CAUSE_ELEMENT_SOURCE
            elif state is None:
                cause = rm2.CAUSE_READER_SOURCE
            elif rcause is not None:
                cause = rcause
            elif qcause is not None:
                cause = qcause
            elif pert["observation_key"] in seen_keys:
                cause = rm2.CAUSE_DUPLICATE_KEY
            else:
                cause = None
            if cause is not None:
                causes.append(cause)
                pairs_excluded[cause] += 1
                continue
            seen_keys.add(pert["observation_key"])
            built.append(
                {
                    "candidate": c,
                    "pert": pert,
                    "rule": rule,
                    "state": state,
                    "reporter": reporter,
                    "eqtl": eq,
                }
            )
        further = [b for b in built if b["reporter"] or b["eqtl"]]
        if not built or not further:
            cause = causes[0] if causes else CAUSE_NO_FURTHER_KIND
            excluded[cause].append(
                {
                    "locus": lid,
                    "chrom": keys[0].chrom,
                    "genes": sorted({k.gene for k in keys}),
                    "causes_of_its_pairs": sorted(set(causes)) or [CAUSE_NO_FURTHER_KIND],
                }
            )
            continue
        ids: list[str] = []
        for b in built:
            for a in (b["pert"], b["rule"], b["state"], *b["reporter"], *b["eqtl"]):
                if a["id"] not in added:
                    added.add(a["id"])
                    m.add(a)
            c = b["candidate"]
            state_counts[c.state] += 1
            if c.state == ce.STATE_NOT_OPEN:
                not_open.append(
                    {
                        "locus": lid,
                        "assertion": b["state"]["id"],
                        "element": c.element["id"] if c.element else None,
                        "cell": c.pair.cell,
                        "biosample": c.biosample,
                        "reading": ce.READING[ce.STATE_NOT_OPEN],
                        "not_closed": ce.NOT_CLOSED,
                    }
                )
            steps = [b["pert"]["id"], b["rule"]["id"], b["state"]["id"]]
            beside = [a["id"] for a in (*b["reporter"], *b["eqtl"])]
            ids += steps + beside
            chains.append(
                {
                    "id": f"chain3|{b['pert']['id']}",
                    "title": (
                        f"{c.pair.chrom}:{c.pair.start}-{c.pair.end} -> {c.element['id']} -> "
                        f"{c.pair.gene} in {c.pair.cell}"
                    ),
                    "steps": steps,
                    "beside": beside,
                    "beside_note": (
                        "read beside the chain, not as a step in it: a reporter reading is of a copy "
                        "of the sequence and an eQTL is an association over donors, so neither "
                        "continues the perturbation's chain nor checks it"
                    ),
                    "kinds_of_evidence_at_this_chain": _kinds_of(b),
                    "derived_claims": [],
                    "stops_at": list(CHAIN_STOPS),
                    "note": CHAIN_NOTE,
                }
            )
        shown_keys = [cands[i].key for i in idx if cands[i].status == rm2.STATUS_ADDABLE]
        tag = "+".join(
            sorted({k for b in built for k in _kinds_of(b)} - {"crispri_perturbation", "reader_dnase"})
        )
        kinds_per_locus[tag or "nothing_beyond_the_two"] += 1
        loci.append(
            {
                "id": lid,
                "chrom": keys[0].chrom,
                "chromosomes": sorted({k.chrom for k in keys}),
                "start": min(k.start for k in shown_keys),
                "end": max(k.end for k in shown_keys),
                "span_over_every_pair_grouped_here": [
                    min(k.start for k in keys),
                    max(k.end for k in keys),
                ],
                "span_note": (
                    "`start` and `end` are the span of the shown pairs. The 1 Mb grouping chains, so "
                    "the group can reach much further than the pairs the map shows at it; "
                    "`span_over_every_pair_grouped_here` is that wider interval, and neither is "
                    "established biological independence"
                ),
                "cells": sorted({k.cell for k in shown_keys}),
                "cells_of_every_pair_grouped_here": sorted({k.cell for k in keys}),
                "genes_measured": sorted({k.gene for k in shown_keys}),
                "elements": sorted({b["candidate"].element["id"] for b in built}),
                "assertions": ids,
                "chains": [f"chain3|{b['pert']['id']}" for b in built],
                "reader_states": sorted({b["candidate"].state for b in built}),
                "kinds_of_evidence": sorted({k for b in built for k in _kinds_of(b)}),
                "kinds_beyond_the_two": tag.split("+") if tag else [],
                "more_kinds_is_not_more_established": MORE_KINDS_IS_NOT_MORE_ESTABLISHED,
                "reporter_readings": sum(len(b["reporter"]) for b in built),
                "eqtl_records": sum(len(b["eqtl"]) for b in built),
                "silent_is_not_a_contradiction": SILENT_IS_NOT_A_CONTRADICTION,
                "pairs_shown": len(built),
                "pairs_in_this_locus_not_shown": {
                    s: sum(1 for i in idx if cands[i].status == s) for s in rm2.STATUSES[1:]
                },
                "already_in_increment_1": g in theirs,
                "already_in_increment_2": True,
            }
        )

    payload = rm.assemble(
        m.entities,
        m.assertions,
        [],
        chains,
        view=VIEW,
        increment=INCREMENT,
        example=RESULT,
        title=TITLE,
        header=HEADER,
        why=_why(chroms, population),
        conventions=_conventions(assembly, margin),
        vocabulary={
            "entity_kinds": list(rm.ENTITY_KINDS),
            "relations": list(rm.RELATIONS),
            "statuses": list(rm.STATUSES),
            "categories": rm.CATEGORIES,
        },
        loci=loci,
        page_size=PAGE,
        counted=counted,
        population=population,
        kinds_of_evidence=_kinds_block(payload_loci=loci, kinds_per_locus=kinds_per_locus),
        not_open_in_reader=rm2._not_open_block(not_open, state_counts),
        cannot_establish=list(CANNOT_ESTABLISH),
        requests={"network": 0, "model": 0, "downloads": 0},
    )
    payload["reconciliation"] = _reconcile(loci, population, excluded, pairs_excluded, counted)
    refuse_unless_reconciled(payload["reconciliation"])
    payload["counts"]["loci"] = len(loci)
    payload["counts"]["observed_by_assay"] = observed_by_assay(payload["assertions"])
    payload["outcomes_by_assay"] = outcomes_by_assay(payload["assertions"], loci)
    payload["sources"] = [files.table[p] for p in sorted(files.table)]
    return payload


def refuse_unless_population_matches(derived: int, counted: int) -> None:
    """Refuse a build over a set of a different size than the one that was counted before it.

    The count was committed first on purpose. If the population re-derived at build time is not the
    one that was counted, the build's own reconciliation would be against a denominator nothing on the
    record carries, so it does not run at all.
    """
    if derived != counted:
        raise rm.RefusedError(
            f"the population re-derived here holds {derived} loci and the committed count named "
            f"{counted}: the build would not be over the set that was counted"
        )


def refuse_unless_reconciled(reconciliation: dict[str, Any]) -> None:
    """Refuse to serve a payload whose built and excluded loci do not sum to the counted population."""
    if not reconciliation["reconciles"]:
        raise rm.RefusedError(
            f"refusing to serve this payload: {reconciliation['loci_built']} loci built plus "
            f"{sum(reconciliation['loci_excluded_by_cause'].values())} excluded is "
            f"{reconciliation['loci_accounted_for']}, not the "
            f"{reconciliation['loci_counted_before_the_build']} the count named"
        )


def _kinds_of(b: dict[str, Any]) -> list[str]:
    out = ["crispri_perturbation", "reader_dnase"]
    if b["reporter"]:
        out.append("lentimpra_sequence")
    if b["eqtl"]:
        out.append("eqtl_native_locus")
    return out


def observed_by_assay(assertions: list[dict[str, Any]]) -> dict[str, int]:
    """`observed` split by what measured it. Four assays are not one experiment."""
    out: Counter[str] = Counter()
    for a in assertions:
        if a["status"] != "observed":
            continue
        for p, label in ASSAY_OF_PREFIX.items():
            if a["id"].startswith(p):
                out[label] += 1
                break
        else:  # pragma: no cover - a new observed kind must name what measured it
            raise rm.RefusedError(f"{a['id']}: observed, but nothing says what measured it")
    return dict(sorted(out.items()))


def _outcome_of(a: dict[str, Any]) -> str:
    """What one assertion's own record says it found, read out of the field that holds it.

    Nothing is thresholded, compared or recomputed here: the CRISPRi outcome is the benchmark's own
    call as `measured` recorded it, the reader state is reader v1's own call, the reporter label is the
    measured layer's declared per-cell aggregation carried on the tile, the eQTL direction is the sign
    of GTEx's own slope, and the predicted arm's direction is the compiler's. A kind whose outcome
    field this does not know is refused rather than pooled into a total that would hide it.
    """
    m = a["measurement"]
    if a["id"].startswith("p2|"):
        return str(m["outcome"])
    if a["id"].startswith("s2|"):
        return str(m["value"])
    if a["id"].startswith("m3|"):
        return str(a["aggregated_label"]["label"])
    if a["id"].startswith(("q3|", "r2|")):
        return str(m["direction"])
    raise rm.RefusedError(f"{a['id']}: nothing on this assertion says what it found")


#: Which assay each assertion id's prefix came from. `r2|` is the compiled rule and is deliberately
#: absent: it is predicted and is never reported beside the measured ones.
ASSAY_OF_PREFIX = {
    "p2|": "CRISPRi perturbation (ENCODE benchmark)",
    "s2|": reader.EVIDENCE,
    "m3|": mpra.EVIDENCE,
    "q3|": eqtl.EVIDENCE,
}

A_COUNT_IS_NOT_A_MEASUREMENT = (
    "a count of observations is not a measurement of what they say, and the two have been confused on "
    "this project before: 198 loci where the model says nothing were quoted as 198 places it is wrong "
    "when 66 of them held a measured decrease, and 48 measured links cleared a floor while only 21 "
    "were answerable. So every count on this map carries what its own observations found, at its own "
    "denominator"
)

OUTCOMES_ARE_NOT_AGREEMENT = (
    "these breakdowns stand side by side and are never crossed. Each is at its own denominator, in its "
    "own context, from its own assay, and no row here is read against another: nothing on this map "
    "counts a locus whose readings point the same way, because two of these kinds measuring different "
    "things in different contexts cannot agree or disagree in the first place"
)


def outcomes_by_assay(assertions: list[dict[str, Any]], payload_loci: list[dict[str, Any]]) -> dict[str, Any]:
    """What each assay's observations on this map actually say, each at its own denominator.

    The map's `counts` block says how many observations each assay contributed. That is a fact about
    where the assays were pointed and says nothing about what any of them found, which is the error
    this block exists to prevent. Every outcome is read off the assertion that carries it, by
    `_outcome_of`, so no rule, threshold or call is applied here that was not already registered.

    Two denominators are reported for each assay, because they answer different questions and have
    been mixed up before: `observations`, one per assertion, and `loci_with_at_least_one_observation`,
    one per locus of the map. A locus can hold observations of more than one outcome, so the locus rows
    are counted independently and need not sum to the map's locus count; the block says so of itself.
    """
    by_locus: dict[str, set[str]] = {}
    for locus in payload_loci:
        for aid in locus["assertions"]:
            by_locus.setdefault(aid, set()).add(locus["id"])

    blocks: dict[str, Any] = {}
    for prefix, label in ASSAY_OF_PREFIX.items():
        mine = [a for a in assertions if a["id"].startswith(prefix)]
        obs: Counter[str] = Counter()
        loci_of: dict[str, set[str]] = {}
        for a in mine:
            o = _outcome_of(a)
            obs[o] += 1
            loci_of.setdefault(o, set()).update(by_locus.get(a["id"], set()))
        blocks[label] = {
            "status": "observed",
            "observations": len(mine),
            "what_they_found": dict(sorted(obs.items())),
            "loci_with_at_least_one_observation": {o: len(ids) for o, ids in sorted(loci_of.items())},
            "loci_rows_do_not_sum": (
                "a locus can carry observations of more than one outcome, so these rows are counted "
                "independently and do not sum to the map's locus count"
            ),
            "outcome_field": OUTCOME_FIELD[prefix],
        }

    pred = [a for a in assertions if a["status"] == "predicted"]
    predicted: Counter[str] = Counter()
    bases = sorted({a["basis"] for a in pred})
    for a in pred:
        predicted[_outcome_of(a)] += 1
    return {
        "measured": blocks,
        "predicted": {
            "status": "predicted",
            "is_a_measurement": False,
            "assertions": len(pred),
            "what_the_compiler_said": dict(sorted(predicted.items())),
            "basis": bases,
            "never_observed": (
                "a predicted quantity is never `observed` and is never counted beside a measured one. "
                "It carries its own basis word for word, and it is not a check on the measurement at "
                "the same locus"
            ),
        },
        "a_count_is_not_a_measurement": A_COUNT_IS_NOT_A_MEASUREMENT,
        "no_outcome_is_read_against_another": OUTCOMES_ARE_NOT_AGREEMENT,
        "the_crispri_arm_holds_only_positives_by_construction": CRISPRI_ARM_IS_POSITIVES_ONLY,
    }


#: Why the CRISPRi outcomes on this map are the shape they are. A reader who takes the breakdown above
#: for the screen's own direction mix would be reading a selection rule as a result.
CRISPRI_ARM_IS_POSITIVES_ONLY = (
    "a pair reaches this map only with the status `all_three_present`, and that status requires the "
    "screen to have measured regulation at the pair. The pairs where it measured none are counted under "
    "`no_measured_rule_the_screen_measured_no_regulation` and are not built, so the CRISPRi breakdown "
    "above is over measured positives only. No rate, direction mix or share of significant effects may "
    "be read off it, for the screen or for anything else: the selection chose them"
)


#: Where each kind's outcome is recorded, named so a reader can go to the field rather than trust this
#: block, and so a kind cannot quietly start being summarised from somewhere else.
OUTCOME_FIELD = {
    "p2|": (
        "measurement.outcome: the benchmark's own significance call and direction for this pair, as "
        "`attribution.measured` recorded it. No threshold is applied here"
    ),
    "s2|": (
        "measurement.value: reader v1's own call at the element's interval. `not_open_in_reader` means "
        "not detected open at the reader's registered call, which is not closed"
    ),
    "m3|": (
        "aggregated_label.label: the measured layer's declared per-cell aggregation over this "
        "element's matched tiles, carried on each of that element's tile readings. It is NOT a reading "
        "of the tile, so an element's label is counted once for each of its tiles in this cell, and "
        "the tile's own observation is its log2(RNA/DNA) value"
    ),
    "q3|": (
        "measurement.direction: the sign of GTEx's own slope for this variant-gene pair in this "
        "tissue, over the retained significant set. No threshold is applied here"
    ),
    "r2|": (
        "measurement.direction: which way the compiler's rule says the element acts on the gene it "
        "named. Predicted, not measured"
    ),
}


def _kinds_block(payload_loci: list[dict[str, Any]], kinds_per_locus: Counter[str]) -> dict[str, Any]:
    return {
        "every_locus_here_carries": [
            "crispri_perturbation: the element silenced in place and the gene's expression read out",
            "reader_dnase: the element's own accessibility in the biosample, at the reader's call",
        ],
        "beyond_those_two": dict(sorted(kinds_per_locus.items())),
        "what_each_kind_observes": {
            "crispri_perturbation": ms.OUTCOME_KIND["crispri"],
            "reader_dnase": ce.OPENNESS_CALL,
            "lentimpra_sequence": ms.OUTCOME_KIND["lentimpra"],
            "eqtl_native_locus": eqtl.EVIDENCE,
        },
        "which_kinds_read_the_native_locus": {
            "crispri_perturbation": True,
            "reader_dnase": True,
            "eqtl_native_locus": True,
            "lentimpra_sequence": False,
        },
        "reporter_is_not_the_locus": REPORTER_NOT_THE_LOCUS,
        "eqtl_is_an_association": EQTL_UNCERTAINTY[0],
        "eqtl_frame": EQTL_FRAME,
        "more_kinds_is_not_more_established": MORE_KINDS_IS_NOT_MORE_ESTABLISHED,
        "absent_kinds": {
            "vista": (
                f"{ms.OUTCOME_KIND['vista']}. No element of this map carries one, so VISTA contributes "
                "nothing here; and a transgenic reporter in a mouse embryo would in any case be a "
                "reading of the sequence and not of the native locus"
            ),
            "satmut": (
                f"{ms.OUTCOME_KIND['satmut']}. No element of this map carries one. "
                f"{ms.SATMUT_CANNOT_DISAGREE}"
            ),
        },
        "loci": len(payload_loci),
    }


CANNOT_ESTABLISH = (
    "nothing here is a verdict, and no quantity on this map was fitted, scored or compared against a "
    "threshold of this increment's own",
    "several kinds of evidence at one locus are not agreement, confirmation or validation of one "
    "another. " + MORE_KINDS_IS_NOT_MORE_ESTABLISHED,
    "a lentiMPRA reading is an observation of a copy of the sequence in a reporter and never of the "
    "native locus after a perturbation of it. " + SILENT_IS_NOT_A_CONTRADICTION,
    "a GTEx cis-eQTL is an association over donors, in a tissue that is not the CRISPRi cell line. "
    "Linkage carries the signal some distance, so the variant need not be causal and the eGene need "
    "not be regulated by the interval it sits in",
    "an interval showing no eQTL here may simply be outside the frame the retained set was distilled "
    "against. " + EQTL_FRAME,
    "the locus convention is imported from `attribution.cell2` and is an operational grouping that "
    "pools cell types, not established biological independence",
    "the compiled rule beside each perturbation is predicted, never observed, and it carries its own "
    "basis word for word; a predicted quantity on this map is not a measurement of anything",
    "the loci here are the ones that happen to carry a further kind of evidence, which is a view of "
    "where the assays were pointed and not of the genome",
)


def _why(chroms: list[str], population: dict[str, Any]) -> dict[str, Any]:
    return {
        "question": (
            "where more than one kind of evidence exists at a locus of the measured perturbation "
            "layer, what does each kind actually say, and what does having several of them not "
            "establish?"
        ),
        "increment_2_stated_limit": INCREMENT_2_LIMIT,
        "what_this_increment_adds": (
            "two further kinds of observation beside increment 2's three assertions: a lentiMPRA tile "
            "reading of a copy of the element's sequence, and a retained GTEx cis-eQTL inside the "
            "element. Increment 2's own three assertions are built by its own makers, imported "
            "unchanged"
        ),
        "selection_rule": (
            "every addable locus of the measured perturbation layer that carries at least one of "
            "those two further kinds. The set was counted before this was built, at "
            f"{population['path']} (sha256 {population['sha256']}), and the build refuses to run over "
            "a set of a different size"
        ),
        "chromosomes": chroms,
        "no_cut_off_is_introduced": (
            "the overlap rule, the locus span, the openness call, the reporter activity threshold and "
            "its label rule, and GTEx's own significance set are all imported from where they were "
            "registered. This increment sets no threshold of its own"
        ),
    }


def _conventions(assembly: dict[str, Any], margin: int) -> dict[str, Any]:
    return {
        "assembly": assembly["assembly"],
        "coordinates": assembly["coordinates"],
        "assembly_basis": ", ".join(assembly["from"]),
        "independent_locus_rule": cell2.INDEPENDENT_LOCUS_RULE,
        "not_biological_independence": (
            "the grouping is operational and pools cell types: it joins two keys on a shared measured "
            "gene or on 1 Mb proximity on one chromosome and never reads the cell. It is not "
            "established biological independence and nothing may cite it as one"
        ),
        "reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
        "attachment": rm2.ATTACHMENT_NOTE,
        "openness_call": ce.OPENNESS_CALL,
        "not_validation": ce.NOT_VALIDATION,
        "not_closed": ce.NOT_CLOSED,
        "reporter_unit": list(ms.REPORTER_UNIT),
        "reporter_label_rule": ms.REPORTER_LABEL_RULE,
        "reporter_threshold_log2": ms.MPRA_ACTIVE,
        "reporter_uncertainty": ms.REPORTER_UNCERTAINTY,
        "eqtl_window_bp": margin,
        "eqtl_window_basis": (
            "the margin the retained GTEx set was distilled under, kept at `s < pos <= e` over the "
            "element widened by it, which is `eqtl.Intervals.at`'s own window"
        ),
        "eqtl_frame": EQTL_FRAME,
        "eqtl_gene_naming": (
            "GTEx names a gene by its Ensembl id and the CRISPRi benchmark by symbol, so an eQTL "
            "assertion's gene entity carries the Ensembl id as given. Two entities naming one gene "
            "under two identifiers are not joined here, and no gene identity across the two kinds is "
            "asserted by this map"
        ),
        "predicted_relation": rm2.PREDICTED_RELATION_NOTE,
        "reader_relation": rm2.READER_RELATION_NOTE,
        "reporter_relation": REPORTER_RELATION_NOTE,
        "eqtl_relation": EQTL_RELATION_NOTE,
        "chain": CHAIN_NOTE,
        "chain_stops": list(CHAIN_STOPS),
    }


def _reconcile(
    loci: list[dict[str, Any]],
    population: dict[str, Any],
    excluded: dict[str, list[dict[str, Any]]],
    pairs_excluded: dict[str, int],
    counted: dict[str, Any],
) -> dict[str, Any]:
    by_cause = {c: len(excluded[c]) for c in CAUSES}
    total = len(loci) + sum(by_cause.values())
    want = population["loci_carrying_at_least_one_further_kind"]
    return {
        "loci_counted_before_the_build": want,
        "counted_at": population["path"],
        "counted_sha256": population["sha256"],
        "loci_built": len(loci),
        "loci_excluded_by_cause": by_cause,
        "loci_accounted_for": total,
        "reconciles": total == want,
        "rule": (
            "loci_built + sum(loci_excluded_by_cause) == loci_counted_before_the_build. `build` "
            "refuses to return a payload where this does not hold, so a locus can never leave the "
            "map without a named cause"
        ),
        "causes": list(CAUSES),
        "pairs_excluded_by_cause": dict(pairs_excluded),
        "pairs_of_the_measured_layer_by_status": counted["by_status"],
        "addable_loci_of_the_whole_layer": population["addable_loci"],
        "share_of_the_addable_loci": (
            round(len(loci) / population["addable_loci"], 4) if population["addable_loci"] else None
        ),
        "what_that_share_is_not": (
            "the share of the addable loci that carry a further kind of evidence, which is a fact "
            "about where the assays were pointed. It is not a share of the genome, of regulatory "
            "elements, or of anything established"
        ),
    }


# --- paging ------------------------------------------------------------------------------------------
#: The blocks that describe the map rather than its loci. Carried on page 1 only, because repeating
#: them on every page is what made increment 2's eleven pages eleven copies of one preamble.
DESCRIBES_THE_MAP = (
    "why",
    "conventions",
    "vocabulary",
    "kinds_of_evidence",
    "cannot_establish",
    "counted",
    "population",
    "not_open_in_reader",
    "sources",
    "result_manifest",
    "code_cleanliness",
)


def page(payload: dict[str, Any], number: int = 1, per: int = PAGE) -> dict[str, Any]:
    """One page of the loci, with the whole map's counts beside it and nothing repeated.

    Every page carries `loci_total`, `assertions_total`, `pages` and the reconciliation, so a page is
    never a truncation. The blocks that describe the map rather than its loci are carried on page 1
    only and named in `described_on_page_1` on the others, so reading the map through does not mean
    reading one preamble eleven times.
    """
    per = max(1, int(per))
    loci = payload["loci"]
    pages = max(1, -(-len(loci) // per))
    number = min(max(1, int(number)), pages)
    shown = loci[(number - 1) * per : number * per]
    want = {aid for x in shown for aid in x["assertions"]}
    wanted_chains = {cid for x in shown for cid in x["chains"]}
    assertions = [a for a in payload["assertions"] if a["id"] in want]
    entities = {e for a in assertions for e in a["entities"]}
    drop = {"loci", "assertions", "entities", "chains", "cycles"}
    if number != 1:
        drop |= set(DESCRIBES_THE_MAP)
    out = {k: v for k, v in payload.items() if k not in drop}
    out.update(
        {
            "loci": shown,
            "assertions": assertions,
            "entities": {k: v for k, v in payload["entities"].items() if k in entities},
            "chains": [c for c in payload["chains"] if c["id"] in wanted_chains],
            "cycles": [c for c in payload["cycles"] if set(c["assertions"]) <= want],
            "page": number,
            "pages": pages,
            "per_page": per,
            "loci_total": len(loci),
            "assertions_total": len(payload["assertions"]),
            "counts": payload["counts"],
            "reconciliation": payload["reconciliation"],
            "paging": (
                f"page {number} of {pages}: loci {(number - 1) * per + 1} to "
                f"{(number - 1) * per + len(shown)} of {len(loci)}. Nothing is truncated; the counts "
                "beside this page are of the whole map"
            ),
        }
    )
    if number != 1:
        out["described_on_page_1"] = sorted(set(DESCRIBES_THE_MAP) & set(payload))
        out["not_repeated"] = (
            "the blocks naming this map's question, conventions, vocabulary, kinds of evidence, "
            "sources and limits are on page 1 and are the same for every page, so they are not "
            "repeated here. The counts and the reconciliation beside this page are of the whole map"
        )
    return out


def summary(payload: dict[str, Any]) -> str:
    c, r = payload["counts"], payload["reconciliation"]
    k = payload["kinds_of_evidence"]
    return "\n".join(
        [
            f"{payload['title']} (response map, increment {payload['increment']})",
            payload["header"],
            f"built {r['loci_built']} loci of {r['loci_counted_before_the_build']} counted "
            f"({r['counted_at']}); reconciles {r['reconciles']}",
            "excluded by cause: "
            + (", ".join(f"{x} {y}" for x, y in r["loci_excluded_by_cause"].items() if y) or "none"),
            f"{c['assertions']} assertions: " + ", ".join(f"{x} {y}" for x, y in c["by_status"].items()),
            "observed by assay: " + ", ".join(f"{x} {y}" for x, y in c["observed_by_assay"].items()),
            "kinds beyond the perturbation and the reader: "
            + ", ".join(f"{x} {y}" for x, y in k["beyond_those_two"].items()),
            "reads the native locus: "
            + ", ".join(
                f"{x} {'yes' if y else 'no'}" for x, y in k["which_kinds_read_the_native_locus"].items()
            ),
            MORE_KINDS_IS_NOT_MORE_ESTABLISHED,
        ]
    )
