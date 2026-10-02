# SPDX-License-Identifier: AGPL-3.0-or-later
"""The cellular control and response map, increment 2: every locus that carries all three inputs.

Increment 1 (`genomeos.response_map`) arranged one example, the beta-like globin genes in K562, by
reference. This increment arranges the same kind of chain for every locus of the measured layer's
perturbation assay that holds all three of the inputs such a chain needs, as
`data/results/response_map_coverage.json` counted them before anything was built:

- *a measured perturbation* - a CRISPRi element-gene pair of the ENCODE benchmark;
- *a compiled rule* - the rule the compiler emits for the element the pair attaches to under the
  measured layer's own overlap rule (`measured.RECIPROCAL_OVERLAP`);
- *a reader state* - `attribution.context_evidence.state_for`, lane-context2's own call, unchanged,
  on the rule's cell and the element's locus; `not_assessable` is not a state.

What it adds is *arrangement*, not evidence. Every assertion is one record of one file, named with
its path, the record key inside it and the file's sha256, and refused when it cannot name all three.
`observed` is only what an experiment measured: the CRISPRi effect and the DNase-seq peak reading.
A compiled rule is `predicted`, however strong. No validation, no evaluator verdict, no new cut-off,
no score, no fit: `response_map.check` is the same gate increment 1 passes, and it is run here.

**What this does not cover.** The perturbation assay holds 473 independent loci under the imported
1 Mb convention, so these loci are 23.26% of one of the measured layer's four assays; the other three
(lentiMPRA, VISTA, saturation mutagenesis) read a sequence or the bases inside it, not the native
locus after a perturbation of it, and no locus of theirs is here. The locus grouping is the
operational convention imported from `attribution.cell2`, **not established biological
independence**. A reader state is never validation: `context_evidence.NOT_VALIDATION` says why, and
`not_open_in_reader` means not detected open at the reader's registered call, never closed.

    uv run --frozen python scripts/response_map_increment2.py
    uv run --frozen python scripts/response_map_increment2.py --chroms chr21 \
        --result response_map_increment2_chr21
"""

from __future__ import annotations

import bisect
import copy
import hashlib
from collections import Counter
from pathlib import Path
from typing import Any, NamedTuple

from genomeos import response_map as rm
from genomeos.attribution import cell2
from genomeos.attribution import compile as cp
from genomeos.attribution import context_evidence as ce
from genomeos.attribution import measured as ms
from genomeos.attribution import targets as tg
from genomeos.genome import reader
from genomeos.results import RESULTS_DIR

INCREMENT = 2
VIEW = rm.VIEW
RESULT = "response_map_increment2"
CHROMS = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y"])

HEADER = (
    "A read-only map of existing evidence and of where it stops, for every locus of the measured "
    "perturbation assay that carries all three inputs. It adds no biological validation or "
    "prediction, and no evaluator verdict."
)

#: The count this increment is built against, and the file that produced it. Lane-rmap2 counted the
#: denominator before anything was built, so the number is not re-derived here: this module
#: enumerates the loci itself and then reconciles its own count against this one, naming every
#: exclusion by cause. A different sha256 is reported as a moved denominator and never adopted in
#: silence.
COVERAGE = "data/results/response_map_coverage.json"
COVERAGE_SHA256 = "ad0978b9269b59545f7f35dae96213f3c65131645f3cd52979fb1e71bf6b9691"
ADDABLE_LOCI = 110
ADDABLE_PAIRS = 195
MEASURED_LAYER_LOCI = 473
#: What increment 1 already showed, under its own unit of one example. It is one of the ADDABLE_LOCI,
#: so the ceiling after this increment is ADDABLE_LOCI and not ADDABLE_LOCI + 1.
LOCI_INCREMENT_1_COVERED = 1

#: The statuses of the count, carried by name so this module and `scripts/response_map_coverage.py`
#: cannot answer differently about which input a pair is missing. `tests/test_response_map2.py` pins
#: each one against that script.
STATUS_ADDABLE = "all_three_present"
STATUS_NO_READER = "no_reader_state"
STATUS_NO_RULE = "no_measured_rule_the_screen_measured_no_regulation"
STATUS_NO_LINK = "no_compiled_link_attaches_to_no_compiled_element"
STATUSES = (STATUS_ADDABLE, STATUS_NO_READER, STATUS_NO_RULE, STATUS_NO_LINK)
STATUS_ORDER = {s: i for i, s in enumerate(STATUSES)}

#: Why a locus with all three inputs can still fail to be built, each cause named rather than a
#: locus disappearing quietly. The build reports a count per cause and refuses a payload whose built
#: loci plus its exclusions do not come back to ADDABLE_LOCI.
CAUSE_ELEMENT_SOURCE = "element_source_file_not_in_this_checkout"
CAUSE_PAIR_SOURCE = "perturbation_source_file_not_in_this_checkout"
CAUSE_READER_SOURCE = "reader_peak_file_not_in_this_checkout"
CAUSE_DUPLICATE_KEY = "two_records_claim_one_observation_key"
CAUSES = (CAUSE_ELEMENT_SOURCE, CAUSE_PAIR_SOURCE, CAUSE_READER_SOURCE, CAUSE_DUPLICATE_KEY)

ATTACHMENT_NOTE = (
    "the measured interval and the compiled element are two intervals, not one: the pair is a "
    "measurement of this element only because it meets the measured layer's own overlap rule "
    f"(RECIPROCAL_OVERLAP = {ms.RECIPROCAL_OVERLAP}, both ways), which is how "
    "`measured.rows` attaches a pair to an element. The overlap is carried on the assertion"
)
PREDICTED_RELATION_NOTE = (
    "the relation names the quantity the scorer predicts, an RNA-seq expression change on deleting "
    "the element; the status says it was predicted and not measured, and no part of this is a check "
    "on the measurement beside it"
)
READER_RELATION_NOTE = (
    "the state is a reading of one ENCODE DNase-seq peak set at the element's own interval, and "
    "nothing more: the assertion says which peak set, at which interval, and what it found"
)
NAMED_INSTEAD = (
    "the compiled rule names a different gene for this element than the screen measured; both are "
    "shown, neither is resolved, and the disagreement is not scored here"
)
CHAIN_NOTE = (
    rm.CHAIN_NOTE + " The step from the measured interval to the compiled element is the measured "
    "layer's overlap rule and not an identity, and the reader state is a consistency check between "
    "two readings of the same ENCODE chromatin, never support for the rule."
)
CHAIN_STOPS = (
    "readout_not_recorded",
    "mechanism_not_established",
    "reader_state_is_not_validation",
    "locus_grouping_is_operational",
)

#: One page of loci the web serves at a time. 110 loci do not fit one screen, so the page is cut into
#: pages server-side and every page carries the whole count beside it; nothing is truncated.
PAGE = 10


class Candidate(NamedTuple):
    """One measured pair, with the status the three-input rule gives it and what it attached to."""

    pair: Any
    status: str
    element: dict[str, Any] | None
    row: dict[str, Any] | None
    overlap: float | None
    state: str
    biosample: str | None
    key: cell2.LocusKey


def attaches_to(
    pair: Any, starts: list[int], elements: list[dict[str, Any]]
) -> tuple[dict[str, Any] | None, float | None]:
    """The first compiled element the measured layer's overlap rule makes this pair a measurement of.

    The same scan as `scripts/response_map_coverage.py`: `measured.REACH` wide on each side, which is
    the longest measured interval any assay holds, so no element that could meet the rule lies
    outside it. The overlap fraction is returned beside the element because the map shows it.
    """
    lo = bisect.bisect_left(starts, pair.start - ms.REACH)
    hi = bisect.bisect_right(starts, pair.end + ms.REACH)
    for e in elements[lo:hi]:
        if ms.measures(e["start"], e["end"], pair.start, pair.end):
            return e, round(ms.reciprocal_overlap(e["start"], e["end"], pair.start, pair.end), 3)
    return None, None


def candidates(
    chroms: list[str],
    results_dir: Path = RESULTS_DIR,
    readers: ce.Readers | None = None,
    table: ce.Mapping | None = None,
) -> tuple[list[Candidate], dict[str, Any]]:
    """Every measured pair of these chromosomes with its three-input status, in the counting order.

    The order is the one `scripts/response_map_coverage.py` grouped in - the chromosomes as given,
    the pairs as `measured.load_crispri` sorts them - because `cell2.group` names a group after the
    index of its first member and the two must agree on which group is which.
    """
    readers = readers if readers is not None else ce.Readers()
    table = table if table is not None else ce.mapping()
    out: list[Candidate] = []
    per_chrom: dict[str, Counter[str]] = {}
    invalid = 0
    for chrom in chroms:
        elements = cp._attributed(chrom, results_dir)
        starts = [e["start"] for e in elements]
        _layer, measured_rows = cp._measured_rows(chrom, elements, results_dir)
        rule_of: dict[tuple[str, str, str], dict[str, Any]] = {}
        for r in measured_rows:
            for gene, _action, _strength, cell, _split in ms.rule_links(r):
                rule_of[(r["id"], gene, cell)] = r
        pairs, bad = ms.load_crispri(chrom)
        invalid += bad
        here: Counter[str] = Counter()
        for p in pairs:
            e, overlap = attaches_to(p, starts, elements)
            state, biosample, row = "", None, None
            if e is None:
                status, locus = STATUS_NO_LINK, (p.chrom, p.start, p.end)
            else:
                locus = (chrom, e["start"], e["end"])
                row = rule_of.get((e["id"], p.gene, p.cell))
                if row is None:
                    status = STATUS_NO_RULE
                else:
                    raw = ce.state_for(cp.context(p.cell), chrom, row["start"], row["end"], readers, table)
                    state, rest = ce.parse_value(raw)
                    biosample = rest[-1] if rest else None
                    status = STATUS_NO_READER if state == ce.STATE_NOT_ASSESSABLE else STATUS_ADDABLE
            out.append(
                Candidate(
                    pair=p,
                    status=status,
                    element=e,
                    row=row,
                    overlap=overlap,
                    state=state,
                    biosample=biosample,
                    key=cell2.LocusKey(p.cell, locus[0], locus[1], locus[2], p.gene),
                )
            )
            here[status] += 1
        per_chrom[chrom] = here
    counted = {
        "pairs": len(out),
        "invalid_rows_the_parser_rejected": invalid,
        "by_status": {s: sum(c[s] for c in per_chrom.values()) for s in STATUSES},
        "per_chromosome": {c: {s: v[s] for s in STATUSES} for c, v in per_chrom.items()},
    }
    return out, counted


# --- sources ---------------------------------------------------------------------------------------
class _Files:
    """The files the assertions name, hashed once each, with whether git holds them at HEAD.

    `response_map._Sources` reads a file whole to quote records out of it; here the records come from
    the measured layer and the compiler, which have already read them, so only the digest and the
    availability are needed. A file absent from this checkout cannot carry a record key and a
    sha256, so an assertion that needs it is excluded by cause instead of shown.

    Whether a file is committed is read from git rather than from its directory: reader v1's peak sets
    live under `data/results/` and are git-ignored all the same, so the path alone would call them
    committed and the assertion would say something untrue about its own source.
    """

    def __init__(self, root: Path) -> None:
        self.root = Path(root)
        self.table: dict[str, dict[str, Any]] = {}
        self.at_head = _tracked_at_head(self.root)

    def committed(self, path: str) -> bool:
        if self.at_head is None:  # no git here: fall back to increment 1's path rule, and say so
            return rm._Sources.committed(path)
        return path in self.at_head

    def hash(self, path: str) -> dict[str, Any]:
        if path not in self.table:
            p = self.root / path
            digest = size = None
            if p.is_file():
                h = hashlib.sha256()
                with p.open("rb") as fh:
                    for chunk in iter(lambda: fh.read(1 << 20), b""):
                        h.update(chunk)
                digest, size = h.hexdigest(), p.stat().st_size
            entry = {
                "path": path,
                "kind": rm.COMMITTED if self.committed(path) else rm.LOCAL,
                "committed_at_head": self.committed(path),
                "committed_read_from": "git ls-tree HEAD" if self.at_head is not None else "the path",
                "available": digest is not None,
                "sha256": digest,
                "bytes": size,
            }
            if digest is None:
                entry["note"] = rm.UNAVAILABLE
            self.table[path] = entry
        return self.table[path]

    def ref(self, path: str, key: str, record_path: list[Any], **extra: Any) -> dict[str, Any] | None:
        """The reference an assertion carries, or None when the file is not in this checkout."""
        s = self.hash(path)
        if not s["available"]:
            return None
        return {
            "path": path,
            "record_key": key,
            "record_path": record_path,
            "sha256": s["sha256"],
            "kind": s["kind"],
            **extra,
        }


def element_source(chrom: str, origin: str, results_dir: Path = RESULTS_DIR) -> tuple[str, str] | None:
    """The file one attributed element was read from, and how, for the run that produced it.

    `attribution.targets` reads three runs per chromosome and marks each element with the one it came
    from. A run writes its elements inline into its own committed result, or writes a summary pointing
    at a local table; both are returned as a repository-relative path so the assertion can name it.
    """
    name = next((n for n, o in tg.ORIGIN.items() if o == origin), None)
    if name is None:
        return None
    from genomeos.results import load_result

    r = load_result(f"{name}_{chrom}", results_dir) or {}
    if "elements" in r:
        return f"{results_dir.as_posix()}/{name}_{chrom}.json", "elements"
    where = r.get("elements_where")
    return (str(where), "the table the result points at") if where else None


def pair_key(pair: Any) -> str:
    """The benchmark row's key inside its file, spelled as increment 1 spells it, with the cell and
    the split, because one (interval, gene, dataset) can be measured in more than one of each."""
    return f"{pair.chrom}:{pair.start}-{pair.end}|{pair.gene}|{pair.dataset}|{pair.cell}|{pair.split}"


# --- the assertions of one locus --------------------------------------------------------------------
_PAIR_FIELDS = (
    "chrom",
    "start",
    "end",
    "gene",
    "cell",
    "dataset",
    "reference",
    "regulated",
    "significant",
    "effect_size",
    "p_adjusted",
    "split",
    "source_file",
    "power_at_effect_size_20",
)
_ELEMENT_FIELDS = ("id", "start", "end", "domain", "origin", "predicted_coding", "verdict_coding")


def _state_entity(m: rm._Map, state: str, biosample: str) -> str:
    """The measured cell state a reader state is an association with, as its own entity."""
    return m.entity(
        f"state:{state}:{biosample}",
        "cell_state",
        state=state,
        biosample=biosample,
        reading=ce.READING[state],
        assay=reader.EVIDENCE,
        is_validation=False,
        not_validation=ce.NOT_VALIDATION,
        **({"not_closed": ce.NOT_CLOSED} if state == ce.STATE_NOT_OPEN else {}),
    )


def _perturbation(
    c: Candidate, m: rm._Map, files: _Files, element_iv: str, basis: str
) -> dict[str, Any] | None:
    """The CRISPRi pair as one assertion: observed, from the benchmark row, or None without a source."""
    p = c.pair
    path = (ms.CRISPRI_KNOWLEDGE / p.source_file).as_posix()
    ref = files.ref(path, pair_key(p), ["row", pair_key(p)], partition=p.split)
    if ref is None:
        return None
    iv = m.interval(p.chrom, p.start, p.end, basis, tested_by=f"{p.assay} ({p.dataset})")
    direction = "decrease" if p.effect_size < 0 else "increase"
    return {
        "id": f"p2|{p.chrom}:{p.start}-{p.end}|{p.gene}|{p.dataset}|{p.cell}|{p.split}",
        **rm._expression_relation(),
        "subject": iv,
        "object": m.gene(p.gene),
        "entities": [element_iv],
        "status": "observed",
        "basis": "measured perturbation effect (CRISPRi screen, ENCODE benchmark)",
        "context": {
            "cell": p.cell,
            "condition": rm.NOT_RECORDED,
            "perturbation": f"CRISPRi of the tested interval ({p.dataset})",
        },
        "measurement": {
            "quantity": "EffectSize: signed effect on target gene expression (the benchmark's column)",
            "value": round(p.effect_size, 4),
            "unit": ms.EFFECT_UNIT,
            "direction": direction,
            "outcome": p.outcome,
            "significant": bool(p.significant),
            "power_at_effect_size_20": p.power_at_effect_size_20,
        },
        "split": p.split,
        "split_basis": ms.HELDOUT_MARK if p.split == ms.HELDOUT else "the benchmark's training file",
        "mechanism": "not established: the change after perturbation may be indirect",
        "attached_to_element": {
            "element": c.element["id"] if c.element else None,
            "entity": element_iv,
            "reciprocal_overlap": c.overlap,
            "rule": ATTACHMENT_NOTE,
        },
        "source": ref,
        "quoted": {k: getattr(p, k) for k in _PAIR_FIELDS},
        "assembly_basis": basis,
        "observation_key": pair_key(p),
    }


def _compiled_rule(
    c: Candidate, m: rm._Map, files: _Files, chrom: str, results_dir: Path, basis: str
) -> dict[str, Any] | None:
    """The compiler's own predicted rule for the element: predicted, never observed."""
    e = c.element
    assert e is not None
    found = element_source(chrom, e.get("origin", ""), results_dir)
    if found is None:
        return None
    path, how = found
    ref = files.ref(path, f"element {e['id']}", ["elements", e["id"]], read_as=how)
    if ref is None:
        return None
    pc = e.get("predicted_coding") or {}
    basis = str(pc.get("basis") or rm.NOT_RECORDED)
    iv = m.interval(
        chrom,
        e["start"],
        e["end"],
        basis,
        registry_ids=[e["id"]],
        domain=e.get("domain", ""),
        compiled_element=e["id"],
    )
    measured_gene = c.pair.gene
    conflicts = []
    if pc.get("gene") and pc["gene"] != measured_gene:
        conflicts.append(
            {
                "kind": "model_named_instead",
                "detail": f"for this element the compiler named {pc['gene']}, not {measured_gene}",
                "model_named": [pc["gene"]],
                "status": "predicted",
                "resolution": "unresolved",
                "note": NAMED_INSTEAD,
                "source": {
                    "path": path,
                    "record_key": f"element {e['id']}",
                    "field": "predicted_coding.gene",
                },
            }
        )
    return {
        "id": f"r2|{chrom}:{e['start']}-{e['end']}|{e['id']}|{pc.get('gene') or 'none'}",
        "relation": "changes_measured_rna",
        "relation_note": PREDICTED_RELATION_NOTE,
        "subject": iv,
        "object": m.gene(pc["gene"]) if pc.get("gene") else None,
        "status": "predicted",
        "basis": basis,
        "context": {
            "cell": rm.NOT_RECORDED,
            "condition": rm.NOT_RECORDED,
            "perturbation": "in silico deletion of the element (the scorer's own perturbation)",
        },
        "measurement": {
            "quantity": "predicted log2 fold change in the gene's expression on deleting the element",
            "value": pc.get("log2_fold_change"),
            "unit": "log2 fold change",
            "direction": pc.get("action"),
            "tissue_or_track": pc.get("tissue"),
            "is_a_measurement": False,
        },
        "mechanism": "not established: a predicted expression change names no mechanism",
        "conflicts": conflicts,
        "source": ref,
        "quoted": {k: copy.deepcopy(e[k]) for k in _ELEMENT_FIELDS if k in e},
        "assembly_basis": basis,
    }


def _reader_state(
    c: Candidate,
    m: rm._Map,
    files: _Files,
    chrom: str,
    readers: ce.Readers,
    element_iv: str,
    basis: str,
) -> dict[str, Any] | None:
    """The element's state in the pair's cell, as reader v1 calls it: observed, and not validation."""
    if c.biosample is None or c.state not in (ce.STATE_OPEN, ce.STATE_NOT_OPEN):
        return None
    path = (RESULTS_DIR / reader.peaks_path(c.biosample, chrom).name).as_posix()
    row = c.row
    assert row is not None
    key = f"{c.biosample} {chrom}:{row['start']}-{row['end']}"
    span = readers.span(c.biosample, chrom)
    ref = files.ref(path, key, ["peaks", chrom, [row["start"], row["end"]]])
    if ref is None:
        return None
    return {
        "id": f"s2|{chrom}:{row['start']}-{row['end']}|{c.biosample}|{c.state}",
        "relation": "associated_with_state",
        "relation_note": READER_RELATION_NOTE,
        "subject": element_iv,
        "object": _state_entity(m, c.state, c.biosample),
        "status": "observed",
        "basis": ce.OPENNESS_CALL,
        "context": {
            "cell": c.pair.cell,
            "condition": rm.NOT_RECORDED,
            "perturbation": "none: an unperturbed accessibility reading",
        },
        "measurement": {
            "quantity": "overlap of the element's interval with a DNase narrowPeak of the biosample",
            "value": c.state,
            "unit": rm.NOT_RECORDED,
            "reading": ce.READING[c.state],
            "measured_span_of_the_biosample": list(span) if span else None,
            "span_call": ce.SPAN_CALL,
        },
        "is_validation": False,
        "uncertainty": [ce.NOT_VALIDATION] + ([ce.NOT_CLOSED] if c.state == ce.STATE_NOT_OPEN else []),
        "cell_mapping": {
            "rule_cell_label": cp.context(c.pair.cell),
            "reader_biosample": c.biosample,
            "term": ce.READER_TERMS.get(c.biosample),
        },
        "source": ref,
        "assembly_basis": basis,
    }


# --- the build --------------------------------------------------------------------------------------
def _locus_id(keys: list[cell2.LocusKey]) -> str:
    chroms = sorted({k.chrom for k in keys})
    span = f"{min(k.start for k in keys)}-{max(k.end for k in keys)}"
    return "locus|" + "+".join(chroms) + ":" + span


def map1_keys(root: Path) -> list[cell2.LocusKey]:
    """Increment 1's own measured perturbation assertions as locus keys, so the two increments agree
    on which locus the map already covered. Read from `response_map.view`, not restated."""
    payload = rm.view(Path(root))
    entities = payload["entities"]
    out = []
    for a in payload["assertions"]:
        if not str(a.get("basis", "")).startswith("measured perturbation effect"):
            continue
        iv = entities[a["subject"]]
        out.append(
            cell2.LocusKey(
                a["context"]["cell"],
                iv["chrom"],
                int(iv["start"]),
                int(iv["end"]),
                entities[a["object"]]["symbol"],
            )
        )
    return out


def build(
    root: Path = Path(),
    chroms: tuple[str, ...] | list[str] = CHROMS,
    results_dir: Path = RESULTS_DIR,
    readers: ce.Readers | None = None,
    table: ce.Mapping | None = None,
) -> dict[str, Any]:
    """The map of every addable locus, checked by `response_map.check` and reconciled against the
    count it was built from. Refuses rather than serve a payload whose loci do not reconcile."""
    chroms = list(chroms)
    readers = readers if readers is not None else ce.Readers()
    table = table if table is not None else ce.mapping()
    cands, counted = candidates(chroms, results_dir, readers, table)

    covered = map1_keys(root)
    groups = cell2.group([c.key for c in cands] + covered)
    mine, theirs = groups[: len(cands)], set(groups[len(cands) :])
    best: dict[int, str] = {}
    for g, c in zip(mine, cands, strict=True):
        if g not in best or STATUS_ORDER[c.status] < STATUS_ORDER[best[g]]:
            best[g] = c.status
    addable = {g for g, s in best.items() if s == STATUS_ADDABLE}
    members: dict[int, list[int]] = {}
    for i, g in enumerate(mine):
        members.setdefault(g, []).append(i)

    files = _Files(root)
    assembly = rm._assembly(rm._Sources(Path(root)))  # GRCh38, half-open, or a refusal
    basis = ", ".join(assembly["from"])
    m = rm._Map(assembly)
    loci: list[dict[str, Any]] = []
    chains: list[dict[str, Any]] = []
    excluded: dict[str, list[dict[str, Any]]] = {c: [] for c in CAUSES}
    pairs_excluded: dict[str, int] = dict.fromkeys(CAUSES, 0)
    seen_keys: set[str] = set()
    state_counts: Counter[str] = Counter()
    not_open: list[dict[str, Any]] = []

    added: set[str] = set()
    for g in sorted(addable):
        idx = members[g]
        keys = [cands[i].key for i in idx]
        lid = _locus_id([cands[i].key for i in idx if cands[i].status == STATUS_ADDABLE])
        built: list[dict[str, Any]] = []
        causes: list[str] = []
        for i in idx:
            c = cands[i]
            if c.status != STATUS_ADDABLE:
                continue
            assert c.element is not None and c.row is not None
            element_iv = rm._iv(c.pair.chrom, c.element["start"], c.element["end"])
            pert = _perturbation(c, m, files, element_iv, basis)
            rule = _compiled_rule(c, m, files, c.pair.chrom, results_dir, basis)
            state = _reader_state(c, m, files, c.pair.chrom, readers, element_iv, basis)
            if pert is None:
                cause = CAUSE_PAIR_SOURCE
            elif rule is None:
                cause = CAUSE_ELEMENT_SOURCE
            elif state is None:
                cause = CAUSE_READER_SOURCE
            elif pert["observation_key"] in seen_keys:
                cause = CAUSE_DUPLICATE_KEY
            else:
                seen_keys.add(pert["observation_key"])
                built.append({"candidate": c, "pert": pert, "rule": rule, "state": state})
                continue
            causes.append(cause)
            pairs_excluded[cause] += 1
        if not built:
            excluded[causes[0]].append(
                {
                    "locus": lid,
                    "chrom": keys[0].chrom,
                    "genes": sorted({k.gene for k in keys}),
                    "causes_of_its_pairs": sorted(set(causes)),
                }
            )
            continue
        ids: list[str] = []
        for b in built:
            for a in (b["pert"], b["rule"], b["state"]):
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
            ids += steps
            chains.append(
                {
                    "id": f"chain|{b['pert']['id']}",
                    "title": (
                        f"{c.pair.chrom}:{c.pair.start}-{c.pair.end} -> {c.element['id']} -> "
                        f"{c.pair.gene} in {c.pair.cell}"
                    ),
                    "steps": steps,
                    "beside": [],
                    "derived_claims": [],
                    "stops_at": list(CHAIN_STOPS),
                    "note": CHAIN_NOTE,
                }
            )
        shown_keys = [cands[i].key for i in idx if cands[i].status == STATUS_ADDABLE]
        loci.append(
            {
                "id": lid,
                "chrom": keys[0].chrom,
                "chromosomes": sorted({k.chrom for k in keys}),
                # the span of the pairs this locus shows, which is what `id` is made of. The grouping
                # also holds pairs that are not shown, and their span is wider: both are reported, so
                # the interval beside a locus is never read as covering more than the map shows.
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
                "chains": [f"chain|{b['pert']['id']}" for b in built],
                "reader_states": sorted({b["candidate"].state for b in built}),
                "pairs_shown": len(built),
                "pairs_in_this_locus_not_shown": {
                    s: sum(1 for i in idx if cands[i].status == s) for s in STATUSES[1:]
                },
                "already_in_increment_1": g in theirs,
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
        title="The loci of the measured perturbation layer that carry all three inputs",
        header=HEADER,
        why=_why(chroms),
        conventions=_conventions(assembly),
        vocabulary={
            "entity_kinds": list(rm.ENTITY_KINDS),
            "relations": list(rm.RELATIONS),
            "statuses": list(rm.STATUSES),
            "categories": rm.CATEGORIES,
        },
        loci=loci,
        page_size=PAGE,
        counted=counted,
        coverage=_coverage(counted, loci, len(best), list(chroms) == list(CHROMS)),
        not_open_in_reader=_not_open_block(not_open, state_counts),
        requests={"network": 0, "model": 0, "downloads": 0},
    )
    payload["reconciliation"] = _reconcile(
        files, chroms, loci, best, addable, excluded, pairs_excluded, counted
    )
    payload["counts"]["loci"] = len(loci)
    payload["counts"]["observed_by_assay"] = _observed_by_assay(payload["assertions"])
    payload["sources"] = [files.table[p] for p in sorted(files.table)]
    return payload


def _observed_by_assay(assertions: list[dict[str, Any]]) -> dict[str, int]:
    """`observed` split by what measured it, so the figure is never read as one experiment's."""
    out: Counter[str] = Counter()
    for a in assertions:
        if a["status"] != "observed":
            continue
        out["CRISPRi perturbation (ENCODE benchmark)" if a["id"].startswith("p2|") else reader.EVIDENCE] += 1
    return dict(sorted(out.items()))


def _not_open_block(rows: list[dict[str, Any]], states: Counter[str]) -> dict[str, Any]:
    return {
        "shown_where_it_applies": rows,
        "shown": len(rows),
        "pairs_by_reader_state": dict(sorted(states.items())),
        "meaning": ce.READING[ce.STATE_NOT_OPEN],
        "never": ce.NOT_CLOSED,
        "not_validation": ce.NOT_VALIDATION,
        "population_it_comes_from": (
            "lane-context2's genome-wide census of all 440,589 compiled rules: 55,084 of the 81,635 "
            "assessable rules assert their element in a cell where the reader does not detect it "
            "open. The rules shown here are the ones a measured perturbation attaches to, which is a "
            "different and much smaller population, so no share of that census may be read off this "
            "map"
        ),
    }


def _coverage(
    counted: dict[str, Any], loci: list[dict[str, Any]], loci_in_scope: int, genome_wide: bool
) -> dict[str, Any]:
    """What the built loci are a share of. A run over fewer chromosomes is a share of its own scope,
    never of the genome-wide 473, so a partial build cannot be quoted as a genome-wide share."""
    total = MEASURED_LAYER_LOCI if genome_wide else loci_in_scope
    return {
        "loci_shown": len(loci),
        "independent_loci_of_the_perturbation_assay": total,
        "independent_loci_of_the_perturbation_assay_genome_wide": MEASURED_LAYER_LOCI,
        "share_of_the_perturbation_assay": round(len(loci) / total, 4) if total else None,
        "share_is_of": "all 24 chromosomes" if genome_wide else "the chromosomes of this run only",
        "assays_in_the_measured_layer": len(ms.ASSAYS),
        "what_it_covers": (
            f"{len(loci)} of the {total} independent loci of the CRISPRi arm "
            + ("genome-wide" if genome_wide else "on the chromosomes of this run")
            + ", under the operational 1 Mb grouping imported from `attribution.cell2`"
        ),
        "what_it_does_not_cover": (
            f"the other {total - len(loci)} loci of the same arm, named by cause in the count this "
            f"was built from; and the measured layer's other three assays "
            f"({', '.join(a for a in ms.ASSAYS if a != 'crispri')}), which read a sequence or the "
            "bases inside it and not the native locus after a perturbation of it, so none of their "
            "loci is here"
        ),
        "not_biological_independence": cell2.INDEPENDENT_LOCUS_RULE,
        "measured_pairs_seen": counted["pairs"],
    }


def _why(chroms: list[str]) -> dict[str, Any]:
    return {
        "chosen_for": (
            "completeness of the chain over a counted population rather than one example: every locus "
            "of the measured perturbation assay that holds a measured perturbation, a compiled rule "
            "for the element it attaches to, and a reader state, each assertion naming its file, its "
            "record key and that file's sha256"
        ),
        "plan": "docs/ROADMAP.md, area I, the response map's second increment",
        "counted_before_it_was_built": {
            "result": "response_map_coverage",
            "path": COVERAGE,
            "key": "addable.independent_loci",
            "what": (
                f"{ADDABLE_LOCI} independent loci carry all three inputs out of {MEASURED_LAYER_LOCI}, "
                f"{LOCI_INCREMENT_1_COVERED} of them already shown by increment 1; the denominator was "
                "counted before anything was built, so this increment adds no selection of its own"
            ),
        },
        "selection_rule": (
            "every locus whose best pair status is "
            f"`{STATUS_ADDABLE}` over the CRISPRi pairs of {len(chroms)} chromosomes, with the status "
            "taken per pair from the three inputs and per locus as the best any of its pairs reaches. "
            "No outcome, effect size, agreement or score takes part in the selection"
        ),
        "consequence": (
            "the population was fixed by a count that read no assertion, and every locus in it is "
            "shown or excluded by a named cause, so nothing here is a choice made after seeing a "
            "result"
        ),
        "it_is_not": [
            "not validation: a reader state is a consistency check between two readings of the same "
            "ENCODE chromatin, and AlphaGenome was trained on those tracks",
            "not a claim that the loci are independent: the grouping is the operational 1 Mb "
            "convention imported from attribution.cell2",
            "not an extension of coverage to the measured layer: it is one of four assays, and part "
            "of that one",
        ],
    }


def _conventions(assembly: dict[str, Any]) -> dict[str, Any]:
    return {
        **assembly,
        "not_recorded": rm.NOT_RECORDED,
        "unavailable": rm.UNAVAILABLE,
        "local": rm.LOCAL,
        "readout": rm.READOUT_NOTE,
        "attachment": ATTACHMENT_NOTE,
        "predicted_relation": PREDICTED_RELATION_NOTE,
        "reader_state": READER_RELATION_NOTE,
        "reader_state_is_not_validation": ce.NOT_VALIDATION,
        "not_open_is_not_closed": ce.NOT_CLOSED,
        "openness_call": ce.OPENNESS_CALL,
        "span_call": ce.SPAN_CALL,
        "locus_grouping": cell2.INDEPENDENT_LOCUS_RULE,
        "overlap_rule": ATTACHMENT_NOTE,
        "chains": CHAIN_NOTE,
        "cycles": rm.CYCLE_NOTE,
        "observed_means": "what an experiment measured; a model output is predicted, never observed",
    }


def _reconcile(
    files: _Files,
    chroms: list[str],
    loci: list[dict[str, Any]],
    best: dict[int, str],
    addable: set[int],
    excluded: dict[str, list[dict[str, Any]]],
    pairs_excluded: dict[str, int],
    counted: dict[str, Any],
) -> dict[str, Any]:
    """The built locus count against the one it was built from, with every exclusion named by cause."""
    genome_wide = list(chroms) == list(CHROMS)
    entry = files.hash(COVERAGE)
    recorded = {"path": COVERAGE, "sha256_recorded": COVERAGE_SHA256, "sha256_found": entry["sha256"]}
    recorded["denominator_file_is_the_one_recorded"] = entry["sha256"] == COVERAGE_SHA256
    recorded["denominator_file_is_committed"] = files.committed(COVERAGE)
    counts = {c: len(v) for c, v in excluded.items()}
    built = len(loci)
    total = built + sum(counts.values())
    expected = ADDABLE_LOCI if genome_wide else len(addable)
    out = {
        "unit": "one independent locus under the imported 1 Mb grouping",
        "denominator": recorded,
        "addable_loci_counted_before_the_build": expected,
        "addable_loci_counted_in_this_run": len(addable),
        "loci_built": built,
        "loci_excluded_by_cause": counts,
        "loci_excluded": excluded,
        "pairs_excluded_by_cause": pairs_excluded,
        "built_plus_excluded": total,
        "reconciles": total == expected and len(addable) == expected,
        "every_exclusion_is_named": sum(counts.values()) == sum(len(v) for v in excluded.values()),
        "already_covered_by_increment_1": sum(1 for x in loci if x["already_in_increment_1"]),
        "newly_covered_by_this_increment": sum(1 for x in loci if not x["already_in_increment_1"]),
        "ceiling_after_this_increment": expected,
        "the_ceiling_is_not_one_more": (
            f"the locus increment 1 covered is one of the {expected} addable loci, not a locus beside "
            f"them, so the map covers at most {expected} loci after this increment and not "
            f"{expected + 1}"
        ),
        "loci_by_best_status": {s: sum(1 for v in best.values() if v == s) for s in STATUSES},
        "pairs_by_status": counted["by_status"],
        "scope": ("all 24 chromosomes" if genome_wide else f"{len(chroms)} of 24: {', '.join(chroms)}"),
    }
    if not out["reconciles"]:
        raise rm.RefusedError(
            f"the built loci do not reconcile: {built} built plus {sum(counts.values())} excluded is "
            f"{total}, against {expected} counted ({len(addable)} addable in this run)"
        )
    return out


def _tracked_at_head(root: Path) -> set[str] | None:
    """Every path git holds at HEAD, in one call, or None where this is not a checkout.

    A file the map reads but nobody has committed is reported as such rather than assumed; the sha256
    pins the bytes either way.
    """
    import subprocess

    try:
        r = subprocess.run(
            ["git", "ls-tree", "-r", "HEAD", "--name-only"],
            cwd=str(root or Path()),
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        return None
    if r.returncode != 0:
        return None
    return set(r.stdout.splitlines())


def page(payload: dict[str, Any], number: int = 1, per: int = PAGE) -> dict[str, Any]:
    """One page of the loci with the whole count beside it, and the assertions those loci cite.

    110 loci and their assertions do not fit one screen, so the server cuts the payload into pages.
    Every page carries `loci_total` and `pages`, so a page is never a truncation: a reader can always
    see how many loci the map holds and which ones this page is.
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
    drop = ("loci", "assertions", "entities", "chains", "cycles")
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
            "paging": (
                f"page {number} of {pages}: loci {(number - 1) * per + 1} to "
                f"{(number - 1) * per + len(shown)} of {len(loci)}. Nothing is truncated; the counts "
                "beside this page are of the whole map"
            ),
        }
    )
    return out


def summary(payload: dict[str, Any]) -> str:
    c, r = payload["counts"], payload["reconciliation"]
    lines = [
        f"{payload['title']} (response map, increment {payload['increment']})",
        payload["header"],
        f"built {r['loci_built']} loci of {r['addable_loci_counted_before_the_build']} counted; "
        f"{r['newly_covered_by_this_increment']} not covered by increment 1; reconciles {r['reconciles']}",
        "excluded by cause: "
        + (", ".join(f"{k} {v}" for k, v in r["loci_excluded_by_cause"].items() if v) or "none"),
        f"{c['assertions']} assertions: " + ", ".join(f"{k} {v}" for k, v in c["by_status"].items()),
        "observed by assay: " + ", ".join(f"{k} {v}" for k, v in c["observed_by_assay"].items()),
        "reader states: "
        + ", ".join(f"{k} {v}" for k, v in payload["not_open_in_reader"]["pairs_by_reader_state"].items()),
        f"not_open_in_reader shown at {payload['not_open_in_reader']['shown']} of "
        f"{sum(payload['not_open_in_reader']['pairs_by_reader_state'].values())} shown pairs",
        payload["coverage"]["what_it_covers"],
        "it does not cover " + payload["coverage"]["what_it_does_not_cover"],
    ]
    return "\n".join(lines)
