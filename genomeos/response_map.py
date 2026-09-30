# SPDX-License-Identifier: AGPL-3.0-or-later
"""The cellular control and response map, increment 1: one example, read-only, by reference.

The map arranges evidence this project already holds, and the places where it stops, for one
example: the beta-globin locus in K562 (docs/ROADMAP.md, row "↳ MAP", commit 304a1c8). It is an
evidence-organisation view. It adds no biological validation or prediction, changes no label,
verdict, schema or result, and runs no model, request or download.

Three kinds of information are kept apart:

- entity identity: a DNA interval (GRCh38, 0-based half-open, as the result manifests state), a
  gene, a protein, a complex, or a measured cellular state;
- functional participation: which process an entity takes part in, under which conditions, tagged
  with responsibility categories (`CATEGORIES`), which overlap, form no hierarchy and are never
  written as block labels;
- evidence status: observed, predicted, inferred or unknown.

Every assertion names its source file, the record key inside it and the file's sha256, and quotes
the record's fields unchanged. A git-ignored cache (`data/knowledge/`, `data/cache/`) may be read
when present; what comes from it is marked "local cache, not committed", and when it is absent the
assertion says "not available in this checkout", never zero or negative.

The rules the builder enforces (`check`), and refuses to serve without:

- a perturbation's measured effect is kept apart from any direct mechanism; the committed records
  say "effect on target gene expression", so whether RNA or protein was read out stays unresolved;
- contact is never regulation, an RNA change is never a protein change, and a coherent path is never
  validation: a chain is an ordered list of assertion ids and derives no claim;
- the v3 stated-context judge stays authoritative: other-cell evidence is shown and never counts as
  support in the claimed cell; historical C4 readings are labelled apart from v3 verdicts;
- set-valued candidates stay set-valued, and an underpowered observation is not negative evidence;
- a feedback loop is a set of assertions with no execution order.

    uv run python -m genomeos.response_map globin_k562            # the payload, as JSON
    uv run python -m genomeos.response_map globin_k562 --summary  # a short reading
"""

from __future__ import annotations

import copy
import csv
import gzip
import hashlib
import json
import re
import sys
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from genomeos import evidence_discovery as evd
from genomeos.attribution.correctness import NOT_ASSESSED, norm_cell, stated
from genomeos.attribution.measured import DECREASE, INCREASE, NULL_INCONCLUSIVE, NULL_INFORMATIVE
from genomeos.genome import hic_contact

EXAMPLES = ("globin_k562",)
VIEW = "response_map"
INCREMENT = 1

HEADER = (
    "A read-only map of existing evidence and of where it stops, for one example. "
    "It adds no biological validation or prediction."
)

ENTITY_KINDS = ("dna_interval", "gene", "rna", "protein", "complex", "cell_state")
RELATIONS = (
    "encodes",
    "binds",
    "physically_contacts",
    "changes_accessibility",
    "changes_measured_rna",
    "changes_measured_protein",
    "participates_in_process",
    "associated_with_state",
)
#: a relation the source leaves open between named alternatives: `relation_candidates` lists them
UNRESOLVED = "unresolved"
STATUSES = ("observed", "predicted", "inferred", "unknown")
CATEGORIES = {
    "detect_changes": "detect changes",
    "relay_and_combine_signals": "relay and combine signals",
    "regulate_gene_responses": "regulate gene responses",
    "control_dna_access_and_organisation": "control DNA access and organisation",
    "produce_and_process_molecules": "produce and process molecules",
    "carry_out_cellular_functions": "carry out cellular functions",
    "maintain_limit_or_change_responses": "maintain, limit or change responses",
}
NOT_RECORDED = "not recorded"
UNAVAILABLE = "not available in this checkout"
LOCAL = "local cache, not committed"
COMMITTED = "committed"

CHAIN_NOTE = (
    "An ordered list of referenced assertions. A path through several relations is not a proven "
    "causal chain, and nothing is claimed beyond its steps."
)
CYCLE_NOTE = (
    "A feedback set: the assertions whose subject and object both lie in one strongly connected set "
    "of entities. It is listed without an execution order and claims nothing beyond its assertions."
)
OTHER_CELL_NOTE = (
    "Other-cell evidence is shown beside the claim and never counts as support in the cell the claim "
    "states; the v3 stated-context verdict is authoritative."
)
UNDERPOWERED_NOTE = "underpowered: not negative evidence, and no conclusion either way"
SET_NOTE = "set-valued: an admitted registry candidate is not an established cause"
SIGNED_BASIS = (
    "the record carries no significance field; genomeos/attribution/crispri_direction.py keeps only signed "
    "pairs (Significant TRUE, EffectSize non-zero, ValidConnection TRUE)"
)
SPLIT_BASIS = "genomeos/attribution/crispri_direction.py fills both sections from the held-out arm"
READOUT_NOTE = (
    "the benchmark records an effect on 'target gene expression'; which molecule the screen read out "
    "(RNA or protein) is not recorded in the committed records, so the relation stays unresolved"
)

# --- the example's sources -------------------------------------------------------------------------
LOCUS = "chr11"
LOCUS_GENES = ("HBB", "HBD", "HBG1", "HBG2", "HBE1")  # the beta-like globin genes the records measure
CELL = "K562"
CRISPRI_DIRECTION = "data/results/crispri_direction.json"
DISCOVERY = str(evd.SOURCE)
V3 = "data/results/attribution_correctness_v3.json"
PROTEOME = "genomeos/lib/data/proteome.json.gz"
ERYTHROCYTE = "data/organisms/human/erythrocyte.bio"
LOCI = "data/results/loci_benchmark.json"
LOCI_ENTRY = "HBB_LCR"
# local caches, git-ignored
CRISPRI_FILES = {
    "training": "data/knowledge/crispri/EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz",
    "heldout": "data/knowledge/crispri/EPCrisprBenchmark_combined_data.heldout_5_cell_types.GRCh38.tsv.gz",
}
CRISPRI_ACCESSION = "EngreitzLab/CRISPR_comparison EPCrisprBenchmark (Gschwind et al. 2025)"
HIC = str(hic_contact.KNOWLEDGE / f"{CELL}_{hic_contact.MATRICES[CELL]}_{hic_contact.RESOLUTION}.json")
HIC_ACCESSION = f"4D Nucleome {hic_contact.MATRICES[CELL]} (K562 in-situ Hi-C, merged replicates)"
COMPILED = "data/knowledge/compiled/noncoding_chr11.bio"
REACTOME = "data/knowledge/ReactomePathways.txt"
HALF_LIVES = "data/cache/rates/schofield2018_TableS2_halflives.xlsx"
HALF_LIFE_SHEET = "Table S2_K562"

WHY = {
    "chosen_for": (
        "traceability rather than completeness: committed CRISPRi decreases of three globin genes in "
        "K562, a model that named a different gene for the same interval, a compiled claim stated in "
        "another cell, and a clear point where the explanation stops"
    ),
    "plan": "docs/ROADMAP.md, row '↳ MAP' (commit 304a1c8), which records the architecture review",
    "examined_before": [
        {
            "result": "crispri_direction",
            "path": CRISPRI_DIRECTION,
            "key": "what_it_named_instead.rows",
            "what": "the Reilly HBG1 decreases and the K562_DC_TAP HBD decreases, with the gene the model "
            "named instead",
        },
        {
            "result": "discovery_review",
            "path": DISCOVERY,
            "key": "changed_observations; prior_exposure",
            "what": "the HBE1, HBG2, HBB and HBD observations, their C4 attachment and audit A's "
            "candidates; the file lists its own prior exposure",
        },
        {
            "result": "attribution_correctness_v3",
            "path": V3,
            "key": "not_assessed",
            "what": "EH38E2941908 -> HBE1 stated in HT1080, not assessed in this context, with the K562 "
            "decrease kept as other-cell evidence",
        },
    ],
    "consequence": (
        "every outcome shown here has been read before: this view is not a fresh test, and it selects "
        "records by a stated rule rather than by their outcomes"
    ),
    "selection_rule": (
        f"every committed record whose measured gene is one of {', '.join(LOCUS_GENES)} on {LOCUS}, "
        f"in {CELL} for measurements and in any stated cell for compiled claims; the protein records "
        "of the same genes; the committed rule that assembles a globin complex; the loci benchmark's "
        f"{LOCI_ENTRY} entry"
    ),
}


class MissingError(FileNotFoundError):
    """A committed source is not in this checkout."""


class RefusedError(ValueError):
    """The map would break one of its rules, so it serves nothing rather than a wrong assertion."""


# --- reading ---------------------------------------------------------------------------------------
class _Sources:
    """Files under `root`, read once each, with their sha256 and whether they are committed."""

    def __init__(self, root: Path) -> None:
        self.root = Path(root)
        self.table: dict[str, dict[str, Any]] = {}
        self._bytes: dict[str, bytes | None] = {}

    @staticmethod
    def committed(path: str) -> bool:
        return not path.startswith(("data/knowledge/", "data/cache/"))

    def raw(self, path: str) -> bytes | None:
        if path not in self._bytes:
            p = self.root / path
            data = p.read_bytes() if p.is_file() else None
            if data is None and self.committed(path):
                raise MissingError(f"{path} is not in this checkout")
            self._bytes[path] = data
            self.table[path] = {
                "path": path,
                "kind": COMMITTED if self.committed(path) else LOCAL,
                "available": data is not None,
                "sha256": hashlib.sha256(data).hexdigest() if data is not None else None,
                "bytes": len(data) if data is not None else None,
            }
            if data is None:
                self.table[path]["note"] = UNAVAILABLE
        return self._bytes[path]

    def json(self, path: str) -> Any:
        data = self.raw(path)
        if data is None:
            return None
        return json.loads(gzip.decompress(data) if path.endswith(".gz") else data)

    def lines(self, path: str) -> list[str] | None:
        data = self.raw(path)
        return None if data is None else data.decode().splitlines()

    def ref(self, path: str, key: str, record_path: list[Any], **extra: Any) -> dict[str, Any]:
        """The reference an assertion carries: file, record key, the path to it, sha256, kind."""
        self.raw(path)
        s = self.table[path]
        out = {
            "path": path,
            "record_key": key,
            "record_path": record_path,
            "sha256": s["sha256"],
            "kind": s["kind"],
        }
        if not s["available"]:
            out["availability"] = UNAVAILABLE
        return {**out, **extra}


def _pick(record: dict[str, Any], keys: Iterable[str]) -> dict[str, Any]:
    return {k: copy.deepcopy(record[k]) for k in keys if k in record}


def _iv(chrom: str, start: int, end: int) -> str:
    return f"interval:{chrom}:{start}-{end}"


class _Map:
    """Entities, assertions and participation records, added once each."""

    def __init__(self, assembly: dict[str, Any]) -> None:
        self.assembly = assembly
        self.entities: dict[str, dict[str, Any]] = {}
        self.assertions: list[dict[str, Any]] = []
        self.participation: list[dict[str, Any]] = []

    def entity(self, eid: str, kind: str, **fields: Any) -> str:
        if kind not in ENTITY_KINDS:
            raise RefusedError(f"entity kind {kind!r} is not one of {ENTITY_KINDS}")
        if eid not in self.entities:
            self.entities[eid] = {"id": eid, "kind": kind, **fields}
        else:
            for k, v in fields.items():  # a later source may add an alias, never replace a field
                self.entities[eid].setdefault(k, v)
        return eid

    def interval(self, chrom: str, start: int, end: int, basis: str, **fields: Any) -> str:
        return self.entity(
            _iv(chrom, start, end),
            "dna_interval",
            chrom=chrom,
            start=start,
            end=end,
            assembly=self.assembly["assembly"],
            coordinates=self.assembly["coordinates"],
            assembly_basis=basis,
            **fields,
        )

    def gene(self, symbol: str) -> str:
        return self.entity(f"gene:{symbol}", "gene", symbol=symbol)

    def add(self, a: dict[str, Any]) -> dict[str, Any]:
        base = {
            "kind": "evidence",
            "relation_candidates": [],
            "relation_note": None,
            "entities": [],
            "context": {"cell": NOT_RECORDED, "condition": NOT_RECORDED, "perturbation": NOT_RECORDED},
            "measurement": None,
            "time_window": NOT_RECORDED,
            "uncertainty": [],
            "conflicts": [],
            "unresolved_candidates": None,
        }
        a = {**base, **a}
        a["entities"] = sorted(
            set(a["entities"]) | {a["subject"]} | ({a["object"]} if a.get("object") else set())
        )
        self.assertions.append(a)
        return a

    def participate(self, p: dict[str, Any]) -> None:
        self.participation.append(p)


def _assembly(src: _Sources) -> dict[str, Any]:
    """GRCh38 and 0-based half-open, as the result manifests state; refuse if they disagree."""
    found = []
    for path in (DISCOVERY, V3):
        m = src.json(path)["result_manifest"]
        found.append((m.get("assembly"), json.dumps(m.get("coordinates"), sort_keys=True), path))
    if len({(a, c) for a, c, _ in found}) != 1 or found[0][0] is None:
        raise RefusedError(f"the result manifests disagree on assembly or coordinates: {found}")
    return {
        "assembly": found[0][0],
        "coordinates": json.loads(found[0][1]),
        "from": [f"{p}#result_manifest" for _, _, p in found],
    }


# --- perturbations ---------------------------------------------------------------------------------
def _expression_relation() -> dict[str, Any]:
    return {
        "relation": UNRESOLVED,
        "relation_candidates": ["changes_measured_rna", "changes_measured_protein"],
        "relation_note": READOUT_NOTE,
    }


def _perturbation_context(dataset: str | None) -> dict[str, Any]:
    return {
        "cell": CELL,
        "condition": NOT_RECORDED,
        "perturbation": "CRISPRi of the tested interval" + (f" ({dataset})" if dataset else ""),
    }


def _crispri_direction(src: _Sources, m: _Map) -> None:
    d = src.json(CRISPRI_DIRECTION)
    basis = (
        f"{CRISPRI_DIRECTION} carries no result manifest; its evidence field names the ENCODE CRISPRi "
        f"benchmark, which {DISCOVERY}'s manifest pins as {m.assembly['assembly']}, half-open"
    )
    for section in ("what_it_named_instead", "answered"):
        rows = d[section]["rows"] if section == "what_it_named_instead" else d[section]
        for i, r in enumerate(rows):
            if r.get("cell") != CELL or r.get("chrom") != LOCUS or r.get("gene") not in LOCUS_GENES:
                continue
            s, e = r["element"]
            iv = m.interval(LOCUS, s, e, basis)
            key = (
                f"{section}[cell={r['cell']},chrom={r['chrom']},element={s}-{e},"
                f"gene={r['gene']},dataset={r['dataset']}]"
            )
            path = [section, "rows", i] if section == "what_it_named_instead" else [section, i]
            direction = {"down": "decrease", "up": "increase"}.get(r.get("measured_sign"), NOT_RECORDED)
            conflicts = []
            named = r.get("model_named_instead")
            if named:
                conflicts.append(
                    {
                        "kind": "model_named_instead",
                        "detail": f"for this interval the model named {', '.join(named)}, not {r['gene']}",
                        "model_named": list(named),
                        "status": "predicted",
                        "resolution": "unresolved",
                        "source": {
                            "path": CRISPRI_DIRECTION,
                            "record_key": key,
                            "field": "model_named_instead",
                        },
                    }
                )
            m.add(
                {
                    "id": f"cd|{r['chrom']}:{s}-{e}|{r['gene']}|{r['dataset']}",
                    **_expression_relation(),
                    "subject": iv,
                    "object": m.gene(r["gene"]),
                    "status": "observed",
                    "basis": "measured perturbation effect (CRISPRi screen, ENCODE benchmark)",
                    "context": _perturbation_context(r["dataset"]),
                    "measurement": {
                        "quantity": "EffectSize: signed effect on target gene expression (the benchmark's "
                        "column)",
                        "value": r["measured"],
                        "direction": direction,
                        "unit": NOT_RECORDED,
                        "significant": True,
                        "significant_basis": SIGNED_BASIS,
                    },
                    "split": "heldout",
                    "split_basis": SPLIT_BASIS,
                    "mechanism": "not established: the change after perturbation may be indirect",
                    "conflicts": conflicts,
                    "source": src.ref(CRISPRI_DIRECTION, key, path),
                    "quoted": copy.deepcopy(r),
                    "assembly_basis": basis,
                    "observation_key": f"{r['chrom']}:{s}-{e}|{r['gene']}|{r['dataset']}|{CELL}",
                }
            )


def _discovery(src: _Sources, m: _Map) -> None:
    d = src.json(DISCOVERY)
    basis = f"{DISCOVERY}#result_manifest"
    for i, r in enumerate(d["changed_observations"]):
        iv = evd.parse_id(r["id"])
        if r["cell"] != CELL or iv["chrom"] != LOCUS or iv["gene"] not in LOCUS_GENES:
            continue
        eid = m.interval(LOCUS, iv["start"], iv["end"], basis)
        outcome = r["outcome"]
        measurement: dict[str, Any] = {
            "quantity": "the benchmark's outcome class for the observation",
            "outcome": outcome,
            "value": "not recorded in this record",
        }
        extra: dict[str, Any] = {}
        if outcome in (DECREASE, INCREASE):
            status = "observed"
            measurement["direction"] = "decrease" if outcome == DECREASE else "increase"
            measurement["significant"] = True
        elif outcome == NULL_INCONCLUSIVE:
            status = "unknown"
            measurement["direction"] = None
            measurement["reading"] = UNDERPOWERED_NOTE
            extra["negative_evidence"] = False
        elif outcome == NULL_INFORMATIVE:
            status = "observed"
            measurement["direction"] = "no change detected at the screen's power"
        else:
            status = "unknown"
            measurement["direction"] = None
        cands = [c["id"] for c in r["candidates"]]
        attached = [{**copy.deepcopy(a), "context": evd.claim_context(a, r["cell"])} for a in r["attached"]]
        m.add(
            {
                "id": f"dr|{r['id']}",
                **_expression_relation(),
                "subject": eid,
                "object": m.gene(iv["gene"]),
                "status": status,
                "basis": "measured perturbation outcome (CRISPRi screen, ENCODE benchmark)",
                "context": _perturbation_context(r["study"]),
                "measurement": measurement,
                "split": r["split"],
                "mechanism": "not established: the change after perturbation may be indirect",
                "unresolved_candidates": {
                    "kind": f"admitted registry candidates under audit A's discovery-only policy "
                    f"({evd.NEW_RULE})",
                    "set": cands,
                    "fraction_covered": {c["id"]: c["fraction_covered"] for c in r["candidates"]},
                    "union_coverage": dict(r["union_coverage"]),
                    "resolution": dict(r["resolution"]),
                    "note": SET_NOTE,
                },
                "historical_c4": {
                    "reading": r["verdict_current_rule"],
                    "rule": evd.VERDICT_RULE,
                    "label": evd.VERDICT_NOTE,
                },
                "discovery_policy": {
                    "status": r["status_under_new_rule"],
                    "rule": f"{evd.NEW_RULE}: audit A's policy, discovery only; it changes no matcher",
                    "attached_claims": attached,
                    "stated_cell_agreement": dict(r["stated_cell_agreement"]),
                },
                "source": src.ref(
                    DISCOVERY, f"changed_observations[id={r['id']}]", ["changed_observations", i]
                ),
                "quoted": copy.deepcopy(r),
                "assembly_basis": basis,
                "observation_key": f"{LOCUS}:{iv['start']}-{iv['end']}|{iv['gene']}|{r['study']}|{r['cell']}",
                **extra,
            }
        )


def _participation_of_effects(m: _Map) -> None:
    """An interval whose perturbation significantly moved a gene: a candidate role, inferred."""
    for a in m.assertions:
        meas = a.get("measurement") or {}
        if a["kind"] != "evidence" or a["status"] != "observed" or not meas.get("significant"):
            continue
        gene = a["object"].split(":", 1)[1]
        m.participate(
            {
                "id": f"P|{a['id']}",
                "entity": a["subject"],
                "process": f"the measured expression of {gene} in {CELL}",
                "conditions": f"{CELL}; CRISPRi of the interval",
                "categories": ["regulate_gene_responses"],
                "status": "inferred",
                "from_assertions": [a["id"]],
                "alternatives": [
                    f"acts on the regulation of {gene} directly",
                    "acts through another gene or process (an indirect effect)",
                ],
                "note": "inferred from a measured perturbation effect; the alternatives are unresolved",
            }
        )


# --- the claim stated in another cell --------------------------------------------------------------
def _local_crispri(src: _Sources) -> dict[str, list[dict[str, Any]]] | None:
    """The local benchmark rows at the locus, by split; None when the cache is absent."""
    out: dict[str, list[dict[str, Any]]] = {}
    for split, path in CRISPRI_FILES.items():
        data = src.raw(path)
        if data is None:
            return None
        text = gzip.decompress(data).decode()
        rows = []
        for n, row in enumerate(csv.DictReader(text.splitlines(), delimiter="\t"), start=2):
            if row["chrom"] == LOCUS and row["measuredGeneSymbol"] in LOCUS_GENES and row["CellType"] == CELL:
                rows.append({**row, "_line": n, "_path": path})
        out[split] = rows
    return out


def _row_key(row: dict[str, Any]) -> str:
    return (
        f"{row['chrom']}:{row['chromStart']}-{row['chromEnd']}|{row['measuredGeneSymbol']}|{row['Dataset']}"
    )


_ROW_FIELDS = (
    "chrom",
    "chromStart",
    "chromEnd",
    "measuredGeneSymbol",
    "Dataset",
    "EffectSize",
    "Significant",
    "PowerAtEffectSize20",
    "startTSS",
    "endTSS",
    "CellType",
    "Reference",
)


def _compiled_rules(src: _Sources, elements: set[str]) -> dict[str, tuple[int, str]] | None:
    """The compiled program's rule line per element (not the `_measured` layer); None when absent."""
    lines = src.lines(COMPILED)
    if lines is None:
        return None
    out: dict[str, tuple[int, str]] = {}
    for n, line in enumerate(lines, start=1):
        if line.startswith("rule "):
            name = line.split(" ", 2)[1]
            if name in elements and name not in out:
                out[name] = (n, line)
    return out


def _claims(src: _Sources, m: _Map, local: dict[str, list[dict[str, Any]]] | None) -> None:
    d = src.json(V3)
    basis = f"{V3}#result_manifest"
    picked = []
    for lst in ("judged", "not_assessed"):
        for i, r in enumerate(d[lst]):
            if r.get("gene") in LOCUS_GENES and str(r.get("locus", "")).startswith(f"{LOCUS}:"):
                picked.append((lst, i, r))
    pin = next((x["sha256"] for x in d["result_manifest"]["inputs"] if x["path"] == COMPILED), None)
    rules = _compiled_rules(src, {r["element"] for _, _, r in picked})
    for lst, i, r in picked:
        chrom, span = r["locus"].split(":")
        s, e = (int(x) for x in span.split("-"))
        iv = m.interval(chrom, s, e, basis, registry_ids=[r["element"]])
        key = f"{lst}[element={r['element']},axis={r['axis']},gene={r['gene']},cell={r['cell']}]"
        if rules is None:
            status, rule_ref = "unknown", {"availability": UNAVAILABLE, "path": COMPILED}
            status_basis = (
                "the claim's evidence kind is recorded in the compiled program, not in this checkout"
            )
        elif r["element"] in rules:
            n, line = rules[r["element"]]
            kind = re.search(r"evidence: (\w+)", line)
            status = kind.group(1) if kind and kind.group(1) in STATUSES else "unknown"
            rule_ref = src.ref(COMPILED, f"line {n}: rule {r['element']}", ["line", n], quoted=line)
            rule_ref["matches_v3_manifest_pin"] = rule_ref["sha256"] == pin
            status_basis = "the compiled program's rule for the element states this evidence kind"
        else:
            status, rule_ref = "unknown", {"path": COMPILED, "note": "no rule for the element in the program"}
            status_basis = "no compiled rule found for the element"
        other = []
        for x in r.get("cross_cell") or []:
            item = {
                **copy.deepcopy(x),
                "counts_as_support_in_stated_cell": False,
                "note": OTHER_CELL_NOTE,
            }
            item["candidate_rows"] = _cross_cell_rows(x, (chrom, s, e), local, m)
            other.append(item)
        m.add(
            {
                "id": f"v3|{r['element']}|{r['axis']}|{r['gene']}|{r['cell']}",
                "kind": "claim",
                "relation": None,
                "claim": {
                    "element": r["element"],
                    "axis": r["axis"],
                    "value": r["value"],
                    "gene": r["gene"],
                    "stated_cell": r["cell"],
                },
                "subject": iv,
                "object": m.gene(r["gene"]),
                "status": status,
                "basis": status_basis,
                "status_from": LOCAL if rules is not None else UNAVAILABLE,
                "compiled_rule": rule_ref,
                "context": {
                    "cell": r["cell"],
                    "condition": NOT_RECORDED,
                    "perturbation": "none: a compiled claim",
                },
                "verdict": {
                    "rule": r["rule"],
                    "verdict": r["verdict"],
                    "reason": r["reason"],
                    "detail": r["detail"],
                    "authoritative": True,
                    "label": "v3 stated-context verdict",
                },
                "support_in_stated_cell": [],
                "other_cell_evidence": other,
                "source": src.ref(V3, key, [lst, i]),
                "quoted": copy.deepcopy(r),
                "assembly_basis": basis,
            }
        )


def _cross_cell_rows(
    ref: dict[str, Any], element: tuple[str, int, int], local: dict[str, list[dict[str, Any]]] | None, m: _Map
) -> dict[str, Any]:
    """The local benchmark rows a v3 cross-cell reference could name: set-valued, never picked."""
    if local is None:
        return {"availability": UNAVAILABLE, "set": None}
    study = str(ref.get("source", "")).removeprefix("crispri:")
    want = {"crispri_decrease": "decrease", "crispri_increase": "increase"}.get(ref.get("kind"))
    chrom, s, e = element
    rows = []
    for row in local.get(ref.get("split"), []):
        if (
            row["Dataset"] != study
            or row["measuredGeneSymbol"] != ref.get("gene")
            or row["CellType"] != ref.get("cell")
        ):
            continue
        if row["Significant"].upper() != "TRUE":
            continue
        if want and ("decrease" if float(row["EffectSize"]) < 0 else "increase") != want:
            continue
        rs, re_ = int(row["chromStart"]), int(row["chromEnd"])
        rows.append(
            {
                "key": _row_key(row),
                "source": {
                    "path": row["_path"],
                    "record_key": _row_key(row),
                    "line": row["_line"],
                    "kind": LOCAL,
                },
                "overlap_bp_with_claim_element": max(0, min(e, re_) - max(s, rs)),
                "assertion": next(
                    (
                        a["id"]
                        for a in m.assertions
                        if a.get("observation_key", "").startswith(_row_key(row) + "|")
                    ),
                    None,
                ),
                "quoted": _pick(row, _ROW_FIELDS),
            }
        )
    return {
        "availability": LOCAL,
        "set": rows,
        "note": (
            "the v3 reference names study, gene, cell and split only, so every local row matching those "
            f"is listed; {len(rows)} match and none is picked"
        ),
    }


# --- proteins and complexes ------------------------------------------------------------------------
def _reactome_names(src: _Sources, ids: set[str]) -> dict[str, tuple[int, str]] | None:
    lines = src.lines(REACTOME)
    if lines is None:
        return None
    out = {}
    for n, line in enumerate(lines, start=1):
        parts = line.split("\t")
        if parts and parts[0] in ids:
            out[parts[0]] = (n, parts[1])
    return out


def _proteins(src: _Sources, m: _Map) -> dict[str, dict[str, Any]]:
    """encodes, function and Reactome participation, and STRING partner sets; returns what each states."""
    d = src.json(PROTEOME)
    labels = d["evidence"]
    annotation = "annotation, not measured in K562"
    ctx = {
        "cell": "not stated by the source (general annotation)",
        "condition": NOT_RECORDED,
        "perturbation": "none",
    }
    pathways = {p for g in LOCUS_GENES for p in d["proteins"][g]["pathways"]}
    names = _reactome_names(src, pathways)
    stated_here: dict[str, dict[str, Any]] = {}
    for g in LOCUS_GENES:
        rec = d["proteins"][g]
        stated_here[g] = {
            "function": rec["function"],
            "pathways": list(rec["pathways"]),
            "pathway_names": None if names is None else [names[p][1] for p in rec["pathways"] if p in names],
            "oxygen_domains": [x for x in rec["domains"] if "oxygen" in x.lower()],
        }
        prot = m.entity(f"protein:{g}", "protein", symbol=g, accession=rec["accession"], name=rec["name"])
        base = {"context": ctx, "status": "inferred", "time_window": "not applicable (annotation)"}
        m.add(
            {
                **base,
                "id": f"prot|{g}|encodes",
                "relation": "encodes",
                "subject": m.gene(g),
                "object": prot,
                "basis": f"{labels['accession']}; {annotation}",
                "source": src.ref(PROTEOME, f"proteins.{g}", ["proteins", g]),
                "quoted": _pick(rec, ("accession", "name", "symbol_match", "chrom")),
            }
        )
        if rec["function"]:
            fid = f"prot|{g}|function"
            m.add(
                {
                    **base,
                    "id": fid,
                    "relation": "participates_in_process",
                    "subject": prot,
                    "object": None,
                    "process": {"text": rec["function"], "source": "UniProt function text, verbatim"},
                    "basis": f"{labels['function']}; {annotation}",
                    "source": src.ref(PROTEOME, f"proteins.{g}.function", ["proteins", g]),
                    "quoted": _pick(rec, ("function",)),
                }
            )
            m.participate(
                {
                    "id": f"P|{fid}",
                    "entity": prot,
                    "process": rec["function"],
                    "conditions": "not stated by the source; not measured in K562",
                    "categories": ["carry_out_cellular_functions"],
                    "status": "inferred",
                    "from_assertions": [fid],
                    "alternatives": [],
                    "note": "a curated annotation, transferred; no committed record measures it in K562",
                }
            )
        for pw in rec["pathways"]:
            name = None if names is None else names.get(pw)
            m.add(
                {
                    **base,
                    "id": f"prot|{g}|reactome|{pw}",
                    "relation": "participates_in_process",
                    "subject": prot,
                    "object": None,
                    "process": {
                        "id": pw,
                        "name": name[1] if name else (UNAVAILABLE if names is None else NOT_RECORDED),
                        "name_source": (
                            {
                                "path": REACTOME,
                                "line": name[0],
                                "kind": LOCAL,
                                "sha256": src.table[REACTOME]["sha256"],
                            }
                            if name
                            else {
                                "path": REACTOME,
                                "availability": UNAVAILABLE if names is None else "id not listed",
                            }
                        ),
                    },
                    "basis": f"{labels['pathways']}; {annotation}",
                    "source": src.ref(PROTEOME, f"proteins.{g}.pathways[{pw}]", ["proteins", g]),
                    "quoted": {"pathways": list(rec["pathways"])},
                }
            )
        if rec["partners"]:
            partners = [m.entity(f"protein:{p}", "protein", symbol=p) for p in rec["partners"]]
            m.add(
                {
                    **base,
                    "id": f"prot|{g}|partners",
                    "relation": "binds",
                    "subject": prot,
                    "object": None,
                    "partners": partners,
                    "status": "predicted",
                    "basis": f"{labels['partners']}; a partner set, not a measured complex or its "
                    "stoichiometry",
                    "unresolved_candidates": {
                        "kind": "partner proteins",
                        "set": partners,
                        "note": "set-valued: which partners form a complex in K562 is not recorded",
                    },
                    "entities": partners,
                    "source": src.ref(PROTEOME, f"proteins.{g}.partners", ["proteins", g]),
                    "quoted": _pick(rec, ("partners",)),
                }
            )
    return stated_here


def _complexes(src: _Sources, m: _Map) -> list[str]:
    """The committed rules that assemble a globin complex, quoted line by line."""
    lines = src.lines(ERYTHROCYTE)
    assert lines is not None  # committed: raw() raised otherwise
    found = []
    for n, line in enumerate(lines, start=1):
        hit = re.match(r"rule (\S+) binds (\S+) \{(.*)\}\s*$", line.strip())
        if not hit:
            continue
        chain, cplx = hit.group(1), hit.group(2)
        decl = next((k for k, x in enumerate(lines, start=1) if x.startswith(f"protein {cplx} ")), None)
        cid = m.entity(f"complex:{cplx}", "complex", name=cplx, declared_at=f"{ERYTHROCYTE}:{decl}")
        pid = m.entity(f"protein:{chain}", "protein", symbol=chain)
        ev = re.search(r'evidence: (\w+) "([^"]*)"', hit.group(3))
        found.append(cplx)
        m.add(
            {
                "id": f"bio|{ERYTHROCYTE}|{n}",
                "relation": "binds",
                "subject": pid,
                "object": cid,
                "status": "inferred",
                "basis": (
                    f"the committed program's rule, evidence {ev.group(1)} ({ev.group(2)})"
                    if ev
                    else "the committed program's rule"
                ),
                "context": {
                    "cell": "erythrocyte (the program models a red blood cell, not K562)",
                    "condition": NOT_RECORDED,
                    "perturbation": "none",
                },
                "time_window": "not applicable (a program rule)",
                "uncertainty": [
                    f"the program names {chain} and {cplx}; it states no fetal or embryonic chain"
                ],
                "source": src.ref(ERYTHROCYTE, f"line {n}: rule {chain} binds {cplx}", ["line", n]),
                "quoted": {"line": line},
            }
        )
    return sorted(set(found))


# --- the loci benchmark's entry ---------------------------------------------------------------------
def _loci(src: _Sources, m: _Map) -> dict[str, Any]:
    d = src.json(LOCI)
    idx, rec = next((i, x) for i, x in enumerate(d["loci"]) if x["locus"] == LOCI_ENTRY)
    exp = rec["expected"]
    s, e = exp["element"]
    basis = f"{LOCI} records no assembly; the interval is read as it stands"
    iv = m.interval(exp["chrom"], s, e, basis, name=LOCI_ENTRY)
    rd = rec["readings"]
    reader = rd["reader"]
    acc = m.entity(
        f"state:dnase_accessible_{CELL}", "cell_state", name=f"DNase-seq accessible chromatin in {CELL}"
    )
    m.add(
        {
            "id": f"loci|{LOCI_ENTRY}|reader",
            "relation": "associated_with_state",
            "subject": iv,
            "object": acc,
            "status": "observed",
            "basis": reader["evidence"],
            "context": {
                "cell": CELL,
                "condition": NOT_RECORDED,
                "perturbation": "none (an unperturbed state)",
            },
            "measurement": {
                "quantity": "the reader's open_in value for K562 (unit not recorded in the record)",
                "value": reader["open_in"].get(CELL),
                "strongest_cell": reader["strongest_cell"],
                "peaks": len(reader.get(f"peaks_{CELL}") or []),
            },
            "uncertainty": ["accessibility is a state of the interval, not a change caused by it"],
            "source": src.ref(
                LOCI, f"loci[locus={LOCI_ENTRY}].readings.reader", ["loci", idx, "readings", "reader"]
            ),
            "quoted": _pick(reader, ("evidence", "strongest_cell", f"peaks_{CELL}"))
            | {"open_in": {CELL: reader["open_in"].get(CELL)}},
            "assembly_basis": basis,
        }
    )
    mpra = rd["mpra"]
    rep = m.entity(
        f"state:lentimpra_reporter_activity_{CELL}",
        "cell_state",
        name=f"lentiMPRA reporter activity in {CELL}",
    )
    m.add(
        {
            "id": f"loci|{LOCI_ENTRY}|mpra",
            "relation": "associated_with_state",
            "subject": iv,
            "object": rep,
            "status": "observed",
            "basis": mpra["evidence"],
            "context": {
                "cell": CELL,
                "condition": NOT_RECORDED,
                "perturbation": "reporter assay (sequence outside its locus)",
            },
            "measurement": {
                "quantity": "log2 RNA/DNA of reporter constructs from the window",
                "tested": mpra["tested"],
                "active": mpra["active"].get(CELL),
                "best": mpra["best"].get(CELL),
            },
            "uncertainty": ["a reporter assay outside the native locus: it names no endogenous target gene"],
            "source": src.ref(
                LOCI, f"loci[locus={LOCI_ENTRY}].readings.mpra", ["loci", idx, "readings", "mpra"]
            ),
            "quoted": _pick(mpra, ("evidence", "tested"))
            | {"active": {CELL: mpra["active"].get(CELL)}, "best": {CELL: mpra["best"].get(CELL)}},
            "assembly_basis": basis,
        }
    )
    dele = rd["deletion"]
    m.add(
        {
            "id": f"loci|{LOCI_ENTRY}|deletion",
            **_expression_relation(),
            "relation_note": "a predicted expression change on deletion; the record names no gene for K562",
            "subject": iv,
            "object": None,
            "status": "predicted",
            "basis": dele["evidence"],
            "context": {
                "cell": CELL,
                "condition": NOT_RECORDED,
                "perturbation": "in-silico deletion of each scored element",
            },
            "measurement": {
                "quantity": "summed_by_cell: per-element predicted log2 fold change of each element's own "
                "coding gene, summed over the scored elements",
                "value": (dele.get("summed_by_cell") or {}).get(CELL),
                "elements_scored": dele["elements_scored"],
            },
            "uncertainty": [
                f"the record's top target {dele['target']} comes from the '{dele['tissue']}' track, not from "
                f"{CELL}, and is not carried into {CELL} here"
            ],
            "source": src.ref(
                LOCI, f"loci[locus={LOCI_ENTRY}].readings.deletion", ["loci", idx, "readings", "deletion"]
            ),
            "quoted": _pick(dele, ("evidence", "elements_scored", "target", "tissue"))
            | {"summed_by_cell": {CELL: (dele.get("summed_by_cell") or {}).get(CELL)}},
            "assembly_basis": basis,
        }
    )
    return {"entry": LOCI_ENTRY, "top_target": dele["target"], "top_target_track": dele["tissue"]}


# --- contact, from the local cache only ----------------------------------------------------------------
def _contacts(src: _Sources, m: _Map, local: dict[str, list[dict[str, Any]]] | None) -> dict[str, Any]:
    """K562 Hi-C contact between each significant interval and its gene's TSS as the benchmark places it."""
    cache = src.json(HIC)
    tss: dict[str, set[int]] = {}
    index: dict[str, dict[str, Any]] = {}
    for rows in (local or {}).values():
        for row in rows:
            tss.setdefault(row["measuredGeneSymbol"], set()).add(int(row["startTSS"]))
            index.setdefault(_row_key(row), row)
    perturbs = [
        a
        for a in m.assertions
        if a["kind"] == "evidence"
        and a.get("observation_key")
        and (a.get("measurement") or {}).get("significant")
    ]
    for a in perturbs:
        key = a["observation_key"].rsplit("|", 1)[0]
        gene = a["object"].split(":", 1)[1]
        chrom, span = key.split("|")[0].split(":")
        s, e = (int(x) for x in span.split("-"))
        cid = f"hic|{a['observation_key']}"
        common = {
            "id": cid,
            "relation": "physically_contacts",
            "subject": a["subject"],
            "context": {"cell": CELL, "condition": NOT_RECORDED, "perturbation": "none (unperturbed cells)"},
            "uncertainty": ["contact is not regulation, and says nothing about direction or effect"],
            "paired_with": a["id"],
        }
        row = index.get(key)
        if cache is None or row is None:
            m.add(
                {
                    **common,
                    "object": m.gene(gene),
                    "status": "unknown",
                    "basis": f"K562 Hi-C contact: {UNAVAILABLE}",
                    "availability": UNAVAILABLE,
                    "measurement": None,
                    "source": {
                        "path": HIC,
                        "kind": LOCAL,
                        "availability": UNAVAILABLE,
                        "sha256": src.table[HIC]["sha256"],
                        "record_key": None,
                        "record_path": None,
                    },
                    "quoted": None,
                }
            )
            continue
        t0, t1 = int(row["startTSS"]), int(row["endTSS"])
        target = m.interval(
            chrom,
            t0,
            t1,
            "the benchmark's startTSS and endTSS (local cache)",
            name=f"{gene} TSS as the benchmark places it",
        )
        mid = (s + e) // 2
        b1, b2 = mid // hic_contact.RESOLUTION, t0 // hic_contact.RESOLUTION
        hkey = f"{hic_contact._chrom_key(chrom)}:{min(b1, b2)}:{max(b1, b2)}"
        value = cache.get(hkey)
        conflicts = []
        if len(tss.get(gene, ())) > 1:
            conflicts.append(
                {
                    "kind": "tss_position",
                    "detail": f"the local benchmark files place the {gene} TSS at "
                    + ", ".join(f"{x:,}" for x in sorted(tss[gene])),
                    "resolution": "unresolved: the contact is read at the position the matching row states",
                }
            )
        m.add(
            {
                **common,
                "object": target,
                "entities": [m.gene(gene)],
                "status": "observed" if value else "unknown",
                "basis": hic_contact.EVIDENCE,
                "availability": LOCAL,
                "measurement": (
                    {
                        "quantity": f"Hi-C contact, {value['norm']}-normalised, observed/expected, "
                        f"{value['binsize']:,} bp bins",
                        "oe": value["oe"],
                        "observed": value["observed"],
                        "expected": value["expected"],
                        "raw": value["raw"],
                        "bins": [min(b1, b2), max(b1, b2)],
                    }
                    if value
                    else None
                ),
                "conflicts": conflicts,
                "source": src.ref(HIC, hkey, [hkey], accession=HIC_ACCESSION),
                "position_source": {
                    "path": row["_path"],
                    "record_key": _row_key(row),
                    "line": row["_line"],
                    "kind": LOCAL,
                    "sha256": src.table[row["_path"]]["sha256"],
                },
                "quoted": copy.deepcopy(value),
            }
        )
    return {"tss_positions": {g: sorted(v) for g, v in sorted(tss.items())}}


# --- stop points and dynamics ---------------------------------------------------------------------------
def _stop_points(
    m: _Map,
    local: dict[str, list[dict[str, Any]]] | None,
    loci: dict[str, Any],
    complexes: list[str],
    proteins: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    by = {a["id"]: a for a in m.assertions}
    perturb = [a for a in m.assertions if a.get("observation_key")]
    sig = [a for a in perturb if (a.get("measurement") or {}).get("significant")]
    under = [a for a in perturb if a.get("negative_evidence") is False]
    groups: dict[tuple[str, str], list[str]] = {}
    for a in sig:
        groups.setdefault((a["subject"], a["measurement"]["direction"]), []).append(a["id"])
    shared = {k: v for k, v in groups.items() if len({by[i]["object"] for i in v}) > 1}
    named = [a["id"] for a in m.assertions if any(c["kind"] == "model_named_instead" for c in a["conflicts"])]
    claims = [a["id"] for a in m.assertions if a["kind"] == "claim"]
    contacts = [a["id"] for a in m.assertions if a["relation"] == "physically_contacts"]

    def says_oxygen(g: str) -> bool | None:
        x = proteins[g]
        if "oxygen" in x["function"].lower():
            return True
        if x["pathway_names"] is None:
            return None  # the pathway names are not in this checkout: not assessable, not "no"
        return any("oxygen" in n.lower() for n in x["pathway_names"])

    fetal = [g for g in ("HBG1", "HBG2", "HBE1") if says_oxygen(g) is False]
    fetal_unknown = [g for g in ("HBG1", "HBG2", "HBE1") if says_oxygen(g) is None]
    oxygen_detail = "; ".join(
        f"{g}: Reactome {', '.join(proteins[g]['pathways']) or 'none'}"
        + (f" ({'; '.join(proteins[g]['pathway_names'])})" if proteins[g]["pathway_names"] else "")
        + f"; function text '{proteins[g]['function']}'; oxygen named only in the InterPro domain "
        + ", ".join(proteins[g]["oxygen_domains"])
        for g in fetal
    )
    if fetal and proteins[fetal[0]]["pathway_names"] is not None:
        oxygen_detail += f" (Reactome names from {REACTOME}, {LOCAL})"
    if fetal_unknown:
        oxygen_detail += (
            f"; for {', '.join(fetal_unknown)} the Reactome names are {UNAVAILABLE}, so this is not assessed"
        )
    out = [
        {
            "id": "readout_not_recorded",
            "what_is_missing": "a protein measurement after these perturbations",
            "detail": READOUT_NOTE + ". No committed record measures globin protein or hemoglobin after them",
            "assertions": [a["id"] for a in perturb],
        },
        {
            "id": "no_phenotype",
            "what_is_missing": "a measured cellular phenotype after these perturbations",
            "detail": "no committed record measures a cellular state (hemoglobin content, "
            "differentiation, growth) "
            "after CRISPRi of these intervals",
            "assertions": [],
        },
        {
            "id": "no_time_course",
            "what_is_missing": "a time course",
            "detail": "each screen reports one sampling point, and its timing is not recorded in the "
            "committed records",
            "assertions": [a["id"] for a in sig],
        },
        {
            "id": "contact_not_committed",
            "what_is_missing": "committed contact",
            "detail": "K562 Hi-C contact exists only as a local cache; in a clean checkout it is not "
            "available. "
            "Where shown, contact is not regulation",
            "assertions": contacts,
        },
        {
            "id": "one_interval_several_genes",
            "what_is_missing": "which gene an interval acts on directly",
            "detail": "; ".join(
                f"{k[0].removeprefix('interval:')} ({k[1]}): "
                + ", ".join(sorted({by[i]["object"].removeprefix("gene:") for i in v}))
                for k, v in sorted(shared.items())
            )
            + ". Unresolved: direct action on each, on one, or through another gene",
            "assertions": sorted(i for v in shared.values() for i in v),
        },
        {
            "id": "model_named_another_gene",
            "what_is_missing": "agreement between the measurement and the model on the target",
            "detail": "for these intervals the model named a different gene; both are shown and neither "
            "overrides the other",
            "assertions": named,
        },
        {
            "id": "stated_cell_mismatch",
            "what_is_missing": "an observation in the cell each compiled claim states",
            "detail": "the claims state cells other than K562; v3 holds them not assessed in this "
            "context, and the K562 "
            "decreases stay other-cell evidence",
            "assertions": claims,
        },
        {
            "id": "underpowered",
            "what_is_missing": "a powered test of these interval-gene pairs",
            "detail": f"{len(under)} observations are underpowered: {UNDERPOWERED_NOTE}",
            "assertions": [a["id"] for a in under],
        },
        {
            "id": "oxygen_transport_not_stated",
            "what_is_missing": "a committed process annotation of oxygen transport for the fetal and "
            "embryonic chains",
            "detail": oxygen_detail + ". A domain name is not a process",
            "genes": fetal,
            "assertions": [
                a["id"] for a in m.assertions if a["id"].startswith(tuple(f"prot|{g}|" for g in fetal))
            ],
        },
        {
            "id": "complex_adult_only",
            "what_is_missing": "a committed complex rule for the fetal or embryonic chains",
            "detail": f"the committed complex rules assemble {', '.join(complexes)} from HBA and HBB only "
            "(a red blood "
            "cell program); they are not extended to HBG1, HBG2 or HBE1. The proteome's partner sets are "
            "STRING "
            "predictions",
            "assertions": [a["id"] for a in m.assertions if a["id"].startswith("bio|")],
        },
        {
            "id": "deletion_track_not_K562",
            "what_is_missing": "a K562 target for the predicted deletion effect",
            "detail": f"the loci benchmark's deletion reading names {loci['top_target']} from the "
            f"'{loci['top_target_track']}' track; for K562 it holds only a sum over elements, naming no gene",
            "assertions": [f"loci|{LOCI_ENTRY}|deletion"],
        },
    ]
    committed_keys = {a["observation_key"].rsplit("|", 1)[0] for a in perturb}
    if local is None:
        extra = {"count": None, "availability": UNAVAILABLE}
    else:
        rows = [r for rs in local.values() for r in rs if r["Significant"].upper() == "TRUE"]
        extra = {
            "count": sum(1 for r in rows if _row_key(r) not in committed_keys),
            "of": len(rows),
            "availability": LOCAL,
        }
    out.append(
        {
            "id": "local_rows_not_assembled",
            "what_is_missing": "the benchmark rows no committed result carries one by one",
            "detail": (
                f"the local benchmark cache holds {extra['count']} of {extra['of']} significant {CELL} "
                "rows for "
                f"{', '.join(LOCUS_GENES)} that no committed record lists; increment 1 does not assemble them"
                if local is not None
                else f"the local benchmark cache is {UNAVAILABLE}"
            ),
            "assertions": [],
            **extra,
        }
    )
    return out


def _half_lives(src: _Sources) -> dict[str, Any]:
    data = src.raw(HALF_LIVES)
    if data is None:
        return {"availability": UNAVAILABLE, "path": HALF_LIVES}
    try:
        import io

        import openpyxl
    except ImportError:
        return {"availability": "local cache present, but no reader for it here", "path": HALF_LIVES}
    wb = openpyxl.load_workbook(io.BytesIO(data), read_only=True)
    rows = wb[HALF_LIFE_SHEET].iter_rows(values_only=True)
    header = next(rows)
    col = header.index("mean_half_life")
    found = {r[0]: r[col] for r in rows if r and r[0] in LOCUS_GENES}
    return {
        "availability": LOCAL,
        "path": HALF_LIVES,
        "sha256": src.table[HALF_LIVES]["sha256"],
        "sheet": HALF_LIFE_SHEET,
        "mean_half_life_h": {g: found.get(g) for g in LOCUS_GENES},
        "absent": [g for g in LOCUS_GENES if g not in found],
        "note": "Schofield et al. 2018, TimeLapse-seq in K562; hours, as genomeos.attribution.bridge "
        "reads the table",
    }


def _dynamics(src: _Sources) -> dict[str, Any]:
    hl = _half_lives(src)
    na = "not assessable with current data"
    return {
        "note": (
            "A static assessment: no simulation is run. No human time course, dose response or persistence "
            "measurement exists locally for this example. A computational ordering of cells is not "
            "elapsed time."
        ),
        "half_lives": hl,
        "rows": [
            {
                "question": "response timing: how soon after the perturbation does the gene's expression "
                "change",
                "measurements_needed": "expression at several times after an acute perturbation with a "
                "known onset",
                "available": "none: the screens report one sampling point, its timing not recorded; "
                "K562 mRNA half-lives are local only (a steady-state decay rate, not a response time)",
                "assessment": na,
            },
            {
                "question": "stimulus strength: how strong a perturbation is needed",
                "measurements_needed": "a graded (dose-response) perturbation of the interval with "
                "expression read at each level",
                "available": "none",
                "assessment": na,
            },
            {
                "question": "persistence: does the change last after the perturbation is removed",
                "measurements_needed": "expression after the perturbation is withdrawn, over time",
                "available": "none",
                "assessment": na,
            },
            {
                "question": "recovery and cell-to-cell variability",
                "measurements_needed": "single-cell expression over time, before, during and after the "
                "perturbation",
                "available": "none: the committed records hold one effect per interval and gene, not a "
                "distribution over cells or time",
                "assessment": na,
            },
            {
                "question": "effects conditional on another perturbation",
                "measurements_needed": "combinatorial perturbations (this interval with another, and each "
                "alone)",
                "available": "none in the committed records",
                "assessment": na,
            },
        ],
    }


# --- structure over the assertions ---------------------------------------------------------------
def _chains(m: _Map) -> list[dict[str, Any]]:
    """Reading paths as ordered lists of assertion ids; each consecutive pair shares an entity."""
    ids = {a["id"] for a in m.assertions}
    specs = [
        (
            "chain|reilly_hbg1_proximal",
            "Reilly: an interval lowers HBG1 expression; the model named HBG2",
            ["cd|chr11:5253147-5253547|HBG1|Reilly", "prot|HBG1|encodes", "prot|HBG1|function"],
        ),
        (
            "chain|reilly_hbg1_distal",
            "Reilly: an interval lowers HBG1 expression; the model named HBE1",
            ["cd|chr11:5275847-5276247|HBG1|Reilly", "prot|HBG1|encodes", "prot|HBG1|function"],
        ),
        (
            "chain|nasser_hbe1",
            "Nasser2021: one interval lowers HBE1 expression (and HBG2's)",
            ["dr|training|chr11:5280370-5281170|HBE1|K562", "prot|HBE1|encodes", "prot|HBE1|function"],
        ),
        (
            "chain|nasser_hbg2",
            "Nasser2021: the same interval lowers HBG2 expression (and HBE1's)",
            ["dr|training|chr11:5280370-5281170|HBG2|K562", "prot|HBG2|encodes", "prot|HBG2|function"],
        ),
    ]
    out = []
    for cid, title, steps in specs:
        if not all(s in ids for s in steps):
            continue
        contact = next((a["id"] for a in m.assertions if a.get("paired_with") == steps[0]), None)
        out.append(
            {
                "id": cid,
                "title": title,
                "steps": steps,
                "beside": [contact] if contact else [],
                "derived_claims": [],
                "stops_at": ["readout_not_recorded", "no_phenotype", "no_time_course"],
                "note": CHAIN_NOTE
                + " No step says that perturbing the interval changes the protein or what it does.",
            }
        )
    return out


def _links(assertions: dict[str, dict[str, Any]], steps: list[str]) -> list[dict[str, Any]]:
    links = []
    for x, y in zip(steps, steps[1:], strict=False):
        shared = sorted(set(assertions[x]["entities"]) & set(assertions[y]["entities"]))
        if not shared:
            raise RefusedError(f"chain steps {x} and {y} share no entity")
        links.append({"from": x, "to": y, "shared_entities": shared})
    return links


def cycles(assertions: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Feedback sets among the assertions: strongly connected sets of entities, with no order."""
    edges: dict[str, set[str]] = {}
    for a in assertions:
        if a.get("kind") == "evidence" and a.get("object"):
            edges.setdefault(a["subject"], set()).add(a["object"])
            edges.setdefault(a["object"], set())
    index: dict[str, int] = {}
    low: dict[str, int] = {}
    stack: list[str] = []
    on: set[str] = set()
    comps: list[set[str]] = []

    def visit(v: str) -> None:  # Tarjan, recursive: these graphs hold tens of nodes
        index[v] = low[v] = len(index)
        stack.append(v)
        on.add(v)
        for w in sorted(edges[v]):
            if w not in index:
                visit(w)
                low[v] = min(low[v], low[w])
            elif w in on:
                low[v] = min(low[v], index[w])
        if low[v] == index[v]:
            comp = set()
            while True:
                w = stack.pop()
                on.discard(w)
                comp.add(w)
                if w == v:
                    break
            comps.append(comp)

    for v in sorted(edges):
        if v not in index:
            visit(v)
    out = []
    for comp in comps:
        members = sorted(
            a["id"]
            for a in assertions
            if a.get("kind") == "evidence" and a.get("object") in comp and a["subject"] in comp
        )
        if len(comp) > 1 or members:
            out.append(
                {
                    "id": "cycle|" + "|".join(sorted(comp)),
                    "entities": sorted(comp),
                    "assertions": members,
                    "execution_order": None,
                    "note": CYCLE_NOTE,
                }
            )
    return out


# --- the rules, enforced -----------------------------------------------------------------------------
_FORBIDDEN_CHAIN_KEYS = {"claim", "claims", "conclusion", "supports", "validated", "mechanism"}


def check(payload: dict[str, Any]) -> None:
    """Refuse a payload that breaks one of the map's rules (see the module docstring)."""
    by = {}
    for a in payload["assertions"]:
        if a["id"] in by:
            raise RefusedError(f"assertion {a['id']} appears twice")
        by[a["id"]] = a
        if a["status"] not in STATUSES:
            raise RefusedError(f"{a['id']}: status {a['status']!r} is not one of {STATUSES}")
        src = a.get("source") or {}
        if not src.get("path") or "record_key" not in src or "sha256" not in src:
            raise RefusedError(f"{a['id']}: every assertion names its source file, record key and sha256")
        if a.get("derived_from"):
            raise RefusedError(f"{a['id']}: an assertion is read from a source, never derived from others")
        if a["kind"] == "evidence":
            rel = a["relation"]
            if rel == UNRESOLVED:
                if not set(a["relation_candidates"]) <= set(RELATIONS) or len(a["relation_candidates"]) < 2:
                    raise RefusedError(f"{a['id']}: an unresolved relation names its alternatives")
            elif rel not in RELATIONS:
                raise RefusedError(f"{a['id']}: relation {rel!r} is not one of {RELATIONS}")
            if rel == "physically_contacts":
                meas = a.get("measurement") or {}
                if meas.get("direction") or a.get("mechanism") or a.get("relation_candidates"):
                    raise RefusedError(f"{a['id']}: contact carries no direction, effect or mechanism")
        elif a["kind"] == "claim":
            stated_cell = a["claim"]["stated_cell"]
            if a["verdict"].get("label") != "v3 stated-context verdict" or not a["verdict"].get(
                "authoritative"
            ):
                raise RefusedError(f"{a['id']}: a claim carries its v3 verdict as authoritative")
            for x in a["other_cell_evidence"]:
                if x.get("counts_as_support_in_stated_cell") is not False:
                    raise RefusedError(
                        f"{a['id']}: other-cell evidence never counts as support in the claimed cell"
                    )
                if stated(x.get("cell", "")) and norm_cell(x["cell"]) == norm_cell(stated_cell):
                    raise RefusedError(f"{a['id']}: evidence in the stated cell is not other-cell evidence")
            for x in a["support_in_stated_cell"]:
                if norm_cell(x.get("cell", "")) != norm_cell(stated_cell):
                    raise RefusedError(
                        f"{a['id']}: support in the stated cell comes from the stated cell only"
                    )
            if a["verdict"]["reason"] == NOT_ASSESSED and a["support_in_stated_cell"]:
                raise RefusedError(
                    f"{a['id']}: a claim not assessed in its context has no support listed there"
                )
        else:
            raise RefusedError(f"{a['id']}: kind {a['kind']!r} is neither evidence nor claim")
        meas = a.get("measurement") or {}
        if meas.get("outcome") == NULL_INCONCLUSIVE and (
            a["status"] != "unknown" or a.get("negative_evidence") is not False
        ):
            raise RefusedError(f"{a['id']}: an underpowered observation is unknown and not negative evidence")
        cands = a.get("unresolved_candidates")
        if cands and cands.get("set") is not None and not isinstance(cands["set"], list):
            raise RefusedError(f"{a['id']}: candidates stay a set")
        if "historical_c4" in a and a["historical_c4"].get("label") != evd.VERDICT_NOTE:
            raise RefusedError(f"{a['id']}: a C4 reading is labelled as historical, apart from v3")
    keys: dict[str, str] = {}
    for a in payload["assertions"]:
        k = a.get("observation_key")
        if k:
            if k in keys:
                raise RefusedError(f"observation {k} is carried by {keys[k]} and {a['id']}")
            keys[k] = a["id"]
    for eid, e in payload["entities"].items():
        if e["kind"] not in ENTITY_KINDS:
            raise RefusedError(f"entity {eid}: kind {e['kind']!r}")
        if {"label", "labels", "role", "roles", "categories"} & set(e):
            raise RefusedError(
                f"entity {eid}: roles and categories live in participation, never on an entity"
            )
    for p in payload["participation"]:
        bad = set(p["categories"]) - set(CATEGORIES)
        if bad:
            raise RefusedError(f"{p['id']}: unknown categories {sorted(bad)}")
        for aid in p["from_assertions"]:
            if aid not in by:
                raise RefusedError(f"{p['id']}: cites a missing assertion {aid}")
            if by[aid]["relation"] == "physically_contacts":
                raise RefusedError(f"{p['id']}: contact never becomes a role")
    for c in payload["chains"]:
        if _FORBIDDEN_CHAIN_KEYS & set(c) or c.get("derived_claims"):
            raise RefusedError(f"{c['id']}: a chain lists assertions and derives no claim")
        for aid in c["steps"] + c.get("beside", []):
            if aid not in by:
                raise RefusedError(f"{c['id']}: cites a missing assertion {aid}")
        if c.get("links") != _links(by, c["steps"]):
            raise RefusedError(f"{c['id']}: its links are not the shared entities of its steps")
    for cy in payload["cycles"]:
        if cy.get("execution_order") is not None:
            raise RefusedError(f"{cy['id']}: a feedback set has no execution order")
        if not set(cy["assertions"]) <= set(by):
            raise RefusedError(f"{cy['id']}: a feedback set is made only of existing assertions")


def _counts(payload: dict[str, Any]) -> dict[str, Any]:
    status: dict[str, int] = dict.fromkeys(STATUSES, 0)
    relation: dict[str, int] = {}
    kinds: dict[str, int] = {}
    source_kind: dict[str, int] = {}
    for a in payload["assertions"]:
        status[a["status"]] += 1
        r = a["relation"] or "claim"
        relation[r] = relation.get(r, 0) + 1
        kinds[a["kind"]] = kinds.get(a["kind"], 0) + 1
        k = a["source"].get("availability") or a["source"]["kind"]
        source_kind[k] = source_kind.get(k, 0) + 1
    return {
        "assertions": len(payload["assertions"]),
        "by_status": status,
        "by_relation": dict(sorted(relation.items())),
        "by_kind": kinds,
        "by_source": source_kind,
        "entities": len(payload["entities"]),
        "participation": len(payload["participation"]),
        "chains": len(payload["chains"]),
        "cycles": len(payload["cycles"]),
        "stop_points": len(payload.get("stop_points") or []),
        "unit": "assertions, each one record of one source",
    }


def assemble(
    entities: dict[str, dict[str, Any]],
    assertions: list[dict[str, Any]],
    participation: list[dict[str, Any]],
    chains: list[dict[str, Any]],
    **rest: Any,
) -> dict[str, Any]:
    """Put a payload together from its parts, link the chains, find feedback sets, and check it."""
    by = {a["id"]: a for a in assertions}
    for c in chains:
        c["links"] = _links(by, c["steps"])
    payload = {
        "entities": entities,
        "assertions": assertions,
        "participation": participation,
        "chains": chains,
        "cycles": cycles(assertions),
        **rest,
    }
    check(payload)
    payload["counts"] = _counts(payload)
    return payload


def build(root: Path, example: str = "globin_k562") -> dict[str, Any]:
    """The map of one example from the files under `root`; refuses rather than serve a broken rule."""
    if example not in EXAMPLES:
        raise KeyError(f"unknown example {example!r}; known: {', '.join(EXAMPLES)}")
    src = _Sources(Path(root))
    assembly = _assembly(src)
    m = _Map(assembly)
    _crispri_direction(src, m)
    _discovery(src, m)
    _participation_of_effects(m)
    local = _local_crispri(src)
    _claims(src, m, local)
    proteins = _proteins(src, m)
    complexes = _complexes(src, m)
    loci = _loci(src, m)
    tss = _contacts(src, m, local)
    stops = _stop_points(m, local, loci, complexes, proteins)
    dynamics = _dynamics(src)
    payload = assemble(
        m.entities,
        m.assertions,
        m.participation,
        _chains(m),
        view=VIEW,
        increment=INCREMENT,
        example=example,
        title="The beta-globin locus in K562",
        header=HEADER,
        why=WHY,
        conventions={
            **assembly,
            "not_recorded": NOT_RECORDED,
            "unavailable": UNAVAILABLE,
            "local": LOCAL,
            "other_cell": OTHER_CELL_NOTE,
            "underpowered": UNDERPOWERED_NOTE,
            "set_valued": SET_NOTE,
            "readout": READOUT_NOTE,
            "historical_c4": evd.VERDICT_NOTE,
            "chains": CHAIN_NOTE,
            "cycles": CYCLE_NOTE,
        },
        vocabulary={
            "entity_kinds": list(ENTITY_KINDS),
            "relations": list(RELATIONS),
            "statuses": list(STATUSES),
            "categories": CATEGORIES,
        },
        stop_points=stops,
        dynamics=dynamics,
        tss_positions=tss["tss_positions"],
        requests={"network": 0, "model": 0, "downloads": 0},
    )
    payload["sources"] = [src.table[p] for p in sorted(src.table)]
    return payload


_CACHE: dict[str, tuple[tuple[Any, ...], dict[str, Any]]] = {}


def _stamp(root: Path) -> tuple[Any, ...]:
    paths = [
        CRISPRI_DIRECTION,
        DISCOVERY,
        V3,
        PROTEOME,
        ERYTHROCYTE,
        LOCI,
        HIC,
        COMPILED,
        REACTOME,
        HALF_LIVES,
        *CRISPRI_FILES.values(),
    ]
    out = []
    for p in paths:
        try:
            st = (Path(root) / p).stat()
            out.append((p, st.st_mtime_ns, st.st_size))
        except FileNotFoundError:
            out.append((p, None, None))
    return tuple(out)


def view(root: Path, example: str = "globin_k562") -> dict[str, Any]:
    """`build`, read once per version of its input files."""
    key = f"{Path(root).resolve()}|{example}"
    stamp = _stamp(root)
    hit = _CACHE.get(key)
    if hit is None or hit[0] != stamp:
        hit = _CACHE[key] = (stamp, build(root, example))
    return copy.deepcopy(hit[1])


def summary(payload: dict[str, Any]) -> str:
    c = payload["counts"]
    lines = [
        f"{payload['title']} (response map, increment {payload['increment']})",
        payload["header"],
        f"{c['assertions']} assertions: " + ", ".join(f"{k} {v}" for k, v in c["by_status"].items()),
        "by source: " + ", ".join(f"{k} {v}" for k, v in sorted(c["by_source"].items())),
        f"{c['entities']} entities, {c['participation']} participation records, {c['chains']} chains, "
        f"{c['cycles']} feedback sets",
        "",
        "where the evidence stops:",
    ]
    for s in payload["stop_points"]:
        lines.append(f"  - {s['id']}: {s['what_is_missing']}")
    lines += ["", "dynamics:"]
    for r in payload["dynamics"]["rows"]:
        lines.append(f"  - {r['question'].split(':')[0]}: {r['assessment']}")
    lines += ["", "sources:"]
    for s in payload["sources"]:
        lines.append(f"  - {s['path']} [{s['kind']}] " + (s["sha256"][:12] if s["sha256"] else UNAVAILABLE))
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(prog="python -m genomeos.response_map", description=__doc__.split("\n")[0])
    ap.add_argument("example", nargs="?", default="globin_k562", choices=EXAMPLES)
    ap.add_argument("--summary", action="store_true", help="a short reading instead of the JSON payload")
    ap.add_argument("--root", default=".", help="the checkout to read (default: the current directory)")
    args = ap.parse_args(argv)
    payload = build(Path(args.root), args.example)
    sys.stdout.write(
        (summary(payload) if args.summary else json.dumps(payload, indent=1, ensure_ascii=False)) + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
