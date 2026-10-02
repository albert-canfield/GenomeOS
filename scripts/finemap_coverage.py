# SPDX-License-Identifier: AGPL-3.0-or-later
"""How many of the 1,505 CRISPRi elements are ASSESSABLE for a fine-mapped eQTL at all.

A re-distillation (registration `eqtl_crispri_frame_registration`, result `eqtl_crispri_frame`)
reported 1,262 of 1,505 compiled elements carrying a retained cis-eQTL and 96 addable loci against a
registered floor of 20. Those hits are GTEx v8 SINGLE-TISSUE SIGNIFICANT variant-gene pairs:
`genomeos.attribution.eqtl.HIT_COLUMNS` carries no PIP, no credible set and no posterior. Significant
variants travel in LD blocks, so an overlap of that kind supports "a variant in LD with an
association lies here" and not "the association is localised here".

This writer asks the prior question, and asks it the way the re-distillation's own registration
imposed one level down: before counting how many elements CARRY a fine-mapped variant, count how
many could have been SEEN to carry one. An element with no fine-mapped variant in range is
UNASSESSED, not an element without one.

EVERYTHING QUANTITATIVE IS IMPORTED, nothing is chosen here:

  * the PIP threshold is 0.5, `genomeos.attribution.executor.DAPG_PIP`, whose registered wording in
    `data/results/executor_replication.json` is "the unit's recurring value is a variant DAP-G
    fine-maps (PIP at or above 0.5) for a gene with GTEx p at or below 1e-5 and one sign across
    tissues";
  * the label for an overlap that fails the fine-mapping test is this project's own existing
    distinction between `E2R_eqtl_replication` and `E3R_eqtl_linked_replication`, the latter being
    "the other significant eQTL values, which are usually linked to a cause rather than causal";
  * the margin is `scripts/eqtl_targets.MARGIN = 500`, reached through
    `response_map_increment3_count.EQTL_MARGIN`;
  * the frame of elements is `response_map2.candidates`'s own attachment, the same call the
    re-distillation's frame was built from.

NO BASELINE AND NO DEEPENING ARM IS BUILT HERE. Both were conditional on coverage supporting a
count, and the count below says it does not.
"""

from __future__ import annotations

import argparse
import bisect
import collections
import gzip
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import response_map_increment3_count as c3  # noqa: E402

from genomeos import manifest as mf  # noqa: E402
from genomeos import response_map2 as rm2  # noqa: E402
from genomeos.attribution import cell2, eqtl  # noqa: E402
from genomeos.attribution import context_evidence as ce  # noqa: E402
from genomeos.attribution import executor as ex  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

NAME = "finemap_coverage"
ENTRY = "scripts/finemap_coverage.py"
OWN_CODE = (ENTRY,)

#: The one threshold, imported. Not chosen here and not restated as a number of this lane's own.
PIP = ex.DAPG_PIP

#: The margin, imported through increment 3's count from `scripts/eqtl_targets.MARGIN`.
MARGIN = c3.EQTL_MARGIN

#: Where the fine-mapped posterior lives at all. The per-element records on disk are what a previous
#: run happened to retain from it; the track itself is remote and is NOT fetched by this writer.
DAPG_TRACK = "https://hgdownload.soe.ucsc.edu/gbdb/hg38/gtex/eQtl/gtexDapg.bb"
DAPG_SOURCE = "experimental: GTEx v8 fine-mapped cis-eQTLs (DAP-G, UCSC gtexEqtlDapg)"

#: The label an overlap that fails the fine-mapping test gets. This project's own distinction,
#: carried in its registered wording rather than named afresh.
LD_TAGGED_LABEL = "LD-tagged association, not localised"
E3R_WORDING = "the other significant eQTL values, which are usually linked to a cause rather than causal"
E2R_WORDING = (
    "the unit's recurring value is a variant DAP-G fine-maps (PIP at or above 0.5) for a gene with "
    "GTEx p at or below 1e-5 and one sign across tissues"
)

#: The distinction this writer applies, quoted from the re-distillation's own registration.
THE_DISTINCTION = {
    "assessed": "how many of the frame's elements could have been seen at all",
    "carrying_a_hit": "how many of them carry at least one retained hit",
    "why": (
        "conflating them is how an unassessed element becomes an element with no eQTL. The first is "
        "a property of the frame and the second is a measurement"
    ),
    "quoted_from": "data/results/eqtl_crispri_frame_registration.json, the_two_figures_kept_apart",
}

PANEL = Path("data/knowledge/human_panel")
EXECUTOR = PANEL / "executor"

#: The executor files that carry a `pip` at all, which is the whole fine-mapped record on disk.
PIP_FILES = (
    "replication_assembled.json.gz",
    "assembled.json.gz",
    "two_instrument_rows.json.gz",
    "run_rows.json.gz",
    "run_rows_E1_mpra-E2_eqtl.json.gz",
    "replication_rows_E2R_eqtl_replication.json.gz",
)


# --------------------------------------------------------------------------- intervals


class Spans:
    """Merged half-open intervals per chromosome, asked whether a window touches one."""

    def __init__(self) -> None:
        self._by: dict[str, list[list[int]]] = {}
        self._starts: dict[str, list[int]] = {}

    def add(self, chrom: str, start: int, end: int) -> None:
        self._by.setdefault(chrom, []).append([start, end])

    def freeze(self) -> Spans:
        for c, rows in self._by.items():
            rows.sort()
            merged: list[list[int]] = []
            for s, e in rows:
                if merged and s <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            self._by[c] = merged
            self._starts[c] = [r[0] for r in merged]
        return self

    def touches(self, chrom: str, lo: int, hi: int) -> bool:
        rows = self._by.get(chrom)
        if not rows:
            return False
        i = bisect.bisect_right(self._starts[chrom], hi) - 1
        while i >= 0 and rows[i][1] > lo:
            if rows[i][0] < hi and rows[i][1] > lo:
                return True
            i -= 1
        return False

    @property
    def spans(self) -> int:
        return sum(len(v) for v in self._by.values())

    @property
    def bases(self) -> int:
        return sum(e - s for v in self._by.values() for s, e in v)

    @property
    def chromosomes(self) -> list[str]:
        return sorted(self._by)


class Points:
    """Sorted positions per chromosome, asked whether a window contains one."""

    def __init__(self) -> None:
        self._by: dict[str, set[int]] = collections.defaultdict(set)
        self._sorted: dict[str, list[int]] = {}

    def add(self, chrom: str, pos: int) -> None:
        self._by[chrom].add(pos)

    def freeze(self) -> Points:
        self._sorted = {c: sorted(v) for c, v in self._by.items()}
        return self

    def inside(self, chrom: str, lo: int, hi: int) -> list[int]:
        p = self._sorted.get(chrom)
        if not p:
            return []
        return p[bisect.bisect_left(p, lo) : bisect.bisect_right(p, hi)]

    @property
    def n(self) -> int:
        return sum(len(v) for v in self._by.values())

    @property
    def by_chromosome(self) -> dict[str, int]:
        return {c: len(v) for c, v in sorted(self._by.items())}


# --------------------------------------------------------------------------- reading


def elements(root: Path = ROOT) -> dict[str, dict[str, Any]]:
    """The CRISPRi element frame: `response_map2.candidates`'s own attachment, as the frame used it."""
    cands, _counted = rm2.candidates(list(rm2.CHROMS), RESULTS_DIR, ce.Readers(), ce.mapping())
    out: dict[str, dict[str, Any]] = {}
    for c in cands:
        if c.element is None:
            continue
        rec = out.setdefault(
            c.element["id"],
            {
                "chrom": c.pair.chrom,
                "start": c.element["start"],
                "end": c.element["end"],
                "genes": set(),
                "statuses": set(),
            },
        )
        rec["genes"].add(c.pair.gene)
        rec["statuses"].add(c.status)
    return out


def fine_mapped() -> tuple[Points, dict[tuple[str, int], set[str]], dict[str, int]]:
    """Every fine-mapped variant at PIP >= 0.5 anywhere on disk, with the genes it is fine-mapped for."""
    pts = Points()
    genes: dict[tuple[str, int], set[str]] = {}
    rows_by_file: dict[str, int] = {}

    def keep(chrom: str, pos: int, gene: str | None, pip: float | None, name: str) -> None:
        if pip is None or pip < PIP:
            return
        pts.add(chrom, int(pos))
        genes.setdefault((chrom, int(pos)), set()).add(gene or "")
        rows_by_file[name] = rows_by_file.get(name, 0) + 1

    def walk(o: Any, name: str, unit: str | None = None) -> None:
        if isinstance(o, dict):
            u = o["unit"] if isinstance(o.get("unit"), str) else unit
            if o.get("pip") is not None and o.get("pos") is not None and u:
                keep(u.split(":")[0], o["pos"], o.get("gene"), o.get("pip"), name)
            v = o.get("variant")
            if isinstance(v, str) and ":" in v and o.get("pip") is not None:
                c, p = v.split(":")[:2]
                if p.isdigit():
                    keep(c, int(p), o.get("gtex_gene") or o.get("gene"), o.get("pip"), name)
            for x in o.values():
                walk(x, name, u)
        elif isinstance(o, list):
            for x in o:
                walk(x, name, unit)

    for f in PIP_FILES:
        p = EXECUTOR / f
        if not p.exists():
            continue
        with gzip.open(p, "rt") as fh:
            walk(json.load(fh), f)
    return pts.freeze(), genes, rows_by_file


def queried_upper_bound() -> tuple[Spans, Points]:
    """The widest region a PIP could ever have been recorded over, from what is on disk.

    DAP-G was read over the human panel's storage units and over MPRAVarDB's tested variants, never
    genome-wide: `genomeos/attribution/human_panel.py` reads `gtex_dapg` for a unit's own interval,
    and the panel's catalogues are the only unit tables on disk. Outside the union below no posterior
    was ever looked up, so no element there can be assessed whatever is really present at it.
    """
    units = Spans()
    for p in sorted(PANEL.glob("chr*/storage_catalogue.json.gz")):
        chrom = p.parent.name
        with gzip.open(p, "rt") as fh:
            for u in json.load(fh):
                units.add(chrom, u["start"], u["end"])
    tested = Points()
    p = EXECUTOR / "two_instrument_rows.json.gz"
    if p.exists():
        with gzip.open(p, "rt") as fh:
            for r in json.load(fh):
                v = r.get("variant", "")
                if ":" in v:
                    c, q = v.split(":")[:2]
                    if q.isdigit():
                        tested.add(c, int(q))
    return units.freeze(), tested.freeze()


def queried_lower_bound() -> Spans:
    """The regions the executor's own records show it looked at, a floor on the queried frame.

    Every unit the executor names in its eQTL records is a unit it read. Units it read and found
    nothing at do not all appear, so this is a floor and not the frame.
    """
    s = Spans()
    for f in ("replication_assembled.json.gz", "assembled.json.gz"):
        p = EXECUTOR / f
        if not p.exists():
            continue
        with gzip.open(p, "rt") as fh:
            d = json.load(fh)

        def walk(o: Any) -> None:
            if isinstance(o, dict):
                u = o.get("unit")
                if isinstance(u, str) and ":" in u and "-" in u:
                    chrom, rest = u.split(":")
                    a, b = rest.split("-")
                    s.add(chrom, int(a), int(b))
                for v in o.values():
                    walk(v)
            elif isinstance(o, list):
                for v in o:
                    walk(v)

        walk(d)
    return s.freeze()


# --------------------------------------------------------------------------- the count


def count(root: Path = ROOT) -> dict[str, Any]:
    els = elements(root)
    pts, var_genes, rows_by_file = fine_mapped()
    units, tested = queried_upper_bound()
    floor_spans = queried_lower_bound()

    by_chrom = collections.Counter(v["chrom"] for v in els.values())
    assessable_upper = assessable_units = assessable_mpra = assessable_floor = 0
    unassessable_by_chrom: collections.Counter[str] = collections.Counter()
    carriers: list[dict[str, Any]] = []
    for eid, v in sorted(els.items()):
        c, lo, hi = v["chrom"], v["start"] - MARGIN, v["end"] + MARGIN
        u = units.touches(c, lo, hi)
        m = bool(tested.inside(c, lo, hi))
        assessable_units += u
        assessable_mpra += m
        assessable_floor += floor_spans.touches(c, lo, hi)
        if u or m:
            assessable_upper += 1
        else:
            unassessable_by_chrom[c] += 1
        hits = pts.inside(c, lo, hi)
        if hits:
            carriers.append(
                {
                    "element": eid,
                    "chrom": c,
                    "linked_genes": sorted(v["genes"]),
                    "fine_mapped": [
                        {"pos": p, "genes": sorted(g for g in var_genes[(c, p)] if g)} for p in hits
                    ],
                }
            )

    symbols = c3.symbol_map({v["chrom"] for v in els.values()})
    matched, unresolved = 0, []
    for rec in carriers:
        linked = set(rec["linked_genes"])
        names: set[str] = set()
        for h in rec["fine_mapped"]:
            for g in h["genes"]:
                if g.startswith("ENSG"):
                    resolved = symbols.get(g.split(".")[0])
                    if not resolved or resolved == {g}:
                        unresolved.append(g)
                        continue
                    names |= resolved
                else:
                    names.add(g)
        rec["fine_mapped_gene_symbols"] = sorted(names)
        rec["gene_is_the_linked_gene"] = bool(names & linked)
        matched += bool(names & linked)

    return {
        "elements": els,
        "n_elements": len(els),
        "elements_by_chromosome": dict(sorted(by_chrom.items())),
        "fine_mapped_points": pts,
        "rows_by_file": rows_by_file,
        "n_fine_mapped": pts.n,
        "fine_mapped_by_chromosome": pts.by_chromosome,
        "units": units,
        "tested": tested,
        "floor_spans": floor_spans,
        "assessable_upper": assessable_upper,
        "assessable_units": assessable_units,
        "assessable_mpra": assessable_mpra,
        "assessable_floor": assessable_floor,
        "unassessable_by_chromosome": dict(sorted(unassessable_by_chrom.items())),
        "carriers": carriers,
        "matched": matched,
        "unresolved_gene_ids": sorted(set(unresolved)),
    }


def payload(k: dict[str, Any]) -> dict[str, Any]:
    n = k["n_elements"]
    return {
        "result": NAME,
        "date": "2026-10-02",
        "question": (
            "how many of the 1,505 compiled CRISPRi elements are ASSESSABLE for a fine-mapped eQTL "
            "at all, counted before and apart from how many carry one"
        ),
        "why": (
            "data/results/eqtl_crispri_frame.json reports 1,262 of 1,505 elements carrying a "
            "retained cis-eQTL and 96 addable loci over a floor of 20. Those hits are GTEx v8 "
            "single-tissue SIGNIFICANT variant-gene pairs and carry no posterior: "
            f"`genomeos.attribution.eqtl.HIT_COLUMNS` is {list(eqtl.HIT_COLUMNS)}. Significant "
            "variants travel in LD blocks, so an element-sized window catches one almost anywhere "
            "and 84% of elements carrying one is close to what density alone gives"
        ),
        "the_distinction_this_applies": THE_DISTINCTION,
        "imported_not_chosen": {
            "pip_threshold": PIP,
            "pip_threshold_source": (
                "genomeos.attribution.executor.DAPG_PIP, whose registered wording in "
                f'data/results/executor_replication.json pre_registration.population is: "{E2R_WORDING}"'
            ),
            "margin_bp": MARGIN,
            "margin_source": (
                "scripts/eqtl_targets.MARGIN, imported through "
                "scripts/response_map_increment3_count.EQTL_MARGIN, not chosen here"
            ),
            "label_for_an_overlap_that_is_not_fine_mapped": LD_TAGGED_LABEL,
            "label_source": (
                "this project's own existing distinction, not a new category: "
                "`E2R_eqtl_replication` is the fine-mapped arm and "
                f'`E3R_eqtl_linked_replication` is "{E3R_WORDING}"'
            ),
            "element_frame": (
                "`response_map2.candidates`'s own attachment at "
                f"`measured.RECIPROCAL_OVERLAP = {ms.RECIPROCAL_OVERLAP}`, the same call the "
                "re-distillation's frame was built from"
            ),
            "locus_convention": (
                "`genomeos.attribution.cell2.group`'s, carrying its wording that it is an "
                "operational grouping and NOT established biological independence"
            ),
        },
        "the_fine_mapped_record_on_disk": {
            "where_the_posterior_lives": DAPG_TRACK,
            "source_note": DAPG_SOURCE,
            "the_track_is_not_on_disk": (
                "no local copy of gtexDapg.bb exists under data/. Nothing was fetched: this writer "
                "makes no network request"
            ),
            "files_that_carry_a_pip": k["rows_by_file"],
            "distinct_variants_at_or_above_the_threshold": k["n_fine_mapped"],
            "by_chromosome": k["fine_mapped_by_chromosome"],
            "what_this_set_is": (
                "not a genome-wide fine-mapping record. It is whatever a previous run retained from "
                "the remote track, and that run read the track only over the human panel's storage "
                "units and over MPRAVarDB's tested variants, one interval at a time "
                "(`genomeos/attribution/human_panel.py` reads `gtex_dapg` for a unit's own interval). "
                "It was assembled for a different population and its coverage of the CRISPRi element "
                "set is the open question this writer answers"
            ),
        },
        "ASSESSABLE_FIRST": {
            "quantity": (
                "compiled CRISPRi elements, widened by the imported margin, that lie inside a region "
                "where a fine-mapped posterior could have been recorded at all"
            ),
            "elements_of_the_frame": n,
            "upper_bound": k["assessable_upper"],
            "upper_bound_rule": (
                "the element window touches a human-panel storage unit OR contains an "
                "MPRAVarDB-tested variant: the union of everywhere DAP-G was ever read. It is an "
                "UPPER bound because only the subset of units the executor shortlisted was actually "
                "queried, not every catalogued unit"
            ),
            "upper_bound_from_storage_units_alone": k["assessable_units"],
            "upper_bound_from_mpravardb_variants_alone": k["assessable_mpra"],
            "lower_bound": k["assessable_floor"],
            "lower_bound_rule": (
                "the element window touches a unit the executor's own eQTL records name, so a unit "
                "it demonstrably read. A floor, because a unit read and found empty need not appear"
            ),
            "unassessed": n - k["assessable_upper"],
            "unassessed_meaning": (
                "no fine-mapped posterior was ever looked up over these elements, so they are "
                "UNASSESSED for fine-mapping and not elements without a fine-mapped variant. This is "
                "the same distinction the re-distillation's registration imposed on itself, applied "
                "one level up"
            ),
            "unassessed_by_chromosome": k["unassessable_by_chromosome"],
            "elements_by_chromosome": k["elements_by_chromosome"],
            "a_bound_is_not_a_denominator": (
                "the assessable count is known only between a floor and a ceiling, so no ratio of "
                "carriers to assessed elements is reportable. That is the finding, not a limitation "
                "of the reporting"
            ),
        },
        "CARRYING_SECOND": {
            "quantity": (
                "elements of the frame with at least one fine-mapped variant at PIP >= "
                f"{PIP} within {MARGIN} bases"
            ),
            "carrying_any_fine_mapped_variant": len(k["carriers"]),
            "carrying_one_whose_gene_is_the_element_linked_gene": k["matched"],
            "gene_rule": (
                "condition (a) of the authorisation: the eQTL's gene must be the element's linked "
                "gene. An eQTL for some other gene inside the window is not evidence about this chain"
            ),
            "gene_ids_that_could_not_be_resolved_to_a_symbol": k["unresolved_gene_ids"],
            "unresolved_note": (
                "`scripts/eqtl_targets.py symbol_map` leaves a gene whose GENCODE record carries no "
                "symbol as its own Ensembl id, so a match count is a floor. Each unresolved id is "
                "named above rather than counted as a non-match"
            ),
            "carriers": k["carriers"],
        },
        "verdict": {
            "coverage": "poor",
            "what_follows": (
                "a report and nothing else. The authorisation made the two further conditions - the "
                "gene match and the matched-window baseline - conditional on coverage supporting a "
                "count. It does not, so no baseline was registered, no baseline was constructed, no "
                "excess was computed and nothing was built"
            ),
            "why_poor": (
                f"at most {k['assessable_upper']} of {n} elements could have been assessed and at "
                f"least {n - k['assessable_upper']} could not. Of the frame, "
                f"{len(k['carriers'])} elements carry a fine-mapped variant in range and "
                f"{k['matched']} carry one for the element's own linked gene"
            ),
            "the_number_that_cannot_be_recomputed": (
                "the floor quantity of the re-distillation is addable loci carrying a hit on a shown "
                "element, reached at 96. Under the gene match no element of the frame carries a "
                "fine-mapped eQTL for its linked gene, so that quantity is 0 under fine-mapping and "
                "no locus count derived from it can be anything else. 96 and 0 are not two readings "
                "of one measurement: 96 counts LD-tagged associations over a frame that was fully "
                "assessed, and 0 counts localised associations over a frame that mostly was not"
            ),
            "what_the_1262_may_be_called": (
                f'{LD_TAGGED_LABEL}. The overlap supports "a variant in LD with an association lies '
                'here" and not "the association is localised here". Nothing here says the 1,262 '
                "is wrong as a count of what it counted"
            ),
            "what_would_answer_it": (
                "the fine-mapped track itself, read against the CRISPRi element set the way the "
                f"re-distillation read the significant set: {DAPG_TRACK}. That is a fetch and a "
                "ruling, not a step this writer takes"
            ),
        },
        "alphagenome_requests": 0,
        "money": "none: every input was already on disk and no network request was made",
    }


def _inputs() -> list[dict[str, Any]]:
    """Everything this writer opens under data/, including what it reaches through a pointer field.

    Declared repository-relative throughout: an absolute path cannot be satisfied from a second
    checkout, and two of them once made a rebuild report 979 of 981 and refuse a verdict.
    """
    from genomeos.genome import default_gencode, reader

    seen: set[str] = set()
    out: list[dict[str, Any]] = []

    def put(p: Path, partition: str | None = None) -> None:
        q = p if not p.is_absolute() else p.relative_to(ROOT)
        key = q.as_posix()
        if key in seen or not q.exists():
            return
        seen.add(key)
        out.append(mf.input_entry(q, partition=partition))

    put(Path(ce.TRACK_METADATA))
    for f in ms.CRISPRI_FILES:
        put(ms.CRISPRI_KNOWLEDGE / f, ms.CRISPRI_SPLIT_OF[f])
    paths: list[Path] = [
        p
        for cell in ce.READER_TERMS
        for chrom in rm2.CHROMS
        if (p := RESULTS_DIR / reader.peaks_path(cell, chrom).name).exists()
    ]
    for g in (
        "budget_chr*.json",
        "budget_axes_chr*.json",
        "variation_chr*.json",
        "duplication_chr*.json",
        "domains_chr*.json",
        "unknown_chr*.json",
        "reader_*_chr*.json",
        "ccres_chr*.bed.gz",
        "enhancer_targets_chr*.json",
        "enhancer_targets_all_chr*.json",
        "constrained_targets_chr*.json",
        "vista_chr*.json",
    ):
        for chrom in rm2.CHROMS:
            paths += sorted(RESULTS_DIR.glob(g.replace("chr*", chrom)))
    for group in c3.OPENED_BESIDE_THE_RESULTS.values():
        paths += [Path(x) for x in group]
    paths += [EXECUTOR / f for f in PIP_FILES]
    for r in ("eqtl_crispri_frame_registration.json", "eqtl_crispri_frame.json", "executor_replication.json"):
        paths.append(Path("data/results") / r)
    paths += sorted(PANEL.glob("chr*/storage_catalogue.json.gz"))
    for chrom in rm2.CHROMS:
        g = default_gencode({chrom})
        if g is not None:
            paths.append(Path(g))
    for q in paths:
        put(q)
    return out


def manifest() -> dict[str, Any]:
    return {
        "sources": [
            {"accession": DAPG_SOURCE, "version": DAPG_TRACK},
            {
                "accession": "EngreitzLab/CRISPR_comparison (Gschwind et al. 2025)",
                "version": "the two local benchmark tables under data/knowledge/crispri",
            },
            {
                "accession": "GTEx v8 single-tissue cis-eQTLs (significant variant-gene pairs), hg38",
                "version": "named only to say what the 1,262 counted; not read by this writer",
            },
        ],
        "inputs": _inputs(),
        "inputs_opened_beside_the_declared_results": {
            "why_they_are_listed": (
                "an input is what the writer opened, not what it said it read. The audit hook in "
                "`genomeos.manifest` records every open under `data/` and `save_result` refuses a "
                "result whose traced reads exceed its declared inputs"
            ),
            "groups": {k: len(v) for k, v in c3.OPENED_BESIDE_THE_RESULTS.items()},
            "the_pointer_case": (
                "`response_map2.candidates` opens the 24 per-chromosome element tables through a "
                "result file's `elements_where` POINTER field - about 1.4 GB that no reading of the "
                "source can name, and how 17 published results read a table none declared"
            ),
        },
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "dapg_pip": PIP,
            "eqtl_margin_bp": MARGIN,
            "reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
            "independent_locus_span": cell2.INDEPENDENT_LOCUS_SPAN,
            "chromosomes": list(rm2.CHROMS),
        },
        "exclusions": [
            "no element of the frame is excluded: it is every element a valid CRISPRi pair attaches to",
            "no fine-mapped variant on disk is excluded other than by the imported PIP threshold",
            "no matched-window baseline was constructed, because coverage did not support a count",
        ],
        "partitions": {
            "frame": "the CRISPRi element set, the same frame the re-distillation indexed against",
            "fine_mapped_record": (
                "whatever a previous run retained from the DAP-G track over the human panel's "
                "storage units and MPRAVarDB's tested variants - a different population"
            ),
        },
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--save", action="store_true", help="write data/results/finemap_coverage.json")
    args = ap.parse_args()
    k = count()  # `genomeos.manifest` installs its one audit hook at import, so every open is traced
    p = payload(k)
    print(json.dumps({kk: vv for kk, vv in p["ASSESSABLE_FIRST"].items() if isinstance(vv, int)}, indent=1))
    print(json.dumps({kk: vv for kk, vv in p["CARRYING_SECOND"].items() if isinstance(vv, int)}, indent=1))
    if args.save:
        out = save_result(NAME, p, manifest=manifest())
        print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
