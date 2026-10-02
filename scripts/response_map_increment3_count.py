# SPDX-License-Identifier: AGPL-3.0-or-later
"""Count, before anything is built, what a third increment of the response map could read.

    uv run --frozen python scripts/response_map_increment3_count.py
    uv run --frozen python scripts/response_map_increment3_count.py --chroms chr21 \
        --result response_map_increment3_count_chr21

Increment 2 built a chain over 110 loci and stated its own limit, which is carried into the result
below word for word: the evidenced chain runs over 110 loci, and that is "23.26% of one of the
measured layer's four assays"; lentiMPRA, VISTA and saturation mutagenesis contributed nothing to it,
because "they read a sequence, or the bases inside it, not the native locus after a perturbation of
it". So increment 2 is one assay deep. The question here is how much of a second, third or fourth kind
of evidence is actually on disk at those loci, and what each of those observations is.

**Counting only.** No increment is built, nothing is added to the map, nothing is scored or fitted, no
cut-off is introduced, 0 model requests, no downloads: every input was already on disk.

**The counting method is imported, not rewritten.** `genomeos.response_map2.candidates` is the one
implementation of the three-input rule and of the order `cell2.group` is given its keys, and
`genomeos.response_map2.map1_keys` supplies increment 1's keys for the same grouping; both come from
lane-rmap2's count at `619148b` by way of the increment-2 builder that reuses it. The 110 loci are
re-derived by that function here rather than read out of the increment-2 result, so this count stands
or falls on reproducing the number it starts from.

**What an observation is, per kind of evidence**, taken from `measured.OUTCOME_KIND` where it was
registered and never restated here. The distinction increment 2 drew is the whole point of this count:
a reporter assay reads a *copy of a sequence* in a reporter, and that is not an observation of the
native locus after a perturbation of it. Such an assay is therefore counted at **zero native-locus
observations**, and the sequence observations it does carry are counted separately and under their own
name. The rule is not bent to admit them.

**GTEx eQTL** is the one further kind that reads the native locus: a significant cis-eQTL is a variant
inside the element whose alleles go with a gene's expression across donors. It is natural variation and
an association, not an experimental perturbation, and it is counted as its own kind. Its retained set
was distilled against the 9,286 intervals of three sampled element sets, so **a locus outside that
frame could not have been seen at all**: the assessable denominator is reported beside the count and
the count is a floor over the frame, never a genome-wide eQTL query.
"""

from __future__ import annotations

import argparse
import bisect
import builtins
import collections
import glob as globmod
import importlib.util
import json
import sys
import time
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos import response_map2 as rm2  # noqa: E402
from genomeos.attribution import cell2, eqtl  # noqa: E402
from genomeos.attribution import compile as cp  # noqa: E402
from genomeos.attribution import context_evidence as ce  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.attribution import targets as tg  # noqa: E402
from genomeos.genome import reader  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

RESULT = "response_map_increment3_count"
CHROMS = rm2.CHROMS
ROOT = Path(__file__).resolve().parents[1]
ENTRY = "scripts/response_map_increment3_count.py"
OWN_CODE = (
    "scripts/response_map_increment3_count.py",
    "tests/test_response_map_increment3.py",
)

#: The counting functions imported rather than written a second time, named in the result so a reader
#: can check that no parallel method was introduced here.
IMPORTED_COUNTING = (
    "genomeos.response_map2.candidates - the three-input status of every CRISPRi pair, in the order "
    "scripts/response_map_coverage.py grouped in (lane-rmap2, 619148b)",
    "genomeos.response_map2.map1_keys - increment 1's locus keys, grouped alongside",
    "genomeos.response_map2.STATUS_ORDER and STATUS_ADDABLE - the best-status-per-locus rule",
    "genomeos.attribution.cell2.group - the independent-locus grouping, imported where registered",
    "genomeos.attribution.compile._attributed and ._measured_rows - the measured layer's own "
    "element and assay enumeration",
    "genomeos.attribution.eqtl.load_hits - the retained GTEx hits, read from local knowledge",
    "scripts/eqtl_targets.py symbol_map - the one Ensembl-id-to-symbol resolution, loaded from its "
    "own file rather than written again here",
)


def symbol_map(chroms: set[str]) -> dict[str, set[str]]:
    """Ensembl gene id to symbol, by `scripts/eqtl_targets.py`'s own resolution, loaded from it.

    Imported from the file rather than copied, so one implementation names GTEx's genes everywhere.
    A gene whose GENCODE record carries no symbol keeps its Ensembl id as its symbol there, so such a
    gene cannot match a benchmark symbol and the match count below is a floor.
    """
    spec = importlib.util.spec_from_file_location("_eqtl_targets", ROOT / "scripts/eqtl_targets.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    out: dict[str, set[str]] = {}
    for gid, sym in mod.symbol_map(chroms).items():
        out.setdefault(gid, set()).add(sym)
    return out


#: Increment 2's stated limit, carried word for word. This is the starting point of the count.
INCREMENT_2_LIMIT = (
    "the evidenced chain runs over 110 loci, and that is \"23.26% of one of the measured layer's four "
    'assays". lentiMPRA, VISTA and saturation mutagenesis contributed NOTHING to it, because "they '
    'read a sequence, or the bases inside it, not the native locus after a perturbation of it"'
)

#: The three element sets the GTEx distillation indexed, from `scripts/eqtl_targets.py`, whose MARGIN
#: is imported below rather than chosen here. Restated as globs so this script does not import a
#: module that runs an argument parser.
EQTL_SETS = {
    "uniform": "data/results/enhancer_targets_chr*.json",
    "constrained": "data/results/constrained_targets_chr*.json",
    "vista": "data/results/vista_chr*.json",
}
EQTL_MARGIN = 500  # scripts/eqtl_targets.MARGIN, the margin the retained hits were distilled under
EQTL_RESULT = "data/results/eqtl_targets.json"

#: Files this count opens that no glob of `data/results` reaches, each with why it is opened. They are
#: declared here because an input is what the writer opened, not what it said it read: the audit hook
#: in `genomeos.manifest` records every open under `data/` and `save_result` refuses a result whose
#: reads exceed its declared inputs. Three of these groups would have gone undeclared without it.
OPENED_BESIDE_THE_RESULTS = {
    # `attribution.targets` reads a run's elements out of the table its result's `elements_where` field
    # points at, so no reading of the source names these: the path comes out of a result file at run
    # time. 1.4 GB over 26 chromosomes, and the reason this project traces reads instead of trusting
    # declarations.
    "the per-chromosome element tables reached through a result file's `elements_where` pointer": [
        f"data/knowledge/alphagenome/all_elements/{c}.json" for c in CHROMS
    ],
    # the measured layer loads all four of its assays per chromosome, not only the perturbation one, so
    # the three reporter and base-level sources are opened even though they contribute no native-locus
    # observation to this count
    "the measured layer's other three assays, opened by `compile._measured_rows`": [
        *[f"data/knowledge/mpra/{a}.bed.gz" for a in ("ENCFF475FKV", "ENCFF769REH", "ENCFF802FUV")],
        "data/knowledge/vista/locus.tsv.gz",
        "data/knowledge/satmut/elements.tsv.gz",
    ],
    # increment 1's map is read whole by `response_map2.map1_keys`, which is how this count knows which
    # locus the map already covered; everything that view reads is opened with it
    "increment 1's own view, read by `response_map2.map1_keys`": [
        "data/results/attribution_correctness_v3.json",
        "data/results/crispri_direction.json",
        "data/results/discovery_review.json",
        "data/results/loci_benchmark.json",
        "data/knowledge/ReactomePathways.txt",
        "data/knowledge/compiled/noncoding_chr11.bio",
        "data/knowledge/hic_contact/K562_4DNFITUOMFUQ_5000.json",
        "data/organisms/human/erythrocyte.bio",
        "data/cache/rates/schofield2018_TableS2_halflives.xlsx",
    ],
}
#: GENCODE is read per chromosome, and only for the chromosomes carrying a retained eQTL hit, so the
#: set is not fixed: it is collected from what was actually opened and declared by name.
GENCODE = "data/reference/gencode_v50_{chrom}.gff3.gz"

#: A gene token for a kind of evidence that names no measured gene. `cell2.group` unions on a shared
#: measured gene or on 1 Mb proximity; a reporter observation names no gene, so a token unique to the
#: element suppresses the gene arm and leaves the proximity arm alone.
NO_MEASURED_GENE = (
    "a reporter or base-level observation names no measured gene, so its locus key carries a token "
    "unique to its element and `cell2.group`'s gene arm cannot fire for it. Only the 1 Mb proximity "
    "arm joins such a key to anything, so the further-locus count below is an UPPER BOUND on loci: "
    "with gene identity the groups could only be fewer, never more"
)

#: Why each kind is or is not an observation of the native locus. Stated per kind, beside the
#: registered wording of what the assay measures, and never softened into "an observation at the
#: locus".
NATIVE_LOCUS_RULE = (
    "an observation OF a native locus is a measurement taken on the chromosome in place: the locus "
    "perturbed and the consequence read out (CRISPRi), the locus's own state read out (the reader's "
    "DNase peak), or the locus's own natural variation read out against expression (a GTEx cis-eQTL). "
    "A reporter assay reads a copy of the sequence outside the locus and is not such an observation, "
    "however well measured it is"
)


class _CacheAudit:
    """Every open of a file under the per-element response cache, recorded by patching `open`.

    Whether a writer reads the per-element cache cannot be seen by reading code:
    `targets.run_elements` reads a path out of a result file's field. So it is recorded here instead
    of asserted, and the count goes into the result whatever it turns out to be.
    """

    def __init__(self, root: Path = tg.ELEMENT_CACHE) -> None:
        #: the path as the cache declares it, reported as such. `root` is resolved, for matching only:
        #: in a rebuild worktree the store is a read-only link, so a resolved path there names the
        #: checkout it points at and the report would differ between machines without differing in
        #: substance.
        self.declared = Path(root).as_posix()
        self.root = Path(root).resolve()
        self.opens: list[str] = []
        self._real = builtins.open

    def __enter__(self) -> _CacheAudit:
        audit = self

        def patched(file, *a, **kw):  # type: ignore[no-untyped-def]
            try:
                p = Path(file).resolve()
                if p == audit.root or audit.root in p.parents:
                    audit.opens.append(p.as_posix())
            except (TypeError, ValueError, OSError):
                pass
            return audit._real(file, *a, **kw)

        builtins.open = patched  # type: ignore[assignment]
        return self

    def __exit__(self, *exc: object) -> None:
        builtins.open = self._real  # type: ignore[assignment]

    def _named(self) -> set[str]:
        """Each open spelled under the declared root, so the list does not carry this machine's paths."""
        out = set()
        for x in self.opens:
            p = Path(x)
            under = p.is_relative_to(self.root)
            out.add((Path(self.declared) / p.relative_to(self.root)).as_posix() if under else p.name)
        return out

    def block(self) -> dict[str, Any]:
        return {
            "per_element_response_cache_root": self.declared,
            "opens": len(self.opens),
            "files": sorted(self._named()),
            "how_this_is_known": (
                "recorded, not inferred: `builtins.open` was patched for the whole run and every open "
                "of a path under the cache root was appended. No grep of the source could see it, "
                "because `targets.run_elements` reads its path out of a result file's field"
            ),
        }


def eqtl_index(root: Path = ROOT) -> tuple[dict[str, list[tuple[int, int]]], int]:
    """The intervals the GTEx distillation indexed, rebuilt with its own margin.

    A retained hit exists only inside one of these, so an element outside them could not have been
    seen whatever eQTLs are really there. Returned so that absence can be reported as unassessed
    rather than as no eQTL.
    """
    out: dict[str, set[tuple[int, int]]] = {}
    for pattern in EQTL_SETS.values():
        for f in sorted(globmod.glob((root / pattern).as_posix())):
            with open(f) as fh:
                d = json.load(fh)
            rows = d.get("rows") if isinstance(d.get("rows"), list) else d.get("elements") or []
            for e in rows:
                if "start" not in e:
                    continue
                out.setdefault(d["chrom"], set()).add(
                    (max(0, e["start"] - EQTL_MARGIN), e["end"] + EQTL_MARGIN)
                )
    n = sum(len(v) for v in out.values())
    return {c: sorted(v) for c, v in out.items()}, n


def eqtl_hits_by_chrom(knowledge: Path) -> dict[str, list[dict[str, Any]]]:
    """Every retained GTEx hit, by chromosome and sorted by position, read a line at a time."""
    out: dict[str, list[dict[str, Any]]] = {}
    for p in sorted(knowledge.glob("hits_*.tsv")):
        with p.open() as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                out.setdefault(f[1], []).append(
                    {
                        "tissue": f[0],
                        "pos": int(f[2]),
                        "gene_id": f[5],
                        "slope": float(f[6]),
                        "pval": float(f[7]),
                    }
                )
    for c in out:
        out[c].sort(key=lambda r: r["pos"])
    return out


def covered_by(index: list[tuple[int, int]], start: int, end: int) -> bool:
    """Whether a half-open interval overlaps any indexed interval at all."""
    if not index:
        return False
    i = bisect.bisect_right([s for s, _ in index], end) - 1
    while i >= 0:
        s, e = index[i]
        if e > start:
            return True
        if start - s > 10_000_000:
            break
        i -= 1
    return False


def hits_in(rows: list[dict[str, Any]], start: int, end: int) -> list[dict[str, Any]]:
    """The retained hits whose variant position falls in the element, widened by the distil margin."""
    pos = [r["pos"] for r in rows]
    lo = bisect.bisect_right(pos, start - EQTL_MARGIN)
    hi = bisect.bisect_right(pos, end + EQTL_MARGIN)
    return rows[lo:hi]


#: Why the outcome breakdown is taken at every denominator, not just the count. Stated because the
#: project has twice reported a count as if it were the thing the count was about.
OUTCOME_BREAKDOWN_RULE = (
    "a count of loci carrying a kind of evidence is not a measurement of what that evidence says. "
    "198 loci where the model says nothing were once quoted as 198 places it is wrong when 66 held a "
    "measured decrease; 48 measured links cleared a floor while only 21 were answerable. So every "
    "count below carries what its observations actually say, at its own denominator"
)

#: Why a silent reporter beside a measured CRISPRi decrease is not a contradiction.
SILENT_IS_NOT_A_CONTRADICTION = (
    "a reporter calling the sequence silent and a CRISPRi screen measuring a decrease at the native "
    "locus are not in conflict: the reporter read a 200 bp copy outside the locus and the screen "
    "perturbed the locus in place. Neither checks the other, and no agreement or disagreement is "
    "counted here"
)


def reporter_outcomes(
    eids: set[str],
    element_row: dict[str, dict[str, Any]],
    element_cells: dict[str, set[str]],
    assay: str,
) -> dict[str, Any]:
    """What the reporter observations on these elements actually say, by the assay's own label rule.

    Reported per element and then split by whether the reporter cell is the cell the CRISPRi
    observation of the same chain was taken in, because a label in another cell is another context.
    """
    by_element: collections.Counter[str] = collections.Counter()
    labels: collections.Counter[str] = collections.Counter()
    same_cell: collections.Counter[str] = collections.Counter()
    cells_seen: collections.Counter[str] = collections.Counter()
    no_overlap = 0
    for eid in sorted(eids):
        blk = (element_row[eid].get("measured") or {}).get(assay) or {}
        per_cell = blk.get("label_by_cell") or {}
        for cell, label in per_cell.items():
            labels[label] += 1
            cells_seen[cell] += 1
        if blk.get("cells_active"):
            by_element["at_least_one_cell_active"] += 1
        elif blk.get("cells_conflicting"):
            by_element["no_cell_active_and_one_conflicting"] += 1
        elif per_cell:
            by_element["every_cell_silent"] += 1
        else:
            by_element["no_cell_label_in_the_record"] += 1
        shared = set(per_cell) & element_cells.get(eid, set())
        if not shared:
            no_overlap += 1
        for cell in sorted(shared):
            same_cell[per_cell[cell]] += 1
    return {
        "elements": len(eids),
        "by_element": dict(by_element),
        "cell_labels": dict(labels),
        "reporter_cells": dict(cells_seen),
        "in_the_same_cell_as_the_chains_crispri_observation": dict(same_cell),
        "elements_whose_reporter_cells_include_none_of_the_chains_cells": no_overlap,
        "label_rule": ms.REPORTER_LABEL_RULE,
        "silent_is_not_a_contradiction": SILENT_IS_NOT_A_CONTRADICTION,
    }


def eqtl_outcomes(
    eqtl_shown: dict[str, list[dict[str, Any]]],
    element_genes: dict[str, set[str]],
    symbols: dict[str, set[str]],
) -> dict[str, Any]:
    """What the retained eQTLs on these elements actually say, element by element.

    Whether the gene the CRISPRi screen measured is among the retained eGenes is counted as a
    DESCRIPTIVE coincidence of two observations in different contexts, never as agreement.
    """
    rows = []
    both = 0
    for eid, hits in sorted(eqtl_shown.items()):
        egenes = sorted({h["gene_id"] for h in hits})
        tissues = sorted({h["tissue"] for h in hits})
        names = set()
        for g in egenes:
            names |= symbols.get(g, set())
        measured = sorted(element_genes.get(eid, set()))
        shared = sorted(set(measured) & names)
        both += bool(shared)
        rows.append(
            {
                "element": eid,
                "retained_eqtls": len(hits),
                "egenes": len(egenes),
                "tissues": len(tissues),
                "slopes_negative": sum(1 for h in hits if h["slope"] < 0),
                "slopes_positive": sum(1 for h in hits if h["slope"] > 0),
                "genes_the_crispri_screen_measured_here": measured,
                "of_those_also_a_retained_egene": shared,
                "egene_symbols_resolved": sorted(names)[:40],
            }
        )
    return {
        "elements": len(eqtl_shown),
        "per_element": rows,
        "elements_where_a_measured_gene_is_also_a_retained_egene": both,
        "that_is_not_agreement": (
            "the eQTL is an association over GTEx donor tissues and the CRISPRi effect is a "
            "perturbation in a cell line. A gene appearing in both is two observations in different "
            "contexts naming the same gene, and is counted descriptively; it is not agreement, "
            "confirmation or validation, and the slope sign is not compared with the CRISPRi direction"
        ),
        "egene_symbol_resolution": (
            "GTEx names genes by Ensembl id and the CRISPRi benchmark by symbol. Symbols are read "
            "from the committed eqtl_targets.json, which resolved them once; an Ensembl id that "
            "result did not resolve stays unresolved and cannot match, so this count is a floor"
        ),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", default=",".join(CHROMS))
    ap.add_argument("--result", default=RESULT)
    args = ap.parse_args()
    chroms = [c for c in args.chroms.split(",") if c]
    name = args.result
    if name == RESULT and tuple(chroms) != tuple(CHROMS):
        raise SystemExit(
            f"{RESULT} is the genome-wide count; a run over {len(chroms)} chromosomes needs --result"
        )
    started = time.time()

    with _CacheAudit() as audit:
        readers = ce.Readers()
        table = ce.mapping()
        cands, counted = rm2.candidates(chroms, RESULTS_DIR, readers, table)

        # --- the 110, re-derived by the imported function, not read out of increment 2's result ----
        covered = rm2.map1_keys(ROOT)
        groups = cell2.group([c.key for c in cands] + covered)
        mine = groups[: len(cands)]
        best: dict[int, str] = {}
        for g, c in zip(mine, cands, strict=True):
            if g not in best or rm2.STATUS_ORDER[c.status] < rm2.STATUS_ORDER[best[g]]:
                best[g] = c.status
        addable = {g for g, s in best.items() if s == rm2.STATUS_ADDABLE}

        # the elements the built chains actually stand on, which is where another kind of evidence
        # has to attach to be evidence about the same chain
        shown: dict[int, set[str]] = collections.defaultdict(set)
        element_row: dict[str, dict[str, Any]] = {}
        element_iv: dict[str, tuple[str, int, int]] = {}
        element_cells: dict[str, set[str]] = collections.defaultdict(set)
        element_genes: dict[str, set[str]] = collections.defaultdict(set)
        for g, c in zip(mine, cands, strict=True):
            if g in addable and c.status == rm2.STATUS_ADDABLE and c.element is not None:
                shown[g].add(c.element["id"])
                element_row[c.element["id"]] = c.row or {}
                element_iv[c.element["id"]] = (c.pair.chrom, c.element["start"], c.element["end"])
                element_cells[c.element["id"]].add(c.pair.cell)
                element_genes[c.element["id"]].add(c.pair.gene)
        crispri_groups = set(mine)

        # --- the measured layer's whole assay inventory, for the further-locus arm ----------------
        assay_elements: collections.Counter[str] = collections.Counter()
        other_keys: dict[str, list[cell2.LocusKey]] = {a: [] for a in ms.ASSAYS}
        for chrom in chroms:
            elements = cp._attributed(chrom, RESULTS_DIR)
            _layer, measured_rows = cp._measured_rows(chrom, elements, RESULTS_DIR)
            for r in measured_rows:
                for a in r["assays"]:
                    assay_elements[a] += 1
                    if a == "crispri":
                        continue
                    other_keys[a].append(cell2.LocusKey(a, chrom, r["start"], r["end"], f"{a}|{r['id']}"))

        # --- eQTL: the retained set, its frame, and what lands on the shown elements ---------------
        index, indexed_n = eqtl_index(ROOT)
        hits = eqtl_hits_by_chrom(eqtl.KNOWLEDGE)
        eqtl_shown: dict[str, list[dict[str, Any]]] = {}
        eqtl_assessable: set[str] = set()
        for eid, (chrom, s, e) in element_iv.items():
            if covered_by(index.get(chrom, []), s - EQTL_MARGIN, e + EQTL_MARGIN):
                eqtl_assessable.add(eid)
            got = hits_in(hits.get(chrom, []), s, e)
            if got:
                eqtl_shown[eid] = got
        eqtl_chroms = {element_iv[eid][0] for eid in eqtl_shown}
        symbols = symbol_map(eqtl_chroms) if eqtl_chroms else {}
        eqtl_loci = {g for g, eids in shown.items() if eids & set(eqtl_shown)}
        eqtl_assessable_loci = {g for g, eids in shown.items() if eids & eqtl_assessable}

        # the further eQTL arm, from the committed result: its element ids and its eGene symbols
        with open(ROOT / EQTL_RESULT) as fh:
            eq = json.load(fh)
        eqtl_keys: list[cell2.LocusKey] = []
        eqtl_further_elements: set[str] = set()
        for blk in eq["sets"].values():
            for r in blk["elements_with_eqtl"]:
                chrom = r["chrom"]
                s, e = r["id"].split(":")[1].split("-")
                eqtl_further_elements.add(r["id"])
                for gene in r["egenes"]:
                    eqtl_keys.append(cell2.LocusKey("gtex_donors", chrom, int(s), int(e), gene))

    # --- further loci, one grouping per kind, against every CRISPRi key ----------------------------
    crispri_keys = [c.key for c in cands]
    further: dict[str, dict[str, Any]] = {}
    for kind, keys in [*other_keys.items(), ("eqtl", eqtl_keys)]:
        if kind == "crispri":
            continue
        if not keys:
            further[kind] = {"observations": 0, "further_independent_loci": 0}
            continue
        g2 = cell2.group([*crispri_keys, *keys])
        theirs = g2[len(crispri_keys) :]
        reached = set(g2[: len(crispri_keys)])
        further[kind] = {
            "observations": len(keys),
            "independent_loci_of_this_kind": len(set(theirs)),
            "of_those_the_crispri_arm_already_reaches": len(set(theirs) & reached),
            "further_independent_loci": len(set(theirs) - reached),
        }

    loci_total = len(set(mine))
    n_addable = len(addable)
    n_shown = len(element_row)

    def on_the_110(assay: str) -> dict[str, Any]:
        eids = {eid for eid, row in element_row.items() if assay in (row.get("assays") or ())}
        loci = {g for g, s in shown.items() if s & eids}
        return {"shown_elements": len(eids), "loci": len(loci), "element_ids": sorted(eids)[:40]}

    kinds: list[dict[str, Any]] = [
        {
            "kind": "crispri",
            "what_the_observation_is": ms.OUTCOME_KIND["crispri"],
            "source": ms.SOURCES["crispri"],
            "an_observation_of_the_native_locus": True,
            "why": "the element is silenced in place and the gene's expression is read out",
            "already_in_increment_2": True,
            "loci_on_the_existing_110": n_addable,
            "observations_increment_2_showed": 195,
        },
        {
            "kind": "reader_dnase",
            "what_the_observation_is": (
                f"chromatin accessibility of the element in the biosample, at {ce.OPENNESS_CALL}"
            ),
            "source": "ENCODE DNase-seq narrowPeak, GRCh38, as reader v1 cached it",
            "an_observation_of_the_native_locus": True,
            "why": "the locus's own state is read out on the chromosome in place",
            "already_in_increment_2": True,
            "loci_on_the_existing_110": n_addable,
            "observations_increment_2_showed": 166,
            "not_one_of_the_four_assays": (
                "the reader is not part of the measured layer's four assays; it is the third input of "
                "the three-input rule, and increment 2 counted its 166 readings apart from the 195 "
                "CRISPRi effects for exactly this reason"
            ),
        },
    ]
    for assay in ("lentimpra", "vista", "satmut"):
        blk = on_the_110(assay)
        kinds.append(
            {
                "kind": assay,
                "what_the_observation_is": ms.OUTCOME_KIND[assay],
                "source": ms.SOURCES[assay],
                "an_observation_of_the_native_locus": False,
                "why_not": (
                    "it reads a sequence, or the bases inside it, not the native locus after a "
                    "perturbation of it - increment 2's own wording, carried unchanged"
                ),
                "native_locus_observations_on_the_existing_110": 0,
                "native_locus_loci_on_the_existing_110": 0,
                "sequence_observations_on_the_existing_110": blk,
                "what_those_sequence_observations_say": reporter_outcomes(
                    {eid for eid, row in element_row.items() if assay in (row.get("assays") or ())},
                    element_row,
                    dict(element_cells),
                    assay,
                ),
                "what_the_sequence_observation_can_be_read_as": (
                    "evidence about the element's SEQUENCE in a reporter, named by what measured it, "
                    "and never as the locus having been perturbed"
                ),
                "further": further.get(assay, {}),
                "elements_carrying_it_genome_wide": assay_elements[assay],
            }
        )
        if assay == "satmut":
            kinds[-1]["cannot_disagree"] = ms.SATMUT_CANNOT_DISAGREE
        if assay == "lentimpra":
            kinds[-1]["label_rule"] = ms.REPORTER_LABEL_RULE
            kinds[-1]["observation_unit"] = list(ms.REPORTER_UNIT)
            kinds[-1]["uncertainty"] = ms.REPORTER_UNCERTAINTY
    kinds.append(
        {
            "kind": "eqtl",
            "what_the_observation_is": (
                "an association, not a perturbation: a significant cis-eQTL is a variant inside the "
                "element whose alleles go with a gene's expression across GTEx donors, in one tissue"
            ),
            "source": eqtl.EVIDENCE,
            "an_observation_of_the_native_locus": True,
            "why": (
                "the variant sits on the chromosome in place and the expression is of the donor's own "
                "gene, so the reading is of the native locus. It is natural variation and an "
                "association, never an experimental perturbation, and association is not mechanism"
            ),
            "loci_on_the_existing_110": len(eqtl_loci),
            "shown_elements_with_a_retained_hit": len(eqtl_shown),
            "shown_elements_the_retained_set_could_have_seen_at_all": len(eqtl_assessable),
            "loci_the_retained_set_could_have_seen_at_all": len(eqtl_assessable_loci),
            "attachment_rule": (
                "a retained hit whose variant position falls inside the compiled element widened by "
                f"{EQTL_MARGIN} bases each side, which is the margin the retained set was distilled "
                "under (scripts/eqtl_targets.MARGIN), imported and not chosen here"
            ),
            "the_frame_is_the_limit": (
                f"the retained set was distilled against {indexed_n:,} intervals of three sampled "
                "element sets (enhancer_targets_chr*, constrained_targets_chr*, vista_chr*), so a "
                "variant outside that frame was discarded before this count existed. "
                f"{len(eqtl_assessable)} of the {n_shown} shown elements lie inside the frame; for "
                "the rest, no eQTL could have been seen whatever is really there, and their absence "
                "is UNASSESSED and not an absence of eQTLs"
            ),
            "what_those_observations_say": eqtl_outcomes(eqtl_shown, dict(element_genes), symbols),
            "further": further.get("eqtl", {}),
            "further_elements_with_an_eqtl_in_the_committed_result": len(eqtl_further_elements),
            "elements_by_set": {k: v["summary"]["with_eqtl"] for k, v in eq["sets"].items()},
        }
    )

    # --- how many kinds of evidence each of the addable loci carries -------------------------------
    # Every addable locus carries CRISPRi and a reader state by construction, which is what made it
    # addable; the cross-tabulation is of what it carries BEYOND those two.
    lenti_loci = {
        g
        for g, eids in shown.items()
        if any("lentimpra" in (element_row[e].get("assays") or ()) for e in eids)
    }
    extra: collections.Counter[str] = collections.Counter()
    for g in sorted(addable):
        tags = []
        if g in lenti_loci:
            tags.append("lentimpra_sequence")
        if g in eqtl_loci:
            tags.append("eqtl_native_locus")
        extra["+".join(tags) if tags else "nothing_beyond_the_two"] += 1
    kinds_per_locus = {
        "every_addable_locus_already_carries": [
            "crispri: a measured perturbation of the native locus",
            "reader_dnase: the locus's own accessibility in the biosample",
        ],
        "beyond_those_two": dict(extra),
        "loci_carrying_at_least_one_further_kind": len(lenti_loci | eqtl_loci),
        "loci_carrying_both_further_kinds": len(lenti_loci & eqtl_loci),
        "loci_carrying_nothing_further": n_addable - len(lenti_loci | eqtl_loci),
        "vista_and_satmut_contribute_to_none": True,
        "a_further_kind_is_not_a_second_opinion": (
            "a lentiMPRA reading is of the element's sequence in a reporter and a GTEx eQTL is an "
            "association over donor tissues. Neither is a second measurement of the CRISPRi "
            "observation, and a locus carrying three kinds is not a locus three times established"
        ),
    }

    native_kinds = [k for k in kinds if k.get("an_observation_of_the_native_locus")]
    payload: dict[str, Any] = {
        "question": (
            "how many of the 110 loci increment 2 built carry observations from the measured layer's "
            "other three assays, and from GTEx eQTLs; how many further loci those kinds of evidence "
            "would add that the CRISPRi arm does not reach; and what each of those observations is"
        ),
        "counting_only": (
            "no increment is built, nothing is added to the map, nothing is scored and nothing is "
            "fitted. No cut-off is introduced: every rule, span and call is imported from where it "
            "was registered. 0 model requests, no downloads, no money: every input was on disk"
        ),
        "starting_point": {
            "increment_2_stated_limit": INCREMENT_2_LIMIT,
            "so_increment_2_is_one_assay_deep": True,
        },
        "chromosomes": chroms,
        "imported_counting_functions": list(IMPORTED_COUNTING),
        "native_locus_rule": NATIVE_LOCUS_RULE,
        "outcome_breakdown_rule": OUTCOME_BREAKDOWN_RULE,
        "no_measured_gene_token": NO_MEASURED_GENE,
        "denominator": {
            "measured_perturbation_pairs": counted["pairs"],
            "invalid_rows_the_parser_rejected": counted["invalid_rows_the_parser_rejected"],
            "independent_loci_of_the_measured_layer": loci_total,
            "independent_locus_rule": cell2.INDEPENDENT_LOCUS_RULE,
            "not_biological_independence": (
                "the grouping is operational and pools cell types: it joins two keys on a shared "
                "measured gene or on 1 Mb proximity on one chromosome and never reads the cell. It is "
                "not established biological independence and nothing may cite it as one"
            ),
            "addable_loci_re_derived_here": n_addable,
            "shown_elements_of_those_loci": n_shown,
            "by_status": counted["by_status"],
            "elements_by_assay_genome_wide": dict(sorted(assay_elements.items())),
        },
        "kinds_of_evidence": kinds,
        "how_many_kinds_reach_a_native_locus": {
            "kinds": [k["kind"] for k in native_kinds],
            "count": len(native_kinds),
            "already_in_increment_2": ["crispri", "reader_dnase"],
            "new_here": [k["kind"] for k in native_kinds if not k.get("already_in_increment_2")],
        },
        "kinds_per_locus": kinds_per_locus,
        "what_a_build_would_gain": {
            "loci_gaining_a_sequence_observation_from_lentimpra": on_the_110("lentimpra")["loci"],
            "loci_gaining_a_native_locus_eqtl_observation": len(eqtl_loci),
            "loci_gaining_nothing_from_vista_or_satmut": n_addable,
            "loci_gaining_at_least_one_further_kind": kinds_per_locus[
                "loci_carrying_at_least_one_further_kind"
            ],
            "a_count_is_not_a_measurement": (
                "these are counts of loci where an observation exists, not of loci where the "
                "observations agree, and nothing here says what any of them shows. The outcome "
                "breakdown of each kind is in its own block above"
            ),
        },
        "reconciliation": {
            "addable_loci": n_addable,
            "addable_loci_increment_2_committed": 110,
            "holds": n_addable == 110,
            "how": (
                "the 110 are re-derived here by the imported counting function and its imported "
                "grouping, not read out of increment 2's result, so this equality is a check on the "
                "method and not a restatement"
            ),
            "shown_elements": n_shown,
            "crispri_groups_cover_every_pair": len(crispri_groups) == loci_total,
        },
        "cannot_establish": [
            "nothing here is a verdict. A second kind of evidence at a locus is not agreement, "
            "confirmation or validation of the first, and this count says nothing about whether any "
            "two of them point the same way",
            "a lentiMPRA, VISTA or saturation-mutagenesis reading is an observation of a sequence in "
            "a reporter. It is NOT an observation of the native locus, it cannot be added to the "
            "evidence for a perturbation of that locus, and the zero above is the rule holding rather "
            "than a gap in the data",
            "the eQTL count is a floor over a sampled frame: only "
            f"{len(eqtl_assessable)} of the {n_shown} shown elements lie inside the "
            f"{indexed_n:,} intervals the retained set was distilled against, so for the remainder "
            "no eQTL could have been seen and their absence is unassessed",
            "a cis-eQTL is an association over donors. Linkage carries the signal some distance, so "
            "the variant inside the element need not be the causal one, and the eGene need not be "
            "regulated by this element",
            "the independent-locus grouping is an operational grouping, not established biological "
            "independence, and it pools cell types",
            "the further-locus counts for the three reporter assays are upper bounds, because a "
            "reporter observation names no measured gene and only the proximity arm of the grouping "
            "can fire for it",
            "GTEx tissues are not the CRISPRi cell lines. An eQTL in a GTEx tissue and a CRISPRi "
            "effect in K562 are observations in different contexts, and no cross-context agreement is "
            "counted, claimed or implied here",
        ],
        "per_element_response_cache": audit.block(),
        "alphagenome_requests": 0,
        "money": "none: every input was already on disk",
        "seconds": round(time.time() - started, 1),
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }
    if not payload["reconciliation"]["holds"]:
        raise SystemExit(
            f"refusing to save: the imported counting function gives {n_addable} addable loci, not "
            "the 110 increment 2 committed, so this count is not of the same population"
        )

    inputs = [mf.input_entry(ce.TRACK_METADATA, partition=None)]
    for f in ms.CRISPRI_FILES:
        p = ms.CRISPRI_KNOWLEDGE / f
        if p.exists():
            inputs.append(mf.input_entry(p, partition=ms.CRISPRI_SPLIT_OF[f]))
    seen = {str(i["path"]) for i in inputs}
    paths = [
        p
        for cell in ce.READER_TERMS
        for chrom in chroms
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
        for chrom in chroms:
            paths += sorted(RESULTS_DIR.glob(g.replace("chr*", chrom)))
    paths.append(Path(EQTL_RESULT))
    paths += sorted(eqtl.KNOWLEDGE.glob("hits_*.tsv"))
    p = eqtl.KNOWLEDGE / "distil_summary.json"
    if p.exists():
        paths.append(p)
    for group in OPENED_BESIDE_THE_RESULTS.values():
        # repository-relative, so `manifest_rebuild.py` can check the bytes its own
        # worktree reads rather than an absolute path outside the linked stores
        paths += [Path(x) for x in group]
    for chrom in sorted(eqtl_chroms):
        paths.append(Path(GENCODE.format(chrom=chrom)))
    for q in paths:
        if q.as_posix() not in seen:
            seen.add(q.as_posix())
            inputs.append(mf.input_entry(q, partition=None))

    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "EngreitzLab/CRISPR_comparison (Gschwind et al. 2025)",
                "version": "the two local benchmark tables under data/knowledge/crispri",
            },
            {
                "accession": "ENCODE DNase-seq narrowPeak, GRCh38, released, 13 biosamples",
                "version": "as reader v1 cached them under data/results/dnase_*_chr*.bed.gz",
            },
            {
                "accession": "GTEx v8 single-tissue cis-eQTLs, significant variant-gene pairs",
                "version": (
                    "the local distillation under data/knowledge/gtex, kept only inside the three "
                    "sampled element sets; see data/knowledge/gtex/distil_summary.json"
                ),
            },
            {
                "accession": "the compiled non-coding programs of the 24 chromosomes",
                "version": (
                    "the compiler's own element and measured-layer enumeration, run in this process "
                    "from the results on disk, not read back from a .bio file"
                ),
            },
        ],
        "inputs": inputs,
        "inputs_opened_beside_the_declared_results": {
            "why_they_are_listed": (
                "an input is what the writer opened, not what it said it read. The audit hook in "
                "`genomeos.manifest` records every open under `data/` and `save_result` refuses a "
                "result whose reads exceed its declared inputs, which is how these came to be named"
            ),
            "groups": {k: len(v) for k, v in OPENED_BESIDE_THE_RESULTS.items()},
            "the_pointer_case": (
                "`attribution.targets` reads a run's elements out of the table its result's "
                "`elements_where` field points at, so the path comes out of a result file at run time "
                "and no reading of the source could name it"
            ),
            "gencode_chromosomes": sorted(eqtl_chroms),
        },
        "inputs_are_recorded_one_file_per_entry": (
            "scripts/manifest_rebuild.py resolves an input by its `path`, and a grouped entry carries "
            "a label there, so a grouped input would read as absent without any file being checked"
        ),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "independent_locus_span": cell2.INDEPENDENT_LOCUS_SPAN,
            "reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
            "openness_call": ce.OPENNESS_CALL,
            "eqtl_margin_bp": EQTL_MARGIN,
            "mpra_active_log2": ms.MPRA_ACTIVE,
            "chromosomes": chroms,
        },
        "exclusions": [
            "no pair is excluded: the denominator is every CRISPRi pair the measured layer parses as "
            "valid on these chromosomes, both arms and every cell type",
            "lentiMPRA, VISTA and saturation mutagenesis are excluded from the native-locus count "
            "under a named cause - they read a sequence, not the native locus after a perturbation - "
            "and their sequence observations are counted separately rather than dropped silently",
            "an eQTL outside the distillation's sampled frame is excluded as unassessed, by name, and "
            "is not counted as an absence of eQTLs",
        ],
        "partitions": {
            "kinds_of_evidence": "one block per kind, with what its observation is",
            "native_locus": "whether the kind reads the native locus or a copy of its sequence",
        },
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }

    path = save_result(name, payload)
    print(f"{name}: {path}")
    print(f"  addable loci re-derived: {n_addable} (increment 2 committed 110)")
    print(f"  shown elements: {n_shown}")
    for k in kinds:
        if k["kind"] in ("crispri", "reader_dnase"):
            continue
        if k["kind"] == "eqtl":
            print(
                f"  eqtl: {k['loci_on_the_existing_110']} loci, "
                f"{k['shown_elements_with_a_retained_hit']} of "
                f"{k['shown_elements_the_retained_set_could_have_seen_at_all']} assessable elements; "
                f"further loci {k['further'].get('further_independent_loci')}"
            )
        else:
            s = k["sequence_observations_on_the_existing_110"]
            print(
                f"  {k['kind']}: native-locus 0; sequence observations on {s['loci']} loci "
                f"({s['shown_elements']} elements); further loci "
                f"{k['further'].get('further_independent_loci')}"
            )
    print(f"  per-element cache opens: {audit.block()['opens']}")


if __name__ == "__main__":
    main()
